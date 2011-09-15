/**********************************************************************
Copyright (C) 2006 by Fredrik Wallner
Some portions Copyright (C) 2006-2007 by Geoffrey Hutchsion
Some portions Copyright (C) 2004 by Chris Morley

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include <openbabel/babelconfig.h>
#include <openbabel/obmolecformat.h>
#include <openbabel/reaction.h>
#include "chemdrawcdx.h"

#include <iostream>
#include <fstream>
#include <map>
#include <list>

//#define debug

#if !defined(__CYGWIN__)
static inline unsigned short bswap_16(unsigned short x) {
  return (x>>8) | (x<<8);
}

static inline unsigned int bswap_32(unsigned int x) {
  return (bswap_16(x&0xffff)<<16) | (bswap_16(x>>16));
}

static inline unsigned long long bswap_64(unsigned long long x) {
  return (((unsigned long long)bswap_32(x&0xffffffffull))<<32) | (bswap_32(x>>32));
}
#endif

// Macs -- need to use Apple macros to deal with Universal binaries correctly
#ifdef __APPLE__
#include <machine/endian.h>
#if BYTE_ORDER == BIG_ENDIAN
#    define READ_INT16(stream,data) \
(stream).read ((char*)&data, sizeof(data)); \
data = bswap_16 (data);
#    define READ_INT32(stream,data) \
(stream).read ((char*)&data, sizeof(data)); \
data = bswap_32 (data);
#else BYTE_ORDER == LITTLE_ENDIAN
#    define READ_INT16(stream,data) \
(stream).read ((char*)&data, sizeof(data));
#    define READ_INT32(stream,data) \
(stream).read ((char*)&data, sizeof(data));
#endif
#else

// Non-Apple systems
// defined in babelconfig.h by autoconf (portable to Solaris, BSD, Linux)
#ifdef WORDS_BIGENDIAN
#    define READ_INT16(stream,data) \
(stream).read ((char*)&data, sizeof(data)); \
data = bswap_16 (data);
#    define READ_INT32(stream,data) \
(stream).read ((char*)&data, sizeof(data)); \
data = bswap_32 (data);
#else
#    define READ_INT16(stream,data) \
(stream).read ((char*)&data, sizeof(data));
#    define READ_INT32(stream,data) \
(stream).read ((char*)&data, sizeof(data));
#endif
// end endian / bigendian issues (on non-Mac systems)
#endif
// end Apple/non-Apple systems

using namespace std;
namespace OpenBabel
{

  class ChemDrawBinaryFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID in the constructor
    ChemDrawBinaryFormat()
    {
      OBConversion::RegisterFormat("cdx",this, "chemical/x-cdx");
    }

    virtual const char* Description() //required
    {
      return
        "ChemDraw binary format\n"
        "Read only.\n";
    };

    //Optional URL where the file format is specified
    virtual const char* SpecificationURL()
    {return "http://www.cambridgesoft.com/services/documentation/sdk/chemdraw/cdx/IntroCDX.htm";};

    virtual const char* GetMIMEType()
    { return "chemical/x-cdx"; };

    /* Flags() can return be any of the following combined by |
       or be omitted if none apply
       NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY */
    virtual unsigned int Flags()
    {
      return READBINARY|NOTWRITABLE;
    };

    ////////////////////////////////////////////////////
    /// Declarations for the "API" interface functions. Definitions are below
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
    OBBase* ReadObject(OBConversion* pConv);
    //  virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);
  private:
    struct cdBond
    {
      UINT32	begin;
      UINT32	end;
      int	order;
      int stereo;

			cdBond() {}
			cdBond(UINT32 bgn, UINT32 e, int ord, int st = 0) : begin(bgn), end (e), order(ord), stereo(st) {}
			~cdBond() {}
    };


    int readFragment(istream *ifs, UINT32 fragmentId,
                     OBMol *pmol, map<UINT32, int> &atoms, list<cdBond> &bonds);
    int readNode(istream *ifs, UINT32 nodeId, OBMol *pmol,
                 map<UINT32, int> &atoms, list<cdBond> &bonds,
                 UINT32 fragmentID);
    int readBond(istream *ifs, UINT32 nodeId, OBMol *pmol, list<cdBond> &bonds);
    int readGeneric(istream *ifs, UINT32 objId);
    const char* getName(istream *ifs, UINT32 size);

  };
	////////////////////////////////////////////////////

  //Make an instance of the format class
  ChemDrawBinaryFormat theChemDrawBinaryFormat;

  /////////////////////////////////////////////////////////////////
  bool ChemDrawBinaryFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
      return false;

    istream& ifs = *pConv->GetInStream();
    char buffer[BUFF_SIZE];
    char errorMsg[BUFF_SIZE];
    UINT16 tag;
    UINT16 size;
    UINT32 id;
    map<UINT32, int> atoms;
    list<cdBond> bonds;
    list<cdBond>::const_iterator bondsIter;
    cdBond prevBond;
    cdBond oneBond;
    //    OBPairData *pd;
    int iInt=0, depth=1;

    if (!ifs.good() || ifs.peek() == EOF)
      return false;

    pmol->BeginModify();

    pmol->SetTitle(pConv->GetTitle());

    if((int)ifs.tellg() == 0)	// Beginning of file
      {
        ifs.read(buffer,kCDX_HeaderStringLen);
        if(strncmp(buffer, kCDX_HeaderString, kCDX_HeaderStringLen) == 0)  // File header
          {
            ifs.seekg (kCDX_HeaderLength - kCDX_HeaderStringLen, ios_base::cur);	// Discard rest of header.
          }
        else
          {
            cout << "Invalid file, no ChemDraw Header" << endl;	// No header, error.
            ifs.seekg(0, ios_base::end);
            return false;
          }
      }
    while(ifs.good())
      {
        READ_INT16 (ifs, tag);
        if(tag & kCDXTag_Object)	// Object
          {
            READ_INT32 (ifs, id);
            snprintf(errorMsg, BUFF_SIZE, "Object ID: %08X in root has type: %04X\n", id, tag);
            obErrorLog.ThrowError(__FUNCTION__, errorMsg, obDebug);

            if(tag == kCDXObj_Fragment)
              {
                if(readFragment(&ifs, id, pmol, atoms, bonds) != 0)
                  {
                    obErrorLog.ThrowError(__FUNCTION__, "Error reading fragment", obWarning);
                    return false;
                  }
                iInt = 0;
                prevBond = cdBond(0,0,0);
                for (bondsIter=bonds.begin(); bondsIter != bonds.end(); bondsIter++)
                  {
                    // printf("Bond between %08X (%d) and %08x (%d)\n", bondsIter->begin, atoms[bondsIter->begin], bondsIter->end, atoms[bondsIter->end]);
                    if(atoms[bondsIter->begin] == -1)
                      {
                        //						printf("Bond starts at a non-atom\n");
                        if(atoms[bondsIter->end] == -1)
                          {
                            //							printf("Bond between two non-atoms pushed back\n");
                            if((prevBond.begin == bondsIter->begin) && (prevBond.end == bondsIter->end))
                              {
                                obErrorLog.ThrowError(__FUNCTION__, "Can't assign all bonds!", obDebug);
                                break;
                              }
                            bonds.push_back(*bondsIter);
                          }
                        else
                          atoms[bondsIter->begin] = atoms[bondsIter->end];
                      }
                    else if(atoms[bondsIter->end] == -1)
                      atoms[bondsIter->end] = atoms[bondsIter->begin];
                    else
                      {
                        int flags = 0;
                        switch (bondsIter->stereo)
                          {
                          case kCDXBondDisplay_WedgedHashBegin:
                          case kCDXBondDisplay_WedgeEnd:
                            flags = OB_HASH_BOND;
                            break;
                          case kCDXBondDisplay_WedgedHashEnd:
                          case kCDXBondDisplay_WedgeBegin:
                            flags = OB_WEDGE_BOND;
                            break;
                          default:;
                          }
                        pmol->AddBond(atoms[bondsIter->begin], atoms[bondsIter->end], bondsIter->order, flags);
                      }
                    prevBond = *bondsIter;
                  }
                break;
              }
            else if ((tag == kCDXObj_Graphic) || (tag == kCDXObj_Text)
                     || (tag == kCDXObj_BracketedGroup)
                     || (tag == kCDXObj_BracketAttachment)
                     || (tag == kCDXObj_CrossingBond)
                     || (tag == kCDXObj_ReactionStep)
                     || (tag == kCDXObj_Curve)
                     || (tag == kCDXObj_EmbeddedObject))
              {	// Objects that can be ignored
                readGeneric(&ifs, id);
              }
            else if((tag == kCDXObj_Page) || ( tag == kCDXObj_Group))
              {	// Objects where the header can be ignored
              }
            else
              {
                depth++;
                snprintf(errorMsg, BUFF_SIZE, "New object in root, type %04X\n", tag);
                obErrorLog.ThrowError(__FUNCTION__, errorMsg, obDebug);
              }
          }
        else if(tag == 0)	// End of object
          {
            if(depth > 1)
              {
                obErrorLog.ThrowError(__FUNCTION__, errorMsg, obDebug);
                depth--;
              }
          }
        else	// Property
          {
            READ_INT16 (ifs ,size);
            switch(tag)
              {
              case kCDXProp_Name:
                pmol->SetTitle(getName(&ifs, size));
                break;
              case kCDXProp_FontTable:
              case kCDXProp_BoundingBox:
              case kCDXProp_Window_Position:
              case kCDXProp_Window_Size:
              case kCDXProp_Atom_ShowQuery:
              case kCDXProp_Atom_ShowStereo:
              case kCDXProp_Atom_ShowAtomNumber:
              case kCDXProp_Bond_ShowQuery:
              case kCDXProp_Bond_ShowStereo:
              case kCDXProp_MacPrintInfo:
              case kCDXProp_DrawingSpaceType:
              case kCDXProp_Width:
              case kCDXProp_Height:
              case kCDXProp_PageOverlap:
              case kCDXProp_Magnification:
              case kCDXProp_ChainAngle:
              case kCDXProp_BondSpacing:
              case kCDXProp_BondSpacingAbs:
              case kCDXProp_BondLength:
              case kCDXProp_BoldWidth:
              case kCDXProp_LineWidth:
              case kCDXProp_MarginWidth:
              case kCDXProp_HashSpacing:
              case kCDXProp_LabelStyle:
              case kCDXProp_CaptionStyle:
              case kCDXProp_CaptionJustification:
              case kCDXProp_LabelJustification:
              case kCDXProp_FractionalWidths:
              case kCDXProp_LabelLineHeight:
              case kCDXProp_CaptionLineHeight:
              case kCDXProp_Window_IsZoomed:
              case kCDXProp_WidthPages:
              case kCDXProp_HeightPages:
              case kCDXProp_HeaderPosition:
              case kCDXProp_FooterPosition:
              case kCDXProp_PrintTrimMarks:
              case kCDXProp_PrintMargins:
              case kCDXProp_ColorTable:
              case kCDXProp_CreationProgram:
              case kCDXProp_Group_Integral:
                ifs.seekg(size, ios_base::cur);
                break;
              default: // some unknown tag
                ifs.seekg(size, ios_base::cur);
                snprintf(errorMsg, BUFF_SIZE, "Root Tag: %04X\tSize: %04X\n", tag, size);
                obErrorLog.ThrowError(__FUNCTION__, errorMsg, obDebug);
                break;
              }
          }
      }

    pmol->SetDimension(2);
    pmol->EndModify();

    return true;
  }

  OBBase* ChemDrawBinaryFormat::ReadObject(OBConversion* pConv)
  {
    istream& ifs = *pConv->GetInStream();
    char buffer[BUFF_SIZE];
    char errorMsg[BUFF_SIZE];
    UINT16 tag;
    UINT16 size;
    UINT32 id;
    map<UINT32, int> atoms;
    list<cdBond> bonds;
    list<cdBond>::const_iterator bondsIter;
    cdBond prevBond;
    cdBond oneBond;
    //    OBPairData *pd;
    int iInt=0, depth=1;

    if((int)ifs.tellg() == 0)	// Beginning of file
      {
        ifs.read(buffer,kCDX_HeaderStringLen);
        if(strncmp(buffer, kCDX_HeaderString, kCDX_HeaderStringLen) == 0)  // File header
          {
            ifs.seekg (kCDX_HeaderLength - kCDX_HeaderStringLen, ios_base::cur);	// Discard rest of header.
          }
        else
          {
            cout << "Invalid file, no ChemDraw Header" << endl;	// No header, error.
            ifs.seekg(0, ios_base::end);
            return (OBBase*)NULL;
          }
      }
    while(ifs.good())
      {
        READ_INT16 (ifs, tag);
        if(tag & kCDXTag_Object)	// Object
          {
            READ_INT32 (ifs, id);
            snprintf(errorMsg, BUFF_SIZE, "Object ID: %08X in root has type: %04X\n", id, tag);
            obErrorLog.ThrowError(__FUNCTION__, errorMsg, obDebug);

            if(tag == kCDXObj_Fragment)
              {
                OBMol *pmol = new OBMol();
                if(pmol==NULL)
                  {
                  	return NULL;
                  }

                pmol->BeginModify();

  	            pmol->SetTitle(pConv->GetTitle());
                if(readFragment(&ifs, id, pmol, atoms, bonds) != 0)
                  {
                    obErrorLog.ThrowError(__FUNCTION__, "Error reading fragment", obWarning);
                    delete pmol;
                    return NULL;
                  }
                iInt = 0;
                prevBond = cdBond(0,0,0);
                for (bondsIter=bonds.begin(); bondsIter != bonds.end(); bondsIter++)
                  {
                    // printf("Bond between %08X (%d) and %08x (%d)\n", bondsIter->begin, atoms[bondsIter->begin], bondsIter->end, atoms[bondsIter->end]);
                    if(atoms[bondsIter->begin] == -1)
                      {
                        //						printf("Bond starts at a non-atom\n");
                        if(atoms[bondsIter->end] == -1)
                          {
                            //							printf("Bond between two non-atoms pushed back\n");
                            if((prevBond.begin == bondsIter->begin) && (prevBond.end == bondsIter->end))
                              {
                                obErrorLog.ThrowError(__FUNCTION__, "Can't assign all bonds!", obDebug);
                                break;
                              }
                            bonds.push_back(*bondsIter);
                          }
                        else
                          atoms[bondsIter->begin] = atoms[bondsIter->end];
                      }
                    else if(atoms[bondsIter->end] == -1)
                      atoms[bondsIter->end] = atoms[bondsIter->begin];
                    else
                      {
                        int flags = 0;
                        switch (bondsIter->stereo)
                          {
                          case kCDXBondDisplay_WedgedHashBegin:
                          case kCDXBondDisplay_WedgeEnd:
                            flags = OB_HASH_BOND;
                            break;
                          case kCDXBondDisplay_WedgedHashEnd:
                          case kCDXBondDisplay_WedgeBegin:
                            flags = OB_WEDGE_BOND;
                            break;
                          default:;
                          }
                        pmol->AddBond(atoms[bondsIter->begin], atoms[bondsIter->end], bondsIter->order, flags);
                      }
                    prevBond = *bondsIter;
                  }
                pmol->SetDimension(2);
                pmol->EndModify();
                return pmol;
              }
            else if (tag == kCDXObj_ReactionScheme)
              {
                readGeneric(&ifs, id);
puts("found a reaction");
                return new OBReaction();
              }
           else if ((tag == kCDXObj_Graphic) || (tag == kCDXObj_Text)
                     || (tag == kCDXObj_BracketedGroup)
                     || (tag == kCDXObj_BracketAttachment)
                     || (tag == kCDXObj_CrossingBond)
                     || (tag == kCDXObj_Curve)
                     || (tag == kCDXObj_EmbeddedObject))
              {	// Objects that can be ignored
                readGeneric(&ifs, id);
              }
            else if((tag == kCDXObj_Page) || ( tag == kCDXObj_Group))
              {	// Objects where the header can be ignored
              }
            else
              {
                depth++;
                snprintf(errorMsg, BUFF_SIZE, "New object in root, type %04X\n", tag);
                obErrorLog.ThrowError(__FUNCTION__, errorMsg, obDebug);
              }
          }
        else if(tag == 0)	// End of object
          {
            if(depth > 1)
              {
                obErrorLog.ThrowError(__FUNCTION__, errorMsg, obDebug);
                depth--;
              }
          }
        else	// Property
          {
            READ_INT16 (ifs ,size);
            switch(tag)
              {
              case kCDXProp_Name:
puts("found name");
                /*pmol->SetTitle(*/puts(getName(&ifs, size));//);
                break;
              case kCDXProp_FontTable:
              case kCDXProp_BoundingBox:
              case kCDXProp_Window_Position:
              case kCDXProp_Window_Size:
              case kCDXProp_Atom_ShowQuery:
              case kCDXProp_Atom_ShowStereo:
              case kCDXProp_Atom_ShowAtomNumber:
              case kCDXProp_Bond_ShowQuery:
              case kCDXProp_Bond_ShowStereo:
              case kCDXProp_MacPrintInfo:
              case kCDXProp_DrawingSpaceType:
              case kCDXProp_Width:
              case kCDXProp_Height:
              case kCDXProp_PageOverlap:
              case kCDXProp_Magnification:
              case kCDXProp_ChainAngle:
              case kCDXProp_BondSpacing:
              case kCDXProp_BondSpacingAbs:
              case kCDXProp_BondLength:
              case kCDXProp_BoldWidth:
              case kCDXProp_LineWidth:
              case kCDXProp_MarginWidth:
              case kCDXProp_HashSpacing:
              case kCDXProp_LabelStyle:
              case kCDXProp_CaptionStyle:
              case kCDXProp_CaptionJustification:
              case kCDXProp_LabelJustification:
              case kCDXProp_FractionalWidths:
              case kCDXProp_LabelLineHeight:
              case kCDXProp_CaptionLineHeight:
              case kCDXProp_Window_IsZoomed:
              case kCDXProp_WidthPages:
              case kCDXProp_HeightPages:
              case kCDXProp_HeaderPosition:
              case kCDXProp_FooterPosition:
              case kCDXProp_PrintTrimMarks:
              case kCDXProp_PrintMargins:
              case kCDXProp_ColorTable:
              case kCDXProp_CreationProgram:
              case kCDXProp_Group_Integral:
                ifs.seekg(size, ios_base::cur);
                break;
              default: // some unknown tag
                ifs.seekg(size, ios_base::cur);
                snprintf(errorMsg, BUFF_SIZE, "Root Tag: %04X\tSize: %04X\n", tag, size);
                obErrorLog.ThrowError(__FUNCTION__, errorMsg, obDebug);
                break;
              }
          }
      }

    return NULL;
  }

  const char* ChemDrawBinaryFormat::getName(istream *ifs, UINT32 size)
  {
    char *buff;
    UINT16 temp;
    UINT32 skipsize;

    //	buff = new char[size];
    READ_INT16 ((*ifs), temp);
    if(temp == 0)
      {
        buff = new char[size-1];
        ifs->read(buff, size-2);
        buff[size-2] = 0;
      }
    else
      {
        skipsize = temp * 10;
        ifs->seekg(skipsize, ios_base::cur);
        buff = new char[size-skipsize-1];
        ifs->read(buff, size-skipsize-2);
        buff[size-skipsize-2] = 0;
      }
    return buff;
  }


  int get2DPosition(istream *ifs, UINT32 size, INT32 &x, INT32 &y)
  {
    if(size != 8) {
      return -1; }

    READ_INT32 (*ifs, y);
    READ_INT32 (*ifs, x);

    // 	printf("2D coordinates - x: %d, y: %d\n", x, y);
    return 0;
  }

  UINT32 getBondStart(istream *ifs, UINT32 size)
  {
    UINT32 atomID;

    if(size != 4)
      return -1;

    READ_INT32 (*ifs, atomID);
    //	printf("Bond starts at atom nr: %08X\n", atomID);
    return atomID;
  }

  UINT32 getBondEnd(istream *ifs, UINT32 size)
  {
    UINT32 atomID;

    if(size != 4)
      return -1;

    READ_INT32 (*ifs, atomID);
    //	printf("Bond ends at atom nr: %08X\n", atomID);
    return atomID;
  }

  int getBondOrder(istream *ifs, UINT32 size)
  {
    UINT16 order;

    if(size != 2)
      return -1;

    READ_INT16 (*ifs, order);
    //	printf("Bond order is: %d\n", order);
    return (int) order;
  }

  int getBondDisplay(istream *ifs, UINT32 size)
  {
    UINT16 stereo;

    if(size != 2)
      return -1;

    READ_INT16 (*ifs, stereo);
    //    printf("Bond stereo is: %d\n", stereo);
    return (int) stereo;
  }

  int getNodeType(istream *ifs, UINT32 size)
  {
    UINT16 nodeType;

    if(size != 2)
      return -1;

    READ_INT16 (*ifs, nodeType);
    //	printf("Node type is: %d\n", nodeType);
    return (int) nodeType;
  }

  int getElement(istream *ifs, UINT32 size, OBAtom &atom)
  {
    UINT16 element;

    if(size != 2)
      return -1;

    READ_INT16 (*ifs, element);
    atom.SetAtomicNum(element);
    //	printf("Atomic number: %d\n", element);
    return 0;
  }

  int getIsotope(istream *ifs, UINT32 size, OBAtom &atom)
  {
    UINT16 isotope;

    if(size != 2)
      return -1;

    READ_INT16 (*ifs, isotope);
    atom.SetIsotope(isotope);
    return 0;
  }

  int getRadical(istream *ifs, UINT32 size, OBAtom &atom)
  {
    int radical;

    ifs->read((char *)&radical, size);
#if __BYTE_ORDER == __BIG_ENDIAN
    radical = radical >> 24;
#endif

    //    cout << " Atomic radical " << size << " " << radical << endl;
    switch (radical)
      {
      case 2:
        atom.SetSpinMultiplicity(2);
        break;
      case 3:
        atom.SetSpinMultiplicity(3);
        break;
      case 1:
      case 0:
      default:
        break; // all singlets (current no way in OB to set diradical singlet)
      }

    return 0;
  }

  int getCharge(istream *ifs, UINT32 size)
  {
    int charge;

    if(size == 4)		// Bug in ChemDraw 8.0, see http://www.cambridgesoft.com/services/documentation/sdk/chemdraw/cdx/properties/Atom_Charge.htm
      {
        READ_INT32 (*ifs, charge);
      }
    else
      if(size == 1)
        {
          ifs->read((char *)&charge, size);
#if __BYTE_ORDER == __BIG_ENDIAN
          charge = charge >> 24;
#endif
        }
      else
        return 0;

    //	printf("Charge: %d\n",charge);
    return charge;
  }

  int getAtomHydrogens(istream *ifs, UINT32 size)
  {
    UINT16 hydrogens;

    if(size != 2)
      return -1;

    READ_INT16 (*ifs, hydrogens);
    //	printf("Number of explicit hydrogens: %d\n", hydrogens);
    return 0;
  }

  int readText(istream *ifs, UINT32 textId)
  {
    char errorMsg[BUFF_SIZE];
    UINT16 tag;
    UINT16 size;
    UINT32 id;
    int depth = 1;

    while(ifs->good())
      {
        READ_INT16 (*ifs, tag);
        if(tag & kCDXTag_Object)	// Object
          {
            depth++;
            READ_INT32 (*ifs, id);
            //			printf("Object ID (in text-object %08X): %08X has type: %04X\n", textId, id, tag);
            snprintf(errorMsg, BUFF_SIZE, "New object in text, type %04X\n", tag);
            obErrorLog.ThrowError(__FUNCTION__, errorMsg, obDebug);
          }
        else if(tag == 0)	// End of object
          {
            depth--;
            //			printf("End of textobject %08X\n", textId);
          }
        else	// Property
          {
            READ_INT16 (*ifs, size);
            //			printf("Tag: %04X\tSize: %04X\n", tag, size);
            //switch(tag)
            //  {
              //default:
              ifs->seekg(size, ios_base::cur);
              //break;
            //  }
          }
        if(depth < 1)
          return 0;
      }
    return -1;
  }


  int ChemDrawBinaryFormat::readNode(istream *ifs, UINT32 nodeId,
                                     OBMol *pmol, map<UINT32, int> &atoms,
                                     list<cdBond> &bonds, UINT32 fragmentID)
  {
    char errorMsg[BUFF_SIZE];
    UINT16 tag;
    UINT16 size;
    UINT32 id;
    int depth = 1;
    OBAtom atom;
    OBAtom *atom2;
    //	OBPairData *data = new OBPairData;   // To hold ChemDraw's indexnr
    INT32 x, y;
    int nodeType = 1;
    //    char strNodeId[20];

    atom.SetAtomicNum(6);
    //	data->SetAttribute("nodeId");
    //	snprintf(strNodeId, 20, "%d", nodeId);
    //	data->SetValue(strNodeId);
    //	atom.SetData(data);
    //	if(atom.HasData("nodeId"))
    //	{
    //		printf("Adding data to atom - attr: %s\tvalue: %s\n", dynamic_cast<OBPairData *> (atom.GetData("nodeId"))->GetAttribute().c_str(), dynamic_cast<OBPairData *> (atom.GetData("nodeId"))->GetValue().c_str());
    //	}
    //	else
    //		printf("No data added\n");

    while(ifs->good())
      {
        READ_INT16 (*ifs, tag);
        if(tag & kCDXTag_Object)	// Object
          {
            depth++;
            READ_INT32 (*ifs, id);

            snprintf(errorMsg, BUFF_SIZE, "Object ID (in node %08X): %08X has type: %04X\n", nodeId, id, tag);
            obErrorLog.ThrowError(__FUNCTION__, errorMsg, obDebug);

            if(tag == kCDXObj_Fragment)
              {
                if(readFragment(ifs, id, pmol, atoms, bonds) != 0)
                  {
                    obErrorLog.ThrowError(__FUNCTION__, "Error reading fragment", obWarning);
                    return false;
                  }
                bonds.push_back(cdBond(nodeId, id, 1));
                depth--;
              }
            else if(tag == kCDXObj_Text)
              {
                readText(ifs, id);
                depth--;
              }
            else if((tag == kCDXObj_ObjectTag))
              {	// Objects to ignore
                readGeneric(ifs, id);
                depth--;
              }
            else if(tag == kCDXObj_Group)
              {	// Objects where we can ignore the header
              }
            else
              {
                snprintf(errorMsg, BUFF_SIZE, "New object in node, type %04X\n", tag);
                obErrorLog.ThrowError(__FUNCTION__, errorMsg, obDebug);
              }
          }
        else if(tag == 0)	// End of object
          {
            depth--;
            //			printf("End of Object, depth=%d\n", depth);
          }
        else	// Property
          {
            READ_INT16 (*ifs, size);
            snprintf(errorMsg, BUFF_SIZE, "Node Tag: %04X\tSize: %04X\n", tag, size);
            obErrorLog.ThrowError(__FUNCTION__, errorMsg, obDebug);

            switch(tag)
              {
              case kCDXProp_Atom_NumHydrogens:
                getAtomHydrogens(ifs, size);
                break;
              case kCDXProp_2DPosition:
                get2DPosition(ifs, size, x, y);
                atom.SetVector((double) x / 500000., (double) y / -500000.0, (double) 0.0);
                break;
              case kCDXProp_Node_Element:
                getElement(ifs, size, atom);
                break;
              case kCDXProp_Atom_Charge:
                atom.SetFormalCharge(getCharge(ifs, size));
                break;
              case kCDXProp_Node_Type:
                nodeType = getNodeType(ifs, size);
                break;
              case kCDXProp_Atom_Isotope:
                getIsotope(ifs, size, atom);
                break;
              case kCDXProp_Atom_Radical:
                getRadical(ifs, size, atom);
                break;
              case kCDXProp_Atom_ElementList:
              case kCDXProp_Atom_AbnormalValence:
              case kCDXProp_Name:
              case kCDXProp_IgnoreWarnings:
              case kCDXProp_ChemicalWarning:
                ifs->seekg(size, ios_base::cur);
                break;  // Maybe use? (any of the above properties)
              case kCDXProp_MarginWidth:
              case kCDXProp_ZOrder:
              case kCDXProp_Atom_GenericNickname:
              case kCDXProp_Atom_CIPStereochemistry: // maybe use for chirality
              case kCDXProp_Atom_Geometry:
              case kCDXProp_Atom_BondOrdering:
              case kCDXProp_Node_LabelDisplay:
              case kCDXProp_LabelStyle:
              case kCDXProp_ForegroundColor:
              case kCDXProp_BackgroundColor:
                ifs->seekg(size, ios_base::cur);
                break;
              default: // some unknown node tag
                ifs->seekg(size, ios_base::cur);
                snprintf(errorMsg, BUFF_SIZE, "Node Tag: %04X\tSize: %04X\n", tag, size);
                obErrorLog.ThrowError(__FUNCTION__, errorMsg, obDebug);
                break;
              }
          }
        if(depth < 1)
          {
						// Check the data before returning...
            switch(nodeType)
              {
              case kCDXNodeType_Fragment:
                atoms[nodeId] = -1;
                break;
              case kCDXNodeType_ExternalConnectionPoint:
                atoms[nodeId] = -1;
                bonds.push_back(cdBond(nodeId,fragmentID,1));
                break;
              default:
                atom2 = pmol->NewAtom();
                atom2->SetVector(atom.GetVector());
                atom2->SetFormalCharge(atom.GetFormalCharge());
                atom2->SetAtomicNum(atom.GetAtomicNum());
                /*
                if (atom.IsClockwise())
                  atom2->SetClockwiseStereo();
                else if (atom.IsAntiClockwise())
                  atom2->SetAntiClockwiseStereo();
                */
                atoms[nodeId] = atom2->GetIdx();
                break;
              }
            return 0; // everything worked correctly
          }
      } // end while(ifs.good())
    return -1; // file ended early
  }

  int ChemDrawBinaryFormat::readBond(istream *ifs, UINT32 bondId,
                                     OBMol *pmol, list<cdBond> &bonds)
  {
    char errorMsg[BUFF_SIZE];
    UINT16 tag;
    UINT16 size;
    UINT32 id, bgnID, endID;
    int depth = 1, order=1, stereo = 0;

    while(ifs->good())
      {
        READ_INT16 (*ifs, tag);
        if(tag & kCDXTag_Object)	// Object
          {
            depth++;
            READ_INT32 (*ifs, id);
            snprintf(errorMsg, BUFF_SIZE, "Object ID (in bond %08X): %08X has type: %04X\n", bondId, id, tag);
            obErrorLog.ThrowError(__FUNCTION__, errorMsg, obDebug);

            if(tag == kCDXObj_Text)
              {
                readText(ifs, id);
                depth--;
              }
            else
              {
                snprintf(errorMsg, BUFF_SIZE, "New object in bond, type %04X\n", tag);
                obErrorLog.ThrowError(__FUNCTION__, errorMsg, obDebug);
              }
          }
        else if(tag == 0)	// End of object
          {
            depth--;
            //			printf("End of Object\n");
          }
        else	// Property
          {
            READ_INT16 (*ifs, size);

            snprintf(errorMsg, BUFF_SIZE, "Bond Tag: %04X\tSize: %04X\n", tag, size);
            obErrorLog.ThrowError(__FUNCTION__, errorMsg, obDebug);

            switch(tag)
              {
              case kCDXProp_Bond_Begin:
                bgnID = getBondStart(ifs, size);
                break;
              case kCDXProp_Bond_End:
                endID = getBondEnd(ifs, size);
                break;
              case kCDXProp_Bond_Order:
                order = getBondOrder(ifs, size);
                switch (order)
                  {
                  case 0xFFFF: // undefined, keep 1 for now
                    order = 1;
                  case 0x0001:
                  case 0x0002:
                    break;
                  case 0x0004:
                    order = 3;
                    break;
                  case 0x0080: // aromatic bond
                    order = 5;
                    break;
                  default: // other cases are just not supported, keep 1
                    order = 1;
                    break;
                  }
                break;
              case kCDXProp_Bond_Display:
                stereo = getBondDisplay(ifs, size);
                break; // Stereo info

              case kCDXProp_BondLength:
              case kCDXProp_Bond_CIPStereochemistry:
              case kCDXProp_Bond_BeginAttach:
              case kCDXProp_Bond_EndAttach:
              case kCDXProp_ChemicalWarning:
                ifs->seekg(size, ios_base::cur);
                break; // Maybe use? (any of the above cases)

              case kCDXProp_IgnoreWarnings:
              case kCDXProp_MarginWidth:
              case kCDXProp_Bond_Display2:
              case kCDXProp_Bond_DoublePosition:
              case kCDXProp_Bond_BondOrdering:
              case kCDXProp_BondSpacing:
              case kCDXProp_LineWidth:
              case kCDXProp_BoldWidth:
              case kCDXProp_HashSpacing:
              case kCDXProp_ZOrder:
              case kCDXProp_ForegroundColor:
              case kCDXProp_BackgroundColor:
              case kCDXProp_LabelStyle:
                ifs->seekg(size, ios_base::cur);
                break;

              default: // some unknown, undocumented property
                ifs->seekg(size, ios_base::cur);
                snprintf(errorMsg, BUFF_SIZE, "Bond Tag: %04X\tSize: %04X\n", tag, size);
                obErrorLog.ThrowError(__FUNCTION__, errorMsg, obDebug);
                break;
              }
          }
        if(depth < 1)
          {
            bonds.push_back(cdBond(bgnID,endID,order,stereo));
            return 0;
          }
      }
    return -1;
  }

  int ChemDrawBinaryFormat::readGeneric(istream *ifs, UINT32 objId)
  {
    char errorMsg[BUFF_SIZE];
    UINT16 tag;
    UINT16 size;
    UINT32 id;
    int depth = 1;

    while(ifs->good())
      {
        READ_INT16 (*ifs, tag);
        if(tag & kCDXTag_Object)	// Object
          {
            depth++;
            READ_INT32 (*ifs, id);

            snprintf(errorMsg, BUFF_SIZE, "Object ID (in generic %08X): %08X has type: %04X\n", objId, id, tag);
            obErrorLog.ThrowError(__FUNCTION__, errorMsg, obDebug);

            if ((tag == kCDXObj_BracketAttachment)
                || (tag == kCDXObj_CrossingBond)
                || (tag == kCDXObj_BracketedGroup)
                || (tag == kCDXObj_Text)
                || (tag == kCDXObj_Fragment))
              {	// Objects to ignore (Fragment might be worth parsing)
                readGeneric(ifs, id);
                depth--;
              }
            else {
              snprintf(errorMsg, BUFF_SIZE, "New object in generic, type %04X\n", tag);
              obErrorLog.ThrowError(__FUNCTION__, errorMsg, obDebug);
            }
          }
        else if(tag == 0)	// End of object
          {
            depth--;
            snprintf(errorMsg, BUFF_SIZE, "End of Object in generic %08X\n", objId);
            obErrorLog.ThrowError(__FUNCTION__, errorMsg, obDebug);
          }
        else	// Property
          {
            READ_INT16 (*ifs, size);

            snprintf(errorMsg, BUFF_SIZE, "Generic Tag: %04X\tSize: %04X\n", tag, size);
            obErrorLog.ThrowError(__FUNCTION__, errorMsg, obDebug);
            //switch(tag)
            //  {
            //  default:
                ifs->seekg(size, ios_base::cur);
            //   break;
            //  }
          }
        if(depth < 1) // matched begin and end tags, so we're good
          {
            return 0;
          }
      } // while reading
    return -1; // file ended while reading
  }

  int ChemDrawBinaryFormat::readFragment(istream *ifs, UINT32 fragmentId,
                                         OBMol *pmol, map<UINT32, int> &atoms,
                                         list<cdBond> &bonds)
  {
    char errorMsg[BUFF_SIZE];
    UINT16 tag;
    UINT16 size;
    UINT32 id;
    //char u32[4];
    int depth = 1;

 cerr<<"Reading "<<pmol<<endl;
   atoms[fragmentId] = -1;
    while(ifs->good())
      {
        READ_INT16 ((*ifs), tag);
        if(tag & kCDXTag_Object)	// Object
          {
            depth++;
            READ_INT32 (*ifs, id);

            snprintf(errorMsg, BUFF_SIZE, "Object ID (in fragment %08X): %08X has type: %04X\n", fragmentId, id, tag);
            obErrorLog.ThrowError(__FUNCTION__, errorMsg, obDebug);

            if(tag == kCDXObj_Fragment)
              {
                if(readFragment(ifs, id, pmol, atoms, bonds) != 0)
                  {
                    obErrorLog.ThrowError(__FUNCTION__, "Error reading fragment", obWarning);
                    return false;
                  }
              }
            else if(tag == kCDXObj_Node)
              {
                readNode(ifs, id, pmol, atoms, bonds, fragmentId);
                depth--;
              }
            else if(tag == kCDXObj_Bond)
              {
                readBond(ifs, id, pmol, bonds);
                depth--;
              }
            else if((tag == kCDXObj_Graphic) || (tag == kCDXObj_Text))
              {	// Objects that can be ignored
                readGeneric(ifs, id);
                depth--;
              }
            //			else if((tag == kCDXObj_Page) || ( tag == kCDXObj_Group))
            //			{	// Objects where the header can be ignored
            //			}
            else
              {
                snprintf(errorMsg, BUFF_SIZE, "New object in fragment, type %04X\n", tag);
                obErrorLog.ThrowError(__FUNCTION__, errorMsg, obDebug);
              }

          }
        else if(tag == 0)	// End of object
          {
            depth--;
          }
        else	// Property
          {
            READ_INT16 ((*ifs), size);
            switch(tag)
              {
              case kCDXProp_Frag_ConnectionOrder:
                ifs->seekg(size, ios_base::cur);
                break;  // Should use
              case kCDXProp_BoundingBox:
                ifs->seekg(size, ios_base::cur);
                break;
              default:
                ifs->seekg(size, ios_base::cur);
                snprintf(errorMsg, BUFF_SIZE, "Fragment Tag: %04X\tSize: %04X\n", tag, size);
                obErrorLog.ThrowError(__FUNCTION__, errorMsg, obDebug);

                break;
              }
          }
        if(depth < 1)
          {
cerr<<"Done reading "<<pmol<<endl;
            return 0; // all begin and end tags matched -- everything good
          }
      } // while reading
    return -1; // file ended before reading was finished
  }

} //namespace OpenBabel

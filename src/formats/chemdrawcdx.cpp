/**********************************************************************
Copyright (C) 2006 by Fredrik Wallner
Some portions Copyright (C) 2006 by Geoffrey Hutchsion
Some portions Copyright (C) 2004 by Chris Morley
 
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include "babelconfig.h"
#include "obmolecformat.h"
#include "chemdrawcdx.h"

#include <iostream>
#include <fstream>
#include <map>
#include <list>

//#define debug

using namespace std;
namespace OpenBabel
{

  class ChemDrawBinaryFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID in the constructor
    ChemDrawBinaryFormat()
    {
      OBConversion::RegisterFormat("cdx",this, "chemical-x-cdx");
    }

    virtual const char* Description() //required
    {
      return
        "ChemDraw binary format\n \
		Read only.\n";
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


    int readFragment(istream *ifs, UINT32 fragmentId, OBMol *pmol, map<UINT32, int> &atoms, list<cdBond> &bonds);
    int readNode(istream *ifs, UINT32 nodeId, OBMol *pmol, map<UINT32, int> &atoms, list<cdBond> &bonds, UINT32 fragmentID);
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
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;

    istream& ifs = *pConv->GetInStream();
    char buffer[BUFF_SIZE];
    UINT16 tag;
    UINT16 size;
    UINT32 id;
    char u32[4];
    map<UINT32, int> atoms;
    list<cdBond> bonds;
    list<cdBond>::const_iterator bondsIter;
    cdBond prevBond;
    cdBond oneBond;
    OBPairData *pd;
    int iInt=0, depth=1;


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
        ifs.read((char *)&tag, sizeof(tag));
        tag = (tag << 8) | (tag >> 8);
        if(tag & kCDXTag_Object)	// Object
          {
            ifs.read(u32, sizeof(u32));
            id = ((long) (u32[3] & 0xff) << 24) | ((long) (u32[2] & 0xff) << 16) | ((long) (u32[1] & 0xff) << 8) | ((long) (u32[0] & 0xff));
#ifdef debug
            printf("Object ID: %08X in root has type: %04X\n", id, tag);
#endif
            if(tag == kCDXObj_Fragment)
              {
                if(readFragment(&ifs, id, pmol, atoms, bonds) != 0)
                  {
                    printf("Error reading fragment\n");
                    return false;
                  }
                //				printf("Number of atoms in molecule: %d\n", pmol->NumNodes());
                iInt = 0;
                /*				FOR_ATOMS_OF_MOL(a, pmol)
                          {
                          printf("Atom #%d\tatomic number #%d\t", ++iInt, a->GetAtomicNum());
                          if(a->HasData(OBGenericDataType::PairData))
                          {
                          pd = dynamic_cast<OBPairData *> (a->GetData("nodeId"));
                          printf("%s #%s\n", dynamic_cast<OBPairData *> (a->GetData("nodeId"))->GetAttribute().c_str(), dynamic_cast<OBPairData *> (a->GetData("nodeId"))->GetValue().c_str());
                          atoms[atol(pd->GetValue().c_str())] = a->GetIdx();
                          }
                          else
                          printf("no nodeId!\n");
                          }
                */
                //				printf("Number of bonds: %d\n", bonds.size());
                //				for (bondsIter=bonds.begin(); bondsIter != bonds.end(); bondsIter++)
                //				{
                //					printf("Bond between %08X (%d) and %08x (%d)\n", bondsIter->begin, atoms[bondsIter->begin], bondsIter->end, atoms[bondsIter->end]);
                //				}
                prevBond = cdBond(0,0,0);
                for (bondsIter=bonds.begin(); bondsIter != bonds.end(); bondsIter++)
                  {
                    //					printf("Bond between %08X (%d) and %08x (%d)\n", bondsIter->begin, atoms[bondsIter->begin], bondsIter->end, atoms[bondsIter->end]);
                    if(atoms[bondsIter->begin] == -1)
                      {
                        //						printf("Bond starts at a non-atom\n");
                        if(atoms[bondsIter->end] == -1)
                          {
                            //							printf("Bond between two non-atoms pushed back\n");
                            if((prevBond.begin == bondsIter->begin) && (prevBond.end == bondsIter->end))
                              {
                                //								printf("Can't assign all bonds!\n");
                                break;
                              }	
                            bonds.push_back(*bondsIter);
                          }
                        else
                          {
                            //							printf("Renumber bond start to real atom %d...", atoms[bondsIter->end]);
                            atoms[bondsIter->begin] = atoms[bondsIter->end];
                            //							printf(" Done\n");
                          }
                      }
                    else if(atoms[bondsIter->end] == -1)
                      {
                        //						printf("Renumber bond end to real atom\n");						
                        atoms[bondsIter->end] = atoms[bondsIter->begin];
                      }
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
            else if((tag == kCDXObj_Graphic) || (tag == kCDXObj_Text) || (tag == kCDXObj_BracketedGroup) || (tag == kCDXObj_BracketAttachment) || (tag == kCDXObj_CrossingBond) || (tag == kCDXObj_ReactionStep) || (tag == kCDXObj_Curve) || (tag == kCDXObj_EmbeddedObject)) 
              {	// Objects that can be ignored
                readGeneric(&ifs, id);
              }
            else if((tag == kCDXObj_Page) || ( tag == kCDXObj_Group))
              {	// Objects where the header can be ignored
              }
            else
              {
                depth++;
                printf("New object in root, type %04X\n", tag);
              }
          }
        else if(tag == 0)	// End of object
          {
            if(depth > 1)
              {
                printf("End object\n");
                depth--;
              }
          }
        else	// Property
          {
            ifs.read((char *)&size, sizeof(size));
            size = (size << 8) | (size >> 8);
            //			printf("Root Tag: %04X\tSize: %04X\n", tag, size);
            switch(tag)
              {
              case kCDXProp_Name: pmol->SetTitle(getName(&ifs, size)); break;
              case kCDXProp_FontTable: ifs.seekg(size, ios_base::cur); break;
              case kCDXProp_BoundingBox: ifs.seekg(size, ios_base::cur); break;
              case kCDXProp_Window_Position: ifs.seekg(size, ios_base::cur); break;
              case kCDXProp_Window_Size: ifs.seekg(size, ios_base::cur); break;
              case kCDXProp_Atom_ShowQuery: ifs.seekg(size, ios_base::cur); break;
              case kCDXProp_Atom_ShowStereo: ifs.seekg(size, ios_base::cur); break;
              case kCDXProp_Atom_ShowAtomNumber: ifs.seekg(size, ios_base::cur); break;
              case kCDXProp_Bond_ShowQuery: ifs.seekg(size, ios_base::cur); break;
              case kCDXProp_Bond_ShowStereo: ifs.seekg(size, ios_base::cur); break;
              case kCDXProp_MacPrintInfo: ifs.seekg(size, ios_base::cur); break;
              case kCDXProp_DrawingSpaceType: ifs.seekg(size, ios_base::cur); break;
              case kCDXProp_Width: ifs.seekg(size, ios_base::cur); break;
              case kCDXProp_Height: ifs.seekg(size, ios_base::cur); break;
              case kCDXProp_PageOverlap: ifs.seekg(size, ios_base::cur); break;
              case kCDXProp_Magnification: ifs.seekg(size, ios_base::cur); break;
              case kCDXProp_ChainAngle: ifs.seekg(size, ios_base::cur); break;
              case kCDXProp_BondSpacing: ifs.seekg(size, ios_base::cur); break;
              case kCDXProp_BondSpacingAbs: ifs.seekg(size, ios_base::cur); break;
              case kCDXProp_BondLength: ifs.seekg(size, ios_base::cur); break;
              case kCDXProp_BoldWidth: ifs.seekg(size, ios_base::cur); break;
              case kCDXProp_LineWidth: ifs.seekg(size, ios_base::cur); break;
              case kCDXProp_MarginWidth: ifs.seekg(size, ios_base::cur); break;
              case kCDXProp_HashSpacing: ifs.seekg(size, ios_base::cur); break;
              case kCDXProp_LabelStyle: ifs.seekg(size, ios_base::cur); break;
              case kCDXProp_CaptionStyle: ifs.seekg(size, ios_base::cur); break;
              case kCDXProp_CaptionJustification: ifs.seekg(size, ios_base::cur); break;
              case kCDXProp_LabelJustification: ifs.seekg(size, ios_base::cur); break;
              case kCDXProp_FractionalWidths: ifs.seekg(size, ios_base::cur); break;
              case kCDXProp_LabelLineHeight: ifs.seekg(size, ios_base::cur); break;
              case kCDXProp_CaptionLineHeight: ifs.seekg(size, ios_base::cur); break;
              case kCDXProp_Window_IsZoomed: ifs.seekg(size, ios_base::cur); break;
              case kCDXProp_WidthPages: ifs.seekg(size, ios_base::cur); break;
              case kCDXProp_HeightPages: ifs.seekg(size, ios_base::cur); break;
              case kCDXProp_HeaderPosition: ifs.seekg(size, ios_base::cur); break;
              case kCDXProp_FooterPosition: ifs.seekg(size, ios_base::cur); break;
              case kCDXProp_PrintTrimMarks: ifs.seekg(size, ios_base::cur); break;
              case kCDXProp_PrintMargins: ifs.seekg(size, ios_base::cur); break;
              case kCDXProp_ColorTable: ifs.seekg(size, ios_base::cur); break;
              case kCDXProp_CreationProgram: ifs.seekg(size, ios_base::cur); break;
              case kCDXProp_Group_Integral: ifs.seekg(size, ios_base::cur); break;
              default: ifs.seekg(size, ios_base::cur); printf("Root Tag: %04X\tSize: %04X\n", tag, size); break;
              }
          }
      }
    //	if(atoms.empty())
    //		printf("Atoms is empty\n");
    //	else
    //		printf("Atoms holds %d entries\n", atoms.size());
    //	for (map<UINT32, int>::iterator im = atoms.begin(); im != atoms.end(); ++im )
    //        printf("\"%08X\"  = %d\n", im->first, im->second);
    //
    //    OBAtom *atom;
    //    vector<OBAtom*>::iterator i;
    //
    //    for(atom = pmol->BeginAtom(i);atom;atom = pmol->NextAtom(i))
    //    {
    //        printf("Atom %d: %9.4f %9.4f    0.0000 %-1s\n",
    //        		atom->GetIdx(),
    //                atom->x(),
    //                atom->y(),
    //                etab.GetSymbol(atom->GetAtomicNum()));
    //    }
    //
    //    OBBond *bond;
    //    vector<OBBond*>::iterator j;
    //
    //    for(bond = pmol->BeginBond(j);bond;bond = pmol->NextBond(j))
    //    {
    //        printf("Bond %d: %3d%3d%3d%3d\n",
    //        		bond->GetIdx(),
    //                bond->GetBeginAtomIdx(),
    //                bond->GetEndAtomIdx(),
    //                bond->GetBO(), bond->GetBO());
    //    }
    //  	

    /** Parse the input stream and use the OpenBabel API to populate the OBMol **/

    // To use an input option
    /*	if(pConv->IsOption("s",OBConversion::INOPTIONS))
        {
        //Code for when -as is specified
        }
    */
    /* If the molecule has other than 3D coordinates for its atoms, it
       is necessary to set the dimension to 0, or 2 */
    pmol->SetDimension(2);

    pmol->EndModify();
    pmol->Center();

    /* For multi-molecule formats, leave the input stream at the start of the
       next molecule, ready for this routine to be called again. 

       /* Return true if ok. Returning false means discard the OBMol and stop
       converting, unless the -e option is set. With a multi-molecule inputstream
       this will skip the current molecule and continue with the next, if SkipObjects()
       has been defined. If it has not, and continuation after errors is still required,
       it is necessary to leave the input stream at the beginning of next object when
       returning false;*/
    return true;
  }

  const char* ChemDrawBinaryFormat::getName(istream *ifs, UINT32 size)
  {
    char *buff;
    UINT16 temp;
    UINT32 skipsize;
		
    //	buff = new char[size];
    ifs->read((char *)&temp, 2);
    temp = (temp << 8) | (temp >> 8);
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
    char *pos;
    //	INT32 x, y;
	
    if(size != 8) {
      return -1; }
	
    pos = new char[size];
    ifs->read(pos, size);
    y = ((long) (pos[3] & 0xff) << 24) | ((long) (pos[2] & 0xff) << 16) | ((long) (pos[1] & 0xff) << 8) | ((long) (pos[0] & 0xff));
    x = ((long) (pos[7] & 0xff) << 24) | ((long) (pos[6] & 0xff) << 16) | ((long) (pos[5] & 0xff) << 8) | ((long) (pos[4] & 0xff));
    
    //	printf("2D coordinates - x: %d, y: %d\n", x, y);
    delete pos;
    return 0;
  }

  UINT32 getBondStart(istream *ifs, UINT32 size)
  {
    char *read;
    UINT32 atomID;
	
    if(size != 4)
      return -1;
	
    read = new char[size];
    ifs->read(read, size);
    atomID = (UINT32)((read[3] & 0xff) << 24) | ((long) (read[2] & 0xff) << 16) | ((long) (read[1] & 0xff) << 8) | ((long) (read[0] & 0xff));
    //	printf("Bond starts at atom nr: %08X\n", atomID);
    delete read;
    return atomID;
  }

  UINT32 getBondEnd(istream *ifs, UINT32 size)
  {
    char *read;
    UINT32 atomID;
	
    if(size != 4)
      return -1;
	
    read = new char[size];
    ifs->read(read, size);
    atomID = (UINT32)((read[3] & 0xff) << 24) | ((long) (read[2] & 0xff) << 16) | ((long) (read[1] & 0xff) << 8) | ((long) (read[0] & 0xff));
    //	printf("Bond ends at atom nr: %08X\n", atomID);
    delete read;
    return atomID;
  }

  int getBondOrder(istream *ifs, UINT32 size)
  {
    UINT16 order;
	
    if(size != 2)
      return -1;
	
    ifs->read((char *)&order, size);
    order = (order << 8) | (order >> 8);
    //	printf("Bond order is: %d\n", order);
    return (int) order;
  }

  int getBondDisplay(istream *ifs, UINT32 size)
  {
    UINT16 stereo;
	
    if(size != 2)
      return -1;
	
    ifs->read((char *)&stereo, size);
    stereo = (stereo << 8) | (stereo >> 8);
    //    printf("Bond stereo is: %d\n", stereo);
    return (int) stereo;
  }

  int getNodeType(istream *ifs, UINT32 size)
  {
    UINT16 nodeType;
	
    if(size != 2)
      return -1;
	
    ifs->read((char *)&nodeType, size);
    nodeType = (nodeType << 8) | (nodeType >> 8);
    //	printf("Node type is: %d\n", nodeType);
    return (int) nodeType;
  }

  int getElement(istream *ifs, UINT32 size, OBAtom &atom)
  {
    UINT16 element;
	
    if(size != 2)
      return -1;
	
    ifs->read((char *)&element, size);
    element = (element << 8) | (element >> 8);
    atom.SetAtomicNum(element);
    //	printf("Atomic number: %d\n", element);
    return 0;
  }

  int getCharge(istream *ifs, UINT32 size)
  {
    int charge;
	
    if(size == 4)		// Bug in ChemDraw 8.0, see http://www.cambridgesoft.com/services/documentation/sdk/chemdraw/cdx/properties/Atom_Charge.htm
      {}
    else
      if(size != 1)
        return 0;
	
    ifs->read((char *)&charge, size);
    charge = charge >> 24;
    //	printf("Charge: %d\n",charge);
    return charge;
  }

  int getAtomHydrogens(istream *ifs, UINT32 size)
  {
    UINT16 hydrogens;
	
    if(size != 2)
      return -1;
	
    ifs->read((char *)&hydrogens, size);
    hydrogens = (hydrogens << 8) | (hydrogens >> 8);
    //	printf("Number of explicit hydrogens: %d\n", hydrogens);
    return 0;
  }

  int readText(istream *ifs, UINT32 textId)
  {
    //	char buffer[1024];
    UINT16 tag;
    UINT16 size;
    UINT32 id;
    char u32[4];
    int depth = 1;

    while(ifs->good())
      {
        ifs->read((char *)&tag, sizeof(tag));
        tag = (tag << 8) | (tag >> 8);
        if(tag & kCDXTag_Object)	// Object
          {
            depth++;
            ifs->read(u32, sizeof(u32));
            id = ((long) (u32[3] & 0xff) << 24) | ((long) (u32[2] & 0xff) << 16) | ((long) (u32[1] & 0xff) << 8) | ((long) (u32[0] & 0xff));
            //			printf("Object ID (in text-object %08X): %08X has type: %04X\n", textId, id, tag);
            printf("New object in text, type %04X\n", tag);
          }
        else if(tag == 0)	// End of object
          {
            depth--;
            //			printf("End of textobject %08X\n", textId);
          }
        else	// Property
          {
            ifs->read((char *)&size, sizeof(size));
            size = (size << 8) | (size >> 8);
            //			printf("Tag: %04X\tSize: %04X\n", tag, size);
            switch(tag)
              {
              default: ifs->seekg(size, ios_base::cur); break;
              }
          }
        if(depth < 1)
          return 0;
      }
    return -1;
  }

		
  int ChemDrawBinaryFormat::readNode(istream *ifs, UINT32 nodeId, OBMol *pmol, map<UINT32, int> &atoms, list<cdBond> &bonds, UINT32 fragmentID)
  {
    //	char buffer[1024];
    UINT16 tag;
    UINT16 size;
    UINT32 id;
    char u32[4];
    int depth = 1;
    OBAtom atom;
    OBAtom *atom2;
    //	OBPairData *data = new OBPairData;   // To hold ChemDraw's indexnr
    INT32 x, y;
    int nodeType = 1;
    char strNodeId[20];

    atom.SetAtomicNum(6);
    //	data->SetAttribute("nodeId");
    //	sprintf(strNodeId, "%d", nodeId);
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
        ifs->read((char *)&tag, sizeof(tag));
        tag = (tag << 8) | (tag >> 8);
        if(tag & kCDXTag_Object)	// Object
          {
            depth++;
            ifs->read(u32, sizeof(u32));
            id = ((long) (u32[3] & 0xff) << 24) | ((long) (u32[2] & 0xff) << 16) | ((long) (u32[1] & 0xff) << 8) | ((long) (u32[0] & 0xff));
#ifdef debug
            printf("Object ID (in node %08X): %08X has type: %04X\n", nodeId, id, tag);
#endif
            if(tag == kCDXObj_Fragment)
              {
                if(readFragment(ifs, id, pmol, atoms, bonds) != 0)
                  {
                    printf("Error reading fragment\n");
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
                printf("New object in node, type %04X\n", tag);
              }			
          }
        else if(tag == 0)	// End of object
          {
            depth--;
            //			printf("End of Object, depth=%d\n", depth);
          }
        else	// Property
          {
            ifs->read((char *)&size, sizeof(size));
            size = (size << 8) | (size >> 8);
#ifdef debug
            printf("Node Tag: %04X\tSize: %04X\n", tag, size);
#endif
            switch(tag)
              {
              case kCDXProp_Atom_NumHydrogens: getAtomHydrogens(ifs, size); break;
              case kCDXProp_2DPosition: get2DPosition(ifs, size, x, y); 
                atom.SetVector((double) x / 500000.0, (double) y / -500000.0, (double) 0.0); break;
              case kCDXProp_Node_Element: getElement(ifs, size, atom); break;
              case kCDXProp_Atom_Charge: atom.SetFormalCharge(getCharge(ifs, size)); break;
              case kCDXProp_Node_Type: nodeType = getNodeType(ifs, size); break;
              case kCDXProp_Atom_Isotope: ifs->seekg(size, ios_base::cur); break;  // Should use
              case kCDXProp_Atom_ElementList: ifs->seekg(size, ios_base::cur); break;  // Should use
              case kCDXProp_Atom_Radical: ifs->seekg(size, ios_base::cur); break;  // Should use
              case kCDXProp_Atom_AbnormalValence: ifs->seekg(size, ios_base::cur); break;  
              case kCDXProp_Name: ifs->seekg(size, ios_base::cur); break;  
              case kCDXProp_IgnoreWarnings: ifs->seekg(size, ios_base::cur); break; 
              case kCDXProp_ChemicalWarning: ifs->seekg(size, ios_base::cur); break;  // Maybe use?
              case kCDXProp_MarginWidth: ifs->seekg(size, ios_base::cur); break;
              case kCDXProp_ZOrder: ifs->seekg(size, ios_base::cur); break;
              case kCDXProp_Atom_GenericNickname: ifs->seekg(size, ios_base::cur); break;
              case kCDXProp_Atom_CIPStereochemistry: ifs->seekg(size, ios_base::cur); break;
              case kCDXProp_Atom_Geometry: ifs->seekg(size, ios_base::cur); break;
              case kCDXProp_Atom_BondOrdering: ifs->seekg(size, ios_base::cur); break;
              case kCDXProp_Node_LabelDisplay: ifs->seekg(size, ios_base::cur); break;
              case kCDXProp_LabelStyle: ifs->seekg(size, ios_base::cur); break;
              case kCDXProp_ForegroundColor: ifs->seekg(size, ios_base::cur); break;
              case kCDXProp_BackgroundColor: ifs->seekg(size, ios_base::cur); break;
              default: ifs->seekg(size, ios_base::cur); printf("Node Tag: %04X\tSize: %04X\n", tag, size); break;
              }
          }
        if(depth < 1)
          {
						// Check the data before returning...
            switch(nodeType)
              {
              case kCDXNodeType_Fragment: atoms[nodeId] = -1; break;
              case kCDXNodeType_ExternalConnectionPoint: atoms[nodeId] = -1; bonds.push_back(cdBond(nodeId,fragmentID,1)); break;
              default: atom2 = pmol->NewAtom(); atom2->SetVector(atom.GetVector()); atom2->SetFormalCharge(atom.GetFormalCharge()); atom2->SetAtomicNum(atom.GetAtomicNum()); atoms[nodeId] = atom2->GetIdx(); /*printf("nodeId: %08X\tatomIdx: %d\n", nodeId, atom2->GetIdx());*/  break;
              }
            return 0;
          }
      }
    return -1;
  }

  int ChemDrawBinaryFormat::readBond(istream *ifs, UINT32 bondId, OBMol *pmol, list<cdBond> &bonds)
  {
    //	char buffer[1024];
    UINT16 tag;
    UINT16 size;
    UINT32 id, bgnID, endID;
    char u32[4];
    int depth = 1, order=1, stereo = 0;

    while(ifs->good())
      {
        ifs->read((char *)&tag, sizeof(tag));
        tag = (tag << 8) | (tag >> 8);
        if(tag & kCDXTag_Object)	// Object
          {
            depth++;
            ifs->read(u32, sizeof(u32));
            id = ((long) (u32[3] & 0xff) << 24) | ((long) (u32[2] & 0xff) << 16) | ((long) (u32[1] & 0xff) << 8) | ((long) (u32[0] & 0xff));
#ifdef debug
            printf("Object ID (in bond %08X): %08X has type: %04X\n", bondId, id, tag);
#endif
            if(tag == kCDXObj_Text)
              {
                readText(ifs, id);
                depth--;
              }
            else
              {
                printf("New object in bond, type %04X\n", tag);
              }
          }
        else if(tag == 0)	// End of object
          {
            depth--;
            //			printf("End of Object\n");
          }
        else	// Property
          {
            ifs->read((char *)&size, sizeof(size));
            size = (size << 8) | (size >> 8);
#ifdef debug
            printf("Bond Tag: %04X\tSize: %04X\n", tag, size);
#endif
            switch(tag)
              {
              case kCDXProp_Bond_Begin: bgnID = getBondStart(ifs, size); break;
              case kCDXProp_Bond_End: endID = getBondEnd(ifs, size); break;
              case kCDXProp_Bond_Order: order = getBondOrder(ifs, size); break;
              case kCDXProp_Bond_Display: stereo = getBondDisplay(ifs, size); break; // Stereo info
              case kCDXProp_Bond_CIPStereochemistry: ifs->seekg(size, ios_base::cur); break; // Maybe use?
              case kCDXProp_Bond_BeginAttach: ifs->seekg(size, ios_base::cur); break; // Maybe use?
              case kCDXProp_Bond_EndAttach: ifs->seekg(size, ios_base::cur); break; // Maybe use?
              case kCDXProp_ChemicalWarning: ifs->seekg(size, ios_base::cur); break; // Maybe use?
              case kCDXProp_IgnoreWarnings: ifs->seekg(size, ios_base::cur); break;
              case kCDXProp_MarginWidth: ifs->seekg(size, ios_base::cur); break;
              case kCDXProp_Bond_Display2: ifs->seekg(size, ios_base::cur); break;
              case kCDXProp_Bond_DoublePosition: ifs->seekg(size, ios_base::cur); break;
              case kCDXProp_Bond_BondOrdering: ifs->seekg(size, ios_base::cur); break;
              case kCDXProp_BondLength: ifs->seekg(size, ios_base::cur); break;
              case kCDXProp_BondSpacing: ifs->seekg(size, ios_base::cur); break;
              case kCDXProp_LineWidth: ifs->seekg(size, ios_base::cur); break;
              case kCDXProp_BoldWidth: ifs->seekg(size, ios_base::cur); break;
              case kCDXProp_HashSpacing: ifs->seekg(size, ios_base::cur); break;
              case kCDXProp_ZOrder: ifs->seekg(size, ios_base::cur); break;
              case kCDXProp_ForegroundColor: ifs->seekg(size, ios_base::cur); break;
              case kCDXProp_BackgroundColor: ifs->seekg(size, ios_base::cur); break;
              case kCDXProp_LabelStyle: ifs->seekg(size, ios_base::cur); break;
              default: ifs->seekg(size, ios_base::cur); printf("Bond Tag: %04X\tSize: %04X\n", tag, size); break;
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
    UINT16 tag;
    UINT16 size;
    UINT32 id;
    char u32[4];
    int depth = 1;

    while(ifs->good())
      {
        ifs->read((char *)&tag, sizeof(tag));
        tag = (tag << 8) | (tag >> 8);
        if(tag & kCDXTag_Object)	// Object
          {
            depth++;
            ifs->read(u32, sizeof(u32));
            id = ((long) (u32[3] & 0xff) << 24) | ((long) (u32[2] & 0xff) << 16) | ((long) (u32[1] & 0xff) << 8) | ((long) (u32[0] & 0xff));
#ifdef debug
            printf("Object ID (in generic %08X): %08X has type: %04X\n", objId, id, tag);
#endif
            //			if(tag == kCDXObj_Text)
            //			{
            //				readText(ifs, id);
            //				depth--;
            //			}
            if((tag == kCDXObj_BracketAttachment) || (tag == kCDXObj_CrossingBond) || (tag == kCDXObj_BracketedGroup) || (tag == kCDXObj_Text)|| (tag == kCDXObj_Fragment))
              {	// Objects to ignore (Fragment could possibly be a problem)
                readGeneric(ifs, id);
                depth--;
              }
            else
              printf("New object in generic, type %04X\n", tag);
          }
        else if(tag == 0)	// End of object
          {
            depth--;
#ifdef debug
            printf("End of Object in generic %08X\n", objId);
#endif
          }
        else	// Property
          {
            ifs->read((char *)&size, sizeof(size));
            size = (size << 8) | (size >> 8);
#ifdef debug
            printf("Generic Tag: %04X\tSize: %04X\n", tag, size);
#endif
            switch(tag)
              {
              default: ifs->seekg(size, ios_base::cur); break;
              }
          }
        if(depth < 1)
          {
            return 0;
          }
      }
    return -1;
  }

  int ChemDrawBinaryFormat::readFragment(istream *ifs, UINT32 fragmentId, OBMol *pmol, map<UINT32, int> &atoms, list<cdBond> &bonds)
  {
    //	char buffer[1024];
    UINT16 tag;
    UINT16 size;
    UINT32 id;
    char u32[4];
    int depth = 1;

    atoms[fragmentId] = -1;
    while(ifs->good())
      {
        ifs->read((char *)&tag, sizeof(tag));
        tag = (tag << 8) | (tag >> 8);
        if(tag & kCDXTag_Object)	// Object
          {
            depth++;
            ifs->read(u32, sizeof(u32));
            id = ((long) (u32[3] & 0xff) << 24) | ((long) (u32[2] & 0xff) << 16) | ((long) (u32[1] & 0xff) << 8) | ((long) (u32[0] & 0xff));
#ifdef debug
            printf("Object ID (in fragment %08X): %08X has type: %04X\n", fragmentId, id, tag);
#endif
            if(tag == kCDXObj_Fragment)
              {
                if(readFragment(ifs, id, pmol, atoms, bonds) != 0)
                  {
                    printf("Error reading fragment\n");
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
                printf("New object in fragment, type %04X\n", tag);
              }

          }
        else if(tag == 0)	// End of object
          {
            depth--;
          }
        else	// Property
          {
            ifs->read((char *)&size, sizeof(size));
            size = (size << 8) | (size >> 8);
            //			printf("Fragment Tag: %04X\tSize: %04X\n", tag, size);
            switch(tag)
              {
              case kCDXProp_Frag_ConnectionOrder: ifs->seekg(size, ios_base::cur); break;  // Should use
              case kCDXProp_BoundingBox: ifs->seekg(size, ios_base::cur); break;
              default: ifs->seekg(size, ios_base::cur); printf("Fragment Tag: %04X\tSize: %04X\n", tag, size); break;
              }
          }
        if(depth < 1)
          {
            //			if(atoms.empty())
            //				printf("In readFragment: Atoms is empty\n");
            //			else
            //				printf("In readFragment: Atoms holds %d entries\n", atoms.size());
            return 0;
          }
      }
    return -1;
  }

} //namespace OpenBabel


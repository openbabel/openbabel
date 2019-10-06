/**********************************************************************
Copyright (C) 2006 by Fredrik Wallner
Some portions Copyright (C) 2006-2007 by Geoffrey Hutchsion
Some portions Copyright (C) 2011 by Chris Morley

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
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/elements.h>
#include <openbabel/reactionfacade.h>
#include <openbabel/stereo/stereo.h>
#include <openbabel/obfunctions.h>
#include <openbabel/reaction.h>
#include <openbabel/tokenst.h>
#include <openbabel/alias.h>
#include <openbabel/text.h>
#include "chemdrawcdx.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <list>


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

//Class which traverse the tree in CDX binary files 
class CDXReader
{
public:
  CDXReader(std::istream& is);
  CDXTag ReadNext(bool objectsOnly=false, int targetDepth=-2);
  void IgnoreObject()          { ReadNext(true, GetDepth()-1); }
  operator bool ()const        { return (bool)ifs; }
  int GetDepth()const          { return depth; }
  int GetLen()const            { return _len;} //length of current property data
  CDXObjectID CurrentID()const { return ids.back(); }
  stringstream& data(); //call this only once for each set of property data

  //Routines to display the structure of a cdx binary file
  OBText* WriteTree(const std::string& filename, unsigned wtoptions);
private:
  bool ParseEnums(std::map<CDXTag, std::string>& enummap, const std::string& filename);
  std::string TagName(std::map<CDXTag, std::string>& enummap, CDXTag tag);

private:
  std::istream& ifs;
  int depth;
  std::vector<CDXObjectID> ids;
  CDXObjectID _tempback;
  std::string _buf;
  UINT16 _len;
  std::stringstream _ss;
};

//**************************************************************
class ChemDrawBinaryXFormat : OBMoleculeFormat
{
public:
  //Register this format type ID in the constructor
  ChemDrawBinaryXFormat()
  {
    OBConversion::RegisterFormat("cdx",this);
  }

  virtual const char* Description() //required
  {
    return
      "ChemDraw binary format\n"
      "Read only\n"
      "The whole file is read in one call.\n"
      "Note that a file may contain a mixture of reactions and\n"
      "molecules.\n"

      "With the -ad option, a human-readable representation of the CDX tree\n"
      "structure is output as an OBText object. Use textformat to view it::\n\n"

      "    obabel input.cdx -otext -ad\n\n"

      "Many reactions in CDX files are not fully specified with reaction data\n"
      "structures, and may not be completely interpreted by this parser.\n\n"

      "Read Options, e.g. -am\n"
      " m read molecules only; no reactions\n"
      " d output CDX tree to OBText object\n"
      " o display only objects in tree output\n";
  }

  virtual const char* SpecificationURL()
  {return "http://www.cambridgesoft.com/services/documentation/sdk/chemdraw/cdx/IntroCDX.htm";}

  virtual const char* GetMIMEType()
  { return "chemical/x-cdx"; };

  virtual unsigned int Flags()
  {
    return READBINARY|NOTWRITABLE;
  }

  ////////////////////////////////////////////////////
  virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);

private:
  enum graphicType {none, equilArrow};
  bool        TopLevelParse(CDXReader& cdxr, OBConversion* pConv,CDXObjectID ContainingGroup);
  bool        DoFragment(CDXReader& cdxr, OBMol* pmol);
  bool        DoFragmentImpl(CDXReader& cdxr, OBMol* pmol,
         map<CDXObjectID, unsigned>& atommap, map<OBBond*, OBStereo::BondDirection>& updown);
  bool        DoReaction(CDXReader& cdxr, OBMol* pReact);
  std::string DoText(CDXReader& cdxr);

  std::vector<OBMol*> LookupMol(CDXObjectID id);
  graphicType         LookupGraphic(CDXObjectID id);
  OBMol*              LookupInMolMap(CDXObjectID id);

private:
  bool readReactions;
  static const bool objectsOnly = true;
  std::map<CDXObjectID, graphicType> _graphicmap;
  std::map<CDXObjectID, OBMol*> _molmap;
  std::map<CDXObjectID, std::vector<CDXObjectID> > _groupmap;
  // In case of chain A -> B -> C, B is both reactant and product
  CDXObjectID _lastProdId;
  typedef std::map<CDXObjectID, std::vector<CDXObjectID> >::iterator GroupMapIterator;
  static const unsigned usedFlag = 1<<30;
};

//******************************************************************
  //Global instance of the format
 ChemDrawBinaryXFormat theChemDrawBinaryXFormat;
//******************************************************************

 /*New CDXformat
Each fragment goes into a new OBMol on the heap.
The CDX id and OBMol* are added to _molmap.
When a reaction is found, the reactant/product/agent CDX ids are looked up in molmap,
and added to an OBReaction (made by deleting pOb if it is a OBMol
and assigning pOb to a new OBReaction. The OBMol is marked as Used.
When the reaction is complete it is output via AddChemObject().
At the end, any OBMol in the map not marked as Used is output as an OBMol.
*/


bool ChemDrawBinaryXFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
{
  _molmap.clear();
  _graphicmap.clear();
  _groupmap.clear();
  OBMol* pmol=NULL;
  bool ok = true;

  CDXReader cdxr(*pConv->GetInStream());
  readReactions = (pConv->IsOption("m",OBConversion::INOPTIONS)==NULL);

  // Write CDX tree only if requested
  if(pConv->IsOption("d",OBConversion::INOPTIONS))
  {
    unsigned wtoptions=0;
    if(pConv->IsOption("o",OBConversion::INOPTIONS))
      wtoptions |= 1; //display objects only
    OBText* pText  = cdxr.WriteTree("chemdrawcdx.h", wtoptions);
    if(pText)
    {
      pConv->AddChemObject(pText);
      return true;
    }
    return false;
  }

  // Normal reading of molecules and reactions
  //Top level parse 
  while(cdxr)
  {
    if(!TopLevelParse(cdxr, pConv, 0))
      return false;
  }

  //At the end, output molecules that have not been used in a reaction
  map<CDXObjectID, OBMol*>::iterator mapiter = _molmap.begin();
  for(; mapiter!=_molmap.end(); ++mapiter)
  {
    pmol = mapiter->second;
    if(!(pmol->GetFlags() & usedFlag) && strcmp(pmol->GetTitle(),"justplus"))
    {
      OBMol* ptmol = static_cast<OBMol*>(pmol->DoTransformations(
                    pConv->GetOptions(OBConversion::GENOPTIONS),pConv));
      if(!ptmol)
        delete pmol;
      else
        if(!pConv->AddChemObject(ptmol))
          return false; //error during writing
    }
  }

  return ok;
}
///////////////////////////////////////////////////////////////////////
bool ChemDrawBinaryXFormat::TopLevelParse
        (CDXReader& cdxr, OBConversion* pConv, CDXObjectID ContainingGroup)
{
  bool ok = true;
  CDXTag tag;
  while((tag = cdxr.ReadNext(objectsOnly)))
  {
    if(tag == kCDXObj_Group)
    {
      CDXObjectID cid = cdxr.CurrentID();
      vector<CDXObjectID> v;
      _groupmap.insert(make_pair(cid,v)); //empty vector as yet
      TopLevelParse(cdxr, pConv, cid );
    }

    else if(tag==kCDXObj_Fragment)
    {
      OBMol* pmol = new OBMol;
      //Save all molecules to the end
      _molmap[cdxr.CurrentID()] = pmol;

      if(ContainingGroup)
      {
        // Add the id of this mol to the group's entry in _groupmap 
        GroupMapIterator gmapiter = _groupmap.find(ContainingGroup);
        if(gmapiter!=_groupmap.end())
          gmapiter->second.push_back(cdxr.CurrentID());
      }
      ok = DoFragment(cdxr, pmol);
    }

    else if(tag == kCDXObj_ReactionStep && readReactions)
    {
      OBMol* pReact = new OBMol;
      pReact->SetIsReaction();
      ok = DoReaction(cdxr, pReact);
      // Output OBReaction and continue 
      if(pReact)
        if(!pConv->AddChemObject(pReact))
          return false; //error during writing
    }

    else if(ok && tag==kCDXObj_Graphic)
    {
      while( (tag = cdxr.ReadNext()) )
      {
        stringstream& ss = cdxr.data();
        if(tag == kCDXProp_Arrow_Type)
        {
          char type1=0;
          UINT16 type2=0;
          if(cdxr.GetLen()==1)
            ss.get(type1);
          else
            READ_INT16(ss,type2);
          if(type1==kCDXArrowType_Equilibrium || type2==kCDXArrowType_Equilibrium)
            _graphicmap[type1+type2] = equilArrow; //save in graphicmap
        }
      }
    }
  }
  return true;
}
///////////////////////////////////////////////////////////////////////
bool ChemDrawBinaryXFormat::DoReaction(CDXReader& cdxr, OBMol* pReact)
{
  CDXTag tag;
  CDXObjectID id;
  OBReactionFacade facade(pReact);
  while( (tag = cdxr.ReadNext()) )
  {
    if(tag ==	kCDXProp_ReactionStep_Reactants)
    {
      stringstream& ss = cdxr.data();
      for(unsigned i=0;i<cdxr.GetLen()/4;++i)//for each reactant id
      {
        READ_INT32(ss,id);
        vector<OBMol*> molvec = LookupMol(id); //id could be a group with several mols
        for(unsigned i=0;i<molvec.size();++i)
          if(strcmp(molvec[i]->GetTitle(),"justplus"))
          {
            facade.AddComponent(molvec[i], REACTANT);
          }
      }
    }
    else if(tag == kCDXProp_ReactionStep_Products)
    {
      stringstream& ss = cdxr.data();
      for(unsigned i=0;i<cdxr.GetLen()/4;++i)//for each product id
      {
        READ_INT32(ss,id);
        vector<OBMol*> molvec = LookupMol(id); //id could be a group with several mols
        for(unsigned i=0;i<molvec.size();++i)
          if(strcmp(molvec[i]->GetTitle(),"justplus"))
          {
            facade.AddComponent(molvec[i], PRODUCT);
            _lastProdId = id;
          }
      }
    }
    else if(tag==kCDXProp_ReactionStep_Arrows)
    {
      READ_INT32(cdxr.data(),id);
      //if(LookupGraphic(id)==equilArrow) // TODO? Store reversibility somehow?
      //  pReact->SetReversible();
    }
  }
  return true;
}
///////////////////////////////////////////////////////////////////////
vector<OBMol*> ChemDrawBinaryXFormat::LookupMol(CDXObjectID id)
{
  vector<OBMol*> molvec;
  //Check whether the id is that of a kCDXObj_Group
  GroupMapIterator gmapiter;
  gmapiter = _groupmap.find(id);
  if(gmapiter != _groupmap.end())
  {
    for(unsigned i=0;i<gmapiter->second.size();++i)
    {
      OBMol* pmmol = LookupInMolMap(gmapiter->second[i]);
      if(pmmol)
        molvec.push_back(pmmol);
    }
  }
  else
  {
    //id is not a group; it must be a fragment
    OBMol* pmmol = LookupInMolMap(id);
    if(pmmol)
      molvec.push_back(pmmol);
  }
  return molvec; 
}

OBMol* ChemDrawBinaryXFormat::LookupInMolMap(CDXObjectID id)
{
  std::map<CDXObjectID, OBMol*>::iterator mapiter;
  mapiter = _molmap.find(id);
  if(mapiter!=_molmap.end())
  {
    //Mark mol as used in a reaction, so that it will not be output independently
    mapiter->second->SetFlags(mapiter->second->GetFlags() | usedFlag);
    return mapiter->second;
  }
  else
  {
    stringstream ss;
    ss << "Reactant or product mol not found id = " << hex << showbase << id; 
    obErrorLog.ThrowError(__FUNCTION__, ss.str(), obError);
    return NULL;
  }
}

////////////////////////////////////////////////////////////////////////
ChemDrawBinaryXFormat::graphicType ChemDrawBinaryXFormat::LookupGraphic(CDXObjectID id)
{
  std::map<CDXObjectID, graphicType>::iterator mapiter;
  mapiter = _graphicmap.find(id);
  if(mapiter != _graphicmap.end())
    return mapiter->second;
  else
    return none;
}

////////////////////////////////////////////////////////////////////////
bool ChemDrawBinaryXFormat::DoFragment(CDXReader& cdxr, OBMol* pmol)
{
  map<OBBond*, OBStereo::BondDirection> updown;
  pmol->SetDimension(2);
  pmol->BeginModify();

  map<CDXObjectID, unsigned> atommap; //key = CDX id; value = OB atom idx

  //The inner workings of DoFragment,since Fragment elements can be nested
  DoFragmentImpl(cdxr, pmol, atommap, updown);

  // use 2D coordinates + hash/wedge to determine stereochemistry
  StereoFrom2D(pmol, &updown);

  pmol->EndModify();

  //Expand any aliases after molecule constructed
  //Need to save aliases in list first and expand later
  vector<OBAtom*> aliasatoms;
  for(int idx=1; idx<=pmol->NumAtoms();++idx)
  {
    OBAtom* pAtom = pmol->GetAtom(idx);
    AliasData* ad = dynamic_cast<AliasData*>(pAtom->GetData(AliasDataType));
    if(ad && !ad->IsExpanded())
      aliasatoms.push_back(pAtom);
  }
  for(vector<OBAtom*>::iterator vit=aliasatoms.begin();
      vit!=aliasatoms.end(); ++vit)
  {
    int idx = (*vit)->GetIdx();
    AliasData* ad = dynamic_cast<AliasData*>((*vit)->GetData(AliasDataType));
    if(ad && !ad->IsExpanded())
      ad->Expand(*pmol, idx); //Make chemically meaningful, if possible.
  }
  return true;
}

bool ChemDrawBinaryXFormat::DoFragmentImpl(CDXReader& cdxr, OBMol* pmol, 
       map<CDXObjectID, unsigned>& atommap, map<OBBond*, OBStereo::BondDirection>& updown)
{
  CDXTag tag;
  std::vector<OBAtom*> handleImplicitCarbons;
  while((tag = cdxr.ReadNext(objectsOnly)))
  {
    if(tag==kCDXObj_Node)
    {
      unsigned nodeID = cdxr.CurrentID();
      bool isAlias=false, hasElement=false;
      bool hasNumHs = false;
      UINT16 atnum=-1, spin=0, numHs=0;
      int x, y, charge=0, iso=0;
      string aliastext;

      //Read all node properties
      while( (tag = cdxr.ReadNext()) )
      {
        switch(tag)
        {
        case kCDXProp_Node_Type:
          UINT16 type;
          READ_INT16(cdxr.data(), type);
          if(type==4 || type==5) //Nickname or fragment
            isAlias = true;
          break;
        case kCDXProp_Node_Element:
          READ_INT16(cdxr.data(), atnum);
          hasElement = true;
          break;
        case kCDXProp_2DPosition:
          {
            stringstream& ss = cdxr.data();
            READ_INT32(ss, y); //yes, this way round
            READ_INT32(ss, x);
          }
            break;
        case kCDXProp_Atom_Charge:
          if(cdxr.GetLen()==1)
            charge = cdxr.data().get();
          else
            READ_INT32(cdxr.data(), charge);
          break;
        case kCDXProp_Atom_Radical:
          READ_INT16(cdxr.data(),spin);
          break;
        case kCDXProp_Atom_Isotope:
          READ_INT16(cdxr.data(),iso);
          break;
        case kCDXProp_Atom_NumHydrogens:
          READ_INT16(cdxr.data(), numHs);
          hasNumHs = true;
          break;
        case kCDXProp_Atom_CIPStereochemistry:
          break;
        case kCDXObj_Text:
          aliastext = DoText(cdxr);
          if(aliastext=="+")
          {
            //This node is not an atom, but dangerous to delete
            pmol->SetTitle("justplus");
          }
          break;
        case kCDXObj_Fragment:
        /* ignore fragment contained in node
        if(isAlias)
          {
            unsigned Idxbefore = pmol->NumAtoms();
            if(DoFragmentImpl(cdxr, pmol, atommap, updown))
              return false;
          }
         */
          //ignore the contents of this node
          cdxr.IgnoreObject();
          //cdxr.ReadNext(objectsOnly, cdxr.GetDepth()-1);
          break;
        default:
          if(tag & kCDXTag_Object) //unhandled object
            while(cdxr.ReadNext());
        }
      }
      //All properties of Node have now been read
      OBAtom* pAtom = pmol->NewAtom();
      pAtom->SetVector(x*1.0e-6, -y*1.0e-6, 0); //inv y axis
      atommap[nodeID] = pmol->NumAtoms();
      if(isAlias || (!aliastext.empty() && atnum==0xffff))
      {
        //Treat text as an alias 
        pAtom->SetAtomicNum(0);
        AliasData* ad = new AliasData();
        ad->SetAlias(aliastext);
        ad->SetOrigin(fileformatInput);
        pAtom->SetData(ad);
      } 
      else
      {
        if(atnum==0xffff)
          atnum = 6; //atoms are C by default
        pAtom->SetAtomicNum(atnum);
        if (hasNumHs)
          pAtom->SetImplicitHCount(numHs);
        else if (atnum==6)
          handleImplicitCarbons.push_back(pAtom);
        pAtom->SetFormalCharge(charge);
        pAtom->SetIsotope(iso);
        pAtom->SetSpinMultiplicity(spin);
      }
    }

    else if(tag==kCDXObj_Bond)
    {
      CDXObjectID bgnID, endID;
      int order=1, bgnIdx, endIdx ;
      UINT16 stereo=0;

      while( (tag = cdxr.ReadNext()) )
      {
        switch(tag)
        {
        case kCDXProp_Bond_Begin:
          READ_INT32(cdxr.data(), bgnID);
          bgnIdx = atommap[bgnID];
          break;
        case kCDXProp_Bond_End:
          READ_INT32(cdxr.data(), endID);
          endIdx = atommap[endID];
          break;
        case kCDXProp_Bond_Order:
          READ_INT16(cdxr.data(), order);
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
          READ_INT16(cdxr.data(), stereo);
        break;
        }
      }

      if(!order || !bgnIdx || !endIdx)
      {
        obErrorLog.ThrowError(__FUNCTION__,"Incorrect bond", obError);
        return false;
      }
      if(stereo==4 || stereo==7 || stereo==10 || stereo==12)
        swap(bgnIdx, endIdx);
      pmol->AddBond(bgnIdx, endIdx, order);
      if(stereo)
      {
        OBBond* pBond = pmol->GetBond(pmol->NumBonds()-1);
        if(stereo==3 || stereo==4)
          pBond->SetHash();
        else if(stereo==6 || stereo==7)
          pBond->SetWedge();
      }
    }
  }
  // Handle 'implicit carbons' by adjusting their valence with
  // implicit hydrognes
  for(vector<OBAtom*>::iterator vit=handleImplicitCarbons.begin();
      vit!=handleImplicitCarbons.end(); ++vit)
    OBAtomAssignTypicalImplicitHydrogens(*vit);

  return true;
}

string ChemDrawBinaryXFormat::DoText(CDXReader& cdxr)
{
  CDXTag tag;
  string text;
  while( (tag=cdxr.ReadNext()) )
  {
    stringstream& ss = cdxr.data();
    switch(tag)
    {
    case kCDXProp_Text:
      UINT16 nStyleRuns;
      READ_INT16(ss,nStyleRuns);
      ss.ignore(nStyleRuns*10);
      ss >> text;
    default:
      if(tag & kCDXTag_Object) //unhandled object
        while(cdxr.ReadNext());      
    }
  }
  return text;
}

//****************************************************************
CDXTag CDXReader::ReadNext(bool objectsOnly, int targetDepth)
{
  //ostringstream treestream;
  CDXTag tag;
  CDXObjectID id;

  while(ifs) 
  {
    READ_INT16(ifs, tag);
    if(tag==0)
    {
      if(depth==0)
      {
        ifs.setstate(ios::eofbit); //ignore everything after end of document
        return 0; //end of document
      }
      --depth;
      _tempback = ids.back(); //needed for WriteTree
      ids.pop_back();
      if(targetDepth<0 || depth == targetDepth)
        return 0; //end of object
    }
    else if(tag & kCDXTag_Object)
    {
      READ_INT32(ifs, id);
      ids.push_back(id);
      ++depth;
      if(targetDepth<0 || depth-1 == targetDepth)
        return tag; //object
    }
    else
    {
      //property
      READ_INT16(ifs, _len);

      if(objectsOnly)
        ifs.ignore(_len);
      else
      {
        //copy property data to buffer
        char* p = new char[_len+1];
        ifs.read(p, _len);
        _buf.assign(p, _len);
        delete[] p;
        return tag; //property
      }
    }
  }
  return 0;
}
/////////////////////////////////////////////////////////////////////

stringstream& CDXReader::data()
{
  _ss.clear();
  _ss.str(_buf);
  return _ss;
}
/////////////////////////////////////////////////////////////////////

CDXReader::CDXReader(std::istream& is) : ifs(is), depth(0)
{
  //ReadHeader
  char buffer[kCDX_HeaderStringLen+1];
  ifs.read(buffer,kCDX_HeaderStringLen);
  buffer[kCDX_HeaderStringLen] = '\0';
  if(strncmp(buffer, kCDX_HeaderString, kCDX_HeaderStringLen) == 0)
    ifs.ignore(kCDX_HeaderLength - kCDX_HeaderStringLen);	// Discard rest of header.
  else
  {
    obErrorLog.ThrowError(__FUNCTION__,"Invalid file, no ChemDraw Header",obError);
    ifs.setstate(ios::eofbit);
    throw;
  }
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Routines to display the structure of a cdx binary file

OBText* CDXReader::WriteTree(const string& filename, unsigned wtoptions)
{  
  const char indentchar = '\t';
  std::map<CDXTag, std::string> enummap;
  ParseEnums(enummap, filename);

  stringstream tss;
  tss << hex << showbase;

  while(*this)
  {
    CDXTag tag = ReadNext();
    if(ifs.eof())
      return new OBText(tss.str()); //normal exit
    if(tag==0 && !(wtoptions &1))
    {
      //Object end
      tss << string(depth,indentchar) << "ObjectEnd " << _tempback << endl;
    } 
    else if(tag & kCDXTag_Object)
    {
      //Object
      tss<<string(depth-1,indentchar) << "Object " << tag
                   << TagName(enummap,tag) << " id=" << ids.back() << endl; 
    }
    else
    {
      //Property
      if(!(wtoptions &1))
      {
        stringstream ss;
        ss << _len;
        tss<<string(depth,indentchar) << "Property  "<< tag << TagName(enummap,tag)
                     << " [" << ss.str() << " bytes] ";
        for(unsigned i=0;i<_len;++i)
        {
          ss.str("");
          ss.fill('0');
          ss.width(8);
          ss << hex << static_cast<unsigned>(_buf[i]) << dec;
          tss << ss.str()[6] << ss.str()[7] << ' ';
        }

        if(tag==0x700 || tag==kCDXProp_CreationProgram || tag==kCDXProp_CreationDate
          || tag==kCDXProp_Name)
        {
          stringstream ss(_buf);
          UINT16 nStyleRuns;
          READ_INT16(ss, nStyleRuns);
          tss << '\"';
          for(unsigned i=2+nStyleRuns*10; i<_len; ++i)
            tss << _buf[i];
          tss << '\"';
        }
        tss << endl;
      }
    }
  }
  return NULL; //error exit
}

///////////////////////////////////////////////////////////////////////
bool CDXReader::ParseEnums(map<CDXTag, string>& enummap, const string& filename)
{
  ifstream ihs;
  if(OpenDatafile(ihs, filename).empty())
  {
    obErrorLog.ThrowError(__FUNCTION__, 
      filename + " needs to be in the *data* directory when displaying the tree.\n" , obError);
    return false;
  }
  ignore(ihs, "enum CDXDatumID");
  string ln;
  vector<string> vec;
  stringstream ss;
  CDXTag tag;
  while(ihs)
  {
    getline(ihs, ln);
    tokenize(vec, ln, " \t,{}");
    if(vec.size()==0 || vec[0]=="//")
      continue; //blank and comment lines
    if(vec[0]==";") //line is }; end of enum
      return true;
    if(vec[0][0]!='k') //only collect enums starting with kCDX
      continue;
    int tagpos = (vec[1]=="=" && vec.size()>4) ? 4 : 2;
    ss.str(vec[tagpos]);
    ss.clear();
    ss >> hex >> tag;
    if(ss)
    {
      if(tag==0x0400 && vec[0]=="kCDXUser_TemporaryEnd")//special case
        continue;
      enummap[tag] = vec[0];
    }
  }
  return false;
}
/////////////////////////////////////////////////////////////////////////

string CDXReader::TagName(map<CDXTag, string>& enummap, CDXTag tag)
{
  string tagname;
  if(!enummap.empty())
  {
    map<CDXTag, std::string>::iterator iter = enummap.find(tag);
    if(iter!=enummap.end())
    {
      tagname=iter->second;
      //Remove prefix, e.g. kCDXProp_
      string::size_type pos = tagname.find('_');
      if(pos!=string::npos)
      {
        tagname.erase(0,pos);
        tagname[0] = ' ';
      }
    }
  }
  return tagname;
}

} //namespace

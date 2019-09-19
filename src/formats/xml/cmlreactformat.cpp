/**********************************************************************
Copyright (C) 2005 by Chris Morley

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include "openbabel/babelconfig.h"
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obiter.h>
#include <openbabel/elements.h>
#include "openbabel/reaction.h"
#include "openbabel/xml.h"
#include "openbabel/kinetics.h"
#include "openbabel/text.h"

using namespace std;

namespace OpenBabel
{

class CMLReactFormat : XMLBaseFormat
{
public:
  CMLReactFormat()
  {
    OBConversion::RegisterFormat("cmlr",this);
    XMLConversion::RegisterXMLFormat( this);
    OBConversion::RegisterOptionParam("l", this);
    XMLConversion::RegisterXMLFormat(this, false,"http://www.xml-cml.org/schema/cml2/react");
  }
  virtual const char* NamespaceURI()const
  {return "http://www.xml-cml.org/schema";}


  const char* Description()
  {
    return
      "CML Reaction format\n"
      "A minimal implementation of the CML Reaction format\n"
      "This implementation uses libxml2.\n"
      "Write Options (e.g. -x1a)\n"
      " 1  output CML1 (rather than CML2)\n"
      " a  output array format for atoms and bonds\n"
      " l  molecules NOT in MoleculeList\n"
      " h  use hydrogenCount for all hydrogens\n"
      " x  omit XML declaration\n"
      " r  omit rate constant data\n"
      " N<prefix> add namespace prefix to elements\n"
      " M  add obr prefix on non-CMLReact elements\n"
      " p  add properties to molecules\n\n"

      "The implementation of this format which reads and writes to and from\n"
      "OBReaction objects is fairly minimal at present. (Currently the only\n"
      "other reaction format in OpenBabel is RXN.) During reading, only the\n"
      "elements <reaction>, <reactant>, <product> and <molecule>  are acted\n"
      "upon (the last through CML). The molecules can be collected together\n"
      "in a list at the start of the file and referenced in the reactant and\n"
      "product via e.g. <molecule ref=\"mol1\">.\n\n"

      "On writing, the list format can be specified with the ``-xl`` option. The\n"
      "list containers are <moleculeList> and <reactionList> and the overall\n"
      "wrapper is <mechanism>. These are non-standard CMLReact element names\n"
      "and would have to be changed (in the code) to <list>,<list> and <cml>\n"
      "if this was unacceptable.\n\n";
  }

  virtual const char* TargetClassDescription()
  {
      return OBReaction::ClassDescription();
  }

  unsigned Flags()
  {
    return 0;
  }
  virtual bool ReadChemObject(OBConversion* pConv);
  virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
  virtual bool WriteChemObject(OBConversion* pConv);
  virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);
  virtual bool DoElement(const string& ElName);
  virtual bool EndElement(const string& ElName);
  virtual const char* EndTag(){ return "/reaction>"; };
  const type_info& GetType()
  {
    return typeid(OBReaction*);
  };

private:
  typedef map<string,obsharedptr<OBMol> > MolMap;
  string AddMolToList(obsharedptr<OBMol> spmol, MolMap& mmap);
  bool WriteRateData(OBReaction* pReact, xmlChar* altprefix);
  void WriteMetadataList(OBReaction& react);

private:
  OBReaction* _preact;
  OBMol* pmol;
  obsharedptr<OBMol> _spmol;
  MolMap IMols; //used on input
  MolMap OMols; //used on output
  int nextmol;
  ostringstream ssout; //temporary output
  bool MolsHaveBeenAdded; //separate OBMols during output
  OBRateData* _pRD; //used on input
  string _text; //used on output to delay text so it can be after reactions
  ostream* _pOut;//used on output to save the original output stream
};

//Make an instance of the format class
CMLReactFormat theCMLReactFormat;

////////////////////////////////////////////////////
bool CMLReactFormat::ReadChemObject(OBConversion* pConv)
{
  //Makes a new OBReaction and new associated OBMols
  OBReaction* pReact = new OBReaction;

  if(pConv->IsFirstInput())
  {
    IMols.clear();
    //add special species
    obsharedptr<OBMol> sp(new OBMol);
    sp.get()->SetTitle("M");
    IMols["M"] = sp;
  }

  bool ret=ReadMolecule(pReact,pConv); //call the "API" read function

  std::string auditMsg = "OpenBabel::Read reaction ";
  std::string description(Description());
  auditMsg += description.substr(0,description.find('\n'));
  obErrorLog.ThrowError(__FUNCTION__,
            auditMsg,
            obAuditMsg);

   //Do transformation and return reaction, if it has either reactants or products
  if(ret && (pReact->NumReactants()!=0 || pReact->NumProducts()!=0)) //Do transformation and return molecule
    return pConv->AddChemObject(pReact->DoTransformations(pConv->GetOptions(OBConversion::GENOPTIONS),pConv))!=0;
  else
  {
    delete pReact;
    pConv->AddChemObject(NULL);
    return false;//don't continue after empty reaction
  }
}

bool CMLReactFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
{
  //Not a molecule: the object converted is a reaction
    _preact = dynamic_cast<OBReaction*>(pOb);
    if(!_preact)
      return false;

    _pxmlConv = XMLConversion::GetDerived(pConv);
    if(!_pxmlConv)
      return false;

    _pRD = NULL; //No rate data yet

    bool ret =_pxmlConv->ReadXML(this,pOb);

    return ret;
}


bool CMLReactFormat::DoElement(const string& name)
{
  double val;

  if(name=="reaction")
  {
    _spmol.reset();
    _preact->SetTitle(_pxmlConv->GetAttribute("id"));
  }

  else if(name=="molecule")
  {
    string reference = _pxmlConv->GetAttribute("ref");
    if(!reference.empty())
    {
      _spmol = IMols[reference];
      pmol = _spmol.get();
      if(!pmol)
      {
        cerr << " Molecule reference \"" << reference <<"\" not found" << endl;
        return false;
      }
    }
    else
    {
      obsharedptr<OBMol> sp(new OBMol);
      OBFormat* pCMLFormat = OBConversion::FindFormat("cml");
      if(!pCMLFormat)
        return false;
      _pxmlConv->_SkipNextRead=true;
      pCMLFormat->ReadMolecule(sp.get(), _pxmlConv);

      //Store smart pointers to all molecules in map
      _spmol=sp;
      AddMolToList(_spmol,IMols);
    }
  }
  else if(name=="rateParameters")
  {
    _pRD = new OBRateData; //to store rate constant data
   _preact->SetData(_pRD);

    string rt = _pxmlConv->GetAttribute("reactionType");
    OBRateData::reaction_type enumrt=OBRateData::ARRHENIUS;
    if(rt=="arrhenius")       enumrt=OBRateData::ARRHENIUS;
    else if(rt=="lindermann") enumrt=OBRateData::LINDERMANN;
    else if(rt=="troe")       enumrt=OBRateData::TROE;
    else if(rt=="sri")        enumrt=OBRateData::SRI;
    else if(rt=="threeBody")  enumrt=OBRateData::THREEBODY;
    else
      obErrorLog.ThrowError(__FUNCTION__, rt + " is not a known reactionType", obWarning);
    _pRD->ReactionType = enumrt;

    if(_pxmlConv->GetAttribute("reversible")=="true")
      _preact->SetReversible();
  }

  else if(_pRD && name=="A")
  {
    if(_pxmlConv->GetContentDouble(val))
      _pRD->SetRate(OBRateData::A, val);
  }
  else if(_pRD && name=="n")
  {
    if(_pxmlConv->GetContentDouble(val))
      _pRD->SetRate(OBRateData::n, val);
  }
  else if(_pRD && name=="E")
  {
    if(_pxmlConv->GetContentDouble(val))
      _pRD->SetRate(OBRateData::E, val);
  }
  else if(_pRD && name=="loA")
  {
    if(_pxmlConv->GetContentDouble(val))
      _pRD->SetLoRate(OBRateData::A, val);
  }
  else if(_pRD && name=="lon")
  {
    if(_pxmlConv->GetContentDouble(val))
      _pRD->SetLoRate(OBRateData::n, val);
  }
  else if(_pRD && name=="loE")
  {
    if(_pxmlConv->GetContentDouble(val))
      _pRD->SetLoRate(OBRateData::E, val);
  }
  else if(_pRD && name=="troeParams")
  {
    string txt(_pxmlConv->GetContent());
    if(!txt.empty())
    {
      stringstream ss(txt);
      for(int i=0;i<4;++i)
      {
        ss >>val;
        _pRD->SetTroeParams(i, val);
      }
    }

  }
  else if(_pRD && name=="eff")
  {
    string ref = _pxmlConv->GetAttribute("ref");
    if(!ref.empty() && _pxmlConv->GetContentDouble(val))
      _pRD->SetEfficiency(ref, val);
  }

  //The end element event would not be called for <element/>, so call it explicitly.
    if(xmlTextReaderIsEmptyElement(reader())==1)
      return EndElement(name);

  return true;
}

bool CMLReactFormat::EndElement(const string& name)
{
  if(name=="reactant")
  {
    if(!_spmol.get())
      return false;
    _preact->AddReactant(_spmol);
  }
  else if(name=="product")
  {
    if(!_spmol.get())
      return false;
    _preact->AddProduct(_spmol);
  }
  else if(name=="reaction")
  {
    _spmol.reset();
    return false;//means stop parsing
  }
  else if(name=="rateParameters")
    _pRD = NULL;

  return true;
}

///////////////////////////////////////////////////////////////////
bool CMLReactFormat::WriteChemObject(OBConversion* pConv)
{
  //WriteChemObject() always deletes the object retrieved by GetChemObject
  OBBase* pOb=pConv->GetChemObject();
  OBReaction* pReact = dynamic_cast<OBReaction*>(pOb);
  if(pReact==NULL)
  {
    //If sent a molecule, add it to the list.
    //Molecules added in this way have precedence over internal molecules.
    //They can therefore have structure and properties, perhaps obtained from
    //more than one place using the -C option. Thermodynamic data can also be
    //upgraded from that in molecules in OBReaction object (and so in Chemkin or
    //CMLReact input file) by including a separate file with updates.

    if(pConv->GetOutputIndex()==1)
    {
      _pOut = pConv->GetOutStream();//Save the original output stream
      OMols.clear();
    }

    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol!=NULL)
    {
      obsharedptr<OBMol> sp(pmol);
      AddMolToList(sp, OMols);
      pConv->SetOutputIndex(-1); //Signals that molecules have been added

      bool ret=true;
      if(pConv->IsLast())
      {
        //Only molecules have been supplied; output them via CMLFormat
        OBFormat* pCMLFormat = pConv->FindFormat("cml");
        if(pCMLFormat==NULL)
        {
          obErrorLog.ThrowError(__FUNCTION__,
            "CML format for molecules is needed by CMLReactformat and is not available\n",obError);
            return false;
        }
        unsigned int n=0;
        MolMap::iterator mapitr;
        for(mapitr=OMols.begin();mapitr!=OMols.end() && ret; ++mapitr)
        {
          pConv->SetOutputIndex(++n);      //we have to increment and
          pConv->SetLast(n==OMols.size()); //setLast manually because we are not using Convert()
          ret = pCMLFormat->WriteMolecule(mapitr->second.get(), pConv);
        }
      }
      return ret;
    }
    else
    //If sent text as an OBText object, output the text up to the insertion point
    //and the rest of the text in _text. The content of _text is output at the end.
    {
      OBText* ptext = dynamic_cast<OBText*>(pOb);
      if(ptext==NULL)
        return false;
      string::size_type pos = 0;
      string frontText(ptext->GetText(pos));
      *_pOut << frontText;//Output(to orig outstream) text up to insertion point
      _text = ptext->GetText(pos);  //Save text after insertion point to be output at the end

      //if any of the text contains an xml declaration, do not write again
      // and also omit <cml>...</cml> wrapper
      if(frontText.find("<?xml ")!=string::npos)
        pConv->AddOption("ReactionsNotStandalone");

      pConv->SetOutputIndex(pConv->GetOutputIndex()-1);//not an output we want to count
      return true;
    }
  }

  bool ret=false;
  ret=WriteMolecule(pReact,pConv);

  std::string auditMsg = "OpenBabel::Write reaction ";
  std::string description(Description());
  auditMsg += description.substr( 0, description.find('\n') );
  obErrorLog.ThrowError(__FUNCTION__, auditMsg, obAuditMsg);

  delete pOb;

  if(pConv->IsLast() && !_text.empty())
  {
    *_pOut << _text; //Last part of an OBText object to original output stream
    _text.clear();
  }

  return ret;
}

////////////////////////////////////////////////////////////////
bool CMLReactFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
  //Badly named function: it's really a reaction, not a molecule.

  _pxmlConv = XMLConversion::GetDerived(pConv,false);
  if(!_pxmlConv)
    return false;

  //Cast output object to the class type need, i.e. OBReaction
  OBReaction* pReact = dynamic_cast<OBReaction*>(pOb);
  if(pReact==NULL)
      return false;

  //Two forms of output are supported:
  //With normal form, if more than one reaction to be output,
  //write <cml> at start and </cml> at end.
  //With list form, use ref attributes when refering to molecules and put
  //actual molecules in a separate list.
  // <cml><moleculeList>...</moleculeList><reactionList>...</reactionList></cml>

  //The following element and attribute names may be output
  //For strict CMLReact compliance use the alternatives:
  static const xmlChar C_REACTIONLIST[] = "reactionList";
  static const xmlChar C_MOLECULELIST[] = "moleculeList";
  static const xmlChar C_LISTWRAPPER[]  = "cml";    //was non-std "mechanism"
  static const xmlChar C_WRAPPER[]      = "cml";
  static const xmlChar C_MOLECULE[]     = "molecule";
  static const xmlChar C_REACTION[]     = "reaction";
  static const xmlChar C_REACTANT[]     = "reactant";
  static const xmlChar C_PRODUCT[]      = "product";
  static const xmlChar C_TS[]           = "me:transitionState";
  static const xmlChar C_REACTANTLIST[] = "reactantList";
  static const xmlChar C_PRODUCTLIST[]  = "productList";
  static const xmlChar C_REF[]          = "ref";
  static const xmlChar C_TITLE[]        = "id"; //"title";
  static const xmlChar C_REVERSIBLE[]   = "reversible";


  bool list = _pxmlConv->IsOption("l")==NULL; //Output with molecules in a separate list

  xmlChar* prefix = BAD_CAST _pxmlConv->IsOption("N");
  xmlChar* altprefix = _pxmlConv->IsOption("M") ? BAD_CAST "obr" : NULL;

  xmlChar* uri=NULL;
  xmlChar* alturi=NULL;

  _pxmlConv->AddOption("MolsNotStandalone",OBConversion::OUTOPTIONS); //inform CMLFormat

//  OBConversion MolConv(*_pxmlConv); //new copy to use to write associated CML molecules
//	MolConv.SetAuxConv(NULL); //temporary until a proper OBConversion copy constructor written
//	MolConv.SetOneObjectOnly();

  OBFormat* pCMLFormat = _pxmlConv->FindFormat("cml");
  if(pCMLFormat==NULL)
  {
    obErrorLog.ThrowError(__FUNCTION__,
      "CML format for molecules is needed by CMLReactformat and is not available\n",obError);
    return false;
  }

  bool MolsHaveBeenAdded=false;
  if(pConv->GetOutputIndex()==0)
  {
    MolsHaveBeenAdded = true;
    pConv->SetOutputIndex(1);//set to proper rather than artificial value
  }

  //For first reaction
  if(pConv->GetOutputIndex()==1) //OBConversion::Convert() is still using original pConv
  {
    _pOut = pConv->GetOutStream();//Save the original output stream
    if(list)
    {
      //Use a temporary stringstream for list format so that moleculeList can be first
      ssout.str("");
      if(!MolsHaveBeenAdded)
        OMols.clear();

      _pxmlConv->SetOutStream(&ssout);
      //Do also for original object; GetDerived() copies the old to the new on subsequent reactions
      pConv->SetOutStream(&ssout);
    }
    if(!_pxmlConv->IsOption("x") && !_pxmlConv->IsOption("ReactionsNotStandalone"))
    {
      xmlTextWriterStartDocument(writer(), NULL, NULL, NULL);
      uri=BAD_CAST NamespaceURI();
      alturi = BAD_CAST "http://www.gaseq.co.uk/obr"; //@todo: better alturi
    }

    if(list)
    {
      if(!pConv->IsOption("ReactionsNotStandalone"))
      {
        xmlTextWriterStartElementNS(writer(), prefix, C_LISTWRAPPER, uri);
      if(altprefix)
        xmlTextWriterWriteAttributeNS(writer(),BAD_CAST "xmlns",altprefix,NULL,alturi);
      }
      xmlTextWriterStartElementNS(writer(), prefix, C_REACTIONLIST, NULL);
    }
    else if(!_pxmlConv->IsLast() && !_pxmlConv->IsOption("ReactionsNotStandalone"))
      xmlTextWriterStartElementNS(writer(), prefix, C_WRAPPER, uri);
    uri=NULL; //not needed again
  }

  xmlTextWriterStartElementNS(writer(), prefix, C_REACTION, NULL);
  if(!pReact->GetTitle().empty())
    xmlTextWriterWriteFormatAttribute(writer(), C_TITLE,"%s", pReact->GetTitle().c_str());
  if(pReact->IsReversible())
    xmlTextWriterWriteFormatAttribute(writer(), C_REVERSIBLE,"%s", "true");

  WriteMetadataList(*pReact);

  xmlTextWriterStartElementNS(writer(), prefix, C_REACTANTLIST, NULL);
  unsigned i;
  for(i=0;i<pReact->NumReactants();i++)
  {
    xmlTextWriterStartElementNS(writer(), prefix, C_REACTANT, NULL);
    // put molecules into map even if they are to be output immediately
    // since the map may already have the same molecules with additional info
    string id = AddMolToList(pReact->GetReactant(i), OMols);
    if(list) //output references
    {
      xmlTextWriterStartElementNS(writer(), prefix, C_MOLECULE, NULL);
      xmlTextWriterWriteFormatAttribute(writer(), C_REF,"%s", id.c_str());
      xmlTextWriterEndElement(writer());//molecule
    }
    else
    {
      //Write reactant in CML format
      pCMLFormat->WriteMolecule(OMols[id].get(), _pxmlConv);
    }
    xmlTextWriterEndElement(writer());//reactant
  }
  xmlTextWriterEndElement(writer());//reactantList

  xmlTextWriterStartElementNS(writer(), prefix, C_PRODUCTLIST, NULL);
  for(i=0;i<pReact->NumProducts();i++)
  {
    xmlTextWriterStartElementNS(writer(), prefix, C_PRODUCT, NULL);
    // put molecules into map even if they are to be output immediately
    // since the map may already have the same molecules with additional info
    string id = AddMolToList(pReact->GetProduct(i), OMols);
    if(list) //output references
    {
      xmlTextWriterStartElementNS(writer(), prefix, C_MOLECULE, NULL);
      xmlTextWriterWriteFormatAttribute(writer(), C_REF,"%s", id.c_str());
      xmlTextWriterEndElement(writer());//molecule
    }
    else
      //Write product in CML format
      pCMLFormat->WriteMolecule(OMols[id].get(), _pxmlConv);
    xmlTextWriterEndElement(writer());//product
 }
  xmlTextWriterEndElement(writer());//productList

  //Output transition state if one is specified (non-standard)
  if(pReact->GetTransitionState().get())
  {
    xmlTextWriterStartElementNS(writer(), prefix, C_TS, NULL);
    string id = AddMolToList(pReact->GetTransitionState(), OMols);
    if(list) //output a reference
    {
      xmlTextWriterStartElementNS(writer(), prefix, C_MOLECULE, NULL);
      xmlTextWriterWriteFormatAttribute(writer(), C_REF,"%s", id.c_str());
      xmlTextWriterEndElement(writer());//molecule
    }
    else
      //Write product in CML format
      pCMLFormat->WriteMolecule(OMols[id].get(), _pxmlConv);
    xmlTextWriterEndElement(writer());//me:transitionState
  }

  //Write reaction rate data, if there is any (non-standard)
  if(!_pxmlConv->IsOption("r"))
    WriteRateData(pReact, altprefix);

  xmlTextWriterEndElement(writer());//reaction

  if(pConv->IsLast()) //After all reactions converted
  {
    if(list)
    {
      xmlTextWriterEndElement(writer());//reactionList

      //output moleculeList
      xmlTextWriterStartElementNS(writer(), prefix, C_MOLECULELIST, NULL);
      MolMap::iterator mapitr;
      for(mapitr=OMols.begin();mapitr!=OMols.end();++mapitr)
      {
        if(strcmp((mapitr->second)->GetTitle(), "M"))
          pCMLFormat->WriteMolecule(mapitr->second.get(), _pxmlConv);
      }
      xmlTextWriterEndElement(writer());//moleculeList

      if(!_pxmlConv->IsOption("ReactionsNotStandalone"))
        xmlTextWriterEndElement(writer());//LISTWRAPPER

      xmlTextWriterEndDocument(writer());
      OutputToStream();

      //Send to real output stream with the moleculeList first
      string::size_type reaclistPos, mollistPos, footerPos;
      const string& s = ssout.str();
      reaclistPos = s.find("<reactionList");
      if(reaclistPos!=string::npos)
        mollistPos = s.find("<moleculeList",reaclistPos+1);
      footerPos = s.find("</cml");
      if(footerPos==string::npos)// cml tag may have been supressed
        footerPos = s.size();
      *_pOut << s.substr(0, reaclistPos)                      //header
            << s.substr(mollistPos, footerPos-mollistPos)    //moleculeList
            << s.substr(reaclistPos, mollistPos-reaclistPos) //reactionList
            << s.substr(footerPos) << endl;                  //footer
    }
    else if(_pxmlConv->GetOutputIndex()>1 && !_pxmlConv->IsOption("ReactionsNotStandalone"))
      xmlTextWriterEndElement(writer());//WRAPPER
    xmlTextWriterEndDocument(writer());
    OutputToStream();

    OMols.clear(); //clean up, delete molecules
  }
  _pxmlConv->RemoveOption("MolsNotStandalone",OBConversion::OUTOPTIONS);
  return true;
}

string CMLReactFormat::AddMolToList(obsharedptr<OBMol> spmol, MolMap& mmap)
{
  //Adds a molecule to the map
  string id = spmol->GetTitle();
  MolMap::iterator mapitr;
  if(id.empty())
  {
    //no id, so make one
    stringstream ssid;
    ssid << "m" << nextmol++;
    id = ssid.str();
      spmol->SetTitle(id);
    mmap[id] = spmol;
  }
  else
  {
    //If id is a filename with a path or an extension, which can
    //happen if no explicit title given to a molecule,
    // remove path and extension and use the filename as id
    string::size_type pos;
    pos = id.find_last_of("/\\:");
    if(pos!=string::npos)
      id.erase(0, pos+1);

    pos = id.rfind('.');
    if(pos!=string::npos)
      id.erase(pos);

    if(!isalpha(id[0])) //since ids have to start with a letter, add "id" to those that don't...
      id = "id" + id;
    spmol->SetTitle(id.c_str());//ensure that molecule title is the same as id

    mapitr = mmap.find(id);
    if(mapitr==mmap.end())
      //not in map; need to add
      mmap[id] = spmol;
    else
    {
      //already in map.
      //Get a molecule with the best bits of both old and new molecules and immediately make a shared_ ptr
      obsharedptr<OBMol> spnew(OBMoleculeFormat::MakeCombinedMolecule(mapitr->second.get(), spmol.get()));
      if(spnew)
      {
        spmol.swap(spnew);
        mapitr->second = spmol; //replace with new molecule
      //The OBMol originally in map will be deleted when local shared_ptr goes out of scope.
      }
    }
  }
  return id;
}



bool CMLReactFormat::WriteRateData(OBReaction* pReact, xmlChar* altprefix)
{
  OBRateData* pRD = static_cast<OBRateData*>(pReact->GetData(RateData));
  if(!pRD || pRD->GetRate(OBRateData::A)==0)
    return false; //nothing written

  static const xmlChar C_REF[]          = "ref";
  //following not strict cml
  static const xmlChar C_RATEPARAMS[]   = "rateParameters";
  static const xmlChar C_A[]            = "A";
  static const xmlChar C_N[]            = "n";
  static const xmlChar C_E[]            = "E";
  static const xmlChar C_LOA[]          = "loA";
  static const xmlChar C_LON[]          = "lon";
  static const xmlChar C_LOE[]          = "loE";
  static const xmlChar C_TROEPARAMS[]   = "troeParams";
  static const xmlChar C_REACTIONTYPE[] = "reactionType";
  static const xmlChar C_REVERSIBLE[]   = "reversible";
  static const xmlChar C_EFF[]          = "eff";

  string rtype("arrhenius");
  switch(pRD->ReactionType)
  {
  case OBRateData::TROE:
    rtype="troe"; break;
  case OBRateData::SRI:
    rtype="sri"; break;
  case OBRateData::LINDERMANN:
    rtype="lindermann"; break;
  case OBRateData::THREEBODY:
    rtype="threeBody";
  }
  xmlTextWriterStartElementNS(writer(), altprefix, C_RATEPARAMS, NULL);
  xmlTextWriterWriteFormatAttribute(writer(), C_REACTIONTYPE,"%s", rtype.c_str());
  if(pReact->IsReversible())
    xmlTextWriterWriteFormatAttribute(writer(), C_REVERSIBLE,"%s", "true");

  xmlTextWriterStartElementNS(writer(), altprefix, C_A, NULL);
  xmlTextWriterWriteFormatString(writer(),"%.3e", pRD->GetRate(OBRateData::A));
  xmlTextWriterEndElement(writer());//A

  xmlTextWriterStartElementNS(writer(), altprefix, C_N, NULL);
  xmlTextWriterWriteFormatString(writer(),"%g", pRD->GetRate(OBRateData::n));
  xmlTextWriterEndElement(writer());//n

  xmlTextWriterStartElementNS(writer(), altprefix, C_E, NULL);
  xmlTextWriterWriteFormatString(writer(),"%g", pRD->GetRate(OBRateData::E));
  xmlTextWriterEndElement(writer());//E

  switch(pRD->ReactionType)
  {
  case OBRateData::TROE:
    xmlTextWriterStartElementNS(writer(), altprefix, C_TROEPARAMS, NULL);
    xmlTextWriterWriteFormatString(writer(),"%g %g %g %g",
      pRD->GetTroeParam(0),pRD->GetTroeParam(1),pRD->GetTroeParam(2),pRD->GetTroeParam(3));
    xmlTextWriterEndElement(writer());
    //fallthrough
  case OBRateData::LINDERMANN:
    xmlTextWriterStartElementNS(writer(), altprefix, C_LOA, NULL);
    xmlTextWriterWriteFormatString(writer(),"%.3e", pRD->GetLoRate(OBRateData::A));
    xmlTextWriterEndElement(writer());//loA

    xmlTextWriterStartElementNS(writer(), altprefix, C_LON, NULL);
    xmlTextWriterWriteFormatString(writer(),"%g", pRD->GetLoRate(OBRateData::n));
    xmlTextWriterEndElement(writer());//lon

    xmlTextWriterStartElementNS(writer(), altprefix, C_LOE, NULL);
    xmlTextWriterWriteFormatString(writer(),"%g", pRD->GetLoRate(OBRateData::E));
    xmlTextWriterEndElement(writer());//loE
  //fallthrough
  case OBRateData::THREEBODY:
    string id;
    double Eff;
    while(pRD->GetNextEff(id,Eff))
    {
      xmlTextWriterStartElementNS(writer(), altprefix, C_EFF, NULL);
      xmlTextWriterWriteFormatAttribute(writer(), C_REF,"%s", id.c_str());
      xmlTextWriterWriteFormatString(writer(),"%g", Eff);
      xmlTextWriterEndElement(writer());//Eff
    }
  }
  xmlTextWriterEndElement(writer());//rateParams
  return true;
}
  void CMLReactFormat::WriteMetadataList(OBReaction& react)
  {
    static const xmlChar C_METADATALIST[] = "metadataList";
    static const xmlChar C_METADATA[]     = "metadata";
    static const xmlChar C_TITLE[]        = "title";
    static const xmlChar C_NAME[]         = "name";
    static const xmlChar C_CONTENT[]      = "content";

    string comment = react.GetComment();
    if(!comment.empty())
    {
      xmlTextWriterStartElement(writer(), C_METADATALIST);
      xmlTextWriterWriteAttributeNS(writer(),BAD_CAST "xmlns",BAD_CAST "dc",NULL,BAD_CAST "http://purl.org/dc/elements/1.1/");
      xmlTextWriterStartElement(writer(), C_METADATA);
      xmlTextWriterWriteAttribute(writer(), C_NAME, BAD_CAST "dc:description");
      xmlTextWriterWriteAttribute(writer(), C_CONTENT, BAD_CAST comment.c_str());
      xmlTextWriterEndElement(writer());
      xmlTextWriterEndElement(writer());
    }
   
  }

} //namespace OpenBabel

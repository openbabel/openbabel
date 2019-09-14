/**********************************************************************
Copyright (C) 2007 by Chris Morley

This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/
#include "openbabel/babelconfig.h"
#include <string>
#include <iomanip>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/elements.h>

#include "openbabel/obconversion.h"
#include "openbabel/reaction.h"

using namespace std;
//using std::tr1::shared_ptr;

namespace OpenBabel
{
  class SmiReactFormat : public OBFormat
  {
  public:
    //Register this format type ID
    SmiReactFormat()
    {
      OBConversion::RegisterFormat("rsmi",this);
    }

    virtual const char* Description()
    {
      return
        "Reaction SMILES format\n"
        "Write Options e.g. -xt\n"
        "  r radicals lower case eg ethyl is Cc\n"
        "\n";

    }

    virtual const char* GetMIMEType()
    { return "chemical/x-daylight-smiles"; }; // not right, need something else

    virtual const char* TargetClassDescription()
    {
      return OBReaction::ClassDescription();
    }

    const type_info& GetType()
    {
      return typeid(OBReaction*);
    }


    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool ReadMolecule(OBBase* pReact, OBConversion* pConv);
    virtual bool WriteMolecule(OBBase* pReact, OBConversion* pConv);

    ////////////////////////////////////////////////////
    /// The "Convert" interface functions
    virtual bool ReadChemObject(OBConversion* pConv)
    {
      //Makes a new OBReaction and new associated OBMols
      OBReaction* pReact = new OBReaction;
      bool ret=ReadMolecule(pReact,pConv); //call the "API" read function

      std::string auditMsg = "OpenBabel::Read reaction ";
      std::string description(Description());
      auditMsg += description.substr(0,description.find('\n'));
      obErrorLog.ThrowError(__FUNCTION__,
          auditMsg,
          obAuditMsg);

      if(ret) //Do transformation and return molecule
        return pConv->AddChemObject(pReact->DoTransformations(pConv->GetOptions(OBConversion::GENOPTIONS),pConv))!=0;
      else
      {
        pConv->AddChemObject(NULL);
        delete pReact;
        pReact=NULL;
        return false;
      }
    }

    virtual bool WriteChemObject(OBConversion* pConv)
    {
      //WriteChemObject() always deletes the object retrieved by GetChemObject
      OBBase* pOb = pConv->GetChemObject();
      OBReaction* pReact = dynamic_cast<OBReaction*>(pOb);
      if(pReact==NULL)
        return false;

      bool ret=false;
      ret=WriteMolecule(pReact,pConv);

      std::string auditMsg = "OpenBabel::Write reaction ";
      std::string description(Description());
      auditMsg += description.substr( 0, description.find('\n') );
      obErrorLog.ThrowError(__FUNCTION__,
          auditMsg,
          obAuditMsg);
      delete pOb;
      return ret;
    }

  };

  //Make an instance of the format class
  SmiReactFormat theSmiReactFormat;

  static bool IsNotEndChar(char t)
  {
    switch (t) {
    case '\0':
    case '\t':
    case ' ':
      return false;
    }
    return true;
  }

  /////////////////////////////////////////////////////////////////
  bool SmiReactFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {
    //It's really a reaction, not a molecule.
    //Doesn't make a new OBReaction object, but does make mew reactant and product OBMols
    OBReaction* pReact = pOb->CastAndClear<OBReaction>();

    istream &ifs = *pConv->GetInStream();
    OBConversion sconv; //Copy
    if(!sconv.SetInFormat("smi"))
    {
      obErrorLog.ThrowError(__FUNCTION__, "Smiles format needed but not found", obError);
      return false;
    }

    string ln, rsmiles, title, s;
    string::size_type pos, pos2;

    //Ignore lines that start with # or /
    while ((ifs && ifs.peek()=='#') || ifs.peek()=='/')
      if(!getline(ifs, ln))
        return false;

    //Get title
    if(!getline(ifs, ln))
      return false;
    pos = ln.find_first_of(" \t");
    if(pos!=string::npos)
    {
      rsmiles = ln.substr(0,pos);
      title = ln.substr(pos+1);
      Trim(title);
      pReact->SetTitle(title);
    }
    else
      rsmiles = ln;

    //Check for illegal characters
    pos = rsmiles.find_first_of(",<\"\'!^&_|{}");
    if(pos!=string::npos)
    {
      obErrorLog.ThrowError(__FUNCTION__,
          rsmiles + " contained a character '" + rsmiles[pos] + "' which is invalid in SMILES", obError);
      return false;
    }

    pos = rsmiles.find('>');
    if(pos==string::npos)
    {
      obErrorLog.ThrowError(__FUNCTION__, "No > in reaction", obError);
      return false;
    }

    vector<OBMol> mols;
    vector<OBMol>::iterator itr;

    //Extract reactants and split into individual molecules
    OBMol jreactants;
    s = rsmiles.substr(0,pos);
    if(pos > 0 && !sconv.ReadString(&jreactants, s))
    {
      obErrorLog.ThrowError(__FUNCTION__, "Cannot read reactant", obError);
      return false;
    }
    mols = jreactants.Separate();
    for(itr=mols.begin();itr!=mols.end();++itr)
      pReact->AddReactant(obsharedptr<OBMol>(new OBMol(*itr)));

    pos2 = rsmiles.find('>', pos+1);
    if(pos2==string::npos)
    {
      obErrorLog.ThrowError(__FUNCTION__, "Only one > in reaction", obError);
      return false;
    }

    //Extract agent (not split into separate molecules)
    if(pos2-pos>1)
    {
      OBMol* pAgent = new OBMol;
      s = rsmiles.substr(pos+1,pos2-pos-1);
      if(!sconv.ReadString(pAgent, s))
      {
        obErrorLog.ThrowError(__FUNCTION__, "Cannot read agent", obError);
        delete pAgent;
        return false;
      }
      pReact->AddAgent(obsharedptr<OBMol>(pAgent));
    }

    //Extract products and split into separate molecules
    OBMol jproducts;
    s = rsmiles.substr(pos2+1);
    if(IsNotEndChar(s[0]) && !sconv.ReadString(&jproducts, s))
    {
      obErrorLog.ThrowError(__FUNCTION__, "Cannot read product", obError);
      return false;
    }
    mols.clear();
    mols = jproducts.Separate();
    for(itr=mols.begin();itr!=mols.end();++itr)
      pReact->AddProduct(obsharedptr<OBMol>(new OBMol(*itr)));

    return true;
  }


  /////////////////////////////////////////////////////////////////
  bool SmiReactFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    //It's really a reaction, not a molecule.
    //Cast output object to the class type need, i.e. OBReaction
    OBReaction* pReact = dynamic_cast<OBReaction*>(pOb);
    if(pReact==NULL)
      return false;
    ostream &ofs = *pConv->GetOutStream();


    OBFormat* pSmiFormat = OBConversion::FindFormat("SMI");
    if(!pSmiFormat)
      return false;
    pConv->AddOption("smilesonly",OBConversion::OUTOPTIONS);//supresses title and new line
    pConv->AddOption("c",OBConversion::OUTOPTIONS);//output atom classes if available

    OBMol jReactants;
    for(int i=0;i<pReact->NumReactants();++i)
      jReactants += *(pReact->GetReactant(i));

    if(!pSmiFormat->WriteMolecule(&jReactants, pConv))
      return false;

    ofs << '>';

    OBMol jAgents;
    for (int i = 0; i<pReact->NumAgents(); ++i)
      jAgents += *(pReact->GetAgent(i));

    if(!pSmiFormat->WriteMolecule(&jAgents, pConv))
      return false;

    ofs << '>';

    OBMol jProducts;
    for(int i=0;i<pReact->NumProducts();++i)
      jProducts += *(pReact->GetProduct(i));

    if(!pSmiFormat->WriteMolecule(&jProducts, pConv))
      return false;

    if(!pReact->GetTitle().empty())
      ofs << '\t' << pReact->GetTitle();

    ofs << endl;

    return true;
  }

} //namespace

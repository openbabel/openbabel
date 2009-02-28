/**********************************************************************
Copyright (C) 2007 by Chris Morley
 
This file is part of the Open Babel project.
For more information, see <http://openbabel.sourceforge.net/>
 
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
#include "openbabel/mol.h"
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
  { return "chemical/x-mdl-rxn"; };

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
      return pConv->AddChemObject(pReact)!=0;
    else
        pConv->AddChemObject(NULL);
    return false;
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

/////////////////////////////////////////////////////////////////
bool SmiReactFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
{
  //It's really a reaction, not a molecule.
  //Doesn't make a new OBReactionObject, but does make mew reactant and product OBMols
  OBReaction* pReact = dynamic_cast<OBReaction*>(pOb);

  istream &ifs = *pConv->GetInStream();

  OBFormat* pSmiFormat = OBConversion::FindFormat("SMI");
  if(!pSmiFormat)
    return false;

  //Read reactant
  shared_ptr<OBMol> spReactant(new OBMol);
  if(!pSmiFormat->ReadMolecule(spReactant.get(), pConv))
  {
    obErrorLog.ThrowError(__FUNCTION__, "Cannot read reactant", obError);
    return false;
  }
  pReact->AddReactant(spReactant);

  char ch;
  if(!ifs.get(ch) || ch!='>')
  {
    obErrorLog.ThrowError(__FUNCTION__, "No > in reaction", obError);
    return false;
  }
  
  //Read >> characters and possibly an agent molecule between them
  if(ifs.get(ch) && ch!='>')
  {
    //there is an agent
    ifs.unget();
    shared_ptr<OBMol> spAgent(new OBMol);
    if(!pSmiFormat->ReadMolecule(spAgent.get(), pConv))
    {
      obErrorLog.ThrowError(__FUNCTION__, "Error in agent molecule", obError);
      return false;
    }
    pReact->AddAgent(spAgent);
    if(!ifs.get(ch) || ch!='>')
    {
      obErrorLog.ThrowError(__FUNCTION__, "The second > is missing", obError);
      return false;
    }
  }

  //Read product
  shared_ptr<OBMol> spProduct(new OBMol);
  if(!pSmiFormat->ReadMolecule(spProduct.get(), pConv))
  {
    obErrorLog.ThrowError(__FUNCTION__, "Cannot read product", obError);
    return false;
  }
  pReact->AddProduct(spProduct);

  //The comment at the end of the line ends up as the title of the product molecule
  string comment = spProduct->GetTitle();
  spProduct->SetTitle("");
  pReact->SetComment(comment);

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

  if(pReact->NumReactants()!=1 || pReact->NumProducts()>1)
    obErrorLog.ThrowError(__FUNCTION__,
      "ReactionSMILES format is only for a single reactant and product", obError);

  shared_ptr<OBMol> spReactant = pReact->GetReactant(0);
  if(!spReactant.get() || spReactant->NumAtoms()==0)
    obErrorLog.ThrowError(__FUNCTION__,"Missing or empty reactant", obWarning);

  if(!pSmiFormat->WriteMolecule(spReactant.get(), pConv))
    return false;

  ofs << '>';

  shared_ptr<OBMol> spAgent = pReact->GetAgent();
  if(spAgent.get())
    if(!pSmiFormat->WriteMolecule(spAgent.get(), pConv))
      return false;

  ofs << '>';

  shared_ptr<OBMol> spProduct = pReact->GetProduct(0);
  if(!spProduct.get() || spProduct->NumAtoms()==0)
    obErrorLog.ThrowError(__FUNCTION__,"Missing or empty product", obWarning);
  if(!pSmiFormat->WriteMolecule(spProduct.get(), pConv))
    return false;

  if(!pReact->GetComment().empty())
    ofs << '\t' << pReact->GetComment();

  ofs << endl;

  return true;
}

} //namespace

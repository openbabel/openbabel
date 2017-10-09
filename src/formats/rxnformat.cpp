/**********************************************************************
Copyright (C) 2004 by Chris Morley

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
#include <typeinfo>
#include "openbabel/mol.h"
#include "openbabel/obconversion.h"
#include "openbabel/reaction.h"

using namespace std;
//using std::tr1::shared_ptr;

namespace OpenBabel
{
class RXNFormat : public OBFormat
{
public:
  //Register this format type ID
  RXNFormat()
  {
      OBConversion::RegisterFormat("rxn",this);
  }

  virtual const char* Description()
  {
      return
        "MDL RXN format\n"
        "The MDL reaction format is used to store information on chemical reactions.\n\n"
        "Output Options, e.g. -xA\n"
        " A  output in Alias form, e.g. Ph, if present\n"
        " G <option> how to handle any agents present\n\n"
        "            One of the following options should be specifed:\n\n"
        "            - agent - Treat as an agent (default). Note that some programs\n"
        "                      may not read agents in RXN files.\n"
        "            - reactant - Treat any agent as a reactant\n"
        "            - product - Treat any agent as a product\n"
        "            - ignore - Ignore any agent\n"
        "            - both - Treat as both a reactant and a product\n\n";
  };

  virtual const char* GetMIMEType()
  { return "chemical/x-mdl-rxn"; };

  virtual const char* TargetClassDescription()
  {
      return OBReaction::ClassDescription();
  };

  const type_info& GetType()
  {
    return typeid(OBReaction*);
  };


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
};

  virtual bool WriteChemObject(OBConversion* pConv)
  {
    //WriteChemObject() always deletes the object retrieved by GetChemObject
    //For RXN NO LONGER deletes the associated molecules which are handled by a smart pointer
    //Cast to the class type need, e.g. OBMol
    OBBase* pOb=pConv->GetChemObject();
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
  };
};

//Make an instance of the format class
RXNFormat theRXNFormat;

/////////////////////////////////////////////////////////////////
bool RXNFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
{
    //It's really a reaction, not a molecule.
    //Doesn't make a new OBReactionObject, but does make mew reactant and product OBMols
   OBReaction* pReact = pOb->CastAndClear<OBReaction>();

    OBFormat* pMolFormat = pConv->FindFormat("MOL");
    if(pMolFormat==NULL || !pReact)
        return false;

    //	OBConversion MolConv(*pConv); //new copy to use to read associated MOL

    istream &ifs = *pConv->GetInStream();
    string ln;
    // When MDLFormat reads the last product it may also read and discard
    // the line with $RXN for the next reaction. But it then sets $RXNread option.
    if(pConv->IsOption("$RXNread"))
      pConv->RemoveOption("$RXNread", OBConversion::OUTOPTIONS);
    else
    {
      if (!getline(ifs,ln))
        return(false);
      if(Trim(ln).find("$RXN")!=0)
        return false; //Has to start with $RXN
    }
    if (!getline(ifs,ln))
      return(false); //reaction title
    pReact->SetTitle(Trim(ln));

    if (!getline(ifs,ln))
      return false; //creator
    if (!getline(ifs, ln))
      return(false); //comment
    pReact->SetComment(Trim(ln));

    int nReactants, nProducts, i;
    ifs >> setw(3) >> nReactants >> setw(3) >> nProducts >> ws;
    if(!ifs) return false;

    if(nReactants + nProducts)
    {
      //Read the first $MOL. The others are read at the end of the previous MOL
      if(!getline(ifs, ln))
        return false;
      if(Trim(ln).find("$MOL")==string::npos)
        return false;
    }

    OBMol* pmol;

    for(i=0;i<nReactants;i++)
    {
      //Read a MOL file	using the same OBConversion object but with a different format
      pmol=new OBMol;
      if(!pMolFormat->ReadMolecule(pmol,pConv))
        obErrorLog.ThrowError(__FUNCTION__, "Failed to read a reactant", obWarning);
      else
      {
        obsharedptr<OBMol> p(pmol);
        pReact->AddReactant(p);
      }
    }

    for(i=0;i<nProducts;i++)
    {
      //Read a MOL file
      pmol=new OBMol;
      if(!pMolFormat->ReadMolecule(pmol,pConv))
        obErrorLog.ThrowError(__FUNCTION__, "Failed to read a product", obWarning);
      else
      {
        //        pReact->products.push_back(pmol);
        obsharedptr<OBMol> p(pmol);
        pReact->AddProduct(p);
      }
    }

    return(true);
}

enum HandleAgent {
  AS_AGENT, IGNORE, AS_REACT, AS_PROD, BOTH_REACT_AND_PROD
};

static HandleAgent ReadAgentOption(const char* t)
{
  if (!t)
    return AS_AGENT; // default
  switch(t[0]) {
  case 'a':
    if (t[1]=='g' && t[2]=='e' && t[3]=='n' && t[4]=='t' && t[5]=='\0')
      return AS_AGENT;
    break;
  case 'i':
    if (t[1]=='g' && t[2]=='n' && t[3]=='o' && t[4]=='r' && t[5]=='e' && t[6]=='\0')
      return IGNORE;
    break;
  case 'r':
    if (t[1]=='e' && t[2]=='a' && t[3]=='c' && t[4]=='t' && t[5]=='a' && t[6]=='n' && t[7]=='t' && t[8]=='\0')
      return AS_REACT;
    break;
  case 'p':
    if (t[1]=='r' && t[2]=='o' && t[3]=='d' && t[4]=='u' && t[5]=='c' && t[6]=='t' && t[7]=='\0')
      return AS_PROD;
    break;
  case 'b':
    if (t[1]=='o' && t[2]=='t' && t[3]=='h' && t[4]=='\0')
      return BOTH_REACT_AND_PROD;
    break;
  }
  return AS_AGENT;
}

static void WriteMolFile(OBMol* pmol, OBConversion* pconv, OBFormat* pformat)
{
  ostream &ofs = *pconv->GetOutStream();
  ofs << "$MOL" << '\n';
  pformat->WriteMolecule(pmol, pconv);
}

/////////////////////////////////////////////////////////////////
bool RXNFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
    //It's really a reaction, not a molecule.
    //Cast output object to the class type need, i.e. OBReaction
    OBReaction* pReact = dynamic_cast<OBReaction*>(pOb);
    if(pReact==NULL)
        return false;

    pConv->AddOption("no$$$$",OBConversion::OUTOPTIONS);

    OBFormat* pMolFormat = pConv->FindFormat("MOL");
    if(pMolFormat==NULL)
    {
      obErrorLog.ThrowError(__FUNCTION__, "MDL MOL format not available", obError);
        return false;
    }

    HandleAgent handleagent = ReadAgentOption(pConv->IsOption("G"));
    obsharedptr<OBMol> agent = pReact->GetAgent();
    bool agentInReactants;
    if (agent && (handleagent==BOTH_REACT_AND_PROD || handleagent==AS_REACT))
      agentInReactants = true;
    else
      agentInReactants = false;
    bool agentInProducts;
    if (agent && (handleagent==BOTH_REACT_AND_PROD || handleagent==AS_PROD))
      agentInProducts = true;
    else
      agentInProducts = false;

    ostream &ofs = *pConv->GetOutStream();

    ofs << "$RXN" << '\n';
    ofs << pReact->GetTitle() << '\n';
    ofs << "      OpenBabel" << '\n';
    ofs << pReact->GetComment() << '\n';

    ofs << setw(3);
    if (agentInReactants)
      ofs << pReact->NumReactants() + 1;
    else
      ofs << pReact->NumReactants();
    ofs << setw(3);
    if (agentInProducts)
      ofs << pReact->NumProducts() + 1;
    else
      ofs << pReact->NumProducts();
    if (agent && handleagent==AS_AGENT)
      ofs << setw(3) << 1;
    ofs << '\n';

    // Write reactants
    for(unsigned int i=0; i<pReact->NumReactants(); i++)
      WriteMolFile(pReact->GetReactant(i).get(), pConv, pMolFormat);
    if (agentInReactants)
      WriteMolFile(pReact->GetAgent().get(), pConv, pMolFormat);

    // Write products
    for(unsigned int i=0; i<pReact->NumProducts(); i++)
      WriteMolFile(pReact->GetProduct(i).get(), pConv, pMolFormat);
    if (agentInProducts)
      WriteMolFile(pReact->GetAgent().get(), pConv, pMolFormat);

    // Write agent out (if treating AS_AGENT)
    if(agent && handleagent==AS_AGENT)
      WriteMolFile(pReact->GetAgent().get(), pConv, pMolFormat);

    return true;
}

} //namespace

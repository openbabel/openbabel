/**********************************************************************
Copyright (C) 2004 by Chris Morley
Copyright (C) 2018 by Noel M. O'Boyle

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
#include <openbabel/babelconfig.h>
#include <string>
#include <iomanip>
#include <typeinfo>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/elements.h>
#include <openbabel/generic.h>

#include <openbabel/obconversion.h>
#include <openbabel/obmolecformat.h>
#include <openbabel/reactionfacade.h>

using namespace std;

namespace OpenBabel
{
class RXNFormat : public OBMoleculeFormat
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
      return OBMol::ClassDescription();
  }


  ////////////////////////////////////////////////////
  /// The "API" interface functions
  virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
  virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

};

//Make an instance of the format class
RXNFormat theRXNFormat;

static bool ParseComponent(const char* t, unsigned int *ans)
{
  const char *p = t;
  while (*p == ' ')
    p++;
  while (p - t < 3) {
    if (*p < '0' || *p > '9')
      return false;
    *ans *= 10;
    *ans += *p - '0';
    p++;
  }
  return true;
}

/////////////////////////////////////////////////////////////////
bool RXNFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
{
    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if (pmol == NULL)
      return false;

    OBFormat* pMolFormat = pConv->FindFormat("MOL");
    if (pMolFormat==NULL)
      return false;

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
      return false; //reaction title
    pmol->SetTitle(Trim(ln));

    if (!getline(ifs,ln))
      return false; //creator
    if (!getline(ifs, ln))
      return false; //comment
    // Originally the comment was added to the reaction via:
    //     pmol->SetComment(Trim(ln));

    if (!getline(ifs, ln))
      return false; // num reactants, products, and optionally agents

    unsigned int nReactants = 0, nProducts = 0, nAgents = 0;
    bool ok = ParseComponent(ln.c_str() + 0, &nReactants);
    if (!ok)
      return false;
    ok = ParseComponent(ln.c_str() + 3, &nProducts);
    if (!ok)
      return false;
    if (ln[6] != '\0') { // optional agents
      ok = ParseComponent(ln.c_str() + 6, &nAgents);
      if (!ok)
        return false;
    }

    if(nReactants + nProducts + nAgents)
    {
      //Read the first $MOL. The others are read at the end of the previous MOL
      if(!getline(ifs, ln))
        return false;
      if(Trim(ln).find("$MOL")==string::npos)
        return false;
    }

    OBReactionFacade rxnfacade(pmol);

    // Note: If we supported it, we could read each of the rxn components directly
    // into the returned OBMol instead of having to do a copy. Unfortunately,
    // this isn't possible at the moment (MOL format will need some work first).
    // Here is some example code to do it:
    //
    //unsigned int old_numatoms = 0;
    //unsigned int compid = 1;
    //for (int i = 0; i<nReactants; i++)
    //{
    //  //Read a MOL file	using the same OBConversion object but with a different format
    //  if (!pMolFormat->ReadMolecule(pmol, pConv))
    //    obErrorLog.ThrowError(__FUNCTION__, "Failed to read a reactant", obWarning);
    //  unsigned int numatoms = pmol->NumAtoms();
    //  for (unsigned int idx = old_numatoms + 1; idx <= numatoms; ++idx) {
    //    OBAtom* atom = pmol->GetAtom(idx);
    //    rxnfacade.SetRole(atom, REACTANT);
    //    rxnfacade.SetComponentId(atom, compid);
    //  }
    //  old_numatoms = numatoms;
    //  compid++;
    //}

    const char* type[3] = {"a reactant", "a product", "an agent"};
    OBReactionRole role;
    unsigned int num_components;
    for(unsigned int N=0; N<3; N++) {
      switch(N) {
      case 0:
        role = REACTANT;
        num_components = nReactants;
        break;
      case 1:
        role = PRODUCT;
        num_components = nProducts;
        break;
      case 2:
        role = AGENT;
        num_components = nAgents;
        break;
      }
      for (int i=0; i<num_components; i++)
      {
        //Read a MOL file	using the same OBConversion object but with a different format
        OBMol mol;
        if (!pMolFormat->ReadMolecule(&mol, pConv)) {
          std::string error = "Failed to read ";
          error += type[N];
          obErrorLog.ThrowError(__FUNCTION__, error, obWarning);
          continue;
        }
        if (mol.NumAtoms() == 0) {
          OBAtom* dummy = mol.NewAtom(); // Treat the empty OBMol as having a single dummy atom
          OBPairData *pd = new OBPairData();
          pd->SetAttribute("rxndummy");
          pd->SetValue("");
          pd->SetOrigin(fileformatInput);
          dummy->SetData(pd);
        }

        rxnfacade.AddComponent(&mol, role);
      }
    }

    pmol->SetIsReaction();
    return true;
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
  // Treat a dummy atom with "rxndummy" as the empty file
  if (pmol->NumAtoms() == 1) {
    OBAtom *atm = pmol->GetFirstAtom();
    if (atm->GetAtomicNum() == 0 && atm->HasData("rxndummy"))
      pmol->DeleteAtom(atm);
  }
  pformat->WriteMolecule(pmol, pconv);
}

static void WriteAgents(OBMol& mol, OBReactionFacade& rxnfacade, OBConversion* pconv, OBFormat* pformat)
{
  for(unsigned int i=0; i<rxnfacade.NumComponents(AGENT); i++) {
    mol.Clear();
    rxnfacade.GetComponent(&mol, AGENT, i);
    WriteMolFile(&mol, pconv, pformat);
  }
}

/////////////////////////////////////////////////////////////////
bool RXNFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if (pmol == NULL || !pmol->IsReaction())
      return false;

    pConv->AddOption("no$$$$",OBConversion::OUTOPTIONS);

    OBFormat* pMolFormat = pConv->FindFormat("MOL");
    if(pMolFormat==NULL)
    {
      obErrorLog.ThrowError(__FUNCTION__, "MDL MOL format not available", obError);
        return false;
    }

    OBReactionFacade rxnfacade(pmol);

    HandleAgent handleagent = ReadAgentOption(pConv->IsOption("G"));
    bool hasAgent = rxnfacade.NumComponents(AGENT) > 0;
    bool agentInReactants, agentInProducts;
    if (hasAgent && (handleagent==BOTH_REACT_AND_PROD || handleagent==AS_REACT))
      agentInReactants = true;
    else
      agentInReactants = false;
    if (hasAgent && (handleagent==BOTH_REACT_AND_PROD || handleagent==AS_PROD))
      agentInProducts = true;
    else
      agentInProducts = false;

    ostream &ofs = *pConv->GetOutStream();

    ofs << "$RXN" << '\n';
    ofs << pmol->GetTitle() << '\n';
    ofs << "      OpenBabel" << '\n';
    //ofs << pReact->GetComment() << '\n';
    ofs << "\n";

    ofs << setw(3);
    if (agentInReactants)
      ofs << rxnfacade.NumComponents(REACTANT) + rxnfacade.NumComponents(AGENT);
    else
      ofs << rxnfacade.NumComponents(REACTANT);
    ofs << setw(3);
    if (agentInProducts)
      ofs << rxnfacade.NumComponents(PRODUCT) + rxnfacade.NumComponents(AGENT);
    else
      ofs << rxnfacade.NumComponents(PRODUCT);
    if (hasAgent && handleagent==AS_AGENT)
      ofs << setw(3) << rxnfacade.NumComponents(AGENT);
    ofs << '\n';

    // Write reactants
    OBMol mol;
    for(unsigned int i=0; i<rxnfacade.NumComponents(REACTANT); i++) {
      mol.Clear();
      rxnfacade.GetComponent(&mol, REACTANT, i);
      WriteMolFile(&mol, pConv, pMolFormat);
    }
    if (agentInReactants)
      WriteAgents(mol, rxnfacade, pConv, pMolFormat);

    // Write products
    for(unsigned int i=0; i<rxnfacade.NumComponents(PRODUCT); i++) {
      mol.Clear();
      rxnfacade.GetComponent(&mol, PRODUCT, i);
      WriteMolFile(&mol, pConv, pMolFormat);
    }
    if (agentInProducts)
      WriteAgents(mol, rxnfacade, pConv, pMolFormat);

    // Write agent out (if treating AS_AGENT)
    if(hasAgent && handleagent==AS_AGENT)
      WriteAgents(mol, rxnfacade, pConv, pMolFormat);
    
    return true;
}

} //namespace

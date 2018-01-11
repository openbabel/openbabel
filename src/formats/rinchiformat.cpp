/**********************************************************************
Copyright (C) 2017 by Noel M. O'Boyle

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
#include <algorithm>
#include "openbabel/mol.h"
#include "openbabel/obconversion.h"
#include "openbabel/reaction.h"

using namespace std;

#define RINCHI_VERSION_STRING "RInChI=1.00.1S/"

namespace OpenBabel
{
  class ReactionInChIFormat : public OBFormat
  {
  public:
    //Register this format type ID
    ReactionInChIFormat()
    {
      OBConversion::RegisterFormat("rinchi",this);
    }

    virtual const char* Description()
    {
      return
        "RInChI\n"
        "The Reaction InChI (or RInChI) is intended to be a unique\n"
        "string that describes a reaction. This may be useful for\n"
        "indexing and searching reaction databases. As with the InChI\n"
        "it is recommended that you always keep the original reaction\n"
        "information and use the RInChI in addition.\n\n"

        "The RInChI format is a hierarchical, layered description of a\n"
        "reaction with different levels based on the Standard InChI\n"
        "representation of each structural component participating in\n"
        "the reaction.\n\n"

        "Write Options e.g. -xe\n"
        "  e Treat this reaction as an equilibrium reaction\n"
        "    Layer 5 of the generated RInChI will have /d=\n"
        "\n";
    }

    virtual const char* TargetClassDescription()
    {
      return OBReaction::ClassDescription();
    }

    virtual unsigned int Flags()
    {
      return NOTREADABLE;
    };

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
      return true;
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

      delete pOb;
      return ret;
    }

  };

  //Make an instance of the format class
  ReactionInChIFormat theReactionInChIFormat;

  /////////////////////////////////////////////////////////////////
  bool ReactionInChIFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {
    return true;
  }

  // Trim off the "InChI=1S/" and the trailing \n
  static std::string TrimInChI(const char *inchi)
  {
    std::string trimmed;
    const char *p = inchi + 9;
    while (true) {
      trimmed += *p;
      p++;
      if (*p == '\n' || *p == '\0')
        break;
    }
    return trimmed;
  }

  /////////////////////////////////////////////////////////////////
  bool ReactionInChIFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    //It's really a reaction, not a molecule.
    //Cast output object to the class type need, i.e. OBReaction
    OBReaction* pReact = dynamic_cast<OBReaction*>(pOb);
    if (pReact == NULL)
      return false;
    ostream &ofs = *pConv->GetOutStream();

    OBFormat* pInChIFormat = OBConversion::FindFormat("inchi");
    if (!pInChIFormat)
      return false;

    bool isEquilibrium = pReact->IsReversible() || pConv->IsOption("e");

    OBConversion inchiconv;
    inchiconv.SetOutFormat(pInChIFormat);
    stringstream ss;
    inchiconv.SetOutStream(&ss);

    //if (pReact->NumReactants() == 0 && pReact->NumProducts() == 0 && pReact->NumAgents() == 0) {
    //  // Special-case the empty RInChI
    //  ofs << RINCHI_VERSION_STRING << "/d+\n";
    //  return true;
    //}

#define REACTANTS 0
#define PRODUCTS 1
#define AGENTS 2

    std::vector<std::vector<std::string> > inchis(3);
    unsigned int nonInchi[3] = { 0, 0, 0 };
    bool hasNonInchi = false;
    for (int part = 0; part <= 2; ++part) {
      unsigned int N;
      switch (part) {
      case REACTANTS: N = pReact->NumReactants(); break;
      case PRODUCTS: N = pReact->NumProducts(); break;
      case AGENTS: N = pReact->NumAgents(); break;
      }
      for (unsigned int i = 0; i < N; ++i) {
        OBMol* mol;
        switch (part) {
        case REACTANTS: mol = &*(pReact->GetReactant(i)); break;
        case PRODUCTS: mol = &*(pReact->GetProduct(i)); break;
        case AGENTS: mol = &*(pReact->GetAgent(i)); break;
        }
        if (mol->NumAtoms() == 0) {
          nonInchi[part]++;
          hasNonInchi = true;
        }
        else {
          bool ok = inchiconv.Write(mol);
          if (!ok) {
            nonInchi[part]++;
            hasNonInchi = true;
          }
          else {
            string inchi = ss.str();
            if (strncmp(inchi.c_str(), "InChI=1S/", 9) != 0)
              return false;
            inchis[part].push_back(TrimInChI(inchi.c_str()));
          }
          ss.str("");
        }
      }
    }

    std::sort(inchis[REACTANTS].begin(), inchis[REACTANTS].end());
    std::sort(inchis[PRODUCTS].begin(), inchis[PRODUCTS].end());
    std::sort(inchis[AGENTS].begin(), inchis[AGENTS].end());

    std::string reactants_string = "";
    const int rsize = inchis[REACTANTS].size();
    for (int i = 0; i < rsize; ++i) {
      if (i > 0)
        reactants_string += '!';
      reactants_string += inchis[REACTANTS][i];
    }
    std::string products_string = "";
    const int psize = inchis[PRODUCTS].size();
    for (int i = 0; i < psize; ++i) {
      if (i > 0)
        products_string += '!';
      products_string += inchis[PRODUCTS][i];
    }

    bool reactants_first = reactants_string <= products_string;

    ofs << RINCHI_VERSION_STRING;
    if (rsize > 0 || psize > 0 || !inchis[AGENTS].empty()) {
      ofs << (reactants_first ? reactants_string : products_string);
      ofs << "<>";
      ofs << (reactants_first ? products_string : reactants_string);
      if (!inchis[AGENTS].empty()) {
        ofs << "<>";
        for (std::vector<std::string>::const_iterator vit = inchis[AGENTS].begin(); vit != inchis[AGENTS].end(); ++vit) {
          if (vit != inchis[AGENTS].begin())
            ofs << '!';
          ofs << *vit;
        }
      }
    }
    ofs << "/d";
    if (isEquilibrium)
      ofs << '=';
    else
      ofs << (reactants_first ? '+' : '-');
    if (hasNonInchi) {
      ofs << "/u" << (reactants_first ? nonInchi[REACTANTS] : nonInchi[PRODUCTS]) << '-'
        << (reactants_first ? nonInchi[PRODUCTS] : nonInchi[REACTANTS]) << '-'
        << nonInchi[AGENTS];
    }

    ofs << '\n';
    return true;
  }

} //namespace

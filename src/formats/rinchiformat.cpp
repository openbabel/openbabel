/**********************************************************************
Copyright (C) 2017,2018 by Noel M. O'Boyle

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
#include <algorithm>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/elements.h>

#include <openbabel/obconversion.h>
#include <openbabel/reactionfacade.h>
#include <openbabel/obmolecformat.h>

using namespace std;

#define RINCHI_VERSION_STRING "RInChI=1.00.1S/"

namespace OpenBabel
{
  class ReactionInChIFormat : public OBMoleculeFormat
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
        "The Reaction InChI\n"
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
      return OBMol::ClassDescription();
    }

    virtual unsigned int Flags()
    {
      return NOTREADABLE;
    }

    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

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
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if (pmol == NULL || !pmol->IsReaction())
      return false;
    ostream &ofs = *pConv->GetOutStream();

    OBFormat* pInChIFormat = OBConversion::FindFormat("inchi");
    if (!pInChIFormat)
      return false;

    bool isEquilibrium = pConv->IsOption("e");

    OBConversion inchiconv;
    inchiconv.SetOutFormat(pInChIFormat);
    stringstream ss;
    inchiconv.SetOutStream(&ss);

#define M_REACTANTS 0
#define M_PRODUCTS 1
#define M_AGENTS 2

    OBReactionFacade facade(pmol);

    std::vector<std::vector<std::string> > inchis(3);
    unsigned int nonInchi[3] = { 0, 0, 0 };
    bool hasNonInchi = false;
    OBMol mol;
    for (int part = 0; part <= 2; ++part) {
      unsigned int N;
      switch (part) {
      case M_REACTANTS: N = facade.NumComponents(REACTANT); break;
      case M_PRODUCTS: N = facade.NumComponents(PRODUCT); break;
      case M_AGENTS: N = facade.NumComponents(AGENT); break;
      }
      for (unsigned int i = 0; i < N; ++i) {
        mol.Clear();
        switch (part) {
        case M_REACTANTS: facade.GetComponent(&mol, REACTANT, i); break;
        case M_PRODUCTS: facade.GetComponent(&mol, PRODUCT, i); break;
        case M_AGENTS: facade.GetComponent(&mol, AGENT, i); break;
        }
        if (mol.NumAtoms() == 1 && mol.GetFirstAtom()->GetAtomicNum() == 0) {
          // This represents an unknown component
          nonInchi[part]++;
          hasNonInchi = true;
        }
        else {
          bool ok = inchiconv.Write(&mol);
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

    std::sort(inchis[M_REACTANTS].begin(), inchis[M_REACTANTS].end());
    std::sort(inchis[M_PRODUCTS].begin(), inchis[M_PRODUCTS].end());
    std::sort(inchis[M_AGENTS].begin(), inchis[M_AGENTS].end());

    std::string reactants_string = "";
    const int rsize = inchis[M_REACTANTS].size();
    for (int i = 0; i < rsize; ++i) {
      if (i > 0)
        reactants_string += '!';
      reactants_string += inchis[M_REACTANTS][i];
    }
    std::string products_string = "";
    const int psize = inchis[M_PRODUCTS].size();
    for (int i = 0; i < psize; ++i) {
      if (i > 0)
        products_string += '!';
      products_string += inchis[M_PRODUCTS][i];
    }

    bool reactants_first = reactants_string <= products_string;

    ofs << RINCHI_VERSION_STRING;
    if (rsize > 0 || psize > 0 || !inchis[M_AGENTS].empty()) {
      ofs << (reactants_first ? reactants_string : products_string);
      ofs << "<>";
      ofs << (reactants_first ? products_string : reactants_string);
      if (!inchis[M_AGENTS].empty()) {
        ofs << "<>";
        for (std::vector<std::string>::const_iterator vit = inchis[M_AGENTS].begin(); vit != inchis[M_AGENTS].end(); ++vit) {
          if (vit != inchis[M_AGENTS].begin())
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
      ofs << "/u" << (reactants_first ? nonInchi[M_REACTANTS] : nonInchi[M_PRODUCTS]) << '-'
        << (reactants_first ? nonInchi[M_PRODUCTS] : nonInchi[M_REACTANTS]) << '-'
        << nonInchi[M_AGENTS];
    }

    ofs << '\n';
    return true;
  }

} //namespace

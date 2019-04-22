/**********************************************************************
Copyright (C) 2012 by Tim Vandermeersch

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
#include <openbabel/elements.h>
#include <openbabel/bond.h>
#include <openbabel/obiter.h>

#include <openbabel/typer.h>
#include <openbabel/stereo/tetrahedral.h>
#include <openbabel/stereo/cistrans.h>

#include "smiley.h"

using namespace std;

namespace OpenBabel
{

  template<typename T>
  void print_vector(const std::string &label, const std::vector<T> &v)
  {
    std::cout << label << ": ";
    for (std::size_t i = 0; i < v.size(); ++i)
      std::cout << v[i] << " ";
    std::cout << std::endl;
  }

  /**
   *
   * Smiley callback.
   *
   */
  struct OpenBabelCallback : public Smiley::CallbackBase
  {
    enum UpDown {
      IsNotUpDown,
      IsUp,
      IsDown
    };

    OpenBabelCallback(OBMol *pmol) : mol(pmol)
    {
    }

    void clear()
    {
      // prepare for new SMILES
    }

    void addAtom(int element, bool aromatic, int isotope, int hCount, int charge, int atomClass)
    {
      // invoked when an atom is completely parsed

      // element: 0, 1, 2, ... (0 is '*' in SMILES)
      // aromatic: true if the atom was lower case c, n, ...
      // isotope: -1, 0, 1, 2, 3, 4 ... (-1 means no isotope specified and is not the same as 0)
      // hCount: -1, 0, 1, 2, ..., 9 (-1 means default valence and is only for unbracketed atoms)
      // charge: -9, -8, ..., -1, 0, 1, ..., 8, 9
      // atomClass: 0, 1, 2, 3, 4, ... (0 means no atom class, specified atom classes should start from 1)

      OBAtom *atom = mol->NewAtom();
      atom->SetAtomicNum(element);

      indices.push_back(mol->NumAtoms());

      if (aromatic)
        atom->SetAromatic();
      //else if (hCount == -1)
      //  atom->ForceImplH();

      if (hCount > -1) {
        if (hCount == 0)
          atom->SetSpinMultiplicity(2);

        for (int i = 0; i < hCount; ++i) {
          OBAtom *hydrogen = mol->NewAtom();
          hydrogen->SetAtomicNum(1);
          mol->AddBond(atom->GetIdx(), hydrogen->GetIdx(), 1);
          upDown.push_back(IsNotUpDown);
        }
      }

      if (isotope > 0)
        atom->SetIsotope(isotope);

      atom->SetFormalCharge(charge);
    }

    void addBond(int source, int target, int order, bool isUp, bool isDown)
    {
      // invoked for each bond once both of it's atoms have been added by
      // calling addAtom(). This ensures that the bond's atom indexes are always valid.

      // source: source atom index starting from 0 (order from addAtom() calls)
      // target: target atom index starting from 0 (order from addAtom() calls)
      // order: 1, 2, 3, 4, 5 (5 means aromatic)
      // isUp: true if bond is single order up bond '/'
      // isDown: true if bond is single order down bond '\'

      //std::cout << "addBond(" << source << ", " << target << ")" << std::endl;

      if (isDown)
        upDown.push_back(IsDown);
      else if (isUp)
        upDown.push_back(IsUp);
      else
        upDown.push_back(IsNotUpDown);

      /*
      if (isUp || isDown) {
        std::cout << "addBond(" << source << ", " << target << ")" << std::endl;
        std::cout << "isUp: " << isUp << ", isDown: " << isDown << std::endl;
      }
      */

      mol->AddBond(indices[source], indices[target], order);
      if (order == 5)
        mol->GetBond(mol->NumBonds() - 1)->SetAromatic();
    }

    void setChiral(int index, Smiley::Chirality chirality, const std::vector<int> &chiralNbrs)
    {
      // invoked at the end of parsing for each chiral atom

      // index: atom index starting from 0
      // chirality: Clockwise, AntiClockwise, TH1, AL2, SP3, TB14, OH26, ...
      // chiralNbrs: atom indices of neighbors, size 4 for TH, AL and SP, size 5 for TB and 6 for OH

      //print_vector("chiralNbrs", chiralNbrs);

      unsigned long center = indices[index] - 1;
      unsigned long from = indices[chiralNbrs[0]] - 1;
      std::vector<unsigned long> refs(chiralNbrs.size() - 1);
      for (std::size_t i = 0; i < refs.size(); ++i)
        if (chiralNbrs[i + 1] == Smiley::implicitHydrogen())
          refs[i] = OBStereo::ImplicitRef;
        else
          refs[i] = indices[chiralNbrs[i + 1]] - 1;

      //std::cout << "center: " << center << std::endl;
      //std::cout << "from: " << from << std::endl;
      //print_vector("refs", refs);

      switch (chirality) {
        case Smiley::Clockwise:
          switch (chiralNbrs.size()) {
            case 4:
              OBTetrahedralStereo *stereo = new OBTetrahedralStereo(mol);
              stereo->SetConfig(OBTetrahedralStereo::Config(center, from, refs));
              mol->SetData(stereo);
              break;
          }
          break;
        case Smiley::AntiClockwise:
          switch (chiralNbrs.size()) {
            case 4:
              OBTetrahedralStereo *stereo = new OBTetrahedralStereo(mol);
              stereo->SetConfig(OBTetrahedralStereo::Config(center, from, refs, OBStereo::AntiClockwise));
              mol->SetData(stereo);
              break;
          }
          break;
      }
    }

    OBMol *mol;
    std::vector<UpDown> upDown;
    std::vector<int> indices;
  };

  /**
   *
   * OpenBabel format.
   *
   */
  class SmileyFormat : public OBMoleculeFormat
  {
    public:
      SmileyFormat()
      {
        OBConversion::RegisterFormat("smy", this);
      }

      virtual const char* Description() //required
      {
        return "SMILES format using Smiley parser\n\n"

"The Smiley parser presents an alternative to the standard SMILES parser\n"
"(:ref:`SMILES_format`). It was written to be strictly compatible with the\n"
"OpenSMILES standard (http://opensmiles.org). In comparison, the standard\n"
"parser is more forgiving to erroneous input, and also supports some extensions\n"
"such as for radicals.\n\n"

"In addition, the Smiley parser returns detailed error messages when problems\n"
"arise parsing or validating the SMILES, whereas the standard parser seldom\n"
"describes the specific problem. For a detailed description of the OpenSMILES\n"
"semantics, the specification should be consulted. In addition to syntactical\n"
"and grammatical correctness, the Smiley parser also verifies some basic\n"
"semantics.\n"
"\n"
"Here are some examples of the errors reported::\n"
"\n"
"   SyntaxError: Bracket atom expression contains invalid trailing characters.\n"
"   F.FB(F)F.[NH2+251][C@@H](CP(c1ccccc1)c1ccccc1)C(C)(C)C 31586112\n"
"                  ^^\n"
"   SyntaxError: Unmatched branch opening.\n"
"   CC(CC\n"
"     ^^^\n"
"   SyntaxError: Unmatched branch closing.\n"
"   CC)CC\n"
"   ^^^\n"
"   SemanticsError: Unmatched ring bond.\n"
"   C1CCC\n"
"   ^\n"
"   SemanticsError: Conflicing ring bonds.\n"
"   C-1CCCCC=1\n"
"\n"
"Hydrogen with Hydrogen Count\n"
"~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
"Hydrogen atoms can not have a hydrogen count. Hydrogen bound to a hydrogen\n"
"atom should be specified by two bracket atom expressions.\n"
"\n"
"Examples::\n"
"\n"
"  [HH]        invalid\n"
"  [HH1]       invalid (same as [HH]\n"
"  [HH3]       invalid\n"
"  [HH0]       valid (same as [H])\n"
"  [H][H]      valid\n"
"\n"
"Unmatched Ring Bond\n"
"~~~~~~~~~~~~~~~~~~~\n"
"Report unmatched ring bonds.\n"
"\n"
"Example::\n"
"\n"
"  C1CCC\n"
"\n"
"Conflicting Ring Bonds\n"
"~~~~~~~~~~~~~~~~~~~~~~\n"
"When the bond type for ring bonds are explicitly specified at both ends,\n"
"these should be the same.\n"
"\n"
"Example::\n"
"\n"
"  C-1CCCCCC=1\n"
"\n"
"Invalid Ring Bonds\n"
"~~~~~~~~~~~~~~~~~~\n"
"There are two types of invalid ring bonds. The first is when two atoms both\n"
"have the same two ring bonds. This would mean adding a parallel edge in the\n"
"graph which is not allowed. The second type is similar but results in a\n"
"self-loop by having a ring bond number twice.\n"
"\n"
"Examples::\n"
"\n"
"  C12CCCC12      parallel bond\n"
"  C11            self-loop bond\n"
"\n"
"Invalid Chiral Valence\n"
"~~~~~~~~~~~~~~~~~~~~~~\n"
"When an atom is specified as being chiral, it should have the correct\n"
"number of neighboring atoms (possibly including an implicit H inside the\n"
"bracket.\n"
"\n"
"The valid valences are::\n"
"\n"
"  Tetrahedral (TH)          : 4\n"
"  Allene (AL)               : 4 (*)\n"
"  Square Planar (SP)        : 4\n"
"  Trigonal Bypiramidal (TB) : 5\n"
"  Octahedral(OH)            : 6\n"
"\n"
"  (*) The chiral atom has only 2 bonds but the neighbor's neighbors are\n"
"      counted: NC(Br)=[C@AL1]=C(F)I\n"
"\n"
"Invalid Chiral Hydrogen Count\n"
"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
"Chiral atoms can only have one hydrogen in their bracket since multiple\n"
"hydrogens would make them not chiral.\n"
"\n"
"Example::\n"
"\n"
"  C[C@H2]F\n\n";
      }

      virtual const char* SpecificationURL()
      {
        return "http://opensmiles.org";
      }

      virtual const char* GetMIMEType()
      {
        return "chemical/x-daylight-smiles";
      }

      virtual unsigned int Flags()
      {
        return NOTWRITABLE;
      }

      virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);

    private:
      void CreateCisTrans(OBMol *mol, const std::vector<OpenBabelCallback::UpDown> &upDown);
      bool AssignNbrAtoms(const std::vector<OpenBabelCallback::UpDown> &upDown, OBAtom *atom,
          unsigned long &above, unsigned long &below);

  };

  ////////////////////////////////////////////////////

  //Make an instance of the format class
  SmileyFormat theSmileyFormat;

  /////////////////////////////////////////////////////////////////




  bool SmileyFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if (!pmol)
      return false;

    // get the input stream
    istream& ifs = *pConv->GetInStream();

    // read the smiles string
    std::string smiles;
    std::getline(ifs, smiles);

    // extract title
    std::size_t space_pos = smiles.find(" ");
    std::size_t tab_pos = smiles.find("\t");
    if (space_pos != std::string::npos && tab_pos != std::string::npos)
      space_pos = std::min(space_pos, tab_pos);
    else if (tab_pos != std::string::npos)
      space_pos = tab_pos;

    if (space_pos != std::string::npos) {
      while (space_pos < smiles.size() && (smiles[space_pos] == ' ' || smiles[space_pos] == '\t'))
        ++space_pos;
      pmol->SetTitle(smiles.substr(space_pos).c_str());
    }

    pmol->BeginModify();
    pmol->SetDimension(0);

    // create callback and parser object
    OpenBabelCallback callback(pmol);
    Smiley::Parser<OpenBabelCallback> parser(callback);

    try {
      parser.parse(smiles);
    } catch (Smiley::Exception &e) {
      if (e.type() == Smiley::Exception::SyntaxError)
        std::cerr << "Syntax";
      else
        std::cerr << "Semantics";
      std::cerr << "Error: " << e.what() << "." << std::endl;
      std::cerr << smiles << std::endl;
      for (std::size_t i = 0; i < e.pos(); ++i)
        std::cerr << " ";
      for (std::size_t i = 0; i < e.length(); ++i)
        std::cerr << "^";
      std::cerr << std::endl;
    }

    pmol->EndModify();

    // handle aromaticity
    pmol->SetAromaticPerceived();

    // create cis/trans stereo objects
    CreateCisTrans(pmol, callback.upDown);
    StereoFrom0D(pmol);

    return true;
  }

  bool SmileyFormat::AssignNbrAtoms(const std::vector<OpenBabelCallback::UpDown> &upDown,
      OBAtom *atom, unsigned long &aboveId, unsigned long &belowId)
  {
    OBAtom *above = 0;
    OBAtom *below = 0;
    OBAtom *unspecified = 0;

    FOR_BONDS_OF_ATOM (bond, atom) {
      if (!bond->IsAromatic() && bond->GetBondOrder() == 2)
        continue;

      OBAtom *nbr = bond->GetNbrAtom(atom);
      //std::cout << "atom: " << atom->GetIndex() << std::endl;
      //std::cout << "nbr: " << nbr->GetIndex() << std::endl;

      switch (upDown[bond->GetIdx()]) {
        case OpenBabelCallback::IsUp:
          if (nbr->GetIndex() < atom->GetIndex() && bond->GetBeginAtomIdx() < bond->GetEndAtomIdx()) {
            //std::cout << "below: " << nbr->GetIndex() << std::endl;
            if (below)
              return false;
            below = nbr;
          } else {
            //std::cout << "above: " << nbr->GetIndex() << std::endl;
            if (above)
              return false;
            above = nbr;
          }
          break;
        case OpenBabelCallback::IsDown:
          if (nbr->GetIndex() < atom->GetIndex() && bond->GetBeginAtomIdx() < bond->GetEndAtomIdx()) {
            if (above)
              return false;
            above = nbr;
          } else {
            if (below)
              return false;
            below = nbr;
          }
          break;
        case OpenBabelCallback::IsNotUpDown:
          //std::cout << "unspecified: " << nbr->GetIndex() << std::endl;
          unspecified = nbr;
          break;
      }
    }

    // at least 1 bond should be specified
    if (!above && !below)
      return true;

    if (above && unspecified)
      below = unspecified;
    else if (below && unspecified)
      above = unspecified;

    aboveId = above ? above->GetId() : OBStereo::ImplicitRef;
    belowId = below ? below->GetId() : OBStereo::ImplicitRef;

    return true;
  }

  void SmileyFormat::CreateCisTrans(OBMol *mol, const std::vector<OpenBabelCallback::UpDown> &upDown)
  {
    FOR_BONDS_OF_MOL (doubleBond, mol) {
      if (doubleBond->GetBondOrder() != 2 || doubleBond->IsAromatic())
        continue;

      OBAtom *source = doubleBond->GetBeginAtom();
      OBAtom *target = doubleBond->GetEndAtom();

      //std::cout << "double bond: " << source->GetIndex() << " " << target->GetIndex() << std::endl;

      // Check that both atoms on the double bond have at least one
      // other neighbor, but not more than two other neighbors;
      int v1 = source->GetExplicitDegree();
      int v2 = target->GetExplicitDegree();
      if (v1 < 2 || v1 > 3 || v2 < 2 || v2 > 3)
        continue;

      unsigned long aboveSource = OBStereo::ImplicitRef;
      unsigned long belowSource = OBStereo::ImplicitRef;
      if (!AssignNbrAtoms(upDown, source, aboveSource, belowSource)) {
        std::cerr << "Invalid cis/trans specification" << std::endl;
        continue;
      }

      if (aboveSource == OBStereo::ImplicitRef && belowSource == OBStereo::ImplicitRef)
        continue;

      unsigned long aboveTarget = OBStereo::ImplicitRef;
      unsigned long belowTarget = OBStereo::ImplicitRef;
      if (!AssignNbrAtoms(upDown, target, aboveTarget, belowTarget)) {
        std::cerr << "Invalid cis/trans specification" << std::endl;
        continue;
      }

      if (aboveTarget == OBStereo::ImplicitRef && belowTarget == OBStereo::ImplicitRef)
        continue;

      //std::cout << "refs: " << aboveSource << " " << aboveTarget << " " << belowTarget << " " << belowSource << std::endl;

      OBCisTransStereo *stereo = new OBCisTransStereo(mol);
      stereo->SetConfig(OBCisTransStereo::Config(source->GetId(), target->GetId(),
            OBStereo::MakeRefs(aboveSource, belowSource, belowTarget, aboveTarget),
            OBStereo::ShapeU));

      mol->SetData(stereo);
    }
  }

} //namespace OpenBabel


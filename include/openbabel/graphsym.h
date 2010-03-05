/**********************************************************************
graphsym.h - Class for handling graph symmetry.

To determine copyright, please analyse the Subversion commit log.

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

#ifndef OB_GRAPHSYM_H
#define OB_GRAPHSYM_H

#include <openbabel/babelconfig.h>
#include <vector>

#ifndef EXTERN
#  define EXTERN extern
#endif

using namespace std;

namespace OpenBabel {

  class OBBitVec;
  class OBMol;
  class OBAtom;
  class OBBond;
  class OBMol;

  class OBAPI OBGraphSym {
    private:
      OBBitVec* _frag_atoms;
      OBMol* _pmol;

      static bool CompareUnsigned(const unsigned int &a,const unsigned int &b);
      static bool ComparePairFirst(const std::pair<OBAtom*,unsigned int> &a,const std::pair<OBAtom*,unsigned int> &b);
      static bool ComparePairSecond(const std::pair<OBAtom*,unsigned int> &a,const std::pair<OBAtom*,unsigned int> &b);
      static bool CompareBondPairSecond(const std::pair<OBBond*,unsigned int> &a,const std::pair<OBBond*,unsigned int> &b);

      unsigned int GetValence(OBAtom *atom);
      unsigned int GetHvyValence(OBAtom *atom);
      unsigned int GetHvyBondSum(OBAtom *atom);

      void FindRingAtoms(OBBitVec &ring_atoms);
      void CreateNewClassVector(std::vector<std::pair<OBAtom*,unsigned int> > &vp1,
                                std::vector<std::pair<OBAtom*,unsigned int> > &vp2);
      void GetGIVector(std::vector<unsigned int> &vid);
      bool GetGTDVector(std::vector<int> &gtd);
      void CountAndRenumberClasses(std::vector<std::pair<OBAtom*,unsigned int> > &vp, unsigned int &count);
      int ExtendInvariants(std::vector<std::pair<OBAtom*, unsigned int> > &symmetry_classes, bool breakChiralTies);
                           
      int CalculateSymmetry(std::vector<unsigned int> &symmetry_classes, bool breakChiralTies);
      void BreakChiralTies(vector<pair<OBAtom*, unsigned int> > &atom_sym_classes);
    public:
      //! Constructor
      OBGraphSym(OBMol* pmol, OBBitVec* frag_atoms = NULL);
      //! Destructor
      virtual ~OBGraphSym();

      static const unsigned int NoSymmetryClass;
      
      /**
       * Calculate the symmetry classes for the molecule. The result will be 
       * stored in @p symmetry_classes.
       *
       * The results in @p symmetry_classes will be ordered by symmetry 
       * classes. Use the OBAtom* pointer in the std::pair to match the atoms 
       * with the right symmetry classes.
       *
       * @return The number of symmetry classes.
       */
      int GetSymmetry(vector<unsigned int> &symmetry_classes, bool breakChiralTies = true);
      void ClearSymmetry();
      /**
       * Calculate the canonical labels for the molecule. The result will be 
       * stored in @p canonical_labels.
       *
       * FIXME
       *
       * @return FIXME
       */
      void CanonicalLabels(vector<unsigned int> &symmetry_classes);
    };

      
      
      
      //&OBBitVec GetFragment();
      //SetFragment(&OBBitVec);

} // namespace OpenBabel

//! \file graphsym.h
//! \brief XXXX

  #endif // OB_GRAPHSYM_H

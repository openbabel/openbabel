/**********************************************************************
builder.h - OBBuilder class.

Copyright (C) 2007-2008 by Tim Vandermeersch
                           <tim.vandermeersch@gmail.com>

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

#ifndef OB_BUILDER_H
#define OB_BUILDER_H

#include <vector>
#include <string>
#include <map>

#include <list>
#include <set>

#include <openbabel/stereo/stereo.h>

namespace OpenBabel
{
  class OBMol;
  class OBAtom;
  class OBSmartsPattern;
  class vector3;
  class OBBitVec;

  //! \class OBBuilder builder.h <openbabel/builder.h>
  //! \brief Class to build 3D structures
  class OBAPI OBBuilder {
    public:

      OBBuilder(): _keeprings(false) {}

      ///@name Call the build algorithm
      //@{
      /*! The mol object contains all connectivity information (atomic numbers, bonds, bond orders, ..)
       *  but no 3D coordinates. Build generates these coordinates and assigns them.
       *  \param mol Molecule with the connectivity (from SMILES for example). The coordinates are also
       *         changed in this mol.
       *  \param stereoWarnings Warn if the stereochemistry is incorrect (default is true)
       */
    bool Build(OBMol &mol, bool stereoWarnings = true);
      //@}

      ///@name Setup build parameters
      //@{
      /*! If the molecule already contains 3D coordinates, if you set KeepRings to true it will use
       *  retain the 3D coordinates of the rings. By default KeepRings is false, and ring conformations
       *  are obtained by lookup in a library of ring conformers. However, since the ring conformer library
       *  is not exhaustive, if the ring system is not found in the library, the resulting 3D structure can
       *  be poor, and require geometry optimisation before it is reasonable. If your starting point is
       *  a 3D structure, you can set KeepRings to true, and the conformation will be taken from the input.
       *  The remaining (acyclic) bonds will still all be built by the builder.
       */
      void SetKeepRings() { _keeprings = true; }
      void UnsetKeepRings() { _keeprings = false; }
      //@}


      //! Used by LoadFragments to check for invalid (all zero coordinates) fragments
      void AddRingFragment(OBSmartsPattern *sp, const std::vector<vector3> &coords);
      //! Load fragment info from file, if is it has not already been done
      void LoadFragments();
      std::vector<vector3> GetFragmentCoord(std::string smiles);

      /*! Get the position for a new neighbour on atom.  Returns
       * non-finite vector if there is no reasonable location.
       *  \param atom Atom for which we want a new neighbour location.
       *  \returns The position for the new atom.
       */
      static vector3 GetNewBondVector(OBAtom *atom);
      static vector3 GetNewBondVector(OBAtom *atom, double length);

      /*! Atoms a and b are part of two fragments that are not connected in mol.
       *  Connect will translate and rotate the fragment that contains b so that
       *  a and b are separated by a bond. This bond is also added.
       *  \param mol The molecule to be modified
       *  \param a Index for atom in fragment that should not be rotated.
       *  \param b Index for atom in fragment that should be rotated.
       *  \param newpos Direction for new bond between a and b
       *  \param bondOrder Bond order of the new bond between a and b.
       *  \returns true if successful or fails when failed (most likely cause
       *  for failing: a and b are in the same fragment, they are connected)
       */
      static bool Connect(OBMol &mol, int a, int b, vector3 &newpos, int bondOrder = 1);
      /*! Atoms a and b are part of two fragments that are not connected in mol.
       *  Connect will translate and rotate the fragment that contains b so that
       *  a and b are separated by a bond. This bond is also added.
       *  \param mol The molecule to be modified
       *  \param a Index for atom in fragment that should not be rotated.
       *  \param b Index for atom in fragment that should be rotated.
       *  \param bondOrder Bond order of the new bond bewtween a and b.
       *  \returns true if successful or fails when failed (most likely cause
       *  for failing: a and b are in the same fragment, they are connected)
       */
      static bool Connect(OBMol &mol, int a, int b, int bondOrder = 1);
      /*! Swap group b, bonded to a with group d, bonded to c. The bonds a-b and b-c cannot be
       *  part of a ring. Atoms a and b will not be moved. Atoms b, d and their connected atoms
       *  (after deleting bonds ab and cd) will be translated/rotated.
       *
       *  Example:
       *  \code
       *    \ /                            /
       *     b                            d
       *      \     /     Swap(a,b,c,d)    \     /
       *       a---x          ---->         a---x
       *      /     \     /                /     \     /
       *     x       c---d                x       c---b
       *                                               \
       *  \endcode
       *
       *
       *  This function can also be used to invert chiral centers if a and c are the same atom.
       *
       *  Example
       *  \code
       *     1                        3
       *     |      Swap(C,1,C,3)     |
       *  2>-C-<3      ----->      2>-C-<1
       *     |                        |
       *     4                        4
       *  \endcode
       */
      static bool Swap(OBMol &mol, int a, int b, int c, int d);
      /*! Atoms a and b must be bonded and this bond cannot be part of a ring. The bond will
       *  be broken and the smiles fragment will be inserted bewteen the two remaining fragments.
       *  The fragment that contains a will not be translated or rotated. Parameters c and d are
       *  the index in the smiles to which atoms a and b will be connected respectivly.
       *
       */
      //bool Insert(OBMol &mol, int a, int b, std::string smiles, int c, int d);
      /*! Correct double bond stereochemistry
       *
       * \returns Success or failure
       */
      static bool CorrectStereoBonds(OBMol &mol);
      /*! Correct stereochemistry at tetrahedral atoms with at least two non-ring
       * bonds. It also works for spiro atoms.
       *
       * \returns Success or failure
       */
      static bool CorrectStereoAtoms(OBMol &mol, bool warn = true);
      /*! Does this atom connect two rings which are not otherwise connected?
      */
      static bool IsSpiroAtom(unsigned long atomId, OBMol &mol);
      /*! Get the fragment to which this atom belongs.
       *  \param atom Atom in the fragment.
       *  \returns The OBBitVec defining the fragment to which a belongs.
       */
      static OBBitVec GetFragment(OBAtom *atom);
      static void AddNbrs(OBBitVec &fragment, OBAtom *atom);

    private:
      //! used to hold the fragments loaded in the constructor
      //static std::map<std::string, double> _torsion;
      static std::vector<std::string> _rigid_fragments;
      static std::vector<std::pair<OBSmartsPattern*, std::vector<vector3> > > _ring_fragments;
      static std::map<std::string, int> _rigid_fragments_index;
      static std::map<std::string, std::vector<vector3> > _rigid_fragments_cache;
      //! Connect a ring fragment to an already matched fragment. Currently only
      //  supports the case where the fragments overlap at a spiro atom only.
      static void ConnectFrags(OBMol &mol, OBMol &workmol, std::vector<int> match, std::vector<vector3> coords,
                               std::vector<int> pivot);
      //! Rotate one of the spiro rings 180 degrees
      static void FlipSpiro(OBMol &mol, int idx);
      static bool FixRingStereo(std::vector<std::pair<OBStereo::Ref, bool> > atomIds,
                                OBMol &mol, OBStereo::Refs &unfixedcenters);
      static void AddRingNbrs(OBBitVec &fragment, OBAtom *atom, OBMol &mol);
      static bool SwapWithVector(OBMol &mol, int a, int b, int c, const vector3 &newlocation);
      bool _keeprings;
  }; // class OBBuilder

}// namespace OpenBabel

#endif   // OB_BUILDER_H

//! \file builder.h
//! \brief Class to build 3D structures

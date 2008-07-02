/**********************************************************************
builder.h - OBBuilder class.
 
Copyright (C) 2007-2008 by Tim Vandermeersch 
                           <tim.vandermeersch@gmail.com>
 
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

#ifndef OB_BUILDER_H
#define OB_BUILDER_H

#include <vector>
#include <string>
#include <map>

#include <list>
#include <set>
#include <openbabel/base.h>
#include <openbabel/mol.h>

namespace OpenBabel
{
  //! \class OBBuilder
  //! \brief Class to build 3D structures
  class OBAPI OBBuilder {
    public:

      //! Load fragment info from file, if is it has not already been done
      void LoadFragments();

      /*! Get the position for a new neighbour on atom.
       *  \param atom Atom for which we want a new neighbour location.
       *  \returns The position for the new atom.
       */  
      static vector3 GetNewBondVector(OBAtom *atom);
      static vector3 GetNewBondVector(OBAtom *atom, double length);
      /*! The mol object contains all connectivity information (atomic numbers, bonds, bond orders, ..) 
       *  but no 3D coordinates. Build generates these coordinates and assigns them.
       *  \param mol Molecule with the connectivity (from smiles for example). The coordinates are also
       *         changed in this mol.
       */    
      bool Build(OBMol &mol);
      /*! Atoms a and b are part of two fragments that are not connected in mol.
       *  Connect will translate and rotate the fragment that contains b so that
       *  a and b are seperated by a bond. This bond is also added.
       *  \param mol The molecule to be modified
       *  \param a Index for atom in fragment that should not be rotated.
       *  \param b Index for atom in fragment that should be rotated.
       *  \param newpos Direction for new bond between a and b
       *  \param bondOrder Bond order of the new bond between a and b.
       *  \returns true if succesful or fails when failed (most likely cause 
       *  for failing: a and b are in the same fragment, they are connected)
       */
      static bool Connect(OBMol &mol, int a, int b, vector3 &newpos, int bondOrder = 1);
      /*! Atoms a and b are part of two fragments that are not connected in mol.
       *  Connect will translate and rotate the fragment that contains b so that
       *  a and b are seperated by a bond. This bond is also added.
       *  \param mol The molecule to be modified
       *  \param a Index for atom in fragment that should not be rotated.
       *  \param b Index for atom in fragment that should be rotated.
       *  \param bondOrder Bond order of the new bond bewtween a and b.
       *  \returns true if succesfull or fails when failed (most likely cause 
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
      /*! Currently only corrects double bond chemistry comming from smiles. (OBBond::IsUp() / OBBond::IsDown())
       */ 
      static void CorrectStereoBonds(OBMol &mol);
      /*! Currently only corrects atom chirality comming from smiles. (OBAtom::IsClockwize() / OBBond::IsAntiClockwise())
       */ 
      static void CorrectStereoAtoms(OBMol &mol);
      /*! Get the fragment to which this atom belongs.
       *  \param atom Atom in the fragment.
       *  \returns The OBBitVec defining the fragment to which a belongs.
       */
      static OBBitVec GetFragment(OBAtom *atom);
      static void AddNbrs(OBBitVec &fragment, OBAtom *atom);
 
    private:
      //! used to hold the fragments loaded in the constructor
      static std::vector<std::pair<OBSmartsPattern*, std::vector<vector3> > > _fragments;
  }; // class OBBuilder

}// namespace OpenBabel

#endif   // OB_BUILDER_H

//! \file builder.h
//! \brief Class to build 3D structures

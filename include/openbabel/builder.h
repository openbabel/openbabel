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
      //! constructor
      OBBuilder();
      //! destructor
      ~OBBuilder();
      /*! Get the position for a new neighbour on atom.
       *  \param atom Atom for which we want a new neighbour location.
       *  \returns The position for the new atom.
       */  
      vector3 GetNewBondVector(OBAtom *atom);
      /*! The mol object contains all connectivity information (atomic numbers, bonds, bond orders, ..) 
       *  but no 3D coordinates. Build generates these coordinates and assigns them.
       *  \param mol Molecule with the connectivity (from smiles for example). The coordinates are also
       *         changed in this mol.
       */    
      bool Build(OBMol &mol);
      /*! Atoms a and b are part of two fragments that are not connected in mol.
       *  Connect will translate and rotate the fragment that contains b so that
       *  a and b are seperated by a bond. This bond is also added.
       *  \param a Index for atom in fragment that should not be rotated.
       *  \param b Index for atom in fragment that should be rotated.
       *  \param bondOrder Bond order of the new bond bewtween a and b.
       *  \returns true if succesfull or fails when failed (most likely cause 
       *  for failing: a and b are in the same fragment, they are connected)
       */
      bool Connect(OBMol &mol, int a, int b, int bondOrder = 1);
    private:
      /*! Atoms a and b are part of two fragments that are not connected in _workMol.
       *  Connect will translate and rotate the fragment containing b so that
       *  a and b are seperated by a bond. This bond is also added.
       *  \param a Index for atom in fragment that should not be rotated.
       *  \param b Index for atom in fragment that should be rotated.
       *  \returns true if succesfull or fails when failed (most likely cause 
       *  for failing: a and b are in the same fragment, they are connected)
       */
      bool Connect(OBAtom *a, OBAtom *b);
      /*! Get the fragment to which atom with index belongs.
       *  \param index Atom index.
       *  \returns The OBBitVec defining the fragment to which a belongs.
       */
      OBBitVec GetFragment(int index);
      void AddNbrs(OBBitVec &fragment, OBAtom *atom);
      //! used to hold the fragments loaded in the constructor
      std::vector<std::pair<OBSmartsPattern*, std::vector<vector3> > > _fragments;
      //! used to hold the molecule we're working on 
      OBMol _workMol;
  }; // class OBBuilder

}// namespace OpenBabel

#endif   // OB_BUILDER_H

//! \file builder.h
//! \brief Class to build 3D structures

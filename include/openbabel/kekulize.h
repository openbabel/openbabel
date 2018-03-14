/**********************************************************************
Copyright (C) 2017 Noel M. O'Boyle

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
#ifndef OB_KEKULIZE_H
#define OB_KEKULIZE_H

namespace OpenBabel
{  /** 
   \brief Kekulize a molecule by assigning bond orders of 1 or 2 to aromatic bonds

   Some file formats describe bond orders as aromatic. Such bonds require kekulization
   before the molecule is returned by the reader. Normally, a user should never need
   to call this function themselves.

   This function takes an OBMol which has atoms and bonds marked as aromatic,
   aromatic bonds whose bond orders have been set to single, and aromaticity set as perceived.
   The function assumes that atoms joined by aromatic bonds have been marked as aromatic
   (if they are so intended).

   The purpose of the function is to set the bond orders of the aromatic bonds to either 1 or 2
   in such a way that the valencies of all of the aromatic atoms are satisfied. Failure to do
   this will result in one or more atoms having unsatisfied valences, indicated by a radical.
   Such a failure can only occur if an atom is incorrectly marked as aromatic, or is correctly
   marked as aromatic but has incorrect valence (e.g. 'n' instead of '[nH]' in SMILES).

   \return Whether kekulization was successful
   **/
  
  OBAPI bool OBKekulize(OBMol *mol);
}

#endif //OB_KEKULIZE_H

//! \file kekulize.h
//! \brief Functions relating to kekulization

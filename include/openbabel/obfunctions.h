/**********************************************************************
obfunctions.h - Various global functions

Copyright (C) 2017 by Noel O'Boyle

This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.
***********************************************************************/

#ifndef OB_FUNCTIONS_H
#define OB_FUNCTIONS_H

#include <openbabel/atom.h>

namespace OpenBabel
{
  /**
  \brief Return the size of the smallest ring in which a bond appears

  This function returns the size of the smallest ring in which a bond appears. The
  search is bounded by the specified bound. A value of 0 is returned if the bond
  is not in a ring or if no ring is found of size less than or equal to the bound.

  Note that alternative algorithms may be more appropriate if you wish to calculate
  this value for all atoms in a molecule.

  \return The size of the smallest ring, or 0
  **/

  OBAPI unsigned int OBBondGetSmallestRingSize(OBBond *bond, unsigned int bound);

  /**
   \brief Return the typical valence of an atom of a particular element

   This function returns the typical valence of an atom given its element, current
   valence (that is, the current sum of the bond orders of its bonds) and formal
   charge.

   This is typically used on atoms that are missing hydrogens, to decide how many
   implicit hydrogens should be assigned (should one have to guess). For example,
   the value 3 is returned for a positively charged carbon with no attached atoms.

   \return A value for the typical valence
   **/
  	
  OBAPI unsigned int GetTypicalValence(unsigned int element, unsigned int bosum, int charge);
    /**
   \brief Assign implicit hydrogens to an OBAtom based on typical valences

   This function uses the return value of GetTypicalValence to determine how many
   implicit hydrogens to assign to an OBAtom. Note that most file formats describe
   exactly how many hydrogens are present, and do not require this function.
   **/
  
  OBAPI void OBAtomAssignTypicalImplicitHydrogens(OBAtom* atom);

} // end namespace OpenBabel

#endif

//! \file obfunctions.h
//! \brief A collection of global functions

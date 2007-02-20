/* -*-C++-*-

**********************************************************************
Copyright (C) 2005-2006, eMolecules, Inc. (www.emolecules.com)
 
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************

+======================================================================
| FILE:         canon.h
| AUTHOR:       Craig A. James
| DESCRIPTION:
|       Declarations for canon.cpp
+======================================================================
*/

// Return vector is indexed from zero, corresponds to "atom->GetIdx()-1"
namespace OpenBabel {

void CanonicalLabels(OBMol *pmol,
                     OBBitVec &frag_atoms,
                     std::vector<unsigned int> &symmetry_classes,
                     std::vector<unsigned int> &canonical_labels);

} // namespace OpenBabel

//! \file canon.h
//! \brief Canonical numbering of SMILES, molecules and fragments

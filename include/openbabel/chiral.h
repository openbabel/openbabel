/**********************************************************************
chiral.h - Detect chiral atoms and molecules.
 
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2005 by Geoffrey R. Hutchison
 
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

#ifndef OB_CHIRAL_H
#define OB_CHIRAL_H

#include <vector>
#include <openbabel/matrix.h>

namespace OpenBabel
{
OBAPI void GraphPotentials(OBMol &mol, std::vector<double> &pot);
OBAPI void construct_g_matrix(OBMol &mol, std::vector<std::vector<double> > &m);
OBAPI void construct_c_matrix(OBMol &mol, std::vector<std::vector<double > > &m);

//! Calculate the signed volume for an atom.
OBAPI double CalcSignedVolume(OBMol &mol,OBAtom*,bool ReZeroZ=true);

OBAPI double signed_volume(const vector3 &a, const vector3 &b, const vector3 &c, const vector3 &d);
OBAPI void GetChirality(OBMol &mol, std::vector<int> &chirality);

//! Calculates parity of a vector of 4 items
OBAPI int GetParity4Ref(std::vector<unsigned int> pref);
OBAPI bool CorrectChirality(OBMol &mol, OBAtom *atm, atomreftype i=input, atomreftype o=output);
}

#endif // OB_CHIRAL_H

//! \file chiral.h
//! \brief Detect chiral atoms and molecules.

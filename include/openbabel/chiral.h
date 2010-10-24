/**********************************************************************
chiral.h - Detect chiral atoms and molecules.

Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2005 by Geoffrey R. Hutchison

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

#ifndef OB_CHIRAL_H
#define OB_CHIRAL_H

#include <vector>
#include <openbabel/matrix.h>

namespace OpenBabel
{

//! @deprecated
OBAPI void GraphPotentials(OBMol &mol, std::vector<double> &pot);
//! @deprecated
OBAPI void construct_g_matrix(OBMol &mol, std::vector<std::vector<double> > &m);
//! @deprecated
OBAPI void construct_c_matrix(OBMol &mol, std::vector<std::vector<double > > &m);

//! @deprecated Use new @ref stereo classes.
OBAPI double CalcSignedVolume(OBMol &mol,OBAtom*,bool ReZeroZ=true);
//! @deprecated Use new @ref stereo classes.
OBAPI double signed_volume(const vector3 &a, const vector3 &b, const vector3 &c, const vector3 &d);
//! @deprecated Use new @ref stereo classes.
OBAPI void GetChirality(OBMol &mol, std::vector<int> &chirality);
//! @deprecated Use new @ref stereo classes.
OBAPI int GetParity4Ref(std::vector<unsigned int> pref);
//! @deprecated Use new @ref stereo classes.
OBAPI bool CorrectChirality(OBMol &mol, OBAtom *atm, atomreftype i=input, atomreftype o=output);

}

#endif // OB_CHIRAL_H

//! \file chiral.h
//! \brief Detect chiral atoms and molecules.

/* -*-C-*-

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
| FILE:         cansmilesformat.h
| AUTHOR:       Craig A. James
| DESCRIPTION:
|       Declaration for cansmilesformat.cpp
+======================================================================
*/

namespace OpenBabel {

void OpenBabel::CreateCansmiString(OBMol &mol, char *buffer, OBBitVec &frag_atoms, bool iso);

}

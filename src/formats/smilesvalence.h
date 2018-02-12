/**********************************************************************
smilesvalence.h - Implement SMILES valence model.

Copyright (C) 2016,2017 by Noel O'Boyle

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

static bool IsOutsideOrganicSubset(unsigned int elem)
{
  switch (elem) {
  case  0: // *
  case  5: // B
  case  6: // C
  case  7: // N
  case 15: // P
  case  8: // O
  case 16: // S
  case  9: // F
  case 17: // Cl
  case 35: // Br
  case 53: // I
    return false;
  default:
    return true;
  }
}

/* Return the implicit Smiles valence for element "elem" with neutral charge
 * and "val" explicit nbrs. When writing SMILES, a return value of 0 indicates a hypervalent structure
 */
static unsigned int SmilesValence(unsigned int elem, unsigned int val, bool reading=true)
{
  switch (elem) {
  case  5: // B
    if (val <= 3) return 3;
    break;
  case  6: // C
    if (val <= 4) return 4;
    break;
  case  7: // N
  case 15: // P
    switch (val) {
    case 0:
    case 1:
    case 2:
    case 3:
      return 3;
    case 4:
    case 5:
      return 5;
    }
    break;
  case  8: // O
    if (val <= 2) return 2;
    break;
  case 16:  // S
    switch (val) {
    case 0:
    case 1:
    case 2:
      return 2;
    case 3:
    case 4:
      return 4;
    case 5:
    case 6:
      return 6;
    }
    break;
  case  9:  // F
  case 17:  // Cl
  case 35:  // Br
  case 53:  // I
    if (val <= 1) return 1;
    break;
  }
  return reading ? val : 0;
}
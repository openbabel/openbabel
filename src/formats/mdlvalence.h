/**********************************************************************
mdlvalence.h - Implement MDL valence model.

Copyright (C) 2012 by NextMove Software

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

/* Return the implicit MDL valence for element "elem" with charge "q".  */
static unsigned int MDLValence(unsigned int elem, int q, unsigned int val)
{
  switch (elem) {
  case  1:  // H
  case  3:  // Li
  case 11:  // Na
  case 19:  // K
  case 37:  // Rb
  case 55:  // Cs
  case 87:  // Fr
    if (q == 0 && val <= 1)
      return 1;
    break;

  case  4:  // Be
  case 12:  // Mg
  case 20:  // Ca
  case 38:  // Sr
  case 56:  // Ba
  case 88:  // Ra
    switch (q) {
    case 0:  if (val <= 2) return 2;  break;
    case 1:  if (val <= 1) return 1;  break;
    }
    break;

  case  5:  // B
    switch (q) {
    case -4:  if (val <= 1) return 1;  break;
    case -3:  if (val <= 2) return 2;  break;
    case -2:  if (val <= 3) return 3;
              if (val <= 5) return 5;  break;
    case -1:  if (val <= 4) return 4;  break;
    case  0:  if (val <= 3) return 3;  break;
    case  1:  if (val <= 2) return 2;  break;
    case  2:  if (val <= 1) return 1;  break;
    }
    break;

  case  6:  // C
    switch (q) {
    case -3:  if (val <= 1) return 1;  break;
    case -2:  if (val <= 2) return 2;  break;
    case -1:  if (val <= 3) return 3;
              if (val <= 5) return 5;  break;
    case  0:  if (val <= 4) return 4;  break;
    case  1:  if (val <= 3) return 3;  break;
    case  2:  if (val <= 2) return 2;  break;
    case  3:  if (val <= 1) return 1;  break;
    }
    break;

  case  7:  // N
    switch (q) {
    case -2:  if (val <= 1) return 1;  break;
    case -1:  if (val <= 2) return 2;  break;
    case  0:  if (val <= 3) return 3;
              if (val <= 5) return 5;  break;
    case  1:  if (val <= 4) return 4;  break;
    case  2:  if (val <= 3) return 3;  break;
    case  3:  if (val <= 2) return 2;  break;
    case  4:  if (val <= 1) return 1;  break;
    }
    break;

  case  8:  // O
    switch (q) {
    case -1:  if (val <= 1) return 1;  break;
    case  0:  if (val <= 2) return 2;  break;
    case  1:  if (val <= 3) return 3;
              if (val <= 5) return 5;  break;
    case  2:  if (val <= 4) return 4;  break;
    case  3:  if (val <= 3) return 3;  break;
    case  4:  if (val <= 2) return 2;  break;
    case  5:  if (val <= 1) return 1;  break;
    }
    break;

  case  9:  // F
    switch (q) {
    case  0:  if (val <= 1) return 1;  break;
    case  1:  if (val <= 2) return 2;  break;
    case  2:  if (val <= 3) return 3;
              if (val <= 5) return 5;  break;
    case  3:  if (val <= 4) return 4;  break;
    case  4:  if (val <= 3) return 3;  break;
    case  5:  if (val <= 2) return 2;  break;
    case  6:  if (val <= 1) return 1;  break;
    }
    break;

  case 13:  // Al
    switch (q) {
    case -4:  if (val <= 1) return 1;
              if (val <= 3) return 3;
              if (val <= 5) return 5;
              if (val <= 7) return 7;  break;
    case -3:  if (val <= 2) return 2;
              if (val <= 4) return 4;
              if (val <= 6) return 6;  break;
    case -2:  if (val <= 3) return 3;
              if (val <= 5) return 5;  break;
    case -1:  if (val <= 4) return 4;  break;
    case  0:  if (val <= 3) return 3;  break;
    case  1:  if (val <= 2) return 2;  break;
    case  2:  if (val <= 1) return 1;  break;
    }
    break;

  case 14:  // Si
    switch (q) {
    case -3:  if (val <= 1) return 1;
              if (val <= 3) return 3;
              if (val <= 5) return 5;
              if (val <= 7) return 7;  break;
    case -2:  if (val <= 2) return 2;
              if (val <= 4) return 4;
              if (val <= 6) return 6;  break;
    case -1:  if (val <= 3) return 3;
              if (val <= 5) return 5;  break;
    case  0:  if (val <= 4) return 4;  break;
    case  1:  if (val <= 3) return 3;  break;
    case  2:  if (val <= 2) return 2;  break;
    case  3:  if (val <= 1) return 1;  break;
    }
    break;

  case 15:  // P
    switch (q) {
    case -2:  if (val <= 1) return 1;
              if (val <= 3) return 3;
              if (val <= 5) return 5;
              if (val <= 7) return 7;  break;
    case -1:  if (val <= 2) return 2;
              if (val <= 4) return 4;
              if (val <= 6) return 6;  break;
    case  0:  if (val <= 3) return 3;
              if (val <= 5) return 5;  break;
    case  1:  if (val <= 4) return 4;  break;
    case  2:  if (val <= 3) return 3;  break;
    case  3:  if (val <= 2) return 2;  break;
    case  4:  if (val <= 1) return 1;  break;
    }
    break;

  case 16:  // S
    switch (q) {
    case -1:  if (val <= 1) return 1;
              if (val <= 3) return 3;
              if (val <= 5) return 5;
              if (val <= 7) return 7;  break;
    case  0:  if (val <= 2) return 2;
              if (val <= 4) return 4;
              if (val <= 6) return 6;  break;
    case  1:  if (val <= 3) return 3;
              if (val <= 5) return 5;  break;
    case  2:  if (val <= 4) return 4;  break;
    case  3:  if (val <= 3) return 3;  break;
    case  4:  if (val <= 2) return 2;  break;
    case  5:  if (val <= 1) return 1;  break;
    }
    break;

  case 17:  // Cl
    switch (q) {
    case  0:  if (val <= 1) return 1;
              if (val <= 3) return 3;
              if (val <= 5) return 5;
              if (val <= 7) return 7;  break;
    case  1:  if (val <= 2) return 2;
              if (val <= 4) return 4;
              if (val <= 6) return 6;  break;
    case  2:  if (val <= 3) return 3;
              if (val <= 5) return 5;  break;
    case  3:  if (val <= 4) return 4;  break;
    case  4:  if (val <= 3) return 3;  break;
    case  5:  if (val <= 2) return 2;  break;
    case  6:  if (val <= 1) return 1;  break;
    }
    break;

  case 31:  // Ga
    switch (q) {
    case -4:  if (val <= 1) return 1;
              if (val <= 3) return 3;
              if (val <= 5) return 5;
              if (val <= 7) return 7;  break;
    case -3:  if (val <= 2) return 2;
              if (val <= 4) return 4;
              if (val <= 6) return 6;  break;
    case -2:  if (val <= 3) return 3;
              if (val <= 5) return 5;  break;
    case -1:  if (val <= 4) return 4;  break;
    case  0:  if (val <= 3) return 3;  break;
    case  2:  if (val <= 1) return 1;  break;
    }
    break;

  case 32:  // Ge
    switch (q) {
    case -3:  if (val <= 1) return 1;
              if (val <= 3) return 3;
              if (val <= 5) return 5;
              if (val <= 7) return 7;  break;
    case -2:  if (val <= 2) return 2;
              if (val <= 4) return 4;
              if (val <= 6) return 6;  break;
    case -1:  if (val <= 3) return 3;
              if (val <= 5) return 5;  break;
    case  0:  if (val <= 4) return 4;  break;
    case  1:  if (val <= 3) return 3;  break;
    case  3:  if (val <= 1) return 1;  break;
    }
    break;

  case 33:  // As
    switch (q) {
    case -2:  if (val <= 1) return 1;
              if (val <= 3) return 3;
              if (val <= 5) return 5;
              if (val <= 7) return 7;  break;
    case -1:  if (val <= 2) return 2;
              if (val <= 4) return 4;
              if (val <= 6) return 6;  break;
    case  0:  if (val <= 3) return 3;
              if (val <= 5) return 5;  break;
    case  1:  if (val <= 4) return 4;  break;
    case  2:  if (val <= 3) return 3;  break;
    case  4:  if (val <= 1) return 1;  break;
    }
    break;

  case 34:  // Se
    switch (q) {
    case -1:  if (val <= 1) return 1;
              if (val <= 3) return 3;
              if (val <= 5) return 5;
              if (val <= 7) return 7;  break;
    case  0:  if (val <= 2) return 2;
              if (val <= 4) return 4;
              if (val <= 6) return 6;  break;
    case  1:  if (val <= 3) return 3;
              if (val <= 5) return 5;  break;
    case  2:  if (val <= 4) return 4;  break;
    case  3:  if (val <= 3) return 3;  break;
    case  5:  if (val <= 1) return 1;  break;
    }
    break;

  case 35:  // Br
    switch (q) {
    case  0:  if (val <= 1) return 1;
              if (val <= 3) return 3;
              if (val <= 5) return 5;
              if (val <= 7) return 7;  break;
    case  1:  if (val <= 2) return 2;
              if (val <= 4) return 4;
              if (val <= 6) return 6;  break;
    case  2:  if (val <= 3) return 3;
              if (val <= 5) return 5;  break;
    case  3:  if (val <= 4) return 4;  break;
    case  4:  if (val <= 3) return 3;  break;
    case  6:  if (val <= 1) return 1;  break;
    }
    break;

  case 49:  // In
    switch (q) {
    case -4:  if (val <= 1) return 1;
              if (val <= 3) return 3;
              if (val <= 5) return 5;
              if (val <= 7) return 7;  break;
    case -3:  if (val <= 2) return 2;
              if (val <= 4) return 4;
              if (val <= 6) return 6;  break;
    case -2:  if (val <= 3) return 3;
              if (val <= 5) return 5;  break;
    case -1:  if (val <= 2) return 2;
              if (val <= 4) return 4;  break;
    case  0:  if (val <= 3) return 3;  break;
    case  2:  if (val <= 1) return 1;  break;
    }
    break;

  case 50:  // Sn
  case 82:  // Pb
    switch (q) {
    case -3:  if (val <= 1) return 1;
              if (val <= 3) return 3;
              if (val <= 5) return 5;
              if (val <= 7) return 7;  break;
    case -2:  if (val <= 2) return 2;
              if (val <= 4) return 4;
              if (val <= 6) return 6;  break;
    case -1:  if (val <= 3) return 3;
              if (val <= 5) return 5;  break;
    case  0:  if (val <= 2) return 2;
              if (val <= 4) return 4;  break;
    case  1:  if (val <= 3) return 3;  break;
    case  3:  if (val <= 1) return 1;  break;
    }
    break;

  case 51:  // Sb
  case 83:  // Bi
    switch (q) {
    case -2:  if (val <= 1) return 1;
              if (val <= 3) return 3;
              if (val <= 5) return 5;
              if (val <= 7) return 7;  break;
    case -1:  if (val <= 2) return 2;
              if (val <= 4) return 4;
              if (val <= 6) return 6;  break;
    case  0:  if (val <= 3) return 3;
              if (val <= 5) return 5;  break;
    case  1:  if (val <= 2) return 2;
              if (val <= 4) return 4;  break;
    case  2:  if (val <= 3) return 3;  break;
    case  4:  if (val <= 1) return 1;  break;
    }
    break;

  case 52:  // Te
  case 84:  // Po
    switch (q) {
    case -1:  if (val <= 1) return 1;
              if (val <= 3) return 3;
              if (val <= 5) return 5;
              if (val <= 7) return 7;  break;
    case  0:  if (val <= 2) return 2;
              if (val <= 4) return 4;
              if (val <= 6) return 6;  break;
    case  1:  if (val <= 3) return 3;
              if (val <= 5) return 5;  break;
    case  2:  if (val <= 2) return 2;
              if (val <= 4) return 4;  break;
    case  3:  if (val <= 3) return 3;  break;
    case  5:  if (val <= 1) return 1;  break;
    }
    break;

  case 53:  // I
  case 85:  // At
    switch (q) {
    case  0:  if (val <= 1) return 1;
              if (val <= 3) return 3;
              if (val <= 5) return 5;
              if (val <= 7) return 7;  break;
    case  1:  if (val <= 2) return 2;
              if (val <= 4) return 4;
              if (val <= 6) return 6;  break;
    case  2:  if (val <= 3) return 3;
              if (val <= 5) return 5;  break;
    case  3:  if (val <= 2) return 2;
              if (val <= 4) return 4;  break;
    case  4:  if (val <= 3) return 3;  break;
    case  6:  if (val <= 1) return 1;  break;
    }
    break;

  case 81:  // Tl
    switch (q) {
    case -4:  if (val <= 1) return 1;
              if (val <= 3) return 3;
              if (val <= 5) return 5;
              if (val <= 7) return 7;  break;
    case -3:  if (val <= 2) return 2;
              if (val <= 4) return 4;
              if (val <= 6) return 6;  break;
    case -2:  if (val <= 3) return 3;
              if (val <= 5) return 5;  break;
    case -1:  if (val <= 2) return 2;
              if (val <= 4) return 4;  break;
    case  0:  if (val <= 1) return 1;
              if (val <= 3) return 3;  break;
    }
    break;

  }
  return val;
}

/* Return the implicit valence for element "elem" with charge "q" when using HYD extension */
static unsigned int HYDValence(unsigned int elem, int q, unsigned int val)
{
  int impval = 0;
  if (elem == 6) {  // C
    impval = 4 - abs(q);
  } else if (elem == 7 || elem == 15) {  // N or P
    impval = 3 + q;
  } else if (elem == 8 || elem == 16) {  // O or S
    impval = 2 + q;
  }
  if (impval < 0) {
    impval = 0;
  }
  if (val > impval) {
    impval = val;
  }
  return impval;  
}


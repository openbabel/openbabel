/**********************************************************************
elements.cpp - Handle element conversion

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

#include <openbabel/elements.h>
#include <vector>

#define NUMELEMENTS 118

// Use the C++ preprocessor to turn a CSV into separate arrays per column
// - FYI, the "#B" means put the B in quotes. In you do it directly (i.e. "B")
//   the macro argument is not expanded and every element is set to the letter B.
// - NUMELEMENTS + 1 is used as there is a dummy atom at the start. Also for
//   direct array lookup using the atomic number (starts at 1).
#define ELEMTABLE(A, B, C, D, E, F, G, H, I, J, K, L, M, N, O) #B
const char* element_symbols[NUMELEMENTS + 1] = {
#include "elementtable.h"
};
#undef ELEMTABLE
#define ELEMTABLE(A, B, C, D, E, F, G, H, I, J, K, L, M, N, O) C
double eneg_allred[NUMELEMENTS + 1] = {
#include "elementtable.h"
};
#undef ELEMTABLE
#define ELEMTABLE(A, B, C, D, E, F, G, H, I, J, K, L, M, N, O) D
double cov_rad[NUMELEMENTS + 1] = {
#include "elementtable.h"
};
#undef ELEMTABLE
#define ELEMTABLE(A, B, C, D, E, F, G, H, I, J, K, L, M, N, O) F
double vdw_rad[NUMELEMENTS + 1] = {
#include "elementtable.h"
};
#undef ELEMTABLE
#define ELEMTABLE(A, B, C, D, E, F, G, H, I, J, K, L, M, N, O) G
unsigned int maxbonds[NUMELEMENTS + 1] = {
#include "elementtable.h"
};
#undef ELEMTABLE
#define ELEMTABLE(A, B, C, D, E, F, G, H, I, J, K, L, M, N, O) H
double atomic_mass[NUMELEMENTS + 1] = {
#include "elementtable.h"
};
#undef ELEMTABLE
#define ELEMTABLE(A, B, C, D, E, F, G, H, I, J, K, L, M, N, O) I
double eneg_pauling[NUMELEMENTS + 1] = {
#include "elementtable.h"
};
#undef ELEMTABLE
#define ELEMTABLE(A, B, C, D, E, F, G, H, I, J, K, L, M, N, O) J
double ionization[NUMELEMENTS + 1] = {
#include "elementtable.h"
};
#undef ELEMTABLE
#define ELEMTABLE(A, B, C, D, E, F, G, H, I, J, K, L, M, N, O) K
double electron_affinity[NUMELEMENTS + 1] = {
#include "elementtable.h"
};
#undef ELEMTABLE
#define ELEMTABLE(A, B, C, D, E, F, G, H, I, J, K, L, M, N, O) {L, M, N}
double rgb[NUMELEMENTS + 1][3] = {
#include "elementtable.h"
};
#undef ELEMTABLE
#define ELEMTABLE(A, B, C, D, E, F, G, H, I, J, K, L, M, N, O) #O
const char* element_name[NUMELEMENTS + 1] = {
#include "elementtable.h"
};


namespace OpenBabel
{
  namespace OBElements {

    const char* GetSymbol(unsigned int atomic_number)
    {
      if (atomic_number > NUMELEMENTS)
        return "";
      return element_symbols[atomic_number];
    }

    const char* GetName(unsigned int atomic_number)
    {
      if (atomic_number > NUMELEMENTS)
        return "";
      return element_name[atomic_number];
    }

    double GetMass(unsigned int atomic_number)
    {
      if (atomic_number > NUMELEMENTS)
        return 0.0;
      return atomic_mass[atomic_number];
    }

    double GetAllredRochowElectroNeg(unsigned int atomic_number)
    {
      if (atomic_number > NUMELEMENTS)
        return 0.0;
      return eneg_allred[atomic_number];
    }

    double GetCovalentRad(unsigned int atomic_number)
    {
      if (atomic_number > NUMELEMENTS)
        return 0.0;
      return cov_rad[atomic_number];
    }

    double GetVdwRad(unsigned int atomic_number)
    {
      if (atomic_number > NUMELEMENTS)
        return 0.0;
      return vdw_rad[atomic_number];
    }

    double GetElectronAffinity(unsigned int atomic_number)
    {
      if (atomic_number > NUMELEMENTS)
        return 0.0;
      return electron_affinity[atomic_number];
    }

    double GetIonization(unsigned int atomic_number)
    {
      if (atomic_number > NUMELEMENTS)
        return 0.0;
      return ionization[atomic_number];
    }

    unsigned int GetMaxBonds(unsigned int atomic_number)
    {
      if (atomic_number > NUMELEMENTS)
        return 0;
      return maxbonds[atomic_number];
    }

    double GetElectroNeg(unsigned int atomic_number)
    {
      if (atomic_number > NUMELEMENTS)
        return 0.0;
      return eneg_pauling[atomic_number];
    }

    void GetRGB(unsigned int atomic_number, double *r, double *g, double *b)
    {
      if (atomic_number > NUMELEMENTS) {
        *r = 0;
        *g = 0;
        *b = 0;
      }
      else {
        double *ans = rgb[atomic_number];
        *r = *ans;
        *g = *(ans + 1);
        *b = *(ans + 2);
      }
    }

    unsigned int GetAtomicNum(const char* ptr)
    {
      switch (ptr[0]) {
      case 'A':
        switch (ptr[1]) {
        case 'c':
          if (ptr[2] == '\0') { // Ac
            return 89; // Actinium
          }
          break;
        case 'g':
          if (ptr[2] == '\0') { // Ag
            return 47; // Silver
          }
          break;
        case 'l':
          if (ptr[2] == '\0') { // Al
            return 13; // Aluminium
          }
          break;
        case 'm':
          if (ptr[2] == '\0') { // Am
            return 95; // Americium
          }
          break;
        case 'r':
          if (ptr[2] == '\0') { // Ar
            return 18; // Argon
          }
          break;
        case 's':
          if (ptr[2] == '\0') { // As
            return 33; // Arsenic
          }
          break;
        case 't':
          if (ptr[2] == '\0') { // At
            return 85; // Astatine
          }
          break;
        case 'u':
          if (ptr[2] == '\0') { // Au
            return 79; // Gold
          }
          break;
        }
        break;
      case 'B':
        switch (ptr[1]) {
        case '\0': // B
          return 5; // Boron
        case 'a':
          if (ptr[2] == '\0') { // Ba
            return 56; // Barium
          }
          break;
        case 'e':
          if (ptr[2] == '\0') { // Be
            return 4; // Beryllium
          }
          break;
        case 'h':
          if (ptr[2] == '\0') { // Bh
            return 107; // Bohrium
          }
          break;
        case 'i':
          if (ptr[2] == '\0') { // Bi
            return 83; // Bismuth
          }
          break;
        case 'k':
          if (ptr[2] == '\0') { // Bk
            return 97; // Berkelium
          }
          break;
        case 'r':
          if (ptr[2] == '\0') { // Br
            return 35; // Bromine
          }
          break;
        }
        break;
      case 'C':
        switch (ptr[1]) {
        case '\0': // C
          return 6; // Carbon
        case 'a':
          if (ptr[2] == '\0') { // Ca
            return 20; // Calcium
          }
          break;
        case 'd':
          if (ptr[2] == '\0') { // Cd
            return 48; // Cadmium
          }
          break;
        case 'e':
          if (ptr[2] == '\0') { // Ce
            return 58; // Cerium
          }
          break;
        case 'f':
          if (ptr[2] == '\0') { // Cf
            return 98; // Californium
          }
          break;
        case 'l':
          if (ptr[2] == '\0') { // Cl
            return 17; // Chlorine
          }
          break;
        case 'm':
          if (ptr[2] == '\0') { // Cm
            return 96; // Curium
          }
          break;
        case 'n':
          if (ptr[2] == '\0') { // Cn
            return 112; // Copernicium
          }
          break;
        case 'o':
          if (ptr[2] == '\0') { // Co
            return 27; // Cobalt
          }
          break;
        case 'r':
          if (ptr[2] == '\0') { // Cr
            return 24; // Chromium
          }
          break;
        case 's':
          if (ptr[2] == '\0') { // Cs
            return 55; // Caesium
          }
          break;
        case 'u':
          if (ptr[2] == '\0') { // Cu
            return 29; // Copper
          }
          break;
        }
        break;
      case 'D':
        switch (ptr[1]) {
        case 'b':
          if (ptr[2] == '\0') { // Db
            return 105; // Dubnium
          }
          break;
        case 's':
          if (ptr[2] == '\0') { // Ds
            return 110; // Darmstadtium
          }
          break;
        case 'y':
          if (ptr[2] == '\0') { // Dy
            return 66; // Dysprosium
          }
          break;
        }
        break;
      case 'E':
        switch (ptr[1]) {
        case 'r':
          if (ptr[2] == '\0') { // Er
            return 68; // Erbium
          }
          break;
        case 's':
          if (ptr[2] == '\0') { // Es
            return 99; // Einsteinium
          }
          break;
        case 'u':
          if (ptr[2] == '\0') { // Eu
            return 63; // Europium
          }
          break;
        }
        break;
      case 'F':
        switch (ptr[1]) {
        case '\0': // F
          return 9; // Fluorine
        case 'e':
          if (ptr[2] == '\0') { // Fe
            return 26; // Iron
          }
          break;
        case 'l':
          if (ptr[2] == '\0') { // Fl
            return 114; // Flerovium
          }
          break;
        case 'm':
          if (ptr[2] == '\0') { // Fm
            return 100; // Fermium
          }
          break;
        case 'r':
          if (ptr[2] == '\0') { // Fr
            return 87; // Francium
          }
          break;
        }
        break;
      case 'G':
        switch (ptr[1]) {
        case 'a':
          if (ptr[2] == '\0') { // Ga
            return 31; // Gallium
          }
          break;
        case 'd':
          if (ptr[2] == '\0') { // Gd
            return 64; // Gadolinium
          }
          break;
        case 'e':
          if (ptr[2] == '\0') { // Ge
            return 32; // Germanium
          }
          break;
        }
        break;
      case 'H':
        switch (ptr[1]) {
        case '\0': // H
          return 1; // Hydrogen
        case 'e':
          if (ptr[2] == '\0') { // He
            return 2; // Helium
          }
          break;
        case 'f':
          if (ptr[2] == '\0') { // Hf
            return 72; // Hafnium
          }
          break;
        case 'g':
          if (ptr[2] == '\0') { // Hg
            return 80; // Mercury
          }
          break;
        case 'o':
          if (ptr[2] == '\0') { // Ho
            return 67; // Holmium
          }
          break;
        case 's':
          if (ptr[2] == '\0') { // Hs
            return 108; // Hassium
          }
          break;
        }
        break;
      case 'I':
        switch (ptr[1]) {
        case '\0': // I
          return 53; // Iodine
        case 'n':
          if (ptr[2] == '\0') { // In
            return 49; // Indium
          }
          break;
        case 'r':
          if (ptr[2] == '\0') { // Ir
            return 77; // Iridium
          }
          break;
        }
        break;
      case 'K':
        switch (ptr[1]) {
        case '\0': // K
          return 19; // Potassium
        case 'r':
          if (ptr[2] == '\0') { // Kr
            return 36; // Krypton
          }
          break;
        }
        break;
      case 'L':
        switch (ptr[1]) {
        case 'a':
          if (ptr[2] == '\0') { // La
            return 57; // Lanthanum
          }
          break;
        case 'i':
          if (ptr[2] == '\0') { // Li
            return 3; // Lithium
          }
          break;
        case 'r':
          if (ptr[2] == '\0') { // Lr
            return 103; // Lawrencium
          }
          break;
        case 'u':
          if (ptr[2] == '\0') { // Lu
            return 71; // Lutetium
          }
          break;
        case 'v':
          if (ptr[2] == '\0') { // Lv
            return 116; // Livermorium
          }
          break;
        }
        break;
      case 'M':
        switch (ptr[1]) {
        case 'c':
          if (ptr[2] == '\0') { // Mc
            return 115; // Moscovium
          }
          break;
        case 'd':
          if (ptr[2] == '\0') { // Md
            return 101; // Mendelevium
          }
          break;
        case 'g':
          if (ptr[2] == '\0') { // Mg
            return 12; // Magnesium
          }
          break;
        case 'n':
          if (ptr[2] == '\0') { // Mn
            return 25; // Manganese
          }
          break;
        case 'o':
          if (ptr[2] == '\0') { // Mo
            return 42; // Molybdenum
          }
          break;
        case 't':
          if (ptr[2] == '\0') { // Mt
            return 109; // Meitnerium
          }
          break;
        }
        break;
      case 'N':
        switch (ptr[1]) {
        case '\0': // N
          return 7; // Nitrogen
        case 'a':
          if (ptr[2] == '\0') { // Na
            return 11; // Sodium
          }
          break;
        case 'b':
          if (ptr[2] == '\0') { // Nb
            return 41; // Niobium
          }
          break;
        case 'd':
          if (ptr[2] == '\0') { // Nd
            return 60; // Neodymium
          }
          break;
        case 'e':
          if (ptr[2] == '\0') { // Ne
            return 10; // Neon
          }
          break;
        case 'h':
          if (ptr[2] == '\0') { // Nh
            return 113; // Nihonium
          }
          break;
        case 'i':
          if (ptr[2] == '\0') { // Ni
            return 28; // Nickel
          }
          break;
        case 'o':
          if (ptr[2] == '\0') { // No
            return 102; // Nobelium
          }
          break;
        case 'p':
          if (ptr[2] == '\0') { // Np
            return 93; // Neptunium
          }
          break;
        }
        break;
      case 'O':
        switch (ptr[1]) {
        case '\0': // O
          return 8; // Oxygen
        case 'g':
          if (ptr[2] == '\0') { // Og
            return 118; // Oganesson
          }
          break;
        case 's':
          if (ptr[2] == '\0') { // Os
            return 76; // Osmium
          }
          break;
        }
        break;
      case 'P':
        switch (ptr[1]) {
        case '\0': // P
          return 15; // Phosphorus
        case 'a':
          if (ptr[2] == '\0') { // Pa
            return 91; // Protactinium
          }
          break;
        case 'b':
          if (ptr[2] == '\0') { // Pb
            return 82; // Lead
          }
          break;
        case 'd':
          if (ptr[2] == '\0') { // Pd
            return 46; // Palladium
          }
          break;
        case 'm':
          if (ptr[2] == '\0') { // Pm
            return 61; // Promethium
          }
          break;
        case 'o':
          if (ptr[2] == '\0') { // Po
            return 84; // Polonium
          }
          break;
        case 'r':
          if (ptr[2] == '\0') { // Pr
            return 59; // Praseodymium
          }
          break;
        case 't':
          if (ptr[2] == '\0') { // Pt
            return 78; // Platinum
          }
          break;
        case 'u':
          if (ptr[2] == '\0') { // Pu
            return 94; // Plutonium
          }
          break;
        }
        break;
      case 'R':
        switch (ptr[1]) {
        case 'a':
          if (ptr[2] == '\0') { // Ra
            return 88; // Radium
          }
          break;
        case 'b':
          if (ptr[2] == '\0') { // Rb
            return 37; // Rubidium
          }
          break;
        case 'e':
          if (ptr[2] == '\0') { // Re
            return 75; // Rhenium
          }
          break;
        case 'f':
          if (ptr[2] == '\0') { // Rf
            return 104; // Rutherfordium
          }
          break;
        case 'g':
          if (ptr[2] == '\0') { // Rg
            return 111; // Roentgenium
          }
          break;
        case 'h':
          if (ptr[2] == '\0') { // Rh
            return 45; // Rhodium
          }
          break;
        case 'n':
          if (ptr[2] == '\0') { // Rn
            return 86; // Radon
          }
          break;
        case 'u':
          if (ptr[2] == '\0') { // Ru
            return 44; // Ruthenium
          }
          break;
        }
        break;
      case 'S':
        switch (ptr[1]) {
        case '\0': // S
          return 16; // Sulfur
        case 'b':
          if (ptr[2] == '\0') { // Sb
            return 51; // Antimony
          }
          break;
        case 'c':
          if (ptr[2] == '\0') { // Sc
            return 21; // Scandium
          }
          break;
        case 'e':
          if (ptr[2] == '\0') { // Se
            return 34; // Selenium
          }
          break;
        case 'g':
          if (ptr[2] == '\0') { // Sg
            return 106; // Seaborgium
          }
          break;
        case 'i':
          if (ptr[2] == '\0') { // Si
            return 14; // Silicon
          }
          break;
        case 'm':
          if (ptr[2] == '\0') { // Sm
            return 62; // Samarium
          }
          break;
        case 'n':
          if (ptr[2] == '\0') { // Sn
            return 50; // Tin
          }
          break;
        case 'r':
          if (ptr[2] == '\0') { // Sr
            return 38; // Strontium
          }
          break;
        }
        break;
      case 'T':
        switch (ptr[1]) {
        case 'a':
          if (ptr[2] == '\0') { // Ta
            return 73; // Tantalum
          }
          break;
        case 'b':
          if (ptr[2] == '\0') { // Tb
            return 65; // Terbium
          }
          break;
        case 'c':
          if (ptr[2] == '\0') { // Tc
            return 43; // Technetium
          }
          break;
        case 'e':
          if (ptr[2] == '\0') { // Te
            return 52; // Tellurium
          }
          break;
        case 'h':
          if (ptr[2] == '\0') { // Th
            return 90; // Thorium
          }
          break;
        case 'i':
          if (ptr[2] == '\0') { // Ti
            return 22; // Titanium
          }
          break;
        case 'l':
          if (ptr[2] == '\0') { // Tl
            return 81; // Thallium
          }
          break;
        case 'm':
          if (ptr[2] == '\0') { // Tm
            return 69; // Thulium
          }
          break;
        case 's':
          if (ptr[2] == '\0') { // Ts
            return 117; // Tennessine
          }
          break;
        }
        break;
      case 'U':
        if (ptr[1] == '\0') { // U
          return 92; // Uranium
        }
        break;
      case 'V':
        if (ptr[1] == '\0') { // V
          return 23; // Vanadium
        }
        break;
      case 'W':
        if (ptr[1] == '\0') { // W
          return 74; // Tungsten
        }
        break;
      case 'X':
        if (ptr[1] == 'e' && ptr[2] == '\0') { // Xe
          return 54; // Xenon
        }
        break;
      case 'Y':
        switch (ptr[1]) {
        case '\0': // Y
          return 39; // Yttrium
        case 'b':
          if (ptr[2] == '\0') { // Yb
            return 70; // Ytterbium
          }
          break;
        }
        break;
      case 'Z':
        switch (ptr[1]) {
        case 'n':
          if (ptr[2] == '\0') { // Zn
            return 30; // Zinc
          }
          break;
        case 'r':
          if (ptr[2] == '\0') { // Zr
            return 40; // Zirconium
          }
          break;
        }
        break;
      }
      return 0;
    }
  } // namespace OBElements
} // namespace OpenBabel
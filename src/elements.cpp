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

//#  CAUTION: Masses have been zero - padded to 9 decimal places               #
//#   (we don't currently have a math package that can deal with sig. fig.)    #
//#   so trailing zeros may not be significant and / or useful                 #
//#                                                                            #
//#  Values are compiled from :                                                #
//#      "The Ame2003 atomic mass evaluation (II)"                             #
//#        by G.Audi, A.H.Wapstra and C.Thibault                               #
//#          Nuclear Physics A729 p. 337 - 676, December 22, 2003.             #
//#  as made available in the mass.mas03round file                             #
//#  (these are values for publication)                                        #
//#                                                                            #
//#  Note that since element 0 often represents a dummy atom, the default      #
//#  mass is 0 amu, *not* the mass of a neutron, accessible as isotope 1       #
    double GetExactMass(unsigned int atomic_num, unsigned int isotope)
    {
      // Note that an isotope value of 0 (the default) should return the mass
      // of the most abundant isotope
      switch(atomic_num) {
      case 0:
        switch (isotope) {
        case 0: return 0.000000000;
        case 1: return 1.008664915;
        }
        break;
      case 1:
        switch (isotope) {
        case 0: case 1: return 1.007825032;
        case 2: return 2.014101778;
        case 3: return 3.016049278;
        case 4: return 4.027810000;
        case 5: return 5.035310000;
        case 6: return 6.044940000;
        case 7: return 7.052750000;
        }
        break;
      case 2:
        switch (isotope) {
        case 3: return 3.016029319;
        case 0: case 4: return 4.002603254;
        case 5: return 5.012220000;
        case 6: return 6.018889100;
        case 7: return 7.028021000;
        case 8: return 8.033922000;
        case 9: return 9.043950000;
        case 10: return 10.052400000;
        }
        break;
      case 3:
        switch (isotope) {
        case 3: return 3.030780000;
        case 4: return 4.027190000;
        case 5: return 5.012540000;
        case 6: return 6.015122795;
        case 0: case 7: return 7.016004550;
        case 8: return 8.022487360;
        case 9: return 9.026789500;
        case 10: return 10.035481000;
        case 11: return 11.043798000;
        case 12: return 12.053780000;
        }
        break;
      case 4:
        switch (isotope) {
        case 5: return 5.040790000;
        case 6: return 6.019726000;
        case 7: return 7.016929830;
        case 8: return 8.005305100;
        case 0: case 9: return 9.012182200;
        case 10: return 10.013533800;
        case 11: return 11.021658000;
        case 12: return 12.026921000;
        case 13: return 13.035690000;
        case 14: return 14.042890000;
        case 15: return 15.053460000;
        case 16: return 16.061920000;
        }
        break;
      case 5:
        switch (isotope) {
        case 6: return 6.046810000;
        case 7: return 7.029920000;
        case 8: return 8.024607200;
        case 9: return 9.013328800;
        case 10: return 10.012937000;
        case 0: case 11: return 11.009305400;
        case 12: return 12.014352100;
        case 13: return 13.017780200;
        case 14: return 14.025404000;
        case 15: return 15.031103000;
        case 16: return 16.039810000;
        case 17: return 17.046990000;
        case 18: return 18.056170000;
        case 19: return 19.063730000;
        }
        break;
      case 6:
        switch (isotope) {
        case 8: return 8.037675000;
        case 9: return 9.031036700;
        case 10: return 10.016853200;
        case 11: return 11.011433600;
        case 0: case 12: return 12.000000000;
        case 13: return 13.003354838;
        case 14: return 14.003241989;
        case 15: return 15.010599300;
        case 16: return 16.014701000;
        case 17: return 17.022586000;
        case 18: return 18.026760000;
        case 19: return 19.034810000;
        case 20: return 20.040320000;
        case 21: return 21.049340000;
        case 22: return 22.057200000;
        }
        break;
      case 7:
        switch (isotope) {
        case 10: return 10.041650000;
        case 11: return 11.026090000;
        case 12: return 12.018613200;
        case 13: return 13.005738610;
        case 0: case 14: return 14.003074005;
        case 15: return 15.000108898;
        case 16: return 16.006101700;
        case 17: return 17.008450000;
        case 18: return 18.014079000;
        case 19: return 19.017029000;
        case 20: return 20.023370000;
        case 21: return 21.027110000;
        case 22: return 22.034390000;
        case 23: return 23.041220000;
        case 24: return 24.051040000;
        case 25: return 25.060660000;
        }
        break;
      case 8:
        switch (isotope) {
        case 12: return 12.034405000;
        case 13: return 13.024812000;
        case 14: return 14.008596250;
        case 15: return 15.003065600;
        case 0: case 16: return 15.994914620;
        case 17: return 16.999131700;
        case 18: return 17.999161000;
        case 19: return 19.003580000;
        case 20: return 20.004076700;
        case 21: return 21.008656000;
        case 22: return 22.009970000;
        case 23: return 23.015690000;
        case 24: return 24.020470000;
        case 25: return 25.029460000;
        case 26: return 26.038340000;
        case 27: return 27.048260000;
        case 28: return 28.057810000;
        }
        break;
      case 9:
        switch (isotope) {
        case 14: return 14.035060000;
        case 15: return 15.018010000;
        case 16: return 16.011466000;
        case 17: return 17.002095240;
        case 18: return 18.000938000;
        case 0: case 19: return 18.998403220;
        case 20: return 19.999981320;
        case 21: return 20.999949000;
        case 22: return 22.002999000;
        case 23: return 23.003570000;
        case 24: return 24.008120000;
        case 25: return 25.012100000;
        case 26: return 26.019620000;
        case 27: return 27.026760000;
        case 28: return 28.035670000;
        case 29: return 29.043260000;
        case 30: return 30.052500000;
        case 31: return 31.060430000;
        }
        break;
      case 10:
        switch (isotope) {
        case 16: return 16.025761000;
        case 17: return 17.017672000;
        case 18: return 18.005708200;
        case 19: return 19.001880200;
        case 0: case 20: return 19.992440175;
        case 21: return 20.993846680;
        case 22: return 21.991385114;
        case 23: return 22.994466900;
        case 24: return 23.993610800;
        case 25: return 24.997737000;
        case 26: return 26.000461000;
        case 27: return 27.007590000;
        case 28: return 28.012070000;
        case 29: return 29.019390000;
        case 30: return 30.024800000;
        case 31: return 31.033110000;
        case 32: return 32.040020000;
        case 33: return 33.049380000;
        case 34: return 34.057030000;
        }
        break;
      case 11:
        switch (isotope) {
        case 18: return 18.025970000;
        case 19: return 19.013877000;
        case 20: return 20.007351000;
        case 21: return 20.997655200;
        case 22: return 21.994436400;
        case 0: case 23: return 22.989769281;
        case 24: return 23.990962780;
        case 25: return 24.989954000;
        case 26: return 25.992633000;
        case 27: return 26.994077000;
        case 28: return 27.998938000;
        case 29: return 29.002861000;
        case 30: return 30.008976000;
        case 31: return 31.013590000;
        case 32: return 32.020470000;
        case 33: return 33.026720000;
        case 34: return 34.035170000;
        case 35: return 35.042490000;
        case 36: return 36.051480000;
        case 37: return 37.059340000;
        }
        break;
      case 12:
        switch (isotope) {
        case 19: return 19.035470000;
        case 20: return 20.018863000;
        case 21: return 21.011713000;
        case 22: return 21.999573800;
        case 23: return 22.994123700;
        case 0: case 24: return 23.985041700;
        case 25: return 24.985836920;
        case 26: return 25.982592929;
        case 27: return 26.984340590;
        case 28: return 27.983876800;
        case 29: return 28.988600000;
        case 30: return 29.990434000;
        case 31: return 30.996546000;
        case 32: return 31.998975000;
        case 33: return 33.005254000;
        case 34: return 34.009460000;
        case 35: return 35.017340000;
        case 36: return 36.023000000;
        case 37: return 37.031400000;
        case 38: return 38.037570000;
        case 39: return 39.046770000;
        case 40: return 40.053930000;
        }
        break;
      case 13:
        switch (isotope) {
        case 21: return 21.028040000;
        case 22: return 22.019520000;
        case 23: return 23.007267000;
        case 24: return 23.999938900;
        case 25: return 24.990428100;
        case 26: return 25.986891690;
        case 0: case 27: return 26.981538630;
        case 28: return 27.981910310;
        case 29: return 28.980445000;
        case 30: return 29.982960000;
        case 31: return 30.983947000;
        case 32: return 31.988120000;
        case 33: return 32.990840000;
        case 34: return 33.996850000;
        case 35: return 34.999860000;
        case 36: return 36.006210000;
        case 37: return 37.010680000;
        case 38: return 38.017230000;
        case 39: return 39.022970000;
        case 40: return 40.031450000;
        case 41: return 41.038330000;
        case 42: return 42.046890000;
        }
        break;
      case 14:
        switch (isotope) {
        case 22: return 22.034530000;
        case 23: return 23.025520000;
        case 24: return 24.011546000;
        case 25: return 25.004106000;
        case 26: return 25.992330000;
        case 27: return 26.986704910;
        case 0: case 28: return 27.976926532;
        case 29: return 28.976494700;
        case 30: return 29.973770170;
        case 31: return 30.975363230;
        case 32: return 31.974148080;
        case 33: return 32.978000000;
        case 34: return 33.978576000;
        case 35: return 34.984580000;
        case 36: return 35.986600000;
        case 37: return 36.992940000;
        case 38: return 37.995630000;
        case 39: return 39.002070000;
        case 40: return 40.005870000;
        case 41: return 41.014560000;
        case 42: return 42.019790000;
        case 43: return 43.028660000;
        case 44: return 44.035260000;
        }
        break;
      case 15:
        switch (isotope) {
        case 24: return 24.034350000;
        case 25: return 25.020260000;
        case 26: return 26.011780000;
        case 27: return 26.999230000;
        case 28: return 27.992315000;
        case 29: return 28.981800600;
        case 30: return 29.978313800;
        case 0: case 31: return 30.973761630;
        case 32: return 31.973907270;
        case 33: return 32.971725500;
        case 34: return 33.973636000;
        case 35: return 34.973314100;
        case 36: return 35.978260000;
        case 37: return 36.979610000;
        case 38: return 37.984160000;
        case 39: return 38.986180000;
        case 40: return 39.991300000;
        case 41: return 40.994340000;
        case 42: return 42.001010000;
        case 43: return 43.006190000;
        case 44: return 44.012990000;
        case 45: return 45.019220000;
        case 46: return 46.027380000;
        }
        break;
      case 16:
        switch (isotope) {
        case 26: return 26.027880000;
        case 27: return 27.018830000;
        case 28: return 28.004370000;
        case 29: return 28.996610000;
        case 30: return 29.984903000;
        case 31: return 30.979554700;
        case 0: case 32: return 31.972071000;
        case 33: return 32.971458760;
        case 34: return 33.967866900;
        case 35: return 34.969032160;
        case 36: return 35.967080760;
        case 37: return 36.971125570;
        case 38: return 37.971163000;
        case 39: return 38.975130000;
        case 40: return 39.975450000;
        case 41: return 40.979580000;
        case 42: return 41.981020000;
        case 43: return 42.987150000;
        case 44: return 43.990210000;
        case 45: return 44.996510000;
        case 46: return 46.000750000;
        case 47: return 47.008590000;
        case 48: return 48.014170000;
        case 49: return 49.023620000;
        }
        break;
      case 17:
        switch (isotope) {
        case 28: return 28.028510000;
        case 29: return 29.014110000;
        case 30: return 30.004770000;
        case 31: return 30.992410000;
        case 32: return 31.985690000;
        case 33: return 32.977451900;
        case 34: return 33.973762820;
        case 0: case 35: return 34.968852680;
        case 36: return 35.968306980;
        case 37: return 36.965902590;
        case 38: return 37.968010430;
        case 39: return 38.968008200;
        case 40: return 39.970420000;
        case 41: return 40.970680000;
        case 42: return 41.973250000;
        case 43: return 42.974050000;
        case 44: return 43.978280000;
        case 45: return 44.980290000;
        case 46: return 45.984210000;
        case 47: return 46.988710000;
        case 48: return 47.994950000;
        case 49: return 49.000320000;
        case 50: return 50.007840000;
        case 51: return 51.014490000;
        }
        break;
      case 18:
        switch (isotope) {
        case 30: return 30.021560000;
        case 31: return 31.012120000;
        case 32: return 31.997638000;
        case 33: return 32.989925700;
        case 34: return 33.980271200;
        case 35: return 34.975257600;
        case 36: return 35.967545106;
        case 37: return 36.966776320;
        case 38: return 37.962732400;
        case 39: return 38.964313000;
        case 0: case 40: return 39.962383123;
        case 41: return 40.964500600;
        case 42: return 41.963046000;
        case 43: return 42.965636000;
        case 44: return 43.964924000;
        case 45: return 44.968040000;
        case 46: return 45.968090000;
        case 47: return 46.972190000;
        case 48: return 47.974540000;
        case 49: return 48.980520000;
        case 50: return 49.984430000;
        case 51: return 50.991630000;
        case 52: return 51.996780000;
        case 53: return 53.004940000;
        }
        break;
      case 19:
        switch (isotope) {
        case 32: return 32.021920000;
        case 33: return 33.007260000;
        case 34: return 33.998410000;
        case 35: return 34.988010000;
        case 36: return 35.981292000;
        case 37: return 36.973375890;
        case 38: return 37.969081200;
        case 0: case 39: return 38.963706680;
        case 40: return 39.963998480;
        case 41: return 40.961825760;
        case 42: return 41.962402810;
        case 43: return 42.960716000;
        case 44: return 43.961560000;
        case 45: return 44.960699000;
        case 46: return 45.961977000;
        case 47: return 46.961678000;
        case 48: return 47.965514000;
        case 49: return 48.967450000;
        case 50: return 49.972780000;
        case 51: return 50.976380000;
        case 52: return 51.982610000;
        case 53: return 52.987120000;
        case 54: return 53.994200000;
        case 55: return 54.999710000;
        }
        break;
      case 20:
        switch (isotope) {
        case 34: return 34.014120000;
        case 35: return 35.004940000;
        case 36: return 35.993090000;
        case 37: return 36.985870000;
        case 38: return 37.976318000;
        case 39: return 38.970719700;
        case 0: case 40: return 39.962590980;
        case 41: return 40.962278060;
        case 42: return 41.958618010;
        case 43: return 42.958766600;
        case 44: return 43.955481800;
        case 45: return 44.956186600;
        case 46: return 45.953692600;
        case 47: return 46.954546000;
        case 48: return 47.952534000;
        case 49: return 48.955674000;
        case 50: return 49.957519000;
        case 51: return 50.961500000;
        case 52: return 51.965100000;
        case 53: return 52.970050000;
        case 54: return 53.974350000;
        case 55: return 54.980550000;
        case 56: return 55.985570000;
        case 57: return 56.992360000;
        }
        break;
      case 21:
        switch (isotope) {
        case 36: return 36.014920000;
        case 37: return 37.003050000;
        case 38: return 37.994700000;
        case 39: return 38.984790000;
        case 40: return 39.977967000;
        case 41: return 40.969251130;
        case 42: return 41.965516430;
        case 43: return 42.961150700;
        case 44: return 43.959402800;
        case 0: case 45: return 44.955911900;
        case 46: return 45.955171900;
        case 47: return 46.952407500;
        case 48: return 47.952231000;
        case 49: return 48.950024000;
        case 50: return 49.952188000;
        case 51: return 50.953603000;
        case 52: return 51.956680000;
        case 53: return 52.959610000;
        case 54: return 53.963260000;
        case 55: return 54.968240000;
        case 56: return 55.972870000;
        case 57: return 56.977790000;
        case 58: return 57.983710000;
        case 59: return 58.989220000;
        case 60: return 59.995710000;
        }
        break;
      case 22:
        switch (isotope) {
        case 38: return 38.009770000;
        case 39: return 39.001610000;
        case 40: return 39.990500000;
        case 41: return 40.983150000;
        case 42: return 41.973031000;
        case 43: return 42.968522000;
        case 44: return 43.959690100;
        case 45: return 44.958125600;
        case 46: return 45.952631600;
        case 47: return 46.951763100;
        case 0: case 48: return 47.947946300;
        case 49: return 48.947870000;
        case 50: return 49.944791200;
        case 51: return 50.946615000;
        case 52: return 51.946897000;
        case 53: return 52.949730000;
        case 54: return 53.951050000;
        case 55: return 54.955270000;
        case 56: return 55.958200000;
        case 57: return 56.963990000;
        case 58: return 57.966970000;
        case 59: return 58.972930000;
        case 60: return 59.976760000;
        case 61: return 60.983200000;
        case 62: return 61.987490000;
        case 63: return 62.994420000;
        }
        break;
      case 23:
        switch (isotope) {
        case 40: return 40.011090000;
        case 41: return 40.999780000;
        case 42: return 41.991230000;
        case 43: return 42.980650000;
        case 44: return 43.974110000;
        case 45: return 44.965776000;
        case 46: return 45.960200500;
        case 47: return 46.954908900;
        case 48: return 47.952253700;
        case 49: return 48.948516100;
        case 50: return 49.947158500;
        case 0: case 51: return 50.943959500;
        case 52: return 51.944775500;
        case 53: return 52.944338000;
        case 54: return 53.946440000;
        case 55: return 54.947230000;
        case 56: return 55.950530000;
        case 57: return 56.952560000;
        case 58: return 57.956830000;
        case 59: return 58.960210000;
        case 60: return 59.965030000;
        case 61: return 60.968480000;
        case 62: return 61.973780000;
        case 63: return 62.977550000;
        case 64: return 63.983470000;
        case 65: return 64.987920000;
        }
        break;
      case 24:
        switch (isotope) {
        case 42: return 42.006430000;
        case 43: return 42.997710000;
        case 44: return 43.985550000;
        case 45: return 44.979640000;
        case 46: return 45.968359000;
        case 47: return 46.962900000;
        case 48: return 47.954032000;
        case 49: return 48.951335700;
        case 50: return 49.946044200;
        case 51: return 50.944767400;
        case 0: case 52: return 51.940507500;
        case 53: return 52.940649400;
        case 54: return 53.938880400;
        case 55: return 54.940839700;
        case 56: return 55.940653100;
        case 57: return 56.943613000;
        case 58: return 57.944350000;
        case 59: return 58.948590000;
        case 60: return 59.950080000;
        case 61: return 60.954720000;
        case 62: return 61.956610000;
        case 63: return 62.961860000;
        case 64: return 63.964410000;
        case 65: return 64.970160000;
        case 66: return 65.973380000;
        case 67: return 66.979550000;
        }
        break;
      case 25:
        switch (isotope) {
        case 44: return 44.006870000;
        case 45: return 44.994510000;
        case 46: return 45.986720000;
        case 47: return 46.976100000;
        case 48: return 47.968520000;
        case 49: return 48.959618000;
        case 50: return 49.954238200;
        case 51: return 50.948210800;
        case 52: return 51.945565500;
        case 53: return 52.941290100;
        case 54: return 53.940358900;
        case 0: case 55: return 54.938045100;
        case 56: return 55.938904900;
        case 57: return 56.938285400;
        case 58: return 57.939980000;
        case 59: return 58.940440000;
        case 60: return 59.942910000;
        case 61: return 60.944650000;
        case 62: return 61.948430000;
        case 63: return 62.950240000;
        case 64: return 63.954250000;
        case 65: return 64.956340000;
        case 66: return 65.961080000;
        case 67: return 66.964140000;
        case 68: return 67.969300000;
        case 69: return 68.972840000;
        }
        break;
      case 26:
        switch (isotope) {
        case 45: return 45.014580000;
        case 46: return 46.000810000;
        case 47: return 46.992890000;
        case 48: return 47.980500000;
        case 49: return 48.973610000;
        case 50: return 49.962990000;
        case 51: return 50.956820000;
        case 52: return 51.948114000;
        case 53: return 52.945307900;
        case 54: return 53.939610500;
        case 55: return 54.938293400;
        case 0: case 56: return 55.934937500;
        case 57: return 56.935394000;
        case 58: return 57.933275600;
        case 59: return 58.934875500;
        case 60: return 59.934072000;
        case 61: return 60.936745000;
        case 62: return 61.936767000;
        case 63: return 62.940370000;
        case 64: return 63.941200000;
        case 65: return 64.945380000;
        case 66: return 65.946780000;
        case 67: return 66.950950000;
        case 68: return 67.953700000;
        case 69: return 68.958780000;
        case 70: return 69.961460000;
        case 71: return 70.966720000;
        case 72: return 71.969620000;
        }
        break;
      case 27:
        switch (isotope) {
        case 47: return 47.011490000;
        case 48: return 48.001760000;
        case 49: return 48.989720000;
        case 50: return 49.981540000;
        case 51: return 50.970720000;
        case 52: return 51.963590000;
        case 53: return 52.954219000;
        case 54: return 53.948459600;
        case 55: return 54.941999000;
        case 56: return 55.939839300;
        case 57: return 56.936291400;
        case 58: return 57.935752800;
        case 0: case 59: return 58.933195000;
        case 60: return 59.933817100;
        case 61: return 60.932475800;
        case 62: return 61.934051000;
        case 63: return 62.933612000;
        case 64: return 63.935810000;
        case 65: return 64.936478000;
        case 66: return 65.939760000;
        case 67: return 66.940890000;
        case 68: return 67.944870000;
        case 69: return 68.946320000;
        case 70: return 69.951000000;
        case 71: return 70.952900000;
        case 72: return 71.957810000;
        case 73: return 72.960240000;
        case 74: return 73.965380000;
        case 75: return 74.968330000;
        }
        break;
      case 28:
        switch (isotope) {
        case 48: return 48.019750000;
        case 49: return 49.009660000;
        case 50: return 49.995930000;
        case 51: return 50.987720000;
        case 52: return 51.975680000;
        case 53: return 52.968470000;
        case 54: return 53.957910000;
        case 55: return 54.951330000;
        case 56: return 55.942132000;
        case 57: return 56.939793500;
        case 0: case 58: return 57.935342900;
        case 59: return 58.934346700;
        case 60: return 59.930786400;
        case 61: return 60.931056000;
        case 62: return 61.928345100;
        case 63: return 62.929669400;
        case 64: return 63.927966000;
        case 65: return 64.930084300;
        case 66: return 65.929139300;
        case 67: return 66.931569000;
        case 68: return 67.931869000;
        case 69: return 68.935610000;
        case 70: return 69.936500000;
        case 71: return 70.940740000;
        case 72: return 71.942090000;
        case 73: return 72.946470000;
        case 74: return 73.948070000;
        case 75: return 74.952870000;
        case 76: return 75.955330000;
        case 77: return 76.960550000;
        case 78: return 77.963180000;
        }
        break;
      case 29:
        switch (isotope) {
        case 52: return 51.997180000;
        case 53: return 52.985550000;
        case 54: return 53.976710000;
        case 55: return 54.966050000;
        case 56: return 55.958560000;
        case 57: return 56.949211000;
        case 58: return 57.944538500;
        case 59: return 58.939498000;
        case 60: return 59.937365000;
        case 61: return 60.933457800;
        case 62: return 61.932584000;
        case 0: case 63: return 62.929597500;
        case 64: return 63.929764200;
        case 65: return 64.927789500;
        case 66: return 65.928868800;
        case 67: return 66.927730300;
        case 68: return 67.929610900;
        case 69: return 68.929429300;
        case 70: return 69.932392300;
        case 71: return 70.932676800;
        case 72: return 71.935820300;
        case 73: return 72.936675000;
        case 74: return 73.939875000;
        case 75: return 74.941900000;
        case 76: return 75.945275000;
        case 77: return 76.947850000;
        case 78: return 77.951960000;
        case 79: return 78.954560000;
        case 80: return 79.960870000;
        }
        break;
      case 30:
        switch (isotope) {
        case 54: return 53.992950000;
        case 55: return 54.983980000;
        case 56: return 55.972380000;
        case 57: return 56.964790000;
        case 58: return 57.954590000;
        case 59: return 58.949260000;
        case 60: return 59.941827000;
        case 61: return 60.939511000;
        case 62: return 61.934330000;
        case 63: return 62.933211000;
        case 0: case 64: return 63.929142000;
        case 65: return 64.929241000;
        case 66: return 65.926033000;
        case 67: return 66.927127000;
        case 68: return 67.924844000;
        case 69: return 68.926550000;
        case 70: return 69.925319000;
        case 71: return 70.927722000;
        case 72: return 71.926858000;
        case 73: return 72.929780000;
        case 74: return 73.929460000;
        case 75: return 74.932940000;
        case 76: return 75.933290000;
        case 77: return 76.936960000;
        case 78: return 77.938440000;
        case 79: return 78.942650000;
        case 80: return 79.944340000;
        case 81: return 80.950480000;
        case 82: return 81.954420000;
        case 83: return 82.961030000;
        }
        break;
      case 31:
        switch (isotope) {
        case 56: return 55.994910000;
        case 57: return 56.982930000;
        case 58: return 57.974250000;
        case 59: return 58.963370000;
        case 60: return 59.957060000;
        case 61: return 60.949450000;
        case 62: return 61.944175000;
        case 63: return 62.939294000;
        case 64: return 63.936838000;
        case 65: return 64.932734000;
        case 66: return 65.931589000;
        case 67: return 66.928201000;
        case 68: return 67.927980000;
        case 0: case 69: return 68.925573000;
        case 70: return 69.926022000;
        case 71: return 70.924701000;
        case 72: return 71.926366000;
        case 73: return 72.925174000;
        case 74: return 73.926946000;
        case 75: return 74.926500000;
        case 76: return 75.928827000;
        case 77: return 76.929154000;
        case 78: return 77.931608000;
        case 79: return 78.932890000;
        case 80: return 79.936520000;
        case 81: return 80.937750000;
        case 82: return 81.942990000;
        case 83: return 82.946980000;
        case 84: return 83.952650000;
        case 85: return 84.957000000;
        case 86: return 85.963120000;
        }
        break;
      case 32:
        switch (isotope) {
        case 58: return 57.991010000;
        case 59: return 58.981750000;
        case 60: return 59.970190000;
        case 61: return 60.963790000;
        case 62: return 61.954650000;
        case 63: return 62.949640000;
        case 64: return 63.941650000;
        case 65: return 64.939440000;
        case 66: return 65.933840000;
        case 67: return 66.932734000;
        case 68: return 67.928094000;
        case 69: return 68.927964000;
        case 70: return 69.924247000;
        case 71: return 70.924951000;
        case 72: return 71.922075000;
        case 73: return 72.923458000;
        case 0: case 74: return 73.921177000;
        case 75: return 74.922858000;
        case 76: return 75.921402000;
        case 77: return 76.923548000;
        case 78: return 77.922853000;
        case 79: return 78.925400000;
        case 80: return 79.925370000;
        case 81: return 80.928820000;
        case 82: return 81.929550000;
        case 83: return 82.934620000;
        case 84: return 83.937470000;
        case 85: return 84.943030000;
        case 86: return 85.946490000;
        case 87: return 86.952510000;
        case 88: return 87.956910000;
        case 89: return 88.963830000;
        }
        break;
      case 33:
        switch (isotope) {
        case 60: return 59.993130000;
        case 61: return 60.980620000;
        case 62: return 61.973200000;
        case 63: return 62.963690000;
        case 64: return 63.957570000;
        case 65: return 64.949560000;
        case 66: return 65.944710000;
        case 67: return 66.939190000;
        case 68: return 67.936770000;
        case 69: return 68.932270000;
        case 70: return 69.930920000;
        case 71: return 70.927112000;
        case 72: return 71.926752000;
        case 73: return 72.923825000;
        case 74: return 73.923928000;
        case 0: case 75: return 74.921596000;
        case 76: return 75.922394000;
        case 77: return 76.920647000;
        case 78: return 77.921827000;
        case 79: return 78.920948000;
        case 80: return 79.922534000;
        case 81: return 80.922132000;
        case 82: return 81.924500000;
        case 83: return 82.924980000;
        case 84: return 83.929060000;
        case 85: return 84.932020000;
        case 86: return 85.936500000;
        case 87: return 86.939900000;
        case 88: return 87.944940000;
        case 89: return 88.949390000;
        case 90: return 89.955500000;
        case 91: return 90.960430000;
        case 92: return 91.966800000;
        }
        break;
      case 34:
        switch (isotope) {
        case 65: return 64.964660000;
        case 66: return 65.955210000;
        case 67: return 66.950090000;
        case 68: return 67.941800000;
        case 69: return 68.939560000;
        case 70: return 69.933390000;
        case 71: return 70.932240000;
        case 72: return 71.927112000;
        case 73: return 72.926765000;
        case 74: return 73.922476000;
        case 75: return 74.922523000;
        case 76: return 75.919213000;
        case 77: return 76.919914000;
        case 78: return 77.917309000;
        case 79: return 78.918499000;
        case 0: case 80: return 79.916521000;
        case 81: return 80.917992000;
        case 82: return 81.916699000;
        case 83: return 82.919118000;
        case 84: return 83.918462000;
        case 85: return 84.922250000;
        case 86: return 85.924272000;
        case 87: return 86.928520000;
        case 88: return 87.931420000;
        case 89: return 88.936450000;
        case 90: return 89.939960000;
        case 91: return 90.945960000;
        case 92: return 91.949920000;
        case 93: return 92.956290000;
        case 94: return 93.960490000;
        }
        break;
      case 35:
        switch (isotope) {
        case 67: return 66.964790000;
        case 68: return 67.958520000;
        case 69: return 68.950110000;
        case 70: return 69.944790000;
        case 71: return 70.938740000;
        case 72: return 71.936640000;
        case 73: return 72.931690000;
        case 74: return 73.929891000;
        case 75: return 74.925776000;
        case 76: return 75.924541000;
        case 77: return 76.921379000;
        case 78: return 77.921146000;
        case 0: case 79: return 78.918337000;
        case 80: return 79.918529000;
        case 81: return 80.916290000;
        case 82: return 81.916804000;
        case 83: return 82.915180000;
        case 84: return 83.916479000;
        case 85: return 84.915608000;
        case 86: return 85.918798000;
        case 87: return 86.920711000;
        case 88: return 87.924070000;
        case 89: return 88.926390000;
        case 90: return 89.930630000;
        case 91: return 90.933970000;
        case 92: return 91.939260000;
        case 93: return 92.943050000;
        case 94: return 93.948680000;
        case 95: return 94.952870000;
        case 96: return 95.958530000;
        case 97: return 96.962800000;
        }
        break;
      case 36:
        switch (isotope) {
        case 69: return 68.965180000;
        case 70: return 69.955260000;
        case 71: return 70.949630000;
        case 72: return 71.942092000;
        case 73: return 72.939289000;
        case 74: return 73.933084000;
        case 75: return 74.930946000;
        case 76: return 75.925910000;
        case 77: return 76.924670000;
        case 78: return 77.920364000;
        case 79: return 78.920082000;
        case 80: return 79.916379000;
        case 81: return 80.916592000;
        case 82: return 81.913483000;
        case 83: return 82.914136000;
        case 0: case 84: return 83.911507000;
        case 85: return 84.912527000;
        case 86: return 85.910610000;
        case 87: return 86.913354000;
        case 88: return 87.914447000;
        case 89: return 88.917630000;
        case 90: return 89.919517000;
        case 91: return 90.923450000;
        case 92: return 91.926156000;
        case 93: return 92.931270000;
        case 94: return 93.934360000;
        case 95: return 94.939840000;
        case 96: return 95.943070000;
        case 97: return 96.948560000;
        case 98: return 97.951910000;
        case 99: return 98.957600000;
        case 100: return 99.961140000;
        }
        break;
      case 37:
        switch (isotope) {
        case 71: return 70.965320000;
        case 72: return 71.959080000;
        case 73: return 72.950560000;
        case 74: return 73.944265000;
        case 75: return 74.938570000;
        case 76: return 75.935072000;
        case 77: return 76.930408000;
        case 78: return 77.928141000;
        case 79: return 78.923989000;
        case 80: return 79.922519000;
        case 81: return 80.918996000;
        case 82: return 81.918208000;
        case 83: return 82.915110000;
        case 84: return 83.914385000;
        case 0: case 85: return 84.911789000;
        case 86: return 85.911167000;
        case 87: return 86.909180000;
        case 88: return 87.911315000;
        case 89: return 88.912278000;
        case 90: return 89.914802000;
        case 91: return 90.916537000;
        case 92: return 91.919729000;
        case 93: return 92.922042000;
        case 94: return 93.926405000;
        case 95: return 94.929303000;
        case 96: return 95.934270000;
        case 97: return 96.937350000;
        case 98: return 97.941790000;
        case 99: return 98.945380000;
        case 100: return 99.949870000;
        case 101: return 100.953200000;
        case 102: return 101.958870000;
        }
        break;
      case 38:
        switch (isotope) {
        case 73: return 72.965970000;
        case 74: return 73.956310000;
        case 75: return 74.949950000;
        case 76: return 75.941770000;
        case 77: return 76.937945000;
        case 78: return 77.932180000;
        case 79: return 78.929708000;
        case 80: return 79.924521000;
        case 81: return 80.923212000;
        case 82: return 81.918402000;
        case 83: return 82.917557000;
        case 84: return 83.913425000;
        case 85: return 84.912933000;
        case 86: return 85.909260000;
        case 87: return 86.908877000;
        case 0: case 88: return 87.905612000;
        case 89: return 88.907450000;
        case 90: return 89.907738000;
        case 91: return 90.910203000;
        case 92: return 91.911038000;
        case 93: return 92.914026000;
        case 94: return 93.915361000;
        case 95: return 94.919359000;
        case 96: return 95.921697000;
        case 97: return 96.926153000;
        case 98: return 97.928453000;
        case 99: return 98.933240000;
        case 100: return 99.935350000;
        case 101: return 100.940520000;
        case 102: return 101.943020000;
        case 103: return 102.948950000;
        case 104: return 103.952330000;
        case 105: return 104.958580000;
        }
        break;
      case 39:
        switch (isotope) {
        case 76: return 75.958450000;
        case 77: return 76.949650000;
        case 78: return 77.943610000;
        case 79: return 78.937350000;
        case 80: return 79.934280000;
        case 81: return 80.929130000;
        case 82: return 81.926790000;
        case 83: return 82.922350000;
        case 84: return 83.920390000;
        case 85: return 84.916433000;
        case 86: return 85.914886000;
        case 87: return 86.910875000;
        case 88: return 87.909501000;
        case 0: case 89: return 88.905848000;
        case 90: return 89.907151000;
        case 91: return 90.907305000;
        case 92: return 91.908949000;
        case 93: return 92.909583000;
        case 94: return 93.911595000;
        case 95: return 94.912821000;
        case 96: return 95.915891000;
        case 97: return 96.918134000;
        case 98: return 97.922203000;
        case 99: return 98.924636000;
        case 100: return 99.927760000;
        case 101: return 100.930310000;
        case 102: return 101.933560000;
        case 103: return 102.936730000;
        case 104: return 103.941050000;
        case 105: return 104.944870000;
        case 106: return 105.949790000;
        case 107: return 106.954140000;
        case 108: return 107.959480000;
        }
        break;
      case 40:
        switch (isotope) {
        case 78: return 77.955230000;
        case 79: return 78.949160000;
        case 80: return 79.940400000;
        case 81: return 80.937210000;
        case 82: return 81.931090000;
        case 83: return 82.928650000;
        case 84: return 83.923250000;
        case 85: return 84.921470000;
        case 86: return 85.916470000;
        case 87: return 86.914816000;
        case 88: return 87.910227000;
        case 89: return 88.908890000;
        case 0: case 90: return 89.904704000;
        case 91: return 90.905645000;
        case 92: return 91.905040000;
        case 93: return 92.906476000;
        case 94: return 93.906315000;
        case 95: return 94.908042000;
        case 96: return 95.908273000;
        case 97: return 96.910953000;
        case 98: return 97.912735000;
        case 99: return 98.916512000;
        case 100: return 99.917760000;
        case 101: return 100.921140000;
        case 102: return 101.922980000;
        case 103: return 102.926600000;
        case 104: return 103.928780000;
        case 105: return 104.933050000;
        case 106: return 105.935910000;
        case 107: return 106.940750000;
        case 108: return 107.943960000;
        case 109: return 108.949240000;
        case 110: return 109.952870000;
        }
        break;
      case 41:
        switch (isotope) {
        case 81: return 80.949030000;
        case 82: return 81.943130000;
        case 83: return 82.936710000;
        case 84: return 83.933570000;
        case 85: return 84.927910000;
        case 86: return 85.925040000;
        case 87: return 86.920360000;
        case 88: return 87.918330000;
        case 89: return 88.913418000;
        case 90: return 89.911265000;
        case 91: return 90.906996000;
        case 92: return 91.907194000;
        case 0: case 93: return 92.906378000;
        case 94: return 93.907283000;
        case 95: return 94.906835000;
        case 96: return 95.908101000;
        case 97: return 96.908098000;
        case 98: return 97.910328000;
        case 99: return 98.911618000;
        case 100: return 99.914182000;
        case 101: return 100.915252000;
        case 102: return 101.918040000;
        case 103: return 102.919140000;
        case 104: return 103.922460000;
        case 105: return 104.923940000;
        case 106: return 105.927970000;
        case 107: return 106.930310000;
        case 108: return 107.934840000;
        case 109: return 108.937630000;
        case 110: return 109.942440000;
        case 111: return 110.945650000;
        case 112: return 111.950830000;
        case 113: return 112.954700000;
        }
        break;
      case 42:
        switch (isotope) {
        case 83: return 82.948740000;
        case 84: return 83.940090000;
        case 85: return 84.936550000;
        case 86: return 85.930700000;
        case 87: return 86.927330000;
        case 88: return 87.921953000;
        case 89: return 88.919480000;
        case 90: return 89.913937000;
        case 91: return 90.911750000;
        case 92: return 91.906811000;
        case 93: return 92.906813000;
        case 94: return 93.905088000;
        case 95: return 94.905842000;
        case 96: return 95.904679000;
        case 97: return 96.906021000;
        case 0: case 98: return 97.905408000;
        case 99: return 98.907711000;
        case 100: return 99.907477000;
        case 101: return 100.910347000;
        case 102: return 101.910297000;
        case 103: return 102.913210000;
        case 104: return 103.913760000;
        case 105: return 104.916970000;
        case 106: return 105.918137000;
        case 107: return 106.921690000;
        case 108: return 107.923450000;
        case 109: return 108.927810000;
        case 110: return 109.929730000;
        case 111: return 110.934410000;
        case 112: return 111.936840000;
        case 113: return 112.941880000;
        case 114: return 113.944920000;
        case 115: return 114.950290000;
        }
        break;
      case 43:
        switch (isotope) {
        case 85: return 84.948830000;
        case 86: return 85.942880000;
        case 87: return 86.936530000;
        case 88: return 87.932680000;
        case 89: return 88.927170000;
        case 90: return 89.923560000;
        case 91: return 90.918430000;
        case 92: return 91.915260000;
        case 93: return 92.910249000;
        case 94: return 93.909657000;
        case 95: return 94.907657000;
        case 96: return 95.907871000;
        case 97: return 96.906365000;
        case 0: case 98: return 97.907216000;
        case 99: return 98.906254000;
        case 100: return 99.907657000;
        case 101: return 100.907315000;
        case 102: return 101.909215000;
        case 103: return 102.909181000;
        case 104: return 103.911450000;
        case 105: return 104.911660000;
        case 106: return 105.914358000;
        case 107: return 106.915080000;
        case 108: return 107.918460000;
        case 109: return 108.919980000;
        case 110: return 109.923820000;
        case 111: return 110.925690000;
        case 112: return 111.929150000;
        case 113: return 112.931590000;
        case 114: return 113.935880000;
        case 115: return 114.938690000;
        case 116: return 115.943370000;
        case 117: return 116.946480000;
        case 118: return 117.951480000;
        }
        break;
      case 44:
        switch (isotope) {
        case 87: return 86.949180000;
        case 88: return 87.940260000;
        case 89: return 88.936110000;
        case 90: return 89.929890000;
        case 91: return 90.926290000;
        case 92: return 91.920120000;
        case 93: return 92.917050000;
        case 94: return 93.911360000;
        case 95: return 94.910413000;
        case 96: return 95.907598000;
        case 97: return 96.907555000;
        case 98: return 97.905287000;
        case 99: return 98.905939000;
        case 100: return 99.904219000;
        case 101: return 100.905582000;
        case 0: case 102: return 101.904349000;
        case 103: return 102.906323000;
        case 104: return 103.905433000;
        case 105: return 104.907753000;
        case 106: return 105.907329000;
        case 107: return 106.909910000;
        case 108: return 107.910170000;
        case 109: return 108.913200000;
        case 110: return 109.914140000;
        case 111: return 110.917700000;
        case 112: return 111.918970000;
        case 113: return 112.922490000;
        case 114: return 113.924280000;
        case 115: return 114.928690000;
        case 116: return 115.930810000;
        case 117: return 116.935580000;
        case 118: return 117.937820000;
        case 119: return 118.942840000;
        case 120: return 119.945310000;
        }
        break;
      case 45:
        switch (isotope) {
        case 89: return 88.948840000;
        case 90: return 89.942870000;
        case 91: return 90.936550000;
        case 92: return 91.931980000;
        case 93: return 92.925740000;
        case 94: return 93.921700000;
        case 95: return 94.915900000;
        case 96: return 95.914461000;
        case 97: return 96.911340000;
        case 98: return 97.910708000;
        case 99: return 98.908132000;
        case 100: return 99.908122000;
        case 101: return 100.906164000;
        case 102: return 101.906843000;
        case 0: case 103: return 102.905504000;
        case 104: return 103.906656000;
        case 105: return 104.905694000;
        case 106: return 105.907287000;
        case 107: return 106.906748000;
        case 108: return 107.908730000;
        case 109: return 108.908737000;
        case 110: return 109.911140000;
        case 111: return 110.911590000;
        case 112: return 111.914390000;
        case 113: return 112.915530000;
        case 114: return 113.918810000;
        case 115: return 114.920330000;
        case 116: return 115.924060000;
        case 117: return 116.925980000;
        case 118: return 117.930070000;
        case 119: return 118.932110000;
        case 120: return 119.936410000;
        case 121: return 120.938720000;
        case 122: return 121.943210000;
        }
        break;
      case 46:
        switch (isotope) {
        case 91: return 90.949110000;
        case 92: return 91.940420000;
        case 93: return 92.935910000;
        case 94: return 93.928770000;
        case 95: return 94.924690000;
        case 96: return 95.918160000;
        case 97: return 96.916480000;
        case 98: return 97.912721000;
        case 99: return 98.911768000;
        case 100: return 99.908506000;
        case 101: return 100.908289000;
        case 102: return 101.905609000;
        case 103: return 102.906087000;
        case 104: return 103.904036000;
        case 105: return 104.905085000;
        case 0: case 106: return 105.903486000;
        case 107: return 106.905133000;
        case 108: return 107.903892000;
        case 109: return 108.905950000;
        case 110: return 109.905153000;
        case 111: return 110.907671000;
        case 112: return 111.907314000;
        case 113: return 112.910150000;
        case 114: return 113.910363000;
        case 115: return 114.913680000;
        case 116: return 115.914160000;
        case 117: return 116.917840000;
        case 118: return 117.918980000;
        case 119: return 118.923110000;
        case 120: return 119.924690000;
        case 121: return 120.928870000;
        case 122: return 121.930550000;
        case 123: return 122.934930000;
        case 124: return 123.936880000;
        }
        break;
      case 47:
        switch (isotope) {
        case 93: return 92.949780000;
        case 94: return 93.942780000;
        case 95: return 94.935480000;
        case 96: return 95.930680000;
        case 97: return 96.923970000;
        case 98: return 97.921570000;
        case 99: return 98.917600000;
        case 100: return 99.916100000;
        case 101: return 100.912800000;
        case 102: return 101.911690000;
        case 103: return 102.908973000;
        case 104: return 103.908629000;
        case 105: return 104.906529000;
        case 106: return 105.906669000;
        case 0: case 107: return 106.905097000;
        case 108: return 107.905956000;
        case 109: return 108.904752000;
        case 110: return 109.906107000;
        case 111: return 110.905291000;
        case 112: return 111.907005000;
        case 113: return 112.906567000;
        case 114: return 113.908804000;
        case 115: return 114.908760000;
        case 116: return 115.911360000;
        case 117: return 116.911680000;
        case 118: return 117.914580000;
        case 119: return 118.915670000;
        case 120: return 119.918790000;
        case 121: return 120.919850000;
        case 122: return 121.923530000;
        case 123: return 122.924900000;
        case 124: return 123.928640000;
        case 125: return 124.930430000;
        case 126: return 125.934500000;
        case 127: return 126.936770000;
        case 128: return 127.941170000;
        case 129: return 128.943690000;
        case 130: return 129.950450000;
        }
        break;
      case 48:
        switch (isotope) {
        case 95: return 94.949870000;
        case 96: return 95.939770000;
        case 97: return 96.934940000;
        case 98: return 97.927400000;
        case 99: return 98.925010000;
        case 100: return 99.920290000;
        case 101: return 100.918680000;
        case 102: return 101.914460000;
        case 103: return 102.913419000;
        case 104: return 103.909849000;
        case 105: return 104.909468000;
        case 106: return 105.906459000;
        case 107: return 106.906618000;
        case 108: return 107.904184000;
        case 109: return 108.904982000;
        case 110: return 109.903002000;
        case 111: return 110.904178000;
        case 112: return 111.902757000;
        case 113: return 112.904401000;
        case 0: case 114: return 113.903358000;
        case 115: return 114.905431000;
        case 116: return 115.904756000;
        case 117: return 116.907219000;
        case 118: return 117.906915000;
        case 119: return 118.909920000;
        case 120: return 119.909850000;
        case 121: return 120.912980000;
        case 122: return 121.913330000;
        case 123: return 122.917000000;
        case 124: return 123.917650000;
        case 125: return 124.921250000;
        case 126: return 125.922350000;
        case 127: return 126.926440000;
        case 128: return 127.927760000;
        case 129: return 128.932150000;
        case 130: return 129.933900000;
        case 131: return 130.940670000;
        case 132: return 131.945550000;
        }
        break;
      case 49:
        switch (isotope) {
        case 97: return 96.949540000;
        case 98: return 97.942140000;
        case 99: return 98.934220000;
        case 100: return 99.931110000;
        case 101: return 100.926340000;
        case 102: return 101.924090000;
        case 103: return 102.919914000;
        case 104: return 103.918300000;
        case 105: return 104.914674000;
        case 106: return 105.913465000;
        case 107: return 106.910295000;
        case 108: return 107.909698000;
        case 109: return 108.907151000;
        case 110: return 109.907165000;
        case 111: return 110.905103000;
        case 112: return 111.905532000;
        case 113: return 112.904058000;
        case 114: return 113.904914000;
        case 0: case 115: return 114.903878000;
        case 116: return 115.905260000;
        case 117: return 116.904514000;
        case 118: return 117.906354000;
        case 119: return 118.905845000;
        case 120: return 119.907960000;
        case 121: return 120.907846000;
        case 122: return 121.910280000;
        case 123: return 122.910438000;
        case 124: return 123.913180000;
        case 125: return 124.913600000;
        case 126: return 125.916460000;
        case 127: return 126.917350000;
        case 128: return 127.920170000;
        case 129: return 128.921700000;
        case 130: return 129.924970000;
        case 131: return 130.926850000;
        case 132: return 131.932990000;
        case 133: return 132.937810000;
        case 134: return 133.944150000;
        case 135: return 134.949330000;
        }
        break;
      case 50:
        switch (isotope) {
        case 99: return 98.949330000;
        case 100: return 99.939040000;
        case 101: return 100.936060000;
        case 102: return 101.930300000;
        case 103: return 102.928100000;
        case 104: return 103.923140000;
        case 105: return 104.921350000;
        case 106: return 105.916880000;
        case 107: return 106.915640000;
        case 108: return 107.911925000;
        case 109: return 108.911283000;
        case 110: return 109.907843000;
        case 111: return 110.907734000;
        case 112: return 111.904818000;
        case 113: return 112.905171000;
        case 114: return 113.902779000;
        case 115: return 114.903342000;
        case 116: return 115.901741000;
        case 117: return 116.902952000;
        case 118: return 117.901603000;
        case 119: return 118.903308000;
        case 0: case 120: return 119.902194000;
        case 121: return 120.904235000;
        case 122: return 121.903439000;
        case 123: return 122.905720000;
        case 124: return 123.905273000;
        case 125: return 124.907784000;
        case 126: return 125.907653000;
        case 127: return 126.910360000;
        case 128: return 127.910537000;
        case 129: return 128.913480000;
        case 130: return 129.913967000;
        case 131: return 130.917000000;
        case 132: return 131.917816000;
        case 133: return 132.923830000;
        case 134: return 133.928290000;
        case 135: return 134.934730000;
        case 136: return 135.939340000;
        case 137: return 136.945990000;
        }
        break;
      case 51:
        switch (isotope) {
        case 103: return 102.939690000;
        case 104: return 103.936470000;
        case 105: return 104.931490000;
        case 106: return 105.928790000;
        case 107: return 106.924150000;
        case 108: return 107.922160000;
        case 109: return 108.918132000;
        case 110: return 109.916750000;
        case 111: return 110.913160000;
        case 112: return 111.912398000;
        case 113: return 112.909372000;
        case 114: return 113.909270000;
        case 115: return 114.906598000;
        case 116: return 115.906794000;
        case 117: return 116.904836000;
        case 118: return 117.905529000;
        case 119: return 118.903942000;
        case 120: return 119.905072000;
        case 0: case 121: return 120.903815000;
        case 122: return 121.905173000;
        case 123: return 122.904214000;
        case 124: return 123.905935000;
        case 125: return 124.905253000;
        case 126: return 125.907250000;
        case 127: return 126.906924000;
        case 128: return 127.909169000;
        case 129: return 128.909148000;
        case 130: return 129.911656000;
        case 131: return 130.911982000;
        case 132: return 131.914467000;
        case 133: return 132.915252000;
        case 134: return 133.920380000;
        case 135: return 134.925170000;
        case 136: return 135.930350000;
        case 137: return 136.935310000;
        case 138: return 137.940790000;
        case 139: return 138.945980000;
        }
        break;
      case 52:
        switch (isotope) {
        case 105: return 104.943640000;
        case 106: return 105.937500000;
        case 107: return 106.935010000;
        case 108: return 107.929440000;
        case 109: return 108.927420000;
        case 110: return 109.922410000;
        case 111: return 110.921110000;
        case 112: return 111.917010000;
        case 113: return 112.915890000;
        case 114: return 113.912090000;
        case 115: return 114.911900000;
        case 116: return 115.908460000;
        case 117: return 116.908645000;
        case 118: return 117.905828000;
        case 119: return 118.906404000;
        case 120: return 119.904020000;
        case 121: return 120.904936000;
        case 122: return 121.903043000;
        case 123: return 122.904270000;
        case 124: return 123.902817000;
        case 125: return 124.904430000;
        case 126: return 125.903311000;
        case 127: return 126.905226000;
        case 128: return 127.904463000;
        case 129: return 128.906598000;
        case 0: case 130: return 129.906224000;
        case 131: return 130.908523000;
        case 132: return 131.908553000;
        case 133: return 132.910955000;
        case 134: return 133.911369000;
        case 135: return 134.916450000;
        case 136: return 135.920100000;
        case 137: return 136.925320000;
        case 138: return 137.929220000;
        case 139: return 138.934730000;
        case 140: return 139.938850000;
        case 141: return 140.944650000;
        case 142: return 141.949080000;
        }
        break;
      case 53:
        switch (isotope) {
        case 108: return 107.943480000;
        case 109: return 108.938150000;
        case 110: return 109.935240000;
        case 111: return 110.930280000;
        case 112: return 111.927970000;
        case 113: return 112.923640000;
        case 114: return 113.921850000;
        case 115: return 114.918050000;
        case 116: return 115.916810000;
        case 117: return 116.913650000;
        case 118: return 117.913074000;
        case 119: return 118.910070000;
        case 120: return 119.910048000;
        case 121: return 120.907367000;
        case 122: return 121.907589000;
        case 123: return 122.905589000;
        case 124: return 123.906209000;
        case 125: return 124.904630000;
        case 126: return 125.905624000;
        case 0: case 127: return 126.904473000;
        case 128: return 127.905809000;
        case 129: return 128.904988000;
        case 130: return 129.906674000;
        case 131: return 130.906124000;
        case 132: return 131.907997000;
        case 133: return 132.907797000;
        case 134: return 133.909744000;
        case 135: return 134.910048000;
        case 136: return 135.914650000;
        case 137: return 136.917871000;
        case 138: return 137.922350000;
        case 139: return 138.926100000;
        case 140: return 139.931000000;
        case 141: return 140.935030000;
        case 142: return 141.940180000;
        case 143: return 142.944560000;
        case 144: return 143.949990000;
        }
        break;
      case 54:
        switch (isotope) {
        case 110: return 109.944280000;
        case 111: return 110.941600000;
        case 112: return 111.935620000;
        case 113: return 112.933340000;
        case 114: return 113.927980000;
        case 115: return 114.926294000;
        case 116: return 115.921581000;
        case 117: return 116.920359000;
        case 118: return 117.916179000;
        case 119: return 118.915411000;
        case 120: return 119.911784000;
        case 121: return 120.911462000;
        case 122: return 121.908368000;
        case 123: return 122.908482000;
        case 124: return 123.905893000;
        case 125: return 124.906395000;
        case 126: return 125.904274000;
        case 127: return 126.905184000;
        case 128: return 127.903531000;
        case 129: return 128.904779000;
        case 130: return 129.903508000;
        case 131: return 130.905082000;
        case 0: case 132: return 131.904153000;
        case 133: return 132.905910000;
        case 134: return 133.905394000;
        case 135: return 134.907227000;
        case 136: return 135.907219000;
        case 137: return 136.911562000;
        case 138: return 137.913950000;
        case 139: return 138.918793000;
        case 140: return 139.921640000;
        case 141: return 140.926650000;
        case 142: return 141.929710000;
        case 143: return 142.935110000;
        case 144: return 143.938510000;
        case 145: return 144.944070000;
        case 146: return 145.947750000;
        case 147: return 146.953560000;
        }
        break;
      case 55:
        switch (isotope) {
        case 112: return 111.950300000;
        case 113: return 112.944490000;
        case 114: return 113.941450000;
        case 115: return 114.935910000;
        case 116: return 115.933370000;
        case 117: return 116.928670000;
        case 118: return 117.926559000;
        case 119: return 118.922377000;
        case 120: return 119.920677000;
        case 121: return 120.917229000;
        case 122: return 121.916110000;
        case 123: return 122.912996000;
        case 124: return 123.912258000;
        case 125: return 124.909728000;
        case 126: return 125.909452000;
        case 127: return 126.907418000;
        case 128: return 127.907749000;
        case 129: return 128.906064000;
        case 130: return 129.906709000;
        case 131: return 130.905464000;
        case 132: return 131.906434000;
        case 0: case 133: return 132.905451000;
        case 134: return 133.906718000;
        case 135: return 134.905977000;
        case 136: return 135.907311000;
        case 137: return 136.907089000;
        case 138: return 137.911017000;
        case 139: return 138.913364000;
        case 140: return 139.917282000;
        case 141: return 140.920046000;
        case 142: return 141.924299000;
        case 143: return 142.927352000;
        case 144: return 143.932077000;
        case 145: return 144.935526000;
        case 146: return 145.940290000;
        case 147: return 146.944160000;
        case 148: return 147.949220000;
        case 149: return 148.952930000;
        case 150: return 149.958170000;
        case 151: return 150.962190000;
        }
        break;
      case 56:
        switch (isotope) {
        case 114: return 113.950680000;
        case 115: return 114.947370000;
        case 116: return 115.941380000;
        case 117: return 116.938500000;
        case 118: return 117.933040000;
        case 119: return 118.930660000;
        case 120: return 119.926040000;
        case 121: return 120.924050000;
        case 122: return 121.919900000;
        case 123: return 122.918781000;
        case 124: return 123.915094000;
        case 125: return 124.914473000;
        case 126: return 125.911250000;
        case 127: return 126.911094000;
        case 128: return 127.908318000;
        case 129: return 128.908679000;
        case 130: return 129.906320000;
        case 131: return 130.906941000;
        case 132: return 131.905061000;
        case 133: return 132.906007000;
        case 134: return 133.904508000;
        case 135: return 134.905688000;
        case 136: return 135.904575000;
        case 137: return 136.905827000;
        case 0: case 138: return 137.905247000;
        case 139: return 138.908841000;
        case 140: return 139.910605000;
        case 141: return 140.914411000;
        case 142: return 141.916453000;
        case 143: return 142.920627000;
        case 144: return 143.922953000;
        case 145: return 144.927630000;
        case 146: return 145.930220000;
        case 147: return 146.934950000;
        case 148: return 147.937720000;
        case 149: return 148.942580000;
        case 150: return 149.945680000;
        case 151: return 150.950810000;
        case 152: return 151.954270000;
        case 153: return 152.959610000;
        }
        break;
      case 57:
        switch (isotope) {
        case 117: return 116.950070000;
        case 118: return 117.946730000;
        case 119: return 118.940990000;
        case 120: return 119.938070000;
        case 121: return 120.933010000;
        case 122: return 121.930710000;
        case 123: return 122.926240000;
        case 124: return 123.924570000;
        case 125: return 124.920816000;
        case 126: return 125.919510000;
        case 127: return 126.916375000;
        case 128: return 127.915590000;
        case 129: return 128.912693000;
        case 130: return 129.912369000;
        case 131: return 130.910070000;
        case 132: return 131.910100000;
        case 133: return 132.908220000;
        case 134: return 133.908514000;
        case 135: return 134.906977000;
        case 136: return 135.907640000;
        case 137: return 136.906494000;
        case 138: return 137.907112000;
        case 0: case 139: return 138.906353000;
        case 140: return 139.909477000;
        case 141: return 140.910962000;
        case 142: return 141.914079000;
        case 143: return 142.916063000;
        case 144: return 143.919600000;
        case 145: return 144.921650000;
        case 146: return 145.925790000;
        case 147: return 146.928240000;
        case 148: return 147.932230000;
        case 149: return 148.934730000;
        case 150: return 149.938770000;
        case 151: return 150.941720000;
        case 152: return 151.946250000;
        case 153: return 152.949620000;
        case 154: return 153.954500000;
        case 155: return 154.958350000;
        }
        break;
      case 58:
        switch (isotope) {
        case 119: return 118.952760000;
        case 120: return 119.946640000;
        case 121: return 120.943420000;
        case 122: return 121.937910000;
        case 123: return 122.935400000;
        case 124: return 123.930410000;
        case 125: return 124.928440000;
        case 126: return 125.923970000;
        case 127: return 126.922730000;
        case 128: return 127.918910000;
        case 129: return 128.918100000;
        case 130: return 129.914740000;
        case 131: return 130.914420000;
        case 132: return 131.911460000;
        case 133: return 132.911515000;
        case 134: return 133.908925000;
        case 135: return 134.909151000;
        case 136: return 135.907172000;
        case 137: return 136.907806000;
        case 138: return 137.905991000;
        case 139: return 138.906653000;
        case 0: case 140: return 139.905438000;
        case 141: return 140.908276000;
        case 142: return 141.909244000;
        case 143: return 142.912386000;
        case 144: return 143.913647000;
        case 145: return 144.917230000;
        case 146: return 145.918760000;
        case 147: return 146.922670000;
        case 148: return 147.924430000;
        case 149: return 148.928400000;
        case 150: return 149.930410000;
        case 151: return 150.933980000;
        case 152: return 151.936540000;
        case 153: return 152.940580000;
        case 154: return 153.943420000;
        case 155: return 154.948040000;
        case 156: return 155.951260000;
        case 157: return 156.956340000;
        }
        break;
      case 59:
        switch (isotope) {
        case 121: return 120.955360000;
        case 122: return 121.951810000;
        case 123: return 122.945960000;
        case 124: return 123.942960000;
        case 125: return 124.937830000;
        case 126: return 125.935310000;
        case 127: return 126.930830000;
        case 128: return 127.928790000;
        case 129: return 128.925100000;
        case 130: return 129.923590000;
        case 131: return 130.920260000;
        case 132: return 131.919260000;
        case 133: return 132.916331000;
        case 134: return 133.915710000;
        case 135: return 134.913112000;
        case 136: return 135.912692000;
        case 137: return 136.910705000;
        case 138: return 137.910755000;
        case 139: return 138.908938000;
        case 140: return 139.909076000;
        case 0: case 141: return 140.907652000;
        case 142: return 141.910044000;
        case 143: return 142.910816000;
        case 144: return 143.913305000;
        case 145: return 144.914512000;
        case 146: return 145.917640000;
        case 147: return 146.918996000;
        case 148: return 147.922135000;
        case 149: return 148.923720000;
        case 150: return 149.926673000;
        case 151: return 150.928319000;
        case 152: return 151.931500000;
        case 153: return 152.933840000;
        case 154: return 153.937520000;
        case 155: return 154.940120000;
        case 156: return 155.944270000;
        case 157: return 156.947430000;
        case 158: return 157.951980000;
        case 159: return 158.955500000;
        }
        break;
      case 60:
        switch (isotope) {
        case 124: return 123.952230000;
        case 125: return 124.948880000;
        case 126: return 125.943220000;
        case 127: return 126.940500000;
        case 128: return 127.935390000;
        case 129: return 128.933190000;
        case 130: return 129.928510000;
        case 131: return 130.927250000;
        case 132: return 131.923321000;
        case 133: return 132.922350000;
        case 134: return 133.918790000;
        case 135: return 134.918181000;
        case 136: return 135.914976000;
        case 137: return 136.914567000;
        case 138: return 137.911950000;
        case 139: return 138.911978000;
        case 140: return 139.909550000;
        case 141: return 140.909610000;
        case 0: case 142: return 141.907723000;
        case 143: return 142.909814000;
        case 144: return 143.910087000;
        case 145: return 144.912573000;
        case 146: return 145.913116000;
        case 147: return 146.916100000;
        case 148: return 147.916893000;
        case 149: return 148.920149000;
        case 150: return 149.920891000;
        case 151: return 150.923829000;
        case 152: return 151.924682000;
        case 153: return 152.927698000;
        case 154: return 153.929480000;
        case 155: return 154.932930000;
        case 156: return 155.935020000;
        case 157: return 156.939030000;
        case 158: return 157.941600000;
        case 159: return 158.946090000;
        case 160: return 159.949090000;
        case 161: return 160.953880000;
        }
        break;
      case 61:
        switch (isotope) {
        case 126: return 125.957520000;
        case 127: return 126.951630000;
        case 128: return 127.948420000;
        case 129: return 128.943160000;
        case 130: return 129.940450000;
        case 131: return 130.935870000;
        case 132: return 131.933750000;
        case 133: return 132.929780000;
        case 134: return 133.928350000;
        case 135: return 134.924880000;
        case 136: return 135.923570000;
        case 137: return 136.920479000;
        case 138: return 137.919548000;
        case 139: return 138.916804000;
        case 140: return 139.916040000;
        case 141: return 140.913555000;
        case 142: return 141.912874000;
        case 143: return 142.910933000;
        case 144: return 143.912591000;
        case 0: case 145: return 144.912749000;
        case 146: return 145.914696000;
        case 147: return 146.915138000;
        case 148: return 147.917475000;
        case 149: return 148.918334000;
        case 150: return 149.920984000;
        case 151: return 150.921207000;
        case 152: return 151.923497000;
        case 153: return 152.924117000;
        case 154: return 153.926460000;
        case 155: return 154.928100000;
        case 156: return 155.931060000;
        case 157: return 156.933040000;
        case 158: return 157.936560000;
        case 159: return 158.938970000;
        case 160: return 159.942990000;
        case 161: return 160.945860000;
        case 162: return 161.950290000;
        case 163: return 162.953680000;
        }
        break;
      case 62:
        switch (isotope) {
        case 128: return 127.958080000;
        case 129: return 128.954640000;
        case 130: return 129.948920000;
        case 131: return 130.946110000;
        case 132: return 131.940690000;
        case 133: return 132.938670000;
        case 134: return 133.933970000;
        case 135: return 134.932520000;
        case 136: return 135.928276000;
        case 137: return 136.926970000;
        case 138: return 137.923244000;
        case 139: return 138.922297000;
        case 140: return 139.918995000;
        case 141: return 140.918476000;
        case 142: return 141.915198000;
        case 143: return 142.914628000;
        case 144: return 143.911999000;
        case 145: return 144.913410000;
        case 146: return 145.913041000;
        case 147: return 146.914897000;
        case 148: return 147.914822000;
        case 149: return 148.917184000;
        case 150: return 149.917275000;
        case 151: return 150.919932000;
        case 0: case 152: return 151.919732000;
        case 153: return 152.922097000;
        case 154: return 153.922209000;
        case 155: return 154.924640000;
        case 156: return 155.925528000;
        case 157: return 156.928360000;
        case 158: return 157.929990000;
        case 159: return 158.933210000;
        case 160: return 159.935140000;
        case 161: return 160.938830000;
        case 162: return 161.941220000;
        case 163: return 162.945360000;
        case 164: return 163.948280000;
        case 165: return 164.952980000;
        }
        break;
      case 63:
        switch (isotope) {
        case 130: return 129.963570000;
        case 131: return 130.957750000;
        case 132: return 131.954370000;
        case 133: return 132.949240000;
        case 134: return 133.946510000;
        case 135: return 134.941820000;
        case 136: return 135.939600000;
        case 137: return 136.935570000;
        case 138: return 137.933710000;
        case 139: return 138.929792000;
        case 140: return 139.928090000;
        case 141: return 140.924931000;
        case 142: return 141.923430000;
        case 143: return 142.920298000;
        case 144: return 143.918817000;
        case 145: return 144.916265000;
        case 146: return 145.917206000;
        case 147: return 146.916746000;
        case 148: return 147.918086000;
        case 149: return 148.917931000;
        case 150: return 149.919702000;
        case 151: return 150.919850000;
        case 152: return 151.921744000;
        case 0: case 153: return 152.921230000;
        case 154: return 153.922979000;
        case 155: return 154.922893000;
        case 156: return 155.924752000;
        case 157: return 156.925424000;
        case 158: return 157.927850000;
        case 159: return 158.929089000;
        case 160: return 159.931970000;
        case 161: return 160.933680000;
        case 162: return 161.937040000;
        case 163: return 162.939210000;
        case 164: return 163.942990000;
        case 165: return 164.945720000;
        case 166: return 165.949970000;
        case 167: return 166.953210000;
        }
        break;
      case 64:
        switch (isotope) {
        case 134: return 133.955370000;
        case 135: return 134.952570000;
        case 136: return 135.947340000;
        case 137: return 136.945020000;
        case 138: return 137.940120000;
        case 139: return 138.938240000;
        case 140: return 139.933670000;
        case 141: return 140.932126000;
        case 142: return 141.928120000;
        case 143: return 142.926750000;
        case 144: return 143.922960000;
        case 145: return 144.921709000;
        case 146: return 145.918311000;
        case 147: return 146.919094000;
        case 148: return 147.918115000;
        case 149: return 148.919341000;
        case 150: return 149.918659000;
        case 151: return 150.920348000;
        case 152: return 151.919791000;
        case 153: return 152.921749000;
        case 154: return 153.920865000;
        case 155: return 154.922622000;
        case 156: return 155.922122000;
        case 157: return 156.923960000;
        case 0: case 158: return 157.924103000;
        case 159: return 158.926388000;
        case 160: return 159.927054000;
        case 161: return 160.929669000;
        case 162: return 161.930985000;
        case 163: return 162.933990000;
        case 164: return 163.935860000;
        case 165: return 164.939380000;
        case 166: return 165.941600000;
        case 167: return 166.945570000;
        case 168: return 167.948360000;
        case 169: return 168.952870000;
        }
        break;
      case 65:
        switch (isotope) {
        case 136: return 135.961380000;
        case 137: return 136.955980000;
        case 138: return 137.953160000;
        case 139: return 138.948290000;
        case 140: return 139.945810000;
        case 141: return 140.941450000;
        case 142: return 141.938740000;
        case 143: return 142.935120000;
        case 144: return 143.933050000;
        case 145: return 144.929270000;
        case 146: return 145.927250000;
        case 147: return 146.924045000;
        case 148: return 147.924272000;
        case 149: return 148.923246000;
        case 150: return 149.923660000;
        case 151: return 150.923103000;
        case 152: return 151.924070000;
        case 153: return 152.923435000;
        case 154: return 153.924680000;
        case 155: return 154.923505000;
        case 156: return 155.924747000;
        case 157: return 156.924024000;
        case 158: return 157.925413000;
        case 0: case 159: return 158.925346000;
        case 160: return 159.927167000;
        case 161: return 160.927569000;
        case 162: return 161.929490000;
        case 163: return 162.930648000;
        case 164: return 163.933350000;
        case 165: return 164.934880000;
        case 166: return 165.937990000;
        case 167: return 166.940050000;
        case 168: return 167.943640000;
        case 169: return 168.946220000;
        case 170: return 169.950250000;
        case 171: return 170.953300000;
        }
        break;
      case 66:
        switch (isotope) {
        case 138: return 137.962490000;
        case 139: return 138.959540000;
        case 140: return 139.954010000;
        case 141: return 140.951350000;
        case 142: return 141.946370000;
        case 143: return 142.943830000;
        case 144: return 143.939250000;
        case 145: return 144.937430000;
        case 146: return 145.932845000;
        case 147: return 146.931092000;
        case 148: return 147.927150000;
        case 149: return 148.927305000;
        case 150: return 149.925585000;
        case 151: return 150.926185000;
        case 152: return 151.924718000;
        case 153: return 152.925765000;
        case 154: return 153.924424000;
        case 155: return 154.925754000;
        case 156: return 155.924283000;
        case 157: return 156.925466000;
        case 158: return 157.924409000;
        case 159: return 158.925739000;
        case 160: return 159.925197000;
        case 161: return 160.926933000;
        case 162: return 161.926798000;
        case 163: return 162.928731000;
        case 0: case 164: return 163.929174000;
        case 165: return 164.931703000;
        case 166: return 165.932806000;
        case 167: return 166.935660000;
        case 168: return 167.937130000;
        case 169: return 168.940310000;
        case 170: return 169.942390000;
        case 171: return 170.946200000;
        case 172: return 171.948760000;
        case 173: return 172.953000000;
        }
        break;
      case 67:
        switch (isotope) {
        case 140: return 139.968540000;
        case 141: return 140.963100000;
        case 142: return 141.959770000;
        case 143: return 142.954610000;
        case 144: return 143.951480000;
        case 145: return 144.947200000;
        case 146: return 145.944640000;
        case 147: return 146.940060000;
        case 148: return 147.937720000;
        case 149: return 148.933775000;
        case 150: return 149.933496000;
        case 151: return 150.931688000;
        case 152: return 151.931714000;
        case 153: return 152.930199000;
        case 154: return 153.930602000;
        case 155: return 154.929103000;
        case 156: return 155.929840000;
        case 157: return 156.928256000;
        case 158: return 157.928941000;
        case 159: return 158.927712000;
        case 160: return 159.928729000;
        case 161: return 160.927855000;
        case 162: return 161.929096000;
        case 163: return 162.928733000;
        case 164: return 163.930233000;
        case 0: case 165: return 164.930322000;
        case 166: return 165.932284000;
        case 167: return 166.933133000;
        case 168: return 167.935520000;
        case 169: return 168.936872000;
        case 170: return 169.939620000;
        case 171: return 170.941470000;
        case 172: return 171.944820000;
        case 173: return 172.947290000;
        case 174: return 173.951150000;
        case 175: return 174.954050000;
        }
        break;
      case 68:
        switch (isotope) {
        case 143: return 142.966340000;
        case 144: return 143.960380000;
        case 145: return 144.957390000;
        case 146: return 145.952000000;
        case 147: return 146.949490000;
        case 148: return 147.944550000;
        case 149: return 148.942310000;
        case 150: return 149.937914000;
        case 151: return 150.937449000;
        case 152: return 151.935050000;
        case 153: return 152.935063000;
        case 154: return 153.932783000;
        case 155: return 154.933209000;
        case 156: return 155.931065000;
        case 157: return 156.931920000;
        case 158: return 157.929893000;
        case 159: return 158.930684000;
        case 160: return 159.929083000;
        case 161: return 160.929995000;
        case 162: return 161.928778000;
        case 163: return 162.930033000;
        case 164: return 163.929200000;
        case 165: return 164.930726000;
        case 0: case 166: return 165.930293000;
        case 167: return 166.932048000;
        case 168: return 167.932370000;
        case 169: return 168.934590000;
        case 170: return 169.935464000;
        case 171: return 170.938029000;
        case 172: return 171.939356000;
        case 173: return 172.942400000;
        case 174: return 173.944230000;
        case 175: return 174.947770000;
        case 176: return 175.950080000;
        case 177: return 176.954050000;
        }
        break;
      case 69:
        switch (isotope) {
        case 145: return 144.970070000;
        case 146: return 145.966430000;
        case 147: return 146.960960000;
        case 148: return 147.957840000;
        case 149: return 148.952720000;
        case 150: return 149.949960000;
        case 151: return 150.945483000;
        case 152: return 151.944420000;
        case 153: return 152.942012000;
        case 154: return 153.941568000;
        case 155: return 154.939199000;
        case 156: return 155.938980000;
        case 157: return 156.936970000;
        case 158: return 157.936980000;
        case 159: return 158.934980000;
        case 160: return 159.935260000;
        case 161: return 160.933550000;
        case 162: return 161.933995000;
        case 163: return 162.932651000;
        case 164: return 163.933560000;
        case 165: return 164.932435000;
        case 166: return 165.933554000;
        case 167: return 166.932851000;
        case 168: return 167.934173000;
        case 0: case 169: return 168.934213000;
        case 170: return 169.935801000;
        case 171: return 170.936429000;
        case 172: return 171.938400000;
        case 173: return 172.939604000;
        case 174: return 173.942170000;
        case 175: return 174.943840000;
        case 176: return 175.946990000;
        case 177: return 176.949040000;
        case 178: return 177.952640000;
        case 179: return 178.955340000;
        }
        break;
      case 70:
        switch (isotope) {
        case 148: return 147.967420000;
        case 149: return 148.964040000;
        case 150: return 149.958420000;
        case 151: return 150.955400000;
        case 152: return 151.950290000;
        case 153: return 152.949480000;
        case 154: return 153.946394000;
        case 155: return 154.945782000;
        case 156: return 155.942818000;
        case 157: return 156.942628000;
        case 158: return 157.939866000;
        case 159: return 158.940050000;
        case 160: return 159.937552000;
        case 161: return 160.937902000;
        case 162: return 161.935768000;
        case 163: return 162.936334000;
        case 164: return 163.934489000;
        case 165: return 164.935280000;
        case 166: return 165.933882000;
        case 167: return 166.934950000;
        case 168: return 167.933897000;
        case 169: return 168.935190000;
        case 170: return 169.934761000;
        case 171: return 170.936325000;
        case 172: return 171.936381000;
        case 173: return 172.938210000;
        case 0: case 174: return 173.938862000;
        case 175: return 174.941276000;
        case 176: return 175.942571000;
        case 177: return 176.945260000;
        case 178: return 177.946647000;
        case 179: return 178.950170000;
        case 180: return 179.952330000;
        case 181: return 180.956150000;
        }
        break;
      case 71:
        switch (isotope) {
        case 150: return 149.973230000;
        case 151: return 150.967580000;
        case 152: return 151.964120000;
        case 153: return 152.958770000;
        case 154: return 153.957520000;
        case 155: return 154.954316000;
        case 156: return 155.953030000;
        case 157: return 156.950098000;
        case 158: return 157.949313000;
        case 159: return 158.946630000;
        case 160: return 159.946030000;
        case 161: return 160.943570000;
        case 162: return 161.943280000;
        case 163: return 162.941180000;
        case 164: return 163.941340000;
        case 165: return 164.939407000;
        case 166: return 165.939860000;
        case 167: return 166.938270000;
        case 168: return 167.938740000;
        case 169: return 168.937651000;
        case 170: return 169.938475000;
        case 171: return 170.937913000;
        case 172: return 171.939086000;
        case 173: return 172.938930000;
        case 174: return 173.940337000;
        case 0: case 175: return 174.940771000;
        case 176: return 175.942686000;
        case 177: return 176.943758000;
        case 178: return 177.945955000;
        case 179: return 178.947327000;
        case 180: return 179.949880000;
        case 181: return 180.951970000;
        case 182: return 181.955040000;
        case 183: return 182.957570000;
        case 184: return 183.960910000;
        }
        break;
      case 72:
        switch (isotope) {
        case 153: return 152.970690000;
        case 154: return 153.964860000;
        case 155: return 154.963390000;
        case 156: return 155.959360000;
        case 157: return 156.958400000;
        case 158: return 157.954799000;
        case 159: return 158.953995000;
        case 160: return 159.950684000;
        case 161: return 160.950275000;
        case 162: return 161.947210000;
        case 163: return 162.947090000;
        case 164: return 163.944367000;
        case 165: return 164.944570000;
        case 166: return 165.942180000;
        case 167: return 166.942600000;
        case 168: return 167.940570000;
        case 169: return 168.941260000;
        case 170: return 169.939610000;
        case 171: return 170.940490000;
        case 172: return 171.939448000;
        case 173: return 172.940510000;
        case 174: return 173.940046000;
        case 175: return 174.941509000;
        case 176: return 175.941408000;
        case 177: return 176.943220000;
        case 178: return 177.943698000;
        case 179: return 178.945816000;
        case 0: case 180: return 179.946550000;
        case 181: return 180.949101000;
        case 182: return 181.950554000;
        case 183: return 182.953530000;
        case 184: return 183.955450000;
        case 185: return 184.958820000;
        case 186: return 185.960890000;
        case 187: return 186.964590000;
        case 188: return 187.966850000;
        }
        break;
      case 73:
        switch (isotope) {
        case 155: return 154.974590000;
        case 156: return 155.972300000;
        case 157: return 156.968190000;
        case 158: return 157.966700000;
        case 159: return 158.963018000;
        case 160: return 159.961490000;
        case 161: return 160.958420000;
        case 162: return 161.957290000;
        case 163: return 162.954330000;
        case 164: return 163.953530000;
        case 165: return 164.950773000;
        case 166: return 165.950510000;
        case 167: return 166.948090000;
        case 168: return 167.948050000;
        case 169: return 168.946010000;
        case 170: return 169.946180000;
        case 171: return 170.944480000;
        case 172: return 171.944900000;
        case 173: return 172.943750000;
        case 174: return 173.944450000;
        case 175: return 174.943740000;
        case 176: return 175.944860000;
        case 177: return 176.944472000;
        case 178: return 177.945778000;
        case 179: return 178.945929000;
        case 180: return 179.947464000;
        case 0: case 181: return 180.947995000;
        case 182: return 181.950151000;
        case 183: return 182.951372000;
        case 184: return 183.954008000;
        case 185: return 184.955559000;
        case 186: return 185.958550000;
        case 187: return 186.960530000;
        case 188: return 187.963700000;
        case 189: return 188.965830000;
        case 190: return 189.969230000;
        }
        break;
      case 74:
        switch (isotope) {
        case 158: return 157.974560000;
        case 159: return 158.972920000;
        case 160: return 159.968480000;
        case 161: return 160.967360000;
        case 162: return 161.963497000;
        case 163: return 162.962520000;
        case 164: return 163.958954000;
        case 165: return 164.958280000;
        case 166: return 165.955027000;
        case 167: return 166.954816000;
        case 168: return 167.951808000;
        case 169: return 168.951779000;
        case 170: return 169.949228000;
        case 171: return 170.949450000;
        case 172: return 171.947290000;
        case 173: return 172.947690000;
        case 174: return 173.946080000;
        case 175: return 174.946720000;
        case 176: return 175.945630000;
        case 177: return 176.946640000;
        case 178: return 177.945876000;
        case 179: return 178.947070000;
        case 180: return 179.946704000;
        case 181: return 180.948197000;
        case 182: return 181.948204000;
        case 183: return 182.950223000;
        case 0: case 184: return 183.950931000;
        case 185: return 184.953419000;
        case 186: return 185.954364000;
        case 187: return 186.957160000;
        case 188: return 187.958489000;
        case 189: return 188.961910000;
        case 190: return 189.963180000;
        case 191: return 190.966600000;
        case 192: return 191.968170000;
        }
        break;
      case 75:
        switch (isotope) {
        case 160: return 159.982120000;
        case 161: return 160.977590000;
        case 162: return 161.976000000;
        case 163: return 162.972081000;
        case 164: return 163.970320000;
        case 165: return 164.967089000;
        case 166: return 165.965810000;
        case 167: return 166.962600000;
        case 168: return 167.961570000;
        case 169: return 168.958790000;
        case 170: return 169.958220000;
        case 171: return 170.955720000;
        case 172: return 171.955420000;
        case 173: return 172.953240000;
        case 174: return 173.953120000;
        case 175: return 174.951380000;
        case 176: return 175.951620000;
        case 177: return 176.950330000;
        case 178: return 177.950990000;
        case 179: return 178.949988000;
        case 180: return 179.950789000;
        case 181: return 180.950068000;
        case 182: return 181.951210000;
        case 183: return 182.950820000;
        case 184: return 183.952521000;
        case 185: return 184.952955000;
        case 186: return 185.954986000;
        case 0: case 187: return 186.955753000;
        case 188: return 187.958114000;
        case 189: return 188.959229000;
        case 190: return 189.961820000;
        case 191: return 190.963125000;
        case 192: return 191.965960000;
        case 193: return 192.967470000;
        case 194: return 193.970420000;
        }
        break;
      case 76:
        switch (isotope) {
        case 162: return 161.984430000;
        case 163: return 162.982690000;
        case 164: return 163.978040000;
        case 165: return 164.976760000;
        case 166: return 165.972691000;
        case 167: return 166.971550000;
        case 168: return 167.967804000;
        case 169: return 168.967019000;
        case 170: return 169.963577000;
        case 171: return 170.963185000;
        case 172: return 171.960023000;
        case 173: return 172.959808000;
        case 174: return 173.957062000;
        case 175: return 174.956946000;
        case 176: return 175.954810000;
        case 177: return 176.954965000;
        case 178: return 177.953251000;
        case 179: return 178.953816000;
        case 180: return 179.952379000;
        case 181: return 180.953240000;
        case 182: return 181.952110000;
        case 183: return 182.953130000;
        case 184: return 183.952489000;
        case 185: return 184.954042000;
        case 186: return 185.953838000;
        case 187: return 186.955750000;
        case 188: return 187.955838000;
        case 189: return 188.958147000;
        case 190: return 189.958447000;
        case 191: return 190.960929000;
        case 0: case 192: return 191.961480000;
        case 193: return 192.964151000;
        case 194: return 193.965182000;
        case 195: return 194.968130000;
        case 196: return 195.969640000;
        }
        break;
      case 77:
        switch (isotope) {
        case 164: return 163.992200000;
        case 165: return 164.987520000;
        case 166: return 165.985820000;
        case 167: return 166.981665000;
        case 168: return 167.979880000;
        case 169: return 168.976295000;
        case 170: return 169.974970000;
        case 171: return 170.971630000;
        case 172: return 171.970460000;
        case 173: return 172.967502000;
        case 174: return 173.966861000;
        case 175: return 174.964113000;
        case 176: return 175.963649000;
        case 177: return 176.961302000;
        case 178: return 177.961082000;
        case 179: return 178.959122000;
        case 180: return 179.959229000;
        case 181: return 180.957625000;
        case 182: return 181.958076000;
        case 183: return 182.956846000;
        case 184: return 183.957480000;
        case 185: return 184.956700000;
        case 186: return 185.957946000;
        case 187: return 186.957363000;
        case 188: return 187.958853000;
        case 189: return 188.958719000;
        case 190: return 189.960546000;
        case 191: return 190.960594000;
        case 192: return 191.962605000;
        case 0: case 193: return 192.962926000;
        case 194: return 193.965078000;
        case 195: return 194.965979000;
        case 196: return 195.968400000;
        case 197: return 196.969653000;
        case 198: return 197.972280000;
        case 199: return 198.973800000;
        }
        break;
      case 78:
        switch (isotope) {
        case 166: return 165.994860000;
        case 167: return 166.992980000;
        case 168: return 167.988150000;
        case 169: return 168.986720000;
        case 170: return 169.982495000;
        case 171: return 170.981240000;
        case 172: return 171.977347000;
        case 173: return 172.976440000;
        case 174: return 173.972819000;
        case 175: return 174.972421000;
        case 176: return 175.968945000;
        case 177: return 176.968469000;
        case 178: return 177.965649000;
        case 179: return 178.965363000;
        case 180: return 179.963031000;
        case 181: return 180.963097000;
        case 182: return 181.961171000;
        case 183: return 182.961597000;
        case 184: return 183.959922000;
        case 185: return 184.960620000;
        case 186: return 185.959351000;
        case 187: return 186.960590000;
        case 188: return 187.959395000;
        case 189: return 188.960834000;
        case 190: return 189.959932000;
        case 191: return 190.961677000;
        case 192: return 191.961038000;
        case 193: return 192.962987000;
        case 194: return 193.962680000;
        case 0: case 195: return 194.964791000;
        case 196: return 195.964951000;
        case 197: return 196.967340000;
        case 198: return 197.967893000;
        case 199: return 198.970593000;
        case 200: return 199.971441000;
        case 201: return 200.974510000;
        case 202: return 201.975740000;
        }
        break;
      case 79:
        switch (isotope) {
        case 169: return 168.998080000;
        case 170: return 169.996120000;
        case 171: return 170.991879000;
        case 172: return 171.990040000;
        case 173: return 172.986237000;
        case 174: return 173.984760000;
        case 175: return 174.981270000;
        case 176: return 175.980100000;
        case 177: return 176.976865000;
        case 178: return 177.976030000;
        case 179: return 178.973213000;
        case 180: return 179.972521000;
        case 181: return 180.970079000;
        case 182: return 181.969618000;
        case 183: return 182.967593000;
        case 184: return 183.967452000;
        case 185: return 184.965789000;
        case 186: return 185.965953000;
        case 187: return 186.964568000;
        case 188: return 187.965324000;
        case 189: return 188.963948000;
        case 190: return 189.964700000;
        case 191: return 190.963700000;
        case 192: return 191.964813000;
        case 193: return 192.964150000;
        case 194: return 193.965365000;
        case 195: return 194.965034000;
        case 196: return 195.966570000;
        case 0: case 197: return 196.966568000;
        case 198: return 197.968242000;
        case 199: return 198.968765000;
        case 200: return 199.970730000;
        case 201: return 200.971657000;
        case 202: return 201.973810000;
        case 203: return 202.975155000;
        case 204: return 203.977720000;
        case 205: return 204.979870000;
        }
        break;
      case 80:
        switch (isotope) {
        case 171: return 171.003760000;
        case 172: return 171.998830000;
        case 173: return 172.997240000;
        case 174: return 173.992864000;
        case 175: return 174.991420000;
        case 176: return 175.987355000;
        case 177: return 176.986280000;
        case 178: return 177.982483000;
        case 179: return 178.981834000;
        case 180: return 179.978266000;
        case 181: return 180.977819000;
        case 182: return 181.974690000;
        case 183: return 182.974450000;
        case 184: return 183.971713000;
        case 185: return 184.971899000;
        case 186: return 185.969362000;
        case 187: return 186.969814000;
        case 188: return 187.967577000;
        case 189: return 188.968190000;
        case 190: return 189.966322000;
        case 191: return 190.967157000;
        case 192: return 191.965634000;
        case 193: return 192.966665000;
        case 194: return 193.965439000;
        case 195: return 194.966720000;
        case 196: return 195.965833000;
        case 197: return 196.967213000;
        case 198: return 197.966769000;
        case 199: return 198.968279000;
        case 200: return 199.968326000;
        case 201: return 200.970302000;
        case 0: case 202: return 201.970643000;
        case 203: return 202.972872000;
        case 204: return 203.973493000;
        case 205: return 204.976073000;
        case 206: return 205.977514000;
        case 207: return 206.982590000;
        case 208: return 207.985940000;
        case 209: return 208.991040000;
        case 210: return 209.994510000;
        }
        break;
      case 81:
        switch (isotope) {
        case 176: return 176.000590000;
        case 177: return 176.996427000;
        case 178: return 177.994900000;
        case 179: return 178.991090000;
        case 180: return 179.989910000;
        case 181: return 180.986257000;
        case 182: return 181.985670000;
        case 183: return 182.982193000;
        case 184: return 183.981870000;
        case 185: return 184.978790000;
        case 186: return 185.978330000;
        case 187: return 186.975906000;
        case 188: return 187.976010000;
        case 189: return 188.973588000;
        case 190: return 189.973880000;
        case 191: return 190.971786000;
        case 192: return 191.972230000;
        case 193: return 192.970670000;
        case 194: return 193.971200000;
        case 195: return 194.969774000;
        case 196: return 195.970481000;
        case 197: return 196.969575000;
        case 198: return 197.970480000;
        case 199: return 198.969880000;
        case 200: return 199.970963000;
        case 201: return 200.970819000;
        case 202: return 201.972106000;
        case 203: return 202.972344000;
        case 204: return 203.973863000;
        case 0: case 205: return 204.974427000;
        case 206: return 205.976110000;
        case 207: return 206.977419000;
        case 208: return 207.982018000;
        case 209: return 208.985359000;
        case 210: return 209.990074000;
        case 211: return 210.993480000;
        case 212: return 211.998230000;
        }
        break;
      case 82:
        switch (isotope) {
        case 178: return 178.003830000;
        case 179: return 179.002150000;
        case 180: return 179.997918000;
        case 181: return 180.996620000;
        case 182: return 181.992672000;
        case 183: return 182.991870000;
        case 184: return 183.988142000;
        case 185: return 184.987610000;
        case 186: return 185.984239000;
        case 187: return 186.983918000;
        case 188: return 187.980874000;
        case 189: return 188.980810000;
        case 190: return 189.978082000;
        case 191: return 190.978270000;
        case 192: return 191.975785000;
        case 193: return 192.976170000;
        case 194: return 193.974012000;
        case 195: return 194.974542000;
        case 196: return 195.972774000;
        case 197: return 196.973431000;
        case 198: return 197.972034000;
        case 199: return 198.972917000;
        case 200: return 199.971827000;
        case 201: return 200.972885000;
        case 202: return 201.972159000;
        case 203: return 202.973391000;
        case 204: return 203.973043000;
        case 205: return 204.974481000;
        case 206: return 205.974465000;
        case 207: return 206.975896000;
        case 0: case 208: return 207.976652000;
        case 209: return 208.981090000;
        case 210: return 209.984188000;
        case 211: return 210.988737000;
        case 212: return 211.991897000;
        case 213: return 212.996581000;
        case 214: return 213.999805000;
        case 215: return 215.004810000;
        }
        break;
      case 83:
        switch (isotope) {
        case 184: return 184.001120000;
        case 185: return 184.997630000;
        case 186: return 185.996600000;
        case 187: return 186.993158000;
        case 188: return 187.992270000;
        case 189: return 188.989200000;
        case 190: return 189.988300000;
        case 191: return 190.985786000;
        case 192: return 191.985460000;
        case 193: return 192.982960000;
        case 194: return 193.982830000;
        case 195: return 194.980651000;
        case 196: return 195.980667000;
        case 197: return 196.978864000;
        case 198: return 197.979210000;
        case 199: return 198.977672000;
        case 200: return 199.978132000;
        case 201: return 200.977009000;
        case 202: return 201.977742000;
        case 203: return 202.976876000;
        case 204: return 203.977813000;
        case 205: return 204.977389000;
        case 206: return 205.978499000;
        case 207: return 206.978470000;
        case 208: return 207.979742000;
        case 0: case 209: return 208.980398000;
        case 210: return 209.984120000;
        case 211: return 210.987269000;
        case 212: return 211.991285000;
        case 213: return 212.994385000;
        case 214: return 213.998712000;
        case 215: return 215.001770000;
        case 216: return 216.006306000;
        case 217: return 217.009470000;
        case 218: return 218.014320000;
        }
        break;
      case 84:
        switch (isotope) {
        case 188: return 187.999422000;
        case 189: return 188.998481000;
        case 190: return 189.995101000;
        case 191: return 190.994574000;
        case 192: return 191.991335000;
        case 193: return 192.991030000;
        case 194: return 193.988186000;
        case 195: return 194.988110000;
        case 196: return 195.985535000;
        case 197: return 196.985660000;
        case 198: return 197.983389000;
        case 199: return 198.983666000;
        case 200: return 199.981799000;
        case 201: return 200.982260000;
        case 202: return 201.980758000;
        case 203: return 202.981420000;
        case 204: return 203.980318000;
        case 205: return 204.981203000;
        case 206: return 205.980481000;
        case 207: return 206.981593000;
        case 208: return 207.981245000;
        case 0: case 209: return 208.982430000;
        case 210: return 209.982873000;
        case 211: return 210.986653000;
        case 212: return 211.988868000;
        case 213: return 212.992857000;
        case 214: return 213.995201000;
        case 215: return 214.999420000;
        case 216: return 216.001915000;
        case 217: return 217.006335000;
        case 218: return 218.008973000;
        case 219: return 219.013740000;
        case 220: return 220.016600000;
        }
        break;
      case 85:
        switch (isotope) {
        case 193: return 192.999840000;
        case 194: return 193.998730000;
        case 195: return 194.996268000;
        case 196: return 195.995790000;
        case 197: return 196.993190000;
        case 198: return 197.992840000;
        case 199: return 198.990530000;
        case 200: return 199.990351000;
        case 201: return 200.988417000;
        case 202: return 201.988630000;
        case 203: return 202.986942000;
        case 204: return 203.987251000;
        case 205: return 204.986074000;
        case 206: return 205.986667000;
        case 207: return 206.985784000;
        case 208: return 207.986590000;
        case 209: return 208.986173000;
        case 0: case 210: return 209.987148000;
        case 211: return 210.987496000;
        case 212: return 211.990745000;
        case 213: return 212.992937000;
        case 214: return 213.996372000;
        case 215: return 214.998653000;
        case 216: return 216.002423000;
        case 217: return 217.004719000;
        case 218: return 218.008694000;
        case 219: return 219.011162000;
        case 220: return 220.015410000;
        case 221: return 221.018050000;
        case 222: return 222.022330000;
        case 223: return 223.025190000;
        }
        break;
      case 86:
        switch (isotope) {
        case 195: return 195.005440000;
        case 196: return 196.002115000;
        case 197: return 197.001580000;
        case 198: return 197.998679000;
        case 199: return 198.998370000;
        case 200: return 199.995699000;
        case 201: return 200.995630000;
        case 202: return 201.993263000;
        case 203: return 202.993387000;
        case 204: return 203.991429000;
        case 205: return 204.991720000;
        case 206: return 205.990214000;
        case 207: return 206.990734000;
        case 208: return 207.989642000;
        case 209: return 208.990415000;
        case 210: return 209.989696000;
        case 211: return 210.990601000;
        case 212: return 211.990704000;
        case 213: return 212.993883000;
        case 214: return 213.995363000;
        case 215: return 214.998745000;
        case 216: return 216.000274000;
        case 217: return 217.003928000;
        case 218: return 218.005601000;
        case 219: return 219.009480000;
        case 220: return 220.011394000;
        case 221: return 221.015537000;
        case 0: case 222: return 222.017577000;
        case 223: return 223.021790000;
        case 224: return 224.024090000;
        case 225: return 225.028440000;
        case 226: return 226.030890000;
        case 227: return 227.035410000;
        case 228: return 228.037990000;
        }
        break;
      case 87:
        switch (isotope) {
        case 199: return 199.007260000;
        case 200: return 200.006570000;
        case 201: return 201.003860000;
        case 202: return 202.003370000;
        case 203: return 203.000925000;
        case 204: return 204.000653000;
        case 205: return 204.998594000;
        case 206: return 205.998670000;
        case 207: return 206.996950000;
        case 208: return 207.997140000;
        case 209: return 208.995954000;
        case 210: return 209.996408000;
        case 211: return 210.995537000;
        case 212: return 211.996202000;
        case 213: return 212.996189000;
        case 214: return 213.998971000;
        case 215: return 215.000341000;
        case 216: return 216.003198000;
        case 217: return 217.004632000;
        case 218: return 218.007578000;
        case 219: return 219.009252000;
        case 220: return 220.012327000;
        case 221: return 221.014255000;
        case 222: return 222.017552000;
        case 0: case 223: return 223.019735000;
        case 224: return 224.023250000;
        case 225: return 225.025570000;
        case 226: return 226.029390000;
        case 227: return 227.031840000;
        case 228: return 228.035730000;
        case 229: return 229.038450000;
        case 230: return 230.042510000;
        case 231: return 231.045440000;
        case 232: return 232.049770000;
        }
        break;
      case 88:
        switch (isotope) {
        case 202: return 202.009890000;
        case 203: return 203.009270000;
        case 204: return 204.006500000;
        case 205: return 205.006270000;
        case 206: return 206.003827000;
        case 207: return 207.003800000;
        case 208: return 208.001840000;
        case 209: return 209.001990000;
        case 210: return 210.000495000;
        case 211: return 211.000898000;
        case 212: return 211.999794000;
        case 213: return 213.000384000;
        case 214: return 214.000108000;
        case 215: return 215.002720000;
        case 216: return 216.003533000;
        case 217: return 217.006320000;
        case 218: return 218.007140000;
        case 219: return 219.010085000;
        case 220: return 220.011028000;
        case 221: return 221.013917000;
        case 222: return 222.015375000;
        case 223: return 223.018502000;
        case 224: return 224.020211000;
        case 225: return 225.023612000;
        case 0: case 226: return 226.025409000;
        case 227: return 227.029177000;
        case 228: return 228.031070000;
        case 229: return 229.034958000;
        case 230: return 230.037056000;
        case 231: return 231.041220000;
        case 232: return 232.043640000;
        case 233: return 233.048060000;
        case 234: return 234.050700000;
        }
        break;
      case 89:
        switch (isotope) {
        case 206: return 206.014500000;
        case 207: return 207.011950000;
        case 208: return 208.011550000;
        case 209: return 209.009490000;
        case 210: return 210.009440000;
        case 211: return 211.007730000;
        case 212: return 212.007810000;
        case 213: return 213.006610000;
        case 214: return 214.006902000;
        case 215: return 215.006454000;
        case 216: return 216.008720000;
        case 217: return 217.009347000;
        case 218: return 218.011640000;
        case 219: return 219.012420000;
        case 220: return 220.014763000;
        case 221: return 221.015590000;
        case 222: return 222.017844000;
        case 223: return 223.019137000;
        case 224: return 224.021723000;
        case 225: return 225.023230000;
        case 226: return 226.026098000;
        case 0: case 227: return 227.027752000;
        case 228: return 228.031021000;
        case 229: return 229.033020000;
        case 230: return 230.036290000;
        case 231: return 231.038560000;
        case 232: return 232.042030000;
        case 233: return 233.044550000;
        case 234: return 234.048420000;
        case 235: return 235.051230000;
        case 236: return 236.055300000;
        }
        break;
      case 90:
        switch (isotope) {
        case 209: return 209.017720000;
        case 210: return 210.015075000;
        case 211: return 211.014930000;
        case 212: return 212.012980000;
        case 213: return 213.013010000;
        case 214: return 214.011500000;
        case 215: return 215.011730000;
        case 216: return 216.011062000;
        case 217: return 217.013114000;
        case 218: return 218.013284000;
        case 219: return 219.015540000;
        case 220: return 220.015748000;
        case 221: return 221.018184000;
        case 222: return 222.018468000;
        case 223: return 223.020811000;
        case 224: return 224.021467000;
        case 225: return 225.023951000;
        case 226: return 226.024903000;
        case 227: return 227.027704000;
        case 228: return 228.028741000;
        case 229: return 229.031762000;
        case 230: return 230.033133000;
        case 231: return 231.036304000;
        case 0: case 232: return 232.038055000;
        case 233: return 233.041581000;
        case 234: return 234.043601000;
        case 235: return 235.047510000;
        case 236: return 236.049870000;
        case 237: return 237.053890000;
        case 238: return 238.056500000;
        }
        break;
      case 91:
        switch (isotope) {
        case 212: return 212.023200000;
        case 213: return 213.021110000;
        case 214: return 214.020920000;
        case 215: return 215.019190000;
        case 216: return 216.019110000;
        case 217: return 217.018320000;
        case 218: return 218.020042000;
        case 219: return 219.019880000;
        case 220: return 220.021880000;
        case 221: return 221.021880000;
        case 222: return 222.023740000;
        case 223: return 223.023960000;
        case 224: return 224.025626000;
        case 225: return 225.026130000;
        case 226: return 226.027948000;
        case 227: return 227.028805000;
        case 228: return 228.031051000;
        case 229: return 229.032096000;
        case 230: return 230.034541000;
        case 0: case 231: return 231.035884000;
        case 232: return 232.038592000;
        case 233: return 233.040247000;
        case 234: return 234.043308000;
        case 235: return 235.045440000;
        case 236: return 236.048680000;
        case 237: return 237.051150000;
        case 238: return 238.054500000;
        case 239: return 239.057260000;
        case 240: return 240.060980000;
        }
        break;
      case 92:
        switch (isotope) {
        case 217: return 217.024370000;
        case 218: return 218.023540000;
        case 219: return 219.024920000;
        case 220: return 220.024720000;
        case 221: return 221.026400000;
        case 222: return 222.026090000;
        case 223: return 223.027740000;
        case 224: return 224.027605000;
        case 225: return 225.029391000;
        case 226: return 226.029339000;
        case 227: return 227.031156000;
        case 228: return 228.031374000;
        case 229: return 229.033506000;
        case 230: return 230.033940000;
        case 231: return 231.036294000;
        case 232: return 232.037156000;
        case 233: return 233.039635000;
        case 234: return 234.040952000;
        case 235: return 235.043929000;
        case 236: return 236.045568000;
        case 237: return 237.048730000;
        case 0: case 238: return 238.050788000;
        case 239: return 239.054293000;
        case 240: return 240.056592000;
        case 241: return 241.060330000;
        case 242: return 242.062930000;
        }
        break;
      case 93:
        switch (isotope) {
        case 225: return 225.033910000;
        case 226: return 226.035150000;
        case 227: return 227.034960000;
        case 228: return 228.036180000;
        case 229: return 229.036260000;
        case 230: return 230.037830000;
        case 231: return 231.038250000;
        case 232: return 232.040110000;
        case 233: return 233.040740000;
        case 234: return 234.042895000;
        case 235: return 235.044063000;
        case 236: return 236.046570000;
        case 0: case 237: return 237.048173000;
        case 238: return 238.050946000;
        case 239: return 239.052939000;
        case 240: return 240.056162000;
        case 241: return 241.058250000;
        case 242: return 242.061640000;
        case 243: return 243.064280000;
        case 244: return 244.067850000;
        }
        break;
      case 94:
        switch (isotope) {
        case 228: return 228.038740000;
        case 229: return 229.040150000;
        case 230: return 230.039650000;
        case 231: return 231.041101000;
        case 232: return 232.041187000;
        case 233: return 233.043000000;
        case 234: return 234.043317000;
        case 235: return 235.045286000;
        case 236: return 236.046058000;
        case 237: return 237.048409000;
        case 238: return 238.049559000;
        case 239: return 239.052163000;
        case 240: return 240.053813000;
        case 241: return 241.056851000;
        case 242: return 242.058742000;
        case 243: return 243.062003000;
        case 0: case 244: return 244.064204000;
        case 245: return 245.067747000;
        case 246: return 246.070205000;
        case 247: return 247.074070000;
        }
        break;
      case 95:
        switch (isotope) {
        case 231: return 231.045560000;
        case 232: return 232.046590000;
        case 233: return 233.046350000;
        case 234: return 234.047810000;
        case 235: return 235.047950000;
        case 236: return 236.049580000;
        case 237: return 237.050000000;
        case 238: return 238.051980000;
        case 239: return 239.053024000;
        case 240: return 240.055300000;
        case 241: return 241.056829000;
        case 242: return 242.059549000;
        case 0: case 243: return 243.061381000;
        case 244: return 244.064284000;
        case 245: return 245.066452000;
        case 246: return 246.069775000;
        case 247: return 247.072090000;
        case 248: return 248.075750000;
        case 249: return 249.078480000;
        }
        break;
      case 96:
        switch (isotope) {
        case 233: return 233.050770000;
        case 234: return 234.050160000;
        case 235: return 235.051430000;
        case 236: return 236.051410000;
        case 237: return 237.052900000;
        case 238: return 238.053030000;
        case 239: return 239.054960000;
        case 240: return 240.055529000;
        case 241: return 241.057653000;
        case 242: return 242.058835000;
        case 243: return 243.061389000;
        case 244: return 244.062752000;
        case 245: return 245.065491000;
        case 246: return 246.067223000;
        case 0: case 247: return 247.070354000;
        case 248: return 248.072349000;
        case 249: return 249.075953000;
        case 250: return 250.078357000;
        case 251: return 251.082285000;
        case 252: return 252.084870000;
        }
        break;
      case 97:
        switch (isotope) {
        case 235: return 235.056580000;
        case 236: return 236.057330000;
        case 237: return 237.057000000;
        case 238: return 238.058280000;
        case 239: return 239.058280000;
        case 240: return 240.059760000;
        case 241: return 241.060230000;
        case 242: return 242.061980000;
        case 243: return 243.063008000;
        case 244: return 244.065181000;
        case 245: return 245.066361000;
        case 246: return 246.068670000;
        case 0: case 247: return 247.070307000;
        case 248: return 248.073090000;
        case 249: return 249.074986000;
        case 250: return 250.078317000;
        case 251: return 251.080760000;
        case 252: return 252.084310000;
        case 253: return 253.086880000;
        case 254: return 254.090600000;
        }
        break;
      case 98:
        switch (isotope) {
        case 237: return 237.062070000;
        case 238: return 238.061410000;
        case 239: return 239.062420000;
        case 240: return 240.062300000;
        case 241: return 241.063730000;
        case 242: return 242.063700000;
        case 243: return 243.065430000;
        case 244: return 244.066001000;
        case 245: return 245.068049000;
        case 246: return 246.068805000;
        case 247: return 247.071001000;
        case 248: return 248.072185000;
        case 249: return 249.074853000;
        case 250: return 250.076406000;
        case 0: case 251: return 251.079587000;
        case 252: return 252.081626000;
        case 253: return 253.085133000;
        case 254: return 254.087323000;
        case 255: return 255.091050000;
        case 256: return 256.093440000;
        }
        break;
      case 99:
        switch (isotope) {
        case 240: return 240.068920000;
        case 241: return 241.068540000;
        case 242: return 242.069750000;
        case 243: return 243.069550000;
        case 244: return 244.070880000;
        case 245: return 245.071320000;
        case 246: return 246.072900000;
        case 247: return 247.073660000;
        case 248: return 248.075470000;
        case 249: return 249.076410000;
        case 250: return 250.078610000;
        case 251: return 251.079992000;
        case 0: case 252: return 252.082980000;
        case 253: return 253.084824000;
        case 254: return 254.088022000;
        case 255: return 255.090273000;
        case 256: return 256.093600000;
        case 257: return 257.095980000;
        case 258: return 258.099520000;
        }
        break;
      case 100:
        switch (isotope) {
        case 242: return 242.073430000;
        case 243: return 243.074350000;
        case 244: return 244.074080000;
        case 245: return 245.075390000;
        case 246: return 246.075300000;
        case 247: return 247.076850000;
        case 248: return 248.077195000;
        case 249: return 249.079030000;
        case 250: return 250.079521000;
        case 251: return 251.081575000;
        case 252: return 252.082467000;
        case 253: return 253.085185000;
        case 254: return 254.086854000;
        case 255: return 255.089962000;
        case 256: return 256.091773000;
        case 0: case 257: return 257.095105000;
        case 258: return 258.097080000;
        case 259: return 259.100600000;
        case 260: return 260.102680000;
        }
        break;
      case 101:
        switch (isotope) {
        case 245: return 245.080830000;
        case 246: return 246.081890000;
        case 247: return 247.081640000;
        case 248: return 248.082820000;
        case 249: return 249.083010000;
        case 250: return 250.084420000;
        case 251: return 251.084840000;
        case 252: return 252.086560000;
        case 253: return 253.087280000;
        case 254: return 254.089660000;
        case 255: return 255.091083000;
        case 256: return 256.094060000;
        case 257: return 257.095541000;
        case 0: case 258: return 258.098431000;
        case 259: return 259.100510000;
        case 260: return 260.103650000;
        case 261: return 261.105720000;
        case 262: return 262.108870000;
        }
        break;
      case 102:
        switch (isotope) {
        case 248: return 248.086600000;
        case 249: return 249.087830000;
        case 250: return 250.087510000;
        case 251: return 251.089010000;
        case 252: return 252.088977000;
        case 253: return 253.090680000;
        case 254: return 254.090955000;
        case 255: return 255.093241000;
        case 256: return 256.094283000;
        case 257: return 257.096877000;
        case 258: return 258.098210000;
        case 0: case 259: return 259.101030000;
        case 260: return 260.102640000;
        case 261: return 261.105750000;
        case 262: return 262.107300000;
        case 263: return 263.110550000;
        case 264: return 264.112350000;
        }
        break;
      case 103:
        switch (isotope) {
        case 251: return 251.094360000;
        case 252: return 252.095370000;
        case 253: return 253.095210000;
        case 254: return 254.096450000;
        case 255: return 255.096680000;
        case 256: return 256.098630000;
        case 257: return 257.099560000;
        case 258: return 258.101810000;
        case 259: return 259.102900000;
        case 260: return 260.105500000;
        case 261: return 261.106880000;
        case 0: case 262: return 262.109630000;
        case 263: return 263.111290000;
        case 264: return 264.114040000;
        case 265: return 265.115840000;
        case 266: return 266.119310000;
        }
        break;
      case 104:
        switch (isotope) {
        case 253: return 253.100690000;
        case 254: return 254.100180000;
        case 255: return 255.101340000;
        case 256: return 256.101166000;
        case 257: return 257.102990000;
        case 258: return 258.103490000;
        case 259: return 259.105640000;
        case 260: return 260.106440000;
        case 0: case 261: return 261.108770000;
        case 262: return 262.109930000;
        case 263: return 263.112550000;
        case 264: return 264.113990000;
        case 265: return 265.116700000;
        case 266: return 266.117960000;
        case 267: return 267.121530000;
        case 268: return 268.123640000;
        }
        break;
      case 105:
        switch (isotope) {
        case 255: return 255.107400000;
        case 256: return 256.108130000;
        case 257: return 257.107720000;
        case 258: return 258.109230000;
        case 259: return 259.109610000;
        case 260: return 260.111300000;
        case 261: return 261.112060000;
        case 0: case 262: return 262.114080000;
        case 263: return 263.114990000;
        case 264: return 264.117400000;
        case 265: return 265.118600000;
        case 266: return 266.121030000;
        case 267: return 267.122380000;
        case 268: return 268.125450000;
        case 269: return 269.127460000;
        case 270: return 270.130710000;
        }
        break;
      case 106:
        switch (isotope) {
        case 258: return 258.113170000;
        case 259: return 259.114500000;
        case 260: return 260.114420000;
        case 261: return 261.116120000;
        case 262: return 262.116400000;
        case 0: case 263: return 263.118320000;
        case 264: return 264.118930000;
        case 265: return 265.121110000;
        case 266: return 266.122070000;
        case 267: return 267.124430000;
        case 268: return 268.125610000;
        case 269: return 269.128760000;
        case 270: return 270.130330000;
        case 271: return 271.133470000;
        case 272: return 272.135160000;
        case 273: return 273.138220000;
        }
        break;
      case 107:
        switch (isotope) {
        case 260: return 260.121970000;
        case 261: return 261.121660000;
        case 262: return 262.122890000;
        case 263: return 263.123040000;
        case 0: case 264: return 264.124600000;
        case 265: return 265.125150000;
        case 266: return 266.126940000;
        case 267: return 267.127650000;
        case 268: return 268.129760000;
        case 269: return 269.130690000;
        case 270: return 270.133620000;
        case 271: return 271.135180000;
        case 272: return 272.138030000;
        case 273: return 273.139620000;
        case 274: return 274.142440000;
        case 275: return 275.144250000;
        }
        break;
      case 108:
        switch (isotope) {
        case 263: return 263.128560000;
        case 264: return 264.128390000;
        case 0: case 265: return 265.130090000;
        case 266: return 266.130100000;
        case 267: return 267.131790000;
        case 268: return 268.132160000;
        case 269: return 269.134060000;
        case 270: return 270.134650000;
        case 271: return 271.137660000;
        case 272: return 272.139050000;
        case 273: return 273.141990000;
        case 274: return 274.143130000;
        case 275: return 275.145950000;
        case 276: return 276.147210000;
        case 277: return 277.149840000;
        }
        break;
      case 109:
        switch (isotope) {
        case 265: return 265.136150000;
        case 266: return 266.137300000;
        case 267: return 267.137310000;
        case 0: case 268: return 268.138730000;
        case 269: return 269.139060000;
        case 270: return 270.140660000;
        case 271: return 271.141140000;
        case 272: return 272.143740000;
        case 273: return 273.144910000;
        case 274: return 274.147490000;
        case 275: return 275.148650000;
        case 276: return 276.151160000;
        case 277: return 277.152420000;
        case 278: return 278.154810000;
        case 279: return 279.156190000;
        }
        break;
      case 110:
        switch (isotope) {
        case 0: case 281: return 281.162061;
        }
        break;
      case 111:
        switch (isotope) {
        case 0: case 280: return 280.164473;
        }
        break;
      case 112:
        switch (isotope) {
        case 0: case 285: return 285.174105;
        }
        break;
      case 113:
        switch (isotope) {
        case 0: case 284: return 284.178080;
        }
        break;
      case 114:
        switch (isotope) {
        case 0: case 289: return 289.187279;
        }
        break;
      case 115:
        switch (isotope) {
        case 0: case 288: return 288.192492;
        }
        break;
      case 116:
        switch (isotope) {
        case 0: case 292: return 292.199786;
        }
        break;
      case 117:
        switch (isotope) {
        case 291: return 291.206560;
        case 0: case 292: return 292.207550;
        }
        break;
      case 118:
        switch (isotope) {
        case 0: case 293: return 293.214670;
        }
        break;

      }
      return 0.0;
    }

  } // namespace OBElements
} // namespace OpenBabel

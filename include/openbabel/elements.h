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
#ifndef OB_ELEMENTS_H
#define OB_ELEMENTS_H

#include <openbabel/babelconfig.h>

namespace OpenBabel
{
  /** \brief Functions and constants for handling the elements of the periodic table and associated data

  Translating element data is a common task given that many file
  formats give either element symbol or atomic number information, but
  not both. The OBElements namespace contains functions to handle this conversion,
  as well as providing information associated with particular elements.
  The following code sample demonstrates the use of members of this namespace:
  \code
  cout << "The symbol for element 6 is " << OBElements::GetSymbol(6) << endl;
  cout << "The atomic number for Sulfur is " << OBElements::GetAtomicNum(16) << endl;
  cout << "The van der Waal radius for Nitrogen is " << OBElements::GetVdwRad(7);
  if (atom->GetAtomicNum() == OBElements::Carbon) {
    // then do something
  }
  \endcode

  Stored information in the OBElementTable includes elemental:
  - symbols
  - covalent radii
  - van der Waal radii
  - expected maximum bonding valence
  - molar mass (by IUPAC recommended atomic masses)
  - isotopic mass
  - electronegativity (Pauling and Allred-Rochow)
  - ionization potential
  - electron affinity
  - RGB colors for visualization programs
  - names (by IUPAC recommendation)
  **/
  namespace OBElements {

    //! \return the element symbol matching the atomic number passed
    OBAPI const char* GetSymbol(unsigned int atomic_number);
    //! \return the name of this element
    OBAPI const char* GetName(unsigned int atomic_number);
    //! \return the average atomic mass for this element.
    //! For exact isotope masses, use GetExactMass()
    OBAPI double GetMass(unsigned int atomic_number);
    //! \return the exact mass of the specified isotope
    //!         or by default (i.e. "isotope 0") the mass of the most abundant isotope
    OBAPI double GetExactMass(unsigned int atomic_number, unsigned int isotope=0);
    //  //! \return the atomic number matching the element symbol
    OBAPI unsigned int GetAtomicNum(const char* ptr);
    //! \return the Allred-Rochow electronegativity for this element
    OBAPI double GetAllredRochowElectroNeg(unsigned int atomic_number);
    //! \return the covalent radius (in Angstrom) for this atomic number
    OBAPI double GetCovalentRad(unsigned int atomic_number);
    //! \return the van der Waals radius (in Angstrom) for this atomic number
    OBAPI double GetVdwRad(unsigned int atomic_number);
    //! \return the electron affinity (in eV) for this element
    OBAPI double GetElectronAffinity(unsigned int atomic_number);
    //! \return the ionization potential (in eV) for this element
    OBAPI double GetIonization(unsigned int atomic_number);
    //! \return the maximum expected number of bonds to this element
    OBAPI unsigned int GetMaxBonds(unsigned int atomic_number);
    //! \return the Pauling electronegativity for this element
    OBAPI double GetElectroNeg(unsigned int atomic_number);
    //! Sets the red, green, and blue color values for this element
    OBAPI void GetRGB(unsigned int atomic_number, double *r, double *g, double *b);
    //! The atomic numbers of the elements
    const unsigned int Hydrogen = 1; // H
    const unsigned int Helium = 2; // He
    const unsigned int Lithium = 3; // Li
    const unsigned int Beryllium = 4; // Be
    const unsigned int Boron = 5; // B
    const unsigned int Carbon = 6; // C
    const unsigned int Nitrogen = 7; // N
    const unsigned int Oxygen = 8; // O
    const unsigned int Fluorine = 9; // F
    const unsigned int Neon = 10; // Ne
    const unsigned int Sodium = 11; // Na
    const unsigned int Magnesium = 12; // Mg
    const unsigned int Aluminium = 13; // Al
    const unsigned int Silicon = 14; // Si
    const unsigned int Phosphorus = 15; // P
    const unsigned int Sulfur = 16; // S
    const unsigned int Chlorine = 17; // Cl
    const unsigned int Argon = 18; // Ar
    const unsigned int Potassium = 19; // K
    const unsigned int Calcium = 20; // Ca
    const unsigned int Scandium = 21; // Sc
    const unsigned int Titanium = 22; // Ti
    const unsigned int Vanadium = 23; // V
    const unsigned int Chromium = 24; // Cr
    const unsigned int Manganese = 25; // Mn
    const unsigned int Iron = 26; // Fe
    const unsigned int Cobalt = 27; // Co
    const unsigned int Nickel = 28; // Ni
    const unsigned int Copper = 29; // Cu
    const unsigned int Zinc = 30; // Zn
    const unsigned int Gallium = 31; // Ga
    const unsigned int Germanium = 32; // Ge
    const unsigned int Arsenic = 33; // As
    const unsigned int Selenium = 34; // Se
    const unsigned int Bromine = 35; // Br
    const unsigned int Krypton = 36; // Kr
    const unsigned int Rubidium = 37; // Rb
    const unsigned int Strontium = 38; // Sr
    const unsigned int Yttrium = 39; // Y
    const unsigned int Zirconium = 40; // Zr
    const unsigned int Niobium = 41; // Nb
    const unsigned int Molybdenum = 42; // Mo
    const unsigned int Technetium = 43; // Tc
    const unsigned int Ruthenium = 44; // Ru
    const unsigned int Rhodium = 45; // Rh
    const unsigned int Palladium = 46; // Pd
    const unsigned int Silver = 47; // Ag
    const unsigned int Cadmium = 48; // Cd
    const unsigned int Indium = 49; // In
    const unsigned int Tin = 50; // Sn
    const unsigned int Antimony = 51; // Sb
    const unsigned int Tellurium = 52; // Te
    const unsigned int Iodine = 53; // I
    const unsigned int Xenon = 54; // Xe
    const unsigned int Caesium = 55; // Cs
    const unsigned int Barium = 56; // Ba
    const unsigned int Lanthanum = 57; // La
    const unsigned int Cerium = 58; // Ce
    const unsigned int Praseodymium = 59; // Pr
    const unsigned int Neodymium = 60; // Nd
    const unsigned int Promethium = 61; // Pm
    const unsigned int Samarium = 62; // Sm
    const unsigned int Europium = 63; // Eu
    const unsigned int Gadolinium = 64; // Gd
    const unsigned int Terbium = 65; // Tb
    const unsigned int Dysprosium = 66; // Dy
    const unsigned int Holmium = 67; // Ho
    const unsigned int Erbium = 68; // Er
    const unsigned int Thulium = 69; // Tm
    const unsigned int Ytterbium = 70; // Yb
    const unsigned int Lutetium = 71; // Lu
    const unsigned int Hafnium = 72; // Hf
    const unsigned int Tantalum = 73; // Ta
    const unsigned int Tungsten = 74; // W
    const unsigned int Rhenium = 75; // Re
    const unsigned int Osmium = 76; // Os
    const unsigned int Iridium = 77; // Ir
    const unsigned int Platinum = 78; // Pt
    const unsigned int Gold = 79; // Au
    const unsigned int Mercury = 80; // Hg
    const unsigned int Thallium = 81; // Tl
    const unsigned int Lead = 82; // Pb
    const unsigned int Bismuth = 83; // Bi
    const unsigned int Polonium = 84; // Po
    const unsigned int Astatine = 85; // At
    const unsigned int Radon = 86; // Rn
    const unsigned int Francium = 87; // Fr
    const unsigned int Radium = 88; // Ra
    const unsigned int Actinium = 89; // Ac
    const unsigned int Thorium = 90; // Th
    const unsigned int Protactinium = 91; // Pa
    const unsigned int Uranium = 92; // U
    const unsigned int Neptunium = 93; // Np
    const unsigned int Plutonium = 94; // Pu
    const unsigned int Americium = 95; // Am
    const unsigned int Curium = 96; // Cm
    const unsigned int Berkelium = 97; // Bk
    const unsigned int Californium = 98; // Cf
    const unsigned int Einsteinium = 99; // Es
    const unsigned int Fermium = 100; // Fm
    const unsigned int Mendelevium = 101; // Md
    const unsigned int Nobelium = 102; // No
    const unsigned int Lawrencium = 103; // Lr
    const unsigned int Rutherfordium = 104; // Rf
    const unsigned int Dubnium = 105; // Db
    const unsigned int Seaborgium = 106; // Sg
    const unsigned int Bohrium = 107; // Bh
    const unsigned int Hassium = 108; // Hs
    const unsigned int Meitnerium = 109; // Mt
    const unsigned int Darmstadtium = 110; // Ds
    const unsigned int Roentgenium = 111; // Rg
    const unsigned int Copernicium = 112; // Cn
    const unsigned int Nihonium = 113; // Nh
    const unsigned int Flerovium = 114; // Fl
    const unsigned int Moscovium = 115; // Mc
    const unsigned int Livermorium = 116; // Lv
    const unsigned int Tennessine = 117; // Ts
    const unsigned int Oganesson = 118; // Og

  }
}

#endif //OB_ELEMENTS_H

//! \file elements.h
//! \brief Functions relating to elements

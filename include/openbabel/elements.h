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
    const unsigned int Dummy = 0;
    const unsigned int Hydrogen = 1;
    const unsigned int H = 1;
    const unsigned int Helium = 2;
    const unsigned int He = 2;
    const unsigned int Lithium = 3;
    const unsigned int Li = 3;
    const unsigned int Beryllium = 4;
    const unsigned int Be = 4;
    const unsigned int Boron = 5;
    const unsigned int B = 5;
    const unsigned int Carbon = 6;
    const unsigned int C = 6;
    const unsigned int Nitrogen = 7;
    const unsigned int N = 7;
    const unsigned int Oxygen = 8;
    const unsigned int O = 8;
    const unsigned int Fluorine = 9;
    const unsigned int F = 9;
    const unsigned int Neon = 10;
    const unsigned int Ne = 10;
    const unsigned int Sodium = 11;
    const unsigned int Na = 11;
    const unsigned int Magnesium = 12;
    const unsigned int Mg = 12;
    const unsigned int Aluminium = 13;
    const unsigned int Al = 13;
    const unsigned int Silicon = 14;
    const unsigned int Si = 14;
    const unsigned int Phosphorus = 15;
    const unsigned int P = 15;
    const unsigned int Sulfur = 16;
    const unsigned int S = 16;
    const unsigned int Chlorine = 17;
    const unsigned int Cl = 17;
    const unsigned int Argon = 18;
    const unsigned int Ar = 18;
    const unsigned int Potassium = 19;
    const unsigned int K = 19;
    const unsigned int Calcium = 20;
    const unsigned int Ca = 20;
    const unsigned int Scandium = 21;
    const unsigned int Sc = 21;
    const unsigned int Titanium = 22;
    const unsigned int Ti = 22;
    const unsigned int Vanadium = 23;
    const unsigned int V = 23;
    const unsigned int Chromium = 24;
    const unsigned int Cr = 24;
    const unsigned int Manganese = 25;
    const unsigned int Mn = 25;
    const unsigned int Iron = 26;
    const unsigned int Fe = 26;
    const unsigned int Cobalt = 27;
    const unsigned int Co = 27;
    const unsigned int Nickel = 28;
    const unsigned int Ni = 28;
    const unsigned int Copper = 29;
    const unsigned int Cu = 29;
    const unsigned int Zinc = 30;
    const unsigned int Zn = 30;
    const unsigned int Gallium = 31;
    const unsigned int Ga = 31;
    const unsigned int Germanium = 32;
    const unsigned int Ge = 32;
    const unsigned int Arsenic = 33;
    const unsigned int As = 33;
    const unsigned int Selenium = 34;
    const unsigned int Se = 34;
    const unsigned int Bromine = 35;
    const unsigned int Br = 35;
    const unsigned int Krypton = 36;
    const unsigned int Kr = 36;
    const unsigned int Rubidium = 37;
    const unsigned int Rb = 37;
    const unsigned int Strontium = 38;
    const unsigned int Sr = 38;
    const unsigned int Yttrium = 39;
    const unsigned int Y = 39;
    const unsigned int Zirconium = 40;
    const unsigned int Zr = 40;
    const unsigned int Niobium = 41;
    const unsigned int Nb = 41;
    const unsigned int Molybdenum = 42;
    const unsigned int Mo = 42;
    const unsigned int Technetium = 43;
    const unsigned int Tc = 43;
    const unsigned int Ruthenium = 44;
    const unsigned int Ru = 44;
    const unsigned int Rhodium = 45;
    const unsigned int Rh = 45;
    const unsigned int Palladium = 46;
    const unsigned int Pd = 46;
    const unsigned int Silver = 47;
    const unsigned int Ag = 47;
    const unsigned int Cadmium = 48;
    const unsigned int Cd = 48;
    const unsigned int Indium = 49;
    const unsigned int In = 49;
    const unsigned int Tin = 50;
    const unsigned int Sn = 50;
    const unsigned int Antimony = 51;
    const unsigned int Sb = 51;
    const unsigned int Tellurium = 52;
    const unsigned int Te = 52;
    const unsigned int Iodine = 53;
    const unsigned int I = 53;
    const unsigned int Xenon = 54;
    const unsigned int Xe = 54;
    const unsigned int Caesium = 55;
    const unsigned int Cs = 55;
    const unsigned int Barium = 56;
    const unsigned int Ba = 56;
    const unsigned int Lanthanum = 57;
    const unsigned int La = 57;
    const unsigned int Cerium = 58;
    const unsigned int Ce = 58;
    const unsigned int Praseodymium = 59;
    const unsigned int Pr = 59;
    const unsigned int Neodymium = 60;
    const unsigned int Nd = 60;
    const unsigned int Promethium = 61;
    const unsigned int Pm = 61;
    const unsigned int Samarium = 62;
    const unsigned int Sm = 62;
    const unsigned int Europium = 63;
    const unsigned int Eu = 63;
    const unsigned int Gadolinium = 64;
    const unsigned int Gd = 64;
    const unsigned int Terbium = 65;
    const unsigned int Tb = 65;
    const unsigned int Dysprosium = 66;
    const unsigned int Dy = 66;
    const unsigned int Holmium = 67;
    const unsigned int Ho = 67;
    const unsigned int Erbium = 68;
    const unsigned int Er = 68;
    const unsigned int Thulium = 69;
    const unsigned int Tm = 69;
    const unsigned int Ytterbium = 70;
    const unsigned int Yb = 70;
    const unsigned int Lutetium = 71;
    const unsigned int Lu = 71;
    const unsigned int Hafnium = 72;
    const unsigned int Hf = 72;
    const unsigned int Tantalum = 73;
    const unsigned int Ta = 73;
    const unsigned int Tungsten = 74;
    const unsigned int W = 74;
    const unsigned int Rhenium = 75;
    const unsigned int Re = 75;
    const unsigned int Osmium = 76;
    const unsigned int Os = 76;
    const unsigned int Iridium = 77;
    const unsigned int Ir = 77;
    const unsigned int Platinum = 78;
    const unsigned int Pt = 78;
    const unsigned int Gold = 79;
    const unsigned int Au = 79;
    const unsigned int Mercury = 80;
    const unsigned int Hg = 80;
    const unsigned int Thallium = 81;
    const unsigned int Tl = 81;
    const unsigned int Lead = 82;
    const unsigned int Pb = 82;
    const unsigned int Bismuth = 83;
    const unsigned int Bi = 83;
    const unsigned int Polonium = 84;
    const unsigned int Po = 84;
    const unsigned int Astatine = 85;
    const unsigned int At = 85;
    const unsigned int Radon = 86;
    const unsigned int Rn = 86;
    const unsigned int Francium = 87;
    const unsigned int Fr = 87;
    const unsigned int Radium = 88;
    const unsigned int Ra = 88;
    const unsigned int Actinium = 89;
    const unsigned int Ac = 89;
    const unsigned int Thorium = 90;
    const unsigned int Th = 90;
    const unsigned int Protactinium = 91;
    const unsigned int Pa = 91;
    const unsigned int Uranium = 92;
    const unsigned int U = 92;
    const unsigned int Neptunium = 93;
    const unsigned int Np = 93;
    const unsigned int Plutonium = 94;
    const unsigned int Pu = 94;
    const unsigned int Americium = 95;
    const unsigned int Am = 95;
    const unsigned int Curium = 96;
    const unsigned int Cm = 96;
    const unsigned int Berkelium = 97;
    const unsigned int Bk = 97;
    const unsigned int Californium = 98;
    const unsigned int Cf = 98;
    const unsigned int Einsteinium = 99;
    const unsigned int Es = 99;
    const unsigned int Fermium = 100;
    const unsigned int Fm = 100;
    const unsigned int Mendelevium = 101;
    const unsigned int Md = 101;
    const unsigned int Nobelium = 102;
    const unsigned int No = 102;
    const unsigned int Lawrencium = 103;
    const unsigned int Lr = 103;
    const unsigned int Rutherfordium = 104;
    const unsigned int Rf = 104;
    const unsigned int Dubnium = 105;
    const unsigned int Db = 105;
    const unsigned int Seaborgium = 106;
    const unsigned int Sg = 106;
    const unsigned int Bohrium = 107;
    const unsigned int Bh = 107;
    const unsigned int Hassium = 108;
    const unsigned int Hs = 108;
    const unsigned int Meitnerium = 109;
    const unsigned int Mt = 109;
    const unsigned int Darmstadtium = 110;
    const unsigned int Ds = 110;
    const unsigned int Roentgenium = 111;
    const unsigned int Rg = 111;
    const unsigned int Copernicium = 112;
    const unsigned int Cn = 112;
    const unsigned int Nihonium = 113;
    const unsigned int Nh = 113;
    const unsigned int Flerovium = 114;
    const unsigned int Fl = 114;
    const unsigned int Moscovium = 115;
    const unsigned int Mc = 115;
    const unsigned int Livermorium = 116;
    const unsigned int Lv = 116;
    const unsigned int Tennessine = 117;
    const unsigned int Ts = 117;
    const unsigned int Oganesson = 118;
    const unsigned int Og = 118;
  }
}

#endif //OB_ELEMENTS_H

//! \file elements.h
//! \brief Functions relating to elements

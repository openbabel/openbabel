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
    enum Element {
        Dummy = 0,
        Hydrogen = 1,
        H = 1,
        Helium = 2,
        He = 2,
        Lithium = 3,
        Li = 3,
        Beryllium = 4,
        Be = 4,
        Boron = 5,
        B = 5,
        Carbon = 6,
        C = 6,
        Nitrogen = 7,
        N = 7,
        Oxygen = 8,
        O = 8,
        Fluorine = 9, 
        F = 9,
        Neon = 10,
        Ne = 10,
        Sodium = 11,
        Na = 11,
        Magnesium = 12,
        Mg = 12,
        Aluminium = 13,
        Al = 13,
        Silicon = 14,
        Si = 14,
        Phosphorus = 15,
        P = 15,
        Sulfur = 16,
        S = 16,
        Chlorine = 17,
        Cl = 17,
        Argon = 18, 
        Ar = 18,
        Potassium = 19, 
        K = 19,
        Calcium = 20,
        Ca = 20,
        Scandium = 21,
        Sc = 21,
        Titanium = 22,
        Ti = 22,
        Vanadium = 23,
        V = 23,
        Chromium = 24,
        Cr = 24,
        Manganese = 25,
        Mn = 25,
        Iron = 26,
        Fe = 26,
        Cobalt = 27,
        Co = 27,
        Nickel = 28,
        Ni = 28,
        Copper = 29,
        Cu = 29,
        Zinc = 30,
        Zn = 30,
        Gallium = 31,
        Ga = 31,
        Germanium = 32,
        Ge = 32,
        Arsenic = 33,
        As = 33,
        Selenium = 34,
        Se = 34,
        Bromine = 35,
        Br = 35,
        Krypton = 36,
        Kr = 36,
        Rubidium = 37,
        Rb = 37,
        Strontium = 38,
        Sr = 38,
        Yttrium = 39,
        Y = 39,
        Zirconium = 40,
        Zr = 40,
        Niobium = 41,
        Nb = 41,
        Molybdenum = 42,
        Mo = 42,
        Technetium = 43,
        Tc = 43,
        Ruthenium = 44,
        Ru = 44,
        Rhodium = 45,
        Rh = 45,
        Palladium = 46,
        Pd = 46,
        Silver = 47,
        Ag = 47,
        Cadmium = 48,
        Cd = 48,
        Indium = 49,
        In = 49,
        Tin = 50,
        Sn = 50,
        Antimony = 51,
        Sb = 51,
        Tellurium = 52,
        Te = 52,
        Iodine = 53,
        I = 53,
        Xenon = 54,
        Xe = 54,
        Caesium = 55,
        Cs = 55,
        Barium = 56,
        Ba = 56,
        Lanthanum = 57,
        La = 57,
        Cerium = 58,
        Ce = 58,
        Praseodymium = 59,
        Pr = 59,
        Neodymium = 60,
        Nd = 60,
        Promethium = 61,
        Pm = 61,
        Samarium = 62,
        Sm = 62,
        Europium = 63,
        Eu = 63,
        Gadolinium = 64,
        Gd = 64,
        Terbium = 65,
        Tb = 65,
        Dysprosium = 66,
        Dy = 66,
        Holmium = 67,
        Ho = 67,
        Erbium = 68,
        Er = 68,
        Thulium = 69,
        Tm = 69,
        Ytterbium = 70,
        Yb = 70,
        Lutetium = 71,
        Lu = 71,
        Hafnium = 72,
        Hf = 72,
        Tantalum = 73,
        Ta = 73,
        Tungsten = 74,
        W = 74,
        Rhenium = 75,
        Re = 75,
        Osmium = 76,
        Os = 76,
        Iridium = 77,
        Ir = 77,
        Platinum = 78,
        Pt = 78,
        Gold = 79,
        Au = 79,
        Mercury = 80,
        Hg = 80,
        Thallium = 81,
        Tl = 81,
        Lead = 82,
        Pb = 82,
        Bismuth = 83,
        Bi = 83,
        Polonium = 84,
        Po = 84,
        Astatine = 85,
        At = 85,
        Radon = 86,
        Rn = 86,
        Francium = 87,
        Fr = 87,
        Radium = 88,
        Ra = 88,
        Actinium = 89,
        Ac = 89,
        Thorium = 90,
        Th = 90,
        Protactinium = 91,
        Pa = 91,
        Uranium = 92,
        U = 92,
        Neptunium = 93,
        Np = 93,
        Plutonium = 94,
        Pu = 94,
        Americium = 95,
        Am = 95,
        Curium = 96,
        Cm = 96,
        Berkelium = 97,
        Bk = 97,
        Californium = 98,
        Cf = 98,
        Einsteinium = 99,
        Es = 99,
        Fermium = 100,
        Fm = 100,
        Mendelevium = 101,
        Md = 101,
        Nobelium = 102,
        No = 102,
        Lawrencium = 103,
        Lr = 103,
        Rutherfordium = 104,
        Rf = 104,
        Dubnium = 105,
        Db = 105,
        Seaborgium = 106,
        Sg = 106,
        Bohrium = 107,
        Bh = 107,
        Hassium = 108,
        Hs = 108,
        Meitnerium = 109,
        Mt = 109,
        Darmstadtium = 110,
        Ds = 110,
        Roentgenium = 111,
        Rg = 111,
        Copernicium = 112,
        Cn = 112,
        Nihonium = 113,
        Nh = 113,
        Flerovium = 114,
        Fl = 114,
        Moscovium = 115,
        Mc = 115,
        Livermorium = 116,
        Lv = 116,
        Tennessine = 117,
        Ts = 117,
        Oganesson = 118,
        Og = 118
    };
  }
}

#endif //OB_ELEMENTS_H

//! \file elements.h
//! \brief Functions relating to elements

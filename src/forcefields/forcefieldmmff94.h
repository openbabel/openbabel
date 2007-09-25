/**********************************************************************
forcefieldmmff94.h - Merck Molecular Force Field (94).
 
Copyright (C) 2006 by Tim Vandermeersch <tim.vandermeersch@gmail.com>
 
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

#include <vector>
#include <string>
#include <map>

#include <openbabel/parsmart.h>
#include <openbabel/forcefield.h>
#include <openbabel/base.h>
#include <openbabel/mol.h>

namespace OpenBabel
{
  class OBFFBondCalculationMMFF94 : public OBFFCalculation
  {
    public:
      OBAtom *a, *b; // atoms of the bond
      double kb, r0, rab, delta;
      int bt; // bondtype (BTIJ)
      
      double GetEnergy();
      vector3 GetGradient(OBAtom *atom);
  };
  
  class OBFFAngleCalculationMMFF94 : public OBFFCalculation
  {
    public:
      OBAtom *a, *b, *c; // atoms of the angle
      double ka, theta0, theta, delta;
      int at; //angletype (ATIJK)
      
      double GetEnergy();
      vector3 GetGradient(OBAtom *atom);
  };
  
  class OBFFStrBndCalculationMMFF94 : public OBFFCalculation
  {
    public:
      OBAtom *a, *b, *c; // atoms of the angle
      double kbaABC, kbaCBA, theta0, theta, rab0, rbc0, rab, rbc, delta_theta, delta_rab, delta_rbc;
      int sbt; //strbndtype (SBTIJK)
      
      double GetEnergy();
      vector3 GetGradient(OBAtom *atom);
  };

  class OBFFTorsionCalculationMMFF94 : public OBFFCalculation
  {
    public:
      OBAtom *a, *b, *c, *d; // atoms of the torsion
      double v1, v2, v3, tor, cosine;
      int tt; //torsiontype (TTIJKL)
      
      double GetEnergy();
      vector3 GetGradient(OBAtom *atom);
  };

 class OBFFOOPCalculationMMFF94 : public OBFFCalculation
  {
    public:
      OBAtom *a, *b, *c, *d; // atoms of the oop angle
      double koop, angle;
      
      double GetEnergy();
  };

  class OBFFVDWCalculationMMFF94 : public OBFFCalculation
  {
    public:
      OBAtom *a, *b; // atoms of the pair
      double rab, epsilon, alpha_a, alpha_b, Na, Nb, Aa, Ab, Ga, Gb;
      double R_AB, R_AB7, erep, eattr, escale;
      int aDA, bDA; // hydrogen donor/acceptor (A=1, D=2, neither=0)

      double GetEnergy();
  };

  class OBFFElectrostaticCalculationMMFF94 : public OBFFCalculation
  {
    public:
      OBAtom *a, *b, *c; // atoms of the angle
      double ka, theta0, thetaabc, delta, delta2;
      int sbt; //strbndtype (SBTIJK)
      
      double GetEnergy();
  };

  // Class OBForceFieldMMFF94
  // class introduction in forcefieldmmff94.cpp
  class OBForceFieldMMFF94: public OBForceField
  {
    protected:
      //! \return Parses the parameter file
      bool ParseParamFile();
      bool ParseParamProp();
      bool ParseParamBond();
      bool ParseParamBndk();
      bool ParseParamAngle();
      bool ParseParamStrBnd();
      bool ParseParamDfsb();
      bool ParseParamOOP();
      bool ParseParamTorsion();
      bool ParseParamVDW();
      bool ParseParamCharge();
      bool ParseParamPbci();
      //! \return Sets atomtypes to MMFF94 in _mol
      bool SetMMFFTypes();
      //! fill OBFFXXXCalculation vectors
      bool SetupCalculations();
      //!  Sets charges to MMFF94 charges in _mol
      bool SetMMFF94Charges();
      //! \return The row of the element atom in the periodic table
      int GetElementRow(OBAtom *atom);
      //! \return The bond type (BTIJ)
      int GetBondType(OBAtom* a, OBAtom* b);
      //! \return The angle type (ATIJK)
      int GetAngleType(OBAtom* a, OBAtom* b, OBAtom *c);
      //! \return The stretch-bend type (SBTIJK)
      int GetStrBndType(OBAtom* a, OBAtom* b, OBAtom *c);
      //! \return The torsion type (TTIJKL)
      int GetTorsionType(OBAtom* a, OBAtom* b, OBAtom *c, OBAtom *d);
      //! \return true if atom a and b are in the same ring
      //bool IsInSameRing(OBAtom* a, OBAtom* b);
      //! \return true if atomtype has sbmb set in mmffprop.par
      bool HasSbmbSet(int atomtype);
      //! \return true if atomtype has arom set in mmffprop.par
      bool HasAromSet(int atomtype);
      //! \return true if atomtype has lin set in mmffprop.par
      bool HasLinSet(int atomtype);
      //! Same as OBForceField::GetParameter, but takes (bond/angle/torsion) type in account.
      OBFFParameter* GetParameterMMFF94(int ffclass, int a, int b, int c, int d, std::vector<OBFFParameter> &parameter);
      //! Returns the negative gradient (force) on atom a
      vector3 GetGradient(OBAtom *a, int terms = OBFF_ENERGY);
      
      // OBFFParameter vectors to contain the parameters
      std::vector<OBFFParameter> _ffbondparams;
      std::vector<OBFFParameter> _ffbndkparams;
      std::vector<OBFFParameter> _ffangleparams;
      std::vector<OBFFParameter> _ffstrbndparams;
      std::vector<OBFFParameter> _ffdfsbparams;
      std::vector<OBFFParameter> _fftorsionparams;
      std::vector<OBFFParameter> _ffoopparams;
      std::vector<OBFFParameter> _ffvdwparams;
      std::vector<OBFFParameter> _ffchgparams;
      std::vector<OBFFParameter> _ffpbciparams;

      // OBFFXXXCalculationYYY vectors to contain the calculations
      std::vector<OBFFBondCalculationMMFF94>          _bondcalculations;
      std::vector<OBFFAngleCalculationMMFF94>         _anglecalculations;
      std::vector<OBFFStrBndCalculationMMFF94>        _strbndcalculations;
      std::vector<OBFFTorsionCalculationMMFF94>       _torsioncalculations;
      std::vector<OBFFOOPCalculationMMFF94>           _oopcalculations;
      std::vector<OBFFVDWCalculationMMFF94>           _vdwcalculations;
      std::vector<OBFFElectrostaticCalculationMMFF94> _electrostaticcalculations;
      
      // parameters from mmffprop.par
      std::vector<int> _sbmb; // single bond - multiple bond
      std::vector<int> _arom; // aromatic
      std::vector<int> _lin; // lineair

    public:
      //! Constructor
      OBForceFieldMMFF94(const char* ID, bool IsDefault=true) : OBForceField(ID, IsDefault)
      {
        ParseParamFile();
      }
      
      //! Destructor
      virtual ~OBForceFieldMMFF94();
      
      //! Assignment
      OBForceFieldMMFF94 &operator = (OBForceFieldMMFF94 &);

      //! Get the description for this force field
      const char* Description()
      { 
        return "Merck Molecular Force Field. (94)";
      }

      //! Get the unit in wich the energy is expressed
      std::string GetUnit() 
      { 
        return std::string("kcal/mol"); 
      }
      //! Setup
      bool Setup(OBMol &mol);
 
      //! Returns total energy
      double Energy();
      //! Returns the bond stretching energy
      double E_Bond();
      //! Returns the angle bending energy
      double E_Angle();
      //! Returns the stretch-bend energy
      double E_StrBnd();
      //! Returns the torsional energy
      double E_Torsion();
      //! Returns the out-of-plane bending energy
      double E_OOP();
      //! Returns the Van der Waals energy (Buckingham potential)
      double E_VDW();
      //! Returns the dipole-dipole interaction energy
      double E_Electrostatic();
      
      
      //! Validate MMFF94 using validation suite
      bool Validate();
      //! Compare and print the numerical and analytical gradients
      bool ValidateGradients();
 
  }; // class OBForceFieldMM2

}// namespace OpenBabel

//! \file forcefieldmmff94.h
//! \brief MMFF94 force field


/**********************************************************************
forcefieldmmff94.h - MMFF94
 
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
      double kb, r0, rab, delta;
      int bt; // bondtype (BTIJ)
      
      void Compute(bool gradients = true);
  };
  
  class OBFFAngleCalculationMMFF94 : public OBFFCalculation
  {
    public:
      double ka, theta0, theta, delta;
      int at; //angletype (ATIJK)
      
      void Compute(bool gradients = true);
  };
  
  class OBFFStrBndCalculationMMFF94 : public OBFFCalculation
  {
    public:
      double kbaABC, kbaCBA, theta0, theta, rab0, rbc0, rab, rbc, delta_theta, delta_rab, delta_rbc;
      int sbt; //strbndtype (SBTIJK)
      
      void Compute(bool gradients = true);
  };

  class OBFFTorsionCalculationMMFF94 : public OBFFCalculation
  {
    public:
      double v1, v2, v3, tor, cosine;
      int tt; //torsiontype (TTIJKL)
      
      void Compute(bool gradients = true);
  };

 class OBFFOOPCalculationMMFF94 : public OBFFCalculation
  {
    public:
      double koop, angle;
      
      void Compute(bool gradients = true);
  };

  class OBFFVDWCalculationMMFF94 : public OBFFCalculation
  {
    public:
      double rab, epsilon, alpha_a, alpha_b, Na, Nb, Aa, Ab, Ga, Gb;
      double R_AB, R_AB7, erep, eattr, escale;
      int aDA, bDA; // hydrogen donor/acceptor (A=1, D=2, neither=0)

      void Compute(bool gradients = true);
  };

  class OBFFElectrostaticCalculationMMFF94 : public OBFFCalculation
  {
    public:
      double qq, rab;
      
      void Compute(bool gradients = true);
  };

  // Class OBForceFieldMMFF94
  // class introduction in forcefieldmmff94.cpp
  class OBForceFieldMMFF94: public OBForceField
  {
    protected:
      //! \return Parses the parameter file
      bool ParseParamFile();
      bool ParseParamProp();
      bool ParseParamDef();
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
      //! detect which rings are aromatic
      bool PerceiveAromatic();
      //! \return Sets atomtypes to MMFF94 in _mol
      bool SetTypes();
      //! fill OBFFXXXCalculation vectors
      bool SetupCalculations();
      //!  Sets formal charges
      bool SetFormalCharges();
      //!  Sets partial charges
      bool SetPartialCharges();
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
      //! \return true if atomtype has sbmb set in mmffprop.par
      bool HasSbmbSet(int atomtype);
      //! \return true if atomtype has pilp set in mmffprop.par
      bool HasPilpSet(int atomtype);
      //! \return true if atomtype has arom set in mmffprop.par
      bool HasAromSet(int atomtype);
      //! \return true if atomtype has lin set in mmffprop.par
      bool HasLinSet(int atomtype);
      //! \return the crd value for the atomtype in mmffprop.par
      int GetCrd(int atomtype);
      //! \return the val value for the atomtype in mmffprop.par
      int GetVal(int atomtype);
      //! \return the mltb value for the atomtype in mmffprop.par
      int GetMltb(int atomtype);
      //! \return the level 2 equivalent atom type for type (mmffdef.par)
      int EqLvl2(int type);
      //! \return the level 3 equivalent atom type for type (mmffdef.par)
      int EqLvl3(int type);
      //! \return the level 4 equivalent atom type for type (mmffdef.par)
      int EqLvl4(int type);
      //! \return the level 5 equivalent atom type for type (mmffdef.par)
      int EqLvl5(int type);
      //! \return the U value for the atom from table X page 631
      double GetUParam(OBAtom* atom);
      //! \return the Z value for the atom from table VI page 628
      double GetZParam(OBAtom* atom);
      //! \return the C value for the atom from table VI page 628
      double GetCParam(OBAtom* atom);
      //! \return the V value for the atom from table X page 631
      double GetVParam(OBAtom* atom);
      //! return the covalent radius from Blom and Haaland, value from etab if not available
      double GetCovalentRadius(OBAtom* a);
      //! return the bond length calculated with a modified version of the Schomaker-Stevenson rule
      double GetRuleBondLength(OBAtom* a, OBAtom* b);
      //! return the bond length from mmffbond.par, if not found, one is calculated with a modified version of the Schomaker-Stevenson rule
      double GetBondLength(OBAtom* a, OBAtom* b);
      
      //! Same as OBForceField::GetParameter, but takes (bond/angle/torsion) type in account and takes 0 as wildcart.
      OBFFParameter* GetParameter1Atom(int a, std::vector<OBFFParameter> &parameter);
      OBFFParameter* GetParameter2Atom(int a, int b, std::vector<OBFFParameter> &parameter);
      OBFFParameter* GetParameter3Atom(int a, int b, int c, std::vector<OBFFParameter> &parameter);
      
      //! Same as OBForceField::GetParameter, but takes (bond/angle/torsion) type in account and takes 0 as wildcart.
      OBFFParameter* GetTypedParameter2Atom(int ffclass, int a, int b, std::vector<OBFFParameter> &parameter);
      OBFFParameter* GetTypedParameter3Atom(int ffclass, int a, int b, int c, std::vector<OBFFParameter> &parameter);
      OBFFParameter* GetTypedParameter4Atom(int ffclass, int a, int b, int c, int d, std::vector<OBFFParameter> &parameter);
      
      
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
      std::vector<OBFFParameter> _ffdefparams;
      std::vector<OBFFParameter> _ffpropparams;

      // OBFFXXXCalculationYYY vectors to contain the calculations
      std::vector<OBFFBondCalculationMMFF94>          _bondcalculations;
      std::vector<OBFFAngleCalculationMMFF94>         _anglecalculations;
      std::vector<OBFFStrBndCalculationMMFF94>        _strbndcalculations;
      std::vector<OBFFTorsionCalculationMMFF94>       _torsioncalculations;
      std::vector<OBFFOOPCalculationMMFF94>           _oopcalculations;
      std::vector<OBFFVDWCalculationMMFF94>           _vdwcalculations;
      std::vector<OBFFElectrostaticCalculationMMFF94> _electrostaticcalculations;
      
    public:
      //! Constructor
      OBForceFieldMMFF94(const char* ID, bool IsDefault=true) : OBForceField(ID, IsDefault)
      {
        _init = false;
      }
      
      //! Destructor
      virtual ~OBForceFieldMMFF94();
      
      //! Assignment
      OBForceFieldMMFF94 &operator = (OBForceFieldMMFF94 &);
      
      //!Clone the current instance. May be desirable in multithreaded environments
      virtual OBForceFieldMMFF94* MakeNewInstance(){ return new OBForceFieldMMFF94(*this); }

      //! Get the description for this force field
      const char* Description()
      { 
        return "MMFF94 force field.";
      }

      //! Get the unit in wich the energy is expressed
      std::string GetUnit() 
      { 
        return std::string("kcal/mol"); 
      }
      //! Setup
      bool Setup(OBMol &mol);
 
      //! Returns total energy
      double Energy(bool gradients = true);
      //! Returns the bond stretching energy
      double E_Bond(bool gradients = true);
      //! Returns the angle bending energy
      double E_Angle(bool gradients = true);
      //! Returns the stretch-bend energy
      double E_StrBnd(bool gradients = true);
      //! Returns the torsional energy
      double E_Torsion(bool gradients = true);
      //! Returns the out-of-plane bending energy
      double E_OOP(bool gradients = true);
      //! Returns the Van der Waals energy (Buckingham potential)
      double E_VDW(bool gradients = true);
      //! Returns the dipole-dipole interaction energy
      double E_Electrostatic(bool gradients = true);
      
      
      //! Validate MMFF94 using validation suite
      bool Validate();
      //! Compare and print the numerical and analytical gradients
      bool ValidateGradients();
 
  }; // class OBForceFieldMM2

}// namespace OpenBabel

//! \file forcefieldmmff94.h
//! \brief MMFF94 force field


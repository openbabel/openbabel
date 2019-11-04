/**********************************************************************
forcefieldmmff94.h - MMFF94

Copyright (C) 2006 by Tim Vandermeersch <tim.vandermeersch@gmail.com>

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

#include <vector>
#include <string>
#include <map>

#include <openbabel/parsmart.h>
#include <openbabel/forcefield.h>
#include <openbabel/base.h>
#include <openbabel/mol.h>

namespace OpenBabel
{
  class OBFFBondCalculationMMFF94 : public OBFFCalculation2
  {
    public:
      int bt; // bondtype (BTIJ)
      double kb, r0, rab, delta;

      template<bool> void Compute();
  };

  class OBFFAngleCalculationMMFF94 : public OBFFCalculation3
  {
    public:
      int at; //angletype (ATIJK)
      bool linear;
      double ka, theta, theta0, delta;

      template<bool> void Compute();
  };

  class OBFFStrBndCalculationMMFF94 : public OBFFCalculation3
  {
    public:
      int sbt; //strbndtype (SBTIJK)
      double kbaABC, kbaCBA, theta0, rab0, rbc0, delta_theta, delta_rab, delta_rbc;
      double theta, rab, rbc;
      double force_ab_a[3], force_ab_b[3], force_bc_b[3], force_bc_c[3];
      double force_abc_a[3], force_abc_b[3], force_abc_c[3];

      template<bool> void Compute();
  };

  class OBFFTorsionCalculationMMFF94 : public OBFFCalculation4
  {
    public:
      int tt; //torsiontype (TTIJKL)
      double v1, v2, v3, tor, cosine;

      template<bool> void Compute();
  };

  class OBFFOOPCalculationMMFF94 : public OBFFCalculation4
  {
    public:
      double koop, angle;

      template<bool> void Compute();
  };

  class OBFFVDWCalculationMMFF94 : public OBFFCalculation2
  {
    public:
      int aDA, bDA; // hydrogen donor/acceptor (A=1, D=2, neither=0)
      double rab, epsilon, alpha_a, alpha_b, Na, Nb, Aa, Ab, Ga, Gb;
      double R_AB, R_AB7/*, erep, erep7, eattr*/;
      int pairIndex; // index into iteration using FOR_PAIRS_OF_MOL(..., _mol)

      template<bool> void Compute();
  };

  class OBFFElectrostaticCalculationMMFF94 : public OBFFCalculation2
  {
    public:
      double qq, rab;
      int pairIndex; // index into iteration using FOR_PAIRS_OF_MOL(..., _mol)

      template<bool> void Compute();
  };

  // Class OBForceFieldMMFF94
  // class introduction in forcefieldmmff94.cpp
  class OBForceFieldMMFF94: public OBForceField
  {
    protected:
      //! \return Parses the parameter file
      bool ParseParamFile();
      bool ParseParamProp(std::string &filename);
      bool ParseParamDef(std::string &filename);
      bool ParseParamBond(std::string &filename);
      bool ParseParamBndk(std::string &filename);
      bool ParseParamAngle(std::string &filename);
      bool ParseParamStrBnd(std::string &filename);
      bool ParseParamDfsb(std::string &filename);
      bool ParseParamOOP(std::string &filename);
      bool ParseParamTorsion(std::string &filename);
      bool ParseParamVDW(std::string &filename);
      bool ParseParamCharge(std::string &filename);
      bool ParseParamPbci(std::string &filename);
      //! detect which rings are aromatic
      bool PerceiveAromatic();
      //! \return Get the MMFF94 atom type for atom
      int GetType(OBAtom *atom);
      //! \return Sets atomtypes to MMFF94 in _mol
      bool SetTypes();
      //! fill OBFFXXXCalculation vectors
      bool SetupCalculations();
      //! Setup pointers in OBFFXXXCalculation vectors
      bool SetupPointers();
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
      //! \return the canonical bond index
      unsigned int GetCXB(int type, int a, int b);
      //! \return the canonical angle index
      unsigned int GetCXA(int type, int a, int b, int c);
      //! \return the canonical stretch-bend index
      unsigned int GetCXS(int type, int a, int b, int c);
      //! \return the canonical out-of-plane index
      unsigned int GetCXO(int a, int b, int c, int d);
      //! \return the canonical torsion index
      unsigned int GetCXT(int type, int a, int b, int c, int d);
      //! \return the canonical bond-charge-increment index
      unsigned int GetCXQ(int type, int a, int b);
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
      OBBitVec			 _ffpropPilp;
      OBBitVec			 _ffpropArom;
      OBBitVec			 _ffpropLin;
      OBBitVec			 _ffpropSbmb;

      // OBFFXXXCalculationYYY vectors to contain the calculations
      std::vector<OBFFBondCalculationMMFF94>          _bondcalculations;
      std::vector<OBFFAngleCalculationMMFF94>         _anglecalculations;
      std::vector<OBFFStrBndCalculationMMFF94>        _strbndcalculations;
      std::vector<OBFFTorsionCalculationMMFF94>       _torsioncalculations;
      std::vector<OBFFOOPCalculationMMFF94>           _oopcalculations;
      std::vector<OBFFVDWCalculationMMFF94>           _vdwcalculations;
      std::vector<OBFFElectrostaticCalculationMMFF94> _electrostaticcalculations;

      bool mmff94s;

    public:
      //! Constructor
      explicit OBForceFieldMMFF94(const char* ID, bool IsDefault=true) : OBForceField(ID, IsDefault)
      {
        _validSetup = false;
        _init = false;
        _rvdw = 7.0;
        _rele = 15.0;
        _epsilon = 1.0; // default electrostatics
        _pairfreq = 15;
        _cutoff = false;
        _linesearch = LineSearchType::Newton2Num;
        _gradientPtr = NULL;
        _grad1 = NULL;
	if (!strncmp(ID, "MMFF94s", 7)) {
          mmff94s = true;
          _parFile = std::string("mmff94s.ff");
	} else {
          mmff94s = false;
          _parFile = std::string("mmff94.ff");
	}
      }

      //! Destructor
      virtual ~OBForceFieldMMFF94();

      //! Assignment
      OBForceFieldMMFF94 &operator = (OBForceFieldMMFF94 &);

      //!Clone the current instance. May be desirable in multithreaded environments
      virtual OBForceFieldMMFF94* MakeNewInstance()
      {
        return new OBForceFieldMMFF94(_id, false);
      }

      //! Get the description for this force field
      const char* Description()
      {
        if (mmff94s)
          return "MMFF94s force field.";
	else
          return "MMFF94 force field.";
      }

      //! Get the unit in which the energy is expressed
      std::string GetUnit()
      {
        return std::string("kcal/mol");
      }

      //! \return that analytical gradients are implemented for MMFF94
      bool HasAnalyticalGradients() { return true; }

      //! Returns total energy
      double Energy(bool gradients = true);
      //! Returns the bond stretching energy
      template<bool> double E_Bond();
      double E_Bond(bool gradients = true)
      {
        return gradients ? E_Bond<true>() : E_Bond<false>();
      }
      //! Returns the angle bending energy
      template<bool> double E_Angle();
      double E_Angle(bool gradients = true)
      {
        return gradients ? E_Angle<true>() : E_Angle<false>();
      }
      //! Returns the stretch-bend energy
      template<bool> double E_StrBnd();
      double E_StrBnd(bool gradients = true)
      {
        return gradients ? E_StrBnd<true>() : E_StrBnd<false>();
      }
      //! Returns the torsional energy
      template<bool> double E_Torsion();
      double E_Torsion(bool gradients = true)
      {
        return gradients ? E_Torsion<true>() : E_Torsion<false>();
      }
      //! Returns the out-of-plane bending energy
      template<bool> double E_OOP();
      double E_OOP(bool gradients = true)
      {
        return gradients ? E_OOP<true>() : E_OOP<false>();
      }
      //! Returns the Van der Waals energy (Buckingham potential)
      template<bool> double E_VDW();
      double E_VDW(bool gradients = true)
      {
        return gradients ? E_VDW<true>() : E_VDW<false>();
      }
      //! Returns the dipole-dipole interaction energy
      template<bool> double E_Electrostatic();
      double E_Electrostatic(bool gradients = true)
      {
        return gradients ? E_Electrostatic<true>() : E_Electrostatic<false>();
      }

      //! Validate MMFF94 using validation suite
      bool Validate();
      //! Compare and print the numerical and analytical gradients
      bool ValidateGradients();

  }; // class OBForceFieldMM2

}// namespace OpenBabel

//! \file forcefieldmmff94.h
//! \brief MMFF94 force field

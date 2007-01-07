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
      double kb, r0, e, rab, delta;
      int bt; // bondtype (BTIJ)
      
      double Result();
  };
  
  class OBFFAngleCalculationMMFF94 : public OBFFCalculation
  {
    public:
      OBAtom *a, *b, *c; // atoms of the angle
      double ka, theta0, e, theta, delta;
      int at; //angletype (ATIJK)
      
      double Result();
  };
  
  class OBFFStrBndCalculationMMFF94 : public OBFFCalculation
  {
    public:
      OBAtom *a, *b, *c; // atoms of the angle
      double kbaABC, kbaCBA, e, theta0, theta, rab0, rbc0, rab, rbc, delta_theta, delta_rab, delta_rbc;
      int sbt; //strbndtype (SBTIJK)
      
      double Result();
  };

  class OBFFTorsionCalculationMMFF94 : public OBFFCalculation
  {
    public:
      OBAtom *a, *b, *c, *d; // atoms of the torsion
      double v1, v2, v3, e, tor, cosine;
      int tt; //torsiontype (TTIJKL)
      
      double Result();
  };

 class OBFFOOPCalculationMMFF94 : public OBFFCalculation
  {
    public:
      OBAtom *a, *b, *c, *d; // atoms of the oop angle
      double koop, e, angle;
      
      double Result();
  };

  class OBFFVDWCalculationMMFF94 : public OBFFCalculation
  {
    public:
      OBAtom *a, *b; // atoms of the pair
      double e, rab, epsilon, alpha_a, alpha_b, Na, Nb, Aa, Ab, Ga, Gb;
      double R_AB, R_AB7, erep, eattr, escale;
      int aDA, bDA; // hydrogen donor/acceptor (A=1, D=2, neither=0)

      double Result();
  };

  class OBFFElectrostaticCalculationMMFF94 : public OBFFCalculation
  {
    public:
      OBAtom *a, *b, *c; // atoms of the angle
      double ka, theta0, e, thetaabc, delta, delta2;
      int sbt; //strbndtype (SBTIJK)
      
      double Result();
  };

  // Class OBForceFieldMMFF94
  // class introduction in forcefieldmmff94.cpp
  class OBAPI OBForceFieldMMFF94: public OBForceField
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
      //! \return Sets atomtypes to MM2 in _mol
      bool SetMMFFTypes();
      //! fill OBFFXXXCalculation vectors
      bool SetupCalculations();
      bool CalcCharges();
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
      bool IsInSameRing(OBAtom* a, OBAtom* b);
      //! \return true if atomtype has sbmb set in mmffprop.par
      bool HasSbmbSet(int atomtype);
      //! \return true if atomtype has arom set in mmffprop.par
      bool HasAromSet(int atomtype);
      //! \return true if atomtype has lin set in mmffprop.par
      bool HasLinSet(int atomtype);
      
      
      bool Validate();
      OBFFParameter* GetParameterMMFF94(int ffclass, int a, int b, int c, int d, std::vector<OBFFParameter> &parameter);

      
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

      std::vector<OBFFBondCalculationMMFF94>          _bondcalculations;
      std::vector<OBFFAngleCalculationMMFF94>         _anglecalculations;
      std::vector<OBFFStrBndCalculationMMFF94>        _strbndcalculations;
      std::vector<OBFFTorsionCalculationMMFF94>       _torsioncalculations;
      std::vector<OBFFOOPCalculationMMFF94>           _oopcalculations;
      std::vector<OBFFVDWCalculationMMFF94>           _vdwcalculations;
      std::vector<OBFFElectrostaticCalculationMMFF94> _electostaticcalculations;
      
      // parameters from mmffprop.par
      std::vector<int> _sbmb; // single bond - multiple bond
      std::vector<int> _arom; // aromatic
      std::vector<int> _lin; // lineair

    public:
      //! Setup
      bool Setup(OBMol &mol);
      //! Constructor
      OBForceFieldMMFF94(std::string ID, bool IsDefault=true) : OBForceField(ID, IsDefault)
      {
        ParseParamFile();
      }
      
      std::string Description()
	{ return "Merck Molecular Force Field. (94)";};

      std::string GetUnit() { return std::string("kcal/mol"); }


      //! Destructor
      virtual ~OBForceFieldMMFF94();
      //! Assignment
      OBForceFieldMMFF94 &operator = (OBForceFieldMMFF94 &);
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

  }; // class OBForceFieldMM2

}// namespace OpenBabel

//! \file forcefieldmmff94.h
//! \brief MMFF94 force field


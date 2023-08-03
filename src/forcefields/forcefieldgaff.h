/**********************************************************************
forcefieldgaff.h - Gaff force field.

Copyright (C) 2009 by Frank Peters <e.a.j.f.peters@tue.nl>
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

#include <openbabel/forcefield.h>
#include <openbabel/base.h>
#include <openbabel/mol.h>

namespace OpenBabel
{

  class OBFFBondCalculationGaff : public OBFFCalculation2
  {
    public:
      double kr, r0, rab, delta;

      template<bool> void Compute();
  };

  class OBFFAngleCalculationGaff : public OBFFCalculation3
  {
    public:
      double kth, theta, theta0, delta;

      template<bool> void Compute();
  };

  class OBFFTorsionCalculationGaff : public OBFFCalculation4
  {
    public:
    double tor, vn_half, gamma, n;

      template<bool> void Compute();
  };

  class OBFFOOPCalculationGaff : public OBFFCalculation4
  {
    public:
    double tor, vn_half, gamma, n;

      template<bool> void Compute();
  };

  class OBFFVDWCalculationGaff : public OBFFCalculation2
  {
    public:
      bool is14, samering;
      double Eab, RVDWab, rab;

      template<bool> void Compute();
  };

  class OBFFElectrostaticCalculationGaff : public OBFFCalculation2
  {
    public:
      double qq, rab;

      template<bool> void Compute();
  };

  // Class OBForceFieldGaff
  // class introduction in forcefieldgaff.cpp
  class OBForceFieldGaff: public OBForceField
  {
    protected:
      //!  Parses the parameter file
      bool ParseParamFile();
      //!  Sets atomtypes to Gaff types in _mol
      bool SetTypes();
      //!  Sets partial charges to Gaff charges in _mol
      bool SetPartialCharges();
      //! fill OBFFXXXCalculation vectors
      bool SetupCalculations();
      //! Setup pointers in OBFFXXXCalculation vectors
      bool SetupPointers();
      //! Calculate Gasteiger charges 'out of order' before atom typing
      bool SetPartialChargesBeforeAtomTyping();
      // GetParameterOOP for improper-dihedrals
      // This specialization is needed because improper-dihedral have different symmetry as dihedrals
      OBFFParameter* GetParameterOOP(const char* a, const char* b, const char* c, const char* d,
        std::vector<OBFFParameter> &parameter);

      // OBFFParameter vectors to contain the parameters
      std::vector<OBFFParameter> _ffpropparams;
      std::vector<OBFFParameter> _ffbondparams;
      std::vector<OBFFParameter> _ffangleparams;
      std::vector<OBFFParameter> _fftorsionparams;
      std::vector<OBFFParameter> _ffoopparams;
      std::vector<OBFFParameter> _ffhbondparams;
      std::vector<OBFFParameter> _ffvdwparams;
      std::vector<OBFFParameter> _ffchargeparams;


      // OBFFXXXCalculationYYY vectors to contain the calculations
      std::vector<OBFFBondCalculationGaff>          _bondcalculations;
      std::vector<OBFFAngleCalculationGaff>         _anglecalculations;
      std::vector<OBFFTorsionCalculationGaff>       _torsioncalculations;
      std::vector<OBFFOOPCalculationGaff>      _oopcalculations;
      std::vector<OBFFVDWCalculationGaff>           _vdwcalculations;
      std::vector<OBFFElectrostaticCalculationGaff> _electrostaticcalculations;

    public:
      //! Constructor
      explicit OBForceFieldGaff(const char* ID, bool IsDefault=true) : OBForceField(ID, IsDefault)
      {
        _validSetup = false;
        _init = false;
        _rvdw = 7.0;
        _rele = 15.0;
        _epsilon = 1.0;
        _pairfreq = 10;
        _cutoff = false;
        _linesearch = LineSearchType::Newton2Num;
      }

      //! Destructor
      virtual ~OBForceFieldGaff();

      //! Assignment
      OBForceFieldGaff &operator = (OBForceFieldGaff &);

      //! Get the description for this force field
      const char* Description()
      {
        return "General Amber Force Field (GAFF).";
      }

      //!Clone the current instance. May be desirable in multithreaded environments
      virtual OBForceFieldGaff* MakeNewInstance()
      {
        return new OBForceFieldGaff(_id, false);
      }

      //! Get the unit in which the energy is expressed
      std::string GetUnit()
      {
        return std::string("kJ/mol");
      }

      //! \return that analytical gradients are implemented for Gaff
      bool HasAnalyticalGradients() { return true; }

      //! Setup
      bool Setup(OBMol &mol);

      //! \return total energy
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
      //! Returns the torsional energy
      template<bool> double E_Torsion();
      double E_Torsion(bool gradients = true)
      {
        return gradients ? E_Torsion<true>() : E_Torsion<false>();
      }
      //! Returns the improper (out-of-plane) bending energy
      template<bool> double E_OOP();
      double E_OOP(bool gradients = true)
      {
        return gradients ? E_OOP<true>() : E_OOP<false>();
      }
      //! Returns the Van der Waals energy
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

      //! Compare and print the numerical and analytical gradients
      bool ValidateGradients();

  }; // class OBForceFieldGaff

}// namespace OpenBabel

//! \file forcefieldgaff.h
//! \brief Gaff force field

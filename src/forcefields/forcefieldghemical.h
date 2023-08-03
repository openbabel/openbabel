/**********************************************************************
forcefieldghemical.h - Ghemical force field.

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
  class OBFFBondCalculationGhemical : public OBFFCalculation2
  {
    public:
      int bt; // bondtype (BTIJ)
      double kb, r0, rab, delta;

      template<bool> void Compute();
  };

  class OBFFAngleCalculationGhemical : public OBFFCalculation3
  {
    public:
      double ka, theta, theta0, delta;

      template<bool> void Compute();
  };

  class OBFFTorsionCalculationGhemical : public OBFFCalculation4
  {
    public:
      int tt; //torsiontype (TTIJKL)
      double V, s, n, tor;
      double k1, k2, k3;

      template<bool> void Compute();
  };

  class OBFFVDWCalculationGhemical : public OBFFCalculation2
  {
    public:
      bool is14, samering;
      double Ra, Rb, kab, rab;
      union {
        double ka, sigma12;
      };
      union {
        double kb, sigma6;
      };

      template<bool> void Compute();
  };

  class OBFFElectrostaticCalculationGhemical : public OBFFCalculation2
  {
    public:
      double qq, rab;

      template<bool> void Compute();
  };

  // Class OBForceFieldGhemical
  // class introduction in forcefieldghemical.cpp
  class OBForceFieldGhemical: public OBForceField
  {
    protected:
      //!  Parses the parameter file
      bool ParseParamFile();
      //!  Sets atomtypes to Ghemical types in _mol
      bool SetTypes();
      //!  Sets partial charges to Ghemical charges in _mol
      bool SetPartialCharges();
      //! fill OBFFXXXCalculation vectors
      bool SetupCalculations();
      //! Setup pointers in OBFFXXXCalculation vectors
      bool SetupPointers();
      //! Same as OBForceField::GetParameter, but takes (bond/angle/torsion) type in account.
      OBFFParameter* GetParameterGhemical(int type, const char* a, const char* b,
          const char* c, const char* d, std::vector<OBFFParameter> &parameter);

      // OBFFParameter vectors to contain the parameters
      std::vector<OBFFParameter> _ffbondparams;
      std::vector<OBFFParameter> _ffangleparams;
      std::vector<OBFFParameter> _fftorsionparams;
      std::vector<OBFFParameter> _ffvdwparams;
      std::vector<OBFFParameter> _ffchargeparams;

      // OBFFXXXCalculationYYY vectors to contain the calculations
      std::vector<OBFFBondCalculationGhemical>          _bondcalculations;
      std::vector<OBFFAngleCalculationGhemical>         _anglecalculations;
      std::vector<OBFFTorsionCalculationGhemical>       _torsioncalculations;
      std::vector<OBFFVDWCalculationGhemical>           _vdwcalculations;
      std::vector<OBFFElectrostaticCalculationGhemical> _electrostaticcalculations;

    public:
      //! Constructor
      explicit OBForceFieldGhemical(const char* ID, bool IsDefault=true) : OBForceField(ID, IsDefault)
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
      virtual ~OBForceFieldGhemical();

      //! Assignment
      OBForceFieldGhemical &operator = (OBForceFieldGhemical &);

      //! Get the description for this force field
      const char* Description()
      {
        return "Ghemical force field.";
      }

      //!Clone the current instance. May be desirable in multithreaded environments
      virtual OBForceFieldGhemical* MakeNewInstance()
      {
        return new OBForceFieldGhemical(_id, false);
      }

      //! Get the unit in which the energy is expressed
      std::string GetUnit()
      {
        return std::string("kJ/mol");
      }

      //! \return that analytical gradients are implemented for Ghemical
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

      //! Compare and print the numerical and analytical gradients
      bool ValidateGradients();

  }; // class OBForceFieldGhemical

}// namespace OpenBabel

//! \file forcefieldghemical.h
//! \brief Ghemical force field

/**********************************************************************
forcefielduff.h - UFF force field.

Copyright (C) 2007 by Geoffrey Hutchison
Some portions Copyright (C) 2006-2007 by Tim Vandermeersch

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
  class OBFFBondCalculationUFF : public OBFFCalculation2
  {
    public:
      double bt; // bond order (e.g., 1.41 for amide)
      double kb, r0, rab, delta;

      template<bool> void Compute();
  };

  class OBFFAngleCalculationUFF : public OBFFCalculation3
  {
    public:
      int at; //angletype (ATIJK)
      bool linear;
      double ka, theta0, theta, delta;
      double c0, c1, c2;
      double zi, zk, rij, rjk, rik;
      double cosT0; // cos theta0
      int coord, n;

      template<bool> void Compute();
  };

  class OBFFTorsionCalculationUFF : public OBFFCalculation4
  {
    public:
      int n;
      double tt; //torsiontype (i.e. b-c bond order)
      double V, tor, cosNPhi0;

      template<bool> void Compute();

  };

  class OBFFOOPCalculationUFF : public OBFFCalculation4
  {
    public:
      double koop, angle;
      double c0, c1, c2;

      template<bool> void Compute();
  };

  class OBFFVDWCalculationUFF : public OBFFCalculation2
  {
    public:
      bool is14, samering;
      double ka, kaSquared, Ra, kb, Rb, kab, rab;

      template<bool> void Compute();
  };

  class OBFFElectrostaticCalculationUFF : public OBFFCalculation2
  {
    public:
      double qq, rab;

      template<bool> void Compute();
  };

  // Class OBForceFieldUFF
  // class introduction in forcefieldUFF.cpp
  class OBForceFieldUFF: public OBForceField
  {
  protected:
    //!  Parses the parameter file
    bool ParseParamFile();
    //!  Sets atomtypes to UFF types in _mol
    bool SetTypes();
    //!  Fill OBFFXXXCalculation vectors
    bool SetupCalculations();
    //! Setup pointers in OBFFXXXCalculation vectors
    bool SetupPointers();
    bool SetupVDWCalculation(OBAtom *a, OBAtom *b, OBFFVDWCalculationUFF &vdwcalc);
    //!  By default, electrostatic terms are disabled
    //!  This is discouraged, since the parameterization is not designed for it
    //!  But if you want, we give you the option.
    bool SetupElectrostatics();
    //! Same as OBForceField::GetParameter, but simpler
    OBFFParameter* GetParameterUFF(std::string a, std::vector<OBFFParameter> &parameter);

    // OBFFParameter vectors to contain the parameters
    std::vector<OBFFParameter> _ffparams;

    // OBFFXXXCalculationYYY vectors to contain the calculations
    std::vector<OBFFBondCalculationUFF>          _bondcalculations;
    std::vector<OBFFAngleCalculationUFF>         _anglecalculations;
    std::vector<OBFFTorsionCalculationUFF>       _torsioncalculations;
    std::vector<OBFFOOPCalculationUFF>           _oopcalculations;
    std::vector<OBFFVDWCalculationUFF>           _vdwcalculations;
    std::vector<OBFFElectrostaticCalculationUFF> _electrostaticcalculations;

  public:
    //! Constructor
    explicit OBForceFieldUFF(const char* ID, bool IsDefault=true) : OBForceField(ID, IsDefault)
    {
      _validSetup = false;
      _init = false;
      _rvdw = 7.0;
      _rele = 15.0;
      _epsilon = 1.0; // electrostatics not used
      _pairfreq = 10;
      _cutoff = false;
      _linesearch = LineSearchType::Newton2Num;
    }

    //! Destructor
    virtual ~OBForceFieldUFF();

     //!Clone the current instance. May be desirable in multithreaded environments
     virtual OBForceFieldUFF* MakeNewInstance()
     {
       return new OBForceFieldUFF(_id, false);
     }

    //! Assignment
    OBForceFieldUFF &operator = (OBForceFieldUFF &);

    //! Get the description for this force field
    const char* Description()
    {
      return "Universal Force Field.";
    }

    //! Get the unit in which the energy is expressed
    std::string GetUnit()
      {
        return std::string("kJ/mol");  // Note that we convert from kcal/mol internally
      }

    //! \return that analytical gradients are implemented for UFF
    bool HasAnalyticalGradients() { return true; }

    //! \return total energy
    double Energy(bool gradients = true);
    //! \return the bond stretching energy
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

    //! Compare and print the numerical and analytical gradients
    bool ValidateGradients();

  }; // class OBForceFieldUFF

}// namespace OpenBabel

//! \file forcefieldUFF.h
//! \brief UFF force field

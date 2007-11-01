/**********************************************************************
forcefielduff.h - UFF force field.
 
Copyright (C) 2007 by Geoffrey Hutchison
Some portions Copyright (C) 2006-2007 by Tim Vandermeersch
 
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

#include <openbabel/forcefield.h>
#include <openbabel/base.h>
#include <openbabel/mol.h>

namespace OpenBabel
{
  class OBFFBondCalculationUFF : public OBFFCalculation
  {
  public:
    double kb, r0, rab, delta;
    int bt; // bondtype
      
    void Compute(bool gradients = true);
  };
  
  class OBFFAngleCalculationUFF : public OBFFCalculation
  {
  public:
    double ka, theta0, theta, delta;
    double c1, c2, c0;
      
    void Compute(bool gradients = true);
  };
  
  class OBFFTorsionCalculationUFF : public OBFFCalculation
  {
  public:
    double V, s, n, tor;
    double k1, k2, k3;
    int tt; //torsiontype
      
    void Compute(bool gradients = true);
  };

  class OBFFOOPCalculationUFF : public OBFFCalculation
  {
  public:
    double koop, angle;
      
    void Compute(bool gradients = true);
  };  
  
  class OBFFVDWCalculationUFF : public OBFFCalculation
  {
  public:
    double ka, Ra, kb, Rb, kab, rab;
    bool is14, samering;

    void Compute(bool gradients = true);
  };

  class OBFFElectrostaticCalculationUFF : public OBFFCalculation
  {
  public:
    double qq, rab;
      
    void Compute(bool gradients = true);
  };

  // Class OBForceFieldUFF
  // class introduction in forcefieldUFF.cpp
  class OBForceFieldUFF: public OBForceField
  {
  protected:
    //!  Parses the parameter file
    bool ParseParamFile();
    //!  Sets atomtypes to UFF types in _mol
    bool SetUFFTypes();
    //! fill OBFFXXXCalculation vectors
    bool SetupCalculations();
    //! Same as OBForceField::GetParameter, but simpler
    OBFFParameter* GetParameterUFF(std::string a, std::vector<OBFFParameter> &parameter);
    //! Returns the negative gradient (force) on atom a
    vector3 GetGradient(OBAtom *a, int terms = OBFF_ENERGY);
      
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
  OBForceFieldUFF(const char* ID, bool IsDefault=true) : OBForceField(ID, IsDefault)
      {
        ParseParamFile();
      }
      
    //! Destructor
    virtual ~OBForceFieldUFF();
      
    //! Assignment
    OBForceFieldUFF &operator = (OBForceFieldUFF &);
      
    //! Get the description for this force field
    const char* Description() 
    { 
      return "Universal Force Field.";
    }
      
    //! Get the unit in wich the energy is expressed
    std::string GetUnit() 
      { 
        return std::string("kJ/mol");  // Note that we convert from kcal/mol internally
      }

    //! Setup
    bool Setup(OBMol &mol);
      
    //! Returns total energy
    double Energy(bool gradients = true);
    //! Returns the bond stretching energy
    double E_Bond(bool gradients = true);
    //! Returns the angle bending energy
    double E_Angle(bool gradients = true);
    //! Returns the torsional energy
    double E_Torsion(bool gradients = true);
    //! Returns the out-of-plane (inversion) energy
    double E_OOP(bool gradients = true);
    //! Returns energy due to Van der Waals interactions
    double E_VDW(bool gradients = true);
    //! Returns energy due to electrostatic interactions
    double E_Electrostatic(bool gradients = true);
      
    //! Compare and print the numerical and analytical gradients
    bool ValidateGradients();

  }; // class OBForceFieldUFF

}// namespace OpenBabel

//! \file forcefieldUFF.h
//! \brief UFF force field

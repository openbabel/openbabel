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
    double bt; // bond order (e.g., 1.41 for amide)
      
    void Compute(bool gradients = true);
  };
  
  class OBFFAngleCalculationUFF : public OBFFCalculation
  {
  public:
    double ka, theta0, theta, delta;
    double c0, c1, c2;
    double zi, zk, rij, rjk, rik;
    double cosT0; // cos theta0
    int coord, n;

    void Compute(bool gradients = true);
  };
  
  class OBFFTorsionCalculationUFF : public OBFFCalculation
  {
  public:
    double V, tor, cosPhi0;
    int n;
    double tt; //torsiontype (i.e. b-c bond order)
      
    void Compute(bool gradients = true);
  };

  class OBFFOOPCalculationUFF : public OBFFCalculation
  {
  public:
    double koop, angle;
    double c0, c1, c2;
      
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
    bool SetTypes();
    //!  Fill OBFFXXXCalculation vectors
    bool SetupCalculations();
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
      OBForceFieldUFF(const char* ID, bool IsDefault=true) : OBForceField(ID, IsDefault)
      {
        _init = false;
      }
      
    //! Destructor
    virtual ~OBForceFieldUFF();
     
     //!Clone the current instance. May be desirable in multithreaded environments
     virtual OBForceFieldUFF* MakeNewInstance(){ return new OBForceFieldUFF(*this); }
 
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
    double E_Bond(bool gradients = true);
    //! \return the angle bending energy
    double E_Angle(bool gradients = true);
    //! \return the torsional energy
    double E_Torsion(bool gradients = true);
    //! \return the out-of-plane (inversion) energy
    double E_OOP(bool gradients = true);
    //! \return energy due to Van der Waals interactions
    double E_VDW(bool gradients = true);
    //! \return energy due to electrostatic interactions
    double E_Electrostatic(bool gradients = true);
      
    //! Compare and print the numerical and analytical gradients
    bool ValidateGradients();

  }; // class OBForceFieldUFF

}// namespace OpenBabel

//! \file forcefieldUFF.h
//! \brief UFF force field

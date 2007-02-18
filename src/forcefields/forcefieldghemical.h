/**********************************************************************
forcefieldghemical.h - Ghemical force field.
 
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

#include <openbabel/forcefield.h>
#include <openbabel/base.h>
#include <openbabel/mol.h>

namespace OpenBabel
{
  class OBFFBondCalculationGhemical : public OBFFCalculation
  {
    public:
      OBAtom *a, *b; // atoms of the bond
      double kb, r0, rab, delta;
      int bt; // bondtype
      
      double GetEnergy();
      vector3 GetGradient(OBAtom *atom);
  };
  
  class OBFFAngleCalculationGhemical : public OBFFCalculation
  {
    public:
      OBAtom *a, *b, *c; // atoms of the angle
      double ka, theta0, theta, delta;
      
      double GetEnergy();
      vector3 GetGradient(OBAtom *atom);
  };
  
  class OBFFTorsionCalculationGhemical : public OBFFCalculation
  {
    public:
      OBAtom *a, *b, *c, *d; // atoms of the torsion
      double V, s, n, tor;
      double k1, k2, k3;
      int tt; //torsiontype
      
      double GetEnergy();
      vector3 GetGradient(OBAtom *atom);
  };

  class OBFFVDWCalculationGhemical : public OBFFCalculation
  {
    public:
      OBAtom *a, *b; // atoms of the pair
      double ka, Ra, kb, Rb, kab, rab, sigma, sigma6, sigma12;

      double GetEnergy();
      vector3 GetGradient(OBAtom *atom);
  };

  class OBFFElectrostaticCalculationGhemical : public OBFFCalculation
  {
    public:
      OBAtom *a, *b; // atoms of the pair
      double qq, rab;
      
      double GetEnergy();
      vector3 GetGradient(OBAtom *atom);
  };


  // Class OBForceFieldGhemical
  // class introduction in forcefieldghemical.cpp
  class OBAPI OBForceFieldGhemical: public OBForceField
  {
    protected:
      //!  Parses the parameter file
      bool ParseParamFile();
      //!  Sets atomtypes to Ghemical types in _mol
      bool SetGhemicalTypes();
      //!  Sets partial charges to Ghemical charges in _mol
      bool SetGhemicalCharges();
      //! fill OBFFXXXCalculation vectors
      bool SetupCalculations();
      //! Same as OBForceField::GetParameter, but takes (bond/angle/torsion) type in account.
      OBFFParameter* GetParameterGhemical(int type, const char* a, const char* b, const char* c, const char* d, std::vector<OBFFParameter> &parameter);

      
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
      OBForceFieldGhemical(std::string ID, bool IsDefault=true) : OBForceField(ID, IsDefault)
      {
        ParseParamFile();
      }
      
      //! Destructor
      virtual ~OBForceFieldGhemical();
      
      //! Assignment
      OBForceFieldGhemical &operator = (OBForceFieldGhemical &);
      
      //! Get the description for this force field
      std::string Description() 
      { 
        return "Ghemical force field.";
      }
      
      //! Get the unit in wich the energy is expressed
      std::string GetUnit() 
      { 
        return std::string("kJ/mol"); 
      }

      //! Setup
      bool Setup(OBMol &mol);
      
      //! Returns total energy
      double Energy();
     //! Returns the bond stretching energy
      double E_Bond();
      //! Returns the angle bending energy
      double E_Angle();
      //! Returns the torsional energy
      double E_Torsion();
      //! Returns energy due to Van der Waals interactions
      double E_VDW();
      //! Returns energy due to electrostatic interactions
      double E_Electrostatic();
      
      //! Returns the negative gradient (force) on atom a
      vector3 GetGradient(OBAtom *a, int terms = OBFF_ENERGY);
      //! Compare and print the numerical and analytical gradients
      bool ValidateGradients();

  }; // class OBForceFieldGhemical

}// namespace OpenBabel

//! \file forcefieldGhemical.h
//! \brief Ghemical force field


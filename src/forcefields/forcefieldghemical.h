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
      OBAtom *a, *b, *c; // atoms of the angle
      double ka, theta0, thetaabc, delta, delta2;
      int sbt; //strbndtype (SBTIJK)
      
      double GetEnergy();
  };


  // Class OBForceFieldMM2
  // class introduction in forcefield.cpp
  class OBAPI OBForceFieldGhemical: public OBForceField
  {
    protected:
      //! \return Parses the parameter file
      bool ParseParamFile();
      //! \return Sets atomtypes to Ghemical in _mol
      bool SetGhemicalTypes();
      //! fill OBFFXXXCalculation vectors
      bool SetupCalculations();
      
      OBFFParameter* GetParameterGhemical(int type, const char* a, const char* b, const char* c, const char* d, 
          std::vector<OBFFParameter> &parameter);

      std::vector<OBFFParameter> _ffbondparams; 
      std::vector<OBFFParameter> _ffangleparams; 
      std::vector<OBFFParameter> _fftorsionparams; 
      std::vector<OBFFParameter> _ffvdwparams;

      std::vector<vector3> forces;

      std::vector<OBFFBondCalculationGhemical>          _bondcalculations;
      std::vector<OBFFAngleCalculationGhemical>         _anglecalculations;
      std::vector<OBFFTorsionCalculationGhemical>       _torsioncalculations;
      std::vector<OBFFVDWCalculationGhemical>           _vdwcalculations;
      std::vector<OBFFElectrostaticCalculationGhemical> _electostaticcalculations;

    public:
      //! Setup
      bool Setup(OBMol &mol);
      //! Constructor
      OBForceFieldGhemical(std::string ID, bool IsDefault=true) : OBForceField(ID, IsDefault)
      {
        ParseParamFile();
      }
      
      std::string Description()
	{ return "Ghemical force field.";};
      
      std::string GetUnit() { return std::string("kcal/mol"); }


      //! Destructor
      virtual ~OBForceFieldGhemical();
      //! Assignment
      OBForceFieldGhemical &operator = (OBForceFieldGhemical &);
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

      vector3 GetGradient(OBAtom *a);
      bool ValidateGradients();

  }; // class OBForceFieldGhemical

}// namespace OpenBabel

//! \file forcefieldGhemical.h
//! \brief Ghemical force field


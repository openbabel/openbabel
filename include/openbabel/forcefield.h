/**********************************************************************
forcefield.h - Handle OBForceField class.
 
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

#ifndef OB_FORCEFIELD_H
#define OB_FORCEFIELD_H

#include <vector>
#include <string>
#include <map>

#include <list>
#include <set>
#include <openbabel/base.h>
#include <openbabel/mol.h>
#include <openbabel/pluginiter.h>

namespace OpenBabel
{
  #define OBFF_LOGLVL_NONE	0 // no output
  #define OBFF_LOGLVL_LOW	1 // SteepestDescent progress... (no output from Energy())
  #define OBFF_LOGLVL_MEDIUM	2 // individual energy terms
  #define OBFF_LOGLVL_HIGH	3 // individual calculations and parameters

  #define IF_OBFF_LOGLVL_LOW    if(loglvl >= OBFF_LOGLVL_LOW)
  #define IF_OBFF_LOGLVL_MEDIUM if(loglvl >= OBFF_LOGLVL_MEDIUM)
  #define IF_OBFF_LOGLVL_HIGH   if(loglvl >= OBFF_LOGLVL_HIGH)

  class OBFFParameter {
    public:
      int         a, b, c, d; //! used to store integer atom types
      std::string _a, _b, _c, _d; //! used to store ascii atom types

      int       ipar1, ipar2, ipar3, ipar4, ipar5;
      double    dpar1, dpar2, dpar3, dpar4, dpar5;

      OBFFParameter& operator=(const OBFFParameter &ai) {
        if (this != &ai) {
          a = ai.a;
	  b = ai.b;
	  c = ai.c;
	  d = ai.d;
	  _a = ai._a;
	  _b = ai._b;
	  _c = ai._c;
	  _d = ai._d;
          ipar1 = ai.ipar1;
          ipar2 = ai.ipar2;
          ipar3 = ai.ipar3;
          ipar4 = ai.ipar4;
          ipar5 = ai.ipar5;
	  dpar1 = ai.dpar1;
	  dpar2 = ai.dpar2;
	  dpar3 = ai.dpar3;
	  dpar4 = ai.dpar4;
	  dpar5 = ai.dpar5;
        }
        
	return *this;
      }

      void clear () {
        a = 0;
	b = 0;
	c = 0;
	d = 0;
        ipar1 = 0;
        ipar2 = 0;
        ipar3 = 0;
        ipar4 = 0;
        ipar5 = 0;
	dpar1 = 0.0f;
	dpar2 = 0.0f;
	dpar3 = 0.0f;
	dpar4 = 0.0f;
	dpar5 = 0.0f;
      }
  };
  
  class OBFFCalculation
  {
    public:
      OBFFCalculation() 
      {
      }
      
      ~OBFFCalculation()
      {
      }

      virtual double Result() { return 0.0f; }
  };

  // Class OBForceField
  // class introduction in forcefield.cpp
  class OBAPI OBForceField
    {
  
    MAKE_PLUGIN(OBForceField)

    public:
      /// Required short description of the force field type.
      virtual std::string Description()=0;
      /// see Energy()
      double StoredEnergy;
      bool EnergyCalculated;
      void SetEnergyCalculated()   { EnergyCalculated = true; } 
      void UnsetEnergyCalculated() { EnergyCalculated = false; } //was bool
      bool IsSetEnergyCalculated() { return EnergyCalculated; } 

      /// Returns a pointer to a forcefield (the default if ID is empty), or NULL if not available
      static OBForceField* FindForceField(const std::string& ID){ return Iter().FindType(ID);}
      
      static OBForceField* FindForceField(const char *ID)
      {
        std::string ffname(ID);
	return FindForceField(ffname);
      }

    protected:
      //! Get the correct OBFFParameter from a OBFFParameter vector.
      //!
      //!  vector<OBFFParameter> parameters;
      //!  
      //!      this vector is filled with entries (as OBFFParameter) from 
      //!      a parameter file.
      //!
      //!  GetParameter(a, 0, 0, 0, parameters); 
      //!  
      //!      returns the first OBFFParameter from vector<OBFFParameter> 
      //!      parameters where: pa = a (pa = parameter.a)
      //!
      //!      use: vdw parameters, ...
      //!
      //!  GetParameter(a, b, 0, 0, parameters);
      //!
      //!      returns the first OBFFParameter from vector<OBFFParameter> 
      //!      parameters where: pa = a & pb = b      (ab)
      //!                    or: pa = b & pb = a      (ba)
      //!      
      //!      use: bond parameters, vdw parameters (pairs), ...
      //!
      //!  GetParameter(a, b, c, 0, parameters);
      //!
      //!      returns the first OBFFParameter from vector<OBFFParameter> 
      //!      parameters where: pa = a & pb = b & pc = c     (abc)
      //!                    or: pa = c & pb = b & pc = a     (cba)
      //!
      //!      use: angle parameters, ...
      //!
      //!  GetParameter(a, b, c, d, parameters);
      //!
      //!      returns the first OBFFParameter from vector<OBFFParameter> 
      //!      parameters where: pa = a & pb = b & pc = c & pd = d    (abcd)
      //!                    or: pa = d & pb = b & pc = c & pd = a    (dbca)
      //!                    or: pa = a & pb = c & pc = b & pd = d    (acbd)
      //!                    or: pa = d & pb = c & pc = b & pd = a    (dcba)
      //!
      //!      use: torsion parameters, ...
      //!
      OBFFParameter* GetParameter(int a, int b, int c, int d, std::vector<OBFFParameter> &parameter);
      OBFFParameter* GetParameter(const char* a, const char* b, const char* c, const char* d, std::vector<OBFFParameter> &parameter);
      
      OBFFParameter* GetParameter(const char* a, const char* b, const char* c, std::vector<OBFFParameter> &parameter)
      { return GetParameter(a, b, c, NULL, parameter); }
      OBFFParameter* GetParameter(const char* a, const char* b, std::vector<OBFFParameter> &parameter)
      { return GetParameter(a, b, NULL, NULL, parameter); }
      OBFFParameter* GetParameter(const char* a, std::vector<OBFFParameter> &parameter)
      { return GetParameter(a, NULL, NULL, NULL, parameter); }

      //! Calculate the energy ('only' called from GetEnergy, programs 
      //! should use GetEnergy(). This way, the Energy is not computed twice
      //! when the coordinates are unchanged. GetEnergy() evaluates the 
      //! value of 'bool EnergyCalculated' 
      //!   
      //!     true: return stored energy
      //!     false: call Energy, set EnergyCalculated = true and 
      //!            return newly calculated
      //!
      //! Each functions wich changes the coordinates of the atoms should 
      //! set EnergyCalculated to false. (SteepestDescent(), ...)
      //!
      //! If a function frequently changes coordinates it may be better
      //! to use Energy() instead of using:
      //!     
      //!     UnsetEnergyCalculated();
      //!     energy = GetEnergy();
      //!
      //! for each change in coordinates.
      virtual double Energy() { return 0.0f; }
 
      //! Get index for vector<OBFFParameter> ...
      int GetParameterIdx(int a, int b, int c, int d, std::vector<OBFFParameter> &parameter);
      //! Calculate the potential energy function derivative numerically with repect 
      //! to the coordinates of atom with index a (this vector is the gradient)
      vector3 NumericalDerivative(int a);

      int get_nbr (OBAtom* atom, int level);
      
      OBMol _mol;

      // ofstream for logfile
      std::ostream* logos;
      int loglvl;

    public:
      // see Energy()
      double GetEnergy()
      {
        if (!IsSetEnergyCalculated()) { 
	  StoredEnergy = Energy();
	  SetEnergyCalculated();
	}

	return StoredEnergy;
      }
     
      virtual std::string GetUnit() { return std::string("au"); }
      virtual double E_Bond() { return 0.0f; }
      virtual double E_Angle() { return 0.0f; }
      virtual double E_StrBnd() { return 0.0f; }
      virtual double E_Torsion() { return 0.0f; }
      virtual double E_OOP() { return 0.0f; }
      virtual double E_VDW() { return 0.0f; }
      virtual double E_Electrostatic() { return 0.0f; }
      virtual bool Setup(OBMol &mol) { return false; }
      bool SetLogFile(std::ostream *pos);
      bool SetLogLevel(int level);
      virtual bool Validate() { return false; }
 

      //! Update coordinates after steepest descent, ...
      void UpdateCoordinates(OBMol &mol);
      //! Generate coordinates for the molecule.
      void GenerateCoordinates();
      //! Perform steepest descent optimalization
      void SteepestDescent(int steps);
      //! Perform conjugate gradients optimalization
      void ConjugateGradients(int steps);
 
   }; // class OBForceField

}// namespace OpenBabel

#endif   // OB_FORCEFIELD_H

//! \file forcefield.h
//! \brief Handle forcefields

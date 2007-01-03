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
  

  // Class OBForceField
  // class introduction in forcefield.cpp
  class OBAPI OBForceField
    {
  
    MAKE_PLUGIN(OBForceField)

    public:
      /// Required short description of the force field type.
      virtual std::string Description()=0;
      /// Returns a pointer to a fingerprint (the default if ID is empty), or NULL if not available
      static OBForceField* FindForceField(const std::string& ID){ return Iter().FindType(ID);}
      
      static OBForceField* FindForceField(const char *ID)
      {
        std::string ffname(ID);
	return FindForceField(ffname);
      }

    protected:
      //! Get the correct OBFFParameter from a OBFFParameter vector
      OBFFParameter* GetParameter(int a, int b, int c, int d, std::vector<OBFFParameter> &parameter);
      OBFFParameter* GetParameter(const char* a, const char* b, const char* c, const char* d, std::vector<OBFFParameter> &parameter);
      OBFFParameter* GetParameter(const char* a, const char* b, const char* c, std::vector<OBFFParameter> &parameter)
      {
        return GetParameter(a, b, c, NULL, parameter);
      }
      OBFFParameter* GetParameter(const char* a, const char* b, std::vector<OBFFParameter> &parameter)
      {
        return GetParameter(a, b, NULL, NULL, parameter);
      }
      OBFFParameter* GetParameter(const char* a, std::vector<OBFFParameter> &parameter)
      {
        return GetParameter(a, NULL, NULL, NULL, parameter);
      }
      //! Get index for std::vector<OBFFParameter> ...
      int GetParameterIdx(int a, int b, int c, int d, std::vector<OBFFParameter> &parameter);
      //! Calculate angle between point d and plane a-b-c (used for out-of-plane-bending)
      double PointPlaneAngle(const vector3 &a, const vector3 &b, const vector3 &c, const vector3 &d);
      //! Get largest absolute value of "std::vector<vector3> forces" (used for steepest descent)
      //double GetFmax(OBMol &mol);
      vector3 NumericalDerivative(int a);
      
      OBMol _mol;

    public:
      //! Returns total energy
      virtual double Energy()
        { std::cerr << "Not a valid force field"; return false;}
      
      virtual double E_Bond() { return 0.0f; }
      virtual double E_Angle() { return 0.0f; }
      virtual double E_StrBnd() { return 0.0f; }
      virtual double E_Torsion() { return 0.0f; }
      virtual double E_OOP() { return 0.0f; }
      virtual double E_VDW() { return 0.0f; }
      virtual double E_Electrostatic() { return 0.0f; }
      virtual bool Setup(OBMol &mol) { return false; }
 

      //! Update coordinates after steepest descent, ...
      void UpdateCoordinates(OBMol &mol);
      //! Perform steepest descent optimalization
      void SteepestDescent(int steps);
 
   }; // class OBForceField

}// namespace OpenBabel

#endif   // OB_FORCEFIELD_H

//! \file forcefield.h
//! \brief Handle forcefields

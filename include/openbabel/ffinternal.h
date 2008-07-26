/**********************************************************************
ffinternal.h - Internal classes for the force field implementations.
 
Copyright (C) 2006-2008 by Tim Vandermeersch
 
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

#ifndef OB_FFINTERNAL_H
#define OB_FFINTERNAL_H

#include <vector>
#include <string>
#include <map>

#include <list>
#include <set>
#include <openbabel/babelconfig.h>
#include <openbabel/base.h>
#include <openbabel/mol.h>
#include <openbabel/plugin.h>
#include <openbabel/grid.h>
#include <openbabel/griddata.h>
#include <float.h>

namespace OpenBabel
{
  //! \class OBFFParameter forcefield.h <openbabel/forcefield.h>
  //! \brief Internal class for OBForceField to hold forcefield parameters
  class OBFPRT OBFFParameter {
  public:
    //! Used to store integer atom types
    int         a, b, c, d;
    //! used to store string atom types
    std::string _a, _b, _c, _d; 
    //! Used to store integer type parameters (bondtypes, multiplicity, ...)
    std::vector<int>    _ipar;
    //! Used to store double type parameters (force constants, bond lengths, angles, ...)
    std::vector<double> _dpar;

    //! Assignment 
    OBFFParameter& operator=(const OBFFParameter &ai) 
      {
        if (this != &ai) {
          a = ai.a;
          b = ai.b;
          c = ai.c;
          d = ai.d;
          _a = ai._a;
          _b = ai._b;
          _c = ai._c;
          _d = ai._d;
          _ipar = ai._ipar;
          _dpar = ai._dpar;
        }
        
        return *this;
      }

    //! Reset the atom types and set all parameters to zero
    void clear () 
    {
      a = b = c = d = 0;
      _ipar.clear();
      _dpar.clear();
    }
  }; // class OBFFParameter
  
  // specific class introductions in forcefieldYYYY.cpp (for YYYY calculations)

  //! \class OBFFCalculation2 forcefield.h <openbabel/forcefield.h>
  //! \brief Internal class for OBForceField to hold energy and gradient calculations on specific force fields
  class OBFPRT OBFFCalculation2
  {
  public:
    //! Used to store the atoms for this OBFFCalculation
    OBAtom *a, *b;
    //! Destructor
    virtual ~OBFFCalculation2() 
    {
    }    
  };
 
  //! \class OBFFCalculation3 forcefield.h <openbabel/forcefield.h>
  //! \brief Internal class for OBForceField to hold energy and gradient calculations on specific force fields
  class OBFPRT OBFFCalculation3: public OBFFCalculation2
  {
  public:
    //! Used to store the atoms for this OBFFCalculation
    OBAtom *c;
    //! Destructor
    virtual ~OBFFCalculation3() 
    {
    }    
  };

  //! \class OBFFCalculation4 forcefield.h <openbabel/forcefield.h>
  //! \brief Internal class for OBForceField to hold energy and gradient calculations on specific force fields
  class OBFPRT OBFFCalculation4: public OBFFCalculation3
  {
  public:
    //! Used to store the atoms for this OBFFCalculation
    OBAtom *d;
    //! Destructor
    virtual ~OBFFCalculation4() 
    {
    }    
  };

  //! \class OBFFConstraint forcefield.h <openbabel/forcefield.h>
  //! \brief Internal class for OBForceField to hold constraints
  //! \since version 2.2
  class OBFPRT OBFFConstraint
  {
  public:
    //! Used to store the contraint energy for this OBFFConstraint
    double factor, constraint_value;
    double rab0, rbc0;
    //! Used to store the contraint type for this OBFFConstraint
    int type, ia, ib, ic, id;
    //! Used to store the atoms for this OBFFCostraint
    OBAtom *a, *b, *c, *d;
    //! Used to store the gradients for this OBFFCalculation
    Eigen::Vector3d Fa, Fb, Fc, Fd;

    //! Constructor
    OBFFConstraint() 
      {
        a = b = c = d = NULL;
        ia = ib = ic = id = 0;
        constraint_value = 0.0;
        factor = 0.0;
      }
    //! Destructor
    ~OBFFConstraint()
      {
      }
      
    Eigen::Vector3d GetGradient(int a) 
    {
      if (a == ia)
        return Fa;
      else if (a == ib)
        return Fb;
      else if (a == ic)
        return Fc;
      else if (a == id)
        return Fd;
      else 
        return  VZero;
    }
  };

}// namespace OpenBabel

#endif   // OB_FORCEFIELD_H

//! \file forcefield.h
//! \brief Handle forcefields

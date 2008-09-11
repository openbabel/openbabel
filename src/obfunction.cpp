/**********************************************************************
obfunction.cpp - 

Copyright (C) 2008 by Tim Vandermeersch
 
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
#include <openbabel/babelconfig.h>

#include <openbabel/obfunction.h>

using namespace std;

namespace OpenBabel
{
  
  //  
  //         f(1) - f(0)
  // f'(0) = -----------      f(1) = f(0+h)
  //              h
  //
  Eigen::Vector3d OBFunction::NumericalDerivative(unsigned int idx) // start at 0!!
  {
    Eigen::Vector3d va, grad;
    double e_orig, e_plus_delta, delta, dx, dy, dz;

    delta = 1.0e-5;
    va = GetPositions()[idx];

    e_orig = Eval(false);
    
    // X direction
    GetPositions()[idx] = Eigen::Vector3d(va.x() + delta, va.y(), va.z());
    e_plus_delta = Eval(false);
    dx = (e_plus_delta - e_orig) / delta;
    
    // Y direction
    GetPositions()[idx] = Eigen::Vector3d(va.x(), va.y() + delta, va.z());
    e_plus_delta = Eval(false);
    dy = (e_plus_delta - e_orig) / delta;
    
    // Z direction
    GetPositions()[idx] = Eigen::Vector3d(va.x(), va.y(), va.z() + delta);
    e_plus_delta = Eval(false);
    dz = (e_plus_delta - e_orig) / delta;

    // reset coordinates to original
    GetPositions()[idx] = Eigen::Vector3d(va.x(), va.y(), va.z());

    grad = Eigen::Vector3d(-dx, -dy, -dz);
    return (grad);
  }
  
  //  
  //         f(2) - 2f(1) + f(0)
  // f'(0) = -------------------      f(1) = f(0+h)
  //                 h^2              f(1) = f(0+2h)
  //
  Eigen::Vector3d OBFunction::NumericalSecondDerivative(unsigned int idx)
  {
    Eigen::Vector3d va, grad;
    double e_0, e_1, e_2, delta, dx, dy, dz;

    delta = 1.0e-5;

    va = GetPositions()[idx];

    // calculate f(0)
    e_0 = Eval(false);
    
    // 
    // X direction
    //
    
    // calculate f(1)
    GetPositions()[idx] = Eigen::Vector3d(va.x() + delta, va.y(), va.z());
    e_1 = Eval(false);
    
    // calculate f(2)
    GetPositions()[idx] = Eigen::Vector3d(va.x() + 2 * delta, va.y(), va.z());
    e_2 = Eval(false);
    
    dx = (e_2 - 2 * e_1 + e_0) / (delta * delta);
    
    // 
    // Y direction
    //
    
    // calculate f(1)
    GetPositions()[idx] = Eigen::Vector3d(va.x(), va.y() + delta, va.z());
    e_1 = Eval(false);
    
    // calculate f(2)
    GetPositions()[idx] = Eigen::Vector3d(va.x(), va.y() + 2 * delta, va.z());
    e_2 = Eval(false);
    
    dy = (e_2 - 2 * e_1 + e_0) / (delta * delta);

    // 
    // Z direction
    //
    
    // calculate f(1)
    GetPositions()[idx] = Eigen::Vector3d(va.x(), va.y(), va.z() + delta);
    e_1 = Eval(false);
    
    // calculate f(2)
    GetPositions()[idx] = Eigen::Vector3d(va.x(), va.y(), va.z() + 2 * delta);
    e_2 = Eval(false);
    
    dz = (e_2 - 2 * e_1 + e_0) / (delta * delta);

    // reset coordinates to original
    GetPositions()[idx] = Eigen::Vector3d(va.x(), va.y(), va.z());

    grad = Eigen::Vector3d(-dx, -dy, -dz);
    return (grad);
  }
  
} // end namespace OpenBabel


//! @file obfunction.cpp
//! @brief Handle OBFunction class

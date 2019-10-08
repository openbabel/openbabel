/**********************************************************************
qtpie.h - A OBChargeModel to handle QTPIE charges

Copyright (C) 2010 by Jiahao Chen <jiahao@mit.edu>

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

#ifndef __QTPIE__H__
#define __QTPIE__H__

#include <openbabel/chargemodel.h>
#include <openbabel/mol.h>

#include <math.h>

#ifdef HAVE_EIGEN

#include <Eigen/LU>
#include <Eigen/SVD>

#define pi 3.1415926
#define sqr(x) ((x)*(x))
#define cube(x) ((x)*(x)*(x))

//// conversion factor from electron volt to Hartree
const double eV = 3.67493245e-2;

/// conversion factor from Angstrom to bohr
const double Angstrom = 1./0.529177249;

/// Determines threshold value of Coulomb interaction beyond which the classical 1/R expression
/// will be used instead of an integral over the Gaussian orbitals
/// This is a pre-screening cutoff.
const double CoulombThreshold = 1.e-9;
const double OverlapThreshold = 1.e-9;

namespace OpenBabel
{

class QTPIECharges : public OBChargeModel
{
public:
  QTPIECharges(void) : OBChargeModel("fake ID", false){};
  QTPIECharges(const char* ID) : OBChargeModel(ID, false){};
  const char* Description(){ return "Assign QTPIE (charge transfer, polarization and equilibration) partial charges (Chen and Martinez, 2007)"; }

  /// \return whether partial charges were successfully assigned to this molecule
  bool ComputeCharges(OBMol &mol);

  double DipoleScalingFactor() { return 4.041; } // fit from regression

private:
  Eigen::Vector3d GetParameters(unsigned int Z, int Q);
  bool solver(Eigen::MatrixXd A, Eigen::VectorXd b, Eigen::VectorXd &x, const double NormThreshold = 1.e-6);
  double CoulombInt(double a, double b, double R);
  double OverlapInt(double a, double b, double R);

  Eigen::MatrixXd Hardness; ///The hardness matrix
  Eigen::VectorXd Electronegativity, Voltage, Charge;
  double ChemicalPotential;

  std::vector<Eigen::Vector3d> _parameters;
  void ParseParamFile();
};
}; //namespace OpenBabel
#endif //HAVE_EIGEN
#endif //__QTPIE_H__

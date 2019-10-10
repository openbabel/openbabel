/**********************************************************************
qeq.h - A OBChargeModel to handle QEq charges

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

#ifndef __QEQ__H_
#define __QEQ__H_

#include <openbabel/chargemodel.h>
#include <openbabel/mol.h>

#include <math.h>

#ifdef HAVE_EIGEN

#include <Eigen/LU>
#include <Eigen/SVD>

//// conversion factor from electron volt to Hartree
const double eV = 3.67493245e-2;

/// conversion factor from Angstrom to bohr
const double Angstrom = 1./0.529177249;

/// Determines threshold value of Coulomb interaction beyond which the classical 1/R expression
/// will be used instead of an integral over the Gaussian orbitals
/// This is a pre-screening cutoff.
const double CoulombThreshold = 1.e-9;

namespace OpenBabel
{

class QEqCharges : public OBChargeModel
{
public:
  QEqCharges(void) : OBChargeModel("fake ID", false){};
  QEqCharges(const char* ID) : OBChargeModel(ID, false){};
  const char* Description(){ return "Assign QEq (charge equilibration) partial charges (Rappe and Goddard, 1991)"; }

  /// \return whether partial charges were successfully assigned to this molecule
  bool ComputeCharges(OBMol &mol);

  double DipoleScalingFactor() { return 1.959; } // fit from regression

private:
  Eigen::Vector3d GetParameters(unsigned int Z, int Q);
  bool solver(Eigen::MatrixXd A, Eigen::VectorXd b, Eigen::VectorXd &x, const double NormThreshold = 1.e-6);
  double CoulombInt(double a, double b, double R);

  Eigen::MatrixXd Hardness; ///The hardness matrix
  Eigen::VectorXd Electronegativity, Voltage, Charge;
  double ChemicalPotential;

  std::vector<Eigen::Vector3d> _parameters;
  void ParseParamFile();
};

}; //namespace OpenBabel
#endif //HAVE_EIGEN
#endif //__QEQ_H__

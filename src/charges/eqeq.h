/**********************************************************************
eqeq.h - A OBChargeModel to handle EQEq charges

Copyright (C) 2013 by Patrick Fuller <patrickfuller@gmail.com>

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

#ifndef __EQEQ__H_
#define __EQEQ__H_

#include <openbabel/chargemodel.h>
#include <openbabel/mol.h>

#include <math.h>

#ifdef HAVE_EIGEN3

#include <Eigen/LU>
#include <Eigen/SVD>

#define TABLE_OF_ELEMENTS_SIZE 84 // Number of atoms in data/eqeqIonizations.txt
#define PI 3.1415926535897932384626433832795 // 32 digits of pi

static bool _paramFileLoaded = false; // Flipped after parameter file is loaded

// The following constants are functionalized in EQeq's original implementation,
// but rarely ever changed in practice.
const double minCellLength = 150.0; // Minimum dimension of supercell used in calculation [A]
const double k = 14.4; // Vacuum permittivity (1/(4 * pi * eps_0)) [A * eV]
const double lambda = 1.2; // Coulomb scaling parameter
const double hi_0 = -2.0; // Electron affinity of hydrogen used in EQeq paper
const double eta = 50; // Ewald splitting parameter

namespace OpenBabel
{

class EQEqCharges : public OBChargeModel
{
public:
  EQEqCharges(void) : OBChargeModel("fake ID", false){};
  EQEqCharges(const char* ID) : OBChargeModel(ID, false){};
  const char* Description(){ return "Assign EQEq (charge equilibration) partial charges."; }

  /// \return whether partial charges were successfully assigned to this molecule
  bool ComputeCharges(OBMol &mol);

private:
  int _chargeCenter[TABLE_OF_ELEMENTS_SIZE + 1]; // Common charge of metallic ions
  double _ionizations[TABLE_OF_ELEMENTS_SIZE + 1][9]; // Electron affinity + 8x ionizations of elements
  bool ParseParamFile();
  double GetNonperiodicJij(double J_i, double J_j, double R_ij, bool isSameAtom);
  double GetPeriodicEwaldJij(double J_i, double J_j, vector3 dx, bool isSameAtom, matrix3x3 unitcell, matrix3x3 fourier, double cellVolume, int numNeighbors[]);
};

}; //namespace OpenBabel
#endif //HAVE_EIGEN3
#endif //__EQEQ_H__

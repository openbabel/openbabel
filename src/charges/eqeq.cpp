/**********************************************************************
eqeq.cpp - A OBChargeModel to handle EQEq charges

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

#ifdef HAVE_EIGEN3

#include "eqeq.h"
#include <openbabel/locale.h>
#include <openbabel/oberror.h>
#include <openbabel/atom.h>
#include <openbabel/generic.h>



using namespace std;

#if defined(_MSC_VER) && _MSC_VER < 1800
// Older MSVC doesn't have error function erfc, use local implementation
#include <openbabel/math/erf.h>
using temperf::erfc;
#endif

namespace OpenBabel
{
  /*! \class EQEqCharges eqeq.h "eqeq.h"

    \brief Assigns partial charges according to the extended charge equilibration (EQEq) model, according to Wilmer et al., J. Phys. Chem. Let. 2012.
  */
  EQEqCharges theEQEqCharges("eqeq"); // Global instance

  //! Loads a file containing nth ionization energies and metal ionic charges
  // TODO Move all relevant data into OBElementTable and eliminate this
  // function / data file.
  bool EQEqCharges::ParseParamFile()
  {
    int atomicNumber, i;
    vector<string> vs;
    char buffer[BUFF_SIZE];

    // open data/eqeqIonizations.txt
    ifstream ifs;
    if (OpenDatafile(ifs, "eqeqIonizations.txt").length() == 0)
    {
      obErrorLog.ThrowError(__FUNCTION__, "Cannot open eqeqIonizations.txt", obError);
      return false;
    }

    // Set the locale for number parsing to avoid locale issues: PR#1785463
    obLocale.SetLocale();
    while (ifs.getline(buffer, BUFF_SIZE))
    {
      if (buffer[0] == '#')
        continue;

      tokenize(vs, buffer);
      if (vs.size() != 12)
      {
        obErrorLog.ThrowError(__FUNCTION__, "Format error in eqeqIonizations.txt. Each data row should have exactly 12 elements.", obError);
        return false;
      }

      // Reads in a line of parameters. The format is:
      //   Atomic Number    Element    Ion Charge    Affinity    Ionizations (x8)
      atomicNumber = atoi(vs[0].c_str());
      _chargeCenter[atomicNumber] = atoi(vs[2].c_str());
      for (i = 0; i < 9; i ++)
        _ionizations[atomicNumber][i] = atof(vs[i + 3].c_str());

      // The electron affinity of hydrogen is a custom-fit parameter
      _ionizations[1][0] = hi_0;
    }
    return true;
  }

  //! \return whether partial charges were successfully assigned to this molecule
  bool EQEqCharges::ComputeCharges(OBMol &mol)
  {
    int i, j, a, c, N = mol.NumAtoms();
    double cellVolume;
    Eigen::VectorXf chi(N), J(N), b(N), x(N);
    Eigen::MatrixXf J_ij(N, N), A(N, N);
    OBUnitCell *obuc;
    matrix3x3 unitcell, fourier;
    vector3 dx;
    int numNeighbors[3];
    OBAtom *atom;

    // If parameters have not yet been loaded, do that
    if (!_paramFileLoaded)
    {
      if (ParseParamFile())
      {
        _paramFileLoaded = true;
      } else
      {
        return false;
      }
    }

    // Calculate atomic properties based around their ionic charge
    for (i = 0; i < N; i++)
    {
      atom = mol.GetAtom(i + 1);
      a = atom->GetAtomicNum();
      c = _chargeCenter[a];

      // Fail if ionization data is missing for any atom in the molecule
      if (_ionizations[a][c + 1] == -1 || _ionizations[a][c] == -1 || a > TABLE_OF_ELEMENTS_SIZE)
      {
        obErrorLog.ThrowError(__FUNCTION__, "Insufficient ionization data for atoms in the given molecule. Update `data/eqeqIonizations.txt` with missing information and re-run this function.", obError);
        return false;
      }

      J(i) = _ionizations[a][c + 1] - _ionizations[a][c];
      chi(i) = 0.5 * (_ionizations[a][c + 1] + _ionizations[a][c]) - (a == 1? 0 : c * J(i));
    }

    // If a unit cell is defined, use the periodic Ewald calculation
    if (mol.HasData(OBGenericDataType::UnitCell))
    {
      // Get unit cell and calculate its Fourier transform + volume
      obuc = (OBUnitCell *) mol.GetData(OBGenericDataType::UnitCell);
      unitcell = obuc->GetCellMatrix();
      fourier = (2 * PI * unitcell.inverse()).transpose();
      cellVolume = obuc->GetCellVolume();

      // Get the number of radial unit cells to use in x, y, and z
      numNeighbors[0] = int(ceil(minCellLength / (2.0 * (obuc->GetA())))) - 1;
      numNeighbors[1] = int(ceil(minCellLength / (2.0 * (obuc->GetB())))) - 1;
      numNeighbors[2] = int(ceil(minCellLength / (2.0 * (obuc->GetC())))) - 1;

      for (i = 0; i < N; i++)
      {
        atom = mol.GetAtom(i + 1);
        for (j = 0; j < N; j++)
        {
          dx = atom->GetVector() - (mol.GetAtom(j + 1))->GetVector();
          J_ij(i, j) = GetPeriodicEwaldJij(J(i), J(j), dx, (i == j), unitcell, fourier, cellVolume, numNeighbors);
        }
      }
    // If no unit cell, use the simplified nonperiodic calculation
    } else
    {
      for (i = 0; i < N; i++)
      {
        atom = mol.GetAtom(i + 1);
        for (j = 0; j < N; j++)
        {
          J_ij(i, j) = GetNonperiodicJij(J(i), J(j), atom->GetDistance(j + 1), (i == j));
        }
        return false;
      }
    }

    // Formulate problem as A x = b, where x is the calculated partial charges
    // First equation is a simple overall balance: sum(Q) = 0
    A.row(0) = Eigen::VectorXf::Ones(N);
    b(0) = 0;

    // Remaining equations are based off of the fact that, at equilibrium, the
    // energy of the system changes equally for a change in any charge:
    //     dE/dQ_1 = dE/dQ_2 = ... = dE/dQ_N
    A.block(1, 0, N - 1, N) = J_ij.block(0, 0, N - 1, N) - J_ij.block(1, 0, N - 1, N);
    b.tail(N - 1) = chi.tail(N - 1) - chi.head(N - 1);

    // The solution is a list of charges in the system
    x = A.colPivHouseholderQr().solve(b);

    // Now we are done calculating, pass all this back to OpenBabel molecule
    mol.SetPartialChargesPerceived();
    OBPairData *dp = new OBPairData;
    dp->SetAttribute("PartialCharges");
    dp->SetValue("EQEq");
    dp->SetOrigin(perceived);
    mol.SetData(dp);

    m_partialCharges.clear();
    m_partialCharges.reserve(N);
    m_formalCharges.clear();
    m_formalCharges.reserve(N);

    for (i = 0; i < N; i ++)
    {
      atom = mol.GetAtom(i + 1);
      atom->SetPartialCharge(x(i));
      m_partialCharges.push_back(x(i));
      m_formalCharges.push_back(atom->GetFormalCharge());
    }

    obErrorLog.ThrowError(__FUNCTION__, "EQEq charges successfully assigned.", obInfo);
    return true;
  }

  //! Calculates a lumped charge coefficient for simplified nonperiodic systems
  double EQEqCharges::GetNonperiodicJij(double J_i, double J_j, double R_ij, bool isSameAtom)
  {
    if (isSameAtom)
      return J_i;

    double a_ij, E_0ij;
    a_ij = sqrt(J_i * J_j) / k;
    E_0ij = exp(-a_ij * a_ij * R_ij * R_ij) * (2 * a_ij - a_ij * a_ij * R_ij - 1 / R_ij); // Orbital overlap term

    return lambda * (k / 2.0) * (1.0 / R_ij + E_0ij);
  }

  //! Calculates a lumped charge coefficient that incorporates neighboring unit cells
  double EQEqCharges::GetPeriodicEwaldJij(double J_i, double J_j, vector3 dx, bool isSameAtom, matrix3x3 unitcell, matrix3x3 fourier, double cellVolume, int numNeighbors[])
  {
    int u, v, w;
    double R_ij, a_ij, h_2, orbital = 0.0, alpha = 0.0, beta = 0.0;
    vector3 uvw, rlv;

    // Calculate unit-cell-independent data
    a_ij = sqrt(J_i * J_j) / k;
    uvw = vector3();

    for (u = -numNeighbors[0]; u <= numNeighbors[0]; u++)
    {
      for (v = -numNeighbors[1]; v <= numNeighbors[1]; v++)
      {
        for (w = -numNeighbors[2]; w <= numNeighbors[2]; w++)
        {
          // If we're iterating over same atom + center unit cell, skip
          if (isSameAtom && u == 0 && v == 0 && w == 0)
            continue;

          // Orbital overlap term summation
          uvw.Set(u, v, w);
          R_ij = (dx + unitcell * uvw).length();
          orbital += exp(-a_ij * a_ij * R_ij * R_ij) * (2 * a_ij - a_ij * a_ij * R_ij - 1 / R_ij);

          // Real-space Coulomb component summation
          alpha += erfc(R_ij / eta) / R_ij;

          // Skip Fourier math below if same unit cell
          if (u == 0 && v == 0 && w == 0)
            continue;

          // K-space component summation (from reciprocal lattice vector)
          rlv = fourier * uvw;
          h_2 = rlv.length_2();
          beta += cos(dot(rlv, dx)) * exp(-0.25 * h_2 * eta * eta) / h_2;
        }
      }
    }
    // Fourier normalization
    beta *= 4.0 * PI / cellVolume;

    // Combine and return final J_ij, with modifier if the atoms are the same
    return lambda * k / 2.0 * (alpha + beta + orbital) + (isSameAtom? J_i - lambda * k / (eta * sqrt(PI)) : 0);
  }

}//namespace

#endif //HAVE_EIGEN3

//! \file eqeq.cpp
//! \brief Assign EQEq partial charges.

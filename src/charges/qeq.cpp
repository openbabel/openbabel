/**********************************************************************
qeq.cpp - A OBChargeModel to handle QEq charges

Copyright (C) 2010 by Jiahao Chen <jiahao@mit.edu>

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

#if (defined(HAVE_EIGEN2))

#include "qeq.h"
#include <openbabel/locale.h>

using namespace std;

#ifdef _MSC_VER
// MSVC doesn't have error function erf, use local implementation
#include <openbabel/math/erf.h>
using temperf::erf;
#endif

namespace OpenBabel
{
  /////////////////////////////////////////////////////////////////
  QEqCharges theQEqCharges("qeq"); //Global instance

  /////////////////////////////////////////////////////////////////

  void QEqCharges::ParseParamFile()
  {
    vector<string> vs;
    char buffer[BUFF_SIZE];

    // open data/qeq.txt
    ifstream ifs;
    if (OpenDatafile(ifs, "qeq.txt").length() == 0) {
      obErrorLog.ThrowError(__FUNCTION__, "Cannot open qeq.txt", obError);
      return;
    }

    // Set the locale for number parsing to avoid locale issues: PR#1785463
    obLocale.SetLocale();
    Vector3d P;

    while (ifs.getline(buffer, BUFF_SIZE)) {
      if (buffer[0] == '#')
        continue;

      tokenize(vs, buffer);
      if (vs.size() < 6)
        continue;

      // Element, n, Xi, Hardness, Slater, Gaussian
      P << atof(vs[2].c_str())*eV, atof(vs[3].c_str())*eV*2, atof(vs[5].c_str());
      _parameters.push_back(P);
    }
  }

  Vector3d QEqCharges::GetParameters(unsigned int Z, int Q)
  {
    //Returns a triple of numbers: electronegativity (in eV), hardness (in eV), and Gaussian exponent (in bohr^-2)

    Vector3d P;
    //For now, completely ignore the formal charge
    if (_parameters.size() == 0)
      ParseParamFile();

    if (Z > 0 && Z < _parameters.size()-1)
      return _parameters[Z - 1];

    P<<0., 1.e10, 1.e10; //Magic value
    return P;
  }

  bool QEqCharges::ComputeCharges(OBMol &mol)
  {

    ///////////////////////////////////////////////////////////////////////////////
    //Some OpenBabel bookkeeping that I copied from the Gasteiger scheme
    mol.SetPartialChargesPerceived();

    OBPairData *dp = new OBPairData;
    dp->SetAttribute("PartialCharges");
    dp->SetValue("QEq");
    dp->SetOrigin(perceived);
    mol.SetData(dp);


    ///////////////////////////////////////////////////////////////////////////////
    //Read in atomic information from OpenBabel molecule and parameterize

    //Read in total number of atoms
    int i, N = mol.NumAtoms();

    Hardness = MatrixXd::Zero(N+1, N+1);
    Voltage = VectorXd::Zero(N+1);
    Electronegativity = VectorXd::Zero(N);
    VectorXd BasisSet = VectorXd::Zero(N);

    Vector3d Parameters;

    FOR_ATOMS_OF_MOL(atom, mol)
      {
       	Parameters = GetParameters(atom->GetAtomicNum(), atom->GetFormalCharge());
        i = atom->GetIdx() - 1;

        if (Parameters[0] == 0.)
          {
            stringstream msg;
            msg << "Some QEq Parameters not found!" << endl
                << "Parameters not found for atom no. " << i+1 << endl
                << "Atom will be ignored in the charge computation.";
            obErrorLog.ThrowError(__FUNCTION__, msg.str(), obError);
          }

        Electronegativity[i] = Parameters[0];
        Hardness(i,i) = Parameters[1];
        BasisSet[i] = Parameters[2];
      }

    //Read in total charge of molecule
    double Z = (double)(mol.GetTotalCharge());

    ///////////////////////////////////////////////////////////////////////////////
    // Now populate integrals

    // Calculate integral pre-screening threshold
    double SmallestGaussianExponent = BasisSet.minCoeff();
    double CoulMaxDistance = 2 * sqrt(-log(CoulombThreshold) / SmallestGaussianExponent);

    int j;
    double R, Coulomb;
    FOR_ATOMS_OF_MOL(atom1, mol)
      {
        i = atom1->GetIdx() - 1;
        FOR_ATOMS_OF_MOL(atom2, mol)
          {
            j = atom2->GetIdx() - 1;
            if (i>j)
              {
                //For some reason, this did _not_ produce the expected pairwise distances
                //(2,1) evaluated to the same distance as (2,0) in H2O: bug?
                // - jiahao@mit.edu 2010-04-20 r3691
                ////R = atom1->GetDistance(atom2)*Angstrom;
                //
                R = atom1->GetDistance(j+1)*Angstrom;
                if (R<CoulMaxDistance)
                  Coulomb = CoulombInt(BasisSet[i], BasisSet[j], R);
                else
                  Coulomb = 1./R;
                Hardness(i,j) = Coulomb;
                Hardness(j,i) = Coulomb;
              }
          }
      }

    Hardness.block(N,0,1,N).setOnes();
    Hardness.block(0,N,N,1).setOnes();

    Voltage.segment(0,N) = Electronegativity;
    Voltage[N] = Z;

    ///////////////////////////////////////////////////////////////////////////////
    // Call linear algebra solver
    //
    bool status = solver(Hardness, Voltage, Charge);
    if (!status)
      {
        obErrorLog.ThrowError(__FUNCTION__, "Charges could not be computed.", obError);
        return false;
      }

    ChemicalPotential = Charge[N];

    //For debugging purposes only
#if 0
    cout << "Ready to solve QEq problem" << endl << endl;
    cout << "Hardness matrix:" << endl << Hardness << endl << endl;
    cout << "Voltage vector:" << endl << Voltage << endl << endl;
    cout << "Charge vector:" << endl << Charge << endl << endl;
#endif

    //Now we are done calculating, pass all this back to OpenBabel molecule
    m_partialCharges.clear();
    m_partialCharges.reserve(mol.NumAtoms());
    m_formalCharges.clear();
    m_formalCharges.reserve(mol.NumAtoms());
    FOR_ATOMS_OF_MOL(atom, mol)
      {
        atom->SetPartialCharge(Charge[atom->GetIdx()-1]);
        m_partialCharges.push_back(Charge[atom->GetIdx() - 1]);
        m_formalCharges.push_back(atom->GetFormalCharge());
      }
    obErrorLog.ThrowError(__FUNCTION__, "QEq charges successfully assigned.", obInfo);
    return true;
  }


  /// Calculates Coulomb integral
  double QEqCharges::CoulombInt(double a, double b, double R)
  {
    double p = sqrt(a * b / (a + b));
    return erf(p * R) / R;
  }

  /// Here's a wrapper around the Eigen solver routine
  bool QEqCharges::solver(MatrixXd A, VectorXd b, VectorXd &x, const double NormThreshold)
  {
    // using a LU factorization
    bool SolverOK = A.lu().solve(b, &x);
    //bool SolverOK = A.svd().solve(b, &x);

    VectorXd resid = A*x - b;
    double resnorm = resid.norm();
    if (IsNan(resnorm) || resnorm > NormThreshold || !SolverOK)
      {
        stringstream msg;
        msg << "Warning, LU solver failed." << endl;
        if (!SolverOK) msg << "Solver returned error." << endl;
        if (IsNan(resnorm)) msg << "NaNs were returned" << endl;
        if (resnorm > NormThreshold) msg << "Residual has norm " << resnorm
                                         << " which exceeds the recommended threshold of " << NormThreshold
                                         << endl;
        msg << "Proceeding with singular value decomposition.";

        obErrorLog.ThrowError(__FUNCTION__, msg.str(), obWarning);

        SolverOK = A.svd().solve(b, &x);
        resid = A*x - b;
        resnorm = resid.norm();

        if (IsNan(resnorm) || !SolverOK)
          {
            obErrorLog.ThrowError(__FUNCTION__, "SVD solver returned an error. Charges may not be reliable!", obError);
            return false;
          }
      }

    stringstream msg_resid;
    msg_resid << "The residual of the solution has norm " << resnorm;
    obErrorLog.ThrowError(__FUNCTION__, msg_resid.str(), obInfo);

    if (resnorm > NormThreshold) {
      stringstream msg_reswarn;
      msg_reswarn << "Warning, the norm of the residual is " << resnorm
                  << "which exceeds the recommended threshold of " << NormThreshold;
      obErrorLog.ThrowError(__FUNCTION__, msg_reswarn.str(), obWarning);
    }
    return true;
  }

}//namespace

#endif //HAVE_EIGEN2

//! \file qeq.cpp
//! \brief Assign QEq partial charges.

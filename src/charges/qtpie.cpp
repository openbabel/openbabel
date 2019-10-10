/**********************************************************************
qtpie.cpp - A OBChargeModel to handle QTPIE charges

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

#ifdef HAVE_EIGEN

#include "qtpie.h"
#include <openbabel/locale.h>
#include <openbabel/oberror.h>
#include <openbabel/generic.h>
#include <openbabel/obiter.h>
#include <openbabel/atom.h>

using namespace std;

#if defined(_MSC_VER) && _MSC_VER < 1800
// Older MSVC doesn't have error function erf, use local implementation
#include <openbabel/math/erf.h>
using temperf::erf;
#endif

namespace OpenBabel
{

  /*! \class QTPIECharges qtpie.h "qtpie.h"

    \brief Assigns partial charges according to the charge transfer with polarization current equilibration
    (QTPIE) model of Chen and Martinez, 2007; 2008.

    The QTPIE model solves for charges by minimizing an energy function of the form

    \f[

    E\left(\mathbf{q}\right)=\mathbf{q}\cdot\boldsymbol{\chi}+\frac{1}{2}\mathbf{q}\boldsymbol{\eta}\mathbf{q}

    \f]

    where \f$\mathbf q\f$ is the vector of atomic charges, \f$\boldsymbol \chi\f$ are atomic electronegativities
    and \f$\boldsymbol \eta \f$ is the matrix of Coulomb interactions with chemical hardnesses along the diagonal,
    such that the overall charge is \f$Q\f$, i.e.

    \f[

    \mathbf{q} \cdot \mathbf{1} = \sum_{i=1}^N q_i = Q

    \f]

	in the current implementation, we assume \f$Q = 0\f$ always.

    The off-diagonal Coulomb interactions are screened using the following integral

    \f[

    \eta_{ij}=\int\frac{\phi_{i}^{2}\left(\mathbf{r}\right)\phi_{j}^{2}\left(\mathbf{r}^{\prime}\right)}
          {\left|\mathbf{r}-\mathbf{r}^{\prime}\right|}d^{3}\mathbf{r}d^{3}\mathbf{r}^{\prime}

    \f]

    where the orbitals \f$ \phi_{i}(\mathbf{r}) \f$ are s-type Gaussian orbitals (GTOs) of the form

    \f[

    \phi_i\left(\mathbf{r}\right)=\left(\frac{2}{\pi\sigma_{i}^{2}}\right)^{3/4}
          \exp\left(-\frac{\left|\mathbf{r}-\mathbf{R}_{i}\right|^{2}}{\sigma_{i}^{2}}\right)

    \f]

    where \f$ \sigma_i \f$ is a Gaussian screening radius and \f$ \mathbf R_i \f$ is the Cartesian coordinate of atom \f$ i \f$.

    The parameters in this model are the atomic electronegativities \$f \chi_i = \chi_i^0 \$f in volts (V),
    hardnesses \$f \eta_{ii} = \eta_i^0 \f$ in V/e, and screening exponents \$f \alpha_i \$f in Angstroms.

    The default parameter set uses that of QEq.

    References:
      J. Chen and T. J. Martínez, Chem. Phys. Lett. 438 (2007): 315-320.
      doi:10.1021/j100161a070

      J. Chen, D. Hundertmark and T. J. Martínez, J. Chem. Phys. 129 (2008): 214113.
      doi:10.1063/1.3021400

      The Gaussian parameters are published in
      J. Chen and T. J. Martínez, Prog. Theor. Chem. Phys. 19 (2009): 397-416.
      doi:10.1007/978-90-481-2596-8_19

      The method of solving the model is given by
      J. Chen and T. J. Martínez, J. Chem. Phys. 131 (2009): 044114.
      doi:10.1063/1.3183167

    \author Jiahao Chen

    \since version 2.3.
  */

/////////////////////////////////////////////////////////////////
QTPIECharges theQTPIECharges("qtpie"); //Global instance

/////////////////////////////////////////////////////////////////

  void QTPIECharges::ParseParamFile()
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
    Eigen::Vector3d P;
    float radius;

    while (ifs.getline(buffer, BUFF_SIZE)) {
      if (buffer[0] == '#')
        continue;

      tokenize(vs, buffer);
      if (vs.size() < 4)
        continue;

      // Reads in a line of parameters
      // The format is:              code converts to
      //   Element
      //   Electronegativity (V) ->  Electronegativity (a.u.)
      //   Hardness (V/e)        ->  Hardness (a.u.)
      //   radius (Angstrom)     ->  Gaussian exponent (bohr^-2)
      radius = atof(vs[3].c_str())*Angstrom;
      P << atof(vs[1].c_str())*eV, atof(vs[2].c_str())*eV, 1.0/(radius*radius);
      _parameters.push_back(P);
    }
  }

  Eigen::Vector3d QTPIECharges::GetParameters(unsigned int Z, int Q)
  {
    //Returns a triple of numbers: electronegativity (in eV), hardness (in eV), and Gaussian exponent (in bohr^-2)

    Eigen::Vector3d P;
    //For now, completely ignore the formal charge
    if (_parameters.size() == 0)
      ParseParamFile();

    if (Z > 0 && Z < _parameters.size()-1) {
      return _parameters[Z - 1];
    }

    P<<0., 1.e10, 1.e10; //Magic value
    return P;
  }

  bool QTPIECharges::ComputeCharges(OBMol &mol)
  {

     ///////////////////////////////////////////////////////////////////////////////
    //Some OpenBabel bookkeeping that I copied from the Gasteiger scheme
    mol.SetPartialChargesPerceived();

    OBPairData *dp = new OBPairData;
    dp->SetAttribute("PartialCharges");
    dp->SetValue("QTPIE");
    dp->SetOrigin(perceived);
    mol.SetData(dp);


    ///////////////////////////////////////////////////////////////////////////////
    //Read in atomic information from OpenBabel molecule and parameterize

    //Read in total number of atoms
    int i, N = mol.NumAtoms();

    Hardness = Eigen::MatrixXd::Zero(N+1, N+1);
    Voltage = Eigen::VectorXd::Zero(N+1);
    Electronegativity = Eigen::VectorXd::Zero(N);
    Eigen::VectorXd BasisSet = Eigen::VectorXd::Zero(N);

    Eigen::Vector3d Parameters;

    FOR_ATOMS_OF_MOL(atom, mol)
    {
       	Parameters = GetParameters(atom->GetAtomicNum(), atom->GetFormalCharge());
	i = atom->GetIdx() - 1;

	if (Parameters[0] == 0.)
        {
		stringstream msg;
		msg << "Some QTPIE Parameters not found!" << endl
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

    if (Z != 0.) obErrorLog.ThrowError(__FUNCTION__, "Warning, total charge on molecule is not zero. QTPIE routine may give nonsense.", obWarning);


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
		//R = atom1->GetDistance(atom2)*Angstrom;
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

    double OverlapMaxDistance = sqrt( log( (pi/cube(2*SmallestGaussianExponent))
			      / sqr(OverlapThreshold)) /SmallestGaussianExponent);

    //This matrix can be sparse, but I didn't get Eigen's SparseMatrix to
    //play well with this - jiahao@mit.edu 2010-04-20 eigen-2.0.12 r3691
    Eigen::MatrixXd Overlap = Eigen::MatrixXd::Zero(N,N);
    double OverlapVal;
    FOR_ATOMS_OF_MOL(atom1, mol)
    {
	i = atom1->GetIdx() - 1;
        FOR_ATOMS_OF_MOL(atom2, mol)
	{
	    j = atom2->GetIdx() - 1;
	    if (i>j)
	    {
	        R = atom1->GetDistance(j+1)*Angstrom;
	        if (R<OverlapMaxDistance)
	        {
		    OverlapVal = OverlapInt(BasisSet[i], BasisSet[j], R);
		    Overlap(i,j) = OverlapVal;
		    Overlap(j,i) = OverlapVal;
	        }
	    }
	}
    }

    // Calculate normalization factors
    Eigen::VectorXd OvNorm(N);
    for (i=0; i<N; i++) OvNorm[i] = 1.0 / (1.0 + Overlap.row(i).sum());

    // Calculate voltages
    double PotentialDiff, Norm, ThisVoltage, ThisOverlap;
    for (i=0; i<N; i++)
    {
	Norm = OvNorm[i];
	ThisVoltage = 0.;
	for (j=0; j<N; j++)
	{
		PotentialDiff = Electronegativity[i] - Electronegativity[j];
		ThisOverlap = Overlap(i,j);
		if (ThisOverlap > OverlapThreshold)
		       ThisVoltage -= PotentialDiff * Norm* ThisOverlap;
	}

	Voltage[i] = ThisVoltage;
    }

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
    cout << "Ready to solve QTPIE problem" << endl << endl;
    cout << "Hardness matrix:" << endl << Hardness << endl << endl;
    cout << "Voltage vector:" << endl << Voltage << endl << endl;
    cout << "Overlap matrix:" << endl << Overlap << endl << endl;
    cout << "Overlap norm:" << endl << OvNorm << endl << endl;
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
    obErrorLog.ThrowError(__FUNCTION__, "QTPIE charges successfully assigned.", obInfo);
    return true;
  }


/// Calculates Coulomb integral
double QTPIECharges::CoulombInt(double a, double b, double R)
{
	double p = sqrt(a * b / (a + b));
	return erf(p * R) / R;
}

double QTPIECharges::OverlapInt(double a, double b, double R)
{
	double p = a + b;
	double q = a * b / p;
	return pow(4*q/p, 0.75) * exp(-q*R*R);
}

/// Here's a wrapper around the Eigen solver routine
bool QTPIECharges::solver(Eigen::MatrixXd A, Eigen::VectorXd b, Eigen::VectorXd &x, const double NormThreshold)
{
    // using a LU factorization
#ifdef HAVE_EIGEN3
    bool SolverOK = true;
    x = A.partialPivLu().solve(b);
#else
    bool SolverOK = A.lu().solve(b, &x);
#endif
    //bool SolverOK = A.svd().solve(b, &x);

    Eigen::VectorXd resid = A*x - b;
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

#ifdef HAVE_EIGEN3
        x = A.jacobiSvd().solve(b);
#else
        SolverOK = A.svd().solve(b, &x);
#endif
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

//! \file qtpie.cpp
//! \brief Assign QTPIE partial charges.

/**********************************************************************
align.cpp - Align two molecules or vectors of vector3

Copyright (C) 2010 by Noel M. O'Boyle

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

#include <openbabel/babelconfig.h>

#include <vector>
#include <climits> // UINT_MAX

#include <openbabel/math/align.h>
#include <openbabel/atom.h>
#include <openbabel/oberror.h>
#include <openbabel/obiter.h>
#include <openbabel/graphsym.h>
#include <openbabel/math/vector3.h>
#include <openbabel/elements.h>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/SVD>
#include <Eigen/LU>

using namespace std;

namespace OpenBabel
{
  extern OBMessageHandler obErrorLog;

  OBAlign::OBAlign(bool includeH, bool symmetry) : _method(OBAlign::Kabsch)
  {
    _ready = false;
    _symmetry = symmetry;
    _includeH = includeH;
    _prefmol = 0;
  }

  OBAlign::OBAlign(const vector<vector3> &ref, const vector<vector3> &target) : _method(OBAlign::Kabsch)
  {
    SetRef(ref);
    SetTarget(target);
    _symmetry = false;
    _prefmol = 0;
  }

  OBAlign::OBAlign(const OBMol &refmol, const OBMol &targetmol, bool includeH, bool symmetry) : _method(OBAlign::Kabsch)
  {
    _symmetry = symmetry;
    _includeH = includeH;
    SetRefMol(refmol);
    SetTargetMol(targetmol);
  }

  void OBAlign::VectorsToMatrix(const vector<vector3> *pcoords, Eigen::MatrixXd &coords) {

    vector<vector3>::size_type N = pcoords->size();
    coords.resize(3, N);

    // Create a 3xN matrix of the coords
    vector<vector3>::const_iterator it;
    vector<vector3>::size_type colm;
    for (colm=0,it=pcoords->begin();colm<N;++colm,++it)
      coords.col(colm) = Eigen::Vector3d( it->AsArray() );
  }

  Eigen::Vector3d OBAlign::MoveToOrigin(Eigen::MatrixXd &coords) {

    vector<vector3>::size_type N = coords.cols();

    // Find the centroid
    Eigen::Vector3d centroid;
    centroid = coords.rowwise().sum() / N;

    // Subtract the centroids
    for (vector<vector3>::size_type i=0; i<N; ++i)
      coords.col(i) -= centroid;
    return centroid;
  }

  void OBAlign::SetRef(const vector<vector3> &ref) {
    _pref = &ref;
    VectorsToMatrix(_pref, _mref);
    _ref_centr = MoveToOrigin(_mref);

    _ready = false;
  }

  void OBAlign::SetTarget(const vector<vector3> &target) {
    _ptarget = &target;
    VectorsToMatrix(_ptarget, _mtarget);
    _target_centr = MoveToOrigin(_mtarget);

    _ready = false;
  }

  void OBAlign::SetRefMol(const OBMol &refmol) {
    _prefmol = &refmol;

    // Set up the BitVec for the hydrogens and store the refmol coords
    _frag_atoms.Clear();
    _frag_atoms.Resize(refmol.NumAtoms() + 1);
    _refmol_coords.resize(0);
    OBAtom* atom;
    int delta = 1;
    _newidx.resize(0);

    for (unsigned int i=1; i<=refmol.NumAtoms(); ++i) {
      atom = refmol.GetAtom(i);
      if (_includeH || atom->GetAtomicNum() != OBElements::Hydrogen) {
        _frag_atoms.SetBitOn(i);
        _newidx.push_back(i - delta);
        _refmol_coords.push_back(atom->GetVector());
      }
      else {
        delta++;
        _newidx.push_back(UINT_MAX);
      }
    }
    SetRef(_refmol_coords);

    if (_symmetry) {
      FindAutomorphisms((OBMol*)&refmol, _aut, _frag_atoms);
    }
  }

  void OBAlign::SetTargetMol(const OBMol &targetmol) {
    _ptargetmol = &targetmol;
    _targetmol_coords.resize(0);
    OBAtom const *atom;
    for (unsigned int i=1; i<=targetmol.NumAtoms(); ++i) {
      atom = targetmol.GetAtom(i);
      if (_includeH || atom->GetAtomicNum() != OBElements::Hydrogen)
        _targetmol_coords.push_back(atom->GetVector());
    }
    SetTarget(_targetmol_coords);
  }

  void OBAlign::SetMethod(OBAlign::AlignMethod method) {
    _method = method;
  }

/* Evaluates the Newton-Raphson correction for the Horn quartic.
   only 11 FLOPs */
  static double eval_horn_NR_corrxn(const vector<double> &c, const double x)
  {
    double x2 = x*x;
    double b = (x2 + c[2])*x;
    double a = b + c[1];

    return((a*x + c[0])/(2.0*x2*x + b + a));
  }

  /* Newton-Raphson root finding */
  static double QCProot(const vector<double> &coeff, double guess, const double delta)
  {
    int             i;
    double          oldg;
    double initialg = guess;

    for (i = 0; i < 50; ++i)
    {
        oldg = guess;
        /* guess -= (eval_horn_quart(coeff, guess) / eval_horn_quart_deriv(coeff, guess)); */
        guess -= eval_horn_NR_corrxn(coeff, guess);

        if (fabs(guess - oldg) < fabs(delta*guess))
            return(guess);
    }

    return initialg + 1.0; // Failed to converge!
  }

  vector<double> CalcQuarticCoeffs(const Eigen::Matrix3d &M)
  {
    vector<double> coeff(4);

    double          Sxx, Sxy, Sxz, Syx, Syy, Syz, Szx, Szy, Szz;
    double          Szz2, Syy2, Sxx2, Sxy2, Syz2, Sxz2, Syx2, Szy2, Szx2,
                    SyzSzymSyySzz2, Sxx2Syy2Szz2Syz2Szy2, Sxy2Sxz2Syx2Szx2,
                    SxzpSzx, SyzpSzy, SxypSyx, SyzmSzy,
                    SxzmSzx, SxymSyx, SxxpSyy, SxxmSyy;

#ifdef HAVE_EIGEN3
    Eigen::MatrixXd M_sqr = M.array().square();
#else
    Eigen::MatrixXd M_sqr = M.cwise().square();
#endif

    Sxx = M(0, 0);
    Sxy = M(1, 0);
    Sxz = M(2, 0);
    Syx = M(0, 1);
    Syy = M(1, 1);
    Syz = M(2, 1);
    Szx = M(0, 2);
    Szy = M(1, 2);
    Szz = M(2, 2);

    Sxx2 = Sxx * Sxx;
    Syy2 = Syy * Syy;
    Szz2 = Szz * Szz;

    Sxy2 = Sxy * Sxy;
    Syz2 = Syz * Syz;
    Sxz2 = Sxz * Sxz;

    Syx2 = Syx * Syx;
    Szy2 = Szy * Szy;
    Szx2 = Szx * Szx;

    SyzSzymSyySzz2 = 2.0*(Syz*Szy - Syy*Szz);
    Sxx2Syy2Szz2Syz2Szy2 = Syy2 + Szz2 - Sxx2 + Syz2 + Szy2;

    /* coeff[4] = 1.0; */
    /* coeff[3] = 0.0; */
    // coeff[2] = -2.0 * (Sxx2 + Syy2 + Szz2 + Sxy2 + Syx2 + Sxz2 + Szx2 + Syz2 + Szy2);
    coeff[2] = -2.0 * M_sqr.sum();
    coeff[1] = 8.0 * (Sxx*Syz*Szy + Syy*Szx*Sxz + Szz*Sxy*Syx - Sxx*Syy*Szz - Syz*Szx*Sxy - Szy*Syx*Sxz);

    SxzpSzx = Sxz+Szx;
    SyzpSzy = Syz+Szy;
    SxypSyx = Sxy+Syx;
    SyzmSzy = Syz-Szy;
    SxzmSzx = Sxz-Szx;
    SxymSyx = Sxy-Syx;
    SxxpSyy = Sxx+Syy;
    SxxmSyy = Sxx-Syy;
    Sxy2Sxz2Syx2Szx2 = Sxy2 + Sxz2 - Syx2 - Szx2;

    coeff[0] = Sxy2Sxz2Syx2Szx2 * Sxy2Sxz2Syx2Szx2
             + (Sxx2Syy2Szz2Syz2Szy2 + SyzSzymSyySzz2) * (Sxx2Syy2Szz2Syz2Szy2 - SyzSzymSyySzz2)
             + (-(SxzpSzx)*(SyzmSzy)+(SxymSyx)*(SxxmSyy-Szz)) * (-(SxzmSzx)*(SyzpSzy)+(SxymSyx)*(SxxmSyy+Szz))
             + (-(SxzpSzx)*(SyzpSzy)-(SxypSyx)*(SxxpSyy-Szz)) * (-(SxzmSzx)*(SyzmSzy)-(SxypSyx)*(SxxpSyy+Szz))
             + (+(SxypSyx)*(SyzpSzy)+(SxzpSzx)*(SxxmSyy+Szz)) * (-(SxymSyx)*(SyzmSzy)+(SxzpSzx)*(SxxpSyy+Szz))
             + (+(SxypSyx)*(SyzmSzy)+(SxzmSzx)*(SxxmSyy-Szz)) * (-(SxymSyx)*(SyzpSzy)+(SxzmSzx)*(SxxpSyy-Szz));

    return coeff;
  }

  void OBAlign::TheobaldAlign(const Eigen::MatrixXd &mtarget)
  {
    // M = B(t) times A (where A, B are N x 3 matrices)
    Eigen::Matrix3d M = mtarget * _mref.transpose();

    // Maximum value for lambda is (Ga + Gb) / 2
    double innerprod = mtarget.squaredNorm() + _mref.squaredNorm();

    vector<double> coeffs = CalcQuarticCoeffs(M);
    double lambdamax = QCProot(coeffs, 0.5 * innerprod, 1e-6);
    if (lambdamax > (0.5 * innerprod))
      _fail = true;
    else {
      double sqrdev = innerprod - (2.0 * lambdamax);
      _rmsd = sqrt(sqrdev / mtarget.cols());
    }
  }

  void OBAlign::SimpleAlign(const Eigen::MatrixXd &mtarget)
  {
    // Covariance matrix C = X times Y(t)
    Eigen::Matrix3d C = _mref * mtarget.transpose();

    // Singular Value Decomposition of C into USV(t)
#ifdef HAVE_EIGEN3
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(C, Eigen::ComputeFullU | Eigen::ComputeFullV);
#else
    Eigen::SVD<Eigen::Matrix3d> svd(C);
#endif

    // Prepare matrix T
    double sign = (C.determinant() > 0) ? 1. : -1.; // Sign of determinant
    Eigen::Matrix3d T = Eigen::Matrix3d::Identity();
    T(2,2) = sign;

    // Optimal rotation matrix, U, is V T U(t)
    _rotMatrix = svd.matrixV() * T * svd.matrixU().transpose();

    // Rotate target using rotMatrix
    _result = _rotMatrix.transpose() * mtarget;

    Eigen::MatrixXd deviation = _result - _mref;
#ifdef HAVE_EIGEN3
    Eigen::MatrixXd sqr = deviation.array().square();
#else
    Eigen::MatrixXd sqr = deviation.cwise().square();
#endif
    double sum = sqr.sum();
    _rmsd = sqrt( sum / sqr.cols() );

  }

  bool OBAlign::Align()
  {
    vector<vector3>::size_type N = _ptarget->size();

    if (_pref->size() != N) {
      obErrorLog.ThrowError(__FUNCTION__, "Cannot align the reference and target as they are of different size" , obError);
      return false;
    }

    if (!_symmetry || _aut.size() == 1) {
      if (_method == OBAlign::Kabsch)
        SimpleAlign(_mtarget);
      else
        TheobaldAlign(_mtarget);
    }
    else {  // Iterate over the automorphisms

      // ...for storing the results from the lowest rmsd to date
      double min_rmsd = DBL_MAX;
      Eigen::MatrixXd result, rotMatrix;

      // Try all of the symmetry-allowed permutations
      OBIsomorphismMapper::Mappings::const_iterator cit;
      Eigen::MatrixXd mtarget(_mtarget.rows(), _mtarget.cols());

      for (unsigned int k = 0; k < _aut.size(); ++k) {
        // Rearrange columns of _mtarget for this permutation
        unsigned int i=0;
        for (unsigned int j=1; j<=_prefmol->NumAtoms(); ++j) {
          if (_frag_atoms.BitIsSet(j)) {
            for (std::size_t l = 0; l < _aut[k].size(); ++l)
              if (_aut[k][l].first == j - 1) {
                mtarget.col(i) = _mtarget.col(_newidx[_aut[k][l].second]);
                break;
              }
            i++;
          }
        }
        if (_method == OBAlign::Kabsch)
          SimpleAlign(mtarget);
        else
          TheobaldAlign(mtarget);
        if (_rmsd < min_rmsd) {
          min_rmsd = _rmsd;
          result = _result;
          rotMatrix = _rotMatrix;
        }
      }

      // Restore the best answer from memory
      _rmsd = min_rmsd;
      _result = result;
      _rotMatrix = rotMatrix;
    }

    _ready = true;
    return true;
  }

  matrix3x3 OBAlign::GetRotMatrix()
  {
    if (!_ready) {
      obErrorLog.ThrowError(__FUNCTION__, "Rotation matrix not available until you call Align()" , obError);
      return matrix3x3();
    }

    // Convert Eigen::Matrix to matrix3x3
    double rot[3][3];
    for (int row=0; row<3; ++row)
       for (int col=0; col<3; ++col)
         rot[col][row] = _rotMatrix(row, col); // Return in form suitable for use in expressions like "result *= rotMatrix";
    matrix3x3 rotmat = matrix3x3(rot);

    return rotmat;
  }

  vector<vector3> OBAlign::GetAlignment() {
    vector<vector3> aligned_coords;
    if (!_ready) {
      obErrorLog.ThrowError(__FUNCTION__, "Alignment not available until you call Align()" , obError);
      return aligned_coords;
    }

    if (!_prefmol || _includeH) {
      // Add back the centroid of the reference and convert to vv3
      Eigen::Vector3d tmp;
      aligned_coords.reserve(_result.cols());
      for (int i=0; i<_result.cols(); ++i) {
        tmp = _result.col(i) + _ref_centr;
        aligned_coords.push_back(vector3(tmp(0), tmp(1), tmp(2)));
      }
    }
    else { // Need to deal with the case where hydrogens were excluded
      vector<vector3> target_coords;
      for (unsigned int i=1; i<=_ptargetmol->NumAtoms(); ++i)
        target_coords.push_back(_ptargetmol->GetAtom(i)->GetVector());
      Eigen::MatrixXd mtarget;
      VectorsToMatrix(&target_coords, mtarget);

      // Subtract the centroid of the non-H atoms
      for (unsigned int i=0; i<mtarget.cols(); ++i)
        mtarget.col(i) -= _target_centr;

      // Rotate
      Eigen::MatrixXd result = mtarget.transpose() * _rotMatrix;
      result.transposeInPlace();

      // Add back the centroid of the reference and convert to vv3
      Eigen::Vector3d tmp;
      aligned_coords.reserve(_result.cols());
      for (int i=0; i<result.cols(); ++i) {
        tmp = result.col(i) + _ref_centr;
        aligned_coords.push_back(vector3(tmp(0), tmp(1), tmp(2)));
      }
    }

    return aligned_coords;
  }

  bool OBAlign::UpdateCoords(OBMol* target) {
    if (!_ready) {
      obErrorLog.ThrowError(__FUNCTION__, "Alignment not available until you call Align()" , obError);
      return false;
    }

    vector<vector3> newcoords = GetAlignment();
    if (newcoords.size() != target->NumAtoms()) {
      obErrorLog.ThrowError(__FUNCTION__, "Cannot update the target molecule with the alignment coordinates as they are of different size" , obError);
      return false;
    }

    int i = 0;
    FOR_ATOMS_OF_MOL(a, *target) {
      a->SetVector(newcoords.at(i));
      i++;
    }

    return true;
  }

  double OBAlign::GetRMSD() {
    if (!_ready) {
      obErrorLog.ThrowError(__FUNCTION__, "RMSD not available until you call Align()" , obError);
      return (double) NULL;
    }

    return _rmsd;
  }

} // namespace OpenBabel

//! \file align.cpp
//! \brief Handle 3D coordinates.

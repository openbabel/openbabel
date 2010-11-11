/**********************************************************************
eem.cpp - A OBChargeModel to handle EEM charges from Bultinck.

Copyright (C) 2005-2010 by Silicos NV

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
#include <openbabel/chargemodel.h>
#include <openbabel/mol.h>
#include <openbabel/molchrg.h>

namespace OpenBabel
{

  class EEMCharges : public OBChargeModel
  {
  public:
    EEMCharges(const char* ID) : OBChargeModel(ID, false){};
    const char* Description(){ return "Assign Electronegativity Equilization Method (EEM) atomic partial charges"; }

    /// \return whether partial charges were successfully assigned to this molecule
    bool ComputeCharges(OBMol &mol);

    double DipoleScalingFactor() { return 1.0; } // fit from regression

  private:
    void _solveMatrix(double**, double*, unsigned int);
    void _luDecompose(double**, std::vector<int>&, unsigned int);
    void _luSolve(double**, std::vector<int>&, double*, unsigned int);
    void _swapRows(double*, unsigned int, unsigned int);
    void _swapRows(double**, unsigned int, unsigned int, unsigned int);
  };

  /////////////////////////////////////////////////////////////////
  EEMCharges theEEMCharges("eem"); //Global instance

  /////////////////////////////////////////////////////////////////

  bool EEMCharges::ComputeCharges(OBMol &mol)
  {
    mol.SetPartialChargesPerceived();

    // Copied from spectrophore.cpp
    // CHI and ETA
    unsigned int _nAtoms = mol.NumAtoms();
    unsigned int dim(_nAtoms + 1);
    std::vector<double> CHI(dim);
    double** ETA = new double*[dim];
    for (unsigned int i = 0; i < dim; ++i)
      {
        ETA[i] = new double[dim];
      }
    double totalCharge(0.0);
    unsigned int i(0);
    unsigned int n;
    double hardness;
    double electronegativity;
    for (OpenBabel::OBMolAtomIter atom(mol); atom; ++atom)
      {
        n = (unsigned int) atom->GetAtomicNum();
        switch (n) // this should move to a parameter file
          {
          case 1:  // H
            hardness = 0.65971;
            electronegativity = 0.20606;
            break;
          case 3:  // Li
            hardness = 0.32966;
            electronegativity = 0.36237;
            break;
          case 5:  // B
            hardness = 0.32966;
            electronegativity = 0.36237;
            break;
          case 6:  // C
            hardness = 0.32966;
            electronegativity = 0.36237;
            break;
          case 7:  // N
            hardness = 0.34519;
            electronegativity = 0.49279;
            break;
          case 8:  // O
            hardness = 0.54428;
            electronegativity = 0.73013;
            break;
          case 9:  // F
            hardness = 0.72664;
            electronegativity = 0.72052;
            break;
          case 11: // Na
            hardness = 0.32966;
            electronegativity = 0.36237;
            break;
          case 12: // Mg
            hardness = 0.32966;
            electronegativity = 0.36237;
            break;
          case 14: // Si
            hardness = 0.32966;
            electronegativity = 0.36237;
            break;
          case 15: // P
            hardness = 0.32966;
            electronegativity = 0.36237;
            break;
          case 16: // S
            hardness = 0.20640;
            electronegativity = 0.62020;
            break;
          case 17: // Cl
            hardness = 0.32966;
            electronegativity = 0.36237;
            break;
          case 19: // K
            hardness = 0.32966;
            electronegativity = 0.36237;
            break;
          case 20: // Ca
            hardness = 0.32966;
            electronegativity = 0.36237;
            break;
          case 26: // Fe
            hardness = 0.32966;
            electronegativity = 0.36237;
            break;
          case 29: // Cu
            hardness = 0.32966;
            electronegativity = 0.36237;
            break;
          case 30: // Zn
            hardness = 0.32966;
            electronegativity = 0.36237;
            break;
          case 35: // Br
            hardness = 0.54554;
            electronegativity = 0.70052;
            break;
          case 53: // I
            hardness = 0.30664;
            electronegativity = 0.68052;
            break;
          default:
            hardness = 0.65971;
            electronegativity = 0.20606;
            break;
          }

        CHI[i] = -electronegativity;
        ETA[i][i] = 2.0 * hardness;

        // Adjust the total molecular charge
        totalCharge += atom->GetFormalCharge();

        // Increment
        ++i;
      }

    // Complete CHI
    CHI[_nAtoms] = totalCharge;

    // Complete ETA
    double d;
    OBAtom *rAtom, *cAtom;
    for (unsigned int r = 0; r < _nAtoms; ++r)
      {
        rAtom = mol.GetAtom(r+1); // Atom index
        for (unsigned int c = r + 1; c < _nAtoms; ++c)
          {
            cAtom = mol.GetAtom(c+1); // Atom index
            ETA[r][c] = 0.529176 / cAtom->GetDistance(rAtom);     // 0.529176: Angstrom to au
            ETA[c][r] = ETA[r][c];
          }
      }
    for (int i = 0; i < dim; ++i)
      {
        ETA[i][_nAtoms] = -1.0;
        ETA[_nAtoms][i] = +1.0;
      }
    ETA[_nAtoms][_nAtoms] = 0.0;

    // Solve the matrix equation
    _solveMatrix(ETA, &(CHI[0]), dim);    // CHI will contain the values

    OBAtom *atom;
    for (unsigned int i = 0; i < _nAtoms; ++i)
      {
        atom = mol.GetAtom(i+1); // atom index issue
        atom->SetPartialCharge(CHI[i]);
      }

    OBChargeModel::FillChargeVectors(mol);

    return true;
  }

  void
  EEMCharges::_solveMatrix(double** A, double* B, unsigned int dim)
  {
    std::vector<int> temp(dim);
    _luDecompose(A, temp, dim);
    _luSolve(A, temp, B, dim);
  }



  void
  EEMCharges::_luDecompose(double** A, std::vector<int>& I, unsigned int dim)
  {
    int i, j, k, kMax, iMax;
    std::vector<double> vScales(dim, 0);
    double maxVal = 0, dummy = 0;
    double * pRowi = NULL;

    // first find the highest pivot element in each row and store it for implicit scaling
    for (i = 0; i < dim; ++i)
      {
        maxVal = 0.0;
        for (j = 0; j < dim; ++j)
          {
            if ((dummy=fabs(A[i][j])) > maxVal)
              {
                maxVal = dummy;
              }
          }
        if (maxVal == 0)
          {
            std::cerr << "EEMCharges: Warning singular matrix..." << std::endl;
          }

        vScales[i] = 1.0 / maxVal;
      }

    std::vector<double> colJ(dim); // variable to store local copy of column

    // loop over columns
    for (j = 0; j < dim; ++j)
      {
        // make a local copy of column j
        for (i = 0; i < dim; ++i) colJ[i] = A[i][j];
        for (i = 0; i < dim; ++i)
          {
            pRowi = A[i];
            dummy = pRowi[j];
            kMax = i < j ? i : j;
            for (k = 0; k < kMax; ++k) dummy -= pRowi[k] * colJ[k];
            colJ[i] = dummy;
            pRowi[j] = colJ[i];
          }

        // search largest pivot element
        maxVal = 0.0;
        iMax = j;
        for (i = j + 1; i < dim; ++i)
          {
            if ((dummy = fabs(colJ[i]) * vScales[i]) >= maxVal)
              {
                maxVal = dummy;
                iMax = i;
              }
          }

        // check if we need to interchange rows
        if (j != iMax) // if current column index is not the maximal row index we need to interchange
          {
            // std::cerr << "Swap rows: " << iMax << " <-> " << j << std::endl;
            _swapRows(A, iMax, j, dim);
            vScales[iMax] = vScales[j];
          }
        // store row index in I
        I[j] = iMax;

        // finally divide by the pivot element
        if (j != dim - 1)
          {
            dummy = 1.0 / A[j][j]; // A.GetValueAt(j,j);
            for (i = j + 1; i < dim; ++i) A[i][j] *= dummy;
          }


      } // next column

    return;
  }



  void
  EEMCharges::_luSolve(double** A, std::vector<int>& I, double* B, unsigned int dim)
  {
    int i, k;

    for (int i = 0; i < dim; ++i) _swapRows(B, i, I[i]);

    // forward substitution pass
    for (k = 0; k < dim; ++k)
      {
        for (i = k+1; i < dim; ++i)
          {
            B[i] -= A[i][k] * B[k];
          }
      }

    // do the backsubstitution
    for (i = dim - 1; i >= 0; --i)
      {
        B[i] /= A[i][i];
        for (k = 0; k < i; ++k)
          {
            B[k] -= A[k][i] * B[i];
          }
      }

    return;
  }

  void
  EEMCharges::_swapRows(double** _pMatrix, unsigned int i, unsigned int j, unsigned int nCols)
  {
    double dummy;
    for (unsigned int k = 0; k < nCols; ++k)         // loop over all columns
      {
        dummy = _pMatrix[i][k];
        _pMatrix[i][k] = _pMatrix[j][k];
        _pMatrix[j][k] = dummy;
      }
    return;
  }

  void
  EEMCharges::_swapRows(double* _pMatrix, unsigned int i, unsigned int j)
  {
    double dummy;
    dummy = _pMatrix[i];
    _pMatrix[i] = _pMatrix[j];
    _pMatrix[j] = dummy;
    return;
  }

}//namespace

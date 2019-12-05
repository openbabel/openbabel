/**********************************************************************
align.h - Align two molecules or vectors of vector3

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

#ifndef OB_ALIGN_H
#define OB_ALIGN_H

#include <openbabel/mol.h>
#include <openbabel/math/vector3.h>
#include <openbabel/math/matrix3x3.h>
#include <openbabel/isomorphism.h>
#include <Eigen/Core>

namespace OpenBabel
{
  /**
   * \class OBAlign align.h <openbabel/math/align.h>
   * \brief Perform a least-squares alignment of two molecules or two vectors of vector3 objects
   *
   * This class may be used to perform a least-squares alignment of two OBMol
   * objects or two vector<vector3> objects. The Kabsch algorithm is used
   * for the alignment.
   *
   * During the alignment, the Target is aligned to the Reference. Note that
   * mutiple alignments to the same Reference will be much faster than multiple
   * alignments to the same Target. When carrying out multiple alignments,
   * a single OBAlign instance should be reused by calling
   * SetTarget() or SetTargetMol() for each additional Target and then calling
   * Align().
   *
   * When aligning molecules, the atoms of the two molecules must be in the same
   * order for the results to be sensible. By default, hydrogens are not
   * included in the least-squares fitting procedure (set @p includeH to
   * true if you wish to include them) and so the resulting RMSD is the
   * heavy-atom RMSD (which is usually what you want).
   *
   * By default, symmetry is taken
   * into account when comparing molecules. For example, if a benzene is flipped
   * by 180 degrees along one of its 2-fold symmetry axes, it will only have an
   * RMSD of 0 (with respect to its original orientation) if symmetry is
   * enabled. To turn off handling of symmetry set @p symmetry to false (this
   * will speed up the alignment).
   *
   * Note that neither the Target nor the Reference
   * are modified by the algorithm. As a result, to update a TargetMol with the
   * new coordinates from the alignment, you need to use UpdateCoords().
   *
   * @since version 2.3
   */
  class OBAPI OBAlign {
  public:
    ///@name Constructors
    //@{
    /**
     * If this constructor is used, the Target and Reference must be
     * set using SetRef/SetRefMol and SetTarget/SetTargetMol before running
     * the alignment.
     */
    OBAlign(bool includeH=false, bool symmetry=true);
    /**
     * Align two molecules.
     */
    OBAlign(const OBMol &refmol, const OBMol &targetmol, bool includeH=false, bool symmetry=true);
    /**
     * Align two vectors of vector3 objects.
     */
    OBAlign(const std::vector<vector3> &ref, const std::vector<vector3> &target);
    //@}

    ///@name Partial Setup
    //@{
    /**
     * Set the Reference (to which the Target will be aligned) in
     * terms of a vector of vector3 objects. Note that it is faster
     * to perform multiple alignments to the same Reference, rather than
     * multiple alignments to the same Target.
     */
    void SetRef(const std::vector<vector3> &ref);
    /**
     * Set the Target (which will be aligned to the Reference) in
     * terms of a vector of vector3 objects.
     */
    void SetTarget(const std::vector<vector3> &target);
    /**
     * Set the Reference Molecule (to which the Target Molecule must
     * be aligned). Note that is faster to perform multiple alignments
     * to the same Reference Molecule, rather than multple alignments
     * to the same Target Molecule.
     */
    void SetRefMol(const OBMol &refmol);
    /**
     * Set the Target Molecule (which will be aligned to the
     * Reference Molecule).
     */
    void SetTargetMol(const OBMol &targetmol);
    //@}

    ///@name Execute the alignment
    //@{
    /**
     * Align the Target to the Reference using a least-squares alignment.
     */
    bool Align();

    enum AlignMethod {
      Kabsch = 0,   // Returns matrix and RMSD
      QCP    = 1    // Returns just RMSD (fast)
    };

    void SetMethod(enum AlignMethod method);
    //@}

    ///@name Access the result of the alignment
    //@{
    /**
     * Return the root mean squared deviation of the target from
     * the reference. This function should only
     * be called after running the alignment (using Align()).
     */
    double GetRMSD();
    /**
     * Return the rotation matrix associated with the least squares
     * alignment. This function should only
     * be called after running the alignment (using Align()).
     * 
     * The following example shows how to use the rotation matrix
     * to rotate all of the atoms in a molecule. 
     * \code
     * matrix3x3 rotmatrix = align.GetRotMatrix();
     * for (unsigned int i = 1; i <= mol.NumAtoms(); ++i) {
     *    vector3 tmpvec = mol.GetAtom(i)->GetVector();
     *    tmpvec *= rotmatrix; //apply the rotation
     *    mol.GetAtom(i)->SetVector(tmpvec);
     * }
     * \endcode
     * Note that if you wish to use the rotation matrix to find the
     * aligned coordinates (that is, the same coordinates returned by
     * GetAlignment()), you should first translate the set of coordinates
     * to the origin by subtracting the centroid, then apply the rotation,
     * and finally add the centroid of the reference coordinates.
     */
    matrix3x3 GetRotMatrix();
    /**
     * Return the actual alignment of the Target to the Reference
     * in terms of a vector of vector3 objects. If you want an OBMol
     * with the aligned coordinates, you should use UpdateCoords() instead.
     * This function should only
     * be called after running the alignment (using Align()).
     */
    std::vector<vector3> GetAlignment();
    /**
     * Set the coordinates of an OBMol to those from the alignment.
     * This function should only
     * be called after running the alignment (using Align()).
     */
    bool UpdateCoords(OBMol* target);
    //@}

  private:
    bool _ready;
    bool _fail;
    bool _symmetry;
    bool _includeH;
    enum AlignMethod _method;
    double _rmsd;
    OBBitVec _frag_atoms;
    Automorphisms _aut;
    const OBMol* _prefmol;
    const OBMol* _ptargetmol;
    Eigen::MatrixXd _rotMatrix;
    Eigen::Vector3d _ref_centr, _target_centr;
    const std::vector<vector3> *_pref;
    const std::vector<vector3> *_ptarget;
    std::vector<vector3> _refmol_coords;
    std::vector<vector3> _targetmol_coords;
    Eigen::MatrixXd _result;
    Eigen::MatrixXd _mref, _mtarget;
    void VectorsToMatrix(const std::vector<vector3> *pcoords, Eigen::MatrixXd &coords);
    Eigen::Vector3d MoveToOrigin(Eigen::MatrixXd &coords);
    void SimpleAlign(const Eigen::MatrixXd &mtarget);
    void TheobaldAlign(const Eigen::MatrixXd &mtarget);
    // Generate a mapping from the permutation map to the index of
    // correct column in _mtarget. Need to handle the fact that the
    // permutation group contains non-fragment atoms.
    // For example, map(213465) will be converted to newidx(102354).
    // If the atom with Idx=3 is not in the fragment, it will be
    // converted to newidx(10X243) instead.
    std::vector<unsigned int> _newidx;
  };
}

#endif // OB_ALIGN_H

/// @file align.h
/// @brief Perform a least-squares alignment of two molecules or two vectors of vector3 objects

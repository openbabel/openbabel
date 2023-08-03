/**********************************************************************
distgeom.h - Distance Geometry generation and sampling

  Copyright (C) 2011 by Tim Vandermeersch
  Copyright (C) 2012 by Geoffrey Hutchison

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

#ifndef OB_DISTGEOM_H
#define OB_DISTGEOM_H

#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>

#include <iostream>

#ifndef OBAPI
  #define OBAPI
#endif

#ifdef HAVE_EIGEN

#include <Eigen/Core>
#include <LBFGS.h>

namespace OpenBabel {

  class DistanceGeometryPrivate;
  class OBCisTransStereo;

  class TetrahedralInfo {
    int c;
    std::vector<unsigned long> nbrs;
    double lb, ub;
    public:
    TetrahedralInfo(int center, std::vector<unsigned long> neighbors,
                    double lower_bound, double upper_bound) :
                    c(center), nbrs(neighbors),
                    lb(lower_bound), ub(upper_bound) {}
    int GetCenter() {
      return c;
    }
    std::vector<unsigned long> GetNeighbors() {
      return nbrs;
    }
    double GetUpperBound() {
      return ub;
    }
    double GetLowerBound() {
      return lb;
    }
  };

  class OBAPI OBDistanceGeometry {
    friend class DistgeomFunc;
  public:
    OBDistanceGeometry();
    OBDistanceGeometry(const OBMol &mol, bool useCurrentGeometry);
    ~OBDistanceGeometry();

    /**
     * Setup this instance with the specified molecule.
     *
     * @param mol The molecule to use
     * @param useCurrentGeom Whether to use the current bond distances and angles
     * for the bounds matrix
     *
     * \return Success or failure
     */
    bool Setup(const OBMol &mol, bool useCurrentGeom = false);

    void Generate();
    void AddConformer();
    void GetConformers(OBMol &mol);

    /**
     * Check if last call to AddConformer was successful.
     *
     * \return Success or failure.
     */
    bool WasSuccessful() const;

    /**
     * Convenience method to set up this molecule, generate a geometry and return it
     *
     * \return Success or failure
     */
    bool GetGeometry(OBMol &mol, bool useCurrentGeom = false);

    //! \return The bounds matrix after setup (e.g., for debugging)
    Eigen::MatrixXf GetBoundsMatrix();
    /**
     * \brief Set the bounds matrix explicitly
     * \return Success or failure (e.g., bounds matrix does not match the number of atoms)
     */
    bool SetBoundsMatrix(const Eigen::MatrixXf bounds);
    float GetUpperBounds(int i, int j);
    float GetLowerBounds(int i, int j);
    unsigned int GetDimension() {return dim;};
    std::vector<TetrahedralInfo>  _stereo;       //!< Internal private data, including stereo info
  private:
    OBMol                     _mol;
    std::vector<OBGenericData*> _vdata;
    DistanceGeometryPrivate  *_d;    //!< Internal private data, including bounds matrix
    Eigen::VectorXd _coord;          // one-dimensional vector containing coordinates of atoms
    std::string input_smiles;

    unsigned int dim;

    bool generateInitialCoords();
    bool firstMinimization();
    bool minimizeFourthDimension();
    
    //! \brief Set the default upper bounds for the constraint matrix
    //! Upper bounds = maximum length of the molecule, or 1/2 the body diagonal in a unit cell
    void SetUpperBounds();
    //! \brief Update the upper and lower bounds based on bonding distances
    //! @param useCurrentGeom Whether to use the current bond distances and angles
    void Set12Bounds(bool useCurrentGeom);
    //! \brief Update the upper and lower bounds based on bonded angles
    //! @param useCurrentGeom Whether to use the current bond distances and angles
    void Set13Bounds(bool useCurrentGeom);
    //! \brief Set the upper and lower bounds for aromatic ring cycles
    void SetAromaticRingBounds();
    //! \brief Update the upper and lower bounds based on bonded torsions
    void Set14Bounds();
    //! \brief Update the upper and lower bounds based on 1-5 connections
    void Set15Bounds();
    //! \returns 0 if not in the same ring, or the ring size otherwise
    int AreInSameRing(OBAtom *a, OBAtom *b);
    //! \brief Self-consistently smooth the bounds matrix using the triangle inequality
    void TriangleSmooth();
    //! \brief Set the lower bounds to retain VdW distances after all 1-X bounds are set
    void SetLowerBounds();

    //! \return The specified cis/trans stereo configuration for this bond. NULL if not specified
    OBCisTransStereo *GetCisTransStereo(OBBond *bond);

    //! \brief Use OBBuilder to attempt to correct stereo constraints
    void CorrectStereoConstraints(double scale = 1.0);
    //! \brief Check that the double bond and atom stereo constraints are met
    //! \return True if all constraints are valid
    bool CheckStereoConstraints();

    //! \return True if the bounds are met
    bool CheckBounds();
  };
  class DistGeomFunc {
    OBDistanceGeometry* const owner;
    public:
      DistGeomFunc(OBDistanceGeometry* owner) : owner(owner) {}
      double operator() (const Eigen::VectorXd& x, Eigen::VectorXd& grad);
  };

  class DistGeomFunc4D {
    OBDistanceGeometry* const owner;
    public:
      DistGeomFunc4D(OBDistanceGeometry* owner) : owner(owner) {}
      double operator() (const Eigen::VectorXd& x, Eigen::VectorXd& grad);
  };
}

#endif

#endif

//! \file distgeom.h
//! \brief Distance Geometry generation and sampling

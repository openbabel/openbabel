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

#ifndef OBAPI
  #define OBAPI
#endif

#ifdef HAVE_EIGEN

#include <Eigen/Core>

namespace OpenBabel {

  class DistanceGeometryPrivate;

  class OBAPI OBDistanceGeometry {
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
     */
    bool Setup(const OBMol &mol, bool useCurrentGeom = false);

    Eigen::MatrixXf GetBoundsMatrix();
    bool SetBoundsMatrix(const Eigen::MatrixXf bounds);

    void AddConformer();
    void GetConformers(OBMol &mol);

  private:
    OBMol                     _mol;
    DistanceGeometryPrivate  *_d;    //!< Internal private data, including bounds matrix

    void SetDefaultBounds();
    void Set12Bounds(bool useCurrentGeom);
    void Set13Bounds(bool useCurrentGeom);
    void Set14Bounds();
    void Set15Bounds();
    //! \returns 0 if not in the same ring, or the ring size otherwise
    int AreInSameRing(OBAtom *a, OBAtom *b);
    //! \brief Self-consistently smooth the bounds matrix using the triangle inequality
    void TriangleSmooth(int iterations = 5);
    bool CheckChiralConstraints();
  };


}

#endif

#endif

//! \file distgeom.h
//! \brief Distance Geometry generation and sampling

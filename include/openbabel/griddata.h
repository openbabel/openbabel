/**********************************************************************
griddata.h - Store grids of data linked to a molecule (e.g. Gaussian cube)

// Molekel - Molecular Visualization Program
// Copyright (C) 2006, 2007 Swiss National Supercomputing Centre (CSCS)

Some Portions Copyright (C) 2008 by Marcus D. Hanwell

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

#ifndef OB_GRIDDATA_H
#define OB_GRIDDATA_H

#include <openbabel/babelconfig.h>
#include <openbabel/math/vector3.h>
#include <openbabel/generic.h>

#include <vector>

namespace OpenBabel {

  // Uses private data pointer for ABI compatibility
  // http://techbase.kde.org/Policies/Library_Code_Policy#D-Pointers
  class GridDataPrivate;

  // Class to store values for generic (non axis aligned) grids like Gaussian cube files
  // Class introduction in griddata.cpp
  class OBAPI OBGridData : public OBGenericData
  {
  public:
    /// Constructor.
    OBGridData();

    /// Destructor
    ~OBGridData();

    /// Unit of measure for grid spacings.
    typedef enum { BOHR, ANGSTROM, OTHER } Unit;

    //! \name Property Methods
    //@{
    /// \return the three axes parallel to the grid edges. The
    /// length of the returned vector is the step along that
    /// direction.
    void GetAxes(double x[3], double y[3], double z[3]) const;
    //! \return X axis direction.
    vector3 GetXAxis() const;
    //! \return Y axis direction.
    vector3 GetYAxis() const;
    //! \return Z axis direction.
    vector3 GetZAxis() const;
    /// \return the three axes parallel to the grid edges
    void GetAxes(vector3 &v1, vector3 &v2, vector3 &v3) const;
    /// \return number of points along the three axes parallel to the grid edges.
    void GetNumberOfPoints(int &nx, int &ny, int &nz) const;
    /// \return total number of points in the grid.
    int GetNumberOfPoints() const;
    /// \return number of points along the three axes parallel to the grid edges.
    void GetNumberOfSteps(int steps[3]) const;
    /// \return grid values as a vector of doubles.
    std::vector< double > GetValues() const;
    /// \return the value at position i, j, k in the grid.
    double GetValue(int i, int j, int k) const;
    /// \return the value at a position in the grid (by interpolation)
    double GetValue (vector3 pos) const;
    /// \return the unit of measure for grid spacings.
    Unit GetUnit() const;
    /// \return the minimum value of all points in the grid.
    double GetMinValue() const;
    /// \return the maximum value of all points in the grid.
    double GetMaxValue() const;
    /// \return vector3 of the origin (i.e., the minimum x, y, and z coords of the grid).
    vector3 GetOriginVector() const;
    /// \param o set to the origin (i.e., the minimum x, y, and z coords of the grid).
    /// \deprecated Will be removed.
    /// \sa GetOriginVector()
    void GetOriginVector(double o[3]) const;
    /// \return The maximum point in the grid.
    vector3 GetMaxVector() const;
    /// \return the unrestricted flag.
    bool GetUnrestricted() const;
    /// \return the number of symmetries.
    int GetNumSymmetries() const;
    //@}


    //! \name Modification Methods
    //@{
    /// Set number of points along the three axes.
    void SetNumberOfPoints(int nx, int ny, int nz);
    /// Set the limits (i.e., the origin point and the axes)
    /// NOTE: You must set the number of points first,
    ///       with SetNumberOfPoints
    ///       so the grid spacing can be calculated
    void SetLimits(const vector3 &origin, const vector3 &x, const vector3 &y,
                   const vector3 &z);
    /// \deprecated Will be removed.
    /// \sa SetLimits(const vector3 &origin, const vector3 &x, const vector3 &y, const vector3 &z)
    void SetLimits(const double origin[3], const double x[3], const double y[3],
                   const double z[3]);
    /// Set an individual value, grid must have been initialised
    bool SetValue(int i, int j, int k, double val);
    /// Set the values, this vector must match the dimensions of the grid
    void SetValues(const std::vector< double >& v);
    /// Set the unit of measure
    void SetUnit(Unit u);
    /// Set the unrestricted flag
    void SetUnrestricted(bool u);
    /// Set the number of symmetries
    void SetNumSymmetries(int s);
    //@}

  private:
    GridDataPrivate *const d;

  };

} // end namespace

#endif /*OBGRIDDATA_H_*/

//! \file griddata.h
//! \brief OBGenericData class to connect numeric grids (e.g., orbitals, electrostatic potential) to molecules

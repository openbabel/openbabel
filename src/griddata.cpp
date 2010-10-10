/**********************************************************************
griddata.cpp - Store grids of data linked to a molecule (e.g. Gaussian cube)

// Molekel - Molecular Visualization Program
// Copyright (C) 2006, 2007 Swiss National Supercomputing Centre (CSCS)

 Some Portions Copyright (c) 2007 by Geoffrey R. Hutchison
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

#include <openbabel/babelconfig.h>

#include <openbabel/griddata.h>
#include <openbabel/mol.h>
#include <openbabel/grid.h>

#include <vector>
#include <algorithm>
#include <limits>

using namespace std;

namespace OpenBabel {

  class GridDataPrivate {
  public:
    GridDataPrivate() {    }

    OBFloatGrid  floatGrid;
    OBGridData::Unit _unit;

    double           _max;
    double           _min;

    bool             _unrestricted;
    int              _symmetries;
  };

  /** \class OBGridData griddata.h <openbabel/griddata.h>
    \brief Store values for numeric grids such as orbitals or electrostatic potential
    \since version 2.2
    \sa OBFloatGrid

    OBGridData facilitates attaching grids and cubes to molecular data. A "grid" is
    data representing some function f(x,y,z), such as a molecule's electrostatic potential
    or molecular orbitals. This need not be a "cube" even though this file format from Gaussian
    is frequently used. Axes need not be identical, and indeed do not need to be orthogonal.

    Open Babel supports reading several types of grid file formats, including Gaussian cube,
    and OpenDX. The latter is notably used by the APBS program for numeric evaluation of molecular
    and protein electrostatic potential.

    \code
    OBGridData *gd = new OBGridData;
    gd->SetAttribute("Example Grid"); // the title of the grid -- e.g., for user display
    vector<int> voxels(3); // the number of voxels in each direction
    vector3 origin; // the beginning x, y, z coordinate of the grid
    vector<vector3> axes; // the xyz displacements for each of the grid axes
    ...
    gd->SetNumberOfPoints(voxels[0], voxels[1], voxels[2]);
    gd->SetLimits(origin, axes[0], axes[1], axes[2]);
    gd->SetUnit(OBGridData::ANGSTROM);
    gd->SetOrigin(fileformatInput); // i.e., is this data from a file or determined by Open Babel

    for (int k = 0; k < voxels[2]; ++k)
      for (int j = 0; j < voxels[1]; ++j)
        for (int i = 0; i < voxels[0]; ++i)
          {
            gd->SetValue(i, j, k,
                         grid[k*voxels[0]*voxels[1] + j*voxels[0] + i]);
          }

    mol->SetData(gd);
    \endcode

    \code
    if (mol->HasData(OBGenericDataType::GridData)) {
      vector<OBGenericData*> grids = mol->GetAllData(OBGenericDataType::GridData)
      // Output the name of the grid
      if (grids[0] != NULL)
        cout << grids[0]->GetAttribute();
    }
    \endcode

  */

  OBGridData::OBGridData() : OBGenericData("GridData", OBGenericDataType::GridData),
    d(new GridDataPrivate)
  {
  }

  OBGridData::~OBGridData()
  {
    delete d;
  }

  void OBGridData::GetAxes( double x[3], double y[3], double z[3] ) const
  {
    vector3 v1, v2, v3;
    v1 = d->floatGrid.GetXAxis();
    v2 = d->floatGrid.GetYAxis();
    v3 = d->floatGrid.GetZAxis();

    x[0] = v1.x(); x[1] = v1.y(), x[2] = v1.z();
    y[0] = v2.x(); y[1] = v2.y(), y[2] = v2.z();
    z[0] = v3.x(); z[1] = v3.y(), z[2] = v3.z();
  }

  vector3 OBGridData::GetXAxis() const
  {
    return d->floatGrid.GetXAxis();
  }

  vector3 OBGridData::GetYAxis() const
  {
    return d->floatGrid.GetYAxis();
  }

  vector3 OBGridData::GetZAxis() const
  {
    return d->floatGrid.GetZAxis();
  }

  void OBGridData::GetAxes( vector3 &v1, vector3 &v2, vector3 &v3 ) const
  {
    v1 = d->floatGrid.GetXAxis();
    v2 = d->floatGrid.GetYAxis();
    v3 = d->floatGrid.GetZAxis();
  }

  void OBGridData::GetNumberOfPoints( int &nx, int &ny, int &nz) const
  {
    nx = d->floatGrid.GetXdim();
    ny = d->floatGrid.GetYdim();
    nz = d->floatGrid.GetZdim();
  }

  int OBGridData::GetNumberOfPoints() const
  {
    return d->floatGrid.GetXdim() * d->floatGrid.GetYdim() * d->floatGrid.GetZdim();
  }

  void OBGridData::GetNumberOfSteps( int steps[ 3 ] ) const
  {
    steps[0] = d->floatGrid.GetXdim() - 1;
    steps[1] = d->floatGrid.GetYdim() - 1;
    steps[2] = d->floatGrid.GetZdim() - 1;
  }

  std::vector< double > OBGridData::GetValues() const
  {
    return d->floatGrid.GetDataVector();
  }

  double OBGridData::GetValue( int i, int j, int k ) const
  {
    return d->floatGrid.GetValue(i, j, k);
  }

  double OBGridData::GetValue(vector3 pos) const
  {
    return d->floatGrid.Interpolate(pos.x(), pos.y(), pos.z());
  }

  OBGridData::Unit OBGridData::GetUnit() const
  {
    return d->_unit;
  }

  double OBGridData::GetMinValue() const
  {
    return d->_min;
  }

  double OBGridData::GetMaxValue() const
  {
    return d->_max;
  }

  void OBGridData::GetOriginVector( double o[ 3 ] ) const
  {
    d->floatGrid.GetMin(o);
  }

  vector3 OBGridData::GetOriginVector() const
  {
    return d->floatGrid.GetMin();
  }

  vector3 OBGridData::GetMaxVector() const
  {
    return d->floatGrid.GetMax();
  }

  bool OBGridData::GetUnrestricted() const
  {
    return d->_unrestricted;
  }

  int OBGridData::GetNumSymmetries() const
  {
    return d->_symmetries;
  }

  void OBGridData::SetUnrestricted( bool u )
  {
    d->_unrestricted = u;
  }

  void OBGridData::SetNumSymmetries( int s )
  {
    d->_symmetries = s;
  }

  void OBGridData::SetNumberOfPoints( int nx, int ny, int nz )
  {
    d->floatGrid.SetNumberOfPoints(nx, ny, nz);
  }

  void OBGridData::SetLimits(const double origin [3], const double x[3],
                             const double y[3], const double z[3])
  {
    d->floatGrid.SetLimits(origin, x, y, z);
  }

  void OBGridData::SetLimits(const vector3 &origin, const vector3 &x,
                             const vector3 &y, const vector3 &z)
  {
    d->floatGrid.SetLimits(origin, x, y, z);
  }

  bool OBGridData::SetValue(int i, int j, int k, double val)
  {
    return d->floatGrid.SetValue(i, j, k, val);
  }

  void OBGridData::SetValues( const std::vector< double >& v )
  {
    d->floatGrid.SetVals(v);
    d->_min = *std::min_element( v.begin(), v.end() );
    d->_max = *std::max_element( v.begin(), v.end() );
  }

  void OBGridData::SetUnit( OBGridData::Unit u )
  {
    d->_unit = u;
  }

} // end namespace

//! \file griddata.cpp
//! \brief OBGenericData class to connect numeric grids (e.g., orbitals, electrostatic potential) to molecules

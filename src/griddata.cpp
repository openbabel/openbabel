/**********************************************************************
griddata.cpp - Store grids of data linked to a molecule (e.g. Gaussian cube)

// Molekel - Molecular Visualization Program
// Copyright (C) 2006, 2007 Swiss National Supercomputing Centre (CSCS)

 Some Portions Copyright (c) 2007 by Geoffrey R. Hutchison
 Some Portions Copyright (C) 2008 by Marcus D. Hanwell

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

    OBFloatGrid _floatGrid;
    OBGridData::Unit _unit;

    double           _max;
    double           _min;
  };

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
    v1 = d->_floatGrid.GetXAxis();
    v2 = d->_floatGrid.GetYAxis();
    v3 = d->_floatGrid.GetZAxis();

    x[0] = v1.x(); x[1] = v1.y(), x[2] = v1.z();
    y[0] = v2.x(); y[1] = v2.y(), y[2] = v2.z();
    z[0] = v3.x(); z[1] = v3.y(), z[2] = v3.z();
  }

  void OBGridData::GetAxes( vector3 &v1, vector3 &v2, vector3 &v3 ) const
  {
    v1 = d->_floatGrid.GetXAxis();
    v2 = d->_floatGrid.GetYAxis();
    v3 = d->_floatGrid.GetZAxis();
  }

  void OBGridData::GetNumberOfPoints( int &nx, int &ny, int &nz) const
  {
    nx = d->_floatGrid.GetXdim();
    ny = d->_floatGrid.GetYdim();
    nz = d->_floatGrid.GetZdim();
  }

  void OBGridData::GetNumberOfSteps( int steps[ 3 ] ) const
  {
    steps[0] = d->_floatGrid.GetXdim() - 1;
    steps[1] = d->_floatGrid.GetYdim() - 1;
    steps[2] = d->_floatGrid.GetZdim() - 1;
  }

  std::vector< double > OBGridData::GetValues() const
  {
    return d->_floatGrid.GetDataVector();
  }

  double OBGridData::GetValue( int i, int j, int k ) const
  {
    return d->_floatGrid.GetValue(i, j, k);
  }

  double OBGridData::GetValue(vector3 pos) const
  {
    return d->_floatGrid.Interpolate(pos.x(), pos.y(), pos.z());
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
    d->_floatGrid.GetMin(o);
  }

  vector3 OBGridData::GetOriginVector() const
  {
    return d->_floatGrid.GetMin();
  }

  void OBGridData::SetNumberOfPoints( int nx, int ny, int nz )
  {
    d->_floatGrid.SetNumberOfPoints(nx, ny, nz);
  }

  void OBGridData::SetLimits(double origin [ 3 ], double x[ 3 ], double y[ 3 ], double z[ 3 ] )
  {
    d->_floatGrid.SetLimits(origin, x, y, z);
  }

  void OBGridData::SetLimits(vector3 &origin, vector3 &x, vector3 &y, vector3 &z)
  {
    d->_floatGrid.SetLimits(origin, x, y, z);
  }

  bool OBGridData::SetValue(int i, int j, int k, double val)
  {
    return d->_floatGrid.SetValue(i, j, k, val);
  }

  void OBGridData::SetValues( const std::vector< double >& v )
  {
    d->_floatGrid.SetVals(v);
    d->_min = *std::min_element( v.begin(), v.end() );
    d->_max = *std::max_element( v.begin(), v.end() );
  }

  void OBGridData::SetUnit( OBGridData::Unit u )
  {
    d->_unit = u;
  }

} // end namespace


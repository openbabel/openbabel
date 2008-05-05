//
// Molekel - Molecular Visualization Program
// Copyright (C) 2006, 2007 Swiss National Supercomputing Centre (CSCS)
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
// MA  02110-1301, USA.
//
// $Author$
// $Date$
// $Revision$
//

#ifndef OBT41DATA_H_
#define OBT41DATA_H_

#include <openbabel/obmolecformat.h>
#include <openbabel/base.h>
#include <vector>
#include <algorithm>
#include <limits>
#include <cassert>
#include <string>
#include <map>
#include <vector>

namespace OpenBabel {

/// Class to store values read from ADF Tape41 files.
class OBT41Data : public OBGenericData
{
public:
    /// Constructor assigns the values of type and attr protected data
    /// This values will be accessed through the GetDataType, HasData methods.
    OBT41Data() : OpenBabel::OBGenericData()
    {
        _type = OpenBabel::OBGenericDataType::CustomData1;
        _attr = "T41Data";
    }

    /// Returns the three axes parallel to the grid edges the
    /// length of the returned vector is the step along that
    /// direction.
    void GetAxes( double x[ 3 ], double y[ 3 ], double z[ 3 ] ) const
    {
        x[ 0 ] = xAxis_[ 0 ]; x[ 1 ] = xAxis_[ 1 ]; x[ 2 ] = xAxis_[ 2 ];
        y[ 0 ] = yAxis_[ 0 ]; y[ 1 ] = yAxis_[ 1 ]; y[ 2 ] = yAxis_[ 2 ];
        z[ 0 ] = zAxis_[ 0 ]; z[ 1 ] = zAxis_[ 1 ]; z[ 2 ] = zAxis_[ 2 ];
    }

    /// Return number of points along the three axes paralled to the
    /// grid edges.
    void GetNumberOfPoints( int& nx, int& ny, int& nz ) const
    {
        nx = nx_;
        ny = ny_;
        nz = nz_;
    }

    /// Return total number of points in grid.
    int GetNumberOfPoints() const
    {
        return nx_ * ny_ * nz_;
    }


    /// Return number of points along the three axes paralled to the
    /// grid edges.
    void GetNumberOfSteps( int steps[ 3 ] ) const
    {
        steps[ 0 ] = nx_ - 1;
        steps[ 1 ] = ny_ - 1;
        steps[ 2 ] = nz_ - 1;
    }

    /// Return grid values as an array of doubles.
    const std::vector< double >& GetValues( const std::string& key ) const
    {
        assert( values_.find( key ) != values_.end() );
        return values_.find( key )->second;
    }

    /// Returns point at position i, j, k in the grid.
    double GetValue( const std::string& key, int i, int j, int k ) const
    {
        assert( values_.find( key ) != values_.end() );
        const int idx = ComputeIndex( i, j, k );
        return values_.find( key )->second[ idx ];
    }

    /// Returns min value.
    double GetMinValue( const std::string& key ) const
    {
        assert( min_.find( key ) != min_.end() );
        return min_.find( key )->second;
    }

    /// Returns max value.
    double GetMaxValue( const std::string& key ) const
    {
        assert( max_.find( key ) != max_.end() );
        return max_.find( key )->second;
    }

    /// Returns origin.
    const double* GetStartPoint() const { return &startPoint_[ 0 ]; }

    /// Returns origin.
    void GetStartPoint( double sp[ 3 ] ) const
    {
        sp[ 0 ] = startPoint_[ 0 ];
        sp[ 1 ] = startPoint_[ 1 ];
        sp[ 2 ] = startPoint_[ 2 ];
    }

    /// Grid unrestricted flag.
    bool GetUnrestricted() const { return unrestricted_; }

    /// Number of symmetries.
    int GetNumSymmetries() const { return numSymmetries_; }


    /// Return grid labels
    typedef std::vector< std::string > GridLabels;
    GridLabels GetGridLabels() const
    {
        GridLabels labels;
        labels.reserve( values_.size() );
        typedef std::map< std::string, std::vector< double > > Grid;
        Grid::const_iterator i = values_.begin();
        const Grid::const_iterator end = values_.end();
        for( ; i != end; ++i ) labels.push_back( i->first );
        return labels;
    }


    /// Reserve data in value vector.
    void Reserve( const std::string& key, int size ) { values_[ key ].reserve( size ); }

    /// Append value to value vector.
    void AddValue( const std::string& key, double v )
    {
        values_[ key ].push_back( v );
        if( v < min_[ key ] ) min_[ key ] = v;
        if( v > max_[ key ] ) max_[ key ] = v;
    }

    /// Set number of points along the three axes.
    void SetNumberOfPoints( int nx, int ny, int nz )
    {
        nx_ = nx;
        ny_ = ny;
        nz_ = nz;
    }

    /// Set grid axes.
    void SetAxes( double x[ 3 ], double y[ 3 ], double z[ 3 ] )
    {
        xAxis_[ 0 ] = x[ 0 ]; xAxis_[ 1 ] = x[ 1 ]; xAxis_[ 2 ] = x[ 2 ];
        yAxis_[ 0 ] = y[ 0 ]; yAxis_[ 1 ] = y[ 1 ]; yAxis_[ 2 ] = y[ 2 ];
        zAxis_[ 0 ] = z[ 0 ]; zAxis_[ 1 ] = z[ 1 ]; zAxis_[ 2 ] = z[ 2 ];
    }

    /// Set value vector.
    void SetValues( const std::string& key, const std::vector< double >& v )
    {
        values_[ key ] = v;
        min_[ key ] = *std::min_element( values_[ key ].begin(), values_[ key ].end() );
        max_[ key ] = *std::max_element( values_[ key ].begin(), values_[ key ].end() );
    }

    /// Set start point.
    void SetStartPoint( double sp[ 3 ] )
    {
        startPoint_[ 0 ] = sp[ 0 ];
        startPoint_[ 1 ] = sp[ 1 ];
        startPoint_[ 2 ] = sp[ 2 ];
    }

    /// Grid unrestricted flag.
    void SetUnrestricted( bool u ) { unrestricted_ = u; }

    /// Number of symmetries.
    void SetNumSymmetries( int s ) { numSymmetries_ = s; }

private:

    // @{ Axes
    double xAxis_[ 3 ];
    double yAxis_[ 3 ];
    double zAxis_[ 3 ];
    // @}
    //@{ Steps along axes
    int nx_;
    int ny_;
    int nz_;
    // @}
    /// Grid unrestricted flag.
    bool unrestricted_;
    /// Number of symmetries.
    int numSymmetries_;
    /// Grid values.
    std::map< std::string, std::vector< double > > values_;
    /// Start point.
    double startPoint_[ 3 ];
    /// Min value.
    std::map< std::string, double > min_;
    /// Max value.
    std::map< std::string, double > max_;
    /// Return vector index given i, j, k grid coordinates.
    int ComputeIndex( int i, int j, int k ) const
    {
        assert( i >= 0 && i < nx_ &&
                j >= 0 && j < ny_ &&
                k >= 0 && k < nz_ &&
                "Grid index out of bounds" );
        return  i + nx_ *( j +  ny_ * k );
    }
};

} // end namespace OpenBabel

#endif /*OBT41DATA_H_*/

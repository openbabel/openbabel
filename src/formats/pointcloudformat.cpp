//
// Point Cloud Plugin for Open Babel
// Copyright (C) 2015 M J Harvey,
// Acellera Ltd
// m.j.harvey ( at ) acellera.com
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


#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/obconversion.h>
#include <openbabel/obmolecformat.h>
#include <openbabel/elements.h>

#include <openbabel/obiter.h>
#include <openbabel/data.h>
#include <openbabel/bitvec.h>

#include <iostream>

#include <stdlib.h>
#include <math.h>
#include <cstdlib>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


using namespace std;

namespace OpenBabel
{


  //==============================================================================
  /// Class to output a point cloud on a surface around a molecule

  class PointCloudFormat : public OpenBabel::OBMoleculeFormat
  {
    public:
      PointCloudFormat()
      {
        OpenBabel::OBConversion::RegisterFormat( "pointcloud", this );
        /*
        OBConversion::RegisterOptionParam("r", this, 1, OBConversion::OUTOPTIONS);
        OBConversion::RegisterOptionParam("d", this, 1, OBConversion::OUTOPTIONS);
        OBConversion::RegisterOptionParam("p", this, 1, OBConversion::OUTOPTIONS);
        OBConversion::RegisterOptionParam("x", this, 1, OBConversion::OUTOPTIONS);
        */
      }

      /// Return description.
      virtual const char* Description() //required
      {
        return
          "Point cloud on VDW surface\n"
          "Generates a point cloud on the VDW surface around the molecule\n\n"

          "The surface location is calculated by adding the probe atom radius\n"
          "(if specified) to the Van der Waal radius of the particular atom multipled\n"
          "by the specified multiple (1.0 if unspecified)."

          "Output is a list of {x,y,z} tuples in Angstrom. Alternatively, if the ``x``\n"
          "option is specified, the :ref:`XYZ_cartesian_coordinates_format` is used\n"
          "instead.\n\n"

          "Write Options, e.g. -xx\n"
          "  r <radii> create a surface for each VDS radius (default 1.0)\n"
          "        A comma-separated list of VDW radius multiples\n"
          "  d <densities> for each surface, specify the point density (default 1.0 Angstrom^2)\n"
          "        A comma-separated list of densities\n"
          "  p <radius> radius of the probe atom in Angstrom (default 0.0)\n"
          "  x output in xyz format\n\n";
      }

      /// Return a specification url, not really a specification since
      /// I couldn't find it but close enough.
      virtual const char* SpecificationURL()
      {
        return "N/A";
      }

      /// Return MIME type, NULL in this case.
      virtual const char* GetMIMEType() { return "chemical/x-pointcloud"; };

      /// Return read/write flag: read only.
      virtual unsigned int Flags()
      {
        return WRITEONEONLY | NOTREADABLE;
      };

      /// Skip to object: used for multi-object file formats.
      virtual int SkipObjects( int n, OpenBabel::OBConversion* pConv ) { return 0; }

      /// Read: always return false.
      virtual bool ReadMolecule( OpenBabel::OBBase*, OpenBabel::OBConversion* )
      {
        return false;
      }

      /// Write.
      virtual bool WriteMolecule( OpenBabel::OBBase* , OpenBabel::OBConversion* );
  };


  // Global variable used to register PointCloud format.
  PointCloudFormat thePointCloudFormat;



  //==============================================================================

  double rv( void ){
    return ((double)std::rand()) / RAND_MAX; // C++11 provides a better <random> but I think this will do here
  }

  vector3 surface_point( double cx, double cy, double cz, double radius ) {
    double  theta = rv() * 2. * M_PI;
    double  phi   = acos( 2. * rv() - 1. );

    double  x = radius * cos(theta) * sin(phi);
    double  y = radius * sin(theta) * sin(phi);
    double  z = radius * cos(phi);

    return vector3( cx+x, cy+y, cz+z );
  }

  // Only add a point if it's more than densit_r away from all others
  bool conditional_add( vector<vector3> &list, vector3 point, double density_r ) {

    double density_r2 = density_r * density_r;
    for( std::vector<vector3>::iterator it = list.begin(); it != list.end(); ++it ) {
      vector3 r = *it - point;
      double r2 = dot(r,r);
      if( r2 < density_r2 ) {
        return false;
      }
    }
    list.push_back( point );
    return true;
  }

  //------------------------------------------------------------------------------
  bool PointCloudFormat::WriteMolecule( OBBase* pOb, OBConversion* pConv )
  {
    OBMol* pmol = dynamic_cast< OBMol* >(pOb);
    if( pmol == 0 ) return false;

    ostream& os = *pConv->GetOutStream();

    const char *radius_list_str  = NULL;
    const char *density_list_str = NULL;
    double probe_radius = 0.;
    bool format_xyz = false;

    if( pConv->IsOption( "r" ) ) {
      radius_list_str = pConv->IsOption("r", OBConversion::OUTOPTIONS);
    }
    if( pConv->IsOption( "d" ) ) {
      density_list_str = pConv->IsOption("d", OBConversion::OUTOPTIONS);
    }
    if( pConv->IsOption( "p" ) ) {
      probe_radius = atof( pConv->IsOption("p", OBConversion::OUTOPTIONS) );
      if( !isfinite(probe_radius) || probe_radius < 0. ) { probe_radius = 0.; }
    }
    if( pConv->IsOption( "x" ) ) {
      format_xyz = true;
    }


    int total_points = 0;
    std::srand(0); // let's be deterministic(FSVO)

    vector< vector3 > filtered_points;

    vector<double> radius_mult_list;
    vector<double> density_list;

    // Extract the any lists of radii and point densities
    if( radius_list_str ) {
      char*a = strdup( radius_list_str );
      const char * x = strtok( a, "," );
      while( x ) {
        double d = atof(x);
        if( isfinite(d) && d>0. ) { radius_mult_list.push_back( d ); }
        x = strtok( NULL, "," );
      }
      free(a);
    }
    if( density_list_str ) {
      char*a = strdup( density_list_str );
      const char * x = strtok( a, "," );
      while( x ) {
        double d = atof(x);
        if( isfinite(d) && d>0. ) { density_list.push_back( d ); }
        x = strtok( NULL, "," );
      }
      free(a);
    }
    if( radius_mult_list.size() == 0 ) { radius_mult_list.push_back( 1.0 ); }
    while( density_list.size() < radius_mult_list.size()) {
      density_list.push_back(1.0);
    }

    // foreach radius generate a point cloud
    for( int j = 0; j < radius_mult_list.size(); j++ ) {

      double radius_mult  = radius_mult_list[j];
      double density      = density_list[j];
      double density_r    = sqrt ( density / M_PI );

      FOR_ATOMS_OF_MOL( a, *pmol )
      {
        vector<vector3> pt;


        // for each atom generate a cloud of points
        // on the atom-centred spherical surface r_VDW * radius_multiplier
        // ensuring each point is >= density_r from all others
        const double* c = a->GetCoordinate();
        double vdwrad   = probe_radius + ( OBElements::GetVdwRad( a->GetAtomicNum() ) * radius_mult );

        // estimate # of points by dividing area of VDW sphere
        // by the required density
        int estimate_number_of_points = (int)  ( 0.6 * (( 4. * M_PI * M_PI * vdwrad * vdwrad ) / ( density )) );
        int count = 0;
        while ( count < estimate_number_of_points ) {
          vector3 p = surface_point( c[0], c[1], c[2], vdwrad );
          if( conditional_add( pt, p, density_r ) ) {
            count++;
            total_points++;
          }
        }

        // now cull any points on that surface that are within r_VDW * radius_multiplier
        // of any other atom
        for( std::vector< vector3 >::iterator it = pt.begin(); it != pt.end(); ++it ) {
          bool exclude = false;
          FOR_ATOMS_OF_MOL( a, *pmol )
          {
            const double* c = a->GetCoordinate();
            double vdwrad   =  probe_radius + ( OBElements::GetVdwRad( a->GetAtomicNum() ) * radius_mult );
            vdwrad *= vdwrad;
            vector3 r = *it - vector3( c[0], c[1], c[2] );
            double r2 = r[0]*r[0] + r[1]*r[1] + r[2] * r[2];
            if( r2 < vdwrad ) { exclude = true; break; }
          }
          if( !exclude ) {
            filtered_points.push_back( *it );
          }
        }
      }
    }

    // Output the final list of points in XYZ format
    if( format_xyz ) {
      os << filtered_points.size() << "\n\n";
    }

    for( std::vector<vector3>::iterator it2 = filtered_points.begin(); it2 != filtered_points.end(); ++it2 ) {
      if( format_xyz ) {
        os << "Xx\t";
      }
      os << (*it2)[0] << "\t" << (*it2)[1] << "\t" << (*it2)[2] << "\n";
    }

    os.flush();
    return true;
  }

}

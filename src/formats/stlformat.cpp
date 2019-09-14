//
// STL Plugin for Open Babel
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
#include <openbabel/elements.h>

#include <openbabel/obconversion.h>
#include <openbabel/obmolecformat.h>

#include <openbabel/obiter.h>
#include <openbabel/data.h>

#include <iostream>

#include <stdlib.h>
#include <math.h>

// uint8_t and uint16_t are not defined for earlier versions of msvc
#if defined(_MSC_VER) && _MSC_VER <= 1600
  typedef unsigned __int8 uint8_t;
  typedef unsigned __int16 uint16_t;
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


using namespace std;

namespace OpenBabel
{


  //==============================================================================
  /// Class to output a point cloud on a surface around a molecule

  class STLFormat : public OpenBabel::OBMoleculeFormat
  {
    public:
      STLFormat()
      {
        OpenBabel::OBConversion::RegisterFormat( "stl", this );
        /*
        OBConversion::RegisterOptionParam("p", this, 1, OBConversion::OUTOPTIONS);
        OBConversion::RegisterOptionParam("s", this, 1, OBConversion::OUTOPTIONS);
        OBConversion::RegisterOptionParam("c", this, 1, OBConversion::OUTOPTIONS);
        */
      }

      /// Return description.
      virtual const char* Description() //required
      {
        return
          "STL 3D-printing format\n"
          "The STereoLithography format developed by 3D Systems\n\n"
          "Write Options, e.g. -xc\n"
          "  p <radius> radius for probe particle (default 0.0 A)\n"
          "  s <scale> scale-factor for VDW radius (default 1.0 A)\n"
          "  c add CPK colours\n\n";
      }

      virtual const char* SpecificationURL()
      {
        return "http://www.ennex.com/~fabbers/StL.asp";
      }

      /// Return MIME type, NULL in this case.
      virtual const char* GetMIMEType() { return "application/sla"; };

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


  // Global variable used to register STL format.
  STLFormat theSTLFormat;


  // Produces triangles over the surface of a sphere


  struct Triangle {
    vector3 a, b, c;
    uint16_t col;
    Triangle(vector3 x, vector3 y, vector3 z, uint16_t colour ) {
      a=x;
      b=y;
      c=z;
      col = colour;
    }
  };

  void map_sphere ( vector<Triangle> &triangles, vector3 origin, double r, uint16_t col )
  {
    vector<vector3> points;

    int  longitude_steps = 144;
    int  latitude_steps  = longitude_steps / 2;
    double theta =  ( 2*M_PI / longitude_steps );
    int p2 = longitude_steps / 2;
    int r2 = latitude_steps / 2;
    for(int y = -r2; y < r2; ++y) {
      double cy = cos(y*theta);
      double cy1 = cos((y+1)*theta);
      double sy = sin(y*theta);
      double sy1 = sin((y+1)*theta);

      for(int i = -p2; i < p2; ++i) {
        double ci = cos(i*theta);
        double si = sin(i*theta);
        points.push_back(vector3( origin[0] + r * ci*cy , origin[1] + r*sy , origin[2] + r * si*cy));
        points.push_back(vector3( origin[0] + r * ci*cy1, origin[1] + r*sy1, origin[2] + r * si*cy1));
      }
    }

    for( int i=0 ; i < points.size()-2; i++ ) {
      // Order the points to obey the 'right-hand rule', so normals all face outwards
      if( i%2 == 0 ) {
        triangles.push_back( Triangle( points[i], points[i+1], points[i+2], col ) );
      }
      else {
        triangles.push_back( Triangle( points[i+2], points[i+1], points[i], col ) );
      }
    }
  }

  static uint16_t stl_colour( int atomicnum ) {
    // CPK colouring
    // STL format: 5 bits per col, bit 15 set
#define COL( R,G,B ) ( 0x800 | (((R)&0x1F)<<10) | (((G)&0x1F)<<5) | ((B)&0x1F) );

    switch( atomicnum ) {
      case 1:  return COL( 0x1F, 0x1F, 0x1F ); // H, WHITE
      case 6:  return COL( 0x04, 0x04, 0x04 ); // C, (almost) BLACK
      case 7:  return COL( 0x12, 0x1A, 0x1F );  // N, SKY BLUE
      case 8:  return COL( 0x1F, 0x00, 0x00 );  // O, RED
      case 16: return COL( 0x1F, 0x1F, 0x00 );  // S, YELLOW
      case 15: return COL( 0x1F, 0x00, 0x1F );  // P, PURPLE
      case 9:  return COL( 0x00, 0x1F, 0x00 );  //F, LIGHT GREEN
      case 17: return COL( 0x00, 0x17, 0x00 );  //Cl, MEDIUM GREEN
      case 35: return COL( 0x00, 0x0F, 0x00 );  //Br, MEDIUM DARK GREEN
      case 53: return COL( 0x00, 0x07, 0x00 );   //I, DARK GREEN
      case 26: // Fe
      case 27: // Co
      case 28: // Ni
      case 29:	return COL(0x18, 0x18, 0x18 ); // Cu, Metals Silver
      default:
                return COL( 0x08, 0x08, 0x08 ); // Default, grey

    }
    return 0x8FFF;
  }

  static void output_stl( ostream &os, vector<Triangle> &triangles, bool cpk_colour ) {
    uint8_t  tok1 = 0;
    uint32_t tok2 = 0;

    if( cpk_colour ) {
      tok1=0xff;
      os.write( "COLOR=", 6 );
      os.write( (const char*) &tok1, 1);
      os.write( (const char*) &tok1, 1);
      os.write( (const char*) &tok1, 1);
      os.write( (const char*) &tok1, 1);

      tok1=0;
      for( int i=0; i < 70; i++ ) {
        os.write( (const char*) &tok1, 1);
      }
    }
    else {
      for( int i=0; i < 80; i++ ) {
        os.write( (const char*) &tok1, 1);
      }
    }

    tok2 = triangles.size();
    os.write( (const char*) &tok2, sizeof(uint32_t) );

    for(vector<Triangle>::iterator i = triangles.begin(); i != triangles.end(); ++i) { // go through all triangles
      float nx,ny,nz;
      //		nx= (i->a.x + i->b.x + i->c.x) / 3.;
      //		ny= (i->a.y + i->b.y + i->c.y) / 3.;
      //		nz= (i->a.z + i->b.z + i->c.z) / 3.;

      //		float r = sqrt( nx*nx + ny*ny + nz*nz );
      //		nx /= r;
      //		ny /= r;
      //		nz /= r;

      //  Turns out we don't have to specify a normal
      nx = ny = nz = 0.f;

      os.write( (const char*) &nx, sizeof(float) );
      os.write( (const char*) &ny, sizeof(float) );
      os.write( (const char*) &nz, sizeof(float) );

      nx = i->a[0];
      ny = i->a[1];
      nz = i->a[2];

      os.write( (const char*) &nx, sizeof(float) );
      os.write( (const char*) &ny, sizeof(float) );
      os.write( (const char*) &nz, sizeof(float) );

      nx = i->b[0];
      ny = i->b[1];
      nz = i->b[2];

      os.write( (const char*) &nx, sizeof(float) );
      os.write( (const char*) &ny, sizeof(float) );
      os.write( (const char*) &nz, sizeof(float) );

      nx = i->c[0];
      ny = i->c[1];
      nz = i->c[2];

      os.write( (const char*) &nx, sizeof(float) );
      os.write( (const char*) &ny, sizeof(float) );
      os.write( (const char*) &nz, sizeof(float) );

      os.write( (const char*) &(i->col), sizeof(uint16_t) );
    }

    os.flush();
  }




  bool STLFormat::WriteMolecule( OBBase* pOb, OBConversion* pConv )
  {
    OBMol* pmol = dynamic_cast< OBMol* >(pOb);
    if( pmol == 0 ) return false;

    ostream& os = *pConv->GetOutStream();

    double probe_radius = 0.;
    double scale_factor = 1.;
    uint16_t col = 0x0; // no colour colour
    bool cpk_colours = false;

    if( pConv->IsOption( "p" ) ) {
      probe_radius = atof( pConv->IsOption("p", OBConversion::OUTOPTIONS) );
      if( !isfinite(probe_radius) || probe_radius < 0. ) { probe_radius = 0.; }
    }
    if( pConv->IsOption( "s" ) ) {
      probe_radius = atof( pConv->IsOption("s", OBConversion::OUTOPTIONS) );
      if( !isfinite(scale_factor) || scale_factor < 0. ) { scale_factor = 0.; }
    }
    if( pConv->IsOption( "c" ) ) {
      cpk_colours = true;
    }


    vector<Triangle> triangles;
    FOR_ATOMS_OF_MOL(a, *pmol) {
      const double *coord = a->GetCoordinate();
      const double vdwrad = scale_factor * OBElements::GetVdwRad( a->GetAtomicNum() ) + probe_radius;
      if( cpk_colours ) {
        col =  stl_colour(  a->GetAtomicNum() ) ;
      }
      map_sphere( triangles, vector3(coord[0], coord[1], coord[2]), vdwrad, col );
    }

    output_stl ( os, triangles, cpk_colours );
    os.flush();
    return true;
  }

}

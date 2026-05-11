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


#include <openbabel/math/vector3.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/elements.h>

#include <openbabel/obconversion.h>
#include <openbabel/obmolecformat.h>

#include <openbabel/obiter.h>
#include <openbabel/data.h>

#include <iostream>
#include <ostream>
#include <vector>
#include <array>
#include <algorithm>

#include <cstdlib>
#include <cmath>

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
      /**
        * @brief Enumeration of supported non-standard color channel orders for STL.
        *        See https://en.wikipedia.org/wiki/STL_(file_format)#Binary.
        */
      enum class ColourOrder : uint8_t { RedGreenBlue, BlueGreenRed };

    public:
      /**
       * @brief Non-standard STL value for no color for triangle face.
       */
      static constexpr uint16_t NoColor{0};

      /**
       * @brief Default colour order.
       */
      static constexpr ColourOrder DefaultColourOrder{ ColourOrder::BlueGreenRed };

      /**
       * @brief Default value for -xc option. (Unit: A)
       */
      static constexpr double DefaultProbeRadius{0.0};

      /**
       * @brief Default scaling factor for -xs option.
       */
      static constexpr double DefaultScalingFactor{1.0};

      /**
       * @brief Default bond radius for --bond-radius option. (Unit: A)
       */
      static constexpr double DefaultBondRadius{0.5};

      /**
       * @brief Default bond spacing for --bond-spacing option. (Unit: A)
       */
      static constexpr double DefaultBondSpacing{0.1};

    public:
      STLFormat()
      {
        OpenBabel::OBConversion::RegisterFormat( "stl", this );
        /*
        OBConversion::RegisterOptionParam("p", this, 1, OBConversion::OUTOPTIONS);
        OBConversion::RegisterOptionParam("s", this, 1, OBConversion::OUTOPTIONS);
        OBConversion::RegisterOptionParam("c", this, 1, OBConversion::OUTOPTIONS);
        */
        OBConversion::RegisterOptionParam("bond-size", this, 1, OBConversion::GENOPTIONS);
        OBConversion::RegisterOptionParam("bond-radius", this, 1, OBConversion::GENOPTIONS);
        OBConversion::RegisterOptionParam("colour-order-rgb", this, 0, OBConversion::GENOPTIONS);
      }

      /// Return description.
      const char* Description() override  // required
      {
        return
          "STL 3D-printing format\n"
          "The STereoLithography format developed by 3D Systems\n\n"
          "Write Options, e.g. -xc\n"
          "  p <radius> radius for probe particle (default 0.0 A)\n"
          "  s <scale> scale-factor for VDW radius and bond radius/spacing (default 1.0)\n"
          "  c add CPK colours\n"
          "  b draw bonds\n"
          "General options,\n"
          "  --bond-radius <radius> radius for bond visualization (default: 0.5 A)\n"
          "  --bond-spacing <spacing> spacing between bond visualizations (default: 0.1 A)\n"
          "  --colour-order-rgb export colours in red green blue (RGB) order (default: BGR order)\n\n";
      }

      const char* SpecificationURL() override
      {
        return "https://www.fabbers.com/tech/STL_Format";
      }

      /// Return MIME type, NULL in this case.
      const char* GetMIMEType() override { return "application/sla"; }

      /// Return read/write flag: read only.
      unsigned int Flags() override
      {
        return WRITEONEONLY | NOTREADABLE;
      }

      /// Skip to object: used for multi-object file formats.
      int SkipObjects(int n, OpenBabel::OBConversion* pConv) override { return 0; }

      /// Read: always return false.
      bool ReadMolecule(OpenBabel::OBBase*, OpenBabel::OBConversion*) override
      {
        return false;
      }

      /// Write.
      bool WriteMolecule(OpenBabel::OBBase*, OpenBabel::OBConversion*) override;
  };


  // Global variable used to register STL format.
  STLFormat theSTLFormat;


  // Produces triangles over the surface of a sphere


  struct Triangle {
    vector3 a, b, c;
    uint16_t col;
    Triangle( vector3 x, vector3 y, vector3 z, uint16_t colour ) : a{x}, b{y}, c{z}, col{colour} {
    }
  };

  /**
   * @brief Helper struct containing the necessary information to visualize
   *        a bond between two atoms.
   */
  struct StlBondModelInfo {
    /**
     * @brief radius of the cylinder visualizing a single bond (Unit: Angstrom).
     *        Scales with the scale factor.
     */
    double radius;

    /**
     * @brief spacing between two bonds - only relevant for bond order > 1 (Unit: Angstrom).
     *        Scales with the scale factor.
     */
    double spacing;

    /**
     * @brief origin of the bond - e.g. the position of one of the two bonding atoms.
     */
    vector3 origin;

    /**
     * @brief direction and length of the bond in 3D space.
     */
    vector3 direction_and_length;

    /**
     * @brief colour for the first half of the bond - for example CPK colour of the start atom.
     */
    uint16_t first_half_colour;

    /**
     * @brief colour for the second half of the bond - for example CPK colour of the end atom.
     */
    uint16_t second_half_colour;

    /**
     * @brief order of the bond - see OBBond for more information.
     *        For visualization purposes this value is clamped between [1, 4].
     */
    unsigned int order;
  };

  static void map_sphere ( vector<Triangle> &triangles, vector3 origin, double r, uint16_t col )
  {
    static constexpr std::size_t LongitudeSteps{ 144 };
    static constexpr std::size_t LatitudeSteps{LongitudeSteps / 2};

    vector<vector3> points;

    double theta =  ( 2*M_PI / LongitudeSteps );
    int p2 = LongitudeSteps / 2;
    int r2 = LatitudeSteps / 2;
    for(int y = -r2; y < r2; ++y) {
      double cy = cos(y*theta);
      double cy1 = cos((y+1)*theta);
      double sy = sin(y*theta);
      double sy1 = sin((y+1)*theta);

      for(int i = -p2; i < p2; ++i) {
        double ci = cos(i*theta);
        double si = sin(i*theta);
        points.emplace_back(origin[0] + r * ci*cy , origin[1] + r*sy , origin[2] + r * si*cy);
        points.emplace_back(origin[0] + r * ci*cy1, origin[1] + r*sy1, origin[2] + r * si*cy1);
      }
    }

    for( int i=0 ; i < points.size()-2; i++ ) {
      // Order the points to obey the 'right-hand rule', so normals all face outwards
      if( i%2 == 0 ) {
        triangles.emplace_back(points[i], points[i+1], points[i+2], col);
      }
      else {
        triangles.emplace_back(points[i+2], points[i+1], points[i], col);
      }
    }
  }

  static void map_cylinder(vector<Triangle> &out, vector3 cylinder_origin,
                           vector3 cylinder_direction_and_length, double radius,
                           uint16_t color) {
    static constexpr uint32_t QuadsPerCylinder{128};
    static constexpr double AngleIncrease{2.0 * M_PI / QuadsPerCylinder};

    // plane for the two circles at the top and bottom of the cylinder
    vector3 plane_v0{};
    cylinder_direction_and_length.createOrthoVector(plane_v0);

    vector3 plane_v1 = cross(plane_v0, cylinder_direction_and_length);
    plane_v1.normalize();

    // remember 2 points as the starting point of the quad
    vector3 base_point = cylinder_origin + radius * plane_v1;
    vector3 top_point = base_point + cylinder_direction_and_length;
    double angle = 0.0;

    for (uint32_t i = 0; i < QuadsPerCylinder; ++i) {
      // compute the 2 next points along the top and bottom circle
      angle += AngleIncrease;
      vector3 next_base_point = cylinder_origin +
                                radius * sin(angle) * plane_v0 +
                                radius * cos(angle) * plane_v1;
      vector3 next_top_point = next_base_point + cylinder_direction_and_length;

      // create quad (2 triangles) from previous and next points
      // note the point order in order to obey the 'right-hand rule'
      out.emplace_back(next_base_point, top_point, base_point, color);
      out.emplace_back(top_point, next_top_point, next_base_point, color);

      // remember next points
      base_point = next_base_point;
      top_point = next_top_point;
    }
  }

  static void model_bond(vector<Triangle>& out, const StlBondModelInfo& bond)
  {
    const vector3 half_bond_vector = bond.direction_and_length * 0.5;
    const bool same_colour = (bond.first_half_colour == bond.second_half_colour);

    // compute how to shift bond visualizations in case of bond order > 1
    vector3 bond_plane_vector{};
    bond.direction_and_length.createOrthoVector(bond_plane_vector);
    const vector3 bond_visualization_offset = (2 * bond.radius + bond.spacing) * bond_plane_vector;
  
    // compute origin of the first bond visualization depending on the bond order
    const unsigned int order = max(unsigned(1), min(unsigned(4), bond.order));
    const double offset_from_spacing = (order - 1) * bond.spacing * 0.5;
    const double offset_from_bond_visualizations = (order - 1) * bond.radius;
    vector3 origin = bond.origin - (offset_from_bond_visualizations + offset_from_spacing) * bond_plane_vector;
  
    // draw as many bond visualizations as the bond order dictates
    for (uint32_t i = 0; i < order; ++i) {
      // colour the bond depending on the first/second half colour
      // if both colours are the same only make one cylinder in order to reduce the triangle count
      if (same_colour) {
        map_cylinder(out, origin, bond.direction_and_length, bond.radius,
                     bond.first_half_colour);
      } else {
        map_cylinder(out, origin, half_bond_vector, bond.radius, bond.first_half_colour);
        map_cylinder(out, origin + half_bond_vector, half_bond_vector, bond.radius,
                     bond.second_half_colour);
      }
      origin += bond_visualization_offset;
    }
  }

  /**
   * @brief Get STL compatible colour for the given atom based on the colours defined
   *        in OBElements::GetRGB().
   *        Note that colour in STL is a non-standard extension allocating 5 bits per
   *        colour channel with bit 15 set to 1.
   */
  static uint16_t stl_colour( unsigned int atomicnum , STLFormat::ColourOrder order ) {
    static constexpr uint8_t MaxColor{0x1F};
    static constexpr uint8_t ColorChannelWidth{5};
    static constexpr uint16_t Bit15SetTo1{0x800};

    array<double, 3> colour{};
    OBElements::GetRGB(atomicnum, &colour[0], &colour[1], &colour[2]);

    // adjust channel order
    if (order == STLFormat::ColourOrder::BlueGreenRed) {
      std::swap(colour[0], colour[2]);
    }

    // compress colour into 15 bits
    uint16_t stl_color{0};
    for (double channel : colour) {
      uint8_t scaled_color = static_cast<uint8_t>(std::round(channel * MaxColor));
      uint8_t clamped_color = std::min(scaled_color, MaxColor);

      stl_color <<= ColorChannelWidth;
      stl_color |= clamped_color;
    }

    return Bit15SetTo1 | stl_color;
  }

  static void write_stl_vector(ostream& out, const vector3& v) {
    array<float, 3> tmp{
      static_cast<float>(v.x()), static_cast<float>(v.y()), static_cast<float>(v.z())
    };
    out.write(reinterpret_cast<const char*>(tmp.data()), tmp.size() * sizeof(float));
  }

  static void write_stl_triangle(ostream& out, const Triangle& t) {
    write_stl_vector(out, t.a);
    write_stl_vector(out, t.b);
    write_stl_vector(out, t.c);
    out.write(reinterpret_cast<const char*>(&t.col), sizeof(t.col));
  }

  static void output_stl( ostream &os, const vector<Triangle> &triangles, bool cpk_colour ) {
    static constexpr size_t StlHeaderSizeInBytes{80};
    static constexpr array<uint8_t, 10> StlGlobalColorHeader{
      'C', 'O', 'L', 'O', 'R', '=', 0xFF, 0xFF, 0xFF, 0xFF
    };
    
    // STL header
    size_t zero_bytes_in_header = StlHeaderSizeInBytes;
    if (cpk_colour) {
      os.write(reinterpret_cast<const char*>(StlGlobalColorHeader.data()), StlGlobalColorHeader.size());
      zero_bytes_in_header -= StlGlobalColorHeader.size();
    }
    for (size_t i = 0; i < zero_bytes_in_header; ++i) {
      os.put(0);
    }

    uint32_t triangle_count = triangles.size();
    os.write( reinterpret_cast<const char*>(&triangle_count), sizeof(uint32_t) );

    // go through all triangles
    for(const auto& t : triangles) {
      // write zero normal - most programs ignore the normals
      write_stl_vector(os, VZero);

      // write triangle
      write_stl_triangle(os, t);
    }

    os.flush();
  }

  bool STLFormat::WriteMolecule( OBBase* pOb, OBConversion* pConv )
  {
    OBMol* pmol = dynamic_cast< OBMol* >(pOb);
    if (pmol == nullptr) return false;

    ostream& os = *pConv->GetOutStream();

    double probe_radius{DefaultProbeRadius};
    double scale_factor{DefaultScalingFactor};
    double bond_radius{DefaultBondRadius};
    double bond_spacing{DefaultBondSpacing};
    uint16_t col{NoColor};
    STLFormat::ColourOrder colour_order{DefaultColourOrder};

    bool cpk_colours{pConv->IsOption("c") != nullptr};
    bool draw_bonds{pConv->IsOption("b") != nullptr};

    if( pConv->IsOption( "p" ) ) {
      probe_radius = atof( pConv->IsOption("p", OBConversion::OUTOPTIONS) );
      if( !isfinite(probe_radius) || probe_radius < 0. ) { probe_radius = 0.; }
    }
    if( pConv->IsOption( "s" ) ) {
      scale_factor = atof( pConv->IsOption("s", OBConversion::OUTOPTIONS) );
      if( !isfinite(scale_factor) || scale_factor < 0. ) { scale_factor = 0.; }
    }
    if (pConv->IsOption("bond-radius", OBConversion::GENOPTIONS)) {
      bond_radius = atof(pConv->IsOption("bond-radius", OBConversion::GENOPTIONS));
      if (!isfinite(bond_radius) || bond_radius < 0.) { bond_radius = 0.; }
    }
    if (pConv->IsOption("bond-spacing", OBConversion::GENOPTIONS)) {
      bond_spacing = atof(pConv->IsOption("bond-spacing", OBConversion::GENOPTIONS));
      if (!isfinite(bond_spacing) || bond_spacing < 0.) { bond_spacing = 0.; }
    }
    if (pConv->IsOption("colour-order-rgb", OBConversion::GENOPTIONS)) {
      colour_order = STLFormat::ColourOrder::RedGreenBlue;
    }

    vector<Triangle> triangles;
    FOR_ATOMS_OF_MOL(a, *pmol) {
      const double vdwrad = scale_factor * OBElements::GetVdwRad( a->GetAtomicNum() ) + probe_radius;
      col = (cpk_colours) ? stl_colour(  a->GetAtomicNum(), colour_order ) : NoColor;
      map_sphere( triangles, a->GetVector(), vdwrad, col );
    }

    if (draw_bonds) {
      // pre-scale bond radius and spacing
      bond_radius *= scale_factor;
      bond_spacing *= scale_factor;

      FOR_BONDBFS_OF_MOL(b, *pmol) {
        const OBAtom& begin = *b->GetBeginAtom();
        const OBAtom& end = *b->GetEndAtom();

        StlBondModelInfo bond_info{};
        bond_info.radius = bond_radius;
        bond_info.spacing = bond_spacing;
        bond_info.origin = begin.GetVector();
        bond_info.direction_and_length = (end.GetVector() - begin.GetVector());
        bond_info.first_half_colour = (cpk_colours) ? stl_colour(begin.GetAtomicNum(), colour_order) : NoColor;
        bond_info.second_half_colour = (cpk_colours) ? stl_colour(end.GetAtomicNum(), colour_order) : NoColor;
        bond_info.order = b->GetBondOrder();

        model_bond(triangles, bond_info);
      }
    }

    output_stl ( os, triangles, cpk_colours );
    os.flush();
    return true;
  }

}

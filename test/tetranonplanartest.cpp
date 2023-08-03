#include "obtest.h"
#include <openbabel/stereo/tetrahedral.h>

using namespace std;
using namespace OpenBabel;

bool hasSameWinding(const OBStereo::Refs &refs1, const OBStereo::Refs &refs2)
{
  OB_REQUIRE( refs1.size() == 3 );
  OB_REQUIRE( refs2.size() == 3 );

  int Ni1 = OBStereo::NumInversions(refs1);
  int Ni2 = OBStereo::NumInversions(refs2);

  return ((Ni1 + Ni2) % 2 == 0);
}

//        2    3
//        -  ///
//        - //
//  0 --- 1
//         \
//          \
//           4


int tetranonplanartest(int argc, char* argv[])
{
  OBTetrahedralStereo::Config cfg;

  // set clockwise, viewing from 1
  cfg.from = 0;
  cfg.center = 1;
  cfg.refs = OBStereo::MakeRefs(2, 3, 4);
  OB_REQUIRE( cfg.from == 0 );
  OB_REQUIRE( cfg.center == 1 );
  OB_REQUIRE( cfg.refs.size() == 3 );
  OB_REQUIRE( cfg.refs[0] == 2 );
  OB_REQUIRE( cfg.refs[1] == 3 );
  OB_REQUIRE( cfg.refs[2] == 4 );

  // test nothing operation
  OBTetrahedralStereo::Config cfg2;
  cfg2 = OBTetraNonPlanarStereo::ToConfig(cfg, 0);
  OB_ASSERT( cfg == cfg2 );
  
  OBTetrahedralStereo::Config cfg3;
  // try viewing from other atom: 2
  cfg2 = OBTetraNonPlanarStereo::ToConfig(cfg, 2);
  OB_ASSERT( cfg2.center == 1 );
  OB_ASSERT( cfg2.from == 2 );
  OB_ASSERT( cfg2.refs.size() == 3 );
  OB_ASSERT( hasSameWinding(cfg2.refs, OBStereo::MakeRefs(3, 0, 4)) );
 
  // try viewing from other atom: 3
  cfg2 = OBTetraNonPlanarStereo::ToConfig(cfg, 3);
  OB_ASSERT( cfg2.center == 1 );
  OB_ASSERT( cfg2.from == 3 );
  OB_ASSERT( cfg2.refs.size() == 3 );
  OB_ASSERT( hasSameWinding(cfg2.refs, OBStereo::MakeRefs(0, 2, 4)) );
 
  // try viewing from other atom: 4
  cfg2 = OBTetraNonPlanarStereo::ToConfig(cfg, 4);
  OB_ASSERT( cfg2.center == 1 );
  OB_ASSERT( cfg2.from == 4 );
  OB_ASSERT( cfg2.refs.size() == 3 );
  OB_ASSERT( hasSameWinding(cfg2.refs, OBStereo::MakeRefs(3, 2, 0)) );

  // try viewing anti-clockwise 
  cfg2 = OBTetraNonPlanarStereo::ToConfig(cfg, 3, OBStereo::AntiClockwise);
  OB_ASSERT( cfg2.center == 1 );
  OB_ASSERT( cfg2.towards == 3 );
  OB_ASSERT( cfg2.refs.size() == 3 );
  OB_ASSERT( hasSameWinding(cfg2.refs, OBStereo::MakeRefs(2, 0, 4)) ); // CW <-> ACW = inversion

  // try viewing towards atom
  cfg2 = OBTetraNonPlanarStereo::ToConfig(cfg, 3, OBStereo::Clockwise, OBStereo::ViewTowards);
  OB_ASSERT( cfg2.center == 1 );
  OB_ASSERT( cfg2.towards == 3 );
  OB_ASSERT( cfg2.refs.size() == 3 );
  OB_ASSERT( hasSameWinding(cfg2.refs, OBStereo::MakeRefs(2, 0, 4)) ); // from <-> towards = inversion
 
  // try viewing towards atom anti-clockwise
  cfg2 = OBTetraNonPlanarStereo::ToConfig(cfg, 3, OBStereo::AntiClockwise, OBStereo::ViewTowards);
  OB_ASSERT( cfg2.center == 1 );
  OB_ASSERT( cfg2.towards == 3 );
  OB_ASSERT( cfg2.refs.size() == 3 );
  OB_ASSERT( hasSameWinding(cfg2.refs, OBStereo::MakeRefs(0, 2, 4)) ); // 2 permutations cancel out

  return 0;
}

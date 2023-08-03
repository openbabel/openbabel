#include "obtest.h"
#include <openbabel/stereo/cistrans.h>

using namespace std;
using namespace OpenBabel;

//        2    3
//        -  ///
//        - //
//  0 --- 1
//         \
//          \
//           4


int tetraplanartest(int argc, char* argv[])
{
  OBCisTransStereo::Config cfg;

  // set clockwise, viewing from 1
  cfg.begin = 0;
  cfg.end = 1;
  cfg.refs = OBStereo::MakeRefs(2, 3, 4, 5);
  OB_REQUIRE( cfg.begin == 0 );
  OB_REQUIRE( cfg.end == 1 );
  OB_REQUIRE( cfg.refs.size() == 4 );
  OB_REQUIRE( cfg.refs[0] == 2 );
  OB_REQUIRE( cfg.refs[1] == 3 );
  OB_REQUIRE( cfg.refs[2] == 4 );
  OB_REQUIRE( cfg.refs[3] == 5 );
  OB_REQUIRE( cfg.shape == OBStereo::ShapeU );

  // test nothing operation
  OBCisTransStereo::Config cfg2;
  cfg2 = OBTetraPlanarStereo::ToConfig(cfg, 2);
  OB_ASSERT( cfg == cfg2 );
  OB_ASSERT( cfg2.refs[0] == cfg.refs[0] );
  OB_ASSERT( cfg2.refs[1] == cfg.refs[1] );
  OB_ASSERT( cfg2.refs[2] == cfg.refs[2] );
  OB_ASSERT( cfg2.refs[3] == cfg.refs[3] );
  
  // try start = 3
  cfg2 = OBTetraPlanarStereo::ToConfig(cfg, 3);
  OB_ASSERT( cfg2.begin == 0 );
  OB_ASSERT( cfg2.end == 1 );
  OB_ASSERT( cfg2.refs.size() == 4 );
  OB_ASSERT( cfg2.refs[0] == 3 );
  OB_ASSERT( cfg2.refs[1] == 4 );
  OB_ASSERT( cfg2.refs[2] == 5 );
  OB_ASSERT( cfg2.refs[3] == 2 );
  OB_ASSERT( cfg2.shape == OBStereo::ShapeU );

  // try start = 5
  cfg2 = OBTetraPlanarStereo::ToConfig(cfg, 5);
  OB_ASSERT( cfg2.begin == 0 );
  OB_ASSERT( cfg2.end == 1 );
  OB_ASSERT( cfg2.refs.size() == 4 );
  OB_ASSERT( cfg2.refs[0] == 5 );
  OB_ASSERT( cfg2.refs[1] == 2 );
  OB_ASSERT( cfg2.refs[2] == 3 );
  OB_ASSERT( cfg2.refs[3] == 4 );
  OB_ASSERT( cfg2.shape == OBStereo::ShapeU );

  // try U -> Z
  OBCisTransStereo::Config shapeZ = OBTetraPlanarStereo::ToConfig(cfg, 2, OBStereo::ShapeZ);
  OB_ASSERT( shapeZ.begin == 0 );
  OB_ASSERT( shapeZ.end == 1 );
  OB_ASSERT( shapeZ.refs.size() == 4 );
  OB_ASSERT( shapeZ.refs[0] == 2 );
  OB_ASSERT( shapeZ.refs[1] == 3 );
  OB_ASSERT( shapeZ.refs[2] == 5 );
  OB_ASSERT( shapeZ.refs[3] == 4 );
  OB_ASSERT( shapeZ.shape == OBStereo::ShapeZ );

  // try U -> 4
  OBCisTransStereo::Config shape4 = OBTetraPlanarStereo::ToConfig(cfg, 2, OBStereo::Shape4);
  OB_ASSERT( shape4.begin == 0 );
  OB_ASSERT( shape4.end == 1 );
  OB_ASSERT( shape4.refs.size() == 4 );
  OB_ASSERT( shape4.refs[0] == 2 );
  OB_ASSERT( shape4.refs[1] == 4 );
  OB_ASSERT( shape4.refs[2] == 3 );
  OB_ASSERT( shape4.refs[3] == 5 );
  OB_ASSERT( shape4.shape == OBStereo::Shape4 );

  // try Z -> U
  OBCisTransStereo::Config shapeU = OBTetraPlanarStereo::ToConfig(shapeZ, 2, OBStereo::ShapeU);
  OB_ASSERT( shapeU.begin == 0 );
  OB_ASSERT( shapeU.end == 1 );
  OB_ASSERT( shapeU.refs.size() == 4 );
  OB_ASSERT( shapeU.refs[0] == 2 );
  OB_ASSERT( shapeU.refs[1] == 3 );
  OB_ASSERT( shapeU.refs[2] == 4 );
  OB_ASSERT( shapeU.refs[3] == 5 );
  OB_ASSERT( shapeU.shape == OBStereo::ShapeU );

  // try 4 -> U
  shapeU = OBTetraPlanarStereo::ToConfig(shape4, 2, OBStereo::ShapeU);
  OB_ASSERT( shapeU.begin == 0 );
  OB_ASSERT( shapeU.end == 1 );
  OB_ASSERT( shapeU.refs.size() == 4 );
  OB_ASSERT( shapeU.refs[0] == 2 );
  OB_ASSERT( shapeU.refs[1] == 3 );
  OB_ASSERT( shapeU.refs[2] == 4 );
  OB_ASSERT( shapeU.refs[3] == 5 );
  OB_ASSERT( shapeU.shape == OBStereo::ShapeU );

  return 0;
}

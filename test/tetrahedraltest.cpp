#include "obtest.h"
#include <openbabel/stereo/tetrahedral.h>

using namespace std;
using namespace OpenBabel;

bool sameWinding(const OBStereo::Refs &refs1, const OBStereo::Refs &refs2)
{
  OB_REQUIRE( refs1.size() == 3 );
  OB_REQUIRE( refs2.size() == 3 );

  int Ni1 = OBStereo::NumInversions(refs1);
  int Ni2 = OBStereo::NumInversions(refs2);

  return ((Ni1 + Ni2) % 2 == 0);
}
  
void test_configStruct()
{
  // reference Config
  OBTetrahedralStereo::Config reference(0, 1, OBStereo::MakeRefs(2, 3, 4), OBStereo::Clockwise, OBStereo::ViewFrom);

  // test copying
  OBTetrahedralStereo::Config referenceCopy = reference;
  OB_ASSERT( reference == referenceCopy );
  OB_ASSERT( referenceCopy == reference );

  // invalid center (chiral) id
  OBTetrahedralStereo::Config invalidCenter(45, 1, OBStereo::MakeRefs(2, 3, 4), OBStereo::Clockwise, OBStereo::ViewFrom);
  OB_ASSERT( reference != invalidCenter );

  // invalid from/towards id
  OBTetrahedralStereo::Config invalidFrom(0, 45, OBStereo::MakeRefs(2, 3, 4), OBStereo::Clockwise, OBStereo::ViewFrom);
  OB_ASSERT( reference != invalidFrom );

  // test other refs
  OBTetrahedralStereo::Config cfg1(0, 1, OBStereo::MakeRefs(2, 4, 3), OBStereo::Clockwise, OBStereo::ViewFrom);
  OB_ASSERT( reference != cfg1 );

  // test anti-clockwise
  OBTetrahedralStereo::Config cfg2(0, 1, OBStereo::MakeRefs(2, 3, 4), OBStereo::AntiClockwise, OBStereo::ViewFrom);
  OB_ASSERT( reference != cfg2 );

  // test viewing towards
  OBTetrahedralStereo::Config cfg3(0, 1, OBStereo::MakeRefs(2, 3, 4), OBStereo::Clockwise, OBStereo::ViewTowards);
  OB_ASSERT( reference != cfg3 );

  // test anti-clockwise + viewing towards
  OBTetrahedralStereo::Config cfg4(0, 1, OBStereo::MakeRefs(2, 3, 4), OBStereo::AntiClockwise, OBStereo::ViewTowards);
  OB_ASSERT( reference == cfg4 );

}

void test_IsValid()
{
  OBTetrahedralStereo ts(0);
  OBTetrahedralStereo::Config cfg, cfgCopy;
  cfg.center = 0;
  cfg.from = 1;
  cfg.refs = OBStereo::MakeRefs(2, 3, 4);

  ts.SetConfig(cfg);
  OB_ASSERT( ts.IsValid() );

  // no center
  cfgCopy = cfg;
  cfgCopy.center = OBStereo::NoRef;
  ts.SetConfig(cfgCopy);
  OB_ASSERT( !ts.IsValid() );

  // no from
  cfgCopy = cfg;
  cfgCopy.from = OBStereo::NoRef;
  ts.SetConfig(cfgCopy);
  OB_ASSERT( !ts.IsValid() );

  // no refs
  cfgCopy = cfg;
  cfgCopy.refs = std::vector<unsigned long>();
  ts.SetConfig(cfgCopy);
  OB_ASSERT( !ts.IsValid() );
}

void test_equalsOperator()
{
  OBTetrahedralStereo ts1(0), ts2(0);
  OBTetrahedralStereo::Config cfg;
  cfg.center = 0;
  cfg.from = 1;
  cfg.refs = OBStereo::MakeRefs(2, 3, 4);

  ts1.SetConfig(cfg);
  ts2.SetConfig(cfg);
  OB_ASSERT( ts1 == ts2 );

  cfg.winding = OBStereo::AntiClockwise;
  ts2.SetConfig(cfg);
  OB_ASSERT( ts1 != ts2 );
}

void test_GetSetConfig()
{
  OBTetrahedralStereo th(0);
  OBTetrahedralStereo::Config cfg;

  // set clockwise, viewing from 1
  OB_ASSERT( !th.IsValid() );
  cfg.refs = OBStereo::MakeRefs(2, 3, 4);
  cfg.from = 1;
  cfg.center = 5;
  th.SetConfig(cfg);
  OB_ASSERT( th.IsValid() );

  OBTetrahedralStereo::Config cfg2 = th.GetConfig();
  OB_ASSERT( cfg2.center == 5 );
  OB_ASSERT( cfg2.from == 1 );
  OB_ASSERT( cfg2.refs.size() == 3 );
  OB_ASSERT( cfg2.refs[0] == 2 );
  OB_ASSERT( cfg2.refs[1] == 3 );
  OB_ASSERT( cfg2.refs[2] == 4 );
  OB_ASSERT( cfg2.winding == OBStereo::Clockwise );
  OB_ASSERT( cfg2.view == OBStereo::ViewFrom );
  OB_ASSERT( cfg == cfg2 );

  // get viewing from 3
  cfg2 = th.GetConfig(3);
  OB_ASSERT( cfg2.from == 3 );
  OB_ASSERT( cfg == cfg2 );

  // get viewing from 3, AC
  cfg2 = th.GetConfig(3, OBStereo::AntiClockwise);
  OB_ASSERT( cfg2.winding == OBStereo::AntiClockwise );
  OB_ASSERT( cfg == cfg2 );

  // get viewing towards 3
  cfg2 = th.GetConfig(3, OBStereo::Clockwise, OBStereo::ViewTowards);
  OB_ASSERT( cfg2.view == OBStereo::ViewTowards );
  OB_ASSERT( cfg == cfg2 );

  // get viewing towards 3, AC
  cfg2 = th.GetConfig(3, OBStereo::AntiClockwise, OBStereo::ViewTowards);
  OB_ASSERT( cfg2.winding == OBStereo::AntiClockwise );
  OB_ASSERT( cfg2.view == OBStereo::ViewTowards );
  OB_ASSERT( cfg == cfg2 );
}

void test_Refs()
{
  OBTetrahedralStereo th(0);
 
  // center 2, view from 1, clockwise 9 4 34
  OBTetrahedralStereo::Config cfg;
  cfg.center = 2;
  cfg.from = 1;
  cfg.refs = OBStereo::MakeRefs(9, 4, 34);
  th.SetConfig(cfg);

  OB_REQUIRE( th.IsValid() );

  // 
  // test viewing from/towards all atoms
  //

  // from/towards 1
  OB_ASSERT( th.GetConfig() == cfg ); // from 1, clockwise
  cfg.winding = OBStereo::AntiClockwise;
  OB_ASSERT( th.GetConfig() != cfg ); // from 1, anti-clockwise
  cfg.view = OBStereo::ViewTowards;
  OB_ASSERT( th.GetConfig() == cfg ); // towards 1, anti-clockwise

  
  // from/towards 9
  OBTetrahedralStereo::Config cfg2;
  cfg2.center = 2;
  cfg2.from = 9;
  cfg2.refs = OBStereo::MakeRefs(4, 1, 34);
  OB_ASSERT( th.GetConfig() == cfg2 ); // from 9, clockwise
  cfg2.winding = OBStereo::AntiClockwise;
  OB_ASSERT( th.GetConfig() != cfg2 ); // from 9, anti-clockwise
  cfg2.view = OBStereo::ViewTowards;
  OB_ASSERT( th.GetConfig() == cfg2 ); // towards 9, anti-clockwise

  // from/towards 4
  OBTetrahedralStereo::Config cfg3;
  cfg3.center = 2;
  cfg3.from = 4;
  cfg3.refs = OBStereo::MakeRefs(1, 9, 34);
  OB_ASSERT( th.GetConfig() == cfg3 ); // from 4, clockwise
  cfg3.winding = OBStereo::AntiClockwise;
  OB_ASSERT( th.GetConfig() != cfg3 ); // from 4, anti-clockwise
  cfg3.view = OBStereo::ViewTowards;
  OB_ASSERT( th.GetConfig() == cfg3 ); // towards 4, anti-clockwise

  // from/towards 34
  OBTetrahedralStereo::Config cfg4;
  cfg4.center = 2;
  cfg4.from = 34;
  cfg4.refs = OBStereo::MakeRefs(4, 9, 1);
  OB_ASSERT( th.GetConfig() == cfg4 ); // from 34, clockwise
  cfg4.winding = OBStereo::AntiClockwise;
  OB_ASSERT( th.GetConfig() != cfg4 ); // from 34, anti-clockwise
  cfg4.view = OBStereo::ViewTowards;
  OB_ASSERT( th.GetConfig() == cfg4 ); // towards 34, anti-clockwise

  //
  // test invalid Config structs
  //
  OBTetrahedralStereo::Config cfg5;
  OB_ASSERT( th.GetConfig() != cfg5 );
  cfg5.center = 2;
  OB_ASSERT( th.GetConfig() != cfg5 );
  cfg5.from = 34;
  OB_ASSERT( th.GetConfig() != cfg5 ); 
  cfg5.refs = OBStereo::MakeRefs(4, 9, 1);
  OB_ASSERT( th.GetConfig() == cfg5 ); 
  cfg5.from = OBStereo::NoRef;
  OB_ASSERT( th.GetConfig() != cfg5 ); 
  cfg5.center = OBStereo::NoRef;
  cfg5.from = 34;
  OB_ASSERT( th.GetConfig() != cfg5 ); 
}

int main()
{
  test_configStruct();
  test_IsValid();
  test_equalsOperator();
  test_GetSetConfig();
  test_Refs();

  return 0;
}

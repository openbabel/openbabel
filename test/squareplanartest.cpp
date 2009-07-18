#include <boost/test/minimal.hpp>
#include <cassert>

#include <iostream>

#include <openbabel/mol.h>
#include "squareplanar.h"

using namespace std;
using namespace OpenBabel;


void testGetType ()
{
  OBMol mol;
  OBSquarePlanarStereo sp(&mol);
  BOOST_REQUIRE( sp.GetType() == OBStereo::SquarePlanar );
}

void testCenter()
{
  OBMol mol;
  OBSquarePlanarStereo sp(&mol);
  sp.SetCenter(42); 
  BOOST_REQUIRE( sp.GetCenter() );
}

void testIsValid()
{
  OBMol mol;
  OBSquarePlanarStereo sp(&mol);
  BOOST_REQUIRE( !sp.IsValid() );
  // set center atom
  sp.SetCenter(5);
  BOOST_REQUIRE( !sp.IsValid() );
  // set reference atoms
  sp.SetRefs( OBStereo::MakeRefs(1, 2, 3, 4) );
  // the object should now be valid
  BOOST_REQUIRE( sp.IsValid() );
}

// test basic ref setting/getting
void testRefs1()
{
  OBMol mol;
  OBSquarePlanarStereo sp(&mol);

  // 1   4     1   4     1---4     1   4
  //  \ /      |   |        /       \ /
  //   C    =  | C |  =    C    =    C
  //  / \      |   |      /         / \
  // 2   3     2---3     2---3     2---3

  // set refs using default U shape
  sp.SetRefs( OBStereo::MakeRefs(1, 2, 3, 4) );

  vector<unsigned long> refs;
  // get refs using default U shape
  refs = sp.GetRefs();
  BOOST_REQUIRE( refs.size() == 4 );
  BOOST_REQUIRE( refs[0] == 1 );
  BOOST_REQUIRE( refs[1] == 2 );
  BOOST_REQUIRE( refs[2] == 3 );
  BOOST_REQUIRE( refs[3] == 4 );

  // get refs using Z shape
  refs = sp.GetRefs(OBStereo::ShapeZ);
  BOOST_REQUIRE( refs.size() == 4 );
  BOOST_REQUIRE( refs[0] == 1 );
  BOOST_REQUIRE( refs[1] == 4 );
  BOOST_REQUIRE( refs[2] == 2 );
  BOOST_REQUIRE( refs[3] == 3 );

  // get refs using 4 shape
  refs = sp.GetRefs(OBStereo::Shape4);
  BOOST_REQUIRE( refs.size() == 4 );
  BOOST_REQUIRE( refs[0] == 1 );
  BOOST_REQUIRE( refs[1] == 3 );
  BOOST_REQUIRE( refs[2] == 2 );
  BOOST_REQUIRE( refs[3] == 4 );

}
 
void testRefs2()
{
  OBMol mol;
  OBSquarePlanarStereo sp(&mol);

  // 1   4     1   4     1---4     1   4
  //  \ /      |   |        /       \ /
  //   C    =  | C |  =    C    =    C
  //  / \      |   |      /         / \
  // 2   3     2---3     2---3     2---3

  // set refs using default U shape
  sp.SetRefs( OBStereo::MakeRefs(1, 2, 3, 4) );

  vector<unsigned long> refs;
  // get refs using default U shape starting from 1
  refs = sp.GetRefs((unsigned long)1);
  BOOST_REQUIRE( refs.size() == 4 );
  BOOST_REQUIRE( refs[0] == 1 );
  BOOST_REQUIRE( refs[1] == 2 );
  BOOST_REQUIRE( refs[2] == 3 );
  BOOST_REQUIRE( refs[3] == 4 );

  // get refs using Z shape
  refs = sp.GetRefs(OBStereo::ShapeZ);
  BOOST_REQUIRE( refs.size() == 4 );
  BOOST_REQUIRE( refs[0] == 1 );
  BOOST_REQUIRE( refs[1] == 4 );
  BOOST_REQUIRE( refs[2] == 2 );
  BOOST_REQUIRE( refs[3] == 3 );

  // get refs using 4 shape
  refs = sp.GetRefs(OBStereo::Shape4);
  BOOST_REQUIRE( refs.size() == 4 );
  BOOST_REQUIRE( refs[0] == 1 );
  BOOST_REQUIRE( refs[1] == 3 );
  BOOST_REQUIRE( refs[2] == 2 );
  BOOST_REQUIRE( refs[3] == 4 );

}
  
void testCisTrans()
{
  OBMol mol;
  OBSquarePlanarStereo sp(&mol);
  // set refs using default U shape
  sp.SetRefs( OBStereo::MakeRefs(1, 2, 3, 4) );

  //
  // Trans
  //

  // test invalid ids
  BOOST_REQUIRE( !sp.IsTrans(43, 3) );
  BOOST_REQUIRE( !sp.IsTrans(1, 1) );
  
  // test real trans refs in all combinations
  BOOST_REQUIRE( sp.IsTrans(1, 3) );
  BOOST_REQUIRE( sp.IsTrans(3, 1) );
  BOOST_REQUIRE( sp.IsTrans(2, 4) );
  BOOST_REQUIRE( sp.IsTrans(4, 2) );

  // test cis atoms, should not be trans...
  BOOST_REQUIRE( !sp.IsTrans(1, 2) );
  BOOST_REQUIRE( !sp.IsTrans(1, 4) );
  BOOST_REQUIRE( !sp.IsTrans(2, 1) );
  BOOST_REQUIRE( !sp.IsTrans(2, 3) );
  BOOST_REQUIRE( !sp.IsTrans(3, 2) );
  BOOST_REQUIRE( !sp.IsTrans(3, 4) );
  BOOST_REQUIRE( !sp.IsTrans(4, 1) );
  BOOST_REQUIRE( !sp.IsTrans(4, 3) );

  // 
  // Cis
  //

  // test invalid ids
  BOOST_REQUIRE( !sp.IsCis(43, 3) );
  BOOST_REQUIRE( !sp.IsCis(1, 1) );
  
  // test real cis refs in all combinations
  BOOST_REQUIRE( sp.IsCis(1, 2) );
  BOOST_REQUIRE( sp.IsCis(2, 1) );
  BOOST_REQUIRE( sp.IsCis(1, 4) );
  BOOST_REQUIRE( sp.IsCis(4, 1) );

  BOOST_REQUIRE( sp.IsCis(2, 1) );
  BOOST_REQUIRE( sp.IsCis(2, 3) );
  BOOST_REQUIRE( sp.IsCis(3, 2) );
  BOOST_REQUIRE( sp.IsCis(3, 4) );
  BOOST_REQUIRE( sp.IsCis(4, 3) );
  BOOST_REQUIRE( sp.IsCis(4, 1) );

  // test trans atoms, should not be cis...
  BOOST_REQUIRE( !sp.IsCis(1, 3) );
  BOOST_REQUIRE( !sp.IsCis(2, 4) );

  // test GetTransRef 
  BOOST_REQUIRE( sp.GetTransRef(1) == 3);
  BOOST_REQUIRE( sp.GetTransRef(3) == 1);
  BOOST_REQUIRE( sp.GetTransRef(2) == 4);
  BOOST_REQUIRE( sp.GetTransRef(4) == 2);
   // test GetCisRef 
  vector<unsigned long> cis = sp.GetCisRefs(1);
  BOOST_REQUIRE( cis.size() == 2 );
  BOOST_REQUIRE( cis[0] == 2 || cis[1] == 2);
  BOOST_REQUIRE( cis[0] == 4 || cis[1] == 4);
  cis = sp.GetCisRefs(4);
  BOOST_REQUIRE( cis.size() == 2 );
  BOOST_REQUIRE( cis[0] == 1 || cis[1] == 1);
  BOOST_REQUIRE( cis[0] == 3 || cis[1] == 3);
  cis = sp.GetCisRefs(2);
  BOOST_REQUIRE( cis.size() == 2 );
  BOOST_REQUIRE( cis[0] == 1 || cis[1] == 1);
  BOOST_REQUIRE( cis[0] == 3 || cis[1] == 3);
 
}

void testCompare()
{
  OBMol mol;
  OBSquarePlanarStereo sp(&mol);
  // set refs using default U shape
  sp.SetRefs( OBStereo::MakeRefs(1, 2, 3, 4) );
  sp.SetCenter(5);
  
  BOOST_REQUIRE( sp.Compare(OBStereo::MakeRefs(1, 2, 3, 4), OBStereo::ShapeU) );
  
  BOOST_REQUIRE( sp.Compare(OBStereo::MakeRefs(1, 2, 4, 3), OBStereo::ShapeZ) );
  BOOST_REQUIRE( !sp.Compare(OBStereo::MakeRefs(1, 2, 3, 4), OBStereo::ShapeZ) );
  
  BOOST_REQUIRE( sp.Compare(OBStereo::MakeRefs(1, 3, 2, 4), OBStereo::Shape4) );
  BOOST_REQUIRE( !sp.Compare(OBStereo::MakeRefs(1, 2, 3, 4), OBStereo::Shape4) );
 
}
  

int test_main(int,char**) 
{
  testGetType();
  testCenter();
  testIsValid();
  testRefs1();
  testRefs2();
  testCisTrans();
  testCompare();
  
  return 0;
}

                

#include "obtest.h"

#include <iostream>

#include <openbabel/mol.h>
#include <openbabel/stereo/squareplanar.h>

using namespace std;
using namespace OpenBabel;

void testGetType ()
{
  OBMol mol;
  OBSquarePlanarStereo sp(&mol);
  OB_ASSERT( sp.GetType() == OBStereo::SquarePlanar );
}

void testCenter()
{
  OBMol mol;
  OBSquarePlanarStereo sp(&mol);
  sp.SetCenter(42); 
  OB_ASSERT( sp.GetCenter() );
}

void testIsValid()
{
  OBMol mol;
  OBSquarePlanarStereo sp(&mol);
  OB_ASSERT( !sp.IsValid() );
  // set center atom
  sp.SetCenter(5);
  OB_ASSERT( !sp.IsValid() );
  // set reference atoms
  sp.SetRefs( OBStereo::MakeRefs(1, 2, 3, 4) );
  // the object should now be valid
  OB_ASSERT( sp.IsValid() );
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
  OB_ASSERT( refs.size() == 4 );
  OB_ASSERT( refs[0] == 1 );
  OB_ASSERT( refs[1] == 2 );
  OB_ASSERT( refs[2] == 3 );
  OB_ASSERT( refs[3] == 4 );

  // get refs using Z shape
  refs = sp.GetRefs(OBStereo::ShapeZ);
  OB_ASSERT( refs.size() == 4 );
  OB_ASSERT( refs[0] == 1 );
  OB_ASSERT( refs[1] == 4 );
  OB_ASSERT( refs[2] == 2 );
  OB_ASSERT( refs[3] == 3 );

  // get refs using 4 shape
  refs = sp.GetRefs(OBStereo::Shape4);
  OB_ASSERT( refs.size() == 4 );
  OB_ASSERT( refs[0] == 1 );
  OB_ASSERT( refs[1] == 3 );
  OB_ASSERT( refs[2] == 2 );
  OB_ASSERT( refs[3] == 4 );

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
  OB_ASSERT( refs.size() == 4 );
  OB_ASSERT( refs[0] == 1 );
  OB_ASSERT( refs[1] == 2 );
  OB_ASSERT( refs[2] == 3 );
  OB_ASSERT( refs[3] == 4 );

  // get refs using Z shape
  refs = sp.GetRefs(OBStereo::ShapeZ);
  OB_ASSERT( refs.size() == 4 );
  OB_ASSERT( refs[0] == 1 );
  OB_ASSERT( refs[1] == 4 );
  OB_ASSERT( refs[2] == 2 );
  OB_ASSERT( refs[3] == 3 );

  // get refs using 4 shape
  refs = sp.GetRefs(OBStereo::Shape4);
  OB_ASSERT( refs.size() == 4 );
  OB_ASSERT( refs[0] == 1 );
  OB_ASSERT( refs[1] == 3 );
  OB_ASSERT( refs[2] == 2 );
  OB_ASSERT( refs[3] == 4 );

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
  OB_ASSERT( !sp.IsTrans(43, 3) );
  OB_ASSERT( !sp.IsTrans(1, 1) );
  
  // test real trans refs in all combinations
  OB_ASSERT( sp.IsTrans(1, 3) );
  OB_ASSERT( sp.IsTrans(3, 1) );
  OB_ASSERT( sp.IsTrans(2, 4) );
  OB_ASSERT( sp.IsTrans(4, 2) );

  // test cis atoms, should not be trans...
  OB_ASSERT( !sp.IsTrans(1, 2) );
  OB_ASSERT( !sp.IsTrans(1, 4) );
  OB_ASSERT( !sp.IsTrans(2, 1) );
  OB_ASSERT( !sp.IsTrans(2, 3) );
  OB_ASSERT( !sp.IsTrans(3, 2) );
  OB_ASSERT( !sp.IsTrans(3, 4) );
  OB_ASSERT( !sp.IsTrans(4, 1) );
  OB_ASSERT( !sp.IsTrans(4, 3) );

  // 
  // Cis
  //

  // test invalid ids
  OB_ASSERT( !sp.IsCis(43, 3) );
  OB_ASSERT( !sp.IsCis(1, 1) );
  
  // test real cis refs in all combinations
  OB_ASSERT( sp.IsCis(1, 2) );
  OB_ASSERT( sp.IsCis(2, 1) );
  OB_ASSERT( sp.IsCis(1, 4) );
  OB_ASSERT( sp.IsCis(4, 1) );

  OB_ASSERT( sp.IsCis(2, 1) );
  OB_ASSERT( sp.IsCis(2, 3) );
  OB_ASSERT( sp.IsCis(3, 2) );
  OB_ASSERT( sp.IsCis(3, 4) );
  OB_ASSERT( sp.IsCis(4, 3) );
  OB_ASSERT( sp.IsCis(4, 1) );

  // test trans atoms, should not be cis...
  OB_ASSERT( !sp.IsCis(1, 3) );
  OB_ASSERT( !sp.IsCis(2, 4) );

  // test GetTransRef 
  OB_ASSERT( sp.GetTransRef(1) == 3);
  OB_ASSERT( sp.GetTransRef(3) == 1);
  OB_ASSERT( sp.GetTransRef(2) == 4);
  OB_ASSERT( sp.GetTransRef(4) == 2);
   // test GetCisRef 
  vector<unsigned long> cis = sp.GetCisRefs(1);
  OB_ASSERT( cis.size() == 2 );
  OB_ASSERT( cis[0] == 2 || cis[1] == 2);
  OB_ASSERT( cis[0] == 4 || cis[1] == 4);
  cis = sp.GetCisRefs(4);
  OB_ASSERT( cis.size() == 2 );
  OB_ASSERT( cis[0] == 1 || cis[1] == 1);
  OB_ASSERT( cis[0] == 3 || cis[1] == 3);
  cis = sp.GetCisRefs(2);
  OB_ASSERT( cis.size() == 2 );
  OB_ASSERT( cis[0] == 1 || cis[1] == 1);
  OB_ASSERT( cis[0] == 3 || cis[1] == 3);
 
}

void testCompare()
{
  OBMol mol;
  OBSquarePlanarStereo sp(&mol);
  // set refs using default U shape
  sp.SetRefs( OBStereo::MakeRefs(1, 2, 3, 4) );
  sp.SetCenter(5);
  
  OB_ASSERT( sp.Compare(OBStereo::MakeRefs(1, 2, 3, 4), OBStereo::ShapeU) );
  
  OB_ASSERT( sp.Compare(OBStereo::MakeRefs(1, 2, 4, 3), OBStereo::ShapeZ) );
  OB_ASSERT( !sp.Compare(OBStereo::MakeRefs(1, 2, 3, 4), OBStereo::ShapeZ) );
  
  OB_ASSERT( sp.Compare(OBStereo::MakeRefs(1, 3, 2, 4), OBStereo::Shape4) );
  OB_ASSERT( !sp.Compare(OBStereo::MakeRefs(1, 2, 3, 4), OBStereo::Shape4) );
 
}
  

int main() 
{
  // Define location of file formats for testing
  #ifdef FORMATDIR
    char env[BUFF_SIZE];
    snprintf(env, BUFF_SIZE, "BABEL_LIBDIR=%s", FORMATDIR);
    putenv(env);
  #endif

  testGetType();
  testCenter();
  testIsValid();
  testRefs1();
  testRefs2();
  testCisTrans();
  testCompare();
  
  cout << "end" << endl;

  return 0;
}

                

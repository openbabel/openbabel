#include "obtest.h"
#include <openbabel/stereo/stereo.h>

using namespace std;
using namespace OpenBabel;

void test_MakeRefs()
{
  OBStereo::Refs refs = OBStereo::MakeRefs(1, 3, 5);
  OB_ASSERT( refs.size() == 3 );
  OB_ASSERT( refs[0] == 1 );
  OB_ASSERT( refs[1] == 3 );
  OB_ASSERT( refs[2] == 5 );

  refs = OBStereo::MakeRefs(1, 3, 5, 7);
  OB_ASSERT( refs.size() == 4 );
  OB_ASSERT( refs[0] == 1 );
  OB_ASSERT( refs[1] == 3 );
  OB_ASSERT( refs[2] == 5 );
  OB_ASSERT( refs[3] == 7 );
}

void test_ContainsSameRefs()
{
  OB_ASSERT( OBStereo::ContainsSameRefs(OBStereo::MakeRefs(1, 2, 3), 
                                        OBStereo::MakeRefs(2, 3, 1)) );

  OB_ASSERT( !OBStereo::ContainsSameRefs(OBStereo::MakeRefs(1, 2, 3), 
                                         OBStereo::MakeRefs(2, 7, 1)) );

  // different sizes
  OB_ASSERT( !OBStereo::ContainsSameRefs(OBStereo::MakeRefs(1, 2, 3), 
                                        OBStereo::MakeRefs(2, 3, 1, 4)) );
}

void test_NumInversions()
{
  OB_ASSERT( OBStereo::NumInversions(OBStereo::MakeRefs(1, 2, 3)) == 0 );
  OB_ASSERT( OBStereo::NumInversions(OBStereo::MakeRefs(1, 3, 2)) == 1 );
  OB_ASSERT( OBStereo::NumInversions(OBStereo::MakeRefs(2, 1, 3)) == 1 );
  OB_ASSERT( OBStereo::NumInversions(OBStereo::MakeRefs(2, 3, 1)) == 2 );
  OB_ASSERT( OBStereo::NumInversions(OBStereo::MakeRefs(3, 1, 2)) == 2 );
  OB_ASSERT( OBStereo::NumInversions(OBStereo::MakeRefs(3, 2, 1)) == 3 );
}

void test_Permutate()
{
  OBStereo::Refs refs = OBStereo::MakeRefs(1, 2, 3);
  OBStereo::Permutate(refs, 1, 2);

  OB_ASSERT( refs[0] == 1 );
  OB_ASSERT( refs[1] == 3 );
  OB_ASSERT( refs[2] == 2 );

  // test out of range indexes for crashes
  OBStereo::Permutate(refs, -1, 1);
  OBStereo::Permutate(refs, 1, -1);
  OBStereo::Permutate(refs, 9, 1);
  OBStereo::Permutate(refs, 1, 9);
}

void test_Permutated()
{
  OBStereo::Refs refs = OBStereo::MakeRefs(1, 2, 3);
  OBStereo::Refs mutant = OBStereo::Permutated(refs, 1, 2);

  OB_ASSERT( refs[0] == 1 );
  OB_ASSERT( refs[1] == 2 );
  OB_ASSERT( refs[2] == 3 );

  OB_ASSERT( mutant[0] == 1 );
  OB_ASSERT( mutant[1] == 3 );
  OB_ASSERT( mutant[2] == 2 );

}

int main()
{
  test_MakeRefs();
  test_ContainsSameRefs();
  test_NumInversions();
  test_Permutate();
  test_Permutated();

  return 0;
}

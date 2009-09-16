#include "obtest.h"
#include <openbabel/stereo/stereo.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/graphsym.h>

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


std::string GetFilename(const std::string &filename)
{
  string path = TESTDATADIR + filename;
  return path;
}

bool doStereoPerception(OBMol &mol, int numTetrahedral, int numCisTrans)
{
  // need to calculate symmetry first
  std::vector<unsigned int> symmetry_classes;
  OpenBabel::OBGraphSym graphsym(&mol);
  graphsym.GetSymmetry(symmetry_classes);

  std::vector<OpenBabel::StereogenicUnit> units = FindStereogenicUnits(&mol, symmetry_classes);

  int tetrahedralCount = 0;
  int cistransCount = 0;
  for (unsigned int i = 0; i < units.size(); ++i) {
    switch (units.at(i).type) {
      case OBStereo::Tetrahedral:
        tetrahedralCount++;
        break;
      case OBStereo::CisTrans:
        cistransCount++;
        break;
    }
  }
  
  OB_COMPARE(tetrahedralCount, numTetrahedral);
  OB_COMPARE(cistransCount, numCisTrans);

  return (tetrahedralCount == numTetrahedral) && (cistransCount == numCisTrans);
}


/**
 * Test detection of stereocenters
 *
 *
 */
void test_StereoPerception()
{
  OBMol mol;
  OBConversion conv;
  OB_ASSERT( conv.SetInFormat("mol") );


  cout << "structure perception1" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/perception1.mol")) );
  OB_ASSERT( doStereoPerception(mol, 0, 0) );

  cout << "structure perception2" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/perception2.mol")) );
  OB_ASSERT( doStereoPerception(mol, 0, 0) );

  cout << "structure perception3" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/perception3.mol")) );
  OB_ASSERT( doStereoPerception(mol, 2, 0) );

  cout << "structure perception4" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/perception4.mol")) );
  OB_ASSERT( doStereoPerception(mol, 0, 0) );

  cout << "structure perception5" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/perception5.mol")) );
  OB_ASSERT( doStereoPerception(mol, 0, 0) );

  cout << "structure perception6" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/perception6.mol")) );
  OB_ASSERT( doStereoPerception(mol, 0, 2) );

  cout << "structure perception7" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/perception7.mol")) );
  OB_ASSERT( doStereoPerception(mol, 0, 0) );

  cout << "structure perception8" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/perception8.mol")) );
  OB_ASSERT( doStereoPerception(mol, 0, 0) );

  cout << "structure perception9" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/perception9.mol")) );
  OB_ASSERT( doStereoPerception(mol, 0, 0) );

  cout << "structure perception10" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/perception10.mol")) );
  OB_ASSERT( doStereoPerception(mol, 1, 1) );

  cout << "structure perception11" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/perception11.mol")) );
  OB_ASSERT( doStereoPerception(mol, 3, 0) );

  cout << "structure perception12" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/perception12.mol")) );
  OB_ASSERT( doStereoPerception(mol, 0, 0) );

  cout << "structure perception13" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/perception13.mol")) );
  OB_ASSERT( doStereoPerception(mol, 0, 0) );

  cout << "structure perception14" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/perception14.mol")) );
  OB_ASSERT( doStereoPerception(mol, 0, 0) );

  cout << "structure perception15" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/perception15.mol")) );
  OB_ASSERT( doStereoPerception(mol, 0, 0) );

  cout << "structure perception16" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/perception16.mol")) );
  OB_ASSERT( doStereoPerception(mol, 1, 2) );


  /*
   * J. Chem. Inf. Comput. Sci., Vol. 33, No. 6, 1993
   *
   * Figure 1. Examples of compounds containing parastereocenters
   */
  cout << "Razinger paper, fig. 1: structure a" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig1_a.mol")) );
  OB_ASSERT( doStereoPerception(mol, 0, 2) );

  cout << "Razinger paper, fig. 1: structure b" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig1_b.mol")) ); 
  OB_ASSERT( doStereoPerception(mol, 2, 2) );

  cout << "Razinger paper, fig. 1: structure c" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig1_c.mol")) ); 
  OB_ASSERT( doStereoPerception(mol, 2, 0) );

  cout << "Razinger paper, fig. 1: structure d" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig1_d.mol")) ); 
  OB_ASSERT( doStereoPerception(mol, 3, 0) );

  cout << "Razinger paper, fig. 1: structure e" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig1_e.mol")) ); 
  OB_ASSERT( doStereoPerception(mol, 6, 0) );

  cout << "Razinger paper, fig. 1: structure f" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig1_f.mol")) ); 
  OB_ASSERT( doStereoPerception(mol, 2, 1) );

  cout << "Razinger paper, fig. 1: structure g" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig1_g.mol")) ); 
  OB_ASSERT( doStereoPerception(mol, 0, 0) );

  cout << "Razinger paper, fig. 1: structure h" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig1_h.mol")) ); 
  OB_ASSERT( doStereoPerception(mol, 2, 0) );

  cout << "Razinger paper, fig. 1: structure i" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig1_i.mol")) ); 
  OB_ASSERT( doStereoPerception(mol, 3, 0) );



//  conv.ReadString(&mol, "CC1CCC(C)CC1");
}




int main()
{
  test_MakeRefs();
  test_ContainsSameRefs();
  test_NumInversions();
  test_Permutate();
  test_Permutated();
  test_StereoPerception();

  return 0;
}

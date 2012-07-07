#include "obtest.h"
#include <openbabel/stereo/tetrahedral.h>
#include <openbabel/stereo/cistrans.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>

using namespace std;
using namespace OpenBabel;

std::string GetFilename(const std::string &filename)
{
  string path = TESTDATADIR + filename;

  return path;
}

std::string test_singleTetrahedral(const std::string &file, 
    const OBTetrahedralStereo::Config &correct)
{
  OBMol mol;
  OBConversion conv;
  conv.SetInFormat("mol");
  conv.SetOutFormat("can");

  const std::string filename = GetFilename(file);
  ifstream ifs(filename.c_str());
  if (!ifs) {
    cerr << "Could not open " << filename << endl;
    return std::string();
  }

  conv.Read(&mol, &ifs);
  ifs.close();

  std::vector<OBGenericData *> stereoData = mol.GetAllData(OBGenericDataType::StereoData);
  OB_ASSERT( stereoData.size() == 1 );

  // compare the stereochemistry
  for (std::vector<OBGenericData*>::iterator data = stereoData.begin(); data != stereoData.end(); ++data) {
    if (((OBStereoBase*)*data)->GetType() == OBStereo::Tetrahedral) {
      OBTetrahedralStereo *ts = dynamic_cast<OBTetrahedralStereo*>(*data);
      OB_ASSERT( ts->GetConfig() == correct );
      if ( ts->GetConfig() != correct ) {
        cout << "found = " << ts->GetConfig() << endl;
        cout << "correct = " << correct << endl;
      }

      OBTetrahedralStereo::Config cfg = ts->GetConfig();
      OB_ASSERT( cfg.specified );
      // change refs 
      OBStereo::Permutate(cfg.refs, 1, 2);
      OB_ASSERT( cfg != correct );
      // change winding
      cfg.winding = (cfg.winding == OBStereo::Clockwise) ? OBStereo::AntiClockwise : OBStereo::Clockwise;
      OB_ASSERT( cfg == correct );
      // change view
      cfg.view = (cfg.view == OBStereo::ViewFrom) ? OBStereo::ViewTowards : OBStereo::ViewFrom;
      OB_ASSERT( cfg != correct );
      cfg.view = (cfg.view == OBStereo::ViewFrom) ? OBStereo::ViewTowards : OBStereo::ViewFrom;
      OB_ASSERT( cfg == correct );
      // change center
      cfg.center = 3994;
      OB_ASSERT( cfg != correct );


    }
  }

  return conv.WriteString(&mol);
}

std::string test_singleCisTrans(const std::string &file, 
    const OBCisTransStereo::Config &correct)
{
  OBMol mol;
  OBConversion conv;
  conv.SetInFormat("mol");
  conv.SetOutFormat("can");

  const std::string filename = GetFilename(file);
  ifstream ifs(filename.c_str());
  if (!ifs) {
    cerr << "Could not open " << filename << endl;
    return std::string();
  }

  conv.Read(&mol, &ifs);
  ifs.close();

  std::vector<OBGenericData *> stereoData = mol.GetAllData(OBGenericDataType::StereoData);
  OB_ASSERT( stereoData.size() == 1 );

  // compare the stereochemistry
  for (std::vector<OBGenericData*>::iterator data = stereoData.begin(); data != stereoData.end(); ++data) {
    if (((OBStereoBase*)*data)->GetType() == OBStereo::CisTrans) {
      OBCisTransStereo *ct = dynamic_cast<OBCisTransStereo*>(*data);
      OB_ASSERT( ct->GetConfig() == correct );
      if ( ct->GetConfig() != correct ) {
        cout << "found = " << ct->GetConfig() << endl;
        cout << "correct = " << correct << endl;
      }
    }
  }

  return conv.WriteString(&mol);
}

std::string readMol(OBMol *pmol, const std::string &file)
{
  OBConversion conv;
  conv.SetInFormat("mol");
  conv.SetOutFormat("can");

  const std::string filename = GetFilename(file);
  ifstream ifs(filename.c_str());
  if (!ifs) {
    cerr << "Could not open " << filename << endl;
    return NULL;
  }

  conv.Read(pmol, &ifs);
  ifs.close();
  return conv.WriteString(pmol);
}

std::string test_singleUnspecifiedTetrahedral(const std::string &file, 
    unsigned long center)
{
  OBMol mol;
  string retval = readMol(&mol, file);
  if (retval.size() > 0) {
    OBStereoFacade stereo(&mol);
    OB_ASSERT( stereo.HasTetrahedralStereo(center) );
    OBTetrahedralStereo *ts = stereo.GetTetrahedralStereo(center);
    OB_ASSERT( ts->GetConfig().specified == false );
  }

  return retval;
}

std::string test_singleUnknownTetrahedral(const std::string &file, 
    unsigned long center)
{
  OBMol mol;
  string retval = readMol(&mol, file);
  if (retval.size() > 0) {
    OBStereoFacade stereo(&mol);
    OB_ASSERT( stereo.HasTetrahedralStereo(center) );
    OBTetrahedralStereo *ts = stereo.GetTetrahedralStereo(center);
    OBTetrahedralStereo::Config cfg = ts->GetConfig();
    OB_ASSERT( cfg.specified == true );
    OB_ASSERT( cfg.winding == OBStereo::UnknownWinding );
  }

  return retval;
}

void test_noStereo(const std::string &file)
{
  OBMol mol;
  OBConversion conv;
  conv.SetInFormat("mol");
  conv.SetOutFormat("can");

  const std::string filename = GetFilename(file);
  ifstream ifs(filename.c_str());
  if (!ifs) {
    cerr << "Could not open " << filename << endl;
    return;
  }

  conv.Read(&mol, &ifs);
  ifs.close();

  // perceive stereochemistry from 3D
  StereoFrom3D(&mol);

  std::vector<OBGenericData *> stereoData = mol.GetAllData(OBGenericDataType::StereoData);
  OB_ASSERT( stereoData.size() == 0 );
} 


int main()
{
  // Define location of file formats for testing
  #ifdef FORMATDIR
    char env[BUFF_SIZE];
    snprintf(env, BUFF_SIZE, "BABEL_LIBDIR=%s", FORMATDIR);
    putenv(env);
  #endif

  // There are several cases to test:

  //////////////////////////////////////////////////////////////////////////////
  // 
  // 1      StereoFrom3D for tetrahedral atoms
  //
  //////////////////////////////////////////////////////////////////////////////

  // 1.1    Input molecule with 4 refs (explicit H)
 
  // 1.1.1  234 == 234  @
  string smiles3D_1 = test_singleTetrahedral("stereo/tetrahedral3D_1.mol",
      OBTetrahedralStereo::Config(1, 0, OBStereo::MakeRefs(2, 3, 4), OBStereo::AntiClockwise));
  // 1.1.2  234 == 234  @@
  string smiles3D_2 = test_singleTetrahedral("stereo/tetrahedral3D_2.mol",
      OBTetrahedralStereo::Config(1, 0, OBStereo::MakeRefs(2, 3, 4), OBStereo::Clockwise));
  
  // 1.1.3 compare using ImplicitRef
  string smiles3D_5 = test_singleTetrahedral("stereo/tetrahedral3D_1.mol",
      OBTetrahedralStereo::Config(1, 0, OBStereo::MakeRefs(2, 3, OBStereo::ImplicitRef), OBStereo::AntiClockwise));
   
  // 1.2    Input molecule with 3 refs (implicit H)
  
  // 1.2.1  23H == 23H  @@
  string smiles3D_3 = test_singleTetrahedral("stereo/tetrahedral3D_3.mol",
      OBTetrahedralStereo::Config(0, 1, OBStereo::MakeRefs(2, 3, OBStereo::ImplicitRef), OBStereo::Clockwise));
  // 1.2.2  23H == 23H  @
  string smiles3D_4 = test_singleTetrahedral("stereo/tetrahedral3D_4.mol",
      OBTetrahedralStereo::Config(0, 1, OBStereo::MakeRefs(2, 3, OBStereo::ImplicitRef), OBStereo::AntiClockwise));
  
  // 1.2.3 compare using explicit id 
  string smiles3D_6 = test_singleTetrahedral("stereo/tetrahedral3D_4.mol",
      OBTetrahedralStereo::Config(0, 1, OBStereo::MakeRefs(2, 3, 4), OBStereo::AntiClockwise));

  OB_ASSERT( smiles3D_1 == smiles3D_2 );
  OB_ASSERT( smiles3D_1 == smiles3D_3 );
  OB_ASSERT( smiles3D_1 == smiles3D_4 );
  OB_ASSERT( smiles3D_1 == smiles3D_5 );
  OB_ASSERT( smiles3D_1 == smiles3D_6 );
  
  cout << smiles3D_1 << endl;

  //////////////////////////////////////////////////////////////////////////////
  // 
  // 2      StereoFrom3D for cis/trans bonds
  //
  //////////////////////////////////////////////////////////////////////////////

  // 2.1    Input molecule with 4 refs (explicit H)
 
  // 2.1.1  F      H   2      5
  //         \    /
  //          C==C       1  3
  //         /    \     
  //        H      F   0      4
  string smiles7 = test_singleCisTrans("stereo/cistrans3D_1.mol",
      OBCisTransStereo::Config(1, 3, OBStereo::MakeRefs(0, 2, 4, 5), OBStereo::ShapeZ));
  // implicit refs
  test_singleCisTrans("stereo/cistrans3D_1.mol", OBCisTransStereo::Config(1, 3, 
        OBStereo::MakeRefs(OBStereo::ImplicitRef, 2, 4, 5), OBStereo::ShapeZ));
  test_singleCisTrans("stereo/cistrans3D_1.mol", OBCisTransStereo::Config(1, 3, 
        OBStereo::MakeRefs(0, 2, 4, OBStereo::ImplicitRef), OBStereo::ShapeZ));
  test_singleCisTrans("stereo/cistrans3D_1.mol", OBCisTransStereo::Config(1, 3, 
        OBStereo::MakeRefs(OBStereo::ImplicitRef, 2, 4, OBStereo::ImplicitRef), OBStereo::ShapeZ));
 
  // 2.1.2  F      F   2      4
  //         \    /
  //          C==C       1  3
  //         /    \     
  //        H      H   0      5
  string smiles8 = test_singleCisTrans("stereo/cistrans3D_4.mol",
      OBCisTransStereo::Config(1, 3, OBStereo::MakeRefs(0, 2, 4, 5), OBStereo::ShapeU));
  // implicit refs
  test_singleCisTrans("stereo/cistrans3D_4.mol", OBCisTransStereo::Config(1, 3, 
        OBStereo::MakeRefs(OBStereo::ImplicitRef, 2, 4, 5), OBStereo::ShapeU));
  test_singleCisTrans("stereo/cistrans3D_4.mol", OBCisTransStereo::Config(1, 3, 
        OBStereo::MakeRefs(0, 2, 4, OBStereo::ImplicitRef), OBStereo::ShapeU));
  test_singleCisTrans("stereo/cistrans3D_4.mol", OBCisTransStereo::Config(1, 3, 
        OBStereo::MakeRefs(OBStereo::ImplicitRef, 2, 4, OBStereo::ImplicitRef), OBStereo::ShapeU));
 
  // 2.2    Input molecule with 3 refs (implicit H)

  // 2.2.1  F          1      
  //         \    
  //          C==C       0  2
  //              \     
  //               F          3
  string smiles9 = test_singleCisTrans("stereo/cistrans3D_2.mol", OBCisTransStereo::Config(0, 2, 
      OBStereo::MakeRefs(1, OBStereo::ImplicitRef, 3, OBStereo::ImplicitRef), OBStereo::ShapeU));

  // 2.2.1  F      F   2      4
  //         \    /
  //          C==C       1  3
  //         /
  //        H          0
  string smiles10 = test_singleCisTrans("stereo/cistrans3D_5.mol", OBCisTransStereo::Config(1, 3, 
      OBStereo::MakeRefs(0, 2, 4, OBStereo::ImplicitRef), OBStereo::ShapeU));

  OB_ASSERT( smiles7 == smiles9 );
  OB_ASSERT( smiles8 == smiles10 );
  
  cout << smiles7 << endl;

  // 3      No stereochemistry

  // 3.1    H      H   3      2
  //         \    /
  //          C==C       0  1
  //              \     
  //               F          4
  test_noStereo("stereo/cistrans3D_3.mol");

  //////////////////////////////////////////////////////////////////////////////
  // 
  // 3      StereoFrom2D for tetrahedral atoms
  //
  //////////////////////////////////////////////////////////////////////////////
  
  // 3.1    Input molecule with 4 refs (2x in plane, one behind plane, one in front of plane)

  // 3.1.1  Input molecule with 4 refs (2x in plane bond, hash & wedge bond)
  //
  //         F   I         C-F : hash (from C to F)
  //          \ /          C-I : wedge (from C to F)
  //           C
  //          / \
  //        Br   Cl
  string smiles2D_1 = test_singleTetrahedral("stereo/tetrahedral2D_1.mol",
      OBTetrahedralStereo::Config(1, 0, OBStereo::MakeRefs(2, 3, 4), OBStereo::AntiClockwise));
 
  // 3.1.2  Input molecule with 4 refs (2x in plane bond, 'real' hash & 'inverted' hash bond)
  //
  //         F   I         C-F : 'real' hash (from C to F)
  //          \ /          C-I : 'inverted' hash (from I to C) = 'real' wedge
  //           C
  //          / \
  //        Br   Cl
  string smiles2D_2 = test_singleTetrahedral("stereo/tetrahedral2D_4.mol",
      OBTetrahedralStereo::Config(1, 0, OBStereo::MakeRefs(2, 3, 4), OBStereo::AntiClockwise));
 
  // 3.1.3  Input molecule with 4 refs (2x in plane bond, 'real' wedge & 'inverted' wedge bond)
  //
  //         F   I         C-F : 'inverted' wedge (from F to C) = 'real' hash
  //          \ /          C-I : 'real' wedge (from C to I)
  //           C
  //          / \
  //        Br   Cl
  string smiles2D_3 = test_singleTetrahedral("stereo/tetrahedral2D_5.mol",
      OBTetrahedralStereo::Config(1, 0, OBStereo::MakeRefs(2, 3, 4), OBStereo::AntiClockwise));

  OB_ASSERT( smiles2D_1 == smiles2D_2 );
  OB_ASSERT( smiles2D_1 == smiles2D_3 );
  
  cout << smiles2D_1 << endl;
  
  // 3.2    Input molecule with 3 refs (2x in plane, one behind plane or in front of plane)
 
  // 3.2.1  Input molecule with 3 refs (2x in plane bond, real wedge bond)
  //
  //           F           C-F : 'real' wedge (from C to F)
  //           |           
  //           C
  //          / \
  //        Br   Cl
  string smiles2D_4 = test_singleTetrahedral("stereo/tetrahedral2D_2.mol",
      OBTetrahedralStereo::Config(1, 0, OBStereo::MakeRefs(2, 3, OBStereo::ImplicitRef), 
      OBStereo::AntiClockwise));
 
  // 3.2.2  Input molecule with 3 refs (2x in plane bond, inverted wedge bond)
  //
  //           F           C-F : 'inverted' wedge (from F to C)
  //           |           
  //           C
  //          / \
  //        Br   Cl
  
  /* No support for inverted bonds in 'tip-only' convention
  string smiles2D_5 = test_singleTetrahedral("stereo/tetrahedral2D_6.mol",
      OBTetrahedralStereo::Config(1, 0, OBStereo::MakeRefs(2, 3, OBStereo::ImplicitRef), 
      OBStereo::Clockwise));*/
 
  // 3.2.3  Input molecule with 3 refs (2x in plane bond, real hash bond)
  //
  //           F           C-F : 'real' hash (from C to F)
  //           |           
  //           C
  //          / \
  //        Br   Cl
  string smiles2D_6 = test_singleTetrahedral("stereo/tetrahedral2D_3.mol",
      OBTetrahedralStereo::Config(1, 0, OBStereo::MakeRefs(2, 3, OBStereo::ImplicitRef), 
      OBStereo::Clockwise));
 
  // 3.2.3  Input molecule with 3 refs (2x in plane bond, real hash bond)
  //
  //           F           C-F : 'inverted' hash (from F to C)
  //           |           
  //           C
  //          / \
  //        Br   Cl
  /* No support for inverted bonds in 'tip-only' convention
  string smiles2D_7 = test_singleTetrahedral("stereo/tetrahedral2D_7.mol",
      OBTetrahedralStereo::Config(1, 0, OBStereo::MakeRefs(2, 3, OBStereo::ImplicitRef), 
      OBStereo::AntiClockwise));

  OB_ASSERT( smiles2D_4 == smiles2D_7 );
  OB_ASSERT( smiles2D_5 == smiles2D_6 );*/

  cout << smiles2D_4 << endl;

 // Input molecule with unknown stereochemistry

  string smiles2D_8 = test_singleUnknownTetrahedral("stereo/tetrahedral2D_8.mol", 1);
  cout << smiles2D_8 << endl;

 // Input molecule with unspecified stereochemistry

  string smiles2D_9 = test_singleUnspecifiedTetrahedral("stereo/tetrahedral2D_9.mol", 1);
  cout << smiles2D_9 << endl;

  //////////////////////////////////////////////////////////////////////////////
  // 
  // 4      StereoFrom2D for cis/trans bonds
  //
  //////////////////////////////////////////////////////////////////////////////
 
  //  C     
  //   \    
  //    C==C 
  //        \
  //         C
  string cistrans2D_1 = test_singleCisTrans("stereo/cistrans2D_1.mol",
      OBCisTransStereo::Config(0, 1, OBStereo::MakeRefs(2, OBStereo::ImplicitRef, 
      3, OBStereo::ImplicitRef), OBStereo::ShapeU));

  //  
  //  C      C
  //   \    / 
  //    C==C
  //        
  string cistrans2D_2 = test_singleCisTrans("stereo/cistrans2D_2.mol",
      OBCisTransStereo::Config(0, 1, OBStereo::MakeRefs(2, OBStereo::ImplicitRef, 
      OBStereo::ImplicitRef, 3), OBStereo::ShapeU));



  return 0;
}

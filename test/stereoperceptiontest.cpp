#include "obtest.h"
#include <openbabel/stereo/tetrahedral.h>
#include <openbabel/stereo/cistrans.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>

using namespace std;
using namespace OpenBabel;

std::string GetFilename(const std::string &filename)
{
  #ifdef TESTDATADIR
    string testdatadir = TESTDATADIR;
    string path = testdatadir + filename;
  #else
    string path = "files/" + filename;
  #endif

  return path;
}

std::string test_singleTetrahedral(const std::string &file, const OBTetrahedralStereo::Config &correct)
{
  OBMol mol;
  OBConversion conv;
  conv.SetInFormat("sdf");
  conv.SetOutFormat("can");

  const std::string filename = GetFilename(file);
  ifstream ifs(filename.c_str());
  if (!ifs) {
    cerr << "Could not open " << filename << endl;
    return std::string();
  }

  conv.Read(&mol, &ifs);
  ifs.close();

  // perceive stereochemistry from 3D
  StereoFrom3D(&mol);

  std::vector<OBGenericData *> stereoData = mol.GetAllData(OBGenericDataType::StereoData);
  OB_ASSERT( stereoData.size() == 1 );

  // compare the stereochemistry
  for (std::vector<OBGenericData*>::iterator data = stereoData.begin(); data != stereoData.end(); ++data) {
    if (((OBStereoBase*)*data)->GetType() == OBStereo::Tetrahedral) {
      OBTetrahedralStereo *ts = dynamic_cast<OBTetrahedralStereo*>(*data);
      OB_ASSERT( ts->GetConfig() == correct );
      if ( ts->GetConfig() != correct ) {
        cout << ts->GetConfig() << endl;
        cout << correct << endl;
      }

      OBTetrahedralStereo::Config cfg = ts->GetConfig();
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


int main()
{
  // There are several cases to test:

  // 1      StereoFrom3D for tetrahedral atoms

  // 1.1    Input molecule with 4 refs (explicit H)
 
  // 1.1.1  234 == 234  @
  string smiles1 = test_singleTetrahedral("tetrahedral1.sdf",
      OBTetrahedralStereo::Config(1, 0, OBStereo::MakeRefs(2, 3, 4), OBStereo::AntiClockwise));
  // 1.1.2  234 == 234  @@
  string smiles2 = test_singleTetrahedral("tetrahedral2.sdf",
      OBTetrahedralStereo::Config(1, 0, OBStereo::MakeRefs(2, 3, 4), OBStereo::Clockwise));
  
  // 1.1.3 compare using ImplicitId
  string smiles5 = test_singleTetrahedral("tetrahedral1.sdf",
      OBTetrahedralStereo::Config(1, 0, OBStereo::MakeRefs(2, 3, OBStereo::ImplicitId), OBStereo::AntiClockwise));
  
  
  // 1.2    Input molecule with 3 refs (implicit H)
  
  // 1.2.1  23H == 23H  @@
  string smiles3 = test_singleTetrahedral("tetrahedral3.sdf",
      OBTetrahedralStereo::Config(0, 1, OBStereo::MakeRefs(2, 3, OBStereo::ImplicitId), OBStereo::Clockwise));
  // 1.2.2  23H == 23H  @
  string smiles4 = test_singleTetrahedral("tetrahedral4.sdf",
      OBTetrahedralStereo::Config(0, 1, OBStereo::MakeRefs(2, 3, OBStereo::ImplicitId), OBStereo::AntiClockwise));
  
  // 1.2.3 compare using explicit id 
  string smiles6 = test_singleTetrahedral("tetrahedral4.sdf",
      OBTetrahedralStereo::Config(0, 1, OBStereo::MakeRefs(2, 3, 4), OBStereo::AntiClockwise));
  

  OB_ASSERT( smiles1 == smiles2 );
  OB_ASSERT( smiles1 == smiles3 );
  OB_ASSERT( smiles1 == smiles4 );
  OB_ASSERT( smiles1 == smiles5 );
  OB_ASSERT( smiles1 == smiles6 );



  return 0;
}

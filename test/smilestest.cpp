#include "obtest.h"

#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/stereo/tetrahedral.h>
#include <openbabel/stereo/cistrans.h>
#include <openbabel/stereo/squareplanar.h>

#include <openbabel/graphsym.h>
#include <openbabel/canon.h>
#include <openbabel/atom.h>

using namespace std;
using namespace OpenBabel;

/*
 * Stereo classes have their own tests. This file tests if the smiles
 * format uses them correctly.
 */

void testTetrahedralStereo1()
{
  cout << "testTetrahedralStereo1()" << endl;
  // read a smiles string
  OBMol mol;
  OBConversion conv;
  OB_REQUIRE( conv.SetInFormat("smi") );
  cout << "smiles: C[C@H](O)N" << endl;
  OB_REQUIRE( conv.ReadString(&mol, "C[C@H](O)N") );

  // get the stereo data
  OB_REQUIRE( mol.HasData(OBGenericDataType::StereoData) );
  std::vector<OBGenericData *> stereoData = mol.GetAllData(OBGenericDataType::StereoData);
  OB_REQUIRE( stereoData.size() == 1 );

  // convert to tetrahedral data
  OB_REQUIRE( ((OBStereoBase*)stereoData[0])->GetType() == OBStereo::Tetrahedral );
  OBTetrahedralStereo *ts = dynamic_cast<OBTetrahedralStereo*>(stereoData[0]);
  OB_REQUIRE( ts );

  // print the configuration
  cout << *ts << endl;

  // construct a valid configuration here
  //
  // C[C@H](O)N
  // 0 1 2  3 4  <- ids
  //
  OBTetrahedralStereo::Config cfg(1, 0, OBStereo::MakeRefs(4, 3, 2), OBStereo::Clockwise);

  // compare stereochemistry
  OB_REQUIRE( ts->GetConfig() == cfg );

  cout << endl;
}

void genericSmilesCanonicalTest(const std::string &smiles)
{
  cout << "Testing generic smiles <-> canonical smiles" << endl;
  // read a smiles string
  OBMol mol;
  OBConversion conv;
  OB_REQUIRE( conv.SetInFormat("smi") );
  OB_REQUIRE( conv.SetOutFormat("can") );
  cout << "smiles: " << smiles << endl;
  // read a smiles string
  OB_REQUIRE( conv.ReadString(&mol, smiles) );

  // store the stereo data for the smiles string using unique symmetry ids
  std::vector<OBTetrahedralStereo::Config> tetrahedral1;
  std::vector<OBCisTransStereo::Config> cistrans1;
  std::vector<OBSquarePlanarStereo::Config> squareplanar1;

  // get the stereo data
  OB_ASSERT( mol.HasData(OBGenericDataType::StereoData) );
  std::vector<OBGenericData *> stereoData = mol.GetAllData(OBGenericDataType::StereoData);

  std::vector<unsigned int> canlbls;
  std::vector<unsigned int> symclasses;
  OBGraphSym gs1(&mol);
  gs1.GetSymmetry(symclasses);
  CanonicalLabels(&mol, symclasses, canlbls);
  cout << "mol.NumAtoms = " << mol.NumAtoms() << endl;
  for (std::vector<OBGenericData*>::iterator data = stereoData.begin(); data != stereoData.end(); ++data) {
    if (((OBStereoBase*)*data)->GetType() == OBStereo::Tetrahedral) {
      // convert to tetrahedral data
      OBTetrahedralStereo *ts = dynamic_cast<OBTetrahedralStereo*>(*data);
      OB_REQUIRE( ts );
      OB_ASSERT( ts->IsValid() );
      if (!ts->IsValid())
        continue;

      OBTetrahedralStereo::Config config = ts->GetConfig();
      // convert atom ids to symmetry ids
     if (mol.GetAtomById(config.center))
        config.center = canlbls.at( mol.GetAtomById(config.center)->GetIdx() - 1 );
      if (mol.GetAtomById(config.from))
        config.from = canlbls.at( mol.GetAtomById(config.from)->GetIdx() - 1 );
      if (mol.GetAtomById(config.refs[0]))
        config.refs[0] = canlbls.at( mol.GetAtomById(config.refs[0])->GetIdx() - 1 );
      if (mol.GetAtomById(config.refs[1]))
        config.refs[1] = canlbls.at( mol.GetAtomById(config.refs[1])->GetIdx() - 1 );
      if (mol.GetAtomById(config.refs[2]))
        config.refs[2] = canlbls.at( mol.GetAtomById(config.refs[2])->GetIdx() - 1 );
      cout << "Config with symmetry ids: " << config << endl;
      tetrahedral1.push_back(config);
    } else
    if (((OBStereoBase*)*data)->GetType() == OBStereo::CisTrans) {
      // convert to tetrahedral data
      OBCisTransStereo *ct = dynamic_cast<OBCisTransStereo*>(*data);
      OB_REQUIRE( ct );
      OB_ASSERT( ct->IsValid() );

      OBCisTransStereo::Config config = ct->GetConfig();
      // convert atom ids to symmetry ids
      config.begin = canlbls.at( mol.GetAtomById(config.begin)->GetIdx() - 1 );
      config.end = canlbls.at( mol.GetAtomById(config.end)->GetIdx() - 1 );
      if (mol.GetAtomById(config.refs[0]))
        config.refs[0] = canlbls.at( mol.GetAtomById(config.refs[0])->GetIdx() - 1 );
      if (mol.GetAtomById(config.refs[1]))
        config.refs[1] = canlbls.at( mol.GetAtomById(config.refs[1])->GetIdx() - 1 );
      if (mol.GetAtomById(config.refs[2]))
        config.refs[2] = canlbls.at( mol.GetAtomById(config.refs[2])->GetIdx() - 1 );
      if (mol.GetAtomById(config.refs[3]))
        config.refs[3] = canlbls.at( mol.GetAtomById(config.refs[3])->GetIdx() - 1 );
      cout << "Config with symmetry ids: " << config << endl;
      cistrans1.push_back(config);
    } else
    if (((OBStereoBase*)*data)->GetType() == OBStereo::SquarePlanar) {
      // convert to tetrahedral data
      OBSquarePlanarStereo *sp = dynamic_cast<OBSquarePlanarStereo*>(*data);
      OB_REQUIRE( sp );
      OB_ASSERT( sp->IsValid() );
      if (!sp->IsValid())
        continue;

      OBSquarePlanarStereo::Config config = sp->GetConfig();
      // convert atom ids to symmetry ids
     if (mol.GetAtomById(config.center))
        config.center = canlbls.at( mol.GetAtomById(config.center)->GetIdx() - 1 );
      if (mol.GetAtomById(config.refs[0]))
        config.refs[0] = canlbls.at( mol.GetAtomById(config.refs[0])->GetIdx() - 1 );
      if (mol.GetAtomById(config.refs[1]))
        config.refs[1] = canlbls.at( mol.GetAtomById(config.refs[1])->GetIdx() - 1 );
      if (mol.GetAtomById(config.refs[2]))
        config.refs[2] = canlbls.at( mol.GetAtomById(config.refs[2])->GetIdx() - 1 );
      if (mol.GetAtomById(config.refs[3]))
        config.refs[3] = canlbls.at( mol.GetAtomById(config.refs[3])->GetIdx() - 1 );
      cout << "Config with symmetry ids: " << config << endl;
      squareplanar1.push_back(config);
    }


  }

  // write to can smiles
  std::string canSmiles = conv.WriteString(&mol);
  cout << "canSmiles: " << canSmiles;
  // read can smiles in again
  OB_REQUIRE( conv.ReadString(&mol, canSmiles) );

  // store the stereo data for the smiles string using unique symmetry ids
  std::vector<OBTetrahedralStereo::Config> tetrahedral2;
  std::vector<OBCisTransStereo::Config> cistrans2;
  std::vector<OBSquarePlanarStereo::Config> squareplanar2;

  // get the stereo data
  OB_ASSERT( mol.HasData(OBGenericDataType::StereoData) );
  stereoData = mol.GetAllData(OBGenericDataType::StereoData);

  OBGraphSym gs2(&mol);
  gs2.GetSymmetry(symclasses);
  CanonicalLabels(&mol, symclasses, canlbls);
  cout << "mol.NumAtoms = " << mol.NumAtoms() << endl;
  for (std::vector<OBGenericData*>::iterator data = stereoData.begin(); data != stereoData.end(); ++data) {
    if (((OBStereoBase*)*data)->GetType() == OBStereo::Tetrahedral) {
      // convert to tetrahedral data
      OBTetrahedralStereo *ts = dynamic_cast<OBTetrahedralStereo*>(*data);
      OB_REQUIRE( ts );
      OB_ASSERT( ts->IsValid() );

      OBTetrahedralStereo::Config config = ts->GetConfig();
      // convert atom ids to symmetry ids
      if (mol.GetAtomById(config.center))
        config.center = canlbls.at( mol.GetAtomById(config.center)->GetIdx() - 1 );
      if (mol.GetAtomById(config.from))
        config.from = canlbls.at( mol.GetAtomById(config.from)->GetIdx() - 1 );
      if (mol.GetAtomById(config.refs[0]))
        config.refs[0] = canlbls.at( mol.GetAtomById(config.refs[0])->GetIdx() - 1 );
      if (mol.GetAtomById(config.refs[1]))
        config.refs[1] = canlbls.at( mol.GetAtomById(config.refs[1])->GetIdx() - 1 );
      if (mol.GetAtomById(config.refs[2]))
        config.refs[2] = canlbls.at( mol.GetAtomById(config.refs[2])->GetIdx() - 1 );
      cout << "Config with symmetry ids: " << config << endl;
      tetrahedral2.push_back(config);
    }
    if (((OBStereoBase*)*data)->GetType() == OBStereo::CisTrans) {
      // convert to tetrahedral data
      OBCisTransStereo *ct = dynamic_cast<OBCisTransStereo*>(*data);
      OB_REQUIRE( ct );
      OB_ASSERT( ct->IsValid() );

      OBCisTransStereo::Config config = ct->GetConfig();
      // convert atom ids to symmetry ids
      config.begin = canlbls.at( mol.GetAtomById(config.begin)->GetIdx() - 1 );
      config.end = canlbls.at( mol.GetAtomById(config.end)->GetIdx() - 1 );
      if (mol.GetAtomById(config.refs[0]))
        config.refs[0] = canlbls.at( mol.GetAtomById(config.refs[0])->GetIdx() - 1 );
      if (mol.GetAtomById(config.refs[1]))
        config.refs[1] = canlbls.at( mol.GetAtomById(config.refs[1])->GetIdx() - 1 );
      if (mol.GetAtomById(config.refs[2]))
        config.refs[2] = canlbls.at( mol.GetAtomById(config.refs[2])->GetIdx() - 1 );
      if (mol.GetAtomById(config.refs[3]))
        config.refs[3] = canlbls.at( mol.GetAtomById(config.refs[3])->GetIdx() - 1 );
      cout << "Config with symmetry ids: " << config << endl;
      cistrans2.push_back(config);
    } else
    if (((OBStereoBase*)*data)->GetType() == OBStereo::SquarePlanar) {
      // convert to tetrahedral data
      OBSquarePlanarStereo *sp = dynamic_cast<OBSquarePlanarStereo*>(*data);
      OB_REQUIRE( sp );
      OB_ASSERT( sp->IsValid() );

      OBSquarePlanarStereo::Config config = sp->GetConfig();
      // convert atom ids to symmetry ids
      if (mol.GetAtomById(config.center))
        config.center = canlbls.at( mol.GetAtomById(config.center)->GetIdx() - 1 );
      if (mol.GetAtomById(config.refs[0]))
        config.refs[0] = canlbls.at( mol.GetAtomById(config.refs[0])->GetIdx() - 1 );
      if (mol.GetAtomById(config.refs[1]))
        config.refs[1] = canlbls.at( mol.GetAtomById(config.refs[1])->GetIdx() - 1 );
      if (mol.GetAtomById(config.refs[2]))
        config.refs[2] = canlbls.at( mol.GetAtomById(config.refs[2])->GetIdx() - 1 );
      if (mol.GetAtomById(config.refs[3]))
        config.refs[3] = canlbls.at( mol.GetAtomById(config.refs[3])->GetIdx() - 1 );
      cout << "Config with symmetry ids: " << config << endl;
      squareplanar2.push_back(config);
    }

  }

  // compare the tetrahedral structs
  OB_ASSERT( tetrahedral1.size() == tetrahedral2.size() );
  for (unsigned int i = 0; i < tetrahedral1.size(); ++i) {
    for (unsigned int j = 0; j < tetrahedral2.size(); ++j) {
      if (tetrahedral1[i].center == tetrahedral2[j].center)
        OB_ASSERT( tetrahedral1[i] == tetrahedral2[j] );
        if ( tetrahedral1[i] != tetrahedral2[j] ) {
          cout << "1 = " << tetrahedral1[i] << endl;
          cout << "2 = " << tetrahedral2[j] << endl;
        }
    }
  }
  // compare the cistrans structs
  OB_ASSERT( cistrans1.size() == cistrans2.size() );
  for (unsigned int i = 0; i < cistrans1.size(); ++i) {
    for (unsigned int j = 0; j < cistrans2.size(); ++j) {
      if ((cistrans1[i].begin == cistrans2[j].begin) && (cistrans1[i].end == cistrans2[j].end))
        OB_ASSERT( cistrans1[i] == cistrans2[j] );
      if ((cistrans1[i].begin == cistrans2[j].end) && (cistrans1[i].end == cistrans2[j].begin))
        OB_ASSERT( cistrans1[i] == cistrans2[j] );
    }
  }
  // compare the square-planar structs
  OB_ASSERT( squareplanar1.size() == squareplanar2.size() );
  for (unsigned int i = 0; i < squareplanar1.size(); ++i) {
    for (unsigned int j = 0; j < squareplanar2.size(); ++j) {
      if (squareplanar1[i].center == squareplanar2[j].center)
        OB_ASSERT( squareplanar1[i] == squareplanar2[j] );
        if ( squareplanar1[i] != squareplanar2[j] ) {
          cout << "1 = " << squareplanar1[i] << endl;
          cout << "2 = " << squareplanar2[j] << endl;
        }
    }
  }

  cout << "." << endl << endl;
}


int smilestest(int argc, char* argv[])
{
  int defaultchoice = 1;
  
  int choice = defaultchoice;

  if (argc > 1) {
    if(sscanf(argv[1], "%d", &choice) != 1) {
      printf("Couldn't parse that input as a number\n");
      return -1;
    }
  }

  // Define location of file formats for testing
  #ifdef FORMATDIR
    char env[BUFF_SIZE];
    snprintf(env, BUFF_SIZE, "BABEL_LIBDIR=%s", FORMATDIR);
    putenv(env);
  #endif

  switch(choice) {
  case 1:
    testTetrahedralStereo1();
    // Tetrahedral
    genericSmilesCanonicalTest("C[C@H](O)N");
    genericSmilesCanonicalTest("Cl[C@@](CCl)(I)Br");
    // CisTrans
    genericSmilesCanonicalTest("Cl/C=C/F");
    break;
  case 2:
    // SquarePlanar
    genericSmilesCanonicalTest("F[Po@SP1](Cl)(Br)I");
    genericSmilesCanonicalTest("F[Po@SP2](Br)(Cl)I");
    genericSmilesCanonicalTest("F[Po@SP3](Cl)(I)Br");
    break;
  case 3:
    // Mixed
    genericSmilesCanonicalTest("CCC[C@@H](O)CC\\C=C\\C=C\\C#CC#C\\C=C\\CO");
    genericSmilesCanonicalTest("OC[C@@H](O1)[C@@H](O)[C@H](O)[C@@H](O)[C@@H](O)1");
    genericSmilesCanonicalTest("OC[C@@H](O1)[C@@H](O)[C@H](O)[C@@H]2[C@@H]1c3c(O)c(OC)c(O)cc3C(=O)O2");
    genericSmilesCanonicalTest("CC(=O)OCCC(/C)=C\\C[C@H](C(C)=C)CCC=C");
    genericSmilesCanonicalTest("CC[C@H](O1)CC[C@@]12CCCO2");
    genericSmilesCanonicalTest("CN1CCC[C@H]1c2cccnc2");
    genericSmilesCanonicalTest("CC(C)[C@@]12C[C@@H]1[C@@H](C)C(=O)C2");
    genericSmilesCanonicalTest("CC(C)[C@H]1CC[C@]([C@@H]2[C@@H]1C=C(COC2=O)C(=O)O)(CCl)O");
    genericSmilesCanonicalTest("C(CS[14CH2][14C@@H]1[14C@H]([14C@H]([14CH](O1)O)O)O)[C@@H](C(=O)O)N");
    genericSmilesCanonicalTest("CCC[C@@H]1C[C@H](N(C1)C)C(=O)NC([C@@H]2[C@@H]([C@@H]([C@H]([C@H](O2)SC)OP(=O)(O)O)O)O)C(C)Cl");
    break;
    // this structure fails but this is an error in genericSmilesCanonicalTest (It should sort the centers). However, this is
    // extensivily tested in other tests.
    //genericSmilesCanonicalTest("O1C=C[C@H]([C@H]1O2)c3c2cc(OC)c4c3OC(=O)C5=C4CCC(=O)5");

    // FAILING: need to fix graphsymtest first!!
    // ring gets converted to aromatic ring, adding H on n (i.e. N -> [nH])
    //genericSmilesCanonicalTest("CC1=CN(C(=O)NC1=O)[C@H]2C[C@@H]([C@H](O2)CNCC3=CC=CC=C3)O");
  default:
    cout << "Test number " << choice << " does not exist!\n";
    return -1;
  }

  return 0;
}



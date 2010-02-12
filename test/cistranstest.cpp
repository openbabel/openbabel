#include "obtest.h"

#include <iostream>

#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/stereo/cistrans.h>

using namespace std;
using namespace OpenBabel;

void test_GetType ()
{
  OBMol mol;
  OBCisTransStereo ct(&mol);
  OB_ASSERT( ct.GetType() == OBStereo::CisTrans );
}

void test_configStruct()
{
  // reference Config
  OBCisTransStereo::Config reference(0, 1, OBStereo::MakeRefs(2, 3, 4, 5), OBStereo::ShapeU);
  
  OB_ASSERT( reference.begin == 0 );
  OB_ASSERT( reference.end == 1 );
  OB_ASSERT( reference.refs.size() == 4 );
  OB_ASSERT( reference.shape == OBStereo::ShapeU );

  // test copying
  OBCisTransStereo::Config referenceCopy = reference;
  OB_ASSERT( reference == referenceCopy );

  // invalid begin id
  OBCisTransStereo::Config invalidBegin(45, 1, OBStereo::MakeRefs(2, 3, 4, 5), OBStereo::ShapeU);
  OB_ASSERT( reference != invalidBegin );

  // invalid end id
  OBCisTransStereo::Config invalidEnd(0, 45, OBStereo::MakeRefs(2, 3, 4, 5), OBStereo::ShapeU);
  OB_ASSERT( reference != invalidEnd );

  // test other refs
  OBCisTransStereo::Config cfg1(0, 1, OBStereo::MakeRefs(2, 4, 3, 5), OBStereo::ShapeU);
  OB_ASSERT( reference != cfg1 );

  // test Z shape == U shape
  OBCisTransStereo::Config cfg2(0, 1, OBStereo::MakeRefs(3, 2, 4, 5), OBStereo::ShapeZ);
  OB_ASSERT( reference == cfg2 );

  // test 4 shape == U shape
  OBCisTransStereo::Config cfg3(0, 1, OBStereo::MakeRefs(2, 4, 3, 5), OBStereo::Shape4);
  OB_ASSERT( reference == cfg3 );

  // test Z shape == U shape
  OB_ASSERT( cfg2 == cfg3 );

}

void test_IsValid()
{
  OBMol mol;
  OBCisTransStereo ct(&mol);
  OBCisTransStereo::Config cfg, cfgCopy;
  cfg.begin = 0;
  cfg.end = 1;
  cfg.refs = OBStereo::MakeRefs(2, 3, 4, 5);

  ct.SetConfig(cfg);
  OB_ASSERT( ct.IsValid() );

  // no begin
  cfgCopy = cfg;
  cfgCopy.begin = OBStereo::NoRef;
  ct.SetConfig(cfgCopy);
  OB_ASSERT( !ct.IsValid() );

  // no end
  cfgCopy = cfg;
  cfgCopy.end = OBStereo::NoRef;
  ct.SetConfig(cfgCopy);
  OB_ASSERT( !ct.IsValid() );

  // no refs
  cfgCopy = cfg;
  cfgCopy.refs = std::vector<unsigned long>();
  ct.SetConfig(cfgCopy);
  OB_ASSERT( !ct.IsValid() );
}

void test_equalsOperator()
{
  OBCisTransStereo ct1(0), ct2(0);
  OBCisTransStereo::Config cfg;
  cfg.begin = 0;
  cfg.end = 1;
  cfg.refs = OBStereo::MakeRefs(2, 3, 4, 5);

  ct1.SetConfig(cfg);
  ct2.SetConfig(cfg);
  OB_ASSERT( ct1 == ct2 );

  cfg.shape = OBStereo::ShapeZ;
  ct2.SetConfig(cfg);
  OB_ASSERT( ct1 != ct2 );

  OBStereo::Permutate(cfg.refs, 0, 1);
  ct2.SetConfig(cfg);
  OB_ASSERT( ct1 == ct2 );
}

void test_GetSetConfig()
{
  OBCisTransStereo ct(0);
  OBCisTransStereo::Config cfg;

  // set clockwise, viewing from 1
  OB_ASSERT( !ct.IsValid() );
  cfg.refs = OBStereo::MakeRefs(2, 3, 4, 5);
  cfg.begin = 0;
  cfg.end = 1;
  ct.SetConfig(cfg);
  OB_ASSERT( ct.IsValid() );

  OBCisTransStereo::Config cfg2 = ct.GetConfig();
  OB_ASSERT( cfg2.begin == 0 );
  OB_ASSERT( cfg2.end == 1 );
  OB_ASSERT( cfg2.refs.size() == 4 );
  OB_ASSERT( cfg2.refs[0] == 2 );
  OB_ASSERT( cfg2.refs[1] == 3 );
  OB_ASSERT( cfg2.refs[2] == 4 );
  OB_ASSERT( cfg2.refs[3] == 5 );
  OB_ASSERT( cfg2.shape == OBStereo::ShapeU );
  OB_ASSERT( cfg == cfg2 );

}

void testRefs()
{
  OBMol mol;
  OBCisTransStereo ct(&mol);

  // 1   4     1   4     1---4     1   4
  //  \ /      |   |        /       \ /
  //   C    =  | C |  =    C    =    C
  //  / \      |   |      /         / \
  // 2   3     2---3     2---3     2---3

  // set refs using default U shape
  OBCisTransStereo::Config cfg(0, 1, OBStereo::MakeRefs(2, 3, 4, 5));
  ct.SetConfig(cfg);

  // get config using Z shape
  OBCisTransStereo::Config cfg2 = ct.GetConfig(OBStereo::ShapeZ);
  OB_ASSERT( cfg2.refs.size() == 4 );
  OB_ASSERT( cfg2.refs[0] == 2 );
  OB_ASSERT( cfg2.refs[1] == 3 );
  OB_ASSERT( cfg2.refs[2] == 5 );
  OB_ASSERT( cfg2.refs[3] == 4 );

  // get config using 4 shape
  cfg2 = ct.GetConfig(OBStereo::Shape4);
  OB_ASSERT( cfg2.refs.size() == 4 );
  OB_ASSERT( cfg2.refs[0] == 2 );
  OB_ASSERT( cfg2.refs[1] == 4 );
  OB_ASSERT( cfg2.refs[2] == 3 );
  OB_ASSERT( cfg2.refs[3] == 5 );
}
  
void test_IsOnSameAtom1()
{
  OBConversion conv;
  conv.SetInFormat("smi");
  OBMol mol;
  //
  // C       I     C      I    0      5
  //  \     /      |      |    |      |
  //   C===C       | C  C |    | 1  3 |
  //  /     \      |      |    |      |
  // Cl      F     Cl-----F    2------4
  //
  OB_ASSERT( conv.ReadString(&mol, "CC(Cl)=C(F)I") );
  OB_ASSERT( mol.NumAtoms() == 6 );
  
  OB_ASSERT( mol.GetAtomById(1) );
  OB_ASSERT( mol.GetAtomById(1)->IsCarbon() );
  OB_ASSERT( mol.GetAtomById(1)->GetValence() == 3);
  OB_ASSERT( mol.GetAtomById(3) );
  OB_ASSERT( mol.GetAtomById(3)->IsCarbon() );
  OB_ASSERT( mol.GetAtomById(3)->GetValence() == 3);

  OB_ASSERT( mol.GetAtomById(4) );
  OB_ASSERT( mol.GetAtomById(4)->GetAtomicNum() == 9 );
  OB_ASSERT( mol.GetAtomById(5) );
  OB_ASSERT( mol.GetAtomById(5)->GetAtomicNum() == 53 );

  OBCisTransStereo ct(&mol);
  // set refs using default U shape
  ct.SetConfig(OBCisTransStereo::Config(1, 3, OBStereo::MakeRefs(0, 2, 4, 5)));

  OB_ASSERT( ct.IsOnSameAtom(0, 2) );
  OB_ASSERT( ct.IsOnSameAtom(2, 0) );
  OB_ASSERT( ct.IsOnSameAtom(4, 5) );
  OB_ASSERT( ct.IsOnSameAtom(5, 4) );
  
  OB_ASSERT( !ct.IsOnSameAtom(0, 5) );
  OB_ASSERT( !ct.IsOnSameAtom(5, 0) );
  OB_ASSERT( !ct.IsOnSameAtom(0, 4) );
  OB_ASSERT( !ct.IsOnSameAtom(4, 0) );
 
  OB_ASSERT( !ct.IsOnSameAtom(2, 5) );
  OB_ASSERT( !ct.IsOnSameAtom(5, 2) );
  OB_ASSERT( !ct.IsOnSameAtom(2, 4) );
  OB_ASSERT( !ct.IsOnSameAtom(4, 2) );
  
}

void test_IsOnSameAtom2()
{
  obErrorLog.SetOutputLevel(obInfo);
  OBConversion conv;
  conv.SetInFormat("smi");
  OBMol mol;
  //
  //  Br      I     C      I    0      5
  //   \     /      |      |    |      |
  //    C===C       | C  C |    | 1  3 |
  //   /     \      |      |    |      |
  //  Cl      F     Cl-----F    2------4
  //
  OB_ASSERT( conv.ReadString(&mol, "[Br]C([Cl])=C(F)I") );
  OB_ASSERT( mol.NumAtoms() == 6 );
  
  OBCisTransStereo ct(&mol);
  // set refs using default U shape
  ct.SetConfig(OBCisTransStereo::Config(1, 3, OBStereo::MakeRefs(0, 2, 4, 5)));
  
  OB_ASSERT( ct.IsOnSameAtom(0, 2) );
  OB_ASSERT( ct.IsOnSameAtom(2, 0) );
  OB_ASSERT( ct.IsOnSameAtom(4, 5) );
  OB_ASSERT( ct.IsOnSameAtom(5, 4) );
  OB_ASSERT( !ct.IsOnSameAtom(0, 5) );
  OB_ASSERT( !ct.IsOnSameAtom(5, 0) );
  OB_ASSERT( !ct.IsOnSameAtom(0, 4) );
  OB_ASSERT( !ct.IsOnSameAtom(4, 0) );

  //
  // Trans
  //

  OB_ASSERT( ct.GetTransRef(0) == 4 );
  OB_ASSERT( ct.GetTransRef(4) == 0 );
  OB_ASSERT( ct.GetTransRef(2) == 5 );
  OB_ASSERT( ct.GetTransRef(5) == 2 );

  OB_ASSERT( ct.GetTransRef(0) != 2 );
  OB_ASSERT( ct.GetTransRef(2) != 0 );
  OB_ASSERT( ct.GetTransRef(4) != 5 );
  OB_ASSERT( ct.GetTransRef(5) != 4 );
  
  OB_ASSERT( ct.IsTrans(0, 4) );
  OB_ASSERT( ct.IsTrans(2, 5) );
  OB_ASSERT( !ct.IsTrans(0, 2) );
  OB_ASSERT( !ct.IsTrans(0, 5) );

  //
  // Cis
  //

  OB_ASSERT( ct.GetCisRef(0) == 5 );
  OB_ASSERT( ct.GetCisRef(5) == 0 );
  OB_ASSERT( ct.GetCisRef(2) == 4 );
  OB_ASSERT( ct.GetCisRef(4) == 2 );

  OB_ASSERT( ct.GetCisRef(2) != 5 );
  OB_ASSERT( ct.GetCisRef(5) != 2 );
  OB_ASSERT( ct.GetCisRef(0) != 4 );
  OB_ASSERT( ct.GetCisRef(4) != 0 );
  
  OB_ASSERT( ct.IsCis(0, 5) );
  OB_ASSERT( ct.IsCis(2, 4) );
  OB_ASSERT( !ct.IsCis(0, 4) );
  OB_ASSERT( !ct.IsCis(2, 5) );

}

void test_CisTrans1()
{
  obErrorLog.SetOutputLevel(obInfo);
  OBConversion conv;
  conv.SetInFormat("smi");
  OBMol mol;
  // 
  //  F       F      F      F     0      3
  //   \     /       |      |     |      |
  //    C===C        | C  C |     | 1  2 |
  //   /     \       |      |     |      |
  // (H)     (H)    (H)----(H)    H------H
  //
  OB_ASSERT( conv.ReadString(&mol, "FC=CF") );
  OBCisTransStereo ct(&mol);
  ct.SetConfig(OBCisTransStereo::Config(1, 2, OBStereo::MakeRefs(0, 
      OBStereo::ImplicitRef, OBStereo::ImplicitRef, 3)));

  OB_ASSERT( ct.IsOnSameAtom(0, OBStereo::ImplicitRef) );
  OB_ASSERT( ct.IsOnSameAtom(3, OBStereo::ImplicitRef) );
  OB_ASSERT( !ct.IsOnSameAtom(OBStereo::ImplicitRef, OBStereo::ImplicitRef) );
  OB_ASSERT( !ct.IsOnSameAtom(0, 3) );

  OB_ASSERT( ct.GetCisRef(0) == 3 );
  OB_ASSERT( ct.GetCisRef(3) == 0 );
  
  OB_ASSERT( ct.GetTransRef(0) == OBStereo::ImplicitRef );
  OB_ASSERT( ct.GetTransRef(3) == OBStereo::ImplicitRef );
}

void test_CisTrans2()
{
  obErrorLog.SetOutputLevel(obInfo);
  OBConversion conv;
  conv.SetInFormat("smi");
  OBMol mol;
  // 
  //  F      (H)     F     (H)   0      H
  //   \     /       |      |    |      |
  //    C===C        | C  C |    | 1  2 |
  //   /     \       |      |    |      |
  // (H)      F     (H)-----F    H------3
  //
  OB_ASSERT( conv.ReadString(&mol, "FC=CF") );
  OBCisTransStereo ct(&mol);
  OBCisTransStereo::Config cfg(1, 2, OBStereo::MakeRefs(0, OBStereo::ImplicitRef, 3, OBStereo::ImplicitRef));
  ct.SetConfig(cfg);

  OB_ASSERT( ct.IsOnSameAtom(0, OBStereo::ImplicitRef) );
  OB_ASSERT( ct.IsOnSameAtom(3, OBStereo::ImplicitRef) );
  OB_ASSERT( !ct.IsOnSameAtom(OBStereo::ImplicitRef, OBStereo::ImplicitRef) );
  OB_ASSERT( !ct.IsOnSameAtom(0, 3) );

  OB_ASSERT( ct.GetTransRef(0) == 3 );
  OB_ASSERT( ct.GetTransRef(3) == 0 );
  
  OB_ASSERT( ct.GetCisRef(0) == OBStereo::ImplicitRef );
  OB_ASSERT( ct.GetCisRef(3) == OBStereo::ImplicitRef );
 

}


int main() 
{
  // Define location of file formats for testing
  #ifdef FORMATDIR
    char env[BUFF_SIZE];
    snprintf(env, BUFF_SIZE, "BABEL_LIBDIR=%s", FORMATDIR);
    putenv(env);
  #endif

  test_GetType();
  test_configStruct();
  test_IsValid();
  test_equalsOperator();
  test_GetSetConfig();
  testRefs();
  test_IsOnSameAtom1();
  test_IsOnSameAtom2();
  test_CisTrans1();
  test_CisTrans2();
  
  return 0;
}

                

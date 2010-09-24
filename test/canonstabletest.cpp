#include "obtest.h"
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>

/*
#include <openbabel/graphsym.h>
#include <openbabel/stereo/tetrahedral.h>
#include <openbabel/canon.h>
*/

#include <iostream>
#include <vector>
#include <algorithm>


using std::cout;
using std::endl;
using namespace OpenBabel;

bool testCanSmiles(const std::string &smiles, const std::string &stable_cansmiles)
{
  cout << " Testing: " << smiles << endl;
  // read a smiles string
  OBMol mol;
  OBConversion canConv, smiConv;
  OB_REQUIRE( canConv.SetInFormat("smi") );
  OB_REQUIRE( canConv.SetOutFormat("can") );
  OB_REQUIRE( smiConv.SetOutFormat("smi") );
  // read a smiles string
  OB_REQUIRE( canConv.ReadString(&mol, smiles) );


  // get can smiles
  std::string cansmiles = canConv.WriteString(&mol, true);
  OB_COMPARE( cansmiles, stable_cansmiles );
  // comapare with ref
  if (cansmiles != stable_cansmiles) {
    cout << " " << cansmiles << endl;
    cout << " " << stable_cansmiles << endl;
    return false;
  }

  return true;
}


int main(int argc, char **argv)
{
  // Define location of file formats for testing
#ifdef FORMATDIR
    char env[BUFF_SIZE];
    snprintf(env, BUFF_SIZE, "BABEL_LIBDIR=%s", FORMATDIR);
    putenv(env);
#endif  

  testCanSmiles("N#CC1[C@@H](C#N)[C@@H]1C#N.N#CC1[C@@H](C#N)[C@@H]1C#N","C1([C@H]([C@@H]1C#N)C#N)C#N.C1([C@H]([C@@H]1C#N)C#N)C#N");
  testCanSmiles("O=C(N[C@@H]1[C@@H]2C[C@H]3C[C@@H](C2)C[C@@H]1C3)c1ccncn1","c1(C(=O)N[C@@H]2[C@@H]3C[C@@H]4C[C@H]2C[C@H](C3)C4)ncncc1");
  testCanSmiles("Clc1nc(Cl)nc(c1)O[C@H]1[C@@H]2C[C@H]3C[C@@H](C2)C[C@@H]1C3","[C@H]1([C@@H]2C[C@@H]3C[C@H]1C[C@H](C2)C3)Oc1nc(nc(c1)Cl)Cl");
  testCanSmiles("OC(=O)[C@@H]1[C@@H](c2ccccc2)[C@H](C(=O)O)[C@@H]1c1ccccc1","[C@H]1([C@H]([C@H]([C@H]1C(=O)O)c1ccccc1)C(=O)O)c1ccccc1");
  testCanSmiles("OC(=O)[C@]12[C@@H]3CC[C@@H](O3)[C@@]1(C(=O)O)[C@H]1CC[C@@H]2O1","[C@]12([C@]([C@H]3O[C@@H]2CC3)([C@H]2O[C@@H]1CC2)C(=O)O)C(=O)O");
  testCanSmiles("Nc1nc(S)nc(c1)C(=O)N[C@@H]1[C@@H]2C[C@H]3C[C@@H](C2)C[C@@H]1C3","c1(C(=O)N[C@@H]2[C@@H]3C[C@@H]4C[C@H]2C[C@H](C3)C4)nc(nc(c1)N)S");
  testCanSmiles("O1[C@@H]2[C@@H]1[C@@H]1O[C@H]1[C@@H]1O[C@H]1[C@@H]1O[C@H]21","[C@@H]12[C@H](O2)[C@H]2[C@@H](O2)[C@H]2[C@H]([C@H]3[C@H]1O3)O2");
  testCanSmiles("c1ccc(cc1)C[C@@]12N=N[C@@](Cc3ccccc3)([C@H]3C=C[C@@H]13)[C@H]1C=C[C@@H]21","[C@@]12([C@@H]3[C@H]([C@]([C@H]4[C@@H]1C=C4)(N=N2)Cc1ccccc1)C=C3)Cc1ccccc1");
  testCanSmiles("c1ccc2O[Cu@]34[N+](=Cc2c1)CCC[O+]3[Cu@]12Oc3ccccc3C=[N+]1CCC[O+]42","[Cu@@]123[O+]4[Cu@]5([O+]1CCC[N+]2=Cc1c(O3)cccc1)[N+](=Cc1c(O5)cccc1)CCC4");
  testCanSmiles("CN1C[C@H]2CN(C)C[C@@H](C1)[C@@]2(O)c1ccccc1","[C@@]1(c2ccccc2)(O)[C@@H]2CN(C[C@H]1CN(C2)C)C");
  testCanSmiles("C1CC[C@H]2N(C1)C1CCCCN1[C@@H]1CCCCN21","[C@@H]12N3[C@H](N4C(N1CCCC2)CCCC4)CCCC3");
  testCanSmiles("O=C1C[C@@H]2[C@@H]3C[C@@H](C3)[C@H](C1)N2C","[C@@H]12N([C@@H]([C@H]3C[C@@H]1C3)CC(=O)C2)C");
  testCanSmiles("[NaH].CCCCCCCCCCCCCCCCCCSC[C@@]12[C@H](COS(=O)(=O)O)[C@@H](COS(=O)(=O)O)[C@@H](c3ccccc13)c1ccccc21",
                "[C@]12(c3c([C@@H](c4c1cccc4)[C@@H]([C@H]2COS(=O)(=O)O)COS(=O)(=O)O)cccc3)CSCCCCCCCCCCCCCCCCCC.[NaH]");
  testCanSmiles("[NaH].Cl.CN1CCCCC1.O=C(Nc1ccc(cc1)S(=O)(=O)c1ccc(cc1)NC(=O)c1cc(ccc1O)S(=O)(=O)O)c1cc(ccc1O)S(=O)(=O)O",
                "S(=O)(=O)(c1cc(c(cc1)O)C(=O)Nc1ccc(cc1)S(=O)(=O)c1ccc(cc1)NC(=O)c1c(ccc(c1)S(=O)(=O)O)O)O.N1(CCCCC1)C.[NaH].Cl");
  testCanSmiles("CC(=O)O.N=C\\1/C=C/C(=C(/c2ccc(N)cc2)\\c2ccc(N)c(C)c2)/C=C1","c1(c(cc(cc1)C(=C1C=CC(=N)C=C1)c1ccc(cc1)N)C)N.C(=O)(O)C");
  testCanSmiles("CCNc1ccc(cc1)/C(=C/1\\C=C/C(=N)/C=C1)/c1ccc(N(CC)CC)c(C)c1","c1(c(cc(cc1)C(=C1C=CC(=N)C=C1)c1ccc(cc1)NCC)C)N(CC)CC");
  testCanSmiles("C[C@H]1CCCCN1/C=C/1\\C(=O)O[C@@](C)(C)OC1=O","C1(OC(=O)C(=CN2[C@H](CCCC2)C)C(=O)O1)(C)C");
  testCanSmiles("N=C\\1/C=C/C(=C(\\c2ccc(N)cc2)/c2ccc(N)c(C)c2)/C=C1","c1(c(cc(cc1)C(=C1C=CC(=N)C=C1)c1ccc(cc1)N)C)N");
  testCanSmiles("N=C\\1/C=C/C(=N\\Cc2ccc(N)cc2)/C=C1","C1(=NCc2ccc(cc2)N)C=CC(=N)C=C1");
  testCanSmiles("[O+]#C[Fe+]1234(C#[O+])(/C=C/CCC/C=C/[Fe+]5678(C#[O+])(C#[O+])[C@H]9C7=C6C5=C89)[C@H]5C3=C2C1=C45.F[B-](F)(F)F",
                "[Fe+]1234(C5=C3C2=C1C45)(/C=C/CCC/C=C/[Fe+]1234(C5=C3C2=C1C45)(C#[O+])C#[O+])(C#[O+])C#[O+].[B-](F)(F)(F)F");
  testCanSmiles("[O+]#C[Fe+]1234(C#[O+])(/C=C/CCCC/C=C/[Fe+]5678(C#[O+])(C#[O+])C9=C6C7=C5[C@@H]89)C5=C3C2=C1[C@H]45.F[B-](F)(F)F",
                "[Fe+]1234(C5=C3C2=C1C45)(/C=C/CCCC/C=C/[Fe+]1234(C5=C3C2=C1C45)(C#[O+])C#[O+])(C#[O+])C#[O+].[B-](F)(F)(F)F");
  testCanSmiles("c1[14cH]cccc1","[14cH]1ccccc1");
  testCanSmiles("[14cH]1[14cH]cccc1","[14cH]1[14cH]cccc1");
  testCanSmiles("[14cH]1[14cH]ccc[14cH]1","[14cH]1[14cH]ccc[14cH]1");
  testCanSmiles("[14cH]1[14cH]cc[14cH][14cH]1","[14cH]1[14cH][14cH]cc[14cH]1");
  testCanSmiles("[14cH]1[14cH]c[14cH][14cH][14cH]1","[14cH]1[14cH][14cH]c[14cH][14cH]1");
  testCanSmiles("[14cH]1[14cH][14cH][14cH][14cH][14cH]1","[14cH]1[14cH][14cH][14cH][14cH][14cH]1");
  testCanSmiles("F[Po@SP1](Cl)(Br)I","[Po@SP1](I)(Br)(Cl)F");
  testCanSmiles("F[Po@SP2](Br)(Cl)I","[Po@SP1](I)(Br)(Cl)F");
  testCanSmiles("F[Po@SP3](Cl)(I)Br","[Po@SP1](I)(Br)(Cl)F");
  testCanSmiles("Br.CC(=O)Nc1ccc(cc1)S(=O)(=O)c1ccc(cc1)NC(=O)c1ccccc1SC(=O)CCCCn1ccccc1","S(=O)(=O)(c1ccc(cc1)NC(=O)C)c1ccc(cc1)NC(=O)c1c(SC(=O)CCCCN2CCCCC2)cccc1.Br");
  testCanSmiles("Br.CC(=O)Nc1cccc(c1)S(=O)(=O)c1cccc(NC(=O)c2ccccc2SC(=O)CCCCn2ccccc2)c1","S(=O)(=O)(c1cc(NC(=O)C)ccc1)c1cc(NC(=O)c2c(SC(=O)CCCCN3CCCCC3)cccc2)ccc1.Br");
  testCanSmiles("Br.COc1ccc(cc1)Cn1cccc(O)c1","c1(ccc(cc1)CN1CC(CCC1)O)OC.Br");
  testCanSmiles("Br.COc1ccc2c(c1)[C@]13CCCC[C@@H]3[C@H](N(C)CC1)[C@@]12CC1","C12(c3c([C@@]45[C@@H]([C@@H]1N(CC4)C)CCCC5)cc(cc3)OC)CC2.Br");
  testCanSmiles("Br.NCCCCNCCCC[C@-]12C3=C4C5=C1[Fe+2]16782345[C-]2C7=C6C1=C82","[Fe+2]12345678([C-]9(C1=C5C6=C29)CCCCNCCCCN)[C-]1C3=C7C8=C41.Br");
  testCanSmiles("Br.NCCCCNCCCC[C@@-]12C3=C4C5=C1[Fe+2]16782345[C-]2C7=C6C1=C82","[Fe+2]12345678([C-]9(C1=C5C6=C29)CCCCNCCCCN)[C-]1C3=C7C8=C41.Br");
  testCanSmiles("Br.NCCCNCCCCNCCCC[C@-]12C3=C4C5=C1[Fe+2]16782345[C-]2C7=C6C1=C82","[C-]12([Fe+2]3456789(C1=C4C5=C23)[C-]1C6=C8C9=C71)CCCCNCCCCNCCCN.Br");
  testCanSmiles("Br.NCCCNCCCCNCCCC[C@@-]12C3=C4C5=C1[Fe+2]16782345[C-]2C7=C6C1=C82","[C-]12([Fe+2]3456789(C1=C4C5=C23)[C-]1C6=C8C9=C71)CCCCNCCCCNCCCN.Br");
  testCanSmiles("Br.NCCCNCCCCNCCCNCCCC[C@-]12C3=C4C5=C1[Fe+2]16782345[C-]2C7=C6C1=C82","[C-]12([Fe+2]3456789(C1=C4C5=C23)[C-]1C6=C8C9=C71)CCCCNCCCNCCCCNCCCN.Br");
  testCanSmiles("Br.NCCCNCCCCNCCCNCCCC[C@@-]12C3=C4C5=C1[Fe+2]16782345[C-]2C7=C6C1=C82","[C-]12([Fe+2]3456789(C1=C4C5=C23)[C-]1C6=C8C9=C71)CCCCNCCCNCCCCNCCCN.Br");
  testCanSmiles("Br.NCCCNCCCC[C@-]12C3=C4C5=C1[Fe+2]16782345[C-]2C7=C6C1=C82","[Fe+2]12345678([C-]9(C1=C5C6=C29)CCCCNCCCN)[C-]1C3=C7C8=C41.Br");
  testCanSmiles("Br.NCCCNCCCC[C@@-]12C3=C4C5=C1[Fe+2]16782345[C-]2C7=C6C1=C82","[Fe+2]12345678([C-]9(C1=C5C6=C29)CCCCNCCCN)[C-]1C3=C7C8=C41.Br");
  testCanSmiles("CC(=O)O.Oc1nc(N)c2nc3c(nc2n1)c1cccc2cccc3c12.OS(=O)(=O)O","S(=O)(=O)(O)O.c12c3c(c4c1[nH]c1c(n4)nc([nH]c1N)O)cccc3ccc2.C(=O)(O)C");
  testCanSmiles("C12=C3C4=C5[C@H]1[Fe]16782345C2=C6C7=C1C82","[Fe]12345678(C9C1=C5C6=C29)C1C3=C7C8=C41");
  testCanSmiles("C12([C@@H]3CCC[C@@H]1CCC[C@@H]2CCC3)OC(=O)C","C12([C@@H]3CCC[C@@H]2CCC[C@H]1CCC3)OC(=O)C");
  testCanSmiles("[Fe]12345678(C9(C1=C4C3=C29)C=O)C1=C6C7=C5C81","[Fe]12345678(C9(C1=C5C6=C29)C=O)C1C3=C7C8=C41");
  testCanSmiles("[Fe+2]12345678([C-]9(C1=C7C6=C29)CCCCNCCCCN)C1=C4C5=C3[C-]81.Br","[Fe+2]12345678([C-]9(C1=C5C6=C29)CCCCNCCCCN)[C-]1C3=C7C8=C41.Br");
  testCanSmiles("O[C@H]1CC[C@@H](O)CC1","[C@H]1(CC[C@@H](CC1)O)O");
  testCanSmiles("c1csc(c1)c1nnc(o1)C12C[C@@H]3C[C@@H](C[C@@H](C3)C2)C1","C12(c3oc(nn3)c3sccc3)C[C@@H]3C[C@H](C2)C[C@H](C1)C3");
  testCanSmiles("CC(=O)OC12[C@@H]3CCC[C@@H]1CCC[C@@H]2CCC3","C12([C@@H]3CCC[C@@H]2CCC[C@H]1CCC3)OC(=O)C");
  testCanSmiles("CC(=O)OC12[C@@H]3CCC[C@@H]1CCC[C@@H]2CCC3","C12([C@@H]3CCC[C@@H]2CCC[C@H]1CCC3)OC(=O)C");
  testCanSmiles("C/C(=N/c1ccncc1)/C(=N/c1ccncc1)/C","C(=N\\c1ccncc1)(\\C(=N\\c1ccncc1)\\C)/C");
  testCanSmiles("CC(=CCN[C@]12C[C@@H]3C[C@@H](C[C@@H](C3)C2)C1)C","[C@]12(C[C@@H]3C[C@H](C2)C[C@H](C1)C3)NCC=C(C)C");
  return 0;
}


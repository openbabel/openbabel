/**********************************************************************
residue.cpp - Unit tests for Open Babel OBResidue class

Copyright (C) 2005-2006 Geoffrey R. Hutchison

This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

// used to set import/export for Cygwin DLLs
#ifdef WIN32
#define USING_OBDLL
#endif

#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/obiter.h>
#include <openbabel/chains.h>
#include <openbabel/residue.h>
#include <cstdlib>
#include <stdio.h>
#include <iostream>

using namespace std;
using namespace OpenBabel;

// helper functions for chains parsing
void CheckValidResidue(OBConversion &conv,
                         const string &test,
                         unsigned int testCount);
void CheckInvalidResidue(OBConversion &conv,
                         const string &test,
                         unsigned int testCount);
void CheckValidDipeptide(OBConversion &conv,
                         const string &test,
                         unsigned int testCount);

int residue(int argc, char* argv[])
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


  cout << "# Unit tests for OBResidue \n";

  cout << "ok 1\n"; // for loading tests

  // OBResidue isolation tests
  OBResidue emptyResidue, testRes1;
  cout << "ok 2\n"; // ctor works

  // chains parser tests

  // PR#1515198
  static const string loopTest1("C1(C(NC(C(N1C(C(NC(C=Cc1ccccc1)=O)C)=O)Cc1ccccc1)=O)Cc1ccccc1)=O");
  OBConversion conv;
  OBMol mol;
  OBFormat *inFormat = conv.FindFormat("SMI");

  conv.SetInFormat(inFormat);
  conv.ReadString(&mol, loopTest1);
  chainsparser.PerceiveChains(mol);

  // OK if it doesn't crash
  cout << "ok 3\n";

  // parse common residues
  unsigned int testCount = 3;
  static const string ala("NC(C)C(O)(=O)");
  CheckValidResidue(conv, ala, ++testCount);
  static const string arg("NC(CCCNC(N)=N)C(O)(=O)");
  CheckValidResidue(conv, arg, ++testCount);
  static const string asn("NC(CC(N)=O)C(O)(=O)");
  CheckValidResidue(conv, asn, ++testCount);
  static const string asp("NC(CC(O)=O)C(O)(=O)");
  CheckValidResidue(conv, asp, ++testCount);
  static const string cys("NC(CS)C(O)(=O)");
  CheckValidResidue(conv, cys, ++testCount);
  static const string glu("NC(CCC(O)=O)C(O)(=O)");
  CheckValidResidue(conv, glu, ++testCount);
  static const string gln("NC(CCC(N)=O)C(O)(=O)");
  CheckValidResidue(conv, gln, ++testCount);
  static const string gly("NC([H])C(O)(=O)");
  CheckValidResidue(conv, gly, ++testCount);
  static const string his("NC(CC1=CNC=N1)C(O)(=O)");
  CheckValidResidue(conv, his, ++testCount);
  static const string ile("NC(C(CC)C)C(O)(=O)");
  CheckValidResidue(conv, ile, ++testCount);
  static const string leu("NC(CC(C)C)C(O)(=O)");
  CheckValidResidue(conv, leu, ++testCount);
  static const string lys("NC(CCCCN)C(O)(=O)");
  CheckValidResidue(conv, lys, ++testCount);
  static const string met("NC(CCSC)C(O)(=O)");
  CheckValidResidue(conv, met, ++testCount);
  static const string phe("NC(CC1=CC=CC=C1)C(O)(=O)");
  CheckValidResidue(conv, phe, ++testCount);
  static const string pro("OC(C1CCCN1)(=O)");
  CheckValidResidue(conv, pro, ++testCount);
  static const string ser("NC(CO)C(O)(=O)");
  CheckValidResidue(conv, ser, ++testCount);
  static const string thr("NC(C(C)O)C(O)(=O)");
  CheckValidResidue(conv, thr, ++testCount);
  static const string trp("NC(CC1=CNC2=C1C=CC=C2)C(O)(=O)");
  CheckValidResidue(conv, trp, ++testCount);
  static const string tyr("NC(CC1=CC=C(O)C=C1)C(O)(=O)");
  CheckValidResidue(conv, tyr, ++testCount);
  static const string val("NC(C(C)C)C(O)(=O)");
  CheckValidResidue(conv, val, ++testCount);

  // nucleics
  static const string a("OC[C@H]1O[C@H](C[C@@H]1O)n1cnc2c(ncnc12)N");
  CheckValidResidue(conv, a, ++testCount);
  static const string g("OC[C@H]1O[C@H](C[C@@H]1O)n1c(nc(cc1)N)=O");
  CheckValidResidue(conv, g, ++testCount);
  static const string c("OC[C@H]1O[C@H](C[C@@H]1O)n1cnc2c([nH]c(nc12)N)=O");
  CheckValidResidue(conv, c, ++testCount);
  static const string t("OC[C@H]1O[C@H](C[C@@H]1O)n1c([nH]c(c(c1)C)=O)=O");
  CheckValidResidue(conv, t, ++testCount);
  static const string u("OC[C@H]1O[C@H]([C@@H]([C@@H]1O)O)n1c([nH]c(cc1)=O)=O");
  CheckValidResidue(conv, u, ++testCount);

  // invalid residues
  static const string benzene("c1ccccc1");
  CheckInvalidResidue(conv, benzene, ++testCount);
  static const string pyrrole("c1cccn[H]1");
  CheckInvalidResidue(conv, pyrrole, ++testCount);
  static const string amine("CC(=O)CCN");
  CheckInvalidResidue(conv, amine, ++testCount);

  // check some dipeptides
  static const string ala_val("NC(C)C(=O)NC(C(C)C)C(=O)O");
  CheckValidDipeptide(conv, ala_val, ++testCount);
  static const string cys_leu("NC(CS)C(=O)NC(CC(C)C)C(=O)O");
  CheckValidDipeptide(conv, cys_leu, ++testCount);

  // the number of tests for "prove"
  cout << "1.." << testCount << "\n";

  return(0);
}

void CheckValidResidue(OBConversion &conv,
                       const string &test,
                       unsigned int testCount)
{
  OBMol mol;

  mol.Clear();
  conv.ReadString(&mol, test);
  chainsparser.PerceiveChains(mol);
  if (mol.NumResidues() != 1)
    {
      cout << "not ok " << testCount << " # expected 1 residue, but found "
           << mol.NumResidues() << '\n';
      cout << "# ";
      FOR_RESIDUES_OF_MOL(res, mol)
        cout << res->GetName() << " ";
      cout << endl;
    }
  else {
    OBResidue *res;
    res = mol.GetResidue(0);
    cout << "ok " << testCount << " # " << res->GetName() << endl;
  }
}


void CheckInvalidResidue(OBConversion &conv,
                         const string &test,
                         unsigned int testCount)
{
  OBMol mol;

  mol.Clear();
  conv.ReadString(&mol, test);
  chainsparser.PerceiveChains(mol);
  if (mol.NumResidues() != 0) {
    OBResidue *res = mol.GetResidue(0);
    if (res->GetName() == "LIG" || res->GetName() == "UNL" ) { // ligand, not residue
      cout << "ok " << testCount << " # found ligand, not residue "
           << test << '\n';
    } else {
      cout << "not ok " << testCount << " # expected 0 residues, found "
           << mol.NumResidues() << '\n';
      cout << "# " << res->GetName() << endl;
    }
  } else
    cout << "ok " << testCount << " # correctly rejected " << test << '\n';
}

void CheckValidDipeptide(OBConversion &conv,
                         const string &test,
                         unsigned int testCount)
{
  OBMol mol;

  mol.Clear();
  conv.ReadString(&mol, test);
  chainsparser.PerceiveChains(mol);
  if (mol.NumResidues() != 2) {
    cout << "not ok " << testCount << " # expected 2 residues, but found "
         << mol.NumResidues() << '\n';
    cout << "# ";
    FOR_RESIDUES_OF_MOL(res, mol)
      cout << res->GetName() << " ";
    cout << endl;
  } else {
    OBResidue *res;
    res = mol.GetResidue(0);
    cout << "ok " << testCount << " # " << res->GetName();
    res = mol.GetResidue(1);
    cout << " " << res->GetName() << '\n';
  }
}

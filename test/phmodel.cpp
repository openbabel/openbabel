/**********************************************************************
phmodel.cpp - Unit tests for Open Babel OBPhModel class

Copyright (C) 2008 Tim Vandermmeersch
 
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

#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/forcefield.h>
#include <openbabel/phmodel.h>
#include <cstdlib>

using namespace std;
using namespace OpenBabel;

int phmodel(int argc, char* argv[])
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

  OBMol mol;
  OBConversion conv;
  conv.SetInFormat("smi");
  conv.SetOutFormat("smi");

  unsigned int test = 0;

  // amine acid COOH pKa = 4.0 (carboxylic acid entry in phmodel.txt)
  // amino acid NH3+ pKa = 10.0 (amine entry in phmodel.txt)

  // 
  // Aspartic acid (sidechain COOH pKa = 3.8)
  //
  cout << "# Asp" << endl;
  conv.ReadString(&mol, "NC(CC(O)=O)C(O)=O");
  mol.SetAutomaticFormalCharge(true);
  mol.AddHydrogens(false, true, 1.0); // NH3+ COOH COOH
  cout << "#pH = 1.0 : ";
  conv.Write(&mol, &cout);
  
  if (mol.NumAtoms() == 17)
    cout << "ok " << ++test << "\n";
  else
    cout << "not ok " << ++test << "\n";

  conv.ReadString(&mol, "NC(CC(O)=O)C(O)=O");
  mol.SetAutomaticFormalCharge(true);
  mol.AddHydrogens(false, true, 3.9); // NH3+ COOH COO-
  cout << "#pH = 3.9 : ";
  conv.Write(&mol, &cout);
 
  if (mol.NumAtoms() == 16)
    cout << "ok " << ++test << "\n";
  else
    cout << "not ok " << ++test << "\n";

  conv.ReadString(&mol, "NC(CC(O)=O)C(O)=O");
  mol.SetAutomaticFormalCharge(true);
  mol.AddHydrogens(false, true, 7.4); // NH3+ COO- COO-
  cout << "#pH = 7.4 : ";
  conv.Write(&mol, &cout);
 
  if (mol.NumAtoms() == 15)
    cout << "ok " << ++test << "\n";
  else
    cout << "not ok " << ++test << "\n";

  conv.ReadString(&mol, "NC(CC(O)=O)C(O)=O");
  mol.SetAutomaticFormalCharge(true);
  mol.AddHydrogens(false, true, 13.0); // NH2 COO- COO-
  cout << "#pH = 13.0 : ";
  conv.Write(&mol, &cout);

  if (mol.NumAtoms() == 14)
    cout << "ok " << ++test << "\n";
  else
    cout << "not ok " << ++test << "\n";

  // 
  // Glutamic acid (sidechain COOH pKa = 4.3)
  //
  cout << "# Glu" << endl;
  conv.ReadString(&mol, "NC(CCC(O)=O)C(O)=O");
  mol.SetAutomaticFormalCharge(true);
  mol.AddHydrogens(false, true, 1.0); // NH3+ COOH COOH
  cout << "#pH = 1.0 : ";
  conv.Write(&mol, &cout);
  
  if (mol.NumAtoms() == 20)
    cout << "ok " << ++test << "\n";
  else
    cout << "not ok " << ++test << "\n";

  conv.ReadString(&mol, "NC(CCC(O)=O)C(O)=O");
  mol.SetAutomaticFormalCharge(true);
  mol.AddHydrogens(false, true, 4.15); // NH3+ COOH COO-
  cout << "#pH = 4.15 : ";
  conv.Write(&mol, &cout);
 
  // known bug
//   if (mol.NumAtoms() == 19)
//     cout << "ok " << ++test << "\n";
//   else
//     cout << "not ok " << ++test << "\n";

  conv.ReadString(&mol, "NC(CCC(O)=O)C(O)=O");
  mol.SetAutomaticFormalCharge(true);
  mol.AddHydrogens(false, true, 7.4); // NH3+ COO- COO-
  cout << "#pH = 7.4 : ";
  conv.Write(&mol, &cout);
 
  if (mol.NumAtoms() == 18)
    cout << "ok " << ++test << "\n";
  else
    cout << "not ok " << ++test << "\n";

  conv.ReadString(&mol, "NC(CCC(O)=O)C(O)=O");
  mol.SetAutomaticFormalCharge(true);
  mol.AddHydrogens(false, true, 13.0); // NH2 COO- COO-
  cout << "#pH = 13.0 : ";
  conv.Write(&mol, &cout);

  if (mol.NumAtoms() == 17)
    cout << "ok " << ++test << "\n";
  else
    cout << "not ok " << ++test << "\n";

  // 
  // Histidine (sidechain nH+ pKa = 6.08)
  //
  cout << "# His" << endl;
  conv.ReadString(&mol, "NC(Cc1nc[nH]c1)C(O)=O");
  mol.SetAutomaticFormalCharge(true);
  mol.AddHydrogens(false, true, 1.0); // NH3+ COOH nH+
  cout << "#pH = 1.0 : ";
  conv.Write(&mol, &cout);
  
  if (mol.NumAtoms() == 22)
    cout << "ok " << ++test << "\n";
  else
    cout << "not ok " << ++test << "\n";
  
  conv.ReadString(&mol, "NC(Cc1nc[nH]c1)C(O)=O");
  mol.SetAutomaticFormalCharge(true);
  mol.AddHydrogens(false, true, 5.0); // NH3+ COO- nH+
  cout << "#pH = 5.0 : ";
  conv.Write(&mol, &cout);
 
  if (mol.NumAtoms() == 21)
    cout << "ok " << ++test << "\n";
  else
    cout << "not ok " << ++test << "\n";

  conv.ReadString(&mol, "NC(Cc1nc[nH]c1)C(O)=O");
  mol.SetAutomaticFormalCharge(true);
  mol.AddHydrogens(false, true, 7.4); // NH3+ COO- n:
  cout << "#pH = 7.4 : ";
  conv.Write(&mol, &cout);
 
  if (mol.NumAtoms() == 20)
    cout << "ok " << ++test << "\n";
  else
    cout << "not ok " << ++test << "\n";

  conv.ReadString(&mol, "NC(Cc1nc[nH]c1)C(O)=O");
  mol.SetAutomaticFormalCharge(true);
  mol.AddHydrogens(false, true, 13.0); // NH2 COO- n:
  cout << "#pH = 13.0 : ";
  conv.Write(&mol, &cout);

  if (mol.NumAtoms() == 19)
    cout << "ok " << ++test << "\n";
  else
    cout << "not ok " << ++test << "\n";

  // 
  // Lysine (sidechain NH3+ pKa = 8.28)
  //
  cout << "# Lys" << endl;
  conv.ReadString(&mol, "NC(CCCCN)C(O)=O");
  mol.SetAutomaticFormalCharge(true);
  mol.AddHydrogens(false, true, 1.0); // NH3+ COOH NH3+
  cout << "#pH = 1.0 : ";
  conv.Write(&mol, &cout);
  
  if (mol.NumAtoms() == 26)
    cout << "ok " << ++test << "\n";
  else
    cout << "not ok " << ++test << "\n";

  conv.ReadString(&mol, "NC(CCCCN)C(O)=O");
  mol.SetAutomaticFormalCharge(true);
  mol.AddHydrogens(false, true, 7.4); // NH3+ COO- NH3+
  cout << "#pH = 7.4 : ";
  conv.Write(&mol, &cout);
 
  if (mol.NumAtoms() == 25)
    cout << "ok " << ++test << "\n";
  else
    cout << "not ok " << ++test << "\n";

  conv.ReadString(&mol, "NC(CCCCN)C(O)=O");
  mol.SetAutomaticFormalCharge(true);
  mol.AddHydrogens(false, true, 9.0); // NH3+ COO- NH2
  cout << "#pH = 9.0 : ";
  conv.Write(&mol, &cout);
  
  // known bug
  //  if (mol.NumAtoms() == 24)
  //    cout << "ok " << ++test << "\n";
  //  else
  //    cout << "not ok " << ++test << "\n";

  conv.ReadString(&mol, "NC(CCCCN)C(O)=O");
  mol.SetAutomaticFormalCharge(true);
  mol.AddHydrogens(false, true, 13.0); // NH2 COO- NH2
  cout << "#pH = 13.0 : ";
  conv.Write(&mol, &cout);

  if (mol.NumAtoms() == 23)
    cout << "ok " << ++test << "\n";
  else
    cout << "not ok " << ++test << "\n";

  // 
  // Tyrosine (sidechain OH pKa = 10.1)
  //
  cout << "# Tyr" << endl;
  conv.ReadString(&mol, "NC(Cc1ccc(O)cc1)C(O)=O");
  mol.SetAutomaticFormalCharge(true);
  mol.AddHydrogens(false, true, 1.0); // NH3+ COOH NH3+
  cout << "#pH = 1.0 : ";
  conv.Write(&mol, &cout);
  
  if (mol.NumAtoms() == 25)
    cout << "ok " << ++test << "\n";
  else
    cout << "not ok " << ++test << "\n";

  conv.ReadString(&mol, "NC(Cc1ccc(O)cc1)C(O)=O");
  mol.SetAutomaticFormalCharge(true);
  mol.AddHydrogens(false, true, 7.4); // NH3+ COO- NH3+
  cout << "#pH = 7.4 : ";
  conv.Write(&mol, &cout);
 
  if (mol.NumAtoms() == 24)
    cout << "ok " << ++test << "\n";
  else
    cout << "not ok " << ++test << "\n";

  conv.ReadString(&mol, "NC(Cc1ccc(O)cc1)C(O)=O");
  mol.SetAutomaticFormalCharge(true);
  mol.AddHydrogens(false, true, 10.05); // NH3+ COO- NH2
  cout << "#pH = 10.05 : ";
  conv.Write(&mol, &cout);
 
  if (mol.NumAtoms() == 23)
    cout << "ok " << ++test << "\n";
  else
    cout << "not ok " << ++test << "\n";

  conv.ReadString(&mol, "NC(Cc1ccc(O)cc1)C(O)=O");
  mol.SetAutomaticFormalCharge(true);
  mol.AddHydrogens(false, true, 13.0); // NH2 COO- NH2
  cout << "#pH = 13.0 : ";
  conv.Write(&mol, &cout);

  if (mol.NumAtoms() == 22)
    cout << "ok " << ++test << "\n";
  else
    cout << "not ok " << ++test << "\n";

  // 
  // Arginine (sidechain =NH2+ pKa = 12.0)
  //
  cout << "# Arg" << endl;
  conv.ReadString(&mol, "NC(CCCNC(N)=N)C(O)=O");
  mol.SetAutomaticFormalCharge(true);
  mol.AddHydrogens(false, true, 1.0); // NH3+ COOH =NH2+
  cout << "#pH = 1.0 : ";
  conv.Write(&mol, &cout);
  
  if (mol.NumAtoms() == 28)
    cout << "ok " << ++test << "\n";
  else
    cout << "not ok " << ++test << "\n";

  conv.ReadString(&mol, "NC(CCCNC(N)=N)C(O)=O");
  mol.SetAutomaticFormalCharge(true);
  mol.AddHydrogens(false, true, 7.4); // NH3+ COO- NH3+
  cout << "#pH = 7.4 : ";
  conv.Write(&mol, &cout);
 
  if (mol.NumAtoms() == 27)
    cout << "ok " << ++test << "\n";
  else
    cout << "not ok " << ++test << "\n";

  conv.ReadString(&mol, "NC(CCCNC(N)=N)C(O)=O");
  mol.SetAutomaticFormalCharge(true);
  mol.AddHydrogens(false, true, 11.0); // NH3+ COO- NH2
  cout << "#pH = 11.0 : ";
  conv.Write(&mol, &cout);
 
  if (mol.NumAtoms() == 26)
    cout << "ok " << ++test << "\n";
  else
    cout << "not ok " << ++test << "\n";

  conv.ReadString(&mol, "NC(CCCNC(N)=N)C(O)=O");
  mol.SetAutomaticFormalCharge(true);
  mol.AddHydrogens(false, true, 13.0); // NH2 COO- NH2
  cout << "#pH = 13.0 : ";
  conv.Write(&mol, &cout);

  if (mol.NumAtoms() == 25)
    cout << "ok " << ++test << "\n";
  else
    cout << "not ok " << ++test << "\n";

  // 
  // Gly-Gly (terminal NH3+, COOH and the amide bond)
  //
  cout << "# Gly-Gly" << endl;
  conv.ReadString(&mol, "NCC(=O)NCC(=O)O");
  mol.SetAutomaticFormalCharge(true);
  mol.AddHydrogens(false, true, 1.0); // NH3+ COOH =NH2+
  cout << "#pH = 1.0 : ";
  conv.Write(&mol, &cout);
  
  if (mol.NumAtoms() == 18)
    cout << "ok " << ++test << "\n";
  else
    cout << "not ok " << ++test << "\n";

  conv.ReadString(&mol, "NCC(=O)NCC(=O)O");
  mol.SetAutomaticFormalCharge(true);
  mol.AddHydrogens(false, true, 7.4); // NH3+ COO- NH3+
  cout << "#pH = 7.4 : ";
  conv.Write(&mol, &cout);
 
  if (mol.NumAtoms() == 17)
    cout << "ok " << ++test << "\n";
  else
    cout << "not ok " << ++test << "\n";

  conv.ReadString(&mol, "NCC(=O)NCC(=O)O");
  mol.SetAutomaticFormalCharge(true);
  mol.AddHydrogens(false, true, 13.0); // NH2 COO- NH2
  cout << "#pH = 13.0 : ";
  conv.Write(&mol, &cout);

  if (mol.NumAtoms() == 16)
    cout << "ok " << ++test << "\n";
  else
    cout << "not ok " << ++test << "\n";

  // the number of tests for "prove"
  cout << "1.." << test << "\n";

//cout << mol.NumAtoms() << endl;  
  return 0;
}

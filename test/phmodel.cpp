/**********************************************************************
phmodel.cpp - Unit tests for Open Babel OBPhModel class

Copyright (C) 2008 Tim Vandermmeersch
 
This file is part of the Open Babel project.
For more information, see <http://openbabel.sourceforge.net/>

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

using namespace std;
using namespace OpenBabel;

int main (int argc, char **argv)
{
  OBMol mol;
  OBConversion conv;
  conv.SetInFormat("smi");
  conv.SetOutFormat("smi");

  // the number of tests for "prove"
  cout << "1..27\n";

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
    cout << "ok 1\n";
  else
    cout << "not ok 1\n";

  conv.ReadString(&mol, "NC(CC(O)=O)C(O)=O");
  mol.SetAutomaticFormalCharge(true);
  mol.AddHydrogens(false, true, 3.9); // NH3+ COOH COO-
  cout << "#pH = 3.9 : ";
  conv.Write(&mol, &cout);
 
  if (mol.NumAtoms() == 16)
    cout << "ok 2\n";
   else
    cout << "not ok 2\n";

  conv.ReadString(&mol, "NC(CC(O)=O)C(O)=O");
  mol.SetAutomaticFormalCharge(true);
  mol.AddHydrogens(false, true, 7.4); // NH3+ COO- COO-
  cout << "#pH = 7.4 : ";
  conv.Write(&mol, &cout);
 
  if (mol.NumAtoms() == 15)
    cout << "ok 3\n";
  else
    cout << "not ok 3\n";

  conv.ReadString(&mol, "NC(CC(O)=O)C(O)=O");
  mol.SetAutomaticFormalCharge(true);
  mol.AddHydrogens(false, true, 13.0); // NH2 COO- COO-
  cout << "#pH = 13.0 : ";
  conv.Write(&mol, &cout);

  if (mol.NumAtoms() == 14)
    cout << "ok 4\n";
  else
    cout << "not ok 4\n";

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
    cout << "ok 5\n";
  else
    cout << "not ok 5\n";

  conv.ReadString(&mol, "NC(CCC(O)=O)C(O)=O");
  mol.SetAutomaticFormalCharge(true);
  mol.AddHydrogens(false, true, 4.15); // NH3+ COOH COO-
  cout << "#pH = 4.15 : ";
  conv.Write(&mol, &cout);
 
  if (mol.NumAtoms() == 19)
    cout << "ok 6\n";
   else
    cout << "not ok 6\n";

  conv.ReadString(&mol, "NC(CCC(O)=O)C(O)=O");
  mol.SetAutomaticFormalCharge(true);
  mol.AddHydrogens(false, true, 7.4); // NH3+ COO- COO-
  cout << "#pH = 7.4 : ";
  conv.Write(&mol, &cout);
 
  if (mol.NumAtoms() == 18)
    cout << "ok 7\n";
  else
    cout << "not ok 7\n";

  conv.ReadString(&mol, "NC(CCC(O)=O)C(O)=O");
  mol.SetAutomaticFormalCharge(true);
  mol.AddHydrogens(false, true, 13.0); // NH2 COO- COO-
  cout << "#pH = 13.0 : ";
  conv.Write(&mol, &cout);

  if (mol.NumAtoms() == 17)
    cout << "ok 8\n";
  else
    cout << "not ok 8\n";

  // 
  // Histidine (sidechain nH+ pKa = 6.08)
  //
  cout << "# His" << endl;
  conv.ReadString(&mol, "NC(Cc1ncnc1)C(O)=O");
  mol.SetAutomaticFormalCharge(true);
  mol.AddHydrogens(false, true, 1.0); // NH3+ COOH nH+
  cout << "#pH = 1.0 : ";
  conv.Write(&mol, &cout);
  
  if (mol.NumAtoms() == 22)
    cout << "ok 9\n";
  else
    cout << "not ok 9\n";
  
  conv.ReadString(&mol, "NC(Cc1ncnc1)C(O)=O");
  mol.SetAutomaticFormalCharge(true);
  mol.AddHydrogens(false, true, 5.0); // NH3+ COO- nH+
  cout << "#pH = 5.0 : ";
  conv.Write(&mol, &cout);
 
  if (mol.NumAtoms() == 21)
    cout << "ok 10\n";
   else
    cout << "not ok 10\n";

  conv.ReadString(&mol, "NC(Cc1ncnc1)C(O)=O");
  mol.SetAutomaticFormalCharge(true);
  mol.AddHydrogens(false, true, 7.4); // NH3+ COO- n:
  cout << "#pH = 7.4 : ";
  conv.Write(&mol, &cout);
 
  if (mol.NumAtoms() == 20)
    cout << "ok 11\n";
  else
    cout << "not ok 11\n";

  conv.ReadString(&mol, "NC(Cc1ncnc1)C(O)=O");
  mol.SetAutomaticFormalCharge(true);
  mol.AddHydrogens(false, true, 13.0); // NH2 COO- n:
  cout << "#pH = 13.0 : ";
  conv.Write(&mol, &cout);

  if (mol.NumAtoms() == 19)
    cout << "ok 12\n";
  else
    cout << "not ok 12\n";

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
    cout << "ok 13\n";
  else
    cout << "not ok 13\n";

  conv.ReadString(&mol, "NC(CCCCN)C(O)=O");
  mol.SetAutomaticFormalCharge(true);
  mol.AddHydrogens(false, true, 7.4); // NH3+ COO- NH3+
  cout << "#pH = 7.4 : ";
  conv.Write(&mol, &cout);
 
  if (mol.NumAtoms() == 25)
    cout << "ok 14\n";
   else
    cout << "not ok 14\n";

  conv.ReadString(&mol, "NC(CCCCN)C(O)=O");
  mol.SetAutomaticFormalCharge(true);
  mol.AddHydrogens(false, true, 9.0); // NH3+ COO- NH2
  cout << "#pH = 9.0 : ";
  conv.Write(&mol, &cout);
 
  if (mol.NumAtoms() == 24)
    cout << "ok 15\n";
  else
    cout << "not ok 15\n";

  conv.ReadString(&mol, "NC(CCCCN)C(O)=O");
  mol.SetAutomaticFormalCharge(true);
  mol.AddHydrogens(false, true, 13.0); // NH2 COO- NH2
  cout << "#pH = 13.0 : ";
  conv.Write(&mol, &cout);

  if (mol.NumAtoms() == 23)
    cout << "ok 16\n";
  else
    cout << "not ok 16\n";

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
    cout << "ok 17\n";
  else
    cout << "not ok 17\n";

  conv.ReadString(&mol, "NC(Cc1ccc(O)cc1)C(O)=O");
  mol.SetAutomaticFormalCharge(true);
  mol.AddHydrogens(false, true, 7.4); // NH3+ COO- NH3+
  cout << "#pH = 7.4 : ";
  conv.Write(&mol, &cout);
 
  if (mol.NumAtoms() == 24)
    cout << "ok 18\n";
   else
    cout << "not ok 18\n";

  conv.ReadString(&mol, "NC(Cc1ccc(O)cc1)C(O)=O");
  mol.SetAutomaticFormalCharge(true);
  mol.AddHydrogens(false, true, 10.05); // NH3+ COO- NH2
  cout << "#pH = 10.05 : ";
  conv.Write(&mol, &cout);
 
  if (mol.NumAtoms() == 23)
    cout << "ok 19\n";
  else
    cout << "not ok 19\n";

  conv.ReadString(&mol, "NC(Cc1ccc(O)cc1)C(O)=O");
  mol.SetAutomaticFormalCharge(true);
  mol.AddHydrogens(false, true, 13.0); // NH2 COO- NH2
  cout << "#pH = 13.0 : ";
  conv.Write(&mol, &cout);

  if (mol.NumAtoms() == 22)
    cout << "ok 20\n";
  else
    cout << "not ok 20\n";

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
    cout << "ok 21\n";
  else
    cout << "not ok 21\n";

  conv.ReadString(&mol, "NC(CCCNC(N)=N)C(O)=O");
  mol.SetAutomaticFormalCharge(true);
  mol.AddHydrogens(false, true, 7.4); // NH3+ COO- NH3+
  cout << "#pH = 7.4 : ";
  conv.Write(&mol, &cout);
 
  if (mol.NumAtoms() == 27)
    cout << "ok 22\n";
   else
    cout << "not ok 22\n";

  conv.ReadString(&mol, "NC(CCCNC(N)=N)C(O)=O");
  mol.SetAutomaticFormalCharge(true);
  mol.AddHydrogens(false, true, 11.0); // NH3+ COO- NH2
  cout << "#pH = 11.0 : ";
  conv.Write(&mol, &cout);
 
  if (mol.NumAtoms() == 26)
    cout << "ok 23\n";
  else
    cout << "not ok 23\n";

  conv.ReadString(&mol, "NC(CCCNC(N)=N)C(O)=O");
  mol.SetAutomaticFormalCharge(true);
  mol.AddHydrogens(false, true, 13.0); // NH2 COO- NH2
  cout << "#pH = 13.0 : ";
  conv.Write(&mol, &cout);

  if (mol.NumAtoms() == 25)
    cout << "ok 24\n";
  else
    cout << "not ok 24\n";

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
    cout << "ok 25\n";
  else
    cout << "not ok 25\n";

  conv.ReadString(&mol, "NCC(=O)NCC(=O)O");
  mol.SetAutomaticFormalCharge(true);
  mol.AddHydrogens(false, true, 7.4); // NH3+ COO- NH3+
  cout << "#pH = 7.4 : ";
  conv.Write(&mol, &cout);
 
  if (mol.NumAtoms() == 17)
    cout << "ok 26\n";
   else
    cout << "not ok 26\n";

  conv.ReadString(&mol, "NCC(=O)NCC(=O)O");
  mol.SetAutomaticFormalCharge(true);
  mol.AddHydrogens(false, true, 13.0); // NH2 COO- NH2
  cout << "#pH = 13.0 : ";
  conv.Write(&mol, &cout);

  if (mol.NumAtoms() == 16)
    cout << "ok 27\n";
  else
    cout << "not ok 27\n";



//cout << mol.NumAtoms() << endl;  

}

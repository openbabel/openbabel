/**********************************************************************
cleartitletest.cpp - OBMol::Clear drops the molecule title.

This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
***********************************************************************/

#include "obtest.h"

#include <openbabel/atom.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>

#include <cstdio>
#include <cstdlib>
#include <string>

using namespace OpenBabel;

static void test_clear_drops_title()
{
  OBMol mol;
  mol.SetTitle("ethanol");
  mol.NewAtom()->SetAtomicNum(6);
  mol.SetIsPatternStructure();
  OB_REQUIRE(mol.NumAtoms() == 1);
  OB_REQUIRE(std::string(mol.GetTitle()) == "ethanol");

  OB_REQUIRE(mol.Clear());
  OB_REQUIRE(mol.NumAtoms() == 0);
  OB_REQUIRE(mol.HasFlag(OB_PATTERN_STRUCTURE));
  OB_REQUIRE(std::string(mol.GetTitle()).empty());
}

static void test_reused_smiles_has_no_leftover_title()
{
  OBConversion conv;
  OB_REQUIRE(conv.SetInFormat("smi"));

  OBMol mol;
  OB_REQUIRE(conv.ReadString(&mol, "CCO ethanol\n"));
  OB_REQUIRE(std::string(mol.GetTitle()) == "ethanol");
  OB_REQUIRE(mol.GetFormula() == "C2H6O");

  OB_REQUIRE(mol.Clear());
  OB_REQUIRE(conv.ReadString(&mol, "c1ccccc1\n"));
  OB_REQUIRE(mol.GetFormula() == "C6H6");
  OB_REQUIRE(std::string(mol.GetTitle()).empty());

  OBMol fresh;
  OB_REQUIRE(conv.ReadString(&fresh, "c1ccccc1\n"));
  OB_REQUIRE(fresh.GetFormula() == "C6H6");
  OB_REQUIRE(std::string(fresh.GetTitle()).empty());
}

int cleartitletest(int argc, char* argv[])
{
  int choice = 1;
  if (argc > 1 && sscanf(argv[1], "%d", &choice) != 1) {
    printf("Couldn't parse that input as a number\n");
    return -1;
  }

#ifdef FORMATDIR
  char env[BUFF_SIZE];
  snprintf(env, BUFF_SIZE, "BABEL_LIBDIR=%s", FORMATDIR);
  putenv(env);
#endif

  if (choice == 1)
    test_clear_drops_title();
  else if (choice == 2)
    test_reused_smiles_has_no_leftover_title();
  else
    return -1;
  return 0;
}

/**********************************************************************
fuzzregresstest.cpp - Regression tests for previously crashing inputs.

Each test case loads a saved fuzzer reproducer (typically from OSS-Fuzz
or a CVE report) and runs it through OBConversion. The test passes as
long as the read returns cleanly. Run under ASAN/UBSAN to catch the
original memory-safety issue if it regresses.

Add new cases by:
  1. Saving the minimized reproducer to test/files/fuzz_regress/
     using the naming convention cve-YYYY-NNNNN.<ext>
  2. Adding a case<N>() function below
  3. Registering it in the switch and incrementing
     fuzzregresstest_parts in test/CMakeLists.txt

Copyright (C) 2026 by the Open Babel project.

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

#include "obtest.h"
#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>

#include <fstream>
#include <sstream>
#include <string>

using namespace std;
using namespace OpenBabel;

static string GetFuzzFile(const string &filename)
{
  return string(TESTDATADIR) + "fuzz_regress/" + filename;
}

// Run a single reproducer through the named input format. We don't care
// whether the read succeeds or fails; we only require that it returns
// without crashing. Skip silently if the corpus file isn't available
// (e.g. submodule not checked out), so the test is non-fatal in that
// environment but still meaningful in CI where the file is present.
static bool RunRepro(const string &cveId, const string &inFormat,
                     const string &filename)
{
  string path = GetFuzzFile(filename);
  ifstream probe(path.c_str());
  if (!probe.good()) {
    cout << "# skip " << cveId << ": corpus file missing (" << path << ")\n";
    return true;
  }

  OBConversion conv;
  if (!conv.SetInFormat(inFormat.c_str())) {
    cout << "# skip " << cveId << ": format " << inFormat
         << " not registered in this build\n";
    return true;
  }

  OBMol mol;
  // ReadFile may return true or false; both are fine. The point is the
  // sanitizer must not fire while parsing this previously-crashing input.
  conv.ReadFile(&mol, path);
  return true;
}

// Read a reproducer in one format and round-trip it through a writer
// in another format. As with RunRepro, success or failure of the read
// or write is irrelevant -- the goal is to ensure neither path crashes
// or trips a sanitizer.
static bool RunReproConvert(const string &caseId, const string &inFormat,
                            const string &outFormat, const string &filename)
{
  string path = GetFuzzFile(filename);
  ifstream probe(path.c_str());
  if (!probe.good()) {
    cout << "# skip " << caseId << ": corpus file missing (" << path << ")\n";
    return true;
  }

  OBConversion conv;
  if (!conv.SetInFormat(inFormat.c_str())) {
    cout << "# skip " << caseId << ": input format " << inFormat
         << " not registered in this build\n";
    return true;
  }
  if (!conv.SetOutFormat(outFormat.c_str())) {
    cout << "# skip " << caseId << ": output format " << outFormat
         << " not registered in this build\n";
    return true;
  }

  OBMol mol;
  if (!conv.ReadFile(&mol, path))
    return true;

  ostringstream out;
  conv.Write(&mol, &out);
  return true;
}

// CVE-2026-2704: heap-buffer-overflow in transform3d::DescribeAsString
// when parsing a CIF with an all-zero row in a space-group transform.
// Fixed in PR #2862.
void caseCVE_2026_2704()
{
  OB_ASSERT(RunRepro("CVE-2026-2704", "cif", "cve-2026-2704.cif"));
}

// CVE-2026-2705: NULL pointer dereference in OBAtom::SetFormalCharge via
// MOL2Format::ReadMolecule with an out-of-range UNITY_ATOM_ATTR id.
// Fixed in PR #2862.
void caseCVE_2026_2705()
{
  OB_ASSERT(RunRepro("CVE-2026-2705", "mol2", "cve-2026-2705.mol2"));
}

// CVE-2026-3408: NULL pointer dereference in
// ChemDrawXMLFormat::EndElement when fragment atom indices fail to
// resolve. Fixed in PR #2862.
void caseCVE_2026_3408()
{
  OB_ASSERT(RunRepro("CVE-2026-3408", "cdxml", "cve-2026-3408.cdxml"));
}

// ANT-2026-00770: crash when writing a star-shaped molecule (one
// central atom bonded to seven peripheral atoms) as MCDL. Read the
// MOL/SDF reproducer and exercise the MCDL writer.
void caseANT_2026_00770()
{
  OB_ASSERT(RunReproConvert("ANT-2026-00770", "sdf", "mcdl",
                            "ant-2026-00770.sdf"));
}

int fuzzregresstest(int argc, char *argv[])
{
  int defaultchoice = 1;
  int choice = defaultchoice;

  if (argc > 1) {
    if (sscanf(argv[1], "%d", &choice) != 1) {
      printf("Couldn't parse that input as a number\n");
      return -1;
    }
  }

#ifdef FORMATDIR
  char env[BUFF_SIZE];
  snprintf(env, BUFF_SIZE, "BABEL_LIBDIR=%s", FORMATDIR);
  putenv(env);
#endif

  switch (choice) {
  case 1:
    caseCVE_2026_2704();
    break;
  case 2:
    caseCVE_2026_2705();
    break;
  case 3:
    caseCVE_2026_3408();
    break;
  case 4:
    caseANT_2026_00770();
    break;
  default:
    cout << "Test number " << choice << " does not exist!\n";
    return -1;
  }

  return 0;
}

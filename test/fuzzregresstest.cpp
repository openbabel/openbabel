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

// Like RunRepro but sets one INOPTION flag before reading.  Needed for
// formats that gate a code path behind a conversion option (e.g. mol2 -c).
static bool RunReproWithInputFlag(const string &cveId, const string &inFormat,
                                  const string &filename, const string &flag)
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
  conv.AddOption(flag.c_str(), OBConversion::INOPTIONS);

  OBMol mol;
  conv.ReadFile(&mol, path);
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

// CVE-2022-46291: out-of-bounds write into a fixed 3-element
// translationVectors[] in GaussianOutputFormat when the orientation
// block contains more than three atomicNum=-2 (Tv) rows.
void caseCVE_2022_46291()
{
  OB_ASSERT(RunRepro("CVE-2022-46291", "g09", "cve-2022-46291.g09"));
}

// CVE-2022-46292: out-of-bounds write into translationVectors[] in
// MOPACFormat when the "UNIT CELL TRANSLATION" block contains more
// than three lattice-vector rows.
void caseCVE_2022_46292()
{
  OB_ASSERT(RunRepro("CVE-2022-46292", "mopout", "cve-2022-46292.out"));
}

// CVE-2022-46293: out-of-bounds write into translationVectors[] in
// MOPACFormat when the "FINAL POINT AND DERIVATIVES" block contains
// more than three Tv-atom Z-component rows.
void caseCVE_2022_46293()
{
  OB_ASSERT(RunRepro("CVE-2022-46293", "mopout", "cve-2022-46293.out"));
}

// CVE-2022-46294: out-of-bounds write into translationVectors[] in
// MOPACCARTFormat when the input contains more than three Tv-element
// atom rows.
void caseCVE_2022_46294()
{
  OB_ASSERT(RunRepro("CVE-2022-46294", "mop", "cve-2022-46294.mop"));
}

// CVE-2022-46295: out-of-bounds write into translationVectors[] in
// MSIFormat when a PeriodicType record is followed by more than
// three lattice-vector lines.
void caseCVE_2022_46295()
{
  OB_ASSERT(RunRepro("CVE-2022-46295", "msi", "cve-2022-46295.msi"));
}

// CVE-2022-46289: out-of-bounds write into the confCoords[] heap
// buffer in OrcaOutputFormat when the "Number of atoms" header
// understates the row count of the following CARTESIAN COORDINATES
// (ANGSTROEM) block.
void caseCVE_2022_46289()
{
  OB_ASSERT(RunRepro("CVE-2022-46289", "orca", "cve-2022-46289.out"));
}

// CVE-2022-46290: out-of-bounds write in OrcaOutputFormat reachable
// via a malformed "Number of atoms" value (e.g. negative) that
// previously skipped the confCoords[] bounds check.
void caseCVE_2022_46290()
{
  OB_ASSERT(RunRepro("CVE-2022-46290", "orca", "cve-2022-46290.out"));
}

// CVE-2022-42885: uninitialized OBResidue* in GROFormat::WriteMolecule
// when the molecule has no residue information (e.g. read from XYZ).
// The write path dereferenced `res` without a null check.
void caseCVE_2022_42885()
{
  OB_ASSERT(RunReproConvert("CVE-2022-42885", "xyz", "gro",
                            "cve-2022-42885.xyz"));
}

// CVE-2022-44451: uninitialized OBAtom* in MSIFormat::ReadMolecule
// when an atom record contains an XYZ line before the ACL line that
// allocates the atom object.
void caseCVE_2022_44451()
{
  OB_ASSERT(RunRepro("CVE-2022-44451", "msi", "cve-2022-44451.msi"));
}

// CVE-2022-46280: uninitialized OBFormat* in PQSFormat::ReadMolecule
// when a "geom file=" line references an external file whose suffix
// does not match any of the recognized =car/=hin/=pdb/=mop patterns,
// leaving pFormat garbage before the dispatch call.
void caseCVE_2022_46280()
{
  OB_ASSERT(RunRepro("CVE-2022-46280", "pqs", "cve-2022-46280.pqs"));
}

// CVE-2022-43467: out-of-bounds write in PQSFormat::ReadMolecule when
// a relative "geom file=" name is long enough that the directory prefix
// copied from the title plus the filename together overflow the 256-byte
// full_coord_path[] buffer in the strcat() call.
void caseCVE_2022_43467()
{
  OB_ASSERT(RunRepro("CVE-2022-43467", "pqs", "cve-2022-43467.pqs"));
}

// CVE-2025-10998: NULL dereference in ChemKinFormat::ReadReactionQualifierLines.
// A qualifier line that tokenizes to empty caused toks[0] to be accessed
// out-of-bounds; the resulting invalid std::string had _M_data()==nullptr,
// so the subsequent strcasecmp(toks[0].c_str(),...) SEGVed at address 0x0.
// Fixed in af4a4212 by adding an upfront toks.empty() early-continue.
void caseCVE_2025_10998()
{
  OB_ASSERT(RunRepro("CVE-2025-10998", "ck", "cve-2025-10998.ck"));
}

// CVE-2025-10997: heap-buffer-overflow in ChemKinFormat::CheckSpecies via
// ReadReactionQualifierLines. A "TS" qualifier line with fewer than two
// tokens caused toks[1] to be accessed out-of-bounds; the resulting garbage
// std::string was passed to IMols.find(), whose map-node traversal then
// tripped a heap-buffer-overflow. Fixed in af4a4212 by adding a
// toks.size()>=2 guard and an upfront toks.empty() early-continue.
void caseCVE_2025_10997()
{
  OB_ASSERT(RunRepro("CVE-2025-10997", "ck", "cve-2025-10997.ck"));
}

// CVE-2025-10995: memcpy-param-overlap in basic_unzip_streambuf::underflow.
// When the get pointer is close enough to the start of _buffer, the putback
// copy destination and source overlap. memcpy has undefined behaviour on
// overlapping ranges; fixed by replacing it with memmove.
void caseCVE_2025_10995()
{
  OB_ASSERT(RunRepro("CVE-2025-10995", "sdf", "cve-2025-10995.sdf.gz"));
}

// CVE-2025-10996: heap-buffer-overflow in OBSmilesParser::ParseSmiles when
// a malformed SMILES string (e.g. from a bracket path) creates atoms that
// have no corresponding _hcount entry. The loop applied implicit valence
// using _hcount[idx-1] without checking whether idx-1 is in range.
// Fixed in b34cd604 by adding an `idx-1 >= _hcount.size()` bounds check.
void caseCVE_2025_10996()
{
  OB_ASSERT(RunRepro("CVE-2025-10996", "smi", "cve-2025-10996.smi"));
}

// CVE-2025-10999: NULL dereference in CacaoFormat::SetHilderbrandt when
// atoms are spaced >10 Å apart. The original code used sum=100.0 as a
// distance² threshold, so no reference atom was ever found and vit[i]->_a
// was left nullptr; the subsequent GetIdx() call then crashed. Fixed by
// using numeric_limits::max() as the initial sum and adding null guards.
void caseCVE_2025_10999()
{
  OB_ASSERT(RunReproConvert("CVE-2025-10999", "xyz", "cacint",
                            "cve-2025-10999.xyz"));
}

// CVE-2025-11000: out-of-bounds read in PQSFormat::ReadMolecule via lowerit().
// When a '=' appeared at position i<4, the original code did
// strncpy(tmp, &s[i-4], 5) reading 5 bytes before the start of the buffer.
// Fixed by replacing strncpy+strcmp with an i<4 guard + strncmp in-place.
void caseCVE_2025_11000()
{
  OB_ASSERT(RunRepro("CVE-2025-11000", "pqs", "cve-2025-11000.pqs"));
}

// CVE-2022-37331: stack-buffer-overflow in GaussianOutputFormat::ReadMolecule
// during the orientation pre-scan. strncpy(coords_type, vs[0], 24) followed
// by strcat(coords_type, " orientation:") overflowed the 25-byte coords_type[]
// buffer when the first token on an "orientation:" line exceeded 10 chars.
// Fixed by replacing coords_type with std::string.
void caseCVE_2022_37331()
{
  OB_ASSERT(RunRepro("CVE-2022-37331", "g09", "cve-2022-37331.g09"));
}

// CVE-2022-41793: heap-buffer-overflow in CSRFormat::PadString when
// mol.GetTitle() is longer than the 100-byte output buffer. strncpy was
// called with strlen(input) as the limit instead of the buffer size,
// so a title longer than 100 chars overflowed into adjacent heap memory
// during WriteCSRHeader.
void caseCVE_2022_41793()
{
  OB_ASSERT(RunReproConvert("CVE-2022-41793", "xyz", "csr",
                            "cve-2022-41793.xyz"));
}

// CVE-2025-10994: heap-use-after-free in GAMESSOutputFormat::ReadMolecule
// when a line matching "ICHARG=" or "MULT " has fewer whitespace-separated
// tokens than expected, causing vs[1] or vs[2] to read freed vector memory.
void caseCVE_2025_10994()
{
  OB_ASSERT(RunRepro("CVE-2025-10994", "gamout", "cve-2025-10994.out"));
}

// CVE-2022-43607: stack-buffer-overflow in MOL2Format::ReadMolecule when
// parsing a UCSF Dock "##########" comment line with an attribute or value
// field longer than 31 chars.  sscanf %[^:] and %s had no width limit,
// overflowing the 32-byte attr[] and val[] stack buffers.
// Requires the -c INOPTION (UCSF Dock comment mode) to reach the sscanf.
void caseCVE_2022_43607()
{
  OB_ASSERT(RunReproWithInputFlag("CVE-2022-43607", "mol2",
                                  "cve-2022-43607.mol2", "c"));
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
  case 5:
    caseCVE_2022_46291();
    break;
  case 6:
    caseCVE_2022_46292();
    break;
  case 7:
    caseCVE_2022_46293();
    break;
  case 8:
    caseCVE_2022_46294();
    break;
  case 9:
    caseCVE_2022_46295();
    break;
  case 10:
    caseCVE_2022_46289();
    break;
  case 11:
    caseCVE_2022_46290();
    break;
  case 12:
    caseCVE_2022_42885();
    break;
  case 13:
    caseCVE_2022_44451();
    break;
  case 14:
    caseCVE_2022_46280();
    break;
  case 15:
    caseCVE_2022_43467();
    break;
  case 16:
    caseCVE_2022_43607();
    break;
  case 17:
    caseCVE_2025_10994();
    break;
  case 18:
    caseCVE_2022_41793();
    break;
  case 19:
    caseCVE_2022_37331();
    break;
  case 20:
    caseCVE_2025_11000();
    break;
  case 21:
    caseCVE_2025_10999();
    break;
  case 22:
    caseCVE_2025_10996();
    break;
  case 23:
    caseCVE_2025_10995();
    break;
  case 24:
    caseCVE_2025_10997();
    break;
  case 25:
    caseCVE_2025_10998();
    break;
  default:
    cout << "Test number " << choice << " does not exist!\n";
    return -1;
  }

  return 0;
}

/**********************************************************************
ketformattest.cpp - Unit tests for the Ketcher KET (JSON) format.

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

#include <openbabel/atom.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/obiter.h>
#include <openbabel/reactionfacade.h>

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>

using namespace std;
using namespace OpenBabel;

namespace {

string testFilePath(const string &name)
{
#ifdef TESTDATADIR
    return string(TESTDATADIR) + "ket/" + name;
#else
    return string("files/ket/") + name;
#endif
}

OBMol readKetFile(const string &path)
{
    OBConversion conv;
    OB_REQUIRE(conv.SetInFormat("ket"));
    OBMol mol;
    OB_REQUIRE(conv.ReadFile(&mol, path));
    return mol;
}

string writeKet(OBMol &mol)
{
    OBConversion conv;
    OB_REQUIRE(conv.SetOutFormat("ket"));
    return conv.WriteString(&mol);
}

string writeKetMinified(OBMol &mol)
{
    OBConversion conv;
    OB_REQUIRE(conv.SetOutFormat("ket"));
    conv.AddOption("m", OBConversion::OUTOPTIONS);
    return conv.WriteString(&mol);
}

OBMol readKetString(const string &text)
{
    OBConversion conv;
    OB_REQUIRE(conv.SetInFormat("ket"));
    OBMol mol;
    OB_REQUIRE(conv.ReadString(&mol, text));
    return mol;
}

}  // namespace

// ---------------------------------------------------------------------
// Test 1: parse a single molecule from a KET file.
// ---------------------------------------------------------------------

void testReadSingleMolecule()
{
    OBMol mol = readKetFile(testFilePath("ethanol.ket"));
    OB_ASSERT(mol.NumAtoms() == 3);
    OB_ASSERT(mol.NumBonds() == 2);

    unsigned int carbons = 0, oxygens = 0;
    FOR_ATOMS_OF_MOL(a, mol) {
        if (a->GetAtomicNum() == 6) ++carbons;
        else if (a->GetAtomicNum() == 8) ++oxygens;
    }
    OB_ASSERT(carbons == 2);
    OB_ASSERT(oxygens == 1);
}

// ---------------------------------------------------------------------
// Test 2: read a KET reaction and verify OBReactionFacade roles.
// ---------------------------------------------------------------------

void testReadReaction()
{
    OBMol mol = readKetFile(testFilePath("reaction.ket"));
    OB_ASSERT(mol.IsReaction());

    OBReactionFacade facade(&mol);
    // reaction.ket: ethene + Br2 -> 1,2-dibromoethane
    OB_ASSERT(facade.NumComponents(REACTANT) == 2);
    OB_ASSERT(facade.NumComponents(PRODUCT) == 1);
}

// ---------------------------------------------------------------------
// Test 3: round-trip a molecule through KET preserving atoms/bonds and
// formal charges across two passes.
// ---------------------------------------------------------------------

void testRoundTripMolecule()
{
    OBMol original = readKetFile(testFilePath("multi.ket"));
    const unsigned int origAtoms = original.NumAtoms();
    const unsigned int origBonds = original.NumBonds();

    int origCharge = 0;
    FOR_ATOMS_OF_MOL(a, original) origCharge += a->GetFormalCharge();

    // First round-trip.
    const string first = writeKet(original);
    OB_REQUIRE(!first.empty());
    OBMol passOne = readKetString(first);
    OB_ASSERT(passOne.NumAtoms() == origAtoms);
    OB_ASSERT(passOne.NumBonds() == origBonds);

    int passOneCharge = 0;
    FOR_ATOMS_OF_MOL(a, passOne) passOneCharge += a->GetFormalCharge();
    OB_ASSERT(passOneCharge == origCharge);

    // Second round-trip — output must remain stable and identical to first.
    const string second = writeKet(passOne);
    OB_ASSERT(second == first);

    OBMol passTwo = readKetString(second);
    OB_ASSERT(passTwo.NumAtoms() == origAtoms);
    OB_ASSERT(passTwo.NumBonds() == origBonds);
}

// ---------------------------------------------------------------------
// Test 4: round-trip a KET reaction. Verifies the document remains a
// reaction with the expected reactant/product counts after re-reading.
// ---------------------------------------------------------------------

void testRoundTripReaction()
{
    OBMol original = readKetFile(testFilePath("reaction.ket"));
    OB_REQUIRE(original.IsReaction());

    const string written = writeKet(original);
    OB_REQUIRE(!written.empty());
    // Sanity-check the writer emitted at least one arrow.
    OB_ASSERT(written.find("\"arrow\"") != string::npos);

    OBMol rt = readKetString(written);
    OB_ASSERT(rt.IsReaction());

    OBReactionFacade facade(&rt);
    OB_ASSERT(facade.NumComponents(REACTANT) == 2);
    OB_ASSERT(facade.NumComponents(PRODUCT) == 1);
}

// ---------------------------------------------------------------------
// Test 5: ref-name preservation — a KET with a non-default $ref ("mol5")
// must keep that name through a round-trip.
// ---------------------------------------------------------------------

void testRefNamePreserved()
{
    OBMol mol = readKetFile(testFilePath("with_connections.ket"));
    const string out = writeKet(mol);
    OB_ASSERT(out.find("\"$ref\":\"mol5\"") != string::npos ||
              out.find("\"$ref\": \"mol5\"") != string::npos);
    OB_ASSERT(out.find("\"mol5\":") != string::npos);
}

// ---------------------------------------------------------------------
// Test 6: root.connections preservation — must survive a full round-trip.
// ---------------------------------------------------------------------

void testConnectionsPreserved()
{
    OBMol mol = readKetFile(testFilePath("with_connections.ket"));
    const string out = writeKet(mol);
    OB_ASSERT(out.find("\"connections\"") != string::npos);
    OB_ASSERT(out.find("monomerId") != string::npos);
    OB_ASSERT(out.find("attachmentPointId") != string::npos);

    // Re-read and re-write: stable output the second time.
    OBMol second = readKetString(out);
    const string out2 = writeKet(second);
    OB_ASSERT(out2 == out);
}

// ---------------------------------------------------------------------
// Test 7: node ordering preservation — a text meta object placed between
// two molecule refs must keep its position relative to them.
// ---------------------------------------------------------------------

void testNodeOrderPreserved()
{
    OBMol mol = readKetFile(testFilePath("with_connections.ket"));
    const string out = writeKet(mol);
    const auto molPos = out.find("\"$ref\":\"mol5\"");
    const auto textPos = out.find("\"type\":\"text\"");
    const auto monomerPos = out.find("\"$ref\":\"monomer0\"");
    const auto molPosAlt = out.find("\"$ref\": \"mol5\"");
    const auto textPosAlt = out.find("\"type\": \"text\"");
    const auto monomerPosAlt = out.find("\"$ref\": \"monomer0\"");
    const auto mp = (molPos != string::npos) ? molPos : molPosAlt;
    const auto tp = (textPos != string::npos) ? textPos : textPosAlt;
    const auto np = (monomerPos != string::npos) ? monomerPos : monomerPosAlt;
    OB_REQUIRE(mp != string::npos && tp != string::npos && np != string::npos);
    OB_ASSERT(mp < tp);  // mol5 first
    OB_ASSERT(tp < np);  // then text, then monomer
}

// ---------------------------------------------------------------------
// Test 8: original reaction arrow geometry must be preserved (not
// synthesized from centroids).
// ---------------------------------------------------------------------

void testOriginalArrowPreserved()
{
    OBMol mol = readKetFile(testFilePath("reaction.ket"));
    OB_REQUIRE(mol.IsReaction());
    const string out = writeKet(mol);
    // The original arrow runs from x=2 to x=5. The synthesizer would emit
    // tail+gap (~1.1) to head-gap (~5.9). So the presence of x=2 and x=5
    // proves we kept the original.
    OB_ASSERT(out.find("\"x\":2.0") != string::npos ||
              out.find("\"x\": 2.0") != string::npos);
    OB_ASSERT(out.find("\"x\":5.0") != string::npos ||
              out.find("\"x\": 5.0") != string::npos);
    // Exactly one arrow — no duplicate synthesis.
    size_t count = 0, pos = 0;
    while ((pos = out.find("\"type\":\"arrow\"", pos)) != string::npos) {
        ++count; ++pos;
    }
    if (count == 0) {
        pos = 0;
        while ((pos = out.find("\"type\": \"arrow\"", pos)) != string::npos) {
            ++count; ++pos;
        }
    }
    OB_ASSERT(count == 1);
}

// Helper: count occurrences of a substring.
static size_t countOccurrences(const string &haystack, const string &needle)
{
    size_t count = 0, pos = 0;
    while ((pos = haystack.find(needle, pos)) != string::npos) {
        ++count; ++pos;
    }
    return count;
}

// ---------------------------------------------------------------------
// Test 9: S-group round-trip with bond indices for a multi-fragment file.
// Verifies that S-group `bonds: [...]` indices are correctly remapped
// back to per-molecule-local KET indices (not left as global).
// ---------------------------------------------------------------------

void testSgroupBondsLocal()
{
    OBMol mol = readKetFile(testFilePath("full_fidelity.ket"));
    const string out = writeKetMinified(mol);

    // Two molecules. mol1's bonds array has exactly one bond (index 0).
    // mol1's S-group must reference "bonds":[0] (local), not the global
    // bond index (which would be 3 after mol0's three bonds).
    const size_t mol0Def = out.find("\"mol0\":");
    const size_t mol1Def = out.find("\"mol1\":");
    OB_REQUIRE(mol0Def != string::npos);
    OB_REQUIRE(mol1Def != string::npos);
    const string mol0Body = out.substr(mol0Def, mol1Def - mol0Def);
    const string mol1Body = out.substr(mol1Def);

    OB_ASSERT(mol0Body.find("\"bonds\":[0,2]") != string::npos);
    OB_ASSERT(mol1Body.find("\"bonds\":[0]") != string::npos);
    OB_ASSERT(mol1Body.find("\"bonds\":[3]") == string::npos);
}

// ---------------------------------------------------------------------
// Test 10: fragment-level KET properties and highlight survive round-trip.
// ---------------------------------------------------------------------

void testPropertiesAndHighlightPreserved()
{
    OBMol mol = readKetFile(testFilePath("full_fidelity.ket"));
    const string out = writeKetMinified(mol);

    // mol0 has both highlight and properties; mol1 has neither. Properties
    // / highlight must NOT leak from mol0 into mol1.
    const size_t mol0Def = out.find("\"mol0\":");
    const size_t mol1Def = out.find("\"mol1\":");
    OB_REQUIRE(mol0Def != string::npos);
    OB_REQUIRE(mol1Def != string::npos);
    const string mol0Body = out.substr(mol0Def, mol1Def - mol0Def);
    const string mol1Body = out.substr(mol1Def);

    OB_ASSERT(mol0Body.find("\"highlight\"") != string::npos);
    OB_ASSERT(mol0Body.find("\"properties\"") != string::npos);
    OB_ASSERT(mol0Body.find("\"key\":\"name\"") != string::npos);
    OB_ASSERT(mol0Body.find("\"value\":\"test-compound\"") != string::npos);
    OB_ASSERT(mol0Body.find("\"entityType\":\"atom\"") != string::npos);
    OB_ASSERT(mol0Body.find("\"entityType\":\"bond\"") != string::npos);

    OB_ASSERT(mol1Body.find("\"highlight\"") == string::npos);
    OB_ASSERT(mol1Body.find("\"properties\"") == string::npos);
}

// ---------------------------------------------------------------------
// Test 11: byte-stable round-trip — for a file that uses only
// natively-supported features, reading then writing twice must produce
// identical output.
// ---------------------------------------------------------------------

void testStableRoundTrip()
{
    OBMol mol = readKetFile(testFilePath("full_fidelity.ket"));
    const string first = writeKet(mol);
    OBMol pass1 = readKetString(first);
    const string second = writeKet(pass1);
    OB_ASSERT(first == second);

    OBMol pass2 = readKetString(second);
    const string third = writeKet(pass2);
    OB_ASSERT(second == third);
}

// ---------------------------------------------------------------------
// Test 12: empty OBMol -> KET emits a structurally valid empty document
// (no dangling refs).
// ---------------------------------------------------------------------

void testEmptyMolEmitsNoDanglingRef()
{
    OBMol mol;  // empty
    const string out = writeKet(mol);
    // No mol* members nor refs.
    OB_ASSERT(out.find("\"$ref\"") == string::npos);
    OB_ASSERT(out.find("\"mol0\"") == string::npos);
    // But the document is still valid: contains a root.nodes (possibly empty).
    OB_ASSERT(out.find("\"root\"") != string::npos);
    OB_ASSERT(out.find("\"nodes\"") != string::npos);

    // Re-read should not fail.
    OBMol roundtrip = readKetString(out);
    OB_ASSERT(roundtrip.NumAtoms() == 0);
}

// ---------------------------------------------------------------------
// Test 13: no duplicate arrow on reaction round-trip — the original arrow
// must be preserved AND not duplicated by the synthesizer.
// ---------------------------------------------------------------------

void testNoDuplicateArrow()
{
    OBMol mol = readKetFile(testFilePath("reaction.ket"));
    const string out = writeKet(mol);
    const size_t arrowCount =
        countOccurrences(out, "\"type\":\"arrow\"") +
        countOccurrences(out, "\"type\": \"arrow\"");
    OB_ASSERT(arrowCount == 1);
    const size_t plusCount =
        countOccurrences(out, "\"type\":\"plus\"") +
        countOccurrences(out, "\"type\": \"plus\"");
    OB_ASSERT(plusCount == 1);  // reaction.ket has 1 original plus
}

// ---------------------------------------------------------------------
// Test 14: aromatic KET bonds remain KET type 4 after round-trip, even
// though the reader kekulizes the OBMol internally.
// ---------------------------------------------------------------------

void testAromaticBondTypePreserved()
{
    OBMol mol = readKetFile(testFilePath("aromatic.ket"));
    const string out = writeKetMinified(mol);

    OB_ASSERT(countOccurrences(out, "\"type\":4") == 6);
    OB_ASSERT(countOccurrences(out, "\"type\":1") == 0);
    OB_ASSERT(countOccurrences(out, "\"type\":2") == 0);

    OBMol second = readKetString(out);
    const string out2 = writeKetMinified(second);
    OB_ASSERT(out2 == out);
}

// ---------------------------------------------------------------------
//  Entry point — the test harness uses cpptests=1..14 (see CMakeLists).
// ---------------------------------------------------------------------

int ketformattest(int argc, char *argv[])
{
#ifdef FORMATDIR
    char env[BUFF_SIZE];
    snprintf(env, BUFF_SIZE, "BABEL_LIBDIR=%s", FORMATDIR);
    putenv(env);
#endif

    int choice = 1;
    if (argc > 1) {
        if (sscanf(argv[1], "%d", &choice) != 1) {
            printf("Couldn't parse that input as a number\n");
            return -1;
        }
    }

    switch (choice) {
    case 1:  testReadSingleMolecule();             break;
    case 2:  testReadReaction();                   break;
    case 3:  testRoundTripMolecule();              break;
    case 4:  testRoundTripReaction();              break;
    case 5:  testRefNamePreserved();               break;
    case 6:  testConnectionsPreserved();           break;
    case 7:  testNodeOrderPreserved();             break;
    case 8:  testOriginalArrowPreserved();         break;
    case 9:  testSgroupBondsLocal();               break;
    case 10: testPropertiesAndHighlightPreserved(); break;
    case 11: testStableRoundTrip();                break;
    case 12: testEmptyMolEmitsNoDanglingRef();     break;
    case 13: testNoDuplicateArrow();               break;
    case 14: testAromaticBondTypePreserved();      break;
    default:
        cout << "Test number " << choice << " does not exist!\n";
        return -1;
    }

    return 0;
}

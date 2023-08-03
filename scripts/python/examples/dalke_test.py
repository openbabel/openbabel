# This file is meant primarily for people who want to see an example
# of how to use part of the OpenBabel API, with a secondary use as a
# set of quick unit tests to make sure there's no strange but obvious
# problem with your OpenBabel setup.

# Regression tests, coverage tests, stress tests, performance tests,
# etc. should not go in this file.


import math
import os
import shutil
import tempfile
import time
import unittest
import warnings

from openbabel import openbabel as ob

def testfile(name):
    return os.path.join("files", name)

class MyTestCase(unittest.TestCase):
    def assertClose(self, val, expect):
        if expect > 0:
            self.assertTrue((expect * 0.9999) < val < (expect * 1.0001), val)
        else:
            self.assertTrue((expect * 1.0001) < val < (expect * 0.9999), val)

    def assertZero(self, val):
        self.assertTrue(abs(val) < 0.00001, val)


# Make a temporary directory for use during the "with" context block.
# When finished, remove the directory.
class TempDir(object):
    def __init__(self):
        self.dirname = None
    def __enter__(self):
        self.dirname = tempfile.mkdtemp(prefix="ob_py_test")
        return self
    def __call__(self, name):
        return os.path.join(self.dirname, name)
    def __exit__(self, exc, value, tb):
        shutil.rmtree(self.dirname)

# Some of the API calls generate log messages. I don't want to
# see them when doing the testing.
class SuppressLogging(object):
    def __enter__(self):
        self.lvl = ob.obErrorLog.GetOutputLevel()
        ob.obErrorLog.SetOutputLevel(-1)
        return self
    def __exit__(self, exc, value, tb):
        ob.obErrorLog.SetOutputLevel(self.lvl)

# The plugin system requires that OBConversion be called first.
# This is done once, and it affects the entire system
# Check for that case now
result = ob.OBPlugin.ListAsString("fingerprints")
assert "FP2" not in result, result
ob.OBConversion()
result = ob.OBPlugin.ListAsString("fingerprints")
assert "FP2" in result, result

_smiles_parser = ob.OBConversion()
_smiles_parser.SetInFormat("smi")
_smiles_parser.SetOutFormat("can")

def parse_smiles(smiles):
    "parse a SMILES into a molecule"
    mol = ob.OBMol()
    _smiles_parser.ReadString(mol, smiles)
    return mol

def cansmiles(mol):
    "return the canonical SMILES for molecule"
    return _smiles_parser.WriteString(mol).strip()

def parse_smarts(smarts):
    "Parse a SMARTS into an ObSmartsPattern"
    pat = ob.OBSmartsPattern()
    assert pat.Init(smarts)
    return pat

def readfile(filename, filetype):
    "Iterate over all molecule records in the given file, with given type"
    if "/" not in filename and "\\" not in filename:
        filename = testfile(filename)
    obconversion = ob.OBConversion()
    obconversion.SetInFormat(filetype)
    mol = ob.OBMol()
    notatend = obconversion.ReadFile(mol, filename)
    while notatend:
        yield mol
        mol.Clear()
        notatend = obconversion.Read(mol)

class TestIO(MyTestCase):
    def test_read_smiles(self):
        it = readfile("FormulaTest.smi", "smi")
        mol = it.next()
        self.assertEqual(mol.GetTitle(), "CH4")
        self.assertTrue(mol.NumAtoms())
        mol = it.next()
        self.assertEqual(mol.GetTitle(), "C atom")
        self.assertTrue(mol.NumAtoms())

    def test_read_sdf(self):
        it = readfile("cantest.sdf", "sdf")
        mol = it.next()
        self.assertEqual(mol.GetTitle(), "8978")
        self.assertEqual(mol.NumBonds(), 64)
        mol = it.next()
        self.assertEqual(mol.GetTitle(), "10617")
        self.assertEqual(mol.NumBonds(), 40)

    def test_read_sdf_gz(self):
        it = readfile("ziptest.sdf.gz", "sdf")
        mol = it.next()
        self.assertEqual(mol.GetTitle(), "ZINC04985529")
        self.assertEqual(mol.NumAtoms(), 49)
        mol = it.next()
        self.assertEqual(mol.GetTitle(), "ZINC01700999")
        self.assertEqual(mol.NumAtoms(), 34)

    def test_write_smiles(self):
        conv = ob.OBConversion()
        conv.SetOutFormat("smi")
        with TempDir() as tempdir:
            mol = parse_smiles("CCO")
            mol.SetTitle("#1")
            conv.WriteFile(mol, tempdir("blah.smi"))
            mol = parse_smiles("[NH4+]")
            mol.SetTitle("mol2")
            conv.Write(mol)
            conv.CloseOutFile()

            lines = open(tempdir("blah.smi"), "U").readlines()
            self.assertTrue(lines[0] == "CCO\t#1\n" or
                              lines[0] == "OCC\t#1\n", repr(lines[0]))
            self.assertTrue(lines[1] == "[NH4+]\tmol2\n", repr(lines[1]))

    def test_write_sdf(self):
        conv = ob.OBConversion()
        conv.SetOutFormat("sdf")
        with TempDir() as tempdir:
            mol = parse_smiles("CCO")
            mol.SetTitle("#1")
            with SuppressLogging():
                # XXX For some reason, this generates the warning
                #   Warning in WriteMolecule No 2D or 3D coordinates exist.
                #   Any stereochemical information will be lost. To generate
                #   2D or 3D coordinates use --gen2D or --gen3d.
                # Since not all users of the API will have a --gen2D/--gen3d option,
                # that's not always going to be useful. Plus, my test cases
                # have no stereochemical information. Oh, and hey - I don't even
                # call WriteMolecule directly
                conv.WriteFile(mol, tempdir("blah.sdf"))
            mol = parse_smiles("[NH4+]")
            mol.SetTitle("mol2")
            conv.Write(mol)
            conv.CloseOutFile()

            titles = []
            atom_counts = []
            for mol in readfile(tempdir("blah.sdf"), "sdf"):
                titles.append(mol.GetTitle())
                atom_counts.append(mol.NumAtoms())
            self.assertEqual(titles, ["#1", "mol2"])
            self.assertEqual(atom_counts, [3, 5])

    def test_write_inchi(self):
        mol = parse_smiles("c1ccccc1O")
        conv = ob.OBConversion()
        conv.SetOutFormat("inchi")
        s = conv.WriteString(mol)
        # Note the newline!
        self.assertEqual(s, "InChI=1S/C6H6O/c7-6-4-2-1-3-5-6/h1-5,7H\n")

    def test_perception_and_canonicalization(self):
        mol = parse_smiles("C1=CC=C(O)C=C1")
        conv = ob.OBConversion()
        # Input does perception. Output is not canonical
        conv.SetOutFormat("smi")
        s = conv.WriteString(mol)
        self.assertEqual(s, "c1ccc(O)cc1\t\n")

        conv = ob.OBConversion()
        # Perception and canonical generation
        conv.SetOutFormat("can")
        s = conv.WriteString(mol)
        self.assertEqual(s, "Oc1ccccc1\t\n")


class TestPlugins(MyTestCase):
    known_types = ["charges", "descriptors", "fingerprints", "forcefields",
                   "formats", "loaders", "ops"]
    def test_known_types(self):
        for name in TestPlugins.known_types:
            s = ob.OBPlugin.ListAsString(name)
            self.assertFalse("not a recognized" in s, s)

            v = ob.vectorString()
            ob.OBPlugin.ListAsVector(name, None, v)
            self.assertTrue(len(v) > 0, list(v))

    def test_as_string(self):
        s = ob.OBPlugin.ListAsString("fingerprints")
        self.assertTrue("FP2" in s, s)
        self.assertTrue("FP3" in s, s)
        self.assertTrue("MACCS" in s, s)

    def test_as_string_unknown_type(self):
        s = ob.OBPlugin.ListAsString("qwerty.shrdlu")
        self.assertTrue("\nfingerprints\n" in s, s)
        self.assertTrue("\nloaders\n" in s, s)

    def test_as_vector(self):
        v = ob.vectorString()
        ob.OBPlugin.ListAsVector("formats", None, v)
        formats = set(v)
        self.assertTrue("smiles -- SMILES format" in formats, formats)

##     def test_list(self):
##         # XXX GRR! To capture requires passing a 3rd argument which is a std:ostream
##         # I can't figure out how to do that in OpenBabel
##         s = ob.OBFingerprint.List("fingerprints").splitlines(True)
##         self.assertEquals(s[0], "FP2    Indexes linear fragments up to 7 atoms.\n")
##         self.assertEquals(s[1], "FP3    SMARTS patterns specified in the file "
##                           "patterns.txt\n")
##         self.assertEquals(s[2], "FP4    SMARTS patterns specified in the file "
##                           "SMARTS_InteLigand.txt\n")
##         self.assertEquals(s[3], "MACCS    SMARTS patterns specified in the file MACCS.txt\n")
##         self.assertEquals(len(s)>=4, True, s)

class TestFingerprints(MyTestCase):
    def test_descriptions(self):
        P = "\nPatternFP is definable"
        for name, expected_description in (
            ("FP2", "Indexes linear fragments up to 7 atoms."),
            ("FP3", "SMARTS patterns specified in the file patterns.txt" + P),
            ("FP4", "SMARTS patterns specified in the file SMARTS_InteLigand.txt" + P),
            ("MACCS", "SMARTS patterns specified in the file MACCS.txt" + P)):
            fingerprinter = ob.OBFingerprint.FindFingerprint(name)
            self.assertFalse(fingerprinter is None)
            self.assertEqual(fingerprinter.GetID(), name)
            self.assertEqual(fingerprinter.Description(), expected_description)
            # Which supported platforms have non-32-bit integers?
            self.assertEqual(fingerprinter.Getbitsperint(), 32)

    # XXX I don't think DescribeBits is accessible from Python
    # XXX What do I do with the result of GetMap?
    # XXX What are "Flags()" for?

    def test_fp_words(self):
        mol = parse_smiles("c1ccccc1O.C#N.[Ge].C1CCC1")
        def next_highest_power_of_two(n):
            i = 8
            while i < n:
                i *= 2
            return i
        for (name, nbits, v0, v1) in ( ("FP2", 1021, 0, 1),
                                       ("FP3", 55, 67108864, 1159170),
                                       ("FP4", 307, 2, 0),
                                       # TODO: change my MACCS.txt so it's correct
                                       # then rerun this test and change to the right answer
                                       ("MACCS", 166, 2097156, 256),
                                       ):
            fingerprinter = ob.OBFingerprint.FindFingerprint(name)
            v = ob.vectorUnsignedInt()
            fingerprinter.GetFingerprint(mol, v)
            size = next_highest_power_of_two(nbits)//32 # bits-per-int
            self.assertEqual(len(v), size)
            self.assertEqual(v[0], v0, (name, v[0], v0))
            self.assertEqual(v[1], v1, (name, v[1], v1))

    def test_fold(self):
        v = ob.vectorUnsignedInt([0x2A, 0x41])
        self.assertEqual(len(v), 2)
        x = ob.OBFingerprint.FindFingerprint("FP2")
        x.Fold(v, 32)
        self.assertEqual(len(v), 1)
        self.assertEqual(v[0], (0x2A | 0x41))

        v = ob.vectorUnsignedInt([0x01, 0x04, 0x20, 0x00])
        self.assertEqual(len(v), 4)
        x.Fold(v, 64)
        self.assertEqual(len(v), 2)
        self.assertEqual(v[0], 0x21)
        self.assertEqual(v[1], 0x04)

    def test_get_set(self):
        v = ob.vectorUnsignedInt([1, 6])
        # XXX Why does GetBit need an actual instance?
        x = ob.OBFingerprint.FindFingerprint("FP2")
        self.assertTrue(x.GetBit(v, 0))
        for i in range(1, 32):
            self.assertFalse(x.GetBit(v, i), i)
        self.assertFalse(x.GetBit(v, 32))
        self.assertTrue(x.GetBit(v, 33))
        self.assertTrue(x.GetBit(v, 34))
        self.assertFalse(x.GetBit(v, 35))

        x.SetBit(v, 35)
        self.assertTrue(x.GetBit(v, 35))

    def test_tanimoto(self):
        v1 = ob.vectorUnsignedInt([0x1, 0x6])
        v2 = ob.vectorUnsignedInt([0x1, 0x7])
        x = ob.OBFingerprint.FindFingerprint("FP2")
        self.assertEqual(x.Tanimoto(v1, v2), (1 + 2) / (1 + 3 + 0.0))

    def test_tanimoto_size_mismatch(self):
        v1 = ob.vectorUnsignedInt([0x1, 0x6])
        v2 = ob.vectorUnsignedInt([1, 2, 0])
        x = ob.OBFingerprint.FindFingerprint("FP2")
        self.assertEqual(x.Tanimoto(v1, v2), -1.0)

    def test_tanimoto_with_no_set_bits(self):
        v1 = ob.vectorUnsignedInt([0, 0, 0, 0])
        x = ob.OBFingerprint.FindFingerprint("FP2")
        # Again, this is an arbitrary decision by toolkit providers
        self.assertEqual(x.Tanimoto(v1, v1), 0.0)


mol_with_many_rings = parse_smiles(
    "C1C3C5C7C9C%11C%13C%15C%17." +
    "C12C34C56C78C9%10C%11%12C%13%14C%15%16C%17%18." * 11 +
    "C2C4C6C8C%10C%12C%14C%16C%18C")

class TestSmarts(MyTestCase):
    def test_4_membered_ring(self):
        pat = ob.OBSmartsPattern()
        self.assertTrue(pat.Init("*1~*~*~*~1"), "failed to Init")
        mol = parse_smiles("C1CCCC1")
        m = pat.Match(mol)
        self.assertFalse(m, "had a match?")

        mol = parse_smiles("C1CCC1")
        m = pat.Match(mol)
        self.assertTrue(m, "no match?")

    def test_is_valid_and_test_empty(self):
        pat = ob.OBSmartsPattern()
        self.assertFalse(pat.IsValid())
        self.assertTrue(pat.Empty())

        pat.Init("CO")
        self.assertTrue(pat.IsValid())
        self.assertFalse(pat.Empty())
        pat = ob.OBSmartsPattern()
        with SuppressLogging():
            # This will send message to the error log.
            self.assertFalse(pat.Init("=O"))
        self.assertFalse(pat.IsValid())
        self.assertTrue(pat.Empty())

    def test_num_atoms_and_bonds(self):
        pat = ob.OBSmartsPattern()
        self.assertEqual(pat.NumAtoms(), 0)
        self.assertEqual(pat.NumBonds(), 0)
        pat.Init("C")
        self.assertEqual(pat.NumAtoms(), 1)
        self.assertEqual(pat.NumBonds(), 0)
        pat.Init("C#N")
        self.assertEqual(pat.NumAtoms(), 2)
        self.assertEqual(pat.NumBonds(), 1)
        pat.Init("c1ccccc1")
        self.assertEqual(pat.NumAtoms(), 6)
        self.assertEqual(pat.NumBonds(), 6)

    def test_basic_match_fails(self):
        mol = parse_smiles("c1ccccc1O")
        pat = parse_smarts("[#7]")
        self.assertFalse(pat.Match(mol))
        self.assertEqual(pat.NumMatches(), 0)

    def test_basic_match_fails_with_single_flag_set(self):
        mol = parse_smiles("c1ccccc1O")
        pat = parse_smarts("[#7]")
        self.assertFalse(pat.Match(mol, True))
        self.assertEqual(pat.NumMatches(), 0)

    def test_basic_match(self):
        mol = parse_smiles("c1ccccc1O")
        pat = parse_smarts("[#6][#8]")
        self.assertTrue(pat.Match(mol))
        results = pat.GetUMapList()
        self.assertEqual(len(results), 1)
        self.assertEqual(results[0], (6, 7))
        self.assertEqual(pat.NumMatches(), 1)

    def test_basic_match_with_two_unique_hits(self):
        mol = parse_smiles("c1ccccc1O")
        pat = parse_smarts("ccO")
        self.assertTrue(pat.Match(mol))
        results = pat.GetUMapList()
        self.assertEqual(len(results), 2)
        results = set(results)
        self.assertEqual(results, set([(5, 6, 7), (1, 6, 7)]))

        self.assertEqual(results, set(pat.GetMapList()))
        self.assertEqual(pat.NumMatches(), 2)

    def test_basic_match_with_one_unique_hit(self):
        mol = parse_smiles("c1ccccc1O")
        pat = parse_smarts("c1ccccc1")
        self.assertTrue(pat.Match(mol))
        results = pat.GetUMapList()
        self.assertEqual(len(results), 1)
        self.assertEqual(pat.NumMatches(), 1)
        self.assertEqual(set(results[0]), set([1, 2, 3, 4, 5, 6]))

    def test_basic_match_with_nonunique_hits(self):
        mol = parse_smiles("c1ccccc1O")
        pat = parse_smarts("c1ccccc1")
        self.assertTrue(pat.Match(mol))
        results = pat.GetMapList()
        self.assertEqual(len(results), 12)
        results = list(map(set, results))
        for i in range(12):
            self.assertEqual(results[0], results[i])

    def test_basic_match_behavior_which_I_did_not_expect(self):
        mol = parse_smiles("c1ccccc1O")
        pat = parse_smarts("c1ccccc1")
        pat.Match(mol)
        self.assertEqual(pat.NumMatches(), 12)
        results = pat.GetUMapList()
        self.assertEqual(pat.NumMatches(), 1)
        results = pat.GetMapList()
        # I really expected these to be 12.
        # It appears the UMapList does an in-place trim.
        # XXX Is that the right/expected behavior?
        self.assertEqual(pat.NumMatches(), 1)
        self.assertEqual(len(results), 1)

        pat.Match(mol)
        # Here they are 12
        results = pat.GetMapList()
        self.assertEqual(pat.NumMatches(), 12)
        results = pat.GetUMapList()
        self.assertEqual(pat.NumMatches(), 1)
        self.assertEqual(len(results), 1)

    def test_basic_match_with_single_flag_set(self):
        # I want something which takes a long time
        mol = mol_with_many_rings
        pat = parse_smarts("C1CCCCCCCCCCCCC1")
        t1 = time.time()
        self.assertTrue(pat.Match(mol))
        t2 = time.time()
        if t2 - t1 > 0.01:
            warnings.warn("test_basic_match_with_single_flag_set took too long")
        self.assertEqual(len(pat.GetMapList()), 1)
        self.assertEqual(len(pat.GetUMapList()), 1)

    def test_has_match(self):
        mol = mol_with_many_rings
        pat = parse_smarts("C1CCCCCCCCCCCCC1")
        t1 = time.time()
        self.assertTrue(pat.HasMatch(mol))
        t2 = time.time()
        if t2-t1 > 0.01:
            warnings.warn("test_has_match took too long")


    def test_vector_match_false(self):
        # Create a  vector< vector<int> >, wherein the results go
        v = ob.vectorvInt()
        mol = parse_smiles("c1ccccc1O")
        pat = parse_smarts("N")
        self.assertEqual(pat.Match(mol, v), 0)
        self.assertEqual(len(v), 0)

    def test_vector_match(self):
        v = ob.vectorvInt()
        mol = parse_smiles("c1ccccc1O")
        pat = parse_smarts("cO")
        self.assertEqual(pat.Match(mol, v), 1)
        self.assertEqual(len(v), 1)
        self.assertEqual(set(v[0]), set([6, 7]))

    def test_vector_match_with_two_hits(self):
        v = ob.vectorvInt()
        mol = parse_smiles("c1ccccc1O")
        pat = parse_smarts("ccO")
        self.assertEqual(pat.Match(mol, v), 1)
        self.assertEqual(len(v), 2)
        results = list(v)
        self.assertTrue((5, 6, 7) in results, results)
        self.assertTrue((1, 6, 7) in results, results)

    def test_vector_match_with_one_unique_hit(self):
        mol = parse_smiles("c1ccccc1O")
        pat = parse_smarts("c1ccccc1")
        v = ob.vectorvInt()
        self.assertTrue(pat.Match(mol, v, ob.OBSmartsPattern.AllUnique))
        self.assertEqual(len(v), 1)
        self.assertEqual(set(v[0]), set([1, 2, 3, 4, 5, 6]))

    def test_vector_match_with_single_hit(self):
        v = ob.vectorvInt()
        mol = parse_smiles("c1ccccc1O")
        pat = parse_smarts("ccO")
        self.assertEqual(pat.Match(mol, v, ob.OBSmartsPattern.Single), 1)
        self.assertEqual(len(v), 1)
        result = v[0]
        self.assertTrue(result == (5, 6, 7) or result == (1, 6, 7), result)

    def test_vector_match_with_all_hits(self):
        mol = parse_smiles("c1ccccc1O")
        pat = parse_smarts("c1ccccc1")
        v = ob.vectorvInt()
        self.assertEqual(pat.Match(mol, v, ob.OBSmartsPattern.All), 1)
        self.assertEqual(len(v), 12)
        expect = set([1, 2, 3, 4, 5, 6])
        for x in v:
            self.assertEqual(set(x), expect)

    def test_bad_smarts(self):
        pat = ob.OBSmartsPattern()
        # This writes an error to the log
        with SuppressLogging():
            self.assertFalse(pat.Init("%"))
        self.assertEqual(pat.NumAtoms(), 0)
        self.assertFalse(pat.IsValid())

    def test_replace_with_bad_smarts(self):
        pat = ob.OBSmartsPattern()
        self.assertTrue(pat.Init("CCCC"))
        self.assertEqual(pat.NumAtoms(), 4)
        # Re-init and verify that there's an overwrite
        # This writes an error to the log
        with SuppressLogging():
            self.assertFalse(pat.Init("Q"))
        self.assertEqual(pat.NumAtoms(), 0)
        self.assertFalse(pat.IsValid())


# The BeginMList/EndMList seems broken in Python XXX

class TestDescriptors(MyTestCase):
    def test_logp(self):
        calc_logp = ob.OBDescriptor.FindType("logP")
        mol = parse_smiles("Oc1ccccc1OC")
        #mol.AddHydrogens() # doesn't change the results
        logp = calc_logp.Predict(mol)

        self.assertTrue(abs(logp - 1.4008) <= 0.0001, logp)

    def test_tpsa(self):
        calc_tpsa = ob.OBDescriptor.FindType("TPSA")
        mol = parse_smiles("Oc1ccccc1OC")
        #mol.AddHydrogens() # doesn't change the results
        tpsa = calc_tpsa.Predict(mol)
        self.assertTrue(abs(tpsa - 29.460) <= 0.001, tpsa)

    def test_mr(self):
        calc_mr = ob.OBDescriptor.FindType("MR")
        mol = parse_smiles("Oc1ccccc1OC")
        #mol.AddHydrogens() # doesn't change the results
        mr = calc_mr.Predict(mol)
        self.assertTrue(abs(mr - 34.957) <= 0.001, mr)

    def test_gotta_try_them_all(self):
        v = ob.vectorString()
        ob.OBDescriptor.ListAsVector("descriptors", None, v)
        mol = parse_smiles("c1ccccc1O")
        for term in v:
            name = term.split()[0]
            prop_calculator = ob.OBDescriptor.FindType(name)
            self.assertFalse(prop_calculator is None, "Could not find " + name)
            prop_calculator.Predict(mol)


def add_atom(mol, atomno):
    atom = mol.NewAtom()
    atom.SetAtomicNum(atomno)
    return atom

class TestMolecule(MyTestCase):
    def test_mol_iteration(self):
        mol = parse_smiles("c12c(O[CH](C1=O)C(C)C)cc1c(c2)ccc(=O)o1")
        element_counts = {}
        for atom in ob.OBMolAtomIter(mol):
            n = atom.GetAtomicNum()
            element_counts[n] = element_counts.get(n, 0) + 1
        self.assertEqual(element_counts[8], 4)

        bond_counts = {}
        for bond in ob.OBMolBondIter(mol):
            n = bond.GetBondOrder()
            if not bond.IsAromatic():
                bond_counts[n] = bond_counts.get(n, 0) + 1
        self.assertEqual(bond_counts[2], 2)

    def test_atom_iteration(self):
        mol = parse_smiles("[U](F)(F)(F)[Cl]")
        atom = mol.GetAtom(1)

        counts = {9: 0, 17: 0}
        for bond in ob.OBAtomBondIter(atom):
            xatom = bond.GetNbrAtom(atom)
            n = xatom.GetAtomicNum()
            counts[n] += 1
        self.assertEqual(counts, {9: 3, 17: 1})

        counts = {9: 0, 17: 0}
        for atom in ob.OBAtomAtomIter(atom):
            n = atom.GetAtomicNum()
            counts[n] += 1
        self.assertEqual(counts, {9: 3, 17: 1})

### XXX By symmetry I thought something like this would work
# It does not since there is no ob.OBBondAtomIter
#    def test_bond_iteration(self):
#        mol = parse_smiles("C#N")
#        elements = []
#        for atom in ob.OBBondAtomIter(bond):
#            elements.append(atom.GetAtomicNum())
#        elements.sort()
#        self.assertEquals(elements, [6, 7])



    # Most people don't do molecule building, so I'm not going to test all the variations
    def test_building_a_molecule(self):
        mol = ob.OBMol()
        C = add_atom(mol, 6)
        N = add_atom(mol, 7)
        # XXX Why can't I do mol.AddBond(C, N, 3)?
        mol.AddBond(C.GetIdx(), N.GetIdx(), 3)
        self.assertEqual(C.ImplicitHydrogenCount(), 1)
        C.IncrementImplicitValence()  # Is this how to increment the implicit hcount?
        self.assertEqual(C.ImplicitHydrogenCount(), 2)
        conv = ob.OBConversion()
        conv.SetOutFormat("can")
        s = conv.WriteString(mol).strip()
        # XXX How does this work when the ImplicitHydrogenCount is 2?? XXX
        self.assertEqual(s, "C#N")

        # I can even add an atom this way. (Why are there 2 ways?)
        O = ob.OBAtom()
        O.SetAtomicNum(8)
        mol.AddAtom(O)
        O.SetImplicitValence(2)

        s = conv.WriteString(mol).strip()
        self.assertEqual(s, "C#N.O")

    def test_molecule_properties(self):
        # Have Cl because the average MW is 35.5 so it's easy to
        # tell the difference between the "average isotopic weight"
        # and "weight of the most common isotope" answers
        mol = parse_smiles("c1ccccc1O.[NH4+].[Cl]")
        # 13? That includes the 4 hydrogens in [NH4+], but
        # not the implicit hydrogen in c1ccccc1O. I don't get it. XXX
        self.assertEqual(mol.NumAtoms(), 13)
        self.assertEqual(mol.NumBonds(), 11)  # includes the -H bonds
        self.assertEqual(mol.NumHvyAtoms(), 9)

        self.assertClose(mol.GetMolWt(), 147.6027)
        self.assertClose(mol.GetMolWt(True), 147.6027)
        self.assertClose(mol.GetMolWt(False), 141.55506)

        self.assertClose(mol.GetExactMass(), 147.0451)
        self.assertClose(mol.GetExactMass(True), 147.0451)
        self.assertClose(mol.GetExactMass(False), 140.998)

        self.assertEqual(mol.GetTotalCharge(), 1)

#    def test_title(self): # tested in the IO module

    def test_get_bond(self):
        mol = parse_smiles("C=O")
        C = mol.GetAtomById(0)
        self.assertEqual(C.GetAtomicNum(), 6)
        O = mol.GetAtomById(1)
        self.assertEqual(O.GetAtomicNum(), 8)
        bond = mol.GetBond(C, O)
        self.assertEqual(bond.GetBondOrder(), 2)

    def test_formula(self):
        mol = parse_smiles("c1ccccc1O.[NH4+]")
        # XXX Leaves out the "+"?
        self.assertEqual(mol.GetFormula(), "C6H10NO")
        # XXX Why are the extra spaces there? "N  1", "O  1" and the terminal " "
        self.assertEqual(mol.GetSpacedFormula(),  "C 6 H 10 N  1 O  1 ")
        self.assertEqual(mol.GetSpacedFormula(0), "C 6 H 10 N  1 O  1 ")
        self.assertEqual(mol.GetSpacedFormula(1), "C 6 H 10 N O ")
        self.assertEqual(mol.GetSpacedFormula(1, '>'), "C>6>H>10>N>O>")
        # It seems that OpenBabel and I have different definitions of "implicit"
        self.assertEqual(mol.GetSpacedFormula(0, ' ', 0), "C 6 H 4 N  1 O  1 ")
        self.assertEqual(mol.GetSpacedFormula(1, ' ', 0), "C 6 H 4 N O ")

    # There's a huge number of properties I've omitted

class TestAtomAndBond(MyTestCase):
    def test_atom_properties(self):
        mol = parse_smiles("[12CH4-]")
        mol.SetTitle("Spam!")
        atom = mol.GetAtom(0)
        self.assertTrue(atom is None, "GetAtom(0)")
        atom = mol.GetAtom(1)
        self.assertTrue(atom is not None, "GetAtom(1)")

        self.assertEqual(atom.GetAtomicNum(), 6)
        self.assertEqual(atom.GetIsotope(), 12)
        self.assertEqual(atom.GetFormalCharge(), -1)
        self.assertEqual(atom.GetImplicitValence(), 4)
        self.assertEqual(atom.GetIdx(), 1)
        self.assertEqual(atom.GetIndex(), 0)
        self.assertEqual(atom.GetSpinMultiplicity(), 0)
        self.assertEqual(atom.GetAtomicMass(), 12.0)
        self.assertEqual(atom.GetExactMass(), 12.0)
        self.assertEqual(atom.GetValence(), 4)
        self.assertEqual(atom.GetHyb(), 3) # sp3
        self.assertEqual(atom.GetHvyValence(), 0)

        self.assertEqual(atom.GetX(), 0.0)
        self.assertEqual(atom.GetY(), 0.0)
        self.assertEqual(atom.GetZ(), 0.0)
        atom.SetVector(1.25, 2.5, 5.125)

        self.assertEqual(atom.x(), 1.25)
        self.assertEqual(atom.y(), 2.5)
        self.assertEqual(atom.z(), 5.125)

        self.assertEqual(atom.ImplicitHydrogenCount(), 0)
        self.assertEqual(atom.ExplicitHydrogenCount(), 4)
        self.assertEqual(atom.MemberOfRingCount(), 0)
        self.assertEqual(atom.MemberOfRingSize(), 0)
        self.assertEqual(atom.CountRingBonds(), 0)

        # *sigh* I don't like all these silly methods. I would
        # rather they be functions.
        self.assertTrue(atom.IsCarbon())

        self.assertFalse(atom.IsHydrogen())
        self.assertTrue(atom.IsCarbon())
        self.assertFalse(atom.IsNitrogen())
        self.assertFalse(atom.IsOxygen())
        self.assertFalse(atom.IsSulfur())
        self.assertFalse(atom.IsPhosphorus())
        self.assertFalse(atom.IsAromatic())
        self.assertFalse(atom.IsInRing())
        for i in range(10):
            self.assertFalse(atom.IsInRingSize(i))
        self.assertFalse(atom.IsNotCorH())
        self.assertFalse(atom.IsCarboxylOxygen())
        self.assertFalse(atom.IsPhosphateOxygen())
        self.assertFalse(atom.IsSulfateOxygen())
        self.assertFalse(atom.IsNitroOxygen())
        self.assertFalse(atom.IsAmideNitrogen())
        self.assertFalse(atom.IsPolarHydrogen())
        self.assertFalse(atom.IsNonPolarHydrogen())
        self.assertFalse(atom.IsAromaticNOxide())
        self.assertFalse(atom.IsChiral())
        self.assertFalse(atom.IsAxial())
        self.assertFalse(atom.IsHbondAcceptor())
        self.assertFalse(atom.IsHbondDonor())
        self.assertFalse(atom.IsHbondDonorH())

        self.assertFalse(atom.HasBondOfOrder(0))
        self.assertTrue(atom.HasBondOfOrder(1))
        self.assertFalse(atom.HasBondOfOrder(2))
        self.assertFalse(atom.HasBondOfOrder(3))

        self.assertEqual(atom.CountBondsOfOrder(1), 4)

        self.assertFalse(atom.HasNonSingleBond())
        self.assertTrue(atom.HasSingleBond())
        self.assertFalse(atom.HasDoubleBond())
        self.assertFalse(atom.HasAromaticBond())

        # In the 15th or 16th main group (N, O, P, S, ...)
        self.assertFalse(atom.IsHeteroatom())

        # Whee! This isn't really accessible to Python. XXX
        # Should I use ctypes to peer into the object?
        # self.assertEquals(atom.GetCoordinate(), ...?)
        v = atom.GetVector()
        self.assertEqual(v.GetX(), 1.25)
        self.assertEqual(v.GetY(), 2.5)
        self.assertEqual(v.GetZ(), 5.125)

        self.assertClose(atom.GetPartialCharge(), -0.25658)

        self.assertTrue(atom.GetParent().GetTitle() == mol.GetTitle(),
                          "parent is mol")

        self.assertFalse(atom.IsAromatic())
        atom.SetAromatic()
        self.assertTrue(atom.IsAromatic())
        atom.UnsetAromatic()

    def test_more_atom_properties(self):
        mol = parse_smiles("Nc1cc(S)ccc1O")
        self.assertTrue(mol.GetAtom(8).CountFreeOxygens())
        self.assertFalse(mol.GetAtom(9).CountFreeOxygens())
        self.assertTrue(mol.GetAtom(9).ImplicitHydrogenCount())

        atom = mol.GetAtom(2)
        self.assertEqual(atom.MemberOfRingCount(), 1)
        self.assertEqual(atom.MemberOfRingSize(), 6)
        self.assertEqual(atom.CountRingBonds(), 2)

        self.assertEqual(atom.BOSum(), 4)
        self.assertTrue(atom.IsAromatic())
        self.assertTrue(atom.IsInRing())

    def test_bond_length(self):
        mol = parse_smiles("C#N")
        C = mol.GetAtom(1)
        N = mol.GetAtom(2)
        # XXX Why do bonds starts from 0 and not 1
        self.assertTrue(mol.GetBond(1) is None)

        bond = mol.GetBond(0)
        self.assertEqual(bond.GetLength(), 0.0)
        N.SetVector(0.0, 1.0, 0.0)
        self.assertEqual(bond.GetLength(), 1.0)
        length = bond.GetEquibLength()
        bond.SetLength(C, length)

        self.assertEqual(C.GetX(), 0.0)
        self.assertEqual(C.GetY(), 0.0)
        self.assertEqual(C.GetZ(), 0.0)
        self.assertEqual(N.GetX(), 0.0)
        self.assertEqual(N.GetY(), length)
        self.assertEqual(N.GetZ(), 0.0)

    def test_bond_neighbor(self):
        mol = parse_smiles("CNS")
        C = mol.GetAtom(1)
        N = mol.GetAtom(2)
        S = mol.GetAtom(3)
        bond = mol.GetBond(C, N)
        self.assertEqual(bond.GetNbrAtom(C).GetIdx(), N.GetIdx())
        self.assertEqual(bond.GetNbrAtom(N).GetIdx(), C.GetIdx())
        # XXX S isn't part of the bond. The docs need to warn about this behavior
        self.assertEqual(bond.GetNbrAtom(S), C)
        self.assertEqual(bond.GetNbrAtomIdx(C), N.GetIdx())
        self.assertEqual(bond.GetNbrAtomIdx(N), C.GetIdx())
        # This is the documented failure condition
        self.assertEqual(bond.GetNbrAtomIdx(S), bond.GetBeginAtomIdx())

    def test_bond_properties(self):
        mol = parse_smiles("C#N")
        bond = mol.GetBond(0)
        self.assertEqual(bond.GetBondOrder(), 3)
        self.assertFalse(bond.IsDouble())
        self.assertTrue(bond.IsTriple())

        bond.SetBondOrder(2)
        self.assertEqual(bond.GetBondOrder(), 2)

        # It looks like OpenBabel tracks the valences and not the
        # hydrogen counts, which is why this works out right.
        # Interesting.
        smiles = cansmiles(mol)
        self.assertEqual(smiles, "C=N")

        self.assertFalse(bond.IsAromatic())
        self.assertEqual(bond.GetIdx(), 0)
        self.assertEqual(bond.GetId(), 0)
        self.assertEqual(bond.GetBeginAtomIdx(), 1)
        self.assertEqual(bond.GetBeginAtom().GetAtomicNum(), 6)
        self.assertEqual(bond.GetEndAtomIdx(), 2)
        self.assertEqual(bond.GetEndAtom().GetAtomicNum(), 7)

        self.assertFalse(bond.IsAmide())
        self.assertFalse(bond.IsPrimaryAmide())
        self.assertFalse(bond.IsRotor())
        self.assertFalse(bond.IsInRing())
        self.assertFalse(bond.IsSecondaryAmide())
        self.assertFalse(bond.IsTertiaryAmide())
        self.assertFalse(bond.IsEster())
        self.assertFalse(bond.IsCarbonyl())
        self.assertFalse(bond.IsSingle())
        self.assertTrue(bond.IsDouble())
        self.assertFalse(bond.IsTriple())
        self.assertFalse(bond.IsClosure())
        self.assertFalse(bond.IsUp())
        self.assertFalse(bond.IsDown())
        self.assertFalse(bond.IsWedge())
        self.assertFalse(bond.IsHash())
        self.assertFalse(bond.IsWedgeOrHash())
        self.assertFalse(bond.IsCisOrTrans())

        ## This returns True, but the test is rather meaningless
        # since there are no coordinates.
        self.assertTrue(bond.IsDoubleBondGeometry())

    def test_more_bond_properties(self):
        mol = parse_smiles("Sc1nccc1")
        bond = mol.GetBond(2)
        self.assertEqual(bond.GetBeginAtom().GetAtomicNum(), 7)
        self.assertEqual(bond.GetEndAtom().GetIdx(), 4)
        self.assertTrue(bond.IsInRing())

    def test_rings(self):
        mol = parse_smiles("C12CNCC3C1.C2CCC3")
        atom = mol.GetAtom(1)
        self.assertTrue(atom.IsInRing())
        self.assertTrue(atom.IsInRingSize(6))
        self.assertTrue(atom.IsInRingSize(7))
        self.assertFalse(atom.IsInRingSize(10))
        self.assertEqual(atom.MemberOfRingCount(), 2)
        self.assertEqual(atom.MemberOfRingSize(), 6)
        self.assertEqual(atom.CountRingBonds(), 3)

        sssr = mol.GetSSSR()
        self.assertEqual(len(sssr), 2)
        ring_info = [(ring.Size(), ring) for ring in sssr]
        ring_info.sort()

        sizes = [x[0] for x in ring_info]
        self.assertEqual(sizes, [6, 7])

        ring = ring_info[0][1]
        self.assertFalse(ring.IsAromatic())
        self.assertEqual(ring.GetType(), "")
        # XXX *which* of the non-carbons is the root? That isn't documented
        idx = ring.GetRootAtom() # Shouldn't that be "Idx"?
        # Since there's only one non-C, it must be the N
        atom = mol.GetAtom(idx)
        self.assertEqual(atom.GetAtomicNum(), 7)
        self.assertTrue(ring.IsMember(atom))
        for bond in ob.OBAtomBondIter(atom):
            self.assertTrue(ring.IsMember(bond))
        self.assertTrue(ring.IsInRing(idx))


        lssr = mol.GetLSSR()
        self.assertEqual(len(lssr), 2)
        sizes = [ring.Size() for ring in lssr]
        sizes.sort()
        self.assertEqual(sizes, [6, 7])

    def test_ring_center_and_normal(self):
        mol = parse_smiles("c1ccccc1")
        R = 1.5
        for i in range(6):
            atom = mol.GetAtom(i+1)
            atom.SetVector(R*math.cos(2*math.pi*i/6),
                           R*math.sin(2*math.pi*i/6),
                           0.0)
        for ring in ob.OBMolRingIter(mol):
            break

        center = ob.vector3()
        norm1 = ob.vector3()
        norm2 = ob.vector3()
        ring.findCenterAndNormal(center, norm1, norm2)
        self.assertZero(center.GetX())
        self.assertZero(center.GetY())
        self.assertZero(center.GetZ())

        self.assertZero(norm1.GetX())
        self.assertZero(norm1.GetY())
        self.assertClose(norm1.GetZ(), 1.0)

        self.assertZero(norm2.GetX())
        self.assertZero(norm2.GetY())
        self.assertClose(norm2.GetZ(), -1.0)

    def test_geometry_calculations(self):
        mol = parse_smiles("CNOS")
        C = mol.GetAtom(1)
        N = mol.GetAtom(2)
        O = mol.GetAtom(3)
        S = mol.GetAtom(4)
        C.SetVector(0.0, 0.0, 0.0)
        N.SetVector(1.0, 0.0, 0.0)
        O.SetVector(1.5, 0.5, 0.0)
        S.SetVector(1.0, 1.0, 1.0)

        self.assertEqual(C.GetDistance(1), 0.0)
        self.assertEqual(C.GetDistance(N), 1.0)

        # XXX This returns degrees?!
        self.assertClose(C.GetAngle(2, 3), 135.0)
        self.assertClose(C.GetAngle(N, O), 135.0)
        self.assertEqual(C.GetAngle(C, O), 0.0)

        self.assertClose(N.SmallestBondAngle(), 135.0)
        self.assertClose(N.AverageBondAngle(), 135.0)

        self.assertClose(N.SmallestBondAngle(), 135.0)
        self.assertClose(N.AverageBondAngle(), 135.0)

        # The molecule also has an angle method, PLUS torsion
        self.assertClose(mol.GetAngle(C, N, O), 135.0)
        self.assertEqual(mol.GetAngle(C, C, O), 0.0)
        self.assertClose(mol.GetTorsion(C, N, O, S), 54.7356)

        self.assertTrue(C.IsConnected(N))
        self.assertFalse(C.IsConnected(O))

        # XXX I don't expect this
        self.assertTrue(C.IsConnected(C))

        self.assertFalse(C.IsOneThree(S))
        self.assertTrue(N.IsOneThree(S))

        self.assertTrue(C.IsOneFour(S))
        # XXX I don't expect this.
        # I think it's a consequence of X.IsConnected(X) == True
        self.assertTrue(C.IsOneFour(O))
        self.assertTrue(C.IsOneFour(C)) # XXX completely suprising!

    def test_HtoMethyl(self):
        mol = parse_smiles("[H]Cl")
        # If I don't move this atom then I get the message
        #  *** Open Babel Warning  in SetLength
        #    Atoms are both at the same location, moving out of the way.
        mol.GetAtom(2).SetVector(1.5, 0, 0)

        atom = mol.GetAtom(1)

        self.assertEqual(atom.GetAtomicNum(), 1)
        # This triggers some debug code which dumps to cerr
        atom.HtoMethyl()
        self.assertEqual(atom.GetAtomicNum(), 6)

    def test_MatchesSMARTS(self):
        # I don't much like this function.
        mol = parse_smiles("CCO")
        atom = mol.GetAtom(1)
        self.assertEqual(atom.MatchesSMARTS("O"), 0)
        self.assertEqual(atom.MatchesSMARTS("OCC"), 0)
        self.assertEqual(atom.MatchesSMARTS("CC"), 1)
    # HasAlphaBetaUnsat



# These values are taken directly from OB's spectrophoretest.cpp
class SpectorphoreTest(MyTestCase):
    def assertWithin_0_001(self, val, expect):
        assert val > 0
        self.assertTrue(abs(val - expect) < 0.001, val)
    def _make_mol(self):
        mol = ob.OBMol()
        def new_atom(eleno):
            a = mol.NewAtom()
            a.SetAtomicNum(eleno)
            return a
        atoms = []
        atoms.append(new_atom(6))
        atoms.append(new_atom(1))
        atoms.append(new_atom(9))
        atoms.append(new_atom(35))
        atoms.append(new_atom(17))
        for atom in atoms[1:]:
            b = mol.NewBond()
            b.SetBegin(atoms[0])
            b.SetEnd(atom)
            b.SetBondOrder(1)

        mol.GetAtom(1).SetVector(-0.013, 1.086, 0.008)
        mol.GetAtom(2).SetVector(0.002, -0.004, 0.002)
        mol.GetAtom(3).SetVector(1.300, 1.570, -0.002)
        mol.GetAtom(4).SetVector(-0.964, 1.737, -1.585)
        mol.GetAtom(5).SetVector(-0.857, 1.667, 1.491)
        return mol

    def test_1(self):
        s = ob.OBSpectrophore()
        s.SetNormalization(ob.OBSpectrophore.NoNormalization)
        s.SetResolution(3.0)
        s.SetAccuracy(ob.OBSpectrophore.AngStepSize20)
        s.SetStereo(ob.OBSpectrophore.NoStereoSpecificProbes)
        r = s.GetSpectrophore(self._make_mol())

        C = self.assertWithin_0_001
        C(r[ 0],  1.599)
        C(r[ 1],  1.577)
        C(r[ 2],  1.170)
        C(r[ 3],  3.761)
        C(r[ 4],  5.175)
        C(r[ 5],  5.781)
        C(r[ 6],  3.797)
        C(r[ 7],  3.713)
        C(r[ 8],  4.651)
        C(r[ 9],  7.737)
        C(r[10],  7.950)
        C(r[11],  4.869)
        C(r[12],  2.708)
        C(r[13],  3.471)
        C(r[14],  6.698)  # XXX The original code has a bug in the upper bound
        C(r[15],  9.486)
        C(r[16],  7.668)
        C(r[17],  8.882)
        C(r[18],  4.900)  # XXX The original code is slightly too high in the upper bound
        C(r[19],  7.479)
        C(r[20],  9.324)
        C(r[21], 10.293)
        C(r[22], 12.956)
        C(r[23], 10.335)
        C(r[24],  4.021)
        C(r[25],  3.814)
        C(r[26],  2.947)
        C(r[27],  6.381)
        C(r[28], 11.004)
        C(r[29],  8.279)
        C(r[30],  6.549)
        C(r[31],  7.136)
        C(r[32],  8.613)
        C(r[33], 13.182)
        C(r[34], 13.744)
        C(r[35],  9.084)
        C(r[36],  0.459)
        C(r[37],  0.642)
        C(r[38],  2.172)
        C(r[39],  2.753)
        C(r[40],  2.348)
        C(r[41],  2.605)
        C(r[42],  1.614)
        C(r[43],  3.166)
        C(r[44],  3.391)
        C(r[45],  3.132)
        C(r[46],  4.105)
        C(r[47],  2.875)

    def test_with_increased_accuracy(self):
        s = ob.OBSpectrophore()
        s.SetNormalization(ob.OBSpectrophore.NoNormalization)
        s.SetResolution(3.0)
        s.SetAccuracy(ob.OBSpectrophore.AngStepSize5)
        s.SetStereo(ob.OBSpectrophore.NoStereoSpecificProbes)
        r = s.GetSpectrophore(self._make_mol())
        C = self.assertWithin_0_001
        C(r[0], 1.6445)
        C(r[12], 2.7245)
        C(r[24], 4.0435)
        C(r[36], 0.4585)

    # Look at spectrophoretest.cpp for many more examples.

class TestForceFields(MyTestCase):
    def test_plugin(self):
        v = ob.vectorString()
        ob.OBPlugin.ListAsVector("forcefields", None, v)
        # Huh. The plugin system uses case-insensitive lookup
        names = [x.split()[0].lower() for x in v]

        self.assertTrue("gaff" in names, names)
        self.assertTrue("mmff94" in names, names)
        self.assertTrue("uff" in names, names)

        pFF1 = ob.OBForceField.FindForceField("GAFF")
        pFF2 = ob.OBForceField.FindForceField("GafF")
        self.assertFalse(pFF1 is None)
        self.assertFalse(pFF2 is None)
        self.assertEqual(pFF1.GetID(), pFF2.GetID())

    def _test_energies(self, plugin_name, expected_results, filename = None):
        pFF = ob.OBForceField.FindForceField(plugin_name)
        self.assertFalse(pFF is None, "Cannot load " + plugin_name)

        if filename is None:
            filename = testfile("forcefield.sdf")

        for i, mol in enumerate(readfile(filename, "sdf")):
            self.assertEqual(pFF.Setup(mol), 1,
                              "Could not set up forcefield on " + mol.GetTitle())
            energy = pFF.Energy(False)

            self.assertClose(energy, expected_results[i])
            self.assertEqual(pFF.ValidateGradients(), 1,
                              "gradients do not validate for molecule " + mol.GetTitle())

            # These are meant to be fast unit tests, and not a validation suite.
            # Therefore I'll only run two tests then exit
            if i == 0:
                break

    # The basis for these tests come from ffghemical.cpp
    def test_ghemical_energy_calculation(self):
        expected_results = list(map(float, open(testfile("ghemicalresults.txt")).readlines()))
        self._test_energies("Ghemical", expected_results)

    # The basis for these tests come from ffgaff.cpp
    def test_gaff_energy_calculation(self):
        expected_results = list(map(float, open(testfile("gaffresults.txt")).readlines()))
        self._test_energies("GAFF", expected_results)

    # These basis for these tests comes from ffmmff94.cpp
    def test_mmff94_energy_calculation(self):
        # XXX The MMFF94 ValidateGradients() dumps output to stdout
        expected_results = list(map(float, open(testfile("mmff94results.txt")).readlines()))
        self._test_energies("MMFF94", expected_results)

        ## Doing this test does not show anything new about the OpenBabel API
        ## and it dumps more useless text (for the purposes of testing) to stdout
        #expected_results = map(float, open(testfile("more-mmff94results.txt")).readlines())
        #self._test_energies("MMFF94", expected_results, testfile("more-mmff94.sdf"))

    # These basis for these tests comes from ffuff.cpp
    def test_uff_energy_calculation(self):
        expected_results = list(map(float, open(testfile("uffresults.txt")).readlines()))
        self._test_energies("UFF", expected_results)


# Does not seem to work. Don't know if I'm doing the wrong thing
#    def test_mmff94_validates(self):
#        pFF = ob.OBForceField.FindForceField("MMFF94")
#        self.assertEquals(pFF.Validate(), True)

if __name__ == "__main__":
    unittest.main()

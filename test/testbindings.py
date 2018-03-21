"""Test OpenBabel Python bindings

On Windows or Linux, you can run these tests at the commandline
in the build folder with:
"C:\Program Files\CMake 2.6\bin\ctest.exe" -C CTestTestfile.cmake
                                           -R pybindtest -VV

The runtime directory is ${CMAKE_SRC_DIR}/test. 

You could also "chdir" into build and run the test file directly:
python ../../test/testbindings.py

In this latter case, you will need to set the environment variables
PYTHONPATH, LD_LIBRARY_PATH, BABEL_LIBDIR and BABEL_DATADIR beforehand.
The CMake script does this automatically.

In both cases, the test file is run directly from the source folder,
and so you can quickly develop the tests and try them out.
"""

import os
import re
import sys
import unittest

here = sys.path[0]
iswin = sys.platform.startswith("win")

try:
    import openbabel as ob
except ImportError:
    ob = None

try:
    import pybel
except ImportError:
    pybel = None

class PythonBindings(unittest.TestCase):
    def setUp(self):
        self.assertTrue(ob is not None, "Failed to import the openbabel module")

class TestPythonBindings(PythonBindings):
    def testSimple(self):
        mol = ob.OBMol()
        conv = ob.OBConversion()
        conv.SetInFormat("smi")
        conv.ReadString(mol, "CC(=O)Cl")
        self.assertAlmostEqual(mol.GetMolWt(), 78.5, 1)
    
class PybelWrapper(PythonBindings):
    def testDummy(self):
        self.assertTrue(pybel is not None, "Failed to import the Pybel module")

class TestSuite(PythonBindings):

    def testInChIIsotopes(self):
        """Ensure that we correctly set and read isotopes in InChIs"""
        with open(os.path.join(here, "inchi", "inchi_isotopes.txt")) as inp:
            for line in inp:
                if line.startswith("#"): continue
                smi, inchi = line.rstrip().split("\t")
                minchi = pybel.readstring("smi", smi).write("inchi").rstrip()
                self.assertEqual(minchi, inchi)
                msmi = pybel.readstring("inchi", minchi).write("smi").rstrip()
                self.assertEqual(msmi, smi)

    def testAsterisk(self):
        """Ensure that asterisk in SMILES is bracketed when needed
        and not otherwise"""
        smis = [ # these don't need brackets for *
                "*", "*C", "C*C",
                # these do need brackets for *
                "[*+]", "[*-]", "[*H]",
                "*[H]", # this one is written as [*H]
                "[1*]", "[*@@](F)(Cl)(Br)I", "[*:1]"]

        for smi in smis:
            needsbracket = "[" in smi
            roundtrip = pybel.readstring("smi", smi).write("smi", opt={"a":True})
            hasbracket = "[" in roundtrip
            self.assertEqual(hasbracket, needsbracket)

    def testSmilesAtomOrder(self):
        """Ensure that SMILES atom order is written correctly"""
        data = [("CC", "1 2"),
                ("O=CCl", "3 2 1")]
        for smi, atomorder in data:
            mol = pybel.readstring("smi", smi)
            mol.write("can", opt={"O": True})
            res = mol.data["SMILES Atom Order"]
            self.assertEqual(res, atomorder)
        mol = pybel.readstring("smi", "CC")
        mol.write("can")
        self.assertFalse("SMILES Atom Order" in mol.data)

    def testOBMolSeparatePreservesAromaticity(self):
        """If the original molecule had aromaticity perceived,
        then the fragments should also.
        """
        smi = "C.c1ccccc1"
        # Two passes: One with aromaticity perceived on the orig mol and
        #             one without
        for N in range(2):
            obmol = pybel.readstring("smi", smi).OBMol
            # Aromaticity is perceived during the last step of reading SMILES
            # so let's unset it here for the first pass
            if N == 0:
                obmol.UnsetAromaticPerceived()
            else:
                self.assertTrue(obmol.HasAromaticPerceived())

            # After separation, is aromaticity the same as the parent?
            mols = obmol.Separate()
            if N == 0:
                self.assertFalse(mols[1].HasAromaticPerceived())
            else:
                self.assertTrue(mols[1].HasAromaticPerceived())

            atom = mols[1].GetAtom(1)
            atom.SetImplicitHCount(0) # mess up the structure
            if N == 0:
                self.assertFalse(atom.IsAromatic())
            else:
                self.assertTrue(atom.IsAromatic())

    def testAtomMapsAfterCopying(self):
        """Copying a molecule should copy the atom maps"""
        smi = "C[CH2:2]O[Cl:6]"
        obmol = pybel.readstring("smi", smi).OBMol
        copy = pybel.ob.OBMol(obmol)
        copysmi = pybel.Molecule(copy).write("smi", opt={"a": True})
        self.assertEqual(copysmi.rstrip(), smi)

    def testAtomMapsAfterAddition(self):
        """Adding two molecules should not mess up the atom maps"""
        data = [
                ("C", "[OH2:2]"),
                ("[OH2:2]", "C"),
                ("[CH4:6]", "[OH2:2]"),
                ("C", "O"),
               ]
        for a, b in data:
            mols = [pybel.readstring("smi", x) for x in [a, b]]
            mols[0].OBMol += mols[1].OBMol
            joined = mols[0].write("smi", opt={"a":True}).rstrip()
            self.assertEqual(joined, "%s.%s" % (a, b))

    def testAtomMapsAfterDeletion(self):
        """Removing atoms/hydrogens should not mess up the atom maps"""
        smis = ["C[NH2:2]", "[CH3:1][NH2:2]"]
        for smi in smis:
            mol = pybel.readstring("smi", smi)
            mol.OBMol.DeleteAtom(mol.OBMol.GetAtom(1))
            self.assertEqual(mol.write("smi", opt={"a":True}).rstrip(), "[NH2:2]")
        smi = "[H]C[NH:2]"
        mol = pybel.readstring("smi", smi)
        mol.removeh()
        self.assertEqual(mol.write("smi", opt={"a":True}).rstrip(), "C[NH:2]")

    def testSquarePlanar(self):
        """Tighten up the parsing of SP stereochemistry in SMILES"""
        good = [
                "C[S@SP1](Cl)(Br)I",
                "C[S@SP2](Cl)(Br)I",
                "C[S@SP3](Cl)(Br)I",
                ]
        bad = [ # raises error
                "C[S@SP0](Cl)(Br)I",
                "C[S@SP4](Cl)(Br)I",
                "C[S@@SP1](Cl)(Br)I",
                "C[S@SP11](Cl)(Br)I",
                "C[S@SO1](Cl)(Br)I",
              ]
        alsobad = [ # just a warning
                "C[S@SP1](Cl)(Br)(F)I",
                "C[S@SP1](Cl)(Br)(F)1CCCC1",
                ]
        for smi in good:
            mol = pybel.readstring("smi", smi)
            self.assertTrue(mol.OBMol.GetData(ob.StereoData))
        for smi in bad:
            self.assertRaises(IOError, pybel.readstring, "smi", smi)
        for smi in alsobad:
            mol = pybel.readstring("smi", smi)
            self.assertTrue(mol.OBMol.GetData(ob.StereoData))

    def testFuzzingTestCases(self):
        """Ensure that fuzzing testcases do not cause crashes"""

        # rejected as invalid smiles
        smis = [r"\0", "&0", "=&",
                "[H][S][S][S@S00]0[S][S@S00H](0[S@S00][S])0n"]
        for smi in smis:
            self.assertRaises(IOError, pybel.readstring, "smi", smi)

        smis = ["c0C[C@H](B)00O0"] # warning and stereo ignored
        for smi in smis:
            pybel.readstring("smi", smi)

    def testImplicitCisDblBond(self):
        """Ensure that dbl bonds in rings of size 8 or less are always
        implicitly cis"""
        smi = "C1/C=C/C"
        for i in range(5): # from size 4 to 8
            ringsize = i + 4
            ringsmi = smi + "1"
            roundtrip = pybel.readstring("smi", ringsmi).write("smi")
            self.assertTrue("/" not in roundtrip)
            smi += "C"
        ringsize = 9
        ringsmi = smi + "1"
        roundtrip = pybel.readstring("smi", ringsmi).write("smi")
        self.assertTrue("/" in roundtrip)

    def testKekulizationOfHypervalents(self):
        # We should support hypervalent aromatic S and N (the latter
        # as we write them)
        data = [("Cs1(=O)ccccn1",
                 "CS1(=O)=NC=CC=C1"),
                ("n1c2-c(c3cccc4cccc2c34)n(=N)c2ccccc12",
                 "n1c2-c(c3cccc4cccc2c34)n(=N)c2ccccc12")]
        for inp, out in data:
            mol = pybel.readstring("smi", inp)
            self.assertEqual(out, mol.write("smi").rstrip())

    def testKekulizationOfcn(self):
        """We were previously not reading 'cn' correctly, or at least how
        Daylight would"""
        mol = pybel.readstring("smi", "cn")
        self.assertEqual("C=N", mol.write("smi").rstrip())

    def testReadingMassDifferenceInMolfiles(self):
        """Previously we were rounding incorrectly when reading the mass diff"""
        template = """
 OpenBabel02181811152D

  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 %2s %2d  0  0  0  0  0  0  0  0  0  0  0
M  END
"""
        # Positive test cases:
        # These are the BIOVIA Draw answers for the first 50 elements for
        # a mass diff of 1
        answers = [2,5,8,10,12,13,15,17,20,21,24,25,28,29,32,33,36,41,40,41,46,49,52,53,56,57,60,60,65,66,71,74,76,80,81,85,86,89,90,92,94,97,99,102,104,107,109,113,116,120,123]
        for idx, answer in enumerate(answers):
            elem = idx + 1
            molfile = template % (ob.GetSymbol(elem), 1)
            mol = pybel.readstring("mol", molfile).OBMol
            iso = mol.GetAtom(1).GetIsotope()
            self.assertEqual(answer, iso)

        # Also test D and T - BIOVIA Draw ignores the mass diff
        for elem, answer in zip("DT", [2, 3]):
            molfile = template % (elem, 1)
            mol = pybel.readstring("mol", molfile).OBMol
            iso = mol.GetAtom(1).GetIsotope()
            self.assertEqual(answer, iso)

        # Negative test cases:
        # Test error message for out-of-range values
        for value in [5, -4]:
            molfile = template % ("C", value)
            mol = pybel.readstring("mol", molfile).OBMol
            iso = mol.GetAtom(1).GetIsotope()
            self.assertEqual(0, iso)

    def testInvalidCharsInSmiles(self):
        """Check that inserting a comma in a SMILES string in various positions
        does not result in a valid SMILES"""
        data = "C, ,C [C,] [,C] [1,C] [C:,1] [C:1,] [CH,] [C+,] [C++,]"
        data += " [C@,](Br)(Cl)(I)F C1,CC1 C%1,1CC%11 C%(1,1)CC%11"
        data += " C=,C"
        for smi in data.split(" "):
            self.assertRaises(IOError, pybel.readstring, "smi", smi)

    def testOBMolAssignTotalChargeToAtoms(self):
        """Run the test cases described in the source code"""
        data = [("[NH4]", +1, "[NH4+]"),
                ("CC(=O)[O]", -1, "CC(=O)[O-]"),
                ("C[CH2]", +1, "C[CH2+]"),
                ("C[CH2]", -1, "C[CH2-]"),
                ("[NH3]CC(=O)[O]", 0, "[NH3+]CC(=O)[O-]"),
                ("S(=O)(=O)([O])[O]", -2, "S(=O)(=O)([O-])[O-]"),
                ("[NH4].[Cl]", 0, "[NH4+].[Cl-]")]
        for smi, charge, ans in data:
            mol = pybel.readstring("smi", smi)
            mol.OBMol.AssignTotalChargeToAtoms(charge)
            self.assertEqual(mol.write("smi").rstrip(), ans)

    def testReadingBenzyne(self):
        """Check that benzyne is read correctly"""
        smi = "c1cccc#c1"
        mol = pybel.readstring("smi", smi)
        self.assertEqual("C1=CC=CC#C1", mol.write("smi").rstrip())

    def testSmilesParsingOfAllElements(self):
        smi = "*[H][He][Li][Be][B][C][N][O][F][Ne][Na][Mg][Al][Si][P][S][Cl][Ar][K][Ca][Sc][Ti][V][Cr][Mn][Fe][Co][Ni][Cu][Zn][Ga][Ge][As][Se][Br][Kr][Rb][Sr][Y][Zr][Nb][Mo][Tc][Ru][Rh][Pd][Ag][Cd][In][Sn][Sb][Te][I][Xe][Cs][Ba][La][Ce][Pr][Nd][Pm][Sm][Eu][Gd][Tb][Dy][Ho][Er][Tm][Yb][Lu][Hf][Ta][W][Re][Os][Ir][Pt][Au][Hg][Tl][Pb][Bi][Po][At][Rn][Fr][Ra][Ac][Th][Pa][U][Np][Pu][Am][Cm][Bk][Cf][Es][Fm][Md][No][Lr][Rf][Db][Sg][Bh][Hs][Mt][Ds][Rg][Cn][Nh][Fl][Mc][Lv][Ts][Og]"
        roundtrip = pybel.readstring("smi", smi).write("smi").rstrip()
        self.assertEqual(roundtrip, smi.replace("[O]", "O").replace("[S]", "S"))
        # aromatic
        # - like C
        elems = ["C", "Si", "Ge", "Sn"]
        for elem in elems:
            mol = pybel.readstring("smi", "c1ccc[%sH]c1" % elem.lower(), opt={"a": True})
            self.assertNotEqual(0, mol.OBMol.NumAtoms())
        # - like N
        elems = ["B", "N", "P", "As", "Sb", "Bi"]
        for elem in elems:
            mol = pybel.readstring("smi", "c1cc[%sH]c1" % elem.lower())
            self.assertNotEqual(0, mol.OBMol.NumAtoms())
        # - like O
        elems = ["O", "S", "Se", "Te"]
        for elem in elems:
            mol = pybel.readstring("smi", "c1cc[%s]c1" % elem.lower())
            self.assertNotEqual(0, mol.OBMol.NumAtoms())
        # Organic subset aromatics
        # - like N
        elems = ["B", "N", "P"]
        for elem in elems:
            mol = pybel.readstring("smi", "c1cc%scc1" % elem.lower())
            self.assertNotEqual(0, mol.OBMol.NumAtoms())
        # - like O
        elems = ["O", "S"]
        for elem in elems:
            mol = pybel.readstring("smi", "c1cc%sc1" % elem.lower())
            self.assertNotEqual(0, mol.OBMol.NumAtoms())

    def testSmilesParsingAndWritingOfLargeIsotopes(self):
        smis = ["[1C]", "[11C]", "[111C]", "[1111C]"]
        for smi in smis:
            mol = pybel.readstring("smi", smi)
            self.assertEqual(mol.write("smi").rstrip(), smi)
        self.assertRaises(IOError, pybel.readstring, "smi", "[11111C]")
        mol = pybel.readstring("smi", "[C]")
        mol.atoms[0].OBAtom.SetIsotope(65535)
        self.assertEqual(mol.write("smi").rstrip(), "[C]")

    def testOBMolSeparatePreservesAtomOrder(self):
        """Originally Separate() preserved DFS order rather
        than atom order"""
        # First test
        smi = "C123.F3.Cl2.Br1"
        mol = pybel.readstring("smi", smi)
        atomicnums = [atom.OBAtom.GetAtomicNum() for atom in mol]
        mols = mol.OBMol.Separate()
        new_atomicnums = [atom.OBAtom.GetAtomicNum() for atom in pybel.Molecule(mols[0])]
        for x, y in zip(atomicnums, new_atomicnums):
            self.assertEqual(x, y) # check that the atoms have not been permuted
        # Second test
        xyz = """6
examples/water_dimer.xyz
O          0.12908       -0.26336        0.64798
H          0.89795        0.28805        0.85518
H          0.10833       -0.20468       -0.33302
O          0.31020        0.07569       -2.07524
H          0.64083       -0.57862       -2.71449
H         -0.26065        0.64232       -2.62218
"""
        mol = pybel.readstring("xyz", xyz)
        mols = mol.OBMol.Separate()
        allatoms = pybel.Molecule(mols[0]).atoms + pybel.Molecule(mols[1]).atoms
        for idx, atom in enumerate(allatoms):
            xcoord = atom.OBAtom.GetX()
            orig_xcoord = mol.OBMol.GetAtom(idx+1).GetX()
            self.assertEqual(xcoord, orig_xcoord)

    def testStereoRefsAfterAddingOBMols(self):
        """The stereo ref for an implicit H ref was being set to 0"""
        smis = ["C", "[C@@H](Br)(Cl)I"]
        mols = [pybel.readstring("smi", smi) for smi in smis]
        # FIXME - does not seem to be possible to work out whether
        # tetrahedral or not from Python?
        stereodata = mols[1].OBMol.GetData(ob.StereoData)
        config = ob.toTetrahedralStereo(stereodata).GetConfig()
        self.assertEqual(config.from_or_towards, 4294967294)
        mols[0].OBMol += mols[1].OBMol
        self.assertEqual(mols[0].write("smi").rstrip(), ".".join(smis))
        stereodata = mols[0].OBMol.GetData(ob.StereoData)
        config = ob.toTetrahedralStereo(stereodata).GetConfig()
        self.assertEqual(config.from_or_towards, 4294967294)

    def testWhetherAllElementsAreSupported(self):
        """Check whether a new element has been correctly added"""
        N = 0
        while ob.GetSymbol(N+1):
            N += 1
            # Is the symbol parsed?
            symbol = ob.GetSymbol(N)
            self.assertEqual(N, ob.GetAtomicNum(symbol))
            # Has an exact mass been set?
            self.assertNotEqual(0.0, ob.GetExactMass(N))
            # Has the symbol been added to the SMILES parser?
            numatoms = pybel.readstring("smi", "[%s]" % symbol).OBMol.NumAtoms()
            self.assertEqual(numatoms, 1)
            # Check whether the element is available as a constant
            self.assertEqual(N, getattr(ob, ob.GetName(N)))

        self.assertTrue(N > 100)

    def testIterators(self):
        """Basic check that at least two iterators are working"""
        mol = pybel.readstring("smi", "c1ccccc1C(=O)Cl")
        atoms = list(ob.OBMolAtomIter(mol.OBMol))
        self.assertEqual(len(atoms), 9)
        elements = [atom.GetAtomicNum() for atom in atoms]
        self.assertEqual(elements, [6,6,6,6,6,6,6,8,17])
        bonds = list(ob.OBMolBondIter(mol.OBMol))
        self.assertEqual(len(bonds), 9)

class Radicals(PythonBindings):
    def testSmilesToMol(self):
        smis = ["C", "[CH3]", "[CH2]", "[CH2]C", "[C]"]
        valences = [0, 3, 2, 3, 15]
        for smi, valence in zip(smis, valences):
            mol = pybel.readstring("smi", smi)
            molfile = mol.write("mol")
            firstcarbon = molfile.split("\n")[4]
            mvalence = int(firstcarbon[48:53])
            self.assertEqual(valence, mvalence)
            # test molfile->smiles
            msmi = pybel.readstring("mol", molfile).write("smi").rstrip()
            self.assertEqual(smi, msmi)
    def testMolfileInvalidValencyToSmiles(self):
        """Should trigger warning about invalid valency field
        Unfortunately - there's no way to test this from Python as it requires
        streams :-/"""
        molfile = """
 OpenBabel07291719012D

  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0 15  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  4  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  END
"""
        mol = pybel.readstring("mol", molfile)
        self.assertEqual("[C]C", mol.write("smi").rstrip())
    def testSettingSpinMult(self):
        """Set spin and read/write it"""
        mol = pybel.readstring("smi", "C")
        mol.atoms[0].OBAtom.SetSpinMultiplicity(2)
        molfile = mol.write("mol")
        self.assertEqual("M  RAD  1   1   2", molfile.split("\n")[5])
        molb = pybel.readstring("mol", molfile)
        self.assertEqual(2, molb.atoms[0].OBAtom.GetSpinMultiplicity())
        self.assertEqual(4, molb.atoms[0].OBAtom.GetImplicitHCount())

class AcceptStereoAsGiven(PythonBindings):
    """Tests to ensure that stereo is accepted as given in input from
    SMILES and Mol files"""
    def testSmiToSmi(self):
        # Should preserve stereo
        tet = "[C@@H](Br)(Br)Br"
        out = pybel.readstring("smi", tet).write("smi")
        self.assertTrue("@" in out)
        cistrans = r"C/C=C(\C)/C"
        out = pybel.readstring("smi", cistrans).write("smi")
        self.assertTrue("/" in out)
        # Should wipe stereo
        out = pybel.readstring("smi", tet, opt={"S": True}).write("smi")
        self.assertFalse("@" in out)
        cistrans = r"C/C=C(\C)/C"
        out = pybel.readstring("smi", cistrans, opt={"S": True}).write("smi")
        self.assertFalse("/" in out)

class AtomClass(PythonBindings):
    """Tests to ensure that refactoring the atom class handling retains
    functionality"""

    def testSMILES(self):
        mol = pybel.readstring("smi", "C[CH3:6]")
        atom = mol.OBMol.GetAtom(2)
        data = atom.GetData("Atom Class")
        self.assertTrue(data)
        self.assertEqual(6, ob.toPairInteger(data).GetGenericValue())

        atom.DeleteData("Atom Class")
        ac = ob.obpairtemplateint()
        ac.SetAttribute("Atom Class")
        ac.SetValue(2)
        mol.OBMol.GetAtom(1).CloneData(ac)
        out = mol.write("smi", opt={"a":True, "n":True, "nonewline":True})
        self.assertEqual("[CH3:2]C", out)

    def testMOL(self):
        """Roundtrip thru MOL file"""
        smi = "C[CH3:6]"
        mol = pybel.readstring("smi", smi)
        molfile = mol.write("mol", opt={"a":True})
        molb = pybel.readstring("mol", molfile)
        out = mol.write("smi", opt={"a":True, "n":True, "nonewline":True})
        self.assertEqual(smi, out)

    def testRGroup(self):
        """[*:1] is converted to R1 in MOL file handling"""
        smi = "[*:6]C"
        mol = pybel.readstring("smi", smi)
        molfile = mol.write("mol")
        self.assertTrue("M  RGP  1   1   6" in molfile)
        molb = pybel.readstring("mol", molfile)
        out = mol.write("smi", opt={"a":True, "n":True, "nonewline":True})
        self.assertEqual(smi, out)

    def testCML(self):
        """OB stores atom classes using _NN at the end of atom ids"""
        smis = ["[CH3:6]C", "[CH3:6][OH:6]",
                "O"+"[CH2:2]"*27+"O"
                ]
        for smi in smis:
            mol = pybel.readstring("smi", smi)
            cml = mol.write("cml")
            molb = pybel.readstring("mol", cml)
            out = mol.write("smi", opt={"a":True, "n":True, "nonewline":True})
            self.assertEqual(smi, out)

    def testTinkerXYZ(self):
        """Atom classes are written out as the atom types (though
        not currently read)"""
        smi = "[CH4:23]"
        mol = pybel.readstring("smi", smi)
        xyz = mol.write("txyz", opt={"c": True})
        lines = xyz.split("\n")
        broken = lines[1].split()
        self.assertEqual("23", broken[-1].rstrip())

    def testDeleteHydrogens(self):
        """Don't suppress a hydrogen with an atom class"""
        smi = "C([H])([H])([H])[H:1]"
        mol = pybel.readstring("smi", smi)
        mol.OBMol.DeleteHydrogens()
        nsmi = mol.write("smi", opt={"a": True, "h": True})
        self.assertEqual("C[H:1]", nsmi.rstrip())

if __name__ == "__main__":
    unittest.main()

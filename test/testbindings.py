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
import sys
import unittest
import itertools

here = sys.path[0]
iswin = sys.platform.startswith("win")

try:
    from openbabel import openbabel as ob
except ImportError:
    ob = None

try:
    from openbabel import pybel
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

    def testAromaticityPreservedOnAtomDeletion(self):
        """Ensure that aromaticity is preserved on atom deleteion"""
        mol = pybel.readstring("smi", "c1ccccc1").OBMol
        mol.DeleteAtom(mol.GetFirstAtom())
        self.assertTrue(mol.GetFirstAtom().IsAromatic())
        mol.SetAromaticPerceived(False)
        self.assertFalse(mol.GetFirstAtom().IsAromatic())

    def testLPStereo(self):
        """Ensure that nitrogen and sulfur can support LP stereo"""
        data = ["[N@@](Cl)(Br)I", "Cl[N@@](Br)I",
                "[S@@](Cl)(Br)I", "Cl[S@@](Br)I"]
        for smi in data:
            mol = pybel.readstring("smi", smi)
            self.assertTrue(mol.OBMol.GetData(ob.StereoData))
            nsmi = mol.write("smi").rstrip()
            self.assertEqual(smi, nsmi)

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

    def testECFP(self):
        data = [
                ("CC", 1, 2),
                ("CCC", 2, 4),
                ("CC(C)C", 2, 4),
                ("CC(C)(C)C", 2, 4),
                ]
        for smi, numA, numB in data:
            mol = pybel.readstring("smi", smi)
            ecfp0 = mol.calcfp("ecfp0").bits
            self.assertEqual(len(ecfp0), numA)
            ecfp2 = mol.calcfp("ecfp2").bits
            self.assertEqual(len(ecfp2), numB)
            for bit in ecfp0:
                self.assertTrue(bit in ecfp2)

    def testOldRingInformationIsWipedOnReperception(self):
        """Previously, the code that identified ring atoms and bonds
        did not set the flags of non-ring atoms. This meant that no
        matter what you did to the structure, once a ring-atom, always a
        ring atom."""
        mol = pybel.readstring("smi", "c1ccccc1")
        atom = mol.atoms[0].OBAtom
        self.assertTrue(atom.IsInRing()) # trigger perception
        mol.OBMol.DeleteAtom(mol.atoms[-1].OBAtom)
        self.assertFalse(atom.IsInRing()) # this used to return True

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
                obmol.SetAromaticPerceived(False)
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

    def testFFGradients(self):
        """Support public access of FF gradients"""
        xyz = """3
water
O          1.02585       -0.07579        0.08189
H          1.99374       -0.04667        0.04572
H          0.74700        0.50628       -0.64089
"""
        mol = pybel.readstring("xyz", xyz)
        ff = pybel._forcefields["mmff94"]

        self.assertTrue(ff.Setup(mol.OBMol))
        energy = ff.Energy() # note: calling GetGradient w/o calling Energy()
                             #       just returns random numbers
        for atom in mol.atoms:
            # this should throw an AttributeError if not available
            grad = ff.GetGradient(atom.OBAtom)
            self.assertNotEqual(0.0, grad.GetX())
            self.assertNotEqual(0.0, grad.GetY())
            self.assertNotEqual(0.0, grad.GetZ())

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

    def testElementsSpecifiedByAtomicNumberInSmiles(self):
        smis = [
                ("[#100]", "[Fm]"),
                ("[#255]", None),
                ("[254#255]", None),
                ("[#6]", -1),
                ("[#256]", -1),
                ("[#10a]", -1)
               ]

        for smi, rt in smis:
            if rt==-1:
                self.assertRaises(IOError, pybel.readstring, "smi", smi)
                continue
            if rt is None:
                rt = smi
            nsmi = pybel.readstring("smi", smi).write("smi").rstrip()
            self.assertEqual(rt, nsmi)

    def testIterators(self):
        """Basic check that at least two iterators are working"""
        mol = pybel.readstring("smi", "c1ccccc1C(=O)Cl")
        atoms = list(ob.OBMolAtomIter(mol.OBMol))
        self.assertEqual(len(atoms), 9)
        elements = [atom.GetAtomicNum() for atom in atoms]
        self.assertEqual(elements, [6,6,6,6,6,6,6,8,17])
        bonds = list(ob.OBMolBondIter(mol.OBMol))
        self.assertEqual(len(bonds), 9)

    def testProper2DofFragments(self):
        """Check for proper handling of fragments in mcdl routines, see issue #1889"""
        mol = pybel.readstring("smi", "[H+].CC[O-].CC[O-]")
        mol.draw(show=False, update=True)
        dists = [
            abs(a.coords[0] - b.coords[0]) + abs(a.coords[1] - b.coords[1])
            for a, b in itertools.combinations(mol.atoms, 2)
        ]
        mindist = min(dists)
        self.assertTrue(mindist > 0.00001)

    def testRegressionBenzene2D(self):
        """Check that benzene is given a correct layout, see #1900"""
        mol = pybel.readstring("smi", "c1ccccc1")
        mol.draw(show=False, update=True)
        benzmol = """
 OpenBabel10161813072D

  6  6  0  0  0  0  0  0  0  0999 V2000
   -0.8660   -0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7321   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7321    1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8660    1.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0000    1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  6  2  0  0  0  0
  1  2  1  0  0  0  0
  2  3  2  0  0  0  0
  3  4  1  0  0  0  0
  4  5  2  0  0  0  0
  5  6  1  0  0  0  0
M  END
"""
        self.assertEqual(mol.write("mol")[24:], benzmol[24:])

    def testTemplates(self):
        """Check for regressions to #1851"""
        smis = [
            "O=C(C1=CN=CS1)N1C2CCC1CN(CC1CC3CCC1O3)C2",
            "O=C(CC1CC1)N1C2CCC1CC(NC(=O)C13CCC(CC1)CC3)C2",
            "O=C([C@@H]1C[C@H]1C1CCC1)N1C2CCC1CN(C(=O)C13CCN(CC1)C3)C2",
            "O=C(CCN1C=CN=C1)N1C2CCC1CN(CC1CC3CCC1C3)C2"
            ]
        for s in smis:
            mol = pybel.readstring("smi", s)
            mol.draw(show=False, update=True)
        assert(True) # Segfaults before...


class NewReactionHandling(PythonBindings):

    def testBasic(self):
        smis = ["C>N>O", "C>N>", ">N>O", ">N>", "C>>", ">>O", ">>"]
        for smi in smis:
            nsmi = pybel.readstring("smi", smi).write("smi").rstrip()
            self.assertEqual(smi, nsmi)
        badsmis = ["C>>N>O", ">>>", "C>N>O>", ">", ">N", "N>"]
        for smi in badsmis:
            self.assertRaises(IOError, pybel.readstring, "smi", smi)

    def testFacade(self):
        parts = "CNO"
        mol = pybel.readstring("smi", ">".join(parts)).OBMol
        facade = ob.OBReactionFacade(mol)

        # basic test - copy out each component
        comps = [ob.REACTANT, ob.AGENT, ob.PRODUCT]
        for part, comp in zip(parts, comps):
            nmol = ob.OBMol()
            facade.GetComponent(nmol, comp, 0)
            nsmi = pybel.Molecule(nmol).write("smi").rstrip()
            self.assertEqual(nsmi, part)
        nmol = ob.OBMol()
        facade.GetComponent(nmol, ob.NO_REACTIONROLE, 0)
        self.assertEqual(0, nmol.NumAtoms())

        # Reassign role
        facade.ReassignComponent(ob.AGENT, 0, ob.NO_REACTIONROLE)
        nsmi = pybel.Molecule(mol).write("smi").rstrip()
        self.assertEqual("C>>O", nsmi)
        # ...and back again
        facade.ReassignComponent(ob.NO_REACTIONROLE, 0, ob.AGENT)

        # Add a new reactant
        molb = pybel.readstring("smi", "S").OBMol
        facade.AddComponent(molb, ob.REACTANT)
        nsmi = pybel.Molecule(mol).write("smi").rstrip()
        self.assertEqual("C.S>N>O", nsmi)

        # ...and copy it back out
        nmol = ob.OBMol()
        facade.GetComponent(nmol, ob.REACTANT, 1)
        nsmi = pybel.Molecule(nmol).write("smi").rstrip()
        self.assertEqual("S", nsmi)
        # ...and also copy the O
        facade.GetComponent(nmol, ob.PRODUCT, 0)
        nmol.SetIsReaction()
        nsmi = pybel.Molecule(nmol).write("smi").rstrip()
        self.assertEqual("S>>O", nsmi)

    def testInvalidRxn(self):
        """IsValid() should flag up invalid reaction data"""
        mol = pybel.readstring("smi", "CC>>O").OBMol
        facade = ob.OBReactionFacade(mol)
        self.assertTrue(facade.IsValid())
        mol.SetIsReaction(False)
        self.assertFalse(facade.IsValid())
        mol.SetIsReaction()
        self.assertTrue(facade.IsValid())

        atom = mol.GetAtom(1)

        facade.SetRole(atom, 4)
        self.assertFalse(facade.IsValid()) # invalid role
        facade.SetRole(atom, ob.REACTANT)
        self.assertTrue(facade.IsValid())

        data = atom.GetData("rxncomp")
        ob.toPairInteger(data).SetValue(-1)
        self.assertFalse(facade.IsValid()) # invalid rxn component id

        atom.DeleteData(data)
        self.assertFalse(atom.HasData("rxncomp"))
        self.assertFalse(facade.IsValid()) # data missing

        newdata = ob.OBPairData()
        newdata.SetAttribute("rxncomp")
        newdata.SetValue("1")
        atom.CloneData(newdata)
        self.assertTrue(atom.HasData("rxncomp"))
        self.assertFalse(facade.IsValid()) # wrong type of data

        # Connected component should not belong to two different
        # rxn components or two different reaction roles
        mol = pybel.readstring("smi", "CC>>O").OBMol
        facade = ob.OBReactionFacade(mol)
        self.assertTrue(facade.IsValid())
        atom = mol.GetAtom(1)
        facade.SetComponentId(atom, 99)
        self.assertFalse(facade.IsValid())
        facade.SetComponentId(atom, 1)
        self.assertTrue(facade.IsValid())
        facade.SetRole(atom, ob.AGENT)
        self.assertFalse(facade.IsValid())


    def testRoundtripThroughRXN(self):
        data = ["C>N>O", "C>>O", "C.N>>O", "C>>O.N",
                "C>>O", ">>O", "C>>", ">N>", ">>"]
        for rsmi in data:
            rxn = pybel.readstring("smi", rsmi).write("rxn")
            mrsmi = pybel.readstring("rxn", rxn).write("smi").rstrip()
            self.assertEqual(mrsmi, rsmi)
        # Test -G option, which changes the treatment of agents
        rsmi = "C>N>O"
        ans = {"agent": "C>N>O",
               "reactant": "C.N>>O",
               "product": "C>>O.N",
               "both": "C.N>>O.N",
               "ignore": "C>>O"}
        for option, result in ans.items():
            rxn = pybel.readstring("smi", rsmi).write("rxn", opt={"G":option})
            mrsmi = pybel.readstring("rxn", rxn).write("smi").rstrip()
            self.assertEqual(mrsmi, result)


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

class OBMolCopySubstructure(PythonBindings):
    """Tests for copying a component of an OBMol"""

    def createBitVec(self, size, bits):
        bv = ob.OBBitVec(size)
        for bit in bits:
            bv.SetBitOn(bit)
        return bv

    def testBasic(self):
        mol = pybel.readstring("smi", "ICBr")
        bv = self.createBitVec(4, (1, 3))
        nmol = ob.OBMol()
        ok = mol.OBMol.CopySubstructure(nmol, bv, None, 0)
        self.assertTrue(ok)
        self.assertEqual(pybel.Molecule(nmol).write("smi").rstrip(), "[I].[Br]")
        bv = self.createBitVec(4, (2,))
        ok = mol.OBMol.CopySubstructure(nmol, bv, None, 0)
        self.assertTrue(ok)
        self.assertEqual(pybel.Molecule(nmol).write("smi").rstrip(), "[I].[Br].[CH2]")

        mol = pybel.readstring("smi", "CCC")
        bv = self.createBitVec(4, (1,))
        bondv = self.createBitVec(2, (1,))
        nmol = ob.OBMol()
        ok = mol.OBMol.CopySubstructure(nmol, bv, bondv, 0)
        self.assertTrue(ok)

    def testNonexistentAtom(self):
        mol = pybel.readstring("smi", "ICBr")
        bv = self.createBitVec(10, (9,))
        nmol = ob.OBMol()
        ok = mol.OBMol.CopySubstructure(nmol, bv)
        self.assertFalse(ok)

    def testOptions(self):
        mol = pybel.readstring("smi", "ICBr")
        bv = self.createBitVec(4, (1, 3))
        ans = ["[I].[Br]", "I.Br", "I*.Br*"]
        ans_atomorder = [[1, 3], [1, 3], [1, 3, 2, 2]]
        ans_bondorder = [ [], [], [0, 1] ]
        for option in range(3):
            nmol = ob.OBMol()
            atomorder = ob.vectorUnsignedInt()
            bondorder = ob.vectorUnsignedInt()
            ok = mol.OBMol.CopySubstructure(nmol, bv, None, option, atomorder, bondorder)
            self.assertTrue(ok)
            self.assertEqual(pybel.Molecule(nmol).write("smi").rstrip(), ans[option])
            self.assertEqual(ans_atomorder[option], list(atomorder))
            self.assertEqual(ans_bondorder[option], list(bondorder))

    def testSpecifyBonds(self):
        mol = pybel.readstring("smi", "ICBr")
        bv = self.createBitVec(4, (1, 2, 3))
        bondbv = self.createBitVec(2, (1,))
        ans = ["I[CH2].[Br]", "IC.Br", "IC*.Br*"]
        ans_atomorder = [[1, 2, 3], [1, 2, 3], [1, 2, 3, 3, 2]]
        ans_bondorder = [ [0], [0], [0, 1, 1] ]
        for option in range(3):
            nmol = ob.OBMol()
            atomorder = ob.vectorUnsignedInt()
            bondorder = ob.vectorUnsignedInt()
            ok = mol.OBMol.CopySubstructure(nmol, bv, bondbv, option, atomorder, bondorder)
            self.assertTrue(ok)
            self.assertEqual(pybel.Molecule(nmol).write("smi").rstrip(), ans[option])
            self.assertEqual(ans_atomorder[option], list(atomorder))
            self.assertEqual(ans_bondorder[option], list(bondorder))

    def testSpecifyAtomsAndBonds(self):
        # Now copy just a subset of atoms too
        mol = pybel.readstring("smi", "ICBr")
        bv = self.createBitVec(4, (1, 3))
        bondbv = self.createBitVec(2, (1,))
        ans = ["[I].[Br]", "I.Br", "I*.Br*"]
        for option in range(3):
            nmol = ob.OBMol()
            ok = mol.OBMol.CopySubstructure(nmol, bv, bondbv, option)
            self.assertTrue(ok)
            self.assertEqual(pybel.Molecule(nmol).write("smi").rstrip(), ans[option])

    def testBondOrders(self):
        mol = pybel.readstring("smi", "O=C=O")
        bv = self.createBitVec(3, (2, 3))
        bondbv = self.createBitVec(2, (1,))
        ans = ["[C].[O]", "C.O", "C(=*)=*.O=*"]
        for option in range(3):
            nmol = ob.OBMol()
            ok = mol.OBMol.CopySubstructure(nmol, bv, bondbv, option)
            self.assertTrue(ok)
            self.assertEqual(pybel.Molecule(nmol).write("smi").rstrip(), ans[option])
        ans = ["[C]=O", "C=O", "C(=O)=*"]
        for option in range(3):
            nmol = ob.OBMol()
            ok = mol.OBMol.CopySubstructure(nmol, bv, None, option)
            self.assertTrue(ok)
            self.assertEqual(pybel.Molecule(nmol).write("smi").rstrip(), ans[option])

    def testStereo(self):
        data = [
                ("FC[C@@](Br)(Cl)I",
                    [((2, 3, 4, 5, 6), None, "C[C@@](Br)(Cl)I"),
                    ((2, 3, 4, 5), None, "CC(Br)Cl"),
                    ((1, 2, 3, 4, 5, 6), (4,), "FCC(Br)Cl.I")]
                ),
                ("[C@@H](Br)(Cl)I",
                    [((1, 2, 3), None, "C(Br)Cl"),
                    ((1, 2, 3, 4), (2,), "C(Br)Cl.I")]
                ),
                ("C[C@@H]1CO1",
                    [((2, 3, 4), None, "C1CO1"),]
                ),
                ("F/C=C/I",
                    [
                     ((1, 2, 3, 4), None, "F/C=C/I"),
                     ((1, 2, 3), None, "FC=C"),
                     ((1, 2, 3, 4), (0,), "F.C=CI"),
                     ((1, 2, 3, 4), (1,), "FC.CI")]
                ),
               ]
        for smi, d in data:
            mol = pybel.readstring("smi", smi)
            for a, b, ans in d:
                nmol = ob.OBMol()
                bv = self.createBitVec(7, a)
                bondbv = None if b is None else self.createBitVec(5, b)
                ok = mol.OBMol.CopySubstructure(nmol, bv, bondbv)
                self.assertTrue(ok)
                if "@" not in ans and "/" not in ans:
                    self.assertFalse(nmol.GetData(ob.StereoData))
                self.assertEqual(pybel.Molecule(nmol).write("smi").rstrip(),
                                 ans)

    def testResidueCopying(self):
        smi = "C[C@@H](C(=O)N[C@@H](CS)C(=O)O)N" # H-Ala-Cys-OH
        mol = pybel.readstring("smi", smi).OBMol
        ob.cvar.chainsparser.PerceiveChains(mol)
        mol.SetChainsPerceived()
        residues = list(ob.OBResidueIter(mol))
        self.assertEqual(len(residues), 2)

        # Copy just the Cys N
        bv = self.createBitVec(mol.NumAtoms() + 1, (5,))
        nmol = ob.OBMol()
        mol.CopySubstructure(nmol, bv)
        self.assertEqual(len(list(ob.OBResidueIter(nmol))), 1)
        pdb = pybel.Molecule(nmol).write("pdb")
        atoms = [line for line in pdb.split("\n") if line.startswith("ATOM")]
        self.assertEqual(len(atoms), 1)
        cysN = "ATOM      1  N   CYS A   2"
        self.assertTrue(atoms[0].startswith(cysN))

        # Copy the Cys N and Ca
        bv = self.createBitVec(mol.NumAtoms() + 1, (5, 6))
        nmol.Clear()
        mol.CopySubstructure(nmol, bv)
        self.assertEqual(len(list(ob.OBResidueIter(nmol))), 1)
        pdb = pybel.Molecule(nmol).write("pdb")
        atoms = [line for line in pdb.split("\n") if line.startswith("ATOM")]
        self.assertEqual(len(atoms), 2)
        self.assertTrue(atoms[0].startswith(cysN))
        cysCa = "ATOM      2  CA  CYS A   2"
        self.assertTrue(atoms[1].startswith(cysCa))

        # Copy the Ala C and Cys N
        bv = self.createBitVec(mol.NumAtoms() + 1, (3, 5))
        nmol.Clear()
        mol.CopySubstructure(nmol, bv)
        self.assertEqual(len(list(ob.OBResidueIter(nmol))), 2)
        pdb = pybel.Molecule(nmol).write("pdb")
        atoms = [line for line in pdb.split("\n") if line.startswith("ATOM")]
        self.assertEqual(len(atoms), 2)
        alaC = "ATOM      1  C   ALA A   1"
        self.assertTrue(atoms[0].startswith(alaC))
        cysN = "ATOM      2  N   CYS A   2"
        self.assertTrue(atoms[1].startswith(cysN))

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
        ac = ob.OBPairInteger()
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

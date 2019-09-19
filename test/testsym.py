"""Test OpenBabel executables from Python

Note: Python bindings not used

On Windows or Linux, you can run these tests at the commandline
in the build folder with:
"C:\Program Files\CMake 2.6\bin\ctest.exe" -C CTestTestfile.cmake
                                           -R pytest -VV

You could also "chdir" into build/test and run the test file directly:
python ../../../test/testsym.py

In both cases, the test file is run directly from the source folder,
and so you can quickly develop the tests and try them out.
"""

import os
import unittest

from testbabel import run_exec, BaseTest

class TestSym(BaseTest):
    """Base class for a series of tests relating to symmetry"""

    def testInChItoSMI(self):
        """Verify that the InChI is read correctly"""
        output, error = run_exec(self.inchi, "obabel -iinchi -ocan")
        self.assertEqual(output.rstrip(), self.cansmi)

    def testSMItoInChI(self):
        """Verify that all molecules give the same InChI"""
        output, error = run_exec("\n".join(self.smiles), "obabel -ismi -oinchi")
        output = "\n".join([x.rstrip() for x in output.split("\n")])
        self.assertEqual(output.rstrip(), "\n".join([self.inchi] * len(self.smiles)))

    def testSMItoCAN(self):
        """Verify that all molecules give the same cansmi"""
        output, error = run_exec("\n".join(self.smiles), "obabel -ismi -ocan")
        output = "\n".join([x.rstrip() for x in output.split("\n")])
        self.assertEqual(output.rstrip(), "\n".join([self.cansmi] * len(self.smiles)))

    def testSMIthruXML(self):
        """Verify that roundtripping through CML preserves stereo"""
        output, error = run_exec("\n".join(self.smiles), "obabel -ismi -O tmp.cml")
        output, error = run_exec(output.rstrip(), "obabel -icml tmp.cml -ocan")
        output = "\n".join([x.rstrip() for x in output.split("\n")])
        self.assertEqual(output.rstrip(), "\n".join([self.cansmi] * len(self.smiles)))
        os.remove("tmp.cml")

class TestTetSym(TestSym):
    """A series of tests relating to tetrahedral symmetry"""

    def setUp(self):
        self.canFindExecutable("obabel")

        # The following all represent the same molecule
        self.cansmi = "C[C@](Br)(Cl)F"
        self.inchi = "InChI=1S/C2H3BrClF/c1-2(3,4)5/h1H3/t2-/m0/s1"
        self.smiles = [
             'C[C@@](Cl)(Br)F',
             'C[C@](Cl)(F)Br',
             'C[C@](Br)(Cl)F',
             'C[C@@](Br)(F)Cl',
             'C[C@@](F)(Cl)Br',
             'C[C@](F)(Br)Cl',
             'Cl[C@](C)(Br)F',
             'Cl[C@@](C)(F)Br',
             'Cl[C@@](Br)(C)F',
             'Cl[C@](Br)(F)C',
             'Cl[C@](F)(C)Br',
             'Cl[C@@](F)(Br)C',
             'Br[C@@](C)(Cl)F',
             'Br[C@](C)(F)Cl',
             'Br[C@](Cl)(C)F',
             'Br[C@@](Cl)(F)C',
             'Br[C@@](F)(C)Cl',
             'Br[C@](F)(Cl)C',
             'F[C@](C)(Cl)Br',
             'F[C@@](C)(Br)Cl',
             'F[C@@](Cl)(C)Br',
             'F[C@](Cl)(Br)C',
             'F[C@](Br)(C)Cl',
             'F[C@@](Br)(Cl)C'
             ]


class TestCisTransSym(TestSym):
    """A series of tests relating to cistrans symmetry"""

    def setUp(self):
        self.canFindExecutable("obabel")

        # The following all represent the same molecule
        self.cansmi = "Cl/C=C/C=C\\Br"
        self.inchi = "InChI=1S/C4H4BrCl/c5-3-1-2-4-6/h1-4H/b3-1-,4-2+"
        self.smiles = [
                "C(=C\C=C/Br)/Cl",
                "Cl/C=C/C=C\Br",
                "Br/C=C\C=C\Cl",
                "C(=C\Cl)/C=C\Br",
                "C(=C\C=C\Cl)\Br",
                "C(=C\Br)\C=C\Cl"
                ]

class TestLonePairTetSym(TestSym):
    """A series of tests relating to tet symmetry involving a lone pair"""

    def setUp(self):
        self.canFindExecutable("obabel")

        # The following all represent the same molecule
        self.cansmi = "C[S@](=O)Cl"
        self.inchi = "InChI=1S/CH3ClOS/c1-4(2)3/h1H3/t4-/m0/s1"
        self.smiles = [
                self.cansmi,
                "O=[S@](Cl)C",
                "O=[S@@](C)Cl",
                "[S@](Cl)(=O)C",
                ]

class TestRingBondCisTransSym(TestSym):
    """A series of tests relating to tet symmetry involving a lone pair"""

    def setUp(self):
        self.canFindExecutable("obabel")

        # The following all represent the same molecule
        self.cansmi = r"I/C=C/1\CN1"
        self.inchi = "InChI=1S/C3H4IN/c4-1-3-2-5-3/h1,5H,2H2/b3-1+"
        self.smiles = [
                self.cansmi,
                r"I/C=C\1/NC1",
                r"I/C=C1NC/1",
                 "I/C=C/1/NC/1",
                ]

class TestConversions(BaseTest):
    """A series of tests relating to file format conversions and symmetry"""

    def setUp(self):
        self.canFindExecutable("obabel")
        self.data = [
('ClC=CF', 'FC=CCl',       'InChI=1S/C2H2ClF/c3-1-2-4/h1-2H'),
('ClC=CF', 'FC=CCl',       'InChI=1S/C2H2ClF/c3-1-2-4/h1-2H'),
('Cl/C=C/F', 'F/C=C/Cl',   'InChI=1S/C2H2ClF/c3-1-2-4/h1-2H/b2-1+'),
(r"Cl/C=C\F", r"F/C=C\Cl", 'InChI=1S/C2H2ClF/c3-1-2-4/h1-2H/b2-1-'),
('Cl[C@@](Br)(F)I', 'F[C@](I)(Br)Cl', 'InChI=1S/CBrClFI/c2-1(3,4)5/t1-/m0/s1'),
('Cl[C@](Br)(F)I', 'F[C@@](I)(Br)Cl',   'InChI=1S/CBrClFI/c2-1(3,4)5/t1-/m1/s1'),
('ClC(Br)(F)I', 'FC(I)(Br)Cl',         'InChI=1S/CBrClFI/c2-1(3,4)5'),
('O=[S@@](Cl)I', "Cl[S@](=O)I", "InChI=1S/ClIOS/c1-4(2)3/t4-/m0/s1"),
('O=[S@](Cl)I', "Cl[S@@](=O)I", "InChI=1S/ClIOS/c1-4(2)3/t4-/m1/s1"),
('O=S(Cl)I', "ClS(=O)I", "InChI=1S/ClIOS/c1-4(2)3"),
(r"IC=C1NC1", r"IC=C1CN1", "InChI=1S/C3H4IN/c4-1-3-2-5-3/h1,5H,2H2"),
(r"I/C=C\1/NC1", r"I/C=C/1\CN1", "InChI=1S/C3H4IN/c4-1-3-2-5-3/h1,5H,2H2/b3-1+"),
(r"I/C=C/1\NC1", r"I/C=C\1/CN1", "InChI=1S/C3H4IN/c4-1-3-2-5-3/h1,5H,2H2/b3-1-"),
]

    def testSMILEStoInChI(self):
        # Tests interconversions between the SMILES on the left versus
        # the InChI on the right.
        # The canonical smiles (in the middle) were derived from the SMILES.
        for smiles, can, inchi in self.data:
            output, error = run_exec(smiles, "obabel -ismi -oinchi")
            self.assertEqual(output.rstrip(), inchi)
            output, error = run_exec(inchi, "obabel -iinchi -ocan")
            self.assertEqual(output.rstrip(), can)

    def parseMDL(self, text):
        lines = text.split("\n")
        broken = lines[3].split()
        Natoms = int(broken[0])
        Nbonds = int(broken[1])
        atoms = []
        for i in range(Natoms):
            broken = lines[i+4].split()
            atoms.append({'parity':int(broken[6])})
        bonds = []
        for i in range(Nbonds):
            broken = lines[i+4+Natoms].split()
            bonds.append({'stereo':int(broken[3])})
        return atoms, bonds

    def testSMILESto2D(self):
        """Test gen2d for some basic cases"""
        for smi, can, inchi in self.data:
            output, error = run_exec(smi, "obabel -ismi --gen2d -omdl")
            output, error = run_exec(output.rstrip(), "obabel -imdl -ocan")
            self.assertEqual(can, output.rstrip())

    def testSMILESto3DMDL(self):
        """Test interconversion between SMILES and 3D MDL"""
        data = [
([0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 3]), # 'ClC=CF'
([0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 3]), # 'ClC=CF'
([0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0]), # 'Cl/C=C/F'
([0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0]), # 'Cl/C=C\\F'
# The bond parities are irrelevant/meaningless for the next two
([0, 0, 0, 0, 1], []), # 'Cl[C@@](Br)(F)I'
([0, 0, 0, 0, 2], []), # 'Cl[C@](Br)(F)I'
([0, 0, 0, 0, 3], [0, 0, 0, 0]), # 'ClC(Br)(F)I'
([0, 0, 0, 1], []), # 'O=[S@@](Cl)I),
([0, 0, 0, 2], []), # 'O=[S@](Cl)I),
([0, 0, 0, 3], []), # 'O=S(Cl)I),
([0]*9, [0]*8 + [3]), #  "IC=C1NC1"
([0]*9, [0]*9), # r"I/C=C\1/NC1"
([0]*9, [0]*9), # r"I/C=C/1\NC1"
]
        for i, (atompar, bondstereo) in enumerate(data):
            smiles, can = self.data[i][0:2]
            output, error = run_exec(smiles, "obabel -ismi -osdf --gen3d")
            atoms, bonds = self.parseMDL(output)
            parities = [atom['parity'] for atom in atoms]
            parities.sort()
            stereos = [bond['stereo'] for bond in bonds]
            stereos.sort()
            self.assertEqual(atompar, parities)
            if bondstereo:
                self.assertEqual(bondstereo, stereos)
            output, error = run_exec(output, "obabel -isdf -as -ocan")
            # "-as" is necessary to identify the unknown stereo
            self.assertEqual(output.rstrip(), can)

    def testXYZtoSMILESand3DMDL(self):
        """Test conversion from XYZ to SMILES and 3D MDL"""
        # Since the XYZ format does not trigger stereo perception,
        # this test makes sure that the SMILES and 3D MDL formats
        # perceive stereo themselves.
        data = [
([0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 3]), # 'ClC=CF'
([0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 3]), # 'ClC=CF'
([0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0]), # 'Cl/C=C/F'
([0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0]), # 'Cl/C=C\\F'
# The bond parities are irrelevant/meaningless for the next two
([0, 0, 0, 0, 1], []), # 'Cl[C@@](Br)(F)I'
([0, 0, 0, 0, 2], []), # 'Cl[C@](Br)(F)I'
([0, 0, 0, 0, 3], [0, 0, 0, 4]), # 'ClC(Br)(F)I'
([0, 0, 0, 1], []), # 'O=[S@@](Cl)I),
([0, 0, 0, 2], []), # 'O=[S@](Cl)I),
([0, 0, 0, 3], []), # 'O=S(Cl)I),
([0]*9, [0]*8 + [3]), #  "IC=C1NC1"
([0]*9, [0]*9), # r"I/C=C\1/NC1"
([0]*9, [0]*9), # r"I/C=C/1\NC1"
]
        for i, (atompar, bondstereo) in enumerate(data):
            if i in [0, 1, 6, 10]: continue # ambiguous stereo is lost in XYZ
            if i in [7, 8, 9]: continue # perception of S=O from XYZ fails

            smiles, can = self.data[i][0:2]
            output, error = run_exec(smiles, "obabel -ismi -oxyz --gen3d")

            canoutput, error = run_exec(output, "obabel -ixyz -ocan")
            self.assertEqual(canoutput.rstrip(), can)

            sdfoutput, error = run_exec(output, "obabel -ixyz -osdf")
            atoms, bonds = self.parseMDL(sdfoutput)
            parities = [atom['parity'] for atom in atoms]
            parities.sort()
            stereos = [bond['stereo'] for bond in bonds]
            stereos.sort()
            self.assertEqual(atompar, parities)
            if bondstereo:
                self.assertEqual(bondstereo, stereos)

    def test2DMDLto0D(self):
        """Test conversion for 2D MDL to CAN and InChI"""
        # The following file was created using RDKit starting from
        # the SMILES strings in data[x][0] below.
        filename = self.getTestFile("testsym_2Dtests.sdf")

        output, error = run_exec("obabel -isdf %s -ocan" % filename)
        for i, smiles in enumerate(output.rstrip().split("\n")):
            self.assertEqual(smiles.rstrip(), self.data[i][1])

        output, error = run_exec("obabel -isdf %s -oinchi" % filename)
        for i, inchi in enumerate(output.rstrip().split("\n")):
            self.assertEqual(inchi.rstrip(), self.data[i][2])

    def test2DMDLto0D_more(self):
        """Test various combinations of stereobonds in 2D perception"""
        filenames = [self.getTestFile(x) for x in
                     ["testsym_2Dtests_more.sdf",
                      "testsym_2Dtests_threeligands.sdf"]]
        # The test files have the correct canonical SMILES string
        # stored in the data field "smiles"

        output, error = run_exec("obabel -isdf %s %s -ocan --append smiles" %
                                 (filenames[0], filenames[1]))
        for line in output.rstrip().split("\n"):
            result, correct_answer = line.split()
            self.assertEqual(result, correct_answer)

    def test2DMDLto2DMDL(self):
        """Make sure that stereo is preserved when writing wedge bonds"""
        filenames = [self.getTestFile(x) for x in
                     ["testsym_2Dtests_more.sdf",
                      "testsym_2Dtests_threeligands.sdf"]]
        # The test files have the correct canonical SMILES string
        # stored in the data field "smiles"

        output, error = run_exec("obabel -isdf %s %s -osdf --append smiles" %
                                 (filenames[0], filenames[1]))
        finaloutput, error = run_exec(output, "obabel -isdf -ocan")
        for line in finaloutput.rstrip().split("\n"):
            result, correct_answer = line.split()
            self.assertEqual(result, correct_answer)

    def testSMILESto0DMDL(self):
        """Test interconversion between SMILES and 0D MDL"""
        data = [
([0, 0, 0, 0, 1], [0, 0, 0, 0]), # 'Cl[C@@](Br)(F)I'
([0, 0, 0, 0, 2], [0, 0, 0, 0]), # 'Cl[C@](Br)(F)I'
([0, 0, 0, 0, 0], [0, 0, 0, 0]), # 'ClC(Br)(F)I'
([0, 0, 0, 0, 3], [0, 0, 0, 0])  # 'ClC(Br)(F)I' option 'S' when reading
]
        for i, (atompar, bondstereo) in enumerate(data):
            if i == 3:
                smiles, can = self.data[6][0:2]
                output, error = run_exec(smiles, "obabel -ismi -osdf -aS")
            else:
                smiles, can = self.data[i + 4][0:2]
                output, error = run_exec(smiles, "obabel -ismi -osdf")
            atoms, bonds = self.parseMDL(output)
            parities = [atom['parity'] for atom in atoms]
            parities.sort()
            stereos = [bond['stereo'] for bond in bonds]
            stereos.sort()
            self.assertEqual(atompar, parities)
            self.assertEqual(bondstereo, stereos)
            output, error = run_exec(output, "obabel -isdf -as -ocan")
            self.assertEqual(output.rstrip(), can)


class TestStereoConversion(BaseTest):
    """Random tests relating to roundtripping stereochemistry"""
    def setUp(self):
        self.canFindExecutable("obabel")
    def testInChIToSMILES_Bug(self):
        """PR#2101034- InChI <-> SMILES conv misrepresents stereo"""
        test_inchi = 'InChI=1S/C10H10/c1-2-3-7-10-8-5-4-6-9-10/h2-9H,1H2/b7-3+'
        output, error = run_exec(test_inchi, "obabel -iinchi -osmi")
        self.assertEqual(output.rstrip(), "C=C/C=C/c1ccccc1")

        test_smiles = "C=C\C=C/c1ccccc1"
        output, error = run_exec(test_smiles, "obabel -ismi -oinchi")
        self.assertEqual(output.rstrip(), "InChI=1S/C10H10/c1-2-3-7-10-8-5-4-6-9-10/h2-9H,1H2/b7-3-")
    def testChiralToLonePair(self):
        """PR#3058701 - Handle stereochemistry at lone pair on S"""
        # Note to self: Need to ensure that roundtripping through the various
        # 2D and 3D formats works. In the meanwhile, this test at least ensures
        # that SMILES reading and writing works fine.
        can = 'C[S@](=O)Cl'
        smiles = [can, '[S@](Cl)(=O)C', 'O=[S@](Cl)C']
        for smile in smiles:
            output, error = run_exec(smile, "obabel -ismi -ocan")
            self.assertEqual(output.rstrip(), can)
        # Check that regular chiral S still work fine
        smi = "[S@](=O)(=N)(C)O"
        output, error = run_exec(smi, "obabel -ismi -osmi")
        self.assertEqual(output.rstrip(), smi)

del TestSym # remove base class to avoid tests

if __name__ == "__main__":
    unittest.main()

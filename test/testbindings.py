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

    def testKekulizationOfcn(self):
        """We were previously not reading 'cn' correctly, or at least how
        Daylight would"""
        mol = pybel.readstring("smi", "cn")
        self.assertEqual("C=N", mol.write("smi").rstrip())

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
        smi = "[*][H][He][Li][Be][B][C][N][O][F][Ne][Na][Mg][Al][Si][P][S][Cl][Ar][K][Ca][Sc][Ti][V][Cr][Mn][Fe][Co][Ni][Cu][Zn][Ga][Ge][As][Se][Br][Kr][Rb][Sr][Y][Zr][Nb][Mo][Tc][Ru][Rh][Pd][Ag][Cd][In][Sn][Sb][Te][I][Xe][Cs][Ba][La][Ce][Pr][Nd][Pm][Sm][Eu][Gd][Tb][Dy][Ho][Er][Tm][Yb][Lu][Hf][Ta][W][Re][Os][Ir][Pt][Au][Hg][Tl][Pb][Bi][Po][At][Rn][Fr][Ra][Ac][Th][Pa][U][Np][Pu][Am][Cm][Bk][Cf][Es][Fm][Md][No][Lr][Rf][Db][Sg][Bh][Hs][Mt][Ds][Rg][Cn][Nh][Fl][Mc][Lv][Ts][Og]"
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
        self.assertEqual(mol.write("smi").rstrip(), "")

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

class TestBFSiter(PythonBindings):
    def testZeroDepth(self):
        mol_txt = """
 OpenBabel09261715552D

 12 12  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  6  1  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  7 12  1  0  0  0  0
  7  8  2  0  0  0  0
  8  9  1  0  0  0  0
  9 10  2  0  0  0  0
 10 11  1  0  0  0  0
 11 12  2  0  0  0  0
M  END
$$$$
"""
        mol = pybel.readstring('sdf', mol_txt)

        fp = []
        for i in range(len(mol.atoms)):
            envs = {}
            for atom, current_depth in pybel.ob.OBMolAtomBFSIter(mol.OBMol, i + 1):
                if current_depth not in envs:
                    envs[current_depth] = []
                envs[current_depth].append(atom.GetIdx())

                if current_depth == 1:
                    self.assertEqual(1, len(envs[current_depth]))


if __name__ == "__main__":
    unittest.main()

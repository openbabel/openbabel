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
        elems = ["N", "P", "As", "Sb", "Bi"]
        for elem in elems:
            mol = pybel.readstring("smi", "c1cc[%sH]c1" % elem.lower(), opt={"a": True})
            self.assertNotEqual(0, mol.OBMol.NumAtoms())
        # - like O
        elems = ["O", "S", "Se", "Te"]
        for elem in elems:
            mol = pybel.readstring("smi", "c1cc[%s]c1" % elem.lower(), opt={"a": True})
            self.assertNotEqual(0, mol.OBMol.NumAtoms())
        # - untested and unsupported at the moment, aromatic B


    def testSmilesParsingAndWritingOfLargeIsotopes(self):
        smis = ["[1C]", "[11C]", "[111C]", "[1111C]"]
        for smi in smis:
            mol = pybel.readstring("smi", smi)
            self.assertEqual(mol.write("smi").rstrip(), smi)
        self.assertRaises(IOError, pybel.readstring, "smi", "[11111C]")
        mol = pybel.readstring("smi", "[C]")
        mol.atoms[0].OBAtom.SetIsotope(65535)
        self.assertEqual(mol.write("smi").rstrip(), "")

if __name__ == "__main__":
    unittest.main()

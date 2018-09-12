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
    
class PybelWrapper(PythonBindings):
    def testDummy(self):
        self.assertTrue(pybel is not None, "Failed to import the Pybel module")
   
if __name__ == "__main__":
    unittest.main()

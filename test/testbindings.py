"""Test OpenBabel Python bindings

On Windows or Linux, you can run these tests at the commandline
in the build folder with:
"C:\Program Files\CMake 2.6\bin\ctest.exe" -C CTestTestfile.cmake
                                           -R pybindtest -VV

You could also "chdir" into build/test and run the test file directly:
python ../../test/testbindings.py

In both cases, the test file is run directly from the source folder,
and so you can quickly develop the tests and try them out.
"""

import os
import re
import sys
import unittest

modulelocation = os.path.join("..", "Release")
sys.path = [modulelocation] + sys.path

try:
    import openbabel as ob
    # We check later to make sure we have imported the right module
except ImportError:
    ob = None

class PythonBindings(unittest.TestCase):
    def setUp(self):
        self.assertTrue(ob is not None, "Failed to import the openbabel module")
        self.assertTrue(os.path.isfile(os.path.join(
            modulelocation, "openbabel.py")), "openbabel module not found")
    def testSimple(self):
        mol = ob.OBMol()
        conv = ob.OBConversion()
        conv.SetInFormat("smi")
        conv.ReadString(mol, "CC(=O)Cl")
        self.assertAlmostEqual(mol.GetMolWt(), 78.5, 1)
    
pybellocation = os.path.join("..", "..", "..", "scripts", "python")
sys.path = [pybellocation] + sys.path

try:
    import pybel
except ImportError:
    pybel = None
    
class PybelWrapper(unittest.TestCase):
    def testDummy(self):
        self.assertTrue(pybel is not None, "Failed to import the Pybel module")
        self.assertTrue(os.path.isfile(os.path.join(
            pybellocation, "pybel.py")), "Pybel module not found")
   
testlocation = os.path.join("..", "..", "..", "scripts", "python", "examples")
sys.path = [testlocation] + sys.path

try:
    import testpybel
except ImportError:
    testpybel = None

def gettests():
    testsuite = []
    for myclass in [PybelWrapper, PythonBindings]:
        suite = unittest.TestLoader().loadTestsFromTestCase(myclass)
        testsuite.append(suite)
    suite = unittest.TestLoader().loadTestsFromName("TestOBPybel", testpybel)
    testsuite.append(suite)
    return testsuite

if __name__ == "__main__":
    unittest.TextTestRunner().run(unittest.TestSuite(gettests()))

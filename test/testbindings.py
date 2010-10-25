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

here = sys.path[0]

iswin = sys.platform.startswith("win")

if iswin:
    modulelocation = os.path.join("..", "Release")
else:
    modulelocation = os.path.join("..", "scripts", "pybuild")
sys.path = [modulelocation] + sys.path

try:
    import openbabel as ob
    # We check later to make sure we have imported the right module
except ImportError:
    ob = None

if iswin:
    pybellocation = os.path.join(here, "..", "scripts", "python")
else:
    pybellocation = os.path.join("..", "scripts")
sys.path = [pybellocation] + sys.path

try:
    if iswin:
        import pybel_py2x as pybel
    else:
        import pybel
except ImportError:
    pybel = None

class PythonBindings(unittest.TestCase):
    def setUp(self):
        self.assertTrue(ob is not None, "Failed to import the openbabel module")
        self.assertTrue(os.path.isfile(os.path.join(
            modulelocation, "openbabel.py")), "openbabel module not found")

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
        self.assertTrue(os.path.isfile(os.path.join(
            pybellocation, "pybel.py")), "Pybel module not found at %s" % pybellocation)
   
if __name__ == "__main__":
    unittest.main()

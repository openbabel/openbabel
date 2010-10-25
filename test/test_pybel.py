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
import sys

from testbindings import *

testlocation = os.path.join(here, "..", "scripts", "python", "examples")
sys.path = [testlocation] + sys.path

try:
    import testpybel
except ImportError:
    testpybel = None

def gettests():
    suite = unittest.TestLoader().loadTestsFromName("TestOBPybelNoDraw", testpybel)
    return suite

if __name__ == "__main__":
    unittest.TextTestRunner().run(unittest.TestSuite(gettests()))

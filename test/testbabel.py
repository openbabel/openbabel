"""Test OpenBabel executables from Python

Note: Python bindings not used

On Windows or Linux, you can run these tests at the commandline
in the build folder with:
"C:\Program Files\CMake 2.6\bin\ctest.exe" -C CTestTestfile.cmake
                                           -R pytest -VV

You could also "chdir" into build/test and run the test file directly:
python ../../test/testbabel.py

In both cases, the test file is run directly from the source folder,
and so you can quickly develop the tests and try them out.
"""

import os
import sys
import unittest

from subprocess import Popen, PIPE

def run_exec(text, commandline):
    """Pipe from stdin to an executable

    Example: run_exec("CC(=O)Cl", "babel -ismi -oinchi")

    Return a tuple (stdout, stderr)
    """
    broken = commandline.split()
    exe = executable(broken[0])
    p = Popen([exe] + broken[1:], 
              stdin=PIPE, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate(text)
    return stdout, stderr

def executable(name):
    """Return the full path to an executable"""
    suffix = ""
    folder = "bin"
    if sys.platform == "win32":
        suffix = ".exe"
        folder = "Release"
    return os.path.join("..", folder, name + suffix)

def log(text):
    """Convenience function for debugging tests

    The log file (log.txt) is created in build/test
    """
    output = open("log.txt", "a")
    print >> output, text
    output.close()

class BaseTest(unittest.TestCase):
    def canFindExecutable(self, name):
        fullpath = executable(name)
        self.assertTrue(os.path.isfile(fullpath),
                        "'%s' executable not found at %s" % (name, fullpath))

class testBabel(BaseTest):
    """A series of tests relating to the Babel executable"""
        
    def testSMItoInChI(self):
        self.canFindExecutable("babel")
        output, error = run_exec("CC(=O)Cl", "babel -ismi -oinchi")
        self.assertEqual(output.rstrip(), "InChI=1S/C2H3ClO/c1-2(3)4/h1H3")

if __name__ == "__main__":
    unittest.main()

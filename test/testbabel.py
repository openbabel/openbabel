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
import re
import sys
import unittest

from subprocess import Popen, PIPE

def run_exec(*args):
    """Run one of OpenBabel's executables

    With two arguments (stdin, commandline) it pipes
    the stdin through the executable.

    Example: run_exec("CC(=O)Cl", "babel -ismi -oinchi")

    Return a tuple (stdout, stderr)
    """
    if len(args) == 2:
        text, commandline = args
    elif len(args) == 1:
        text, commandline = "", args[0]
    else:
        raise Exception("One or two arguments expected")
    
    broken = commandline.encode().split()
    exe = executable(broken[0].decode())
    # Note that bufsize = -1 means default buffering
    # Without this, it's unbuffered and it takes 10x longer on MacOSX
    if text:
        p = Popen([exe] + broken[1:], 
                  stdin=PIPE, stdout=PIPE, stderr=PIPE, bufsize=1)
        stdout, stderr = p.communicate(text.encode())
    else:
        p = Popen([exe] + broken[1:], 
                  stdout=PIPE, stderr=PIPE, bufsize=-1)
        stdout, stderr = p.communicate()
    
    return stdout.decode(), stderr.decode()

def executable(name):
    """Return the full path to an executable"""
    suffix = ""
    folder = "bin"
    if sys.platform == "win32":
        suffix = ".exe"
        folder = os.path.join(folder, "Release")
    return os.path.join("..", folder, name + suffix)

def log(text):
    """Convenience function for debugging tests

    The log file (log.txt) is created in build/test
    """
    output = open("log.txt", "a")
    print >> output, text
    output.close()

class BaseTest(unittest.TestCase):
    """A base class for test classes that adds additional
    test methods"""
    
    def canFindExecutable(self, name):
        fullpath = executable(name)
        self.assertTrue(os.path.isfile(fullpath),
                        "'%s' executable not found at %s" % (name, fullpath))
        
    def canFindFile(self, filename):
        self.assertTrue(os.path.isfile(filename),
                        "Cannot find the file '%s'" % filename)

    def getTestFile(self, filename):
        here = sys.path[0]
        fullpath = os.path.join(here, "files", filename)
        self.canFindFile(fullpath)
        return fullpath

    def assertConverted(self, stderr, N):
        """Assert that N molecules were converted."""
        pat = "(-?\d.*) molecule(?:s?) converted"
        lines = stderr.split("\n")
        convertedline =  [line for line in lines if re.match(pat, line)]
        if len(convertedline) == 0:
            self.fail("Cannot find the number of molecules converted")
        conversion_no = int(re.findall(pat, convertedline[0])[0])
        self.assertEqual(N, conversion_no,
                         "Number of molecules converted is %d "
                         "but should be %d" % (conversion_no, N))        

class testOBabel(BaseTest):
    """A series of tests relating to the obabel executable"""
        
    def testSMItoInChI(self):
        self.canFindExecutable("obabel")
        output, error = run_exec("CC(=O)Cl", "obabel -ismi -oinchi")
        self.assertEqual(output.rstrip(), "InChI=1S/C2H3ClO/c1-2(3)4/h1H3")

    def testRSMItoRSMI(self):
        # Check possible combinations of missing rxn components
        data = ["O>N>S", "O>>S", "O>N>", "O>>",
                ">N>S", ">>S", ">>"]
        for rsmi in data:
            output, error = run_exec('obabel -:%s -irsmi -orsmi' % rsmi)
            self.assertEqual(output.rstrip(), rsmi)
        # Check handling of invalid rxn components
        data = ["Noel>N>S", "O>Noel>S", "O>N>Noel"]
        errors = ["reactant", "agent", "product"]
        for rsmi, error in zip(data, errors):
            output, errormsg = run_exec('obabel -:%s -irsmi -orsmi' % rsmi)
            self.assertTrue(error in errormsg)

    def sort(self, rsmi):
        # TODO: Change OBMol.Separate to preserve the order. This
        # function shouldn't be necessary.
        tmp = rsmi.split(">")
        r = ".".join(sorted(tmp[0].split(".")))
        p = ".".join(sorted(tmp[2].split(".")))
        return "%s>%s>%s" % (r, tmp[1], p)

    def testRingClosures(self):
        # Test positives
        data = ["c1ccccc1", "c%11ccccc%11", "c%(1)ccccc%(1)", "c%(51)ccccc%51",
                "c%(99999)ccccc%(99999)"]
        for smi in data:
            output, error = run_exec("obabel -:%s -osmi" % smi)
            self.assertEqual("c1ccccc1", output.rstrip())
        # Test negatives
        data = ["c%1ccccc%1", "c%a1cccc%a1", "c%(a1)ccccc%(a1)",
                "c%(000001)ccccc%(000001)", "c%(51)ccccc%(15)"]
        for smi in data:
            output, error = run_exec("obabel -:%s -osmi" % smi)
            self.assertTrue("0 molecules converted" in error)
        # Now test writing of %(NNN) notation
        output, error = run_exec("obabel %s -osmi" % self.getTestFile("102Uridine.smi"))
        self.assertTrue("%(100)" in output)

    def testRoundtripThroughRXN(self):
        self.canFindExecutable("obabel")
        data = ["C>N>O", "C>>O", "C.N>>O", "C>>O.N",
                "C>>O", ">>O", "C>>", ">N>", ">>"]
        for rsmi in data:
            output, error = run_exec("obabel -irsmi -:%s -orxn" % rsmi)
            moutput, error = run_exec(output, "obabel -irxn -orsmi")
            self.assertEqual(self.sort(moutput.rstrip()), self.sort(rsmi))
        rsmi = "C>N>O"
        ans = {"agent": "C>N>O",
               "reactant": "C.N>>O",
               "product": "C>>O.N",
               "both": "C.N>>O.N",
               "ignore": "C>>O"}
        for option, result in ans.iteritems():
            output, error = run_exec("obabel -irsmi -:%s -orxn -xG %s" %
                                     (rsmi, option))
            moutput, error = run_exec(output, "obabel -irxn -orsmi")
            self.assertEqual(self.sort(moutput.rstrip()), self.sort(result))

if __name__ == "__main__":
    unittest.main()

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

import unittest

from testbabel import run_exec, executable, log, BaseTest

class TestKekuleIsotope(BaseTest):
    """A series of tests relating to aromaticity/kekule"""

    def testSMItoCAN(self):
        """PR#1842055- bad isotope canonicalization"""
        self.canFindExecutable("babel")

        # A series of isotopamers, and their canonical forms
        self.smiles = [
            'c1ccccc1',
            '[14cH]1ccccc1',
            '[14cH]1[14cH]cccc1',
            '[14cH]1[14cH][14cH]ccc1',
            '[14cH]1[14cH][14cH][14cH]cc1',
            '[14cH]1[14cH][14cH][14cH][14cH]c1',
            '[14cH]1[14cH][14cH][14cH][14cH][14cH]1'
             ]
        self.cansmis = [
            'c1ccccc1',
            'c1[14cH]cccc1',
            '[14cH]1[14cH]cccc1',
            '[14cH]1[14cH]ccc[14cH]1',
            '[14cH]1[14cH]cc[14cH][14cH]1',
            '[14cH]1[14cH]c[14cH][14cH][14cH]1',
            '[14cH]1[14cH][14cH][14cH][14cH][14cH]1',
            ]
        for i in range(0, len(self.smiles)):
            output, error = run_exec(self.smiles[i], "babel -ismi -ocan")
            self.assertEqual(output.rstrip(), self.cansmis[i])

if __name__ == "__main__":
    testsuite = []
    for myclass in [TestKekuleIsotope]:
        suite = unittest.TestLoader().loadTestsFromTestCase(myclass)
        testsuite.append(suite)
    unittest.TextTestRunner().run(unittest.TestSuite(testsuite))

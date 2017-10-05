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

class TestReactionInChIWriter(BaseTest):
    """A series of tests relating to writing Reaction InChI"""

    def testRSMItoRINCHI(self):
        data = [
                ("C>N>O", "RInChI=1.00.1S/CH4/h1H4<>H2O/h1H2<>H3N/h1H3/d+"),
                ("O>N>C", "RInChI=1.00.1S/CH4/h1H4<>H2O/h1H2<>H3N/h1H3/d-"),
                # Example: esterification of acetic acid
                ("OCC.CC(=O)O>S(=O)(=O)(O)O>CC(=O)OCC.O", "RInChI=1.00.1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)!C2H6O/c1-2-3/h3H,2H2,1H3<>C4H8O2/c1-3-6-4(2)5/h3H2,1-2H3!H2O/h1H2<>H2O4S/c1-5(2,3)4/h(H2,1,2,3,4)/d+"),
                # Example: alkaline ring opening
                ("CC[C@]1(C)O[C@H]1C.[OH-]>>CC[C@](C)(O)[C@@H](C)O", "RInChI=1.00.1S/C6H12O/c1-4-6(3)5(2)7-6/h5H,4H2,1-3H3/t5-,6-/m0/s1!H2O/h1H2/p-1<>C6H14O2/c1-4-6(3,8)5(2)7/h5,7-8H,4H2,1-3H3/t5-,6+/m1/s1/d+"), 
                ]
        for rsmi, rinchi in data:
            output, error = run_exec('obabel -:%s -irsmi -orinchi' % rsmi)
            self.assertEqual(output.rstrip(), rinchi)

if __name__ == "__main__":
    unittest.main()
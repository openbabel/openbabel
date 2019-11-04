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

import os
import glob
import unittest

from testbabel import run_exec, executable, log, BaseTest

here = os.path.dirname(__file__)

class TestReactionInChIWriter(BaseTest):
    """A series of tests relating to writing Reaction InChI"""

    def testRSMItoRINCHI(self):
        data = [
                ("C>N>O", "RInChI=1.00.1S/CH4/h1H4<>H2O/h1H2<>H3N/h1H3/d+"),
                ("O>N>C", "RInChI=1.00.1S/CH4/h1H4<>H2O/h1H2<>H3N/h1H3/d-"),
                ("O>>C", "RInChI=1.00.1S/CH4/h1H4<>H2O/h1H2/d-"),
                # The following is assumed to be d+ by analogy with
                # the empty reaction which is d+
                ("O>>O", "RInChI=1.00.1S/H2O/h1H2<>H2O/h1H2/d+"),
                # Example: esterification of acetic acid
                ("OCC.CC(=O)O>S(=O)(=O)(O)O>CC(=O)OCC.O", "RInChI=1.00.1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)!C2H6O/c1-2-3/h3H,2H2,1H3<>C4H8O2/c1-3-6-4(2)5/h3H2,1-2H3!H2O/h1H2<>H2O4S/c1-5(2,3)4/h(H2,1,2,3,4)/d+"),
                # Example: alkaline ring opening
                ("CC[C@]1(C)O[C@H]1C.[OH-]>>CC[C@](C)(O)[C@@H](C)O", "RInChI=1.00.1S/C6H12O/c1-4-6(3)5(2)7-6/h5H,4H2,1-3H3/t5-,6-/m0/s1!H2O/h1H2/p-1<>C6H14O2/c1-4-6(3,8)5(2)7/h5,7-8H,4H2,1-3H3/t5-,6+/m1/s1/d+"),
                # Partial reactions
                (">>C1CC=C(O)CC1", "RInChI=1.00.1S/<>C6H10O/c7-6-4-2-1-3-5-6/h4,7H,1-3,5H2/d+"),
                ("C1CC=C(O)CC1>>", "RInChI=1.00.1S/<>C6H10O/c7-6-4-2-1-3-5-6/h4,7H,1-3,5H2/d-"),
                # The empty reaction
                (">>", "RInChI=1.00.1S//d+"),
                # Test 'no-structure'
                ("c1ccccc1C=C>>*", "RInChI=1.00.1S/<>C8H8/c1-2-8-6-4-3-5-7-8/h2-7H,1H2/d-/u1-0-0"),
                ("*>>C1CC=C(O)CC1", "RInChI=1.00.1S/<>C6H10O/c7-6-4-2-1-3-5-6/h4,7H,1-3,5H2/d+/u1-0-0"),
                ("O>*>C", "RInChI=1.00.1S/CH4/h1H4<>H2O/h1H2/d-/u0-0-1"),
                ("*.O>>C", "RInChI=1.00.1S/CH4/h1H4<>H2O/h1H2/d-/u0-1-0"),
                # Empty except for 'no-structures' (assumed)
                ("*>*>*", "RInChI=1.00.1S//d+/u1-1-1"),
                ]
        for eqm in [False, True]:
            for rsmi, rinchi in data:
                if eqm:
                    output, error = run_exec('obabel -:%s -ismi -orinchi -xe' % rsmi)
                    ans = rinchi.replace("/d-", "/d=").replace("/d+", "/d=")
                    self.assertEqual(output.rstrip(), ans)
                else:
                    output, error = run_exec('obabel -:%s -ismi -orinchi' % rsmi)
                    self.assertEqual(output.rstrip(), rinchi)

    def testRInChIOfficialExamples(self):
        """These test RXN to RInChI using the examples in the RInChI distrib"""
        for rxnfile in glob.glob(os.path.join(here, "rinchi", "*.rxn")):
            dirname, fname = os.path.split(rxnfile)
            output, error = run_exec('obabel %s -orinchi' % rxnfile)
            with open(os.path.join(dirname, fname.split(".")[0]+".txt")) as inp:
                ans = inp.readlines()[0]
            self.assertEqual(output.rstrip(), ans.rstrip())

if __name__ == "__main__":
    unittest.main()

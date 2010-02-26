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

class TestSym(BaseTest):
    """Base class for a series of tests relating to symmetry"""

    def testInChItoSMI(self):
        """Verify that the InChI is read correctly"""
        output, error = run_exec(self.inchi, "babel -iinchi -ocan")
        self.assertEqual(output.rstrip(), self.cansmi)

    def testSMItoInChI(self):
        """Verify that all molecules give the same InChI"""
        for smi in self.smiles:
            output, error = run_exec(smi, "babel -ismi -oinchi")
            self.assertEqual(output.rstrip(), self.inchi)

    def testSMItoCAN(self):
        """Verify that all molecules give the same cansmi"""
        for smi in self.smiles:
            output, error = run_exec(smi, "babel -ismi -ocan")
            self.assertEqual(output.rstrip(), self.cansmi)            

    def testSMIthruSDF(self):
        """Verify that roundtripping through SDF preserves stereo"""
        for smi in self.smiles:
            output, error = run_exec(smi, "babel -ismi -osdf")
            output, error = run_exec(output.rstrip(), "babel -isdf -ocan")
            self.assertEqual(output.rstrip(), self.cansmi)

    def testSMIthruXML(self):
        """Verify that roundtripping through CML preserves stereo"""
        for smi in self.smiles:
            output, error = run_exec(smi, "babel -ismi -ocml")
            output, error = run_exec(output.rstrip(), "babel -icml -ocan")
            self.assertEqual(output.rstrip(), self.cansmi)

class TestTetSym(TestSym):
    """A series of tests relating to tetrahedral symmetry"""

    def setUp(self):
        self.canFindExecutable("babel")

        # The following all represent the same molecule
        self.cansmi = "C[C@@](F)(Cl)Br"
        self.inchi = "InChI=1S/C2H3BrClF/c1-2(3,4)5/h1H3/t2-/m0/s1"
        self.smiles = [
             'C[C@@](Cl)(Br)F',
             'C[C@](Cl)(F)Br',
             'C[C@](Br)(Cl)F',
             'C[C@@](Br)(F)Cl',
             'C[C@@](F)(Cl)Br',
             'C[C@](F)(Br)Cl',
             'Cl[C@](C)(Br)F',
             'Cl[C@@](C)(F)Br',
             'Cl[C@@](Br)(C)F',
             'Cl[C@](Br)(F)C',
             'Cl[C@](F)(C)Br',
             'Cl[C@@](F)(Br)C',
             'Br[C@@](C)(Cl)F',
             'Br[C@](C)(F)Cl',
             'Br[C@](Cl)(C)F',
             'Br[C@@](Cl)(F)C',
             'Br[C@@](F)(C)Cl',
             'Br[C@](F)(Cl)C',
             'F[C@](C)(Cl)Br',
             'F[C@@](C)(Br)Cl',
             'F[C@@](Cl)(C)Br',
             'F[C@](Cl)(Br)C',
             'F[C@](Br)(C)Cl',
             'F[C@@](Br)(Cl)C'
             ]


class TestCisTransSym(TestSym):
    """A series of tests relating to cistrans symmetry"""

    def setUp(self):
        self.canFindExecutable("babel")

        # The following all represent the same molecule
        self.cansmi = "Cl/C=C/C=C\Br"
        self.inchi = "InChI=1S/C4H4BrCl/c5-3-1-2-4-6/h1-4H/b3-1-,4-2+"
        self.smiles = [
                "C(=C\C=C/Br)/Cl",
                "Cl/C=C/C=C\Br", 
                "Br/C=C\C=C\Cl",
                "C(=C\Cl)/C=C\Br",
                "C(=C\C=C\Cl)\Br",
                "C(=C\Br)\C=C\Cl"
                ]

if __name__ == "__main__":
    testsuite = []
    for myclass in [TestCisTransSym, TestTetSym]:
        suite = unittest.TestLoader().loadTestsFromTestCase(myclass)
        testsuite.append(suite)
    unittest.TextTestRunner().run(unittest.TestSuite(testsuite))

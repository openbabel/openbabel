"""Test OpenBabel executables from Python

Note: Python bindings not used

On Windows or Linux, you can run these tests at the commandline
in the build folder with:
"C:\Program Files\CMake 2.6\bin\ctest.exe" -C CTestTestfile.cmake
                                           -R pytest -VV

You could also "chdir" into build/test and run the test file directly:
python ../../../test/testdistgeom.py

In both cases, the test file is run directly from the source folder,
and so you can quickly develop the tests and try them out.
"""

import unittest

from testbabel import run_exec, executable, log, BaseTest

class TestDistanceGeomStereo(BaseTest):
    """A series of tests relating to 3D stereo"""

    def testSMItoSMI(self):
        """
        Some initial tests based on previous geometry stereo bugs
        (i.e., OBBuilder can't handle correctly)
        """
        self.canFindExecutable("obabel")

        # A series of aromatic strings, which should convert to themselves
        self.smiles = [
            'c1ccccc1',  # benzene
            'C#C', # triple bond
            'CC=CC', # butene unspecified
            'C/C=C\\C',  # Z-butene
            'C/C=C/C',  # E-butene
            'NC(Br)(O)C',
            'N[C@](Br)(O)C',
            'N[C@@](Br)(O)C',
            'CCC[C@@H]([C@H](CC(C)C)C)C',
            'C1CC[C@H]2[C@@H](C1)CCCC2',  # cis-decalin
            'C1CC[C@@H]2[C@@H](C1)CCCC2',  # trans-decalin
            '[C@H]1(NC[C@H]2[C@H]1N2)OC',
            'Clc1cccc(Cl)c1\C=N\NC(=O)c1cccs1',
            'O=C1NC(=S)S\C1=C/c1ccco1',
            'S=C1NC(=O)/C(=C/c2ccco2)/S1',
            'O=C1NC(=S)N\C1=C\c1ccncc1',
            'S=C1NC(=O)C(=C)N1',
            'CC(=O)N\N=C\c1ccncc1',
            'N/N=c/1\sc2c(n1C)cccc2',
            'OCCN/C=C\\1/C(=NN(C1=O)c1ccccc1)C',
            'Cc1ccc(o1)/C=C/C=O',
            # disabled to make test run faster:
            #'CCCNC1=C(C)C(=O)C2=C(C1=O)[C@@H](COC(=O)N)[C@]1(N2C[C@H]2[C@H]1N2)OC',
            #'CN([C@H]1C(=O)C(=C([C@]2([C@@H]1C[C@@H]1Cc3c(C(=C1C2=O)O)c(O)ccc3N(C)C)O)O)C(=O)N)C',
            #'CN([C@@H]1C(=O)C(=C([C@@]2([C@H]1C[C@@H]1C(=C(O)c3c([C@@]1(C)O)c(Cl)ccc3O)C2=O)O)O)C(=O)N)C',
            #'C[C@@H](CC(=O)OC[C@@]12CC[C@H]3[C@@]([C@@H]2C[C@@H](O1)C1=CC(=O)O[C@H]1O)(C)CC[C@@H]1[C@]3(C)CCCC1(C)C)O',
            #'CC(=O)OC[C@@]12CC[C@H]3[C@@]([C@@H]2C[C@H](O1)C1=CC(=O)O[C@H]1O)(C)CC[C@@H]1[C@]3(C)CCCC1(C)C'
            ]
        for smi in self.smiles:
            # generate a canonical SMILES in case the ordering changes
            canSMI, error = run_exec(smi, "obabel -ismi -ocan")
            # generate a mol2 (any 3D format without implicit hydrogens)
            mol2, error = run_exec(smi, "obabel -ismi -osdf --gen3d dg")
            # now check if it matches the previous canonical SMILES
            output, error = run_exec(mol2, "obabel -isdf -ocan")

            self.assertEqual(output.split('\t')[0].rstrip(), canSMI.rstrip())


if __name__ == "__main__":
    testsuite = []
    for myclass in [TestDistanceGeomStereo]:
        suite = unittest.TestLoader().loadTestsFromTestCase(myclass)
        testsuite.append(suite)
    unittest.TextTestRunner().run(unittest.TestSuite(testsuite))

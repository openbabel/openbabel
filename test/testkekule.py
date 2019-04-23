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

from testbabel import run_exec, executable, BaseTest

class TestKekuleAssignment(BaseTest):
    """A series of tests relating to aromaticity/kekule"""

    def testSMItoSMI(self):
        """
        PR#2705497 aromatic - kekule conversion issue
        PR#1445453 SMILES aromaticity fails on 4-valent N+ atoms
        PR#1814248 Aromaticity munged by SMILES input
        PR#1761638 Error in Aromaticity / Kekulize
        PR#2948661 Trunk fails aromaticity
        """
        self.canFindExecutable("obabel")

        # A series of aromatic strings, which should convert to themselves
        self.smiles = [
            'c12-c3c(cc(N)cc3)Cc1cccc2',
            'c1(=O)n(c2c(c(=O)o1)cccc2)CC(=O)OCC',
            'c1n[nH]c(=S)[nH]1',
            'O=c1[nH]ccc2nc3oc4ccccc4c(=O)c3cc12',
            'c1nc2sccn2c1',
            'c1[nH+]cnc2[nH]cnc12',
            'c1onc(c2ccccc2Cl)c1',
            'c1ccc2[nH]c3ccc4cc[nH+]cc4c3c2c1',
            '[nH]1c2ccccc2c2c3C(=O)NCc3c3c4ccccc4[nH]c3c12',
            'c1c(C)c2C=c3[n-]c(=Cc4[nH]c(C=c5[n-]c(=Cc1[nH]2)c(C)c5C=C)c(C)c4CCC(=O)O)c(CCC(=O)O)c3',
            'C1=C2CCC(=Cc3ccc([nH]3)C=c3ccc(=Cc4ccc1[nH]4)[nH]3)N2',
            'c1(NC(=O)C2CC2)nc2-c3c(cccc3)CCc2cn1',
            'O=C1N(CCCC)C(=O)NC2C1C1N(N2)CCN1',
            'Cn1cccnc1=O',
            'O=c1n(C)c(=O)nc2-c1c1n([nH]2)cc[nH]1',
            'Cn1ccn2c1nc1c2c(=O)n(C)c(=O)n1C'
            ]
        for i in range(0, len(self.smiles)):
            output, error = run_exec(self.smiles[i], "obabel -ismi -osmi")
            self.assertEqual(output.rstrip(), self.smiles[i])

class TestKekuleIsotope(BaseTest):
    """A series of tests relating to aromaticity/kekule"""

    def testSMItoCAN(self):
        """PR#1842055- bad isotope canonicalization"""
        self.canFindExecutable("obabel")

        # A series of isotopamers, and their canonical forms
        self.smiles = [
            'c1ccccc1',
            'c1[14cH]cccc1',
            'c1[14cH][14cH]ccc1',
            'c1[14cH][14cH][14cH]cc1',
            'c1[14cH][14cH][14cH][14cH]c1',
            'c1[14cH][14cH][14cH][14cH][14cH]1',
            '[14cH]1[14cH][14cH][14cH][14cH][14cH]1',
            ]
        self.cansmis = [
            'c1ccccc1',
            '[14cH]1ccccc1',
            '[14cH]1[14cH]cccc1',
            '[14cH]1[14cH]ccc[14cH]1',
            '[14cH]1[14cH][14cH]cc[14cH]1',
            '[14cH]1[14cH][14cH]c[14cH][14cH]1',
            '[14cH]1[14cH][14cH][14cH][14cH][14cH]1',
            '[14cH]1[14cH][14cH][14cH][14cH][14cH]1',
            ]
        for i in range(0, len(self.smiles)):
            output, error = run_exec(self.smiles[i], "obabel -ismi -ocan")
            self.assertEqual(output.rstrip(), self.cansmis[i])

class TestKekuleCrashers(BaseTest):
    """A series of tests which caused crashes"""

    def testXYZtoXYZ(self):
        """PR#2956135- crash in kekulize"""
        self.canFindExecutable("obabel")

        # A series of isotopamers, and their canonical forms
        self.xyz = """39
crash.gamout
C         -0.31501       -0.05904        0.00332
C          0.47846        1.04480        0.28483
N          1.83248        1.00566        0.38129
C          2.42024       -0.19701        0.19943
C          1.71953       -1.35983       -0.07408
C          0.33143       -1.28624       -0.17253
H         -0.24386       -2.18519       -0.38442
H          2.23907       -2.30258       -0.20635
H          3.50341       -0.20251        0.28526
H          0.06258        2.02989        0.46235
C         -1.79310       -0.00135       -0.09779
O         -2.46156       -1.02575       -0.18756
N         -2.41033        1.21816       -0.10797
H         -1.94649        2.11169       -0.16198
H         -3.41687        1.20930       -0.22381
C          0.26924        0.47947       -3.35313
C         -0.44373       -0.03140       -4.37204
H         -0.48212        0.45752       -5.34214
H         -1.00864       -0.95087       -4.26294
C          1.00998        1.73774       -3.59280
O          1.09874        2.34143       -4.64607
O          1.61155        2.19132       -2.48091
H          2.04906        3.01525       -2.77997
C          0.34525       -0.16868       -2.00668
H         -0.05656        0.49552       -1.23462
H          1.38253       -0.41344       -1.75510
H         -0.22965       -1.10058       -1.97049
C          0.26203        0.28311        3.26498
C          1.38529        0.78393        2.72190
H          1.56449        0.74957        1.65029
H          2.15564        1.24022        3.33413
C         -0.73494       -0.32012        2.35309
O         -0.68045       -0.35957        1.13769
O         -1.76945       -0.85942        3.01857
H         -2.33520       -1.21968        2.30453
C         -0.01692        0.31604        4.73465
H         -0.11977       -0.69915        5.13149
H         -0.93979        0.86887        4.93917
H          0.78936        0.80651        5.29109
"""
        output, error = run_exec(self.xyz, "obabel -ixyz -oxyz")
        self.assertConverted(error, 1)

if __name__ == "__main__":
    testsuite = []
    for myclass in [TestKekuleAssignment, TestKekuleIsotope, TestKekuleCrashers]:
        suite = unittest.TestLoader().loadTestsFromTestCase(myclass)
        testsuite.append(suite)
    unittest.TextTestRunner().run(unittest.TestSuite(testsuite))

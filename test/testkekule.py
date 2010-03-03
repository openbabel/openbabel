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

    def testXYZtoXYZ(self):
        """PR#2956135- crash in kekulize"""
        self.canFindExecutable("babel")

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
        output, error = run_exec(self.xyz, "babel -ixyz -oxyz")
        self.assertConverted(error, 1)

if __name__ == "__main__":
    testsuite = []
    for myclass in [TestKekuleIsotope]:
        suite = unittest.TestLoader().loadTestsFromTestCase(myclass)
        testsuite.append(suite)
    unittest.TextTestRunner().run(unittest.TestSuite(testsuite))

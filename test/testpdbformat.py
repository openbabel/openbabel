"""Test OpenBabel executables from Python

Note: Python bindings not used

On Windows or Linux, you can run these tests at the commandline
in the build folder with:
"C:\Program Files\CMake 2.6\bin\ctest.exe" -C CTestTestfile.cmake
                                           -R pytest -VV

You could also "chdir" into build/test and run the test file directly:
python ../../test/testpdbformat.py


In both cases, the test file is run directly from the source folder,
and so you can quickly develop the tests and try them out.
"""

import unittest

from testbabel import run_exec, executable, BaseTest

class TestPDBFormat(BaseTest):
    """A series of tests relating to PDB"""

    def testInsertionCodes(self):
        """
        Testing a PDB entry with insertion codes to distinguish residues
        upon conversion to FASTA.
        """
        self.canFindExecutable("obabel")

        self.entryPDBwithInsertioncodes="""ATOM    406  N   VAL L  29      58.041  17.797  48.254  1.00  0.00           N
ATOM    407  CA  VAL L  29      57.124  18.088  47.170  1.00  0.00           C
ATOM    408  C   VAL L  29      55.739  17.571  47.538  1.00  0.00           C
ATOM    409  O   VAL L  29      55.535  16.362  47.550  1.00  0.00           O
ATOM    410  CB  VAL L  29      57.580  17.456  45.842  1.00  0.00           C
ATOM    411  CG1 VAL L  29      56.571  17.743  44.741  1.00  0.00           C
ATOM    412  CG2 VAL L  29      58.957  17.973  45.450  1.00  0.00           C
ATOM    413  H   VAL L  29      58.603  16.959  48.212  1.00  0.00           H
ATOM    414  HA  VAL L  29      57.012  19.163  47.024  1.00  0.00           H
ATOM    415  HB  VAL L  29      57.674  16.378  45.977  1.00  0.00           H
ATOM    416 1HG1 VAL L  29      56.909  17.289  43.809  1.00  0.00           H
ATOM    417 2HG1 VAL L  29      55.603  17.327  45.016  1.00  0.00           H
ATOM    418 3HG1 VAL L  29      56.479  18.821  44.604  1.00  0.00           H
ATOM    419 1HG2 VAL L  29      59.263  17.515  44.510  1.00  0.00           H
ATOM    420 2HG2 VAL L  29      58.917  19.055  45.331  1.00  0.00           H
ATOM    421 3HG2 VAL L  29      59.676  17.719  46.229  1.00  0.00           H
ATOM    422  N   SER L  30      54.838  18.500  47.837  1.00  0.00           N
ATOM    423  CA  SER L  30      53.494  18.162  48.273  1.00  0.00           C
ATOM    424  C   SER L  30      52.725  17.364  47.221  1.00  0.00           C
ATOM    425  O   SER L  30      52.723  17.697  46.056  1.00  0.00           O
ATOM    426  CB  SER L  30      52.734  19.429  48.610  1.00  0.00           C
ATOM    427  OG  SER L  30      51.403  19.143  48.941  1.00  0.00           O
ATOM    428  H   SER L  30      55.100  19.472  47.757  1.00  0.00           H
ATOM    429  HA  SER L  30      53.471  17.585  49.199  1.00  0.00           H
ATOM    430 1HB  SER L  30      53.219  19.934  49.445  1.00  0.00           H
ATOM    431 2HB  SER L  30      52.761  20.107  47.758  1.00  0.00           H
ATOM    432  HG  SER L  30      50.919  19.965  48.828  1.00  0.00           H
ATOM    433  N   SER L  30A     52.170  16.303  47.698  1.00  0.00           N
ATOM    434  CA  SER L  30A     51.329  15.409  46.920  1.00  0.00           C
ATOM    435  C   SER L  30A     52.015  14.812  45.685  1.00  0.00           C
ATOM    436  O   SER L  30A     51.350  14.366  44.764  1.00  0.00           O
ATOM    437  CB  SER L  30A     50.082  16.156  46.488  1.00  0.00           C
ATOM    438  OG  SER L  30A     49.348  16.592  47.599  1.00  0.00           O
ATOM    439  H   SER L  30A     52.421  16.046  48.642  1.00  0.00           H
ATOM    440  HA  SER L  30A     50.943  14.567  47.497  1.00  0.00           H
ATOM    441 1HB  SER L  30A     50.364  17.013  45.876  1.00  0.00           H
ATOM    442 2HB  SER L  30A     49.463  15.505  45.873  1.00  0.00           H
ATOM    443  HG  SER L  30A     49.931  17.176  48.090  1.00  0.00           H
ATOM    444  N   SER L  31      53.347  14.792  45.683  1.00  0.00           N
ATOM    445  CA  SER L  31      54.094  14.259  44.549  1.00  0.00           C
ATOM    446  C   SER L  31      53.734  14.959  43.242  1.00  0.00           C
ATOM    447  O   SER L  31      53.703  14.356  42.179  1.00  0.00           O
ATOM    448  CB  SER L  31      53.835  12.771  44.418  1.00  0.00           C
ATOM    449  OG  SER L  31      54.240  12.087  45.572  1.00  0.00           O
ATOM    450  H   SER L  31      53.852  15.150  46.480  1.00  0.00           H
ATOM    451  HA  SER L  31      55.175  14.292  44.689  1.00  0.00           H
ATOM    452 1HB  SER L  31      52.774  12.600  44.243  1.00  0.00           H
ATOM    453 2HB  SER L  31      54.375  12.383  43.555  1.00  0.00           H
ATOM    454  HG  SER L  31      53.773  11.248  45.560  1.00  0.00           H
ATOM    455  N   TYR L  32      53.460  16.259  43.402  1.00  0.00           N
ATOM    456  CA  TYR L  32      53.176  17.161  42.301  1.00  0.00           C
ATOM    457  C   TYR L  32      54.489  17.641  41.668  1.00  0.00           C
ATOM    458  O   TYR L  32      54.910  18.762  41.892  1.00  0.00           O
ATOM    459  CB  TYR L  32      52.342  18.352  42.780  1.00  0.00           C
ATOM    460  CG  TYR L  32      50.880  18.031  42.990  1.00  0.00           C
ATOM    461  CD1 TYR L  32      50.294  16.936  42.371  1.00  0.00           C
ATOM    462  CD2 TYR L  32      50.089  18.824  43.807  1.00  0.00           C
ATOM    463  CE1 TYR L  32      48.958  16.639  42.559  1.00  0.00           C
ATOM    464  CE2 TYR L  32      48.751  18.535  44.002  1.00  0.00           C
ATOM    465  CZ  TYR L  32      48.190  17.441  43.376  1.00  0.00           C
ATOM    466  OH  TYR L  32      46.859  17.150  43.569  1.00  0.00           O
ATOM    467  H   TYR L  32      53.456  16.618  44.347  1.00  0.00           H
ATOM    468  HA  TYR L  32      52.651  16.625  41.509  1.00  0.00           H
ATOM    469 1HB  TYR L  32      52.778  18.693  43.721  1.00  0.00           H
ATOM    470 2HB  TYR L  32      52.439  19.136  42.030  1.00  0.00           H
ATOM    471  HD1 TYR L  32      50.908  16.305  41.727  1.00  0.00           H
ATOM    472  HD2 TYR L  32      50.537  19.687  44.299  1.00  0.00           H
ATOM    473  HE1 TYR L  32      48.512  15.775  42.066  1.00  0.00           H
ATOM    474  HE2 TYR L  32      48.145  19.172  44.648  1.00  0.00           H
ATOM    475  HH  TYR L  32      46.462  17.658  44.280  1.00  0.00           H
"""
        output, error = run_exec(self.entryPDBwithInsertioncodes,
                                     "obabel -ipdb -ofasta")
        self.assertEqual(output.rstrip().rsplit("\n",1)[1], "VSSSY")

if __name__ == "__main__":
    testsuite = []
    for myclass in [TestPDBFormat]:
        suite = unittest.TestLoader().loadTestsFromTestCase(myclass)
        testsuite.append(suite)
    unittest.TextTestRunner().run(unittest.TestSuite(testsuite))

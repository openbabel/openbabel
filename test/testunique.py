"""Test OpenBabel executables from Python

Note: Python bindings not used

On Windows or Linux, you can run these tests at the commandline
in the build folder with:
"C:\Program Files\CMake 2.6\bin\ctest.exe" -C CTestTestfile.cmake
                                           -R pytest -VV

You could also "chdir" into build/test and run the test file directly:
python ../../../test/testunique.py

In both cases, the test file is run directly from the source folder,
and so you can quickly develop the tests and try them out.
"""

import unittest

from testbabel import run_exec, BaseTest

class TestUnique(BaseTest):
    """A series of tests relating to obabel --unique"""

    def setUp(self):
        self.canFindExecutable("obabel")
        self.smiles = """C	methane
COC	dimethyl ether
OC(C)CC	2-butanol
O[C@H](C)CC	R-2-butanol
O[C@@H](C)CC	S-2-butanol
O([CH3])[CH3]	DME
FC=CF	difluoroethene
F/C=C/F	trans-difluoroethene
F/C=C\F	cis-difluoroethene
C(N)N	diamino methane
C(N)[NH3+]	protonated diamino methane
CCO	ethanol
CC[18O][2H]	ethanol-18OD
C([2H])([2H])([2H])[2H]	deuteromethane"""

    def testFindDups(self):
        """Look for duplicates using --unique"""

        params = [("", 13), ("/formula", 5),
                  ("/connect", 6), ("/nostereo", 9),
                  ("/nosp3", 11), ("/noEZ", 11),
                  ("/nochg", 12), ("/noiso", 11),
                  ("cansmi", 13), ("cansmiNS", 7)]

        for param in params:
            output, error = run_exec(self.smiles,
                                     "obabel -ismi -osmi --unique %s" % param[0])
            self.assertConverted(error, param[1])

if __name__ == "__main__":
    unittest.main()

"""Test OpenBabel executables from Python

Note: Python bindings not used

On Windows or Linux, you can run these tests at the commandline
in the build folder with:
"C:\Program Files\CMake 2.6\bin\ctest.exe" -C CTestTestfile.cmake
                                           -R pytest -VV

You could also "chdir" into build/test and run the test file directly:
python ../../../test/testfastsearch.py

In both cases, the test file is run directly from the source folder,
and so you can quickly develop the tests and try them out.
"""

import unittest

from testbabel import run_exec, BaseTest

class TestSym(BaseTest):
    """A series of tests relating to fastsearch functionality"""

    def setUp(self):
        self.canFindExecutable("obabel")

    def testSingleHit(self):
        """PR#2955101 - Difficulty reading from a fastsearch index"""

        smiles = """C12(C(N(C(=O)C)c3c2cccc3)=O)Nc2c(ccc(c2N1)OCCCC)OCCCC
n1c([nH]c(cc1c1ccccc1)=O)c1ccc(cc1)Br
n1c(nc2c(c1N(C)C)cccc2)c1c(O)cccc1
C1(/[CH]2[CH]3\C(=C4/CC(C)(C)NC(C4)(C)C)C=C[CH]3[CH]1C=C2)=C1/CC(C)(C)NC(C1)(C)C
n1c(c2ccc(C(=O)O)cc2)ccc(c1)CCCCC
N1(C(CN(CC1=O)C(=O)C1CCCCC1)=O)CCc1ccccc1
S(N1[CH](c2ccccc2C=C1)C#N)(c1ccc(cc1)C)(=O)=O
c12c(c(OC)c3c(c1OC)occ3)ccc(o2)=O
c12c(O[CH](C1=O)C(C)C)cc1c(c2)ccc(=O)o1
c12[C]3([C@H]4([N@@](CCc1c1ccccc1[nH]2)C[C@H](C=C4CC)C3))C(=O)OC"""

        outputfile = open("ten.smi", "w")
        outputfile.write(smiles)
        outputfile.close()

        output, error = run_exec("obabel ten.smi -O ten.fs")
        self.canFindFile("ten.fs")
        self.assertConverted(error, 10)

        query = "Nc2nc(c1ccccc1)nc3ccccc23"
        output, error = run_exec("obabel ten.fs -ifs -s %s -osmi" % query)
        self.assertConverted(error, 1)

        output, error = run_exec("obabel ten.fs -ifs -s %s -at 0.5 -aa -osmi" % query)
        self.assertConverted(error, 1)



if __name__ == "__main__":
    unittest.main()

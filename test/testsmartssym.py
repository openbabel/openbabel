"""Test OpenBabel executables from Python

Note: Python bindings not used

On Windows or Linux, you can run these tests at the commandline
in the build folder with:
"C:\Program Files\CMake 2.6\bin\ctest.exe" -C CTestTestfile.cmake
                                           -R pytest -VV

You could also "chdir" into build/test and run the test file directly:
python ../../../test/testsmartssym.py

In both cases, the test file is run directly from the source folder,
and so you can quickly develop the tests and try them out.
"""

import unittest

from testbabel import run_exec, BaseTest

def checkmatch(query, molecules):
    result = []
    for smi in molecules:
        output, error = run_exec("obabel -:%s -s%s -osmi" % (smi, query))
        result.append(output.strip() != "")
    return result

def fastcheckmatch(query, molecules):
    """May fail where Open Babel does not output the input query, e.g.
    [C@@]([H])(Br)(Cl)I is output as [C@@H](Br)(Cl)I"""
    output, error = run_exec("\n".join(molecules), "obabel -ismi -s%s -osmi" % query)
    converted = [x.rstrip() for x in output.split("\n")]
    results = [smi in converted for smi in molecules]
    return results

class TestSmartsSym(BaseTest):
    """Base class for a series of tests relating to symmetry"""

    def testSelfMatch(self):
        """Verify that a molecule matches itself"""
        data = [
                '[C@@](F)(Br)(Cl)I',
                '[C@](F)(Br)(Cl)I',
                'F[C@](Br)(Cl)I',
                '[C@H](Br)(Cl)I',
                'Br[C@H](Cl)I',
                '[C@]1(Br)(Cl)NC1',
                '[C@@]1(Br)(Cl)NC1',
                'Br[C@]1(Cl)NC1',
                'C1N[C@]1(Cl)Br',
                'F[C@]1(Br)N[C@]1(Br)Cl',
                '[C@H]1(Cl)NC1'
            ]
        for smi in data:
            output, error = run_exec("obabel -:%s -s%s -osmi" % (smi, smi))
            self.assertEqual(output.rstrip(), smi)

    def testTetStereo(self):
        data = ['[C@@](F)(Br)(Cl)I',
                '[C@](F)(Br)(Cl)I',
                'F[C@](Br)(Cl)I',
                'F[C@@](Br)(Cl)I',
                'C(F)(Br)(Cl)I',
                'FC(Br)(Cl)I']
        self.assertEqual(fastcheckmatch(data[0], data[0:6]),
                         [True, False, False, True, False, False])
        self.assertEqual(fastcheckmatch(data[2], data[0:6]),
                         [False, True, True, False, False, False])
        self.assertEqual(fastcheckmatch(data[4], data[0:6]), [True]*6)

    def testTetStereoImplicitH(self):
        data = ['[C@H](Br)(Cl)I',
                '[C@@H](Br)(Cl)I',
                'Br[C@H](Cl)I',
                'Br[C@@H](Cl)I',
                'BrC(Cl)I',
                'BrC([H])(Cl)I',
                'Br[C@@]([H])(Cl)I'
                ]
        self.assertEqual(checkmatch(data[0], data[0:7]),
                         [True, False, False, True, False, False, True])
        self.assertEqual(checkmatch(data[2], data[0:7]),
                         [False, True, True, False, False, False, False])
        self.assertEqual(checkmatch(data[4], data[0:7]), [True]*7)
        self.assertEqual(checkmatch(data[6], data[0:7]),
                         [True, False, False, True, False, False, True])

    def testRingClosures(self):
        data = ['[C@]1(Br)(Cl)NC1',
                '[C@@]1(Br)(Cl)NC1',
                'Br[C@]1(Cl)NC1',
                'Br[C@@]1(Cl)NC1',
                'C1N[C@]1(Cl)Br',
                'C1NC1(Cl)Br']
        self.assertEqual(fastcheckmatch(data[0], data[0:6]),
                         [True, False, False, True, True, False])
        self.assertEqual(fastcheckmatch(data[2], data[0:6]),
                         [False, True, True, False, False, False])
        self.assertEqual(fastcheckmatch(data[5], data[0:6]), [True]*6)


if __name__ == "__main__":
    testsuite = []
    allclasses = [TestSmartsSym]
    for myclass in allclasses:
        suite = unittest.TestLoader().loadTestsFromTestCase(myclass)
        testsuite.append(suite)
    unittest.TextTestRunner().run(unittest.TestSuite(testsuite))

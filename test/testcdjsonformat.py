"""Test ChemDoodle JSON format using the OpenBabel Python bindings

On Windows or Linux, you can run these tests at the commandline
in the build folder with:
"C:\Program Files\CMake 2.6\bin\ctest.exe" -C CTestTestfile.cmake
                                           -R pybindtest -VV

You could also "chdir" into build/test and run the test file directly:
python ../../test/testcdjsonformat.py

In both cases, the test file is run directly from the source folder,
and so you can quickly develop the tests and try them out.
"""

import json
import os
import unittest
from testbindings import PybelWrapper, pybel


filedir = os.path.join(os.path.dirname(__file__), 'cdjson')


class TestCdJsonFormat(PybelWrapper):
    """Test ChemDoodle JSON format."""

    # def test_read_invalid(self):
    #     """Test reading a file that is not valid JSON."""
    #     # Expected error log: JSON parse error at offset 25: Invalid escape character in string.
    #     mols = list(pybel.readfile("cdjson", os.path.join(filedir, 'invalid.json')))
    #     self.assertEqual(mols, [])

    # def test_read_not_object(self):
    #     """Test reading a JSON file that doesn't have a root object."""
    #     # Expected error log: JSON file should be a single object
    #     mols = list(pybel.readfile("cdjson", os.path.join(filedir, 'notobject.json')))
    #     self.assertEqual(mols, [])

    def test_read_empty(self):
        """Test reading a file with an empty molecules array."""
        mols = list(pybel.readfile("cdjson", os.path.join(filedir, 'empty.json')))
        self.assertEqual(mols, [])

    def test_read_proton(self):
        """Test reading a file with a single hydrogen atom (and no bonds)."""
        mols = list(pybel.readfile("cdjson", os.path.join(filedir, 'proton.json')))
        self.assertEqual(len(mols), 1)
        self.assertEqual(mols[0].OBMol.NumAtoms(), 1)
        self.assertEqual(mols[0].OBMol.GetAtom(1).GetFormalCharge(), 1)
        self.assertEqual(mols[0].OBMol.NumBonds(), 0)

    def test_read_multiple(self):
        """Test reading a file with multiple molecules."""
        mols = list(pybel.readfile("cdjson", os.path.join(filedir, 'butane.json')))
        self.assertEqual(len(mols), 2)
        self.assertEqual(mols[0].OBMol.NumAtoms(), 4)
        self.assertEqual(mols[1].OBMol.NumAtoms(), 4)

    def test_read_atoms(self):
        """Test reading atoms."""
        mols = list(pybel.readfile("cdjson", os.path.join(filedir, 'butane.json')))
        self.assertEqual([a.atomicnum for a in mols[0].atoms], [6, 6, 6, 6])
        self.assertEqual([a.formalcharge for a in mols[0].atoms], [0, 0, 0, 0])

    def test_read_bonds(self):
        """Test reading bonds."""
        mols = list(pybel.readfile("cdjson", os.path.join(filedir, 'butane.json')))
        self.assertEqual(mols[0].OBMol.NumBonds(), 3)
        self.assertEqual([mols[0].OBMol.GetBond(i).GetBondOrder() for i in range(0, 3)], [1, 1, 1])

    def test_write_atoms(self):
        """Test writing atoms."""
        mols = list(pybel.readfile("cdjson", os.path.join(filedir, 'butane.json')))
        output = json.loads(mols[0].write('cdjson'))
        self.assertEqual(len(output['m'][0]['a']), 4)
        for a in output['m'][0]['a']:
            self.assertTrue('x' in a)
            self.assertTrue('y' in a)

    def test_write_bonds(self):
        """Test writing bonds."""
        mols = list(pybel.readfile("cdjson", os.path.join(filedir, 'butane.json')))
        output = json.loads(mols[0].write('cdjson'))
        self.assertEqual(output['m'][0]['b'], [{'b': 0, 'e': 1}, {'b': 0, 'e': 2}, {'b': 0, 'e': 3}])

    def test_write_minified(self):
        """Test writing minified output."""
        mols = list(pybel.readfile("cdjson", os.path.join(filedir, 'butane.json')))
        output = mols[0].write('cdjson', opt={'m': None})
        self.assertTrue('\n' not in output)
        output = mols[0].write('cdjson')
        self.assertTrue('\n' in output)


if __name__ == "__main__":
    unittest.main()

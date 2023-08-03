"""Test PubChem JSON format using the OpenBabel Python bindings

On Windows or Linux, you can run these tests at the commandline
in the build folder with:
"C:\Program Files\CMake 2.6\bin\ctest.exe" -C CTestTestfile.cmake
                                           -R pybindtest -VV

You could also "chdir" into build/test and run the test file directly:
python ../../test/testpcjsonformat.py

In both cases, the test file is run directly from the source folder,
and so you can quickly develop the tests and try them out.
"""

import json
import os
import unittest
from testbindings import PybelWrapper, pybel


filedir = os.path.join(os.path.dirname(__file__), 'pcjson')


class TestPcJsonFormat(PybelWrapper):
    """Test PubChem JSON format."""

    # def test_read_invalid(self):
    #     """Test reading a file that is not valid JSON."""
    #     # Expected error log: JSON parse error at offset 37: Invalid escape character in string.
    #     mols = list(pybel.readfile("pcjson", os.path.join(filedir, 'invalid.json')))
    #     self.assertEqual(mols, [])
    #
    # def test_read_not_object(self):
    #     """Test reading a JSON file that doesn't have a root object."""
    #     # Expected error log: JSON file should be a single object
    #     mols = list(pybel.readfile("pcjson", os.path.join(filedir, 'notobject.json')))
    #     self.assertEqual(mols, [])

    def test_read_empty(self):
        """Test reading a file with an empty compounds array."""
        mols = list(pybel.readfile("pcjson", os.path.join(filedir, 'empty.json')))
        self.assertEqual(mols, [])

    def test_read_proton(self):
        """Test reading a file with a single hydrogen atom (and no bonds)."""
        mols = list(pybel.readfile("pcjson", os.path.join(filedir, 'CID_1038_2D.json')))
        self.assertEqual(len(mols), 1)
        self.assertEqual(mols[0].OBMol.NumAtoms(), 1)
        self.assertEqual(mols[0].OBMol.GetAtom(1).GetFormalCharge(), 1)
        self.assertEqual(mols[0].OBMol.NumBonds(), 0)

    def test_read_sid(self):
        """Test reading a PubChem substance SID."""
        mols = list(pybel.readfile("pcjson", os.path.join(filedir, 'SID_348356775_2D.json')))
        self.assertEqual(len(mols), 1)
        self.assertEqual(mols[0].data['sid'], '348356775')
        self.assertEqual(mols[0].title, '348356775')

    def test_read_cid(self):
        """Test reading a PubChem compound CID."""
        mols = list(pybel.readfile("pcjson", os.path.join(filedir, 'CID_2244_2D.json')))
        self.assertEqual(len(mols), 1)
        self.assertEqual(mols[0].data['cid'], '2244')
        self.assertEqual(mols[0].title, '2244')

    def test_read_atoms(self):
        """Test reading atoms for a PubChem compound."""
        mols = list(pybel.readfile("pcjson", os.path.join(filedir, 'CID_6857552_2D.json')))
        self.assertEqual(len(mols), 1)
        self.assertEqual(len(mols[0].atoms), 14)
        self.assertEqual([a.atomicnum for a in mols[0].atoms], [8, 8, 8, 7, 6, 6, 6, 1, 1, 1, 1, 1, 1, 1])
        self.assertEqual([a.formalcharge for a in mols[0].atoms], [0, -1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

    def test_read_bonds(self):
        """Test reading bonds for a PubChem compound."""
        mols = list(pybel.readfile("pcjson", os.path.join(filedir, 'CID_6857552_2D.json')))
        self.assertEqual(len(mols), 1)
        self.assertEqual(mols[0].OBMol.NumBonds(), 13)
        self.assertEqual([mols[0].OBMol.GetBond(i).GetBondOrder() for i in range(0, 13)], [1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1])

    def test_write_cid(self):
        """Test writing a PubChem compound CID."""
        mols = list(pybel.readfile("pcjson", os.path.join(filedir, 'CID_2244_2D.json')))
        output = json.loads(mols[0].write('pcjson'))
        self.assertEqual(len(output['PC_Compounds']), 1)
        self.assertTrue('id' in output['PC_Compounds'][0])
        self.assertEqual(output['PC_Compounds'][0]['id']['id']['cid'], '2244')

    def test_write_atoms(self):
        """Test writing atoms for a PubChem compound."""
        mols = list(pybel.readfile("pcjson", os.path.join(filedir, 'CID_6857552_2D.json')))
        output = json.loads(mols[0].write('pcjson'))
        self.assertEqual(output['PC_Compounds'][0]['atoms']['aids'], list(range(1, 15)))
        self.assertEqual(output['PC_Compounds'][0]['atoms']['element'], [8, 8, 8, 7, 6, 6, 6, 1, 1, 1, 1, 1, 1, 1])
        self.assertEqual(output['PC_Compounds'][0]['atoms']['charge'][0], {'aid': 2, 'value': -1})
        self.assertEqual(output['PC_Compounds'][0]['atoms']['charge'][1], {'aid': 4, 'value': 1})

    def test_write_minified(self):
        """Test writing minified output."""
        mols = list(pybel.readfile("pcjson", os.path.join(filedir, 'CID_6857552_2D.json')))
        output = mols[0].write('pcjson', opt={'m': None})
        self.assertTrue('\n' not in output)
        output = mols[0].write('pcjson')
        self.assertTrue('\n' in output)

    def test_write_complex_bonds(self):
        """Test writing complex bonds."""
        mols = list(pybel.readfile("pcjson", os.path.join(filedir, 'CID_368849_2D.json'), opt={'s': None}))
        output = json.loads(mols[0].write('pcjson', opt={'w': None}))
        self.assertEqual(output['PC_Compounds'][0]['bonds']['order'][:4], [6, 6, 6, 6])

    def test_write_charge(self):
        """Test writing molecule charge."""
        mols = list(pybel.readfile("pcjson", os.path.join(filedir, 'CID_1038_2D.json')))
        output = json.loads(mols[0].write('pcjson'))
        self.assertEqual(output['PC_Compounds'][0]['charge'], 1)

    def test_write_stereo_tetrahedral(self):
        """Test writing tetrahedral stereochemistry."""
        mols = list(pybel.readfile("pcjson", os.path.join(filedir, 'CID_6137_2D.json')))
        output = json.loads(mols[0].write('pcjson'))
        self.assertEqual(output['PC_Compounds'][0]['stereo'], [
            {'tetrahedral': {'above': 12, 'below': 5, 'bottom': 8, 'center': 6, 'parity': 1, 'top': 4, 'type': 1}}
        ])

    def test_write_stereo_cistrans(self):
        """Test writing cis-trans stereochemistry."""
        # (Superfluous?) parity value is not set to same/opposite (1/2). unknown (255) is valid though.
        # Cis
        mols = list(pybel.readfile("pcjson", os.path.join(filedir, 'CID_643833_2D.json')))
        output = json.loads(mols[0].write('pcjson'))
        self.assertEqual(output['PC_Compounds'][0]['stereo'], [
            {'planar': {'type': 1, 'ltop': 1, 'left': 3, 'right': 4, 'rbottom': 6, 'lbottom': 5, 'rtop': 2, 'parity': 255}}
        ])
        # Trans
        mols = list(pybel.readfile("pcjson", os.path.join(filedir, 'CID_638186_2D.json')))
        output = json.loads(mols[0].write('pcjson'))
        self.assertEqual(output['PC_Compounds'][0]['stereo'], [
            {'planar': {'type': 1, 'ltop': 2, 'left': 3, 'right': 4, 'rbottom': 1, 'lbottom': 5, 'rtop': 6, 'parity': 255}}
        ])
        # Any
        mols = list(pybel.readfile("pcjson", os.path.join(filedir, 'CID_10900_2D.json')))
        output = json.loads(mols[0].write('pcjson'))
        self.assertEqual(output['PC_Compounds'][0]['stereo'][0]['planar']['left'], 3)
        self.assertEqual(output['PC_Compounds'][0]['stereo'][0]['planar']['right'], 4)
        self.assertEqual(output['PC_Compounds'][0]['stereo'][0]['planar']['left'], 3)

    def test_read(self):
        """Test reading a PubChem JSON file."""

        pcjson_input = """
{
  "PC_Compounds": [
    {
      "atoms": {
        "aid": [1,2,3,4,5,6,7,8,9,10,11,12,13],
        "element": [8,8,7,6,6,6,1,1,1,1,1,1,1]
      },
      "bonds": {
        "aid1": [1,1,2,4,3,3,4,4,4,5,5,5],
        "aid2": [6,13,6,3,11,12,5,6,7,8,9,10],
        "order": [1,1,2,1,1,1,1,1,1,1,1,1]
      },
      "charge": 0,
      "coords": [
        {
          "aid": [1,2,3,4,5,6,7,8,9,10,11,12,13],
          "conformers": [
            {
              "style": {"aid1": [4],"aid2": [7],"annotation": [6]},
              "x": [
                5.1350,
                4.2690,
                2.53690,
                3.4030,
                3.4030,
                4.2690,
                3.4030,
                2.7830,
                3.4030,
                4.0230,
                2.0,
                2.53690,
                5.6720
              ],
              "y": [
                -0.250,
                1.250,
                0.250,
                -0.250,
                -1.250,
                0.250,
                0.370,
                -1.250,
                -1.870,
                -1.250,
                -0.060,
                0.870,
                0.060
              ]
            }
          ]
        }
      ],
      "id": {"id": {"cid": 71080}},
      "stereo": [
        {
          "tetrahedral": {
            "above": 6,
            "below": 5,
            "bottom": 7,
            "center": 4,
            "parity": 1,
            "top": 3,
            "type": 1
          }
        }
      ]
    }
  ]
}
        """
        mol = pybel.readstring("pcjson", pcjson_input)
        self.assertEqual(mol.OBMol.NumAtoms(), 13)

    def test_old_string_format(self):
        """Test reading a PubChem JSON file using the old string format."""

        pcjson_input2 = """
{
  "PC_Compounds": [
    {
      "atoms": {"aid": [1,2,3,4,5,6],"element": ["cl","cl","c","c","h","h"]},
      "bonds": {"aid1": [1,2,3,3,4],"aid2": [4,3,4,5,6],"order": ["single","single","double","single","single"]},
      "charge": 0,
      "coords": [
        {
          "aid": [1,2,3,4,5,6],
          "conformers": [
            {
              "x": [4.59810,2.0,2.8660,3.7320,2.8660,3.7320],
              "y": [-0.250,0.250,-0.250,0.250,-0.870,0.870]
            }
          ]
        }
      ],
      "id": {"id": {"cid": 638186}},
      "stereo": [
        {
          "planar": {
            "left": 3,
            "ltop": 2,
            "lbottom": 5,
            "right": 4,
            "rtop": 6,
            "rbottom": 1,
            "parity": "opposite",
            "type": "planar"
          }
        }
      ]
    }
  ]
}
        """
        mol = pybel.readstring("pcjson", pcjson_input2)
        self.assertEqual(mol.OBMol.NumAtoms(), 6)


if __name__ == "__main__":
    unittest.main()

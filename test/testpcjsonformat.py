"""Example test that uses the OpenBabel Python bindings

On Windows or Linux, you can run these tests at the commandline
in the build folder with:
"C:\Program Files\CMake 2.6\bin\ctest.exe" -C CTestTestfile.cmake
                                           -R pybindtest -VV

You could also "chdir" into build/test and run the test file directly:
python ../../test/testpcjsonformat.py

In both cases, the test file is run directly from the source folder,
and so you can quickly develop the tests and try them out.
"""

import unittest
from testbindings import PybelWrapper, pybel


class TestPcJsonFormat(PybelWrapper):
    """Test PubChem JSON format."""

    def testRead(self):
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

    def testOldStringFormat(self):
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

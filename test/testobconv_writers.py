"""Test WriteFile() and WriteString() for each of the supported formats

These tests started as a way to verify that the OBConversion Index is
correctly reset to 1 for each format.

"""

import os
import sys
import unittest
import tempfile
from openbabel import openbabel as ob
import re

# Set the following to enable a workaround so the tests work on older
# versions of Open Babel.
ENABLE_WORKAROUND = 0


# Some of the formats embed the version in the output
VERSION = ob.OBReleaseVersion()

# Most of the tests use an OBMol made from this phenol structure

_default_conv = ob.OBConversion()
_default_conv.SetInAndOutFormats("sdf", "smi")

PHENOL_SDF = """\
phenol
 OpenBabel01151914482D

  7  7  0  0  0  0  0  0  0  0999 V2000
    1.5846   -0.0249    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5703    0.9755    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4295    1.4882    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3031    1.0004    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3175   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0005    0.0051    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  6  2  0  0  0  0
  1  2  1  0  0  0  0
  2  3  2  0  0  0  0
  3  4  1  0  0  0  0
  4  5  2  0  0  0  0
  5  6  1  0  0  0  0
  6  7  1  0  0  0  0
M  END
$$$$
"""

# Some of the tests use a reaction.

_rxn_conv = ob.OBConversion()
_rxn_conv.SetInAndOutFormats("rxn", "rxn")
_alchemy_mol = ob.OBMol()

ALCHEMY_RXN = """\
$RXN
lead_to_gold
      OpenBabel

  1  1
$MOL

 OpenBabel01151916222D

  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000   10.0000    0.0000 Pb  0  0  0  0  0 15  0  0  0  0  0  0
M  END
$MOL

 OpenBabel01151916222D

  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000   11.0000    0.0000 Au  0  0  0  0  0  0  0  0  0  0  0  0
M  END
"""
if not _rxn_conv.ReadString(_alchemy_mol, ALCHEMY_RXN):
    if ENABLE_WORKAROUND:
        # For some reason this record fails under Open Babel 2.4.1
        sys.stderr.write("Unable to parse RXN record? Reaction tests will fail.\n")
    else:
        raise AssertionError("Cannot parse RXN record")

# Some of the tests pass in a SMILES string
_smi_conv = ob.OBConversion()
_smi_conv.SetInAndOutFormats("smi", "smi")

def get_mol(test_case, mol):
    if mol is None:
        # Always make a new molecule so the tests don't
        # interfere with each other
        mol = ob.OBMol()
        assert _default_conv.ReadString(mol, PHENOL_SDF)
        return mol
    
    if isinstance(mol, str):
        # Parse it as a SMILES string with optional title
        terms = mol.split(None, 1)
        if len(terms) == 1:
            smiles = terms[0]
            title = "unknown"
        elif len(terms) == 2:
            smiles = terms[0]
            title = terms[1]
        else:
            raise AssertionError(mol)
            
        mol = ob.OBMol()
        if not _smi_conv.ReadString(mol, smiles):
            test_case.assert_("Cannot parse SMILES %r" % (smiles,))
        mol.SetTitle(title)
        return mol

    # Must have passed in a molecule. Return it.
    return mol

# Create a new OBConversion for the given format.
# Optionally pass in the options to set.
def get_converter(test_case, output_format, options=None):
    conv = ob.OBConversion()
    if not conv.SetInAndOutFormats("smi", output_format):
        test_case.assert_("Cannot set output format %r" % (output_format,))
    
    if options:
        # Can pass in a dictionary ...
        if isinstance(options, dict):
            for k, v in options.items():
                conv.AddOption(k, ob.OBConversion.OUTOPTIONS, v)
        else:
            # ... or an iterable
            for k in options:
                conv.AddOption(k, ob.OBConversion.OUTOPTIONS)

    if ENABLE_WORKAROUND:
        conv.SetOutputIndex(1)
    return conv

def save_to_pasteboard(text):
    # This test suite was developed on a Mac.
    # This code copies the text to the paste buffer,
    # which I can then use as the expected text.
    import subprocess
    p = subprocess.Popen(["pbcopy"],
                         stdin=subprocess.PIPE)
    p.stdin.write(text.encode("utf8"))
    p.stdin.close()
    p.wait()

def test_write_string(test_case, mol, conv, expected_output, normalize):
    output = conv.WriteString(mol)
    ### Debugging output
    if 0:
        print("===")
        print(output)
        print("===")
    if 0:
        save_to_pasteboard(text)

    # Apply normalizations to both sides
    if normalize is not None:
        output = normalize(output)
        expected_output = normalize(expected_output)
        
    test_case.assertMultiLineEqual(output.replace("\r\n", "\n"), expected_output.replace("\r\n", "\n"))

if type(u"") == type(""):
    # Python 3
    def test_binary_write_string(test_case, mol, conv, expected_output, normalize):
        # I think 'surrogateescape' is the right way to handle this
        output = conv.WriteString(mol).encode("utf8", "surrogateescape")
        if normalize:
            output = normalize(output)
            expected_output = normalize(expected_output)
        ## print("===", repr(output))
        test_case.assertEqual(output, expected_output)
else:
    # Python 2
    def test_binary_write_string(test_case, mol, conv, expected_output, normalize):
        output = conv.WriteString(mol)
##        print("===", repr(output))
        if normalize:
            output = normalize(output)
            expected_output = normalize(expected_output)
        test_case.assertEqual(output, expected_output)

def test_write_file(test_case, mol, conv, expected_output, normalize):
    temp_file_object = tempfile.NamedTemporaryFile(delete=False) # we will delete it manually
    temp_filename = temp_file_object.name
    if os.name == 'nt':
        temp_file_object.close() # Can't write to open file on Windows so we have to close it (but this could lead to a race condition if someone else uses the same temporary file name)
    try:
        test_case.assertTrue(conv.WriteFile(mol, temp_filename))
        conv.CloseOutFile() # we can't delete it on Windows otherwise
        with open(temp_filename) as f:
            output = f.read()
    finally:
        temp_file_object.close()
        os.remove(temp_filename)

    if 0:
        save_to_pasteboard(output)
        
    if normalize is not None:
        output = normalize(output)
        expected_output = normalize(expected_output)
    test_case.assertMultiLineEqual(output.replace("\r\n", "\n"), expected_output.replace("\r\n", "\n"))
    
def test_binary_write_file(test_case, mol, conv, expected_output, normalize):
    temp_file_object = tempfile.NamedTemporaryFile()
    temp_filename = temp_file_object.name
    try:
        test_case.assertTrue(conv.WriteFile(mol, temp_filename))
        with open(temp_filename, "rb") as f:
            output = f.read()
    finally:
        temp_file_object.close()

    ## print("==", repr(output))
    if normalize is not None:
        output = normalize(output)
        expected_output = normalize(expected_output)
    
    test_case.assertEqual(output, expected_output)

def test_write_multi_file(test_case, mols, conv, expected_output, normalize):
    temp_file_object = tempfile.NamedTemporaryFile()
    temp_filename = temp_file_object.name
    n = len(mols)
    test_case.assertGreater(n, 0, "must have at least one molecule")
    last = n-1
    
    try:
        for i, mol in enumerate(mols):
            conv.SetLast(i == last)
            if i == 0:
                test_case.assertTrue(conv.WriteFile(mol, temp_filename))
            else:
                test_case.assertTrue(conv.Write(mol))
        with open(temp_filename) as f:
            output = f.read()
    finally:
        temp_file_object.close()

    if 0:
        save_to_pasteboard(output)
        
    if normalize is not None:
        output = normalize(output)
        expected_output = normalize(expected_output)
    test_case.assertMultiLineEqual(output, expected_output)
    

class WriteMixin(object):
    def assertWriters(self, output_format, expected_output, options=None, mol=None, normalize=None):
        mol = get_mol(self, mol)
        conv = get_converter(self, output_format, options)
        test_write_string(self, mol, conv, expected_output, normalize)
        test_write_file(self, mol, conv, expected_output, normalize)
        
    def assertWriteString(self, output_format, expected_output, options=None, mol=None, normalize=None):
        mol = get_mol(self, mol)
        conv = get_converter(self, output_format, options)
        test_write_string(self, mol, conv, expected_output, normalize)
        
    def assertWriteFile(self, output_format, expected_output, options=None, mol=None, normalize=None):
        mol = get_mol(self, mol)
        conv = get_converter(self, output_format, options)
        test_write_file(self, mol, conv, expected_output, normalize)

    # Write 1 or more molecule to a file
    def assertWriteMultiFile(self, output_format, expected_output, options=None, mols=None, normalize=None):
        if mols is None:
            # Get two of the default molecules
            mols = [get_mol(self, None), get_mol(self, None)]
            
        conv = get_converter(self, output_format, options)
        test_write_multi_file(self, mols, conv, expected_output, normalize)
        
    def assertBinaryWriters(self, output_format, expected_output, options=None, mol=None, normalize=None):
        mol = get_mol(self, mol)
        conv = get_converter(self, output_format, options)
        test_binary_write_string(self, mol, conv, expected_output, normalize)
        test_binary_write_file(self, mol, conv, expected_output, normalize)
        
    def assertBinaryWriteString(self, output_format, expected_output, options=None, mol=None, normalize=None):
        mol = get_mol(self, mol)
        conv = get_converter(self, output_format, options)
        test_binary_write_string(self, mol, conv, expected_output, normalize)
        
    def assertBinaryWriteFile(self, output_format, expected_output, options=None, mol=None, normalize=None):
        mol = get_mol(self, mol)
        conv = get_converter(self, output_format, options)
        test_binary_write_file(self, mol, conv, expected_output, normalize)
        
    

# acesin -- ACES input format [Write-only]
class TestACES(unittest.TestCase, WriteMixin):
    fmt = "acesin"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
phenol
  C        1.58460       -0.02490        0.00000
  C        1.57030        0.97550        0.00000
  C        2.42950        1.48820        0.00000
  C        3.30310        1.00040        0.00000
  C        3.31750       -0.00000        0.00000
  C        0.00000        0.00000        0.00000
  O       -1.00050        0.00510        0.00000

*ACES2(__ADD_SETUP_HERE__)

""")

# adf -- ADF cartesian input format [Write-only]
class TestADF(unittest.TestCase, WriteMixin):
    fmt = "adf"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
TITLE phenol

CHARGE 0  0

Number of atoms
 7

ATOMS Cartesian
C          1.58460       -0.02490        0.00000
C          1.57030        0.97550        0.00000
C          2.42950        1.48820        0.00000
C          3.30310        1.00040        0.00000
C          3.31750       -0.00000        0.00000
C          0.00000        0.00000        0.00000
O         -1.00050        0.00510        0.00000
End

Basis
End

Geometry
End


""")

# alc -- Alchemy format
class TestALC(unittest.TestCase, WriteMixin):
    fmt = "alc"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
    7 ATOMS,     7 BONDS,     0 CHARGES
    1 C2      1.5846  -0.0249   0.0000     0.0000
    2 C2      1.5703   0.9755   0.0000     0.0000
    3 C2      2.4295   1.4882   0.0000     0.0000
    4 C2      3.3031   1.0004   0.0000     0.0000
    5 C2      3.3175  -0.0000   0.0000     0.0000
    6 C2      0.0000   0.0000   0.0000     0.0000
    7 O3     -1.0005   0.0051   0.0000     0.0000
    1     1     6  DOUBLE
    2     1     2  SINGLE
    3     2     3  DOUBLE
    4     3     4  SINGLE
    5     4     5  DOUBLE
    6     5     6  SINGLE
    7     6     7  SINGLE
""")

## # ascii -- ASCII format [Write-only]
## # XXX Doesn't look good
## class TestASCII(unittest.TestCase, WriteMixin):
##     fmt = "ascii"
##     maxDiff = None
##     def test_default(self):
##         self.assertWriters(self.fmt, """\
## """)

# bgf -- MSI BGF format
class TestBGF(unittest.TestCase, WriteMixin):
    fmt = "bgf"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
BIOGRF 200
DESCRP phenol
FORCEFIELD DREIDING  
FORMAT ATOM   (a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,3f10.5,1x,a5,i3,i2,1x,f8.5)
HETATM     1 C1    RES A   444   1.58460  -0.02490   0.00000 C_R    3 0  0.04203
HETATM     2 C2    RES A   444   1.57030   0.97550   0.00000 C_R    3 0  0.00328
HETATM     3 C3    RES A   444   2.42950   1.48820   0.00000 C_R    3 0  0.00021
HETATM     4 C4    RES A   444   3.30310   1.00040   0.00000 C_R    3 0  0.00328
HETATM     5 C5    RES A   444   3.31750  -0.00000   0.00000 C_R    3 0  0.04203
HETATM     6 C6    RES A   444   0.00000   0.00000   0.00000 C_R    3 0  0.19575
HETATM     7 O7    RES A   444  -1.00050   0.00510   0.00000 O_3    2 0 -0.28657
FORMAT CONECT (a6,12i6)

CONECT     1     6     2
ORDER      1     2     1
CONECT     2     1     3
ORDER      2     1     2
CONECT     3     2     4
ORDER      3     2     1
CONECT     4     3     5
ORDER      4     1     2
CONECT     5     4     6
ORDER      5     2     1
CONECT     6     1     5     7
ORDER      6     2     1     1
CONECT     7     6
ORDER      7     1
END
""")

# box -- Dock 3.5 Box format
class TestBOX(unittest.TestCase, WriteMixin):
    fmt = "box"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
HEADER    CORNERS OF BOX
REMARK    CENTER (X Y Z)           1.159      0.732      0.000
REMARK    DIMENSIONS (X Y Z)       6.318      3.513      2.000
ATOM      1  DUA BOX     1      -2.000  -1.025  -1.000
ATOM      2  DUA BOX     1       4.317  -1.025  -1.000
ATOM      3  DUA BOX     1       4.317  -1.025   1.000
ATOM      4  DUA BOX     1      -2.000  -1.025   1.000
ATOM      5  DUA BOX     1      -2.000   2.488  -1.000
ATOM      6  DUA BOX     1       4.317   2.488  -1.000
ATOM      7  DUA BOX     1       4.317   2.488   1.000
ATOM      8  DUA BOX     1      -2.000   2.488   1.000
CONECT    1    2    4    5
CONECT    2    1    3    6
CONECT    3    2    4    7
CONECT    4    1    3    8
CONECT    5    1    6    8
CONECT    6    2    5    7
CONECT    7    3    6    8
CONECT    8    4    5    7
""")

# bs -- Ball and Stick format
class TestBS(unittest.TestCase, WriteMixin):
    fmt = "bs"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
phenol
7
C     1.5846   -0.0249    0.0000     6     2
C     1.5703    0.9755    0.0000     1     3
C     2.4295    1.4882    0.0000     2     4
C     3.3031    1.0004    0.0000     3     5
C     3.3175   -0.0000    0.0000     4     6
C     0.0000    0.0000    0.0000     1     5     7
O    -1.0005    0.0051    0.0000     6
""")

# c3d1 -- Chem3D Cartesian 1 format
class TestC3D1(unittest.TestCase, WriteMixin):
    fmt = "c3d1"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
7
C   1       1.5846   -0.0249    0.0000     2     6     2
C   2       1.5703    0.9755    0.0000     2     1     3
C   3       2.4295    1.4882    0.0000     2     2     4
C   4       3.3031    1.0004    0.0000     2     3     5
C   5       3.3175   -0.0000    0.0000     2     4     6
C   6       0.0000    0.0000    0.0000     2     1     5     7
O   7      -1.0005    0.0051    0.0000     6     6
""")

# c3d2 -- Chem3D Cartesian 2 format
class TestC3D2(unittest.TestCase, WriteMixin):
    fmt = "c3d2"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
7
C   1       1.5846   -0.0249    0.0000     2     6     2
C   2       1.5703    0.9755    0.0000     2     1     3
C   3       2.4295    1.4882    0.0000     2     2     4
C   4       3.3031    1.0004    0.0000     2     3     5
C   5       3.3175   -0.0000    0.0000     2     4     6
C   6       0.0000    0.0000    0.0000     2     1     5     7
O   7      -1.0005    0.0051    0.0000    82     6
""")

# cac -- CAChe MolStruct format [Write-only]
class TestCAC(unittest.TestCase, WriteMixin):
    fmt = "cac"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
molstruct88_Apr_30_1993_11:02:29 <molecule> 0x1d00
Written by Molecular Editor on <date>
Using data dictionary         9/9/93  4:47 AM
Version 6
local_transform
0.100000 0.000000 0.000000 0.000000
0.000000 0.100000 0.000000 0.000000
0.000000 0.000000 0.100000 0.000000
0.000000 0.000000 0.000000 1.000000
object_class atom
property xyz_coordinates MoleculeEditor angstrom 6 3 FLOAT
property anum MoleculeEditor unit 0 1 INTEGER
property sym MoleculeEditor noUnit 0 2 STRING
property chrg MoleculeEditor charge_au 0 1 INTEGER
property rflag MoleculeEditor noUnit 0 1 HEX
ID xyz_coordinates             anum sym	chrg rflag
  1   1.584600  -0.024900   0.000000  6  C  0 0x7052
  2   1.570300   0.975500   0.000000  6  C  0 0x7052
  3   2.429500   1.488200   0.000000  6  C  0 0x7052
  4   3.303100   1.000400   0.000000  6  C  0 0x7052
  5   3.317500  -0.000000   0.000000  6  C  0 0x7052
  6   0.000000   0.000000   0.000000  6  C  0 0x7052
  7  -1.000500   0.005100   0.000000  8  O  0 0x7052
property_flags:
object_class bond
property rflag MoleculeEditor noUnit 0 1 HEX
property type MoleculeEditor noUnit 0 1 NAME
property bond_order MoleculeEditor noUnit 4 1 FLOAT
ID rflag type bond_order
  1 0x7005 double
  2 0x7005 single
  3 0x7005 double
  4 0x7005 single
  5 0x7005 double
  6 0x7005 single
  7 0x7005 single
property_flags:
object_class connector
property dflag MoleculeEditor noUnit 0 1 HEX
property objCls1 MoleculeEditor noUnit 0 1 NAME
property objCls2 MoleculeEditor noUnit 0 1 NAME
property objID1 MoleculeEditor noUnit 0 1 INTEGER
property objID2 MoleculeEditor noUnit 0 1 INTEGER
ID dflag objCls1 objCls2 objID1 objID2
  1 0xa1 atom bond 1 1
  2 0xa1 atom bond 6 1
  3 0xa1 atom bond 1 2
  4 0xa1 atom bond 2 2
  5 0xa1 atom bond 2 3
  6 0xa1 atom bond 3 3
  7 0xa1 atom bond 3 4
  8 0xa1 atom bond 4 4
  9 0xa1 atom bond 4 5
 10 0xa1 atom bond 5 5
 11 0xa1 atom bond 5 6
 12 0xa1 atom bond 6 6
 13 0xa1 atom bond 6 7
 14 0xa1 atom bond 7 7
property_flags:
""")

# caccrt -- Cacao Cartesian format
class TestCACCRT(unittest.TestCase, WriteMixin):
    fmt = "caccrt"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
phenol
  7   DIST  0  0  0
CELL 1.,1.,1.,90.,90.,90.
 C  1.5846, -0.0249,  0.0000
 C  1.5703,  0.9755,  0.0000
 C  2.4295,  1.4882,  0.0000
 C  3.3031,  1.0004,  0.0000
 C  3.3175, -0.0000,  0.0000
 C  0.0000,  0.0000,  0.0000
 O -1.0005,  0.0051,  0.0000
""")

# cache -- CAChe MolStruct format [Write-only]
class TestCACHE(unittest.TestCase, WriteMixin):
    fmt = "cache"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
molstruct88_Apr_30_1993_11:02:29 <molecule> 0x1d00
Written by Molecular Editor on <date>
Using data dictionary         9/9/93  4:47 AM
Version 6
local_transform
0.100000 0.000000 0.000000 0.000000
0.000000 0.100000 0.000000 0.000000
0.000000 0.000000 0.100000 0.000000
0.000000 0.000000 0.000000 1.000000
object_class atom
property xyz_coordinates MoleculeEditor angstrom 6 3 FLOAT
property anum MoleculeEditor unit 0 1 INTEGER
property sym MoleculeEditor noUnit 0 2 STRING
property chrg MoleculeEditor charge_au 0 1 INTEGER
property rflag MoleculeEditor noUnit 0 1 HEX
ID xyz_coordinates             anum sym	chrg rflag
  1   1.584600  -0.024900   0.000000  6  C  0 0x7052
  2   1.570300   0.975500   0.000000  6  C  0 0x7052
  3   2.429500   1.488200   0.000000  6  C  0 0x7052
  4   3.303100   1.000400   0.000000  6  C  0 0x7052
  5   3.317500  -0.000000   0.000000  6  C  0 0x7052
  6   0.000000   0.000000   0.000000  6  C  0 0x7052
  7  -1.000500   0.005100   0.000000  8  O  0 0x7052
property_flags:
object_class bond
property rflag MoleculeEditor noUnit 0 1 HEX
property type MoleculeEditor noUnit 0 1 NAME
property bond_order MoleculeEditor noUnit 4 1 FLOAT
ID rflag type bond_order
  1 0x7005 double
  2 0x7005 single
  3 0x7005 double
  4 0x7005 single
  5 0x7005 double
  6 0x7005 single
  7 0x7005 single
property_flags:
object_class connector
property dflag MoleculeEditor noUnit 0 1 HEX
property objCls1 MoleculeEditor noUnit 0 1 NAME
property objCls2 MoleculeEditor noUnit 0 1 NAME
property objID1 MoleculeEditor noUnit 0 1 INTEGER
property objID2 MoleculeEditor noUnit 0 1 INTEGER
ID dflag objCls1 objCls2 objID1 objID2
  1 0xa1 atom bond 1 1
  2 0xa1 atom bond 6 1
  3 0xa1 atom bond 1 2
  4 0xa1 atom bond 2 2
  5 0xa1 atom bond 2 3
  6 0xa1 atom bond 3 3
  7 0xa1 atom bond 3 4
  8 0xa1 atom bond 4 4
  9 0xa1 atom bond 4 5
 10 0xa1 atom bond 5 5
 11 0xa1 atom bond 5 6
 12 0xa1 atom bond 6 6
 13 0xa1 atom bond 6 7
 14 0xa1 atom bond 7 7
property_flags:
""")

# cacint -- Cacao Internal format [Write-only]
class TestCACINT(unittest.TestCase, WriteMixin):
    fmt = "cacint"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
 # TITLE
  EL
0.,0.,0., C
 1,2, C  1.001, 90.000,269.181
 2,3, C  1.001,120.006, 90.000
 3,4, C  1.001,119.997, 49.094
 4,5, C  1.001,120.003, -0.000
 1,6, C  1.585, 90.000, -0.000
 6,7, O  1.001,179.392,180.000
""")

# can -- Canonical SMILES format
class TestCAN(unittest.TestCase, WriteMixin):
    fmt = "can"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
Oc1ccccc1\tphenol
""")

# cdjson -- ChemDoodle JSON
class TestCDJSON(unittest.TestCase, WriteMixin):
    fmt = "cdjson"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
{
  "m": [
    {
      "a": [
        {
          "x": 31.692,
          "y": -0.498
        },
        {
          "x": 31.406,
          "y": 19.51
        },
        {
          "x": 48.59,
          "y": 29.764
        },
        {
          "x": 66.062,
          "y": 20.008
        },
        {
          "x": 66.35,
          "y": -0.0
        },
        {
          "x": 0.0,
          "y": 0.0
        },
        {
          "x": -20.009999999999999,
          "y": 0.10200000000000001,
          "l": 8
        }
      ],
      "b": [
        {
          "b": 0,
          "e": 5,
          "o": 2
        },
        {
          "b": 0,
          "e": 1
        },
        {
          "b": 1,
          "e": 2,
          "o": 2
        },
        {
          "b": 2,
          "e": 3
        },
        {
          "b": 3,
          "e": 4,
          "o": 2
        },
        {
          "b": 4,
          "e": 5
        },
        {
          "b": 5,
          "e": 6
        }
      ]
    }
  ]
}""")

## # cdxml -- ChemDraw CDXML format
## XXX fails on an unpatched system
## class TestCDXML(unittest.TestCase, WriteMixin):
##     fmt = "cdxml"
##     maxDiff = None
##     def test_default(self):
##         self.assertWriters(self.fmt, """\
## <fragment>
##  <n id="1" p="0.000000 0.000000"/>
##  <n id="2" p="0.000000 0.000000"/>
##  <n id="3" p="0.000000 0.000000"/>
##  <n id="4" p="0.000000 0.000000"/>
##  <n id="5" p="0.000000 0.000000"/>
##  <n id="6" p="0.000000 0.000000"/>
##  <n id="7" p="0.000000 0.000000" Element="8"/>
##  <b B="1" E="2"/>
##  <b B="2" E="3" Order="2"/>
##  <b B="3" E="4"/>
##  <b B="4" E="5" Order="2"/>
##  <b B="5" E="6"/>
##  <b B="1" E="6" Order="2"/>
##  <b B="6" E="7"/>
## </fragment>
## """)

# cht -- Chemtool format [Write-only]
class TestCHT(unittest.TestCase, WriteMixin):
    fmt = "cht"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
Chemtool Version 1.4
geometry 165 55
bonds 7
79	-1	0	0	1
79	-1	79	49	0
79	49	121	74	1
121	74	165	50	0
165	50	166	0	1
166	0	0	0	0
0	0	-50	0	0
atoms 1
-50	0	O	-1
splines 0
""")

# cif -- Crystallographic Information File
class TestCIF(unittest.TestCase, WriteMixin):
    fmt = "cif"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
# CIF file generated by openbabel %(VERSION)s, see https://openbabel.org
data_I
_chemical_name_common 'phenol'
loop_
    _atom_site_label
    _atom_site_type_symbol
    _atom_site_fract_x
    _atom_site_fract_y
    _atom_site_fract_z
    _atom_site_occupancy
    C0       C        1.58460   -0.02490    0.00000    1.000
    C1       C        1.57030    0.97550    0.00000    1.000
    C2       C        2.42950    1.48820    0.00000    1.000
    C3       C        3.30310    1.00040    0.00000    1.000
    C4       C        3.31750   -0.00000    0.00000    1.000
    C5       C        0.00000    0.00000    0.00000    1.000
    O6       O       -1.00050    0.00510    0.00000    1.000
""" % dict(VERSION=VERSION))

## # ck -- ChemKin format
## XXX I don't know why this fails
## class TestCK(unittest.TestCase, WriteMixin):
##     fmt = "ck"
##     maxDiff = None
##     def test_default(self):
##         self.assertWriters(self.fmt, """\
## """)

# cml -- Chemical Markup Language
class TestCML(unittest.TestCase, WriteMixin):
    fmt = "cml"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
<?xml version="1.0"?>
<molecule id="phenol" xmlns="http://www.xml-cml.org/schema">
 <atomArray>
  <atom id="a1" elementType="C" hydrogenCount="1" x2="1.584600" y2="-0.024900"/>
  <atom id="a2" elementType="C" hydrogenCount="1" x2="1.570300" y2="0.975500"/>
  <atom id="a3" elementType="C" hydrogenCount="1" x2="2.429500" y2="1.488200"/>
  <atom id="a4" elementType="C" hydrogenCount="1" x2="3.303100" y2="1.000400"/>
  <atom id="a5" elementType="C" hydrogenCount="1" x2="3.317500" y2="-0.000000"/>
  <atom id="a6" elementType="C" hydrogenCount="0" x2="0.000000" y2="0.000000"/>
  <atom id="a7" elementType="O" hydrogenCount="1" x2="-1.000500" y2="0.005100"/>
 </atomArray>
 <bondArray>
  <bond atomRefs2="a1 a6" order="2"/>
  <bond atomRefs2="a1 a2" order="1"/>
  <bond atomRefs2="a2 a3" order="2"/>
  <bond atomRefs2="a3 a4" order="1"/>
  <bond atomRefs2="a4 a5" order="2"/>
  <bond atomRefs2="a5 a6" order="1"/>
  <bond atomRefs2="a6 a7" order="1"/>
 </bondArray>
</molecule>
""")

    def test_multimol_default(self):
        # Write two phenols.
        # When there are 2 or more molecules then each molecule
        # is wrapped in a <cml> element.
        self.assertWriteMultiFile("cml", """\
<?xml version="1.0"?>
<cml xmlns="http://www.xml-cml.org/schema">
 <molecule id="phenol">
  <atomArray>
   <atom id="a1" elementType="C" hydrogenCount="1" x2="1.584600" y2="-0.024900"/>
   <atom id="a2" elementType="C" hydrogenCount="1" x2="1.570300" y2="0.975500"/>
   <atom id="a3" elementType="C" hydrogenCount="1" x2="2.429500" y2="1.488200"/>
   <atom id="a4" elementType="C" hydrogenCount="1" x2="3.303100" y2="1.000400"/>
   <atom id="a5" elementType="C" hydrogenCount="1" x2="3.317500" y2="-0.000000"/>
   <atom id="a6" elementType="C" hydrogenCount="0" x2="0.000000" y2="0.000000"/>
   <atom id="a7" elementType="O" hydrogenCount="1" x2="-1.000500" y2="0.005100"/>
  </atomArray>
  <bondArray>
   <bond atomRefs2="a1 a6" order="2"/>
   <bond atomRefs2="a1 a2" order="1"/>
   <bond atomRefs2="a2 a3" order="2"/>
   <bond atomRefs2="a3 a4" order="1"/>
   <bond atomRefs2="a4 a5" order="2"/>
   <bond atomRefs2="a5 a6" order="1"/>
   <bond atomRefs2="a6 a7" order="1"/>
  </bondArray>
 </molecule>
 <molecule id="phenol">
  <atomArray>
   <atom id="a1" elementType="C" hydrogenCount="1" x2="1.584600" y2="-0.024900"/>
   <atom id="a2" elementType="C" hydrogenCount="1" x2="1.570300" y2="0.975500"/>
   <atom id="a3" elementType="C" hydrogenCount="1" x2="2.429500" y2="1.488200"/>
   <atom id="a4" elementType="C" hydrogenCount="1" x2="3.303100" y2="1.000400"/>
   <atom id="a5" elementType="C" hydrogenCount="1" x2="3.317500" y2="-0.000000"/>
   <atom id="a6" elementType="C" hydrogenCount="0" x2="0.000000" y2="0.000000"/>
   <atom id="a7" elementType="O" hydrogenCount="1" x2="-1.000500" y2="0.005100"/>
  </atomArray>
  <bondArray>
   <bond atomRefs2="a1 a6" order="2"/>
   <bond atomRefs2="a1 a2" order="1"/>
   <bond atomRefs2="a2 a3" order="2"/>
   <bond atomRefs2="a3 a4" order="1"/>
   <bond atomRefs2="a4 a5" order="2"/>
   <bond atomRefs2="a5 a6" order="1"/>
   <bond atomRefs2="a6 a7" order="1"/>
  </bondArray>
 </molecule>
</cml>
""")
        


## # cmlr -- CML Reaction format
## XXX I don't know why the result is the empty string
## class TestCMLR(unittest.TestCase, WriteMixin):
##     fmt = "cmlr"
##     maxDiff = None
##     def test_default(self):
##         self.assertWriters(self.fmt, """\
## """, mol=_alchemy_mol)

# com -- Gaussian 98/03 Input [Write-only]
class TestCOM(unittest.TestCase, WriteMixin):
    fmt = "com"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
!Put Keywords Here, check Charge and Multiplicity.
#

 phenol

0  1
C           1.58460        -0.02490         0.00000
C           1.57030         0.97550         0.00000
C           2.42950         1.48820         0.00000
C           3.30310         1.00040         0.00000
C           3.31750        -0.00000         0.00000
C           0.00000         0.00000         0.00000
O          -1.00050         0.00510         0.00000

""")

## # confabreport -- Confab report format [Write-only]
## XXX no conformations
## class TestCONFABREPORT(unittest.TestCase, WriteMixin):
##     fmt = "confabreport"
##     maxDiff = None
##     def test_default(self):
##         self.assertWriters(self.fmt, """\
## """)

# CONFIG -- DL-POLY CONFIG
class TestCONFIG(unittest.TestCase, WriteMixin):
    fmt = "CONFIG"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
phenol
         0         0
       C         1         6
   1.584600000000000   -0.024900000000000    0.000000000000000
       C         2         6
   1.570300000000000    0.975500000000000    0.000000000000000
       C         3         6
   2.429500000000000    1.488200000000000    0.000000000000000
       C         4         6
   3.303100000000000    1.000400000000000    0.000000000000000
       C         5         6
   3.317500000000000   -0.000000000000000    0.000000000000000
       C         6         6
   0.000000000000000    0.000000000000000    0.000000000000000
       O         7         8
  -1.000500000000000    0.005100000000000    0.000000000000000
""")

# CONTCAR -- VASP format
class TestCONTCAR(unittest.TestCase, WriteMixin):
    fmt = "CONTCAR"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
phenol
1.000 
0.0  0.0  0.0
0.0  0.0  0.0
0.0  0.0  0.0
C   O   
6   1   
Cartesian
     1.5846000000000000085     -0.0248999999999999985      0.0000000000000000000
     1.5703000000000000291      0.9755000000000000338      0.0000000000000000000
     2.4294999999999999929      1.4881999999999999673      0.0000000000000000000
     3.3031000000000001471      1.0003999999999999559      0.0000000000000000000
     3.3174999999999998934     -0.0000000000000000000      0.0000000000000000000
     0.0000000000000000000      0.0000000000000000000      0.0000000000000000000
    -1.0004999999999999449      0.0051000000000000004      0.0000000000000000000
""")

# CONTFF -- MDFF format
class TestCONTFF(unittest.TestCase, WriteMixin):
    fmt = "CONTFF"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
7
phenol
0.0  0.0  0.0
0.0  0.0  0.0
0.0  0.0  0.0
2
C   O   
6   1   
Cartesian
C        1.5846000000000000085     -0.0248999999999999985      0.0000000000000000000
C        1.5703000000000000291      0.9755000000000000338      0.0000000000000000000
C        2.4294999999999999929      1.4881999999999999673      0.0000000000000000000
C        3.3031000000000001471      1.0003999999999999559      0.0000000000000000000
C        3.3174999999999998934     -0.0000000000000000000      0.0000000000000000000
C        0.0000000000000000000      0.0000000000000000000      0.0000000000000000000
O       -1.0004999999999999449      0.0051000000000000004      0.0000000000000000000
""")

## # copy -- Copy raw text [Write-only]
## XXX "Not a valid output format"
## class TestCOPY(unittest.TestCase, WriteMixin):
##     fmt = "copy"
##     maxDiff = None
##     def test_default(self):
##         self.assertWriters(self.fmt, """\
## """)

# crk2d -- Chemical Resource Kit diagram(2D)
class TestCRK2D(unittest.TestCase, WriteMixin):
    fmt = "crk2d"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
<Property Type="DiagramStructure">
 <Structure2D>
  <Group Charge="0" Spin="0">
   <Atom ID="1">
    <X>1.5846</X>
    <Y>-0.0249</Y>
    <Z>0</Z>
    <Element>C</Element>
   </Atom>
   <Atom ID="2">
    <X>1.5703</X>
    <Y>0.9755</Y>
    <Z>0</Z>
    <Element>C</Element>
   </Atom>
   <Atom ID="3">
    <X>2.4295</X>
    <Y>1.4882</Y>
    <Z>0</Z>
    <Element>C</Element>
   </Atom>
   <Atom ID="4">
    <X>3.3031</X>
    <Y>1.0004</Y>
    <Z>0</Z>
    <Element>C</Element>
   </Atom>
   <Atom ID="5">
    <X>3.3175</X>
    <Y>-0</Y>
    <Z>0</Z>
    <Element>C</Element>
   </Atom>
   <Atom ID="6">
    <X>0</X>
    <Y>0</Y>
    <Z>0</Z>
    <Element>C</Element>
   </Atom>
   <Atom ID="7">
    <X>-1.0005</X>
    <Y>0.0051</Y>
    <Z>0</Z>
    <Element>O</Element>
   </Atom>
   <Bond>
    <From>1</From>
    <To>6</To>
    <Order>1.5</Order>
    <Style>0</Style>
   </Bond>
   <Bond>
    <From>1</From>
    <To>2</To>
    <Order>1.5</Order>
    <Style>0</Style>
   </Bond>
   <Bond>
    <From>2</From>
    <To>3</To>
    <Order>1.5</Order>
    <Style>0</Style>
   </Bond>
   <Bond>
    <From>3</From>
    <To>4</To>
    <Order>1.5</Order>
    <Style>0</Style>
   </Bond>
   <Bond>
    <From>4</From>
    <To>5</To>
    <Order>1.5</Order>
    <Style>0</Style>
   </Bond>
   <Bond>
    <From>5</From>
    <To>6</To>
    <Order>1.5</Order>
    <Style>0</Style>
   </Bond>
   <Bond>
    <From>6</From>
    <To>7</To>
    <Order>1</Order>
    <Style>0</Style>
   </Bond>
  </Group>
 </Structure2D>
</Property>
""")

# crk3d -- Chemical Resource Kit 3D format
class TestCRK3D(unittest.TestCase, WriteMixin):
    fmt = "crk3d"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
<Property Type="ModelStructure">
 <Structure3D>
  <Group Charge="0" Spin="0">
   <Atom ID="1">
    <X>1.5846</X>
    <Y>-0.0249</Y>
    <Z>0</Z>
    <Element>C</Element>
   </Atom>
   <Atom ID="2">
    <X>1.5703</X>
    <Y>0.9755</Y>
    <Z>0</Z>
    <Element>C</Element>
   </Atom>
   <Atom ID="3">
    <X>2.4295</X>
    <Y>1.4882</Y>
    <Z>0</Z>
    <Element>C</Element>
   </Atom>
   <Atom ID="4">
    <X>3.3031</X>
    <Y>1.0004</Y>
    <Z>0</Z>
    <Element>C</Element>
   </Atom>
   <Atom ID="5">
    <X>3.3175</X>
    <Y>-0</Y>
    <Z>0</Z>
    <Element>C</Element>
   </Atom>
   <Atom ID="6">
    <X>0</X>
    <Y>0</Y>
    <Z>0</Z>
    <Element>C</Element>
   </Atom>
   <Atom ID="7">
    <X>-1.0005</X>
    <Y>0.0051</Y>
    <Z>0</Z>
    <Element>O</Element>
   </Atom>
   <Bond>
    <From>1</From>
    <To>6</To>
    <Order>1.5</Order>
    <Style>0</Style>
   </Bond>
   <Bond>
    <From>1</From>
    <To>2</To>
    <Order>1.5</Order>
    <Style>0</Style>
   </Bond>
   <Bond>
    <From>2</From>
    <To>3</To>
    <Order>1.5</Order>
    <Style>0</Style>
   </Bond>
   <Bond>
    <From>3</From>
    <To>4</To>
    <Order>1.5</Order>
    <Style>0</Style>
   </Bond>
   <Bond>
    <From>4</From>
    <To>5</To>
    <Order>1.5</Order>
    <Style>0</Style>
   </Bond>
   <Bond>
    <From>5</From>
    <To>6</To>
    <Order>1.5</Order>
    <Style>0</Style>
   </Bond>
   <Bond>
    <From>6</From>
    <To>7</To>
    <Order>1</Order>
    <Style>0</Style>
   </Bond>
  </Group>
 </Structure3D>
</Property>
""")

# csr -- Accelrys/MSI Quanta CSR format [Write-only]
class TestCSR(unittest.TestCase, WriteMixin):
    fmt = "csr"
    maxDiff = None
    def test_default(self):
        self.assertBinaryWriters(self.fmt, b'\x04\x00\x00\x00V33 \x04\x00\x00\x00\x08\x00\x00\x00\x07\x00\x00\x00\x01\x00\x00\x00\x08\x00\x00\x00d\x00\x00\x00phenol                                                                                             \x00d\x00\x00\x00\x04\x00\x00\x00\x07\x00\x00\x00\x04\x00\x00\x00\\\x00\x00\x00\x01\x00\x00\x00\x05\x17+j0\xad\x04\xc0phenol:1                                                                       \x00\\\x00\x00\x008\x00\x00\x00\x98\xdd\x93\x87\x85Z\xf9?r\x8a\x8e\xe4\xf2\x1f\xf9?V\x0e-\xb2\x9do\x03@?W[\xb1\xbfl\n@\n\xd7\xa3p=\x8a\n@\x00\x00\x00\x00\x00\x00\x00\x005^\xbaI\x0c\x02\xf0\xbf8\x00\x00\x008\x00\x00\x00V}\xae\xb6b\x7f\x99\xbf\x9e\xef\xa7\xc6K7\xef?\xe4\x83\x9e\xcd\xaa\xcf\xf7?\xc4\xb1.n\xa3\x01\xf0?\x00\x00\x00\x00\x00\x00\x00\x80\x00\x00\x00\x00\x00\x00\x00\x00\x88\x85Z\xd3\xbc\xe3t?8\x00\x00\x008\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x008\x00\x00\x00')

# cssr -- CSD CSSR format [Write-only]
class TestCSSR(unittest.TestCase, WriteMixin):
    fmt = "cssr"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
 REFERENCE STRUCTURE = 00000   A,B,C =   1.000   1.000   1.000
   ALPHA,BETA,GAMMA =  90.000  90.000  90.000    SPGR =    P1
   7   1 phenol


   1 C1      1.58460  -0.02490   0.00000    6   2   0   0   0   0   0   0   0.042   1
   2 C2      1.57030   0.97550   0.00000    1   3   0   0   0   0   0   0   0.003   1
   3 C3      2.42950   1.48820   0.00000    2   4   0   0   0   0   0   0   0.000   1
   4 C4      3.30310   1.00040   0.00000    3   5   0   0   0   0   0   0   0.003   1
   5 C5      3.31750  -0.00000   0.00000    4   6   0   0   0   0   0   0   0.042   1
   6 C6      0.00000   0.00000   0.00000    1   5   7   0   0   0   0   0   0.196   1
   7 O1     -1.00050   0.00510   0.00000    6   0   0   0   0   0   0   0  -0.287   1
""")

# ct -- ChemDraw Connection Table format
class TestCT(unittest.TestCase, WriteMixin):
    fmt = "ct"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
phenol
 7 7
    1.5846   -0.0249    0.0000 C
    1.5703    0.9755    0.0000 C
    2.4295    1.4882    0.0000 C
    3.3031    1.0004    0.0000 C
    3.3175   -0.0000    0.0000 C
    0.0000    0.0000    0.0000 C
   -1.0005    0.0051    0.0000 O
  1  6  2  2
  1  2  1  1
  2  3  2  2
  3  4  1  1
  4  5  2  2
  5  6  1  1
  6  7  1  1
""")

## # cub -- Gaussian cube format
## XXX "The molecule has no grid."
## class TestCUB(unittest.TestCase, WriteMixin):
##     fmt = "cub"
##     maxDiff = None
##     def test_default(self):
##         self.assertWriters(self.fmt, """\
## """)

## # cube -- Gaussian cube format
## XXX "The molecule has no grid."
## class TestCUBE(unittest.TestCase, WriteMixin):
##     fmt = "cube"
##     maxDiff = None
##     def test_default(self):
##         self.assertWriters(self.fmt, """\
## """)

# dalmol -- DALTON input format
class TestDALMOL(unittest.TestCase, WriteMixin):
    fmt = "dalmol"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
BASIS
6-31G*
phenol
Generated by Open Babel. Check overall charge below.
AtomTypes=2 Charge=0 NoSymmetry Angstrom
Charge=6.0 Atoms=6
C             1.5846000000   -0.0249000000    0.0000000000 
C             1.5703000000    0.9755000000    0.0000000000 
C             2.4295000000    1.4882000000    0.0000000000 
C             3.3031000000    1.0004000000    0.0000000000 
C             3.3175000000   -0.0000000000    0.0000000000 
C             0.0000000000    0.0000000000    0.0000000000 
Charge=8.0 Atoms=1
O            -1.0005000000    0.0051000000    0.0000000000 
""")

# dmol -- DMol3 coordinates format
class TestDMOL(unittest.TestCase, WriteMixin):
    fmt = "dmol"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
$coordinates
C             2.99445980216940   -0.04705417712610    0.00000000000000
C             2.96743672052670    1.84342770226950    0.00000000000000
C             4.59108929027550    2.81229021682980    0.00000000000000
C             6.24195391426590    1.89048187939560    0.00000000000000
C             6.26916596850750   -0.00000000000000    0.00000000000000
C             0.00000000000000    0.00000000000000    0.00000000000000
O            -1.89067085199450    0.00963760254390    0.00000000000000
$end
""")

## # dx -- OpenDX cube format for APBS
## XXX "The molecule has no grid."
## class TestDX(unittest.TestCase, WriteMixin):
##     fmt = "dx"
##     maxDiff = None
##     def test_default(self):
##         self.assertWriters(self.fmt, """\
## """)

# ent -- Protein Data Bank format
class TestENT(unittest.TestCase, WriteMixin):
    fmt = "ent"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
COMPND    phenol 
AUTHOR    GENERATED BY OPEN BABEL %(VERSION)s
HETATM    1  C   UNL     1       1.585  -0.025   0.000  1.00  0.00           C  
HETATM    2  C   UNL     1       1.570   0.976   0.000  1.00  0.00           C  
HETATM    3  C   UNL     1       2.429   1.488   0.000  1.00  0.00           C  
HETATM    4  C   UNL     1       3.303   1.000   0.000  1.00  0.00           C  
HETATM    5  C   UNL     1       3.317   0.000   0.000  1.00  0.00           C  
HETATM    6  C   UNL     1       0.000   0.000   0.000  1.00  0.00           C  
HETATM    7  O   UNL     1      -1.000   0.005   0.000  1.00  0.00           O  
CONECT    1    6    6    2                                            
CONECT    2    1    3    3                                            
CONECT    3    2    2    4                                            
CONECT    4    3    5    5                                            
CONECT    5    4    4    6                                            
CONECT    6    1    1    5    7                                       
CONECT    7    6                                                      
MASTER        0    0    0    0    0    0    0    0    7    0    7    0
END
""" % dict(VERSION=VERSION))

# exyz -- Extended XYZ cartesian coordinates format
class TestEXYZ(unittest.TestCase, WriteMixin):
    fmt = "exyz"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
7
phenol %PBC
   C        1.58460        -0.02490         0.00000
   C        1.57030         0.97550         0.00000
   C        2.42950         1.48820         0.00000
   C        3.30310         1.00040         0.00000
   C        3.31750        -0.00000         0.00000
   C        0.00000         0.00000         0.00000
   O       -1.00050         0.00510         0.00000

Vector1        1.00000         0.00000         0.00000
Vector2        0.00000         1.00000         0.00000
Vector3        0.00000         0.00000         1.00000
Offset         0.00000         0.00000         0.00000
""")

## # fa -- FASTA format
## XXX need a protein
## class TestFA(unittest.TestCase, WriteMixin):
##     fmt = "fa"
##     maxDiff = None
##     def test_default(self):
##         self.assertWriters(self.fmt, """\
## """)

## # fasta -- FASTA format
## XXX need a protein
## class TestFASTA(unittest.TestCase, WriteMixin):
##     fmt = "fasta"
##     maxDiff = None
##     def test_default(self):
##         self.assertWriters(self.fmt, """\
## """)

# feat -- Feature format
class TestFEAT(unittest.TestCase, WriteMixin):
    fmt = "feat"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
7
phenol
C    1.58460  -0.02490   0.00000 
C    1.57030   0.97550   0.00000 
C    2.42950   1.48820   0.00000 
C    3.30310   1.00040   0.00000 
C    3.31750  -0.00000   0.00000 
C    0.00000   0.00000   0.00000 
O   -1.00050   0.00510   0.00000 
""")

# fh -- Fenske-Hall Z-Matrix format [Write-only]
class TestFH(unittest.TestCase, WriteMixin):
    fmt = "fh"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\

7
C   1
C   1 1.001
C   2 1.001  1 120.006
C   3 1.001  2 119.997  1  -0.0
C   4 1.001  3 120.003  2  -0.0
C   1 1.585  2  88.281  3 180.0
O   2 2.748  1  70.139  3 180.0
""")

# fhiaims -- FHIaims XYZ format
class TestFHIAIMS(unittest.TestCase, WriteMixin):
    fmt = "fhiaims"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
#
# phenol
# Generated by Open Babel %(VERSION)s
#
atom                   1.58460       -0.02490        0.00000  C
atom                   1.57030        0.97550        0.00000  C
atom                   2.42950        1.48820        0.00000  C
atom                   3.30310        1.00040        0.00000  C
atom                   3.31750       -0.00000        0.00000  C
atom                   0.00000        0.00000        0.00000  C
atom                  -1.00050        0.00510        0.00000  O
""" % dict(VERSION=VERSION))

# fix -- SMILES FIX format [Write-only]
class TestFIX(unittest.TestCase, WriteMixin):
    fmt = "fix"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
c1ccccc1O
""")

# fps -- FPS text fingerprint format (Dalke) [Write-only]
_fps_date_pat = re.compile("#date=[0-9T:-]+")
_fps_type_version_pat = re.compile("(#type=[^/]+/)[0-9A-Za-z.]+")
_fps_software_version_pat = re.compile("(#software=OpenBabel/)[0-9A-Za-z.]+")
def normalize_fps(content):
    content = _fps_date_pat.sub("#date=Right now", content)
    content = _fps_type_version_pat.sub(r"\1test", content)
    content = _fps_software_version_pat.sub(r"\1test", content)
    return content

class TestFPS(unittest.TestCase, WriteMixin):
    fmt = "fps"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
#FPS1
#num_bits=1021
#type=OpenBabel-FP2/1
#software=OpenBabel/2.4.90
#source=
#date=2019-01-15T15:10:11
0000000000000000000002000000000000000000000000000000000000000000000000000000000000000008000000000000020000000000000000000000000008000000000000000000000002000000008000000000000040080000000000000000000000000002000000000000000000020000000000200800000000000000\tphenol
""", normalize=normalize_fps)

    def test_multimol_default(self):
        # Test that the header is written once, rather than once per molecule.
        phenol = get_mol(self, None)
        ethane = get_mol(self, "CC ethane")
        temp_file_object = tempfile.NamedTemporaryFile(suffix=".fps")
        filename = temp_file_object.name
        
        conv = get_converter(self, "fps")
        self.assertTrue(conv.WriteFile(phenol, filename))
        self.assertEqual(conv.GetOutputIndex(), 1)
        self.assertTrue(conv.Write(ethane))
        self.assertEqual(conv.GetOutputIndex(), 2)
        conv.CloseOutFile()

        with open(filename) as f:
            # Ensure there is a header
            line = f.readline()
            self.assertEqual(line[:5], "#FPS1", line)
            # Skip the rest of the header
            ids = []
            inheader = True
            for line in f:
                if line[:1] != "#":
                    inheader = False
                    self.assertEqual(line.count("\t"), 1, "Wrong number of fields?: %r" % (line,))
                    hex_fp, mid = line.rstrip("\n").split("\t", 1)
                    ids.append(mid)
                elif not inheader:
                    self.fail("Second header?: %r" % (line,))
            if inheader:
                self.fail("Reached end of file too early, after: %r" % (line,))

            self.assertEqual(ids, ["phenol", "ethane"])
        

    def test_MACCS(self):
        self.assertWriters(self.fmt, """\
#FPS1
#num_bits=166
#type=OpenBabel-MACCS/1
#software=OpenBabel/2.4.90
#source=
#date=2019-01-15T15:10:11
00000000000000000000000000000140004480101e\tphenol
""", normalize=normalize_fps, options={"f": "MACCS"})

    def test_FP2(self):
        self.assertWriters(self.fmt, """\
#FPS1
#num_bits=1021
#type=OpenBabel-FP2/1
#software=OpenBabel/2.4.90
#source=
#date=2019-01-15T15:10:11
0000000000000000000002000000000000000000000000000000000000000000000000000000000000000008000000000000020000000000000000000000000008000000000000000000000002000000008000000000000040080000000000000000000000000002000000000000000000020000000000200800000000000000\tphenol
""", normalize=normalize_fps, options={"f": "FP2"})

    def test_FP3(self):
        self.assertWriters(self.fmt, """\
#FPS1
#num_bits=55
#type=OpenBabel-FP3/1
#software=OpenBabel/2.4.90
#source=
#date=2019-01-15T15:10:11
0000000402b001\tphenol
""", normalize=normalize_fps, options={"f": "FP3"})

    def test_FP4(self):
        self.assertWriters(self.fmt, """\
#FPS1
#num_bits=307
#type=OpenBabel-FP4/1
#software=OpenBabel/2.4.90
#source=
#date=2019-01-15T15:10:11
000000000000000000000000000000000000000000010000000000000000000000000200400000\tphenol
""", normalize=normalize_fps, options={"f": "FP4"})

        
# fpt -- Fingerprint format [Write-only]
class TestFPT(unittest.TestCase, WriteMixin):
    fmt = "fpt"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
>phenol   12 bits set 
00000000 00000008 20000000 00000200 00000000 00000000 
02000000 00000000 00000000 00000840 00000000 00008000 
00000002 00000000 00000000 00000008 00000000 00000000 
00000000 00020000 00000000 08000000 00000000 00000000 
00000000 00000000 00000000 00000000 00000000 00020000 
00000000 00000000 
""")

    def test_MACCS(self):
        self.assertWriters(self.fmt, """\
>phenol   10 bits set 
00000000 00000000 0000001e 10804400 40010000 00000000 
00000000 00000000 
""", options={"f": "MACCS"})
        
    def test_describe_set_MACCS_bits(self):
        self.assertWriters(self.fmt, """\
>phenol
113: Onot%A%A\t127: A$A!O > 1 (&...) *2\t139: OH\t143: A$A!O\t152: OC(C)C\t157: C-O\t162: Aromatic\t163: 6M Ring\t164: O\t165: Ring	
""", options={"f": "MACCS", "s": None})


# fract -- Free Form Fractional format
class TestFRACT(unittest.TestCase, WriteMixin):
    fmt = "fract"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
phenol
   1.00000   1.00000   1.00000  90.00000  90.00000  90.00000
C    1.58460  -0.02490   0.00000
C    1.57030   0.97550   0.00000
C    2.42950   1.48820   0.00000
C    3.30310   1.00040   0.00000
C    3.31750  -0.00000   0.00000
C    0.00000   0.00000   0.00000
O   -1.00050   0.00510   0.00000

""")

## # fs -- Fastsearch format
## XXX "Not a valid output forma"
## class TestFS(unittest.TestCase, WriteMixin):
##     fmt = "fs"
##     maxDiff = None
##     def test_default(self):
##         self.assertWriters(self.fmt, """\
## """)

## # fsa -- FASTA format
## XXX need a protein structure
## class TestFSA(unittest.TestCase, WriteMixin):
##     fmt = "fsa"
##     maxDiff = None
##     def test_default(self):
##         self.assertWriters(self.fmt, """\
## """)

# gamin -- GAMESS Input
class TestGAMIN(unittest.TestCase, WriteMixin):
    fmt = "gamin"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
 $CONTRL COORD=CART UNITS=ANGS $END

 $DATA
phenol
C1
C      6.0      1.5846000000   -0.0249000000    0.0000000000 
C      6.0      1.5703000000    0.9755000000    0.0000000000 
C      6.0      2.4295000000    1.4882000000    0.0000000000 
C      6.0      3.3031000000    1.0004000000    0.0000000000 
C      6.0      3.3175000000   -0.0000000000    0.0000000000 
C      6.0      0.0000000000    0.0000000000    0.0000000000 
O      8.0     -1.0005000000    0.0051000000    0.0000000000 
 $END


""")

# gau -- Gaussian 98/03 Input [Write-only]
class TestGAU(unittest.TestCase, WriteMixin):
    fmt = "gau"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
!Put Keywords Here, check Charge and Multiplicity.
#

 phenol

0  1
C           1.58460        -0.02490         0.00000
C           1.57030         0.97550         0.00000
C           2.42950         1.48820         0.00000
C           3.30310         1.00040         0.00000
C           3.31750        -0.00000         0.00000
C           0.00000         0.00000         0.00000
O          -1.00050         0.00510         0.00000

""")

# gjc -- Gaussian 98/03 Input [Write-only]
class TestGJC(unittest.TestCase, WriteMixin):
    fmt = "gjc"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
!Put Keywords Here, check Charge and Multiplicity.
#

 phenol

0  1
C           1.58460        -0.02490         0.00000
C           1.57030         0.97550         0.00000
C           2.42950         1.48820         0.00000
C           3.30310         1.00040         0.00000
C           3.31750        -0.00000         0.00000
C           0.00000         0.00000         0.00000
O          -1.00050         0.00510         0.00000

""")

# gjf -- Gaussian 98/03 Input [Write-only]
class TestGJF(unittest.TestCase, WriteMixin):
    fmt = "gjf"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
!Put Keywords Here, check Charge and Multiplicity.
#

 phenol

0  1
C           1.58460        -0.02490         0.00000
C           1.57030         0.97550         0.00000
C           2.42950         1.48820         0.00000
C           3.30310         1.00040         0.00000
C           3.31750        -0.00000         0.00000
C           0.00000         0.00000         0.00000
O          -1.00050         0.00510         0.00000

""")

# gpr -- Ghemical format
class TestGPR(unittest.TestCase, WriteMixin):
    fmt = "gpr"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
!Header gpr 100
!Info 1
!Atoms 7
0 6
1 6
2 6
3 6
4 6
5 6
6 8
!Bonds 7
0 5 C
0 1 C
1 2 C
2 3 C
3 4 C
4 5 C
5 6 S
!Coord
0 0.15846 -0.00249 0
1 0.15703 0.09755 0
2 0.24295 0.14882 0
3 0.33031 0.10004 0
4 0.33175 -0 0
5 0 0 0
6 -0.10005 0.00051 0
!Charges
0 0.0420281
1 0.00328151
2 0.000205843
3 0.00328151
4 0.0420281
5 0.195745
6 -0.28657
!End
""")

# gr96 -- GROMOS96 format [Write-only]
class TestGR96(unittest.TestCase, WriteMixin):
    fmt = "gr96"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
#GENERATED BY OPEN BABEL %(VERSION)s
TITLE
phenol
END
POSITION
    1   UNL    C       1         1.58460        -0.02490         0.00000
    1   UNL    C       2         1.57030         0.97550         0.00000
    1   UNL    C       3         2.42950         1.48820         0.00000
    1   UNL    C       4         3.30310         1.00040         0.00000
    1   UNL    C       5         3.31750        -0.00000         0.00000
    1   UNL    C       6         0.00000         0.00000         0.00000
    1   UNL    O       7        -1.00050         0.00510         0.00000
END
""" % dict(VERSION=VERSION))

# gro -- GRO format
class TestGRO(unittest.TestCase, WriteMixin):
    fmt = "gro"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
phenol
7
    1UNL      C    1   0.158  -0.002   0.000
    1UNL      C    2   0.157   0.098   0.000
    1UNL      C    3   0.243   0.149   0.000
    1UNL      C    4   0.330   0.100   0.000
    1UNL      C    5   0.332  -0.000   0.000
    1UNL      C    6   0.000   0.000   0.000
    1UNL      O    7  -0.100   0.001   0.000
   0.00000   0.00000   0.00000
""")

# gukin -- GAMESS-UK Input
class TestGUKIN(unittest.TestCase, WriteMixin):
    fmt = "gukin"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
title
phenol

#
# NB: Class I directives (e.g. memory, multiplicity, charge etc) go here
#
# For more information see: http://www.cfs.dl.ac.uk/docs/index.shtml
#

geometry angstrom
     1.58460000     -0.02490000      0.00000000   6   C
     1.57030000      0.97550000      0.00000000   6   C
     2.42950000      1.48820000      0.00000000   6   C
     3.30310000      1.00040000      0.00000000   6   C
     3.31750000     -0.00000000      0.00000000   6   C
     0.00000000      0.00000000      0.00000000   6   C
    -1.00050000      0.00510000      0.00000000   8   O
end


basis 6-31G

#
# NB: Class II directives go here
#
# To perform a dft calculation with b3lyp and medium quadrature uncomment the below
# dft b3lyp
# dft quadrature medium
#

runtype scf

enter
""")

## # gukout -- GAMESS-UK Output
## XXX "Not a valid output format"
## class TestGUKOUT(unittest.TestCase, WriteMixin):
##     fmt = "gukout"
##     maxDiff = None
##     def test_default(self):
##         self.assertWriters(self.fmt, """\
## """)

# gzmat -- Gaussian Z-Matrix Input
class TestGZMAT(unittest.TestCase, WriteMixin):
    fmt = "gzmat"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
!Put Keywords Here, check Charge and Multiplicity.
#

 phenol

0  1
C
C  1  r2
C  2  r3  1  a3
C  3  r4  2  a4  1  d4
C  4  r5  3  a5  2  d5
C  1  r6  2  a6  3  d6
O  2  r7  1  a7  3  d7
Variables:
r2= 1.0005
r3= 1.0005
a3= 120.01
r4= 1.0006
a4= 120.00
d4=  -0.00
r5= 1.0005
a5= 120.00
d5=  -0.00
r6= 1.5848
a6=  88.28
d6= 180.00
r7= 2.7479
a7=  70.14
d7= 180.00

""")

# hin -- HyperChem HIN format
class TestHIN(unittest.TestCase, WriteMixin):
    fmt = "hin"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
mol 1 "phenol"
atom 1 - C   **  -  0.04203  1.58460  -0.02490   0.00000 2 6 a 2 a 
atom 2 - C   **  -  0.00328  1.57030   0.97550   0.00000 2 1 a 3 a 
atom 3 - C   **  -  0.00021  2.42950   1.48820   0.00000 2 2 a 4 a 
atom 4 - C   **  -  0.00328  3.30310   1.00040   0.00000 2 3 a 5 a 
atom 5 - C   **  -  0.04203  3.31750  -0.00000   0.00000 2 4 a 6 a 
atom 6 - C   **  -  0.19575  0.00000   0.00000   0.00000 3 1 a 5 a 7 s 
atom 7 - O   **  - -0.28657 -1.00050   0.00510   0.00000 1 6 s 
endmol 1
""")

# inchi -- InChI format
class TestINCHI(unittest.TestCase, WriteMixin):
    fmt = "inchi"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
InChI=1S/C6H6O/c7-6-4-2-1-3-5-6/h1-5,7H
""")

# inchikey -- InChIKey [Write-only]
class TestINCHIKEY(unittest.TestCase, WriteMixin):
    fmt = "inchikey"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
ISWSIDIOOBJBQZ-UHFFFAOYSA-N
""")

# inp -- GAMESS Input
class TestINP(unittest.TestCase, WriteMixin):
    fmt = "inp"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
 $CONTRL COORD=CART UNITS=ANGS $END

 $DATA
phenol
C1
C      6.0      1.5846000000   -0.0249000000    0.0000000000 
C      6.0      1.5703000000    0.9755000000    0.0000000000 
C      6.0      2.4295000000    1.4882000000    0.0000000000 
C      6.0      3.3031000000    1.0004000000    0.0000000000 
C      6.0      3.3175000000   -0.0000000000    0.0000000000 
C      6.0      0.0000000000    0.0000000000    0.0000000000 
O      8.0     -1.0005000000    0.0051000000    0.0000000000 
 $END


""")

# jin -- Jaguar input format
class TestJIN(unittest.TestCase, WriteMixin):
    fmt = "jin"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
phenol

&gen
&
&zmat
  C1      1.5846000    -0.0249000     0.0000000
  C2      1.5703000     0.9755000     0.0000000
  C3      2.4295000     1.4882000     0.0000000
  C4      3.3031000     1.0004000     0.0000000
  C5      3.3175000    -0.0000000     0.0000000
  C6      0.0000000     0.0000000     0.0000000
  O7     -1.0005000     0.0051000     0.0000000
&
""")

# k -- Compare molecules using InChI [Write-only]
class TestK(unittest.TestCase, WriteMixin):
    fmt = "k"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
InChI=1S/C6H6O/c7-6-4-2-1-3-5-6/h1-5,7H phenol
""")

# lmpdat -- The LAMMPS data format [Write-only]
class TestLMPDAT(unittest.TestCase, WriteMixin):
    fmt = "lmpdat"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
LAMMPS data file generated by OpenBabel
7 atoms
7 bonds
8 angles
8 dihedrals
0 impropers
2 atom types
2 bond types
2 angle types
2 dihedral types
0 improper types
  -1.50050    3.81750 xlo xhi
  -1.50050    3.81750 ylo yhi
  -1.50050    3.81750 zlo zhi



Masses

1 12.0107 # C
2 15.9994 # O


Atoms

1       1    1    0.00000    1.58460   -0.02490    0.00000 #   C
2       1    1    0.00000    1.57030    0.97550    0.00000 #   C
3       1    1    0.00000    2.42950    1.48820    0.00000 #   C
4       1    1    0.00000    3.30310    1.00040    0.00000 #   C
5       1    1    0.00000    3.31750   -0.00000    0.00000 #   C
6       1    1    0.00000    0.00000    0.00000    0.00000 #   C
7       1    2   -0.82000   -1.00050    0.00510    0.00000 #   O


Bonds

1       1    1    6 #  C: C
2       1    1    2 #  C: C
3       1    2    3 #  C: C
4       1    3    4 #  C: C
5       1    4    5 #  C: C
6       1    5    6 #  C: C
7       2    7    6 #  O: C


Angles

1       1    6    1    2 #  C: C: C
2       1    3    2    1 #  C: C: C
3       1    4    3    2 #  C: C: C
4       1    5    4    3 #  C: C: C
5       1    6    5    4 #  C: C: C
6       1    5    6    1 #  C: C: C
7       2    7    6    1 #  O: C: C
8       2    7    6    5 #  O: C: C


Dihedrals

1       1    5    1    6    2 #  C: C: C: C
2       2    7    1    6    2 #  O: C: C: C
3       1    3    1    2    6 #  C: C: C: C
4       1    4    2    3    1 #  C: C: C: C
5       1    5    3    4    2 #  C: C: C: C
6       1    6    4    5    3 #  C: C: C: C
7       1    1    5    6    4 #  C: C: C: C
8       2    7    5    6    4 #  O: C: C: C
""")

## # lpmd -- LPMD format
## XXX "The original file doesn't have the information about the unitcell"
## class TestLPMD(unittest.TestCase, WriteMixin):
##     fmt = "lpmd"
##     maxDiff = None
##     def test_default(self):
##         self.assertWriters(self.fmt, """\
## """)

# mcdl -- MCDL format
class TestMCDL(unittest.TestCase, WriteMixin):
    fmt = "mcdl"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
C;5CH;OH[2,3,7;4;5;6;6]{CN:}phenol}
""")

# mcif -- Macromolecular Crystallographic Info
class TestMCIF(unittest.TestCase, WriteMixin):
    fmt = "mcif"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
# --------------------------------------------------------------------------
#
# CIF file generated by openbabel %(VERSION)s http://openbabel.org/
# to comply with the Macromolecular CIF Dictionary  (cif_mm.dic) version  2.0.11 http://mmcif.pdb.org/
# The contents of this file were derived from 
#
#---------------------------------------------------------------------------

data_PHENOL

###########
## ENTRY ##
###########

_entry.id	PHENOL

##############
## CHEMICAL ##
##############

_chemical.entry_id	PHENOL
_chemical.name_common	'phenol'

######################
## CHEMICAL FORMULA ##
######################

_chemical_formula.entry_id	PHENOL
_chemical_formula.structural	'C6H6O'

###############
## ATOM_SITE ##
###############

loop_
_atom_site.id
_atom_site.type_symbol
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
	1	C	1.5846	-0.0249	0
	2	C	1.5703	0.9755	0
	3	C	2.4295	1.4882	0
	4	C	3.3031	1.0004	0
	5	C	3.3175	-0	0
	6	C	0	0	0
	7	O	-1.0005	0.0051	0

""" % dict(VERSION=VERSION))

# MDFF -- MDFF format
class TestMDFF(unittest.TestCase, WriteMixin):
    fmt = "MDFF"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
7
phenol
0.0  0.0  0.0
0.0  0.0  0.0
0.0  0.0  0.0
2
C   O   
6   1   
Cartesian
C        1.5846000000000000085     -0.0248999999999999985      0.0000000000000000000
C        1.5703000000000000291      0.9755000000000000338      0.0000000000000000000
C        2.4294999999999999929      1.4881999999999999673      0.0000000000000000000
C        3.3031000000000001471      1.0003999999999999559      0.0000000000000000000
C        3.3174999999999998934     -0.0000000000000000000      0.0000000000000000000
C        0.0000000000000000000      0.0000000000000000000      0.0000000000000000000
O       -1.0004999999999999449      0.0051000000000000004      0.0000000000000000000
""")

# Normalize MDL formats by removing the timestamp from the string
_sd_timestamp_pat_u = re.compile(u"OpenBabel\\d\\d\\d\\d\\d\\d\\d\\d\\d\\d")
_sd_timestamp_pat_b = re.compile(b"OpenBabel\\d\\d\\d\\d\\d\\d\\d\\d\\d\\d")
def normalize_sd_timestamp(data):
    if isinstance(data, type(b"")):
        return _sd_timestamp_pat_b.sub(b"OpenBabel2020202020", data)
    else:
        return _sd_timestamp_pat_u.sub(u"OpenBabel2020202020", data)

        
# mdl -- MDL MOL format
class TestMDL(unittest.TestCase, WriteMixin):
    fmt = "mdl"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
phenol
 OpenBabel01151915422D

  7  7  0  0  0  0  0  0  0  0999 V2000
    1.5846   -0.0249    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5703    0.9755    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4295    1.4882    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3031    1.0004    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3175   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0005    0.0051    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  6  2  0  0  0  0
  1  2  1  0  0  0  0
  2  3  2  0  0  0  0
  3  4  1  0  0  0  0
  4  5  2  0  0  0  0
  5  6  1  0  0  0  0
  6  7  1  0  0  0  0
M  END
""", normalize=normalize_sd_timestamp)

# ml2 -- Sybyl Mol2 format
class TestML2(unittest.TestCase, WriteMixin):
    fmt = "ml2"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
@<TRIPOS>MOLECULE
phenol
 7 7 0 0 0
SMALL
GASTEIGER

@<TRIPOS>ATOM
      1 C           1.5846   -0.0249    0.0000 C.ar    1  UNL1        0.0420
      2 C           1.5703    0.9755    0.0000 C.ar    1  UNL1        0.0033
      3 C           2.4295    1.4882    0.0000 C.ar    1  UNL1        0.0002
      4 C           3.3031    1.0004    0.0000 C.ar    1  UNL1        0.0033
      5 C           3.3175   -0.0000    0.0000 C.ar    1  UNL1        0.0420
      6 C           0.0000    0.0000    0.0000 C.ar    1  UNL1        0.1957
      7 O          -1.0005    0.0051    0.0000 O.3     1  UNL1       -0.2866
@<TRIPOS>BOND
     1     1     6   ar
     2     1     2   ar
     3     2     3   ar
     4     3     4   ar
     5     4     5   ar
     6     5     6   ar
     7     6     7    1
""")

# mmcif -- Macromolecular Crystallographic Info
class TestMMCIF(unittest.TestCase, WriteMixin):
    fmt = "mmcif"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
# --------------------------------------------------------------------------
#
# CIF file generated by openbabel %(VERSION)s http://openbabel.org/
# to comply with the Macromolecular CIF Dictionary  (cif_mm.dic) version  2.0.11 http://mmcif.pdb.org/
# The contents of this file were derived from 
#
#---------------------------------------------------------------------------

data_PHENOL

###########
## ENTRY ##
###########

_entry.id	PHENOL

##############
## CHEMICAL ##
##############

_chemical.entry_id	PHENOL
_chemical.name_common	'phenol'

######################
## CHEMICAL FORMULA ##
######################

_chemical_formula.entry_id	PHENOL
_chemical_formula.structural	'C6H6O'

###############
## ATOM_SITE ##
###############

loop_
_atom_site.id
_atom_site.type_symbol
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
	1	C	1.5846	-0.0249	0
	2	C	1.5703	0.9755	0
	3	C	2.4295	1.4882	0
	4	C	3.3031	1.0004	0
	5	C	3.3175	-0	0
	6	C	0	0	0
	7	O	-1.0005	0.0051	0

""" % dict(VERSION=VERSION))

# mmd -- MacroModel format
class TestMMD(unittest.TestCase, WriteMixin):
    fmt = "mmd"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
     7 phenol      E =   0.000 KJ/mol
   2     6 2     2 1     0 0     0 0     0 0     0 0    1.584600   -0.024900    0.000000     0     0  0.04203 
   2     1 1     3 2     0 0     0 0     0 0     0 0    1.570300    0.975500    0.000000     0     0  0.00328 
   2     2 2     4 1     0 0     0 0     0 0     0 0    2.429500    1.488200    0.000000     0     0  0.00021 
   2     3 1     5 2     0 0     0 0     0 0     0 0    3.303100    1.000400    0.000000     0     0  0.00328 
   2     4 2     6 1     0 0     0 0     0 0     0 0    3.317500   -0.000000    0.000000     0     0  0.04203 
   2     1 2     5 1     7 1     0 0     0 0     0 0    0.000000    0.000000    0.000000     0     0  0.19575 
  16     6 1     0 0     0 0     0 0     0 0     0 0   -1.000500    0.005100    0.000000     0     0 -0.28657 
""")

# mmod -- MacroModel format
class TestMMOD(unittest.TestCase, WriteMixin):
    fmt = "mmod"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
     7 phenol      E =   0.000 KJ/mol
   2     6 2     2 1     0 0     0 0     0 0     0 0    1.584600   -0.024900    0.000000     0     0  0.04203 
   2     1 1     3 2     0 0     0 0     0 0     0 0    1.570300    0.975500    0.000000     0     0  0.00328 
   2     2 2     4 1     0 0     0 0     0 0     0 0    2.429500    1.488200    0.000000     0     0  0.00021 
   2     3 1     5 2     0 0     0 0     0 0     0 0    3.303100    1.000400    0.000000     0     0  0.00328 
   2     4 2     6 1     0 0     0 0     0 0     0 0    3.317500   -0.000000    0.000000     0     0  0.04203 
   2     1 2     5 1     7 1     0 0     0 0     0 0    0.000000    0.000000    0.000000     0     0  0.19575 
  16     6 1     0 0     0 0     0 0     0 0     0 0   -1.000500    0.005100    0.000000     0     0 -0.28657 
""")

# mna -- Multilevel Neighborhoods of Atoms (MNA) [Write-only]
class TestMNA(unittest.TestCase, WriteMixin):
    fmt = "mna"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
# The contents of this file were derived from 
# Title = phenol
C(C(CC-O)C(CC-H)-H(C))
C(C(CC-H)C(CC-H)-H(C))
C(C(CC-H)C(CC-H)-H(C))
C(C(CC-H)C(CC-H)-H(C))
C(C(CC-H)C(CC-O)-H(C))
C(C(CC-H)C(CC-H)-O(C-H))
-O(C(CC-O)-H(-O))
-H(C(CC-H))
-H(C(CC-H))
-H(C(CC-H))
-H(C(CC-H))
-H(C(CC-H))
-H(-O(C-H))
""")

# mol -- MDL MOL format
class TestMOL(unittest.TestCase, WriteMixin):
    fmt = "mol"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
phenol
 OpenBabel01151915412D

  7  7  0  0  0  0  0  0  0  0999 V2000
    1.5846   -0.0249    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5703    0.9755    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4295    1.4882    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3031    1.0004    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3175   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0005    0.0051    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  6  2  0  0  0  0
  1  2  1  0  0  0  0
  2  3  2  0  0  0  0
  3  4  1  0  0  0  0
  4  5  2  0  0  0  0
  5  6  1  0  0  0  0
  6  7  1  0  0  0  0
M  END
""", normalize=normalize_sd_timestamp)

# mol2 -- Sybyl Mol2 format
class TestMOL2(unittest.TestCase, WriteMixin):
    fmt = "mol2"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
@<TRIPOS>MOLECULE
phenol
 7 7 0 0 0
SMALL
GASTEIGER

@<TRIPOS>ATOM
      1 C           1.5846   -0.0249    0.0000 C.ar    1  UNL1        0.0420
      2 C           1.5703    0.9755    0.0000 C.ar    1  UNL1        0.0033
      3 C           2.4295    1.4882    0.0000 C.ar    1  UNL1        0.0002
      4 C           3.3031    1.0004    0.0000 C.ar    1  UNL1        0.0033
      5 C           3.3175   -0.0000    0.0000 C.ar    1  UNL1        0.0420
      6 C           0.0000    0.0000    0.0000 C.ar    1  UNL1        0.1957
      7 O          -1.0005    0.0051    0.0000 O.3     1  UNL1       -0.2866
@<TRIPOS>BOND
     1     1     6   ar
     2     1     2   ar
     3     2     3   ar
     4     3     4   ar
     5     4     5   ar
     6     5     6   ar
     7     6     7    1
""")

# mold -- Molden format
class TestMOLD(unittest.TestCase, WriteMixin):
    fmt = "mold"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
[Molden Format]
[Atoms] Angs
 C     1  6     1.584600    -0.024900     0.000000
 C     2  6     1.570300     0.975500     0.000000
 C     3  6     2.429500     1.488200     0.000000
 C     4  6     3.303100     1.000400     0.000000
 C     5  6     3.317500    -0.000000     0.000000
 C     6  6     0.000000     0.000000     0.000000
 O     7  8    -1.000500     0.005100     0.000000
""")

# molden -- Molden format
class TestMOLDEN(unittest.TestCase, WriteMixin):
    fmt = "molden"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
[Molden Format]
[Atoms] Angs
 C     1  6     1.584600    -0.024900     0.000000
 C     2  6     1.570300     0.975500     0.000000
 C     3  6     2.429500     1.488200     0.000000
 C     4  6     3.303100     1.000400     0.000000
 C     5  6     3.317500    -0.000000     0.000000
 C     6  6     0.000000     0.000000     0.000000
 O     7  8    -1.000500     0.005100     0.000000
""")

# molf -- Molden format
class TestMOLF(unittest.TestCase, WriteMixin):
    fmt = "molf"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
[Molden Format]
[Atoms] Angs
 C     1  6     1.584600    -0.024900     0.000000
 C     2  6     1.570300     0.975500     0.000000
 C     3  6     2.429500     1.488200     0.000000
 C     4  6     3.303100     1.000400     0.000000
 C     5  6     3.317500    -0.000000     0.000000
 C     6  6     0.000000     0.000000     0.000000
 O     7  8    -1.000500     0.005100     0.000000
""")

# molreport -- Open Babel molecule report [Write-only]
class TestMOLREPORT(unittest.TestCase, WriteMixin):
    fmt = "molreport"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
TITLE: phenol
FORMULA: C6H6O
MASS: 94.1112
ATOM:         1   C TYPE: Car    HYB:  2 CHARGE:   0.0420
ATOM:         2   C TYPE: Car    HYB:  2 CHARGE:   0.0033
ATOM:         3   C TYPE: Car    HYB:  2 CHARGE:   0.0002
ATOM:         4   C TYPE: Car    HYB:  2 CHARGE:   0.0033
ATOM:         5   C TYPE: Car    HYB:  2 CHARGE:   0.0420
ATOM:         6   C TYPE: Car    HYB:  2 CHARGE:   0.1957
ATOM:         7   O TYPE: O3     HYB:  2 CHARGE:  -0.2866
BOND:         0 START:         1 END:         6 ORDER:   2
BOND:         1 START:         1 END:         2 ORDER:   1
BOND:         2 START:         2 END:         3 ORDER:   2
BOND:         3 START:         3 END:         4 ORDER:   1
BOND:         4 START:         4 END:         5 ORDER:   2
BOND:         5 START:         5 END:         6 ORDER:   1
BOND:         6 START:         6 END:         7 ORDER:   1
""")

# mop -- MOPAC Cartesian format
class TestMOP(unittest.TestCase, WriteMixin):
    fmt = "mop"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
PUT KEYWORDS HERE
phenol

C   1.58460 1 -0.02490 1  0.00000 1
C   1.57030 1  0.97550 1  0.00000 1
C   2.42950 1  1.48820 1  0.00000 1
C   3.30310 1  1.00040 1  0.00000 1
C   3.31750 1 -0.00000 1  0.00000 1
C   0.00000 1  0.00000 1  0.00000 1
O  -1.00050 1  0.00510 1  0.00000 1
""")

# mopcrt -- MOPAC Cartesian format
class TestMOPCRT(unittest.TestCase, WriteMixin):
    fmt = "mopcrt"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
PUT KEYWORDS HERE
phenol

C   1.58460 1 -0.02490 1  0.00000 1
C   1.57030 1  0.97550 1  0.00000 1
C   2.42950 1  1.48820 1  0.00000 1
C   3.30310 1  1.00040 1  0.00000 1
C   3.31750 1 -0.00000 1  0.00000 1
C   0.00000 1  0.00000 1  0.00000 1
O  -1.00050 1  0.00510 1  0.00000 1
""")

# mopin -- MOPAC Internal
class TestMOPIN(unittest.TestCase, WriteMixin):
    fmt = "mopin"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
PUT KEYWORDS HERE
phenol

C    0.000000  1    0.000000  1    0.000000  1     0   0   0
C    1.000502  1    0.000000  1    0.000000  1     1   0   0
C    1.000543  1  120.006337  1    0.000000  1     2   1   0
C    1.000563  1  119.996638  1   -0.000000  1     3   2   1
C    1.000504  1  120.002751  1   -0.000000  1     4   3   2
C    1.584796  1   88.280797  1  180.000000  1     1   2   3
O    2.747852  1   70.138927  1  180.000000  1     2   1   3
""")

# mp -- Molpro input format [Write-only]
class TestMP(unittest.TestCase, WriteMixin):
    fmt = "mp"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
*** phenol
!file,2,INSERT WAVEFUNCTION FILE LOCATION HERE
!memory,INSERT MEMORY HERE
!basis,INSERT BASIS SET HERE

geomtyp=xyz
geometry={
7
Geometry specification:
  C,        1.58460,       -0.02490,        0.00000
  C,        1.57030,        0.97550,        0.00000
  C,        2.42950,        1.48820,        0.00000
  C,        3.30310,        1.00040,        0.00000
  C,        3.31750,       -0.00000,        0.00000
  C,        0.00000,        0.00000,        0.00000
  O,       -1.00050,        0.00510,        0.00000
}

!INSERT QM METHODS HERE
!hf
---
""")

# mpc -- MOPAC Cartesian format
class TestMPC(unittest.TestCase, WriteMixin):
    fmt = "mpc"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
PUT KEYWORDS HERE
phenol

C   1.58460 1 -0.02490 1  0.00000 1
C   1.57030 1  0.97550 1  0.00000 1
C   2.42950 1  1.48820 1  0.00000 1
C   3.30310 1  1.00040 1  0.00000 1
C   3.31750 1 -0.00000 1  0.00000 1
C   0.00000 1  0.00000 1  0.00000 1
O  -1.00050 1  0.00510 1  0.00000 1
""")

# mpd -- MolPrint2D format [Write-only]
class TestMPD(unittest.TestCase, WriteMixin):
    fmt = "mpd"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
phenol	3;1-2-3;2-2-3;2-1-8;	3;1-2-3;2-2-3;	3;1-2-3;2-2-3;	3;1-2-3;2-2-3;	3;1-2-3;2-2-3;2-1-8;	3;1-2-3;1-1-8;2-2-3;	8;1-1-3;2-2-3;\t
""")

# mpqcin -- MPQC simplified input format [Write-only]
class TestMPQCIN(unittest.TestCase, WriteMixin):
    fmt = "mpqcin"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
% phenol

molecule:
   C   1.58460  -0.02490   0.00000 
   C   1.57030   0.97550   0.00000 
   C   2.42950   1.48820   0.00000 
   C   3.30310   1.00040   0.00000 
   C   3.31750  -0.00000   0.00000 
   C   0.00000   0.00000   0.00000 
   O  -1.00050   0.00510   0.00000 



""")

# mrv -- Chemical Markup Language
class TestMRV(unittest.TestCase, WriteMixin):
    fmt = "mrv"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
<?xml version="1.0"?>
<molecule id="phenol" xmlns="http://www.xml-cml.org/schema">
 <atomArray>
  <atom id="a1" elementType="C" hydrogenCount="1" x2="1.584600" y2="-0.024900"/>
  <atom id="a2" elementType="C" hydrogenCount="1" x2="1.570300" y2="0.975500"/>
  <atom id="a3" elementType="C" hydrogenCount="1" x2="2.429500" y2="1.488200"/>
  <atom id="a4" elementType="C" hydrogenCount="1" x2="3.303100" y2="1.000400"/>
  <atom id="a5" elementType="C" hydrogenCount="1" x2="3.317500" y2="-0.000000"/>
  <atom id="a6" elementType="C" hydrogenCount="0" x2="0.000000" y2="0.000000"/>
  <atom id="a7" elementType="O" hydrogenCount="1" x2="-1.000500" y2="0.005100"/>
 </atomArray>
 <bondArray>
  <bond atomRefs2="a1 a6" order="2"/>
  <bond atomRefs2="a1 a2" order="1"/>
  <bond atomRefs2="a2 a3" order="2"/>
  <bond atomRefs2="a3 a4" order="1"/>
  <bond atomRefs2="a4 a5" order="2"/>
  <bond atomRefs2="a5 a6" order="1"/>
  <bond atomRefs2="a6 a7" order="1"/>
 </bondArray>
</molecule>
""")

# msms -- M.F. Sanner's MSMS input format [Write-only]
class TestMSMS(unittest.TestCase, WriteMixin):
    fmt = "msms"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
1.5846	-0.0249	0	1.7
1.5703	0.9755	0	1.7
2.4295	1.4882	0	1.7
3.3031	1.0004	0	1.7
3.3175	-0	0	1.7
0	0	0	1.7
-1.0005	0.0051	0	1.52
""")

# nul -- Outputs nothing [Write-only]
class TestNUL(unittest.TestCase, WriteMixin):
    fmt = "nul"
    maxDiff = None
    def test_default(self):
        # Why can't I write this to a file?
        self.assertWriteString(self.fmt, "")

# nw -- NWChem input format [Write-only]
class TestNW(unittest.TestCase, WriteMixin):
    fmt = "nw"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
start molecule

title 
 phenol

geometry units angstroms print xyz autosym
  C        1.58460       -0.02490        0.00000
  C        1.57030        0.97550        0.00000
  C        2.42950        1.48820        0.00000
  C        3.30310        1.00040        0.00000
  C        3.31750       -0.00000        0.00000
  C        0.00000        0.00000        0.00000
  O       -1.00050        0.00510        0.00000
end
""")

# orcainp -- ORCA input format [Write-only]
class TestORCAINP(unittest.TestCase, WriteMixin):
    fmt = "orcainp"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
# ORCA input file
# phenol
! insert inline commands here 
* xyz 0 1
   C        1.58460        -0.02490         0.00000
   C        1.57030         0.97550         0.00000
   C        2.42950         1.48820         0.00000
   C        3.30310         1.00040         0.00000
   C        3.31750        -0.00000         0.00000
   C        0.00000         0.00000         0.00000
   O       -1.00050         0.00510         0.00000
*
""")

# outmol -- DMol3 coordinates format
class TestOUTMOL(unittest.TestCase, WriteMixin):
    fmt = "outmol"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
$coordinates
C             2.99445980216940   -0.04705417712610    0.00000000000000
C             2.96743672052670    1.84342770226950    0.00000000000000
C             4.59108929027550    2.81229021682980    0.00000000000000
C             6.24195391426590    1.89048187939560    0.00000000000000
C             6.26916596850750   -0.00000000000000    0.00000000000000
C             0.00000000000000    0.00000000000000    0.00000000000000
O            -1.89067085199450    0.00963760254390    0.00000000000000
$end
""")

# paint -- Painter format [Write-only]
class TestPAINT(unittest.TestCase, WriteMixin):
    fmt = "paint"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
NewCanvas 202.1 122.8
SetPenColor 0.0 0.0 0.0 1.0 (rgba)
SetPenColor 0.0 0.0 0.0 1.0 (rgba)
SetPenColor 0.0 0.0 0.0 1.0 (rgba)
SetPenColor 0.0 0.0 0.0 1.0 (rgba)
SetPenColor 0.0 0.0 0.0 1.0 (rgba)
SetPenColor 0.0 0.0 0.0 1.0 (rgba)
SetPenColor 0.0 0.0 0.0 1.0 (rgba)
DrawLine 68.3 82.1 to 53.0 82.0
SetPenColor 0.0 0.0 0.0 1.0 (rgba)
DrawLine 137.0 40.0 to 161.7 53.8
SetPenColor 0.0 0.0 0.0 1.0 (rgba)
DrawLine 161.7 53.8 to 162.1 82.1
DrawLine 154.5 59.9 to 154.8 76.2
SetPenColor 0.0 0.0 0.0 1.0 (rgba)
DrawLine 162.1 82.1 to 68.3 82.1
SetPenColor 0.0 0.0 0.0 1.0 (rgba)
DrawLine 68.3 82.1 to 113.1 82.8
DrawLine 74.4 75.0 to 107.2 75.5
SetPenColor 0.0 0.0 0.0 1.0 (rgba)
DrawLine 113.1 82.8 to 112.7 54.5
SetPenColor 0.0 0.0 0.0 1.0 (rgba)
DrawLine 112.7 54.5 to 137.0 40.0
DrawLine 121.5 57.6 to 135.5 49.3
SetPenColor 0.4 0.4 0.4 1.0 (rgba)
SetPenColor 0.4 0.4 0.4 1.0 (rgba)
SetPenColor 0.4 0.4 0.4 1.0 (rgba)
SetPenColor 0.4 0.4 0.4 1.0 (rgba)
SetPenColor 0.4 0.4 0.4 1.0 (rgba)
SetPenColor 0.4 0.4 0.4 1.0 (rgba)
SetPenColor 1.0 0.1 0.1 1.0 (rgba)
SetFontSize 16
SetFontSize 16
SetFontSize 16
SetFontSize 16
DrawText 40.0 81.9 "HO"
""")

# pcjson -- PubChem JSON
class TestPCJSON(unittest.TestCase, WriteMixin):
    fmt = "pcjson"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
{
  "PC_Compounds": [
    {
      "atoms": {
        "aids": [
          1,
          2,
          3,
          4,
          5,
          6,
          7,
          8,
          9,
          10,
          11,
          12,
          13
        ],
        "element": [
          6,
          6,
          6,
          6,
          6,
          6,
          8,
          1,
          1,
          1,
          1,
          1,
          1
        ]
      },
      "bonds": {
        "aid1": [
          1,
          1,
          2,
          3,
          4,
          5,
          6,
          1,
          2,
          3,
          4,
          5,
          7
        ],
        "aid2": [
          6,
          2,
          3,
          4,
          5,
          6,
          7,
          8,
          9,
          10,
          11,
          12,
          13
        ],
        "order": [
          2,
          1,
          2,
          1,
          2,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1
        ]
      },
      "coords": [
        {
          "type": [
            1
          ],
          "aids": [
            1,
            2,
            3,
            4,
            5,
            6,
            7,
            8,
            9,
            10,
            11,
            12,
            13
          ],
          "conformers": [
            {
              "x": [
                1.5846,
                1.5703,
                2.4295,
                3.3031,
                3.3175,
                0.0,
                -1.0005,
                2.313816216007316,
                0.669250157347277,
                2.4146659588503769,
                4.189331679349326,
                4.052466878708012,
                -1.4648575597102012
              ],
              "y": [
                1.5846,
                1.5703,
                2.4295,
                3.3031,
                3.3175,
                0.0,
                -1.0005,
                2.313816216007316,
                0.669250157347277,
                2.4146659588503769,
                4.189331679349326,
                4.052466878708012,
                -1.4648575597102012
              ],
              "style": {
                "annotation": [
                  8,
                  8,
                  8,
                  8,
                  8,
                  8
                ],
                "aid1": [
                  1,
                  1,
                  2,
                  3,
                  4,
                  5
                ],
                "aid2": [
                  6,
                  2,
                  3,
                  4,
                  5,
                  6
                ]
              }
            }
          ]
        }
      ],
      "charge": 0
    }
  ]
}""")

# pcm -- PCModel Format
class TestPCM(unittest.TestCase, WriteMixin):
    fmt = "pcm"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
{PCM phenol
NA 7
ATOMTYPES 1
AT 1,40:1.5846,-0.0249,0 B 6,2 2,1 C 0.0420281
AT 2,40:1.5703,0.9755,0 B 1,1 3,2 C 0.00328151
AT 3,40:2.4295,1.4882,0 B 2,2 4,1 C 0.000205843
AT 4,40:3.3031,1.0004,0 B 3,1 5,2 C 0.00328151
AT 5,40:3.3175,-0,0 B 4,2 6,1 C 0.0420281
AT 6,40:0,0,0 B 1,2 5,1 7,1 C 0.195745
AT 7,6:-1.0005,0.0051,0 B 6,1 C -0.28657
}
""")

# pdb -- Protein Data Bank format
class TestPDB(unittest.TestCase, WriteMixin):
    fmt = "pdb"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
COMPND    phenol 
AUTHOR    GENERATED BY OPEN BABEL %(VERSION)s
HETATM    1  C   UNL     1       1.585  -0.025   0.000  1.00  0.00           C  
HETATM    2  C   UNL     1       1.570   0.976   0.000  1.00  0.00           C  
HETATM    3  C   UNL     1       2.429   1.488   0.000  1.00  0.00           C  
HETATM    4  C   UNL     1       3.303   1.000   0.000  1.00  0.00           C  
HETATM    5  C   UNL     1       3.317   0.000   0.000  1.00  0.00           C  
HETATM    6  C   UNL     1       0.000   0.000   0.000  1.00  0.00           C  
HETATM    7  O   UNL     1      -1.000   0.005   0.000  1.00  0.00           O  
CONECT    1    6    6    2                                            
CONECT    2    1    3    3                                            
CONECT    3    2    2    4                                            
CONECT    4    3    5    5                                            
CONECT    5    4    4    6                                            
CONECT    6    1    1    5    7                                       
CONECT    7    6                                                      
MASTER        0    0    0    0    0    0    0    0    7    0    7    0
END
""" % dict(VERSION=VERSION))

# pdbqt -- AutoDock PDBQT format
class TestPDBQT(unittest.TestCase, WriteMixin):
    fmt = "pdbqt"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
REMARK  Name = phenol
REMARK  0 active torsions:
REMARK  status: ('A' for Active; 'I' for Inactive)
REMARK                            x       y       z     vdW  Elec       q    Type
REMARK                         _______ _______ _______ _____ _____    ______ ____
ROOT
ATOM      1  C   UNL     1       1.585  -0.025   0.000  0.00  0.00    +0.000 A 
ATOM      2  C   UNL     1       1.570   0.976   0.000  0.00  0.00    +0.000 A 
ATOM      3  C   UNL     1       2.429   1.488   0.000  0.00  0.00    +0.000 A 
ATOM      4  C   UNL     1       3.303   1.000   0.000  0.00  0.00    +0.000 A 
ATOM      5  C   UNL     1       3.317   0.000   0.000  0.00  0.00    +0.000 A 
ATOM      6  C   UNL     1       0.000   0.000   0.000  0.00  0.00    +0.000 A 
ATOM      7  O   UNL     1      -1.000   0.005   0.000  0.00  0.00    +0.000 OA
ENDROOT
TORSDOF 0
""")

# png -- PNG 2D depiction
class TestPNG(unittest.TestCase, WriteMixin):
    fmt = "png"
    maxDiff = None
    def test_default(self):
        # This doesn't seem to work for a string?
        self.assertBinaryWriteFile(self.fmt, b"\x00\x00\x00\x18tEXtsmiles\x00c1ccccc1O\tphenol\nt\x82e\xc6")

## # pointcloud -- Point cloud on VDW surface [Write-only]
## class TestPOINTCLOUD(unittest.TestCase, WriteMixin):
##     fmt = "pointcloud"
##     maxDiff = None
##     def test_default(self):
##         self.assertWriters(self.fmt, """\
## 1.60287	0.366606	-1.6542
## 1.84042	0.0911282	1.67663
## 2.16155	0.0920641	-1.59482
## 2.24669	-1.17133	-1.06646
## 1.48693	-1.10539	-1.30882
## 1.06863	-1.56851	-0.490946
## 0.950812	-0.774568	1.38792
## 1.53266	-0.423934	1.65169
## 1.02418	0.0825103	1.60137
## 1.35345	-0.521534	-1.60932
## 2.03809	-1.53972	-0.624231
## 1.81315	-1.28599	1.11688
## 1.61398	-1.72433	-0.0329242
## 2.00427	-0.832704	1.43573
## 1.95914	-0.558616	-1.56999
## 0.856928	-0.923254	-1.24638
## 1.17686	-1.28285	1.06833
## 1.20053	1.21326	1.64218
## 2.21578	0.539197	-1.51096
## 0.503037	2.05994	-0.758254
## 0.915434	2.03276	-1.15903
## 1.49448	1.98584	-1.36509
## 0.677649	1.86202	1.14335
## 1.70666	1.01844	-1.69398
## 0.606207	0.940736	1.39976
## 0.800302	2.48089	-0.175754
## 2.30163	0.492979	1.45682
## 0.359758	2.12949	0.304793
## 0.172569	1.70773	0.6326
## 1.08273	2.32878	0.906044
## 1.20281	1.7605	1.46244
## 1.09029	2.41042	-0.774986
## 1.06925	0.58422	1.57666
## 0.819971	1.51412	-1.4272
## 2.24158	1.24937	1.67262
## 2.268	0.913847	-1.59187
## 3.21823	2.70485	-0.88751
## 1.69832	2.84941	-0.708855
## 3.37734	2.88568	0.196553
## 1.81476	2.94375	0.627274
## 2.45416	3.04162	-0.690131
## 1.39655	2.82919	-0.15739
## 2.78789	2.69877	1.13845
## 2.58101	2.10658	-1.57628
## 2.62307	1.99508	1.61109
## 1.91432	2.54924	1.22424
## 2.63768	3.10226	0.491385
## 3.37994	2.32048	1.13753
## 1.97981	2.37337	-1.37994
## 2.81402	1.57517	-1.65366
## 3.09355	2.34502	-1.30954
## 2.51748	2.67703	-1.212
## 4.50657	1.72936	0.954082
## 4.10815	1.38645	-1.44667
## 4.16586	2.45354	0.184494
## 2.76711	0.762946	-1.59572
## 3.58896	2.09447	-1.26937
## 3.63994	2.15084	1.20542
## 3.61046	1.63061	1.54867
## 4.21658	0.640438	1.3878
## 2.76058	0.831469	1.60223
## 4.01515	1.15557	1.53588
## 3.22636	1.37375	-1.65672
## 4.97539	0.760939	0.190012
## 4.38158	2.08244	-0.745703
## 4.67339	1.3979	-0.924286
## 3.37937	0.813165	-1.68794
## 4.8542	1.5958	-0.359974
## 4.14815	2.19582	0.864209
## 3.54144	0.223295	-1.67033
## 2.67516	-1.48112	0.532625
## 2.75609	-0.36393	1.56281
## 2.97253	0.177718	-1.65512
## 4.10254	-0.86271	-1.23671
## 4.88406	0.422935	-0.506982
## 3.01875	-1.26774	-1.09252
## 4.36976	-0.410129	1.27065
## 4.79486	-0.579758	0.609343
## 3.8726	-0.707977	1.44244
## 3.14802	-0.752984	-1.51469
## 4.35612	0.135931	-1.33895
## 3.23252	-1.56672	-0.654333
## 4.93256	-0.246193	-0.470066
## 4.67656	0.479203	0.901842
## 4.8206	-0.763065	-0.220075
## 4.53051	-0.970282	-0.69076
## 4.00273	-1.33313	-0.802011
## 3.30219	-1.04063	1.34419
## 2.64751	-0.872718	1.29595
## 4.3737	-1.32296	-0.155576
## 2.61384	-1.541	-0.142047
## 2.48304	-0.950605	-1.13579
## 3.98815	-0.0974654	1.55908
## 3.33228	-1.6999	0.0114309
## 2.64073	-1.39937	-0.688291
## 3.2404	0.0494107	1.69753
## 3.73107	-1.30226	1.01148
## 4.58433	-0.411858	-1.05618
## 3.77852	-0.367714	-1.59444
## 0.0660381	-1.67619	0.275722
## -0.518499	-1.13932	1.15026
## -0.235893	-0.285261	-1.65921
## -0.586135	-0.629111	1.46651
## 0.0181033	1.66149	0.359348
## 0.12942	-0.831931	-1.47687
## -0.472954	1.38793	-0.860204
## 0.21481	-1.33883	1.02537
## -0.460881	-1.00007	-1.29516
## -0.491812	0.843394	1.39169
## -0.54576	1.37757	0.833342
## -0.234615	0.801952	-1.48048
## -0.074239	0.0507306	1.69762
## -0.682068	-1.41573	-0.648448
## 0.0462573	-1.27611	-1.12223
## 0.488512	-0.166156	-1.6198
## 0.321571	-0.918684	1.39378
## 0.393719	0.468739	1.58596
## -0.529335	1.59791	-0.237701
## 0.443801	-0.371712	1.5984
## 0.307199	0.62341	-1.55145
## 0.186585	1.13753	1.24948
## -0.591614	-1.50216	0.532456
## -0.622001	0.198299	1.56965
## -0.909119	-0.181918	-1.50568
## -2.36154	0.2186	-0.642182
## -1.91762	-1.05198	0.593188
## -1.71575	0.531615	-1.23353
## -1.1743	-1.48946	0.215575
## -1.50771	0.884715	1.13111
## -2.1386	0.765019	0.661547
## -1.33334	1.37418	0.570295
## -1.14925	-1.23227	-0.870164
## -2.45608	0.397831	0.193536
## -1.73352	-0.0174981	1.33138
## -2.47167	-0.313266	0.211401
## -1.82559	1.25998	0.234329
## -2.19613	-0.430075	-0.831562
## -1.66563	-0.758896	-1.13327
## -2.21639	-0.345173	0.842217
## -1.3875	-0.466092	1.39234
## -1.71358	1.25851	-0.480506
## -1.08675	0.687493	-1.35547
## -2.16224	-0.973077	-0.0626776
## -1.63096	-0.890019	1.05436
## -1.51062	-0.138403	-1.42464
## -1.65999	-1.36268	-0.068253
## -0.877431	0.951618	1.18294
## -0.871612	-0.21869	1.4979
## """)

# POSCAR -- VASP format
class TestPOSCAR(unittest.TestCase, WriteMixin):
    fmt = "POSCAR"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
phenol
1.000 
0.0  0.0  0.0
0.0  0.0  0.0
0.0  0.0  0.0
C   O   
6   1   
Cartesian
     1.5846000000000000085     -0.0248999999999999985      0.0000000000000000000
     1.5703000000000000291      0.9755000000000000338      0.0000000000000000000
     2.4294999999999999929      1.4881999999999999673      0.0000000000000000000
     3.3031000000000001471      1.0003999999999999559      0.0000000000000000000
     3.3174999999999998934     -0.0000000000000000000      0.0000000000000000000
     0.0000000000000000000      0.0000000000000000000      0.0000000000000000000
    -1.0004999999999999449      0.0051000000000000004      0.0000000000000000000
""")

# POSFF -- MDFF format
class TestPOSFF(unittest.TestCase, WriteMixin):
    fmt = "POSFF"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
7
phenol
0.0  0.0  0.0
0.0  0.0  0.0
0.0  0.0  0.0
2
C   O   
6   1   
Cartesian
C        1.5846000000000000085     -0.0248999999999999985      0.0000000000000000000
C        1.5703000000000000291      0.9755000000000000338      0.0000000000000000000
C        2.4294999999999999929      1.4881999999999999673      0.0000000000000000000
C        3.3031000000000001471      1.0003999999999999559      0.0000000000000000000
C        3.3174999999999998934     -0.0000000000000000000      0.0000000000000000000
C        0.0000000000000000000      0.0000000000000000000      0.0000000000000000000
O       -1.0004999999999999449      0.0051000000000000004      0.0000000000000000000
""")

# pov -- POV-Ray input format [Write-only]
_pov_date = re.compile("//Date: [A-Za-z0-9 :]*")
def normalize_pov_date(content):
    return _pov_date.sub("//Date: Somewhere in time", content)

## class TestPOV(unittest.TestCase, WriteMixin):
## XXX Does not work on unpatched system
##     fmt = "pov"
##     maxDiff = None
##     def test_default(self):
##         self.assertWriters(self.fmt, """\
## //Povray v3 code generated by Open Babel
## //Author: Steffen Reith <streit@streit.cc>
## //Update (2010): Noel O'Boyle and Steven Wathen
## //Date: Tue Jan 15 15:28:49 CET 2019

## //Set some global parameters for display options
## #declare BAS = true;
## #declare TRANS = false;

## #include "colors.inc"

## // create a regular point light source
## light_source {
##   <3.60064,3.49204,-8>
##   color rgb <1,1,1>    // light's color
## }

## // set a color of the background (sky)
## background { color rgb <0.95 0.95 0.95> }

## // perspective (default) camera
## camera {
##   location  <1.60064,0.492043,-10>
##   look_at   <1.60064,0.492043,0>
##   right     x*image_width/image_height
## }

## //Include header for povray
## #include "babel_povray3.inc"

## //Use PovRay3.6
## #version 3.6;

## //Print name of molecule while rendering
## #render "\\b\\b phenol\\n\\n"

## //Coodinates of atoms 1 - 7
## #declare mol_0_pos_1 = <1.5846,-0.0249,0>;
## #declare mol_0_pos_2 = <1.5703,0.9755,0>;
## #declare mol_0_pos_3 = <2.4295,1.4882,0>;
## #declare mol_0_pos_4 = <3.3031,1.0004,0>;
## #declare mol_0_pos_5 = <3.3175,-0,0>;
## #declare mol_0_pos_6 = <0,0,0>;
## #declare mol_0_pos_7 = <-1.0005,0.0051,0>;

## //Povray-description of atoms 1 - 7
## #declare mol_0_atom1 = object {
## 	  Atom_C
## 	  translate mol_0_pos_1
## 	 }
## #declare mol_0_atom2 = object {
## 	  Atom_C
## 	  translate mol_0_pos_2
## 	 }
## #declare mol_0_atom3 = object {
## 	  Atom_C
## 	  translate mol_0_pos_3
## 	 }
## #declare mol_0_atom4 = object {
## 	  Atom_C
## 	  translate mol_0_pos_4
## 	 }
## #declare mol_0_atom5 = object {
## 	  Atom_C
## 	  translate mol_0_pos_5
## 	 }
## #declare mol_0_atom6 = object {
## 	  Atom_C
## 	  translate mol_0_pos_6
## 	 }
## #declare mol_0_atom7 = object {
## 	  Atom_O
## 	  translate mol_0_pos_7
## 	 }

## //Povray-description of bonds 1 - 7
## #if (BAS)
## #declare mol_0_bond0 = object {
## 	  bond_2
## 	  scale <1.5848,1.0000,1.0000>
## 	  rotate <0.0000,0.0000,0.900257>
## 	  rotate <0.0000,-180,0.0000>
## 	  translate mol_0_pos_1
## 	 }
## #declare mol_0_bond1 = object {
## 	  bond_1
## 	  scale <1.0005,1.0000,1.0000>
## 	  rotate <0.0000,0.0000,89.1811>
## 	  rotate <0.0000,-180,0.0000>
## 	  translate mol_0_pos_1
## 	 }
## #declare mol_0_bond2 = object {
## 	  bond_2
## 	  scale <1.00054,1.0000,1.0000>
## 	  rotate <0.0000,0.0000,30.8253>
## 	  translate mol_0_pos_2
## 	 }
## #declare mol_0_bond3 = object {
## 	  bond_1
## 	  scale <1.00056,1.0000,1.0000>
## 	  rotate <0.0000,0.0000,-29.1781>
## 	  translate mol_0_pos_3
## 	 }
## #declare mol_0_bond4 = object {
## 	  bond_2
## 	  scale <1.0005,1.0000,1.0000>
## 	  rotate <0.0000,0.0000,-89.1753>
## 	  translate mol_0_pos_4
## 	 }
## #declare mol_0_bond5 = object {
## 	  bond_1
## 	  scale <3.3175,1.0000,1.0000>
## 	  rotate <0.0000,-180,0.0000>
## 	  translate mol_0_pos_5
## 	 }
## #declare mol_0_bond6 = object {
## 	  bond_1
## 	  scale <1.00051,1.0000,1.0000>
## 	  rotate <0.0000,0.0000,0.29206>
## 	  rotate <0.0000,-180,0.0000>
## 	  translate mol_0_pos_6
## 	 }
## #end //(BAS-Bonds)

## #if (CST)
## #declare mol_0_bond0 = object {
## 	  union {
## 	   object {
## 	    bond_2
## 	    pigment{color Color_Car}
## 	    scale <0.792398,1.0000,1.0000>
## 	    rotate <0.0000,0.0000,0.900257>
## 	    rotate <0.0000,-180,0.0000>
## 	    translate mol_0_pos_1
## 	   }
## 	   object {
## 	    bond_2
## 	    pigment{color Color_Car}
## 	    scale <0.792398,1.0000,1.0000>
## 	    rotate <0.0000,0.0000,180.9>
## 	    rotate <0.0000,-180,0.0000>
## 	    translate mol_0_pos_6
## 	   }
## 	  }
## 	 }

## #declare mol_0_bond1 = object {
## 	  union {
## 	   object {
## 	    bond_1
## 	    pigment{color Color_Car}
## 	    scale <0.500251,1.0000,1.0000>
## 	    rotate <0.0000,0.0000,89.1811>
## 	    rotate <0.0000,-180,0.0000>
## 	    translate mol_0_pos_1
## 	   }
## 	   object {
## 	    bond_1
## 	    pigment{color Color_Car}
## 	    scale <0.500251,1.0000,1.0000>
## 	    rotate <0.0000,0.0000,269.181>
## 	    rotate <0.0000,-180,0.0000>
## 	    translate mol_0_pos_2
## 	   }
## 	  }
## 	 }

## #declare mol_0_bond2 = object {
## 	  union {
## 	   object {
## 	    bond_2
## 	    pigment{color Color_Car}
## 	    scale <0.500271,1.0000,1.0000>
## 	    rotate <0.0000,0.0000,30.8253>
## 	    translate mol_0_pos_2
## 	   }
## 	   object {
## 	    bond_2
## 	    pigment{color Color_Car}
## 	    scale <0.500271,1.0000,1.0000>
## 	    rotate <0.0000,0.0000,210.825>
## 	    translate mol_0_pos_3
## 	   }
## 	  }
## 	 }

## #declare mol_0_bond3 = object {
## 	  union {
## 	   object {
## 	    bond_1
## 	    pigment{color Color_Car}
## 	    scale <0.500281,1.0000,1.0000>
## 	    rotate <0.0000,0.0000,-29.1781>
## 	    translate mol_0_pos_3
## 	   }
## 	   object {
## 	    bond_1
## 	    pigment{color Color_Car}
## 	    scale <0.500281,1.0000,1.0000>
## 	    rotate <0.0000,0.0000,150.822>
## 	    translate mol_0_pos_4
## 	   }
## 	  }
## 	 }

## #declare mol_0_bond4 = object {
## 	  union {
## 	   object {
## 	    bond_2
## 	    pigment{color Color_Car}
## 	    scale <0.500252,1.0000,1.0000>
## 	    rotate <0.0000,0.0000,-89.1753>
## 	    translate mol_0_pos_4
## 	   }
## 	   object {
## 	    bond_2
## 	    pigment{color Color_Car}
## 	    scale <0.500252,1.0000,1.0000>
## 	    rotate <0.0000,0.0000,90.8247>
## 	    translate mol_0_pos_5
## 	   }
## 	  }
## 	 }

## #declare mol_0_bond5 = object {
## 	  union {
## 	   object {
## 	    bond_1
## 	    pigment{color Color_Car}
## 	    scale <1.65875,1.0000,1.0000>
## 	    rotate <0.0000,-180,0.0000>
## 	    translate mol_0_pos_5
## 	   }
## 	   object {
## 	    bond_1
## 	    pigment{color Color_Car}
## 	    scale <1.65875,1.0000,1.0000>
## 	    rotate <0.0000,0.0000,180>
## 	    rotate <0.0000,-180,0.0000>
## 	    translate mol_0_pos_6
## 	   }
## 	  }
## 	 }

## #declare mol_0_bond6 = object {
## 	  union {
## 	   object {
## 	    bond_1
## 	    pigment{color Color_Car}
## 	    scale <0.500256,1.0000,1.0000>
## 	    rotate <0.0000,0.0000,0.29206>
## 	    rotate <0.0000,-180,0.0000>
## 	    translate mol_0_pos_6
## 	   }
## 	   object {
## 	    bond_1
## 	    pigment{color Color_O3}
## 	    scale <0.500256,1.0000,1.0000>
## 	    rotate <0.0000,0.0000,180.292>
## 	    rotate <0.0000,-180,0.0000>
## 	    translate mol_0_pos_7
## 	   }
## 	  }
## 	 }

## #end // (CST-Bonds)


## //All atoms of molecule mol_0
## #ifdef (TRANS)
## #declare mol_0_atoms = merge {
## #else
## #declare mol_0_atoms = union {
## #end //(End of TRANS)
## 	  object{mol_0_atom1}
## 	  object{mol_0_atom2}
## 	  object{mol_0_atom3}
## 	  object{mol_0_atom4}
## 	  object{mol_0_atom5}
## 	  object{mol_0_atom6}
## 	  object{mol_0_atom7}
## 	 }

## //Bonds only needed for ball and sticks or capped sticks models
## #if (BAS | CST)
## #declare mol_0_bonds = union {
## 	  object{mol_0_bond0}
## 	  object{mol_0_bond1}
## 	  object{mol_0_bond2}
## 	  object{mol_0_bond3}
## 	  object{mol_0_bond4}
## 	  object{mol_0_bond5}
## 	  object{mol_0_bond6}
## 	 }
## #end


## //Definition of molecule mol_0
## #if (SPF)
## #declare mol_0 = object{
## 	  mol_0_atoms
## #else
## #declare mol_0 = union {
## 	  object{mol_0_atoms}
## #if (BAS | CST)//(Not really needed at moment!)
## #if (TRANS)
## 	  difference {
## 	   object{mol_0_bonds}
## 	   object{mol_0_atoms}
## 	  }
## #else
## 	  object{mol_0_bonds}
## #end //(End of TRANS)
## #end //(End of (BAS|CST))
## #end //(End of SPF)
## //	  bounded_by {
## //	   box {
## //	    <-4.0005,-3.0249,-3>
## //	    <6.3175,4.4882,3>
## 	 }

## //Center of molecule mol_0 (bounding box)
## #declare mol_0_center = <-1.1585,-0.73165,-0>;

## mol_0
## """, normalize=normalize_pov_date)

# pqr -- PQR format
class TestPQR(unittest.TestCase, WriteMixin):
    fmt = "pqr"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
COMPND    phenol 
AUTHOR    GENERATED BY OPEN BABEL %(VERSION)s
HETATM    1  C   UNL     1       1.585  -0.025   0.000  0.04202806   1.700  C  
HETATM    2  C   UNL     1       1.570   0.976   0.000  0.00328151   1.700  C  
HETATM    3  C   UNL     1       2.429   1.488   0.000  0.00020584   1.700  C  
HETATM    4  C   UNL     1       3.303   1.000   0.000  0.00328151   1.700  C  
HETATM    5  C   UNL     1       3.317   0.000   0.000  0.04202806   1.700  C  
HETATM    6  C   UNL     1       0.000   0.000   0.000  0.19574524   1.700  C  
HETATM    7  O   UNL     1      -1.000   0.005   0.000 -0.28657022   1.520  O  
CONECT    1    6    2                                                 
CONECT    2    1    3                                                 
CONECT    3    2    4                                                 
CONECT    4    3    5                                                 
CONECT    5    4    6                                                 
CONECT    6    1    5    7                                            
CONECT    7    6                                                      
MASTER        0    0    0    0    0    0    0    0    7    0    7    0
END
""" % dict(VERSION=VERSION))

# pqs -- Parallel Quantum Solutions format
class TestPQS(unittest.TestCase, WriteMixin):
    fmt = "pqs"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
TEXT=phenol
GEOM=PQS
C             1.584600    -0.024900     0.000000
C             1.570300     0.975500     0.000000
C             2.429500     1.488200     0.000000
C             3.303100     1.000400     0.000000
C             3.317500    -0.000000     0.000000
C             0.000000     0.000000     0.000000
O            -1.000500     0.005100     0.000000
""")

# qcin -- Q-Chem input format [Write-only]
class TestQCIN(unittest.TestCase, WriteMixin):
    fmt = "qcin"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
$comment
phenol
$end

$molecule
0 1
6 1.5846 -0.0249 0
6 1.5703 0.9755 0
6 2.4295 1.4882 0
6 3.3031 1.0004 0
6 3.3175 -0 0
6 0 0 0
8 -1.0005 0.0051 0
$end

$rem

$end
""")

# report -- Open Babel report format [Write-only]
class TestREPORT(unittest.TestCase, WriteMixin):
    fmt = "report"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
FILENAME: phenol
FORMULA: C6H6O
MASS: 94.1112
EXACT MASS: 94.0418648
INTERATOMIC DISTANCES

              C   1      C   2      C   3      C   4      C   5      C   6
              ------------------------------------------------------------------
   C   1    0.0000 
   C   2    1.0005     0.0000 
   C   3    1.7330     1.0005     0.0000 
   C   4    2.0011     1.7330     1.0006     0.0000 
   C   5    1.7331     2.0011     1.7330     1.0005     0.0000 
   C   6    1.5848     1.8486     2.8491     3.4513     3.3175     0.0000 
   O   7    2.5853     2.7479     3.7369     4.4172     4.3180     1.0005 

              O   7
              -----------
   O   7    0.0000 



ATOMIC CHARGES
   C   1    0.0420280594
   C   2    0.0032815119
   C   3    0.0002058430
   C   4    0.0032815119
   C   5    0.0420280594
   C   6    0.1957452352
   O   7   -0.2865702207


BOND ANGLES
   2    1    6  Car  Car  Car     88.281
   1    2    3  Car  Car  Car    120.006
   2    3    4  Car  Car  Car    119.997
   3    4    5  Car  Car  Car    120.003
   4    5    6  Car  Car  Car     89.175
   1    6    5  Car  Car  Car      0.900
   1    6    7  Car  Car   O3    179.392
   5    6    7  Car  Car   O3    179.708


TORSION ANGLES
   2    1    6    5     -0.000
   2    1    6    7   -180.000
   6    1    2    3    180.000
   1    2    3    4      0.000
   2    3    4    5      0.000
   3    4    5    6     -0.000
   4    5    6    1   -180.000
   4    5    6    7     -0.000


""")

# rinchi -- RInChI [Write-only]
class TestRINCHI(unittest.TestCase, WriteMixin):
    fmt = "rinchi"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
RInChI=1.00.1S/Au<>Pb/d-
""", mol=_alchemy_mol)

## # rsmi -- Reaction SMILES format
## XXX I don't know why this fails
## class TestRSMI(unittest.TestCase, WriteMixin):
##     fmt = "rsmi"
##     def test_default(self):
##         self.assertWriters(self.fmt, """\
## """, mol="[Pb]>>[Au]")

# rxn -- MDL RXN format
class TestRXN(unittest.TestCase, WriteMixin):
    fmt = "rxn"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, ALCHEMY_RXN,
                           mol=_alchemy_mol, normalize=normalize_sd_timestamp)

# sd -- MDL MOL format
class TestSD(unittest.TestCase, WriteMixin):
    fmt = "sd"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
phenol
 OpenBabel01151915132D

  7  7  0  0  0  0  0  0  0  0999 V2000
    1.5846   -0.0249    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5703    0.9755    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4295    1.4882    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3031    1.0004    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3175   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0005    0.0051    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  6  2  0  0  0  0
  1  2  1  0  0  0  0
  2  3  2  0  0  0  0
  3  4  1  0  0  0  0
  4  5  2  0  0  0  0
  5  6  1  0  0  0  0
  6  7  1  0  0  0  0
M  END
$$$$
""", normalize=normalize_sd_timestamp)

# sdf -- MDL MOL format
class TestSDF(unittest.TestCase, WriteMixin):
    fmt = "sdf"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
phenol
 OpenBabel01151915132D

  7  7  0  0  0  0  0  0  0  0999 V2000
    1.5846   -0.0249    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5703    0.9755    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4295    1.4882    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3031    1.0004    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3175   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0005    0.0051    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  6  2  0  0  0  0
  1  2  1  0  0  0  0
  2  3  2  0  0  0  0
  3  4  1  0  0  0  0
  4  5  2  0  0  0  0
  5  6  1  0  0  0  0
  6  7  1  0  0  0  0
M  END
$$$$
""", normalize=normalize_sd_timestamp)


class _BaseSmiles(object):
    def test_default(self):
        self.assertWriters("smi", "c1ccccc1O\tphenol\n")
        
    def test_kekule(self):
        self.assertWriters("smi", "C1C=CC=CC=1O\tphenol\n", options=["k"])
        
    ## def test_explicit_hydrogens(self):
    ##     self.assertWriters("smi", "C1C=CC=CC=1O\n", options=["h"])
        
    def test_no_molecule_name(self):
        self.assertWriters("smi", "C1C=CC=CC=1O\n", options=["k", "n"])
    def test_molecule_name_only(self):
        self.assertWriters("smi", "phenol\n", options=["t"])
        
    def test_append_coordinates(self):
        self.assertWriters("smi", "c1ccccc1O\tphenol\t1.5846,-0.0249,1.5703,0.9755,2.4295,1.4882,3.3031,1.0004,3.3175,-0.0000,0.0000,0.0000,-1.0005,0.0051\n", options=["x"])
        
    def test_reuse_ring_closures(self):
        self.assertWriters("smi", "c1ccccc1c1ccccc1\tblah\n", mol="c1ccccc1c1ccccc1 blah")
    def test_do_not_reuse_ring_closures(self):
        self.assertWriters("smi", "c1ccccc1c2ccccc2\tblah\n", mol="c1ccccc1c1ccccc1 blah",
                               options=["R"])
    
    def test_fragment_smiles(self):
        self.assertWriters("smi", "P=N\n", mol="P=N-C=O blah2",
                               options={"n": None, "F": "1 2"})
        
    def test_atom_priority_order(self):
        self.assertWriters("smi", "Oc1ccccc1\tphenol\n",
                               options={"o": "7-6-5-4-3-2-1"})
    def test_first_atom(self):
        self.assertWriters("smi", "c1cc(ccc1)O\tphenol\n",
                               options={"f": "2"})
    def test_last_atom(self):
        self.assertWriters("smi", "c1c(cccc1O)\tphenol\n",
                               options={"l": "2"})
        
    def test_disable_isomeric(self):
        self.assertWriters("smi", "C[C@]12CCC(=O)[C@@]1(C)CCCC2O\tXYZ\n",
                               options=[], mol="C[C@]12CCC(=O)[C@@]1(C)CCCC2O XYZ")
        self.assertWriters("smi", "CC12CCC(=O)C1(C)CCCC2O\tXYZ\n",
                               options=["i"], mol="C[C@]12CCC(=O)[C@@]1(C)CCCC2O XYZ")
        

# smi -- SMILES format
class TestSMI(unittest.TestCase, WriteMixin, _BaseSmiles):
    fmt = "smi"
    maxDiff = None

# smiles -- SMILES format
class TestSMILES(unittest.TestCase, WriteMixin,  _BaseSmiles):
    fmt = "smiles"
    maxDiff = None

## # stl -- STL 3D-printing format [Write-only]
## XXX the output is far too extensive to test here
## class TestSTL(unittest.TestCase, WriteMixin):
##     fmt = "stl"
##     maxDiff = None
##     def test_default(self):
##         self.assertBinaryWriteFile(self.fmt, """\
## """)

# svg -- SVG 2D depiction [Write-only]
class TestSVG(unittest.TestCase, WriteMixin):
    fmt = "svg"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
<?xml version="1.0"?>
<svg version="1.1" id="topsvg"
xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink"
xmlns:cml="http://www.xml-cml.org/schema" x="0" y="0" width="200px" height="200px" viewBox="0 0 100 100">
<title>phenol - Open Babel Depiction</title>
<rect x="0" y="0" width="100" height="100" fill="white"/>
<g transform="translate(0,0)">
<svg width="100" height="100" x="0" y="0" viewBox="0 0 202.065 122.773"
font-family="sans-serif" stroke="rgb(0,0,0)" stroke-width="2"  stroke-linecap="round">
<line x1="68.3" y1="82.1" x2="53.0" y2="82.0" opacity="1.0" stroke="rgb(0,0,0)"  stroke-width="2.0"/>
<line x1="137.0" y1="40.0" x2="161.7" y2="53.8" opacity="1.0" stroke="rgb(0,0,0)"  stroke-width="2.0"/>
<line x1="161.7" y1="53.8" x2="162.1" y2="82.1" opacity="1.0" stroke="rgb(0,0,0)"  stroke-width="2.0"/>
<line x1="154.5" y1="59.9" x2="154.8" y2="76.2" opacity="1.0" stroke="rgb(0,0,0)"  stroke-width="2.0"/>
<line x1="162.1" y1="82.1" x2="68.3" y2="82.1" opacity="1.0" stroke="rgb(0,0,0)"  stroke-width="2.0"/>
<line x1="68.3" y1="82.1" x2="113.1" y2="82.8" opacity="1.0" stroke="rgb(0,0,0)"  stroke-width="2.0"/>
<line x1="74.4" y1="75.0" x2="107.2" y2="75.5" opacity="1.0" stroke="rgb(0,0,0)"  stroke-width="2.0"/>
<line x1="113.1" y1="82.8" x2="112.7" y2="54.5" opacity="1.0" stroke="rgb(0,0,0)"  stroke-width="2.0"/>
<line x1="112.7" y1="54.5" x2="137.0" y2="40.0" opacity="1.0" stroke="rgb(0,0,0)"  stroke-width="2.0"/>
<line x1="121.5" y1="57.6" x2="135.5" y2="49.3" opacity="1.0" stroke="rgb(0,0,0)"  stroke-width="2.0"/>
<text x="22.000000" y="89.925427" fill="rgb(255,12,12)" stroke-width="0" font-weight="bold" font-size="16" >HO</text>
</svg>
</g>
<text font-size="18.000000" fill ="black" font-family="sans-serif"
x="10.000000" y="20.000000" >phenol</text>
</svg>
""")

# sy2 -- Sybyl Mol2 format
class TestSY2(unittest.TestCase, WriteMixin):
    fmt = "sy2"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
@<TRIPOS>MOLECULE
phenol
 7 7 0 0 0
SMALL
GASTEIGER

@<TRIPOS>ATOM
      1 C           1.5846   -0.0249    0.0000 C.ar    1  UNL1        0.0420
      2 C           1.5703    0.9755    0.0000 C.ar    1  UNL1        0.0033
      3 C           2.4295    1.4882    0.0000 C.ar    1  UNL1        0.0002
      4 C           3.3031    1.0004    0.0000 C.ar    1  UNL1        0.0033
      5 C           3.3175   -0.0000    0.0000 C.ar    1  UNL1        0.0420
      6 C           0.0000    0.0000    0.0000 C.ar    1  UNL1        0.1957
      7 O          -1.0005    0.0051    0.0000 O.3     1  UNL1       -0.2866
@<TRIPOS>BOND
     1     1     6   ar
     2     1     2   ar
     3     2     3   ar
     4     3     4   ar
     5     4     5   ar
     6     5     6   ar
     7     6     7    1
""")

## # tdd -- Thermo format
## XXX need thermo data
## class TestTDD(unittest.TestCase, WriteMixin):
##     fmt = "tdd"
##     maxDiff = None
##     def test_default(self):
##         self.assertWriters(self.fmt, """\
## """)

## # text -- Read and write raw text
## XXX Not valid output format?
## class TestTEXT(unittest.TestCase, WriteMixin):
##     fmt = "text"
##     maxDiff = None
##     def test_default(self):
##         self.assertWriters(self.fmt, """\
## """)

## # therm -- Thermo format
## XXX need thermo data
## class TestTHERM(unittest.TestCase, WriteMixin):
##     fmt = "therm"
##     maxDiff = None
##     def test_default(self):
##         self.assertWriters(self.fmt, """\
## """)

# tmol -- TurboMole Coordinate format
class TestTMOL(unittest.TestCase, WriteMixin):
    fmt = "tmol"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
$coord
    2.99446001766484     -0.04705418051234      0.00000000000000      c
    2.96743693407743      1.84342783493125      0.00000000000000      c
    4.59108962067193      2.81229041921546      0.00000000000000      c
    6.24195436346633      1.89048201544359      0.00000000000000      c
    6.26916641966623     -0.00000000000000      0.00000000000000      c
    0.00000000000000      0.00000000000000      0.00000000000000      c
   -1.89067098805609      0.00963760323747      0.00000000000000      o
$end
""")

# txt -- Title format
class TestTXT(unittest.TestCase, WriteMixin):
    fmt = "txt"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
phenol
""")

# txyz -- Tinker XYZ format
class TestTXYZ(unittest.TestCase, WriteMixin):
    fmt = "txyz"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
     7 phenol                 MM2 parameters
     1  C      1.584600   -0.024900    0.000000     2     6     2
     2  C      1.570300    0.975500    0.000000     2     1     3
     3  C      2.429500    1.488200    0.000000     2     2     4
     4  C      3.303100    1.000400    0.000000     2     3     5
     5  C      3.317500   -0.000000    0.000000     2     4     6
     6  C      0.000000    0.000000    0.000000     2     1     5     7
     7  O     -1.000500    0.005100    0.000000     6     6
""")

# unixyz -- UniChem XYZ format
class TestUNIXYZ(unittest.TestCase, WriteMixin):
    fmt = "unixyz"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
phenol
7
  6        1.58460       -0.02490        0.00000
  6        1.57030        0.97550        0.00000
  6        2.42950        1.48820        0.00000
  6        3.30310        1.00040        0.00000
  6        3.31750       -0.00000        0.00000
  6        0.00000        0.00000        0.00000
  8       -1.00050        0.00510        0.00000
""")

# VASP -- VASP format
class TestVASP(unittest.TestCase, WriteMixin):
    fmt = "VASP"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
phenol
1.000 
0.0  0.0  0.0
0.0  0.0  0.0
0.0  0.0  0.0
C   O   
6   1   
Cartesian
     1.5846000000000000085     -0.0248999999999999985      0.0000000000000000000
     1.5703000000000000291      0.9755000000000000338      0.0000000000000000000
     2.4294999999999999929      1.4881999999999999673      0.0000000000000000000
     3.3031000000000001471      1.0003999999999999559      0.0000000000000000000
     3.3174999999999998934     -0.0000000000000000000      0.0000000000000000000
     0.0000000000000000000      0.0000000000000000000      0.0000000000000000000
    -1.0004999999999999449      0.0051000000000000004      0.0000000000000000000
""")

# vmol -- ViewMol format
class TestVMOL(unittest.TestCase, WriteMixin):
    fmt = "vmol"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
$title
phenol
$coord 1.0
      1.58460000000000     -0.02490000000000      0.00000000000000 C
      1.57030000000000      0.97550000000000      0.00000000000000 C
      2.42950000000000      1.48820000000000      0.00000000000000 C
      3.30310000000000      1.00040000000000      0.00000000000000 C
      3.31750000000000     -0.00000000000000      0.00000000000000 C
      0.00000000000000      0.00000000000000      0.00000000000000 C
     -1.00050000000000      0.00510000000000      0.00000000000000 O
$end
""")

# xed -- XED format [Write-only]
class TestXED(unittest.TestCase, WriteMixin):
    fmt = "xed"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
     0.000         7         7
File conversion by Open Babel
       1       6       1       2       2       3       3       4       4       5
       5       6       6       7
     6       1.584600      -0.024900       0.000000     3      0.0000
     6       1.570300       0.975500       0.000000     3      0.0000
     6       2.429500       1.488200       0.000000     3      0.0000
     6       3.303100       1.000400       0.000000     3      0.0000
     6       3.317500      -0.000000       0.000000     3      0.0000
     6       0.000000       0.000000       0.000000     3      0.0000
     8      -1.000500       0.005100       0.000000    10      0.0000
    1         0.0000    0         0.0000
""")

# xyz -- XYZ cartesian coordinates format
class TestXYZ(unittest.TestCase, WriteMixin):
    fmt = "xyz"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
7
phenol
C          1.58460       -0.02490        0.00000
C          1.57030        0.97550        0.00000
C          2.42950        1.48820        0.00000
C          3.30310        1.00040        0.00000
C          3.31750       -0.00000        0.00000
C          0.00000        0.00000        0.00000
O         -1.00050        0.00510        0.00000
""")

# yob -- YASARA.org YOB format
class TestYOB(unittest.TestCase, WriteMixin):
    fmt = "yob"
    maxDiff = None
    def test_default(self):
        self.assertBinaryWriters(self.fmt, b'YMOB\x90\x00\x00\x00\x06\x00\x00\x00\x88\x00\x00\x00\x00\x00\x00\x00\x00\x00\xf0?\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\xf0?\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\xf0?\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\xf0?\xff\xff\xff\x7f\x08\x00\x00\x00$\x01\x00\x00\x07\x00\x00\x00\x01\x00\x00\x00\x06\x00\x00\x00\x02\x04\x06@\x04\x95\xfd\xffF\xf6\xff\xff\x00\x00\x00\x00\x05\x00\x00\x02\x01\x00\x00\x01\x03\x00\x00\x00C   UNK    1\x02\x04\x06@\x9a\x9a\xfd\xff\x0e}\x01\x00\x00\x00\x00\x00\x00\x00\x00\x01\x02\x00\x00\x02\x03\x00\x00\x00C   UNK    1\x02\x04\x06@\xfaJ\xfc\xffTE\x02\x00\x00\x00\x00\x00\x01\x00\x00\x02\x03\x00\x00\x01\x03\x00\x00\x00C   UNK    1\x02\x04\x06@\xba\xf5\xfa\xff\xc8\x86\x01\x00\x00\x00\x00\x00\x02\x00\x00\x01\x04\x00\x00\x02\x03\x00\x00\x00C   UNK    1\x02\x04\x06@\x1a\xf0\xfa\xff\x00\x00\x00\x00\x00\x00\x00\x00\x03\x00\x00\x02\x05\x00\x00\x01\x03\x00\x00\x00C   UNK    1\x03\x04\x06@\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x02\x04\x00\x00\x01\x06\x00\x00\x01\x03\x00\x00\x00C   UNK    1\x01\x04\x08@\xd2\x86\x01\x00\xfe\x01\x00\x00\x00\x00\x00\x00\x05\x00\x00\x01\x03\x00\x00\x00O   UNK    1')

# zin -- ZINDO input format [Write-only]
class TestZIN(unittest.TestCase, WriteMixin):
    fmt = "zin"
    maxDiff = None
    def test_default(self):
        self.assertWriters(self.fmt, """\
 $TITLEI

   phenol

 $END

 $CONTRL

 SCFTYP       ROHF   RUNTYP       CI   ENTTYP     COORD
 UNITS        ANGS   INTTYP        1   IAPX           3

 NOP = 1 
 NDT = 1 
 FOP(1) =  29  1.000000
 NAT             7   NEL          30   MULT           1
 IPRINT         -1   ITMAX       100

! ***** BASIS SET AND C. I. SIZE INFORMATION ***** 

 DYNAL(1) =     0    0    7    0    0 1200   40

 INTFA(1) =   1.000000 1.267000  0.680000  1.000000  1.000000 

! ***** OUTPUT FILE NAME ***** 

   ONAME =  zindo 

 $END

 $DATAIN 

  1.584600 -0.024900  0.000000    6
  1.570300  0.975500  0.000000    6
  2.429500  1.488200  0.000000    6
  3.303100  1.000400  0.000000    6
  3.317500 -0.000000  0.000000    6
  0.000000  0.000000  0.000000    6
 -1.000500  0.005100  0.000000    8



 $END 

 $CIINPU

! ***** C. I. SPECIFICATION *****

    2    1   25    1    0    0    0    1   10    1   10
  -60000.0 0.0000000

    1   15   15   16
   21    7   16   16   26

 $END 
""")

        

if __name__ == "__main__":
    unittest.main()

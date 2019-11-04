from __future__ import print_function
from openbabel import openbabel as ob

# Initialize the OBConversion object
conv = ob.OBConversion()
if not conv.SetInFormat('smi'):
  print('could not find smiles format')

# Read the smiles string
mol = ob.OBMol()
if not conv.ReadString(mol, 'CCCC'):
  print('could not read the smiles string')

# ... Use OBMol object ...

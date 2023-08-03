from __future__ import print_function
from openbabel import openbabel

mol = openbabel.OBMol()
print('Should print 0 (atoms)')
print(mol.NumAtoms())

a = mol.NewAtom()
b = mol.NewAtom()
mol.AddBond(1, 2, 1)
print('Should print 2 (atoms)')
print(mol.NumAtoms())
print('Should print 1 (bond)')
print(mol.NumBonds())

obConversion = openbabel.OBConversion()
obConversion.SetInAndOutFormats("smi", "mdl")

mol.Clear()
obConversion.ReadString(mol, "C1=CC=CS1")

print('Should print 5 (atoms)')
print(mol.NumAtoms())

mol.AddHydrogens()
print('Should print 9 (atoms) after adding hydrogens')
print(mol.NumAtoms())

outMDL = obConversion.WriteString(mol)

obConversion.WriteFile(mol, 'temp.mdl')

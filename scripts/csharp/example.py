import clr
clr.AddReference("OBDotNet.dll")
import OpenBabel as ob
conv = ob.OBConversion()
conv.SetInFormat("smi")
mol = ob.OBMol()
conv.ReadString(mol, "CC=O")
print mol.GetMolWt()

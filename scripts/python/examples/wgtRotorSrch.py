######################################################################
#
# wgtRotorSrch.py: weighted rotor search (conformer search)
#
######################################################################

from __future__ import print_function
from openbabel import openbabel
import sys

# Make sure we have a filename
try:
  filename = sys.argv[1]
except:
  print("Usage: python energy.py filename")
  sys.exit(1)

# Read the file.
mol = openbabel.OBMol()
conv = openbabel.OBConversion()
format = conv.FormatFromExt(filename)
conv.SetInAndOutFormats(format, format)
conv.ReadFile(mol, filename)

# Find the MMFF94 force field.
ff = openbabel.OBForceField.FindForceField("MMFF94")
if ff == 0:
  print("Could not find forcefield")

# Set the log level to low since we only want to print out the conformer search
# steps and not all individual interactions for each call to Energy()
ff.SetLogLevel(openbabel.OBFF_LOGLVL_LOW)
# python specific, python doesn't have std::ostream so the SetLogFile()
# function is replaced by SetLogToStdOut and SetLogToStdErr in the SWIG
# interface file
ff.SetLogToStdErr()

# Setup the molecule. This assigns atoms types, charges and parameters
if ff.Setup(mol) == 0:
  print("Could not setup forcefield")

# Weighted rotor search: generate 25 conformers, optimize each conformer for
# 500 steps.
ff.WeightedRotorSearch(25, 500)

# Get all the coordinates back from the force field. The best conformer is also
# set in mol
ff.GetConformers(mol)
# Write the best conformer back to the file
conv.WriteFile(mol, filename)
# Other conformers can also be written by calling:
# mol.SetConformer(0)
# conv.Write(mol)
# mol.SetConformer(1)
# conv.Write(mol)
# ...

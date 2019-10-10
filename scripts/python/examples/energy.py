######################################################################
#
# energy.py: calculate the energy of a molecule
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

# Set the log level to high since we want to print out individual
# interactions.
ff.SetLogLevel(openbabel.OBFF_LOGLVL_HIGH)
# python specific, python doesn't have std::ostream so the SetLogFile()
# function is replaced by SetLogToStdOut and SetLogToStdErr in the SWIG
# interface file
ff.SetLogToStdErr()

# Setup the molecule. This assigns atoms types, charges and parameters
if ff.Setup(mol) == 0:
  print("Could not setup forcefield")

# Calculate the energy
ff.Energy()

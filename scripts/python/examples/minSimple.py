######################################################################
#
# minSimple.py: minimize a molecule
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

# Set the log level to low since we only want to print out the minimization
# steps and not all individual interactions for each call to Energy()
ff.SetLogLevel(openbabel.OBFF_LOGLVL_LOW)
# python specific, python doesn't have std::ostream so the SetLogFile()
# function is replaced by SetLogToStdOut and SetLogToStdErr in the SWIG
# interface file
ff.SetLogToStdErr()

# Setup the molecule. This assigns atoms types, charges and parameters
if ff.Setup(mol) == 0:
  print("Could not setup forcefield")

# Minimize using steepest descent for 2000 steps
ff.SteepestDescent(2000)

# Get the coordinates back from the force field.
ff.GetCoordinates(mol)
# Write the minimized molecule back to the file
conv.WriteFile(mol, filename)

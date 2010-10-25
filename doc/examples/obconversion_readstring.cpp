#include <openbabel/obconversion.h>
#include <openbabel/mol.h>

using namespace OpenBabel;

int main()
{
  // Create an OBConversion object.
  OBConversion conv;
  // Set the input format.
  if (!conv.SetInFormat("smi")) {
    // Handle error.
    return 1;
  }

  // Create the OBMol object.
  OBMol mol;

  // Read the smiles string.
  if (conv.ReadString(&mol, "CCCC")) {
    // Handle error.
    return 1;
  }

  // ...Use OBMol object...
}

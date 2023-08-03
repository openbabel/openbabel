#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/shared_ptr.h>
#include <openbabel/forcefield.h>

#include <iostream>

using namespace OpenBabel;

// Helper function to read molecule from file
shared_ptr<OBMol> GetMol(const std::string &filename)
{
  // Create the OBMol object.
  shared_ptr<OBMol> mol(new OBMol);

  // Create the OBConversion object.
  OBConversion conv;
  OBFormat *format = conv.FormatFromExt(filename.c_str());
  if (!format || !conv.SetInFormat(format)) {
    std::cout << "Could not find input format for file " << filename << std::endl;
    return mol;
  }

  // Open the file.
  std::ifstream ifs(filename.c_str());
  if (!ifs) {
    std::cout << "Could not open " << filename << " for reading." << std::endl;
    return mol;
  }
  // Read the molecule.
  if (!conv.Read(mol.get(), &ifs)) {
    std::cout << "Could not read molecule from file " << filename << std::endl;
    return mol;
  }

  return mol;
}

int main(int argc, char **argv)
{
  if (argc < 2) {
    std::cout << "Usage: " << argv[0] << " <filename>" << std::endl;
    return 1;
  }

  // Read the file.
  shared_ptr<OBMol> mol = GetMol(argv[1]);

  // Get the forcefield.
  OBForceField *ff = OBForceField::FindType("MMFF94");
  if (!ff) {
    std::cout << "Could not find forcefield." << std::endl;
    return 1;
  }

  // Setup the forcefield.
  if (!ff->Setup(mol)) {
    std::cout << "Could not setup forcefield." << std::endl;
    return 1;
  }

  // Print the enegy and unit.
  std::cout << ff->Energy() << " " << ff->GetUnit() << std::endl;

  return 0;
}

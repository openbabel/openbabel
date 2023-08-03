#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/shared_ptr.h>
#include <openbabel/conformersearch.h>

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

  // Read the file
  shared_ptr<OBMol> mol = GetMol(argv[1]);

  // Create the OBConformerSearch object
  OBConformerSearch cs;
  
  // Setup 
  std::cout << "Setting up conformer searching..." << std::endl
            << "   conformers:  30" << std::endl
            << "   children:    5" << std::endl
            << "   mutability:  5" << std::endl
            << "   convergence: 25" << std::endl;
  cs.Setup(*mol.get(),
           30, // numConformers
           5, // numChildren
           5, // mutability
           25); // convergence

  // Perform searching
  cs.Search();

  // Print the rotor keys
  RotorKeys keys = cs.GetRotorKeys();
  for (RotorKeys::iterator key = keys.begin(); key != keys.end(); ++key) {
    for (unsigned int i = 1; i < key->size(); ++i)
      std::cout << key->at(i) << " ";
    std::cout << std::endl;
  }

  // Get the conformers
  cs.GetConformers(*mol.get());

  std::cout << mol->NumConformers() << std::endl;
  OBConversion conv;
  conv.SetOutFormat("sdf");
  for (unsigned int c = 0; c < mol->NumConformers(); ++c) {
    mol->SetConformer(c);
    conv.Write(mol.get(), &std::cerr);
  }

}

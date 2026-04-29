#include <openbabel/babelconfig.h>
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/depict/svgpainter.h>
#include <openbabel/depict/depict.h>

#include <iostream>
#include <sstream>

using namespace OpenBabel;

extern "C" int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
  if (size < 1) return 0;

  // Use the first byte to decide some options
  uint8_t options = data[0];
  const char* fuzzed_data = reinterpret_cast<const char*>(data + 1);
  size_t fuzzed_size = size - 1;

  std::string input(fuzzed_data, fuzzed_size);
  
  OBConversion conv;
  // We can try different input formats, but SMILES is most versatile for random strings
  if (!conv.SetInFormat("SMI")) {
    return 0;
  }

  OBMol mol;
  if (conv.ReadString(&mol, input)) {
    // Successfully read a molecule, now try to depict it
    
    std::stringstream oss;
    std::set<ColorGradient> gradients;
    SVGPainter painter(oss, &gradients);
    OBDepict depictor(&painter);
    
    if (options & 0x02) depictor.SetAliasMode();
    if (options & 0x04) depictor.SetBondLength(10.0);
    
    depictor.DrawMolecule(&mol);

    // Set some options based on the first byte
    if (options & 0x01) depictor.AddAtomLabels(OBDepict::AtomSymmetryClass);
    // Destructor of SVGPainter finishes the SVG
  }

  return 0;
}

#include <openbabel/babelconfig.h>
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
 
#include <openbabel/depict/svgpainter.h>
#include <openbabel/depict/depict.h>
 
#include <iostream>
using namespace std;
 
using namespace OpenBabel;
 
int main(int argc, char **argv)
{
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " filename\n"
      "Writes an SVG file containing a 2D depiction for each molecule." << std::endl;
    return -1;
  }
 
  // read a molecule
  OBConversion conv;
  OBFormat *format = conv.FormatFromExt(argv[1]);
  conv.SetInFormat(format);
 
  std::string filename(argv[1]);
  std::ifstream ifs;
  ifs.open(filename.c_str());
  if (!ifs) {
    std::cerr << "Could not open " << filename << std::endl;
    return -1;
  }
 
  string::size_type pos = filename.find_last_of(".");
  if (pos != string::npos)
    filename = filename.substr(0, pos);

  unsigned int count = 1;
  OBMol mol;
  while (conv.Read(&mol, &ifs)) {
    if(!mol.Has2D(true)) {
      cerr << "Molecule " << count << ", " << mol.GetTitle()
           << ", does not have 2D coordinates" << endl;
      continue;
    }
    std::stringstream ss;
    ss << filename << count << ".svg";
    ofstream ofs(ss.str().c_str());
    if(ofs)
    {
      SVGPainter painter(ofs);     
      OBDepict depictor(&painter);
      depictor.DrawMolecule(&mol);   // <svg...> written here
      //depictor.AddAtomLabels(OBDepict::AtomSymmetryClass);
    }                                // </svg> written here 
    //Note that the svg is not complete until the destructor of the SVGPainter has been called
    ofs.close();
    count++;
  }
 
  return 0;
}
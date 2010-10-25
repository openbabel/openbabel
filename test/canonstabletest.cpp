#include "obtest.h"
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>

/*
#include <openbabel/graphsym.h>
#include <openbabel/stereo/tetrahedral.h>
#include <openbabel/canon.h>
*/

#include <iostream>
#include <vector>
#include <algorithm>

std::string GetFilename(const std::string &filename)
{
  std::string path = TESTDATADIR + filename;
  return path;
}


using std::cout;
using std::endl;
using namespace OpenBabel;

bool testCanSmiles(const std::string &smiles, const std::string &stable_cansmiles)
{
  cout << " Testing: " << smiles << endl;
  // read a smiles string
  OBMol mol;
  OBConversion canConv, smiConv;
  OB_REQUIRE( canConv.SetInFormat("smi") );
  OB_REQUIRE( canConv.SetOutFormat("can") );
  OB_REQUIRE( smiConv.SetOutFormat("smi") );
  // read a smiles string
  OB_REQUIRE( canConv.ReadString(&mol, smiles) );


  // get can smiles
  std::string cansmiles = canConv.WriteString(&mol, true);
  OB_COMPARE( cansmiles, stable_cansmiles );
  // comapare with ref
  if (cansmiles != stable_cansmiles) {
    cout << " " << cansmiles << endl;
    cout << " " << stable_cansmiles << endl;
    return false;
  }

  return true;
}


int main(int argc, char **argv)
{
  // Define location of file formats for testing
#ifdef FORMATDIR
    char env[BUFF_SIZE];
    snprintf(env, BUFF_SIZE, "BABEL_LIBDIR=%s", FORMATDIR);
    putenv(env);
#endif  

  std::ifstream ifs(GetFilename("canonstable.can").c_str());
  OB_REQUIRE( ifs );


  OBMol mol;
  OBConversion conv;
  conv.SetInFormat("smi");
  conv.SetOutFormat("can");

  std::string line;
  while (std::getline(ifs, line)) {
    OB_REQUIRE( conv.ReadString(&mol, line.c_str()) );

    std::vector<OBAtom*> atoms;
    FOR_ATOMS_OF_MOL(atom, mol)
      atoms.push_back(&*atom);

    for (int i = 0; i < 5; ++i) {
      // shuffle the atoms
      std::random_shuffle(atoms.begin(), atoms.end());
      mol.RenumberAtoms(atoms);

      // get can smiles
      mol.SetTitle("");
      std::string cansmi = conv.WriteString(&mol, true);
      // comapare with ref
      if (cansmi != line) {
        cout << "ref = " << line << endl;
        cout << "can = " << cansmi << endl;
        OB_ASSERT( cansmi == line );
      }
    }
  }
 
  return 0;
}


#include "obtest.h"
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/atom.h>
#include <openbabel/obiter.h>

#include <iostream>
#include <vector>
#include <algorithm>

using std::cout;
using std::endl;
using namespace OpenBabel;


int canonstabletest(int argc, char *argv[])
{

  // Define location of file formats for testing
#ifdef FORMATDIR
    char env[BUFF_SIZE];
    snprintf(env, BUFF_SIZE, "BABEL_LIBDIR=%s", FORMATDIR);
    putenv(env);
#endif  

  std::ifstream ifs(OBTestUtil::GetFilename("canonstable.can").c_str());
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


#include <iostream>
#include <fstream>
#include <algorithm>

#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/graphsym.h>
#include <openbabel/builder.h>
#include <openbabel/parsmart.h>
#include <openbabel/elements.h>

using namespace std;
using namespace OpenBabel;

int main(int argc, char** argv) {
  if (argc != 2) {
    cerr << "Usage: " << argv[0] << " <fragment list>" << endl;
    exit(EXIT_FAILURE);
  }

  ifstream ifs(argv[1]);
  if(!ifs) {
    cerr << "Falied to open " << argv[1] << endl;
    exit(EXIT_FAILURE);
  }

  char buffer[BUFF_SIZE];
  vector<string> vs;
  while (ifs.getline(buffer, BUFF_SIZE)) {
    int pos = ifs.tellg();
    if (buffer[0] != '#') {
      tokenize(vs, buffer);
      if (vs.size() == 1) { // SMARTS pattern
        cout << vs[0] << " " << pos << endl;
      }
    }
    pos = ifs.tellg();
  }
}


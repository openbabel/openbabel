#include <openbabel/mol.h>
#include <openbabel/obconversion.h>

#include <iostream>
#include <fstream>
#include <iterator>
#include <string>
#include <map>
#include <vector>

using namespace std;
using namespace OpenBabel;

int main(int argc, char *argv[])
{
  // turn off slow sync with C-style output (we don't use it anyway).
  std::ios::sync_with_stdio(false);

  OBConversion conv;
  conv.SetOptions("O", conv.OUTOPTIONS);
  conv.SetOutFormat("can");

  ifstream ifs;
  OBMol mol;

  if (argc < 2)
  {
    cerr << "Usage: obtorsion <file>" << endl;
    return(-1);
  }

  map<string, vector<double> > torsion_angles;
  for (int i = 1; i < argc; i++) {
    cerr << " Reading file " << argv[i] << endl;

    OBFormat *inFormat = conv.FormatFromExt(argv[i]);
    if(inFormat==NULL || !conv.SetInFormat(inFormat))
    {
      cerr << " Cannot read file format for " << argv[i] << endl;
      continue; // try next file
    }

    ifs.open(argv[i]);

    if (!ifs)
    {
      cerr << "Cannot read input file: " << argv[i] << endl;
      continue;
    }


    while(ifs.peek() != EOF && ifs.good())
    {
      conv.Read(&mol, &ifs);
      mol.DeleteHydrogens();

      FOR_BONDS_OF_MOL(bond, mol) {
        if(bond->IsRotor()) {
          OBBitVec atomsToCopy;
          OBAtom *atom = bond->GetBeginAtom();
          FOR_NBORS_OF_ATOM(a, &*atom) {
            atomsToCopy.SetBitOn(a->GetIdx());
          }
          atom = bond->GetEndAtom();
          FOR_NBORS_OF_ATOM(a, &*atom) {
            atomsToCopy.SetBitOn(a->GetIdx());
          }
          OBMol mol_copy;
          mol.CopySubstructure(mol_copy, &atomsToCopy);
          string smiles = conv.WriteString(&mol_copy, true);

          OBAtom* b = bond->GetBeginAtom();
          OBAtom* c = bond->GetEndAtom();
          OBAtom* a = NULL;
          FOR_NBORS_OF_ATOM(t, &*b) {
            a = &*t;
            if(a != c)
              break;
          }
          OBAtom* d = NULL;
          FOR_NBORS_OF_ATOM(t, &*c) {
            d = &*t;
            if(d != b)
              break;
          }
          double angle = mol.GetTorsion(a, b, c, d);
          if(torsion_angles.count(smiles) == 0) {
            vector<double> t(36, 0);
            torsion_angles[smiles] = t;
          }
          torsion_angles[smiles][static_cast<size_t>(angle+180)/10] += 1;
        }
      }
    }
  }

  map<string, vector<double> >::iterator it;
  for(it=torsion_angles.begin(); it!=torsion_angles.end(); it++) {
    vector<double>::iterator j;
    int max_index = distance(it->second.begin(), max_element(it->second.begin(), it->second.end()));
    cout << it->first << " " << max_index * 10 - 180 << endl;
  }
}

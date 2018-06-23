/**********************************************************************
obfragment - Generate coordinate database of ring fragments

Copyright (C) 2007 Geoffrey R. Hutchison
 
This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

// used to set import/export for Cygwin DLLs
#ifdef WIN32
#define USING_OBDLL
#endif

#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/graphsym.h>
#include <openbabel/builder.h>
#include <openbabel/parsmart.h>
#include <openbabel/data.h>
#include <openbabel/math/align.h>

#if !HAVE_STRNCASECMP
extern "C" int strncasecmp(const char *s1, const char *s2, size_t n);
#endif

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <queue>
#include <map>

using namespace std;
using namespace OpenBabel;

class UnionFind {
  vector<int> table;
  int root(int x) {
    return table[x] < 0 ? x : table[x] = root(table[x]);
  }
public:
  UnionFind(size_t size) : table(size, -1) { }
  void unite(int x, int y) {
    x = root(x);
    y = root(y);
    if (x != y) {
      if (table[y] < table[x]) swap(x, y);
      if (table[x] == table[y]) table[x]--;
      table[y] = x;
    }
  }
  bool isRoot(int x) const {
    return table[x] < 0;
  }
};

int main(int argc, char *argv[])
{
  // turn off slow sync with C-style output (we don't use it anyway).
  std::ios::sync_with_stdio(false);

  OBConversion conv, conv_sdf;
  conv.SetOptions("O", conv.OUTOPTIONS);
  OBFormat *inFormat, *canFormat;
  OBMol mol;
  ifstream ifs;
  vector<OBMol> fragments;
  unsigned int fragmentCount = 0; // track how many in library -- give a running count
  map<string, int> index; // index of cansmi
  string currentCAN, currentSMARTS;
  unsigned int size;
  OBAtom *atom;
  OBBond *bond;
  char buffer[BUFF_SIZE];
  map<string, vector<OBMol> > fragment_list;

  canFormat = conv.FindFormat("can"); // Canonical SMILES format
  conv.SetOutFormat(canFormat);
  canFormat = conv.FindFormat("sdf");
  conv_sdf.SetOutFormat(canFormat);

  if (argc < 2)
  {
    cout << "Usage: obfragment <file>" << endl;
    return(-1);
  }

  for (int i = 1; i < argc; i++) {
    cerr << " Reading file " << argv[i] << endl;

    inFormat = conv.FormatFromExt(argv[i]);
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

      size = mol.NumAtoms();
      OBBitVec atomsToCopy(size+1);
      for (unsigned int i = 1; i <= size; ++i) {
        atom = mol.GetAtom(i);
        atomsToCopy.SetBitOn(atom->GetIdx());
      }

      size = mol.NumBonds();
      OBBitVec bondsToExclude(size);
      for (unsigned int i = 0; i < size; ++i) {
        bond = mol.GetBond(i);
        if (bond->IsRotor()) {
          bondsToExclude.SetBitOn(bond->GetIdx());
        }
      }

      OBMol mol_copy;
      mol.CopySubstructure(mol_copy, &atomsToCopy, &bondsToExclude);
      fragments = mol_copy.Separate(); // Copies each disconnected fragment as a separate
      for (unsigned int i = 0; i < fragments.size(); ++i)
      {
        if (fragments[i].NumAtoms() < 3) // too small to care
          continue;

        string smiles = conv.WriteString(&fragments[i], true);

        if (fragment_list.count(smiles) == 0) {
          vector<OBMol> v;
          v.push_back(fragments[i]);
          fragment_list[smiles] = v;
        } else {
          fragment_list[smiles].push_back(fragments[i]);
        }
      }
    } // while reading molecules (in this file)
    ifs.close();
    ifs.clear();
  } // while reading files
 
  map<string, vector<OBMol> >::iterator i;
  for (i = fragment_list.begin(); i != fragment_list.end(); i++) {
    // OK, now retrieve the canonical SMILES ordering for the fragment
    UnionFind uf(i->second.size());
    cout << i->first << ": " << i->second.size() << endl;
    vector<OBMol>::iterator j, j2;
    for (j = i->second.begin(); j != i->second.end(); j++) {
      for (j2 = i->second.begin(); j2 != i->second.end(); j2++) {
        //if(j == j2) continue;
        OBAlign aln(*j, *j2);
        aln.Align();
        if (aln.GetRMSD() < 2.0) {
          uf.unite(j - i->second.begin(), j2 - i->second.begin());
        }
      }
    }
    for (size_t idx = 0; idx < i->second.size(); idx++) {
      if (uf.isRoot(idx)) {
        OBPairData *pd = dynamic_cast<OBPairData*>(i->second[idx].GetData("SMILES Atom Order"));
        istringstream iss(pd->GetValue());
        vector<unsigned int> canonical_order;
        canonical_order.clear();
        copy(istream_iterator<unsigned int>(iss),
            istream_iterator<unsigned int>(),
            back_inserter<vector<unsigned int> >(canonical_order));

        // Write out an XYZ-style file with the CANSMI as the title
        cout << i->first << '\n'; // endl causes a flush

        unsigned int order;
        OBAtom *atom;

        i->second[idx].Center(); // Translate to the center of all coordinates
        i->second[idx].ToInertialFrame(); // Translate all conformers to the inertial frame-of-reference.

        for (unsigned int index = 0; index < canonical_order.size(); ++index) {
          order = canonical_order[index];
          atom = i->second[idx].GetAtom(order);

          snprintf(buffer, BUFF_SIZE, "%d %9.3f %9.3f %9.3f\n",
              atom->GetAtomicNum(),
              atom->x(), atom->y(), atom->z());
          cout << buffer;
        }
      }
    }
  }
}

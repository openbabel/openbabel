/**********************************************************************
  obfragment - Generate coordinate database of ring fragments

  Copyright (C) 2007 Geoffrey R. Hutchison

  This file is part of the Open Babel project.
  For more information, see <http://openbabel.sourceforge.net/>

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
#include <algorithm>
#include <queue>
#include <map>

using namespace std;
using namespace OpenBabel;

int main(int argc, char *argv[])
{
  // turn off slow sync with C-style output (we don't use it anyway).
  std::ios::sync_with_stdio(false);

  OBConversion conv;
  conv.SetOptions("O", conv.OUTOPTIONS);
  OBFormat *inFormat, *canFormat;
  ifstream ifs;
  OBMol mol;
  OBAtom *atom;
  OBBond *bond;
  char buffer[BUFF_SIZE];
  map<string, vector<OBMol> > fragment_list;
  map<string, int> fragment_count;
  vector<pair<int, string> > fragment_size;

  canFormat = conv.FindFormat("can"); // Canonical SMILES format
  conv.SetOutFormat(canFormat);

  if (argc < 2)
  {
    cerr << "Usage: obfragment <file>" << endl;
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

      unsigned int size = mol.NumAtoms();
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
      vector<OBMol> fragments = mol_copy.Separate(); // Copies each disconnected fragment as a separate
      for (unsigned int i = 0; i < fragments.size(); ++i) {
        if (fragments[i].NumAtoms() < 5) // too small to care
          continue;

        string smiles = conv.WriteString(&fragments[i], true);

        if (fragment_count.count(smiles) == 0) {
          vector<OBMol> v;
          v.push_back(fragments[i]);
          fragment_list[smiles] = v;
          fragment_count[smiles] = 1;
        } else {
          fragment_count[smiles]++;
          vector<OBMol> &f = fragment_list[smiles];
          bool isDifferent = true;
          for (vector<OBMol>::iterator j = f.begin(); j != f.end(); j++) {
            OBAlign aln(fragments[i], *j, false, false);
            aln.Align();
            if (aln.GetRMSD() < 5.0) {
              isDifferent = false;
              break;
            }
          }
          if(isDifferent)
            fragment_list[smiles].push_back(fragments[i]);
        }
      } 
    } // while reading molecules (in this file)
    ifs.close();
    ifs.clear();
  } // // while reading files

  // sort fragments by the number of molecules
  for (map<string, vector<OBMol> >::iterator i = fragment_list.begin(); i != fragment_list.end(); i++) {
    fragment_size.push_back(make_pair<int, string>(i->second[0].NumAtoms(), i->first));
  }
  sort(fragment_size.rbegin(), fragment_size.rend());

  for (vector<pair<int, string> >::iterator i = fragment_size.begin(); i != fragment_size.end(); i++) {
    if (fragment_count[i->second] < 3) continue;
    // OK, now retrieve the canonical SMILES ordering for the fragment
    vector<OBMol>& fragments = fragment_list[i->second];
    vector<OBMol> t;
    for (size_t idx = 0; idx < fragments.size(); idx++) {
      t.push_back(fragments[idx]);
      OBPairData *pd = dynamic_cast<OBPairData*>(fragments[idx].GetData("SMILES Atom Order"));
      istringstream iss(pd->GetValue());
      vector<unsigned int> canonical_order;
      canonical_order.clear();
      copy(istream_iterator<unsigned int>(iss),
          istream_iterator<unsigned int>(),
          back_inserter<vector<unsigned int> >(canonical_order));

      // Write out an XYZ-style file with the CANSMI as the title
      cout << i->second << '\n'; // endl causes a flush

      unsigned int order;
      OBAtom *atom;

      fragments[idx].Center(); // Translate to the center of all coordinates
      fragments[idx].ToInertialFrame(); // Translate all conformers to the inertial frame-of-reference.

      for (unsigned int index = 0; index < canonical_order.size(); ++index) {
        order = canonical_order[index];
        atom = fragments[idx].GetAtom(order);

        snprintf(buffer, BUFF_SIZE, "%d %9.3f %9.3f %9.3f\n",
            atom->GetAtomicNum(),
            atom->x(), atom->y(), atom->z());
        cout << buffer;
      }
    }
  }
}

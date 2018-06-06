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

string RewriteSMILES(const string smiles);

int main(int argc,char *argv[])
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
              
            currentCAN = conv.WriteString(&fragments[i], true);
            currentSMARTS = currentCAN;
            std::stringstream stitle;
            stitle << currentSMARTS << "\t" << fragments[i].NumAtoms() << "\t" << fragments[i].NumBonds();
            std::string title = stitle.str();
            stitle << "\t" << currentCAN;
            if (index.find(title) != index.end()) { // already got this
              index[title] += 1; // add to the count for bookkeeping
              continue;
            }

            index[title] = 1; // don't ever write this ring fragment again

            // OK, now retrieve the canonical SMILES ordering for the fragment
            OBPairData *pd = dynamic_cast<OBPairData*>(fragments[i].GetData("SMILES Atom Order"));
            /*cout << "Canonical order " << pd->GetValue() << "\n";*/
            istringstream iss(pd->GetValue());
            vector<unsigned int> canonical_order;
            canonical_order.clear();
            copy(istream_iterator<unsigned int>(iss),
                 istream_iterator<unsigned int>(),
                 back_inserter<vector<unsigned int> >(canonical_order));

            // Write out an XYZ-style file with the CANSMI as the title
            //cout << fragments[i].NumAtoms() << '\n';
            //cout << stitle.str() << '\n'; // endl causes a flush
            cout << currentSMARTS << '\n'; // endl causes a flush

            unsigned int order;
            OBAtom *atom;

            fragments[i].Center(); // Translate to the center of all coordinates
            fragments[i].ToInertialFrame(); // Translate all conformers to the inertial frame-of-reference.

            for (unsigned int index = 0; index < canonical_order.size(); 
                 ++index) {
              order = canonical_order[index];
              atom = fragments[i].GetAtom(order);
              
              snprintf(buffer, BUFF_SIZE, "%d %9.3f %9.3f %9.3f\n",
                       atom->GetAtomicNum(),
                       atom->x(), atom->y(), atom->z());
              cout << buffer;
            }
          }
        fragments.clear();
        if (index.size() > fragmentCount) {
          fragmentCount = index.size();
          cerr << " Fragments: " << fragmentCount << endl;
        }

      } // while reading molecules (in this file)
    ifs.close();
    ifs.clear();
  } // while reading files

  // loop through the map and output frequencies
  map<string, int>::const_iterator indexItr;
  for (indexItr = index.begin(); indexItr != index.end(); ++indexItr) {
      cerr << (*indexItr).second << " INDEX " << (*indexItr).first << "\n";
  }

  /*priority_queue<pair<int, string> > freq;
  int total = 0;
  map<string, int>::const_iterator indexItr;
  for (indexItr = index.begin(); indexItr != index.end(); ++indexItr) {
      freq.push(make_pair((*indexItr).second, (*indexItr).first));
      total += (*indexItr).second;
  }
  for(int i = 0; i < 25; ++i) {
      if(freq.empty()) break;
      pair<int, string> f = freq.top();
      freq.pop();
      stringstream ss(f.second);
      string smarts;
      ss >> smarts;
      cerr << f.second << "\t" << 100.0*f.first/total << "\n";
  }*/
    
  return(0);
}

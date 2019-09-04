/**********************************************************************
obsort -- sort the fragment database

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

#if !HAVE_STRNCASECMP
extern "C" int strncasecmp(const char *s1, const char *s2, size_t n);
#endif

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <algorithm>

using namespace std;
using namespace OpenBabel;

bool CompareSMILES(const string &a, const string &b)
{ 
  // Originally a.length() > b.length()
  // Now it uses the actual number of atoms and number of bonds
  vector<string> va;
  tokenize(va, a.c_str());
  vector<string> vb;
  tokenize(vb, b.c_str());

  // Sort based on the number of atoms, then the number of bonds, then the length
  if (atoi(va[1].c_str()) > atoi(vb[1].c_str()))
    return true;
  else if (atoi(va[1].c_str()) == atoi(vb[1].c_str())) {
    if (atoi(va[2].c_str()) > atoi(vb[2].c_str()))
      return true;
    else if (atoi(va[2].c_str()) == atoi(vb[2].c_str()))
      return a.length() > b.length();
  }
  else
    return false;
}

int main(int argc,char *argv[])
{
  // turn off slow sync with C-style output (we don't use it anyway).
  std::ios::sync_with_stdio(false);

  if (argc != 3)    {
    cout << "Usage: obsort <fragments> <freq>" << endl;
    return(-1);
  }
  
  ifstream ifs(argv[1]);
  if (!ifs)    {
    cout << "Bail out! Cannot read input file!" << endl;
    return(-1);
  }
  
  OBConversion conv(&ifs, &cout);
  OBFormat* pFormat;
  
  pFormat = conv.FormatFromExt(argv[1]);
  if ( pFormat == NULL )    {
    cout << "Bail out! Cannot read file format!" << endl;
    return(-1);
  }
  
  if (! conv.SetInAndOutFormats(pFormat, pFormat))    {
    cout << "Bail out! File format isn't loaded" << endl;
    return (-1);
  }
  
  ifstream freq(argv[2]);
  if (!freq) {
    cout << "Bail out! Cannot read list of frequencies!" << endl;
  }
  
  OBMol *mol;
  OBAtom *atom;
  OBSmartsPattern sp;
  map<string, OBMol*> molIndex; // index of molecules

  // First read in the molecules and hash them
  int molCount = 0;
  string temp, title;
  vector<string> vs;
  while(ifs.peek() != EOF && ifs.good())    {
    mol = new OBMol;
    conv.Read(mol);
    
    /*if (!mol->Has3D())
      continue; // invalid coordinates*/
    
    temp = std::string(mol->GetTitle());
    
    //cout << "Title" << temp << endl;
    OBMol clone;
    clone = OBMol(*mol);
    temp = Trim(temp);
    if (!sp.Init(temp)) // not a valid SMARTS pattern
      continue;

    tokenize(vs, temp);
    title = vs[0] + "\t" + vs[1] + "\t" + vs[2];
    molIndex[title] = mol;

    // Give some progress
    molCount++;
    if (molCount % 500 == 0)
      cerr << "Count: " << molCount << endl;
  } // while reading molecules
  
  // Now read in the list of frequencies and SMILES
  // (We sort the SMILES by length -- it's already sorted by frequency)
  vector<string> smilesIndex;
  string line;
  char buffer[BUFF_SIZE];
  do {
    getline(freq, line); // Count INDEX CanSMI
    tokenize(vs, line);
    string title = vs[2] + "\t" + vs[3] + "\t" + vs[4];
    if (molIndex.find(title) == molIndex.end())
      continue;
    
    smilesIndex.push_back(title);
  }  while (freq.peek() != EOF && freq.good());
  sort(smilesIndex.begin(), smilesIndex.end(), CompareSMILES);
   
  // Finally, we loop through the index of SMILES and write the molecules
  for (vector<string>::iterator j = smilesIndex.begin(); j != smilesIndex.end(); ++j) {
    mol = molIndex[(*j)];
    if (!mol) // Shouldn't happen, but it's good defensive programming
      continue;
      
    mol->Center();
//    cout << mol->NumAtoms() << "\n";
    tokenize(vs, *j);
    cout << vs[0] << "\n";
    for(unsigned int i = 1;i <= mol->NumAtoms(); i++)
    {
      atom = mol->GetAtom(i);
      snprintf(buffer, BUFF_SIZE, "%7.3f%7.3f%7.3f\n",
               atom->x(), atom->y(), atom->z());
      cout << buffer;
    }

//    delete(mol);
    mol = NULL;
  }

  return(0);
}

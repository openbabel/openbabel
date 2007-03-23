/**********************************************************************
 inchiwrite.cpp - tests for writing InChIs

 This file is part of the Open Babel project.
 For more information, see <http://openbabel.sourceforge.net/>

 Copyright (C) 2007 by Chris Morley

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


#include <iostream>
#include <fstream>
#include <sstream>
#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>

using namespace std;
namespace OpenBabel{
extern string GetInChI(istream& is);
}
using namespace OpenBabel;

int main(int argc, char* argv[])
{
  if(argc !=3)
  {
    cout << "Usage:  inchiwrite Testmols.xxx Results.txt\n"
    "Testmols.xxx can be any format but must have the appropriate extension.\n"
    "Results.txt is a text file with the correct the correct InChIs.\n"
    "It can contain other text interspersed." << endl;
    return 1;
  }
  ifstream ifs(argv[1]);
  if(!ifs)
    cerr << "Cannot open " << argv[1];
  ifstream results(argv[2]);
  if(!results)
    cerr << "Cannot open " << argv[2];

  stringstream ssout;
  OBConversion conv(&ifs, &ssout);
  OBFormat* pFormat = conv.FormatFromExt(argv[1]);
  if(!pFormat)
  {
    cerr << "Did not recognize the Format of " << argv[1] << endl;
    return 1;
  }
  conv.SetInFormat(pFormat);
  if(!conv.SetOutFormat("inchi"))
  {  
    cerr << "Did not find InChI Format" << endl;
    return 1;
  }

  conv.AddOption("w",OBConversion::OUTOPTIONS); //suppress routine warnings
  int count = conv.Convert();
  cout << count << " molecules converted\n" << endl;
  if(count<=0)
    return 2;

  //Compare ssout with the correct InChIs
  int nfail=0;
  int n, skipped=0;
  for(n=1;ssout;++n)
  {
    string correct = GetInChI(results);
    string generated;
    ssout >> generated;
    if(correct.empty() || generated.empty())
      break;
    if(correct.size()<9)//skip molecules with an empty InChI
    {
      skipped++;
      continue;
    }
    if(generated != correct)
    {
      cout << "Mismatch in molecule #" << n << ". generated / correct\n"
        << generated << '\n' << correct << '\n' << endl;
      nfail++;
    }
  }
  cout << n-1-skipped<< " molecules compared" << endl;
  cout << nfail << " failures" << endl;
  return 0;
}

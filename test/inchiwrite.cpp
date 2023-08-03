/**********************************************************************
 inchiwrite.cpp - tests for writing InChIs

 This file is part of the Open Babel project.
 For more information, see <http://openbabel.org/>

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

#include <cstdlib>
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

  // Define location of file formats for testing
  #ifdef FORMATDIR
    char env[BUFF_SIZE];
    snprintf(env, BUFF_SIZE, "BABEL_LIBDIR=%s", FORMATDIR);
    putenv(env);
  #endif

  ifstream ifs(argv[1]);
  if(!ifs) { 
    cout << "# Cannot open " << argv[1];
    return 1;
  }
  ifstream results(argv[2]);
  if(!results) {
    cout << "# Cannot open " << argv[2];
    return 1;
  }

  stringstream ssout;
  OBConversion conv(&ifs, &ssout);
  OBFormat* pFormat = conv.FormatFromExt(argv[1]);
  if(!pFormat)
  {
    cout << "# Skipped. Did not recognize the format of " << argv[1] << endl;
    return 1;
  }
  conv.SetInFormat(pFormat);
  if(!conv.SetOutFormat("inchi"))
  {  
    cout << "# Skipped. Did not find InChI Format" << endl;
    return 1;
  }

  conv.AddOption("w",OBConversion::OUTOPTIONS); //suppress routine warnings

  /*In OB version 2.2.0 a set of 'recommended' InChI options is used by default.
    The compiled Windows program wInChI.exe, which is used to generate the results file
    was not available for InChI version 1.02, the beta of which is used in OB.
    So the following option turns off all the InChI options. This needs to be revisted
    when InChI 1.02 is released.
  */
  conv.AddOption("n",OBConversion::OUTOPTIONS);

  int count = conv.Convert();
  if(count<=0) {
    cout << "# Skipped. Did not convert any molecules" << endl;
    return 2;
  }

  cout << "1.." << count << "\n # molecules converted\n";

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
      cout << "ok " << n << " # ignored empty InChI\n";
      continue;
    }

    if(generated != correct)
    {
      cout << "Not ok " << n << " # Mismatch in molecule #" 
           << n << ". generated / correct\n"
           << "# " << generated << "\n# " << correct << '\n';
      nfail++;
    }
    else
      cout << "ok " << n << " # InChI matched expected\n";
  }
  cout << n-1-skipped<< " molecules compared" << endl;
  cout << "# " << skipped << " # empty molecules skipped" << endl;
  cout << nfail << " failures" << endl;
  return 0;
}

/**********************************************************************
obtautomer - Enumerate tautomer smiles and canonical tautomer smiles

Copyright (C) 2011 Tim Vandermeersch
 
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
#include <cstdlib>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/tautomer.h>

using namespace std;

/**
 * TautomerFunctor to print out smiles for each found tautomer.
 */
class Functor : public OpenBabel::TautomerFunctor
{
  public:
    bool foundTautomers;

    Functor() : foundTautomers(false) 
    {
    }
  
    void operator()(OpenBabel::OBMol *mol)
    {
      foundTautomers = true;
      OpenBabel::OBConversion conv;
      conv.SetOutFormat("can");
      std::cout << conv.WriteString(mol);
    }
};


int main(int argc,char **argv)
{
  char *program_name= argv[0];
  int c;
  char *FileIn = NULL;

  if (argc < 2) {
    string err = "Usage: ";
    err += program_name;
    err += " [-c] <filename>\n";
    err += "-c: Canonical tautomer only\n";
    cerr << err; //Why not do directly
    exit(-1);
  } else {
    FileIn  = argv[argc-1];
  }

  // Find Input filetype
  OpenBabel::OBConversion conv;
  OpenBabel::OBFormat *format = conv.FormatFromExt(FileIn);
    
  if (!format || !conv.SetInFormat(format)) {
    cerr << program_name << ": cannot read input format!" << endl;
    exit (-1);
  }
    
  if (!conv.SetOutFormat("can")) {
    cerr << program_name << ": cannot find output format!" << endl;
    exit (-1);
  }

  ifstream ifs;

  // Read the file
  ifs.open(FileIn);
  if (!ifs) {
    cerr << program_name << ": cannot read input file!" << endl;
    exit (-1);
  }
  
  OpenBabel::OBMol mol;


  for (c = 1;; ++c) {
      mol.Clear();
      conv.Read(&mol, &ifs);
      if (mol.Empty())
        break;

      if (std::string(argv[1]) == "-c") {
        CanonicalTautomer(&mol);
        std::cout << conv.WriteString(&mol);
      } else {
        Functor f;
        EnumerateTautomers(&mol, f);

        if (!f.foundTautomers)
          std::cout << conv.WriteString(&mol);
      }
      
  } // end for loop
  
  return(0);
}


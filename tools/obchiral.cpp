/**********************************************************************
OBChiral
-Lists chiral centers in a molecule
Copyright 2005 - Nick England
 
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
#include <openbabel/chiral.h>
#ifndef _MSC_VER
  #include <unistd.h>
#endif

using namespace std;
using namespace OpenBabel;

// PROTOTYPES /////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
//! \brief Compute some properties easy to access from open babel
//


int main(int argc,char **argv)
{
  char *program_name= argv[0];
  char *FileIn = NULL;

  if (argc != 2)
    {
      cout << "Usage: " << program_name << " <filename>" << endl;
      exit(-1);
    }
  else
    {
      FileIn  = argv[1];
      //   const char* p = strrchr(FileIn,'.');
    }

  // Find Input filetype
  OBConversion conv;
  OBFormat *format = conv.FormatFromExt(FileIn);
    
  if (!format || !conv.SetInAndOutFormats(format, format))
    {
      cerr << program_name << ": cannot read input format!" << endl;
      exit (-1);
    }

  ifstream ifs;

  // Read the file
  ifs.open(FileIn);
  if (!ifs)
    {
      cerr << program_name << ": cannot read input file!" << endl;
      exit (-1);
    }

  OBMol mol;
  OBAtom *atom;

  for (int c=1;;++c) // big for loop (replace with do while?)
    {
      mol.Clear();
      conv.Read(&mol, &ifs);
      if (mol.Empty())
        break;
      cout << "Molecule "<< c << ": " << mol.GetTitle() << endl;
      //mol.FindChiralCenters(); // labels all chiral atoms
      vector<OBAtom*>::iterator i; // iterate over all atoms
      for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
        {
          if(!atom->IsChiral())continue; // aborts if atom isn't chiral
          cout << "Atom " << atom->GetIdx() << " Is Chiral ";
          cout << atom->GetType()<<endl;
        
          OBChiralData* cd = (OBChiralData*)atom->GetData(OBGenericDataType::ChiralData);
        
          if (cd){
            vector<unsigned int> x=cd->GetAtom4Refs(input);
            size_t n=0;
            cout <<"Atom4refs:";
            for (n=0;n<x.size();++n)
              cout <<" "<<x[n];
            cout <<endl;
          }
          else{cd=new OBChiralData;atom->SetData(cd);}
          vector<unsigned int> _output;
          unsigned int n;
          for(n=1;n<5;++n) _output.push_back(n);
          cd->SetAtom4Refs(_output,output);
          /* // MOLV3000 uses 1234 unless an H then 123H
             if (atom->GetHvyValence()==3)
             {
             OBAtom *nbr;
             int Hid=1000;// max Atom ID +1 should be used here
             vector<unsigned int> nbr_atms;
             vector<OBBond*>::iterator i;
             for (nbr = atom->BeginNbrAtom(i);nbr;nbr = atom->NextNbrAtom(i))
             {
             if (nbr->IsHydrogen()){Hid=nbr->GetIdx();continue;}
             nbr_atms.push_back(nbr->GetIdx());
             }
             sort(nbr_atms.begin(),nbr_atms.end());
             nbr_atms.push_back(Hid);
             OBChiralData* cd=(OBChiralData*)atom->GetData(OBGenericDataType::ChiralData);
             cd->SetAtom4Refs(nbr_atms,output);   
             } 
             else if (atom->GetHvyValence()==4)
             {
             OBChiralData* cd=(OBChiralData*)atom->GetData(OBGenericDataType::ChiralData);
             vector<unsigned int> nbr_atms;
             int n;
             for(n=1;n<5;++n)nbr_atms.push_back(n);
             cd->SetAtom4Refs(nbr_atms,output); 
             } */
    /* FIXME          
          if (!mol.HasNonZeroCoords())
            {
              cout << "Calcing 0D chirality "<< CorrectChirality(mol,atom)<<endl;
            }
          else {
            cout << "Volume= "<< CalcSignedVolume(mol,atom) << endl;
            OBChiralData* cd=(OBChiralData*)atom->GetData(OBGenericDataType::ChiralData);
            size_t n;
            vector<unsigned int> refs=cd->GetAtom4Refs(output);
            cout<<"Atom refs=";
            for(n=0;n<refs.size();++n)cout<<" "<<refs[n];
            cout<<endl;
          }
          cout << "Clockwise? " << atom->IsClockwise() << endl;
          */
        } // end iterating over atoms

    } // end big for loop

  return(0);
} // end main


/* obchiral man page*/
/** \page obchiral print molecular chirality information
*
* \n
* \par SYNOPSIS
*
* \b obchiral \<filename\>
*
* \par DESCRIPTION
*
* The obchiral program is a tool to print the chirality information
* for all molecules in a file. It also serves as example code for using the
* Open Babel library (libopenbabel).
* 
* \par EXAMPLES
*
*   obchiral pyridines.sdf
*
* \par AUTHORS
*
* The obchiral program was contributed by \b Nick \b England.
*
* Open Babel is currently maintained by \b Geoff \b Hutchison, \b Chris \b Morley and \b Michael \b Banck.
*
* For more contributors to Open Babel, see http://openbabel.org/THANKS.shtml
*
* \par COPYRIGHT
*  Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
*  Some portions Copyright (C) 2001-2005 by Geoffrey R. Hutchison \n \n
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation version 2 of the License.\n \n
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
* \par SEE ALSO
*   The web pages for Open Babel can be found at: http://openbabel.org/ \n
**/

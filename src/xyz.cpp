/**********************************************************************
Copyright (C) 2000 by OpenEye Scientific Software, Inc.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include <stdlib.h>
#include "mol.h"

using namespace std;

namespace OpenBabel {

bool ReadXYZ(istream &ifs,OBMol &mol,const char *title)
{
  char buffer[BUFF_SIZE];

  // Read the first line, which should contain the number of atoms in
  // the molecule.
  if (!ifs.getline(buffer,BUFF_SIZE)) {
    cerr << "WARNING: Problems reading an XYZ file, method 'bool ReadXYZ(istream &,OBMol &,const char *)'" << endl
	 << "  Could not read the first line, file error." << endl;
    return(false);
  }
  unsigned int natoms;	// [ejk] assumed natoms could not be -ve
  if (sscanf(buffer,"%d", &natoms) == 0) {
    cerr << "WARNING: Problems reading an XYZ file, method 'bool ReadXYZ(istream &,OBMol &,const char *)'" << endl
	 << "  Problems reading the first line. The first line must contain the number of atoms." << endl
	 << "  OpenBabel found the line '" << buffer << "'" << endl
	 << "  which could not be interpreted as a number." << endl;
    return(false);
  }
  if (!natoms) {
    cerr << "WARNING: Problems reading an XYZ file, method 'bool ReadXYZ(istream &,OBMol &,const char *)'" << endl
	 << "  Problems reading the first line. The first line must contain the number of atoms." << endl
	 << "  OpenBabel found the number '0' which obviously does not work." << endl;
    return(false);
  }
  mol.ReserveAtoms(natoms);
  
  // The next line contains a title string for the molecule. Use this
  // as the title for the molecule if the line is not
  // empty. Otherwise, use the title given by the calling function.
  if (!ifs.getline(buffer,BUFF_SIZE)) {
    cerr << "WARNING: Problems reading an XYZ file, method 'bool ReadXYZ(istream &,OBMol &,const char *)'" << endl
	 << "  Could not read the second line, file error." << endl;
    return(false);
  }
  if (strlen(buffer) == 0)
    mol.SetTitle(buffer);
  else
    mol.SetTitle(title);
  
  // The next lines contain four items each, separated by white
  // spaces: the atom type, and the coordinates of the atom
  vector<string> vs;
  for (unsigned int i = 1; i <= natoms; i ++) {
    if (!ifs.getline(buffer,BUFF_SIZE)) {
      cerr << "WARNING: Problems reading an XYZ file, method 'bool ReadXYZ(istream &,OBMol &,const char *)'" << endl
	   << "  Could not read line #" << i+2 << ", file error." << endl
	   << "  According to line one, there should be " << natoms << " atoms, and therefore " << natoms+2 << " lines in the file" << endl;
      return(false);
    }
    tokenize(vs,buffer);
    if (vs.size() != 4) {
      cerr << "WARNING: Problems reading an XYZ file, method 'bool ReadXYZ(istream &,OBMol &,const char *)'" << endl
	   << "  Could not read line #" << i+2 << "." << endl
	   << "  OpenBabel found the line '" << buffer << "'" << endl
	   << "  According to the specifications, this line should contain exactly 4 entries, separated by white space." << endl
	   << "  However, OpenBabel found " << vs.size() << " items." << endl;
      return(false);
    }
    
    // Atom Type: get the atomic number from the element table, using
    // the first entry in the currently read line. If the entry makes
    // sense, set the atomic number and leave the atomic type open
    // (the type is then later faulted in when atom->GetType() is
    // called). If the entry does not make sense to use, set the atom
    // type manually, assuming that the author of the xyz-file had
    // something "special" in mind.
    OBAtom *atom  = mol.NewAtom();
    int atomicNum = etab.GetAtomicNum(vs[0].c_str());
    atom->SetAtomicNum(etab.GetAtomicNum(vs[0].c_str())); //set atomic number, or '0' if the atom type is not recognized
    if (atomicNum == 0)
      atom->SetType(vs[0]);

    // Read the atom coordinates
    char *endptr;
    double x = strtod((char*)vs[1].c_str(),&endptr);
    if (endptr == (char*)vs[1].c_str()) {
      cerr << "WARNING: Problems reading an XYZ file, method 'bool ReadXYZ(istream &,OBMol &,const char *)'" << endl
	   << "  Could not read line #" << i+2 << "." << endl
	   << "  OpenBabel found the line '" << buffer << "'" << endl
	   << "  According to the specifications, this line should contain exactly 4 entries, separated by white space." << endl
	   << "  Item #0: Atom Type, Items #1-#3 cartesian coordinates in Angstrom." << endl
	   << "  OpenBabel could not interpret item #1 as a number." << endl;
      return(false);
    }
    double y = strtod((char*)vs[2].c_str(),&endptr);
    if (endptr == (char*)vs[2].c_str()) {
      cerr << "WARNING: Problems reading an XYZ file, method 'bool ReadXYZ(istream &,OBMol &,const char *)'" << endl
	   << "  Could not read line #" << i+2 << "." << endl
	   << "  OpenBabel found the line '" << buffer << "'" << endl
	   << "  According to the specifications, this line should contain exactly 4 entries, separated by white space." << endl
	   << "  Item #0: Atom Type, Items #1-#3 cartesian coordinates in Angstrom." << endl
	   << "  OpenBabel could not interpret item #2 as a number." << endl;
      return(false);
    }
    double z = strtod((char*)vs[3].c_str(),&endptr);
    if (endptr == (char*)vs[3].c_str()) {
      cerr << "WARNING: Problems reading an XYZ file, method 'bool ReadXYZ(istream &,OBMol &,const char *)'" << endl
	   << "  Could not read line #" << i+2 << "." << endl
	   << "  OpenBabel found the line '" << buffer << "'" << endl
	   << "  According to the specifications, this line should contain exactly 4 entries, separated by white space." << endl
	   << "  Item #0: Atom Type, Items #1-#3 cartesian coordinates in Angstrom." << endl
	   << "  OpenBabel could not interpret item #3 as a number." << endl;
      return(false);
    }
    atom->SetVector(x,y,z); //set coordinates
  }

  mol.ConnectTheDots();
  mol.PerceiveBondOrders();
  
  return(true);
}
  
  
bool WriteXYZ(ostream &ofs,OBMol &mol)
{
  unsigned int i;
  char buffer[BUFF_SIZE];
  
  sprintf(buffer,"%d", mol.NumAtoms());
  ofs << buffer << endl;
  sprintf(buffer,"%s\tEnergy: %15.7f", mol.GetTitle(), mol.GetEnergy());
  ofs << buffer << endl;

  OBAtom *atom;
  string str,str1;
  for(i = 1;i <= mol.NumAtoms(); i++)
  {
    atom = mol.GetAtom(i);
    sprintf(buffer,"%3s%15.5f%15.5f%15.5f",
	    etab.GetSymbol(atom->GetAtomicNum()),
	    atom->GetX(),
	    atom->GetY(),
	    atom->GetZ());
    ofs << buffer << endl;
  }

  return(true);
}

}

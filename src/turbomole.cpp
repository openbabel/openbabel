/**********************************************************************
Copyright (C) 2004 by Michael Banck

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

#include <stdlib.h>
#include "mol.h"

using namespace std;

namespace OpenBabel {

void ReadTurbomoleCoords(istream &, OBMol &);

bool ReadTurbomole(istream &ifs, OBMol &mol, const char *title)
{
	string buff;
	string::size_type i;

	if (!getline(ifs, buff)) {
		return(false);
	}
	if (buff != "$coord") {
		cerr << "No coordinates found" << endl;
    		return(false);
	}
	mol.SetTitle(title);
	ReadTurbomoleCoords(ifs, mol);
	// FIXME: handle rest of file
	mol.ConnectTheDots();
	mol.PerceiveBondOrders();
  
	return(true);
}
  
void ReadTurbomoleCoords(istream &ifs, OBMol &mol)
{
	string buff;
	vector<string> vs;
	OBAtom *atom;
	double x, y, z;

	/* Each line contains four items each, separated by white
	   spaces: the coordinates of the atom and the atom type */
	while (getline(ifs, buff)) {
		tokenize(vs, buff); 
		if (vs.size() != 4) {
			// Assume the coordinate block is over
			return; // FIXME: Throw?  
		}
		atom = mol.NewAtom();
		//set atomic number, or '0' if the atom type is not recognized
		atom->SetAtomicNum(etab.GetAtomicNum(vs[3].c_str())); 
		if (etab.GetAtomicNum(vs[3].c_str()) == 0)
			atom->SetType(vs[3]);
		x = atof(vs[0].c_str());
		y = atof(vs[1].c_str());
		z = atof(vs[2].c_str());
		atom->SetVector(x, y, z); 
	}
}
  
bool WriteTurbomole(ostream &ofs,OBMol &mol)
{
	OBAtom *atom;
	vector<OBNodeBase*>::iterator i;

	ofs << "$coord" << endl;
	for (atom = mol.BeginAtom(i); atom; atom = mol.NextAtom(i)) {
		ofs.setf(ios_base::fixed, ios_base::floatfield);
		ofs.precision(14); 
		ofs.width(22);
		ofs << atom->GetX();
		ofs.width(22);
		ofs << atom->GetY(); 
		ofs.width(22);
		ofs << atom->GetZ();
		ofs.width(3);
		ofs << etab.GetSymbol(atom->GetAtomicNum()) << endl;
		ofs.setf(ios_base::fmtflags(0), ios_base::floatfield);
	}
	ofs << "$end" << endl;

	return(true);
}

} // namespace OpenBabel

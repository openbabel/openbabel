/**********************************************************************
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (c) 2001-2003 by Geoffrey R. Hutchison

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include "mol.h"

using namespace std;

namespace OpenBabel
{

bool WriteCache(ostream &ofs,OBMol &mol)
{ 
  char type_name[10];
  char buffer[BUFF_SIZE];
  
  mol.Kekulize();

  ofs << "molstruct88_Apr_30_1993_11:02:29 <molecule> 0x1d00" << endl;
  ofs << "Written by Molecular Editor on <date>" << endl;
  ofs << "Using data dictionary         9/9/93  4:47 AM" << endl;
  ofs << "Version 6" << endl;
  ofs << "local_transform" << endl;
  ofs << "0.100000 0.000000 0.000000 0.000000" << endl;
  ofs << "0.000000 0.100000 0.000000 0.000000" << endl;
  ofs << "0.000000 0.000000 0.100000 0.000000" << endl;
  ofs << "0.000000 0.000000 0.000000 1.000000" << endl;
  ofs << "object_class atom" << endl;
  ofs << "property xyz_coordinates MoleculeEditor angstrom 6 3 FLOAT" << endl;
  ofs << "property anum MoleculeEditor unit 0 1 INTEGER" << endl;
  ofs << "property sym MoleculeEditor noUnit 0 2 STRING" << endl;
  ofs << "property chrg MoleculeEditor charge_au 0 1 INTEGER" << endl;
  ofs << "property rflag MoleculeEditor noUnit 0 1 HEX" << endl;
  ofs << "ID xyz_coordinates             anum sym	chrg rflag" << endl;

  OBAtom *atom;
  vector<OBNodeBase*>::iterator i;
  for(atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
  {
    strcpy(type_name,etab.GetSymbol(atom->GetAtomicNum()));
    
    sprintf(buffer,"%3d %10.6f %10.6f %10.6f %2d %2s %2d 0x7052",
            atom->GetIdx(),
	    atom->x(),
	    atom->y(),
	    atom->z(),
	    atom->GetAtomicNum(),
	    type_name,
	    atom->GetFormalCharge());
    ofs << buffer << endl;
  }

  ofs << "property_flags:" << endl;
  ofs << "object_class bond" << endl;
  ofs << "property rflag MoleculeEditor noUnit 0 1 HEX" << endl;
  ofs << "property type MoleculeEditor noUnit 0 1 NAME" << endl;
  ofs << "property bond_order MoleculeEditor noUnit 4 1 FLOAT" << endl;
  ofs << "ID rflag type bond_order" << endl;

  char bstr[10];
  OBBond *bond;
  vector<OBEdgeBase*>::iterator j;
  for (bond = mol.BeginBond(j);bond;bond = mol.NextBond(j))
    {
      switch (bond->GetBO())
	{
	case 1: strcpy(bstr,"single"); break;
	case 2: strcpy(bstr,"double"); break;
	case 3: strcpy(bstr,"triple"); break;
	default: strcpy(bstr,"weak");
	}

      sprintf(buffer,"%3d 0x7005 %s", bond->GetIdx()+1,bstr);
      ofs << buffer << endl;
    }

  ofs << "property_flags:" << endl;
  ofs << "object_class connector" << endl;
  ofs << "property dflag MoleculeEditor noUnit 0 1 HEX" << endl;
  ofs << "property objCls1 MoleculeEditor noUnit 0 1 NAME" << endl;
  ofs << "property objCls2 MoleculeEditor noUnit 0 1 NAME" << endl;
  ofs << "property objID1 MoleculeEditor noUnit 0 1 INTEGER" << endl;
  ofs << "property objID2 MoleculeEditor noUnit 0 1 INTEGER" << endl;
  ofs << "ID dflag objCls1 objCls2 objID1 objID2" << endl;


  int k;
  for (bond = mol.BeginBond(j),k=1;bond;bond = mol.NextBond(j))
    {
      sprintf(buffer,"%3d 0xa1 atom bond %d %d",
	      k++,bond->GetBeginAtomIdx(),bond->GetIdx()+1);
      ofs << buffer << endl;
      sprintf(buffer,"%3d 0xa1 atom bond %d %d",
	      k++,bond->GetEndAtomIdx(),bond->GetIdx()+1);
      ofs << buffer << endl;
    }

  sprintf(buffer,"property_flags:"); ofs << buffer << endl;
  return(true);
}


}

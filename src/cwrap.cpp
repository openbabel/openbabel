/**********************************************************************
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include "mol.h"
#include "obutil.h"
#include "typer.h"
#include "binary.h"

extern "C" 
{
#include "cwrap.h"
}

using namespace std;
using namespace OpenBabel;

char AtomicNumToOBSymbol(int atomno);

long int ob_make_dbase(char* fname)
{
  return((long int)new OBBinaryDBase (fname));
}

int ob_dbase_get_number_of_molecules(long int dbhandle)
{
  return(((OBBinaryDBase*)dbhandle)->Size());
}

int ob_dbase_get_cmol(int idx,long int mhandle,long int dbhandle)
{
  ((OBBinaryDBase*)dbhandle)->GetMolecule(*((OBMol*)mhandle),idx);

  return(1);
}

long int ob_make_cmol()
{
  return((long int)new OBMol);
}

void ob_delete_cmol(long int handle)
{
  delete (OBMol*)handle;
  handle = 0;
}

int ob_get_cmol_atom_number(int *natoms,long int handle)
{
  *natoms = ((OBMol*)handle)->NumAtoms();
  return(1);
}

int ob_get_cmol_bond_number(int *nbonds,long int handle)
{
  *nbonds = ((OBMol*)handle)->NumBonds();
  return(1);
}

int ob_get_cmol_coordinates(float **c,long int handle)
{
  *c = ((OBMol*)handle)->GetConformer(0);
  return(1);
}

int ob_get_cmol_element(char *c,long int handle)
{
  OBAtom *atom;
  OBMol *mol = (OBMol*)handle;
  vector<OBNodeBase*>::iterator i;

  for(atom = mol->BeginAtom(i);atom;atom = mol->NextAtom(i))
    c[atom->GetIdx()-1] = AtomicNumToOBSymbol(atom->GetAtomicNum());

  return(1);
}

int ob_get_cmol_name(char *name,long int handle)
{
  strcpy(name,((OBMol*)handle)->GetTitle());
  return(1);
}

int ob_get_cmol_conformer_number(int *cnum,long int handle)
{
  *cnum = ((OBMol*)handle)->NumConformers();
  return(1);
}

int ob_get_cmol_conformer(float **c,int num,long int handle)
{
  *c = ((OBMol*)handle)->GetConformer(num);
  return(1);
}


char AtomicNumToOBSymbol(int atomno)
{
  char ele;
  switch (atomno)
    {
    case 1: ele = 'H'; break; //hydrogen
    case 3: ele = 'T'; break; //lithium
    case 5: ele = 'B'; break; //boron
    case 6: ele = 'C'; break; //carbon
    case 7: ele = 'N'; break; //nitrogen
    case 8: ele = 'O'; break; //oxygen
    case 9: ele = 'F'; break; //florine
    case 11: ele = 'D'; break; //sodium
    case 12: ele = 'M'; break; //magnesium
    case 14: ele = 'G'; break; //silicon
    case 15: ele = 'P'; break; //phosphorus
    case 16: ele = 'S'; break; //sulfur
    case 17: ele = 'L'; break; //chlorine
    case 19: ele = 'K'; break; //potassium
    case 25: ele = 'A'; break; //manganese
    case 26: ele = 'E'; break; //iron
    case 29: ele = 'U'; break; //copper
    case 34: ele = 'E'; break; //selenium
    case 35: ele = 'R'; break; //bromine
    case 53: ele = 'I'; break; //iodine
    default: ele = 'Z'; break; //dummy
    }
  return(ele);
}


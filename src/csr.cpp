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

namespace OpenBabel
{

static bool FirstTime = true;
static int MolCount = 1;

static void WriteSize(int,ostream&);
static char *PadString(char*,int);
static void WriteCSRHeader(ostream&,OBMol&);
static void WriteCSRCoords(ostream&,OBMol&);

bool WriteCSR(ostream &ofs,OBMol &mol)
{
  if (FirstTime)
  {
    WriteCSRHeader(ofs,mol);
    FirstTime = false;
  }
  
  WriteCSRCoords(ofs,mol);
  MolCount++;
  
  return(true);
}

void WriteCSRHeader(ostream &ofs,OBMol &mol)
{
  char *molnames;
  int nmol, natom;

  molnames = PadString((char*)mol.GetTitle(),100);

  nmol = 1;
  natom = mol.NumAtoms();

  WriteSize(4*sizeof(char),ofs);
  ofs.write("V33 ",strlen("V33 ")*sizeof(char));
  WriteSize(4*sizeof(char),ofs);
  
  WriteSize(2*sizeof(int),ofs);
  ofs.write((char*)&natom,sizeof(int));
  ofs.write((char*)&nmol,sizeof(int));
  WriteSize(2*sizeof(int),ofs);
  
  WriteSize(100*sizeof(char),ofs);
  ofs.write(molnames,100*sizeof(char));
  WriteSize(100*sizeof(char),ofs);
  
  WriteSize(sizeof(int),ofs);
  ofs.write((char*)&natom,sizeof(int));
  WriteSize(sizeof(int),ofs);

  delete [] molnames;
}

void WriteCSRCoords(ostream &ofs,OBMol &mol)
{
  int the_size,jconf;
  float x,y,z,energy;
  char title[100];
  char *tag;

  the_size = sizeof(int) + sizeof(float) + (80 * sizeof(char));
  
  jconf = 1;
  energy = -2.584565f;

  sprintf(title,"%s:%d",mol.GetTitle(),MolCount);
  tag = PadString(title,80);

  WriteSize(the_size,ofs);
  ofs.write((char*)&jconf,sizeof(int));
  ofs.write((char*)&energy,sizeof(float));
  ofs.write(tag,80*sizeof(char));
  WriteSize(the_size,ofs);

  WriteSize(mol.NumAtoms()*sizeof(float),ofs);

  OBAtom *atom;
  vector<OBNodeBase*>::iterator i;
  for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
  {
    x = atom->x();
    ofs.write((char*)&x,sizeof(float));
  }
  WriteSize(mol.NumAtoms()*sizeof(float),ofs);

  WriteSize(mol.NumAtoms()*sizeof(float),ofs);
  for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
  {
    y = atom->y();
    ofs.write((char*)&y,sizeof(float));
  }
  WriteSize(mol.NumAtoms()*sizeof(float),ofs);

  WriteSize(mol.NumAtoms()*sizeof(float),ofs);
  for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
  {
    z = atom->z();
    ofs.write((char*)&z,sizeof(float));
  }
  WriteSize(mol.NumAtoms()*sizeof(float),ofs);

  delete [] tag;
}

void WriteSize(int size,ostream &ofs)
{
  ofs.write((char*)&size,sizeof(int));
}

char *PadString(char *input, int size)
{
  unsigned int i;
  char *output;
  
  output = new char [size];
  for (i = 0; i < (unsigned)size; i++)
    output[i] = ' ';
  for (i = 0; i < strlen(input); i++)
    output[i] = input[i];
  return(output);
}

}

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

using namespace std;

namespace OpenBabel
{

static bool WriteChem3d(ostream &ofs,OBMol &mol, char *mol_typ);

bool WriteChem3d2(ostream &ofs,OBMol &mol)
{
  return(WriteChem3d(ofs,mol,"C3D"));
}

bool WriteChem3d1(ostream &ofs,OBMol &mol)
{
  return(WriteChem3d(ofs,mol,"MM2"));
}

bool WriteMmads(ostream &ofs,OBMol &mol)
{
  return(WriteChem3d(ofs,mol,"MMADS"));
}

static bool WriteChem3d(ostream &ofs,OBMol &mol, char *mol_typ)
{ 
  int atnum;
  int type_num;
  char buffer[BUFF_SIZE],type_name[10],ele_type[10];

  sprintf(buffer,"%d",mol.NumAtoms()); ofs << buffer;
  if (EQ(mol_typ,"MMADS"))
  {
    sprintf(buffer," %s",mol.GetTitle()); ofs << buffer;
    ttab.SetToType("MM2");
  }
  else ttab.SetToType(mol_typ);
  ofs << endl;

  ttab.SetFromType("INT"); 

  OBAtom *atom,*nbr;
  vector<OBNodeBase*>::iterator i;
  vector<OBEdgeBase*>::iterator j;

  for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
  {
    if (!ttab.Translate(type_name,atom->GetType()))
    {
      sprintf(buffer,"Unable to assign %s type to atom %d type = %s\n",
	      mol_typ,atom->GetIdx(),atom->GetType());
      ThrowError(buffer);
      atnum = atom->GetAtomicNum();
      type_num = atnum * 10 + atom->GetValence();
      sprintf(type_name,"%d",type_num);
    }
    strcpy(ele_type,etab.GetSymbol(atom->GetAtomicNum()));
    sprintf(buffer,"%-3s %-5d %8.4f  %8.4f  %8.4f %5s",
	    ele_type,
	    atom->GetIdx(),
	    atom->x(),
	    atom->y(),
	    atom->z(),
	    type_name); ofs << buffer;

    for (nbr = atom->BeginNbrAtom(j);nbr;nbr = atom->NextNbrAtom(j))
      {
	sprintf(buffer,"%6d",nbr->GetIdx()); 
	ofs << buffer;
      }
    ofs << endl;
  }
  return(true);
}

static bool ReadChem3d(istream&,OBMol&,bool,char*);

bool ReadMmads(istream &ifs,OBMol &mol,char *title)
{
  return(ReadChem3d(ifs,mol,true,"MM2"));
}

bool ReadChem3d1(istream &ifs,OBMol &mol,char *title)
{
  return(ReadChem3d(ifs,mol,false,"MM2"));
}

bool ReadChem3d2(istream &ifs,OBMol &mol,char *title)
{
  return(ReadChem3d(ifs,mol,false,"C3D"));
}

bool ReadChem3d(istream &ifs,OBMol &mol,bool mmads,char *type_key)
{
  char buffer[BUFF_SIZE];
  int natoms,i;
  char tmp[10],tmp1[10];
  char atomic_type[10];
  float exponent = 0.0f;
  float divisor = 1.0f;
  float Alpha,Beta,Gamma,A,B,C;
  bool has_fractional = false, has_divisor = false;
  Matrix3x3 m;
  
  vector<string> vs;
  ifs.getline(buffer,BUFF_SIZE);
  tokenize(vs,buffer);

  if (mmads)
  {
    if (vs.empty()) return(false);
    natoms = atoi((char*)vs[0].c_str());
    if (vs.size() == 2) mol.SetTitle(vs[1]);
  }
  else
    {
    switch(vs.size())
      {
      case 7 :
	sscanf(buffer,"%d%f%f%f%f%f%f",
	       &natoms,&Alpha,&Beta,&Gamma,&A,&B,&C);
	m.FillOrth(Alpha,Beta,Gamma,A,B,C);
	has_fractional = true;
	break;
      case 8 :
	sscanf(buffer,"%d%f%f%f%f%f%f%f",
	       &natoms,&Alpha,&Beta,&Gamma,&A,&B,&C,&exponent);
	m.FillOrth(Alpha,Beta,Gamma,A,B,C);
	has_fractional = true;
	has_divisor = true;
	break;   
      default :
	sscanf(buffer,"%d",&natoms);
	break;
      }
    }

  if (!natoms) return(false);
  divisor = pow(10.0,exponent);
  mol.ReserveAtoms(natoms);

  ttab.SetToType("INT");
  ttab.SetFromType(type_key);

  OBAtom *atom;
  float x,y,z;
  Vector v;

  unsigned int k;
  for (i = 1; i <= natoms; i++)
  {
    ifs.getline(buffer,BUFF_SIZE);
    sscanf(buffer,"%s%*d%f%f%f%s",
	   atomic_type,
	   &x,
	   &y,
	   &z,
	   tmp);
    v.Set(x,y,z);
    if (has_fractional) v *= m;
    if (has_divisor) v/= divisor;
      
    tokenize(vs,buffer);
    if (vs.empty()) return(false);

    atom = mol.NewAtom();
    ttab.Translate(tmp1,tmp);
    atom->SetType(tmp1);
    atom->SetVector(v);
    atom->SetAtomicNum(etab.GetAtomicNum(atomic_type));

    for (k = 6;k < vs.size(); k++)
      mol.AddBond(atom->GetIdx(),atoi((char*)vs[k].c_str()),1);
  }
  
  //assign_bond_order(mol);
  
  return(true);
}

}

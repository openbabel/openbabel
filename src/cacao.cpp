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

bool ReadCaccrt(istream &ifs,OBMol &mol,char *title)
{
  char buffer[BUFF_SIZE];
  int natoms;
  float A,B,C,Alpha,Beta,Gamma;
  Matrix3x3 m;
  
  ifs.getline(buffer,BUFF_SIZE); mol.SetTitle(buffer);
  ifs.getline(buffer,BUFF_SIZE); sscanf(buffer,"%d",&natoms);
  
  while (ifs.getline(buffer,BUFF_SIZE) &&!EQn(buffer,"CELL",4));

  if (!EQn(buffer,"CELL",4)) return(false);
  vector<string> vs;
  tokenize(vs,buffer," \n\t,");
  if (vs.size() != 7) return(false);

  //parse cell values
  A = atof((char*)vs[1].c_str());
  B = atof((char*)vs[2].c_str());
  C = atof((char*)vs[3].c_str());
  Alpha = atof((char*)vs[4].c_str());
  Beta  = atof((char*)vs[5].c_str());
  Gamma = atof((char*)vs[6].c_str());
  
  m.FillOrth(Alpha,Beta,Gamma,A,B,C);

  int i;
  float x,y,z;
  char type[10];
  OBAtom *atom;
  Vector v;

  for (i = 1; i <= natoms;i++)
    {
      if (!ifs.getline(buffer,BUFF_SIZE)) return(false);
      tokenize(vs,buffer," \n\t,");
      if (vs.size() < 4) return(false);
      atom = mol.NewAtom();

      x = atof((char*)vs[1].c_str());
      y = atof((char*)vs[2].c_str());
      z = atof((char*)vs[3].c_str());
      v.Set(x,y,z); v *= m;

      strcpy(type,vs[0].c_str());
      atom->SetAtomicNum(etab.GetAtomicNum(type));
      atom->SetVector(v);
    }

  //result = assign_radii(mol);
  //result = assign_bonds(mol);
  //result = assign_types(mol);
  //assign_bond_order(mol);

  return(true);
}


bool WriteCaccrt(ostream &ofs,OBMol &mol)
{
  OBAtom *atom;
  char type_name[10],buffer[BUFF_SIZE];
  vector<OBNodeBase*>::iterator i;

  sprintf(buffer,"%s\n",mol.GetTitle());
  sprintf(buffer,"%3d   DIST  0  0  0\n",mol.NumAtoms());
  sprintf(buffer,"CELL 1.,1.,1.,90.,90.,90.\n");

  for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
  {
    sprintf(buffer,"%2s %7.4f, %7.4f, %7.4f \n",
            etab.GetSymbol(atom->GetAtomicNum()),
            atom->x(),
            atom->y(),
            atom->z());
  }

  return(true);
}

static void SetHilderbrandt(OBMol&,vector<OBInternalCoord*>&);

bool WriteCacaoInternal(ostream &ofs,OBMol &mol)
{
  unsigned int i;
  Vector v;
  char tmptype[10],buffer[BUFF_SIZE];

  if (mol.Empty()) return(false);

  //translate first atom to origin
  v = mol.GetAtom(1)->GetVector(); v *= -1.0f; mol.Translate(v);

  vector<OBInternalCoord*> vit;
  SetHilderbrandt(mol,vit);
  strcpy(tmptype,etab.GetSymbol(mol.GetAtom(1)->GetAtomicNum()));

  sprintf(buffer," # TITLE"); ofs << buffer << endl;
  sprintf(buffer,"%3d  0DIST  0  0  0",mol.NumAtoms()); ofs << buffer << endl;
  sprintf(buffer,"  EL"); ofs << buffer << endl;
  sprintf(buffer,"0.,0.,0., %s",tmptype); ofs << buffer << endl;
  for (i = 2; i <= mol.NumAtoms(); i++)
  {
    strcpy(tmptype,etab.GetSymbol(mol.GetAtom(i)->GetAtomicNum()));

    if (vit[i]->_tor < 0.0) vit[i]->_tor += 360.0;
    sprintf(buffer,"%2d,%d,%2s%7.3f,%7.3f,%7.3f",
	    vit[i]->_a->GetIdx(),i,tmptype,
	    vit[i]->_dst,
	    vit[i]->_ang,
	    vit[i]->_tor);
    ofs << buffer << endl;
  }

  vector<OBInternalCoord*>::iterator j;
  for (j = vit.begin();j != vit.end();j++) 
    if (*j)
      delete *j;

  return(true);
}

void SetHilderbrandt(OBMol &mol,vector<OBInternalCoord*> &vit)
{
  float sum,r;

  OBAtom dummy1,dummy2;
  dummy1.SetVector(0.0f,0.0f,1.0f);
  dummy2.SetVector(1.0f,0.0f,0.0f);

  OBAtom *atom,*a1,*a2,*ref;
  vector<OBNodeBase*>::iterator ai;
  vit.push_back((OBInternalCoord*)NULL);
  for (atom = mol.BeginAtom(ai);atom;atom = mol.NextAtom(ai))
    vit.push_back(new OBInternalCoord (atom));

  vit[1]->_a = &dummy1;
  vit[1]->_b = &dummy2;
  vit[2]->_b = &dummy1;
  vit[2]->_c = &dummy2;
  vit[3]->_c = &dummy1;

  unsigned int i,j;
  for (i = 2;i <= mol.NumAtoms();i++)
    {
      ref = (OBAtom*)NULL;
      a1 = mol.GetAtom(i);
      sum = 100.0f;
      for (j = 1;j < i;j++)
	{
	  a2 = mol.GetAtom(j);
	  r = (a1->GetVector()-a2->GetVector()).length_2();
	  if ((r < sum) && (vit[j]->_a != a2) && (vit[j]->_b != a2))
	    {
	      sum = r;
	      ref = a2;
	    }
	}
      vit[i]->_a = ref;
    }

  for (i = 3;i <= mol.NumAtoms();i++)
    vit[i]->_b = vit[vit[i]->_a->GetIdx()]->_a;

  for (i = 4;i <= mol.NumAtoms();i++)
    {
      if (vit[i]->_b->GetIdx())
	vit[i]->_c = vit[vit[i]->_b->GetIdx()]->_b;
      else
	vit[i]->_c = &dummy1;
    }

  OBAtom *a,*b,*c;
  Vector v1,v2;
  for (i = 2;i <= mol.NumAtoms();i++)
    {
      atom = mol.GetAtom(i);
      a = vit[i]->_a; b = vit[i]->_b; c = vit[i]->_c;
      v1 = atom->GetVector() - a->GetVector();
      v2 = b->GetVector() - a->GetVector();
      vit[i]->_ang = VectorAngle(v1,v2);
      vit[i]->_tor = CalcTorsionAngle(atom->GetVector(),
				      a->GetVector(),
				      b->GetVector(),
				      c->GetVector());
      vit[i]->_dst = (vit[i]->_a->GetVector() - atom->GetVector()).length();
    }

}

}

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

#include <list>
#include "mol.h"
#include "obutil.h"

using namespace OpenBabel;

#include "matrix.h"
#include "chiral.h"

void OBMol::FindChiralCenters()
{
  if (HasChiralityPerceived()) return;
  SetChiralityPerceived();

  //do quick test to see if there are any possible chiral centers
  bool mayHaveChiralCenter=false;
  OBAtom *atom,*nbr;
  vector<OBNodeBase*>::iterator i;
  for (atom = BeginAtom(i);atom;atom = NextAtom(i))
    if (atom->GetHyb() == 3 && atom->GetHvyValence() >= 3)
      {
	mayHaveChiralCenter=true;
	break;
      }

  if (!mayHaveChiralCenter) return;

  OBBond *bond;
  vector<OBEdgeBase*>::iterator j;
  for (bond = BeginBond(j);bond;bond = NextBond(j))
    if (bond->IsWedge() || bond->IsHash())
	(bond->GetBeginAtom())->SetChiral();

#define INTMETHOD

#ifdef INTMETHOD
	vector<unsigned int> vgid;
	GetGIDVector(vgid);
	vector<unsigned int> tlist;
	vector<unsigned int>::iterator k;
#else //use Golender floating point method
	vector<float> gp;  
	GraphPotentials(*this,gp);
	vector<float> tlist;
	vector<float>::iterator k;
#endif

  bool ischiral;
  for (atom = BeginAtom(i);atom;atom = NextAtom(i))
    if (atom->GetHyb() == 3 && atom->GetHvyValence() >= 3 && !atom->IsChiral())
      {
		tlist.clear();
		ischiral = true;

		for (nbr = atom->BeginNbrAtom(j);nbr;nbr = atom->NextNbrAtom(j))
		{
			for (k = tlist.begin();k != tlist.end();k++)
#ifdef INTMETHOD
				if (vgid[nbr->GetIdx()-1] == *k)
#else
				if (fabs(gp[nbr->GetIdx()]-*k) < 0.001)
#endif
					ischiral = false;

#ifdef INTMETHOD
			if (ischiral) tlist.push_back(vgid[nbr->GetIdx()-1]);
#else
			if (ischiral) tlist.push_back(gp[nbr->GetIdx()]);
#endif
			else break;
		}

	if (ischiral) atom->SetChiral();
      }
}

void GetChirality(OBMol &mol, vector<int> &chirality)
{
  chirality.resize(mol.NumAtoms()+1);
  fill(chirality.begin(),chirality.end(),0);

  OBAtom *atom;
  vector<OBNodeBase*>::iterator i;
  for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
    if (atom->IsChiral())
      {
	float sv = CalcSignedVolume(mol,atom);
	if (sv < 0.0f)	chirality[atom->GetIdx()-1] = -1;
	else if (sv > 0.0)  chirality[atom->GetIdx()-1] = 1;
    }
}

// Calculate the signed volume for an atom.  If the atom has a valence of 3
// the coordinates of an attached hydrogen are calculated

float CalcSignedVolume(OBMol &mol,OBAtom *atm)
{
  Vector tmp_crd;
  vector<int> nbr_atms;
  vector<Vector> nbr_crds;
  float hbrad = etab.CorrectedBondRad(1,0);
	
  if (atm->GetHvyValence() < 3)
    {
      cerr << "Cannot calculate a signed volume for an atom with a heavy atom valence of " << atm->GetHvyValence() << endl;
      exit(0);
    }

  // Create a vector with the coordinates of the neighbor atoms
  OBAtom *nbr;
  vector<OBEdgeBase*>::iterator bint;
  for (nbr = atm->BeginNbrAtom(bint);nbr;nbr = atm->NextNbrAtom(bint))
    {
      nbr_atms.push_back(nbr->GetIdx());		
    }
  // sort the neighbor atoms to insure a consistent ordering
  sort(nbr_atms.begin(),nbr_atms.end());
  for (unsigned int i = 0; i < nbr_atms.size(); i++)
    {
      OBAtom *tmp_atm = mol.GetAtom(nbr_atms[i]);
      nbr_crds.push_back(tmp_atm->GetVector());
    }	

  // If we have three heavy atoms we need to calculate the position of the fourth
  if (atm->GetHvyValence() == 3)
    {		
      float bondlen = hbrad+etab.CorrectedBondRad(atm->GetAtomicNum(),atm->GetHyb());
      atm->GetNewBondVector(tmp_crd,bondlen);
      nbr_crds.push_back(tmp_crd);
    }

  return(signed_volume(nbr_crds[0],nbr_crds[1],nbr_crds[2],nbr_crds[3]));
}

// calculate a signed volume given a set of 4 coordinates
float signed_volume(const Vector &a, const Vector &b, const Vector &c, const Vector &d)
{
  Vector A,B,C;
  A = b-a;
  B = c-a;
  C = d-a;
  Matrix3x3 m(A,B,C);
  return(m.determinant());
}

// Calculate the Graph Potentials of a molecule
// based on 
// V.E. and Rozenblit, A.B. Golender
// Logical and Combinatorial Algorithms for Drug Design
// for an example see
// Walters, W. P., Yalkowsky, S. H., JCICS, 1996, 36(5), 1015-1017

void GraphPotentials(OBMol &mol, vector<float> &pot)
{
  float det;

  vector<vector<float> > g,c,h;
  construct_g_matrix(mol,g);
  invert_matrix(g,det); 
  construct_c_matrix(mol,c);
  mult_matrix(h,g,c);
  pot.resize(mol.NumAtoms()+1);

  for (unsigned int i = 0; i < mol.NumAtoms();i++) pot[i+1] = h[i][0];
}


// Construct the matrix G, which puts each atoms valence+1
// on the diagonal and and -1 on the off diagonal if two
// atoms are connected.

void construct_g_matrix(OBMol &mol, vector<vector<float> > &m)
{
  unsigned int i,j;

  OBAtom *atm1,*atm2;
  vector<OBNodeBase*>::iterator aint,bint;

  m.resize(mol.NumAtoms());
  for (i = 0; i < m.size(); i++)
    m[i].resize(mol.NumAtoms());

  for (atm1 = mol.BeginAtom(aint),i=0;atm1;atm1 = mol.NextAtom(aint),i++)
    for (atm2 = mol.BeginAtom(bint),j=0;atm2;atm2 = mol.NextAtom(bint),j++)
      {
	if (i == j)
	  {
	    m[i][j] = atm1->GetValence() + 1;
	    m[i][j] += (float)atm1->GetAtomicNum()/10.0;
	    m[i][j] += (float)atm1->GetHyb()/100.0;
	  }
	else
	  {
	    if (atm1->IsConnected(atm2))
	      m[i][j] = -1;
	    else
	      m[i][j] = 0;
	  }
      }
}

// Construct the matrix C, which is simply a column vector
// consisting of the valence for each atom
void construct_c_matrix(OBMol &mol,vector<vector<float > > &m)
{
  unsigned int i;
  OBAtom *atm1;
  vector<OBNodeBase*>::iterator aint;
  
  m.resize(mol.NumAtoms());
  for (i = 0; i < m.size(); i++)
    m[i].resize(1);
  for (atm1 = mol.BeginAtom(aint),i=0;atm1;atm1 = mol.NextAtom(aint),i++)
    {
      m[i][0] = atm1->GetValence();
    }
}





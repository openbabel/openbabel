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
#include "molchrg.h"

using namespace std;
namespace OpenBabel {

bool OBGastChrg::AssignPartialCharges(OBMol &mol)
{
  //InitialPartialCharges(mol);
  
  OBAtom *atom;
  vector<OBNodeBase*>::iterator i;

  GSVResize(mol.NumAtoms()+1);
  float a,b,c;
  for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
    {
      if (!GasteigerSigmaChi(atom,a,b,c)) return(false);
      _gsv[atom->GetIdx()]->SetValues(a,b,c,atom->GetPartialCharge());
    }

  float alpha,charge,denom;
  unsigned j;
  int iter;
  OBBond *bond;
  OBAtom *src,*dst;
  vector<OBEdgeBase*>::iterator k;
  alpha = 1.0;
  for(iter = 0;iter < MX_GASTEIGER_ITERS;iter++) 
    {
      alpha *= MX_GASTEIGER_DAMP;
      
      for( j=1;j < _gsv.size();j++)
	{
	  charge = _gsv[j]->q;
	  _gsv[j]->chi = (_gsv[j]->c*charge+_gsv[j]->b)*charge+_gsv[j]->a;
	}

      for (bond = mol.BeginBond(k);bond;bond = mol.NextBond(k))
	{
	  src = bond->GetBeginAtom();
	  dst = bond->GetEndAtom();

	  if (_gsv[src->GetIdx()]->chi >= _gsv[dst->GetIdx()]->chi)
	    {
	      if (dst->IsHydrogen()) denom = float(MX_GASTEIGER_DENOM);
	      else                   denom = _gsv[dst->GetIdx()]->denom;
	    }
	  else
	    {
	      if (src->IsHydrogen()) denom = float(MX_GASTEIGER_DENOM);
	      else                   denom = _gsv[src->GetIdx()]->denom;
	    }

	  charge = (_gsv[src->GetIdx()]->chi - _gsv[dst->GetIdx()]->chi)/denom;
	  _gsv[src->GetIdx()]->q -= alpha*charge;
	  _gsv[dst->GetIdx()]->q += alpha*charge;
	}
    }

  for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
      atom->SetPartialCharge(_gsv[atom->GetIdx()]->q);

  return(true);
}

void OBGastChrg::InitialPartialCharges(OBMol &mol)
{
  OBAtom *atom;
  vector<OBNodeBase*>::iterator i;
  
  for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
    {
      if (atom->IsCarboxylOxygen()) atom->SetPartialCharge(-0.500f);
      else if (atom->IsPhosphateOxygen() &&
	       atom->GetHvyValence() == 1) atom->SetPartialCharge(-0.666f);
      else if (atom->IsSulfateOxygen())   atom->SetPartialCharge(-0.500f);
      else atom->SetPartialCharge((float)atom->GetFormalCharge());
    }
}

bool OBGastChrg::GasteigerSigmaChi(OBAtom *atom,float &a,float &b,float &c )
{
  int count;
  float val[3] = {0.0,0.0,0.0};

  switch(atom->GetAtomicNum()) 
    {
    case 1: //H
      val[0] = 0.37f;val[1] = 7.17f;val[2] = 12.85f; 
      break;
    case 6: //C
      if (atom->GetHyb() == 3)	{val[0] = 0.68f;val[1] = 7.98f;val[2] = 19.04f;}
      if (atom->GetHyb() == 2)	{val[0] = 0.98f;val[1] = 8.79f;val[2] = 19.62f;}
      if (atom->GetHyb() == 1)	{val[0] = 1.67f;val[1] = 10.39f;val[2] = 20.57f;}
      break;
    case 7: //N
      if (atom->GetHyb() == 3)
	{
	  if (atom->GetValence() == 4 || atom->GetFormalCharge())
	    {val[0] = 0.0f;val[1] = 0.0f;val[2] = 23.72f;}
	  else
	    {val[0] = 2.08f;val[1] = 11.54f;val[2] = 23.72f;}
	}

      if (atom->GetHyb() == 2)
	{
	  if (EQ(atom->GetType(),"Npl") || EQ(atom->GetType(),"Nam"))
	    {val[0] = 2.46f;val[1] = 12.32f;val[2] = 24.86f;}
	  else
	    {val[0] = 2.57f;val[1] = 12.87f;val[2] = 24.87f;}
	}

      if (atom->GetHyb() == 1)	{val[0] = 3.71f;val[1] = 15.68f;val[2] = 27.11f;}
      break;
    case 8: //O
      if (atom->GetHyb() == 3) {val[0] = 2.65f;val[1] = 14.18f;val[2] = 28.49f;}
      if (atom->GetHyb() == 2) {val[0] = 3.75f;val[1] = 17.07f;val[2] = 31.33f;}
      break;
    case 9: //F
      val[0] = 3.12f;val[1] = 14.66f;val[2] = 30.82f;
      break;
    case 15: //P
      val[0] = 1.62f;val[1] = 8.90f;val[2] = 18.10f;
      break;
    case 16: //S
      count = atom->CountFreeOxygens();
      if (count == 0 || count == 1) 
	{val[0] = 2.39f;val[1] = 10.14f;val[2] = 20.65f;}
      if (count > 1) {val[0] = 2.39f;val[1] = 12.00f;val[2] = 24.00f;}
      /*S2? if (count == 0) {val[0] = 2.72;val[1] = 10.88;val[2] = 21.69;}*/
      break;
    case 17: //Cl
      val[0] = 2.66f;val[1] = 11.00f;val[2] = 22.04f;
      break;
    case 35: //Br
      val[0] = 2.77f;val[1] = 10.08f;val[2] = 19.71f;
      break;
    case 53: //I
      val[0] = 2.90f;val[1] = 9.90f;val[2] = 18.82f;
      break;
    case 13: //Al
      val[0] = 1.06f;val[1] = 5.47f;val[2] = 11.65f;
	break;
    }

  if (val[2] != 0.0)
    {
      a = val[1];
      b = (val[2]-val[0])/2;
      c = (val[2]+val[0])/2 - val[1];
    }
  else
    return(false);

  return(true);
}

void OBGastChrg::GSVResize(int size)
{
  vector <GasteigerState*>::iterator i;
  for (i = _gsv.begin();i != _gsv.end();i++) delete *i;
  _gsv.clear();

  for (int j = 0;j < size;j++) _gsv.push_back(new GasteigerState);
}

OBGastChrg::~OBGastChrg()
{
  vector <GasteigerState*>::iterator i;
  for (i = _gsv.begin();i != _gsv.end();i++) delete *i;
}

GasteigerState::GasteigerState()
{
  a = 0.0;
  b = 0.0;
  c = 0.0;
  denom = 0.0;
  chi = 0.0;
  q = 0.0;
}

}

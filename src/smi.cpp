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
#include "typer.h"
#include "smi.h"

//#define KEKULE

using namespace std;

namespace OpenBabel {

bool WriteSmiles(ostream &ofs,OBMol &mol,const char *title)
{
  char buffer[BUFF_SIZE],tmp[BUFF_SIZE];

  OBMol2Smi m2s;

  m2s.Init();
  m2s.CorrectAromaticAmineCharge(mol);
  m2s.CreateSmiString(mol,buffer);

  strcpy(tmp,(title) ? title:mol.GetTitle());
  ofs << buffer << ' ' << tmp << endl;
  return(true);
}

void OBMol2Smi::CreateSmiString(OBMol &mol,char *buffer)
{
  OBAtom *atom;
  OBSmiNode *root;
  buffer[0] = '\0';
  vector<OBNodeBase*>::iterator i;

  for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
//    if ((!atom->IsHydrogen() || atom->GetValence() == 0) && !_uatoms[atom->GetIdx()])
    if (!atom->IsHydrogen() && !_uatoms[atom->GetIdx()])
      if (!atom->IsChiral()) //don't use chiral atoms as root node
	{
	  //clear out closures in case structure is dot disconnected
	  _vclose.clear();
	  _atmorder.clear();
	  _storder.clear();
	  _vopen.clear();
	  //dot disconnected structure
	  if (strlen(buffer) > 0) strcat(buffer,"."); 
	  root = new OBSmiNode (atom);
	  BuildTree(root);
	  FindClosureBonds(mol);
	  if (mol.Has2D()) AssignCisTrans(root);
	  ToSmilesString(root,buffer);
	  delete root;
	}
}

bool OBMol2Smi::BuildTree(OBSmiNode *node)
{
  vector<OBEdgeBase*>::iterator i;
  OBAtom *nbr,*atom = node->GetAtom();
  
  _uatoms.SetBitOn(atom->GetIdx()); //mark the atom as visited
  _atmorder.push_back(atom->GetIdx()); //store the atom ordering
  _storder.push_back(atom->GetIdx()); //store the atom ordering for stereo

  for (nbr = atom->BeginNbrAtom(i);nbr;nbr = atom->NextNbrAtom(i))
    if (!nbr->IsHydrogen() && !_uatoms[nbr->GetIdx()])
      {
			_ubonds.SetBitOn((*i)->GetIdx());
			OBSmiNode *next = new OBSmiNode (nbr);
			next->SetParent(atom);
			node->SetNextNode(next,(OBBond*)*i);
			BuildTree(next);
	}

  return(true);
}

void OBMol2Smi::ToSmilesString(OBSmiNode *node,char *buffer)
{
  char tmpbuf[10];
  OBAtom *atom = node->GetAtom();

  //write the current atom to the string
  GetSmilesElement(node,tmpbuf);
  strcat(buffer,tmpbuf);

    //handle ring closures here
  vector<pair<int,OBBond*> > vc = GetClosureDigits(atom);
  if (!vc.empty())
    {
      vector<pair<int,OBBond*> >::iterator i;
      for (i = vc.begin();i != vc.end();i++)
	{
	  if (i->second)
	    {
	      if (i->second->IsUp())                              strcat(buffer,"\\");
	      if (i->second->IsDown())                            strcat(buffer,"/");
#ifndef KEKULE
	      if (i->second->GetBO() == 2 && !i->second->IsAromatic()) strcat(buffer,"=");
#else
	      if (i->second->GetBO() == 2)                        strcat(buffer,"=");
#endif
	      if (i->second->GetBO() == 3)                        strcat(buffer,"#");
	    }
	      
	  if (i->first > 9) strcat(buffer,"%");
	  sprintf(tmpbuf,"%d",i->first); 
	  strcat(buffer,tmpbuf);
	}
    }

  //follow path to child atoms
  OBBond *bond;
  for (int i = 0;i < node->Size();i++)
    {
      bond = node->GetNextBond(i);
      if (i+1 < node->Size())                        strcat(buffer,"(");
      if (bond->IsUp())                              strcat(buffer,"\\");
      if (bond->IsDown())                            strcat(buffer,"/");
#ifndef KEKULE
      if (bond->GetBO() == 2 && !bond->IsAromatic()) strcat(buffer,"=");
#else
      if (bond->GetBO() == 2)                        strcat(buffer,"=");
#endif
      if (bond->GetBO() == 3)                        strcat(buffer,"#");

      ToSmilesString(node->GetNextNode(i),buffer);
      if (i+1 < node->Size()) strcat(buffer,")");
    }
}

void OBMol2Smi::GetClosureAtoms(OBAtom *atom,vector<OBNodeBase*> &va)
{

//look through closure list for start atom
  vector<OBEdgeBase*>::iterator i;
  for (i = _vclose.begin();i != _vclose.end();i++)
	  if (*i)
  {
	  if (((OBBond*)*i)->GetBeginAtom() == atom)
	    va.push_back(((OBBond*)*i)->GetEndAtom());
	  if (((OBBond*)*i)->GetEndAtom() == atom)
	    va.push_back(((OBBond*)*i)->GetBeginAtom());
  }

  OBAtom *nbr;
  vector<pair<OBAtom*,pair<int,int> > >::iterator j;
  for (j = _vopen.begin();j != _vopen.end();j++)
	  for (nbr = atom->BeginNbrAtom(i);nbr;nbr = atom->NextNbrAtom(i))
		  if (nbr == j->first)
			va.push_back(nbr);  
}

vector<pair<int,OBBond*> > OBMol2Smi::GetClosureDigits(OBAtom *atom)
{
  vector<pair<int,OBBond*> > vc; vc.clear();

  //look through closure list for start atom
  int idx,bo;
  OBBond *bond;
  vector<OBEdgeBase*>::iterator i;
  for (i = _vclose.begin();i != _vclose.end();i++)
    if ((bond=(OBBond*)*i))
      if (bond->GetBeginAtom() == atom || bond->GetEndAtom() == atom)
	{
	  idx = GetUnusedIndex();
	  vc.push_back(pair<int,OBBond*> (idx,bond));
	  bo = (bond->IsAromatic())? 1 : bond->GetBO();
	  _vopen.push_back(pair<OBAtom*,pair<int,int> >
			   (bond->GetNbrAtom(atom),pair<int,int>(idx,bo)));
	  *i = NULL;//remove bond from closure list
	}

  //try to complete closures
  if (!_vopen.empty())
    {
      vector<pair<OBAtom*,pair<int,int> > >::iterator j;
      for (j = _vopen.begin();j != _vopen.end();)
	if (j->first == atom)
	  {
	    vc.push_back(pair<int,OBBond*> (j->second.first,(OBBond*)NULL));
	    _vopen.erase(j);
	    j = _vopen.begin();
	  }
	else j++;
    }

  return(vc);
}

void OBMol2Smi::FindClosureBonds(OBMol &mol)
{
  //find closure bonds
  OBAtom *a1,*a2;
  OBBond *bond;
  vector<OBEdgeBase*>::iterator i;
  OBBitVec bv;
  bv.FromVecInt(_storder);

  for (bond = mol.BeginBond(i);bond;bond = mol.NextBond(i))
    if (!_ubonds[bond->GetIdx()] && bv[bond->GetBeginAtomIdx()])
      {
	a1 = bond->GetBeginAtom();
	a2 = bond->GetEndAtom();
	if (!a1->IsHydrogen() && !a2->IsHydrogen()) 
	  _vclose.push_back(bond);
      }

  vector<OBEdgeBase*>::reverse_iterator j;
  vector<int>::iterator k;

  //modify _order to reflect ring closures
  for (j = _vclose.rbegin();j != _vclose.rend();j++)
    {
      bond = (OBBond*)*j;
      a1 = a2 = NULL;

      for (k = _storder.begin();k != _storder.end();k++)
	if (bond->GetBeginAtomIdx() == static_cast<unsigned int>(*k) || 
	    bond->GetEndAtomIdx() == static_cast<unsigned int>(*k))
	  if (!a1) a1 = mol.GetAtom(*k);
	  else if (!a2) 
	    {
	      a2 = mol.GetAtom(*k);
	      _storder.erase(k);
	      break;
	    }

      for (k = _storder.begin();k != _storder.end();k++)
		if (a1->GetIdx() == static_cast<unsigned int>(*k))
		{
			k++;
			if (k != _storder.end()) _storder.insert(k,a2->GetIdx());
			else                     _storder.push_back(a2->GetIdx());
			break;
		}
    }
}

int OBMol2Smi::GetUnusedIndex()
{
  int idx=1;

  vector<pair<OBAtom*,pair<int,int> > >::iterator j;
  for (j = _vopen.begin();j != _vopen.end();)
    if (j->second.first == idx)
      {
	idx++; //increment idx and start over if digit is already used
	j = _vopen.begin();
      }
    else j++;

  return(idx);
}

void OBMol2Smi::CorrectAromaticAmineCharge(OBMol &mol)
{
  OBAtom *atom;
  vector<OBNodeBase*>::iterator i;

  _aromNH.clear();
  _aromNH.resize(mol.NumAtoms()+1);

  for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
    if (atom->IsNitrogen() && atom->IsAromatic())
      if (atom->GetHvyValence() == 2)
	{
	  if (atom->GetValence() == 3 || atom->GetImplicitValence() == 3)
	    _aromNH[atom->GetIdx()] = true;
	}
}

void OBMol2Smi::AssignCisTrans(OBSmiNode *node)
{
  //traverse the tree searching for acyclic olefins - if it
  //has at least one heavy atom attachment on each end assign stereochem

  OBBond *bond;
  for (int i = 0;i < node->Size();i++)
    {
      bond = node->GetNextBond(i);
      if (bond->GetBO() == 2 && !bond->IsInRing())
	{
	  OBAtom *b = node->GetAtom();
	  OBAtom *c = bond->GetNbrAtom(b);

	//skip allenes
	  if (b->GetHyb() == 1 || c->GetHyb() == 1) continue;

	  if (b->GetHvyValence() > 1 && c->GetHvyValence() > 1)
	    {
	      OBAtom *a,*d;
	      vector<OBEdgeBase*>::iterator j,k;

	      //look for bond with assigned stereo as in poly-ene
	      for (a = b->BeginNbrAtom(j);a;a = b->NextNbrAtom(j))
		if (((OBBond*)*j)->IsUp() ||((OBBond*)*j)->IsDown())
		  break;

	      if (!a)
	      for (a = b->BeginNbrAtom(j);a;a = b->NextNbrAtom(j))
		if (a != c && !a->IsHydrogen())
		  break;
	      for (d = c->BeginNbrAtom(k);d;d = c->NextNbrAtom(k))
		if (d != b && !d->IsHydrogen())
		  break;
	      obAssert(a); obAssert(d);
	      
	      if (((OBBond*)*j)->IsUp() || ((OBBond*)*j)->IsDown()) //stereo already assigned
		{
		  if (fabs(CalcTorsionAngle(a->GetVector(),b->GetVector(),
				       c->GetVector(),d->GetVector())) > 10.0)
		    if (((OBBond*)*j)->IsUp()) ((OBBond*)*k)->SetUp();
		    else              ((OBBond*)*k)->SetDown();
		  else
		    if (((OBBond*)*j)->IsUp()) ((OBBond*)*k)->SetDown();
		    else              ((OBBond*)*k)->SetUp();
		}
	      else //assign stereo to both ends
		{
		  ((OBBond*)*j)->SetUp();
		  if (fabs(CalcTorsionAngle(a->GetVector(),b->GetVector(),
				       c->GetVector(),d->GetVector())) > 10.0)
		    ((OBBond*)*k)->SetUp();
		  else
		    ((OBBond*)*k)->SetDown();
		}
	    }
	}
      AssignCisTrans(node->GetNextNode(i));
    }
}

void OBMol2Smi::Init()
  {
  _vclose.clear();
  _atmorder.clear();
  _storder.clear();
  _aromNH.clear(); 
  _uatoms.Clear();
  _ubonds.Clear();
  _vopen.clear();
}

bool OBMol2Smi::GetSmilesElement(OBSmiNode *node,char *element)
{
  //***handle reference atom stuff here and return***
  char symbol[10];
  bool bracketElement = false;
  bool normalValence = true;

  OBAtom *atom = node->GetAtom();

  int bosum = atom->KBOSum();

  switch (atom->GetAtomicNum())
    {
    case 0: break;
    case 5: /*bracketElement = !(normalValence = (bosum == 3)); break;*/ break;
    case 6: break;
    case 7: 
				if (atom->IsAromatic() && atom->GetHvyValence() == 2 && atom->GetImplicitValence() == 3)
				{
					bracketElement = !(normalValence = false); break;
				}
				else
					bracketElement = !(normalValence = (bosum == 3 || bosum == 5)); break;
    case 8: break;
    case 9: break;
    case 15: break;
    case 16: 
			bracketElement = !(normalValence = (bosum == 2 || bosum == 4 || bosum == 6)); break;
    case 17: break;
    case 35: break;
    case 53: break;
    
    default:
      bracketElement = true;
    }

  if (atom->GetHvyValence() > 2 && atom->IsChiral())
    if (((OBMol*)atom->GetParent())->HasNonZeroCoords() || atom->HasChiralitySpecified())
      bracketElement = true;

  if (atom->GetFormalCharge() != 0) //bracket charged elements
    bracketElement = true;

  if (!bracketElement)
    {
      if (!atom->GetAtomicNum())
	{
	  bool external = false;
	  vector<pair<int,pair<OBAtom *,OBBond *> > > *externalBonds = 
			(vector<pair<int,pair<OBAtom *,OBBond *> > > *)((OBMol*)atom->GetParent())->GetData("extBonds");
	  vector<pair<int,pair<OBAtom *,OBBond *> > >::iterator externalBond;

	  if (externalBonds)
	    for(externalBond = externalBonds->begin();externalBond != externalBonds->end();externalBond++)
	      {
		if (externalBond->second.first == atom)
		  {
		    external = true;
		    strcpy(symbol,"&");
		    OBBond *bond = externalBond->second.second;
		    if (bond->IsUp())                              strcat(symbol,"\\");
		    if (bond->IsDown())                            strcat(symbol,"/");
#ifndef KEKULE
		    if (bond->GetBO() == 2 && !bond->IsAromatic()) strcat(symbol,"=");
		    if (bond->GetBO() == 2 && bond->IsAromatic())  strcat(symbol,";");
#else
		    if (bond->GetBO() == 2)                        strcat(symbol,"=");
#endif
		    if (bond->GetBO() == 3)                        strcat(symbol,"#");
		    sprintf(symbol,"%s%d",symbol,externalBond->first);
		    break;
		  }
	      }

	  if(!external) strcpy(symbol,"*");
	}
      else
	{
	  strcpy(symbol,etab.GetSymbol(atom->GetAtomicNum()));
#ifndef KEKULE
	  if (atom->IsAromatic()) symbol[0] = tolower(symbol[0]);
#endif
	}
      strcpy(element,symbol);

      return(true);
    }

  strcpy(element,"[");
  if (!atom->GetAtomicNum()) strcpy(symbol,"*");
  else
    {
      strcpy(symbol,etab.GetSymbol(atom->GetAtomicNum()));
#ifndef KEKULE
      if (atom->IsAromatic()) symbol[0] = tolower(symbol[0]);
#endif
    }
  strcat(element,symbol);
  
	  //if (atom->IsCarbon() && atom->GetHvyValence() > 2 && atom->IsChiral())
  if (atom->GetHvyValence() > 2 && atom->IsChiral())
    {
      char stereo[5];
      if (GetChiralStereo(node,stereo))
		  strcat(element,stereo);
    }

  //add extra hydrogens
//  if (!normalValence && atom->ImplicitHydrogenCount())
  if (atom->ImplicitHydrogenCount())
    {
      strcat(element,"H");
      if (atom->ImplicitHydrogenCount() > 1)
	{
	  char tcount[10];
	  sprintf(tcount,"%d",atom->ImplicitHydrogenCount());
	  strcat(element,tcount);
	}
    }

  //cat charge on the end
  if (atom->GetFormalCharge() != 0)
    {

/*
		if (atom->ImplicitHydrogenCount())
		{
			cerr << "imp = " << atom->GetAtomicNum() << ' ' << atom->GetImplicitValence() << endl;
      strcat(element,"H");
      if (atom->ImplicitHydrogenCount() > 1)
			{
				char tcount[10];
				sprintf(tcount,"%d",atom->ImplicitHydrogenCount());
				strcat(element,tcount);
			}
		}
	*/	
			if (atom->GetFormalCharge() > 0) strcat(element,"+");
			else                             strcat(element,"-");

      if (abs(atom->GetFormalCharge()) > 1)
			{
				char tcharge[10];
				sprintf(tcharge,"%d",abs(atom->GetFormalCharge()));
				strcat(element,tcharge);
			}
    }

  strcat(element,"]");

  return(true);
}

bool OBMol2Smi::GetChiralStereo(OBSmiNode *node,char *stereo)
{
  bool is2D=false;
  float torsion;
  OBAtom *a,*b,*c,*d,hydrogen;

  b = node->GetAtom();
  OBMol *mol = (OBMol*)b->GetParent();

  if (!mol->HasNonZeroCoords()) //must have come in from smiles string
    {
      if (!b->HasChiralitySpecified()) return(false);
      if (b->IsClockwise())          strcpy(stereo,"@@");
      else if (b->IsAntiClockwise()) strcpy(stereo,"@");
      else return(false);
      //if (b->GetHvyValence() == 3) strcat(stereo,"H");
      return(true);
    }

  //give peudo Z coords if mol is 2D
  if (!mol->Has3D())
    {
      vector3 v,vz(0.0,0.0,1.0);
      is2D = true;
	  OBAtom *nbr;
      OBBond *bond;
      vector<OBEdgeBase*>::iterator i;
      for (bond = b->BeginBond(i);bond;bond = b->NextBond(i))
	  {
		  nbr = bond->GetEndAtom();
		  if (nbr != b)
		  {
			  v = nbr->GetVector();
			  if (bond->IsWedge()) v += vz;
			  else
				  if (bond->IsHash()) v -= vz;

			  nbr->SetVector(v);
		  }
		  else
		  {
			  nbr = bond->GetBeginAtom();
			  v = nbr->GetVector();
			  if (bond->IsWedge()) v -= vz;
			  else
				  if (bond->IsHash()) v += vz;

			  nbr->SetVector(v);
		  }
	}
    }

  c = d = NULL;
  a = node->GetParent();
  obAssert(a); //chiral atom can't be used as root node - must have parent
      
  if (b->GetHvyValence() == 3) //must have attached hydrogen
    {
      if (b->GetValence() == 4)//has explicit hydrogen
	{
	  vector<OBEdgeBase*>::iterator i;
	  for (c = b->BeginNbrAtom(i);c;c = b->NextNbrAtom(i))
	    if (c->IsHydrogen())
	      break;
	  obAssert(c);
	}
      else  //implicit hydrogen
	{
	  vector3 v;
	  b->GetNewBondVector(v,1.0);
	  hydrogen.SetVector(v);
	  c = &hydrogen;
	}
    }

  //get connected atoms in order
  OBAtom *nbr;
  vector<int>::iterator j;

  //try to get neighbors that are closure atoms in the order they appear in the string
  vector<OBNodeBase*> va;
  GetClosureAtoms(b,va);
  if (!va.empty())
  {
	  vector<OBNodeBase*>::iterator k;
	  for (k = va.begin();k != va.end();k++)
		  if (*k != a)
		  {
			if (!c) c = (OBAtom*)*k;
			else if (!d) d = (OBAtom*)*k;
		  }
  }

  for (j = _storder.begin();j != _storder.end();j++)
    {
      nbr = mol->GetAtom(*j);
      if (!b->IsConnected(nbr)) continue;
      if (nbr == a || nbr == b || nbr == c) continue;
      if (!c)      c = nbr;
      else if (!d) d = nbr;
    }

  torsion = CalcTorsionAngle(a->GetVector(),b->GetVector(), 
			     c->GetVector(),d->GetVector());

  strcpy(stereo,(torsion<0.0)?"@":"@@");
  //if (b->GetHvyValence() == 3) strcat(stereo,"H");

  //re-zero psuedo-coords
  if (is2D)
    {
      vector3 v;
      OBAtom *atom;
      vector<OBNodeBase*>::iterator k;
      for (atom = mol->BeginAtom(k);atom;atom = mol->NextAtom(k))
	  {
		  v = atom->GetVector();
		  v.SetZ(0.0);
		  atom->SetVector(v);
	  }
    }

  return(true);
}

bool WriteFixFile(ostream &ofs,OBMol &mol)
{
  char buffer[BUFF_SIZE];
  OBMol2Smi m2s;

  m2s.Init();
  //m2s.AssignCisTrans(mol);
  m2s.CorrectAromaticAmineCharge(mol);
  m2s.CreateSmiString(mol,buffer);

  OBAtom *atom;
  vector<int>::iterator i;
  vector<int> order = m2s.GetOutputOrder();
  ofs << buffer << endl;
  
  int j;
  for (j = 0;j < mol.NumConformers();j++)
    {
      mol.SetConformer(j);
      for (i = order.begin();i != order.end();i++)
	{
	  atom = mol.GetAtom(*i);
	  sprintf(buffer,"%9.3f %9.3f %9.3f",atom->GetX(),atom->GetY(),atom->GetZ());
	  ofs << buffer<< endl;
	}
    }
  return(true);
}

OBSmiNode::OBSmiNode(OBAtom *atom) 
{
  _atom = atom;
  _parent = NULL;
  _nextnode.clear();
  _nextbond.clear();
}

void OBSmiNode::SetNextNode(OBSmiNode *node,OBBond *bond) 
{
  _nextnode.push_back(node);
  _nextbond.push_back(bond);
}

OBSmiNode::~OBSmiNode()
{
  vector<OBSmiNode*>::iterator i;
  for (i = _nextnode.begin();i != _nextnode.end();i++)
    delete (*i);
}


bool WriteTheSmiles(OBMol & mol,char *out)
{
  char buffer[2*BUFF_SIZE];

  OBMol2Smi m2s;

  m2s.Init();
  m2s.CorrectAromaticAmineCharge(mol);
  m2s.CreateSmiString(mol,buffer);

  strcpy(out,buffer);
  return(true);

}

}

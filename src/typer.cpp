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

#include "babelconfig.h"
#include "mol.h"
#include "typer.h"
#include "atomtyp.h"
#include "aromatic.h"

#ifdef WIN32
#pragma warning (disable : 4786)
#endif

using namespace std;

namespace OpenBabel {

OBAromaticTyper  aromtyper;
OBAtomTyper      atomtyper;

  /*! \class OBAtomTyper
      \brief Assigns atom types, hybridization, implicit valence and formal charges

  The OBAtomTyper class is designed to read in a list of atom typing
  rules and apply them to molecules. The code that performs atom
  typing is not usually used directly as atom typing, hybridization
  assignment, implicit valence assignment and charge are all done
  automatically when their corresponding values are requested of
  atoms:
\code
  atom->GetType();
  atom->GetFormalCharge();
  atom->GetHyb();
\endcode
*/
OBAtomTyper::OBAtomTyper()
{
  _init = false;
  _dir = DATADIR;
  _envvar = "BABEL_DATADIR";
  _filename = "atomtyp.txt";
  _subdir = "data";
  _dataptr = AtomTypeData;
}

void OBAtomTyper::ParseLine(const char *buffer)
{
  vector<string> vs;
  OBSmartsPattern *sp;

  if (EQn(buffer,"INTHYB",6))
    {
      tokenize(vs,buffer);
      if (vs.empty() || vs.size() < 3) return;
      sp = new OBSmartsPattern;
      if (sp->Init(vs[1]))
	_vinthyb.push_back(pair<OBSmartsPattern*,int> (sp,atoi((char*)vs[2].c_str())));
      else {delete sp; sp = NULL;}
    }
  else if (EQn(buffer,"IMPVAL",6))
    {
      tokenize(vs,buffer);
      if (vs.empty() || vs.size() < 3) return;
      sp = new OBSmartsPattern;
      if (sp->Init(vs[1]))
	_vimpval.push_back(pair<OBSmartsPattern*,int> (sp,atoi((char*)vs[2].c_str())));
      else {delete sp; sp = NULL;}
    }
  else if (EQn(buffer,"EXTTYP",6))
    {
      tokenize(vs,buffer);
      if (vs.empty() || vs.size() < 3) return;
      sp = new OBSmartsPattern;
      if (sp->Init(vs[1]))
	_vexttyp.push_back(pair<OBSmartsPattern*,string> (sp,vs[2]));
      else {delete sp; sp = NULL;}
    }
}

OBAtomTyper::~OBAtomTyper()
{
  vector<pair<OBSmartsPattern*,int> >::iterator i;
  for (i = _vinthyb.begin();i != _vinthyb.end();i++) {delete i->first; i->first = NULL;}
  for (i = _vimpval.begin();i != _vimpval.end();i++) {delete i->first; i->first = NULL;}

  vector<pair<OBSmartsPattern*,string> >::iterator j;
  for (j = _vexttyp.begin();j != _vexttyp.end();j++) {delete j->first; j->first = NULL;}
  
}

void OBAtomTyper::AssignTypes(OBMol &mol)
{
  if (!_init) Init();

  mol.SetAtomTypesPerceived();

  vector<vector<int> >::iterator j;
  vector<pair<OBSmartsPattern*,string> >::iterator i;

  for (i = _vexttyp.begin();i != _vexttyp.end();i++)
    if (i->first->Match(mol))
      {
		_mlist = i->first->GetMapList();
		for (j = _mlist.begin();j != _mlist.end();j++)
			mol.GetAtom((*j)[0])->SetType(i->second);
      }
}

void OBAtomTyper::AssignHyb(OBMol &mol)
{
  if (!_init) Init();

  aromtyper.AssignAromaticFlags(mol);

  mol.SetHybridizationPerceived();

  OBAtom *atom;
  vector<OBNodeBase*>::iterator k;
  for (atom = mol.BeginAtom(k);atom;atom = mol.NextAtom(k)) atom->SetHyb(0);

  vector<vector<int> >::iterator j;
  vector<pair<OBSmartsPattern*,int> >::iterator i;

  for (i = _vinthyb.begin();i != _vinthyb.end();i++)
    if (i->first->Match(mol))
      {
		_mlist = i->first->GetMapList();
		for (j = _mlist.begin();j != _mlist.end();j++)
			mol.GetAtom((*j)[0])->SetHyb(i->second);
      }
}

void OBAtomTyper::AssignImplicitValence(OBMol &mol)
{
  if (!_init) Init();

  mol.SetImplicitValencePerceived();

  OBAtom *atom;
  vector<OBNodeBase*>::iterator k;
  for (atom = mol.BeginAtom(k);atom;atom = mol.NextAtom(k))
    atom->SetImplicitValence(atom->GetValence());

  vector<vector<int> >::iterator j;
  vector<pair<OBSmartsPattern*,int> >::iterator i;

  for (i = _vimpval.begin();i != _vimpval.end();i++)
    if (i->first->Match(mol))
    {
		_mlist = i->first->GetMapList();
		for (j = _mlist.begin();j != _mlist.end();j++)
			mol.GetAtom((*j)[0])->SetImplicitValence(i->second);
	}

  if (!mol.HasAromaticCorrected()) CorrectAromaticNitrogens(mol);

  for (atom = mol.BeginAtom(k);atom;atom = mol.NextAtom(k))
    if (atom->GetImplicitValence() < atom->GetValence())
      atom->SetImplicitValence(atom->GetValence());
}


void OBAtomTyper::CorrectAromaticNitrogens(OBMol &mol)
{
  if (!_init) Init();

  if (mol.HasAromaticCorrected()) return;
  mol.SetAromaticCorrected();

	return;

  int j;
  OBAtom *atom,*nbr,*a;
  vector<OBNodeBase*>::iterator i,m;
  vector<OBEdgeBase*>::iterator k;
  OBBitVec curr,used,next;
  vector<OBNodeBase*> v_N,v_OS;

  vector<bool> _aromNH;
  _aromNH.clear();
  _aromNH.resize(mol.NumAtoms()+1);

  for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
    if (atom->IsAromatic() && !atom->IsCarbon() && !used[atom->GetIdx()])
    {
      vector<OBNodeBase*> rsys;
      rsys.push_back(atom);
      curr |= atom->GetIdx();
      used |= atom->GetIdx();
      
      while (!curr.Empty())
	{
	  next.Clear();
	  for (j = curr.NextBit(0);j != curr.EndBit();j = curr.NextBit(j))
	    {
	      atom = mol.GetAtom(j);
	      for (nbr = atom->BeginNbrAtom(k);nbr;nbr = atom->NextNbrAtom(k))
		{
		  if (!(*k)->IsAromatic()) continue;
		  if (used[nbr->GetIdx()]) continue;

		  rsys.push_back(nbr);
		  next |= nbr->GetIdx();
		  used |= nbr->GetIdx();
		}
	    }

	  curr = next;
	}

      //count number of electrons in the ring system
      v_N.clear();
      v_OS.clear();
      int nelectrons = 0;
      
      for (m = rsys.begin();m != rsys.end();m++)
	{
	  a = (OBAtom *)*m;

	  int hasExoDoubleBond = false;
	  for (nbr = a->BeginNbrAtom(k);nbr;nbr = a->NextNbrAtom(k))
	    if ((*k)->GetBO() == 2 && !(*k)->IsInRing())
	      if (nbr->IsOxygen() || nbr->IsSulfur() || nbr->IsNitrogen())
		hasExoDoubleBond = true;

	  if (a->IsCarbon() && hasExoDoubleBond) continue;

	  if (a->IsCarbon()) nelectrons++;
	  else
	    if (a->IsOxygen() || a->IsSulfur())
	      {
		v_OS.push_back(a);
		nelectrons += 2;
	      }
	    else
	      if (a->IsNitrogen())
		{
		  v_N.push_back(a); //store nitrogens
		  nelectrons++;
		}
	}
      
      //calculate what the number of electrons should be for aromaticity
      int naromatic = 2+4*((int)((double)(rsys.size()-2)*0.25+0.5));

      if (nelectrons > naromatic) //try to give one electron back to O or S
	for (m = v_OS.begin();m != v_OS.end();m++)
	  {
	    if (naromatic == nelectrons) break;
	    if ((*m)->GetValence() == 2 && (*m)->GetHvyValence() == 2)
	      {
		nelectrons--;
		((OBAtom*)*m)->SetFormalCharge(1);
	      }
	  }

      if (v_N.empty()) continue; //no nitrogens found in ring

      //check for protonated nitrogens
      for (m = v_N.begin();m != v_N.end();m++)
	{
	  if (naromatic == nelectrons) break;
	  if (((OBAtom*)*m)->GetValence() == 3 && 
	      ((OBAtom*)*m)->GetHvyValence() == 2)
	    {nelectrons++;_aromNH[((OBAtom*)*m)->GetIdx()] = true;}
	  if (((OBAtom*)*m)->GetFormalCharge() == -1)
	    {nelectrons++;_aromNH[((OBAtom*)*m)->GetIdx()] = false;}
	}

       //charge up tert nitrogens
      for (m = v_N.begin();m != v_N.end();m++)
	if (((OBAtom*)*m)->GetHvyValence() == 3 
	    && ((OBAtom*)*m)->BOSum() < 5)
	  ((OBAtom*)*m)->SetFormalCharge(1);

      //try to uncharge nitrogens first
      for (m = v_N.begin();m != v_N.end();m++)
	{
	  if (((OBAtom*)*m)->BOSum() > 4) continue; //skip n=O
	  if (naromatic == nelectrons) break;
	  if (((OBAtom*)*m)->GetHvyValence() == 3)
	    {nelectrons++;((OBAtom*)*m)->SetFormalCharge(0);}
	}

      if (naromatic == nelectrons) continue;
      
      //try to protonate amides next
      for (m = v_N.begin();m != v_N.end();m++)
	{
	  if (naromatic == nelectrons) break;
	  if (((OBAtom*)*m)->IsAmideNitrogen() && 
	      ((OBAtom*)*m)->GetValence() == 2 && 
	      ((OBAtom*)*m)->GetHvyValence() == 2 &&
	      !_aromNH[((OBAtom*)*m)->GetIdx()])
	    {nelectrons++;_aromNH[((OBAtom*)*m)->GetIdx()] = true;}
	}
	  
      if (naromatic == nelectrons) continue;

      //protonate amines in 5 membered rings first - try to match kekule
      for (m = v_N.begin();m != v_N.end();m++)
	{
	  if (naromatic == nelectrons) break;
	  if (((OBAtom*)*m)->GetValence() == 2 &&
	      !_aromNH[((OBAtom*)*m)->GetIdx()] &&
	      ((OBAtom*)*m)->IsInRingSize(5) &&
	      ((OBAtom*)*m)->BOSum() == 2)
	    {nelectrons++;_aromNH[((OBAtom*)*m)->GetIdx()] = true;}
	}

      if (naromatic == nelectrons) continue;

      //protonate amines in 5 membered rings first - no kekule restriction
      for (m = v_N.begin();m != v_N.end();m++)
	{
	  if (naromatic == nelectrons) break;
	  if ((*m)->GetValence() == 2 && !_aromNH[(*m)->GetIdx()] &&
	      (*m)->IsInRingSize(5))
	    {nelectrons++; _aromNH[(*m)->GetIdx()] = true;}
	}

      if (naromatic == nelectrons) continue;

      //then others -- try to find an atom w/o a double bond first
      for (m = v_N.begin();m != v_N.end();m++)
	{
	  if (naromatic == nelectrons) break;
	  if (((OBAtom*)*m)->GetHvyValence() == 2 &&
	      !_aromNH[((OBAtom*)*m)->GetIdx()] && 
	      !((OBAtom*)*m)->HasDoubleBond())
	    {nelectrons++; _aromNH[(*m)->GetIdx()] = true;}
	}
      
      for (m = v_N.begin();m != v_N.end();m++)
	{
	  if (naromatic == nelectrons) break;
	  if ((*m)->GetHvyValence() == 2 && !_aromNH[(*m)->GetIdx()])
	    {nelectrons++; _aromNH[(*m)->GetIdx()] = true;}
	}
    }

  for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
    if (_aromNH[atom->GetIdx()] && atom->GetValence() == 2)
	atom->SetImplicitValence(3);
}

OBAromaticTyper::OBAromaticTyper()
{
  _init = false;
  _dir = DATADIR;
  _envvar = "BABEL_DATADIR";
  _filename = "aromatic.txt";
  _subdir = "data";
  _dataptr = AromaticData;
}

void OBAromaticTyper::ParseLine(const char *buffer)
{
  OBSmartsPattern *sp;
  char temp_buffer[BUFF_SIZE];

  if (buffer[0] == '#') return;
  vector<string> vs;
  tokenize(vs,buffer);
  if (!vs.empty() && vs.size() == 3)
    {
      strcpy(temp_buffer,vs[0].c_str());
      sp = new OBSmartsPattern();
      if (sp->Init(temp_buffer)) 
	{
	  _vsp.push_back(sp);
	  _verange.push_back(pair<int,int> 
			     (atoi((char*)vs[1].c_str()),
			      atoi((char*)vs[2].c_str())));
	}
      else
	{delete sp; sp = NULL;}
    }
}

OBAromaticTyper::~OBAromaticTyper()
{
  vector<OBSmartsPattern*>::iterator i;
   for (i = _vsp.begin();i != _vsp.end();i++)
     {delete *i; *i = NULL;}
}

void OBAromaticTyper::AssignAromaticFlags(OBMol &mol)
{
  if (!_init) Init();

  if (mol.HasAromaticPerceived()) return;
  mol.SetAromaticPerceived();

  _vpa.clear();   _vpa.resize(mol.NumAtoms()+1);
  _velec.clear(); _velec.resize(mol.NumAtoms()+1);
  _root.clear();  _root.resize(mol.NumAtoms()+1);
  
  OBBond *bond;
  OBAtom *atom;
  vector<OBNodeBase*>::iterator i;
  vector<OBEdgeBase*>::iterator j;

  //unset all aromatic flags
  for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i)) 
    atom->UnsetAromatic();
  for (bond = mol.BeginBond(j);bond;bond = mol.NextBond(j)) 
    bond->UnsetAromatic();
	
  int idx;
  vector<vector<int> >::iterator m;
  vector<OBSmartsPattern*>::iterator k;

  //mark atoms as potentially aromatic
  for (idx=0,k = _vsp.begin();k != _vsp.end();k++,idx++)
    if ((*k)->Match(mol))
    {
      _mlist = (*k)->GetMapList();
      for (m = _mlist.begin();m != _mlist.end();m++)
	  {
		  _vpa[(*m)[0]] = true;
		  _velec[(*m)[0]] = _verange[idx];
	  }
    }

	//sanity check - exclude all 4 substituted atoms and sp centers
  for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
  {
		if (atom->GetImplicitValence() > 3)
		{
			_vpa[atom->GetIdx()] = false;
			continue;
		}

		switch(atom->GetAtomicNum()) //phosphorus and sulfur may be initially typed as sp3
		{
		case 6:
		case 7:
		case 8:
			if (atom->GetHyb() != 2)
				_vpa[atom->GetIdx()] = false;
			break;
		}
  }

  //propagate potentially aromatic atoms
  for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
    if (_vpa[atom->GetIdx()])
      PropagatePotentialAromatic(atom);

  //select root atoms
  for (bond = mol.BeginBond(j);bond;bond = mol.NextBond(j)) 
    if (bond->IsClosure())
      _root[bond->GetBeginAtomIdx()] = true;

  ExcludeSmallRing(mol); //remove 3 membered rings from consideration

  //loop over root atoms and look for 5-6 membered aromatic rings
  _visit.clear(); _visit.resize(mol.NumAtoms()+1);
  for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
    if (_root[atom->GetIdx()])
      CheckAromaticity(atom,6);
	
	for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
    if (_root[atom->GetIdx()])
      CheckAromaticity(atom,20);

  //for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
  //	  if (atom->IsAromatic())
  //		  cerr << "aro = " <<atom->GetIdx()  << endl;

  //for (bond = mol.BeginBond(j);bond;bond = mol.NextBond(j)) 
  //if (bond->IsAromatic())
  //cerr << bond->GetIdx() << ' ' << bond->IsAromatic() << endl;
}

bool OBAromaticTyper::TraverseCycle(OBAtom *root,OBAtom *atom,OBBond *prev,pair<int,int> &er,int depth)
{
  if (atom == root)
    {
      int i;
      for (i = er.first;i <= er.second;i++)
				if (i%4 == 2 && i > 2)
				  return(true);

      return(false);
    }
  
  if (!depth || !_vpa[atom->GetIdx()] || _visit[atom->GetIdx()]) return(false);
  
  bool result = false;

  depth--;
  er.first  += _velec[atom->GetIdx()].first;
  er.second += _velec[atom->GetIdx()].second;

  _visit[atom->GetIdx()] = true;
  OBAtom *nbr;
  vector<OBEdgeBase*>::iterator i;
  for (nbr = atom->BeginNbrAtom(i);nbr;nbr = atom->NextNbrAtom(i))
    if (*i != prev && (*i)->IsInRing() && _vpa[nbr->GetIdx()])
      {
	if (TraverseCycle(root,nbr,(OBBond*)(*i),er,depth))
	  {
	    result = true;
	    ((OBBond*) *i)->SetAromatic();
	  }
      }

  _visit[atom->GetIdx()] = false;
  if (result) atom->SetAromatic();

  er.first  -= _velec[atom->GetIdx()].first;
  er.second -= _velec[atom->GetIdx()].second;

  return(result);
}

void OBAromaticTyper::CheckAromaticity(OBAtom *atom,int depth)
{
  OBAtom *nbr;
  vector<OBEdgeBase*>::iterator i;

  pair<int,int> erange; 
  for (nbr = atom->BeginNbrAtom(i);nbr;nbr = atom->NextNbrAtom(i))
    if ((*i)->IsInRing() && !(*i)->IsAromatic())
      {
	erange = _velec[atom->GetIdx()];
	
	if (TraverseCycle(atom,nbr,(OBBond *)*i,erange,depth-1))
	  {
	    atom->SetAromatic();
	    ((OBBond*) *i)->SetAromatic();
	  }
      }
}

void OBAromaticTyper::PropagatePotentialAromatic(OBAtom *atom)
{
  int count = 0;
  OBAtom *nbr;
  vector<OBEdgeBase*>::iterator i;
  
  for (nbr = atom->BeginNbrAtom(i);nbr;nbr = atom->NextNbrAtom(i))
    if ((*i)->IsInRing() && _vpa[nbr->GetIdx()])
      count++;

  if (count < 2)
    {
      _vpa[atom->GetIdx()] = false;
      if (count == 1)
	for (nbr = atom->BeginNbrAtom(i);nbr;nbr = atom->NextNbrAtom(i))
	  if ((*i)->IsInRing() && _vpa[nbr->GetIdx()])
	    PropagatePotentialAromatic(nbr);
    }
}

void OBAromaticTyper::ExcludeSmallRing(OBMol &mol)
{
  OBAtom *atom,*nbr1,*nbr2;
  vector<OBNodeBase*>::iterator i;
  vector<OBEdgeBase*>::iterator j,k;

  for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
    if (_root[atom->GetIdx()])
      for (nbr1 = atom->BeginNbrAtom(j);nbr1;nbr1 = atom->NextNbrAtom(j))
	if ((*j)->IsInRing() && _vpa[nbr1->GetIdx()])
	  for (nbr2 = nbr1->BeginNbrAtom(k);nbr2;nbr2 = nbr1->NextNbrAtom(k))
	    if (nbr2 != atom && (*k)->IsInRing() && _vpa[nbr2->GetIdx()])
	      if (atom->IsConnected(nbr2))
		_root[atom->GetIdx()] = false;
}

} //namespace OpenBabel;



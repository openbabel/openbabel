/**********************************************************************
Copyright (C) 2003 by Geoffrey R. Hutchison

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include "mol.h"
#include "bondtyper.h"
#include "bondtyp.h"

#ifdef WIN32
#pragma warning (disable : 4786)
#endif

using namespace std;

namespace OpenBabel {

OBBondTyper  bondtyper;

  /*! \class OBBondTyper
      \brief Assigns bond types for file formats without bond information

  The OBBondTyper class is designed to read in a list of bond typing
  rules and apply them to molecules.
*/
OBBondTyper::OBBondTyper()
{
  _init = false;
  _dir = BABEL_DATADIR;
  _envvar = "BABEL_DATADIR";
  _filename = "bondtyp.txt";
  _subdir = "data";
  _dataptr = BondTypeData;
}

void OBBondTyper::ParseLine(const char *buffer)
{
  vector<string> vs;
  vector<int>    bovector;
  OBSmartsPattern *sp;

  if (buffer[0] != '#')
    {
      tokenize(vs,buffer);
      // Make sure we actually have a SMARTS pattern plus at least one triple
      // and make sure we have the correct number of integers
      if (vs.empty() || vs.size() < 4 || (vs.size() % 3 == 1))
    {
	  cerr << " Error in OBBondTyper. Pattern is incorrect, found "
	       << vs.size() << " tokens." << endl;
	  cerr << " Buffer is: " << buffer << endl;
	  return;
    }
      sp = new OBSmartsPattern;
      if (sp->Init(vs[0]))
	{
	  for (int i; i <= vs.size() ; i++)
	    {
	      bovector.push_back( atoi((char *)vs[i].c_str()) );
	    }
	  _fgbonds.push_back(pair<OBSmartsPattern*,vector<int> >
			     (sp, bovector));
	}
      else {delete sp; sp = NULL;}
    }
}

OBBondTyper::~OBBondTyper()
{
  vector<pair<OBSmartsPattern*, vector<int> > >::iterator i;
  for (i = _fgbonds.begin();i != _fgbonds.end();i++) {delete i->first; i->first = NULL;}
}


/*! This method adds single bonds between all atoms
  closer than their combined atomic covalent radii,
  then "cleans up" making sure bonded atoms are not
  closer than 0.4A and the atom does not exceed its valence. */
void OBBondTyper::ConnectTheDots(OBMol &mol)
{
  if (Empty()) return;
//  if (!Has3D()) return; // not useful on 2D structures

  int j,k,max;
  bool unset = false;
  OBAtom *atom,*nbr;
  vector<OBNodeBase*>::iterator i;
  vector<pair<OBAtom*,double> > zsortedAtoms;
  vector<double> rad; 
  vector<int> zsorted;

  double *c = new double [NumAtoms()*3];
  rad.resize(_natoms);

  for (j = 0, atom = BeginAtom(i) ; atom ; atom = NextAtom(i), j++)
    {
      (atom->GetVector()).Get(&c[j*3]);
      pair<OBAtom*,double> entry(atom, atom->GetVector().z());
      zsortedAtoms.push_back(entry);
    }
  sort(zsortedAtoms.begin(), zsortedAtoms.end(), SortAtomZ);

  max = zsortedAtoms.size();

  for ( j = 0 ; j < max ; j++ )
  {
      atom   = zsortedAtoms[j].first;
      rad[j] = etab.GetCovalentRad(atom->GetAtomicNum());
      zsorted.push_back(atom->GetIdx()-1);
  }

  int idx1, idx2;
  double d2,cutoff,zd;
  for (j = 0 ; j < max ; j++)
  {
    idx1 = zsorted[j];
    for (k = j + 1 ; k < max ; k++ )
      {
        idx2 = zsorted[k];

	// bonded if closer than elemental Rcov + tolerance
        cutoff = SQUARE(rad[j] + rad[k] + 0.45); 

        zd  = SQUARE(c[idx1*3+2] - c[idx2*3+2]);
	if (zd > 25.0 ) break; // bigger than max cutoff

        d2  = SQUARE(c[idx1*3]   - c[idx2*3]);
        d2 += SQUARE(c[idx1*3+1] - c[idx2*3+1]);
        d2 += zd;
  
        if (d2 > cutoff) continue;
	if (d2 < 0.40) continue;
  
        atom = GetAtom(idx1+1); 
        nbr  = GetAtom(idx2+1);
  
        if (atom->IsConnected(nbr)) continue;
        if (atom->IsHydrogen() && nbr->IsHydrogen()) continue;

        AddBond(idx1+1,idx2+1,1);
      }
  }

  // If between BeginModify and EndModify, coord pointers are NULL
  // setup molecule to handle current coordinates

  if (_c == NULL)
  {
	  _c = c;
	  for (atom = BeginAtom(i);atom;atom = NextAtom(i))
		  atom->SetCoordPtr(&_c);
	  _vconf.push_back(c);
	  unset = true;
  }

  // Cleanup -- delete long bonds that exceed max valence
  OBBond *maxbond, *bond;
  double maxlength;
  vector<OBEdgeBase*>::iterator l;
  for (atom = BeginAtom(i);atom;atom = NextAtom(i))
    {
      while (atom->BOSum() > static_cast<unsigned int>(etab.GetMaxBonds(atom->GetAtomicNum())) || atom->SmallestBondAngle() < 45.0)
	{
	  maxbond = atom->BeginBond(l);
	  maxlength = maxbond->GetLength();
	  for (bond = atom->BeginBond(l);bond;bond = atom->NextBond(l))
	    {
	      if (bond->GetLength() > maxlength)
		{
		  maxbond = bond;
		  maxlength = bond->GetLength();
		}
	    }
	  DeleteBond(maxbond);
	}
    }

  if (unset)
  {
	  _c = NULL;
	  for (atom = BeginAtom(i);atom;atom = NextAtom(i))
		  atom->ClearCoordPtr();
	  _vconf.resize(_vconf.size()-1);
  }
  
  delete [] c;
}

/*! This method uses bond angles and geometries from current
  connectivity to guess atom types and then filling empty valences
  with multiple bonds. It currently has a pass to detect some
  frequent functional groups. It still needs a pass to detect aromatic
  rings to "clean up." */
void OBBondTyper::PerceiveBondOrders(OBMol &mol)
{
  if (Empty()) return;
  //  if (!Has3D()) return; // not useful on 2D structures

  OBAtom *atom, *b, *c;
  vector3 v1, v2;
  int angles;
  double degrees;
  vector<OBNodeBase*>::iterator i;
  vector<OBEdgeBase*>::iterator j,k;
  
  //  BeginModify();

  // Pass 1: Assign estimated hybridization based on avg. angles
  for (atom = BeginAtom(i);atom;atom = NextAtom(i))
    {
      degrees = 0.0;
      angles = 0;
      for (b = atom->BeginNbrAtom(j); b; b = atom->NextNbrAtom(j))
	{
	  k = j;
	  for (c = atom->NextNbrAtom(k); c; c = atom->NextNbrAtom(k))
	    {
	      v1 = b->GetVector() - atom->GetVector();
	      v2 = c->GetVector() - atom->GetVector();	
	      degrees += vectorAngle(v1, v2);
	      angles++;
	    }
	}
      if (angles >= 1)
	{
	  if ((degrees / angles) > 155.0)
	    atom->SetHyb(1);
	  else if ((degrees / angles) <= 155.0 && (degrees / angles) > 115)
	    atom->SetHyb(2);
	}
    } // pass 1

  // Make sure upcoming calls to GetHyb() don't kill these temporary values
  SetHybridizationPerceived();
  
  // Pass 2: look for 5-member rings with torsions <= 7.5 degrees
  //         and 6-member rings with torsions <= 12 degrees
  //         (set all atoms with at least two bonds to sp2)

  vector<OBRing*> rlist;
  vector<OBRing*>::iterator ringit;
  vector<int> path;
  double torsions = 0.0;
  int ringAtom;

  if (!HasSSSRPerceived())
    FindSSSR();
  rlist = GetSSSR();
  for (ringit = rlist.begin(); ringit != rlist.end(); ringit++)
    {
      if ((*ringit)->PathSize() == 5)
	{
	  path = (*ringit)->_path;
	  torsions = 
	    GetTorsion(path[0], path[1], path[2], path[3]) +
	    GetTorsion(path[1], path[2], path[3], path[4]) +
	    GetTorsion(path[2], path[3], path[4], path[0]) +
	    GetTorsion(path[3], path[4], path[0], path[1]) +
	    GetTorsion(path[4], path[0], path[1], path[2]) / 5.0;
	  if (torsions <= 7.5)
	    {
	      for (ringAtom = 0; ringAtom != path.size(); ringAtom++)
		{
		  b = GetAtom(path[ringAtom]);
		  if (b->GetValence() == 2 || b->GetValence() == 3)
		    b->SetHyb(2);
		}
	    }
	}
      else if ((*ringit)->PathSize() == 6)
	{
	  path = (*ringit)->_path;
	  torsions = 
	    GetTorsion(path[0], path[1], path[2], path[3]) +
	    GetTorsion(path[1], path[2], path[3], path[4]) +
	    GetTorsion(path[2], path[3], path[4], path[5]) +
	    GetTorsion(path[3], path[4], path[5], path[0]) +
	    GetTorsion(path[4], path[5], path[0], path[1]) +
	    GetTorsion(path[5], path[0], path[1], path[2]) / 6.0;
	  if (torsions <= 12.0)
	    {
	      for (ringAtom = 0; ringAtom != path.size(); ringAtom++)
		{
		  b = GetAtom(path[ringAtom]);
		  if (b->GetValence() == 2 || b->GetValence() == 3)
		    b->SetHyb(2);
		}
	    }
	}
  }

  // Pass 3: "Antialiasing" If an atom marked as sp hybrid isn't 
  //          bonded to another or an sp2 hybrid isn't bonded 
  //          to another (or terminal atoms in both cases)
  //          mark them to a lower hybridization for now
  bool openNbr;
  for (atom = BeginAtom(i);atom;atom = NextAtom(i))
    {
      if (atom->GetHyb() == 2 || atom->GetHyb() == 1)
	{
	  openNbr = false;
	  for (b = atom->BeginNbrAtom(j); b; b = atom->NextNbrAtom(j))
	    {
	      if (b->GetHyb() < 3 || b->GetValence() == 1)
		{
		  openNbr = true;
		  break;
		}
	    }
	  if (!openNbr && atom->GetHyb() == 2)
	    atom->SetHyb(3);
	  else if (!openNbr && atom->GetHyb() == 1)
	    atom->SetHyb(2);
	}
    } // pass 3

  // Pass 4: Check for known functional group patterns and assign bonds
  //         to the canonical form
  //      Currently we have explicit code to do this, but a "bond typer"
  //      is in progress to make it simpler to test and debug.
  OBBond *b1,*b2,*b3;
  OBAtom *a1,*a2,*a3,*a4;
  vector<vector<int> > mlist;
  vector<vector<int> >::iterator l;

  // azide -N=N=N
  OBSmartsPattern azide; azide.Init("[#7D2][#7D2][#7D1]");
  if (azide.Match(*this))
    {
      mlist = azide.GetUMapList();
      for (l = mlist.begin(); l != mlist.end(); l++)
        {
          a1 = GetAtom((*l)[0]);
          a2 = GetAtom((*l)[1]);
          a3 = GetAtom((*l)[2]);
          b1 = a2->GetBond(a1); b2 = a2->GetBond(a3);

          if (!b1 || !b2) continue;
          b1->SetBO(2);
          b2->SetBO(2);
        }
    } // Azide

  // Nitro -NO2
  OBSmartsPattern nitro; nitro.Init("[#8D1][#7D3][#8D1]");
  if (nitro.Match(*this))
    {
      mlist = nitro.GetUMapList();
      for (l = mlist.begin(); l != mlist.end(); l++)
        {
          a1 = GetAtom((*l)[0]);
          a2 = GetAtom((*l)[1]);
          a3 = GetAtom((*l)[2]);
          b1 = a1->GetBond(a2); b2 = a2->GetBond(a3);

          if (!b1 || !b2) continue;
          b1->SetBO(2);
          b2->SetBO(2);
        }      
    } // Nitro

  // Sulphone -SO2-
  OBSmartsPattern sulphone; sulphone.Init("[#16D4]([#8D1])([#8D1])(*)(*)");
  if (sulphone.Match(*this))
    {
      mlist = sulphone.GetUMapList();
      for (l = mlist.begin(); l != mlist.end(); l++)
        {
          a1 = GetAtom((*l)[0]);
          a2 = GetAtom((*l)[1]);
          a3 = GetAtom((*l)[2]);
	  a4 = GetAtom((*l)[3]);
          b1 = a1->GetBond(a2); b2 = a1->GetBond(a3); b3 = a1->GetBond(a4);

          if (!b1 || !b2 || !b3) continue;
          b1->SetBO(2);
          b2->SetBO(2);
	  b3->SetBO(1);
	  a4 = GetAtom((*l)[4]);
	  b3 = a1->GetBond(a4);
	  if (!b3) continue;
	  b3->SetBO(1);
        }
    } // Sulfone

  // Pass 5: Check for aromatic rings and assign bonds as appropriate
  // This is just a quick and dirty approximation that marks everything
  //  as potentially aromatic

  // This doesn't work perfectly, but it's pretty decent.
  //  Need to have a list of SMARTS patterns for common rings
  //  which would "break ties" on complicated multi-ring systems
  // (Most of the current problems lie in the interface with the 
  //   Kekulize code anyway, not in marking everything as potentially aromatic)

  bool typed; // has this ring been typed?
  unsigned int loop, loopSize;
  for (ringit = rlist.begin(); ringit != rlist.end(); ringit++)
    {
      typed = false;
      loopSize = (*ringit)->PathSize();
      if (loopSize == 5 || loopSize == 6)
	{
	  path = (*ringit)->_path;
	  for(loop = 0; loop < loopSize; loop++)
	    {
	      atom = GetAtom(path[loop]);
	      if(atom->HasBondOfOrder(2) || atom->HasBondOfOrder(3)
		 || atom->GetHyb() != 2)
		{
		  typed = true;
		  break;
		}
	    }

	  if (!typed)
	    for(loop = 0; loop < loopSize; loop++)
	      {
		//		cout << " set aromatic " << path[loop] << endl;
		(GetBond(path[loop], path[(loop+1) % loopSize]))->SetBO(5);
		(GetBond(path[loop], path[(loop+1) % loopSize]))->UnsetKekule();
	      }
	}
    }
  _flags &= (~(OB_KEKULE_MOL));
  //  cout << " calling Kekulize " << endl;
  Kekulize();

  // Pass 6: Assign remaining bond types, ordered by atom electronegativity
  vector<pair<OBAtom*,double> > sortedAtoms;
  vector<double> rad; 
  vector<int> sorted;
  int iter, max;
  double maxElNeg, shortestBond, currentElNeg;

  for (atom = BeginAtom(i) ; atom ; atom = NextAtom(i))
    {
      pair<OBAtom*,double> entry(atom, etab.GetElectroNeg(atom->GetAtomicNum()));
      sortedAtoms.push_back(entry);
    }
  sort(sortedAtoms.begin(), sortedAtoms.end(), SortAtomZ);

  max = sortedAtoms.size();

  for (iter = 0 ; iter < max ; iter++ )
  {
      atom = sortedAtoms[iter].first;
      if ( (atom->GetHyb() == 1 || atom->GetValence() == 1)
	   && atom->BOSum() + 2  <= static_cast<unsigned int>(etab.GetMaxBonds(atom->GetAtomicNum())) )
	{
	  // loop through the neighbors looking for a hybrid or terminal atom
	  // (and pick the one with highest electronegativity first)
	  // *or* pick a neighbor that's a terminal atom

	  if (atom->HasNonSingleBond() || 
	      (atom->GetAtomicNum() == 7 && atom->BOSum() + 2 > 3))
	    continue;

	  maxElNeg = 0.0;
	  shortestBond = 5000.0;
	  c = NULL;
	  for (b = atom->BeginNbrAtom(j); b; b = atom->NextNbrAtom(j))
	    {
	      currentElNeg = etab.GetElectroNeg(b->GetAtomicNum());
	      if ( (b->GetHyb() == 1 || b->GetValence() == 1)
		   && b->BOSum() + 2 <= static_cast<unsigned int>(etab.GetMaxBonds(b->GetAtomicNum()))
		   && (currentElNeg > maxElNeg ||
		       (IsNear(currentElNeg,maxElNeg)
			&& (atom->GetBond(b))->GetLength() < shortestBond)) )
		{
		  if (b->HasNonSingleBond() || 
		      (b->GetAtomicNum() == 7 && b->BOSum() + 2 > 3))
		    continue;

		  shortestBond = (atom->GetBond(b))->GetLength();
		  maxElNeg = etab.GetElectroNeg(b->GetAtomicNum());
		  c = b; // save this atom for later use
		}
	    }
	  if (c)
	    (atom->GetBond(c))->SetBO(3);
	}
      else if ( (atom->GetHyb() == 2 || atom->GetValence() == 1)
		&& atom->BOSum() + 1 <= static_cast<unsigned int>(etab.GetMaxBonds(atom->GetAtomicNum())) )
	{
	  // as above
	  if (atom->HasNonSingleBond() ||
	      (atom->GetAtomicNum() == 7 && atom->BOSum() + 1 > 3))
	    continue;

	  maxElNeg = 0.0;
	  shortestBond = 5000.0;
	  c = NULL;
	  for (b = atom->BeginNbrAtom(j); b; b = atom->NextNbrAtom(j))
	    {
	      currentElNeg = etab.GetElectroNeg(b->GetAtomicNum());
	      if ( (b->GetHyb() == 2 || b->GetValence() == 1)
		   && b->BOSum() + 1 <= static_cast<unsigned int>(etab.GetMaxBonds(b->GetAtomicNum()))
		   && (currentElNeg > maxElNeg ||
		       (IsNear(currentElNeg,maxElNeg)
                       && (atom->GetBond(b))->GetLength() < shortestBond)) )
		{
		  if (b->HasNonSingleBond() ||
		      (b->GetAtomicNum() == 7 && b->BOSum() + 1 > 3))
		    continue;

		  shortestBond = (atom->GetBond(b))->GetLength();
		  maxElNeg = etab.GetElectroNeg(b->GetAtomicNum());
		  c = b; // save this atom for later use
		}
	    }
	  if (c)
	    (atom->GetBond(c))->SetBO(2);
	}
 } // pass 6

  // Now let the atom typer go to work again
  _flags &= (~(OB_HYBRID_MOL));
  _flags &= (~(OB_KEKULE_MOL));
  _flags &= (~(OB_AROMATIC_MOL));
  _flags &= (~(OB_ATOMTYPES_MOL));
  _flags &= (~(OB_IMPVAL_MOL));
  //  EndModify(true); // "nuke" perceived data
}

} //namespace OpenBabel;



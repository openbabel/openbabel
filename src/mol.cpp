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
#include "binary.h"
#include "phmodel.h"

using namespace std;

namespace OpenBabel {

extern bool SwabInt;
extern OBPhModel phmodel;

//
// OBMol member functions
//

bool SortVVInt(const vector<int> &a,const vector<int> &b)
{
  return(a.size() > b.size());
}

bool SortAtomZ(const pair<OBAtom*,float> &a, const pair<OBAtom*,float> &b)
{
  return (a.second < b.second);
}

float OBMol::GetTorsion(int a,int b,int c,int d)
{
  return(CalcTorsionAngle(((OBAtom*)_vatom[a-1])->GetVector(),
                          ((OBAtom*)_vatom[b-1])->GetVector(),
                          ((OBAtom*)_vatom[c-1])->GetVector(),
                          ((OBAtom*)_vatom[d-1])->GetVector()));
}

void OBMol::SetTorsion(OBAtom *a,OBAtom *b,OBAtom *c, OBAtom *d, float ang)
{
  vector<int> tor;
  vector<int> atoms;

  tor.push_back(a->GetCIdx());
  tor.push_back(b->GetCIdx());
  tor.push_back(c->GetCIdx());
  tor.push_back(d->GetCIdx());

  FindChildren(atoms, b->GetIdx(), c->GetIdx());
  int j;
  for (j = 0 ; (unsigned)j < atoms.size() ; j++ )
    atoms[j] = (atoms[j] - 1) * 3;

  float v1x,v1y,v1z,v2x,v2y,v2z,v3x,v3y,v3z;
  float c1x,c1y,c1z,c2x,c2y,c2z,c3x,c3y,c3z;
  float c1mag,c2mag,radang,costheta,m[9];
  float x,y,z,mag,rotang,sn,cs,t,tx,ty,tz;

  //calculate the torsion angle

  v1x = _c[tor[0]]   - _c[tor[1]];   v2x = _c[tor[1]]   - _c[tor[2]];
  v1y = _c[tor[0]+1] - _c[tor[1]+1]; v2y = _c[tor[1]+1] - _c[tor[2]+1];
  v1z = _c[tor[0]+2] - _c[tor[1]+2]; v2z = _c[tor[1]+2] - _c[tor[2]+2];
  v3x = _c[tor[2]]   - _c[tor[3]];
  v3y = _c[tor[2]+1] - _c[tor[3]+1];
  v3z = _c[tor[2]+2] - _c[tor[3]+2];

  c1x = v1y*v2z - v1z*v2y;   c2x = v2y*v3z - v2z*v3y;
  c1y = -v1x*v2z + v1z*v2x;  c2y = -v2x*v3z + v2z*v3x;
  c1z = v1x*v2y - v1y*v2x;   c2z = v2x*v3y - v2y*v3x;
  c3x = c1y*c2z - c1z*c2y;
  c3y = -c1x*c2z + c1z*c2x;
  c3z = c1x*c2y - c1y*c2x; 

  c1mag = SQUARE(c1x)+SQUARE(c1y)+SQUARE(c1z);
  c2mag = SQUARE(c2x)+SQUARE(c2y)+SQUARE(c2z);
  if (c1mag*c2mag < 0.01) costheta = 1.0; //avoid div by zero error
  else costheta = (c1x*c2x + c1y*c2y + c1z*c2z)/(sqrt(c1mag*c2mag));

  if (costheta < -0.999999) costheta = -0.999999f;
  if (costheta >  0.999999) costheta =  0.999999f;
                                  
  if ((v2x*c3x + v2y*c3y + v2z*c3z) > 0.0) radang = -acos(costheta);
  else                                     radang = acos(costheta);

  //
  // now we have the torsion angle (radang) - set up the rot matrix
  //

  //find the difference between current and requested
  rotang = ang - radang; 

  sn = sin(rotang); cs = cos(rotang);t = 1 - cs;
  //normalize the rotation vector
  mag = sqrt(SQUARE(v2x)+SQUARE(v2y)+SQUARE(v2z));
  x = v2x/mag; y = v2y/mag; z = v2z/mag;

  //set up the rotation matrix
  m[0]= t*x*x + cs;     m[1] = t*x*y + sn*z;  m[2] = t*x*z - sn*y;
  m[3] = t*x*y - sn*z;  m[4] = t*y*y + cs;    m[5] = t*y*z + sn*x;
  m[6] = t*x*z + sn*y;  m[7] = t*y*z - sn*x;  m[8] = t*z*z + cs;

  //
  //now the matrix is set - time to rotate the atoms
  //
  tx = _c[tor[1]];ty = _c[tor[1]+1];tz = _c[tor[1]+2];
  vector<int>::iterator i;
  for (i = atoms.begin(),j=*i;i != atoms.end();i++,j=*i)
    {
      _c[j] -= tx;_c[j+1] -= ty;_c[j+2]-= tz;
      x = _c[j]*m[0] + _c[j+1]*m[1] + _c[j+2]*m[2];
      y = _c[j]*m[3] + _c[j+1]*m[4] + _c[j+2]*m[5];
      z = _c[j]*m[6] + _c[j+1]*m[7] + _c[j+2]*m[8];
      _c[j] = x; _c[j+1] = y; _c[j+2] = z;
      _c[j] += tx;_c[j+1] += ty;_c[j+2] += tz;
    }
}


float OBMol::GetTorsion(OBAtom *a,OBAtom *b,OBAtom *c,OBAtom *d)
{
  return(CalcTorsionAngle(a->GetVector(),
                          b->GetVector(),
                          c->GetVector(),
                          d->GetVector()));
}

void OBMol::ContigFragList(vector<vector<int> >&cfl) 
     //each vector<int> contains the atom numbers of a contig fragment
     //the vectors are sorted by size from largest to smallest
{
  int j;
  OBAtom *atom;
  OBBond *bond;
  vector<OBNodeBase*>::iterator i;
  vector<OBEdgeBase*>::iterator k;
  OBBitVec used,curr,next,frag;
  vector<int> tmp;

  used.Resize(NumAtoms()+1);
  curr.Resize(NumAtoms()+1);
  next.Resize(NumAtoms()+1);
  frag.Resize(NumAtoms()+1);

  while ((unsigned)used.CountBits() < NumAtoms())
    {
      curr.Clear();frag.Clear();
      for (atom = BeginAtom(i);atom;atom = NextAtom(i))
          if (!used.BitIsOn(atom->GetIdx())) { curr.SetBitOn(atom->GetIdx()); break;}

      frag |= curr;
      while (!curr.IsEmpty())
        {
          next.Clear();
          for (j = curr.NextBit(-1);j != curr.EndBit();j = curr.NextBit(j))
            {
              atom = GetAtom(j);
              for (bond = atom->BeginBond(k);bond;bond = atom->NextBond(k))
              if (!used.BitIsOn(bond->GetNbrAtomIdx(atom)))
              next.SetBitOn(bond->GetNbrAtomIdx(atom));
            }

          used |= curr;
          used |= next;
          frag |= next;
          curr = next;
        }

      tmp.clear();
      frag.ToVecInt(tmp);
      cfl.push_back(tmp);
    }

  sort(cfl.begin(),cfl.end(),SortVVInt);
}

void OBMol::FindLargestFragment(OBBitVec &lf) 
     //each vector<int> contains the atom numbers of a contig fragment
     //the vectors are sorted by size from largest to smallest
{
  int j;
  OBAtom *atom;
  OBBond *bond;
  vector<OBNodeBase*>::iterator i;
  vector<OBEdgeBase*>::iterator k;
  OBBitVec used,curr,next,frag;

  lf.Clear();
  while ((unsigned)used.CountBits() < NumAtoms())
    {
      curr.Clear();frag.Clear();
      for (atom = BeginAtom(i);atom;atom = NextAtom(i))
          if (!used.BitIsOn(atom->GetIdx())) { curr.SetBitOn(atom->GetIdx()); break;}

      frag |= curr;
      while (!curr.IsEmpty())
        {
          next.Clear();
          for (j = curr.NextBit(-1);j != curr.EndBit();j = curr.NextBit(j))
            {
              atom = GetAtom(j);
              for (bond = atom->BeginBond(k);bond;bond = atom->NextBond(k))
              if (!used.BitIsOn(bond->GetNbrAtomIdx(atom))) next.SetBitOn(bond->GetNbrAtomIdx(atom));
            }

          used |= curr;
          used |= next;
          frag |= next;
          curr = next;
        }

      if (lf.Empty() || lf.CountBits() < frag.CountBits()) lf = frag;
    }
}


void OBMol::FindChildren(vector<OBAtom*> &children,OBAtom *bgn,OBAtom *end)
     //locates all atoms for which there exists a path to 'second'
     //without going through 'first' 
     //children does not include 'second'
{
  OBBitVec used,curr,next;

  used |= bgn->GetIdx();
  used |= end->GetIdx();
  curr |= end->GetIdx();
  children.clear();

  int i;
  OBAtom *atom,*nbr;
  vector<OBEdgeBase*>::iterator j;

  for (;;)
    {
      next.Clear();
      for (i = curr.NextBit(-1);i != curr.EndBit();i = curr.NextBit(i))
        {
          atom = GetAtom(i);
          for (nbr = atom->BeginNbrAtom(j);nbr;nbr = atom->NextNbrAtom(j))
              if (!used[nbr->GetIdx()])
                {
                  children.push_back(nbr);
                  next |= nbr->GetIdx();
                  used |= nbr->GetIdx();
                }
        }
      if (next.Empty()) break;
      curr = next;
    }
}

void OBMol::FindChildren(vector<int> &children,int first,int second)
     //locates all atoms for which there exists a path to 'end'
     //without going through 'bgn' 
     //children does not include 'end'
{
  int i;
  OBBitVec used,curr,next;

  used.SetBitOn(first);
  used.SetBitOn(second);
  curr.SetBitOn(second);

  OBAtom *atom;
  OBBond *bond;
  vector<OBEdgeBase*>::iterator j;

  while (!curr.IsEmpty())
    {
      next.Clear();
      for (i = curr.NextBit(-1);i != curr.EndBit();i = curr.NextBit(i))
        {
          atom = GetAtom(i);
          for (j = atom->BeginBonds(),bond=(OBBond *)*j;
	       j != atom->EndBonds();j++,bond=(OBBond *)*j)
          if (!used.BitIsOn(bond->GetNbrAtomIdx(atom)))
          next.SetBitOn(bond->GetNbrAtomIdx(atom));
        }

      used |= next;
      curr = next;
    }

  used.SetBitOff(first);used.SetBitOff(second);
  used.ToVecInt(children);
}

/*!
**\brief calculates the graph theoretical distance of each atom
** vector is indexed from zero
*/

bool OBMol::GetGTDVector(vector<int> &gtd)
     //calculates the graph theoretical distance for every atom 
     //and puts it into gtd
{
  gtd.clear();
  gtd.resize(NumAtoms());
  
  int gtdcount,natom;
  OBBitVec used,curr,next;
  OBAtom *atom,*atom1;
  OBBond *bond;
  vector<OBNodeBase*>::iterator i;
  vector<OBEdgeBase*>::iterator j;

  next.Clear();

  for (atom = BeginAtom(i);atom;atom = NextAtom(i))
    {
      gtdcount = 0;
      used.Clear();curr.Clear();
      used.SetBitOn(atom->GetIdx());
      curr.SetBitOn(atom->GetIdx());

      while (!curr.IsEmpty())
        {
          next.Clear();
          for (natom = curr.NextBit(-1);natom != curr.EndBit();natom = curr.NextBit(natom))
            {
              atom1 = GetAtom(natom);
              for (bond = atom1->BeginBond(j);bond;bond = atom1->NextBond(j))
                  if (!used.BitIsOn(bond->GetNbrAtomIdx(atom1)) && !curr.BitIsOn(bond->GetNbrAtomIdx(atom1)))
                      if (!(bond->GetNbrAtom(atom1))->IsHydrogen()) next.SetBitOn(bond->GetNbrAtomIdx(atom1));
            }

          used |= next;
          curr = next;
          gtdcount++;
        }
      gtd[atom->GetIdx()-1] = gtdcount;
    }
  return(true);
}

/*!
**\brief calculates a set of graph invariant index using
** graph theoretical distance, number of connected heavy atoms,
** aromatic boolean, ring boolean, atomic number, and 
** summation of bond orders connected to the atom
** vector is indexed from zero
*/
void OBMol::GetGIVector(vector<unsigned int> &vid)
{
	vid.clear();
	vid.resize(NumAtoms()+1);

	vector<int> v;
	GetGTDVector(v);

	int i;
	OBAtom *atom;
	vector<OBNodeBase*>::iterator j;
	for (i=0,atom = BeginAtom(j);atom;atom = NextAtom(j),i++)
	{
		vid[i] =  (unsigned int)v[i];
		vid[i] += (unsigned int)(atom->GetHvyValence()*100);
		vid[i] += (unsigned int)(((atom->IsAromatic()) ? 1 : 0)*1000);
		vid[i] += (unsigned int)(((atom->IsInRing()) ? 1 : 0)*10000);
		vid[i] += (unsigned int)(atom->GetAtomicNum()*100000);
		vid[i] += (unsigned int)(atom->GetImplicitValence()*10000000);
	}
}

static bool OBComparePairSecond(const pair<OBAtom*,unsigned int> &a,const pair<OBAtom*,unsigned int> &b)
{
	return(a.second < b.second);
}

static bool OBComparePairFirst(const pair<OBAtom*,unsigned int> &a,const pair<OBAtom*,unsigned int> &b)
{
	return(a.first->GetIdx() < b.first->GetIdx());
}

static void ClassCount(vector<pair<OBAtom*,unsigned int> > &vp,int &count)
//counts the number of unique symmetry classes in a list
{
	count = 0;
	unsigned int id;
	vector<pair<OBAtom*,unsigned int> >::iterator k;
	sort(vp.begin(),vp.end(),OBComparePairSecond);

	for (k = vp.begin();k != vp.end();k++)
	{
		if (k == vp.begin()) 
		{
			id = k->second;
			k->second = count = 0;
		} 
		else 
			if (k->second != id)
			{
				id = k->second;
				k->second = ++count;
			}
			else k->second = count;
	}
	count++;
}

static void	CreateNewClassVector(vector<pair<OBAtom*,unsigned int> > &vp1,vector<pair<OBAtom*,unsigned int> > &vp2)
//creates a new vector of symmetry classes base on an existing vector
//helper routine to GetGIDVector
{
	int m,id;
	OBAtom *nbr;
	vector<OBEdgeBase*>::iterator j;
	vector<unsigned int>::iterator k;
	vector<pair<OBAtom*,unsigned int> >::iterator i;
	sort(vp1.begin(),vp1.end(),OBComparePairFirst);
	vp2.clear();
	for (i = vp1.begin();i != vp1.end();i++)
	{
		vector<unsigned int> vtmp;
		for (nbr = i->first->BeginNbrAtom(j);nbr;nbr = i->first->NextNbrAtom(j))
			vtmp.push_back(vp1[nbr->GetIdx()-1].second);
		sort(vtmp.begin(),vtmp.end(),OBCompareUnsigned);
		for (id=i->second,m=100,k = vtmp.begin();k != vtmp.end();k++,m*=100)
			id += *k * m;

		vp2.push_back(pair<OBAtom*,unsigned int> (i->first,id));
	}
}

/*!
**\brief Calculates a set of symmetry identifiers for a mol.
** Atoms with the same sym ident are symmetrically equivalent
** vector is indexed from zero
*/
void OBMol::GetGIDVector(vector<unsigned int> &vgid)
{
	vector<unsigned int> vgi;
	GetGIVector(vgi);  //get vector of graph invariants

	int i;
	OBAtom *atom;
	vector<OBNodeBase*>::iterator j;
	vector<pair<OBAtom*,unsigned int> > vp1,vp2;
	for (i=0,atom = BeginAtom(j);atom;atom = NextAtom(j),i++)
		vp1.push_back(pair<OBAtom*,unsigned int> (atom,vgi[i]));

	int nclass1,nclass2; //number of classes
	ClassCount(vp1,nclass1);
	
	if (nclass1 < NumAtoms())
	{
		for (i = 0;i < 100;i++) //sanity check - shouldn't ever hit this number
		{
			CreateNewClassVector(vp1,vp2);
			ClassCount(vp2,nclass2);
			vp1 = vp2;
			if (nclass1 == nclass2) break;
			nclass1 = nclass2;
		}
	}

	vgid.clear();
	sort(vp1.begin(),vp1.end(),OBComparePairFirst);
	vector<pair<OBAtom*,unsigned int> >::iterator k;
	for (k = vp1.begin();k != vp1.end();k++)
		vgid.push_back(k->second);
}

unsigned int OBMol::NumHvyAtoms()
{
  OBAtom *atom;
  vector<OBNodeBase*>::iterator(i);
  unsigned int count = 0;

  for(atom = this->BeginAtom(i);atom;atom = this->NextAtom(i))
    {
      if(!atom->IsHydrogen()) count++;
    }

  return(count);
}

unsigned int OBMol::NumRotors()
{
  OBBond *bond;
  vector<OBEdgeBase*>::iterator i;

  unsigned int count = 0;
  for (bond = BeginBond(i);bond;bond = NextBond(i))
    if (bond->IsRotor())
      count++;

  return(count);
}

OBAtom *OBMol::GetAtom(int idx)
     //0 < idx <= NumAtoms
     //returns a pointer to the atom after a safety check
{
  if ((unsigned)idx < 1 || (unsigned)idx > NumAtoms())
    {
      cerr << idx << ' ' << NumAtoms() << endl;

      obAssert(false);
      ThrowError("Requested Atom Out of Range");
      return((OBAtom*)NULL);
    }

  return((OBAtom*)_vatom[idx-1]);
}

OBAtom *OBMol::GetFirstAtom()
{
  return((_vatom.empty()) ? (OBAtom*)NULL : (OBAtom*)_vatom[0]);
}

OBBond *OBMol::GetBond(int idx)
     //0 <= idx < NumBonds
     //returns a pointer to the bond after a safety check
{
  if (idx < 0 || (unsigned)idx >= NumBonds())
    {
      ThrowError("Requested Bond Out of Range");
      return((OBBond*)NULL);
    }

  return((OBBond*)_vbond[idx]);
}

OBBond *OBMol::GetBond(int bgn, int end)
{
  return(GetBond(GetAtom(bgn),GetAtom(end)));
}

OBBond *OBMol::GetBond(OBAtom *bgn,OBAtom *end)
{
  OBAtom *nbr;
  vector<OBEdgeBase*>::iterator i;
  
  for (nbr = bgn->BeginNbrAtom(i);nbr;nbr = bgn->NextNbrAtom(i))
    if (nbr == end)
      return((OBBond *)*i);

  //obAssert(false); //should never get here

  return(NULL); //just to keep the SGI compiler happy
}

OBResidue *OBMol::GetResidue(int idx)
{
    if (idx < 0 || (unsigned)idx >= NumResidues())
    {
        ThrowError("Requested Residue Out of Range");
        return((OBResidue*)NULL);
    }

    return (_residue[idx]);
}

vector<OBRing*> &OBMol::GetSSSR()
{
	if (!HasSSSRPerceived()) 
		FindSSSR();

	if (!HasData(obRingData)) 
		SetData(new OBRingData);

	OBRingData *rd = (OBRingData *) GetData(obRingData);
	return(rd->GetData());
}

float OBMol::GetMolWt()
{
  float molwt=0.0;
  OBAtom *atom;
  vector<OBNodeBase*>::iterator i;

  for (atom = BeginAtom(i);atom;atom = NextAtom(i))
    molwt += etab.GetMass(atom->GetAtomicNum());

  return(molwt);
}

OBMol &OBMol::operator=(const OBMol &source)
     //only atom and bond info is copied from src to dest
     //Conformers are now copied also, MM 2/7/01
     //Rotamers and residue information are copied, MM 4-27-01
{
  OBMol &src = (OBMol &)source;
  vector<OBNodeBase*>::iterator i;
  vector<OBEdgeBase*>::iterator j;
  OBAtom *atom;
  OBBond *bond;
 
  Clear();
  BeginModify();

  _vatom.reserve(src.NumAtoms());
  _vbond.reserve(src.NumBonds());

  for (atom = src.BeginAtom(i);atom;atom = src.NextAtom(i)) AddAtom(*atom);
  for (bond = src.BeginBond(j);bond;bond = src.NextBond(j)) AddBond(*bond);

  _itype = src.GetInputType();
  _otype = src.GetOutputType();

  this->_title  = src.GetTitle();
  this->_energy = src.GetEnergy();

  EndModify();

  //Copy Residue information
  unsigned int NumRes = src.NumResidues();
  if (NumRes) {
      unsigned int k;
      OBResidue *src_res=NULL;
      OBResidue *res=NULL;
      OBAtom *src_atom=NULL;
      OBAtom *atom=NULL;
      vector<OBAtom*>::iterator ii;
      for (k=0 ; k<NumRes ; k++) {
          res = NewResidue();
          src_res = src.GetResidue(k);
          res->SetName(src_res->GetName());
          res->SetNum(src_res->GetNum());
          res->SetChain(src_res->GetChain());
          res->SetChainNum(src_res->GetChainNum());
          for (src_atom=src_res->BeginAtom(ii) ; src_atom ; src_atom=src_res->NextAtom(ii)) {
              atom = GetAtom(src_atom->GetIdx());
              res->AddAtom(atom);
              res->SetAtomID(atom,src_res->GetAtomID(src_atom));
              res->SetHetAtom(atom,src_res->IsHetAtom(src_atom));
              res->SetSerialNum(atom,src_res->GetSerialNum(src_atom));
            }
        }
    }

  //Copy conformer information
  if (src.NumConformers() > 1) {
      int k,l;
      vector<float*> conf;
      float* xyz = NULL;
      for (k=0 ; k<src.NumConformers() ; k++) {
          xyz = new float [3*src.NumAtoms()];
          for (l=0 ; l<(int) (3*src.NumAtoms()) ; l++) xyz[l] = src.GetConformer(k)[l];
          conf.push_back(xyz);
        }
      SetConformers(conf);
    } 

  //Copy rotamer list
  OBRotamerList *rml = (OBRotamerList *)src.GetData(obRotamerList);
  //if (rml) {cout << "DEBUG : OBMol assignment operator.  Source HAS RotamerList" << endl;}
  //else     {cout << "DEBUG : OBMol assignment operator.  Source does NOT have RotamerList" << endl;}
  if (rml && rml->NumAtoms() == src.NumAtoms()) {
      //Destroy old rotamer list if necessary
      if ((OBRotamerList *)GetData(obRotamerList)) {
          DeleteData(obRotamerList);
        }

      //Set base coordinates
      OBRotamerList *cp_rml = new OBRotamerList;
      unsigned int k,l;
      vector<float*> bc; 
      float *c=NULL;
      float *cc=NULL;
      for (k=0 ; k<rml->NumBaseCoordinateSets() ; k++) {
          c = new float [3*rml->NumAtoms()];
          cc = rml->GetBaseCoordinateSet(k);
          for (l=0 ; l<3*rml->NumAtoms() ; l++) c[l] = cc[l];
          bc.push_back(c);
        }
      if (rml->NumBaseCoordinateSets()) cp_rml->SetBaseCoordinateSets(bc,rml->NumAtoms());

      //Set reference array
      unsigned char *ref = new unsigned char [rml->NumRotors()*4];
      if (ref)
	{
	  rml->GetReferenceArray(ref);
	  cp_rml->Setup((*this),ref,rml->NumRotors());
	  delete [] ref;
	}

      //Set Rotamers
      unsigned char *rotamers = new unsigned char [(rml->NumRotors()+1)*rml->NumRotamers()];
      if (rotamers)
	{
	  vector<unsigned char*>::iterator kk;
	  unsigned int idx=0;
	  for (kk = rml->BeginRotamer();kk != rml->EndRotamer();kk++)
	    {
	      memcpy(&rotamers[idx],(const unsigned char*)*kk,sizeof(unsigned char)*(rml->NumRotors()+1));
	      idx += sizeof(unsigned char)*(rml->NumRotors()+1);
	    }
	  cp_rml->AddRotamers(rotamers,rml->NumRotamers());
	  delete [] rotamers;
	}
      SetData(cp_rml); 
    }

  
  //Copy pose information
  _pose = src._pose;
  if (!_pose.empty()) SetPose(src.CurrentPoseIndex());

  return(*this);
}

OBMol &OBMol::operator+=(const OBMol &source)
{
  OBMol &src = (OBMol &)source;
  vector<OBNodeBase*>::iterator i;
  vector<OBEdgeBase*>::iterator j;
  OBAtom *atom;
  OBBond *bond;

  BeginModify();

  int prevatms = NumAtoms();

  _title += "_" + string(src.GetTitle());
  
  for (atom = src.BeginAtom(i) ; atom ; atom = src.NextAtom(i))
    AddAtom(*atom);
  for (bond = src.BeginBond(j) ; bond ; bond = src.NextBond(j))
    AddBond(bond->GetBeginAtomIdx() + prevatms, bond->GetEndAtomIdx() + prevatms, bond->GetBO());

  EndModify();

  return(*this);
}

bool OBMol::Clear()
     //clear all the info in a molecule
{
  vector<OBNodeBase*>::iterator i;
  vector<OBEdgeBase*>::iterator j;
  for (i = _vatom.begin();i != _vatom.end();i++) {DestroyAtom(*i); *i = NULL;}    
  for (j = _vbond.begin();j != _vbond.end();j++) {DestroyBond(*j); *j = NULL;}

  _natoms = _nbonds = 0;

  //Delete residues
  unsigned int ii;
  for (ii=0 ; ii<_residue.size() ; ii++) {
      delete _residue[ii];
    }
  _residue.clear();

  //clear out the multiconformer data
  vector<float*>::iterator k;
  for (k = _vconf.begin();k != _vconf.end();k++)
    delete [] *k;
  _vconf.clear();

  //Clear out the pose data
  DeletePoses();

  if (!_vdata.empty()) //clean up generic data
  {
	  vector<OBGenericData*>::iterator m;
	  for (m = _vdata.begin();m != _vdata.end();m++) delete *m;
	  _vdata.clear();
  }

  _c = (float*) NULL;
  _flags = 0;
  _mod = 0;

  return(true);
}

void OBMol::BeginAccess(void)
{
  if (_access == 0) UnCompress();
  _access++;
}

void OBMol::EndAccess(void)
{
  _access--;
  if (_access == 0) Compress();
}

void OBMol::BeginModify()
{
  //suck coordinates from _c into _v for each atom
  if (!_mod && !Empty())
    {
      OBAtom *atom;
      vector<OBNodeBase*>::iterator i;
      for (atom = BeginAtom(i);atom;atom = NextAtom(i))
        {
          atom->SetVector();
          atom->ClearCoordPtr();
        }

      vector<float*>::iterator j;
      for (j = _vconf.begin();j != _vconf.end();j++) delete [] *j;

      _c = NULL;
      _vconf.clear();
      
      DeletePoses();

      //Destroy rotamer list if necessary
      if ((OBRotamerList *)GetData("RotamerList")) {
          delete (OBRotamerList *)GetData("RotamerList");
          DeleteData(obRotamerList);
        }
    }

  _mod++;
}

void OBMol::EndModify(bool nukePerceivedData)
{
  if (_mod == 0)
    {
      obAssert(false);
      ThrowError("_mod is negative - EndModify() called too many times");
      return;
    }

  _mod--;

  if (_mod) return;

  if (nukePerceivedData) _flags = 0;
  _c = NULL;

  /*
    leave generic data alone for now - just nuke it on clear()
  if (HasData("Comment")) delete [] (char*)GetData("Comment");
  _vdata.clear();
  */

  if (Empty()) return;

  //if atoms present convert coords into array
  float *c = new float [NumAtoms()*3];
  _c = c;

  int idx;
  OBAtom *atom;
  vector<OBNodeBase*>::iterator j;
  for (idx=0,atom = BeginAtom(j);atom;atom = NextAtom(j),idx++)
    {
      atom->SetIdx(idx+1);
      (atom->GetVector()).Get(&_c[idx*3]);
      atom->SetCoordPtr(&_c);
    }
  _vconf.push_back(c);

  //kekulize structure
  SetAromaticPerceived();
  Kekulize();
  UnsetAromaticPerceived();

  //    for (atom = BeginAtom(j);atom;atom = NextAtom(j))
  //      atom->UnsetAromatic();

  //    OBBond *bond;
  //    vector<OBEdgeBase*>::iterator k;
  //    for (bond = BeginBond(k);bond;bond = NextBond(k))
  //      bond->UnsetAromatic();

  //UnsetImplicitValencePerceived();
}

OBAtom *OBMol::CreateAtom(void)
{
  return new OBAtom;
}

OBBond *OBMol::CreateBond(void)
{
  return new OBBond;
}

void OBMol::DestroyAtom(OBNodeBase *atom)
{
  if (atom)
    {
      delete atom;
      atom = NULL;
    }
}

void OBMol::DestroyBond(OBEdgeBase *bond)
{
  if (bond)
    {    
      delete bond;
      bond = NULL;
    }
}

OBAtom *OBMol::NewAtom()
     //add an atom to a molecule
     //also checks bond_queue for any bonds that should be made to the new atom
{
  BeginModify();

  OBAtom *obatom = CreateAtom();
  obatom->SetIdx(_natoms+1);
  obatom->SetParent(this);


#define OBAtomIncrement 100
  if (_vatom.empty() || _natoms+1 >= (signed)_vatom.size())
    {
      _vatom.resize(_natoms+OBAtomIncrement);
      vector<OBNodeBase*>::iterator j;
      for (j = _vatom.begin(),j+=(_natoms+1);j != _vatom.end();j++) 
		  *j = (OBNodeBase*)NULL;
    }
#undef OBAtomIncrement 

  _vatom[_natoms] = obatom;
  _natoms++;

  if (HasData(obVirtualBondData))
    {
      /*add bonds that have been queued*/
	  OBVirtualBond *vb;
	  vector<OBGenericData*> verase;
	  vector<OBGenericData*>::iterator i;
	  for (i = BeginData();i != EndData();i++)
		  if ((*i)->GetDataType() == obVirtualBondData)
		  {
			  vb = (OBVirtualBond*)*i;
			  if (vb->GetBgn() > _natoms || vb->GetEnd() > _natoms) continue;
			  if (obatom->GetIdx() == vb->GetBgn() || obatom->GetIdx() == vb->GetEnd())
			  {
				  AddBond(vb->GetBgn(),vb->GetEnd(),vb->GetOrder());
				  verase.push_back(*i);
			  }
		  }

	   if (!verase.empty()) DeleteData(verase);
    }
    
  EndModify();

  return(obatom);
}

OBResidue *OBMol::NewResidue()
{
    OBResidue *obresidue = new OBResidue;
    obresidue->SetIdx(_residue.size());
    _residue.push_back(obresidue);
    return(obresidue);
}

bool OBMol::AddAtom(OBAtom &atom)
     //add an atom to a molecule
     //also checks bond_queue for any bonds that should be made to the new atom
{
  BeginModify();

    OBAtom *obatom = CreateAtom();
    *obatom = atom;
    obatom->SetIdx(_natoms+1);
    obatom->SetParent(this);


#define OBAtomIncrement 100
    if (_vatom.empty() || _natoms+1 >= (signed)_vatom.size())
      {
        _vatom.resize(_natoms+OBAtomIncrement);
        vector<OBNodeBase*>::iterator j;
        for (j = _vatom.begin(),j+=(_natoms+1);j != _vatom.end();j++) 
			*j = (OBNodeBase*)NULL;
      }
#undef OBAtomIncrement 

    _vatom[_natoms] = (OBNodeBase*)obatom;
    _natoms++;

  if (HasData(obVirtualBondData))
    {
      /*add bonds that have been queued*/
	  OBVirtualBond *vb;
	  vector<OBGenericData*> verase;
	  vector<OBGenericData*>::iterator i;
	  for (i = BeginData();i != EndData();i++)
		  if ((*i)->GetDataType() == obVirtualBondData)
		  {
			  vb = (OBVirtualBond*)*i;
			  if (vb->GetBgn() > _natoms || vb->GetEnd() > _natoms) continue;
			  if (obatom->GetIdx() == vb->GetBgn() || obatom->GetIdx() == vb->GetEnd())
			  {
				  AddBond(vb->GetBgn(),vb->GetEnd(),vb->GetOrder());
				  verase.push_back(*i);
			  }
		  }

	   if (!verase.empty()) DeleteData(verase);
    }
      
    EndModify();

    return(true);
}

bool OBMol::InsertAtom(OBAtom &atom)
{
    BeginModify();

    OBAtom *obatom = CreateAtom();
    *obatom = atom;
    obatom->SetIdx(_natoms+1);
    obatom->SetParent(this);


#define OBAtomIncrement 100
    if (_vatom.empty() || _natoms+1 >= (signed)_vatom.size())
      {
        _vatom.resize(_natoms+OBAtomIncrement);
        vector<OBNodeBase*>::iterator j;
        for (j = _vatom.begin(),j+=(_natoms+1);j != _vatom.end();j++) 
			*j = (OBNodeBase*)NULL;
      }
#undef OBAtomIncrement 

    _vatom[_natoms] = (OBNodeBase*)obatom;
    _natoms++;

  if (HasData(obVirtualBondData))
    {
      /*add bonds that have been queued*/
	  OBVirtualBond *vb;
	  vector<OBGenericData*> verase;
	  vector<OBGenericData*>::iterator i;
	  for (i = BeginData();i != EndData();i++)
		  if ((*i)->GetDataType() == obVirtualBondData)
		  {
			  vb = (OBVirtualBond*)*i;
			  if (vb->GetBgn() > _natoms || vb->GetEnd() > _natoms) continue;
			  if (obatom->GetIdx() == vb->GetBgn() || obatom->GetIdx() == vb->GetEnd())
			  {
				  AddBond(vb->GetBgn(),vb->GetEnd(),vb->GetOrder());
				  verase.push_back(*i);
			  }
		  }

	   if (!verase.empty()) DeleteData(verase);
    }
    
    EndModify();

    return(true);
}

bool OBMol::AddResidue(OBResidue &residue)
{
    BeginModify();

    OBResidue *obresidue = new OBResidue;
    *obresidue = residue;

    obresidue->SetIdx(_residue.size());

    _residue.push_back(obresidue);

    EndModify();

    return(true);
}

bool OBMol::StripSalts()
     //deletes all atoms except for the largest contiguous fragment
{
  vector<vector<int> > cfl;
  vector<vector<int> >::iterator i,max;

  ContigFragList(cfl);
  if (cfl.empty() || cfl.size() == 1)  return(false);

  max = cfl.begin();
  for (i = cfl.begin();i != cfl.end();i++)
      if ((*max).size() < (*i).size()) max = i;

  vector<int>::iterator j;
  vector<OBNodeBase*> delatoms;
  for (i = cfl.begin();i != cfl.end();i++)
      if (i != max)
          for (j = (*i).begin();j != (*i).end();j++) delatoms.push_back(GetAtom(*j));

  if (!delatoms.empty())
    {
      int tmpflags = _flags & (~(OB_SSSR_MOL));
      BeginModify();
      vector<OBNodeBase*>::iterator k;
      for (k = delatoms.begin();k != delatoms.end();k++) DeleteAtom((OBAtom*)*k);
      EndModify();
      _flags = tmpflags;
    }

  return(true);
}

bool OBMol::DeleteNonPolarHydrogens()
{
  OBAtom *atom;
  vector<OBNodeBase*>::iterator i;
  vector<OBNodeBase*> delatoms;

  for (atom = BeginAtom(i);atom;atom = NextAtom(i))
    if (atom->IsNonPolarHydrogen())
      delatoms.push_back(atom);   
  
  if (delatoms.empty()) return(true);

/*
  IncrementMod();

  int idx1,idx2;
  vector<float*>::iterator j;
  for (idx1=0,idx2=0,atom = BeginAtom(i);atom;atom = NextAtom(i),idx1++)
    if (!atom->IsHydrogen())
      {
        for (j = _vconf.begin();j != _vconf.end();j++)
           memcpy((char*)&((*j)[idx2*3]),(char*)&((*j)[idx1*3]),sizeof(float)*3);
        idx2++;
      }
*/

  for (i = delatoms.begin();i != delatoms.end();i++)  DeleteAtom((OBAtom *)*i);

  DecrementMod();

  return(true);
}

bool OBMol::DeleteHydrogens()
{
  OBAtom *atom,*nbr;
  vector<OBNodeBase*>::iterator i;
  vector<OBNodeBase*> delatoms,va;

  for (atom = BeginAtom(i);atom;atom = NextAtom(i))
    if (atom->IsHydrogen())
      delatoms.push_back(atom);
  
  if (delatoms.empty()) return(true);

/* decide whether these flags need to be reset
      _flags &= (~(OB_ATOMTYPES_MOL));
      _flags &= (~(OB_HYBRID_MOL));
      _flags &= (~(OB_PCHARGE_MOL)); 
      _flags &= (~(OB_IMPVAL_MOL));
*/

  //find bonds to delete
  vector<OBEdgeBase*> vdb;
  vector<OBEdgeBase*>::iterator j;
  for (atom = BeginAtom(i);atom;atom = NextAtom(i))
      if (!atom->IsHydrogen())
          for (nbr = atom->BeginNbrAtom(j);nbr;nbr = atom->NextNbrAtom(j))
              if (nbr->IsHydrogen())
                  vdb.push_back(*j);

  IncrementMod();
  for (j = vdb.begin();j != vdb.end();j++) DeleteBond((OBBond *)*j); //delete bonds
  DecrementMod();
   
  int idx1,idx2;
  vector<float*>::iterator k;
  for (idx1=0,idx2=0,atom = BeginAtom(i);atom;atom = NextAtom(i),idx1++)
    if (!atom->IsHydrogen())
      {
        //Update conformer coordinates
        for (k = _vconf.begin();k != _vconf.end();k++)
            memcpy((char*)&((*k)[idx2*3]),(char*)&((*k)[idx1*3]),sizeof(float)*3);

        //Update pose coordinates if necessary
        if (_xyz_pose) {
            memcpy((char*)&((_xyz_pose)[idx2*3]),(char*)&((_xyz_pose)[idx1*3]),sizeof(float)*3);
          }  

        idx2++;
        va.push_back(atom);
      }

  for (i = delatoms.begin();i != delatoms.end();i++)
  {
    DestroyAtom(*i);
    _natoms--;
  }

  _vatom.clear();
  for (i = va.begin();i != va.end();i++) _vatom.push_back((OBNodeBase*)*i);

  //_atom = va;
  //_atom.resize(_atom.size()+1);
  //_atom[_atom.size()-1] = NULL;
  _natoms = va.size();

//reset all the indices to the atoms
  for (idx1=1,atom = BeginAtom(i);atom;atom = NextAtom(i),idx1++)
      atom->SetIdx(idx1);

  return(true);
}

bool OBMol::DeleteHydrogens(OBAtom *atom)
//deletes all hydrogens attached to the atom passed to the function
{
  OBAtom *nbr;
  vector<OBNodeBase*>::iterator i;
  vector<OBEdgeBase*>::iterator k;
  vector<OBNodeBase*> delatoms;

  for (nbr = atom->BeginNbrAtom(k);nbr;nbr = atom->NextNbrAtom(k))
    if (nbr->IsHydrogen())
        delatoms.push_back(nbr);

  if (delatoms.empty())
  return(true);

  IncrementMod();
  for (i = delatoms.begin();i != delatoms.end();i++)
    DeleteHydrogen((OBAtom*)*i);
  DecrementMod();

  return(true);
}


bool OBMol::DeleteHydrogen(OBAtom *atom)
//deletes the hydrogen atom passed to the function
{
  //find bonds to delete
  OBAtom *nbr;
  vector<OBEdgeBase*> vdb;
  vector<OBEdgeBase*>::iterator j;
  for (nbr = atom->BeginNbrAtom(j);nbr;nbr = atom->NextNbrAtom(j)) vdb.push_back(*j);

  IncrementMod();
  for (j = vdb.begin();j != vdb.end();j++) DeleteBond((OBBond*)*j); //delete bonds
  DecrementMod();

  int idx;
  if (atom->GetIdx() != (int) NumAtoms())
    {
      idx = atom->GetCIdx();
      int size = NumAtoms()-atom->GetIdx();
      vector<float*>::iterator k;
      for (k = _vconf.begin();k != _vconf.end();k++)
          memmove((char*)&(*k)[idx],(char*)&(*k)[idx+3],sizeof(float)*3*size);

      //Update pose coordinates if necessary
      if (_xyz_pose) {
          memmove((char*)&(_xyz_pose)[idx],(char*)&(_xyz_pose)[idx+3],sizeof(float)*3*size);
        }  

    }
  
  _vatom.erase(_vatom.begin()+(atom->GetIdx()-1));
  DestroyAtom(atom);
  _natoms--;

//reset all the indices to the atoms
  vector<OBNodeBase*>::iterator i;
  for (idx=1,atom = BeginAtom(i);atom;atom = NextAtom(i),idx++)
    atom->SetIdx(idx);

  return(true);
}

bool OBMol::AddHydrogens(bool polaronly,bool correctForPH)
{
  if (!IsCorrectedForPH() && correctForPH) CorrectForPH();

  if (HasHydrogensAdded()) return(true);
  SetHydrogensAdded();

  //count up number of hydrogens to add
  OBAtom *atom,*h;
  int hcount,count=0;
  vector<pair<OBAtom*,int> > vhadd;
  vector<OBNodeBase*>::iterator i;
  for (atom = BeginAtom(i);atom;atom = NextAtom(i))
    {
      if (polaronly && !(atom->IsNitrogen() || atom->IsOxygen() || 
                         atom->IsSulfur() || atom->IsPhosphorus()))
      continue;

      hcount = atom->GetImplicitValence() - atom->GetValence();
      if (hcount < 0) hcount = 0;
      if (hcount)
        {
          vhadd.push_back(pair<OBAtom*,int>(atom,hcount));
          count += hcount;
        }
    }

  if (count == 0) return(true);
  bool hasCoords = HasNonZeroCoords();

  //realloc memory in coordinate arrays for new hydrogens
  float *tmpf;
  vector<float*>::iterator j;
  for (j = _vconf.begin();j != _vconf.end();j++)
    {
      tmpf = new float [(NumAtoms()+count)*3];
      memset(tmpf,'\0',sizeof(float)*(NumAtoms()+count)*3);
      if (hasCoords) memcpy(tmpf,(*j),sizeof(float)*NumAtoms()*3);
      delete []*j;
      *j = tmpf;
    }

  //Just delete the pose coordinate array.  It will automatically
  //be reallocated when SetPose(unsigned int) is called.
  if (_xyz_pose) {delete [] _xyz_pose; _xyz_pose = NULL;}

  IncrementMod();

  int m,n;
  Vector v;
  vector<pair<OBAtom*,int> >::iterator k;
  float hbrad = etab.CorrectedBondRad(1,0);


  for (k = vhadd.begin();k != vhadd.end();k++)
    {
      atom = k->first;
      float bondlen = hbrad+etab.CorrectedBondRad(atom->GetAtomicNum(),atom->GetHyb());
      for (m = 0;m < k->second;m++)
        {
          for (n = 0;n < NumConformers();n++)
            {
              SetConformer(n);
              if (hasCoords)
                {
                  atom->GetNewBondVector(v,bondlen);
                  _c[(NumAtoms())*3]   = v.x();
                  _c[(NumAtoms())*3+1] = v.y();
                  _c[(NumAtoms())*3+2] = v.z();
                }
              else memset((char*)&_c[NumAtoms()*3],'\0',sizeof(float)*3);
            }
          h = NewAtom(); h->SetType("H"); h->SetAtomicNum(1);
          AddBond(atom->GetIdx(),h->GetIdx(),1);
          h->SetCoordPtr(&_c);
        }
    }

  DecrementMod();
  SetConformer(0);

  //reset atom type and partial charge flags
  _flags &= (~(OB_PCHARGE_MOL|OB_ATOMTYPES_MOL));

  return(true);
}

bool OBMol::AddPolarHydrogens()
{
  return(AddHydrogens(true));
}

bool OBMol::AddHydrogens(OBAtom *atom)
{
  OBAtom *h;

  //count up number of hydrogens to add
  int hcount,count=0;
  vector<pair<OBAtom*,int> > vhadd;

  hcount = atom->GetImplicitValence() - atom->GetValence();

  if (hcount < 0) hcount = 0;
  if (hcount)
    {
      vhadd.push_back(pair<OBAtom*,int>(atom,hcount));
      count += hcount;
    }

  if (count == 0) return(true);

  //realloc memory in coordinate arrays for new hydroges
  float *tmpf;
  vector<float*>::iterator j;
  for (j = _vconf.begin();j != _vconf.end();j++)
    {
      tmpf = new float [(NumAtoms()+count)*3+10];
      memcpy(tmpf,(*j),sizeof(float)*NumAtoms()*3);
      delete []*j;
      *j = tmpf;
    }

  //Just delete the pose coordinate array.  It will automatically
  //be reallocated when SetPose(unsigned int) is called.
  if (_xyz_pose) {delete [] _xyz_pose; _xyz_pose = NULL;}

  IncrementMod();

  int m,n;
  Vector v;
  vector<pair<OBAtom*,int> >::iterator k;
  float hbrad = etab.CorrectedBondRad(1,0);

  for (k = vhadd.begin();k != vhadd.end();k++)
    {
      atom = k->first;
      float bondlen = hbrad+etab.CorrectedBondRad(atom->GetAtomicNum(),atom->GetHyb());
      for (m = 0;m < k->second;m++)
        {
          for (n = 0;n < NumConformers();n++)
            {
              SetConformer(n);
              atom->GetNewBondVector(v,bondlen);
              _c[(NumAtoms())*3]   = v.x();
              _c[(NumAtoms())*3+1] = v.y();
              _c[(NumAtoms())*3+2] = v.z();
            }
          h = NewAtom(); h->SetType("H"); h->SetAtomicNum(1);
          AddBond(atom->GetIdx(),h->GetIdx(),1);
          h->SetCoordPtr(&_c);
        }
    }

  DecrementMod();
  SetConformer(0);

  //reset atom type and partial charge flags
  //_flags &= (~(OB_PCHARGE_MOL|OB_ATOMTYPES_MOL));

  return(true);
}

bool OBMol::CorrectForPH()
{
  if (IsCorrectedForPH()) return(true);
  phmodel.CorrectForPH(*this);

  return(true);
}

static void ResetVisit(OBMol &mol,vector<int> &visit,int depth)
{
  OBBond *bond;
  vector<OBEdgeBase*>::iterator i;
  
  for (bond = mol.BeginBond(i);bond;bond = mol.NextBond(i))
    if (bond->IsAromatic() && visit[bond->GetIdx()] >= depth)
      visit[bond->GetIdx()] = 0;
}

static int ValenceSum(OBAtom *atom)
{
  int count = atom->GetImplicitValence();

  OBBond *bond;
  vector<OBEdgeBase*>::iterator i;
  for (bond = atom->BeginBond(i);bond;bond = atom->NextBond(i))
    if (bond->IsKDouble())
      count++;

  return(count);
}

static bool KekulePropagate(OBAtom *atom,vector<int> &visit,vector<int> &ival,int depth)
{
  int count = 0;
  OBBond *bond;
  vector<OBEdgeBase*>::iterator i;
  for (bond = atom->BeginBond(i);bond;bond = atom->NextBond(i))
    if (!visit[bond->GetIdx()])
      count++;

  if (!count)
    return(ValenceSum(atom) == ival[atom->GetIdx()]);

  bool result = true;
  OBAtom *nbr;

  if (ValenceSum(atom) >= ival[atom->GetIdx()])
    {
      for (nbr = atom->BeginNbrAtom(i);nbr;nbr = atom->NextNbrAtom(i))
          if (nbr->IsAromatic() && !visit[(*i)->GetIdx()])
            {
              visit[(*i)->GetIdx()] = depth;
              ((OBBond*)*i)->SetKSingle();
              result = KekulePropagate(nbr,visit,ival,depth);
              if (result) break;
//            if (!result) break;
            }
    }
  else if (count == 1)
    for (nbr = atom->BeginNbrAtom(i);nbr;nbr = atom->NextNbrAtom(i))
      if (nbr->IsAromatic() && !visit[(*i)->GetIdx()])
        {
          visit[(*i)->GetIdx()] = depth;
          ((OBBond*)*i)->SetKDouble();
          result = KekulePropagate(nbr,visit,ival,depth);
		  //break;
          if (result) break;
        }
  return(result);
}

int GetCurrentValence(OBAtom *atom)
{
  int count = atom->GetImplicitValence();

  OBBond *bond;
  vector<OBEdgeBase*>::iterator i;
  for (bond = atom->BeginBond(i);bond;bond = atom->NextBond(i))
  if (bond->IsKDouble()) count++;
  else if (bond->IsKTriple()) count += 2;
  return(count);
}

bool ExpandKekule(OBMol &mol, vector<OBNodeBase*> &va,
                  vector<OBNodeBase*>::iterator i,
		  vector<int> &maxv,bool secondpass)
{
	if (i == va.end())
    {
      //check to see that the ideal valence has been achieved for all atoms
      vector<OBNodeBase*>::iterator j;
      for (j = va.begin();j != va.end();j++)
        {
          //let erroneously aromatic carboxylates pass
          if (((OBAtom*)*j)->IsOxygen() && ((OBAtom*)*j)->GetValence() == 1)
	    continue;
          if (GetCurrentValence((OBAtom*)*j) != maxv[(*j)->GetIdx()])
	    return(false);
        }
      return(true);
    }

  //jump to next atom in list if current atom doesn't have any attached
  //aromatic bonds
  OBBond *bond;
  OBAtom *atom = (OBAtom*)*i;
  vector<OBEdgeBase*>::iterator j;
  bool done = true;
  for (bond = atom->BeginBond(j);bond;bond = atom->NextBond(j))
      if (bond->GetBO() == 5)
        {
          done = false;
          break;
        }
  if (done)
	  return(ExpandKekule(mol,va,i+1,maxv,secondpass));

  //store list of attached aromatic atoms
  OBAtom *nbr;
  vector<OBEdgeBase*> vb;
  for (nbr = atom->BeginNbrAtom(j);nbr;nbr = atom->NextNbrAtom(j))
      if ((*j)->GetBO() == 5)
        {
          vb.push_back(*j);
          ((OBBond *)*j)->SetBO(1);
          ((OBBond *)*j)->SetKSingle();
        }

  //try setting a double bond
  if (GetCurrentValence(atom) < maxv[atom->GetIdx()]) 
    {
	  for (j = vb.begin();j != vb.end();j++)
        {
          nbr = ((OBBond *)*j)->GetNbrAtom(atom);
          if (GetCurrentValence(nbr) <= maxv[nbr->GetIdx()])
            {
              ((OBBond*)*j)->SetKDouble();
              ((OBBond*)*j)->SetBO(2);
              if (ExpandKekule(mol,va,i+1,maxv,secondpass)) return(true);
              ((OBBond*)*j)->SetKSingle();
              ((OBBond*)*j)->SetBO(1);
            }
        }

      if (secondpass && atom->IsNitrogen() && atom->GetFormalCharge() == 0 && 
		  atom->GetImplicitValence() == 2)
        {
		  atom->IncrementImplicitValence();
		  if (ExpandKekule(mol,va,i+1,maxv,secondpass)) return(true);
		  atom->DecrementImplicitValence();
        }
    }
  else  //full valence - no double bond to set
    {
      if (ExpandKekule(mol,va,i+1,maxv,secondpass)) return(true);

	  bool trycharge = false;
	  if (secondpass && atom->GetFormalCharge() == 0)
	  {
		  if (atom->IsNitrogen() && atom->GetHvyValence() == 3) trycharge = true;
		  if (atom->IsOxygen() && atom->GetHvyValence() == 2) trycharge = true;
		  if (atom->IsSulfur() && atom->GetHvyValence() == 2) trycharge = true;
	  }

	  if (trycharge) //attempt to charge up O,N,S to make a valid kekule form
	  {
		  maxv[atom->GetIdx()]++;
		  atom->SetFormalCharge(1);
	  	  for (j = vb.begin();j != vb.end();j++)
		  {
			  nbr = ((OBBond*)*j)->GetNbrAtom(atom);
			  if (GetCurrentValence(nbr) <= maxv[nbr->GetIdx()])
			  {
				  ((OBBond*)*j)->SetKDouble();
				  ((OBBond*)*j)->SetBO(2);
				  if (ExpandKekule(mol,va,i+1,maxv,secondpass)) return(true);
				  ((OBBond*)*j)->SetKSingle();
				  ((OBBond*)*j)->SetBO(1);
			  }
		  }
		  maxv[atom->GetIdx()]--;
		  atom->SetFormalCharge(0);
	  }

      if (secondpass && atom->IsNitrogen() && atom->GetFormalCharge() == 0 && 
		  atom->GetImplicitValence() == 2) //try protonating the nitrogen
        {
		  atom->IncrementImplicitValence();
		  if (ExpandKekule(mol,va,i+1,maxv,secondpass)) return(true);
		  atom->DecrementImplicitValence();
        }
    }

  //failed to find a valid solution - reset attached bonds
  for (j = vb.begin();j != vb.end();j++) 
    {
      ((OBBond*)*j)->SetKSingle();
      ((OBBond*)*j)->SetBO(5);
    }

  return(false);
}

void CorrectBadResonanceForm(OBMol &mol)
{
  string s;
  OBBond *b1,*b2,*b3;
  OBAtom *a1,*a2,*a3,*a4;
  vector<vector<int> > mlist;
  vector<vector<int> >::iterator i;

  OBSmartsPattern acid; acid.Init("[oD1]c[oD1]");

  //carboxylic acid
  if (acid.Match(mol))
    {
      mlist = acid.GetUMapList();
      for (i = mlist.begin();i != mlist.end();i++)
        {
          a1 = mol.GetAtom((*i)[0]);
          a2 = mol.GetAtom((*i)[1]);
          a3 = mol.GetAtom((*i)[2]);
          b1 = a2->GetBond(a1); b2 = a2->GetBond(a3);
          if (!b1 || !b2) continue;
          b1->SetKDouble();
          b2->SetKSingle();
        }
    }

  //phosphonic acid
  OBSmartsPattern phosphate; phosphate.Init("[p]([oD1])([oD1])([oD1])[#6,#8]");
  if (phosphate.Match(mol))
    {
      mlist = phosphate.GetUMapList();
      for (i = mlist.begin();i != mlist.end();i++)
        {
          a1 = mol.GetAtom((*i)[0]);
          a2 = mol.GetAtom((*i)[1]);
          a3 = mol.GetAtom((*i)[2]);
          a4 = mol.GetAtom((*i)[3]);
          b1 = a1->GetBond(a2); b2 = a1->GetBond(a3); b3 = a1->GetBond(a4);

          if (!b1 || !b2 || !b3) continue;
          b1->SetKDouble();
          b2->SetKSingle();
          b3->SetKSingle();
        }
    }

  //amidene and guanidine
  OBSmartsPattern amidene; amidene.Init("[nD1]c([nD1])*");
  if (amidene.Match(mol))
    {
      mlist = amidene.GetUMapList();
      for (i = mlist.begin();i != mlist.end();i++)
        {
          a1 = mol.GetAtom((*i)[0]);
          a2 = mol.GetAtom((*i)[1]);
          a3 = mol.GetAtom((*i)[2]);
          b1 = a2->GetBond(a1); b2 = a2->GetBond(a3);
          if (!b1 || !b2) continue;
          b1->SetKDouble();
          b2->SetKSingle();
        }
    }
}

bool OBMol::PerceiveKekuleBonds()
{
  if (HasKekulePerceived())  return(true);
  SetKekulePerceived();

  OBBond *bond;
  vector<OBEdgeBase*>::iterator i;

  //initialize kekule bonds
  bool done = true;
  bool badResonanceForm = false;
  vector<bool> varo; varo.resize(NumAtoms()+1,false);
  for (bond = BeginBond(i);bond;bond = NextBond(i))
  switch (bond->GetBO())
    {
      case 2: bond->SetKDouble(); break;
      case 3: bond->SetKTriple(); break;
      case 5: 

      bond->SetKSingle(); 
      if (bond->IsInRing())
        {
          varo[bond->GetBeginAtomIdx()] = true;
          varo[bond->GetEndAtomIdx()]   = true;
          done = false;
        }
      else badResonanceForm = true;

      break;

      default: bond->SetKSingle(); break;
    }

  if (badResonanceForm) CorrectBadResonanceForm(*this);

  if (done) return(true);

  //set the maximum valence for each aromatic atom
  OBAtom *atom,*nbr;
  vector<OBNodeBase*>::iterator j,k;
  vector<int> maxv; maxv.resize(NumAtoms()+1);

  for (atom = BeginAtom(j);atom;atom = NextAtom(j))
  if (varo[atom->GetIdx()])
    {
      switch (atom->GetAtomicNum())
        {
          case 6:  
            maxv[atom->GetIdx()] = 4; 
          break;
          case 8:  
          case 16: 
          case 34: 
          case 52: maxv[atom->GetIdx()] = 2; break;
          case 7:  
          case 15: 
          case 33: maxv[atom->GetIdx()] = 3; break;
        }
      //correct valence for formal charges
      if (atom->IsCarbon()) maxv[atom->GetIdx()] -= abs(atom->GetFormalCharge());
      else                  maxv[atom->GetIdx()] += atom->GetFormalCharge();

      if (atom->IsNitrogen() || atom->IsSulfur())
          for (nbr = atom->BeginNbrAtom(i);nbr;nbr = atom->NextNbrAtom(i))
              if (nbr->IsOxygen() && (*i)->GetBO() == 2)
                  maxv[atom->GetIdx()] += 2;
  }

  bool result = true;
  vector<bool> used; used.resize(NumAtoms()+1);
  vector<OBNodeBase*> va,curr,next;
  for (atom = BeginAtom(j);atom;atom = NextAtom(j))
      if (varo[atom->GetIdx()] && !used[atom->GetIdx()])
        {
          va.clear(); va.push_back(atom);
          curr.clear();
          curr.push_back(atom);
          used[atom->GetIdx()] = true;

          for (;!curr.empty();)
            {
              next.clear();
              for (k = curr.begin();k != curr.end();k++)
                  for (nbr = ((OBAtom*)*k)->BeginNbrAtom(i);nbr;nbr = ((OBAtom*)*k)->NextNbrAtom(i))
                      if (varo[nbr->GetIdx()] && !used[nbr->GetIdx()])
                        {
                          used[nbr->GetIdx()] = true;
                          next.push_back(nbr);
                          va.push_back(nbr);
                        }
              curr = next;
            }

          //try it first without protonating aromatic nitrogens
          if (!ExpandKekule(*this,va,va.begin(),maxv,false) && 
			  !ExpandKekule(*this,va,va.begin(),maxv,true)) 
			  result = false;
        }

  if (!result)
    { 
      cerr << "Kekulization Error = " << GetTitle() << endl;
      //exit(0);
    }
  
  return(result);
}

bool OBMol::Kekulize()
{
  OBBond *bond;
  vector<OBEdgeBase*>::iterator i;
  if (NumAtoms() > 255) return(false);

  for (bond = BeginBond(i);bond;bond = NextBond(i))
    if (bond->IsKSingle())      bond->SetBO(1);
    else if (bond->IsKDouble()) bond->SetBO(2);
    else if (bond->IsKTriple()) bond->SetBO(3);

  return(true);
}

bool OBMol::DeleteAtom(OBAtom *atom)
{
  if (atom->IsHydrogen()) return(DeleteHydrogen(atom));

  BeginModify(); 
  //don't need to do anything with coordinates b/c 
  //BeginModify() blows away coordinates

  //find bonds to delete
  OBAtom *nbr;
  vector<OBEdgeBase*> vdb;
  vector<OBEdgeBase*>::iterator j;
  for (nbr = atom->BeginNbrAtom(j);nbr;nbr = atom->NextNbrAtom(j))
      vdb.push_back(*j);

  for (j = vdb.begin();j != vdb.end();j++) DeleteBond((OBBond *)*j); //delete bonds

  _vatom.erase(_vatom.begin()+(atom->GetIdx()-1));
  DestroyAtom(atom);
  _natoms--;

  //reset all the indices to the atoms
  int idx;
  vector<OBNodeBase*>::iterator i;
  for (idx=1,atom = BeginAtom(i);atom;atom = NextAtom(i),idx++)
    atom->SetIdx(idx);

  EndModify();

  return(true);
}

bool OBMol::DeleteResidue(OBResidue *residue)
{
    unsigned short idx = residue->GetIdx();
    for ( unsigned short i = idx ; i < _residue.size() ; i++ )
        _residue[i]->SetIdx(i-1);

    _residue.erase(_residue.begin() + idx);

    if (residue)
      {
	delete residue;
	residue = NULL;
      }

    return(true);
}

bool OBMol::AddBond(int first,int second,int order,int stereo,int insertpos)
{
  BeginModify();

  if ((unsigned)first <= NumAtoms() && (unsigned)second <= NumAtoms()
      && !GetBond(first, second))
    //atoms exist and bond doesn't
    {
      OBBond *bond = CreateBond();
      if (!bond)
        {
          EndModify();
          return(false);
        }

      OBAtom *bgn,*end;
      bgn = GetAtom(first);
      end = GetAtom(second);
      if (!bgn || !end)
        {
          ThrowError("Unable to add bond - invalid atom index");
          return(false);
        }
      bond->Set(_nbonds,bgn,end,order,stereo);
      bond->SetParent(this);

      //set aromatic flags if it has the appropriate order
      if (order == 5)
        {
          bond->SetAromatic();
          bgn->SetAromatic();
          end->SetAromatic();
        }

#define OBBondIncrement 100
      if (_vbond.empty() || _nbonds+1 >= (signed)_vbond.size()) 
        {
          _vbond.resize(_nbonds+OBBondIncrement);
          vector<OBEdgeBase*>::iterator i;
          for (i = _vbond.begin(),i+=(_nbonds+1);i != _vbond.end();i++) 
			  *i = (OBEdgeBase*)NULL;
        }
#undef  OBBondIncrement

      _vbond[_nbonds] = (OBEdgeBase*)bond;
      _nbonds++;

      if (insertpos == -1)
        {
          bgn->AddBond(bond);
          end->AddBond(bond);
        }
      else
        {
          if (insertpos >= bgn->GetValence()) bgn->AddBond(bond);
          else //need to insert the bond for the connectivity order to be preserved
            {    //otherwise stereochemistry gets screwed up
              vector<OBEdgeBase*>::iterator bi;
              bgn->BeginNbrAtom(bi);
              bi += insertpos;
              bgn->InsertBond(bi,bond);
            }
          end->AddBond(bond);
        }
    }
  else //at least one atom doesn't exist yet - add to bond_q
	  SetData(new OBVirtualBond(first,second,order,stereo));

  EndModify();
  return(true);
}

bool OBMol::AddBond(OBBond &bond)
{
  return(AddBond(bond.GetBeginAtomIdx(),
                 bond.GetEndAtomIdx(),
                 bond.GetBO(),
                 bond.GetFlags()));
}

bool OBMol::DeleteBond(OBBond *bond)
{
  BeginModify();

  (bond->GetBeginAtom())->DeleteBond(bond);
  (bond->GetEndAtom())->DeleteBond(bond);
  _vbond.erase(_vbond.begin() + bond->GetIdx());

  DestroyBond(bond);

  vector<OBEdgeBase*>::iterator i; int j;
  for (bond = BeginBond(i),j=0;bond;bond = NextBond(i),j++)
    bond->SetIdx(j);

  _nbonds--;
  EndModify();
  return(true);
}

void OBMol::Align(OBAtom *a1,OBAtom *a2,Vector &p1,Vector &p2)
  // aligns atom a on p1 and atom b along p1->p2 vector
{
  vector<int> children;
  
  //find which atoms to rotate
  FindChildren(children,a1->GetIdx(),a2->GetIdx());
  children.push_back(a2->GetIdx());

  //find the rotation vector and angle
  Vector v1,v2,v3;
  v1 = p2 - p1;
  v2 = a2->GetVector() - a1->GetVector();
  v3 = cross(v1,v2);
  float angle = VectorAngle(v1,v2);

  //find the rotation matrix
  Matrix3x3 m;
  m.RotAboutAxisByAngle(v3,angle);

  //rotate atoms
  Vector v;
  OBAtom *atom;
  vector<int>::iterator i;
  for (i = children.begin();i != children.end();i++)
    {
      atom = GetAtom(*i);
      v = atom->GetVector();
      v -= a1->GetVector();
      v *= m;   //rotate the point
      v += p1;  //translate the vector
      atom->SetVector(v);
    }
  //set a1 = p1
  a1->SetVector(p1); 
}


void jacobi3x3(float a[3][3], float v[3][3])
{
#define MAX_SWEEPS 50
  float onorm, dnorm;
  float b, dma, q, t, c, s,d[3];
  float atemp, vtemp, dtemp;
  int i, j, k, l;

  memset((char*)d,'\0',sizeof(float)*3);
  
  for (j = 0; j < 3; j++) 
    {
      for (i = 0; i < 3; i++) v[i][j] = 0.0;

      v[j][j] = 1.0;
      d[j] = a[j][j];
    }
  
  for (l = 1; l <= MAX_SWEEPS; l++) 
    {
      dnorm = 0.0;
      onorm = 0.0;
      for (j = 0; j < 3; j++) 
        {
          dnorm = dnorm + fabs(d[j]);
          for (i = 0; i <= j - 1; i++) 
            {
              onorm = onorm + fabs(a[i][j]);
            }
        }
      
      if((onorm/dnorm) <= 1.0e-12) goto Exit_now;
      for (j = 1; j < 3; j++) 
        {
          for (i = 0; i <= j - 1; i++) 
            {
              b = a[i][j];
              if(fabs(b) > 0.0) 
                {
                  dma = d[j] - d[i];
                  if((fabs(dma) + fabs(b)) <=  fabs(dma)) 
                    {
                      t = b / dma;
                    }
                  else 
                    {
                      q = 0.5 * dma / b;
                      t = 1.0/(fabs(q) + sqrt(1.0+q*q));
                      if(q < 0.0) 
                        {
                          t = -t;
                        }
                    }
                  c = 1.0/sqrt(t * t + 1.0);
                  s = t * c;
                  a[i][j] = 0.0;
                  for (k = 0; k <= i-1; k++) 
                    {
                      atemp = c * a[k][i] - s * a[k][j];
                      a[k][j] = s * a[k][i] + c * a[k][j];
                      a[k][i] = atemp;
                    }
                  for (k = i+1; k <= j-1; k++) 
                    {
                      atemp = c * a[i][k] - s * a[k][j];
                      a[k][j] = s * a[i][k] + c * a[k][j];
                      a[i][k] = atemp;
                    }
                  for (k = j+1; k < 3; k++) 
                    {
                      atemp = c * a[i][k] - s * a[j][k];
                      a[j][k] = s * a[i][k] + c * a[j][k];
                      a[i][k] = atemp;
                    }
                  for (k = 0; k < 3; k++) 
                   {
                      vtemp = c * v[k][i] - s * v[k][j];
                      v[k][j] = s * v[k][i] + c * v[k][j];
                      v[k][i] = vtemp;
                   }
                  dtemp = c*c*d[i] + s*s*d[j] - 2.0*c*s*b;
                  d[j] = s*s*d[i] + c*c*d[j] +  2.0*c*s*b;
                  d[i] = dtemp;
                }  /* end if */
            } /* end for i */
        } /* end for j */
    } /* end for l */
  
Exit_now:
  
  /* max_sweeps = l;*/
  
  for (j = 0; j < 3-1; j++) 
    {
      k = j;
      dtemp = d[k];
      for (i = j+1; i < 3; i++) 
        {
          if(d[i] < dtemp) 
           {
             k = i;
             dtemp = d[k];
           }
        }
      
      if(k > j) 
        {
          d[k] = d[j];
          d[j] = dtemp;
          for (i = 0; i < 3 ; i++) 
            {
              dtemp = v[i][k];
              v[i][k] = v[i][j];
              v[i][j] = dtemp;
            }
      }
    }
#undef MAX_SWEEPS
}

void OBMol::ToInertialFrame()
{
  float m[9];
  for (int i = 0;i < NumConformers();i++)
    ToInertialFrame(i,m);
}

void OBMol::ToInertialFrame(int conf,float *rmat)
{
  unsigned int i;
  int count=0;
  float x,y,z;
  float center[3],m[3][3];

  for (i = 0;i < 3;i++) memset(&m[i],'\0',sizeof(float)*3);
  memset(center,'\0',sizeof(float)*3);

  SetConformer(conf);
  OBAtom *atom;
  vector<OBNodeBase*>::iterator j;
  //find center of mass
  for (atom = BeginAtom(j);atom;atom = NextAtom(j))
    if (!atom->IsHydrogen())
      {center[0] += atom->x(); center[1] += atom->y();center[2] += atom->z();count++;}

  center[0] /= (float)count; center[1] /= (float)count; center[2] /= (float)count;

  //calculate inertial tensor
  for (atom = BeginAtom(j);atom;atom = NextAtom(j))
    if (!atom->IsHydrogen())
    {
      x = atom->x()-center[0]; y = atom->y()-center[1]; z = atom->z()-center[2];

      m[0][0] += y*y + z*z; m[0][1] -= x*y;       m[0][2] -= x*z;       
      m[1][0] -= x*y;       m[1][1] += x*x + z*z; m[1][2] -= y*z;       
      m[2][0] -= x*z;       m[2][1] -= y*z;       m[2][2] += x*x + y*y; 
    }

  /* find rotation matrix for moment of inertia */
  ob_make_rmat(m,rmat);

  /*
  float v[3][3];
  jacobi3x3(m,v);

  Vector v1,v2,v3,r1,r2;
  r1.Set(v[0][0],v[1][0],v[2][0]);
  r2.Set(v[0][1],v[1][1],v[2][1]);
  
  v3 = cross(r1,r2); v3 = v3.normalize();
  v2 = cross(v3,r1); v2 = v2.normalize();
  v1 = cross(v2,v3); v1 = v1.normalize();

  float rmat[9]; 
  rmat[0] = v1.x(); rmat[1] = v1.y(); rmat[2] = v1.z();
  rmat[3] = v2.x(); rmat[4] = v2.y(); rmat[5] = v2.z();
  rmat[6] = v3.x(); rmat[7] = v3.y(); rmat[8] = v3.z();
  */

  /* rotate all coordinates */
  float *c = GetConformer(conf);
  for(i=0; i < NumAtoms();i++)
    {
      x = c[i*3]-center[0]; y = c[i*3+1]-center[1]; z = c[i*3+2]-center[2];
      c[i*3]   = x*rmat[0] + y*rmat[1] + z*rmat[2];
      c[i*3+1] = x*rmat[3] + y*rmat[4] + z*rmat[5];
      c[i*3+2] = x*rmat[6] + y*rmat[7] + z*rmat[8];
  }
}

OBMol::OBMol(io_type itype,io_type otype)
{
  _natoms = _nbonds = 0;
  _mod = 0;
  _access = 0;
  _energy = 0.0f;
  _itype = itype;
  _otype = otype;
  _vatom.clear();
  _vbond.clear();
  _vdata.clear();
  _title = "";
  _c = (float*)NULL;
  _flags = 0;
  _vconf.clear();
  _xyz_pose = NULL;
  _pose.clear();
  _cur_pose_idx=0;
  _autoPartialCharge = true;
  _autoFormalCharge = true;
  _compressed = false;
}

OBMol::OBMol(const OBMol &mol)
{
  _natoms = _nbonds = 0;
  _mod = 0;
  _access = 0;
  _vatom.clear();
  _vbond.clear();
  _vdata.clear();
  _title = "";
  _c = (float*)NULL;
  _flags = 0;
  _vconf.clear();
  _xyz_pose = NULL;
  _pose.clear();
  _cur_pose_idx=0;
  _autoPartialCharge = true;
  _autoFormalCharge = true;
  _compressed = false;
  *this = mol;
}

OBMol::~OBMol()
{
  OBAtom    *atom;
  OBBond    *bond;
  OBResidue *residue;
  vector<OBNodeBase*>::iterator i;
  vector<OBEdgeBase*>::iterator j;
  vector<OBResidue*>::iterator r;
  for (atom = BeginAtom(i);atom;atom = NextAtom(i)) DestroyAtom(atom);
  for (bond = BeginBond(j);bond;bond = NextBond(j)) DestroyBond(bond);
  for (residue = BeginResidue(r);residue;residue = NextResidue(r)) delete residue;

  //clear out the multiconformer data
  vector<float*>::iterator k;
  for (k = _vconf.begin();k != _vconf.end();k++)
    delete [] *k;
  _vconf.clear();

  //Deallocate the pose coordinate array if necessary
  if (_xyz_pose) {delete [] _xyz_pose; _xyz_pose = NULL;}

  if (!_vdata.empty())
  {
    vector<OBGenericData*>::iterator m;
    for (m = _vdata.begin();m != _vdata.end();m++) delete *m;
      _vdata.clear();
  }
}

bool OBMol::HasData(string &s)
     //returns true if the generic attribute/value pair exists
{
  if (_vdata.empty()) return(false);

    vector<OBGenericData*>::iterator i;

    for (i = _vdata.begin();i != _vdata.end();i++)
        if ((*i)->GetAttribute() == s)
            return(true);
    
    return(false);
}

bool OBMol::HasData(const char *s)
     //returns true if the generic attribute/value pair exists
{
  if (_vdata.empty()) return(false);

    vector<OBGenericData*>::iterator i;

    for (i = _vdata.begin();i != _vdata.end();i++)
        if ((*i)->GetAttribute() == s)
            return(true);
    
    return(false);
}


bool OBMol::HasData(obDataType dt)
     //returns true if the generic attribute/value pair exists
{
  if (_vdata.empty()) return(false);

    vector<OBGenericData*>::iterator i;

    for (i = _vdata.begin();i != _vdata.end();i++)
        if ((*i)->GetDataType() == dt)
            return(true);
    
    return(false);
}

OBGenericData *OBMol::GetData(string &s)
     //returns the value given an attribute
{
    vector<OBGenericData*>::iterator i;

    for (i = _vdata.begin();i != _vdata.end();i++)
		if ((*i)->GetAttribute() == s)
            return(*i);

    return(NULL);
}

OBGenericData *OBMol::GetData(const char *s)
     //returns the value given an attribute
{
    vector<OBGenericData*>::iterator i;

    for (i = _vdata.begin();i != _vdata.end();i++)
		if ((*i)->GetAttribute() == s)
            return(*i);

    return(NULL);
}

OBGenericData *OBMol::GetData(obDataType dt)
{
    vector<OBGenericData*>::iterator i;
    for (i = _vdata.begin();i != _vdata.end();i++)
        if ((*i)->GetDataType() == dt)
            return(*i);
    return(NULL);
}

void OBMol::DeleteData(obDataType dt)
{
  vector<OBGenericData*> vdata;
  vector<OBGenericData*>::iterator i;
    for (i = _vdata.begin();i != _vdata.end();i++)
        if ((*i)->GetDataType() == dt) delete *i;
        else vdata.push_back(*i);
  _vdata = vdata;
}

void OBMol::DeleteData(vector<OBGenericData*> &vg)
{
  vector<OBGenericData*> vdata;
  vector<OBGenericData*>::iterator i,j;

  bool del;
  for (i = _vdata.begin();i != _vdata.end();i++)
  {
	  del = false;
	  for (j = vg.begin();j != vg.end();j++)
		  if (*i == *j)
		  {
			  del = true;
			  break;
		  }
	   if (del) delete *i;
	   else     vdata.push_back(*i);
  }
  _vdata = vdata;
}

void OBMol::DeleteData(OBGenericData *gd)
{
  vector<OBGenericData*>::iterator i;
  for (i = _vdata.begin();i != _vdata.end();i++)
	  if (*i == gd)
	  {
		delete *i;
		_vdata.erase(i);
	  }

}

bool OBMol::HasNonZeroCoords()
{
  OBAtom *atom;
  vector<OBNodeBase*>::iterator i;
  
  for (atom = BeginAtom(i);atom;atom = NextAtom(i))
    if (atom->GetVector() != VZero)
      return(true);

    return(false);
}

bool OBMol::Has2D()
{
  bool hasX,hasY;
  OBAtom *atom;
  vector<OBNodeBase*>::iterator i;

  hasX = hasY = false;
  for (atom = BeginAtom(i);atom;atom = NextAtom(i))
    {
      if (!hasX && atom->x() != 0.0) hasX = true;
      if (!hasY && atom->y() != 0.0) hasY = true;

      if (hasX && hasY) return(true);
    }
  return(false);
}

bool OBMol::Has3D()
{
  bool hasX,hasY,hasZ;
  OBAtom *atom;
  vector<OBNodeBase*>::iterator i;

  hasX = hasY = hasZ = false;
  if (this->_c == NULL) return(false);
  for (atom = BeginAtom(i);atom;atom = NextAtom(i))
    {
      if (!hasX && atom->x() != 0.0) hasX = true;
      if (!hasY && atom->y() != 0.0) hasY = true;
      if (!hasZ && atom->z() != 0.0) hasZ = true;

      if (hasX && hasY && hasZ) return(true);
    }
  return(false);
}

bool OBMol::IsChiral()
{
  OBAtom *atom;
  vector<OBNodeBase*>::iterator i;

  for (atom = BeginAtom(i);atom;atom = NextAtom(i))
    if ((atom->IsCarbon() || atom->IsNitrogen()) && atom->GetHvyValence() > 2 && atom->IsChiral()) 
        return(true);

  return(false);
}

void OBMol::RenumberAtoms(vector<OBNodeBase*> &v)
{
  if (Empty()) return;

  OBAtom *atom;
  vector<OBNodeBase*> va;
  vector<OBNodeBase*>::iterator i;
  
  va = v;

  if (!va.empty() && va.size() < NumAtoms())
    //make sure all atoms are represented in the vector
    {
      OBBitVec bv;
      for (i = va.begin();i != va.end();i++) bv |= (*i)->GetIdx();
     
      for (atom = BeginAtom(i);atom;atom = NextAtom(i))
          if (!bv[atom->GetIdx()])
              va.push_back(atom);
    }

  int j,k;
  float *c;
  float *ctmp = new float [NumAtoms()*3];

  for (j = 0;j < NumConformers();j++)
    {
      c = GetConformer(j);
      for (k=0,i = va.begin();i != va.end();i++,k++)
          memcpy((char*)&ctmp[k*3],(char*)&c[((OBAtom*)*i)->GetCIdx()],sizeof(float)*3);
      memcpy((char*)c,(char*)ctmp,sizeof(float)*3*NumAtoms());
    }

  for (k=1,i = va.begin();i != va.end();i++,k++) (*i)->SetIdx(k);

  delete [] ctmp;

  _vatom.clear();
  for (i = va.begin();i != va.end();i++) _vatom.push_back(*i);
}

#ifdef REMOVE_LATER
bool CompareBonds(const OBEdgeBase *a,const OBEdgeBase *b) 
{
  if (a->GetBgn()->GetIdx() == b->GetBgn()->GetIdx())
    return(a->GetEnd()->GetIdx() < b->GetEnd()->GetIdx());

  //return((a->GetBgn())->GetIdx() < (b->GetBgn())->GetIdx());
}
#endif

bool WriteTitles(ostream &ofs, OBMol &mol)
{
	ofs << mol.GetTitle() << endl;
	return true;
}

/*
void OBMol::ConnectTheDots(void)
  //use inter-atomic distances to identify bonds
{
  if (Empty()) return;

  int j,k;
  OBAtom *atom,*nbr;
  vector<OBNodeBase*>::iterator i;
  vector<float> rad; 

  float *c = new float [NumAtoms()*3];
  rad.resize(_natoms);

  for (j=0,atom = BeginAtom(i);atom;atom = NextAtom(i),j++)
    {
      (atom->GetVector()).Get(&c[j*3]);
      rad[j] = etab.CorrectedBondRad(atom->GetAtomicNum(),3);
      rad[j] *= 1.10f;
    }

  float d2,cutoff;
  for (j = 0;j < _natoms;j++)
    for (k = j+1;k < _natoms;k++)
      {
        cutoff = rad[j]+rad[k]; 
        cutoff *= cutoff;

        d2 = SQUARE(c[j*3]-c[k*3]);
        d2 += SQUARE(c[j*3+1]-c[k*3+1]);
        d2 += SQUARE(c[j*3+2]-c[k*3+2]);
        if (d2 > cutoff) continue;
        atom = GetAtom(j+1); nbr = GetAtom(k+1);
        if (atom->IsConnected(nbr)) continue;
        if (atom->IsHydrogen() && nbr->IsHydrogen()) continue;
        AddBond(j+1,k+1,1);
      }

  delete [] c;
}
*/
void OBMol::ConnectTheDots(void)
  //use inter-atomic distances to identify bonds
{
  if (Empty()) return;

  int j,k,max;
  bool unset = false;
  OBAtom *atom,*nbr;
  vector<OBNodeBase*>::iterator i;
  vector<pair<OBAtom*,float> > zsortedAtoms;
  vector<float> rad; 
  vector<int> zsorted;

  float *c = new float [NumAtoms()*3];
  rad.resize(_natoms);

  for (j = 0, atom = BeginAtom(i) ; atom ; atom = NextAtom(i), j++)
    {
      (atom->GetVector()).Get(&c[j*3]);
      pair<OBAtom*,float> entry(atom, atom->GetVector().z());
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
  float d2,cutoff,zd;
  for (j = 0 ; j < max ; j++)
  {
    idx1 = zsorted[j];
    for (k = j + 1 ; k < max ; k++ )
      {
        idx2 = zsorted[k];

	// bonded if closer than elemental Rcov + tolerance
        cutoff = SQUARE(rad[j] + rad[k] + 0.45f); 

        zd  = SQUARE(c[idx1*3+2] - c[idx2*3+2]);
	if (zd > 25.0f ) break; // bigger than max cutoff

        d2  = SQUARE(c[idx1*3]   - c[idx2*3]);
        d2 += SQUARE(c[idx1*3+1] - c[idx2*3+1]);
        d2 += zd;
  
        if (d2 > cutoff) continue;
	if (d2 < 0.40f) continue;
  
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
  float maxlength;
  vector<OBEdgeBase*>::iterator l;
  for (atom = BeginAtom(i);atom;atom = NextAtom(i))
    {
      while (atom->BOSum() > etab.GetMaxBonds(atom->GetAtomicNum()))
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

// From current connectivity (but little or no atom types)
// assign multiple bonds to fill out valences, etc.
void OBMol::PerceiveBondOrders()
{
  OBAtom *atom, *b, *c;
  Vector v1, v2;
  int angles;
  float degrees;
  vector<OBNodeBase*>::iterator i;
  vector<OBEdgeBase*>::iterator j,k;
  
  // Pass 1: Assign estimated hybridization based on avg. angles
  for (atom = BeginAtom(i);atom;atom = NextAtom(i))
    {
      degrees = 0.0f;
      angles = 0;
      for (b = atom->BeginNbrAtom(j); b; b = atom->NextNbrAtom(j))
	{
	  k = j;
	  for (c = atom->NextNbrAtom(k); c; c = atom->NextNbrAtom(k))
	    {
	      v1 = b->GetVector() - atom->GetVector();
	      v2 = c->GetVector() - atom->GetVector();	
	      degrees += VectorAngle(v1, v2);
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
  // (GRH): This isn't necessary because flat atoms will have ~120 angles
  //        on average (exactly flat atoms = 360 deg/3 bonds => 120.0)

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

  // Sulphone -SO2-
  OBSmartsPattern sulphone; sulphone.Init("[#16D4]([#8D1])([#8D1])(*)(*)");
  if (sulphone.Match(*this))
    {
      cerr << " got match! " << endl;
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

  // Pass 6: Assign remaining bond types, ordered by atom electronegativity
  vector<pair<OBAtom*,float> > sortedAtoms;
  vector<float> rad; 
  vector<int> sorted;
  int iter, max;
  float maxElNeg, shortestBond, currentElNeg;

  for (atom = BeginAtom(i) ; atom ; atom = NextAtom(i))
    {
      pair<OBAtom*,float> entry(atom, etab.GetElectroNeg(atom->GetAtomicNum()));
      sortedAtoms.push_back(entry);
    }
  sort(sortedAtoms.begin(), sortedAtoms.end(), SortAtomZ);

  max = sortedAtoms.size();

  for (iter = 0 ; iter < max ; iter++ )
  {
      atom = sortedAtoms[iter].first;
      if ( (atom->GetHyb() == 1 || atom->GetValence() == 1)
	   && atom->BOSum() < etab.GetMaxBonds(atom->GetAtomicNum()) )
	{
	  // loop through the neighbors looking for a hybrid or terminal atom
	  // (and pick the one with highest electronegativity first)
	  // *or* pick a neighbor that's a terminal atom
	  if (atom->HasNonSingleBond())
	    continue;

	  maxElNeg = 0.0f;
	  shortestBond = 5000.0f;
	  c = NULL;
	  for (b = atom->BeginNbrAtom(j); b; b = atom->NextNbrAtom(j))
	    {
	      currentElNeg = etab.GetElectroNeg(b->GetAtomicNum());
	      if ( (b->GetHyb() == 1 || b->GetValence() == 1)
		   && b->BOSum() < etab.GetMaxBonds(b->GetAtomicNum())
		   && (currentElNeg > maxElNeg ||
		       (currentElNeg == maxElNeg
			&& (atom->GetBond(b))->GetLength() < shortestBond)) )
		{
		  if (atom->HasNonSingleBond())
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
		&& atom->BOSum() < etab.GetMaxBonds(atom->GetAtomicNum()) )
	{
	  // as above
	  if (atom->HasNonSingleBond())
	    continue;
	  
	  maxElNeg = 0.0f;
	  shortestBond = 5000.0f;
	  c = NULL;
	  for (b = atom->BeginNbrAtom(j); b; b = atom->NextNbrAtom(j))
	    {
	      currentElNeg = etab.GetElectroNeg(b->GetAtomicNum());
	      if ( (b->GetHyb() == 2 || b->GetValence() == 1)
		   && b->BOSum() < etab.GetMaxBonds(b->GetAtomicNum())
		   && (currentElNeg > maxElNeg ||
                      (currentElNeg == maxElNeg
                       && (atom->GetBond(b))->GetLength() < shortestBond)) )
		{
		  if (atom->HasNonSingleBond())
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
}

void OBMol::Center()
{
  int j,size;
  float *c,x,y,z,fsize;
  
  size = NumAtoms();
  fsize = -1.0f/(float)NumAtoms();

  vector<float*>::iterator i;
  for (i = _vconf.begin();i != _vconf.end();i++)
    {
      c = *i;
      x = y = z = 0.0f;
      for (j = 0;j < size;j++) {x += c[j*3]; y += c[j*3+1]; z += c[j*3+2];}
      x *= fsize;
      y *= fsize;
      z *= fsize;

      for (j = 0;j < size;j++) {c[j*3]+=x; c[j*3+1]+=y; c[j*3+2]+=z;}
   }

}

Vector OBMol::Center(int nconf)
{
  SetConformer(nconf);

  OBAtom *atom;
  vector<OBNodeBase*>::iterator i;

  float x=0.0,y=0.0,z=0.0;
  for (atom = BeginAtom(i);atom;atom = NextAtom(i))
    {x += atom->x();y += atom->y();z += atom->z();}

  x /= (float)NumAtoms();
  y /= (float)NumAtoms();
  z /= (float)NumAtoms();
  
  Vector vtmp;
  Vector v(x,y,z);

  for (atom = BeginAtom(i);atom;atom = NextAtom(i))
    {
      vtmp = atom->GetVector() - v;
      atom->SetVector(vtmp);
    }

  return(v);
}


void OBMol::Translate(const Vector &v)
{
  for (int i = 0;i < NumConformers();i++) 
    Translate(v,i);
}

void OBMol::Translate(const Vector &v,int nconf)
{
  int i,size;
  float x,y,z;
  float *c = (nconf == OB_CURRENT_CONFORMER)? _c : GetConformer(nconf);

  x = v.x(); y = v.y(); z = v.z();
  size = NumAtoms();
  for (i = 0;i < size;i++)
    {
      c[i*3  ] += x;
      c[i*3+1] += y;
      c[i*3+2] += z;
    }
}

void OBMol::Rotate(const float u[3][3])
{
  int i,j,k;
  float m[9];
  for (k=0,i = 0;i < 3;i++)
    for (j = 0;j < 3;j++)
      m[k++] = u[i][j];

  for (i = 0;i < NumConformers();i++) 
    Rotate(m,i);
}

void OBMol::Rotate(const float m[9])
{
  for (int i = 0;i < NumConformers();i++) 
    Rotate(m,i);
}

void OBMol::Rotate(const float m[9],int nconf)
{
  int i,size;
  float x,y,z;
  float *c = (nconf == OB_CURRENT_CONFORMER)? _c : GetConformer(nconf);
  
  size = NumAtoms();
  for (i = 0;i < size;i++)
    {
      x = c[i*3  ];
      y = c[i*3+1];
      z = c[i*3+2];
      c[i*3  ] = m[0]*x + m[1]*y + m[2]*z;
      c[i*3+1] = m[3]*x + m[4]*y + m[5]*z;
      c[i*3+2] = m[6]*x + m[7]*y + m[8]*z;
    }
}


void OBMol::SetConformers(vector<float*> &v)
{
  vector<float*>::iterator i;
  for (i = _vconf.begin();i != _vconf.end();i++)
    delete [] *i;

  _vconf = v;
  _c = (_vconf.empty()) ? NULL : _vconf[0];

}

void OBMol::CopyConformer(float *c,int idx)
{
  obAssert(!_vconf.empty() && (unsigned)idx < _vconf.size());
  memcpy((char*)_vconf[idx],(char*)c,sizeof(float)*3*NumAtoms());
}

void OBMol::CopyConformer(double *c,int idx)
{
  obAssert(!_vconf.empty() && (unsigned)idx < _vconf.size());

  unsigned int i;
  for (i = 0;i < NumAtoms();i++)
    {
      _vconf[idx][i*3  ] = (float)c[i*3  ];
      _vconf[idx][i*3+1] = (float)c[i*3+1];
      _vconf[idx][i*3+2] = (float)c[i*3+2];
    }
}

void OBMol::DeleteConformer(int idx)
{
  if (idx < 0 || idx >= (signed)_vconf.size()) return;

  delete [] _vconf[idx];
  _vconf.erase((_vconf.begin()+idx));
}

bool OBMol::Compress(void)
{
  int size = 0;
  unsigned char buf[100000];

  if (!_compressed && NumAtoms() < 256)
    {
      WriteBinary(buf, size, *this);
		
      if (size)
        {
          _compressed = true;

		  OBCompressData *cd = new OBCompressData;
		  cd->SetData(buf,size);

          Clear();
		  SetData(cd);
        }
      else
          _compressed = false;
    }

  return _compressed;
}

bool OBMol::UnCompress(void)
{
  int size = 0;
  unsigned char *buf;
  unsigned char *data;

  if (_compressed)
    {
	  OBCompressData *cd = (OBCompressData*)GetData(obCompressData);
      data = cd->GetData();

      if (data != NULL)
        {
		  size = cd->GetSize();

          if (size > 0)
            {
              _compressed = false;

              buf = new unsigned char[size];

              memcpy(buf, data, size);

              Clear();
              ReadBinary(buf,*this,size);
			  delete [] buf;
              return true;
            }
        }
    }

  return false;
}

OBAtom *OBMol::BeginAtom(vector<OBNodeBase*>::iterator &i) 
{
	i = _vatom.begin();
	return((i == _vatom.end()) ? (OBAtom*)NULL : (OBAtom*)*i);
}

OBAtom *OBMol::NextAtom(vector<OBNodeBase*>::iterator &i) 
{
	i++;
	return((i == _vatom.end()) ? (OBAtom*)NULL : (OBAtom*)*i);
}

OBBond *OBMol::BeginBond(vector<OBEdgeBase*>::iterator &i) 
{
	i = _vbond.begin();
	return((i == _vbond.end()) ? (OBBond*)NULL : (OBBond*)*i);
}

OBBond *OBMol::NextBond(vector<OBEdgeBase*>::iterator &i) 
{
	i++;
	return((i == _vbond.end()) ? (OBBond*)NULL : (OBBond*)*i);
}


//////Pose Member functions of OBMol/////////
/*!
**\brief Deletes all pose information for the OBMol
*/
void OBMol::DeletePoses()
  {
    //If there are no poses don't do anything
    if (_pose.empty()) return;

    //If the atom coordinate array is pointing to the pose
    //array change it to point to the 1st conformer
    if (_c == _xyz_pose && _c != NULL) {
        if (!_vconf.empty()) _c = _vconf[0];
        else _c = NULL;
      }

    //Free the pose coordinate array
    if (_xyz_pose) {delete [] _xyz_pose; _xyz_pose = NULL;}

    _pose.clear();

    _cur_pose_idx=0;

  }


/*!
**\brief Deletes a specified pose
**\param i Number of the pose to be deleted
*/
void OBMol::DeletePose(unsigned int i)
  {
    //Check that a valid pose is being deleted.
    if (i >= NumPoses()) return;

    //If this is the last pose just call DeletePoses.
    if (NumPoses() == 1) {DeletePoses(); return;} 

    //Delete the pose
    _pose.erase(_pose.begin()+i);

    return;
  }

/*!
**\brief Adds an additional pose to the pose array.
**\param pose An OBPose object holding the pose to be added
*/
void OBMol::AddPose(OBPose& pose)
  {
    //Check that the pose references a valid conformer
    if (pose.ConformerNum() >= (unsigned int) NumConformers()) {
        ThrowError("WARNING! Pose does not reference a valid conformer");
        cerr << "WARNING! OBMol::AddPose(OBPose)  ";
        cerr << "Pose references invalid conformer";
        cerr << endl;
        return;
      }
    _pose.push_back(pose);
  }

/*!
**\brief Sets the poses of the OBMol
**\param poses A vector of poses
*/
void OBMol::SetPoses(vector<OBPose>& poses)
  {
    //Check that all the poses reference valid conformers
    unsigned int i;
    for (i=0 ; i<poses.size() ; i++) {
        if (poses[i].ConformerNum() >= (unsigned int) NumConformers()) {
            ThrowError("WARNING! Poses do not reference valid conformers");
            cerr << "WARNING! OBMol::SetPoses(vector<OBPose>&)  ";
            cerr << "Pose references invalid conformer";
            cerr << endl;
            return;
          }
      }

    //Set the poses
    _pose = poses;

    //If the atom coordinate array was looking at poses have it look at the first new pose
    if (!_pose.empty()) {
        if (_c == _xyz_pose) SetPose(0);
      }
    else if (!_vconf.empty()) SetConformer(0);
    else _c = NULL;
  }

/*!
**\brief Sets the OBMol's atomic coordiates to the specified pose
**\param i The number of the pose
*/
void OBMol::SetPose(unsigned int i)
  {
    //check that the pose is valid
    if (i >= NumPoses()) {
        ThrowError("WARNING! Invalid pose specified");
        cerr << "WARNING! OBMol::SetPose(unsigned int) ";
        cerr << "Invalid pose specified!" << endl;
        return;
      }

    //Check that the pose references a valid conformer
    if (_pose[i].ConformerNum() >= (unsigned int) NumConformers()) {
        ThrowError("WARNING! Pose references invalid conformer");
        cerr << "WARNING! OBMol::SetPose(unsigned int) ";
        cerr << "Pose references invalid conformer" << endl;
        return;
      }

    //Make sure the pose coordinate array has memory
    if (_xyz_pose==NULL) _xyz_pose = new float [3*NumAtoms()];

    //Generate coordinates for the pose
    unsigned int ii;
    float *xyz_conf = GetConformer(_pose[i].ConformerNum());
    for (ii=0 ; ii<3*NumAtoms() ; ii++) _xyz_pose[ii] = xyz_conf[ii];
    _pose[i].CoordTrans().Transform(_xyz_pose,NumAtoms());

    //Point the atom coordinate pointer to the coordinates of the pose
    _c = _xyz_pose;

    //Remember current pose
    _cur_pose_idx = i;

    return;
  }

/*!
**\brief Returns the number of the current pose
**\return The number of the current pose.  If no poses are
**present 0 is returned.
*/
unsigned int OBMol::CurrentPoseIndex()
  {
    if (_pose.empty()) return 0;
    return _cur_pose_idx;
  }

/*!
**\brief Returns the coordinates of a pose
**\param i The number of the pose desired.
**\param xyz An array (preallocated) that is filled
**with the coordinates of the pose. 
*/
void OBMol::GetPoseCoordinates(unsigned int i, float *xyz)
  {
    //check that the pose is valid
    if (i >= NumPoses()) {
        ThrowError("WARNING! Invalid pose specified");
        cerr << "WARNING! OBMol::GetPoseCoordinates(unsigned int) ";
        cerr << "Invalid pose specified!" << endl;
        return;
      }
 
    //Check that the pose references a valid conformer
    if (_pose[i].ConformerNum() >= (unsigned int) NumConformers()) {
        ThrowError("WARNING! Pose references invalid conformer");
        cerr << "WARNING! OBMol::GetPoseCoordinates(unsigned int) ";
        cerr << "Pose references invalid conformer" << endl;
        return;
      }

    //Check that xyz is not NULL
    if (xyz==NULL) return;
 
    //Generate coordinates for the pose
    unsigned int ii;
    float *xyz_conf = GetConformer(_pose[i].ConformerNum());
    for (ii=0 ; ii<3*NumAtoms() ; ii++) xyz[ii] = xyz_conf[ii];
    _pose[i].CoordTrans().Transform(xyz,NumAtoms());
 
    return;
  }

/*!
**\brief Provied access to the pose object
**\param i The number of the pose desired.
**\return A reference to the i'th OBPose.
*/
OBPose& OBMol::GetPose(unsigned int i)
  {
/*
    //check that the pose is valid
    if (i >= NumPoses()) {
        ThrowError("WARNING! Invalid pose specified");
        cerr << "WARNING! OBMol::GetPoseCoordinates(unsigned int) ";
        cerr << "Invalid pose specified!" << endl;
        OBPose pose;
        return pose;
      }
*/
    return _pose[i];
  }

/*!
**\brief Converts the poses of the OBMol into conformers.
**This destroys the original conformer list as well as the
**original pose information.
*/
void OBMol::ChangePosesToConformers()
  {
    //If there aren't any poses don't do anything
    if (_pose.empty()) return;

    //Generate the coordinates of all the poses
    unsigned int i;
    vector<float*> dconf;
    dconf.resize(NumPoses());
    for (i=0 ; i<NumPoses() ; i++) {
        dconf[i] = new float [3*NumAtoms()];
        GetPoseCoordinates(i,dconf[i]);
      }

    //Now that we have the coordinates clear the pose info for OBMol
    DeletePoses();

    //Assign the pose coordinates to the conformers
    SetConformers(dconf);

    return;
  }
//////////End pose member functions of OBMol////////////////

//
// OBResidue member functions
//

OBResidue::OBResidue() 
{
    _chainnum = 0;
    _resnum   = 0;
    _resname  = "";
}

OBResidue::~OBResidue()
{
    vector<OBAtom*>::iterator a;
    for ( a = _atoms.begin() ; a != _atoms.end() ; a++ )
        (*a)->SetResidue(NULL);
    _atoms.clear();
}

OBResidue &OBResidue::operator=(const OBResidue &src)
     //copy residue information
{
    _chainnum = src._chainnum;
    _resnum   = src._resnum;
    _resname  = src._resname;
    _atomid   = src._atomid;
    _hetatm   = src._hetatm;
    _sernum   = src._sernum;
    return(*this);
}

void OBResidue::AddAtom(OBAtom *atom)
{
    if (atom != NULL)
    {
        atom->SetResidue(this);

        _atoms.push_back(atom);
        _atomid.push_back("");
        _hetatm.push_back(false);
        _sernum.push_back(0);
    }
}

void OBResidue::InsertAtom(OBAtom *atom)
{
    if (atom != NULL)
    {
        atom->SetResidue(this);

        _atoms.push_back(atom);
        _atomid.resize(_atoms.size());
        _hetatm.resize(_atoms.size());
        _sernum.resize(_atoms.size());
    }
}

void OBResidue::RemoveAtom(OBAtom *atom)
{
    if (atom != NULL)
    {
        int idx = GetIndex(atom);
        if (idx >= 0)
        {
            atom->SetResidue(NULL);

            _atoms.erase(_atoms.begin() + idx);
            _atomid.erase(_atomid.begin() + idx);
            _hetatm.erase(_hetatm.begin() + idx);
            _sernum.erase(_sernum.begin() + idx);
        }
    }
}

void OBResidue::Clear(void)
{
  for (unsigned int i = 0 ; i < _atoms.size() ; i++ )
    _atoms[i]->SetResidue(NULL);

    _chainnum = 0;
    _resnum   = 0;
    _resname  = "";

    _atoms.clear();
    _atomid.clear();
    _hetatm.clear();
    _sernum.clear();
}

}

/**********************************************************************
Copyright (C) 1998, 1999, 2000, 2001, 2002
OpenEye Scientific Software, Inc.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include "binary.h"
#include "fileformat.h"
#include "mol.h"
#include "obutil.h"
#include "rotor.h"

#define OB_TITLE_SIZE     254
#define OB_BINARY_SETWORD 32

using namespace std;

namespace OpenBabel
{

//test byte ordering
static int SINT = 0x00000001;
static unsigned char *STPTR = (unsigned char*)&SINT;
const bool SwabInt = (STPTR[0]!=0);

void SetRotorToAngle(float *c,OBAtom **ref,float ang,vector<int> atoms);

int Swab(int i)
{
  unsigned char tmp[4],c;
  memcpy(tmp,(char*)&i,sizeof(int));
  c = tmp[0]; tmp[0] = tmp[3]; tmp[3] = c;
  c = tmp[1]; tmp[1] = tmp[2]; tmp[2] = c;
  memcpy((char*)&i,tmp,sizeof(int));
  return(i);
}

OBRotamerList::~OBRotamerList()
{
  vector<unsigned char*>::iterator i;
  for (i = _vrotamer.begin();i != _vrotamer.end();i++)
    delete [] *i;

  vector<pair<OBAtom**,vector<int> > >::iterator j;
  for (j = _vrotor.begin();j != _vrotor.end();j++)
    delete [] j->first;

  //Delete the interal base coordinate list
  unsigned int k;
  for (k=0 ; k<_c.size() ; k++) delete [] _c[k];
}

void OBRotamerList::GetReferenceArray(unsigned char *ref)
{
  int j;
  vector<pair<OBAtom**,vector<int> > >::iterator i;
  for (j=0,i = _vrotor.begin();i != _vrotor.end();i++)
    {
      ref[j++] = (unsigned char)(i->first[0])->GetIdx();
      ref[j++] = (unsigned char)(i->first[1])->GetIdx();
      ref[j++] = (unsigned char)(i->first[2])->GetIdx();
      ref[j++] = (unsigned char)(i->first[3])->GetIdx();
    }
}

void OBRotamerList::Setup(OBMol &mol,OBRotorList &rl)
{
  //clear the old stuff out if necessary
  _vres.clear();
  vector<unsigned char*>::iterator j;
  for (j = _vrotamer.begin();j != _vrotamer.end();j++) delete [] *j;
  _vrotamer.clear();

  vector<pair<OBAtom**,vector<int> > >::iterator k;
  for (k = _vrotor.begin();k != _vrotor.end();k++)
    delete [] k->first;
  _vrotor.clear();

  //create the new list
  OBRotor *rotor;
  vector<OBRotor*>::iterator i;
  vector<int> children;

  int ref[4];
  OBAtom **atomlist;
  for (rotor = rl.BeginRotor(i);rotor;rotor = rl.NextRotor(i))
    {
      atomlist = new OBAtom* [4];
      rotor->GetDihedralAtoms(ref);
      atomlist[0] = mol.GetAtom(ref[0]); atomlist[1] = mol.GetAtom(ref[1]);
      atomlist[2] = mol.GetAtom(ref[2]); atomlist[3] = mol.GetAtom(ref[3]);
      mol.FindChildren(children,ref[1],ref[2]);
      _vrotor.push_back(pair<OBAtom**,vector<int> > (atomlist,children));
      _vres.push_back(rotor->GetResolution());
    }

  vector<float>::iterator n;
  vector<vector<float> >::iterator m;
  for (m = _vres.begin();m != _vres.end();m++)
    for (n = m->begin();n != m->end();n++)
      *n *= RAD_TO_DEG;
}

void OBRotamerList::Setup(OBMol &mol,unsigned char *ref,int nrotors)
{
  //clear the old stuff out if necessary
  _vres.clear();
  vector<unsigned char*>::iterator j;
  for (j = _vrotamer.begin();j != _vrotamer.end();j++) delete [] *j;
  _vrotamer.clear();

  vector<pair<OBAtom**,vector<int> > >::iterator k;
  for (k = _vrotor.begin();k != _vrotor.end();k++)
    delete [] k->first;
  _vrotor.clear();

  //create the new list
  int i;
  vector<int> children;

  int refatoms[4];
  OBAtom **atomlist;
  for (i = 0;i < nrotors;i++)
    {
      atomlist = new OBAtom* [4];
      refatoms[0] = (int)ref[i*4  ];
      refatoms[1] = (int)ref[i*4+1];
      refatoms[2] = (int)ref[i*4+2];
      refatoms[3] = (int)ref[i*4+3];
      mol.FindChildren(children,refatoms[1],refatoms[2]);
      atomlist[0] = mol.GetAtom(refatoms[0]); 
      atomlist[1] = mol.GetAtom(refatoms[1]);
      atomlist[2] = mol.GetAtom(refatoms[2]);
      atomlist[3] = mol.GetAtom(refatoms[3]);
      _vrotor.push_back(pair<OBAtom**,vector<int> > (atomlist,children));
    }

}

void OBRotamerList::AddRotamer(float *c)
{
  int idx,size;
  float angle,res=255.0f/360.0f;
  Vector v1,v2,v3,v4;

  unsigned char *rot = new unsigned char [_vrotor.size()+1];
  rot[0] = (char) 0;

  vector<pair<OBAtom**,vector<int> > >::iterator i;
  for (size=1,i = _vrotor.begin();i != _vrotor.end();i++,size++)
    {
      idx = (i->first[0])->GetCIdx(); v1.Set(c[idx],c[idx+1],c[idx+2]);
      idx = (i->first[1])->GetCIdx(); v2.Set(c[idx],c[idx+1],c[idx+2]);
      idx = (i->first[2])->GetCIdx(); v3.Set(c[idx],c[idx+1],c[idx+2]);
      idx = (i->first[3])->GetCIdx(); v4.Set(c[idx],c[idx+1],c[idx+2]);

      angle = CalcTorsionAngle(v1,v2,v3,v4);
      while (angle < 0.0f)   angle += 360.0f;
      while (angle > 360.0f) angle -= 360.0f;
      rot[size] = (unsigned char)rint(angle*res);
    }

  _vrotamer.push_back(rot);
}

void OBRotamerList::AddRotamer(int *arr)
{
  unsigned int i;
  float angle,res=255.0f/360.0f;

  unsigned char *rot = new unsigned char [_vrotor.size()+1];
  rot[0] = (unsigned char)arr[0];

  for (i = 0;i < _vrotor.size();i++)
    {
      angle = _vres[i][arr[i+1]];
      while (angle < 0.0f)   angle += 360.0f;
      while (angle > 360.0f) angle -= 360.0f;
      rot[i+1] = (unsigned char)rint(angle*res);
    }
  _vrotamer.push_back(rot);
}

void OBRotamerList::AddRotamer(unsigned char *arr)
{
  unsigned int i;
  float angle,res=255.0f/360.0f;

  unsigned char *rot = new unsigned char [_vrotor.size()+1];
  rot[0] = (unsigned char)arr[0];

  for (i = 0;i < _vrotor.size();i++)
    {
      angle = _vres[i][(int)arr[i+1]];
      while (angle < 0.0f)   angle += 360.0f;
      while (angle > 360.0f) angle -= 360.0f;
      rot[i+1] = (unsigned char)rint(angle*res);
    }
  _vrotamer.push_back(rot);
}

void OBRotamerList::AddRotamers(unsigned char *arr,int nrotamers)
{
  unsigned int size;
  int i;

  size = (unsigned int)_vrotor.size()+1;
  for (i = 0;i < nrotamers;i++)
    {
      unsigned char *rot = new unsigned char [size];
      memcpy(rot,&arr[i*size],sizeof(char)*size);
      _vrotamer.push_back(rot);
    }
}

void OBRotamerList::ExpandConformerList(OBMol &mol,vector<float*> &clist)
{
  unsigned int j;
  float angle,invres=360.0f/255.0f;
  unsigned char *conf;
  vector<float*> tmpclist;
  vector<unsigned char*>::iterator i;

  for (i = _vrotamer.begin();i != _vrotamer.end();i++)
    {
      conf = *i;
      float *c = new float [mol.NumAtoms()*3];
      memcpy(c,clist[(int)conf[0]],sizeof(float)*mol.NumAtoms()*3);

      for (j = 0;j < _vrotor.size();j++)
	{
	  angle = invres*((float)conf[j+1]);
	  if (angle > 180.0) angle -= 360.0;
	  SetRotorToAngle(c,_vrotor[j].first,angle,_vrotor[j].second);
	}
      tmpclist.push_back(c);
    }

  //transfer the conf list
  vector<float*>::iterator k;
  for (k = clist.begin();k != clist.end();k++)
    delete [] *k;
  clist = tmpclist;
}

//Create a conformer list using the internal base set of coordinates
vector<float*> OBRotamerList::CreateConformerList(OBMol& mol)
  {
    unsigned int j;
    float angle,invres=360.0f/255.0f;
    unsigned char *conf;
    vector<float*> tmpclist;
    vector<unsigned char*>::iterator i;
   
    for (i = _vrotamer.begin();i != _vrotamer.end();i++)
      {
        conf = *i;
        float *c = new float [mol.NumAtoms()*3];
        memcpy(c,_c[(int)conf[0]],sizeof(float)*mol.NumAtoms()*3);
   
        for (j = 0;j < _vrotor.size();j++)
          {
            angle = invres*((float)conf[j+1]);
            if (angle > 180.0) angle -= 360.0;
            SetRotorToAngle(c,_vrotor[j].first,angle,_vrotor[j].second);
          }
        tmpclist.push_back(c);
      }
   
    return tmpclist;
  }

//Copies the coordinates in bc, NOT the pointers, into the object
void OBRotamerList::SetBaseCoordinateSets(vector<float*> bc, unsigned int N)
  {
    unsigned int i,j;

    //Clear out old data
    for (i=0 ; i<_c.size() ; i++) delete [] _c[i];
    _c.clear();

    //Copy new data
    float *c = NULL;
    float *cc= NULL;
    for (i=0 ; i<bc.size() ; i++) {
        c = new float [3*N];
        cc = bc[i];
        for (j=0 ; j<3*N ; j++) c[j] = cc[j];
        _c.push_back(c);
      }
    _NBaseCoords = N;
  }



int PackCoordinate(float c[3],float max[3])
{
  int tmp;
  tmp  = ((int)(c[0]*max[0])) << 20;
  tmp |= ((int)(c[1]*max[1])) << 10;
  tmp |= ((int)(c[2]*max[2]));
  return(tmp);
}

void UnpackCoordinate(float c[3],float max[3],int tmp)
{
  c[0] = (float)(tmp>>20);            c[0] *= max[0];
  c[1] = (float)((tmp&0xffc00)>>10);  c[1] *= max[1];
  c[2] = (float)(tmp&0x3ff);          c[2] *= max[2];
}

bool WriteBinary(ostream &ofs, OBMol &mol)
{
    int    size = 0;
    string buf;

    mol.SetOutputType(OEBINARY);
    WriteBinary(buf, size, mol);

    int tmp = size;
    if (SwabInt) 
        tmp = Swab(tmp);

    ofs.write((char*)&tmp, sizeof(int));
    ofs.write((char*)buf.data(), size);

    return true;
}

bool WriteBinary(string &buf, int &size, OBMol &mol)
{
    vector<float*>::iterator j;
    unsigned int k, sz;
    int idx, m, tmp;

    idx = 0;

    //read title first
    const char *title = mol.GetTitle();
    
    int len = (int) strlen(title);
    if (len > OB_TITLE_SIZE) 
        len = OB_TITLE_SIZE;
    
    if (len > 0)
    {
        buf += (char) len;
        buf.insert(++idx, title, len);
        idx += len;
    }
    else
    {
        buf += (char) 0;
        idx++;
    }

    tmp = (mol.NumAtoms() << 16) | mol.NumBonds();
    if (SwabInt) 
        tmp = Swab(tmp);
    
    buf.insert(idx, (const char *) &tmp, sizeof(int));
    idx += (int) sizeof(int);

    OBAtom *atom;
    vector<OBNodeBase*>::iterator atomIter;
    for (atom = mol.BeginAtom(atomIter); atom; atom = mol.NextAtom(atomIter))
      {
        buf += (unsigned char) atom->GetAtomicNum();
        idx++;
      }

    OBBond *bond;
    vector<OBEdgeBase*>::iterator bondIter;
    for (bond = mol.BeginBond(bondIter); bond; bond = mol.NextBond(bondIter))
      {
        buf += (unsigned char) bond->GetBeginAtomIdx();
        buf += (unsigned char) bond->GetEndAtomIdx();
        buf += (unsigned char) bond->GetBO();
        idx += 3;
      }

    //Write out conformers and coordinates
    OBRotamerList *rml = static_cast<OBRotamerList*>(mol.GetData(obRotamerList));

    //find min and max
    int   imin[3], imax[3];
    float  min[3] = { 10E10f, 10E10f, 10E10f};
    float  max[3] = {-10E10f,-10E10f,-10E10f};

    vector<float*> clist;

    //If we have a rotamer list with internally stored base coordinates
    //have clist point to them rather than the molecules conformer coordinates

    if (rml && ((rml) ? rml->NumRotamers() : 0) && rml->NumBaseCoordinateSets()) 
    {
        if (rml->NumAtoms() == mol.NumAtoms()) //Error check, these should match
        {
            for ( k = 0 ; k < rml->NumBaseCoordinateSets() ; k++ ) 
                clist.push_back(rml->GetBaseCoordinateSet(k));
        }
        else 
            clist = mol.GetConformers();
    }
    else 
        clist = mol.GetConformers();

    sz = mol.NumAtoms() * 3;
    for ( j = clist.begin() ; j != clist.end() ; j++ )
    {
        float *coords = *j;
        for ( k = 0 ; k < sz ; k += 3 )
        {
            if (coords[k  ] < min[0]) min[0] = coords[k  ];
            if (coords[k+1] < min[1]) min[1] = coords[k+1];
            if (coords[k+2] < min[2]) min[2] = coords[k+2];
            if (coords[k  ] > max[0]) max[0] = coords[k  ];
            if (coords[k+1] > max[1]) max[1] = coords[k+1];
            if (coords[k+2] > max[2]) max[2] = coords[k+2];
        }
    }

    //store integer versions of min and max
    for ( k = 0 ; k < 3 ; k++ )
    {
        max[k] -= min[k];
        imin[k] = (int) (1000000.0f * min[k]);
        imax[k] = (int) (1000000.0f * max[k]);

        if (SwabInt)
        {
            imin[k] = Swab(imin[k]);
            imax[k] = Swab(imax[k]);
        }
    }

    //write the min and max
    buf.insert(idx, (const char *) imin, sizeof(int) * 3);
    idx += (int) sizeof(int) * 3;

    buf.insert(idx, (const char *) imax, sizeof(int) * 3);
    idx += (int) sizeof(int) * 3;

    //quantize max for packing coordinates
    for ( k = 0 ; k < 3 ; k++ ) 
        max[k] = (fabs(max[k]) > 0.01) ? 1023.0f / max[k] : 0.0;

    //write the number of confs and rotamers
    tmp = (int) clist.size();
    if (SwabInt) 
        tmp = Swab(tmp);
    
    buf.insert(idx, (const char *) &tmp, sizeof(int));
    idx += (int) sizeof(int);

    tmp = (rml) ? rml->NumRotamers() : 0;
    if (SwabInt) 
        tmp = Swab(tmp);
    
    buf.insert(idx, (const char *) &tmp, sizeof(int));
    idx += (int) sizeof(int);

    //The third boolean expression in the next if statement is an error check to make
    //sure than the number of atoms in the OBRotamerLists internal coordinates are the
    //same as the molecules IF we are using the OBRotamerLists internal coordinates.
    //If we are using the OBRotamerList's internal coordinates but the number of atoms
    //in those coordinates don't match the number of atoms in the molecule then  the
    //rotamer list is incorrect for the molecule and we just default to writing the
    //molecules conformers.  I put this in because it strikes me as a very easy error
    //to make if the molecule is modified in any way and the user is not aware that
    //he/she is responsible for correctly updating the OBRotamerList MM 4/20/01

    if (rml && ((rml) ? rml->NumRotamers() : 0) && 
        (rml->NumBaseCoordinateSets() == 0 || rml->NumAtoms() == mol.NumAtoms()) ) 
        // Store conformers as torsion list
    {        
        // Write base coordinates

        float tc[3];
        sz = mol.NumAtoms() * 3;
        for ( j = clist.begin() ; j != clist.end() ; j++ )
        {
            float *coords = *j;
            for ( k = 0 ; k < sz ; k += 3 )
            {
                tc[0] = coords[k  ] - min[0];
                tc[1] = coords[k+1] - min[1];
                tc[2] = coords[k+2] - min[2];

                tmp = PackCoordinate(tc, max);
                if (SwabInt) 
                    tmp = Swab(tmp);

                buf.insert(idx, (const char *) &tmp, sizeof(int));
                idx += (int) sizeof(int);
            }
        }

        //Write rotors
        tmp = rml->NumRotors(); 
        if (SwabInt) 
            tmp = Swab(tmp);

        buf.insert(idx, (const char *) &tmp, sizeof(int));
        idx += (int) sizeof(int);

        unsigned char *ref = new unsigned char [rml->NumRotors() * 4];

        rml->GetReferenceArray(ref);

        buf.insert(idx, (const char *) ref, rml->NumRotors() * 4);
        idx += (int) rml->NumRotors() * 4;

        delete [] ref;

        //Write Rotamers
        vector<unsigned char*>::iterator r;
        sz = rml->NumRotors() + 1;
        for ( r = rml->BeginRotamer() ; r != rml->EndRotamer() ; r++ )
        {
            buf.insert(idx, (const char *) *r, sz);
            idx += sz;
        }
    }
    else if (mol.NumConformers() > 1) // Store conformers as coordinates
    {
        // Write coordinate conformers 

        float tc[3];
        sz = mol.NumAtoms() * 3;
        for ( j = clist.begin() ; j != clist.end() ; j++ )
        {
            float *coords = *j;
            for ( k = 0 ; k < sz ; k += 3 )
            {
                tc[0] = coords[k  ] - min[0];
                tc[1] = coords[k+1] - min[1];
                tc[2] = coords[k+2] - min[2];
                
                tmp = PackCoordinate(tc, max);
                if (SwabInt) 
                    tmp = Swab(tmp);

                buf.insert(idx, (const char *) &tmp, sizeof(int));
                idx += (int) sizeof(int);
            }
        }
    }
    else //must be storing single-conformer structure
    {
        //write the coordinates
        float coord[3];
	OBAtom *atom2;
	vector<OBNodeBase*>::iterator i;
	for (atom2 = mol.BeginAtom(i); atom2; atom2 = mol.NextAtom(i))
        {
            (atom2->GetVector()).Get(coord);
            coord[0] -= min[0];
            coord[1] -= min[1];
            coord[2] -= min[2];

            tmp = PackCoordinate(coord,max);
            if (SwabInt) 
                tmp = Swab(tmp);
            
            buf.insert(idx, (const char *) &tmp, sizeof(int));
            idx += (int) sizeof(int);
        }
    }

    if (mol.NumAtoms()) //set bits on for aromatic atoms
    {
        int nwords = mol.NumAtoms() / OB_BINARY_SETWORD;
        if (mol.NumAtoms() % OB_BINARY_SETWORD) 
            nwords++;
    
        unsigned int *arobits = new unsigned int [nwords];
        memset((char*)arobits, '\0', sizeof(int)*nwords);

        int word, bit;
	OBAtom *atom3;
	for (atom3 = mol.BeginAtom(atomIter); atom3; atom3 = mol.NextAtom(atomIter))
        {
            if (atom3->IsAromatic())
            {
                word           = (atom3->GetIdx()-1) / OB_BINARY_SETWORD;
                bit            = (atom3->GetIdx()-1) % OB_BINARY_SETWORD;
                arobits[word] |= (1<<bit);
            }
        }

        if (SwabInt)
            for ( m = 0; m < nwords ; m++)
                arobits[m] = Swab(arobits[m]);

        buf.insert(idx, (const char *) arobits, sizeof(int) * nwords);
        idx += (int) sizeof(int) * nwords;

        delete [] arobits;
    }

    if (mol.NumBonds()) //set bits on for aromatic bonds
    {
        int nwords = mol.NumBonds() / OB_BINARY_SETWORD;
        if (mol.NumBonds() % OB_BINARY_SETWORD) 
            nwords++;

        unsigned int *arobits = new unsigned int [nwords];
        memset((char*)arobits, '\0', sizeof(int)*nwords);

        int word, bit;
	OBBond *bond2;
	for (bond2 = mol.BeginBond(bondIter); bond2; bond2 = mol.NextBond(bondIter))
        {
            if (bond2->IsAromatic())
            {
                word           = (bond2->GetIdx()) / OB_BINARY_SETWORD;
                bit            = (bond2->GetIdx()) % OB_BINARY_SETWORD;
                arobits[word] |= (1<<bit);
            }
        }

        if (SwabInt)
            for ( m = 0 ; m < nwords ; m++)
                arobits[m] = Swab(arobits[m]);

        buf.insert(idx, (const char *) arobits, sizeof(int) * nwords);
        idx += (int) sizeof(int) * nwords;

        delete [] arobits;
    }

    size = idx;
    return true;
}


bool ReadBinary(istream &ifs, OBMol &mol)
{
    int size = 0;
    if (!ifs.read((char*)&size,sizeof(int))) 
        return false;

    if (SwabInt) 
        size = Swab(size);

    if (size > 0)
    {
        unsigned char *buf = new unsigned char[size];

        if (!ifs.read((char*)buf,(int)sizeof(char)*size)) 
            return(false);

        ReadBinary(buf,mol,size);

        delete [] buf;
        return true;
    }
    else
        return false;
}


bool ReadBinary(istream &ifs, unsigned char **bin)
{
    obAssert(bin != NULL);

    int temp = 0;
    if (!ifs.read((char*)&temp,sizeof(int)))
        return false;

    int size = (SwabInt) ? Swab(temp) : temp;

    if (size > 0)
    {
        unsigned char *buf = new unsigned char[sizeof(int) + size];

        memcpy(buf, &temp, sizeof(int));
        if (!ifs.read((char*)&buf[sizeof(int)], (int)sizeof(char)*size))
            return false;

        *bin = buf;        
        return true;
    }
    else
    {
        *bin = NULL;
        return false;
    }
}


bool ReadBinary(unsigned char *buf, OBMol &mol, int size)
{
  int i,j,k,idx,natoms,nbonds,tmp;
  char title[OB_TITLE_SIZE+1]; 
  idx = 0;

  //read title
  i = (int)buf[0];
  idx += (int)sizeof(char);
  if (i > 0)
    {
      memcpy(title,&buf[idx],sizeof(char)*i);
      title[i] = '\0';
      idx += i;
    }
  else
    strcpy(title,"****");

  //readnumber of atoms and bonds
  memcpy(&tmp,(unsigned char*)&buf[idx],sizeof(int)); idx += (int)sizeof(int);
  if (SwabInt) tmp = Swab(tmp);

  natoms = (tmp >> 16);
  nbonds = tmp & 0xffff;

  unsigned char *anum = new unsigned char [natoms];
  memcpy((unsigned char*)anum,
	 &buf[idx],sizeof(unsigned char)*natoms); 
  idx += (int)sizeof(unsigned char)*natoms;
  
  mol.BeginModify();
  //read atom data
  OBAtom atom;
  for (i = 0;i < natoms;i++)
    {
      atom.SetAtomicNum((int)anum[i]);
      atom.SetType(etab.GetSymbol((int)anum[i]));
      mol.AddAtom(atom);
      atom.Clear();
    }
  delete [] anum;

  //read bond data
  int start,end,order;
  unsigned char *bnd = new unsigned char [nbonds*3];
  memcpy((unsigned char*)bnd,
	 &buf[idx],sizeof(unsigned char)*3*nbonds); 
  idx += (int)sizeof(unsigned char)*3*nbonds;
  for (i = 0;i < nbonds;i++)
    {
      start = bnd[i*3  ];
      end   = bnd[i*3+1];
      order = bnd[i*3+2];
      mol.AddBond(start,end,order);
    }
  delete [] bnd;

  mol.EndModify();

  //read the min and max
  int imin[3],imax[3];
  float min[3],max[3];
  memcpy((char*)imin,&buf[idx],sizeof(int)*3); idx += (int)sizeof(int)*3;
  memcpy((char*)imax,&buf[idx],sizeof(int)*3); idx += (int)sizeof(int)*3;

  //unpack min and max
  for (i = 0;i < 3;i++)
    {
      if (SwabInt)
	{
	  imin[i] = Swab(imin[i]);
	  imax[i] = Swab(imax[i]);
	}
      min[i] = (float) imin[i]/1000000.0f;
      max[i] = (float) imax[i]/1000000.0f;
      max[i] = (fabs(max[i])> 0.01) ? max[i]/1023.0f:0.0; 
    }

  //read conformer information if available
  int nconfs,rotmrs;
  memcpy((char*)&nconfs,&buf[idx],sizeof(int)); idx += (int)sizeof(int);
  memcpy((char*)&rotmrs,&buf[idx],sizeof(int)); idx += (int)sizeof(int);
  if (SwabInt)
    {
      nconfs = Swab(nconfs);
      rotmrs = Swab(rotmrs);
    }

  if (nconfs == 1 && !rotmrs) //only a single conformer
    {
      Vector v;
      int *tmpi = new int [natoms];
      memcpy((char*)tmpi,&buf[idx],sizeof(int)*natoms); idx += (int)sizeof(int)*natoms;
      float coord[3];
      for (i = 0;i < natoms;i++)
	{
	  if (SwabInt) tmpi[i] = Swab(tmpi[i]);
	  UnpackCoordinate(coord,max,tmpi[i]);
	  for (j = 0;j < 3;j++) coord[j] += min[j];
	  v.Set(coord);
	  (mol.GetAtom(i+1))->SetVector(v);
	}
      delete [] tmpi;
    }
  else
    {
      int *tmpi = new int [natoms];
      vector<float*> cltmp;
      for (i = 0;i < nconfs;i++)
	{
	  memcpy((char*)tmpi,(unsigned char*)&buf[idx],sizeof(int)*natoms); 
	  idx += (int)sizeof(int)*natoms;
	  float *coord = new float [mol.NumAtoms()*3];
	  for (j = 0;j < natoms;j++) 
	    {
	      if (SwabInt) tmpi[j] = Swab(tmpi[j]);
	      UnpackCoordinate(&coord[j*3],max,tmpi[j]);
	      for (k = 0;k < 3;k++) coord[j*3+k] += min[k];
	    }
	  cltmp.push_back(coord);
	}

      delete [] tmpi;

      if (!rotmrs)
	mol.SetConformers(cltmp);
      else
	{
       
	  int nrotors;
	  OBRotamerList *rml = new OBRotamerList;
	  memcpy((char*)&nrotors,&buf[idx],sizeof(int)); idx += (int)sizeof(int);
	  if (SwabInt) nrotors = Swab(nrotors);

	  unsigned char *ref = new unsigned char [nrotors*4];
	  memcpy((unsigned char*)ref,&buf[idx],sizeof(unsigned char)*nrotors*4);
	  idx += (int)sizeof(unsigned char)*nrotors*4;
	  rml->Setup(mol,ref,nrotors);
	  delete [] ref;
	  
	  unsigned char *rotamers = new unsigned char [(nrotors+1)*rotmrs];
	  memcpy((unsigned char*)rotamers,
		 &buf[idx],sizeof(unsigned char)*(nrotors+1)*rotmrs);
	  idx += (int)sizeof(unsigned char)*(nrotors+1)*rotmrs;
	  rml->AddRotamers(rotamers,rotmrs);
	  delete [] rotamers;
	 
          //Copy the base coordinate list into the OBRotamerList object
          rml->SetBaseCoordinateSets(cltmp,mol.NumAtoms());
 
	  //expand rotamer information to a conformer list
	  rml->ExpandConformerList(mol,cltmp);
	  mol.SetConformers(cltmp);

          //Add the OBRotamerList to the molecule as user data
          mol.SetData(rml);
	}
    }

  mol.SetTitle(title);

  if (idx >= size) return(true);
  
  int nwords;
  unsigned int *arobits;

  if (mol.NumAtoms()) //set bits on for aromatic atoms
    {
      nwords = mol.NumAtoms()/OB_BINARY_SETWORD;
      if (mol.NumAtoms()%OB_BINARY_SETWORD) nwords++;
      arobits = new unsigned int [nwords];

      memcpy((unsigned char*)arobits,&buf[idx],sizeof(int)*nwords);
      idx += (int)sizeof(int)*nwords;

      if (SwabInt)
	for (i = 0;i < nwords;i++) arobits[i] =  Swab(arobits[i]);

      for (i = 0;i < (signed)mol.NumAtoms();i++)
	if ((arobits[i/OB_BINARY_SETWORD]>>(i%OB_BINARY_SETWORD))&1)
	  mol.GetAtom(i+1)->SetAromatic();

      delete [] arobits;
    }

  if (mol.NumBonds()) //set bits on for aromatic atoms
    {
      nwords = (mol.NumBonds()/OB_BINARY_SETWORD);
      if (mol.NumBonds()%OB_BINARY_SETWORD) nwords++;
      arobits = new unsigned int [nwords];

      memcpy((unsigned char*)arobits,&buf[idx],sizeof(int)*nwords);
      idx += (int)sizeof(int)*nwords;

      if (SwabInt)
	for (i = 0;i < nwords;i++) arobits[i] =  Swab(arobits[i]);

      for (i = 0;i < (signed)mol.NumBonds();i++)
	if ((arobits[i/OB_BINARY_SETWORD]>>(i%OB_BINARY_SETWORD))&1)
	  mol.GetBond(i)->SetAromatic();

      delete [] arobits;
    }

  mol.SetAromaticPerceived();

  return(true);
}

void SetRotorToAngle(float *c,OBAtom **ref,float ang,vector<int> atoms)
     //this function will rotate the coordinates of 'atoms'
     //such that tor == ang - atoms in 'tor' should be ordered such 
     //that the 3rd atom is the pivot around which atoms rotate
     //ang is in degrees
{
  float v1x,v1y,v1z,v2x,v2y,v2z,v3x,v3y,v3z;
  float c1x,c1y,c1z,c2x,c2y,c2z,c3x,c3y,c3z;
  float c1mag,c2mag,radang,costheta,m[9];
  float x,y,z,mag,rotang,sn,cs,t,tx,ty,tz;

  int tor[4];
  tor[0] = ref[0]->GetCIdx();
  tor[1] = ref[1]->GetCIdx();
  tor[2] = ref[2]->GetCIdx();
  tor[3] = ref[3]->GetCIdx();

  //
  //calculate the torsion angle
  //
  v1x = c[tor[0]]   - c[tor[1]];   v2x = c[tor[1]]   - c[tor[2]];
  v1y = c[tor[0]+1] - c[tor[1]+1]; v2y = c[tor[1]+1] - c[tor[2]+1];
  v1z = c[tor[0]+2] - c[tor[1]+2]; v2z = c[tor[1]+2] - c[tor[2]+2];
  v3x = c[tor[2]]   - c[tor[3]];
  v3y = c[tor[2]+1] - c[tor[3]+1];
  v3z = c[tor[2]+2] - c[tor[3]+2];

  c1x = v1y*v2z - v1z*v2y;   c2x = v2y*v3z - v2z*v3y;
  c1y = -v1x*v2z + v1z*v2x;  c2y = -v2x*v3z + v2z*v3x;
  c1z = v1x*v2y - v1y*v2x;   c2z = v2x*v3y - v2y*v3x;
  c3x = c1y*c2z - c1z*c2y;
  c3y = -c1x*c2z + c1z*c2x;
  c3z = c1x*c2y - c1y*c2x; 
  
  c1mag = c1x*c1x + c1y*c1y + c1z*c1z;
  c2mag = c2x*c2x + c2y*c2y + c2z*c2z;
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
  rotang = (DEG_TO_RAD*ang) - radang; 

  sn = sin(rotang); cs = cos(rotang);t = 1 - cs;
  //normalize the rotation vector
  mag = sqrt(v2x*v2x + v2y*v2y + v2z*v2z);
  x = v2x/mag; y = v2y/mag; z = v2z/mag;
  
  //set up the rotation matrix
  m[0]= t*x*x + cs;     m[1] = t*x*y + sn*z;  m[2] = t*x*z - sn*y;
  m[3] = t*x*y - sn*z;  m[4] = t*y*y + cs;    m[5] = t*y*z + sn*x;
  m[6] = t*x*z + sn*y;  m[7] = t*y*z - sn*x;  m[8] = t*z*z + cs;

  //
  //now the matrix is set - time to rotate the atoms
  //
  tx = c[tor[1]];ty = c[tor[1]+1];tz = c[tor[1]+2];
  vector<int>::iterator i;int j;
  for (i = atoms.begin();i != atoms.end();i++)
    {
      j = ((*i)-1)*3;
      c[j] -= tx;c[j+1] -= ty;c[j+2]-= tz;
      x = c[j]*m[0] + c[j+1]*m[1] + c[j+2]*m[2];
      y = c[j]*m[3] + c[j+1]*m[4] + c[j+2]*m[5];
      z = c[j]*m[6] + c[j+1]*m[7] + c[j+2]*m[8];
      c[j] = x; c[j+1] = y; c[j+2] = z;
      c[j] += tx;c[j+1] += ty;c[j+2] += tz;
    }
}

//OBBinaryDBase class - facilitates random access to OBBinary files

OBBinaryDBase::OBBinaryDBase(const char *fname)
{
    int size;
    std::streampos pos;

    _ifs.open(fname);
    if (!_ifs)
    {
        exit(0);
    }

    for (;;)
    {
        pos = _ifs.tellg();
        
        if (!_ifs.read((char*)&size,sizeof(int))) 
            break;
        
        if (SwabInt) 
            size = Swab(size);
        
        if (!_ifs.seekg(size, ostream::cur))
            break;
        
        _vpos.push_back(pos);
    }

    _ifs.close();
    _ifs.open(fname);
    if (!_ifs) 
    {
        exit(0);
    }
}

OBBinaryDBase::OBBinaryDBase(string &fname)
{
    int size;
    std::streampos pos;

    _ifs.open(fname.c_str());
    if (!_ifs)
    {
        exit(0);
    }

    for (;;)
    {
        pos = _ifs.tellg();
        
        if (!_ifs.read((char*)&size,sizeof(int))) 
            break;
        
        if (SwabInt) 
            size = Swab(size);
        
        if (!_ifs.seekg(size, ostream::cur))
            break;
        
        _vpos.push_back(pos);
    }

    _ifs.close();
    _ifs.open(fname.c_str());
    if (!_ifs)
    {
        exit(0);
    }
}

int OBBinaryDBase::Size()
{
    return (int) _vpos.size();
}

void OBBinaryDBase::GetMolecule(OBMol &mol,int idx)
{
    OBFileFormat translator;
    mol.Clear();
    mol.SetInputType(OEBINARY);
    _ifs.seekg(_vpos[idx]);
    translator.ReadMolecule(_ifs,mol);
}

OBBinaryDBase::~OBBinaryDBase(void)
{
}


} //namespace OpenEye


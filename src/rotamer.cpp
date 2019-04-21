/**********************************************************************
rotamer.cpp - Handle rotamer list data.

Copyright (C) 1998, 1999, 2000-2002 OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2006 by Geoffrey R. Hutchison

This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/
#include <openbabel/babelconfig.h>

#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/ring.h>
#include <openbabel/obiter.h>
#include <openbabel/rotamer.h>

#define OB_TITLE_SIZE     254
#define OB_BINARY_SETWORD 32

using namespace std;

namespace OpenBabel
{

  /** \class OBRotamerList rotamer.h <openbabel/rotamer.h>

      A high-level class for rotamer / conformer generation, intended mainly
      for use with the related class OBRotorList and the OBRotorRules database

      Rotamers represent conformational isomers formed simply by rotation of
      dihedral angles. On the other hand, conformers may include geometric
      relaxation (i.e., slight modification of bond lengths, bond angles, etc.)

      The following shows an example of generating 2 conformers using different
      rotor states. Similar code could be used for systematic or Monte Carlo
      conformer sampling when combined with energy evaluation (molecular
      mechanics or otherwise).

      \code
      OBRotorList rl; // used to sample all rotatable bonds via the OBRotorRules data
      // If you want to "fix" any particular atoms (i.e., freeze them in space)
      // then set up an OBBitVec of the fixed atoms and call
      // rl.SetFixAtoms(bitvec);
      rl.Setup(mol);

      // How many rotatable bonds are there?
      cerr << " Number of rotors: " << rl.Size() << endl;

      // indexed from 1, rotorKey[0] = 0
      std::vector<int> rotorKey(rl.Size() + 1, 0);

      // each entry represents the configuration of a rotor
      // e.g. indexes into OBRotor::GetResolution() -- the different angles
      //   to sample for a rotamer search
      for (unsigned int i = 0; i < rl.Size() + 1; ++i)
      rotorKey[i] = 0; // could be anything from 0 .. OBRotor->GetResolution().size()
      // -1 is for no rotation

      // The OBRotamerList can generate conformations (i.e., coordinate sets)
      OBRotamerList rotamers;
      rotamers.SetBaseCoordinateSets(mol);
      rotamers.Setup(mol, rl);

      rotamers.AddRotamer(rotorKey);
      rotorKey[1] = 2; // switch one rotor
      rotamers.AddRotamer(rotorKey);

      rotamers.ExpandConformerList(mol, mol.GetConformers());

      // change the molecule conformation
      mol.SetConformer(0); // rotorKey 0, 0, ...
      conv.Write(&mol);

      mol.SetConformer(1); // rotorKey 0, 2, ...

      \endcode

  **/

  //test byte ordering
  static int SINT = 0x00000001;
  static unsigned char *STPTR = (unsigned char*)&SINT;
  const bool SwabInt = (STPTR[0]!=0);

#if !HAVE_RINT
  inline double rint(double x)
  {
    return ( (x < 0.0) ? ceil(x-0.5) : floor(x+0.5));
  }
#endif

  void SetRotorToAngle(double *c,OBAtom **ref,double ang,vector<int> atoms);

  int Swab(int i)
  {
    unsigned char tmp[4],c;
    memcpy(tmp,(char*)&i,sizeof(int));
    c = tmp[0];
    tmp[0] = tmp[3];
    tmp[3] = c;
    c = tmp[1];
    tmp[1] = tmp[2];
    tmp[2] = c;
    memcpy((char*)&i,tmp,sizeof(int));
    return(i);
  }

  OBGenericData* OBRotamerList::Clone(OBBase* newparent) const
  {
    //Since the class contains OBAtom pointers, the new copy use
    //these from the new molecule, newparent
    OBMol* newmol = static_cast<OBMol*>(newparent);

    OBRotamerList *new_rml = new OBRotamerList;
    new_rml->_attr = _attr;
    new_rml->_type = _type;

    //Set base coordinates
    unsigned int k,l;
    vector<double*> bc;
    double *c=NULL;
    double *cc=NULL;
    for (k=0 ; k<NumBaseCoordinateSets() ; ++k)
      {
        c = new double [3*NumAtoms()];
        cc = GetBaseCoordinateSet(k);
        for (l=0 ; l<3*NumAtoms() ; ++l)
					c[l] = cc[l];
        bc.push_back(c);
      }
    if (NumBaseCoordinateSets())
			new_rml->SetBaseCoordinateSets(bc,NumAtoms());

    //Set reference array
    unsigned char *ref = new unsigned char [NumRotors()*4];
    if (ref)
      {
        GetReferenceArray(ref);
        new_rml->Setup(*newmol,ref,NumRotors());
        delete [] ref;
      }

    //Set Rotamers
    unsigned char *rotamers = new unsigned char [(NumRotors()+1)*NumRotamers()];
    if (rotamers)
      {
        vector<unsigned char*>::const_iterator kk;
        unsigned int idx=0;
        for (kk = _vrotamer.begin();kk != _vrotamer.end();++kk)
          {
            memcpy(&rotamers[idx],(const unsigned char*)*kk,sizeof(unsigned char)*(NumRotors()+1));
            idx += sizeof(unsigned char)*(NumRotors()+1);
          }
        new_rml->AddRotamers(rotamers,NumRotamers());
        delete [] rotamers;
      }
    return new_rml;
  }

  OBRotamerList::~OBRotamerList()
  {
    vector<unsigned char*>::iterator i;
    for (i = _vrotamer.begin();i != _vrotamer.end();++i)
      delete [] *i;

    vector<pair<OBAtom**,vector<int> > >::iterator j;
    for (j = _vrotor.begin();j != _vrotor.end();++j)
      delete [] j->first;

    //Delete the interal base coordinate list
    unsigned int k;
    for (k=0 ; k<_c.size() ; ++k)
      delete [] _c[k];
  }

  void OBRotamerList::GetReferenceArray(unsigned char *ref)const
  {
    int j;
		vector<pair<OBAtom**,vector<int> > >::const_iterator i;
    for (j=0,i = _vrotor.begin();i != _vrotor.end();++i)
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
    for (j = _vrotamer.begin();j != _vrotamer.end();++j)
      delete [] *j;
    _vrotamer.clear();

    vector<pair<OBAtom**,vector<int> > >::iterator k;
    for (k = _vrotor.begin();k != _vrotor.end();++k)
      delete [] k->first;
    _vrotor.clear();
    _vrings.clear();
    _vringTors.clear();

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
        atomlist[0] = mol.GetAtom(ref[0]);
        atomlist[1] = mol.GetAtom(ref[1]);
        atomlist[2] = mol.GetAtom(ref[2]);
        atomlist[3] = mol.GetAtom(ref[3]);
        mol.FindChildren(children,ref[1],ref[2]);
        _vrotor.push_back(pair<OBAtom**,vector<int> > (atomlist,children));
        _vres.push_back(rotor->GetResolution());
      }

    // if the rotor list has ring bonds, build up an index
    if (rl.HasRingRotors()){
      // go through rings
      // for each step of the path, see if there's a matching rotor
      vector<int> path;
      int pSize;
      vector<double> ringTorsions;
      vector<int> ringRotors;
      FOR_RINGS_OF_MOL(r, mol)
      {
        ringTorsions.clear();
        ringRotors.clear();

        pSize = r->Size();
        if (pSize < 4)
          continue; // not rotatable

        path = r->_path;
        for (int j = 0; j < pSize; ++j) {
          double torsion = mol.GetTorsion(path[(j + pSize - 1) % pSize],
                                       path[(j + pSize) % pSize],
                                       path[(j + pSize + 1) % pSize],
                                       path[(j + pSize + 2) % pSize]);
          ringTorsions.push_back(torsion);

          // now check to see if any of these things are rotors
          int rotorIndex = -1; // not a rotor
          OBBond *bond = mol.GetBond(path[(j + pSize) % pSize], path[(j + pSize + 1) % pSize]);
          for (rotor = rl.BeginRotor(i);rotor;rotor = rl.NextRotor(i))
            {
              if (bond != rotor->GetBond())
                continue; // no match at all

              // Central bond matches, make sure 1..4 atoms are in the path
              rotor->GetDihedralAtoms(ref);
              if ( (ref[0] == path[(j + pSize - 1) % pSize] &&
                    ref[3] == path[(j + pSize + 2) % pSize])
                   ||
                   (ref[3] == path[(j + pSize - 1) % pSize] &&
                    ref[0] == path[(j + pSize + 2) % pSize]) ) {
                rotorIndex = rotor->GetIdx();
              }
            } // end checking all the rotors
          ringRotors.push_back(rotorIndex); // could be -1 if it's not rotatable
        }

        _vringTors.push_back(ringTorsions);
        _vrings.push_back(ringRotors);
      } // finished with the rings
    } // if (ring rotors)

    vector<double>::iterator n;
    vector<vector<double> >::iterator m;
    for (m = _vres.begin();m != _vres.end();++m)
      for (n = m->begin();n != m->end();++n)
        *n *= RAD_TO_DEG;
  }

  void OBRotamerList::Setup(OBMol &mol,unsigned char *ref,int nrotors)
  {
    //clear the old stuff out if necessary
    _vres.clear();
    vector<unsigned char*>::iterator j;
    for (j = _vrotamer.begin();j != _vrotamer.end();++j)
      delete [] *j;
    _vrotamer.clear();

    vector<pair<OBAtom**,vector<int> > >::iterator k;
    for (k = _vrotor.begin();k != _vrotor.end();++k)
      delete [] k->first;
    _vrotor.clear();
    _vrings.clear();
    _vringTors.clear();

    //create the new list
    int i;
    vector<int> children;

    int refatoms[4];
    OBAtom **atomlist;
    for (i = 0; i < nrotors; ++i)
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

  void OBRotamerList::AddRotamer(double *c)
  {
    int idx,size;
    double angle,res=255.0/360.0;
    vector3 v1,v2,v3,v4;

    unsigned char *rot = new unsigned char [_vrotor.size()+1];
    rot[0] = (char) 0;

    vector<pair<OBAtom**,vector<int> > >::iterator i;
    for (size=1,i = _vrotor.begin();i != _vrotor.end();++i,++size)
      {
        idx = (i->first[0])->GetCoordinateIdx();
        v1.Set(c[idx],c[idx+1],c[idx+2]);
        idx = (i->first[1])->GetCoordinateIdx();
        v2.Set(c[idx],c[idx+1],c[idx+2]);
        idx = (i->first[2])->GetCoordinateIdx();
        v3.Set(c[idx],c[idx+1],c[idx+2]);
        idx = (i->first[3])->GetCoordinateIdx();
        v4.Set(c[idx],c[idx+1],c[idx+2]);

        angle = CalcTorsionAngle(v1,v2,v3,v4);
        while (angle < 0.0)
          angle += 360.0;
        while (angle > 360.0)
          angle -= 360.0;
        rot[size] = (unsigned char)rint(angle*res);
      }

    _vrotamer.push_back(rot);
  }

  void OBRotamerList::AddRotamer(int *arr)
  {
    unsigned int i;
    double angle,res=255.0/360.0;

    unsigned char *rot = new unsigned char [_vrotor.size()+1];
    rot[0] = (unsigned char)arr[0];

    for (i = 0;i < _vrotor.size();++i)
      {
        angle = _vres[i][arr[i+1]];
        while (angle < 0.0)
          angle += 360.0;
        while (angle > 360.0)
          angle -= 360.0;
        rot[i+1] = (unsigned char)rint(angle*res);
      }
    _vrotamer.push_back(rot);
  }

  void OBRotamerList::AddRotamer(std::vector<int> arr)
  {
    unsigned int i;
    double angle,res=255.0/360.0;

    if (arr.size() != (_vrotor.size() + 1))
      return; // wrong size key

    // gotta check for weird ring torsion combinations
    if (_vrings.size()) {
      // go through each ring and update the possible torsions
      for (unsigned int j = 0; j < _vrings.size(); ++j) {
        vector<int> path = _vrings[j];
        double torsionSum = 0.0;

        // go around the loop and add up the torsions
        for (unsigned int i = 0; i < path.size(); ++i) {
          if (path[i] == -1) {
            // not a rotor, use the fixed value
            torsionSum += _vringTors[j][i];
            continue;
          }

          // what angle are we trying to use with this key?
          angle = _vres[ path[i] ][arr[ path[i]+1 ]]*res;
          while (angle < 0.0)
          	angle += 360.0;
        	while (angle > 360.0)
          	angle -= 360.0;

          // update the ring torsion for this setting
          _vringTors[j][i] = angle;
          torsionSum += angle;
        }

        // if the sum of the ring torsions is not ~0, bad move
        if (fabs(torsionSum) > 45.0) {
          //          cerr << " Bad move! " << fabs(torsionSum) << endl;
          return; // don't make the move
        }
        //        cerr << " Good move!" << endl;
      }
    }

    unsigned char *rot = new unsigned char [_vrotor.size()+1];
    rot[0] = (unsigned char)arr[0];

    for (i = 0;i < _vrotor.size();++i)
      {
        angle = _vres[i][arr[i+1]];
        while (angle < 0.0)
          angle += 360.0;
        while (angle > 360.0)
          angle -= 360.0;
        rot[i+1] = (unsigned char)rint(angle*res);
      }
    _vrotamer.push_back(rot);
  }

  void OBRotamerList::AddRotamer(unsigned char *arr)
  {
    unsigned int i;
    double angle,res=255.0/360.0;

    unsigned char *rot = new unsigned char [_vrotor.size()+1];
    rot[0] = (unsigned char)arr[0];

    for (i = 0;i < _vrotor.size();++i)
      {
        angle = _vres[i][(int)arr[i+1]];
        while (angle < 0.0)
          angle += 360.0;
        while (angle > 360.0)
          angle -= 360.0;
        rot[i+1] = (unsigned char)rint(angle*res);
      }
    _vrotamer.push_back(rot);
  }

  void OBRotamerList::AddRotamers(unsigned char *arr,int nrotamers)
  {
    unsigned int size;
    int i;

    size = (unsigned int)_vrotor.size()+1;
    for (i = 0;i < nrotamers;++i)
      {
        unsigned char *rot = new unsigned char [size];
        memcpy(rot,&arr[i*size],sizeof(char)*size);
        _vrotamer.push_back(rot);
      }
  }

  void OBRotamerList::ExpandConformerList(OBMol &mol,vector<double*> &clist)
  {
    vector<double*> tmpclist = CreateConformerList(mol);

    //transfer the conf list
    vector<double*>::iterator k;
    for (k = clist.begin();k != clist.end();++k)
      delete [] *k;
    clist = tmpclist;
  }

  //! Create a conformer list using the internal base set of coordinates
  vector<double*> OBRotamerList::CreateConformerList(OBMol& mol)
  {
    unsigned int j;
    double angle,invres=360.0/255.0;
    unsigned char *conf;
    vector<double*> tmpclist;
    vector<unsigned char*>::iterator i;

    for (i = _vrotamer.begin();i != _vrotamer.end();++i)
      {
        conf = *i;
        double *c = new double [mol.NumAtoms()*3];
        memcpy(c,_c[(int)conf[0]],sizeof(double)*mol.NumAtoms()*3);

        for (j = 0;j < _vrotor.size();++j)
          {
            angle = invres*((double)conf[j+1]);
            if (angle > 180.0)
              angle -= 360.0;
            SetRotorToAngle(c,_vrotor[j].first,angle,_vrotor[j].second);
          }
        tmpclist.push_back(c);
      }

    return tmpclist;
  }

  //! Change the current coordinate set
  //! \since version 2.2
  bool OBRotamerList::SetCurrentCoordinates(OBMol &mol, std::vector<int> arr)
  {
    double angle;

    if (arr.size() != (_vrotor.size() + 1))
      return false; // wrong size key

    // gotta check for weird ring torsion combinations
    if (_vrings.size()) {
      // go through each ring and update the possible torsions
      for (unsigned int j = 0; j < _vrings.size(); ++j) {
        vector<int> path = _vrings[j];
        double torsionSum = 0.0;

        // go around the loop and add up the torsions
        for (unsigned int i = 0; i < path.size(); ++i) {
          if (path[i] == -1) {
            // not a rotor, use the fixed value
            torsionSum += _vringTors[j][i];
            continue;
          }

          // what angle are we trying to use with this key?
          angle = _vres[ path[i] ][arr[ path[i]+1 ]];
          while (angle < 0.0)
          	angle += 360.0;
        	while (angle > 360.0)
          	angle -= 360.0;

          // update the ring torsion for this setting
          _vringTors[j][i] = angle;
          torsionSum += angle;
        }

        // if the sum of the ring torsions is not ~0, bad move
        if (fabs(torsionSum) > 45.0) {
          //          cerr << " Bad move!" << endl;
          return false; // don't make the move
        }
      }
    }

    double *c = mol.GetCoordinates();
    for (unsigned int i = 0;i < _vrotor.size();++i)
      {
				if (arr[i+1] == -1) // skip this rotor
					continue;
				else {
        	angle = _vres[i][arr[i+1]];
        	while (angle < 0.0)
          	angle += 360.0;
        	while (angle > 360.0)
          	angle -= 360.0;
	        SetRotorToAngle(c,_vrotor[i].first,angle,_vrotor[i].second);
				} // set an angle
      } // for rotors

    return true;
  }

  void OBRotamerList::SetBaseCoordinateSets(OBMol& mol)
  {
    SetBaseCoordinateSets(mol.GetConformers(), mol.NumAtoms());
  }

  //Copies the coordinates in bc, NOT the pointers, into the object
  void OBRotamerList::SetBaseCoordinateSets(vector<double*> bc, unsigned int N)
  {
    unsigned int i,j;

    //Clear out old data
    for (i=0 ; i<_c.size() ; ++i)
      delete [] _c[i];
    _c.clear();

    //Copy new data
    double *c = NULL;
    double *cc= NULL;
    for (i=0 ; i<bc.size() ; ++i)
      {
        c = new double [3*N];
        cc = bc[i];
        for (j=0 ; j<3*N ; ++j)
          c[j] = cc[j];
        _c.push_back(c);
      }
    _NBaseCoords = N;
  }

  //! Rotate the coordinates of 'atoms'
  //! such that tor == ang.
  //! Atoms in 'tor' should be ordered such that the 3rd atom is
  //! the pivot around which atoms rotate (ang is in degrees)
  //! \todo This code is identical to OBMol::SetTorsion() and should be combined
  void SetRotorToAngle(double *c, OBAtom **ref,double ang,vector<int> atoms)
  {
    double v1x,v1y,v1z,v2x,v2y,v2z,v3x,v3y,v3z;
    double c1x,c1y,c1z,c2x,c2y,c2z,c3x,c3y,c3z;
    double c1mag,c2mag,radang,costheta,m[9];
    double x,y,z,mag,rotang,sn,cs,t,tx,ty,tz;

    int tor[4];
    tor[0] = ref[0]->GetCoordinateIdx();
    tor[1] = ref[1]->GetCoordinateIdx();
    tor[2] = ref[2]->GetCoordinateIdx();
    tor[3] = ref[3]->GetCoordinateIdx();

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

    if (costheta < -0.999999) costheta = -0.999999;
    if (costheta >  0.999999) costheta =  0.999999;

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
    if (mag < 0.1) mag = 0.1; // avoid divide by zero error
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
    for (i = atoms.begin();i != atoms.end();++i)
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

  int PackCoordinate(double c[3],double max[3])
  {
    int tmp;
    double cf;
    cf = c[0];
    tmp  = ((int)(cf*max[0])) << 20;
    cf = c[1];
    tmp |= ((int)(cf*max[1])) << 10;
    cf = c[2];
    tmp |= ((int)(cf*max[2]));
    return(tmp);
  }

  void UnpackCoordinate(double c[3],double max[3],int tmp)
  {
    double cf;
    cf = (double)(tmp>>20);
    c[0] = cf;
    c[0] *= max[0];
    cf = (double)((tmp&0xffc00)>>10);
    c[1] = cf;
    c[1] *= max[1];
    cf = (double)(tmp&0x3ff);
    c[2] = cf;
    c[2] *= max[2];
  }

} //namespace OpenBabel

//! \file rotamer.cpp
//! \brief Handle rotamer list data.

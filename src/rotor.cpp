/**********************************************************************
rotor.cpp - Rotate dihedral angles according to rotor rules.

Copyright (C) 1998-2000 by OpenEye Scientific Software, Inc.
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
#include <openbabel/bond.h>
#include <openbabel/ring.h>
#include <openbabel/obiter.h>
#include <openbabel/oberror.h>
#include <openbabel/rotor.h>
#include <openbabel/graphsym.h>
#include <openbabel/elements.h>

#include <set>
#include <assert.h>

// private data headers with default parameters
#include "torlib.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


using namespace std;
namespace OpenBabel
{

  //! Default step resolution for a dihedral angle (in degrees)
#define OB_DEFAULT_DELTA 15.0
  static bool GetDFFVector(OBMol&,vector<int>&,OBBitVec&);
  static bool CompareRotor(const pair<OBBond*,int>&,const pair<OBBond*,int>&);


  //**************************************
  //**** OBRotorList Member Functions ****
  //**************************************

  bool OBRotorList::Setup(OBMol &mol, bool sampleRingBonds)
  {
    Clear();
    // find the rotatable bonds
    FindRotors(mol, sampleRingBonds);
    if (!Size())
      return(false);

    // set the atoms that should be evaluated when this rotor changes
    SetEvalAtoms(mol);
    AssignTorVals(mol);

    OBRotor *rotor;
    vector<OBRotor*>::iterator i;
    for (rotor = BeginRotor(i);rotor;rotor = NextRotor(i))
      if (!rotor->Size())
        {
          int ref[4];
          rotor->GetDihedralAtoms(ref);
          char buffer[BUFF_SIZE];
          snprintf(buffer, BUFF_SIZE, "The rotor has no associated torsion values -> %d %d %d %d",
		  ref[0],ref[1],ref[2],ref[3]);
          obErrorLog.ThrowError(__FUNCTION__, buffer, obDebug);
        }

    // Reduce the number of torsions to be checked through symmetry considerations
    if (_removesym)
      RemoveSymVals(mol);

    return(true);
  }

  bool OBRotorList::FindRotors(OBMol &mol, bool sampleRingBonds)
  {
    // Find ring atoms & bonds
    // This function will set OBBond::IsRotor().
    mol.FindRingAtomsAndBonds();

    obErrorLog.ThrowError(__FUNCTION__,
                          "Ran OpenBabel::FindRotors", obAuditMsg);

    //
    // Score the bonds using the graph theoretical distance (GTD).
    // The GTD is the distance from atom i to every other atom j.
    // Atoms on the "inside" of the molecule will have a lower GTD
    // value than atoms on the "outside"
    //
    // The scoring will rank "inside" bonds first.
    //
    vector<int> gtd;
    mol.GetGTDVector(gtd);
    // compute the scores
    vector<OBBond*>::iterator i;
    vector<pair<OBBond*,int> > vtmp;
    for (OBBond *bond = mol.BeginBond(i);bond;bond = mol.NextBond(i)) {
      // check if the bond is "rotatable"
      if (bond->IsRotor(sampleRingBonds)) {
        // check if the bond is fixed (using deprecated fixed atoms or new fixed bonds)
        if ((HasFixedAtoms() || HasFixedBonds()) && IsFixedBond(bond))
          continue;

        if (bond->IsInRing()) {
          //otherwise mark that we have them and add it to the pile
          _ringRotors = true;
        }

        int score = gtd[bond->GetBeginAtomIdx()-1] + gtd[bond->GetEndAtomIdx()-1];
        // compute the GTD bond score as sum of atom GTD scores
        vtmp.push_back(pair<OBBond*,int> (bond,score));
      }
    }

    // sort the rotatable bonds by GTD score
    sort(vtmp.begin(),vtmp.end(),CompareRotor);

    // create rotors for the bonds
    int count = 0;
    vector<pair<OBBond*,int> >::iterator j;
    for (j = vtmp.begin(); j != vtmp.end(); ++j, ++count) {
      OBRotor *rotor = new OBRotor;
      rotor->SetBond((*j).first);
      rotor->SetIdx(count);
      rotor->SetNumCoords(mol.NumAtoms()*3);
      _rotor.push_back(rotor);
    }

    return true;
  }

  bool OBRotorList::IsFixedBond(OBBond *bond)
  {
    if (_fixedatoms.IsEmpty() && _fixedbonds.IsEmpty())
      return false;

    // new fixed bonds
    if (!_fixedbonds.IsEmpty()) {
      return _fixedbonds.BitIsSet(bond->GetIdx());
    }

    if (_fixedatoms.IsEmpty())
      return false;

    // deprecated fixed atoms
    OBAtom *a1,*a2,*a3;
    vector<OBBond*>::iterator i;

    a1 = bond->GetBeginAtom();
    a2 = bond->GetEndAtom();
    if (!_fixedatoms[a1->GetIdx()] || !_fixedatoms[a2->GetIdx()])
      return(false);

    bool isfixed=false;
    for (a3 = a1->BeginNbrAtom(i);a3;a3 = a1->NextNbrAtom(i))
      if (a3 != a2 && _fixedatoms[a3->GetIdx()])
        {
          isfixed = true;
          break;
        }

    if (!isfixed)
      return(false);

    isfixed = false;
    for (a3 = a2->BeginNbrAtom(i);a3;a3 = a2->NextNbrAtom(i))
      if (a3 != a1 && _fixedatoms[a3->GetIdx()])
        {
          isfixed = true;
          break;
        }

    return(isfixed);
  }

  bool GetDFFVector(OBMol &mol,vector<int> &dffv,OBBitVec &bv)
  {
    dffv.clear();
    dffv.resize(mol.NumAtoms());

    int dffcount,natom;
    OBBitVec used,curr,next;
    OBAtom *atom,*atom1;
    OBBond *bond;
    vector<OBAtom*>::iterator i;
    vector<OBBond*>::iterator j;

    next.Clear();

    for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
      {
        if (bv[atom->GetIdx()])
          {
            dffv[atom->GetIdx()-1] = 0;
            continue;
          }

        dffcount = 0;
        used.Clear();
        curr.Clear();
        used.SetBitOn(atom->GetIdx());
        curr.SetBitOn(atom->GetIdx());

        while (!curr.IsEmpty() && (bv&curr).IsEmpty())
          {
            next.Clear();
            for (natom = curr.NextBit(-1);natom != curr.EndBit();natom = curr.NextBit(natom))
              {
                atom1 = mol.GetAtom(natom);
                for (bond = atom1->BeginBond(j);bond;bond = atom1->NextBond(j))
                  if (!used.BitIsSet(bond->GetNbrAtomIdx(atom1)) &&
                      !curr.BitIsSet(bond->GetNbrAtomIdx(atom1)))
                      if (bond->GetNbrAtom(atom1)->GetAtomicNum() != OBElements::Hydrogen)
                      next.SetBitOn(bond->GetNbrAtomIdx(atom1));
              }

            used |= next;
            curr = next;
            dffcount++;
          }

        dffv[atom->GetIdx()-1] = dffcount;
      }

    return(true);
  }

  void OBRotorList::RemoveSymVals(OBMol &mol)
  {
    OBGraphSym gs(&mol);
    vector<unsigned int> sym_classes;
    gs.GetSymmetry(sym_classes);

    OBRotor *rotor;
    vector<OBRotor*>::iterator i;
    std::set<unsigned int> syms;
    for (rotor = BeginRotor(i);rotor;rotor = NextRotor(i)) {
      OBBond* bond = rotor->GetBond();
      OBAtom* end = bond->GetEndAtom();
      OBAtom* begin = bond->GetBeginAtom();
      int N_fold_symmetry = 1;
      for (int here=0; here <= 1; ++here) { // Try each side of the bond in turn

        OBAtom *this_side, *other_side;
        if (here == 0) {
          this_side = begin; other_side = end;
        }
        else {
          this_side = end; other_side = begin;
        }

        for (unsigned int hyb=2; hyb<=3; ++hyb) { // sp2 and sp3 carbons, with explicit Hs
          if (this_side->GetAtomicNum() == 6 && this_side->GetHyb() == hyb && this_side->GetExplicitDegree() == (hyb + 1) ) {
            syms.clear();
            FOR_NBORS_OF_ATOM(nbr, this_side) {
              if ( &(*nbr) == other_side ) continue;
              syms.insert(sym_classes[nbr->GetIdx() - 1]);
            }
            if (syms.size() == 1) // All of the rotated atoms have the same symmetry class
              N_fold_symmetry *= hyb;
                      }
                  }
              }

      if (N_fold_symmetry  > 1) {
        size_t old_size = rotor->Size();
        rotor->RemoveSymTorsionValues(N_fold_symmetry);
        if (!_quiet) {
          cout << "...." << N_fold_symmetry << "-fold symmetry at rotor between " <<
                 begin->GetIdx() << " and " << end->GetIdx();
          cout << " - reduced from " << old_size << " to " << rotor->Size() << endl;
                  }
              }
      }
  }

  bool OBRotorList::SetEvalAtoms(OBMol &mol)
  {
    int j;
    OBBond *bond;
    OBAtom *a1,*a2;
    OBRotor *rotor;
    vector<OBRotor*>::iterator i;
    OBBitVec eval,curr,next;
    vector<OBBond*>::iterator k;

    for (rotor = BeginRotor(i);rotor;rotor = NextRotor(i))
      {
        bond = rotor->GetBond();
        curr.Clear();
        eval.Clear();
        curr.SetBitOn(bond->GetBeginAtomIdx());
        curr.SetBitOn(bond->GetEndAtomIdx());
        eval |= curr;

        //follow all non-rotor bonds and add atoms to eval list
        for (;!curr.IsEmpty();)
          {
            next.Clear();
            for (j = curr.NextBit(0);j != curr.EndBit();j = curr.NextBit(j))
              {
                a1 = mol.GetAtom(j);
                for (a2 = a1->BeginNbrAtom(k);a2;a2 = a1->NextNbrAtom(k))
                  if (!eval[a2->GetIdx()])
                    if (!((OBBond*)*k)->IsRotor(_ringRotors)||((HasFixedAtoms()||HasFixedBonds())&&IsFixedBond((OBBond*)*k)))
                      {
                        next.SetBitOn(a2->GetIdx());
                        eval.SetBitOn(a2->GetIdx());
                      }
              }
            curr = next;
          }

        //add atoms alpha to eval list
        next.Clear();
        for (j = eval.NextBit(0);j != eval.EndBit();j = eval.NextBit(j))
          {
            a1 = mol.GetAtom(j);
            for (a2 = a1->BeginNbrAtom(k);a2;a2 = a1->NextNbrAtom(k))
              next.SetBitOn(a2->GetIdx());
          }
        eval |= next;
        rotor->SetEvalAtoms(eval);
      }

    return(true);
  }

  bool OBRotorList::AssignTorVals(OBMol &mol)
  {
    vector<OBRotor*>::iterator i;
    for (i = _rotor.begin(); i != _rotor.end(); ++i) {
      OBRotor *rotor = *i;
      OBBond *bond = rotor->GetBond();

      // query the rotor database
      int ref[4];
      vector<double> angles;
      double delta;
      _rr.GetRotorIncrements(mol, bond, ref, angles, delta);
      rotor->SetTorsionValues(angles);
      rotor->SetDelta(delta);

      // Find the smallest set of atoms to rotate. There are two candidate sets,
      // one on either side of the bond. If the first tried set size plus one is
      // larger than half of the number of atoms, the other set is smaller and the
      // rotor atom indexes are inverted (.i.e. a-b-c-d -> d-c-b-a).
      vector<int> atoms;
      vector<int>::iterator j;
      // find all atoms for which there is a path to ref[2] without goinf through ref[1]
      mol.FindChildren(atoms, ref[1], ref[2]);
      if (atoms.size() + 1 > mol.NumAtoms() / 2) {
        atoms.clear();
        // select the other smaller set
        mol.FindChildren(atoms, ref[2], ref[1]);
        // invert the rotor
        swap(ref[0],ref[3]);
        swap(ref[1],ref[2]);
      }

      // translate the rotate atom indexes to coordinate indexes (i.e. from 0, multiplied by 3)
      for (j = atoms.begin(); j != atoms.end(); ++j)
        *j = ((*j) - 1) * 3;
      // set the rotate atoms and dihedral atom indexes
      rotor->SetRotAtoms(atoms);
      rotor->SetDihedralAtoms(ref);
    }

    return true;
  }

  bool OBRotorList::SetRotAtoms(OBMol &mol)
  {
    OBRotor *rotor;
    vector<int> rotatoms;
    vector<OBRotor*>::iterator i;

    int ref[4];
    vector<int>::iterator j;
    for (rotor = BeginRotor(i);rotor;rotor = NextRotor(i))
      {
        rotor->GetDihedralAtoms(ref);

        mol.FindChildren(rotatoms,ref[1],ref[2]);
        if (rotatoms.size()+1 > mol.NumAtoms()/2)
          {
            rotatoms.clear();
            mol.FindChildren(rotatoms,ref[2],ref[1]);
            swap(ref[0],ref[3]);
            swap(ref[1],ref[2]);
          }

        for (j = rotatoms.begin();j != rotatoms.end();++j)
          *j = ((*j)-1)*3;
        rotor->SetRotAtoms(rotatoms);
        rotor->SetDihedralAtoms(ref);
      }

    return(true);
  }

  void OBRotorList::SetRotAtomsByFix(OBMol &mol)
  {
    int ref[4];
    OBRotor *rotor;
    vector<int> rotatoms;
    vector<int>::iterator j;
    vector<OBRotor*>::iterator i;

    GetDFFVector(mol,_dffv,_fixedatoms);

    for (rotor = BeginRotor(i);rotor;rotor = NextRotor(i))
      {
        rotatoms.clear();
        rotor->GetDihedralAtoms(ref);

        if (_fixedatoms[ref[1]] && _fixedatoms[ref[2]])
          {
            if (!_fixedatoms[ref[0]])
              {
                swap(ref[0],ref[3]);
                swap(ref[1],ref[2]);
                mol.FindChildren(rotatoms,ref[1],ref[2]);
                for (j = rotatoms.begin();j != rotatoms.end();++j)
                  *j = ((*j)-1)*3;
                rotor->SetRotAtoms(rotatoms);
                rotor->SetDihedralAtoms(ref);
              }
          }
        else
          if (_dffv[ref[1]-1] > _dffv[ref[2]-1])
            {
              swap(ref[0],ref[3]);
              swap(ref[1],ref[2]);
              mol.FindChildren(rotatoms,ref[1],ref[2]);
              for (j = rotatoms.begin();j != rotatoms.end();++j)
                *j = ((*j)-1)*3;
              rotor->SetRotAtoms(rotatoms);
              rotor->SetDihedralAtoms(ref);
            }
      }
  }

  OBRotorList::OBRotorList()
  {
    _rotor.clear();
    _quiet = true;
    _removesym = true;
    _ringRotors = false;
  }

  OBRotorList::~OBRotorList()
  {
    vector<OBRotor*>::iterator i;
    for (i = _rotor.begin();i != _rotor.end();++i)
      delete *i;
  }

  void OBRotorList::Clear()
  {
    vector<OBRotor*>::iterator i;
    for (i = _rotor.begin();i != _rotor.end();++i)
      delete *i;
    _rotor.clear();
    _ringRotors = false;
    //_fix.Clear();
  }

  bool CompareRotor(const pair<OBBond*,int> &a,const pair<OBBond*,int> &b)
  {
    //return(a.second > b.second); //outside->in
    return(a.second < b.second);   //inside->out
  }

  //**********************************
  //**** OBRotor Member Functions ****
  //**********************************

  OBRotor::OBRotor()
  {
  }

  void OBRotor::SetRings()
  {
    _rings.clear();
    if (_bond == NULL)
      return; // nothing to do

    vector<OBRing*> rlist;
    vector<OBRing*>::iterator i;

    OBMol *mol = _bond->GetParent();

    if (mol == NULL)
      return; // nothing to do

    rlist = mol->GetSSSR();
    for (i = rlist.begin();i != rlist.end();++i) {
      if ((*i)->IsMember(_bond))
        _rings.push_back(*i);
    }
  }

  double OBRotor::CalcTorsion(double *c)
  {
    double v1x,v1y,v1z,v2x,v2y,v2z,v3x,v3y,v3z;
    double c1x,c1y,c1z,c2x,c2y,c2z,c3x,c3y,c3z;
    double c1mag,c2mag,ang,costheta;

    //
    //calculate the torsion angle
    //
    v1x = c[_torsion[0]] -  c[_torsion[1]];
    v1y = c[_torsion[0]+1] - c[_torsion[1]+1];
    v1z = c[_torsion[0]+2] - c[_torsion[1]+2];
    v2x = c[_torsion[1]] - c[_torsion[2]];
    v2y = c[_torsion[1]+1] - c[_torsion[2]+1];
    v2z = c[_torsion[1]+2] - c[_torsion[2]+2];
    v3x = c[_torsion[2]]   - c[_torsion[3]];
    v3y = c[_torsion[2]+1] - c[_torsion[3]+1];
    v3z = c[_torsion[2]+2] - c[_torsion[3]+2];

    c1x = v1y*v2z - v1z*v2y;
    c2x = v2y*v3z - v2z*v3y;
    c1y = -v1x*v2z + v1z*v2x;
    c2y = -v2x*v3z + v2z*v3x;
    c1z = v1x*v2y - v1y*v2x;
    c2z = v2x*v3y - v2y*v3x;
    c3x = c1y*c2z - c1z*c2y;
    c3y = -c1x*c2z + c1z*c2x;
    c3z = c1x*c2y - c1y*c2x;

    c1mag = SQUARE(c1x)+SQUARE(c1y)+SQUARE(c1z);
    c2mag = SQUARE(c2x)+SQUARE(c2y)+SQUARE(c2z);
    if (c1mag*c2mag < 0.01)
      costheta = 1.0; //avoid div by zero error
    else
      costheta = (c1x*c2x + c1y*c2y + c1z*c2z)/(sqrt(c1mag*c2mag));

    if (costheta < -0.9999999)
      costheta = -0.9999999;
    if (costheta >  0.9999999)
      costheta =  0.9999999;

    if ((v2x*c3x + v2y*c3y + v2z*c3z) > 0.0)
      ang = -acos(costheta);
    else
      ang = acos(costheta);

    return(ang);
  }

  double OBRotor::CalcBondLength(double *c)
  {
    // compute the difference
    double dx, dy, dz;
    dx = c[_torsion[1]  ] - c[_torsion[2]  ];
    dy = c[_torsion[1]+1] - c[_torsion[2]+1];
    dz = c[_torsion[1]+2] - c[_torsion[2]+2];
    // compute the length
    return sqrt(SQUARE(dx) + SQUARE(dy) + SQUARE(dz));
  }


  void OBRotor::SetRotor(double *c,int idx,int prev)
  {
    double ang,sn,cs,t,dx,dy,dz,mag;

    if (prev == -1)
      ang = _torsionAngles[idx] - CalcTorsion(c);
    else
      ang = _torsionAngles[idx] - _torsionAngles[prev];

    sn = sin(ang);
    cs = cos(ang);
    t = 1 - cs;

    dx = c[_torsion[1]]   - c[_torsion[2]];
    dy = c[_torsion[1]+1] - c[_torsion[2]+1];
    dz = c[_torsion[1]+2] - c[_torsion[2]+2];
    mag = sqrt(SQUARE(dx) + SQUARE(dy) + SQUARE(dz));

    Set(c,sn,cs,t,1.0/mag);
  }

  // used in combination with Set(double *coordinates, int idx)
  void OBRotor::Precompute(double *coordinates)
  {
    _imag = 1.0 / CalcBondLength(coordinates);
    _refang = CalcTorsion(coordinates);
  }

  // used in combination with Precompute(double *coordinates)
  void OBRotor::Set(double *coordinates, int idx)
  {
    double ang,sn,cs,t;

    // compute the rotation angle
    ang = _torsionAngles[idx] - _refang;
    // compute some values for the rotation matrix
    sn = sin(ang);
    cs = cos(ang);
    t = 1 - cs;

    double x,y,z,tx,ty,tz,m[9];

    // compute the bond vector
    x = coordinates[_torsion[1]]   - coordinates[_torsion[2]];
    y = coordinates[_torsion[1]+1] - coordinates[_torsion[2]+1];
    z = coordinates[_torsion[1]+2] - coordinates[_torsion[2]+2];

    // normalize the bond vector
    x *= _imag;
    y *= _imag;
    z *= _imag;

    // set up the rotation matrix
    tx = t*x;
    ty = t*y;
    tz = t*z;
    m[0]= tx*x + cs;
    m[1] = tx*y + sn*z;
    m[2] = tx*z - sn*y;
    m[3] = tx*y - sn*z;
    m[4] = ty*y + cs;
    m[5] = ty*z + sn*x;
    m[6] = tx*z + sn*y;
    m[7] = ty*z - sn*x;
    m[8] = tz*z + cs;

    //
    //now the matrix is set - time to rotate the atoms
    //
    tx = coordinates[_torsion[1]];
    ty = coordinates[_torsion[1]+1];
    tz = coordinates[_torsion[1]+2];
    unsigned int i, j;
    for (i = 0; i < _rotatoms.size(); ++i)
      {
        j = _rotatoms[i];
        coordinates[j] -= tx;
        coordinates[j+1] -= ty;
        coordinates[j+2]-= tz;
        x = coordinates[j]*m[0] + coordinates[j+1]*m[1] + coordinates[j+2]*m[2];
        y = coordinates[j]*m[3] + coordinates[j+1]*m[4] + coordinates[j+2]*m[5];
        z = coordinates[j]*m[6] + coordinates[j+1]*m[7] + coordinates[j+2]*m[8];
        coordinates[j] = x+tx;
        coordinates[j+1] = y+ty;
        coordinates[j+2] = z+tz;
      }
  }

  // used in combination with Set
  void OBRotor::Precalc(vector<double*> &cv)
  {
    double *c,ang;
    vector<double*>::iterator i;
    vector<double>::iterator j;
    vector<double> cs,sn,t;
    for (i = cv.begin();i != cv.end();++i)
      {
        c = *i;
        cs.clear();
        sn.clear();
        t.clear();
        ang = CalcTorsion(c);

        for (j = _torsionAngles.begin();j != _torsionAngles.end();++j)
          {
            cs.push_back(cos(*j-ang));
            sn.push_back(sin(*j-ang));
            t.push_back(1 - cos(*j-ang));
          }

        _cs.push_back(cs);
        _sn.push_back(sn);
        _t.push_back(t);
        _invmag.push_back(1.0/CalcBondLength(c));
      }
  }


  void OBRotor::Set(double *c,double sn,double cs,double t,double invmag)
  {
    double x,y,z,tx,ty,tz,m[9];

    x = c[_torsion[1]]   - c[_torsion[2]];
    y = c[_torsion[1]+1] - c[_torsion[2]+1];
    z = c[_torsion[1]+2] - c[_torsion[2]+2];

    //normalize the rotation vector

    x *= invmag;
    y *= invmag;
    z *= invmag;

    //set up the rotation matrix
    tx = t*x;
    ty = t*y;
    tz = t*z;
    m[0]= tx*x + cs;
    m[1] = tx*y + sn*z;
    m[2] = tx*z - sn*y;
    m[3] = tx*y - sn*z;
    m[4] = ty*y + cs;
    m[5] = ty*z + sn*x;
    m[6] = tx*z + sn*y;
    m[7] = ty*z - sn*x;
    m[8] = tz*z + cs;

    //
    //now the matrix is set - time to rotate the atoms
    //
    tx = c[_torsion[1]];
    ty = c[_torsion[1]+1];
    tz = c[_torsion[1]+2];
    unsigned int i, j;
    for (i = 0; i < _rotatoms.size(); ++i)
      {
        j = _rotatoms[i];
        c[j] -= tx;
        c[j+1] -= ty;
        c[j+2]-= tz;
        x = c[j]*m[0] + c[j+1]*m[1] + c[j+2]*m[2];
        y = c[j]*m[3] + c[j+1]*m[4] + c[j+2]*m[5];
        z = c[j]*m[6] + c[j+1]*m[7] + c[j+2]*m[8];
        c[j] = x+tx;
        c[j+1] = y+ty;
        c[j+2] = z+tz;
      }
  }

  //! Remove all torsions angles between 0 and 360/fold
  void OBRotor::RemoveSymTorsionValues(int fold)
  {
    vector<double>::iterator i;
    vector<double> tv;
    if (_torsionAngles.size() == 1)
      return;

    for (i = _torsionAngles.begin();i != _torsionAngles.end();++i)
      if (*i >= 0.0 && *i < 2.0*M_PI / fold)
            tv.push_back(*i);

    if (tv.empty())
      return;
    _torsionAngles = tv;
  }

  void OBRotor::SetDihedralAtoms(std::vector<int> &ref)
  {
    if (ref.size() != 4)
      return;
    // copy indexes starting from 1
    _ref.resize(4);
    for (int i = 0;i < 4;++i)
      _ref[i] = ref[i];
    _torsion.resize(4);
    // convert the indexes (start from 0, multiplied by 3) for easy access to coordinates
    _torsion[0] = (ref[0]-1)*3;
    _torsion[1] = (ref[1]-1)*3;
    _torsion[2] = (ref[2]-1)*3;
    _torsion[3] = (ref[3]-1)*3;
  }

  void OBRotor::SetDihedralAtoms(int ref[4])
  {
    // copy indexes starting from 1
    _ref.resize(4);
    for (int i = 0;i < 4;++i)
      _ref[i] = ref[i];
    _torsion.resize(4);
    // convert the indexes (start from 0, multiplied by 3) for easy access to coordinates
    _torsion[0] = (ref[0]-1)*3;
    _torsion[1] = (ref[1]-1)*3;
    _torsion[2] = (ref[2]-1)*3;
    _torsion[3] = (ref[3]-1)*3;
  }

  void OBRotor::SetRotAtoms(vector<int> &vi)
  {
    _rotatoms = vi;
  }

  //***************************************
  //**** OBRotorRules Member functions ****
  //***************************************
  OBRotorRules::OBRotorRules()
  {
    _quiet=false;
    _init = false;
    _dir = BABEL_DATADIR;
    _envvar = "BABEL_DATADIR";
    _filename = "torlib.txt";
    _subdir = "data";
    _dataptr = TorsionDefaults;
  }

  void OBRotorRules::ParseLine(const char *buffer)
  {
    int i;
    int ref[4];
    double delta;
    vector<double> vals;
    vector<string> vs;
    vector<string>::iterator j;
    char temp_buffer[BUFF_SIZE];

    if (buffer[0] == '#')
      return;
    tokenize(vs,buffer);
    if (vs.empty())
      return;

    if (EQn(buffer,"SP3-SP3",7))
      {
        _sp3sp3.clear();
        //        assert (vs.size() > 1);
        for (j = vs.begin(),++j;j != vs.end();++j)
          _sp3sp3.push_back(DEG_TO_RAD*atof(j->c_str()));
        return;
      }

    if (EQn(buffer,"SP3-SP2",7))
      {
        _sp3sp2.clear();
        //        assert(vs.size() > 1);
        for (j = vs.begin(),++j;j != vs.end();++j)
          _sp3sp2.push_back(DEG_TO_RAD*atof(j->c_str()));
        return;
      }

    if (EQn(buffer,"SP2-SP2",7))
      {
        _sp2sp2.clear();
        //        assert(vs.size() > 1);
        for (j = vs.begin(),++j;j != vs.end();++j)
          _sp2sp2.push_back(DEG_TO_RAD*atof(j->c_str()));
        return;
      }

    if (vs.size() > 5)
      {
        strncpy(temp_buffer,vs[0].c_str(), sizeof(temp_buffer) - 1);
        temp_buffer[sizeof(temp_buffer) - 1] = '\0';
        //reference atoms
        for (i = 0;i < 4;++i)
          ref[i] = atoi(vs[i+1].c_str())-1;
        //possible torsions
        vals.clear();
        delta = OB_DEFAULT_DELTA;
        for (i = 5;(unsigned)i < vs.size();++i)
          {
            if (i == (signed)(vs.size()-2) && vs[i] == "Delta")
              {
                delta = atof(vs[i+1].c_str());
                i += 2;
              }
            else
              vals.push_back(DEG_TO_RAD*atof(vs[i].c_str()));
          }

        if (vals.empty())
          {
            string err = "The following rule has no associated torsions: ";
            err += vs[0];
            obErrorLog.ThrowError(__FUNCTION__, err, obDebug);
          }
        OBRotorRule *rr = new OBRotorRule (temp_buffer,ref,vals,delta);
        if (rr->IsValid())
          _vr.push_back(rr);
        else
          delete rr;
      }

  }

  void OBRotorRules::GetRotorIncrements(OBMol &mol,OBBond *bond,
                                        int ref[4],vector<double> &vals,double &delta)
  {
    if (!_init)
      Init();

    vals.clear();
    vector<pair<int,int> > vpr;
    vpr.push_back(pair<int,int> (0,bond->GetBeginAtomIdx()));
    vpr.push_back(pair<int,int> (0,bond->GetEndAtomIdx()));

    delta = OB_DEFAULT_DELTA;

    int j;
    OBSmartsPattern *sp;
    vector<vector<int> > map;
    vector<OBRotorRule*>::iterator i;
    for (i = _vr.begin();i != _vr.end();++i)
      {
        sp = (*i)->GetSmartsPattern();
        (*i)->GetReferenceAtoms(ref);
        vpr[0].first = ref[1];
        vpr[1].first = ref[2];

        if (!sp->RestrictedMatch(mol,vpr,true))
          {
            swap(vpr[0].first,vpr[1].first);
            if (!sp->RestrictedMatch(mol,vpr,true))
              continue;
          }

        map = sp->GetMapList();
        for (j = 0;j < 4;++j)
          ref[j] = map[0][ref[j]];
        vals = (*i)->GetTorsionVals();
        delta = (*i)->GetDelta();

        OBAtom *a1,*a2,*a3,*a4,*r;
        a1 = mol.GetAtom(ref[0]);
        a4 = mol.GetAtom(ref[3]);
        if (a1->GetAtomicNum() == OBElements::Hydrogen && a4->GetAtomicNum() == OBElements::Hydrogen)
          continue; //don't allow hydrogens at both ends
        if (a1->GetAtomicNum() == OBElements::Hydrogen || a4->GetAtomicNum() == OBElements::Hydrogen) //need a heavy atom reference - can use hydrogen
          {
            bool swapped = false;
            a2 = mol.GetAtom(ref[1]);
            a3 = mol.GetAtom(ref[2]);
            if (a4->GetAtomicNum() == OBElements::Hydrogen)
              {
                swap(a1,a4);
                swap(a2,a3);
                swapped = true;
              }

            vector<OBBond*>::iterator k;
            for (r = a2->BeginNbrAtom(k);r;r = a2->NextNbrAtom(k))
              if (r->GetAtomicNum() != OBElements::Hydrogen && r != a3)
                break;

            if (!r)
              continue; //unable to find reference heavy atom
            //			cerr << "r = " << r->GetIdx() << endl;

            double t1 = mol.GetTorsion(a1,a2,a3,a4);
            double t2 = mol.GetTorsion(r,a2,a3,a4);
            double diff = t2 - t1;
            if (diff > 180.0)
              diff -= 360.0;
            if (diff < -180.0)
              diff += 360.0;
            diff *= DEG_TO_RAD;

            vector<double>::iterator m;
            for (m = vals.begin();m != vals.end();++m)
              {
                *m += diff;
                if (*m < M_PI)
                  *m += 2.0*M_PI;
                if (*m > M_PI)
                  *m -= 2.0*M_PI;
              }

            if (swapped)
              ref[3] = r->GetIdx();
            else
              ref[0] = r->GetIdx();

            /*
              mol.SetTorsion(r,a2,a3,a4,vals[0]);
              cerr << "test = " << (vals[0]-diff)*RAD_TO_DEG << ' ';
              cerr << mol.GetTorsion(a1,a2,a3,a4) <<  ' ';
              cerr << mol.GetTorsion(r,a2,a3,a4) << endl;
            */
          }

        char buffer[BUFF_SIZE];
        if (!_quiet)
          {
            snprintf(buffer,BUFF_SIZE,"%3d%3d%3d%3d %s",
                     ref[0],ref[1],ref[2],ref[3],
                     ((*i)->GetSmartsString()).c_str());
            obErrorLog.ThrowError(__FUNCTION__, buffer, obDebug);
          }
        return;
      }

    //***didn't match any rules - assign based on hybridization***
    OBAtom *a1,*a2,*a3,*a4;
    a2 = bond->GetBeginAtom();
    a3 = bond->GetEndAtom();
    vector<OBBond*>::iterator k;

    for (a1 = a2->BeginNbrAtom(k);a1;a1 = a2->NextNbrAtom(k))
      if (a1->GetAtomicNum() != OBElements::Hydrogen && a1 != a3)
        break;
    for (a4 = a3->BeginNbrAtom(k);a4;a4 = a3->NextNbrAtom(k))
      if (a4->GetAtomicNum() != OBElements::Hydrogen && a4 != a2)
        break;

    ref[0] = a1->GetIdx();
    ref[1] = a2->GetIdx();
    ref[2] = a3->GetIdx();
    ref[3] = a4->GetIdx();

    if (a2->GetHyb() == 3 && a3->GetHyb() == 3) //sp3-sp3
      {
        vals = _sp3sp3;

        if (!_quiet)
          {
            char buffer[BUFF_SIZE];
            snprintf(buffer,BUFF_SIZE,"%3d%3d%3d%3d %s",
                     ref[0],ref[1],ref[2],ref[3],"sp3-sp3");
            obErrorLog.ThrowError(__FUNCTION__, buffer, obDebug);
          }
      }
    else
      if (a2->GetHyb() == 2 && a3->GetHyb() == 2) //sp2-sp2
        {
          vals = _sp2sp2;

          if (!_quiet)
            {
              char buffer[BUFF_SIZE];
              snprintf(buffer,BUFF_SIZE,"%3d%3d%3d%3d %s",
                       ref[0],ref[1],ref[2],ref[3],"sp2-sp2");
              obErrorLog.ThrowError(__FUNCTION__, buffer, obDebug);
            }
        }
      else //must be sp2-sp3
        {
          vals = _sp3sp2;

          if (!_quiet)
            {
              char buffer[BUFF_SIZE];
              snprintf(buffer,BUFF_SIZE,"%3d%3d%3d%3d %s",
                       ref[0],ref[1],ref[2],ref[3],"sp2-sp3");
              obErrorLog.ThrowError(__FUNCTION__, buffer, obDebug);
            }
        }
  }

  OBRotorRules::~OBRotorRules()
  {
    vector<OBRotorRule*>::iterator i;
    for (i = _vr.begin();i != _vr.end();++i)
      delete (*i);
  }

#undef OB_DEFAULT_DELTA
}

//! \file rotor.cpp
//! \brief Rotate dihedral angles according to rotor rules.

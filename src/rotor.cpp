/**********************************************************************
rotor.cpp - Rotate dihedral angles according to rotor rules.
 
Copyright (C) 1998-2000 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2006 by Geoffrey R. Hutchison
Some portions Copyright (C) 2008 by Tim Vandermeersch
 
This file is part of the Open Babel project.
For more information, see <http://openbabel.sourceforge.net/>
 
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
#include <openbabel/rotor.h>

// private data headers with default parameters
#include "torlib.h"

using namespace std;
namespace OpenBabel
{

  static bool GetDFFVector(OBMol&,vector<int>&,OBBitVec&);
  static bool CompareRotor(const pair<OBBond*,int>&,const pair<OBBond*,int>&);


  //**************************************
  //**** OBRotorList Member Functions ****
  //**************************************

  bool OBRotorList::Setup(OBMol &mol)
  {
    Clear();
    FindRotors(mol);
    if (!Size())
      return(false);

    AssignTorVals(mol);

    OBRotor *rotor;
    vector<OBRotor*>::iterator i;
    for (rotor = BeginRotor(i);rotor;rotor = NextRotor(i))
      if (!rotor->Size())
        {
          unsigned int ref[4];
          rotor->GetDihedralAtoms(ref);
          char buffer[BUFF_SIZE];
          snprintf(buffer, BUFF_SIZE, "The rotor has no associated torsion values -> %d %d %d %d",
		  ref[0],ref[1],ref[2],ref[3]);
          obErrorLog.ThrowError(__FUNCTION__, buffer, obDebug);
        }

    mol.DeleteData(OBGenericDataType::TorsionData);

    return(true);
  }

  bool OBRotorList::FindRotors(OBMol &mol)
  {
    mol.FindRingAtomsAndBonds();
    vector<int> gtd;
    mol.GetGTDVector(gtd);

    obErrorLog.ThrowError(__FUNCTION__,
                          "Ran OpenBabel::FindRotors", obAuditMsg);

    OBBond *bond;
    vector<OBBond*>::iterator i;
    vector<pair<OBBond*,int> > vtmp;

    int score;
    for (bond = mol.BeginBond(i);bond;bond = mol.NextBond(i))
      if (bond->IsRotor())
        {
          if (HasFixedAtoms() && IsFixedBond(bond))
            continue;
          score = gtd[bond->GetBeginAtomIdx()-1] +
            gtd[bond->GetEndAtomIdx()-1];
          vtmp.push_back(pair<OBBond*,int> (bond,score));
        }

    sort(vtmp.begin(),vtmp.end(),CompareRotor);

    OBRotor *rotor;
    int count;
    vector<pair<OBBond*,int> >::iterator j;
    for (j = vtmp.begin(),count=0;j != vtmp.end();++j,++count)
      {
        rotor = new OBRotor;
        rotor->SetBond((*j).first);
        rotor->SetIdx(count);
        rotor->SetNumCoords(mol.NumAtoms()*3);
        m_rotors.push_back(rotor);
      }

    return(true);
  }

  bool OBRotorList::IsFixedBond(OBBond *bond)
  {
    OBAtom *a1,*a2,*a3;
    vector<OBBond*>::iterator i;

    a1 = bond->GetBeginAtom();
    a2 = bond->GetEndAtom();
    if (!_fix[a1->GetIdx()] || !_fix[a2->GetIdx()])
      return(false);

    bool isfixed=false;
    for (a3 = a1->BeginNbrAtom(i);a3;a3 = a1->NextNbrAtom(i))
      if (a3 != a2 && _fix[a3->GetIdx()])
        {
          isfixed = true;
          break;
        }

    if (!isfixed)
      return(false);

    isfixed = false;
    for (a3 = a2->BeginNbrAtom(i);a3;a3 = a2->NextNbrAtom(i))
      if (a3 != a1 && _fix[a3->GetIdx()])
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

        while (!curr.IsEmpty() && (bv&curr).Empty())
          {
            next.Clear();
            for (natom = curr.NextBit(-1);natom != curr.EndBit();natom = curr.NextBit(natom))
              {
                atom1 = mol.GetAtom(natom);
                for (bond = atom1->BeginBond(j);bond;bond = atom1->NextBond(j))
                  if (!used.BitIsOn(bond->GetNbrAtomIdx(atom1)) &&
                      !curr.BitIsOn(bond->GetNbrAtomIdx(atom1)))
                    if (!(bond->GetNbrAtom(atom1))->IsHydrogen())
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


  static double MinimumPairRMS(OBMol&,double*,double*,bool &);

  void OBRotorList::RemoveSymVals(OBMol &mol)
  {
    double *c,*c1,*c2;
    c1 = new double [mol.NumAtoms()*3];
    c2 = new double [mol.NumAtoms()*3];
    c = mol.GetCoordinates();
    bool one2one;
    double cutoff = 0.20;

    OBRotor *rotor;
    vector<OBRotor*>::iterator i;
    for (rotor = BeginRotor(i);rotor;rotor = NextRotor(i))
      {
        //look for 2-fold symmetry about a bond
        memcpy(c1,c,sizeof(double)*mol.NumAtoms()*3);
        memcpy(c2,c,sizeof(double)*mol.NumAtoms()*3);
        rotor->SetToAngle(c1,(double)(0.0*DEG_TO_RAD));
        rotor->SetToAngle(c2,(double)(180.0*DEG_TO_RAD));

        if (MinimumPairRMS(mol,c1,c2,one2one) <cutoff && !one2one)
          {
            rotor->RemoveSymTorsionValues(2);
            OBBond *bond = rotor->GetBond();

            if (!_quiet)
              {
                stringstream errorMsg;
                errorMsg << "symmetry found = " << ' ';
                errorMsg << bond->GetBeginAtomIdx() << ' ' << bond->GetEndAtomIdx() << ' ' ;
                errorMsg << "rms = " << ' ';
                errorMsg << MinimumPairRMS(mol,c1,c2,one2one) << endl;
                obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obDebug);
              }
            continue;
          }

        //look for 3-fold symmetry about a bond
        memcpy(c1,c,sizeof(double)*mol.NumAtoms()*3);
        memcpy(c2,c,sizeof(double)*mol.NumAtoms()*3);
        rotor->SetToAngle(c1,(double)(0.0*DEG_TO_RAD));
        rotor->SetToAngle(c2,(double)(120.0*DEG_TO_RAD));

        if (MinimumPairRMS(mol,c1,c2,one2one) <cutoff && !one2one)
          {
            rotor->RemoveSymTorsionValues(3);
            OBBond *bond = rotor->GetBond();

            if (!_quiet)
              {
                stringstream errorMsg;
                errorMsg << "3-fold symmetry found = " << ' ';
                errorMsg << bond->GetBeginAtomIdx() << ' ' << bond->GetEndAtomIdx() << ' ' ;
                errorMsg << "rms = " << ' ';
                errorMsg << MinimumPairRMS(mol,c1,c2,one2one) << endl;
                obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obDebug);
              }
          }
      }

    delete [] c1;
    delete [] c2;

    //pattern based duplicate removal
    unsigned int ref[4];
    vector<vector<int> > mlist;
    vector<vector<int> >::iterator k;
    vector<pair<OBSmartsPattern*,pair<int,int> > >::iterator j;
    for (j = _vsym2.begin();j != _vsym2.end();++j)
      if (j->first->Match(mol))
        {
          mlist = j->first->GetUMapList();

          for (k = mlist.begin();k != mlist.end();++k)
            for (rotor = BeginRotor(i);rotor;rotor = NextRotor(i))
              {
                rotor->GetDihedralAtoms(ref);
                if (((*k)[j->second.first] == ref[1] && (*k)[j->second.second] == ref[2]) ||
                    ((*k)[j->second.first] == ref[2] && (*k)[j->second.second] == ref[1]))
                  {
                    rotor->RemoveSymTorsionValues(2);
                    if (!_quiet)
                      {
                        stringstream errorMsg;
                        errorMsg << "2-fold pattern-based symmetry found = " << ' ';
                        errorMsg << ref[1] << ' ' << ref[2] << endl;
                        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obDebug);
                      }
                  }
              }
        }

    for (j = _vsym3.begin();j != _vsym3.end();++j)
      if (j->first->Match(mol))
        {
          mlist = j->first->GetUMapList();

          for (k = mlist.begin();k != mlist.end();++k)
            for (rotor = BeginRotor(i);rotor;rotor = NextRotor(i))
              {
                rotor->GetDihedralAtoms(ref);
                if (((*k)[j->second.first] == ref[1] && (*k)[j->second.second] == ref[2]) ||
                    ((*k)[j->second.first] == ref[2] && (*k)[j->second.second] == ref[1]))
                  {
                    rotor->RemoveSymTorsionValues(3);
                    if (!_quiet)
                      {
                        stringstream errorMsg;
                        errorMsg << "3-fold pattern-based symmetry found = " << ' ';
                        errorMsg << ref[1] << ' ' << ref[2] << endl;
                        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obDebug);
                      }
                  }
              }
        }
  }

  static double MinimumPairRMS(OBMol &mol,double *a,double *b,bool &one2one)
  {
    unsigned int i,j,k=0;
    double min,tmp,d_2 = 0.0;
    OBBitVec bset;
    one2one = true;
    vector<OBAtom*> _atom;
    _atom.resize(mol.NumAtoms());
    for (i = 0;i < mol.NumAtoms();++i)
      _atom[i] = mol.GetAtom(i+1);

    for (i = 0;i < mol.NumAtoms();++i)
      {
        min = 10E10;
        for (j = 0;j < mol.NumAtoms();++j)
          if ((_atom[i])->GetAtomicNum() == (_atom[j])->GetAtomicNum() &&
              (_atom[i])->GetHyb()       == (_atom[j])->GetHyb())
            if (!bset[j])
              {
                tmp = SQUARE(a[3*i]-b[3*j]) +
                  SQUARE(a[3*i+1]-b[3*j+1]) +
                  SQUARE(a[3*i+2]-b[3*j+2]);
                if (tmp < min)
                  {
                    k = j;
                    min = tmp;
                  }
              }

        if (i != j)
          one2one = false;
        bset.SetBitOn(k);
        d_2 += min;
      }

    d_2 /= (double)mol.NumAtoms();

    return(sqrt(d_2));
  }

  bool OBRotorList::AssignTorVals(OBMol &mol)
  {
    OBBond *bond;
    OBRotor *rotor;
    vector<OBRotor*>::iterator i;

    unsigned int ref[4];
    vector<double> vals;
    vector<int> itmp1;
    vector<int>::iterator j;
    for (i = m_rotors.begin();i != m_rotors.end();++i)
      {
        rotor=*i;
        bond = rotor->GetBond();
        _rr.GetRotorIncrements(mol, bond, ref, vals);
        rotor->SetTorsionValues(vals);

        mol.FindChildren(itmp1,ref[1],ref[2]);
        if (itmp1.size()+1 > mol.NumAtoms()/2)
          {
            itmp1.clear();
            mol.FindChildren(itmp1,ref[2],ref[1]);
            swap(ref[0],ref[3]);
            swap(ref[1],ref[2]);
          }

        rotor->SetRotAtoms(itmp1);
        rotor->SetDihedralAtoms(ref);
      }

    return(true);
  }

  bool OBRotorList::SetRotAtoms(OBMol &mol)
  {
    OBRotor *rotor;
    vector<int> rotatoms;
    vector<OBRotor*>::iterator i;
    vector<int>::iterator j;

    unsigned int ref[4];
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

        rotor->SetRotAtoms(rotatoms);
        rotor->SetDihedralAtoms(ref);
      }

    return(true);
  }

  void OBRotorList::SetRotAtomsByFix(OBMol &mol)
  {
    unsigned int ref[4];
    OBRotor *rotor;
    vector<int> rotatoms;
    vector<int>::iterator j;
    vector<OBRotor*>::iterator i;

    GetDFFVector(mol,_dffv,_fix);

    for (rotor = BeginRotor(i);rotor;rotor = NextRotor(i))
      {
        rotatoms.clear();
        rotor->GetDihedralAtoms(ref);
        if (_fix[ref[1]] && _fix[ref[2]])
          {
            if (!_fix[ref[0]])
              {
                swap(ref[0],ref[3]);
                swap(ref[1],ref[2]);
                mol.FindChildren(rotatoms,ref[1],ref[2]);
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
              rotor->SetRotAtoms(rotatoms);
              rotor->SetDihedralAtoms(ref);
            }
      }
  }

  OBRotorList::OBRotorList()
  {
    m_rotors.clear();
    _quiet=false;

    //para-disub benzene
    OBSmartsPattern *sp;
    sp = new OBSmartsPattern;
    sp->Init("*c1[cD2][cD2]c(*)[cD2][cD2]1");
    _vsym2.push_back(pair<OBSmartsPattern*,pair<int,int> > (sp,pair<int,int> (0,1)));

    //piperidine amide
    sp = new OBSmartsPattern;
    sp->Init("O=CN1[CD2][CD2][CD2][CD2][CD2]1");
    _vsym2.push_back(pair<OBSmartsPattern*,pair<int,int> > (sp,pair<int,int> (1,2)));

    //terminal phosphate
    sp = new OBSmartsPattern;
    sp->Init("[#8D2][#15,#16](~[#8D1])(~[#8D1])~[#8D1]");
    _vsym3.push_back(pair<OBSmartsPattern*,pair<int,int> > (sp,pair<int,int> (0,1)));

  }

  OBRotorList::~OBRotorList()
  {
    vector<OBRotor*>::iterator i;
    for (i = m_rotors.begin();i != m_rotors.end();++i)
      delete *i;

    vector<pair<OBSmartsPattern*,pair<int,int> > >::iterator j;
    for (j = _vsym2.begin();j != _vsym2.end();++j)
      delete j->first;

    for (j = _vsym3.begin();j != _vsym3.end();++j)
      delete j->first;
  }

  void OBRotorList::Clear()
  {
    vector<OBRotor*>::iterator i;
    for (i = m_rotors.begin();i != m_rotors.end();++i)
      delete *i;
    m_rotors.clear();
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

  OBRotor::~OBRotor()
  {
  }

  void OBRotor::SetToAngle(double *c, double setang)
  {
    SetTorsion(c, m_ref, setang, m_rotatoms);
  }

  void OBRotor::RemoveSymTorsionValues(int fold)
  {
    vector<double>::iterator i;
    vector<double> tv;

    // only one torsion value
    if (m_torsions.size() == 1)
      return;

    for (i = m_torsions.begin();i != m_torsions.end();++i) {
      if (*i >= 0.0) {
        // for 2-fold symmetry, we only keep torsion values < 180
        if (fold == 2 && *i < DEG_TO_RAD*180.0)
          tv.push_back(*i);
        // for 3-fold symmetry, we only keep torsion values < 120
        if (fold == 3 && *i < DEG_TO_RAD*120.0)
          tv.push_back(*i);
      }
    }

    if (tv.empty())
      return;

    m_torsions = tv;
  }

  void OBRotor::SetDihedralAtoms(unsigned int ref[4])
  {
    for (int i = 0;i < 4;++i)
      m_ref[i] = ref[i];

    m_cidx.resize(4);
    m_cidx[0] = (ref[0] - 1) * 3;
    m_cidx[1] = (ref[1] - 1) * 3;
    m_cidx[2] = (ref[2] - 1) * 3;
    m_cidx[3] = (ref[3] - 1) * 3;
  }

  //***************************************
  //**** OBRotorRules Member functions ****
  //***************************************
  
  OBRotorRules::OBRotorRules()
  {
    m_quiet=false;
    _init = false;
    _dir = BABEL_DATADIR;
    _envvar = "BABEL_DATADIR";
    _filename = "torlib.txt";
    _subdir = "data";
    _dataptr = TorsionDefaults;
  }

  void OBRotorRules::ParseLine(const char *buffer)
  {
    unsigned int ref[4];
    vector<double> vals;
    vector<string> vs;
    vector<string>::iterator j;
    char temp_buffer[BUFF_SIZE];

    // ignore comments
    if (buffer[0] == '#')
      return;

    tokenize(vs,buffer);
    // ignore empty lines
    if (vs.empty())
      return;

    // default sp3-sp3 rule
    if (EQn(buffer, "SP3-SP3", 7)) {
      m_sp3sp3.clear();
      for (j = vs.begin(),j++;j != vs.end();++j)
        m_sp3sp3.push_back(DEG_TO_RAD*atof(j->c_str()));
      
      return;
    }

    // default sp3-sp2 rule
    if (EQn(buffer, "SP3-SP2", 7)) {
      m_sp3sp2.clear();
      for (j = vs.begin(),j++;j != vs.end();++j)
        m_sp3sp2.push_back(DEG_TO_RAD*atof(j->c_str()));
      
      return;
    }

    // default sp3-sp2 rule
    if (EQn(buffer,"SP2-SP2",7)) {
      m_sp2sp2.clear();
      for (j = vs.begin(),j++;j != vs.end();++j)
        m_sp2sp2.push_back(DEG_TO_RAD*atof(j->c_str()));
      
      return;
    }

    // we need at least a smarts pattern + 4 atom refs + 1 angle
    if (vs.size() > 5) {
      strncpy(temp_buffer, vs[0].c_str(), sizeof(temp_buffer) - 1);
      temp_buffer[sizeof(temp_buffer) - 1] = '\0';
        
      //reference atoms
      for (int i = 0;i < 4;++i)
        ref[i] = atoi(vs[i+1].c_str())-1;
      
      //possible torsions
      vals.clear();
      for (unsigned int i = 5; i < vs.size(); ++i) {
        if (i == (vs.size()-2) && vs[i] == "Delta") {
          // $SMARTS $a1 $a2 $a3 $a4 Delta $delta
          double delta = atof(vs[i+1].c_str());
          
          OBRotorRule *rr = new OBRotorRule (temp_buffer, ref, delta);
          if (rr->IsValid())
            m_vr.push_back(rr);
          else
            delete rr;
          
          return; // return here

        } else {
          // $SMARTS $a1 $a2 $a3 $a4 $tor1 $tor2 $tor3 ...
          vals.push_back(DEG_TO_RAD * atof(vs[i].c_str()));
        }
      }

      if (vals.empty()) {
        string err = "The following rule has no associated torsions: ";
        err += vs[0];
        obErrorLog.ThrowError(__FUNCTION__, err, obDebug);
      }
      
      // creat a new OBRotorRule and save it if valid
      OBRotorRule *rr = new OBRotorRule (temp_buffer, ref, vals);
      if (rr->IsValid())
        m_vr.push_back(rr);
      else
        delete rr;
    }
  }

  void OBRotorRules::GetRotorIncrements(OBMol &mol, OBBond *bond, unsigned int ref[4], vector<double> &vals)
  {
    if (!_init)
      Init();

    vals.clear();
    vector<pair<int, int> > vpr;
    vpr.push_back(pair<int, int>(0, bond->GetBeginAtomIdx()));
    vpr.push_back(pair<int, int>(0, bond->GetEndAtomIdx()));

    int j;
    OBSmartsPattern *sp;
    vector<vector<int> > map;
    vector<OBRotorRule*>::iterator i;
    for (i = m_vr.begin(); i != m_vr.end(); ++i) { // for each OBRotorRule
      sp = (*i)->GetSmartsPattern();
      (*i)->GetReferenceAtoms(ref);
      vpr[0].first = ref[1];
      vpr[1].first = ref[2];

      if (!sp->RestrictedMatch(mol, vpr, true)) {
        swap(vpr[0].first, vpr[1].first);
        if (!sp->RestrictedMatch(mol, vpr, true))
          continue;
      }

      map = sp->GetMapList();
      for (j = 0;j < 4;++j)
        ref[j] = map[0][ref[j]];
      vals = (*i)->GetTorsionVals();

      OBAtom *a1,*a2,*a3,*a4,*r;
      a1 = mol.GetAtom(ref[0]);
      a4 = mol.GetAtom(ref[3]);
      if (a1->IsHydrogen() && a4->IsHydrogen())
        continue; //don't allow hydrogens at both ends
      if (a1->IsHydrogen() || a4->IsHydrogen()) { //need a heavy atom reference - can use hydrogen
        bool swapped = false;
        a2 = mol.GetAtom(ref[1]);
        a3 = mol.GetAtom(ref[2]);
        if (a4->IsHydrogen()) {
          swap(a1, a4);
          swap(a2, a3);
          swapped = true;
        }

        vector<OBBond*>::iterator k;
        for (r = a2->BeginNbrAtom(k);r;r = a2->NextNbrAtom(k))
          if (!r->IsHydrogen() && r != a3)
            break;

        if (!r)
          continue; //unable to find reference heavy atom
                    //cerr << "r = " << r->GetIdx() << endl;

        double t1 = mol.GetTorsion(a1,a2,a3,a4);
        double t2 = mol.GetTorsion(r,a2,a3,a4);
        double diff = t2 - t1;
        if (diff > 180.0)
          diff -= 360.0;
        if (diff < -180.0)
          diff += 360.0;
        diff *= DEG_TO_RAD;

        vector<double>::iterator m;
        for (m = vals.begin(); m != vals.end(); ++m) {
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
      } // if (a1->IsHydrogen() || a4->IsHydrogen()) 

      char buffer[BUFF_SIZE];
      if (!m_quiet) {
        snprintf(buffer,BUFF_SIZE,"%3d%3d%3d%3d %s",
            ref[0], ref[1], ref[2], ref[3],
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
      if (!a1->IsHydrogen() && a1 != a3)
        break;
    for (a4 = a3->BeginNbrAtom(k);a4;a4 = a3->NextNbrAtom(k))
      if (!a4->IsHydrogen() && a4 != a2)
        break;

    ref[0] = a1->GetIdx();
    ref[1] = a2->GetIdx();
    ref[2] = a3->GetIdx();
    ref[3] = a4->GetIdx();

    if (a2->GetHyb() == 3 && a3->GetHyb() == 3) { //sp3-sp3
      vals = m_sp3sp3;

      if (!m_quiet) {
        char buffer[BUFF_SIZE];
        snprintf(buffer, BUFF_SIZE, "%3d%3d%3d%3d %s",
            ref[0], ref[1], ref[2], ref[3], "sp3-sp3");
        obErrorLog.ThrowError(__FUNCTION__, buffer, obDebug);
      }
    } else if (a2->GetHyb() == 2 && a3->GetHyb() == 2) { //sp2-sp2
      vals = m_sp2sp2;

      if (!m_quiet) {
        char buffer[BUFF_SIZE];
        snprintf(buffer, BUFF_SIZE, "%3d%3d%3d%3d %s",
            ref[0], ref[1], ref[2], ref[3], "sp2-sp2");
        obErrorLog.ThrowError(__FUNCTION__, buffer, obDebug);
      }
    } else { //must be sp2-sp3
      vals = m_sp3sp2;

      if (!m_quiet) {
        char buffer[BUFF_SIZE];
        snprintf(buffer, BUFF_SIZE, "%3d%3d%3d%3d %s",
            ref[0], ref[1], ref[2], ref[3], "sp2-sp3");
        obErrorLog.ThrowError(__FUNCTION__, buffer, obDebug);
      }
    }
  }

  OBRotorRules::~OBRotorRules()
  {
    // delte the rotor rules
    vector<OBRotorRule*>::iterator i;
    for (i = m_vr.begin();i != m_vr.end();++i)
      delete (*i);
  }

}

//! \file rotor.cpp
//! \brief Rotate dihedral angles according to rotor rules.

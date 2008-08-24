/**********************************************************************
phmodel.cpp - Read pH rules and assign charges.
 
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2005 by Geoffrey R. Hutchison
 
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
#include <openbabel/phmodel.h>

// private data header with default parameters
#include "phmodeldata.h"

#ifdef WIN32
#pragma warning (disable : 4786)
#endif

using namespace std;

namespace OpenBabel
{

  // Global OBPhModel for assigning formal charges and hydrogen addition rules
  OBPhModel phmodel;
  extern OBAtomTyper atomtyper;

  OBPhModel::OBPhModel()
  {
    _init = false;
    _dir = BABEL_DATADIR;
    _envvar = "BABEL_DATADIR";
    _filename = "phmodel.txt";
    _subdir = "data";
    _dataptr = PhModelData;
  }

OBPhModel::~OBPhModel()
{
    vector<OBChemTsfm*>::iterator k;
    for (k = _vtsfm.begin();k != _vtsfm.end();++k)
        delete *k;

    vector<pair<OBSmartsPattern*,vector<double> > >::iterator m;
    for (m = _vschrg.begin();m != _vschrg.end();++m)
        delete m->first;
}

void OBPhModel::ParseLine(const char *buffer)
{
    vector<string> vs;
    OBSmartsPattern *sp;

    if (buffer[0] == '#')
        return;

    if (EQn(buffer,"TRANSFORM",7))
    {
        tokenize(vs,buffer);
        if (vs.empty() || vs.size() < 5)
	  {
	    obErrorLog.ThrowError(__FUNCTION__, " Could not parse line in phmodel table from phmodel.txt", obInfo);
            return;
	  }

        OBChemTsfm *tsfm = new OBChemTsfm;
        if (!tsfm->Init(vs[1],vs[3]))
        {
            delete tsfm;
            tsfm = NULL;
	    obErrorLog.ThrowError(__FUNCTION__, " Could not parse line in phmodel table from phmodel.txt", obInfo);
            return;
        }

        _vtsfm.push_back(tsfm);
        _vpKa.push_back(atof(vs[4].c_str()));
    }
    else if (EQn(buffer,"SEEDCHARGE",10))
    {
        tokenize(vs,buffer);
        if (vs.empty() || vs.size() < 2)
	  {
	    obErrorLog.ThrowError(__FUNCTION__, " Could not parse line in phmodel table from phmodel.txt", obInfo);
            return;
	  }

        sp = new OBSmartsPattern;
        if (!sp->Init(vs[1]) || (vs.size()-2) != sp->NumAtoms())
        {
            delete sp;
            sp = NULL;
	    obErrorLog.ThrowError(__FUNCTION__, " Could not parse line in phmodel table from phmodel.txt", obInfo);
            return;
        }

        vector<double> vf;
        vector<string>::iterator i;
        for (i = vs.begin()+2;i != vs.end();++i)
            vf.push_back(atof((char*)i->c_str()));

        _vschrg.push_back(pair<OBSmartsPattern*,vector<double> > (sp,vf));
    }
}

void OBPhModel::AssignSeedPartialCharge(OBMol &mol)
{
    if (!_init)
        Init();

    mol.SetPartialChargesPerceived();
    if (!mol.AutomaticPartialCharge())
        return;

    vector<pair<OBSmartsPattern*,vector<double> > >::iterator i;
    for (i = _vschrg.begin();i != _vschrg.end();++i)
        if (i->first->Match(mol))
        {
            _mlist = i->first->GetUMapList();
            unsigned int k;
            vector<vector<int> >::iterator j;

            for (j = _mlist.begin();j != _mlist.end();++j)
                for (k = 0;k < j->size();++k)
                    mol.GetAtom((*j)[k])->SetPartialCharge(i->second[k]);
        }
}

void OBPhModel::CorrectForPH(OBMol &mol, double pH)
{
    if (!_init)
        Init();
    if (mol.IsCorrectedForPH())
        return;
    if (mol.GetDimension() > 0 && !mol.AutomaticFormalCharge())
        return;

    mol.SetCorrectedForPH();

    obErrorLog.ThrowError(__FUNCTION__,
                          "Ran OpenBabel::CorrectForPH", obAuditMsg);

    mol.DeleteHydrogens();

    for (unsigned int i = 0; i < _vtsfm.size(); ++i) {
      
      if (_vpKa[i] > 1E+9) { 
	// always apply when pKa is > 1e+9
        _vtsfm[i]->Apply(mol);
      } else {
        // 10^(pKa - pH) = [HA] / [A-]
        //
        // > 1 : [HA] > [A-]
        // < 1 : [HA] < [A-]
        if (_vtsfm[i]->IsAcid()) {
          //cout << "IsAcid == " << _vtsfm[i]->IsAcid() << endl;
          //cout << "pKa == " << _vpKa[i] << endl;
          //cout << "pow(10, _vpKa[i] - pH) == " << pow(10, _vpKa[i] - pH) << endl;
          if (pow(10, _vpKa[i] - pH) < 1.0) {
            //cout << "APPLY!!" << endl;
            _vtsfm[i]->Apply(mol);
	  }
        }

        // 10^(pKa - pH) = [BH+] / [B:]
        //
        // > 1 : [BH+] > [B:]
        // < 1 : [BH+] < [B:]
        if (_vtsfm[i]->IsBase()) {
          //cout << "IsBase == " << _vtsfm[i]->IsBase() << endl;
          //cout << "pKa == " << _vpKa[i] << endl;
          //cout << "pow(10, _vpKa[i] - pH) == " << pow(10, _vpKa[i] - pH) << endl;
          if (pow(10, _vpKa[i] - pH) > 1.0) {
            //cout << "APPLY!!" << endl;
            _vtsfm[i]->Apply(mol);
	  }
        }
      }
    }

    atomtyper.CorrectAromaticNitrogens(mol);
}


  // Portions of this documentation adapted from the JOELib docs, written by
  // Joerg Wegner
  /** @class OBChemTsfm phmodel.h <openbabel/phmodel.h>
      @brief SMARTS based structural modification (chemical transformation)
 
      Transformation of chemical structures can be used for pH value correction
      (i.e. via OBPhModel and OBMol::CorrectForPH()). The OBChemTsfm class
      defines SMARTS based TRANSFORM patterns to delete atoms, change atom types,
      atom formal charges, and bond types.

      For storing and converting chemical reaction files, use the OBReaction class.
   */
  bool OBChemTsfm::Init(std::string &bgn, std::string &end)
  {
    if (!m_bgn.Init(bgn))
      return false;
    
    if (!end.empty())
        if (!m_end.Init(end))
            return false;

    //find atoms to be deleted
    unsigned int i, j;
    int vb;
    bool found;
    for (i = 0;i < m_bgn.NumAtoms();++i) {
      if ((vb = m_bgn.GetVectorBinding(i))) {
        found = false;
        for (j = 0;j < m_end.NumAtoms();++j) {
          if (vb == m_end.GetVectorBinding(j)) {
            found = true;
            break;
          }
        }

        if (!found)
	  m_vadel.push_back(i);
      }
    }

    //find elements to be changed
    int ele;
    for (i = 0;i < m_bgn.NumAtoms();++i) {
      if ((vb = m_bgn.GetVectorBinding(i)) != 0) {
        ele = m_bgn.GetAtomicNum(i);
        
        for (j = 0;j < m_end.NumAtoms();++j) {
          if (vb == m_end.GetVectorBinding(j)) {
            
            if (ele != m_end.GetAtomicNum(j)) {
              m_vele.push_back(pair<int,int> (i,m_end.GetAtomicNum(j)));
              break;
            }
          }
        }

      }
    }

    //find charges to modify
    int chrg;
    for (i = 0;i < m_bgn.NumAtoms();++i) {
      if ((vb = m_bgn.GetVectorBinding(i))) {
        chrg = m_bgn.GetCharge(i);
        
        for (j = 0;j < m_end.NumAtoms();++j) {
          if (vb == m_end.GetVectorBinding(j)) {
            if (chrg != m_end.GetCharge(j))
              m_vchrg.push_back(pair<int,int> (i,m_end.GetCharge(j)));
          }
        }
      
      }
    }

    //find bonds to be modified
    //find bonds to be modified
    int bsrc,bdst,bord,bvb1,bvb2;
    int esrc,edst,eord,evb1,evb2;

    for (i = 0;i < m_bgn.NumBonds();++i) {
      m_bgn.GetBond(bsrc,bdst,bord,i);
      bvb1 = m_bgn.GetVectorBinding(bsrc);
      bvb2 = m_bgn.GetVectorBinding(bdst);
      if (!bvb1 || !bvb2)
        continue;

      for (j = 0;j < m_end.NumBonds();++j) {
        m_end.GetBond(esrc,edst,eord,j);
        evb1 = m_end.GetVectorBinding(esrc);
        evb2 = m_end.GetVectorBinding(edst);
        if ((bvb1 == evb1 && bvb2 == evb2) || (bvb1 == evb2 && bvb2 == evb1)) {
          if (bord == eord)
            break; //nothing to modify if bond orders identical
          m_vbond.push_back(pair<pair<int,int>,int> (pair<int,int> (bsrc,bdst),eord));
          break;
        }
      }
    }

    //make sure there is some kind of transform to do here
    if (m_vadel.empty() && m_vchrg.empty() && m_vbond.empty())
      return false;

    return true;
  }

bool OBChemTsfm::Apply(OBMol &mol)
{
    if (!m_bgn.Match(mol))
        return(false);

    vector<vector<int> > mlist = m_bgn.GetUMapList();

    obErrorLog.ThrowError(__FUNCTION__,
                          "Ran OpenBabel::OBChemTransform", obAuditMsg);

    if (!m_vchrg.empty()) //modify charges
    {
        vector<vector<int> >::iterator i;
        vector<pair<int,int> >::iterator j;

        for (i = mlist.begin();i != mlist.end();++i)
            for (j = m_vchrg.begin();j != m_vchrg.end();++j)
                if (j->first < (signed)i->size()) //goof proofing
                    mol.GetAtom((*i)[j->first])->SetFormalCharge(j->second);

        mol.UnsetImplicitValencePerceived();
    }

    if (!m_vbond.empty()) //modify bond orders
    {
        OBBond *bond;
        vector<vector<int> >::iterator i;
        vector<pair<pair<int,int>,int> >::iterator j;
        for (i = mlist.begin();i != mlist.end();++i)
            for (j = m_vbond.begin();j != m_vbond.end();++j)
            {
                bond = mol.GetBond((*i)[j->first.first],(*i)[j->first.second]);
                if (!bond)
                {
		  obErrorLog.ThrowError(__FUNCTION__, "unable to find bond", obDebug);
                    continue;
                }

                bond->SetBondOrder(j->second);
            }
    }

    if (!m_vadel.empty() || !m_vele.empty()) //delete atoms and change elements
    {
        vector<int>::iterator j;
        vector<vector<int> >::iterator i;

        if (!m_vele.empty())
        {
            vector<pair<int,int> >::iterator k;
            for (i = mlist.begin();i != mlist.end();++i)
                for (k = m_vele.begin();k != m_vele.end();++k)
                    mol.GetAtom((*i)[k->first])->SetAtomicNum(k->second);
        }

        //make sure same atom isn't deleted twice
        vector<bool> vda;
        vector<OBAtom*> vdel;
        vda.resize(mol.NumAtoms()+1,false);
        for (i = mlist.begin();i != mlist.end();++i)
            for (j = m_vadel.begin();j != m_vadel.end();++j)
                if (!vda[(*i)[*j]])
                {
                    vda[(*i)[*j]] = true;
                    vdel.push_back(mol.GetAtom((*i)[*j]));
                }

        vector<OBAtom*>::iterator k;
        for (k = vdel.begin();k != vdel.end();++k)
            mol.DeleteAtom((OBAtom*)*k);
    }

    return(true);
}

bool OBChemTsfm::IsAcid()
{
  //cout << m_bgn.GetSMARTS() << " >> " << m_end.GetSMARTS() << endl;
  //cout << "  NumAtoms = " << m_end.NumAtoms() << endl;

  if (m_bgn.NumAtoms() > m_end.NumAtoms())  // O=CO[#1:1] >> O=CO
    return true;

  for (unsigned int i = 0; i < m_end.NumAtoms(); ++i) {
    //cout << "    m_end(" << i << ")  " << m_end.GetCharge(i) << endl;
    //cout << "    m_bgn(" << i << ")  " << m_bgn.GetCharge(i) << endl;
    if (m_end.GetCharge(i) < 0)
      return true;
  }
  
  return false;
}

bool OBChemTsfm::IsBase()
{
  //cout << m_bgn.GetSMARTS() << " >> " << m_end.GetSMARTS() << endl;
  //cout << "  NumAtoms = " << m_end.NumAtoms() << endl;

  for (unsigned int i = 0; i < m_end.NumAtoms(); ++i) {
    //cout << "    m_end(" << i << ")  " << m_end.GetCharge(i) << endl;
    //cout << "    m_bgn(" << i << ")  " << m_bgn.GetCharge(i) << endl;
    if (m_end.GetCharge(i) > 0)
      return true;
  }

  return false;
}

} //namespace OpenBabel

//! \file phmodel.cpp
//! \brief Read pH rules and assign charges.

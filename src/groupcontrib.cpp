/**********************************************************************
logp.cpp - Handle logP prediction algorithms.
 
Copyright (C) 2007 by Tim Vandermeersch
 
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

#include <openbabel/groupcontrib.h>

using namespace std;

namespace OpenBabel
{
  OBGroupContrib::OBGroupContrib()
  {
  }

  OBGroupContrib::~OBGroupContrib()
  {
  }

  bool OBGroupContrib::ParseFile(const char *filename)
  {
    OBSmartsPattern *sp;
    
    // open data/psa.txt
    ifstream ifs;

    if (OpenDatafile(ifs, filename).length() == 0) {
      obErrorLog.ThrowError(__FUNCTION__, " Could not find contribution data file.", obError);
      return false;
    }

    vector<string> vs;
    bool heavy = false;
    
    char buffer[80];
    while (ifs.getline(buffer, 80)) {
      if (EQn(buffer, "#", 1)) continue;
      if (EQn(buffer, ";heavy", 6))
        heavy = true;
      else if (EQn(buffer, ";", 1)) continue;

	
      tokenize(vs, buffer);
      if (vs.size() < 2)
        continue;
      
      sp = new OBSmartsPattern;
      if (sp->Init(vs[0])) {
        if (heavy)
          _contribsHeavy.push_back(pair<OBSmartsPattern*, double> (sp, atof(vs[1].c_str())));
	else
          _contribsHydrogen.push_back(pair<OBSmartsPattern*, double> (sp, atof(vs[1].c_str())));
      } else {
        delete sp;
        sp = NULL;
        obErrorLog.ThrowError(__FUNCTION__, " Could not parse SMARTS from contribution data file", obInfo);
	return false;
      }
    }

    return true;
  }
  
 
  double OBGroupContrib::GroupContributions(OBMol &mol)
  {
    vector<vector<int> > _mlist; //!< match list for atom typing
    vector<vector<int> >::iterator j;
    vector<pair<OBSmartsPattern*, double> >::iterator i;
    
    vector<double> atomValues(mol.NumAtoms(), 0.0);

    // atom contributions
    //cout << "atom contributions:" << endl;
    for (i = _contribsHeavy.begin();i != _contribsHeavy.end();++i) {
      if (i->first->Match(mol)) {
        _mlist = i->first->GetMapList();
        for (j = _mlist.begin();j != _mlist.end();++j) {
	  atomValues[(*j)[0] - 1] = i->second;
	  //cout << (*j)[0] << " = " << i->first->GetSMARTS() << " : " << i->second << endl;
        }
      }
    }
    
    vector<double> hydrogenValues(mol.NumAtoms(), 0.0);
    //hydrogenValues.resize(mol.NumAtoms());   
    
    // hydrogen contributions
    //cout << "hydrogen contributions:" << endl;
    for (i = _contribsHydrogen.begin();i != _contribsHydrogen.end();++i) {
      if (i->first->Match(mol)) {
        _mlist = i->first->GetMapList();
        for (j = _mlist.begin();j != _mlist.end();++j) {
	  int Hcount = mol.GetAtom((*j)[0])->GetValence() - mol.GetAtom((*j)[0])->GetHvyValence();
	  hydrogenValues[(*j)[0] - 1] = i->second * Hcount;
	  //cout << (*j)[0] << " = " << i->first->GetSMARTS() << " : " << i->second << endl;
        }
      }
    }

    // total atomic and hydrogen contribution
    double total = 0.0;

    for (int index = 0; index < mol.NumAtoms(); index++) {
      if (mol.GetAtom(index+1)->IsHydrogen())
        continue;

      total += atomValues[index];
      total += hydrogenValues[index];
    }
   
    /*
    FOR_ATOMS_OF_MOL (a, mol)
      cout << "hydrogens on atom " << a->GetIdx() << ": " << a->GetValence() - a->GetHvyValence() << endl;
    for (int index = 0; index < mol.NumAtoms(); index++)
      cout << "atom " << index << ": " << atomValues[index] << endl;
    for (int index = 0; index < mol.NumAtoms(); index++)
      cout << "hydrogen " << index << ": " << hydrogenValues[index] << endl;
    */

    return total;
  }
  
  OBLogP::OBLogP()
  {
    ParseFile("logp.txt");
  }

  OBLogP::~OBLogP()
  {
  }
  
  double OBLogP::Predict(OBMol &mol)
  {
    return GroupContributions(mol);
  }

  OBPSA::OBPSA()
  {
    ParseFile("psa.txt");
  }

  OBPSA::~OBPSA()
  {
  }
  
  double OBPSA::Predict(OBMol &mol)
  {
    return GroupContributions(mol);
  }

  OBMR::OBMR()
  {
    ParseFile("mr.txt");
  }

  OBMR::~OBMR()
  {
  }
  
  double OBMR::Predict(OBMol &mol)
  {
    return GroupContributions(mol);
  }

} // end namespace OpenBabel

//! \file logp.cpp
//! \brief Handle logP prediction algorithms.

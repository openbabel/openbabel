/**********************************************************************
groupcontrib.cpp - Handle logP, PSA, MR, and other group-based predictions

Copyright (C) 2007      by Tim Vandermeersch
              2001-2007 by Stephen Jelfs
              2001-2007 by Joerg Kurt Wegner, me@cheminformatics.eu

Original version: JOELib2, http://joelib.sf.net
 
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
  /** \class OBGroupContrib groupcontrib.h <openbabel/groupcontrib.h>
      \brief Handle group contribution algorithms.
 
      This is the base class for calculations that use the JOELib2 contribution 
      algorithm. See the derived OBPSA, OBLogP, OBMR classes for more 
      information on how to use these classes.
    */

  OBGroupContrib::OBGroupContrib()
  {
  }

  OBGroupContrib::~OBGroupContrib()
  {
  }

  bool OBGroupContrib::ParseFile(const char *filename)
  {
    OBSmartsPattern *sp;
    
    // open data file
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
    vector<vector<int> > _mlist; // match list for atom typing
    vector<vector<int> >::iterator j;
    vector<pair<OBSmartsPattern*, double> >::iterator i;
    
    vector<double> atomValues(mol.NumAtoms(), 0.0);

    OBMol tmpmol;
    tmpmol = mol;

    tmpmol.ConvertDativeBonds();

    // atom contributions
    //cout << "atom contributions:" << endl;
    for (i = _contribsHeavy.begin();i != _contribsHeavy.end();++i) {
      if (i->first->Match(tmpmol)) {
        _mlist = i->first->GetMapList();
        for (j = _mlist.begin();j != _mlist.end();++j) {
	  atomValues[(*j)[0] - 1] = i->second;
	  //cout << (*j)[0] << " = " << i->first->GetSMARTS() << " : " << i->second << endl;
        }
      }
    }
    
    vector<double> hydrogenValues(tmpmol.NumAtoms(), 0.0);
    //hydrogenValues.resize(tmpmol.NumAtoms());   
    
    // hydrogen contributions
    //cout << "hydrogen contributions:" << endl;
    for (i = _contribsHydrogen.begin();i != _contribsHydrogen.end();++i) {
      if (i->first->Match(tmpmol)) {
        _mlist = i->first->GetMapList();
        for (j = _mlist.begin();j != _mlist.end();++j) {
	  int Hcount = tmpmol.GetAtom((*j)[0])->GetValence() - tmpmol.GetAtom((*j)[0])->GetHvyValence();
	  hydrogenValues[(*j)[0] - 1] = i->second * Hcount;
	  //cout << (*j)[0] << " = " << i->first->GetSMARTS() << " : " << i->second << endl;
        }
      }
    }

    // total atomic and hydrogen contribution
    double total = 0.0;

    for (int index = 0; index < tmpmol.NumAtoms(); index++) {
      if (tmpmol.GetAtom(index+1)->IsHydrogen())
        continue;

      total += atomValues[index];
      total += hydrogenValues[index];
    }
   
    /*
    FOR_ATOMS_OF_MOL (a, tmpmol)
      cout << "hydrogens on atom " << a->GetIdx() << ": " << a->GetValence() - a->GetHvyValence() << endl;
    for (int index = 0; index < tmpmol.NumAtoms(); index++)
      cout << "atom " << index << ": " << atomValues[index] << endl;
    for (int index = 0; index < tmpmol.NumAtoms(); index++)
      cout << "hydrogen " << index << ": " << hydrogenValues[index] << endl;
    */

    return total;
  }
  
  /** \class OBLogP groupcontrib.h <openbabel/groupcontrib.h>
      \brief calculate the LogP (octanol/water partition coefficient).
 
      This class uses the JOELib2 group contribution algorithm to calculate 
      the logP (octanol/water partition coefficient) of a molecule.

      example:
      \code
      #include <openbabel/groupcontrib.h>
      #include <openbabel/mol.h>

      OBMol mol;
      OBLogP logP;
      
      cout << "logP = " << logP.Predict(mol) << endl;
      \endcode
   */

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

  /** \class OBPSA groupcontrib.h <openbabel/groupcontrib.h>
      \brief calculate the TPSA (topological polar surface area).
 
      This class uses the JOELib2 group contribution algorithm to calculate 
      the TPSA (Topological Polar Surface Area) of a molecule.

      example:
      \code
      #include <openbabel/groupcontrib.h>
      #include <openbabel/mol.h>

      OBMol mol;
      OBLogP psa;
      
      cout << "TPSA = " << psa.Predict(mol) << endl;
      \endcode
   */

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
  
  /** \class OBMR groupcontrib.h <openbabel/groupcontrib.h>
      \brief calculate the MR (molar refractivity).
 
      This class uses the JOELib2 group contribution algorithm to calculate 
      the MR (Molar Refractivity) of a molecule.

      example:
      \code
      #include <openbabel/groupcontrib.h>
      #include <openbabel/mol.h>

      OBMol mol;
      OBLogP mr;
      
      cout << "MR = " << mr.Predict(mol) << endl;
      \endcode
   */

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

//! \file groupcontrib.cpp
//! \brief Handle logP, PSA and other group-based prediction algorithms.

/**********************************************************************
groupcontrib.cpp - Handle logP, PSA, MR, and other group-based predictions

Copyright (C) 2007      by Tim Vandermeersch
              2001-2007 by Stephen Jelfs
              2001-2007 by Joerg Kurt Wegner, me@cheminformatics.eu
              2007      by Chris Morley

Original version: JOELib2, http://joelib.sf.net

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
#include <vector>
#include <utility>
#include <cstdlib>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/oberror.h>
#include <openbabel/parsmart.h>
#include <openbabel/bitvec.h>
#include <openbabel/groupcontrib.h>
#include <openbabel/locale.h>
#include <openbabel/elements.h>

using namespace std;

namespace OpenBabel
{

  const char* OBGroupContrib::Description()
  {
   //Adds name of datafile containing SMARTS strings to the description
    static string txt;
    txt =  _descr;
    txt += "\n Datafile: ";
    txt += _filename;
    txt += "\nOBGroupContrib is definable";
    return txt.c_str();
  }

  bool OBGroupContrib::ParseFile()
  {
    OBSmartsPattern *sp;

    // open data file
    ifstream ifs;

    if (OpenDatafile(ifs, _filename).length() == 0) {
      obErrorLog.ThrowError(__FUNCTION__, " Could not find contribution data file.", obError);
      return false;
    }

    // Set the locale for number parsing to avoid locale issues: PR#1785463
    obLocale.SetLocale();

    vector<string> vs;
    bool heavy = false;
    string ln;
    while(getline(ifs,ln)){
      if(ln[0]=='#') continue;
      if(ln.find(";heavy")!=string::npos)
        heavy=true;
      if(ln.find(";debug")!=string::npos)
        _debug = true;
      if(ln[0]==';') continue;
      tokenize(vs, ln);

      if (vs.size() < 2)
        continue;

      sp = new OBSmartsPattern;//causes non-serious memory leak.
      // Could be cured by copying OBSmartsPattern rather than a pointer in vectors
      if (sp->Init(vs[0]))
      {
        if (heavy)
          _contribsHeavy.push_back(pair<OBSmartsPattern*, double> (sp, atof(vs[1].c_str())));
        else
          _contribsHydrogen.push_back(pair<OBSmartsPattern*, double> (sp, atof(vs[1].c_str())));
      }
      else
      {
        delete sp;
        sp = NULL;
        obErrorLog.ThrowError(__FUNCTION__, " Could not parse SMARTS from contribution data file", obInfo);

        // return the locale to the original one
        obLocale.RestoreLocale();

        return false;
      }
    }

    // return the locale to the original one
    obLocale.RestoreLocale();
    return true;
  }


  double OBGroupContrib::Predict(OBBase* pOb, string* param)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(!pmol)
      return 0.0;

    //Need to add hydrogens, so do this to a copy to leave original unchanged
    OBMol mol(*pmol);
    mol.AddHydrogens(false, false);

    //Read in data, unless it has already been done.
    if(_contribsHeavy.empty() && _contribsHydrogen.empty())
      ParseFile();

    vector<vector<int> > _mlist; // match list for atom typing
    vector<vector<int> >::iterator j;
    vector<pair<OBSmartsPattern*, double> >::iterator i;

    stringstream debugMessage;
    OBBitVec seenHeavy(mol.NumAtoms() + 1);
    OBBitVec seenHydrogen(mol.NumAtoms() + 1);
    vector<double> atomValues(mol.NumAtoms(), 0.0);

    OBMol tmpmol;
    tmpmol = mol;
    tmpmol.ConvertDativeBonds();

    // atom contributions
    if (_debug) debugMessage << "Heavy atom contributions:" << endl;
    for (i = _contribsHeavy.begin();i != _contribsHeavy.end();++i) {
      if (i->first->Match(tmpmol)) {
        _mlist = i->first->GetMapList();
        for (j = _mlist.begin();j != _mlist.end();++j) {
          atomValues[(*j)[0] - 1] = i->second;
          seenHeavy.SetBitOn((*j)[0]);
	        if (_debug)
            debugMessage << (*j)[0] << " = " << i->first->GetSMARTS() << " : " << i->second << endl;
        }
      }
    }

    vector<double> hydrogenValues(tmpmol.NumAtoms(), 0.0);

    // Hydrogen contributions - note that matches to hydrogens themselves are ignored
    if (_debug) debugMessage << "  Hydrogen contributions:" << endl;
    for (i = _contribsHydrogen.begin();i != _contribsHydrogen.end();++i) {
      if (i->first->Match(tmpmol)) {
        _mlist = i->first->GetMapList();
        for (j = _mlist.begin();j != _mlist.end();++j) {
          if (tmpmol.GetAtom((*j)[0])->GetAtomicNum() == OBElements::Hydrogen)
            continue;
          int Hcount = tmpmol.GetAtom((*j)[0])->GetExplicitDegree() - tmpmol.GetAtom((*j)[0])->GetHvyDegree();
          hydrogenValues[(*j)[0] - 1] = i->second * Hcount;
          seenHydrogen.SetBitOn((*j)[0]);
          if (_debug)
            debugMessage << (*j)[0] << " = " << i->first->GetSMARTS() << " : " << i->second << " Hcount " << Hcount << endl;
        }
      }
    }

    // total atomic and hydrogen contribution
    double total = 0.0;

    if (_debug)
      debugMessage << "  Final contributions:\n";
    for (unsigned int index = 0; index < tmpmol.NumAtoms(); index++) {
      if (tmpmol.GetAtom(index + 1)->GetAtomicNum() == OBElements::Hydrogen)
        continue;
      total += atomValues[index];
      total += hydrogenValues[index];
      if (_debug) {
        debugMessage << index+1 << " = " << atomValues[index] << " ";
        if (!seenHeavy.BitIsSet(index + 1)) debugMessage << "un";
        debugMessage << "matched...";
        int Hcount = tmpmol.GetAtom(index + 1)->GetExplicitDegree() - tmpmol.GetAtom(index + 1)->GetHvyDegree();
        debugMessage << "   " << Hcount << " hydrogens = " << hydrogenValues[index] << " ";
        if (!seenHydrogen.BitIsSet(index + 1)) debugMessage << "un";
        debugMessage << "matched\n";
      }
    }

    if (_debug)
      obErrorLog.ThrowError(__FUNCTION__, debugMessage.str(), obWarning);

    return total;
  }

  //******************************************************
  // Make global instances for descriptors which are all calculated
  // from group contibutions in the same way but with different data.

  // LogP (octanol/water partition coefficient)
  OBGroupContrib thelogP("logP", "logp.txt",
    "octanol/water partition coefficient");

  // TPSA (topological polar surface area)
  OBGroupContrib theTPSA("TPSA", "psa.txt",
    "topological polar surface area");

  // MR (molar refractivity)
  OBGroupContrib theMR("MR", "mr.txt",
    "molar refractivity");

//! \file groupcontrib.cpp
//! \brief Handle logP, PSA and other group-based prediction algorithms.

  }//namespace

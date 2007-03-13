/**********************************************************************
psa.cpp - Handle polar surface area (PSA) prediction algorithms.
 
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

#include <openbabel/psa.h>

using namespace std;

namespace OpenBabel
{
  OBPSA::OBPSA()
  {
    OBSmartsPattern *sp;
    
    // open data/psa.txt
    ifstream ifs;

    if (OpenDatafile(ifs, "psa.txt").length() == 0) {
      obErrorLog.ThrowError(__FUNCTION__, " Could not find PSA data file.", obError);
      return;
    }


    vector<string> vs;
    
    char buffer[80];
    while (ifs.getline(buffer, 80)) {
      if (EQn(buffer, "#", 1)) continue;
      if (EQn(buffer, ";", 1)) continue;

	
      tokenize(vs, buffer);
      if (vs.size() < 2)
        continue;
      
      sp = new OBSmartsPattern;
      if (sp->Init(vs[0])) {
        _PSAcontribs.push_back(pair<OBSmartsPattern*, double> (sp, atof(vs[1].c_str())));
      } else {
        delete sp;
        sp = NULL;
        obErrorLog.ThrowError(__FUNCTION__, " Could not parse SMARTS from psa.txt", obInfo);
      }
    }
  }
  
  OBPSA::~OBPSA()
  {
  }
  
  double OBPSA::GroupContributions(OBMol &mol)
  {
    vector<vector<int> > _mlist; //!< match list for atom typing
    vector<vector<int> >::iterator j;
    vector<pair<OBSmartsPattern*, double> >::iterator i;
    
    vector<double> atomValues(mol.NumAtoms(), 0.0);

    // atom contributions
    for (i = _PSAcontribs.begin();i != _PSAcontribs.end();++i) {
      if (i->first->Match(mol)) {
        _mlist = i->first->GetMapList();
        for (j = _mlist.begin();j != _mlist.end();++j) {
	  atomValues[(*j)[0] - 1] = i->second;
        }
      }
    }
    
    // total atomic and hydrogen contribution
    double total = 0.0;

    for (int index = 0; index < mol.NumAtoms(); index++) {
      if (mol.GetAtom(index+1)->IsHydrogen())
        continue;

      total += atomValues[index];
    }

    return total;
  }
 

} // end namespace OpenBabel

//! \file psa.cpp
//! \brief Handle PSA prediction algorithms.

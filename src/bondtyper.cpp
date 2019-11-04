/**********************************************************************
bondtyper.cpp - Bond typer to perceive connectivity and bond orders/types.

Copyright (C) 2003-2006 by Geoffrey R. Hutchison

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
#include <cstdlib>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/oberror.h>
#include <openbabel/bondtyper.h>
#include <openbabel/elements.h>

// data header with default parameters
#include "bondtyp.h"

using namespace std;

namespace OpenBabel
{
  extern OBMessageHandler obErrorLog;


  //! Global OBBondTyper for perception of bond order assignment.
#if __cplusplus >= 201103L
  thread_local //this is required for correct multi-threading
#endif
	OBBondTyper  bondtyper;

  /*! \class OBBondTyper bondtyper.h <openbabel/bondtyper.cpp>
    \brief Assigns bond types for file formats without bond information

    The OBBondTyper class is designed to read in a list of bond typing
    rules and apply them to molecules. It is called from the
    OBMol::PerceiveBondOrders() method.
  */
  OBBondTyper::OBBondTyper()
  {
    _init = false;
    _dir = BABEL_DATADIR;
    _envvar = "BABEL_DATADIR";
    _filename = "bondtyp.txt";
    _subdir = "data";
    _dataptr = BondTypeData;
  }

  void OBBondTyper::ParseLine(const char *buffer)
  {
    vector<string> vs;
    vector<int>    bovector;
    OBSmartsPattern *sp;

    if (buffer[0] != '#')
      {
        tokenize(vs,buffer);
        // Make sure we actually have a SMARTS pattern plus at least one triple
        // and make sure we have the correct number of integers
        if (vs.size() < 4)
          return; // just ignore empty (or short lines)
        else if (vs.size() >= 4 && (vs.size() % 3 != 1))
          {
            stringstream errorMsg;
            errorMsg << " Error in OBBondTyper. Pattern is incorrect, found "
                     << vs.size() << " tokens." << endl;
            errorMsg << " Buffer is: " << buffer << endl;
            obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
            return;
          }

        sp = new OBSmartsPattern;
        if (sp->Init(vs[0]))
          {
            for (unsigned int i = 1; i < vs.size() ; ++i)
              {
                bovector.push_back( atoi((char *)vs[i].c_str()) );
              }
            _fgbonds.push_back(pair<OBSmartsPattern*,vector<int> >
                               (sp, bovector));
          }
        else
          {
            delete sp;
            sp = NULL;
          }
      }
  }

  OBBondTyper::~OBBondTyper()
  {
    vector<pair<OBSmartsPattern*, vector<int> > >::iterator i;
    for (i = _fgbonds.begin();i != _fgbonds.end();++i)
      {
        delete i->first;
        i->first = NULL;
      }
  }

  void OBBondTyper::AssignFunctionalGroupBonds(OBMol &mol)
  {
    if (!_init)
      Init();

    OBSmartsPattern *currentPattern;
    OBBond *b1, *b2;
    OBAtom *a1,*a2, *a3;
    double angle, dist1, dist2;
    vector<int> assignments;
    vector<vector<int> > mlist;
    vector<vector<int> >::iterator matches, l;
    vector<pair<OBSmartsPattern*, vector<int> > >::iterator i;
    unsigned int j;

    // Loop through for all the functional groups and assign bond orders
    for (i = _fgbonds.begin();i != _fgbonds.end();++i)
      {
        currentPattern = i->first;
        assignments = i->second;

        if (currentPattern && currentPattern->Match(mol))
          {
            mlist = currentPattern->GetUMapList();
            for (matches = mlist.begin(); matches != mlist.end(); ++matches)
              {
                // Now loop through the bonds to assign from _fgbonds
                for (j = 0; j < assignments.size(); j += 3)
                  {
                    // along the assignments vector: atomID1 atomID2 bondOrder
                    a1 = mol.GetAtom((*matches)[ assignments[j] ]);
                    a2 = mol.GetAtom((*matches)[ assignments[j+1 ] ]);
                    if (!a1 || !a2) continue;

                    b1 = a1->GetBond(a2);

                    if (!b1) continue;
                    b1->SetBondOrder(assignments[j+2]);
                  } // bond order assignments
              } // each match
          } // current pattern matches

      } // for(functional groups)

    // FG with distance and/or bond criteria
    // Carbonyl oxygen C=O (O must be neutral)
    OBSmartsPattern carbo; carbo.Init("[#8D1;!-][#6](*)(*)");

    if (carbo.Match(mol))
      {
        mlist = carbo.GetUMapList();
        for (l = mlist.begin(); l != mlist.end(); ++l)
          {
            a1 = mol.GetAtom((*l)[0]);
            a2 = mol.GetAtom((*l)[1]);

            angle = a2->AverageBondAngle();
            dist1 = a1->GetDistance(a2);

            // carbonyl geometries ?
            if (angle > 115 && angle < 150 && dist1 < 1.28) {

              if ( !a1->HasDoubleBond() ) {// no double bond already assigned
                b1 = a1->GetBond(a2);

                if (!b1 ) continue;
                b1->SetBondOrder(2);
              }
            }
          }
      } // Carbonyl oxygen

    // thione C=S
    OBSmartsPattern thione; thione.Init("[#16D1][#6](*)(*)");

    if (thione.Match(mol))
      {
        mlist = thione.GetUMapList();
        for (l = mlist.begin(); l != mlist.end(); ++l)
          {
            a1 = mol.GetAtom((*l)[0]);
            a2 = mol.GetAtom((*l)[1]);

            angle = a2->AverageBondAngle();
            dist1 = a1->GetDistance(a2);

            // thione geometries ?
            if (angle > 115 && angle < 150 && dist1 < 1.72) {

              if ( !a1->HasDoubleBond() ) {// no double bond already assigned
                b1 = a1->GetBond(a2);

                if (!b1 ) continue;
                b1->SetBondOrder(2);
              }
            }
          }
      } // thione

    // Isocyanate N=C=O or Isothiocyanate
    bool dist1OK;
    OBSmartsPattern isocyanate; isocyanate.Init("[#8,#16;D1][#6D2][#7D2]");
    if (isocyanate.Match(mol))
      {
        mlist = isocyanate.GetUMapList();
        for (l = mlist.begin(); l != mlist.end(); ++l)
          {
            a1 = mol.GetAtom((*l)[0]);
            a2 = mol.GetAtom((*l)[1]);
            a3 = mol.GetAtom((*l)[2]);

            angle = a2->AverageBondAngle();
            dist1 = a1->GetDistance(a2);
            dist2 = a2->GetDistance(a3);

            // isocyanate geometry or Isothiocyanate geometry ?
            if (a1->GetAtomicNum() == OBElements::Oxygen)
              dist1OK =  dist1 < 1.28;
            else
              dist1OK =  dist1 < 1.72;

            if (angle > 150 && dist1OK && dist2 < 1.34) {

              b1 = a1->GetBond(a2);
              b2 = a2->GetBond(a3);
              if (!b1 || !b2) continue;
              b1->SetBondOrder(2);
              b2->SetBondOrder(2);

            }

          }
      } // Isocyanate

    // oxime C=S
    OBSmartsPattern oxime; oxime.Init("[#6D3][#7D2][#8D2]");

    if (oxime.Match(mol))
      {
        mlist = oxime.GetUMapList();
        for (l = mlist.begin(); l != mlist.end(); ++l)
          {
            a1 = mol.GetAtom((*l)[0]);
            a2 = mol.GetAtom((*l)[1]);

            angle = a2->AverageBondAngle();
            dist1 = a1->GetDistance(a2);

            // thione geometries ?
            if (angle > 110 && angle < 150 && dist1 < 1.4) {

              if ( !a1->HasDoubleBond() ) {// no double bond already assigned
                b1 = a1->GetBond(a2);

                if (!b1 ) continue;
                b1->SetBondOrder(2);
              }
            }
          }
      } // oxime

    // oxido-n+ (e.g., pyridine-N-oxide)
    OBSmartsPattern oxidopyr; oxidopyr.Init("[#8D1][#7D3r6]");
    if (oxidopyr.Match(mol))
      {
        mlist = oxidopyr.GetUMapList();
        for (l = mlist.begin(); l != mlist.end(); ++l)
          {
            a1 = mol.GetAtom((*l)[0]);
            a2 = mol.GetAtom((*l)[1]);

            angle = a2->AverageBondAngle();
            dist1 = a1->GetDistance(a2);

            if (angle > 110 && angle < 150 && dist1 < 1.35) {
              a1->SetFormalCharge(-1); // oxygen
              a2->SetFormalCharge(+1); // nitrogen
            }
          }
      }

  }

} //namespace OpenBabel;

//! \file bondtyper.cpp
//! \brief Bond typer to perceive connectivity and bond orders/types.
//! \todo Needs to add aromatic ring bond order assignment.
//!   Eventually need to migrate OBMol::PerceiveBondOrders(),
//!   OBMol::ConnectTheDots(), and possibly some of the Kekulize routines

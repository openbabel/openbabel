/**********************************************************************
patty.cpp - Programmable atom typer.

Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
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
#include <openbabel/oberror.h>
#include <openbabel/patty.h>

#include <cstring>
#include <cstdlib>

// Simple programmable atom typer
// WPW - 070199
// Usage is in sample main below

using namespace std;

namespace OpenBabel
{

  /*! \class patty patty.h <openbabel/patty.h>
    \brief Programmable Atom Typer

    \deprecated This code is currently not used by the Open Babel
    library. Instead, OBAtomTyper and OBAromaticTyper are used. Unless
    there is interest in retaining this independent class, it will be
    removed in the future.

    Patty stands for programmable atom typer. The patty class was kindly
    donated by W. Patrick Walters. The patty class provides a more
    flexible means for atom typing than the OBAtomTyper. The behavior of
    patty is similar to the OBAtomTyper in that rules apply only to the
    first atom in the SMARTS pattern. The patty class can read any free
    format ASCII file which contains SMARTS patterns associated with user
    defined atom type. The following is an example of a valid patty rule
    \code
    O=C hbacceptor
    \endcode
    The following is a code sample that demonstrates the use of patty
    class:
    \code
    OBMol mol;

    string rulefile = "rules.txt";
    patty p;
    p.read_rules(p);
    vector<string> type;
    p.assign_types(mol,type);
    for (int i = 1;i <= mol.NumAtoms();++i)
       cout << "atom number " << i << " was given a type " << type[i] << endl;
    \endcode
    The array indices in the vector<string> into which the result values
    are placed match the corresponding atom numbers. Since atoms are
    numbered beginning from one, the first element in the vector<string>
    is empty, and the values are placed in [1...mol.NumAtoms()].
  */

  void patty::read_rules(const string &infile)
  {
    ifstream ifs, ifs1, *ifsP;
    vector<string> vs;
    char buffer[BUFF_SIZE];
    char tmp_str[BUFF_SIZE];
    string patty_dir;
    OBSmartsPattern *sp;

    ifs.open(infile.c_str());
    ifsP= &ifs;
    if (!ifs)
      {
        if (getenv("BABEL_DATADIR") == NULL)
          {
            stringstream errorMsg;
            errorMsg << "The BABEL_DATADIR environment variable is not defined" << endl;
            errorMsg << "Please define it so the program can find " << infile << endl;
            obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
            //            exit(0);
          }
        else
          patty_dir = getenv("BABEL_DATADIR");
        patty_dir += FILE_SEP_CHAR;
        patty_dir += infile;
        ifs1.open(patty_dir.c_str());
        ifsP= &ifs1;
        //     if (!ifs1)
        //    {
        //     cerr << "Could not open " << patty_dir << endl;
        //    exit(0);
        // }
      }

    if (!ifsP)
      {
        stringstream errorMsg;
        errorMsg << "Could not open " << patty_dir << endl;
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
        //        exit(0);
      }
    while (ifsP->getline(buffer,BUFF_SIZE))
      {
        if (buffer[0] != '#')
          {
            tokenize(vs,buffer," \t\n");
            if (vs.size() >= 2)
              {
                strncpy(tmp_str,vs[0].c_str(), sizeof(tmp_str) - 1);
                tmp_str[sizeof(tmp_str) - 1] = '\0';
                sp = new OBSmartsPattern;
                sp->Init(tmp_str);
                _sp.push_back(sp);
                smarts.push_back(vs[0]);
                typ.push_back(vs[1]);
              }
          }
      }
  }

  void patty::assign_rules(std::vector<std::string> &rules)
  {
    vector<string> vs;
    char buffer[BUFF_SIZE];
    char tmp_str[BUFF_SIZE];
    unsigned int i;
    OBSmartsPattern *sp;

    for ( i = 0 ; i < rules.size() ; i++ )
      {
        strncpy(buffer, rules[i].c_str(), BUFF_SIZE - 1); // leave space for null termination
        if (buffer[0] != '#')
          {
            tokenize(vs,buffer," \t\n");
            if (vs.size() >= 2)
              {
                strncpy(tmp_str,vs[0].c_str(), sizeof(tmp_str) - 1);
                tmp_str[sizeof(tmp_str) - 1] = '\0';
                sp = new OBSmartsPattern;
                sp->Init(tmp_str);
                _sp.push_back(sp);
                smarts.push_back(vs[0]);
                typ.push_back(vs[1]);
              }
          }
      }
  }


  void patty::assign_types(OBMol &mol, std::vector<std::string> &atm_typ)
  {
    atm_typ.resize(mol.NumAtoms()+1);

    obErrorLog.ThrowError(__FUNCTION__,
                          "Ran OpenBabel::PATTY::AssignTypes", obAuditMsg);

    for (unsigned int i = 0; i < _sp.size(); ++i)
      {
        _sp[i]->Match(mol);
        vector<vector<int> > match = _sp[i]->GetMapList();
        //vector<vector<int> >& match = _sp[i]->GetMapList();
        if (!match.empty())
          {
            if (debug)
              {
                stringstream errorMsg;
                errorMsg << typ[i] << " " << smarts[i] << " matched ";
                obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obDebug);
              }

            for (unsigned int j = 0; j < match.size(); ++j)
              {
                if (debug)
                  {
                    stringstream errorMsg;
                    errorMsg << match[j][0] << " ";
                    obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obDebug);
                  }
                atm_typ[match[j][0]] = typ[i];
              }
          }
      }
  }

  void patty::assign_types(OBMol &mol,vector<int> &atm_typ)
  {
    atm_typ.resize(mol.NumAtoms()+1);

    obErrorLog.ThrowError(__FUNCTION__,
                          "Ran OpenBabel::PATTY::AssignTypes", obAuditMsg);

    for (unsigned int i = 0; i < _sp.size(); ++i)
      {
        _sp[i]->Match(mol);
        vector<vector<int> > match = _sp[i]->GetMapList();
        //vector<vector<int> >& match = _sp[i]->GetMapList();
        if (!match.empty())
          {
            if (debug)
              {
                stringstream errorMsg;
                errorMsg << typ[i] << " " << smarts[i] << " matched " ;
                obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obDebug);
              }

            for (unsigned int j = 0; j < match.size(); ++j)
              {
                if (debug)
                  {
                    stringstream errorMsg;
                    errorMsg << match[j][0] << " ";
                    obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obDebug);
                  }
                atm_typ[match[j][0]] = type_to_int(typ[i]);
              }
          }
      }
  }


  int patty::type_to_int(const string &type, bool failOnUndefined)
  {
    int result;

    switch(toupper(type.c_str()[0]))
      {
      case 'C' : // CAT - CATION
        result = PT_CATION;
        break;
      case 'A' :
        if (toupper(type.c_str()[1]) == 'N') // ANI - ANION
          result = PT_ANION;
        else
          result = PT_ACCEPTOR;
        break;
      case 'P' : // POL - POLAR
        result = PT_POLAR;
        break;
      case 'D' : // DON - DONOR
        result = PT_DONOR;
        break;
      case 'H' : // HYD - HYDROPHOBIC
        result = PT_HYDROPHOBIC;
        break;
      case 'M' : // Metal
        result = PT_METAL;
        break;
      case 'O' : // OTH - OTHER
        result = PT_OTHER;
        break;
      default :
        // This was added by Brian,
        // Behavior will fail if type is undefined
        if (failOnUndefined)
          {
            stringstream errorMsg;
            errorMsg << "Unable to find type of feature passed in " << endl;
            errorMsg << "Feature passed in is " << type << endl;
            obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
          }
        result = 7;
      }
    return(result);
  }

  //! return null if the type does not exist, the type position otherwise
  //! the first position start at 1
  int patty::Istype(const std::string &type)
  {
    for(unsigned int pos=0; pos < typ.size(); ++pos)
      {
        if(typ[pos] == type)
          return (pos + 1);
      }

    return (0);
  }

}

//! \file patty.cpp
//! \brief Programmable atom typer.

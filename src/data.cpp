/**********************************************************************
data.cpp - Global data and resource file parsers.

Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2008 by Geoffrey R. Hutchison

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

#ifdef WIN32
#pragma warning (disable : 4786)
#endif
#include <cstdlib>
#include <openbabel/babelconfig.h>
#include <openbabel/data.h>
#include <openbabel/data_utilities.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/locale.h>
#include <openbabel/oberror.h>
#include <openbabel/elements.h>

// data headers with default parameters
#include "types.h"
#include "resdata.h"


#if !HAVE_STRNCASECMP
extern "C" int strncasecmp(const char *s1, const char *s2, size_t n);
#endif

using namespace std;

namespace OpenBabel
{
  // Initialize the globals (declared in data.h)
  OBTypeTable ttab;
  OBResidueData resdat;

  OBAtomicHeatOfFormationTable::OBAtomicHeatOfFormationTable(void)
  {
    _init = false;
    _dir = BABEL_DATADIR;
    _envvar = "BABEL_DATADIR";
    _filename = "atomization-energies.txt";
    _subdir = "data";
    Init();
  }

  static double UnitNameToConversionFactor(const char* unit) {
    const char* p = unit;
    switch(p[0]) {
    case 'e':
      if (p[1]=='V' && p[2]=='\0')
        return ELECTRONVOLT_TO_KCALPERMOL; // eV
      if (p[1]=='l' && p[2]=='e' && p[3]=='c' && p[4]=='t' && p[5]=='r' && p[6]=='o' && p[7]=='n' &&
          p[8]=='v' && p[9]=='o' && p[10]=='l' && p[11]=='t' && p[12]=='\0')
        return ELECTRONVOLT_TO_KCALPERMOL; // electronvolt
      break;
    case 'k':
      if (p[1]=='J' && p[2]=='/' && p[3]=='m' && p[4]=='o' && p[5]=='l' && p[6]=='\0')
        return KJPERMOL_TO_KCALPERMOL; // kJ/mol
      if (p[1]=='c' && p[2]=='a' && p[3]=='l' && p[4]=='/' && p[5]=='m' && p[6]=='o' && p[7]=='l' && p[8]=='\0')
        return 1.0; // kcal/mol
      break;
    case 'H':
      if (p[1]=='a' && p[2]=='r' && p[3]=='t' && p[4]=='r' && p[5]=='e' && p[6]=='e' && p[7]=='\0')
        return HARTEE_TO_KCALPERMOL; // Hartree
      break;
    case 'J':
      if (p[1]=='/' && p[2]=='m' && p[3]=='o' && p[4]=='l' && p[5]==' ' && p[6]=='K' && p[7]=='\0')
        return KJPERMOL_TO_KCALPERMOL; // J/mol K
      break;
    case 'R':
      if (p[1]=='y' && p[2]=='d' && p[3]=='b' && p[4]=='e' && p[5]=='r' && p[6]=='g' && p[7]=='\0')
        return RYDBERG_TO_KCALPERMOL; // Rydberg
      break;
    }

    std::stringstream errorMsg;
    errorMsg << "WARNING: Unknown energy unit in thermochemistry file\n";
    obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obWarning);

    return 1.0;
  }

  void OBAtomicHeatOfFormationTable::ParseLine(const char *line)
  {
    char *ptr;
    vector<string> vs;
    OBAtomHOF *oba;

    ptr = const_cast<char*>( strchr(line,'#'));
    if (NULL != ptr)
      ptr[0] = '\0';
    if (strlen(line) > 0)
      {
        tokenize(vs,line,"|");
        if (vs.size() >= 8)
          {
              oba = new OBAtomHOF(vs[0],
                                  atoi(vs[1].c_str()),
                                  vs[2],
                                  vs[3],
                                  atof(vs[4].c_str()),
                                  atof(vs[5].c_str()),
                                  atoi(vs[6].c_str()),
                                  vs[7]);
            _atomhof.push_back(*oba);
          }
      }
  }

  int OBAtomicHeatOfFormationTable::GetHeatOfFormation(std::string elem,
                                                       int charge,
                                                       std::string meth,
                                                       double T,
                                                       double *dhof0,
                                                       double *dhofT,
                                                       double *S0T)
  {
    int    found;
    double Ttol = 0.05; /* Kelvin */
    double Vmodel, Vdhf, S0, HexpT;
    std::vector<OBAtomHOF>::iterator it;
    char desc[128];

    found = 0;
    Vmodel = Vdhf = S0 = HexpT = 0;
    snprintf(desc,sizeof(desc),"%s(0K)",meth.c_str());

    for(it = _atomhof.begin(); it != _atomhof.end(); ++it)
    {
        if ((0 == it->Element().compare(elem)) &&
            (it->Charge() == charge))
        {
            double eFac = UnitNameToConversionFactor(it->Unit().c_str());
            if (fabs(T - it->T()) < Ttol)
            {
                if (0 == it->Method().compare("exp"))
                {
                    if (0 == it->Desc().compare("H(0)-H(T)"))
                    {
                        HexpT += it->Value()*eFac;
                        found++;
                    }
                    else if (0 == it->Desc().compare("S0(T)"))
                    {
                        S0 += it->Value();
                        found++;
                    }
                }
            }
            else if (0 == it->T()) 
            {
                if ((0 == it->Method().compare(meth)) &&
                    (0 == it->Desc().compare(desc)))
                {
                    Vmodel += it->Value()*eFac;
                    found++;
                }
                if (0 == it->Method().compare("exp"))
                {
                    if (0 == it->Desc().compare("DHf(T)"))
                    {
                        Vdhf += it->Value()*eFac;
                        found++;
                    }
                }
            }
        }
    }

    if (found == 4)
    {
        *dhof0 = Vdhf-Vmodel;
        *dhofT = Vdhf-Vmodel-HexpT;
        *S0T   = -S0/4.184;
        return 1;
    }
    return 0;
  }

  /** \class OBTypeTable data.h <openbabel/data.h>
      \brief Atom Type Translation Table

      Molecular file formats frequently store information about atoms in an
      atom type field. Some formats store only the element for each atom,
      while others include hybridization and local environments, such as the
      Sybyl mol2 atom type field. The OBTypeTable class acts as a translation
      table to convert atom types between a number of different molecular
      file formats. The constructor for OBTypeTable automatically reads the
      text file types.txt. An instance of
      OBTypeTable (ttab) is declared external in data.cpp and is referenced as
      extern OBTypeTable ttab in mol.h.  The following code demonstrates how
      to use the OBTypeTable class to translate the internal representation
      of atom types in an OBMol Internal to Sybyl Mol2 atom types.

      \code
      ttab.SetFromType("INT");
      ttab.SetToType("SYB");
      OBAtom *atom;
      vector<OBAtom*>::iterator i;
      string src,dst;
      for (atom = mol.BeginAtom(i);atom;atom = mol.EndAtom(i))
      {
         src = atom->GetType();
         ttab.Translate(dst,src);
         cout << "atom number " << atom->GetIdx() << "has mol2 type " << dst << endl;
      }
      \endcode

      Current atom types include (defined in the top line of the data file types.txt):
      - INT (Open Babel internal codes)
      - ATN (atomic numbers)
      - HYB (hybridization)
      - MMD (MacroModel)
      - MM2 (MM2 force field)
      - XYZ (element symbols from XYZ file format)
      - ALC (Alchemy)
      - HAD (H added)
      - MCML (MacMolecule)
      - C3D (Chem3D)
      - SYB (Sybyl mol2)
      - MOL (Sybyl mol)
      - MAP (Gasteiger partial charge types)
      - DRE (Dreiding)
      - XED (XED format)
      - DOK (Dock)
      - M3D (Molecular Arts M3D)
      - SBN (Sybyl descriptor types for MPD files)
      - PCM (PC Model)
  */

  OBTypeTable::OBTypeTable()
  {
    _init = false;
    _dir = BABEL_DATADIR;
    _envvar = "BABEL_DATADIR";
    _filename = "types.txt";
    _subdir = "data";
    _dataptr = TypesData;
    _linecount = 0;
    _from = _to = -1;
  }

  void OBTypeTable::ParseLine(const char *buffer)
  {
    if (buffer[0] == '#')
      return; // just a comment line

    if (_linecount == 0) {
      tokenize(_colnames,buffer);
      _ncols = _colnames.size();
    }
    else
      {
        vector<string> vc;
        tokenize(vc,buffer);
        if (vc.size() == (unsigned)_ncols)
          _table.push_back(vc);
        else
          {
            stringstream errorMsg;
            errorMsg << " Could not parse line in type translation table types.txt -- incorect number of columns";
            errorMsg << " found " << vc.size() << " expected " << _ncols << ".";
            obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
          }
      }
    _linecount++;
  }

  bool OBTypeTable::SetFromType(const char* from)
  {
    if (!_init)
      Init();

    string tmp = from;

    unsigned int i;
    for (i = 0;i < _colnames.size();++i)
      if (tmp == _colnames[i])
        {
          _from = i;
          return(true);
        }

    obErrorLog.ThrowError(__FUNCTION__, "Requested type column not found", obInfo);

    return(false);
  }

  bool OBTypeTable::SetToType(const char* to)
  {
    if (!_init)
      Init();

    string tmp = to;

    unsigned int i;
    for (i = 0;i < _colnames.size();++i)
      if (tmp == _colnames[i])
        {
          _to = i;
          return(true);
        }

    obErrorLog.ThrowError(__FUNCTION__, "Requested type column not found", obInfo);

    return(false);
  }

  //! Translates atom types (to, from), checking for size of destination
  //!  string and null-terminating as needed
  //! \deprecated Because there is no guarantee on the length of an atom type
  //!  you should consider using std::string instead
  bool OBTypeTable::Translate(char *to, const char *from)
  {
    if (!_init)
      Init();

    bool rval;
    string sto,sfrom;
    sfrom = from;
    rval = Translate(sto,sfrom);
    strncpy(to,(char*)sto.c_str(), OBATOM_TYPE_LEN - 1);
    to[OBATOM_TYPE_LEN - 1] = '\0';

    return(rval);
  }

  bool OBTypeTable::Translate(string &to, const string &from)
  {
    if (!_init)
      Init();

    if (from == "")
      return(false);

    if (_from >= 0 && _to >= 0 &&
        _from < (signed)_table.size() && _to < (signed)_table.size())
      {
        vector<vector<string> >::iterator i;
        for (i = _table.begin();i != _table.end();++i)
          if ((signed)(*i).size() > _from &&  (*i)[_from] == from)
            {
              to = (*i)[_to];
              return(true);
            }
      }

    // Throw an error, copy the string and return false
    obErrorLog.ThrowError(__FUNCTION__, "Cannot perform atom type translation: table cannot find requested types.", obWarning);
    to = from;
    return(false);
  }

  std::string OBTypeTable::Translate(const string &from)
  {
    if (!_init)
      Init();

    if (from.empty())
      return("");

    if (_from >= 0 && _to >= 0 &&
        _from < (signed)_table.size() && _to < (signed)_table.size())
      {
        vector<vector<string> >::iterator i;
        for (i = _table.begin();i != _table.end();++i)
          if ((signed)(*i).size() > _from &&  (*i)[_from] == from)
            {
              return (*i)[_to];
            }
      }

    // Throw an error, copy the string and return false
    obErrorLog.ThrowError(__FUNCTION__, "Cannot perform atom type translation: table cannot find requested types.", obWarning);
    return("");
  }

  std::string OBTypeTable::GetFromType()
  {
    if (!_init)
      Init();

    if (_from > 0 && _from < (signed)_table.size())
      return( _colnames[_from] );
    else
      return( _colnames[0] );
  }

  std::string OBTypeTable::GetToType()
  {
    if (!_init)
      Init();

    if (_to > 0 && _to < (signed)_table.size())
      return( _colnames[_to] );
    else
      return( _colnames[0] );
  }

  void Toupper(string &s)
  {
    unsigned int i;
    for (i = 0;i < s.size();++i)
      s[i] = toupper(s[i]);
  }

  void Tolower(string &s)
  {
    unsigned int i;
    for (i = 0;i < s.size();++i)
      s[i] = tolower(s[i]);
  }

  ///////////////////////////////////////////////////////////////////////
  OBResidueData::OBResidueData()
  {
    _init = false;
    _dir = BABEL_DATADIR;
    _envvar = "BABEL_DATADIR";
    _filename = "resdata.txt";
    _subdir = "data";
    _dataptr = ResidueData;
  }

  bool OBResidueData::AssignBonds(OBMol &mol)
  {
    if (!_init)
      Init();

    OBAtom *a1,*a2;
    OBResidue *r1,*r2;
    vector<OBAtom*>::iterator i,j;
    vector3 v;

    int bo;
    string skipres = ""; // Residue Number to skip
    string rname = "";
    //assign residue bonds
    for (a1 = mol.BeginAtom(i);a1;a1 = mol.NextAtom(i))
      {
        r1 = a1->GetResidue();
        if (r1 == NULL) // atoms may not have residues
          continue;

        if (skipres.length() && r1->GetNumString() == skipres)
          continue;

        if (r1->GetName() != rname)
          {
            skipres = SetResName(r1->GetName()) ? "" : r1->GetNumString();
            rname = r1->GetName();
          }
        //assign bonds for each atom
        for (j=i,a2 = mol.NextAtom(j);a2;a2 = mol.NextAtom(j))
          {
            r2 = a2->GetResidue();
            if (r2 == NULL) // atoms may not have residues
              continue;

            if (r1->GetNumString() != r2->GetNumString())
              break;
            if (r1->GetName() != r2->GetName())
              break;
            if (r1->GetChain() != r2->GetChain())
              break; // Fixes PR#2889763 - Fabian

            if ((bo = LookupBO(r1->GetAtomID(a1),r2->GetAtomID(a2))))
              {
                // Suggested by Liu Zhiguo 2007-08-13
                // for predefined residues, don't perceive connection
                // by distance
                //                v = a1->GetVector() - a2->GetVector();
                //                if (v.length_2() < 3.5) //check by distance
                  mol.AddBond(a1->GetIdx(),a2->GetIdx(),bo);
              }
          }
      }

    int hyb;
    string type;

    //types and hybridization
    rname = ""; // name of current residue
    skipres = ""; // don't skip any residues right now
    for (a1 = mol.BeginAtom(i);a1;a1 = mol.NextAtom(i))
      {
        if (a1->GetAtomicNum() == OBElements::Oxygen && !a1->GetExplicitDegree())
          {
            a1->SetType("O3");
            continue;
          }
        if (a1->GetAtomicNum() == OBElements::Hydrogen)
          {
            a1->SetType("H");
            continue;
          }

        //***valence rule for O-
        if (a1->GetAtomicNum() == OBElements::Oxygen && a1->GetExplicitDegree() == 1)
          {
            OBBond *bond;
            bond = (OBBond*)*(a1->BeginBonds());
            if (bond->GetBondOrder() == 2)
              {
                a1->SetType("O2");
                a1->SetHyb(2);
              }
            else if (bond->GetBondOrder() == 1)
              {
                // Leave the protonation/deprotonation to phmodel.txt
                a1->SetType("O3");
                a1->SetHyb(3);
                // PR#3203039 -- Fix from Magnus Lundborg
                //                a1->SetFormalCharge(0);
              }
            continue;
          }

        r1 = a1->GetResidue();
        if (r1 == NULL) continue; // atoms may not have residues
        if (skipres.length() && r1->GetNumString() == skipres)
          continue;

        if (r1->GetName() != rname)
          {
            // if SetResName fails, skip this residue
            skipres = SetResName(r1->GetName()) ? "" : r1->GetNumString();
            rname = r1->GetName();
          }

        if (LookupType(r1->GetAtomID(a1),type,hyb))
          {
            a1->SetType(type);
            a1->SetHyb(hyb);
          }
        else // try to figure it out by bond order ???
          {}
      }

    return(true);
  }

  void OBResidueData::ParseLine(const char *buffer)
  {
    int bo;
    string s;
    vector<string> vs;

    if (buffer[0] == '#')
      return;

    tokenize(vs,buffer);
    if (!vs.empty())
      {
        if (vs[0] == "BOND")
          {
            s = (vs[1] < vs[2]) ? vs[1] + " " + vs[2] :
              vs[2] + " " + vs[1];
            bo = atoi(vs[3].c_str());
            _vtmp.push_back(pair<string,int> (s,bo));
          }

        if (vs[0] == "ATOM" && vs.size() == 4)
          {
            _vatmtmp.push_back(vs[1]);
            _vatmtmp.push_back(vs[2]);
            _vatmtmp.push_back(vs[3]);
          }

        if (vs[0] == "RES")
          _resname.push_back(vs[1]);

        if (vs[0]== "END")
          {
            _resatoms.push_back(_vatmtmp);
            _resbonds.push_back(_vtmp);
            _vtmp.clear();
            _vatmtmp.clear();
          }
      }
  }

  bool OBResidueData::SetResName(const string &s)
  {
    if (!_init)
      Init();

    unsigned int i;

    for (i = 0;i < _resname.size();++i)
      if (_resname[i] == s)
        {
          _resnum = i;
          return(true);
        }

    _resnum = -1;
    return(false);
  }

  int OBResidueData::LookupBO(const string &s)
  {
    if (_resnum == -1)
      return(0);

    unsigned int i;
    for (i = 0;i < _resbonds[_resnum].size();++i)
      if (_resbonds[_resnum][i].first == s)
        return(_resbonds[_resnum][i].second);

    return(0);
  }

  int OBResidueData::LookupBO(const string &s1, const string &s2)
  {
    if (_resnum == -1)
      return(0);
    string s;

    s = (s1 < s2) ? s1 + " " + s2 : s2 + " " + s1;

    unsigned int i;
    for (i = 0;i < _resbonds[_resnum].size();++i)
      if (_resbonds[_resnum][i].first == s)
        return(_resbonds[_resnum][i].second);

    return(0);
  }

  bool OBResidueData::LookupType(const string &atmid,string &type,int &hyb)
  {
    if (_resnum == -1)
      return(false);

    string s;
    vector<string>::iterator i;

    for (i = _resatoms[_resnum].begin();i != _resatoms[_resnum].end();i+=3)
      if (atmid == *i)
        {
          ++i;
          type = *i;
          ++i;
          hyb = atoi((*i).c_str());
          return(true);
        }

    return(false);
  }

  void OBGlobalDataBase::Init()
  {
    if (_init)
      return;
    _init = true;

    ifstream ifs;
    char charBuffer[BUFF_SIZE];

    // Set the locale for number parsing to avoid locale issues: PR#1785463
    obLocale.SetLocale();

    // Check return value from OpenDatafile
    // Suggestion from Zhiguo Liu
    string fn_open = OpenDatafile(ifs, _filename, _envvar);
    
    // Check _subdir directory
    if (fn_open == "")
      string fn_open = OpenDatafile(ifs, _filename, _subdir);

    if (fn_open != "" && (ifs))
      {
        while(ifs.getline(charBuffer,BUFF_SIZE))
          ParseLine(charBuffer);
      }

    else
      // If all else fails, use the compiled in values
      if (_dataptr)
        {
          obErrorLog.ThrowError(__FUNCTION__, "Cannot open " + _filename + " defaulting to compiled data.", obDebug);

          const char *p1,*p2;
          for (p1 = p2 = _dataptr;*p2 != '\0';++p2)
            if (*p2 == '\n')
              {
                strncpy(charBuffer, p1, (p2 - p1));
                charBuffer[(p2 - p1)] = '\0';
                ParseLine(charBuffer);
                p1 = ++p2;
              }
        }
      else
        {
          string s = "Unable to open data file '";
          s += _filename;
          s += "'";
          obErrorLog.ThrowError(__FUNCTION__, s, obWarning);
        }

    // return the locale to the original one
    obLocale.RestoreLocale();

    if (ifs)
      ifs.close();

    if (GetSize() == 0)
      {
        string s = "Cannot initialize database '";
        s += _filename;
        s += "' which may cause further errors.";
        obErrorLog.ThrowError(__FUNCTION__, s, obWarning);
      }

  }

} // end namespace OpenBabel

//! \file data.cpp
//! \brief Global data and resource file parsers.

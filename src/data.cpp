/**********************************************************************
data.cpp - Global data and resource file parsers.
 
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (c) 2001-2003 by Geoffrey R. Hutchison
 
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

#ifdef WIN32
#pragma warning (disable : 4786)
#endif

#ifndef BUFF_SIZE
#define BUFF_SIZE 1024
#endif

#include "babelconfig.h"
#include "data.h"
#include "element.h"
#include "types.h"
// #include "extable.h"
#include "isotope.h"
#include "mol.h"

using namespace std;

namespace OpenBabel
{

OBElementTable   etab;
OBTypeTable      ttab;
OBIsotopeTable   isotab;

/** \class OBElementTable
    \brief Periodic Table of the Elements
 
    Translating element data is a common task given that many file
  formats give either element symbol or atomic number information, but
  not both. The OBElementTable class facilitates conversion between
  textual and numeric element information. An instance of the
  OBElementTable class (etab) is declared as external in data.cpp. Source
  files that include the header file mol.h automatically have an extern
  definition to etab. The following code sample demonstrates the use
  of the OBElementTable class:
\code
 cout << "The symbol for element 6 is " << etab.GetSymbol(6) << endl;
 cout << "The atomic number for Sulfur is " << etab.GetAtomicNum(16) << endl;
 cout << "The van der Waal radius for Nitrogen is " << etab.GetVdwRad(7);
\endcode
 
  Stored information in the OBElementTable includes atomic:
   - symbols
   - van der Waal radii
   - covalent radii
   - bond order radii (deprecated -- use covalent instead)
   - expected maximum bonding valence
   - molar mass (by IUPAC recommended atomic masses)
   - electronegativity
 
*/

OBElementTable::OBElementTable()
{
    _init = false;
    _dir = BABEL_DATADIR;
    _envvar = "BABEL_DATADIR";
    _filename = "element.txt";
    _subdir = "data";
    _dataptr = ElementData;
}

OBElementTable::~OBElementTable()
{
    vector<OBElement*>::iterator i;
    for (i = _element.begin();i != _element.end();i++)
        delete *i;
}

void OBElementTable::ParseLine(const char *buffer)
{
    int num,maxbonds;
    char symbol[3];
    double Rbo,Rcov,Rvdw,mass, elNeg;

    if (buffer[0] != '#') // skip comment line (at the top)
    {
        // Ignore RGB columns
        sscanf(buffer,"%d %s %lf %lf %lf %d %lf %lf %*lf %*lf %*lf",
               &num,
               symbol,
               &Rcov,
               &Rbo,
               &Rvdw,
               &maxbonds,
               &mass,
               &elNeg);

        OBElement *ele = new OBElement(num,symbol,Rcov,Rbo,Rvdw,maxbonds,mass,elNeg);
        _element.push_back(ele);
    }
}

char *OBElementTable::GetSymbol(int atomicnum)
{
    if (!_init)
        Init();

    if (atomicnum < 0 || atomicnum > static_cast<int>(_element.size()))
        return("\0");

    return(_element[atomicnum]->GetSymbol());
}

int OBElementTable::GetMaxBonds(int atomicnum)
{
    if (!_init)
        Init();

    if (atomicnum < 0 || atomicnum > static_cast<int>(_element.size()))
        return(0);

    return(_element[atomicnum]->GetMaxBonds());
}

double OBElementTable::GetElectroNeg(int atomicnum)
{
    if (!_init)
        Init();

    if (atomicnum < 0 || atomicnum > static_cast<int>(_element.size()))
        return(0.0);

    return(_element[atomicnum]->GetElectroNeg());
}

double OBElementTable::GetVdwRad(int atomicnum)
{
    if (!_init)
        Init();

    if (atomicnum < 0 || atomicnum > static_cast<int>(_element.size()))
        return(0.0);

    return(_element[atomicnum]->GetVdwRad());
}

double OBElementTable::GetBORad(int atomicnum)
{
    if (!_init)
        Init();

    if (atomicnum < 0 || atomicnum > static_cast<int>(_element.size()))
        return(0.0);

    return(_element[atomicnum]->GetBoRad());
}

double OBElementTable::CorrectedBondRad(int atomicnum, int hyb)
{
    double rad;
    if (!_init)
        Init();

    if (atomicnum < 0 || atomicnum > static_cast<int>(_element.size()))
        return(1.0);

    rad = _element[atomicnum]->GetCovalentRad();

    if (hyb == 2)
        rad *= 0.95;
    else if (hyb == 1)
        rad *= 0.90;

    return(rad);
}

double OBElementTable::CorrectedVdwRad(int atomicnum, int hyb)
{
    double rad;
    if (!_init)
        Init();

    if (atomicnum < 0 || atomicnum > static_cast<int>(_element.size()))
        return(1.95);

    rad = _element[atomicnum]->GetVdwRad();

    if (hyb == 2)
        rad *= 0.95;
    else if (hyb == 1)
        rad *= 0.90;

    return(rad);
}

double OBElementTable::GetCovalentRad(int atomicnum)
{
    if (!_init)
        Init();

    if (atomicnum < 0 || atomicnum > static_cast<int>(_element.size()))
        return(0.0);

    return(_element[atomicnum]->GetCovalentRad());
}

double OBElementTable::GetMass(int atomicnum)
{
    if (!_init)
        Init();

    if (atomicnum < 0 || atomicnum > static_cast<int>(_element.size()))
        return(0.0);

    return(_element[atomicnum]->GetMass());
}

int OBElementTable::GetAtomicNum(const char *sym, unsigned short iso)
{
    if (!_init)
        Init();

    vector<OBElement*>::iterator i;
    for (i = _element.begin();i != _element.end();i++)
        if (!strncasecmp(sym,(*i)->GetSymbol(),2))
            return((*i)->GetAtomicNum());
    if (strcasecmp(sym, "D") == 0)
    {
        iso = 2;
        return(1);
    }
    else if (strcasecmp(sym, "T") == 0)
    {
        iso = 3;
        return(1);
    }

    return(0);
}

/** \class OBIsotopeTable
    \brief Table of atomic isotope masses
 
*/

OBIsotopeTable::OBIsotopeTable()
{
    _init = false;
    _dir = BABEL_DATADIR;
    _envvar = "BABEL_DATADIR";
    _filename = "isotope.txt";
    _subdir = "data";
    _dataptr = IsotopeData;
}

void OBIsotopeTable::ParseLine(const char *buffer)
{
    unsigned int atomicNum;
    unsigned int i;
    vector<string> vs;

    pair <unsigned int, double> entry;
    vector <pair <unsigned int, double> > row;

    if (buffer[0] != '#') // skip comment line (at the top)
    {
        tokenize(vs,buffer);
        if (vs.size() > 3) // atomic number, 0, most abundant mass (...)
        {
            atomicNum = atoi(vs[0].c_str());
            for (i = 1; i < vs.size() - 1; i += 2) // make sure i+1 still exists
            {
                entry.first = atoi(vs[i].c_str()); // isotope
                entry.second = atof(vs[i + 1].c_str()); // exact mass
                row.push_back(entry);
            }
            _isotopes.push_back(row);
        }
    }
}

double	OBIsotopeTable::GetExactMass(const unsigned int ele,
                                    const unsigned int isotope)
{
    if (!_init)
        Init();

    if (ele > _isotopes.size())
        return 0.0;

    unsigned int iso;
    for (iso = 0; iso < _isotopes[ele].size(); iso++)
        if (isotope == _isotopes[ele][iso].first)
            return _isotopes[ele][iso].second;

    return 0.0;
}

/** \class OBTypeTable
    \brief Atom Type Translation Table
 
Molecular file formats frequently store information about atoms in an
atom type field. Some formats store only the element for each atom,
while others include hybridization and local environments, such as the
Sybyl mol2 atom type field. The OBTypeTable class acts as a translation
table to convert atom types between a number of different molecular
file formats. The constructor for OBTypeTable automatically reads the
text file types.txt. Just as OBElementTable, an instance of
OBTypeTable (ttab) is declared external in data.cpp and is referenced as
extern OBTypeTable ttab in mol.h.  The following code demonstrates how
to use the OBTypeTable class to translate the internal representation
of atom types in an OBMol Internal to Sybyl Mol2 atom types.
 
\code
ttab.SetFromType("INT");
ttab.SetToType("SYB");
OEAtom *atom;
vector<OEAtom*>::iterator i;
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
- MMD
- MM2 (MM2 force field)
- XYZ (element symbols from XYZ file format)
- ALC (Alchemy file)
- HAD
- MCML
- C3D (Chem3D)
- SYB (Sybyl mol2)
- MOL
- MAP
- DRE
- XED (XED format)
- DOK (Dock)
- M3D
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

    if (_linecount == 0)
        sscanf(buffer,"%d%d",&_ncols,&_nrows);
    else if (_linecount == 1)
        tokenize(_colnames,buffer);
    else
    {
        vector<string> vc;
        tokenize(vc,buffer);
        if (vc.size() == (unsigned)_nrows)
            _table.push_back(vc);
    }
    _linecount++;
}

bool OBTypeTable::SetFromType(char* from)
{
    if (!_init)
        Init();

    string tmp = from;

    unsigned int i;
    for (i = 0;i < _colnames.size();i++)
        if (tmp == _colnames[i])
        {
            _from = i;
            return(true);
        }

    obErrorLog.ThrowError(__FUNCTION__, "Requested type column not found", obInfo);

    return(false);
}

bool OBTypeTable::SetToType(char* to)
{
    if (!_init)
        Init();

    string tmp = to;

    unsigned int i;
    for (i = 0;i < _colnames.size();i++)
        if (tmp == _colnames[i])
        {
            _to = i;
            return(true);
        }

    obErrorLog.ThrowError(__FUNCTION__, "Requested type column not found", obInfo);

    return(false);
}

bool OBTypeTable::Translate(char *to, char *from)
{
    if (!_init)
        Init();

    bool rval;
    string sto,sfrom;
    sfrom = from;
    rval = Translate(sto,sfrom);
    strcpy(to,(char*)sto.c_str());

    return(rval);
}

bool OBTypeTable::Translate(string &to,string &from)
{
    if (!_init)
        Init();

    if (from == "")
        return(false);

    vector<vector<string> >::iterator i;
    for (i = _table.begin();i != _table.end();i++)
        if ((signed)(*i).size() > _from &&  (*i)[_from] == from)
        {
            to = (*i)[_to];
            return(true);
        }

    to = from;
    return(false);
}

void Toupper(string &s)
{
    unsigned int i;
    for (i = 0;i < s.size();i++)
        s[i] = toupper(s[i]);
}

void Tolower(string &s)
{
    unsigned int i;
    for (i = 0;i < s.size();i++)
        s[i] = tolower(s[i]);
}

void OBGlobalDataBase::Init()
{
    if (_init)
        return;
    _init = true;

    char buffer[BUFF_SIZE],subbuffer[BUFF_SIZE];
    ifstream ifs1, ifs2, ifs3, ifs4, *ifsP;
    // First, look for an environment variable
    if (getenv(_envvar.c_str()) != NULL)
    {
        strcpy(buffer,getenv(_envvar.c_str()));
        strcat(buffer,FILE_SEP_CHAR);

        if (!_subdir.empty())
        {
            strcpy(subbuffer,buffer);
            strcat(subbuffer,_subdir.c_str());
            strcat(subbuffer,FILE_SEP_CHAR);
        }

        strcat(buffer,(char*)_filename.c_str());
        strcat(subbuffer,(char*)_filename.c_str());

        ifs1.open(subbuffer);
        ifsP= &ifs1;
        if (!(*ifsP))
        {
            ifs2.open(buffer);
            ifsP = &ifs2;
        }
    }
    // Then, check the configured data directory
    else // if (!(*ifsP))
    {
        strcpy(buffer,_dir.c_str());
        strcat(buffer,FILE_SEP_CHAR);

	strcpy(subbuffer,buffer);
	strcat(subbuffer,BABEL_VERSION);
	strcat(subbuffer,FILE_SEP_CHAR);
        strcat(subbuffer,(char*)_filename.c_str());

        strcat(buffer,(char*)_filename.c_str());

        ifs3.open(subbuffer);
        ifsP= &ifs3;
        if (!(*ifsP))
        {
            ifs4.open(buffer);
            ifsP = &ifs4;
        }
    }

    if ((*ifsP))
    {
        while(ifsP->getline(buffer,BUFF_SIZE))
            ParseLine(buffer);
    }

    else
        // If all else fails, use the compiled in values
        if (_dataptr)
        {
            const char *p1,*p2;
            for (p1 = p2 = _dataptr;*p2 != '\0';p2++)
                if (*p2 == '\n')
                {
                    strncpy(buffer, p1, (p2 - p1));
                    buffer[(p2 - p1)] = '\0';
                    ParseLine(buffer);
                    p1 = ++p2;
                }
        }
        else
        {
            string s = "Unable to open data file '";
            s += _filename;
            s += "'";
	    obErrorLog.ThrowError(__FUNCTION__, s, obInfo);
        }

    if (ifs1)
        ifs1.close();
    if (ifs2)
        ifs2.close();
    if (ifs3)
        ifs3.close();
}

}


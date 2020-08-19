/**********************************************************************
data_utilities.cpp - Global data and resource file parsers.

Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Copyright (C) 2015 by David van der Spoel

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
#include <cstring>
#include <fstream>
#include <iostream>
#include <openbabel/babelconfig.h>
#include <openbabel/data.h>
#include <openbabel/data_utilities.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/generic.h>
#include <openbabel/locale.h>

using std::string;
using std::vector;

namespace OpenBabel {

bool extract_thermochemistry(OpenBabel::OBMol  &mol,
                             bool    bVerbose,
                             int    *Nsymm,
                             int     Nrotbonds,
                             double  dBdT,
                             double *temperature,
                             double *DeltaHf0,
                             double *DeltaHfT,
                             double *DeltaGfT,
                             double *DeltaSfT,
                             double *S0T,
                             double *CVT,
                             double *CPT,
                             std::vector<double> &Scomponents,
                             double *ZPVE)
{
    enum kkTYPE {kkDH, kkDG, kkDS, kkS0, kkCV, kkSt, kkSr, kkSv, kkZP};
    typedef struct {
        std::string term;
        kkTYPE kk;
    } energy_unit;
    double St = 0, Sr = 0, Sv = 0, Sconf = 0, Ssymm = 0;
    double Rgas      = 1.9872041; 
    int    RotSymNum = 1;
    OpenBabel::OBRotationData* rd;
    
    rd = (OpenBabel::OBRotationData*)mol.GetData("RotationData");
    if (nullptr != rd)
    {
        RotSymNum = rd->GetSymmetryNumber();
        if (bVerbose)
        {
            printf("Found symmetry number %d in input file.\n", RotSymNum);
        }
    }
    else if (bVerbose)
    {
        printf("Using default symmetry number %d\n", RotSymNum);
    }
    if ((*Nsymm > 0) && (*Nsymm != RotSymNum))
    {
        // Rgas in cal/mol K http://en.wikipedia.org/wiki/Gas_constant
        Ssymm = -Rgas*log((1.0* *Nsymm)/RotSymNum);
        RotSymNum = *Nsymm;
        if (bVerbose)
        {
            printf("Changing symmetry number to %d\n", RotSymNum);
        }
    }
    else if (*Nsymm == 0)
    {
        *Nsymm = RotSymNum;
    }
    if (Nrotbonds > 0) 
    {
        Sconf = Rgas*Nrotbonds*log(3.0);
    }
    energy_unit eu[] = {
        { "zpe",        kkZP },
        { "DeltaHform", kkDH },
        { "DeltaGform", kkDG },
        { "DeltaSform", kkDS },
        { "S0",         kkS0 },
        { "cv",         kkCV },
        { "Strans",     kkSt },
        { "Srot",       kkSr },
        { "Svib",       kkSv }
    };
#define NEU (sizeof(eu)/sizeof(eu[0]))
    int found = 0;
    std::vector<OpenBabel::OBGenericData*> obdata = mol.GetData();
    for(std::vector<OpenBabel::OBGenericData*>::iterator j = obdata.begin(); (j<obdata.end()); ++j)
    {
        std::string term  = (*j)->GetAttribute();
        double value = atof((*j)->GetValue().c_str());
        double T     = 0;
        {
            size_t lh = term.find("(");
            size_t rh = term.find("K)");
            double TT = atof(term.substr(lh+1,rh-lh-1).c_str());
            if (0 != TT)
            {
                if (0 == T)
                {
                    T            = TT;
                    *temperature = TT;
                }
                else
                {
                    std::cerr << "Different T in the input file, found " << T << " before and now " << TT << ". Output maybe inconsistent." << std::endl;
                    T = TT;
                }
            }
        }
        for(unsigned int i = 0; (i<NEU); i++)
        {
            if (strstr(term.c_str(), eu[i].term.c_str()) != nullptr)
            {
                switch (eu[i].kk)
                {
		            case kkZP:
		                {
		                    *ZPVE = value;
		                }
		                break;
                case kkDH:
                    if (0 == T)
                    {
                        *DeltaHf0 = value;
                    }
                    else
                    {
                        *DeltaHfT = value;
                    }
                    found ++;
                    break;
                case kkDG:
                    *DeltaGfT = value - T*(Ssymm+Sconf)/1000;
                    found ++;
                    break;
                case kkDS:
                    *DeltaSfT = value + Ssymm + Sconf;
                    found ++;
                    break;
                case kkS0:
                    *S0T = value + Ssymm + Sconf;
                    found ++;
                    break;
                case kkSt:
                    St = value;
                    found ++;
                    break;
                case kkSr:
                    Sr = value;
                    found ++;
                    break;
                case kkSv:
                    Sv = value;
                    found ++;
                    break;
                case kkCV:
                    *CVT = value;
                    found++;
                    break;
                default:
                    break;
                }
            }
        }
    }
    double P   = 16.605/4.184; // Convert pressure to kcal/mol
    *CPT       = *CVT + Rgas + (2*P*dBdT + pow(P*dBdT, 2.0)/Rgas);

    Scomponents.push_back(St);
    Scomponents.push_back(Sr);
    Scomponents.push_back(Sv);
    Scomponents.push_back(Ssymm);
    Scomponents.push_back(Sconf);
    if (bVerbose && (Ssymm != 0))
    {
        printf("Applyied symmetry correction to free energy of %g kcal/mol\n",
               -(*temperature*Ssymm)/1000);
    }
    if (bVerbose && (Sconf != 0))
    {
        printf("Applyied conformational correction to free energy of %g kcal/mol\n",
               -(*temperature*Sconf)/1000);
    }
    return (found == 9);
}

// class OBTranslator
OBTranslator::OBTranslator() {
	if (!ttab._init)
		ttab.Init();
	_from = _to = -1;
}

OBTranslator::OBTranslator(const char* from, const char* to) {
	OBTranslator();
	SetFromType(from);
	SetToType(to);
}

bool OBTranslator::SetFromType(const char* from)
{
	string tmp = from;

	unsigned int i;
	for (i = 0;i < ttab._colnames.size();++i)
		if (tmp == ttab._colnames[i])
		{
			_from = i;
			return(true);
		}

	obErrorLog.ThrowError(__FUNCTION__, "Requested type column not found", obInfo);

	return(false);
}

bool OBTranslator::SetToType(const char* to)
{
	string tmp = to;

	unsigned int i;
	for (i = 0;i < ttab._colnames.size();++i)
		if (tmp == ttab._colnames[i])
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
bool OBTranslator::Translate(char *to, const char *from) const
{
	bool rval;
	string sto,sfrom;
	sfrom = from;
	rval = Translate(sto,sfrom);
	strncpy(to,(char*)sto.c_str(), OBATOM_TYPE_LEN - 1);
	to[OBATOM_TYPE_LEN - 1] = '\0';

	return(rval);
}

bool OBTranslator::Translate(string &to, const string &from) const
{
	using std::vector;

	if (from == "")
		return(false);

	if (_from >= 0 && _to >= 0 &&
		_from < (signed)ttab._table.size() && _to < (signed)ttab._table.size())
	{
		vector<vector<string> >::iterator i;
		for (i = ttab._table.begin();i != ttab._table.end();++i)
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

std::string OBTranslator::Translate(const string &from) const
{
	using std::vector;

	if (from.empty())
		return("");

	if (_from >= 0 && _to >= 0 &&
		_from < (signed)ttab._table.size() && _to < (signed)ttab._table.size())
	{
		vector<vector<string> >::iterator i;
		for (i = ttab._table.begin();i != ttab._table.end();++i)
			if ((signed)(*i).size() > _from &&  (*i)[_from] == from)
			{
				return (*i)[_to];
			}
	}

	// Throw an error, copy the string and return false
	obErrorLog.ThrowError(__FUNCTION__, "Cannot perform atom type translation: table cannot find requested types.", obWarning);
	return("");
}

std::string OBTranslator::GetFromType() const
{
	if (_from > 0 && _from < (signed)ttab._table.size())
		return( ttab._colnames[_from] );
	else
		return( ttab._colnames[0] );
}

std::string OBTranslator::GetToType() const
{
	if (_to > 0 && _to < (signed)ttab._table.size())
		return( ttab._colnames[_to] );
	else
		return( ttab._colnames[0] );
}
// End class OBTranslator

// class OBResidueObserver
bool OBResidueObserver::SetResName(const string &s)
{
	if (!resdat._init)
		resdat.Init();

	unsigned int i;

	for (i = 0;i < resdat._resname.size();++i)
		if (resdat._resname[i] == s)
		{
			_resnum = i;
			return(true);
		}

	_resnum = -1;
	return(false);
}

int OBResidueObserver::LookupBO(const string &s)
{
	if (_resnum == -1)
		return(0);

	unsigned int i;
	for (i = 0;i < resdat._resbonds[_resnum].size();++i)
		if (resdat._resbonds[_resnum][i].first == s)
			return(resdat._resbonds[_resnum][i].second);

	return(0);
}

int OBResidueObserver::LookupBO(const string &s1, const string &s2)
{
	if (_resnum == -1)
		return(0);
	string s;

	s = (s1 < s2) ? s1 + " " + s2 : s2 + " " + s1;

	unsigned int i;
	for (i = 0;i < resdat._resbonds[_resnum].size();++i)
		if (resdat._resbonds[_resnum][i].first == s)
			return(resdat._resbonds[_resnum][i].second);

	return(0);
}

bool OBResidueObserver::LookupType(const string &atmid,string &type,int &hyb)
{
	if (_resnum == -1)
		return(false);

	string s;
	vector<string>::iterator i;

	for (i = resdat._resatoms[_resnum].begin();i != resdat._resatoms[_resnum].end();i+=3)
		if (atmid == *i)
		{
			++i;
			type = *i;
			++i;
			hyb = atoi((*i).c_str());
			return(true);
		}

	return(false);
} // End class OBResidueObserver
}

//! \file data_utilities.cpp
//! \brief Data related tools and utilities

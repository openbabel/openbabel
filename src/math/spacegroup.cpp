/**********************************************************************
spacegroup.cpp - Handle Space Groups.

Copyright (C) 2007-2011 by Jean Bréfort

This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License as
published by the Free Software Foundation; either version 2 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/
#include <openbabel/babelconfig.h>

#include <openbabel/math/spacegroup.h>
#include <openbabel/data.h>
#include <openbabel/obutil.h>
#include <iostream>
#include <map>
#include <set>
#include <vector>
#include <locale>

#include <cstdarg>
#include <cstdlib>

using namespace std;

namespace OpenBabel
{
  /** Function to remove whitespaces from a string, returning
   * a new string
   */
  std::string RemoveWhiteSpaceUnderscore(const string &in){
    std::string out=in;
    for(std::string::iterator pos=out.begin();pos!=out.end();){
      if( ((char)(*pos)==' ') ||((char)(*pos)=='_'))  pos=out.erase(pos);
      else ++pos;
    }
    return out;
  }

  class SpaceGroups: public OBGlobalDataBase
  {
  public:
    SpaceGroups();
    virtual ~SpaceGroups();

		void ParseLine(const char*);
    size_t GetSize() { return sgs.size();}
    bool Inited() { return _init;}

    map<string, const SpaceGroup*> sgbn;
    vector< list<const SpaceGroup*> > sgbi;
    set<SpaceGroup*> sgs;
  };

  static SpaceGroups _SpaceGroups;

  SpaceGroups::SpaceGroups()
  {
    sgbi.assign(230, list<const SpaceGroup*>());
    _dir = BABEL_DATADIR;
    _envvar = "BABEL_DATADIR";
    _filename = "space-groups.txt";
    _subdir = "data";
  }

  SpaceGroups::~SpaceGroups()
  {
    set<SpaceGroup*>::iterator i, end = sgs.end();
    for (i = sgs.begin(); i != end; ++i)
      delete (*i);
  }

  enum
    {
      SPACE_GROUP_ID,
      SPACE_GROUP_HALL,
      SPACE_GROUP_HM,
      SPACE_GROUP_TRANSFORM
    };

  void SpaceGroups::ParseLine(const char* line)
  {
    static SpaceGroup *group = nullptr;
    static int step = SPACE_GROUP_ID;
    static string HMs;
    switch (step)
      {
      case SPACE_GROUP_ID:
        group = new SpaceGroup();
        group->SetId(atoi (line));
        step++;
        break;
      case SPACE_GROUP_HALL:
        group->SetHallName(line);
        step++;
        break;
      case SPACE_GROUP_HM:
        {
          string linestr = std::string(line);
          std::string::size_type idx = linestr.find(',');
          if (idx != std::string::npos)
            {
              std::string alt = linestr.substr(0, idx);
              if (alt.length() > 0 && _SpaceGroups.sgbn[alt] == nullptr)
                _SpaceGroups.sgbn[alt] = group;
              std::string stripped_HM=RemoveWhiteSpaceUnderscore(alt);
              if (stripped_HM.length() > 0 && _SpaceGroups.sgbn[stripped_HM] == nullptr)
                _SpaceGroups.sgbn[stripped_HM] = group;
              group->SetHMName(linestr.substr(idx+1, std::string::npos).c_str());
            }
          else
            group->SetHMName(line);
          step++;
          break;
        }
      case SPACE_GROUP_TRANSFORM:
        if (strlen(line) == 0)
          {
            step = SPACE_GROUP_ID;
            if (HMs.length() > 0)
              group->RegisterSpaceGroup(1, HMs.c_str());
            else
              group->RegisterSpaceGroup();
            group = nullptr;
            HMs.clear();
		      }
        else
          group->AddTransform(line);
        break;
      }
  }

  SpaceGroup::SpaceGroup():
    HEXAGONAL_ORIGIN(10), m_HM(""),m_Hall(""),m_id(0),m_OriginAlternative(0)
  {
  }

  SpaceGroup::~SpaceGroup()
  {
    list<transform3d*>::iterator i, end = m_transforms.end();
    for (i = m_transforms.begin(); i != end; ++i)
      delete *i;
  }

  void SpaceGroup::SetHMName(const char *name_in)
  {
    string name = std::string(name_in);
    std::string::size_type idx = name.find(':');
    if (idx != std::string::npos)
      {
        std::string origin = name.substr(idx + 1, std::string::npos);
        if (origin == "H")
          {
            m_OriginAlternative = HEXAGONAL_ORIGIN;
          } else {
            m_OriginAlternative = atoi (origin.c_str());
          }
      }
    m_HM = name;
  }

  /*!
   */
  void SpaceGroup::AddTransform(const string &s)
  {
    matrix3x3 m;
    vector3 v;
    locale cLocale("C");

    if (s.find(',') != string::npos)
      {
        string s1 = RemoveWhiteSpaceUnderscore(s);
        istringstream iss(s1);
        iss.imbue(cLocale);

        string row;
        int i;
        size_t j;
        bool neg;
        double *t;
        for (i = 0; i < 3; i++)
          {
            getline(iss, row, ',');
            j = 0;
            neg = false;
            while (j < row.length())
              {
                switch (row[j])
                  {
                  case '0':
                  case '.': // anticipating something like 0.5 or .3333
                    {
                      char *end;
                      switch (i)
                        {
                        case 0:
                          t = &v.x();
                          break;
                        case 1:
                          t = &v.y();
                          break;
                        case 2:
                          t = &v.z();
                          break;
                        }
                      *t = strtod(row.c_str() + j, &end);
                      j = end - row.c_str() - 1;
                      if (neg)
                        *t = - *t;
                      break;
                    }
                  case '1':
                  case '2':
                  case '3':
                  case '4':
                  case '5':
                  case '6':
                  case '7':
                  case '8':
                  case '9':
                    if (j+2 < row.length() && row[j+1] == '/')
                      {
                        double *t = nullptr;
                        switch (i)
                          {
                          case 0:
                            t = &v.x();
                            break;
                          case 1:
                            t = &v.y();
                            break;
                          case 2:
                            t = &v.z();
                            break;
                          }
                        *t = ((double) (row[j] - '0')) / (row[j+2] - '0');
                        if (neg)
                          *t = - *t;

                        j +=2;
                      }
                    break;
                  case '-':
                    neg = true;
                    break;
                  case '+':
                    neg = false;
                    break;
                  case 'X':
                  case 'x':
                    m(i, 0) = (neg)? -1.: 1.;
                  break;
                  case 'Y':
                  case 'y':
                    m(i, 1) = (neg)? -1.: 1.;
                  break;
                  case 'Z':
                  case 'z':
                    m(i, 2) = (neg)? -1.: 1.;
                  break;
                  }
                j++;
              }
          }
      }
    else if (s.find(' ') != string::npos)
      {
        istringstream iss(s);
        iss.imbue(cLocale);
        /* supposing the string is a list of at least 12 float values. If there are
           16, the last four are 0., 0., 0. and 1. and are not needed */
        iss >> m(0,0) >> m(0,1) >> m(0,2) >> v.x();
        iss >> m(1,0) >> m(1,1) >> m(1,2) >> v.y();
        iss >> m(2,0) >> m(2,1) >> m(2,2) >> v.z();
      }
		if (v.x() < 0)
			v.x() += 1.;
		else if (v.x() >= 1.)
			v.x() -= 1.;
		if (v.y() < 0)
			v.y() += 1.;
		else if (v.y() >= 1.)
			v.y() -= 1.;
		if (v.z() < 0)
			v.z() += 1.;
		else if (v.z() >= 1.)
			v.z() -= 1.;

    // only push_back unique transformations
    transform3dIterator i, iend = m_transforms.end();
    transform3d* candidate = new transform3d (m, v);
    bool transform_exists = false;

    for (i = m_transforms.begin(); i!= iend; i++)
    {
      if (candidate->DescribeAsString() == (*i)->DescribeAsString())
      {
        transform_exists = true;
        break;
      }
    }

    if (!transform_exists)
      m_transforms.push_back (candidate);
  }

  /*!
   */
  list<vector3> SpaceGroup::Transform(const vector3 &v) const
  {
    static double prec = 2e-5;
    list<vector3> res;
    transform3dIterator i, iend = m_transforms.end();
    for (i = m_transforms.begin(); i!= iend; i++)
      {
        vector3 t;
        t = *(*i) * v;
        if (t.x() < 0.)
          t.x() += 1.;
        if (t.x() >= 1.)
          t.x() -= 1.;
        if (t.y() < 0.)
          t.y() += 1.;
        if (t.y() >= 1.)
          t.y() -= 1.;
        if (t.z() < 0.)
          t.z() += 1.;
        if (t.z() >= 1.)
          t.z() -= 1.;
        list<vector3>::iterator j, jend = res.end();
        bool duplicate = false;
        for (j = res.begin(); j != jend; ++j)
          if (fabs(t.x() - (*j).x()) < prec &&
              fabs(t.y() - (*j).y()) < prec &&
              fabs(t.z() - (*j).z()) < prec)
            {
              duplicate = true;
              break;
            }
        if (!duplicate)
          res.push_back (t);
      }
    return res;
  }

  /*!
   */
  transform3d const * SpaceGroup::BeginTransform(transform3dIterator &i) const
  {
    i = m_transforms.begin ();
    return (i == m_transforms.end())? static_cast<transform3d*>(nullptr): *i++;
  }

  /*!
   */
  transform3d const * SpaceGroup::NextTransform(transform3dIterator &i) const
  {
    return (i == m_transforms.end())? static_cast<transform3d*>(nullptr): *i++;
  }

  /*!
   */
  const SpaceGroup * SpaceGroup::GetSpaceGroup (char const *name)
  {
    return GetSpaceGroup(std::string(name)); // let's only use one method
  }

  /*!
   */
  const SpaceGroup * SpaceGroup::GetSpaceGroup (const string &name_in)
  {
    if (!_SpaceGroups.Inited())
      _SpaceGroups.Init();

    // This needs to be more forgiving
    // First, try it without removing the white space
    const SpaceGroup *match = _SpaceGroups.sgbn.find(name_in) != _SpaceGroups.sgbn.end() ? _SpaceGroups.sgbn[name_in] : nullptr;
    if (match) return match;

    // If a match wasn't found, remove the white space and try again
    string name = RemoveWhiteSpaceUnderscore(name_in);
    match = _SpaceGroups.sgbn.find(name) != _SpaceGroups.sgbn.end() ? _SpaceGroups.sgbn[name] : nullptr;

    if (!match) {
      // Try another search, e.g. Fm-3m instead of Fm3m
      string search = name;
      bool hasMirror = (name.find('m') != string::npos || name.find('d') != string::npos || name.find('n') != string::npos || name.find('c') != string::npos);
      if (name.find('4') != string::npos && hasMirror && name.find('-') == string::npos) {
        search.insert(name.find('4'), "-");
      } else if (name.find('3') != string::npos && hasMirror && name.find('-') == string::npos) {
        search.insert(name.find('3'), "-");
      } else if (name.find('6') != string::npos && hasMirror && name.find('-') == string::npos) {
        search.insert(name.find('6'), "-");
      }

      match = _SpaceGroups.sgbn.find(search) != _SpaceGroups.sgbn.end() ? _SpaceGroups.sgbn[search] : nullptr;
    }

    return (match);
  }

  /*!
   */
  const SpaceGroup * SpaceGroup::GetSpaceGroup (unsigned id)
  {
    if (!_SpaceGroups.Inited())
      _SpaceGroups.Init();
    return (id > 0 && id <= 230)? _SpaceGroups.sgbi[id - 1].front() : nullptr;
  }

  /*!
   */
  void SpaceGroup::RegisterSpaceGroup (int nb, ...)
  {
    _SpaceGroups.sgs.insert(this);
    if (m_id > 0 && m_id <= 230)
      _SpaceGroups.sgbi[m_id - 1].push_back(this);
    if (m_HM.length() > 0)
	  {
        if (m_OriginAlternative != 0)
		  {
            char a = '0' + m_OriginAlternative;
            std::string nm = m_HM + ':' + a;
            if (_SpaceGroups.sgbn[nm] == nullptr)
              _SpaceGroups.sgbn[nm] = this;
            // Also use the symbol stripped from whitespaces as key
            std::string stripped_HM=RemoveWhiteSpaceUnderscore(nm);
            if (stripped_HM.length() > 0 && _SpaceGroups.sgbn[nm] == nullptr)
              _SpaceGroups.sgbn[nm] = this;
		  }
        if ((m_OriginAlternative & 1) == 0 && _SpaceGroups.sgbn[m_HM] == nullptr)
          _SpaceGroups.sgbn[m_HM] = this;
	  }
    // Also use the HM symbol stripped from whitespaces as key
	  std::string stripped_HM=RemoveWhiteSpaceUnderscore(m_HM);
    if (stripped_HM.length() > 0 && _SpaceGroups.sgbn[stripped_HM] == nullptr)
      _SpaceGroups.sgbn[stripped_HM] = this;
    if (m_Hall.length() > 0 && _SpaceGroups.sgbn[m_Hall] == nullptr)
      _SpaceGroups.sgbn[m_Hall] = this;
    if (nb == 0)
      return;
    va_list args;
    va_start(args, nb);
    string name;
    for (int i = 0; i < nb; i++)
      {
        name=va_arg(args, const char *);
        if (name.length() > 0 && _SpaceGroups.sgbn[name] == nullptr)
          _SpaceGroups.sgbn[name] = this;
      }
    va_end(args);
  }

  /*!
   */
  bool SpaceGroup::operator ==(const SpaceGroup &sg) const
  {
    if (m_transforms.size() != sg.m_transforms.size())
      return false;
    set<string> s0, s1;
    list<transform3d*>::const_iterator i, iend;
    iend = m_transforms.end();
    for (i = m_transforms.begin(); i != iend; ++i)
      s0.insert((*i)->DescribeAsString());
    iend = sg.m_transforms.end();
    for (i = sg.m_transforms.begin(); i != iend; ++i)
      s1.insert((*i)->DescribeAsString());
    if (s0.size() != s1.size())
      return false;
    set<string>::iterator j, jend = s0.end();
    for (j = s0.begin(); j != jend; ++j)
      if (s1.find(*j) == s1.end())
        return false;
    return true;
  }

  /*!
   */
  bool SpaceGroup::IsValid() const
  {
    if (!m_transforms.size())
      return false;
    list<transform3d*>::const_iterator i, iend = m_transforms.end();
    map <string, transform3d*>T;
    for (i = m_transforms.begin(); i != iend; ++i)
      {
        if (T.find((*i)->DescribeAsString()) != T.end())
          {
            cerr << "Duplicated transform: " << (*i)->DescribeAsString() << endl;
            return false;
          }
        T[(*i)->DescribeAsString()] = *i;
      }
		// calculate all products and check if they are in the group
		map <string, transform3d*>::iterator j, k, end = T.end();
    string s;
    bool has_inverse;
		for (j = T.begin(); j != end; ++j)
		  {
        has_inverse = false;
        for (k = T.begin(); k != end; ++k)
          {
            s = (*(*j).second * *(*k).second).DescribeAsString();
            if (T.find(s) == end)
              {
                cerr << "Invalid transform: " << (*j).first << " * " << (*k).first << " = " << s << endl;
                return false;
              }
            if (!has_inverse && s == "x,y,z")
              has_inverse = true;
          }
        if (!has_inverse)
          {
            cerr << "Transform with no inverse: " << (*j).first << endl;
            return false;
          }
      }
    return true;
  }

  /*!
   */
  const SpaceGroup * SpaceGroup::Find (SpaceGroup* group)
  {
    const SpaceGroup *found = nullptr;
    if (group->m_Hall.length() > 0 && _SpaceGroups.sgbn.find(group->m_Hall)!=_SpaceGroups.sgbn.end())
      {
        found = _SpaceGroups.sgbn[group->m_Hall];
        if (!found)
          obErrorLog.ThrowError(__FUNCTION__, "Unknown space group (Hall symbol:"+group->m_Hall+") error, please file a bug report.", obError);
        if (group->m_transforms.size() && *found  != *group)
          {
            unsigned id = group->GetId();
            if (id != 3 && id != 68) // these groups have duplicates
              {
                obErrorLog.ThrowError(__FUNCTION__, "Space group error (Hall symbol and list of transforms do not match), please file a bug report.", obWarning);
                return found;
              }
          }
        else
        /* even if there is an error (this should not occur) return the found group, since
           Hall names are secure */
          return found;
      }
    // Identify from the HM symbol, after removing all whitespaces or underscore (which are valid separators in
    // old CIF files)
    std::string stripped_hm=RemoveWhiteSpaceUnderscore(group->m_HM);
    if (stripped_hm.length() > 0 &&
        _SpaceGroups.sgbn.find(stripped_hm)!=_SpaceGroups.sgbn.end() &&
        (found = _SpaceGroups.sgbn[stripped_hm]))
      {
        if (*found == *group){
          found = _SpaceGroups.sgbn[found->GetHallName()];
          return found;
        }
        if (group->m_transforms.size())
          {// If transforms (symmetry operations) are listed, make sure they match the tabulated ones
            list<const SpaceGroup*>::const_iterator i, end = _SpaceGroups.sgbi[found->m_id - 1].end();
            for (i = _SpaceGroups.sgbi[found->m_id - 1].begin(); i!= end; ++i)
              if ((**i) == *group)
                return *i;
            obErrorLog.ThrowError(__FUNCTION__, "Unknown space group error (H-M symbol:"+group->m_HM+"), cannot match the list of transforms, please file a bug report.", obError);
            return nullptr;
          }
        else if (group->m_transforms.size() == 0)
          {// No transforms (symmetry operations) are listed, warn if HM symbol can match several spacegroups
            int n = 0;
            list<const SpaceGroup*>::const_iterator i, end = _SpaceGroups.sgbi[group->m_id].end();
            for (i = _SpaceGroups.sgbi[group->m_id].begin(); i!= end; ++i)
              if (RemoveWhiteSpaceUnderscore((*i)->m_HM) == stripped_hm)
                n++;
            if (n > 1)
              obErrorLog.ThrowError(__FUNCTION__, "Ambiguous space group: HM symbol corresponds to several space groups.", obWarning);
            return found;
          }
        /* even if there is an error (this should not occur) return the found group, since
           Hall names are secure */
      }
    else if (group->m_id > 0 && group->m_id <= 230)
      {
        if (group->m_transforms.size())
          {
            list<const SpaceGroup*>::const_iterator i, end = _SpaceGroups.sgbi[group->m_id - 1].end();
            for (i = _SpaceGroups.sgbi[group->m_id - 1].begin(); i!= end; ++i)
              if ((**i) == *group)
                return *i;
          }
        else if (group->m_transforms.size() == 0)
          {
            if (_SpaceGroups.sgbi[group->m_id - 1].size() > 1)
              obErrorLog.ThrowError(__FUNCTION__, "Ambiguous space group: sg number corresponds to several space groups.", obWarning);
            return _SpaceGroups.sgbi[group->m_id - 1].front();
          }
      }
    // If we are there, we need to make a hard search through the whole collection
    if (!group->IsValid())
      {
        obErrorLog.ThrowError(__FUNCTION__, "Unknown space group (HM:"+group->m_HM+",Hall:"+group->m_Hall
                                           +") with incomplete or wrong definition.", obWarning);
        return nullptr;
      }
    set<SpaceGroup*>::iterator i, end = _SpaceGroups.sgs.end();
    for (i = _SpaceGroups.sgs.begin(); i != end; ++i)
      if (**i == *group)
        return *i;
    obErrorLog.ThrowError(__FUNCTION__, "Unknown space group error, please file a bug report.", obWarning);
    return nullptr;
	}
}

//! \file spacegroup.cpp
//! \brief Handle Crystallographic Space Groups

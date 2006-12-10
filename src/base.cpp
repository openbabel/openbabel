/**********************************************************************
base.cpp - Base classes to build a graph
 
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2006 by Geoffrey R. Hutchison
 
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
#include <openbabel/base.h>

using namespace std;

//! Global namespace for all Open Babel code
namespace OpenBabel
{
  bool OBBase::HasData(const string &s)
  {
    if (_vdata.empty())
      return(false);

    vector<OBGenericData*>::iterator i;

    for (i = _vdata.begin();i != _vdata.end();++i)
      if ((*i)->GetAttribute() == s)
        return(true);

    return(false);
  }

  bool OBBase::HasData(const char *s)
    //returns true if the generic attribute/value pair exists
  {
    if (_vdata.empty())
      return(false);

    vector<OBGenericData*>::iterator i;

    for (i = _vdata.begin();i != _vdata.end();++i)
      if ((*i)->GetAttribute() == s)
        return(true);

    return(false);
  }


  bool OBBase::HasData(const unsigned int dt)
    //returns true if the generic attribute/value pair exists
  {
    if (_vdata.empty())
      return(false);

    vector<OBGenericData*>::iterator i;

    for (i = _vdata.begin();i != _vdata.end();++i)
      if ((*i)->GetDataType() == dt)
        return(true);

    return(false);
  }

  //! Returns the value given an attribute name
  OBGenericData *OBBase::GetData(const string &s)
  {
    vector<OBGenericData*>::iterator i;

    for (i = _vdata.begin();i != _vdata.end();++i)
      if ((*i)->GetAttribute() == s)
        return(*i);

    return(NULL);
  }

  //! Returns the value given an attribute name
  OBGenericData *OBBase::GetData(const char *s)
  {
    vector<OBGenericData*>::iterator i;

    for (i = _vdata.begin();i != _vdata.end();++i)
      if ((*i)->GetAttribute() == s)
        return(*i);

    return(NULL);
  }

  OBGenericData *OBBase::GetData(const unsigned int dt)
  {
    vector<OBGenericData*>::iterator i;
    for (i = _vdata.begin();i != _vdata.end();++i)
      if ((*i)->GetDataType() == dt)
        return(*i);
    return(NULL);
  }

  void OBBase::DeleteData(unsigned int dt)
  {
    vector<OBGenericData*> vdata;
    vector<OBGenericData*>::iterator i;
    for (i = _vdata.begin();i != _vdata.end();++i)
      if ((*i)->GetDataType() == dt)
        delete *i;
      else
        vdata.push_back(*i);
    _vdata = vdata;
  }

  void OBBase::DeleteData(vector<OBGenericData*> &vg)
  {
    vector<OBGenericData*> vdata;
    vector<OBGenericData*>::iterator i,j;

    bool del;
    for (i = _vdata.begin();i != _vdata.end();++i)
      {
        del = false;
        for (j = vg.begin();j != vg.end();++j)
          if (*i == *j)
            {
              del = true;
              break;
            }
        if (del)
          delete *i;
        else
          vdata.push_back(*i);
      }
    _vdata = vdata;
  }

  void OBBase::DeleteData(OBGenericData *gd)
  {
    vector<OBGenericData*>::iterator i;
    for (i = _vdata.begin();i != _vdata.end();++i)
      if (*i == gd)
        {
          delete *i;
          _vdata.erase(i);
        }

  }

  /*! \mainpage

  \section intro Introduction and Key Modules
 
  Open Babel is a full chemical software toolbox. In addition to converting
  file formats, it offers a complete programming library for developing 
  chemistry software. The library is written primarily in C++ and also offers
  interfaces to other languages (e.g., Perl and Python) using essentially
  the same API.

  The heart of Open Babel lies in the \link OpenBabel::OBMol OBMol\endlink, 
  \link OpenBabel::OBAtom OBAtom\endlink, and 
  \link OpenBabel::OBBond OBBond\endlink classes,
  which handle operations on atoms, bonds and molecules. Newcomers should 
  start with looking at the \link OpenBabel::OBMol OBMol\endlink class, 
  designed to store the basic information
  in a molecule and to perceive information about a molecule.

  One of the key philosophies in the code is that transformations and
  automatic perception of properties are performed in a <a href="http://en.wikipedia.org/wiki/Lazy_evaluation">"lazy"</a>
  manner. That is, until you call for partial atomic charges, no
  charges are calculated. This ensures faster transformations of
  chemical data -- properties that are not needed for your code will
  typically not be calculated. When such data is needed, appropriate
  routines are called, and a "flag" is set (e.g., via OBMol::SetFlag
  or OBAtom::SetFlag etc.) so that the code is only run once.

  Arbitrary custom data and text descriptors can be stored in any atom,
  bond, molecule, or residue using the \link OpenBabel::OBGenericData
  OBGenericData\endlink or \link OpenBabel::OBPairData
  OBPairData\endlink classes.

  Conversion between various chemical file formats is accomplished through
  the \link OpenBabel::OBConversion OBConversion\endlink and \link 
  OpenBabel::OBFormat OBFormat\endlink classes, often through use of the \link 
  OpenBabel::OBMoleculeFormat OBMoleculeFormat\endlink subclass which is designed
  for easy read/write access to one or more \link OpenBabel::OBMol OBMol\endlink
  objects. The philosophy of the file format codes is to parse as much
  chemical information from a given file as possible (no data left
  behind) and ideally any perception or transformations will occur when
  writing to some other format later.

  \section other Other Resources

  Open Babel is a community project. In addition to this API documentation,
  the website offers a variety of up-to-date and useful information for
  developing with the library.

  For more resources, check out the <a href="http://openbabel.sourceforge.net/wiki/Develop">Developing with Open Babel</a> section,
  particularly, further code examples in the various <a href="http://openbabel.sourceforge.net/wiki/Developer:Tutorial">developer tutorials</a>.

  */

} // namespace OpenBabel

//! \file base.cpp
//! \brief Implementation of base classes.

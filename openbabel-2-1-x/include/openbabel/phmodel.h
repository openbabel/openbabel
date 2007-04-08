/**********************************************************************
phmodel.h - Read pH rules and assign charges.
 
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2005 by Geoffrey R. Hutchison
 
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

#ifndef OB_PHMODEL_H
#define OB_PHMODEL_H

#include <openbabel/parsmart.h>
#include <openbabel/data.h>

namespace OpenBabel
{

// class introduction in phmodel.cpp
class OBAPI OBChemTsfm
{
    std::vector<int>                            _vadel;
    std::vector<std::pair<int,int> >            _vele;
    std::vector<std::pair<int,int> >            _vchrg;
    std::vector<std::pair<int,int> >            _vbdel;
    std::vector<std::pair<std::pair<int,int>,int> >  _vbond;
    OBSmartsPattern _bgn,_end;
public:
    OBChemTsfm()    {}
    ~OBChemTsfm()   {}
    //! Initialize this transformation with the supplied SMARTS patterns
    bool Init(std::string&start, std::string &end);
    //! Apply this transformation to all matches in the supplied OBMol
    bool Apply(OBMol&);
};

//! \brief Corrections for pH used by OBMol::CorrectForPH()
class OBPhModel : public OBGlobalDataBase
{
    std::vector<std::vector<int> >                      _mlist;
    std::vector<OBChemTsfm*>                            _vtsfm;
    std::vector<std::pair<OBSmartsPattern*,std::vector<double> > > _vschrg;
public:
    OBPhModel();
    ~OBPhModel();

    void ParseLine(const char*);
    //! \return the number of chemical transformations
    unsigned int GetSize()                 { return _vtsfm.size();}
    void AssignSeedPartialCharge(OBMol&);
    void CorrectForPH(OBMol&);
};



} //namespace OpenBabel

#endif // OB_PHMODEL_H

//! \file phmodel.h
//! \brief Read pH rules and assign charges.

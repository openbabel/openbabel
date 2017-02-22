/**********************************************************************
phmodel.h - Read pH rules and assign charges.

Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2005 by Geoffrey R. Hutchison

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
    /*! Is this transformation an acid dissociation?
     *  \code
     *      Ka
     *  HA ----> A(-)         (the H(+) will be deleted)
     *  \endcode
     *
     *  IsAcid() will check the charge in the end SMARTS pattern.
     *  \return true if the charge is less than 0 (-1).
     */
    bool IsAcid();
    /*! Is this a transformation to the conjugated acid from a base?
     *  \code
     *      Ka
     *  HA ----> A(-)         (the H(+) will be deleted)
     *  \endcode
     *
     *  IsBase() will check the charge in the end SMARTS pattern.
     *  \return true if the charge is higher than 0 (+1).
     */
    bool IsBase();
};

/*! \brief Corrections for pH used by OBMol::CorrectForPH()
 *
 *  The data/phmodel.txt file contains transformations which are applied
 *  to correct the charges for a given pH. This function uses the
 *  Henderson-Hasselbalch equation to calculate which species (protonated/
 *  unprotonated) is present in the highest concentration at the given pH.
 *
 *  For acids an entry would look like:
 *  \code
 *  # carboxylic acid
 *  O=C[OD1:1] >> O=C[O-:1]	4.0
 *  \endcode
 *
 *  The 4.0 is the pKa for the dissociation [HA] -> [H+] + [A-]. To
 *  calculate [HA]/[A-] we use:
 *  \code
 *  [HA] / [A-] = 10^(pKa - pH)
 *
 *  [HA]/[A-] > 1  :  [HA] > [A-]
 *  [HA]/[A-] < 1  :  [A-] > [HA]
 *  \endcode
 *
 *  For a base, an entry would look be:
 *  \code
 *  # methyl amine
 *  C[N:1] >> C[N+:1]	10.7
 *  \endcode
 *
 *  Here, the 10.7 is the pKa for the dissociation [BH+] -> [H+] + [B:]. To
 *  calculate [BH+]/[B:] we use:
 *  \code
 *  [BH+] / [B:] = 10^(pKa - pH)
 *
 *  [BH+]/[B:] > 1  :  [BH+] > [B:]
 *  [BH+]/[B:] < 1  :  [B:] > [BH+]
 *  \endcode
 *
 *  The transformations are all applied (if needed at the specified pH value) in
 *  the same order they are found in data/phmodel.txt.
 */
class OBAPI OBPhModel : public OBGlobalDataBase
{
    std::vector<OBChemTsfm*>                            _vtsfm;
    std::vector<double>                                 _vpKa;
    std::vector<std::pair<OBSmartsPattern*,std::vector<double> > > _vschrg;
public:
    OBPhModel();
    ~OBPhModel();

    void ParseLine(const char*);
    //! \return the number of chemical transformations
    size_t GetSize()                 { return _vtsfm.size();}
    void AssignSeedPartialCharge(OBMol&);
    //void CorrectForPH(OBMol&);
    void CorrectForPH(OBMol&, double pH  = 7.4 );
};



} //namespace OpenBabel

#endif // OB_PHMODEL_H

//! \file phmodel.h
//! \brief Read pH rules and assign charges.


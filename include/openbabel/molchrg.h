/**********************************************************************
molchrg.h - Assign Gasteiger partial charges.

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

#ifndef OB_MOLCHRG_H
#define OB_MOLCHRG_H

#include <openbabel/babelconfig.h>
#include <vector>

namespace OpenBabel
{

  // Foward declarations
  class OBMol;
  class OBAtom;

//! \class GasteigerState molchrg.h <openbabel/molchrg.h>
//! \brief Helper class for OBGastChrg which stores the Gasteiger states of a given atom
class OBAPI GasteigerState
{
public:
    GasteigerState();
    ~GasteigerState() {}

    void SetValues(double _a,double _b,double _c,double _q)
    {
        a = _a;
        b = _b;
        c = _c;
        denom=a+b+c;
        q = _q;
    }

    double a, b, c;
    double denom;
    double chi;
    double q;
};

// class introduction in molchrg.cpp
class OBAPI OBGastChrg
{
  std::vector <GasteigerState*> _gsv; //!< vector of internal GasteigerState (for each atom)

    void InitialPartialCharges(OBMol &);
    bool GasteigerSigmaChi(OBAtom *,double &,double &,double &);

public:
    OBGastChrg()    {}
    ~OBGastChrg();

    //! Assign partial charges to this OBMol, setting each OBAtom partial charge
    bool AssignPartialCharges(OBMol &);
    //! Resize the internal GasteigerState vector to match the size of the molecule
    void GSVResize(int);
};

}

#define OB_GASTEIGER_DENOM  20.02
#define OB_GASTEIGER_DAMP   0.5
#define OB_GASTEIGER_ITERS  6

#endif // OB_MOLCHRG_H

//! \file molchrg.h
//! \brief Assign Gasteiger partial charges.

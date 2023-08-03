/**********************************************************************
patty.h - Programmable atom typer.

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

#ifndef OB_PATTY_H
#define OB_PATTY_H

#include <openbabel/parsmart.h>

namespace OpenBabel
{
#define PT_CATION      1
#define PT_ANION       2
#define PT_ACCEPTOR    3
#define PT_POLAR       4
#define PT_DONOR       5
#define PT_HYDROPHOBIC 6
#define PT_OTHER       7
#define PT_METAL	   8

// class introduction in patty.cpp
class OBAPI patty
{
    std::vector<OBSmartsPattern*> _sp;
    std::vector<std::string> smarts;
    std::vector<std::string> typ;
    bool debug;

public :

    patty()
    {
        debug = false;
    }
    patty(char *s)
    {
        debug = false;
        read_rules(std::string(s));
    }

    patty(const std::string &s)
    {
        debug = false;
        read_rules(s);
    }

    ~patty()
    {
        std::vector<OBSmartsPattern*>::iterator i;
        for (i = _sp.begin();i != _sp.end();++i)
            delete *i;
    }
    void debug_on()
    {
        debug = true;
    }
    void debug_off()
    {
        debug = false;
    }
    void read_rules(const std::string &infile);
    void assign_rules(std::vector<std::string> &rules);
    void assign_types(OBMol &mol,std::vector<std::string> &atm_typ);
    void assign_types(OBMol &mol,std::vector<int> &atm_typ);
    int type_to_int(const std::string &type, bool failOnUndefined= false);
    int Istype(const std::string &type);//!< return atom type index, 0 otherwise
};

} // end namespace OpenBabel

#endif // OB_PATTY_H

//! \file patty.h
//! \brief Programmable atom typer. (abbreviated P.At.Ty.)

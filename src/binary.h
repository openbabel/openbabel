/**********************************************************************
binary.h - Handle rotamer list and binary base files.

Copyright (C) 1998-2000 by OpenEye Scientific Software, Inc.
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

#ifndef OB_BINARY_H
#define OB_BINARY_H

#include <vector>
#include <map>

#if HAVE_FSTREAM
#include <fstream>
#elif HAVE_FSTREAM_H
#include <fstream.h>
#endif

#include "mol.h"
#include "rotor.h"
#include "generic.h"

namespace OpenBabel
{

class OBRotamerList : public OBGenericData
{
    unsigned int                         _NBaseCoords;
    std::vector<double*>                       _c;
    std::vector<std::vector<double> >               _vres;
    std::vector<unsigned char*>               _vrotamer;
    std::vector<std::pair<OBAtom**,std::vector<int> > > _vrotor;

    OBRotamerList(const OBRotamerList &cpy) : OBGenericData(cpy) {}
    OBRotamerList& operator =(const OBRotamerList &) { return *this; }

public:
    OBRotamerList() {_NBaseCoords=0; _type=obRotamerList; _attr="RotamerList";}
    ~OBRotamerList();
    void Setup(OBMol&,OBRotorList&);
    void Setup(OBMol&,unsigned char*,int);
    unsigned int NumRotors()   const { return (unsigned int)_vrotor.size();   }
    unsigned int NumRotamers() const { return (unsigned int)_vrotamer.size(); }
    void AddRotamer(double*);
    void AddRotamer(int *arr);
    void AddRotamer(unsigned char *arr);
    void AddRotamers(unsigned char*,int);
    void GetReferenceArray(unsigned char*);
    void ExpandConformerList(OBMol&,std::vector<double*>&);
    std::vector<unsigned char*>::iterator BeginRotamer() {return _vrotamer.begin();}
    std::vector<unsigned char*>::iterator EndRotamer()   {return _vrotamer.end();  }

    // Support for internal storage of base coordinate sets that
    // rotamers operate on

    //! Create a conformer list using the internal base set of coordinates
    std::vector<double*> CreateConformerList(OBMol& mol);

    //! Copies the mol's conformers (the coordinates, NOT the pointers)
    //! into the object as base coordinates
    void SetBaseCoordinateSets(OBMol& mol) {
        SetBaseCoordinateSets(mol.GetConformers(),mol.NumAtoms());
    } 

    //! Copies the coordinates in bc, NOT the pointers, into the object
    void SetBaseCoordinateSets(std::vector<double*> bc, unsigned int N); 

    unsigned int NumBaseCoordinateSets() const {
        return (unsigned int)_c.size();
    }

    //! Get a pointer to a specific base pointer
    double *GetBaseCoordinateSet(unsigned int i) {
        return (i<_c.size()) ? _c[i] : NULL;
    }

    unsigned int NumAtoms() const { return _NBaseCoords; }
};

//! \brief Store a compressed binary database of OBMol as OEBinary files
class OBBinaryDBase
{
    std::ifstream		 _ifs;
    std::vector<std::streampos>  _vpos;

public:

     OBBinaryDBase(const char*);
     OBBinaryDBase(std::string&);
    ~OBBinaryDBase();

    int  Size();
    void GetMolecule(OBMol&, int);
};

int Swab(int);

}

#endif // OB_BINARY_H


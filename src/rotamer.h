/**********************************************************************
rotamer.h - Handle rotamer list data.
 
Copyright (C) 1998-2000 by OpenEye Scientific Software, Inc.
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

#ifndef OB_ROTAMER_H
#define OB_ROTAMER_H

#include <vector>
#include <map>

#include "mol.h"
#include "rotor.h"
#include "generic.h"

namespace OpenBabel
{

  //! \brief Supports a set of rotamer coordinate sets for some number of potentially rotatable bonds
  // Further class introduction in rotamer.cpp
 class OBAPI OBRotamerList : public OBGenericData
  {
    //! Number of atoms in the base coordinate set (i.e., OBMol::NumAtoms())
    unsigned int                         _NBaseCoords;
    //! Base coordinate sets (i.e., existing conformers to be modified)
    std::vector<double*>                 _c;
    //! Individual bond rotors (from an OBRotor object or other)
    std::vector<std::pair<OBAtom**,std::vector<int> > > _vrotor;
    //! \brief Index of each rotor's different sampling states ("resolution")
    //! Usually from OBRotor::GetResolution()
    std::vector<std::vector<double> >    _vres;
    //! Individual rotamer states (i.e., the array of rotor settings)
    std::vector<unsigned char*>          _vrotamer;

    /*Because contains OBAtom*, these aren't meaningful without knowing the parent molecule
      OBRotamerList(const OBRotamerList &cpy) : OBGenericData(cpy)
      {}
      OBRotamerList& operator =(const OBRotamerList &)
      {
      return *this;
      }
		*/

  public:
    OBRotamerList()
      {
        _NBaseCoords=0;
        _type= OBGenericDataType::RotamerList;
        _attr="RotamerList";
      }
		virtual OBGenericData* Clone(OBBase* parent) const;

    ~OBRotamerList();
    //! Set up a rotamer list based on an already created OBRotorList
    void Setup(OBMol&,OBRotorList&);
    //! Set up a rotamer list based on the supplied reference atoms and the number of rotors
    //! \param mol The molecule to evaluate
    //! \param ref An array of the 4 dihedral atoms for each rotor
    //! \param nrotors The number of rotors (i.e., the size of ref / 4)
    void Setup(OBMol &mol,unsigned char*ref,int nrotors);
    //! \return the number of rotatable bonds considered
    unsigned int NumRotors()   const
    {
      return (unsigned int)_vrotor.size();
    }
    //! \return the number of rotamer (conformation) coordinate sets
    unsigned int NumRotamers() const
    {
      return (unsigned int)_vrotamer.size();
    }
    //! Add a rotamer to the list based on the supplied coordinate set as a double*
    void AddRotamer(double*);
    //! Add a rotamer to the list based on @p key as a configuration of the individual rotor bonds
    void AddRotamer(int *key);
    //! Add a rotamer to the list based on @p key as a configuration of the individual rotor bonds
    void AddRotamer(unsigned char *key);
    //! Add @p nconf rotamers based on @p as an array of configurations much like AddRotamer()
    void AddRotamers(unsigned char *arr,int nconf);
    void GetReferenceArray(unsigned char*) const;

    //! \name Iterator methods
    //@{
    std::vector<unsigned char*>::iterator BeginRotamer()
      {
        return _vrotamer.begin();
      }
    std::vector<unsigned char*>::iterator EndRotamer()
      {
        return _vrotamer.end();
      }
    //@}

    //! Create a conformer list using the internal base set of coordinates
    //! \return The set of coordinates by rotating the bonds in each rotamer
    std::vector<double*> CreateConformerList(OBMol& mol);

    //! Create a conformer list using the internal base set of coordinates
    //! Returns the set of coordinates as a reference in @p confs
    void ExpandConformerList(OBMol&mol,std::vector<double*>& confs);

    //! Copies the mol's conformers (the coordinates, NOT the pointers)
    //! into the object as base coordinates
    void SetBaseCoordinateSets(OBMol& mol)
    {
      SetBaseCoordinateSets(mol.GetConformers(), mol.NumAtoms());
    }

    //! Copies the coordinates in bc, NOT the pointers, into the object
    void SetBaseCoordinateSets(std::vector<double*> bc, unsigned int N);

    unsigned int NumBaseCoordinateSets() const
    {
      return (unsigned int)_c.size();
    }

    //! Get a pointer to a specific base pointer
    double *GetBaseCoordinateSet(unsigned int i) const
    {
      return (i<_c.size()) ? _c[i] : NULL;
    }

    unsigned int NumAtoms() const
    {
      return _NBaseCoords;
    }
  };

  int Swab(int);

}

#endif // OB_ROTAMER_H

//! \file rotamer.h
//! \brief Handle rotamer list data.

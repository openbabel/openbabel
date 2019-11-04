/**********************************************************************
rotamer.h - Handle rotamer list data.

Copyright (C) 1998-2000 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2006 by Geoffrey R. Hutchison

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

#ifndef OB_ROTAMER_H
#define OB_ROTAMER_H

#include <vector>
#include <map>

#include <openbabel/rotor.h>
#include <openbabel/generic.h>

namespace OpenBabel
{
  class OBMol;

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
    //! Rotors in rings
    std::vector<std::vector<int> >       _vrings;
    //! Dihedral angles of ring bonds
    std::vector<std::vector<double> >    _vringTors;

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
    void AddRotamer(std::vector<int> key);
    //! Add a rotamer to the list based on @p key as a configuration of the individual rotor bonds
    void AddRotamer(unsigned char *key);
    //! Add @p nconf rotamers based on @p as an array of configurations much like AddRotamer()
    void AddRotamers(unsigned char *arr,int nconf);
    //! \return A reference array (as used by AddRotamer() as a configuration of the individual rotor bonds
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

    //! \brief Create a conformer list using the internal base set of coordinates
    //! \return The set of coordinates by rotating the bonds in each rotamer
    std::vector<double*> CreateConformerList(OBMol& mol);

    //! \brief Create a conformer list using the internal base set of coordinates
    //! \return The set of coordinates as a reference in @p confs
    void ExpandConformerList(OBMol&mol,std::vector<double*>& confs);

    bool SetCurrentCoordinates(OBMol &mol, std::vector<int> arr);

    //! \brief Copies the mol's conformers (the coordinates, NOT the pointers)
    //! into the object as base coordinates
    void SetBaseCoordinateSets(OBMol& mol);

    //! Copies the coordinates in bc, NOT the pointers, into this object
    /** \param bc The conformer set for the molecule
        \param N  The number of atoms in the molecule
     **/
    void SetBaseCoordinateSets(std::vector<double*> bc, unsigned int N);

    //! \return The number of "base" coordinate sets (i.e., the number of conformers in the base OBMol)
    unsigned int NumBaseCoordinateSets() const
    {
      return static_cast<unsigned int> (_c.size());
    }

    //! Get a pointer to a specific base pointer (i.e., specific conformer)
    double *GetBaseCoordinateSet(unsigned int i) const
    {
      return (i<_c.size()) ? _c[i] : NULL;
    }

    //! \return The number of atoms in the base OBMol
    unsigned int NumAtoms() const
    {
      return _NBaseCoords;
    }
  };

  //! Swap Byte instruction (i.e., handle transfers between endian forms)
  int Swab(int);

}

#endif // OB_ROTAMER_H

//! \file rotamer.h
//! \brief Handle rotamer list data.

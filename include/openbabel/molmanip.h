/**********************************************************************
molmanip.h - Manipulate the geometries of molecules. Declaration of OBMolManip

Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2006 by Geoffrey R. Hutchison
Some portions Copyright (C) 2003 by Michael Banck
Some portions Copyright (C) 2017 by Torsten Sachse

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

#ifndef OB_MOLMANIP_H
#define OB_MOLMANIP_H

#include <stdexcept>
#include <queue>
#include <openbabel/mol.h>

//The following macro can be used to perform sanity checks
//if desired. Not used so far but example given below.
#define MOLMANIPSANITY(fallback)
//if(_manip==NULL){ \
//    std::cerr << "ERROR in method of OBMolManip: empty molecule pointer" << std::endl; \
//    return fallback; }

namespace OpenBabel
{

  class OBMol;

  // Class OBMolManip

  // class introduction in molmanip.cpp
 class OBAPI OBMolManip
  {
  protected:
    OBMol*        _manip; //!< The molecule that can be manipulated

  public:

    //! \name Initialization methods
    //@{
    //Disallow the empty constructor OBMolManip();
    //! Assignment constructor, directly assigns the molecule to be manipulated
    OBMolManip(OBMol* mol) : _manip(mol) {
        if(_manip==NULL){ throw std::runtime_error(
            "ERROR in constructor of OBMolManip: empty molecule pointer.");}
    }
    //! Destructor
    virtual ~OBMolManip();
    //! Get the pointer to the current molecule
    OBMol* GetMolecule(){return _manip;}
    //@}

    //! \name Molecule modification methods
    //@{
    //! Rotate a molecule using only doubles (only current conformer)
    void Rotate(const double axis[3], const double angle);
    //! Rotate using a vector3 object as axis (only current conformer)
    void Rotate(const vector3 axis, const double angle);
    //! Rotate only the atoms in the given vector (only current conformer)
    void Rotate(std::vector<OBAtom*> &atom_vec, const matrix3x3 m);
    //! Rotate a molecule around one of its main axes (only current conformer)
    void Rotate(const int main_axis_nr, const double angle);
    //! Rorate only the atoms in the given vector around their main axis (only current conformer)
    void Rotate(std::vector<OBAtom*> &atom_vec, const int main_axis_nr, const double angle);
    //! Rotate all conformers using the supplied matrix @p u (a 3x3 array of double)
    void Rotate(const double u[3][3]);
    //! Rotate all conformers using the supplied matrix @p m (a linear 3x3 row-major array of double)
    void Rotate(const double m[9]);
    //! Rotate a specific conformer @p nconf using the supplied rotation matrix @p m
    void Rotate(const double m[9],int nconf);
    //! Translate a molecule using only doubles (only current conformer)
    void Translate(const double vector[3]);
    //! Translate only the atoms in the given vector (only current conformer)
    void Translate(std::vector<OBAtom*> &atom_vec, const double vector[3]);
    //! Translates all conformers in the molecule by the supplied vector
    void Translate(const vector3 &v);
    //! Translates one conformer in the molecule by the supplied vector
    void Translate(const vector3 &v, int conf);
    //! Align a molecule with a plane as good as possible
    //! The plane is given by the poitn p which will be the molecule's new centre
    //! and two axes v1 and v2 which are the third main axis and the second main axis respectively
    void Align(const double p[3], const double v1[3], const double v2[3]);
    //! Align a molecule with a plane as good as possible
    //! Use vector3 objects as axes
    void Align(const vector3 p, const vector3 v1, const vector3 v2);
    //! Align only the atoms in the given vector with their main axes
    void Align(std::vector<OBAtom*> &atom_vec, const vector3 p, const vector3 v1, const vector3 v2);
    //! Set the dihedral angle defined by the four atom indices to angle
    //! A cis-configurations corresponds to angle=0 whereas a trans configuration corresponds to angle=180
    bool SetDihedralAngle(const int idxa1, const int idxa2, const int idxa3, const int idxa4, const double angle);
    double GetDihedralAngle(const int idxa1, const int idxa2, const int idxa3, const int idxa4);
    //! Set the angle defined by three atoms to angle
    //! Atom 1 will be kept fixed in its position
    bool SetAngle(const int idxa1, const int idxa2, const int idxa3, const double angle);
    double GetAngle(const int idxa1, const int idxa2, const int idxa3);
    //! Partition the molecule using a plane in Hessian normal form
    //! Only keep those atoms on the side of the plane in which @a direction points (WARNING: bonds will not be copied over)
    void PartMolecule(OBMol &dest, const double direction[3], const double point[3]);
    //! Mirrors the molecule at a point (inversion, if normal=(0,0,0)) or plane in Hessian normal form
    void Mirror(const double normal[3], const double point[3], bool center_it=false);
    //! Mirrors the molecule using vector3 objects
    void Mirror(const vector3 normal, const vector3 point, bool center_it=false);
    //! Mirror only the atoms in the given vector
    void Mirror(std::vector<OBAtom*> &atom_vec, const vector3 normal, const vector3 point, bool center_it=false);
    //! Bonds another molecule to the main one using the provided atom indices.
    //! The indices i1 and i2 are for the current molecule (i.e., this) and m1 and m2 for the to-be-connected one.
    //! Everything connected to the atoms indexed by i2 and m2 will be cut off. The atoms indexed by i1 and m1
    //! will be part of the new molecule, the atoms indexed by i2 and m2 will not.
    //! This automatically rotated 'mol' so that the bonds i1-i2 and m1-m2 are parallel to each other.
    bool BondMolecule(OBMol &mol, int i1, int i2, int m1, int m2);
    //! Cleave off part of a molecule that is connected to the atom indexed by i2 but leaving present
    //! everything that is connected to the atom indexed by i1.
    bool CleaveOff(int i1, int i2);
    //! Append the atoms and bonds stored in one molecule to another. This does not automatically bond
    //! those two bits together.
    void AppendMolecule(OBMol &mol);
    //! Add a new set of coordinates @p f as a new conformer
    //! You can also force to deep-copy the data
    void AddConformer(double *f, bool deepcopy=false);
    //! Delete a set of conformers in molecule
    void DeleteConformers(int start_idx, int end_idx);
    //@}

  };

} // end namespace OpenBabel

#endif // OB_MOLMANIP_H

//! \file molmanip.h
//! \brief Manipulate the geometries of molecules. Declaration of OBMolManip

/**********************************************************************
cpcomplex.h - Handle and define Cp (cyclopentadienyl) structures.

Copyright (C) 2023 by Jesus N. M.

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

#ifndef OB_CPCOMPLEX_H
#define OB_CPCOMPLEX_H


#include <openbabel/babelconfig.h>

#include <openbabel/mol.h>
#include <openbabel/math/vector3.h>
#include <openbabel/oberror.h>
#include <openbabel/bond.h>



namespace OpenBabel {

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define deg2rads(deg) ((deg) * M_PI / 180.0)

#define PERSPECTIVE_DEG (75.0)

    class BranchBlock {
        unsigned int _idx;                          //!< Block identifier
        std::vector<unsigned int> vidx_atoms;       //!< Vector idx of the atoms that are part of the block.

    public:
        //! Default constructor
        BranchBlock() {
            _idx = 0;
            vidx_atoms.clear();
        }
        //! Destructor
        ~BranchBlock() {}
        //! \return the size of the block (number of atoms in the block)
        int Size(){ return(vidx_atoms.empty() ? 0 : vidx_atoms.size()); }
        //! \return the block identifier
        unsigned int GetIdx() { return(_idx); }
        //! Set the block identifier
        void SetIdx(int idx) { _idx = idx; }
        //! Add an atom's idx to the block
        void AddAtom(int i){ vidx_atoms.push_back(i); }
        //! \return the idx of the atom at position @p i. Zero based access.
        unsigned int GetAtomIdx(int i)
        {
            if ((unsigned)i < 0 || (unsigned)i >= vidx_atoms.size())
            {
                obErrorLog.ThrowError(__FUNCTION__, "Requested AtomBlock Out of Range", obDebug);
            }

            return(vidx_atoms[i]);
        }

        //! \return Whether the @p idx exists within the atoms already inserted in the block
        bool HasAtom(int idx)
        {
            return ((std::find(vidx_atoms.begin(), vidx_atoms.end(), idx) != vidx_atoms.end())
                ? true : false);
        }

        //! \return whether or not it appears to be a Cp block
        // Cp will be possible if all the elements in the block are carbons up to that point and have a bond with an ogm metal.
        bool IsPossibleCp(OBMol &mol) {
            OBAtom* atom;
            bool test = false;
            for (int i = 0; i < vidx_atoms.size(); i++) {
                atom = mol.GetAtom(vidx_atoms[i]);
                if (!atom->IsCarbon()) {
                    return false;
                }
                else { //Si es un carbono, miramos si tiene algun bond con un Metal ogm
                    OBBond* bond;
                    OBBondIterator i;
                    for (bond = atom->BeginBond(i); bond && !test; bond = atom->NextBond(i)) {
                        test = bond->GetNbrAtom(atom)->IsOgmMetal();
                    }
                    if (!test) //Si no tiene ningun enlace con un ogm metal, devolvemos false
                        return false;
                    test = false; //Reseteamos el bool para el siguiente carbono
                }
            }
            return true; //Si llegamos aqui es porque todos los carbonos de la branch tienen al menos 1 enlace con metal
        }
    };

    



	class CpComplex {
		protected:
			OBMol* _parent;                         //!< Parent molecule
			unsigned int _idx;                      //!< Cp identifier within the molecule
			unsigned int metal_idx;                 //!< Atom idx of central metal
			std::vector<OBAtom*> _cpAtoms;          //!< Atoms for the carbons of the Cp structure
			std::vector<unsigned int> idx_carbons;  //!< Atom indexes for the carbons of the Cp structure
			vector3 center;                         //!< Cp center, for normal bond connection with metal atom, and aromatic circle position
            std::vector<vector3> circlePath;        //!< Coordinates for the cp circle (needed to achieve a perspective circunference)
			double radius;                          //!< Cp's aromatic circle radius
			unsigned int dummy_idx;                 //!< Dummy central atom idx


    public:
        //! Constructor
        CpComplex() {
            _parent = nullptr;
            _idx = 0;
            metal_idx = 0;
            _cpAtoms.clear();
            idx_carbons.clear();
            center.Set(0.0, 0.0, 0.0);
            radius = 0.0;
            dummy_idx = 0;
        }

        
        //! \name Methods to modify internal information
        //@{
        //! Attach an OBMol @p ptr as the parent container for this Cp
        void SetParent(OBMol* ptr) { _parent = ptr; }
        //! Set the center point of the Cp, sprecified by @p _v. It is equidistant to every carbon in th Cp, as they are disposed in a regular polygon
        void SetCentroid(vector3& _v) { center = _v; }
        //! Set the center point of the Cp, sprecified by @p v_x, v_y, v_z. It is equidistant to every carbon in th Cp, as they are disposed in a regular polygon
        void SetCentroid(const double v_x, const double v_y, const double v_z) { center.Set(v_x, v_y, v_z); }
        //! Set the radius of the Cp circle
        void SetRadius(double r) { radius = r; }
        //! Set the Cp identifier
        void SetIdx(int idx) { _idx = idx; }
        //! Set the idx of the central metal to which this Cp is attached
        void SetMetalIdx(int midx) { metal_idx = midx; }
        //! Dummy atom is created to make a perpendicular bond between the metal and the Cp drawing
        //! Set the atom idx of the dummy atom created for this Cp 
        void SetDummyIdx(int idx) { dummy_idx = idx; }
        //! Set the point of the Cp circle at index @p i to the coordinates specified by @p _v
        void SetCircleCoord(unsigned int i, vector3 _v);
        //! Set the point of the Cp circle at index @p i to the coordinates specified by @p _vx, _vy, _vz
        void SetCircleCoord(unsigned int i, double _vx, double _vy, double _vz = 0.0);
        //@}


        //! \name Methods to retrieve information
        //@{
        //! \return number of carbon atoms in the cp
        unsigned int GetCarbonsSize() { return idx_carbons.size(); }
        //! \return the molecule which contains this Cp, or NULL if none exists
        OBMol* GetParent() { return((OBMol*)_parent); }
        //! \return dummy atom idx for this Cp structure, or 0 if none exists
        unsigned int GetDummyIdx()   const { return((int)dummy_idx); }
        //! \return Cp identifier
        unsigned int GetIdx()   const { return((int)_idx); }
        //! \return Central metal identifier
        unsigned int GetMetalIdx()   const { return((int)metal_idx); }
        //! \return the whole contanier of carbon idx
        const std::vector<unsigned int>& GetIdxCarbons() { return idx_carbons; }
        //! \return the centroid of this Cp in a coordinate vector
        vector3& GetCentroid() { return center; };
        //! \return the radius of the Cp circle
        double GetRadius() { return radius; }
        //! \return the coordinate vector for the Cp circle point at position @p i in the container. Zero based access method to vector
        vector3 GetCircleCoord(unsigned int i);
        //! \return the number of points of the Cp circle
        int GetCirclePathSize() const { return circlePath.size(); }
        //! \return the whole container of coordinates of the Cp circle
        std::vector<vector3> GetCircleCoords() const { return circlePath; }
        //! \return carbon idx at position @p i in tha container. Zero based access method to vector
        unsigned int GetCarbonIdx(int i) const
        {
            if ((unsigned)i < 0 || (unsigned)i >= idx_carbons.size())
            {
                obErrorLog.ThrowError(__FUNCTION__, "Requested CarbonIdx Out of Range", obDebug);
            }

            return(idx_carbons[i]);
        }
        //@}


        //! \name Addition of data for a Cp
        //@{
        //! Adds a new atom idx to this Cp
        void AddIdxCarbon(int idx) { idx_carbons.push_back(idx); }
        //! Adds a new atom to this Cp
        void AddCpAtom(OBAtom* atom) { _cpAtoms.push_back(atom); }
        //! Adds a new point to the coordinate vector that forms the Cp circle
        void AddCircleCoord(vector3 _v) { circlePath.push_back(_v); }
        //@}


        //! \name Iteration methods
        //@{
        //! Set the iterator to the beginning of the Cp atom list
        //! \return the first atom, or NULL if none exist
        OBAtom* CpComplex::BeginAtomCp(OBAtomIterator& i)
        {
            i = _cpAtoms.begin();
            return i == _cpAtoms.end() ? nullptr : (OBAtom*)*i;
        }
        //! Advance the iterator to the next atom in the Cp
        //! \return the next first atom record, or NULL if none exist
        OBAtom* CpComplex::NextAtomCp(OBAtomIterator& i)
        {
            ++i;
            return i == _cpAtoms.end() ? nullptr : (OBAtom*)*i;
        }
        //@}


        //! \name Other operations
        //@{
        //! Calculate and set the centroid of this Cp, taking into consideration all atoms stored in _cpAtoms
        void FindCentroid();
        //! Equivalence operator
        bool operator==(const CpComplex* other) const { return (GetIdx() == other->GetIdx()); }
        //@}



    }; // class CpComplex

}// namespace OpenBabel

#endif   // OB_CPCOMPLEX_H

//! \file cpcomplex.h
//! \brief Handle Cp structures
























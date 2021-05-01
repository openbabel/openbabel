/**********************************************************************
atom.h - Handle OBAtom class.

Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2006 by Geoffrey R. Hutchison
Some portions Copyright (C) 2003 by Michael Banck

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

#ifndef OB_ATOM_H
#define OB_ATOM_H

#include <openbabel/babelconfig.h>

#ifndef OB_EXTERN
#  define OB_EXTERN extern
#endif

#include <vector>
#include <string>

#include <openbabel/base.h>
#include <openbabel/residue.h>
#include <openbabel/math/vector3.h>

namespace OpenBabel
{

  class OBBond;
  class OBMol;

  //! OBNodeBase is declared for backwards-compatibility with 2.0 and earlier code
  typedef OBAtom OBNodeBase;
  //! A standard iterator over a vector of bonds
  typedef std::vector<OBBond*>::iterator OBBondIterator;
  //! A standard iterator over a vector of atoms
  typedef std::vector<OBAtom*>::iterator OBAtomIterator;

  //ATOM Property Macros (flags)
  //! Atom is in a 4-membered ring
#define OB_4RING_ATOM     (1<<1)
  //! Atom is in a 3-membered ring
#define OB_3RING_ATOM     (1<<2)
  //! Atom is aromatic
#define OB_AROMATIC_ATOM  (1<<3)
  //! Atom is in a ring
#define OB_RING_ATOM      (1<<4)
  //! Atom is an electron donor
#define OB_DONOR_ATOM     (1<<7)
  //! Atom is an electron acceptor
#define OB_ACCEPTOR_ATOM  (1<<8)

#define SET_OR_UNSET_FLAG(X) \
  if (value) SetFlag(X); \
  else     UnsetFlag(X);

  // Class OBAtom
  // class introduction in atom.cpp
 #define OBATOM_TYPE_LEN 6
 class OBAPI OBAtom: public OBBase
    {
    protected:
      unsigned char                 _ele;       //!< atomic number (type unsigned char to minimize space -- allows for 0..255 elements)
      unsigned char                 _imph;      //!< number of implicit hydrogens
      char                          _type[OBATOM_TYPE_LEN];   //!< atomic type
      short                         _fcharge;   //!< formal charge
      unsigned short                _isotope;   //!< isotope (0 = most abundant)
      short                         _spinmultiplicity;//!< atomic spin, e.g., 2 for radical  1 or 3 for carbene

      unsigned int                  _idx;       //!< unique node index (GetIdx(), SetIdx())
      OBMol                        *_parent;    //!< parent molecule (if any)
      std::vector<OBBond*>          _vbond;     //!< bonds to this atom -- assumed to be one of the endpoints

      unsigned int                  _cidx;      //!< index into coordinate array
      unsigned short                _hyb;       //!< hybridization
      unsigned short                _flags;     //!< bitwise flags (e.g. aromaticity)
      double                        _pcharge;   //!< partial charge
      double                      **_c;         //!< coordinate array in double*
      mutable vector3               _v;         //!< coordinate vector
      OBResidue                    *_residue;   //!< parent residue (if applicable)

      unsigned long                 _id;        //!< unique id

      //! \return All flags
      int  GetFlag() const    {  return(_flags);  }
      //! Sets the bitwise @p flag
      void SetFlag(int flag)  { _flags |= flag;   }
      //! Unsets the bitwise @p flag
      void UnsetFlag(int flag) { _flags &= (~(flag)); }
      //! \return True of the atom has the @p flag
      bool HasFlag(int flag)  {  return((_flags & flag) ? true : false); }

    public:
      enum StereoFlag {

      };


       //! Used internally by graph traversal algorithms
      bool Visit;

      //! Constructor
      OBAtom();
      //! Destructor
      virtual ~OBAtom();
      //! Assignment
      OBAtom &operator = (OBAtom &);
      //! Equivalence
      bool operator==(const OBAtom * other) const {  return (GetIdx() == other->GetIdx()); }
      //! Duplicate another atom. Copies all information with the exception of index
      //! \since version 2.2
      void Duplicate(OBAtom *);
      //! Clear all data. Calls OBBase::Clear() to handle any generic data.
      //! \return True if successful.
      bool Clear();

      //! \name Methods to set atomic information
      //@{
      //! Set atom index (i.e., in an OBMol)
      void SetIdx(int idx)    { _idx = idx; _cidx = (idx-1)*3; }
      void SetId(unsigned long id) { _id = id; }
      //! Set atom hybridization (i.e., 1 = sp, 2 = sp2, 3 = sp3 ...)
      void SetHyb(int hyb)    { _hyb = hyb; }
      //! Set atomic number
      void SetAtomicNum(int atomicnum)    { _ele = (char)atomicnum; }
      //! Set isotope number (actual atomic weight is tabulated automatically, 0 = most abundant)
      void SetIsotope(unsigned int iso);
      //! Set the implicit hydrogen count to @p val
      void SetImplicitHCount(unsigned int val)    { _imph = (unsigned char)val; }
      //! Set the formal charge of the atom to @p fcharge
      void SetFormalCharge(int fcharge)   { _fcharge = fcharge; }
      //! Set the atomic spin to @p spin. See _spinmultiplicity
      void SetSpinMultiplicity(short spin){ _spinmultiplicity = spin; }
      //! Set the atomic type symbol (see OBTypeTable and OBAtomTyper for more)
      void SetType(const char *type);
      //! Set the atomic type symbol (see OBTypeTable and OBAtomTyper for more)
      void SetType(const std::string &type);
      //! Set the partial charge to @p pcharge
      void SetPartialCharge(double pcharge){ _pcharge = pcharge; }
      //! Set the coordinate vector for this atom to @p v as a vector3
      void SetVector(const vector3 &v);
      //! Set the coordinate vector for this atom based on @p x @p y & @p z
      void SetVector(const double x,const double y,const double z);
      //! Set the position of this atom from a pointer-driven array of coordinates
      void SetCoordPtr(double **c)        { _c = c; _cidx = (GetIdx()-1)*3; }
      //! Set the position of this atom based on the internal pointer array (i.e. from SetCoordPtr() )
      void SetVector();
      //! Attach an OBResidue @p res as containing this atom
      void SetResidue(OBResidue *res)     { _residue=res; }
      //! Attach an OBMol @p ptr as the parent container for this atom
      void SetParent(OBMol *ptr)          { _parent=ptr; }
      //! Mark atom as being aromatic
      void SetAromatic(bool value=true)                  { SET_OR_UNSET_FLAG(OB_AROMATIC_ATOM); }
      //! Mark an atom as belonging to at least one ring
      void SetInRing(bool value=true)         { SET_OR_UNSET_FLAG(OB_RING_ATOM); }
      //! Clear the internal coordinate pointer
      void ClearCoordPtr()     { _c = nullptr; _cidx=0; }
      //@}

      //! \name Methods to retrieve atomic information
      //@{
      //! \return the formal charge for this atom
      int          GetFormalCharge()  const { return(_fcharge);    }
      //! \return the atomic number for this atom
      unsigned int GetAtomicNum()     const { return((unsigned int)_ele); }
      //! \return the isotope for this atom, if specified, or 0 for unspecified
      unsigned short int GetIsotope() const { return(_isotope);    }
      //! \return the atomic spin, e.g., 0 (default) for normal atoms - note that this value is a convention,
      //!   2 for radical  1 or 3 for carbene
      int          GetSpinMultiplicity() const { return(_spinmultiplicity); }
      //! \return the atomic mass of this atom given by standard IUPAC
      //!  average molar mass
      double     GetAtomicMass()    const;
      //! \return the atomic mass of this atom given by the isotope
      //! (default of 0 gives the most abundant isotope)
      double     GetExactMass()     const;
      //! \return the internal atom index (e.g., inside an OBMol)
      unsigned int GetIdx()           const { return((int)_idx);  }
      unsigned int GetIndex() const { return _idx - 1; }
      unsigned long GetId() const { return _id; }
      //! \return the index into a pointer-driven array as used by
      //!   GetCoordPtr() or SetCoordPtr()
      unsigned int GetCoordinateIdx() const { return((int)_cidx); }
      //! \return The number of explicit bonds to this atom
      unsigned int GetExplicitDegree() const { return (unsigned int)_vbond.size(); }
      //! \return The total number of bonds to this atom including bonds to implicit hydrogens
      unsigned int GetTotalDegree() const { return (unsigned int)(_vbond.size() + _imph); }
      //! \return The sum of the bond orders of the explicit bonds to this atom
      unsigned int GetExplicitValence() const;
      //! \return The sum of the bond orders of all bonds to this atom including bonds to implicit hydrogens
      unsigned int GetTotalValence() const;
      //! \return The hybridization of this atom: 1 for sp, 2 for sp2, 3 for sp3, 4 for sq. planar, 5 for trig. bipy, 6 for octahedral
      unsigned int GetHyb()             const;
      //! \return The number of implicit hydrogens attached to this atom
      unsigned char GetImplicitHCount() const { return _imph; };
      //! \return The number of non-hydrogens connected to this atom
      unsigned int GetHvyDegree()      const;
      //! \return The number of heteroatoms connected to an atom
      unsigned int GetHeteroDegree()   const;
      //! \return the atomic type (e.g., for molecular mechanics)
      char        *GetType();

      //! \return the x coordinate
      double      GetX() const   {        return(x());    }
      //! \return the y coordinate
      double      GetY() const  {        return(y());    }
      //! \return the z coordinate
      double      GetZ() const  {        return(z());    }

      // These methods check to see if there is a coordinate pointer
      // or an internal vector (e.g., SetCoordPtr())
      //! \return the x coordinate
      double      x() const {
        if (_c)            return((*_c)[_cidx]);
        else               return _v.x();
      }
      //! \return the y coordinate
      double      y() const {
        if (_c)            return((*_c)[_cidx+1]);
        else               return _v.y();
      }
      //! \return the z coordinate
      double      z() const {
        if (_c)            return((*_c)[_cidx+2]);
        else               return _v.z();
      }
      //! \return the coordinates as a double* or NULL if none.
      //!
      //! See SetCoordPtr() for more. If no coordinate pointer is used
      //! (e.g., only vector3), NULL will be returned.
      double     *GetCoordinate(){
        if (_c)          return(&(*_c)[_cidx]);
        else             return nullptr;
      }
      //! \return the coordinates as a vector3 object
      vector3   &GetVector();
      //! \return the coordinates as a vector3 object
      const vector3   &GetVector() const;
      //! \return the partial charge of this atom, calculating a Gasteiger charge if needed
      double     GetPartialCharge();
      //! \return the residue which contains this atom, or NULL if none exists
      OBResidue *GetResidue();
      //! \return the molecule which contains this atom, or NULL if none exists
      OBMol     *GetParent()        {return((OBMol*)_parent);}
      //! Create a vector for a new bond from this atom, with length given by the supplied parameter
      //! \return success or failure
      bool       GetNewBondVector(vector3 &v,double length);
      //! \return the OBBond object between this atom and that supplied,
      //! or NULL if the two atoms are not bonded
      OBBond    *GetBond(OBAtom *);
      //@}

      //! \name Iterator methods
      //@{
      //! \return An iterator to the beginning of the bonds to this atom
      OBBondIterator BeginBonds()
        { return(_vbond.begin()); }
      //! \return An iterator to the end of the bonds to this atom
      OBBondIterator EndBonds()
        { return(_vbond.end());   }
      //! Set the iterator @p i to the beginning of the bonds
      //! \return The first bond to this atom (or NULL if none exist)
      OBBond *BeginBond(OBBondIterator &i);
      //! Increment the iterator @p i
      //! \return The next bond to this atom (or NULL if none exist)
      OBBond *NextBond(OBBondIterator &i);
      //! Set the iterator @p i to the beginning of the bonds
      //! \return The first neighboring atom (or NULL if none exist)
      OBAtom *BeginNbrAtom(OBBondIterator &i);
      //! Increment the iterator @p i
      //! \return The next neighboring atom (or NULL if none exist)
      OBAtom *NextNbrAtom(OBBondIterator &i);
      //@}

      //! \return the distance to the atom defined by OBMol::GetAtom()
      double GetDistance(int index);
      //! \return the distance to the supplied OBAtom
      double GetDistance(OBAtom*);
      //! \return the distance to the coordinates of the supplied vector3
      //! \since version 2.4
      double GetDistance(vector3* v);
      //! \return the angle defined by this atom -> b (vertex) -> c
      double GetAngle(int b, int c);
      //! \return the angle defined by this atom -> b (vertex) -> c
      double GetAngle(OBAtom *b, OBAtom *c);

      //! \name Addition of residue/bond info. for an atom
      //@{

      //! If no residue has been set for this atom, create a new one
      void NewResidue()
        {
          if (!_residue)
            _residue = new OBResidue;
        }
      //! Add (set) the residue for this atom
      void AddResidue(OBResidue *res) { SetResidue(res); }
      //! Delete any residue associated with this atom
      void DeleteResidue(){
        if (_residue) {
          delete _residue;
          _residue = nullptr; // Make sure to clear that a residue existed
        }
      }
      //! Add a bond to the internal list. Does not update the bond.
      void AddBond(OBBond *bond) { _vbond.push_back(bond); }
      //! \brief Insert @p bond into the internal list at the position from @p i
      //! Does not modify the bond
      void InsertBond(OBBondIterator &i, OBBond *bond)
        {
          _vbond.insert(i, bond);
        }
      //! Find @p bond and remove it from the internal list. Does not update the bond.
      bool DeleteBond(OBBond* bond);
      //! Clear all bonding information in this atom (does not delete them)
      void ClearBond() {_vbond.clear();}
      //@}

      //! \name Builder utilities
      //@{

      //! \brief If this is a hydrogen atom, transform into a methyl group
      //! \return success or failure
      bool HtoMethyl();
      //! Change the hybridization of this atom and modify the geometry accordingly
      //! \return success or failure
      bool SetHybAndGeom(int);
      //@}

      //! \name Property information
      //@{
      //! \return The number of oxygen atoms connected that only have one heavy valence
      unsigned int  CountFreeOxygens()      const;
      //! \return The number of sulfur atoms connected that only have one heavy valence
      //! \since version 2.4
      unsigned int  CountFreeSulfurs()      const;
      //! \return The number of hydrogens explicitly bound to this atom, optionally excluding D,T and isotope explicitly set to 1
      unsigned int  ExplicitHydrogenCount(bool ExcludeIsotopes=false) const;
      //! \return The number of rings that contain this atom
      unsigned int  MemberOfRingCount()     const;
      //! \return The size of the smallest ring that contains this atom (0 if not in a ring)
      unsigned int  MemberOfRingSize()	  const;
      //! \return The number of explicit ring connections to this atom
      unsigned int  CountRingBonds() const;
      //! \return The smallest angle of bonds to this atom
      double	  SmallestBondAngle();
      //! \return The average angle of bonds to this atom
      double	  AverageBondAngle();
      /** Lewis acid/base vacancies for this atom
       *  \return A pair of integers, where first is acid count and second is base count
       *  \since version 2.3
       */
      std::pair<int, int> LewisAcidBaseCounts() const;
      //! \return Is there any residue information?
      bool HasResidue()    { return(_residue != nullptr);    }
      //! \return Is this a HETATM in a residue (returns false if not in a residue)
      //! \since version 2.4
      bool IsHetAtom() {
        if (_residue == nullptr)
          return false;
        else
          return _residue->IsHetAtom(this);
      }
      //! \return Is the specified element, as specified by atom number (see OBElement namespace)?
      bool IsElement(const unsigned int e) const {
           return e == _ele;
      }
      //! \return Is the atom aromatic?
      bool IsAromatic()      const;
      //! \return Is the atom in a ring?
      bool IsInRing()        const;
      //! \return Is the atom in a ring of a given size?
      bool IsInRingSize(int) const;
      //! \return Is this atom an element in the 15th or 16th main groups
      //!  (i.e., N, O, P, S ...) ?
      bool IsHeteroatom();
      //! \return Is this atom directly connected to the supplied OBAtom?
      bool IsConnected(OBAtom*);
      //! \return Is this atom related to the supplied OBAtom in
      //!  a 1,3 bonding pattern?
      bool IsOneThree(OBAtom*);
      //! \return Is this atom related to the supplied OBAtom in
      //!  a 1,4 bonding pattern?
      bool IsOneFour(OBAtom*);
      //! \return Is this atom an oxygen in a carboxyl (-CO2 or CO2H) group?
      bool IsCarboxylOxygen();
      //! \return Is this atom an oxygen in a phosphate (R-PO3) group?
      bool IsPhosphateOxygen();
      //! \return Is this atom an oxygen in a sulfate (-SO3) group?
      bool IsSulfateOxygen();
      //! \return Is this atom an oxygen in a nitro (-NO2) group?
      bool IsNitroOxygen();
      //! \return Is this atom a nitrogen in an amide (-C(=O)NR2) group?
      bool IsAmideNitrogen();
      //! \return Is this atom a hydrogen connected to a polar atom
      //!  (i.e., N, O, P, S)
      bool IsPolarHydrogen();
      //! \return Is this atom a hydrogen connected to a non-polar atom
      //!  (i.e., C)
      bool IsNonPolarHydrogen();
      //! \return Is this atom an aromatic nitrogen with at least one
      //!  double bond to an oxygen atom
      bool IsAromaticNOxide();
      //! \return Is this atom chiral?
      bool IsChiral();
      //! \return Is the atom part of a periodic unit cell?
      bool IsPeriodic() const;
      //! \return Is this atom an axial atom in a ring
      bool IsAxial();
      //! \return Is this atom a hydrogen-bond acceptor  (considering also atom surrounding)
      bool IsHbondAcceptor();
      //! \return Is this atom a hydrogen-bond acceptor (old function)?
      bool IsHbondAcceptorSimple();
            //! \return Is this atom a hydrogen-bond donor?
      bool IsHbondDonor();
      //! \return Is this a hydrogen atom attached to a hydrogen-bond donor?
      bool IsHbondDonorH();
      //! \return Is this atom a metal?
      //! \since version 2.4
      bool IsMetal();
      //! \return Whether a neighboring atom (alpha) has an unsaturated bond
      //!   to a third atom (beta).
      //! \param includePandS Whether to include phosphorus and sulfur neighbors
      //! in this determination (or to exclude them)
      bool HasAlphaBetaUnsat(bool includePandS=true);
      //! \return Whether this atom is connected to any bond with order == @p bo
      bool HasBondOfOrder(unsigned int bo);
      //! \return The count of bonds connected to this atom with order == @p bo
      int  CountBondsOfOrder(unsigned int bo);
      //! \return The maximum bond order for this atom
      int  HighestBondOrder();
      //! \return Whether this atom is connected to any bond with order >1
      bool HasNonSingleBond();
      //! \return Does this atom have a single bond
      bool HasSingleBond()    {        return(HasBondOfOrder(1));    }
      //! \return Does this atom have a double bond
      bool HasDoubleBond()    {        return(HasBondOfOrder(2));    }
      //! \return Does this atom have an aromatic bond
      bool HasAromaticBond()  {        return(HasBondOfOrder(5));    }
      //! \return Whether this atom matches the first atom in a given SMARTS pattern
      bool MatchesSMARTS(const char *);
      //@}

    }; // class OBAtom

}// namespace OpenBabel

#endif   // OB_ATOM_H

//! \file atom.h
//! \brief Handle atoms

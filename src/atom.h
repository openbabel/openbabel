/**********************************************************************
atom.h - Handle OBAtom class.
 
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2006 by Geoffrey R. Hutchison
Some portions Copyright (C) 2003 by Michael Banck
 
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

#ifndef OB_ATOM_H
#define OB_ATOM_H

#include "babelconfig.h"

#ifndef EXTERN
#  define EXTERN extern
#endif

#include <vector>
#include <string>
#include <map>

#include "base.h"
#include "residue.h"

namespace OpenBabel
{

  class OBBond;
  class OBMol;

  //! OBNodeBase is declared for backwards-compatibility with 2.0 and earlier code
  typedef OBAtom OBNodeBase;

  //ATOM Property Macros (flags)
  //! Atom is in a 4-membered ring
#define OB_4RING_ATOM     (1<<1)
  //! Atom is in a 3-membered ring
#define OB_3RING_ATOM     (1<<2)
  //! Atom is aromatic
#define OB_AROMATIC_ATOM  (1<<3)
  //! Atom is in a ring
#define OB_RING_ATOM      (1<<4)
  //! Atom has clockwise SMILES chiral stereochemistry (i.e., "@@")
#define OB_CSTEREO_ATOM   (1<<5)
  //! Atom has anticlockwise SMILES chiral stereochemistry (i.e., "@")
#define OB_ACSTEREO_ATOM  (1<<6)
  //! Atom is an electron donor
#define OB_DONOR_ATOM     (1<<7)
  //! Atom is an electron acceptor
#define OB_ACCEPTOR_ATOM  (1<<8)
  //! Atom is chiral
#define OB_CHIRAL_ATOM    (1<<9)
  //! Atom has + chiral volume
#define OB_POS_CHIRAL_ATOM (1<<10)
  //! Atom has - chiral volume
#define OB_NEG_CHIRAL_ATOM (1<<11)
  //! Atom has no hydrogen attached. Temporary use only during SMILES input
#define OB_ATOM_HAS_NO_H   (1<<12)
  // 13-16 currently unused

  // Class OBAtom
  // class introduction in atom.cpp
 class OBAPI OBAtom: public OBBase
    {
    protected:
      char                          _ele;       //!< atomic number (type char to minimize space -- allows for 0..255 elements)
      char                          _impval;    //!< implicit valence
      char                          _type[6];   //!< atomic type
      short                         _fcharge;   //!< formal charge
      unsigned short                _isotope;   //!< isotope (0 = most abundant)
      short                         _spinmultiplicity;//!< atomic spin, e.g., 2 for radical  1 or 3 for carbene

      unsigned int                  _idx;       //!< unique node index (GetIdx(), SetIdx())
      OBMol                        *_parent;    //!< parent molecule (if any)
      std::vector<OBBond*>          _vbond;     //!< bonds to this atom -- assumed to be one of the endpoints

      unsigned short                _cidx;      //!< index into coordinate array
      unsigned short                _hyb;       //!< hybridization
      unsigned short                _flags;     //!< bitwise flags (e.g. aromaticity)
      double                        _pcharge;  //!< partial charge
      double                      **_c;        //!< coordinate array in double*
      vector3                       _v;         //!< coordinate vector
      OBResidue                    *_residue;   //!< parent residue (if applicable)

      int  GetFlag() const    {  return(_flags);  }
      //! Sets the bitwise @p flag
      void SetFlag(int flag)  { _flags |= flag;   }
      //! \return Returns true of the atom has the @p flag
      bool HasFlag(int flag)  {  return((_flags & flag) ? true : false); }

    public:
       //! Used internally by graph traversal algorithms
      bool Visit;

      //! Constructor
      OBAtom();
      //! Destructor
      virtual ~OBAtom();
      //! Assignment
      OBAtom &operator = (OBAtom &);
      //! Clear all data
      void Clear();

      //! \name Methods to set atomic information
      //@{
      //! Set atom index (i.e., in an OBMol)
      void SetIdx(int idx)    { _idx = idx; _cidx = (idx-1)*3; }
      //! Set atom hybridization (i.e., 1 = sp, 2 = sp2, 3 = sp3 ...)
      void SetHyb(int hyb)    { _hyb = hyb; }
      //! Set atomic number
      void SetAtomicNum(int atomicnum)    { _ele = (char)atomicnum; }
      //! Set isotope number (actual atomic weight is tabulated automatically, 0 = most abundant)
      void SetIsotope(unsigned int iso);
      //! Set the implicit valence to @p val
      void SetImplicitValence(int val)    { _impval = (char)val; }
      //! Increase the implicit valence by one
      void IncrementImplicitValence()     { _impval++; }
      //! Decrease the implicit valence by one
      void DecrementImplicitValence()     { _impval--; }
      //! Set the formal charge of the atom to @p fcharge
      void SetFormalCharge(int fcharge)   { _fcharge = fcharge; }
      //! Set the atomic spin to @p spin. See _spinmultiplicity
      void SetSpinMultiplicity(short spin){ _spinmultiplicity = spin; }
      void SetType(char *type);
      void SetType(std::string &type);
      //! Set the patrical charge to @p pcharge
      void SetPartialCharge(double pcharge){ _pcharge = pcharge; }
      void SetVector(vector3 &v);
      void SetVector(const double x,const double y,const double z);
      //! Set the position of this atom from a pointer-driven array of coordinates
      void SetCoordPtr(double **c)        { _c = c; _cidx = (GetIdx()-1)*3; }
      //! Set the position of this atom based on the internal pointer array (i.e. from SetCoordPtr() )
      void SetVector();
      void SetResidue(OBResidue *res)     { _residue=res; }
      void SetParent(OBMol *ptr)          { _parent=ptr; }
      void SetAromatic()                  { SetFlag(OB_AROMATIC_ATOM); }
      void UnsetAromatic()                { _flags &= (~(OB_AROMATIC_ATOM)); }
      //! Mark atom as having SMILES clockwise stereochemistry (i.e., "@@")
      void SetClockwiseStereo()           { SetFlag(OB_CSTEREO_ATOM|OB_CHIRAL_ATOM); }
      //! Mark atom as having SMILES anticlockwise stereochemistry (i.e., "@")
      void SetAntiClockwiseStereo()       { SetFlag(OB_ACSTEREO_ATOM|OB_CHIRAL_ATOM); }
      //! Mark an atom as having + chiral volume
      void SetPositiveStereo() { SetFlag(OB_POS_CHIRAL_ATOM|OB_CHIRAL_ATOM); }
      //! Mark an atom as having - chiral volume
      void SetNegativeStereo() { SetFlag(OB_NEG_CHIRAL_ATOM|OB_CHIRAL_ATOM); }
      //! Clear all stereochemistry information
      void UnsetStereo()
        {
          _flags &= ~(OB_ACSTEREO_ATOM);
          _flags &= ~(OB_CSTEREO_ATOM);
          _flags &= ~(OB_POS_CHIRAL_ATOM);
          _flags &= ~(OB_NEG_CHIRAL_ATOM);
          _flags &= ~(OB_CHIRAL_ATOM);
        }
      //! Mark an atom as belonging to at least one ring
      void SetInRing()         { SetFlag(OB_RING_ATOM); }
      //! Mark an atom as being chiral with unknown stereochemistry
      void SetChiral()         { SetFlag(OB_CHIRAL_ATOM); }
      //! Clear the internal coordinate pointer
      void ClearCoordPtr()     { _c = NULL; _cidx=0; }
      //@}

      //! \name Methods to retrieve atomic information
      //@{
      //int        GetStereo()        const { return((int)_stereo);}
      //! \return the formal charge for this atom
      int          GetFormalCharge()  const { return(_fcharge);    }
      //! \return the atomic number for this atom
      unsigned int GetAtomicNum()     const { return((unsigned int)_ele); }
      //! \return the isotope for this atom, if specified, or 0 for unspecified
      unsigned short int GetIsotope() const { return(_isotope);    }
      //! \return the atomic spin, e.g., 0 (default) for singlet,
      //!   2 for radical  1 or 3 for carbene
      int          GetSpinMultiplicity() const { return(_spinmultiplicity); }
      //! \return the atomic mass of this atom given by standard IUPAC
      //!  average molar mass
      double     GetAtomicMass()    const;
      //! \return the atomic mass of given by the isotope
      //! (default of 0 gives the most abundant isotope)
      double     GetExactMass()     const;
      //! \return the internal atom index (e.g., inside an OBMol)
      unsigned int GetIdx()           const { return((int)_idx);  }
      //! \return the index into a pointer-driven array as used by
      //!   GetCoordPtr() or SetCoordPtr()
      unsigned int GetCoordinateIdx() const { return((int)_cidx); }
      //! \deprecated Use GetCoordinateIdx() instead
      unsigned int GetCIdx()          const { return((int)_cidx); }
      //! \return The current number of explicit connections
      unsigned int GetValence()       const
        {
          return((_vbond.empty()) ? 0 : _vbond.size());
        }
      //! \return The hybridization of this atom (i.e. 1 for sp, 2 for sp2, 3 for sp3)
      unsigned int GetHyb()             const;
      //! \return The implicit valence of this atom type (i.e. maximum number of connections expected)
      unsigned int GetImplicitValence() const;
      //! \return The number of non-hydrogens connected to this atom
      unsigned int GetHvyValence()      const;
      //! \return The number of heteroatoms connected to an atom
      unsigned int GetHeteroValence()   const;
      //! \return the atomic type (e.g., for molecular mechanics)
      char        *GetType();

      //! \return the x coordinate
      double      GetX()    {        return(x());    }
      //! \return the x coordinate
      double      x()
        {
          if (_c)
            return((*_c)[_cidx]);
          else
            return _v.x();
        }
      //! \return the y coordinate
      double      GetY()    {        return(y());    }
      //! \return the y coordinate
      double      y()
        {
          if (_c)
            return((*_c)[_cidx+1]);
          else
            return _v.y();
        }
      //! \return the z coordinate
      double      GetZ()    {        return(z());    }
      //! \return the z coordinate
      double      z()
        {
          if (_c)
            return((*_c)[_cidx+2]);
          else
            return _v.z();
        }
      //! \return the coordinates as a double*
      double     *GetCoordinate()
        {
          if (_c)
            return(&(*_c)[_cidx]);
          else
            return NULL;
        }
      //! Return the coordinates as a vector3 object
      vector3   &GetVector();
      //! Return the partial charge of this atom, calculating a Gasteiger charge if needed
      double     GetPartialCharge();
      //! Return the residue which contains this atom, or NULL if none exists
      OBResidue *GetResidue();
      //! Return the molecule which contains this atom, or NULL if none exists
      OBMol     *GetParent()        {return((OBMol*)_parent);}
      //! Create a vector for a new bond from this atom, with length given by the supplied parameter
      //! \return success or failure
      bool       GetNewBondVector(vector3 &v,double length);
      //! Return the OBBond object between this atom and that supplied,
      //! or NULL if the two atoms are not bonded
      OBBond    *GetBond(OBAtom *);
      //! Return a pointer to the "next" atom (by atom index) in the
      //! parent OBMol, or NULL if no such atom exists.
      //! \deprecated Use any of the other iterator methods. This
      //! method will be removed in the future.
      OBAtom    *GetNextAtom();
      //@}

      //! \name Iterator methods
      //@{
      std::vector<OBBond*>::iterator BeginBonds()
        { return(_vbond.begin()); }
      std::vector<OBBond*>::iterator EndBonds()
        { return(_vbond.end());   }
      OBBond *BeginBond(std::vector<OBBond*>::iterator &i);
      OBBond *NextBond(std::vector<OBBond*>::iterator &i);
      OBAtom *BeginNbrAtom(std::vector<OBBond*>::iterator &);
      OBAtom *NextNbrAtom(std::vector<OBBond*>::iterator &);
      //@}

      //! \return the distance to the atom defined by OBMol::GetAtom()
      double GetDistance(int index);
      //! \return the distance to the supplied OBAtom
      double GetDistance(OBAtom*);
      //! \return the angle defined by this atom -> b (vertex) -> c
      double GetAngle(int b, int c);
      //! \return the angle defined by this atom -> b (vertex) -> c
      double GetAngle(OBAtom *b, OBAtom *c);

      //! \name Addition of residue/bond info. for an atom
      //@{
      void NewResidue()
        {
          if (!_residue)
            _residue = new OBResidue;
        }
      void DeleteResidue()
        {
          if (_residue)
            delete _residue;
        }
      void AddBond(OBBond *bond)
        {
          _vbond.push_back(bond);
        }
      void InsertBond(std::vector<OBBond*>::iterator &i, OBBond *bond)
        {
          _vbond.insert(i, bond);
        }
      bool DeleteBond(OBBond*);
      void ClearBond() {_vbond.clear();}
      //@}

      //! \name Requests for atomic property information
      //@{
      //! \return The number of oxygen atoms connected that only have one heavy valence
      unsigned int  CountFreeOxygens()      const;
      //! \return The number of hydrogens needed to fill the implicit valence of this atom
      unsigned int  ImplicitHydrogenCount() const;
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
      //! \return The sum of the bond orders of the bonds to the atom (i.e. double bond = 2...)
      unsigned int  BOSum()                 const;
      //! \return The sum of the bond orders of bonds to the atom, considering only KDouble, KTriple bonds
      unsigned int  KBOSum()                const;
      //@}

      //! \name Builder utilities
      //@{
      //! If this is a hydrogen atom, transform into a methyl group
      //! \return success or failure
      bool HtoMethyl();
      //! Change the hybridization of this atom and modify the geometry accordingly
      //! \return success or failure
      bool SetHybAndGeom(int);
      //! Mark that atom has no hydrogens attached
      void ForceNoH() {SetFlag(OB_ATOM_HAS_NO_H);}
      //! \return if atom has been marked as having
      //!  no hydrogens attached
      bool HasNoHForced() {return HasFlag(OB_ATOM_HAS_NO_H);}
      //@}

      //! \name Property information
      //@{
      //! \return Is there any residue information?
      bool HasResidue()    { return(_residue != NULL);    }
      //! \return Is the atom hydrogen?
      bool IsHydrogen()    { return(GetAtomicNum() == 1); }
      //! \return Is the atom carbon?
      bool IsCarbon()      { return(GetAtomicNum() == 6); }
      //! \return Is the atom nitrogen?
      bool IsNitrogen()    { return(GetAtomicNum() == 7); }
      //! \return Is the atom oxygen?
      bool IsOxygen()      { return(GetAtomicNum() == 8); }
      //! \return Is the atom sulfur?
      bool IsSulfur()      { return(GetAtomicNum() == 16);}
      //! \return Is the atom phosphorus?
      bool IsPhosphorus()  { return(GetAtomicNum() == 15);}
      //! \return Is the atom aromatic?
      bool IsAromatic()      const;
      //! \return Is the atom in a ring?
      bool IsInRing()        const;
      //! \return Is the atom in a ring of a given size?
      bool IsInRingSize(int) const;
      //! \return Is this atom an element in the 15th or 16th main groups
      //!  (i.e., N, O, P, S ...) ?
      bool IsHeteroatom();
      //! \return Is this atom any element except carbon or hydrogen?
      bool IsNotCorH();
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
      //! \return Is this atom an axial atom in a ring
      bool IsAxial();
      //! \return Does this atom have SMILES-specified clockwise "@@" stereochemistry?
      bool IsClockwise()         { return(HasFlag(OB_CSTEREO_ATOM));  }
      //! \return Does this atom have SMILES-specified anticlockwise "@" stereochemistry?
      bool IsAntiClockwise()     { return(HasFlag(OB_ACSTEREO_ATOM)); }
      //! \return Does this atom have a positive chiral volume?
      bool IsPositiveStereo() { return(HasFlag(OB_POS_CHIRAL_ATOM)); }
      //! \return Does this atom have a negative chiral volume?
      bool IsNegativeStereo() { return(HasFlag(OB_NEG_CHIRAL_ATOM)); }
      //! \return Does this atom have SMILES-specified stereochemistry?
      bool HasChiralitySpecified()
        { return(HasFlag(OB_CSTEREO_ATOM|OB_ACSTEREO_ATOM)); }
      //! \return Does this atom have a specified chiral volume?
      bool HasChiralVolume()
        { return(HasFlag(OB_POS_CHIRAL_ATOM|OB_NEG_CHIRAL_ATOM)); }
      //! \return Is this atom a hydrogen-bond acceptor (receptor)?
      bool IsHbondAcceptor();
      //! \return Is this atom a hydrogen-bond donor?
      bool IsHbondDonor();
      //! \return Is this a hydrogen atom attached to a hydrogen-bond donor?
      bool IsHbondDonorH();
      bool HasAlphaBetaUnsat(bool includePandS=true);
      bool HasBondOfOrder(unsigned int);
      int  CountBondsOfOrder(unsigned int);
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

  //! A standard iterator over a vector of atoms
  typedef std::vector<OBAtom*>::iterator OBAtomIterator;

}// namespace OpenBabel

#endif   // OB_ATOM_H

//! \file atom.h
//! \brief Handle atoms

/**********************************************************************
mol.h - Handle molecules. Declarations of OBMol, OBAtom, OBBond, OBResidue.
        (the main header for Open Babel)
 
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2005 by Geoffrey R. Hutchison
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

#ifndef OB_MOL_H
#define OB_MOL_H

#include "babelconfig.h"

#ifndef EXTERN
#  define EXTERN extern
#endif

#include <math.h>

#include <algorithm>
#include <vector>
#include <string>
#include <map>

#if HAVE_IOSTREAM
#include <iostream>
#elif HAVE_IOSTREAM_H
#include <iostream.h>
#endif

#if HAVE_FSTREAM
#include <fstream>
#elif HAVE_FSTREAM_H
#include <fstream.h>
#endif

#include "base.h"
#include "data.h"
#include "chains.h"
#include "math/vector3.h"
#include "bitvec.h"
#include "ring.h"
#include "generic.h"
#include "typer.h"
#include "oberror.h"
#include "obiter.h"

namespace OpenBabel
{

  class OBAtom;
  class OBBond;
  class OBMol;
  class OBInternalCoord;

  // Class OBResidue
  // class introduction in residue.cpp
  class OBAPI OBResidue : public OBBase
    {
    public:

      //! Constructor
      OBResidue(void);
      //! Copy constructor
      OBResidue(const OBResidue &);
      //! Destructor
      virtual ~OBResidue(void);

      OBResidue &operator=(const OBResidue &);

      void    AddAtom(OBAtom *atom);
      void    InsertAtom(OBAtom *atom);
      void    RemoveAtom(OBAtom *atom);
      void    Clear(void);

      void    SetName(const std::string &resname);
      void    SetNum(const unsigned int resnum);
      void    SetChain(const char chain);
      void    SetChainNum(const unsigned int chainnum);
      void    SetIdx(const unsigned int idx);

      void    SetAtomID(OBAtom *atom, const std::string &id);
      void    SetHetAtom(OBAtom *atom, bool hetatm);
      //! Set the atomic serial number for a given atom (see OBSerialNums)
      void    SetSerialNum(OBAtom *atom, unsigned int sernum);

      std::string    GetName(void)                  const;
      unsigned int   GetNum(void)                   const;
      unsigned int   GetNumAtoms()                  const;
      char           GetChain(void)                 const;
      unsigned int   GetChainNum(void)              const;
      unsigned int   GetIdx(void)                   const;
      unsigned int   GetResKey(void)                const;

      std::vector<OBAtom*> GetAtoms(void)           const;
      std::vector<OBBond*> GetBonds(bool exterior= true)const;

      std::string    GetAtomID(OBAtom *atom)        const;
      //! \return the serial number of the supplied atom (uses OBSerialNums)
      unsigned       GetSerialNum(OBAtom *atom)     const;

      bool           GetAminoAcidProperty(int)      const;
      bool           GetAtomProperty(OBAtom *, int) const;
      bool           GetResidueProperty(int)        const;

      bool           IsHetAtom(OBAtom *atom)        const;
      bool           IsResidueType(int)             const;

      //! \deprecated Use FOR_ATOMS_OF_RESIDUE and OBResidueAtomIter instead
      OBAtom *BeginAtom(std::vector<OBAtom*>::iterator &i);
      //! \deprecated Use FOR_ATOMS_OF_RESIDUE and OBResidueAtomIter instead
      OBAtom *NextAtom(std::vector<OBAtom*>::iterator &i);

    protected: // members

      unsigned int              _idx;   //!< Residue index (i.e., internal index in an OBMol)
      char                              _chain; //!< Chain ID
      unsigned int              _aakey; //!< Amino Acid key ID -- see SetResidueKeys()
      unsigned int              _reskey;//!< Residue key ID -- see SetResidueKeys()
      unsigned int              _resnum;//!< Residue number (i.e., in file)
      std::string                 _resname;//!< Residue text name

      std::vector<bool>           _hetatm;//!< Is a given atom a HETAM
      std::vector<std::string>    _atomid;//!< Residue atom text IDs
      std::vector<OBAtom*>        _atoms; //!< List of OBAtom in this residue
      std::vector<unsigned int>   _sernum;//!< List of serial numbers
      //    std::vector<OBGenericData*> _vdata; //!< Custom data
    }; // OBResidue


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
  class OBAPI OBAtom : public OBNodeBase
    {
    protected:
      char                          _ele;       //!< atomic number (type char to minimize space -- allows for 0..255 elements)
      char                          _impval;    //!< implicit valence
      char                          _type[6];   //!< atomic type
      short                         _fcharge;   //!< formal charge
      unsigned short                _isotope;   //!< isotope (0 = most abundant)
      short                           _spinmultiplicity;//!< atomic spin, e.g., 2 for radical  1 or 3 for carbene

      //unsigned short int          _idx;       //!< index in parent (inherited)
      unsigned short            _cidx;          //!< index into coordinate array
      unsigned short                _hyb;       //!< hybridization
      unsigned short                _flags;     //!< bitwise flags (e.g. aromaticity)
      double                         _pcharge;  //!< partial charge
      double                       **_c;        //!< coordinate array in double*
      vector3                       _v;         //!< coordinate vector
      OBResidue                    *_residue;   //!< parent residue (if applicable)
      //OBMol                      *_parent;    //!< parent molecule (inherited)
      //vector<OBBond*>             _bond;      //!< connections (inherited)
      //   std::vector<OBGenericData*>   _vdata;//!< custom data

      int  GetFlag() const    {  return(_flags);  }
      void SetFlag(int flag)  { _flags |= flag;   }
      bool HasFlag(int flag)  {  return((_flags & flag) ? true : false); }

    public:

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
      void SetImplicitValence(int val)    { _impval = (char)val; }
      void IncrementImplicitValence()     { _impval++; }
      void DecrementImplicitValence()     { _impval--; }
      void SetFormalCharge(int fcharge)   { _fcharge = fcharge; }
      void SetSpinMultiplicity(short spin){ _spinmultiplicity = spin; }
      void SetType(char *type);
      void SetType(std::string &type);
      void SetPartialCharge(double pcharge){ _pcharge = pcharge; }
      void SetVector(vector3 &v);
      void SetVector(const double x,const double y,const double z);
      //! Set the position of this atom from a pointer-driven array of coordinates
      void SetCoordPtr(double **c)        { _c = c; _cidx = (GetIdx()-1)*3; }
      //! Set the position of this atom based on the internal pointer array (i.e. from SetCoordPtr() )
      void SetVector();
      void SetResidue(OBResidue *res)     { _residue=res; }
      //  void SetParent(OBMol *ptr)      { _parent=ptr; } // inherited
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
      int          GetFormalCharge()  const { return(_fcharge);    }
      unsigned int GetAtomicNum()     const { return((unsigned int)_ele); }
      unsigned short int GetIsotope() const { return(_isotope);    }
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
      //! The current number of explicit connections
      unsigned int GetValence()       const
        {
          return((_vbond.empty()) ? 0 : _vbond.size());
        }
      //! The hybridization of this atom (i.e. 1 for sp, 2 for sp2, 3 for sp3)
      unsigned int GetHyb()             const;
      //! The implicit valence of this atom type (i.e. maximum number of connections expected)
      unsigned int GetImplicitValence() const;
      //! The number of non-hydrogens connected to this atom
      unsigned int GetHvyValence()      const;
      //! The number of heteroatoms connected to an atom
      unsigned int GetHeteroValence()   const;
      char        *GetType();

      //! The x coordinate
      double      GetX()    {        return(x());    }
      //! The y coordinate
      double      GetY()    {        return(y());    }
      //! The z coordinate
      double      GetZ()    {        return(z());    }
      double      x()
        {
          if (_c)
            return((*_c)[_cidx]);
          else
            return _v.x();
        }
      double      y()
        {
          if (_c)
            return((*_c)[_cidx+1]);
          else
            return _v.y();
        }
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
      //! \return the coordinates as a vector3 object
      vector3   &GetVector();
      //! \return the partial charge of this atom, calculating a Gasteiger charge if needed
      double     GetPartialCharge();
      OBResidue *GetResidue();
      //OBMol   *GetParent()        {return((OBMol*)_parent);}
      //! Create a vector for a new bond from this atom, with length given by the supplied parameter
      bool       GetNewBondVector(vector3 &v,double length);
      OBBond    *GetBond(OBAtom *);
      OBAtom    *GetNextAtom();
      //@}

      //! \name Iterator methods
      //@{
      //! \deprecated Use FOR_BONDS_OF_ATOM and OBAtomBondIter instead
      std::vector<OBEdgeBase*>::iterator BeginBonds()
        { return(_vbond.begin()); }
      //! \deprecated Use FOR_BONDS_OF_ATOM and OBAtomBondIter instead
      std::vector<OBEdgeBase*>::iterator EndBonds()
        { return(_vbond.end());   }
      //! \deprecated Use FOR_BONDS_OF_ATOM and OBAtomBondIter instead
      OBBond *BeginBond(std::vector<OBEdgeBase*>::iterator &i);
      //! \deprecated Use FOR_BONDS_OF_ATOM and OBAtomBondIter instead
      OBBond *NextBond(std::vector<OBEdgeBase*>::iterator &i);
      //! \deprecated Use FOR_NBORS_OF_ATOM and OBAtomAtomIter instead
      OBAtom *BeginNbrAtom(std::vector<OBEdgeBase*>::iterator &);
      //! \deprecated Use FOR_NBORS_OF_ATOM and OBAtomAtomIter instead
      OBAtom *NextNbrAtom(std::vector<OBEdgeBase*>::iterator &);
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
          _vbond.push_back((OBEdgeBase*)bond);
        }
      void InsertBond(std::vector<OBEdgeBase*>::iterator &i, OBBond *bond)
        {
          _vbond.insert(i, (OBEdgeBase*)bond);
        }
      bool DeleteBond(OBBond*);
      void ClearBond() {_vbond.clear();}
      //@}

      //! \name Requests for atomic property information
      //@{
      //! The number of oxygen atoms connected that only have one heavy valence
      unsigned int  CountFreeOxygens()      const;
      //! The number of hydrogens needed to fill the implicit valence of this atom
      unsigned int  ImplicitHydrogenCount() const;
      //! The number of hydrogens explicitly bound to this atom, optionally excluding D,T and isotope explicitly set to 1
      unsigned int  ExplicitHydrogenCount(bool ExcludeIsotopes=false) const;
      //! The number of rings that contain this atom
      unsigned int  MemberOfRingCount()     const;
      //! The size of the smallest ring that contains this atom (0 if not in a ring)
      unsigned int  MemberOfRingSize()	  const;
      //! \return The number of explicit ring connections to this atom
      unsigned int  CountRingBonds() const;
      //! The smallest angle of bonds to this atom
      double	  SmallestBondAngle();
      //! The average angle of bonds to this atom
      double	  AverageBondAngle();
      //! The sum of the bond orders of the bonds to the atom (i.e. double bond = 2...)
      unsigned int  BOSum()                 const;
      //! The sum of the bond orders of bonds to the atom, considering only KDouble, KTriple bonds
      unsigned int  KBOSum()                const;
      //@}

      //! \name Builder utilities
      //@{
      //! If this is a hydrogen atom, transform into a methyl group
      bool HtoMethyl();
      //! Change the hybridization of this atom and modify the geometry accordingly
      bool SetHybAndGeom(int);
      //! Mark that atom has no hydrogens attached
      void ForceNoH() {SetFlag(OB_ATOM_HAS_NO_H);}
      //! \return Return true if atom has been marked as having
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
      bool IsSulfur()      { return(GetAtomicNum() == 16);}
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
      //! \return Is this atom connected to the supplied OBAtom?
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
      bool IsAromaticNOxide();
      //! Is this atom chiral?
      bool IsChiral();
      bool IsAxial();
      //! Does this atom have SMILES-specified clockwise "@@" stereochemistry?
      bool IsClockwise()         { return(HasFlag(OB_CSTEREO_ATOM));  }
      //! Does this atom have SMILES-specified anticlockwise "@" stereochemistry?
      bool IsAntiClockwise()     { return(HasFlag(OB_ACSTEREO_ATOM)); }
      //! Does this atom have a positive chiral volume?
      bool IsPositiveStereo() { return(HasFlag(OB_POS_CHIRAL_ATOM)); }
      //! Does this atom have a negative chiral volume?
      bool IsNegativeStereo() { return(HasFlag(OB_NEG_CHIRAL_ATOM)); }
      //! Does this atom have SMILES-specified stereochemistry?
      bool HasChiralitySpecified()
        { return(HasFlag(OB_CSTEREO_ATOM|OB_ACSTEREO_ATOM)); }
      //! Does this atom have a specified chiral volume?
      bool HasChiralVolume()
        { return(HasFlag(OB_POS_CHIRAL_ATOM|OB_NEG_CHIRAL_ATOM)); }
      //! Is this atom a hydrogen-bond acceptor (receptor)?
      bool IsHbondAcceptor();
      //! Is this atom a hydrogen-bond donor?
      bool IsHbondDonor();
      //! Is this a hydrogen atom attached to a hydrogen-bond donor?
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


  // Class OBBond

  //BOND Property Macros (flags)
  //! An aromatic bond (regardless of bond order)
#define OB_AROMATIC_BOND  (1<<1)
  //! A solid black wedge in 2D representations -- i.e., "up" from the 2D plane
#define OB_WEDGE_BOND     (1<<2)
  //! A dashed "hash" bond in 2D representations -- i.e., "down" from the 2D plane
#define OB_HASH_BOND      (1<<3)
  //! A bond in a ring
#define OB_RING_BOND      (1<<4)
  //! The "upper" bond in a double bond cis/trans isomer (i.e., "/" in SMILES)
#define OB_TORUP_BOND     (1<<5)
  //! The "down" bond in a double bond cis/trans isomer (i.e., "\" in SMILES)
#define OB_TORDOWN_BOND   (1<<6)
  //! A Kekule single bond
#define OB_KSINGLE_BOND   (1<<7)
  //! A Kekule double bond
#define OB_KDOUBLE_BOND   (1<<8)
  //! A Kekule triple bond
#define OB_KTRIPLE_BOND   (1<<9)
#define OB_CLOSURE_BOND   (1<<10)
  // 11-16 currently unused

  // class introduction in bond.cpp
  class OBAPI OBBond : public OBEdgeBase
    {
    protected:
      char                          _order; //!< Bond order (1, 2, 3, 5=aromatic)
      unsigned short int            _flags; //!< Any flags for this bond
      //OBAtom                     *_bgn;   //!< Not needed, inherited from OBEdgeBase
      //OBAtom                     *_end;   //!< Not needed, inherited from OBEdgeBase
      //OBMol                      *_parent;//!< Not needed, inherited from OBEdgeBase
      //unsigned short int          _idx;   //!< Not needed, inherited from OBEdgeBase
      //    std::vector<OBGenericData*>   _vdata; //!< Generic data for custom information

      bool HasFlag(int flag)    { return((_flags & flag) != 0); }
      void SetFlag(int flag)    { _flags |= flag;               }
      void UnsetFlag(int flag)  { _flags &= (~(flag));          }

    public:
      //! Constructor
      OBBond();
      //! Destructor
      virtual ~OBBond();

      //! \name Bond modification methods
      //@{
      void SetIdx(int idx)
        {
          _idx = idx;
        }
      void SetBO(int order);
      void SetBegin(OBAtom *begin)
        {
          _bgn = begin;
        }
      void SetEnd(OBAtom *end)
        {
          _end = end;
        }
      // void SetParent(OBMol *ptr)               {_parent=ptr;} // (inherited)
      void SetLength(OBAtom*,double);
      void Set(int,OBAtom*,OBAtom*,int,int);
      void SetKSingle();
      void SetKDouble();
      void SetKTriple();
      void SetAromatic()    { SetFlag(OB_AROMATIC_BOND); }
      void SetHash()        { SetFlag(OB_HASH_BOND);     }
      void SetWedge()       { SetFlag(OB_WEDGE_BOND);    }
      void SetUp()          { SetFlag(OB_TORUP_BOND);   UnsetFlag(OB_TORDOWN_BOND); }
      void SetDown()        { SetFlag(OB_TORDOWN_BOND); UnsetFlag(OB_TORUP_BOND);   }
      void SetInRing()      { SetFlag(OB_RING_BOND);     }
      void SetClosure()     { SetFlag(OB_CLOSURE_BOND);  }

      void UnsetHash()      { UnsetFlag(OB_HASH_BOND);    }
      void UnsetWedge()     { UnsetFlag(OB_WEDGE_BOND);   }
      void UnsetUp()        { UnsetFlag(OB_TORUP_BOND);   }
      void UnsetDown()      { UnsetFlag(OB_TORDOWN_BOND); }
      void UnsetAromatic()  { UnsetFlag(OB_AROMATIC_BOND);}
      void UnsetKekule()
        {
          _flags &= (~(OB_KSINGLE_BOND|OB_KDOUBLE_BOND|OB_KTRIPLE_BOND));
        }
      //@}

      //! \name bond data request methods
      //@{
      unsigned int     GetBO()            const { return((int)_order); }
      unsigned int     GetBondOrder()     const { return((int)_order); }
      unsigned int     GetFlags()         const { return(_flags);      }
      unsigned int     GetBeginAtomIdx()  const { return(_bgn->GetIdx()); }
      unsigned int     GetEndAtomIdx()    const { return(_end->GetIdx()); }
      OBAtom *GetBeginAtom()    { return((OBAtom*)_bgn);    }
      OBAtom *GetEndAtom()      { return((OBAtom*)_end);    }
      OBAtom *GetNbrAtom(OBAtom *ptr)
        {
          return((ptr != _bgn)? (OBAtom*)_bgn : (OBAtom*)_end);
        }
      // OBMol  *GetParent()                 {return(_parent);}  // (inherited)
      double   GetEquibLength();
      double   GetLength();
      int     GetNbrAtomIdx(OBAtom *ptr)
        {
          return((ptr!=_bgn)?_bgn->GetIdx():_end->GetIdx());
        }
      //@}

      //! \name property request methods
      //@{
      //! \return Is the bond aromatic? 
      //!  (Note that the two atoms of the bond may be aromatic, 
      //!   but not the bond)
      bool IsAromatic() const;
      //! \return Is the bond part of a ring?
      bool IsInRing() const;
      //! Is the bond a rotatable bond?
      //!  Currently, this function classifies any bond with at least one heavy
      //!  atom, no sp-hybrid atoms (e.g., a triple bond somewhere) not in a ring
      //!  as a potential rotor. No other bond typing is attempted.
      bool IsRotor();
      bool IsAmide();
      bool IsPrimaryAmide();
      bool IsSecondaryAmide();
      bool IsEster();
      bool IsCarbonyl();
      bool IsSingle();
      bool IsDouble();
      bool IsTriple();
      bool IsKSingle();
      bool IsKDouble();
      bool IsKTriple();
      bool IsClosure();
      //! \return whether this is the "upper" bond in a double bond cis/trans
      //!   isomer (i.e., "/" in SMILES)
      bool IsUp()    {    return(HasFlag(OB_TORUP_BOND));    }
      //! \return whether this is the "lower" bond in a double bond cis/trans
      //!   isomer (i.e., "\" in SMILES)
      bool IsDown()  {    return(HasFlag(OB_TORDOWN_BOND));  }
      bool IsWedge() {    return(HasFlag(OB_WEDGE_BOND));    }
      bool IsHash()  {    return(HasFlag(OB_HASH_BOND));     }
      //! \return whether the geometry around this bond looks unsaturated
      bool IsDoubleBondGeometry();
      //@}

    }; // class OBBond


  // Class OBMol

  //MOL Property Macros (flags) -- 32+ bits
#define OB_SSSR_MOL              (1<<1)
#define OB_RINGFLAGS_MOL         (1<<2)
#define OB_AROMATIC_MOL          (1<<3)
#define OB_ATOMTYPES_MOL         (1<<4)
#define OB_CHIRALITY_MOL         (1<<5)
#define OB_PCHARGE_MOL           (1<<6)
#define OB_HYBRID_MOL            (1<<8)
#define OB_IMPVAL_MOL            (1<<9)
#define OB_KEKULE_MOL            (1<<10)
#define OB_CLOSURE_MOL           (1<<11)
#define OB_H_ADDED_MOL           (1<<12)
#define OB_PH_CORRECTED_MOL      (1<<13)
#define OB_AROM_CORRECTED_MOL    (1<<14)
#define OB_CHAINS_MOL            (1<<15)
#define OB_TCHARGE_MOL		 (1<<16)
#define OB_TSPIN_MOL             (1<<17)
  // flags 18-32 unspecified
#define OB_CURRENT_CONFORMER	 -1

  // class introduction in mol.cpp
  class OBAPI OBMol : public OBGraphBase
    {
    protected:
      int                           _flags;	//!< bitfield of flags
      bool                          _autoPartialCharge; //!< Assign partial charges automatically
      bool                          _autoFormalCharge; //!< Assign formal charges automatically
      std::string                   _title;	//!< Molecule title
      //vector<OBAtom*>             _atom;	//!< not needed (inherited)
      //vector<OBBond*>             _bond;	//!< not needed (inherited)
      unsigned short int            _dimension;   //!< Dimensionality of coordinates
      double                        _energy;      //!< Molecular heat of formation (if applicable)
      int				  _totalCharge; //!< Total charge on the molecule
      unsigned int                  _totalSpin;   //!< Total spin on the molecule (if not specified, assumes lowest possible spin)
      double                       *_c;	        //!< coordinate array
      std::vector<double*>          _vconf;       //!< vector of conformers
      unsigned short int            _natoms;      //!< Number of atoms
      unsigned short int            _nbonds;      //!< Number of bonds
      std::vector<OBResidue*>       _residue;     //!< Residue information (if applicable)
      std::vector<OBInternalCoord*> _internals;   //!< Internal Coordinates (if applicable)
      //    std::vector<OBGenericData*>   _vdata;       //!< Custom data -- see OBGenericData class for more
      unsigned short int            _mod;	        //!< Number of nested calls to BeginModify()

      bool  HasFlag(int flag)    { return((_flags & flag) ? true : false); }
      void  SetFlag(int flag)    { _flags |= flag; }

      //! \name Internal Kekulization routines -- see kekulize.cpp and NewPerceiveKekuleBonds()
      //@{
      void start_kekulize(std::vector <OBAtom*> &cycle, std::vector<int> &electron);
      int expand_kekulize(OBAtom *atom1, OBAtom *atom2, std::vector<int> &currentState, std::vector<int> &initState, std::vector<int> &bcurrentState, std::vector<int> &binitState, std::vector<bool> &mark);
      int getorden(OBAtom *atom);
      void expandcycle(OBAtom *atom, OBBitVec &avisit);
      //@}

    public:

      //! \name Initialization and data (re)size methods
      //@{
      //! Constructor
      OBMol();
      //! Copy constructor, copies atoms,bonds and OBGenericData
      OBMol(const OBMol &);
      //! Destructor
      virtual ~OBMol();
      //! Assignment, copies atoms,bonds and OBGenericData
      OBMol &operator=(const OBMol &mol);      
      //! Copies atoms and bonds but not OBGenericData
      OBMol &operator+=(const OBMol &mol);

      void ReserveAtoms(int natoms)
        {
          if (natoms && _mod)
            _vatom.reserve(natoms);
        }
      virtual OBAtom *CreateAtom(void);
      virtual OBBond *CreateBond(void);
      virtual void DestroyAtom(OBNodeBase*);
      virtual void DestroyBond(OBEdgeBase*);
      bool AddAtom(OBAtom&);
      bool AddBond(int,int,int,int flags=0,int insertpos=-1);
      bool AddBond(OBBond&);
      bool AddResidue(OBResidue&);
      bool InsertAtom(OBAtom &);
      bool DeleteAtom(OBAtom*);
      bool DeleteBond(OBBond*);
      bool DeleteResidue(OBResidue*);
      OBAtom    *NewAtom();
      OBResidue *NewResidue();
      //@}

      //! \name Molecule modification methods
      //@{
      //! Call when making many modifications -- clears conformer/rotomer data.
      virtual void BeginModify(void);
      //! Call when done with modificaions -- re-perceive data as needed.
      virtual void EndModify(bool nukePerceivedData=true);
      int GetMod()
        {
          return(_mod);
        }
      void IncrementMod()
        {
          _mod++;
        }
      void DecrementMod()
        {
          _mod--;
        }
      //@}

      //! \name Data retrieval methods
      //@{
      int          GetFlags()               { return(_flags); }
      //! \return the title of this molecule (often the filename)
      const char  *GetTitle() const         { return(_title.c_str()); }
      //! \return the number of atoms (i.e. OBAtom children)
      unsigned int NumAtoms() const         {  return(_natoms); }
      //! \return the number of bonds (i.e. OBBond children)
      unsigned int NumBonds() const         {  return(_nbonds); }
      //! \return the number of non-hydrogen atoms
      unsigned int NumHvyAtoms();
      //! \return the number of residues (i.e. OBResidue substituents)
      unsigned int NumResidues() const      { return(_residue.size()); }
      //! \return the number of rotatble bonds. See OBBond::IsRotor() for details
      unsigned int NumRotors();
    
      OBAtom      *GetAtom(int);
      OBAtom      *GetFirstAtom();
      OBBond      *GetBond(int);
      OBBond      *GetBond(int, int);
      OBBond      *GetBond(OBAtom*,OBAtom*);
      OBResidue   *GetResidue(int);
      std::vector<OBInternalCoord*> GetInternalCoord();
      //! \return the dihedral angle between the four atoms supplied a1-a2-a3-a4)
      double       GetTorsion(int,int,int,int);
      //! \return the dihedral angle between the four atoms supplied a1-a2-a3-a4)
      double       GetTorsion(OBAtom*,OBAtom*,OBAtom*,OBAtom*);
      //! \return the stochoimetric formula (e.g., C4H6O)
      std::string  GetFormula();
      //! \return the stochoimetric formula in spaced format e.g. C 4 H 6 O 1
      std::string OBMol::GetSpacedFormula(int ones=0, const char* sp=" ");
      //! \return the heat of formation for this molecule (in kcal/mol)
      double       GetEnergy() const { return(_energy); }
      //! \return the standard molar mass given by IUPAC atomic masses (amu)
      double       GetMolWt();
      //! \return the mass given by isotopes (or most abundant isotope, if not specified)
      double	 GetExactMass();
      //! \return the total charge on this molecule (i.e., 0 = neutral, +1, -1...)
      int		 GetTotalCharge();
      //! \return the total spin on this molecule (i.e., 1 = singlet, 2 = doublet...)
      unsigned int GetTotalSpinMultiplicity();
      //! \return the dimensionality of coordinates (i.e., 0 = unknown or no coord, 2=2D, 3=3D)
      unsigned short int GetDimension() const { return _dimension; }
      double      *GetCoordinates() { return(_c); }
      //! \return the Smallest Set of Smallest Rings has been run (see OBRing class
      std::vector<OBRing*> &GetSSSR();
      //! Get the current flag for whether formal charges are set with pH correction
      bool AutomaticFormalCharge()   { return(_autoFormalCharge);  }
      //! Get the current flag for whether partial charges are auto-determined
      bool AutomaticPartialCharge()  { return(_autoPartialCharge); }
      //@}


      //! \name Data modification methods
      //@{
      void   SetTitle(const char *title);
      void   SetTitle(std::string &title);
      //! Set the stochiometric formula for this molecule
      void   SetFormula(std::string molFormula);
      //! Set the heat of formation for this molecule (in kcal/mol)
      void   SetEnergy(double energy) { _energy = energy; }
      //! Set the dimension of this molecule (i.e., 0, 1 , 2, 3)
      void   SetDimension(unsigned short int d) { _dimension = d; }
      void   SetTotalCharge(int charge);
      void   SetTotalSpinMultiplicity(unsigned int spin);
      void   SetInternalCoord(std::vector<OBInternalCoord*> int_coord)
        { _internals = int_coord; }
      //! Set the flag for determining automatic formal charges with pH (default=true)
      void SetAutomaticFormalCharge(bool val)
        { _autoFormalCharge=val;  }
      //! Set the flag for determining partial charges automatically (default=true)
      void SetAutomaticPartialCharge(bool val)
        { _autoPartialCharge=val; }

      //! Mark that aromaticity has been perceived for this molecule (see OBAromaticTyper)
      void   SetAromaticPerceived()    { SetFlag(OB_AROMATIC_MOL);    }
      //! Mark that Smallest Set of Smallest Rings has been run (see OBRing class)
      void   SetSSSRPerceived()        { SetFlag(OB_SSSR_MOL);        }
      //! Mark that rings have been perceived (see OBRing class for details)
      void   SetRingAtomsAndBondsPerceived(){SetFlag(OB_RINGFLAGS_MOL);}
      //! Mark that atom types have been perceived (see OBAtomTyper for details)
      void   SetAtomTypesPerceived()   { SetFlag(OB_ATOMTYPES_MOL);   }
      //! Mark that chains and residues have been perceived (see OBChainsParser)
      void   SetChainsPerceived()      { SetFlag(OB_CHAINS_MOL);      }
      //! Mark that chirality has been perceived
      void   SetChiralityPerceived()   { SetFlag(OB_CHIRALITY_MOL);   }
      //! Mark that partial charges have been assigned
      void   SetPartialChargesPerceived(){ SetFlag(OB_PCHARGE_MOL);   }
      void   SetHybridizationPerceived() { SetFlag(OB_HYBRID_MOL);    }
      void   SetImplicitValencePerceived(){ SetFlag(OB_IMPVAL_MOL);   }
      void   SetKekulePerceived()      { SetFlag(OB_KEKULE_MOL);      }
      void   SetClosureBondsPerceived(){ SetFlag(OB_CLOSURE_MOL);     }
      void   SetHydrogensAdded()       { SetFlag(OB_H_ADDED_MOL);     }
      void   SetCorrectedForPH()       { SetFlag(OB_PH_CORRECTED_MOL);}
      void   SetAromaticCorrected()    { SetFlag(OB_AROM_CORRECTED_MOL);}
      void   SetSpinMultiplicityAssigned(){ SetFlag(OB_TSPIN_MOL);    }
      void   SetFlags(int flags)       { _flags = flags;              }

      void   UnsetAromaticPerceived()  { _flags &= (~(OB_AROMATIC_MOL));   }
      void   UnsetPartialChargesPerceived(){ _flags &= (~(OB_PCHARGE_MOL));}
      void   UnsetImplicitValencePerceived(){_flags &= (~(OB_IMPVAL_MOL)); }
      void   UnsetFlag(int flag)       { _flags &= (~(flag));              }

      //! \name Molecule modification methods
      //@{
      // Description in transform.cpp
      virtual OBBase*    DoTransformations(const std::map<std::string,std::string>* pOptions);
      static const char* ClassDescription();
      //! Clear all information from a molecule
      bool Clear();
      //! Renumber the atoms of this molecule according to the order in the supplied vector
      void RenumberAtoms(std::vector<OBNodeBase*>&);
      //! Translate one conformer and rotate by a rotation matrix (which is returned) to the inertial frame-of-reference
      void ToInertialFrame(int conf, double *rmat);
      //! Translate all conformers to the inertial frame-of-reference
      void ToInertialFrame();
      //! Translates all conformers in the molecule by the supplied vector
      void Translate(const vector3 &v);
      //! Translates one conformer in the molecule by the supplied vector
      void Translate(const vector3 &v, int conf);
      void Rotate(const double u[3][3]);
      void Rotate(const double m[9]);
      void Rotate(const double m[9],int nconf);
      //! Translate to the center of all coordinates (for this conformer)
      void Center();
      //! Transform to standard Kekule bond structure (presumably from an aromatic form)
      bool Kekulize();
      bool PerceiveKekuleBonds();

      void NewPerceiveKekuleBonds();

      bool DeleteHydrogen(OBAtom*);
      bool DeleteHydrogens();
      bool DeleteHydrogens(OBAtom*);
      bool DeleteNonPolarHydrogens();
      bool AddHydrogens(bool polaronly=false,bool correctForPH=true);
      bool AddHydrogens(OBAtom*);
      bool AddPolarHydrogens();

      //! Deletes all atoms except for the largest contiguous fragment
      bool StripSalts();
      //! Converts the charged form of coordinate bonds, e.g.[N+]([O-])=O to N(=O)=O 
      bool ConvertDativeBonds();

      bool CorrectForPH();
      bool AssignSpinMultiplicity();
      vector3 Center(int nconf);
      //! Set the torsion defined by these atoms, rotating bonded neighbors
      void SetTorsion(OBAtom*,OBAtom*,OBAtom*,OBAtom*,double);
      //@}

      //! \name Molecule utilities and perception methods
      //@{
      //! Find Smallest Set of Smallest Rings (see OBRing class for more details)
      void FindSSSR();
      void FindRingAtomsAndBonds();
      void FindChiralCenters();
      void FindChildren(std::vector<int> &,int,int);
      void FindChildren(std::vector<OBAtom*>&,OBAtom*,OBAtom*);
      void FindLargestFragment(OBBitVec &);
      //! Sort a list of contig fragments by size from largest to smallest
      //! Each vector<int> contains the atom numbers of a contig fragment
      void ContigFragList(std::vector<std::vector<int> >&);
      //! Aligns atom a on p1 and atom b along p1->p2 vector
      void Align(OBAtom*,OBAtom*,vector3&,vector3&);
      //! Adds single bonds based on atom proximity
      void ConnectTheDots();
      //! Attempts to perceive multiple bonds based on geometries
      void PerceiveBondOrders();
      void FindTorsions();
      // documented in mol.cpp: graph-theoretical distance for each atom
      bool         GetGTDVector(std::vector<int> &);
      // documented in mol.cpp: graph-invariant index for each atom
      void         GetGIVector(std::vector<unsigned int> &);
      // documented in mol.cpp: calculate symmetry-unique identifiers
      void         GetGIDVector(std::vector<unsigned int> &);
      //@}

      //! \name Methods to check for existence of properties
      //@{
      //! Are there non-zero coordinates in two dimensions (i.e. X and Y)?
      bool Has2D();
      //! Are there non-zero coordinates in all three dimensions (i.e. X, Y, Z)?
      bool Has3D();
      //! Are there any non-zero coordinates?
      bool HasNonZeroCoords();
      bool HasAromaticPerceived()     { return(HasFlag(OB_AROMATIC_MOL)); }
      bool HasSSSRPerceived()         { return(HasFlag(OB_SSSR_MOL));     }
      bool HasRingAtomsAndBondsPerceived(){return(HasFlag(OB_RINGFLAGS_MOL));}
      bool HasAtomTypesPerceived()    { return(HasFlag(OB_ATOMTYPES_MOL));}
      bool HasChiralityPerceived()    { return(HasFlag(OB_CHIRALITY_MOL));}
      bool HasPartialChargesPerceived() { return(HasFlag(OB_PCHARGE_MOL));}
      bool HasHybridizationPerceived() { return(HasFlag(OB_HYBRID_MOL));  }
      bool HasImplicitValencePerceived() { return(HasFlag(OB_IMPVAL_MOL));}
      bool HasKekulePerceived() { return(HasFlag(OB_KEKULE_MOL));         }
      bool HasClosureBondsPerceived() { return(HasFlag(OB_CLOSURE_MOL));  }
      bool HasChainsPerceived() { return(HasFlag(OB_CHAINS_MOL));         }
      bool HasHydrogensAdded() { return(HasFlag(OB_H_ADDED_MOL));         }
      bool HasAromaticCorrected() { return(HasFlag(OB_AROM_CORRECTED_MOL));}
      bool IsCorrectedForPH() { return(HasFlag(OB_PH_CORRECTED_MOL));     }
      bool HasSpinMultiplicityAssigned() { return(HasFlag(OB_TSPIN_MOL)); }
      //! Is this molecule chiral?
      bool IsChiral();
      //! Are there any atoms in this molecule?
      bool Empty()                       { return(_natoms == 0);          }
      //@}

      //! \name Multiple conformer member functions
      //@{
      int     NumConformers()    { return((_vconf.empty())?0:_vconf.size()); }
      void    SetConformers(std::vector<double*> &v);
      void    AddConformer(double *f)    {  _vconf.push_back(f);    }
      void    SetConformer(int i)        {  _c = _vconf[i];         }
      void    CopyConformer(double*,int);
      void    DeleteConformer(int);
      double  *GetConformer(int i)       {  return(_vconf[i]);      }
      double  *BeginConformer(std::vector<double*>::iterator&i)
        { i = _vconf.begin();
        return((i == _vconf.end()) ? NULL:*i); }
      double  *NextConformer(std::vector<double*>::iterator&i)
        { i++;
        return((i == _vconf.end()) ? NULL:*i); }
      std::vector<double*> &GetConformers() {   return(_vconf);     }
      //@}

      //! \name Iterator methods
      //@{
      //! \deprecated Use FOR_ATOMS_OF_MOL and OBMolAtomIter instead
      OBAtom *BeginAtom(std::vector<OBNodeBase*>::iterator &i);
      //! \deprecated Use FOR_ATOMS_OF_MOL and OBMolAtomIter instead
      OBAtom *NextAtom(std::vector<OBNodeBase*>::iterator &i);
      //! \deprecated Use FOR_BONDS_OF_MOL and OBMolBondIter instead
      OBBond *BeginBond(std::vector<OBEdgeBase*>::iterator &i);
      //! \deprecated Use FOR_BONDS_OF_MOL and OBMolBondIter instead
      OBBond *NextBond(std::vector<OBEdgeBase*>::iterator &i);
      //! \deprecated Use FOR_RESIDUES_OF_MOL and OBResidueIter instead
      OBResidue *BeginResidue(std::vector<OBResidue*>::iterator &i)
        {
          i = _residue.begin();
          return((i == _residue.end()) ? NULL:*i);
        }
      //! \deprecated Use FOR_RESIDUES_OF_MOL and OBResidueIter instead
      OBResidue *NextResidue(std::vector<OBResidue*>::iterator &i)
        {
          i++;
          return((i == _residue.end()) ? NULL:*i);
        }
      OBInternalCoord *BeginInternalCoord(std::vector<OBInternalCoord*>::iterator &i)
        {
          i = _internals.begin();
          return((i == _internals.end()) ? NULL:*i);
        }
      OBInternalCoord *NextInternalCoord(std::vector<OBInternalCoord*>::iterator &i)
        {
          i++;
          return((i == _internals.end()) ? NULL:*i);
        }
      //@}

      //  Removed with OBConversion framework -- see OBConversion class instead
      //! \name Convenience functions for I/O
      //@{
      // friend std::ostream&       operator<< ( std::ostream&, OBMol& ) ;
      // friend std::istream&       operator>> ( std::istream&, OBMol& ) ;
      //@}
    };

  //! \brief Used to transform from z-matrix to cartesian coordinates.
  class OBAPI OBInternalCoord
    {
    public:
      //class members
      OBAtom *_a,*_b,*_c;
      double   _dst,_ang,_tor;
      //! Constructor
      OBInternalCoord(OBAtom *a=(OBAtom*)NULL,
                      OBAtom *b=(OBAtom*)NULL,
                      OBAtom *c=(OBAtom*)NULL)
        {
          _a = a;
          _b = b;
          _c = c;
          _dst = _ang = _tor = 0.0;
        }
    };

  //function prototypes

  OBAPI bool tokenize(std::vector<std::string>&, const char *buf, const char *delimstr=" \t\n");
  OBAPI bool tokenize(std::vector<std::string>&, std::string&, const char *delimstr=" \t\n", int limit=-1);
  //! remove leading and trailing whitespace from a string
  OBAPI std::string& Trim(std::string& txt);
  //! \deprecated -- use OBMessageHandler class instead
  OBAPI void ThrowError(char *str);
  //! \deprecated -- use OBMessageHandler class instead
  OBAPI void ThrowError(std::string &str);
  OBAPI void CartesianToInternal(std::vector<OBInternalCoord*>&,OBMol&);
  OBAPI void InternalToCartesian(std::vector<OBInternalCoord*>&,OBMol&);
  OBAPI std::string NewExtension(std::string&,char*);
  // Now handled by OBConversion class
  // OBAPI bool SetInputType(OBMol&,std::string&);
  // OBAPI bool SetOutputType(OBMol&,std::string&);

  //global definitions
  //! Global OBElementTable for element properties
  EXTERN  OBElementTable   etab;
  //! Global OBTypeTable for translating between different atom types
  //! (e.g., Sybyl <-> MM2)
  EXTERN  OBTypeTable      ttab;
  //! Global OBIsotopeTable for isotope properties
  EXTERN  OBIsotopeTable   isotab;
  //! Global OBAromaticTyper for detecting aromatic atoms and bonds
  EXTERN  OBAromaticTyper  aromtyper;
  //! Global OBAtomTyper for marking internal valence, hybridization,
  //!  and atom types (for internal and external use)
  EXTERN  OBAtomTyper      atomtyper;
  //! Global OBChainsParser for detecting macromolecular chains and residues
  EXTERN  OBChainsParser   chainsparser;
  //! Global OBMessageHandler error handler
  OBERROR extern  OBMessageHandler obErrorLog;
  //! Global OBResidueData biomolecule residue database
  EXTERN  OBResidueData    resdat;

  //Utility Macros

#ifndef BUFF_SIZE
#define BUFF_SIZE 32768
#endif

#ifndef EQ
#define EQ(a,b) (!strcmp((a), (b)))
#endif

#ifndef EQn
#define EQn(a,b,n) (!strncmp((a), (b), (n)))
#endif

#ifndef SQUARE
#define SQUARE(x) ((x)*(x))
#endif

#ifndef IsUnsatType
#define IsUnsatType(x)  (EQ(x,"Car") || EQ(x,"C2") || EQ(x,"Sox") || EQ(x,"Sac") || EQ(x,"Pac") || EQ(x,"So2"))
#endif

#ifndef __KCC
  extern "C"
  {
    OBAPI void  get_rmat(double*,double*,double*,int);
    OBAPI void  ob_make_rmat(double mat[3][3],double rmat[9]);
    OBAPI void  qtrfit (double *r,double *f,int size,double u[3][3]);
    OBAPI double superimpose(double*,double*,int);
  }
#else
  OBAPI void get_rmat(double*,double*,double*,int);
  OBAPI void ob_make_rmat(double mat[3][3],double rmat[9]);
  OBAPI void qtrfit (double *r,double *f,int size,double u[3][3]);
  OBAPI double superimpose(double*,double*,int);
#endif // __KCC

} // end namespace OpenBabel

#endif // OB_MOL_H

//! \file
//! \brief Handle molecules. Declarations of OBMol, OBAtom, OBBond, OBResidue.
//!        (the main header for Open Babel)

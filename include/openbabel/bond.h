/**********************************************************************
bond.h - Handle OBBond class.

Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2006 by Geoffrey R. Hutchison
Some portions Copyright (C) 2003 by Michael Banck
Some portions Copyright (C) 2008 by Tim Vandermeersch

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

#ifndef OB_BOND_H
#define OB_BOND_H

#include <openbabel/babelconfig.h>

#ifndef EXTERN
#  define EXTERN extern
#endif

#include <openbabel/base.h>
#include <openbabel/atom.h>

namespace OpenBabel
{
  class OBAtom;
  class OBRing;

  //BOND Property Macros (flags)
  //! An aromatic bond (regardless of bond order)
#define OB_AROMATIC_BOND  (1<<1)
  //! A solid black wedge in 2D representations -- i.e., "up" from the 2D plane
#define OB_WEDGE_BOND     (1<<2)
  //! A dashed "hash" bond in 2D representations -- i.e., "down" from the 2D plane
#define OB_HASH_BOND      (1<<3)
  //! A bond in a ring
#define OB_RING_BOND      (1<<4)
  //! A bond which "closes" a ring when walking the molecular graph
#define OB_CLOSURE_BOND   (1<<10)
  // 11-16 currently unused
#define OB_WEDGE_OR_HASH_BOND     (1<<11)

#define SET_OR_UNSET_FLAG(X) \
  if (value) SetFlag(X); \
  else     UnsetFlag(X);


  class OBAPI OBBond: public OBBase
  {
    protected:
      unsigned int                _idx;   //!< Unique edge index used by GetIdx() and SetIdx()
      OBMol                      *_parent;//!< The molecule which contains me (if any)
      OBAtom                     *_bgn;   //!< I connect one node
      OBAtom                     *_end;   //!< to another node
      char                        _order; //!< Bond order (1, 2, 3, 5=aromatic)
      unsigned short int          _flags; //!< Any flags for this bond
      unsigned long                 _id;        //!< unique id
      //OBBondPrivate * const d;

      /**
      * @return True id the @p flag is set.
       */
      bool HasFlag(int flag) const { return ((_flags & flag) != 0); }
      /**
       * Sets the bitwise @p flag
       */
      void SetFlag(int flag) { _flags |= flag; }
      /**
       * Unsets the bitwise @p flag
       */
      void UnsetFlag(int flag) { _flags &= (~(flag)); }

    public:
      enum Flag {
        Aromatic = (1<<1), //!< An aromatic bond (regardless of bond order)
        Ring     = (1<<4), //!< A bond in a ring
        Closure  = (1<<10) //!< A bond which "closes" a ring when walking the molecular graph
      };
      enum StereoFlag {
        Wedge       = (1<<2),  //!< A solid black wedge in 2D representations -- i.e., "up" from the 2D plane
        Hash        = (1<<3),  //!< A dashed "hash" bond in 2D representations -- i.e., "down" from the 2D plane
        WedgeOrHash = (1<<11), //!< The bond is either wedge or hash, this is a seperate flag!
        CisOrTrans  = (1<<12)  //!< Indicates the 2D/3D coordinates are accidently cis/trans.
      };
      //! Whether this bond has been visited by a graph algorithm
      /** \deprecated Use OBBitVec objects instead to be fully thread-safe. **/
      bool Visit;

      //! Constructor
      OBBond();
      //! Destructor
      virtual ~OBBond();

      //! \name Bond modification methods
      //@{
      //! Set the internal bond index
      /** \warning This will not update the index in the parent OBMol.
          Intended mainly for internal use. Use with care. **/
      void SetIdx(int idx)        {          _idx = idx;        }
      void SetId(unsigned long id) { _id = id; }
      //! Set the bond order to @p order (i.e., 1 = single, 2 = double, 5 = aromatic)
      void SetBondOrder(int order);
      //! Set the beginning atom of this bond to @p begin. Does not update @p begin.
      void SetBegin(OBAtom *begin){          _bgn = begin;      }
      //! Set the ending atom of this bond to @p end. Does not update @p end.
      void SetEnd(OBAtom *end)    {          _end = end;        }
      //! Set the parent molecule to @p ptr. Does not update parent.
      void SetParent(OBMol *ptr)  {        _parent= ptr;        }
      //! Change the bond length to @p length, while keeping @p fixed stationary
      void SetLength(OBAtom *fixed,double length);
      //! Change the bond length to @p length, moving both atoms halfway
      //! \since version 2.2
      void SetLength(double length);
      //! Set the main bond information (i.e., when creating a bond)
      void Set(int index, OBAtom* begin,OBAtom* end,int order,int flags);
      //! Mark that this bond is aromatic. Does not update atoms or validate.
      void SetAromatic(bool value=true)    { SET_OR_UNSET_FLAG(OB_AROMATIC_BOND); }
      /**
       * Mark that this bond has 2D "wedge" notation (i.e., goes in a positive
       * Z direction from the beginning to end atoms)
       */
      void SetWedge(bool value=true) { SET_OR_UNSET_FLAG(Wedge); }
      /**
       * Mark that this bond has 2D "hash" notation (i.e., goes in a negative
       * Z direction from the beginning to end atoms)
       */
      void SetHash(bool value=true) { SET_OR_UNSET_FLAG(Hash); }
      /**
       * Set the WedgeOrHash flag on a bond (??)
       */
      void SetWedgeOrHash(bool value=true) { SET_OR_UNSET_FLAG(WedgeOrHash); }
      //! Mark that this bond is in a ring. Primarily for internal use.
      void SetInRing(bool value=true) { SET_OR_UNSET_FLAG(OB_RING_BOND); }
      //! Mark that this bond indicates a ring closure when walking the molecule
      /** \warning This is for internal use only. All closure bonds are marked
          automatically by lazy evaluation when requesting
          OBBond::IsClosure() **/
      void SetClosure(bool value=true)     { SET_OR_UNSET_FLAG(OB_CLOSURE_BOND); }
      //@}

      //! \name Bond data request methods
      //@{
      //! \return The unique bond index in a molecule.
      unsigned int     GetIdx()           const { return(_idx);  }
      unsigned long GetId()           const { return _id; }
      //! \return The bond order for the bond
      unsigned int     GetBondOrder()     const { return(_order); }
      //! \return The set of property flags defined for this bond.
      unsigned int     GetFlags()         const { return(_flags);      }
      //! \return The atom index for the end atom in this bond (from OBAtom::GetIdx()
      unsigned int     GetBeginAtomIdx()  const
        { return (_bgn ? _bgn->GetIdx() : 0); }
      //! \return The atom index for the end atom in this bond (from OBAtom::GetIdx()
      unsigned int     GetEndAtomIdx()    const
        { return (_end ? _end->GetIdx() : 0); }
      //! \return The "beginning" atom for this bond
      OBAtom *GetBeginAtom()    { return(_bgn);    }
      const OBAtom *GetBeginAtom() const
        { return(_bgn);    }
      //! \return The "end" atom for this bond
      OBAtom *GetEndAtom()      { return(_end);    }
      const OBAtom *GetEndAtom() const
        { return(_end);    }
      //! \return The neighboring atom to @p ptr (i.e., the end if @p ptr is the start)
      /** \warning If @p ptr is not part of the bond, the beginning atom
          will always be returned **/
      OBAtom *GetNbrAtom(OBAtom *ptr)
        {
          return((ptr != _bgn)? _bgn : _end);
        }
      //! \return The enclosing OBMol for this bond, or NULL if none is defined.
      OBMol  *GetParent()                 {return(_parent);}
      //! \return The expected "equilibrium" length based on the covalent radii and bond order
      /** Length is given in Angstroms **/
      double  GetEquibLength() const;
      //! \return The current length of this bond in Angstroms
      double  GetLength() const;
      //! \return The index to the neighboring atom of @p ptr (i.e., the end if @p ptr is the start)
      /** \warning If @p ptr is not part of the bond, the beginning atom
          index will always be returned **/
      unsigned int     GetNbrAtomIdx(OBAtom *ptr)
        {
          if (ptr!=_bgn)
            return (_bgn ? _bgn->GetIdx() : 0);
          else
            return (_end ? _end->GetIdx() : 0);
        }
      //! Find the smallest ring containing this bond (returns a NULL pointer if none exists)
      OBRing* FindSmallestRing() const;
      //@}

      //! \name property request methods
      //@{
      //! \return Is the bond aromatic?
      //!  (Note that the two atoms of the bond may be aromatic,
      //!   but not the bond)
      bool IsAromatic() const;
      //! \return Is the bond part of a ring?
      bool IsInRing() const;
      //! \return Is the bond a rotatable bond?
      /**  Currently, this function classifies any bond with at least one heavy
           atom, no sp-hybrid atoms (e.g., a triple bond somewhere) not in a ring
           as a potential rotor if includeRingsBonds is false.  If true, rotors in
           rings with more than 3 atoms may be included. No other bond typing is attempted.
           For more detailed rotor detection, check the OBRotorList and
           OBRotorRules classes **/
      bool IsRotor(bool includeRingBonds=false);
      /** \return Is the bond an amide link (i.e., between a carbonyl C and a N)?
           No distinction is made between primary, secondary, and tertiary amides. **/
      bool IsAmide();
      /** \return Is the bond a primary amide (i.e., between carbonyl C and a NH2)?
           In versions prior to 2.3, this function incorrectly identified secondary amides. **/
      bool IsPrimaryAmide();
      /** \return Is the bond a secondary amide (i.e., between a carbonyl C and a NH1)?
           In versions prior to 2.3, this function incorrectly identified tertiary amides. **/
      bool IsSecondaryAmide();
      //! \return Is the bond a teriary amide (i.e., between a carbonyl C and a NH0)?
      //!  \since version 2.3.
      bool IsTertiaryAmide();
      //! \return Is the bond an ester link (i.e., between a carbonyl C and an O)?
      bool IsEster();
      //! \return Is the bond a carbonyl C=O?
      bool IsCarbonyl();
      //! \return Does this bond "close" a ring when walking the molecular graph?
      bool IsClosure();
      /** \return Whether this bond is a "wedge" in 2D representations
          (i.e., goes in a positive Z direction from the beginning to end atoms) **/
      bool IsWedge() {    return(HasFlag(OB_WEDGE_BOND));    }
      /** \return Whether this bond is a "hash" in 2D representations
          (i.e., goes in a negative Z direction from the beginning to end atoms) **/
      bool IsHash()  {    return(HasFlag(OB_HASH_BOND));     }
      /**
       * @return True if this bond is either a wedge or hash.
       * @note: This is a seperate bond type
       * @since version 2.3
       */
      bool IsWedgeOrHash() const { return(HasFlag(WedgeOrHash)); }
      /**
       * @return True if this bond is either a cis or trans.
       * @since version 2.3
       */
      bool IsCisOrTrans() const { return(HasFlag(CisOrTrans)); }

      //! \return whether the geometry around this bond "looks" unsaturated
      bool IsDoubleBondGeometry();
      //@}

    }; // class OBBond

  //! A standard iterator over a vector of bonds
  typedef std::vector<OBBond*>::iterator OBBondIterator;

}// namespace OpenBabel

#endif   // OB_BOND_H

//! @file bond.h
//! @brief Handle bonds

/**********************************************************************
bond.h - Handle OBBond class.
 
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2006 by Geoffrey R. Hutchison
Some portions Copyright (C) 2003 by Michael Banck
Some portions Copyright (C) 2008 by Tim Vandermeersch
 
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

#ifndef OB_BOND_H
#define OB_BOND_H

#include <openbabel/babelconfig.h>

#ifndef EXTERN
#  define EXTERN extern
#endif

#include <openbabel/base.h>

namespace OpenBabel
{
  class OBAtom;
  class OBMol;

  //! OBBond flags (Aromatic, InRing, ...)
  namespace OBBondFlag {
    enum {
      //! An aromatic bond (regardless of bond order)
      Aromatic = (1<<1),
      //! A solid black wedge in 2D representations -- i.e., "up" from the 2D plane
      Wedge    = (1<<2),
      //! A dashed "hash" bond in 2D representations -- i.e., "down" from the 2D plane
      Hash     = (1<<4),
      //! A bond in a ring
      Ring     = (1<<5),
      //! The "upper" bond in a double bond cis/trans isomer (i.e., "/" in SMILES)
      Up       = (1<<6),
      //! The "down" bond in a double bond cis/trans isomer (i.e., "\" in SMILES)
      Down     = (1<<7),
      //! A Kekule single bond
      Single   = (1<<8),
      //! A Kekule double bond
      Double   = (1<<9),
      //! A Kekule triple bond
      Triple   = (1<<10),
      //! A bond which "closes" a ring when walking the molecular graph
      Closure  = (1<<11)
      // 11-16 currently unused
    };
  };

  // class introduction in bond.cpp
  class OBBondPrivate;
  class OBAPI OBBond: public OBBase
  {
    protected:
      // Some protected data declared in the class itself so we can access it 
      // through inline functions
      unsigned int                m_idx;   //!< Unique edge index used by GetIdx() and SetIdx()
      OBMol                      *m_parent;//!< The molecule which contains me (if any)
      OBAtom                     *m_bgn;   //!< I connect one node
      OBAtom                     *m_end;   //!< to another node
      char                        m_order; //!< Bond order (1, 2, 3, 5=aromatic)
      unsigned short int          m_flags; //!< Any flags for this bond
      //! The private d pointer for which the content can be changed in 
      //! futute versions without breaking binary compatibility.
      OBBondPrivate * const d;
     
      /** 
       * @return True id the @p flag is set.
       */
      bool HasFlag(int flag);
      /** 
       * Sets the bitwise @p flag
       */
      void SetFlag(int flag);
      /** 
       * Unsets the bitwise @p flag
       */
      void UnsetFlag(int flag);
    

    public:
      //! Whether this bond has been visited by a graph algorithm
      /** \deprecated Use OBBitVec objects instead to be fully thread-safe. **/
      //bool Visit;

      //! Constructor
      OBBond();
      //! Destructor
      virtual ~OBBond();

      //! \name Bond modification methods
      //@{
      //! Set the internal bond index
      /** \warning This will not update the index in the parent OBMol.
          Intended mainly for internal use. Use with care. **/
      void SetIdx(int idx);
      //! Set the bond order to @p order (i.e., 1 = single, 2 = double, 5 = aromatic)
      /** \deprecated Use SetBondOrder() instead. **/
      void SetBO(int order);
      //! Set the bond order to @p order (i.e., 1 = single, 2 = double, 5 = aromatic)
      void SetBondOrder(int order);
      //! Set the beginning atom of this bond to @p begin. Does not update @p begin.
      void SetBegin(OBAtom *begin);
      //! Set the ending atom of this bond to @p end. Does not update @p end.
      void SetEnd(OBAtom *end);    
      //! Set the parent molecule to @p ptr. Does not update parent.
      void SetParent(OBMol *ptr);
      //! Change the bond length to @p length, while keeping @p fixed stationary
      void SetLength(OBAtom *fixed,double length);
      //! Change the bond length to @p length, moving both atoms halfway
      //! \since version 2.2
      void SetLength(double length);
      //! Set the main bond information (i.e., when creating a bond)
      void Set(int index, OBAtom* begin, OBAtom* end, int order, int flags);
      //! \deprecated Use SetBondOrder() instead
      void SetKSingle();
      //! \deprecated Use SetBondOrder() instead
      void SetKDouble();
      //! \deprecated Use SetBondOrder() instead
      void SetKTriple();
      //! Mark that this bond is aromatic. Does not update atoms or validate.
      void SetAromatic()    { SetFlag(OBBondFlag::Aromatic); }
      //! Mark that this bond has 2D "hash" notation (i.e., goes in a negative Z direction from the beginning to end atoms)
      void SetHash()        { SetFlag(OBBondFlag::Hash);     }
      //! Mark that this bond has 2D "wedge" notation (i.e., goes in a positive Z direction from the beginning to end atoms)
      void SetWedge()       { SetFlag(OBBondFlag::Wedge);    }
      //! Mark that this bond has an "up" torsion for double-bond stereochem (i.e., "/" in SMILES notation
      void SetUp()          { SetFlag(OBBondFlag::Up); UnsetFlag(OBBondFlag::Down); }
      //! Mark that this bond has an "down" torsion for double-bond stereochem (i.e., "\" in SMILES notation
      void SetDown()        { SetFlag(OBBondFlag::Down); UnsetFlag(OBBondFlag::Up);   }
      //! Mark that this bond is in a ring. Primarily for internal use.
      void SetInRing()      { SetFlag(OBBondFlag::Ring);     }
      //! Mark that this bond indicates a ring closure when walking the molecule
      /** \warning This is for internal use only. All closure bonds are marked
          automatically by lazy evaluation when requesting 
          OBBond::IsClosure() **/
      void SetClosure()     { SetFlag(OBBondFlag::Closure);  }
      //! Clear any indication of 2D "hash" notation from SetHash()
      void UnsetHash()      { UnsetFlag(OBBondFlag::Hash);    }
      //! Clear any indication of 2D "wedge" notation from SetWedge()
      void UnsetWedge()     { UnsetFlag(OBBondFlag::Wedge);   }
      //! Clear any indication of "/" double bond stereochemistry from SetUp()
      void UnsetUp()        { UnsetFlag(OBBondFlag::Up);   }
      //! Clear any indication of "\" double bond stereochemistry from SetDown()
      void UnsetDown()      { UnsetFlag(OBBondFlag::Down); }
      //! Clear all aromaticity information for the bond
      void UnsetAromatic()  { UnsetFlag(OBBondFlag::Aromatic);}
      //! Clear all Kekule information for the bond
      void UnsetKekule();
      //@}

      //! \name Bond data request methods
      //@{
      //! \return The unique bond index in a molecule.
      unsigned int GetIdx() const;
      //! \return The bond order for the bond
      /** \deprecated Use GetBondOrder() as this method may be removed. **/
      unsigned int GetBO() const;
      //! \return The bond order for the bond
      unsigned int GetBondOrder() const;
      //! \return The set of property flags defined for this bond.
      unsigned int GetFlags() const;
      //! \return The atom index for the end atom in this bond (from OBAtom::GetIdx()
      unsigned int GetBeginAtomIdx() const; 
      //! \return The atom index for the end atom in this bond (from OBAtom::GetIdx()
      unsigned int GetEndAtomIdx() const;
      //! \return The "beginning" atom for this bond
      OBAtom *GetBeginAtom();
      const OBAtom *GetBeginAtom() const; 
      //! \return The "end" atom for this bond
      OBAtom *GetEndAtom();
      const OBAtom *GetEndAtom() const;
      //! \return The neighboring atom to @p ptr (i.e., the end if @p ptr is the start)
      /** \warning If @p ptr is not part of the bond, the beginning atom
          will always be returned **/
      OBAtom *GetNbrAtom(OBAtom *ptr);
      //! \return The enclosing OBMol for this bond, or NULL if none is defined.
      OBMol  *GetParent();
      //! \return The expected "equilibrium" length based on the covalent radii and bond order
      /** Length is given in Angstroms **/
      double  GetEquibLength() const;
      //! \return The current length of this bond in Angstroms
      double  GetLength() const;
      //! \return The index to the neighboring atom of @p ptr (i.e., the end if @p ptr is the start)
      /** \warning If @p ptr is not part of the bond, the beginning atom
          index will always be returned **/
      unsigned int GetNbrAtomIdx(OBAtom *ptr);
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
           as a potential rotor. No other bond typing is attempted.
           For more detailed rotor detection, check the OBRotorList and 
           OBRotorRules classes **/
      bool IsRotor();
      //! \return Is the bond an amide link (i.e., between a carbonyl C and a N)
      bool IsAmide();
      //! \return Is the bond an amide (i.e., between carbonyl C and a NH group)
      bool IsPrimaryAmide();
      //! \return Is the bond an amide between a carbonyl C and a N with no hydrogens
      bool IsSecondaryAmide();
      //! \return Is the bond an ester link (i.e., between a carbonyl C and an O)
      bool IsEster();
      //! \return Is the bond a carbonyl C=O?
      bool IsCarbonyl();
      //! \return Is the bond a single bond?
      bool IsSingle();
      //! \return Is the bond is a double bond?
      bool IsDouble();
      //! \return Is the bond is a triple bond?
      bool IsTriple();
      //! \deprecated Use IsSingle() instead
      bool IsKSingle();
      //! \deprecated Use IsDouble() instead
      bool IsKDouble();
      //! \deprecated Use IsTriple() instead
      bool IsKTriple();
      //! \return Does this bond "close" a ring when walking the molecular graph?
      bool IsClosure();
      /** \return Whether this is the "upper" bond in a double bond cis/trans
          isomer (i.e., "/" in SMILES) **/
      bool IsUp()    {    return(HasFlag(OBBondFlag::Up));    }
      /** \return Whether this is the "lower" bond in a double bond cis/trans
          isomer (i.e., "\" in SMILES) **/
      bool IsDown()  {    return(HasFlag(OBBondFlag::Down));  }
      /** \return Whether this bond is a "wedge" in 2D representations
          (i.e., goes in a positive Z direction from the beginning to end atoms) **/
      bool IsWedge() {    return(HasFlag(OBBondFlag::Wedge));    }
      /** \return Whether this bond is a "hash" in 2D representations
          (i.e., goes in a negative Z direction from the beginning to end atoms) **/
      bool IsHash()  {    return(HasFlag(OBBondFlag::Hash));     }
      //! \return whether the geometry around this bond "looks" unsaturated
      bool IsDoubleBondGeometry();
      //@}

    }; // class OBBond

  //! A standard iterator over a vector of bonds
  typedef std::vector<OBBond*>::iterator OBBondIterator;

}// namespace OpenBabel

#endif   // OB_BOND_H

//! \file bond.h
//! \brief Handle bonds

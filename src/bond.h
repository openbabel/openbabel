/**********************************************************************
bond.h - Handle OBBond class.
 
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

#ifndef OB_BOND_H
#define OB_BOND_H

#include "babelconfig.h"

#ifndef EXTERN
#  define EXTERN extern
#endif

#include "base.h"
#include "atom.h"

namespace OpenBabel
{

  class OBAtom;

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
      //! \return Return true of the bond is a single bond
      bool IsSingle();
      //! \return Return true of the bond is a double bond
      bool IsDouble();
      //! \return Return true of the bond is a tripple bond
      bool IsTriple();
      //! \deprecated Use IsSingle() instead
      bool IsKSingle();
      //! \deprecated Use IsDouble() instead
      bool IsKDouble();
      //! \deprecated Use IsTriple() instead
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

}// namespace OpenBabel

#endif   // OB_BOND_H

//! \file bond.h
//! \brief Handle bonds

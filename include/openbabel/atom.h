/**********************************************************************
atom.h - Handle OBAtom class.
 
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

#ifndef OB_ATOM_H
#define OB_ATOM_H

#include <openbabel/babelconfig.h>

#ifndef EXTERN
#  define EXTERN extern
#endif

#include <openbabel/base.h>
#include <openbabel/math/vector.h>

#include <vector>
#include <string>

namespace OpenBabel
{
  // forward declarations...
  class OBAtom;
  class OBBond;
  class OBResidue;
  class OBMol;

  //! A standard iterator over a vector of bonds
  typedef std::vector<OBBond*>::iterator OBBondIterator;
  //! A standard iterator over a vector of atoms
  typedef std::vector<OBAtom*>::iterator OBAtomIterator;

  namespace OBAtomFlag {
    enum {
      //! Atom is in a 4-membered ring
      Ring4         = (1<<1),
      //! Atom is in a 3-membered ring
      Ring3         = (1<<2),
      //! Atom is aromatic
      AromaticAtom  = (1<<3),
      //! Atom is in a ring
      Ring          = (1<<4),
      //! Atom has clockwise SMILES chiral stereochemistry (i.e., "@@")
      StereoC       = (1<<5),
      //! Atom has anticlockwise SMILES chiral stereochemistry (i.e., "@")
      StereoAC      = (1<<6),
      //! Atom is an electron donor
      Donor         = (1<<7),
      //! Atom is an electron acceptor
      Acceptor      = (1<<8),
      //! Atom is chiral
      Chiral        = (1<<9),
      //! Atom has + chiral volume
      ChiralPos     = (1<<10),
      //! Atom has - chiral volume
      ChiralNeg     = (1<<11),
      //! Atom has no hydrogen attached. Temporary use only during input of some formats
      HasNoH        = (1<<12),
      //! Atom is not hydrogen deficient. (for SMILES input)
      NotHDeficient = (1<<13)
    };
  };

  /// @addtogroup core Core classes 
  //@{

  class OBAtomPrivate;
  class OBAPI OBAtom: public OBBase
  {
    private:   
      // Some protected data declared in the class itself so we can access it 
      // through inline functions
      unsigned int              m_idx;  //!< unique node index (GetIdx(), SetIdx())
      unsigned int              m_cidx; //!< index into coordinate array
      double                  **m_cptr; //!< coordinate array in double*
      mutable Eigen::Vector3d   m_pos;  //!< coordinate vector
      //! The private d pointer for which the content can be changed in 
      //! futute versions without breaking binary compatibility.
      OBAtomPrivate * const d;

      /** 
       * @return All flags
       */
      int  GetFlag() const;
      /** 
       * Sets The bitwise @p flag
       */
      void SetFlag(int flag);
      /** 
       * @return True of the atom has the @p flag
       */
      bool HasFlag(int flag);

    public:
      /** 
       * Constructor.
       */
      OBAtom();
      /**
       * Destructor.
       */
      virtual ~OBAtom();
      /** 
       * Assignment operator.
       */
      OBAtom& operator=(OBAtom &);
      /** 
       * Duplicate another atom. Copies all information with the exception of index
       * @since version 2.2
       */
      void Duplicate(OBAtom *);
      /** 
       * Clear all data. Calls OBBase::Clear() to handle any generic data.
       * @return True if successful.
       */
      bool Clear();

      //! @name Methods to set atomic information
      //@{
      /** 
       * Set atom index (i.e., in an OBMol).
       */
      void SetIdx(unsigned int idx);
      /** 
       * Set atom hybridization. 
       * @param hyb The hybridization (i.e., 1 = sp, 2 = sp2, 3 = sp3 ...).
       */
      void SetHyb(int hyb);
      /** 
       * Set atomic number.
       * @param atomicnum The atomic number (i.e., 1 = Hydrogen, 2 = Helium, ...).
       */
      void SetAtomicNum(int atomicnum);
      /** 
       * Set isotope number (actual atomic weight is tabulated automatically, 0 = most abundant).
       */
      void SetIsotope(unsigned int iso);
      /** 
       * Set implicit valence to @p val.
       */
      void SetImplicitValence(int val);
      /** 
       * Increase the implicit valence by one.
       */
      void IncrementImplicitValence();
      /** 
       * Decrease the implicit valence by one.
       */
      void DecrementImplicitValence();
      /** 
       * Set the formal charge of the atom to @p fcharge.
       */
      void SetFormalCharge(int fcharge);
      /** 
       * Set the atomic spin to @p spin. See d->spinmultiplicity.
       */
      void SetSpinMultiplicity(short spin);
      /** 
       * Set the atomic type symbol (see OBTypeTable and OBAtomTyper for more).
       * @warning OpenBabel does not attempt to interpret this type, use 
       * SetAtomicNum() to change an atom's element.
       */
      void SetType(const char *type);
      /** 
       * Set the atomic type symbol (see OBTypeTable and OBAtomTyper for more).
       * @warning OpenBabel does not attempt to interpret this type, use 
       * SetAtomicNum() to change an atom's element.
       */
      void SetType(const std::string &type);
      /** 
       * Set the partial charge to @p pcharge.
       */
      void SetPartialCharge(double pcharge);
      /** 
       * Set the coordinate vector for this atom to @p v as a Eigen::Vector3d.
       */
      void SetVector(const Eigen::Vector3d &v);
      /** 
       * Set the coordinate vector for this atom based on @p x @p y & @p z.
       */
      void SetVector(const double x, const double y, const double z);
      /** 
       * Set the position of this atom from a pointer-driven array of coordinates.
       */
      void SetCoordPtr(double **c);
      /** 
       * Set the position of this atom based on the internal pointer array (i.e. from SetCoordPtr() ).
       */
      void SetVector();
      /** 
       * Attach an OBResidue @p res as containing this atom.
       */
      void SetResidue(OBResidue *res);
      /** 
       * Attach an OBMol @p ptr as the parent container for this atom.
       */
      void SetParent(OBMol *ptr);
      /** 
       * Mark atom as being aromatic.
       */
      void SetAromatic() { SetFlag(OBAtomFlag::AromaticAtom); }
      /** 
       * Clear aromatic information from the atom.
       */
      void UnsetAromatic();
      /** 
       * Mark atom as having SMILES clockwise stereochemistry (i.e., "@@").
       */
      void SetClockwiseStereo() { SetFlag(OBAtomFlag::StereoC | OBAtomFlag::Chiral); }
      /** 
       * Mark atom as having SMILES anticlockwise stereochemistry (i.e., "@").
       */
      void SetAntiClockwiseStereo() { SetFlag(OBAtomFlag::StereoAC | OBAtomFlag::Chiral); }
      /** 
       * Mark an atom as having + chiral volume.
       */
      void SetPositiveStereo() { SetFlag(OBAtomFlag::ChiralPos | OBAtomFlag::Chiral); }
      /** 
       * Mark an atom as having - chiral volume.
       */
      void SetNegativeStereo() { SetFlag(OBAtomFlag::ChiralNeg | OBAtomFlag::Chiral); }
      /** 
       * Clear all stereochemistry information.
       */
      void UnsetStereo();
      /** 
       * Mark an atom as belonging to at least one ring.
       */
      void SetInRing() { SetFlag(OBAtomFlag::Ring); }
      /** 
       * Mark an atom as being chiral with unknown stereochemistry.
       */
      void SetChiral() { SetFlag(OBAtomFlag::Chiral); }
      /** 
       * Clear the internal coordinate pointer.
       */
      void ClearCoordPtr();
      //@}

      //! \name Methods to retrieve atomic information
      //@{
      /** 
       * @return The formal charge for this atom.
       */
      int GetFormalCharge() const;
      /** 
       * @return The atomic number for this atom.
       */
      unsigned int GetAtomicNum() const;
      /** 
       * @return The isotope for this atom, if specified, or 0 for unspecified.
       */
      unsigned short int GetIsotope() const;
      /** 
       * @return The atomic spin, e.g., 0 (default) for normal atoms - note that 
       * this value is a convention, 2 for radical 1 or 3 for carbene.
       */
      int GetSpinMultiplicity() const;
      /** 
       * @return The atomic mass of this atom given by standard IUPAC
       * average molar mass.
       */
      double GetAtomicMass() const;
      /** 
       * @return The atomic mass of given by the isotope
       * (default of 0 gives the most abundant isotope).
       */
      double GetExactMass() const;
      /** 
       * @return The internal atom index (e.g., inside an OBMol).
       */
      unsigned int GetIdx() const { return m_idx; } 
      /** 
       * @return The index into a pointer-driven array as used by
       * GetCoordPtr() or SetCoordPtr().
       */
      unsigned int GetCoordinateIdx() const { return m_cidx; }
      /** 
       * @return The current number of explicit connections.
       */
      unsigned int GetValence() const;
      /** 
       * @return The hybridization of this atom (i.e. 1 for sp, 2 for sp2, 3 for sp3).
       */
      unsigned int GetHyb() const;
      /** 
       * @return The implicit valence of this atom type (i.e. maximum number of connections expected).
       */
      unsigned int GetImplicitValence() const;
      /** 
       * @return The number of non-hydrogens connected to this atom.
       */
      unsigned int GetHvyValence() const;
      /** 
       * @return The number of heteroatoms connected to an atom.
       */
      unsigned int GetHeteroValence() const;
      /** 
       * @return The atomic type (e.g., for molecular mechanics).
       */
      char* GetType();
      /**
       * These methods check to see if there is a coordinate pointer
       * or an internal vector (e.g., SetCoordPtr()).
       * @return The x coordinate.
       */
      double GetX() const { return x(); }
      /**
       * These methods check to see if there is a coordinate pointer
       * or an internal vector (e.g., SetCoordPtr()).
       * @return The y coordinate.
       */
      double GetY() const { return y(); }
      /**
       * These methods check to see if there is a coordinate pointer
       * or an internal vector (e.g., SetCoordPtr()).
       * @return The z coordinate.
       */
      double GetZ() const { return z(); }
      /**
       * These methods check to see if there is a coordinate pointer
       * or an internal vector (e.g., SetCoordPtr()).
       * @return The x coordinate.
       */
      double x() const;
      /**
       * These methods check to see if there is a coordinate pointer
       * or an internal vector (e.g., SetCoordPtr()).
       * @return The y coordinate.
       */
      double y() const;
      /**
       * These methods check to see if there is a coordinate pointer
       * or an internal vector (e.g., SetCoordPtr()).
       * @return The z coordinate.
       */
      double z() const;
      /** 
       * See SetCoordPtr() for more. If no coordinate pointer is used
       * (e.g., only Eigen::Vector3d), NULL will be returned.
       * @return The coordinates as a double* or NULL if none.
       */
      double* GetCoordinate() 
      { 
        if (m_cptr) return(&(*m_cptr)[m_cidx]);
        else     return NULL;
      }
      /** 
       * See SetCoordPtr() for more. If no coordinate pointer is used
       * (e.g., only Eigen::Vector3d), the internal Eigen::Vector3d object 
       * will be returned.
       * @return The coordinates as a Eigen::Vector3d object.
       */
      Eigen::Vector3d& GetVector()
      {
        if (!m_cptr)
          return m_pos;

        m_pos = Eigen::Vector3d( (*m_cptr)[m_cidx], (*m_cptr)[m_cidx+1], (*m_cptr)[m_cidx+2] );
        return m_pos;
      }
      /** 
       * See SetCoordPtr() for more. If no coordinate pointer is used
       * (e.g., only Eigen::Vector3d), the internal Eigen::Vector3d object 
       * will be returned.
       * @return The coordinates as a Eigen::Vector3d object.
       */
      const Eigen::Vector3d& GetVector() const
      {
        if (!m_cptr)
          return m_pos;

        m_pos = Eigen::Vector3d( (*m_cptr)[m_cidx], (*m_cptr)[m_cidx+1], (*m_cptr)[m_cidx+2] );
        return m_pos;
      }
      /** 
       * @return the partial charge of this atom, calculating a Gasteiger charge if needed.
       */
      double GetPartialCharge();
      //! \return the residue which contains this atom, or NULL if none exists.
      OBResidue* GetResidue();
      /** 
       * @param perception implies whether chain perception should occur.
       * @return the residue which contains this atom, or NULL if none exists.
       */
      OBResidue* GetResidue(bool perception);
      /** 
       * @return the molecule which contains this atom, or NULL if none exists.
       */
      OBMol* GetParent();
      /** 
       * @return the OBBond object between this atom and that supplied,
       * or NULL if the two atoms are not bonded.
       */
      OBBond* GetBond(OBAtom *);
      //@}

      //! \name Iterator methods
      //@{
      /** 
       * @return An iterator to the beginning of the bonds to this atom.
       */
      OBBondIterator BeginBonds();
      /** 
       * @return An iterator to the end of the bonds to this atom.
       */
      OBBondIterator EndBonds();
      /** 
       * Set the iterator @p i to the beginning of the bonds.
       * @return The first bond to this atom (or NULL if none exist).
       */
      OBBond* BeginBond(OBBondIterator &i);
      /** 
       * Increment the iterator @p i.
       * @return The next bond to this atom (or NULL if none exist).
       */
      OBBond* NextBond(OBBondIterator &i);
      /** 
       * Set the iterator @p i to the beginning of the bonds.
       * @return The first neighboring atom (or NULL if none exist).
       */
      OBAtom* BeginNbrAtom(OBBondIterator &i);
      /** 
       * Increment the iterator @p i.
       * @return The next neighboring atom (or NULL if none exist).
       */
      OBAtom* NextNbrAtom(OBBondIterator &i);
      //@}

      /** 
       * @return The distance to the atom defined by OBMol::GetAtom().
       */
      double GetDistance(unsigned int index);
      /** 
       * @return the distance to the supplied OBAtom.
       */
      double GetDistance(OBAtom*);
      /** 
       * @return the angle (in degrees) defined by this atom -> b (vertex) -> c.
       */
      double GetAngle(int b, int c);
      /** 
       * @return the angle (in degrees) defined by this atom -> b (vertex) -> c.
       */
      double GetAngle(OBAtom *b, OBAtom *c);

      //! @name Addition of residue/bond info. for an atom
      //@{
      /**
       * If no residue has been set for this atom, create a new one.
       */
      void NewResidue();
      /** 
       * Add (set) the residue for this atom.
       */
      void AddResidue(OBResidue *res);
      /** 
       * Delete any residue associated with this atom.
       */
      void DeleteResidue();
      /** 
       * Add a bond to the internal list. Does not update the bond.
       */
      void AddBond(OBBond *bond);
      /** 
       * Insert @p bond into the internal list at the position from @p i.
       * Does not update the bond.
       */
      void InsertBond(OBBondIterator &i, OBBond *bond);
      /** 
       * Find @p bond and remove it from the internal list. Does not update the bond.
       */
      bool DeleteBond(OBBond* bond);
      /** 
       * Clear all bonding information in this atom (does not delete them).
       */
      void ClearBond();
      //@}

      //! @name Builder utilities
      //@{
      /** 
       * If this is a hydrogen atom, transform into a methyl group.
       * @return Success (true) or failure (flase).
       */
      bool HtoMethyl();
      /** 
       * Change the hybridization of this atom and modify the geometry accordingly
       * @param hyb The hybridization (i.e., 1 = sp, 2 = sp2, 3 = sp3 ...).
       * @return Success (true) or failure (flase).
       */
      bool SetHybAndGeom(int hyb);
      /**
       * Mark that atom has no hydrogens attached.
       */
      void ForceNoH() { SetFlag(OBAtomFlag::HasNoH); }
      /** 
       * @return True if atom has been marked as having no hydrogens attached.
       */
      bool HasNoHForced() { return HasFlag(OBAtomFlag::HasNoH); }
      /** 
       * Mark that atom is not hydrogen deficient (For SMILES input).
       * @since version 2.2
       */
      void ForceImplH() { SetFlag(OBAtomFlag::NotHDeficient); }
      /** 
       * @return if atom has been marked as having no hydrogens attached.
       * @since version 2.2
       */
      bool HasImplHForced() { return HasFlag(OBAtomFlag::NotHDeficient); }
      //@}

      //! \name Property information
      //@{
      /** 
       * @return The number of oxygen atoms connected that only have one heavy valence.
       */
      unsigned int CountFreeOxygens() const;
      /** 
       * @return The number of hydrogens needed to fill the implicit valence of this atom.
       */
      unsigned int ImplicitHydrogenCount() const;
      /** 
       * @return The number of hydrogens explicitly bound to this atom, optionally 
       * excluding D,T and isotope explicitly set to 1.
       */
      unsigned int ExplicitHydrogenCount(bool ExcludeIsotopes = false) const;
      /** 
       * @return The number of rings that contain this atom.
       */
      unsigned int MemberOfRingCount() const;
      /** 
       * @return The size of the smallest ring that contains this atom (0 if not in a ring).
       */
      unsigned int MemberOfRingSize() const;
      /** 
       * @return The number of explicit ring connections to this atom.
       */
      unsigned int CountRingBonds() const;
      /** 
       * @return The smallest angle of bonds to this atom.
       */
      double SmallestBondAngle();
      /** 
       * @return The average angle of bonds to this atom.
       */
      double AverageBondAngle();
      /** 
       * @return The sum of the bond orders of the bonds to the atom (i.e. double bond = 2...).
       */
      unsigned int BOSum() const;
      /** 
       * @return True if there is any residue information.
       */
      bool HasResidue();
      /** 
       * @return Is the atom hydrogen?
       */
      bool IsHydrogen() { return (GetAtomicNum() == 1); }
      /** 
       * @return Is the atom carbon?
       */
      bool IsCarbon() { return (GetAtomicNum() == 6); }
      /** 
       * @return Is the atom nitrogen?
       */
      bool IsNitrogen() { return (GetAtomicNum() == 7); }
      /** 
       * @return Is the atom oxygen?
       */
      bool IsOxygen() { return(GetAtomicNum() == 8); }
      /** 
       * @return Is the atom sulfur?
       */
      bool IsSulfur() { return(GetAtomicNum() == 16); }
      /** 
       * @return Is the atom phosphorus?
       */
      bool IsPhosphorus() { return(GetAtomicNum() == 15); }
      /** 
       * @return Is the atom aromatic?
       */
      bool IsAromatic() const;
      /** 
       * @return Is the atom in a ring?
       */
      bool IsInRing() const;
      /** 
       * @return Is the atom in a ring of a given size?
       */
      bool IsInRingSize(int) const;
      /** 
       * @return Is this atom an element in the 15th or 16th main groups
       * (i.e., N, O, P, S ...) ?
       */
      bool IsHeteroatom();
      /** 
       * @return Is this atom any element except carbon or hydrogen?
       */
      bool IsNotCorH();
      /** 
       * @return Is this atom directly connected to the supplied OBAtom?
       */
      bool IsConnected(OBAtom*);
      /** 
       * @return Is this atom related to the supplied OBAtom in 
       * a 1,3 bonding pattern?
       */
      bool IsOneThree(OBAtom*);
      /** 
       * @return Is this atom related to the supplied OBAtom in
       * a 1,4 bonding pattern?
       */
      bool IsOneFour(OBAtom*);
      /** 
       * @return Is this atom an oxygen in a carboxyl (-CO2 or CO2H) group?
       */
      bool IsCarboxylOxygen();
      /** 
       * @return Is this atom an oxygen in a phosphate (R-PO3) group? 
       */
      bool IsPhosphateOxygen();
      /** 
       * @return Is this atom an oxygen in a sulfate (-SO3) group?
       */
      bool IsSulfateOxygen();
      /** 
       * @return Is this atom an oxygen in a nitro (-NO2) group?
       */
      bool IsNitroOxygen();
      /** 
       * @return Is this atom a nitrogen in an amide (-C(=O)NR2) group?
       */
      bool IsAmideNitrogen();
      /** 
       * @return Is this atom a hydrogen connected to a polar atom
       * (i.e., N, O, P, S)
       */
      bool IsPolarHydrogen();
      /** 
       * @return Is this atom a hydrogen connected to a non-polar atom
       * (i.e., C).
       */
      bool IsNonPolarHydrogen();
      /** 
       * @return Is this atom an aromatic nitrogen with at least one
       * double bond to an oxygen atom.
       */
      bool IsAromaticNOxide();
      /** 
       * @return Is this atom chiral?
       */
      bool IsChiral();
      /** 
       * @return Is this atom an axial atom in a ring?
       */
      bool IsAxial();
      /** 
       * @return Does this atom have SMILES-specified clockwise "@@" stereochemistry?
       */
      bool IsClockwise() { return(HasFlag(OBAtomFlag::StereoC));  }
      /** 
       * @return Does this atom have SMILES-specified anticlockwise "@" stereochemistry?
       */
      bool IsAntiClockwise()  { return(HasFlag(OBAtomFlag::StereoAC)); }
      /** 
       * @return Does this atom have a positive chiral volume?
       */
      bool IsPositiveStereo() { return(HasFlag(OBAtomFlag::ChiralPos)); }
      /** 
       * @return Does this atom have a negative chiral volume?
       */
      bool IsNegativeStereo() { return(HasFlag(OBAtomFlag::ChiralNeg)); }
      /** 
       * @return Does this atom have SMILES-specified stereochemistry?
       */
      bool HasChiralitySpecified()
      { return(HasFlag(OBAtomFlag::StereoC | OBAtomFlag::StereoAC)); }
      /** 
       * @return Does this atom have a specified chiral volume?
       */
      bool HasChiralVolume()
      { return(HasFlag(OBAtomFlag::ChiralPos | OBAtomFlag::ChiralNeg)); }
      /** 
       * @return Is this atom a hydrogen-bond acceptor (receptor)?
       */
      bool IsHbondAcceptor();
      /** 
       * @return Is this atom a hydrogen-bond donor?
       */
      bool IsHbondDonor();
      /** 
       * @return Is this a hydrogen atom attached to a hydrogen-bond donor?
       */
      bool IsHbondDonorH();
      /**
       * @param includePandS Whether to include phosphorus and sulfur neighbors
       * in this determination (or to exclude them).
       * @return Whether a neighboring atom (alpha) has an unsaturated bond
       * to a third atom (beta).
       */
      bool HasAlphaBetaUnsat(bool includePandS = true);
      /** 
       * @return Whether this atom is connected to any bond with order == @p bo.
       */
      bool HasBondOfOrder(unsigned int bo);
      /** 
       * @return The count of bonds connected to this atom with order == @p bo.
       */
      int CountBondsOfOrder(unsigned int bo);
      /** 
       * @return Whether this atom is connected to any bond with order >1.
       */
      bool HasNonSingleBond();
      /** 
       * @return Does this atom have a single bond.
       */
      bool HasSingleBond()   { return HasBondOfOrder(1); }
      /** 
       * @return Does this atom have a double bond.
       */
      bool HasDoubleBond()   { return HasBondOfOrder(2); }
      /** 
       * @return Does this atom have an aromatic bond.
       */
      bool HasAromaticBond() { return HasBondOfOrder(5); }
      /** 
       * @return Whether this atom matches the first atom in a given SMARTS pattern.
       */
      bool MatchesSMARTS(const char *);
      //@}
      
      //! @todo remove this functions...
      unsigned int KBOSum() const;
  
  }; // class OBAtom
  
  //@} group

}// namespace OpenBabel

#endif   // OB_ATOM_H

//! @file atom.h
//! @brief Handle atoms

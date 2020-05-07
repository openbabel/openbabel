/**********************************************************************
residue.h - Defines for residue properties, names, etc.

Copyright (C) 2001, 2002  OpenEye Scientific Software, Inc.
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

/**********************************************************************
Global arrays Residue, ElemDesc and function GetResidueNumber were
obtained in part or whole from RasMol2 by Roger Sayle.
***********************************************************************/

#ifndef OB_RESIDUE_H
#define OB_RESIDUE_H

#include <openbabel/babelconfig.h>

#ifndef OB_EXTERN
#  define OB_EXTERN extern
#endif

#include <vector>
#include <string>

#include <openbabel/base.h>

namespace OpenBabel {

  class OBAtom;
  //! A standard iterator over a vector of atoms
  typedef std::vector<OBAtom*>::iterator OBAtomIterator;
  class OBBond;
  //! A standard iterator over a vector of bonds
  typedef std::vector<OBBond*>::iterator OBBondIterator;

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

    //! Add @p atom to this residue. Updates the atom via OBAtom::SetResidue()
    void    AddAtom(OBAtom *atom);
    //! Add @p atom to this residue. Updates the atom via OBAtom::SetResidue()
    void    InsertAtom(OBAtom *atom);
    //! Remove @p atom from this residue and update the atom.
    void    RemoveAtom(OBAtom *atom);
    //! Clear any and all data associated with this residue. Updates all atoms
    //!  included in the residue, as well as calling OBBase::Clear() for any
    //!  generic data.
    //! \return Whether the call was successful.
    bool    Clear();

    //! \brief Set the name of this residue (e.g., "ALA"). Use 3-char PDB standard names.
    //! http://www.rcsb.org/pdb/file_formats/pdb/pdbguide2.2/part_79.html
    //! MODRES records for modified residues:
    //! http://www.rcsb.org/pdb/file_formats/pdb/pdbguide2.2/part_36.html
    void    SetName(const std::string &resname);
    //! Set the residue number (in the sequence)
    void    SetNum(const unsigned int resnum);
    void    SetNum(const std::string  resnum);
    //! Set the chain ID for this residue
    void    SetChain(const char chain);
    //! Set the chain number for this residue
    void    SetChainNum(const unsigned int chainnum);
    //! Set the internal index of this residue in the parent OBMol.
    //! Intended mostly for internal use
    void    SetIdx(const unsigned int idx);
    //! Set  PDB insertion code information for this residue. This allows
    //! consecutive residues to have the same number. Some communities
    //! that work in a well-conserved structural world use this, e.g.
    //! for immunoglobulins.
    void    SetInsertionCode(const char insertioncode);

    //! Set the character code ID for an ATOM record for the supplied atom
    //! This does nothing if the supplied atom is not found in the residue
    void    SetAtomID(OBAtom *atom, const std::string &id);
    void    SetHetAtom(OBAtom *atom, bool hetatm);
    //! Set the atomic serial number for a given atom (see OBSerialNums)
    void    SetSerialNum(OBAtom *atom, unsigned int sernum);

    //! \return The residue name
    std::string    GetName(void)                  const;
    //! \return The residue number (in the sequence)
    int    GetNum(void);
    std::string     GetNumString(void);
    //! \return The number of atoms in this residue
    unsigned int   GetNumAtoms()                  const;
    //! \return The ID of the chain which includes this residue
    char           GetChain(void)                 const;
    //! \return The number of the chain which includes this residue
    unsigned int   GetChainNum(void)              const;
    //! \return The internal index of this residue in the parent OBMol
    unsigned int   GetIdx(void)                   const;
    //! \return The residue key (i.e., an entry in the OBResidueIndex namespace)
    unsigned int   GetResKey(void)                const;

    //! \return a vector of all atoms in this residue
    std::vector<OBAtom*> GetAtoms(void)           const;
    //! \return all bonds in this residue. @p exterior includes bonds to atoms
    //!  outside this residue (default is true)
    std::vector<OBBond*> GetBonds(bool exterior= true)const;

    //! \return the atom ID (character code) for the supplied atom or ""
    //!  if the atom is not found in this residue
    std::string    GetAtomID(OBAtom *atom)        const;
    //! \return the serial number of the supplied atom (uses OBSerialNums)
    unsigned       GetSerialNum(OBAtom *atom)     const;
    //! \return The Insertion Code (i.e., an extra position motivated by a
    //! multiple sequence alignment against a template with defined numbers)
    char           GetInsertionCode(void)	  const;

    //! \return Whether this residue has the supplied amino acid property
    //!  defined from the OBAminoAcidProperty namespace
    bool           GetAminoAcidProperty(int)      const;
    //! \return Whether atom @p a has the supplied residue atom property
    //!  defined from the OBResidueAtomProperty namespace
    bool           GetAtomProperty(OBAtom *a, int) const;
    //! \return Whether this residue has the supplied property
    //!  defined from the OBResidueProperty namespace
    bool           GetResidueProperty(int)        const;

    //! \return If the given atom is a HETATM record
    bool           IsHetAtom(OBAtom *atom)        const;
    //! \return If this residue matches the supplied @p restype
    //! Set by SetResidueKeys()
    bool           IsResidueType(int)             const;

    //! \name Iterator methods
    //@{
    //! \return An iterator to the beginning of the atom list in this residue
    OBAtomIterator BeginAtoms()   { return _atoms.begin(); }
    //! \return An iterator to the end of the atom list in this residue
    OBAtomIterator EndAtoms()     { return _atoms.end();   }
    //! Set the iterator @p i to the beginning of the atom list in this residue
    //! \return The first atom (or NULL if none exist)
    OBAtom *BeginAtom(std::vector<OBAtom*>::iterator &i);
    //! Increment the iterator @p i
    //! \return The next atom (or NULL if none exist)
    OBAtom *NextAtom(std::vector<OBAtom*>::iterator &i);
    //@}

  protected: // members

    unsigned int              _idx;   //!< Residue index (i.e., internal index in an OBMol)
    char                      _chain; //!< Chain ID
    unsigned int              _aakey; //!< Amino Acid key ID -- see SetResidueKeys()
    unsigned int              _reskey;//!< Residue key ID -- see SetResidueKeys()
    std::string               _resnum;//!< Residue number (i.e., in file) 23, 1B, etc.
    std::string               _resname;//!<Residue text name
    char                _insertioncode;//!<PBB insertion code

    std::vector<bool>         _hetatm;//!< Is a given atom a HETAM
    std::vector<std::string>  _atomid;//!< Residue atom text IDs
    std::vector<OBAtom*>      _atoms; //!< List of OBAtom in this residue
    std::vector<unsigned int> _sernum;//!< List of serial numbers
    // Now in OBBase
    //    std::vector<OBGenericData*> _vdata; //!< Custom data
  }; // OBResidue


  ///////////////////////////////////////////////////////////////////////////////
  // Global Definitions
  ///////////////////////////////////////////////////////////////////////////////

#define MAXSETNO 40
#define MAXELEM  29
#define MAXRES   54

  ///////////////////////////////////////////////////////////////////////////////
  // Amino Acid Definitions
  ///////////////////////////////////////////////////////////////////////////////

#define AA_ALA (1<<1)
#define AA_GLY (1<<2)
#define AA_LEU (1<<3)
#define AA_SER (1<<4)
#define AA_VAL (1<<5)
#define AA_THR (1<<6)
#define AA_LYS (1<<7)
#define AA_ASP (1<<8)
#define AA_ILE (1<<9)
#define AA_ASN (1<<10)
#define AA_GLU (1<<11)
#define AA_PRO (1<<12)
#define AA_ARG (1<<13)
#define AA_PHE (1<<14)
#define AA_GLN (1<<15)
#define AA_TYR (1<<16)
#define AA_HIS (1<<17)
#define AA_CYS (1<<18)
#define AA_MET (1<<19)
#define AA_TRP (1<<20)

  /////////////////////////////////////////////////////////////////////////////
  // Amino Acid Property Definitions
  /////////////////////////////////////////////////////////////////////////////
#define IS_ACIDIC(x)      ((x) & ((AA_ASP)|(AA_GLU)))
#define IS_ACYCLIC(x)     ((x) & ((AA_ALA)|(AA_GLY)|(AA_LEU)|(AA_SER)|  \
                                  (AA_VAL)|(AA_THR)|(AA_LYS)|(AA_ASP)|  \
                                  (AA_ILE)|(AA_ASN)|(AA_GLU)|(AA_GLN)|  \
                                  (AA_CYS)|(AA_MET)))
#define IS_ALIPHATIC(x)   ((x) & ((AA_ALA)|(AA_GLY)|(AA_ILE)|(AA_LEU)|  \
                                  (AA_VAL)))
#define IS_AROMATIC(x)    ((x) & ((AA_HIS)|(AA_PHE)|(AA_TRP)|(AA_TYR)))
#define IS_BASIC(x)       ((x) & ((AA_ARG)|(AA_HIS)|(AA_LYS)))
#define IS_BURIED(x)      ((x) & ((AA_ALA)|(AA_CYS)|(AA_ILE)|(AA_LEU)|  \
                                  (AA_MET)|(AA_PHE)|(AA_TRP)|(AA_VAL)))
#define IS_CHARGED(x)     ((x) & ((AA_ASP)|(AA_GLU)|(AA_ARG)|(AA_HIS)|  \
                                  (AA_LYS)))
#define IS_CYCLIC(x)      ((x) & ((AA_HIS)|(AA_PHE)|(AA_PRO)|(AA_TRP)|  \
                                  (AA_TYR)))
#define IS_HYDROPHOBIC(x) ((x) & ((AA_ALA)|(AA_LEU)|(AA_VAL)|(AA_ILE)|  \
                                  (AA_PRO)|(AA_PHE)|(AA_MET)|(AA_TRP)))
#define IS_LARGE(x)       ((x) & ((AA_ARG)|(AA_PHE)|(AA_GLN)|(AA_TYR)|  \
                                  (AA_HIS)|(AA_LEU)|(AA_LYS)|(AA_ILE)|  \
                                  (AA_GLU)|(AA_MET)|(AA_TRP)))
#define IS_MEDIUM(x)      ((x) & ((AA_VAL)|(AA_THR)|(AA_ASP)|(AA_ASN)|  \
                                  (AA_PRO)|(AA_CYS)))
#define IS_NEGATIVE(x)    ((x) & ((AA_ASP)|(AA_GLU)))
#define IS_NEUTRAL(x)     ((x) & ((AA_ALA)|(AA_GLY)|(AA_LEU)|(AA_SER)|  \
                                  (AA_VAL)|(AA_THR)|(AA_PHE)|(AA_GLN)|  \
                                  (AA_TYR)|(AA_HIS)|(AA_CYS)|(AA_MET)|  \
                                  (AA_TRP)|(AA_ILE)|(AA_ASN)|(AA_PRO)))
#define IS_POLAR(x)       ((x) & ((AA_ASP)|(AA_ILE)|(AA_ASN)|(AA_GLU)|  \
                                  (AA_SER)|(AA_THR)|(AA_ARG)|(AA_GLN)|  \
                                  (AA_CYS)|(AA_HIS)))
#define IS_POSITIVE(x)    ((x) & ((AA_ARG)|(AA_HIS)|(AA_LYS)))
#define IS_SMALL(x)       ((x) & ((AA_ALA)|(AA_GLY)|(AA_SER)))
#define IS_SURFACE(x)     ((x) & ((AA_THR)|(AA_LYS)|(AA_ASP)|(AA_ILE)|  \
                                  (AA_ASN)|(AA_GLU)|(AA_PRO)|(AA_ARG)|  \
                                  (AA_GLY)|(AA_SER)|(AA_GLN)|(AA_TYR)|  \
                                  (AA_HIS)))

  //! Residue property definitions
  namespace OBAminoAcidProperty
  {
	enum
    {
      ACIDIC      =  0,
      ACYCLIC     =  1,
      ALIPHATIC   =  2,
      AROMATIC    =  3,
      BASIC       =  4,
      BURIED      =  5,
      CHARGED     =  6,
      CYCLIC      =  7,
      HYDROPHOBIC =  8,
      LARGE       =  9,
      MEDIUM      = 10,
      NEGATIVE    = 11,
      NEUTRAL     = 12,
      POLAR       = 13,
      POSITIVE    = 14,
      SMALL       = 15,
      SURFACE     = 16
    };
  }

  //! Residue atom properties
  namespace OBResidueAtomProperty
  {
    enum
    {
      ALPHA_CARBON     = 0,
      AMINO_BACKBONE   = 1,
      BACKBONE         = 2,
      CYSTEINE_SULPHUR = 3,
      LIGAND           = 4,
      NUCLEIC_BACKBONE = 5,
      SHAPELY_BACKBONE = 6,
      SHAPELY_SPECIAL  = 7,
      SIDECHAIN        = 8,
      SUGAR_PHOSPHATE  = 9
    };
  }

  //! Residue names (index into Residue[] array)
  // some of these are invalid or troublesome in scripting interfaces
  // so they are removed by the #ifndef SWIG parts
  // (otherwise ignore them for C++ use)
  namespace OBResidueIndex
  {
    
    enum
    {
      ALA   =  0,
      GLY   =  1,
      LEU   =  2,
      SER   =  3,
      VAL   =  4,
#ifndef SWIGPERL
      THR   =  5,
#endif
      LYS   =  6,
      ASP   =  7,
      ILE   =  8,
      ASN   =  9,
      GLU   = 10,
      PRO   = 11,
      ARG   = 12,
      PHE   = 13,
      GLN   = 14,
      TYR   = 15,
      HIS   = 16,
      CYS   = 17,
      MET   = 18,
      TRP   = 19,
      ASX   = 20,
      GLX   = 21,
      PCA   = 22,
      HYP   = 23,
      A     = 24,
      C     = 25,
      G     = 26,
      T     = 27,
      U     = 28,
      UPLUS = 29,
      I     = 30,
      _1MA  = 32,
      _5MC  = 32,
      OMC   = 33,
      _1MG  = 34,
      _2MG  = 35,
      M2G   = 36,
      _7MG  = 37,
      OMG   = 38,
      YG    = 39,
      H2U   = 40,
      _5MU  = 41,
      PSU   = 42,
      UNK   = 43,
      ACE   = 44,
      FOR   = 45,
      HOH   = 46,
      DOD   = 47,
      SO4   = 48,
      PO4   = 49,
      NAD   = 50,
      COA   = 51,
      NAP   = 52,
      NDP   = 53
    };
  }

  //! Residue types.
  namespace OBResidueProperty
  {
    enum
     {
      AMINO        = 0,
      AMINO_NUCLEO = 1,
      COENZYME     = 2,
      ION          = 3,
      NUCLEO       = 4,
      PROTEIN      = 5,
      PURINE       = 6,
      PYRIMIDINE   = 7,
      SOLVENT      = 8,
      WATER        = 9
    };
  }

} // end namespace OpenBabel

#endif

//! \file residue.h
//! \brief Defines for residue properties, names, etc.

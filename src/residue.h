/**********************************************************************
residue.h - Defines for residue properties, names, etc.
 
Copyright (C) 2001, 2002  OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2006 by Geoffrey R. Hutchison
 
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

/**********************************************************************
Global arrays Residue, ElemDesc and function GetResidueNumber were
obtained in part or whole from RasMol2 by Roger Sayle.
***********************************************************************/

#ifndef OB_RESIDUE_H
#define OB_RESIDUE_H

#include "babelconfig.h"

#include <vector>
#include <string>

#include "base.h"

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

    OBAtomIterator BeginAtoms()   { return _atoms.begin(); }
    OBAtomIterator EndAtoms()     { return _atoms.end();   }
    OBAtom *BeginAtom(std::vector<OBAtom*>::iterator &i);
    OBAtom *NextAtom(std::vector<OBAtom*>::iterator &i);

  protected: // members

    unsigned int              _idx;   //!< Residue index (i.e., internal index in an OBMol)
    char                      _chain; //!< Chain ID
    unsigned int              _aakey; //!< Amino Acid key ID -- see SetResidueKeys()
    unsigned int              _reskey;//!< Residue key ID -- see SetResidueKeys()
    unsigned int              _resnum;//!< Residue number (i.e., in file)
    std::string               _resname;//!< Residue text name

    std::vector<bool>         _hetatm;//!< Is a given atom a HETAM
    std::vector<std::string>  _atomid;//!< Residue atom text IDs
    std::vector<OBAtom*>      _atoms; //!< List of OBAtom in this residue
    std::vector<unsigned int> _sernum;//!< List of serial numbers
    //    std::vector<OBGenericData*> _vdata; //!< Custom data
  }; // OBResidue

  //! A standard iterator over a vector of residues
  typedef std::vector<OBResidue*>::iterator OBResidueIterator;

  ///////////////////////////////////////////////////////////////////////////////
  // Global Definitions
  ///////////////////////////////////////////////////////////////////////////////

#define MAXSETNO 40
#define MAXELEM  1024
#define MINELEM  29
#define MAXRES   100
#define MINRES   54

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
    static const unsigned int ACIDIC      =  0;
    static const unsigned int ACYCLIC     =  1;
    static const unsigned int ALIPHATIC   =  2;
    static const unsigned int AROMATIC    =  3;
    static const unsigned int BASIC       =  4;
    static const unsigned int BURIED      =  5;
    static const unsigned int CHARGED     =  6;
    static const unsigned int CYCLIC      =  7;
    static const unsigned int HYDROPHOBIC =  8;
    static const unsigned int LARGE       =  9;
    static const unsigned int MEDIUM      = 10;
    static const unsigned int NEGATIVE    = 11;
    static const unsigned int NEUTRAL     = 12;
    static const unsigned int POLAR       = 13;
    static const unsigned int POSITIVE    = 14;
    static const unsigned int SMALL       = 15;
    static const unsigned int SURFACE     = 16;
  }

  //! Residue atom properties
  namespace OBResidueAtomProperty
  {
    static const unsigned int ALPHA_CARBON     = 0;
    static const unsigned int AMINO_BACKBONE   = 1;
    static const unsigned int BACKBONE         = 2;
    static const unsigned int CYSTEINE_SULPHUR = 3;
    static const unsigned int LIGAND           = 4;
    static const unsigned int NUCLEIC_BACKBONE = 5;
    static const unsigned int SHAPELY_BACKBONE = 6;
    static const unsigned int SHAPELY_SPECIAL  = 7;
    static const unsigned int SIDECHAIN        = 8;
    static const unsigned int SUGAR_PHOSPHATE  = 9;
  }

  //! Residue names
  namespace OBResidueIndex
  {
    static const unsigned int ALA   =  0;
    static const unsigned int GLY   =  1;
    static const unsigned int LEU   =  2;
    static const unsigned int SER   =  3;
    static const unsigned int VAL   =  4;
#ifndef SWIGPERL
    static const unsigned int THR   =  5;
#endif
    static const unsigned int LYS   =  6;
    static const unsigned int ASP   =  7;
    static const unsigned int ILE   =  8;
    static const unsigned int ASN   =  9;
    static const unsigned int GLU   = 10;
    static const unsigned int PRO   = 11;
    static const unsigned int ARG   = 12;
    static const unsigned int PHE   = 13;
    static const unsigned int GLN   = 14;
    static const unsigned int TYR   = 15;
    static const unsigned int HIS   = 16;
    static const unsigned int CYS   = 17;
    static const unsigned int MET   = 18;
    static const unsigned int TRP   = 19;
    static const unsigned int ASX   = 20;
    static const unsigned int GLX   = 21;
    static const unsigned int PCA   = 22;
    static const unsigned int HYP   = 23;
    static const unsigned int A     = 24;
    static const unsigned int C     = 25;
    static const unsigned int G     = 26;
    static const unsigned int T     = 27;
    static const unsigned int U     = 28;
    static const unsigned int UPLUS = 29;
    static const unsigned int I     = 30;
    static const unsigned int _1MA  = 31;
    static const unsigned int _5MC  = 32;
    static const unsigned int OMC   = 33;
    static const unsigned int _1MG  = 34;
    static const unsigned int _2MG  = 35;
    static const unsigned int M2G   = 36;
    static const unsigned int _7MG  = 37;
    static const unsigned int OMG   = 38;
    static const unsigned int YG    = 39;
    static const unsigned int H2U   = 40;
    static const unsigned int _5MU  = 41;
    static const unsigned int PSU   = 42;
    static const unsigned int UNK   = 43;
    static const unsigned int ACE   = 44;
    static const unsigned int FOR   = 45;
    static const unsigned int HOH   = 46;
    static const unsigned int DOD   = 47;
    static const unsigned int SO4   = 48;
    static const unsigned int PO4   = 49;
    static const unsigned int NAD   = 50;
    static const unsigned int COA   = 51;
    static const unsigned int NAP   = 52;
    static const unsigned int NDP   = 53;
  }

  //! Residue types.
  namespace OBResidueProperty
  {
    static const unsigned int AMINO        = 0;
    static const unsigned int AMINO_NUCLEO = 1;
    static const unsigned int COENZYME     = 2;
    static const unsigned int ION          = 3;
    static const unsigned int NUCLEO       = 4;
    static const unsigned int PROTEIN      = 5;
    static const unsigned int PURINE       = 6;
    static const unsigned int PYRIMIDINE   = 7;
    static const unsigned int SOLVENT      = 8;
    static const unsigned int WATER        = 9;
  }

  ////////////////////////////////////////////////////////////////////////////////
  // Global Variables
  ////////////////////////////////////////////////////////////////////////////////

  static char Residue[MAXRES][4] = {
    /*===============*/
    /*  Amino Acids  */
    /*===============*/

    /* Ordered by Cumulative Frequency in Brookhaven *
     * Protein Databank, December 1991               */

    "ALA", /* 8.4% */     "GLY", /* 8.3% */
    "LEU", /* 8.0% */     "SER", /* 7.5% */
    "VAL", /* 7.1% */     "THR", /* 6.4% */
    "LYS", /* 5.8% */     "ASP", /* 5.5% */
    "ILE", /* 5.2% */     "ASN", /* 4.9% */
    "GLU", /* 4.9% */     "PRO", /* 4.4% */
    "ARG", /* 3.8% */     "PHE", /* 3.7% */
    "GLN", /* 3.5% */     "TYR", /* 3.5% */
    "HIS", /* 2.3% */     "CYS", /* 2.0% */
    "MET", /* 1.8% */     "TRP", /* 1.4% */

    "ASX", "GLX", "PCA", "HYP",

    /*===================*/
    /*  DNA Nucleotides  */
    /*===================*/
    "  A", "  C", "  G", "  T",

    /*===================*/
    /*  RNA Nucleotides  */
    /*===================*/
    "  U", " +U", "  I", "1MA",
    "5MC", "OMC", "1MG", "2MG",
    "M2G", "7MG", "OMG", " YG",
    "H2U", "5MU", "PSU",

    /*=================*/
    /*  Miscellaneous  */
    /*=================*/
    "UNK", "ACE", "FOR", "HOH",
    "DOD", "SO4", "PO4", "NAD",
    "COA", "NAP", "NDP"
  };

  /* Avoid SGI Compiler Warnings! */
  static char ElemDesc[MAXELEM][4] = {
    { ' ', 'N', ' ', ' ' },  /* 0*/
    { ' ', 'C', 'A', ' ' },  /* 1*/
    { ' ', 'C', ' ', ' ' },  /* 2*/
    { ' ', 'O', ' ', ' ' },  /* 3*/   /* 0-3   Amino Acid Backbone    */
    { ' ', 'C', '\'', ' ' }, /* 4*/
    { ' ', 'O', 'T', ' ' },  /* 5*/
    { ' ', 'S', ' ', ' ' },  /* 6*/
    { ' ', 'P', ' ', ' ' },  /* 7*/   /* 4-7   Shapely Amino Backbone */
    { ' ', 'O', '1', 'P' },  /* 8*/
    { ' ', 'O', '2', 'P' },  /* 9*/
    { ' ', 'O', '5', '*' },  /*10*/
    { ' ', 'C', '5', '*' },  /*11*/
    { ' ', 'C', '4', '*' },  /*12*/
    { ' ', 'O', '4', '*' },  /*13*/
    { ' ', 'C', '3', '*' },  /*14*/
    { ' ', 'O', '3', '*' },  /*15*/
    { ' ', 'C', '2', '*' },  /*16*/
    { ' ', 'O', '2', '*' },  /*17*/
    { ' ', 'C', '1', '*' },  /*18*/   /* 7-18  Nucleic Acid Backbone  */
    { ' ', 'C', 'A', '2' },  /*19*/   /* 19    Shapely Special        */
    { ' ', 'S', 'G', ' ' },  /*20*/   /* 20    Cysteine Sulphur       */
    { ' ', 'N', '1', ' ' },  /*21*/
    { ' ', 'N', '2', ' ' },  /*22*/
    { ' ', 'N', '3', ' ' },  /*23*/
    { ' ', 'N', '4', ' ' },  /*24*/
    { ' ', 'N', '6', ' ' },  /*25*/
    { ' ', 'O', '2', ' ' },  /*26*/
    { ' ', 'O', '4', ' ' },  /*27*/
    { ' ', 'O', '6', ' ' }   /*28*/   /* 21-28 Nucleic Acid H-Bonding */
  };

  static unsigned int ResNo  = MINRES;
  static unsigned int ElemNo = MINELEM;

} // end namespace OpenBabel

#endif

//! \file residue.h
//! \brief Defines for residue properties, names, etc.

/*
 * International Chemical Identifier (InChI)
 * Version 1
 * Software version 1.07
 * April 30, 2024
 *
 * MIT License
 *
 * Copyright (c) 2024 IUPAC and InChI Trust
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 * The InChI library and programs are free software developed under the
 * auspices of the International Union of Pure and Applied Chemistry (IUPAC).
 * Originally developed at NIST.
 * Modifications and additions by IUPAC and the InChI Trust.
 * Some portions of code were developed/changed by external contributors
 * (either contractor or volunteer) which are listed in the file
 * 'External-contributors' included in this distribution.
 *
 * info@inchi-trust.org
 *
 */

#ifndef _INPDEF_H_
#define _INPDEF_H_

/* input/output format */

#include "ichidrp.h"
#include "mode.h"
#include "mol_fmt.h"

#define CLOSING_STARRED_SRU_IS_A_MUST 1
#define ALLOW_CLOSING_SRU_VIA_HIGHER_ORDER_BOND 1
#define ALLOW_CLOSING_SRU_VIA_DIRADICAL 1

#define CLOSING_SRU_NOT_APPLICABLE 0
#define CLOSING_SRU_RING 1
#define CLOSING_SRU_HIGHER_ORDER_BOND 2
#define CLOSING_SRU_DIRADICAL 3
#define CLOSING_SRU_RING_OPENED 4
#define CLOSING_SRU_MOVED_BRACKETS 5

#define bDrawingLabelLeftShift endpoint /**< for drawing only */
typedef S_SHORT ST_CAP_FLOW;

/* inp_ATOM::at_type */
#define ATT_NONE 0x0000
#define ATT_ACIDIC_CO 0x0001
#define ATT_ACIDIC_S 0x0002
#define ATT_OO 0x0004
#define ATT_ZOO 0x0008
#define ATT_NO 0x0010
#define ATT_N_O 0x0020
#define ATT_ATOM_N 0x0040
#define ATT_ATOM_P 0x0080
#define ATT_OTHER_NEG_O 0x0100
#define ATT_OTHER_ZO 0x0200 /**<  -Z=O or =Z=O */
#define ATT_OH_MINUS 0x0400 /**< OH(-), O=O,S,Se,Te */
#define ATT_O_PLUS 0x0800   /**< -OH2(+), =OH(+), -OH(+)-, OH3(+), =O(+)-, etc; O=O,S,Se,Te */
#define ATT_PROTON 0x1000
#define ATT_HalAnion 0x2000
#define ATT_HalAcid 0x4000
#if (FIX_NP_MINUS_BUG == 1)
#define ATT_NP_MINUS_V23 0x8000 /**< =N(-) or =P(-) where = previously was triple */
#endif

#define AT_FLAG_ISO_H_POINT 0x01 /**< may have isotopic H */

#define PERIODIC_NUMBER_H 1

#ifndef NUMH
#define NUM_ISO_H(AT, N) (AT[N].num_iso_H[0] + AT[N].num_iso_H[1] + AT[N].num_iso_H[2])
#define NUMH(AT, N) (AT[N].num_H + NUM_ISO_H(AT, N))
#endif

#define FlagSC_0D 1 /**< bUsed0DParity */
#define FlagSB_0D 2 /**< bUsed0DParity */

#define SB_PARITY_FLAG 0x38 /**< mask for disconnected metal parity if it is different */
#define SB_PARITY_SHFT 3    /**< number of right shift bits to get disconnected metal parity */
#define SB_PARITY_MASK 0x07
#define SB_PARITY_1(X) (X & SB_PARITY_MASK)                       /**< refers to connected structure */
#define SB_PARITY_2(X) (((X) >> SB_PARITY_SHFT) & SB_PARITY_MASK) /**< refers to connected structure */

/**
 * @brief Structure describing an input atom
 *
 * @param elname        Name of chemical element
 * @param el_number     Number of the element in the Periodic Table
 * @param neighbor      Positions (from 0) of the neighbors in the inp_ATOM array
 * @param orig_at_number  Original atom number
 * @param orig_compt_at_numb Original atom number within the component before terminal H removal
 * @param bond_stereo   Bond stereo information (1=Up,4=Either,6=Down; this atom is at the pointing wedge, negative => on the opposite side; 3=Either double bond)
 * @param bond_type     Bond type (1..4; 4="aromatic", should be discouraged on input)
 * @param valence       Valence (for most chemists, it is coordination number, CN; number of bonds = number of neighbors)
 * @param chem_bonds_valence Chemical bonds valence (for most chemists, it is what usually called valence; sum of bond types (type 4 needs special treatment))
 * @param num_H         Number of implicit hydrogens, including D and T
 * @param num_iso_H     Number of implicit 1H, 2H(D), 3H(T) < 16
 * @param iso_atw_diff  Isotopic atomic weight difference (=0 => natural isotopic abundances; >0 => (mass) - (mass of the most abundant isotope) + 1; <0 => (mass) - (mass of the most abundant
 * isotope))
 * @param charge        Charge
 * @param radical       Radical type (RADICAL_SINGLET, RADICAL_DOUBLET, or RADICAL_TRIPLET)
 * @param bAmbiguousStereo  Flag indicating ambiguous stereochemistry
 * @param cFlags        Additional flags (e.g., AT_FLAG_ISO_H_POINT)
 * @param at_type       Atom type (ATT_NONE, ATT_ACIDIC, etc.)
 * @param component     Number of the structure component > 0
 * @param endpoint      ID of a tautomeric group
 * @param c_point       ID of a positive charge group
 * @param x, y, z      3D coordinates
 * @param bUsed0DParity Bit flags for 0D parities (bit=1 => stereobond; bit=2 => stereocenter)
 * @param p_parity      0D tetrahedral parity
 * @param p_orig_at_num Original atom numbers of neighbors for tetrahedral stereochemistry
 * @param sb_ord        Stereo bond/neighbor ordering number, starts from 0
 * @param sn_ord        Ord. num. of the neighbor adjacent to the SB; starts from 0; -1 means removed explicit H
 * @param sb_parity    Stereo bond parities
 * @param sn_orig_at_num Original atom numbers of neighbors for stereobonds
 * @param bCutVertex    Flag indicating if the atom is a cut vertex (part of ring system analysis)
 * @param nRingSystem  Ring system number (part of ring system analysis)
 * @param nNumAtInRingSystem Number of atoms in the ring system (part of ring system analysis)
 * @param nBlockSystem Block system number (part of ring system analysis)
 * @param nDistanceFromTerminal Distance from terminal atom or ring system (terminal atom or ring system has 1, next has 2, etc.)
 */
typedef struct tagInputAtom
{
    char elname[ATOM_EL_LEN];         /* name of chemical element                                 */
    U_CHAR el_number;                 /* number of the element in the Periodic Table              */
    AT_NUMB neighbor[MAXVAL];         /* positions (from 0) of the neighbors in the inp_ATOM array*/
    AT_NUMB orig_at_number;           /* original atom number                                     */
    AT_NUMB orig_compt_at_numb;       /* original atom number within the component                */
                                      /* before terminal H removal                                */
    S_CHAR bond_stereo[MAXVAL];       /* 1=Up,4=Either,6=Down; this atom is at the pointing wedge,*/
                                      /*   negative => on the opposite side; 3=Either double bond */
    U_CHAR bond_type[MAXVAL];         /* 1..4; 4="aromatic", should be discouraged on input       */
    S_CHAR valence;                   /* for most chemists, it is coordination number, CN;        */
                                      /* number of bonds = number of neighbors					*/
    S_CHAR chem_bonds_valence;        /* for most chemists, it is what usually called valence;    */
                                      /* sum of bond types (type 4 needs special treatment)       */
    S_CHAR num_H;                     /* number of implicit hydrogens, including D and T          */
    S_CHAR num_iso_H[NUM_H_ISOTOPES]; /* number of implicit 1H, 2H(D), 3H(T) < 16                */
    S_CHAR iso_atw_diff;              /* =0 => natural isotopic abundances                        */
                                      /* >0 => (mass) - (mass of the most abundant isotope) + 1   */
                                      /* <0 => (mass) - (mass of the most abundant isotope)       */
    S_CHAR charge;                    /* charge                                                   */
    S_CHAR radical;                   /* RADICAL_SINGLET, RADICAL_DOUBLET, or RADICAL_TRIPLET     */
    S_CHAR bAmbiguousStereo;
    S_CHAR cFlags;     /* AT_FLAG_ISO_H_POINT                                      */
    AT_NUMB at_type;   /* ATT_NONE, ATT_ACIDIC                                     */
    AT_NUMB component; /* number of the structure component > 0                    */
    AT_NUMB endpoint;  /* id of a tautomeric group                                 */
    AT_NUMB c_point;   /* id of a positive charge group                            */
    double x;
    double y;
    double z;
    /* 0D parities (originally were used with CML) */
    S_CHAR bUsed0DParity; /* bit=1 => stereobond; bit=2 => stereocenter               */
    /* 0D tetrahedral parity */
    S_CHAR p_parity;
    AT_NUMB p_orig_at_num[MAX_NUM_STEREO_ATOM_NEIGH];
    /* 0D bond parities */
    S_CHAR sb_ord[MAX_NUM_STEREO_BONDS]; /* stereo bond/neighbor ordering number, starts from 0  */
    /* neighbors on both sides of stereobond have same sign=> trans/T/E, diff. signs => cis/C/Z         */
    S_CHAR sn_ord[MAX_NUM_STEREO_BONDS]; /* ord. num. of the neighbor adjacent to the SB;        */
                                         /* starts from 0;                                       */
                                         /* -1 means removed explicit H                          */
    /* neighbors on both sides of stereobond have same parity => trans/T/E/2, diff. parities => cis/C/Z/1 */
    S_CHAR sb_parity[MAX_NUM_STEREO_BONDS];
    AT_NUMB sn_orig_at_num[MAX_NUM_STEREO_BONDS]; /* orig. at number of sn_ord[] neighbors        */

#if (FIND_RING_SYSTEMS == 1)
    S_CHAR bCutVertex;
    AT_NUMB nRingSystem;
    AT_NUMB nNumAtInRingSystem;
    AT_NUMB nBlockSystem;

#if (FIND_RINS_SYSTEMS_DISTANCES == 1)
    AT_NUMB nDistanceFromTerminal; /* terminal atom or ring system has 1, next has 2, etc. */
#endif

#endif
} inp_ATOM;

/*  v. 1.05 extensions: extended input supporting V3000; polymers   */
/*  OAD stands for OrigAtData                                       */

/* v. 1.05 Polymer stuff                                            */

/* Polymer representation type */
#define NO_POLYMER -1
#define POLYMER_REPRESENTATION_SOURCE_BASED 1
#define POLYMER_REPRESENTATION_STRUCTURE_BASED 2
#define POLYMER_REPRESENTATION_MIXED 3
#define POLYMER_REPRESENTATION_UNRECOGNIZED 4

#define ALLOW_MIXED_SRU_AND_MON 1 /**< allow simultaneous presence of source-based and structure-based units in embedding copolymer unit   */

#define POLYMER_STY_NON 0
#define POLYMER_STY_SRU 1 /**<  structure-based unit or copolymer subunit                        */
#define POLYMER_STY_MON 2 /**< source-based    polymer unit or copolymer subunit                    */
#define POLYMER_STY_COP 3 /**< copolymer unit embedding >1 subunits (may be SRU, MON, MER, CRO, MOD; others not supported yet)                               */
#define POLYMER_STY_MOD 4 /**< copolymer subunit only, designates chemical modification of SRU    */
#define POLYMER_STY_CRO 5 /**< copolymer subunit only, designates cross-linked version of SRU    */
#define POLYMER_STY_MER 6 /**< copolymer subunit only, source-based with no homopolymerize        */

#define POLYMER_SST_NON 0
#define POLYMER_SST_ALT 1
#define POLYMER_SST_RAN 2
#define POLYMER_SST_BLK 3

#define POLYMER_CONN_NON 0
#define POLYMER_CONN_HT 1
#define POLYMER_CONN_HH 2
#define POLYMER_CONN_EU 3

/**
 * @brief Structure describing a polymer unit
 *
 * @param id             It is what is called 'Sgroup number' in CTFILE
 * @param type           Type as by MDL format (STY)
 * @param subtype        Subtype as by MDL format (SST)
 * @param conn           Connection scheme  as by MDL format (SCN)
 * @param label          It is what is called 'unique Sgroup identifier' CTFILE
 * @param na             Number of atoms in the unit
 * @param nb             Number of bonds in the unit
 * @param cyclizable     =1 if frame shift via CRU 'cyclization' is applicable
 * @param cyclized       =1 if CRU already was frame_shift' treated
 * @param xbr1           Bracket ends coordinates (SDI)
 * @param xbr2           Bracket ends coordinates (SDI)
 * @param smt            Sgroup Subscript (SMT)
 * @param representation
 * @param cap1           Cap atom 1
 * @param end_atom1      End atom 1
 * @param end_atom2      End atom 2
 * @param cap2           Cap atom 2
 * @param cap1_is_undef  =1 if cap1 is undefined (star)
 * @param cap2_is_undef  =1 if cap2 is undefined (star)
 * @param alist          List of atoms in the unit (SAL)
 * @param blist          Bonds in the unit as list [atom1, atom2; atom1, atom2,..] for crossing bonds (S)
 * @param maxbkbonds     Max (allocd) number of frame_shift involved bonds
 * @param nbkbonds       Number of bkbonds in CRU main chain bkbonds are [breakable at possible frame shift] backbone bonds in between "left" end_atom1 and "right" end_atom2 ends (they may be only
 * single-order non-intraring bonds)
 * @param bkbonds        List of [breakable] backbone bonds [(a1,a2), (a3,a4),...]
 */
typedef struct OAD_PolymerUnit
{
    int id;         /* it is what is called 'Sgroup number' in CTFILE           */
    int type;       /* type as by MDL format (STY)                              */
    int subtype;    /* subtype as by MDL format (SST)                           */
    int conn;       /* connection scheme  as by MDL format (SCN)                */
    int label;      /* it is what is called 'unique Sgroup identifier' CTFILE   */
    int na;         /* number of atoms in the unit                              */
    int nb;         /* number of bonds in the unit                              */
    int cyclizable; /* =1 if frame shift via CRU 'cyclization' is applicable    */
    int cyclized;   /* =1 if CRU already was frame_shift' treated                 */
    double xbr1[4]; /* bracket ends coordinates (SDI)                           */
    double xbr2[4]; /* bracket ends coordinates (SDI)                           */
    char smt[80];   /* Sgroup Subscript (SMT)                                   */
    int representation;
    /* CRU structure is:            cap1-[-end_atom1/\/\/\end_atom2-]-cap2              */
    int cap1;
    int end_atom1;
    int end_atom2;
    int cap2;
    int cap1_is_undef;
    int cap2_is_undef;
    int *alist;     /* list of atoms in the unit (SAL)                          */
    int *blist;     /* bonds in the unit as list [atom1, atom2; atom1, atom2,..]*/
                    /* for crossing bonds (S)                                   */
    int maxbkbonds; /* max (allocd) number of frame_shift involved bonds		*/
    int nbkbonds;   /* number of bkbonds in CRU main chain						*/
                    /* bkbonds are [breakable at possible frame shift] backbone */
                    /* bonds in between "left" end_atom1 and "right" end_atom2  */
                    /* ends (they may be only single-order non-intraring bonds) */
    int **bkbonds;  /* list of [breakable] backbone bonds [(a1,a2), (a3,a4),...}*/
} OAD_PolymerUnit;

/**
 * @brief Structure describing a polymer
 *
 * @param units          array of pointers to units
 * @param n              number of units
 * @param n_pzz          number of polymeric Zz atoms
 * @param pzz            polymeric Zz atoms
 * @param really_do_frame_shift  flag indicating whether to do frame shift
 * @param frame_shift_scheme    frame shift analysis scheme
 * @param treat          treatment flag
 * @param representation representation type
 * @param is_in_reconn   flag indicating whether polymer is in reconnection mode
 * @param edit_repeats   flag indicating whether repeats need to be edited (-1 unknown, to be checked; 0 no, no edits required; 1 yes, edits are necessary)
 */
typedef struct OAD_Polymer
{
    OAD_PolymerUnit **units; /* array of pointers to units               */
    int n;
    int n_pzz; /* number of polymeric Zz atoms             */
               /* paired undefined-nature caps/stars around*/
               /* CRY's (must be even)                     */
    int *pzz;  /* polymeric Zz atoms                       */
    int really_do_frame_shift;
    int frame_shift_scheme; /* frame shift analysis scheme              */
    int treat;
    int representation;
    int is_in_reconn;
    int edit_repeats; /*  -1 unknown, to be checked               */
                      /*   0 no, no edits required                */
                      /*   1 yes, edits are necessary             */
} OAD_Polymer;

/**
 * @brief Structure describing atom properties for OAD
 *
 * @param erank         Rank of element; 2 - C, >2 - rank of heteroatom in chain, O > S > Se > Te > N ...., Rule 4
 * @param ring_erank    0 - not ring or just carbocycle, >2 - rank of senior heteroatom in this cycle, according to Rule 2 ( N > O >... )
 * @param ring_num      Ring number
 * @param ring_size     0 or ring system size
 *                      that is:
 *                      ring_erank != 0    heterocycle of ring_size
 *                      ring_erank==0 && ring_size>0    carbocycle of ring_size
 */
typedef struct OAD_AtProps
{
    int erank;      /* rank of element; 2 - C, >2 - rank of heteroatom in chain,    */
                    /* O > S > Se > Te > N ...., Rule 4                             */
    int ring_erank; /* 0 - not ring or just carbocycle,
                       >2 - rank of senior heteroatom in this cycle                 */
                    /* according to Rule 2 ( N > O >... )                           */
    int ring_num;
    int ring_size; /* 0 or ring system size                                        */
                   /* that is:                                                     */
                   /* ring_erank != 0    heterocycle of ring_size                  */
                   /* ring_erank==0 && ring_size>0    carbocycle of ring_size      */
} OAD_AtProps;

/**
 * @brief Structure describing v. 1.05 extended input supporting V3000; polymers
 *
 * @param n_non_star_atoms   Number of non-star atoms
 * @param n_star_atoms       Number of star atoms
 * @param atom_index_orig    index as supplied for atoms
 * @param atom_index_fin     = index or -1 for star atom
 * @param n_sgroups         currently, we do not use this.
 * @param n_3d_constraints   currently, we do not use this.
 * @param n_collections     number of collections
 * @param n_non_haptic_bonds number of non-haptic bonds
 * @param n_haptic_bonds     number of haptic bonds
 * @param lists_haptic_bonds  haptic_bonds[i] is pointer to int* which contains:
 *                           bond_type, non-star atom number, nendpts,
 *                           then endpts themselves
 * @param n_steabs          number of enhanced stereo absolute configurations
 * @param lists_steabs      steabs[k][0] - not used
 *                           steabs[k][1] -  number of members in collection
 *                           steabs[k][2..] - member atom numbers
 * @param n_sterel          number of enhanced stereo relative configurations
 * @param lists_sterel      sterel[k][0] - n from "STERELn" tag
 *                           sterel[k][1] -  number of members in collection
 *                           sterel[k][2..] - member atom numbers
 * @param n_sterac          number of enhanced stereo racemic configurations
 * @param lists_sterac      sterac[k][0] - n from "STERACn" tag
 *                           sterac[k][1] -  number of members in collection
 */
/* Extended input supports v. 1.05 extensions: V3000; polymers  */
typedef struct OAD_V3000
{
    int n_non_star_atoms;
    int n_star_atoms;
    int *atom_index_orig; /* index as supplied for atoms                                  */
    int *atom_index_fin;  /* = index or -1 for star atom                                  */
    int n_sgroups;        /* currently, we do not use this.                               */
    int n_3d_constraints; /* currently, we do not use this.                               */
    int n_collections;
    int n_non_haptic_bonds;
    int n_haptic_bonds;
    int **lists_haptic_bonds;
    /* haptic_bonds[i] is pointer to int* which contains:           */
    /* bond_type, non-star atom number, nendpts,                    */
    /* then endpts themselves                                       */
    /* Enhanced stereo */
    int n_steabs;
    int **lists_steabs; /* steabs[k][0] - not used                                      */
                        /* steabs[k][1] -  number of members in collection              */
                        /* steabs[k][2..] - member atom numbers                         */
    int n_sterel;
    int **lists_sterel; /* sterel[k][0] - n from "STERELn" tag                          */
                        /* sterel[k][1] -  number of members in collection              */
                        /* sterel[k][2..] - member atom numbers                         */
    int n_sterac;
    int **lists_sterac; /* sterac[k][0] - n from "STERACn" tag                          */
                        /* sterac[k][1] -  number of members in collection              */
                        /* sterac[k][2..] - member atom numbers                      */
} OAD_V3000;

/**
 * @brief Structure describing original atom data
 *
 * @param at               Array of input atoms
 * @param num_dimensions   Number of dimensions (2D or 3D)
 * @param num_inp_bonds    Number of input bonds
 * @param num_inp_atoms    Number of input atoms
 * @param num_components   Number of components in the structure
 * @param bDisconnectSalts Flag indicating whether salt disconnection is possible
 * @param bDisconnectCoord Flag indicating whether coordinate disconnection is needed
 * @param nCurAtLen        Array with current atom counts per component
 * @param nOldCompNumber   Array with old component numbers
 * @param nNumEquSets      Number of found component equivalence sets
 * @param nEquLabels       Array with equivalence labels for atoms
 * @param nSortedOrder     Array with sorted order of components
 * @param bSavedInINCHI_LIB Array indicating if saved in InChI library
 * @param bPreprocessed    Array indicating if preprocessed
 * @param szCoord          Array of molecular coordinates
 * @param polymer         Pointer to polymer data structure
 * @param v3000           Pointer to V3000 data structure
 * @param valid_polymer    Flag indicating if polymer data is valid
 * @param n_zy            Number of non-polymeric pseudoatoms (Zy)
 */
typedef struct tagOrigAtom
{
    /* Initially filled out by CreateOrigInpDataFromMolfile()                               */
    /* may be changed by disconnecting salts and disconnecting metals                       */
    inp_ATOM *at;
    int num_dimensions;
    int num_inp_bonds;
    int num_inp_atoms;
    /* may be changed by disconnecting salts and disconnecting metals                       */
    int num_components;   /* set by MarkDisconnectedComponents and disconn. metals*/
    int bDisconnectSalts; /* whether salt disconnection is possible               */
    int bDisconnectCoord; /* 0 if no disconnection needed                         */
                          /* else (NumImplicitH to disconnect)+1                  */
#if (bRELEASE_VERSION == 0)
    int bExtract;
#endif
    AT_NUMB *nCurAtLen;      /* has max_num_components elements                      */
    AT_NUMB *nOldCompNumber; /* 0 or component number in previous numbering          */
    int nNumEquSets;         /* number of found component equivalence sets           */
    AT_NUMB *nEquLabels;     /* num_inp_atoms elements, value>0 marks atoms          */
                             /* in the set #value                                    */
    AT_NUMB *nSortedOrder;   /* num_components elements, values = 1..num_components; */
                             /* only if num_components > 1                           */
    int bSavedInINCHI_LIB[INCHI_NUM];
    int bPreprocessed[INCHI_NUM];
    MOL_COORD *szCoord;
    /* v. 1.05 extensions                                                                   */
    OAD_Polymer *polymer;
    OAD_V3000 *v3000;
    int valid_polymer;
    int n_zy; /* number of non-polymeric pseudoatoms (Zy)             */

} ORIG_ATOM_DATA;

/**
 * @brief Structure describing the original structure
 *
 * @param num_atoms   Number of atoms
 * @param szAtoms     String representation of atoms
 * @param szBonds     String representation of bonds
 * @param szCoord     String representation of coordinates
 * @param polymer     Pointer to polymer data structure (uses pointer copy from orig_inp_data, do not free after use!)
 * @param v3000       Pointer to V3000 data structure (uses pointer copy from orig_inp_data, do not free after use!)
 * @param n_zy        Number of non-polymeric pseudoatoms (Zy)
 */
typedef struct tagOriginalStruct
{
    int num_atoms;
    char *szAtoms;
    char *szBonds;
    char *szCoord;
    /* v. 1.05 extensions                                                                   */
    OAD_Polymer *polymer; /* uses pointer copy from orig_inp_data, do not free after use! */
    OAD_V3000 *v3000;     /* uses pointer copy from orig_inp_data, do not free after use! */
    int n_zy;             /* number of non-polymeric pseudoatoms (Zy) */

} ORIG_STRUCT;

/**
 * @brief Structure describing atom parameters for drawing
 *
 * @param at_string                Atom string for drawing
 * @param DrawingLabelLeftShift    Left shift for drawing label
 * @param DrawingLabelLength       Length of drawing label
 * @param nCanonNbr                Canonical number; if zero then do not use all data for the atom
 * @param nCanonEquNbr             Canonical equivalence number
 * @param nTautGroupCanonNbr       Tautomeric group canonical number
 * @param nTautGroupEquNbr         Tautomeric group equivalence number
 * @param cFlags                   Additional flags (e.g., AT_FLAG_ISO_H_POINT)
 * @param nDebugData               Debug data (only if DISPLAY_DEBUG_DATA is defined)
 * @param cHighlightTheAtom        Flag indicating whether to highlight the atom
 * @param cStereoCenterParity      Stereo center parity
 * @param cStereoBondParity        Stereo bond parities
 * @param cStereoBondWarning       Stereo bond warnings
 * @param cStereoBondNumber        Stereo bond numbers
 */
typedef struct tagAtomParmsForDrawing
{
    char at_string[ATOM_INFO_LEN];
    int DrawingLabelLeftShift;
    int DrawingLabelLength;
    AT_NUMB nCanonNbr; /* if zero then do not use all data for the atom                */
    AT_NUMB nCanonEquNbr;
    AT_NUMB nTautGroupCanonNbr;
    AT_NUMB nTautGroupEquNbr;
    S_CHAR cFlags; /* AT_FLAG_ISO_H_POINT                                          */
#ifdef DISPLAY_DEBUG_DATA
    int nDebugData;
#endif
    S_CHAR cHighlightTheAtom;
    S_CHAR cStereoCenterParity;
    S_CHAR cStereoBondParity[MAX_STEREO_BONDS];
    S_CHAR cStereoBondWarning[MAX_STEREO_BONDS];
    S_CHAR cStereoBondNumber[MAX_STEREO_BONDS];
} inf_ATOM;

#define INF_STEREO_ABS 0x0001
#define INF_STEREO_REL 0x0002
#define INF_STEREO_RAC 0x0004
#define INF_STEREO_NORM 0x0008
#define INF_STEREO_INV 0x0010
#define INF_STEREO 0x0020
#define INF_STEREO_ABS_REL_RAC (INF_STEREO_ABS | INF_STEREO_REL | INF_STEREO_RAC)
#define INF_STEREO_NORM_INV (INF_STEREO_NORM | INF_STEREO_INV)

#define MAX_LEN_REMOVED_PROTONS 128

/**
 * @brief Structure describing information about atoms for InChI generation
 *
 * @param at              Array of atoms for InChI generation
 * @param num_at          Number of atoms
 * @param StereoFlags     Stereo flags
 * @param num_components  Number of components
 * @param pStereoFlags    Pointer to stereo flags per component
 * @param nNumRemovedProtons Number of removed protons
 * @param num_removed_iso_H  Number of removed/exchangeable isotopic hydrogens
 * @param num_iso_H       Array with number of removed isotopic hydrogens per isotope
 * @param szRemovedProtons String with removed protons
 */
typedef struct tagInfoAtomData
{
    inf_ATOM *at;
    int num_at;
    AT_NUMB StereoFlags;
    AT_NUMB num_components;
    AT_NUMB *pStereoFlags;

    int nNumRemovedProtons;
    int num_removed_iso_H;           /* number of exchangable isotopic H     */
    NUM_H num_iso_H[NUM_H_ISOTOPES]; /* number of exchangable isotopic H     */
    char szRemovedProtons[MAX_LEN_REMOVED_PROTONS];
} INF_ATOM_DATA;

/**
 * @brief Structure describing input atom data for InChI generation
 *
 * @param at                   Array of input atoms
 * @param at_fixed_bonds       Array of input atoms with fixed bonds (tautomeric case, added or removed H)
 * @param num_at               Number of atoms
 * @param num_removed_H        Number of removed hydrogens
 * @param num_bonds            Number of bonds
 * @param num_isotopic         Number of isotopic atoms
 * @param bExists              Flag indicating if the atom data exists
 * @param bDeleted             Flag indicating if the atom data is deleted
 * @param bHasIsotopicLayer    Flag indicating if the atom data has isotopic layer
 * @param bTautomeric          Flag indicating if the atom data is tautomeric
 * @param bTautPreprocessed    Flag indicating if the atom data is tautomeric and preprocessed
 * @param nNumRemovedProtons   Number of removed protons
 * @param nNumRemovedProtonsIsotopic Array with number of removed isotopic protons per isotope
 * @param num_iso_H            Array with number of isotopic hydrogens per isotope
 * @param bTautFlags           Tautomeric flags
 * @param bTautFlagsDone       Tautomeric flags done
 * @param bNormalizationFlags  Normalization flags
 */
typedef struct tagInputAtomData
{
    inp_ATOM *at;
    inp_ATOM *at_fixed_bonds; /* tautomeric case, added or removed H  */
    int num_at;
    int num_removed_H;
    int num_bonds;
    int num_isotopic;
    int bExists;
    int bDeleted;
    int bHasIsotopicLayer;
    int bTautomeric;
    int bTautPreprocessed;
    int nNumRemovedProtons;
    NUM_H nNumRemovedProtonsIsotopic[NUM_H_ISOTOPES];
    /* isotopic composition of removed protons, not included in num_iso_H[] */
    NUM_H num_iso_H[NUM_H_ISOTOPES];
    /* isotopic H on tautomerIC atoms and those in nIsotopicEndpointAtomNumber */
    INCHI_MODE bTautFlags;
    INCHI_MODE bTautFlagsDone;
    INCHI_MODE bNormalizationFlags;
} INP_ATOM_DATA;

typedef INP_ATOM_DATA INP_ATOM_DATA2[TAUT_NUM];

/**
 * @brief Structure describing normalization and canonicalization flags
 *
 * @param bTautFlags            Tautomeric flags
 * @param bTautFlagsDone        Tautomeric flags done
 * @param bNormalizationFlags   Normalization flags
 * @param nCanonFlags           Canonicalization flags
 */
typedef struct tagNormCanonFlags
{
    INCHI_MODE bTautFlags[INCHI_NUM][TAUT_NUM];
    INCHI_MODE bTautFlagsDone[INCHI_NUM][TAUT_NUM];
    INCHI_MODE bNormalizationFlags[INCHI_NUM][TAUT_NUM];
    int nCanonFlags[INCHI_NUM][TAUT_NUM];
} NORM_CANON_FLAGS;

/**
 * @brief Structure describing composite atom data for InChI generation
 *
 * @param at                   Array of input atoms
 * @param num_at               Number of atoms
 * @param num_removed_H        Number of removed hydrogens
 * @param num_bonds            Number of bonds
 * @param num_isotopic         Number of isotopic atoms
 * @param bExists              Flag indicating if the atom data exists
 * @param bDeleted             Flag indicating if the atom data is deleted
 * @param bHasIsotopicLayer    Flag indicating if the atom data has isotopic layer
 * @param bTautomeric          Flag indicating if the atom data is tautomeric
 * @param nNumRemovedProtons   Number of removed protons
 * @param nNumRemovedProtonsIsotopic Array with number of removed isotopic protons per isotope
 * @param num_iso_H            Array with number of isotopic hydrogens per isotope
 * @param nOffsetAtAndH        Array with offsets for atoms and hydrogens
 * @param num_components       Number of components
 */
typedef struct tagCompositeAtomData
{
    inp_ATOM *at;
    int num_at;
    int num_removed_H;
    int num_bonds;
    int num_isotopic;
    int bExists;
    int bDeleted; /* unused */
    int bHasIsotopicLayer;
    int bTautomeric;
    int nNumRemovedProtons;
    NUM_H nNumRemovedProtonsIsotopic[NUM_H_ISOTOPES];
    /* isotopic composition of removed protons, not included in num_iso_H[] */
    NUM_H num_iso_H[NUM_H_ISOTOPES];
    /* isotopic H on tautomeric atoms and those                             */
    /* in nIsotopicEndpointAtomNumber                                       */
    AT_NUMB *nOffsetAtAndH;
    int num_components;
} COMP_ATOM_DATA;

/*
typedef COMP_ATOM_DATA COMP_ATOM_DATA3[TAUT_NUM+1];
*/

#define ADD_LEN_STRUCT_FPTRS 100 /* allocation increments                                  */

typedef long INCHI_FPTR;

/**
 * @brief Structure describing file pointers for structures
 *
 * @param fptr          Array of file pointers to structures (input: fptr[cur_fptr] =file ptr to the struct to read; output: fptr[cur_fptr+1]=file ptr to the next struct or EOF)
 * @param len_fptr      Allocated length of fptr
 * @param cur_fptr      Current file pointer index (input: k-1 to read the kth struct, k = 1, 2, 3,...; left unchanged; struct number := cur_fptr+1)
 * @param max_fptr      Length of the filled out portion of fptr
 */
typedef struct tagStructFptrs
{
    INCHI_FPTR *fptr; /* input:  fptr[cur_fptr]  =file ptr to the struct to read      */
                      /* output: fptr[cur_fptr+1]=file ptr to the next struct or EOF  */
    int len_fptr;     /* allocated length of fptr                                     */
    int cur_fptr;     /* input: k-1 to read the kth struct, k = 1, 2, 3,...;          */
                      /* left unchanged; struct number := cur_fptr+1                  */
    int max_fptr;     /* length of the filled out portion of fptr                     */
} STRUCT_FPTRS;

/* forward declaration */
struct tagCANON_GLOBALS;

#ifndef __ICHITIME_H__
struct tagInchiTime;
struct tagINCHI_CLOCK;
int bInchiTimeIsOver(struct tagINCHI_CLOCK *ic, struct tagInchiTime *TickEnd);
#endif

#define FLAG_INP_AT_CHIRAL 1
#define FLAG_INP_AT_NONCHIRAL 2
#define FLAG_SET_INP_AT_CHIRAL 4
#define FLAG_SET_INP_AT_NONCHIRAL 8
#define FLAG_SET_INP_LARGE_MOLS 16

#ifndef COMPILE_ALL_CPP
#ifdef __cplusplus
extern "C"
{
#endif
#endif

    /**
     * @brief Create original input data from molfile
     *
     * @param inp_file Input file stream
     * @param orig_at_data Pointer to original atom data structure to be filled
     * @param bMergeAllInputStructures Flag indicating whether to merge all input structures
     * @param bGetOrigCoord Flag indicating whether to get original coordinates
     * @param bDoNotAddH Flag indicating whether to avoid adding hydrogens
     * @param treat_polymers Flag indicating how to treat polymers
     * @param treat_NPZz Flag indicating how to treat non-polymeric Zz atoms
     * @param pSdfLabel Pointer to SDF label
     * @param pSdfValue Pointer to SDF value
     * @param lSdfId Pointer to SDF ID
     * @param lMolfileNumber Pointer to molfile number
     * @param pInpAtomFlags Pointer to input atom flags
     * @param err Pointer to error code
     * @param pStrErr Pointer to error string
     * @param bNoWarnings Flag indicating whether to suppress warnings
     *
     * @return int Status code
     */
    int CreateOrigInpDataFromMolfile(INCHI_IOSTREAM *inp_file, ORIG_ATOM_DATA *orig_at_data, int bMergeAllInputStructures, int bGetOrigCoord, int bDoNotAddH, int treat_polymers, int treat_NPZz,
                                     const char *pSdfLabel, char *pSdfValue, unsigned long *lSdfId, long *lMolfileNumber, INCHI_MODE *pInpAtomFlags, int *err, char *pStrErr, int bNoWarnings);

    /**
     * @brief Convert InChI representation to original atom data
     *
     * @param infile Input file stream
     * @param orig_at_data Pointer to original atom data structure to be filled
     * @param bMergeAllInputStructures Flag indicating whether to merge all input structures
     * @param bGetOrigCoord Flag indicating whether to get original coordinates
     * @param bDoNotAddH Flag indicating whether to avoid adding hydrogens
     * @param vABParityUnknown Flag indicating whether to treat vAB parity as unknown
     * @param nInputType Input type
     * @param pSdfLabel Pointer to SDF label
     * @param pSdfValue Pointer to SDF value
     * @param lSdfId Pointer to SDF ID
     * @param pInpAtomFlags Pointer to input atom flags
     * @param err Pointer to error code
     * @param pStrErr Pointer to error string
     * @return int Status code
     */
    int InchiToOrigAtom(INCHI_IOSTREAM *infile, ORIG_ATOM_DATA *orig_at_data, int bMergeAllInputStructures, int bGetOrigCoord, int bDoNotAddH, int vABParityUnknown, INPUT_TYPE nInputType, char *pSdfLabel,
                        char *pSdfValue, unsigned long *lSdfId, INCHI_MODE *pInpAtomFlags, int *err, char *pStrErr);

    /**
     * @brief Mark disconnected components in the original atom data
     *
     * @param orig_at_data Pointer to original atom data structure
     * @param bProcessOldCompNumbers Flag indicating whether to process old component numbers
     * @return int Status code
     */
    int MarkDisconnectedComponents(ORIG_ATOM_DATA *orig_at_data, int bProcessOldCompNumbers);

    /**
     * @brief Disconnect salts in the original atom data
     *
     * @param orig_inp_data Pointer to original atom data structure
     * @param bDisconnect Flag indicating whether to disconnect salts
     * @return int Status code
     */
    int DisconnectSalts(ORIG_ATOM_DATA *orig_inp_data, int bDisconnect);

    /**
     * @brief Disconnect metals in the original atom data
     *
     * @param orig_inp_data Pointer to original atom data structure
     * @param bCheckMetalValence Flag indicating whether to check metal valence
     * @param bTautFlagsDone Pointer to tautomeric flags done
     * @return int Status code
     */
    int DisconnectMetals(ORIG_ATOM_DATA *orig_inp_data, int bCheckMetalValence, INCHI_MODE *bTautFlagsDone);

    /**
     * @brief Check if metals may be disconnected in the original atom data
     *
     * @param orig_inp_data Pointer to original atom data structure
     * @param bCheckMetalValence Flag indicating whether to check metal valence
     * @param bTautFlagsDone Pointer to tautomeric flags done
     * @return int Status code
     */
    int bMayDisconnectMetals(ORIG_ATOM_DATA *orig_inp_data, int bCheckMetalValence, INCHI_MODE *bTautFlagsDone);

    /**
     * @brief Check if there is a metal atom in the original atom data
     *
     * @param orig_inp_data Pointer to original atom data structure
     * @return int Status code
     */
    int bHasMetalAtom(ORIG_ATOM_DATA *orig_inp_data);

    /**
     * @brief Fix adjacent radicals in the input atoms
     *
     * @param num_inp_atoms Number of input atoms
     * @param at Array of input atoms
     * @return int Status code
     */
    int FixAdjacentRadicals(int num_inp_atoms, inp_ATOM *at); /* FIX_ADJ_RAD == 1 */

    /**
     * @brief Fix odd things in the input atoms
     *
     * @param num_atoms Number of input atoms
     * @param at Array of input atoms
     * @param bFixBug Flag indicating whether to fix bugs
     * @param bFixNonUniformDraw Flag indicating whether to fix non-uniform drawing
     * @return int Status code
     */
    int fix_odd_things(int num_atoms, inp_ATOM *at, int bFixBug, int bFixNonUniformDraw);

    /**
     * @brief Post-fix odd things in the input atoms (does nothing, returns 0)
     *
     * @param num_atoms Number of input atoms
     * @param at Array of input atoms
     * @return int Status code
     */
    int post_fix_odd_things(int num_atoms, inp_ATOM *at);

    /**
     * @brief Remove ion pairs from the input atoms
     *
     * @param num_atoms Number of input atoms
     * @param at Array of input atoms
     * @return int Status code
     */
    int remove_ion_pairs(int num_atoms, inp_ATOM *at);

    /**
     * @brief Check if a feature is found in the input atoms ???
     *
     * @param at Array of input atoms
     * @param num_atoms Number of input atoms
     * @return int Status code
     */
    int bFoundFeature(inp_ATOM *at, int num_atoms);

    /**
     * @brief Unmark ring systems in the input atoms
     *
     * @param at Array of input atoms
     * @param num_atoms Number of input atoms
     * @return int Status code
     */
    int UnMarkRingSystemsInp(inp_ATOM *at, int num_atoms);

    /**
     * @brief Write original atom data to SD file
     *
     * @param inp_at_data Pointer to original atom data structure
     * @param fcb File stream to write to
     * @param name Name of the structure
     * @param comment Comment to include
     * @param bChiralFlag Flag indicating for chiral information
     * @param bAtomsDT Flag indicating how H isotopes (D/T) are represented in the SD output
     * @param szLabel SDF label
     * @param szValue SDF value
     * @return int Status code
     */
    int OrigAtData_WriteToSDfile(const ORIG_ATOM_DATA *inp_at_data, INCHI_IOSTREAM *fcb, const char *name, const char *comment, int bChiralFlag, int bAtomsDT, const char *szLabel, const char *szValue);

    /**
     * @brief Free input atom array
     *
     * @param at Pointer to input atom array
     */
    void FreeInpAtom(inp_ATOM **at);

    /**
     * @brief Free information atom array (same as FreeInpAtom)
     *
     * @param at Pointer to information atom array
     */
    void FreeInfAtom(inf_ATOM **at);

    /**
     * @brief Free original atom data
     *
     * @param orig_at_data Pointer to original atom data structure
     */
    void FreeOrigAtData(ORIG_ATOM_DATA *orig_at_data);

    /**
     * @brief Free extended original atom data (polymers and V3000)
     *
     * @param pd Pointer to polymer data structure
     * @param v3k Pointer to V3000 data structure
     */
    void FreeExtOrigAtData(OAD_Polymer *pd, OAD_V3000 *v3k);

    /**
     * @brief Free input atom data
     *
     * @param inp_at_data Pointer to input atom data structure
     */
    void FreeInpAtomData(INP_ATOM_DATA *inp_at_data);

    /**
     * @brief Free composite atom data
     *
     * @param inp_at_data Pointer to composite atom data structure
     */
    void FreeCompAtomData(COMP_ATOM_DATA *inp_at_data);

    /**
     * @brief Free information atom data
     *
     * @param inf_at_data Pointer to information atom data structure
     */
    void FreeInfoAtomData(INF_ATOM_DATA *inf_at_data);

    /*
        ORIG_ATOM_DATA
            functions
    */

    /* OAD_Edit */

    /**
     * @brief Structure describing structure edits
     *
     * @param del_atom    Array of deleted atom indices
     * @param del_bond    Array of deleted bond indices
     * @param new_bond    Array of new bond indices
     * @param mod_bond    Array of modified bond indices
     * @param mod_coord   Array of modified coordinate indices
     * @param del_side_chains Flag indicating whether to delete side chains
     */
    typedef struct tagOAD_StructureEdits
    {
        INT_ARRAY *del_atom;
        INT_ARRAY *del_bond;
        INT_ARRAY *new_bond;
        INT_ARRAY *mod_bond;
        INT_ARRAY *mod_coord;
        int del_side_chains;
    } OAD_StructureEdits;

    /**
     * @brief Initialize structure edits
     *
     * @param ed Pointer to structure edits
     * @return int Error status code
     */
    int OAD_StructureEdits_Init(OAD_StructureEdits *ed);

    /**
     * @brief Clear structure edits
     *
     * @param ed Pointer to structure edits
     */
    void OAD_StructureEdits_Clear(OAD_StructureEdits *ed);

    /**
     * @brief Debug print structure edits
     *
     * @param ed Pointer to structure edits
     */
    void OAD_StructureEdits_DebugPrint(OAD_StructureEdits *ed);

#if (RING2CHAIN == 1)

    /**
     * @brief Convert rings to chains in the original atom data
     *
     * @param ic Pointer to INCHI_CLOCK structure
     * @param pCG Pointer to CANON_GLOBALS structure
     * @param orig_inp_data Pointer to original atom data structure
     * @return int Status code or the number of ring cuts made
     */
    int Ring2Chain(struct tagINCHI_CLOCK *ic, struct tagCANON_GLOBALS *pCG, ORIG_ATOM_DATA *orig_inp_data);
#endif

#if (UNDERIVATIZE == 1)
    /**
     * @brief Underivatize the original atom data (main underivatization procedure)
     *
     * @param ic Pointer to INCHI_CLOCK structure
     * @param pCG Pointer to CANON_GLOBALS structure
     * @param orig_inp_data Pointer to original atom data structure
     * @param bOutputSdf Flag indicating whether to output SDF
     * @param bOutputReport Flag indicating whether to output report
     * @param pSdfValue Pointer to SDF value
     * @return int Status code
     */
    int OAD_Edit_Underivatize(struct tagINCHI_CLOCK *ic, struct tagCANON_GLOBALS *pCG, ORIG_ATOM_DATA *orig_inp_data, int bOutputSdf, int bOutputReport, char *pSdfValue);

    /**
     * @brief Merge components and recreate original atom data
     *
     * @param orig_OrigAtomData Pointer to original atom data structure
     * @param curr_InpAtomData Pointer to current input atom data structure
     * @param num_components Number of components
     * @param errcode Pointer to error code
     */
    void OAD_Edit_MergeComponentsAndRecreateOAD(ORIG_ATOM_DATA *orig_OrigAtomData, INP_ATOM_DATA *curr_InpAtomData, int num_components, int *errcode);

#endif

    /* OAD misc. */

    /**
     * @brief Validate polymer and pseudo-element data in the original atom data
     *
     * @param orig_at_data Pointer to original atom data structure
     * @param treat_polymers Flag indicating how to treat polymers
     * @param bNPZz Flag indicating whether to treat non-polymeric Zz atoms
     * @param pStrErr Pointer to error string
     * @param bNoWarnings Flag indicating whether to suppress warnings
     * @return int Status code
     */
    int OAD_ValidatePolymerAndPseudoElementData(ORIG_ATOM_DATA *orig_at_data, int treat_polymers, int bNPZz, char *pStrErr, int bNoWarnings);

    /* OAD_Polymer */
    /**
     * @brief Validate and sort out pseudo-element atoms in the original atom data
     *
     * @param orig_at_data Pointer to original atom data structure
     * @param treat_polymers Flag indicating how to treat polymers (POLYMER_NO, POLYMER_MODERN, POLYMER_LEGACY, POLYMER_LEGACY_PLUS)
     * @param use_zz Flag indicating whether to use Zz atoms
     * @param err Pointer to error code
     * @param pStrErr Pointer to error string
     */
    void OAD_ValidateAndSortOutPseudoElementAtoms(ORIG_ATOM_DATA *orig_at_data, int treat_polymers, int use_zz, int *err, char *pStrErr);

    /**
     * @brief Get polymer representation type
     *
     * @param p Pointer to polymer data structure
     * @return int Polymer representation type (NO_POLYMER, POLYMER_REPRESENTATION_SOURCE_BASED, POLYMER_REPRESENTATION_STRUCTURE_BASED, POLYMER_REPRESENTATION_MIXED,POLYMER_REPRESENTATION_UNRECOGNIZED)
     */
    int OAD_Polymer_GetRepresentation(OAD_Polymer *p);

    /**
     * @brief Cyclize closeable units in the original atom data
     *
     * @param orig_at_data Pointer to original atom data structure
     * @param use_zz Flag indicating whether to use Zz atoms
     * @param pStrErr Pointer to error string
     * @param bNoWarnings Flag indicating whether to suppress warnings
     * @return int Status code
     */
    int OAD_Polymer_CyclizeCloseableUnits(ORIG_ATOM_DATA *orig_at_data, int use_zz, char *pStrErr, int bNoWarnings);

    /**
     * @brief Find ring systems in the polymer data
     *
     * @param pd Pointer to polymer data structure
     * @param at Array of input atoms
     * @param nat Number of input atoms
     * @param num_inp_bonds Pointer to number of input bonds
     * @param num_ring_sys Pointer to number of ring systems found
     * @param size_ring_sys Pointer to size of ring systems found
     * @param start Starting atom index
     * @return int Status code
     */
    int OAD_Polymer_FindRingSystems(OAD_Polymer *pd, inp_ATOM *at, int nat, int *num_inp_bonds, int *num_ring_sys, int *size_ring_sys, int start);

    /**
     * @brief Set atom properties for the polymer data
     *
     * @param pd Pointer to polymer data structure
     * @param at Array of input atoms
     * @param nat Number of input atoms
     * @param num_inp_bonds Pointer to number of input bonds
     * @param aprops Pointer to atom properties structure
     * @param cano_nums Pointer to canonical numbers array
     */
    void OAD_Polymer_SetAtProps(OAD_Polymer *pd, inp_ATOM *at, int nat, int *num_inp_bonds, OAD_AtProps *aprops, int *cano_nums);

    /**
     * @brief Compare backbone bonds seniority (for sorting SRU cyclizing bonds (PS=='frame-shift') in descending order)
     *
     * @param b1 Pointer to first backbone bond
     * @param b2 Pointer to second backbone bond
     * @param aprops Pointer to atom properties structure
     * @return int Comparison result (-1 if b1 < b2, 0 if equal, 1 if b1 > b2)
     */
    int OAD_Polymer_CompareBackboneBondsSeniority(int *b1, int *b2, OAD_AtProps *aprops);

    /**
     * @brief Compare ranks of two atoms (compare seniority of two atoms in polymer SRU)
     *
     * @param atom1 Index of the first atom
     * @param atom2 Index of the second atom
     * @param aprops Pointer to atom properties structure
     * @return int Comparison result
     */
    int OAD_Polymer_CompareRanksOfTwoAtoms(int atom1, int atom2, OAD_AtProps *aprops);

    /**
     * @brief Check if the first atom has a lower rank than the second atom (check seniority of two atoms in polymer SRU)
     *
     * @param atom1 Index of the first atom
     * @param atom2 Index of the second atom
     * @param aprops Pointer to atom properties structure
     * @return int 1 if atom1 has a lower rank than atom2, 0 otherwise
     */
    int OAD_Polymer_IsFirstAtomRankLower(int atom1, int atom2, OAD_AtProps *aprops);

    /**
     * @brief Set reopening details for a polymer unit
     *
     * @param u Pointer to polymer unit
     * @param at Array of input atoms
     * @return int Number of breakable bonds
     */
    int OAD_PolymerUnit_SetReopeningDetails(OAD_PolymerUnit *u, inp_ATOM *at);

    /**
     * @brief Sort backbone bonds and set seniors for a polymer unit
     *
     * @param u Pointer to polymer unit
     * @param at Array of input atoms
     * @param aprops Pointer to atom properties structure
     * @param senior_bond Pointer to array to store senior bond indices
     */
    void OAD_PolymerUnit_SortBackboneBondsAndSetSeniors(OAD_PolymerUnit *u, inp_ATOM *at, OAD_AtProps *aprops, int *senior_bond);

    /**
     * @brief Find backbones in the original atom data
     *
     * @param where_to_look Pointer to original atom data structure
     * @param composite_norm_data Pointer to composite atom data structure
     * @param err Pointer to error code
     * @param pStrErr Pointer to error string
     */
    void OAD_Polymer_FindBackbones(ORIG_ATOM_DATA *where_to_look, COMP_ATOM_DATA *composite_norm_data, int *err, char *pStrErr);

    /**
     * @brief Prepare working set for polymer processing (replace original atom numbers in polymer data with (canonical num + 1))
     *
     * @param p Pointer to polymer data structure
     * @param cano_nums Pointer to canonical numbers array
     * @param compnt_nums Pointer to component numbers array
     * @param units2 Pointer to array of polymer units (a copy of original polymer units with updated atom numbers)
     * @param unum Pointer to number of polymer units
     * @return int Status code
     */
    int OAD_Polymer_PrepareWorkingSet(OAD_Polymer *p, int *cano_nums, int *compnt_nums, OAD_PolymerUnit **units2, int *unum);

    /**
     * @brief Free polymer data structure
     *
     * @param p Pointer to polymer data structure
     */
    void OAD_Polymer_Free(OAD_Polymer *p);

    /**
     * @brief Debug trace for polymer data structure (print the whole polymer data)
     *
     * @param p Pointer to polymer data structure
     */
    void OAD_Polymer_DebugTrace(OAD_Polymer *p);

    /**
     * @brief Smartly reopen cyclized units in the original atom data (open pre-cyclized CRUs appropriately (i.e., make frame shift))
     *
     * @param p Pointer to polymer data structure
     * @param at Pointer to input atom data
     * @param nat Number of input atoms
     * @param num_inp_bonds Pointer to number of input bonds
     */
    void OAD_Polymer_SmartReopenCyclizedUnits(OAD_Polymer *p, inp_ATOM *at, int nat, int *num_inp_bonds);

    /**
     * @brief Prepare CRU edits for polymer data structure (prepare CRU fold edits as suggested by the strings with preliminary generated interim (1.05+ flavoured) InChI and AuxInfo)
     *
     * @param orig_at_data Pointer to original atom data structure
     * @param sinchi_noedits Pointer to InChI string without edits
     * @param saux_noedits Pointer to auxiliary InChI string without edits
     * @param sinchi Pointer to InChI string with edits
     * @param saux Pointer to auxiliary InChI string with edits
     * @param ed Pointer to structure edits to be filled
     * @return int Status code
     */
    int OAD_Polymer_PrepareFoldCRUEdits(ORIG_ATOM_DATA *orig_at_data, char *sinchi_noedits, char *saux_noedits, char *sinchi, char *saux, OAD_StructureEdits *ed);

    /**
     * @brief Prepare frame shift edits for polymer data structure (prepare frame shift edits as suggested by the strings with preliminary generated interim (1.05+ flavoured) InChI and AuxInfo)
     *
     * @param orig_at_data Pointer to original atom data structure
     * @param sinchi Pointer to InChI string with edits
     * @param saux Pointer to auxiliary InChI string with edits
     * @param ed Pointer to structure edits to be filled
     * @return int Status code
     */
    int OAD_Polymer_PrepareFrameShiftEdits(ORIG_ATOM_DATA *orig_at_data, char *sinchi, char *saux, OAD_StructureEdits *ed);

    /* OAD_PolymerUnit */

    /**
     * @brief Create a new polymer unit
     *
     * @param maxatoms Maximum number of atoms
     * @param maxbonds Maximum number of bonds
     * @param id Unique identifier for the polymer unit
     * @param label Label for the polymer unit
     * @param type Type of the polymer unit
     * @param subtype Subtype of the polymer unit
     * @param conn Connection type
     * @param smt SMARTS pattern?
     * @param na Number of atoms
     * @param alist Pointer to the atom list
     * @param nb Number of bonds
     * @param blist Pointer to the bond list
     * @param nbkbonds Number of backbone bonds
     * @param bkbonds Pointer to the backbone bond list
     * @return OAD_PolymerUnit* Pointer to the new polymer unit
     */
    OAD_PolymerUnit *OAD_PolymerUnit_New(int maxatoms, int maxbonds, int id, int label, int type, int subtype, int conn, char *smt, int na, INT_ARRAY *alist, int nb, INT_ARRAY *blist, int nbkbonds,
                                         int **bkbonds);

    /**
     * @brief Create a copy of a polymer unit
     *
     * @param u Pointer to the original polymer unit
     * @return OAD_PolymerUnit* Pointer to the new polymer unit
     */
    OAD_PolymerUnit *OAD_PolymerUnit_CreateCopy(OAD_PolymerUnit *u);

    /**
     * @brief Free a polymer unit
     *
     * @param unit Pointer to the polymer unit to free
     */
    void OAD_PolymerUnit_Free(OAD_PolymerUnit *unit);

    /**
     * @brief Debug trace for a polymer unit (print the whole polymer unit)
     *
     * @param unit Pointer to the polymer unit
     */
    void OAD_PolymerUnit_DebugTrace(OAD_PolymerUnit *unit);

    /**
     * @brief Find the ends and caps of a polymer unit
     *
     * @param unit Pointer to the polymer unit
     * @param orig_at_data Pointer to the original atom data
     * @param end1 Pointer to the first end atom
     * @param cap1 Pointer to the first cap atom
     * @param cap1_is_undef Pointer to the first cap atom's undefined status
     * @param end2 Pointer to the second end atom
     * @param cap2 Pointer to the second cap atom
     * @param cap2_is_undef Pointer to the second cap atom's undefined status
     * @param err Pointer to the error status
     * @param pStrErr Pointer to the error string
     */
    void OAD_PolymerUnit_FindEndsAndCaps(OAD_PolymerUnit *unit, ORIG_ATOM_DATA *orig_at_data, int *end1, int *cap1, int *cap1_is_undef, int *end2, int *cap2, int *cap2_is_undef, int *err, char *pStrErr);

    /**
     * @brief Set the ends and caps of a polymer unit
     *
     * @param unit Pointer to the polymer unit
     * @param orig_at_data Pointer to the original atom data
     * @param err Pointer to the error status
     * @param pStrErr Pointer to the error string
     */
    void OAD_PolymerUnit_SetEndsAndCaps(OAD_PolymerUnit *unit, ORIG_ATOM_DATA *orig_at_data, int *err, char *pStrErr);

    /**
     * @brief Unlink caps and connect end atoms of a polymer unit
     *
     * @param unit Pointer to the polymer unit
     * @param orig_at_data Pointer to the original atom data
     * @param err Pointer to the error status
     * @param pStrErr Pointer to the error string
     */
    void OAD_PolymerUnit_UnlinkCapsAndConnectEndAtoms(OAD_PolymerUnit *unit, ORIG_ATOM_DATA *orig_at_data, int *err, char *pStrErr);

    /**
     * @brief Prepare a polymer unit for frame shifting
     *
     * @param unit Pointer to the polymer unit
     * @param orig_at_data Pointer to the original atom data
     * @param err Pointer to the error status
     * @param pStrErr Pointer to the error string
     */
    void OAD_PolymerUnit_PrepareToFrameShift(OAD_PolymerUnit *unit, ORIG_ATOM_DATA *orig_at_data, int *err, char *pStrErr);

    /**
     * @brief Find reachable atoms from a given start atom
     *
     * @param orig_at_data Pointer to the original atom data
     * @param start_atom The starting atom
     * @param nforbidden_bonds The number of forbidden bonds
     * @param forbidden_bonds Array of forbidden bonds
     * @param n_reachable Pointer to the number of reachable atoms
     * @param reachable Pointer to the array of reachable atoms
     * @param err Pointer to the error status
     * @param pStrErr Pointer to the error string
     * @return int Status code
     */
    int OAD_CollectReachableAtoms(ORIG_ATOM_DATA *orig_at_data, int start_atom, int nforbidden_bonds, int *forbidden_bonds, int *n_reachable, int *reachable, int *err, char *pStrErr);

    /**
     * @brief Collect backbone atoms from a polymer unit
     *
     * @param at_data Pointer to the original atom data
     * @param na The number of atoms
     * @param alist Array of atom indices
     * @param end_atom1 The first end atom
     * @param end_atom2 The second end atom
     * @param nbkatoms Pointer to the number of backbone atoms
     * @param bkatoms Pointer to the array of backbone atoms
     * @param err Pointer to the error status
     * @param pStrErr Pointer to the error string
     */
    void OAD_CollectBackboneAtoms(ORIG_ATOM_DATA *at_data, int na, int *alist, int end_atom1, int end_atom2, int *nbkatoms, int *bkatoms, int *err, char *pStrErr);

    /**
     * @brief Collect backbone bonds from a polymer unit
     *
     * @param at_data Pointer to the original atom data
     * @param na The number of atoms
     * @param alist Array of atom indices
     * @param end_atom1 The first end atom
     * @param end_atom2 The second end atom
     * @param nbkbonds Pointer to the number of backbone bonds
     * @param bkbonds Pointer to the array of backbone bonds
     * @param err Pointer to the error status
     * @param pStrErr Pointer to the error string
     */
    void OAD_CollectBackboneBonds(ORIG_ATOM_DATA *at_data, int na, int *alist, int end_atom1, int end_atom2, int *nbkbonds, int **bkbonds, int *err, char *pStrErr);

    /**
     * @brief Collect backbone bonds from a polymer unit and delist intra-ring backbone bonds
     *
     * @param unit Pointer to the polymer unit
     * @param orig_at_data Pointer to the original atom data
     * @param err Pointer to the error status
     * @param pStrErr Pointer to the error string
     */
    void OAD_PolymerUnit_DelistIntraRingBackboneBonds(OAD_PolymerUnit *unit, ORIG_ATOM_DATA *orig_at_data, int *err, char *pStrErr);

    /**
     * @brief Collect high-order backbone bonds from a polymer unit
     *
     * @param unit Pointer to the polymer unit
     * @param orig_at_data Pointer to the original atom data
     * @param composite_norm_data Pointer to the composite atom data
     * @param err Pointer to the error status
     * @param pStrErr Pointer to the error string
     */
    void OAD_PolymerUnit_DelistHighOrderBackboneBonds(OAD_PolymerUnit *unit, ORIG_ATOM_DATA *orig_at_data, COMP_ATOM_DATA *composite_norm_data, int *err, char *pStrErr);

    /**
     * @brief Collect backbone bonds from a polymer unit and sort them
     *
     * @param u Pointer to the polymer unit
     * @param aprops Pointer to the atom properties
     * @param bnum Pointer to the number of backbone bonds
     */
    void OAD_PolymerUnit_SortBackboneBonds(OAD_PolymerUnit *u, OAD_AtProps *aprops, int *bnum);

    /**
     * @brief Collect backbone atoms from a polymer unit
     *
     * @param u Pointer to the polymer unit
     * @param at Pointer to the atom data
     * @param aprops Pointer to the atom properties
     * @param nat Pointer to the number of atoms
     * @param num_inp_bonds Pointer to the number of input bonds
     */
    void OAD_PolymerUnit_ReopenCyclized(OAD_PolymerUnit *u, inp_ATOM *at, OAD_AtProps *aprops, int nat, int *num_inp_bonds);

    /**
     * @brief Compare two polymer units based on their atom lists
     *
     * @param u1 Pointer to the first polymer unit
     * @param u2 Pointer to the second polymer unit
     * @return int comparison result (-1, 0, 1)
     */
    int OAD_PolymerUnit_CompareAtomLists(OAD_PolymerUnit *u1, OAD_PolymerUnit *u2);

    /**
     * @brief Compare two polymer units based on their atom lists with modifications
     *
     * @param u1 Pointer to the first polymer unit
     * @param u2 Pointer to the second polymer unit
     * @return int comparison result (-1, 0, 1)
     */
    int OAD_PolymerUnit_CompareAtomListsMod(OAD_PolymerUnit *u1, OAD_PolymerUnit *u2);

    /**
     * @brief Order bond atoms and bonds themselves in a polymer unit
     *
     * @param u Pointer to the polymer unit
     * @param n_stars Pointer to the number of stars
     * @param stars Pointer to the array of star atoms
     * @return int
     */
    int OAD_PolymerUnit_OrderBondAtomsAndBondsThemselves(OAD_PolymerUnit *u, int n_stars, int *stars);

    /**
     * @brief Check if a polymer unit has a specific atom (e.g., metal)
     *
     * @param u Pointer to the polymer unit
     * @param at Pointer to the atom data
     * @return int Status code (1 if an atom is found, 0 otherwise)
     */
    int OAD_PolymerUnit_HasMetal(OAD_PolymerUnit *u, inp_ATOM *at);

    /**
     * @brief Fix unknown stereo bonds in the input atoms
     *
     * @param at Pointer to the atom data
     * @param num_at Pointer to the number of atoms
     * @return int Number of fixed stereo bonds
     */
    int FixUnkn0DStereoBonds(inp_ATOM *at, int num_at);

    /**
     * @brief Create/allocate a Inf Atom structure
     *
     * @param num_atoms Number of atoms
     * @return inf_ATOM* Pointer to the created Inf Atom structure
     */
    inf_ATOM *CreateInfAtom(int num_atoms);

    /**
     * @brief Create/allocate a Inp Atom structure
     *
     * @param num_atoms Number of atoms
     * @return inp_ATOM* Pointer to the created Inp Atom structure
     */
    inp_ATOM *CreateInpAtom(int num_atoms);

    /**
     * @brief Create a Info Atom Data structure
     *
     * @param inf_at_data Pointer to the Info Atom Data structure
     * @param num_atoms Number of atoms
     * @param num_components Number of components
     * @return int Status code (1 if successful, 0 otherwise)
     */
    int CreateInfoAtomData(INF_ATOM_DATA *inf_at_data, int num_atoms, int num_components);

    /**
     * @brief Allocate memory for Info Atom Data structure
     *
     * @param inf_at_data Pointer to the Info Atom Data structure
     * @param num_atoms Number of atoms
     * @param num_components Number of components
     * @return int Status code (1 if successful, 0 otherwise)
     */
    int AllocateInfoAtomData(INF_ATOM_DATA *inf_at_data, int num_atoms, int num_components);

    /**
     * @brief Duplicate Info Atom Data structure
     *
     * @param inf_at_data_to Pointer to the destination Info Atom Data structure
     * @param inf_at_data_from Pointer to the source Info Atom Data structure
     * @return int Status code (1 if successful, 0 otherwise)
     */
    int DuplicateInfoAtomData(INF_ATOM_DATA *inf_at_data_to, const INF_ATOM_DATA *inf_at_data_from);

    /**
     * @brief Create a Inp Atom Data structure
     *
     * @param inp_at_data Pointer to the Inp Atom Data structure
     * @param num_atoms Number of atoms
     * @param create_at_fixed_bonds Flag to create fixed bonds
     * @return int Status code (1 if successful, 0 otherwise)
     */
    int CreateInpAtomData(INP_ATOM_DATA *inp_at_data, int num_atoms, int create_at_fixed_bonds);

    /**
     * @brief Create a Comp Atom Data structure
     *
     * @param inp_at_data Pointer to the Comp Atom Data structure
     * @param num_atoms Number of atoms
     * @param num_components Number of components
     * @param bIntermediateTaut Flag for intermediate tautomer
     * @return int Status code (1 if successful, 0 otherwise)
     */
    int CreateCompAtomData(COMP_ATOM_DATA *inp_at_data, int num_atoms, int num_components, int bIntermediateTaut);

#ifndef COMPILE_ANSI_ONLY
    int DisplayInputStructure(char *szOutputString, inp_ATOM *at, INF_ATOM_DATA *inf_at_data, int num_at, DRAW_PARMS *dp);
#endif

    /**
     * @brief Print the file name
     *
     * @param fmt Format string
     * @param out_file Output file pointer
     * @param szFname File name string
     */
    void PrintFileName(const char *fmt, FILE *out_file, const char *szFname);

    /**
     * @brief Sleep for a specified amount of time
     *
     * @param ms Milliseconds to sleep
     */
    void MySleep(unsigned long ms);

    /**
     * @brief Reconcile all CML bond parities in the input atoms
     *
     * @param at Pointer to the atom data
     * @param num_atoms Number of atoms
     * @param bDisconnected Flag for disconnected atoms
     * @return int Status code
     */
    int ReconcileAllCmlBondParities(inp_ATOM *at, int num_atoms, int bDisconnected);

#ifndef COMPILE_ALL_CPP
#ifdef __cplusplus
}
#endif
#endif

#endif /* _INPDEF_H_ */

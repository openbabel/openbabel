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

#ifndef _MOL_FMT_H_
#define _MOL_FMT_H_

#include <stdio.h>

#include "ichisize.h"
#include "mode.h" /* djb-rwth: necessary header file */

/*
    Data structures and constants
*/

/*************************** read MOL file V2000.************************/
/* ref: A.Dalby et al, "Description of Several Chemical Structure
 * File Formats Used by Computer Programs Developed at Molecular
 * Design Limited", J. Chem. Inf. Comput. Sci., 1992, 32, 244-255.
 */

/* ************************** read MOL file V3000.************************ */
/* http://download.accelrys.com/freeware/ctfile-formats/CTFile-formats.zip
 * Last accessed 2013-06-11
 */

/*-----------*/
/* CONSTANTS */
/*-----------*/

#define SD_FMT_END_OF_DATA "$$$$"

#define MOL_FMT_INPLINELEN 204 /* add cr, lf, double zero termination */
#ifndef MOL_FMT_MAXLINELEN
#define MOL_FMT_MAXLINELEN 200
#endif

#define MOL_FMT_PRESENT 1
#define MOL_FMT_ABSENT 0

/* configuration */
#define MOL_FMT_QUERY MOL_FMT_ABSENT
#define MOL_FMT_CPSS MOL_FMT_ABSENT
#define MOL_FMT_REACT MOL_FMT_ABSENT

#define MOL_FMT_STRING_DATA 'S'
#define MOL_FMT_CHAR_INT_DATA 'C'
#define MOL_FMT_SHORT_INT_DATA 'N'
#define MOL_FMT_LONG_INT_DATA 'L'
#define MOL_FMT_DOUBLE_DATA 'D'
#define MOL_FMT_FLOAT_DATA 'F'
#define MOL_FMT_JUMP_TO_RIGHT 'J'
#define MOL_FMT_INT_DATA 'I'

#define MOL_FMT_MAX_VALUE_LEN 32 /* max length of string containing a numerical value */

#define MOL_FMT_M_STY_NON 0 /**< None */
#define MOL_FMT_M_STY_SRU 1 /**< Structure repeating unit */
#define MOL_FMT_M_STY_MON 2 /**< Monomer */
#define MOL_FMT_M_STY_COP 3 /**< Copolymer */
#define MOL_FMT_M_STY_MOD 4 /**< Modification */
#define MOL_FMT_M_STY_CRO 5 /**< Crosslink */
#define MOL_FMT_M_STY_MER 6 /**< Mer type */

#define MOL_FMT_M_SST_NON 0 /**< None */
#define MOL_FMT_M_SST_ALT 1 /**< Alternating */
#define MOL_FMT_M_SST_RAN 2 /**< Random */
#define MOL_FMT_M_SST_BLK 3 /**< Block */

#define MOL_FMT_M_CONN_NON 0
#define MOL_FMT_M_CONN_HT 1
#define MOL_FMT_M_CONN_HH 2
#define MOL_FMT_M_CONN_EU 3

/* V3000 specific constants */
#define MOL_FMT_V3000_STENON -1
#define MOL_FMT_V3000_STEABS 1
#define MOL_FMT_V3000_STEREL 2
#define MOL_FMT_V3000_STERAC 3

/* provisional limits for V3000 */
#define MOL_FMT_V3000_INPLINELEN 32004 /* add cr, lf, double zero termination */
#ifndef MOL_FMT_V3000_MAXLINELEN
#define MOL_FMT_V3000_MAXLINELEN 32000
#endif
#define MOL_FMT_V3000_MAXFIELDLEN 4096

/*#ifdef TARGET_EXE_USING_API*/
#ifndef ISOTOPIC_SHIFT_FLAG
#define ISOTOPIC_SHIFT_FLAG 10000 /* add to isotopic mass if isotopic_mass = */
#endif
/*#endif*/

/*-------------------*/
/* SIMPLE DATA TYPES */
/*-------------------*/

#ifndef INCHI_US_CHAR_DEF
typedef signed char S_CHAR;
typedef unsigned char U_CHAR;
#define INCHI_US_CHAR_DEF
#endif

#ifndef LEN_COORD
#define LEN_COORD 10
#endif
#ifndef NUM_COORD
#define NUM_COORD 3
#endif

/*-----------------*/
/* DATA STRUCTURES */
/*-----------------*/
/**
 * @brief NUM_LISTS - Dynamically growing array of numeric lists
 * @param lists Pointer to the array of integer lists.
 * @param allocated Amount of memory allocated.
 * @param used Amount of memory used.
 * @param increment Amount to increment when expanding.
 */
typedef struct A_NUM_LISTS
{
    int **lists;
    int allocated;
    int used;
    int increment;
} NUM_LISTS;

/**
 * @brief Allocates memory for a specified number of numeric lists.
 * @param num_lists Pointer to the NUM_LISTS structure.
 * @param nlists Number of lists to allocate.
 * @return 0 on success,
 *         -1 on failure.
 */
int NumLists_Alloc(NUM_LISTS *num_lists, int nlists);

/**
 * @brief Reallocates memory for the numeric lists.
 * @param num_lists Pointer to the NUM_LISTS structure.
 * @return 0 on success,
 *         -1 on failure.
 */
int NumLists_ReAlloc(NUM_LISTS *num_lists);

/**
 * @brief Push new item to the end of array.
 * @param num_lists Pointer to the NUM_LISTS structure.
 * @param list Pointer to the list of integers to append.
 * @return 0 on success,
 *         -1 on failure.
 */
int NumLists_Append(NUM_LISTS *num_lists, int *list);

/**
 * @brief Frees the memory allocated for the numeric list.
 * @param num_lists Pointer to the NUM_LISTS structure.
 */
void NumLists_Free(NUM_LISTS *num_lists);

/**
 * @brief INT_ARRAY - Dynamically growing array of int
 * @param item Pointer to the array of integers.
 * @param allocated Amount of memory allocated.
 * @param used Amount of memory used.
 * @param increment Amount to increment when expanding.
 */
typedef struct tagINT_ARRAY
{
    int *item;
    int allocated;
    int used;
    int increment;
} INT_ARRAY;

/**
 * @brief Initializes an INT_ARRAY structure.
 *
 * @param items Pointer to the INT_ARRAY structure to initialize.
 * @param nitems Initial number of items to allocate.
 * @return 0 on success,
 *         -1 on failure.
 */
int IntArray_Alloc(INT_ARRAY *items, int nitems);

/**
 * @brief Reallocates memory for the INT_ARRAY structure.
 *
 * @param items Pointer to the INT_ARRAY structure to reallocate.
 * @return 0 on success,
 *         -1 on failure.
 */
int IntArray_ReAlloc(INT_ARRAY *items);

/**
 * @brief Appends a new item to the INT_ARRAY structure.
 *
 * @param items Pointer to the INT_ARRAY structure to modify.
 * @param new_item The new item to append.
 * @return 0 on success,
 *         -1 on failure.
 */
int IntArray_Append(INT_ARRAY *items, int new_item);

/**
 * @brief Appends a new item to the INT_ARRAY structure if it is not already present.
 *
 * @param items Pointer to the INT_ARRAY structure to modify.
 * @param new_item The new item to append.
 * @return 0 on success,
 *         -1 on failure.
 */
int IntArray_AppendIfAbsent(INT_ARRAY *items, int new_item);

/**
 * @brief Resets the INT_ARRAY structure to its initial state.
 *
 * @param items Pointer to the INT_ARRAY structure to reset.
 */
void IntArray_Reset(INT_ARRAY *items);

/**
 * @brief Frees the memory allocated for the INT_ARRAY structure.
 *
 * @param items Pointer to the INT_ARRAY structure to free.
 */
void IntArray_Free(INT_ARRAY *items);

/**
 * @brief Prints the contents of the INT_ARRAY structure for debugging purposes.
 *
 * @param items Pointer to the INT_ARRAY structure to print.
 */
void IntArray_DebugPrint(INT_ARRAY *items);

/**
 * @brief Data structure for Sgroup data (substance group data).
 *
 * @param id it is what is called 'Sgroup number' in CTFile
 * @param type Sgroup type: SUP = superatom, MUL = multiple group, SRU = SRU type, MON = monomer, MER = Mer type, COP = copolymer, CRO = crosslink, MOD = modification, GRA = graft, COM = component, MIX = mixture, FOR = formulation, DAT = data Sgroup, ANY = any polymer, GEN = generic
 * @param subtype (SST)
 * @param conn (SCN) Connectivity
 * @param label  what is called 'unique Sgroup identifier' in CTFile (SLB)
 * @param xbr1 bracket ends coordinates (SDI)
 * @param xbr2 bracket ends coordinates (SDI)
 * @param smt Sgroup Subscript (SMT)
 * @param alist list of atom indices (AL)
 * @param blist list of bond indices (BL)
 */
typedef struct A_MOL_FMT_SGROUP
{
    int id;
    int type;
    int subtype;
    int conn;
    int label;
    double xbr1[4];
    double xbr2[4];
    char smt[80];
    INT_ARRAY alist;
    INT_ARRAY blist;
} MOL_FMT_SGROUP;

/**
 * @brief  Allocate new array Sgroup.
 *
 * @param sgroup Pointer to the Sgroup structure to allocate.
 * @param id Sgroup ID.
 * @param type Sgroup type.
 * @return 0 on success, -1 on failure.
 */
int MolFmtSgroup_Create(MOL_FMT_SGROUP **sgroup,
                        int id,
                        int type);

/**
 * @brief Frees the memory allocated for the Sgroup.
 *
 * @param sgroup Pointer to the Sgroup structure to free.
 */
void MolFmtSgroup_Free(MOL_FMT_SGROUP *sgroup);

/**
 * @brief MOL_FMT_SGROUPS is a dynamically growing array of pointers to MOL_FMT_SGROUP objects.
 *
 * @param group Pointer to the growable array of pointers to MOL_FMT_SGROUPs
 * @param allocated Number of allocated objects
 * @param used Number of used objects
 * @param increment Array expansion increment
 */
typedef struct A_MOL_FMT_SGROUPS
{
    MOL_FMT_SGROUP **group;
    int allocated;
    int used;
    int increment;
} MOL_FMT_SGROUPS;

/**
 * @brief Allocates memory for a specified number of Sgroup objects.
 *
 * @param items Pointer to the MOL_FMT_SGROUPS structure.
 * @param nitems Number of Sgroup objects to allocate.
 * @return 0 on success, -1 on failure.
 */
int MolFmtSgroups_Alloc(MOL_FMT_SGROUPS *items, int nitems);

/**
 * @brief Expand array of Sgroups.
 *
 * @param items Pointer to the MOL_FMT_SGROUPS structure.
 * @return 0 on success, -1 on failure.
 */
int MolFmtSgroups_ReAlloc(MOL_FMT_SGROUPS *items);

/**
 * @brief Appends a new Sgroup to the array.
 *
 * @param items Pointer to the MOL_FMT_SGROUPS structure.
 * @param id Sgroup ID.
 * @param type Sgroup type.
 * @return 0 on success, -1 on failure.
 */
int MolFmtSgroups_Append(MOL_FMT_SGROUPS *items, int id, int type);

/**
 * @brief Frees the memory allocated for the MOL_FMT_SGROUPS list.
 *
 * @param items Pointer to the MOL_FMT_SGROUPS structure.
 */
void MolFmtSgroups_Free(MOL_FMT_SGROUPS *items);

/**
 * @brief Gets the index of a Sgroup by its ID.
 *
 * @param id Sgroup ID.
 * @param items Pointer to the MOL_FMT_SGROUPS structure.
 * @return Index of the Sgroup on success, -1 on failure.
 */
int MolFmtSgroups_GetIndexBySgroupId(int id, MOL_FMT_SGROUPS *items);

/**
 * @brief Data structure for MOL file header block (3 lines).
 *
 * @param molname Name of molecule (up to 80 characters)
 * @param line2 Second line of the header (the whole line2, up to 80 chars)
 * @param user_initls User initials (2 bytes; char)
 * @param prog_name Program name (8 bytes; char)
 * @param month Month (2 bytes; integral)
 * @param day Day (2 bytes; integral)
 * @param year Year (2 bytes; integral)
 * @param hour Hour (2 bytes; integral)
 * @param minute Minute (2 bytes; integral )
 * @param dim_code Dimensional code (2 bytes, dimensional code, char)
 * @param scaling_factor1 Scaling factor 1 (2 bytes, I2)
 * @param scaling_factor2 Scaling factor 2 (10 bytes, F10.5 )
 * @param energy Energy (10 bytes, F10.5  )
 * @param internal_regno Internal registration number (6 bytes, integral)
 * @param comment Comment (Line #3, up to 80 characters)
 */
typedef struct A_MOL_FMT_HEADER_BLOCK
{
    char molname[MOL_FMT_MAXLINELEN + 1];
    char line2[MOL_FMT_MAXLINELEN + 1];
    char user_initls[3];
    char prog_name[9];
    char month;
    char day;
    char year;
    char hour;
    char minute;
    char dim_code[3];
    short scaling_factor1;
    double scaling_factor2;
    double energy;
    long internal_regno;
    char comment[81];
} MOL_FMT_HEADER_BLOCK;

/**
 * @brief Data structure for atom representation in the MOL format.
 *
 * @param fx x coordinate: F10.5;       Generic
 * @param fy y coordinate: F10.5;       Generic
 * @param fz z coordinate: F10.5;       Generic
 * @param symbol Element symbol (aaa; up to 6 characters)
 * @param mass_difference Mass difference (dd;  (M_ISO); Generic: -3..+4 otherwise 0 or 127=most abund. isotope )
 * @param charge Formal charge (M CHG)
 * @param radical Radical status (M RAD)
 * @param stereo_parity Stereochemical parity (M STY)
 * @param H_count_plus_1 Hydrogen count plus one (hhh;         Query; Hn means >= n H; H0 means no H)
 * @param stereo_care bbb; Query: 0=ignore; 1=must match
 * @param valence vvv: <vvv: 0=no marking; (1..14)=(1..14); 15=zero valence.Number of bonds includes bonds to impl. H's>
 * @param H0_designator HHH: CPSS
 * @param reaction_component_type rrr: CPSS: 1=reactant, 2=product, 3=intermediate
 * @param reaction_component_num iii: CPSS: 0 to (n-1)
 * @param atom_atom_mapping_num mmm:        Reaction: 1..255
 * @param cInversionRetentionFlag nnn: 1=inverted, 2=retained config.; 0=property not applied
 * @param exact_change_flag eee
 * @param my_n_impH number of implicit H calculated for adding H to strings in STDATA
 * @param display_tom Do not hide element's name (applies to C 7-25-98 DCh)
 * @param atom_aliased_flag Do not remove charge/radical/isotope if it is in the alias. 9-3-99 DCh
 */
typedef struct A_MOL_FMT_ATOM
{
    double fx;      /* F10.5;       Generic             */
    double fy;      /* F10.5;       Generic             */
    double fz;      /* F10.5;       Generic             */
    char symbol[6]; /* aaa;         Generic             */

    S_CHAR mass_difference; /* dd;  (M_ISO)                     */
                            /*   Generic: -3..+4 otherwise 0 or */
                            /*   127=most abund. isotope        */
    S_CHAR charge;          /* ccc; (M CHG),                    */
                            /*  Generic: 1=+3, 2=+2,3=+1,       */
                            /*  4=doublet,5=-1,6=-2,7=-3        */
    char radical;           /*      (M RAD)                     */
    char stereo_parity;     /*   sss;         Generic           */
#if (MOL_FMT_QUERY == MOL_FMT_PRESENT)
    char H_count_plus_1; /* hhh;         Query;              */
                         /*   Hn means >= n H;               */
                         /*   H0 means no H                  */
    char stereo_care;    /* bbb;         Query: 0=ignore;    */
                         /*   1=must match                   */
#endif
    char valence; /* vvv:                             */
                  /*   0=no marking; (1..14)=(1..14); */
                  /*   15=zero valence.Number of bonds*/
                  /*   includes bonds to impl. H's    */

#if (MOL_FMT_CPSS == MOL_FMT_PRESENT)
    char H0_designator;           /* HHH:         CPSS                */
    char reaction_component_type; /* rrr:                             */
                                  /* CPSS: 1=reactant,                */
                                  /*       2=product,                 */
                                  /*       3=intermediate             */
    char reaction_component_num;  /* iii:         CPSS: 0 to (n-1)    */
#endif
#if (MOL_FMT_REACT == MOL_FMT_PRESENT)
    short atom_atom_mapping_num;  /* mmm:        Reaction: 1..255     */
    char cInversionRetentionFlag; /* nnn:                             */
                                  /*   1=inverted, 2=retained config.;*/
                                  /*   0=property not applied         */
#endif
#if (MOL_FMT_REACT == MOL_FMT_PRESENT || MOL_FMT_QUERY == MOL_FMT_PRESENT)
    char exact_change_flag; /* eee                              */
#endif
    char my_n_impH;         /* number of implicit H calculated  */
                            /* for adding H to strings in STDATA*/
    char display_tom;       /* Do not hide element's name       */
                            /*   (applies to C 7-25-98 DCh      */
    char atom_aliased_flag; /* Do not remove                    */
                            /* charge/radical/isotope if it     */
                            /*   is in the alias. 9-3-99 DCh    */
} MOL_FMT_ATOM;

/**
 * @brief Data structure for bond representation in the MOL format.
 *
 * @param atnum1 First atom number.
 * @param atnum2 Second atom number.
 * @param bond_type Type of bond (single, double, triple, etc.).
 * @param bond_stereo Stereo information for the bond.
 * @param bond_topology Topology information for the bond.
 * @param react_center_status Reaction center status information.
 */
typedef struct A_MOL_FMT_BOND
{
    short atnum1;     /* 111: First atom number: Generic  */
    short atnum2;     /* 222: Second atom number:Generic  */
    char bond_type;   /* ttt:                             */
                      /* 1,2,3=single, double, triple;    */
                      /* 4=aromatic;                      */
                      /* 5=single or double;              */
                      /* 6=single or aromatic;            */
                      /* 7=double or aromatic;            */
                      /* 8=any.                           */
                      /* Values 4-8 are for               */
                      /* SSS queries only                 */
    char bond_stereo; /* sss:                             */
                      /* Single bonds:                    */
                      /*   0=not stereo, 1=up,            */
                      /*   4=either, 6=down               */
                      /* Double bonds:                    */
                      /*   0=use x,y,z to                 */
                      /*     determine cis/trans,         */
                      /*   3=cis or trans (either)        */
                      /*   xxx:     not used              */

#if (MOL_FMT_QUERY == MOL_FMT_PRESENT)
    char bond_topology; /* rrr:                             */
                        /* 0=either, 1=ring, 2=chain:       */
                        /* SSS queries only                 */
#endif
#if (MOL_FMT_REACT == MOL_FMT_PRESENT)
    char react_center_status; /* ccc:                             */
                              /* 0 = unmarked,                    */
                              /* 1 = a center,                    */
                              /* -1 = not a center;               */
                              /* Additional:                      */
                              /* 2 = no charge,                   */
                              /* 4 = bond made/broken,            */
                              /* 8 = bond order changes           */
                              /*  12=4+8; 5=4+1, 9=8+1, 13=12+1   */
                              /*  12=4+8; 5=4+1, 9=8+1, 13=12+1   */
                              /*  are also possible               */
#endif
} MOL_FMT_BOND;

/**
 * @brief Data structure for V3000 representation in the MOL format.
 *
 * @param n_non_star_atoms Number of non-star atoms.
 * @param n_star_atoms Number of star atoms.
 * @param atom_index_orig Original atom indices as supplied.
 * @param atom_index_fin Final atom indices, with -1 for star atoms.
 * @param n_sgroups Number of S-groups.
 * @param n_3d_constraints Number of 3D constraints.
 * @param n_collections Number of collections.
 * @param n_non_haptic_bonds Number of non-haptic bonds.
 * @param n_haptic_bonds Number of haptic bonds.
 * @param haptic_bonds Pointer to the list of haptic bonds (int* contains bond type, non-star atom number, nendpts, then endpts themselves).
 * @param n_steabs Number of absolute stereo groups.
 * @param steabs Pointer to the list of absolute stereo groups (e.g. R and S).
 * @param n_sterel Number of relative stereo groups.
 * @param sterel Pointer to the list of relative stereo groups (OR - describes the orientation of groups relative to each other, such as in cis- and trans-isomers).
 * @param n_sterac Number of racemic stereo groups.
 * @param sterac Pointer to the list of racemic stereo groups (AND - equal 50:50 mixture of two enantiomers of a chiral molecule).
 */
typedef struct A_MOL_FMT_v3000
{
    int n_non_star_atoms;
    int n_star_atoms;
    int *atom_index_orig; /* index as supplied for atoms                      */
    int *atom_index_fin;  /* = index or -1 for star atom                      */
    int n_sgroups;        /* currently, we do not use this.                   */
    int n_3d_constraints; /* currently, we do not use this.                   */
    int n_collections;
    int n_non_haptic_bonds;
    int n_haptic_bonds;
    NUM_LISTS *haptic_bonds; /* haptic_bonds[i] is ptr to int* which contains    */
                             /* bond_type, non-star atom number,                 */
                             /* nendpts, then endpts themselves                  */
    /* Enhanced stereo */
    int n_steabs;
    NUM_LISTS *steabs; /* steabs[k][0] - not used                          */
                       /* steabs[k][1] -  number of members in collection  */
                       /* steabs[k][2..] - member atom numbers             */
    int n_sterel;
    NUM_LISTS *sterel; /* sterel[k][0] - n from "STERELn" tag              */
                       /* sterel[k][1] -  number of members in collection  */
                       /* sterel[k][2..] - member atom numbers             */
    int n_sterac;
    NUM_LISTS *sterac; /* sterac[k][0] - n from "STERACn" tag              */
                       /* sterac[k][1] -  number of members in collection  */
                       /* sterac[k][2..] - member atom numbers          */
} MOL_FMT_v3000;

/**
 * @brief Connection table data structure.
 *
 * @param n_atoms Number of atoms in the molecule.
 * @param n_bonds Number of bonds in the molecule.
 * @param n_atom_lists Number of atom lists.
 * @param chiral_flag Chiral flag indicating the presence of chiral centers.
 * @param n_stext_entries Number of stereo text entries.
 * @param n_reaction_components_plus_1 Number of reaction components.
 * @param n_reactants Number of reactants.
 * @param n_products Number of products.
 * @param n_intermediates Number of intermediates.
 * @param n_property_lines Number of property lines.
 * @param follow_inchi_1_treating_iso_mass Flag indicating whether to follow InChI-1 treating isotopic mass.
 * @param version_string Version string indicating the format version.
 * @param atoms Pointer to the array of atom block data structure.
 * @param bonds Pointer to the array of bond block data structure.
 * @param coords Pointer to the array of coordinate data structure.
 * @param sgroups Growable array of pointers to Sgroup objects.
 * @param v3000 Pointer to the V3000 specific data structure.
 */
typedef struct A_MOL_FMT_CTAB
{
    /* Line #1: Counts line */
    int n_atoms; /* int accounts for possible V3000. Was: aaa; <= 255; Generic */
    int n_bonds; /* int accounts for possible V3000. Was: bbb; <= 255; Generic */
#if (MOL_FMT_QUERY == MOL_FMT_PRESENT)
    short n_atom_lists; /* lll; <=  30; Query                               */
#endif

    char chiral_flag;      /* ccc; 0 or 1; Generic                             */
    short n_stext_entries; /* sss;         CPSS                                */
#if (MOL_FMT_CPSS == MOL_FMT_PRESENT)
    short n_reaction_components_plus_1; /* xxx;         CPSS                    */
    short n_reactants;                  /* rrr;         CPSS                                */
    short n_products;                   /* ppp;         CPSS                                */
    short n_intermediates;              /* iii;         CPSS                                */
#endif
    short n_property_lines; /* mmm;         Generic                             */
    short follow_inchi_1_treating_iso_mass;
    char version_string[7]; /* vvvvvv;      Generic; 'V2000'                    */
    MOL_FMT_ATOM *atoms;    /* The Atom Block                                   */
    MOL_FMT_BOND *bonds;
    MOL_COORD *coords;
    MOL_FMT_SGROUPS sgroups; /*    growable array of pointers to Sgroup objects  */
    MOL_FMT_v3000 *v3000;
} MOL_FMT_CTAB;

/**
 * @brief Data structure for a MOL file
 *
 * @param hdr Header block data structure.
 * @param ctab Connection table data structure.
 */
typedef struct A_MOL_FMT_DATA
{
    MOL_FMT_HEADER_BLOCK hdr;
    MOL_FMT_CTAB ctab;
} MOL_FMT_DATA;

/*
    Functions
*/

/**
 * @brief Reads header lines and connection table block from input SD or MOL file, ignore STEXT block, queries, and 3D features
 * @param inp_file Input file
 * @param OnlyHeaderBlock Output header data structure, return MOL_FMT_DATA will point to it
 * @param OnlyCTab Output connection table data structure, return MOL_FMT_DATA will point to it
 * @param bGetOrigCoord Flag to get original coordinates
 * @param treat_polymers Flag to treat polymers
 * @param err Error code
 * @param pStrErr Error string
 * @param bNoWarnings Flag to show warnings
 * @return MOL_FMT_DATA* returns mol file data structure, includes e.g. header block, connection table, ...
 */
static MOL_FMT_DATA* MolfileReadDataLines(INCHI_IOSTREAM* inp_file,
    MOL_FMT_HEADER_BLOCK* OnlyHeaderBlock,
    MOL_FMT_CTAB* OnlyCTab,
    int bGetOrigCoord,
    int treat_polymers,
    int* err, char* pStrErr, int bNoWarnings);

/**
 * @brief Reads header lines from input MOL file. A MOL file can have 3 header lines: (1) the name of the molecule, (2) details about the software used and (3) a comment line
 *
 * @param hdr Output header data structure
 * @param inp_file Input file
 * @param pStrErr Error string
 * @return * int Error code, retuns 0 - no error, 1 - error: can't read header block name, 3 - error: can't read header block 2 line, 7 - error: cant' read header block comment line
 */
static int MolfileReadHeaderLines(MOL_FMT_HEADER_BLOCK* hdr, INCHI_IOSTREAM* inp_file, char* pStrErr);

/**
 * @brief Reads counts line from input MOL file, includes information about the number of atoms, bonds, and atom lists, the chiral flag setting, and the Ctab version.
 *
 * @param ctab Connection table data structure
 * @param inp_file Input file
 * @param pStrErr Error string
 * @return int Error code, returns 0 - no error, 3 - error: can't read counts line, -1 - error: out of RAM
 */
static int MolfileReadCountsLine(MOL_FMT_CTAB* ctab, INCHI_IOSTREAM* inp_file, char* pStrErr);

/**
 * @brief Reads an atom block from input MOL file (V2000).
 * @param ctab Connection table data structure
 * @param inp_file Input file
 * @param err Error code
 * @param pStrErr Error string
 * @return int Error code, returns 0 - no error, 4 - error: can't interpret atom block, 5 - error: can't interpret second half of atom block?
 */
static int MolfileReadAtomsBlock(MOL_FMT_CTAB* ctab, INCHI_IOSTREAM* inp_file,
    int err, char* pStrErr);

/**
 * @brief Reads a bond block from input MOL file (V2000)
 *
 * @param ctab Connection table data structure
 * @param inp_file Input file
 * @param err Error code
 * @param pStrErr Error string
 * @return int Error code
 */
static int MolfileReadBondsBlock(MOL_FMT_CTAB* ctab, INCHI_IOSTREAM* inp_file,
    int err, char* pStrErr);

/**
 * @brief Reads a substance text block
 *
 * @param ctab Connection table data structure
 * @param inp_file Input file
 * @param err Error code
 * @param pStrErr Error string
 * @return int Error code
 */

/**
 * @brief Read MOL file data in to data structures.
 *
 * @param inp_file Pointer to the input file stream.
 * @param OnlyHeaderBlock Pointer to the header block to read.
 * @param OnlyCTab Pointer to the connection table to read.
 * @param bGetOrigCoord Flag indicating whether to get original coordinates.
 * @param treat_polymers Flag indicating whether to treat polymers.
 * @param pseudos_allowed Flag indicating whether pseudo atoms are allowed.
 * @param pname Pointer to the name buffer.
 * @param lname Length of the name buffer.
 * @param Id Pointer to the ID buffer.
 * @param pSdfLabel Pointer to the SDF label buffer.
 * @param pSdfValue Pointer to the SDF value buffer.
 * @param err Pointer to the error code.
 * @param pStrErr Pointer to the error string buffer.
 * @param bNoWarnings Flag indicating whether to suppress warnings.
 * @return MOL_FMT_DATA*
 */
MOL_FMT_DATA *ReadMolfile(INCHI_IOSTREAM *inp_file,
                          MOL_FMT_HEADER_BLOCK *OnlyHeaderBlock,
                          MOL_FMT_CTAB *OnlyCTab,
                          int bGetOrigCoord,
                          int treat_polymers,
                          int pseudos_allowed,
                          char *pname,
                          int lname,
                          unsigned long *Id,
                          const char *pSdfLabel,
                          char *pSdfValue,
                          int *err,
                          char *pStrErr,
                          int bNoWarnings);

/**
 * @brief Read a string from a source buffer into a destination buffer.
 *
 * @param dest Pointer to the destination buffer.
 * @param source Pointer to the source buffer.
 * @param len Length of the string to read.
 * @param first_space Pointer to a variable that will receive the address of the first space character in the source buffer.
 * @return number of actually processed bytes excluding zero terminator.
 */
int MolfileStrnread(char *dest,
                    char *source,
                    int len,
                    char **first_space);

/**
 * @brief Extract the 'data' in the MOL file field at given text position 'line_ptr'.
 *
 * @param data Pointer to the destination buffer.
 * @param field_len Length of the field to read. For MOL_FMT_STRING_DATA does not include
 *                  trailing zero, that is actual length of the string pointed by 'data'
 *                  should be at least field_len+1 bytes.
 *                  For numerical data 'field_len' is length of input data field
 *                  For numerical integral data field_len <= 0 means read up to first
 *                  non-numeric character as strtod() does ("free format")
 * @param data_type Type of the data to read. E.g. MOL_FMT_STRING_DATA, MOL_FMT_CHAR_INT_DATA, etc.
 * @param line_ptr Pointer to a variable that will receive the address of the next line.
 * @return for MOL_FMT_STRING_DATA: number of bytes excluding trailing zero
 *         for all others:  1=success; 0 = empty; -1 = error
 * @details on exit *line_ptr points to the next byte after the last entered
 */
int MolfileReadField(void *data,
                     int field_len,
                     int data_type,
                     char **line_ptr);

/**
 * @brief Extract the MOL file number from the header name line like "Structure #22"
 *
 * @param pHdr Pointer to the header block.
 * @return The structure number.
 */
long MolfileExtractStrucNum(MOL_FMT_HEADER_BLOCK *pHdr);

/**
 * @brief Check if the MOL file has no chemical structure.
 *
 * @param mfdata Pointer to the MOL file data structure.
 * @return 1 if no chemical structure, 0 otherwise.
 */
int MolfileHasNoChemStruc(MOL_FMT_DATA *mfdata);

/**
 * @brief Copy MOL-formatted data of SDF record or Molfile to another file
 *
 * @param inp_file Pointer to the input file stream.
 * @param fPtrStart File pointer position to start copying from.
 * @param fPtrEnd File pointer position to stop copying at.
 * @param outfile Pointer to the output file stream.
 * @param num Structure number to write in the header line.
 * @return last position of the output file stream.
 */
int MolfileSaveCopy(INCHI_IOSTREAM *inp_file,
                    long fPtrStart,
                    long fPtrEnd,
                    FILE *outfile,
                    long num);

/**
 * @brief Get xyz dimensionality and normalization factors in the MOL file
 *
 * @param mfdata Pointer to the MOL file data structure.
 * @param find_norm_factors Flag indicating whether to find normalization/scaling factors.
 * @param x0 Pointer to the minimum x-coordinate of the molecule.
 * @param y0 Pointer to the minimum y-coordinate of the molecule.
 * @param z0 Pointer to the minimum z-coordinate of the molecule.
 * @param xmin Pointer to the x-coordinate (is set to 0).
 * @param ymin Pointer to the y-coordinate (is set to 0).
 * @param zmin Pointer to the z-coordinate (is set to 0).
 * @param scaler Pointer to the scaling factor.
 * @param err Pointer to the error code.
 * @param pStrErr Pointer to the error string buffer.
 * @return number of dimensions: 0, 2, 3
 */
int MolfileGetXYZDimAndNormFactors(MOL_FMT_DATA *mfdata,
                                   int find_norm_factors,
                                   double *x0,
                                   double *y0,
                                   double *z0,
                                   double *xmin,
                                   double *ymin,
                                   double *zmin,
                                   double *scaler,
                                   int *err,
                                   char *pStrErr);

/**
 * @brief Free MOL file data structure.
 *
 * @param mfdata Pointer to the MOL file data structure to free.
 * @return NULL pointer.
 */
MOL_FMT_DATA *FreeMolfileData(MOL_FMT_DATA *mfdata);

/*
    V3000 Molfile
*/

/**
 * @brief Initialize V3000 connection table in MOL file data structure.
 *
 * @param ctab Pointer to the connection table data structure.
 * @param pStrErr Pointer to the error string buffer.
 * @return 0 on success, -1 on failure.
 */
int MolfileV3000Init(MOL_FMT_CTAB *ctab,
                     char *pStrErr);

/**
 *  @brief Read V3000 head (begin and counts line) of in MOL file data structure from input file stream.
 *
 * @param ctab Pointer to the connection table data structure.
 * @param inp_file Pointer to the input file stream.
 * @param pStrErr Pointer to the error string buffer.
 * @return 0 on success, -1 on failure.
 */
int MolfileV3000ReadCTABBeginAndCountsLine(MOL_FMT_CTAB *ctab,
                                           INCHI_IOSTREAM *inp_file,
                                           char *pStrErr);

/**
 * @brief Read V3000 atoms block in MOL file data structure from input file stream.
 *
 * @param ctab Pointer to the connection table data structure.
 * @param inp_file Pointer to the input file stream.
 * @param err Error code.
 * @param pStrErr Pointer to the error string buffer.
 * @return 0 on success, -1 on failure.
 */
int MolfileV3000ReadAtomsBlock(MOL_FMT_CTAB *ctab,
                               INCHI_IOSTREAM *inp_file,
                               int err,
                               char *pStrErr);

/**
 * @brief Read V3000 bonds block in MOL file data structure from input file stream.
 *
 * @param ctab Pointer to the connection table data structure.
 * @param inp_file Pointer to the input file stream.
 * @param err Error code.
 * @param pStrErr Pointer to the error string buffer.
 * @return 0: Success (no error, bonds block read correctly),
 *         1: Error: No V3000 Bond block start marker ("BEGIN BOND" missing),
 *         2: Error: Cannot read V3000 bond block line (input line missing or unreadable),
 *         4: Error: Cannot interpret V3000 bond block line (parsing failure)
 *         Negative values: If the end-of-data marker ($$$$) is encountered, err is set to -abs(err) (e.g., -1, -2, -4)
 *         Other values: If the end marker ("END BOND") is missing, err is set to 1
 */
int MolfileV3000ReadBondsBlock(MOL_FMT_CTAB *ctab,
                               INCHI_IOSTREAM *inp_file,
                               int err,
                               char *pStrErr);

/**
 * @brief Read V3000 tail (haptic bonds, stereo collections, Sgroups, 3D constraints, collections, end line) in to the MOL file data structure from input file stream.
 *
 * @param ctab Pointer to the connection table data structure.
 * @param inp_file Pointer to the input file stream.
 * @param err Error code.
 * @param pStrErr Pointer to the error string buffer.
 * @return 0: Success (no error encountered, tail of CTAB read correctly)
 *         1: Error (e.g., missing "END CTAB" marker, or other parsing errors)
 *         7: Error in reading or interpreting V3000 collection lines
 *         71, 77, etc.: Error codes from sub-blocks (e.g., SGroup, 3DBlock, Collection) plus 70, indicating which sub-block failed
 *         Other positive values: Error codes from called functions (e.g., MolfileV3000ReadSGroup, MolfileV3000Read3DBlock, MolfileV3000ReadCollections) plus 70
 *         Negative values: If an end-of-data marker ($$$$) is encountered during an error, the error code may be returned as a negative value (e.g., -1, -7, -71, etc.)
 */
int MolfileV3000ReadTailOfCTAB(MOL_FMT_CTAB *ctab,
                               INCHI_IOSTREAM *inp_file,
                               int err,
                               char *pStrErr);

/**
 * @brief Read V3000 haptic bond information from the input line.
 *
 * @param ctab Pointer to the connection table data structure.
 * @param line_ptr Pointer to the line buffer.
 * @param num_list Pointer to the list of numbers.
 * @param pStrErr Pointer to the error string buffer.
 * @return Positive integer (nread > 0): Success; number of bytes/fields read from the haptic bond block.
 *         0: No data read (rare, usually means nothing was parsed).
 *        -1: Error; parsing failed at any step (missing '(', invalid count, allocation failure, field read error, missing "ATTACH=ALL", etc.).
 */
int MolfileV3000ReadHapticBond(MOL_FMT_CTAB *ctab,
                               char **line_ptr,
                               int **num_list,
                               char *pStrErr);

/**
 * @brief Read V3000 stereo collection information from the input line.
 *
 * @param ctab Pointer to the connection table data structure.
 * @param line_ptr Pointer to the line buffer.
 * @param num_list Pointer to the list of numbers.
 * @param pStrErr Pointer to the error string buffer.
 * @return Positive integer (nread > 0): Success; number of bytes/fields read from the stereo collection block.
 *         0: No data read (rare, usually means nothing was parsed).
 *        -1: Error; parsing failed at any step (missing '(', invalid count, allocation failure, field read error, missing "ATTACH=ALL", etc.).
 */
int MolfileV3000ReadStereoCollection(MOL_FMT_CTAB *ctab,
                                     char **line_ptr,
                                     int **num_list,
                                     char *pStrErr);

/**
 * @brief Read V3000 Sgroup information from the input file stream into the MOL file data structure.
 *
 * @param ctab Pointer to the connection table data structure.
 * @param inp_file Pointer to the input file stream.
 * @param err Error code.
 * @param pStrErr Pointer to the error string buffer.
 * @return 0: Success (no error encountered, SGroup read correctly)
 *         1: Error (e.g., missing "END SGROUP" marker, or other parsing errors)
 */
int MolfileV3000ReadSGroup(MOL_FMT_CTAB *ctab,
                           INCHI_IOSTREAM *inp_file,
                           int err,
                           char *pStrErr);

/**
 * @brief Read V3000 3D constraints block in MOL file data structure from input file stream.
 *
 * @param ctab Pointer to the connection table data structure.
 * @param inp_file Pointer to the input file stream.
 * @param err Error code.
 * @param pStrErr Pointer to the error string buffer.
 * @return 0: Success (no error encountered, 3D block read correctly)
 *         1: Error (e.g., missing "END 3D" marker, or other parsing errors)
 */
int MolfileV3000Read3DBlock(MOL_FMT_CTAB *ctab,
                            INCHI_IOSTREAM *inp_file,
                            int err,
                            char *pStrErr);

/**
 * @brief Read V3000 collections block in MOL file data structure from input file stream.
 *
 * @param ctab Pointer to the connection table data structure.
 * @param inp_file Pointer to the input file stream.
 * @param err Error code.
 * @param pStrErr Pointer to the error string buffer.
 * @return 0: Success (collections block read and parsed correctly)
 *         7: Error in reading or interpreting V3000 collection lines (parsing failure, unexpected format, missing required fields, etc.)
 */
int MolfileV3000ReadCollections(MOL_FMT_CTAB *ctab,
                                INCHI_IOSTREAM *inp_file,
                                int err,
                                char *pStrErr);
/*    Clean V3000 stuff */

/**
 * @brief Free memory allocated for V3000 specific data in the MOL file data structure.
 *
 * @param v3000 Pointer to the V3000 specific data structure to free.
 * @return 0.
 */
int DeleteMolfileV3000Info(MOL_FMT_v3000 *v3000);

/**
 * @brief Read a line from the input file stream, handling V3000 line continuations: Extended version of inchi_fgetsLf which is able of reading
 *        concatenated lines (ending with '-') of V3000 Molfile. Also removes "M  V30 " prefix" and normalizes the rest of string.
 *
 * @param line Pointer to the buffer to store the read line.
 * @param inp_stream Pointer to the input file stream.
 * @return Pointer to the read line, or NULL on end of file or error.
 */
char *inchi_fgetsLf_V3000(char *line,
                          INCHI_IOSTREAM *inp_stream);

/**
 * @brief Get a V3000 input line and store it in a string buffer.
 *
 * @param buf Pointer to the string buffer to store the line.
 * @param inp_stream Pointer to the input file stream.
 * @return Length of buffer stored on success, -1 on failure.
 */
int get_V3000_input_line_to_strbuf(INCHI_IOS_STRING *buf,
                                   INCHI_IOSTREAM *inp_stream);

/**
 * @brief Extract the 'data' in specified mol file field at given text position 'line_ptr'
 *
 * @param data Pointer to the destination buffer.
 * @param data_type Type of the data to read. E.g. MOL_FMT_STRING_DATA, MOL_FMT_CHAR_INT_DATA, etc.
 * @param line_ptr Pointer to the line buffer.
 * @return for MOL_FMT_STRING_DATA: number of bytes excluding trailing zero
 *         for all others:  1=success; 0 = empty
 */
int MolfileV3000ReadField(void *data,
                          int data_type,
                          char **line_ptr);
/**
 * @brief Read keyword from the specified line.
 *
 * @param key Pointer to the buffer to store the keyword.
 * @param line_ptr Pointer to the line buffer.
 * @return Length of the keyword on success, -1 on failure.
 */
int MolfileV3000ReadKeyword(char *key,
                            char **line_ptr);

/*
    SDF
*/

/**
 * @brief Skip extra data in SDF file after the MOL file data.
 *
 * @param inp_file Pointer to the input file stream.
 * @param CAS_num Pointer to the variable to store the CAS number.
 * @param comment Pointer to the buffer to store the comment.
 * @param lcomment Length of the comment buffer.
 * @param name Pointer to the buffer to store the name.
 * @param lname Length of the name buffer.
 * @param prev_err Previous error code.
 * @param pSdfLabel Pointer to the SDF label buffer.
 * @param pSdfValue Pointer to the SDF value buffer.
 * @param pStrErr Pointer to the error string buffer.
 * @param bNoWarnings Flag indicating whether to suppress warnings.
 * @return 0: Success; extra SDF data was skipped without errors.
 *         3: Unexpected SData header line (unexpected contents or format).
 *         5: Only blank lines encountered (all lines empty).
 *         9: Error occurred, but successfully bypassed to the next structure ($$$$ marker found); non-fatal, warning issued.
 *         Other values: The function may propagate the value of prev_err if set before entering the loop.
 */
int SDFileSkipExtraData(INCHI_IOSTREAM *inp_file,
                        unsigned long *CAS_num,
                        char *comment,
                        int lcomment,
                        char *name,
                        int lname,
                        int prev_err,
                        const char *pSdfLabel,
                        char *pSdfValue,
                        char *pStrErr,
                        int bNoWarnings);

/**
 * @brief Identify if the given line matches the specified SDF label.
 *
 * @param inp_line Pointer to the input line.
 * @param pSdfLabel Pointer to the SDF label to match.
 * @return SDF_DATA_HEADER_USER: The label matches the user-specified label (pSdfLabel).
 *         SDF_DATA_HEADER_NAME: The label is "NAME".
 *         SDF_DATA_HEADER_COMMENT: The label is "COMMENT".
 *         SDF_DATA_HEADER_CAS: The label is "CAS".
 *         SDF_DATA_HEADER: The label does not match any of the above (default case).
 */
int SDFileIdentifyLabel(char *inp_line, const char *pSdfLabel);

/**
 * @brief Extract CAS number from the given line.
 *
 * @param line Pointer to the line containing the CAS number.
 * @return The extracted CAS number, or 0 if not found or invalid.
 */
unsigned long SDFileExtractCASNo(char *line);

#endif /* _MOL_FMT_H_ */

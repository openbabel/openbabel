/*
 * International Chemical Identifier (InChI)
 * Version 1
 * Software version 1.02
 * October 31, 2008
 * Developed at NIST
 *
 * The InChI library and programs are free software developed under the
 * auspices of the International Union of Pure and Applied Chemistry (IUPAC);
 * you can redistribute this software and/or modify it under the terms of 
 * the GNU Lesser General Public License as published by the Free Software 
 * Foundation:
 * http://www.opensource.org/licenses/lgpl-license.php
 */


#ifndef __READMOL_H__
#define __READMOL_H__

/*************** read MOL file V2000.************************/
/* ref: A.Dalby et al, "Description of Several Chemical Structure
 * File Formats Used by Computer Programs Developed at Molecular
 * Design Limited", J. Chem. Inf. Comput. Sci., 1992, 32, 244-255.
 */

/*
#define MOLFILEINPLINELEN   84  // add cr, lf, double zero termination
#ifndef MOLFILEMAXLINELEN
#define MOLFILEMAXLINELEN   80
#endif
*/

#define MOLFILEINPLINELEN   204  /* add cr, lf, double zero termination */
#ifndef MOLFILEMAXLINELEN
#define MOLFILEMAXLINELEN   200
#endif


#define MOL_PRESENT 1
#define MOL_ABSENT  0

/* configuration */
#define MOL_QUERY   MOL_ABSENT
#define MOL_CPSS    MOL_ABSENT
#define MOL_REACT   MOL_ABSENT

#define MOL_STRING_DATA      'S'
#define MOL_CHAR_INT_DATA    'C'
#define MOL_SHORT_INT_DATA   'N'
#define MOL_LONG_INT_DATA    'L'
#define MOL_DOUBLE_DATA      'D'
#define MOL_FLOAT_DATA       'F'
#define MOL_JUMP_TO_RIGHT    'J'
#define MOL_MAX_VALUE_LEN    32    /* max length of string containing a numerical value */


#define SDF_END_OF_DATA "$$$$"

/****************************************************************************/
typedef struct tagMOL_HEADER_BLOCK {
    /* Line #1 */
    char    szMoleculeName[MOLFILEMAXLINELEN+1]; /* up to 80 characters */
    /* Line #2: optional */
    char    szMoleculeLine2[MOLFILEMAXLINELEN+1]; /* the whole line2 -- up to 80 characters */
    char    szUserInitials[3];      /* 2 bytes; char */
    char    szProgramName[9];       /* 8 bytes; char */
    char    cMonth;                 /* 2 bytes; integral */
    char    cDay;                   /* 2 bytes; integral */
    char    cYear;                  /* 2 bytes; integral */
    char    cHour;                  /* 2 bytes; integral */
    char    cMinute;                /* 2 bytes; integral */
    char    szDimCode[3];           /* 2 bytes: dimensional code; char */
    short   nScalingFactor1;        /* 2 bytes;  I2    */
    double  dScalingFactor2;        /* 10 bytes, F10.5 */
    double  dEnergy;                /* 10 bytes, F10.5 */
    long    lInternalRegistryNumber;/* 6 bytes, integral */
    /* Line #3: comment */
    char    szComment[81];
}MOL_HEADER_BLOCK;

/****************************************************************************/
typedef struct tagMOL_ATOM {
    double   fX;                       /* F10.5;       Generic */
    double   fY;                       /* F10.5;       Generic */
    double   fZ;                       /* F10.5;       Generic */
    char     szAtomSymbol[6];          /* aaa;         Generic */ /* changed from 4 to 6 to match STDATA */
    S_CHAR   cMassDifference;   /* dd;  (M_ISO) Generic: -3..+4 otherwise 0 or 127=most abund. isotope */
    S_CHAR   cCharge;           /* ccc; (M CHG), Generic: 1=+3,2=+2,3=+1,4=doublet,5=-1,6=-2,7=-3 */
    char     cRadical;                 /*      (M RAD) */
    char     cStereoParity;            /* sss;         Generic */
#if ( MOL_QUERY == MOL_PRESENT )
    char     cH_countPlus1;            /* hhh;         Query; Hn means >= n H; H0 means no H */
    char     cStereoCare;              /* bbb;         Query: 0=ignore; 1=must match */
#endif
    char     cValence;                 /* vvv:         0=no marking; (1..14)=(1..14); 15=zero valence */
                                       /*              number of bonds including bonds to implies H's */
#if ( MOL_CPSS == MOL_PRESENT )
    char     cH0_designator;           /* HHH:         CPSS */
    char     cReactionComponentType;   /* rrr:         CPSS: 1=reactant, 2=product, 3=intermediate */
    char     cReactionComponentNumber; /* iii:         CPSS: 0 to (n-1) */
#endif

#if ( MOL_REACT == MOL_PRESENT )
    short    nAtomAtomMappingNumber;   /* mmm:         Reaction: 1..255 */
    char     cInversionRetentionFlag;  /* nnn:         1=inverted,2=retained config.; 0=property not applied */
#endif
#if ( MOL_REACT == MOL_PRESENT || MOL_QUERY == MOL_PRESENT )
    char     cExactChargeFlag;         /* eee:         1=charge on atom must match exactly, 0=property not applied */
#endif
    char cMyNumImpH;                 /* number of implicit H calculated for adding H to strings in STDATA */
    char cDisplayAtom;               /* Do not hide element's name ( applies to C 7-25-98 DCh */
    char cAtomAliasedFlag;           /* Do not remove charge/radical/isotope if it is in the alias. 9-3-99 DCh */
} MOL_ATOM;

/****************************************************************************/
typedef struct tagMOL_BONDS {
    short nAtomNo1;                  /* 111:     First atom number:  Generic */
    short nAtomNo2;                  /* 222:     Second atom number: Generic */
    char  cBondType;                 /* ttt:     1,2,3=single, double, triple; 4=aromatic; 5=single or double */
                                     /*          6=single or aromatic, 7=double or aromatic, 8=any.           */
                                     /*          values 4-8 are for SSS queries only */
    char cBondStereo;                /* sss:     Single bonds: 0=not stereo, 1=up, 4=either, 6=down           */
                                     /*          Double bonds: 0=use x,y,z to determine cis/trans             */
                                     /*                        3=cis or trans (either)                        */
                                     /* xxx:     not used */
#if ( MOL_QUERY == MOL_PRESENT )
    char cBondTopology;              /* rrr:     0=either, 1=ring, 2=chain: SSS queries only */
#endif
#if ( MOL_REACT == MOL_PRESENT )
    char cReactingCenterStatus;      /* ccc:     0=unmarked, 1=a center, -1=not a center; Additional: */
                                     /*          2=no charge,4=bond made/broken,8=bond order changes  */
                                     /*          12=4+8; 5=4+1, 9=8+1, 13=12+1 are also possible      */
#endif
} MOL_BONDS;
/****************************************************************************/
typedef struct tagMOL_CTAB {
    /* Line #1: Counts line */
    short nNumberOfAtoms;                  /* aaa; <= 255; Generic */
    short nNumberOfBonds;                  /* bbb; <= 255; Generic */
#if ( MOL_QUERY == MOL_PRESENT )
    short nNumberOfAtomsLists;              /* lll; <=  30; Query   */
#endif
                                            /* fff; Obsolete */
    char cChiralFlag;                       /* ccc; 0 or 1; Generic */
    short nNumberOfStextEntries;            /* sss;         CPSS    */
#if ( MOL_CPSS == MOL_PRESENT )
    short nNumberOfReactionComponentsPlus1; /* xxx;         CPSS    */
    short nNumberOfReactants;               /* rrr;         CPSS    */
    short nNumberOfProducts;                /* ppp;         CPSS    */
    short nNumberOfIntermediates;           /* iii;         CPSS    */
#endif
    short nNumberOfPropertyLines;            /* mmm;         Generic */
    char  csCurrentCtabVersion[7];          /* vvvvvv;      Generic; 'V2000' */
    /* The Atom Block */
    MOL_ATOM  *MolAtom;
    MOL_BONDS *MolBond;
    MOL_COORD *szCoord;
} MOL_CTAB;

typedef struct tagMOL_DATA {
    MOL_HEADER_BLOCK hdr;
    MOL_CTAB         ctab;
} MOL_DATA;

#ifndef INCHI_ALL_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif


MOL_DATA* delete_mol_data( MOL_DATA* mol_data );
MOL_DATA* read_sdfile_segment(FILE* inp, MOL_HEADER_BLOCK *OnlyHeaderBlock, MOL_CTAB *OnlyCtab,
                              int bGetOrigCoord,
                              char *pname, int lname,
                              long *Id, const char *pSdfLabel, char *pSdfValue,
                              int *err, char *pStrErr );

#ifndef INCHI_ALL_CPP
#ifdef __cplusplus
}
#endif
#endif

#endif /*__READMOL_H__*/

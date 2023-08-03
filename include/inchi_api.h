/*
 * International Chemical Identifier (InChI)
 * Version 1
 * Software version 1.04
 * September 9, 2011
 *
 * The InChI library and programs are free software developed under the
 * auspices of the International Union of Pure and Applied Chemistry (IUPAC).
 * Originally developed at NIST. Modifications and additions by IUPAC 
 * and the InChI Trust.
 *
 * IUPAC/InChI-Trust Licence No.1.0 for the 
 * International Chemical Identifier (InChI) Software version 1.04
 * Copyright (C) IUPAC and InChI Trust Limited
 * 
 * This library is free software; you can redistribute it and/or modify it 
 * under the terms of the IUPAC/InChI Trust InChI Licence No.1.0, 
 * or any later version.
 * 
 * Please note that this library is distributed WITHOUT ANY WARRANTIES 
 * whatsoever, whether expressed or implied.  See the IUPAC/InChI Trust 
 * Licence for the International Chemical Identifier (InChI) Software 
 * version 1.04, October 2011 ("IUPAC/InChI-Trust InChI Licence No.1.0") 
 * for more details.
 * 
 * You should have received a copy of the IUPAC/InChI Trust InChI 
 * Licence No. 1.0 with this library; if not, please write to:
 * 
 * The InChI Trust
 * c/o FIZ CHEMIE Berlin
 *
 * Franklinstrasse 11
 * 10587 Berlin
 * GERMANY
 *
 * or email to: ulrich@inchi-trust.org.
 * 
 */


#ifndef __INHCH_API_H__
#define __INHCH_API_H__




/*^^^ Post-1.02b fix - thanks to David Foss */
#ifndef FIND_RING_SYSTEMS
#define FIND_RING_SYSTEMS 1
#endif
#ifndef FIND_RINS_SYSTEMS_DISTANCES
#define FIND_RINS_SYSTEMS_DISTANCES 0
#endif


/* radical definitions */
typedef enum tagINCHIRadical {
   INCHI_RADICAL_NONE    = 0,
   INCHI_RADICAL_SINGLET = 1,
   INCHI_RADICAL_DOUBLET = 2,
   INCHI_RADICAL_TRIPLET = 3
} inchi_Radical;

/* bond type definitions */
typedef enum tagINCHIBondType {
   INCHI_BOND_TYPE_NONE    =  0,
   INCHI_BOND_TYPE_SINGLE  =  1,
   INCHI_BOND_TYPE_DOUBLE  =  2,
   INCHI_BOND_TYPE_TRIPLE  =  3,
   INCHI_BOND_TYPE_ALTERN  =  4  /* avoid by all means */
} inchi_BondType;
/* 2D stereo definitions */
typedef enum tagINCHIBondStereo2D {
   /* stereocenter-related; positive: the sharp end points to this atom  */
   INCHI_BOND_STEREO_NONE           =  0,
   INCHI_BOND_STEREO_SINGLE_1UP     =  1,
   INCHI_BOND_STEREO_SINGLE_1EITHER =  4,
   INCHI_BOND_STEREO_SINGLE_1DOWN   =  6,
   /* stereocenter-related; negative: the sharp end points to the opposite atom  */
   INCHI_BOND_STEREO_SINGLE_2UP     = -1,
   INCHI_BOND_STEREO_SINGLE_2EITHER = -4,
   INCHI_BOND_STEREO_SINGLE_2DOWN   = -6,
   /* stereobond-related */
   INCHI_BOND_STEREO_DOUBLE_EITHER  =  3 /* unknown stereobond geometry */
} inchi_BondStereo2D;

/*************************************************************************
 * Notes on using INCHI_BOND_STEREO_SINGLE_*  from inchi_BondStereo2D    *
 *                                                                       *
 * These stereo markings are used by InChI to characterize a stereogenic *
 * atom if and only if all neighbors of this atom have same z-coordinate *
 * as this atom (that is, in case of 2D fragment).                       *
 * The only exception is INCHI_BOND_STEREO_SINGLE_?EITHER marking which  *
 * always assigns to the atom an "unknown" parity (u).                   *
 *                                                                       *
 * Note the behavior which is default for InChI software v.1.04/03/02std *
 * (at -NEWPSOFF option is not supplied) 2D stereo interpretation:       *
 * only bonds that have sharp end pointing to the stereogenic atom are   *
 * considered as being out of plane and only sharp ends of               *
 * INCHI_BOND_STEREO_SINGLE_?EITHER bonds are considered to determine    *
 * whether the stereochemistry is unknown.                               *
 *************************************************************************/

/* sizes definitions */
#define MAXVAL                   20 /* max number of bonds per atom                 */
#define ATOM_EL_LEN               6 /* length of ASCIIZ element symbol field        */
#define NUM_H_ISOTOPES            3 /* number of hydrogen isotopes: protium, D, T   */
#define ISOTOPIC_SHIFT_FLAG   10000 /* add to isotopic mass if isotopic_mass =      */
                                    /* (isotopic mass - average atomic mass)        */
#define ISOTOPIC_SHIFT_MAX      100 /* max abs(isotopic mass - average atomic mass) */

#ifndef INCHI_US_CHAR_DEF
typedef signed char   S_CHAR;
typedef unsigned char U_CHAR;
#define INCHI_US_CHAR_DEF
#endif

#ifndef INCHI_US_SHORT_DEF
typedef signed short   S_SHORT;
typedef unsigned short U_SHORT;
#define INCHI_US_SHORT_DEF
#endif

typedef  S_SHORT AT_NUM;            /* atom number; starts from 0 */

/*************************************************
 *
 *
 *  A T O M S   a n d   C O N N E C T I V I T Y
 *
 *
 *************************************************/

typedef struct tagInchiAtom {
    /* atom coordinates */
    double x;
    double y;
    double z;
    /* connectivity */
    AT_NUM  neighbor[MAXVAL];     /* adjacency list: ordering numbers of */
                                  /*            the adjacent atoms, >= 0 */
    S_CHAR  bond_type[MAXVAL];    /* inchi_BondType */
    /* 2D stereo */
    S_CHAR  bond_stereo[MAXVAL];  /* inchi_BondStereo2D; negative if the */
                                  /* sharp end points to opposite atom */
    /* other atom properties */
    char    elname[ATOM_EL_LEN];  /* zero-terminated chemical element name:*/
                                  /* "H", "Si", etc. */
    AT_NUM  num_bonds;            /* number of neighbors, bond types and bond*/
                                  /* stereo in the adjacency list */
    S_CHAR  num_iso_H[NUM_H_ISOTOPES+1]; /* implicit hydrogen atoms */
                                  /* [0]: number of implicit non-isotopic H
                                       (exception: num_iso_H[0]=-1 means INCHI
                                       adds implicit H automatically),
                                     [1]: number of implicit isotopic 1H (protium),
                                     [2]: number of implicit 2H (deuterium),
                                     [3]: number of implicit 3H (tritium) */
    AT_NUM  isotopic_mass;        /* 0 => non-isotopic; isotopic mass or  */
                                  /* ISOTOPIC_SHIFT_FLAG + mass - (average atomic mass) */
    S_CHAR  radical;              /* inchi_Radical */
    S_CHAR  charge;               /* positive or negative; 0 => no charge */
}inchi_Atom;

/*******************************************************************
 * Notes: 1. Atom ordering numbers (i, k, and atom[i].neighbor[j] below)
 *           start from zero; max. ordering number is (num_atoms-1).
 *        2. inchi_Atom atom[i] is connected to the atom[atom[i].neighbor[j]]
 *           by a bond that has type atom[i].bond_type[j] and 2D stereo type
 *           atom[i].bond_stereo[j] (in case of no stereo
 *           atom[i].bond_stereo[j] = INCHI_BOND_STEREO_NONE)
 *           Index j is in the range 0 <= j <= (atom[i].num_bonds-1)
 *        3. Any connection (represented by atom[i].neighbor[j],
 *           atom[i].bond_type[j], and atom[i].bond_stereo[j])
 *           should be present in one or both adjacency list:
 *             if k = atom[i].neighbor[j] then i may or may not be present in
 *           atom[k].neighbor[] list. For example, the adjacency lists may be
 *           populated with only such neighbors that atom[i].neighbor[j] < i
 *           All elements of an adjacency list must be different, that is,
 *           a bond must be specified in an adjacency list only once.
 *        4. in Molfiles usually
 *           (number of implicit H) = Valence - SUM(bond_type[])
 *        5. Seemingly illogical order of the inchi_Atom members was
 *           chosen in an attempt to avoid alignment problems when
 *           accessing inchi_Atom from unrelated to C programming
 *           languages such as Visual Basic.
 *******************************************************************/

/*******************************************************************
    0D Stereo Parity and Type definitions
 *******************************************************************
            Note:
            =====
            o Below #A is the ordering number of atom A, starting from 0
            o See parity values corresponding to 'o', 'e', and 'u' in
              inchi_StereoParity0D definition below)

           =============================================
            stereogenic bond >A=B< or cumulene >A=C=C=B<
           =============================================

                                 neighbor[4]  : {#X,#A,#B,#Y} in this order
     X                           central_atom : NO_ATOM
      \            X      Y      type         : INCHI_StereoType_DoubleBond
       A==B         \    /
           \         A==B
            Y

    parity= 'e'    parity= 'o'   unknown parity = 'u'

    Limitations:
    ============
    o Atoms A and B in cumulenes MUST be connected by a chain of double bonds;
      atoms A and B in a stereogenic 'double bond' may be connected by a double,
      single, or alternating bond.
    o One atom may belong to up to 3 stereogenic bonds (i.g. in a fused
      aromatic structure).
    o Multiple stereogenic bonds incident to any given atom should
      either all except possibly one have (possibly different) defined
      parities ('o' or 'e') or should all have an unknown parity 'u'.

      Note on parities of alternating stereobonds
      ===========================================
                                                     D--E
      In large rings  (see Fig. 1, all              //   \\
      atoms are C) all alternating bonds         B--C      F--G
      are treated as stereogenic.              //              \\
      To avoid "undefined" bond parities      A                  H
      for bonds BC, DE, FG, HI, JK, LM, AN     \               /
      it is recommended to mark them with       N==M       J==I
      parities.                                     \     /
                                                      L==K    Fig. 1
      Such a marking will make
      the stereochemical layer unambiguous
      and it will be different from the          B--C      F--G
      stereochemical layer of the second       //   \\   //    \\
      structure (Fig. 2).                     A      D--E        H
                                               \               /
                                                N==M       J==I
      By default, double and alternating            \     /
      bonds in 8-member and greater rings             L==K    Fig. 2
      are treated by InChI as stereogenic.


           =============================================
            tetrahedral atom
           =============================================

   4 neighbors

            X                    neighbor[4] : {#W, #X, #Y, #Z}
            |                    central_atom: #A
         W--A--Y                 type        : INCHI_StereoType_Tetrahedral
            |
            Z
   parity: if (X,Y,Z) are clockwize when seen from W then parity is 'e' otherwise 'o'
   Example (see AXYZW above): if W is above the plane XYZ then parity = 'e'

   3 neighbors

              Y          Y       neighbor[4] : {#A, #X, #Y, #Z}
             /          /        central_atom: #A
         X--A  (e.g. O=S   )     type        : INCHI_StereoType_Tetrahedral
             \          \
              Z          Z

   parity: if (X,Y,Z) are clockwize when seen from A then parity is 'e',
                                                          otherwise 'o'
   unknown parity = 'u'
   Example (see AXYZ above): if A is above the plane XYZ then parity = 'e'
   This approach may be used also in case of an implicit H attached to A.

           =============================================
            allene
           =============================================

       X       Y                 neighbor[4]  : {#X,#A,#B,#Y}
        \     /                  central_atom : #C
         A=C=B                   type         : INCHI_StereoType_Allene

                                      Y      X
                                      |      |
     when seen from A along A=C=B:  X-A    Y-A

                          parity:   'e'    'o'

   parity: if A, B, Y are clockwise when seen from X then parity is 'e',
                                                          otherwise 'o'
   unknown parity = 'u'
   Example (see XACBY above): if X on the diagram is above the plane ABY
                                                      then parity is 'o'

   Limitations
   ===========
   o Atoms A and B in allenes MUST be connected by a chain of double bonds;

   ==============================================
   Note. Correspondence to CML 0D stereo parities
   ==============================================
   a list of 4 atoms corresponds to CML atomRefs4

   tetrahedral atom
   ================
       CML atomParity > 0 <=> INCHI_PARITY_EVEN
       CML atomParity < 0 <=> INCHI_PARITY_ODD

                                    | 1   1   1   1  |  where xW is x-coordinate of
                                    | xW  xX  xY  xZ |  atom W, etc. (xyz is a
       CML atomParity = determinant | yW  yX  yY  yZ |  'right-handed' Cartesian
                                    | zW  zX  xY  zZ |  coordinate system)

   allene (not yet defined in CML)
   ===============================
   the parity corresponds to the sign of the following determinant
   in exactly same way as for tetrahedral atoms:

       | 1   1   1   1  |  where bonds and neighbor[4] array are
       | xX  xA  xB  xY |  same as defined above for allenes
       | yX  yA  yB  yY |  Obviously, the parity is same for
       | zX  zA  xB  zY |  {#X,#A,#B,#Y} and {#Y,#B,#A,#X}
                           because of the even number of column permutations.

   stereogenic double bond and (not yet defined in CML) cumulenes
   ==============================================================
       CML 'C' (cis)      <=> INCHI_PARITY_ODD
       CML 'T' (trans)    <=> INCHI_PARITY_EVEN


   How InChI uses 0D parities
   ==========================

   1. 0D parities are used if all atom coordinates are zeroes.

   In addition to that:

   2. 0D parities are used for Stereobonds, Allenes, or Cumulenes if:

   2a. A bond to the end-atom is shorter than MIN_BOND_LEN=0.000001
   2b. A ratio of two bond lengths to the end-atom is smaller than MIN_SINE=0.03
   2c. In case of a linear fragment X-A=B end-atom A is treated as satisfying 2a-b

       0D parities are used if 2a or 2b or 2c applies to one or both end-atoms.

   3. 0D parities are used for Tetrahedral Atoms if at least one of 3a-c is true:

   3a. One of bonds to the central atom is shorter than MIN_BOND_LEN=0.000001
   3b. A ratio of two bond lengths to the central atom is smaller than MIN_SINE=0.03
   3c. The four neighbors are almost in one plane or the central atom and
       its only 3 explicit neighbors are almost in one plane

   Notes on 0D parities and 'undefined' stereogenic elements
   =========================================================

   If 0D parity is to be used according to 1-3 but    CH3     CH3
   has not been provided then the corresponding         \    /
   stereogenic element is considered 'undefined'.        C=CH
                                                        /
   For example, if in the structure (Fig. 3)           H
   the explicit H has been moved so that it                Fig. 3
   has same coordinates as atom >C= (that is,
   the length of the bond H-C became zero)
   then the double bond is assigned 'undefined'       CH3      CH3
   parity which by default is omitted from the          \     /
   Identifier.                                           CH=CH

   However, the structure on Fig. 4 will have double        Fig. 4
   bond parity 'o' and its parity in the Identifier is (-).

   Notes on 0D parities in structures containing metals
   ====================================================
   Since InChI disconnects bonds to metals the 0D parities upon the
   disconnection may change in several different ways:

   1) previously non-stereogenic bond may become stereogenic:

         \     /                            \     /  
          CH==CH          disconnection      CH==CH  
           \ /               ======>                 
            M                                  M     

     before the disconnection:    after the disconnection:
     atoms C have valence=5 and   the double bond may become
     the double bond is not       stereogenic
     recognized as stereogenic

   2) previously stereogenic bond may become non-stereogenic:

       M                           M(+)       
        \    /                             / 
         N==C      disconnection    (-)N==C  
             \        ======>              \ 

   3) Oddball structures, usually resulting from projecting 3D
      structures on the plane, may contain fragment like that
      depicted on Fig. 5:

              M   A                      M   A   
              |\ /   B                      /   B 
              | X   /     disconnection    /   /  
              |/ \ /         ======>      /   /   
              C===C                      C===C    
             Fig. 5
     (X stands for bond intersection)
    
     A-C=C-B parity is              A-C=C-B parity is
     trans (e)                      cis (o) or undefined
     because the bond               because C valence = 3,
     orientation is same            not 4.
     as on Fig, 6 below:

          A       M
           \     /     Removal of M from the structure
            C===C      on Fig. 5 changes the geometry from trans
           /     \     to cis. 
          M'      B    Removal of M and M' from the structure
          Fig. 6       on Fig. 6 does not change the A-C=C-B
                       geometry: it is trans.

   To resolve the problem InChI API accepts the second parity
   corresponding to the metal-disconnected structure.
   To store both bond parities use left shift by 3 bits:

   inchi_Stereo0D::parity = ParityOfConnected | (ParityOfDisconnected<<3)

   In case when only disconnected structure parity exists set
   ParityOfConnected = INCHI_PARITY_UNDEFINED.
   This is the only case when INCHI_PARITY_UNDEFINED parity
   may be fed to the InChI.

   In cases when the bond parity in a disconnected structure exists and
   differs from the parity in the connected structure the atoms A and B
   should be non-metals.

****************************************************************************/

#define NO_ATOM          (-1) /* non-existent (central) atom */

/* 0D parity types */
typedef enum tagINCHIStereoType0D {
   INCHI_StereoType_None        = 0,
   INCHI_StereoType_DoubleBond  = 1,
   INCHI_StereoType_Tetrahedral = 2,
   INCHI_StereoType_Allene      = 3
} inchi_StereoType0D;

/* 0D parities */
typedef enum tagINCHIStereoParity0D {
   INCHI_PARITY_NONE      = 0,
   INCHI_PARITY_ODD       = 1,  /* 'o' */
   INCHI_PARITY_EVEN      = 2,  /* 'e' */
   INCHI_PARITY_UNKNOWN   = 3,  /* 'u' */ /* (see also readinch.c)
                                          used in: Extract0DParities, INChITo_Atom  */
   INCHI_PARITY_UNDEFINED = 4   /* '?' -- should not be used; however, see Note above */
} inchi_StereoParity0D;


/*************************************************
 *
 *
 *  0D - S T E R E O  (if no coordinates given)
 *
 *
 *************************************************/


typedef struct tagINCHIStereo0D {
    AT_NUM  neighbor[4];    /* 4 atoms always */
    AT_NUM  central_atom;   /* central tetrahedral atom or a central */
                            /* atom of allene; otherwise NO_ATOM */
    S_CHAR  type;           /* inchi_StereoType0D */
    S_CHAR  parity;         /* inchi_StereoParity0D: may be a combination of two parities: */
                            /* ParityOfConnected | (ParityOfDisconnected << 3), see Note above */
}inchi_Stereo0D;




/*************************************************
 *
 *
 *  I N C h I    D L L     I n p u t
 *
 *
 *************************************************/


/* Structure -> InChI, GetINCHI() / GetStdINCHI() */
typedef struct tagINCHI_Input {
    /* the caller is responsible for the data allocation and deallocation */
    inchi_Atom     *atom;         /* array of num_atoms elements */
    inchi_Stereo0D *stereo0D;     /* array of num_stereo0D 0D stereo elements or NULL */
    char           *szOptions;    /* InChI options: space-delimited; each is preceded by */
                                  /* '/' or '-' depending on OS and compiler */
    AT_NUM          num_atoms;    /* number of atoms in the structure < 1024 */
    AT_NUM          num_stereo0D; /* number of 0D stereo elements */
}inchi_Input;

/* InChI -> Structure, GetStructFromINCHI()/GetStructFromStdINCHI() */
typedef struct tagINCHI_InputINCHI {
    /* the caller is responsible for the data allocation and deallocation */
    char *szInChI;     /* InChI ASCIIZ string to be converted to a strucure */
    char *szOptions;   /* InChI options: space-delimited; each is preceded by */
                       /* '/' or '-' depending on OS and compiler */
} inchi_InputINCHI;


/*************************************************
 *
 *
 *  I N C h I     D L L     O u t p u t
 *
 *
 *************************************************/

/* Structure -> InChI */
typedef struct tagINCHI_Output {
    /* zero-terminated C-strings allocated by GetStdINCHI() */
    /* to deallocate all of them call FreeStdINCHI() (see below) */
    char *szInChI;     /* InChI ASCIIZ string */
    char *szAuxInfo;   /* Aux info ASCIIZ string */
    char *szMessage;   /* Error/warning ASCIIZ message */
    char *szLog;       /* log-file ASCIIZ string, contains a human-readable list */
                       /* of recognized options and possibly an Error/warning message */
} inchi_Output;

/* InChI -> Structure */
typedef struct tagINCHI_OutputStruct {
    /* 4 pointers are allocated by GetStructFromINCHI()/GetStructFromStdINCHI()     */
    /* to deallocate all of them call FreeStructFromStdINCHI()/FreeStructFromStdINCHI() */
    inchi_Atom     *atom;         /* array of num_atoms elements */
    inchi_Stereo0D *stereo0D;     /* array of num_stereo0D 0D stereo elements or NULL */
    AT_NUM          num_atoms;    /* number of atoms in the structure < 1024 */
    AT_NUM          num_stereo0D; /* number of 0D stereo elements */
    char           *szMessage;    /* Error/warning ASCIIZ message */
    char           *szLog;        /* log-file ASCIIZ string, contains a human-readable list */
                                  /* of recognized options and possibly an Error/warning message */
    unsigned long  WarningFlags[2][2]; /* warnings, see INCHIDIFF in inchicmp.h */
                                       /* [x][y]: x=0 => Reconnected if present in InChI otherwise Disconnected/Normal
                                                  x=1 => Disconnected layer if Reconnected layer is present
                                                  y=1 => Main layer or Mobile-H
                                                  y=0 => Fixed-H layer
                                        */
}inchi_OutputStruct;




/*************************************************
 *
 *
 *  I N C h I      D L L     I n t e r f a c e
 *
 *
 *************************************************/

#if (defined( _WIN32 ) && defined( _MSC_VER ) && defined(BUILD_LINK_AS_DLL) )
    /* Win32 & MS VC ++, compile and link as a DLL */
    #ifdef _USRDLL
        /* InChI library dll */
        #define INCHI_API __declspec(dllexport)
        #define EXPIMP_TEMPLATE
        #define INCHI_DECL __stdcall
     #else
        /* calling the InChI dll program */
        #define INCHI_API __declspec(dllimport)
        #define EXPIMP_TEMPLATE extern
        #define INCHI_DECL __stdcall
     #endif
#else
    /* create a statically linked InChI library or link to an executable */
    #define INCHI_API
    #define EXPIMP_TEMPLATE
    #define INCHI_DECL
#endif

/*^^^ Return codes for 
        GetINCHI 
        GetStdINCHI 
        Get_inchi_Input_FromAuxInfo 
        Get_std_inchi_Input_FromAuxInfo 
        GetStructFromINCHI
        GetStructFromStdINCHI
*/
typedef enum tagRetValGetINCHI {
    inchi_Ret_SKIP    = -2, /* not used in InChI library */
    inchi_Ret_EOF     = -1, /* no structural data has been provided */
    inchi_Ret_OKAY    =  0, /* Success; no errors or warnings */
    inchi_Ret_WARNING =  1, /* Success; warning(s) issued */
    inchi_Ret_ERROR   =  2, /* Error: no InChI has been created */
    inchi_Ret_FATAL   =  3, /* Severe error: no InChI has been created (typically, memory allocation failure) */
    inchi_Ret_UNKNOWN =  4, /* Unknown program error */
    inchi_Ret_BUSY    =  5  /* Previuos call to InChI has not returned yet */

} RetValGetINCHI;

/*^^^ Return codes for CheckINCHI */
typedef enum tagRetValCheckINCHI 
{
    INCHI_VALID_STANDARD     =   0,
    INCHI_VALID_NON_STANDARD =  -1,    
    INCHI_INVALID_PREFIX     =   1,
    INCHI_INVALID_VERSION    =   2,
    INCHI_INVALID_LAYOUT     =   3,
    INCHI_FAIL_I2I           =   4
    
} RetValCheckINCHI;



/* to compile all InChI code as a C++ code #define COMPILE_ALL_CPP */
#ifndef COMPILE_ALL_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif


/*^^^ InChI PREFIX */
#define INCHI_STRING_PREFIX "InChI="
#define LEN_INCHI_STRING_PREFIX 6

/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Format:
    
    Standard InChI starts with: InChI=1S/
    Non-standard one with:      InChI=1/
    Empty std InChI:            InChI=1S//
    Empty InChI:                InChI=1//
                                AuxInfo=1//
 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/




/* EXPORTED FUNCTIONS */



/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
GetINCHI / GetStdINCHI


    inchi_Input is created by the user; strings in inchi_Output are allocated and deallocated by InChI
    inchi_Output does not need to be initilized out to zeroes; see FreeNCHI()/FreeSTDINCHI() on how to deallocate it


    Valid options for GetINCHI:
    (use - instead of / for O.S. other than MS Windows)

    Structure perception (compatible with stdInChI)                     
        /NEWPSOFF   /DoNotAddH   /SNon
    Stereo interpretation (lead to generation of non-standard InChI)
        /SRel /SRac /SUCF /ChiralFlagON /ChiralFlagOFF
    InChI creation options (lead to generation of non-standard InChI)
        /SUU /SLUUD   /FixedH  /RecMet  /KET /15T

    
    GetINCHI produces standard InChI if no InChI creation/stereo modification options 
    are specified. Inveresely, if any of SUU/SLUUD/RecMet/FixedH/Ket/15T/SRel/SRac/SUCF 
    options are specified, generated InChI will be non-standard one. 
 
    
    GetStdINCHI produces standard InChI only. 
    The valid structure perception options are:
        /NEWPSOFF   /DoNotAddH   /SNon

    
    Other options are:    
        /AuxNone    Omit auxiliary information (default: Include)
        /Wnumber    Set time-out per structure in seconds; W0 means unlimited
                    In InChI library the default value is unlimited
        /OutputSDF  Output SDfile instead of InChI
        /WarnOnEmptyStructure 
                    Warn and produce empty InChI for empty structure
        /SaveOpt    Save custom InChI creation options (non-standard InChI)

 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
EXPIMP_TEMPLATE INCHI_API int INCHI_DECL GetINCHI( inchi_Input *inp, inchi_Output *out );
EXPIMP_TEMPLATE INCHI_API int INCHI_DECL GetStdINCHI( inchi_Input *inp, inchi_Output *out );



/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
FreeINCHI / FreeStdINCHI

    should be called to deallocate char* pointers 
    obtained from each GetINCHI /GetStdINCHI call                         

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
EXPIMP_TEMPLATE INCHI_API void INCHI_DECL FreeINCHI ( inchi_Output *out );
EXPIMP_TEMPLATE INCHI_API void INCHI_DECL FreeStdINCHI ( inchi_Output *out );



/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
GetStringLength
    
    helper: get string length

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
EXPIMP_TEMPLATE INCHI_API int INCHI_DECL GetStringLength( char *p );


/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
GetStructFromINCHI / GetStructFromStdINCHI

    inchi_Inputinchi_InputINCHI is created by the user; pointers in inchi_OutputStruct are allocated and deallocated by InChI 
    inchi_OutputStruct does not need to be initilized out to zeroes; see FreeStructFromStdINCHI() on how to deallocate it 
    Option /Inchi2Struct is not needed for GetStructFromINCHI()/GetStructFromStdINCHI()

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
EXPIMP_TEMPLATE INCHI_API int INCHI_DECL GetStructFromINCHI( inchi_InputINCHI *inpInChI, inchi_OutputStruct *outStruct );
EXPIMP_TEMPLATE INCHI_API int INCHI_DECL GetStructFromStdINCHI( inchi_InputINCHI *inpInChI, inchi_OutputStruct *outStruct );


/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
FreeStructFromINCHI / FreeStructFromStdINCHI

    should be called to deallocate pointers obtained from each 
    GetStructFromStdINCHI / GetStructFromINCHI

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
EXPIMP_TEMPLATE INCHI_API void INCHI_DECL FreeStructFromINCHI( inchi_OutputStruct *out );
EXPIMP_TEMPLATE INCHI_API void INCHI_DECL FreeStructFromStdINCHI( inchi_OutputStruct *out );


/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
GetINCHIfromINCHI

    GetINCHIfromINCHI does same as -InChI2InChI option: converts InChI into InChI for validation purposes    
    It may also be used to filter out specific layers. For instance, /Snon would remove stereochemical layer 
    Omitting /FixedH and/or /RecMet would remove Fixed-H or Reconnected layers                               
    To keep all InChI layers use options string "/FixedH /RecMet"; option /InChI2InChI is not needed         
    inchi_InputINCHI is created by the user; strings in inchi_Output are allocated and deallocated by InChI  
    inchi_Output does not need to be initilized out to zeroes; see FreeINCHI() on how to deallocate it  
                                                                                                        
    Note: there is no explicit tool to conversion from/to standard InChI                                   

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
EXPIMP_TEMPLATE INCHI_API int INCHI_DECL GetINCHIfromINCHI( inchi_InputINCHI *inpInChI, inchi_Output *out );

#ifndef COMPILE_ALL_CPP
#ifdef __cplusplus
}
#endif
#endif


/*****************************************************************
 *
 *
 *   C o n v e r s i o n:   InChI  AuxInfo string => inchi_Input
 *
 *
 *****************************************************************/

#ifndef STR_ERR_LEN
#define STR_ERR_LEN     256
#endif

typedef struct tagInchiInpData {
    inchi_Input *pInp;    /* a pointer to pInp that has all items 0 or NULL */
    int          bChiral; /* 1 => the structure was marked as chiral, 2=> not chiral, 0=> not marked */
    char         szErrMsg[STR_ERR_LEN];
} InchiInpData;

/* to compile all InChI code as a C++ code #define COMPILE_ALL_CPP */
#ifndef COMPILE_ALL_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif



/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Get_inchi_Input_FromAuxInfo / Get_std_inchi_Input_FromAuxInfo

Input:
    szInchiAuxInfo: contains ASCIIZ string of InChI output for a single
                   structure or only the AuxInfo line
    bDoNotAddH:    if 0 then InChI will be allowed to add implicit H
    bDiffUnkUndfStereo
                   if not 0, use different labels for unknown and undefined stereo
    pInchiInp:     should have a valid pointer pInchiInp->pInp to an empty 
                   (all members = 0) inchi_Input structure 

Output:
    pInchiInp:     The following members of pInp may be filled during the call:
                   atom, num_atoms, stereo0D, num_stereo0D
    Return value:  see RetValGetINCHI

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
EXPIMP_TEMPLATE INCHI_API int INCHI_DECL Get_inchi_Input_FromAuxInfo( 
                                                        char *szInchiAuxInfo, 
                                                        int bDoNotAddH, 
                                                        int bDiffUnkUndfStereo,
                                                        InchiInpData *pInchiInp );
EXPIMP_TEMPLATE INCHI_API int INCHI_DECL Get_std_inchi_Input_FromAuxInfo( char *szInchiAuxInfo, 
                                                        int bDoNotAddH, 
                                                        InchiInpData *pInchiInp );



/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Free_inchi_Input / Free_std_inchi_Input

    To deallocate and write zeroes into the changed members of pInchiInp->pInp call
    Free_std_inchi_Input( inchi_Input *pInp )

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
EXPIMP_TEMPLATE INCHI_API void INCHI_DECL Free_inchi_Input( inchi_Input *pInp );
EXPIMP_TEMPLATE INCHI_API void INCHI_DECL Free_std_inchi_Input( inchi_Input *pInp );



/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
CheckINCHI

Check if the string represents valid InChI/standard InChI.          
Input:
    szINCHI     source InChI
    strict      if 0, just briefly check for proper layout (prefix, version, etc.)
                The result may not be strict.
                If not 0, try to perform InChI2InChI conversion and 
                returns success if a resulting InChI string exactly match source.
                The result may be 'false alarm' due to imperfectness of conversion.
Returns:
    success/errors codes

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
EXPIMP_TEMPLATE INCHI_API int INCHI_DECL CheckINCHI(const char *szINCHI, const int strict);    


#ifndef COMPILE_ALL_CPP
#ifdef __cplusplus
}
#endif
#endif



/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    
                                InChIKey API
                                

    InChIKey description



    The InChIKey is a character signature based on a hash code of the InChI string.        
    Standard InChIKey is produced out of standard InChI.
    Non-standard InChIKey is produced out of non-standard InChI.
           
                    AAAAAAAAAAAAAA-BBBBBBBBCD-P


    InChIKey layout is as follows:

    AAAAAAAAAAAAAA
        First block (14 letters)
        Encodes molecular skeleton (connectivity)

    BBBBBBBB
        Second block (8 letters)
        Encodes tautomers, stereochemistry, isotopomers, reconnected layer 
    C
        'S' for standard
        'N' for non-standard
    D
        InChI version ('A' for 1)
    P - (de)protonation flag
        Protonization encoding:
        N 0
        O +1 P +2 Q +3 R +4 S +5 T +6 U +7 V +8 W +9 X +10 Y +11 Z +12
        M -1 L-2 K -3 J -4 I -5 H -6 G -7 F -8 E -9 D -10 C -11 B -12 
        A < -12 or > +12


    All symbols except delimiter (dash, that is, minus) are uppercase English 
    letters representing a "base 26" encoding.
    
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/


/*^^^ Return codes for key generation procedure */
#define INCHIKEY_OK 0
#define INCHIKEY_UNKNOWN_ERROR 1
#define INCHIKEY_EMPTY_INPUT 2
#define INCHIKEY_INVALID_INCHI_PREFIX 3
#define INCHIKEY_NOT_ENOUGH_MEMORY 4
#define INCHIKEY_INVALID_INCHI 20
#define INCHIKEY_INVALID_STD_INCHI 21


/*^^^ Return codes for CheckINCHIKey */
typedef enum tagRetValGetINCHIKey
{
    INCHIKEY_VALID_STANDARD     =   0,
    INCHIKEY_VALID_NON_STANDARD =  -1,
    INCHIKEY_INVALID_LENGTH     =   1,
    INCHIKEY_INVALID_LAYOUT     =   2,
    INCHIKEY_INVALID_VERSION    =   3
} RetValCheckINCHIKeyv;



/* EXPORTED FUNCTIONS */



/* To compile all InChI code as a C++ code #define COMPILE_ALL_CPP */

#ifndef COMPILE_ALL_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif



/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
GetINCHIKeyFromINCHI

Calculate InChIKey by InChI string. 

Input:
        szINCHISource
            source InChI string 
        xtra1
            =1 calculate hash extension (up to 256 bits; 1st block)
        xtra2
            =1 calculate hash extension (up to 256 bits; 2nd block)

Output:
        szINCHIKey
            InChIKey string 
            The user-supplied buffer szINCHIKey should be at least 28 bytes long.
        szXtra1
            hash extension (up to 256 bits; 1st block) string 
            Caller should allocate space for 64 characters + trailing NULL
        szXtra2
            hash extension (up to 256 bits; 2nd block) string 
            Caller should allocate space for 64 characters + trailing NULL

Returns:
        success/errors codes

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
EXPIMP_TEMPLATE INCHI_API int INCHI_DECL GetINCHIKeyFromINCHI(const char* szINCHISource, 
                                                              const int xtra1,
                                                              const int xtra2,
                                                              char* szINCHIKey, 
                                                              char* szXtra1, 
                                                              char* szXtra2);



/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
GetStdINCHIKeyFromStdINCHI

    "Standard" counterpart

    For compatibility with v. 1.02std, no extra hash calculation is allowed.
    To calculate extra hash(es), use GetINCHIKeyFromINCHI with stdInChI as input.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
EXPIMP_TEMPLATE INCHI_API int INCHI_DECL GetStdINCHIKeyFromStdINCHI(const char* szINCHISource, 
                                                                    char* szINCHIKey);


/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
CheckINCHIKey

Check if the string represents valid InChIKey.          
Input:
        szINCHIKey
            source InChIKey string 
Returns:
        success/errors codes

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
EXPIMP_TEMPLATE INCHI_API int INCHI_DECL CheckINCHIKey(const char *szINCHIKey);    

#ifndef COMPILE_ALL_CPP
#ifdef __cplusplus
}
#endif
#endif



/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    
                          Modularized InChI generation API



    Note. Functions with STDINCHIGEN prefix are 
    retained for compatibility with v. 1.02std

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/




/*^^^ Data structures holding intermediate (normalization) results */

#ifndef MAX_NUM_STEREO_ATOM_NEIGH
#define MAX_NUM_STEREO_ATOM_NEIGH 4
#endif
#ifndef MAX_NUM_STEREO_BONDS
#define MAX_NUM_STEREO_BONDS      3
#endif
#ifndef INCHI_NUM
#define INCHI_NUM            2    /* = array size; member indexes: */
#endif


typedef unsigned short AT_NUMBR;    
typedef signed short NUM_HS; 
typedef unsigned long INCHI_MODES;        

typedef struct tagNormAtom 
{
    char          elname[ATOM_EL_LEN];      /* chem. element name */
    U_CHAR        el_number;                /* number of the element in the Periodic Table */
    AT_NUMBR      neighbor[MAXVAL];         /* positions (from 0) of the neighbors in the NORM_ATOM array */
    AT_NUMBR      orig_at_number;           /* original atom number, starts from 1 */
    AT_NUMBR      orig_compt_at_numb;       /* atom number within a component before terminal H removal */
    S_CHAR        bond_stereo[MAXVAL];      /* 1=Up,4=Either,6=Down (this atom is at the pointing wedge)
                                               negative => on the opposite side of the wedge; 3=Either double bond  */
    U_CHAR        bond_type[MAXVAL];        /* 1=single, 2=double, 3=triple, 4=1/2 (bond order is 1 or 2) */
                                            /* 5=1/2/3, 6=1/3, 7=2/3, 8=tautomeric, 9=1/2 non-stereogenic */

    S_CHAR        valence;                  /* number of bonds = number of neighbors not greater than MAXVAL */
    S_CHAR        chem_bonds_valence;       /* sum of bond types (1,2,3); type 4 needs special treatment */
    S_CHAR        num_H;                    /* number of adjacent implicit hydrogen atoms including D and T    */
    S_CHAR        num_iso_H[NUM_H_ISOTOPES];/* number of adjacent implicit 1H(protium), 2H(D), 3H(T) < 16 */
    S_CHAR        iso_atw_diff;             /* =0 => natural isotopic abundances  */
                                            /* >0 => (isotopic mass) - (rounded average atomic mass) + 1 */
                                            /* <0 => (isotopic mass) - (rounded average atomic mass) */
    S_CHAR        charge;                   /* charge */
    S_CHAR        radical;                  /* RADICAL_SINGLET, RADICAL_DOUBLET, or RADICAL_TRIPLET */
    S_CHAR        bAmbiguousStereo;         /* flag of detected stereo ambiguity */
    S_CHAR        cFlags;                   /* AT_FLAG_ISO_H_POINT: atom may have exchangeable isotopic H */
    AT_NUMBR      at_type;                  /* ATT_NONE, ATT_ACIDIC, etc. See InChI normalization code */
    AT_NUMBR      component;                /* number of the structure component > 0 */
    AT_NUMBR      endpoint;                 /* id of a tautomeric group */
    AT_NUMBR      c_point;                  /* id of a positive charge group */
    double        x;                        /* x coordinate */
    double        y;                        /* y coordinate */
    double        z;                        /* x coordinate */
    /*---------  0D parities ----------*/
    S_CHAR        bUsed0DParity;            /* bit=1 => stereobond; bit=2 => stereocenter */
    /*-----  tetrahedral stereo parity */
    S_CHAR        p_parity;                 /* tetrahedral (sp3) cml parity */
    AT_NUMBR      p_orig_at_num[MAX_NUM_STEREO_ATOM_NEIGH]; /* orig_at_number of each neighbor > 0; 0=> no neighbor */
    /*----- stereo bond (SB) parities */
    S_CHAR        sb_ord[MAX_NUM_STEREO_BONDS];  /* neighbor[] index of another end of this SB, starts from 0 */
    S_CHAR        sn_ord[MAX_NUM_STEREO_BONDS];  /* neighbor[] index of a bond that is not this SB; starts from 0;
                                                  -1 means the neighbor is a removed explicit H */
    /* atoms on both ends of a stereobond have same parity => trans/T/E/2, diff. parities => cis/C/Z/1 */
    S_CHAR        sb_parity[MAX_NUM_STEREO_BONDS];       /* parities of stereobonds (sp2) incident to this atom */	
    AT_NUMBR      sn_orig_at_num[MAX_NUM_STEREO_BONDS];  /* orig_at_number of sn_ord[] neighbor > 0 */

#if ( FIND_RING_SYSTEMS == 1 )
    S_CHAR  bCutVertex;                    /* is the atom a cut-vertex or not */
    AT_NUMBR nRingSystem;                  /* starts from 1; number of a ring system */
    AT_NUMBR nNumAtInRingSystem;           /* number of atoms in a ring system to which this at belongs */
    AT_NUMBR nBlockSystem;                 /* ambiguous if the atom is a cut-vertex: better apply this to bonds */

#if ( FIND_RINS_SYSTEMS_DISTANCES == 1 )
    AT_NUMBR nDistanceFromTerminal;        /* not used */
#endif

#endif
} NORM_ATOM;



typedef struct tagNormAtomData
{
    NORM_ATOM *at;                  /* atom list */
    NORM_ATOM *at_fixed_bonds;      /* atom list with added or removed protons only */
    int       num_at;               /* number of atoms except removed terminal H */
    int       num_removed_H;        /* number of removed H; at[] has (num_at+num_removed_H) elements */
    int       num_bonds;
    int       num_isotopic;         /* number of isotopic atoms */
    int       bExists;              /* for internal use */
    int       bDeleted;             /* for internal use */
    int       bHasIsotopicLayer;
    int       bTautomeric;
    int       bTautPreprocessed;    /* for internal use */
    int       nNumRemovedProtons;
    NUM_HS    nNumRemovedProtonsIsotopic[NUM_H_ISOTOPES]; 
                                    /* isotopic composition of removed protons, not included in num_iso_H[] */
    NUM_HS    num_iso_H[NUM_H_ISOTOPES]; 
                                    /* isotopic H on tautomeric atoms and those 
                                       in nIsotopicEndpointAtomNumber */
    INCHI_MODES bTautFlags;         /* for internal use */
    INCHI_MODES bTautFlagsDone;     /* for internal use */
    INCHI_MODES bNormalizationFlags;/* for internal use */
} NORM_ATOMS;


typedef struct tagINCHIGEN_DATA 
{
    
    char          pStrErrStruct[STR_ERR_LEN]; /* intermediate log (warning/error report) */
    int           num_components[INCHI_NUM];  /* number of allocated INChI, INChI_Aux data structures */
                                              /* index=0 => disconnected, 1 => reconnected structure */
    
    /*^^^ The results of normalization stage */
    /*^^^ for each member of pair disconnected/reconnected structures: */
    NORM_ATOMS   *NormAtomsNontaut[INCHI_NUM]; 
    NORM_ATOMS   *NormAtomsTaut[INCHI_NUM];

} INCHIGEN_DATA;


/*^^^ InChI Generator Handle */

typedef void* INCHIGEN_HANDLE;




/* EXPORTED FUNCTIONS */



/* to compile all InChI code as a C++ code #define COMPILE_ALL_CPP */
#ifndef COMPILE_ALL_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif



/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
INCHIGEN_Create / STDINCHIGEN_Create

    InChI Generator: create generator
    Returns handle of generator object or NULL on failure

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
EXPIMP_TEMPLATE INCHI_API 
INCHIGEN_HANDLE INCHI_DECL INCHIGEN_Create(void);
EXPIMP_TEMPLATE INCHI_API 
INCHIGEN_HANDLE INCHI_DECL STDINCHIGEN_Create(void);



/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
INCHIGEN_Setup / STDINCHIGEN_Setup

    InChI Generator: setup

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
EXPIMP_TEMPLATE INCHI_API int INCHI_DECL INCHIGEN_Setup(INCHIGEN_HANDLE HGen, 
                                                        INCHIGEN_DATA * pGenData, 
                                                        inchi_Input * pInp);
EXPIMP_TEMPLATE INCHI_API int INCHI_DECL STDINCHIGEN_Setup(INCHIGEN_HANDLE HGen, 
                                                           INCHIGEN_DATA * pGenData, 
                                                           inchi_Input * pInp);



/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
INCHIGEN_DoNormalization / STDINCHIGEN_DoNormalization

    InChI Generator: structure normalization stage

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
EXPIMP_TEMPLATE INCHI_API int INCHI_DECL INCHIGEN_DoNormalization(INCHIGEN_HANDLE HGen, 
                                                                     INCHIGEN_DATA * pGenData);
EXPIMP_TEMPLATE INCHI_API int INCHI_DECL STDINCHIGEN_DoNormalization(INCHIGEN_HANDLE HGen, 
                                                                     INCHIGEN_DATA * pGenData);



/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
INCHIGEN_DoCanonicalization / STDINCHIGEN_DoCanonicalization

    InChI Generator: structure canonicalization stage

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
EXPIMP_TEMPLATE INCHI_API int INCHI_DECL INCHIGEN_DoCanonicalization
                                (INCHIGEN_HANDLE HGen, INCHIGEN_DATA * pGenData);
EXPIMP_TEMPLATE INCHI_API int INCHI_DECL STDINCHIGEN_DoCanonicalization
                                (INCHIGEN_HANDLE HGen, INCHIGEN_DATA * pGenData);



/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
INCHIGEN_DoSerialization / STDINCHIGEN_DoSerialization

    InChI Generator: InChI serialization stage

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
EXPIMP_TEMPLATE INCHI_API int INCHI_DECL INCHIGEN_DoSerialization(INCHIGEN_HANDLE HGen, 
                                                                  INCHIGEN_DATA * pGenData, 
                                                                  inchi_Output * pResults);
EXPIMP_TEMPLATE INCHI_API int INCHI_DECL STDINCHIGEN_DoSerialization(INCHIGEN_HANDLE HGen, 
                                                                     INCHIGEN_DATA * pGenData, 
                                                                     inchi_Output * pResults);



/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
INCHIGEN_DoSerialization / STDINCHIGEN_DoSerialization

    InChI Generator: reset stage (use before get next structure)

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
EXPIMP_TEMPLATE INCHI_API 
void INCHI_DECL INCHIGEN_Reset(INCHIGEN_HANDLE HGen, 
                               INCHIGEN_DATA * pGenData, 
                               inchi_Output * pResults);
EXPIMP_TEMPLATE INCHI_API 
void INCHI_DECL STDINCHIGEN_Reset(INCHIGEN_HANDLE HGen, 
                               INCHIGEN_DATA * pGenData, 
                               inchi_Output * pResults);



/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
INCHIGEN_DoSerialization / STDINCHIGEN_DoSerialization

    InChI Generator: destroy generator

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
EXPIMP_TEMPLATE INCHI_API void INCHI_DECL INCHIGEN_Destroy(INCHIGEN_HANDLE HGen);
EXPIMP_TEMPLATE INCHI_API void INCHI_DECL STDINCHIGEN_Destroy(INCHIGEN_HANDLE HGen);





#ifndef COMPILE_ALL_CPP
#ifdef __cplusplus
}
#endif
#endif



/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Prototypes for C calling conventions:

    int  GetINCHI( inchi_Input *inp, inchi_Output *out );
    int  GetStdINCHI( inchi_Input *inp, inchi_Output *out );
    void FreeINCHI( inchi_Output *out );
    void FreeStdINCHI( inchi_Output *out );
    int  GetStringLength( char *p );
    int  Get_inchi_Input_FromAuxInfo
     ( char *szInchiAuxInfo, int bDoNotAddH, int bDiffUnkUndfStereo, InchiInpData *pInchiInp );
    int  Get_std_inchi_Input_FromAuxInfo
     ( char *szInchiAuxInfo, int bDoNotAddH, int bDiffUnkUndfStereo,InchiInpData *pInchiInp );
    void Free_inchi_Input( inchi_Input *pInp );
    void Free_std_inchi_Input( inchi_Input *pInp );
    int GetStructFromINCHI( inchi_InputINCHI *inpInChI, inchi_OutputStruct *outStruct );
    int GetStructFromStdINCHI( inchi_InputINCHI *inpInChI, inchi_OutputStruct *outStruct );
    void FreeStructFromStdINCHI( inchi_OutputStruct *out );

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/


/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Win32 Dumpbin export information

    ordinal hint RVA      name
cdecl 
          1    0 000B7EAB CheckINCHI
          2    1 000B7221 CheckINCHIKey
          3    2 000B7C62 FreeINCHI
          4    3 000B7B04 FreeStdINCHI
          5    4 000B72B7 FreeStructFromINCHI
          6    5 000B7BC2 FreeStructFromStdINCHI
          7    6 000B7E33 Free_inchi_Input
          8    7 000B7C58 Free_std_inchi_Input
          9    8 000B727B GetINCHI
         10    9 000B75B4 GetINCHIKeyFromINCHI
         11    A 000B757D GetINCHIfromINCHI
         12    B 000B8211 GetStdINCHI
         13    C 000B7F0A GetStdINCHIKeyFromStdINCHI
         14    D 000B77CB GetStringLength
         15    E 000B7CA3 GetStructFromINCHI
         16    F 000B778A GetStructFromStdINCHI
         17   10 000B7DAC Get_inchi_Input_FromAuxInfo
         18   11 000B7D6B Get_std_inchi_Input_FromAuxInfo
         19   12 000B7D6B INCHIGEN_Create
         20   13 000B7F2D INCHIGEN_Destroy
         21   14 000B7F23 INCHIGEN_DoCanonicalization
         22   15 000B7F23 INCHIGEN_DoNormalization
         23   16 000B714A INCHIGEN_DoSerialization
         24   17 000B7FCD INCHIGEN_Reset
         25   18 000B7FCD INCHIGEN_Setup
         26   19 000B7EA6 STDINCHIGEN_Create
         27   1A 000B7EA6 STDINCHIGEN_Destroy
         28   1B 000B711D STDINCHIGEN_DoCanonicalization
         29   1C 000B7073 STDINCHIGEN_DoNormalization
         30   1D 000B7FC3 STDINCHIGEN_DoSerialization
         31   1E 000B7668 STDINCHIGEN_Reset
         32   1F 000B7438 STDINCHIGEN_Setup
__stdcall or PASCAL
         33   20 000B7DFC _CheckINCHI@8
         34   21 000B7802 _CheckINCHIKey@4
         35   22 000B7F73 _FreeINCHI@4
         36   23 000B7F82 _FreeStdINCHI@4
         37   24 000B75E1 _FreeStructFromINCHI@4
         38   25 000B7B81 _FreeStructFromStdINCHI@4
         39   26 000B7B86 _Free_inchi_Input@4
         40   27 000B7A96 _Free_std_inchi_Input@4
         41   28 000B7B5E _GetINCHI@8
         42   29 000B7285 _GetINCHIKeyFromINCHI@24
         43   2A 000B758C _GetINCHIfromINCHI@8
         44   2B 000B7CDA _GetStdINCHI@8
         45   2C 000B7979 _GetStdINCHIKeyFromStdINCHI@8
         46   2D 000B7BA4 _GetStringLength@4
         47   2E 000B70A5 _GetStructFromINCHI@8
         48   2F 000B79B0 _GetStructFromStdINCHI@8
         49   30 000B8022 _Get_inchi_Input_FromAuxInfo@16
         50   31 000B76E0 _Get_std_inchi_Input_FromAuxInfo@12
         51   32 000B7230 _INCHIGEN_Create@0
         52   33 000B760E _INCHIGEN_Destroy@4
         53   34 000B7087 _INCHIGEN_DoCanonicalization@8
         54   35 000B70B4 _INCHIGEN_DoNormalization@8
         55   36 000B72D5 _INCHIGEN_DoSerialization@12
         56   37 000B7FE1 _INCHIGEN_Reset@12
         57   38 000B7163 _INCHIGEN_Setup@12
         58   39 000B7159 _STDINCHIGEN_Create@0
         59   3A 000B78A7 _STDINCHIGEN_Destroy@4
         60   3B 000B72F3 _STDINCHIGEN_DoCanonicalization@8
         61   3C 000B737A _STDINCHIGEN_DoNormalization@8
         62   3D 000B7B72 _STDINCHIGEN_DoSerialization@12
         63   3E 000B7654 _STDINCHIGEN_Reset@12
         64   3F 000B75FF _STDINCHIGEN_Setup@12


    Note. Currently there is no callback function for aborting, progress, etc.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/

#endif /* __INHCH_API_H__ */

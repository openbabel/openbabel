/*
 * International Chemical Identifier (InChI)
 * Version 1
 * Software version 1.02-beta
 * August 23, 2007
 * Developed at NIST
 *
 * The InChI library and programs are free software developed under the
 * auspices of the International Union of Pure and Applied Chemistry (IUPAC);
 * you can redistribute this software and/or modify it under the terms of 
 * the GNU Lesser General Public License as published by the Free Software 
 * Foundation:
 * http://www.opensource.org/licenses/lgpl-license.php
 */


/*
 The code in this #include file reads InChI AuxInfo
*/

/****************************************************************************/
#define MIN_BOND_LENGTH   (1.0e-6)
#define INCHI_LINE_LEN   512 /*1024*/ /*256*/ 
#define INCHI_LINE_ADD   384  /*128*/  /*64*/
/* Note: (INCHI_LINE_LEN - INCHI_LINE_ADD) > (length of the longest item: szCoord) = 33 */
/*****************************************************************************/

#if( defined(INCHI_LIBRARY) || defined(INCHI_MAIN) )

#define AB_MAX_WELL_DEFINED_PARITY inchi_max(INCHI_PARITY_ODD, INCHI_PARITY_EVEN) /* 1, 2 => well defined parities, uncluding 'unknown' */
#define AB_MIN_WELL_DEFINED_PARITY inchi_min(INCHI_PARITY_ODD, INCHI_PARITY_EVEN) /* min(INCHI_PARITY_ODD, INCHI_PARITY_EVEN) */
#define ATOM_PARITY_WELL_DEF(X)     (AB_MIN_WELL_DEFINED_PARITY <= (X) && (X) <= AB_MAX_WELL_DEFINED_PARITY)

#define inchi_NUMH2(AT,CUR_AT) ((AT[CUR_AT].num_iso_H[0]>0?AT[CUR_AT].num_iso_H[0]:0) +AT[CUR_AT].num_iso_H[1]+AT[CUR_AT].num_iso_H[2]+AT[CUR_AT].num_iso_H[3])

#define SB_PARITY_FLAG  0x38 /* disconnected structure has undef. parity */
#define SB_PARITY_SHFT  3
#define SB_PARITY_MASK  0x07
#define SB_PARITY_1(X) (X & SB_PARITY_MASK)  /* refers to connected structure */
#define SB_PARITY_2(X) (((X) >> SB_PARITY_SHFT) & SB_PARITY_MASK) /* refers to connected structure */


int      str_fgetc( INCHI_FILE *f );
char    *str_fgets( char *szLine, int len, INCHI_FILE *f );
int      my_fgets( char *szLine, int len, INCHI_FILE *f, int *bTooLongLine );

#endif


#ifdef INCHI_LIBRARY

void            FreeInchi_Atom( inchi_Atom **at );
inchi_Atom     *CreateInchi_Atom( int num_atoms );
void            FreeInchi_Input( inchi_Input *inp_at_data );
S_SHORT        *is_in_the_slist( S_SHORT *pathAtom, S_SHORT nNextAtom, int nPathLen );
int             is_element_a_metal( char szEl[] );

char    *str_fgetsTab( char *szLine, int len, INCHI_FILE *f );
int      my_fgetsTab( char *szLine, int len, INCHI_FILE *f, int *bTooLongLine );
int      my_fgetsTab1( char *szLine, int len, INCHI_FILE *f, int *bTooLongLine );

#endif

#ifndef INCHI_MAIN

void            FreeInchi_Stereo0D( inchi_Stereo0D **stereo0D );
inchi_Stereo0D *CreateInchi_Stereo0D( int num_stereo0D );
int Extract0DParities( inp_ATOM *at, int nNumAtoms, inchi_Stereo0D *stereo0D,
                       int num_stereo0D, char *pStrErr, int *err );

#endif


#if( defined(INCHI_LIBRARY) || defined(INCHI_MAIN) )
/*******************************************************************/
int str_fgetc( INCHI_FILE *f )
{
#ifdef INCHI_LIBRARY
    if ( f->nPtr < f->nUsedLength ) {
        return (int)f->pStr[f->nPtr++];
    }
    return EOF;
#else
    return fgetc( f );
#endif
}
/*******************************************************************/
char *str_fgets( char *szLine, int len, INCHI_FILE *f )
{
    int  length=0, c=0;
    if ( -- len < 0 ) {
        return NULL;
    }
    while ( length < len && EOF != (c = str_fgetc( f )) ) {
        szLine[length++] = (char)c;
        if ( c == '\n' )
            break;
    }
    if ( !length && EOF == c ) {
        return NULL;
    }
    szLine[length] = '\0';
    return szLine;
}
/*******************************************************************/
int my_fgets( char *szLine, int len, INCHI_FILE *f, int *bTooLongLine )
{
    int  length;
    char *p;
    do {
        p = str_fgets( szLine, len-1, f );
        if ( !p ) {
            *bTooLongLine = 0;
            return -1; /* end of file or cannot read */
        }
        szLine[len-1] = '\0';
        /*
        *bTooLongLine = !strchr( szLine, '\n' );
        */
        p = strchr( szLine, '\n' );
        *bTooLongLine = ( !p && ((int)strlen(szLine)) == len-2 );
        LtrimRtrim( szLine, &length );
    } while ( !length );
    return length;
}

#endif


#ifdef INCHI_LIBRARY
/******************************************************************************************************/
void FreeInchi_Atom( inchi_Atom **at )
{
    if ( at && *at ) {
        inchi_free( *at );
        *at = NULL;
    }
}
/******************************************************************************************************/
inchi_Atom *CreateInchi_Atom( int num_atoms )
{
   inchi_Atom *p = (inchi_Atom* ) inchi_calloc(num_atoms, sizeof(inchi_Atom) );
   return p;
}
/******************************************************************************************************/
void FreeInchi_Input( inchi_Input *inp_at_data )
{
    FreeInchi_Atom( &inp_at_data->atom );
    FreeInchi_Stereo0D( &inp_at_data->stereo0D );
    memset( inp_at_data, 0, sizeof(*inp_at_data) );
}
/*************************************************************************/
S_SHORT *is_in_the_slist( S_SHORT *pathAtom, S_SHORT nNextAtom, int nPathLen )
{
    for ( ; nPathLen && *pathAtom != nNextAtom; nPathLen--,  pathAtom++ )
        ;
    return nPathLen? pathAtom : NULL;
}
/************************************************/
int is_element_a_metal( char szEl[] )
{
    static const char szMetals[] = "K;V;Y;W;U;"
        "Li;Be;Na;Mg;Al;Ca;Sc;Ti;Cr;Mn;Fe;Co;Ni;Cu;Zn;Ga;Rb;Sr;Zr;"
        "Nb;Mo;Tc;Ru;Rh;Pd;Ag;Cd;In;Sn;Sb;Cs;Ba;La;Ce;Pr;Nd;Pm;Sm;"
        "Eu;Gd;Tb;Dy;Ho;Er;Tm;Yb;Lu;Hf;Ta;Re;Os;Ir;Pt;Au;Hg;Tl;Pb;"
        "Bi;Po;Fr;Ra;Ac;Th;Pa;Np;Pu;Am;Cm;Bk;Cf;Es;Fm;Md;No;Lr;Rf;";
    const int len = strlen(szEl);
    const char *p;

    if ( 0 < len && len <= 2 &&
         isalpha( UCINT szEl[0] ) && isupper( szEl[0] ) &&
         (p = strstr(szMetals, szEl) ) && p[len] == ';' ) {

            return 1; /*return AtType_Metal;*/
    }
    return 0;
}
/*******************************************************************/
/* read up to len or tab or LF; if empty read next until finds non-empty line   */
/* remove leading and trailing white spaces; keep zero termination */
/*******************************************************************/
char *str_fgetsTab( char *szLine, int len, INCHI_FILE *f )
{
    int  length=0, c=0;
    if ( --len < 0 ) {
        return NULL;
    }
    while ( length < len && EOF != (c = str_fgetc( f )) ) {
        if ( c == '\t' )
            c = '\n';
        szLine[length++] = (char)c;
        if ( c == '\n' )
            break;
    }
    if ( !length && EOF == c ) {
        return NULL;
    }
    szLine[length] = '\0';
    return szLine;
}
/*******************************************************************/
/* read up to len or tab or LF; if empty read next until finds non-empty line   */
/* remove leading and trailing white spaces; keep zero termination */
/*******************************************************************/
int my_fgetsTab( char *szLine, int len, INCHI_FILE *f, int *bTooLongLine )
{
    int  length;
    char *p;
    do {
        p = str_fgetsTab( szLine, len-1, f );
        if ( !p ) {
            *bTooLongLine = 0;
            return -1; /* end of file or cannot read */
        }
        szLine[len-1] = '\0';
        /*
        *bTooLongLine = !strchr( szLine, '\n' );
        */
        p = strchr( szLine, '\n' );
        *bTooLongLine = ( !p && ((int)strlen(szLine)) == len-2 );
        LtrimRtrim( szLine, &length );
    } while ( !length );
    return length;
}
/*******************************************************************/
int my_fgetsTab1( char *szLine, int len, INCHI_FILE *f, int *bTooLongLine )
{
    int  length;
    char *p;
    /*do {*/
        p = str_fgetsTab( szLine, len-1, f );
        if ( !p ) {
            *bTooLongLine = 0;
            return -1; /* end of file or cannot read */
        }
        szLine[len-1] = '\0';
        /*
        *bTooLongLine = !strchr( szLine, '\n' );
        */
        p = strchr( szLine, '\n' );
        *bTooLongLine = ( !p && ((int)strlen(szLine)) == len-2 );
        LtrimRtrim( szLine, &length );
    /*} while ( !length );*/
    return length;
}

#endif


#ifndef INCHI_MAIN
/******************************************************************************************************/
inchi_Stereo0D *CreateInchi_Stereo0D( int num_stereo0D )
{
   return (inchi_Stereo0D* ) inchi_calloc(num_stereo0D, sizeof(inchi_Stereo0D) );
}
/******************************************************************************************************/
void FreeInchi_Stereo0D( inchi_Stereo0D **stereo0D )
{
    if ( stereo0D && *stereo0D ) {
        inchi_free( *stereo0D );
        *stereo0D = NULL;
    }
}
#endif


#if( defined(INCHI_LIBRARY) || defined(INCHI_MAIN) )

#ifdef INCHI_LIBRARY
#define INChITo_Atom        ll_INChIToInchi_Atom
#else
#define INChITo_Atom        ee_INChIToIbChI_Atom
#define FindToken           e_FindToken
#define LoadLine            e_LoadLine
#endif

#define AT_NUM_BONDS(AT)    (AT).num_bonds
#define ATOM_NUMBER         AT_NUM
#define IN_NEIGH_LIST       is_in_the_slist
#define INPUT_FILE          INCHI_FILE
#define Create_Atom         CreateInchi_Atom 
#define AT_BONDS_VAL(AT,I)  AT[I].num_iso_H[0]
#define ISOLATED_ATOM       (-15)
#define NUM_ISO_Hk(AT,I,K)  AT[I].num_iso_H[K+1]
#define IS_METAL_ATOM(AT,I) is_element_a_metal( AT[I].elname )

#else 

#define inchi_Atom          inp_ATOM
#define AT_NUM_BONDS(AT)    (AT).valence
#define ATOM_NUMBER         AT_NUMB
#define IN_NEIGH_LIST       is_in_the_list
#define inchi_NUMH2(AT,N)   NUMH(AT,N)
#define INChITo_Atom        cc_INChIToInpAtom
#define INPUT_FILE          FILE
#define Create_Atom         CreateInpAtom 
#define AT_BONDS_VAL(AT,I)  AT[I].chem_bonds_valence
#define ISOLATED_ATOM       15
#define NUM_ISO_Hk(AT,I,K)  AT[I].num_iso_H[K]
#define IS_METAL_ATOM(AT,I) is_el_a_metal( AT[I].el_number )

#endif

/*****************************************************************************/
/* local prototypes */
char *FindToken( INPUT_FILE *inp_molfile, int *bTooLongLine, const char *sToken, int lToken,
                        char *szLine, int nLenLine, char *p, int *res );
char *LoadLine( INPUT_FILE *inp_molfile, int *bTooLongLine, int *bItemIsOver, char **s,
                        char *szLine, int nLenLine, int nMinLen2Load, char *p, int *res ); 


int INChITo_Atom(INPUT_FILE *inp_molfile, MOL_COORD **szCoord,
                      inchi_Stereo0D **stereo0D, int *num_stereo0D,
                      int bDoNotAddH, INPUT_TYPE nInputType, inchi_Atom **at,
                      int max_num_at,
                      int *num_dimensions, int *num_bonds, char *pSdfLabel, char *pSdfValue,
                      long *Id, INCHI_MODE *pInpAtomFlags, int *err, char *pStrErr );


#if( defined(INCHI_LIBRARY) || defined(INCHI_MAIN) )
/*****************************************************************************/
int INChIToInchi_Atom ( INCHI_FILE *inp_molfile, inchi_Stereo0D **stereo0D, int *num_stereo0D,
                      int bDoNotAddH, INPUT_TYPE nInputType, inchi_Atom **at,
                      int max_num_at,
                      int *num_dimensions, int *num_bonds, char *pSdfLabel, char *pSdfValue,
                      long *Id, INCHI_MODE *pInpAtomFlags, int *err, char *pStrErr );

int INChIToInchi_Atom ( INCHI_FILE *inp_molfile, inchi_Stereo0D **stereo0D, int *num_stereo0D,
                      int bDoNotAddH, INPUT_TYPE nInputType, inchi_Atom **at,
                      int max_num_at,
                      int *num_dimensions, int *num_bonds, char *pSdfLabel, char *pSdfValue,
                      long *Id, INCHI_MODE *pInpAtomFlags, int *err, char *pStrErr )
{
    return INChITo_Atom ( inp_molfile, NULL, stereo0D, num_stereo0D,
                          bDoNotAddH, nInputType, at, max_num_at,
                          num_dimensions, num_bonds, pSdfLabel, pSdfValue,
                          Id, pInpAtomFlags, err, pStrErr );
}

#else

/*****************************************************************************/
int INChIToInpAtom (  FILE *inp_molfile, MOL_COORD **szCoord,
                      int bDoNotAddH, INPUT_TYPE nInputType, inp_ATOM **at,
                      int max_num_at,
                      int *num_dimensions, int *num_bonds, char *pSdfLabel, char *pSdfValue,
                      long *Id, INCHI_MODE *pInpAtomFlags, int *err, char *pStrErr );

int INChIToInpAtom (  FILE *inp_molfile, MOL_COORD **szCoord,
                      int bDoNotAddH, INPUT_TYPE nInputType, inp_ATOM **at,
                      int max_num_at,
                      int *num_dimensions, int *num_bonds, char *pSdfLabel, char *pSdfValue,
                      long *Id, INCHI_MODE *pInpAtomFlags, int *err, char *pStrErr )
{
    return INChITo_Atom ( inp_molfile, szCoord, NULL, NULL,
                          bDoNotAddH, nInputType, at, max_num_at,
                          num_dimensions, num_bonds, pSdfLabel, pSdfValue,
                          Id, pInpAtomFlags, err, pStrErr );
}

#endif

/*****************************************************************************/
char *FindToken( INPUT_FILE *inp_molfile, int *bTooLongLine, const char *sToken, int lToken,
                        char *szLine, int nLenLine, char *p, int *res )
{
    char *q;
    int   res2;
                
    while ( !(q = strstr( p, sToken ) ) ) {
        if ( (q = strrchr( p, '/' )) && (q + lToken > szLine + *res) ) {
            *res -= q - szLine; /* res = the length of the szLine to be left in */
            memmove( szLine, q, *res + 1);
        } else {
            *res = 0;
        }
        if ( !*bTooLongLine || 
             0 > (res2 = my_fgetsTab1( szLine + *res, nLenLine - *res - 1,
                                       inp_molfile, bTooLongLine ) ) ) {
            /* the line is over or end of file */
            return NULL;
        } else {
            *res += res2;
            p = szLine;
        }
    }

    return q + lToken;
}
/*****************************************************************************/
char *LoadLine( INPUT_FILE *inp_molfile, int *bTooLongLine, int *bItemIsOver, char **s,
                        char *szLine, int nLenLine, int nMinLen2Load, char *p, int *res ) 
{
    int pos = p - szLine, res2;
    if ( !*bItemIsOver && nLenLine - (*res - pos) > nMinLen2Load ) {
        /* load the next portion if possible */
        if ( pos ) {
            *res -= pos;
            memmove( szLine, p, *res+1 );
            p = szLine;
            if ( *s ) {
                *s -= pos;
            }
            pos = 0;
        }
        res2 = my_fgetsTab1( szLine + *res, nLenLine - *res - 1, inp_molfile, bTooLongLine );
        if ( res2 > 0 ) {
            *bItemIsOver = ( (*s = strchr( p + *res, '/') ) || !*bTooLongLine );
            *res += res2;
        } else {
            *bItemIsOver = 1;
        }
    }
    return p;
}
/*****************************************************************************/
int INChITo_Atom(INPUT_FILE *inp_molfile, MOL_COORD **szCoord,
                      inchi_Stereo0D **stereo0D, int *num_stereo0D,
                      int bDoNotAddH, INPUT_TYPE nInputType, inchi_Atom **at,
                      int max_num_at,
                      int *num_dimensions, int *num_bonds, char *pSdfLabel, char *pSdfValue,
                      long *Id, INCHI_MODE *pInpAtomFlags, int *err, char *pStrErr )
{
    int      num_atoms = 0, bFindNext = 0, len, bHeaderRead, bItemIsOver, bErrorMsg, bRestoreInfo;
    int      bFatal = 0, num_struct = 0;
    int      i, k, k2, res, bond_type, bond_stereo1, bond_stereo2, bond_char, neigh, bond_parity, bond_parityNM;
    int      bTooLongLine, res2, bTooLongLine2, pos, hlen, hk;
    long     longID;
    char     szLine[INCHI_LINE_LEN], szNextLine[INCHI_LINE_ADD], *p, *q, *s, parity;
    int      b2D=0, b3D=0, b23D, nNumBonds = 0, bNonZeroXYZ, bNonMetal;
    int      len_stereo0D = 0, max_len_stereo0D = 0;
    inchi_Stereo0D  *atom_stereo0D = NULL;
    inchi_Atom      *atom          = NULL;
    MOL_COORD       *pszCoord      = NULL;
    INCHI_MODE InpAtomFlags = 0; /* 0 or FLAG_INP_AT_NONCHIRAL or FLAG_INP_AT_CHIRAL */
    static const char szIsoH[] = "hdt";
    /* plain tags */
    static const char sStructHdrPln[]         = "Structure:";
    static const char sStructHdrPlnNoLblVal[] = " is missing";
    static char sStructHdrPlnAuxStart[64] =""; /*"$1.1Beta/";*/
    static int  lenStructHdrPlnAuxStart = 0;
    static const char sStructHdrPlnRevAt[]    = "/rA:";
    static const char sStructHdrPlnRevBn[]    = "/rB:";
    static const char sStructHdrPlnRevXYZ[]   = "/rC:";
    const  char *sToken;
    int  lToken;
    if ( !lenStructHdrPlnAuxStart ) {
        lenStructHdrPlnAuxStart = sprintf( sStructHdrPlnAuxStart, "AuxInfo=" );
    }

    if ( at ) {
        if ( *at && max_num_at ) {
            memset( *at, 0, max_num_at * sizeof(**at) );
        }
#if( defined(INCHI_LIBRARY) || defined(INCHI_MAIN) )
        if ( stereo0D && num_stereo0D ) {
            if ( *stereo0D && *num_stereo0D ) {
                max_len_stereo0D = *num_stereo0D;
                memset( *stereo0D, 0, max_len_stereo0D * sizeof( **stereo0D ));
            } else {
                max_len_stereo0D = 0;
            }
        }
#else
        if ( szCoord && *szCoord ) {
            inchi_free( *szCoord );
            *szCoord = NULL;
        }
#endif
    } else {
        bFindNext = 1;
    }
    bHeaderRead = bErrorMsg = bRestoreInfo = 0;
    *num_dimensions = *num_bonds = 0;

    /*************************************************************/
    /*   extract reversibility info from plain text INChI format */
    /*************************************************************/
    if ( nInputType == INPUT_INCHI_PLAIN ) {
        bHeaderRead = hk = 0;
        while ( 0 < (res = my_fgetsTab( szLine, sizeof(szLine)-1, inp_molfile, &bTooLongLine ) ) ) {

            /********************* find and interpret structure header ************/
            if ( !bTooLongLine &&
                 (hlen=sizeof(sStructHdrPln)-1, !memcmp(szLine, sStructHdrPln, hlen)) ) {
                p = szLine + hlen;
                longID = 0;
                num_atoms = 0;
                /* structure number */
                longID = strtol( p, &q, 10 );
                if ( q && q[0] == '.' && q[1] == ' ' ) {
                    p = q+2;
                }
                p = p + strspn( p, " \n\r" );
                
                if ( pSdfLabel ) {
                    pSdfLabel[0] = '\0';
                }
                if ( pSdfValue ) {
                    pSdfValue[0] = '\0';
                }

                if ( *p ) {
                    /* has label name */
                    /*p ++;*/
                    if ( q = strchr( p, '=' ) ) {
                        /* '=' separates label name from the value */
                        len = inchi_min( q-p+1, MAX_SDF_HEADER-1);
                        if ( pSdfLabel ) {
                            mystrncpy( pSdfLabel, p, len );
                            LtrimRtrim( pSdfLabel, &len );
                        }
                        p = q+1;
                        q = p + (int)strlen( p );
                        if ( q-p > 0 ) {
                            len = inchi_min( q-p+1, MAX_SDF_VALUE-1);
                            if ( pSdfValue ) {
                                mystrncpy( pSdfValue, p, len );
                            }
                            p = q;
                        }

                    } else
                    if ( q = strstr( p, sStructHdrPlnNoLblVal ) ) {
                        len = inchi_min( q-p+1, MAX_SDF_HEADER-1);
                        if ( pSdfLabel ) {
                            mystrncpy( pSdfLabel, p, len );
                        }
                        p = q+1;
                    }
                }
                if ( Id )
                    *Id = longID;

                bHeaderRead = 1;
                bErrorMsg = bRestoreInfo = 0;
            } else
            if ( !memcmp( szLine, sStructHdrPlnAuxStart, lenStructHdrPlnAuxStart) ) {
                /* found the header of the AuxInfo, read AuxInfo head of the line */
                if ( !bHeaderRead ) {
                    longID = 0;
                    if ( Id )
                        *Id = longID;
                    if ( pSdfLabel ) {
                        pSdfLabel[0] = '\0';
                    }
                    if ( pSdfValue ) {
                        pSdfValue[0] = '\0';
                    }
                }
                bHeaderRead = 0;
                /* check for empty "AuxInfo=ver//" */
                p = strchr( szLine + lenStructHdrPlnAuxStart, '/' );
                if ( p && p[1] == '/' && (!p[2] || '\n' == p[2]) ) {
                    goto bypass_end_of_INChI_plain;
                }
                /***************** search for atoms block (plain) **********************/
                p = szLine;
                sToken = sStructHdrPlnRevAt;
                lToken = sizeof(sStructHdrPlnRevAt)-1;
                /* search for sToken in the line; load next segments of the line if sToken has not found */
                p = FindToken( inp_molfile, &bTooLongLine, sToken, lToken,
                               szLine, sizeof(szLine), p, &res );
                if ( !p ) {
                    *err      = INCHI_INP_ERROR_ERR;
                    num_atoms = INCHI_INP_ERROR_RET;
                    MOLFILE_ERR_SET (*err, 0, "Missing atom data");
                    goto bypass_end_of_INChI_plain;
                } else {
                    /* atoms block started */
                    i = 0;
                    res2 = bTooLongLine2 = -1;
                    bItemIsOver = (s = strchr( p, '/') ) || !bTooLongLine;
                    while ( 1 ) {
                        p = LoadLine( inp_molfile, &bTooLongLine, &bItemIsOver, &s,
                                      szLine, sizeof(szLine), INCHI_LINE_ADD, p, &res );
                        if ( !i ) {
                            /* allocate atom */
                            num_atoms = strtol( p, &q, 10 );
                            if ( !num_atoms || !q || !*q ) {
                                num_atoms = 0; /* no atom data */
                                goto bypass_end_of_INChI_plain;
                            }
                            p = q;
                            /* Molfile chirality flag */
                            switch( *p ) {
                            case 'c':
                                InpAtomFlags |= FLAG_INP_AT_CHIRAL;
                                p ++;
                                break;
                            case 'n':
                                InpAtomFlags |= FLAG_INP_AT_NONCHIRAL;
                                p ++;
                                break;
                            }
                            if ( at && *at ) {
                                if ( num_atoms > max_num_at ) {
                                    inchi_free( *at );
                                    *at = NULL;
                                } else {
                                    memset( *at, 0, max_num_at * sizeof( **at ) );
                                    atom = *at;
                                }
                            }
                            if ( !at || !*at ) {
                                atom = Create_Atom( num_atoms+1 );
                                if ( !atom ) {
                                    num_atoms = INCHI_INP_FATAL_RET; /* was -1; error */
                                    *err      = INCHI_INP_FATAL_ERR;
                                    MOLFILE_ERR_SET (*err, 0, "Out of RAM");
                                    goto bypass_end_of_INChI_plain;
                                }
                            }
                            if ( stereo0D && *stereo0D ) {
                                if ( num_atoms > max_len_stereo0D ) {
                                    FreeInchi_Stereo0D( stereo0D );
                                } else {
                                    memset( *stereo0D, 0, max_len_stereo0D * sizeof( **stereo0D ) );
                                    atom_stereo0D = *stereo0D;
                                }
                            }
                            if ( !stereo0D || !*stereo0D ) {
                                max_len_stereo0D = num_atoms+1;
                                atom_stereo0D = CreateInchi_Stereo0D( max_len_stereo0D );
                                if ( !atom_stereo0D ) {
                                    num_atoms = INCHI_INP_FATAL_RET; /* fatal error: cannot allocate */
                                    *err      = INCHI_INP_FATAL_ERR;
                                    MOLFILE_ERR_SET (*err, 0, "Out of RAM");
                                    goto bypass_end_of_INChI_plain;
                                }
                            }
                        }
                        /* element, first char */
                        if ( !isalpha( UCINT *p ) || !isupper( UCINT *p ) || i >= num_atoms ) {
                            break; /* end of atoms block */
                        }
                        atom[i].elname[0] = *p ++;
                        /* element, second char */
                        if ( isalpha( UCINT *p ) && islower( UCINT *p ) ) {
                            atom[i].elname[1] = *p ++;
                        }
    #if( defined(INCHI_LIBRARY) || defined(INCHI_MAIN) )
    #else
                        atom[i].el_number = get_periodic_table_number( atom[i].elname );
    #endif
                        /* bonds' valence + number of non-isotopic H */
                        if ( isdigit( UCINT *p ) ) {
                            AT_BONDS_VAL(atom,i) = (char)strtol( p, &q, 10 );
                            if ( !AT_BONDS_VAL(atom,i) )
                                AT_BONDS_VAL(atom,i) = ISOLATED_ATOM; /* same convention as in MOLfile, found zero bonds valence */
                            p = q;
                        }
                        /* charge */
                        atom[i].charge = (*p == '+')? 1 : (*p == '-')? -1 : 0;
                        if ( atom[i].charge ) {
                            p ++;
                            if ( isdigit( UCINT *p ) ) {
                                atom[i].charge *= (S_CHAR)(strtol( p, &q, 10 ) & CHAR_MASK);
                                p = q;
                            }
                        }
                        /* radical */
                        if ( *p == '.' ) {
                            p ++;
                            if ( isdigit( UCINT *p ) ) {
                                atom[i].radical = (S_CHAR)strtol( p, &q, 10 );
                                p = q;
                            }
                        }
                        /* isotopic mass */
                        if ( *p == 'i' ) {
                            p ++;
                            if ( isdigit( UCINT *p ) ) {
                                int mw = strtol( p, &q, 10 );
                                p = q;
    #if( defined(INCHI_LIBRARY) || defined(INCHI_MAIN) )
                                atom[i].isotopic_mass = mw;
    #else
                                mw -= get_atw_from_elnum( atom[i].el_number );
                                if ( mw >= 0 )
                                    mw ++;
                                atom[i].iso_atw_diff = mw;
    #endif
                            }
                        }
                        /* parity */
                        switch( *p ) {
                        case 'o':
                            parity = INCHI_PARITY_ODD;
                            p ++;
                            break;
                        case 'e':
                            parity = INCHI_PARITY_EVEN;
                            p ++;
                            break;
                        case 'u':
                            parity = INCHI_PARITY_UNKNOWN;
                            p ++;
                            break;
                        case '?':
                            parity = INCHI_PARITY_UNDEFINED;
                            p ++;
                            break;
                        default:
                            parity = 0;
                            break;
                        }
                        if ( parity ) {
                            atom_stereo0D[len_stereo0D].central_atom = i;
                            atom_stereo0D[len_stereo0D].parity       = parity;
                            atom_stereo0D[len_stereo0D].type         = INCHI_StereoType_Tetrahedral;
                            len_stereo0D ++;
                        }
                        /* isotopic h, d, t */
                        for ( k = 0; k < NUM_H_ISOTOPES; k ++ ) {
                            if ( *p == szIsoH[k] ) {
                                NUM_ISO_Hk(atom,i,k) = 1;
                                p ++;
                                if ( isdigit( UCINT *p ) ) {
                                    NUM_ISO_Hk(atom,i,k) = (char)strtol( p, &q, 10 );
                                    p = q;
                                }
                            }
                        }
                        i ++;
                    }
                    if ( !bItemIsOver || i != num_atoms || s && p != s ) {
                        num_atoms = INCHI_INP_ERROR_RET; /* error */
                        *err      = INCHI_INP_ERROR_ERR;
                        MOLFILE_ERR_SET (*err, 0, "Wrong number of atoms");
                        goto bypass_end_of_INChI_plain;
                    }
                }
                /***************** search for bonds block (plain) and read it *****************/
                /*p = szLine;*/
                sToken = sStructHdrPlnRevBn;
                lToken = sizeof(sStructHdrPlnRevBn)-1;
                /* search for sToken in the line; load next segments of the line if sToken has not found */
                p = FindToken( inp_molfile, &bTooLongLine, sToken, lToken,
                               szLine, sizeof(szLine), p, &res );
                if ( !p ) {
                    num_atoms = INCHI_INP_ERROR_RET; /* error */
                    *err      = INCHI_INP_ERROR_ERR;
                    MOLFILE_ERR_SET (*err, 0, "Missing bonds data");
                    goto bypass_end_of_INChI_plain;
                } else {
                    /* bonds block started */
                    i = 1;
                    res2 = bTooLongLine2 = -1;
                    bItemIsOver = (s = strchr( p, '/') ) || !bTooLongLine;
                    if ( 1 == num_atoms ) {
                        /* needed because the next '/' may be still out of szLine */
                        p = LoadLine( inp_molfile, &bTooLongLine, &bItemIsOver, &s,
                                      szLine, sizeof(szLine), INCHI_LINE_ADD, p, &res );
                    }
                    while ( i < num_atoms ) {
                        p = LoadLine( inp_molfile, &bTooLongLine, &bItemIsOver, &s,
                                      szLine, sizeof(szLine), INCHI_LINE_ADD, p, &res );
                        if ( i >= num_atoms || s && p >= s ) {
                            break; /* end of bonds (plain) */
                        }
                        /* bond, first char */
                        if ( *p == ';' ) {
                            p ++;
                            i ++;
                            continue;
                        }
                        if ( !isalpha( UCINT *p ) ) {
                            num_atoms = INCHI_INP_ERROR_RET; /* error */
                            *err      = INCHI_INP_ERROR_ERR;
                            MOLFILE_ERR_SET (*err, 0, "Wrong bonds data");
                            goto bypass_end_of_INChI_plain;
                        }
                        bond_char = *p ++;
                        /* bond parity */
                        switch( *p ) {
                        case '-':
                            bond_parity = INCHI_PARITY_ODD;
                            p ++;
                            break;
                        case '+':
                            bond_parity = INCHI_PARITY_EVEN;
                            p ++;
                            break;
                        case 'u':
                            bond_parity = INCHI_PARITY_UNKNOWN;
                            p ++;
                            break;
                        case '?':
                            bond_parity = INCHI_PARITY_UNDEFINED;
                            p ++;
                            break;
                        default:
                            bond_parity = 0;
                            break;
                        }
                        if ( bond_parity ) {
                            switch( *p ) {
                            case '-':
                                bond_parityNM = INCHI_PARITY_ODD;
                                p ++;
                                break;
                            case '+':
                                bond_parityNM = INCHI_PARITY_EVEN;
                                p ++;
                                break;
                            case 'u':
                                bond_parityNM = INCHI_PARITY_UNKNOWN;
                                p ++;
                                break;
                            case '?':
                                bond_parityNM = INCHI_PARITY_UNDEFINED;
                                p ++;
                                break;
                            default:
                                bond_parityNM = 0;
                                break;
                            }
                        } else {
                            bond_parityNM = 0;
                        }

                        /* neighbor of the current atom */
                        if ( !isdigit( UCINT *p ) ) {
                            num_atoms = INCHI_INP_ERROR_RET; /* error */
                            *err      = INCHI_INP_ERROR_ERR;
                            MOLFILE_ERR_SET (*err, 0, "Wrong bonds data");
                            goto bypass_end_of_INChI_plain;
                        }
                        neigh = (int)strtol( p, &q, 10 )-1;

                        if ( i >= num_atoms || neigh >= num_atoms ) {
                            num_atoms = INCHI_INP_ERROR_RET; /* error */
                            *err      = INCHI_INP_ERROR_ERR;
                            MOLFILE_ERR_SET (*err, 0, "Bond to nonexistent atom");
                            goto bypass_end_of_INChI_plain;
                        }
                        p = q;
                        bond_stereo1 = bond_stereo2 = 0;

                        /* bond type & 2D stereo */
                        switch( bond_char ) {
                        case 'v':
                            bond_type    = INCHI_BOND_TYPE_SINGLE;
                            bond_stereo1 = INCHI_BOND_STEREO_SINGLE_1EITHER;
                            bond_stereo2 = INCHI_BOND_STEREO_SINGLE_2EITHER;
                            break;
                        case 'V':
                            bond_type    = INCHI_BOND_TYPE_SINGLE;
                            bond_stereo1 = INCHI_BOND_STEREO_SINGLE_2EITHER;
                            bond_stereo2 = INCHI_BOND_STEREO_SINGLE_1EITHER;
                            break;
                        case 'w':
                            bond_type    = INCHI_BOND_TYPE_DOUBLE;
                            bond_stereo1 = 
                            bond_stereo2 = INCHI_BOND_STEREO_DOUBLE_EITHER;
                            break;
                        case 's':
                            bond_type    = INCHI_BOND_TYPE_SINGLE;
                            break;
                        case 'd':
                            bond_type    = INCHI_BOND_TYPE_DOUBLE;
                            break;
                        case 't':
                            bond_type    = INCHI_BOND_TYPE_TRIPLE;
                            break;
                        case 'a':
                            bond_type    = INCHI_BOND_TYPE_ALTERN;
                            break;
                        case 'p':
                            bond_type    =  INCHI_BOND_TYPE_SINGLE;
                            bond_stereo1 =  INCHI_BOND_STEREO_SINGLE_1UP;
                            bond_stereo2 =  INCHI_BOND_STEREO_SINGLE_2UP;
                            break;
                        case 'P':
                            bond_type    =  INCHI_BOND_TYPE_SINGLE;
                            bond_stereo1 =  INCHI_BOND_STEREO_SINGLE_2UP;
                            bond_stereo2 =  INCHI_BOND_STEREO_SINGLE_1UP;
                            break;
                        case 'n':
                            bond_type    =  INCHI_BOND_TYPE_SINGLE;
                            bond_stereo1 =  INCHI_BOND_STEREO_SINGLE_1DOWN;
                            bond_stereo2 =  INCHI_BOND_STEREO_SINGLE_2DOWN;
                            break;
                        case 'N':
                            bond_type    =  INCHI_BOND_TYPE_SINGLE;
                            bond_stereo1 =  INCHI_BOND_STEREO_SINGLE_2DOWN;
                            bond_stereo2 =  INCHI_BOND_STEREO_SINGLE_1DOWN;
                            break;
                        default:
                            num_atoms = INCHI_INP_ERROR_RET; /* error */
                            *err      = INCHI_INP_ERROR_ERR;
                            MOLFILE_ERR_SET (*err, 0, "Wrong bond type");
                            goto bypass_end_of_INChI_plain;
                        }
                        k = AT_NUM_BONDS(atom[i]) ++;
                        atom[i].bond_type[k]   = bond_type;
                        atom[i].bond_stereo[k] = bond_stereo1;
                        atom[i].neighbor[k]    = (ATOM_NUMBER)neigh;
                        k2 = AT_NUM_BONDS(atom[neigh]) ++;
                        atom[neigh].bond_type[k2]   = bond_type;
                        atom[neigh].bond_stereo[k2] = bond_stereo2;
                        atom[neigh].neighbor[k2]    = (ATOM_NUMBER)i;
                        bond_parity |= (bond_parityNM << SB_PARITY_SHFT);

                        if ( bond_parity ) {
                            if ( max_len_stereo0D <= len_stereo0D ) {
                                /* realloc atom_Stereo0D */
                                inchi_Stereo0D *new_atom_stereo0D = CreateInchi_Stereo0D( max_len_stereo0D+num_atoms );
                                if ( !new_atom_stereo0D ) {
                                    num_atoms = INCHI_INP_FATAL_RET; /* fatal error: cannot allocate */
                                    *err      = INCHI_INP_FATAL_ERR;
                                    MOLFILE_ERR_SET (*err, 0, "Out of RAM");
                                    goto bypass_end_of_INChI_plain;
                                }
                                memcpy( new_atom_stereo0D, atom_stereo0D, len_stereo0D * sizeof(*atom_stereo0D) );
                                FreeInchi_Stereo0D( &atom_stereo0D );
                                atom_stereo0D = new_atom_stereo0D;
                                max_len_stereo0D += num_atoms;
                            }
                            /* (a) i may be allene endpoint and     neigh = allene middle point or
                               (b) i may be allene middle point and neigh = allene endpoint
                               !!!!! CURRENTLY ONLY (b) IS ALLOWED !!!!!
                            */
                            atom_stereo0D[len_stereo0D].neighbor[1] = neigh; /* neigh < i */
                            atom_stereo0D[len_stereo0D].neighbor[2] = i;
                            atom_stereo0D[len_stereo0D].parity      = bond_parity;
                            atom_stereo0D[len_stereo0D].type        = INCHI_StereoType_DoubleBond; /* incl allenes & cumulenes */
                            len_stereo0D ++;
                        }
                    }
                    if ( !bItemIsOver || i != num_atoms || s && p != s ) {
                        num_atoms = INCHI_INP_ERROR_RET; /* error */
                        *err      = INCHI_INP_ERROR_ERR;
                        MOLFILE_ERR_SET (*err, 0, "Wrong number of bonds");
                        goto bypass_end_of_INChI_plain;
                    }
                }
                /***************** search for coordinates block (plain) **********************/
                /*p = szLine;*/
                sToken = sStructHdrPlnRevXYZ;
                lToken = sizeof(sStructHdrPlnRevXYZ)-1;
                /* search for sToken in the line; load next segments of the line if sToken has not found */
                p = FindToken( inp_molfile, &bTooLongLine, sToken, lToken,
                               szLine, sizeof(szLine), p, &res );
                if ( !p ) {
                    num_atoms = INCHI_INP_ERROR_RET; /* error */
                    *err      = INCHI_INP_ERROR_ERR;
                    MOLFILE_ERR_SET (*err, 0, "Missing atom coordinates data");
                    goto bypass_end_of_INChI_plain;
                } else {
                    /* coordinates block started */
                    if ( pszCoord = (MOL_COORD*)inchi_malloc(inchi_max(num_atoms,1) * sizeof(MOL_COORD)) ) {
                        memset( pszCoord, ' ', inchi_max(num_atoms,1) * sizeof(MOL_COORD));
                    } else {
                        num_atoms = INCHI_INP_FATAL_RET; /* allocation error */
                        *err      = INCHI_INP_FATAL_ERR;
                        MOLFILE_ERR_SET (*err, 0, "Out of RAM");
                        goto bypass_end_of_INChI_plain;
                    }
                    i = 0;
                    res2 = bTooLongLine2 = -1;
                    bItemIsOver = (s = strchr( p, '/') ) || !bTooLongLine;
                    while ( i < num_atoms ) {
                        p = LoadLine( inp_molfile, &bTooLongLine, &bItemIsOver, &s,
                                      szLine, sizeof(szLine), INCHI_LINE_ADD, p, &res );
                        if ( i >= num_atoms || s && p >= s ) {
                            break; /* end of bonds (plain) */
                        }

                        /* coord, first char */
                        if ( *p == ';' ) {
                            for ( k = 0; k < NUM_COORD; k ++ ) {
                                pszCoord[i][LEN_COORD*k + 4] = '0';
                            }
                            p ++;
                            i ++;
                            continue;
                        }
                        for ( k = 0; k < 3; k ++ ) {
                            double xyz;
                            bNonZeroXYZ = 0;
                            if ( *p == ';' ) {
                                pszCoord[i][LEN_COORD*k + 4] = '0';
                                xyz = 0.0;
                            } else
                            if ( *p == ',' ) {
                                /* empty */
                                pszCoord[i][LEN_COORD*k + 4] = '0';
                                xyz = 0.0;
                                p ++;
                            } else {
                                xyz = strtod( p, &q );
                                bNonZeroXYZ = fabs(xyz) > MIN_BOND_LENGTH;
                                if ( q != NULL ) {
                                    memcpy( pszCoord[i]+LEN_COORD*k, p, q-p );
                                    if ( *q == ',' )
                                        q ++;
                                    p = q;
                                } else {
                                    pszCoord[i][LEN_COORD*k + 4] = '0';
                                }
                            }
                            switch( k ) {
                            case 0:
                                atom[i].x = xyz;
                                b2D |= bNonZeroXYZ;
                                break;
                            case 1:
                                atom[i].y = xyz;
                                b2D |= bNonZeroXYZ;
                                break;
                            case 2:
                                b3D |= bNonZeroXYZ;
                                atom[i].z = xyz;
                                break;
                            }
                        }
                        if ( *p == ';' ) {
                            p ++; /* end of this triple of coordinates */
                            i ++;
                        } else {
                            num_atoms = INCHI_INP_ERROR_RET; /* error in input data: atoms, bonds & coord must be present together */
                            *err      = INCHI_INP_ERROR_ERR;
                            MOLFILE_ERR_SET (*err, 0, "Wrong atom coordinates data");
                            goto bypass_end_of_INChI_plain;
                        }
                    }
                    if ( !bItemIsOver || s && p != s || i != num_atoms ) {
                        num_atoms = INCHI_INP_ERROR_RET; /* error */
                        *err      = INCHI_INP_ERROR_ERR;
                        MOLFILE_ERR_SET (*err, 0, "Wrong number of coordinates");
                        goto bypass_end_of_INChI_plain;
                    }
                } /* end of coordinates */
                /* set special valences and implicit H (xml) */
                b23D = b2D | b3D;
                b2D = b3D = 0;
                if ( at ) {
                    if ( !*at ) {
                        int a1, a2, n1, n2, valence;
                        int chem_bonds_valence;
                        int    nX=0, nY=0, nZ=0, nXYZ;
                        *at = atom;
                        /* special valences */
                        for ( bNonMetal = 0; bNonMetal < 1; bNonMetal ++ ) {
                            for ( a1 = 0; a1 < num_atoms; a1 ++ ) {
                                int num_bond_type[MAX_INPUT_BOND_TYPE - MIN_INPUT_BOND_TYPE + 1];
#if( defined(INCHI_LIBRARY) || defined(INCHI_MAIN) )
#else
                                int bHasMetalNeighbor=0;
#endif
                                memset( num_bond_type, 0, sizeof(num_bond_type) );

                                valence = AT_BONDS_VAL(atom, a1); /*  save atom valence if available */
                                AT_BONDS_VAL(atom, a1) = 0;
#if( defined(INCHI_LIBRARY) || defined(INCHI_MAIN) )
#else
                                atom[a1].orig_at_number = a1+1;
#endif
                                nX = nY = nZ = 0;
                                for ( n1 = 0; n1 < AT_NUM_BONDS(atom[a1]); n1 ++ ) {
                                    bond_type = atom[a1].bond_type[n1] - MIN_INPUT_BOND_TYPE;
                                    if (  bond_type < 0 || bond_type > MAX_INPUT_BOND_TYPE - MIN_INPUT_BOND_TYPE ) {
                                        bond_type = 0;
                                        MOLFILE_ERR_SET (*err, 0, "Unknown bond type in InChI aux assigned as a single bond");
                                    }

                                    num_bond_type[ bond_type ] ++;
                                    nNumBonds ++;
                                    if ( b23D ) {
                                        neigh = atom[a1].neighbor[n1];
                                        nX |= (fabs(atom[a1].x - atom[neigh].x) > MIN_BOND_LENGTH);
                                        nY |= (fabs(atom[a1].y - atom[neigh].y) > MIN_BOND_LENGTH);
                                        nZ |= (fabs(atom[a1].z - atom[neigh].z) > MIN_BOND_LENGTH);
                                    }
                                }
                                chem_bonds_valence = 0;
                                for ( n1 = 0; MIN_INPUT_BOND_TYPE + n1 <= 3 && MIN_INPUT_BOND_TYPE + n1 <= MAX_INPUT_BOND_TYPE; n1 ++ ) {
                                    chem_bonds_valence += (MIN_INPUT_BOND_TYPE + n1) * num_bond_type[n1];
                                }
                                if ( MIN_INPUT_BOND_TYPE <= INCHI_BOND_TYPE_ALTERN && INCHI_BOND_TYPE_ALTERN <= MAX_INPUT_BOND_TYPE &&
                                     ( n2 = num_bond_type[INCHI_BOND_TYPE_ALTERN-MIN_INPUT_BOND_TYPE] ) ) {
                                    /* accept input aromatic bonds for now */
                                    switch ( n2 ) {
                                    case 2:
                                        chem_bonds_valence += 3;  /* =A- */
                                        break;
                                    case 3:
                                        chem_bonds_valence += 4;  /* =A< */
                                        break;
                                    default:
                                        /*  if 1 or >= 4 aromatic bonds then replace such bonds with single bonds */
                                        for ( n1 = 0; n1 < AT_NUM_BONDS(atom[a1]); n1 ++ ) {
                                            if ( atom[a1].bond_type[n1] == INCHI_BOND_TYPE_ALTERN ) {
                                                ATOM_NUMBER *p1;
                                                a2 = atom[a1].neighbor[n1];
                                                p1 = IN_NEIGH_LIST( atom[a2].neighbor, (ATOM_NUMBER)a1, AT_NUM_BONDS(atom[a2]) );
                                                if ( p1 ) {
                                                    atom[a1].bond_type[n1] = 
                                                    atom[a2].bond_type[p1-atom[a2].neighbor] = INCHI_BOND_TYPE_SINGLE;
                                                } else {
                                                    *err = -2;  /*  Program error */
                                                    MOLFILE_ERR_SET (*err, 0, "Program error interpreting InChI aux");
                                                    num_atoms = 0;
                                                    goto bypass_end_of_INChI_plain; /*  no structure */
                                                }
                                            }
                                        }
                                        chem_bonds_valence += n2;
                                        *err |= 32; /*  Unrecognized aromatic bond(s) replaced with single */
                                        MOLFILE_ERR_SET (*err, 0, "Atom has 1 or more than 3 aromatic bonds");
                                        break;
                                    }
                                }
#if( defined(INCHI_LIBRARY) || defined(INCHI_MAIN) )
                                /*************************************************************************************
                                 *
                                 *  Set number of hydrogen atoms
                                 */
                                {
                                    int num_iso_H;
                                    num_iso_H = atom[a1].num_iso_H[1] + atom[a1].num_iso_H[2] + atom[a1].num_iso_H[3];
                                    if ( valence == ISOLATED_ATOM ) {
                                        atom[a1].num_iso_H[0] = 0;
                                    } else
                                    if ( valence && valence >= chem_bonds_valence ) {
                                        atom[a1].num_iso_H[0] = valence - chem_bonds_valence;
                                    } else
                                    if ( valence || bDoNotAddH ) {
                                        atom[a1].num_iso_H[0] = 0;
                                    } else
                                    if ( !bDoNotAddH ) {
                                        atom[a1].num_iso_H[0] = -1; /* auto add H */
                                    }
                                }
#else                                
                                /* added 2006-07-19 to process aromatic bonds same way as from molfile */
                                if ( n2 && !valence ) {
                                    int num_H = NUMH(atom, a1); /* only isotopic */
                                    int chem_valence = chem_bonds_valence;
                                    int bUnusualValenceArom = 
                                        detect_unusual_el_valence( (int)atom[a1].el_number, atom[a1].charge,
                                                                    atom[a1].radical, chem_valence,
                                                                    num_H, atom[a1].valence );
                                    int bUnusualValenceNoArom = 
                                        detect_unusual_el_valence( (int)atom[a1].el_number, atom[a1].charge,
                                                                    atom[a1].radical, chem_valence-1,
                                                                    num_H, atom[a1].valence );
#if ( CHECK_AROMBOND2ALT == 1 )
                                    if ( bUnusualValenceArom && !bUnusualValenceNoArom && 0 == nBondsValToMetal( atom, a1) )
#else
                                    if ( bUnusualValenceArom && !bUnusualValenceNoArom )
#endif                     
                                    {
                                        /* typically NH in 5-member aromatic ring */
                                        chem_bonds_valence --;
                                    }
                                } else
                                if ( n2 && valence ) {
                                    /* atom has aromatic bonds AND the chemical valence is known */
                                    int num_H = NUMH(atom, a1);
                                    int chem_valence = chem_bonds_valence + num_H;
                                    if ( valence == chem_valence-1 ) {
                                        /* typically NH in 5-member aromatic ring */
                                        chem_bonds_valence --;
                                    }
                                }

                                atom[a1].chem_bonds_valence = chem_bonds_valence;
                                atom[a1].num_H = get_num_H( atom[a1].elname, atom[a1].num_H, atom[a1].num_iso_H, atom[a1].charge, atom[a1].radical,
                                                          atom[a1].chem_bonds_valence,
                                                          valence,
                                                          0, bDoNotAddH, bHasMetalNeighbor );
#endif
                            }
                        }
                        nNumBonds /= 2;
                        if ( b23D && nNumBonds ) {
                            nXYZ = nX+nY+nZ;
                            b2D  = (nXYZ > 0);
                            b3D  = (nXYZ == 3);
                            *num_dimensions = b3D? 3 : b2D? 2 : 0;
                            *num_bonds = nNumBonds;
                        }
                        /*======= 0D parities =================================*/
#if( defined(INCHI_LIBRARY) || defined(INCHI_MAIN) )
                        if ( len_stereo0D > 0 && atom_stereo0D && stereo0D ) {
                            *stereo0D     = atom_stereo0D;
                            *num_stereo0D = len_stereo0D;
                        } else {
                            FreeInchi_Stereo0D( &atom_stereo0D );
                            *num_stereo0D = len_stereo0D = 0;
                        }
#endif
                        for ( i = 0; i < len_stereo0D; i ++ ) {
                            ATOM_NUMBER *p1, *p2;
                            int     sb_ord_from_a1 = -1, sb_ord_from_a2 = -1, bEnd1 = 0, bEnd2 = 0;
                            switch( atom_stereo0D[i].type ) {

                            case INCHI_StereoType_Tetrahedral:
                                a1 = atom_stereo0D[i].central_atom;
                                if ( atom_stereo0D[i].parity && (AT_NUM_BONDS(atom[a1]) == 3 || AT_NUM_BONDS(atom[a1]) == 4) ) {
                                    int ii, kk = 0;
                                    if ( AT_NUM_BONDS(atom[a1]) == 3 ) {
                                        atom_stereo0D[i].neighbor[kk++] = a1;
                                    }
                                    for ( ii = 0; ii < AT_NUM_BONDS(atom[a1]); ii ++ ) {
                                        atom_stereo0D[i].neighbor[kk++] = atom[a1].neighbor[ii];
                                    }
                                }
                            
                            break;

                            case INCHI_StereoType_DoubleBond:
#define MAX_CHAIN_LEN 20
                                a1 = atom_stereo0D[i].neighbor[1];
                                a2 = atom_stereo0D[i].neighbor[2];
                                p1 = IN_NEIGH_LIST( atom[a1].neighbor, (ATOM_NUMBER)a2, AT_NUM_BONDS(atom[a1]) );
                                p2 = IN_NEIGH_LIST( atom[a2].neighbor, (ATOM_NUMBER)a1, AT_NUM_BONDS(atom[a2]) );
                                if ( !p1 || !p2 ) {
                                    atom_stereo0D[i].type = INCHI_StereoType_None;
                                    atom_stereo0D[i].central_atom = NO_ATOM;
                                    atom_stereo0D[i].neighbor[0] =
                                    atom_stereo0D[i].neighbor[3] = -1;
                                    *err |= 64; /* Error in cumulene stereo */
                                    MOLFILE_ERR_SET (*err, 0, "0D stereobond not recognized");
                                    break;
                                }
                                /* streobond, allene, or cumulene */

                                sb_ord_from_a1 = p1 - atom[a1].neighbor;
                                sb_ord_from_a2 = p2 - atom[a2].neighbor;
                                
                                if (  AT_NUM_BONDS(atom[a1]) == 2 &&
                                      atom[a1].bond_type[0] + atom[a1].bond_type[1] == 2*INCHI_BOND_TYPE_DOUBLE &&
                                      0 == inchi_NUMH2(atom, a1) &&
                                     (AT_NUM_BONDS(atom[a2]) != 2 ||
                                      atom[a2].bond_type[0] + atom[a2].bond_type[1] != 2*INCHI_BOND_TYPE_DOUBLE ) ) {
                                    bEnd2 = 1; /* a2 is the end-atom, a1 is middle atom */   
                                }
                                if (  AT_NUM_BONDS(atom[a2]) == 2 &&
                                      atom[a2].bond_type[0] + atom[a2].bond_type[1] == 2*INCHI_BOND_TYPE_DOUBLE &&
                                      0 == inchi_NUMH2(atom, a2) &&
                                     (AT_NUM_BONDS(atom[a1]) != 2 ||
                                      atom[a1].bond_type[0] + atom[a1].bond_type[1] != 2*INCHI_BOND_TYPE_DOUBLE ) ) {
                                    bEnd1 = 1; /* a1 is the end-atom, a2 is middle atom */   
                                }
                                
                                if ( bEnd2 + bEnd1 == 1 ) {
                                    /* allene or cumulene */
                                    ATOM_NUMBER  chain[MAX_CHAIN_LEN+1], prev, cur, next;
                                    if ( bEnd2 && !bEnd1 ) {
                                        cur = a1;
                                        a1 = a2;
                                        a2 = cur;
                                        sb_ord_from_a1 = sb_ord_from_a2;
                                    }
                                    sb_ord_from_a2 = -1;
                                    cur  = a1;
                                    next = a2;
                                    len = 0;
                                    chain[len++] = cur;
                                    chain[len++] = next;
                                    while ( len < MAX_CHAIN_LEN ) { /* arbitrary very high upper limit to prevent infinite loop */
                                        prev = cur;
                                        cur  = next;
                                            /* follow double bond path && avoid going back */
                                        if ( AT_NUM_BONDS(atom[cur]) == 2 &&
                                             atom[cur].bond_type[0]+atom[cur].bond_type[1] == 2*INCHI_BOND_TYPE_DOUBLE &&
                                             0 == inchi_NUMH2(atom, cur) ) {
                                            next     = atom[cur].neighbor[atom[cur].neighbor[0] == prev];
                                            chain[len++] = next;
                                        } else {
                                            break;
                                        }
                                    }
                                    if ( len > 2 &&
                                         (p2 = IN_NEIGH_LIST( atom[cur].neighbor, (ATOM_NUMBER)prev, AT_NUM_BONDS(atom[cur]))) ) {
                                        sb_ord_from_a2 = p2 - atom[cur].neighbor;
                                        a2 = cur;
                                        /* by design we need to pick up the first non-stereo-bond-neighbor as "sn"-atom */
                                        atom_stereo0D[i].neighbor[0] = atom[a1].neighbor[sb_ord_from_a1 == 0];
                                        atom_stereo0D[i].neighbor[1] = a1;
                                        atom_stereo0D[i].neighbor[2] = a2;
                                        atom_stereo0D[i].neighbor[3] = atom[a2].neighbor[sb_ord_from_a2 == 0];
                                        if ( len % 2 ) {
                                            atom_stereo0D[i].central_atom = chain[len/2];
                                            atom_stereo0D[i].type         = INCHI_StereoType_Allene; 
                                        } else {
                                            atom_stereo0D[i].central_atom = NO_ATOM;
                                        }
                                    } else {
                                        /* error */
                                        atom_stereo0D[i].type = INCHI_StereoType_None;
                                        atom_stereo0D[i].central_atom = NO_ATOM;
                                        atom_stereo0D[i].neighbor[0] =
                                        atom_stereo0D[i].neighbor[3] = -1;
                                        *err |= 64; /* Error in cumulene stereo */
                                        MOLFILE_ERR_SET (*err, 0, "Cumulene stereo not recognized (0D)");

                                    }
#undef MAX_CHAIN_LEN 
                                } else {
                                    /****** a normal possibly stereogenic bond -- not an allene or cumulene *******/
                                    /* by design we need to pick up the first non-stereo-bond-neighbor as "sn"-atom */
                                    sb_ord_from_a1 = p1 - atom[a1].neighbor;
                                    sb_ord_from_a2 = p2 - atom[a2].neighbor;
                                    atom_stereo0D[i].neighbor[0] = atom[a1].neighbor[p1 == atom[a1].neighbor];
                                    atom_stereo0D[i].neighbor[3] = atom[a2].neighbor[p2 == atom[a2].neighbor];
                                    atom_stereo0D[i].central_atom = NO_ATOM;
                                }
                                if ( atom_stereo0D[i].type != INCHI_StereoType_None &&
                                     sb_ord_from_a1 >= 0 && sb_ord_from_a2 >= 0 &&
                                     ATOM_PARITY_WELL_DEF( SB_PARITY_2(atom_stereo0D[i].parity) ) ) {
                                    /* Detected well-defined disconnected stereo
                                     * locate first non-metal neighbors */
                                    int    a, n, j, /* k,*/ sb_ord, cur_neigh, min_neigh;
                                    for ( k = 0; k < 2; k ++ ) {
                                        a      = k? atom_stereo0D[i].neighbor[2] : atom_stereo0D[i].neighbor[1];
                                        sb_ord = k? sb_ord_from_a2 : sb_ord_from_a1;
                                        min_neigh = num_atoms;
                                        for ( n  = j = 0; j < AT_NUM_BONDS(atom[a]); j ++ ) {
                                            cur_neigh = atom[a].neighbor[j];
                                            if ( j != sb_ord && !IS_METAL_ATOM(atom, cur_neigh) ) {
                                                min_neigh = inchi_min( cur_neigh, min_neigh );
                                            }
                                        }
                                        if ( min_neigh < num_atoms ) {
                                            atom_stereo0D[i].neighbor[k?3:0] = min_neigh;
                                        } else {
                                            MOLFILE_ERR_SET (*err, 0, "Cannot find non-metal stereobond neighor (0D)");
                                        }
                                    }
                                }

                                break;
                            }
                        }
                        /* end of 0D parities extraction */
/*exit_cycle:;*/
                    }
#if( defined(INCHI_LIBRARY) || defined(INCHI_MAIN) )
#else
                    /* transfer atom_stereo0D[] to atom[] */
                    if ( len_stereo0D ) {
                        Extract0DParities( atom, num_atoms, atom_stereo0D, len_stereo0D, pStrErr, err );
                    }
#endif
                    if ( pInpAtomFlags ) {
                        /* save chirality flag */
                        *pInpAtomFlags |= InpAtomFlags;
                    }
                } else
                if ( atom ) {
                    inchi_free( atom );
                    atom = NULL;
                }
#if( defined(INCHI_LIBRARY) || defined(INCHI_MAIN) )
#else
#if( FIX_READ_AUX_MEM_LEAK == 1 )
                /* 2005-08-04 avoid memory leak */
                if ( atom_stereo0D && !(stereo0D && *stereo0D == atom_stereo0D) ) {
                    FreeInchi_Stereo0D( &atom_stereo0D );
                }
#endif
                if ( szCoord ) {
                    *szCoord = pszCoord;
                    pszCoord = NULL;
                } else
#endif
                if ( pszCoord ) {
                    inchi_free( pszCoord );
                    pszCoord = NULL;
                }
                goto bypass_end_of_INChI_plain;
                /*return num_atoms;*/
            }
        }
        if ( atom_stereo0D ) {
            FreeInchi_Stereo0D( &atom_stereo0D );
        }
        /* end of struct. reading cycle */
        if ( res <= 0 ) {
            if ( *err == INCHI_INP_ERROR_ERR ) {
                return num_atoms;
            }
            *err = INCHI_INP_EOF_ERR;
            return INCHI_INP_EOF_RET; /* no more data */
        }
bypass_end_of_INChI_plain:
        /* cleanup */
        if ( num_atoms == INCHI_INP_ERROR_RET && atom_stereo0D ) {
            if ( stereo0D && *stereo0D == atom_stereo0D ) {
                *stereo0D     = NULL;
                *num_stereo0D = 0;
            }
            FreeInchi_Stereo0D( &atom_stereo0D );
        }
        while ( bTooLongLine && 
                0 < my_fgetsTab1( szLine, sizeof(szLine)-1, inp_molfile, &bTooLongLine ) ) {
            ;
        }
#if( defined(INCHI_LIBRARY) || defined(INCHI_MAIN) )
        /* cleanup */
        if ( !*at ) {
            if ( atom ) {
                inchi_free( atom );
                atom = NULL;
            }
            if ( pszCoord ) {
                inchi_free( pszCoord );
                pszCoord = NULL;
            }
        }
#endif

        return num_atoms;
    }
    
    /***********************************************************/
    /*   extract reversibility info from xml text INChI format */
    /*                                                         */
    /*   OBSOLETE CODE because InChI output in XML             */
    /*      does not exist anymore. Unsupported.               */
    /*                                                         */
    /***********************************************************/
    if ( nInputType == INPUT_INCHI_XML ) {
        /* xml tags */
        static const char sStructHdrXml[]         = "<structure";
        static const char sStructHdrXmlEnd[]      = "</structure";
        static const char sStructHdrXmlNumber[]   = "number=\"";
        static const char sStructHdrXmlIdName[]   = "id.name=\"";
        static const char sStructHdrXmlIdValue[]  = "id.value=\"";
#if( SPECIAL_BUILD == 1 )
        static const char sStructMsgXmlErr[]      = "<message type=\"error (no MoChI)\" value=\"";
#else
        static const char sStructMsgXmlErr[]      = "<message type=\"error (no InChI)\" value=\"";
#endif
        static const char sStructMsgXmlErrFatal[] = "<message type=\"fatal (aborted)\" value=\"";
        static const char sStructRevXmlRevHdr[]   = "<reversibility>";
        static const char sStructRevXmlRevAt[]    = "<atoms>";
        static const char sStructRevXmlRevAtEnd[] = "</atoms>";
        static const char sStructRevXmlRevBn[]    = "<bonds>";
        static const char sStructRevXmlRevBnEnd[] = "</bonds>";
        static const char sStructRevXmlRevXYZ[]   = "<xyz>";
        static const char sStructRevXmlRevXYZEnd[]= "</xyz>";
        static const char sStructAuxXml[]         = "<identifier.auxiliary-info";
        static const char sStructAuxXmlEnd[]      = "</identifier.auxiliary-info";
        int         bInTheAuxInfo           = 0;

        while ( 0 < (res = my_fgets( szLine, sizeof(szLine)-1, inp_molfile, &bTooLongLine ) ) ) {
        
            /********************* find and interpret structure header ************/
            if ( !memcmp(szLine, sStructHdrXml, sizeof(sStructHdrXml)-1) ) {
                num_struct = 1;
                p = szLine + sizeof(sStructHdrXml)-1;
                longID = 0;
                num_atoms = 0;
                /* structure number */
                if ( q = strstr( p, sStructHdrXmlNumber ) ) {
                    p = q + sizeof(sStructHdrXmlNumber)-1;
                    longID = strtol( p, &q, 10);
                    if ( q && *q == '\"' )
                        p = q+1;
                }
                if ( pSdfLabel ) {
                    pSdfLabel[0] = '\0';
                }
                if ( pSdfValue ) {
                    pSdfValue[0] = '\0';
                }
                /* pSdfLabel */
                if ( q = strstr( p, sStructHdrXmlIdName ) ) {
                    p = q + sizeof(sStructHdrXmlIdName)-1;
                    q = strchr( p, '\"' );
                    if ( q ) {
                        len = inchi_min( q-p+1, MAX_SDF_HEADER-1);
                        if ( pSdfLabel ) {
                            mystrncpy( pSdfLabel, p, len );
                        }
                        p = q+1;
                    }
                }
                /* pSdfValue */
                if ( q = strstr( p, sStructHdrXmlIdValue ) ) {
                    p = q + sizeof(sStructHdrXmlIdValue)-1;
                    q = strchr( p, '\"' );
                    if ( q ) {
                        len = inchi_min( q-p+1, MAX_SDF_VALUE-1);
                        if ( pSdfValue ) {
                            mystrncpy( pSdfValue, p, len );
                        }
                        p = q+1;
                    }
                }
                if ( Id )
                    *Id = longID;
                bHeaderRead = 1;
                bErrorMsg = bRestoreInfo = 0;
            } else
            if ( bHeaderRead && (bFatal=0, len=sizeof(sStructMsgXmlErr)-1,      !memcmp(szLine, sStructMsgXmlErr, len)) ||
                 bHeaderRead && (len=sizeof(sStructMsgXmlErrFatal)-1, !memcmp(szLine, sStructMsgXmlErrFatal, len))&&(bFatal=1)) {
                p = szLine+len;
                q = strchr( p, '\"' );
                if ( q && !bFindNext ) {
                    int c;
                    bErrorMsg = 1;
                    pStrErr[0] = '\0';
                    c = *q;
                    *q = '\0';
                    MOLFILE_ERR_SET (*err, 0, p);
                    *q = c;
                }
                *err      = bFatal? INCHI_INP_FATAL_ERR : INCHI_INP_ERROR_ERR;
                num_atoms = bFatal? INCHI_INP_FATAL_RET : INCHI_INP_ERROR_RET;
                goto bypass_end_of_INChI;
            } else
            if ( bHeaderRead && !memcmp(szLine, sStructAuxXml, sizeof(sStructAuxXml)-1) ) {
                bInTheAuxInfo = 1;
            } else
            if ( bHeaderRead && !memcmp(szLine, sStructAuxXmlEnd, sizeof(sStructAuxXmlEnd)-1) ) {
                *err      = INCHI_INP_ERROR_ERR;
                num_atoms = INCHI_INP_ERROR_RET;
                MOLFILE_ERR_SET (*err, 0, "Missing reversibility info" );
                goto bypass_end_of_INChI; /* reversibility info not found */
            } else
            if ( bHeaderRead && bInTheAuxInfo && !memcmp(szLine, sStructRevXmlRevHdr, sizeof(sStructRevXmlRevHdr)-1) ) {
                /***********************  atoms xml ***************************/
                num_struct = 1;
                res = my_fgets( szLine, sizeof(szLine)-1, inp_molfile, &bTooLongLine );
                if ( res <= 0 ) {
                    num_atoms = INCHI_INP_EOF_RET; /* no data, probably end of file */
                    *err      = INCHI_INP_EOF_ERR;
                    goto bypass_end_of_INChI;
                }
                if ( memcmp(szLine, sStructRevXmlRevAt, sizeof(sStructRevXmlRevAt)-1) ) {
                    bHeaderRead = 0; /* invalid reversibility info; look for another header */
                    continue;
                }
                /* read (the head of) the atoms line */
                res = my_fgets( szLine, sizeof(szLine)-1, inp_molfile, &bTooLongLine );
                if ( res <= 0 ) {
                    num_atoms = INCHI_INP_EOF_RET; /* no data */
                    *err      = INCHI_INP_EOF_ERR;
                    goto bypass_end_of_INChI;
                }
                p = szLine;
                num_atoms = strtol( p, &q, 10 );
                if ( !num_atoms || !q || !*q ) {
                    num_atoms = INCHI_INP_EOF_RET; /* no atom data */
                    *err      = INCHI_INP_EOF_ERR;
                    goto bypass_end_of_INChI;
                }
                p = q;
                /* Molfile chirality flag */
                switch( *p ) {
                case 'c':
                    InpAtomFlags |= FLAG_INP_AT_CHIRAL;
                    p ++;
                    break;
                case 'n':
                    InpAtomFlags |= FLAG_INP_AT_NONCHIRAL;
                    p ++;
                    break;
                }
                if ( at && *at ) {
                    if ( num_atoms > max_num_at ) {
                        inchi_free( *at );
                        *at = NULL;
                    } else {
                        memset( *at, 0, max_num_at * sizeof( **at ) );
                        atom = *at;
                    }
                }
                if ( !at || !*at ) {
                    atom = Create_Atom( num_atoms+1 );
                    if ( !atom ) {
                        num_atoms = INCHI_INP_FATAL_RET; /* fatal error: cannot allocate */
                        *err      = INCHI_INP_FATAL_ERR;
                        MOLFILE_ERR_SET (*err, 0, "Out of RAM");
                        goto bypass_end_of_INChI;
                    }
                }
                if ( stereo0D && *stereo0D ) {
                    if ( num_atoms > max_len_stereo0D ) {
                        FreeInchi_Stereo0D( stereo0D );
                    } else {
                        memset( *stereo0D, 0, max_len_stereo0D * sizeof( **stereo0D ) );
                        atom_stereo0D = *stereo0D;
                    }
                }
                if ( !stereo0D || !*stereo0D ) {
                    max_len_stereo0D = num_atoms+1;
                    atom_stereo0D = CreateInchi_Stereo0D( max_len_stereo0D );
                    if ( !atom_stereo0D ) {
                        num_atoms = INCHI_INP_FATAL_RET; /* fatal error: cannot allocate */
                        *err      = INCHI_INP_FATAL_ERR;
                        MOLFILE_ERR_SET (*err, 0, "Out of RAM");
                        goto bypass_end_of_INChI;
                    }
                }

                i = 0;
                bItemIsOver = 0;
                res2 = bTooLongLine2 = -1;
                
                /* read all atoms xml */
                while ( i < num_atoms ) {
                    pos = p - szLine;
                    if ( !bItemIsOver && (int)sizeof(szLine)-res + pos > (int)sizeof(szNextLine) ) {
                        /* load next line if possible */
                        res2 = my_fgets( szNextLine, sizeof(szNextLine)-1, inp_molfile, &bTooLongLine2 );
                        if ( res2 > 0 && memcmp(szNextLine, sStructRevXmlRevAtEnd, sizeof(sStructRevXmlRevAtEnd)-1) ) {
                            if ( pos ) {
                                res -= pos;  /* number of chars left to process in szLine */
                                memmove( szLine, p, res*sizeof(szLine[0]) ); /* move them to the start of the line */
                            }
                            memcpy( szLine+res, szNextLine, (res2+1)*sizeof(szNextLine[0]) );
                            res += res2;
                            szLine[res] = '\0';
                            bTooLongLine = bTooLongLine2;
                            p = szLine;
                        } else {
                            bItemIsOver = 1;
                        }
                    }
                    /* element, first char */
                    if ( !isalpha( UCINT *p ) || !isupper( UCINT *p ) || i >= num_atoms ) {
                        bHeaderRead = 0; /* wrong atom data */
                        num_atoms = INCHI_INP_ERROR_RET; /* was 0, error */
                        *err      = INCHI_INP_ERROR_ERR;     /* 40 */
                        MOLFILE_ERR_SET (*err, 0, "Wrong atoms data");
                        goto bypass_end_of_INChI;
                    }
                    atom[i].elname[0] = *p ++;
                    /* element, second char */
                    if ( isalpha( UCINT *p ) && islower( UCINT *p ) ) {
                        atom[i].elname[1] = *p ++;
                    }
#if( defined(INCHI_LIBRARY) || defined(INCHI_MAIN) )
#else
                    atom[i].el_number = get_periodic_table_number( atom[i].elname );
#endif
                    /* bonds' valence */
                    if ( isdigit( UCINT *p ) ) {
                        AT_BONDS_VAL(atom,i) = (char)strtol( p, &q, 10 );
                        if ( !AT_BONDS_VAL(atom,i) )
                            AT_BONDS_VAL(atom,i) = ISOLATED_ATOM; /* same convention as in MOLfile, found zero bonds valence */
                        p = q;
                    }
                    /* charge */
                    atom[i].charge = (*p == '+')? 1 : (*p == '-')? -1 : 0;
                    if ( atom[i].charge ) {
                        p ++;
                        if ( isdigit( UCINT *p ) ) {
                            atom[i].charge *= (S_CHAR)(strtol( p, &q, 10 ) & CHAR_MASK);
                            p = q;
                        }
                    }
                    /* radical */
                    if ( *p == '.' ) {
                        p ++;
                        if ( isdigit( UCINT *p ) ) {
                            atom[i].radical = (S_CHAR)strtol( p, &q, 10 );
                            p = q;
                        }
                    }
                    /* isotopic mass */
                    if ( *p == 'i' ) {
                        p ++;
                        if ( isdigit( UCINT *p ) ) {
                            int mw = strtol( p, &q, 10 );
                            p = q;
#if( defined(INCHI_LIBRARY) || defined(INCHI_MAIN) )
                            atom[i].isotopic_mass = mw;
#else
                            mw -= get_atw_from_elnum( atom[i].el_number );
                            if ( mw >= 0 )
                                mw ++;
                            atom[i].iso_atw_diff = mw;
#endif
                        }
                    }
                    /* parity */
                    switch( *p ) {
                    case 'o':
                        parity = INCHI_PARITY_ODD;
                        p ++;
                        break;
                    case 'e':
                        parity = INCHI_PARITY_EVEN;
                        p ++;
                        break;
                    case 'u':
                        parity = INCHI_PARITY_UNKNOWN;
                        p ++;
                        break;
                    case '?':
                        parity = INCHI_PARITY_UNDEFINED;
                        p ++;
                        break;
                    default:
                        parity = 0;
                        break;
                    }
                    if ( parity ) {
                        atom_stereo0D[len_stereo0D].central_atom = i;
                        atom_stereo0D[len_stereo0D].parity       = parity;
                        atom_stereo0D[len_stereo0D].type         = INCHI_StereoType_Tetrahedral;
                        len_stereo0D ++;
                    }
                    /* isotopic h, d, t */
                    for ( k = 0; k < NUM_H_ISOTOPES; k ++ ) {
                        if ( *p == szIsoH[k] ) {
                            NUM_ISO_Hk(atom,i,k) = 1;
                            p ++;
                            if ( isdigit( UCINT *p ) ) {
                                NUM_ISO_Hk(atom,i,k) = (char)strtol( p, &q, 10 );
                                p = q;
                            }
                        }
                    }
                    i ++;
                }
                if ( !bItemIsOver || p - szLine != res || i != num_atoms ) {
                    num_atoms = INCHI_INP_ERROR_RET; /* error */
                    *err      = INCHI_INP_ERROR_ERR;
                    MOLFILE_ERR_SET (*err, 0, "Wrong number of atoms");
                    goto bypass_end_of_INChI;
                }
                /********************** bonds xml ****************************/
                res = my_fgets( szLine, sizeof(szLine)-1, inp_molfile, &bTooLongLine );
                if ( res <= 0 ) {
                    num_atoms = 0; /* no data */
                    goto bypass_end_of_INChI;
                }
                if ( memcmp(szLine, sStructRevXmlRevBn, sizeof(sStructRevXmlRevBn)-1) ) {
                    bHeaderRead = 0; /* invalid reversibility info; look for another header */
                    continue;
                }
                /* read (the head of) the xml bonds line */
                res = my_fgets( szLine, sizeof(szLine)-1, inp_molfile, &bTooLongLine );
                if ( res <= 0 ) {
                    num_atoms = INCHI_INP_ERROR_RET; /* was 0; error: no data -- eof? */
                    *err      = INCHI_INP_ERROR_ERR;
                    goto bypass_end_of_INChI;
                }
                i = 1;
                bItemIsOver = 0;
                res2 = bTooLongLine2 = -1;
                p = szLine;
                if ( !memcmp(szLine, sStructRevXmlRevBnEnd, sizeof(sStructRevXmlRevBnEnd)-1) ) {
                    /* empty bonds section */
                    res = 0;
                    bItemIsOver = 1;
                }
                /* read all bonds (xml), starting from atom 1 (not 0) */
                while ( i < num_atoms ) {
                    pos = p - szLine;
                    if ( !bItemIsOver && (int)sizeof(szLine)-res + pos > (int)sizeof(szNextLine) ) {
                        /* load next line if possible */
                        res2 = my_fgets( szNextLine, sizeof(szNextLine)-1, inp_molfile, &bTooLongLine2 );
                        if ( res2 > 0 && memcmp(szNextLine, sStructRevXmlRevBnEnd, sizeof(sStructRevXmlRevBnEnd)-1) ) {
                            if ( pos ) {
                                res -= pos;  /* number of chars left to process in szLine */
                                memmove( szLine, p, res*sizeof(szLine[0]) ); /* move them to the start of the line */
                            }
                            memcpy( szLine+res, szNextLine, (res2+1)*sizeof(szNextLine[0]) );
                            res += res2;
                            szLine[res] = '\0';
                            bTooLongLine = bTooLongLine2;
                            p = szLine;
                        } else {
                            bItemIsOver = 1;
                        }
                    }
                    if ( i >= num_atoms ) {
                        break;
                    }
                    /* bond, first char */
                    if ( *p == ';' ) {
                        p ++;
                        i ++;
                        continue;
                    }
                    if ( !isalpha( UCINT *p ) ) {
                        num_atoms = INCHI_INP_ERROR_RET; /* error in input data */
                        *err      = INCHI_INP_ERROR_ERR;
                        MOLFILE_ERR_SET (*err, 0, "Wrong bonds data");
                        goto bypass_end_of_INChI;
                    }
                    bond_char = *p ++;
                    /* bond parity */
                    switch( *p ) {
                    case '-':
                        bond_parity = INCHI_PARITY_ODD;
                        p ++;
                        break;
                    case '+':
                        bond_parity = INCHI_PARITY_EVEN;
                        p ++;
                        break;
                    case 'u':
                        bond_parity = INCHI_PARITY_UNKNOWN;
                        p ++;
                        break;
                    case '?':
                        bond_parity = INCHI_PARITY_UNDEFINED;
                        p ++;
                        break;
                    default:
                        bond_parity = 0;
                        break;
                    }
                    if ( bond_parity ) {
                        switch( *p ) {
                        case '-':
                            bond_parityNM = INCHI_PARITY_ODD;
                            p ++;
                            break;
                        case '+':
                            bond_parityNM = INCHI_PARITY_EVEN;
                            p ++;
                            break;
                        case 'u':
                            bond_parityNM = INCHI_PARITY_UNKNOWN;
                            p ++;
                            break;
                        case '?':
                            bond_parityNM = INCHI_PARITY_UNDEFINED;
                            p ++;
                            break;
                        default:
                            bond_parityNM = 0;
                            break;
                        }
                    } else {
                        bond_parityNM = 0;
                    }

                    /* neighbor of the current atom */
                    if ( !isdigit( UCINT *p ) ) {
                        num_atoms = INCHI_INP_ERROR_RET; /* error in input data */
                        *err      = INCHI_INP_ERROR_ERR;
                        MOLFILE_ERR_SET (*err, 0, "Wrong bonds data");
                        goto bypass_end_of_INChI;
                    }
                    neigh = (int)strtol( p, &q, 10 )-1;

                    if ( i >= num_atoms || neigh >= num_atoms ) {
                        num_atoms = INCHI_INP_ERROR_RET; /* error in input data */
                        *err      = INCHI_INP_ERROR_ERR;
                        MOLFILE_ERR_SET (*err, 0, "Bond to nonexistent atom");
                        goto bypass_end_of_INChI;
                    }
                    p = q;
                    bond_stereo1 = bond_stereo2 = 0;

                    /* bond type & 2D stereo */
                    switch( bond_char ) {
                    case 'v':
                        bond_type    = INCHI_BOND_TYPE_SINGLE;
                        bond_stereo1 = INCHI_BOND_STEREO_SINGLE_1EITHER;
                        bond_stereo2 = INCHI_BOND_STEREO_SINGLE_2EITHER;
                        break;
                    case 'V':
                        bond_type    = INCHI_BOND_TYPE_SINGLE;
                        bond_stereo1 = INCHI_BOND_STEREO_SINGLE_2EITHER;
                        bond_stereo2 = INCHI_BOND_STEREO_SINGLE_1EITHER;
                        break;
                    case 'w':
                        bond_type    = INCHI_BOND_TYPE_DOUBLE;
                        bond_stereo1 = 
                        bond_stereo2 = INCHI_BOND_STEREO_DOUBLE_EITHER;
                        break;
                    case 's':
                        bond_type    = INCHI_BOND_TYPE_SINGLE;
                        break;
                    case 'd':
                        bond_type    = INCHI_BOND_TYPE_DOUBLE;
                        break;
                    case 't':
                        bond_type    = INCHI_BOND_TYPE_TRIPLE;
                        break;
                    case 'a':
                        bond_type    = INCHI_BOND_TYPE_ALTERN;
                        break;
                    case 'p':
                        bond_type    =  INCHI_BOND_TYPE_SINGLE;
                        bond_stereo1 =  INCHI_BOND_STEREO_SINGLE_1UP;
                        bond_stereo2 =  INCHI_BOND_STEREO_SINGLE_2UP;
                        break;
                    case 'P':
                        bond_type    =  INCHI_BOND_TYPE_SINGLE;
                        bond_stereo1 =  INCHI_BOND_STEREO_SINGLE_2UP;
                        bond_stereo2 =  INCHI_BOND_STEREO_SINGLE_1UP;
                        break;
                    case 'n':
                        bond_type    =  INCHI_BOND_TYPE_SINGLE;
                        bond_stereo1 =  INCHI_BOND_STEREO_SINGLE_1DOWN;
                        bond_stereo2 =  INCHI_BOND_STEREO_SINGLE_2DOWN;
                        break;
                    case 'N':
                        bond_type    =  INCHI_BOND_TYPE_SINGLE;
                        bond_stereo1 =  INCHI_BOND_STEREO_SINGLE_2DOWN;
                        bond_stereo2 =  INCHI_BOND_STEREO_SINGLE_1DOWN;
                        break;
                    default:
                        num_atoms = INCHI_INP_ERROR_RET; /* error */
                        *err      = INCHI_INP_ERROR_ERR;
                        MOLFILE_ERR_SET (*err, 0, "Wrong bond type");
                        goto bypass_end_of_INChI;
                    }
                    k = AT_NUM_BONDS(atom[i]) ++;
                    atom[i].bond_type[k]   = bond_type;
                    atom[i].bond_stereo[k] = bond_stereo1;
                    atom[i].neighbor[k]    = (ATOM_NUMBER)neigh;
                    k2 = AT_NUM_BONDS(atom[neigh]) ++;
                    atom[neigh].bond_type[k2]   = bond_type;
                    atom[neigh].bond_stereo[k2] = bond_stereo2;
                    atom[neigh].neighbor[k2]    = (ATOM_NUMBER)i;
                    bond_parity |= (bond_parityNM << SB_PARITY_SHFT);

                    if ( bond_parity ) {
                        if ( max_len_stereo0D <= len_stereo0D ) {
                            /* realloc atom_Stereo0D */
                            inchi_Stereo0D *new_atom_stereo0D = CreateInchi_Stereo0D( max_len_stereo0D+num_atoms );
                            if ( !new_atom_stereo0D ) {
                                num_atoms = INCHI_INP_FATAL_RET; /* fatal error: cannot allocate */
                                *err      = INCHI_INP_FATAL_ERR;
                                MOLFILE_ERR_SET (*err, 0, "Out of RAM");
                                goto bypass_end_of_INChI;
                            }
                            memcpy( new_atom_stereo0D, atom_stereo0D, len_stereo0D * sizeof(*atom_stereo0D) );
                            FreeInchi_Stereo0D( &atom_stereo0D );
                            atom_stereo0D = new_atom_stereo0D;
                            max_len_stereo0D += num_atoms;
                        }
                        /* (a) i may be allene endpoint and     neigh = allene middle point or
                           (b) i may be allene middle point and neigh = allene endpoint
                           !!!!! CURRENTLY ONLY (b) IS ALLOWED !!!!!
                        */
                        atom_stereo0D[len_stereo0D].neighbor[1] = neigh; /* neigh < i */
                        atom_stereo0D[len_stereo0D].neighbor[2] = i;
                        atom_stereo0D[len_stereo0D].parity      = bond_parity;
                        atom_stereo0D[len_stereo0D].type        = INCHI_StereoType_DoubleBond; /* incl allenes & cumulenes */
                        len_stereo0D ++;
                    }
                }
                if ( !bItemIsOver || p - szLine != res || i != num_atoms ) {
                    num_atoms = INCHI_INP_ERROR_RET; /* error in input data */
                    *err      = INCHI_INP_ERROR_ERR;
                    MOLFILE_ERR_SET (*err, 0, "Wrong number of bonds");
                    goto bypass_end_of_INChI;
                }
                /********************** coordinates xml ****************************/
                if ( pszCoord = (MOL_COORD*)inchi_malloc(inchi_max(num_atoms,1) * sizeof(MOL_COORD)) ) {
                    memset( pszCoord, ' ', inchi_max(num_atoms,1) * sizeof(MOL_COORD));
                    res = my_fgets( szLine, sizeof(szLine)-1, inp_molfile, &bTooLongLine );
                    if ( res <= 0  ||
                         /* compare the header */
                         memcmp(szLine, sStructRevXmlRevXYZ, sizeof(sStructRevXmlRevXYZ)-1) ||
                         /* read (the head of) the coordinates (xml) line */
                         0 >= (res = my_fgets( szLine, sizeof(szLine)-1, inp_molfile, &bTooLongLine ))) {
                        num_atoms = INCHI_INP_ERROR_RET; /* error in input data: atoms, bonds & coord must be present together */
                        *err      = INCHI_INP_ERROR_ERR;
                        MOLFILE_ERR_SET (*err, 0, "Missing atom coordinates data");
                        goto bypass_end_of_INChI;
                    }
                    i = 0;
                    bItemIsOver = 0;
                    res2 = bTooLongLine2 = -1;
                    p = szLine;
                    if ( !memcmp(szLine, sStructRevXmlRevXYZEnd, sizeof(sStructRevXmlRevXYZEnd)-1) ) {
                        /* empty bonds section */
                        res = 0;
                        bItemIsOver = 1;
                    }
                    /* read all coordinates (xml), starting from atom 1 (not 0) */
                    while ( i < num_atoms ) {
                        pos = p - szLine;
                        if ( !bItemIsOver && (int)sizeof(szLine)-res + pos > (int)sizeof(szNextLine) ) {
                            /* load next line if possible */
                            res2 = my_fgets( szNextLine, sizeof(szNextLine)-1, inp_molfile, &bTooLongLine2 );
                            if ( res2 > 0 && memcmp(szNextLine, sStructRevXmlRevXYZEnd, sizeof(sStructRevXmlRevXYZEnd)-1) ) {
                                if ( pos ) {
                                    res -= pos;  /* number of chars left to process in szLine */
                                    memmove( szLine, p, res*sizeof(szLine[0]) ); /* move them to the start of the line */
                                }
                                memcpy( szLine+res, szNextLine, (res2+1)*sizeof(szNextLine[0]) );
                                res += res2;
                                szLine[res] = '\0';
                                bTooLongLine = bTooLongLine2;
                                p = szLine;
                            } else {
                                bItemIsOver = 1;
                            }
                        }
                        /* coord, first char */
                        if ( *p == ';' ) {
                            for ( k = 0; k < NUM_COORD; k ++ ) {
                                pszCoord[i][LEN_COORD*k + 4] = '0';
                            }
                            p ++;
                            i ++;
                            continue;
                        }
                        for ( k = 0; k < 3; k ++ ) {
                            double xyz;
                            bNonZeroXYZ = 0;
                            if ( *p == ';' ) {
                                pszCoord[i][LEN_COORD*k + 4] = '0';
                                xyz = 0.0;
                            } else
                            if ( *p == ',' ) {
                                /* empty */
                                pszCoord[i][LEN_COORD*k + 4] = '0';
                                xyz = 0.0;
                                p ++;
                            } else {
                                xyz = strtod( p, &q );
                                bNonZeroXYZ = fabs(xyz) > MIN_BOND_LENGTH;
                                if ( q != NULL ) {
                                    memcpy( pszCoord[i]+LEN_COORD*k, p, q-p );
                                    if ( *q == ',' )
                                        q ++;
                                    p = q;
                                } else {
                                    pszCoord[i][LEN_COORD*k + 4] = '0';
                                }
                            }
                            switch( k ) {
                            case 0:
                                atom[i].x = xyz;
                                b2D |= bNonZeroXYZ;
                                break;
                            case 1:
                                atom[i].y = xyz;
                                b2D |= bNonZeroXYZ;
                                break;
                            case 2:
                                b3D |= bNonZeroXYZ;
                                atom[i].z = xyz;
                                break;
                            }
                        }
                        if ( *p == ';' ) {
                            p ++;  /* end of this triple of coordinates */
                            i ++;
                        } else {
                            num_atoms = INCHI_INP_ERROR_RET; /* error in input data: atoms, bonds & coord must be present together */
                            *err      = INCHI_INP_ERROR_ERR;
                            MOLFILE_ERR_SET (*err, 0, "Wrong atom coordinates data");
                            goto bypass_end_of_INChI;
                        }
                    }
                    if ( !bItemIsOver || p - szLine != res || i != num_atoms ) {
                        num_atoms = INCHI_INP_ERROR_RET; /* error in input data: atoms, bonds & coord must be present together */
                        *err      = INCHI_INP_ERROR_ERR;
                        MOLFILE_ERR_SET (*err, 0, "Wrong number of coordinates");
                        goto bypass_end_of_INChI;
                    }
                } else { /* allocation failed */
                    num_atoms = INCHI_INP_FATAL_RET; 
                    *err      = INCHI_INP_FATAL_ERR;
                    MOLFILE_ERR_SET (*err, 0, "Out of RAM");
                    goto bypass_end_of_INChI;
                }
                /* set special valences and implicit H (xml) */
                b23D = b2D | b3D;
                b2D = b3D = 0;
                if ( at ) {
                    if ( !*at ) {
                        int a1, a2, n1, n2, valence;
                        int chem_bonds_valence;
                        int    nX=0, nY=0, nZ=0, nXYZ;
                        *at = atom;
                        /* special valences */
                        for ( bNonMetal = 0; bNonMetal < 1 /*2*/; bNonMetal ++ ) {
                            for ( a1 = 0; a1 < num_atoms; a1 ++ ) {
                                int num_bond_type[MAX_INPUT_BOND_TYPE - MIN_INPUT_BOND_TYPE + 1];
#if( defined(INCHI_LIBRARY) || defined(INCHI_MAIN) )
#else
                                int bHasMetalNeighbor=0;
#endif
                                memset( num_bond_type, 0, sizeof(num_bond_type) );

                                valence = AT_BONDS_VAL(atom, a1); /*  save atom valence if available */
                                AT_BONDS_VAL(atom, a1) = 0;
#if( defined(INCHI_LIBRARY) || defined(INCHI_MAIN) )
#else
                                atom[a1].orig_at_number = a1+1;
#endif
                                nX = nY = nZ = 0;
                                for ( n1 = 0; n1 < AT_NUM_BONDS(atom[a1]); n1 ++ ) {
                                    bond_type = atom[a1].bond_type[n1] - MIN_INPUT_BOND_TYPE;
                                    if (  bond_type < 0 || bond_type > MAX_INPUT_BOND_TYPE - MIN_INPUT_BOND_TYPE ) {
                                        bond_type = 0; /* cannot happen */
                                        MOLFILE_ERR_SET (*err, 0, "Unknown bond type in InChI aux assigned as a single bond");
                                    }

                                    num_bond_type[ bond_type ] ++;
                                    nNumBonds ++;
                                    if ( b23D ) {
                                        neigh = atom[a1].neighbor[n1];
                                        nX |= (fabs(atom[a1].x - atom[neigh].x) > MIN_BOND_LENGTH);
                                        nY |= (fabs(atom[a1].y - atom[neigh].y) > MIN_BOND_LENGTH);
                                        nZ |= (fabs(atom[a1].z - atom[neigh].z) > MIN_BOND_LENGTH);
                                    }
                                }
                                chem_bonds_valence = 0;
                                for ( n1 = 0; MIN_INPUT_BOND_TYPE + n1 <= 3 && MIN_INPUT_BOND_TYPE + n1 <= MAX_INPUT_BOND_TYPE; n1 ++ ) {
                                    chem_bonds_valence += (MIN_INPUT_BOND_TYPE + n1) * num_bond_type[n1];
                                }
                                if ( MIN_INPUT_BOND_TYPE <= INCHI_BOND_TYPE_ALTERN && INCHI_BOND_TYPE_ALTERN <= MAX_INPUT_BOND_TYPE &&
                                     ( n2 = num_bond_type[INCHI_BOND_TYPE_ALTERN-MIN_INPUT_BOND_TYPE] ) ) {
                                    /* accept input aromatic bonds for now */
                                    switch ( n2 ) {
                                    case 2:
                                        chem_bonds_valence += 3;  /* =A- */
                                        break;
                                    case 3:
                                        chem_bonds_valence += 4;  /* =A< */
                                        break;
                                    default:
                                        /*  if 1 or >= 4 aromatic bonds then replace such bonds with single bonds */
                                        for ( n1 = 0; n1 < AT_NUM_BONDS(atom[a1]); n1 ++ ) {
                                            if ( atom[a1].bond_type[n1] == INCHI_BOND_TYPE_ALTERN ) {
                                                ATOM_NUMBER *p1;
                                                a2 = atom[a1].neighbor[n1];
                                                p1 = IN_NEIGH_LIST( atom[a2].neighbor, (ATOM_NUMBER)a1, AT_NUM_BONDS(atom[a2]) );
                                                if ( p1 ) {
                                                    atom[a1].bond_type[n1] = 
                                                    atom[a2].bond_type[p1-atom[a2].neighbor] = INCHI_BOND_TYPE_SINGLE;
                                                } else {
                                                    *err = -2;  /*  Program error */
                                                    MOLFILE_ERR_SET (*err, 0, "Program error interpreting InChI aux");
                                                    num_atoms = 0;
                                                    goto bypass_end_of_INChI; /*  no structure */
                                                }
                                            }
                                        }
                                        chem_bonds_valence += n2;
                                        *err |= 32; /*  Unrecognized aromatic bond(s) replaced with single */
                                        MOLFILE_ERR_SET (*err, 0, "Atom has 1 or more than 3 aromatic bonds");
                                        break;
                                    }
                                }

#if( defined(INCHI_LIBRARY) || defined(INCHI_MAIN) )
                                /*************************************************************************************
                                 *
                                 *  Set number of hydrogen atoms
                                 */
                                {
                                    int num_iso_H;
                                    num_iso_H = atom[a1].num_iso_H[1] + atom[a1].num_iso_H[2] + atom[a1].num_iso_H[3];
                                    if ( valence == ISOLATED_ATOM ) {
                                        atom[a1].num_iso_H[0] = 0;
                                    } else
                                    if ( valence && valence >= chem_bonds_valence ) {
                                        atom[a1].num_iso_H[0] = valence - chem_bonds_valence;
                                    } else
                                    if ( valence || bDoNotAddH ) {
                                        atom[a1].num_iso_H[0] = 0;
                                    } else
                                    if ( !bDoNotAddH ) {
                                        atom[a1].num_iso_H[0] = -1; /* auto add H */
                                    }
                                }
#else                                
                                /* added 2006-07-19 to process aromatic bonds same way as from molfile */
                                if ( n2 && !valence ) {
                                    int num_H = NUMH(atom, a1); /* only isotopic */
                                    int chem_valence = chem_bonds_valence;
                                    int bUnusualValenceArom = 
                                        detect_unusual_el_valence( (int)atom[a1].el_number, atom[a1].charge,
                                                                    atom[a1].radical, chem_valence,
                                                                    num_H, atom[a1].valence );
                                    int bUnusualValenceNoArom = 
                                        detect_unusual_el_valence( (int)atom[a1].el_number, atom[a1].charge,
                                                                    atom[a1].radical, chem_valence-1,
                                                                    num_H, atom[a1].valence );
#if ( CHECK_AROMBOND2ALT == 1 )
                                    if ( bUnusualValenceArom && !bUnusualValenceNoArom && 0 == nBondsValToMetal( atom, a1) )
#else
                                    if ( bUnusualValenceArom && !bUnusualValenceNoArom )
#endif                     
                                    {
                                        /* typically NH in 5-member aromatic ring */
                                        chem_bonds_valence --;
                                    }
                                } else
                                if ( n2 && valence ) {
                                    /* atom has aromatic bonds AND the chemical valence is known */
                                    int num_H = NUMH(atom, a1);
                                    int chem_valence = chem_bonds_valence + num_H;
                                    if ( valence == chem_valence-1 ) {
                                        /* typically NH in 5-member aromatic ring */
                                        chem_bonds_valence --;
                                    }
                                }

                                atom[a1].chem_bonds_valence = chem_bonds_valence;
                                atom[a1].num_H = get_num_H( atom[a1].elname, atom[a1].num_H, atom[a1].num_iso_H, atom[a1].charge, atom[a1].radical,
                                                          atom[a1].chem_bonds_valence,
                                                          valence,
                                                          0, bDoNotAddH, bHasMetalNeighbor );
#endif
                            }
                        }
                        nNumBonds /= 2;
                        if ( b23D && nNumBonds ) {
                            nXYZ = nX+nY+nZ;
                            b2D  = (nXYZ > 0);
                            b3D  = (nXYZ == 3);
                            *num_dimensions = b3D? 3 : b2D? 2 : 0;
                            *num_bonds = nNumBonds;
                        }
                        /*======= 0D parities =================================*/
#if( defined(INCHI_LIBRARY) || defined(INCHI_MAIN) )
                        if ( len_stereo0D > 0 && atom_stereo0D && stereo0D ) {
                            *stereo0D     = atom_stereo0D;
                            *num_stereo0D = len_stereo0D;
                        } else {
                            FreeInchi_Stereo0D( &atom_stereo0D );
                            *num_stereo0D = len_stereo0D = 0;
                        }
#endif
                        for ( i = 0; i < len_stereo0D; i ++ ) {
                            ATOM_NUMBER *p1, *p2;
                            int     sb_ord_from_a1 = -1, sb_ord_from_a2 = -1, bEnd1 = 0, bEnd2 = 0;
                            switch( atom_stereo0D[i].type ) {

                            case INCHI_StereoType_Tetrahedral:
                                a1 = atom_stereo0D[i].central_atom;
                                if ( atom_stereo0D[i].parity && (AT_NUM_BONDS(atom[a1]) == 3 || AT_NUM_BONDS(atom[a1]) == 4) ) {
                                    int ii, kk = 0;
                                    if ( AT_NUM_BONDS(atom[a1]) == 3 ) {
                                        atom_stereo0D[i].neighbor[kk++] = a1;
                                    }
                                    for ( ii = 0; ii < AT_NUM_BONDS(atom[a1]); ii ++ ) {
                                        atom_stereo0D[i].neighbor[kk++] = atom[a1].neighbor[ii];
                                    }
                                }
                            
                            break;

                            case INCHI_StereoType_DoubleBond:
#define MAX_CHAIN_LEN 20
                                a1 = atom_stereo0D[i].neighbor[1];
                                a2 = atom_stereo0D[i].neighbor[2];
                                p1 = IN_NEIGH_LIST( atom[a1].neighbor, (ATOM_NUMBER)a2, AT_NUM_BONDS(atom[a1]) );
                                p2 = IN_NEIGH_LIST( atom[a2].neighbor, (ATOM_NUMBER)a1, AT_NUM_BONDS(atom[a2]) );
                                if ( !p1 || !p2 ) {
                                    atom_stereo0D[i].type = INCHI_StereoType_None;
                                    atom_stereo0D[i].central_atom = NO_ATOM;
                                    atom_stereo0D[i].neighbor[0] =
                                    atom_stereo0D[i].neighbor[3] = -1;
                                    *err |= 64; /* Error in cumulene stereo */
                                    MOLFILE_ERR_SET (*err, 0, "0D stereobond not recognized");
                                    break;
                                }
                                /* streobond, allene, or cumulene */

                                sb_ord_from_a1 = p1 - atom[a1].neighbor;
                                sb_ord_from_a2 = p2 - atom[a2].neighbor;
                                
                                if (  AT_NUM_BONDS(atom[a1]) == 2 &&
                                      atom[a1].bond_type[0] + atom[a1].bond_type[1] == 2*INCHI_BOND_TYPE_DOUBLE &&
                                      0 == inchi_NUMH2(atom, a1) &&
                                     (AT_NUM_BONDS(atom[a2]) != 2 ||
                                      atom[a2].bond_type[0] + atom[a2].bond_type[1] != 2*INCHI_BOND_TYPE_DOUBLE ) ) {
                                    bEnd2 = 1; /* a2 is the end-atom, a1 is middle atom */   
                                }
                                if (  AT_NUM_BONDS(atom[a2]) == 2 &&
                                      atom[a2].bond_type[0] + atom[a2].bond_type[1] == 2*INCHI_BOND_TYPE_DOUBLE &&
                                      0 == inchi_NUMH2(atom, a2) &&
                                     (AT_NUM_BONDS(atom[a1]) != 2 ||
                                      atom[a1].bond_type[0] + atom[a1].bond_type[1] != 2*INCHI_BOND_TYPE_DOUBLE ) ) {
                                    bEnd1 = 1; /* a1 is the end-atom, a2 is middle atom */   
                                }
                                
                                if ( bEnd2 + bEnd1 == 1 ) {
                                    /* allene or cumulene */
                                    ATOM_NUMBER  chain[MAX_CHAIN_LEN+1], prev, cur, next;
                                    if ( bEnd2 && !bEnd1 ) {
                                        cur = a1;
                                        a1 = a2;
                                        a2 = cur;
                                        sb_ord_from_a1 = sb_ord_from_a2;
                                    }
                                    sb_ord_from_a2 = -1;
                                    cur  = a1;
                                    next = a2;
                                    len = 0;
                                    chain[len++] = cur;
                                    chain[len++] = next;
                                    while ( len < MAX_CHAIN_LEN ) { /* arbitrary very high upper limit to prevent infinite loop */
                                        prev = cur;
                                        cur  = next;
                                            /* follow double bond path && avoid going back */
                                        if ( AT_NUM_BONDS(atom[cur]) == 2 &&
                                             atom[cur].bond_type[0]+atom[cur].bond_type[1] == 2*INCHI_BOND_TYPE_DOUBLE &&
                                             0 == inchi_NUMH2(atom, cur) ) {
                                            next     = atom[cur].neighbor[atom[cur].neighbor[0] == prev];
                                            chain[len++] = next;
                                        } else {
                                            break;
                                        }
                                    }
                                    if ( len > 2 &&
                                         (p2 = IN_NEIGH_LIST( atom[cur].neighbor, (ATOM_NUMBER)prev, AT_NUM_BONDS(atom[cur]))) ) {
                                        sb_ord_from_a2 = p2 - atom[cur].neighbor;
                                        a2 = cur;
                                        /* by design we need to pick up the first non-stereo-bond-neighbor as "sn"-atom */
                                        atom_stereo0D[i].neighbor[0] = atom[a1].neighbor[sb_ord_from_a1 == 0];
                                        atom_stereo0D[i].neighbor[1] = a1;
                                        atom_stereo0D[i].neighbor[2] = a2;
                                        atom_stereo0D[i].neighbor[3] = atom[a2].neighbor[sb_ord_from_a2 == 0];
                                        if ( len % 2 ) {
                                            atom_stereo0D[i].central_atom = chain[len/2];
                                            atom_stereo0D[i].type         = INCHI_StereoType_Allene; 
                                        } else {
                                            atom_stereo0D[i].central_atom = NO_ATOM;
                                        }
                                    } else {
                                        /* error */
                                        atom_stereo0D[i].type = INCHI_StereoType_None;
                                        atom_stereo0D[i].central_atom = NO_ATOM;
                                        atom_stereo0D[i].neighbor[0] =
                                        atom_stereo0D[i].neighbor[3] = -1;
                                        *err |= 64; /* Error in cumulene stereo */
                                        MOLFILE_ERR_SET (*err, 0, "Cumulene stereo not recognized (0D)");

                                    }
#undef MAX_CHAIN_LEN 
                                } else {
                                    /****** a normal possibly stereogenic bond -- not an allene or cumulene *******/
                                    /* by design we need to pick up the first non-stereo-bond-neighbor as "sn"-atom */
                                    sb_ord_from_a1 = p1 - atom[a1].neighbor;
                                    sb_ord_from_a2 = p2 - atom[a2].neighbor;
                                    atom_stereo0D[i].neighbor[0] = atom[a1].neighbor[p1 == atom[a1].neighbor];
                                    atom_stereo0D[i].neighbor[3] = atom[a2].neighbor[p2 == atom[a2].neighbor];
                                    atom_stereo0D[i].central_atom = NO_ATOM;
                                }
                                if ( atom_stereo0D[i].type != INCHI_StereoType_None &&
                                     sb_ord_from_a1 >= 0 && sb_ord_from_a2 >= 0 &&
                                     ATOM_PARITY_WELL_DEF( SB_PARITY_2(atom_stereo0D[i].parity) ) ) {
                                    /* Detected well-defined disconnected stereo
                                     * locate first non-metal neighbors */
                                    int    a, n, j, /* k,*/ sb_ord, cur_neigh, min_neigh;
                                    for ( k = 0; k < 2; k ++ ) {
                                        a      = k? atom_stereo0D[i].neighbor[2] : atom_stereo0D[i].neighbor[1];
                                        sb_ord = k? sb_ord_from_a2 : sb_ord_from_a1;
                                        min_neigh = num_atoms;
                                        for ( n  = j = 0; j < AT_NUM_BONDS(atom[a]); j ++ ) {
                                            cur_neigh = atom[a].neighbor[j];
                                            if ( j != sb_ord && !IS_METAL_ATOM(atom, cur_neigh) ) {
                                                min_neigh = inchi_min( cur_neigh, min_neigh );
                                            }
                                        }
                                        if ( min_neigh < num_atoms ) {
                                            atom_stereo0D[i].neighbor[k?3:0] = min_neigh;
                                        } else {
                                            MOLFILE_ERR_SET (*err, 0, "Cannot find non-metal stereobond neighor (0D)");
                                        }
                                    }
                                }

                                break;
                            }
                        }
                        /* end of 0D parities extraction */
/*exit_cycle:;*/
                    }
#if( defined(INCHI_LIBRARY) || defined(INCHI_MAIN) )
#else
                    /* transfer atom_stereo0D[] to atom[] */
                    if ( len_stereo0D ) {
                        Extract0DParities( atom, num_atoms, atom_stereo0D, len_stereo0D, pStrErr, err );
                    }
#endif
                    if ( pInpAtomFlags ) {
                        /* save chirality flag */
                        *pInpAtomFlags |= InpAtomFlags;
                    }
                } else
                if ( atom ) {
                    inchi_free( atom );
                    atom = NULL;
                }
#if( defined(INCHI_LIBRARY) || defined(INCHI_MAIN) )
#else
                if ( szCoord ) {
                    *szCoord = pszCoord;
                    pszCoord = NULL;
                } else
#endif
                if ( pszCoord ) {
                    inchi_free( pszCoord );
                }
                goto bypass_end_of_INChI;
                /*return num_atoms;*/
            }
        }
        if ( atom_stereo0D ) {
            FreeInchi_Stereo0D( &atom_stereo0D );
        }
        /* end of struct. reading cycle, code never used? */
        if ( res <= 0 ) {
            if ( *err == INCHI_INP_ERROR_ERR ) {
                return num_atoms;
            }
            *err = INCHI_INP_EOF_ERR;
            return INCHI_INP_EOF_RET; /* no more data */
        }
bypass_end_of_INChI:
        /* cleanup */
        if ( num_atoms == INCHI_INP_ERROR_RET && atom_stereo0D ) {
            if ( stereo0D && *stereo0D == atom_stereo0D ) {
                *stereo0D     = NULL;
                *num_stereo0D = 0;
            }
            FreeInchi_Stereo0D( &atom_stereo0D );
        }
        if ( !memcmp(szLine, sStructHdrXmlEnd, sizeof(sStructHdrXmlEnd)-1) )
            num_struct --;
        if ( !memcmp(szLine, sStructHdrXml, sizeof(sStructHdrXml)-1) )
            num_struct ++;

        while ( num_struct > 0 && 0 < my_fgets( szLine, sizeof(szLine)-1, inp_molfile, &bTooLongLine ) ) {
            if ( !memcmp(szLine, sStructHdrXmlEnd, sizeof(sStructHdrXmlEnd)-1) )
                num_struct --;
            else
            if ( !memcmp(szLine, sStructHdrXml, sizeof(sStructHdrXml)-1) )
                num_struct ++;
        }
        return num_atoms;

    }
    
    return num_atoms;

#undef AT_NUM_BONDS
#undef ATOM_NUMBER
#undef IN_NEIGH_LIST
#undef inchi_NUMH2

#if( defined(INCHI_LIBRARY) || defined(INCHI_MAIN) )
#else 
#undef inchi_Atom
#endif

#undef AT_NUM_BONDS 
#undef ATOM_NUMBER        
#undef IN_NEIGH_LIST      
#undef inchi_NUMH2 
#undef INChITo_Atom       
#undef MoreParms          
#undef INPUT_FILE         
#undef Create_Atom        
#undef AT_BONDS_VAL
#undef ISOLATED_ATOM      
#undef NUM_ISO_Hk
#undef IS_METAL_ATOM


}
#ifdef INCHI_MAIN

/**********************************************************************************/
int INChIToInchi_Input( INCHI_FILE *inp_molfile, inchi_Input *orig_at_data, int bMergeAllInputStructures,
                       int bDoNotAddH, INPUT_TYPE nInputType,
                       char *pSdfLabel, char *pSdfValue, long *lSdfId, INCHI_MODE *pInpAtomFlags,
                       int *err, char *pStrErr )
{
    /* inp_ATOM       *at = NULL; */
    int             num_dimensions_new;
    int             num_inp_bonds_new;
    int             num_inp_atoms_new;
    int             num_inp_0D_new;
    inchi_Atom     *at_new     = NULL;
    inchi_Atom     *at_old     = NULL;
    inchi_Stereo0D *stereo0D_new = NULL;
    inchi_Stereo0D *stereo0D_old = NULL;
    int             nNumAtoms  = 0, nNumStereo0D = 0;
    MOL_COORD      *szCoordNew = NULL;
    MOL_COORD      *szCoordOld = NULL;
    int            i, j;

    if ( pStrErr ) {
        pStrErr[0] = '\0';
    }

    /*FreeOrigAtData( orig_at_data );*/
    if ( lSdfId )
        *lSdfId = 0;
    do {
        
        at_old       = orig_at_data? orig_at_data->atom      : NULL; /*  save pointer to the previous allocation */
        stereo0D_old = orig_at_data? orig_at_data->stereo0D  : NULL;
        szCoordOld = NULL;
        num_inp_atoms_new =
            INChIToInchi_Atom( inp_molfile, orig_at_data? &stereo0D_new:NULL, &num_inp_0D_new,
                          bDoNotAddH, nInputType, orig_at_data? &at_new:NULL, MAX_ATOMS,
                          &num_dimensions_new, &num_inp_bonds_new,
                          pSdfLabel, pSdfValue, lSdfId, pInpAtomFlags, err, pStrErr );
        if ( num_inp_atoms_new <= 0 && !*err ) {
            MOLFILE_ERR_SET (*err, 0, "Empty structure");
            *err = 98;
        } else
        if ( orig_at_data && !num_inp_atoms_new && 10 < *err && *err < 20 && orig_at_data->num_atoms > 0 && bMergeAllInputStructures ) {
            *err = 0; /* end of file */
            break;
        } else
        if ( num_inp_atoms_new > 0 && orig_at_data ) {
            /*  merge pOrigDataTmp + orig_at_data => pOrigDataTmp; */
            nNumAtoms    = num_inp_atoms_new + orig_at_data->num_atoms;
            nNumStereo0D = num_inp_0D_new    + orig_at_data->num_stereo0D;
            if ( nNumAtoms >= MAX_ATOMS ) {
                MOLFILE_ERR_SET (*err, 0, "Too many atoms");
                *err = 70;
                orig_at_data->num_atoms = -1;
            } else
            if ( !at_old ) {
                /* the first structure */
                orig_at_data->atom         = at_new;            at_new            = NULL;
                orig_at_data->num_atoms    = num_inp_atoms_new; num_inp_atoms_new = 0;
                orig_at_data->stereo0D     = stereo0D_new;      stereo0D_new      = NULL;
                orig_at_data->num_stereo0D = num_inp_0D_new;    num_inp_0D_new    = 0;
            } else
            if ( orig_at_data->atom = CreateInchi_Atom( nNumAtoms ) ) {
                /*  switch at_new <--> orig_at_data->at; */
                if ( orig_at_data->num_atoms ) {
                    memcpy( orig_at_data->atom, at_old, orig_at_data->num_atoms * sizeof(orig_at_data->atom[0]) );
                    /*  adjust numbering in the newly read structure */
                    for ( i = 0; i < num_inp_atoms_new; i ++ ) {
                        for ( j = 0; j < at_new[i].num_bonds; j ++ ) {
                            at_new[i].neighbor[j] += orig_at_data->num_atoms;
                        }
                    }
                }
                FreeInchi_Atom( &at_old );
                /*  copy newly read structure */
                memcpy( orig_at_data->atom + orig_at_data->num_atoms,
                        at_new,
                        num_inp_atoms_new * sizeof(orig_at_data->atom[0]) );
                /*  cpy newly read 0D stereo */
                if ( num_inp_0D_new > 0 && stereo0D_new ) {
                    if ( orig_at_data->stereo0D = CreateInchi_Stereo0D( nNumStereo0D ) ) {
                        memcpy( orig_at_data->stereo0D, stereo0D_old, orig_at_data->num_stereo0D * sizeof(orig_at_data->stereo0D[0]) );
                        /*  adjust numbering in the newly read structure */
                        for ( i = 0; i < num_inp_0D_new; i ++ ) {
                            if ( stereo0D_new[i].central_atom >= 0 ) {
                                stereo0D_new[i].central_atom += orig_at_data->num_atoms;
                            }
                            for ( j = 0; j < 4; j ++ ) {
                                stereo0D_new[i].neighbor[j] += orig_at_data->num_atoms;
                            }
                        }
                        FreeInchi_Stereo0D( &stereo0D_old );
                        memcpy( orig_at_data->stereo0D+orig_at_data->num_stereo0D,
                                stereo0D_new,
                                num_inp_0D_new * sizeof(orig_at_data->stereo0D[0]) );
                    } else {
                        num_inp_0D_new = 0;
                        MOLFILE_ERR_SET (*err, 0, "Out of RAM");
                        *err = -1;
                    }
                } else {
                    num_inp_0D_new = 0;
                }
                /* update lengths */
                orig_at_data->num_atoms    += num_inp_atoms_new;
                orig_at_data->num_stereo0D += num_inp_0D_new;
            } else {
                MOLFILE_ERR_SET (*err, 0, "Out of RAM");
                *err = -1;
            }
        } else
        if ( num_inp_atoms_new > 0 ) {
            nNumAtoms += num_inp_atoms_new;
        }
        FreeInchi_Atom( &at_new );
        num_inp_atoms_new = 0;
        FreeInchi_Stereo0D( &stereo0D_new );
        num_inp_0D_new = 0;

    } while ( !*err && bMergeAllInputStructures );
    /*
    if ( !*err ) {
        orig_at_data->num_components =
            MarkDisconnectedComponents( orig_at_data );
        if ( orig_at_data->num_components == 0 ) {
            MOLFILE_ERR_SET (*err, 0, "No components found");
            *err = 99;
        }
        if ( orig_at_data->num_components < 0 ) {
            MOLFILE_ERR_SET (*err, 0, "Too many components");
            *err = 99;
        }
    }
    */
    if ( szCoordNew ) {
        inchi_free( szCoordNew );
    }
    if ( at_new ) {
        inchi_free( at_new );
    }
    /*
    if ( !*err ) {
        if ( ReconcileAllCmlBondParities( orig_at_data->atom, orig_at_data->num_atoms ) ) {
            MOLFILE_ERR_SET (*err, 0, "Cannot reconcile stereobond parities");
            if (!orig_at_data->num_atoms) {
                *err = 1;
            }
        }
    }
    */
    if ( *err ) {
        FreeInchi_Input( orig_at_data );
    }
    if ( *err && !(10 < *err && *err < 20) && pStrErr && !pStrErr[0] ) {
        MOLFILE_ERR_SET (*err, 0, "Unknown error");  /*   <BRKPT> */
    }
    return orig_at_data? orig_at_data->num_atoms : nNumAtoms;
}

#endif

#ifndef INCHI_MAIN
#undef AB_MAX_WELL_DEFINED_PARITY
#undef AB_MIN_WELL_DEFINED_PARITY
#include "extr_ct.h"
/****************************************************************************************/
int Extract0DParities( inp_ATOM *at, int nNumAtoms, inchi_Stereo0D *stereo0D,
                       int num_stereo0D, char *pStrErr, int *err )
{
    if ( stereo0D && num_stereo0D > 0 ) {
        int i0D, a2, k, k_prev, type, j, j1, j2, len, parity, parityNM;
        int sb_ord_from_i1, sb_ord_from_i2, sn_ord_from_i1, sn_ord_from_i2;
        AT_NUMB i1n, i2n, i1, i2;
        for ( i0D = 0; i0D < num_stereo0D; i0D ++ ) {
            parity   = (stereo0D[i0D].parity & SB_PARITY_MASK);
            parityNM = (stereo0D[i0D].parity & SB_PARITY_FLAG) >> SB_PARITY_SHFT;
            if ( parity == INCHI_PARITY_NONE ||
                 parity != INCHI_PARITY_ODD && parity != INCHI_PARITY_EVEN &&
                 parity != INCHI_PARITY_UNKNOWN && parity != INCHI_PARITY_UNDEFINED ) {
                char szTemp[16];
                sprintf( szTemp, "#%d", i0D+1 );
                MOLFILE_ERR_SET (*err, 0, "Wrong 0D stereo descriptor(s):");
                MOLFILE_ERR_SET (*err, 0, szTemp);
                continue; /* warning */
            }
            type   = stereo0D[i0D].type;
            a2     = stereo0D[i0D].central_atom; /* central atom or -1 */
            j   = -1;
            len =  0;
            sb_ord_from_i1 = sb_ord_from_i2 = sn_ord_from_i1 = sn_ord_from_i2 = -1;
            i1n = i2n = i1 = i2 = MAX_ATOMS+1;

            if ( (type == INCHI_StereoType_Tetrahedral ||
                  type == INCHI_StereoType_Allene ) &&
                  0 <= a2 && a2 < nNumAtoms ||
                  type == INCHI_StereoType_DoubleBond &&
                  a2 == NO_ATOM) {
                /* test the quadruplet */
                for ( j = 0, k_prev = -1; j < 4; j ++, k_prev = k ) {
                    k = stereo0D[i0D].neighbor[j];
                    if ( k < 0 || k >= nNumAtoms || k_prev == k )
                        break;
                    /* tetrahedral atom connectivity test */
                    if ( type == INCHI_StereoType_Tetrahedral &&
                         k != a2 &&
                         !is_in_the_list( at[a2].neighbor, (AT_NUMB)k, at[a2].valence) ) {
                        break;
                    }
                    /* Double bond, Cumulene and allene are tested in the next if() */
                }
            }
            /* find in the adjacency lists the double bond neighbor that leads to the opposite atom */
            if ( j == 4 && (type == INCHI_StereoType_Allene ||
                            type == INCHI_StereoType_DoubleBond) ) {
                AT_NUMB *p1 = NULL, *p2 = NULL, *q1 = NULL, *q2 = NULL;
                i1n = (AT_NUMB)stereo0D[i0D].neighbor[0];
                i1  = (AT_NUMB)stereo0D[i0D].neighbor[1];
                i2  = (AT_NUMB)stereo0D[i0D].neighbor[2];
                i2n = (AT_NUMB)stereo0D[i0D].neighbor[3];
                /* find q1 and q2 */
                if ( !(q1 = is_in_the_list( at[i1].neighbor, i1n, at[i1].valence)) ||
                     !(q2 = is_in_the_list( at[i2].neighbor, i2n, at[i2].valence)) ) {
                    j = -2; /* error flag */
                } else
                /* allene or cumulene; follow double bonds from i1 to i2 */
                if ( !(p1 = is_in_the_list( at[i1].neighbor, i2, at[i1].valence)) ) {
                    /* at[i1] and at[i2] are not connected: can be only allene or cumulene */
                    AT_NUMB prev, cur, next;
                    int     num_dbond, i, next_ord, half_len;

                    cur = next = i1;
                    len = half_len = 0;
                    while ( len < 20 ) { /* arbitrary very high upper limit to prevent infinite loop */
                        prev = cur;
                        cur  = next;
                        for ( i = 0, num_dbond = 0; i < at[cur].valence; i ++ ) {
                            /* follow double bond path && avoid going back */
                            if ( at[cur].bond_type[i] == BOND_TYPE_DOUBLE &&
                                 prev != at[cur].neighbor[i] ) {
                                next = at[cur].neighbor[i];
                                next_ord = i;
                                num_dbond ++;
                            }
                        }
                        if ( num_dbond == 1 && next != i1 ) {
                            len ++;
                            if ( len == 1 ) {
                                sb_ord_from_i1 = next_ord;
                            }
                            if ( type == INCHI_StereoType_Allene && next == (AT_NUMB)a2 ) {
                                half_len = len;
                            }
                        } else {
                            break;
                        }
                    }
                    if ( cur == i2 && prev != cur && 0 == num_dbond && len > 1 &&
                         (p2 = is_in_the_list( at[i2].neighbor, prev, at[i2].valence)) &&
                         (type != INCHI_StereoType_Allene || len == 2*half_len )) {
                        sb_ord_from_i2 = p2 - at[i2].neighbor;
                        sn_ord_from_i1 = q1 - at[i1].neighbor;
                        sn_ord_from_i2 = q2 - at[i2].neighbor;
                    } else {
                        j = -5; /* error flag */
                    }
                } else
                /* allene must have been already processed, otherwise error */
                if ( type == INCHI_StereoType_Allene ) {
                    /* error: atoms #1 and #2 of allene are connected */
                    j = -3; /* error flag */
                } else
                /* double bond only; the bond type is not checked because at the end
                   of the normalization it may happen to be alternating */
                if ( type == INCHI_StereoType_DoubleBond &&
                     (p2 = is_in_the_list( at[i2].neighbor, i1, at[i2].valence) ) ) {
                    sb_ord_from_i1 = p1 - at[i1].neighbor;
                    sb_ord_from_i2 = p2 - at[i2].neighbor;
                    sn_ord_from_i1 = q1 - at[i1].neighbor;
                    sn_ord_from_i2 = q2 - at[i2].neighbor;
                } else {
                    j = -4; /* error flag */
                }
            }
            if ( j != 4 ) {
                char szTemp[16];
                sprintf( szTemp, "#%d", i0D+1 );
                MOLFILE_ERR_SET (*err, 0, "Wrong 0D stereo descriptor(s):");
                MOLFILE_ERR_SET (*err, 0, szTemp);
                continue; /* error */
            }

            switch ( type ) {
            case INCHI_StereoType_None:
                continue;
            case INCHI_StereoType_DoubleBond:
            case INCHI_StereoType_Allene:
                for ( j1 = 0; j1 < MAX_NUM_STEREO_BONDS && at[i1].sb_parity[j1]; j1 ++ )
                    ;
                for ( j2 = 0; j2 < MAX_NUM_STEREO_BONDS && at[i2].sb_parity[j2]; j2 ++ )
                    ;
                if ( j1 < MAX_NUM_STEREO_BONDS && j2 < MAX_NUM_STEREO_BONDS &&
                     sb_ord_from_i1 >= 0 && sb_ord_from_i2 >= 0 &&
                     sn_ord_from_i1 >= 0 && sn_ord_from_i2 >= 0) {

                    switch( parity ) {
                    case INCHI_PARITY_ODD:
                        at[i1].sb_parity[j1] = AB_PARITY_ODD;
                        at[i2].sb_parity[j2] = AB_PARITY_EVEN;
                        break;
                    case INCHI_PARITY_EVEN:
                        at[i1].sb_parity[j1] = AB_PARITY_ODD;
                        at[i2].sb_parity[j2] = AB_PARITY_ODD;
                        break;
                    case INCHI_PARITY_UNKNOWN:
                        at[i1].sb_parity[j1] = AB_PARITY_UNKN;
                        at[i2].sb_parity[j2] = AB_PARITY_UNKN;
                        break;
                    case INCHI_PARITY_UNDEFINED:
                        at[i1].sb_parity[j1] = AB_PARITY_UNDF;
                        at[i2].sb_parity[j2] = AB_PARITY_UNDF;
                        break;
                    default:
                        at[i1].sb_parity[j1] = AB_PARITY_NONE;
                        at[i2].sb_parity[j2] = AB_PARITY_NONE;
                    }

                    switch( parityNM ) {
                    case INCHI_PARITY_ODD:
                        at[i1].sb_parity[j1] |= AB_PARITY_ODD << SB_PARITY_SHFT;
                        at[i2].sb_parity[j2] |= AB_PARITY_EVEN << SB_PARITY_SHFT;
                        break;
                    case INCHI_PARITY_EVEN:
                        at[i1].sb_parity[j1] |= AB_PARITY_ODD << SB_PARITY_SHFT;
                        at[i2].sb_parity[j2] |= AB_PARITY_ODD << SB_PARITY_SHFT;
                        break;
                    case INCHI_PARITY_UNKNOWN:
                        at[i1].sb_parity[j1] |= AB_PARITY_UNKN << SB_PARITY_SHFT;
                        at[i2].sb_parity[j2] |= AB_PARITY_UNKN << SB_PARITY_SHFT;
                        break;
                    case INCHI_PARITY_UNDEFINED:
                        at[i1].sb_parity[j1] |= AB_PARITY_UNDF << SB_PARITY_SHFT;
                        at[i2].sb_parity[j2] |= AB_PARITY_UNDF << SB_PARITY_SHFT;
                        break;
                    default:
                        break;
                    }

                    at[i1].sb_ord[j1]         = sb_ord_from_i1;
                    at[i1].sn_ord[j1]         = sn_ord_from_i1;
                    at[i1].sn_orig_at_num[j1] = at[i1n].orig_at_number;

                    at[i2].sb_ord[j2]         = sb_ord_from_i2;
                    at[i2].sn_ord[j2]         = sn_ord_from_i2;
                    at[i2].sn_orig_at_num[j2] = at[i2n].orig_at_number;
                }
                break;
            case INCHI_StereoType_Tetrahedral:
                    switch( parity ) {
                    case INCHI_PARITY_ODD:
                        at[a2].p_parity = AB_PARITY_ODD;
                        break;
                    case INCHI_PARITY_EVEN:
                        at[a2].p_parity = AB_PARITY_EVEN;
                        break;
                    case INCHI_PARITY_UNKNOWN:
                        at[a2].p_parity = AB_PARITY_UNKN;
                        break;
                    case INCHI_PARITY_UNDEFINED:
                        at[a2].p_parity = AB_PARITY_UNDF;
                        break;
                    default:
                        continue;
                    }
                for ( j = 0; j < 4; j ++ ) {
                    k = stereo0D[i0D].neighbor[j];
                    at[a2].p_orig_at_num[j] = at[k].orig_at_number;
                }
                break;
            default:
                break;
            }
        }
        /* take care of Unknown stereobonds:                                     */
        /* copy their Unknown stereo descriptors to at->bond_stereo (2005-03-01) */
        FixUnkn0DStereoBonds(at, nNumAtoms);

#ifdef INCHI_LIBRARY

        if ( k = ReconcileAllCmlBondParities( at, nNumAtoms, 0 ) ) {
            char szErrCode[16];
            sprintf( szErrCode, "%d", k);
            AddMOLfileError( pStrErr, "0D Parities Reconciliation failed:" );
            AddMOLfileError( pStrErr, szErrCode );
        }

#endif

    }
    return 0;
}

#endif


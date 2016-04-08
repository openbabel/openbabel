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
 * IUPAC/InChI-Trust Licence for the International Chemical Identifier (InChI)
 * Software version 1.0.
 * Copyright (C) IUPAC and InChI Trust Limited
 *
 * This library is free software; you can redistribute it and/or modify it under the
 * terms of the IUPAC/InChI Trust Licence for the International Chemical Identifier
 * (InChI) Software version 1.0; either version 1.0 of the License, or
 * (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the IUPAC/InChI Trust Licence for the International Chemical Identifier (InChI)
 * Software version 1.0 for more details.
 *
 * You should have received a copy of the IUPAC/InChI Trust Licence for the
 * International Chemical Identifier (InChI) Software version 1.0 along with
 * this library; if not, write to:
 *
 * The InChI Trust
 * c/o FIZ CHEMIE Berlin
 * Franklinstrasse 11
 * 10587 Berlin
 * GERMANY
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "mode.h"

#include "inpdef.h"
#include "ichi.h"
#include "strutil.h"
#include "util.h"
#include "extr_ct.h"
#include "ichitaut.h"
#include "ichinorm.h"
#include "ichicant.h"
#include "ichicano.h"
#include "ichicomn.h"
#include "ichister.h"

#include "ichicomp.h"
#include "ichimain.h"
#include "ichimake.h"

#include "ichi_io.h"

int PrintXmlStartTag(char *pStr,
                     int indent, int bEnd, const char *tag,
                     const char *l1, int v1, const char *l2, int v2,
                     const char *l3, int v3, const char *l4, int v4,
                     const char *l5, int v5, const char *l6, int v6);
int Needs2addXmlEntityRefs(const char *s );
int AddXmlEntityRefs(const char *p, char *d );
#if ( TEST_RENUMB_ATOMS == 1 ) /*  { */
int CompareStereoINChI( INChI_Stereo *s1, INChI_Stereo *s2 );
#endif

int str_LineStart(const char *tag, char *tag2, int val2, char *pStr, int ind );
int str_LineEnd(const char *tag, int tot_len, int nStrLen,
                int *bOverflow, char *pStr, int ind, int bPlainTextTags );
int CleanOrigCoord(MOL_COORD szCoord, int delim );
int WriteOrigCoord(int num_inp_atoms, MOL_COORD *szMolCoord, int *i,
                   char *szBuf, int buf_len);
int WriteOrigAtoms(int num_inp_atoms, inp_ATOM *at, int *i,
                   char *szBuf, int buf_len,
                   STRUCT_DATA *sd);
int WriteOrigBonds(int num_inp_atoms, inp_ATOM *at, int *i,
                   char *szBuf, int buf_len,
                   STRUCT_DATA *sd);

void GetSaveOptLetters(unsigned char save_opt_bits, char* let1, char* let2);


char VER_STRING[64];

const char sCompDelim[]       = ";"; /* component delimiter */
const char sIdenticalValues[] = "*"; /* identical component */
const char x_space[]          = "                  ";

/* xml output: words & additional tags */
const char x_inchi[]          = INCHI_NAME;
const char x_inchi_ver[]      = "version"; /* "InChI.version"; */
const char x_curr_ver[]       = INCHI_VERSION;

const char x_structure[]      = "structure";
const char x_number[]         = "number";
const char x_header[]         = "id.name";
const char x_value[]          = "id.value";

const char x_empty[]          = "";

const char x_type[]           = "type";

const char x_message[]        = "message";
const char x_text[]           = "value";

const char x_ferr[]           = "fatal (aborted)";
const char x_err[]            = "error (no InChI)";
const char x_warn[]           = "warning";

const char x_basic[]          = "identifier";
const char x_tautomeric[]     = "mobile-H";
const char x_reconnected[]    = "reconnected";

const char x_ver[]            = "version";

const char x_type_alpha[]     = "alpha";
const char x_type_numer[]     = "numeric";
const char x_type_predec[]    = "sct";
const char x_type_normal[]    = "normal";
const char x_type_short[]     = "compressed";
const char x_basic_layer[]    = "basic";

const char x_aux_basic[]      = "identifier.auxiliary-info";
const char x_aux_comm[]       = "!-- This section is NOT a part of the identifier, it is not unique --";

const char x_ign_uu_sp2[]     = "omit_undef_dbond";
const char x_ign_uu_sp3[]     = "omit_undef_sp3";

const char x_line_opening[]   = "<";
const char x_line_closing[]   = "</";
const char x_close_line[]     = ">";

const char x_abs[]            = "1";
const char x_rel[]            = "2";
const char x_rac[]            = "3";

#define MAX_TAG_LEN 64

typedef struct tagInchiTag
{
    const char *szPlainLabel;
    const char *szPlainComment;
    const char *szXmlLabel;
    int  bAlwaysOutput;
} INCHI_TAG;

/* identifier */
const INCHI_TAG IdentLbl[] =
{
                                                                  /* prefixes: may be combined in this order */
/* IL_FIXH_ORD, */    { "/",   "fixed_H",        "fixed-H",        0 }, /* fixed H */
/* IL_ISOT_ORD, */    { "/",   "isotopic",       "isotopic",       0 }, /* isotopic */
/* IL_STER_ORD, */    { "/",   "stereo",         "stereo",         0 }, /* stereo */
                                                                   /* items */
/* IL_VERS_ORD, */    { "" ,   "version",        "version",        1 },
/* IL_FML__ORD, */    { "/",   "formula",        "formula",        1 }, /* basic part formula */
/* IL_CONN_ORD, */    { "/c",  "connections",    "connections",    1 },
/* IL_ALLH_ORD, */    { "/h",  "H_atoms",        "H",              1 },
/* IL_CHRG_ORD, */    { "/q",  "charge",         "charge",         1 },
/* IL_PROT_ORD, */    { "/p",  "protons",        "protons",        0 },
                                                                    /* stereo */
/* IL_DBND_ORD, */    { "/b",  "dbond",          "dbond",          0 },
/* IL_SP3S_ORD, */    { "/t",  "sp3",            "sp3",            0 },
/* IL_INVS_ORD, */    { "/m",  "sp3:inverted",   "abs.inverted",   0 }, /* mirrored */
/* IL_TYPS_ORD, */    { "/s",  "type (1=abs, 2=rel, 3=rac)", "type",           0 }, /* stereo type */
                                                                    /* isotopic */
/* IL_ATMS_ORD, */    { "/i",  "atoms",          "atoms",          1 },
                                                                    /* isotopic mobile H only */
/* IL_XCGA_ORD, */    { "/h",  "exchangeable_H", "H-isotopic",     1 },
                                                                    /* fixed H only */
/* IL_FMLF_ORD, */    { "/f",  "formula",        "formula",        1 }, /* fixed H formula */
/* IL_HFIX_ORD, */    { "/h",  "H_fixed" ,       "H-fixed" ,       1 }, /* fixed-H */
/* IL_TRNS_ORD, */    { "/o",  "transposition",  "transposition",  0 }, /* order */
/* IL_REC__ORD, */    { "/r",  "reconnected bond(s) to metal(s) formula",  "formula",  0 }
};

/*

  Parsing plain text InChI (FML is a chemical formula)
  ========================

  1.12Beta/FML       /i      /f[FML]  /i   [/o] /rFML      /i      /f[FML]  /i   [/o] end
          |          |       |        |         |          |       |        |         |
Labels    | chqpbtms | hbtms | hqbtms | btms    | chqpbtms | hbtms | hqbtms | btms    |
inside:   |          |       |        |         |          |       |        |         |
          | non-iso- | iso-  | fix-   | iso-    | non-iso- | iso-  | fix-   | iso-    |
meaning:  | topic    | topic | ed H   | topic   | topic    | topic | ed H   | topic   |
          |----------+-------+--------+---------|----------+-------+--------+---------|
          |        mobile-H  |   fixed-H        |        mobile-H  |   fixed-H        |
          |----------+-------+--------+---------|----------+-------+--------+---------|
          |                                     |                                     |
          |     normal  or disconected metal    |      reconnected bonds to metal     |
          |_____________________________________|_____________________________________|

  meanings of h:

       /h   -  immobile H & mobile H group(s)
     /i/h   -  exchangeable isotopic H (common)
     /f/h   -  fixed-H
     /f/i/h -  never happens

*/

typedef enum tagIdentLblOrd
{
    IL_FIXH_ORD,
    IL_ISOT_ORD,
    IL_STER_ORD,

    IL_VERS_ORD,
    IL_FML__ORD,
    IL_CONN_ORD,
    IL_ALLH_ORD,
    IL_CHRG_ORD,
    IL_PROT_ORD,

    IL_DBND_ORD,
    IL_SP3S_ORD,
    IL_INVS_ORD,
    IL_TYPS_ORD,

    IL_ATMS_ORD,

    IL_XCGA_ORD,

    IL_FMLF_ORD,
    IL_HFIX_ORD,
    IL_TRNS_ORD,
    IL_REC__ORD,

    IL_MAX_ORD /* max number of tags */
} IDENT_LBL_ORD;

typedef enum tagIdentLblBit
{
    IL_FIXH = 1 << IL_FIXH_ORD,
    IL_ISOT = 1 << IL_ISOT_ORD,
    IL_STER = 1 << IL_STER_ORD,

    IL_VERS = 1 << IL_VERS_ORD,
    IL_FML_ = 1 << IL_FML__ORD,
    IL_CONN = 1 << IL_CONN_ORD,
    IL_ALLH = 1 << IL_ALLH_ORD,
    IL_CHRG = 1 << IL_CHRG_ORD,
    IL_PROT = 1 << IL_PROT_ORD,

    IL_DBND = 1 << IL_DBND_ORD,
    IL_SP3S = 1 << IL_SP3S_ORD,
    IL_INVS = 1 << IL_INVS_ORD,
    IL_TYPS = 1 << IL_TYPS_ORD,

    IL_ATMS = 1 << IL_ATMS_ORD,

    IL_XCGA = 1 << IL_XCGA_ORD,

    IL_FMLF = 1 << IL_FMLF_ORD,
    IL_HFIX = 1 << IL_HFIX_ORD,
    IL_TRNS = 1 << IL_TRNS_ORD,
    IL_REC_ = 1 << IL_REC__ORD

} IDENT_LBL_BIT;

/* aux info */
const INCHI_TAG AuxLbl[] =
{
                                                                    /* prefixes may be combined in this order */
/* AL_FIXH_ORD, */    { "/",     "fixed_H",                "fixed-H",             0 }, /* fixed-H */
/* AL_ISOT_ORD, */    { "/",     "isotopic",               "isotopic",            0 }, /* isotopic */
/* AL_STER_ORD, */    { "/",     "abs_stereo_inverted",    "stereo.abs.inverted", 0 }, /* inv abs sp3 stereo */
/* AL_REVR_ORD, */    { "/",     "reversibility",          "reversibility",       0 }, /* reversibility */
                                                                       /* items */
/* AL_VERS_ORD, */    { "",      "version",                "version",             1 },
/* AL_NORM_ORD, */    { "/",     "normalization_type",     "norm-type",           1 },
/* AL_ANBR_ORD, */    { "/N:",   "original_atom_numbers",  "atom.orig-nbr",       1 },
/* AL_AEQU_ORD, */    { "/E:",   "atom_equivalence",       "atom.equivalence",    0 },
/* AL_GEQU_ORD, */    { "/gE:",  "group_equivalence",      "group.equivalence",   0 },
                                                                       /* inv abs sp3 stereo */
/* AL_SP3I_ORD, */    { "/it:",  "sp3",                    "sp3",                 0 },
/* AL_SP3N_ORD, */    { "/iN:",  "original_atom_numbers",  "atom.orig-nbr",       0 },

/* AL_CRV__ORD, */    { "/CRV:", "charge_radical_valence", "charges-rad-val",     0 },
                                                                       /* reversibility */
/* AL_ATMR_ORD, */    { "/rA:",  "atoms",                  "atoms",               0 },
/* AL_BNDR_ORD, */    { "/rB:",  "bonds",                  "bonds",               0 },
/* AL_XYZR_ORD, */    { "/rC:",  "xyz",                    "xyz",                 0 },
                                                                       /* fixed-H only */
/* AL_FIXN_ORD, */    { "/F:",   "original_atom_numbers",  "atom.orig-nbr",       1 },
                                                                       /* isotopic only */
/* AL_ISON_ORD, */    { "/I:",   "original_atom_numbers",  "atom.orig-nbr",       1 },

/* AL_REC__ORD, */    { "/R:",  "reconnected bond(s) to metal(s) part",  "",      1 }

};

typedef enum tagAuxLblOrd
{
    AL_FIXH_ORD,
    AL_ISOT_ORD,
    AL_STER_ORD,
    AL_REVR_ORD,

    AL_VERS_ORD,
    AL_NORM_ORD,
    AL_ANBR_ORD,
    AL_AEQU_ORD,
    AL_GEQU_ORD,

    AL_SP3I_ORD,
    AL_SP3N_ORD,

    AL_CRV__ORD,

    AL_ATMR_ORD,
    AL_BNDR_ORD,
    AL_XYZR_ORD,

    AL_FIXN_ORD,

    AL_ISON_ORD,

    AL_REC__ORD,

    AL_MAX_ORD   /* max number of tags */

} AUX_LBL_ORD;


typedef enum tagAuxLblBit
{
    AL_FIXH = 1 << AL_FIXH_ORD,
    AL_ISOT = 1 << AL_ISOT_ORD,
    AL_STER = 1 << AL_STER_ORD,
    AL_REVR = 1 << AL_REVR_ORD,

    AL_VERS = 1 << AL_VERS_ORD,
    AL_NORM = 1 << AL_NORM_ORD,
    AL_ANBR = 1 << AL_ANBR_ORD,
    AL_AEQU = 1 << AL_AEQU_ORD,
    AL_GEQU = 1 << AL_GEQU_ORD,

    AL_SP3I = 1 << AL_SP3I_ORD,
    AL_SP3N = 1 << AL_SP3N_ORD,

    AL_CRV_ = 1 << AL_CRV__ORD,

    AL_ATMR = 1 << AL_ATMR_ORD,
    AL_BNDR = 1 << AL_BNDR_ORD,
    AL_XYZR = 1 << AL_XYZR_ORD,

    AL_FIXN = 1 << AL_FIXN_ORD,

    AL_ISON = 1 << AL_ISON_ORD,

    AL_REC_ = 1 << AL_REC__ORD

} AUX_LBL_BIT;

const int MAX_TAG_NUM = inchi_max((int)IL_MAX_ORD, (int)AL_MAX_ORD);

char *szGetTag(const INCHI_TAG *Tag, int nTag, int bTag, char *szTag, int *bAlways);

#define SP(N)        (x_space+sizeof(x_space)-1-(N))
/**********************************************************************************************/
typedef struct tagXmlEntityRef
{
    char nChar;
    const char *pRef;
} X_REF;
const X_REF xmlRef[] = { {'<', "&lt;"}, {'&', "&amp;"}, {'>', "&gt;"}, {'"', "&quot;"}, {'\'', "&apos;"}, {0, NULL}, };
const char szRefChars[sizeof(xmlRef)/sizeof(xmlRef[0])] = {'<', '&', '>', '"', '\'', '\0' };
/**********************************************************************************************/
int PrintXmlStartTag(char *pStr, int indent, int bEnd, const char *tag,
                     const char *l1, int v1, const char *l2, int v2,
                     const char *l3, int v3, const char *l4, int v4,
                     const char *l5, int v5, const char *l6, int v6)
{
    int len=0;
    if ( tag ) {
        len += sprintf( pStr+len, "%s<%s", SP(indent), tag);
    }
    if ( l1 ) {
        len += sprintf( pStr+len, " %s=\"%d\"", l1, v1);
    }
    if ( l2 ) {
        len += sprintf( pStr+len, " %s=\"%d\"", l2, v2);
    }
    if ( l3 ) {
        len += sprintf( pStr+len, " %s=\"%d\"", l3, v3);
    }
    if ( l4 ) {
        len += sprintf( pStr+len, " %s=\"%d\"", l4, v4);
    }
    if ( l5 ) {
        len += sprintf( pStr+len, " %s=\"%d\"", l5, v5);
    }
    if ( l6 ) {
        len += sprintf( pStr+len, " %s=\"%d\"", l6, v6);
    }
    if ( (bEnd & 3) ) {
        len += sprintf( pStr+len, "%s%s", (bEnd & 1)?"/":"", (bEnd & 2)?">":"");
    }
    return len;
}

/**********************************************************************************************/
int Needs2addXmlEntityRefs( const char *s )
{
    int len = 0;
    const X_REF *q = xmlRef, *r;
    const char  *p;
    if ( s && *s ) {
        for ( q = xmlRef, len = 0; q->nChar; q ++ ) {
            for ( p = s; (p = strchr( p, q->nChar )); p ++ ) {
                if ( q->nChar == '&' ) {
                    for ( r = xmlRef; r->nChar; r ++ ) {
                        if ( !memcmp( p, r->pRef, strlen(r->pRef) ) )
                            goto DoNotSubstitute;
                    }
                }
                len += strlen(q->pRef)-1;
DoNotSubstitute:;
            }
        }
        if ( len ) {
            len += strlen( s );
        }
    }
    return len;
}

/**********************************************************************************************/
int AddXmlEntityRefs( const char *p, char *d )
{
    int len_d, n;
    const X_REF *q = xmlRef, *r;

    len_d = 0;
    while ( *p ) {
        n = strcspn( p, szRefChars );
        if ( n > 0 ) {
            /*  first n characters of p do not contain referenced chars; copy them */
            strncpy( d+len_d, p, n ); /*  does not have zero termination */
            len_d += n;  /*  new destination length */
            p += n;      /*  position of the referenced char in the source */
        }
        if ( *p ) {
            if ( *p == '&' ) {
                for ( r = xmlRef; r->nChar; r ++ ) {
                    if ( !memcmp( p, r->pRef, strlen(r->pRef) ) ) {
                        d[len_d++] = *p;
                        goto DoNotSubstitute;
                    }
                }
            }
            q = xmlRef + (strchr( szRefChars,  UCINT *p) - szRefChars);
            strcpy( d+len_d, q->pRef );   /*  add entity reference and zero termination */
            len_d += strlen( d + len_d ); /*  new destination length */
DoNotSubstitute:
            p ++;
        } else {
            d[len_d] = '\0'; /*  add zero termination */
        }

    }
    return len_d;
}

/**********************************************************************************************/
int OutputINChIXmlRootStartTag( INCHI_IOSTREAM *output_file )
{
    char pStr[128];
    sprintf( pStr, "<%s %s=\"%s\">", x_inchi, x_inchi_ver, x_curr_ver );
    inchi_ios_print_nodisplay( output_file, "%s\n", pStr );
    return 0;
}

/**********************************************************************************************/
int OutputINChIXmlRootEndTag( INCHI_IOSTREAM *output_file )
{
    char pStr[128];
    sprintf( pStr, "</%s>", x_inchi );
    inchi_ios_print_nodisplay( output_file, "%s\n", pStr );
    return 0;
}

/**********************************************************************************************/
int OutputINChIXmlStructStartTag( INCHI_IOSTREAM *output_file, char *pStr, int ind /* indent*/,
                                 int nStrLen, int bNoStructLabels,
                                 int num_input_struct, const char *szSdfLabel, const char *szSdfValue )
{
    char szBuf[64];
    int nEstLen1;
    int nEstLen2;
    int ret = 0;
    int tot_len;
    char *pSdfLabel = NULL, *pSdfValue = NULL, *p;
    /*  substitute special characters (see szRefChars[]) with xml Entity References */
    int   len;
    if ( bNoStructLabels ) {
        /* no labela at all */
        inchi_ios_print( output_file, "%s\n", "" );   /*  empty line */
        tot_len = 0;
        tot_len += sprintf(pStr+tot_len, "%s<%s", SP(ind), x_structure);
        tot_len += sprintf(pStr+tot_len, ">" );
        inchi_ios_print( output_file, "%s\n", pStr );
        ret = 1; /*  success */
    } else
    if ( !(szSdfLabel && szSdfLabel[0]) && !(szSdfValue && szSdfValue[0]) ) {
        /* only structure number if present */
        inchi_ios_print( output_file, "%s\n", "" );   /*  empty line */
        tot_len = 0;
        tot_len += sprintf(pStr+tot_len, "%s<%s", SP(ind), x_structure);
        if ( num_input_struct > 0 ) {
            tot_len += sprintf(pStr+tot_len, " %s=\"%d\"", x_number, num_input_struct);
        }
        tot_len += sprintf(pStr+tot_len, ">" );
        inchi_ios_print( output_file, "%s\n", pStr );
        ret = 1; /*  success */
    } else {
        if ( (len = Needs2addXmlEntityRefs( szSdfLabel )) ) {
            if ( (p = (char*) inchi_malloc( len+1 )) ) {
                AddXmlEntityRefs( szSdfLabel, p );
                szSdfLabel = pSdfLabel = p;
            }
        }
        if ( (len = Needs2addXmlEntityRefs( szSdfValue )) ) {
            if ( (p = (char*) inchi_malloc( len+1 )) ) {
                AddXmlEntityRefs( szSdfValue, p );
                szSdfValue = pSdfValue = p;
            }
        }
        nEstLen1 = ind + 1 + sizeof(x_structure)-1
                       + 1 + sizeof(x_number)-1 + 1 + sprintf(szBuf,"\"%d\"", num_input_struct)  + 2;
        nEstLen2 = 1 + sizeof(x_header)-1 + 1 + 2 + (szSdfLabel? strlen(szSdfLabel):0)
                 + 1 + sizeof(x_value) -1 + 1 + 2 + (szSdfValue? strlen(szSdfValue):0) + 2;
        if ( nEstLen1 <= nStrLen ) {
            inchi_ios_print( output_file, "%s\n", "" );   /*  empty line */
            tot_len = 0;
            tot_len += sprintf(pStr+tot_len, "%s<%s", SP(ind), x_structure);
            tot_len += sprintf(pStr+tot_len, " %s=\"%d\"", x_number, num_input_struct);
            if ( nEstLen1 + nEstLen2 <= nStrLen ) {
                tot_len += sprintf(pStr+tot_len, " %s=\"%s\"", x_header, szSdfLabel? szSdfLabel:x_empty);
                tot_len += sprintf(pStr+tot_len, " %s=\"%s\"", x_value, szSdfValue? szSdfValue:x_empty);
            }
            tot_len += sprintf(pStr+tot_len, ">" );
            inchi_ios_print( output_file, "%s\n", pStr );
            ret = 1; /*  success */
        }
        if ( pSdfValue ) {
            inchi_free ( pSdfValue );
        }
        if ( pSdfLabel ) {
            inchi_free( pSdfLabel );
        }
    }
    return ret;  /*  0 => Buffer overflow */
}

/**********************************************************************************************/
int OutputINChIXmlStructEndTag( INCHI_IOSTREAM *output_file, char *pStr, int nStrLen, int ind )
{
    if ( output_file && pStr )
    {
        int nEstLen1 = ind + 1 + 1 + sizeof(x_structure)-1 + 2;
        if ( nEstLen1 <= nStrLen )
        {
            sprintf(pStr, "%s</%s>", SP(ind), x_structure);
            inchi_ios_print( output_file, "%s\n", pStr );
            return 1;
        }
    }
    return 0;
}

/**********************************************************************************************/
int OutputINChIXmlError( INCHI_IOSTREAM *output_file, char *pStr, int nStrLen, int ind,
                        /*int nErrorNumber,*/ char *pErrorText, int bError )
{
    /* char szBuf[64]; */
    const char *pErr;
    char *pNewErrorText=NULL, *szErrorText = pErrorText;
    int nEstLen, len=0, ret = 0;

    switch( bError ) {
    case _IS_WARNING:
        pErr = x_warn;
        break;
    case _IS_ERROR:
        pErr = x_err;
        break;
    default: /*  _IS_FATAL */
        pErr = x_ferr;
        break;
    }

#if ( ENTITY_REFS_IN_XML_MESSAGES == 1 )
    /*  insert xml entity references if necessary */
    if ( (len = Needs2addXmlEntityRefs( szErrorText )) ) {
        if ( (pNewErrorText = (char*) inchi_malloc( len+1 )) ) {
            AddXmlEntityRefs( szErrorText, pNewErrorText );
            szErrorText = pNewErrorText;
        }
    }
#else
    szErrorText = pErrorText;
#endif


    nEstLen = ind + 1 + sizeof(x_message)-1
                  + 1 + sizeof(x_type)-1 + 1 + 1 + strlen(pErr)-1
                  /* + 1 + sizeof(x_code)-1 + 1 +     sprintf(szBuf, "%d", nErrorNumber) */
                  + 1 + sizeof(x_text)-1 + 1 + 1 + strlen(szErrorText) + 2;
    if ( nEstLen <= nStrLen ) {
        /*
        sprintf( pStr, "%s<%s %s=\"%s\" %s=\"%d\" %s=\"%s\"/>",
                 SP(ind), x_message, x_type, pErr, x_code, nErrorNumber, x_text, szErrorText );
        */
        sprintf( pStr, "%s<%s %s=\"%s\" %s=\"%s\"/>",
                 SP(ind), x_message, x_type, pErr, x_text, szErrorText );
        inchi_ios_print( output_file, "%s\n", pStr );
        /*
        pErrorText[0] = '\0'; // do not repeat same output
        */
        ret = 1;
    }
    if ( pNewErrorText )
        inchi_free( pNewErrorText );
    return ret;

}

/**********************************************************************************************/
int OutputINChIPlainError( INCHI_IOSTREAM *output_file, char *pStr, int nStrLen,
                           char *pErrorText, int bError )
{
    /* char szBuf[64]; */
    const char *pErr;
    char *pNewErrorText=NULL, *szErrorText = pErrorText;
    int nEstLen, ret = 0;

    switch( bError ) {
    case _IS_WARNING:
        pErr = x_warn;
        break;
    case _IS_ERROR:
        pErr = x_err;
        break;
    default: /*  _IS_FATAL */
        pErr = x_ferr;
        break;
    }
                  /* <%s: >, x_message */
    nEstLen =     sizeof(x_message)-1 + 1 + 1
                  /* <%s=\"%s\">, x_type, pErr */
                  + sizeof(x_type)-1 + 1 + 1 + strlen(pErr) + 1
                  /* < %s=\"%s\"\n>, x_text, szErrorText */
                  + 1 + sizeof(x_text)-1 + 1 + 1 + strlen(szErrorText) + 1 + 1;
    if ( nEstLen < nStrLen ) {
        sprintf( pStr, "%s: %s=\"%s\" %s=\"%s\"",
                 x_message, x_type, pErr, x_text, szErrorText );
        inchi_ios_print( output_file, "%s\n", pStr );
        ret = 1;
    }
    if ( pNewErrorText )
        inchi_free( pNewErrorText );
    return ret;

}

/**************************************************************************/

#ifndef OUT_TN    /* defined in mode.h; quoted here for reference purposes only */

#define OUT_N1              0    /* non-tautomeric only */
#define OUT_T1              1    /* tautomeric if present otherwise non-tautomeric */
#define OUT_NT              2    /* only non-taut representations of tautomeric */
#define OUT_TN              3    /* tautomeric if present otherwise non-tautomeric;
                                    sepatately output non-taut representations of tautomeric if present */
/* OUT_TN = OUT_T1 + OUT_NT */
#endif


/******************************************************************/
const char *EquString( int EquVal )
{
    int bFrom = EquVal & (iiSTEREO | iiSTEREO_INV | iiNUMB | iiEQU );
    int bType = EquVal & (iitISO   | iitNONTAUT );
    int bEq2  = EquVal & (iiEq2NONTAUT | iiEq2ISO | iiEq2INV );
    const char *r = "";

#if ( FIX_EMPTY_LAYER_BUG == 1 )
    int bEmpty= EquVal & iiEmpty;
    if ( bEmpty ) {
        r = "e";
        return r;
    }
#endif

    switch ( bFrom ) {

    case iiSTEREO:  /* ------------ Stereo --------------------*/
        switch ( bType ) {
        case iitISO:  /* iso main stereo =... */
            switch( bEq2 ) {
            case 0:
                r = "m";            /* iso main stereo = main stereo */
                break;
            default:
                r = "??";           /* should not happen */
                break;
            }
            break;
        case iitNONTAUT: /* non-taut stereo =... */
            switch( bEq2 ) {
            case 0:
                r = "m";            /* non-taut stereo = main stereo */
                break;
            default:
                r = "??";           /* should not happen */
                break;
            }
            break;
        case (iitNONTAUT | iitISO): /* iso non-taut stereo = ... */
            switch( bEq2 ) {
            case 0:
                r = "m";            /* iso non-taut stereo = main stereo */
                break;
            case iiEq2ISO:
                r = "M";            /* iso non-taut stereo = main iso stereo */
                break;
            case iiEq2NONTAUT:
                r = "n";            /* iso non-taut stereo = non-taut stereo */
                break;
            default:
                r = "??";           /* should not happen */
                break;
            }
            break;
        default:
            r = "??";           /* should not happen */
            break;
        }
        break;

    case iiSTEREO_INV: /*---------- Inverted Aux Stereo ------*/
        if ( bEq2 & iiEq2INV ) { /* stereo = Inverted(another stereo) */
            bEq2 &= ~iiEq2INV;
            switch( bType ) {
            case 0: /* main = ...*/
                switch( bEq2 ) {
                case 0:
                    r = "im";       /* main         = Inv(main) */
                    break;
                case iiEq2ISO:
                    r = "iM";       /* main         = Inv(main iso) */
                    break;
                case iiEq2NONTAUT:
                    r = "in";       /* maim         = Inv(non-taut) */
                    break;
                case (iiEq2NONTAUT | iiEq2ISO):
                    r = "iN";       /* maim         = Inv(non-taut iso ) */
                    break;
                default:
                    r = "??";           /* should not happen */
                    break;
                }
                break;
            case iitISO: /* main iso = ...*/
                switch( bEq2 ) {
                case 0:
                    r = "im";       /* main iso     = Inv(main) */
                    break;
                case iiEq2ISO:
                    r = "iM";       /* main iso     = Inv(main iso) */
                    break;
                case iiEq2NONTAUT:
                    r = "in";       /* maim iso     = Inv(non-taut) */
                    break;
                case (iiEq2NONTAUT | iiEq2ISO):
                    r = "iN";       /* maim         = Inv(non-taut iso ) */
                    break;
                default:
                    r = "??";           /* should not happen */
                    break;
                }
                break;
            case iitNONTAUT: /* non-taut = ... */
                switch( bEq2 ) {
                case 0:
                    r = "im";       /* non-taut     = Inv(main) */
                    break;
                case iiEq2ISO:
                    r = "iM";       /* non-taut     = Inv(main iso) */
                    break;
                case iiEq2NONTAUT:
                    r = "in";       /* non-taut     = Inv(non-taut) */
                    break;
                case (iiEq2NONTAUT | iiEq2ISO):
                    r = "iN";       /* non-taut     = Inv(non-taut iso ) */
                    break;
                default:
                    r = "??";           /* should not happen */
                    break;
                }
                break;
            case (iitNONTAUT | iitISO):
                switch( bEq2 ) {
                case 0:
                    r = "im";       /* non-taut iso = Inv(main) */
                    break;
                case iiEq2ISO:
                    r = "iM";       /* non-taut iso = Inv(main iso) */
                    break;
                case iiEq2NONTAUT:
                    r = "in";       /* non-taut iso = Inv(non-taut) */
                    break;
                case (iiEq2NONTAUT | iiEq2ISO):
                    r = "iN";       /* non-taut iso = Inv(non-taut iso ) */
                    break;
                default:
                    r = "??";           /* should not happen */
                }
                break;
            default:
                r = "??";           /* should not happen */
                break;
            }

        } else {  /* Inv stereo = another (non-inverted) stereo */

            switch( bType ) {
            case iitISO: /* main iso = ...*/
                switch( bEq2 ) {
                case 0:
                    r = "m";       /* main         = (inverted aux) main */
                    break;
                default:
                    r = "??";           /* should not happen */
                    break;
                }
                break;
            case iitNONTAUT: /* non-taut = ... */
                switch( bEq2 ) {
                case 0:
                    r = "m";       /* non-taut     = (inverted aux) main */
                    break;
                default:
                    r = "??";           /* should not happen */
                    break;
                }
                break;
            case (iitNONTAUT | iitISO): /* non-taut iso = ...*/
                switch( bEq2 ) {
                case 0:
                    r = "m";        /* non-taut iso  = (inverted aux) main */
                    break;
                case iiEq2ISO:
                    r = "M";       /* non-taut iso  = (inverted aux) main iso */
                    break;
                case iiEq2NONTAUT:
                    r = "n";       /* non-taut iso  = (inverted aux) non-taut */
                    break;
                default:
                    r = "??";           /* should not happen */
                    break;
                }
                break;
            default:
                r = "??";           /* should not happen */
                break;
            }
        }
        break;

    case ( iiNUMB | iiSTEREO_INV): /*------------- Inv Stereo Numbering ------------*/
        switch( bType ) {
        case 0: /* inv stereo numb main = ...*/
            switch( bEq2 ) {
            case 0:
                r = "m";       /* inv stereo numb main     = main numb */
                break;
            default:
                r = "??";      /* should not happen */
                break;
            }
            break;
        case iitISO: /* inv stereo iso numb main = ...*/
            switch( bEq2 ) {
            case 0:
                r = "m";       /* inv stereo iso numb main = main numb  */
                break;
            case iiEq2INV:
                r = "im";      /* inv stereo iso numb main = InvStereo(main) numb */
                break;
            case iiEq2ISO:
                r = "M";      /* inv stereo iso numb main = isotopic main numb */
                break;
            default:
                r = "??";      /* should not happen */
                break;
            }
            break;
        case iitNONTAUT: /* inv stereo numb non-taut = ... */
            switch( bEq2 ) {
            case 0:
                r = "m";       /* inv stereo numb non-taut = main numb */
                break;
            case iiEq2NONTAUT:
                r = "n";       /* inv stereo numb non-taut = non-taut numb */
                break;
            case iiEq2INV:
                r = "im";      /* inv stereo numb non-taut =  InvStereo(main) numb  */
                break;
            default:
                r = "??";      /* should not happen */
                break;
            }
            break;
        case (iitNONTAUT | iitISO): /* inv stereo numb non-taut iso = ... */
            switch( bEq2 ) {
            case 0:
                r = "m";       /* inv stereo numb non-taut iso = main numb */
                break;
            case iiEq2ISO:
                r = "M";       /* inv stereo numb non-taut iso = main numb iso */
                break;
            case (iiEq2ISO | iiEq2INV):
                r = "iM";       /* inv stereo numb non-taut iso = InvStereo(main iso) numb */
                break;
            case iiEq2NONTAUT:
                r = "n";       /* inv stereo numb non-taut iso = non-taut numb */
                break;
            case (iiEq2NONTAUT | iiEq2ISO):
                r = "N";       /* inv stereo numb non-taut iso = non-taut iso numb */
                break;
            case iiEq2INV:
                r = "im";      /* inv stereo numb non-taut iso = InvStereo(main) numb */
                break;
            case (iiEq2NONTAUT | iiEq2INV):
                r = "in";      /* inv stereo numb non-taut iso = InvStereo(non-taut) numb ) */
                break;
            default:
                r = "??";           /* should not happen  */
                break;
            }
            break;
        default:
            r = "??";           /* should not happen */
            break;
        }
        break;

    case iiNUMB:           /*------------- Canonical Numbering ------------*/
        switch( bType ) {
        case 0:         /* numb main = ...*/
            r = "??";      /* should not happen */
            break;
        case iitISO:     /* iso numb main = ...*/
            switch( bEq2 ) {
            case 0:
                r = "m";       /* iso numb main = main numb  */
                break;
            default:
                r = "??";      /* should not happen */
            }
            break;
        case iitNONTAUT: /* numb non-taut = ... */
            switch( bEq2 ) {
            case 0:
                r = "m";       /* numb non-taut = main numb */
                break;
            default:
                r = "??";      /* should not happen */
            }
            break;
        case (iitNONTAUT | iitISO): /* numb non-taut iso = ... */
            switch( bEq2 ) {
            case 0:
                r = "m";       /* numb non-taut iso = main numb */
                break;
            case iiEq2ISO:
                r = "M";       /* numb non-taut iso = main numb iso */
                break;
            case iiEq2NONTAUT:
                r = "n";       /* numb non-taut iso = non-taut numb */
                break;
            default:
                r = "??";           /* should not happen */
                break;
            }
            break;
        default:
            r = "??";           /* should not happen */
            break;
        }
        break;

    case iiEQU:         /*------------- Atom Equivalence ------------*/
        switch( bType ) {
        case 0:         /* equivalence main = ...*/
            r = "??";      /* should not happen */
            break;
        case iitISO:     /* equivalence main iso = ...*/
            switch( bEq2 ) {
            case 0:
                r = "m";       /* equivalence main = main equ  */
                break;
            default:
                r = "??";      /* should not happen */
                break;
            }
            break;
        case iitNONTAUT: /* equivalence non-taut = ... */
            switch( bEq2 ) {
            case 0:
                r = "m";       /* equivalence non-taut = main equ */
                break;
            default:
                r = "??";      /* should not happen */
                break;
            }
            break;
        case (iitNONTAUT | iitISO): /*  equivalence non-taut iso = ... */
            switch( bEq2 ) {
            case 0:
                r = "m";       /* equivalence non-taut iso = main equ */
                break;
            case iiEq2ISO:
                r = "M";       /* equivalence non-taut iso = main iso equ */
                break;
            case iiEq2NONTAUT:
                r = "n";       /* equivalence non-taut iso = non-taut equ */
                break;
            default:
                r = "??";      /* should not happen */
                break;
            }
            break;
        default:
            r = "??";          /* should not happen */
            break;
        }
        break;
    default:
        r = "??";      /* should not happen */
        break;
    }
    return r;
}

/**********************************************************************************************/

#define OUT_NONTAUT  OUT_NN  /* was OUT_NT until 2004-04-07 */


/**********************************************************************************************/
int OutputINChI2(char *pStr, int nStrLen,
                 INCHI_SORT *pINChISortTautAndNonTaut2[][TAUT_NUM],
                 int iINChI,
                 ORIG_STRUCT *pOrigStruct,
                 int bDisconnectedCoord, int bOutputType, int bINChIOutputOptions,
                 int bXml, int bAbcNumbers,int bCtPredecessors, int bNoStructLabels,
                 int num_components2[],
                 int num_non_taut2[], int num_taut2[],
                 INCHI_IOSTREAM *output_file, INCHI_IOSTREAM *log_file,
                 int num_input_struct,
                 const char *szSdfLabel, const char *szSdfValue, long lSdfId,
                 int *pSortPrintINChIFlags,
                 unsigned char save_opt_bits)
{
    int bINChIOutputOptions0 = bINChIOutputOptions & ~(INCHI_OUT_XML | INCHI_OUT_PLAIN_TEXT | INCHI_OUT_PLAIN_TEXT_COMMENTS);
    int bINChIOutputOptionsCur;
    int bCurOption, ret, i;

    ret = 0;

    for ( i = 0; i < 3; i ++ )
    {
        switch( i )
        {
        case 0:
            bCurOption = INCHI_OUT_XML;
            break;
        case 1:
            bCurOption = INCHI_OUT_PLAIN_TEXT;
            break;
        case 2:
            bCurOption = INCHI_OUT_PLAIN_TEXT_COMMENTS;
            break;
        default:
            continue;
        }
        if ( bINChIOutputOptions & bCurOption )
        {
            bINChIOutputOptionsCur = bINChIOutputOptions0 | bCurOption;
            if ( i != 1 )
            {
                bINChIOutputOptionsCur  &= ~INCHI_OUT_TABBED_OUTPUT;
            }
            ret |= OutputINChI1( pStr, nStrLen,
                                 pINChISortTautAndNonTaut2,
                                 iINChI,
                                 pOrigStruct,
                                 bDisconnectedCoord, bOutputType, bINChIOutputOptionsCur,
                                 bXml, bAbcNumbers,bCtPredecessors, bNoStructLabels,
                                 num_components2,
                                 num_non_taut2, num_taut2,
                                 output_file, log_file,
                                 num_input_struct,
                                 szSdfLabel, szSdfValue, lSdfId,
                                 pSortPrintINChIFlags,
                                 save_opt_bits);
        }
    }

    return ret;
}

/**********************************************************************************/
char *szGetTag( const INCHI_TAG *Tag, int nTag, int bTag, char *szTag, int *bAlways )
{
    int i, j, bit, num, len;
    if ( 0 < nTag && nTag < 3 ) {
        /* no plain text comments: pick up the last tag */
        for ( i = 0, j = -1, bit = 1; i < MAX_TAG_NUM; i ++, bit <<= 1 ) {
            if ( bTag & bit ) {
                j = i;
            }
        }
        if ( j >= 0 ) {
            strcpy( szTag, nTag == 1? Tag[j].szXmlLabel : nTag == 2? Tag[j].szPlainLabel : "???" );
            if ( nTag != 2 ) {
                *bAlways = Tag[j].bAlwaysOutput;
            }
            return szTag;
        }
    } else
    if ( nTag == 3 ) {
        /* plain text with comments */
        szTag[0] = '{';
        szTag[1] = '\0';
        for ( i = 0, j = -1, bit = 1, num=0; i < MAX_TAG_NUM; i ++, bit <<= 1 ) {
            if ( bTag & bit ) {
                j = i;
                if ( num ++ ) {
                    strcat( szTag, ":" );
                }
                strcat( szTag, Tag[i].szPlainComment );
            }
        }
        if ( num ) {
            strcat( szTag, "}" );
            num = strlen( Tag[j].szPlainLabel );
            len = strlen( szTag );
            if ( len ) {
                memmove( szTag + num, szTag, len+1 );
                memcpy( szTag, Tag[j].szPlainLabel, num );
            } else {
                strcpy ( szTag, Tag[j].szPlainLabel );
            }
            *bAlways = Tag[j].bAlwaysOutput;
        } else {
            strcpy( szTag, "???" );
        }
        return szTag;
    }
    strcpy( szTag, "???" );
    return szTag;
}


/***************************************************************************************/
/*  sorting in descending order: return -1 if *p1 > *p2, return +1 if *p1 < *p2               */
/***************************************************************************************/
int OutputINChI1(char *pStr, int nStrLen,
                 INCHI_SORT *pINChISortTautAndNonTaut2[][TAUT_NUM],
                 int iINChI,
                 ORIG_STRUCT *pOrigStruct,
                 int bDisconnectedCoord, int bOutputType, int bINChIOutputOptions,
                 int bXml, int bAbcNumbers,int bCtPredecessors, int bNoStructLabels,
                 int num_components2[], int num_non_taut2[], int num_taut2[],
                 INCHI_IOSTREAM *output_file, INCHI_IOSTREAM *log_file,
                 int num_input_struct,
                 const char *szSdfLabel, const char *szSdfValue, long lSdfId,
                 int *pSortPrintINChIFlags,
                 unsigned char save_opt_bits)
{
/*
  bINChIOutputOptions bits:

    INCHI_OUT_NO_AUX_INFO           0x0001    do not output Aux Info
    INCHI_OUT_SHORT_AUX_INFO        0x0002    output short version of Aux Info
    INCHI_OUT_ONLY_AUX_INFO         0x0004    output only Aux Info
    INCHI_OUT_EMBED_REC             0x0008    embed reconnected INChI into disconnected INChI

*/

    /*int ATOM_MODE = ((bAbcNumbers?2:0)|5|(bCtPredecessors?8:0));*/
    int ATOM_MODE = ((bAbcNumbers?CT_MODE_ABC_NUMBERS:0)
                    | CT_MODE_ATOM_COUNTS
                    | CT_MODE_NO_ORPHANS
#if ( EQL_H_NUM_TOGETHER == 1 )
                    | CT_MODE_EQL_H_TOGETHER
#endif
#if ( ABC_CT_NUM_CLOSURES == 1 )
                    | (bAbcNumbers && bCtPredecessors? CT_MODE_ABC_NUM_CLOSURES:0)
#endif
                    | (bCtPredecessors?CT_MODE_PREDECESSORS:0));

    int TAUT_MODE = (bAbcNumbers?CT_MODE_ABC_NUMBERS:0);
    char sDifSegs[DIFL_LENGTH][DIFS_LENGTH];
    /* bOutputType =
         TAUT_YES  => tautomeric only (if no tautomeric components then no output;
         TAUT_NON  => only non-tautomeric output (if no non-taut present then no output;
         TAUT_BOTH => tautomeric and non-tautomeric */

    int  i, j, ii, jj, /*ii2, jj2,*/ tot_len, tot_len2, bOverflow, bEmbeddedOutputCalled=0;
    int  bIsotopic, bTautIsoHNum, bTautIsoAt, bHasIsotopicAtoms[TAUT_NUM];
    int  bStereoSp2[TAUT_NUM], bStereoSp3[TAUT_NUM];
    int  bIsotopicStereoSp2[TAUT_NUM], bIsotopicStereoSp3[TAUT_NUM];
    int  bStereoAbsInverted[TAUT_NUM], bIsotopicStereoAbsInverted[TAUT_NUM];
    int  bStereoAbs[TAUT_NUM], bIsotopicStereoAbs[TAUT_NUM];
    int  bAtomEqu[TAUT_NUM], bTautEqu[TAUT_NUM], bIsotopicAtomEqu[TAUT_NUM], bIsotopicTautEqu[TAUT_NUM];
    int  bInvStereo[TAUT_NUM], bInvIsotopicStereo[TAUT_NUM];
    int  bInvStereoOrigNumb[TAUT_NUM], bInvIsotopicStereoOrigNumb[TAUT_NUM], bIsotopicOrigNumb[TAUT_NUM];
    int  bTautomeric, bNonTautomeric, bTautomericAcid, bHardAddRemProton, iCurTautMode;
    int  bRequestedRacemicStereo=0, bRequestedRelativeStereo = 0, bRelRac;
    int  bRacemicStereo[TAUT_NUM], bRelativeStereo[TAUT_NUM];
    int  bIsotopicRacemicStereo[TAUT_NUM], bIsotopicRelativeStereo[TAUT_NUM];
    int  bChargesRadVal[TAUT_NUM], bOrigCoord[TAUT_NUM];
    int  bIgn_UU_Sp3[TAUT_NUM], bIgn_UU_Sp2[TAUT_NUM];
    int  bIgn_UU_Sp3_Iso[TAUT_NUM], bIgn_UU_Sp2_Iso[TAUT_NUM];
    int  ind, inc, bNonTautIsIdenticalToTaut = 1;
    int  bNonTautNonIsoIdentifierNotEmpty = 0, bNonTautIsoIdentifierNotEmpty = 0;
    INCHI_SORT   **pINChISortTautAndNonTaut = pINChISortTautAndNonTaut2[iINChI];
    INCHI_SORT   *pINChISort =pINChISortTautAndNonTaut[TAUT_YES];
    INCHI_SORT   *pINChISort2=pINChISortTautAndNonTaut[TAUT_YES];
    INCHI_SORT   *is, *is2;
    INChI        *pINChI /*, *pINChI2*/;
    INChI_Aux    *pINChI_Aux = NULL;


    int  ret = 0; /*  0=>failed, 1=>success */
    int  bOutType = bOutputType; /* ??? */
    int  nTag;
    int  bTautomericOutputAllowed, bSecondNonTautPass;
    int  num_components = num_components2[iINChI];
    int  num_comp[TAUT_NUM], max_num_comp;
    int  num_iso_H[NUM_H_ISOTOPES], bHasIsoH;
    int  nNumRemovedProtons, nNumMovedProtons;
    int  bTautAndNonTaut, bTautIsNonTaut;

    int  bAlways           = 0;
    int  bUseMulipliers    = 1;
    int  bOmitRepetitions  = 1;
    int  bPlainTextTags    = 2;  /* 0 => no plain tags, 1=> plain text tags, 2=>plaintext tags without consecutive // */
    int  bPlainText        = 0 != (bINChIOutputOptions & (INCHI_OUT_PLAIN_TEXT | INCHI_OUT_PLAIN_TEXT_COMMENTS));
    int  bPlainTextCommnts = 0 != (bINChIOutputOptions & INCHI_OUT_PLAIN_TEXT_COMMENTS);
    int  bPlainTabbedOutput;
    int  bTag1, bTag2, bTag3, bFhTag; /* tag bits */
    int  nCurINChISegment, nSegmAction;
    char szTag1[MAX_TAG_LEN], szTag2[MAX_TAG_LEN], szTag3[MAX_TAG_LEN];
    const char *pLF, *pTAB;

    /*^^^ 15 April, 2008 */
    int bFixTranspChargeBug = 0;
#if ( FIX_TRANSPOSITION_CHARGE_BUG == 1 ) /* 2008-01-02 */
    if ( INCHI_OUT_FIX_TRANSPOSITION_CHARGE_BUG & bINChIOutputOptions )
        bFixTranspChargeBug = 1;
#endif
    /*^^^ 15 April, 2008 */

    bXml                   = 0 != (bINChIOutputOptions & INCHI_OUT_XML);
    nTag  = bPlainTextCommnts? 3 : bPlainText? 2 : bXml? 1 : 0; /* tag type */
    ind                    = bXml? 1 : -1;
    inc                    = bXml? 1 : -1;
    pLF                    = bPlainTextCommnts? "\n" : "\0";
    bFhTag                 = 0;
    bPlainTabbedOutput     = 0 != (bINChIOutputOptions & INCHI_OUT_TABBED_OUTPUT) &&
                             bPlainText && !bXml && !bPlainTextCommnts;
#if ( !defined(TARGET_API_LIB) && !defined(TARGET_LIB_FOR_WINCHI) )
    pTAB                   = bPlainTabbedOutput? "\t" : "\n";
#else
    pTAB                   = "\n";
#endif


    memset( sDifSegs, DIFV_BOTH_EMPTY, sizeof(sDifSegs) );

    if ( !pStr ) {
        inchi_ios_eprint(log_file,
            "Cannot allocate output buffer. No output for structure #%d.%s%s%s%s\n",
            num_input_struct, SDF_LBL_VAL(szSdfLabel, szSdfValue));
        return ret;
    }

    bSecondNonTautPass = 0;
/* -- commented out to allow empty InChI --
    if (!num_components )
    {
        return 0;
    }
*/

    /* init version string */
    if ( !VER_STRING[0] )
    {
        strcpy(VER_STRING,  "(V");
        strcat(VER_STRING,  INCHI_VERSION);
        strcat(VER_STRING,  ")");
    }
    for ( i = 0; i < TAUT_NUM; i ++ )
    {
        bHasIsotopicAtoms[i]      = num_comp[i]                  =
        bStereoSp2[i]             = bStereoSp3[i]                =
        bIsotopicStereoSp2[i]     = bIsotopicStereoSp3[i]        =
        bIsotopicOrigNumb[i]      =
        bStereoAbs[i]             = bIsotopicStereoAbs[i]         =
        bStereoAbsInverted[i]     = bIsotopicStereoAbsInverted[i] =
        bRacemicStereo[i]         = bRelativeStereo[i]            =
        bIsotopicRacemicStereo[i] = bIsotopicRelativeStereo[i]    =
        bAtomEqu[i]               = bTautEqu[i]                   =
        bIsotopicAtomEqu[i]       = bIsotopicTautEqu[i]           =
        bInvStereo[i]             = bInvIsotopicStereo[i]         =
        bInvStereoOrigNumb[i]     = bInvIsotopicStereoOrigNumb[i] =
        bIgn_UU_Sp3[i]            = bIgn_UU_Sp2[i]                =
        bIgn_UU_Sp3_Iso[i]        = bIgn_UU_Sp2_Iso[i]            =
        bChargesRadVal[i]         = bOrigCoord[i]                 = 0;
    }

    /*  find if it is isotopic */
    bIsotopic       = bTautomeric = bNonTautomeric = bTautomericAcid =
                      bHardAddRemProton = bTautIsoHNum = bTautIsoAt = 0;
    bTautAndNonTaut = bTautIsNonTaut = 0;
    /*
         x = bStereo, bStereoSp2, bStereoSp3, bStereoAbsInverted,
             bIsotopicStereo, bIsotopicStereoSp2, bIsotopicStereoSp3, bIsotopicStereoAbsInverted

         OUT_N1: x[TAUT_NON] refers to non-tautomeric only
         OUT_T1: x[TAUT_YES] refers to tautomeric if exists otherwise non-tautomeric
         OUT_NT: x[TAUT_NON] refers to non-taut representations of tautomeric
         OUT_TN: x[TAUT_YES] refers to tautomeric if exists otherwise non-tautomeric
                 x[TAUT_NON] refers to non-taut representations of tautomeric
     */

    memset( num_iso_H, 0, sizeof(num_iso_H) );
    nNumRemovedProtons = 0;
    nNumMovedProtons   = 0;
    bHasIsoH           = 0;
    bTautomericOutputAllowed = (bOutType==OUT_T1 || bOutType== OUT_TN);
    pINChISort=pINChISortTautAndNonTaut[bTautomericOutputAllowed? TAUT_YES : TAUT_NON];
    is  = pINChISort;
    is2 = (bOutType== OUT_TN)? pINChISortTautAndNonTaut[TAUT_NON] : NULL;

    for ( i = 0, is2 = pINChISortTautAndNonTaut[TAUT_NON]; i < num_components; i ++, is ++, is2? is2++:NULL )
    {
        CompINChILayers( is, is2, sDifSegs, bFixTranspChargeBug );
        bNonTautIsIdenticalToTaut = bNonTautIsIdenticalToTaut && !CompINChITautVsNonTaut(is, is2, 1);
        if ( is && (pINChI_Aux = is->pINChI_Aux[TAUT_YES]) )
        {
            for ( j = 0; j < NUM_H_ISOTOPES; j ++ )
            {
                bHasIsoH     += abs(pINChI_Aux->nNumRemovedIsotopicH[j]);
                num_iso_H[j] += pINChI_Aux->nNumRemovedIsotopicH[j];
            }
            nNumRemovedProtons += pINChI_Aux->nNumRemovedProtons;
            nNumMovedProtons   += abs(pINChI_Aux->nNumRemovedProtons);
        }
        if ( bTautomericOutputAllowed )
        {
            /* check for removed isotopic H */
            for ( j = TAUT_YES; j < TAUT_NUM; j ++ )
            {
                switch ( bOutType ) {
                case OUT_N1: /* x[TAUT_NON]: non-tautomeric only -- never happens */
                    jj = GET_II(bOutType,is);
                    if ( jj != j )
                        continue;
                    ii = TAUT_NON;
                    break;
                case OUT_T1: /* x[TAUT_YES]: tautomeric if present otherwise non-tautomeric */
                    jj = GET_II(bOutType,is);
                    if ( jj != j )
                        continue;
                    ii = TAUT_YES;
                    break;
                case OUT_NT: /* x[TAUT_NON]: only non-taut representations of tautomeric -- never happens */
                    jj = GET_II(bOutType,is);
                    if ( jj != j )
                        continue;
                    ii = TAUT_NON;
                    break;
                /* main path of control flow */
                case OUT_TN: /* x[TAUT_YES]: tautomeric if present otherwise non-tautomeric;
                              * x[TAUT_NON]: non-taut only if tautomeric is present */
                    jj = ( j == TAUT_YES )? GET_II(OUT_T1,is) : ( j == TAUT_NON )? GET_II(OUT_NT,is) : -1;
                    if ( jj == TAUT_YES )
                    {
                        /* Fix12 */
                        if ( is->pINChI[jj]->lenTautomer > 0 )
                        {
                            bTautAndNonTaut += (!is->pINChI[jj]->bDeleted && HAS_N(is));
                        } else
                        {
                            bTautIsNonTaut ++;
                        }
                    }
                    if ( jj < 0 )
                        continue;
                    ii = j;
                    break;
                default:
                    continue;
                }
                if ( jj != j )
                    continue;
                if ( (pINChI = is->pINChI[jj]) && pINChI->nNumberOfAtoms > 0 && (pINChI_Aux = is->pINChI_Aux[jj]) )
                {
                    bTautIsoHNum += (pINChI_Aux->nNumRemovedIsotopicH[0] +
                                     pINChI_Aux->nNumRemovedIsotopicH[1] +
                                     pINChI_Aux->nNumRemovedIsotopicH[2]);
                    bTautIsoAt   += (pINChI->nNumberOfIsotopicAtoms>0 || pINChI->nNumberOfIsotopicTGroups > 0 );
                }
            }
        }
    }
    sDifSegs[DIFL_M ][DIFS_p_PROTONS] = nNumRemovedProtons? DIFV_NEQ2PRECED : DIFV_BOTH_EMPTY;
    sDifSegs[DIFL_MI][DIFS_h_H_ATOMS] = bHasIsoH?           DIFV_NEQ2PRECED : DIFV_BOTH_EMPTY;

    MarkUnusedAndEmptyLayers( sDifSegs );



    bNonTautIsIdenticalToTaut = bNonTautIsIdenticalToTaut && !bTautIsoHNum;
    /*********************************************************************************************/
    for ( i = 0, is = pINChISort; i < num_components; i ++, is ++ )
    {
        int bCurIso, bCurStereo, bCurIsoStereo, bCurHasIsoStereo /* Fix14 */, bCurTaut /*, bCurTaut2*/;
        int bCompExists, bCurIsoHPos, bCurIsoHStereo;
        int bCurStereoSp2, bCurIsoStereoSp2, bCurStereoSp3, bCurIsoStereoSp3, bCurIsoStereoSp3Inv;
        int bCurRacemic, bCurRelative, bCurIsoRacemic, bCurIsoRelative;
        bCompExists = 0;
        for ( j = TAUT_NON; j < TAUT_NUM; j ++ )
        {
            switch ( bOutType ) {
            case OUT_N1: /* x[TAUT_NON]: non-tautomeric only */
                jj = GET_II(bOutType,is);
                if ( jj != j )
                    continue;
                ii = TAUT_NON;
                break;
            case OUT_T1: /* x[TAUT_YES]: tautomeric if present otherwise non-tautomeric */
                jj = GET_II(bOutType,is);
                if ( jj != j )
                    continue;
                ii = TAUT_YES;
                break;
            case OUT_NT: /* x[TAUT_NON]: only non-taut representations of tautomeric */
                jj = GET_II(bOutType,is);
                if ( jj != j )
                    continue;
                ii = TAUT_NON;
                break;
            /* main control flow comes here: requested both mobile and fixed H results */
            case OUT_TN: /* x[TAUT_YES]: tautomeric if present otherwise non-tautomeric;
                          * x[TAUT_NON]: non-taut only if tautomeric is present */
                jj = ( j == TAUT_YES )? GET_II(OUT_T1,is) : ( j == TAUT_NON )? GET_II(OUT_NT,is) : -1;
                if ( jj < 0 )
                {
                    /* Fix12 */
                    if ( bTautAndNonTaut && bTautIsNonTaut &&
                         j == TAUT_NON && 0 <= (jj = GET_II(OUT_T1,is)) &&
                         !is->pINChI[jj]->bDeleted && !is->pINChI[jj]->lenTautomer )
                    {
                        ; /* the requested non-tautomeric component is in tautomeric position
                             (is->pINChI[TAUT_YES]);
                             process it also as non-tautomeric if Fixed-H layer was requested */
                    }
                    else
                    {
                        continue;
                    }
                }
                ii = j; /* ii is what we wanted; jj is what we found (0 = TAUT_NON: fixed_H, 1 = TAUT_YES: mobile_H) */
                /* -- not used 2004-09-16 ---
                if ( is2 ) {
                    jj2 = ( j == TAUT_YES )? GET_II(OUT_T1,is2) : ( j == TAUT_NON )? GET_II(OUT_NT,is2) : -1;
                    if ( jj2 >= 0 ) {
                        ii2 = j;
                    } else {
                        ii2 = -1;
                    }
                } else {
                    jj2 = ii2 = -1;
                }
                -----------------------------*/
                break;
            default:
                continue;
            }
            if ( (pINChI = is->pINChI[jj]) && pINChI->nNumberOfAtoms > 0 )
            {
                /*pINChI_Aux = is->pINChI_Aux[jj];*/
                bCompExists ++;
                bCurTaut            = (pINChI->lenTautomer > 0);
                bCurIso             = (pINChI->nNumberOfIsotopicAtoms>0 || pINChI->nNumberOfIsotopicTGroups > 0 );
                bCurIsoHPos         = ((pINChI->nPossibleLocationsOfIsotopicH && pINChI->nPossibleLocationsOfIsotopicH[0] > 1) || pINChI->lenTautomer > 1);
                /* present isotopic H + their possible positions AND/OR isotopic atoms */
                bCurIsoHStereo      = (bCurIsoHPos && (bTautIsoHNum || bTautIsoAt)) || bCurIso;
                if ( jj == j && pINChI->bDeleted )
                {
                    num_comp[j] --;
                    if ( bCurTaut )
                    {
                        bTautomeric        |= 1; /* tautomeric representation is present */
                        bNonTautomeric     |= HAS_N(is);
                    }
                    bIsotopic              |= bCurIso;
                    continue; /* deleted H(+) in tautomeric representation */
                }
                bCurStereoSp2       = pINChI->Stereo && (pINChI->Stereo->nNumberOfStereoBonds > 0);
                bCurHasIsoStereo    =
                bCurStereoSp3       = pINChI->Stereo && (pINChI->Stereo->nNumberOfStereoCenters > 0 );
                bCurIsoStereoSp2    = bCurIsoHStereo && pINChI->StereoIsotopic && (pINChI->StereoIsotopic->nNumberOfStereoBonds > 0);
                bCurIsoStereoSp3    = bCurIsoHStereo && pINChI->StereoIsotopic && (pINChI->StereoIsotopic->nNumberOfStereoCenters > 0);
                bCurIsoStereoSp3Inv = bCurIsoStereoSp3 && pINChI->StereoIsotopic->nCompInv2Abs; /* inversion changes sp3 stereo */
                bRequestedRacemicStereo         |= (0!=(pINChI->nFlags & INCHI_FLAG_RAC_STEREO));

                bRequestedRelativeStereo        |= (0!=(pINChI->nFlags & INCHI_FLAG_REL_STEREO));
                /* check whether isotopic stereo is same as non-isotopic; if same than do not output isotopic stereo */
                if ( bCurStereoSp2 && bCurIsoStereoSp2 )
                {
                    bCurIsoStereoSp2 = !Eql_INChI_Stereo( pINChI->Stereo, EQL_SP2, pINChI->StereoIsotopic, EQL_SP2, 0 );
                }
                if ( bCurStereoSp3 && bCurIsoStereoSp3 )
                {
                    /* bCurIsoStereoSp3=0 means (iso stereo sp3) = (non-iso stereo sp3) or (iso stereo sp3) = Inv(non-iso stereo sp3) */
                    bCurIsoStereoSp3 = !Eql_INChI_Stereo( pINChI->Stereo, EQL_SP3, pINChI->StereoIsotopic, EQL_SP3,
                              (pINChI->nFlags & INCHI_FLAG_RAC_STEREO) || (pINChI->nFlags & INCHI_FLAG_REL_STEREO) );
                    if ( !bCurIsoStereoSp3 ) {
                        /* inversion changes iso sp3 differently from non-iso sp3 Fix11 */
                        bCurIsoStereoSp3Inv &= (pINChI->StereoIsotopic->nCompInv2Abs != pINChI->Stereo->nCompInv2Abs);
                    }
                }

                bCurRelative        =  bRequestedRelativeStereo && bCurStereoSp3;
#if ( REL_RAC_STEREO_IGN_1_SC == 1 )
                bCurRelative        =  bCurRelative &&
                                      (pINChI->Stereo->nNumberOfStereoCenters > 1 ) &&
                                      (pINChI->Stereo->nCompInv2Abs != 0) &&
#endif



                bCurIsoRelative     = bRequestedRelativeStereo && (bCurIsoStereoSp3 || bCurIsoStereoSp3Inv);
#if ( REL_RAC_STEREO_IGN_1_SC == 1 )
                bCurIsoRelative     = bCurIsoRelative &&
                                      (pINChI->StereoIsotopic->nNumberOfStereoCenters > 1 ) &&
                                      (pINChI->StereoIsotopic->nCompInv2Abs != 0) &&
#endif


                bCurRacemic         = bRequestedRacemicStereo && bCurStereoSp3;
#if ( REL_RAC_STEREO_IGN_1_SC == 1 )
                bCurRacemic         = bCurRacemic &&
                                      (pINChI->Stereo->nCompInv2Abs != 0) &&
                                      (pINChI->Stereo->nNumberOfStereoCenters > 0 ) ?
                                      pINChI->Stereo->nNumberOfStereoCenters : 0;
#endif

                bCurIsoRacemic      = bRequestedRacemicStereo && (bCurIsoStereoSp3 || bCurIsoStereoSp3Inv);
#if ( REL_RAC_STEREO_IGN_1_SC == 1 )
                bCurIsoRacemic      = bCurIsoRacemic &
                                      (pINChI->StereoIsotopic->nCompInv2Abs != 0) &&
                                      (pINChI->StereoIsotopic->nNumberOfStereoCenters > 0 ) ?
                                       pINChI->StereoIsotopic->nNumberOfStereoCenters : 0;
#endif
                if ( bRequestedRelativeStereo )
                {
                    bCurStereoSp3     = bCurRelative || (bCurStereoSp3 && (pINChI->Stereo->nNumberOfStereoCenters > 1 )); /* Fix11 */
                    bCurIsoStereoSp3  = bCurIsoRelative   ? bCurIsoStereoSp3 : 0;
                }
                else
                {
                    if ( bRequestedRacemicStereo )
                    {
                        bCurStereoSp3     = bCurRacemic    > 1 || (bCurStereoSp3 && (pINChI->Stereo->nNumberOfStereoCenters > 1 )); /* Fix11 */
                        bCurIsoStereoSp3  = bCurIsoRacemic > 1? bCurIsoStereoSp3 : 0;
                    }
                }
                bCurStereo          = bCurStereoSp2    || bCurStereoSp3;
                bCurIsoStereo       = bCurIsoStereoSp2 || bCurIsoStereoSp3;

                bIsotopic              |= bCurIso;
                bHasIsotopicAtoms[ii]  |= bCurIso;
                bStereoSp2[ii]         |= bCurStereoSp2;
                bStereoSp3[ii]         |= bCurStereoSp3;
                bIgn_UU_Sp3[ii]        |= !bCurStereoSp3 && (pINChI->nFlags & INCHI_FLAG_SC_IGN_ALL_UU);
                bIgn_UU_Sp2[ii]        |= !bCurStereoSp2 && (pINChI->nFlags & INCHI_FLAG_SB_IGN_ALL_UU);
                bIsotopicStereoSp2[ii] |= bCurIsoStereoSp2;
                bIsotopicStereoSp3[ii] |= bCurIsoStereoSp3;
                bIgn_UU_Sp3_Iso[ii]    |= !bCurIsoStereoSp3 && (pINChI->nFlags & INCHI_FLAG_SC_IGN_ALL_ISO_UU);
                bIgn_UU_Sp2_Iso[ii]    |= !bCurIsoStereoSp2 && (pINChI->nFlags & INCHI_FLAG_SB_IGN_ALL_ISO_UU);
                bStereoAbs[ii]                  |= bCurStereoSp3 && (pINChI->Stereo->nCompInv2Abs != 0);
                bStereoAbsInverted[ii]          |= bCurStereoSp3 && (pINChI->Stereo->nCompInv2Abs < 0);
                /* Fix08: missing isotopic inverted flag if isotopic = inverted non-isotopic */
                bIsotopicStereoAbsInverted[ii]  |= (bCurIsoStereoSp3 && (pINChI->StereoIsotopic->nCompInv2Abs < 0)) ||
                                                   (!bCurIsoStereoSp3  && pINChI->StereoIsotopic  && pINChI->Stereo &&
                                                   pINChI->StereoIsotopic->nCompInv2Abs &&
                                                   pINChI->StereoIsotopic->nCompInv2Abs != pINChI->Stereo->nCompInv2Abs);
                /* Fix 11: missing /s1 if only isotopic stereo is inverted */
                bIsotopicStereoAbs[ii]          |= (bCurIsoStereoSp3 && (pINChI->StereoIsotopic->nCompInv2Abs != 0)) ||
                                                   (!bCurIsoStereoSp3  && pINChI->StereoIsotopic  && pINChI->Stereo &&
                                                   pINChI->StereoIsotopic->nCompInv2Abs &&
                                                   pINChI->StereoIsotopic->nCompInv2Abs != pINChI->Stereo->nCompInv2Abs);

                bRelativeStereo[ii]             |= bCurRelative;
                bIsotopicRelativeStereo[ii]     |= bCurIsoRelative;
                bRacemicStereo[ii]              |= bCurRacemic;
                bIsotopicRacemicStereo[ii]      |= bCurIsoRacemic;

                bTautomericAcid                 |= (0!=(pINChI->nFlags & INCHI_FLAG_ACID_TAUT));
                bHardAddRemProton               |= (0!=(pINChI->nFlags & INCHI_FLAG_HARD_ADD_REM_PROTON));
                if ( bCurTaut )
                {
                    bTautomeric        |= 1; /* tautomeric representation is present */
                    /* does tautomeric structure have also a non-tautomeric repesentation? */
                    bNonTautomeric     |= HAS_N(is);
                }

                /* auxiliary info */
                if ( !(bINChIOutputOptions & INCHI_OUT_NO_AUX_INFO) && (pINChI_Aux = is->pINChI_Aux[jj]) )
                {
                    /* detect presence of constitutional equivalence onfo */
                    int bCurEqu, bCurTautEqu=0, bCurIsoEqu=0, bCurIsoTautEqu=0; /* Fix15-disabled */
                    bAtomEqu[ii] |= (bCurEqu = bHasEquString( pINChI_Aux->nConstitEquNumbers,
                                                   pINChI_Aux->nNumberOfAtoms));
                    if ( bCurTaut )
                    {
                        bTautEqu[ii] |= (bCurTautEqu = bHasEquString( pINChI_Aux->nConstitEquTGroupNumbers,
                                                       pINChI_Aux->nNumberOfTGroups));
                    }
                    if ( bCurIso )
                    {
                        bIsotopicAtomEqu[ii] |= (bCurIsoEqu = bHasEquString( pINChI_Aux->nConstitEquIsotopicNumbers,
                                                               pINChI_Aux->nNumberOfAtoms)) /*|| bCurEqu*/;
                        if ( bCurTaut )
                        {
                            bIsotopicTautEqu[ii] |= (bCurIsoTautEqu = bHasEquString( pINChI_Aux->nConstitEquIsotopicTGroupNumbers,
                                                                   pINChI_Aux->nNumberOfTGroups)) /*|| bCurTautEqu*/;
                        }
                        /* non-zero if isotopic numbering for inverted isotopic stereo is different */
                        bIsotopicOrigNumb[ii] |= bCurHasIsoStereo && /* Fix14 */
                                                 pINChI_Aux->nOrigAtNosInCanonOrdInv &&
                                                 pINChI_Aux->nIsotopicOrigAtNosInCanonOrd &&
                            (0 != memcmp( pINChI_Aux->nOrigAtNosInCanonOrdInv,
                                          pINChI_Aux->nIsotopicOrigAtNosInCanonOrd,
                                          sizeof(pINChI_Aux->nOrigAtNosInCanonOrdInv[0])
                                          * pINChI_Aux->nNumberOfAtoms));

                    }
                    /* inverted stereo */
                    if ( bCurStereoSp3 && pINChI->Stereo->nCompInv2Abs )
                    {
                        bInvStereo[ii]         |= 1;
                        bInvStereoOrigNumb[ii] |= pINChI_Aux->nOrigAtNosInCanonOrd &&
                                                  pINChI_Aux->nOrigAtNosInCanonOrdInv &&
                            (0 != memcmp( pINChI_Aux->nOrigAtNosInCanonOrd,
                                          pINChI_Aux->nOrigAtNosInCanonOrdInv,
                                          sizeof(pINChI_Aux->nOrigAtNosInCanonOrd[0])
                                          * pINChI_Aux->nNumberOfAtoms));
                    }
                    /* inverted isotopic stereo */
                    if ( bCurIsoStereoSp3 && pINChI->StereoIsotopic->nCompInv2Abs )
                    {
                        bInvIsotopicStereo[ii]         |= 1;
                        bInvIsotopicStereoOrigNumb[ii] |= pINChI_Aux->nIsotopicOrigAtNosInCanonOrd &&
                                                          pINChI_Aux->nIsotopicOrigAtNosInCanonOrdInv &&
                            (0 != memcmp( pINChI_Aux->nIsotopicOrigAtNosInCanonOrd,
                                          pINChI_Aux->nIsotopicOrigAtNosInCanonOrdInv,
                                          sizeof(pINChI_Aux->nIsotopicOrigAtNosInCanonOrd[0])
                                          * pINChI_Aux->nNumberOfAtoms));
                    }
                    if ( pINChI_Aux->OrigInfo && bHasOrigInfo(pINChI_Aux->OrigInfo, pINChI_Aux->nNumberOfAtoms) )
                    {
                        bChargesRadVal[ii] |= 1;
                    }
                }
            }
        }
        if ( bCompExists )
        {
            for ( j = TAUT_NON; j < TAUT_NUM; j ++ )
            {
                num_comp[j] ++;
            }
        }
    }
    if ( bTautomeric /*&& bTautomericAcid*/ ) /* "&& bTautomericAcid" commented out 2004-06-02 */
    {
        bTautomeric += bTautomericAcid; /* long-range tautomerism */
        bTautomeric += (bHardAddRemProton? 4 : 0);
    }
    if ( bRequestedRacemicStereo || bRequestedRelativeStereo )
    {
        /* do not output inverted stereo info */
        for ( i = 0; i < TAUT_NUM; i ++ )
        {
            /* Fix11 */
            bStereoAbsInverted[i] =
            bStereoAbs[i]         =
            bInvStereo[i]         =
            bInvStereoOrigNumb[i] =  0;
            /* bIsotopicRelativeStereo[i]=0 may happen because iso stereo is same or inverted non-iso stereo */
            bIsotopicStereoAbsInverted[i] =
            bIsotopicStereoAbs[i]         =
            bInvIsotopicStereo[i]         =
            bInvIsotopicStereoOrigNumb[i] = 0;
            /* -- commented out: Fix11--
            if ( bRacemicStereo[i] || bRelativeStereo[i] )
            {
                bStereoAbsInverted[i] =
                bStereoAbs[i]         =
                bInvStereo[i]         =
                bInvStereoOrigNumb[i] =  0;
            }
            if ( bIsotopicRacemicStereo[i] || bIsotopicRelativeStereo[i] )
            {
                bIsotopicStereoAbsInverted[i] =
                bIsotopicStereoAbs[i]         =
                bInvIsotopicStereo[i]         =
                bInvIsotopicStereoOrigNumb[i] = 0;
            }
            */
        }
    }


    iCurTautMode = bOutType == OUT_N1? TAUT_NON:  /* only non-taut */
                   bOutType == OUT_T1? TAUT_YES:  /* tautomeric if present, otherwise non-tautomeric */
                   bOutType == OUT_NT? TAUT_NON:  /* only non-taut representations of tautomeric */
                   bOutType == OUT_TN? TAUT_YES:  /* tautomeric if present otherwise non-tautomeric; */
                                             -1;   /* separately output non-taut representations of tautomeric if present */

    if ( iCurTautMode < 0 )
    {
        return 0;  /* error */
    }

    if ( bXml )
    {
        ind += inc* (1+iINChI);
    }

    bOverflow = 0;

    num_components = num_comp[iCurTautMode];

    max_num_comp   = inchi_max(num_comp[TAUT_NON], num_comp[TAUT_YES]);

    if ( bINChIOutputOptions & INCHI_OUT_ONLY_AUX_INFO )
    {
        goto output_aux_info;
    }

    nCurINChISegment = DIFL_M;

    /******************************************
     *
     *  Structure (Compound) Header
     *
     ******************************************/
    if ( bXml )
    {
        /* -- moved to the line above goto output_aux_info;
        ind += inc* (1+iINChI);
        */
        /*  basic title, version */
        if ( INCHI_BAS == iINChI )
        {
            inchi_ios_print( output_file, "\n" );   /*  empty line */
        }
        tot_len = sprintf(pStr, "%s<%s %s=\"%s\"",
            SP(ind), x_basic, x_ver, x_curr_ver);
        if ( INCHI_REC == iINChI || (INCHI_BAS == iINChI && bDisconnectedCoord) )
        {
            tot_len += sprintf(pStr+tot_len, " %s=\"%d\"", x_reconnected, iINChI );
        }
        if ( bAbcNumbers || bCtPredecessors )
        {
            const char *pNumber = "";
            const char *pDelim  = "";
            const char *pCtType = "";
            if ( bAbcNumbers && bCtPredecessors )
            {
                pNumber = x_type_short;
            }
            else
            {
                pNumber = bAbcNumbers? x_type_alpha : x_type_numer;
                pDelim  = (bAbcNumbers && bCtPredecessors)? "-":"";
                pCtType = bCtPredecessors? x_type_predec:"";
            }
            /*  type */
            tot_len += sprintf(pStr+tot_len, " %s=\"%s%s%s\"", x_type, pNumber, pDelim, pCtType);
        }
        sprintf(pStr+tot_len,">");
        inchi_ios_print( output_file, "%s\n", pStr );
        ind += inc;
    }
    else
    if ( INCHI_BAS == iINChI )
    {
        /* eliminate empty line in plain text output */
        if ( bNoStructLabels )
        {
            ;
/* -- removed empty line before InChI ---
#ifndef TARGET_API_LIB
            inchi_ios_print( output_file, "\n" );
#else
            ;
#endif
*/
        }
        else
        {
            if ( !(szSdfLabel && szSdfLabel[0]) && !(szSdfValue && szSdfValue[0]) )
            {
                tot_len = sprintf( pStr, "%sStructure: %d", pLF, num_input_struct );
                inchi_ios_print( output_file, "%s%s", pStr, pTAB );
            }
            else
            {
                tot_len = sprintf( pStr, "%sStructure: %d.%s%s%s%s",
                                         pLF,
                                        num_input_struct,
                                        SDF_LBL_VAL(szSdfLabel, szSdfValue) );
                if ( lSdfId )
                {
                    tot_len --;
                    tot_len += sprintf( pStr + tot_len, ":%ld", lSdfId );
                }
                inchi_ios_print( output_file, "%s%s", pStr, pTAB );
            }
        }
        /* inchi_ios_print( output_file, "%s%s=%s", pLF, (FLAG_SORT_PRINT_ReChI_PREFIX & *pSortPrintINChIFlags)? INCHI_REC_NAME : INCHI_NAME, pLF ); */
        inchi_ios_print( output_file, "%s%s=%s", pLF, INCHI_NAME, pLF );
    }


    /*****************************************************
     *
     * version  (10-29-2003)
     *
     ****************************************************/
    if ( INCHI_BAS == iINChI || !(bINChIOutputOptions & INCHI_OUT_EMBED_REC) /* || !bXml */)
    {
        /* xml: only if the first or not embedded; plain: always */
        szGetTag( IdentLbl, nTag,  bTag1 = IL_VERS, szTag1, &bAlways );
        tot_len = str_LineStart( szTag1, NULL, 0, pStr, ind );
        tot_len += sprintf(pStr + tot_len, "%s", x_curr_ver);

        /* 10-17-2008 Add 'standard' flag if necessary */
        if ( bINChIOutputOptions & INCHI_OUT_STDINCHI )
            tot_len += sprintf(pStr + tot_len, "S");

        /*if ( bXml ) {*/  /* avoid leading slash in plain output */
            if ( str_LineEnd( szTag1, tot_len, nStrLen, &bOverflow, pStr, bXml? 0 : -1, bPlainTextTags ) )
                goto exit_function;
        /*}*/
        inchi_ios_print( output_file, "%s%s", pStr, pLF );
    }


    /*****************************************************
     *
     * atoms, connection tables and tautomeric info
     *
     ****************************************************/
    /******************* constitution: dot-disconnected Hill formulas: <formula> */
    if ( num_components2[0] || num_components2[1] )
    {
        szGetTag( IdentLbl, nTag,  bTag1 = INCHI_REC == iINChI? IL_REC_ : IL_FML_, szTag1, &bAlways );
        tot_len = str_LineStart( szTag1, NULL, 0, pStr, ind );
        tot_len = str_HillFormula(pINChISort, pStr, nStrLen, tot_len,
                                  &bOverflow, bOutType, num_components, bUseMulipliers);

        if ( str_LineEnd( szTag1, tot_len, nStrLen, &bOverflow, pStr, bXml? 0 : -1, 1 ) )
            goto exit_function;
        inchi_ios_print( output_file, "%s%s", pStr, pLF );
    }
    /****************  semicolon/dot-disconnected connection tables */
    szGetTag( IdentLbl, nTag,  bTag1 = IL_CONN, szTag1, &bAlways );
    tot_len  = str_LineStart( szTag1, NULL, 0, pStr, ind );
    tot_len2 = str_Connections(pINChISort, pStr, nStrLen, tot_len,
                              &bOverflow, bOutType, ATOM_MODE, num_components, bUseMulipliers);
    /* current version does not output empty (";;;;") connectivity */
    if ( tot_len != tot_len2 /*|| !bXml*/ ) { /* 2004-06-30: never output empty connection table */
        tot_len = tot_len2;
        if ( str_LineEnd( szTag1, tot_len, nStrLen, &bOverflow, pStr, bXml? 0 : -2, bPlainTextTags ) )
            goto exit_function; /* pStr overfow */
        inchi_ios_print( output_file, "%s%s", pStr, pLF );
    }
    /************** hydrogen atoms; do not output empty */
    if ( INCHI_SEGM_FILL == INChI_SegmentAction( sDifSegs[nCurINChISegment][DIFS_h_H_ATOMS] ) )
    {
        szGetTag( IdentLbl, nTag,  bTag1 = IL_ALLH, szTag1, &bAlways );
        tot_len  = str_LineStart( szTag1, NULL, 0, pStr, ind );
        tot_len2 = str_H_atoms(pINChISort, pStr, nStrLen, tot_len,
                              &bOverflow, bOutType, ATOM_MODE, TAUT_MODE,
                              num_components, bUseMulipliers);
        if ( tot_len != tot_len2 /*|| !bXml*/ ) { /* 2004-06-21: never output empty */
            tot_len = tot_len2;
            if ( str_LineEnd( szTag1, tot_len, nStrLen, &bOverflow, pStr, bXml? 0 : -2, 1 ) )
                goto exit_function;
            inchi_ios_print( output_file, "%s%s", pStr, pLF );
        }
    }


    bFhTag = 0;


repeat_INChI_output:

    /*****************************************************
     * charge
     */

    nSegmAction = INChI_SegmentAction( sDifSegs[nCurINChISegment][DIFS_q_CHARGE] );
    if ( nSegmAction )
    {
        szGetTag( IdentLbl, nTag,  bTag1 = IL_CHRG | bFhTag, szTag1, &bAlways );
        tot_len = str_LineStart( szTag1, NULL, 0, pStr, ind );
        if ( INCHI_SEGM_FILL == nSegmAction )
        {
            tot_len = str_Charge2(pINChISort, pINChISort2, pStr, nStrLen, tot_len,
                                  &bOverflow, bOutType, num_components,
                                  bSecondNonTautPass, bOmitRepetitions, bUseMulipliers);
            bNonTautNonIsoIdentifierNotEmpty += bSecondNonTautPass;
        }
        if ( str_LineEnd( szTag1, tot_len, nStrLen, &bOverflow, pStr, bXml? 0 : -nSegmAction, bPlainTextTags ) )
            goto exit_function;
        inchi_ios_print( output_file, "%s%s", pStr, pLF );
    }
    else
    {
        if ( !bXml )
        {
            if ( bPlainTextTags == 1 ) inchi_ios_print( output_file, "/" );
        }
    }

    /*****************************************************
     * removed protons
     */
    if ( iCurTautMode == TAUT_YES && !bSecondNonTautPass )
    {

        nSegmAction = INChI_SegmentAction( sDifSegs[nCurINChISegment][DIFS_p_PROTONS] );
        if ( nSegmAction )
        {
            szGetTag( IdentLbl, nTag,  bTag1 = IL_PROT | bFhTag, szTag1, &bAlways );
            tot_len = str_LineStart( szTag1, NULL, 0, pStr, ind );
            tot_len += sprintf( pStr + tot_len, "%+d", nNumRemovedProtons );
            if ( str_LineEnd( szTag1, tot_len, nStrLen, &bOverflow, pStr, bXml? 0 : -nSegmAction, bPlainTextTags ) )
                goto exit_function;
            inchi_ios_print( output_file, "%s%s", pStr, pLF );
        }
        else
        {
            if ( !bXml )
            {
                if ( bPlainTextTags == 1 ) inchi_ios_print( output_file, "/" );
            }
        }

    }


    /**************************************************
     *
     *    non-isotopic stereo
     */

    {
        int i;
        i = INChI_SegmentAction( sDifSegs[nCurINChISegment][DIFS_t_SATOMS] );
        i = i;
    }

    if ( INChI_SegmentAction( sDifSegs[nCurINChISegment][DIFS_b_SBONDS] ) ||
         INChI_SegmentAction( sDifSegs[nCurINChISegment][DIFS_t_SATOMS] ) ||
         INChI_SegmentAction( sDifSegs[nCurINChISegment][DIFS_m_SP3INV] ) ||
         INChI_SegmentAction( sDifSegs[nCurINChISegment][DIFS_s_STYPE] ) )
    {
        /*  stereo */
        szGetTag( IdentLbl, nTag,  bTag1 = IL_STER | bFhTag, szTag1, &bAlways );
        if ( bXml )
        {
            str_LineStart( szTag1, NULL, 0, pStr, ind );
            inchi_ios_print( output_file, "%s\n", pStr );
            ind += inc;
        }

        /*  sp2 */
        /*if ( bStereoSp2[iCurTautMode]  )*/
        if ( (nSegmAction = INChI_SegmentAction( sDifSegs[nCurINChISegment][DIFS_b_SBONDS] )) )
        {
            szGetTag( IdentLbl, nTag,  bTag2 = bTag1 | IL_DBND, szTag2, &bAlways );
            tot_len = str_LineStart( szTag2, NULL, 0, pStr, ind );
            if ( INCHI_SEGM_FILL == nSegmAction )
            {
                tot_len = str_Sp2(pINChISort, pINChISort2, pStr, nStrLen, tot_len,
                                   &bOverflow, bOutType, TAUT_MODE, num_components,
                                   bSecondNonTautPass, bOmitRepetitions, bUseMulipliers);
                bNonTautNonIsoIdentifierNotEmpty += bSecondNonTautPass;
            }
            if ( str_LineEnd( szTag2, tot_len, nStrLen, &bOverflow, pStr, bXml? 0 : -nSegmAction, bPlainTextTags ) )
                goto exit_function;
            inchi_ios_print( output_file, "%s%s", pStr, pLF );
        }
        else
        {
            if ( !bXml )
            {
                if ( bPlainTextTags == 1 ) inchi_ios_print( output_file, "/" ); /* sp2 */
            }
        }

        /*  sp3 */
        /*if ( bStereoSp3[iCurTautMode]  )*/
        if ( (nSegmAction = INChI_SegmentAction( sDifSegs[nCurINChISegment][DIFS_t_SATOMS] )) )
        {
            bRelRac     = bRelativeStereo[iCurTautMode] || bRacemicStereo[iCurTautMode];
            szGetTag( IdentLbl, nTag,  bTag2 = bTag1 | IL_SP3S, szTag2, &bAlways );
            tot_len = str_LineStart( szTag2, NULL, 0, pStr, ind );
            if ( INCHI_SEGM_FILL == nSegmAction ) {
                tot_len = str_Sp3(pINChISort, pINChISort2, pStr, nStrLen, tot_len,
                                  &bOverflow, bOutType, TAUT_MODE, num_components, bRelRac,
                                  bSecondNonTautPass, bOmitRepetitions, bUseMulipliers);
                bNonTautNonIsoIdentifierNotEmpty += bSecondNonTautPass;
            }

            if (str_LineEnd( szTag2, tot_len, nStrLen, &bOverflow, pStr, bXml? 0 : -nSegmAction, bPlainTextTags ))
                goto exit_function;
            inchi_ios_print( output_file, "%s%s", pStr, pLF );
        }
        else
        {
            if ( !bXml )
            {
                if ( bPlainTextTags == 1 ) inchi_ios_print( output_file, "/" ); /* sp3 */
            }
        }

        /* bStereoAbsInverted[iCurTautMode]  */

        /* if ( bStereoAbs[iCurTautMode]  ) */
        if ( (nSegmAction = INChI_SegmentAction( sDifSegs[nCurINChISegment][DIFS_m_SP3INV] )) )
        {
            szGetTag( IdentLbl, nTag,  bTag2 = bTag1 | IL_INVS, szTag2, &bAlways );
            tot_len = str_LineStart( szTag2, NULL, 0, pStr, ind );
            if ( INCHI_SEGM_FILL == nSegmAction ) {
                tot_len = str_StereoAbsInv(pINChISort, pStr, nStrLen, tot_len,
                                           &bOverflow, bOutType, num_components);
                bNonTautNonIsoIdentifierNotEmpty += bSecondNonTautPass;
            }

            if (str_LineEnd( szTag2, tot_len, nStrLen, &bOverflow, pStr, bXml? 0 : -nSegmAction, bPlainTextTags ))
                goto exit_function;
            inchi_ios_print( output_file, "%s%s", pStr, pLF );
        }
        else
        {
            if ( !bXml )
            {
                if ( bPlainTextTags == 1 ) inchi_ios_print( output_file, "/" ); /* stereo-abs-inv */
            }
        }

        /* stereo type */
        /*if ( bRacemicStereo[iCurTautMode] || bRelativeStereo[iCurTautMode] || bStereoAbs[iCurTautMode] )*/
        if ( (nSegmAction = INChI_SegmentAction( sDifSegs[nCurINChISegment][DIFS_s_STYPE] )) )
        {
            const char *p_stereo = bRelativeStereo[iCurTautMode]? x_rel :
                                   bRacemicStereo[iCurTautMode] ? x_rac : x_abs;
            szGetTag( IdentLbl, nTag,  bTag2 = bTag1 | IL_TYPS, szTag2, &bAlways );
            tot_len = str_LineStart( szTag2, NULL, 0, pStr, ind );
            if ( INCHI_SEGM_FILL == nSegmAction ) {
                tot_len += MakeDelim( p_stereo, pStr + tot_len, nStrLen-tot_len, &bOverflow);
                bNonTautNonIsoIdentifierNotEmpty += bSecondNonTautPass;
            }
            if (str_LineEnd( szTag2, tot_len, nStrLen, &bOverflow, pStr, bXml? 0 : -nSegmAction, bPlainTextTags ))
                goto exit_function;
            inchi_ios_print( output_file, "%s%s", pStr, pLF );
        }
        if ( !bXml )
        {
            if ( bPlainTextTags == 1 ) inchi_ios_print( output_file, "/" );  /* no abs, inv or racemic stereo */
        }

        if ( bXml )
        {
            /* close stereo */
            ind -= inc;
            if ( str_LineEnd( szTag1, 0, nStrLen, &bOverflow, pStr, ind, bPlainTextTags ) )
                goto exit_function;
            inchi_ios_print( output_file, "%s", pStr );
        }
    }
    else
    {
        if ( !bXml )
        {
            if ( bPlainTextTags == 1 ) inchi_ios_print( output_file, "////" ); /* sp3, sp2, abs-inv, stereo.type */
        }
    }


    /****************************************************
     *
     *  Isotopic canonical results
     *
     ****************************************************/
    nCurINChISegment ++; /* switch from M to MI or from F to FI */

    /*if ( bIsotopic || !bSecondNonTautPass && bHasIsoH )*/
    if ( INChI_SegmentAction( sDifSegs[nCurINChISegment][DIFS_i_IATOMS] ) )
    {
        /*  isotopic #1:  composition -- atoms -- do not output in xml if empty */
        szGetTag( IdentLbl, nTag,  bTag1 = IL_ISOT | bFhTag, szTag1, &bAlways );
        if ( bXml )
        {
            str_LineStart( szTag1, NULL, 0, pStr, ind );
            inchi_ios_print( output_file, "%s\n", pStr );
            ind += inc;
        }
        /* isotopic atoms without mobile H.
         * Fixed 2004-06-15: always output if not bXml. Note:
         * Previous condition if( bHasIsotopicAtoms[iCurTautMode] || bIsotopic && !bXml)
         * did not optput /i in case of only mobile isotopic H
         */
        if ( (nSegmAction = INChI_SegmentAction( sDifSegs[nCurINChISegment][DIFS_i_IATOMS] )) )
        {
            szGetTag( IdentLbl, nTag,  bTag2 = bTag1 | IL_ATMS, szTag2, &bAlways );
            tot_len = str_LineStart( szTag2, NULL, 0, pStr, ind );
            /*if ( bHasIsotopicAtoms[iCurTautMode] )*/
            if ( INCHI_SEGM_FILL == nSegmAction )
            {
                tot_len2 = str_IsoAtoms(pINChISort, pINChISort2, pStr, nStrLen, tot_len,
                                       &bOverflow, bOutType, TAUT_MODE, num_components, bAbcNumbers,
                                       bSecondNonTautPass, bOmitRepetitions, bUseMulipliers);
                bNonTautIsoIdentifierNotEmpty += bSecondNonTautPass;
            }
            else
            {
                tot_len2 = tot_len;
            }

            tot_len = tot_len2;
            if ( str_LineEnd( szTag2, tot_len, nStrLen, &bOverflow, pStr, bXml? 0 : -nSegmAction, bPlainTextTags ) )
                goto exit_function;
            inchi_ios_print( output_file, "%s%s", pStr, pLF );

        }

        /*  isotopic #1a:  composition -- exchangeable isotopic H (mobile H only) */
        /*if ( !bSecondNonTautPass && bHasIsoH )*/
        if ( (nSegmAction = INChI_SegmentAction( sDifSegs[nCurINChISegment][DIFS_h_H_ATOMS] )) )
        {
            szGetTag( IdentLbl, nTag,  bTag2 = bTag1 | IL_XCGA, szTag2, &bAlways );
            tot_len = str_LineStart( szTag2, NULL, 0, pStr, ind );
            tot_len += MakeIsoHString( num_iso_H, pStr + tot_len, nStrLen-tot_len, TAUT_MODE, &bOverflow);
            bNonTautIsoIdentifierNotEmpty += bSecondNonTautPass;
            if ( str_LineEnd( szTag2, tot_len, nStrLen, &bOverflow, pStr, bXml? 0 : -nSegmAction, bPlainTextTags ) )
                goto exit_function;
            inchi_ios_print( output_file, "%s%s", pStr, pLF );
        }

        /***************************************************
         *
         *       Isotopic stereo
         *
         ***************************************************/

        /*if ( bIsotopicStereo[iCurTautMode] )*/
        if ( INChI_SegmentAction( sDifSegs[nCurINChISegment][DIFS_b_SBONDS] ) ||
             INChI_SegmentAction( sDifSegs[nCurINChISegment][DIFS_t_SATOMS] ) ||
             INChI_SegmentAction( sDifSegs[nCurINChISegment][DIFS_m_SP3INV] ) ||
             INChI_SegmentAction( sDifSegs[nCurINChISegment][DIFS_s_STYPE] ) )
        {
            /*  stereo */
            szGetTag( IdentLbl, nTag,  bTag2 = bTag1 | IL_STER, szTag2, &bAlways );
            if ( bXml )
            {
                str_LineStart( szTag2, NULL, 0, pStr, ind );
                inchi_ios_print( output_file, "%s\n", pStr );
                ind += inc;
            }

            /************************
              isotopic #2:  sp2
             ************************/
            /*if ( bIsotopicStereoSp2[iCurTautMode]  )*/
            if ( (nSegmAction = INChI_SegmentAction( sDifSegs[nCurINChISegment][DIFS_b_SBONDS] )) )
            {
                szGetTag( IdentLbl, nTag,  bTag3 = bTag2 | IL_DBND, szTag3, &bAlways );
                tot_len = str_LineStart( szTag3, NULL, 0, pStr, ind );
                if ( INCHI_SEGM_FILL == nSegmAction ) {
                    tot_len = str_IsoSp2(pINChISort, pINChISort2, pStr, nStrLen, tot_len,
                                         &bOverflow, bOutType, TAUT_MODE, num_components,
                                         bSecondNonTautPass, bOmitRepetitions, bUseMulipliers);
                    bNonTautIsoIdentifierNotEmpty += bSecondNonTautPass;
                }
                if ( str_LineEnd( szTag3, tot_len, nStrLen, &bOverflow, pStr, bXml? 0 : -nSegmAction, bPlainTextTags ) )
                    goto exit_function;
                inchi_ios_print( output_file, "%s%s", pStr, pLF );
            }
            else
            {
                if ( !bXml )
                {
                    if ( bPlainTextTags == 1 ) inchi_ios_print( output_file, "/" ); /* iso sp2 */
                }
            }

            /************************
              isotopic #3:  sp3
             ************************/
            /*if ( bIsotopicStereoSp3[iCurTautMode]  )*/
            if ( (nSegmAction = INChI_SegmentAction( sDifSegs[nCurINChISegment][DIFS_t_SATOMS] )) )
            {
                bRelRac = bIsotopicRelativeStereo[iCurTautMode] || bIsotopicRacemicStereo[iCurTautMode];
                szGetTag( IdentLbl, nTag,  bTag3 = bTag2 | IL_SP3S, szTag3, &bAlways );
                tot_len = str_LineStart( szTag3, NULL, 0, pStr, ind );
                if ( INCHI_SEGM_FILL == nSegmAction ) {
                    tot_len = str_IsoSp3(pINChISort, pINChISort2, pStr, nStrLen, tot_len,
                                         &bOverflow, bOutType, TAUT_MODE, num_components, bRelRac,
                                         bSecondNonTautPass, bOmitRepetitions, bUseMulipliers);
                    bNonTautIsoIdentifierNotEmpty += bSecondNonTautPass;
                }
                if ( str_LineEnd( szTag3, tot_len, nStrLen, &bOverflow, pStr, bXml? 0 : -nSegmAction, bPlainTextTags ) )
                    goto exit_function;
                inchi_ios_print( output_file, "%s%s", pStr, pLF );
            } else
            {
                if ( !bXml )
                {
                    if ( bPlainTextTags == 1 )
                        inchi_ios_print( output_file, "/" ); /* iso-sp3 */
                }
            }

            /* isotopic #4: abs inverted */
            if ( (nSegmAction = INChI_SegmentAction( sDifSegs[nCurINChISegment][DIFS_m_SP3INV] )) )
            {
                szGetTag( IdentLbl, nTag,  bTag3 = bTag2 | IL_INVS, szTag3, &bAlways );
                tot_len = str_LineStart( szTag3, NULL, 0, pStr, ind );
                if ( INCHI_SEGM_FILL == nSegmAction ) {
                    tot_len = str_IsoStereoAbsInv(pINChISort, pStr, nStrLen, tot_len,
                                                  &bOverflow, bOutType, num_components);
                    bNonTautIsoIdentifierNotEmpty += bSecondNonTautPass;
                }
                if ( str_LineEnd( szTag3, tot_len, nStrLen, &bOverflow, pStr, bXml? 0 : -nSegmAction, bPlainTextTags ) )
                    goto exit_function;
                inchi_ios_print( output_file, "%s%s", pStr, pLF );
            }
            else
            {
                if ( !bXml )
                {
                    if ( bPlainTextTags == 1 )
                        inchi_ios_print( output_file, "/" );
                }
            }

            /* isotopic #5: stereo type. Do not output if it has already been output in non-iso */
            if ( (nSegmAction = INChI_SegmentAction( sDifSegs[nCurINChISegment][DIFS_s_STYPE] )) )
            {
                const char *p_stereo = bIsotopicRelativeStereo[iCurTautMode]? x_rel :
                                       bIsotopicRacemicStereo[iCurTautMode] ? x_rac : x_abs;
                szGetTag( IdentLbl, nTag,  bTag3 = bTag2 | IL_TYPS, szTag3, &bAlways );
                tot_len = str_LineStart( szTag3, NULL, 0, pStr, ind );
                if ( INCHI_SEGM_FILL == nSegmAction ) {
                    tot_len += MakeDelim( p_stereo, pStr + tot_len, nStrLen-tot_len, &bOverflow);
                    bNonTautIsoIdentifierNotEmpty += bSecondNonTautPass;
                }
                if ( str_LineEnd( szTag3, tot_len, nStrLen, &bOverflow, pStr, bXml? 0 : -nSegmAction, bPlainTextTags ) )
                    goto exit_function;
                inchi_ios_print( output_file, "%s%s", pStr, pLF );
            }
            if ( !bXml )
            {
                if ( bPlainTextTags == 1 )
                    inchi_ios_print( output_file, "/" );  /* no abs, inv or racemic stereo */
            }
            if ( bXml )
            {
                /************************
                  close isotopic stereo
                 ************************/
                ind -= inc;
                if ( str_LineEnd( szTag2, 0, nStrLen, &bOverflow, pStr, ind, bPlainTextTags ) )
                    goto exit_function;
                inchi_ios_print( output_file, "%s", pStr );
            }
        }
        else
        {
            if ( !bXml )
            {
                /* no isotopic stereo */
                if ( bPlainTextTags == 1 )
                    inchi_ios_print( output_file, "////" ); /* sp3, sp2, abs-inv, stereo.type */
            }
        }

        /*  close isotopic */
        if ( bXml )
        {
            ind -= inc;
            if ( str_LineEnd( szTag1, 0, nStrLen, &bOverflow, pStr, ind, bPlainTextTags ) )
                goto exit_function;
            inchi_ios_print( output_file, "%s", pStr );
        }
    }
    else
    {
        if ( !bXml )
        {
            if ( bPlainTextTags == 1 )
                inchi_ios_print( output_file, "///" ); /* isotopic composition, sp2, sp3 */
            if ( bPlainTextTags == 1 )
                inchi_ios_print( output_file, "//" );   /* inv or racemic stereo */
        }
    }

#if ( CANON_FIXH_TRANS == 1 )
    if ( bOutType == OUT_NONTAUT && bOutputType == OUT_TN && bSecondNonTautPass &&
         INCHI_SEGM_FILL == INChI_SegmentAction( sDifSegs[DIFL_F][DIFS_o_TRANSP] ))
    {
        /* find and print non-tautomeric components transposition, if non-trivial */
        AT_NUMB *nTrans_n, *nTrans_s;

        if ( 0 < bin_AuxTautTrans(pINChISort, pINChISort2, &nTrans_n, &nTrans_s, bOutType,  num_components) )
        {
            /* a non-trivial transposition does exist; output start tag */
            szGetTag( IdentLbl, nTag,  bTag1 = IL_TRNS | bFhTag, szTag1, &bAlways );
            tot_len = str_LineStart( szTag1, NULL, 0, pStr, ind );
            /* print the transposition, cycle after cycle */
            tot_len = str_AuxTautTrans(nTrans_n, nTrans_s, pStr, nStrLen, tot_len,
                                       &bOverflow, TAUT_MODE, num_components);
            bNonTautIsoIdentifierNotEmpty += bSecondNonTautPass;
            if ( str_LineEnd( szTag1, tot_len, nStrLen, &bOverflow, pStr, bXml? 0 : -1, bPlainTextTags ) )
                goto exit_function;
            inchi_ios_print( output_file, "%s%s", pStr, pLF );
             /* detected transposition */
            *pSortPrintINChIFlags |= (INCHI_BAS == iINChI)? FLAG_SORT_PRINT_TRANSPOS_BAS :
                                                            FLAG_SORT_PRINT_TRANSPOS_REC;
        }
        else
        {
            if ( !bXml )
            {
                if ( bPlainTextTags == 1 )
                    inchi_ios_print( output_file, "/" );
            }
        }
    }
#endif


    /**************************************************************
      At this point the INChI part of the output has been done.
      If this INChI is tautomeric and non-tautomeric results exist
      then we need to output non-tautomeric data:
         fixed H,
         stereo,
         isotopic
         isotopic stereo
    ***************************************************************/
    if ( bOutType == OUT_TN && !bSecondNonTautPass &&
         bNonTautIsIdenticalToTaut && bTautomeric && bNonTautomeric )
    {
            /* Fixed-H layer is empty in the Identifier */
            *pSortPrintINChIFlags |= (INCHI_BAS == iINChI)? FLAG_SORT_PRINT_NO_NFIX_H_BAS :
                                                            FLAG_SORT_PRINT_NO_NFIX_H_REC;
            *pSortPrintINChIFlags |= (INCHI_BAS == iINChI)? FLAG_SORT_PRINT_NO_IFIX_H_BAS :
                                                            FLAG_SORT_PRINT_NO_IFIX_H_REC;
    }

    if ( bOutType == OUT_TN && !bNonTautIsIdenticalToTaut /* added 2004-10-04 Fix16 */
#ifdef OLD_ITEM_DISCOVERY
                            && bTautomeric && bNonTautomeric
#endif
                            && INChI_SegmentAction( sDifSegs[DIFL_F][DIFS_f_FORMULA] )
                       /* special case: removed isolated H(+): */
                       /* || iCurTautMode == TAUT_YES && num_comp[TAUT_YES] < num_comp[TAUT_NON] &&
                             0 < num_comp[TAUT_NON]*/
       )
    {
        /* add the second (non-tautomeric) output */
        bOutType     = OUT_NONTAUT;    /* pick up only non-tautomeric representation of tautomeric */
        iCurTautMode = TAUT_NON;
        pINChISort    = pINChISortTautAndNonTaut[TAUT_NON];
        bSecondNonTautPass = 1;
        nCurINChISegment   = DIFL_F;
        num_components = num_comp[iCurTautMode]; /* number of components could change due to removal of isolated H(+) from tautomeric */
        bFhTag = IL_FIXH;
        szGetTag( IdentLbl, nTag,  bTag1 = bFhTag, szTag1, &bAlways );
        if ( bXml )
        {
            /* open non-tautomeric */
            str_LineStart( szTag1, NULL, 0, pStr, ind );
            inchi_ios_print( output_file, "%s\n", pStr );
            ind += inc;
        }
        /***** constitution non-taut: dot-disconnected Hill formulas: <formula> -- only if different */
        szGetTag( IdentLbl, nTag,  bTag1 = IL_FMLF | bFhTag, szTag1, &bAlways );
        tot_len = str_LineStart( szTag1, NULL, 0, pStr, ind );
        nSegmAction = INChI_SegmentAction( sDifSegs[nCurINChISegment][DIFS_f_FORMULA] );
        if ( INCHI_SEGM_FILL == nSegmAction )
        {
            tot_len2 = str_HillFormula2(pINChISort, pINChISort2, pStr, nStrLen, tot_len,
                                      &bOverflow, bOutType, num_components, bUseMulipliers);
            bNonTautNonIsoIdentifierNotEmpty += bSecondNonTautPass;
        }
        else
        {
            tot_len2 = tot_len;
        }
        tot_len = tot_len2;
        if ( str_LineEnd( szTag1, tot_len, nStrLen, &bOverflow, pStr, bXml? 0 : -nSegmAction, bPlainTextTags ) )
            goto exit_function;
        inchi_ios_print( output_file, "%s%s", pStr, pLF );

        nSegmAction = INChI_SegmentAction( sDifSegs[nCurINChISegment][DIFS_h_H_ATOMS] );
        if ( INCHI_SEGM_FILL == nSegmAction )
        {
            szGetTag( IdentLbl, nTag,  bTag1 = IL_HFIX | bFhTag, szTag1, &bAlways );
            tot_len = str_LineStart( szTag1, NULL, 0, pStr, ind ); /* open H-fixed */
            /* output the second non-tautomeric item: fixed H -- do not output in xml if empty */
            tot_len2 = str_FixedH_atoms(pINChISort, pStr, nStrLen, tot_len,
                                       &bOverflow, bOutType, ATOM_MODE, num_components, bUseMulipliers);
            tot_len = tot_len2;
            if ( str_LineEnd( szTag1, tot_len, nStrLen, &bOverflow, pStr, bXml? 0 : -nSegmAction, bPlainTextTags ) )
                goto exit_function;
            inchi_ios_print( output_file, "%s%s", pStr, pLF );
            bNonTautNonIsoIdentifierNotEmpty += bSecondNonTautPass;
        }
        goto repeat_INChI_output;
    }
    else
    {
        if ( bOutType == OUT_NONTAUT && bOutputType == OUT_TN && bSecondNonTautPass /* && bTautomeric && bNonTautomeric*/ )
        {
            /* the second (non-taut) output has been done; restore variables */
            bOutType           = OUT_TN;
            iCurTautMode       = TAUT_YES;
            pINChISort         = pINChISortTautAndNonTaut[TAUT_YES];
            bSecondNonTautPass = 0;
            num_components     = num_comp[iCurTautMode];
            if ( !bNonTautNonIsoIdentifierNotEmpty )
            {
                /* Fixed-H layer is empty in the Identifier */
                *pSortPrintINChIFlags |= (INCHI_BAS == iINChI)? FLAG_SORT_PRINT_NO_NFIX_H_BAS :
                                                            FLAG_SORT_PRINT_NO_NFIX_H_REC;
            }
            if ( !bNonTautIsoIdentifierNotEmpty )
            {
                /* Fixed-H layer is empty in the Identifier */
                *pSortPrintINChIFlags |= (INCHI_BAS == iINChI)? FLAG_SORT_PRINT_NO_IFIX_H_BAS :
                                                                FLAG_SORT_PRINT_NO_IFIX_H_REC;
            }
            if ( bXml )
            {
                /*  close non-tautomeric */
                ind -= inc;
                szGetTag( IdentLbl, nTag,  bTag1 = bFhTag, szTag1, &bAlways );
                if ( str_LineEnd( szTag1, 0, nStrLen, &bOverflow, pStr, ind, bPlainTextTags ) )
                    goto exit_function;
                inchi_ios_print( output_file, "%s", pStr );
            }
            bFhTag             = 0;
        }
    }


    /************************************************
     * output INChI     of the reconnected structure *
     ************************************************/
    bEmbeddedOutputCalled = 0;
    if ( bDisconnectedCoord && INCHI_BAS == iINChI &&
         (bINChIOutputOptions & INCHI_OUT_EMBED_REC) && num_components2[INCHI_REC] )
    {
        int nRet;
        bEmbeddedOutputCalled = 1;

        if ( !bXml )
        {
             /* output blank line before /R: in case of bPlainTextCommnts=1 */
            inchi_ios_print( output_file, "%s", pLF );
        }
        /* end of disconnected INChI output */

        nRet = OutputINChI1( pStr, nStrLen,
                             pINChISortTautAndNonTaut2,
                             INCHI_REC, NULL,
                             0 /*bDisconnectedCoord*/, bOutputType,
                             bINChIOutputOptions | INCHI_OUT_NO_AUX_INFO,
                             bXml, bAbcNumbers, bCtPredecessors, bNoStructLabels,
                             num_components2, num_non_taut2, num_taut2,
                             output_file, log_file,
                             num_input_struct,
                             szSdfLabel, szSdfValue, lSdfId,
                             pSortPrintINChIFlags,
                             save_opt_bits);
        if ( !nRet )
            goto exit_function; /* error */
    }

    if ( bXml )
    {
        /*  close INChI identifier (basic) */
        ind -= inc;
        if ( str_LineEnd( x_basic, 0, nStrLen, &bOverflow, pStr, ind, bPlainTextTags ) )
            goto exit_function;
        inchi_ios_print( output_file, "%s", pStr );
    }
    else
    {
        /* save InChI creation options if requested ...*/
        if ( !bEmbeddedOutputCalled)
        {
            if ( bINChIOutputOptions & INCHI_OUT_SAVEOPT )
            {
                /* ... and not std-InChI output */
                if ( 0 == (bINChIOutputOptions & INCHI_OUT_STDINCHI) )
                {
                    char let1, let2;
                    GetSaveOptLetters(save_opt_bits, &let1, &let2);
                    inchi_ios_print( output_file, "\\%c%c", let1, let2 );
                }
            }
        }

        if ( !bEmbeddedOutputCalled && !bPlainTextCommnts )
        { /* plain text comment earlier ended with LF */
            inchi_ios_print( output_file, "%s%s",
                                (!num_components2[0] && !num_components2[1])? "//":"", /* empty InChI=// */
                                (bINChIOutputOptions & INCHI_OUT_NO_AUX_INFO)? "\n" : pTAB );
        /* end of INChI= output */
        }


    }

output_aux_info:

    bFhTag = 0;

    if( !(bINChIOutputOptions & INCHI_OUT_NO_AUX_INFO) )
    {
        /* output aux info */

        /*************************************************************
         *
         *   Aux info non-isotopic
         *
         *************************************************************/

        num_components = num_comp[iCurTautMode];
        if ( bXml )
        {
            /*  aux. info header */
            /*  empty line if INChI output has been printed */
            if ( !(bINChIOutputOptions & INCHI_OUT_ONLY_AUX_INFO) )
            {
                inchi_ios_print( output_file, "\n" );
            }
            /*  basic.aux-info title, version */
            tot_len = sprintf(pStr, "%s<%s %s=\"%s\"",
                SP(ind), x_aux_basic, x_ver, x_curr_ver );
            if ( INCHI_REC == iINChI || (INCHI_BAS == iINChI && bDisconnectedCoord) )
            {
                tot_len += sprintf(pStr+tot_len, " %s=\"%d\"", x_reconnected, iINChI );
            }
            if ( bAbcNumbers )
            {
                /*  type */
                const char *pNumber = x_type_short;
                tot_len += sprintf(pStr+tot_len, " %s=\"%s\"", x_type, pNumber);
            }

            sprintf(pStr+tot_len,">");
            inchi_ios_print( output_file, "%s\n", pStr );
            ind += inc;
            if ( !(bINChIOutputOptions & INCHI_OUT_ONLY_AUX_INFO) )
            {
                /*  comment */
                tot_len = sprintf( pStr, "%s<%s>", SP(ind), x_aux_comm );
                inchi_ios_print( output_file, "%s\n", pStr );
            }
        }
        else
        {
            if ( INCHI_BAS == iINChI )
            {
                tot_len = sprintf( pStr, "AuxInfo=" ); /* in wINChI window, separate INChI: from AuxInfo: with blank line */
                inchi_ios_print( output_file, "%s%s%s",
                                          /* blank line before AuxInfo in winchi window unless it is an annotation */
                                          (bINChIOutputOptions & INCHI_OUT_WINCHI_WINDOW) ? "\n":"",
                                          pStr,
                                          pLF);
                szGetTag( AuxLbl, nTag,  bTag1 = AL_VERS, szTag1, &bAlways );
                tot_len = str_LineStart( szTag1, NULL, 0, pStr, ind );
                tot_len += sprintf(pStr + tot_len, "%s", x_curr_ver);
                /* avoid leading slash in plain output */
                if ( str_LineEnd( szTag1, tot_len, nStrLen, &bOverflow, pStr, bXml? 0 : -1, bPlainTextTags ) )
                    goto exit_function;
                inchi_ios_print( output_file, "%s%s", pStr, pLF );
            }
            else
            {
                if ( INCHI_REC == iINChI )
                {
                    szGetTag( AuxLbl, nTag,  bTag1 = AL_REC_, szTag1, &bAlways );
                    inchi_ios_print( output_file, "%s%s", szTag1, pLF );
                }
            }

        }
        /* normalization type */
        if ( num_components2[0] || num_components2[1] )
        {
            szGetTag( AuxLbl, nTag,  bTag1 = AL_NORM, szTag1, &bAlways );
            tot_len = str_LineStart( szTag1, NULL, 0, pStr, ind );
            tot_len += sprintf( pStr + tot_len, "%d", (bTautomeric && bTautomericOutputAllowed)? bTautomeric : 0);
            if ( str_LineEnd( szTag1, tot_len, nStrLen, &bOverflow, pStr, bXml? 0 : -1, bPlainTextTags ) )
                goto exit_function;
            inchi_ios_print( output_file, "%s%s", pStr, pLF );
        }



repeat_INChI_Aux_output:

        /**************************************************************
         *   Original atom numbers in order of canonical numbers
         **************************************************************/
        if ( num_components2[0] || num_components2[1] )
        {
            szGetTag( AuxLbl, nTag,  bTag1 = (bSecondNonTautPass? AL_FIXN : AL_ANBR) | bFhTag, szTag1, &bAlways );
            tot_len = str_LineStart( szTag1, NULL, 0, pStr, ind );
            /* original numbering output */
            tot_len = str_AuxNumb(pINChISort, pINChISort2, pStr, nStrLen, tot_len,
                                  &bOverflow, bOutType, TAUT_MODE, num_components,
                                  bSecondNonTautPass, bOmitRepetitions);

            if ( str_LineEnd( szTag1, tot_len, nStrLen, &bOverflow, pStr, bXml? 0 : -1, bPlainTextTags ) )
                goto exit_function;
            inchi_ios_print( output_file, "%s%s", pStr, pLF );
        }
        /**********************************************
         *   Symmetry numbers (constit. equivalence)
         **********************************************/
        if ( bAtomEqu[iCurTautMode] )
        {
            /*  aux equ atoms */
            /* 1. Compare to tautomeric equivalence (in case of second, non-taut, pass only) */
            /* 2. Compare to the previous component if (1) failed to find equivalence */
            szGetTag( AuxLbl, nTag,  bTag1 = AL_AEQU | bFhTag, szTag1, &bAlways );
            tot_len = str_LineStart( szTag1, NULL, 0, pStr, ind );
            tot_len = str_AuxEqu(pINChISort, pINChISort2, pStr, nStrLen, tot_len,
                                 &bOverflow, bOutType, TAUT_MODE, num_components,
                                 bSecondNonTautPass, bOmitRepetitions, bUseMulipliers);

            if ( str_LineEnd( szTag1, tot_len, nStrLen, &bOverflow, pStr, bXml? 0 : -1, bPlainTextTags ) )
                goto exit_function;
            inchi_ios_print( output_file, "%s%s", pStr, pLF );
        }
        else
        {
            if ( !bXml )
            {
                if ( bPlainTextTags == 1 )
                    inchi_ios_print( output_file, "/" );
            }
        }
        /*****************************************************
         *    Tautomeric groups equivalence
         *****************************************************/
        if ( bTautomericOutputAllowed && bTautomeric && bTautEqu[iCurTautMode] && !bSecondNonTautPass )
        {
            /*****************************************************
             *    Tautomeric groups constitutional equivalence
             */
            /*  aux tgroup equ */
            szGetTag( AuxLbl, nTag,  bTag1 = AL_GEQU | bFhTag, szTag1, &bAlways );
            tot_len = str_LineStart( szTag1, NULL, 0, pStr, ind );
            tot_len = str_AuxTgroupEqu(pINChISort, pStr, nStrLen, tot_len,
                                       &bOverflow, bOutType, TAUT_MODE,
                                       num_components, bUseMulipliers);
            if ( str_LineEnd( szTag1, tot_len, nStrLen, &bOverflow, pStr, bXml? 0 : -1, bPlainTextTags ) )
                goto exit_function;
            inchi_ios_print( output_file, "%s", pStr );
        }
        else
        {
            if ( !bXml && bTautomericOutputAllowed && bTautomeric )
            {
                if ( bPlainTextTags == 1 )
                    inchi_ios_print( output_file, "/" );
            }
        }

        /****************************************************
         * Inverted stereo -- sp3 only + canonical numbering
         ****************************************************/
        if ( bInvStereo[iCurTautMode] )
        {
            szGetTag( AuxLbl, nTag,  bTag1 = AL_STER | bFhTag, szTag1, &bAlways );
            if ( bXml )
            {
                /***************************
                     inv stereo start  tag
                ****************************/
                str_LineStart( szTag1, NULL, 0, pStr, ind );
                inchi_ios_print( output_file, "%s\n", pStr );
                ind += inc;
            }
            /****************************
                 inverted sp3 start tag
            *****************************/
            szGetTag( AuxLbl, nTag,  bTag2 = bTag1 | AL_SP3I, szTag2, &bAlways );
            tot_len = str_LineStart( szTag2, NULL, 0, pStr, ind );
            tot_len = str_AuxInvSp3(pINChISort, pINChISort2, pStr, nStrLen, tot_len,
                                    &bOverflow, bOutType, TAUT_MODE, num_components,
                                    bSecondNonTautPass, bOmitRepetitions, bUseMulipliers);
            if ( str_LineEnd( szTag2, tot_len, nStrLen, &bOverflow, pStr, bXml? 0 : -1, bPlainTextTags ) )
                goto exit_function;
            inchi_ios_print( output_file, "%s%s", pStr, pLF );

            /*************************************
              inverted sp3  canonical numbering
            **************************************/
            if ( bInvStereoOrigNumb[iCurTautMode] )
            {
                szGetTag( AuxLbl, nTag,  bTag2 = bTag1 | AL_SP3N, szTag2, &bAlways );
                tot_len = str_LineStart( szTag2, NULL, 0, pStr, ind );
                tot_len = str_AuxInvSp3Numb(pINChISort, pINChISort2, pStr, nStrLen, tot_len,
                                            &bOverflow, bOutType, TAUT_MODE, num_components,
                                            bSecondNonTautPass, bOmitRepetitions);
                if ( str_LineEnd( szTag2, tot_len, nStrLen, &bOverflow, pStr, bXml? 0 : -1, bPlainTextTags ) )
                    goto exit_function;
                inchi_ios_print( output_file, "%s%s", pStr, pLF );
            }
            else
            {
                if ( !bXml )
                {
                    if ( bPlainTextTags == 1 ) inchi_ios_print( output_file, "/" );
                }
            }

            if ( bXml )
            {
                /* close sp3 inv */
                ind -= inc;
                if ( str_LineEnd( szTag1, 0, nStrLen, &bOverflow, pStr, ind, bPlainTextTags ) )
                    goto exit_function;
                inchi_ios_print( output_file, "%s%s", pStr, pLF );
            }
        }
        else
        {
            if ( !bXml )
            {
                if ( bPlainTextTags == 1 )
                    inchi_ios_print( output_file, "//" );
            } /* Inverted stereo -- sp3 only + canonical numbering */
        }


        /* omitted undefined/unknown non-isotopic stereo */
        if ( bXml )
        {
            if ( bIgn_UU_Sp2[iCurTautMode] || bIgn_UU_Sp3[iCurTautMode] )
            {
                /* <stereo omit_undef_dbond="1" omit_undef_sp3="1"/> */
                szGetTag( IdentLbl, nTag,  bTag1 = IL_STER, szTag1, &bAlways );
                tot_len = PrintXmlStartTag( pStr, ind, 3, szTag1,
                                 (bIgn_UU_Sp2[iCurTautMode])? x_ign_uu_sp2 : NULL, 1,
                                 (bIgn_UU_Sp3[iCurTautMode])? x_ign_uu_sp3 : NULL, 1,
                                 NULL, 0, NULL, 0, NULL, 0, NULL, 0 );
                inchi_ios_print( output_file, "%s\n", pStr );
            }
        }

        /***************************************************************
         *
         *  Additional information: charges, radicals,
         *                          special valences, coordinates
         *
         ***************************************************************/
        /**************************************************************
         *
         *   Aux info isotopic
         *
         **************************************************************/
repeat_INChI_Aux_Iso_output:
        /* if InChI Fixed-H isotopic is empty then do not output corresponding AuxInfo */
        i =  bSecondNonTautPass &&
             (*pSortPrintINChIFlags & ((INCHI_BAS == iINChI)? FLAG_SORT_PRINT_NO_IFIX_H_BAS :
                                                              FLAG_SORT_PRINT_NO_IFIX_H_REC ));

        if ( bIsotopic && !i &&
                          (bIsotopicOrigNumb[iCurTautMode] ||
                           bIsotopicAtomEqu[iCurTautMode] ||
                           (bTautomericOutputAllowed && bTautomeric && bIsotopicTautEqu[iCurTautMode]) ||
                           bInvIsotopicStereo[iCurTautMode] ||
                           (bXml && ( bIgn_UU_Sp3_Iso[iCurTautMode] || bIgn_UU_Sp2_Iso[iCurTautMode] )) ) )
        {
            /*************************************/
            /*   isotopic aux info header        */
            /*************************************/
            szGetTag( AuxLbl, nTag,  bTag1 = AL_ISOT | bFhTag, szTag1, &bAlways );
            if ( bXml )
            {
                str_LineStart( szTag1, NULL, 0, pStr, ind );
                inchi_ios_print( output_file, "%s\n", pStr );
                ind += inc;
            }
            else
            {
                pStr[tot_len = 0] = '\0';
            }
            /*****************************************************************
             *   Original atom numbers in order of isotopic canonical numbers
             *****************************************************************/
            szGetTag( AuxLbl, nTag,  bTag2 = bTag1 | AL_ISON, szTag2, &bAlways );
            if ( bIsotopicOrigNumb[iCurTautMode] )
            {
                tot_len = str_LineStart( szTag2, NULL, 0, pStr, ind );
                tot_len = str_AuxIsoNumb(pINChISort, pINChISort2, pStr, nStrLen, tot_len,
                                         &bOverflow, bOutType, TAUT_MODE, num_components,
                                         bSecondNonTautPass, bOmitRepetitions);
                if ( str_LineEnd( szTag2, tot_len, nStrLen, &bOverflow, pStr, bXml? 0 : -1, bPlainTextTags ) )
                    goto exit_function;
                inchi_ios_print( output_file, "%s%s", pStr, pLF );
            }
            else
            {
                if ( !bXml )
                {
                    /*if ( bPlainTextTags == 1 ) inchi_ios_print( output_file, "/" );*/
                    inchi_ios_print( output_file, "%s%s", szTag2, pLF ); /* mark isotopic output */
                }
            }

            /*************************/
            /*  Isotopic symmetry    */
            /*************************/
            if ( bIsotopicAtomEqu[iCurTautMode] )
            {
                /*  atoms */
                szGetTag( AuxLbl, nTag,  bTag2 = bTag1 | AL_AEQU, szTag2, &bAlways );
                tot_len = str_LineStart( szTag2, NULL, 0, pStr, ind );
                tot_len = str_AuxIsoEqu(pINChISort, pINChISort2, pStr, nStrLen, tot_len,
                                        &bOverflow, bOutType, TAUT_MODE, num_components,
                                        bSecondNonTautPass, bOmitRepetitions, bUseMulipliers);
                if ( str_LineEnd( szTag2, tot_len, nStrLen, &bOverflow, pStr, bXml? 0 : -2/*was -1: Fix15*/, bPlainTextTags ) )
                    goto exit_function;
                inchi_ios_print( output_file, "%s%s", pStr, pLF );
            }
            else
            {
                if ( !bXml )
                {
                    if ( bPlainTextTags == 1 )
                        inchi_ios_print( output_file, "/" );
                }
            }

            /********************************/
            /*  Tautomeric groups, isotopic */
            /********************************/
            if ( bTautomericOutputAllowed && bTautomeric && bIsotopicTautEqu[iCurTautMode] )
            {
                /********************************************/
                /*  Isotopic tautomeric groups equivalence */
                /********************************************/
                szGetTag( AuxLbl, nTag,  bTag2 = bTag1 | AL_GEQU, szTag2, &bAlways );
                tot_len = str_LineStart( szTag2, NULL, 0, pStr, ind );
                tot_len = str_AuxIsoTgroupEqu(pINChISort, pStr, nStrLen, tot_len,
                                              &bOverflow, bOutType, TAUT_MODE, num_components,
                                              bOmitRepetitions, bUseMulipliers);
                if ( str_LineEnd( szTag2, tot_len, nStrLen, &bOverflow, pStr, bXml? 0 : -2/*was -1: Fix15*/, bPlainTextTags ) )
                    goto exit_function;
                inchi_ios_print( output_file, "%s%s", pStr, pLF );
            }
            else
            {
                if ( !bXml && bTautomericOutputAllowed && bTautomeric )
                {
                    if ( bPlainTextTags == 1 )
                        inchi_ios_print( output_file, "/" );
                }
            }

            /*************************************
             * Isotopic inverted stereo
             *************************************/
            if ( bInvIsotopicStereo[iCurTautMode] )
            {
                szGetTag( AuxLbl, nTag,  bTag2 = bTag1 | AL_STER, szTag2, &bAlways );
                if ( bXml )
                {
                    /************************************
                         inv isotopic stereo start  tag
                    *************************************/
                    str_LineStart( szTag2, NULL, 0, pStr, ind );
                    inchi_ios_print( output_file, "%s\n", pStr );
                    ind += inc;
                }
                /*************************************
                     inverted isotopic sp3 start tag
                **************************************/
                szGetTag( AuxLbl, nTag,  bTag3 = bTag2 | AL_SP3I, szTag3, &bAlways );
                tot_len = str_LineStart( szTag3, NULL, 0, pStr, ind );
                tot_len = str_AuxInvIsoSp3(pINChISort, pINChISort2, pStr, nStrLen, tot_len,
                                           &bOverflow, bOutType, TAUT_MODE, num_components,
                                           bSecondNonTautPass, bOmitRepetitions, bUseMulipliers);
                if ( str_LineEnd( szTag3, tot_len, nStrLen, &bOverflow, pStr, bXml? 0 : -1, bPlainTextTags ) )
                    goto exit_function;
                inchi_ios_print( output_file, "%s", pStr );
                /*********************************************
                  inverted isotopic sp3  canonical numbering
                **********************************************/
                if ( bInvIsotopicStereoOrigNumb[iCurTautMode] )
                {
                    szGetTag( AuxLbl, nTag,  bTag3 = bTag2 | AL_SP3N, szTag3, &bAlways );
                    tot_len = str_LineStart( szTag3, NULL, 0, pStr, ind );
                    tot_len = str_AuxInvIsoSp3Numb(pINChISort, pINChISort2, pStr, nStrLen, tot_len,
                                                   &bOverflow, bOutType, TAUT_MODE, num_components,
                                                   bSecondNonTautPass, bOmitRepetitions);
                    if ( str_LineEnd( szTag3, tot_len, nStrLen, &bOverflow, pStr, bXml? 0 : -1, bPlainTextTags ) )
                        goto exit_function;
                    inchi_ios_print( output_file, "%s%s", pStr, pLF );
                }
                else
                {
                    if ( !bXml )
                    {
                        if ( bPlainTextTags == 1 )
                            inchi_ios_print( output_file, "/" );
                    }
                }
                if ( bXml )
                {
                /* close sp3 inv */
                    ind -= inc;
                    if ( str_LineEnd( szTag2, 0, nStrLen, &bOverflow, pStr, ind, bPlainTextTags ) )
                        goto exit_function;
                    inchi_ios_print( output_file, "%s", pStr );
                }
            }
            else
            {
                if ( !bXml )
                {
                    if ( bPlainTextTags == 1 )
                        inchi_ios_print( output_file, "//" );
                }
            }

            /* totally omitted undefined/unknown isotopic stereo */
            if ( bXml )
            {
                if ( bIgn_UU_Sp3_Iso[iCurTautMode] || bIgn_UU_Sp2_Iso[iCurTautMode]  )
                {
                    /* <stereo omit_undef_dbond="1" omit_undef_sp3="1"/> */
                    szGetTag( IdentLbl, nTag,  bTag1 = IL_STER, szTag1, &bAlways );
                    tot_len = PrintXmlStartTag( pStr, ind, 3, szTag1,
                                     (bIgn_UU_Sp2_Iso[iCurTautMode])? x_ign_uu_sp2 : NULL, 1,
                                     (bIgn_UU_Sp3_Iso[iCurTautMode])? x_ign_uu_sp3 : NULL, 1,
                                     NULL, 0, NULL, 0, NULL, 0, NULL, 0 );
                    inchi_ios_print( output_file, "%s\n", pStr );
                }
            }


            if ( bXml )
            {
                /*****************  close isotopic ***********************/
                ind -= inc;
                if ( str_LineEnd( szTag1, 0, nStrLen, &bOverflow, pStr, ind, bPlainTextTags ) )
                    goto exit_function;
                inchi_ios_print( output_file, "%s", pStr );
            }
        } /* Aux info isotopic */


#if ( CANON_FIXH_TRANS != 1 )
        if ( bSecondNonTautPass )
        {
            /* find and print non-tautomeric components transposition, if non-trivial */
            AT_NUMB *nTrans_n, *nTrans_s;
            if ( 0 < bin_AuxTautTrans(pINChISort, pINChISort2, &nTrans_n, &nTrans_s, bOutType,  num_components) ) {
                /* a non-trivial transposition does exist; output start tag */
                tot_len = str_LineStart( tag=x_aux_trans, NULL, 0, pStr, ind );
                /* print the transposition, cycle after cycle */
                tot_len = str_AuxTautTrans(nTrans_n, nTrans_s, pStr, nStrLen, tot_len,
                                           &bOverflow, TAUT_MODE, num_components);
                if ( str_LineEnd( bXml? tag:p_aux_at_inv_nbr, tot_len, nStrLen, &bOverflow, pStr, bXml? 0 : -1, bPlainTextTags ) )
                    goto exit_function;
                inchi_ios_print( output_file, "%s", pStr );
                /* detected transposition */
                *pSortPrintINChIFlags |= (INCHI_BAS == iINChI)? FLAG_SORT_PRINT_TRANSPOS_BAS :
                                                                FLAG_SORT_PRINT_TRANSPOS_REC;
            }  else
            if ( !bXml ) {
                if ( bPlainTextTags == 1 ) inchi_ios_print( output_file, "/" );
            }
        }
#endif

        /**************************************************************
          At this point the INChI_Aux part of the output has been completed.
          If this INChI is tautomeric and non-tautomeric results exist
          then we need to output non-tautomeric auxilialy data
          (same as above excluding tautomeric information)
          Currently this is enabled for xml output only
        ***************************************************************/

        if ( bOutType == OUT_TN && bTautomeric && bNonTautomeric &&
            /* Check whether the Fixed-H layer is empty */
            (*pSortPrintINChIFlags & ((INCHI_BAS == iINChI)? FLAG_SORT_PRINT_NO_NFIX_H_BAS :
                                                             FLAG_SORT_PRINT_NO_NFIX_H_REC )) &&
            (*pSortPrintINChIFlags & ((INCHI_BAS == iINChI)? FLAG_SORT_PRINT_NO_IFIX_H_BAS :
                                                             FLAG_SORT_PRINT_NO_IFIX_H_REC ))
              )
        {
            bNonTautomeric = 0; /* bNonTautIdentifierNotEmpty == 0 => no fixed H info 02-10-2995 */
        }

        if ( bOutType == OUT_TN && bTautomeric && bNonTautomeric )
        {
            /* add the second (non-tautomeric) output */
            bOutType     = OUT_NONTAUT;
            iCurTautMode = TAUT_NON;
            pINChISort    = pINChISortTautAndNonTaut[TAUT_NON];
            bSecondNonTautPass = 1;
            num_components = num_comp[iCurTautMode];
            bFhTag = AL_FIXH;
            if ( bXml )
            {
                szGetTag( AuxLbl, nTag,  bTag1 = bFhTag, szTag1, &bAlways );
                str_LineStart( szTag1, NULL, 0, pStr, ind );
                inchi_ios_print( output_file, "%s\n", pStr );
                ind += inc;
            }
            else
            {
                pStr[tot_len=0] = '\0';
            }

            /* if InChI Fixed-H isotopic is empty then do not output corresponding AuxInfo */
            if ( !(*pSortPrintINChIFlags &
                    ((INCHI_BAS == iINChI)? FLAG_SORT_PRINT_NO_NFIX_H_BAS :
                                            FLAG_SORT_PRINT_NO_NFIX_H_REC ))
               )
            {
                goto repeat_INChI_Aux_output;
            }
            else
            {
                goto repeat_INChI_Aux_Iso_output;
            }
        }
        else
        {
            if ( bOutType == OUT_NONTAUT && bOutputType == OUT_TN && bTautomeric && bNonTautomeric )
            {
                /* the second (non-taut) output has been done; restore variables */
                bOutType           = OUT_TN;
                iCurTautMode       = TAUT_YES;
                pINChISort         = pINChISortTautAndNonTaut[TAUT_YES];
                bSecondNonTautPass = 0;
                /* set correct num components for the reversibility info 02-10-2005 */
                num_components     = num_comp[iCurTautMode];
                if ( bXml )
                {
                    /*  close non-tautomeric */
                    szGetTag( AuxLbl, nTag,  bTag1 = bFhTag, szTag1, &bAlways );
                    ind -= inc;
                    if ( str_LineEnd( szTag1, 0, nStrLen, &bOverflow, pStr, ind, bPlainTextTags ) )
                        goto exit_function;
                    inchi_ios_print( output_file, "%s", pStr );
                }
                bFhTag = 0;
            }
        }


        /***************************************/
        /* charges, radicals, unusual valences */
        /***************************************/
        if ( !bSecondNonTautPass && bChargesRadVal[iCurTautMode] )
        {
            /*  aux equ atoms */
            /* 1. Compare to tautomeric equivalence (in case of second, non-taut, pass only) */
            /* 2. Compare to the previous component if (1) failed to find equivalence */
            szGetTag( AuxLbl, nTag,  bTag1 = AL_CRV_ | bFhTag, szTag1, &bAlways );
            tot_len = str_LineStart( szTag1, NULL, 0, pStr, ind );
            tot_len = str_AuxChargeRadVal(pINChISort, pStr, nStrLen, tot_len,
                                          &bOverflow, bOutType, TAUT_MODE,
                                          num_components, bUseMulipliers);
            if ( str_LineEnd( szTag1, tot_len, nStrLen, &bOverflow, pStr, bXml? 0 : -1, bPlainTextTags ) )
                goto exit_function;
            inchi_ios_print( output_file, "%s%s", pStr, pLF );
        }

        /* output the original input structure -- quick fix */
        if ( !bSecondNonTautPass && pOrigStruct && pOrigStruct->num_atoms &&
             pOrigStruct->szAtoms && pOrigStruct->szBonds && pOrigStruct->szCoord )
        {
            int length, cur_pos, line_len, last_pos, nMaxLineLen;
            char *p;
            nMaxLineLen = inchi_min( 80, nStrLen ); /* restrict line length to 80 characters */
            /**********************
               reversibility info
             **********************/
            szGetTag( AuxLbl, nTag,  bTag1 = AL_REVR | bFhTag, szTag1, &bAlways );
            if ( bXml )
            {
                str_LineStart( szTag1, NULL, 0, pStr, ind );
                inchi_ios_print( output_file, "%s\n", pStr );
                ind += inc;
            }
            /*  === atoms === */
            szGetTag( AuxLbl, nTag,  bTag2 = bTag1 | AL_ATMR, szTag2, &bAlways );
            if ( bXml )
            {
                str_LineStart( szTag2, NULL, 0, pStr, ind );
                inchi_ios_print( output_file, "%s\n", pStr );
                ind += inc;
                /* first line indent */
                strcpy( pStr, SP(ind));
                tot_len = ind;
            }
            else
            {
                pStr[tot_len = 0] = '\0';
                inchi_ios_print( output_file, "%s%s", szTag2, pStr );
            }
            p = pOrigStruct->szAtoms;
            length = strlen( p );
            line_len = nMaxLineLen - tot_len;
            for ( cur_pos = 0; cur_pos < length; cur_pos = last_pos )
            {
                if ( length - cur_pos >= line_len )
                {
                    last_pos = cur_pos + line_len;
                    /* search backward for the nearest first atom letter (always uppercase) */
                    while ( cur_pos < last_pos && !isupper( UCINT p[last_pos] ) ) {
                        last_pos --;
                    }
                }
                else
                {
                    last_pos = length;
                }
                if ( last_pos > cur_pos )
                {
                    memcpy( pStr + tot_len, p+cur_pos, last_pos - cur_pos );
                    pStr[tot_len + last_pos - cur_pos] = '\0';
                    inchi_ios_print( output_file, "%s%s", pStr, !bXml && bPlainTextTags? "" : "\n" );
                }
                else
                {
                    break;
                }
            }
            if ( bXml )
            {
                ind -= inc;
                pStr[0] = '\0';
                if ( str_LineEnd( szTag2, 0, nMaxLineLen, &bOverflow, pStr, ind, bPlainTextTags ) )
                    goto exit_function;
                inchi_ios_print( output_file, "%s", pStr );
            }
            else
            {
                if ( pLF[0] )
                {
                    inchi_ios_print( output_file, "%s", pLF );
                }
            }


            /*  === bonds === */
            szGetTag( AuxLbl, nTag,  bTag2 = bTag1 | AL_BNDR, szTag2, &bAlways );
            if ( bXml )
            {
                str_LineStart( szTag2, NULL, 0, pStr, ind );
                inchi_ios_print( output_file, "%s\n", pStr );
                ind += inc;
                /* first line indent */
                strcpy( pStr, SP(ind));
                tot_len = ind;
            }
            else
            {
                pStr[tot_len = 0] = '\0';
                inchi_ios_print( output_file, "%s%s", szTag2, pStr );
            }

            p = pOrigStruct->szBonds;
            length = strlen( p );
            line_len = nMaxLineLen - tot_len;
            for ( cur_pos = 0; cur_pos < length; cur_pos = last_pos )
            {
                if ( length - cur_pos >= line_len )
                {
                    last_pos = cur_pos + line_len - 1;
                    /* search backward for the nearest first bond delimiter ";" */
                    while ( cur_pos < last_pos && p[last_pos] != ';' )
                    {
                        last_pos --;
                    }
                    if ( cur_pos < last_pos )
                    {
                        last_pos ++; /* include ';' at the end of the line */
                    }
                }
                else
                {
                    last_pos = length;
                }
                if ( last_pos > cur_pos )
                {
                    memcpy( pStr + tot_len, p+cur_pos, last_pos - cur_pos );
                    pStr[tot_len + last_pos - cur_pos] = '\0';
                    inchi_ios_print( output_file, "%s%s", pStr, !bXml && bPlainTextTags? "" : "\n" );
                }
                else
                {
                    break;
                }
            }

            if ( bXml )
            {
                ind -= inc;
                pStr[0] = '\0';
                if ( str_LineEnd( szTag2, 0, nMaxLineLen, &bOverflow, pStr, ind, bPlainTextTags ) )
                    goto exit_function;
                inchi_ios_print( output_file, "%s", pStr );
            }
            else
            {
                if ( pLF[0] )
                {
                    inchi_ios_print( output_file, "%s", pLF );
                }
            }


            /*  === coordinates === */
            szGetTag( AuxLbl, nTag,  bTag2 = bTag1 | AL_XYZR, szTag2, &bAlways );
            if ( bXml )
            {
                str_LineStart( szTag2, NULL, 0, pStr, ind );
                inchi_ios_print( output_file, "%s\n", pStr );
                ind += inc;
                /* first line indent */
                strcpy( pStr, SP(ind));
                tot_len = ind;
            }
            else
            {
                pStr[tot_len = 0] = '\0';
                inchi_ios_print( output_file, "%s%s", szTag2, pStr );
            }

            p = pOrigStruct->szCoord;
            length = strlen( p );
            line_len = nMaxLineLen - tot_len;
            for ( cur_pos = 0; cur_pos < length; cur_pos = last_pos )
            {
                if ( length - cur_pos >= line_len )
                {
                    last_pos = cur_pos + line_len - 1;
                    /* search backward for the nearest first coord. delimiter ";" */
                    while ( cur_pos < last_pos && p[last_pos] != ';' )
                    {
                        last_pos --;
                    }
                    if ( cur_pos < last_pos )
                    {
                        last_pos ++; /* include ';' at the end of the line */
                    }
                }
                else
                {
                    last_pos = length;
                }
                if ( last_pos > cur_pos )
                {
                    memcpy( pStr + tot_len, p+cur_pos, last_pos - cur_pos );
                    pStr[tot_len + last_pos - cur_pos] = '\0';
                    inchi_ios_print( output_file, "%s%s", pStr, !bXml && bPlainTextTags? "" : "\n" );
                }
                else
                {
                    break;
                }
            }

            if ( bXml )
            {
                ind -= inc;
                pStr[0] = '\0';
                if ( str_LineEnd( szTag2, 0, nMaxLineLen, &bOverflow, pStr, ind, bPlainTextTags ) )
                    goto exit_function;
                inchi_ios_print( output_file, "%s", pStr );
            }
            else
            {
                if ( pLF[0] )
                {
                    inchi_ios_print( output_file, "%s", pLF );
                }
            }

            if ( bXml )
            {
                /***************************
                  close reversibility info
                 ***************************/
                ind -= inc;
                if ( str_LineEnd( szTag1, 0, nStrLen, &bOverflow, pStr, ind, bPlainTextTags ) )
                    goto exit_function;
                inchi_ios_print( output_file, "%s", pStr );
            }
        }



        /************************************************
         * output INChI_Aux of the reconnected structure *
         ************************************************/
        bEmbeddedOutputCalled = 0;
        if ( bDisconnectedCoord && INCHI_BAS == iINChI && (bINChIOutputOptions & INCHI_OUT_EMBED_REC) &&
             num_components2[INCHI_REC] && !(bINChIOutputOptions & INCHI_OUT_NO_AUX_INFO) )
        {
            int nRet;
            bEmbeddedOutputCalled = 1;
            if ( !bXml )
            {
                inchi_ios_print( output_file, "%s", pLF );
            }

            nRet = OutputINChI1(pStr, nStrLen,
                                pINChISortTautAndNonTaut2,
                                INCHI_REC,
                                NULL,
                                0 /*bDisconnectedCoord*/, bOutputType,
                                INCHI_OUT_ONLY_AUX_INFO | bINChIOutputOptions,
                                bXml, bAbcNumbers, bCtPredecessors, bNoStructLabels,
                                num_components2,
                                num_non_taut2, num_taut2,
                                output_file, log_file,
                                num_input_struct,
                                szSdfLabel, szSdfValue, lSdfId,
                                pSortPrintINChIFlags,
                                save_opt_bits);
            if ( !nRet )
                goto exit_function; /* error */
        }

        /* close INChI_Aux */
        if ( bXml )
        {
                ind -= inc;
            if ( str_LineEnd( x_aux_basic, 0, nStrLen, &bOverflow, pStr, ind, bPlainTextTags ) )
                goto exit_function;
            inchi_ios_print( output_file, "%s", pStr );
        }
        else
        {
            if ( !bEmbeddedOutputCalled && !bPlainTextCommnts )
            {
                inchi_ios_print( output_file, "%s\n", (!num_components2[0] && !num_components2[1])? "//":"" );
                /* plain text comment earlier ended with LF */
            }
        }


        /* in wINChI window, separate AuxInfo: from InChIKey: with blank line */
        inchi_ios_print( output_file, "%s",
                        (bINChIOutputOptions & INCHI_OUT_WINCHI_WINDOW) ? "\n":"");




    } /* end of output aux info */


    ret = 1;
exit_function:

    if ( bOverflow )
    {
        strcpy( pStr, "Output buffer overflow");
        if ( bXml )
        {
            OutputINChIXmlError( output_file, pStr, nStrLen, ind /*, 0*/ /* err number */, pStr, _IS_FATAL );
        }
        else
        {
            inchi_ios_print( output_file, "\nFATAL ERROR: %s\n", pStr );
        }
    }

    /* inchi_free( pStr ); */
    return ret;


} /* OutputINChI1 */



/***************************************************************/
int str_LineStart( const char *tag, char *tag2, int val2, char *pStr, int ind )
{
    int tot_len = 0;
    if ( ind >= 0 ) {
        if ( ind > 0 ) {
            /* xml: indent */
            memset( pStr + tot_len, ' ', ind );
            tot_len += ind;
        }
        /* xml: tag */
        strcpy( pStr + tot_len, x_line_opening );
        strcat( pStr + tot_len, tag );
        if ( tag2 ) {
            tot_len += strlen(pStr + tot_len);
            tot_len += sprintf( pStr + tot_len, " %s=\"%d\"%s", tag2, val2, x_close_line );

        } else {
            strcat( pStr + tot_len, x_close_line );
            tot_len += strlen(pStr + tot_len);
        }
    } else {
        pStr[tot_len] = '\0';
    }
    return tot_len;
}
/***************************************************************/
int str_LineEnd( const char *tag, int tot_len, int nStrLen, int *bOverflow, char *pStr, int ind, int bPlainTextTags )
{
    static int  add_tag_len = sizeof(x_line_closing)-1 + sizeof(x_close_line)-1;
    int tag_len;
    /* check buffer overflow */
    if ( *bOverflow )
        return 1;
    if ( ind >= 0 ) {  /* xml */
        tag_len = ind + add_tag_len + strlen(tag);
        if ( tot_len + tag_len < nStrLen - 2  ) {
            /* output "   </tag>\n" */
            tot_len += sprintf( pStr + tot_len, "%s%s%s%s\n", SP(ind), x_line_closing, tag, x_close_line );
        } else {
            *bOverflow += 1;
            return 1;
        }
    } else { /* plain */
        pStr[tot_len] = '\0'; /* add zero termination 2004-04-26 */
        /* insert plain text tag if:
           (a) pStr has non-zero length, or
           (b) ind < -1
        */

        if ( pStr[0] || ind < -1 ) {
            tag_len = bPlainTextTags? strlen( tag ):0;
            if ( tot_len + tag_len  < nStrLen - 2 ) {
                if ( tag_len > 0 ) {
                    /* insert plain text tag */
                    memmove( pStr+tag_len, pStr, tot_len + 1 );
                    memcpy( pStr, tag, tag_len );
                }
            } else {
                *bOverflow += 1;
                return 1;
            }
        }/* else
        if ( bPlainTextTags == 1 ) {
            strcpy( pStr, "/" );
        }*/
    }
    return 0;
}



/**********************************************************************************************/
int CleanOrigCoord( MOL_COORD szCoord, int delim )
{
#define MIN_BOND_LENGTH   (1.0e-6)
    char szVal[LEN_COORD+1];
    MOL_COORD szBuf;
    char *q;
    int len, last, fst, dec_pnt, num_zer=0, len_buf = 0, e;
    int k, i;
    double coord;

    for ( k = 0; k < NUM_COORD*LEN_COORD; k += LEN_COORD ) {
        memcpy( szVal, szCoord+k, LEN_COORD );
        szVal[LEN_COORD] = '\0';
        LtrimRtrim(szVal, &len);
        coord = strtod(szVal, &q);
        if ( MIN_BOND_LENGTH > fabs(coord)  ) {
            strcpy( szVal, "0" );
            len = 1;
            num_zer ++;
        } else {
            len = q - szVal;
            /* last = (last mantissa digit position + 1)  */
            if ( (q = strchr(szVal, 'e')) || (q = strchr(szVal, 'E')) ||
                 (q = strchr(szVal, 'd')) || (q = strchr(szVal, 'D')) ) {
                /* floating point */
                last = q - szVal;
                /* remove (+) and leading zeroes from the exponent */
                e = (int)strtol( szVal+last+1, &q, 10 ); /* exponent */
                if ( e ) {
                    /* new exp; update the length */
                    len = last+1+sprintf( szVal+last+1, "%d", e ); /* print exp without leading zeroes and '+' */
                } else {
                    /* exponent is zero */
                    len = last;
                }
            } else {
                last = len;
            }
            /* fst = (first mantissa digit); fst=1 if the sign is present, otherwise 0 */
            fst = (szVal[0]!='.' && !isdigit( UCINT szVal[0] ));
            /* dec_pnt = (decimal point position) or last */
            if ( (q = strchr(szVal, '.')) ) {
                dec_pnt = q - szVal;
            } else {
                dec_pnt = last;
            }
            last -= 1; /* last mantissa digit position */
            /* remove trailing zeroes in the range dec_pnt+1..last-1 */
            for ( i = last; dec_pnt < i &&  '0' == szVal[i]; i -- )
                ;
            if ( i == dec_pnt ) {
                i --; /* remove decimal point, too */
            }
            if ( i < last ) {
                memmove( szVal+i+1, szVal+last+1, len-last );
                len -= last-i;
            }
            /* remove leading zeroes */
            for ( i = fst; i < len && '0' == szVal[i]; i ++ )
                ;
            if ( i > fst ) {
                memmove( szVal + fst, szVal+i, len-fst );
                len -= i-fst;
            }
        }
        if ( len_buf )
            szBuf[len_buf++] = delim;
        memcpy( szBuf + len_buf, szVal, len ); /* does not copy zero termination*/
        len_buf += len;
    }
    /* zero termination */
    if ( len_buf < (int)sizeof(MOL_COORD) ) {
        memset( szBuf+len_buf, 0, sizeof(MOL_COORD) - len_buf);
    }
    memcpy( szCoord, szBuf, sizeof(MOL_COORD) );
    return num_zer;
#undef MIN_BOND_LENGTH
}

/******************************************************************************************/
int WriteOrigCoord( int num_inp_atoms, MOL_COORD *szMolCoord, int *i, char *szBuf, int buf_len )
{

    int j, num_zer, len, cur_len;
    char *p;
    MOL_COORD szCurCoord;
    cur_len = 0;
    for ( j = *i; j < num_inp_atoms; ) {
        memcpy( szCurCoord, szMolCoord[j], sizeof(szCurCoord));
        num_zer = CleanOrigCoord( szCurCoord, ',' );
        if ( NUM_COORD == num_zer ) {
            len = 0;
        } else {
            if ( (p = (char *)memchr( szCurCoord, '\0', sizeof(szCurCoord))) ) {
                len = p - szCurCoord;
            } else {
                len = sizeof(szCurCoord);
            }
        }
        if ( len + cur_len + 1 < buf_len ) {
            if ( len ) {
                memcpy( szBuf + cur_len, szCurCoord, len * sizeof(szBuf[0]) );
            }
            szBuf[cur_len += len] = ';';
            cur_len ++;
            j ++;
        } else {
            break;
        }
    }
    szBuf[cur_len] = '\0';
    *i = j; /* next item */
    return cur_len;
}
/******************************************************************************************/
/*
  number of atoms
  [c|n]              chiral/nonchiral

  Element
  #valence
  +/-[charge>1]
  .#rad  (#rad=1, 2, 3: singlet, doulet, triplet)
  [.]i#iso_mass
  [.]{o|e|u|?} atom parity = {1:2:3:4}
  [.]h[#of 1H>1]
  [.]d[#of 2H>1]
  [.]t[#of 3H>1]

  Note: . occurs only once and only if radical or 1-character element
*/
int WriteOrigAtoms( int num_inp_atoms, inp_ATOM *at, int *i, char *szBuf, int buf_len, STRUCT_DATA *sd)
{
    int j, k, n, len, len0, cur_len, val, bonds_val, mw, parity, num_trans, is_ok, b_self;
    static char szIsoH[] = "hdt";
    char szCurAtom[32];
    AT_NUMB nNeighOrder[MAXVAL], neigh;

    cur_len = 0;
    if ( 0 == *i ) {
        cur_len = sprintf( szBuf, "%d%s", num_inp_atoms,
                                  (sd->bChiralFlag & FLAG_INP_AT_CHIRAL)?    "c" :
                                  (sd->bChiralFlag & FLAG_INP_AT_NONCHIRAL)? "n" : "" );
    }
    for ( j = *i; j < num_inp_atoms; ) {
        /* tetrahedral parity treatment */
        parity    = 0;
        num_trans = 0;
        if ( at[j].p_parity ) {
            /* verify neighbors */
            is_ok  = 1;
            b_self = 0;
            for ( n = 0, k = 0; n < MAX_NUM_STEREO_ATOM_NEIGH; n ++ ) {
                neigh = at[j].p_orig_at_num[n]-1;
                if ( is_in_the_list( at[j].neighbor, neigh, at[j].valence ) &&
                     at[neigh].orig_at_number  ==  at[j].p_orig_at_num[n] ) {
                    /* real neighbor */
                    nNeighOrder[k ++] = at[j].p_orig_at_num[n];
                } else
                if ( (int)neigh == j && at[neigh].orig_at_number  ==  at[j].p_orig_at_num[n] ) {
                    /* central atom is a neighbor */
                    num_trans = n; /* move this neighbor to 0 position permutation parity */
                    b_self ++;
                } else {
                    is_ok = 0;
                    break;
                }
            }
            if ( is_ok && b_self <= 1 && b_self + k == MAX_NUM_STEREO_ATOM_NEIGH ) {
                num_trans += insertions_sort( nNeighOrder, k, sizeof(nNeighOrder[0]), comp_AT_RANK );
                if ( ATOM_PARITY_WELL_DEF( at[j].p_parity ) ) {
                    parity = 2 - (num_trans + at[j].p_parity) % 2;
                } else
                if ( ATOM_PARITY_ILL_DEF( at[j].p_parity ) ) {
                    parity = at[j].p_parity;
                } else {
                    ; /* invalid atom parity */
                }
            } else {
                ;/* add error message here */
            }
        }

        len = len0 = strlen( at[j].elname );
        memcpy( szCurAtom, at[j].elname, len );
        bonds_val = nBondsValenceInpAt( at+j, NULL, NULL );

        if ( (val=needed_unusual_el_valence( at[j].el_number, at[j].charge, at[j].radical,
                                 at[j].chem_bonds_valence, bonds_val, at[j].num_H, at[j].valence )) ||
             at[j].charge || at[j].radical || at[j].iso_atw_diff || NUM_ISO_H(at,j) || parity ) {
            /* valence */
            if ( val ) {
                len += sprintf( szCurAtom + len, "%d", val > 0? val : 0 );
            }
            /* charge */
            if ( (val = at[j].charge) ) {
                szCurAtom[len++] = val>0? '+' : '-';
                if ( (val = abs(val)) > 1 ) {
                    len += sprintf( szCurAtom + len, "%d", val );
                }
            }
            /* radical */
            if ( (val = at[j].radical) ) {
                len += sprintf(szCurAtom + len, ".%d", val);
            }
            /* isotopic shift */
            if ( (val = at[j].iso_atw_diff) ) {
                mw = get_atw_from_elnum( at[j].el_number );
                if ( val == 1 )
                    val = mw;
                else
                if ( val > 0 )
                    val = mw + val -1;
                else
                    val = mw + val;
                len += sprintf( szCurAtom + len, "%si%d", len == len0? ".":"", val );
            }
            /* parity */
            if ( parity ) {
                len += sprintf( szCurAtom + len, "%s%s", len == len0? ".":"",
                                parity == AB_PARITY_ODD?   "o" :
                                parity == AB_PARITY_EVEN?  "e" :
                                parity == AB_PARITY_UNKN?  "u" :
                                parity == AB_PARITY_UNDF?  "?" : "" );
            }
            /* implicit isotopic H */
            if ( NUM_ISO_H(at,j) ) {
                for ( k = 0; k < NUM_H_ISOTOPES; k ++ ) {
                    if ( (val = at[j].num_iso_H[k]) ) {
                        len += sprintf( szCurAtom + len, "%s%c", len == len0? ".":"", szIsoH[k] );
                        if ( val > 1 ) {
                            len += sprintf(szCurAtom + len, "%d", val);
                        }
                    }
                }
            }
        }
        if ( len + cur_len < buf_len ) {
            memcpy( szBuf + cur_len, szCurAtom, len );
            cur_len += len;
            j ++;
        } else {
            break;
        }
        szBuf[cur_len] = '\0';
        *i = j;

    }
    return cur_len;
}
/******************************************************************************************/
/*
 <bonds> bpA;bpAbpA... </bonds>

b = bond type:
=============
w = undefined stereo, double
s = single
d = double
t = triple
a = aromatic
p = up from the current atom to the neighbor
P = uP from the neighbor to the current atom
v = undefined stereo Either, single from the current atom to the neighbor
V = undefined stereo Either, single from the neighbor to the current atom
n = down from the current atom to the neighbor
N = dowN from the neighbor to the current atom

p = bond parity:
================
- = odd
+ = even
u = unknown
? = undefined
  = no parity (empty)


A = neighbor orig. atom number
===============
neighbor orig. atom number < number of the current atom
Number of the current atom: 2 until first ";", 3 until 2nd ";", etc.

*/

/************************************************************************************/
/* output bonds in ascending order of the neighboring atom original numbers */
int WriteOrigBonds( int num_inp_atoms, inp_ATOM *at, int *i, char *szBuf, int buf_len, STRUCT_DATA *sd)
{
    int j, k, k2, kk, len, cur_len, j2=0, bond_stereo, bond_char, bond_parity, bond_parityNM, num_trans;
    char szCurBonds[7*MAXVAL+2]; /* num_neigh*(1 byte bond type + 2 bytes for bond parity up to 4 digits per neighbor number) + at the end one ';' */
    AT_RANK nNeighOrder[MAXVAL];
    int  chain_len, pnxt_atom, pinxt2cur, pinxt_sb_parity_ord;
    int  chain_len2, pnxt_atom2, pinxt2cur2, pinxt_sb_parity_ord2, m1, m2;
    int  pcur_atom, picur2nxt, picur_sb_parity_ord;

    cur_len = 0;
    for ( j = *i; j < num_inp_atoms; ) {
        len = 0;
        if ( at[j].valence > 1 ) {
            for ( k = 0; k < at[j].valence; k ++ ) {
                nNeighOrder[k] = k;
            }
            pn_RankForSort = at[j].neighbor;
            num_trans = insertions_sort( nNeighOrder, at[j].valence, sizeof(nNeighOrder[0]), CompRank );
        } else {
            num_trans = 0;
            nNeighOrder[0] = 0;
        }
        for ( kk = 0; kk < at[j].valence; kk ++ ) {
            k = nNeighOrder[kk];
            j2 = at[j].neighbor[k];
            bond_parity = 0;
            bond_parityNM = 0;
            if ( j2 < j ) {
                bond_stereo = at[j].bond_stereo[k];
                switch( at[j].bond_type[k] ) {
                case BOND_TYPE_SINGLE:
                    switch( bond_stereo ) {
                    case  STEREO_SNGL_UP:
                        bond_char = 'p';
                        break;
                    case -STEREO_SNGL_UP:
                        bond_char = 'P';
                        break;
                    case  STEREO_SNGL_DOWN:
                        bond_char = 'n';
                        break;
                    case -STEREO_SNGL_DOWN:
                        bond_char = 'N';
                        break;
#if ( FIX_EITHER_STEREO_IN_AUX_INFO == 1 )
                    case  STEREO_SNGL_EITHER:
                        bond_char = 'v';
                        break;
                    case -STEREO_SNGL_EITHER:
                        bond_char = 'V';
                        break;
#else
                    case  STEREO_SNGL_EITHER:
                    case -STEREO_SNGL_EITHER:
                        bond_char = 'v';
                        break;
#endif
                    default:
                        bond_char = 's';
                        break;
                    }
                    break;
                case BOND_TYPE_DOUBLE:
                    switch( bond_stereo ) {
                    case  STEREO_DBLE_EITHER:
                    case -STEREO_DBLE_EITHER:
                        bond_char = 'w';
                        break;
                    default:
                        bond_char = 'd';
                        break;
                    }
                    break;
                case BOND_TYPE_TRIPLE:
                    bond_char = 't';
                    break;
                case BOND_TYPE_ALTERN:
                    bond_char = 'a';
                    break;
                default:
                    bond_char = 's';
                    break;
                }
                /* check for allene/cumulene */
                k2 = is_in_the_list( at[j2].neighbor, (AT_NUMB)j, at[j2].valence ) - at[j2].neighbor;
                chain_len = chain_len2 = 0;
                if ( at[j].sb_parity[0] ) {
                    for ( m1 = 0; m1 < MAX_NUM_STEREO_BONDS && at[j].sb_parity[m1]; m1 ++ ) {
                        if ( k == at[j].sb_ord[m1] ) {
                            chain_len = get_opposite_sb_atom( at, j, k,
                                          &pnxt_atom, &pinxt2cur, &pinxt_sb_parity_ord );
                            break;
                        }
                    }
                }
                if ( at[j2].sb_parity[0] ) {
                    for ( m2 = 0; m2 < MAX_NUM_STEREO_BONDS && at[j2].sb_parity[m2]; m2 ++ ) {
                        if ( k2 == at[j2].sb_ord[m2] ) {
                            chain_len2 = get_opposite_sb_atom( at, j2, k2,
                                           &pnxt_atom2, &pinxt2cur2, &pinxt_sb_parity_ord2 );
                            break;
                        }
                    }
                }
                if ( (chain_len == 1 && chain_len2 == 1) ||  /* regular stereobond */
                     (chain_len  > 1 && j  > pnxt_atom) ) {  /* j  is a cumulene endpoint */
                    int m;
                    pcur_atom = j;  /* pcur_atom > pnxt_atom */
                    picur2nxt = k;
                    picur_sb_parity_ord = -1;
                    for ( m = 0; m < MAX_NUM_STEREO_BONDS && at[pcur_atom].sb_parity[m]; m ++ ) {
                        if ( at[pcur_atom].sb_ord[m] == k ) {
                            picur_sb_parity_ord = m;
                            break;
                        }
                    }
                    chain_len2 = 0;
                } else
                if ( chain_len2 > 1 && j2 > pnxt_atom2  ) { /* j2 is a cumulene endpoint */
                    int m;
                    pcur_atom = j2;
                    picur2nxt = k2;
                    pnxt_atom = pnxt_atom2;
                    pinxt2cur = pinxt2cur2;
                    pinxt_sb_parity_ord = pinxt_sb_parity_ord2;
                    picur_sb_parity_ord = -1;
                    for ( m = 0; m < MAX_NUM_STEREO_BONDS && at[pcur_atom].sb_parity[m]; m ++ ) {
                        if ( at[pcur_atom].sb_ord[m] == k2 )
                            picur_sb_parity_ord = m;
                    }
                    chain_len  = chain_len2;
                    chain_len2 = 0;
                } else {
                    chain_len = chain_len2 = 0;
                }

                /*len += sprintf( szCurBonds + len, "%c%d", bond_char, val+1);*/
                if ( chain_len ) {
                    /* both atoms belong to a stereo bond */
                    int kc;
                    int p1, p2, p1NM, p2NM, neigh, neigh1, neigh2, bHasMetal, bWellDef;
                    int     bNeighSwitched1, bNeighSwitched2;

                    p1   = SB_PARITY_1(at[pcur_atom].sb_parity[picur_sb_parity_ord]);
                    p1NM = SB_PARITY_2(at[pcur_atom].sb_parity[picur_sb_parity_ord]);
                    p2   = SB_PARITY_1(at[pnxt_atom].sb_parity[pinxt_sb_parity_ord]);
                    p2NM = SB_PARITY_2(at[pnxt_atom].sb_parity[pinxt_sb_parity_ord]);

                    bWellDef  = ATOM_PARITY_WELL_DEF(p1)   && ATOM_PARITY_WELL_DEF(p2);
                    bHasMetal = ATOM_PARITY_WELL_DEF(p1NM) && ATOM_PARITY_WELL_DEF(p2NM);

                    bNeighSwitched1 = bNeighSwitched2 = 0;

                    if ( bWellDef || bHasMetal ) {

                        neigh1  = num_inp_atoms;
                        for ( kc = 0; kc < at[pcur_atom].valence; kc ++ ) {
                            if ( kc == picur2nxt )
                                continue;
                            neigh = at[pcur_atom].neighbor[kc];
                            if ( bHasMetal && is_el_a_metal( at[neigh].el_number ) )
                                continue;
                            if ( neigh < neigh1 )
                                neigh1 = neigh;
                        }
                        if ( neigh1 < num_inp_atoms ) {
                             bNeighSwitched1 = (neigh1 != at[pcur_atom].neighbor[(int)at[pcur_atom].sn_ord[picur_sb_parity_ord]]);
                        } else {
                            AddMOLfileError(sd->pStrErrStruct, "Cannot find 0D stereobond neighbor");
                            /*
                            sd->nStructReadError =  99;
                            sd->nErrorType = _IS_ERROR;
                            */

                        }

                        neigh2  = num_inp_atoms;
                        for ( kc = 0; kc < at[pnxt_atom].valence; kc ++ ) {
                            if ( kc == pinxt2cur )
                                continue;
                            neigh = at[pnxt_atom].neighbor[kc];
                            if ( bHasMetal && is_el_a_metal( at[neigh].el_number ) )
                                continue;
                            if ( neigh < neigh2 )
                                neigh2 = neigh;
                        }
                        if ( neigh2 < num_inp_atoms ) {
                             bNeighSwitched2 = (neigh2 != at[pnxt_atom].neighbor[(int)at[pnxt_atom].sn_ord[pinxt_sb_parity_ord]]);
                        } else {
                            AddMOLfileError(sd->pStrErrStruct, "Cannot find 0D stereobond neighbor");
                            /*
                            sd->nStructReadError =  99;
                            sd->nErrorType = _IS_ERROR;
                            */

                        }

                        if ( neigh1 < num_inp_atoms && neigh2 < num_inp_atoms ) {
                            if ( ATOM_PARITY_WELL_DEF(p1) && ATOM_PARITY_WELL_DEF(p2) ) {
                                bond_parity = 2 - (p1 + p2 + bNeighSwitched1 + bNeighSwitched2) % 2;
                            } else {
                                bond_parity = inchi_min( p1, p2 );
                            }

                            if ( bHasMetal ) {
                                bond_parityNM = 2 - (p1NM + p2NM + bNeighSwitched1 + bNeighSwitched2) % 2;
                            } else
                            if ( p1NM && p2NM ) {
                                bond_parityNM = inchi_min( p1NM, p2NM );
                            }
                        }
                    } else {
                        if ( p1 && p2 ) {
                            bond_parity = inchi_min( p1, p2 );
                        }
                        if ( p1NM && p2NM ) {
                            bond_parityNM = inchi_min( p1NM, p2NM );
                        }
                        if ( bond_parityNM && !bond_parity ) {
                            bond_parity = AB_PARITY_UNDF;
                        }
                    }
                }
                len += sprintf( szCurBonds + len, "%c%s%s%d",

                                                  bond_char,

                                                  (bond_parity == AB_PARITY_ODD)?  "-" :
                                                  (bond_parity == AB_PARITY_EVEN)? "+" :
                                                  (bond_parity == AB_PARITY_UNKN)? "u" :
                                                  (bond_parity == AB_PARITY_UNDF)? "?" : "",

                                                  (bond_parityNM == AB_PARITY_ODD)?  "-" :
                                                  (bond_parityNM == AB_PARITY_EVEN)? "+" :
                                                  (bond_parityNM == AB_PARITY_UNKN)? "u" :
                                                  (bond_parityNM == AB_PARITY_UNDF)? "?" : "",

                                                  j2+1);
            }
        }
        if ( len + cur_len + 2 < buf_len ) {
            memcpy( szBuf + cur_len, szCurBonds, len );
            cur_len += len;
            szBuf[ cur_len ++ ] = ';';
            j ++;
        } else {
            break;
        }
    }
    szBuf[cur_len] = '\0';
    *i = num_inp_atoms>0? j : 0;
    return cur_len;
}


#define ORIG_STR_BUFLEN (7*MAXVAL+2)  /* > 7*MAXVAL+2 = 142 */
/******************************************************************************************/
int FillOutOrigStruct( ORIG_ATOM_DATA *orig_inp_data, ORIG_STRUCT *pOrigStruct, STRUCT_DATA *sd )
{
    char szBuf[ORIG_STR_BUFLEN];
    int  i, len, len_coord, len_atoms, len_bonds;
    /* coordinates */
    len_coord = i = 0;

    if (orig_inp_data->szCoord) {

        while ( (len = WriteOrigCoord( orig_inp_data->num_inp_atoms,
                                      orig_inp_data->szCoord, &i, szBuf, sizeof(szBuf) )) ) {
            len_coord += len;
        }
        pOrigStruct->szCoord = (char*) inchi_malloc( (len_coord + 1)*sizeof(pOrigStruct->szCoord[0]) );
        i = 0;
        if ( pOrigStruct->szCoord &&
             len_coord == WriteOrigCoord( orig_inp_data->num_inp_atoms,
                                      orig_inp_data->szCoord, &i, pOrigStruct->szCoord, len_coord+1 ) &&
             i == orig_inp_data->num_inp_atoms ) {
            /* success */
            if ( orig_inp_data->szCoord ) {
                inchi_free( orig_inp_data->szCoord );
                orig_inp_data->szCoord = NULL;
            }
        } else {
            return -1;
        }

    }

    /* atoms */
    len_atoms = i = 0;
    while ( (len = WriteOrigAtoms( orig_inp_data->num_inp_atoms,
                                  orig_inp_data->at, &i, szBuf, sizeof(szBuf), sd)) ) {
        len_atoms += len;
        if ( !orig_inp_data->num_inp_atoms )
            break;
    }
    pOrigStruct->szAtoms = (char*) inchi_malloc( (len_atoms + 1)*sizeof(pOrigStruct->szAtoms[0]) );
    i = 0;
    if ( pOrigStruct->szAtoms &&
         len_atoms == WriteOrigAtoms( orig_inp_data->num_inp_atoms,
                                  orig_inp_data->at, &i, pOrigStruct->szAtoms, len_atoms+1, sd ) &&
         i == orig_inp_data->num_inp_atoms ) {
        ; /* success */
    } else {
        return -1;
    }
    /* bonds */
    len_bonds = 0;
    i = 1;
    while ( (len = WriteOrigBonds( orig_inp_data->num_inp_atoms,
                                  orig_inp_data->at, &i, szBuf, sizeof(szBuf), NULL)) ) {
        len_bonds += len;
        if ( !orig_inp_data->num_inp_atoms )
            break;
    }
    pOrigStruct->szBonds = (char*) inchi_malloc( (len_bonds + 2)*sizeof(pOrigStruct->szBonds[0]) );
    i = 1;
    if ( pOrigStruct->szBonds &&
         len_bonds == WriteOrigBonds( orig_inp_data->num_inp_atoms,
                                  orig_inp_data->at, &i, pOrigStruct->szBonds, len_bonds+2, sd ) &&
         i == orig_inp_data->num_inp_atoms ) {
        ; /* success */
    } else {
        return -1;
    }
    pOrigStruct->num_atoms = orig_inp_data->num_inp_atoms;
    return 0;
}
/*****************************************************************/
void FreeOrigStruct(  ORIG_STRUCT *pOrigStruct)
{
    if ( pOrigStruct ) {
        if ( pOrigStruct->szAtoms )
            inchi_free( pOrigStruct->szAtoms );
        if ( pOrigStruct->szBonds )
            inchi_free( pOrigStruct->szBonds );
        if ( pOrigStruct->szCoord )
            inchi_free( pOrigStruct->szCoord );
        memset( pOrigStruct, 0, sizeof(*pOrigStruct) );

    }
}



/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Get the two letters encoding the saved InChI creation options.

The first one encodes RecMet/FixedH/SUU/SLUUD options.
Each of options is a binary switch {ON,OFF}, so it totals to 2*2*2*2=16 values
which are encoded by capital letters �A� through �P�.

The second character encodes experimental (InChI 1 extension) options KET and 15T.
Each of these options is a binary switch ON/OFF, so there are 2*2=4 combinations,
currently encoded by �A� through �D�.
Note that anything but 'A' here would indicate "extended" InChI 1 Also, there is a
reservation for future needs: the 2nd memo char may accommodate two more ON/OFF
binary options (at 26-base encoding).
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
void GetSaveOptLetters(unsigned char save_opt_bits, char* let1, char* let2)
{
const char a2p[]="ABCDEFGHIJKLMNOP";
    /* SaveOptBits layout: {unused|unused|Ket|15T|RecMet|FixedH|SUU|SLUUD} */
    *let1 = a2p [ (size_t) ( save_opt_bits & 0x0f ) ];
    *let2 = a2p [ (size_t) ( (save_opt_bits & 0x30) >> 4 ) ];
}

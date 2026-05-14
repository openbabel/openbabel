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

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "mode.h"

#include "ichister.h"
#include "ichimain.h"
#include "ichimake.h"
#include "ichi_io.h"

#include "bcf_s.h"

#include "logging.h"                        /*(@nnuk : Nauman Ullah Khan) :: Needed for logging functionality*/

/*
    Local functions
*/

/* djb-rwth: removing redundant code */
static int str_LineEnd( const char *tag,
                        int *bOverflow,
                        INCHI_IOS_STRING *buf,
                        int ind,
                        int bPlainTextTags );
static int CleanOrigCoord( MOL_COORD szCoord, int delim );
static int WriteOrigCoord( int num_inp_atoms,
                           MOL_COORD *szMolCoord,
                           int *i,
                           char *szBuf,
                           int buf_len );
static int WriteOrigAtoms( CANON_GLOBALS *pCG,
                           int num_inp_atoms,
                           inp_ATOM *at,
                           int *i,
                           char *szBuf,
                           int buf_len,
                           STRUCT_DATA *sd );
static int WriteOrigBonds( CANON_GLOBALS *pCG,
                           int num_inp_atoms,
                           inp_ATOM *at,
                           int *i,
                           char *szBuf,
                           int buf_len,
                           STRUCT_DATA *sd );
static void GetSaveOptLetters( unsigned char save_opt_bits,
                               char* let1,
                               char* let2 );
static int OutputINCHI_VersionAndKind( INCHI_IOSTREAM *out_file,
                                       INCHI_IOS_STRING *strbuf,
                                       int bINChIOutputOptions,
                                       int is_beta,
                                       char *pLF,
                                       char *pTAB );
static int OutputINCHI_MainLayerFormula( CANON_GLOBALS *pCG,
                                         INCHI_IOSTREAM *out_file,
                                         INCHI_IOS_STRING *strbuf,
                                         int num_components2[],
                                         int *INCHI_basic_or_INCHI_reconnected,
                                         INCHI_OUT_CTL *io,
                                         char *pLF,
                                         char *pTAB );
static int OutputINCHI_MainLayerConnections( CANON_GLOBALS *pCG,
                                             INCHI_IOSTREAM *out_file,
                                             INCHI_IOS_STRING *strbuf,
                                             int num_components2[],
                                             int *INCHI_basic_or_INCHI_reconnected,
                                             INCHI_OUT_CTL *io,
                                             char *pLF,
                                             char *pTAB );
static int OutputINCHI_MainLayerHydrogens( CANON_GLOBALS *pCG,
                                           INCHI_IOSTREAM *out_file,
                                           INCHI_IOS_STRING *strbuf,
                                           int num_components2[],
                                           int *INCHI_basic_or_INCHI_reconnected,
                                           INCHI_OUT_CTL *io,
                                           char *pLF,
                                           char *pTAB );
static int OutputINCHI_ChargeAndRemovedAddedProtonsLayers( CANON_GLOBALS *pCG,
                                                           INCHI_IOSTREAM *out_file,
                                                           INCHI_IOS_STRING *strbuf,
                                                           INCHI_OUT_CTL *io,
                                                           char *pLF,
                                                           char *pTAB );
static int OutputINCHI_StereoLayer( CANON_GLOBALS *pCG,
                                    INCHI_IOSTREAM *out_file,
                                    INCHI_IOS_STRING *strbuf,
                                    INCHI_OUT_CTL *io,
                                    char *pLF,
                                    char *pTAB );
static int OutputINCHI_IsotopicLayer( CANON_GLOBALS *pCG,
                                      INCHI_IOSTREAM *out_file,
                                      INCHI_IOS_STRING *strbuf,
                                      int *INCHI_basic_or_INCHI_reconnected,
                                      INCHI_OUT_CTL *io,
                                      char *pLF,
                                      char *pTAB );
static int OutputINCHI_FixedHLayerWithSublayers( CANON_GLOBALS *pCG,
                                                 INCHI_IOSTREAM *out_file,
                                                 INCHI_IOS_STRING *strbuf,
                                                 int *INCHI_basic_or_INCHI_reconnected,
                                                 INCHI_OUT_CTL *io,
                                                 char *pLF,
                                                 char *pTAB,
                                                 int *then_goto_repeat );
static int OutputINCHI_PolymerLayer( CANON_GLOBALS *pCG, INCHI_IOSTREAM *out_file, INCHI_IOS_STRING *strbuf,
                                     int *INCHI_basic_or_INCHI_reconnected,
                                     ORIG_ATOM_DATA *orig_inp_data, ORIG_STRUCT *pOrigStruct,
                                     INCHI_OUT_CTL *io, char *pLF, char *pTAB );
static int OutputINCHI_PolymerLayer_SingleUnit( OAD_PolymerUnit *u,
                                                int bPolymers,
                                                int total_star_atoms, 
                                                int *n_used_stars,
                                                OAD_AtProps *aprops, 
                                                int *cano_nums,
                                                ORIG_ATOM_DATA *orig_inp_data, 
                                                ORIG_STRUCT *pOrigStruct,
                                                INCHI_IOS_STRING *strbuf );
static int OutputAUXINFO_HeaderAndNormalization_type( CANON_GLOBALS *pCG,
                                                      INCHI_IOSTREAM *out_file,
                                                      INCHI_IOS_STRING *strbuf,
                                                      int bINChIOutputOptions,
                                                      int *INCHI_basic_or_INCHI_reconnected,
                                                      int num_components2[],
                                                      INCHI_OUT_CTL *io,
                                                      char *pLF,
                                                      char *pTAB );
static int OutputAUXINFO_OriginalNumbersAndEquivalenceClasses( CANON_GLOBALS *pCG,
                                                               INCHI_IOSTREAM *out_file,
                                                               INCHI_IOS_STRING *strbuf,
                                                               int num_components2[],
                                                               INCHI_OUT_CTL *io,
                                                               char *pLF,
                                                               char *pTAB );
static int OutputAUXINFO_TautomericGroupsEquivalence( CANON_GLOBALS         *pCG,
                                                      INCHI_IOSTREAM        *out_file,
                                                      INCHI_IOS_STRING *strbuf,
                                                      INCHI_OUT_CTL   *io );
static int OutputAUXINFO_Stereo( CANON_GLOBALS *pCG,
                                  INCHI_IOSTREAM *out_file,
                                  INCHI_IOS_STRING *strbuf,
                                  INCHI_OUT_CTL *io,
                                  char *pLF,
                                  char *pTAB );
static int OutputAUXINFO_IsotopicInfo( CANON_GLOBALS *pCG,
                                       INCHI_IOSTREAM *out_file,
                                       INCHI_IOS_STRING *strbuf,
                                       int *INCHI_basic_or_INCHI_reconnected,
                                       INCHI_OUT_CTL *io,
                                       char *pLF, char *pTAB );
/* djb-rwth: removing redundant code */
static int OutputAUXINFO_ChargesRadicalsAndUnusualValences( CANON_GLOBALS *pCG,
                                                            INCHI_IOSTREAM *out_file,
                                                            INCHI_IOS_STRING *strbuf,
                                                            INCHI_OUT_CTL *io,
                                                            char *pLF,
                                                            char *pTAB );
static int OutputAUXINFO_ReversibilityInfo( CANON_GLOBALS *pCG,
                                            INCHI_IOSTREAM *out_file,
                                            INCHI_IOS_STRING *strbuf,
                                            ORIG_STRUCT *pOrigStruct,
                                            INCHI_OUT_CTL *io,
                                            char *pLF,
                                            char *pTAB );

static int OutputAUXINFO_PolymerInfo( CANON_GLOBALS *pCG,
                                      INCHI_IOSTREAM *out_file,
                                      INCHI_IOS_STRING *strbuf,
                                      ORIG_STRUCT *pOrigStruct,
                                      INCHI_OUT_CTL *io,
                                      char *pLF,
                                      char *pTAB );

static int InternallyGetCanoNumsAndComponentNums( CANON_GLOBALS         *pCG,
                                                  INCHI_IOS_STRING *strbuf,
                                                  INCHI_OUT_CTL   *io,
                                                  int                   nat,
                                                  int                   *cano_nums,
                                                  int                   *compnt_nums );

static int  CountPseudoElementInFormula( const char *pseudo, char *s );
static int  IsBondAtomNumsLesser( int *bond1, int* bond2 );

static void inchi_sort_int_pair_ascending( int* a, int* b );

/* djb-rwth: removing redundant code */

static void MergeZzInStrHillFormulaComponent( char *s );

/*
    Local constants
*/
const char sCompDelim[] = ";"; /* component delimiter */
const char sIdenticalValues[] = "*"; /* identical component */
const char x_space[] = "                  ";


/*
    Output: words & additional tags
*/
const char x_inchi[] = INCHI_NAME;
const char x_inchi_ver[] = "version"; /* "InChI.version"; */
const char x_curr_ver[] = INCHI_VERSION;
const char x_structure[] = "structure";
const char x_number[] = "number";
const char x_header[] = "id.name";
const char x_value[] = "id.value";
const char x_empty[] = "";
const char x_type[] = "type";
const char x_message[] = "message";
const char x_text[] = "value";
const char x_ferr[] = "fatal (aborted)";
const char x_err[] = "error (no InChI)";
const char x_warn[] = "warning";
const char x_basic[] = "identifier";
const char x_tautomeric[] = "mobile-H";
const char x_reconnected[] = "reconnected";
const char x_ver[] = "version";
const char x_type_alpha[] = "alpha";
const char x_type_numer[] = "numeric";
const char x_type_predec[] = "sct";
const char x_type_normal[] = "normal";
const char x_type_short[] = "compressed";
const char x_basic_layer[] = "basic";
const char x_aux_basic[] = "identifier.auxiliary-info";
const char x_aux_comm[] = "!-- This section is NOT a part of the identifier, it is not unique --";
const char x_ign_uu_sp2[] = "omit_undef_dbond";
const char x_ign_uu_sp3[] = "omit_undef_sp3";
const char x_line_opening[] = "<";
const char x_line_closing[] = "</";
const char x_close_line[] = ">";
const char x_abs[] = "1";
const char x_rel[] = "2";
const char x_rac[] = "3";

typedef struct tagInchiTag
{
    const char *szPlainLabel;
    const char *szPlainComment;
    const char *szXmlLabel;
    int  bAlwaysOutput;
} INCHI_TAG;


/*
    Identifier
*/
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



/*
    Aux Info constants
*/
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

/* const int MAX_TAG_NUM = inchi_max((short)IL_MAX_ORD, (short)AL_MAX_ORD); */ /* djb-rwth: fixing MSVC warning C5287 */

char *szGetTag( const INCHI_TAG *Tag, int nTag, int bTag, char *szTag, int *bAlways, short tag_flag ); /* djb-rwth: fixing GHI #160 */

#define SP(N)        (x_space+sizeof(x_space)-1-(N))



#define NOT_YET_I2I_FOR_POLYMERS 40


/****************************************************************************
  Print error message (plain text)
****************************************************************************/
int OutputINChIPlainError( INCHI_IOSTREAM *out_file,
                           char           *pErrorText,
                           int            bError )
{
    /* char szBuf[64]; */
    const char *pErr;
    char *szErrorText = pErrorText;
    int ret = 0; /* djb-rwth: removing redundant variables */

    switch (bError)
    {
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

    /* djb-rwth: removing redundant code */

    inchi_ios_print( out_file,
                     "%s: %s=\"%s\" %s=\"%s\"",
                     x_message, x_type, pErr, x_text, szErrorText );
#ifdef TARGET_LIB_FOR_WINCHI
    inchi_ios_print( out_file, "\n" );
#endif
    ret = 1;

    return ret;
}


#ifndef OUT_TN    /* defined in mode.h; quoted here for reference purposes only */

#define OUT_N1              0    /* non-tautomeric only */
#define OUT_T1              1    /* tautomeric if present otherwise non-tautomeric */
#define OUT_NT              2    /* only non-taut representations of tautomeric */
#define OUT_TN              3    /* tautomeric if present otherwise non-tautomeric;
                                    sepatately output non-taut representations of tautomeric if present */
/* OUT_TN = OUT_T1 + OUT_NT */
#endif


/****************************************************************************
 Calculate equivalence mark (used to check for repeating (sub)layer(s) )
****************************************************************************/
const char *EquString( int EquVal )
{
    int bFrom = EquVal & ( iiSTEREO | iiSTEREO_INV | iiNUMB | iiEQU );
    int bType = EquVal & ( iitISO | iitNONTAUT );
    int bEq2 = EquVal & ( iiEq2NONTAUT | iiEq2ISO | iiEq2INV );
    const char *r = "";

#if ( FIX_EMPTY_LAYER_BUG == 1 )
    int bEmpty = EquVal & iiEmpty;
    if (bEmpty)
    {
        r = "e";
        return r;
    }
#endif

    switch (bFrom)
    {

        case iiSTEREO:  /* ------------ Stereo --------------------*/
            switch (bType)
            {
                case iitISO:  /* iso main stereo =... */
                    switch (bEq2)
                    {
                        case 0:
                            r = "m";            /* iso main stereo = main stereo */
                            break;
                        default:
                            r = "??";           /* should not happen */
                            break;
                    }
                    break;
                case iitNONTAUT: /* non-taut stereo =... */
                    switch (bEq2)
                    {
                        case 0:
                            r = "m";            /* non-taut stereo = main stereo */
                            break;
                        default:
                            r = "??";           /* should not happen */
                            break;
                    }
                    break;
                case ( iitNONTAUT | iitISO ): /* iso non-taut stereo = ... */
                    switch (bEq2)
                    {
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
            if (bEq2 & iiEq2INV)
            { /* stereo = Inverted(another stereo) */
                bEq2 &= ~iiEq2INV;
                switch (bType)
                {
                    case 0: /* main = ...*/
                        switch (bEq2)
                        {
                            case 0:
                                r = "im";       /* main         = Inv(main) */
                                break;
                            case iiEq2ISO:
                                r = "iM";       /* main         = Inv(main iso) */
                                break;
                            case iiEq2NONTAUT:
                                r = "in";       /* maim         = Inv(non-taut) */
                                break;
                            case ( iiEq2NONTAUT | iiEq2ISO ):
                                r = "iN";       /* maim         = Inv(non-taut iso ) */
                                break;
                            default:
                                r = "??";           /* should not happen */
                                break;
                        }
                        break;
                    case iitISO: /* main iso = ...*/
                        switch (bEq2)
                        {
                            case 0:
                                r = "im";       /* main iso     = Inv(main) */
                                break;
                            case iiEq2ISO:
                                r = "iM";       /* main iso     = Inv(main iso) */
                                break;
                            case iiEq2NONTAUT:
                                r = "in";       /* maim iso     = Inv(non-taut) */
                                break;
                            case ( iiEq2NONTAUT | iiEq2ISO ):
                                r = "iN";       /* maim         = Inv(non-taut iso ) */
                                break;
                            default:
                                r = "??";           /* should not happen */
                                break;
                        }
                        break;
                    case iitNONTAUT: /* non-taut = ... */
                        switch (bEq2)
                        {
                            case 0:
                                r = "im";       /* non-taut     = Inv(main) */
                                break;
                            case iiEq2ISO:
                                r = "iM";       /* non-taut     = Inv(main iso) */
                                break;
                            case iiEq2NONTAUT:
                                r = "in";       /* non-taut     = Inv(non-taut) */
                                break;
                            case ( iiEq2NONTAUT | iiEq2ISO ):
                                r = "iN";       /* non-taut     = Inv(non-taut iso ) */
                                break;
                            default:
                                r = "??";           /* should not happen */
                                break;
                        }
                        break;
                    case ( iitNONTAUT | iitISO ):
                        switch (bEq2)
                        {
                            case 0:
                                r = "im";       /* non-taut iso = Inv(main) */
                                break;
                            case iiEq2ISO:
                                r = "iM";       /* non-taut iso = Inv(main iso) */
                                break;
                            case iiEq2NONTAUT:
                                r = "in";       /* non-taut iso = Inv(non-taut) */
                                break;
                            case ( iiEq2NONTAUT | iiEq2ISO ):
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
            }
            else
            {  /* Inv stereo = another (non-inverted) stereo */

                switch (bType)
                {
                    case iitISO: /* main iso = ...*/
                        switch (bEq2)
                        {
                            case 0:
                                r = "m";       /* main         = (inverted aux) main */
                                break;
                            default:
                                r = "??";           /* should not happen */
                                break;
                        }
                        break;
                    case iitNONTAUT: /* non-taut = ... */
                        switch (bEq2)
                        {
                            case 0:
                                r = "m";       /* non-taut     = (inverted aux) main */
                                break;
                            default:
                                r = "??";           /* should not happen */
                                break;
                        }
                        break;
                    case ( iitNONTAUT | iitISO ): /* non-taut iso = ...*/
                        switch (bEq2)
                        {
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

        case ( iiNUMB | iiSTEREO_INV ): /*------------- Inv Stereo Numbering ------------*/
            switch (bType)
            {
                case 0: /* inv stereo numb main = ...*/
                    switch (bEq2)
                    {
                        case 0:
                            r = "m";       /* inv stereo numb main     = main numb */
                            break;
                        default:
                            r = "??";      /* should not happen */
                            break;
                    }
                    break;
                case iitISO: /* inv stereo iso numb main = ...*/
                    switch (bEq2)
                    {
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
                    switch (bEq2)
                    {
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
                case ( iitNONTAUT | iitISO ): /* inv stereo numb non-taut iso = ... */
                    switch (bEq2)
                    {
                        case 0:
                            r = "m";       /* inv stereo numb non-taut iso = main numb */
                            break;
                        case iiEq2ISO:
                            r = "M";       /* inv stereo numb non-taut iso = main numb iso */
                            break;
                        case ( iiEq2ISO | iiEq2INV ):
                            r = "iM";       /* inv stereo numb non-taut iso = InvStereo(main iso) numb */
                            break;
                        case iiEq2NONTAUT:
                            r = "n";       /* inv stereo numb non-taut iso = non-taut numb */
                            break;
                        case ( iiEq2NONTAUT | iiEq2ISO ):
                            r = "N";       /* inv stereo numb non-taut iso = non-taut iso numb */
                            break;
                        case iiEq2INV:
                            r = "im";      /* inv stereo numb non-taut iso = InvStereo(main) numb */
                            break;
                        case ( iiEq2NONTAUT | iiEq2INV ):
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
            switch (bType)
            {
                case 0:         /* numb main = ...*/
                    r = "??";      /* should not happen */
                    break;
                case iitISO:     /* iso numb main = ...*/
                    switch (bEq2)
                    {
                        case 0:
                            r = "m";       /* iso numb main = main numb  */
                            break;
                        default:
                            r = "??";      /* should not happen */
                    }
                    break;
                case iitNONTAUT: /* numb non-taut = ... */
                    switch (bEq2)
                    {
                        case 0:
                            r = "m";       /* numb non-taut = main numb */
                            break;
                        default:
                            r = "??";      /* should not happen */
                    }
                    break;
                case ( iitNONTAUT | iitISO ): /* numb non-taut iso = ... */
                    switch (bEq2)
                    {
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
            switch (bType)
            {
                case 0:         /* equivalence main = ...*/
                    r = "??";      /* should not happen */
                    break;
                case iitISO:     /* equivalence main iso = ...*/
                    switch (bEq2)
                    {
                        case 0:
                            r = "m";       /* equivalence main = main equ  */
                            break;
                        default:
                            r = "??";      /* should not happen */
                            break;
                    }
                    break;
                case iitNONTAUT: /* equivalence non-taut = ... */
                    switch (bEq2)
                    {
                        case 0:
                            r = "m";       /* equivalence non-taut = main equ */
                            break;
                        default:
                            r = "??";      /* should not happen */
                            break;
                    }
                    break;
                case ( iitNONTAUT | iitISO ): /*  equivalence non-taut iso = ... */
                    switch (bEq2)
                    {
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

#define OUT_NONTAUT  OUT_NN  /* was OUT_NT until 2004-04-07 */

/****************************************************************************
  OutputINChI2( ... ) is called from SortAndPrintINChI( ... )
****************************************************************************/
int OutputINChI2( CANON_GLOBALS     *pCG,
                  INCHI_IOS_STRING  *strbuf,
                  INCHI_SORT        *pINChISortTautAndNonTaut2[][TAUT_NUM],
                  int               INCHI_basic_or_INCHI_reconnected,
                  ORIG_ATOM_DATA    *orig_inp_data,
                  ORIG_STRUCT       *pOrigStruct,
                  INPUT_PARMS       *ip,
                  int               bDisconnectedCoord,
                  int               bOutputType,
                  int               bINChIOutputOptions,
                  int               num_components2[],
                  int               num_non_taut2[],
                  int               num_taut2[],
                  INCHI_IOSTREAM    *output_file,
                  INCHI_IOSTREAM    *log_file,
                  int               num_input_struct,
                  int               *pSortPrintINChIFlags,
                  unsigned char     save_opt_bits )
{
    int bINChIOutputOptions0 = bINChIOutputOptions & ~( INCHI_OUT_PLAIN_TEXT | INCHI_OUT_PLAIN_TEXT_COMMENTS );
    int bINChIOutputOptionsCur;
    int bCurOption, ret, i;

    ret = 0;

    for (i = 0; i < 3; i++)
    {
        switch (i)
        {
            case 1:
                bCurOption = INCHI_OUT_PLAIN_TEXT;
                break;
            case 2:
                bCurOption = INCHI_OUT_PLAIN_TEXT_COMMENTS;
                break;
            default:
                continue;
        }
        if (bINChIOutputOptions & bCurOption)
        {
            bINChIOutputOptionsCur = bINChIOutputOptions0 | bCurOption;
            if (i != 1)
            {
                bINChIOutputOptionsCur &= ~INCHI_OUT_TABBED_OUTPUT;
            }
            ret |= OutputINChI1( pCG,
                                 strbuf,
                                 pINChISortTautAndNonTaut2,
                                 INCHI_basic_or_INCHI_reconnected,
                                 orig_inp_data,
                                 pOrigStruct,
                                 ip,
                                 bDisconnectedCoord,
                                 bOutputType,
                                 bINChIOutputOptionsCur,
                                 num_components2,
                                 num_non_taut2,
                                 num_taut2,
                                 output_file,
                                 log_file,
                                 num_input_struct,
                                 pSortPrintINChIFlags,
                                 save_opt_bits );
        }
    }

    return ret;
}


/*                                                          */
/*    OutputINChI1( ... )                                   */
/*                                                          */
/*    Main actual worker which serializes InChI to string.  */
/*                                                          */
/*    Called from OutputINChI2( ... ) and from itself       */
/*                                                          */
int OutputINChI1( CANON_GLOBALS *pCG,
                  INCHI_IOS_STRING *strbuf,
                  INCHI_SORT *pINChISortTautAndNonTaut2[][TAUT_NUM],
                  int INCHI_basic_or_INCHI_reconnected,
                  ORIG_ATOM_DATA *orig_inp_data,
                  ORIG_STRUCT *pOrigStruct,
                  INPUT_PARMS *ip,
                  int bDisconnectedCoord,
                  int bOutputType,
                  int bINChIOutputOptions,
                  int num_components2[],
                  int num_non_taut2[],
                  int num_taut2[],
                  INCHI_IOSTREAM *out_file,
                  INCHI_IOSTREAM *log_file,
                  int num_input_struct,
                  int *pSortPrintINChIFlags,
                  unsigned char save_opt_bits )
{
    INCHI_OUT_CTL io;
    /*
        bINChIOutputOptions bits:
          INCHI_OUT_NO_AUX_INFO           0x0001    do not output Aux Info
          INCHI_OUT_SHORT_AUX_INFO        0x0002    output short version of Aux Info
          INCHI_OUT_ONLY_AUX_INFO         0x0004    output only Aux Info
          INCHI_OUT_EMBED_REC             0x0008    embed reconnected INChI into disconnected INChI

        bOutputType =
         TAUT_YES  => tautomeric only (if no tautomeric components then no output;
         TAUT_NON  => only non-tautomeric output (if no non-taut present then no output;
         TAUT_BOTH => tautomeric and non-tautomeric
    */
    int  i, j, ii, jj, /*ii2, jj2,*/ bEmbeddedOutputCalled = 0;
    int  bTautIsoHNum, bTautIsoAt, bHasIsotopicAtoms[TAUT_NUM];
    int  bStereoSp2[TAUT_NUM], bStereoSp3[TAUT_NUM];
    int  bIsotopicStereoSp2[TAUT_NUM], bIsotopicStereoSp3[TAUT_NUM];
    int  bStereoAbsInverted[TAUT_NUM], bIsotopicStereoAbsInverted[TAUT_NUM];
    int  bStereoAbs[TAUT_NUM], bIsotopicStereoAbs[TAUT_NUM];
    int  bTautomericAcid, bHardAddRemProton;
    int  bRequestedRacemicStereo = 0, bRequestedRelativeStereo = 0;
    int  npass = 0; /* djb-rwth: removing redundant variables */

    INCHI_SORT   *is, *is2;
    INChI        *pINChI /*, *pINChI2*/;
    INChI_Aux    *pINChI_Aux = NULL;

    int  ret = 0;        /*  0 failed, 1 success */
    int  intermediate_result = 0;
    int then_goto_repeat = 0;
    /* djb-rwth: removing redundant variables */
    int  bHasIsoH;
    /* djb-rwth: removing redundant variables */
    int  bTautAndNonTaut, bTautIsNonTaut;
    int nAtomsAllComp1, nAtomsAllComp2;    /* v. 1.05 Total atoms in all components */

    int  bPlainText = 0 != ( bINChIOutputOptions & ( INCHI_OUT_PLAIN_TEXT | INCHI_OUT_PLAIN_TEXT_COMMENTS ) );
    int  bPlainTextCommnts = 0 != ( bINChIOutputOptions & INCHI_OUT_PLAIN_TEXT_COMMENTS );

    char *pLF, *pTAB;
#ifdef TARGET_LIB_FOR_WINCHI
    int silent = 1;
#endif

    int bFixTranspChargeBug = 0;
#if ( FIX_TRANSPOSITION_CHARGE_BUG == 1 ) /* 2008-01-02 */
    if (INCHI_OUT_FIX_TRANSPOSITION_CHARGE_BUG & bINChIOutputOptions)
        bFixTranspChargeBug = 1;
#endif

    io.bAbcNumbers = ip->bAbcNumbers;

    io.ATOM_MODE = ( ( io.bAbcNumbers ? CT_MODE_ABC_NUMBERS : 0 )
                    | CT_MODE_ATOM_COUNTS
                    | CT_MODE_NO_ORPHANS
#if ( EQL_H_NUM_TOGETHER == 1 )
                    | CT_MODE_EQL_H_TOGETHER
#endif
#if ( ABC_CT_NUM_CLOSURES == 1 )
                    | ( io.bAbcNumbers && ip->bCtPredecessors ? CT_MODE_ABC_NUM_CLOSURES : 0 )
#endif
                    | ( ip->bCtPredecessors ? CT_MODE_PREDECESSORS : 0 ) );

    io.TAUT_MODE = ( io.bAbcNumbers ? CT_MODE_ABC_NUMBERS : 0 );
    io.pSortPrintINChIFlags = pSortPrintINChIFlags;
    io.num_components = num_components2[INCHI_basic_or_INCHI_reconnected];
    io.pINChISortTautAndNonTaut = pINChISortTautAndNonTaut2[INCHI_basic_or_INCHI_reconnected];
    io.pINChISort = io.pINChISortTautAndNonTaut[TAUT_YES];
    io.pINChISort2 = io.pINChISortTautAndNonTaut[TAUT_YES];
    io.bAlways = 0;
    io.bUseMulipliers = 1;
    io.bOmitRepetitions = 1;
    io.bPlainTextTags = 2;  /* 0 => no plain tags, 1=> plain text tags, 2=>plaintext tags without consecutive // */
    io.bOutputType = bOutputType; /* remains constant */
    io.bOutType = bOutputType;    /* will change! */
    io.bOverflow = 0;
    io.bSecondNonTautPass = 0;
    io.bNonTautIsoIdentifierNotEmpty = 0;
    io.bNonTautNonIsoIdentifierNotEmpty = 0;
    io.bNonTautIsIdenticalToTaut = 1;
    io.bFhTag = 0;
    io.nTag = bPlainTextCommnts ? 3 : bPlainText ? 2 : 0; /* tag type */


    if (NULL==orig_inp_data)
    {
        /*intermediate_result = 1;
        goto exit_function;*/
        io.n_zy     = 0;
        io.n_pzz    = 0;
        io.n_pzz    = 0;
    }
    else
    {
        io.n_zy = orig_inp_data->n_zy;
        io.n_pzz = 0;
        if (orig_inp_data->polymer)
        {
            io.n_pzz = orig_inp_data->polymer->n_pzz;
        }
    }
    

    io.bPolymers = ip->bPolymers;

    /* djb-rwth: removing redundant code */

#ifdef TARGET_LIB_FOR_WINCHI
    /* @@@ From now on we will go with silent output; it ends on @@@ below */
    silent = 1;
#endif

    /* Analyze layers, make adjustments and fixes, etc. */

    set_line_separators( bINChIOutputOptions, &pLF, &pTAB );
    memset( io.sDifSegs, DIFV_BOTH_EMPTY, sizeof( io.sDifSegs ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    if (!strbuf || !( strbuf->pStr ) || strbuf->nAllocatedLength <= 0)
    {
        inchi_ios_eprint( log_file, "Cannot allocate output buffer. No output for structure #%d.%s%s%s%s\n",
                         num_input_struct, SDF_LBL_VAL( ip->pSdfLabel, ip->pSdfValue ) );
        return ret;
    }

    /*    -- commented out to allow empty InChI --
    if (!io.num_components ) return 0;
    */

    for (i = 0; i < TAUT_NUM; i++)
    {
        bHasIsotopicAtoms[i] =
            io.num_comp[i] =
            bStereoSp2[i] =
            bStereoSp3[i] =
            bIsotopicStereoSp2[i] =
            bIsotopicStereoSp3[i] =
            io.bIsotopicOrigNumb[i] =
            bStereoAbs[i] =
            bIsotopicStereoAbs[i] =
            bStereoAbsInverted[i] =
            bIsotopicStereoAbsInverted[i] =
            io.bRacemicStereo[i] =
            io.bRelativeStereo[i] =
            io.bIsotopicRacemicStereo[i] =
            io.bIsotopicRelativeStereo[i] =
            io.bAtomEqu[i] =
            io.bTautEqu[i] =
            io.bIsotopicAtomEqu[i] =
            io.bIsotopicTautEqu[i] =
            io.bInvStereo[i] =
            io.bInvIsotopicStereo[i] =
            io.bInvStereoOrigNumb[i] =
            io.bInvIsotopicStereoOrigNumb[i] =
            io.bIgn_UU_Sp3[i] =
            io.bIgn_UU_Sp2[i] =
            io.bIgn_UU_Sp3_Iso[i] =
            io.bIgn_UU_Sp2_Iso[i] =
            io.bChargesRadVal[i] =
            io.bOrigCoord[i] = 0;
    }

    /*    Find if it is isotopic */
        io.bIsotopic =
        io.bTautomeric =
        io.bNonTautomeric =
        bTautomericAcid =
        bHardAddRemProton =
        bTautIsoHNum =
        bTautIsoAt =
        bTautAndNonTaut =
        bTautIsNonTaut = 0;

        /*
             x = bStereo, bStereoSp2, bStereoSp3, bStereoAbsInverted,
                 bIsotopicStereo, bIsotopicStereoSp2, bIsotopicStereoSp3, bIsotopicStereoAbsInverted

             OUT_N1: x[TAUT_NON] refers to non-tautomeric only
             OUT_T1: x[TAUT_YES] refers to tautomeric if exists otherwise non-tautomeric
             OUT_NT: x[TAUT_NON] refers to non-taut representations of tautomeric
             OUT_TN: x[TAUT_YES] refers to tautomeric if exists otherwise non-tautomeric
                     x[TAUT_NON] refers to non-taut representations of tautomeric
         */

    memset( io.num_iso_H, 0, sizeof( io.num_iso_H ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    io.nNumRemovedProtons = 0;
    /* djb-rwth: removing redundant code */
    bHasIsoH = 0;
    io.bTautomericOutputAllowed = ( io.bOutType == OUT_T1 || io.bOutType == OUT_TN );
    io.pINChISort = io.pINChISortTautAndNonTaut[io.bTautomericOutputAllowed ? TAUT_YES : TAUT_NON];
    is = io.pINChISort;
    /* djb-rwth: removing redundant variables/code */


    for (i = 0, is2 = io.pINChISortTautAndNonTaut[TAUT_NON]; i < io.num_components; i++, is++, is2 ? is2++ : NULL)
    {

        CompINChILayers( is, is2, io.sDifSegs, bFixTranspChargeBug );

        io.bNonTautIsIdenticalToTaut = io.bNonTautIsIdenticalToTaut && !CompINChITautVsNonTaut( is, is2, 1 );

        if (is && ( pINChI_Aux = is->pINChI_Aux[TAUT_YES] ))
        {
            for (j = 0; j < NUM_H_ISOTOPES; j++)
            {
                bHasIsoH += abs( pINChI_Aux->nNumRemovedIsotopicH[j] );
                io.num_iso_H[j] += pINChI_Aux->nNumRemovedIsotopicH[j];
            }
            io.nNumRemovedProtons += pINChI_Aux->nNumRemovedProtons;
            /* djb-rwth: removing redundant code */
        }

        if (io.bTautomericOutputAllowed && is) /* djb-rwth: fixing a NULL pointer dereference */
        {
            /* Check for removed isotopic H */
            for (j = TAUT_YES; j < TAUT_NUM; j++)
            {
                switch (io.bOutType)
                {
                    case OUT_N1: /* x[TAUT_NON]: non-tautomeric only -- never happens */
                        jj = GET_II( io.bOutType, is );
                        if (jj != j)
                            continue;
                        /* djb-rwth: removing redundant code */
                        break;
                    case OUT_T1: /* x[TAUT_YES]: tautomeric if present otherwise non-tautomeric */
                        jj = GET_II( io.bOutType, is );
                        if (jj != j)
                            continue;
                        /* djb-rwth: removing redundant code */
                        break;
                    case OUT_NT: /* x[TAUT_NON]: only non-taut representations of tautomeric -- never happens */
                        jj = GET_II( io.bOutType, is );
                        if (jj != j)
                            continue;
                        /* djb-rwth: removing redundant code */
                        break;
                    /* main path of control flow */
                    case OUT_TN: /* x[TAUT_YES]: tautomeric if present otherwise non-tautomeric;
                                  * x[TAUT_NON]: non-taut only if tautomeric is present */
                        jj = ( j == TAUT_YES ) ? GET_II( OUT_T1, is ) : ( j == TAUT_NON ) ? GET_II( OUT_NT, is ) : -1;
                        if (jj == TAUT_YES)
                        {
                            /* Fix12 */
                            if (is->pINChI[jj]->lenTautomer > 0)
                            {
                                bTautAndNonTaut += ( !is->pINChI[jj]->bDeleted && HAS_N( is ) );
                            }
                            else
                            {
                                bTautIsNonTaut++;
                            }
                        }
                        if (jj < 0)
                            continue;
                        /* djb-rwth: removing redundant code */
                        break;
                    default:
                        continue;
                }
                if (jj != j)
                    continue;
                if (( pINChI = is->pINChI[jj] ) && pINChI->nNumberOfAtoms > 0 && ( pINChI_Aux = is->pINChI_Aux[jj] ))
                {
                    bTautIsoHNum += ( pINChI_Aux->nNumRemovedIsotopicH[0] +
                                     pINChI_Aux->nNumRemovedIsotopicH[1] +
                                     pINChI_Aux->nNumRemovedIsotopicH[2] );
                    bTautIsoAt += ( pINChI->nNumberOfIsotopicAtoms > 0 || pINChI->nNumberOfIsotopicTGroups > 0 );
                }
            }
        }
    }

    io.sDifSegs[DIFL_M][DIFS_p_PROTONS] = io.nNumRemovedProtons ? DIFV_NEQ2PRECED : DIFV_BOTH_EMPTY;
    io.sDifSegs[DIFL_MI][DIFS_h_H_ATOMS] = bHasIsoH ? DIFV_NEQ2PRECED : DIFV_BOTH_EMPTY;

    MarkUnusedAndEmptyLayers( io.sDifSegs );

    io.bNonTautIsIdenticalToTaut = io.bNonTautIsIdenticalToTaut && !bTautIsoHNum;
    nAtomsAllComp1 = nAtomsAllComp2 = 0;

    for (i = 0, is = io.pINChISort; i < io.num_components; i++, is++)
    {
        int bCurIso, bCurHasIsoStereo /* Fix14 */, bCurTaut /*, bCurTaut2*/; /* djb-rwth: removing redundant variables */
        int bCompExists, bCurIsoHPos, bCurIsoHStereo;
        int bCurStereoSp2, bCurIsoStereoSp2, bCurStereoSp3, bCurIsoStereoSp3, bCurIsoStereoSp3Inv;
        int bCurRacemic, bCurRelative, bCurIsoRacemic, bCurIsoRelative;
        bCompExists = 0;

        for (j = TAUT_NON; j < TAUT_NUM; j++)
        {
            switch (io.bOutType)
            {
                case OUT_N1:
                    /* x[TAUT_NON]: non-tautomeric only */
                    jj = GET_II( io.bOutType, is );
                    if (jj != j)
                        continue;
                    ii = TAUT_NON;
                    break;
                case OUT_T1:
                    /* x[TAUT_YES]: tautomeric if present otherwise non-tautomeric */
                    jj = GET_II( io.bOutType, is );
                    if (jj != j)
                        continue;
                    ii = TAUT_YES;
                    break;
                case OUT_NT:
                    /* x[TAUT_NON]: only non-taut representations of tautomeric */
                    jj = GET_II( io.bOutType, is );
                    if (jj != j)
                        continue;
                    ii = TAUT_NON;
                    break;
                /* main control flow comes here: requested both mobile and fixed H results */
                case OUT_TN:
                    /* x[TAUT_YES]: tautomeric if present otherwise non-tautomeric; */
                    /* x[TAUT_NON]: non-taut only if tautomeric is present          */
                    jj = ( j == TAUT_YES ) ? GET_II( OUT_T1, is ) : ( j == TAUT_NON ) ? GET_II( OUT_NT, is ) : -1;
                    if (jj < 0)
                    {
                        /* Fix12 */
                        if (bTautAndNonTaut && bTautIsNonTaut &&
                             j == TAUT_NON && 0 <= ( jj = GET_II( OUT_T1, is ) ) &&
                             !is->pINChI[jj]->bDeleted && !is->pINChI[jj]->lenTautomer)
                        {
                            ; /* the requested non-tautomeric component is in tautomeric position   */
                              /*   (is->pINChI[TAUT_YES]);                                          */
                              /*   process it also as non-tautomeric if Fixed-H layer was requested */
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

            if (( pINChI = is->pINChI[jj] ) && pINChI->nNumberOfAtoms > 0)
            {
                /*pINChI_Aux = is->pINChI_Aux[jj];*/
                bCompExists++;

                if (j == TAUT_NON)
                    nAtomsAllComp1 += pINChI->nNumberOfAtoms;
                else if (j == TAUT_YES)
                    nAtomsAllComp2 += pINChI->nNumberOfAtoms;


                bCurTaut = ( pINChI->lenTautomer > 0 );
                bCurIso = ( pINChI->nNumberOfIsotopicAtoms > 0 || pINChI->nNumberOfIsotopicTGroups > 0 );
                bCurIsoHPos = ( (pINChI->nPossibleLocationsOfIsotopicH && pINChI->nPossibleLocationsOfIsotopicH[0] > 1) || pINChI->lenTautomer > 1 ); /* djb-rwth: addressing LLVM warning */
                /* present isotopic H + their possible positions AND/OR isotopic atoms */
                bCurIsoHStereo = (bCurIsoHPos && ( bTautIsoHNum || bTautIsoAt )) || bCurIso; /* djb-rwth: addressing LLVM warning */
                if (jj == j && pINChI->bDeleted)
                {
                    io.num_comp[j] --;
                    if (bCurTaut)
                    {
                        io.bTautomeric |= 1; /* tautomeric representation is present */
                        io.bNonTautomeric |= HAS_N( is );
                    }
                    io.bIsotopic |= bCurIso;
                    continue; /* deleted H(+) in tautomeric representation */
                }

                bCurStereoSp2 = pINChI->Stereo && ( pINChI->Stereo->nNumberOfStereoBonds > 0 );

                bCurHasIsoStereo =
                    bCurStereoSp3 = pINChI->Stereo && ( pINChI->Stereo->nNumberOfStereoCenters > 0 );

                bCurIsoStereoSp2 = bCurIsoHStereo && pINChI->StereoIsotopic && ( pINChI->StereoIsotopic->nNumberOfStereoBonds > 0 );
                bCurIsoStereoSp3 = bCurIsoHStereo && pINChI->StereoIsotopic && ( pINChI->StereoIsotopic->nNumberOfStereoCenters > 0 );
                bCurIsoStereoSp3Inv = bCurIsoStereoSp3 && pINChI->StereoIsotopic->nCompInv2Abs; /* inversion changes sp3 stereo */
                bRequestedRacemicStereo |= ( 0 != ( pINChI->nFlags & INCHI_FLAG_RAC_STEREO ) );
                bRequestedRelativeStereo |= ( 0 != ( pINChI->nFlags & INCHI_FLAG_REL_STEREO ) );

                /* Check whether isotopic stereo is same as non-isotopic; if same than do not output isotopic stereo */
                if (bCurStereoSp2 && bCurIsoStereoSp2)
                {
                    bCurIsoStereoSp2 = !Eql_INChI_Stereo( pINChI->Stereo, EQL_SP2, pINChI->StereoIsotopic, EQL_SP2, 0 );
                }
                if (bCurStereoSp3 && bCurIsoStereoSp3)
                {
                    /* bCurIsoStereoSp3=0 means (iso stereo sp3) = (non-iso stereo sp3) or (iso stereo sp3) = Inv(non-iso stereo sp3) */
                    bCurIsoStereoSp3 = !Eql_INChI_Stereo( pINChI->Stereo, EQL_SP3, pINChI->StereoIsotopic, EQL_SP3,
                        ( pINChI->nFlags & INCHI_FLAG_RAC_STEREO ) || ( pINChI->nFlags & INCHI_FLAG_REL_STEREO ) );
                    if (!bCurIsoStereoSp3)
                    {
                        /* Inversion changes iso sp3 differently from non-iso sp3 Fix11 */
                        bCurIsoStereoSp3Inv &= ( pINChI->StereoIsotopic->nCompInv2Abs != pINChI->Stereo->nCompInv2Abs );
                    }
                }

                bCurRelative = bRequestedRelativeStereo && bCurStereoSp3;
#if ( REL_RAC_STEREO_IGN_1_SC == 1 )
                bCurRelative = bCurRelative &&
                    ( pINChI->Stereo->nNumberOfStereoCenters > 1 ) &&
                    ( pINChI->Stereo->nCompInv2Abs != 0 ) &&
#endif



                    bCurIsoRelative = bRequestedRelativeStereo && ( bCurIsoStereoSp3 || bCurIsoStereoSp3Inv );
#if ( REL_RAC_STEREO_IGN_1_SC == 1 )
                bCurIsoRelative = bCurIsoRelative &&
                    ( pINChI->StereoIsotopic->nNumberOfStereoCenters > 1 ) &&
                    ( pINChI->StereoIsotopic->nCompInv2Abs != 0 ) &&
#endif


                    bCurRacemic = bRequestedRacemicStereo && bCurStereoSp3;
#if ( REL_RAC_STEREO_IGN_1_SC == 1 )
                bCurRacemic = bCurRacemic &&
                    ( pINChI->Stereo->nCompInv2Abs != 0 ) &&
                    ( pINChI->Stereo->nNumberOfStereoCenters > 0 ) ?
                    pINChI->Stereo->nNumberOfStereoCenters : 0;
#endif

                bCurIsoRacemic = bRequestedRacemicStereo && ( bCurIsoStereoSp3 || bCurIsoStereoSp3Inv );
#if ( REL_RAC_STEREO_IGN_1_SC == 1 )
                bCurIsoRacemic = bCurIsoRacemic &
                    ( pINChI->StereoIsotopic->nCompInv2Abs != 0 ) &&
                    ( pINChI->StereoIsotopic->nNumberOfStereoCenters > 0 ) ?
                    pINChI->StereoIsotopic->nNumberOfStereoCenters : 0;
#endif
                if (bRequestedRelativeStereo)
                {
                    bCurStereoSp3 = bCurRelative || (bCurStereoSp3 && ( pINChI->Stereo->nNumberOfStereoCenters > 1 )); /* Fix11 */ /* djb-rwth: addressing LLVM warning */
                    bCurIsoStereoSp3 = bCurIsoRelative ? bCurIsoStereoSp3 : 0;
                }
                else
                {
                    if (bRequestedRacemicStereo)
                    {
                        bCurStereoSp3 = bCurRacemic > 1 || (bCurStereoSp3 && ( pINChI->Stereo->nNumberOfStereoCenters > 1 )); /* Fix11 */ /* djb-rwth: addressing LLVM warning */
                        bCurIsoStereoSp3 = bCurIsoRacemic > 1 ? bCurIsoStereoSp3 : 0;
                    }
                }
                /* djb-rwth: removing redundant code */

                io.bIsotopic |= bCurIso;
                bHasIsotopicAtoms[ii] |= bCurIso;
                bStereoSp2[ii] |= bCurStereoSp2;
                bStereoSp3[ii] |= bCurStereoSp3;
                io.bIgn_UU_Sp3[ii] |= !bCurStereoSp3 && ( pINChI->nFlags & INCHI_FLAG_SC_IGN_ALL_UU );
                io.bIgn_UU_Sp2[ii] |= !bCurStereoSp2 && ( pINChI->nFlags & INCHI_FLAG_SB_IGN_ALL_UU );
                bIsotopicStereoSp2[ii] |= bCurIsoStereoSp2;
                bIsotopicStereoSp3[ii] |= bCurIsoStereoSp3;
                io.bIgn_UU_Sp3_Iso[ii] |= !bCurIsoStereoSp3 && ( pINChI->nFlags & INCHI_FLAG_SC_IGN_ALL_ISO_UU );
                io.bIgn_UU_Sp2_Iso[ii] |= !bCurIsoStereoSp2 && ( pINChI->nFlags & INCHI_FLAG_SB_IGN_ALL_ISO_UU );
                bStereoAbs[ii] |= bCurStereoSp3 && ( pINChI->Stereo->nCompInv2Abs != 0 );

                bStereoAbsInverted[ii] |= bCurStereoSp3 && ( pINChI->Stereo->nCompInv2Abs < 0 );

                /* Fix08: missing isotopic inverted flag if isotopic = inverted non-isotopic */
                bIsotopicStereoAbsInverted[ii] |= (bCurIsoStereoSp3 && ( pINChI->StereoIsotopic->nCompInv2Abs < 0 )) ||
                    (!bCurIsoStereoSp3  && pINChI->StereoIsotopic  && pINChI->Stereo &&
                    pINChI->StereoIsotopic->nCompInv2Abs &&
                    pINChI->StereoIsotopic->nCompInv2Abs != pINChI->Stereo->nCompInv2Abs); /* djb-rwth: addressing LLVM warnings */

                /* Fix 11: missing /s1 if only isotopic stereo is inverted */
                bIsotopicStereoAbs[ii] |= (bCurIsoStereoSp3 && ( pINChI->StereoIsotopic->nCompInv2Abs != 0 )) ||
                    (!bCurIsoStereoSp3  && pINChI->StereoIsotopic  && pINChI->Stereo &&
                    pINChI->StereoIsotopic->nCompInv2Abs &&
                    pINChI->StereoIsotopic->nCompInv2Abs != pINChI->Stereo->nCompInv2Abs); /* djb-rwth: addressing LLVM warnings */

                io.bRelativeStereo[ii] |= bCurRelative;
                io.bIsotopicRelativeStereo[ii] |= bCurIsoRelative;
                io.bRacemicStereo[ii] |= bCurRacemic;
                io.bIsotopicRacemicStereo[ii] |= bCurIsoRacemic;


                bTautomericAcid |= ( 0 != ( pINChI->nFlags & INCHI_FLAG_ACID_TAUT ) );
                bHardAddRemProton |= ( 0 != ( pINChI->nFlags & INCHI_FLAG_HARD_ADD_REM_PROTON ) );
                if (bCurTaut)
                {
                    io.bTautomeric |= 1; /* tautomeric representation is present */
                    /* does tautomeric structure have also a non-tautomeric repesentation? */
                    io.bNonTautomeric |= HAS_N( is );
                }

                /* Auxiliary info */
                if (!( bINChIOutputOptions & INCHI_OUT_NO_AUX_INFO ) && ( pINChI_Aux = is->pINChI_Aux[jj] ))
                {
                    /* detect presence of constitutional equivalence onfo */
                    int bCurEqu, bCurTautEqu = 0, bCurIsoEqu = 0, bCurIsoTautEqu = 0; /* Fix15-disabled */
                    io.bAtomEqu[ii] |= ( bCurEqu = bHasEquString( pINChI_Aux->nConstitEquNumbers,
                        pINChI_Aux->nNumberOfAtoms ) ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
                    if (bCurTaut)
                    {
                        io.bTautEqu[ii] |= ( bCurTautEqu = bHasEquString( pINChI_Aux->nConstitEquTGroupNumbers,
                            pINChI_Aux->nNumberOfTGroups ) ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
                    }
                    if (bCurIso)
                    {
                        io.bIsotopicAtomEqu[ii] |= ( bCurIsoEqu = bHasEquString( pINChI_Aux->nConstitEquIsotopicNumbers,
                            pINChI_Aux->nNumberOfAtoms ) ) /*|| bCurEqu*/; /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
                        if (bCurTaut)
                        {
                            io.bIsotopicTautEqu[ii] |= ( bCurIsoTautEqu = bHasEquString( pINChI_Aux->nConstitEquIsotopicTGroupNumbers,
                                pINChI_Aux->nNumberOfTGroups ) ) /*|| bCurTautEqu*/; /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
                        }
                        /* non-zero if isotopic numbering for inverted isotopic stereo is different */
                        io.bIsotopicOrigNumb[ii] |= bCurHasIsoStereo && /* Fix14 */
                            pINChI_Aux->nOrigAtNosInCanonOrdInv &&
                            pINChI_Aux->nIsotopicOrigAtNosInCanonOrd &&
                            ( 0 != memcmp( pINChI_Aux->nOrigAtNosInCanonOrdInv,
                                pINChI_Aux->nIsotopicOrigAtNosInCanonOrd,
                                sizeof( pINChI_Aux->nOrigAtNosInCanonOrdInv[0] ) * pINChI_Aux->nNumberOfAtoms ) );
                    }
                    /* Inverted stereo */
                    if (bCurStereoSp3 && pINChI->Stereo->nCompInv2Abs)
                    {
                        io.bInvStereo[ii] |= 1;
                        io.bInvStereoOrigNumb[ii] |= pINChI_Aux->nOrigAtNosInCanonOrd &&
                            pINChI_Aux->nOrigAtNosInCanonOrdInv &&
                            ( 0 != memcmp( pINChI_Aux->nOrigAtNosInCanonOrd,
                                pINChI_Aux->nOrigAtNosInCanonOrdInv,
                                sizeof( pINChI_Aux->nOrigAtNosInCanonOrd[0] ) * pINChI_Aux->nNumberOfAtoms ) );
                    }

                    /* Inverted isotopic stereo */
                    if (bCurIsoStereoSp3 && pINChI->StereoIsotopic->nCompInv2Abs)
                    {
                        io.bInvIsotopicStereo[ii] |= 1;

                        io.bInvIsotopicStereoOrigNumb[ii]
                            |= pINChI_Aux->nIsotopicOrigAtNosInCanonOrd &&
                            pINChI_Aux->nIsotopicOrigAtNosInCanonOrdInv &&
                            ( 0 != memcmp( pINChI_Aux->nIsotopicOrigAtNosInCanonOrd,
                                pINChI_Aux->nIsotopicOrigAtNosInCanonOrdInv,
                                sizeof( pINChI_Aux->nIsotopicOrigAtNosInCanonOrd[0] ) * pINChI_Aux->nNumberOfAtoms ) );
                    }

                    if (pINChI_Aux->OrigInfo && bHasOrigInfo( pINChI_Aux->OrigInfo, pINChI_Aux->nNumberOfAtoms ))
                    {
                        io.bChargesRadVal[ii] |= 1;
                    }
                }
            }
        }
        if (bCompExists)
        {
            for (j = TAUT_NON; j < TAUT_NUM; j++)
            {
                io.num_comp[j] ++;
            }
        }
    }
    if (io.bTautomeric /*&& bTautomericAcid*/) /* "&& bTautomericAcid" commented out 2004-06-02 */
    {
        io.bTautomeric += bTautomericAcid; /* long-range tautomerism */
        io.bTautomeric += ( bHardAddRemProton ? 4 : 0 );
    }
    if (bRequestedRacemicStereo || bRequestedRelativeStereo)
    {
        /* do not output inverted stereo info */
        for (i = 0; i < TAUT_NUM; i++)
        {
            /* Fix11 */
            bStereoAbsInverted[i] =
                bStereoAbs[i] =
                io.bInvStereo[i] =
                io.bInvStereoOrigNumb[i] = 0;
                /* io.bIsotopicRelativeStereo[i]=0 may happen because iso stereo is same or inverted non-iso stereo */
            bIsotopicStereoAbsInverted[i] =
                bIsotopicStereoAbs[i] =
                io.bInvIsotopicStereo[i] =
                io.bInvIsotopicStereoOrigNumb[i] = 0;
        }
    }


    io.iCurTautMode = io.bOutType == OUT_N1 ? TAUT_NON :  /*  only non-taut */

        io.bOutType == OUT_T1 ? TAUT_YES :      /*  tautomeric if present, otherwise non-tautomeric     */
        io.bOutType == OUT_NT ? TAUT_NON :      /*  only non-taut representations of tautomeric         */
        io.bOutType == OUT_TN ? TAUT_YES :       /*  tautomeric if present otherwise non-tautomeric;     */
        -1; /*  separately output non-taut representations of tautomeric if present */

    if (io.iCurTautMode < 0)
    {
        return 0;  /* error */
    }


    /* Now print out */

    io.bOverflow = 0;
    io.num_components = io.num_comp[io.iCurTautMode];
    /* djb-rwth: removing redundant code */

    if (bINChIOutputOptions & INCHI_OUT_ONLY_AUX_INFO)
    {
        goto output_aux_info;
    }

    io.nCurINChISegment = DIFL_M;

    /* InChI output: version and kind */
    if (INCHI_basic_or_INCHI_reconnected == INCHI_BAS || !( bINChIOutputOptions & INCHI_OUT_EMBED_REC ))
    {
        int is_beta = 0;
        int nAtomsAllComp = inchi_max( nAtomsAllComp1, nAtomsAllComp2 );

        if (nAtomsAllComp > NORMALLY_ALLOWED_INP_MAX_ATOMS)
        {
            /* v. 1.05 for LargeMolecules */
            is_beta = 1;
        }
        if (pOrigStruct && pOrigStruct->polymer)
        {
            /* v. 1.05 for Polymers */
            is_beta = 1;
        }
        /* specifically put 'B' for empty structure InChI    */
        /* if "Polymers" or "LargeMolecules" requested        */
        else if (!pOrigStruct && ( ip->bLargeMolecules || ip->bPolymers ))
        {
            is_beta = 1;
        }
        else if (pOrigStruct && pOrigStruct->n_zy)
        {
            is_beta = 1;
        }

        OutputINCHI_VersionAndKind( out_file, strbuf, bINChIOutputOptions, is_beta, pLF, pTAB );
    }



    /* InChI output: atoms */
    intermediate_result = OutputINCHI_MainLayerFormula( pCG, out_file, strbuf,
                                                        num_components2,
                                                        &INCHI_basic_or_INCHI_reconnected,
                                                        &io, pLF, pTAB );
    if (intermediate_result != 0)
        goto exit_function;

    /* InChI output: connection table */
    intermediate_result = OutputINCHI_MainLayerConnections( pCG, out_file, strbuf, num_components2,
                                                                &INCHI_basic_or_INCHI_reconnected,
                                                                &io, pLF, pTAB );
    if (intermediate_result != 0)
        goto exit_function;

    /* InChI output: hydrogens (with tautomeric info) */
    intermediate_result = OutputINCHI_MainLayerHydrogens( pCG, out_file, strbuf, num_components2,
                                                              &INCHI_basic_or_INCHI_reconnected,
                                                              &io, pLF, pTAB );
    if (intermediate_result != 0)
        goto exit_function;

    io.bFhTag = 0;
    npass = 0;


repeat_INChI_output:

    /* InChI output: charge and  removed protons */
    intermediate_result = OutputINCHI_ChargeAndRemovedAddedProtonsLayers( pCG, out_file, strbuf,
                                                                              &io, pLF, pTAB );
    if (intermediate_result != 0)
        goto exit_function;

    /* InChI output: polymer layer */
    if (npass == 0)
    {
        intermediate_result = OutputINCHI_PolymerLayer( pCG, out_file, strbuf,
                                                        &INCHI_basic_or_INCHI_reconnected,
                                                        orig_inp_data, pOrigStruct,
                                                        &io, pLF, pTAB );
        if (intermediate_result != 0)
            goto exit_function;
    }

    /* InChI output: stereo (non-isotopic) */
    intermediate_result = OutputINCHI_StereoLayer( pCG, out_file, strbuf, &io, pLF, pTAB );
    if (intermediate_result != 0)
        goto exit_function;


    /* Switch from M to MI or from F to FI */
    io.nCurINChISegment++;

    /* InChI output: isotopic */
    intermediate_result = OutputINCHI_IsotopicLayer( pCG, out_file, strbuf,
                                                         &INCHI_basic_or_INCHI_reconnected,
                                                         &io, pLF, pTAB );

    if (intermediate_result != 0)
    {
        goto exit_function;
    }



    /*
        At this point the INChI part of the output has been done.
        If this INChI is tautomeric and non-tautomeric results exist,
        then we need to output non-tautomeric data:
            fixed H
            stereo
            isotopic
            isotopic stereo
    */


    /* InChI output: FixedH and sublayers */
    intermediate_result = OutputINCHI_FixedHLayerWithSublayers( pCG, out_file, strbuf,
                                                                &INCHI_basic_or_INCHI_reconnected,
                                                                &io, pLF, pTAB,
                                                                &then_goto_repeat );
    if (intermediate_result != 0)
    {
        goto exit_function;
    }
    if (then_goto_repeat)
    {
        npass++;
        goto repeat_INChI_output;
    }


    /*
        InChI output:  reconnected structure
    */

    bEmbeddedOutputCalled = 0;
    if (bDisconnectedCoord && INCHI_basic_or_INCHI_reconnected == INCHI_BAS &&
        ( bINChIOutputOptions & INCHI_OUT_EMBED_REC ) && num_components2[INCHI_REC])
    {
        int nRet;
        bEmbeddedOutputCalled = 1;

        /* output blank line before /R: in case of bPlainTextCommnts=1 */
        inchi_ios_print_nodisplay( out_file, "%s", pLF );
        /* end of disconnected INChI output */

        nRet = OutputINChI1( pCG,
                             strbuf,
                             pINChISortTautAndNonTaut2,
                             INCHI_REC,
                             orig_inp_data,
                             pOrigStruct,
                             ip,
                             0 /*bDisconnectedCoord*/,
                             bOutputType,
                             bINChIOutputOptions | INCHI_OUT_NO_AUX_INFO,
                             num_components2,
                             num_non_taut2,
                             num_taut2,
                             out_file,
                             log_file,
                             num_input_struct,
                             pSortPrintINChIFlags,
                             save_opt_bits );

        if (!nRet)
        {
            goto exit_function; /* error */
        }
    }

    /* InChI output: save InChI creation options if requested */
    if (!bEmbeddedOutputCalled &&
        ( bINChIOutputOptions & INCHI_OUT_SAVEOPT ) &&
        ( 0 == ( bINChIOutputOptions & INCHI_OUT_STDINCHI ) )    /* not std-InChI output */
        )
    {
        char let1, let2;
        GetSaveOptLetters( save_opt_bits, &let1, &let2 );
        inchi_ios_print_nodisplay( out_file, "\\%c%c", let1, let2 );
    }
    if (!bEmbeddedOutputCalled && !bPlainTextCommnts)
    { /* plain text comment earlier ended with LF */
        inchi_ios_print_nodisplay( out_file, "%s%s",
            ( !num_components2[0] && !num_components2[1] ) ? "//" : "", /* empty InChI=// */
            ( bINChIOutputOptions & INCHI_OUT_NO_AUX_INFO ) ? "\n" : pTAB );
/* end of INChI= output */
    }

    inchi_strbuf_reset( strbuf );

#ifdef TARGET_LIB_FOR_WINCHI
    /* @@@ Here we end up with silent output: display previously hidden output */
    if (inchi_ios_flush_not_displayed( out_file ) != -1)
        silent = 0;
#endif



output_aux_info:

    /*  Output Aux Info */

    io.bFhTag = 0;
    if (( bINChIOutputOptions & INCHI_OUT_NO_AUX_INFO ) == 0)
    {

        io.num_components = io.num_comp[io.iCurTautMode];

        /* AuxInfo: header and normalization type */
        intermediate_result = OutputAUXINFO_HeaderAndNormalization_type( pCG, out_file, strbuf,
                                                                         bINChIOutputOptions,
                                                                         &INCHI_basic_or_INCHI_reconnected,
                                                                         num_components2,
                                                                         &io, pLF, pTAB );
        if (intermediate_result != 0)
            goto exit_function;


    repeat_INChI_Aux_output:

        /* AuxInfo: original atom numbers and symmetry numbers (constit. equivalence /E: )    */
        intermediate_result = OutputAUXINFO_OriginalNumbersAndEquivalenceClasses( pCG, out_file, strbuf,
                                                                                  num_components2,
                                                                                  &io, pLF, pTAB );
        if (intermediate_result != 0)
            goto exit_function;

        /* AuxInfo: tautomeric groups equivalence */
        intermediate_result = OutputAUXINFO_TautomericGroupsEquivalence( pCG, out_file, strbuf, &io );
        if (intermediate_result != 0)
            goto exit_function;

        /* AuxInfo: stereo data */
        intermediate_result = OutputAUXINFO_Stereo( pCG, out_file, strbuf, &io, pLF, pTAB );
        if (intermediate_result != 0)
            goto exit_function;

    repeat_INChI_Aux_Iso_output:
            /* AuxInfo: isotopic info */
        intermediate_result = OutputAUXINFO_IsotopicInfo( pCG, out_file, strbuf,
                                                          &INCHI_basic_or_INCHI_reconnected,
                                                          &io, pLF, pTAB );
        if (intermediate_result != 0)
        {
            goto exit_function;
        }


        /*
          At this point the INChI_Aux part of the output has been completed.
          If this INChI is tautomeric and non-tautomeric results exist,
          then we need to output non-tautomeric auxilialy data
          (same as above excluding tautomeric information).
          Currently, this is enabled for xml output only
        */

        if (io.bOutType == OUT_TN && io.bTautomeric && io.bNonTautomeric &&
            /* Check whether the Fixed-H layer is empty */
            ( *pSortPrintINChIFlags & ( ( INCHI_basic_or_INCHI_reconnected == INCHI_BAS ) ? FLAG_SORT_PRINT_NO_NFIX_H_BAS :
                FLAG_SORT_PRINT_NO_NFIX_H_REC ) ) &&
                ( *pSortPrintINChIFlags & ( ( INCHI_basic_or_INCHI_reconnected == INCHI_BAS ) ? FLAG_SORT_PRINT_NO_IFIX_H_BAS :
                    FLAG_SORT_PRINT_NO_IFIX_H_REC ) )
              )
        {
            io.bNonTautomeric = 0; /* bNonTautIdentifierNotEmpty == 0 => no fixed H info 02-10-2995 */
        }

        if (io.bOutType == OUT_TN && io.bTautomeric && io.bNonTautomeric)
        {
            /* add the second (non-tautomeric) output */
            io.bOutType = OUT_NONTAUT;
            io.iCurTautMode = TAUT_NON;
            io.pINChISort = io.pINChISortTautAndNonTaut[TAUT_NON];
            io.bSecondNonTautPass = 1;
            io.num_components = io.num_comp[io.iCurTautMode];
            io.bFhTag = AL_FIXH;
            inchi_strbuf_reset( strbuf ); /*pStr[io.tot_len=0] = '\0';*/

            /* if InChI Fixed-H isotopic is empty then do not output corresponding AuxInfo */
            if (!( *pSortPrintINChIFlags &
                ( ( INCHI_basic_or_INCHI_reconnected == INCHI_BAS ) ? FLAG_SORT_PRINT_NO_NFIX_H_BAS :
                    FLAG_SORT_PRINT_NO_NFIX_H_REC ) )
               )
            {
                npass++;
                goto repeat_INChI_Aux_output;
            }
            else
            {
                npass++;
                goto repeat_INChI_Aux_Iso_output;
            }
        }
        else
        {
            if (io.bOutType == OUT_NONTAUT && io.bOutputType == OUT_TN && io.bTautomeric && io.bNonTautomeric)
            {
                /* the second (non-taut) output has been done; restore variables */
                io.bOutType = OUT_TN;
                io.iCurTautMode = TAUT_YES;
                io.pINChISort = io.pINChISortTautAndNonTaut[TAUT_YES];
                io.bSecondNonTautPass = 0;
                /* set correct num components for the reversibility info 02-10-2005 */
                io.num_components = io.num_comp[io.iCurTautMode];
                io.bFhTag = 0;
            }
        }

        /*    Charges, radicals, unusual valences */
        intermediate_result = OutputAUXINFO_ChargesRadicalsAndUnusualValences( pCG, out_file, strbuf, &io, pLF, pTAB );
        if (intermediate_result != 0)
        {
            goto exit_function;
        }


        /* Output the original input structure -- quick fix */
        intermediate_result = OutputAUXINFO_ReversibilityInfo( pCG, out_file, strbuf, pOrigStruct, &io, pLF, pTAB );
        if (intermediate_result != 0)
        {
            goto exit_function;
        }

        /* Output polymeric Aux Info */
        intermediate_result = OutputAUXINFO_PolymerInfo( pCG, out_file, strbuf, pOrigStruct, &io, pLF, pTAB );
        if (intermediate_result != 0)
        {
            goto exit_function;
        }

        /*
            Output INChI_Aux of the reconnected structure
        */

        bEmbeddedOutputCalled = 0;
        if (bDisconnectedCoord && INCHI_basic_or_INCHI_reconnected == INCHI_BAS && ( bINChIOutputOptions & INCHI_OUT_EMBED_REC ) &&
             num_components2[INCHI_REC] && !( bINChIOutputOptions & INCHI_OUT_NO_AUX_INFO ))
        {
            int nRet;
            bEmbeddedOutputCalled = 1;
            inchi_ios_print( out_file, "%s", pLF );

            nRet = OutputINChI1( pCG,
                                 strbuf,
                                 pINChISortTautAndNonTaut2,
                                 INCHI_REC,
                                 NULL,
                                 NULL,
                                 ip,
                                 0 /*bDisconnectedCoord*/,
                                 bOutputType,
                                 INCHI_OUT_ONLY_AUX_INFO | bINChIOutputOptions,
                                 num_components2,
                                 num_non_taut2,
                                 num_taut2,
                                 out_file,
                                 log_file,
                                 num_input_struct,
                                 pSortPrintINChIFlags,
                                 save_opt_bits );

            if (!nRet)
            {
                goto exit_function; /* error */
            }
        }

        /* Close INChI_Aux */
        if (!bEmbeddedOutputCalled && !bPlainTextCommnts)
        {
            inchi_ios_print( out_file, "%s\n", ( !num_components2[0] && !num_components2[1] ) ? "//" : "" );
            /* plain text comment earlier ended with LF */
        }

        /* in wINChI window, separate AuxInfo: from InChIKey: with blank line */
        inchi_ios_print( out_file, "%s",
            ( bINChIOutputOptions & INCHI_OUT_WINCHI_WINDOW ) ? "\n" : "" );
    } /* end of output AuxInfo */

    ret = 1;


exit_function:


#ifdef TARGET_LIB_FOR_WINCHI
    /* @@@ If for any error we get here silent, display previously hidden output */
    if (silent)
    {
     /*
        if ( !inchi_ios_flush_not_displayed( out_file ) != -1  )
            silent = 0;
    */
        silent = 0;
    }
#endif

    if (io.bOverflow)
    {
        inchi_ios_print( out_file, "\nFATAL ERROR: Output buffer overflow\n" );
    }

    if (intermediate_result)
    {
        ret = 0;
        inchi_ios_eprint( log_file, "InChI serialization error for structure #%d.%s%s%s%s\n",
                                    num_input_struct, SDF_LBL_VAL( ip->pSdfLabel, ip->pSdfValue ) );
    }

    return ret;
} /* OutputINChI1 */


/****************************************************************************/
char *szGetTag( const INCHI_TAG *Tag,
                int             nTag,
                int             bTag,
                char            *szTag,
                int             *bAlways,
                short           tag_flag)
{
    int i, j, bit, num, len;
    const int MAX_TAG_NUM = tag_flag ? (int)IL_MAX_ORD : (int)AL_MAX_ORD; /* djb-rwth: fixing GHI #160 */
    if (0 < nTag && nTag < 3)
    {
        /* no plain text comments: pick up the last tag */
        for (i = 0, j = -1, bit = 1; i < MAX_TAG_NUM; i++, bit <<= 1)
        {
            if (bTag & bit)
            {
                j = i;
            }
        }
        if (j >= 0)
        {
#if USE_BCF
            int stl1, stl2, dstsz;
            stl1 = strlen(Tag[j].szXmlLabel) + 1;
            stl2 = strlen(Tag[j].szPlainLabel) + 1;
            dstsz = max_3(stl1, stl2, 5);
            strcpy_s( szTag, dstsz, nTag == 1 ? Tag[j].szXmlLabel : nTag == 2 ? Tag[j].szPlainLabel : "???" ); /* djb-rwth: function replaced with its safe C11 variant */
#else
            strcpy(szTag, nTag == 1 ? Tag[j].szXmlLabel : nTag == 2 ? Tag[j].szPlainLabel : "???"); /* djb-rwth: addressing coverity ID #499488 -- when nTag == 2, the "???" is avoided, which is correct */
#endif
            if (nTag != 2)
            {
                *bAlways = Tag[j].bAlwaysOutput;
            }
            return szTag;
        }
    }
    else
        if (nTag == 3)
        {
            /* plain text with comments */
            szTag[0] = '{';
            szTag[1] = '\0';
            for (i = 0, j = -1, bit = 1, num = 0; i < MAX_TAG_NUM; i++, bit <<= 1)
            {
                if (bTag & bit)
                {
                    j = i;
                    if (num++)
                    {
                        strcat(szTag, ":");
                    }
                    strcat(szTag, Tag[i].szPlainComment);
                }
            }
            if (num)
            {
                strcat(szTag, "}");
                num = (int) strlen( Tag[j].szPlainLabel );
                len = (int) strlen( szTag );
                if (len)
                {
                    memmove(szTag + num, szTag, (long long)len + 1); /* djb-rwth: cast operator added */
                    memcpy(szTag, Tag[j].szPlainLabel, num);
                }
                else
                {
                    strcpy(szTag, Tag[j].szPlainLabel);
                }
                *bAlways = Tag[j].bAlwaysOutput;
            }
            else
            {
                strcpy(szTag, "???");
            }
            return szTag;
        }

    strcpy(szTag, "???");
    return szTag;
}


/* djb-rwth: removing redundant code */


/****************************************************************************
    str_LineEnd( ... )

    First, checks if buffer overflow; then:
    if ind < 0 (common usage, plain text output)
        sets terminating '\0' in pStr,
        optionally adds leading tag (e.g., '/' or "/c" )
    *obsolete* if ind >=0 XML output

****************************************************************************/
int str_LineEnd( const char       *tag,
                 int              *bOverflow,
                 INCHI_IOS_STRING *buf,
                 int               ind,
                 int               bPlainTextTags )
{
    /* djb-rwth: removing redundant variables */
    int tag_len;

    /* check buffer overflow */
    if (*bOverflow)
    {
        return 1;
    }

    if (ind < 0)
    {
        /* Plain text */
        /* insert plain text tag if:
           (a) pStr has non-zero length, or
           (b) ind < -1
        */
        if (buf->pStr[0] || ind < -1)
        {
            tag_len = bPlainTextTags ? (int) strlen( tag ) : 0;
            if (tag_len > 0)
            {
                int n_added = tag_len + 2 + 2;
                inchi_strbuf_update( buf, n_added );

                memmove(buf->pStr + tag_len, buf->pStr, (long long)buf->nUsedLength + 1); /* djb-rwth: cast operator added */
                /* NB: trailing 0 is also memmoved */
                memcpy(buf->pStr, tag, tag_len);

                /* to be sure...  */
                buf->nUsedLength = strlen( buf->pStr );
            }
        }
    }

    return 0;
}


/****************************************************************************/
int CleanOrigCoord( MOL_COORD szCoord, int delim )
{
#define MIN_BOND_LENGTH   (1.0e-6)
    char szVal[LEN_COORD + 1];
    MOL_COORD szBuf;
    char *q;
    int len, last, fst, dec_pnt, num_zer = 0, len_buf = 0, e;
    int k, i;
    double coord;

    for (k = 0; k < NUM_COORD*LEN_COORD; k += LEN_COORD)
    {
        memcpy(szVal, szCoord + k, LEN_COORD);
        szVal[LEN_COORD] = '\0';
        lrtrim( szVal, &len );
        coord = strtod( szVal, &q );
        if (MIN_BOND_LENGTH > fabs( coord ))
        {
            strcpy(szVal, "0");
            len = 1;
            num_zer++;
        }
        else
        {
            len = (int) ( q - szVal );
            /* last = (last mantissa digit position + 1)  */
            if (( q = strchr( szVal, 'e' ) ) || ( q = strchr( szVal, 'E' ) ) ||
                ( q = strchr( szVal, 'd' ) ) || ( q = strchr( szVal, 'D' ) ))
            {
                /* floating point */
                last = q - szVal;
                /* remove (+) and leading zeroes from the exponent */
                e = (int) strtol( szVal + last + 1, &q, 10 ); /* exponent */
                if (e)
                {
                    /* new exp; update the length */
                    len = last + 1 + sprintf(szVal + last + 1, "%d", e); /* print exp without leading zeroes and '+' */
                }
                else
                {
                    /* exponent is zero */
                    len = last;
                }
            }
            else
            {
                last = len;
            }
            /* fst = (first mantissa digit); fst=1 if the sign is present, otherwise 0 */
            fst = ( szVal[0] != '.' && !isdigit( UCINT szVal[0] ) );
            /* dec_pnt = (decimal point position) or last */
            if ((q = strchr( szVal, '.' ))) /* djb-rwth: addressing LLVM warning */
            {
                dec_pnt = (int) ( q - szVal );
            }
            else
            {
                dec_pnt = last;
            }
            last -= 1; /* last mantissa digit position */
            /* remove trailing zeroes in the range dec_pnt+1..last-1 */
            for (i = last; dec_pnt < i && '0' == szVal[i]; i--)
                ;
            if (i == dec_pnt)
            {
                i--; /* remove decimal point, too */
            }
            if (i < last)
            {
                memmove(szVal + i + 1, szVal + last + 1, (long long)len - (long long)last); /* djb-rwth: cast operator added */
                len -= last - i;
            }
            /* remove leading zeroes */
            for (i = fst; i < len && '0' == szVal[i]; i++)
            {
                ;
            }
            if ((i > fst) && (len - fst <= LEN_COORD + 1 - i) && (len - fst <= LEN_COORD + 1 - fst)) /* djb-rwth: fixing GHI #138 */
            {
                memmove(szVal + fst, szVal + i, (long long)len - (long long)fst); /* djb-rwth: cast operator added */
                len -= i - fst;
            }
        }
        if (len_buf && (len_buf < (int)sizeof(MOL_COORD)))
        {
#pragma warning (push)
#pragma warning (disable: 6386)
            szBuf[len_buf++] = delim;
#pragma warning (pop)
        }
        if (len_buf >= (int)sizeof(MOL_COORD)) /* djb-rwth: fixing coverity ID #499520 */
        {
            len_buf = (int)sizeof(MOL_COORD) - 1;
            len = 0;
        }
        memcpy(szBuf + len_buf, szVal, len); /* does not copy zero termination*/
        len_buf += len;
    }
    /* zero termination */
    if (len_buf < ( int )sizeof( MOL_COORD ))
    {
        memset( szBuf + len_buf, 0, sizeof( MOL_COORD ) - len_buf ); /* djb-rwth: memset_s C11/Annex K variant? */
    }

    memcpy(szCoord, szBuf, sizeof(MOL_COORD));

    return num_zer;
#undef MIN_BOND_LENGTH
}


/****************************************************************************/
int WriteOrigCoord( int       num_inp_atoms,
                    MOL_COORD *szMolCoord,
                    int       *i,
                    char      *szBuf,
                    int       buf_len )
{

    int j, num_zer, len, cur_len;
    char *p;
    MOL_COORD szCurCoord;
    cur_len = 0;
    for (j = *i; j < num_inp_atoms; )
    {
        memcpy(szCurCoord, szMolCoord[j], sizeof(szCurCoord));
        num_zer = CleanOrigCoord( szCurCoord, ',' );
        if (NUM_COORD == num_zer)
        {
            len = 0;
        }
        else
        {
            if ((p = (char *) memchr( szCurCoord, '\0', sizeof( szCurCoord ) ))) /* djb-rwth: addressing LLVM warning */
            {
                len = (int) ( p - szCurCoord );
            }
            else
            {
                len = sizeof( szCurCoord );
            }
        }
        if (len + cur_len + 1 < buf_len)
        {
            if (len)
            {
                memcpy(szBuf + cur_len, szCurCoord, len * sizeof(szBuf[0]));
            }
            szBuf[cur_len += len] = ';';
            cur_len++;
            j++;
        }
        else
        {
            break;
        }
    }
    szBuf[cur_len] = '\0';
    *i = j; /* next item */

    return cur_len;
}


/****************************************************************************
  WriteOrigAtoms

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
****************************************************************************/
int WriteOrigAtoms( CANON_GLOBALS *pCG,
                    int           num_inp_atoms,
                    inp_ATOM      *at,
                    int           *i,
                    char          *szBuf,
                    int           buf_len,
                    STRUCT_DATA   *sd )
{
    int j, k, n, len, len0, cur_len, val, bonds_val, mw, parity, num_trans, is_ok, b_self;
    static const char szIsoH[] = "hdt";
    char szCurAtom[32];
    AT_NUMB nNeighOrder[MAXVAL], neigh;

    cur_len = 0;
    if (0 == *i)
    {
        cur_len = sprintf(szBuf, "%d%s", num_inp_atoms,
            (sd->bChiralFlag & FLAG_INP_AT_CHIRAL) ? "c" :
            (sd->bChiralFlag & FLAG_INP_AT_NONCHIRAL) ? "n" : "");
    }
    for (j = *i; j < num_inp_atoms; )
    {
        /* tetrahedral parity treatment */
        parity = 0;
        num_trans = 0;
        if (at[j].p_parity)
        {
            /* verify neighbors */
            is_ok = 1;
            b_self = 0;
            for (n = 0, k = 0; n < MAX_NUM_STEREO_ATOM_NEIGH; n++)
            {
                neigh = at[j].p_orig_at_num[n] - 1;
                if (is_in_the_list( at[j].neighbor, neigh, at[j].valence ) &&
                     at[neigh].orig_at_number == at[j].p_orig_at_num[n])
                {
                    /* real neighbor */
                    nNeighOrder[k++] = at[j].p_orig_at_num[n];
                }
                else
                {
                    if ((int) neigh == j && at[neigh].orig_at_number == at[j].p_orig_at_num[n])
                    {
                        /* central atom is a neighbor */
                        num_trans = n; /* move this neighbor to 0 position permutation parity */
                        b_self++;
                    }
                    else
                    {
                        is_ok = 0;
                        break;
                    }
                }
            }
            if (is_ok && b_self <= 1 && b_self + k == MAX_NUM_STEREO_ATOM_NEIGH)
            {
                num_trans += insertions_sort( pCG, nNeighOrder, k, sizeof( nNeighOrder[0] ), comp_AT_RANK );
                if (ATOM_PARITY_WELL_DEF( at[j].p_parity ))
                {
                    parity = 2 - ( num_trans + at[j].p_parity ) % 2;
                }
                else
                {
                    if (ATOM_PARITY_ILL_DEF( at[j].p_parity ))
                    {
                        parity = at[j].p_parity;
                    }
                    else
                    {
                        ; /* invalid atom parity */
                    }
                }
            }
            else
            {
                ;/* add error message here */
            }
        }

        len = len0 = (int) strlen( at[j].elname );

        memcpy(szCurAtom, at[j].elname, len);
        bonds_val = nBondsValenceInpAt( at + j, NULL, NULL );

        if (( val = needed_unusual_el_valence( at[j].el_number, at[j].charge, at[j].radical,
            at[j].chem_bonds_valence, bonds_val, at[j].num_H, at[j].valence ) ) ||
             at[j].charge || at[j].radical || at[j].iso_atw_diff || NUM_ISO_H( at, j ) || parity)
        {
            /* valence */
            if (val)
            {
                len += sprintf(szCurAtom + len, "%d", val > 0 ? val : 0);
            }
            /* charge */
            if ((val = at[j].charge)) /* djb-rwth: addressing LLVM warning */
            {
                szCurAtom[len++] = val > 0 ? '+' : '-';
                if (( val = abs( val ) ) > 1)
                {
                    len += sprintf(szCurAtom + len, "%d", val);
                }
            }
            /* radical */
            if ((val = at[j].radical)) /* djb-rwth: addressing LLVM warning */
            {
                len += sprintf(szCurAtom + len, ".%d", val);
            }
            /* isotopic shift */
            if ((val = at[j].iso_atw_diff)) /* djb-rwth: addressing LLVM warning */
            {
                mw = get_atomic_mass_from_elnum( at[j].el_number );
                if (val == 1)
                    val = mw;
                else
                    if (val > 0)
                        val = mw + val - 1;
                    else
                        val = mw + val;

                len += sprintf(szCurAtom + len, "%si%d", len == len0 ? "." : "", val);
            }
            /* parity */
            if (parity)
            {
                len += sprintf(szCurAtom + len, "%s%s", len == len0 ? "." : "",
                    parity == AB_PARITY_ODD ? "o" :
                    parity == AB_PARITY_EVEN ? "e" :
                    parity == AB_PARITY_UNKN ? "u" :
                    parity == AB_PARITY_UNDF ? "?" : "");
            }
            /* implicit isotopic H */
            if (NUM_ISO_H( at, j ))
            {
                for (k = 0; k < NUM_H_ISOTOPES; k++)
                {
                    if ((val = at[j].num_iso_H[k])) /* djb-rwth: addressing LLVM warning */
                    {
                        len += sprintf(szCurAtom + len, "%s%c", len == len0 ? "." : "", szIsoH[k]);
                        if (val > 1)
                        {
                            len += sprintf(szCurAtom + len, "%d", val);
                        }
                    }
                }
            }
        }

        if (len + cur_len < buf_len)
        {
            memcpy(szBuf + cur_len, szCurAtom, len);
            cur_len += len;
            j++;
        }
        else
        {
            break;
        }
        szBuf[cur_len] = '\0';
        *i = j;
    }

    return cur_len;
}


/****************************************************************************
 WriteOrigBonds( ... )

    Output bonds in ascending order of the neighboring atom original numbers

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
****************************************************************************/
int WriteOrigBonds( CANON_GLOBALS *pCG,
                    int           num_inp_atoms,
                    inp_ATOM      *at,
                    int           *i,
                    char          *szBuf,
                    int           buf_len,
                    STRUCT_DATA   *sd )
{
    int j, k, k2, kk, len, cur_len, j2 = 0, bond_stereo, bond_char, bond_parity, bond_parityNM, num_trans; /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    char szCurBonds[7 * MAXVAL + 2]; /* num_neigh*(1 byte bond type + 2 bytes for bond parity up to 4 digits per neighbor number) + at the end one ';' */
    AT_RANK nNeighOrder[MAXVAL];
    int  chain_len, pnxt_atom, pinxt2cur, pinxt_sb_parity_ord;
    int  chain_len2, pnxt_atom2, pinxt2cur2, pinxt_sb_parity_ord2, m1, m2;
    int  pcur_atom, picur2nxt, picur_sb_parity_ord;

    cur_len = 0;
    for (j = *i; j < num_inp_atoms; )
    {
        len = 0;
        if (at[j].valence >= 1) /* djb-rwth: changing condition to avoid garbage values */
        {
            for (k = 0; k < at[j].valence; k++)
            {
                nNeighOrder[k] = k;
            }
            pCG->m_pn_RankForSort = at[j].neighbor;
            num_trans = insertions_sort( pCG, nNeighOrder, at[j].valence, sizeof( nNeighOrder[0] ), CompRank ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
        }
        else
        {
            num_trans = 0; /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
            nNeighOrder[0] = 0;
        }
        for (kk = 0; kk < at[j].valence; kk++) 
        {
            k = nNeighOrder[kk];
            j2 = at[j].neighbor[k];
            bond_parity = 0;
            bond_parityNM = 0;
            if (j2 < j)
            {
                bond_stereo = at[j].bond_stereo[k];
                switch (at[j].bond_type[k])
                {
                    case BOND_TYPE_SINGLE:
                        switch (bond_stereo)
                        {
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
                        switch (bond_stereo)
                        {
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
                k2 = (int) ( is_in_the_list( at[j2].neighbor, (AT_NUMB) j, at[j2].valence ) - at[j2].neighbor );
                chain_len = chain_len2 = 0;
                if (at[j].sb_parity[0])
                {
                    for (m1 = 0; m1 < MAX_NUM_STEREO_BONDS && at[j].sb_parity[m1]; m1++)
                    {
                        if (k == at[j].sb_ord[m1])
                        {
                            chain_len = get_opposite_sb_atom( at, j, k,
                                          &pnxt_atom, &pinxt2cur, &pinxt_sb_parity_ord );
                            break;
                        }
                    }
                }
                if (at[j2].sb_parity[0])
                {
                    for (m2 = 0; m2 < MAX_NUM_STEREO_BONDS && at[j2].sb_parity[m2]; m2++)
                    {
                        if (k2 == at[j2].sb_ord[m2])
                        {
                            chain_len2 = get_opposite_sb_atom( at, j2, k2,
                                           &pnxt_atom2, &pinxt2cur2, &pinxt_sb_parity_ord2 );
                            break;
                        }
                    }
                }
                if ((chain_len == 1 && chain_len2 == 1) ||  /* regular stereobond */
                     (chain_len > 1 && j > pnxt_atom)) /* djb-rwth: addressing LLVM warnings */
                {
                    /* j  is a cumulene endpoint */
                    int m;
                    pcur_atom = j;  /* pcur_atom > pnxt_atom */
                    picur2nxt = k;
                    picur_sb_parity_ord = -1;
                    for (m = 0; m < MAX_NUM_STEREO_BONDS && at[pcur_atom].sb_parity[m]; m++)
                    {
                        if (at[pcur_atom].sb_ord[m] == k)
                        {
                            picur_sb_parity_ord = m;
                            break;
                        }
                    }
                    /* djb-rwth: removing redundant code */
                }
                else
                {
                    if (chain_len2 > 1 && j2 > pnxt_atom2)
                    { /* j2 is a cumulene endpoint */
                        int m;
                        pcur_atom = j2;
                        picur2nxt = k2;
                        pnxt_atom = pnxt_atom2;
                        pinxt2cur = pinxt2cur2;
                        pinxt_sb_parity_ord = pinxt_sb_parity_ord2;
                        picur_sb_parity_ord = -1;
                        for (m = 0; m < MAX_NUM_STEREO_BONDS && at[pcur_atom].sb_parity[m]; m++)
                        {
                            if (at[pcur_atom].sb_ord[m] == k2)
                                picur_sb_parity_ord = m;
                        }
                        chain_len = chain_len2;
                        /* djb-rwth: removing redundant code */
                    }
                    else
                    {
                        chain_len = 0; /* djb-rwth: removing redundant code */
                    }
                }
                /*len += sprintf( szCurBonds + len, "%c%d", bond_char, val+1);*/
                if (chain_len)
                {
                    /* both atoms belong to a stereo bond */
                    int kc;
                    int p1 = 0, p2, p1NM = 0, p2NM, neigh, neigh1, neigh2, bHasMetal, bWellDef; /* djb-rwth: initialising p1 and p1NM */
                    int     bNeighSwitched1, bNeighSwitched2;

                    /* djb-rwth: avoiding buffer overrun as picur_sb_parity_ord == -1 is possible */
                    if (picur_sb_parity_ord >= 0)
                    {
                        p1 = SB_PARITY_1( at[pcur_atom].sb_parity[picur_sb_parity_ord] );
                        p1NM = SB_PARITY_2( at[pcur_atom].sb_parity[picur_sb_parity_ord] );
                    }

                    p2 = SB_PARITY_1( at[pnxt_atom].sb_parity[pinxt_sb_parity_ord] );
                    p2NM = SB_PARITY_2( at[pnxt_atom].sb_parity[pinxt_sb_parity_ord] );

                    bWellDef = ATOM_PARITY_WELL_DEF( p1 ) && ATOM_PARITY_WELL_DEF( p2 );
                    bHasMetal = ATOM_PARITY_WELL_DEF( p1NM ) && ATOM_PARITY_WELL_DEF( p2NM );

                    bNeighSwitched1 = bNeighSwitched2 = 0;

                    if (bWellDef || bHasMetal)
                    {

                        neigh1 = num_inp_atoms;
                        for (kc = 0; kc < at[pcur_atom].valence; kc++)
                        {
                            if (kc == picur2nxt)
                                continue;
                            neigh = at[pcur_atom].neighbor[kc];
                            if (bHasMetal && is_el_a_metal( at[neigh].el_number ))
                                continue;
                            if (neigh < neigh1)
                                neigh1 = neigh;
                        }
                        if (neigh1 < num_inp_atoms)
                        {
                            bNeighSwitched1 = ( neigh1 != at[pcur_atom].neighbor[(int) at[pcur_atom].sn_ord[picur_sb_parity_ord]] );
                        }
                        else
                        {
                            AddErrorMessage( sd->pStrErrStruct, "Cannot find 0D stereobond neighbor" );
                            /*
                            sd->nStructReadError =  99;
                            sd->nErrorType = _IS_ERROR;
                            */
                        }

                        neigh2 = num_inp_atoms;
                        for (kc = 0; kc < at[pnxt_atom].valence; kc++)
                        {
                            if (kc == pinxt2cur)
                                continue;
                            neigh = at[pnxt_atom].neighbor[kc];
                            if (bHasMetal && is_el_a_metal( at[neigh].el_number ))
                                continue;
                            if (neigh < neigh2)
                                neigh2 = neigh;
                        }
                        if (neigh2 < num_inp_atoms)
                        {
                            bNeighSwitched2 = ( neigh2 != at[pnxt_atom].neighbor[(int) at[pnxt_atom].sn_ord[pinxt_sb_parity_ord]] );
                        }
                        else
                        {
                            AddErrorMessage( sd->pStrErrStruct, "Cannot find 0D stereobond neighbor" );
                            /*
                            sd->nStructReadError =  99;
                            sd->nErrorType = _IS_ERROR;
                            */
                        }

                        if (neigh1 < num_inp_atoms && neigh2 < num_inp_atoms)
                        {
                            if (ATOM_PARITY_WELL_DEF( p1 ) && ATOM_PARITY_WELL_DEF( p2 ))
                            {
                                bond_parity = 2 - ( p1 + p2 + bNeighSwitched1 + bNeighSwitched2 ) % 2;
                            }
                            else
                            {
                                bond_parity = inchi_min( p1, p2 );
                            }

                            if (bHasMetal)
                            {
                                bond_parityNM = 2 - ( p1NM + p2NM + bNeighSwitched1 + bNeighSwitched2 ) % 2;
                            }
                            else
                            {
                                if (p1NM && p2NM)
                                {
                                    bond_parityNM = inchi_min( p1NM, p2NM );
                                }
                            }
                        }
                    }
                    else
                    {
                        if (p1 && p2)
                        {
                            bond_parity = inchi_min( p1, p2 );
                        }
                        if (p1NM && p2NM)
                        {
                            bond_parityNM = inchi_min( p1NM, p2NM );
                        }
                        if (bond_parityNM && !bond_parity)
                        {
                            bond_parity = AB_PARITY_UNDF;
                        }
                    }
                }

                len += sprintf(szCurBonds + len, "%c%s%s%d",
                    bond_char,

                    (bond_parity == AB_PARITY_ODD) ? "-" :
                    (bond_parity == AB_PARITY_EVEN) ? "+" :
                    (bond_parity == AB_PARITY_UNKN) ? "u" :
                    (bond_parity == AB_PARITY_UNDF) ? "?" : "",

                    (bond_parityNM == AB_PARITY_ODD) ? "-" :
                    (bond_parityNM == AB_PARITY_EVEN) ? "+" :
                    (bond_parityNM == AB_PARITY_UNKN) ? "u" :
                    (bond_parityNM == AB_PARITY_UNDF) ? "?" : "",

                    j2 + 1);
            }
        }
        if (len + cur_len + 2 < buf_len)
        {
            memcpy(szBuf + cur_len, szCurBonds, len);
            cur_len += len;
            szBuf[cur_len++] = ';';
            j++;
        }
        else
        {
            break;
        }
    }
    szBuf[cur_len] = '\0';
    *i = num_inp_atoms > 0 ? j : 0;

    return cur_len;
}


#define ORIG_STR_BUFLEN (7*MAXVAL+2)    /* > 7*MAXVAL+2 = 142 */

/****************************************************************************
 Fill out original input structure
 ****************************************************************************/
int OrigStruct_FillOut( CANON_GLOBALS *pCG,
                       ORIG_ATOM_DATA *orig_inp_data,
                       ORIG_STRUCT    *pOrigStruct,
                       STRUCT_DATA    *sd )
{
    char szBuf[ORIG_STR_BUFLEN];
    int  i, len, len_coord, len_atoms, len_bonds;

    pOrigStruct->polymer = NULL;
    pOrigStruct->v3000 = NULL;

    pOrigStruct->n_zy = orig_inp_data->n_zy;
    /* Coordinates */
    len_coord = i = 0;

    if (orig_inp_data->szCoord)
    {

        while ((len = WriteOrigCoord( orig_inp_data->num_inp_atoms,
            orig_inp_data->szCoord, &i, szBuf, sizeof( szBuf ) ))) /* djb-rwth: addressing LLVM warning */
        {
            len_coord += len;
        }
        pOrigStruct->szCoord = (char*) inchi_malloc( ( (long long)len_coord + 1 ) * sizeof( pOrigStruct->szCoord[0] ) ); /* djb-rwth: cast operator added */
        i = 0;
        if (pOrigStruct->szCoord &&
             len_coord == WriteOrigCoord( orig_inp_data->num_inp_atoms,
                 orig_inp_data->szCoord, &i, pOrigStruct->szCoord, len_coord + 1 ) &&
             i == orig_inp_data->num_inp_atoms)
        {
            /* success */
            if (orig_inp_data->szCoord)
            {
                inchi_free( orig_inp_data->szCoord );
                orig_inp_data->szCoord = NULL;
            }
        }
        else
        {
            return -1;
        }
    }

    /* Atoms */
    len_atoms = i = 0;
    while ((len = WriteOrigAtoms( pCG, orig_inp_data->num_inp_atoms,
        orig_inp_data->at, &i, szBuf, sizeof( szBuf ), sd ))) /* djb-rwth: addressing LLVM warning */
    {
        len_atoms += len;
        if (!orig_inp_data->num_inp_atoms)
            break;
    }
    pOrigStruct->szAtoms = (char*) inchi_malloc( ( (long long)len_atoms + 1 ) * sizeof( pOrigStruct->szAtoms[0] ) ); /* djb-rwth: cast operator added */
    i = 0;
    if (pOrigStruct->szAtoms &&
         len_atoms == WriteOrigAtoms( pCG, orig_inp_data->num_inp_atoms,
             orig_inp_data->at, &i, pOrigStruct->szAtoms, len_atoms + 1, sd ) &&
         i == orig_inp_data->num_inp_atoms)
    {
        ; /* success */
    }
    else
    {
        return -1;
    }

    /* Bonds */
    len_bonds = 0;
    i = 1;
    while ((len = WriteOrigBonds( pCG, orig_inp_data->num_inp_atoms,
#if ( FIX_CURE53_ISSUE_OOB_ALREADY_HAVE_THIS_MESSAGE==1 )
        orig_inp_data->at, &i, szBuf, sizeof(szBuf), sd))) /* djb-rwth: addressing LLVM warning */
#else
        orig_inp_data->at, &i, szBuf, sizeof(szBuf), NULL)))
#endif
    {
        len_bonds += len;
        if (!orig_inp_data->num_inp_atoms)
        {
            break;
        }
    }

    pOrigStruct->szBonds = (char*) inchi_malloc( ( (long long)len_bonds + 2 ) * sizeof( pOrigStruct->szBonds[0] ) ); /* djb-rwth: cast operator added */
    i = 1;

    if (pOrigStruct->szBonds &&
         len_bonds == WriteOrigBonds( pCG, orig_inp_data->num_inp_atoms,
             orig_inp_data->at, &i, pOrigStruct->szBonds, len_bonds + 2, sd ) &&
         i == orig_inp_data->num_inp_atoms)
    {
        ; /* success */
    }
    else
    {
        return -1;
    }
    pOrigStruct->num_atoms = orig_inp_data->num_inp_atoms;

    /* Extensions of v. 1.05 */
    if (orig_inp_data->polymer != NULL
         && orig_inp_data->polymer->n > 0
         && orig_inp_data->valid_polymer)
    {
        pOrigStruct->polymer = orig_inp_data->polymer;
                                /* pointer copy, do not free after use! */
    }
    if (orig_inp_data->v3000 != NULL)
    {
        pOrigStruct->v3000 = orig_inp_data->v3000;
                                /* pointer copy, do not free after use! */
    }

    return 0;
}


/****************************************************************************/
void OrigStruct_Free( ORIG_STRUCT *pOrigStruct )
{
    if (pOrigStruct)
    {
        if (pOrigStruct->szAtoms)
        {
            inchi_free( pOrigStruct->szAtoms );
        }
        if (pOrigStruct->szBonds)
        {
            inchi_free( pOrigStruct->szBonds );
        }
        if (pOrigStruct->szCoord)
        {
            inchi_free( pOrigStruct->szCoord );
        }

        /* For

            OAD_Polymer *polymer;
            OAD_V3000    *v3000;

            we used shallow (pointer) copy of analogs from orig_inp_data, so do not free these here */

        /*memset( pOrigStruct, 0, sizeof(*pOrigStruct) );*/
        pOrigStruct->szAtoms = NULL;
        pOrigStruct->szBonds = NULL;
        pOrigStruct->szCoord = NULL;
    }
}


/****************************************************************************
    GetSaveOptLetters

    Get the two letters encoding the saved InChI creation options.

    The first one encodes RecMet/FixedH/SUU/SLUUD options.
    Each of options is a binary switch {ON,OFF}, so it totals to 2*2*2*2=16 values
    which are encoded by capital letters 'A' through 'P'.

    The second character encodes experimental (InChI 1 extension)
    options KET and 15T.
    Each of these options is a binary switch ON/OFF, so there are 2*2=4 combinations,
    currently encoded by 'A' through 'D'.
    Note that anything but 'A' here would indicate "extended" InChI 1
    Also, there is a reservation for future needs: the 2nd memo char
    may accommodate two more ON/OFF
****************************************************************************/
void GetSaveOptLetters( unsigned char save_opt_bits, char* let1, char* let2 )
{
    const char a2p[] = "ABCDEFGHIJKLMNOP";
    /* SaveOptBits layout: {unused|unused|Ket|15T|RecMet|FixedH|SUU|SLUUD} */
    *let1 = a2p[(size_t) ( save_opt_bits & 0x0f )];
    *let2 = a2p[(size_t) ( ( save_opt_bits & 0x30 ) >> 4 )];
}


/****************************************************************************
Set line separators dependent on requested output mode
****************************************************************************/
void set_line_separators( int bINChIOutputOptions, char **pLF, char **pTAB )
{
    int  bPlainTextCommnts = 0 != ( bINChIOutputOptions & INCHI_OUT_PLAIN_TEXT_COMMENTS );

    *pLF = (char *)(bPlainTextCommnts ? "\n" : "\0");

#if ( !defined(TARGET_API_LIB) && !defined(TARGET_LIB_FOR_WINCHI) )
    {
        int  bPlainText = 0 != ( bINChIOutputOptions & ( INCHI_OUT_PLAIN_TEXT | INCHI_OUT_PLAIN_TEXT_COMMENTS ) );
        int  bPlainTabbedOutput = 0 != ( bINChIOutputOptions & INCHI_OUT_TABBED_OUTPUT ) &&
            bPlainText && !bPlainTextCommnts;

        *pTAB = bPlainTabbedOutput ? "\t" : "\n";
    }
#else
    *pTAB = "\n";
#endif

    return;
}


/****************************************************************************
Output InChI: InChI version and kind
****************************************************************************/
int OutputINCHI_VersionAndKind( INCHI_IOSTREAM   *out_file,
                                INCHI_IOS_STRING *strbuf,
                                int              bINChIOutputOptions,
                                int              is_beta,
                                char             *pLF,
                                char             *pTAB )
{
    inchi_ios_print_nodisplay( out_file, "%s%s=%s", pLF, INCHI_NAME, pLF );

    inchi_strbuf_reset( strbuf );
    inchi_strbuf_printf( strbuf, "%s", x_curr_ver );

    /* - add 'Beta' flag if applicable */
    if (is_beta)
    {
        inchi_strbuf_printf( strbuf, "B" );
    }
    /* - add 'Standard' flag if applicable */
    else if (bINChIOutputOptions & INCHI_OUT_STDINCHI)
    {
        inchi_strbuf_printf( strbuf, "S" );
    }

    inchi_ios_print_nodisplay( out_file, "%s%s", strbuf->pStr, pLF );

    return 0;
}


/****************************************************************************
Output InChI: main layer - formula, connections and hydrogens
(incl. tautomeric info == mobile H)
***************************************************************************/
int OutputINCHI_MainLayerFormula( CANON_GLOBALS    *pCG,
                                  INCHI_IOSTREAM   *out_file,
                                  INCHI_IOS_STRING *strbuf,
                                  int              num_components2[],
                                  int              *INCHI_basic_or_INCHI_reconnected,
                                  INCHI_OUT_CTL    *io,
                                  char             *pLF,
                                  char             *pTAB )
{

    /* constitution ( dot-disconnected Hill formulas: <formula> ) */

    if (num_components2[0] || num_components2[1])
    {
        szGetTag( IdentLbl, io->nTag, io->bTag1 = *INCHI_basic_or_INCHI_reconnected == INCHI_REC ? IL_REC_ : IL_FML_, io->szTag1, &io->bAlways, 1 );
        inchi_strbuf_reset( strbuf );
        io->tot_len = str_HillFormula( io->pINChISort, strbuf, &io->bOverflow, io->bOutType,
                                   io->num_components, io->bUseMulipliers );

        if (str_LineEnd( io->szTag1, &io->bOverflow, strbuf, -1, 1 ))
        {
            return 1;
        }
        if (io->n_pzz > 0 && io->n_zy > 0)
        {
            int retm = MergeZzInHillFormula(strbuf);
            if (0 != retm)
            {
                return -1;
            }
        }
        inchi_ios_print_nodisplay( out_file, "%s%s", strbuf->pStr, pLF );
    }

    LOG_NO_ARGS("\n#################### (L3318:ichiprt1.c) ##########################\n");
    LOG_MULT_ARGS("This is the Chemical formula : %s\n", strbuf->pStr);
    LOG_NO_ARGS("####################################################################\n");

    return 0;
}


/****************************************************************************/
int OutputINCHI_MainLayerConnections( CANON_GLOBALS    *pCG,
                                      INCHI_IOSTREAM   *out_file,
                                      INCHI_IOS_STRING *strbuf,
                                      int              num_components2[],
                                      int              *INCHI_basic_or_INCHI_reconnected,
                                      INCHI_OUT_CTL    *io,
                                      char             *pLF,
                                      char             *pTAB )
{
    /* connections ( semicolon/dot-disconnected connection tables ) */

    szGetTag( IdentLbl, io->nTag, io->bTag1 = IL_CONN, io->szTag1, &io->bAlways, 1 );
    inchi_strbuf_reset( strbuf );
    io->tot_len = 0;
    io->tot_len2 = str_Connections( pCG, io->pINChISort, strbuf, &io->bOverflow, io->bOutType,
                                    io->ATOM_MODE, io->num_components, io->bUseMulipliers );

    /* current version does not output empty (";;;;") connectivity */

    if (io->tot_len != io->tot_len2)
    { /* 2004-06-30: never output empty connection table */
        io->tot_len = io->tot_len2;
        if (str_LineEnd( io->szTag1, &io->bOverflow, strbuf, -2, io->bPlainTextTags ))
        {
            return 1; /* pStr overfow */
        }
        inchi_ios_print_nodisplay( out_file, "%s%s", strbuf->pStr, pLF );
    }

    LOG_NO_ARGS("\n##################### (L3357:ichiprt1.c) #########################\n");
    LOG_MULT_ARGS("This is the Connection Layer : %s\n", strbuf->pStr);
    LOG_NO_ARGS("####################################################################\n");

    return 0;
}


/****************************************************************************/
int OutputINCHI_MainLayerHydrogens( CANON_GLOBALS    *pCG,
                                    INCHI_IOSTREAM   *out_file,
                                    INCHI_IOS_STRING *strbuf,
                                    int              num_components2[],
                                    int              *INCHI_basic_or_INCHI_reconnected,
                                    INCHI_OUT_CTL    *io,
                                    char             *pLF,
                                    char             *pTAB )
{

    /* hydrogen atoms (do not output empty) */

    if (INCHI_SEGM_FILL == INChI_SegmentAction( io->sDifSegs[io->nCurINChISegment][DIFS_h_H_ATOMS] ))
    {
        szGetTag( IdentLbl, io->nTag, io->bTag1 = IL_ALLH, io->szTag1, &io->bAlways, 1 );
        inchi_strbuf_reset( strbuf );
        io->tot_len = 0;
        io->tot_len2 = str_H_atoms( io->pINChISort, strbuf, &io->bOverflow, io->bOutType,
                                io->ATOM_MODE, io->TAUT_MODE,
                                io->num_components, io->bUseMulipliers );
        if (io->tot_len != io->tot_len2)
        { /* 2004-06-21: never output empty */
            io->tot_len = io->tot_len2;
            if (str_LineEnd( io->szTag1, &io->bOverflow, strbuf, -2, 1 ))
            {
                return 1;
            }
            inchi_ios_print_nodisplay( out_file, "%s%s", strbuf->pStr, pLF );
        }
    }

    LOG_NO_ARGS("\n###################### (L3396:ichiprt1.c) ########################\n");
    LOG_MULT_ARGS("This is the Hydrogen Layer : %s\n", strbuf->pStr);
    LOG_NO_ARGS("####################################################################\n");

    return 0;
}


/****************************************************************************
Output InChI: charge and  removed protons layers
****************************************************************************/
int OutputINCHI_ChargeAndRemovedAddedProtonsLayers( CANON_GLOBALS    *pCG,
                                                    INCHI_IOSTREAM   *out_file,
                                                    INCHI_IOS_STRING *strbuf,
                                                    INCHI_OUT_CTL    *io,
                                                    char             *pLF,
                                                    char             *pTAB )
{

    /* charge  */

    io->nSegmAction = INChI_SegmentAction( io->sDifSegs[io->nCurINChISegment][DIFS_q_CHARGE] );
    if (io->nSegmAction)
    {
        szGetTag( IdentLbl, io->nTag, io->bTag1 = IL_CHRG | io->bFhTag, io->szTag1, &io->bAlways, 1 );
        inchi_strbuf_reset( strbuf );
        io->tot_len = 0;
        if (INCHI_SEGM_FILL == io->nSegmAction)
        {
            io->tot_len = str_Charge2( io->pINChISort, io->pINChISort2,
                                   strbuf, &io->bOverflow, io->bOutType, io->num_components,
                                   io->bSecondNonTautPass, io->bOmitRepetitions, io->bUseMulipliers );
            io->bNonTautNonIsoIdentifierNotEmpty += io->bSecondNonTautPass;
        }
        if (str_LineEnd( io->szTag1, &io->bOverflow, strbuf, -io->nSegmAction, io->bPlainTextTags ))
        {
            return 1;
        }
        inchi_ios_print_nodisplay( out_file, "%s%s", strbuf->pStr, pLF );
    }

    /* removed protons */

    if (io->iCurTautMode == TAUT_YES && !io->bSecondNonTautPass)
    {

        io->nSegmAction = INChI_SegmentAction( io->sDifSegs[io->nCurINChISegment][DIFS_p_PROTONS] );
        if (io->nSegmAction)
        {
            szGetTag( IdentLbl, io->nTag, io->bTag1 = IL_PROT | io->bFhTag, io->szTag1, &io->bAlways, 1 );
            inchi_strbuf_reset( strbuf );
            io->tot_len = 0;
            inchi_strbuf_printf( strbuf, "%+d", io->nNumRemovedProtons );
            if (str_LineEnd( io->szTag1, &io->bOverflow, strbuf, -io->nSegmAction, io->bPlainTextTags ))
            {
                return 1;
            }
            inchi_ios_print_nodisplay( out_file, "%s%s", strbuf->pStr, pLF );
        }
        else
        {
            if (io->bPlainTextTags == 1) inchi_ios_print_nodisplay( out_file, "/" );
        }
    }

    return 0;
}


/****************************************************************************
Output InChI: stereo layer with sublayers
****************************************************************************/
int OutputINCHI_StereoLayer( CANON_GLOBALS    *pCG,
                             INCHI_IOSTREAM   *out_file,
                             INCHI_IOS_STRING *strbuf,
                             INCHI_OUT_CTL    *io,
                             char             *pLF,
                             char             *pTAB )
{

    {
        int i;
        i = INChI_SegmentAction( io->sDifSegs[io->nCurINChISegment][DIFS_t_SATOMS] ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
        /* djb-rwth: removing redundant code */
    }

    if (INChI_SegmentAction( io->sDifSegs[io->nCurINChISegment][DIFS_b_SBONDS] ) ||
         INChI_SegmentAction( io->sDifSegs[io->nCurINChISegment][DIFS_t_SATOMS] ) ||
         INChI_SegmentAction( io->sDifSegs[io->nCurINChISegment][DIFS_m_SP3INV] ) ||
         INChI_SegmentAction( io->sDifSegs[io->nCurINChISegment][DIFS_s_STYPE] ))
    {

        /*  stereo */

        szGetTag( IdentLbl, io->nTag, io->bTag1 = IL_STER | io->bFhTag, io->szTag1, &io->bAlways, 1 );

        /*  sp2 */

        /*if ( bStereoSp2[io->iCurTautMode]  )*/
        if ((io->nSegmAction = INChI_SegmentAction( io->sDifSegs[io->nCurINChISegment][DIFS_b_SBONDS] ))) /* djb-rwth: addressing LLVM warning */
        {
            szGetTag( IdentLbl, io->nTag, io->bTag2 = io->bTag1 | IL_DBND, io->szTag2, &io->bAlways, 1 );
            inchi_strbuf_reset( strbuf );
            io->tot_len = 0;
            if (INCHI_SEGM_FILL == io->nSegmAction)
            {
                io->tot_len = str_Sp2( io->pINChISort, io->pINChISort2, strbuf, &io->bOverflow,
                                        io->bOutType, io->TAUT_MODE, io->num_components,
                                        io->bSecondNonTautPass, io->bOmitRepetitions, io->bUseMulipliers );

                io->bNonTautNonIsoIdentifierNotEmpty += io->bSecondNonTautPass;
            }

            if (str_LineEnd( io->szTag2, &io->bOverflow, strbuf, -io->nSegmAction, io->bPlainTextTags ))
            {
                return 1;
            }
            inchi_ios_print_nodisplay( out_file, "%s%s", strbuf->pStr, pLF );
        }
        else
        {
            if (io->bPlainTextTags == 1)
            {
                inchi_ios_print_nodisplay( out_file, "/" ); /* sp2 */
            }
        }

        /*  sp3 */

        /*if ( bStereoSp3[io->iCurTautMode]  )*/
        if ((io->nSegmAction = INChI_SegmentAction( io->sDifSegs[io->nCurINChISegment][DIFS_t_SATOMS] ))) /* djb-rwth: addressing LLVM warning */
        {
            io->bRelRac = io->bRelativeStereo[io->iCurTautMode] || io->bRacemicStereo[io->iCurTautMode];
            szGetTag( IdentLbl, io->nTag, io->bTag2 = io->bTag1 | IL_SP3S, io->szTag2, &io->bAlways, 1 );
            inchi_strbuf_reset( strbuf );
            io->tot_len = 0;
            if (INCHI_SEGM_FILL == io->nSegmAction)
            {
                io->tot_len = str_Sp3( io->pINChISort, io->pINChISort2, strbuf, &io->bOverflow,
                                       io->bOutType, io->TAUT_MODE, io->num_components, io->bRelRac,
                                   io->bSecondNonTautPass, io->bOmitRepetitions, io->bUseMulipliers );

                io->bNonTautNonIsoIdentifierNotEmpty += io->bSecondNonTautPass;
            }

            if (str_LineEnd( io->szTag2, &io->bOverflow, strbuf, -io->nSegmAction, io->bPlainTextTags ))
            {
                return 2;
            }
            inchi_ios_print_nodisplay( out_file, "%s%s", strbuf->pStr, pLF );
        }
        else
        {
            if (io->bPlainTextTags == 1) inchi_ios_print_nodisplay( out_file, "/" ); /* sp3 */
        }

        /* bStereoAbsInverted[io->iCurTautMode]  */

        /* if ( bStereoAbs[io->iCurTautMode]  ) */
        if ((io->nSegmAction = INChI_SegmentAction( io->sDifSegs[io->nCurINChISegment][DIFS_m_SP3INV] ))) /* djb-rwth: addressing LLVM warning */
        {
            szGetTag( IdentLbl, io->nTag, io->bTag2 = io->bTag1 | IL_INVS, io->szTag2, &io->bAlways, 1 );
            inchi_strbuf_reset( strbuf ); io->tot_len = 0;
            if (INCHI_SEGM_FILL == io->nSegmAction)
            {
                io->tot_len = str_StereoAbsInv( io->pINChISort, strbuf,
                                            &io->bOverflow, io->bOutType, io->num_components );
                io->bNonTautNonIsoIdentifierNotEmpty += io->bSecondNonTautPass;
            }

            if (str_LineEnd( io->szTag2, &io->bOverflow, strbuf, -io->nSegmAction, io->bPlainTextTags ))
            {
                return 3;
            }
            inchi_ios_print_nodisplay( out_file, "%s%s", strbuf->pStr, pLF );
        }
        else
        {
            if (io->bPlainTextTags == 1)
            {
                inchi_ios_print_nodisplay( out_file, "/" ); /* stereo-abs-inv */
            }
        }

        /* stereo type */

        /*if ( io->bRacemicStereo[io->iCurTautMode] || io->bRelativeStereo[io->iCurTautMode] || bStereoAbs[io->iCurTautMode] )*/
        if ((io->nSegmAction = INChI_SegmentAction( io->sDifSegs[io->nCurINChISegment][DIFS_s_STYPE] ))) /* djb-rwth: addressing LLVM warning */
        {
            const char *p_stereo = io->bRelativeStereo[io->iCurTautMode] ? x_rel :
                io->bRacemicStereo[io->iCurTautMode] ? x_rac : x_abs;
            szGetTag( IdentLbl, io->nTag, io->bTag2 = io->bTag1 | IL_TYPS, io->szTag2, &io->bAlways, 1 );
            inchi_strbuf_reset( strbuf ); io->tot_len = 0;
            if (INCHI_SEGM_FILL == io->nSegmAction)
            {
                ( io->tot_len ) += MakeDelim( p_stereo, strbuf, &io->bOverflow );
                io->bNonTautNonIsoIdentifierNotEmpty += io->bSecondNonTautPass;
            }
            if (str_LineEnd( io->szTag2, &io->bOverflow, strbuf, -io->nSegmAction, io->bPlainTextTags ))
            {
                return 1;
            }
            inchi_ios_print_nodisplay( out_file, "%s%s", strbuf->pStr, pLF );
        }
        if (io->bPlainTextTags == 1)
        {
            inchi_ios_print_nodisplay( out_file, "/" );  /* no abs, inv or racemic stereo */
        }
    }
    else
    {
        if (io->bPlainTextTags == 1) inchi_ios_print_nodisplay( out_file, "////" ); /* sp3, sp2, abs-inv, stereo.type */
    }

    return 0;
}


/****************************************************************************
Output InChI: isotopic layer and sublayers  ****************************************************************************/
int OutputINCHI_IsotopicLayer( CANON_GLOBALS    *pCG,
                               INCHI_IOSTREAM   *out_file,
                               INCHI_IOS_STRING *strbuf,
                               int              *INCHI_basic_or_INCHI_reconnected,
                               INCHI_OUT_CTL    *io,
                               char             *pLF,
                               char             *pTAB )
{

    if (INChI_SegmentAction( io->sDifSegs[io->nCurINChISegment][DIFS_i_IATOMS] ))
    {
        /*  isotopic #1:  composition -- atoms -- do not output in xml if empty */
        szGetTag( IdentLbl, io->nTag, io->bTag1 = IL_ISOT | io->bFhTag, io->szTag1, &io->bAlways, 1 );
        /* isotopic atoms without mobile H.
         * Fixed 2004-06-15: always output if not bXml. Note:
         * Previous condition if( bHasIsotopicAtoms[io->iCurTautMode] || bIsotopic && !bXml)
         * did not optput /i in case of only mobile isotopic H
         */
        if ((io->nSegmAction = INChI_SegmentAction( io->sDifSegs[io->nCurINChISegment][DIFS_i_IATOMS] ))) /* djb-rwth: addressing LLVM warning */
        {
            szGetTag( IdentLbl, io->nTag, io->bTag2 = io->bTag1 | IL_ATMS, io->szTag2, &io->bAlways, 1 );
            inchi_strbuf_reset( strbuf );
            io->tot_len = 0;
            /*if ( bHasIsotopicAtoms[io->iCurTautMode] )*/
            if (INCHI_SEGM_FILL == io->nSegmAction)
            {
                io->tot_len2 = str_IsoAtoms( io->pINChISort, io->pINChISort2, strbuf, &io->bOverflow,
                                             io->bOutType, io->TAUT_MODE, io->num_components, io->bAbcNumbers,
                                             io->bSecondNonTautPass, io->bOmitRepetitions, io->bUseMulipliers );
                io->bNonTautIsoIdentifierNotEmpty += io->bSecondNonTautPass;
            }
            else
            {
                io->tot_len2 = io->tot_len;
            }

            io->tot_len = io->tot_len2;
            if (str_LineEnd( io->szTag2, &io->bOverflow, strbuf, -io->nSegmAction, io->bPlainTextTags ))
            {
                return 1;
            }
            inchi_ios_print_nodisplay( out_file, "%s%s", strbuf->pStr, pLF );
        }

        /*  isotopic #1a:  composition -- exchangeable isotopic H (mobile H only) */
        /*if ( !io->bSecondNonTautPass && bHasIsoH )*/
        if ((io->nSegmAction = INChI_SegmentAction( io->sDifSegs[io->nCurINChISegment][DIFS_h_H_ATOMS] ))) /* djb-rwth: addressing LLVM warning */
        {
            szGetTag( IdentLbl, io->nTag, io->bTag2 = io->bTag1 | IL_XCGA, io->szTag2, &io->bAlways, 1 );
            inchi_strbuf_reset( strbuf );
            io->tot_len = 0;
            ( io->tot_len ) += MakeIsoHString( io->num_iso_H, strbuf, io->TAUT_MODE, &io->bOverflow );
            io->bNonTautIsoIdentifierNotEmpty += io->bSecondNonTautPass;
            if (str_LineEnd( io->szTag2, &io->bOverflow, strbuf, -io->nSegmAction, io->bPlainTextTags ))
            {
                return 2;
            }
            inchi_ios_print_nodisplay( out_file, "%s%s", strbuf->pStr, pLF );
        }

        /***************************************************
         *
         *       Isotopic stereo
         *
         ***************************************************/

        /*if ( bIsotopicStereo[io->iCurTautMode] )*/
        if (INChI_SegmentAction( io->sDifSegs[io->nCurINChISegment][DIFS_b_SBONDS] ) ||
             INChI_SegmentAction( io->sDifSegs[io->nCurINChISegment][DIFS_t_SATOMS] ) ||
             INChI_SegmentAction( io->sDifSegs[io->nCurINChISegment][DIFS_m_SP3INV] ) ||
             INChI_SegmentAction( io->sDifSegs[io->nCurINChISegment][DIFS_s_STYPE] ))
        {
            /*  stereo */
            szGetTag( IdentLbl, io->nTag, io->bTag2 = io->bTag1 | IL_STER, io->szTag2, &io->bAlways, 1 );

            /************************
              isotopic #2:  sp2
             ************************/
            /*if ( bIsotopicStereoSp2[io->iCurTautMode]  )*/
            if ((io->nSegmAction = INChI_SegmentAction( io->sDifSegs[io->nCurINChISegment][DIFS_b_SBONDS] ))) /* djb-rwth: addressing LLVM warning */
            {
                szGetTag( IdentLbl, io->nTag, io->bTag3 = io->bTag2 | IL_DBND, io->szTag3, &io->bAlways, 1 );
                inchi_strbuf_reset( strbuf );
                io->tot_len = 0;
                if (INCHI_SEGM_FILL == io->nSegmAction)
                {
                    io->tot_len = str_IsoSp2( io->pINChISort, io->pINChISort2, strbuf, &io->bOverflow,
                                              io->bOutType, io->TAUT_MODE, io->num_components,
                                          io->bSecondNonTautPass, io->bOmitRepetitions, io->bUseMulipliers );
                    io->bNonTautIsoIdentifierNotEmpty += io->bSecondNonTautPass;
                }
                if (str_LineEnd( io->szTag3, &io->bOverflow, strbuf, -io->nSegmAction, io->bPlainTextTags ))
                {
                    return 3;
                }
                inchi_ios_print_nodisplay( out_file, "%s%s", strbuf->pStr, pLF );
            }
            else
            {
                if (io->bPlainTextTags == 1) inchi_ios_print_nodisplay( out_file, "/" ); /* iso sp2 */
            }

            /************************
              isotopic #3:  sp3
             ************************/
            /*if ( bIsotopicStereoSp3[io->iCurTautMode]  )*/
            if ((io->nSegmAction = INChI_SegmentAction( io->sDifSegs[io->nCurINChISegment][DIFS_t_SATOMS] ))) /* djb-rwth: addressing LLVM warning */
            {
                io->bRelRac = io->bIsotopicRelativeStereo[io->iCurTautMode] || io->bIsotopicRacemicStereo[io->iCurTautMode];

                szGetTag( IdentLbl, io->nTag, io->bTag3 = io->bTag2 | IL_SP3S, io->szTag3, &io->bAlways, 1 );
                inchi_strbuf_reset( strbuf );
                io->tot_len = 0;
                if (INCHI_SEGM_FILL == io->nSegmAction)
                {
                    io->tot_len = str_IsoSp3( io->pINChISort, io->pINChISort2, strbuf, &io->bOverflow,
                                              io->bOutType, io->TAUT_MODE, io->num_components, io->bRelRac,
                                              io->bSecondNonTautPass, io->bOmitRepetitions, io->bUseMulipliers );
                    io->bNonTautIsoIdentifierNotEmpty += io->bSecondNonTautPass;
                }
                if (str_LineEnd( io->szTag3, &io->bOverflow, strbuf, -io->nSegmAction, io->bPlainTextTags ))
                {
                    return 5;
                }
                inchi_ios_print_nodisplay( out_file, "%s%s", strbuf->pStr, pLF );
            }
            else
            {
                if (io->bPlainTextTags == 1)
                {
                    inchi_ios_print_nodisplay( out_file, "/" ); /* iso-sp3 */
                }
            }

            /* isotopic #4: abs inverted */
            if ((io->nSegmAction = INChI_SegmentAction( io->sDifSegs[io->nCurINChISegment][DIFS_m_SP3INV] ))) /* djb-rwth: addressing LLVM warning */
            {
                szGetTag( IdentLbl, io->nTag, io->bTag3 = io->bTag2 | IL_INVS, io->szTag3, &io->bAlways, 1 );
                inchi_strbuf_reset( strbuf );
                io->tot_len = 0;
                if (INCHI_SEGM_FILL == io->nSegmAction)
                {
                    io->tot_len = str_IsoStereoAbsInv( io->pINChISort, strbuf,
                                                   &io->bOverflow, io->bOutType, io->num_components );
                    io->bNonTautIsoIdentifierNotEmpty += io->bSecondNonTautPass;
                }
                if (str_LineEnd( io->szTag3, &io->bOverflow, strbuf, -io->nSegmAction, io->bPlainTextTags ))
                {
                    return 5;
                }
                inchi_ios_print_nodisplay( out_file, "%s%s", strbuf->pStr, pLF );
            }
            else
            {
                if (io->bPlainTextTags == 1)
                {
                    inchi_ios_print_nodisplay( out_file, "/" );
                }
            }

            /* isotopic #5: stereo type. Do not output if it has already been output in non-iso */
            if ((io->nSegmAction = INChI_SegmentAction( io->sDifSegs[io->nCurINChISegment][DIFS_s_STYPE] ))) /* djb-rwth: addressing LLVM warning */
            {
                const char *p_stereo = io->bIsotopicRelativeStereo[io->iCurTautMode] ? x_rel :
                    io->bIsotopicRacemicStereo[io->iCurTautMode] ? x_rac : x_abs;
                szGetTag( IdentLbl, io->nTag, io->bTag3 = io->bTag2 | IL_TYPS, io->szTag3, &io->bAlways, 1 );
                inchi_strbuf_reset( strbuf );
                io->tot_len = 0;
                if (INCHI_SEGM_FILL == io->nSegmAction)
                {
                    io->tot_len += MakeDelim( p_stereo, strbuf, &io->bOverflow );
                    io->bNonTautIsoIdentifierNotEmpty += io->bSecondNonTautPass;
                }
                if (str_LineEnd( io->szTag3, &io->bOverflow, strbuf, -io->nSegmAction, io->bPlainTextTags ))
                {
                    return 6;
                }
                inchi_ios_print_nodisplay( out_file, "%s%s", strbuf->pStr, pLF );
            }
            if (io->bPlainTextTags == 1)
            {
                inchi_ios_print_nodisplay( out_file, "/" );  /* no abs, inv or racemic stereo */
            }
        }
        else
        {
            /* no isotopic stereo */
            if (io->bPlainTextTags == 1)
            {
                inchi_ios_print_nodisplay( out_file, "////" ); /* sp3, sp2, abs-inv, stereo.type */
            }
        }
    }
    else
    {
        if (io->bPlainTextTags == 1)
        {
            inchi_ios_print_nodisplay( out_file, "///" ); /* isotopic composition, sp2, sp3 */
        }
        if (io->bPlainTextTags == 1)
        {
            inchi_ios_print_nodisplay( out_file, "//" );   /* inv or racemic stereo */
        }
    }

#if ( CANON_FIXH_TRANS == 1 )
    if (io->bOutType == OUT_NONTAUT && io->bOutputType == OUT_TN && io->bSecondNonTautPass &&
         INCHI_SEGM_FILL == INChI_SegmentAction( io->sDifSegs[DIFL_F][DIFS_o_TRANSP] ))
    {
        /* find and print non-tautomeric components transposition, if non-trivial */
        AT_NUMB *nTrans_n, *nTrans_s;

        if (0 < bin_AuxTautTrans( io->pINChISort, io->pINChISort2, &nTrans_n, &nTrans_s, io->bOutType, io->num_components ))
        {
            /* a non-trivial transposition does exist; output start tag */
            szGetTag( IdentLbl, io->nTag, io->bTag1 = IL_TRNS | io->bFhTag, io->szTag1, &io->bAlways, 1 );
            inchi_strbuf_reset( strbuf );
            io->tot_len = 0;
            /* print the transposition, cycle after cycle */
            io->tot_len = str_AuxTautTrans( pCG, nTrans_n, nTrans_s, strbuf,
                                            &io->bOverflow, io->TAUT_MODE, io->num_components );
            io->bNonTautIsoIdentifierNotEmpty += io->bSecondNonTautPass;
            if (str_LineEnd( io->szTag1, &io->bOverflow, strbuf, -1, io->bPlainTextTags ))
            {
                return 7;
            }
            inchi_ios_print_nodisplay( out_file, "%s%s", strbuf->pStr, pLF );
             /* detected transposition */
            ( *io->pSortPrintINChIFlags ) |=
                ( *INCHI_basic_or_INCHI_reconnected == INCHI_BAS ) ? FLAG_SORT_PRINT_TRANSPOS_BAS : FLAG_SORT_PRINT_TRANSPOS_REC;
        }
        else
        {
            if (io->bPlainTextTags == 1)
            {
                inchi_ios_print_nodisplay( out_file, "/" );
            }
        }
    }
#endif

    return 0;
}


/****************************************************************************
Output InChI: FixedH layer and related sublayers
****************************************************************************/
int OutputINCHI_FixedHLayerWithSublayers( CANON_GLOBALS    *pCG,
                                          INCHI_IOSTREAM   *out_file,
                                          INCHI_IOS_STRING *strbuf,
                                          int              *INCHI_basic_or_INCHI_reconnected,
                                          INCHI_OUT_CTL    *io,
                                          char             *pLF,
                                          char             *pTAB,
                                          int              *then_goto_repeat )
{

    *then_goto_repeat = 0;

    if (io->bOutType == OUT_TN &&
         !( io->bSecondNonTautPass ) &&
         io->bNonTautIsIdenticalToTaut &&
         io->bTautomeric &&
         io->bNonTautomeric)
    {
            /* Fixed-H layer is empty in the Identifier */
        ( *io->pSortPrintINChIFlags ) |=
            ( *INCHI_basic_or_INCHI_reconnected == INCHI_BAS ) ? FLAG_SORT_PRINT_NO_NFIX_H_BAS :
            FLAG_SORT_PRINT_NO_NFIX_H_REC;
        ( *io->pSortPrintINChIFlags ) |=
            ( *INCHI_basic_or_INCHI_reconnected == INCHI_BAS ) ? FLAG_SORT_PRINT_NO_IFIX_H_BAS :
            FLAG_SORT_PRINT_NO_IFIX_H_REC;
    }

    if (io->bOutType == OUT_TN &&
         !io->bNonTautIsIdenticalToTaut && /* added 2004-10-04 Fix16 */
#ifdef OLD_ITEM_DISCOVERY
         io->bTautomeric &&
         io->bNonTautomeric &&
#endif
         INChI_SegmentAction( io->sDifSegs[DIFL_F][DIFS_f_FORMULA] )
                                    /* special case: removed isolated H(+): */
                                    /* || io->iCurTautMode == TAUT_YES && num_comp[TAUT_YES] < num_comp[TAUT_NON] &&
                                        0 < num_comp[TAUT_NON]*/
       )

    {
        /* add the second (non-tautomeric) output */
        io->bOutType = OUT_NONTAUT;    /* pick up only non-tautomeric representation of tautomeric */
        io->iCurTautMode = TAUT_NON;
        io->pINChISort = io->pINChISortTautAndNonTaut[TAUT_NON];
        io->bSecondNonTautPass = 1;
        io->nCurINChISegment = DIFL_F;
        io->num_components = io->num_comp[io->iCurTautMode]; /* number of components could change due to removal of isolated H(+) from tautomeric */
        io->bFhTag = IL_FIXH;
        szGetTag( IdentLbl, io->nTag, io->bTag1 = io->bFhTag, io->szTag1, &io->bAlways, 1 );
        /***** constitution non-taut: dot-disconnected Hill formulas: <formula> -- only if different */
        szGetTag( IdentLbl, io->nTag, io->bTag1 = IL_FMLF | io->bFhTag, io->szTag1, &io->bAlways, 1 );
        inchi_strbuf_reset( strbuf ); io->tot_len = 0;
        io->nSegmAction = INChI_SegmentAction( io->sDifSegs[io->nCurINChISegment][DIFS_f_FORMULA] );
        if (INCHI_SEGM_FILL == io->nSegmAction)
        {
            io->tot_len2 = str_HillFormula2( io->pINChISort, io->pINChISort2,
                                             strbuf, &io->bOverflow, io->bOutType,
                                             io->num_components, io->bUseMulipliers );
            if (io->n_pzz > 0 && io->n_zy > 0)
            {
                MergeZzInHillFormula(strbuf);
            }
            io->bNonTautNonIsoIdentifierNotEmpty += io->bSecondNonTautPass;
        }
        else
        {
            io->tot_len2 = io->tot_len;
        }
        io->tot_len = io->tot_len2;
        if (str_LineEnd( io->szTag1, &io->bOverflow, strbuf, -io->nSegmAction, io->bPlainTextTags ))
        {
            return 1;
        }
        inchi_ios_print_nodisplay( out_file, "%s%s", strbuf->pStr, pLF );

        io->nSegmAction = INChI_SegmentAction( io->sDifSegs[io->nCurINChISegment][DIFS_h_H_ATOMS] );

        if (INCHI_SEGM_FILL == io->nSegmAction)
        {
            szGetTag( IdentLbl, io->nTag, io->bTag1 = IL_HFIX | io->bFhTag, io->szTag1, &io->bAlways, 1 );
            inchi_strbuf_reset( strbuf ); io->tot_len = 0; /* open H-fixed */
            /* output the second non-tautomeric item: fixed H -- do not output in xml if empty */
            io->tot_len2 = str_FixedH_atoms( io->pINChISort, strbuf,
                                             &io->bOverflow, io->bOutType, io->ATOM_MODE,
                                             io->num_components, io->bUseMulipliers );
            io->tot_len = io->tot_len2;
            if (str_LineEnd( io->szTag1, &io->bOverflow, strbuf, -io->nSegmAction, io->bPlainTextTags ))
            {
                return 2;
            }
            inchi_ios_print_nodisplay( out_file, "%s%s", strbuf->pStr, pLF );
            io->bNonTautNonIsoIdentifierNotEmpty += io->bSecondNonTautPass;
        }
        *then_goto_repeat = 1;
        return 0;
    }

    else
    {
        if (io->bOutType == OUT_NONTAUT && io->bOutputType == OUT_TN && io->bSecondNonTautPass /* && io->bTautomeric && io->bNonTautomeric*/)
        {
            /* the second (non-taut) output has been done; restore variables */
            io->bOutType = OUT_TN;
            io->iCurTautMode = TAUT_YES;
            io->pINChISort = io->pINChISortTautAndNonTaut[TAUT_YES];
            io->bSecondNonTautPass = 0;
            io->num_components = io->num_comp[io->iCurTautMode];
            if (!io->bNonTautNonIsoIdentifierNotEmpty)
            {
                /* Fixed-H layer is empty in the Identifier */
                ( *io->pSortPrintINChIFlags ) |= ( *INCHI_basic_or_INCHI_reconnected == INCHI_BAS ) ? FLAG_SORT_PRINT_NO_NFIX_H_BAS :
                    FLAG_SORT_PRINT_NO_NFIX_H_REC;
            }
            if (!io->bNonTautIsoIdentifierNotEmpty)
            {
                /* Fixed-H layer is empty in the Identifier */
                ( *io->pSortPrintINChIFlags ) |= ( *INCHI_basic_or_INCHI_reconnected == INCHI_BAS ) ? FLAG_SORT_PRINT_NO_IFIX_H_BAS :
                    FLAG_SORT_PRINT_NO_IFIX_H_REC;
            }
            io->bFhTag = 0;
        }
    }

    return 0;
}


/****************************************************************************
Output InChI: polymer layer
****************************************************************************/
static int OutputINCHI_PolymerLayer( CANON_GLOBALS *pCG,
                                     INCHI_IOSTREAM *out_file,
                                     INCHI_IOS_STRING *strbuf,
                                     int *INCHI_basic_or_INCHI_reconnected,
                                     ORIG_ATOM_DATA *orig_inp_data,
                                     ORIG_STRUCT *pOrigStruct,
                                     INCHI_OUT_CTL *io,
                                     char *pLF,
                                     char *pTAB )
{
    int i, err = 0;
    int nunits2 = 0;
    int n_used_stars = 0;
    int *cano_nums = NULL, *compnt_nums = NULL, *unum = NULL, *old_stars = NULL;
    OAD_PolymerUnit *u = NULL;
    OAD_PolymerUnit **units2 = NULL;
    OAD_Polymer *p = NULL;
    OAD_AtProps *aprops = NULL;
    int nat,num_inp_bonds;
    inp_ATOM    *at = NULL;
    int is_inchi2inchi = 0;

    if (!orig_inp_data)
    {
        goto exit_function;
    }

    at = orig_inp_data->at;
    nat = orig_inp_data->num_inp_atoms;
    num_inp_bonds = orig_inp_data->num_inp_bonds;

    
    if (pOrigStruct && !pOrigStruct->polymer)
    {
        return 0;
    }

    if (pOrigStruct)
    {
        p = pOrigStruct->polymer;
        is_inchi2inchi = !pOrigStruct->szAtoms && !pOrigStruct->szBonds && !pOrigStruct->szCoord;

        if (is_inchi2inchi)
        {
            err = NOT_YET_I2I_FOR_POLYMERS;
            goto exit_function;
        }

        /*OAD_Polymer_DebugTrace( p );*/

        /* Get canonical numbers and numbers-of-components for each original atom */
        cano_nums = (int*)inchi_calloc((long long)pOrigStruct->num_atoms + 1, sizeof(int)); /* djb-rwth: cast operator added */
        if (!cano_nums)
        {
            err = 1;
            goto exit_function;
        }
        compnt_nums = (int*)inchi_calloc((long long)pOrigStruct->num_atoms + 1, sizeof(int)); /* djb-rwth: cast operator added */
        if (!compnt_nums)
        {
            err = 2;
            goto exit_function;
        }
        err = InternallyGetCanoNumsAndComponentNums(pCG,
            strbuf,
            io,
            pOrigStruct->num_atoms,
            cano_nums,
            compnt_nums);
        if (err != 0)
        {
            err = 3;
            goto exit_function;
        }


        /* Set atom properties for sorting */
        aprops = (OAD_AtProps*)inchi_calloc((long long)nat + 1, sizeof(OAD_AtProps)); /* djb-rwth: cast operator added */
        /* nat + 1: add extra element for possibe 1-based indexing */
        if (!aprops)
        {
            /* djb-rwth: avoiding memory leak */
            if (cano_nums)
            {
                inchi_free(cano_nums);
            }
            if (compnt_nums)
            {
                inchi_free(compnt_nums);
            }
            if (aprops)
            {
                inchi_free(aprops);
            }
            return 0;
        }

        /* Note that aprops[] is in orig_atoms domain (0-based) and      */
        /* u (from units) are in cano_nums domain (1-based)              */
        /* Supply non-NULL cano_nums to adjust the domains (base will be adjusted at place) */
        OAD_Polymer_SetAtProps(p, at, nat, &num_inp_bonds, aprops, cano_nums);


        /* Make a working copy of polymer units data: units2 is a copy        */
        /* of original polymer units (p->units) with atomic numbers changed    */
        /* to curr canonical ones; atoms in alists sorted; atoms in blists    */
        /* and blists themselves sorted                                     */
        units2 = (OAD_PolymerUnit**)inchi_calloc(p->n, sizeof(OAD_PolymerUnit*));

        if (NULL == units2)
        {
            err = 3;
            goto exit_function;
        }
        memset(units2, 0, sizeof(*units2)); /* djb-rwth: memset_s C11/Annex K variant? */

        old_stars = (int*)inchi_calloc(pOrigStruct->polymer->n_pzz, sizeof(int));
        if (NULL == old_stars)
        {
            err = 3;
            goto exit_function;
        }
        for (i = 0; i < pOrigStruct->polymer->n_pzz; i++)
        {
            old_stars[i] = pOrigStruct->polymer->pzz[i];
        }


        for (i = 0; i < p->n; i++)
        {
            units2[i] = OAD_PolymerUnit_CreateCopy(p->units[i]);
            if (NULL == units2[i]) /* djb-rwth: unresolved issue -- revision required? -- units2 properly allocated, and loop index well defined */
            {
                err = 4;
                goto exit_function;
            }
            nunits2 = i + 1;
        }

        /* unum contains numbers of units (0..p->n) as they go  */
        /* when sorted by alist's in lexicographic order        */
        unum = (int*)inchi_calloc(p->n, sizeof(int));
        if (NULL == unum)
        {
            err = 4;
            goto exit_function;
        }


        err = OAD_Polymer_PrepareWorkingSet(p, cano_nums, compnt_nums, units2, unum);

        if (err != 0)
        {
            err = 5;
            goto exit_function;
        }

        /* Prepare polymer substring */

        /* Mark layer beginning */
        inchi_strbuf_printf(strbuf, "%s", "/z");

        /* Print polymer units data */
        n_used_stars = 0;
        for (i = 0; i < p->n; i++)
        {
            /* For each unit u ... */
            u = units2[unum[i]];
            /* djb-rwth: addressing coverity ID #499574 -- all NULL checks already done above */
            err = OutputINCHI_PolymerLayer_SingleUnit(u,
                io->bPolymers,
                pOrigStruct->polymer->n_pzz,
                &n_used_stars, aprops,
                cano_nums,
                orig_inp_data,
                pOrigStruct, strbuf);
            if (err)
            {
                goto exit_function;
            }
            if (i < p->n - 1)
            {
                inchi_strbuf_printf(strbuf, ";");
            }
        }
        inchi_ios_print_nodisplay(out_file, "%s%s", strbuf->pStr, pLF);

        LOG_NO_ARGS("\n******************* (L4184:ichiprt1.c) ********************\n");
        LOG_MULT_ARGS("Polymer Layer start: %s\n", strbuf->pStr);
        LOG_NO_ARGS("\n***********************************************************\n");

    exit_function:
        if (cano_nums)
        {
            inchi_free(cano_nums);
        }
        if (compnt_nums)
        {
            inchi_free(compnt_nums);
        }
        if (aprops)
        {
            inchi_free(aprops);
        }
        if (unum)
        {
            inchi_free(unum);
        }
        if (units2)
        {
            for (i = 0; i < nunits2; i++)
            {
                OAD_PolymerUnit_Free(units2[i]);
            }
            inchi_free(units2);
        }
        if (old_stars)
        {
            for (i = 0; i < pOrigStruct->polymer->n_pzz; i++)
                pOrigStruct->polymer->pzz[i] = old_stars[i];
            inchi_free(old_stars);
        }

    }
    return err;
}


/****************************************************************************
Output InChI: polymer layer, single CRU data
****************************************************************************/
static int OutputINCHI_PolymerLayer_SingleUnit( OAD_PolymerUnit *u,
                                                int bPolymers,
                                                int total_star_atoms,
                                                int *n_used_stars,
                                                OAD_AtProps *aprops,
                                                int *cano_nums,
                                                ORIG_ATOM_DATA *orig_inp_data,
                                                ORIG_STRUCT *pOrigStruct,
                                                INCHI_IOS_STRING *strbuf )
{
    int j, k, tmp, a1 = 0, a2 = 0, a3 = 0, a4 = 0, b, curr_star_num;
    int err = 0;
    OAD_Polymer *p = orig_inp_data->polymer;
    inp_ATOM    *at = orig_inp_data->at;

    /* print unit type and subtype */
    inchi_strbuf_printf( strbuf, "%-d%-d%-d-", u->type, u->subtype, u->conn );

    /* print unit atoms */
    print_sequence_of_nums_compressing_ranges( u->na, u->alist, strbuf );

    /* Print the crossing bonds or frame-shiftable pattern */
    if (u->nb > 2)
    {
        /* not supported yet, too many bonds in SBL */
        err = 12;
        goto exit_function;
    }

    /* Print crossing bonds "(cap1-partner1,cap2-partner2)"    */
    if (u->nb == 2 && ( !u->cyclizable || !u->cyclized ))
    {
        int swap = 0;
        a1 = u->blist[0];
        a2 = u->blist[1];
        a3 = u->blist[2];
        a4 = u->blist[3];
        
        if (is_in_the_ilist( u->alist, a1, u->na ))
        {
            tmp = a2;
            a2 = a1;
            a1 = tmp;
        }
        if (is_in_the_ilist( u->alist, a3, u->na ))
        {
            tmp = a4;
            a4 = a3;
            a3 = tmp;
        }

        /* Always print first the crossing bond pointing to more senior CRU end ("head")    */
        if (bPolymers==POLYMERS_LEGACY)
        {
            /* old, v. 1.05 */
            /* The first printed is the crossing bond with higher canonical number of the cap */
            swap = a3 < a1;
        }
        else
        {
            /* new in v. 1.06 */
            /* The first printed is the crossing bond pointing to more senior CRU end ("head")    */
            swap = (OAD_Polymer_IsFirstAtomRankLower(a2, a4, aprops) == 1);
        }
        
        if (swap)
        {
            inchi_strbuf_printf( strbuf, "(%-d-%-d,%-d-%-d)", a3, a4, a1, a2 );
        }
        else
        {
            inchi_strbuf_printf( strbuf, "(%-d-%-d,%-d-%-d)", a1, a2, a3, a4 );
        }
    }

    else if (u->nb <= 2 && ( u->cyclizable || u->nbkbonds > 0 ))
    {
        /* Print frame-shiftable pattern "cap1,cap2-(b1a1,b1a2, b2a1,b2a2, ... )"     */
        /* where b1, b2, ... are CRU bonds potentially invilved in frame shift          */

        /*  Get actual star atoms numbers from all-units pool                           */
        /*  NB: ordered according to already established in units2/unum                 */
        if (u->cap1 > 0 || u->cap2 > 0)
        {
            int n_expl_H = 0, pos = 0;
            char *sza = pOrigStruct->szAtoms;
            while (sza[pos])
            {
                if (sza[pos] == 'H')
                {
                    if (isupper( UCINT sza[pos + 1] ) || !sza[pos + 1])        /* if ( next_c is Uppercase or NUL ) */
                    {
                        n_expl_H++;
                    }
                    if (!sza[pos + 1])
                    {
                        break;
                    }
                }
                pos++;
            }

            if (u->cap1 > 0)
            {
                curr_star_num = pOrigStruct->num_atoms - n_expl_H - total_star_atoms + *n_used_stars + 1;
                if (curr_star_num > pOrigStruct->num_atoms)
                {
                    err = 11; goto exit_function;
                }
                a1 = curr_star_num;
                (*n_used_stars)++;
            }
            if (u->cap2 > 0)
            {
                curr_star_num = pOrigStruct->num_atoms - n_expl_H - total_star_atoms + *n_used_stars + 1;
                if (curr_star_num > pOrigStruct->num_atoms)
                {
                    err = 11; goto exit_function;
                }
                a2 = curr_star_num;
                (*n_used_stars)++;
            }
        }
        /* a1 and a2 are number of star atoms associated (but actually  */
        /* disconnected at this moment ) with SRU head and tail atoms   */
        inchi_strbuf_printf( strbuf, "(%-d,%-d-", a1, a2 );

        if (u->cyclizable == CLOSING_SRU_DIRADICAL)
        {
            inchi_strbuf_printf( strbuf, "%-d)", u->end_atom1 );
        }
        else if (u->cyclizable == CLOSING_SRU_HIGHER_ORDER_BOND)
        {
            a3 = u->end_atom1;
            a4 = u->end_atom2;
            inchi_sort_int_pair_ascending( &a3, &a4 );
            /* if ( a3 > a4 )   { tmp = a4; a4 = a3;  a3 = tmp;                  }*/
            inchi_strbuf_printf( strbuf, "%-d.%-d)", a3, a4 );
        }
        else if (u->cyclizable == CLOSING_SRU_RING)
        {
            if (u->nbkbonds == 0)
            {
                /* last resort */
                a3 = u->end_atom1;
                a4 = u->end_atom2;
                inchi_sort_int_pair_ascending( &a3, &a4 );
                /* if ( a3 > a4 ) { tmp = a4; a4 = a3; a3 = tmp; } */
                inchi_strbuf_printf( strbuf, "%-d,%-d)", a3, a4 );
            }
            else
            {
                /* Sort all backbone bonds in min-at-number order */
                for (b = 1; b < u->nbkbonds; b++)
                {
                    int *tmp_psbond = u->bkbonds[b];
                    j = b - 1;
                    while (j >= 0 && IsBondAtomNumsLesser( u->bkbonds[j], tmp_psbond ) > 0)
                    {
                        u->bkbonds[j + 1] = u->bkbonds[j];
                        j--;
                    }
                    u->bkbonds[j + 1] = tmp_psbond;
                }

                if (p->treat==POLYMERS_MODERN || p->treat==POLYMERS_LEGACY_PLUS)
                {
                    /* was #if DO_POLYMER_FRAME_SHIFT_AT_STRUCT_TO_INCHI_CONVERSION==1 */
                    /* was: OAD_PolymerUnit_ReorderPolymerFrameShiftLinks( u, orig_inp_data, aprops, cano_nums );
                    OAD_Polymer_DebugTrace(p);*/

                    /*	Find senior link and move it to the beginning of list
                        If necessary, swap atoms in link so that the first
                        points to SRU head and the next to SRU tail

                        The net result of ordering is as follows:
                        at1,at2,  at3,at4,  at5,at6, ...
                        here
                        at1, at2 is the most senior bond, and at1 is more senior than at2
                        all other pairs at3,at4,  at5,at6, ... are sorted just in increasing
                        order of first number in pair, then second one, e.g.: at3<at4; at3<=at5 
                        (and if at3==at5 then at4<at6) 
                    */

                    if (p->frame_shift_scheme != FSS_NONE && u->nbkbonds >= 1 && u->cap1 >= 1 && u->cap2 >= 1)
                    {
                        if (OAD_PolymerUnit_SetReopeningDetails(u, at))
                        {
                            /* Find senior backbone bond and move it to the beginning of list */
                            {
                                int senior_bond, bond0at1, bond0at2;
                                OAD_PolymerUnit_SortBackboneBondsAndSetSeniors(u, at, aprops, &senior_bond);
                                if (senior_bond)
                                {
                                    bond0at1 = u->bkbonds[0][0];
                                    bond0at2 = u->bkbonds[0][1];
                                    u->bkbonds[0][0] = u->bkbonds[senior_bond][0];
                                    u->bkbonds[0][1] = u->bkbonds[senior_bond][1];
                                    u->end_atom1 = u->bkbonds[0][0];
                                    u->end_atom2 = u->bkbonds[0][1];
                                    u->bkbonds[senior_bond][0] = bond0at1;
                                    u->bkbonds[senior_bond][1] = bond0at2;
                                }
                            }
                        }
                    }
                    /* p->really_do_frame_shift = 0; */
                }

                for (k = 0; k < u->nbkbonds; k++)
                {
                    a3 = u->bkbonds[k][0]; a4 = u->bkbonds[k][1];
                    /*if ( a3 > a4 )    { tmp = a4; a4 = a3; a3 = tmp; }*/
                    inchi_strbuf_printf( strbuf, "%-d,%-d%-c", a3, a4, k == u->nbkbonds - 1 ? ')' : ',' );
                }
            }
        }
    }

exit_function:

    return err;
}


/****************************************************************************
Output AuxInfo: header and normalization type
****************************************************************************/
int OutputAUXINFO_HeaderAndNormalization_type( CANON_GLOBALS    *pCG,
                                               INCHI_IOSTREAM   *out_file,
                                               INCHI_IOS_STRING *strbuf,
                                               int              bINChIOutputOptions,
                                               int              *INCHI_basic_or_INCHI_reconnected,
                                               int              num_components2[],
                                               INCHI_OUT_CTL    *io,
                                               char             *pLF,
                                               char             *pTAB )
{
    /* AuxInfo header  */
    if (*INCHI_basic_or_INCHI_reconnected == INCHI_BAS)
    {
        inchi_strbuf_printf( strbuf, "AuxInfo=" ); /* in wINChI window, separate INChI: from AuxInfo: with blank line */
        inchi_ios_print( out_file, "%s%s%s",
                                  /* blank line before AuxInfo in winchi window unless it is an annotation */
            ( bINChIOutputOptions & INCHI_OUT_WINCHI_WINDOW ) ? "\n" : "",
                                  strbuf->pStr, pLF );
        szGetTag( AuxLbl, io->nTag, io->bTag1 = AL_VERS, io->szTag1, &io->bAlways, 0 );
        inchi_strbuf_reset( strbuf ); io->tot_len = 0;
        inchi_strbuf_printf( strbuf, "%s", x_curr_ver );
        /* avoid leading slash in plain output */
        if (str_LineEnd( io->szTag1, &io->bOverflow, strbuf, -1, io->bPlainTextTags ))
        {
            return 1;
        }
        inchi_ios_print( out_file, "%s%s", strbuf->pStr, pLF );
    }
    else
    {
        if (*INCHI_basic_or_INCHI_reconnected == INCHI_REC)
        {
            szGetTag( AuxLbl, io->nTag, io->bTag1 = AL_REC_, io->szTag1, &io->bAlways, 0 );
            inchi_ios_print( out_file, "%s%s", io->szTag1, pLF );
        }
    }

    /* AuxInfo normalization type */
    if (num_components2[0] || num_components2[1])
    {
        szGetTag( AuxLbl, io->nTag, io->bTag1 = AL_NORM, io->szTag1, &io->bAlways, 0 );
        inchi_strbuf_reset( strbuf ); io->tot_len = 0;
        inchi_strbuf_printf( strbuf, "%d", ( io->bTautomeric && io->bTautomericOutputAllowed ) ? io->bTautomeric : 0 );
        if (str_LineEnd( io->szTag1, &io->bOverflow, strbuf, -1, io->bPlainTextTags ))
        {
            return 1;
        }
        inchi_ios_print( out_file, "%s%s", strbuf->pStr, pLF );
    }

    return 0;
}


/****************************************************************************
Output AuxInfo: original atom numbers and symmetry numbers (constit. equivalence /E: )
****************************************************************************/
int OutputAUXINFO_OriginalNumbersAndEquivalenceClasses( CANON_GLOBALS    *pCG,
                                                        INCHI_IOSTREAM   *out_file,
                                                        INCHI_IOS_STRING *strbuf,
                                                        int              num_components2[],
                                                        INCHI_OUT_CTL   *io,
                                                        char            *pLF,
                                                        char            *pTAB )
{
    /* Original atom numbers in order of canonical numbers */
    if (num_components2[0] || num_components2[1])
    {
        szGetTag( AuxLbl, io->nTag,
                 io->bTag1 = ( io->bSecondNonTautPass ? AL_FIXN : AL_ANBR ) | io->bFhTag, io->szTag1, &io->bAlways, 0 );
        inchi_strbuf_reset( strbuf );
        io->tot_len = 0;
        /* Original numbering output */
        io->tot_len = str_AuxNumb( pCG, io->pINChISort, io->pINChISort2,
                                   strbuf, &io->bOverflow,
                                   io->bOutType, io->TAUT_MODE, io->num_components,
                                   io->bSecondNonTautPass, io->bOmitRepetitions );

        if (str_LineEnd( io->szTag1, &io->bOverflow, strbuf, -1, io->bPlainTextTags ))
        {
            return 1;
        }
        inchi_ios_print( out_file, "%s%s", strbuf->pStr, pLF );
    }

    /*
        Symmetry numbers (constit. equivalence)    /E:
    */
    if (io->bAtomEqu[io->iCurTautMode])
    {
        /*  aux equ atoms */
        /* 1. Compare to tautomeric equivalence (in case of second, non-taut, pass only) */
        /* 2. Compare to the previous component if (1) failed to find equivalence */
        szGetTag( AuxLbl, io->nTag, io->bTag1 = AL_AEQU | io->bFhTag, io->szTag1, &io->bAlways, 0 );
        inchi_strbuf_reset( strbuf );
        io->tot_len = 0;
        io->tot_len = str_AuxEqu( io->pINChISort, io->pINChISort2,
                              strbuf, &io->bOverflow, io->bOutType, io->TAUT_MODE,
                              io->num_components, io->bSecondNonTautPass,
                              io->bOmitRepetitions, io->bUseMulipliers );

        if (str_LineEnd( io->szTag1, &io->bOverflow, strbuf, -1, io->bPlainTextTags ))
        {
            return 1;
        }
        inchi_ios_print( out_file, "%s%s", strbuf->pStr, pLF );
    }
    else
    {
        if (io->bPlainTextTags == 1)
        {
            inchi_ios_print( out_file, "/" );
        }
    }

    return 0;
}


/****************************************************************************
Output AuxInfo: tautomeric groups equivalence
****************************************************************************/
int OutputAUXINFO_TautomericGroupsEquivalence( CANON_GLOBALS    *pCG,
                                               INCHI_IOSTREAM   *out_file,
                                               INCHI_IOS_STRING *strbuf,
                                               INCHI_OUT_CTL    *io )
{
    if (io->bTautomericOutputAllowed && io->bTautomeric && io->bTautEqu[io->iCurTautMode] && !io->bSecondNonTautPass)
    {
        /*-- Tautomeric groups constitutional equivalence */

        /*-- aux tgroup equ */
        szGetTag( AuxLbl, io->nTag, io->bTag1 = AL_GEQU | io->bFhTag, io->szTag1, &io->bAlways, 0 );
        inchi_strbuf_reset( strbuf ); io->tot_len = 0;
        io->tot_len = str_AuxTgroupEqu( io->pINChISort,
                                    strbuf, &io->bOverflow, io->bOutType, io->TAUT_MODE,
                                    io->num_components, io->bUseMulipliers );
        if (str_LineEnd( io->szTag1, &io->bOverflow, strbuf, -1, io->bPlainTextTags ))
        {
            return 1;
        }
        inchi_ios_print( out_file, "%s", strbuf->pStr );
    }
    else
    {
        if (io->bTautomericOutputAllowed && io->bTautomeric)
        {
            if (io->bPlainTextTags == 1)
            {
                inchi_ios_print( out_file, "/" );
            }
        }
    }

    return 0;
}


/****************************************************************************
Output AuxInfo: stereo info
****************************************************************************/
int OutputAUXINFO_Stereo( CANON_GLOBALS     *pCG,
                           INCHI_IOSTREAM   *out_file,
                           INCHI_IOS_STRING *strbuf,
                           INCHI_OUT_CTL   *io,
                           char             *pLF,
                           char             *pTAB )
{
    /*--    Inverted stereo -- sp3 only + canonical numbering
    */
    if (io->bInvStereo[io->iCurTautMode])
    {
        szGetTag( AuxLbl, io->nTag, io->bTag1 = AL_STER | io->bFhTag, io->szTag1, &io->bAlways, 0 );
        /*-- inverted sp3 start tag */
        szGetTag( AuxLbl, io->nTag, io->bTag2 = io->bTag1 | AL_SP3I, io->szTag2, &io->bAlways, 0 );
        inchi_strbuf_reset( strbuf ); io->tot_len = 0;
        io->tot_len = str_AuxInvSp3( io->pINChISort, io->pINChISort2, strbuf,
                                 &io->bOverflow, io->bOutType, io->TAUT_MODE, io->num_components,
                                 io->bSecondNonTautPass, io->bOmitRepetitions, io->bUseMulipliers );
        if (str_LineEnd( io->szTag2, &io->bOverflow, strbuf, -1, io->bPlainTextTags ))
            return 1;
        inchi_ios_print( out_file, "%s%s", strbuf->pStr, pLF );

        /*-- inverted sp3  canonical numbering */
        if (io->bInvStereoOrigNumb[io->iCurTautMode])
        {
            szGetTag( AuxLbl, io->nTag, io->bTag2 = io->bTag1 | AL_SP3N, io->szTag2, &io->bAlways, 0 );
            inchi_strbuf_reset( strbuf ); io->tot_len = 0;

            io->tot_len = str_AuxInvSp3Numb( pCG, io->pINChISort, io->pINChISort2,
                                         strbuf, &io->bOverflow, io->bOutType,
                                         io->TAUT_MODE, io->num_components,
                                         io->bSecondNonTautPass, io->bOmitRepetitions );

            if (str_LineEnd( io->szTag2, &io->bOverflow, strbuf, -1, io->bPlainTextTags ))
            {
                return 1;
            }
            inchi_ios_print( out_file, "%s%s", strbuf->pStr, pLF );
        }
        else
        {
            if (io->bPlainTextTags == 1) inchi_ios_print( out_file, "/" );
        }
    }
    else
    {
        if (io->bPlainTextTags == 1)
        {
            inchi_ios_print( out_file, "//" );
        }
        /* Inverted stereo -- sp3 only + canonical numbering */
    }

    return 0;
}


/****************************************************************************
Output AuxInfo: isotopic info
****************************************************************************/
int OutputAUXINFO_IsotopicInfo( CANON_GLOBALS    *pCG,
                                INCHI_IOSTREAM   *out_file,
                                INCHI_IOS_STRING *strbuf,
                                int              *INCHI_basic_or_INCHI_reconnected,
                                INCHI_OUT_CTL    *io,
                                char             *pLF,
                                char             *pTAB )
{
    int i;

    /* if InChI Fixed-H isotopic is empty, then do not output corresponding AuxInfo */

    i = io->bSecondNonTautPass &&
        ( *io->pSortPrintINChIFlags & ( ( *INCHI_basic_or_INCHI_reconnected == INCHI_BAS ) ? FLAG_SORT_PRINT_NO_IFIX_H_BAS :
            FLAG_SORT_PRINT_NO_IFIX_H_REC ) );

    if (io->bIsotopic && !i &&
        ( io->bIsotopicOrigNumb[io->iCurTautMode] ||
            io->bIsotopicAtomEqu[io->iCurTautMode] ||
            (io->bTautomericOutputAllowed && io->bTautomeric && io->bIsotopicTautEqu[io->iCurTautMode]) ||
            (io->bInvIsotopicStereo[io->iCurTautMode]
            && ( io->bIgn_UU_Sp3_Iso[io->iCurTautMode])) || io->bIgn_UU_Sp2_Iso[io->iCurTautMode] ) ) /* djb-rwth: addressing LLVM warnings */
    {
        /*-- isotopic aux info header */
        szGetTag( AuxLbl, io->nTag, io->bTag1 = AL_ISOT | io->bFhTag, io->szTag1, &io->bAlways, 0 );
        inchi_strbuf_reset( strbuf ); /* pStr[io->tot_len = 0] = '\0'; */
        /*-- Original atom numbers in order of isotopic canonical numbers */
        szGetTag( AuxLbl, io->nTag, io->bTag2 = io->bTag1 | AL_ISON, io->szTag2, &io->bAlways, 0 );
        if (io->bIsotopicOrigNumb[io->iCurTautMode])
        {
            inchi_strbuf_reset( strbuf );
            io->tot_len = 0;
            io->tot_len = str_AuxIsoNumb( pCG, io->pINChISort, io->pINChISort2,
                                      strbuf, &io->bOverflow, io->bOutType,
                                      io->TAUT_MODE, io->num_components,
                                      io->bSecondNonTautPass, io->bOmitRepetitions );
            if (str_LineEnd( io->szTag2, &io->bOverflow, strbuf, -1, io->bPlainTextTags ))
            {
                return 1;
            }
            inchi_ios_print( out_file, "%s%s", strbuf->pStr, pLF );
        }
        else
        {
            /*if ( io->bPlainTextTags == 1 ) inchi_ios_print( out_file, "/" );*/
            inchi_ios_print( out_file, "%s%s", io->szTag2, pLF ); /* mark isotopic output */
        }

        /*-- Isotopic symmetry */
        if (io->bIsotopicAtomEqu[io->iCurTautMode])
        {
            /*-- atoms */
            szGetTag( AuxLbl, io->nTag, io->bTag2 = io->bTag1 | AL_AEQU, io->szTag2, &io->bAlways, 0 );
            inchi_strbuf_reset( strbuf ); io->tot_len = 0;
            io->tot_len = str_AuxIsoEqu( io->pINChISort, io->pINChISort2,
                                     strbuf,
                                     &io->bOverflow, io->bOutType, io->TAUT_MODE, io->num_components,
                                     io->bSecondNonTautPass, io->bOmitRepetitions, io->bUseMulipliers );
            if (str_LineEnd( io->szTag2, &io->bOverflow, strbuf, -2/*was -1: Fix15*/, io->bPlainTextTags ))
            {
                return 1;
            }
            inchi_ios_print( out_file, "%s%s", strbuf->pStr, pLF );
        }
        else
        {
            if (io->bPlainTextTags == 1)
            {
                inchi_ios_print( out_file, "/" );
            }
        }

        /*-- Tautomeric groups, isotopic */
        if (io->bTautomericOutputAllowed && io->bTautomeric && io->bIsotopicTautEqu[io->iCurTautMode])
        {
            /*-- Isotopic tautomeric groups equivalence */
            szGetTag( AuxLbl, io->nTag, io->bTag2 = io->bTag1 | AL_GEQU, io->szTag2, &io->bAlways, 0 );
            inchi_strbuf_reset( strbuf ); io->tot_len = 0;
            io->tot_len = str_AuxIsoTgroupEqu( io->pINChISort,
                                           strbuf, &io->bOverflow,
                                           io->bOutType, io->TAUT_MODE, io->num_components,
                                           io->bOmitRepetitions, io->bUseMulipliers );
            if (str_LineEnd( io->szTag2, &io->bOverflow, strbuf, -2/*was -1: Fix15*/, io->bPlainTextTags ))
            {
                return 1;
            }
            inchi_ios_print( out_file, "%s%s", strbuf->pStr, pLF );
        }
        else
        {
            if (io->bTautomericOutputAllowed && io->bTautomeric)
            {
                if (io->bPlainTextTags == 1)
                {
                    inchi_ios_print( out_file, "/" );
                }
            }
        }
        /*-- Isotopic inverted stereo */
        if (io->bInvIsotopicStereo[io->iCurTautMode])
        {
            szGetTag( AuxLbl, io->nTag, io->bTag2 = io->bTag1 | AL_STER, io->szTag2, &io->bAlways, 0 );
            /*-- inverted isotopic sp3 start tag */
            szGetTag( AuxLbl, io->nTag, io->bTag3 = io->bTag2 | AL_SP3I, io->szTag3, &io->bAlways, 0 );
            inchi_strbuf_reset( strbuf ); io->tot_len = 0;
            io->tot_len = str_AuxInvIsoSp3( io->pINChISort, io->pINChISort2,
                                        strbuf, &io->bOverflow,
                                        io->bOutType, io->TAUT_MODE, io->num_components,
                                        io->bSecondNonTautPass, io->bOmitRepetitions, io->bUseMulipliers );
            if (str_LineEnd( io->szTag3, &io->bOverflow, strbuf, -1, io->bPlainTextTags ))
            {
                return 1;
            }
            inchi_ios_print( out_file, "%s", strbuf->pStr );
            /*-- inverted isotopic sp3  canonical numbering */
            if (io->bInvIsotopicStereoOrigNumb[io->iCurTautMode])
            {
                szGetTag( AuxLbl, io->nTag, io->bTag3 = io->bTag2 | AL_SP3N, io->szTag3, &io->bAlways, 0 );
                inchi_strbuf_reset( strbuf ); io->tot_len = 0;
                io->tot_len = str_AuxInvIsoSp3Numb( pCG, io->pINChISort, io->pINChISort2,
                                                strbuf, &io->bOverflow,
                                                io->bOutType, io->TAUT_MODE,
                                                io->num_components,
                                                io->bSecondNonTautPass,
                                                io->bOmitRepetitions );

                if (str_LineEnd( io->szTag3, &io->bOverflow, strbuf, -1, io->bPlainTextTags ))
                {
                    return 1;
                }
                inchi_ios_print( out_file, "%s%s", strbuf->pStr, pLF );
            }
            else
            {
                if (io->bPlainTextTags == 1)
                {
                    inchi_ios_print( out_file, "/" );
                }
            }
        }
        else
        {
            if (io->bPlainTextTags == 1)
            {
                inchi_ios_print( out_file, "//" );
            }
        }
        /*-- totally omitted undefined/unknown isotopic stereo */
    } /* Aux info isotopic */

    return 0;
}


/****************************************************************************
Output AuxInfo: charges, radicals, unusual valences
****************************************************************************/
int OutputAUXINFO_ChargesRadicalsAndUnusualValences( CANON_GLOBALS    *pCG,
                                                     INCHI_IOSTREAM   *out_file,
                                                     INCHI_IOS_STRING *strbuf,
                                                     INCHI_OUT_CTL    *io,
                                                     char             *pLF,
                                                     char             *pTAB )
{
    if (!io->bSecondNonTautPass && io->bChargesRadVal[io->iCurTautMode])
    {
        /*  aux equ atoms */
        /* 1. Compare to tautomeric equivalence (in case of second, non-taut, pass only) */
        /* 2. Compare to the previous component if (1) failed to find equivalence */
        szGetTag( AuxLbl, io->nTag, io->bTag1 = AL_CRV_ | io->bFhTag, io->szTag1, &io->bAlways, 0 );

        inchi_strbuf_reset( strbuf );
        io->tot_len = 0;

        io->tot_len = str_AuxChargeRadVal( io->pINChISort, strbuf,
                                           &io->bOverflow, io->bOutType, io->TAUT_MODE,
                                           io->num_components, io->bUseMulipliers );

        if (str_LineEnd( io->szTag1, &io->bOverflow, strbuf, -1, io->bPlainTextTags ))
        {
            return 1;
        }

        inchi_ios_print( out_file, "%s%s", strbuf->pStr, pLF );
    }

    return 0;
}


/****************************************************************************
Output AuxInfo: reversibility info (to restore orig. structure)
****************************************************************************/
int OutputAUXINFO_ReversibilityInfo( CANON_GLOBALS    *pCG,
                                     INCHI_IOSTREAM   *out_file,
                                     INCHI_IOS_STRING *strbuf,
                                     ORIG_STRUCT      *pOrigStruct,
                                     INCHI_OUT_CTL    *io,
                                     char             *pLF,
                                     char             *pTAB )
{
    if (!io->bSecondNonTautPass &&
         pOrigStruct && pOrigStruct->num_atoms &&
         pOrigStruct->szAtoms
         && pOrigStruct->szBonds
         && pOrigStruct->szCoord)
    {
        int length, cur_pos, line_len, last_pos, nMaxLineLen;
        char *p;
        nMaxLineLen = inchi_min( 80, strbuf->nAllocatedLength ); /* restrict line length to 80 characters */

        szGetTag( AuxLbl, io->nTag, io->bTag1 = AL_REVR | io->bFhTag, io->szTag1, &io->bAlways, 0 );

        /* Atoms /A: */
        szGetTag( AuxLbl, io->nTag, io->bTag2 = io->bTag1 | AL_ATMR, io->szTag2, &io->bAlways, 0 );
        inchi_strbuf_reset( strbuf );
        inchi_ios_print( out_file, "%s%s", io->szTag2, strbuf->pStr );
        p = pOrigStruct->szAtoms;
        length = (int) strlen( p );
        io->tot_len = strbuf->nUsedLength;
        line_len = nMaxLineLen - io->tot_len;
        for (cur_pos = 0; cur_pos < length; cur_pos = last_pos)
        {
            if (length - cur_pos >= line_len)
            {
                last_pos = cur_pos + line_len;
                /* search backward for the nearest first atom letter (always uppercase) */
                while (cur_pos < last_pos && !isupper( UCINT p[last_pos] ))
                {
                    last_pos--;
                }
            }
            else
            {
                last_pos = length;
            }
            if (last_pos > cur_pos)
            {
                memcpy(strbuf->pStr + strbuf->nUsedLength, p + cur_pos, (long long)last_pos - (long long)cur_pos); /* djb-rwth: cast operators added */
                strbuf->pStr[strbuf->nUsedLength + last_pos - cur_pos] = '\0';
                /*strbuf->nUsedLength = strbuf->nUsedLength + last_pos - cur_pos;*/

                if (1) /* always show "Zy" as "Zz" */
                {
                    char *pzy, *pstart=strbuf->pStr + strbuf->nUsedLength;
                    while ((pzy = strstr( pstart, "Zy" ))) /* djb-rwth: addressing LLVM warning */
                    {
                        *(++pzy) = 'z';
                        pstart = pzy;
                    }
                }

                inchi_ios_print( out_file, "%s%s", strbuf->pStr, io->bPlainTextTags ? "" : "\n" );
            }
            else
            {
                break;
            }
        }
        if (pLF[0])
        {
            inchi_ios_print( out_file, "%s", pLF );
        }

        inchi_strbuf_reset( strbuf );

        /* Bonds /B: */
        szGetTag( AuxLbl, io->nTag, io->bTag2 = io->bTag1 | AL_BNDR, io->szTag2, &io->bAlways, 0 );
        inchi_strbuf_reset( strbuf );
        inchi_ios_print( out_file, "%s%s", io->szTag2, strbuf->pStr );

        p = pOrigStruct->szBonds;
        length = (int) strlen( p );
        line_len = nMaxLineLen - io->tot_len;
        for (cur_pos = 0; cur_pos < length; cur_pos = last_pos)
        {
            if (length - cur_pos >= line_len)
            {
                last_pos = cur_pos + line_len - 1;
                /* search backward for the nearest first bond delimiter ";" */
                while (cur_pos < last_pos && p[last_pos] != ';')
                {
                    last_pos--;
                }
                if (cur_pos < last_pos)
                {
                    last_pos++; /* include ';' at the end of the line */
                }
            }
            else
            {
                last_pos = length;
            }
            if (last_pos > cur_pos)
            {
                memcpy(strbuf->pStr, p + cur_pos, (long long)last_pos - (long long)cur_pos); /* djb-rwth: cast operators added */
                strbuf->pStr[last_pos - cur_pos] = '\0';
                strbuf->nUsedLength = last_pos - cur_pos;
                inchi_ios_print( out_file, "%s%s", strbuf->pStr, io->bPlainTextTags ? "" : "\n" );
                inchi_strbuf_reset( strbuf );
            }
            else
            {
                break;
            }
        }
        if (pLF[0])
        {
            inchi_ios_print( out_file, "%s", pLF );
        }

        /* Coordinates /C:    */
        szGetTag( AuxLbl, io->nTag, io->bTag2 = io->bTag1 | AL_XYZR, io->szTag2, &io->bAlways, 0 );
        inchi_strbuf_reset( strbuf );
        inchi_ios_print( out_file, "%s%s", io->szTag2, strbuf->pStr );

        p = pOrigStruct->szCoord;
        length = (int) strlen( p );
        line_len = nMaxLineLen - io->tot_len;
        for (cur_pos = 0; cur_pos < length; cur_pos = last_pos)
        {
            if (length - cur_pos >= line_len)
            {
                last_pos = cur_pos + line_len - 1;
                /* search backward for the nearest first coord. delimiter ";" */
                while (cur_pos < last_pos && p[last_pos] != ';')
                {
                    last_pos--;
                }
                if (cur_pos < last_pos)
                {
                    last_pos++; /* include ';' at the end of the line */
                }
            }
            else
            {
                last_pos = length;
            }
            if (last_pos > cur_pos)
            {
                memcpy(strbuf->pStr, p + cur_pos, (long long)last_pos - (long long)cur_pos); /* djb-rwth: cast operator added */
                strbuf->pStr[last_pos - cur_pos] = '\0';
                strbuf->nUsedLength = last_pos - cur_pos;
                inchi_ios_print( out_file, "%s%s", strbuf->pStr, io->bPlainTextTags ? "" : "\n" );
                inchi_strbuf_reset( strbuf );
            }
            else
            {
                break;
            }
        }

        if (pLF[0])
        {
            inchi_ios_print( out_file, "%s", pLF );
        }
    }

    return 0;
}


/****************************************************************************/
int OutputAUXINFO_PolymerInfo( CANON_GLOBALS    *pCG,
                               INCHI_IOSTREAM   *out_file,
                               INCHI_IOS_STRING *strbuf,
                               ORIG_STRUCT      *pOrigStruct,
                               INCHI_OUT_CTL    *io,
                               char             *pLF,
                               char             *pTAB )
{
    int k, i, q;
    OAD_Polymer *p;
    OAD_PolymerUnit *u;


    if (!pOrigStruct)
    {
        return 0;
    }
    p = pOrigStruct->polymer;
    if (!p)
    {
        return 0;
    }

    inchi_strbuf_reset( strbuf );

    inchi_ios_print( out_file, "/Z:" );


    /* Print polymer units data */
    for (i = 0; i < p->n; i++)
    {
        /* For each unit u ... */
        u = p->units[i];

        /* print kinds of unit */
        inchi_strbuf_printf( strbuf, "%-d%-d%-d-", u->type, u->subtype, u->conn );
        inchi_strbuf_printf( strbuf, "%-s-", u->smt[0] ? u->smt : "n" );

        /* Print unit atoms */
        print_sequence_of_nums_compressing_ranges( u->na, u->alist, strbuf );

        /* Print bonds from unit to otside */
        if (u->nb > 0)
        {
            inchi_strbuf_printf( strbuf, "(" );
            for (k = 0; k < 2 * u->nb - 1; k++)
            {
                inchi_strbuf_printf( strbuf, "%-d,", u->blist[k] );
            }
            inchi_strbuf_printf( strbuf, "%-d)", u->blist[2 * u->nb - 1] );
        }

        if (fabs( -fabs( u->xbr1[0] ) + 777777.777 ) > 1.e-7)
        {
            inchi_strbuf_printf( strbuf, "[" );
            for (q = 0; q < 3; q++)
            {
                inchi_strbuf_printf( strbuf, "%-f,", u->xbr1[q] );
            }
            inchi_strbuf_printf( strbuf, "%-f]", u->xbr1[3] );
        }
        if (fabs( -fabs( u->xbr2[0] ) + 777777.777 ) > 1.e-7)
        {
            inchi_strbuf_printf( strbuf, "[" );
            for (q = 0; q < 3; q++)
            {
                inchi_strbuf_printf( strbuf, "%-f,", u->xbr2[q] );
            }
            inchi_strbuf_printf( strbuf, "%-f]", u->xbr2[3] );
        }

        if (i < p->n - 1)
        {
            inchi_strbuf_printf( strbuf, ";" );
        }
    }

    inchi_ios_print( out_file, "%s%s", strbuf->pStr, pLF );

    return 0;
}


/****************************************************************************/
int IsBondAtomNumsLesser( int *bond1, int* bond2 )
{
    int min1 = inchi_min( bond1[0], bond1[1] );
    int min2 = inchi_min( bond2[0], bond2[1] );
    int max1 = inchi_max( bond1[0], bond1[1] );
    int max2 = inchi_max( bond2[0], bond2[1] );

    if (min1 < min2)
    {
        return -1;
    }
    if (min1 > min2)
    {
        return 1;
    }
    if (min1 == min2)
    {
        if (max1 < max2)
        {
            return -1;
        }
        if (max1 > max2)
        {
            return 1;
        }
    }

    return 0;
}


/****************************************************************************/
void EditINCHI_HidePolymerZz(INCHI_IOSTREAM *out, int n_pzz, int n_zy)
{
    char *s = out->s.pStr, *s0, *buf = NULL;
    char prev_layer_symbol = '0';
    int i, j, ii, nzz,
        nzz1 = 0, nslash = 0, ncopied = 0,
        start = 0, skip = 0, is_in_z_layer = 0,
        eol_was_consumed = 0, pre_eol = 0,
        nonprt_sym = 0, nonprt_prev = 0;

    if (n_zy > 0) 
    {
        /* We have some placeholder pseudo atoms which should not be removed below (if anyway they are allowed) */
        if (n_pzz == 0)
        {
            /* Have nothing to remove, just exit */
            return;
        }
        /* Have both polymer-related and placeholder pseudo atoms */
        if (n_pzz < 2)
        {
            /* Something strange, should not arrive here, cowardly exit */
            return;
        }
    }

    /* Ensure that polymeric layer is present */
    if (!strstr(s, "/z"))
    {
        return;
    }
    s0 = strstr(s, "InChI=1B/");
    if (!s0)
    {
        return;
    }

#if 0
    nzz1 = CountPseudoElementInFormula("Zz", s0 + strlen("InChI=1B/"));
    if (nzz1 == 0)
    {
        return;
    }
    if (nzz1 != (n_pzz + n_zy))
    {
        /* Something strange, should not arrive here, cowardly exit */
        return;
    }
#endif
    nzz1 = n_pzz;

    /* OK, we must hide n_pzz Zz's*/
    buf = (char *) inchi_calloc( (long long)out->s.nUsedLength + 1, sizeof( char ) ); /* djb-rwth: cast operator added */
    if (!buf)
    {
        return;
    }

    /* Consume '\n' temporarily */
    if (s[out->s.nUsedLength - 1] == '\n')
    {
        s[out->s.nUsedLength - 1] = '\0';
        out->s.nUsedLength--;
        eol_was_consumed = 1;
    }

    start = s0 - s;
    nzz = nzz1;
    is_in_z_layer = skip = 0;
    for (i = start; i < out->s.nUsedLength; i++)
    {
        pre_eol = (i == out->s.nUsedLength - 1);
        nonprt_sym = s[i] == '\n' || s[i] == '\r' || s[i] == '\t';

        if (!skip)
        {
            buf[ncopied] = s[i];
            ncopied++;
            if (nonprt_sym && nonprt_prev)
            {
                continue;
            }
        }
        nonprt_prev = nonprt_sym;

        if (is_in_z_layer && !skip)
        {
            if (s[i] == '(')
            {
                /* Software version 1.07 : skip pattern "(cap,cap-bkbonds)" but not "(cap-end, cap-end)" */
                const char *q;
                const char *p = out->s.pStr + i + 1;
                AT_NUMB ia = (AT_NUMB) inchi_strtol(p, &q, 10); /* make compiler happy: */ /* djb-rwth: removing redundant code; ignoring LLVM warning: variable used to store function return value */
                if (*q != '-')
                {
                    skip = 1; 
                }
            }
        }
        else if (is_in_z_layer && skip)
        {
            if (s[i] == '-')
            {
                skip = 0;
            }
        }

        if (s[i] == '/' || pre_eol || nonprt_sym )
        {
            if (is_in_z_layer)
            {
                is_in_z_layer = 0;
            }

            if (s[i] == '/')
            {
                nslash++;
            }

            if (nslash == 2 ||
                ( nslash == 1 && pre_eol ) ||
                ( prev_layer_symbol == 'f' )
                )
            {
                if (nzz)
                {
                    /* eat Zz's */
                    ii = i;
                    if (pre_eol)
                    {
                        ii = i + 1;
                    }
                    if (s[ii - 1] == 'z' && s[ii - 2] == 'Z')
                    {
                        ncopied -= 2;
                        for (j = ii - 3; j >= 0; j--)
                        {
                            if (s[j] == '.')
                            {
                                break;
                            }
                            ncopied--;
                        }
                        ncopied--;
                        if (!pre_eol)
                        {
                            buf[ncopied - 1] = '/';
                        }
                        else
                        {
                            buf[ncopied - 1] = '\0';
                        }
                    }
                }
            }
            else if (nslash > 2 || pre_eol || s[i] == '\n')
            {
                if (nzz)
                {
                    if (prev_layer_symbol != 'p' &&
                         prev_layer_symbol != 's' &&
                         prev_layer_symbol != 'f' &&
                         prev_layer_symbol != 'z'
                        )
                    {
                        /* eat nzz last ; if any */
                        int n_eaten = 0;
                        char eatable = ';';
                        if (prev_layer_symbol == 'm')
                        {
                            eatable = '.';
                        }
                        ii = i;
                        if (pre_eol)
                        {
                            ii = i + 1;
                        }
                        for (j = ii - 1; j >= 0; j--)
                        {
                            if (s[j] == eatable  && n_eaten < nzz)
                            {
                                ncopied--;
                                n_eaten++;
                            }
                            else
                            {
                                break;
                            }
                        }
                        if (!pre_eol)
                        {
                            if (s[i] == '\n')
                            {
                                buf[ncopied - 1] = '\n';
                            }
                            else
                            {
                                buf[ncopied - 1] = '/';
                            }
                        }
                        else
                        {
                            buf[ncopied] = '\0';
                            break;
                        }
                    }
                }
            }

            prev_layer_symbol = s[i + 1];
            if (s[i + 1] == 'r')
            {
                if (i + 3 < out->s.nUsedLength && s[i + 3] == ':')
                {
                    /* ra: rB; rC: in AuxInfo */
                    ;
                }
                else
                {
                    /* reconnected layer */
                    nslash = 1;
                    /* nzz = nzz2;    paranoidal */
                }
            }
            else if (s[i + 1] == 'z')
            {
                is_in_z_layer = 1;
            }
        }
    }

    out->s.nUsedLength = 0;
    inchi_ios_print_nodisplay( out, "%s%s", buf, eol_was_consumed ? "\n" : "" );
    inchi_free( buf );

    return;
}


/****************************************************************************/
int CountPseudoElementInFormula( const char *pseudo, char *s ) /* djb-rwth: ignoring LLVM warning: function used */
{
    int npseudo=0, mult=1, index=1, new_component=1;
    const char *p, *q;
    char prev = '/';

    /*
        format is 
        [sequence of] [.[int_mult[Zz[int_index]]]]
    */

    if (!s)
    {
        return 0;
    }
    
    p = s; 
    while (*p)
    /*for (p = s ; *p; p++)*/
    {
        if (*p =='/')
        {
            /* end of formula layer */
            break;
        }
        else if (*p == '.')
        {
            /* start of new component subformula */
            new_component = 1;
            prev = *p;
            p++;
            continue;
        }

        if (new_component)
        {
            new_component = 0;
            if (isdigit(*p))
            {
                mult = (int)inchi_strtol(p, &q, 10);
                p = q;
                prev = *q--; 
                continue;
            }
            else
            {
                mult = 1;
            }
            if (!mult)
            {
                break; 
            }
        }
        else if (*p== pseudo[1] && prev== pseudo[0])
        {
            q = p++;
            if ( *q && isdigit(*q))
            {
                index = (int)inchi_strtol(q, &p, 10);
            }
            else
            {
                index = 1;
                p = q;
            }
            npseudo+= mult*index;
        }
        prev = *p;
        p++;
    }

    return npseudo;
}


/****************************************************************************
    Get canonical numbers and component numbers for each original atom
    NB:  cano_nums[orig_num]   - in orig_nums order,
         compnt_nums[cano_num] - in cano_nums  order
    NB': orig_nums   are 1-based
         cano_nums   are 0-based
         compnt_nums are 1-based
****************************************************************************/
int InternallyGetCanoNumsAndComponentNums( CANON_GLOBALS    *pCG,
                                           INCHI_IOS_STRING *strbuf,
                                           INCHI_OUT_CTL    *io,
                                           int              nat,
                                           int              *cano_nums,
                                           int              *compnt_nums )
{
    int orig_num, cano_num, icompnt, i, k, ndigit, err;
    char c, cnum[8];

    if (!cano_nums || !compnt_nums || !strbuf->pStr)
    {
        return 1;
    }

    inchi_strbuf_reset( strbuf );
    io->tot_len = str_AuxNumb( pCG, io->pINChISort, io->pINChISort2,
                               strbuf, &io->bOverflow, io->bOutType,
                               io->TAUT_MODE, io->num_components,
                               io->bSecondNonTautPass, io->bOmitRepetitions );
    for (i = 0; i < nat; i++)
    {
        compnt_nums[i] = -1;
        cano_nums[i + 1] = -1;
    }

    ndigit = 0;
    err = 0;
    /* djb-rwth: removing redundant code */
    icompnt = 1;
    cano_num = 0;
    for (k = 0; k <= strbuf->nUsedLength; k++)
    {
        c = strbuf->pStr[k];
        if (c == ',' || c == ';' || c == '\0')
        {
            cnum[ndigit] = '\0';
            orig_num = atoi( cnum );
            cano_nums[orig_num] = cano_num;
            compnt_nums[cano_num] = icompnt;
            cnum[0] = '\0';
            ndigit = 0;
            cano_num++;
            if (c == ';')
            {
                icompnt++;
            }
            if (c == '\0')
            {
                break;
            }
            continue;
        }
        else if (isdigit( c ))
        {
            cnum[ndigit] = c;
            ndigit++;
        }
        else
        {
            err = 2;
            goto exit_function;
        }
    }

exit_function:
    inchi_strbuf_reset( strbuf );

    return err;
}


/***************************************************************************/
int MergeZzInHillFormula(INCHI_IOS_STRING *strbuf)
{
    char *p, *scopy = NULL, *stmp=NULL, *pend=NULL, *p0 = NULL; /* djb-rwth: removing redundant variables */
    size_t sublen; /* djb-rwth: removing redundant variables */

    if (!strbuf->pStr || strbuf->nUsedLength < 1)
    {
        return 0;
    }
    scopy = (char *)inchi_calloc((long long)strbuf->nAllocatedLength+1, sizeof(char)); /* djb-rwth: cast operator added */
    if (!scopy)
    {
        inchi_free(scopy); /* djb-rwth: avoiding memory leak */
        return -1; /* failed */
    }    
    memcpy(scopy, strbuf->pStr, strbuf->nAllocatedLength);
    stmp = (char *)inchi_calloc((long long)strbuf->nAllocatedLength + 1, sizeof(char)); /* djb-rwth: cast operator added */
    if (!stmp)
    {
        inchi_free(scopy); /* djb-rwth: avoiding memory leak */
        return -1; /* failed */
    }

    inchi_strbuf_reset(strbuf);
    p0 = scopy;
    p = p0;
    do 
    {
        /* djb-rwth: removing redundant code */
        pend = strchr(p, '.');
        if (!pend)
        {
            pend = strchr(p, '\0');
        }
        sublen = pend - p;
        memcpy(stmp, p, sublen);
        stmp[sublen] = '\0';
        MergeZzInStrHillFormulaComponent(stmp);
        if (stmp)
        {
            inchi_strbuf_printf(strbuf, "%-s%-c", stmp, *pend);
        }
        /* djb-rwth: removing redundant code */
    } while ( *pend && (p=pend+1));


    if (scopy)
    {
        inchi_free(scopy);
    }
    if (stmp)
    {
        inchi_free(stmp);
    }

    return 0;
}


/***************************************************************************/
void MergeZzInStrHillFormulaComponent(char *s)
{
    char *pz = strstr(s, "Zz");
    if (pz)
    {
        int n = 1, offset;
        char *pz2 = pz + 2, *pd = pz + 2;
        if (*pd && (isdigit(*pd)))
        {
            n = strtol(pd, &pz2, 10);
        }
        pz2 = strstr(pz2, "Zz");
        if (pz2)
        {
            int n2 = 1;
            char *pd2 = pz2 + 2;
            if (*pd2 && (isdigit(*pd2)))
            {
                n2 = strtol(pd2, &pz2, 10);
            }
            n += n2;
            offset = (int)(pd - s);
            sprintf(s + offset, "%d", n);
        }
    }
    return;
}


/****************************************************************************/
static void inchi_sort_int_pair_ascending(int *a, int *b)
{
    int tmp;
    if (*a > *b)
    {
        tmp = *b;
        *b = *a;
        *a = tmp;
    }

    return;
}



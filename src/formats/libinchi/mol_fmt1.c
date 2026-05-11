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
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>

#include "mode.h"
#include "mol_fmt.h"

#include "ichierr.h"
#include "util.h"
#include "ichi_io.h"
#include "strutil.h"

#include "bcf_s.h"

/*
    MolFile related procedures - 1


*/



/*
    Local functions; static
*/
/*    MOL V2000 */
static int MolfileReadSTextBlock( MOL_FMT_CTAB* ctab, INCHI_IOSTREAM *inp_file,
                                  int err, char *pStrErr );
static int MolfileReadPropBlock( MOL_FMT_CTAB* ctab, MOL_FMT_HEADER_BLOCK *pHdr,
                                 INCHI_IOSTREAM *inp_file,
                                 int treat_polymers,
                                 int err, char *pStrErr,
                                 int bNoWarnings );
static int MolfileReadSgroupOfPolymer( MOL_FMT_CTAB* ctab,
                                       MOL_FMT_HEADER_BLOCK *pHdr,
                                       INCHI_IOSTREAM *inp_file,
                                       char line[MOL_FMT_INPLINELEN],
                                       char *szType, char *p,
                                       int err, char *pStrErr );
static int MolfileTreatPseudoElementAtoms( MOL_FMT_CTAB* ctab,
                                           int pseudos_allowed,
                                           int *err,
                                           char *pStrErr );


/****************************************************************************
 Read MOL-format data from SD or MOL file
****************************************************************************/
MOL_FMT_DATA* ReadMolfile( INCHI_IOSTREAM *inp_file,
                           MOL_FMT_HEADER_BLOCK *OnlyHeaderBlock,
                           MOL_FMT_CTAB *OnlyCTab,
                           int bGetOrigCoord,
                           int treat_polymers,
                           int treat_NPZz,
                           char *pname,
                           int lname,
                           unsigned long *Id,
                           const char *pSdfLabel,
                           char *pSdfValue,
                           int *err,
                           char *pStrErr,
                           int bNoWarnings )
{
    MOL_FMT_DATA* mfdata;

    if (pname && lname)
    {
        pname[0] = '\0';
    }
    if (Id)
    {
        *Id = 0LU;  /* ignore for now */
    }

    mfdata = MolfileReadDataLines( inp_file, OnlyHeaderBlock, OnlyCTab,
                                   bGetOrigCoord, treat_polymers,
                                   err, pStrErr, bNoWarnings );

    if (*err < 0)
    {
        /* read OK, and end of data encountered */
        *err = -*err;
    }
    else
    {
        /* unnecessary extra data may have present in SDF; skip them for now */
        int ret_skip_extras = SDFileSkipExtraData( inp_file, Id, NULL, 0,
                                                   pname, lname, *err,
                                                   pSdfLabel, pSdfValue,
                                                   pStrErr, bNoWarnings);


        if (ret_skip_extras)
        {
            /* important to continue to the next good structure */
            *err = ret_skip_extras;
        }
    }

    /*  Treat star/Zz atoms if present.
        This first-line treating affects only Molfile input.
        No direct input in API calls, no InChI strings.
        These will be processed on checking internal data structs.
    */
    if (mfdata)
    {
        int nzz = 0; /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
        int pseudos_allowed = (treat_NPZz == 1) || (treat_polymers != POLYMERS_NO);		
        nzz = MolfileTreatPseudoElementAtoms( &mfdata->ctab,
                                              pseudos_allowed,
                                              err,
                                              pStrErr ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    }

    return mfdata;
}


/****************************************************************************
 Read data lines completely ingnore STEXT block, queries, and 3D features
****************************************************************************/
MOL_FMT_DATA * MolfileReadDataLines( INCHI_IOSTREAM *inp_file,
                                     MOL_FMT_HEADER_BLOCK *OnlyHeaderBlock,
                                     MOL_FMT_CTAB *OnlyCTab,
                                     int bGetOrigCoord,
                                     int treat_polymers,
                                     int *err,
                                     char *pStrErr,
                                     int bNoWarnings )
{
    int n_alloc_atoms;
    MOL_FMT_CTAB  ctab, *pCtab = NULL;
    MOL_FMT_HEADER_BLOCK *pHdr = NULL;
    MOL_FMT_DATA *mfdata = NULL;
    int retcode, prevcode;
    int data_ended = 0, should_read_all = 0;

    if (!OnlyHeaderBlock)
    {
        should_read_all = 1;
    }

    /* djb-rwth: removing redundant code */
    /* djb-rwth: addressing coverity ID #499502 -- TREAT_ERR_AND_FIN properly used in all cases */
    *err = 0;

    if (should_read_all)
    {

        mfdata = (MOL_FMT_DATA*) inchi_calloc( 1, sizeof( MOL_FMT_DATA ) );
        if (!mfdata)
        {
            retcode = 1;
            AddErrorMessage( pStrErr, "Out of RAM" );
            goto err_fin;
        }

        pHdr = &mfdata->hdr;
        pCtab = &mfdata->ctab;
    }
    else
    {
        pHdr = OnlyHeaderBlock;
        pCtab = OnlyCTab ? OnlyCTab : &ctab;
        memset( pHdr, 0, sizeof( MOL_FMT_HEADER_BLOCK ) ); /* djb-rwth: memset_s C11/Annex K variant? */
        memset( pCtab, 0, sizeof( MOL_FMT_CTAB ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    }

    pCtab->bonds = NULL;
    pCtab->atoms = NULL;
    pCtab->coords = NULL;
    pCtab->v3000 = NULL;

    /* Header lines */
    retcode = MolfileReadHeaderLines( pHdr, inp_file, pStrErr );
    if (retcode)
    {
        /*  most likely end of file */
        retcode += 10;
        goto err_fin;
    }

    /* Read counts line and also check if we deal with V3000 Molfile */
    retcode = MolfileReadCountsLine( pCtab, inp_file, pStrErr );
    if (retcode)
    {
        retcode += 20;
        goto err_fin;
    }

    if (pCtab->v3000)
    {

        /*    In V3000, the previously read counts should be neglected
            and new counts line (1st in Ctab) should be read preceded
            by "M  V30 BEGIN CTAB"                                    */

        retcode = MolfileV3000ReadCTABBeginAndCountsLine( pCtab, inp_file, pStrErr );
        if (retcode)
        {
            retcode += 70;
            TREAT_ERR_AND_FIN( retcode, 1, err_fin, pStrErr );
        }

        retcode = MolfileV3000Init( pCtab, pStrErr );
        if (retcode)
        {
            retcode += 70;
            TREAT_ERR_AND_FIN( retcode, 1, err_fin, pStrErr );
        }
    }

    /* Atomic block */
    if (should_read_all)
    {
        n_alloc_atoms = inchi_max( mfdata->ctab.n_atoms, 1 );
        mfdata->ctab.atoms = (MOL_FMT_ATOM*)
            inchi_calloc( n_alloc_atoms, sizeof( MOL_FMT_ATOM ) );
        if (!mfdata->ctab.atoms)
        {
            retcode = 2;
            TREAT_ERR_AND_FIN( retcode, 2, err_fin, "Out of RAM" );
        }

        if (bGetOrigCoord)
        {
            mfdata->ctab.coords = (MOL_COORD*)
                inchi_calloc( n_alloc_atoms, sizeof( MOL_COORD ) );
            if (!mfdata->ctab.coords)
            {
                retcode = 2;
                TREAT_ERR_AND_FIN( retcode, 2, err_fin, "Out of RAM" );
            }
        }
    }

    if (!pCtab->v3000)
    {
        retcode = MolfileReadAtomsBlock( pCtab, inp_file, retcode, pStrErr );
    }
    else
    {
        retcode = MolfileV3000ReadAtomsBlock( pCtab, inp_file, retcode, pStrErr );
    }
    if (retcode)
    {
        if (retcode < 0)
        {
            retcode = -retcode;
            data_ended = 1;
        }
        retcode += 30;
        /* goto err_fin; */
    }

    /* Bonds block */
    if (should_read_all && retcode < 30)
    {
        if (!data_ended)
        {
            int n_alloc_bonds = inchi_max( mfdata->ctab.n_bonds, 1 );
            mfdata->ctab.bonds = (MOL_FMT_BOND *)
                inchi_calloc( n_alloc_bonds, sizeof( MOL_FMT_BOND ) );
            if (!mfdata->ctab.bonds)
            {
                /* can't allocate bonds structure */
                retcode = 3;
                TREAT_ERR_AND_FIN( retcode, 3, err_fin, "Out of RAM" );
            }
        }
    }

    prevcode = retcode;
    if (!data_ended)
    {
        if (!pCtab->v3000)
        {
            retcode = MolfileReadBondsBlock( pCtab, inp_file, retcode, pStrErr );
        }
        else
        {
            retcode = MolfileV3000ReadBondsBlock( pCtab, inp_file, retcode, pStrErr );
        }
        if (retcode)
        {
            if (retcode < 0)
            {
                retcode = -retcode;
                data_ended = 1;
            }
            retcode = prevcode ? prevcode
                : retcode + 40;
        }


        /* SGroup, 3D, link line(s), collections, END CTAB */
        if (pCtab->v3000)
        {
            retcode = MolfileV3000ReadTailOfCTAB( pCtab, inp_file, retcode, pStrErr );
            if (retcode)
            {
                if (retcode < 0)
                {
                    retcode = -retcode;
                    data_ended = 1;
                }
                retcode = prevcode ? prevcode : retcode + 70;
            }
        }
    }
    prevcode = retcode;

    /* SText */
    if (!data_ended)
    {
        retcode = MolfileReadSTextBlock( pCtab, inp_file, retcode, pStrErr );
        if (retcode)
        {
            retcode = prevcode ? prevcode : retcode + 50;
        }
    }
    prevcode = retcode;

    /* Prop block */
    if (!data_ended)
    {
        retcode = MolfileReadPropBlock( pCtab, pHdr, inp_file, treat_polymers, retcode, pStrErr, bNoWarnings );

        if (retcode)
        {
            if (retcode < 0)
            {
                retcode = -retcode;
                data_ended = 1;
            }
            retcode = prevcode ? prevcode : retcode + 60;
        }
    }

    /* Check that all valences are in allowed range ( <=MAXVAL, currently 20 )          */
    if (1)
#ifdef TARGET_LIB_FOR_WINCHI
        if (pCtab && pCtab->atoms)
#endif
    {

        int i;
        for (i = 0; i < pCtab->n_atoms; i++)
        {
            if (pCtab->atoms) /* djb-rwth: fixing a NULL pointer dereference */
            {
                if (pCtab->atoms[i].valence > MAXVAL)
                {
                    retcode = 70 + 9;
                    TREAT_ERR( retcode, 0, "Too large input atomic valence" );
                    break;
                }
            }
        }
#if ( FIX_CURE53_ISSUE_NULL_DEREFERENCE_MAKE_A_COPY_OF_T_GROUP_INFO==1 || defined(FIX_IMPOSSIBLE_H_ISOTOPE_BUG) )
        /* Do not eat H isotopes other than [H,D,T] */
        for (i = 0; i < pCtab->n_atoms; i++)
        {
            if (pCtab->atoms) /* djb-rwth: fixing a NULL pointer dereference */
            {
                int dmass = pCtab->atoms[i].mass_difference;
                if ((!strcmp(pCtab->atoms[i].symbol, "H") && dmass != 0 && dmass != 1 && dmass != 2 && dmass != 127) ||
                    (!strcmp(pCtab->atoms[i].symbol, "D") && dmass != 0 && dmass != 1 && dmass != -1) ||
                    (!strcmp(pCtab->atoms[i].symbol, "T") && dmass != 0 && dmass != -1 && dmass != -2))
                {
                    retcode = 70 + 8;
                    TREAT_ERR(retcode, 0, "Unacceptable isotope of hydrogen");
                    break;
                }
            }
        }
#endif
    }

err_fin:
    *err = data_ended ? -retcode : retcode;

    if (should_read_all)
    {
        if (retcode)
        {
            mfdata = FreeMolfileData( mfdata );        /* delete all results */
        }
        return mfdata;
    }
    else
    {
        if (retcode)
        {
            return NULL;
        }
        else
        {
            return (MOL_FMT_DATA*) OnlyHeaderBlock;
        }
    }
}


/****************************************************************************
 Read Molfile header
****************************************************************************/
int MolfileReadHeaderLines( MOL_FMT_HEADER_BLOCK *hdr,
                            INCHI_IOSTREAM *inp_file,
                            char *pStrErr )
{
/* All input lines can have up to 80 characters */
/* Header Block */

    char line[MOL_FMT_INPLINELEN]; /* + cr +lf +zero termination + reserve */
    int  err = 0, len; /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    const int  line_len = sizeof( line );
    char *p;


        /* Header line #1: name */
    p = inchi_fgetsLf( line, line_len, inp_file );
    if (!p)
    {
        /* can't read the input file line */
        err = 1;
        /* AddErrorMessage( pStrErr, "Can't read header block name line" ); */
        goto err_fin;
    }

    remove_one_lf( line );

    /* -- Disabled to relax strictness: allow > 80 chars names. */
    /*
    if ( line[MOL_FMT_MAXLINELEN] )
    {
        //err = 2;             too long line
        goto err_fin;
    }
    */

    len = MolfileReadField( hdr->molname,
                            sizeof( hdr->molname ) - 1,
                            MOL_FMT_STRING_DATA,
                            &p ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */

    /* Header line #2 */
    p = inchi_fgetsLf( line, line_len, inp_file );
    if (!p)
    {
        /* can't read the input file line */
        err = 3;
        /* AddErrorMessage( pStrErr, "Can't read header block line 2" ); */
        goto err_fin;
    }

    remove_one_lf( line );

    /* -- Disabled to relax strictness: allow > 80 chars names. */
    /*
    if ( line[MOL_FMT_MAXLINELEN] )
    {
        err = 4;             too long input file line
        goto err_fin;
    }
    */

    len = MolfileReadField( hdr->user_initls,
                            sizeof( hdr->user_initls ) - 1,
                            MOL_FMT_STRING_DATA,
                            &p ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    len = MolfileReadField( hdr->prog_name,
                            sizeof( hdr->prog_name ) - 1,
                            MOL_FMT_STRING_DATA,
                            &p ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */

    /*------------ Relax strictness -----------------------*/
    len = MolfileReadField( &hdr->month, 2, MOL_FMT_CHAR_INT_DATA, &p ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    len = MolfileReadField( &hdr->day, 2, MOL_FMT_CHAR_INT_DATA, &p ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    len = MolfileReadField( &hdr->year, 2, MOL_FMT_CHAR_INT_DATA, &p ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    len = MolfileReadField( &hdr->hour, 2, MOL_FMT_CHAR_INT_DATA, &p ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    len = MolfileReadField( &hdr->minute, 2, MOL_FMT_CHAR_INT_DATA, &p ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    len = MolfileReadField( hdr->dim_code, sizeof( hdr->dim_code ) - 1, MOL_FMT_STRING_DATA, &p ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    len = MolfileReadField( &hdr->scaling_factor1, 2, MOL_FMT_SHORT_INT_DATA, &p ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    len = MolfileReadField( &hdr->scaling_factor2, 10, MOL_FMT_DOUBLE_DATA, &p ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    len = MolfileReadField( &hdr->energy, 12, MOL_FMT_DOUBLE_DATA, &p ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    len = MolfileReadField( &hdr->internal_regno, 6, MOL_FMT_LONG_INT_DATA, &p ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */

    /* Save the whole line 2 */
    p = line;
    len = MolfileReadField( hdr->line2,
                            sizeof( hdr->line2 ) - 1,
                            MOL_FMT_STRING_DATA,
                            &p ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */

    /* Header line #3: comment */
    p = inchi_fgetsLf( line, line_len, inp_file );

    if (!p)
    {
        err = 7;             /* can't read the line */
        /* AddErrorMessage( pStrErr, "Can't read header block comment line" ); */
        goto err_fin;
    }
    remove_one_lf( line );
    /* -- Disabled to relax strictness: allow > 80 chars comments.
    if ( line[MOL_FMT_MAXLINELEN] ){
        err = 8;             too long line
        goto err_fin;
    }
    */
    len = MolfileReadField( hdr->comment,
                            sizeof( hdr->comment ) - 1,
                            MOL_FMT_STRING_DATA,
                            &p ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */

err_fin:

    return err;
}


/****************************************************************************
 Read counts line
****************************************************************************/
int MolfileReadCountsLine( MOL_FMT_CTAB* ctab,
                           INCHI_IOSTREAM *inp_file,
                           char *pStrErr )
{
    char *p;
    char line[MOL_FMT_INPLINELEN];
    const int line_len = sizeof( line );
    int   err = 0, len; /* djb-rwth: ignoring LLVM warning: variable used to store function return value */

    /* djb-rwth: addressing coverity ID #499502 -- TREAT_ERR properly used in all cases */

    p = inchi_fgetsLf( line, line_len, inp_file );

    if (!p)
    {
        TREAT_ERR_AND_FIN( err, 1, err_fin, "Cannot read counts line" );
        /* can't read the input file line */
    }

    remove_one_lf( line );

    if (line[MOL_FMT_MAXLINELEN])
    {
        TREAT_ERR( err, 0, "Too long counts line" );  /* too long input file line */
    }

    if (0 > MolfileReadField( &ctab->n_atoms, 3, MOL_FMT_SHORT_INT_DATA, &p ) /* V2000 only: short int */
         || 0 > MolfileReadField( &ctab->n_bonds, 3, MOL_FMT_SHORT_INT_DATA, &p ) /* V2000 only: short int */

#if ( MOL_FMT_QUERY == MOL_FMT_PRESENT )
         || 0 > MolfileReadField( &ctab->n_atom_lists, 3, MOL_FMT_SHORT_INT_DATA, &p )
#else
         || 0 > MolfileReadField( NULL, 3, MOL_FMT_JUMP_TO_RIGHT, &p )
#endif

         || 0 > MolfileReadField( NULL, /*obsolete*/      3, MOL_FMT_JUMP_TO_RIGHT, &p )
         || 0 > MolfileReadField( &ctab->chiral_flag, 3, MOL_FMT_CHAR_INT_DATA, &p )
         || 0 > MolfileReadField( &ctab->n_stext_entries, 3, MOL_FMT_SHORT_INT_DATA, &p )

#if ( MOL_FMT_CPSS == MOL_FMT_PRESENT )
         || 0 > MolfileReadField( &ctab->n_reaction_components_plus_1, 3, MOL_FMT_SHORT_INT_DATA, &p )
         || 0 > MolfileReadField( &ctab->n_reactants, 3, MOL_FMT_SHORT_INT_DATA, &p )
         || 0 > MolfileReadField( &ctab->n_products, 3, MOL_FMT_SHORT_INT_DATA, &p )
         || 0 > MolfileReadField( &ctab->n_intermediates, 3, MOL_FMT_SHORT_INT_DATA, &p )
#else
         || 0 > MolfileReadField( NULL, 3, MOL_FMT_JUMP_TO_RIGHT, &p )
         || 0 > MolfileReadField( NULL, 3, MOL_FMT_JUMP_TO_RIGHT, &p )
         || 0 > MolfileReadField( NULL, 3, MOL_FMT_JUMP_TO_RIGHT, &p )
         || 0 > MolfileReadField( NULL, 3, MOL_FMT_JUMP_TO_RIGHT, &p )
#endif

         || 0 > MolfileReadField( &ctab->n_property_lines, 3, MOL_FMT_SHORT_INT_DATA, &p ))
    {
        err = 3;  /* can't interpret counts line */
        TREAT_ERR( err, 3, "Cannot interpret counts line:" );  /* too long input file line */
        dotify_non_printable_chars( line );
        AddErrorMessage( pStrErr, line );
        goto err_fin;
    }

    /* Get CTFile version (V2000 or other) */
    len = MolfileReadField( ctab->version_string,
                            sizeof( ctab->version_string ) - 1,
                            MOL_FMT_STRING_DATA,
                            &p ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */

    /* Allocate additional space if V3000 */
    if (!strcmp( ctab->version_string, "V3000" ))
    {
        ctab->v3000 = (MOL_FMT_v3000*) inchi_calloc( 1, sizeof( MOL_FMT_v3000 ) );

        if (!ctab->v3000)
        {
            AddErrorMessage( pStrErr, "Out of RAM" );
            return -1;
        }
    }
    else
        ctab->v3000 = NULL; /* paranoia */

    /* Polymer Sgroups */
    MolFmtSgroups_Alloc( &( ctab->sgroups ), 1 );

err_fin:

    return err;
}


/****************************************************************************
 Read V2000 atomic block
****************************************************************************/
int MolfileReadAtomsBlock( MOL_FMT_CTAB* ctab,
                           INCHI_IOSTREAM *inp_file,
                           int err,
                           char *pStrErr )
{
    char *p;
    char line[MOL_FMT_INPLINELEN];
    const int line_len = sizeof( line );
    int i;
    S_SHORT chg;
    static const S_SHORT charge_val[] = { 0, 3, 2, 1, 'R', -1, -2, -3 };

    /* djb-rwth: addressing coverity ID #499580 -- TREAT_ERR properly used in all cases */
    for (i = 0; i < ctab->n_atoms; i++)
    {
        p = inchi_fgetsLf( line, line_len, inp_file );

        if (!p)
        {
            if (!err)
            {
                TREAT_ERR( err, 2, "Cannot read atom block line" );
            }
            break;
        }

        remove_one_lf( line );


        if (line[MOL_FMT_MAXLINELEN])
        {
            TREAT_ERR( err, 0, "Too long atom block line" );
        }
        if (err)
        {
            if (!strcmp( line, SD_FMT_END_OF_DATA ))
            {
                err = -abs( err );
                break;
            }
            continue; /* bypass the rest of the Atom block */
        }

        if (NULL != ctab->coords)
        {
            mystrncpy( ctab->coords[i], p, 31 ); /* original coordinates */
        }

        if (NULL != ctab->atoms)
        {
            if (0 > MolfileReadField( &ctab->atoms[i].fx, 10, MOL_FMT_DOUBLE_DATA, &p )
                || 0 > MolfileReadField( &ctab->atoms[i].fy, 10, MOL_FMT_DOUBLE_DATA, &p )
                || 0 > MolfileReadField( &ctab->atoms[i].fz, 10, MOL_FMT_DOUBLE_DATA, &p )
                || 0 > MolfileReadField( NULL, /* undescribed in article*/    1, MOL_FMT_JUMP_TO_RIGHT, &p )
                || 0 == MolfileReadField( &ctab->atoms[i].symbol, 3, MOL_FMT_STRING_DATA, &p ) /* was sizeof(ctab->atoms[0].symbol)-1 */
                || 0 > MolfileReadField( &ctab->atoms[i].mass_difference, 2, MOL_FMT_CHAR_INT_DATA, &p )
                || 0 > MolfileReadField( &ctab->atoms[i].charge, 3, MOL_FMT_CHAR_INT_DATA, &p )
                || 0 > MolfileReadField( &ctab->atoms[i].stereo_parity, 3, MOL_FMT_CHAR_INT_DATA, &p )

#if ( MOL_FMT_QUERY == MOL_FMT_PRESENT )
                || 0 > MolfileReadField( &ctab->atoms[i].H_count_plus_1, 3, MOL_FMT_CHAR_INT_DATA, &p )
                || 0 > MolfileReadField( &ctab->atoms[i].stereo_care, 3, MOL_FMT_CHAR_INT_DATA, &p )
#else
                || 0 > MolfileReadField( NULL, 3, MOL_FMT_JUMP_TO_RIGHT, &p )
                || 0 > MolfileReadField( NULL, 3, MOL_FMT_JUMP_TO_RIGHT, &p )
#endif

                || 0 > MolfileReadField( &ctab->atoms[i].valence, 3, MOL_FMT_CHAR_INT_DATA, &p ))
            {

                err = 4;
                TREAT_ERR( err, 4, "Cannot interpret atom block line:" );
                dotify_non_printable_chars( line );
                AddErrorMessage( pStrErr, line );

                if (!strcmp( line, SD_FMT_END_OF_DATA ))
                {
                    err = -abs( err );
                    break;
                }
                continue; /* can't interpret a first half of atom block line */
            }


            if (2 == strlen( ctab->atoms[i].symbol ) && isupper( UCINT ctab->atoms[i].symbol[1] ))
            {
                ctab->atoms[i].symbol[1] = (char) tolower( UCINT ctab->atoms[i].symbol[1] ); /* 5-4-99 DCh*/
            }

            if (( chg = (S_SHORT) ctab->atoms[i].charge ) < 0 || chg >= (int) ( sizeof( charge_val ) / sizeof( charge_val[0] ) ))
            {
                /* ctab->atoms[i].charge = 0; */ /* error; ignore for now */
                ctab->atoms[i].charge = (S_CHAR) ( 4 - chg ); /*  allow greater charges to accommodate NCI structures. 8-20-2002 */
                ctab->atoms[i].radical = 0;
            }
            else if ('R' == ( chg = charge_val[chg] ))
            {
                ctab->atoms[i].charge = 0;
                ctab->atoms[i].radical = RADICAL_DOUBLET;
            }
            else
            {
                ctab->atoms[i].charge = (S_CHAR) chg; /* actual charge value */
                ctab->atoms[i].radical = 0;
            }

            if (

#if ( MOL_FMT_CPSS == MOL_FMT_PRESENT )
                   0 > MolfileReadField( &ctab->atoms[i].H0_designator, 3, MOL_FMT_CHAR_INT_DATA, &p )
                || 0 > MolfileReadField( &ctab->atoms[i].reaction_component_type, 3, MOL_FMT_CHAR_INT_DATA, &p )
                || 0 > MolfileReadField( &ctab->atoms[i].reaction_component_num, 3, MOL_FMT_CHAR_INT_DATA, &p )
#else
                   0 > MolfileReadField( NULL, 3, MOL_FMT_JUMP_TO_RIGHT, &p )
                || 0 > MolfileReadField( NULL, 3, MOL_FMT_JUMP_TO_RIGHT, &p )
                || 0 > MolfileReadField( NULL, 3, MOL_FMT_JUMP_TO_RIGHT, &p )
#endif

#if ( MOL_FMT_REACT == MOL_FMT_PRESENT )
                || 0 > MolfileReadField( &ctab->atoms[i].atom_atom_mapping_num, 3, MOL_FMT_SHORT_INT_DATA, &p )
                || 0 > MolfileReadField( &ctab->atoms[i].reaction_component_type, 3, MOL_FMT_CHAR_INT_DATA, &p )
#else
                || 0 > MolfileReadField( NULL, 3, MOL_FMT_JUMP_TO_RIGHT, &p )
                || 0 > MolfileReadField( NULL, 3, MOL_FMT_JUMP_TO_RIGHT, &p )
#endif

#if ( MOL_FMT_REACT == MOL_FMT_PRESENT || MOL_FMT_QUERY == MOL_FMT_PRESENT )
                || 0 > MolfileReadField( &ctab->atoms[i].exact_change_flag, 3, MOL_FMT_CHAR_INT_DATA, &p )
#else
                || 0 > MolfileReadField( NULL, 3, MOL_FMT_JUMP_TO_RIGHT, &p )
#endif

            )
            {
                err = 5; /* can't interpret a second half of atom block line */

                TREAT_ERR( err, 5, "Cannot interpret atom block line:" );
                dotify_non_printable_chars( line );
                AddErrorMessage( pStrErr, line );

                if (!strcmp( line, SD_FMT_END_OF_DATA ))
                {
                    err = -abs( err );
                    break;
                }
                continue;
            }
        }
    }

/* err_fin: */

    return err;
}


/****************************************************************************
 Read V2000 bonds block
****************************************************************************/
int MolfileReadBondsBlock( MOL_FMT_CTAB* ctab,
                           INCHI_IOSTREAM *inp_file,
                           int err,
                           char *pStrErr )
{
    char *p;
    char line[MOL_FMT_INPLINELEN];
    const int line_len = sizeof( line );
    int i;

#if 0
    if (NULL == ctab->bonds)
    {
        err = 1;
        goto err_fin;    /*internal error: memory has not been allocated for bonds structure*/
    }
#endif 

    /* djb-rwth: addressing coverity ID #499538 -- TREAT_ERR properly used in all cases */
    for (i = 0; i < ctab->n_bonds; i++)
    {
        p = inchi_fgetsLf( line, line_len, inp_file );

        if (!p)
        {
            if (!err)
            {
                TREAT_ERR( err, 2, "Cannot read bond block line" );
            }
            break;
        }

        remove_one_lf( line );

        if (line[MOL_FMT_MAXLINELEN])
        {
            err = err ? err : 3;             /* too long input file line */
        }

        if (err)
        {
            if (!strcmp( line, SD_FMT_END_OF_DATA ))
            {
                err = -abs( err );
                break;
            }
            continue;
        }

        if (ctab->bonds)
        {

            if (0 > MolfileReadField( &ctab->bonds[i].atnum1, 3, MOL_FMT_SHORT_INT_DATA, &p )
                || 0 > MolfileReadField( &ctab->bonds[i].atnum2, 3, MOL_FMT_SHORT_INT_DATA, &p )
                || 0 > MolfileReadField( &ctab->bonds[i].bond_type, 3, MOL_FMT_CHAR_INT_DATA, &p )
                || 0 > MolfileReadField( &ctab->bonds[i].bond_stereo, 3, MOL_FMT_CHAR_INT_DATA, &p )

#if ( MOL_FMT_QUERY == MOL_FMT_PRESENT )
                || 0 > MolfileReadField( &ctab->bonds[i].cBondTopology, 3, MOL_FMT_CHAR_INT_DATA, &p ) /* ring/chain */
#else
                || 0 > MolfileReadField( NULL, 3, MOL_FMT_JUMP_TO_RIGHT, &p )
#endif

#if ( MOL_FMT_REACT == MOL_FMT_PRESENT )
                || 0 > MolfileReadField( &ctab->bonds[i].cReactingCenterStatus, 3, MOL_FMT_CHAR_INT_DATA, &p )
#else
                || 0 > MolfileReadField( NULL, 3, MOL_FMT_JUMP_TO_RIGHT, &p )
#endif

                )
            {

                if (!err)
                {
                    /* can't interpret bonds block line */
                    TREAT_ERR( err, 4, "Cannot interpret bond block line:" );
                    dotify_non_printable_chars( line );
                    AddErrorMessage( pStrErr, line );
                }
                if (!strcmp( line, SD_FMT_END_OF_DATA ))
                {
                    err = -abs( err );
                    break;
                }
            }
        }
    }

/* err_fin: */

    return err;
}


/****************************************************************************
 Read SText
****************************************************************************/
int MolfileReadSTextBlock( MOL_FMT_CTAB* ctab,
                           INCHI_IOSTREAM *inp_file,
                           int err,
                           char *pStrErr )
{
    /* just pass by all stext enties without attemp to interpret them */
    char *p;
    char line[MOL_FMT_INPLINELEN];
    const int line_len = sizeof( line );
    S_SHORT i;

    for (i = 0; i < 2 * ctab->n_stext_entries; i++)
    {
        p = inchi_fgetsLf( line, line_len, inp_file );
        if (!p)
        {
            if (!err)
            {
                TREAT_ERR_AND_FIN( err, 2, err_fin, "Cannot read STEXT block line" ); /* djb-rwth: addressing coverity ID #499517 -- TREAT_ERR_AND_FIN properly used */
            }
            break;
            /* can't read the input file line */
        }

        /*
        remove_one_lf( line );
        if ( line[MOL_FMT_MAXLINELEN] ){
            TREAT_ERR( err, 2, "Warning: Too long STEXT block line");
            too long input file line
        }
        */
    }

err_fin:

    return err;
}



/****************************************************************************
 Read properties block
****************************************************************************/
int MolfileReadPropBlock( MOL_FMT_CTAB* ctab,
                          MOL_FMT_HEADER_BLOCK *pHdr,
                          INCHI_IOSTREAM *inp_file,
                          int treat_polymers,
                          int err,
                          char *pStrErr,
                          int bNoWarnings )
{
    enum { MULTI_LINE_MODE_NO_MODE, MULTI_LINE_MODE_ISIS_ALIAS };
    char *p;
    char line[MOL_FMT_INPLINELEN];
    const int line_len = sizeof( line );
    int   nMultiLineMode = MULTI_LINE_MODE_NO_MODE, nAtomNumber = 0;
    S_SHORT i, j;
    char  charM[2];
    char  szBlank[3];
    char  szType[4];
    S_SHORT  skip_lines = 0;
    S_SHORT  num_entries;
    S_SHORT  num_atoms = ctab->n_atoms;

    int  charge_encountered = 0;
    int  radical_encountered = 0;
    int  isotope_encountered = 0;
    int     polymer_occurred = 0;

    szType[0] = '\0'; /* djb-rwth: adding zero termination */
    int debug_polymers = 0;
#if ( DEBUG_POLYMERS == 1 )
    debug_polymers = 1;
#elif  ( DEBUG_POLYMERS == 2 )
    debug_polymers = 2;
#endif


    for (i = 0; ctab->version_string[0] ? 1 : ( i < ctab->n_property_lines ); i++)
    {
        /* the last line should be M END */

        /* ctab->version_string[0] == 0:
              exactly ctab->n_property_lines lines including M END */
        /* ctab->version_string[0] != 0:
              read until M END line was encountered */

        p = inchi_fgetsLf( line, line_len, inp_file );
        /* djb-rwth: addressing coverity ID #499577 -- TREAT_ERR properly used in all cases */
        if (!p)
        {
            if (!err)
            {
                TREAT_ERR( err, 2, "Cannot read properties block line" );
            }
            goto err_fin;
        }

        remove_one_lf( line );

        if (line[MOL_FMT_MAXLINELEN])
        {
            TREAT_ERR( err, 3, "Too long properties block line" );
            continue;
        }

        if (skip_lines > 0)
        {
            skip_lines--;
            continue;
        }

        /* alias */
        if (nMultiLineMode == MULTI_LINE_MODE_ISIS_ALIAS && nAtomNumber)
        {
            int  len;

            nMultiLineMode = MULTI_LINE_MODE_NO_MODE;
            if (0 >= ( len = normalize_string( p ) ))
            {
                nAtomNumber = 0;
                continue;
            }

            if (0 < len && len < (int) ( sizeof( ctab->atoms->symbol ) ))
            {
                int  nCharge, nRad;

                MOL_FMT_ATOM*  atom = ctab->atoms + nAtomNumber - 1;
                /* ctab->atoms[nAtomNumber-1].atom_aliased_flag = 1; */
                /*  extract radicals & charges */

                extract_charges_and_radicals( p, &nRad, &nCharge );

                /*  Aliased atom cannot have charge, radical & mass difference */
                /*  in the atom table or "M  CHG", "M  RAD", "M  ISO" */
                /* if ( nCharge ) */
                atom->charge = (S_CHAR) nCharge;
                /* if ( nRad ) */
                atom->radical = (char) nRad;

                if (1 == len && 'D' == p[0])
                {
                    /*  H isotope */
                    p[0] = 'H';
                    atom->mass_difference = 1;
                }
                else
                {
                    if (1 == len && 'T' == p[0])
                    {
                        /*  H isotope */

                        p[0] = 'H';
                        atom->mass_difference = 2;
                    }
                    else
                    {
                        atom->mass_difference = 0;
                    }
                }
                if (strlen( p ) < sizeof( ctab->atoms[0].symbol ))
                {
                    strcpy(atom->symbol, p);
                }
                else
                {
                    strcpy(atom->symbol, "???");
                }
                atom->atom_aliased_flag++;
            }
            /* else if( 0 < len  )
            {
               ^^^ Just too long alias name.
                     For consistency with parsing {H,D,T}-containing alias names, this
                     would result in issuing error rather than ignoring... However,
                     for compatibility reasons, the 'ignore' behavior remained intact.
            }
            */

            skip_lines = 0;
            nAtomNumber = 0;
            continue;
        }

        if (1 != MolfileReadField( charM, sizeof( charM ) - 1, MOL_FMT_STRING_DATA, &p )
            || 0 != MolfileReadField( szBlank, sizeof( szBlank ) - 1, MOL_FMT_STRING_DATA, &p ) /* must contain 0 bytes */
            || 0 >= MolfileReadField( szType, sizeof( szType ) - 1, MOL_FMT_STRING_DATA, &p ) /* must contain 3 bytes */
        )
        {
            if (!strcmp( line, SD_FMT_END_OF_DATA ))
            {
                err = err ? -abs( err ) : -4;
                break;
            }
            continue;  /* ignore because cannot recognize */
        }

        if (charM[0] == 'V')
        {
            skip_lines = 0;   /* ISIS/Desktop Atom Value: one-line property */
            continue;
        }

        if (charM[0] == 'G')
        {
            skip_lines = 1;   /* ISIS/Desktop Group abbreviation: two-line property */
            continue;
        }

        if (charM[0] == 'A')
        {
            if (NULL != ctab->atoms &&
                 0 < ( nAtomNumber = (int) strtol( szType, NULL, 10 ) ) &&
                 nAtomNumber <= ctab->n_atoms)
            {
                /* Atom Alias [ISIS/Desktop] two-line property */
                nMultiLineMode = MULTI_LINE_MODE_ISIS_ALIAS;
                continue;
            }
            else
            {
                nAtomNumber = 0;
                skip_lines = 1;
                continue;
            }
        }

        if (charM[0] == 'S' && !strcmp( szType, "SKP" ))
        {  /* skip lines */
            if (0 >= MolfileReadField( &skip_lines, 3, MOL_FMT_SHORT_INT_DATA, &p ))
            {
                skip_lines = 0;
            }
            continue;
        }

        if (charM[0] != 'M')
        {
            /* cannot recognize a line */
            continue;
        }

        if (!strcmp( szType, "REG" ))
        {
            int len;
            p = p + strspn( p, " " );
            len = strcspn( p, " " );
            len = inchi_min( len, MOL_FMT_MAX_VALUE_LEN );
            MolfileReadField( &pHdr->internal_regno, len, MOL_FMT_LONG_INT_DATA, &p );
            continue;
        }

        if (!strcmp( szType, "END" ))
        {
            if (ctab->version_string[0])
            {
                break;  /* end of property lines */
            }
            continue;
        }

        if (NULL == ctab->atoms)
        {
            continue; /* ignore because the user requested to bypass all this stuff */
        }

        /* Charge: Generic */
        if (!strcmp( szType, "CHG" ) &&
             0 < MolfileReadField( &num_entries, 3, MOL_FMT_SHORT_INT_DATA, &p ) &&
             1 <= num_entries && num_entries <= 8)
        {

            S_SHORT atoms[8];
            S_SHORT charges[8];

            if (!charge_encountered && !radical_encountered)
            {
                /* first charge or radical record clears all Atom Block */
                /* entered charge and radical data to zeroes            */
                charge_encountered = -1;
            }
            for (j = 0; j < num_entries; j++)
            {
                if (0 > MolfileReadField( &atoms[j], 0, MOL_FMT_SHORT_INT_DATA, &p ) ||
                     0 > MolfileReadField( &charges[j], 0, MOL_FMT_SHORT_INT_DATA, &p ) ||
                     atoms[j] <= 0 || atoms[j] > num_atoms ||
                     charges[j] < -15 || charges[j]  > 15)
                {
                    goto charge_error;
                }
            }
            if (charge_encountered == -1)
            {
                for (j = 0; j < num_atoms; j++)
                {
                    if (!ctab->atoms[j].atom_aliased_flag) /* do not clear aliased atoms.*/
                    {
                        ctab->atoms[j].charge = ctab->atoms[j].radical = '\0';
                    }
                }
                charge_encountered = 1;
            }
            for (j = 0; j < num_entries; j++)
            {
                if (!ctab->atoms[atoms[j] - 1].atom_aliased_flag) /* do not change aliased atoms.*/
                {
                    ctab->atoms[atoms[j] - 1].charge = (S_CHAR) charges[j];
                }
            }
            continue;
        charge_error:
            TREAT_ERR( err, 0, "Charge not recognized:" );
            dotify_non_printable_chars( line );
            AddErrorMessage( pStrErr, line );
            continue; /* ignore for now */
        }

        /* Radical: Generic */
        if (!strcmp( szType, "RAD" ) &&
             0 < MolfileReadField( &num_entries, 3, MOL_FMT_SHORT_INT_DATA, &p ) &&
             1 <= num_entries && num_entries <= 8)
        {

            S_SHORT atoms[8];
            S_SHORT radicals[8];

            if (!charge_encountered && !radical_encountered)
            {
                /* first charge or radical record clears all Atom Block */
                /* entered charge and radical data to zeroes            */
                radical_encountered = -1;
            }
            for (j = 0; j < num_entries; j++)
            {
                if (0 > MolfileReadField( &atoms[j], 0, MOL_FMT_SHORT_INT_DATA, &p ) ||
                     0 > MolfileReadField( &radicals[j], 0, MOL_FMT_SHORT_INT_DATA, &p ) ||
                     atoms[j] <= 0 || atoms[j] > num_atoms ||
                     radicals[j] < 0 || radicals[j]  > 3)
                {
                    goto radical_error;
                }
            }
            if (radical_encountered == -1)
            {
                for (j = 0; j < num_atoms; j++)
                {
                    if (!ctab->atoms[j].atom_aliased_flag)  /* do not clear aliased atoms. 5-3-99 DCh */
                        ctab->atoms[j].charge = ctab->atoms[j].radical = '\0';
                }
                radical_encountered = 1;
            }
            for (j = 0; j < num_entries; j++)
            {
                if (!ctab->atoms[atoms[j] - 1].atom_aliased_flag)
                {
                    /* do not change aliased atoms. 5-3-99 DCh */
                    ctab->atoms[atoms[j] - 1].radical = (S_CHAR) radicals[j];
                }
            }
            continue;
        radical_error:
            TREAT_ERR( err, 0, "Radical not recognized:" );
            dotify_non_printable_chars( line );
            AddErrorMessage( pStrErr, line );
            continue; /* ignore error for now */
        }

        /* Isotope: Generic */
        if (!strcmp( szType, "ISO" ) &&
             0 < MolfileReadField( &num_entries, 3, MOL_FMT_SHORT_INT_DATA, &p ) &&
             1 <= num_entries && num_entries <= 8)
        {

            S_SHORT atoms[8];
            S_SHORT iso_mass[8]; /*  contains istotope mass number, not difference. 7-14-00 DCh. */

            if (!isotope_encountered)
            {
                /* first charge or radical record clears all Atom Block */
                /* entered charge and radical data to zeroes            */
                isotope_encountered = -1;
            }
            for (j = 0; j < num_entries; j++)
            {
                if (0 > MolfileReadField( &atoms[j], 0, MOL_FMT_SHORT_INT_DATA, &p ) ||
                     0 > MolfileReadField( &iso_mass[j], 0, MOL_FMT_SHORT_INT_DATA, &p ) ||
                     atoms[j] <= 0 || atoms[j] > num_atoms
                     /*|| iso_mass[j] < -18 || iso_mass[j]  > 12*/)
                {
                    /* goto isotope_error; */
                    atoms[j] = -1; /*  flag error */
                    TREAT_ERR( err, 0, "Isotopic data not recognized:" );
                    dotify_non_printable_chars( line );
                    AddErrorMessage( pStrErr, line );
                    continue; /* ignore isotopic error for now */
                }
            }

            if (isotope_encountered == -1)
            {
                for (j = 0; j < num_atoms; j++)
                {
                    /*if ( !ctab->atoms[j].atom_aliased_flag )*/  /* clear even aliased atoms */
                    ctab->atoms[j].mass_difference = 0;
                }
                isotope_encountered = 1;
            }
            for (j = 0; j < num_entries; j++)
            {
                if (atoms[j] <= 0)
                {
                    continue; /* ignore isotopic error for now */
                }

                if (1 /* !ctab->atoms[atoms[j]-1].atom_aliased_flag */)
                {
                    char *at = ctab->atoms[atoms[j] - 1].symbol;
                    if (at[1] || (at[0] != 'D' && at[0] != 'T')) /* djb-rwth: addressing LLVM warning */
                    {  /*  D & T cannot have ISO */
                        /*  need atomic weight to calculate isotope difference. 7-14-00 DCh. */

                        int  atw, atw_diff;
                        /*
                        NB: According to CTFile specification, difference should be in
                        [-18; +12] range, not in [-19; +19] as is checked below. */
                        if (( atw = get_atomic_mass( at ) ) && abs( atw_diff = (int) iso_mass[j] - atw ) < 20)
                        {
                            ctab->atoms[atoms[j] - 1].mass_difference = (char) ( atw_diff ? atw_diff : ZERO_ATW_DIFF );
                        }
                    }
                }
            }
            continue;
        }

        /* Sgroup, polymeric */
        if ( !strcmp( szType, "STY" ) ||
             !strcmp( szType, "SST" ) ||
             !strcmp( szType, "SLB" ) ||
             !strcmp( szType, "SCN" ) ||
             !strcmp( szType, "SAL" ) ||
             !strcmp( szType, "SBL" ) ||
             !strcmp( szType, "SDI" ) ||
             !strcmp( szType, "SMT" ) ||
             !strcmp( szType, "SBT" ))
        {
            int result;
            if (!treat_polymers)
            {
                polymer_occurred = 1;
                continue;
            }
            result = MolfileReadSgroupOfPolymer( ctab, pHdr, inp_file, line, szType, p, err, pStrErr );
            if (result != 0)
            {
                TREAT_ERR( err, result, "Could not interpret Molfile polymer data:" );
                dotify_non_printable_chars( line );
                AddErrorMessage( pStrErr, line );
                continue;
            }
        }
    }

err_fin:
    if (!treat_polymers && polymer_occurred)
    {
        /* for compatibility reasons, inchi-1 by default
        ignores polymer related lines (as v. 1.04 did)    */
        if (!bNoWarnings)
        {
            WarningMessage( pStrErr, "Ignore polymer data" );
        }
    }

    if (( ctab->sgroups.used > 0 ) && ( debug_polymers > 1 ))
    {
        ITRACE_( "\n* THE FOLLOWING %-d POLYMER SGROUP(S) WERE RECOGNISED *\n", ctab->sgroups.used );
        for (i = 0; i < ctab->sgroups.used; i++)
        {
            char *sty[] = { "NON", "SRU", "MON", "COP", "MOD", "XL", "MER" };
            char *sst[] = { "NON", "ALT", "RAN", "BLK" };
            char *con[] = { "NON", "HT", "HH", "EU" };

            ITRACE_( "\n* GROUP %-d\n", i );
            ITRACE_( "* \tindex=%-d\n", ctab->sgroups.group[i]->id );
            ITRACE_( "* \ttype=%-d %-s\n", ctab->sgroups.group[i]->type, sty[ctab->sgroups.group[i]->type] );
            if (ctab->sgroups.group[i]->subtype > -1)
                ITRACE_( "* \tsubtype=%-d %-s\n", ctab->sgroups.group[i]->subtype, sst[ctab->sgroups.group[i]->subtype] );
            if (ctab->sgroups.group[i]->conn)
                ITRACE_( "* \tconnection_type=%-d %-s\n", ctab->sgroups.group[i]->conn,
                                                                 con[ctab->sgroups.group[i]->conn] );
            ITRACE_( "* \tlabel=%-d\n", ctab->sgroups.group[i]->label );
            ITRACE_( "* \t%-d atoms:\t", ctab->sgroups.group[i]->alist.used );
            IntArray_DebugPrint( &( ctab->sgroups.group[i]->alist ) );
            ITRACE_( "* \t%-d bonds:\t", ctab->sgroups.group[i]->blist.used );
            IntArray_DebugPrint( &( ctab->sgroups.group[i]->blist ) );
            ITRACE_( "\n" );
        }
    }

    return err;
}


/****************************************************************************
 Parse polymer SGroups of Molfile
****************************************************************************/
int MolfileReadSgroupOfPolymer( MOL_FMT_CTAB* ctab,
                                MOL_FMT_HEADER_BLOCK *pHdr,
                                INCHI_IOSTREAM *inp_file,
                                char line[MOL_FMT_INPLINELEN],
                                char *szType,
                                char *p,
                                int err,
                                char *pStrErr )
{
    S_SHORT  num_entries;
    S_SHORT  num_atoms = ctab->n_atoms;
    S_SHORT  num_bonds = ctab->n_bonds;
    int j, index = -1, ret;
    char  stmp[4], stmplong[81];
    S_SHORT sg_nums[8], sg_num = -1, sg_atoms[15], sg_bonds[15], tmp;
    int q, fail = 0, len;
    /* djb-rwth: removing redundant variables */
#if ( DEBUG_POLYMERS == 1 )
    debug_polymers = 1;
#elif  ( DEBUG_POLYMERS == 2 )
    debug_polymers = 2;
#endif
    /* djb-rwth: removing redundant code */

    /* Check for possible lead codes */

    /*    STY - Sgroup type
                Polymer-related recognized types are:
                SRU = SRU type,
                MON = monomer,
                COP = copolymer,
                MER = Mer type,
                MOD
                CRO
    */
    if ( !strcmp( szType, "STY" )
         && 0 < MolfileReadField( &num_entries, 3, MOL_FMT_SHORT_INT_DATA, &p )
         && 1 <= num_entries && num_entries <= 8 )
    {
        for (j = 0; j < num_entries; j++)
        {
            fail = 0 > MolfileReadField( &sg_nums[j], 0, MOL_FMT_SHORT_INT_DATA, &p ) ||
                0 > MolfileReadField( stmp, 4, MOL_FMT_STRING_DATA, &p );
            if (!fail)
            {
                int type = MOL_FMT_M_STY_NON;

                lrtrim( stmp, &len );

                if (!strcmp( stmp, "SRU" ))
                {
                    type = MOL_FMT_M_STY_SRU;
                }
                else if (!strcmp( stmp, "MON" ))
                {
                    type = MOL_FMT_M_STY_MON;
                }
                else if (!strcmp( stmp, "COP" ))
                {
                    type = MOL_FMT_M_STY_COP;
                }
                else if (!strcmp( stmp, "MOD" ))
                {
                    type = MOL_FMT_M_STY_MOD;
                }
                else if (!strcmp( stmp, "CRO" ))
                {
                    type = MOL_FMT_M_STY_CRO;
                }
                else if (!strcmp( stmp, "MER" ))
                {
                    type = MOL_FMT_M_STY_MER;
                }
                else
                {
                    fail = 1;
                }
                if (!fail)
                {
                    index = MolFmtSgroups_GetIndexBySgroupId( sg_nums[j], &( ctab->sgroups ) );
                    if (-1 == index)
                    {
                        ret = MolFmtSgroups_Append( &ctab->sgroups, sg_nums[j], type );
                        if (0 != ret)
                            fail = 1;
                        else
                            index = ctab->sgroups.used - 1;
                    }
                    else
                    {
                        ctab->sgroups.group[index]->type = type;
                    }
                    if (!fail)
                    {
                        for (q = 0; q < 4; q++)
                        {
                            ctab->sgroups.group[index]->xbr1[q] = -777777.777;
                            ctab->sgroups.group[index]->xbr2[q] = -777777.777;
                        }
                        ctab->sgroups.group[index]->smt[0] = '\0';
                    }
                }
            }
            if (fail)
            {
                err = 5;
                goto err_exit;
            }
        }
    }
    /*
    SST - Polymer Sgroup subtypes:
                ALT = alternating,
                RAN = random,
                BLK = block
    */
    else if (!strcmp( szType, "SST" ) &&
              0 < MolfileReadField( &num_entries, 3, MOL_FMT_SHORT_INT_DATA, &p ) &&
              1 <= num_entries && num_entries <= 8)
    {
        for (j = 0; j < num_entries; j++)
        {
            fail = 0 > MolfileReadField( &sg_nums[j], 0, MOL_FMT_SHORT_INT_DATA, &p ) ||
                0 > MolfileReadField( stmp, 4, MOL_FMT_STRING_DATA, &p );
            if (!fail)
            {
                index = MolFmtSgroups_GetIndexBySgroupId( sg_nums[j], &( ctab->sgroups ) );
                if (-1 == index)
                {
                    fail = 1;
                }
            }
            if (!fail)
            {
                ctab->sgroups.group[index]->subtype = MOL_FMT_M_SST_NON;
                lrtrim( stmp, &len );
                if (!strcmp( stmp, "ALT" ))
                {
                    ctab->sgroups.group[index]->subtype = MOL_FMT_M_SST_ALT;
                }
                else if (!strcmp( stmp, "RAN" ))
                {
                    ctab->sgroups.group[index]->subtype = MOL_FMT_M_SST_RAN;
                }
                else if (!strcmp( stmp, "BLO" ))
                {
                    ctab->sgroups.group[index]->subtype = MOL_FMT_M_SST_BLK;
                }
                else if (!strcmp( stmp, "BLK" ))
                {
                    ctab->sgroups.group[index]->subtype = MOL_FMT_M_SST_BLK;
                }
                else
                {
                    fail = 1;
                    break;
                }
            }
        }
        if (fail)
        {
            err = 6;
            goto err_exit;
        }
    }
    /*    SLB - Sgroup Labels */
    else if (!strcmp( szType, "SLB" ) &&
              0 < MolfileReadField( &num_entries, 3, MOL_FMT_SHORT_INT_DATA, &p ) &&
              1 <= num_entries && num_entries <= 8)
    {
        for (j = 0; j < num_entries; j++)
        {
            fail = 0 > MolfileReadField( &sg_nums[j], 0, MOL_FMT_SHORT_INT_DATA, &p )
                   || 0 > MolfileReadField( &tmp, 0, MOL_FMT_SHORT_INT_DATA, &p );
            if (!fail)
            {
                index = MolFmtSgroups_GetIndexBySgroupId( sg_nums[j], &( ctab->sgroups ) );
                if (-1 == index)
                {
                    fail = 1;
                }
            }
            if (!fail)
            {
                ctab->sgroups.group[index]->label = tmp;
            }
        }
        if (fail)
        {
            err = 7;
            goto err_exit;
        }
    }
    /*    SCN - Sgroup Connectivity
                    HH = head-to-head,
                    HT = head-to-tail,
                    EU = either unknown.
                    Left justified.
    */
    else if (!strcmp( szType, "SCN" ) &&
              0 < MolfileReadField( &num_entries, 3, MOL_FMT_SHORT_INT_DATA, &p ) &&
              1 <= num_entries && num_entries <= 8)
    {
        for (j = 0; j < num_entries; j++)
        {
            fail = 0 > MolfileReadField( &sg_nums[j], 0, MOL_FMT_SHORT_INT_DATA, &p )
                   || 0 > MolfileReadField( stmp, 4, MOL_FMT_STRING_DATA, &p );
            if (!fail)
            {
                index = MolFmtSgroups_GetIndexBySgroupId( sg_nums[j], &( ctab->sgroups ) );
                if (-1 == index)
                {
                    fail = 1;
                }
            }
            if (!fail)
            {
                ctab->sgroups.group[index]->conn = MOL_FMT_M_CONN_NON;
                lrtrim( stmp, &len );
                if (!strcmp( stmp, "HT" ))
                {
                    ctab->sgroups.group[index]->conn = MOL_FMT_M_CONN_HT;
                }
                else if (!strcmp( stmp, "HH" ))
                {
                    ctab->sgroups.group[index]->conn = MOL_FMT_M_CONN_HH;
                }
                else if (!strcmp( stmp, "EU" ))
                {
                    ctab->sgroups.group[index]->conn = MOL_FMT_M_CONN_EU;
                }
                else
                {
                    fail = 1;       /* NB: we do not allow explicit different abbreviation - but note that  */
                                    /* totally skipping SCN line is allowed ("EU" will be inserted further) */
                }
            }
            if (fail)
            {
                err = 8;
                goto err_exit;
            }
        }
    }
    /* SAL - Sgroup atoms list  */
    else if (!strcmp( szType, "SAL" ))
    {
        if (0 < MolfileReadField( &sg_num, 4, MOL_FMT_SHORT_INT_DATA, &p ) &&
              0 < MolfileReadField( &num_entries, 3, MOL_FMT_SHORT_INT_DATA, &p ))
        {
            index = MolFmtSgroups_GetIndexBySgroupId( sg_num, &( ctab->sgroups ) );
            if (-1 == index)
            {
                fail = 1;
            }
            if (!fail)
            {
                for (j = 0; j < num_entries; j++)
                {
                    if (0 > MolfileReadField( &sg_atoms[j], 0, MOL_FMT_SHORT_INT_DATA, &p ) ||
                        sg_atoms[j] <= 0 || sg_atoms[j] > num_atoms)
                    {
                        fail = 1;
                        break;
                    }
                }
            }
            if (!fail)
            {
                for (j = 0; j < num_entries; j++)
                {
                    if (0 != IntArray_Append( &( ctab->sgroups.group[index]->alist ), sg_atoms[j] ))
                    {
                        fail = 1;
                        break;
                    }
                }
            }
        }
        else
        {
            fail = 1;
        }
        if (fail)
        {
            err = 9;
            goto err_exit;
        }
    }
    /*    SBL - Sgroup bonds list */
    else if (!strcmp( szType, "SBL" ))
    {
        if (0 < MolfileReadField( &sg_num, 4, MOL_FMT_SHORT_INT_DATA, &p ) &&
            0 < MolfileReadField( &num_entries, 3, MOL_FMT_SHORT_INT_DATA, &p ))
        {
            index = MolFmtSgroups_GetIndexBySgroupId( sg_num, &( ctab->sgroups ) );
            if (-1 == index)
            {
                fail = 1;
            }
            if (!fail)
            {
                for (j = 0; j < num_entries; j++)
                {
                    if (0 > MolfileReadField( &sg_bonds[j], 0, MOL_FMT_SHORT_INT_DATA, &p ) ||
                         sg_bonds[j] <= 0 || sg_bonds[j] > num_bonds)
                    {
                        fail = 1;
                        break;
                    }
                }
            }
            if (!fail)
            {
                for (j = 0; j < num_entries; j++)
                {
                    if (0 != IntArray_Append( &( ctab->sgroups.group[index]->blist ), sg_bonds[j] ))
                    {
                        fail = 1;
                        break;
                    }
                }
            }
        }
        else
        {
            fail = 1;
        }
        if (fail)
        {
            err = 10;
            goto err_exit;
        }
    }
    /* SDI */
    else if (!strcmp( szType, "SDI" ))
    {
        double x[4];

        if (0 < MolfileReadField( &sg_num, 4, MOL_FMT_SHORT_INT_DATA, &p ) &&
            0 < MolfileReadField( &num_entries, 3, MOL_FMT_SHORT_INT_DATA, &p ))
        {
            index = MolFmtSgroups_GetIndexBySgroupId( sg_num, &( ctab->sgroups ) );
            if (-1 == index)
            {
                fail = 1;
            }
            else if (num_entries != 4)
            {
                fail = 1;
            }
            if (!fail)
            {
                for (j = 0; j < num_entries; j++)
                {
                    if (0 > MolfileReadField( &x[j], 0, MOL_FMT_DOUBLE_DATA, &p ))
                    {
                        fail = 1;
                        break;
                    }
                }
            }
            if (!fail)
            {
                if (fabs( -fabs( ctab->sgroups.group[index]->xbr1[0] ) + 777777.777 ) < 1.e-7)    /* brkt1 coords not yet here */
                {
                    for (q = 0; q < 4; q++)
                    {
                        ctab->sgroups.group[index]->xbr1[q] = x[q];
                    }
                }
                else
                {
                    for (q = 0; q < 4; q++)
                    {
                        ctab->sgroups.group[index]->xbr2[q] = x[q];
                    }
                }
            }
        }
        else
        {
            fail = 1;
        }
        if (fail)
        {
            err = 11;
            goto err_exit;
        }
    }
    /* SMT - Sgroup Subscript */
    else if (!strcmp( szType, "SMT" ))
    {
        index = -1;
        if (0 < MolfileReadField( &sg_num, 4, MOL_FMT_SHORT_INT_DATA, &p ) &&
            0 < MolfileReadField( stmplong, 80, MOL_FMT_STRING_DATA, &p ))
        {
            index = MolFmtSgroups_GetIndexBySgroupId( sg_num, &( ctab->sgroups ) );
        }
        if (-1 == index)
        {
            fail = 1;
        }
        if (!fail)
        {
            lrtrim( stmplong, &len );
            strcpy(ctab->sgroups.group[index]->smt, stmplong);
        }
        if (fail)
        {
            err = 11;
            goto err_exit;
        }
    }


    ITRACE_( "\n" );
    return 0;

err_exit:
    MolFmtSgroups_Free( &ctab->sgroups );

    return err; /* ignore polymeric error for now */
}


/****************************************************************************

 Pseudoelement treatment

    If allowed:   "*" ===> "Zz"

    If disabled:   fills err and pStrErr

****************************************************************************/
static int MolfileTreatPseudoElementAtoms( MOL_FMT_CTAB* ctab,
                                           int pseudos_allowed,
                                           int *err,
                                           char *pStrErr )
{
    int i, nzz = 0;

    /* djb-rwth: addressing coverity ID #499499 -- TREAT_ERR properly used in all cases */

    for (i = 0; i < ctab->n_atoms; i++)
    {
        int is_zz = 0, is_star = 0;

        /* Zy is specifically disabled */
        if (!strcmp(ctab->atoms[i].symbol, "Zy"))
        {
            TREAT_ERR( *err, ( 70 + 6 ), "Invalid element(s):" );
            TREAT_ERR( *err, ( 70 + 6 ), ctab->atoms[i].symbol );
        }

        is_star = !strcmp( ctab->atoms[i].symbol, "*" );
        if (!is_star)
        {
            is_zz = !strcmp( ctab->atoms[i].symbol, "Zz" );
        }

        if (is_star || is_zz)
        {
            nzz++;
            if (0== pseudos_allowed)
            {
                /* Pseudoelements totally disabled */
                TREAT_ERR( *err, ( 70 + 6 ), "Invalid element(s):" );
                TREAT_ERR( *err, ( 70 + 6 ), ctab->atoms[i].symbol );
            }
            else
            {
                /* That's allowed, if it's star then convert to Zz */
                if (is_star)
                {
                    mystrncpy( ctab->atoms[i].symbol, "Zz", sizeof( "Zz" ) );
                }
            }
        }
    }

    return nzz;
}

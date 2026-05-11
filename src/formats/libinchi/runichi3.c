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


/*
    Pre-processing related functions

*/

#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <ctype.h>


#include "mode.h"
#include "ichitime.h"
#ifndef COMPILE_ANSI_ONLY
#include <conio.h>
#endif
#include "ichimain.h"
#include "ichi_io.h"
#include "mol_fmt.h"
#include "inchi_api.h"
#include "readinch.h"
#ifdef TARGET_LIB_FOR_WINCHI
#include "../../../IChI_lib/src/ichi_lib.h"
#include "inchi_api.h"
#else
#include "inchi_gui.h"
#endif
#include "readinch.h"
#include "inpdef.h"
#include "ichi_io.h"

#include "bcf_s.h"

/* Local prototypes */
static int OrigAtData_bCheckUnusualValences( ORIG_ATOM_DATA *orig_at_data,
                                             int bAddIsoH,
                                             char *pStrErrStruct,
                                             int bNoWarnings);

static void OAD_PolymerUnit_RemoveLinkFromCRUChain( int at1, int at2, int *nbonds, int **bonds );

/****************************************************************************
  Check inp_ATOM's for unusual valence
****************************************************************************/
int OrigAtData_bCheckUnusualValences( ORIG_ATOM_DATA *orig_at_data,
                                      int bAddIsoH,
                                      char *pStrErrStruct,
                                      int bNoWarnings)
{
    int i, val, num_found = 0;
    char msg[32];
    int len, num_H;

    int already_here = ( orig_at_data && orig_at_data->num_inp_atoms > 0 );

    inp_ATOM *at = already_here ? orig_at_data->at : NULL;

    if (at)
    {
        for (i = 0, num_found = 0; i < orig_at_data->num_inp_atoms; i++)
        {
            num_H = bAddIsoH ? NUMH( at, i ) : at[i].num_H;

            val = detect_unusual_el_valence( at[i].el_number,
                                             at[i].charge,
                                             at[i].radical,
                                             at[i].chem_bonds_valence,
                                             num_H,
                                             at[i].valence );
            if (val)
            {
                num_found++;
                if (!bNoWarnings)
                {
                    WarningMessage( pStrErrStruct, "Accepted unusual valence(s):" );
                }
                len = sprintf(msg, "%s", at[i].elname);
                if (at[i].charge)
                {
                    len += sprintf(msg + len, "%+d", at[i].charge);
                }
                if (at[i].radical)
                {
                    len += sprintf(msg + len, ",%s", at[i].radical == RADICAL_SINGLET ? "s" :
                        at[i].radical == RADICAL_DOUBLET ? "d" :
                        at[i].radical == RADICAL_TRIPLET ? "t" : "?");
                }
                len += sprintf(msg + len, "(%d)", val);
                if (!bNoWarnings)
                {
                    WarningMessage( pStrErrStruct, msg );
                }
            }
        }
    }

    return num_found;
}


/****************************************************************************
  Make a copy of ORIG_ATOM_DATA
****************************************************************************/
int OrigAtData_Duplicate( ORIG_ATOM_DATA *new_orig_atom,
                          ORIG_ATOM_DATA *orig_atom )
{
    inp_ATOM  *at = NULL;
    AT_NUMB   *nCurAtLen = NULL;
    AT_NUMB   *nOldCompNumber = NULL;
    int k, m, nn;
    int orig_nat = orig_atom->num_inp_atoms;

    int ret = -1; /* fail; 0 - OK */

    if (new_orig_atom->at &&
         new_orig_atom->num_inp_atoms >= orig_nat)
    {
        at = new_orig_atom->at;
    }
    else
    {
        at = (inp_ATOM *) inchi_calloc( (long long)orig_nat + 1, sizeof( at[0] ) ); /* djb-rwth: cast operator added */
        if (!at)
        {
            goto exit_function;
        }
    }

    if (new_orig_atom->nOldCompNumber &&
         new_orig_atom->num_components >= orig_atom->num_components)
    {
        nCurAtLen = new_orig_atom->nCurAtLen;
    }
    else
    {
        nCurAtLen = (AT_NUMB *) inchi_calloc( (long long)orig_atom->num_components + 1, sizeof( nCurAtLen[0] ) ); /* djb-rwth: cast operator added */
        if (!nCurAtLen)
        {
            goto exit_function;
        }
    }

    if (new_orig_atom->nCurAtLen &&
         new_orig_atom->num_components >= orig_atom->num_components)
    {
        nOldCompNumber = new_orig_atom->nOldCompNumber;
    }
    else
    {
        nOldCompNumber = (AT_NUMB *) inchi_calloc( (long long)orig_atom->num_components + 1,
                                                  sizeof( nOldCompNumber[0] ) ); /* djb-rwth: cast operator added */
        if (!nOldCompNumber)
        {
            goto exit_function;
        }
    }

    if (at && nCurAtLen && nOldCompNumber)
    {
        /* Copy */
        if (orig_atom->at)
        {
            memcpy(at, orig_atom->at,
                orig_nat * sizeof(new_orig_atom->at[0]));
        }
        if (orig_atom->nCurAtLen)
        {
            memcpy(nCurAtLen, orig_atom->nCurAtLen,
                orig_atom->num_components * sizeof(nCurAtLen[0]));
        }
        if (orig_atom->nOldCompNumber)
        {
            memcpy(nOldCompNumber, orig_atom->nOldCompNumber,
                orig_atom->num_components * sizeof(nOldCompNumber[0]));
        }

        /* Deallocate */
        if (new_orig_atom->at && new_orig_atom->at != at)
        {
            inchi_free( new_orig_atom->at );
        }
        if (new_orig_atom->nCurAtLen && new_orig_atom->nCurAtLen != nCurAtLen)
        {
            inchi_free( new_orig_atom->nCurAtLen );
        }
        if (new_orig_atom->nOldCompNumber &&
             new_orig_atom->nOldCompNumber != nOldCompNumber)
        {
            inchi_free( new_orig_atom->nOldCompNumber );
        }

        *new_orig_atom = *orig_atom;
        new_orig_atom->at = at;
        new_orig_atom->nCurAtLen = nCurAtLen;
        new_orig_atom->nOldCompNumber = nOldCompNumber;

        /* Data that are not to be copied */
        new_orig_atom->nNumEquSets = 0;
        memset( new_orig_atom->bSavedInINCHI_LIB, 0, sizeof( new_orig_atom->bSavedInINCHI_LIB ) ); /* djb-rwth: memset_s C11/Annex K variant? */
        memset( new_orig_atom->bPreprocessed, 0, sizeof( new_orig_atom->bPreprocessed ) ); /* djb-rwth: memset_s C11/Annex K variant? */



        new_orig_atom->szCoord = NULL; 
        if (orig_atom->szCoord)
        {
            new_orig_atom->szCoord = (MOL_COORD *) inchi_calloc(orig_nat, sizeof(new_orig_atom->szCoord[0]));
            if (!new_orig_atom->szCoord)
            {
                goto exit_function;
            }
            memcpy(new_orig_atom->szCoord, orig_atom->szCoord, orig_nat * sizeof(new_orig_atom->szCoord[0]));
        }
        

        /* Arrays that are not to be copied */
        
        new_orig_atom->nEquLabels = NULL;
        new_orig_atom->nSortedOrder = NULL;

        new_orig_atom->polymer = NULL;
        if (orig_atom->polymer)
        {
            /* Polymer stuff -- deep copy */
            OAD_Polymer *oldp = orig_atom->polymer;
            OAD_Polymer *newp = NULL;

            newp = (OAD_Polymer *) inchi_calloc( 1, sizeof( OAD_Polymer ) );
            if (!newp)
            {
                inchi_free(newp); /* djb-rwth: avoiding memory leak */
                goto exit_function;
            }
            memcpy(newp, orig_atom->polymer, sizeof(OAD_Polymer));
            newp->units = (OAD_PolymerUnit**) inchi_calloc( newp->n, sizeof(OAD_PolymerUnit*) ); /* djb-rwth: inchi_calloc must return OAD_PolymerUnit** */
            if (!newp->units)
            {
                inchi_free(newp); /* djb-rwth: avoiding memory leak */
                goto exit_function;
            }
            for (k = 0; k < orig_atom->polymer->n; k++)
            {
                newp->units[k] = OAD_PolymerUnit_CreateCopy( orig_atom->polymer->units[k] );
            }
            if (oldp->n_pzz > 0)
            {
                newp->n_pzz = oldp->n_pzz;
                newp->pzz = (int *) inchi_calloc( newp->n_pzz, sizeof( int ) );
                if (!newp->pzz)
                {
                    inchi_free(newp->units); /* djb-rwth: fixing coverity ID #499546 */
                    inchi_free(newp); /* djb-rwth: avoiding memory leak */
                    goto exit_function;
                }
                memcpy(newp->pzz, oldp->pzz, newp->n_pzz * sizeof(oldp->pzz[0]));
            }
            new_orig_atom->polymer = newp;
        }

        new_orig_atom->v3000 = NULL;
        if (orig_atom->v3000)
        {
            /* V3000 features -- deep copy */
            OAD_V3000 *new_v3000 = NULL;
            new_v3000 = (OAD_V3000 *) inchi_calloc( 1, sizeof( OAD_V3000 ) );
            if (!new_v3000)
            {
                inchi_free(new_v3000); /* djb-rwth: avoiding memory leak */
                goto exit_function;
            }
            memcpy(new_v3000, orig_atom->v3000, sizeof(OAD_V3000));
            if (orig_atom->v3000->atom_index_orig)
            {
                new_v3000->atom_index_orig = (int *) inchi_calloc( orig_nat, sizeof( int ) );
                /* if ( NULL==new_v3000->atom_index_orig ) {TREAT_ERR( err, 9001, "Out of RAM"); goto exit_function; } */
                if (!new_v3000->atom_index_orig)
                {
                    inchi_free(new_v3000); /* djb-rwth: avoiding memory leak */
                    goto exit_function;
                }
                memcpy(new_v3000->atom_index_orig, orig_atom->v3000->atom_index_orig, orig_nat * sizeof(int));
            }
            if (orig_atom->v3000->atom_index_fin)
            {
                new_v3000->atom_index_fin = (int *) inchi_calloc( orig_nat, sizeof( int ) );
                /* if ( NULL==new_v3000->atom_index_fin ) {TREAT_ERR( err, 9001, "Out of RAM"); goto exit_function; } */
                if (!new_v3000->atom_index_fin)
                {
                    inchi_free(new_v3000); /* djb-rwth: avoiding memory leak */
                    goto exit_function;
                }
                memcpy(new_v3000->atom_index_fin, orig_atom->v3000->atom_index_fin, orig_nat * sizeof(int));
            }
            if (orig_atom->v3000->n_haptic_bonds && orig_atom->v3000->lists_haptic_bonds)
            {
                new_v3000->lists_haptic_bonds = (int **) inchi_calloc( orig_atom->v3000->n_haptic_bonds, sizeof( int* ) );
                /* if ( NULL==new_v3000->lists_haptic_bonds ) { TREAT_ERR( err, 9001, "Out of RAM"); goto exit_function; }*/
                for (m = 0; m < orig_atom->v3000->n_haptic_bonds; m++)
                {
                    int *lst = NULL;
                    int *old_lst = orig_atom->v3000->lists_haptic_bonds[m];
                    nn = old_lst[2] + 3;
                    lst = new_v3000->lists_haptic_bonds[m] = (int *) inchi_calloc( nn, sizeof( int ) );
                    if (!lst)
                    {
                        inchi_free(new_v3000->lists_haptic_bonds); /* djb-rwth: fixing coverity ID #499504 */
                        inchi_free(new_v3000->atom_index_orig); /* djb-rwth: fixing coverity ID #499540 */
                        inchi_free(new_v3000->atom_index_fin); /* djb-rwth: fixing coverity ID #499613 */
                        inchi_free(new_v3000); /* djb-rwth: avoiding memory leak */
                        goto exit_function;
                    }
                    memcpy(lst, old_lst, nn * sizeof(int));
                }
            }
            if (orig_atom->v3000->n_steabs && orig_atom->v3000->lists_steabs)
            {
                new_v3000->lists_steabs = (int **) inchi_calloc( orig_atom->v3000->n_steabs, sizeof( int* ) );
                /* if ( NULL==new_v3000->lists_steabs ) { TREAT_ERR( err, 9001, "Out of RAM"); goto exit_function; }*/
                for (m = 0; m < orig_atom->v3000->n_steabs; m++)
                {
                    int *lst = NULL;
                    int *old_lst = orig_atom->v3000->lists_steabs[m];
                    nn = old_lst[1] + 2;
                    lst = new_v3000->lists_steabs[m] = (int *) inchi_calloc( nn, sizeof( int ) );
                    if (!lst)
                    {
                        inchi_free(new_v3000->lists_haptic_bonds); /* djb-rwth: fixing coverity ID #499504 */
                        inchi_free(new_v3000->lists_steabs); /* djb-rwth: fixing coverity ID #499543 */
                        inchi_free(new_v3000->atom_index_orig); /* djb-rwth: fixing coverity ID #499540 */
                        inchi_free(new_v3000->atom_index_fin); /* djb-rwth: fixing coverity ID #499613 */
                        inchi_free(new_v3000); /* djb-rwth: avoiding memory leak */
                        goto exit_function;
                    }
                    memcpy(lst, old_lst, nn * sizeof(int));
                }
            }
            if (orig_atom->v3000->n_sterel && orig_atom->v3000->lists_sterel)
            {
                new_v3000->lists_sterel = (int **) inchi_calloc( orig_atom->v3000->n_sterel, sizeof( int* ) );
                if (!new_v3000)
                {
                    inchi_free(new_v3000); /* djb-rwth: avoiding memory leak */
                    goto exit_function;
                }
                /* if ( NULL==new_v3000->lists_sterel ) { TREAT_ERR( err, 9001, "Out of RAM"); goto exit_function; }*/
                for (m = 0; m < orig_atom->v3000->n_sterel; m++)
                {
                    int *lst = NULL;
                    int *old_lst = orig_atom->v3000->lists_sterel[m];
                    nn = old_lst[1] + 2;
                    if (new_v3000->lists_sterel) /* djb-rwth: fixing a NULL pointer dereference */
                        lst = new_v3000->lists_sterel[m] = (int *) inchi_calloc( nn, sizeof( int ) );
                    if (!lst)
                    {
                        inchi_free(new_v3000->lists_haptic_bonds); /* djb-rwth: fixing coverity ID #499504 */
                        inchi_free(new_v3000->lists_steabs); /* djb-rwth: fixing coverity ID #499543 */
                        inchi_free(new_v3000->lists_sterel); /* djb-rwth: fixing coverity ID #499504 */
                        inchi_free(new_v3000->atom_index_orig); /* djb-rwth: fixing coverity ID #499540 */
                        inchi_free(new_v3000->atom_index_fin); /* djb-rwth: fixing coverity ID #499613 */
                        inchi_free(new_v3000); /* djb-rwth: avoiding memory leak */
                        goto exit_function;
                    }
                    memcpy(lst, old_lst, nn * sizeof(int));
                }
            }
            if (orig_atom->v3000->n_sterac && orig_atom->v3000->lists_sterac)
            {
                new_v3000->lists_sterac = (int **) inchi_calloc( orig_atom->v3000->n_sterac, sizeof( int* ) );
                /* if ( NULL==new_v3000->lists_sterac ) { TREAT_ERR( err, 9001, "Out of RAM"); goto exit_function; }*/
                if (new_v3000->lists_sterac) /* djb-rwth: fixing a NULL pointer dereference */
                {
                    for (m = 0; m < orig_atom->v3000->n_sterac; m++)
                    {
                        int* lst = NULL;
                        int* old_lst = orig_atom->v3000->lists_sterac[m];
                        nn = old_lst[1] + 2;
                        lst = new_v3000->lists_sterac[m] = (int*)inchi_calloc(nn, sizeof(int));
                        if (!lst)
                        {
                            inchi_free(new_v3000->lists_haptic_bonds); /* djb-rwth: fixing coverity ID #499504 */
                            inchi_free(new_v3000->lists_steabs); /* djb-rwth: fixing coverity ID #499543 */
                            inchi_free(new_v3000->lists_sterel); /* djb-rwth: fixing coverity ID #499504 */
                            inchi_free(new_v3000->lists_sterac); /* djb-rwth: fixing coverity ID #499575 */
                            inchi_free(new_v3000->atom_index_orig); /* djb-rwth: fixing coverity ID #499540 */
                            inchi_free(new_v3000->atom_index_fin); /* djb-rwth: fixing coverity ID #499613 */
                            inchi_free(new_v3000); /* djb-rwth: avoiding memory leak */
                            goto exit_function;
                        }
                        memcpy(lst, old_lst, nn * sizeof(int));
                    }
                }
            }

            new_orig_atom->v3000 = new_v3000;
        }

        /* Success */
        ret = 0;
    }

exit_function:

    if (ret != 0)
    {
        if (at && new_orig_atom->at != at)
            inchi_free( at );
        if (nCurAtLen && new_orig_atom->nCurAtLen != nCurAtLen)
            inchi_free( nCurAtLen );
        if (nOldCompNumber && new_orig_atom->nOldCompNumber != nOldCompNumber)
            inchi_free( nOldCompNumber );
    }

    return ret;
}


/****************************************************************************
  Preprocess the whole structure

    The plan is as follows.

    1.    Copy orig_inp_data --> prep_inp_data (then work with the latter)

    2.    Fix odd things in prep_inp_data

            Find whether the structure can be disconnected or is a salt
                - check if needs salt disconnection
                - check if needs metal disconnection

    3.    If ( orig_inp_data->bDisconnectSalts ) then
            disconnect salts in prep_inp_data

            Mark the (disconnected) components in prep_inp_data

            Detect isotopic H on heteroatoms (necessary condition
            for global isotopic tautomerism)

    4.    Detect unusual valences (should be called before metal disconnection)

    5.    Create metal-disconnected structure if applicable.
            - save reconnected structure in prep_inp_data+1 if requested
            - make Disconnected structure in prep_inp_data
****************************************************************************/
int PreprocessOneStructure( struct tagINCHI_CLOCK *ic,
                            STRUCT_DATA           *sd,
                            INPUT_PARMS           *ip,
                            ORIG_ATOM_DATA        *orig_inp_data,
                            ORIG_ATOM_DATA        *prep_inp_data )
{
    int i;
    INCHI_MODE bTautFlags = 0;
    INCHI_MODE bTautFlagsDone = 0;

    /* 1. Copy orig_inp_data --> prep_inp_data */

    if (0 > OrigAtData_Duplicate( prep_inp_data, orig_inp_data ))
    {
        AddErrorMessage( sd->pStrErrStruct, "Out of RAM" );
        sd->nStructReadError = 99;
        sd->nErrorType = _IS_FATAL;
        goto exit_function;
    }

#if ( bRELEASE_VERSION == 0 && (EXTR_HAS_METAL_ATOM & (EXTR_MASK | EXTR_FLAG) ) )
    if (bHasMetalAtom( orig_inp_data ))
    {
        sd->bExtract |= EXTR_HAS_METAL_ATOM;
    }
#endif

    /* 2. Fix odd things in prep_inp_data            */

    if (0 < fix_odd_things( prep_inp_data->num_inp_atoms, prep_inp_data->at, /*0*/ip->bTautFlags & TG_FLAG_FIX_SP3_BUG, ip->bFixNonUniformDraw ))
    {
        /* changed 2010-03-17 DT */
        if (!ip->bNoWarnings)
        {
            WarningMessage( sd->pStrErrStruct, "Charges were rearranged" );
        }
        if (sd->nErrorType < _IS_WARNING)
        {
            sd->nErrorType = _IS_WARNING;
        }
        sd->bTautFlagsDone[INCHI_BAS] |= TG_FLAG_FIX_ODD_THINGS_DONE;
    }

#if ( FIX_ADJ_RAD == 1 )
    if (ip->bTautFlags & TG_FLAG_FIX_ADJ_RADICALS)
    {
        if (0 < FixAdjacentRadicals( prep_inp_data->num_inp_atoms, prep_inp_data->at ))
        {
            sd->bTautFlagsDone[INCHI_BAS] |= TG_FLAG_FIX_ADJ_RADICALS_DONE;
        }
    }
#endif

#if ( bRELEASE_VERSION == 0 && (EXTR_FLAGS & EXTR_HAS_FEATURE) )
    if (bFoundFeature( prep_inp_data->at, prep_inp_data->num_inp_atoms ))
    {
        sd->bExtract |= EXTR_HAS_FEATURE;
    }
#endif


    /* Find whether the structure can be disconnected or is a salt */


    /* Needs salt disconnection? */

    if (ip->bTautFlags & TG_FLAG_DISCONNECT_SALTS)
    {
        prep_inp_data->bDisconnectSalts = ( 0 < DisconnectSalts( prep_inp_data, 0 ) );
    }
    else
    {
        prep_inp_data->bDisconnectSalts = 0;
    }

    /* Needs metal disconnection? */

    if (ip->bTautFlags & TG_FLAG_DISCONNECT_COORD)
    {
        i = ( 0 != ( ip->bTautFlags & TG_FLAG_CHECK_VALENCE_COORD ) );
        bMayDisconnectMetals( prep_inp_data, i, &bTautFlagsDone ); /* changes prep_inp_data->bDisconnectCoord */
        sd->bTautFlagsDone[INCHI_BAS] |= bTautFlagsDone; /* whether any disconnection has been rejected because of the metal proper valence */

#if ( bRELEASE_VERSION == 0 )
        if (i && ( bTautFlagsDone & TG_FLAG_CHECK_VALENCE_COORD_DONE ))
        {
            sd->bExtract |= EXTR_METAL_WAS_NOT_DISCONNECTED;
        }
#endif
    }
    else
    {
        prep_inp_data->bDisconnectCoord = 0;
    }
    orig_inp_data->bDisconnectSalts = prep_inp_data->bDisconnectSalts;
    orig_inp_data->bDisconnectCoord = prep_inp_data->bDisconnectCoord;

    /* 3. if( orig_inp_data->bDisconnectSalts ) then
          disconnect salts in prep_inp_data    */

    if (( ip->bTautFlags & TG_FLAG_DISCONNECT_SALTS ) && prep_inp_data->bDisconnectSalts &&
         0 < ( i = DisconnectSalts( prep_inp_data, 1 ) )) /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    {
        if (!ip->bNoWarnings)
        {
            WarningMessage( sd->pStrErrStruct, "Salt was disconnected" );
        }
        sd->bTautFlagsDone[INCHI_BAS] |= TG_FLAG_DISCONNECT_SALTS_DONE;
        if (sd->nErrorType < _IS_WARNING)
        {
            sd->nErrorType = _IS_WARNING;
        }
        if ((i = ReconcileAllCmlBondParities( prep_inp_data->at, prep_inp_data->num_inp_atoms, 0 ))) /* djb-rwth: addressing LLVM warning */
        {
            char szErrCode[16];
            sprintf(szErrCode, "%d", i);
            AddErrorMessage( sd->pStrErrStruct, "0D Parities Reconciliation failed:" );
            AddErrorMessage( sd->pStrErrStruct, szErrCode );
        }

#if ( bRELEASE_VERSION == 0 )
        sd->bExtract |= EXTR_SALT_WAS_DISCONNECTED;
#endif
    }
    else
    {
        prep_inp_data->bDisconnectSalts = 0;
    }

    /*  Mark the (disconnected) components in prep_inp_data    */

    prep_inp_data->num_components = MarkDisconnectedComponents( prep_inp_data, 0 );

    if (prep_inp_data->num_components < 0)
    {
        AddErrorMessage( sd->pStrErrStruct, "Out of RAM" );
        sd->nStructReadError = 99;
        sd->nErrorType = _IS_FATAL;
        goto exit_function;
    }

    /* Detect isotopic H on heteroatoms -- necessary condition
       for global isotopic tautomerism */

    if ((i = bNumHeterAtomHasIsotopicH( prep_inp_data->at, prep_inp_data->num_inp_atoms ))) /* djb-rwth: addressing LLVM warning */
    {
        if (i & 1)
        {
            sd->bTautFlagsDone[INCHI_BAS] |= TG_FLAG_FOUND_ISOTOPIC_H_DONE;
        }
        if (i & 2)
        {
            sd->bTautFlagsDone[INCHI_BAS] |= TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE;
        }
    }

    /* 4a. Detect unusual valences                                              */

    if (OrigAtData_bCheckUnusualValences( prep_inp_data, 1, sd->pStrErrStruct, ip->bNoWarnings ))
    {
#if ( bRELEASE_VERSION == 0 )
        sd->bExtract |= EXTR_UNUSUAL_VALENCES;
#else
        ;
#endif
    }


    /*    5. if( orig_inp_data->bDisconnectCoord ) then
              -- copy prep_inp_data --> prep_inp_data+1
              -- disconnect metals in prep_inp_data            */

    if (prep_inp_data->bDisconnectCoord)
    {

        prep_inp_data->num_components = MarkDisconnectedComponents( prep_inp_data, 0 );
        if (prep_inp_data->num_components < 0)
        {
            AddErrorMessage( sd->pStrErrStruct, "Out of RAM" );
            sd->nStructReadError = 99;
            sd->nErrorType = _IS_FATAL;
            goto exit_function;
        }

        /* Save reconnected structure in prep_inp_data+1 if requested */
        if (0 != ( ip->bTautFlags & TG_FLAG_RECONNECT_COORD ))
        {
            if (0 > OrigAtData_Duplicate( prep_inp_data + 1, prep_inp_data ))
            {
                AddErrorMessage( sd->pStrErrStruct, "Out of RAM" );
                sd->nStructReadError = 99;
                sd->nErrorType = _IS_FATAL;
                goto exit_function;
            }
            sd->bTautFlags[INCHI_REC] = sd->bTautFlags[INCHI_BAS];
            sd->bTautFlagsDone[INCHI_REC] = sd->bTautFlagsDone[INCHI_BAS];
            {
                /* Remove "parity undefined in disconnected structure" flag from reconnected structure */
                int k, m; /* djb-rwth: removing redundant variables */
                inp_ATOM *at = ( prep_inp_data + 1 )->at;
                int       num_at = ( prep_inp_data + 1 )->num_inp_atoms;
                for (k = 0; k < num_at; k++)
                {
                    for (m = 0; m < MAX_NUM_STEREO_BONDS && at[k].sb_parity[m]; m++) /* djb-rwth: removing redundant code */
                    {
                        at[k].sb_parity[m] &= SB_PARITY_MASK;
                    }
                }
            }
        }

        /* Make disconnected structure in prep_inp_data */
        i = ( 0 != ( ip->bTautFlags & TG_FLAG_CHECK_VALENCE_COORD ) );

        /*    prep_inp_data->bDisconnectCoord > 1 means add
                prep_inp_data->bDisconnectCoord-1 explicit H atoms    */
        if (0 < ( i = DisconnectMetals( prep_inp_data, i, &bTautFlagsDone ) ))
        {
            if (!ip->bNoWarnings)
            {
                WarningMessage( sd->pStrErrStruct, "Metal was disconnected" );
            }
            sd->bTautFlagsDone[INCHI_BAS] |= TG_FLAG_DISCONNECT_COORD_DONE;
            if (sd->nErrorType < _IS_WARNING)
            {
                sd->nErrorType = _IS_WARNING;
            }

#if ( bRELEASE_VERSION == 0 )
            sd->bExtract |= EXTR_METAL_WAS_DISCONNECTED;
#endif

            /* last parm=1 means find link between unchanged by Metal Disconnection components */
            prep_inp_data->num_components = MarkDisconnectedComponents( prep_inp_data, 1 );

            if (prep_inp_data->num_components < 0)
            {
                AddErrorMessage( sd->pStrErrStruct, "Out of RAM" );
                sd->nStructReadError = 99;
                sd->nErrorType = _IS_FATAL;
                goto exit_function;
            }

            {
                /* Set parities for the disconnected structure */
                int k, m, p;
                inp_ATOM *at = ( prep_inp_data )->at;
                int       num_at = ( prep_inp_data )->num_inp_atoms;
                for (k = 0; k < num_at; k++)
                {
                    for (m = 0; m < MAX_NUM_STEREO_BONDS && ( p = at[k].sb_parity[m] ); m++)
                    {
                        if (p & SB_PARITY_FLAG)
                        {
                            at[k].sb_parity[m] = ( p >> SB_PARITY_SHFT ) & SB_PARITY_MASK;
                        }
                    }
                }
            }

            if ((i = ReconcileAllCmlBondParities( prep_inp_data->at, prep_inp_data->num_inp_atoms, 1 ))) /* djb-rwth: addressing LLVM warning */
            {
                char szErrCode[16];
                sprintf(szErrCode, "%d", i);
                AddErrorMessage( sd->pStrErrStruct, "0D Parities Reconciliation failed:" );
                AddErrorMessage( sd->pStrErrStruct, szErrCode );
            }

#if ( REMOVE_ION_PAIRS_DISC_STRU == 1 )
            if (0 < remove_ion_pairs( prep_inp_data->num_inp_atoms, prep_inp_data->at ))
            {
                if (!ip->bNoWarnings)
                {
                    WarningMessage( sd->pStrErrStruct, "Charges were rearranged" );
                }
                if (sd->nErrorType < _IS_WARNING)
                {
                    sd->nErrorType = _IS_WARNING;
                }
                sd->bTautFlagsDone[INCHI_REC] |= TG_FLAG_FIX_ODD_THINGS_DONE;
                sd->bTautFlagsDone[INCHI_BAS] |= TG_FLAG_FIX_ODD_THINGS_DONE;
            }
#endif

            /*
              if prep_inp_data->nOldCompNumber[i] = iINChI+1 > 0 then
              component #(i+1) in prep_inp_data is identical to component #(iINChI+1) in prep_inp_data+1
            */
        }
        else if (i < 0)
        {
            AddErrorMessage( sd->pStrErrStruct, "Cannot disconnect metal error" );
            sd->nStructReadError = i;
            sd->nErrorType = _IS_ERROR;
            goto exit_function;
        }
    }
    else
    {
        /* Remove "disconnected structure parities" from the structure */
        int k, m; /* djb-rwth: removing redundant variables */
        inp_ATOM *at = ( prep_inp_data )->at;
        int       num_at = ( prep_inp_data )->num_inp_atoms;
        for (k = 0; k < num_at; k++)
        {
            for (m = 0; m < MAX_NUM_STEREO_BONDS && at[k].sb_parity[m]; m++) /* djb-rwth: removing redundant code */
            {
                at[k].sb_parity[m] &= SB_PARITY_MASK;
            }
        }
    }

exit_function:

if (sd->nErrorType < _IS_ERROR && prep_inp_data)
    {
        if (0 < post_fix_odd_things( prep_inp_data->num_inp_atoms, prep_inp_data->at ))
        {
            if (!ip->bNoWarnings)
            {
                WarningMessage( sd->pStrErrStruct, "Charges were rearranged" );
            }
            if (sd->nErrorType < _IS_WARNING)
            {
                sd->nErrorType = _IS_WARNING;
            }
            sd->bTautFlagsDone[INCHI_BAS] |= TG_FLAG_FIX_ODD_THINGS_DONE;
        }
        if (( sd->bTautFlagsDone[INCHI_BAS] & TG_FLAG_DISCONNECT_COORD_DONE ) &&
            ( prep_inp_data + 1 )->at && ( prep_inp_data + 1 )->num_inp_atoms > 0)
        {
            if (0 < post_fix_odd_things( ( prep_inp_data + 1 )->num_inp_atoms, ( prep_inp_data + 1 )->at ))
            {
                if (!ip->bNoWarnings)
                {
                    WarningMessage( sd->pStrErrStruct, "Charges were rearranged" );
                }
                if (sd->nErrorType < _IS_WARNING)
                {
                    sd->nErrorType = _IS_WARNING;
                }
                sd->bTautFlagsDone[INCHI_REC] |= TG_FLAG_FIX_ODD_THINGS_DONE;
                sd->bTautFlagsDone[INCHI_BAS] |= TG_FLAG_FIX_ODD_THINGS_DONE;
            }
        }
    }

    sd->bTautFlags[INCHI_BAS] |= bTautFlags;  /* TG_FLAG_CHECK_VALENCE_COORD_DONE, TG_FLAG_MOVE_CHARGE_COORD_DONE */
    sd->bTautFlagsDone[INCHI_BAS] |= bTautFlagsDone;  /* TG_FLAG_CHECK_VALENCE_COORD_DONE, TG_FLAG_MOVE_CHARGE_COORD_DONE */

    return sd->nErrorType;
}


#ifndef TARGET_API_LIB


/****************************************************************************/
int CreateCompositeNormAtom( COMP_ATOM_DATA  *composite_norm_data,
                            INP_ATOM_DATA2  *all_inp_norm_data,
                            int             num_components )
{
    int i, j, jj, k, n, m, tot_num_at, tot_num_H, cur_num_at, cur_num_H; /* djb-rwth: removing redundant variables */
    int num_comp[TAUT_NUM + 1], num_taut[TAUT_NUM + 1], num_del[TAUT_NUM + 1], num_at[TAUT_NUM + 1], num_inp_at[TAUT_NUM + 1];
    int ret = 0, indicator = 1;
    inp_ATOM *at, *at_from;
    memset( num_comp, 0, sizeof( num_comp ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    memset( num_taut, 0, sizeof( num_taut ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    memset( num_del, 0, sizeof( num_taut ) ); /* djb-rwth: memset_s C11/Annex K variant? */

    /* count taut and non-taut components */
    for (j = 0; j < TAUT_NUM; j++)
    {
        num_comp[j] = num_taut[j] = 0;
        for (i = 0; i < num_components; i++)
        {
            if (all_inp_norm_data[i][j].bExists)
            {
                num_del[j] += ( 0 != all_inp_norm_data[i][j].bDeleted );
                num_comp[j] ++;
                num_taut[j] += ( 0 != all_inp_norm_data[i][j].bTautomeric );
            }
        }
    }

    /* count intermediate taut structure components */
    if (num_comp[TAUT_YES] > num_del[TAUT_YES] && num_taut[TAUT_YES])
    {
        /*
        num_comp[TAUT_INI] = num_comp[TAUT_YES] - num_del[TAUT_YES];
        */

        for (i = 0, j = TAUT_YES; i < num_components; i++)
        {
            if (all_inp_norm_data[i][j].bExists &&
                ( all_inp_norm_data[i][j].bDeleted ||
                    (all_inp_norm_data[i][j].bTautomeric &&
                    all_inp_norm_data[i][j].at_fixed_bonds &&
                    all_inp_norm_data[i][j].bTautPreprocessed) )) /* djb-rwth: addressing LLVM warning */
            {
                num_comp[TAUT_INI] ++;
            }
        }
    }

    /* count atoms and allocate composite atom data */
    for (jj = 0; jj <= TAUT_INI; jj++)
    {
        num_at[jj] = num_inp_at[jj] = 0;
        j = inchi_min( jj, TAUT_YES );
        if (num_comp[jj])
        {
            for (i = 0; i < num_components; i++)
            {
                if (all_inp_norm_data[i][j].bDeleted)
                {
                    continue;
                }
                /* find k = the normaized structure index */
                if (jj == TAUT_INI)
                {
                    if (all_inp_norm_data[i][j].bExists &&
                         all_inp_norm_data[i][j].at_fixed_bonds)
                    {
                        k = j;
                    }
                    else
                        if (all_inp_norm_data[i][ALT_TAUT( j )].bExists && !all_inp_norm_data[i][ALT_TAUT( j )].bDeleted &&
                             !all_inp_norm_data[i][j].bDeleted)
                        {
                            k = ALT_TAUT( j );
                        }
                        else
                        {
                            if (all_inp_norm_data[i][j].bExists)
                            {
                                k = j;
                            }
                            else
                            {
                                continue;
                            }
                        }
                }
                else
                {
                    if (all_inp_norm_data[i][j].bExists)
                    {
                        k = j;
                    }
                    else
                    {
                        if (all_inp_norm_data[i][ALT_TAUT( j )].bExists && !all_inp_norm_data[i][ALT_TAUT( j )].bDeleted)
                        {
                            k = ALT_TAUT( j );
                        }
                        else
                        {
                            continue;
                        }
                    }
                }
                num_inp_at[jj] += all_inp_norm_data[i][k].num_at; /* all atoms including terminal H */
                num_at[jj] += all_inp_norm_data[i][k].num_at - all_inp_norm_data[i][k].num_removed_H;
            }
            if (num_inp_at[jj])
            {
                if (!CreateCompAtomData( composite_norm_data + jj, num_inp_at[jj], num_components, jj == TAUT_INI ))
                {
                    goto exit_error;
                }
                composite_norm_data[jj].num_removed_H = num_inp_at[jj] - num_at[jj];
            }
        }
    }

    /* fill out composite atom */
    for (jj = 0; jj <= TAUT_INI; jj++, indicator <<= 1)
    {
        j = inchi_min( jj, TAUT_YES );
        if (num_comp[jj])
        {
            tot_num_at = 0;
            tot_num_H = 0;
            for (i = 0; i < num_components; i++)
            {
                if (all_inp_norm_data[i][j].bDeleted)
                {
                    composite_norm_data[jj].nNumRemovedProtons += all_inp_norm_data[i][j].nNumRemovedProtons;
                    for (n = 0; n < NUM_H_ISOTOPES; n++)
                    {
                        composite_norm_data[jj].nNumRemovedProtonsIsotopic[n] += all_inp_norm_data[i][j].nNumRemovedProtonsIsotopic[n];
                    }
                    continue;
                }
                /* djb-rwth: removing redundant code */
                /* find k = the normaized structure index */
                if (jj == TAUT_INI)
                {
                    if (all_inp_norm_data[i][j].bExists && all_inp_norm_data[i][j].at_fixed_bonds)
                    {
                        k = j;
                    }
                    else
                    {
                        if (all_inp_norm_data[i][ALT_TAUT( j )].bExists)
                        {
                            k = ALT_TAUT( j );
                        }
                        else
                        {
                            if (all_inp_norm_data[i][j].bExists && !all_inp_norm_data[i][ALT_TAUT( j )].bDeleted)
                            {
                                k = j;
                            }
                            else
                            {
                                continue;
                            }
                        }
                    }
                }
                else
                {
                    if (all_inp_norm_data[i][j].bExists)
                    {
                        k = j;
                    }
                    else
                    {
                        if (all_inp_norm_data[i][ALT_TAUT( j )].bExists && !all_inp_norm_data[i][ALT_TAUT( j )].bDeleted)
                        {
                            k = ALT_TAUT( j );
                        }
                        else
                        {
                            continue;
                        }
                    }
                }
                /* copy main atoms */
                cur_num_H = all_inp_norm_data[i][k].num_removed_H;       /* number of terminal H atoms */
                cur_num_at = all_inp_norm_data[i][k].num_at - cur_num_H;  /* number of all but explicit terminal H atoms */

                if (( tot_num_at + cur_num_at ) > num_at[jj] ||
                    ( num_at[jj] + tot_num_H + cur_num_H ) > num_inp_at[jj])
                {
                    goto exit_error; /* miscount */
                }
                at = composite_norm_data[jj].at + tot_num_at; /* points to the 1st destination atom */
                at_from = ( jj == TAUT_INI && k == TAUT_YES && all_inp_norm_data[i][k].at_fixed_bonds ) ?
                    all_inp_norm_data[i][k].at_fixed_bonds : all_inp_norm_data[i][k].at;
                memcpy(at, at_from, sizeof(composite_norm_data[0].at[0])* cur_num_at); /* copy atoms except terminal H */
                /* shift neighbors of main atoms */
                for (n = 0; n < cur_num_at; n++, at++)
                {
                    for (m = 0; m < at->valence; m++)
                    {
                        at->neighbor[m] += tot_num_at;
                    }
                }
                /* copy explicit H */
                if (cur_num_H)
                {
                    at = composite_norm_data[jj].at + num_at[jj] + tot_num_H; /* points to the 1st destination atom */
                    memcpy(at, at_from + cur_num_at,
                        sizeof(composite_norm_data[0].at[0])* cur_num_H);
                    /* shift neighbors of explicit H atoms */
                    for (n = 0; n < cur_num_H; n++, at++)
                    {
                        for (m = 0; m < at->valence; m++)
                        {
                            at->neighbor[m] += tot_num_at;
                        }
                    }
                }
                /* composite counts */
                composite_norm_data[jj].bHasIsotopicLayer |= all_inp_norm_data[i][k].bHasIsotopicLayer;
                composite_norm_data[jj].num_isotopic += all_inp_norm_data[i][k].num_isotopic;
                composite_norm_data[jj].num_bonds += all_inp_norm_data[i][k].num_bonds;
                composite_norm_data[jj].bTautomeric += ( j == jj ) && all_inp_norm_data[i][k].bTautomeric;
                composite_norm_data[jj].nNumRemovedProtons += all_inp_norm_data[i][k].nNumRemovedProtons;
                for (n = 0; n < NUM_H_ISOTOPES; n++)
                {
                    composite_norm_data[jj].nNumRemovedProtonsIsotopic[n] += all_inp_norm_data[i][k].nNumRemovedProtonsIsotopic[n];
                    composite_norm_data[jj].num_iso_H[n] += all_inp_norm_data[i][k].num_iso_H[n];
                }
                /*
                composite_norm_data[j].num_at            += cur_num_at + cur_num_H;
                composite_norm_data[j].num_removed_H     += cur_num_H;
                */
                /* total count */
                tot_num_at += cur_num_at;
                tot_num_H += cur_num_H;
                /* offset for the next component */
                if (composite_norm_data[jj].nOffsetAtAndH)
                {
                    composite_norm_data[jj].nOffsetAtAndH[2 * i] = tot_num_at;
                    composite_norm_data[jj].nOffsetAtAndH[2 * i + 1] = num_at[jj] + tot_num_H;
                }
            }
            if (tot_num_at != num_at[jj] ||
                 num_at[jj] + tot_num_H != num_inp_at[jj])
            {
                goto exit_error; /* miscount */
            }
            composite_norm_data[jj].bExists = ( tot_num_at > 0 );
            ret |= indicator;
        }
    }
    return ret;

exit_error:

    return ret;
}
#endif


/****************************************************************************/
void OrigAtData_DebugTrace( ORIG_ATOM_DATA* d )
{
    int i, k;

    ITRACE_( "\n\n*********************************************************************\n* ORIG_ATOM_DATA @ 0x%p", d );
    ITRACE_( "\n*  num_inp_atoms = %-d\n*  num_inp_bonds = %-d\n*  num_dimensions = %-d\n*  num_components = %-d",
            d->num_inp_atoms, d->num_inp_bonds, d->num_dimensions, d->num_components );
    ITRACE_( "\n*  ATOMS" );
    for (i = 0; i < d->num_inp_atoms; i++)
    {
        ITRACE_( "\n*    #%-5d %s%-d ( charge %-d, rad %-d nH %-d val %-d) [%-f %-f %-f]",
            i, d->at[i].elname, d->at[i].orig_at_number, d->at[i].charge, d->at[i].radical, d->at[i].num_H, d->at[i].valence,
            d->at[i].x, d->at[i].y, d->at[i].z );
        if (d->at[i].valence > 0)
        {
            ITRACE_( "\n            bonds to     " );
            for (k = 0; k < d->at[i].valence; k++)
            {
                int nbr = d->at[i].neighbor[k]; /* djb-rwth: ignoring LLVM warning: possible presence of global variables */
                ITRACE_( "%s%-3d ", d->at[nbr].elname, nbr + 1 );
            }
        }
        if (d->at[i].valence > 0)
        {
            ITRACE_( "\n            bond types   " );
            for (k = 0; k < d->at[i].valence; k++)
                ITRACE_( "%-3d ", d->at[i].bond_type[k] );
        }
    }
    /*OAD_Polymer_DebugTrace( d->polymer );*/
    ITRACE_( "\n* V3000 INFO @ 0x%-p", d->v3000 );
    ITRACE_( "\n*\n" );
    if (d->v3000)
    {
        ITRACE_( "\n*  n_star_atoms = %-d\n*  n_haptic_bonds = %-d\n*  n_collections = %-d",
            d->v3000->n_star_atoms, d->v3000->n_haptic_bonds, d->v3000->n_collections );
    }
    ITRACE_( "\n*\n* End ORIG_ATOM_DATA\n*********************************************************************\n" );

    return;
}



/*
    Polymer related procedures
*/




/****************************************************************************
 Create a new OAD_PolymerUnit
****************************************************************************/
OAD_PolymerUnit* OAD_PolymerUnit_New( int       maxatoms,
                                      int       maxbonds,
                                      int       id,
                                      int       label,
                                      int       type,
                                      int       subtype,
                                      int       conn,
                                      char      *smt,
                                      int       na,
                                      INT_ARRAY *alist,
                                      int       nb,
                                      INT_ARRAY *blist,
                                      int       nbkbonds,
                                      int       **bkbonds )
{
    int k, err = 0;
    OAD_PolymerUnit *u2 = NULL;

    u2 = (OAD_PolymerUnit*) inchi_calloc( 1, sizeof( OAD_PolymerUnit ) );
    if (NULL == u2)
    {
        err = 1;
        goto exit_function;
    }

    u2->id = id;
    u2->label = label;
    u2->type = type;
    u2->subtype = subtype;
    u2->conn = conn;
    u2->na = na;
    u2->nb = nb;
    u2->cyclizable = CLOSING_SRU_NOT_APPLICABLE;
    u2->cyclized = 0;
    for (k = 0; k < 4; k++)
    {
        u2->xbr1[k] = 0.0;
        u2->xbr2[k] = 0.0;
    }
    strcpy(u2->smt, smt);
    u2->cap1 = -1;
    u2->end_atom1 = -1;
    u2->cap2 = -1;
    u2->end_atom2 = -1;
    u2->maxbkbonds = maxbonds;
    u2->nbkbonds = nbkbonds;
    u2->cap1_is_undef = 0;
    u2->cap2_is_undef = 0;

    u2->alist = NULL;
    if (na > 0 || maxatoms > 0)
    {
        u2->alist	= (int *) inchi_calloc( na > 0 ? na : maxatoms, sizeof( int ) );
        if (!u2->alist )
        {
            err = 2;
            goto exit_function;
        }
        for (k = 0; k < na; k++)
        {
            u2->alist[k]	= alist->item[k];
        }
    }
    u2->blist = NULL;
    if (nb > 0 || maxbonds > 0)
    {
        u2->blist = (int *) inchi_calloc( nb > 0 ? 2 * nb : 2 * maxbonds, sizeof( int ) );
        if (!u2->blist)
        {
            err = 3;
            goto exit_function;
        }
        if (blist)
        {
            for (k = 0; k < 2 * nb; k++)
            {
                u2->blist[k]	= blist->item[k];
            }
        }
        
    }
    u2->bkbonds = NULL;
    
exit_function:

    if (err)
    {
        OAD_PolymerUnit_Free( u2 );
        return NULL;
    }

    return u2;
}


/****************************************************************************
 Create a copy of OAD_PolymerUnit
****************************************************************************/
OAD_PolymerUnit* OAD_PolymerUnit_CreateCopy( OAD_PolymerUnit *u )
{
    int k, err = 0;
    OAD_PolymerUnit *u2 = NULL;

    u2 = (OAD_PolymerUnit*) inchi_calloc( 1, sizeof( OAD_PolymerUnit ) );
    if (NULL == u2)
    {
        err = 1;
        goto exit_function;
    }
    u2->id = u->id;
    u2->type = u->type;
    u2->subtype = u->subtype;
    u2->conn = u->conn;
    u2->label = u->label;
    u2->na = u->na;
    u2->nb = u->nb;
    u2->cyclizable = u->cyclizable;
    u2->cyclized = u->cyclized;
    u2->cap1_is_undef = u->cap1_is_undef;
    u2->cap2_is_undef = u->cap2_is_undef;

    for (k = 0; k < 4; k++)
    {
        u2->xbr1[k] = u->xbr1[k];
        u2->xbr2[k] = u->xbr2[k];
    }

    strcpy(u2->smt, u->smt);

    u2->cap1 = u->cap1;
    u2->end_atom1 = u->end_atom1;
    u2->cap2 = u->cap2;
    u2->end_atom2 = u->end_atom2;
    u2->nbkbonds = u->nbkbonds;
    u2->maxbkbonds = inchi_max( u->maxbkbonds, u->nbkbonds );

    u2->alist	= (int *) inchi_calloc( u2->na, sizeof( int ) );
    if (!u2->alist)
    {
        err = 2;
        goto exit_function;
    }
    for (k = 0; k < u2->na; k++)
    {
        u2->alist[k]	= u->alist[k];
    }

    u2->blist	= (int *) inchi_calloc( 2 * (long long)u2->nb, sizeof( int ) ); /* djb-rwth: cast operator added */
    if (!u2->blist)
    {
        err = 2;
        goto exit_function;
    }
    for (k = 0; k < 2 * u2->nb; k++)
    {
        u2->blist[k]	= u->blist[k];
    }

    err = imat_new( u2->maxbkbonds, 2, &( u2->bkbonds ) );
    if (!err)
    {
        for (k = 0; k < u2->nbkbonds; k++)
        {
            u2->bkbonds[k][0] = u->bkbonds[k][0];
            u2->bkbonds[k][1] = u->bkbonds[k][1];
        }
    }

exit_function:
    if (err)
    {
        OAD_PolymerUnit_Free( u2 );
        return NULL;
    }

    return u2;
}


/****************************************************************************/
void OAD_PolymerUnit_Free( OAD_PolymerUnit *unit )
{

    ITRACE_( "\n************** About to free OAD_PolymerUnit @ %-p\n", unit );
    OAD_PolymerUnit_DebugTrace( unit );

    if (unit)
    {
        if (unit->alist)
        {
            inchi_free( unit->alist );
            unit->alist = NULL;
        }
        if (unit->blist)
        {
            inchi_free( unit->blist );
            unit->blist = NULL;
        }
        if (unit->bkbonds)
        {
            imat_free( unit->maxbkbonds, unit->bkbonds );
            unit->bkbonds = NULL;
        }
    }

    inchi_free( unit );

    return;
}


/****************************************************************************
 Compare two polymer units, modified lexicographic order
    Modification: unit with smaller alist always go first
****************************************************************************/
int  OAD_PolymerUnit_CompareAtomListsMod( OAD_PolymerUnit *u1,
                                          OAD_PolymerUnit *u2 )
{
    int i;
    int n1 = u1->na;
    int n2 = u2->na;
    int n = n1;
    if (n1 < n2)    return -1;
    if (n1 > n2)    return 1;
    /* n1 == n2 == n */
    for (i = 0; i < n; i++)
    {
        if (u1->alist[i] < u2->alist[i])    return -1;
        if (u1->alist[i] > u2->alist[i])    return    1;
    }

    return 0;
}


/****************************************************************************
 Compare two polymer units, lexicographic order
****************************************************************************/
int  OAD_PolymerUnit_CompareAtomLists( OAD_PolymerUnit *u1,
                                       OAD_PolymerUnit *u2 )
{
    int i;
    int n1 = u1->na;
    int n2 = u2->na;
    int n = inchi_min( n1, n2 );

    for (i = 0; i < n; i++)
    {
        if (u1->alist[i] < u2->alist[i])
        {
            return -1;
        }
        if (u1->alist[i] > u2->alist[i])
        {
            return 1;
        }
    }

    if (n1 < n2)
    {
        return -1;
    }

    if (n1 > n2)
    {
        return    1;
    }

    return 0;
}


/****************************************************************************
 Sort SRU bond lists atoms and bonds themselves
****************************************************************************/
int  OAD_PolymerUnit_OrderBondAtomsAndBondsThemselves( OAD_PolymerUnit  *u,
                                                       int n_star_atoms,
                                                       int *star_atoms )
{
    int k;

    /* Sort bond atoms */
    for (k = 0; k < u->nb; k++)
    {
        /* Place not-in-unit bond end to first place */
        int a1 = u->blist[2 * k];
        int a2 = u->blist[2 * k + 1];
        int a1_is_not_in_alist = 0;
        int a1_is_star_atom = 0;
        int a2_is_not_in_alist = 0;
        int a2_is_star_atom = 0;

        if (!is_in_the_ilist( u->alist, a1, u->na ))
        {
            a1_is_not_in_alist = 1;
        }
        if (is_in_the_ilist( star_atoms, a1, n_star_atoms ))
        {
            a1_is_star_atom = 1;
        }

        if (!is_in_the_ilist( u->alist, a2, u->na ))
        {
            a2_is_not_in_alist = 1;
        }
        if (is_in_the_ilist( star_atoms, a2, n_star_atoms ))
        {
            a2_is_star_atom = 1;
        }

        if (( a1_is_not_in_alist || a1_is_star_atom ) &&
            ( a2_is_not_in_alist || a2_is_star_atom ))
        {
            /* Both the ends are out of unit: the crossing bond is invalid */
            return 1;
        }
        /* If a2 is star atom or non-star external to the current unit, swap(a2,a1) */
        if (a2_is_star_atom || a2_is_not_in_alist)
        {
            u->blist[2 * k] = a2;
            u->blist[2 * k + 1] = a1;
        }
    }

    /* Sort bond themselves
        for now, consider only the simplest cases of 2 bonds
    */
    if (u->nb == 2)            /* two bonds in SBL */
    {
        int b1a1 = u->blist[0];
        int b1a2 = u->blist[1];
        int b2a1 = u->blist[2];
        int b2a2 = u->blist[3];
        if (b1a1 > b2a1)
        {
            /* swap */
            u->blist[0] = b2a1; u->blist[1] = b2a2;
            u->blist[2] = b1a1; u->blist[3] = b1a2;
        }
    }

    /* for single or no bonds, do nothing
    else
        ;
    */

    return 0;
}


/****************************************************************************
 Parse pseudoelement and polymer data
 (unit, types, subtypes, connections, etc.)
****************************************************************************/
int OAD_ValidatePolymerAndPseudoElementData( ORIG_ATOM_DATA *orig_at_data,
                                             int treat_polymers,
                                             int bNPZz,
                                             char *pStrErr,
                                             int bNoWarnings)
{
    int i, k, kk, type, subtype, representation, err = 0;
    int nat = orig_at_data->num_inp_atoms;
    int nsgroups = 0;
    OAD_PolymerUnit* u = NULL;
    OAD_Polymer *pd = orig_at_data->polymer;


    /* Assign polymer type and subunits type and check polymer data for consistency */
    /* djb-rwth: addressing coverity ID #499497 -- TREAT_ERR properly used in all cases */
    
    orig_at_data->valid_polymer = 0;
    if (treat_polymers && pd)
    {
        orig_at_data->valid_polymer = 1;
    }
    if (orig_at_data->valid_polymer)
    {
        nsgroups = pd->n;
    }
    if (nsgroups == 1)
    {
        /* Check if copolymer */
        type = pd->units[0]->type;
        if (type == POLYMER_STY_COP)
        {
            TREAT_ERR( err, 9001, "Copolymer must contain more than one unit" );
            goto exit_function;
        }
        /* Check if copolymer subtype */
        subtype = pd->units[0]->subtype;
        if (subtype == POLYMER_SST_RAN || subtype == POLYMER_SST_ALT || subtype == POLYMER_SST_BLK)
        {
            TREAT_ERR( err, 9002, "Single polymer unit may not be RAN/ALT/BLO" );
            goto exit_function;
        }
    }

    /* For each CRU */
    for (i = 0; i < nsgroups; i++)
    {
        /* Check if unit data makes sense */
        u = pd->units[i];
        if (u->nb != 0 && u->nb != 2)
        {
            TREAT_ERR( err, 9003, "Number of crossing bonds in polymer unit is not 0 or 2" );
            goto exit_function;
        }
        if (u->na < 1)
        {
            TREAT_ERR( err, 9004, "Empty polymer unit" );
            goto exit_function;
        }
        if (u->na > nat)
        {
            TREAT_ERR( err, 9005, "Too large polymer unit" );
            goto exit_function;
        }
        for (k = 0; k < u->na; k++)
        {
            int atom = u->alist[k];
            if (atom < 1 || atom > nat)
            {
                TREAT_ERR( err, 9006, "Invalid atom number in polymer unit" );
                goto exit_function;
            }
            /* was not accounting for COP ...
            if (is_in_the_ilist( pd->pzz, atom, pd->n_pzz ))
            {
                TREAT_ERR( err, 9007, "Star atom inside polymer unit" );
                goto exit_function;
            }
            */
        }

        OAD_PolymerUnit_SetEndsAndCaps( u, orig_at_data, &err, pStrErr );
            /*	Reveal and store CRU caps and ends('stars and partners')
                Also set `unit->cap1_is_undef`, `unit->cap2_is_undef`, `unit->cyclizable` 
            */
        if (err)
        {
            goto exit_function;
        }

        
        /* Set possibly missing unit parameters */
        u->nbkbonds = 0;
        u->cyclizable = CLOSING_SRU_NOT_APPLICABLE;
        u->cyclized = 0;
    }
    

    OAD_ValidateAndSortOutPseudoElementAtoms( orig_at_data, treat_polymers, bNPZz, &err, pStrErr );
    /* Here we:
                Make more polymer and pseudoatom data checks
                Convert both "*" and "Zz" temporarily to "Zy" (polymer-unrelated interal pseudoatoms)
                If applicable, check each CRU and back-convert "Zy" to "Zz" (polymer-related 
                pseudoelement atoms) if they are for valid bi-undef-end CRU
    */
    
    if (err)
    {
        /* already treated TREAT_ERR( err, 9040, "Improper pseudoelement atoms" ); */
        goto exit_function;
    }


    /* Make more polymer and pseudoatom data checks */

    /* Check if non-polymer-related Zz/star atoms enabled */
    if (orig_at_data->n_zy > 0 && bNPZz==0 )
    {
        TREAT_ERR( err, 9, "Non-polymer-related Zz/star atoms are not allowed" );
        goto exit_function;
    }

    if (!pd || !orig_at_data->valid_polymer )
    {
        goto exit_function;
    }

    if (pd->n_pzz > 0)
    {
        /* Allocate memory for polymer-related pseudoatoms */
        if ( pd->treat == POLYMERS_NO)
        {
            TREAT_ERR( err, 9, "Pseudoelement endgroups are not allowed" );
            goto exit_function;
        }
        if (pd->pzz)
        {
            inchi_free(pd->pzz);
            pd->pzz = NULL;
        }
        pd->pzz = (int *) inchi_calloc( pd->n_pzz, sizeof( int ) );
        if (!pd->pzz)
        {
            TREAT_ERR( err, 9010, "Not enough memory" );
            goto exit_function;
        }
        kk = 0;
        for (k = 0; k < nat; k++)
        {
            if (!strcmp( orig_at_data->at[k].elname, "Zz" ))
            {
                pd->pzz[kk++] = k + 1; /* djb-rwth: buffer overrun avoided implicitly */
            }
        }
    }

    /* Check copolymers and ensure that COP includes > 1 SRU */
    for (i = 0; i < pd->n; i++)
    {
        u = pd->units[i];

        if ( u->type == POLYMER_STY_COP ||
             u->type == POLYMER_STY_SRU     /* what drawn as 'SRU' [xyz]n may be actually copolymer [xyz]co */
           )
        {
            int j, in_units = 0;

            if (u->nb > 0)
            {
                /* crossing bonds present, either valid SRU or invalid copolymer */
                if (u->type == POLYMER_STY_COP)
                {
                    TREAT_ERR( err, 9026, "Polymer COP unit contains bracket-crossing bonds, not supported" );
                    goto exit_function;
                }
                else
                {
                    continue;
                }
            }
            /* now we have no crossing bonds units */
            for (j = 0; j < pd->n; j++)
            {
                if (pd->units[j]->type == POLYMER_STY_COP)
                {
                    continue;
                }
                if (is_ilist_inside( pd->units[j]->alist, pd->units[j]->na, pd->units[i]->alist, pd->units[i]->na ))
                {
                    in_units++;
                    if (in_units == 2)
                    {
                        break;
                    }
                }
            }
            if (in_units > 1)
            {
                if (u->type != POLYMER_STY_COP)
                {
                    u->type = POLYMER_STY_COP;
                    if (!bNoWarnings)
                    {
                        WarningMessage( pStrErr, "Convert multiple-subunits unit to copolymer" );
                    }
                }
            }
            else    /* in_units <= 1)*/
            {
                if (u->type == POLYMER_STY_COP)
                {
                    TREAT_ERR( err, 9027, "Polymer COP unit contains a single SRU instead of multiple" );
                    goto exit_function;
                }
            }
        }
    }

    representation = OAD_Polymer_GetRepresentation( pd );

    /* Make more polymer data checks and perform some corrections*/
    if (representation == POLYMER_REPRESENTATION_SOURCE_BASED)
    {
        for (i = 0; i < nsgroups; i++)
        {
            /* Replace source-based 'SRU' with 'MON' */
            if (pd->units[i]->type == POLYMER_STY_SRU)
            {
                pd->units[i]->type = POLYMER_STY_MON;
                if (!bNoWarnings)
                {
                    WarningMessage( pStrErr, "Converted src-based polymer unit type to MON" );
                }
            }
            if (pd->units[i]->type == POLYMER_STY_COP)
            {
                /* Set missing copolymer subtype to RAN */
                if (pd->units[i]->subtype == POLYMER_SST_NON)
                {
                    pd->units[i]->subtype = POLYMER_SST_RAN;
                    if (!bNoWarnings)
                    {
                        WarningMessage( pStrErr, "Set missing copolymer subtype to RAN" );
                    }
                }
            }
            /* Suppress connectivity ("HH", "HT", "EU") */
            if (pd->units[i]->conn != POLYMER_CONN_NON)
            {
                pd->units[i]->conn = POLYMER_CONN_NON;
                if (!bNoWarnings)
                {
                    WarningMessage( pStrErr, "Ignore connection pattern for src-based polymer unit" );
                }
            }
        }
    }

#ifdef ALLOW_MIXED_SRU_AND_MON
    else if (representation == POLYMER_REPRESENTATION_STRUCTURE_BASED ||
              representation == POLYMER_REPRESENTATION_MIXED)
#else
    else if (representation == POLYMER_REPRESENTATION_STRUCTURE_BASED)
#endif
    {
        for (i = 0; i < nsgroups; i++)
        {
            int a1, a2, a1_is_not_in_alist, a1_is_star_atom, a2_is_not_in_alist, a2_is_star_atom;

            u = pd->units[i];

            /*    SRU that is copolymer unit embedding other SRU's */
            if (u->nb == 0)
            {
                if (u->type == POLYMER_STY_COP)
                {
                    ;
                }
                else if (u->type == POLYMER_STY_SRU)
                {
                    u->type = POLYMER_STY_COP;
                    if (!bNoWarnings)
                    {
                        WarningMessage( pStrErr, "Set copolymer embedding unit mark to COP" );
                    }
                }
            }
            if (u->type == POLYMER_STY_COP)
            {
                u->cyclizable = CLOSING_SRU_NOT_APPLICABLE;
                /* Set possibly missing copolymer subtype to RAN */
                if (u->subtype == POLYMER_SST_NON)
                {
                    u->subtype = POLYMER_SST_RAN;
                    if (!bNoWarnings)
                    {
                        WarningMessage( pStrErr, "Set missing copolymer subtype to RAN" );
                    }
                }
                continue;
            }

#ifdef ALLOW_MIXED_SRU_AND_MON
            if (u->type == POLYMER_STY_MON)
            {
                continue;
            }
#endif
            /*    SRU with endgroups or stars. Check it. */
            for (k = 0; k < u->nb; k++)
            {
                /* Check that there are no H end groups */
                a1 = u->blist[2 * k]; a2 = u->blist[2 * k + 1];
                if (!strcmp( orig_at_data->at[a1 - 1].elname, "H" ) ||
                     !strcmp( orig_at_data->at[a1 - 1].elname, "D" ) ||
                     !strcmp( orig_at_data->at[a1 - 1].elname, "T" ))
                {
                    TREAT_ERR( err, 9030, "H as polymer end group is not supported" );
                    goto exit_function;
                }
                if ( !strcmp( orig_at_data->at[a2 - 1].elname, "H" ) ||
                     !strcmp( orig_at_data->at[a2 - 1].elname, "D" ) ||
                     !strcmp( orig_at_data->at[a2 - 1].elname, "T" ))
                {
                    TREAT_ERR( err, 9031, "H as polymer end group is not supported" );
                    goto exit_function;
                }
                /* Ensure that caps of polymer unit lie outside it */
                a1_is_not_in_alist = a1_is_star_atom = 0;
                a2_is_not_in_alist = a2_is_star_atom = 0;
                if (!is_in_the_ilist( u->alist, a1, u->na ))
                {
                    a1_is_not_in_alist = 1;
                }
                if (is_in_the_ilist( pd->pzz, a1, pd->n_pzz ))
                {
                    a1_is_star_atom = 1;
                }
                if (!is_in_the_ilist( u->alist, a2, u->na ))
                {
                    a2_is_not_in_alist = 1;
                }
                if (is_in_the_ilist( pd->pzz, a2, pd->n_pzz ))
                {
                    a2_is_star_atom = 1;
                }
                if (( a1_is_not_in_alist || a1_is_star_atom ) &&
                    ( a2_is_not_in_alist || a2_is_star_atom ))
                {
                    TREAT_ERR( err, 9032, "Caps of polymer unit lie inside it" );
                    goto exit_function;
                }
            }

            if (u->type == POLYMER_STY_SRU || u->type == POLYMER_STY_MOD ||
                 u->type == POLYMER_STY_CRO || u->type == POLYMER_STY_MER)
            {
                /* If SRU connection is missing, set to default ('either') */
                if (u->conn == POLYMER_CONN_NON)
                {
                    if (!bNoWarnings)
                    {
                        WarningMessage( pStrErr, "Set missing copolymer unit connection to EU" );
                    }
                    u->conn = POLYMER_CONN_EU;
                }

                if (u->cap1 && u->cap2)
                {
                    /* Set SRU closure type */
                    if (u->na == 1)
                    {
#ifdef ALLOW_CLOSING_SRU_VIA_DIRADICAL
                        u->cyclizable = CLOSING_SRU_DIRADICAL;
#else
                        u->cyclizable = CLOSING_SRU_NOT_APPLICABLE;
#ifdef  CLOSING_STARRED_SRU_IS_A_MUST
                        TREAT_ERR( err, 9029, "Could not perform SRU closure" );
                        goto exit_function;
#endif
#endif
                    }
                    else if (u->na == 2)
                    {

#ifdef ALLOW_CLOSING_SRU_VIA_HIGHER_ORDER_BOND
                        u->cyclizable = CLOSING_SRU_HIGHER_ORDER_BOND;
#else
                        u->cyclizable = CLOSING_SRU_NOT_APPLICABLE;
#ifdef  CLOSING_STARRED_SRU_IS_A_MUST
                        TREAT_ERR( err, 9029, "Could not perform SRU closure" );
                        goto exit_function;
#endif
#endif
                    }
                    else
                    {
                        u->cyclizable = CLOSING_SRU_RING;
                    }
                }

                if (u->conn != POLYMER_CONN_HT)
                {
                    /* frame shift/SRU cyclization is for head-to-tail connections only */
                    u->cyclizable = CLOSING_SRU_NOT_APPLICABLE;
                }

                if (u->cyclizable != CLOSING_SRU_NOT_APPLICABLE)
                {
                    /* Allocate PS (frame-shiftable) bonds */
                    if (u->bkbonds)
                    {
                        imat_free(u->maxbkbonds, u->bkbonds);
                        u->bkbonds = NULL;
                    }
                    u->maxbkbonds = orig_at_data->num_inp_bonds + 2;
                    err = imat_new( u->maxbkbonds, 2, &( u->bkbonds ) );
                    if (err)
                    {
                        TREAT_ERR( err, 9034, "Not enough memory (polymers)" );
                        goto exit_function;
                    }
                }
            }

        }
    }
    else
    {
       TREAT_ERR( err, 9035, "Invalid kind of polymer representation" );
       goto exit_function;
    }

    orig_at_data->valid_polymer = 1;

exit_function:
    if (err)
    {
        orig_at_data->valid_polymer = 0;
    }

    return err;
}


/****************************************************************************/
int UnMarkRingSystemsInp( inp_ATOM *at, int num_atoms )
{
    int i;
    for (i = 0; i < num_atoms; i++)
    {
        at[i].bCutVertex = 0;
        at[i].nRingSystem = 0;
        at[i].nNumAtInRingSystem = 0;
        at[i].nBlockSystem = 0;
    }

    return 0;
}


/****************************************************************************
 Preprocess OAD_Polymer (NB: frame shift is invoked from here)
****************************************************************************/
int OAD_Polymer_CyclizeCloseableUnits( ORIG_ATOM_DATA *orig_at_data,
                                       int treat_polymers,
                                       char *pStrErr,
                                       int bNoWarnings )
{
    int i, err = 0; /* djb-rwth: removing redundant variables */

    for (i = 0; i < orig_at_data->polymer->n; i++)
    {
        OAD_PolymerUnit *unit = orig_at_data->polymer->units[i];

        if (!unit->cyclizable)
        {
            continue;
        }

        /* Find stars and their partners */
        OAD_PolymerUnit_SetEndsAndCaps( unit, orig_at_data, &err, pStrErr );
            /*	Reveal and store CRU caps and ends('stars and partners')
                Also set `unit->cap1_is_undef`, `unit->cap2_is_undef`, `unit->cyclizable` 
            */
        if (err)
        {
            break;
        }
        if (!unit->cyclizable)
        {
            continue;
        }


        if (OAD_PolymerUnit_HasMetal( unit, orig_at_data->at ))
        {
            /*unit->cyclizable = CLOSING_SRU_NOT_APPLICABLE;*/
            if (unit->cyclizable == CLOSING_SRU_RING)
            {
                /*unit->cyclizable = CLOSING_SRU_HIGHER_ORDER_BOND;*/
                if (!bNoWarnings)
                {
                    WarningMessage( pStrErr, "Frame shift in metallated polymer unit may be missed" );
                }
            }
        }

        /* Now remove bonds to cap ("star atoms and cyclize a SRU */
        OAD_PolymerUnit_UnlinkCapsAndConnectEndAtoms( unit, orig_at_data, &err, pStrErr );

        if (err)
        {
            break;
        }
        if (!unit->cyclizable)
        {
            continue;
        }

        /* djb-rwth: removing redundant code */
    }

    /*
    if ( ncyclized )
    {
        if (!bNoWarnings)
        {
            WarningMessage( pStrErr, "Made provision for frame shift in polymer unit(s)" );
        }
    }
    */

    return err;
}


/****************************************************************************
 Check if SRU contains metal
****************************************************************************/
int OAD_PolymerUnit_HasMetal( OAD_PolymerUnit *u, inp_ATOM *at )
{
    int i;
    for (i = 0; i < u->na; i++)
    {
        if (is_el_a_metal( at[u->alist[i] - 1].el_number ))
        {
            return 1;
        }
    }

    return 0;
}


/****************************************************************************/
void OAD_Polymer_Free( OAD_Polymer *pd )
{
    if (pd)
    {
        if (pd->pzz)
        {
            inchi_free( pd->pzz );
            pd->pzz = NULL;
            pd->n_pzz = 0;
        }
        if (pd->n && pd->units)
        {
            int k;
            for (k = 0; k < pd->n; k++)
            {
                OAD_PolymerUnit_Free( pd->units[k] );
            }
            inchi_free( pd->units );
            pd->units = NULL;
            pd->n = 0;
        }
        inchi_free( pd );
        pd = NULL;
    }

    return;
}


/****************************************************************************/
void OAD_PolymerUnit_UnlinkCapsAndConnectEndAtoms( OAD_PolymerUnit *unit,
                                                   ORIG_ATOM_DATA  *orig_inp_data,
                                                   int             *err,
                                                   char            *pStrErr )
{
    int bond_type, bond_stereo;

    *err = 0;
    if (!unit->cyclizable)
    {
        return;
    }

    if (unit->cyclizable == CLOSING_SRU_RING)
    {
        /* Disconnect both star atoms */
        OrigAtData_RemoveBond( unit->cap1 - 1, unit->end_atom1 - 1, orig_inp_data->at,
                               &bond_type, &bond_stereo, &orig_inp_data->num_inp_bonds );

        OrigAtData_RemoveBond( unit->cap2 - 1, unit->end_atom2 - 1, orig_inp_data->at,
                               &bond_type, &bond_stereo, &orig_inp_data->num_inp_bonds );

        OrigAtData_AddSingleStereolessBond( unit->end_atom1 - 1, unit->end_atom2 - 1,
                                            orig_inp_data->at, &orig_inp_data->num_inp_bonds );
    }

    else if (unit->cyclizable == CLOSING_SRU_HIGHER_ORDER_BOND)
    {
        int elevated; /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
        elevated = OrigAtData_IncreaseBondOrder( unit->end_atom1 - 1, unit->end_atom2 - 1, orig_inp_data->at ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
#if 0
/* the bond may already be broken at metal disconnection, so ignore the result here */
        if (!elevated)
        {
            /* *err = 1; */
            WarningMessage( pStrErr, "SRU closure via higher order bond failed" );
            unit->cyclizable = CLOSING_SRU_NOT_APPLICABLE;
            return;
        }
#endif
        OrigAtData_RemoveBond( unit->cap1 - 1, unit->end_atom1 - 1, orig_inp_data->at,
                               &bond_type, &bond_stereo, &orig_inp_data->num_inp_bonds );
        OrigAtData_RemoveBond( unit->cap2 - 1, unit->end_atom2 - 1, orig_inp_data->at,
                               &bond_type, &bond_stereo, &orig_inp_data->num_inp_bonds );
    }

    else if (unit->cyclizable == CLOSING_SRU_DIRADICAL)
    {
        orig_inp_data->at[unit->end_atom1 - 1].radical = RADICAL_TRIPLET;
        OrigAtData_RemoveBond( unit->cap1 - 1, unit->end_atom1 - 1, orig_inp_data->at,
                               &bond_type, &bond_stereo, &orig_inp_data->num_inp_bonds );
        OrigAtData_RemoveBond( unit->cap2 - 1, unit->end_atom2 - 1, orig_inp_data->at,
                                &bond_type, &bond_stereo, &orig_inp_data->num_inp_bonds );
    }

    if (!*err)
    {
        unit->cyclized = 1;
    }
    
    return;
}


/****************************************************************************/
void OAD_PolymerUnit_FindEndsAndCaps( OAD_PolymerUnit *unit,
                                      ORIG_ATOM_DATA *orig_at_data,
                                      int *end1,
                                      int *cap1,
                                      int *cap1_is_star,
                                      int *end2,
                                      int *cap2,
                                      int *cap2_is_star,
                                      int *err,
                                      char *pStrErr )
{
    int i, j, i_inside = 0, j_inside = 0;
    int num_atoms = orig_at_data->num_inp_atoms;

    *end1 = *end2 = *cap1 = *cap2 = 0;
    *cap1_is_star = *cap2_is_star = 0;
    *err = 0;

    if (!unit->blist || unit->nb < 1)
    {
        return;
    }
    /* Left crossing bond */
    i = unit->blist[0];
    j = unit->blist[1];
    i_inside = (NULL != is_in_the_ilist(unit->alist, i, unit->na));
    j_inside = (NULL != is_in_the_ilist(unit->alist, j, unit->na));
    if (i_inside && j_inside)
    {
        TREAT_ERR(*err, 9032, "Polymer CRU cap(s) lie inside CRU");
        return;
    }
    if (i_inside)
    {
        *end1 = i;
        *cap1 = j;
    }
    else
    {
        *end1 = j;
        *cap1 = i;
    }
    if (!strcmp(orig_at_data->at[*cap1 - 1].elname, "Zz"))
    {
        *cap1_is_star = 1;
    }
    /* Right crossing bond */
    i = unit->blist[2];
    j = unit->blist[3];
    i_inside = NULL != is_in_the_ilist(unit->alist, i, unit->na);
    j_inside = NULL != is_in_the_ilist(unit->alist, j, unit->na);
    if (i_inside && j_inside)
    {
        TREAT_ERR(*err, 9032, "Polymer CRU cap(s) lie inside CRU");
    }
    if (i_inside)
    {
        *end2 = i;
        *cap2 = j;
    }
    else
    {
        *end2 = j;
        *cap2 = i;
    }
    if (!strcmp(orig_at_data->at[*cap2 - 1].elname, "Zz"))
    {
        *cap2_is_star = 1;
    }
    /* Checks */
    if (*end1 <= 0 || *end1 > num_atoms || *cap1 <= 0 || *cap1 > num_atoms)
    {
        TREAT_ERR(*err, 9090, "Invalid polymer CRU crossing bond");
        return;
    }
    if (*end2 <= 0 || *end2 > num_atoms || *cap2 <= 0 || *cap2 > num_atoms)
    {
        TREAT_ERR(*err, 9091, "Invalid polymer CRU crossing bond");
        return;
    }
    if  ( *cap1 == *cap2 ) /* || (*end1 == *end2 && unit->na>1)  */
    {
        TREAT_ERR(*err, 9090, "Invalid polymer CRU surrounding");
        return;
    }

    /* Paranoia, 2020-05-22 */
    unit->end_atom1 = *end1;
    unit->end_atom2 = *end2;
    unit->cap1 = *cap1;
    unit->cap2 = *cap2;

    *err = 0;
    return;
}
    
    
/****************************************************************************
 Reveal and store CRU caps and ends ('stars and partners')
****************************************************************************/
void OAD_PolymerUnit_SetEndsAndCaps( OAD_PolymerUnit *unit,
                                     ORIG_ATOM_DATA  *orig_at_data,
                                     int             *err,
                                     char            *pStrErr )
{
    int k;

    unit->cyclizable = CLOSING_SRU_NOT_APPLICABLE;
    unit->end_atom1 = unit->end_atom2 = unit->cap1 = unit->cap2 = -1;
    unit->cap1_is_undef = unit->cap2_is_undef = 0;

    OAD_PolymerUnit_FindEndsAndCaps( unit, orig_at_data,
                                     &unit->end_atom1, &unit->cap1, &unit->cap1_is_undef,
                                     &unit->end_atom2, &unit->cap2, &unit->cap2_is_undef,
                                     err, pStrErr);

    if (*err)
    {
        goto exit_function;
    }

#if ( defined(DEBUG_POLYMERS) && ( DEBUG_POLYMERS != 0 ) )
    ITRACE_( "Cap-end_atom pairs (numbers are from 1) are: %-d-%-d and %-d-%-d\n",
        unit->cap1, unit->end_atom1, unit->cap2, unit->end_atom2 );
#endif

    if (!unit->cap1_is_undef && !unit->cap2_is_undef)
    {
        goto exit_function;
    }

    /* The rest is applicable only to *---SRU---* case */

    /* Stars are separated by one atom - that's not error but do nothing */
    if (unit->end_atom1 == unit->end_atom2)
    {
#ifdef ALLOW_CLOSING_SRU_VIA_DIRADICAL
        unit->cyclizable = CLOSING_SRU_DIRADICAL;
#else
        unit->cyclizable = CLOSING_SRU_NOT_APPLICABLE;
#endif
        goto exit_function;
    }

    /* Stars are separated by two atoms - that's not error but do nothing */
    for (k = 0; k < orig_at_data->at[unit->end_atom1 - 1].valence; k++)
    {
        if (orig_at_data->at[unit->end_atom1 - 1].neighbor[k] == unit->end_atom2 - 1)
        {
#ifdef ALLOW_CLOSING_SRU_VIA_HIGHER_ORDER_BOND
            unit->cyclizable = CLOSING_SRU_HIGHER_ORDER_BOND;
#else
            unit->cyclizable = CLOSING_SRU_NOT_APPLICABLE;
#endif
            goto exit_function;
        }
    }

    unit->cyclizable = CLOSING_SRU_RING;

exit_function:

    return;
}


/****************************************************************************
  Replace original atom numbers in polymer data with (canonical num + 1)
        Then prepare:
            units2    a copy of original polymer units (p->units) with
                      atomic numbers changed to curr canonical ones;
                      atoms in alists sorted; atoms in blists
                      and blists themselves are sorted
            unum      numbers of units (0..p->n) as they go when
                      sorted by alist's in lexicographic orders
****************************************************************************/
int OAD_Polymer_PrepareWorkingSet( OAD_Polymer     *p,
                                   int             *cano_nums,
                                   int             *compnt_nums,
                                   OAD_PolymerUnit **units2,       /* allocd by caller, to be filled */
                                   int             *unum )         /* allocd by caller, to be filled */

{
    int i, k, err = 0, cano_num1 = -1, cano_num2 = -1;
    OAD_PolymerUnit *u;

    /*OAD_Polymer_DebugTrace( p );*/

    /*  Replace original atom numbers in polymer data with canonical plus 1.                    */
    /*  Note that we use 'cano1 nums', that is, 1-based (InChI internal 'cano nums' are 0-based)*/
    /*  Also remove from the list atoms who mapped to cano number 0  ( == -1 + 1_offset ),      */
    /*  they are explicit H's which have already been deleted.                                  */
    for (k = 0; k < p->n_pzz; k++)
    {
        cano_num1 = cano_nums[p->pzz[k]] + 1;
        if (cano_num1 == 0)
        {
            /* we shouldn't arrive here */
            err = 10;
            goto exit_function;
        }
        p->pzz[k] = cano_num1;
    }

    for (i = 0; i < p->n; i++)
    {
        int na_new = -1;
        u = units2[i];

        for (k = 0; k < u->na; k++)
        {
            cano_num1 = cano_nums[u->alist[k]] + 1;
            if (cano_num1 == 0)
            {
                continue;
            }
            u->alist[++na_new] = cano_num1;
        }
        u->na = na_new + 1;
        for (k = 0; k < 2 * u->nb; k++)
        {
            cano_num1 = cano_nums[u->blist[k]] + 1;
            if (cano_num1 == 0)
            {
                /* Can not proceed further as one of PU crossing bond ends
                   leads to explicit H which has been removed already       */
                err = 11;
                goto exit_function;
            }
            u->blist[k] = cano_num1;
        }

        cano_num1 = cano_nums[u->cap1] + 1;
        if (cano_num1 == 0)
        {
            err = 11;
            goto exit_function;
        }
        u->cap1 = cano_num1;

        cano_num1 = cano_nums[u->cap2] + 1;
        if (cano_num1 == 0)
        {
            err = 11;
            goto exit_function;
        }
        u->cap2 = cano_num1;

        cano_num1 = cano_nums[u->end_atom1] + 1;
        if (cano_num1 == 0)
        {
            err = 11;
            goto exit_function;
        }
        u->end_atom1 = cano_num1;

        cano_num1 = cano_nums[u->end_atom2] + 1;
        if (cano_num1 == 0)
        {
            err = 11;
            goto exit_function;
        }
        u->end_atom2 = cano_num1;

        for (k = 0; k < u->nbkbonds; k++)
        {
            cano_num1 = cano_nums[u->bkbonds[k][0]] + 1;
            if (cano_num1 == 0)
            {
                continue;
            }
            cano_num2 = cano_nums[u->bkbonds[k][1]] + 1;
            if (cano_num2 == 0)
            {
                continue;
            }
            u->bkbonds[k][0] = inchi_min( cano_num1, cano_num2 );
            u->bkbonds[k][1] = inchi_max( cano_num1, cano_num2 );
        }
    }

    /* Sort the atoms and the bonds in all units */
    for (i = 0; i < p->n; i++)
    {
        u = units2[i];

        /* sort atoms (alist) */
        iisort( u->alist, u->na );

        /*ITRACE_( "\n*** Polymer unit %-d : ( ", i );
        for (k = 0; k < u->na - 1; k++)
        {
            ITRACE_( "%-d-", u->alist[k] );
        }
        ITRACE_( "%-d )\n", u->alist[u->na - 1] );*/

        /* Sort bonds (blist) */
        err = OAD_PolymerUnit_OrderBondAtomsAndBondsThemselves( u, p->n_pzz, p->pzz );
        if (err)
        {
            /* crossing bonds in blist are invalid */
            err = 12;
            goto exit_function;
        }

        /* Check each unit for >1 connected components */
#if 0
        {
            int icompnt;
            icompnt = compnt_nums[u->alist[0] - 1];
            for (k = 1; k < u->na; k++)
            {
                if (compnt_nums[u->alist[k] - 1] != icompnt)
                {
                    u->disjoint = 1;
                    break;
                }
            }
    }
#endif

    }

    /* Sort all units in modified alist's lexicographic order 
    (modification is: longer list always go first )			*/
    for (i = 0; i < p->n; i++)
    {
        unum[i] = i;
    }
    for (i = 1; i < p->n; i++)
    {
        int tmp = unum[i];
        int j = i - 1;
        while (j >= 0 && OAD_PolymerUnit_CompareAtomListsMod( units2[unum[j]], units2[tmp] ) > 0)
        /*while ( j >= 0 &&    OAD_PolymerUnit_CompareAtomLists( units2[ unum[j] ], units2[ tmp ] ) > 0  )*/
        {
            unum[j + 1] = unum[j];
            j--;
        }
        unum[j + 1] = tmp;
    }

exit_function:

    return err;
}


/****************************************************************************
  Helper for cyclizing CRU. NB: 0-based
****************************************************************************/
int  OrigAtData_RemoveHalfBond( int      this_atom,
                                int      other_atom,
                                inp_ATOM *at,
                                int      *bond_type,
                                int      *bond_stereo )
{
    int k, kk;
    /* djb-rwth: fixing oss-fuzz issues #68286, #30342 */
    if (at && (this_atom >= 0) && (other_atom >= 0))
    {
        inp_ATOM* a = &(at[this_atom]);
        if (a)
        {
            for (k = 0; k < a->valence; k++)
            {
                if (a->neighbor[k] != other_atom)
                {
                    continue;
                }

                *bond_type = a->bond_type[k];
                *bond_stereo = a->bond_stereo[k];

                a->neighbor[k] = a->bond_type[k] = a->bond_stereo[k] = 0;

                for (kk = k + 1; kk < a->valence; kk++)
                {
                    a->neighbor[kk - 1] = a->neighbor[kk];
                    a->bond_type[kk - 1] = a->bond_type[kk];
                    a->bond_stereo[kk - 1] = a->bond_stereo[kk];
                }
                for (kk = a->valence - 1; kk < MAXVAL; kk++)
                {
                    a->neighbor[kk] = 0;
                    a->bond_type[kk] = (U_CHAR)0;
                    a->bond_stereo[kk] = (S_CHAR)0;
                }
                return 1;
            } /* k */
        }
    }

    return 0;
}


/****************************************************************************/
int  OrigAtData_RemoveAtom(ORIG_ATOM_DATA *orig_at_data, int iatom)
{

    if (0)
    {
        return 1;
    }

    return 0;
}


/****************************************************************************/
int  OrigAtData_RemoveBond( int      this_atom,
                            int      other_atom,
                            inp_ATOM *at,
                            int      *bond_type,
                            int      *bond_stereo,
                            int      *num_inp_bonds )
{
    int del = 0;
    
    if (at && (this_atom >= 0) && (other_atom >= 0)) /* djb-rwth: fixing oss-fuzz issue #68329, #68286 */
    {
        del = OrigAtData_RemoveHalfBond(this_atom, other_atom, at, bond_type, bond_stereo);
        del += OrigAtData_RemoveHalfBond(other_atom, this_atom, at, bond_type, bond_stereo);

        if (del == 2)
        {
            (*num_inp_bonds)--;
            at[this_atom].valence--;
            at[this_atom].chem_bonds_valence -= *bond_type;
            at[other_atom].valence--;
            at[other_atom].chem_bonds_valence -= *bond_type;
            return 1;
        }
    }

    return 0;
}


/****************************************************************************/
int  OrigAtData_AddBond( int        this_atom,
                         int        other_atom,
                         inp_ATOM   *at,
                         int        bond_type,
                         int        bond_stereo,
                         int        *num_bonds )
{
    if (at)
    {
        /* djb-rwth: fixing oss-fuzz issue #68286 */
        int i, k, already_here;
        inp_ATOM* a = &(at[this_atom]);

        if (at[this_atom].valence >= MAXVAL ||
            at[other_atom].valence >= MAXVAL)
        {
            return 0;
        }

        if (bond_type != INCHI_BOND_TYPE_DOUBLE && bond_type != INCHI_BOND_TYPE_TRIPLE)
        {
            bond_type = INCHI_BOND_TYPE_SINGLE;
        }

        k = a->valence;
        already_here = 0;
        for (i = 0; i < k; i++)
        {
            if (a->neighbor[i] == other_atom)
            {
                already_here = 1; break;
            }
        }

        if (!already_here)
        {
            a->neighbor[k] = other_atom;
            a->bond_type[k] = (U_CHAR)bond_type;
            a->bond_stereo[k] = (S_CHAR)bond_stereo;
            a->chem_bonds_valence += bond_type;
            a->valence++;
        }

        a = &(at[other_atom]);
        k = a->valence;
        already_here = 0;
        for (i = 0; i < k; i++)
        {
            if (a->neighbor[i] == this_atom)
            {
                already_here = 1; break;
            }
        }

        if (!already_here && (k < MAXVAL)) /* djb-rwth: condition added to prevent buffer overrun */
        {
            a->neighbor[k] = this_atom;
            a->bond_type[k] = (U_CHAR)bond_type;
            a->bond_stereo[k] = (S_CHAR)bond_stereo;
            a->chem_bonds_valence += bond_type;
            a->valence++;
        }

        (*num_bonds)++;

        return 1;
    }
    else
    {
        return 0;
    }
}


/****************************************************************************/
int  OrigAtData_AddSingleStereolessBond( int      this_atom,
                                         int      other_atom,
                                         inp_ATOM *at,
                                         int      *num_bonds )
{
    return OrigAtData_AddBond( this_atom, other_atom, at, INCHI_BOND_TYPE_SINGLE, 0, num_bonds );
}


/****************************************************************************/
int  OrigAtData_IncreaseBondOrder( int this_atom, int other_atom, inp_ATOM *at )
{
    int i, k, n_up = 0;
    inp_ATOM *a;

    if (at[this_atom].valence >= MAXVAL ||
         at[other_atom].valence >= MAXVAL)
    {
        return 0;
    }

    a = &( at[this_atom] );
    if (a->chem_bonds_valence > MAXVAL - 1)
    {
        return 0;
    }

    k = a->valence;
    for (i = 0; i < k; i++)
    {
        if (a->neighbor[i] != other_atom)
            continue;
        if (a->bond_type[i] > 3)
            return 0;
        a->bond_type[i]++;
        a->chem_bonds_valence++;
        n_up++;
        break;
    }

    a = &( at[other_atom] );
    if (a->chem_bonds_valence > MAXVAL - 1)
    {
        return 0;
    }
    k = a->valence;

    for (i = 0; i < k; i++)
    {
        if (a->neighbor[i] != this_atom)
        {
            continue;
        }
        if (a->bond_type[i] > 3)
        {
            return 0;
        }
        a->bond_type[i]++;
        a->chem_bonds_valence++;
        n_up++;
        break;
    }

    return n_up;
}


/****************************************************************************/
int  OrigAtData_DecreaseBondOrder( int      this_atom,
                                   int      other_atom,
                                   inp_ATOM *at )
{
    int i, k, n_dn = 0;
    inp_ATOM *a;

    a = &( at[this_atom] );
    if (a->chem_bonds_valence > MAXVAL - 1)
    {
        return 0;
    }

    k = a->valence;
    for (i = 0; i < k; i++)
    {
        if (a->neighbor[i] != other_atom)
        {
            continue;
        }
        if (a->bond_type[i] < 2)
        {
            return 0;
        }
        a->bond_type[i]--;
        a->chem_bonds_valence--;
        n_dn++;
        break;
    }

    a = &( at[other_atom] );
    k = a->valence;
    for (i = 0; i < k; i++)
    {
        if (a->neighbor[i] != this_atom)
        {
            continue;
        }
        if (a->bond_type[i] < 2)
        {
            return 0;
        }
        a->bond_type[i]--;
        a->chem_bonds_valence--;
        n_dn++;
        break;
    }

    return n_dn;
}


/****************************************************************************
 Collect bonds and optionally atoms of fragment
****************************************************************************/
void OAD_CollectFragmentBondsAndAtoms(	ORIG_ATOM_DATA  *orig_at_data,
                                        int nforbidden,			/* number of edges forbidden for traversal		*/
                                        int *forbidden_orig,	/* atom nums of forbidden edges					*/
                                                                /* [edge1at1,edge1at2, edge2at1, edge2at2, ... ] */
                                        int *n_fragbonds,
                                        int **fragbonds,
                                        int *n_fragatoms,
                                        int *fragatoms,
                                        int *err,
                                        char *pStrErr)
{
    int i;
    int	max_atoms = orig_at_data->num_inp_atoms;
    int *atnums = NULL;
    subgraf *sg = NULL;
    subgraf_pathfinder *spf = NULL;

    *err = 0;
    atnums = (int *)inchi_calloc(max_atoms, sizeof(int));
    if (!atnums)
    {
        TREAT_ERR(*err, 9045, "Not enough memory");
        goto exit_function;
    }
    for (i = 0; i < max_atoms; i++)
    {
        atnums[i] = orig_at_data->at[i].orig_at_number; /* i+1 normally*/
    }

    sg = subgraf_new(orig_at_data, max_atoms, atnums);
    if (!sg)
    {
        TREAT_ERR(*err, 9045, "Not enough memory");
        goto exit_function;
    }
    spf = subgraf_pathfinder_new(sg, orig_at_data, 0, 0); /* start = end = 0th node */
    if (!spf)
    {
        TREAT_ERR(*err, 9045, "Not enough memory");
        goto exit_function;
    }

    spf->seen[0] = spf->start; 
    spf->nseen = 1;
    *n_fragbonds = 0;
    *n_fragatoms = 0;

    subgraf_pathfinder_run(	spf, 
                            nforbidden, forbidden_orig, /* this corrects cinnectivity of subgraf... */
                            n_fragbonds, fragbonds, 
                            n_fragatoms, fragatoms);


exit_function:
    subgraf_free(sg);
    subgraf_pathfinder_free(spf);
    if (atnums)
    {
        inchi_free(atnums);
    }
    return;
}


/****************************************************************************
 For all CRUs detect the bonds potentially involved in frame shift
****************************************************************************/
void OAD_Polymer_FindBackbones( ORIG_ATOM_DATA *at_data,
                                COMP_ATOM_DATA *composite_norm_data,
                                int            *err,
                                char           *pStrErr )
{
    int i;

    *err = 0;
    for (i = 0; i < at_data->polymer->n; i++)
    {
        if (!at_data->polymer->units[i]->cyclizable)
        {
            continue;
        }

        OAD_CollectBackboneBonds( at_data, 
                                  at_data->polymer->units[i]->na,
                                  at_data->polymer->units[i]->alist,
                                  at_data->polymer->units[i]->end_atom1,
                                  at_data->polymer->units[i]->end_atom2,
                                  &(at_data->polymer->units[i]->nbkbonds),
                                  at_data->polymer->units[i]->bkbonds,
                                  err, pStrErr );
        if (*err)
        {
            at_data->polymer->units[i]->cyclizable = CLOSING_SRU_NOT_APPLICABLE;
            continue;
        }
        if (at_data->polymer->units[i]->nbkbonds < 1)
        {
            continue;
        }
        if (at_data->polymer->units[i]->nbkbonds == 1)
        {
            /*  Special case: we got only one bond between end_atom1 and   */
            /*  end_atom2 (this may be the result of metal disconnection)  */
            continue;
        }

        OAD_PolymerUnit_DelistIntraRingBackboneBonds( at_data->polymer->units[i], at_data, err, pStrErr );
        if (*err)
        {
            continue;
        }

        OAD_PolymerUnit_DelistHighOrderBackboneBonds( at_data->polymer->units[i],
                                                      at_data, composite_norm_data,
                                                      err, pStrErr );
        if (*err)
        {
            continue;
        }
        if (at_data->polymer->units[i]->nbkbonds == 0)
        {
            /* We already cyclized frame-shiftable unit and preprocessed it (in 'prep_inp_data').                       */
            /* Despite that, now we discovered that there are no bonds eligible for frame shift                         */
            /* (as all potentially eligible in-unit bonds are either in-ring or alternate ones).                        */
            /* We can not simply restore original connections as the structure may have been already heavily touched.   */
            /* The most viable action is to hold a single frame-shift bond (between original caps of CRU ends).   */
            /* It is for sure will be converted to original bonds to star atoms on possible inchi2struct.               */
            at_data->polymer->units[i]->cyclizable = 1;
            at_data->polymer->units[i]->nbkbonds = 1;
            at_data->polymer->units[i]->bkbonds[0][0] = at_data->polymer->units[i]->end_atom1;
            at_data->polymer->units[i]->bkbonds[0][1] = at_data->polymer->units[i]->end_atom2;
        }
    }

    return;
}


/****************************************************************************
Collect all backbone atoms - of main chain(s), side chains being ignored
****************************************************************************/
void OAD_CollectBackboneAtoms(ORIG_ATOM_DATA  *at_data,
                              int na,
                              int *alist,
                              int end_atom1,
                              int end_atom2,
                              int *nbkatoms,
                              int *bkatoms,
                              int *err,
                              char *pStrErr)
{
    int start = 0, end = 0;
    subgraf *sg = NULL;
    subgraf_pathfinder *spf = NULL;
    int nbkbonds=0;
    int **bkbonds=NULL;		/* list of [breakable] backbone bonds	[ (a1,a2), (a3,a4), ... ]	*/
    int maxbkbonds;

    *nbkatoms = 0;
    maxbkbonds = at_data->num_inp_bonds + 2;
    *err = imat_new(maxbkbonds, 2, &(bkbonds));
    /* djb-rwth: addressing coverity ID #499570 -- TREAT_ERR properly used in all cases */
    if (*err)
    {
        TREAT_ERR(*err, 9034, "Not enough memory (polymers)");
        goto exit_function;
    }
    nbkbonds = 0;
    sg = subgraf_new(at_data, na, alist);
    if (!sg)
    {
        TREAT_ERR(*err, 9037, "Not enough memory (polymers)");
        goto exit_function;
    }

    start = sg->orig2node[end_atom1]; end = sg->orig2node[end_atom2];
    if (start > end)
    {
        int tmp = end;
        end = start;
        start = tmp;
    }
    spf = subgraf_pathfinder_new(sg, at_data, start, end);
    if (!spf)
    {
        TREAT_ERR(*err, 9039, "Not enough memory (polymers)");
        goto exit_function;
    }

    spf->seen[0] = spf->start; spf->nseen = 1;    
    nbkbonds = 0;
    *nbkatoms = 0;
    
    subgraf_pathfinder_run(spf, 0, NULL, &nbkbonds, bkbonds, nbkatoms, bkatoms);

    subgraf_free(sg);
    subgraf_pathfinder_free(spf);
    *err = 0;


exit_function:
    if (bkbonds)
    {
        imat_free(maxbkbonds, bkbonds);
        bkbonds = NULL;
    }	

    return;
}


/****************************************************************************
Collect all atoms reachable from start_atom
****************************************************************************/
int OAD_CollectReachableAtoms( ORIG_ATOM_DATA  *orig_at_data,
                               int start_atom,
                               int nforbidden_bonds,
                               int *forbidden_bond_atoms,
                               int *n_reachable,
                               int *reachable,
                               int *err,
                               char *pStrErr)
{
    int iatom, natnums, max_atoms, j, ret = _IS_OKAY;
    int max_n_reachable = *n_reachable;
    int *atnums = NULL;
    subgraf *sg = NULL;
    subgraf_pathfinder *spf = NULL;

    /* djb-rwth: removing redundant code */
    max_atoms = orig_at_data->num_inp_atoms;
    iatom = start_atom - 1;
    *n_reachable = 0;

    atnums = (int *)inchi_calloc(max_atoms, sizeof(int));
    if (!atnums)
    {
        ret = _IS_ERROR;
        goto exit_function;
    }
    for (j = 0; j < max_atoms; j++)
    {
        atnums[j] = orig_at_data->at[j].orig_at_number;
    }
    sg = subgraf_new(orig_at_data, max_atoms, atnums);
    if (!sg)
    {
        ret = _IS_ERROR;
        goto exit_function;
    }
    spf = subgraf_pathfinder_new(sg, orig_at_data, iatom, iatom);
    if (!spf)
    {
        ret = _IS_ERROR;
        goto exit_function;
    }

    /* move from orig# to node# */
    spf->start = iatom;
    for (j = 0; j < nforbidden_bonds; j++)
    {
        forbidden_bond_atoms[2 * j]      = sg->orig2node[forbidden_bond_atoms[2 * j] ];
        forbidden_bond_atoms[2 * j + 1]  = sg->orig2node[forbidden_bond_atoms[2 * j + 1] ];
    }

    /*memset(atnums, -1, max_atoms * sizeof(int));*/
    for (j = 0; j < max_atoms; j++)
    {
        atnums[j] = -1;
    }

    spf->nseen = 0;
    natnums = subgraf_pathfinder_collect_all(spf, nforbidden_bonds, forbidden_bond_atoms, atnums);

    if (natnums)
    {
        if (natnums > max_n_reachable)
        {
            ret = _IS_ERROR;
            goto exit_function;
        }

        for (j = 0; j < natnums && j < max_atoms; j++) /* djb-rwth: fixing buffer overruns */
        {
            reachable[(*n_reachable)++ ] = atnums[j];
        }
    }

exit_function:
    subgraf_free(sg);
    subgraf_pathfinder_free(spf);
    if (atnums)
    {
        inchi_free(atnums);
    }
    
    return ret;
}


/****************************************************************************
 Collect all backbone bonds - of main chain(s), side chains being ignored
 (for polymer CRU, these are the bonds potentially involved in frame shift)
****************************************************************************/
void OAD_CollectBackboneBonds(ORIG_ATOM_DATA  *at_data,
                              int na, 
                              int *alist,
                              int end_atom1,
                              int end_atom2,
                              int *nbkbonds,
                              int **bkbonds,
                              int *err,
                              char *pStrErr )
{
    int start = 0, end = 0, dummy;
    subgraf *sg = NULL;
    subgraf_pathfinder *spf = NULL;
    /* Establish subgraph for na atoms of the alist */
    *nbkbonds = 0;
    sg = subgraf_new( at_data, na, alist );
    if (!sg)
    {
        TREAT_ERR( *err, 9037, "Not enough memory (polymers)" );
        /* unit->cyclizable = CLOSING_SRU_NOT_APPLICABLE; */
        return;
    }
    start = sg->orig2node[end_atom1]; 
    end = sg->orig2node[end_atom2];
#if 0
    if (start > end)
    {
        int tmp = end;
        end = start;
        start = tmp;
    }
#endif
    spf = subgraf_pathfinder_new( sg, at_data, start, end );
    if (!spf)
    {
        TREAT_ERR( *err, 9039, "Not enough memory (polymers)" ); /* djb-rwth: addressing coverity ID #499539 -- TREAT_ERR properly used */
        /*unit->cyclizable = CLOSING_SRU_NOT_APPLICABLE;*/
        return;
    }
    spf->seen[0] = spf->start; 
    spf->nseen = 1;    
    *nbkbonds = 0;
    subgraf_pathfinder_run( spf, 0, NULL,
                            nbkbonds,
                            bkbonds, /* we will collect backbone CRU bonds here */
                            &dummy,
                            NULL );

    subgraf_free( sg );
    subgraf_pathfinder_free( spf );
    *err = 0;

    return;
}


/****************************************************************************
 Remove intra-ring bonds from the list of frame-shiftable bonds
****************************************************************************/
void OAD_PolymerUnit_DelistIntraRingBackboneBonds( OAD_PolymerUnit *unit,
                                                   ORIG_ATOM_DATA  *at_data,
                                                   int             *err,
                                                   char            *pStrErr )
{
    int nrings = 0;
    int *num_ring_sys = NULL;

    if (!unit)
    {
        return;
    }
    if (unit->nbkbonds < 1)
    {
        return;
    }

    /* Establish ring systems assignments for all related atoms */

    *err = 1;
    num_ring_sys = (int *) inchi_calloc( (long long)at_data->num_inp_atoms + 1, sizeof( int ) ); /* djb-rwth: cast operator added */
    if (!num_ring_sys)
    {
        goto exit_function;
    }

    *err = 0;

    nrings = OAD_Polymer_FindRingSystems( at_data->polymer, at_data->at, at_data->num_inp_atoms, &at_data->num_inp_bonds,
                                          num_ring_sys, NULL, unit->end_atom1 - 1 ); /* NB: start dfs within connected compt! */

    if (nrings == 0)
    {
        goto exit_function;
    }
    else
    {
        int at1, at2, j = 0;
repeatj:
        at1 = unit->bkbonds[j][0];
        at2 = unit->bkbonds[j][1];
        if (( num_ring_sys[at1] == num_ring_sys[at2] ) && ( num_ring_sys[at1] != -1 ))
        {
            OAD_PolymerUnit_RemoveLinkFromCRUChain( at1, at2, &unit->nbkbonds, unit->bkbonds );
        }
        else
        {
            ++j;
        }
        if (j < unit->nbkbonds)
        {
            goto repeatj;
        }
    }

exit_function:
    if (num_ring_sys)
    {
        inchi_free( num_ring_sys );
    }

    return;
}


/****************************************************************************
 Find ring systems (exclude possible cyclizing bonds) in all polymer SRU's
****************************************************************************/
int  OAD_Polymer_FindRingSystems( OAD_Polymer *pd,
                                  inp_ATOM    *at,
                                  int         nat,
                                  int         *num_inp_bonds,
                                  int         *num_ring_sys,
                                  int         *size_ring_sys,
                                  int         start )
{
    int i, j, nrings = 0, bond_type, bond_stereo;

    if (NULL == num_ring_sys)
    {
        return 0;
    }

    /* Remove polymer SRU 'cyclizing' bonds if any */
    for (j = 0; j < pd->n; j++)
    {
        if (pd->units[j]->cyclized)
        {
            OrigAtData_RemoveBond( pd->units[j]->end_atom1 - 1,
                                   pd->units[j]->end_atom2 - 1,
                                   at, &bond_type, &bond_stereo,
                                   num_inp_bonds );
        }
    }

    MarkRingSystemsInp( at, nat, start ); /*0 );*/

    for (i = 0; i <= nat; i++)
    {
        num_ring_sys[i] = -1;
    }
    for (i = 0; i < nat; i++)
    {
        if (at[i].nNumAtInRingSystem > 2)
        {
            int atnum = at[i].orig_at_number;
            num_ring_sys[atnum] = at[i].nRingSystem;
            if (NULL != size_ring_sys)
            {
                size_ring_sys[atnum] = at[i].nNumAtInRingSystem;
            }
        }
    }

    UnMarkRingSystemsInp( at, nat );

    for (i = 0; i < nat; i++)
    {
        if (num_ring_sys[i] > -1)
        {
            nrings++;
        }
    }

    /* Restore polymer SRU 'cyclizing' bonds if applicable */
    for (j = 0; j < pd->n; j++)
    {
        if (pd->units[j]->cyclized)
        {
            OrigAtData_AddSingleStereolessBond( pd->units[j]->end_atom1 - 1,
                                                pd->units[j]->end_atom2 - 1,
                                                at,
                                                num_inp_bonds );
        }
    }

    return nrings;
}


/****************************************************************************
  Fill atomic properties (OrigAtData ) necessary
  to calc seniority in polymer SRUs
****************************************************************************/
void OAD_Polymer_SetAtProps( OAD_Polymer *pd,
                             inp_ATOM    *at,
                             int         nat,
                             int         *num_inp_bonds,
                             OAD_AtProps *aprops,
                             int         *cano_nums )
{
/*  Max rank for in-ring atom is 216 which is achieved for N (element number 7 in Periodic system & erank_rule2[] ),*/
/*    then goes O with rank 215 (element number 8), and so on... lowest rank is 1 for H .                           */
/*                                                                                                                  */
/*  This follows to IUPAC rule 2 [Pure Appl. Chem., Vol. 74, No. 10, 2002, p. 1926] which states:                   */
/*  a. a ring or ring system containing nitrogen;                                                                   */
/*  b. a ring or ring system containing the heteroatom occurring earliest in the order given in Rule 4;             */
/*  ( which is     O > S > Se > Te > N > P > As > Sb > Bi > Si > Ge > Sn > Pb > B > Hg )                            */
/*  ...                                                                                                             */

    int erank_rule2[] = { 0,1,198,197,196,202,2,216,215,191,190,189,188,187,206,210,214,183,182,181,180,179,178,177,176,
                          175,174,173,172,171,170,169,205,209,213,165,164,163,162,161,160,159,158,157,156,155,154,153,152,
                          151,204,208,212,147,146,145,144,143,142,141,140,139,138,137,136,135,134,133,132,131,130,129,128,
                          127,126,125,124,123,122,121,201,119,203,207,116,115,114,113,112,111,110,109,108,107,106,105,104,
                          103,102,101,100,99,98,97,96,95,94,93,92,91,90,89,88,87,86,85,84,83,82,81 };



    /*  Max rank for chain atom is 215 which is achieved for O (element number 8 in Periodic system & erank_rule4[] ),  */
    /*  then goes N with rank 212 (element number 8), and so on... lowest rank is 1 for H .                             */
    /*                                                                                                                  */
    /*  This follows to IUPAC rule 4 [Pure Appl. Chem., Vol. 74, No. 10, 2002, p. 1927] which states:                   */
    /*  O > S > Se > Te > N > P > As > Sb > Bi > Si > Ge > Sn > Pb > B > Hg                                             */
    /*  Note: Other heteroatoms may be placed within this order as indicated by their positions in the                  */
    /*  periodic table [5].                                                                                             */

    int erank_rule4[] = { 0,1,198,197,196,202,2,211,215,191,190,189,188,187,206,210,214,183,182,181,180,179,178,177,176,
                          175,174,173,172,171,170,169,205,209,213,165,164,163,162,161,160,159,158,157,156,155,154,153,152,
                          151,204,208,212,147,146,145,144,143,142,141,140,139,138,137,136,135,134,133,132,131,130,129,128,
                          127,126,125,124,123,122,121,201,119,203,207,116,115,114,113,112,111,110,109,108,107,106,105,104,
                          103,102,101,100,99,98,97,96,95,94,93,92,91,90,89,88,87,86,85,84,83,82,81 };


    int i, j, k, nrings = 0;
    int a1, a2, dummy = 0, bond_type = -1, bond_stereo = -1;
    int *num_ring_sys = NULL, *size_ring_sys = NULL;
    /* djb-rwth: fixing oss-fuzz issue #68112 */
    int err2_len = sizeof(erank_rule2) / sizeof(erank_rule2[0]);
    int err4_len = sizeof(erank_rule4) / sizeof(erank_rule4[0]);

    if ((NULL == aprops) || !at || !pd) /* djb-rwth: fixing oss-fuzz issue #68329, #68286 */
    {
        return;
    }

    /* Establish element ranks for atoms */
    for (k = 0; k < nat; k++)
    {
        int atnum = at[k].orig_at_number, index = k;
        U_CHAR err4_ind = at[k].el_number;
        if (cano_nums)
        {
            index = cano_nums[atnum];
        }
        if (index >= 0 && err4_ind < err4_len)
        {
            aprops[index].erank = erank_rule4[err4_ind];
            aprops[index].ring_erank = 0;
            aprops[index].ring_size = 0;
            aprops[index].ring_num = -1;
        }
        else
        {
            /* deleted H's go here */
            ;
        }

    }

    /* Establish ring systems assignments for atoms */
    num_ring_sys = (int *) inchi_calloc( (long long)nat + 1, sizeof( int ) ); /* djb-rwth: cast operator added */
    if (NULL == num_ring_sys)
    {
        goto exit_function;
    }
    size_ring_sys = (int *) inchi_calloc( (long long)nat + 1, sizeof( int ) ); /* djb-rwth: cast operator added */
    if (NULL == size_ring_sys)
    {
        goto exit_function;
    }

    /* Note that we get here on the way of InChI2Struct conversion.            */
    /* Break temporarily any of (actually, the first) SRU 'cyclizing' bonds    */
    for (j = 0; j < pd->n; j++)
    {
        if (pd->units[j]->na > 2 && pd->units[j]->nbkbonds > 0 &&
             pd->units[j]->cyclized == 0 &&
             pd->units[j]->cyclizable == CLOSING_SRU_RING)
        {
            a1 = pd->units[j]->bkbonds[0][0] - 1;
            a2 = pd->units[j]->bkbonds[0][1] - 1;
            OrigAtData_RemoveBond( a1, a2, at, &bond_type, &bond_stereo, &dummy );
        }
    }

    nrings = OAD_Polymer_FindRingSystems( pd, at, nat, num_inp_bonds, num_ring_sys, size_ring_sys, 0 );

    /* Immediately restore just broken bond(s) */
    for (j = 0; j < pd->n; j++)
    {
        if (pd->units[j]->na > 2 &&
             pd->units[j]->nbkbonds > 0 &&
             pd->units[j]->cyclized == 0 &&
             pd->units[j]->cyclizable == CLOSING_SRU_RING)
        {
            a1 = pd->units[j]->bkbonds[0][0] - 1;
            a2 = pd->units[j]->bkbonds[0][1] - 1;
            /* OrigAtData_AddSingleStereolessBond( a1, a2, at, &dummy ); */
            OrigAtData_AddBond( a1, a2, at, bond_type, bond_stereo, &dummy );
        }
    }

    if (nrings)
    {
        int max_ring_num = 0;
        /* SRU contains ring[s], proceed with them following (not totally) the IUPAC guidelines */
        for (k = 0; k < nat; k++)
        {
            /* Browse 0-based original atoms, go to 1-based cano nums domain if cano_nums mapping is suppied */
            int atnum = at[k].orig_at_number, index = k;
            if (cano_nums)
            {
                index = cano_nums[atnum] + 1;
            }
            if (num_ring_sys[atnum] >= 0)
            {
                aprops[index].ring_num = num_ring_sys[atnum];  /* temporarily */

                if (max_ring_num < aprops[index].ring_num)
                    max_ring_num = aprops[index].ring_num;          /* NB: OAD_Polymer_FindRingSystems may return num_ring_sys[]  */
                                                                    /* which is not a list of consecutive numbers                 */

                aprops[index].ring_size = size_ring_sys[atnum];     /* Size of ring system which includes the atom k .            */
                                                                    /* It is used as an additional score for in-ring              */
                                                                    /* atoms' prioritizing (instead of criteria in                */
                                                                    /* 2c-2h of IUPAC rule 2 which deal with ring sizes).         */
            }
        }

        for (i = 0; i <= max_ring_num; i++)
        {
            int erank, max_erank = 0;
            for (k = 0; k < nat; k++)
            {
                int atnum = at[k].orig_at_number, index = k;
                if (cano_nums)
                {
                    index = cano_nums[atnum] + 1;
                }
                if (aprops[index].ring_num == i)
                {
                    erank = erank_rule2[at[k].el_number];
                    if (erank > max_erank)
                        max_erank = erank;
                }
            }
            for (k = 0; k < nat; k++)
            {
                int atnum = at[k].orig_at_number, index = k;
                if (cano_nums)
                {
                    index = cano_nums[atnum] + 1;
                }
                if (aprops[index].ring_num == i)
                {
                    if (aprops[index].ring_size > 2)
                    {
                        aprops[index].ring_erank = max_erank;
                    }
                }
            }
        }
    }

exit_function:
    if (num_ring_sys)
    {
        inchi_free( num_ring_sys );
    }
    if (size_ring_sys)
    {
        inchi_free( size_ring_sys );
    }

    return;
}


/****************************************************************************
  Exclude higher order bonds from list of bkbonds
****************************************************************************/
void OAD_PolymerUnit_DelistHighOrderBackboneBonds( OAD_PolymerUnit *unit,
                                                   ORIG_ATOM_DATA  *orig_at_data,
                                                   COMP_ATOM_DATA  *composite_norm_data,
                                                   int             *err,
                                                   char            *pStrErr )
{
    int at1, at2, j = 0, k, check_taut = 0, remove; /* djb-rwth: removing redundant variables/code */

    int *orig_num = NULL, *curr_num = NULL;
    int bond_is_untouchable = 0, btype;  /* DT: moved from below 2024-09-01 DT */

    if (unit->na < 2)
    {
        return;
    }
    if (unit->nb < 2)
    {
        return;
    }
    if (unit->nbkbonds < 1)
    {
        return;
    }

    /* Care on tautomeric bonds */
    if (composite_norm_data)
    {
        check_taut = 1;
        orig_num = (int *) inchi_calloc( (long long)orig_at_data->num_inp_atoms + 2, sizeof( int ) ); /* djb-rwth: cast operator added */
        curr_num = (int *) inchi_calloc( (long long)orig_at_data->num_inp_atoms + 2, sizeof( int ) ); /* djb-rwth: cast operator added */
        if (orig_num && curr_num)
        {
            check_taut = 1;
            CompAtomData_GetNumMapping( composite_norm_data, orig_num, curr_num );
        }
    }
repeatj:
    remove = 0;
    at1 = unit->bkbonds[j][0];
    at2 = unit->bkbonds[j][1];
    /* djb-rwth: removing redundant code */
    for (k = 0; k < orig_at_data->at[at1 - 1].valence; k++)
    {
        if (orig_at_data->at[at1 - 1].neighbor[k] != at2 - 1)
            continue;
        /* djb-rwth: removing redundant code */
    }
    /*if ( border > 1 ) */
    /* djb-rwth: removing redundant code */
    bond_is_untouchable = 0;
    if (check_taut && composite_norm_data && composite_norm_data->at && curr_num) /* djb-rwth: fixing a NULL pointer dereference */
    {
        for (k = 0; k < composite_norm_data->at[curr_num[at1]].valence; k++)
        {
            if (composite_norm_data->at[curr_num[at1]].neighbor[k] != curr_num[at2])
            {
                continue;
            }
            btype = composite_norm_data->at[curr_num[at1]].bond_type[k];
            bond_is_untouchable = ( btype == BOND_TAUTOM ); /*|| btype == BOND_ALTERN );*/
            break;
        }
    }
    if (bond_is_untouchable)
    {
        remove = 1;
    }

    if (remove)
    {
        OAD_PolymerUnit_RemoveLinkFromCRUChain( at1, at2, &unit->nbkbonds, unit->bkbonds );
    }
    else
    {
        ++j;
    }

    if (j < unit->nbkbonds)
    {
        goto repeatj;
    }

    if (orig_num)
    {
        inchi_free( orig_num );
    }
    if (curr_num)
    {
        inchi_free( curr_num );
    }

    return;
}


/****************************************************************************
 Remove bond (at1, at2)
****************************************************************************/
void OAD_PolymerUnit_RemoveLinkFromCRUChain( int at1,
                                             int at2,
                                             int *nbonds,
                                             int **bonds )
{
    int p, q;
#if 0
    if (at1 > at2)
    {
        int tmp = at1;
        at1 = at2;
        at2 = tmp;
    }
#endif
    for (p = 0; p < *nbonds; p++)
    {
        if (bonds[p][0] == at1 && bonds[p][1] == at2)
        {
            for (q = p + 1; q < *nbonds; q++)
            {
                bonds[q - 1][0] = bonds[q][0];
                bonds[q - 1][1] = bonds[q][1];
            }
            ( *nbonds )--;
            break;
        }
    }

    return;
}


/****************************************************************************
 Debug print polymer data for a given SRU
****************************************************************************/
void OAD_PolymerUnit_DebugTrace( OAD_PolymerUnit *u )
{
    char *conn = "ABSENT", *typ = "ABSENT", *styp = "ABSENT";

    if (!u)
    {
        return;
    }

    if (u->conn == 1)
    {
        conn = "HT"; /* djb-rwth: ignoring LLVM warning: possible presence of global variables */
    }
    else if (u->conn == 2)
    {
        conn = "HH"; /* djb-rwth: ignoring LLVM warning: possible presence of global variables */
    }
    else if (u->conn == 3)
    {
        conn = "EU"; /* djb-rwth: ignoring LLVM warning: possible presence of global variables */
    }

    if (u->type == 0)
    {
        typ = "NONE"; /* djb-rwth: ignoring LLVM warning: possible presence of global variables */
    }
    else if (u->type == 1)
    {
        typ = "SRU"; /* djb-rwth: ignoring LLVM warning: possible presence of global variables */
    }
    else if (u->type == 2)
    {
        typ = "MON"; /* djb-rwth: ignoring LLVM warning: possible presence of global variables */
    }
    else if (u->type == 3)
    {
        typ = "COP"; /* djb-rwth: ignoring LLVM warning: possible presence of global variables */
    }
    else if (u->type == 4)
    {
        typ = "MOD"; /* djb-rwth: ignoring LLVM warning: possible presence of global variables */
    }
    else if (u->type == 5)
    {
        typ = "MER"; /* djb-rwth: ignoring LLVM warning: possible presence of global variables */
    }

    if (u->subtype == 1)
    {
        styp = "ALT"; /* djb-rwth: ignoring LLVM warning: possible presence of global variables */
    }
    else if (u->subtype == 2)
    {
        styp = "RAN"; /* djb-rwth: ignoring LLVM warning: possible presence of global variables */
    }
    else if (u->subtype == 3)
    {
        styp = "BLK"; /* djb-rwth: ignoring LLVM warning: possible presence of global variables */
    }

    {
        int i, k;
        int na, nb;

        ITRACE_("\n\nPOLYMER UNIT @ %-p", u);

        ITRACE_( "\n\tid=%-d   label=%-d   type=%-s   subtype=%-s   conn=%-s   subscr='%-s'\n",
                u->id, u->label, typ, styp, conn, u->smt );

        ITRACE_( "\tBracket1 coords: %-f, %-f, %-f, %-f\n", u->xbr1[0], u->xbr1[1], u->xbr1[2], u->xbr1[3] );
        ITRACE_( "\tBracket2 coords: %-f, %-f, %-f, %-f\n", u->xbr2[0], u->xbr2[1], u->xbr2[2], u->xbr2[3] );

        na = u->na;
        ITRACE_( "\t%-d atoms { ", na );
        for (k = 0; k < na - 1; k++)
        {
            ITRACE_( " %-d, ", u->alist[k] );
        }
        ITRACE_( " %-d }\n", u->alist[na - 1] );

        nb = u->nb;
        ITRACE_( "\t%-d bonds crossing unit borders { ", nb );

        for (k = 0; k < nb; k++)
        {
            ITRACE_( " %-d-%-d ", u->blist[2 * k], u->blist[2 * k + 1] );
        }
        ITRACE_( "}\n" );

        ITRACE_( "\tCRU caps and end atoms { " );

        ITRACE_( "*%-d-[-%-d(end1) ... ", u->cap1, u->end_atom1 );
        ITRACE_( "%-d(end2)-]-*%-d", u->end_atom2, u->cap2 );
        ITRACE_( " }\n" );

        ITRACE_( "\tBackbone bonds (may include 'artificially cyclizing' one) : %-d bonds ", u->nbkbonds );
        if (u->nbkbonds)
        {
            ITRACE_(" { ");
            for (i = 0; i < u->nbkbonds; i++)
            {
                ITRACE_( "(%-d, %-d)  ", u->bkbonds[i][0], u->bkbonds[i][1] );
            }
            ITRACE_("}\n");
        }
    }
    

    return;
}


/****************************************************************************
 Debug print the whole polymer data
****************************************************************************/
void OAD_Polymer_DebugTrace( OAD_Polymer *p )
{
    int i;

    if (!p)
    {
        return;
    }

    ITRACE_( "\n\n* POLYMER INFO @ %-p (%-d group(s))", p, p->n );
    ITRACE_( "\n\n* %-d star atoms: ", p->n_pzz );
    for (i = 0; i < p->n_pzz; i++)
    {
        ITRACE_( " %-d", p->pzz[i] );
    }

    for (i = 0; i < p->n; i++)
    {
        ITRACE_( "\n* Polymer unit %-d", i );
        OAD_PolymerUnit_DebugTrace( p->units[i] );
    }
    ITRACE_( "\n* Really-do-PS = %-d", p->really_do_frame_shift );
    ITRACE_( "\n* Frame_shift_scheme = %-d", p->frame_shift_scheme );
    ITRACE_( "\n* Edit-repeats = %-d", p->edit_repeats );
    ITRACE_( "\n* End POLYMER INFO\n" );
    return;

}


/****************************************************************************/
int  OAD_Polymer_GetRepresentation( OAD_Polymer *p )
{
    int i, n_source_based_units = 0, n_structure_based_units = 0;

    if (!p)
    {
        return NO_POLYMER;
    }

    for (i = 0; i < p->n; i++)
    {
        if (p->units[i]->nb == 2 || p->units[i]->nbkbonds > 0 ||
            ( ( p->units[i]->cap1 > 0 ) && ( p->units[i]->cap2 > 0 ) ))
        {
            p->units[i]->representation = POLYMER_REPRESENTATION_STRUCTURE_BASED;
            n_structure_based_units++;
        }
        else if (p->units[i]->nb == 0)
        {
            p->units[i]->representation = POLYMER_REPRESENTATION_SOURCE_BASED;
            n_source_based_units++;
        }
    }
    if (p->n == n_source_based_units)
    {
        return POLYMER_REPRESENTATION_SOURCE_BASED;
    }
    else if (p->n == n_structure_based_units)
    {
        return POLYMER_REPRESENTATION_STRUCTURE_BASED;
    }
    else if (n_source_based_units        &&
              n_structure_based_units &&
              ( n_source_based_units + n_structure_based_units ) == p->n)
    {
        /* TODO: check if SRU/MON are embedded in a single COP (is this check really necessary? ??) */
        return POLYMER_REPRESENTATION_MIXED;
    }
#if 0
    else if (p->n == ( n_source_based_units + n_structure_based_units ))
    {
        /*
            Structure based presentation may include no-crossing bond units
            which only serve as embedding for ( >1 ) structure-based SRU's.
            The code below accounts for this.
        */
        if (n_source_based_units < n_structure_based_units)
        {
            int j, atom, atom_is_shared_with_struct_based_unit = 0;
            for (i = 0; i < p->n; i++)
            {
                int k;
                if (p->units[i]->representation != POLYMER_REPRESENTATION_SOURCE_BASED)
                    continue;
                for (k = 0; k < p->units[i]->na; k++)
                {
                    atom = p->units[i]->alist[k];
                    if (is_in_the_ilist( p->pzz, atom, p->n_pzz ))
                        continue;
                    atom_is_shared_with_struct_based_unit = 0;
                    for (j = 0; j < p->n; j++)
                    {
                        if (p->units[j]->representation != POLYMER_REPRESENTATION_STRUCTURE_BASED)
                            continue;
                        if (is_in_the_ilist( p->units[j]->alist, atom, p->units[j]->na ))
                        {
                            atom_is_shared_with_struct_based_unit = 1;
                            break;
                        }
                    }
                    if (!atom_is_shared_with_struct_based_unit)
                        break;
                }
                if (!atom_is_shared_with_struct_based_unit)
                    break;
            }
            if (atom_is_shared_with_struct_based_unit)
                return POLYMER_REPRESENTATION_STRUCTURE_BASED;
        }
        return POLYMER_REPRESENTATION_MIXED;
    }
#endif

    return POLYMER_REPRESENTATION_UNRECOGNIZED;
}


/****************************************************************************
 Open pre-cyclized CRUs appropriately (i.e., make frame shift)
****************************************************************************/
void OAD_Polymer_SmartReopenCyclizedUnits( OAD_Polymer *p,
                                           inp_ATOM    *at,
                                           int         nat,
                                           int         *num_inp_bonds )
{
    int i;
    /* djb-rwth: fixing oss-fuzz issue #68329 */
    OAD_AtProps *aprops = (OAD_AtProps*)inchi_calloc((long long)nat + 1, sizeof(OAD_AtProps)); /* djb-rwth: cast operator added */
    /* nat + 1: add extra element for possibe 1-based indexing */

    if (!p)
    {
        inchi_free(aprops); /* djb-rwth: avoiding memory leak */
        return;
    }
    if (p->n < 1)
    {
        inchi_free(aprops); /* djb-rwth: avoiding memory leak */
        return;
    }
    if (!p->really_do_frame_shift)
    {
        inchi_free(aprops); /* djb-rwth: avoiding memory leak */
        return;
    }
    /* djb-rwth: fixing oss-fuzz issue #68329 */
    if (nat <= 0)
    {
        inchi_free(aprops); /* djb-rwth: avoiding memory leak */
        return;
    }

    /*ITRACE_( "\n\n*********************************************************************\n* ENTERING OAD_Polymer_SmartReopenCyclizedUnits()" );
    OAD_Polymer_DebugTrace( p );*/

    /* Set atom properties for sorting */
    if (!aprops || !at) /* djb-rwth: fixing oss-fuzz issue #68329, #68286 */
    {
        inchi_free(aprops); /* djb-rwth: avoiding memory leak */
        return;
    }
    OAD_Polymer_SetAtProps( p, at, nat, num_inp_bonds, aprops, NULL ); /* NULL as we alredy are in 1-based cano_nums while at i2s/i2i*/
    for (i = 0; i < p->n; i++)
    {
        if (p->units[i]) /* djb-rwth: fixing oss-fuzz issue #68329 */
        {
            OAD_PolymerUnit *u = p->units[i];
            if (p->frame_shift_scheme == FSS_NONE)
            {
                continue;
            }
            if ( /* !u->cyclizable || u->cyclized  || */
                u->nbkbonds < 1 ||
                u->cap1 < 1 || u->cap2 < 1 ||
                u->cap1 > nat || u->cap2 > nat)
            {
                continue;
            }
            if (OAD_PolymerUnit_SetReopeningDetails(u, at))
            {
                int senior_bond;
                OAD_PolymerUnit_SortBackboneBondsAndSetSeniors(u, at, aprops, &senior_bond);
            }
            OAD_PolymerUnit_ReopenCyclized(u, at, aprops, nat, num_inp_bonds);
        }
    }

    p->really_do_frame_shift = 0;
    inchi_free( aprops );

    return;
}


/****************************************************************************/
void OAD_PolymerUnit_ReopenCyclized( OAD_PolymerUnit *u,
                                     inp_ATOM        *at,
                                     OAD_AtProps     *aprops,
                                     int             nat,
                                     int             *num_inp_bonds )
{
    int bond_type, bond_stereo;

    if (u->cyclizable == CLOSING_SRU_RING)
    {
        /* Decyclize artificially introduced bond */
        OrigAtData_RemoveBond( u->end_atom1 - 1, u->end_atom2 - 1,
                               at, &bond_type, &bond_stereo, num_inp_bonds );
    }
    else if (u->cyclizable == CLOSING_SRU_HIGHER_ORDER_BOND)
    {
        OrigAtData_DecreaseBondOrder( u->end_atom1 - 1, u->end_atom2 - 1, at );
    }
    else if (u->cyclizable == CLOSING_SRU_DIRADICAL)
    {
        if (at[u->end_atom1 - 1].radical == RADICAL_TRIPLET)
        {
            at[u->end_atom1 - 1].radical = 0;
        }
    }

    /* Add explicitly connections to star atoms */
    OrigAtData_AddSingleStereolessBond( u->cap1 - 1, u->end_atom1 - 1,
                                        at, num_inp_bonds );
    OrigAtData_AddSingleStereolessBond( u->cap2 - 1, u->end_atom2 - 1,
                                        at, num_inp_bonds );

    /* Create crossing bonds */
    u->nb = 2;
    u->nbkbonds = 0;
    if (!u->blist)
    {
        u->blist = (int *) inchi_calloc( 2 * (long long)u->nb, sizeof( int ) ); /* djb-rwth: cast operator added */
    }
    if (!u->blist)
    {
        return;
    }

    u->blist[0] = u->cap1;
    u->blist[1] = u->end_atom1;
    u->blist[2] = u->cap2;
    u->blist[3] = u->end_atom2;

    return;
}


/****************************************************************************/
int OAD_PolymerUnit_SetReopeningDetails( OAD_PolymerUnit *u, inp_ATOM *at )
{
    int k;

    /* Check reopening  type */

    /* Caps are separated by one atom - that's not error but do nothing */
    if (u->nbkbonds == 0)
    {
        ;
    }
    else if (u->nbkbonds == 1)
    {
        u->end_atom1 = u->bkbonds[0][0];
        u->end_atom2 = u->bkbonds[0][1];

        if (u->end_atom1 == u->end_atom2)
        {
#ifdef ALLOW_CLOSING_SRU_VIA_DIRADICAL
            u->cyclizable = CLOSING_SRU_DIRADICAL;
#else
            u->cyclizable = CLOSING_SRU_NOT_APPLICABLE;
#endif
        }
        else
        {
            /* If caps are separated by two atoms - that's not error but do nothing */
            for (k = 0; k < at[u->end_atom1 - 1].valence; k++)
            {
                if (at[u->end_atom1 - 1].neighbor[k] == u->end_atom2 - 1)
                {
                    if (at[u->end_atom1 - 1].bond_type[k] > 1)
#ifdef ALLOW_CLOSING_SRU_VIA_HIGHER_ORDER_BOND
                        u->cyclizable = CLOSING_SRU_HIGHER_ORDER_BOND;
#else
/*                  u->cyclizable = CLOSING_SRU_NOT_APPLICABLE;*/
#endif
break;
                }
            }
        }
    } /*    */

    return u->nbkbonds;
}



/****************************************************************************/
void OAD_PolymerUnit_SortBackboneBondsAndSetSeniors( OAD_PolymerUnit *u,
                                                     inp_ATOM        *at,
                                                     OAD_AtProps     *aprops,
                                                     int             *senior_bond )
{
    int j, *bnum = NULL;

    *senior_bond = 0;

    /* Sort backbone (== frame shiftable) bonds if necessary */
    if (u->nbkbonds > 1)
    {
        bnum = (int *) inchi_calloc( u->nbkbonds, sizeof( int ) );
        if (bnum)
        {
            for (j = 0; j < u->nbkbonds; j++)
            {
                bnum[j] = j;
            }
            OAD_PolymerUnit_SortBackboneBonds( u, aprops, bnum );
            *senior_bond = bnum[0];
            inchi_free( bnum );
        }
    }

    /* v. 1.05+ : place senior atom the first ("left") in the senior bond */
    if (OAD_Polymer_IsFirstAtomRankLower( u->bkbonds[*senior_bond][0], u->bkbonds[*senior_bond][1], aprops ) == 1)
    {
        int tmp = u->bkbonds[*senior_bond][0];
        u->bkbonds[*senior_bond][0] = u->bkbonds[*senior_bond][1];
        u->bkbonds[*senior_bond][1] = tmp;
    }

    u->end_atom1 = u->bkbonds[*senior_bond][0];
    u->end_atom2 = u->bkbonds[*senior_bond][1];

    return;
}


/****************************************************************************/
void OAD_PolymerUnit_SortBackboneBonds( OAD_PolymerUnit *u,
                                        OAD_AtProps     *aprops,
                                        int             *bnum )
{
    int i, j, tmp;
    int n = u->nbkbonds;
    if (NULL == bnum)
    {
        return;
    }
    for (i = 1; i < n; i++)
    {
        tmp = bnum[i];
        j = i - 1;
        while (j >= 0 && OAD_Polymer_CompareBackboneBondsSeniority( u->bkbonds[bnum[j]], u->bkbonds[tmp], aprops ) > 0)
        {
            bnum[j + 1] = bnum[j];
            j--;
        }
        bnum[j + 1] = tmp;
    }

    return;
}


/****************************************************************************
 For sorting SRU cyclizing bonds (PS=='frame-shift') in descending order
 In general:  favor greater max-rank end
              if max ends are the same, favor lesser min-rank end
****************************************************************************/
int  OAD_Polymer_CompareBackboneBondsSeniority( int* b1, int* b2, OAD_AtProps *aprops )
{
    int b1min, b1max, b2min, b2max, tmp, cmp = 0;

    /* Find min and max ext-ranked ends of the both bonds */
    b1max = b1[0]; b1min = b1[1];
    b2max = b2[0]; b2min = b2[1];
    if (OAD_Polymer_IsFirstAtomRankLower( b1min, b1max, aprops ) == -1)
    {
        tmp = b1max;
        b1max = b1min;
        b1min = tmp;
    }
    if (OAD_Polymer_IsFirstAtomRankLower( b2min, b2max, aprops ) == -1)
    {
        tmp = b2max;
        b2max = b2min;
        b2min = tmp;
    }

    /* Compare bonds' seniority */

    /* First, favor the bond which has greater ext-rank end
       NB: the result may be 0, that is, equal max ext. ranks */

    cmp = OAD_Polymer_CompareRanksOfTwoAtoms( b1max, b2max, aprops );
    if (cmp == 1)
    {
        return   1;        /* rank(b1max) < rank(b2max), so bond2 is senior */
    }
    else if (cmp == -1)
    {
        return  -1;        /* rank(b1max) > rank(b2max), so bond1 is senior */
    }

    /* Max ends are of the same rank, so favor the bond with lesser min-rank end
       NB: the result may NOT be 0, that is, the case is always resolved    */

    cmp = OAD_Polymer_CompareRanksOfTwoAtoms( b1min, b2min, aprops ); /*OAD_Polymer_IsFirstAtomRankLower( b1min, b2min, aprops );*/

    if (cmp == 1)
    {
        return  -1;         /* rank(b1min) < rank(b2min), so bond1 is senior */
    }
    else if (cmp == -1)
    {
        return   1;         /* rank(b1min) > rank(b2min), so bond2 is senior */
    }

    /* Min ends are of the same rank. Here is the time to compare directly
       which canonical number is larger of max-ends ... */

    if (b1max < b2max)
    {
        return  1;
    }
    if (b1max > b2max)
    {
        return -1;
    }

    /* ... they are the same, so compare which canonical number is larger for min-ends ... */

    if (b1min < b2min)
    {
        return -1;          /* b1min < b2min, so bond1 is senior */
    }
    if (b1min > b2min)
    {
        return  1;          /* b1min > b2min, so bond2 is senior */
    }

    return 0;    /* we should not reach there */
}


/****************************************************************************
 Compare seniority of two atoms in polymer SRU
 loosely following IUPAC guidelines.
 NB: no last resort check here, so 0 (=='same seniority') may be returned
****************************************************************************/
int OAD_Polymer_CompareRanksOfTwoAtoms( int atom1, int atom2, OAD_AtProps *aprops )
{
    const int HETEROCYC = 3, HETEROAT = 2, CARBOCYC = 1, CARBOAT = 0;
        /* NB: Carbon's rank is always 2, next to the lowest */

    int a1 = atom1 - 1;
    int a2 = atom2 - 1;
    int a1typ = CARBOAT;
    int a2typ = CARBOAT;

    /* djb-rwth: fixing oss-fuzz issue #69501, #68277 */
    if ((a1 < 0) || (a2 < 0))
    {
        return 0;
    }

    if (aprops[a1].ring_size > 2)
    {
        if (aprops[a1].ring_erank <= 2)
        {
            a1typ = CARBOCYC;
        }
        else
        {
            a1typ = HETEROCYC;
        }
    }
    else
    {
        if (aprops[a1].erank == 2)
        {
            a1typ = CARBOAT;
        }
        else
        {
            a1typ = HETEROAT;
        }
    }

    if (aprops[a2].ring_size > 2)
    {
        if (aprops[a2].ring_erank <= 2)
        {
            a2typ = CARBOCYC;
        }
        else
        {
            a2typ = HETEROCYC;
        }
    }
    else
    {
        if (aprops[a2].erank == 2)
        {
            a2typ = CARBOAT;
        }
        else
        {
            a2typ = HETEROAT;
        }
    }

    /* Compare */

    /*
        Follow IUPAC Rule 1
            'The basic order of seniority of subunits is:
                heterocyclic rings and ring systems > heteroatom chains >
                    > carbocyclic rings and ring systems > acyclic carbon chains'
    */

    if (a1typ == HETEROCYC && a2typ == HETEROCYC)   /* a1 and a2 are HETEROCYC */
    {
        /* Try resolving by senior-heteroatom ring */
        if (aprops[a1].ring_erank < aprops[a2].ring_erank)
        {
            return  1;
        }
        if (aprops[a1].ring_erank > aprops[a2].ring_erank)
        {
            return -1;
        }
        /* Same senior-heteroatom rings, try resolving by total ring size */
        if (aprops[a1].ring_size < aprops[a2].ring_size)
        {
            return  1;
        }
        if (aprops[a1].ring_size > aprops[a2].ring_size)
        {
            return -1;
        }
        /* Could not resolve... */
        return 0;
    }
    else if (a1typ == HETEROCYC)
    {
        return -1;  /* a1 is HETEROCYC, a2 is any other (==junior) */
    }
    else if (a2typ == HETEROCYC)
    {
        return  1;  /* a2 is HETEROCYC, a1 is any other (==junior) */
    }

    /* HETEROCYC left out */

    if (a1typ == HETEROAT && a2typ == HETEROAT)  /* a1 and a2 are HETEROAT */
    {
        if (aprops[a1].erank < aprops[a2].erank)
        {
            return  1;
        }
        if (aprops[a1].erank > aprops[a2].erank)
        {
            return -1;
        }
        /* Could not resolve... */
        return 0;
    }
    else if (a1typ == HETEROAT)
    {
        return -1;  /* a1 is HETEROAT, a2 is any other (==junior) */
    }
    else if (a2typ == HETEROAT)
    {
        return  1;  /* a2 is HETEROAT, a1 is any other (==junior) */
    }

    /* HETEROAT left out */
    if (a1typ == CARBOCYC && a2typ == CARBOCYC) /* a1 and a2 are CARBOCYC */
    {
        /* Same senior-atom (C) ring, try resolving by total ring size */
        if (aprops[a1].ring_size < aprops[a2].ring_size)
        {
            return  1;
        }
        if (aprops[a1].ring_size > aprops[a2].ring_size)
        {
            return -1;
        }
        /* Could not resolve... */
        return 0;
    }
    else if (a1typ == CARBOCYC)
    {
        return -1;
    }
    else if (a2typ == CARBOCYC)
    {
        return  1;
    }

    return 0;        /* 0 means unresolved. It is legal here */
}


/****************************************************************************
 Compare seniority of two atoms in polymer SRU
 loosely following IUPAC guidelines.
 Always return non-zero result, that is, resolve the case.
****************************************************************************/
int OAD_Polymer_IsFirstAtomRankLower( int atom1, int atom2, OAD_AtProps *aprops )
{
    /* Compare ext-ranks */
    int result = OAD_Polymer_CompareRanksOfTwoAtoms( atom1, atom2, aprops );

    if (result)
    {
        return result;
    }

    /* Could not resolve who is junior by extended-ranks...             */
    /* As a last resort, simply check which canonical number is lesser  */
    if (atom1 < atom2)
    {
        return  1;
    }
    if (atom1 > atom2)
    {
        return -1;
    }

    /* should not reach there */
    return 0;
}


/****************************************************************************/
void OAD_ValidateAndSortOutPseudoElementAtoms( ORIG_ATOM_DATA *orig_at_data,
                                               int treat_polymers,
                                               int bNPZz,
                                               int *err,
                                               char *pStrErr )
{
    int i, k, n_pseudo = 0;
    int nsgroups = 0, nzz = 0;
    int nat = orig_at_data->num_inp_atoms;
    OAD_Polymer *pd = orig_at_data->polymer;
    OAD_PolymerUnit* u = NULL;

    int pseudos_allowed = (bNPZz == 1) || (treat_polymers != POLYMERS_NO);

    for (k = 0; k < nat; k++)
    {
        int is_zz = 0, is_star = 0, is_zy=0;


        /* Though "Zy" is present in our internal Periodic Table,
           _input_ "Zy" atoms are prohibited.
        */
        if (!strcmp( orig_at_data->at[k].elname, "Zy" ))
        {
            is_zy = 1;
#if 0
            Disabled 2020-04-07
            TREAT_ERR(*err, (70 + 5), "Invalid element(s):");
            TREAT_ERR(*err, (70 + 5), orig_at_data->at[k].elname);
            continue;
#endif 
        }
        is_star = !strcmp( orig_at_data->at[k].elname, "*" );
        if (!is_star)
        {
            is_zz = !strcmp( orig_at_data->at[k].elname, "Zz" );
        }

        if (is_star || is_zz || is_zy)
        {
            n_pseudo++;
            if (0==pseudos_allowed)
            {
                /* That's an error */
                TREAT_ERR( *err, ( 70 + 5 ), "Invalid element(s):" );
                TREAT_ERR( *err, ( 70 + 5 ), orig_at_data->at[k].elname );
                continue;
            }

            /* Now check if valid pseudoelement atom */

            /* Should be strictly univalent and single-bonded */
            /* Should not have isotopic enrichment */
            if ( orig_at_data->at[k].valence > 1          ||
                 orig_at_data->at[k].chem_bonds_valence>1 /* || orig_at_data->at[k].iso_atw_diff != 0    */
                                                           )
            {
                TREAT_ERR( *err, ( 70 + 7 ), "Invalid pseudo element(s) bonding" );
                /*TREAT_ERR( *err, ( 70 + 7 ), orig_at_data->at[k].elname );*/
                continue;
            }


            /* Now convert both "*" and "Zz" temporarily to "Zy" */
            mystrncpy( orig_at_data->at[k].elname, "Zy", sizeof( "Zy" ) );
        }
    }

    orig_at_data->n_zy = 0;
    nzz = 0;
    if (orig_at_data->valid_polymer)
    {
        nsgroups = pd->n;
        /* If applicable, check each CRU and back-convert "Zy" to "Zz" (polymer-related pseudoelement atoms)
           if they are from valid paired CRU crossing bond out-of-bracket caps */
        for (i = 0; i < nsgroups; i++)
        {
            u = pd->units[i];
            if (u)
            {
                if (u->cap1_is_undef && u->cap2_is_undef)
                {
                    /* valid pair: CRU is capped with two undefined-nature atoms, call them "Zz", finally */
                    strcpy( orig_at_data->at[u->cap1 - 1].elname, "Zz" );
                    strcpy( orig_at_data->at[u->cap2 - 1].elname, "Zz" );
                    nzz+= 2;
                }
            }
        }
        orig_at_data->polymer->n_pzz = nzz;
    }

    orig_at_data->n_zy = n_pseudo - nzz;
    if (orig_at_data->n_zy)
    {
        /* Have non-polymer-related pseudoelement atoms */
        if (0==bNPZz)
        {
            TREAT_ERR( *err, ( 70 + 4 ), "Polymer-unrelated pseudoatoms are not allowed" );
        }
    }

    if (*err)
    {
        orig_at_data->valid_polymer = 0;
    }

}


/****************************************************************************/
int Inp_Atom_GetBondType(inp_ATOM *at, int iatom1, int iatom2)
{
    int i;

    for (i = 0; i < at[iatom1].valence; i++)
    {
        if (at[iatom1].neighbor[i] == iatom2)
        {
                return at[iatom1].bond_type[i];
        }
    }

    return -1;
}
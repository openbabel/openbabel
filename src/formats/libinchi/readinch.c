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
#include <limits.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <ctype.h>

#include "mode.h"

#include "ichierr.h"
#include "extr_ct.h"
#include "ichi_io.h"

#include "inchi_api.h"
#include "readinch.h"

#include "bcf_s.h"

#define NO_ATOM          (-1) /* non-existent (central) atom */


#ifndef AB_MAX_WELL_DEFINED_PARITY
#define AB_MAX_WELL_DEFINED_PARITY inchi_max(INCHI_PARITY_ODD, INCHI_PARITY_EVEN) /* 1, 2 => well defined parities, uncluding 'unknown' */
#endif

#ifndef AB_MIN_WELL_DEFINED_PARITY
#define AB_MIN_WELL_DEFINED_PARITY inchi_min(INCHI_PARITY_ODD, INCHI_PARITY_EVEN) /* min(INCHI_PARITY_ODD, INCHI_PARITY_EVEN) */
#endif


#if ( defined(TARGET_API_LIB) || defined(TARGET_EXE_USING_API) )

#ifndef AB_PARITY_UNKN
#define AB_PARITY_UNKN   3  /* 3 => user marked as unknown parity */
#endif
#ifndef AB_PARITY_UNDF
#define AB_PARITY_UNDF   4  /* 4 => parity cannot be defined because of symmetry or not well defined geometry */
#endif

#define ATOM_PARITY_WELL_DEF(X)     (AB_MIN_WELL_DEFINED_PARITY <= (X) && (X) <= AB_MAX_WELL_DEFINED_PARITY)

#define SB_PARITY_FLAG  0x38 /* disconnected structure has undef. parity */

#define SB_PARITY_SHFT  3

#define SB_PARITY_MASK  0x07

#define SB_PARITY_1(X) (X & SB_PARITY_MASK)  /* refers to connected structure */

#define SB_PARITY_2(X) (((X) >> SB_PARITY_SHFT) & SB_PARITY_MASK) /* refers to connected structure */



#endif /* #if ( defined(TARGET_API_LIB) || defined(TARGET_EXE_USING_API) )  */


#if ( defined( TARGET_LIB_FOR_WINCHI ) || defined(TARGET_EXE_STANDALONE) )
int Extract0DParities( inp_ATOM *at,
                       int nNumAtoms,
                       inchi_Stereo0D *stereo0D,
                       int num_stereo0D,
                       char *pStrErr,
                       int *err,
                       int vABParityUnknown );
#endif


void find_and_interpret_structure_header( char *szLine,
                                          char *pSdfLabel,
                                          char *pSdfValue,
                                          unsigned long *Id,
                                          int hlen,
                                          ReadINCHI_CtlData *ir );



/****************************************************************************/
inchi_Stereo0D * CreateInchi_Stereo0D( int num_stereo0D )
{
    return (inchi_Stereo0D*) inchi_calloc( num_stereo0D, sizeof( inchi_Stereo0D ) );
}


/****************************************************************************/
void FreeInchi_Stereo0D( inchi_Stereo0D **stereo0D )
{
    if (stereo0D && *stereo0D)
    {
        inchi_free( *stereo0D );
        *stereo0D = NULL;
    }
}


/****************************************************************************/
int Extract0DParities( inp_ATOM *at,
                       int nNumAtoms,
                       inchi_Stereo0D *stereo0D,
                       int num_stereo0D,
                       char *pStrErr,
                       int *err,
                       int vABParityUnknown )
{

    /*
        vABParityUnknown holds actual value of an internal constant signifying
        unknown parity: either the same as for undefined parity (default==standard)
        or a specific one (non-std; requested by SLUUD switch).
    */
    if (stereo0D && num_stereo0D > 0)
    {
        int i0D, a2, k, k_prev, type, j, j1, j2, len, parity, parityNM;
        int sb_ord_from_i1, sb_ord_from_i2, sn_ord_from_i1, sn_ord_from_i2;
        AT_NUMB i1n, i2n, i1, i2;

        for (i0D = 0; i0D < num_stereo0D; i0D++)
        {
            parity = ( stereo0D[i0D].parity & SB_PARITY_MASK );
            parityNM = ( stereo0D[i0D].parity & SB_PARITY_FLAG ) >> SB_PARITY_SHFT;

            if (parity == INCHI_PARITY_NONE ||
                 (parity != INCHI_PARITY_ODD && parity != INCHI_PARITY_EVEN &&
                 parity != INCHI_PARITY_UNKNOWN && parity != INCHI_PARITY_UNDEFINED)) /* djb-rwth: addressing LLVM warning */
            {
                char szTemp[16];
                sprintf(szTemp, "#%d", i0D + 1);
                TREAT_ERR( *err, 0, "Wrong 0D stereo descriptor(s):" );
                TREAT_ERR( *err, 0, szTemp );
                continue; /* warning */
            }

            type = stereo0D[i0D].type;
            a2 = stereo0D[i0D].central_atom; /* central atom or -1 */
            j = -1;
            /* djb-rwth: removing redundant code */
            sb_ord_from_i1 = sb_ord_from_i2 = sn_ord_from_i1 = sn_ord_from_i2 = -1;
            i1n = i2n = i1 = i2 = MAX_ATOMS + 1;

            if (( type == INCHI_StereoType_Tetrahedral ||
                ((type == INCHI_StereoType_Allene ) &&
                  0 <= a2 && a2 < nNumAtoms)) ||
                  (type == INCHI_StereoType_DoubleBond &&
                  a2 == NO_ATOM)) /* djb-rwth: addressing LLVM warning */
            {
                /* test the quadruplet */
                for (j = 0, k_prev = -1; j < 4; j++, k_prev = k)
                {
                    k = stereo0D[i0D].neighbor[j];
                    if (k < 0 || k >= nNumAtoms || k_prev == k)
                        break;
                    /* tetrahedral atom connectivity test */
                    if (type == INCHI_StereoType_Tetrahedral &&
                         k != a2 &&
                         !is_in_the_list( at[a2].neighbor, (AT_NUMB) k, at[a2].valence ))
                    {
                        break;
                    }
                    /* Double bond, Cumulene and allene are tested in the next if() */
                }
            }

            /* Find in the adjacency lists the double bond neighbor that leads to the opposite atom */
            if (j == 4 && ( type == INCHI_StereoType_Allene ||
                type == INCHI_StereoType_DoubleBond ))
            {
                AT_NUMB *p1 = NULL, *p2 = NULL, *q1 = NULL, *q2 = NULL;
                i1n = (AT_NUMB) stereo0D[i0D].neighbor[0];
                i1 = (AT_NUMB) stereo0D[i0D].neighbor[1];
                i2 = (AT_NUMB) stereo0D[i0D].neighbor[2];
                i2n = (AT_NUMB) stereo0D[i0D].neighbor[3];

                /* find q1 and q2 */
                if (!( q1 = is_in_the_list( at[i1].neighbor, i1n, at[i1].valence ) ) ||
                     !( q2 = is_in_the_list( at[i2].neighbor, i2n, at[i2].valence ) ))
                {
                    j = -2; /* error flag */
                }
                else
                {
                    /* allene or cumulene; follow double bonds from i1 to i2 */
                    if (!( p1 = is_in_the_list( at[i1].neighbor, i2, at[i1].valence ) ))
                    {
                        /* at[i1] and at[i2] are not connected: can be only allene or cumulene */

                        AT_NUMB prev, cur, next;
                        int     num_dbond, i, next_ord, half_len;

                        cur = next = i1;
                        len = half_len = 0;
                        while (len < 20)
                        {
                            /* arbitrary very high upper limit to prevent infinite loop */
                            prev = cur;
                            cur = next;

                            for (i = 0, num_dbond = 0; i < at[cur].valence; i++)
                            {
                                /* follow double bond path && avoid going back */
                                if (at[cur].bond_type[i] == BOND_TYPE_DOUBLE &&
                                     prev != at[cur].neighbor[i])
                                {
                                    next = at[cur].neighbor[i];
                                    next_ord = i;
                                    num_dbond++;
                                }
                            }

                            if (num_dbond == 1 && next != i1)
                            {
                                len++;
                                if (len == 1)
                                    sb_ord_from_i1 = next_ord;

                                if (type == INCHI_StereoType_Allene && next == (AT_NUMB) a2)
                                    half_len = len;
                            }
                            else
                                break;
                        }

                        if (cur == i2 && prev != cur && 0 == num_dbond && len > 1 &&
                            ( p2 = is_in_the_list( at[i2].neighbor, prev, at[i2].valence ) ) &&
                            ( type != INCHI_StereoType_Allene || len == 2 * half_len ))
                        {
                            sb_ord_from_i2 = p2 - at[i2].neighbor;
                            sn_ord_from_i1 = q1 - at[i1].neighbor;
                            sn_ord_from_i2 = q2 - at[i2].neighbor;
                        }
                        else
                        {
                            j = -5; /* error flag */
                        }
                    }
                    else
                    {
                        /* allene must have been already processed, otherwise error */
                        if (type == INCHI_StereoType_Allene)
                        {
                            /* error: atoms #1 and #2 of allene are connected */
                            j = -3; /* error flag */
                        }
                        else
                        {
                            /* double bond only; the bond type is not checked because at the end
                               of the normalization it may happen to be alternating */
                            if (type == INCHI_StereoType_DoubleBond &&
                                ( p2 = is_in_the_list( at[i2].neighbor, i1, at[i2].valence ) ))
                            {
                                sb_ord_from_i1 = p1 - at[i1].neighbor;
                                sb_ord_from_i2 = p2 - at[i2].neighbor;
                                sn_ord_from_i1 = q1 - at[i1].neighbor;
                                sn_ord_from_i2 = q2 - at[i2].neighbor;
                            }
                            else
                            {
                                j = -4; /* error flag */
                            }
                        }
                    }
                }
            }

            if (j != 4)
            {
                char szTemp[16];
                sprintf(szTemp, "#%d", i0D + 1);
                TREAT_ERR( *err, 0, "Wrong 0D stereo descriptor(s):" );
                TREAT_ERR( *err, 0, szTemp );
                continue; /* error */
            }

            switch (type)
            {
                case INCHI_StereoType_None:
                    continue;
                case INCHI_StereoType_DoubleBond:
                case INCHI_StereoType_Allene:
                    for (j1 = 0; j1 < MAX_NUM_STEREO_BONDS && at[i1].sb_parity[j1]; j1++)
                    {
                        ;
                    }
                    for (j2 = 0; j2 < MAX_NUM_STEREO_BONDS && at[i2].sb_parity[j2]; j2++)
                    {
                        ;
                    }
                    if (j1 < MAX_NUM_STEREO_BONDS && j2 < MAX_NUM_STEREO_BONDS &&
                         sb_ord_from_i1 >= 0 && sb_ord_from_i2 >= 0 &&
                         sn_ord_from_i1 >= 0 && sn_ord_from_i2 >= 0)
                    {
                        switch (parity)
                        {
                            case INCHI_PARITY_ODD:
                                at[i1].sb_parity[j1] = AB_PARITY_ODD;
                                at[i2].sb_parity[j2] = AB_PARITY_EVEN;
                                break;
                            case INCHI_PARITY_EVEN:
                                at[i1].sb_parity[j1] = AB_PARITY_ODD;
                                at[i2].sb_parity[j2] = AB_PARITY_ODD;
                                break;
                            case INCHI_PARITY_UNDEFINED:
                                at[i1].sb_parity[j1] = AB_PARITY_UNDF;
                                at[i2].sb_parity[j2] = AB_PARITY_UNDF;
                                break;
                            default:
                                if (parity == INCHI_PARITY_UNKNOWN)
                                {
                                    at[i1].sb_parity[j1] = vABParityUnknown;
                                    at[i2].sb_parity[j2] = vABParityUnknown;
                                }
                                else
                                {
                                    at[i1].sb_parity[j1] = AB_PARITY_NONE;
                                    at[i2].sb_parity[j2] = AB_PARITY_NONE;
                                }
                                break;
                        }
                        switch (parityNM)
                        {
                            case INCHI_PARITY_ODD:
                                at[i1].sb_parity[j1] |= AB_PARITY_ODD << SB_PARITY_SHFT;
                                at[i2].sb_parity[j2] |= AB_PARITY_EVEN << SB_PARITY_SHFT;
                                break;
                            case INCHI_PARITY_EVEN:
                                at[i1].sb_parity[j1] |= AB_PARITY_ODD << SB_PARITY_SHFT;
                                at[i2].sb_parity[j2] |= AB_PARITY_ODD << SB_PARITY_SHFT;
                                break;
                            case INCHI_PARITY_UNDEFINED:
                                at[i1].sb_parity[j1] |= AB_PARITY_UNDF << SB_PARITY_SHFT;
                                at[i2].sb_parity[j2] |= AB_PARITY_UNDF << SB_PARITY_SHFT;
                                break;
                            default:
                                if (parityNM == INCHI_PARITY_UNKNOWN)
                                {
                                    at[i1].sb_parity[j1] |= vABParityUnknown << SB_PARITY_SHFT;
                                    at[i2].sb_parity[j2] |= vABParityUnknown << SB_PARITY_SHFT;
                                }
                                break;
                        }
                        at[i1].sb_ord[j1] = sb_ord_from_i1;
                        at[i1].sn_ord[j1] = sn_ord_from_i1;
                        at[i1].sn_orig_at_num[j1] = at[i1n].orig_at_number;

                        at[i2].sb_ord[j2] = sb_ord_from_i2;
                        at[i2].sn_ord[j2] = sn_ord_from_i2;
                        at[i2].sn_orig_at_num[j2] = at[i2n].orig_at_number;
                    }
                    break;
                case INCHI_StereoType_Tetrahedral:
                    switch (parity)
                    {
                        case INCHI_PARITY_ODD:
                            at[a2].p_parity = AB_PARITY_ODD;
                            break;
                        case INCHI_PARITY_EVEN:
                            at[a2].p_parity = AB_PARITY_EVEN;
                            break;
                        case INCHI_PARITY_UNDEFINED:
                            at[a2].p_parity = AB_PARITY_UNDF;
                            break;
                        default:
                            if (parity == INCHI_PARITY_UNKNOWN)
                            {
                                at[a2].p_parity = vABParityUnknown;
                                break;
                            }
                            else
                            {
                                continue;
                            }
                    }
                    for (j = 0; j < 4; j++)
                    {
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
        /* Note: to this stage, unk/undef set to what was requested              */
        /*( through vABParityUnknown )  (2009-12-12)                             */
        FixUnkn0DStereoBonds( at, nNumAtoms );
#ifdef TARGET_API_LIB

        if ((k = ReconcileAllCmlBondParities( at, nNumAtoms, 0 ))) /* djb-rwth: addressing LLVM warning */
        {
            char szErrCode[16];
            sprintf( szErrCode, "%d", k );
            AddErrorMessage( pStrErr, "0D Parities Reconciliation failed:" );
            AddErrorMessage( pStrErr, szErrCode );
        }

#endif
    }

    return 0;
}


/****************************************************************************/
char* FindToken( INCHI_IOSTREAM *inp_file,
                 int *bTooLongLine,
                 const char *sToken,
                 int lToken,
                 char *szLine,
                 int nLenLine,
                 char *p,
                 int *res )
{
    char *q;
    int   res2;

    while (!( q = strstr( p, sToken ) ))
    {
        if (( q = strrchr( p, '/' ) ) && ( q + lToken > szLine + *res ))
        {
            *res -= q - szLine; /* res = the length of the szLine to be left in */
            memmove(szLine, q, (long long)*res + 1); /* djb-rwth: cast operator added */
        }
        else
        {
            *res = 0;
        }

        res2 = inchi_ios_getsTab1( szLine + *res, nLenLine - *res - 1,
                                   inp_file, bTooLongLine );

        if (!*bTooLongLine || 0 > res2)
        {
            /* the line is over or end of file */
            return NULL;
        }
        else
        {
            *res += res2;
            p = szLine;
        }
    }

    return q + lToken;
}


/****************************************************************************/
char *LoadLine( INCHI_IOSTREAM *inp_file,
                int *bTooLongLine,
                int *bItemIsOver,
                char **s,
                char *szLine,
                int nLenLine,
                int nMinLen2Load,
                char *p,
                int *res )
{

    int pos = p - szLine, res2;

    if (!*bItemIsOver && nLenLine - ( *res - pos ) > nMinLen2Load)
    {
        /* load the next portion if possible */

        if (pos)
        {
            *res -= pos;
            memmove(szLine, p, (long long)*res + 1); /* djb-rwth: cast operator added */
            p = szLine;
            if (*s)
            {
                *s -= pos;
            }

            /* djb-rwth: removing redundant code */
        }

        res2 = inchi_ios_getsTab1( szLine + *res,
                                   nLenLine - *res - 1,
                                   inp_file, bTooLongLine );

        if (res2 > 0)
        {
            *bItemIsOver = ( ( *s = strchr( p + *res, '/' ) ) || !*bTooLongLine );
            *res += res2;
        }
        else
        {
            *bItemIsOver = 1;
        }
    }

    return p;
}


/*****************************************************************************/
#define AT_BONDS_VAL(AT,I)  AT[I].chem_bonds_valence
#define ISOLATED_ATOM       15
#define NUM_ISO_Hk(AT,I,K)  AT[I].num_iso_H[K]
#define inchi_NUMH2(AT,N)   NUMH(AT,N)
#define AT_NUM_BONDS(AT)    (AT).valence
#define IS_METAL_ATOM(AT,I) is_el_a_metal( AT[I].el_number )


/****************************************************************************/
char szLine_i2i[INCHI_LINE_LEN]; /* djb-rwth: placed as a global variable to avoid function buffer issues */
int InchiToInpAtom( INCHI_IOSTREAM *inp_file,
                     MOL_COORD **szCoord,
                     int bDoNotAddH,
                     int vABParityUnknown,
                     INPUT_TYPE nInputType,
                     inp_ATOM **at,
                     int max_num_at,
                     int *num_dimensions,
                     int *num_bonds,
                     char *pSdfLabel,
                     char *pSdfValue,
                     unsigned long *Id,
                     INCHI_MODE *pInpAtomFlags,
                     int *err,
                     char *pStrErr )
{
    int      num_atoms = 0, bItemIsOver; /* djb-rwth: removing redundant variables */
    int      i, k, k2, res, bond_type, bond_stereo1, bond_stereo2, bond_char, neigh, bond_parity, bond_parityNM;
    /* djb-rwth: removing redundant variables */
    char     *p, *q, *s, parity;
    int      b2D = 0, b3D = 0, b23D, nNumBonds = 0, bNonZeroXYZ, bNonMetal;
    int      len_stereo0D = 0, max_len_stereo0D = 0;
    inp_ATOM    *atom = NULL;
    MOL_COORD    *pszCoord = NULL;
    INCHI_MODE    InpAtomFlags = 0; /* 0 or FLAG_INP_AT_NONCHIRAL or FLAG_INP_AT_CHIRAL */
    inchi_Stereo0D   *atom_stereo0D = NULL;
    static const char szIsoH[] = "hdt";
    /* plain tags */
    static const char sStructHdrPln[] = "Structure:";
    static char sStructHdrPlnAuxStart[64] = ""; /*"$1.1Beta/";*/
    static int  lenStructHdrPlnAuxStart = 0;
    static const char sStructHdrPlnRevAt[] = "/rA:";
    static const char sStructHdrPlnRevBn[] = "/rB:";
    static const char sStructHdrPlnRevXYZ[] = "/rC:";
    const  char *sToken;
    int  lToken, len, hlen;
    
    ReadINCHI_CtlData ir;

    if (!lenStructHdrPlnAuxStart)
    {
        lenStructHdrPlnAuxStart = sprintf(sStructHdrPlnAuxStart, "AuxInfo=");
    }

    if (at)
    {
        if (*at && max_num_at)
            memset(*at, 0, max_num_at * sizeof(**at)); /* djb-rwth: memset_s C11/Annex K variant? */
        if (szCoord && *szCoord)
        {
            inchi_free(*szCoord);
            *szCoord = NULL;
        }
    }
    /* djb-rwth: removing redundant code */

    ir.bHeaderRead = ir.bErrorMsg = ir.bRestoreInfo = 0;
    *num_dimensions = *num_bonds = 0;


    if (nInputType != INPUT_INCHI_PLAIN)
    {
        return num_atoms;
    }
            

    /*
        Extract reversibility info from plain text INChI format
    */

    ir.bHeaderRead = 0; /* djb-rwth: removing redundant code */
    while (0 < (res = inchi_ios_getsTab(szLine_i2i, sizeof(szLine_i2i) - 1, inp_file, &ir.bTooLongLine)))
    {

        if (!ir.bTooLongLine &&
            (hlen = sizeof(sStructHdrPln) - 1, !memcmp(szLine_i2i, sStructHdrPln, hlen)))

        {
            num_atoms = 0;
            find_and_interpret_structure_header(szLine_i2i, pSdfLabel, pSdfValue,
                Id, hlen, &ir);
        }

        else if (!memcmp(szLine_i2i, sStructHdrPlnAuxStart, lenStructHdrPlnAuxStart))
        {
            /* Reject to deal with polymers for now */
            if (strstr(szLine_i2i, "/Z:"))
            {
                *err = INCHI_INP_ERROR_ERR;
                num_atoms = INCHI_INP_ERROR_RET;
                TREAT_ERR(*err, 0, "Reading polymer AuxInfo is not supported yet");
                goto bypass_end_of_INChI_plain;
            }

            /* Found the header of the AuxInfo, read AuxInfo head of the line */
            if (!ir.bHeaderRead)
            {
                ir.ulongID = 0LU;
                if (Id)
                {
                    *Id = ir.ulongID;
                }
                if (pSdfLabel)
                {
                    pSdfLabel[0] = '\0';
                }
                if (pSdfValue)
                {
                    pSdfValue[0] = '\0';
                }
            }

            ir.bHeaderRead = 0;

            /* Check for empty "AuxInfo=ver//" */
            p = strchr(szLine_i2i + lenStructHdrPlnAuxStart, '/');

            if (p && p[1] == '/' && (!p[2] || '\n' == p[2]))
            {
                goto bypass_end_of_INChI_plain;
            }

            /*
                Search for atoms block (plain)
            */

            p = szLine_i2i;
            sToken = sStructHdrPlnRevAt;
            lToken = sizeof(sStructHdrPlnRevAt) - 1;

            /* Search for sToken in the line; load next segments of the line if sToken has not found */

            p = FindToken(inp_file, &ir.bTooLongLine, sToken, lToken,
                szLine_i2i, sizeof(szLine_i2i), p, &res);

            if (!p)
            {
                *err = INCHI_INP_ERROR_ERR;
                num_atoms = INCHI_INP_ERROR_RET;
                TREAT_ERR(*err, 0, "Missing atom data");
                goto bypass_end_of_INChI_plain;
            }
            else
            {
                /* atoms block started */
                i = 0;
                /* djb-rwth: removing redundant code */
                bItemIsOver = (s = strchr(p, '/')) || !ir.bTooLongLine;
                while (1)
                {

                    p = LoadLine(inp_file, &ir.bTooLongLine, &bItemIsOver, &s,
                        szLine_i2i, sizeof(szLine_i2i), INCHI_LINE_ADD, p, &res);

                    if (!i)
                    {
                        /* allocate atom */
                        num_atoms = strtol(p, &q, 10);

                        if (!num_atoms || !q || !*q)
                        {
                            num_atoms = 0; /* no atom data */
                            goto bypass_end_of_INChI_plain;
                        }
                        p = q;

                        /* Molfile chirality flag */
                        switch (*p)
                        {
                        case 'c':
                            InpAtomFlags |= FLAG_INP_AT_CHIRAL;
                            p++;
                            break;
                        case 'n':
                            InpAtomFlags |= FLAG_INP_AT_NONCHIRAL;
                            p++;
                            break;
                        }

                        if (at && *at)
                        {
                            if (num_atoms > max_num_at)
                            {
                                inchi_free(*at);
                                *at = NULL;
                            }
                            else
                            {
                                memset(*at, 0, max_num_at * sizeof(**at)); /* djb-rwth: memset_s C11/Annex K variant? */
                                atom = *at;
                            }
                        }

                        if (!at || !*at)
                        {

                            atom = CreateInpAtom(num_atoms + 1);

                            if (!atom)
                            {
                                num_atoms = INCHI_INP_FATAL_RET; /* was -1; error */
                                *err = INCHI_INP_FATAL_ERR;
                                TREAT_ERR(*err, 0, "Out of RAM");
                                goto bypass_end_of_INChI_plain;
                            }
                        }

                        {
                            max_len_stereo0D = num_atoms + 1;

                            atom_stereo0D = CreateInchi_Stereo0D(max_len_stereo0D);

                            if (!atom_stereo0D)
                            {
                                num_atoms = INCHI_INP_FATAL_RET; /* fatal error: cannot allocate */
                                *err = INCHI_INP_FATAL_ERR;
                                TREAT_ERR(*err, 0, "Out of RAM");
                                goto bypass_end_of_INChI_plain;
                            }
                        }
                    }

                    /* element, first char */
                    if (!isalpha(UCINT * p) || !isupper(UCINT * p) || i >= num_atoms)
                    {
                        break; /* end of atoms block */
                    }

                    atom[i].elname[0] = *p++;

                    /* element, second char */
                    if (isalpha(UCINT * p) && islower(UCINT * p))
                    {
                        atom[i].elname[1] = *p++;
                    }

                    atom[i].el_number = get_periodic_table_number(atom[i].elname);

                    /* bonds' valence + number of non-isotopic H */
                    if (isdigit(UCINT * p))
                    {
                        AT_BONDS_VAL(atom, i) = (char)strtol(p, &q, 10);
                        if (!AT_BONDS_VAL(atom, i))
                            AT_BONDS_VAL(atom, i) = ISOLATED_ATOM; /* same convention as in MOLfile, found zero bonds valence */
                        p = q;
                    }

                    /* charge */
                    atom[i].charge = (*p == '+') ? 1 : (*p == '-') ? -1 : 0;
                    if (atom[i].charge)
                    {
                        p++;
                        if (isdigit(UCINT * p))
                        {
                            atom[i].charge *= (S_CHAR)(strtol(p, &q, 10) & CHAR_MASK);
                            p = q;
                        }
                    }

                    /* radical */
                    if (*p == '.')
                    {
                        p++;
                        if (isdigit(UCINT * p))
                        {
                            atom[i].radical = (S_CHAR)strtol(p, &q, 10);
                            p = q;
                        }
                    }

                    /* isotopic mass */
                    if (*p == 'i')
                    {
                        p++;
                        if (isdigit(UCINT * p))
                        {
                            int mw = strtol(p, &q, 10);
                            p = q;
                            mw -= get_atomic_mass_from_elnum(atom[i].el_number);
                            if (mw >= 0)
                                mw++;
                            atom[i].iso_atw_diff = mw;
                        }
                    }

                    /* parity */
                    switch (*p)
                    {
                    case 'o':
                        parity = INCHI_PARITY_ODD;
                        p++;
                        break;
                    case 'e':
                        parity = INCHI_PARITY_EVEN;
                        p++;
                        break;
                    case 'u':
                        parity = INCHI_PARITY_UNKNOWN;
                        p++;
                        break;
                    case '?':
                        parity = INCHI_PARITY_UNDEFINED;
                        p++;
                        break;
                    default:
                        parity = 0;
                        break;
                    }

                    if (parity)
                    {
                        atom_stereo0D[len_stereo0D].central_atom = i;
                        atom_stereo0D[len_stereo0D].parity = parity;
                        atom_stereo0D[len_stereo0D].type = INCHI_StereoType_Tetrahedral;
                        len_stereo0D++;
                    }

                    /* isotopic h, d, t */
                    for (k = 0; k < NUM_H_ISOTOPES; k++)
                    {
                        if (*p == szIsoH[k])
                        {
                            NUM_ISO_Hk(atom, i, k) = 1;
                            p++;
                            if (isdigit(UCINT * p))
                            {
                                NUM_ISO_Hk(atom, i, k) = (char)strtol(p, &q, 10);
                                p = q;
                            }
                        }
                    }

                    i++;
                }

                if (!bItemIsOver || i != num_atoms || (s && p != s)) /* djb-rwth: addressing LLVM warning */
                {
                    num_atoms = INCHI_INP_ERROR_RET; /* error */
                    *err = INCHI_INP_ERROR_ERR;
                    TREAT_ERR(*err, 0, "Wrong number of atoms");
                    goto bypass_end_of_INChI_plain;
                }
            }

            /*
                Search for bonds block (plain) and read it
            */

            /*p = szLine;*/
            sToken = sStructHdrPlnRevBn;
            lToken = sizeof(sStructHdrPlnRevBn) - 1;

            /* Search for sToken in the line; load next segments of the line if sToken has not found */
            p = FindToken(inp_file, &ir.bTooLongLine, sToken, lToken, szLine_i2i, sizeof(szLine_i2i), p, &res);

            if (!p)
            {
                num_atoms = INCHI_INP_ERROR_RET; /* error */
                *err = INCHI_INP_ERROR_ERR;
                TREAT_ERR(*err, 0, "Missing bonds data");
                goto bypass_end_of_INChI_plain;
            }
            else
            {
                /* bonds block started */

                i = 1;

                /* djb-rwth: removing redundant code */

                bItemIsOver = (s = strchr(p, '/')) || !ir.bTooLongLine;

                if (1 == num_atoms)
                {
                    /* needed because the next '/' may be still out of szLine */

                    p = LoadLine(inp_file, &ir.bTooLongLine, &bItemIsOver, &s,
                        szLine_i2i, sizeof(szLine_i2i), INCHI_LINE_ADD, p, &res);
                }

                while (i < num_atoms)
                {

                    p = LoadLine(inp_file, &ir.bTooLongLine, &bItemIsOver, &s,
                        szLine_i2i, sizeof(szLine_i2i), INCHI_LINE_ADD, p, &res);

                    if (i >= num_atoms || (s && p >= s)) /* djb-rwth: addressing LLVM warning */
                    {
                        break; /* end of bonds (plain) */
                    }

                    /* bond, first char */
                    if (*p == ';')
                    {
                        p++;
                        i++;
                        continue;
                    }

                    if (!isalpha(UCINT * p))
                    {
                        num_atoms = INCHI_INP_ERROR_RET; /* error */
                        *err = INCHI_INP_ERROR_ERR;
                        TREAT_ERR(*err, 0, "Wrong bonds data");
                        goto bypass_end_of_INChI_plain;
                    }

                    bond_char = *p++;

                    /* bond parity */
                    switch (*p)
                    {
                    case '-':
                        bond_parity = INCHI_PARITY_ODD;
                        p++;
                        break;
                    case '+':
                        bond_parity = INCHI_PARITY_EVEN;
                        p++;
                        break;
                    case 'u':
                        bond_parity = INCHI_PARITY_UNKNOWN;
                        p++;
                        break;
                    case '?':
                        bond_parity = INCHI_PARITY_UNDEFINED;
                        p++;
                        break;
                    default:
                        bond_parity = 0;
                        break;
                    }

                    if (bond_parity)
                    {
                        switch (*p)
                        {
                        case '-':
                            bond_parityNM = INCHI_PARITY_ODD;
                            p++;
                            break;
                        case '+':
                            bond_parityNM = INCHI_PARITY_EVEN;
                            p++;
                            break;
                        case 'u':
                            bond_parityNM = INCHI_PARITY_UNKNOWN;
                            p++;
                            break;
                        case '?':
                            bond_parityNM = INCHI_PARITY_UNDEFINED;
                            p++;
                            break;
                        default:
                            bond_parityNM = 0;
                            break;
                        }
                    }
                    else
                    {
                        bond_parityNM = 0;
                    }

                    /* neighbor of the current atom */
                    if (!isdigit(UCINT * p))
                    {
                        num_atoms = INCHI_INP_ERROR_RET; /* error */
                        *err = INCHI_INP_ERROR_ERR;
                        TREAT_ERR(*err, 0, "Wrong bonds data");
                        goto bypass_end_of_INChI_plain;
                    }

                    neigh = (int)strtol(p, &q, 10) - 1;

#if ( FIX_CURE53_ISSUE_HEAP_BUFFER_OVERFLOW_INCHITOINPATOM==1 )
                    if (i >= num_atoms || neigh >= num_atoms || neigh < 0)
                    {
#else
                    if (i >= num_atoms || neigh >= num_atoms) {
#endif
                        num_atoms = INCHI_INP_ERROR_RET; /* error */
                        *err = INCHI_INP_ERROR_ERR;
                        TREAT_ERR(*err, 0, "Bond to nonexistent atom");
                        goto bypass_end_of_INChI_plain;
                }

                    p = q;
                    bond_stereo1 = bond_stereo2 = 0;

                    /* bond type & 2D stereo */
                    switch (bond_char)
                    {
                    case 'v':
                        bond_type = INCHI_BOND_TYPE_SINGLE;
                        bond_stereo1 = INCHI_BOND_STEREO_SINGLE_1EITHER;
                        bond_stereo2 = INCHI_BOND_STEREO_SINGLE_2EITHER;
                        break;
                    case 'V':
                        bond_type = INCHI_BOND_TYPE_SINGLE;
                        bond_stereo1 = INCHI_BOND_STEREO_SINGLE_2EITHER;
                        bond_stereo2 = INCHI_BOND_STEREO_SINGLE_1EITHER;
                        break;
                    case 'w':
                        bond_type = INCHI_BOND_TYPE_DOUBLE;
                        bond_stereo1 =
                            bond_stereo2 = INCHI_BOND_STEREO_DOUBLE_EITHER;
                        break;
                    case 's':
                        bond_type = INCHI_BOND_TYPE_SINGLE;
                        break;
                    case 'd':
                        bond_type = INCHI_BOND_TYPE_DOUBLE;
                        break;
                    case 't':
                        bond_type = INCHI_BOND_TYPE_TRIPLE;
                        break;
                    case 'a':
                        bond_type = INCHI_BOND_TYPE_ALTERN;
                        break;
                    case 'p':
                        bond_type = INCHI_BOND_TYPE_SINGLE;
                        bond_stereo1 = INCHI_BOND_STEREO_SINGLE_1UP;
                        bond_stereo2 = INCHI_BOND_STEREO_SINGLE_2UP;
                        break;
                    case 'P':
                        bond_type = INCHI_BOND_TYPE_SINGLE;
                        bond_stereo1 = INCHI_BOND_STEREO_SINGLE_2UP;
                        bond_stereo2 = INCHI_BOND_STEREO_SINGLE_1UP;
                        break;
                    case 'n':
                        bond_type = INCHI_BOND_TYPE_SINGLE;
                        bond_stereo1 = INCHI_BOND_STEREO_SINGLE_1DOWN;
                        bond_stereo2 = INCHI_BOND_STEREO_SINGLE_2DOWN;
                        break;
                    case 'N':
                        bond_type = INCHI_BOND_TYPE_SINGLE;
                        bond_stereo1 = INCHI_BOND_STEREO_SINGLE_2DOWN;
                        bond_stereo2 = INCHI_BOND_STEREO_SINGLE_1DOWN;
                        break;
                    default:
                        num_atoms = INCHI_INP_ERROR_RET; /* error */
                        *err = INCHI_INP_ERROR_ERR;
                        TREAT_ERR(*err, 0, "Wrong bond type");
                        goto bypass_end_of_INChI_plain;
                    }

                    k = AT_NUM_BONDS(atom[i])++; /* AT_NUM_BONDS(AT)  ==>  (AT).valence */

                    atom[i].bond_type[k] = bond_type;
                    atom[i].bond_stereo[k] = bond_stereo1;
                    atom[i].neighbor[k] = (AT_NUMB)neigh;

                    k2 = AT_NUM_BONDS(atom[neigh])++; /* AT_NUM_BONDS(AT)  ==>  (AT).valence */
                    atom[neigh].bond_type[k2] = bond_type;
                    atom[neigh].bond_stereo[k2] = bond_stereo2;
                    atom[neigh].neighbor[k2] = (AT_NUMB)i;

                    bond_parity |= (bond_parityNM << SB_PARITY_SHFT);

                    if (bond_parity)
                    {
                        if (max_len_stereo0D <= len_stereo0D)
                        {
                            /* realloc atom_Stereo0D */

                            inchi_Stereo0D* new_atom_stereo0D = CreateInchi_Stereo0D(max_len_stereo0D + num_atoms);

                            if (!new_atom_stereo0D)
                            {
                                num_atoms = INCHI_INP_FATAL_RET; /* fatal error: cannot allocate */
                                *err = INCHI_INP_FATAL_ERR;
                                TREAT_ERR(*err, 0, "Out of RAM");
                                goto bypass_end_of_INChI_plain;
                            }

                            memcpy(new_atom_stereo0D, atom_stereo0D, len_stereo0D * sizeof(*atom_stereo0D));
                            FreeInchi_Stereo0D(&atom_stereo0D);
                            atom_stereo0D = new_atom_stereo0D;
                            max_len_stereo0D += num_atoms;
                        }

                        /* (a) i may be allene endpoint and     neigh = allene middle point or
                            (b) i may be allene middle point and neigh = allene endpoint
                            !!!!! CURRENTLY ONLY (b) IS ALLOWED !!!!!
                        */

                        atom_stereo0D[len_stereo0D].neighbor[1] = neigh; /* neigh < i */
                        atom_stereo0D[len_stereo0D].neighbor[2] = i;
                        atom_stereo0D[len_stereo0D].parity = bond_parity;
                        atom_stereo0D[len_stereo0D].type = INCHI_StereoType_DoubleBond; /* incl allenes & cumulenes */
                        len_stereo0D++;
                    }
            }

                if (!bItemIsOver || i != num_atoms || (s && p != s)) /* djb-rwth: addressing LLVM warning */
                {
                    num_atoms = INCHI_INP_ERROR_RET; /* error */
                    *err = INCHI_INP_ERROR_ERR;
                    TREAT_ERR(*err, 0, "Wrong number of bonds");
                    goto bypass_end_of_INChI_plain;
                }
        }

            /*
                Search for coordinates block (plain)
            */

            /*p = szLine;*/
            sToken = sStructHdrPlnRevXYZ;
            lToken = sizeof(sStructHdrPlnRevXYZ) - 1;

            /* search for sToken in the line; load next segments of the line if sToken has not found */
            p = FindToken(inp_file, &ir.bTooLongLine, sToken, lToken, szLine_i2i, sizeof(szLine_i2i), p, &res);

            if (!p)
            {
                num_atoms = INCHI_INP_ERROR_RET; /* error */
                *err = INCHI_INP_ERROR_ERR;
                TREAT_ERR(*err, 0, "Missing atom coordinates data");
                goto bypass_end_of_INChI_plain;
            }
            else
            {
                /* Coordinates block started */
                if ((pszCoord = (MOL_COORD*)inchi_malloc(inchi_max(num_atoms, 1) * sizeof(MOL_COORD)))) /* djb-rwth: addressing LLVM warning */
                {
                    memset(pszCoord, ' ', inchi_max(num_atoms, 1) * sizeof(MOL_COORD)); /* djb-rwth: memset_s C11/Annex K variant? */
                }
                else
                {
                    num_atoms = INCHI_INP_FATAL_RET; /* allocation error */
                    *err = INCHI_INP_FATAL_ERR;
                    TREAT_ERR(*err, 0, "Out of RAM");
                    goto bypass_end_of_INChI_plain;
                }

                i = 0;
                /* djb-rwth: removing redundant code */
                bItemIsOver = (s = strchr(p, '/')) || !ir.bTooLongLine;

                while (i < num_atoms)
                {

                    p = LoadLine(inp_file, &ir.bTooLongLine, &bItemIsOver, &s,
                        szLine_i2i, sizeof(szLine_i2i), INCHI_LINE_ADD, p, &res);

                    if (i >= num_atoms || (s && p >= s)) /* djb-rwth: addressing LLVM warning */
                    {
                        break; /* end of bonds (plain) */
                    }

                    /* coord, first char */
                    if (*p == ';')
                    {
                        for (k = 0; k < NUM_COORD; k++)
                        {
                            pszCoord[i][LEN_COORD * k + 4] = '0';
                        }
                        p++;
                        i++;
                        continue;
                    }

                    for (k = 0; k < 3; k++)
                    {
                        double xyz;
                        bNonZeroXYZ = 0;
                        if (*p == ';')
                        {
                            pszCoord[i][LEN_COORD * k + 4] = '0';
                            xyz = 0.0;
                        }
                        else
                        {
                            if (*p == ',')
                            {
                                /* empty */
                                pszCoord[i][LEN_COORD * k + 4] = '0';
                                xyz = 0.0;
                                p++;
                            }
                            else
                            {
                                xyz = strtod(p, &q);
                                bNonZeroXYZ = fabs(xyz) > MIN_BOND_LENGTH;
                                if (q != NULL)
                                {
                                    memcpy(pszCoord[i] + LEN_COORD * (long long)k, p, q - p); /* djb-rwth: cast operator added */
                                    if (*q == ',')
                                        q++;
                                    p = q;
                                }
                                else
                                    pszCoord[i][LEN_COORD * k + 4] = '0';
                            }
                        }

                        switch (k)
                        {
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

                    if (*p == ';')
                    {
                        p++; /* end of this triple of coordinates */
                        i++;
                    }
                    else
                    {
                        num_atoms = INCHI_INP_ERROR_RET; /* error in input data: atoms, bonds & coord must be present together */
                        *err = INCHI_INP_ERROR_ERR;
                        TREAT_ERR(*err, 0, "Wrong atom coordinates data");
                        goto bypass_end_of_INChI_plain;
                    }
                }

                if (!bItemIsOver || (s && p != s) || i != num_atoms) /* djb-rwth: addressing LLVM warning */
                {
                    num_atoms = INCHI_INP_ERROR_RET; /* error */
                    *err = INCHI_INP_ERROR_ERR;
                    TREAT_ERR(*err, 0, "Wrong number of coordinates");
                    goto bypass_end_of_INChI_plain;
                }
            } /* end of coordinates */

            /*
                Set special valences and implicit H (xml)
            */

            b23D = b2D | b3D;
            /* djb-rwth: removing redundant code */
            if (at)
            {
                if (!*at)
                {
                    int a1, a2, n1, n2, valence;
                    int chem_bonds_valence;
                    int    nX = 0, nY = 0, nZ = 0, nXYZ;
                    *at = atom;

                    /* special valences */

                    for (bNonMetal = 0; bNonMetal < 1; bNonMetal++)
                    {
                        for (a1 = 0; a1 < num_atoms; a1++)
                        {
                            int num_bond_type[MAX_INPUT_BOND_TYPE - MIN_INPUT_BOND_TYPE + 1];
                            int bHasMetalNeighbor = 0;

                            memset(num_bond_type, 0, sizeof(num_bond_type)); /* djb-rwth: memset_s C11/Annex K variant? */

                            valence = AT_BONDS_VAL(atom, a1); /*  save atom valence if available */
                            AT_BONDS_VAL(atom, a1) = 0;

                            atom[a1].orig_at_number = a1 + 1;

                            nX = nY = nZ = 0;

                            for (n1 = 0; n1 < AT_NUM_BONDS(atom[a1]); n1++) /*AT_NUM_BONDS(AT)  ==>  (AT).valence */
                            {
                                bond_type = atom[a1].bond_type[n1] - MIN_INPUT_BOND_TYPE;
                                if (bond_type < 0 || bond_type > MAX_INPUT_BOND_TYPE - MIN_INPUT_BOND_TYPE)
                                {
                                    bond_type = 0;
                                    TREAT_ERR(*err, 0, "Unknown bond type in InChI aux assigned as a single bond");
                                }

                                num_bond_type[bond_type] ++;
                                nNumBonds++;
                                if (b23D)
                                {
                                    neigh = atom[a1].neighbor[n1];
                                    nX |= (fabs(atom[a1].x - atom[neigh].x) > MIN_BOND_LENGTH);
                                    nY |= (fabs(atom[a1].y - atom[neigh].y) > MIN_BOND_LENGTH);
                                    nZ |= (fabs(atom[a1].z - atom[neigh].z) > MIN_BOND_LENGTH);
                                }
                            }

                            chem_bonds_valence = 0;
                            for (n1 = 0; MIN_INPUT_BOND_TYPE + n1 <= 3 && MIN_INPUT_BOND_TYPE + n1 <= MAX_INPUT_BOND_TYPE; n1++)
                            {
                                chem_bonds_valence += (MIN_INPUT_BOND_TYPE + n1) * num_bond_type[n1];
                            }

                            if (MIN_INPUT_BOND_TYPE <= INCHI_BOND_TYPE_ALTERN && INCHI_BOND_TYPE_ALTERN <= MAX_INPUT_BOND_TYPE &&
                                (n2 = num_bond_type[INCHI_BOND_TYPE_ALTERN - MIN_INPUT_BOND_TYPE]))
                            {
                                /* accept input aromatic bonds for now */
                                switch (n2)
                                {
                                case 2:
                                    chem_bonds_valence += 3;  /* =A- */
                                    break;

                                case 3:
                                    chem_bonds_valence += 4;  /* =A< */
                                    break;

                                default:
                                    /*  if 1 or >= 4 aromatic bonds then replace such bonds with single bonds */
                                    for (n1 = 0; n1 < AT_NUM_BONDS(atom[a1]); n1++) /* AT_NUM_BONDS(AT)  ==>  (AT).valence */
                                    {
                                        if (atom[a1].bond_type[n1] == INCHI_BOND_TYPE_ALTERN)
                                        {
                                            AT_NUMB* p1;
                                            a2 = atom[a1].neighbor[n1];
                                            p1 = is_in_the_list(atom[a2].neighbor, (AT_NUMB)a1, AT_NUM_BONDS(atom[a2])); /*AT_NUM_BONDS(AT)  ==>  (AT).valence*/
                                            if (p1)
                                            {
                                                atom[a1].bond_type[n1] =
                                                    atom[a2].bond_type[p1 - atom[a2].neighbor] = INCHI_BOND_TYPE_SINGLE;
                                            }
                                            else
                                            {
                                                *err = -2;  /*  Program error */
                                                TREAT_ERR(*err, 0, "Program error interpreting InChI aux");
                                                num_atoms = INCHI_INP_FATAL_RET;
                                                goto bypass_end_of_INChI_plain; /*  no structure */
                                            }
                                        }
                                    }

                                    chem_bonds_valence += n2;
                                    *err |= 32; /*  Unrecognized aromatic bond(s) replaced with single */
                                    TREAT_ERR(*err, 0, "Atom has 1 or more than 3 aromatic bonds");
                                    break;
                                }
                            }

                            /* added 2006-07-19 to process aromatic bonds same way as from molfile */
                            if (n2 && !valence)
                            {
                                int num_H = NUMH(atom, a1); /* only isotopic */
                                int chem_valence = chem_bonds_valence;
                                int bUnusualValenceArom =
                                    detect_unusual_el_valence((int)atom[a1].el_number, atom[a1].charge,
                                        atom[a1].radical, chem_valence,
                                        num_H, atom[a1].valence);
                                int bUnusualValenceNoArom =
                                    detect_unusual_el_valence((int)atom[a1].el_number, atom[a1].charge,
                                        atom[a1].radical, chem_valence - 1,
                                        num_H, atom[a1].valence);

#if ( CHECK_AROMBOND2ALT == 1 )
                                if (bUnusualValenceArom && !bUnusualValenceNoArom && 0 == nBondsValToMetal(atom, a1))
#else
                                if (bUnusualValenceArom && !bUnusualValenceNoArom)
#endif

                                {
                                    /* typically NH in 5-member aromatic ring */
                                    chem_bonds_valence--;
                                }
                            }
                            else if (n2 && valence)
                            {
                                /* atom has aromatic bonds AND the chemical valence is known */
                                int num_H = NUMH(atom, a1);
                                int chem_valence = chem_bonds_valence + num_H;
                                if (valence == chem_valence - 1)
                                {
                                    /* typically NH in 5-member aromatic ring */
                                    chem_bonds_valence--;
                                }
                            }

                            atom[a1].chem_bonds_valence = chem_bonds_valence;

                            atom[a1].num_H = get_num_H(atom[a1].elname,
                                atom[a1].num_H,
                                atom[a1].num_iso_H,
                                atom[a1].charge,
                                atom[a1].radical,
                                atom[a1].chem_bonds_valence,
                                valence,
                                0,
                                bDoNotAddH,
                                bHasMetalNeighbor);
                        }
                    }

                    nNumBonds /= 2;

                    if (b23D && nNumBonds)
                    {
                        nXYZ = nX + nY + nZ;
                        b2D = (nXYZ > 0);
                        b3D = (nXYZ == 3);
                        *num_dimensions = b3D ? 3 : b2D ? 2 : 0;
                        *num_bonds = nNumBonds;
                    }

                    /*======= 0D parities =================================*/

                    for (i = 0; i < len_stereo0D; i++)
                    {
                        AT_NUMB* p1, * p2;
                        int     sb_ord_from_a1 = -1, sb_ord_from_a2 = -1, bEnd1 = 0, bEnd2 = 0;

                        switch (atom_stereo0D[i].type)
                        {

                        case INCHI_StereoType_Tetrahedral:
                            a1 = atom_stereo0D[i].central_atom;
                            if (atom_stereo0D[i].parity && (AT_NUM_BONDS(atom[a1]) == 3 || AT_NUM_BONDS(atom[a1]) == 4))
                            {
                                int ii, kk = 0;
                                if (AT_NUM_BONDS(atom[a1]) == 3)
                                    atom_stereo0D[i].neighbor[kk++] = a1;
                                for (ii = 0; ii < AT_NUM_BONDS(atom[a1]); ii++)
                                    atom_stereo0D[i].neighbor[kk++] = atom[a1].neighbor[ii];
                            }

                            break;

                        case INCHI_StereoType_DoubleBond:
#define MAX_CHAIN_LEN 20
                            a1 = atom_stereo0D[i].neighbor[1];
                            a2 = atom_stereo0D[i].neighbor[2];
                            p1 = is_in_the_list(atom[a1].neighbor, (AT_NUMB)a2, AT_NUM_BONDS(atom[a1]));
                            p2 = is_in_the_list(atom[a2].neighbor, (AT_NUMB)a1, AT_NUM_BONDS(atom[a2]));
                            if (!p1 || !p2)
                            {
                                atom_stereo0D[i].type = INCHI_StereoType_None;
                                atom_stereo0D[i].central_atom = NO_ATOM;
                                atom_stereo0D[i].neighbor[0] =
                                    atom_stereo0D[i].neighbor[3] = -1;
                                *err |= 64; /* Error in cumulene stereo */
                                TREAT_ERR(*err, 0, "0D stereobond not recognized");
                                break;
                            }

                            /* streobond, allene, or cumulene */

                            sb_ord_from_a1 = p1 - atom[a1].neighbor;
                            sb_ord_from_a2 = p2 - atom[a2].neighbor;

                            if (AT_NUM_BONDS(atom[a1]) == 2 &&
                                atom[a1].bond_type[0] + atom[a1].bond_type[1] == 2 * INCHI_BOND_TYPE_DOUBLE &&
                                0 == inchi_NUMH2(atom, a1) &&
                                (AT_NUM_BONDS(atom[a2]) != 2 ||
                                    atom[a2].bond_type[0] + atom[a2].bond_type[1] != 2 * INCHI_BOND_TYPE_DOUBLE))
                            {
                                bEnd2 = 1; /* a2 is the end-atom, a1 is middle atom */
                            }

                            if (AT_NUM_BONDS(atom[a2]) == 2 &&
                                atom[a2].bond_type[0] + atom[a2].bond_type[1] == 2 * INCHI_BOND_TYPE_DOUBLE &&
                                0 == inchi_NUMH2(atom, a2) &&
                                (AT_NUM_BONDS(atom[a1]) != 2 ||
                                    atom[a1].bond_type[0] + atom[a1].bond_type[1] != 2 * INCHI_BOND_TYPE_DOUBLE))
                            {
                                bEnd1 = 1; /* a1 is the end-atom, a2 is middle atom */
                            }

                            if (bEnd2 + bEnd1 == 1)
                            {
                                /* allene or cumulene */

                                AT_NUMB  chain[MAX_CHAIN_LEN + 1], prev, cur, next;

                                if (bEnd2 && !bEnd1)
                                {
                                    cur = a1;
                                    a1 = a2;
                                    a2 = cur;
                                    sb_ord_from_a1 = sb_ord_from_a2;
                                }

                                sb_ord_from_a2 = -1;
                                cur = a1;
                                next = a2;
                                len = 0;
                                chain[len++] = cur;
                                chain[len++] = next;

                                while (len < MAX_CHAIN_LEN)
                                {
                                    /* arbitrary very high upper limit to prevent infinite loop */

                                    prev = cur;
                                    cur = next;
                                    /* follow double bond path && avoid going back */
                                    if (AT_NUM_BONDS(atom[cur]) == 2 &&
                                        atom[cur].bond_type[0] + atom[cur].bond_type[1] == 2 * INCHI_BOND_TYPE_DOUBLE &&
                                        0 == inchi_NUMH2(atom, cur))
                                    {
                                        next = atom[cur].neighbor[atom[cur].neighbor[0] == prev];
                                        chain[len++] = next;
                                    }
                                    else
                                    {
                                        break;
                                    }
                                }
                                if (len > 2 &&
                                    (p2 = is_in_the_list(atom[cur].neighbor, (AT_NUMB)prev, AT_NUM_BONDS(atom[cur]))))
                                {
                                    sb_ord_from_a2 = p2 - atom[cur].neighbor;
                                    a2 = cur;
                                    /* by design we need to pick up the first non-stereo-bond-neighbor as "sn"-atom */
                                    atom_stereo0D[i].neighbor[0] = atom[a1].neighbor[sb_ord_from_a1 == 0];
                                    atom_stereo0D[i].neighbor[1] = a1;
                                    atom_stereo0D[i].neighbor[2] = a2;
                                    atom_stereo0D[i].neighbor[3] = atom[a2].neighbor[sb_ord_from_a2 == 0];

                                    if (len % 2)
                                    {
                                        atom_stereo0D[i].central_atom = chain[len / 2];
                                        atom_stereo0D[i].type = INCHI_StereoType_Allene;
                                    }
                                    else
                                    {
                                        atom_stereo0D[i].central_atom = NO_ATOM;
                                    }
                                }
                                else
                                {
                                    /* error */
                                    atom_stereo0D[i].type = INCHI_StereoType_None;
                                    atom_stereo0D[i].central_atom = NO_ATOM;
                                    atom_stereo0D[i].neighbor[0] =
                                        atom_stereo0D[i].neighbor[3] = -1;
                                    *err |= 64; /* Error in cumulene stereo */
                                    TREAT_ERR(*err, 0, "Cumulene stereo not recognized (0D)");
                                }
#undef MAX_CHAIN_LEN
                            }
                            else
                            {
                                /****** a normal possibly stereogenic bond -- not an allene or cumulene *******/
                                /* by design we need to pick up the first non-stereo-bond-neighbor as "sn"-atom */
                                sb_ord_from_a1 = p1 - atom[a1].neighbor;
                                sb_ord_from_a2 = p2 - atom[a2].neighbor;
                                atom_stereo0D[i].neighbor[0] = atom[a1].neighbor[p1 == atom[a1].neighbor];
                                atom_stereo0D[i].neighbor[3] = atom[a2].neighbor[p2 == atom[a2].neighbor];
                                atom_stereo0D[i].central_atom = NO_ATOM;
                            }

                            if (atom_stereo0D[i].type != INCHI_StereoType_None &&
                                sb_ord_from_a1 >= 0 && sb_ord_from_a2 >= 0 &&
                                ATOM_PARITY_WELL_DEF(SB_PARITY_2(atom_stereo0D[i].parity)))
                            {
                                /* Detected well-defined disconnected stereo
                                    * locate first non-metal neighbors */

                                int    a, j, /* k,*/ sb_ord, cur_neigh, min_neigh; /* djb-rwth: removing redundant variables */

                                for (k = 0; k < 2; k++)
                                {
                                    a = k ? atom_stereo0D[i].neighbor[2] : atom_stereo0D[i].neighbor[1];
                                    sb_ord = k ? sb_ord_from_a2 : sb_ord_from_a1;
                                    min_neigh = num_atoms;
                                    for (j = 0; j < AT_NUM_BONDS(atom[a]); j++) /* djb-rwth: removing redundant code */
                                    {
                                        cur_neigh = atom[a].neighbor[j];
                                        if (j != sb_ord && !IS_METAL_ATOM(atom, cur_neigh))
                                        {
                                            min_neigh = inchi_min(cur_neigh, min_neigh);
                                        }
                                    }
                                    if (min_neigh < num_atoms)
                                    {
                                        atom_stereo0D[i].neighbor[k ? 3 : 0] = min_neigh;
                                    }
                                    else
                                    {
                                        TREAT_ERR(*err, 0, "Cannot find non-metal stereobond neighor (0D)");
                                    }
                                }
                            }

                            break;
                        }
                    }
                    /* end of 0D parities extraction */
/*exit_cycle:;*/
                }

                /* Transfer atom_stereo0D[] to atom[] */
                if (len_stereo0D)
                {
                    Extract0DParities(atom, num_atoms, atom_stereo0D, len_stereo0D,
                        pStrErr, err, vABParityUnknown);
                }

                if (pInpAtomFlags)
                {
                    /* save chirality flag */
                    *pInpAtomFlags |= InpAtomFlags;
                }
            }
            else if (atom)
            {
                inchi_free(atom);
                atom = NULL;
            }

#if ( FIX_READ_AUX_MEM_LEAK == 1 )
            /* 2005-08-04 avoid memory leak */
            if (atom_stereo0D) /* && !(stereo0D && *stereo0D == atom_stereo0D) ) */
            {
                FreeInchi_Stereo0D(&atom_stereo0D);
            }
#endif

            if (szCoord)
            {
                *szCoord = pszCoord;
                pszCoord = NULL;
            }
            else if (pszCoord)
            {
                inchi_free(pszCoord);
                pszCoord = NULL;
            }

            goto bypass_end_of_INChI_plain;
            /*return num_atoms;*/
    }
}    /* while ( 0 < (res = inchi_ios_getsTab( szLine, sizeof(szLine)-1, inp_file, &ir.bTooLongLine ) ) )  */

/* End of structure reading cycle */
    if (atom_stereo0D)
        FreeInchi_Stereo0D(&atom_stereo0D);
    if (res <= 0)
    {
        if (*err == INCHI_INP_ERROR_ERR)
        {
            return num_atoms;
        }
        *err = INCHI_INP_EOF_ERR;

        return INCHI_INP_EOF_RET; /* no more data */
    }

bypass_end_of_INChI_plain:
    /* Cleanup */
    if (num_atoms == INCHI_INP_ERROR_RET || num_atoms == INCHI_INP_FATAL_RET)
    {
        if (atom_stereo0D) /* djb-rwth: avoiding memory leak */
        {
            FreeInchi_Stereo0D(&atom_stereo0D);
        }

        if (atom) /* djb-rwth: fixing coverity ID #499615 */
        {
            inchi_free(atom);
        }

        if (pszCoord) /* djb-rwth: fixing coverity ID #499571 */
        {
            inchi_free(pszCoord);
        }
    }

    while (ir.bTooLongLine &&
        0 < inchi_ios_getsTab1(szLine_i2i, sizeof(szLine_i2i) - 1, inp_file, &ir.bTooLongLine))
    {
        ;
    }

    return num_atoms;

#undef AT_NUM_BONDS
#undef AT_NUMB
#undef is_in_the_list
#undef inchi_NUMH2

#undef MoreParms
#undef INPUT_FILE
#undef CreateInpAtom
#undef AT_BONDS_VAL
#undef ISOLATED_ATOM
#undef NUM_ISO_Hk
#undef IS_METAL_ATOM
}



/****************************************************************************/
void find_and_interpret_structure_header( char *szLine,
                                          char *pSdfLabel,
                                          char *pSdfValue,
                                          unsigned long *Id,
                                          int hlen,
                                          ReadINCHI_CtlData *ir )
{
    int len;
    char *p, *q;
    static const char sStructHdrPlnNoLblVal[] = " is missing";


    p = szLine + hlen;
    ir->ulongID = 0LU;

    /* structure number */
    ir->ulongID = strtoul( p, &q, 10 );
    if (q && q[0] == '.' && q[1] == ' ')
    {
        p = q + 2;
    }
    p = p + strspn( p, " \n\r" );

    if (pSdfLabel)
    {
        pSdfLabel[0] = '\0';
    }
    if (pSdfValue)
    {
        pSdfValue[0] = '\0';
    }

    if (*p)
    {
        /* has label name */

        /*p ++;*/
        if ((q = strchr( p, '=' ))) /* djb-rwth: addressing LLVM warning */
        {

            /* '=' separates label name from the value */
            len = inchi_min( q - p + 1, MAX_SDF_HEADER - 1 );

            if (pSdfLabel)
            {
                mystrncpy( pSdfLabel, p, len );
                lrtrim( pSdfLabel, &len );
            }

            p = q + 1;
            q = p + (int) strlen( p );

            if (q - p > 0)
            {
                len = inchi_min( q - p + 1, MAX_SDF_VALUE - 1 );
                if (pSdfValue)
                {
                    mystrncpy( pSdfValue, p, len );
                }
                /* djb-rwth: removing redundant code */
            }
        }
        else if ((q = strstr( p, sStructHdrPlnNoLblVal ))) /* djb-rwth: addressing LLVM warning */
        {
            len = inchi_min( q - p + 1, MAX_SDF_HEADER - 1 );
            if (pSdfLabel)
            {
                mystrncpy( pSdfLabel, p, len );
            }
            /* djb-rwth: removing redundant code */
        }
    }

    if (Id)
    {
        *Id = ir->ulongID;
    }

    ir->bHeaderRead = 1;
    ir->bErrorMsg = ir->bRestoreInfo = 0;

    return;
}

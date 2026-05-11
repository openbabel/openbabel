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

#pragma warning( disable : 4706 4127 4514 4100 4786 4996 4244 4267 )

#include <string.h>

#include "mode.h"
#include "ichierr.h"

#include "bcf_s.h"

static int already_have_this_message( char *prev_messages, const char *new_message );


/****************************************************************************/
const char *ErrMsg( int nErrorCode )
{
    const char *p;
    static char szErrMsg[64];
    switch (nErrorCode)
    {
        case 0:                      p = "";                      break;
        case CT_OVERFLOW:            p = "ARRAY OVERFLOW";        break;
        case CT_LEN_MISMATCH:        p = "LENGTH_MISMATCH";       break;
        case CT_OUT_OF_RAM:          p = "Out of RAM";            break;
        case CT_RANKING_ERR:         p = "RANKING_ERR";           break;
        case CT_ISOCOUNT_ERR:        p = "ISOCOUNT_ERR";          break;
        case CT_TAUCOUNT_ERR:        p = "TAUCOUNT_ERR";          break;
        case CT_ISOTAUCOUNT_ERR:     p = "ISOTAUCOUNT_ERR";       break;
        case CT_MAPCOUNT_ERR:        p = "MAPCOUNT_ERR";          break;
        case CT_TIMEOUT_ERR:         p = "Time limit exceeded";   break;
        case CT_ISO_H_ERR:           p = "ISO_H_ERR";             break;
        case CT_STEREOCOUNT_ERR:     p = "STEREOCOUNT_ERR";       break;
        case CT_ATOMCOUNT_ERR:       p = "ATOMCOUNT_ERR";         break;
        case CT_STEREOBOND_ERROR:    p = "STEREOBOND_ERR";        break;
        case CT_USER_QUIT_ERR:       p = "User requested termination"; break;
        case CT_REMOVE_STEREO_ERR:   p = "REMOVE_STEREO_ERR";     break;
        case CT_CALC_STEREO_ERR:     p = "CALC_STEREO_ERR";       break;
        case CT_STEREO_CANON_ERR:    p = "STEREO_CANON_ERR";      break;
        case CT_CANON_ERR:           p = "CANON_ERR";             break;
        case CT_WRONG_FORMULA:       p = "Wrong or missing chemical formula";  break;
        /*case CT_CANON_ERR2:          p = "CT_CANON_ERR2";         break;*/
        case CT_UNKNOWN_ERR:         p = "UNKNOWN_ERR";           break;
        case BNS_RADICAL_ERR:        p = "Cannot process free radical center"; break;
        case BNS_ALTBOND_ERR:        p = "Cannot process aromatic bonds";      break;
        /* v. 1.05 */
        case BNS_TIMEOUT:             p = "Structure normalization timeout";      break;

        default:
            if (nErrorCode > CT_UNKNOWN_ERR)
            {
                sprintf(szErrMsg, "No description(%d)", nErrorCode);
                p = szErrMsg;
            }
            else
            {
                sprintf(szErrMsg, "UNKNOWN_ERR(%d)", CT_UNKNOWN_ERR - nErrorCode);
                p = szErrMsg;
            }
            break;
    }

    return p;
}


/****************************************************************************/
int AddErrorMessage( char *all_messages, const char *new_message )
{
    int len_all, len;

    if (!all_messages)
    {
        return 0;
    }
    if (!new_message)
    {
        return 0;
    }
    if (!new_message[0])
    {
        return 0;
    }
    if (already_have_this_message( all_messages, new_message ))
    {
        return 1;
    }

    len_all = (int) strlen( all_messages );
    len = (int) strlen( new_message );

    if (len_all + len + 2 * ( len_all > 0 ) < STR_ERR_LEN)
    {
        /* enough room... add message and return */
        if (len_all > 0)
        {
            if (all_messages[len_all - 1] != ':')
            {
                strcat(all_messages, ";");
            }
            strcat(all_messages, " ");
        }
        strcat(all_messages, new_message);
        return 1;
    }

    /*  not enough room... add no-room marker if not yet added */
    if (strstr( all_messages, "..." ))
    {
        return 0;
    }
    if (len_all + 3 < STR_ERR_LEN)
    {
        strcat(all_messages, "...");
    }

    return 0;
}


/****************************************************************************/
int already_have_this_message( char *prev_messages, const char *new_message )
{
    int have = 0;

    char *p = strstr( prev_messages, new_message );

    if (p)
    {
        have = ( p == prev_messages || (*( p - 1 ) == ' ' && ( *( p - 2 ) == ';' || *( p - 2 ) == ':' )) ); /* djb-rwth: addressing LLVM warning */
        if (have)
        {
            int len_prev = (int) strlen( prev_messages );
            int len = (int) strlen( new_message );
            have = ( p + len == prev_messages + len_prev || (p[len] == ';' && p[len + 1] == ' ') || (p[len - 1] == ':' && p[len] == ' ') ); /* djb-rwth: addressing LLVM warning */
        }
    }

    return have;
}

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
    InChIKey: calculation of hash for InChI string

    Uses truncated SHA-256 function.
    SHA-256 implementation: Copyright (C) Brainspark B.V.,
                            see files sha2.c, sha2.h.
*/

#ifdef _MSC_VER
#if _MSC_VER > 1000
#pragma warning( disable : 4996 )
#endif
#endif

#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "sha2.h"
#include "ikey_base26.h"

#include "mode.h"
#include "inchi_api.h"
#include "util.h"

#include "bcf_s.h"

/*    Local options */

#define INCHIKEY_DEBUG 0 /* 2 */

#define INCHIKEY_FLAG_OK 0
#define INCHIKEY_NOT_VALID_FLAG 1

enum { MINOUTLENGTH = 256 };



/*    Local functions */

void fprint_digest( FILE* fw, const char *header, unsigned char *a );


/*
    EXPORTED FUNCTIONS
*/


/****************************************************************************/
EXPIMP_TEMPLATE INCHI_API
int INCHI_DECL GetStdINCHIKeyFromStdINCHI( const char* szINCHISource,
                                           char* szINCHIKey )
{
    if (strlen( szINCHISource ) < LEN_INCHI_STRING_PREFIX + 3)
    {
        return INCHIKEY_INVALID_STD_INCHI;
    }
    if (szINCHISource[LEN_INCHI_STRING_PREFIX + 1] != 'S')
    {
        return INCHIKEY_INVALID_STD_INCHI;
    }

    return GetINCHIKeyFromINCHI( szINCHISource,
                                 0, 0,
                                 szINCHIKey,
                                 (char *) NULL, (char *) NULL );
}


/****************************************************************************/
EXPIMP_TEMPLATE INCHI_API
int INCHI_DECL GetINCHIKeyFromINCHI( const char* szINCHISource,
                                     const int xtra1,
                                     const int xtra2,
                                     char* szINCHIKey,
                                     char* szXtra1,
                                     char* szXtra2 )
{
    int ret = INCHIKEY_OK;
    int ret1 = INCHIKEY_OK;
    int cn;
    size_t slen, i, j, jproto = 0, ncp, pos_slash1 = 0;
    char *str = NULL, *smajor = NULL, *sminor = NULL,
        *sproto = NULL, *stmp = NULL, tmp[MINOUTLENGTH];
    unsigned char
        digest_major[32], digest_minor[32];

    char flagstd = 'S', /* standard key */
        flagnonstd = 'N', /* non-standard key */
        flagexptl = 'B', /* experimental ('beta') key */
        flagver = 'A', /* InChI v. 1 */
        flagproto = 'N'; /* no [de]protonization , by default */
    int  nprotons;
    /*
    Protonization encoding:
    N 0
    O +1 P +2 Q +3 R +4 S +5 T +6 U +7 V +8 W +9 X +10 Y +11 Z +12
    M -1 L-2 K -3 J -4 I -5 H -6 G -7 F -8 E -9 D -10 C -11 B -12
    A < -12 or > +12
    */
    static const char *pplus = "OPQRSTUVWXYZ";
    static const char *pminus = "MLKJIHGFEDCB";
    int is_stdinchi = 0;    /* 0 -  non-standard,
                               1    standard
                              -1    experimental ('beta') */



    if (NULL != szXtra1) /* Software version 1.06 added check to fix bug with NULL szXtra, thanks to WDI */
    {
        szXtra1[0] = '\0';
    }
    if (NULL != szXtra2)
    {
        szXtra2[0] = '\0';
    }

    /* Check if input is a valid InChI string */

    /* .. non-empty */
    if (szINCHISource == NULL)
    {
        return INCHIKEY_EMPTY_INPUT;
    }

    slen = strlen( szINCHISource );

    /* .. has valid prefix */
    if (slen < LEN_INCHI_STRING_PREFIX + 3)
    {
        return INCHIKEY_INVALID_INCHI_PREFIX;
    }
    if (memcmp( szINCHISource, INCHI_STRING_PREFIX, LEN_INCHI_STRING_PREFIX ))
    {
        return INCHIKEY_INVALID_INCHI_PREFIX;
    }

    /* .. has InChI version 1 */
    /* if (!isdigit(szINCHISource[LEN_INCHI_STRING_PREFIX]) )  */
    if (szINCHISource[LEN_INCHI_STRING_PREFIX] != '1')
        return INCHIKEY_INVALID_INCHI_PREFIX;

    /* .. optionally has a 'standard' flag character */
    pos_slash1 = LEN_INCHI_STRING_PREFIX + 1;
    if (szINCHISource[pos_slash1] == 'S')
    {
        /* Standard InChI ==> standard InChIKey */
        is_stdinchi = 1;
        pos_slash1++;
    }
    else if (szINCHISource[pos_slash1] == 'B')
    {
        /* v. 1.05 Experimental ('beta') InChI ==> corresponding InChIKey */
        is_stdinchi = -1;
        pos_slash1++;
    }

    /* .. has trailing slash in the right place */
    if (szINCHISource[pos_slash1] != '/')
    {
        return INCHIKEY_INVALID_INCHI_PREFIX;
    }

    /* .. the rest of source string contains at least one a..Z0.9 or slash */
    if (!isalnum( szINCHISource[pos_slash1 + 1] ))
    {
        int valid_empty_inchi = szINCHISource[pos_slash1 + 1] == '/';
        int allowed_and_valid_err_inchi = 0;
        if (!valid_empty_inchi)
        {
            allowed_and_valid_err_inchi = szINCHISource[pos_slash1 + 1] == '?';
            if (!allowed_and_valid_err_inchi)
            {
                return INCHIKEY_INVALID_INCHI;
            }
        }
    }

    /* Ok. Will use a local copy of the source. */

    extract_inchi_substring( &str, szINCHISource, slen );
    if (NULL == str)
    {
        ret = INCHIKEY_NOT_ENOUGH_MEMORY;
        goto fin;
    }
    slen = strlen( str );


    /* Make buffers. */
    smajor = (char*) inchi_calloc( slen + 1, sizeof( char ) );
    if (NULL == smajor)
    {
        ret = INCHIKEY_NOT_ENOUGH_MEMORY; goto fin;
    }
    sminor = (char*) inchi_calloc( 2 * slen + 2, sizeof( char ) ); /* we may double the length ... */
    if (NULL == sminor)
    {
        ret = INCHIKEY_NOT_ENOUGH_MEMORY; goto fin;
    }
    stmp = (char*) inchi_calloc( slen + 1, sizeof( char ) );
    if (NULL == stmp)
    {
        ret = INCHIKEY_NOT_ENOUGH_MEMORY; goto fin;
    }
    sproto = (char*) inchi_calloc( slen + 1, sizeof( char ) );
    if (NULL == sproto)
    {
        ret = INCHIKEY_NOT_ENOUGH_MEMORY; goto fin;
    }

    szINCHIKey[0] = '\0';

    /* Extract the major block */
    smajor[0] = '\0';
    for (j = pos_slash1 + 1; j < slen - 1; j++)
    {
        if (str[j] == '/')
        {
            cn = str[j + 1];
            switch (cn)
            {
                /* anything allowed from a major part */
                case 'c':
                case 'h':
                case 'q':   continue;

                /* "/p"; protons now go to to special string, not to minor hash */
                case 'p':
                    jproto = j;
                    continue;

        /* "/f",  "/r" : may not occur in stdInChI */
                case 'f':
                case 'r':
                    if (is_stdinchi==1)
                    {
                        ret = INCHIKEY_INVALID_STD_INCHI;
                        goto fin;
                    }
                    break;

        /* anything allowed from a minor part */
                default:
                    break;
            }
            break;
        }
    }
    j++;
    if (j == slen)
    {
        j++;
    }
    else
    {
        j--;
    }

    if (jproto)
    {
        ncp = jproto - pos_slash1 - 1;
    }
    else
    {
        ncp = j - pos_slash1 - 1;
    }


    /* Trim 'InChI=1[S]/' */
    memcpy(smajor, str + pos_slash1 + 1, ncp * sizeof(str[0]));
    smajor[ncp] = '\0';


    /* Treat protonization */
    if (jproto)
    {
        /* 2009-01-07 fix bug/typo: assigned incorrect length to the protonation segment of
           source string ( was sproto[ncp]='\0'; should be sproto[lenproto]='\0'; )  */
        int lenproto = j - (int) jproto;
        if (lenproto < 3)
        {
            /* empty "/p", should not occur */
            ret = INCHIKEY_INVALID_INCHI;
            goto fin;
        }

        memcpy(sproto, str + pos_slash1 + ncp + 1, lenproto * sizeof(str[0]));   
        sproto[lenproto] = '\0';

        nprotons = strtol( sproto + 2, NULL, 10 );

        if (nprotons > 0)
        {
            if (nprotons > 12)
            {
                flagproto = 'A';
            }
            else
            {
                flagproto = pplus[nprotons - 1];
            }
        }
        else if (nprotons < 0)
        {
            if (nprotons < -12)
            {
                flagproto = 'A';
            }
            else
            {
                flagproto = pminus[-nprotons - 1];
            }
        }
        else
        {
            /* should never occur */
            ret = INCHIKEY_INVALID_STD_INCHI;
            goto fin;
        }
    }

    /* Extract the minor block. */

    if (j != slen + 1)    /* check that something exists at right.*/
    {
        ncp = slen - j;
        memcpy(sminor, str + j, (ncp) * sizeof(str[0]));
        sminor[ncp] = '\0';
    }
    else
    {
        sminor[0] = '\0';
    }


#if INCHIKEY_DEBUG
    ITRACE_( "Source:  {%-s}\n", str );
    ITRACE_( "SMajor:  {%-s}\n", smajor );
    ITRACE_( "SMinor:  {%-s}\n", sminor );
    ITRACE_( "SProto:  {%-s}\n", sproto );
#endif

    /* Compute and compose the InChIKey string. */

    /* Major hash sub-string. */
    for (i = 0; i < 32; i++)
    {
        digest_major[i] = 0;
    }

    sha2_csum( (unsigned char *) smajor, (int) strlen( smajor ), digest_major );

    sprintf(tmp, "%-.3s%-.3s%-.3s%-.3s%-.2s",
        base26_triplet_1(digest_major), base26_triplet_2(digest_major),
        base26_triplet_3(digest_major), base26_triplet_4(digest_major),
        base26_dublet_for_bits_56_to_64(digest_major));
    strcat(szINCHIKey, tmp);
#if (INCHIKEY_DEBUG>1)
    fprint_digest( stderr, "Major hash, full SHA-256", digest_major );
#endif


    /* Minor hash sub-string. */
    for (i = 0; i < 32; i++)
    {
        digest_minor[i] = 0;
    }
    slen = strlen( sminor );
    if (( slen > 0 ) && ( slen < 255 ))
    {
        strcpy(stmp, sminor);
        strcpy(sminor + slen, stmp);
    }

    sha2_csum( (unsigned char *) sminor, (int) strlen( sminor ), digest_minor );

#if (INCHIKEY_DEBUG>1)
    fprint_digest( stderr, "Minor hash, full SHA-256", digest_minor );
#endif

    strcat(szINCHIKey, "-");
    sprintf(tmp, "%-.3s%-.3s%-.2s",
        base26_triplet_1(digest_minor),
        base26_triplet_2(digest_minor),
        base26_dublet_for_bits_28_to_36(digest_minor));
    strcat(szINCHIKey, tmp);
    /* Append a standard/non-standard flag */
    slen = strlen( szINCHIKey );
    if (is_stdinchi == 1)
    {
        szINCHIKey[slen] = flagstd;
    }
    else if (is_stdinchi == -1)
    {
        szINCHIKey[slen] = flagexptl;
    }
    else
    {
        szINCHIKey[slen] = flagnonstd;
    }

    /* Append InChI v.1 flag */
    szINCHIKey[slen + 1] = flagver;

    /* Append dash  */
    szINCHIKey[slen + 2] = '-';

    /* Append protonization flag */
    szINCHIKey[slen + 3] = flagproto;
    szINCHIKey[slen + 4] = '\0';


#if INCHIKEY_DEBUG
    ITRACE_( "szINCHIKey:  {%-s}\n", szINCHIKey );
#endif

    /* Hash extensions */
    if (xtra1 && szXtra1)
    {
        get_xtra_hash_major_hex( digest_major, szXtra1 );
#if INCHIKEY_DEBUG
        fprintf( stderr, "XHash1=%-s\n", szXtra1 );
        fprintf( stderr, "j=%-d\n", j );
#endif
    }
    if (xtra2 && szXtra2)
    {
        get_xtra_hash_minor_hex( digest_minor, szXtra2 );
#if INCHIKEY_DEBUG
        fprintf( stderr, "XHash2=%-s\n", szXtra2 );
        fprintf( stderr, "j=%-d\n", j );
#endif
    }


fin:
    /* djb-rwth: fixing GH issue #87 (thanks to Ricardo Rodriguez) and oss-fuzz issue #66746 */
    if (NULL != str)
    {
        inchi_free( str );
    }
    if (NULL != smajor)
    {
        inchi_free( smajor );
    }
    if (NULL != sminor)
    {
        inchi_free( sminor );
        sminor = NULL;
    }
    if (NULL != stmp)
    {
        inchi_free( stmp );
    }
    if (NULL != sproto)
    {
        inchi_free( sproto );
    }
    if (( ret == INCHIKEY_OK ) && ( ret1 != INCHIKEY_OK ))
    {
        ret = ret1;
    }

    return ret;
}


/****************************************************************************
Check if the string represents valid InChIKey.
****************************************************************************/
EXPIMP_TEMPLATE INCHI_API int INCHI_DECL CheckINCHIKey( const char *szINCHIKey )
{
    size_t slen, j;


    slen = strlen( szINCHIKey );


    /* Proper length is 27  */
    if (slen != 27)
    {
        return INCHIKEY_INVALID_LENGTH;
    }

    /* Should have dash in 14-th position */
    if (szINCHIKey[14] != '-')
    {
        return INCHIKEY_INVALID_LAYOUT;
    }
    /* Should have dash in 25-th position */
    if (szINCHIKey[25] != '-')
    {
        return INCHIKEY_INVALID_LAYOUT;
    }

    /* All other should be uppercase */
    for (j = 0; j < 14; j++)
    {
        if (!isbase26( szINCHIKey[j] ))
        {
            return INCHIKEY_INVALID_LAYOUT; /* first block */
        }
    }
    for (j = 15; j < 25; j++)
    {
        if (!isbase26( szINCHIKey[j] ))
        {
            return INCHIKEY_INVALID_LAYOUT; /* second block */
        }
    }

    if (!isbase26( szINCHIKey[26] ))
    {
        return INCHIKEY_INVALID_LAYOUT; /* (de)protonation flag */
    }


    /* No 'E' may appear in 0,3,6,and 9 positions of the 1st block ... */
    for (j = 0; j < 10; j += 3)
    {
        if (szINCHIKey[j] == 'E')
        {
            return INCHIKEY_INVALID_LAYOUT;
        }
    }
    /* ... and 0 and 3 pos. of the second block. */
    for (j = 15; j < 19; j += 3)
    {
        if (szINCHIKey[j] == 'E')
        {
            return INCHIKEY_INVALID_LAYOUT;
        }
    }

    /* Check for version (only 1 allowed) */
    if (szINCHIKey[24] != 'A')
    {
        return INCHIKEY_INVALID_VERSION;
    }

    /* Check for standard-ness */
    if (szINCHIKey[23] == 'S')
    {
        return INCHIKEY_VALID_STANDARD;
    }
    else if (szINCHIKey[23] == 'N')
    {
        return INCHIKEY_VALID_NON_STANDARD;
    }
    else
    {
        return INCHIKEY_INVALID_LAYOUT;
    }

}


/****************************************************************************/
void fprint_digest( FILE* fw, const char *header, unsigned char *a )
{
    size_t i, bytelen = 32;
    fprintf( fw, "%s\n", header );
    for (i = 0; i < bytelen; i++)
    {
        fprintf( fw, "%02x ", a[i] );
    }
    fprintf( fw, "\n" );
}


/****************************************************************************/

#if ( defined( _WIN32 ) && defined( _MSC_VER ) && _MSC_VER >= 800 && defined(_USRDLL) && defined(BUILD_LINK_AS_DLL) )
    /* Win32 & MS VC ++, compile and link as a DLL */
/*********************************************************/
/*   C calling conventions export from Win32 dll         */
/*********************************************************/
/* prototypes */
#ifndef COMPILE_ALL_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif

    int cdecl_CheckINCHIKey( const char *szINCHIKey );
    int cdecl_GetINCHIKeyFromINCHI( const char* szINCHISource,
                                   const int xtra1, const int xtra2,
                                   char* szINCHIKey, char* szXtra1, char* szXtra2 );
    int cdecl_GetStdINCHIKeyFromStdINCHI( const char* szINCHISource, char* szINCHIKey );

#ifndef COMPILE_ALL_CPP
#ifdef __cplusplus
}
#endif
#endif

/* implementation */
/* libinchi.def provides export without cdecl_ prefixes */

/****************************************************************************/
int cdecl_GetStdINCHIKeyFromStdINCHI( const char* szINCHISource, char* szINCHIKey )
{
    return GetStdINCHIKeyFromStdINCHI( szINCHISource, szINCHIKey );
}


/****************************************************************************/
int cdecl_CheckINCHIKey( const char *szINCHIKey )
{
    return CheckINCHIKey( szINCHIKey );
}


/****************************************************************************/
int cdecl_GetINCHIKeyFromINCHI( const char* szINCHISource,
                               const int xtra1, const int xtra2,
                               char* szINCHIKey, char* szXtra1, char* szXtra2 )
{
    return GetINCHIKeyFromINCHI( szINCHISource, xtra1, xtra2, szINCHIKey, szXtra1, szXtra2 );
}
#endif



#if ( defined(__GNUC__) && __GNUC__ >= 3 && defined(__MINGW32__) && defined(_WIN32) )
#include <windows.h>


/****************************************************************************
   Pacal calling conventions export from Win32 dll
****************************************************************************/
#ifndef COMPILE_ALL_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif
/* prototypes */
    int PASCAL pasc_CheckINCHIKey( const char *szINCHIKey );
#ifndef COMPILE_ALL_CPP
#ifdef __cplusplus
}
#endif
#endif

/* implementation */
/* libinchi.def provides export without PASCAL pasc_ prefixes */

/****************************************************************************/
int PASCAL pasc_CheckINCHIKey( const char *szINCHIKey )
{
    return CheckINCHIKey( szINCHIKey );
}

#endif

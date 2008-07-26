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
 
 
/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    InChIKey: calculation of hash for InChI string

    Uses truncated SHA-256 function.
    SHA-256 implementation: Copyright (C) 2003-2006  Christophe Devine, 
                            see files sha2.c, sha2.h.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/

#ifdef _MSC_VER
#if _MSC_VER > 1000
#pragma warning( disable : 4996 )
#endif
#endif

#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "sha2.h"
#include "ikey_base26.h"

#include "ichisize.h"
#include "mode.h"
#include "inchi_api.h"



/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Local options 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/

#define INCHIKEY_DEBUG 0 /* 2 */

#define INCHIKEY_FLAG_OK 0
#define INCHIKEY_NOT_VALID_FLAG 1

/*^^^ InChI PREFIX */
#define INCHI_STRING_PREFIX "InChI="
#define LEN_INCHI_STRING_PREFIX 6

enum { MINOUTLENGTH=256 };  



/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                            Local functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/

const char get_inchikey_flag_char(const char *szINCHISource);

int decode_inchikey_flag_char(const char c,
                              int *version, 
                              int *is_stereo,
                              int *is_fixH,
                              int *is_isotopic);

void fprint_digest(FILE* fw, const char *header, unsigned char *a);





/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                            EXPORTED FUNCTIONS
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/



/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Check if the string represents valid InChIKey.          
Checks both proper letter layout and match of check-character.

Input:
        szINCHIKey
            InChIKey string 
Returns:
        Success/errors codes
                     0  key is valid 
                   !=0  invalid, possible values are:
                        INCHIKEY_INVALID_LENGTH
                        INCHIKEY_INVALID_LAYOUT
                        INCHIKEY_INVALID_CHECKSUM
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
EXPIMP_TEMPLATE INCHI_API int INCHI_DECL CheckINCHIKey(const char *szINCHIKey)
{
size_t slen, j;
int checksum, checkchar;
int keytype = INCHIKEY_VALID;
char str[255];


    slen = strlen(szINCHIKey);

    
    /*^^^ Proper length is 25  */
    if (slen != 25) 
        return INCHIKEY_INVALID_LENGTH;


    
    /*^^^ Should have dash in 14-th position */   

    if (szINCHIKey[14] !='-')   
        return INCHIKEY_INVALID_LAYOUT;
    
    /*^^^ All other should be uppercase */  

    for (j = 0;  j < 14; j++)   
        if ( !isbase26(szINCHIKey[j]) ) 
            return INCHIKEY_INVALID_LAYOUT; /* first block */
    for (j = 15; j < 25; j++)   
        if ( !isbase26(szINCHIKey[j]) ) 
            return INCHIKEY_INVALID_LAYOUT; /* second block */
    
    
    /*^^^ No 'E' may appear in 0,3,6,and 9 positions of the 1st block ... */ 

    for (j=0; j <10; j+=3)  
        if (szINCHIKey[j]=='E')        
            return INCHIKEY_INVALID_LAYOUT;
    
    /*^^^ ... and 0 and 3 pos. of the second block. */
    for (j=15; j <19; j+=3) 
        if (szINCHIKey[j]=='E')        
            return INCHIKEY_INVALID_LAYOUT;
    



    /*^^^ Calculate checksum and compare with the check-character */

    memset(str, 0, sizeof(str));
    for (j = 0;   j < 14; j++)      str[j] =    szINCHIKey[j];  /* first block, 12 letters */
    for (j = 15;  j < 24; j++)      str[j-1] =  szINCHIKey[j];  /* second block, 6 letters */
                                    str[23] =   '\0';
    checksum  = base26_checksum(str);
    checkchar = szINCHIKey[24]; 
    if (checksum!= checkchar)       
        return INCHIKEY_INVALID_CHECKSUM;


    /*^^^ All tests passed. */
    return keytype;
}







/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Calculate InChIKey by InChI string. 

Input:
        szINCHISource
            source InChI string 
        INCHIKeyType
            InChIKey type, always = 1 (rerserved for future)

Output:
        szINCHIKey
            InChIKey string 

NB: the only attempt to check if input szINCHISource reperesents valid InChI identifier is as follows:
    string should
        - starts with "InChI=1/"
        - then contains at least one a..Z0..9 or '/'

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
EXPIMP_TEMPLATE INCHI_API int INCHI_DECL GetINCHIKeyFromINCHI(const char* szINCHISource, 
                                                              char* szINCHIKey)
{
int ret = INCHIKEY_OK;
int cn;
char flagchar, checkchar;
size_t slen, i, j, ncp, pos_slash1=0;
char *str = NULL, *smajor = NULL, *sminor = NULL, 
     *stmp = NULL, tmp[MINOUTLENGTH]; 
unsigned char 
    digest_major[32], digest_minor[32];     



    /*^^^ Check the requested key length */




    /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    Check if InChI string by very minimal requirements:
        - starts with INCHI_STRING_PREFIX,
        - then contains digit (1 presently now),
        - then contains '/',                
        - then contains at least one a..Z0.9 or slash.  
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/    

    if (szINCHISource==NULL)
        return INCHIKEY_EMPTY_INPUT;
    
    slen = strlen(szINCHISource);
    
    if (slen<LEN_INCHI_STRING_PREFIX+3)                         
        return INCHIKEY_NOT_INCHI_INPUT;    
    
    if (memcmp(szINCHISource,INCHI_STRING_PREFIX,LEN_INCHI_STRING_PREFIX))    
        return INCHIKEY_NOT_INCHI_INPUT;
    
    /* if (!isdigit(szINCHISource[LEN_INCHI_STRING_PREFIX]) )  */
    if ( szINCHISource[LEN_INCHI_STRING_PREFIX] != '1' )  
        return INCHIKEY_NOT_INCHI_INPUT;
    
    if (szINCHISource[LEN_INCHI_STRING_PREFIX+1]!='/')            
        return INCHIKEY_NOT_INCHI_INPUT;
    
    if (!isalnum(szINCHISource[LEN_INCHI_STRING_PREFIX+2]) &&
        (szINCHISource[LEN_INCHI_STRING_PREFIX+2]!='/')     ) 
        return INCHIKEY_NOT_INCHI_INPUT;    
    

    /*^^^ Use local copy of the source. */
    
    str = calloc( slen+1, sizeof(char));
    if (NULL==str)
    { 
        ret = INCHIKEY_NOT_ENOUGH_MEMORY; 
        goto fin; 
    }

    strcpy(str, szINCHISource);


    /*^^^ Trim it. */
    
    for (j = slen-1; j > LEN_INCHI_STRING_PREFIX; j--) 
    { 
        cn = str[j]; 
        switch (cn) 
        { 
            case '\r': 
            case '\n':  continue; 
            default:    break; 
        } 
        break; 
    }
    j++;
    str[j]='\0'; 
    slen = strlen(str);

    
    /*^^^ Make buffers. */

    smajor = calloc( slen+1, sizeof(char));
    if (NULL==smajor)
        { ret = INCHIKEY_NOT_ENOUGH_MEMORY; goto fin; }
    sminor = calloc( 2*slen + 2, sizeof(char)); /* we may double the length ... */
    if (NULL==sminor)
        { ret = INCHIKEY_NOT_ENOUGH_MEMORY; goto fin; }
    stmp = calloc( slen+1, sizeof(char));
    if (NULL==stmp)
        { ret = INCHIKEY_NOT_ENOUGH_MEMORY; goto fin; }


    
    szINCHIKey[0] = '\0';

    
    /*^^^ Extract the major block. */
    
    smajor[0]  = '\0';

    pos_slash1=0;
    for (j = LEN_INCHI_STRING_PREFIX; j < slen; j++)
    {
        if (str[j]=='/')
        {
            pos_slash1 = j; 
            break;
        }
    }
    if (pos_slash1==0) 
    { 
        ret = INCHIKEY_NOT_ENOUGH_MEMORY; 
        goto fin; 
    }


    for (j = pos_slash1 + 1; j < slen-1; j++)
    {
        if (str[j]=='/')    /* case 'p': protons go to minor hash */
        {
            cn = str[j+1];
            switch (cn) 
            { 
                case 'c': case 'h': case 'q':   continue; 
                default:                        break; 
            }
            break;
        }
    }       
    j++;
    if (j==slen) 
        j++;
    else 
        j--;

    ncp=0;
    ncp = j-pos_slash1-1; 
    
    
    /*^^^ Trim 'InChi=', keep '1/'  (may be '2/' :) */
    
    memcpy(smajor,str+pos_slash1+1, ncp*sizeof(str[0]));
    smajor[ncp]='\0';



    /*^^^ Extract the minor block. */
    
    if (j != slen+1)    /*^^^ check that something exists at right.*/
    {
        ncp = slen-j; 
        memcpy(sminor,str+j, (ncp)*sizeof(str[0]));
        sminor[ncp]='\0';
    }
    else 
        sminor[0]='\0'; 
            
#if INCHIKEY_DEBUG 
    fprintf(stdout,"Source:  {%-s}\n",str);
    fprintf(stdout,"SMajor:  {%-s}\n",smajor);
    fprintf(stdout,"SMinor:  {%-s}\n",sminor);
#endif      



    /*^^^ Compute and compose the InChIKey string. */

        
    
    /*^^^ Major hash sub-string. */
    
    for( i = 0; i < 32; i++ ) 
        digest_major[i] = 0;            



    sha2_csum( (unsigned char *) smajor, (int) strlen(smajor), digest_major );
    
    

#if (INCHIKEY_DEBUG>1)
        fprint_digest(stderr, "Major hash, full SHA-256",digest_major);         
#endif  

    
    sprintf(tmp,"%-.3s%-.3s%-.3s%-.3s%-.2s", 
                base26_triplet_1(digest_major), base26_triplet_2(digest_major),
                base26_triplet_3(digest_major), base26_triplet_4(digest_major),
                base26_dublet_for_bits_56_to_64(digest_major));
    strcat(szINCHIKey, tmp);    




    /*^^^ Minor hash sub-string. */

    for( i = 0; i < 32; i++ ) 
        digest_minor[i] = 0; 

    slen = strlen(sminor); 
    if ((slen>0)&&(slen<255)) 
    { 
        strcpy(stmp, sminor); 
        strcpy(sminor+slen,stmp); 
    }



    sha2_csum( (unsigned char *) sminor, (int) strlen(sminor), digest_minor );



#if (INCHIKEY_DEBUG>1)  
        fprint_digest(stderr, "Minor hash, full SHA-256",digest_minor);
#endif

    
    strcat(szINCHIKey, "-");
    sprintf(tmp,"%-.3s%-.3s%-.2s", 
        base26_triplet_1(digest_minor), 
        base26_triplet_2(digest_minor), 
        base26_dublet_for_bits_28_to_36(digest_minor));  
    strcat(szINCHIKey, tmp); 



    /*^^^ Append flag character */

    slen = strlen(szINCHIKey);
    flagchar = get_inchikey_flag_char(szINCHISource);
    if (flagchar=='Z')
    { 
        ret = INCHIKEY_ERROR_IN_FLAG_CHAR; 
        goto fin; 
    }

    szINCHIKey[slen] = flagchar; 
    szINCHIKey[slen+1] = '\0';


    /*^^^ Calculate check character and insert it before the flag character */      
    
    checkchar = base26_checksum(szINCHIKey);
    szINCHIKey[slen+1] = checkchar; 
    szINCHIKey[slen+2] = '\0';



fin:if (NULL!=str)      free(str);
    if (NULL!=smajor)   free(smajor);
    if (NULL!=sminor)   free(sminor);
    if (NULL!=stmp)     free(stmp);

    return ret;
}









/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
void fprint_digest(FILE* fw, const char *header, unsigned char *a)
{
size_t i, bytelen = 32;
    fprintf(fw,"%s\n", header);
    for( i = 0; i < bytelen; i++ )
        fprintf(fw,"%02x ", a[i]);
    fprintf(fw,"\n" );
}


/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Calculate flag character for InChIKey string    .
Return:     flag character or 'Z' on errors.

NB: caller is responsible to check that the string is appropriate one.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
const char get_inchikey_flag_char(const char *szINCHISource)
{
/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        We have a base-26 letter A..Z to capture flags.
        
        Note that 2*2*2*8 = 24 so we may place in the letter (3 bits + 1 trit). 

        trit 0 - 0 if InChI version 1.0, otherwise 1 or 2
        bit  1 - 1 if InChI contains stereo information, otherwise 0
        bit  2 - 1 if InChI contains fixed H layer, otherwise 0
        bit  3 - 1 if InChI contains isotopic information, otherwise 0
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
    size_t slen, j;
    unsigned char flag = 0x00;
    unsigned char cver,c,c1; 


    slen = strlen(szINCHISource);
    if (slen<LEN_INCHI_STRING_PREFIX+3)                         
        return 'Z';


    /*^^^ Version number (assume integer) */    

    cver = szINCHISource[LEN_INCHI_STRING_PREFIX];
    if (!isdigit(cver)) 
        return 'Z';

    
    /*^^^ Search for the layers */
    
    for (j=0; j<slen-1; j++) 
    {
        c = szINCHISource[j]; 
        c1 = szINCHISource[j+1];
        if (c=='/')
        {
            switch (c1)
            {
                /*^^^   Have stereo? */
                case 't': 
                case 'b': 
                case 'm': 
                case 's':   flag |= 0x01; 
                            break;
                /*^^^ Have fix-H? */
                case 'f':
                             flag |= 0x02; 
                             break;
                /*^^^ Have isotopic ? */
                case 'i':   flag |= 0x04; 
                            break;
                default:    break;
            }
        }
    }

    switch (cver)
    {
        case '1':   return "ABCDEFGH"[flag];
        case '2':   return "IJKLMNOP"[flag];
        case '3':   return "QRSTUVWX"[flag];
        default :   break;
    }

    return 'Z'; /* error */
}




/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Decode flag character of InChI key string
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
int decode_inchikey_flag_char(const char c,
                              int *version, 
                              int *is_stereo,
                              int *is_fixH,
                              int *is_isotopic)
{
    unsigned char flag;
    unsigned char a='A', i='I', q='Q';

    if (!isbase26(c))   return INCHIKEY_NOT_VALID_FLAG;
    if (c>'X')          return INCHIKEY_NOT_VALID_FLAG;


    *is_stereo = *is_fixH = *is_isotopic = 0;   
    
    if (c>=q)       { *version = 3; flag = c - q; }
    else if (c>=i)  { *version = 2; flag = c - i; }
    else            { *version = 1; flag = c - a; }
    *is_stereo   =  flag&0x01;
    *is_fixH     =  (flag&0x02)  >> 1;
    *is_isotopic =  (flag&0x04)  >> 2;
                                                 
    return INCHIKEY_FLAG_OK;

}



/********************************************************/
int cdecl_CheckINCHIKey(const char* szINCHIKey)
{
    return CheckINCHIKey(szINCHIKey);
}
int cdecl_GetINCHIKeyFromINCHI(const char* szINCHISource, char* szINCHIKey)
{
    return GetINCHIKeyFromINCHI(szINCHISource, szINCHIKey);
}

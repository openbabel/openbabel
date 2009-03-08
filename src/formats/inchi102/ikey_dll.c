/*
 * International Chemical Identifier (InChI)
 * Version 1
 * Software version 1.02
 * January 10, 2009
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

#include "util.h"




/*^^^ Post-1.02b - added check of source InChI string validity (API function GetInChIKey) 
      through calling inchi2inchi and comparing the result with source.
      This option significantly slows down key generation but adds reliability.      
      Comment the #define below if source InChI strings are absolutely trusted
      (as it takes place at cInChI build where InChI's are just generated internally). 
*/
#ifdef INCHI_LIBRARY
#ifndef ONLY_STDINCHI_STDKEY
#define INCHIKEY_CHECK_INCHI_WITH_I2I 1
#endif
#endif


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
Return code for CheckINCHIKey102b
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
#define INCHIKEY_VALID 0
#define INCHIKEY_INVALID_LENGTH 1
#define INCHIKEY_INVALID_LAYOUT 2
#define INCHIKEY_INVALID_CHECKSUM 3


/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                            Local functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/

static int INCHI_DECL GetINCHIKeyFromINCHI(const char* szINCHISource, char* szINCHIKey);

static int INCHI_DECL CheckINCHIKey(const char *szINCHIKey);
static int cdecl_CheckINCHIKey(const char* szINCHIKey);

static char get_inchikey102b_flag_char(const char *szINCHISource);

static int decode_inchikey102b_flag_char(const char c,
                              int *version, 
                              int *is_stereo,
                              int *is_fixH,
                              int *is_isotopic);

static int cdecl_CheckINCHIKey102b(const char* szINCHIKey); 
static int INCHI_DECL GetINCHIKey102bFromINCHI(const char* szINCHISource, char* szINCHIKey);
static int cdecl_GetINCHIKey102bFromINCHI(const char* szINCHISource, char* szINCHIKey); 

void fprint_digest(FILE* fw, const char *header, unsigned char *a);
void extract_trimmed_inchi(char ** buf, const char *str, size_t slen);



/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                            EXPORTED FUNCTIONS
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/

EXPIMP_TEMPLATE INCHI_API int INCHI_DECL GetStdINCHIKeyFromStdINCHI(const char* szINCHISource, 
                                                                 char* szINCHIKey)
{
#ifdef ONLY_STDINCHI_STDKEY
	return GetINCHIKeyFromINCHI(szINCHISource, szINCHIKey);
#else
	/* to be filled */
	return inchi_Ret_ERROR;
#endif
}

static int INCHI_DECL GetINCHIKeyFromINCHI(const char* szINCHISource, char* szINCHIKey)
{
int ret = INCHIKEY_OK;
int ret1 = INCHIKEY_OK;
int cn;
size_t slen, i, j, jproto=0, ncp, pos_slash1=0;
char *str = NULL, *smajor = NULL, *sminor = NULL, 
     *sproto=NULL,
     *stmp = NULL, tmp[MINOUTLENGTH]; 
unsigned char 
    digest_major[32], digest_minor[32];     

char flagstd = 'S', /* standard key */
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


/*^^^ Currently disable check of InChI string via call to inchi2inchi 
      (until testing of non-std InChI options ^^^*/
#if 0
        inchi_InputINCHI    inchi_inp;
        inchi_Output        inchi_out;
        int ret_i2i;
#endif


#ifdef ONLY_STDINCHI_STDKEY

    /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    Check if std InChI string, by very minimal requirements:
        - starts with INCHI_STRING_PREFIX,
        - then contains digit (1 presently now),
        - then contains 'S',                
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
    
    if (szINCHISource[LEN_INCHI_STRING_PREFIX+1]!='S')            
        return INCHIKEY_INVALID_STD_INCHI;

    if (szINCHISource[LEN_INCHI_STRING_PREFIX+2]!='/')            
        return INCHIKEY_NOT_INCHI_INPUT;

    if (!isalnum(szINCHISource[LEN_INCHI_STRING_PREFIX+3]) &&
        (szINCHISource[LEN_INCHI_STRING_PREFIX+3]!='/')     ) 
        return INCHIKEY_NOT_INCHI_INPUT;    
    

    /*^^^ Use local copy of the source. */
    

    extract_trimmed_inchi(&str, szINCHISource, slen);


    if (NULL==str)
    { 
        ret = INCHIKEY_NOT_ENOUGH_MEMORY; 
        goto fin; 
    }
    slen = strlen(str);

/*^^^ Currently disable check of InChI string via call to inchi2inchi 
      (until testing of non-std InChI options ^^^*/
#if 0
    inchi_inp.szInChI = str;

#ifdef _WIN32
    inchi_inp.szOptions  = "/FixedH /RecMet /SPXYZ /SAsXYZ /FB /FB2 /FNUD /NEWPS";
#else
    inchi_inp.szOptions  = "-FixedH -RecMet -SPXYZ -SAsXYZ -FB -FB2 -FNUD -NEWPS";
#endif


    ret_i2i = GetINCHIfromINCHI(&inchi_inp, &inchi_out);

    if ( ((ret_i2i!=inchi_Ret_OKAY) && (ret_i2i!=inchi_Ret_WARNING)) || !inchi_out.szInChI )
    {
        ret = INCHIKEY_INCHI_REVERSAL_FAIL;
        goto fin; 
    }
    else
    {
        if (strcmp(inchi_inp.szInChI, inchi_out.szInChI))
        {
            /* will calculate key but set error/warning return code */
            ret1 = INCHIKEY_INCHI_REVERSED_NOT_THE_SAME;
            /* goto fin;  */
        }
    }
#endif
    
    /*^^^ Make buffers. */

    smajor = (char*) inchi_calloc( slen+1, sizeof(char)); 
    if (NULL==smajor)
        { ret = INCHIKEY_NOT_ENOUGH_MEMORY; goto fin; }
    sminor = (char*) inchi_calloc( 2*slen + 2, sizeof(char)); /* we may double the length ... */ 
    if (NULL==sminor)
        { ret = INCHIKEY_NOT_ENOUGH_MEMORY; goto fin; }
    stmp = (char*) inchi_calloc( slen+1, sizeof(char)); 
    if (NULL==stmp)
        { ret = INCHIKEY_NOT_ENOUGH_MEMORY; goto fin; }
    sproto = (char*) inchi_calloc( slen+1, sizeof(char)); 
    if (NULL==sproto)
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
        if (str[j]=='/')    
        {
            cn = str[j+1];
            switch (cn) 
            { 
                /* anything allowed from a major part */
                case 'c': case 'h': case 'q':   continue; 
                
                /* "/p"; protons now go to to special string, not to minor hash */
                case 'p':						jproto = j;
                                                continue;
                
                /* "/f",  "/r" : may not occur */
                case 'f': case 'r':				ret = INCHIKEY_INVALID_STD_INCHI; 
                                                goto fin;
                
                /* anything allowed from a minor part */
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

    if (jproto)
        ncp = jproto - pos_slash1 - 1; 
    else
        ncp = j - pos_slash1 - 1; 
    
    
    /*^^^ Trim 'InChI=X/' */
    
    memcpy(smajor,str+pos_slash1+1, ncp*sizeof(str[0]));
    smajor[ncp]='\0';
    
    /* Treat protonization */
    if (jproto)
    {
        /* 2009-01-07 fix bug/typo: assigned incorrect length to the protonation segment of 
        /* source string ( was sproto[ncp]='\0'; should be sproto[lenproto]='\0'; )  */
        int lenproto = j - jproto;
        if (lenproto<3)
        {	
            /* empty "/p", should not occur */
            ret = INCHIKEY_INVALID_STD_INCHI; 
            goto fin;
        }
        
        memcpy(sproto,str+pos_slash1+ncp+1, lenproto*sizeof(str[0]));
        sproto[lenproto]='\0';
       
        nprotons = strtol( sproto+2, NULL, 10 );        
        
        if (nprotons > 0)
        {
            if (nprotons > 12) flagproto = 'A';
            else			   flagproto = pplus[nprotons-1];
        }
        else if (nprotons < 0)
        {
            if (nprotons < -12) flagproto = 'A';
            else				flagproto = pminus[-nprotons-1];
        }
        else
        {
            /* should never occur */
            ret = INCHIKEY_INVALID_STD_INCHI; 
            goto fin;
        }
    }



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
    fprintf(stdout,"SProto:  {%-s}\n",sproto);
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



    /*^^^ Append "standard" flag */
    slen = strlen(szINCHIKey);
    szINCHIKey[slen] = flagstd; 
    /*^^^ Append InChI v.1 flag */
    szINCHIKey[slen+1] = flagver; 
    /*^^^ Append dash  */
    szINCHIKey[slen+2] = '-'; 
    /*^^^ Append protonization flag */
    szINCHIKey[slen+3] = flagproto; 
    szINCHIKey[slen+4] = '\0';


#if INCHIKEY_DEBUG 
    fprintf(stdout,"szINCHIKey:  {%-s}\n",szINCHIKey);
#endif

fin:if (NULL!=str)      inchi_free(str);
    if (NULL!=smajor)   inchi_free(smajor);
    if (NULL!=sminor)   inchi_free(sminor);
    if (NULL!=stmp)     inchi_free(stmp);
    if (NULL!=sproto)   inchi_free(sproto);

/*^^^ Currently disable check of InChI string via call to inchi2inchi 
      (until testing of non-std InChI options ^^^*/
#if 0
#ifdef INCHIKEY_CHECK_INCHI_WITH_I2I
    FreeINCHI (&inchi_out);
#endif
#endif

    if ( (ret==INCHIKEY_OK) && (ret1!=INCHIKEY_OK) )
        ret = ret1;
    return ret;


#else

    /* to be filled with non-std InChI/Key */
    return INCHIKEY_UNKNOWN_ERROR;

#endif /* ifdef ONLY_STDINCHI_STDKEY */
}


/********************************************************/
int cdecl_GetStdINCHIKeyFromStdINCHI(const char* szINCHISource, char* szINCHIKey); 
int cdecl_GetStdINCHIKeyFromStdINCHI(const char* szINCHISource, char* szINCHIKey)
{
    return GetStdINCHIKeyFromStdINCHI(szINCHISource, szINCHIKey);
}


/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Check if the string represents valid InChIKey.          
TO BE FILLED
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
static int INCHI_DECL CheckINCHIKey(const char *szINCHIKey)
{

            return INCHIKEY_VALID;
}

static int cdecl_CheckINCHIKey(const char* szINCHIKey)
{
    return CheckINCHIKey(szINCHIKey);
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






/*^^^ Post-1.02b - more correct trimming of InChI string */
/*
    According to
    http://info-uri.info/registry/OAIHandler?verb=GetRecord&metadataPrefix=reg&identifier=info:inchi/

    An InChI identifier may contain the following characters:

    A-Z
    a-z
    0-9
    ()*+,-./;=?@

    
    Here we consider any character not conforming this specification as a whitespace
    which marks the end of the InChI string.
    For example:
    "InChI=1/Ar%"
    "InChI=1/Ar\n"
    "InChI=1/Ar\r\t"
    all will be trimmed to
    "InChI=1/Ar"

*/
void extract_trimmed_inchi(char ** buf, const char *str, size_t slen)
{
size_t i; 
int c;  
    
    *buf = NULL;

    
    for (i=0; i<slen; i++)
    {
        c = (int) (unsigned char) str[i];
        
        if (c>='A' && c<='Z')   continue; 
        if (c>='a' && c<='z')   continue; 
        if (c>='0' && c<='9')   continue;             
        switch (c) 
        { 
            case '(': case ')': 
            case '*': case '+': 
            case ',': case '-': 
            case '.': case '/': 
            case ';': case '=': 
            case '?': case '@': continue;             
            
            default:            break; 
        }         
        break; 
    }


    *buf = (char*) inchi_calloc(i+1, sizeof(char));   
    memcpy(*buf, str, i);
    (*buf)[i] = '\0';

    return;
}




/*^^^ v.102b heritage follows */




/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Calculate flag character for InChIKey string    .
Return:     flag character or 'Z' on errors.

NB: caller is responsible to check that the string is appropriate one.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
static char get_inchikey102b_flag_char(const char *szINCHISource)
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
static int decode_inchikey102b_flag_char(const char c,
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
static int INCHI_DECL CheckINCHIKey102b(const char *szINCHIKey)
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
static int INCHI_DECL GetINCHIKey102bFromINCHI(const char* szINCHISource, 
                                                              char* szINCHIKey)
{
int ret = INCHIKEY_OK;
int ret1 = INCHIKEY_OK;
int cn;
char flagchar, checkchar;
size_t slen, i, j, ncp, pos_slash1=0;
char *str = NULL, *smajor = NULL, *sminor = NULL, 
     *stmp = NULL, tmp[MINOUTLENGTH]; 
unsigned char 
    digest_major[32], digest_minor[32];     

/*^^^ Post-1.02b - check InChI string supplied to GetInChIKey via call to inchi2inchi */
#ifdef INCHIKEY_CHECK_INCHI_WITH_I2I
        inchi_InputINCHI    inchi_inp;
        inchi_Output        inchi_out;
        int ret_i2i;
#endif



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
    

    extract_trimmed_inchi(&str, szINCHISource, slen);


    if (NULL==str)
    { 
        ret = INCHIKEY_NOT_ENOUGH_MEMORY; 
        goto fin; 
    }
    slen = strlen(str);

/*^^^ Post-1.02b - check InChI string supplied to GetInChIKey via call to inchi2inchi */
#ifdef INCHIKEY_CHECK_INCHI_WITH_I2I
    inchi_inp.szInChI = str;

#ifdef _WIN32
    inchi_inp.szOptions  = "/FixedH /RecMet /SPXYZ /SAsXYZ /FB /FB2 /FNUD /NEWPS";
#else
    inchi_inp.szOptions  = "-FixedH -RecMet -SPXYZ -SAsXYZ -FB -FB2 -FNUD -NEWPS";
#endif


    ret_i2i = GetINCHIfromINCHI(&inchi_inp, &inchi_out);

    if ( ((ret_i2i!=inchi_Ret_OKAY) && (ret_i2i!=inchi_Ret_WARNING)) || !inchi_out.szInChI )
    {
        ret = INCHIKEY_INCHI_REVERSAL_FAIL;
        goto fin; 
    }
    else
    {
        if (strcmp(inchi_inp.szInChI, inchi_out.szInChI))
        {
            /* will calculate key but set error/warning return code */
            ret1 = INCHIKEY_INCHI_REVERSED_NOT_THE_SAME;
            /* goto fin;  */
        }
    }
#endif
    
    /*^^^ Make buffers. */

    smajor = (char*) inchi_calloc( slen+1, sizeof(char)); 
    if (NULL==smajor)
        { ret = INCHIKEY_NOT_ENOUGH_MEMORY; goto fin; }
    sminor = (char*) inchi_calloc( 2*slen + 2, sizeof(char)); /* we may double the length ... */ 
    if (NULL==sminor)
        { ret = INCHIKEY_NOT_ENOUGH_MEMORY; goto fin; }
    stmp = (char*) inchi_calloc( slen+1, sizeof(char)); 
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
    
    
    /*^^^ Trim 'InChI=X/' */
    
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
    flagchar = get_inchikey102b_flag_char(szINCHISource);
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



fin:if (NULL!=str)      inchi_free(str);
    if (NULL!=smajor)   inchi_free(smajor);
    if (NULL!=sminor)   inchi_free(sminor);
    if (NULL!=stmp)     inchi_free(stmp);

/*^^^ Post-1.02b - check InChI string supplied to GetInChIKey via call to inchi2inchi */
#ifdef INCHIKEY_CHECK_INCHI_WITH_I2I
    FreeINCHI (&inchi_out);
#endif

    if ( (ret==INCHIKEY_OK) && (ret1!=INCHIKEY_OK) )
        ret = ret1;
    return ret;
}







/********************************************************/
static int cdecl_CheckINCHIKey102b(const char* szINCHIKey)
{
    return CheckINCHIKey102b(szINCHIKey);
}
static int cdecl_GetINCHIKey102bFromINCHI(const char* szINCHISource, char* szINCHIKey)
{
    return GetINCHIKey102bFromINCHI(szINCHISource, szINCHIKey);
}


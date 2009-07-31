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

                                    InChI - API 1.02-beta

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdarg.h>
#include <errno.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include "inpdef.h"
#include "ichi.h"
#include "strutil.h"
#include "util.h"
#include "ichierr.h"
#include "ichimain.h"
#include "extr_ct.h"

#include "ichicomp.h"
#include "ichitaut.h"
#include "ichinorm.h"


#include "ichisize.h"
#include "mode.h"
#include "inchi_api.h"

#include "inchi_dll_a.h" /* not inchi_api.h as it hides internal data types */





/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Local prototypes.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
int parse_options_string ( char *cmd, const char *argv[], int maxargs );

int ExtractOneStructure( STRUCT_DATA *sd, INPUT_PARMS *ip, char *szTitle,
         inchi_Input *pInp, INCHI_FILE *log_file, INCHI_FILE *output_file, FILE *prb_file,
         ORIG_ATOM_DATA *orig_inp_data, long *num_inp, char *pStr, int nStrLen );



int NormOneStructureINChI(INCHIGEN_DATA *pGenData, void * HGen, 
                          int iINChI, FILE *inp_file);
int CanonOneStructureINChI(INCHIGEN_DATA *pGenData, void * HGen, 
                           int iINChI, FILE *inp_file);




/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
InChI Generator: create generator
Returns handle of generator object or NULL on failure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
void * INCHI_DECL INCHIGEN_Create(void)
{
INCHIGEN_CONTROL * HGen = NULL;


    HGen = inchi_malloc( sizeof(INCHIGEN_CONTROL) );
    
    if (!HGen) 
        return (void *) NULL;

    memset(HGen, 0, sizeof(INCHIGEN_CONTROL));

    /*^^^ Set/init aliases */
    memset(&(HGen->InpParms), 0, sizeof(INPUT_PARMS));
    memset(&(HGen->StructData), 0, sizeof(STRUCT_DATA) );

    HGen->ulTotalProcessingTime = 0;
    HGen->num_err = 0;
    HGen->num_inp = 0;
    HGen->szTitle[0] = '\0';
    
    HGen->pStr = (char*) inchi_malloc(PSTR_BUFFER_SIZE);
    if (!HGen->pStr) 
    { 
        free(HGen);
        return (void *) NULL;   
    }
    HGen->pStr[0] = '\0';
    memset(HGen->inchi_file, 0, sizeof(HGen->inchi_file));

    memset(&(HGen->OrigInpData), 0, sizeof( HGen->OrigInpData ) );
    memset(&(HGen->PrepInpData[0]), 0, 2*sizeof( HGen->PrepInpData[0] ) );

    memset(HGen->pINChI,     0, sizeof(HGen->pINChI) );
    memset(HGen->pINChI_Aux, 0, sizeof(HGen->pINChI_Aux) );


    return (void *) HGen;
}

/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
InChI Generator: initialization stage (accepts a specific structure)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/

int INCHI_DECL INCHIGEN_Setup(void * _HGen, INCHIGEN_DATA * pGenData, inchi_Input * pInp)
                                   

{
int retcode = inchi_Ret_OKAY;

INCHIGEN_CONTROL *HGen = _HGen;

ORIG_ATOM_DATA *orig_inp_data = &(HGen->OrigInpData);
STRUCT_DATA *sd = &(HGen->StructData);
INPUT_PARMS *ip = &(HGen->InpParms);
INCHI_FILE *log_file = HGen->inchi_file+1;  
FILE *prb_file=NULL;

const char *argv[INCHI_MAX_NUM_ARG+1];
int   argc;
char *szOptions = NULL;
char szSdfDataValue[MAX_SDF_VALUE+1];
int bReleaseVersion = bRELEASE_VERSION;
unsigned long  ulDisplTime = 0;    /*  infinite, milliseconds */
int p;





    /*^^^ Make allocs/inits */

    if (!pGenData) 
    { 
        retcode = _IS_ERROR; 
        goto ret; 
    }    
    memset(pGenData, 0, sizeof(*pGenData));
    




    /*^^^ Parse 'command-line' options and fill internal INPUT_PARMS structure */

    if ( pInp && pInp->szOptions ) 
    {        
        szOptions = (char*)inchi_malloc( strlen(pInp->szOptions) + 1 );
        if (!szOptions) 
            return _IS_FATAL;   /*^^^ Not enough memory.... */
        else
        {
            /*^^^ Parse. */
            strcpy( szOptions, pInp->szOptions );
            argc = parse_options_string ( szOptions, argv, INCHI_MAX_NUM_ARG );
        } 
    } 
    else    
    {
        /*^^^ Got NULL options string or NULL 'pInp', will use defaults. */
        argc = 1;
        argv[0] = "";
        argv[1] = NULL;
    }


    if ( argc == 1
#ifdef INCHI_LIBRARY
        && (!pInp || pInp->num_atoms <= 0 || !pInp->atom)
#endif        
        || argc==2 && ( argv[1][0]==INCHI_OPTION_PREFX ) &&
        (!strcmp(argv[1]+1, "?") || !stricmp(argv[1]+1, "help") ) ) 
    {
    
        HelpCommandLineParms(log_file);
        memset( log_file, 0, sizeof(*log_file) );
        return _IS_EOF;
    }


    memset( szSdfDataValue    , 0, sizeof( szSdfDataValue    ) );


    /*^^^ Decrypt. */
    /*^^^ NB: ReadCommandLineParms resides in _header_ file, ichiparm.h */

    retcode = ReadCommandLineParms( argc, argv, ip, szSdfDataValue, &ulDisplTime, bReleaseVersion, log_file );    

    if (szOptions) 
        inchi_free( szOptions );

    /* INChI DLL specific */
    ip->bNoStructLabels = 1;

    if ( 0 > retcode)   goto ret;
    
    if ( ip->bNoStructLabels ) 
    {
        ip->pSdfLabel = NULL;
        ip->pSdfValue = NULL;
    } 
    else
    if ( ip->nInputType == INPUT_INCHI_XML || ip->nInputType == INPUT_INCHI_PLAIN  || ip->nInputType == INPUT_CMLFILE ) 
    {
        /* the input may contain both the header and the label of the structure */
        if ( !ip->pSdfLabel ) 
            ip->pSdfLabel  = ip->szSdfDataHeader;
        if ( !ip->pSdfValue )
            ip->pSdfValue  = szSdfDataValue;
    }

    
    if (retcode!=inchi_Ret_OKAY) goto ret; 

    PrintInputParms( log_file, ip);    



    /*^^^ Extract the structure */

    retcode = ExtractOneStructure(  sd, ip, 
                                    HGen->szTitle, pInp, 
                                    log_file, 
                                    HGen->inchi_file, /* output_file */
                                    prb_file,
                                    orig_inp_data, 
                                    &(HGen->num_inp), 
                                    HGen->pStr, 
                                    PSTR_BUFFER_SIZE );


ret:switch (retcode) 
    { 
        case _IS_OKAY   : retcode = inchi_Ret_OKAY   ; HGen->init_passed = 1; break;    /* Success; break; no errors or warnings */
        
        case _IS_ERROR  : (HGen->num_err)++;  retcode = inchi_Ret_ERROR  ; break; 
                                                                /* Error: no INChI has been created */
        case _IS_FATAL  : (HGen->num_err)++;  retcode = inchi_Ret_FATAL  ; break; 
                                                                /* Severe error: no INChI has been created 
                                                                (typically; break; memory allocation failed) */
        case _IS_SKIP   : retcode = inchi_Ret_SKIP   ; break;   /* not used in INChI dll */
        case _IS_EOF    : retcode = inchi_Ret_EOF    ; break;   /* no structural data has been provided */
        case _IS_WARNING: retcode = inchi_Ret_WARNING; HGen->init_passed = 1; break;    /* Success; break; warning(s) issued */
        case _IS_UNKNOWN:
        default         : retcode = inchi_Ret_UNKNOWN; break;   /* Unlnown program error */
    }

    strcpy(pGenData->pStrErrStruct, sd->pStrErrStruct);
    for (p=0; p < INCHI_NUM; p++)  
        pGenData->num_components[p] = sd->num_components[p];

    return retcode;
}




/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Get normalized form of the structure.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
int INCHI_DECL INCHIGEN_DoNormalization( void * _HGen, INCHIGEN_DATA *pGenData)
{
int nRet=0, nRet1=0;
/* int maxINChI=0; */

FILE *inp_file = NULL;


INCHIGEN_CONTROL * HGen = _HGen;
INPUT_PARMS *ip = &(HGen->InpParms);
STRUCT_DATA *sd = &(HGen->StructData);
NORM_CANON_FLAGS *pncFlags = &(HGen->ncFlags);
INCHI_FILE *output_file = HGen->inchi_file;
ORIG_ATOM_DATA *orig_inp_data = &(HGen->OrigInpData);
ORIG_STRUCT      *pOrigStruct = NULL;

int k;

#if ( RING2CHAIN == 1 || UNDERIVATIZE == 1 )
int ret1=0, ret2=0;
#endif




/*^^^ Set debug output */
#if (TRACE_MEMORY_LEAKS == 1)

    _CrtSetDbgFlag(_CRTDBG_CHECK_ALWAYS_DF | _CRTDBG_LEAK_CHECK_DF | _CRTDBG_ALLOC_MEM_DF);
/* for execution outside the VC++ debugger uncomment one of the following two */
#ifdef MY_REPORT_FILE 
   _CrtSetReportMode( _CRT_WARN, _CRTDBG_MODE_FILE );
   _CrtSetReportFile( _CRT_WARN, MY_REPORT_FILE );
   _CrtSetReportMode( _CRT_ERROR, _CRTDBG_MODE_FILE );
   _CrtSetReportFile( _CRT_ERROR, MY_REPORT_FILE );
   _CrtSetReportMode( _CRT_ASSERT, _CRTDBG_MODE_FILE );
   _CrtSetReportFile( _CRT_ASSERT, MY_REPORT_FILE );
#else
    _CrtSetReportMode(_CRT_WARN | _CRT_ERROR, _CRTDBG_MODE_DEBUG);
#endif 
#if ( !defined(__STDC__) || __STDC__ != 1 )
    /* turn on floating point exceptions */
    {
        /* Get the default control word. */
        int cw = _controlfp( 0,0 );
        /* Set the exception masks OFF, turn exceptions on. */
        /*cw &=~(EM_OVERFLOW|EM_UNDERFLOW|EM_INEXACT|EM_ZERODIVIDE|EM_DENORMAL);*/
        cw &=~(EM_OVERFLOW|EM_UNDERFLOW|EM_ZERODIVIDE|EM_DENORMAL);
        /* Set the control word. */
        _controlfp( cw, MCW_EM );
    }
#endif /*^^^ ( !defined(__STDC__) || __STDC__ != 1 ) */

#endif /*^^^ (TRACE_MEMORY_LEAKS == 1) */



    if (HGen->init_passed==0)
    {
        AddMOLfileError(sd->pStrErrStruct, "InChI generator not initialized");
        sd->nStructReadError =  99;
        sd->nErrorType = _IS_ERROR;
        nRet = _IS_ERROR;
        goto exit_function;
    }

    sd->bUserQuitComponent = 0;
    sd->bUserQuitComponentDisplay = 0;
    memset( HGen->composite_norm_data, 0, sizeof(HGen->composite_norm_data) );
    memset( pncFlags, 0, sizeof(*pncFlags) );

    /* for testing only */
#if( REMOVE_ION_PAIRS_ORIG_STRU == 1 )
    fix_odd_things( orig_inp_data->num_inp_atoms, orig_inp_data->at, 0 );
#endif
#if( UNDERIVATIZE == 1 )  /***** post v.1 feature *****/
    if ( ip->bUnderivatize && 0 > (ret2=underivatize( orig_inp_data )) ) 
    {
        long num_inp2 = HGen->num_inp;
        AddMOLfileError(sd->pStrErrStruct, "Underivatization error");
        sd->nStructReadError =  99;
        sd->nErrorType = _IS_ERROR;
        nRet = _IS_ERROR;
        TreatReadTheStructureErrors( sd, ip, LOG_MASK_ALL, inp_file, log_file, output_file, prb_file,
                                        prep_inp_data, &num_inp2, HGen->pStr, PSTR_BUFFER_SIZE);
        goto exit_function; /* output only if derivatives found */
    }
#endif /* UNDERIVATIZE == 1 */
#if( RING2CHAIN == 1 )  /***** post v.1 feature *****/
    if ( ip->bRing2Chain && 0 > (ret1 = Ring2Chain( orig_inp_data )) ) 
    {
        long num_inp2 = HGen->num_inp;
        AddMOLfileError(sd->pStrErrStruct, "Ring to chain error");
        sd->nStructReadError =  99;
        sd->nErrorType = _IS_ERROR;
        nRet = _IS_ERROR;
        TreatReadTheStructureErrors( sd, ip, LOG_MASK_ALL, inp_file, log_file, output_file, prb_file,
                                        prep_inp_data, &num_inp2, HGen->pStr, PSTR_BUFFER_SIZE);
        goto exit_function; /* output only if derivatives found */
    }
#endif /* RING2CHAIN == 1 */
#if ( RING2CHAIN == 1 || UNDERIVATIZE == 1 )  /***** post v.1 feature *****/
    if ( ip->bIngnoreUnchanged && !ret1 && !ret2 ) 
    {
        goto exit_function; /* output only if derivatives or ring/chain found */
    }
#endif /* RING2CHAIN == 1 || UNDERIVATIZE == 1 */


    /***** output MOLfile ***************/
    if ( ip->bINChIOutputOptions & INCHI_OUT_SDFILE_ONLY  ) 
    {
        char szNumber[32];
        int ret1a=0, ret2a=0; /* for derivatives and ring-chain */
        ret1a = sprintf(szNumber, "Structure #%ld", HGen->num_inp);
        ret2a = WriteOrigAtomDataToSDfile( orig_inp_data, output_file, szNumber, NULL,
                (sd->bChiralFlag & FLAG_INP_AT_CHIRAL)? 1:0,
                (ip->bINChIOutputOptions & INCHI_OUT_SDFILE_ATOMS_DT)? 1:0, ip->pSdfLabel, ip->pSdfValue );
        goto exit_function;
    }

    /******* create full reversibility information **************/
    if ( !(ip->bINChIOutputOptions & (INCHI_OUT_NO_AUX_INFO | INCHI_OUT_SHORT_AUX_INFO)) ) 
    {
        pOrigStruct = &(HGen->OrigStruct);
        memset( pOrigStruct, 0, sizeof(*pOrigStruct));
        if ( FillOutOrigStruct( orig_inp_data, pOrigStruct, sd ) ) 
        {
            AddMOLfileError(sd->pStrErrStruct, "Cannot interpret reversibility information");
            sd->nStructReadError =  99;
            sd->nErrorType = _IS_ERROR;
            nRet = _IS_ERROR;
        }
    }


    sd->bUserQuit  = 0;
    if (sd->bUserQuit) goto exit_function;



    /*^^^ Normalize the whole disconnected or original structure */
    if ( nRet != _IS_FATAL && nRet != _IS_ERROR ) {
        nRet1 = NormOneStructureINChI( pGenData, HGen, INCHI_BAS, inp_file);
        nRet = inchi_max(nRet, nRet1);
    }
    
/*
    if ( nRet != _IS_FATAL && nRet != _IS_ERROR )
        maxINChI = 1;
*/

    if ( nRet != _IS_FATAL && nRet != _IS_ERROR &&
        (sd->bTautFlagsDone[INCHI_BAS] & TG_FLAG_DISCONNECT_COORD_DONE) &&
        (ip->bTautFlags               & TG_FLAG_RECONNECT_COORD)          ) 
    {
        /* Normalize  the whole reconnected structure */
        nRet1 = NormOneStructureINChI( pGenData, HGen, INCHI_REC, inp_file);
        nRet = inchi_max(nRet, nRet1);
/*
        if ( nRet != _IS_FATAL && nRet != _IS_ERROR )
                maxINChI = 2;
*/
    }


    

exit_function:


    if ( nRet != _IS_FATAL && nRet != _IS_ERROR )
        HGen->norm_passed = 1;

    
    for (k=0; k < INCHI_NUM; k++)  
        pGenData->num_components[k] = sd->num_components[k];


    /*^^^ issue normalization warnings */
    if ( nRet != _IS_FATAL && nRet != _IS_ERROR )
    {
        int ic, istruct, itaut, nc[2];
        int warn_prot=0, warn_neutr=0;
        INP_ATOM_DATA *inp_norm_data[TAUT_NUM]; /*  = { &InpNormAtData, &InpNormTautData }; */
        nc[0] = pGenData->num_components[0];
        nc[1] = pGenData->num_components[1];
        for (istruct=0; istruct<2; istruct++)
        {
            if (nc[istruct]>0)
            {
                for (ic=0; ic < nc[istruct]; ic++)
                {
                    inp_norm_data[0] = &(HGen->InpNormAtData[istruct][ic]);
                    inp_norm_data[1] = &(HGen->InpNormTautData[istruct][ic]);
                    for (itaut=0;itaut<2;itaut++)
                    {
                        if (NULL!=inp_norm_data[itaut])
                        {
                            if ( inp_norm_data[itaut]->bTautomeric )
                            {
                                if (inp_norm_data[itaut]->bNormalizationFlags & (FLAG_NORM_CONSIDER_TAUT &~FLAG_PROTON_CHARGE_CANCEL) ) 
                                if (warn_prot==0)
                                {
                                    warn_prot++;
                                    AddMOLfileError(sd->pStrErrStruct, "Proton(s) added/removed");
                                }
                                if (inp_norm_data[itaut]->bNormalizationFlags & FLAG_PROTON_CHARGE_CANCEL )
                                if (warn_neutr==0)
                                {
                                    warn_neutr++;
                                    AddMOLfileError(sd->pStrErrStruct, "Charges neutralized");
                                }

                            }
                        }
                    } /* itaut */
                
                }
            }
        }

    }
    strcpy(pGenData->pStrErrStruct, sd->pStrErrStruct);
    make_norm_atoms_from_inp_atoms(pGenData, HGen);
    
    return nRet;
}




/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Get canonicalized form of the structure.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
int INCHI_DECL INCHIGEN_DoCanonicalization(void * _HGen, INCHIGEN_DATA *pGenData )
{
int nRet = 0, nRet1 /*, maxINChI=0*/;

FILE *inp_file = NULL;

INCHIGEN_CONTROL * HGen = _HGen;


STRUCT_DATA *sd = &(HGen->StructData);
INPUT_PARMS *ip = &(HGen->InpParms);
INCHI_FILE *output_file = HGen->inchi_file, *log_file = HGen->inchi_file+1;
FILE *prb_file = NULL;    
ORIG_ATOM_DATA *prep_inp_data = &(HGen->PrepInpData[0]);
int k;


/*^^^ Set debug output */
#if (TRACE_MEMORY_LEAKS == 1)

    _CrtSetDbgFlag(_CRTDBG_CHECK_ALWAYS_DF | _CRTDBG_LEAK_CHECK_DF | _CRTDBG_ALLOC_MEM_DF);
    /* for execution outside the VC++ debugger uncomment one of the following two */
#ifdef MY_REPORT_FILE 
   _CrtSetReportMode( _CRT_WARN, _CRTDBG_MODE_FILE );
   _CrtSetReportFile( _CRT_WARN, MY_REPORT_FILE );
   _CrtSetReportMode( _CRT_ERROR, _CRTDBG_MODE_FILE );
   _CrtSetReportFile( _CRT_ERROR, MY_REPORT_FILE );
   _CrtSetReportMode( _CRT_ASSERT, _CRTDBG_MODE_FILE );
   _CrtSetReportFile( _CRT_ASSERT, MY_REPORT_FILE );
#else
    _CrtSetReportMode(_CRT_WARN | _CRT_ERROR, _CRTDBG_MODE_DEBUG);
#endif 
#if ( !defined(__STDC__) || __STDC__ != 1 )
    /* turn on floating point exceptions */
    {
        /* Get the default control word. */
        int cw = _controlfp( 0,0 );
        /* Set the exception masks OFF, turn exceptions on. */
        /*cw &=~(EM_OVERFLOW|EM_UNDERFLOW|EM_INEXACT|EM_ZERODIVIDE|EM_DENORMAL);*/
        cw &=~(EM_OVERFLOW|EM_UNDERFLOW|EM_ZERODIVIDE|EM_DENORMAL);
        /* Set the control word. */
        _controlfp( cw, MCW_EM );
    }
#endif /*^^^ ( !defined(__STDC__) || __STDC__ != 1 ) */

#endif /*^^^ (TRACE_MEMORY_LEAKS == 1) */



    if (HGen->norm_passed==0)
    {
        AddMOLfileError(sd->pStrErrStruct, "Got non-normalized structure");
        sd->nStructReadError =  99;
        sd->nErrorType = _IS_ERROR;
        nRet = _IS_ERROR;
        goto exit_function;
    }

    sd->bUserQuit  = 0;
    if (sd->bUserQuit) goto exit_function;

   /* create INChI for each connected component of the structure and optionally display them */
   /* output INChI for the whole structure */
    


    /* create INChI for each connected component of the structure and optionally display them */
    /* create INChI for the whole disconnected or original structure */
    if ( nRet != _IS_FATAL && nRet != _IS_ERROR ) {
        nRet1 = CanonOneStructureINChI(pGenData, HGen, INCHI_BAS,inp_file);
        nRet = inchi_max(nRet, nRet1);
    }
    
/*
    if ( nRet != _IS_FATAL && nRet != _IS_ERROR )
        maxINChI = 1;
*/

    if ( nRet != _IS_FATAL && nRet != _IS_ERROR &&
        (sd->bTautFlagsDone[INCHI_BAS] & TG_FLAG_DISCONNECT_COORD_DONE) &&
        (ip->bTautFlags               & TG_FLAG_RECONNECT_COORD)          ) 
    {
        /* create INChI for the whole reconnected structure */
        nRet1 = CanonOneStructureINChI( pGenData, HGen,  INCHI_REC,
                                       inp_file);
        nRet = inchi_max(nRet, nRet1);
/*
        if ( nRet != _IS_FATAL && nRet != _IS_ERROR )
                maxINChI = 2;
*/
    }

    if (nRet != _IS_FATAL && nRet != _IS_ERROR) 
    {
        if (  (sd->bChiralFlag & FLAG_INP_AT_CHIRAL) &&
              (ip->nMode & REQ_MODE_STEREO) &&
              !(ip->nMode & (REQ_MODE_RELATIVE_STEREO | REQ_MODE_RACEMIC_STEREO)) &&
              !bIsStructChiral( HGen->pINChI, sd->num_components ) ) 
        {
            AddMOLfileError(sd->pStrErrStruct, "Not chiral");
        }
        
        /*************************************/
        /*       Output err/warn messages    */
        /*************************************/
        if ( /*!sd->nErrorCode &&*/ !sd->bUserQuitComponent && !sd->bUserQuit ) 
        {
            /*  if successful then returns 0, otherwise returns _IS_FATAL */
            /*  extract the structure if requested */
            nRet1 = TreatCreateINChIWarning(sd, ip, prep_inp_data, HGen->num_inp,
                                 inp_file, log_file, output_file, prb_file, HGen->pStr, PSTR_BUFFER_SIZE);
            nRet = inchi_max(nRet, nRet1);
        }
    }
    



    switch (nRet) 
    { 
        case _IS_SKIP   : nRet = inchi_Ret_SKIP   ; break; /* not used in INChI dll */
        case _IS_EOF    : nRet = inchi_Ret_EOF    ; break; /* no structural data has been provided */
        case _IS_OKAY   : nRet = inchi_Ret_OKAY   ; HGen->canon_passed = 1; break; 
                                                    /* Success; break; no errors or warnings */
        case _IS_WARNING: nRet = inchi_Ret_WARNING; HGen->canon_passed = 1; break; 
                                                    /* Success; break; warning(s) issued */
        case _IS_ERROR  : nRet = inchi_Ret_ERROR  ; break; /* Error: no INChI has been created */
        case _IS_FATAL  : nRet = inchi_Ret_FATAL  ; break; /* Severe error: no INChI has been created (typically; break; memory allocation failed) */
        case _IS_UNKNOWN:
        default         : nRet = inchi_Ret_UNKNOWN; break; /* Unknown program error */
    }
exit_function:    

    strcpy(pGenData->pStrErrStruct, sd->pStrErrStruct);
    for (k=0; k < INCHI_NUM; k++)  
        pGenData->num_components[k] = sd->num_components[k];


    return nRet;

} /*^^^ INCHIGEN_DoCanonicalization */





/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Get serialized form (InChI string).
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
int INCHI_DECL INCHIGEN_DoSerialization(void * _HGen,
                                  INCHIGEN_DATA * pGenData, 
                                  inchi_Output * pResults)
{
int nRet=0, nRet1=0, i, k;


FILE *inp_file = NULL;
FILE *prb_file=NULL;

INCHIGEN_CONTROL * HGen = _HGen;

INPUT_PARMS *ip = &(HGen->InpParms);
INCHI_FILE *output_file = HGen->inchi_file, *log_file = HGen->inchi_file+1;    
STRUCT_DATA *sd = &(HGen->StructData);
NORM_CANON_FLAGS *pncFlags = &(HGen->ncFlags);
ORIG_ATOM_DATA *orig_inp_data = &(HGen->OrigInpData);
ORIG_ATOM_DATA *prep_inp_data = &(HGen->PrepInpData[0]);
ORIG_STRUCT      *pOrigStruct = &(HGen->OrigStruct);
int bSortPrintINChIFlags   = 0;

int retcode = 0;

        pResults->szLog = log_file->pStr;

/*^^^ Set debug output */

#if (TRACE_MEMORY_LEAKS == 1)

    _CrtSetDbgFlag(_CRTDBG_CHECK_ALWAYS_DF | _CRTDBG_LEAK_CHECK_DF | _CRTDBG_ALLOC_MEM_DF);

/* for execution outside the VC++ debugger uncomment one of the following two */
#ifdef MY_REPORT_FILE 
   _CrtSetReportMode( _CRT_WARN, _CRTDBG_MODE_FILE );
   _CrtSetReportFile( _CRT_WARN, MY_REPORT_FILE );
   _CrtSetReportMode( _CRT_ERROR, _CRTDBG_MODE_FILE );
   _CrtSetReportFile( _CRT_ERROR, MY_REPORT_FILE );
   _CrtSetReportMode( _CRT_ASSERT, _CRTDBG_MODE_FILE );
   _CrtSetReportFile( _CRT_ASSERT, MY_REPORT_FILE );
#else
    _CrtSetReportMode(_CRT_WARN | _CRT_ERROR, _CRTDBG_MODE_DEBUG);
#endif 
#if ( !defined(__STDC__) || __STDC__ != 1 )
    /* turn on floating point exceptions */
    {
        /* Get the default control word. */
        int cw = _controlfp( 0,0 );
        /* Set the exception masks OFF, turn exceptions on. */
        /*cw &=~(EM_OVERFLOW|EM_UNDERFLOW|EM_INEXACT|EM_ZERODIVIDE|EM_DENORMAL);*/
        cw &=~(EM_OVERFLOW|EM_UNDERFLOW|EM_ZERODIVIDE|EM_DENORMAL);
        /* Set the control word. */
        _controlfp( cw, MCW_EM );
    }
#endif /*^^^ ( !defined(__STDC__) || __STDC__ != 1 ) */

#endif /*^^^ (TRACE_MEMORY_LEAKS == 1) */


/*****************************/


    if (HGen->canon_passed==0)
    {
        AddMOLfileError(sd->pStrErrStruct, "Got non-canonicalized structure");
        sd->nStructReadError =  99;
        sd->nErrorType = _IS_ERROR;
        retcode = _IS_ERROR;
        goto frees;
    }


    /************************************************/
    /*  sort and print INChI for the whole structure */
    /************************************************/
    nRet = SortAndPrintINChI(output_file, HGen->pStr, PSTR_BUFFER_SIZE, log_file, ip, orig_inp_data, prep_inp_data,
                                 HGen->composite_norm_data, pOrigStruct,
                                 sd->num_components, sd->num_non_taut, sd->num_taut,
                                 sd->bTautFlags, sd->bTautFlagsDone, pncFlags, HGen->num_inp,
                                 HGen->pINChI, HGen->pINChI_Aux, &bSortPrintINChIFlags );


    /* XML struct end tag */
    if ( (ip->bINChIOutputOptions & INCHI_OUT_XML) && sd->bXmlStructStarted > 0 ) 
    {
        if ( !OutputINChIXmlStructEndTag( output_file, HGen->pStr, PSTR_BUFFER_SIZE, 1 ) ) 
        {
            my_fprintf( log_file, "Cannot create end xml tag for structure #%ld.%s%s%s%s Terminating.\n", 
                            HGen->num_inp, SDF_LBL_VAL(ip->pSdfLabel,ip->pSdfValue) );
            sd->bXmlStructStarted = -1; /*  do not repeat same message */
            nRet = _IS_FATAL;
        } 
        else 
        {
            sd->bXmlStructStarted = 0; /*  do not continue xml output for this structure */
        }
    }
    if ( nRet != _IS_FATAL && nRet != _IS_ERROR ) 
    {
        /* Special mode: extract all good MOLfiles into the problem file
        * Do not extract any MOLfile that could not be processed (option /PGO)
        */
        if ( prb_file && 0L <= sd->fPtrStart && sd->fPtrStart < sd->fPtrEnd && ip->bSaveAllGoodStructsAsProblem ) {
            CopyMOLfile(inp_file, sd->fPtrStart, sd->fPtrEnd, prb_file, 0);
        }
#if( /*bRELEASE_VERSION != 1 &&*/ EXTR_FLAGS == EXTR_TRANSPOSITION_EXAMPLES && EXTR_MASK == EXTR_FLAGS )
        else
        if ( prb_file && (bSortPrintINChIFlags & 
             ( FLAG_SORT_PRINT_TRANSPOS_BAS | FLAG_SORT_PRINT_TRANSPOS_REC ) )
            ) 
        {
            CopyMOLfile(inp_file, sd->fPtrStart, sd->fPtrEnd, prb_file, 0);
        }
#endif
    }
    
    for ( i = 0; i < INCHI_NUM; i ++ ) 
    {
        for ( k = 0; k < TAUT_NUM+1; k ++ ) 
        {
            FreeCompAtomData( &(HGen->composite_norm_data[i][k]) );
        }
    }
/*****************************/        

    /*^^^ Prepare output message(s). */

    
    /*^^^ Error/warning. */
    if ( sd->pStrErrStruct[0] ) 
        if ( pGenData && (pResults->szMessage = (char *)inchi_malloc( strlen(sd->pStrErrStruct) + 1 )) ) 
            strcpy( pResults->szMessage, sd->pStrErrStruct );
    
    /*^^^ InChI, AuxInfo  (go to  pResults->szInChI, pResults->szAuxInfo) */
    if ( output_file->pStr && output_file->nUsedLength > 0 && pGenData ) 
    {
        char *p;
        pResults->szInChI   = output_file->pStr;
        pResults->szAuxInfo = NULL;
        for ( p = strchr(pResults->szInChI, '\n'); p; p = strchr(p+1, '\n') ) 
        {
            if ( !memcmp( p, "\nAuxInfo", 8 ) ) 
            {
                *p = '\0';            /* remove LF after INChI */
                pResults->szAuxInfo = p+1; /* save pointer to AuxInfo */
            } 
            else
                if ( pResults->szAuxInfo || !p[1]) 
                {   
                    /* remove LF after aux info or from the last char */
                    *p = '\0';
                    break;
                }
        }
        output_file->pStr = NULL;
    }


    /*^^^ Log message. */
    if ( log_file->pStr && log_file->nUsedLength > 0 ) 
    {
        while ( log_file->nUsedLength && '\n' == log_file->pStr[log_file->nUsedLength-1] ) 
            log_file->pStr[-- log_file->nUsedLength]  = '\0'; /* remove last LF */
        if ( pGenData ) 
        {
            pResults->szLog = log_file->pStr;
            log_file->pStr = NULL;
        }
    }
    
    if ( output_file->pStr )    {inchi_free( output_file->pStr ); output_file->pStr  = NULL;}
    if ( log_file->pStr    )    {inchi_free( log_file->pStr ); log_file->pStr = NULL;}
    

        
        

    HGen->ulTotalProcessingTime += sd->ulStructTime;
    nRet = inchi_max(nRet, nRet1);
    
    switch ( nRet ) 
    {
        case _IS_FATAL:     
        case _IS_ERROR:     HGen->num_err ++;
    }



    /*^^^ XML-related. */
    if ( (ip->bINChIOutputOptions & INCHI_OUT_XML) && sd->bXmlStructStarted > 0 ) 
    {
        if ( !OutputINChIXmlStructEndTag( output_file, HGen->pStr, PSTR_BUFFER_SIZE, 1 ) ) 
        {
            my_fprintf( log_file, "Cannot create end xml tag for structure #%d.%s%s%s%s Terminating.\n", 
                HGen->num_inp, SDF_LBL_VAL(ip->pSdfLabel,ip->pSdfValue) );
            sd->bXmlStructStarted = -1; /*  do not repeat same message */
        }
    }
    if ( (ip->bINChIOutputOptions & INCHI_OUT_XML) && ip->bXmlStarted ) 
    {
        OutputINChIXmlRootEndTag( output_file );
        ip->bXmlStarted = 0;
    }


frees:    
    /*^^^ Free all. */

    /*  free INChI memory */
/*     FreeAllINChIArrays(HGen->pINChI, HGen->pINChI_Aux, sd->num_components ); 
*/
#if( ADD_CMLPP == 1 )
        /* BILLY 8/6/04 */
        /* free CML memory */
        FreeCml ();
        FreeCmlDoc( 1 );
#endif
    for ( i = 0; i < MAX_NUM_PATHS; i ++ ) 
    {
        if ( ip->path[i] ) 
        {
            inchi_free( (char*) ip->path[i] ); /*  cast deliberately discards 'const' qualifier */
            ip->path[i] = NULL;
        }
    }
    SetBitFree( );


    strcpy(pGenData->pStrErrStruct, sd->pStrErrStruct);
    for (k=0; k < INCHI_NUM; k++)  
        pGenData->num_components[k] = sd->num_components[k];


    return retcode;
}


/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

InChI Generator: reset stage (use before get next structure)

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
void INCHI_DECL INCHIGEN_Reset(void * _HGen,
                                 INCHIGEN_DATA * pGenData, 
                                 inchi_Output * pResults)
{
int i, k, nc;
INCHIGEN_CONTROL * HGen = _HGen;


    if ( pResults->szInChI )    inchi_free( pResults->szInChI );
    if ( pResults->szLog )      inchi_free( pResults->szLog );
    if ( pResults->szMessage )  inchi_free( pResults->szMessage );
    

    /* Free all data associated with components of disconn/conn structures */

    if (NULL!=HGen)
    {

        memset(HGen->inchi_file, 0, sizeof(HGen->inchi_file));

        if (HGen->pStr )            
            memset(HGen->pStr, 0, sizeof( (*HGen->pStr) ));


        for ( i = 0; i < MAX_NUM_PATHS; i ++ ) 
        {
            if ( HGen->InpParms.path[i] ) 
            {
                inchi_free( (char*) HGen->InpParms.path[i] ); /*  cast deliberately discards 'const' qualifier */
                HGen->InpParms.path[i] = NULL;
            }
        }
        memset(&(HGen->InpParms), 0, sizeof(INPUT_PARMS));

        FreeOrigAtData( &(HGen->OrigInpData) );
        memset(&(HGen->OrigInpData), 0, sizeof( HGen->OrigInpData ) );
    
    
        FreeOrigAtData( &(HGen->PrepInpData[0]));
        FreeOrigAtData( &(HGen->PrepInpData[1]));
        memset(&(HGen->PrepInpData[0]), 0, 2*sizeof( HGen->PrepInpData[0] ) );  

        FreeOrigStruct( &(HGen->OrigStruct));
        memset(&(HGen->OrigStruct), 0, sizeof( HGen->OrigStruct ) );


        for ( i = 0; i < INCHI_NUM; i ++ ) 
            for ( k = 0; k < TAUT_NUM+1; k ++ ) 
                FreeCompAtomData( &(HGen->composite_norm_data[i][k]) );

        for ( k = 0;  k < INCHI_NUM; k++) 
        {
            nc = HGen->StructData.num_components[k];

            if ( HGen->InpCurAtData[k] )
            {
                for ( i = 0;  i < nc; i ++ ) 
                    FreeInpAtomData( &(HGen->InpCurAtData[k][i]) );
                inchi_free(HGen->InpCurAtData[k]);
                HGen->InpCurAtData[k] = NULL;
            }


            if ( HGen->cti[k] )
            {
                if ( (HGen->cti[k])->at[TAUT_YES] )
                {
                    inchi_free( (HGen->cti[k])->at[TAUT_YES] );
                    (HGen->cti[k])->at[TAUT_YES] = NULL;
                }
    
                if ( (HGen->cti[k])->at[TAUT_NON] )
                {
                    inchi_free( (HGen->cti[k])->at[TAUT_NON] );
                    (HGen->cti[k])->at[TAUT_NON] = NULL;
                }

                if (&((HGen->cti[k])->vt_group_info)) 
                    free_t_group_info(&((HGen->cti[k])->vt_group_info));
    
                if (&((HGen->cti[k])->vt_group_info_orig)) 
                    free_t_group_info(&((HGen->cti[k])->vt_group_info_orig));

                inchi_free(HGen->cti[k]);
                HGen->cti[k] = NULL;
            }

        }


        for ( k = 0;  k < INCHI_NUM; k++) 
        {
            nc = HGen->StructData.num_components[k];
            if ( HGen->InpNormAtData[k] )
            {
                for ( i = 0;  i < nc; i ++ ) 
                    FreeInpAtomData( &(HGen->InpNormAtData[k][i]) );
                inchi_free(HGen->InpNormAtData[k]);
                HGen->InpNormAtData[k] = NULL;
            }

            if ( HGen->InpNormTautData[k] )
            {
                for ( i = 0;  i < nc; i ++ ) 
                    FreeInpAtomData( &(HGen->InpNormTautData[k][i]) );
                inchi_free(HGen->InpNormTautData[k]);
                HGen->InpNormTautData[k] = NULL;
            }


            if ( pGenData->NormAtomsTaut[k] )
            {
                /* 
                for ( i = 0;  i < nc; i ++ ) 
                    FreeInpAtomData( &(pGenData->NormAtomsTaut[k][i]) );
                */
                inchi_free(pGenData->NormAtomsTaut[k]);
                pGenData->NormAtomsTaut[k] = NULL;
            }
            if ( pGenData->NormAtomsNontaut[k] )
            {
                /* for ( i = 0;  i < nc; i ++ ) 
                    FreeInpAtomData( &(pGenData->NormAtomsNontaut[k][i]) );
                */
                inchi_free(pGenData->NormAtomsNontaut[k]);
                pGenData->NormAtomsNontaut[k] = NULL;
            }


        }

        /*  free INChI memory */
        FreeAllINChIArrays( HGen->pINChI, HGen->pINChI_Aux, HGen->StructData.num_components );
        memset(HGen->pINChI,     0, sizeof(HGen->pINChI) );
        memset(HGen->pINChI_Aux, 0, sizeof(HGen->pINChI_Aux) );
        
    }

    memset(&(HGen->StructData), 0, sizeof(STRUCT_DATA) );
    HGen->szTitle[0] = '\0';

    memset( pResults, 0, sizeof(*pResults) );
    memset( pGenData , 0, sizeof(*pGenData) );

    return;
}


/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
InChI Generator: destroy generator
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
EXPIMP_TEMPLATE INCHI_API 
    void INCHI_DECL INCHIGEN_Destroy(void * _HGen)
{
INCHIGEN_CONTROL * HGen = _HGen;

    if (NULL!=HGen)
    {

        if ( HGen->pStr )           
            inchi_free(HGen->pStr);

        inchi_free(HGen);

    }


}





/********************************************************************/

#if( defined( _WIN32 ) && defined( _MSC_VER ) && _MSC_VER >= 800 && defined(_USRDLL) && defined(INCHI_LINK_AS_DLL) )
    /* Win32 & MS VC ++, compile and link as a DLL */
/*********************************************************/
/*   C calling conventions export from Win32 dll         */
/*********************************************************/
/* prototypes */
#ifndef INCHI_ALL_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif
void * cdecl_INCHIGEN_Create(void);
int  cdecl_INCHIGEN_Setup( void * HGen, INCHIGEN_DATA * pGenData, inchi_Input * pInp);
int  cdecl_INCHIGEN_DoNormalization( void * HGen, INCHIGEN_DATA *pGenData );
int  cdecl_INCHIGEN_DoCanonicalization( void * HGen, INCHIGEN_DATA *pGenData );
int  cdecl_INCHIGEN_DoSerialization(void * HGen, INCHIGEN_DATA * pGenData, inchi_Output * pResults );
void cdecl_INCHIGEN_Reset( void * HGen, INCHIGEN_DATA * pGenData, inchi_Output * pResults);
void cdecl_INCHIGEN_dESTROY( void * HGen );
#ifndef INCHI_ALL_CPP
#ifdef __cplusplus
}
#endif
#endif

/* implementation */
/* vc6_libinchi.def provides export withou cdecl_ prefixes */

/********************************************************/
void * cdecl_INCHIGEN_Create(void)
{
    return INCHIGEN_Create( );
}
/********************************************************/
int cdecl_INCHIGEN_Setup( void * HGen, INCHIGEN_DATA * pGenData, inchi_Input * pInp)
{
    return INCHIGEN_Setup( HGen, pGenData, pInp );
}
/********************************************************/
int cdecl_INCHIGEN_DoNormalization( void * HGen, INCHIGEN_DATA *pGenData )
{
    return INCHIGEN_DoNormalization( HGen, pGenData );
}
/********************************************************/
int cdecl_INCHIGEN_DoCanonicalization( void * HGen, INCHIGEN_DATA *pGenData )
{
    return INCHIGEN_DoCanonicalization( HGen, pGenData );
}
/********************************************************/
int cdecl_INCHIGEN_DoSerialization( void * HGen, INCHIGEN_DATA * pGenData, inchi_Output * pResults )
{
    return INCHIGEN_DoSerialization( HGen, pGenData, pResults );
}
/********************************************************/
void cdecl_INCHIGEN_Reset( void * HGen , INCHIGEN_DATA *pGenData, inchi_Output *pResults)
{
    INCHIGEN_Reset( HGen, pGenData, pResults );
}
/********************************************************/
void cdecl_INCHIGEN_Destroy( void * HGen )
{
    INCHIGEN_Destroy( HGen );
}

#endif

#if( defined(__GNUC__) && __GNUC__ >= 3 && defined(__MINGW32__) && defined(_WIN32) )
#include <windows.h>
/*********************************************************/
/*   Pacal calling conventions export from Win32 dll     */
/*********************************************************/
#ifndef INCHI_ALL_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif
/* prototypes */
/********************************************************/
void * PASCAL pasc_INCHIGEN_Create(void);
int  PASCAL pasc_INCHIGEN_Setup( void * HGen, INCHIGEN_DATA * pGenData, inchi_Input * pInp );
int  PASCAL pasc_INCHIGEN_DoNormalization( void * HGen, INCHIGEN_DATA *pGenData );
int  PASCAL pasc_INCHIGEN_DoCanonicalization( void * _HGen, INCHIGEN_DATA *pGenData );
int  PASCAL pasc_INCHIGEN_DoSerialization( void * HGen, INCHIGEN_DATA * pGenData, inchi_Output * pResults );
void PASCAL pasc_INCHIGEN_Reset( void * HGen, INCHIGEN_DATA *pGenData, inchi_Output *pResults );
void PASCAL pasc_INCHIGEN_Destroy(void * HGen);
#ifndef INCHI_ALL_CPP
#ifdef __cplusplus
}
#endif
#endif

/* implementation */
/* vc6_libinchi.def provides export without PASCAL pasc_ prefixes */
/********************************************************/
void * PASCAL pasc_INCHIGEN_Create(void)
{
    return INCHIGEN_Create( );
}
/********************************************************/
int PASCAL pasc_INCHIGEN_Setup( void * HGen, INCHIGEN_DATA * pGenData, inchi_Input * pInp)
{
    return INCHIGEN_Setup( HGen, pGenData, pInp);
}
/********************************************************/
int PASCAL pasc_INCHIGEN_DoNormalization( void * HGen, INCHIGEN_DATA *pGenData )
{
    return INCHIGEN_DoNormalization( HGen, pGenData );
}
/********************************************************/
int PASCAL pasc_INCHIGEN_DoCanonicalization( void * HGen, INCHIGEN_DATA *pGenData )
{
    return INCHIGEN_DoCanonicalization( HGen, pGenData );
}
/********************************************************/
int PASCAL pasc_INCHIGEN_DoSerialization( void * HGen, INCHIGEN_DATA * pGenData, inchi_Output * pResults )
{
    return INCHIGEN_DoSerialization( HGen, pGenData, pResults );
}
/********************************************************/
void PASCAL pasc_INCHIGEN_Reset( void * HGen, INCHIGEN_DATA *pGenData,  inchi_Output *pResults)
{
    INCHIGEN_Reset( HGen, pGenData, pResults );
}
/********************************************************/
void PASCAL pasc_INCHIGEN_Destroy(void * HGen)
{
    INCHIGEN_Destroy( HGen);
}
#endif 










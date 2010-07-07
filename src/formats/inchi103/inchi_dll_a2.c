/*
 * International Chemical Identifier (InChI)
 * Version 1
 * Software version 1.03
 * May 9, 2010
 *
 * Originally developed at NIST
 * Modifications and additions by IUPAC and the InChI Trust
 *
 * The InChI library and programs are free software developed under the
 * auspices of the International Union of Pure and Applied Chemistry (IUPAC);
 * you can redistribute this software and/or modify it under the terms of 
 * the GNU Lesser General Public License as published by the Free Software 
 * Foundation:
 * http://www.opensource.org/licenses/lgpl-2.1.php
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdarg.h>
/* #include <varargs.h> */
#include <errno.h>
#include <limits.h>


/*^^^  */

#include "mode.h"      

#include "ichitime.h"


#include "inpdef.h"
#include "ichi.h"
#include "strutil.h"
#include "util.h"
#include "ichidrp.h"
#include "ichierr.h"
#include "ichimain.h"
#include "extr_ct.h"
#include "ichitaut.h"
#include "ichi_io.h"

#include "ichinorm.h"
#include "ichicant.h"
#include "ichicano.h"
#include "ichicomn.h"
#include "ichimake.h"
#include "ichister.h"
/*^^^ */


#ifdef INCHI_LIB
#include "ichi_lib.h"
#endif

#include "ichicomp.h"

#if( ADD_CMLPP == 1 )
#include "readcml.hpp"
#include "debug.h"
#endif
/*^^^  */

/* for DisplayTheWholeStructure() */

#define COMP_ORIG_0_MAIN  0x0001
#define COMP_ORIG_0_RECN  0x0002
#define COMP_PREP_0_MAIN  0x0004
#define COMP_PREP_0_RECN  0x0008
#define COMP_ORIG_1_MAIN  0x0010
#define COMP_ORIG_1_RECN  0x0020

#include "ichisize.h"
#include "mode.h"
#include "inchi_api.h"
#include "inchi_dll_a.h" /* not inchi_api.h as it hides internal data types */


int  GetProcessingWarnings(INChI *cur_INChI[], INP_ATOM_DATA **inp_norm_data, STRUCT_DATA *sd);

/*^^^ */


int inp2spATOM( inp_ATOM *inp_at, int num_inp_at, sp_ATOM *at );
int CheckCanonNumberingCorrectness(
                 int num_atoms, int num_at_tg,
                 sp_ATOM *at, CANON_STAT *pCS, int bTautomeric,
                 char *pStrErrStruct );

int CreateCompositeNormAtom( COMP_ATOM_DATA *composite_norm_data, 
                             INP_ATOM_DATA2 *all_inp_norm_data,
                             int num_components);

int CopyLinearCTStereoToINChIStereo( INChI_Stereo *Stereo,
           AT_STEREO_CARB *LinearCTStereoCarb, int nLenLinearCTStereoCarb,
           AT_STEREO_DBLE *LinearCTStereoDble, int nLenLinearCTStereoDble
           , AT_NUMB *pCanonOrd, AT_RANK *pCanonRank, sp_ATOM *at, int bIsotopic
           , AT_STEREO_CARB *LinearCTStereoCarbInv
           , AT_STEREO_DBLE *LinearCTStereoDbleInv
           , AT_NUMB *pCanonOrdInv, AT_RANK *pCanonRankInv );
int MarkAmbiguousStereo( sp_ATOM *at, inp_ATOM *norm_at, int bIsotopic, AT_NUMB *pCanonOrd,
           AT_STEREO_CARB *LinearCTStereoCarb, int nLenLinearCTStereoCarb,
           AT_STEREO_DBLE *LinearCTStereoDble, int nLenLinearCTStereoDble );

INCHI_MODE UnmarkAllUndefinedUnknownStereo( INChI_Stereo *Stereo, INCHI_MODE nUserMode );
int FillOutINChIReducedWarn( INChI *pINChI, INChI_Aux *pINChI_Aux,
                 int num_atoms, int num_at_tg, int num_removed_H,
                 sp_ATOM *at, inp_ATOM *norm_at, CANON_STAT *pCS, int bTautomeric,
                 INCHI_MODE nUserMode, char *pStrErrStruct );

int NormOneStructureINChI(INCHIGEN_DATA *gendata, INCHIGEN_CONTROL *genctl, 
                          int iINChI, INCHI_IOSTREAM *inp_file);
int CanonOneStructureINChI(INCHIGEN_CONTROL *genctl, int iINChI, INCHI_IOSTREAM *inp_file);

int NormOneComponentINChI(INCHIGEN_CONTROL *genctl, int iINChI, int i);
int CanonOneComponentINChI(INCHIGEN_CONTROL *genctl, int iINChI, int i);
int  Normalization_step( INChI **ppINChI, INChI_Aux **ppINChI_Aux, 
                        inp_ATOM *inp_at, INP_ATOM_DATA *out_norm_data[2],
                        int num_inp_at,
                        INCHI_MODE *pbTautFlags, 
                        INCHI_MODE *pbTautFlagsDone,
                        COMPONENT_TREAT_INFO *cti);
int  Canonicalization_step( INChI **ppINChI, INChI_Aux **ppINChI_Aux, 
                        INP_ATOM_DATA *out_norm_data[2],
                        struct tagInchiTime *ulMaxTime, 
                        T_GROUP_INFO *ti_out, 
                        char *pStrErrStruct,
                        COMPONENT_TREAT_INFO *cti);





/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
int NormOneStructureINChI(INCHIGEN_DATA *gendata, 
                          INCHIGEN_CONTROL *genctl, 
                          int iINChI,
                          INCHI_IOSTREAM *inp_file)
{
int k, i, j, nRet = 0;
int nStrLen = PSTR_BUFFER_SIZE;


STRUCT_DATA *sd = &(genctl->StructData);

INPUT_PARMS *ip = &(genctl->InpParms);
ORIG_ATOM_DATA *prep_inp_data = &(genctl->PrepInpData[0]);
ORIG_ATOM_DATA *orig_inp_data = &(genctl->OrigInpData);
INCHI_IOSTREAM *output_file = genctl->inchi_file, *log_file = genctl->inchi_file+1;   
INCHI_IOSTREAM prbstr, *prb_file=&prbstr;

PINChI2 **pINChI2 = genctl->pINChI; 
PINChI_Aux2 **pINChI_Aux2 = genctl->pINChI_Aux;
NORM_CANON_FLAGS *pncFlags = &(genctl->ncFlags);    
    

INP_ATOM_DATA *inp_cur_data = NULL;

long num_inp = genctl->num_inp;
char *pStr = genctl->pStr;
    
INP_ATOM_DATA *inp_norm_data[TAUT_NUM]; /*  = { &InpNormAtData, &InpNormTautData }; */
ORIG_ATOM_DATA *cur_prep_inp_data = prep_inp_data + iINChI;

inchiTime      ulTStart;

    /*^^^ To save intermediate data... */

    COMP_ATOM_DATA *composite_norm_data = genctl->composite_norm_data[iINChI];
    INP_ATOM_DATA2 *all_inp_norm_data = NULL;
    memset( composite_norm_data+TAUT_NON, 0, sizeof( composite_norm_data[0] ) );
    memset( composite_norm_data+TAUT_YES, 0, sizeof( composite_norm_data[0] ) );
    memset( composite_norm_data+TAUT_INI, 0, sizeof( composite_norm_data[0] ) );

    inchi_ios_init(prb_file, INCHI_IOSTREAM_FILE, NULL);

    /*
        if ( orig_inp_data is NOT empty AND
             prep_inp_data[0] IS empty ) then:

            1. copy orig_inp_data --> prep_inp_data[0]
            2. fix odd things in prep_inp_data[0]
            3. if( orig_inp_data->bDisconnectSalts ) then
                  -- disconnect salts in prep_inp_data[0]
            4. move protons to neutralize charges on heteroatoms
            5. if( orig_inp_data->bDisconnectCoord ) then
                  -- copy prep_inp_data[0] --> prep_inp_data[1]
                  -- disconnect metals in prep_inp_data[0]

            [ This all is done in PreprocessOneStructure() ]

        iINChI = 0
        =========
        (normal/disconnected layer)

            1. normalize prep_inp_data[0] in inp_norm_data[0,1]
            2. create INChI[ iINChI ] out of inp_norm_data[0,1]


        iINChI = 1 AND orig_inp_data->bDisconnectCoord > 0
        =================================================
        (reconnected layer)

            1. normalize prep_inp_data[1] in inp_norm_data[0,1]
            2. create INChI[ iINChI ] out of inp_norm_data[0,1]

    */



    ip->msec_LeftTime = ip->msec_MaxTime; /* start timeout countdown for each component */




    if ( ip->bAllowEmptyStructure && !orig_inp_data->at && !orig_inp_data->num_inp_atoms ) {
        ;
    } else
    if ( !orig_inp_data->at || !orig_inp_data->num_inp_atoms ) {
        return 0; /* nothing to do */
    }
    if ( iINChI == 1 && orig_inp_data->bDisconnectCoord <= 0 ) {
        return 0;
    }

   /* m = iINChI; */ /* orig_inp_data index */

    if ( iINChI != INCHI_BAS && iINChI != INCHI_REC ) {
        AddMOLfileError(sd->pStrErrStruct, "Fatal undetermined program error");
        sd->nStructReadError =  97;
        nRet = sd->nErrorType = _IS_FATAL;
        goto exit_function;
    }

    /*******************************************************************
     *                                                                 *
     *                                                                 *
     *  Whole structure preprocessing: 1st step of the normalization   *
     *                                                                 *
     *  Happen only on the first call to CreateOneStructureINChI()      *
     *                                                                 *
     *                                                                 *
     *******************************************************************/

    if ( (!prep_inp_data->at || !prep_inp_data->num_inp_atoms) && orig_inp_data->num_inp_atoms > 0 ) 
    {
        /* the structure has not been preprocessed */
        if ( ip->msec_MaxTime ) 
        {
            InchiTimeGet( &ulTStart );
        }

        PreprocessOneStructure( sd, ip, orig_inp_data, prep_inp_data );
        pncFlags->bTautFlags[iINChI][TAUT_YES] =
        pncFlags->bTautFlags[iINChI][TAUT_NON] = sd->bTautFlags[INCHI_BAS] | ip->bTautFlags;
        pncFlags->bTautFlagsDone[iINChI][TAUT_YES] =
        pncFlags->bTautFlagsDone[iINChI][TAUT_NON] = sd->bTautFlagsDone[INCHI_BAS] | ip->bTautFlagsDone;

        switch (sd->nErrorType) 
        {
        case _IS_ERROR:
        case _IS_FATAL:
            /* error message */
            nRet = TreatReadTheStructureErrors( sd, ip, LOG_MASK_ALL, inp_file, log_file, output_file, prb_file,
                                                prep_inp_data, &num_inp, pStr, nStrLen );
            goto exit_function;
        }
    }

    /*^^^ To save intermediate data... */
    if ( prep_inp_data[iINChI].num_components > 1) 
    {
        all_inp_norm_data = (INP_ATOM_DATA2 *)inchi_calloc( prep_inp_data[iINChI].num_components, sizeof(all_inp_norm_data[0]));
    }



    /* allocate pINChI[iINChI] and pINChI_Aux2[iINChI] -- arrays of pointers to INChI and INChI_Aux */
    /* assign values to sd->num_components[]                                                  */

    MYREALLOC2(PINChI2, PINChI_Aux2, pINChI2[iINChI], pINChI_Aux2[iINChI], sd->num_components[iINChI], cur_prep_inp_data->num_components, k);
    if ( k ) 
    {
        AddMOLfileError(sd->pStrErrStruct, "Cannot allocate output data. Terminating");
        sd->nStructReadError =  99;
        sd->nErrorType = _IS_FATAL;
        goto exit_function;
    }


    /*^^^ Allocate */

    /*^^^ visible */
    gendata->NormAtomsNontaut[iINChI] = (NORM_ATOMS *)inchi_calloc( sd->num_components[iINChI], sizeof(NORM_ATOMS));  
    gendata->NormAtomsTaut[iINChI] = (NORM_ATOMS *)inchi_calloc( sd->num_components[iINChI], sizeof(NORM_ATOMS));  
    /*^^^ invisible */
    genctl->InpNormAtData[iINChI] = (INP_ATOM_DATA *)inchi_calloc( sd->num_components[iINChI], sizeof(INP_ATOM_DATA)); 
    genctl->InpNormTautData[iINChI] = (INP_ATOM_DATA *)inchi_calloc( sd->num_components[iINChI], sizeof(INP_ATOM_DATA)); 
    genctl->InpCurAtData[iINChI] = (INP_ATOM_DATA *)inchi_calloc( sd->num_components[iINChI], sizeof(INP_ATOM_DATA)); 
    genctl->cti[iINChI] = (COMPONENT_TREAT_INFO *)inchi_calloc( sd->num_components[iINChI], sizeof(COMPONENT_TREAT_INFO)); 
    memset (genctl->cti[iINChI], 0, sd->num_components[iINChI]*sizeof(COMPONENT_TREAT_INFO));




    /*^^^ Second normalization step - component by component */


    for ( i = 0, nRet = 0; !sd->bUserQuitComponent && i < cur_prep_inp_data->num_components; i ++ ) 
    {

        if (ip->msec_MaxTime) InchiTimeGet( &ulTStart );

        inp_cur_data = &(genctl->InpCurAtData[iINChI][i]);
        
        /*  a) allocate memory and extract current component */
        nRet = GetOneComponent( sd, ip, log_file, output_file, inp_cur_data, 
                                cur_prep_inp_data, i, num_inp, pStr, nStrLen );
        
        if (ip->msec_MaxTime)   ip->msec_LeftTime -= InchiTimeElapsed( &ulTStart );
        
        switch (nRet)           
        { 
            case _IS_ERROR: 
            case _IS_FATAL: 
                goto exit_cycle; 
        }



        /*  c) Create the component's INChI ( copies ip->bTautFlags into sd->bTautFlags)*/

        inp_norm_data[TAUT_NON] = &(genctl->InpNormAtData[iINChI][i]);
        memset( inp_norm_data[TAUT_NON], 0, sizeof( *inp_norm_data[0] ) );
        inp_norm_data[TAUT_YES] = &(genctl->InpNormTautData[iINChI][i]);
        memset( inp_norm_data[TAUT_YES], 0, sizeof( *inp_norm_data[0] ) );
        

        nRet = NormOneComponentINChI(genctl, iINChI, i);

        /*^^^ To save intermediate data... */
        if ( all_inp_norm_data ) 
        {
            for ( j = 0; j < TAUT_NUM; j ++ ) 
            {
                if ( inp_norm_data[j]->bExists ) 
                {
                    all_inp_norm_data[i][j] = *inp_norm_data[j];
                    memset( inp_norm_data[j], 0, sizeof(*inp_norm_data[0]) );
                }
            }
        }


        if (nRet) 
        {
            nRet = TreatCreateOneComponentINChIError(sd, ip, cur_prep_inp_data, i, num_inp,
                                 inp_file, log_file, output_file, prb_file,pStr, nStrLen );
            break;
        }
    }


exit_cycle:

    switch ( nRet ) 
    {
    case _IS_FATAL:
    case _IS_ERROR: break;
    default:

        /*^^^ To save intermediate data... */
        if ( all_inp_norm_data ) 
        {
             int res = CreateCompositeNormAtom( composite_norm_data, 
                                                all_inp_norm_data, 
                                                prep_inp_data[iINChI].num_components);
        }
        break;
    }


    /*^^^ When saving intermediate data - avoid memory leaks in case of error */
    if ( all_inp_norm_data ) 
    {
        for ( i = 0; i < prep_inp_data[iINChI].num_components; i ++ ) 
        {
            for ( k = 0; k < TAUT_NUM; k ++ ) 
            {
                FreeInpAtomData( &all_inp_norm_data[i][k] );
            }
        }
        inchi_free( all_inp_norm_data );
        all_inp_norm_data = NULL;
    }







exit_function:

    return nRet;
}




/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
int CanonOneStructureINChI(INCHIGEN_CONTROL *genctl, int iINChI, INCHI_IOSTREAM *inp_file)
{
int i, /*m,*/ nRet = 0;

STRUCT_DATA *sd = &(genctl->StructData);

INPUT_PARMS *ip = &(genctl->InpParms);
INCHI_IOSTREAM *output_file = genctl->inchi_file, *log_file = genctl->inchi_file+1;
INCHI_IOSTREAM prbstr, *prb_file=&prbstr; 
 
ORIG_ATOM_DATA *prep_inp_data = &(genctl->PrepInpData[0]);
PINChI2 **pINChI2 = genctl->pINChI; 
PINChI_Aux2 **pINChI_Aux2 = genctl->pINChI_Aux;

long num_inp = genctl->num_inp;
char *pStr = genctl->pStr;
int nStrLen = PSTR_BUFFER_SIZE;



    INP_ATOM_DATA *inp_cur_data = NULL;

    INP_ATOM_DATA *inp_norm_data[TAUT_NUM]; /*  = { &InpNormAtData, &InpNormTautData }; */

    ORIG_ATOM_DATA *cur_prep_inp_data = prep_inp_data + iINChI;

    inchiTime      ulTStart;

    inchi_ios_init(prb_file, INCHI_IOSTREAM_FILE, NULL);

    for (i = 0; i < TAUT_NUM; i ++) /* initialize in case no InChI to generate 2008-12-23 DT */
        inp_norm_data[i]=NULL;

    /**************************************************************************/
    /*                                                                        */
    /*                                                                        */
    /*   M A I N   C Y C L E:   P R O C E S S    C O M P O N E N T S          */
    /*                                                                        */
    /*                                                                        */
    /*                     O N E   B Y   O N E                                */
    /*                                                                        */
    /*                                                                        */
    /**************************************************************************/

    for ( i = 0, nRet = 0; !sd->bUserQuitComponent && i < cur_prep_inp_data->num_components; i ++ ) 
    {

        if ( ip->msec_MaxTime ) 
            InchiTimeGet( &ulTStart );




        /*****************************************************/
        /*  a) allocate memory and extract current component */
        /*****************************************************/

        inp_cur_data = &(genctl->InpCurAtData[iINChI][i]);


        nRet = GetOneComponent( sd, ip, log_file, output_file, inp_cur_data, cur_prep_inp_data, 
                                i, num_inp, pStr, nStrLen );
        
        if ( ip->msec_MaxTime ) 
            ip->msec_LeftTime -= InchiTimeElapsed( &ulTStart );
        
        switch ( nRet ) { case _IS_ERROR: case _IS_FATAL: goto exit_cycle; }

#ifndef INCHI_LIBRARY
        /*  console request: Display the component? */
        if ( ip->bDisplay && inp_file != stdin ) 
        {
            if ( user_quit("Enter=Display Component, Esc=Stop ?", ip->ulDisplTime) ) 
            {
                sd->bUserQuitComponent = 1;
                break;
            }
        }
#endif




        /*******************************************************************************/
        /*                                                                             */
        /*                      C A N O N I C A L I Z A T I O N                        */
        /*                                                                             */
        /*         (both tautomeric and non-tautomeric if requested)                   */
        /*                                                                             */
        /*******************************************************************************/
        /*  c) Create the component's INChI ( copies ip->bTautFlags into sd->bTautFlags)*/
        /*******************************************************************************/

        inp_norm_data[TAUT_NON] = &(genctl->InpNormAtData[iINChI][i]);
        inp_norm_data[TAUT_YES] = &(genctl->InpNormTautData[iINChI][i]);

        nRet = CanonOneComponentINChI(genctl, iINChI, i);



        
        if ( nRet ) 
        {
            nRet = TreatCreateOneComponentINChIError(sd, ip, cur_prep_inp_data, i, num_inp,
                                 inp_file, log_file, output_file, prb_file,pStr, nStrLen );
            break;
        }
    }
    /**************************************************************************/
    /*                                                                        */
    /*                                                                        */
    /*   E N D   O F   T H E    M A I N   C Y C L E   P R O C E S S I N G     */
    /*                                                                        */
    /*          C O M P O N E N T S    O N E   B Y   O N E                    */
    /*                                                                        */
    /*                                                                        */
    /**************************************************************************/


exit_cycle:

    switch ( nRet ) 
    {
    case _IS_FATAL:
    case _IS_ERROR: break;
    default:

        break;
    }




    for (i = 0; i < TAUT_NUM; i ++) 
        FreeInpAtomData( inp_norm_data[i] );



    return nRet;
}




/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
int NormOneComponentINChI(INCHIGEN_CONTROL * genctl, int iINChI, int i)
{
STRUCT_DATA *sd = &(genctl->StructData);
INPUT_PARMS *ip = &(genctl->InpParms);
PINChI2 **pINChI2 = genctl->pINChI; 
PINChI_Aux2 **pINChI_Aux2 = genctl->pINChI_Aux;
NORM_CANON_FLAGS *pncFlags = &(genctl->ncFlags);

    inchiTime     ulTStart, ulTEnd, *pulTEnd = NULL;
    int           k, num_at, ret = 0;
    int           bOrigCoord;
    INCHI_MODE     bTautFlags     = ip->bTautFlags;
    INCHI_MODE     bTautFlagsDone = (ip->bTautFlagsDone | sd->bTautFlagsDone[INCHI_BAS]);
    long          lElapsedTime;
    /*
    PINChI2     *pINChI     = pINChI2[iINChI];
    PINChI_Aux2 *pINChI_Aux = pINChI_Aux2[iINChI];
    */

PINChI2     *pINChI     = NULL;
PINChI_Aux2 *pINChI_Aux = NULL;

INChI       *cur_INChI[TAUT_NUM];
INChI_Aux   *cur_INChI_Aux[TAUT_NUM];
    

/* pINChI2[m=iINChI-1][j< prep_inp_data[m].num_components][TAUT_NON] */

INP_ATOM_DATA *inp_norm_data[TAUT_NUM]; /*  = { &InpNormAtData, &InpNormTautData }; */
INP_ATOM_DATA *inp_cur_data = NULL;

COMPONENT_TREAT_INFO *cti = NULL;



    inp_cur_data = &(genctl->InpCurAtData[iINChI][i]);
    
    cti = &(genctl->cti[iINChI][i]);

    inp_norm_data[TAUT_NON] = &(genctl->InpNormAtData[iINChI][i]);
    inp_norm_data[TAUT_YES] = &(genctl->InpNormTautData[iINChI][i]);
    

    pINChI     = pINChI2[iINChI];
    pINChI_Aux = pINChI_Aux2[iINChI];
    
    
    InchiTimeGet( &ulTStart );
    bOrigCoord = !(ip->bINChIOutputOptions & (INCHI_OUT_NO_AUX_INFO | INCHI_OUT_SHORT_AUX_INFO));


    for ( k = 0; k < TAUT_NUM; k ++ ) 
    {
        cur_INChI[k] = pINChI[i][k];
        cur_INChI_Aux[k] = pINChI_Aux[i][k];
    }
    
    
    /*  allocate memory for non-tautimeric (k=0) and tautomeric (k=1) results */
    for ( k = 0; k < TAUT_NUM; k ++ ) 
    {
        int nAllocMode = (k==TAUT_YES? REQ_MODE_TAUT:0) |
                         (bTautFlagsDone & ( TG_FLAG_FOUND_ISOTOPIC_H_DONE |
                                             TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE ))?
                         (ip->nMode & REQ_MODE_ISO):0;

        if ( k==TAUT_NON && (ip->nMode & REQ_MODE_BASIC ) ||
             k==TAUT_YES && (ip->nMode & REQ_MODE_TAUT )     ) 
        {
            /*  alloc INChI and INChI_Aux */
            cur_INChI[k]     = Alloc_INChI( inp_cur_data->at, inp_cur_data->num_at, &inp_cur_data->num_bonds,
                                          &inp_cur_data->num_isotopic, nAllocMode );
            cur_INChI_Aux[k] = Alloc_INChI_Aux( inp_cur_data->num_at,
                                          inp_cur_data->num_isotopic, nAllocMode, bOrigCoord );
            if ( cur_INChI_Aux[k] ) 
                cur_INChI_Aux[k]->bIsIsotopic = inp_cur_data->num_isotopic;
            /*  alloc memory for the output structure: non-tautomeric and tautomeric (for displaying) */
            CreateInpAtomData( inp_norm_data[k], inp_cur_data->num_at, k );
        } 
        else 
            FreeInpAtomData( inp_norm_data[k] );
        
    }
    
    lElapsedTime = InchiTimeElapsed( &ulTStart );
    if ( ip->msec_MaxTime ) 
        ip->msec_LeftTime -= lElapsedTime;
    sd->ulStructTime += lElapsedTime;


    /******************************************************
     *
     *  Get one component canonical numberings, etc.
     *
     ******************************************************/

    /*
     * Create_INChI() return value:
     * num_at <= 0: error code
     * num_at >  0: number of atoms (excluding terminal hydrogen atoms)
     * inp_norm_data[0] => non-tautomeric, inp_norm_data[1] => tautomeric
     */
    InchiTimeGet( &ulTStart );
    if ( ip->msec_MaxTime ) 
    {
        ulTEnd = ulTStart;
        pulTEnd = &ulTEnd;
        if ( ip->msec_LeftTime > 0 ) 
            InchiTimeAddMsec( pulTEnd, ip->msec_LeftTime );
    }


    cti->nUserMode = ip->nMode;

    /* vABParityUnknown holds actual value of an internal constant signifying       */
    /* unknown parity: either the same as for undefined parity (default==standard)  */
    /*  or a specific one (non-std; requested by SLUUD switch).                     */
    cti->vABParityUnknown = AB_PARITY_UNDF;
    if ( 0 != ( ip->nMode & REQ_MODE_DIFF_UU_STEREO) ) 
    {
        /* Make labels for unknown and undefined stereo different */
        cti->vABParityUnknown = AB_PARITY_UNKN;
    }

    num_at = Normalization_step( cur_INChI, cur_INChI_Aux, inp_cur_data->at,
                          inp_norm_data, inp_cur_data->num_at,
                          &bTautFlags, &bTautFlagsDone, cti);

    SetConnectedComponentNumber( inp_cur_data->at, inp_cur_data->num_at, i+1 ); /*  normalization alters structure component number */

    for ( k = 0; k < TAUT_NUM; k ++ ) 
    {
        if ( cur_INChI_Aux[k] && cur_INChI_Aux[k]->nNumberOfAtoms > 0 ) 
        {
            pncFlags->bNormalizationFlags[iINChI][k] |= cur_INChI_Aux[k]->bNormalizationFlags;
            pncFlags->bTautFlags[iINChI][k]          |= cur_INChI_Aux[k]->bTautFlags;
            pncFlags->bTautFlagsDone[iINChI][k]      |= cur_INChI_Aux[k]->bTautFlagsDone;
            pncFlags->nCanonFlags[iINChI][k]         |= cur_INChI_Aux[k]->nCanonFlags;
        }
    }

    /*  Detect errors */
    if ( num_at < 0 ) 
        sd->nErrorCode = num_at;
    else if ( num_at == 0 ) 
        sd->nErrorCode = -1;
    else if ( cur_INChI[TAUT_NON] && cur_INChI[TAUT_NON]->nErrorCode ) 
        /*  non-tautomeric error */
        sd->nErrorCode = cur_INChI[TAUT_NON]->nErrorCode;
    else if ( cur_INChI[TAUT_YES] && cur_INChI[TAUT_YES]->nErrorCode ) 
        /*  tautomeric error */
        sd->nErrorCode = cur_INChI[TAUT_YES]->nErrorCode;
    
    /*  detect and store stereo warnings */
    if ( !sd->nErrorCode ) 
        GetProcessingWarnings(cur_INChI, inp_norm_data, sd);
    

    lElapsedTime = InchiTimeElapsed( &ulTStart );
    if ( ip->msec_MaxTime ) 
        ip->msec_LeftTime -= lElapsedTime;
    
    sd->ulStructTime += lElapsedTime;
#ifndef INCHI_LIBRARY
    /*  Display the results */
    if ( ip->bDisplay )
        eat_keyboard_input();
#endif
    /*  a) No matter what happened save the allocated INChI pointers */
    /*  save the INChI of the current component */

    InchiTimeGet( &ulTStart );
    for ( k = 0; k < TAUT_NUM; k ++ ) 
    {
        pINChI[i][k]     = cur_INChI[k];
        pINChI_Aux[i][k] = cur_INChI_Aux[k];

        cur_INChI[k]     = NULL;
        cur_INChI_Aux[k] = NULL;
    }

    /*  b) Count one component structure and/or INChI results only if there was no error */
    /*     Set inp_norm_data[j]->num_removed_H = number of removed explicit H           */

    if ( !sd->nErrorCode ) 
    {

        /*  find where the current processed structure is located */
        int cur_is_in_non_taut = (pINChI[i][TAUT_NON] && pINChI[i][TAUT_NON]->nNumberOfAtoms>0);
        int cur_is_in_taut     = (pINChI[i][TAUT_YES] && pINChI[i][TAUT_YES]->nNumberOfAtoms>0);
        int cur_is_non_taut = cur_is_in_non_taut && 0 == pINChI[i][TAUT_NON]->lenTautomer ||
                              cur_is_in_taut     && 0 == pINChI[i][TAUT_YES]->lenTautomer;
        int cur_is_taut     = cur_is_in_taut     && 0 <  pINChI[i][TAUT_YES]->lenTautomer;

        if ( cur_is_non_taut + cur_is_taut ) 
        {
            /*  count tautomeric and non-tautomeric components of the structures */
            int j1 = cur_is_in_non_taut? TAUT_NON:TAUT_YES;
            int j2 = cur_is_in_taut?     TAUT_YES:TAUT_NON;
            int j;
            sd->num_non_taut[iINChI] += cur_is_non_taut;
            sd->num_taut[iINChI]     += cur_is_taut;
            for ( j = j1; j <= j2; j ++ ) 
            {
                int bIsotopic = (pINChI[i][j]->nNumberOfIsotopicAtoms ||
                                 pINChI[i][j]->nNumberOfIsotopicTGroups ||
                                 pINChI[i][j]->nPossibleLocationsOfIsotopicH && pINChI[i][j]->nPossibleLocationsOfIsotopicH[0]>1);
                if ( j == TAUT_YES ) {
                    bIsotopic |= (0 < pINChI_Aux[i][j]->nNumRemovedIsotopicH[0] + 
                                      pINChI_Aux[i][j]->nNumRemovedIsotopicH[1] +
                                      pINChI_Aux[i][j]->nNumRemovedIsotopicH[2]);
                }
                inp_norm_data[j]->bExists = 1; /*  j=0: non-taut exists, j=1: taut exists */
                inp_norm_data[j]->bHasIsotopicLayer = bIsotopic;
            }
        }
    }

    if ( sd->nErrorCode==CT_OUT_OF_RAM || sd->nErrorCode==CT_USER_QUIT_ERR ) 
        ret = _IS_FATAL;
    else if ( sd->nErrorCode ) 
        ret = _IS_ERROR;
    
    lElapsedTime = InchiTimeElapsed( &ulTStart );
    if ( ip->msec_MaxTime ) 
        ip->msec_LeftTime -= lElapsedTime;
    
    sd->ulStructTime += lElapsedTime;
    return ret;
}



/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
int CanonOneComponentINChI(INCHIGEN_CONTROL *genctl, int iINChI, int i)
{
STRUCT_DATA *sd = &(genctl->StructData);
INPUT_PARMS *ip = &(genctl->InpParms);
PINChI2 **pINChI2 = genctl->pINChI; 
PINChI_Aux2 **pINChI_Aux2 = genctl->pINChI_Aux;
NORM_CANON_FLAGS *pncFlags = &(genctl->ncFlags);

    inchiTime     ulTStart, ulTEnd, *pulTEnd = NULL;
    int           k, num_at, ret = 0;
    INChI       *cur_INChI[TAUT_NUM];
    INChI_Aux   *cur_INChI_Aux[TAUT_NUM];
    long          lElapsedTime;
    /*
    PINChI2     *pINChI     = pINChI2[iINChI];
    PINChI_Aux2 *pINChI_Aux = pINChI_Aux2[iINChI];
    */

PINChI2     *pINChI     = NULL;
PINChI_Aux2 *pINChI_Aux = NULL;
    
INP_ATOM_DATA *inp_norm_data[TAUT_NUM]; /*  = { &InpNormAtData, &InpNormTautData }; */
INP_ATOM_DATA *inp_cur_data = NULL;

COMPONENT_TREAT_INFO *cti = NULL;


    inp_cur_data = &(genctl->InpCurAtData[iINChI][i]);
    
    cti = &(genctl->cti[iINChI][i]);

    inp_norm_data[TAUT_NON] = &(genctl->InpNormAtData[iINChI][i]);
    inp_norm_data[TAUT_YES] = &(genctl->InpNormTautData[iINChI][i]);
    
    

    pINChI     = pINChI2[iINChI];
    pINChI_Aux = pINChI_Aux2[iINChI];

    InchiTimeGet( &ulTStart );

    for ( k = 0; k < TAUT_NUM; k ++ ) 
    {
        cur_INChI[k] = pINChI[i][k];
        cur_INChI_Aux[k] = pINChI_Aux[i][k];
    }

    
    
    lElapsedTime = InchiTimeElapsed( &ulTStart );
    if ( ip->msec_MaxTime ) {
        ip->msec_LeftTime -= lElapsedTime;
    }
    sd->ulStructTime += lElapsedTime;


    /******************************************************
     *
     *  Get one component canonical numberings, etc.
     *
     ******************************************************/

    /*
     * Create_INChI() return value:
     * num_at <= 0: error code
     * num_at >  0: number of atoms (excluding terminal hydrogen atoms)
     * inp_norm_data[0] => non-tautomeric, inp_norm_data[1] => tautomeric
     */
    InchiTimeGet( &ulTStart );
    if ( ip->msec_MaxTime ) {
        ulTEnd = ulTStart;
        pulTEnd = &ulTEnd;
        if ( ip->msec_LeftTime > 0 ) {
            InchiTimeAddMsec( pulTEnd, ip->msec_LeftTime );
        }
    }
    num_at = Canonicalization_step(cur_INChI, cur_INChI_Aux, inp_norm_data, 
                                   pulTEnd, NULL, sd->pStrErrStruct, cti);

    num_at = cti->num_atoms;

    SetConnectedComponentNumber( inp_cur_data->at, inp_cur_data->num_at, i+1 ); /*  normalization alters structure component number */
    for ( k = 0; k < TAUT_NUM; k ++ ) 
    {
        if ( cur_INChI_Aux[k] && cur_INChI_Aux[k]->nNumberOfAtoms > 0 ) 
        {
            pncFlags->bNormalizationFlags[iINChI][k] |= cur_INChI_Aux[k]->bNormalizationFlags;
            pncFlags->bTautFlags[iINChI][k]          |= cur_INChI_Aux[k]->bTautFlags;
            pncFlags->bTautFlagsDone[iINChI][k]      |= cur_INChI_Aux[k]->bTautFlagsDone;
            pncFlags->nCanonFlags[iINChI][k]         |= cur_INChI_Aux[k]->nCanonFlags;
        }
    }

    /*  Detect errors */
    if ( num_at < 0 ) 
        sd->nErrorCode = num_at;
    else if ( num_at == 0 ) 
        sd->nErrorCode = -1;
    else if ( cur_INChI[TAUT_NON] && cur_INChI[TAUT_NON]->nErrorCode ) 
        /*  non-tautomeric error */
        sd->nErrorCode = cur_INChI[TAUT_NON]->nErrorCode;
    else if ( cur_INChI[TAUT_YES] && cur_INChI[TAUT_YES]->nErrorCode ) 
        /*  tautomeric error */
        sd->nErrorCode = cur_INChI[TAUT_YES]->nErrorCode;
    
    /*  detect and store stereo warnings */
    if ( !sd->nErrorCode ) 
        GetProcessingWarnings(cur_INChI, inp_norm_data, sd);
    

    lElapsedTime = InchiTimeElapsed( &ulTStart );
    if ( ip->msec_MaxTime ) 
        ip->msec_LeftTime -= lElapsedTime;
    
    sd->ulStructTime += lElapsedTime;
#ifndef INCHI_LIBRARY
    /*  Display the results */
    if ( ip->bDisplay )
        eat_keyboard_input();
#endif
    /*  a) No matter what happened save the allocated INChI pointers */
    /*  save the INChI of the current component */

    InchiTimeGet( &ulTStart );
    for ( k = 0; k < TAUT_NUM; k ++ ) 
    {
        pINChI[i][k]     = cur_INChI[k];
        pINChI_Aux[i][k] = cur_INChI_Aux[k];

        cur_INChI[k]     = NULL;
        cur_INChI_Aux[k] = NULL;
    }

    /*  b) Count one component structure and/or INChI results only if there was no error */
    /*     Set inp_norm_data[j]->num_removed_H = number of removed explicit H           */

    if ( !sd->nErrorCode ) 
    {

        /*  find where the current processed structure is located */
        int cur_is_in_non_taut = (pINChI[i][TAUT_NON] && pINChI[i][TAUT_NON]->nNumberOfAtoms>0);
        int cur_is_in_taut     = (pINChI[i][TAUT_YES] && pINChI[i][TAUT_YES]->nNumberOfAtoms>0);
        int cur_is_non_taut = cur_is_in_non_taut && 0 == pINChI[i][TAUT_NON]->lenTautomer ||
                              cur_is_in_taut     && 0 == pINChI[i][TAUT_YES]->lenTautomer;
        int cur_is_taut     = cur_is_in_taut     && 0 <  pINChI[i][TAUT_YES]->lenTautomer;

        if ( cur_is_non_taut + cur_is_taut ) 
        {
            /*  count tautomeric and non-tautomeric components of the structures */
            int j1 = cur_is_in_non_taut? TAUT_NON:TAUT_YES;
            int j2 = cur_is_in_taut?     TAUT_YES:TAUT_NON;
            int j;
            sd->num_non_taut[iINChI] += cur_is_non_taut;
            sd->num_taut[iINChI]     += cur_is_taut;
            for ( j = j1; j <= j2; j ++ ) 
            {
                int bIsotopic = (pINChI[i][j]->nNumberOfIsotopicAtoms ||
                                 pINChI[i][j]->nNumberOfIsotopicTGroups ||
                                 pINChI[i][j]->nPossibleLocationsOfIsotopicH && pINChI[i][j]->nPossibleLocationsOfIsotopicH[0]>1);
                if ( j == TAUT_YES ) {
                    bIsotopic |= (0 < pINChI_Aux[i][j]->nNumRemovedIsotopicH[0] + 
                                      pINChI_Aux[i][j]->nNumRemovedIsotopicH[1] +
                                      pINChI_Aux[i][j]->nNumRemovedIsotopicH[2]);
                }
                inp_norm_data[j]->bExists = 1; /*  j=0: non-taut exists, j=1: taut exists */
                inp_norm_data[j]->bHasIsotopicLayer = bIsotopic;
            }
        }
    }

    if ( sd->nErrorCode==CT_OUT_OF_RAM || sd->nErrorCode==CT_USER_QUIT_ERR ) 
        ret = _IS_FATAL;
    else if ( sd->nErrorCode ) 
        ret = _IS_ERROR;
    
    lElapsedTime = InchiTimeElapsed( &ulTStart );
    if ( ip->msec_MaxTime ) 
        ip->msec_LeftTime -= lElapsedTime;
    
    sd->ulStructTime += lElapsedTime;
    return ret;
}





/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
int  Normalization_step( INChI **ppINChI, 
                        INChI_Aux **ppINChI_Aux, 
                        inp_ATOM *inp_at, 
                        INP_ATOM_DATA *out_norm_data[2],
                        int num_inp_at, 
                        INCHI_MODE *pbTautFlags, 
                        INCHI_MODE *pbTautFlagsDone,
                        COMPONENT_TREAT_INFO *z)
{

int i,ret=0;



T_GROUP_INFO * /*const*/  t_group_info        = &(z->vt_group_info);                
T_GROUP_INFO * /*const*/  t_group_info_orig   = &(z->vt_group_info_orig);
    
    
    BCN *pBCN = &(z->Bcn);

    /*^^^ */
    z->fix_isofixedh = 0;
    z->fix_termhchrg = 0;
    /*^^^ */
    #if( FIX_ISO_FIXEDH_BUG == 1 )
        if (TG_FLAG_FIX_ISO_FIXEDH_BUG & *pbTautFlags)
            z->fix_isofixedh = 1;
    #endif
#if( FIX_TERM_H_CHRG_BUG == 1 )
    if (TG_FLAG_FIX_TERM_H_CHRG_BUG & *pbTautFlags)
            z->fix_termhchrg = 1;
    #endif


    z->bPointedEdgeStereo = ((TG_FLAG_POINTED_EDGE_STEREO & *pbTautFlags)? PES_BIT_POINT_EDGE_STEREO:0)
                                  | ((TG_FLAG_PHOSPHINE_STEREO    & *pbTautFlags)? PES_BIT_PHOSPHINE_STEREO :0)
                                  | ((TG_FLAG_ARSINE_STEREO       & *pbTautFlags)? PES_BIT_ARSINE_STEREO    :0)
                                  | ((TG_FLAG_FIX_SP3_BUG         & *pbTautFlags)? PES_BIT_FIX_SP3_BUG      :0);
    z->bTautFlags         = (*pbTautFlags     & (~(INCHI_MODE)TG_FLAG_ALL_TAUTOMERIC) );
    z->bTautFlagsDone     = (*pbTautFlagsDone /*& (~(INCHI_MODE)TG_FLAG_ALL_TAUTOMERIC) */);

    z->out_at = NULL;       /*, *norm_at_fixed_bonds[TAUT_NUM]; */ /*  = {out_norm_nontaut_at, out_norm_taut_at} ; */

    /* Init: internal structs */

    memset( z->s, 0, sizeof(z->s) );    

    if ( pBCN ) memset( pBCN, 0, sizeof( pBCN[0] ) );    

    memset( t_group_info, 0, sizeof(*t_group_info) );
    memset( t_group_info_orig, 0, sizeof(*t_group_info_orig) );


    /* Allocate: at[] */
    
    for ( i = 0; i < TAUT_NUM; i ++ ) 
    {
        if ( out_norm_data[i]->at ) 
        {
            z->at[i] = 
                    (sp_ATOM  *) inchi_malloc( num_inp_at * sizeof(*(z->at[0])) );
            
            if ( !z->at[i] ) ret = -1;
        } 
        else 
            z->at[i] = NULL;

    }

    if ( !out_norm_data[TAUT_NON]->at && !out_norm_data[TAUT_YES]->at || !inp_at || ret ) 
    {
        ret = -1;
        goto exit_function;
    }

    /* the first struct to process: tautomeric if exists else non-tautomeric */
    
    z->out_at = out_norm_data[TAUT_YES]->at? out_norm_data[TAUT_YES]->at : out_norm_data[TAUT_NON]->at;
    
    /* copy the input structure to be normalized to the buffer for the normalization data */
    
    memcpy( z->out_at, inp_at, num_inp_at*sizeof(z->out_at[0]) );

    
    /*  tautomeric groups setting */
    
    t_group_info->bIgnoreIsotopic = 0;   /*  include tautomeric group isotopic info in MarkTautomerGroups() */
    t_group_info->bTautFlags      = *pbTautFlags;
    t_group_info->bTautFlagsDone  = *pbTautFlagsDone;


    /*  Preprocess the structure; here THE NUMBER OF ATOMS MAY BE REDUCED */
    /*  ??? Ambiguity: H-D may become HD or DH (that is, H+implicit D or D+implicit H) */
    
    if ( TG_FLAG_H_ALREADY_REMOVED & z->bTautFlags ) 
    {
        INP_ATOM_DATA *out_norm_data1 = out_norm_data[TAUT_YES]->at? out_norm_data[TAUT_YES] :
                                        out_norm_data[TAUT_NON]->at? out_norm_data[TAUT_NON] : NULL;
        if ( out_norm_data1 ) 
        {
            z->num_at_tg     =
            z->num_atoms     = out_norm_data1->num_at - out_norm_data1->num_removed_H;
            z->num_deleted_H = out_norm_data1->num_removed_H;
            t_group_info->tni.nNumRemovedExplicitH = z->num_deleted_H;
        } 
        else 
        {
            ret = -1;
            goto exit_function;
        }
    } 
    else 
    {
        z->num_at_tg = z->num_atoms = remove_terminal_HDT( num_inp_at, z->out_at, z->fix_termhchrg);
        
        z->num_deleted_H = num_inp_at - z->num_atoms;
        t_group_info->tni.nNumRemovedExplicitH = z->num_deleted_H;

        add_DT_to_num_H( z->num_atoms, z->out_at );
    }
    
    
    
    /* fix_odd_things( z->num_atoms, z->out_at );*/


#if( FIND_RING_SYSTEMS == 1 )
    MarkRingSystemsInp( z->out_at, z->num_atoms, 0 );
#endif

    /*  duplicate the preprocessed structure so that all supplied out_norm_data[]->at buffers are filled */
    if ( z->out_at != out_norm_data[TAUT_YES]->at && out_norm_data[TAUT_YES]->at ) 
        memcpy( out_norm_data[TAUT_YES]->at, z->out_at, num_inp_at*sizeof(z->out_at[0]) );
    
    if ( out_norm_data[TAUT_YES]->at_fixed_bonds && out_norm_data[TAUT_YES]->at ) 
        memcpy( out_norm_data[TAUT_YES]->at_fixed_bonds, z->out_at, num_inp_at*sizeof(z->out_at[0]) );
    
    if ( z->out_at != out_norm_data[TAUT_NON]->at && out_norm_data[TAUT_NON]->at ) 
        memcpy( out_norm_data[TAUT_NON]->at, z->out_at, num_inp_at*sizeof(z->out_at[0]) );
    

    /*******************************************************************************
     * ??? not true ??? duplicate inp_at and keep inp_at[] unchanged after terminal hydrogens removal
     * set stereo parities in taut_at[], non_taut_at[]
     * obtain max. lenghts of the name stereo parts
     * Ignore absence/presence of isotopic stereo for now
     * mark isotopic atoms
     *******************************************************************************/
    
    if ( out_norm_data[TAUT_YES]->at && z->at[TAUT_YES] ) 
    {
    
        /* final normalization of possibly tautomeric structure */

        ret = 
                mark_alt_bonds_and_taut_groups ( out_norm_data[TAUT_YES]->at, 
                                                 out_norm_data[TAUT_YES]->at_fixed_bonds, 
                                                 z->num_atoms,
                                                 t_group_info, 
                                                 NULL, NULL );
        
        if ( ret < 0 ) 
            goto exit_function;/*  out of RAM or other normalization problem */
        
        z->num_taut_at = ret; /* number of atoms without removed H? */
        z->num_deleted_H_taut = t_group_info->tni.nNumRemovedExplicitH;
        out_norm_data[TAUT_YES]->num_at              = z->num_atoms + z->num_deleted_H_taut; /* protons might have been removed */
        out_norm_data[TAUT_YES]->num_removed_H       = z->num_deleted_H_taut;
        out_norm_data[TAUT_YES]->nNumRemovedProtons += t_group_info->tni.nNumRemovedProtons;
        
        for ( i = 0; i < NUM_H_ISOTOPES; i ++ ) 
        {
            out_norm_data[TAUT_YES]->nNumRemovedProtonsIsotopic[i] += t_group_info->tni.nNumRemovedProtonsIsotopic[i] /*+ t_group_info->num_iso_H[i]*/;
            out_norm_data[TAUT_YES]->num_iso_H[i]                  += t_group_info->num_iso_H[i];
        }

        /* mark deleted isolated tautomeric H(+) */
        
        if ( z->num_taut_at == 1 && out_norm_data[TAUT_YES]->at[0].at_type == ATT_PROTON &&
             t_group_info && t_group_info->tni.nNumRemovedProtons == 1 ) 
        {
            out_norm_data[TAUT_YES]->bDeleted = 1;
            
            FreeInpAtom( &out_norm_data[TAUT_YES]->at_fixed_bonds );
        } 
        else if ( (t_group_info->tni.bNormalizationFlags & FLAG_NORM_CONSIDER_TAUT) &&
                   out_norm_data[TAUT_YES]->at_fixed_bonds) 
        {
             out_norm_data[TAUT_YES]->bTautPreprocessed = 1;
        }

        out_norm_data[TAUT_YES]->bTautFlags     = *pbTautFlags     = t_group_info->bTautFlags;
        out_norm_data[TAUT_YES]->bTautFlagsDone = *pbTautFlagsDone = t_group_info->bTautFlagsDone;
        out_norm_data[TAUT_YES]->bNormalizationFlags = t_group_info->tni.bNormalizationFlags;
        
        /* create internal sp_ATOM at[] out of out_norm_data[]->at */
        
        inp2spATOM( out_norm_data[TAUT_YES]->at, num_inp_at, z->at[TAUT_YES] );

        
        
        /* set stereo parities to at[]; nUserMode: accept alt. stereo bonds, min ring size */
        
        ret = 
                set_stereo_parity( out_norm_data[TAUT_YES]->at, z->at[TAUT_YES], z->num_taut_at, z->num_deleted_H_taut,
                                 &(z->s[TAUT_YES].nMaxNumStereoAtoms), 
                                 &(z->s[TAUT_YES].nMaxNumStereoBonds), z->nUserMode,
                                 z->bPointedEdgeStereo, z->vABParityUnknown );
        
        if ( RETURNED_ERROR(ret) ) goto exit_function; /*  stereo bond error */
        

        z->s[TAUT_YES].bMayHaveStereo    = (z->s[TAUT_YES].nMaxNumStereoAtoms || z->s[TAUT_YES].nMaxNumStereoBonds);

        /* 
         * mark isotopic atoms and atoms that have non-tautomeric
         * isotopic terminal hydrogen atoms 1H, 2H(D), 3H(T)
         */

        z->s[TAUT_YES].num_isotopic_atoms = 
            
                    set_atom_iso_sort_keys( z->num_taut_at, z->at[TAUT_YES], t_group_info,
                                            &(z->s[TAUT_YES].bHasIsotopicTautGroups) );
        

        /**************************************************************************
         *  prepare tautomeric (if no tautomerism found then prepare non-tautomeric)
         *  structure for canonicalizaton:
         **************************************************************************
         *   remove t-groups that have no H,
         *   remove charges from t-groups if requested
         *   renumber t-groups and find final t_group_info->num_t_groups
         *   add to t-groups lists of endpoints tgroup->nEndpointAtomNumber[]
         *   calculate length of the t-group part of the connection table
         **************************************************************************/
        
        z->s[TAUT_YES].nLenLinearCTTautomer = 
                
                    CountTautomerGroups( z->at[TAUT_YES], z->num_taut_at, t_group_info );
        

        if ( RETURNED_ERROR(z->s[TAUT_YES].nLenLinearCTTautomer) ) 
        { 
            /* added error treatment 9-11-2003 */
            ret = z->s[TAUT_YES].nLenLinearCTTautomer;
            goto exit_function;
            /*  error has happened; no breakpoint here
            z->s[TAUT_YES].nLenLinearCTTautomer = 0;
            */
        } 
        else if ( z->s[TAUT_YES].nLenLinearCTTautomer > 0 ) 
        {
            z->num_at_tg = z->num_taut_at+t_group_info->num_t_groups;

            /*  ??? -not true- create t_group_info_orig for multiple calls with atom renumbering */

            make_a_copy_of_t_group_info( t_group_info_orig /* dest*/, t_group_info /* source*/ );

            /*  mark isotopic tautomer groups: calculate t_group->iWeight */
            z->s[TAUT_YES].nLenLinearCTIsotopicTautomer=set_tautomer_iso_sort_keys( t_group_info );
            if ( z->s[TAUT_YES].nLenLinearCTIsotopicTautomer < 0 ) 
            {
                /* ??? -error cannot happen- error has happened; no breakpoint here */
                z->s[TAUT_YES].nLenLinearCTIsotopicTautomer = 0;
            }
            out_norm_data[TAUT_YES]->bTautomeric = z->s[TAUT_YES].nLenLinearCTTautomer;
        }
        
        /*  new variable: z->s[TAUT_YES].nLenCT introduced 7-22-2002 */

        GetCanonLengths( z->num_taut_at, z->at[TAUT_YES], &(z->s[TAUT_YES]), t_group_info );


    } /* end of: final normalization of possibly tautomeric structure */
    


    if ( out_norm_data[TAUT_NON]->at && out_norm_data[TAUT_YES]->at && z->at[TAUT_NON] && !z->s[TAUT_YES].nLenLinearCTTautomer ) 
    {
        /* the structure is non-tautomeric: use tautomeric treatment results only for it */

        inchi_free( z->at[TAUT_NON] );
    
        z->at[TAUT_NON] = NULL;
    } 

    else if ( !out_norm_data[TAUT_NON]->at && out_norm_data[TAUT_YES]->at &&
         !z->at[TAUT_NON] && z->at[TAUT_YES] && !z->s[TAUT_YES].nLenLinearCTTautomer ) 
    {
        /* requested tautomeric; found non-tautomeric; it is located in out_norm_data[TAUT_YES]->at */
    
        out_norm_data[TAUT_YES]->bTautomeric = 0;
    }   

    else if ( out_norm_data[TAUT_NON]->at && z->at[TAUT_NON] ) 
    {
        /* the structure needs non-tautomeric treatment: final normalization of non-tautomeric structure */

        ret = 
                    mark_alt_bonds_and_taut_groups (out_norm_data[TAUT_NON]->at, 
                                                    NULL, 
                                                    z->num_atoms, 
                                                    NULL, 
                                                    &(z->bTautFlags), &(z->bTautFlagsDone) );
        
        if ( ret < 0 ) goto exit_function;  /*  out of RAM or other normalization problem */
        
        out_norm_data[TAUT_NON]->num_at        = z->num_atoms + z->num_deleted_H;
        out_norm_data[TAUT_NON]->num_removed_H = z->num_deleted_H;
        out_norm_data[TAUT_NON]->bTautFlags     = *pbTautFlags;
        out_norm_data[TAUT_NON]->bTautFlagsDone = *pbTautFlagsDone;
        out_norm_data[TAUT_NON]->bNormalizationFlags = 0;
        
        /* create internal sp_ATOM at[] out of out_norm_data[]->at */
        
        inp2spATOM( out_norm_data[TAUT_NON]->at, num_inp_at, z->at[TAUT_NON] );
        
        /* set stereo parities to at[]; nUserMode: accept alt. stereo bonds, min ring size */
        
        ret = 
                set_stereo_parity( out_norm_data[TAUT_NON]->at, 
                                   z->at[TAUT_NON], 
                                   z->num_atoms, z->num_deleted_H,
                                   &(z->s[TAUT_NON].nMaxNumStereoAtoms), 
                                   &(z->s[TAUT_NON].nMaxNumStereoBonds), 
                                   z->nUserMode,
                                   z->bPointedEdgeStereo, z->vABParityUnknown );

        if ( RETURNED_ERROR( ret ) ) goto exit_function; /*  stereo bond error */
        

        z->s[TAUT_NON].bMayHaveStereo = (z->s[TAUT_NON].nMaxNumStereoAtoms || z->s[TAUT_NON].nMaxNumStereoBonds);
        
        /* 
         * mark isotopic atoms and atoms that have non-tautomeric
         * isotopic terminal hydrogen atoms 1H, 2H(D), 3H(T)
         */
        
        z->s[TAUT_NON].num_isotopic_atoms =         
                    
                set_atom_iso_sort_keys( z->num_atoms, z->at[TAUT_NON], NULL, NULL );


        GetCanonLengths( z->num_atoms, z->at[TAUT_NON], &(z->s[TAUT_NON]), NULL);
        
        
        out_norm_data[TAUT_NON]->bTautomeric = 0;

    } /* the structure needs non-tautomeric treatment: final normalization of non-tautomeric structure */ 




    /**********************************************************/

    


    /*  common */
    z->bMayHaveStereo        = z->s[TAUT_YES].bMayHaveStereo || z->s[TAUT_NON].bMayHaveStereo;
    z->bHasIsotopicAtoms     = z->s[TAUT_NON].num_isotopic_atoms > 0 || z->s[TAUT_NON].bHasIsotopicTautGroups > 0 ||
                            z->s[TAUT_YES].num_isotopic_atoms > 0 || z->s[TAUT_YES].bHasIsotopicTautGroups > 0 ;
/*^^^ */
    if (z->fix_isofixedh)  /* 2008-03-21 DT */
        z->bHasIsotopicAtoms     = z->bHasIsotopicAtoms ||
                                z->s[TAUT_YES].nLenLinearCTTautomer > 0 && t_group_info &&
                                (0 < NUM_H_ISOTOPES && t_group_info->tni.nNumRemovedProtonsIsotopic[0] ||
                                 1 < NUM_H_ISOTOPES && t_group_info->tni.nNumRemovedProtonsIsotopic[1] ||
                                 2 < NUM_H_ISOTOPES && t_group_info->tni.nNumRemovedProtonsIsotopic[2]) ;
/*^^^ */
        z->bHasIsotopicAtoms     = z->bHasIsotopicAtoms ||
                                z->s[TAUT_YES].nLenIsotopicEndpoints > 1 && t_group_info &&
                                (t_group_info->bTautFlagsDone & (TG_FLAG_FOUND_ISOTOPIC_H_DONE|TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE));






    /* Set mode */

    /*  default mode */
    
    if ( !(z->nUserMode & REQ_MODE_DEFAULT) ) 
    {
        z->nUserMode |= REQ_MODE_DEFAULT;
    }
    
    
    /*  adjust the mode to the reality */
    
    if ( ( z->nUserMode & REQ_MODE_ISO ) && !z->bHasIsotopicAtoms ) 
    {
        z->nUserMode ^= REQ_MODE_ISO;
        z->nUserMode |= REQ_MODE_NON_ISO;  /*  at least one is needed */
    }
    
    if ( (z->nUserMode & REQ_MODE_STEREO) && ( z->nUserMode & REQ_MODE_ISO ) ) 
    {
        z->nUserMode |= REQ_MODE_ISO_STEREO;
    }
    
    if ( (z->nUserMode & REQ_MODE_STEREO) && !( z->nUserMode & REQ_MODE_NON_ISO ) ) 
    {
        z->nUserMode ^= REQ_MODE_STEREO;
    }
    
    if ( !z->bMayHaveStereo ) 
    {
        if ( z->nUserMode & REQ_MODE_STEREO )
            z->nUserMode ^= REQ_MODE_STEREO;
        if ( z->nUserMode & REQ_MODE_ISO_STEREO )
            z->nUserMode ^= REQ_MODE_ISO_STEREO;
    }

    if ( (z->nUserMode & REQ_MODE_BASIC) && (!out_norm_data[TAUT_NON]->at || !ppINChI[TAUT_NON] || !ppINChI_Aux[TAUT_NON] || !z->at[TAUT_NON]) ) 
    {
        z->nUserMode ^= REQ_MODE_BASIC;
    }
    if ( (z->nUserMode & REQ_MODE_TAUT) && (!out_norm_data[TAUT_YES]->at || !ppINChI[TAUT_YES] || !ppINChI_Aux[TAUT_YES] || !z->at[TAUT_YES]) ) 
    {
        z->nUserMode ^= REQ_MODE_TAUT;
    }
        
    
    /* Set n1, n2 according to the mode */

    switch ((int)z->nUserMode & (REQ_MODE_BASIC | REQ_MODE_TAUT)) 
    {
    case REQ_MODE_BASIC:
        z->n1 = TAUT_NON;
        z->n2 = TAUT_NON;
        break;
    case REQ_MODE_TAUT:
        z->n1 = TAUT_YES;
        z->n2 = TAUT_YES;
        break;
    case (REQ_MODE_BASIC | REQ_MODE_TAUT):
        z->n1 = TAUT_NON;
        z->n2 = TAUT_YES;
        break;
    default:
        /*  program error: inconsistent nUserMode or missing taut/non-taut allocation */ /*   <BRKPT> */
        ret = -3;
        goto exit_function; 
    }


    if ( ret == 0 ) 
        ret = z->num_atoms;
    
    /*  treat the results later */

exit_function:

    return ret;

}






/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
int  Canonicalization_step( INChI **ppINChI, 
                           INChI_Aux **ppINChI_Aux, 
                           INP_ATOM_DATA *out_norm_data[2],
                           struct tagInchiTime *ulMaxTime, 
                           T_GROUP_INFO *ti_out, 
                           char *pStrErrStruct,
                           COMPONENT_TREAT_INFO *z)
{
int i,ret=0, ret2=0;



T_GROUP_INFO * /*const*/  t_group_info        = &(z->vt_group_info);                
T_GROUP_INFO * /*const*/  t_group_info_orig   = &(z->vt_group_info_orig);
    
CANON_STAT  CS, CS2;
CANON_STAT *pCS  = &CS;
CANON_STAT *pCS2 = &CS2;  /*  save all allocations to avoid memory leaks in case Canon_INChI() removes the pointer */
    
    BCN *pBCN = &(z->Bcn);


    INChI     *pINChI=NULL;         /* added initialization 2006-03 */
    INChI_Aux *pINChI_Aux=NULL;     /* added initialization 2006-03 */



    /************************************************************
     *                                                          *
     *       Obtain all non-stereo canonical numberings         *
     *                                                          *
     ************************************************************/

    if ( (z->nUserMode & REQ_MODE_NON_ISO) && !(z->nUserMode & REQ_MODE_ISO) ) 
    {

        /* added for special non-isotopic test mode 2004-10-04 */
        if ( t_group_info ) 
        {
            t_group_info->bIgnoreIsotopic = 1;
            if ( t_group_info->nIsotopicEndpointAtomNumber ) 
            {
                t_group_info->nIsotopicEndpointAtomNumber[0] = inchi_min(1, t_group_info->nIsotopicEndpointAtomNumber[0]);
            }
            memset( t_group_info->num_iso_H, 0, sizeof(t_group_info->num_iso_H) );
            memset ( t_group_info->tni.nNumRemovedProtonsIsotopic, 0, sizeof(t_group_info->tni.nNumRemovedProtonsIsotopic));
            t_group_info->bTautFlagsDone &= ~(TG_FLAG_FOUND_ISOTOPIC_H_DONE|TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE);
        }

        for ( i = 0; i < TAUT_NUM; i ++ ) 
        {
            z->s[i].bHasIsotopicTautGroups = 0;
            z->s[i].bIgnoreIsotopic = 1;
            z->s[i].nLenIsotopic = 0;
            z->s[i].nLenIsotopicEndpoints = 0;
            z->s[i].nLenLinearCTIsotopicTautomer = 0;
            z->s[i].num_isotopic_atoms = 0;
        }
        z->bHasIsotopicAtoms = 0;
    }

    ret = GetBaseCanonRanking( z->num_atoms, z->num_at_tg, z->at, t_group_info, z->s, pBCN, ulMaxTime, z->fix_isofixedh );

    if ( ret < 0 ) 
        goto exit_function; /*  program error */
    

    /* added for special non-isotopic test mode 2004-10-04 */
    if ( !pBCN->ftcn[z->n1].PartitionCt.Rank ) 
        z->n1 = ALT_TAUT(z->n1);
    
    if ( !pBCN->ftcn[z->n2].PartitionCt.Rank ) 
        z->n2 = ALT_TAUT(z->n2);
    
    if ( z->n1 > z->n2 ) 
    {
        ret = CT_TAUCOUNT_ERR;
        goto exit_function; /*  program error */
    }


    /************************************************************
     *                                                          *
     *       Obtain stereo canonical numberings                 *
     *                                                          *
     ************************************************************/

    for ( i = z->n2; i >= z->n1 && !RETURNED_ERROR( ret ); i -- ) 
    {


        memset( pCS, 0, sizeof(*pCS) );

        switch (i) 
        {
        
        case TAUT_NON: /*  non-tautomeric */
        
            z->nMode  = 0;            
            z->nMode  = (z->s[i].nLenLinearCTTautomer == 0)? CANON_MODE_CT:CANON_MODE_TAUT;
            z->nMode |= (z->bHasIsotopicAtoms && (z->nUserMode & REQ_MODE_ISO))? CANON_MODE_ISO:0;
            z->nMode |= (z->s[TAUT_NON].bMayHaveStereo && (z->nUserMode & REQ_MODE_STEREO) )? CANON_MODE_STEREO:0;
            z->nMode |= (z->bHasIsotopicAtoms && z->s[TAUT_NON].bMayHaveStereo && (z->nUserMode & REQ_MODE_ISO_STEREO))? CANON_MODE_ISO_STEREO:0;
            z->nMode |= (z->nUserMode & REQ_MODE_NOEQ_STEREO   )? CMODE_NOEQ_STEREO    : 0;
            z->nMode |= (z->nUserMode & REQ_MODE_REDNDNT_STEREO)? CMODE_REDNDNT_STEREO : 0;
            z->nMode |= (z->nUserMode & REQ_MODE_NO_ALT_SBONDS )? CMODE_NO_ALT_SBONDS  : 0;

            /* 2010-01-12 */
            z->nMode |= (z->vABParityUnknown==AB_PARITY_UNDF)? 0 : REQ_MODE_DIFF_UU_STEREO;

            if ( (z->nMode & CANON_MODE_STEREO)     == CANON_MODE_STEREO ||
                 (z->nMode & CANON_MODE_ISO_STEREO) == CANON_MODE_ISO_STEREO ) 
            {
                z->nMode |= (z->nUserMode & REQ_MODE_RELATIVE_STEREO)? CMODE_RELATIVE_STEREO: 0;
                z->nMode |= (z->nUserMode & REQ_MODE_RACEMIC_STEREO )? CMODE_RACEMIC_STEREO : 0;
                z->nMode |= (z->nUserMode & REQ_MODE_SC_IGN_ALL_UU  )? CMODE_SC_IGN_ALL_UU  : 0;
                z->nMode |= (z->nUserMode & REQ_MODE_SB_IGN_ALL_UU  )? CMODE_SB_IGN_ALL_UU  : 0;
            }
            
            if ( ret= AllocateCS( pCS, z->num_atoms, z->num_atoms, z->s[TAUT_NON].nLenCT, z->s[TAUT_NON].nLenCTAtOnly,
                             z->s[TAUT_NON].nLenLinearCTStereoDble, z->s[TAUT_NON].nMaxNumStereoBonds,
                             z->s[TAUT_NON].nLenLinearCTStereoCarb, z->s[TAUT_NON].nMaxNumStereoAtoms,
                             0, 0, z->s[TAUT_NON].nLenIsotopic, z->nMode, pBCN ) ) 
            {
                goto exit_function;
            }
            
            *pCS2 = *pCS;
            break;


        case TAUT_YES: /*  tautomeric */

            z->nMode  = 0;
            z->nMode  = (z->s[i].nLenLinearCTTautomer == 0)? CANON_MODE_CT:CANON_MODE_TAUT;
            z->nMode |= (z->bHasIsotopicAtoms && (z->nUserMode & REQ_MODE_ISO) )? CANON_MODE_ISO:0;
            z->nMode |= (z->s[TAUT_YES].bMayHaveStereo && (z->nUserMode & REQ_MODE_STEREO) )? CANON_MODE_STEREO:0;
            z->nMode |= (z->bHasIsotopicAtoms && z->s[TAUT_YES].bMayHaveStereo && (z->nUserMode & REQ_MODE_ISO_STEREO))? CANON_MODE_ISO_STEREO:0;
            z->nMode |= (z->nUserMode & REQ_MODE_NOEQ_STEREO   )? CMODE_NOEQ_STEREO    : 0;
            z->nMode |= (z->nUserMode & REQ_MODE_REDNDNT_STEREO)? CMODE_REDNDNT_STEREO : 0;
            z->nMode |= (z->nUserMode & REQ_MODE_NO_ALT_SBONDS )? CMODE_NO_ALT_SBONDS  : 0;

            /* 2010-01-12 */
            z->nMode |= (z->vABParityUnknown==AB_PARITY_UNDF)? 0 : REQ_MODE_DIFF_UU_STEREO;

            if ( (z->nMode & CANON_MODE_STEREO)     == CANON_MODE_STEREO ||
                 (z->nMode & CANON_MODE_ISO_STEREO) == CANON_MODE_ISO_STEREO ) 
            {
                z->nMode |= (z->nUserMode & REQ_MODE_RELATIVE_STEREO)? CMODE_RELATIVE_STEREO: 0;
                z->nMode |= (z->nUserMode & REQ_MODE_RACEMIC_STEREO )? CMODE_RACEMIC_STEREO : 0;
                z->nMode |= (z->nUserMode & REQ_MODE_SC_IGN_ALL_UU  )? CMODE_SC_IGN_ALL_UU  : 0;
                z->nMode |= (z->nUserMode & REQ_MODE_SB_IGN_ALL_UU  )? CMODE_SB_IGN_ALL_UU  : 0;
            }
            
            if ( ret= AllocateCS( pCS, z->num_atoms, z->num_at_tg, z->s[TAUT_YES].nLenCT, z->s[TAUT_YES].nLenCTAtOnly,
                             z->s[TAUT_YES].nLenLinearCTStereoDble, z->s[TAUT_YES].nMaxNumStereoBonds,
                             z->s[TAUT_YES].nLenLinearCTStereoCarb, z->s[TAUT_YES].nMaxNumStereoAtoms,
                             z->s[TAUT_YES].nLenLinearCTTautomer, z->s[TAUT_YES].nLenLinearCTIsotopicTautomer,
                             z->s[TAUT_YES].nLenIsotopic, z->nMode, pBCN ) ) 
            {
                goto exit_function;
            }
            
            *pCS2 = *pCS;
            break;


        } /* switch () */


        /*  settings */
        pCS->lNumDecreasedCT               = -1;                             
        pCS->bDoubleBondSquare             =  DOUBLE_BOND_NEIGH_LIST? 2:0;  /*  2 => special mode */
        pCS->bIgnoreIsotopic               =  !((z->s[TAUT_NON].num_isotopic_atoms ||        
                                                z->s[TAUT_YES].num_isotopic_atoms ||    
                                                z->s[TAUT_YES].bHasIsotopicTautGroups) ||
                                               (z->nUserMode & REQ_MODE_NON_ISO) ||
                                               !(z->nUserMode & REQ_MODE_ISO));
        
        if ( (z->nUserMode & REQ_MODE_NON_ISO) && !(z->nUserMode & REQ_MODE_ISO) ) 
            pCS->bIgnoreIsotopic = 1; /* 10-04-2004 */
        
        
        if ( i == TAUT_YES ) 
        {  
            /* tautomeric */
            pCS->t_group_info                  = t_group_info; /*  ??? make a copy or reuse ???  */
            pCS->t_group_info->bIgnoreIsotopic = !(z->s[TAUT_YES].bHasIsotopicTautGroups ||
                                                   (z->nUserMode & REQ_MODE_NON_ISO) ||
                                                   !(z->nUserMode & REQ_MODE_ISO));
            if ( (z->nUserMode & REQ_MODE_NON_ISO) && !(z->nUserMode & REQ_MODE_ISO) ) 
                pCS->t_group_info->bIgnoreIsotopic = 1; /* 10-04-2004 */
            
        }
        pCS->ulTimeOutTime  = pBCN->ulTimeOutTime;
        /*=========== Obsolete Mode Bits (bit 0 is Least Significant Bit) ===========
         *
         *  Mode      Bits       Description                                
         *   '0' c    0          Only one connection table canonicalization 
         *   '1' C    1          Recalculate CT using fixed nSymmRank       
         *   '2' i    1|2        Isotopic canonicalization (internal)       
         *   '3' I    1|2|4      Isotopic canonicalization (output)
         *   '4' s    1|8        Stereo canonicalization                    
         *   '5' S    1|2|4|16   Stereo isotopic canonicalization           
         *   '6' A    1|2|4|8|16 Output All
         */

        /***************************************
           The last canonicalization step
         ***************************************/
        if ( pBCN ) 
        {
            /* USE_CANON2 == 1 */
            pCS->NeighList  = NULL;
            pCS->pBCN       = pBCN;
            ret = Canon_INChI( z->num_atoms, i?z->num_at_tg:z->num_atoms, 
                               z->at[i], pCS, z->nMode, i);
        } 
        else 
        {
            /* old way */
            pCS->NeighList  = CreateNeighList( z->num_atoms, i?z->num_at_tg:z->num_atoms, z->at[i], pCS->bDoubleBondSquare, pCS->t_group_info );
            pCS->pBCN       = NULL;
            ret = Canon_INChI( z->num_atoms, i?z->num_at_tg:z->num_atoms, 
                               z->at[i], pCS, z->nMode, i);
        }

        pINChI     = ppINChI[i];      /* pointers to already allocated still empty InChI */
        pINChI_Aux = ppINChI_Aux[i];
        
        if ( ret <= 0 ) 
        {
            /***************************************/
            /*  failure in Canon_INChI()            */
            /***************************************/
            pINChI->nErrorCode     = ret;
            pINChI_Aux->nErrorCode = ret;
        } 
        else 
        {
            /***************************************/
            /*  success Canon_INChI()               */
            /*  save canonicalization results in   */
            /*  pINChI and pINChI_Aux                */
            /***************************************/
            pINChI->nErrorCode     = 0;
            pINChI_Aux->nErrorCode = 0;
            pINChI->bDeleted = pINChI_Aux->bDeleted = out_norm_data[i]->bDeleted;
            pINChI_Aux->nCanonFlags         = pCS->nCanonFlags;
            pINChI_Aux->bTautFlags          = out_norm_data[i]->bTautFlags;
            pINChI_Aux->bTautFlagsDone      = out_norm_data[i]->bTautFlagsDone;
            pINChI_Aux->bNormalizationFlags = out_norm_data[i]->bNormalizationFlags;
            /*  may return an error or a warning */
            ret = FillOutINChIReducedWarn( pINChI, pINChI_Aux,
                               z->num_atoms, i?z->num_at_tg:z->num_atoms,
                               i?z->num_deleted_H_taut:z->num_deleted_H, z->at[i], 
                               out_norm_data[i]->at, pCS, i, z->nUserMode, pStrErrStruct );
            if ( RETURNED_ERROR( ret ) ) 
            {
                /* failure in FillOutINChI() */
                pINChI->nErrorCode      = ret;
                pINChI_Aux->nErrorCode  = ret;
            } 
            else 
            {
                /****************************/
                /* success in FillOutINChI() */
                /****************************/

                /* mark non-tautomeric representation as having another, tautomeric representation */
                if ( pINChI_Aux && z->s[TAUT_YES].nLenLinearCTTautomer ) 
                    pINChI_Aux->bIsTautomeric = z->s[TAUT_YES].nLenLinearCTTautomer;
                

                ret2 = CheckCanonNumberingCorrectness(z->num_atoms, 
                                                      i?z->num_at_tg:z->num_atoms,
                                                      z->at[i], pCS, i, pStrErrStruct );
                if (ret2) 
                {
                    pINChI->nErrorCode      = ret2;
                    pINChI_Aux->nErrorCode  = ret2;
                    ret = ret2;
                }
            } /* success in FillOutINChI */
        } /* success Canon_INChI */

        FreeNeighList( pCS->NeighList );
        DeAllocateCS( pCS2 );

        pINChI = NULL;      /* avoid dangling pointers */
        pINChI_Aux = NULL;  /* avoid dangling pointers */
    
    } /* for ( i = z->n2; i >= z->n1 && !RETURNED_ERROR( ret ); i -- )  */



    if ( ret == 0 ) 
        ret = z->num_atoms;


exit_function:

    DeAllocBCN( pBCN );
    
    if ( z->at[TAUT_YES] )
    {
        inchi_free( z->at[TAUT_YES] );
        z->at[TAUT_YES] = NULL;
    }
    
    if ( z->at[TAUT_NON] )
    {
        inchi_free( z->at[TAUT_NON] );
        z->at[TAUT_NON] = NULL;
    }

    if ( ti_out ) 
        *ti_out = *t_group_info;
    else 
    {
        free_t_group_info( t_group_info );
        t_group_info = NULL;
    }
    
    free_t_group_info( t_group_info_orig );

    return ret;
}




/****************************************************************************/
int CreateCompositeNormAtom(COMP_ATOM_DATA *composite_norm_data, INP_ATOM_DATA2 *all_inp_norm_data, int num_components)
{
    int i, j, jj, k, n, m, tot_num_at, tot_num_H, cur_num_at, cur_num_H, nNumRemovedProtons;
    int num_comp[TAUT_NUM+1], num_taut[TAUT_NUM+1], num_del[TAUT_NUM+1], num_at[TAUT_NUM+1], num_inp_at[TAUT_NUM+1];
    int ret = 0, indicator = 1;
    inp_ATOM *at, *at_from;
    memset( num_comp, 0, sizeof(num_comp) );
    memset( num_taut, 0, sizeof(num_taut) );
    memset( num_del, 0, sizeof(num_taut) );
    /* count taut and non-taut components */
    for ( j = 0; j < TAUT_NUM; j ++ ) {
        num_comp[j] = num_taut[j] = 0;
        for ( i = 0; i < num_components; i ++ ) {
            if ( all_inp_norm_data[i][j].bExists ) {
                num_del[j]  += (0 != all_inp_norm_data[i][j].bDeleted );
                num_comp[j] ++;
                num_taut[j] += (0 != all_inp_norm_data[i][j].bTautomeric);
            }
        }
    }
    /* count intermediate taut structure components */
    if ( num_comp[TAUT_YES] > num_del[TAUT_YES] && num_taut[TAUT_YES] ) {
        /*
        num_comp[TAUT_INI] = num_comp[TAUT_YES] - num_del[TAUT_YES];
        */

        for ( i = 0, j=TAUT_YES; i < num_components; i ++ ) {
            if ( all_inp_norm_data[i][j].bExists &&
                (all_inp_norm_data[i][j].bDeleted ||
                 all_inp_norm_data[i][j].bTautomeric &&
                 all_inp_norm_data[i][j].at_fixed_bonds &&
                 all_inp_norm_data[i][j].bTautPreprocessed) ) {
                num_comp[TAUT_INI] ++;
            }
        }

    }
    /* count atoms and allocate composite atom data */
    for ( jj = 0; jj <= TAUT_INI; jj ++ ) {
        num_at[jj] = num_inp_at[jj] = 0;
        j = inchi_min (jj, TAUT_YES);
        if ( num_comp[jj] ) {
            for ( i = 0; i < num_components; i ++ ) {
                if ( all_inp_norm_data[i][j].bDeleted )
                    continue;
                /* find k = the normaized structure index */
                if ( jj == TAUT_INI ) {
                    if ( all_inp_norm_data[i][j].bExists &&
                         all_inp_norm_data[i][j].at_fixed_bonds ) {
                        k = j;
                    } else
                    if ( all_inp_norm_data[i][ALT_TAUT(j)].bExists && !all_inp_norm_data[i][ALT_TAUT(j)].bDeleted &&
                         !all_inp_norm_data[i][j].bDeleted  ) {
                        k = ALT_TAUT(j);
                    } else
                    if ( all_inp_norm_data[i][j].bExists ) {
                        k = j;
                    } else {
                        continue;
                    }
                } else {
                    if ( all_inp_norm_data[i][j].bExists ) {
                        k = j;
                    } else
                    if ( all_inp_norm_data[i][ALT_TAUT(j)].bExists && !all_inp_norm_data[i][ALT_TAUT(j)].bDeleted) {
                        k = ALT_TAUT(j);
                    } else {
                        continue;
                    }
                }
                num_inp_at[jj] += all_inp_norm_data[i][k].num_at; /* all atoms including terminal H */
                num_at[jj]     += all_inp_norm_data[i][k].num_at - all_inp_norm_data[i][k].num_removed_H;
            }
            if ( num_inp_at[jj] ) {
                if ( !CreateCompAtomData( composite_norm_data+jj, num_inp_at[jj], num_components, jj == TAUT_INI ) )
                    goto exit_error;
                composite_norm_data[jj].num_removed_H = num_inp_at[jj] - num_at[jj];
            }
        }
    }
    /* fill out composite atom */
    for ( jj = 0; jj <= TAUT_INI; jj ++, indicator <<= 1 ) {
        j = inchi_min (jj, TAUT_YES);
        if ( num_comp[jj] ) {
            tot_num_at = 0;
            tot_num_H = 0;
            for ( i = 0; i < num_components; i ++ ) {
                if ( all_inp_norm_data[i][j].bDeleted ) {
                    composite_norm_data[jj].nNumRemovedProtons += all_inp_norm_data[i][j].nNumRemovedProtons;
                    for ( n = 0; n < NUM_H_ISOTOPES; n ++ ) {
                        composite_norm_data[jj].nNumRemovedProtonsIsotopic[n] += all_inp_norm_data[i][j].nNumRemovedProtonsIsotopic[n];
                    }
                    continue;
                }
                nNumRemovedProtons = 0;
                k = TAUT_NUM;
                /* find k = the normaized structure index */
                if ( jj == TAUT_INI ) {
                    if ( all_inp_norm_data[i][j].bExists && all_inp_norm_data[i][j].at_fixed_bonds ) {
                        k = j;
                    } else
                    if ( all_inp_norm_data[i][ALT_TAUT(j)].bExists ) {
                        k = ALT_TAUT(j);
                    } else
                    if ( all_inp_norm_data[i][j].bExists && !all_inp_norm_data[i][ALT_TAUT(j)].bDeleted ) {
                        k = j;
                    } else {
                        continue;
                    }
                } else {
                    if ( all_inp_norm_data[i][j].bExists ) {
                        k = j;
                    } else
                    if ( all_inp_norm_data[i][ALT_TAUT(j)].bExists && !all_inp_norm_data[i][ALT_TAUT(j)].bDeleted ) {
                        k = ALT_TAUT(j);
                    } else {
                        continue;
                    }
                }
                /* copy main atoms */
                cur_num_H  = all_inp_norm_data[i][k].num_removed_H;       /* number of terminal H atoms */
                cur_num_at = all_inp_norm_data[i][k].num_at - cur_num_H;  /* number of all but explicit terminal H atoms */

                if ( (tot_num_at + cur_num_at) > num_at[jj] ||
                     (num_at[jj] + tot_num_H + cur_num_H) > num_inp_at[jj] ) {
                    goto exit_error; /* miscount */
                }
                at      = composite_norm_data[jj].at+tot_num_at; /* points to the 1st destination atom */
                at_from = (jj == TAUT_INI && k == TAUT_YES && all_inp_norm_data[i][k].at_fixed_bonds)?
                          all_inp_norm_data[i][k].at_fixed_bonds : all_inp_norm_data[i][k].at;
                memcpy( at, at_from, sizeof(composite_norm_data[0].at[0]) * cur_num_at ); /* copy atoms except terminal H */
                /* shift neighbors of main atoms */
                for ( n = 0; n < cur_num_at; n ++, at ++ ) {
                    for ( m = 0; m < at->valence; m ++ ) {
                        at->neighbor[m] += tot_num_at;
                    }
                }
                /* copy explicit H */
                if ( cur_num_H ) {
                    at = composite_norm_data[jj].at+num_at[jj]+tot_num_H; /* points to the 1st destination atom */
                    memcpy( at, at_from+cur_num_at,
                            sizeof(composite_norm_data[0].at[0]) * cur_num_H );
                    /* shift neighbors of explicit H atoms */
                    for ( n = 0; n < cur_num_H; n ++, at ++ ) {
                        for ( m = 0; m < at->valence; m ++ ) {
                            at->neighbor[m] += tot_num_at;
                        }
                    }
                }
                /* composite counts */
                composite_norm_data[jj].bHasIsotopicLayer   |= all_inp_norm_data[i][k].bHasIsotopicLayer;
                composite_norm_data[jj].num_isotopic        += all_inp_norm_data[i][k].num_isotopic;
                composite_norm_data[jj].num_bonds           += all_inp_norm_data[i][k].num_bonds;
                composite_norm_data[jj].bTautomeric         += (j == jj) && all_inp_norm_data[i][k].bTautomeric;
                composite_norm_data[jj].nNumRemovedProtons  += all_inp_norm_data[i][k].nNumRemovedProtons;
                for ( n = 0; n < NUM_H_ISOTOPES; n ++ ) {
                    composite_norm_data[jj].nNumRemovedProtonsIsotopic[n] += all_inp_norm_data[i][k].nNumRemovedProtonsIsotopic[n];
                    composite_norm_data[jj].num_iso_H[n]                  += all_inp_norm_data[i][k].num_iso_H[n];
                }
                /*
                composite_norm_data[j].num_at            += cur_num_at + cur_num_H;
                composite_norm_data[j].num_removed_H     += cur_num_H;
                */
                /* total count */
                tot_num_at += cur_num_at;
                tot_num_H += cur_num_H;
                /* offset for the next component */
                if (  composite_norm_data[jj].nOffsetAtAndH ) {
                    composite_norm_data[jj].nOffsetAtAndH[2*i]   = tot_num_at;
                    composite_norm_data[jj].nOffsetAtAndH[2*i+1] = num_at[jj]+tot_num_H;
                }
            }
            if ( tot_num_at != num_at[jj] ||
                 num_at[jj] + tot_num_H  != num_inp_at[jj] ) {
                goto exit_error; /* miscount */
            }
            composite_norm_data[jj].bExists       = (tot_num_at>0);
            ret |= indicator;
        }
    }
    return ret;





exit_error:
    return ret;
}

int CreateCompAtomData( COMP_ATOM_DATA *inp_at_data, int num_atoms, int num_components, int bIntermediateTaut )
{
    FreeCompAtomData( inp_at_data );
    if ( (inp_at_data->at = CreateInpAtom( num_atoms )) &&
         (num_components <= 1 || bIntermediateTaut ||
            (inp_at_data->nOffsetAtAndH = (AT_NUMB*)inchi_calloc(sizeof(inp_at_data->nOffsetAtAndH[0]), 2*(num_components+1))))) {

        inp_at_data->num_at = num_atoms;
        inp_at_data->num_components = (num_components>1)? num_components : 0;
        return 1;
    }
    FreeCompAtomData( inp_at_data );
    return 0;
}


/**********************************************************************************************/
int FillOutINChIReducedWarn( INChI *pINChI, INChI_Aux *pINChI_Aux,
                 int num_atoms, int num_at_tg, int num_removed_H,
                 sp_ATOM *at, inp_ATOM *norm_at, CANON_STAT *pCS, int bTautomeric,
                 INCHI_MODE nUserMode, char *pStrErrStruct )
{
    int i, j, m, n, g, len, ii, ret=0;

    AT_NUMB   *pSymmRank, *pOrigNosInCanonOrd, *pConstitEquNumb, *pCanonOrd=NULL, *pCanonOrdInv=NULL, *pCanonOrdTaut;
    T_GROUP_INFO     *t_group_info = pCS->t_group_info;
    T_GROUP *t_group;
    int nErrorCode = 0;
    AT_NUMB *pCanonRank, *pCanonRankInv; /* canonical ranks of the atoms or tautomeric groups */
    AT_NUMB *pCanonRankAtoms=NULL, *pSortOrd = NULL;
    AT_RANK nMinOrd;
    INChI_Stereo *Stereo;
    int          bUseNumberingInv = 0, bUseIsotopicNumberingInv = 0;
    INCHI_MODE    nStereoUnmarkMode;

    /*AT_NUMB  *pCanonOrdNonIso = NULL, *pCanonOrdIso = NULL;*/
    /*AT_NUMB  *nOrigAtNosInCanonOrdNonIso = NULL, *nOrigAtNosInCanonOrdIso = NULL;*/

    /*  Check for warnings */
    if ( pCS->nLenLinearCTStereoCarb < 0 || pCS->nLenLinearCTStereoDble  < 0 ||
         pCS->nLenCanonOrdStereo    < 0 || pCS->nLenCanonOrdStereoTaut < 0) {
        nErrorCode     |= WARN_FAILED_STEREO;
    }
    if ( pCS->nLenLinearCTIsotopic < 0  || pCS->nLenLinearCTIsotopicTautomer < 0 ||
         pCS->nLenCanonOrdIsotopic < 0 || pCS->nLenCanonOrdIsotopicTaut    < 0  ) {
        nErrorCode     |= WARN_FAILED_ISOTOPIC;
    }
    if ( pCS->nLenLinearCTIsotopicStereoCarb < 0 || pCS->nLenLinearCTIsotopicStereoDble  < 0 ||
         pCS->nLenCanonOrdIsotopicStereo    < 0 || pCS->nLenCanonOrdIsotopicStereoTaut < 0) {
        nErrorCode     |= WARN_FAILED_ISOTOPIC_STEREO;
    }
    pCanonRankAtoms = (AT_NUMB *)inchi_calloc( num_at_tg+1, sizeof(pCanonRankAtoms[0]) );
    pSortOrd        = (AT_NUMB *)inchi_calloc( num_at_tg+1, sizeof(pSortOrd[0]) ); /*  must have more than num_atoms */

    if ( !pCanonRankAtoms || !pSortOrd ) {
        nErrorCode = 0;
        ret = CT_OUT_OF_RAM;  /*   <BRKPT> */
        pINChI->nErrorCode = pINChI_Aux->nErrorCode = CT_OUT_OF_RAM;
        goto exit_function;
    }

    /*  total charge */
    for ( i = 0, n = 0; i < num_atoms+num_removed_H; i ++ ) {
        n += at[i].charge;
    }
    pINChI->nTotalCharge = n;

    /*  number of atoms */
    pINChI->nNumberOfAtoms     = num_atoms;
    pINChI_Aux->nNumberOfAtoms = num_atoms;

    /* removed protons and detachable isotopic H */
    if ( bTautomeric && t_group_info ) {
        pINChI_Aux->nNumRemovedProtons = t_group_info->tni.nNumRemovedProtons;
        for ( i = 0; i < NUM_H_ISOTOPES; i ++ ) {
            pINChI_Aux->nNumRemovedIsotopicH[i] = t_group_info->num_iso_H[i] 
                                               + t_group_info->tni.nNumRemovedProtonsIsotopic[i];
        }
        if ( pINChI_Aux->bNormalizationFlags & FLAG_FORCE_SALT_TAUT ) {
            pINChI->nFlags |= INCHI_FLAG_HARD_ADD_REM_PROTON;
        }
/*^^^
        if ( pINChI_Aux->bNormalizationFlags & (FLAG_NORM_CONSIDER_TAUT &~FLAG_PROTON_CHARGE_CANCEL) ) {
            AddMOLfileError(pStrErrStruct, "Proton(s) added/removed");
        }

        if ( pINChI_Aux->bNormalizationFlags & FLAG_PROTON_CHARGE_CANCEL ) {
            AddMOLfileError(pStrErrStruct, "Charges neutralized");
        }
^^^*/
    }

    /* abs or rel stereo may establish one of two canonical numberings */
    if ( (pCS->nLenLinearCTStereoCarb > 0 || pCS->nLenLinearCTStereoDble > 0) &&
          pCS->nLenCanonOrdStereo > 0 &&
         (pCS->LinearCTStereoCarb && pCS->LinearCTStereoCarbInv ||
          pCS->LinearCTStereoDble && pCS->LinearCTStereoDbleInv) &&
          pCS->nCanonOrdStereo    && pCS->nCanonOrdStereoInv
       ) {

        pCanonRank    = pCanonRankAtoms;
        pCanonOrd     = pCS->nCanonOrdStereo;
        pCanonRankInv = pSortOrd;
        pCanonOrdInv  = pCS->nCanonOrdStereoInv;
        Stereo        = pINChI->Stereo;
        for ( i = 0; i < num_at_tg; i ++ ) {
            pCanonRankInv[pCanonOrdInv[i]] = 
            pCanonRank[pCanonOrd[i]]       = (AT_NUMB)(i+1);
        }
        /********************************************************************/
        /* copy stereo bonds and stereo centers; compare Inv and Abs stereo */
        /********************************************************************/
        nErrorCode = CopyLinearCTStereoToINChIStereo( Stereo,
                           pCS->LinearCTStereoCarb, pCS->nLenLinearCTStereoCarb,
                           pCS->LinearCTStereoDble, pCS->nLenLinearCTStereoDble
                           , pCanonOrd, pCanonRank, at, 0 /* non-isotopic */
                           , pCS->LinearCTStereoCarbInv
                           , pCS->LinearCTStereoDbleInv
                           , pCanonOrdInv, pCanonRankInv ); 

        if ( Stereo->t_parityInv && Stereo->nNumberInv ) {
            if ( nUserMode & REQ_MODE_RELATIVE_STEREO ) {
                pINChI->nFlags |= INCHI_FLAG_REL_STEREO;
            }
            if ( nUserMode & REQ_MODE_RACEMIC_STEREO ) {
                pINChI->nFlags |= INCHI_FLAG_RAC_STEREO;
            }
            if ( Stereo->nCompInv2Abs ) {
                if ( Stereo->nCompInv2Abs == -1 ) {
                    /* switch pointers in Stereo so that the stereo becomes the smallest (relative)  */
                    /* flag Stereo->nCompInv2Abs == -1 will keep track of this exchange */
                    AT_NUMB    *nNumberInv  = Stereo->nNumberInv;
                    S_CHAR     *t_parityInv = Stereo->t_parityInv;
                    Stereo->nNumberInv  = Stereo->nNumber;
                    Stereo->t_parityInv = Stereo->t_parity;
                    Stereo->nNumber     = nNumberInv;
                    Stereo->t_parity    = t_parityInv;
                    /* switch pointers to set rel. stereo to pINChI_Aux->nOrigAtNosInCanonOrd
                                       and inv. stereo to pINChI_Aux->nOrigAtNosInCanonOrdInv */
                    switch_ptrs( &pCanonRank, &pCanonRankInv );
                    switch_ptrs( &pCanonOrd,  &pCanonOrdInv  );
                    bUseNumberingInv    = 1; /* use inverted stereo numbering instead of normal */
                }
            }
        }

        for ( i = 0; i < num_atoms; i ++ ) {
            pINChI_Aux->nOrigAtNosInCanonOrdInv[i] = at[pCanonOrdInv[i]].orig_at_number;
            pINChI_Aux->nOrigAtNosInCanonOrd[i]    = at[pCanonOrd[i]].orig_at_number;
        }
        if ( bUseNumberingInv ) {
            /* switch ptrs back to avoid confusion */
            switch_ptrs( &pCanonRank, &pCanonRankInv );
            switch_ptrs( &pCanonOrd,  &pCanonOrdInv  );
            /* save inverted stereo ranks & order because it represents the smallest (relative) */
            memcpy( pCanonRank, pCanonRankInv, num_at_tg * sizeof(pCanonRank[0]) );
            /* change pCS->nCanonOrdStereo[] to inverted: */
            memcpy( pCanonOrd,  pCanonOrdInv, num_at_tg * sizeof(pCanonOrd[0]) );
        }
        pCanonRankInv = NULL;
        pCanonOrdInv  = NULL;
        pOrigNosInCanonOrd = NULL;

    } else { /*------------------------------ no stereo */

        pCanonOrd            = pCS->nLenCanonOrdStereo > 0? pCS->nCanonOrdStereo :
                               pCS->nLenCanonOrd       > 0? pCS->nCanonOrd : NULL;
        pCanonRank           = pCanonRankAtoms;
        pOrigNosInCanonOrd   = pINChI_Aux->nOrigAtNosInCanonOrd;
        if ( pCanonOrd && pCanonRank ) {
            for ( i = 0; i < num_atoms; i ++ ) {
                pCanonRank[pCanonOrd[i]]       = (AT_NUMB)(i+1);
                pOrigNosInCanonOrd[i]          = at[pCanonOrd[i]].orig_at_number;
            }
            for ( ; i < num_at_tg; i ++ ) {
                pCanonRank[pCanonOrd[i]]       = (AT_NUMB)(i+1);
            }
        }
    }
    /*pCanonOrdNonIso = pCanonOrd;*/  /* save for aux info */


    if ( pINChI_Aux->OrigInfo ) {
        /* charges, radicals, valences */
        for ( i = 0; i < num_atoms; i ++ ) {
            ii = pCanonOrd[i];
            if ( norm_at[ii].valence || norm_at[ii].num_H ) {
                pINChI_Aux->OrigInfo[i].cCharge  = norm_at[ii].charge;
                pINChI_Aux->OrigInfo[i].cRadical = (norm_at[ii].radical==RADICAL_SINGLET)? 0 :
                                                  (norm_at[ii].radical==RADICAL_DOUBLET)? 1 :
                                                  (norm_at[ii].radical==RADICAL_TRIPLET)? 2 :
                                                  norm_at[ii].radical? 3 : 0 ;
                pINChI_Aux->OrigInfo[i].cUnusualValence = 
                    get_unusual_el_valence( norm_at[ii].el_number, norm_at[ii].charge, norm_at[ii].radical,
                                            norm_at[ii].chem_bonds_valence, norm_at[ii].num_H, norm_at[ii].valence );
            } else {
                /* charge of a single atom component is in the INChI; valence = 0 is standard */
                pINChI_Aux->OrigInfo[i].cRadical = (norm_at[ii].radical==RADICAL_SINGLET)? 0 :
                                                  (norm_at[ii].radical==RADICAL_DOUBLET)? 1 :
                                                  (norm_at[ii].radical==RADICAL_TRIPLET)? 2 :
                                                  norm_at[ii].radical? 3 : 0 ;
            }

        }
    }

    /* non-isotopic canonical numbers and equivalence of atoms (Aux) */
    pConstitEquNumb      = pINChI_Aux->nConstitEquNumbers;  /*  contitutional equivalence */
    pSymmRank            = pCS->nSymmRank;
    if ( pCanonOrd && pCanonRank && pSymmRank && pConstitEquNumb ) {
        for ( i = 0; i < num_atoms; i ++ ) {
            pConstitEquNumb[i]       = pSymmRank[pCanonOrd[i]]; /*  constit. equ. ranks in order of canonical numbers */
            pSortOrd[i]              = i;
        }
        for ( ; i < num_at_tg; i ++ ) {
            pSortOrd[i]              = MAX_ATOMS; /* for debugging only */
        }
        pn_RankForSort  = pConstitEquNumb;
        qsort( pSortOrd, num_atoms, sizeof(pSortOrd[0]), CompRanksOrd );
        for ( i = 0, nMinOrd = pSortOrd[0], j = 1; j <= num_atoms; j ++ ) {
            if ( j == num_atoms || pConstitEquNumb[pSortOrd[i]] != pConstitEquNumb[pSortOrd[j]] ) {
                nMinOrd ++;
                if ( j - i > 1 ) {
                    /*  found a sequence of equivalent atoms: i..j-1 */
                    while ( i < j ) {
                        pConstitEquNumb[pSortOrd[i++]] = nMinOrd; /*  = min. canon. rank in the group of equ. atoms */
                    }
                    /*  at this point j == i */
                } else {
                    pConstitEquNumb[pSortOrd[i++]] = 0; /*  means the atom is not equivalent to any other */
                }
                nMinOrd = pSortOrd[j]; /*  at the end j = num_atoms */
            }
        }
    } else {
        nErrorCode  |= ERR_NO_CANON_RESULTS;
        ret = -1;  /*  program error; no breakpoint here */
        goto exit_function;
    }
    /*  atomic numbers from the Periodic Table */
    for ( i = 0; i < num_atoms; i ++ ) {
        pINChI->nAtom[i] = (int)at[pCanonOrd[i]].el_number;
    }
    /*  connection table: atoms only (before 7-29-2003 pCS->LinearCT2 contained non-isotopic CT) */
    if ( pCS->nLenLinearCTAtOnly <= 0 || !pCS->LinearCT || !pINChI->nConnTable ) {
        nErrorCode  |= ERR_NO_CANON_RESULTS;
        ret = -2;
        goto exit_function;
    }
    memcpy( pINChI->nConnTable, pCS->LinearCT, sizeof(pINChI->nConnTable[0])*pCS->nLenLinearCTAtOnly);
    pINChI->lenConnTable = pCS->nLenLinearCTAtOnly;
    
    /*  tautomeric group(s) canonical representation */
    len = 0;
    if ( bTautomeric && 0 < (n = SortTautomerGroupsAndEndpoints( t_group_info, num_atoms, num_at_tg, pCanonRank )) ) {
        /* SortTautomerGroupsAndEndpoints() produces canonically ordered t-groups */
        pINChI->nFlags |= (t_group_info->bTautFlagsDone & TG_FLAG_ALL_SALT_DONE)? INCHI_FLAG_ACID_TAUT : 0;
        /*  number of tautomeric groups */
        pINChI->nTautomer[len ++] = (AT_NUMB)n;
        /* store each tautomeric group, one by one */
        for ( i = 0; i < n; i ++ ) {
            g = (int)t_group_info->tGroupNumber[i]; /* original group numbers in sorted order */
            t_group = t_group_info->t_group + g;    /* pointer to the tautomeric group */
            /*  NumAt+INCHI_T_NUM_MOVABLE (group length excluding this number) */
            pINChI->nTautomer[len ++]     = t_group->nNumEndpoints+INCHI_T_NUM_MOVABLE;
            /*  Num(H), Num(-) */
            for ( j = 0; j < INCHI_T_NUM_MOVABLE && j < T_NUM_NO_ISOTOPIC; j ++ )
                pINChI->nTautomer[len ++]     = t_group->num[j];
            for ( j = T_NUM_NO_ISOTOPIC; j < INCHI_T_NUM_MOVABLE; j ++ )
                pINChI->nTautomer[len ++]     = 0; /* should not happen */
            /* tautomeric group endpoint canonical numbers, pre-sorted in ascending order */
            for ( j  = (int)t_group->nFirstEndpointAtNoPos,
                  m  = j + (int)t_group->nNumEndpoints; j < m; j ++ ) {
                pINChI->nTautomer[len ++] = pCanonRank[(int)t_group_info->nEndpointAtomNumber[j]]; /*  At[j] */
            }
        }
        pINChI->lenTautomer = len;
        pINChI_Aux->nNumberOfTGroups = n;
    } else {
        pINChI->lenTautomer = 0;
        pINChI_Aux->nNumberOfTGroups = 0;
        if ( t_group_info && ((t_group_info->tni.bNormalizationFlags & FLAG_NORM_CONSIDER_TAUT) ||
                              t_group_info->nNumIsotopicEndpoints>1 &&
                              (t_group_info->bTautFlagsDone & (TG_FLAG_FOUND_ISOTOPIC_H_DONE | TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE)))
           ) {
            /* only protons (re)moved or added */
            pINChI->lenTautomer  = 1;
            pINChI->nTautomer[0] = 0;
        }
    }

    /*  number of H (excluding tautomeric) */
    if ( pCS->nNum_H ) {
        for ( i = 0; i < num_atoms; i ++ ) {
            pINChI->nNum_H[i] = pCS->nNum_H[i];
        } 
    }
    /*  number of fixed H (tautomeric H in non-tautomeric representation) */
    if ( pCS->nNum_H_fixed && !pINChI->lenTautomer ) {
        for ( i = 0; i < num_atoms; i ++ ) {
            pINChI->nNum_H_fixed[i]  = pCS->nNum_H_fixed[i];
            pINChI->nNum_H[i]       += pCS->nNum_H_fixed[i];
        } 
    }

    /***********************************************************
     *  tautomeric group(s) numbering and symmetry;
     *  should not depend on switching to rel. stereo numbering
     */
    if ( pINChI->lenTautomer && (n=pINChI_Aux->nNumberOfTGroups) ) {
        pCanonOrdTaut   = pCS->nLenCanonOrdStereoTaut > 0? pCS->nCanonOrdStereoTaut :
                          pCS->nLenCanonOrdTaut       > 0? pCS->nCanonOrdTaut : NULL;
        pConstitEquNumb = pINChI_Aux->nConstitEquTGroupNumbers;
        pSymmRank       = pCS->nSymmRankTaut;
        if ( pCanonOrdTaut && pSymmRank && pConstitEquNumb ) {
            for ( i = 0; i < n; i ++ ) {
                pConstitEquNumb[i]       = pSymmRank[pCanonOrdTaut[i]];
                pSortOrd[i]              = i;
            }
            pn_RankForSort  = pConstitEquNumb;
            qsort( pSortOrd, n, sizeof(pSortOrd[0]), CompRanksOrd );
            for ( i = 0, nMinOrd = pSortOrd[0], j = 1; j <= n; j ++ ) {
                if ( j == n || pConstitEquNumb[pSortOrd[i]] != pConstitEquNumb[pSortOrd[j]] ) {
                    nMinOrd ++; /* make is start from 1, not from zero */
                    if ( j - i > 1 ) {
                        /*  found a sequence of more than one equivalent t-groups: i..j-1 */
                        while ( i < j ) {
                            pConstitEquNumb[pSortOrd[i++]] = nMinOrd;
                        }
                    } else {
                        pConstitEquNumb[pSortOrd[i++]] = 0;
                    }
                    nMinOrd = pSortOrd[j]; /*  at the end j == n */
                }
            }
        }
    }

    /*  Allocate and fill Hill formula */
    if ( !(pINChI->szHillFormula = AllocateAndFillHillFormula( pINChI ) ) ) {
        nErrorCode = 0;
        ret = CT_WRONG_FORMULA; /* CT_OUT_OF_RAM;*/  /*   <BRKPT> */
        pINChI->nErrorCode = pINChI_Aux->nErrorCode = ret;
        goto exit_function;
    }

    if ( nStereoUnmarkMode = UnmarkAllUndefinedUnknownStereo( pINChI->Stereo, nUserMode ) ) {
        pINChI->nFlags |= (nStereoUnmarkMode & REQ_MODE_SC_IGN_ALL_UU)? INCHI_FLAG_SC_IGN_ALL_UU : 0;    
        pINChI->nFlags |= (nStereoUnmarkMode & REQ_MODE_SB_IGN_ALL_UU)? INCHI_FLAG_SB_IGN_ALL_UU : 0;
        if ( (nStereoUnmarkMode & REQ_MODE_SC_IGN_ALL_UU) ||
             (nStereoUnmarkMode & REQ_MODE_SB_IGN_ALL_UU) ) {
             AddMOLfileError(pStrErrStruct, "Omitted undefined stereo"); 
        }
    }

    /*************************/
    /* mark ambiguous stereo */
    /*************************/
    MarkAmbiguousStereo( at, norm_at, 0 /* non-isotopic */, pCanonOrd,
           pCS->LinearCTStereoCarb, pCS->nLenLinearCTStereoCarb,
           pCS->LinearCTStereoDble, pCS->nLenLinearCTStereoDble );


    /************************************************************************
     *
     *  isotopic part
     */
    /* abs or rel stereo may establish one of two canonical numberings */
    if ( (pCS->nLenLinearCTIsotopicStereoCarb > 0 || pCS->nLenLinearCTIsotopicStereoDble > 0) &&
          pCS->nLenCanonOrdIsotopicStereo > 0 &&
         (pCS->LinearCTIsotopicStereoCarb && pCS->LinearCTIsotopicStereoCarbInv ||
          pCS->LinearCTIsotopicStereoDble && pCS->LinearCTIsotopicStereoDbleInv) &&
          pCS->nCanonOrdIsotopicStereo    && pCS->nCanonOrdIsotopicStereoInv
          ) {
        /* found isotopic stereo */
        pCanonRank    = pCanonRankAtoms;
        pCanonOrd     = pCS->nCanonOrdIsotopicStereo;
        pCanonRankInv = pSortOrd;
        pCanonOrdInv  = pCS->nCanonOrdIsotopicStereoInv;
        Stereo        = pINChI->StereoIsotopic;
        for ( i = 0; i < num_at_tg; i ++ ) {
            pCanonRankInv[pCanonOrdInv[i]] =
            pCanonRank[pCanonOrd[i]]       = (AT_NUMB)(i+1);
        }
        /********************************************************************/
        /* copy stereo bonds and stereo centers; compare Inv and Abs stereo */
        /********************************************************************/
        nErrorCode = CopyLinearCTStereoToINChIStereo( Stereo,
                           pCS->LinearCTIsotopicStereoCarb, pCS->nLenLinearCTIsotopicStereoCarb,
                           pCS->LinearCTIsotopicStereoDble, pCS->nLenLinearCTIsotopicStereoDble
                           , pCanonOrd, pCanonRank, at, 1 /* isotopic */
                           , pCS->LinearCTIsotopicStereoCarbInv
                           , pCS->LinearCTIsotopicStereoDbleInv
                           , pCanonOrdInv, pCanonRankInv ); 

        if ( Stereo->t_parityInv && Stereo->nNumberInv ) {
            if ( nUserMode & REQ_MODE_RELATIVE_STEREO ) {
                pINChI->nFlags |= INCHI_FLAG_REL_STEREO;
            }
            if ( nUserMode & REQ_MODE_RACEMIC_STEREO ) {
                pINChI->nFlags |= INCHI_FLAG_RAC_STEREO;
            }
            if ( Stereo->nCompInv2Abs ) {
                if ( Stereo->nCompInv2Abs == -1 ) {
                    /* switch pointers so that the stereo becomes the smallest (relative)  */
                    /* flag Stereo->nCompInv2Abs == -1 will keep track of this exchange */
                    AT_NUMB    *nNumberInv  = Stereo->nNumberInv;
                    S_CHAR     *t_parityInv = Stereo->t_parityInv;
                    Stereo->nNumberInv  = Stereo->nNumber;
                    Stereo->t_parityInv = Stereo->t_parity;
                    Stereo->nNumber     = nNumberInv;
                    Stereo->t_parity    = t_parityInv;
                    switch_ptrs( &pCanonRank, &pCanonRankInv );
                    switch_ptrs( &pCanonOrd,  &pCanonOrdInv  );
                    bUseIsotopicNumberingInv    = 1;
                }
            }
        }

        for ( i = 0; i < num_atoms; i ++ ) {
            pINChI_Aux->nIsotopicOrigAtNosInCanonOrdInv[i] = at[pCanonOrdInv[i]].orig_at_number;
            pINChI_Aux->nIsotopicOrigAtNosInCanonOrd[i]    = at[pCanonOrd[i]].orig_at_number;
        }
        if ( bUseIsotopicNumberingInv ) {
            switch_ptrs( &pCanonRank, &pCanonRankInv );
            switch_ptrs( &pCanonOrd,  &pCanonOrdInv  );
            memcpy( pCanonRank, pCanonRankInv, num_at_tg * sizeof(pCanonRank[0]) );
            memcpy( pCanonOrd,  pCanonOrdInv, num_at_tg * sizeof(pCanonOrd[0]) );
        }
        pCanonRankInv = NULL;
        pCanonOrdInv  = NULL;
        pOrigNosInCanonOrd = NULL;

    } else {
        /* no isotopic stereo */
        pCanonOrd = pCS->nLenCanonOrdIsotopicStereo > 0? pCS->nCanonOrdIsotopicStereo :
                    pCS->nLenCanonOrdIsotopic       > 0? pCS->nCanonOrdIsotopic : NULL;
        pCanonRank           = pCanonRankAtoms;
        pOrigNosInCanonOrd   = pINChI_Aux->nIsotopicOrigAtNosInCanonOrd;
        if ( pCanonOrd && pCanonRank ) {
            for ( i = 0; i < num_atoms; i ++ ) { /* Fix13 -- out of bounds */
                pCanonRank[pCanonOrd[i]]       = (AT_NUMB)(i+1);
                pOrigNosInCanonOrd[i]          = at[pCanonOrd[i]].orig_at_number;
            }
            for ( ; i < num_at_tg; i ++ ) { /* Fix13 -- out of bounds */
                pCanonRank[pCanonOrd[i]]       = (AT_NUMB)(i+1);
            }
        }
    }
    /*pCanonOrdIso = pCanonOrd;*/

    pConstitEquNumb      = pINChI_Aux->nConstitEquIsotopicNumbers;
    pSymmRank            = pCS->nSymmRankIsotopic;
    if ( pCanonOrd && pCanonRank && pConstitEquNumb && pSymmRank ) {
        for ( i = 0; i < num_atoms; i ++ ) {
            pConstitEquNumb[i]       = pSymmRank[pCanonOrd[i]];
            pSortOrd[i]              = i;
        }
        for ( ; i < num_at_tg; i ++ ) {
            pSortOrd[i]              = i;
        }
        pn_RankForSort  = pConstitEquNumb;
        qsort( pSortOrd, num_atoms, sizeof(pSortOrd[0]), CompRanksOrd );
        for ( i = 0, nMinOrd = pSortOrd[0], j = 1; j <= num_atoms; j ++ ) {
            if ( j == num_atoms || pConstitEquNumb[pSortOrd[i]] != pConstitEquNumb[pSortOrd[j]] ) {
                nMinOrd ++;
                if ( j - i > 1 ) {
                    /*  found a sequence of equivalent atoms: i..j-1 */
                    while ( i < j ) {
                        pConstitEquNumb[pSortOrd[i++]] = nMinOrd;
                    }
                } else {
                    pConstitEquNumb[pSortOrd[i++]] = 0; /* nMinOrd; */
                }
                nMinOrd = pSortOrd[j];
            }
        }
    } else {
        goto exit_function; /*  no isotopic info available */
    }
    /*  isotopic atoms */
    n = pINChI->nNumberOfIsotopicAtoms = pCS->nLenLinearCTIsotopic;
    for ( i = 0; i < n; i ++ ) {
        pINChI->IsotopicAtom[i].nAtomNumber    = pCS->LinearCTIsotopic[i].at_num;
        pINChI->IsotopicAtom[i].nIsoDifference = pCS->LinearCTIsotopic[i].iso_atw_diff;
        pINChI->IsotopicAtom[i].nNum_H         = pCS->LinearCTIsotopic[i].num_1H;
        pINChI->IsotopicAtom[i].nNum_D         = pCS->LinearCTIsotopic[i].num_D;
        pINChI->IsotopicAtom[i].nNum_T         = pCS->LinearCTIsotopic[i].num_T;
    }
    /*  isotopic tautomeric groups */
    n = pINChI->nNumberOfIsotopicTGroups = pCS->nLenLinearCTIsotopicTautomer;
    for ( i = 0; i < n; i ++ ) {
        pINChI->IsotopicTGroup[i].nTGroupNumber = pCS->LinearCTIsotopicTautomer[i].tgroup_num;
        pINChI->IsotopicTGroup[i].nNum_H        = pCS->LinearCTIsotopicTautomer[i].num[2];
        pINChI->IsotopicTGroup[i].nNum_D        = pCS->LinearCTIsotopicTautomer[i].num[1];
        pINChI->IsotopicTGroup[i].nNum_T        = pCS->LinearCTIsotopicTautomer[i].num[0];
    }
    /* atoms that may exchange isotopic H-atoms */
    if ( pCS->nExchgIsoH && pINChI->nPossibleLocationsOfIsotopicH ) {
        for ( i = 0, j = 1; i < num_atoms; i ++ ) {
            if ( pCS->nExchgIsoH[i] ) {
                pINChI->nPossibleLocationsOfIsotopicH[j++] = (AT_NUMB)(i+1); /* canonical number */
            }
        }
        pINChI->nPossibleLocationsOfIsotopicH[0] = (AT_NUMB)j; /* length including the 0th element */
    }

    if ( nStereoUnmarkMode = UnmarkAllUndefinedUnknownStereo( pINChI->StereoIsotopic, nUserMode ) ) {
        pINChI->nFlags |= (nStereoUnmarkMode & REQ_MODE_SC_IGN_ALL_UU)? INCHI_FLAG_SC_IGN_ALL_ISO_UU : 0;    
        pINChI->nFlags |= (nStereoUnmarkMode & REQ_MODE_SB_IGN_ALL_UU)? INCHI_FLAG_SC_IGN_ALL_ISO_UU : 0;    
        if ( (nStereoUnmarkMode & REQ_MODE_SC_IGN_ALL_UU) ||
             (nStereoUnmarkMode & REQ_MODE_SB_IGN_ALL_UU) ) {
             AddMOLfileError(pStrErrStruct, "Omitted undefined stereo"); 
        }
    }
    /* mark ambiguous stereo */
    MarkAmbiguousStereo( at, norm_at, 1 /* isotopic */, pCanonOrd,
           pCS->LinearCTIsotopicStereoCarb, pCS->nLenLinearCTIsotopicStereoCarb,
           pCS->LinearCTIsotopicStereoDble, pCS->nLenLinearCTIsotopicStereoDble );

    /***********************************************************
     *  isotopic tautomeric group(s) numbering and symmetry;
     *  should not depend on switching to rel. stereo numbering
     */
    if ( pINChI->lenTautomer && pINChI_Aux->nConstitEquIsotopicTGroupNumbers && pCS->nSymmRankIsotopicTaut &&
         (pCS->nLenLinearCTIsotopic || pCS->nLenLinearCTIsotopicTautomer) &&
          t_group_info && t_group_info->num_t_groups > 0 ) {
        n = t_group_info->num_t_groups;
        pCanonOrdTaut   = pCS->nLenCanonOrdIsotopicStereoTaut > 0? 
                              (n=pCS->nLenCanonOrdIsotopicStereoTaut, pCS->nCanonOrdIsotopicStereoTaut) :
                          pCS->nLenCanonOrdIsotopicTaut       > 0?
                              (n=pCS->nLenCanonOrdIsotopicTaut,pCS->nCanonOrdIsotopicTaut) : (n=0,(AT_RANK*)NULL);
        pConstitEquNumb = pINChI_Aux->nConstitEquIsotopicTGroupNumbers;
        pSymmRank       = pCS->nSymmRankIsotopicTaut;
        if ( pCanonOrdTaut && pSymmRank && pConstitEquNumb && n > 0 ) {
            for ( i = 0; i < n; i ++ ) {
                pConstitEquNumb[i]       = pSymmRank[pCanonOrdTaut[i]];
                pSortOrd[i]              = i;
            }
            pn_RankForSort  = pConstitEquNumb;
            qsort( pSortOrd, n, sizeof(pSortOrd[0]), CompRanksOrd );
            for ( i = 0, nMinOrd = pSortOrd[0], j = 1; j <= n; j ++ ) {
                if ( j == n || pConstitEquNumb[pSortOrd[i]] != pConstitEquNumb[pSortOrd[j]] ) {
                    nMinOrd ++;
                    if ( j - i > 1 ) {
                        /*  found a sequence of equivalent t-groups: i..j-1 */
                        while ( i < j ) {
                            pConstitEquNumb[pSortOrd[i++]] = nMinOrd;
                        }
                    } else {
                        pConstitEquNumb[pSortOrd[i++]] = 0; /*  nMinOrd; */
                    }
                    nMinOrd = pSortOrd[j]; /*  at the end j = n */
                }
            }
        }
    }


exit_function:
    if ( pCanonRankAtoms )
        inchi_free( pCanonRankAtoms );
    if ( pSortOrd )
        inchi_free( pSortOrd );

    pINChI->nErrorCode     |= nErrorCode;
    pINChI_Aux->nErrorCode |= nErrorCode;

    return ret;
}

/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
void make_norm_atoms_from_inp_atoms(INCHIGEN_DATA *gendata, INCHIGEN_CONTROL *genctl)
{
/*^^^ TODO: make a full copy (with allocs) of atom arrays */
size_t t1;
int k;

    for ( k = 0;  k < INCHI_NUM; k++) 
    {
        if (NULL!=genctl->InpNormAtData[k])
        {
            t1 = genctl->StructData.num_components[k] * sizeof(NORM_ATOMS);
            memcpy(gendata->NormAtomsNontaut[k], genctl->InpNormAtData[k], t1);
        }
    
        if (NULL!=genctl->InpNormTautData[k])
        {
            t1 = genctl->StructData.num_components[k] * sizeof(NORM_ATOMS);
            memcpy(gendata->NormAtomsTaut[k], genctl->InpNormTautData[k], t1);
        }
    }

}

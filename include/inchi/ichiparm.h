/*
 * International Chemical Identifier (InChI)
 * Version 1
 * Software version 1.04
 * September 9, 2011
 *
 * The InChI library and programs are free software developed under the
 * auspices of the International Union of Pure and Applied Chemistry (IUPAC).
 * Originally developed at NIST. Modifications and additions by IUPAC 
 * and the InChI Trust.
 *
 * IUPAC/InChI-Trust Licence for the International Chemical Identifier (InChI) 
 * Software version 1.0.
 * Copyright (C) IUPAC and InChI Trust Limited
 * 
 * This library is free software; you can redistribute it and/or modify it under the 
 * terms of the IUPAC/InChI Trust Licence for the International Chemical Identifier 
 * (InChI) Software version 1.0; either version 1.0 of the License, or 
 * (at your option) any later version.
 * 
 * This library is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
 * See the IUPAC/InChI Trust Licence for the International Chemical Identifier (InChI) 
 * Software version 1.0 for more details.
 * 
 * You should have received a copy of the IUPAC/InChI Trust Licence for the 
 * International Chemical Identifier (InChI) Software version 1.0 along with 
 * this library; if not, write to:
 * 
 * The InChI Trust
 * c/o FIZ CHEMIE Berlin
 * Franklinstrasse 11
 * 10587 Berlin
 * GERMANY
 * 
 */


int DetectInputINChIFileType(FILE **inp_file, INPUT_PARMS *ip, const char *fmode);
void HelpCommandLineParmsReduced(INCHI_IOSTREAM *f);




/*****************************************************************************************/
int ReadCommandLineParms(int argc, const char *argv[], INPUT_PARMS *ip, 
                         char *szSdfDataValue, unsigned long *ulDisplTime, 
                         int bReleaseVersion, INCHI_IOSTREAM *log_file)
{
    int           i, k, c;
    const char    *q;
    unsigned long ul;
    int           nFontSize    = -9;
    int           nMode        = 0;
    int           nReleaseMode = nMode | (REQ_MODE_BASIC | REQ_MODE_TAUT | REQ_MODE_ISO | REQ_MODE_STEREO);


#if ( MIN_SB_RING_SIZE > 0 )
    int           nMinDbRinSize = MIN_SB_RING_SIZE, mdbr=0;
#endif

    /*int           bNotRecognized=0;*/
    char          szNameSuffix[32], szOutputPath[512];
    int           bNameSuffix, bOutputPath;
    int           bMergeAllInputStructures;
    int           bDisconnectSalts       = (DISCONNECT_SALTS==1);
    int           bDoNotAddH             = 0;



    int           bVer1Options           = 1;
    int           bReconnectCoord        = (RECONNECT_METALS==1);
    int           bDisconnectCoord       = (DISCONNECT_METALS==1);

    

#ifdef TARGET_LIB_FOR_WINCHI
/*	int           bVer1Options           = 0; 
    int           bReconnectCoord        = 1;
    int           bDisconnectCoord       = 1; */
    int			  is_gui				 = 1;
    int           bINChIOutputOptions    = INCHI_OUT_EMBED_REC; /* embed reconnected & output full aux info */
    int           bCompareComponents     = CMP_COMPONENTS;
#else
/*  int           bVer1Options           = 1;
    int           bReconnectCoord        = (RECONNECT_METALS==1);
    int           bDisconnectCoord       = (DISCONNECT_METALS==1); */
    int			  is_gui				 = 0;
    int           bINChIOutputOptions     = ((EMBED_REC_METALS_INCHI==1)? INCHI_OUT_EMBED_REC   : 0);
                                            /*| INCHI_OUT_NO_AUX_INFO INCHI_OUT_SHORT_AUX_INFO*/
    int           bCompareComponents     = 0;
#endif

    int           bDisconnectCoordChkVal = (CHECK_METAL_VALENCE == 1);
    int           bMovePositiveCharges   = (MOVE_CHARGES==1);
    int           bAcidTautomerism       = (DISCONNECT_SALTS==1)? (TEST_REMOVE_S_ATOMS==1? 2:1):0;
    int           bUnchargedAcidTaut     = (CHARGED_SALTS_ONLY==0);
    int           bMergeSaltTGroups      = (DISCONNECT_SALTS==1);
    int           bDisplayCompositeResults = 0;

#define VER103_DEFAULT_MODE    (REQ_MODE_TAUT | REQ_MODE_ISO | REQ_MODE_STEREO |\
                                REQ_MODE_SB_IGN_ALL_UU | REQ_MODE_SC_IGN_ALL_UU)

    INCHI_MODE     bVer1DefaultMode = VER103_DEFAULT_MODE;
   
    const char   *ext[MAX_NUM_PATHS+1];
    const char   *pArg;
    double        t;
    int           bRecognizedOption;
    int           bDisplay = 0;
    int           bNoStructLabels = 0;
    int           bOutputMolfileOnly = 0;
    int           bOutputMolfileDT = 0;
    int           bOutputMolfileSplit = 0;
    int           bForcedChiralFlag = 0;

#if ( READ_INCHI_STRING == 1 )
    int           bDisplayIfRestoreWarnings = 0;
#endif

#ifdef TARGET_LIB_FOR_WINCHI
    int           bXml = INCHI_OUT_XML;
#else
    int           bXml = INCHI_OUT_PLAIN_TEXT;
#endif
    int bTgFlagVariableProtons = 1;
    int bTgFlagHardAddRenProtons = 1;
#ifdef STEREO_WEDGE_ONLY
    int bPointedEdgeStereo = STEREO_WEDGE_ONLY; /*   NEWPS TG_FLAG_POINTED_EDGE_STEREO*/
#endif
#if ( FIX_ADJ_RAD == 1 )
    int bFixAdjacentRad = 0;
#endif
    int bAddPhosphineStereo = 1;
    int bAddArsineStereo    = 1;
    int bFixSp3bug          = 1; 
    int bFixFB2             = 1; 
    int bKetoEnolTaut       = 0;
    int b15TautNonRing      = 0;
    int bStdFormat          = 1;
    int bHashKey            = 0;
    int bHashXtra1          = 0;
    int bHashXtra2          = 0;

    ext[0] = ".mol";
    ext[1] = bVer1Options? ".txt" : ".ich";
    ext[2] = ".log";
    ext[3] = ".prb";
    ext[4] = "";

#if ( MAX_NUM_PATHS < 4 )
  #error Wrong initialization
#endif



    /*  init table parms */
    memset(ip, 0, sizeof(*ip));
    
    /* default are standard InChI generation options */
    bVer1DefaultMode   &=		~REQ_MODE_BASIC;	/* "FIXEDH - OFF" */
    bReconnectCoord     =		0;					/* "RECMET - OFF" */
    bPointedEdgeStereo  =		1;					/* "NEWPS" */
    ip->bFixNonUniformDraw =	1;					/* "FNUD" */
    


#ifndef COMPILE_ANSI_ONLY
    strcpy( ip->tdp.ReqShownFoundTxt[ilSHOWN], "Shown" );
    ip->dp.sdp.tdp = &ip->tdp;
    ip->dp.pdp     = &ip->pdp;
#endif

    memset( szNameSuffix, 0, sizeof(szNameSuffix) );
    bNameSuffix = 0;
    memset(szOutputPath, 0, sizeof(szOutputPath) );
    bOutputPath = 0;
    bMergeAllInputStructures = 0;

    *ulDisplTime  = 0;
#ifdef TARGET_API_LIB
    ip->msec_MaxTime = 0;      /*  milliseconds, default = unlimited in libinchi */
#else
    ip->msec_MaxTime = 60000;  /*  milliseconds, default = 60 sec */
#endif


    if ( bReleaseVersion ) 
    {
        /*  normal */
        /*ip->bINChIOutputOptions |= INCHI_OUT_PLAIN_TEXT;*/
        /*ip->bXml = 1;*/
        ip->bAbcNumbers = 0;
        ip->bCtPredecessors = 0;
        /*
         -- testing --
        ip->bXml = 0;
        ip->bAbcNumbers = 1;
        */
    } 
    else 
    {
        /*bXml = INCHI_OUT_PLAIN_TEXT;*/
        nReleaseMode = 0;
    }

    if ( bVer1Options ) 
    {
        bNameSuffix = 1;
        szNameSuffix[0] = '\0';
    }


    /* Parse */

    for ( i = 1; i < argc; i ++ ) 
    {
    
        /* if ( !(bVer1Options & 1) && INCHI_OPTION_PREFX == argv[i][0] && INCHI_OPTION_PREFX != argv[i][1] ) */
        if ( is_gui && INCHI_OPTION_PREFX == argv[i][0] && INCHI_OPTION_PREFX != argv[i][1] ) 
        {
        /*=== parsing TARGET_LIB_FOR_WINCHI GUI (and v. 0.9xx Beta)options ===*/

            pArg = argv[i]+1;
            

            /*--- Input options ---*/
            if ( !stricmp( pArg, "INPAUX" ))
            {
                if (INPUT_NONE == ip->nInputType) 
                    ip->nInputType = INPUT_INCHI_PLAIN;
            }



            else if ( INPUT_NONE == ip->nInputType &&
                    (!memicmp( pArg, "SDF", 3 )) &&
                    ( pArg[3] == ':' ) ) 
            {
                k = 0;
                mystrncpy( ip->szSdfDataHeader, pArg+4, MAX_SDF_HEADER+1 );
                LtrimRtrim( ip->szSdfDataHeader, &k );
                if ( k ) {
                    ip->pSdfLabel  = ip->szSdfDataHeader;
                    ip->pSdfValue  = szSdfDataValue;
                    ip->nInputType = INPUT_SDFILE;
                } else {
                    ip->pSdfLabel  = NULL;
                    ip->pSdfValue  = NULL;
                    ip->nInputType = INPUT_MOLFILE;
                }
            } 

            else if ( INPUT_NONE == ip->nInputType && !stricmp( pArg, "MOL" ) ) 
            {
                ip->nInputType = INPUT_MOLFILE;
            }             
            else if ( INPUT_NONE == ip->nInputType && !stricmp( pArg, "SDF" ) ) 
            {
                ip->nInputType = INPUT_MOLFILE;
            }             
#if ( ADD_CMLPP == 1 )
            else if ( INPUT_NONE == ip->nInputType && !stricmp( pArg, "CML" )  ) {
                 /* CMLfile label */
                ip->nInputType = INPUT_CMLFILE;
            } 
#endif            

            else if ( !memicmp( pArg, "START:", 6 ) ) 
            {
                ip->first_struct_number = strtol(pArg+6, NULL, 10);
            } 
            else if ( !memicmp( pArg, "END:", 4 ) ) 
            {
                ip->last_struct_number = strtol(pArg+4, NULL, 10);
            }
#ifdef BUILD_WITH_ENG_OPTIONS
            else if ( !memicmp( pArg, "RSB:", 4 ) ) 
            {
                mdbr = (int)strtol(pArg+4, NULL, 10);
            } else
            if ( !memicmp( pArg, "DISCONSALT:", 11 ) ) {
                bDisconnectSalts = (0 != strtol(pArg+11, NULL, 10));
            } else
            if ( !memicmp( pArg, "DISCONMETAL:", 12 ) ) {
                bDisconnectCoord = (0 != strtol(pArg+12, NULL, 10));
            } else
            if ( !memicmp( pArg, "RECONMETAL:", 11 ) ) {
                bReconnectCoord = (0 != strtol(pArg+11, NULL, 10));
            } else
            if ( !memicmp( pArg, "DISCONMETALCHKVAL:", 18 ) ) {
                bDisconnectCoordChkVal = (0 != strtol(pArg+18, NULL, 10));
            } else
            if ( !memicmp( pArg, "MOVEPOS:", 8 ) ) {
                bMovePositiveCharges = (0 != strtol(pArg+8, NULL, 10));
            } else
            if ( !memicmp( pArg, "MERGESALTTG:", 12 ) ) {
                bMergeSaltTGroups = (0 != strtol(pArg+12, NULL, 10));
            } else
            if ( !memicmp( pArg, "UNCHARGEDACIDS:", 15) ) {
                bUnchargedAcidTaut = (0 != strtol(pArg+15, NULL, 16));;
            } else
            if ( !memicmp( pArg, "ACIDTAUT:", 9 ) ) {
                bAcidTautomerism = c = (int)strtol(pArg+9, NULL, 10);
                if ( 0 <= c && c <= 2 )  bAcidTautomerism = c;
                /*else bNotRecognized = 2*bReleaseVersion;*/
            }
#endif 

            /*--- Output options ---*/
#if ( !defined(TARGET_API_LIB) && !defined(TARGET_LIB_FOR_WINCHI) )
            else if ( !stricmp( pArg, "Tabbed" ) ) 
            {
                bXml |=  INCHI_OUT_TABBED_OUTPUT;
            } 
#endif
            else if ( !stricmp( pArg, "NOLABELS" ) ) 
            {
                 bNoStructLabels = 1;
            } 
            else if ( !stricmp( pArg, "SAVEOPT" ) ) 
            {
                 bINChIOutputOptions |= INCHI_OUT_SAVEOPT;
            } 
            else if ( !stricmp( pArg, "AUXNONE" ) ) 
            {    
                /* no aux. info */
                bINChIOutputOptions |= INCHI_OUT_NO_AUX_INFO;  /* no aux info */
                bINChIOutputOptions &= ~INCHI_OUT_SHORT_AUX_INFO;
            } 


#if ( defined(BUILD_WITH_ENG_OPTIONS) || defined(TARGET_LIB_FOR_WINCHI) )
            else if ( !stricmp( pArg, "SDFID" ) ) 
            {
                ip->bGetSdfileId = 1;
            }             
            else if ( !stricmp( pArg, "XML" ) ) 
            {
                bXml &= ~INCHI_OUT_PLAIN_TEXT;
                bXml |=  INCHI_OUT_XML;
                /*bNotRecognized = 2*bReleaseVersion;*/
            }             
            else if ( !stricmp( pArg, "PLAIN" ) ) 
            {
                bXml |=  INCHI_OUT_PLAIN_TEXT;
                bXml &= ~INCHI_OUT_XML;
            } 
            else if ( !stricmp( pArg, "ANNPLAIN" ) ) 
            {
                bXml |=   INCHI_OUT_PLAIN_TEXT_COMMENTS;
                bXml &=  ~INCHI_OUT_XML_TEXT_COMMENTS;
            }             
            else if ( !stricmp( pArg, "ANNXML" ) ) 
            {
                bXml |=   INCHI_OUT_XML_TEXT_COMMENTS;
                bXml &=  ~INCHI_OUT_PLAIN_TEXT_COMMENTS;
            } 
            else if ( !memicmp( pArg, "AUXINFO:", 8 ) && isdigit(UCINT pArg[8]) ) 
            {
                k = strtol(pArg+8, NULL, 10);
                if ( k == 0 ) 
                {
                    bINChIOutputOptions |= INCHI_OUT_NO_AUX_INFO;  /* no aux info */
                    bINChIOutputOptions &= ~INCHI_OUT_SHORT_AUX_INFO;
                } 
                else if ( k == 1 ) 
                {
                    bINChIOutputOptions &= ~(INCHI_OUT_NO_AUX_INFO | INCHI_OUT_SHORT_AUX_INFO); /* include full aux info */
                } 
                else if ( k == 2 ) 
                {
                    bINChIOutputOptions &= ~INCHI_OUT_NO_AUX_INFO; /* include short aux info */
                    bINChIOutputOptions |= INCHI_OUT_SHORT_AUX_INFO;
                } 
                else 
                {
                    bINChIOutputOptions = k;  /* override everything */
                }
            }             
            else if ( !stricmp( pArg, "MERGE" ) ) 
            {
                bMergeAllInputStructures = 1;
            } 
            
            else if ( !stricmp( pArg, "PGO" ) ) 
            {
                ip->bSaveAllGoodStructsAsProblem = 1;
            }             
            else if ( !stricmp( pArg, "DCR" ) ) 
            {
                bDisplayCompositeResults = 1;
            } 
            
            else if ( !stricmp( pArg, "DSB" ) ) 
            {
                nMode |= REQ_MODE_NO_ALT_SBONDS;
            }             
            else if ( !stricmp( pArg, "NOHDR" ) ) 
            {
                 bNoStructLabels = 1;
            } 
            
            else if ( !stricmp( pArg, "NoVarH" ) ) 
            {
                 bTgFlagVariableProtons = 0;
            } 
#endif /* BUILD_WITH_ENG_OPTIONS */


            /*--- All modes (std and non-std InChI) structure perception options ---*/
            else if ( !stricmp( pArg, "SNON" ) ) 
            {
                bVer1DefaultMode &= ~REQ_MODE_STEREO; 
                nMode &= ~(REQ_MODE_RACEMIC_STEREO | REQ_MODE_RELATIVE_STEREO | REQ_MODE_CHIR_FLG_STEREO);
            } 
            else if ( !stricmp( pArg, "NEWPSOFF" ) ) 
            {
                 bPointedEdgeStereo = 0;
            } 
            else if ( !stricmp( pArg, "DONOTADDH" ) ) 
            {
                bDoNotAddH = 1;
            } 

            /* Non-standard */
#ifndef USE_STDINCHI_API
            
            /* Non-std InChI structure perception options */           
            else if ( !stricmp( pArg, "SREL" ) ) 
            {
                if ( nMode & REQ_MODE_RACEMIC_STEREO ) 
                {
                    nMode ^= REQ_MODE_RACEMIC_STEREO;
                }
                if ( nMode & REQ_MODE_CHIR_FLG_STEREO ) 
                {
                    nMode ^= REQ_MODE_CHIR_FLG_STEREO;
                }
                nMode |= REQ_MODE_RELATIVE_STEREO;
                nMode |= REQ_MODE_STEREO;
                bStdFormat = 0;
            }             
            else if ( !stricmp( pArg, "SRAC" ) ) 
            {
                if ( nMode & REQ_MODE_RELATIVE_STEREO ) 
                {
                    nMode ^= REQ_MODE_RELATIVE_STEREO;
                }
                if ( nMode & REQ_MODE_CHIR_FLG_STEREO ) 
                {
                    nMode ^= REQ_MODE_CHIR_FLG_STEREO;
                }
                nMode |= REQ_MODE_RACEMIC_STEREO;
                nMode |= REQ_MODE_STEREO;
                bStdFormat = 0;
            } 
            else if ( !stricmp( pArg, "SUCF" ) ) 
            {
                if ( nMode & REQ_MODE_RELATIVE_STEREO ) 
                {
                    nMode ^= REQ_MODE_RELATIVE_STEREO;
                }
                if ( nMode & REQ_MODE_RACEMIC_STEREO ) 
                {
                    nMode ^= REQ_MODE_RACEMIC_STEREO;
                }
                nMode |= REQ_MODE_CHIR_FLG_STEREO; /* stereo defined by the Chiral flag */
                nMode |= REQ_MODE_STEREO;
                bStdFormat = 0;
            }             
            else if ( !stricmp( pArg, "ChiralFlagON" ) ) 
            { 
                /* used only with /SUCF */
                /* NB: do not toggle off bStdFormat! (if necessary SUCF will do) */
                bForcedChiralFlag &= ~FLAG_SET_INP_AT_NONCHIRAL;
                bForcedChiralFlag |= FLAG_SET_INP_AT_CHIRAL;
            }             
            else if ( !stricmp( pArg, "ChiralFlagOFF" ) ) 
            { 
                /* used only with /SUCF */
                /* NB: do not toggle off bStdFormat! (if necessary SUCF will do) */
                bForcedChiralFlag &= ~FLAG_SET_INP_AT_CHIRAL;
                bForcedChiralFlag |= FLAG_SET_INP_AT_NONCHIRAL;
            }             

            /*--- Non-std InChI creation options ---*/

            else if ( !stricmp( pArg, "SUU" ) ) 
            {       
                /* include omitted undef/unknown stereo */
                bVer1DefaultMode &= ~(REQ_MODE_SB_IGN_ALL_UU | REQ_MODE_SC_IGN_ALL_UU);
                bStdFormat = 0;
            } 
            else if ( !stricmp( pArg, "SLUUD" ) ) 
            {       
                /* make labels for unknown and undefined stereo different */
                bVer1DefaultMode |= REQ_MODE_DIFF_UU_STEREO;
                bStdFormat = 0;
            } 
            /* FixedH */
            else if ( !stricmp( pArg, "FIXEDH" ) ) 
            {
                bVer1DefaultMode |= REQ_MODE_BASIC;  /* tautomeric */
                bStdFormat = 0;
            }             
            /* RecMet */
            else if ( !stricmp( pArg, "RECMET" ) ) 
            {    
                /* reconnect metals */
                bReconnectCoord = 1;
                bStdFormat = 0;
            } 
#if ( KETO_ENOL_TAUT == 1 )
            else if ( !stricmp( pArg, "KET" ) ) 
            {
                bKetoEnolTaut = 1;
                bStdFormat = 0;
            } 
#endif
#if ( TAUT_15_NON_RING == 1 )
            else if ( !stricmp( pArg, "15T" ) ) 
            {
                b15TautNonRing = 1;
                bStdFormat = 0;
            } 
#endif

#endif
            /*--- Generation options ---*/

            /* InChIKey/InChI hash */
            else if ( !stricmp( pArg, "Key" ) ) 
            {
                bHashKey = 1;
            }            
            else if ( !stricmp( pArg, "XHash1" ) ) 
            {
                bHashXtra1 = 1;
            }            
            else if ( !stricmp( pArg, "XHash2" ) ) 
            {
                bHashXtra2 = 1;
            }            


            /*--- (engineering) Undo bug/draw fixes options ---*/
#ifdef BUILD_WITH_ENG_OPTIONS
            else if ( !stricmp( pArg, "FixSp3bugOFF" ) ) 
            {
                 bFixSp3bug = 0;
                 bStdFormat = 0;
            }             
            else if ( !stricmp( pArg, "FBOFF" ) ) 
            {
                 bFixSp3bug = 0;
                 bStdFormat = 0;
            }             
            else  if ( !stricmp( pArg, "FB2OFF" ) ) 
            {
                 bFixFB2 = 0; 
                 bStdFormat = 0;
            }             
            else if ( !stricmp( pArg, "SPXYZOFF" ) ) 
            {
                 bAddPhosphineStereo = 0;
                 bStdFormat = 0;
            } 
            else if ( !stricmp( pArg, "SASXYZOFF" ) ) 
            {
                 bAddArsineStereo = 0;
                 bStdFormat = 0;
            } 
            else if ( !stricmp( pArg, "FNUDOFF" ) ) 
            {
                 ip->bFixNonUniformDraw = 0;
                 bStdFormat = 0;
            }            

            /*--- (hidden) Old structure-perception and InChI creation options ---*/
            /*--- (engineering) Old structure-perception and InChI creation options ---*/

            else if ( !stricmp( pArg, "NOUUSB" ) ) 
            {
                nMode |= REQ_MODE_SB_IGN_ALL_UU;
                bStdFormat = 0;
            }             
            else if ( !stricmp( pArg, "NOUUSC" ) ) 
            {
                nMode |= REQ_MODE_SC_IGN_ALL_UU;
                bStdFormat = 0;
            } 

#if ( FIX_ADJ_RAD == 1 )
            else if ( !stricmp( pArg, "FixRad" ) ) 
            {
                bFixAdjacentRad = 1;
                bStdFormat = 0;
            } 
#endif
#if ( UNDERIVATIZE == 1 )
            else if ( !stricmp( pArg, "DoDRV" ) ) 
            {
                ip->bUnderivatize = 1;
                bStdFormat = 0;
            } 
#endif
#if ( RING2CHAIN == 1 )
            else if ( !stricmp( pArg, "DoR2C" ) ) 
            {
                ip->bRing2Chain = 1;
                bStdFormat = 0;
            } 
#endif
#if ( RING2CHAIN == 1 || UNDERIVATIZE == 1 )
            else if ( !stricmp( pArg, "DoneOnly" ) ) 
            {
                ip->bIngnoreUnchanged = 1;
                bStdFormat = 0;
            } 
#endif
            else if ( !stricmp( pArg, "NoADP" ) ) 
            {
                 bTgFlagHardAddRenProtons = 0;
                 bStdFormat = 0;
            } 
            else if ( !memicmp( pArg, "DISCONSALT:", 11 ) ) 
            {
                bDisconnectSalts = (0 != strtol(pArg+11, NULL, 10));
                bStdFormat = 0;
            }             
            else if ( !memicmp( pArg, "DISCONMETAL:", 12 ) ) 
            {
                bDisconnectCoord = (0 != strtol(pArg+12, NULL, 10));
                bStdFormat = 0;
            }             
            else if ( !memicmp( pArg, "RECONMETAL:", 11 ) ) 
            {
                bReconnectCoord = (0 != strtol(pArg+11, NULL, 10));
                bStdFormat = 0;
            }             
            else if ( !memicmp( pArg, "DISCONMETALCHKVAL:", 18 ) ) 
            {
                bDisconnectCoordChkVal = (0 != strtol(pArg+18, NULL, 10));
                bStdFormat = 0;
            }             
            else if ( !memicmp( pArg, "MOVEPOS:", 8 ) ) 
            {
                bMovePositiveCharges = (0 != strtol(pArg+8, NULL, 10));
                bStdFormat = 0;
            }             
            else if ( !memicmp( pArg, "MERGESALTTG:", 12 ) ) 
            {
                bMergeSaltTGroups = (0 != strtol(pArg+12, NULL, 10));
                bStdFormat = 0;
            }             
            else if ( !memicmp( pArg, "UNCHARGEDACIDS:", 15) ) 
            {
                bUnchargedAcidTaut = (0 != strtol(pArg+15, NULL, 16));
                bStdFormat = 0;
            }             
            else if ( !memicmp( pArg, "ACIDTAUT:", 9 ) ) 
            {
                bAcidTautomerism = c = (int)strtol(pArg+9, NULL, 10);
                if ( 0 <= c && c <= 2 )  bAcidTautomerism = c;
                /*else bNotRecognized = 2*bReleaseVersion;*/
                bStdFormat = 0;
            }             

            /*--- (hidden) Old output and other options ---*/

            else if ( !memicmp( pArg, "O:", 2 ) ) 
            {
                bNameSuffix = 1;
                strncpy(szNameSuffix, pArg+2, sizeof(szNameSuffix)-1);
            }             
            else if ( !memicmp( pArg, "OP:", 3 ) ) 
            {
                bOutputPath = 1;
                strncpy(szOutputPath, pArg+3, sizeof(szOutputPath)-1);
            }             
            else if ( !stricmp( pArg, "ALT" ) ) 
            {
                ip->bAbcNumbers = 1;
                bStdFormat = 0;
            }             
            else if ( !stricmp( pArg, "SCT" ) ) 
            {
                ip->bCtPredecessors = 1;
                bStdFormat = 0;
            }             
            else if ( !stricmp( pArg, "CMP" ) ) 
            {
                bCompareComponents = CMP_COMPONENTS;
            }             
            else if ( !stricmp( pArg, "CMPNONISO" ) ) 
            {
                bCompareComponents = CMP_COMPONENTS | CMP_COMPONENTS_NONISO;
            } 
            else if ( !stricmp( pArg, "PW" ) ) 
            {
                ip->bSaveWarningStructsAsProblem = 1;
            }

#endif /* BUILD_WITH_ENG_OPTIONS */


            else 
            {
                /*for ( k = 0; c=pArg[k]; k ++ )*/
                k = 0;
                c=pArg[k]; /* prohibit multiple option concatenations, strict syntax check 2008-11-05 DT  */ 
                {
                    c = toupper( c );
                    switch ( c ) 
                    {
                    case 'D':
                        bDisplay |= 1;
                        if ( (pArg[k+1] == 'C' || pArg[k+1] == 'c') && !pArg[k+2] ) 
                        {
                            bDisplay |= 1;
                            k++;
                            ip->bDisplayEachComponentINChI = 1;
                        } 
                        else if ( !pArg[k+1] ) 
                        {
                            bDisplay |= 1;
                        }
                        break;
                    case 'W':
                        if ( pArg[k+1] == 'D' ) 
                        {  
                            /* restore Display Time functionality */
                            c = 'D';
                            k ++;
                        }
                        t = strtod( pArg+k+1, (char**)&q ); /*  cast deliberately discards 'const' qualifier */
                        if ( q > pArg+k+1 && errno == ERANGE || t < 0.0 || t*1000.0 > (double)ULONG_MAX) 
                        {
                            ul = 0;
                        } 
                        else 
                        {
                            ul = (unsigned long)(t*1000.0);
                        }
                        if ( /*q > pArg+k &&*/ !*q ) 
                        {
                            k = q - pArg - 1; /* k will be incremented by the for() cycle */
                            switch( c ) 
                            {
                                case 'D':
                                    *ulDisplTime = ul;
                                    break;
                                case 'W':
                                    ip->msec_MaxTime = ul;
                                    break;
                            }
                        }
                        break;
                    case 'F':
                        c =  (int)strtol( pArg+k+1, (char**)&q, 10 ); /*  cast deliberately discards 'const' qualifier */
                        if ( q > pArg+k && !*q ) 
                        {
                            k = q - pArg - 1;
                            if ( abs(c) > 5 ) 
                            {
                                nFontSize = -c;  /* font size 5 or less is too small */
                            }
                        }
                        break;
                    default:
                        if ( !pArg[k+1] ) 
                        {
                            switch ( c ) 
                            {
                            case 'B':
                                nMode |= REQ_MODE_BASIC;
                                nReleaseMode = 0;
                                /*bNotRecognized = bReleaseVersion;*/
                                bStdFormat = 0;
                                break;
                            case 'T':
                                nMode |= REQ_MODE_TAUT;
                                nReleaseMode = 0;
                                /*bNotRecognized = bReleaseVersion;*/
                                break;
                            case 'I':
                                nMode |= REQ_MODE_ISO;
                                nReleaseMode = 0;
                                /*bNotRecognized = bReleaseVersion;*/
                                break;
                            case 'N':
                                nMode |= REQ_MODE_NON_ISO;
                                nReleaseMode = 0;
                                bStdFormat = 0;
                                /*bNotRecognized = bReleaseVersion;*/
                                break;
                            case 'S':
                                nMode |= REQ_MODE_STEREO;
                                nReleaseMode = 0;
                                /*bNotRecognized = bReleaseVersion;*/
                                break;
                            case 'E':
                                if ( nReleaseMode & REQ_MODE_STEREO ) 
                                {
                                    nReleaseMode ^= REQ_MODE_STEREO;
                                    bStdFormat = 0;
                                }
                                break;
#ifndef TARGET_LIB_FOR_WINCHI
                            default:
                                inchi_ios_eprint(log_file, "Unrecognized option: \"%c\".\n", c);

#endif
                            }
                        } 


#ifndef TARGET_LIB_FOR_WINCHI
                        else 
                        {
                            inchi_ios_eprint(log_file, "Unrecognized option: \"%c\".\n", c);
                        }
#endif
                    }
                    /*
                    if ( bNotRecognized && bNotRecognized == bReleaseVersion ) {
                        inchi_ios_eprint(stderr, "Unrecognized option: \"%c\".\n", c);
                        bNotRecognized = 0;
                    }
                    */
                }
            }


            /*
            if ( bNotRecognized && bNotRecognized == 2*bReleaseVersion ) 
            {
               inchi_ios_eprint(stderr, "Unrecognized option: \"%s\".\n", argv[i]);
               bNotRecognized = 0;
            }
            */

        } 
        /*=== end of parsing TARGET_LIB_FOR_WINCHI GUI (and v. 0.9xx Beta)options ===*/


                

        else if ( (bVer1Options & 1) && INCHI_OPTION_PREFX == argv[i][0] && argv[i][1] ) 
        {
        /*=== parsing stand-alone/library InChI options ===*/


            pArg = argv[i] + 1;

            bRecognizedOption = 2;
            bVer1Options += 2;
            /* always on: REQ_MODE_TAUT | REQ_MODE_ISO | REQ_MODE_STEREO */
#ifdef CML_DEBUG
            printf ("1 argv %d %s\n", i, argv[i]);
#endif

            
            /*--- Input options ---*/
            if ( !stricmp( pArg, "STDIO" ) ) 
            {
                bNameSuffix = 0;
            }             
            else if ( !stricmp( pArg, "INPAUX" ))
            {
                if (INPUT_NONE == ip->nInputType) 
                    ip->nInputType = INPUT_INCHI_PLAIN;
            }             
            else if ( /* INPUT_NONE == ip->nInputType &&*/
                    !memicmp( pArg, "SDF:", 4 )  ) 
            {
                 /* SDfile label */
                k = 0;
                mystrncpy( ip->szSdfDataHeader, pArg+4, MAX_SDF_HEADER+1 );
                LtrimRtrim( ip->szSdfDataHeader, &k );
                if ( k ) 
                {
                    ip->pSdfLabel  = ip->szSdfDataHeader;
                    ip->pSdfValue  = szSdfDataValue;
                    if ( INPUT_NONE == ip->nInputType ) 
                    {
                        ip->nInputType = INPUT_SDFILE;
                    }
                } 
                else 
                {
                    ip->pSdfLabel  = NULL;
                    ip->pSdfValue  = NULL;
                    if ( INPUT_NONE == ip->nInputType ) 
                    {
                        ip->nInputType = INPUT_MOLFILE;
                    }
                }
            }

#if ( ADD_CMLPP == 1 )
            else if ( INPUT_NONE == ip->nInputType && !stricmp( pArg, "CML" )  ) 
            {
                /* CMLfile label */
                ip->nInputType = INPUT_CMLFILE;
            } 
#endif
            else if ( !memicmp( pArg, "START:", 6 ) ) 
            {
                ip->first_struct_number = strtol(pArg+6, NULL, 10);
            }             
            else if ( !memicmp( pArg, "END:", 4 ) ) 
            {
                ip->last_struct_number = strtol(pArg+4, NULL, 10);
            }             
#ifdef BUILD_WITH_ENG_OPTIONS
            else if ( !memicmp( pArg, "RSB:", 4 )) 
            {
                mdbr = (int)strtol(pArg+4, NULL, 10);
            } 
#endif


            /*--- Output options ---*/
#if ( !defined(TARGET_API_LIB) && !defined(TARGET_LIB_FOR_WINCHI) )
            else if ( !stricmp( pArg, "Tabbed" ) ) 
            {
                bXml |=  INCHI_OUT_TABBED_OUTPUT;
            } 
#endif
            else if ( !stricmp( pArg, "NOLABELS" ) ) 
            {
                 bNoStructLabels = 1;
            } 
            else if ( !stricmp( pArg, "SAVEOPT" ) ) 
            {
                 bINChIOutputOptions |= INCHI_OUT_SAVEOPT;
            } 
            else if ( !stricmp( pArg, "AUXNONE" ) ) 
            {    
                /* no aux. info */
                bINChIOutputOptions |= INCHI_OUT_NO_AUX_INFO;  /* no aux info */
                bINChIOutputOptions &= ~INCHI_OUT_SHORT_AUX_INFO;
            } 
            else if ( !stricmp( pArg, "OUTPUTSDF" ) ) 
            {  
                /* output SDfile */
                bOutputMolfileOnly = 1;
            }             
            else if ( !stricmp( pArg, "SdfAtomsDT" ) ) 
            {  
                /* output isotopes H as D and T in SDfile */
                bOutputMolfileDT = 1;
            } 
            else if ( !stricmp( pArg, "D" ) ) 
            {          
                /* display the structures */
                bDisplay |= 1;
            } 
            else if ( !memicmp( pArg, "F", 1 ) && (c =  (int)strtol( pArg+1, (char**)&q, 10 ), q > pArg+1) ) 
            {
                nFontSize = -c;                      /* struct. display font size */
            } 
            else if ( !stricmp( pArg, "EQU" ) ) 
            {
                bCompareComponents = CMP_COMPONENTS;
            } 


            /*--- All modes (std and non-std InChI) structure perception options ---*/
            else if ( !stricmp( pArg, "SNON" ) ) 
            {
                bVer1DefaultMode &= ~REQ_MODE_STEREO; /* no stereo */
                nMode &= ~(REQ_MODE_RACEMIC_STEREO | REQ_MODE_RELATIVE_STEREO | REQ_MODE_CHIR_FLG_STEREO);
            } 
            else if ( !stricmp( pArg, "NEWPSOFF" ) ) 
            {
                 bPointedEdgeStereo = 0;
            } 
            else if ( !stricmp( pArg, "DONOTADDH" ) ) 
            {
                bDoNotAddH = 1;
            }


            /* Non-standard */
#ifndef USE_STDINCHI_API

            /* Non-std InChI structure perception options */

            else if ( !stricmp( pArg, "SREL" ) ) 
            {
                bVer1DefaultMode |= REQ_MODE_STEREO;  /* relative stereo */
                nMode &= ~(REQ_MODE_RACEMIC_STEREO | REQ_MODE_CHIR_FLG_STEREO);
                nMode |= REQ_MODE_RELATIVE_STEREO;
                bStdFormat = 0;
            }             
            else  if ( !stricmp( pArg, "SRAC" ) ) 
            {
                /* REQ_MODE_CHIR_FLG_STEREO */
                bVer1DefaultMode |= REQ_MODE_STEREO;  /* racemic stereo */
                nMode &= ~(REQ_MODE_RELATIVE_STEREO | REQ_MODE_CHIR_FLG_STEREO);
                nMode |= REQ_MODE_RACEMIC_STEREO;
                bStdFormat = 0;
            }             
            else if ( !stricmp( pArg, "SUCF" ) ) 
            {
                bVer1DefaultMode |= REQ_MODE_STEREO;  /* stereo defined by the Chiral flag */
                nMode &= ~(REQ_MODE_RELATIVE_STEREO | REQ_MODE_RACEMIC_STEREO);
                nMode |= REQ_MODE_CHIR_FLG_STEREO;
                bStdFormat = 0;
            }             
            else if ( !stricmp( pArg, "ChiralFlagON" ) ) 
            { 
                /* used only with /SUCF */
                /* NB: do not toggle off bStdFormat! (if necessary SUCF will do) */
                bForcedChiralFlag &= ~FLAG_SET_INP_AT_NONCHIRAL;
                bForcedChiralFlag |= FLAG_SET_INP_AT_CHIRAL;
            }             
            else if ( !stricmp( pArg, "ChiralFlagOFF" ) ) 
            { 
                /* used only with /SUCF */
                /* NB: do not toggle off bStdFormat! (if necessary SUCF will do) */
                bForcedChiralFlag &= ~FLAG_SET_INP_AT_CHIRAL;
                bForcedChiralFlag |= FLAG_SET_INP_AT_NONCHIRAL;
            }             


            /*--- Non-std InChI creation options ---*/

            /* Stereo */
            else if ( !stricmp( pArg, "SUU" ) ) 
            {       
                /* include omitted undef/unknown stereo */
                bVer1DefaultMode &= ~(REQ_MODE_SB_IGN_ALL_UU | REQ_MODE_SC_IGN_ALL_UU);
                bStdFormat = 0;
            } 
            else if ( !stricmp( pArg, "SLUUD" ) ) 
            {       
               /* Make labels for unknown and undefined stereo different */
                bVer1DefaultMode |= REQ_MODE_DIFF_UU_STEREO;
                bStdFormat = 0;
            } 
            /* FixedH */
            else if ( !stricmp( pArg, "FIXEDH" ) ) 
            {
                bVer1DefaultMode |= REQ_MODE_BASIC;  /* tautomeric */
                bStdFormat = 0;
            }             
            /* RecMet */
            else if ( !stricmp( pArg, "RECMET" ) ) 
            {    
                /* reconnect metals */
                bReconnectCoord = 1;
                bStdFormat = 0;
            } 
#if ( KETO_ENOL_TAUT == 1 )
            else if ( !stricmp( pArg, "KET" ) ) 
            {
                bKetoEnolTaut = 1;
                bStdFormat = 0;
            } 
#endif
#if ( TAUT_15_NON_RING == 1 )
            else if ( !stricmp( pArg, "15T" ) ) 
            {
                b15TautNonRing = 1;
                bStdFormat = 0;
            } 
#endif

#endif

            /*--- Generation options ---*/
            else if ( !memicmp( pArg, "W", 1 ) && (t = strtod( pArg+1, (char**)&q ), q > pArg+1) ) 
            {
                if ( errno == ERANGE || t < 0.0 || t*1000.0 > (double)ULONG_MAX)  
                {
                    ul = 0;
                } 
                else 
                {
                    ul = (unsigned long)(t*1000.0);  /* max. time per structure */
                }
                ip->msec_MaxTime = ul;
            }             
            else if ( !stricmp( pArg, "WarnOnEmptyStructure" ) ) 
            {
                 ip->bAllowEmptyStructure = 1;
            } 

            /* InChIKey/InChI hash */
            else if ( !stricmp( pArg, "Key" ) ) 
            {
                bHashKey = 1;
            }            
            else if ( !stricmp( pArg, "XHash1" ) ) 
            {
                bHashXtra1 = 1;
            }            
            else if ( !stricmp( pArg, "XHash2" ) ) 
            {
                bHashXtra2 = 1;
            }            


            /*--- Conversion modes ---*/           


#if ( READ_INCHI_STRING == 1 )
            else if ( !stricmp( pArg, "InChI2InChI" )  ) 
            {
                 /* Read InChI Identifiers and output InChI Identifiers */
                ip->nInputType = INPUT_INCHI;
                ip->bReadInChIOptions |= READ_INCHI_OUTPUT_INCHI;
                ip->bReadInChIOptions &= ~READ_INCHI_TO_STRUCTURE;
            } 
            else if ( !stricmp( pArg, "InChI2Struct" )  ) 
            {
                 /* Split InChI Identifiers into components */
                ip->bReadInChIOptions |= READ_INCHI_TO_STRUCTURE;
                ip->bReadInChIOptions &= ~READ_INCHI_OUTPUT_INCHI;
                ip->nInputType = INPUT_INCHI;
            } 
#ifdef BUILD_WITH_ENG_OPTIONS
            else if ( !stricmp( pArg, "KeepBalanceP" )  ) 
            {
                 /* When spliting InChI Identifiers into components: */
                 /* If MobileH output then add p to each component;  */
                 /* Otherwise add one more component containing balance */
                 /* of protons and exchangeable isotopic H */
                ip->bReadInChIOptions |= READ_INCHI_KEEP_BALANCE_P;
                bStdFormat = 0;
            } 
#endif
#endif

          
            
            /*--- (engineering) Undo bug/draw fixes options ---*/
#ifdef BUILD_WITH_ENG_OPTIONS
            else if ( !stricmp( pArg, "FixSp3bugOFF" ) ) 
            {
                 bFixSp3bug = 0;
                 bStdFormat = 0;
            }             
            else if ( !stricmp( pArg, "FBOFF" ) ) 
            {
                 bFixSp3bug = 0;
                 bStdFormat = 0;
            }             
            else if ( !stricmp( pArg, "FB2OFF" ) ) 
            {
                 bFixFB2 = 0;
                 bStdFormat = 0;
            }             
            else if ( !stricmp( pArg, "SPXYZOFF" ) ) 
            {
                 bAddPhosphineStereo = 0;
                 bStdFormat = 0;
            } 
            else if ( !stricmp( pArg, "SASXYZOFF" ) ) 
            {
                 bAddArsineStereo = 0;
                 bStdFormat = 0;
            } 
            else if ( !stricmp( pArg, "FNUDOFF" ) ) 
            {
                 ip->bFixNonUniformDraw = 0;
                 bStdFormat = 0;
            } 

            /*--- (engineering) Old structure-perception and InChI creation options ---*/
#if ( FIX_ADJ_RAD == 1 )
            else if ( !stricmp( pArg, "FixRad" ) ) 
            {
                bFixAdjacentRad = 1;
                bStdFormat = 0;
            } 
#endif
#if ( UNDERIVATIZE == 1 )
            else if ( !stricmp( pArg, "DoDRV" ) ) 
            {
                ip->bUnderivatize = 1;
                bStdFormat = 0;
            } 
#endif
#if ( RING2CHAIN == 1 )
            else if ( !stricmp( pArg, "DoR2C" ) ) 
            {
                ip->bRing2Chain = 1;
                bStdFormat = 0;
            } 
#endif
#if ( RING2CHAIN == 1 || UNDERIVATIZE == 1 )
            else if ( !stricmp( pArg, "DoneOnly" ) ) 
            {
                ip->bIngnoreUnchanged = 1;
                bStdFormat = 0;
            } 
#endif

            else if ( !memicmp( pArg, "MOVEPOS:", 8 ) ) 
            {
                bMovePositiveCharges = (0 != strtol(pArg+8, NULL, 10));
                bStdFormat = 0;
            }             

            else if ( !stricmp( pArg, "NoADP" ) ) 
            {
                 bTgFlagHardAddRenProtons = 0;
                 bStdFormat = 0;
            } 
            /* Tautomer perception off */
            else if ( !stricmp( pArg, "EXACT" ) ) 
            {
                bVer1DefaultMode |= REQ_MODE_BASIC;
                bStdFormat = 0;
            } 
            else if ( !stricmp( pArg, "ONLYRECSALT" ) ) 
            {  
                /* do not disconnect salts */
                bDisconnectSalts = 0;
                bStdFormat = 0;
            } 
            else if ( !stricmp( pArg, "ONLYEXACT" ) || !stricmp( pArg, "ONLYFIXEDH" ) ) 
            {
                bVer1DefaultMode |=  REQ_MODE_BASIC;
                bVer1DefaultMode &= ~REQ_MODE_TAUT;
                bStdFormat = 0;
            }             
            else if ( !stricmp( pArg, "ONLYNONISO" ) ) 
            {
                bVer1DefaultMode |=  REQ_MODE_NON_ISO;
                bVer1DefaultMode &= ~REQ_MODE_ISO;
                bStdFormat = 0;
            }             
            else if ( !stricmp( pArg, "TAUT" ) ) 
            {
                bVer1DefaultMode &= ~REQ_MODE_BASIC;
                bVer1DefaultMode |= REQ_MODE_TAUT;
            }             
            else if ( !stricmp( pArg, "ONLYRECMET" ) ) 
            {  
                /* do not disconnect metals */
                bDisconnectCoord = 0;
                bStdFormat = 0;
            }

            /*--- (hidden) Old output and other options ---*/

            else if ( !stricmp( pArg, "SdfSplit" ) ) 
            {  
                /* Split single Molfiles into disconnected components */
                bOutputMolfileSplit = 1;
            } 
            else if ( !stricmp( pArg, "DCR" ) ) 
            {
                bDisplayCompositeResults = 1;
            } 
            else if ( !stricmp( pArg, "AUXFULL" ) || !stricmp( pArg, "AUXMAX" ) ) 
            {     
                /* full aux info */
                bINChIOutputOptions &= ~(INCHI_OUT_NO_AUX_INFO | INCHI_OUT_SHORT_AUX_INFO); /* include short aux info */
            }             
            else if ( !stricmp( pArg, "AUXMIN" ) ) 
            {     
                /* minimal aux info */
                bINChIOutputOptions &= ~INCHI_OUT_NO_AUX_INFO; /* include short aux info */
                bINChIOutputOptions |= INCHI_OUT_SHORT_AUX_INFO;
            } 
#if ( READ_INCHI_STRING == 1 )
            else if ( !stricmp( pArg, "DDSRC" ) ) 
            {
                bDisplayIfRestoreWarnings = 1;  /* InChI->Structure debugging: Display Debug Structure Restore Components */
            }
#endif
            else if ( !stricmp( pArg, "NoVarH" ) ) 
            {
                 bTgFlagVariableProtons = 0;
            }             
            else if ( !stricmp( pArg, "FULL" ) ) 
            {
                bVer1DefaultMode       = VER103_DEFAULT_MODE;
                nMode                = 0;
                bReconnectCoord      = 1;            /* full output */
                bINChIOutputOptions   = ((EMBED_REC_METALS_INCHI==1)? INCHI_OUT_EMBED_REC   : 0) | INCHI_OUT_SHORT_AUX_INFO;
                ip->bCtPredecessors  = 0;
                ip->bAbcNumbers      = 0;
                bXml                 |=  INCHI_OUT_PLAIN_TEXT | INCHI_OUT_PLAIN_TEXT_COMMENTS;
                bXml                 &= ~(INCHI_OUT_XML | INCHI_OUT_XML_TEXT_COMMENTS);
            }             
            else if ( !stricmp( pArg, "MIN" ) ) 
            {
                bVer1DefaultMode     = VER103_DEFAULT_MODE;
                nMode                = 0;
                bReconnectCoord      = 1;            /* minimal output */
                bINChIOutputOptions   = ((EMBED_REC_METALS_INCHI==1)? INCHI_OUT_EMBED_REC   : 0) | INCHI_OUT_NO_AUX_INFO;            /* minimal compressed output */
                ip->bCtPredecessors  = 1;
                ip->bAbcNumbers      = 1;
                bXml                |= INCHI_OUT_PLAIN_TEXT | INCHI_OUT_PLAIN_TEXT_COMMENTS;
                bXml                &= ~(INCHI_OUT_XML | INCHI_OUT_XML_TEXT_COMMENTS);
            } 
            else if ( !stricmp( pArg, "COMPRESS" ) ) 
            {
                ip->bAbcNumbers = 1;
                ip->bCtPredecessors = 1;             /* compressed output */
            } 

#if ( READ_INCHI_STRING == 1 )
            else if ( !stricmp( pArg, "InChI2InChI" )  ) 
            {
                 /* Read InChI Identifiers and output InChI Identifiers */
                ip->nInputType = INPUT_INCHI;
                ip->bReadInChIOptions |= READ_INCHI_OUTPUT_INCHI;
                ip->bReadInChIOptions &= ~READ_INCHI_TO_STRUCTURE;
            } 

           else if ( !stricmp( pArg, "SplitInChI" )  ) 
            {
                 /* Split InChI Identifiers into components */
                ip->bReadInChIOptions |= READ_INCHI_SPLIT_OUTPUT;
            } 
#endif
            else if ( !stricmp( pArg, "MOLFILENUMBER" ) ) 
            {
                ip->bGetMolfileNumber |= 1;
            } 
            else if ( !stricmp( pArg, "OutputPLAIN" ) ) 
            {
                bXml |=  INCHI_OUT_PLAIN_TEXT;
                bXml &= ~INCHI_OUT_XML;
            }             
            else if ( !stricmp( pArg, "OutputXML" ) ) 
            {
                bXml |= INCHI_OUT_XML;
                bXml &= ~INCHI_OUT_PLAIN_TEXT;
            }             
            else if ( !stricmp( pArg, "OutputANNPLAIN" ) ) 
            {
                bXml |=   INCHI_OUT_PLAIN_TEXT_COMMENTS;
                bXml &=  ~INCHI_OUT_XML_TEXT_COMMENTS;
                bXml |=   INCHI_OUT_WINCHI_WINDOW; /* debug */
            }             
            else if ( !stricmp( pArg, "OutputANNXML" ) ) 
            {
                bXml |=   INCHI_OUT_XML_TEXT_COMMENTS;
                bXml &=  ~INCHI_OUT_PLAIN_TEXT_COMMENTS;
            } else
            if ( !stricmp( pArg, "ONLYEXACT" ) || !stricmp( pArg, "ONLYFIXEDH" ) ) {
                bVer1DefaultMode |=  REQ_MODE_BASIC;
                bVer1DefaultMode &= ~REQ_MODE_TAUT;
            } else
            if ( !stricmp( pArg, "ONLYNONISO" ) ) {
                bVer1DefaultMode |=  REQ_MODE_NON_ISO;
                bVer1DefaultMode &= ~REQ_MODE_ISO;
            } else
            if ( !stricmp( pArg, "TAUT" ) ) {
                bVer1DefaultMode &= ~REQ_MODE_BASIC;
                bVer1DefaultMode |= REQ_MODE_TAUT;
            } else
            if ( !stricmp( pArg, "ONLYRECMET" ) ) {  /* do not disconnect metals */
                bDisconnectCoord = 0;
            } else
            if ( !stricmp( pArg, "ONLYRECSALT" ) ) {  /* do not disconnect salts */
                bDisconnectSalts = 0;
            } else
            if ( !memicmp( pArg, "MOVEPOS:", 8 ) ) {   /* added -- 2010-03-01 DT */
                bMovePositiveCharges = (0 != strtol(pArg+8, NULL, 10));
            } else 
            if ( !memicmp( pArg, "RSB:", 4 )) {
                mdbr = (int)strtol(pArg+4, NULL, 10);
            } else
            if ( !stricmp( pArg, "EQU" ) ) {
                bCompareComponents = CMP_COMPONENTS;
            } else
            if ( !stricmp( pArg, "EQUNONISO" ) ) 
            {
                bCompareComponents = CMP_COMPONENTS | CMP_COMPONENTS_NONISO;
            } 
            else if ( !memicmp( pArg, "OP:", 3 ) ) 
            {
                bOutputPath = 1;
                strncpy(szOutputPath, pArg+3, sizeof(szOutputPath)-1);
            }                         
#endif /* BUILD_WITH_ENG_OPTIONS */            



            /* Display unrecognized option */
            else
            {
                bRecognizedOption = 0;
#ifndef TARGET_LIB_FOR_WINCHI
                inchi_ios_eprint(log_file, "Unrecognized option: \"%s\".\n", pArg);
#endif
            }    
            bVer1Options |= bRecognizedOption;
        

        }  
        /*=== end of parsing stand-alone/library InChI options ===*/



        else if ( ip->num_paths < MAX_NUM_PATHS ) 
        {
            char *sz;
            if ( argv[i] && argv[i][0] ) 
            {
                if ( sz = (char*) inchi_malloc( (strlen(argv[i]) + 1)*sizeof(sz[0])) ) 
                {
                    strcpy( sz, argv[i] );
                }
#ifdef CML_DEBUG
                printf ("1 path %d argv %d %s\n", ip -> num_paths, i, argv [i]);
#endif
                   ip->path[ip->num_paths++] = sz;
            }
        }
    


    } /*  for ( i = 1; i < argc; i ++ )  */



    if ( bHashKey != 0 )
    /* 	Suppress InChIKey calculation if: 
            compressed output OR Inchi2Struct OR Inchi2Inchi */
    {

        if ( (ip->bAbcNumbers ==1) && (ip->bCtPredecessors == 1) )
        {
#ifndef TARGET_LIB_FOR_WINCHI
            inchi_ios_eprint(log_file, "Terminating: generation of InChIKey is not available with 'Compress' option\n");
            return -1;
#endif
            bHashKey = 0;
        }
        if ( ip->nInputType == INPUT_INCHI )
        {
#ifndef TARGET_LIB_FOR_WINCHI
            inchi_ios_eprint(log_file, "Terminating: generation of InChIKey is not available in InChI conversion mode\n");
            return -1;
#endif
            bHashKey = 0;
        }
        else 
        if ( bOutputMolfileOnly == 1 )
        {
#ifndef TARGET_LIB_FOR_WINCHI
            inchi_ios_eprint(log_file, "Terminating: generation of InChIKey is not available with 'OutputSDF' option\n");
            return -1;
#endif
            bHashKey = 0;
        }

    }


    if ( bNameSuffix || bOutputPath ) 
    {
        const char *p = NULL;
        char       *r = NULL;
        char       *sz;
        int  len;
        const char szNUL[] = "NUL"; /* fix for AMD processor: use const char[] instead of just "NUL" constant 2008-11-5 DT */

        /*  find the 1st path */
        for ( i = 0; i < MAX_NUM_PATHS; i ++ ) 
        {
            if ( !p && ip->path[i] && ip->path[i][0] ) 
            {
                p = ip->path[i];
                break;
            }
        }
        /* fix output path */
        if ( bOutputPath && szOutputPath[0] && p ) 
        {
            /* remove last slash */
            len = strlen(szOutputPath);
            if ( len > 0 && szOutputPath[len-1] != INCHI_PATH_DELIM ) 
            {
                szOutputPath[len++] = INCHI_PATH_DELIM;
                szOutputPath[len]   = '\0';
            }
            if ( len > 0 && (r = (char *)strrchr( p, INCHI_PATH_DELIM ) ) && r[1] ) 
            {
                strcat( szOutputPath, r+1 );
                p = szOutputPath;
            }
        }        /*  add missing paths */
        for ( i = 0; p && i < MAX_NUM_PATHS; i ++ ) 
        { 
            /* fix for AMD processor: changed order 2008-11-5 DT */
            if ( !ip->path[i] || !ip->path[i][0] ) 
            {
                len = strlen( p ) + strlen(szNameSuffix) + strlen( ext[i] );
                if ( sz = (char*) inchi_malloc( (len+1)*sizeof(sz[0]) ) ) 
                {
                    strcpy( sz, p );
                    strcat( sz, szNameSuffix );
                    strcat( sz, ext[i] );
                    ip->num_paths++;
                    ip->path[i] =sz;
                }
             } else
            if ( !stricmp( ip->path[i], szNUL ) ) 
            {
                inchi_free( (char *)ip->path[i] ); /* cast deliberately const qualifier */
                ip->path[i] = NULL;
            }
        }
    }




#if ( READ_INCHI_STRING == 1 )
    if ( INPUT_INCHI == ip->nInputType ) 
    {
        bCompareComponents                 = 0;
        /*bDisplayCompositeResults           = 0;*/
#if ( I2S_MODIFY_OUTPUT == 1 )
        if ( !(ip->bReadInChIOptions & READ_INCHI_TO_STRUCTURE ) )
#endif
        {
        bOutputMolfileOnly                 = 0;
        /*bNoStructLabels                    = 1;*/
        bINChIOutputOptions  |= INCHI_OUT_NO_AUX_INFO;
        bINChIOutputOptions  &= ~INCHI_OUT_SHORT_AUX_INFO;
        bINChIOutputOptions  &= ~INCHI_OUT_ONLY_AUX_INFO;
        }
        ip->bDisplayIfRestoreWarnings = bDisplayIfRestoreWarnings;
        if ( !(bINChIOutputOptions &
            
             (INCHI_OUT_SDFILE_ONLY          |    /* not in bINChIOutputOptions yet */
              INCHI_OUT_XML                  |    /* not in bINChIOutputOptions yet */
              INCHI_OUT_PLAIN_TEXT           |    /* not in bINChIOutputOptions yet */
              INCHI_OUT_PLAIN_TEXT_COMMENTS  |    /* not in bINChIOutputOptions yet */
              INCHI_OUT_XML_TEXT_COMMENTS         /* not in bINChIOutputOptions yet */
                                             ) )
#if ( I2S_MODIFY_OUTPUT == 1 )
              && !bOutputMolfileOnly
              && !(bXml & (INCHI_OUT_XML | INCHI_OUT_XML_TEXT_COMMENTS | INCHI_OUT_PLAIN_TEXT | INCHI_OUT_XML_TEXT_COMMENTS))
#endif
           ) {
            bINChIOutputOptions |= INCHI_OUT_PLAIN_TEXT;
        }
    }
#endif


    if ( bVer1Options ) 
    {
        nMode |= bVer1DefaultMode;
    } 
    else if ( bReleaseVersion ) 
    {
        nMode |= nReleaseMode;
    }

#if ( defined(COMPILE_ANSI_ONLY) || defined(TARGET_LIB_FOR_WINCHI) )
    if ( bCompareComponents && !(bDisplay & 1) ) {
        bCompareComponents = 0;
    }
#endif
    /*  Save original options */
    /* nOrigMode = nMode; */
#ifndef COMPILE_ANSI_ONLY
    ip->dp.sdp.nFontSize         = nFontSize;
    ip->dp.sdp.ulDisplTime       = *ulDisplTime;
    ip->bDisplay                 = bDisplay;
#ifdef TARGET_LIB_FOR_WINCHI
    ip->bDisplayCompositeResults = bDisplay;
#else
    ip->bDisplayCompositeResults = bDisplayCompositeResults;
#endif
#else
    ip->bDisplayEachComponentINChI = 0;
    bCompareComponents            = 0;
#endif
    ip->bMergeAllInputStructures = bMergeAllInputStructures;
    ip->bDoNotAddH               = bDoNotAddH;
    /*  set default options */
    if ( !nMode || nMode == REQ_MODE_STEREO ) {
        /*  requested all output */
        nMode |= (REQ_MODE_BASIC | REQ_MODE_TAUT | REQ_MODE_ISO | REQ_MODE_NON_ISO | REQ_MODE_STEREO);
    } else {
        if ( !(nMode & (REQ_MODE_BASIC | REQ_MODE_TAUT)) ) {
            nMode |= (REQ_MODE_BASIC | REQ_MODE_TAUT);
        }
        if ( (nMode & REQ_MODE_STEREO) && !(nMode & (REQ_MODE_ISO | REQ_MODE_NON_ISO)) ) {
            nMode |= (REQ_MODE_ISO | REQ_MODE_NON_ISO);
        }
    }
    /*  if the user requested isotopic then unconditionally add non-isotopic output. */
    if ( nMode & REQ_MODE_ISO ) {
        nMode |= REQ_MODE_NON_ISO;
    }
#if ( MIN_SB_RING_SIZE > 0 )
    if ( mdbr ) {
        nMinDbRinSize = mdbr;
    }
    nMode |= (nMinDbRinSize << REQ_MODE_MIN_SB_RING_SHFT) & REQ_MODE_MIN_SB_RING_MASK;
#endif
    /*  input file */
    if ( ip->nInputType == INPUT_NONE && ip->num_paths > 0 ) {
        ip->nInputType = INPUT_MOLFILE; /*  default */
#if ( ADD_CMLPP == 1 )
        {
            const char *p;
            if ( ip->path[0] && ( p = strrchr(ip->path[0], '.' ) ) &&
                 !stricmp( p, ".cml") ) {
                ip->nInputType = INPUT_CMLFILE;
            }
        }
#endif
    }
    ip->nMode = nMode;
    if ( (bCompareComponents & CMP_COMPONENTS) && (nMode & REQ_MODE_BASIC) ) {
        bCompareComponents |= CMP_COMPONENTS_NONTAUT; /* compare non-tautomeric */
    }
    ip->bCompareComponents = bCompareComponents;

    ip->bINChIOutputOptions = bINChIOutputOptions | (bOutputMolfileOnly? INCHI_OUT_SDFILE_ONLY : 0);
    if ( bOutputMolfileOnly ) {
        bXml &= ~(INCHI_OUT_XML                 | INCHI_OUT_PLAIN_TEXT |
                  INCHI_OUT_PLAIN_TEXT_COMMENTS | INCHI_OUT_XML_TEXT_COMMENTS | INCHI_OUT_TABBED_OUTPUT);
#if ( SDF_OUTPUT_DT == 1 )
        ip->bINChIOutputOptions |= bOutputMolfileDT?    INCHI_OUT_SDFILE_ATOMS_DT : 0;
        ip->bINChIOutputOptions |= bOutputMolfileSplit? INCHI_OUT_SDFILE_SPLIT : 0;
#endif
    }
    if ( bXml & INCHI_OUT_XML ) {
        bXml &= ~(INCHI_OUT_PLAIN_TEXT | INCHI_OUT_XML_TEXT_COMMENTS | INCHI_OUT_TABBED_OUTPUT);
    }
#ifdef TARGET_LIB_FOR_WINCHI
    if ( !(bDisplay & 1) ) {
        bXml &= ~(INCHI_OUT_PLAIN_TEXT_COMMENTS | INCHI_OUT_XML_TEXT_COMMENTS); /* do not output comments in wINChI text file results */
    } else {
        bXml |= INCHI_OUT_WINCHI_WINDOW;
    }
#endif
    ip->bINChIOutputOptions |= bXml;
    ip->bNoStructLabels     = bNoStructLabels;

    if ( bForcedChiralFlag ) {
        ip->bChiralFlag = bForcedChiralFlag;
    }

    /*******************************************/
    /*       tautomeric/salts settings         */
    /*******************************************/

    ip->bTautFlags     = 0;   /* initialize */
    ip->bTautFlagsDone = 0;   /* initialize */

    /* find regular tautomerism */
    ip->bTautFlags |= TG_FLAG_TEST_TAUT__ATOMS;
    /* disconnect salts */
    ip->bTautFlags |= bDisconnectSalts?         TG_FLAG_DISCONNECT_SALTS    : 0;
    /* if possible find long-range H/(-) taut. on =C-OH, >C=O    */
    ip->bTautFlags |= bAcidTautomerism?         TG_FLAG_TEST_TAUT__SALTS    : 0;
    /* allow long-range movement of N(+), P(+) charges           */
    ip->bTautFlags |= bMovePositiveCharges?     TG_FLAG_MOVE_POS_CHARGES    : 0;
    /* multi-attachement long-range H/(-) taut. on =C-OH, >C=O   */
    ip->bTautFlags |= (bAcidTautomerism > 1)?   TG_FLAG_TEST_TAUT2_SALTS    : 0;
    /* (debug) allow to find long-range H-only tautomerism on =C-OH, >C=O */
    ip->bTautFlags |= (bUnchargedAcidTaut==1)?  TG_FLAG_ALLOW_NO_NEGTV_O    : 0;
    /* merge =C-OH and >C=O containing t-groups and other =C-OH groups */
    ip->bTautFlags |= bMergeSaltTGroups?        TG_FLAG_MERGE_TAUT_SALTS    : 0;
    ip->bTautFlags |= bDisconnectCoord?         TG_FLAG_DISCONNECT_COORD    : 0;
    ip->bTautFlags |=(bDisconnectCoord &&
                      bReconnectCoord)?         TG_FLAG_RECONNECT_COORD     : 0;
    ip->bTautFlags |= bDisconnectCoordChkVal?   TG_FLAG_CHECK_VALENCE_COORD : 0;
    ip->bTautFlags |= bTgFlagVariableProtons?   TG_FLAG_VARIABLE_PROTONS     : 0;
    ip->bTautFlags |= bTgFlagHardAddRenProtons? TG_FLAG_HARD_ADD_REM_PROTONS : 0;
    ip->bTautFlags |= bKetoEnolTaut?            TG_FLAG_KETO_ENOL_TAUT : 0;
    ip->bTautFlags |= b15TautNonRing?           TG_FLAG_1_5_TAUT : 0;

#ifdef STEREO_WEDGE_ONLY
    ip->bTautFlags  |= bPointedEdgeStereo?      TG_FLAG_POINTED_EDGE_STEREO  : 0;
#endif
#if ( FIX_ADJ_RAD == 1 )
    ip->bTautFlags  |= bFixAdjacentRad?         TG_FLAG_FIX_ADJ_RADICALS : 0;
#endif
    ip->bTautFlags  |= bAddPhosphineStereo?     TG_FLAG_PHOSPHINE_STEREO : 0;
    ip->bTautFlags  |= bAddArsineStereo?        TG_FLAG_ARSINE_STEREO : 0;
    ip->bTautFlags  |= bFixSp3bug?              TG_FLAG_FIX_SP3_BUG   : 0;


    if (bFixFB2)
    {
#if ( FIX_ISO_FIXEDH_BUG == 1 )
        ip->bTautFlags  |= TG_FLAG_FIX_ISO_FIXEDH_BUG; /* accomodate FIX_ISO_FIXEDH_BUG */
#endif
#if ( FIX_TERM_H_CHRG_BUG == 1 )
        ip->bTautFlags  |= TG_FLAG_FIX_TERM_H_CHRG_BUG; /* accomodate FIX_TERM_H_CHRG_BUG */
#endif
#if ( FIX_TRANSPOSITION_CHARGE_BUG == 1 )
        ip->bINChIOutputOptions |= INCHI_OUT_FIX_TRANSPOSITION_CHARGE_BUG;
#endif
    }




    if ( !ip->nInputType ) 
        ip->nInputType = INPUT_MOLFILE;

    /* Check if /SNon requested turn OFF SUU/SLUUD */
    if ( ! (ip->nMode & REQ_MODE_STEREO) )
    {
        ip->nMode &= ~REQ_MODE_DIFF_UU_STEREO;
        ip->nMode &= ~(REQ_MODE_SB_IGN_ALL_UU | REQ_MODE_SC_IGN_ALL_UU);
    }


    /* Standard InChI ? */
    if (bStdFormat)  
    {
        ip->bINChIOutputOptions |= INCHI_OUT_STDINCHI;
    }

    /* InChIKey ? */
    if ( !bHashKey )
        ip->bCalcInChIHash = INCHIHASH_NONE;
    else
        ip->bCalcInChIHash = INCHIHASH_KEY;        
    
    /* Extension(s) to hash (in non-std mode only) ? */
    if ( !bHashKey )
    {
        if ( (bHashXtra1!=0) || (bHashXtra2!=0) )
            inchi_ios_eprint(log_file,
                            "Hash extension(s) not generated: InChIKey not requested");
    }
    else
    {
        if ( bHashXtra1 )
        {
            if ( bHashXtra2 )   ip->bCalcInChIHash = INCHIHASH_KEY_XTRA1_XTRA2;
            else                ip->bCalcInChIHash = INCHIHASH_KEY_XTRA1;
        }
        else if ( bHashXtra2 )
        {
            ip->bCalcInChIHash = INCHIHASH_KEY_XTRA2;
        }
    }


    return 0;
}



/*******************************************************************/
int PrintInputParms( INCHI_IOSTREAM *log_file, INPUT_PARMS *ip )
{
INCHI_MODE nMode = ip->nMode;
int k;
int bStdFormat = 1;
int first=1;
           
     
#ifdef TARGET_LIB_FOR_WINCHI
    int bInChI2Struct = 0; /* as of December 2008, winchi-1 does not convert InChI to structure */
#else
    int bInChI2Struct = (ip->bReadInChIOptions & READ_INCHI_TO_STRUCTURE) && ip->nInputType == INPUT_INCHI;
#endif
            
    if ( !(ip->bINChIOutputOptions & INCHI_OUT_STDINCHI) )
        bStdFormat = 0;


    /* some stereo */
    if ( ! (nMode & REQ_MODE_STEREO) ) 
    {
        inchi_ios_eprint( log_file, "Using specific structure perception features:\n");
        first = 0;
        inchi_ios_eprint( log_file, "  Stereo OFF\n");
    }
    else
    {
        if ( ! (TG_FLAG_POINTED_EDGE_STEREO & ip->bTautFlags) ) 
        {
            if (first)
            {
                inchi_ios_eprint( log_file, "Using specific structure perception features:\n");
                first = 0;
            }
            inchi_ios_eprint( log_file, "  Both ends of wedge point to stereocenters\n");
        }
    }
    if ( ip->bDoNotAddH ) 
    {
        if (first)
            inchi_ios_eprint( log_file, "Using specific structure perception features:\n");
        inchi_ios_eprint( log_file, "  Do not add H\n");
    }


    /*  Generation/conversion indicator */
    if (bStdFormat)
    {
        if ( !(ip->bINChIOutputOptions & INCHI_OUT_SDFILE_ONLY) && !bInChI2Struct )
            inchi_ios_eprint( log_file, "Generating standard InChI\n" );

#if ( !defined(TARGET_API_LIB) && !defined(TARGET_LIB_FOR_WINCHI) && !defined(TARGET_EXE_USING_API) )
        /* effective only in command line program InChI or stdInChI */
        else if ( bInChI2Struct )             
            inchi_ios_eprint( log_file, "Converting InChI(s) to structure(s) in %s\n",
                              (ip->bINChIOutputOptions & INCHI_OUT_SDFILE_ONLY)? 
                              "MOL format" : "aux. info format" );
#endif
    }
    else
        inchi_ios_eprint( log_file, "Generating non-standard InChI with the options: \n" );

    
    /* SDfile output */

    if ( ip->bINChIOutputOptions & INCHI_OUT_SDFILE_ONLY ) 
    {
        inchi_ios_eprint( log_file, 
                         "Output SDfile only without stereochemical information and atom coordinates%s\n",
                         (ip->bINChIOutputOptions & INCHI_OUT_SDFILE_ATOMS_DT)? 
                         "\n(write H isotopes as D, T)":"" );
    }
    
    /* Fixed/Mobile H */
    if (!bStdFormat)
    {
        if( (nMode & (REQ_MODE_BASIC | REQ_MODE_TAUT )) == (REQ_MODE_BASIC | REQ_MODE_TAUT) ) 
            inchi_ios_eprint( log_file, "  Mobile H Perception OFF (include FixedH layer)\n" );
        else if( (nMode & (REQ_MODE_BASIC | REQ_MODE_TAUT )) == (REQ_MODE_TAUT) ) 
            inchi_ios_eprint( log_file, "  Mobile H Perception ON  (omit FixedH layer)\n" );
        else if( (nMode & (REQ_MODE_BASIC | REQ_MODE_TAUT )) == (REQ_MODE_BASIC) ) 
            inchi_ios_eprint( log_file, "  Mobile H ignored\n" );
        else 
            inchi_ios_eprint( log_file, "  Undefined Mobile H mode\n" );

        if ( (ip->bTautFlags & TG_FLAG_VARIABLE_PROTONS) ) 
         if ( !(ip->bTautFlags & TG_FLAG_HARD_ADD_REM_PROTONS) ) 
            inchi_ios_eprint( log_file, "  Disabled Aggressive (De)protonation\n" );

#if ( FIND_RING_SYSTEMS != 1 )
        inchi_ios_eprint( log_file, "  %s5-, 6-, 7-memb. ring taut. ignored\n", i?"; ":"");
#endif

        /* RecMet */
        if ( ip->bTautFlags & TG_FLAG_DISCONNECT_COORD ) 
        {
            if ( ip->bTautFlags & TG_FLAG_RECONNECT_COORD ) 
                inchi_ios_eprint( log_file, "  Include bonds to metals\n");
            else
                inchi_ios_eprint( log_file, "  Do not reconnect metals (omit RecMet layer)\n");
        }
        else
            inchi_ios_eprint( log_file, "  Do not disconnect metals\n");

        /* isotopic - always ON, output disabled. 09-17-2009*/
        /*
        if ( nMode & REQ_MODE_ISO ) 
            inchi_ios_eprint( log_file, "  Isotopic ON\n");
        else if ( nMode & REQ_MODE_NON_ISO ) 
            inchi_ios_eprint( log_file, "  Isotopic OFF\n");
        */

#if ( FIX_ADJ_RAD == 1 )
    if ( ip->bTautFlags & TG_FLAG_FIX_ADJ_RADICALS ) 
        inchi_ios_eprint( log_file, "Fix Adjacent Radicals\n" );
#endif

        /*  stereo */
        if ( nMode & REQ_MODE_STEREO ) 
        {
            inchi_ios_eprint( log_file,  "  %s%s%s%sStereo ON\n",
                     ( nMode & REQ_MODE_NOEQ_STEREO )?     "Slow ":"",
                     ( nMode & REQ_MODE_REDNDNT_STEREO )?  "Redund. ":"",
                     ( nMode & REQ_MODE_NO_ALT_SBONDS)?    "No AltBond ":"",

                     ( nMode & REQ_MODE_RACEMIC_STEREO)?   "Racemic " :
                     ( nMode &  REQ_MODE_RELATIVE_STEREO)? "Relative " :
                     ( nMode &  REQ_MODE_CHIR_FLG_STEREO)? "Chiral Flag " : "Absolute " );
            if ( 0 == (nMode & (REQ_MODE_SB_IGN_ALL_UU | REQ_MODE_SC_IGN_ALL_UU)) ) 
                inchi_ios_eprint( log_file, "  Include undefined/unknown stereogenic centers and bonds\n");
            else if ( REQ_MODE_SC_IGN_ALL_UU == (nMode & (REQ_MODE_SB_IGN_ALL_UU | REQ_MODE_SC_IGN_ALL_UU)) ) 
                inchi_ios_eprint( log_file, "  Omit undefined/unknown stereogenic centers\n");
            else if ( REQ_MODE_SB_IGN_ALL_UU == (nMode & (REQ_MODE_SB_IGN_ALL_UU | REQ_MODE_SC_IGN_ALL_UU)) ) 
                inchi_ios_eprint( log_file, "  Omit undefined/unknown stereogenic bonds\n");
            else
                /*case REQ_MODE_SB_IGN_ALL_UU | REQ_MODE_SC_IGN_ALL_UU*/
                inchi_ios_eprint( log_file, "  Omit undefined/unknown stereogenic centers and bonds\n");

            if ( 0 != (nMode & REQ_MODE_DIFF_UU_STEREO) )
            {
                inchi_ios_eprint( log_file, "  Make labels for unknown and undefined stereo different\n");
            }


#if ( defined(MIN_SB_RING_SIZE) && MIN_SB_RING_SIZE > 0 )
            k = (ip->nMode & REQ_MODE_MIN_SB_RING_MASK) >> REQ_MODE_MIN_SB_RING_SHFT;
            if ( bRELEASE_VERSION != 1 || k != MIN_SB_RING_SIZE ) 
            {
                if ( k >= 3 ) 
                    inchi_ios_eprint( log_file, "  Min. stereobond ring size: %d\n", k );
                else 
                    inchi_ios_eprint( log_file, "  Min. stereobond ring size: NONE\n" );
            }
#endif
        } /* stereo */

    } /* !bStdFormat */

  

    if ( !bStdFormat )
    {
        if ( TG_FLAG_KETO_ENOL_TAUT & ip->bTautFlags) 
            inchi_ios_eprint( log_file, "  Account for keto-enol tautomerism\n");
        else
            inchi_ios_eprint( log_file, "  Do not account for keto-enol tautomerism\n");    
        if ( TG_FLAG_1_5_TAUT & ip->bTautFlags) 
            inchi_ios_eprint( log_file, "  Account for 1,5-tautomerism\n");
        else
            inchi_ios_eprint( log_file, "  Do not account for 1,5-tautomerism\n");
#ifdef BUILD_WITH_ENG_OPTIONS
        if ( TG_FLAG_PHOSPHINE_STEREO & ip->bTautFlags ) 
            inchi_ios_eprint( log_file, "  Include phosphine stereochemistry\n");
        else
            inchi_ios_eprint( log_file, "  Do not include phosphine stereochemistry\n");
        if ( TG_FLAG_ARSINE_STEREO & ip->bTautFlags ) 
            inchi_ios_eprint( log_file, "  Include arsine stereochemistry\n");
        else
            inchi_ios_eprint( log_file, "  Do not include arsine stereochemistry\n");
        if ( ! (TG_FLAG_FIX_SP3_BUG & ip->bTautFlags) ) 
                inchi_ios_eprint( log_file, "  Turned OFF fix of bug leading to missing or undefined sp3 parity\n");
        if ( !(TG_FLAG_FIX_ISO_FIXEDH_BUG & ip->bTautFlags) )
            inchi_ios_eprint( log_file, "  Turned OFF bug-fixes found after v.1.02b release\n");
        if ( !(ip->bFixNonUniformDraw) ) 
            inchi_ios_eprint( log_file, "  Turned OFF fixes of non-uniform drawing issues\n");
        if ( ! (TG_FLAG_MOVE_POS_CHARGES & ip->bTautFlags) ) 
                inchi_ios_eprint( log_file, "  MovePos turned OFF\n");
#endif
    } /* !bStdFormat */

    if ( ip->bCalcInChIHash != INCHIHASH_NONE ) 
    {
        if (bStdFormat)
            inchi_ios_eprint( log_file, "Generating standard InChIKey\n");
        else
            inchi_ios_eprint( log_file, "Generating InChIKey\n");
        if ( ip->bCalcInChIHash == INCHIHASH_KEY_XTRA1 )
            inchi_ios_eprint( log_file, "Generating hash extension (1st block)\n");
        else if ( ip->bCalcInChIHash == INCHIHASH_KEY_XTRA2 )
            inchi_ios_eprint( log_file, "Generating hash extension (2nd block)\n");
        else if ( ip->bCalcInChIHash == INCHIHASH_KEY_XTRA1_XTRA2 )
            inchi_ios_eprint( log_file, "Generating hash extension (two blocks)\n");
    }

    if ( ip->bINChIOutputOptions & INCHI_OUT_SAVEOPT )
    {
        inchi_ios_eprint( log_file, "Saving InChI creation options" );
        if ( bStdFormat )
        {
            inchi_ios_eprint( log_file, " suppressed for standard InChI" );
            /* NB: actual suppression takes place on InChI serialization */
            /* (as on e.g. Inchi2Inchi conversion it may appear that we create non-std */
            /*  InChI instead of standard one) */
        }
        inchi_ios_eprint( log_file, "\n" );

    }


    if ( ip->bAllowEmptyStructure ) 
        inchi_ios_eprint( log_file, "Issue warning on empty structure\n" );




    /* Input */
    if ( ip->nInputType ) 
    {
        inchi_ios_eprint( log_file, "Input format: %s",
            ip->nInputType == INPUT_MOLFILE?     "MOLfile"       :
            ip->nInputType == INPUT_SDFILE?      "SDfile"        :
            ip->nInputType == INPUT_CMLFILE?     "CMLfile"       :
#if ( READ_INCHI_STRING == 1 )
            ip->nInputType == INPUT_INCHI?       "InChI (plain identifier)" :
#endif
            ip->nInputType == INPUT_INCHI_XML?   "InChI AuxInfo (xml)"   :
            ip->nInputType == INPUT_INCHI_PLAIN? "InChI AuxInfo (plain)" : "Unknown" );
        if ( (ip->nInputType == INPUT_MOLFILE || ip->nInputType == INPUT_SDFILE) &&
             ip->bGetMolfileNumber ) 
            inchi_ios_eprint( log_file, "  (attempting to read Molfile number)" );
        inchi_ios_eprint( log_file, "\n");
    }

    if ( ip->szSdfDataHeader[0] && ip->nInputType != INPUT_SDFILE ) 
        inchi_ios_eprint( log_file, "  SDfile data header: \"%s\"\n", ip->szSdfDataHeader);

    /* Output */

    inchi_ios_eprint( log_file, "Output format: %s%s\n", 
        (ip->bINChIOutputOptions & INCHI_OUT_PLAIN_TEXT)?  "Plain text" :
        (ip->bINChIOutputOptions & INCHI_OUT_XML)?         "XML":
        ((ip->bINChIOutputOptions & INCHI_OUT_SDFILE_ONLY) && bInChI2Struct )? "SDfile only (without stereochemical info and atom coordinates)" :
        ((ip->bINChIOutputOptions & INCHI_OUT_SDFILE_ONLY) && !bInChI2Struct)? "SDfile only" : "Unknown",

        ((ip->bINChIOutputOptions & INCHI_OUT_PLAIN_TEXT) &&
        (ip->bINChIOutputOptions & INCHI_OUT_TABBED_OUTPUT))? ", tabbed":"");

#if ( bRELEASE_VERSION == 1 )
    if ( ip->bCtPredecessors || ip->bAbcNumbers ) 
    {
        if ( ip->bCtPredecessors && ip->bAbcNumbers ) 
            inchi_ios_eprint( log_file, "Representation: Compressed\n");
        else 
            inchi_ios_eprint( log_file, "Connection table: %s, %s\n",
                ip->bCtPredecessors? "Predecessor_numbers(closures)":"Canon_numbers(branching, ring closures)",
                ip->bAbcNumbers?     "Shorter alternative":"Numerical");
    }
#else
    if ( (bRELEASE_VERSION != 1) || ip->bCtPredecessors || ip->bAbcNumbers ) 
        inchi_ios_eprint( log_file, "Connection table: %s, %s\n",
            ip->bCtPredecessors? "Predecessor_numbers(closures)":"Canon_numbers(branching, ring closures)",
            ip->bAbcNumbers?     "Shorter alternative":"Numerical");
    else 
        inchi_ios_eprint( log_file, "Representation: Numerical");
#endif

    if ( !(ip->bINChIOutputOptions & INCHI_OUT_SDFILE_ONLY) ) 
    {
        if( ip->bINChIOutputOptions & INCHI_OUT_NO_AUX_INFO ) 
            inchi_ios_eprint( log_file, "Aux. info suppressed\n");
        else if ( ip->bINChIOutputOptions & INCHI_OUT_SHORT_AUX_INFO ) 
            inchi_ios_eprint( log_file, "Minimal Aux. info\n");
        else 
            inchi_ios_eprint( log_file, "Full Aux. info\n");
    }



    if ( ip->msec_MaxTime ) 
    {
        unsigned long seconds = ip->msec_MaxTime/1000;
        unsigned long milliseconds = (ip->msec_MaxTime%1000);
        inchi_ios_eprint( log_file, "Timeout per structure: %lu.%03lu sec\n", seconds, milliseconds);
    } 
    else 
        inchi_ios_eprint( log_file, "No timeout");
    inchi_ios_eprint( log_file, "Up to %d atoms per structure\n", MAX_ATOMS);
    if ( ip->first_struct_number > 1 ) 
        inchi_ios_eprint( log_file, "Skipping %ld structure%s\n", ip->first_struct_number-1, ip->first_struct_number==2? "":"s" );
    if ( ip->last_struct_number > 0 ) 
        inchi_ios_eprint( log_file, "Terminate after structure #%ld\n", ip->last_struct_number );
    if ( ip->bSaveWarningStructsAsProblem && ip->path[3] && ip->path[3][0] ) 
        inchi_ios_eprint( log_file, "Saving warning structures into the problem file\n");
    if ( ip->bSaveAllGoodStructsAsProblem && ip->path[3] && ip->path[3][0] ) 
        inchi_ios_eprint( log_file, "Saving only all good structures into the problem file\n");
    /*  Report debug modes */
#if ( bRELEASE_VERSION != 1 )
    inchi_ios_eprint( log_file, "Release version = NO\n");
#endif


#if ( TRACE_MEMORY_LEAKS == 1 && defined(_DEBUG) )
    inchi_ios_eprint( log_file, "Tracing memory leaks (SLOW)\n");
#endif

    inchi_ios_eprint( log_file, "\n" );



#if ( bRELEASE_VERSION != 1 )
#if ( FIND_RING_SYSTEMS == 1 )
    inchi_ios_eprint( log_file, "Find ring systems=Y\nTautomers:\n" );
    inchi_ios_eprint( log_file, " 4-pyridinol=%s\n", TAUT_4PYRIDINOL_RINGS==1? "Y":"N");
    inchi_ios_eprint( log_file, " pyrazole=%s\n", TAUT_PYRAZOLE_RINGS==1? "Y":"N");
    inchi_ios_eprint( log_file, " tropolone=%s\n", TAUT_TROPOLONE_7==1? "Y":"N");
    inchi_ios_eprint( log_file, " tropolone-5=%s\n", TAUT_TROPOLONE_5==1? "Y":"N");
    inchi_ios_eprint( log_file, "Only chain attachments to tautomeric rings=%s\n", TAUT_RINGS_ATTACH_CHAIN==1? "Y":"N");
#endif
    if ( ip->bGetSdfileId ) 
        inchi_ios_eprint( log_file, "Extracting SDfile IDs\n");

    inchi_ios_eprint( log_file, "\nDbg: MOVE_CHARGES=%d\n",
                           0!=(ip->bTautFlags&TG_FLAG_MOVE_POS_CHARGES));
    inchi_ios_eprint( log_file, "     REPLACE_ALT_WITH_TAUT=%d; NEUTRALIZE_ENDPOINTS=%d; BNS_PROTECT_FROM_TAUT=%d\n",
                                  REPLACE_ALT_WITH_TAUT,    NEUTRALIZE_ENDPOINTS, BNS_PROTECT_FROM_TAUT);
    inchi_ios_eprint( log_file, "     DISCONNECT_SALTS=%d;   TEST_TAUT_SALTS=%d;    TEST_TAUT2_SALTS=%d\n",
                                  0!=(ip->bTautFlags&TG_FLAG_DISCONNECT_SALTS),
                                  0!=(ip->bTautFlags&TG_FLAG_TEST_TAUT__SALTS),
                                  0!=(ip->bTautFlags&TG_FLAG_TEST_TAUT2_SALTS));

    inchi_ios_eprint( log_file, "     CHARGED_ACID_TAUT_ONLY=%d MERGE_TAUT_SALTS=%d\n",
                                  0==(ip->bTautFlags&TG_FLAG_ALLOW_NO_NEGTV_O),
                                  0!=(ip->bTautFlags&TG_FLAG_MERGE_TAUT_SALTS));
    inchi_ios_eprint( log_file, "     DISCONNECT_COORD=%d\n", 0!=(ip->bTautFlags&TG_FLAG_DISCONNECT_COORD) );
#if ( TEST_RENUMB_ATOMS == 1 )
    inchi_ios_eprint( log_file, "\nDbg: TEST_RENUMB_ATOMS=%d; TEST_RENUMB_NEIGH=%d; TEST_RENUMB_SWITCH=%d\n",
                                  TEST_RENUMB_ATOMS,    TEST_RENUMB_NEIGH,    TEST_RENUMB_SWITCH );
    inchi_ios_eprint( log_file, "     TEST_RENUMB_ATOMS_SAVE_LONGEST=%d\n",
                                  TEST_RENUMB_ATOMS_SAVE_LONGEST);
#endif

#endif /* ( bRELEASE_VERSION != 1 ) */


    return 0;
}



/************************************************************************************/
void HelpCommandLineParms( INCHI_IOSTREAM *f )
{
    if ( !f )
        return;

#if ( bRELEASE_VERSION == 1 )

    inchi_ios_print_nodisplay( f, 
#ifdef TARGET_EXE_USING_API
        "%s ver %s%s.\n\nUsage:\ninchi_main inputFile [outputFile [logFile [problemFile]]] [%coption[ %coption...]]\n", 
        INCHI_NAME, INCHI_VERSION, TARGET_ID_STRING, 
        INCHI_OPTION_PREFX, INCHI_OPTION_PREFX); 
#else
        "%s ver %s%s.\n\nUsage:\ninchi-%s inputFile [outputFile [logFile [problemFile]]] [%coption [%coption...]]\n", 	     
        INCHI_NAME, INCHI_VERSION, TARGET_ID_STRING, 
        INCHI_VERSION, INCHI_OPTION_PREFX, INCHI_OPTION_PREFX); 
#if ( BUILD_WITH_AMI == 1 )
    inchi_ios_print_nodisplay( f, 
        "inchi-%s inputFiles... %cAMI [%coption[ %coption...]]\n", 	     
        INCHI_VERSION, INCHI_OPTION_PREFX, INCHI_OPTION_PREFX, INCHI_OPTION_PREFX); 
#endif
#endif
         
    inchi_ios_print_nodisplay( f, "\nOptions:\n");
    
    inchi_ios_print_nodisplay( f, "\nInput\n");
    inchi_ios_print_nodisplay( f, "  STDIO       Use standard input/output streams\n");
    inchi_ios_print_nodisplay( f, "  InpAux      Input structures in %s default aux. info format\n              (for use with STDIO)\n", INCHI_NAME);
    inchi_ios_print_nodisplay( f, "  SDF:DataHeader Read from the input SDfile the ID under this DataHeader\n");
#if ( ADD_CMLPP == 1 )
    inchi_ios_print_nodisplay( f, "  CML         Input in CML format (default for input file .CML extension)\n");
#endif
/*
    inchi_ios_print_nodisplay( f, "  START:n     Skip structures up to n-th one\n");
    inchi_ios_print_nodisplay( f, "  END:m       Skip structures after m-th one\n");
*/
#if ( BUILD_WITH_AMI == 1 )
    inchi_ios_print_nodisplay( f, "  AMI         Allow multiple input files (wildcards supported)\n");
#endif

    inchi_ios_print_nodisplay( f, "Output\n");
    inchi_ios_print_nodisplay( f, "  AuxNone     Omit auxiliary information (default: Include)\n");
    inchi_ios_print_nodisplay( f, "  SaveOpt     Save custom InChI creation options (non-standard InChI)\n");
    inchi_ios_print_nodisplay( f, "  NoLabels    Omit structure number, DataHeader and ID from %s output\n", INCHI_NAME);
    inchi_ios_print_nodisplay( f, "  Tabbed      Separate structure number, %s, and AuxInfo with tabs\n", INCHI_NAME);
    
    /*inchi_ios_print_nodisplay( f, "  Compress    Compressed output\n"); */
    /*inchi_ios_print_nodisplay( f, "    FULL        Standard set of options for Full Verbose Output\n");*/
    /*inchi_ios_print_nodisplay( f, "    MIN         Standard set of options for Minimal Concise Output\n");*/
#if ( defined(_WIN32) && defined(_MSC_VER) && !defined(COMPILE_ANSI_ONLY) && !defined(TARGET_API_LIB) )
    inchi_ios_print_nodisplay( f, "  D           Display the structures\n");
    inchi_ios_print_nodisplay( f, "  EQU         Display sets of identical components\n");
    inchi_ios_print_nodisplay( f, "  Fnumber     Set display Font size in number of points\n");
#endif
    /*inchi_ios_print_nodisplay( f, "    PLAIN       Plain text output (Default: XML format)\n");*/
    inchi_ios_print_nodisplay( f, "  OutputSDF   Convert %s created with default aux. info to SDfile\n", INCHI_NAME);
#if ( SDF_OUTPUT_DT == 1 )
    inchi_ios_print_nodisplay( f, "  SdfAtomsDT  Output Hydrogen Isotopes to SDfile as Atoms D and T\n");
#endif
#if ( BUILD_WITH_AMI == 1 )
    inchi_ios_print_nodisplay( f, "  AMIOutStd   Write output to stdout (in AMI mode)\n");
    inchi_ios_print_nodisplay( f, "  AMILogStd   Write log to stderr (in AMI mode)\n");
    inchi_ios_print_nodisplay( f, "  AMIPrbNone  Suppress creation of problem files (in AMI mode)\n");
#endif
    inchi_ios_print_nodisplay( f, "Structure perception\n");
    inchi_ios_print_nodisplay( f, "  SNon        Exclude stereo (default: include absolute stereo)\n");
    inchi_ios_print_nodisplay( f, "  NEWPSOFF    Both ends of wedge point to stereocenters (default: a narrow end)\n");
    inchi_ios_print_nodisplay( f, "  DoNotAddH   All H are explicit (default: add H according to usual valences)\n");
    
#ifndef USE_STDINCHI_API
    inchi_ios_print_nodisplay( f, "Stereo perception modifiers (non-standard InChI)\n");
    inchi_ios_print_nodisplay( f, "  SRel        Relative stereo\n");
    inchi_ios_print_nodisplay( f, "  SRac        Racemic stereo\n");
    inchi_ios_print_nodisplay( f, "  SUCF        Use Chiral Flag: On means Absolute stereo, Off - Relative\n"); 

    inchi_ios_print_nodisplay( f, "Customizing InChI creation (non-standard InChI)\n");
    inchi_ios_print_nodisplay( f, "  SUU         Always include omitted unknown/undefined stereo\n");
    inchi_ios_print_nodisplay( f, "  SLUUD       Make labels for unknown and undefined stereo different\n");
    inchi_ios_print_nodisplay( f, "  RecMet      Include reconnected metals results\n");
    inchi_ios_print_nodisplay( f, "  FixedH      Include Fixed H layer\n");
    inchi_ios_print_nodisplay( f, "  KET         Account for keto-enol tautomerism (experimental)\n");
    inchi_ios_print_nodisplay( f, "  15T         Account for 1,5-tautomerism (experimental)\n");
#endif

    inchi_ios_print_nodisplay( f, "Generation\n");
    inchi_ios_print_nodisplay( f, "  Wnumber     Set time-out per structure in seconds; W0 means unlimited\n");
    inchi_ios_print_nodisplay( f, "  WarnOnEmptyStructure Warn and produce empty %s for empty structure\n", INCHI_NAME);
    inchi_ios_print_nodisplay( f, "  Key         Generate InChIKey\n");
    inchi_ios_print_nodisplay( f, "  XHash1      Generate hash extension (to 256 bits) for 1st block of InChIKey\n");
    inchi_ios_print_nodisplay( f, "  XHash2      Generate hash extension (to 256 bits) for 2nd block of InChIKey\n");


    inchi_ios_print_nodisplay( f, "Conversion\n");
#ifdef TARGET_EXE_USING_API
    inchi_ios_print_nodisplay( f, "  InChI2InChI  Test mode: Mol/SDfile->%s->%s\n", INCHI_NAME, INCHI_NAME);
    inchi_ios_print_nodisplay( f, "  InChI2Struct Test mode: Mol/SDfile->%s->Structure->(%s+AuxInfo)\n", INCHI_NAME, INCHI_NAME);
#else
    inchi_ios_print_nodisplay( f, "  InChI2InChI  Convert %s string(s) into %s string(s)\n", INCHI_NAME, INCHI_NAME); 
    inchi_ios_print_nodisplay( f, "  InChI2Struct Convert InChI string(s) to structure(s) in InChI aux.info format\n");
#endif

#ifdef BUILD_WITH_ENG_OPTIONS
#if 0
    inchi_ios_print_nodisplay( f, "Engineering options (for testing only)\n");
    inchi_ios_print_nodisplay( f, "  NoADP       Disable Aggressive Deprotonation\n");
#if ( FIX_ADJ_RAD == 1 )
    inchi_ios_print_nodisplay( f, "  FixRad      Fix Adjacent Radicals\n");
#endif
    inchi_ios_print_nodisplay( f, "  SPXYZOFF    Do not include Phosphines Stereochemistry\n");
    inchi_ios_print_nodisplay( f, "  SAsXYZOFF   Do not include Arsines Stereochemistry\n");
    inchi_ios_print_nodisplay( f, "  FBOFF       Do not fix bug leading to missing or undefined sp3 parity\n" );
    inchi_ios_print_nodisplay( f, "  FB2OFF      Do not fix bugs found after v.1.02b release\n" );
    inchi_ios_print_nodisplay( f, "  FNUDOFF     Do not fix non-uniform drawing issues\n" );
#endif
#endif


/*
breleaseversion<1
#else
    inchi_ios_print_nodisplay( f, "%s ver %s. Special testing version 12-12-2002.\n", INCHI_NAME, INCHI_VERSION);
    inchi_ios_print_nodisplay( f, "\nUsage:\ncINChI09b inputFile [outputFile [logFile [problemFile]]] [%coption[ %coption...]]\n", INCHI_OPTION_PREFX, INCHI_OPTION_PREFX);
    inchi_ios_print_nodisplay( f, "\nOptions:\n");
    inchi_ios_print_nodisplay( f, "\tB        Basic\n");
    inchi_ios_print_nodisplay( f, "\tT        basic Tautomeric\n");
    inchi_ios_print_nodisplay( f, "\tI        Isotopic\n");
    inchi_ios_print_nodisplay( f, "\tN        Non-isotopic\n");
    inchi_ios_print_nodisplay( f, "\tS        Stereo\n");
    inchi_ios_print_nodisplay( f, "\tE        Exclude Stereo\n");
    inchi_ios_print_nodisplay( f, "\tD        Display the structures\n");
    inchi_ios_print_nodisplay( f, "\tALT      produce shorter ALTernative representation (Abc)\n");
    inchi_ios_print_nodisplay( f, "\tSCT      produce shorter connection table representation\n");
    inchi_ios_print_nodisplay( f, "\tXML      output in xml format\n");
    inchi_ios_print_nodisplay( f, "\tPLAIN    output in plain format\n");
    inchi_ios_print_nodisplay( f, "\tMERGE    Merge all MOLfiles from the input file into one compound\n");
    inchi_ios_print_nodisplay( f, "\tWnumber  time-out per structure in seconds, W0 means unlimited\n");
    inchi_ios_print_nodisplay( f, "\tFnumber  set display Font size, points\n");
    inchi_ios_print_nodisplay( f, "\tSREL     Relative Stereo\n");
    inchi_ios_print_nodisplay( f, "\tSRAC     Racemic Stereo\n");
    inchi_ios_print_nodisplay( f, "\tNOUUSB   Omit stereobonds if all are unknown/undefined\n");
    inchi_ios_print_nodisplay( f, "\tNOUUSC   Omit stereocenters if all are unknown/undefined\n");
    inchi_ios_print_nodisplay( f, "\tSS       Slow Stereo: do not use stereo equivalence\n");
    inchi_ios_print_nodisplay( f, "\tRS       Do not test for Redundant Stereo elements\n");
    inchi_ios_print_nodisplay( f, "\tPW       Save warning structures in the problems file\n");
    inchi_ios_print_nodisplay( f, "\tPGO      Save only all good structures in the problems file\n");
    inchi_ios_print_nodisplay( f, "\tDSB      Double Stereo Bonds only (ignore alternating bonds stereo)\n");
    inchi_ios_print_nodisplay( f, "\tRSB:n    Min Ring Size for detecting for Stereo Bonds (n=1 => all)\n");
    inchi_ios_print_nodisplay( f, "\tAUXINFO:0          do not output auxiliary information (default:1)\n");
    inchi_ios_print_nodisplay( f, "\tDISCONSALT:0       do not disconnect salts (default:1)\n");
    inchi_ios_print_nodisplay( f, "\tDISCONMETAL:0      do not disconnect metals (default:1)\n");
    inchi_ios_print_nodisplay( f, "\tDISCONMETALCHKVAL:1 do not disconnect if typical valence (default:0)\n");
    inchi_ios_print_nodisplay( f, "\tRECONMETAL:0       do not reconnect metals (default:1)\n");
    inchi_ios_print_nodisplay( f, "\tMOVEPOS:0          do not check moveable positive charges (default:1)\n");
    inchi_ios_print_nodisplay( f, "\tACIDTAUT:n         n=1: one H/(-) tautomerism, 2: more (deflt), 0:none\n");
    inchi_ios_print_nodisplay( f, "\tMERGESALTTG:1      merge salt t-groups (default), 0: do not merge\n");
    inchi_ios_print_nodisplay( f, "\tUNCHARGEDACIDS:1   Apply salt (acid) tautomerism in neutral species\n");
    inchi_ios_print_nodisplay( f, "\tO:[suffix]         Open all 4 files adding suffix to the inputFile name\n");
    inchi_ios_print_nodisplay( f, "\tOP:outputpath      Set output path\n");
    inchi_ios_print_nodisplay( f, "\tMOL                input file is a MOLfile (default)\n");
    inchi_ios_print_nodisplay( f, "\tSDF[:DataHeader]   Include SDfile data for the header into the results\n");
    inchi_ios_print_nodisplay( f, "\tSDFID              extract CAS r.n. in addition to requested SDfile data\n");
    inchi_ios_print_nodisplay( f, "\tSTART:number       Start at the given structure ordering number\n");
    inchi_ios_print_nodisplay( f, "\tEND:number         Terminate after the given structure ordering number\n");
*/
#endif

}




/*^^^ */
/************************************************************************************/
void HelpCommandLineParmsReduced( INCHI_IOSTREAM *f )
{
    if ( !f )
        return;

#if ( bRELEASE_VERSION == 1 )


    /*^^^ */
     inchi_ios_print_nodisplay( f, 
#ifdef TARGET_EXE_USING_API
         "%s ver %s%s.\n\nUsage:\ninchi_main inputFile [outputFile [logFile [problemFile]]] [%coption[ %coption...]]\n", 
         INCHI_NAME, INCHI_VERSION, TARGET_ID_STRING, 
         INCHI_OPTION_PREFX, INCHI_OPTION_PREFX); 
#else
         "%s ver %s%s.\n\nUsage:\nc%s-%s inputFile [outputFile [logFile [problemFile]]] [%coption[ %coption...]]\n", 
         INCHI_NAME, INCHI_VERSION, TARGET_ID_STRING, 
         INCHI_NAME, INCHI_VERSION, INCHI_OPTION_PREFX, INCHI_OPTION_PREFX); 		 
#endif
    /*^^^ */    

    inchi_ios_print_nodisplay( f, "\nOptions:\n");


    inchi_ios_print_nodisplay( f, "\nInput\n");
    inchi_ios_print_nodisplay( f, "  STDIO       Use standard input/output streams\n");
    inchi_ios_print_nodisplay( f, "  InpAux      Input structures in %s default aux. info format\n              (for use with STDIO)\n", INCHI_NAME);
    inchi_ios_print_nodisplay( f, "  SDF:DataHeader Read from the input SDfile the ID under this DataHeader\n");
#if ( ADD_CMLPP == 1 )
    inchi_ios_print_nodisplay( f, "  CML         Input in CML format (default for input file .CML extension)\n");
#endif
/*
    inchi_ios_print_nodisplay( f, "  START:n     Skip structures up to n-th one\n");
    inchi_ios_print_nodisplay( f, "  END:m       Skip structures after m-th one\n");
*/
#if ( BUILD_WITH_AMI == 1 )
    inchi_ios_print_nodisplay( f, "  AMI         Allow multiple input files (wildcards supported)\n");
#endif

    inchi_ios_print_nodisplay( f, "Output\n");
    inchi_ios_print_nodisplay( f, "  AuxNone     Omit auxiliary information (default: Include)\n");
    inchi_ios_print_nodisplay( f, "  SaveOpt     Save custom InChI creation options (non-standard InChI)\n");
    inchi_ios_print_nodisplay( f, "  NoLabels    Omit structure number, DataHeader and ID from %s output\n", INCHI_NAME);
    inchi_ios_print_nodisplay( f, "  Tabbed      Separate structure number, %s, and AuxInfo with tabs\n", INCHI_NAME);
    /* inchi_ios_print_nodisplay( f, "  Compress    Compressed output\n"); */
#if ( defined(_WIN32) && defined(_MSC_VER) && !defined(COMPILE_ANSI_ONLY) && !defined(TARGET_API_LIB) )
    inchi_ios_print_nodisplay( f, "  D           Display the structures\n");
    inchi_ios_print_nodisplay( f, "  EQU         Display sets of identical components\n");
    inchi_ios_print_nodisplay( f, "  Fnumber     Set display Font size in number of points\n");
#endif
    inchi_ios_print_nodisplay( f, "  OutputSDF   Convert %s created with default aux. info to SDfile\n", INCHI_NAME);
#if ( SDF_OUTPUT_DT == 1 )
    inchi_ios_print_nodisplay( f, "  SdfAtomsDT  Output Hydrogen Isotopes to SDfile as Atoms D and T\n");
#endif
#if ( BUILD_WITH_AMI == 1 )
    inchi_ios_print_nodisplay( f, "  AMIOutStd   Write output to stdout (in AMI mode)\n");
    inchi_ios_print_nodisplay( f, "  AMILogStd   Write log to stderr (in AMI mode)\n");
    inchi_ios_print_nodisplay( f, "  AMIPrbNone  Suppress creation of problem files (in AMI mode)\n");
#endif
    inchi_ios_print_nodisplay( f, "Structure perception\n");
    inchi_ios_print_nodisplay( f, "  SNon        Exclude stereo (default: include absolute stereo)\n");
    inchi_ios_print_nodisplay( f, "  NEWPSOFF    Both ends of wedge point to stereocenters (default: a narrow end)\n");
    inchi_ios_print_nodisplay( f, "  DoNotAddH   All H are explicit (default: add H according to usual valences)\n");
    
#ifndef USE_STDINCHI_API
    inchi_ios_print_nodisplay( f, "Stereo perception modifiers (non-standard InChI)\n");
    inchi_ios_print_nodisplay( f, "  SRel        Relative stereo\n");
    inchi_ios_print_nodisplay( f, "  SRac        Racemic stereo\n");
    inchi_ios_print_nodisplay( f, "  SUCF        Use Chiral Flag: On means Absolute stereo, Off - Relative\n"); 

    inchi_ios_print_nodisplay( f, "Customizing InChI creation (non-standard InChI)\n");
    inchi_ios_print_nodisplay( f, "  SUU         Always include omitted unknown/undefined stereo\n");
    inchi_ios_print_nodisplay( f, "  SLUUD       Make labels for unknown and undefined stereo different\n");
    inchi_ios_print_nodisplay( f, "  RecMet      Include reconnected metals results\n");
    inchi_ios_print_nodisplay( f, "  FixedH      Include Fixed H layer\n");
    inchi_ios_print_nodisplay( f, "  KET         Account for keto-enol tautomerism (experimental)\n");
    inchi_ios_print_nodisplay( f, "  15T         Account for 1,5-tautomerism (experimental)\n");
#endif

    inchi_ios_print_nodisplay( f, "Generation\n");
    inchi_ios_print_nodisplay( f, "  Wnumber     Set time-out per structure in seconds; W0 means unlimited\n");
    inchi_ios_print_nodisplay( f, "  WarnOnEmptyStructure Warn and produce empty %s for empty structure\n", INCHI_NAME);
    inchi_ios_print_nodisplay( f, "  Key         Generate InChIKey\n");
    inchi_ios_print_nodisplay( f, "  XHash1      Generate hash extension (to 256 bits) for 1st block of InChIKey\n");
    inchi_ios_print_nodisplay( f, "  XHash2      Generate hash extension (to 256 bits) for 2nd block of InChIKey\n");

#ifdef BUILD_WITH_ENG_OPTIONS
#if 0
    inchi_ios_print_nodisplay( f, "Engineering options (for testing only)\n");
    inchi_ios_print_nodisplay( f, "  NoADP       Disable Aggressive Deprotonation\n");
#if ( FIX_ADJ_RAD == 1 )
    inchi_ios_print_nodisplay( f, "  FixRad      Fix Adjacent Radicals\n");
#endif
    inchi_ios_print_nodisplay( f, "  SPXYZOFF    Do not include Phosphines Stereochemistry\n");
    inchi_ios_print_nodisplay( f, "  SAsXYZOFF   Do not include Arsines Stereochemistry\n");
    inchi_ios_print_nodisplay( f, "  FBOFF       Do not fix bug leading to missing or undefined sp3 parity\n" );
    inchi_ios_print_nodisplay( f, "  FB2OFF      Do not fix bugs found after v.1.02b release\n" );
    inchi_ios_print_nodisplay( f, "  FNUDOFF     Do not fix non-uniform drawing issues\n" );
#endif
#endif


      
#endif /* #if ( bRELEASE_VERSION == 1 ) */
}




#define fprintf2 inchi_fprintf

#ifndef TARGET_API_LIB
/************************************************************************************/
int OpenFiles( FILE **inp_file, FILE **output_file, FILE **log_file, FILE **prb_file, INPUT_PARMS *ip )
{
/*
  -- Files --
  ip->path[0] => Input
  ip->path[1] => Output (INChI)
  ip->path[2] => Log
  ip->path[3] => Problem structures
  ip->path[4] => Errors file (ACD Labs)

*/
    /*  logfile -- open as early as possible */
    if ( !ip->path[2] || !ip->path[2][0] ) {
        fprintf2( stderr, "%s version %s%s%s\n", INCHI_NAME, INCHI_VERSION, TARGET_ID_STRING, bRELEASE_VERSION? "":""); /* (Pre-release, for evaluation purposes only)" ); */
        fprintf2( stderr, "Log file not specified. Using standard error output.\n");
        *log_file = stderr;
    } else
    if ( !(*log_file = fopen( ip->path[2], "w" ) ) ) {
        fprintf2( stderr, "%s version %s%s%s\n", INCHI_NAME, INCHI_VERSION, TARGET_ID_STRING, bRELEASE_VERSION? "":""); /* (Pre-release, for evaluation purposes only)" );*/
        fprintf2( stderr, "Cannot open log file '%s'. Using standard error output.\n", ip->path[2] );
        *log_file = stderr;
    } else {
        fprintf2( *log_file, "%s version %s%s%s\n", INCHI_NAME, INCHI_VERSION, TARGET_ID_STRING, bRELEASE_VERSION? "":""); /* (Pre-release, for evaluation purposes only)" );*/
        fprintf2( *log_file, "Opened log file '%s'\n", ip->path[2] );
    }
    /* input file */
    if ( (ip->nInputType == INPUT_MOLFILE || ip->nInputType == INPUT_SDFILE ||
         ip->nInputType == INPUT_CMLFILE  || ip->nInputType == INPUT_INCHI  ||
         ip->nInputType == INPUT_INCHI_PLAIN ) && ip->num_paths > 0 ) 
    {
        const char *fmode = NULL;
#if ( defined(_MSC_VER)&&defined(_WIN32) || defined(__BORLANDC__)&&defined(__WIN32__) || defined(__GNUC__)&&defined(__MINGW32__)&&defined(_WIN32) )
        /* compilers that definitely allow fopen "rb" (binary read) mode */
        fmode = "rb";
        if ( !ip->path[0] || !ip->path[0][0] || !(*inp_file = fopen( ip->path[0], "rb" ) ) ) 
        {
            fprintf2( *log_file, "Cannot open input file '%s'. Terminating.\n", ip->path[0]? ip->path[0] : "<No name>" );
            goto exit_function;
        } 
        else 
        {
            if ( ip->nInputType == INPUT_CMLFILE ) 
            {
                int c;
#ifdef CML_DEBUG
                printf ("cr %d lf %d ret %d\n", (int) '\r', (int) '\f', (int) '\n');
#endif
                /* read up to the end of the first line */
                while( (c = fgetc( *inp_file )) && c != EOF && c != '\n' && c != '\r' )
                    ;
                if ( c == '\r' || c == EOF ) 
                {
                    /* text file contains CR; close and reopen as "text" */
                    fclose( *inp_file );
                    if ( !(*inp_file = fopen( ip->path[0], "r" ) ) ) 
                    {
                        fprintf2( *log_file, "Cannot open input file '%s' (2nd attempt). Terminating.\n", ip->path[0] );
                        goto exit_function;
                    }
                    fprintf2( *log_file, "Opened input file '%s'\n", ip->path[0] );
                    fmode = "r";
                } 
                else 
                {
                    fclose( *inp_file );
                    if ( !(*inp_file = fopen( ip->path[0], "rb" ) ) ) 
                    {
                        fprintf2( *log_file, "Cannot open input file '%s' (2nd attempt). Terminating.\n", ip->path[0] );
                        goto exit_function;
                    }
                    fprintf2( *log_file, "Opened input file '%s': no CR.\n", ip->path[0] );
                    fmode = "rb";
                }
            } /* CML */
            else
                fprintf2( *log_file, "Opened input file '%s'\n", ip->path[0] );

        }
#else
        if ( !ip->path[0] || !ip->path[0][0] || !(*inp_file = fopen( ip->path[0], "r" ) ) ) {
            fprintf2( *log_file, "Cannot open input file '%s'. Terminating.\n", ip->path[0]? ip->path[0] : "<No Name>" );
            goto exit_function;
        } else {
            fprintf2( *log_file, "Opened input file '%s'\n", ip->path[0] );
        }
        fmode = "r";
#endif
        DetectInputINChIFileType( inp_file, ip, fmode );
    } else
    if ( (ip->nInputType != INPUT_MOLFILE && ip->nInputType != INPUT_SDFILE && ip->nInputType != INPUT_CMLFILE && ip->nInputType != INPUT_INCHI 
        /*^^^ post-1.02b */
        && ip->nInputType != INPUT_INCHI_PLAIN
        )) {
        fprintf2( *log_file, "Input file type not specified. Terminating.\n");
        goto exit_function;
    } else {
        fprintf2( *log_file, "Input file not specified. Using standard input.\n");
        *inp_file = stdin;
    }
    /*  output file */
    if ( !ip->path[1] || !ip->path[1][0] ) {
        fprintf2( *log_file, "Output file not specified. Using standard output.\n");
        *output_file = stdout;
    } else {
        if ( !(*output_file = fopen( ip->path[1], "w" ) ) ) {
            fprintf2( *log_file, "Cannot open output file '%s'. Terminating.\n", ip->path[1] );
            goto exit_function;
        } else {
             fprintf2( *log_file, "Opened output file '%s'\n", ip->path[1] );
            if ( (ip->bINChIOutputOptions & (INCHI_OUT_PLAIN_TEXT)) &&
                  *inp_file != stdin &&
                  !(ip->bINChIOutputOptions & INCHI_OUT_SDFILE_ONLY) &&
                  !ip->bNoStructLabels &&
                  !(ip->bINChIOutputOptions & INCHI_OUT_TABBED_OUTPUT)) {
                 PrintFileName( "* Input_File: \"%s\"\n", *output_file, ip->path[0] );
             }
        }
    }
    /*  problem file */
    if ( ip->path[3] && ip->path[3][0] ) {
        const char *fmode = "w";
#if ( defined(_MSC_VER)&&defined(_WIN32) || defined(__BORLANDC__)&&defined(__WIN32__) || defined(__GNUC__)&&defined(__MINGW32__)&&defined(_WIN32) )
        if ( ip->nInputType != INPUT_CMLFILE ) {
            fmode = "wb";
        }
#endif
        if ( !(*prb_file = fopen( ip->path[3], fmode ) ) ) {
            fprintf2( *log_file, "Cannot open problem file '%s'. Terminating.\n", ip->path[3] );
            goto exit_function;
        } else {
             fprintf2( *log_file, "Opened problem file '%s'\n", ip->path[3] );
        }
    }
    return 1;  /*  success */

exit_function:
    return 0; /*  failed */

}
#define NUM_VERSIONS 7
#define LEN_VERSIONS 64
/*******************************************************************/
static int bMatchOnePrefix( int len, char *str, int lenPrefix[],
                            char strPrefix[][LEN_VERSIONS], int numPrefix);
/*******************************************************************/
static int bMatchOnePrefix( int len, char *str, int lenPrefix[],
                            char strPrefix[][LEN_VERSIONS], int numPrefix)
{
    int i;
    for ( i = 0; i < numPrefix; i ++ ) {
        if ( len >= lenPrefix[i] &&
             !memcmp( str, strPrefix[i], lenPrefix[i] ) ) { 
            return 1;
        }
    }
    return 0;
}
/*******************************************************************/
int DetectInputINChIFileType( FILE **inp_file, INPUT_PARMS *ip, const char *fmode )
{
    char szLine[256], ret = 0;
    static char szPlnVersion[NUM_VERSIONS][LEN_VERSIONS]; /* = "INChI:1.1Beta/";*/
    static int  lenPlnVersion[NUM_VERSIONS];
    static char szPlnAuxVer[NUM_VERSIONS][LEN_VERSIONS]; /* = "AuxInfo:1.1Beta/";*/
    static int  lenPlnAuxVer[NUM_VERSIONS];
    static char szXmlVersion[NUM_VERSIONS][LEN_VERSIONS]; /* = "<INChI version=\"1.1Beta\">";*/
    static int  lenXmlVersion[NUM_VERSIONS];
    static char szXmlStruct[LEN_VERSIONS]  = "<structure";
    static int  lenXmlStruct;
    static char szXmlIdentVer[NUM_VERSIONS][LEN_VERSIONS]; /*= "<identifier version=\"1.1Beta\"";*/
    static int  lenXmlIdentVer[NUM_VERSIONS];
    static char szXmlMsgError[LEN_VERSIONS];
    static int  lenXmlMsgError = 0;
    static char szXmlMsgFatal[LEN_VERSIONS]= "<message type=\"fatal (aborted)\"";
    static int  lenXmlMsgFatal;
    static int  bInitilized = 0;
    int  bINChI_plain = 0, bINChI_xml = 0, len, i;
    if ( ip->nInputType == INPUT_INCHI_XML || ip->nInputType == INPUT_INCHI_PLAIN || ip->nInputType == INPUT_INCHI ) {
        return 1;
    }
    if ( !bInitilized ) {
        lenPlnVersion[0]  = sprintf( szPlnVersion[0],  "%s=%s/", INCHI_NAME, INCHI_VERSION );
        lenPlnVersion[1]  = sprintf( szPlnVersion[1],  "INChI=1.12Beta/" );
        lenPlnVersion[2]  = sprintf( szPlnVersion[2],  "INChI=1.0RC/" );
        lenPlnVersion[3]  = sprintf( szPlnVersion[3],  "InChI=1.0RC/" );
        lenPlnVersion[4]  = sprintf( szPlnVersion[4],  "InChI=1/" );
        lenPlnVersion[5]  = sprintf( szPlnVersion[5],  "MoChI=1a/" );
        lenPlnVersion[6]  = sprintf( szPlnVersion[6],  "InChI=1S/" );
        lenPlnAuxVer[0]   = sprintf( szPlnAuxVer[0],   "AuxInfo=%s/", INCHI_VERSION );
        lenPlnAuxVer[1]   = sprintf( szPlnAuxVer[1],   "AuxInfo=1.12Beta/" );
        lenPlnAuxVer[2]   = sprintf( szPlnAuxVer[2],   "AuxInfo=1.0RC/" );
        lenPlnAuxVer[3]   = sprintf( szPlnAuxVer[3],   "AuxInfo=1.0RC/" );
        lenPlnAuxVer[4]   = sprintf( szPlnAuxVer[4],   "AuxInfo=1/" );
        lenPlnAuxVer[5]   = sprintf( szPlnAuxVer[5],   "AuxInfo=1a/" );
        lenPlnAuxVer[6]   = sprintf( szPlnAuxVer[6],   "AuxInfo=1/" );
        lenXmlVersion[0]  = sprintf( szXmlVersion[0],  "<%s version=\"%s\">", INCHI_NAME, INCHI_VERSION );
        lenXmlVersion[1]  = sprintf( szXmlVersion[1],  "<INChI version=\"1.12Beta\">" );
        lenXmlVersion[2]  = sprintf( szXmlVersion[2],  "<INChI version=\"1.0RC\">" );
        lenXmlVersion[3]  = sprintf( szXmlVersion[3],  "<InChI version=\"1.0RC\">" );
        lenXmlVersion[4]  = sprintf( szXmlVersion[4],  "<InChI version=\"1\">" );
        lenXmlVersion[5]  = sprintf( szXmlVersion[5],  "<MoChI version=\"1a\">" );
        lenXmlVersion[6]  = sprintf( szXmlVersion[6],  "<InChI version=\"1S\">" );
        lenXmlIdentVer[0] = sprintf( szXmlIdentVer[0], "<identifier version=\"%s\"", INCHI_VERSION );
        lenXmlIdentVer[1] = sprintf( szXmlIdentVer[1], "<identifier version=\"1.12Beta\"" );
        lenXmlIdentVer[2] = sprintf( szXmlIdentVer[2], "<identifier version=\"1.0RC\"" );
        lenXmlIdentVer[3] = sprintf( szXmlIdentVer[3], "<identifier version=\"1.0RC\"" );
        lenXmlIdentVer[4] = sprintf( szXmlIdentVer[4], "<identifier version=\"1\"" );
        lenXmlIdentVer[5] = sprintf( szXmlIdentVer[5], "<identifier version=\"1a\"" );
        lenXmlIdentVer[6] = sprintf( szXmlIdentVer[6], "<identifier version=\"1S\"" );
        lenXmlMsgError    = sprintf( szXmlMsgError,    "<message type=\"error (no %s)\"", INCHI_NAME );
        lenXmlStruct      = strlen(szXmlStruct);
        lenXmlMsgFatal    = strlen(szXmlMsgFatal);
#if ( FIX_DALKE_BUGS == 1 )
        bInitilized = 1;
#endif
    }
    for ( i = 0; i < 4; i ++ ) {
        len = inchi_fgetsLfTab( szLine, sizeof(szLine)-1, *inp_file );
        if ( len < 0 )
            break;
        if ( bMatchOnePrefix( len, szLine, lenPlnVersion, szPlnVersion, NUM_VERSIONS ) ||
             bMatchOnePrefix( len, szLine, lenPlnAuxVer, szPlnAuxVer, NUM_VERSIONS ) ) {
            bINChI_plain ++;
        } else
        if ( bMatchOnePrefix( len, szLine, lenXmlVersion, szXmlVersion, NUM_VERSIONS ) ||
             bMatchOnePrefix( len, szLine, &lenXmlStruct, &szXmlStruct, 1 ) ||
             bMatchOnePrefix( len, szLine, lenXmlIdentVer, szXmlIdentVer, NUM_VERSIONS ) ||
             bMatchOnePrefix( len, szLine, &lenXmlMsgError, &szXmlMsgError, 1 ) ||
             bMatchOnePrefix( len, szLine, &lenXmlMsgFatal, &szXmlMsgFatal, 1 ) ) {
            bINChI_xml ++;
        }
    }
    if ( bINChI_plain >= 2 && !bINChI_xml ) {
        ip->nInputType = INPUT_INCHI_PLAIN;
        ret = 1;
    } else
    if ( !bINChI_plain && bINChI_xml >= 3 ) {
        ip->nInputType = INPUT_INCHI_XML;
        ret = 1;
    }
/*exit_function:*/
    fclose ( *inp_file );
    *inp_file = fopen( ip->path[0], fmode );
    return ret;
}
#undef NUM_VERSIONS
#undef LEN_VERSIONS

#endif /* TARGET_API_LIB */

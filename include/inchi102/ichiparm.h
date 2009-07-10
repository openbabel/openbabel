/*
 * International Chemical Identifier (InChI)
 * Version 1
 * Software version 1.02
 * November 30, 2008
 * Developed at NIST
 *
 * The InChI library and programs are free software developed under the
 * auspices of the International Union of Pure and Applied Chemistry (IUPAC);
 * you can redistribute this software and/or modify it under the terms of 
 * the GNU Lesser General Public License as published by the Free Software 
 * Foundation:
 * http://www.opensource.org/licenses/lgpl-license.php
 */


int DetectInputINChIFileType( FILE **inp_file, INPUT_PARMS *ip, const char *fmode );
void HelpCommandLineParmsReduced( INCHI_IOSTREAM *f );

/*****************************************************************************************/
int ReadCommandLineParms( int argc, const char *argv[], INPUT_PARMS *ip, char *szSdfDataValue,
                          unsigned long *ulDisplTime, int bReleaseVersion, INCHI_IOSTREAM *log_file )
{
    int i, k, c;
    const char *q;
    unsigned long ul;
    int           nFontSize    = -9;
    int           nMode        = 0;
    int           nReleaseMode = nMode | (REQ_MODE_BASIC | REQ_MODE_TAUT | REQ_MODE_ISO | REQ_MODE_STEREO);
#if( MIN_SB_RING_SIZE > 0 )
    int           nMinDbRinSize = MIN_SB_RING_SIZE, mdbr=0;
#endif
    /*int           bNotRecognized=0;*/
    char          szNameSuffix[32], szOutputPath[512];
    int           bNameSuffix, bOutputPath;
    int           bMergeAllInputStructures;
    int           bDisconnectSalts       = (DISCONNECT_SALTS==1);
    int           bDoNotAddH             = 0;

#ifdef INCHI_LIB
    int           bVer1Options           = 0;
    int           bReconnectCoord        = 1;
    int           bDisconnectCoord       = 1;
    int           bINChIOutputOptions     = INCHI_OUT_EMBED_REC; /* embed reconnected & output full aux info */
    int           bCompareComponents     = CMP_COMPONENTS;
#else
    int           bVer1Options           = 1;
    int           bReconnectCoord        = (RECONNECT_METALS==1);
    int           bDisconnectCoord       = (DISCONNECT_METALS==1);
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

#define VER100_DEFAULT_MODE    (REQ_MODE_TAUT | REQ_MODE_ISO | REQ_MODE_STEREO |\
                                REQ_MODE_SB_IGN_ALL_UU | REQ_MODE_SC_IGN_ALL_UU)

    INCHI_MODE     bVer1DefaultMode       = VER100_DEFAULT_MODE;
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
#if( READ_INCHI_STRING == 1 )
    int           bDisplayIfRestoreWarnings = 0;
#endif
#ifdef INCHI_LIB
    int           bXml = INCHI_OUT_XML;
#else
    int           bXml = INCHI_OUT_PLAIN_TEXT;
#endif
    int bTgFlagVariableProtons = 1;
    int bTgFlagHardAddRenProtons = 1;
#ifdef STEREO_WEDGE_ONLY
    int bPointedEdgeStereo = STEREO_WEDGE_ONLY; /*   NEWPS TG_FLAG_POINTED_EDGE_STEREO*/
#endif
#if( FIX_ADJ_RAD == 1 )
    int bFixAdjacentRad = 0;
#endif
#if( ADD_PHOSPHINE_STEREO == 1 )
    int bAddPhosphineStereo = 0;
#endif
#if( ADD_ARSINE_STEREO == 1 )
    int bAddArsineStereo = 0;
#endif

    int bFixSp3bug        = 0; 
    /* post v.1 features */
    int bKetoEnolTaut     = 0;
    int b15TautNonRing    = 0;

    int bFixFB2           = 0; 

    ext[0] = ".mol";
    ext[1] = bVer1Options? ".txt" : ".ich";
    ext[2] = ".log";
    ext[3] = ".prb";
    ext[4] = "";

#if( MAX_NUM_PATHS < 4 )
  #error Wrong initialization
#endif

    /*  init table parms */
    memset ( ip, 0, sizeof(*ip) );

#ifdef INCLUDE_V102_FEATURES
    
    ip->bFixNonUniformDraw = 0;         /* "FNUDOFF" */

#ifdef STDINCHI_AS_DEFAULT
    /* set defaults according to standard InChI generation options */
    bVer1DefaultMode &= ~REQ_MODE_BASIC;/* "FIXEDHOFF" */
    bReconnectCoord = 0;                /* "RECMETOFF" */
    bPointedEdgeStereo = 1;             /* "NEWPS" */
    bAddPhosphineStereo = 1;            /* "SPXYZ" */
    bAddArsineStereo = 1;               /* "SASXYZ" */ 
    bFixSp3bug = 1;                     /* "FB" */	
    bFixFB2 = 1;						/* "FB2" */
    ip->bFixNonUniformDraw = 1;         /* "FNUD" */
#endif

#endif

#ifndef INCHI_ANSI_ONLY
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
#ifdef INCHI_LIBRARY
    ip->msec_MaxTime = 0;      /*  milliseconds, default = unlimited in libinchi */
#else
    ip->msec_MaxTime = 60000;  /*  milliseconds, default = 60 sec */
#endif


    if ( bReleaseVersion ) {

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
    } else {
        /*bXml = INCHI_OUT_PLAIN_TEXT;*/
        nReleaseMode = 0;
    }

    if ( bVer1Options ) {
        bNameSuffix = 1;
        szNameSuffix[0] = '\0';
    }


#if( ACD_LABS_VERSION == 1 ) /* { */
/*
  -- Files --
  ip->path[0] => Input               -I
  ip->path[1] => Output (INChI)       -O
  ip->path[2] => Log                 -L
  ip->path[3] => Problem structures
  ip->path[4] => Errors file (ACD(   -E

*/
        for ( i = 1; i < argc; i ++ ) {
            if ( *argv[i] !memicmp( argv[i], "-I", 2) ) {
                ip->num_paths += !ip->path[0];
                ip->path[0] = _strdup( argv[i] + 2 );
            } else
            if ( *argv[i] !memicmp( argv[i], "-O", 2) ) {
                ip->num_paths += !ip->path[1];
                ip->path[1] = _strdup( argv[i] + 2 );
            } else
            if ( *argv[i] !memicmp( argv[i], "-L", 2) ) {
                ip->num_paths += !ip->path[2];
                ip->path[2] = _strdup( argv[i] + 2 );
            } else
            if ( *argv[i] !memicmp( argv[i], "-E", 2) ) {
                ip->path[4] = _strdup( argv[i] + 2 );
            } else
            if ( *argv[i] !stricmp( argv[i], "-Z" ) ) {
                sprintf( stdout, "%s version %s\n", INCHI_NAME, INCHI_VERSION );
                return -1;
            }
        }

#else /* } ACD_LABS_VERSION { */

    for ( i = 1; i < argc; i ++ ) {
        if ( bVer1Options && INCHI_OPTION_PREFX == argv[i][0] && ':' == argv[i][1] && !argv[i][2] ) {
            bVer1Options &= ~1;  /* turn off ver 1 options mode on "/:" (Win) or "-:" (Linux) */
        } else
        if ( !(bVer1Options & 1) && INCHI_OPTION_PREFX == argv[i][0] && INCHI_OPTION_PREFX != argv[i][1] )
        {
            /*============== Version 0.9xx Beta options & INCHI_LIB ==================*/
            /***************************************************************/
            /*                                                             */
            /*       Version 0.9xx Beta  and INCHI_LIB (GUI) options       */
            /*                                                             */
            /***************************************************************/

            pArg = argv[i]+1;

            /* parameter */
            if ( INPUT_NONE == ip->nInputType && !stricmp( pArg, "MOL" ) ) {
                ip->nInputType = INPUT_MOLFILE;
            } else
            if ( INPUT_NONE == ip->nInputType && !stricmp( pArg, "SDF" ) ) {
                ip->nInputType = INPUT_MOLFILE;
            } else
            if ( INPUT_NONE == ip->nInputType &&
                (!memicmp( pArg, "SDF", 3 )) &&
                ( pArg[3] == ':' ) ) {
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
/*^^^ post-1.02b */
            } else 
            if ( !stricmp( pArg, "INPAUX" ))
            {
                if (INPUT_NONE == ip->nInputType) 
                    ip->nInputType = INPUT_INCHI_PLAIN;
            } else
#if( ADD_CMLPP == 1 )
            if ( INPUT_NONE == ip->nInputType && !stricmp( pArg, "CML" )  ) {
                 /* CMLfile label */
                ip->nInputType = INPUT_CMLFILE;
            } else
#endif
            if ( !memicmp( pArg, "START:", 6 ) ) {
                ip->first_struct_number = strtol(pArg+6, NULL, 10);
            } else
            if ( !memicmp( pArg, "END:", 4 ) ) {
                ip->last_struct_number = strtol(pArg+4, NULL, 10);
            } else /*  RSB: */
            if ( !memicmp( pArg, "RSB:", 4 ) ) {
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
            } else
            if ( !memicmp( pArg, "O:", 2 ) ) {
                bNameSuffix = 1;
                strncpy(szNameSuffix, pArg+2, sizeof(szNameSuffix)-1);
            } else
            if ( !memicmp( pArg, "OP:", 3 ) ) {
                bOutputPath = 1;
                strncpy(szOutputPath, pArg+3, sizeof(szOutputPath)-1);
            } else
            if ( !stricmp( pArg, "ALT" ) ) {
                ip->bAbcNumbers = 1;
            } else
            if ( !stricmp( pArg, "SCT" ) ) {
                ip->bCtPredecessors = 1;
            } else

            if ( !stricmp( pArg, "CMP" ) ) {
                bCompareComponents = CMP_COMPONENTS;
            } else
            if ( !stricmp( pArg, "CMPNONISO" ) ) {
                bCompareComponents = CMP_COMPONENTS | CMP_COMPONENTS_NONISO;
            } else
            if ( !stricmp( pArg, "SREL" ) ) {
                if ( nMode & REQ_MODE_RACEMIC_STEREO ) {
                    nMode ^= REQ_MODE_RACEMIC_STEREO;
                }
                if ( nMode & REQ_MODE_CHIR_FLG_STEREO ) {
                    nMode ^= REQ_MODE_CHIR_FLG_STEREO;
                }
                nMode |= REQ_MODE_RELATIVE_STEREO;
                nMode |= REQ_MODE_STEREO;
            } else
            if ( !stricmp( pArg, "SRAC" ) ) {
                if ( nMode & REQ_MODE_RELATIVE_STEREO ) {
                    nMode ^= REQ_MODE_RELATIVE_STEREO;
                }
                if ( nMode & REQ_MODE_CHIR_FLG_STEREO ) {
                    nMode ^= REQ_MODE_CHIR_FLG_STEREO;
                }
                nMode |= REQ_MODE_RACEMIC_STEREO;
                nMode |= REQ_MODE_STEREO;
            } else
            if ( !stricmp( pArg, "SUCF" ) ) {
                if ( nMode & REQ_MODE_RELATIVE_STEREO ) {
                    nMode ^= REQ_MODE_RELATIVE_STEREO;
                }
                if ( nMode & REQ_MODE_RACEMIC_STEREO ) {
                    nMode ^= REQ_MODE_RACEMIC_STEREO;
                }
                nMode |= REQ_MODE_CHIR_FLG_STEREO;    /* stereo defined by the Chiral flag */
                nMode |= REQ_MODE_STEREO;
            } else
            if ( !stricmp( pArg, "ChiralFlagON" ) ) { /* used only with /SUCF */
                bForcedChiralFlag &= ~FLAG_SET_INP_AT_NONCHIRAL;
                bForcedChiralFlag |= FLAG_SET_INP_AT_CHIRAL;
            } else
            if ( !stricmp( pArg, "ChiralFlagOFF" ) ) { /* used only with /SUCF */
                bForcedChiralFlag &= ~FLAG_SET_INP_AT_CHIRAL;
                bForcedChiralFlag |= FLAG_SET_INP_AT_NONCHIRAL;
            } else
            if ( !stricmp( pArg, "NOUUSB" ) ) {
                nMode |= REQ_MODE_SB_IGN_ALL_UU;
            } else
            if ( !stricmp( pArg, "NOUUSC" ) ) {
                nMode |= REQ_MODE_SC_IGN_ALL_UU;
            } else
            if ( !stricmp( pArg, "SDFID" ) ) {
                ip->bGetSdfileId = 1;
            } else
            if ( !stricmp( pArg, "XML" ) ) {
                bXml &= ~INCHI_OUT_PLAIN_TEXT;
                bXml |=  INCHI_OUT_XML;
                /*bNotRecognized = 2*bReleaseVersion;*/
            } else
            if ( !stricmp( pArg, "PLAIN" ) ) {
                bXml |=  INCHI_OUT_PLAIN_TEXT;
                bXml &= ~INCHI_OUT_XML;
            } else
#if( !defined(INCHI_LIBRARY) && !defined(INCHI_LIB) )
            if ( !stricmp( pArg, "Tabbed" ) ) {
                bXml |=  INCHI_OUT_TABBED_OUTPUT;
            } else
#endif
            if ( !stricmp( pArg, "ANNPLAIN" ) ) {
                bXml |=   INCHI_OUT_PLAIN_TEXT_COMMENTS;
                bXml &=  ~INCHI_OUT_XML_TEXT_COMMENTS;
            } else
            if ( !stricmp( pArg, "ANNXML" ) ) {
                bXml |=   INCHI_OUT_XML_TEXT_COMMENTS;
                bXml &=  ~INCHI_OUT_PLAIN_TEXT_COMMENTS;
            } else
            if ( !stricmp( pArg, "DONOTADDH" ) ) {
                bDoNotAddH = 1;
            } else
            if ( !memicmp( pArg, "AUXINFO:", 8 ) && isdigit(UCINT pArg[8]) ) {
                k = strtol(pArg+8, NULL, 10);
                if ( k == 0 ) {
                    bINChIOutputOptions |= INCHI_OUT_NO_AUX_INFO;  /* no aux info */
                    bINChIOutputOptions &= ~INCHI_OUT_SHORT_AUX_INFO;
                } else
                if ( k == 1 ) {
                    bINChIOutputOptions &= ~(INCHI_OUT_NO_AUX_INFO | INCHI_OUT_SHORT_AUX_INFO); /* include full aux info */
                } else
                if ( k == 2 ) {
                    bINChIOutputOptions &= ~INCHI_OUT_NO_AUX_INFO; /* include short aux info */
                    bINChIOutputOptions |= INCHI_OUT_SHORT_AUX_INFO;
                } else {
                    bINChIOutputOptions = k;  /* override everything */
                }
            } else
            if ( !stricmp( pArg, "MERGE" ) ) {
                bMergeAllInputStructures = 1;
            } else
            if ( !stricmp( pArg, "PGO" ) ) {
                ip->bSaveAllGoodStructsAsProblem = 1;
            } else
            if ( !stricmp( pArg, "DCR" ) ) {
                bDisplayCompositeResults = 1;
            } else
            if ( !stricmp( pArg, "DSB" ) ) {
                nMode |= REQ_MODE_NO_ALT_SBONDS;
            } else
            if ( !stricmp( pArg, "NOHDR" ) ) {
                 bNoStructLabels = 1;
            } else
            if ( !stricmp( pArg, "NoVarH" ) ) {
                 bTgFlagVariableProtons = 0;
            } else

            #ifndef ONLY_STDINCHI_STDKEY
            if ( !stricmp( pArg, "NoADP" ) ) {
                 bTgFlagHardAddRenProtons = 0;
            } else
            if ( !stricmp( pArg, "FixSp3bug" ) ) {
                 bFixSp3bug = 1;
            } else
            if ( !stricmp( pArg, "FB" ) ) {
                 bFixSp3bug = 1; /* fix all known v1 bugs */
            } else
            /*^^^ */
            if ( !stricmp( pArg, "FB2" ) ) {
                 bFixFB2 = 1;
            } else
            /*^^^ */
            if ( !stricmp( pArg, "FNUD" ) ) {
                 ip->bFixNonUniformDraw = 1;
            } else
#ifdef STEREO_WEDGE_ONLY
            if ( !stricmp( pArg, "NEWPS" ) ) {
                 bPointedEdgeStereo = 1;
            } else
#endif
#if( ADD_PHOSPHINE_STEREO == 1 )
            if ( !stricmp( pArg, "SPXYZ" ) ) {
                 bAddPhosphineStereo = 1;
            } else
#endif
#if( ADD_ARSINE_STEREO == 1 )
            if ( !stricmp( pArg, "SASXYZ" ) ) {
                 bAddArsineStereo = 1;
            } else
#endif
            #endif /* ifndef ONLY_STDINCHI_STDKEY */

            if ( !stricmp( pArg, "Key" ) ) {
                ip->bCalcInChIKey = 1;
            } else

/*^^^ */
#ifdef INCLUDE_V102_FEATURES 

            /* Check if turning OFF entire "standard" InChI generation shortcut is requested */
            #ifndef ONLY_STDINCHI_STDKEY
            if ( !stricmp( pArg, "STDINCHIOFF" ) ) {              
                bVer1DefaultMode &= ~REQ_MODE_BASIC;/* FIXEDHOFF */
                bReconnectCoord = 0;				/* RECMETOFF */
                bAddPhosphineStereo = 0;			/* "SPXYZOFF" */
                bAddArsineStereo = 0;				/* "SASXYZOFF" */
                bPointedEdgeStereo = 0;				/* "NEWPSOFF" */
                bFixSp3bug = 0;						/* "FBOFF" */
                bFixFB2 = 0;						/* "FB2OFF" */
                ip->bFixNonUniformDraw = 0;			/* "FNUDOFF" */
            } else
            #endif

            /* Check if "standard" InChI generation is requested via STDINCHI shortcut */
            #ifndef ONLY_STDINCHI_STDKEY
            if ( !stricmp( pArg, "STDINCHI" ) ) 
            {
                /* Set defaults according to "standard" InChI generation options */
                
                /* "SABS" (it is default but STDINCHI may me placed after e.g. SRAC) */
                bVer1DefaultMode |= REQ_MODE_STEREO;  
                nMode &= ~(REQ_MODE_RACEMIC_STEREO | REQ_MODE_RELATIVE_STEREO | REQ_MODE_CHIR_FLG_STEREO);
                bVer1DefaultMode |= REQ_MODE_BASIC;	/* "FIXEDH" */
                bReconnectCoord = 1;                /* "RECMET" */
                bPointedEdgeStereo = 1;             /* "NEWPS" */
                bAddPhosphineStereo = 1;            /* "SPXYZ" */
                bAddArsineStereo = 1;               /* "SASXYZ" */ 
                bFixSp3bug = 1;                     /* "FB" */				
                bFixFB2 = 1;						/* "FB2" */
                ip->bFixNonUniformDraw = 1;         /* "FNUD" */
            } else
            #endif

/*#endif*/ /* STDINCHI_AS_DEFAULT */

            /* Check if turning OFF "standard" InChI generation options is requested */
            if ( !stricmp( pArg, "NEWPSOFF" ) ) {
                 bPointedEdgeStereo = 0;
            } else
            #ifndef ONLY_STDINCHI_STDKEY
            if ( !stricmp( pArg, "FIXEDHOFF" ) ) {
                 bVer1DefaultMode &= ~REQ_MODE_BASIC; /* only taut. */
            } else
            if ( !stricmp( pArg, "RECMETOFF" ) ) {
                 bReconnectCoord = 0;    
            } else
            if ( !stricmp( pArg, "SPXYZOFF" ) ) {
                 bAddPhosphineStereo = 0;
            } else
            if ( !stricmp( pArg, "SASXYZOFF" ) ) {
                 bAddArsineStereo = 0;
            } else
            if ( !stricmp( pArg, "FixSp3bugOFF" ) ) {
                 bFixSp3bug = 0;
            } else 
            if ( !stricmp( pArg, "FBOFF" ) ) {
                 bFixSp3bug = 0; 
            } else 
            if ( !stricmp( pArg, "FB2OFF" ) ) {
                 bFixFB2 = 0; 
            } else
            if ( !stricmp( pArg, "FNUDOFF" ) ) {
                 ip->bFixNonUniformDraw = 0;
            } else
            #endif
                
#endif /* INCLUDE_V102_FEATURES */


            if ( !stricmp( pArg, "PW" ) ) {
                ip->bSaveWarningStructsAsProblem = 1;
            } else {
                /*for ( k = 0; c=pArg[k]; k ++ )*/
                k = 0;
                c=pArg[k]; /* prohibit multiple option concatenations, strict syntax check 2008-11-05 DT  */ 
                {
                    c = toupper( c );
                    switch ( c ) {
                    case 'D':
                        bDisplay |= 1;
                        if ( (pArg[k+1] == 'C' || pArg[k+1] == 'c') && !pArg[k+2] ) {
                            bDisplay |= 1;
                            k++;
                            ip->bDisplayEachComponentINChI = 1;
                        } else
                        if ( !pArg[k+1] ) {
                            bDisplay |= 1;
                        }
                        break;
                    case 'W':
                        if ( pArg[k+1] == 'D' ) {  /* restore Display Time functionality */
                            c = 'D';
                            k ++;
                        }
                        t = strtod( pArg+k+1, (char**)&q ); /*  cast deliberately discards 'const' qualifier */
                        if ( q > pArg+k+1 && errno == ERANGE || t < 0.0 || t*1000.0 > (double)ULONG_MAX) {
                            ul = 0;
                        } else {
                            ul = (unsigned long)(t*1000.0);
                        }
                        if ( /*q > pArg+k &&*/ !*q ) {
                            k = q - pArg - 1; /* k will be incremented by the for() cycle */
                            switch( c ) {
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
                        if ( q > pArg+k && !*q ) {
                            k = q - pArg - 1;
                            if ( abs(c) > 5 ) {
                                nFontSize = -c;  /* font size 5 or less is too small */
                            }
                        }
                        break;
                    default:
                        if ( !pArg[k+1] ) {
                            switch ( c ) {
                            case 'B':
                                nMode |= REQ_MODE_BASIC;
                                nReleaseMode = 0;
                                /*bNotRecognized = bReleaseVersion;*/
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
                                /*bNotRecognized = bReleaseVersion;*/
                                break;
                            case 'S':
                                nMode |= REQ_MODE_STEREO;
                                nReleaseMode = 0;
                                /*bNotRecognized = bReleaseVersion;*/
                                break;
                            case 'E':
                                if ( nReleaseMode & REQ_MODE_STEREO ) {
                                    nReleaseMode ^= REQ_MODE_STEREO;
                                }
                                break;
#ifndef INCHI_LIB
                            default:
#ifdef ONLY_STDINCHI_STDKEY
                                inchi_ios_eprint(log_file, "Terminating: there is no option \"%c\" in standard InChI software library\n.", c);
                                return -1;
#else
                                inchi_ios_eprint(log_file, "Unrecognized option: \"%c\".\n", c);
#endif

#endif
                            }
                        } 
#ifndef INCHI_LIB
                        else 
                        {
#ifdef ONLY_STDINCHI_STDKEY
                            inchi_ios_eprint(log_file, "Terminating: there is no option \"%c\" in standard InChI software library\n", c);
                            return -1;
#else
                            inchi_ios_eprint(log_file, "Unrecognized option: \"%c\".\n", c);
#endif
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
            if ( bNotRecognized && bNotRecognized == 2*bReleaseVersion ) {
               inchi_ios_eprint(stderr, "Unrecognized option: \"%s\".\n", argv[i]);
               bNotRecognized = 0;
            }
            */
        } else
        if ( (bVer1Options & 1) && INCHI_OPTION_PREFX == argv[i][0] && argv[i][1] ) {
            /***************************************************************/
            /*                                                             */
            /*                    Version 1.00 Beta options                */
            /*                                                             */
            /***************************************************************/
            bRecognizedOption = 2;
            pArg = argv[i] + 1;
#ifdef CML_DEBUG
            printf ("1 argv %d %s\n", i, argv[i]);
#endif

            bVer1Options += 2;
            /* always on: REQ_MODE_TAUT | REQ_MODE_ISO | REQ_MODE_STEREO */
            
            #ifndef ONLY_STDINCHI_STDKEY
            if ( !stricmp( pArg, "FIXEDH" ) ) {
                bVer1DefaultMode |= REQ_MODE_BASIC;  /* If tautomeric then tautomeric only */
            } else
            #endif

            if ( !stricmp( pArg, "SNON" ) ) {
                bVer1DefaultMode &= ~REQ_MODE_STEREO; /* no stereo */
                nMode &= ~(REQ_MODE_RACEMIC_STEREO | REQ_MODE_RELATIVE_STEREO | REQ_MODE_CHIR_FLG_STEREO);
            } else

            #ifndef ONLY_STDINCHI_STDKEY
                if ( (!stricmp( pArg, "SABS" ) )) {
                bVer1DefaultMode |= REQ_MODE_STEREO;  /* abs stereo (default) */
                nMode &= ~(REQ_MODE_RACEMIC_STEREO | REQ_MODE_RELATIVE_STEREO | REQ_MODE_CHIR_FLG_STEREO);
            } else
            if ( !stricmp( pArg, "SREL" ) ) {
                bVer1DefaultMode |= REQ_MODE_STEREO;  /* relative stereo */
                nMode &= ~(REQ_MODE_RACEMIC_STEREO | REQ_MODE_CHIR_FLG_STEREO);
                nMode |= REQ_MODE_RELATIVE_STEREO;
            } else  /* REQ_MODE_CHIR_FLG_STEREO */
            if ( !stricmp( pArg, "SRAC" ) ) {
                bVer1DefaultMode |= REQ_MODE_STEREO;  /* racemic stereo */
                nMode &= ~(REQ_MODE_RELATIVE_STEREO | REQ_MODE_CHIR_FLG_STEREO);
                nMode |= REQ_MODE_RACEMIC_STEREO;
            } else
            if ( !stricmp( pArg, "SUCF" ) ) {
                bVer1DefaultMode |= REQ_MODE_STEREO;  /* stereo defined by the Chiral flag */
                nMode &= ~(REQ_MODE_RELATIVE_STEREO | REQ_MODE_RACEMIC_STEREO);
                nMode |= REQ_MODE_CHIR_FLG_STEREO;
            } else
            if ( !stricmp( pArg, "ChiralFlagON" ) ) { /* used only with /SUCF */
                bForcedChiralFlag &= ~FLAG_SET_INP_AT_NONCHIRAL;
                bForcedChiralFlag |= FLAG_SET_INP_AT_CHIRAL;
            } else
            if ( !stricmp( pArg, "ChiralFlagOFF" ) ) { /* used only with /SUCF */
                bForcedChiralFlag &= ~FLAG_SET_INP_AT_CHIRAL;
                bForcedChiralFlag |= FLAG_SET_INP_AT_NONCHIRAL;
            } else
            if ( !stricmp( pArg, "SUU" ) ) {       /* include omitted undef/unkn stereo */
                bVer1DefaultMode &= ~(REQ_MODE_SB_IGN_ALL_UU | REQ_MODE_SC_IGN_ALL_UU);
            } else
            if ( !stricmp( pArg, "FixSp3bug" ) ) {
                 bFixSp3bug = 1;
            } else
            if ( !stricmp( pArg, "FB" ) ) {
                 bFixSp3bug = 1; /* fix all known v1 bugs */
            } else
            /*^^^ */
            if ( !stricmp( pArg, "FB2" ) ) {
                 bFixFB2 = 1; 
            } else
            if ( !stricmp( pArg, "FNUD" ) ) {
                 ip->bFixNonUniformDraw = 1;
            } else
#ifdef STEREO_WEDGE_ONLY
            if ( !stricmp( pArg, "NEWPS" ) ) {
                 bPointedEdgeStereo = 1;
            } else
#endif
#if( ADD_PHOSPHINE_STEREO == 1 )
            if ( !stricmp( pArg, "SPXYZ" ) ) {
                 bAddPhosphineStereo = 1;
            } else
#endif
#if( ADD_ARSINE_STEREO == 1 )
            if ( !stricmp( pArg, "SASXYZ" ) ) {
                 bAddArsineStereo = 1;
            } else
#endif
            #endif /* ifndef ONLY_STDINCHI_STDKEY */

            if ( !stricmp( pArg, "Key" ) ) {
                ip->bCalcInChIKey = 1;
            } else


/*^^^ */
#ifdef INCLUDE_V102_FEATURES 
            /* Check if turning OFF entire "standard" InChI generation shortcut is requested */
            #ifndef ONLY_STDINCHI_STDKEY
            if ( !stricmp( pArg, "STDINCHIOFF" ) ) {              
                bVer1DefaultMode &= ~REQ_MODE_BASIC;/* FIXEDHOFF */
                bReconnectCoord = 0;				/* RECMETOFF */
                bAddPhosphineStereo = 0;			/* "SPXYZOFF" */
                bAddArsineStereo = 0;				/* "SASXYZOFF" */
                bPointedEdgeStereo = 0;				/* "NEWPSOFF" */
                bFixSp3bug = 0;						/* "FBOFF" */
                bFixFB2 = 0;						/* "FB2OFF" */
                ip->bFixNonUniformDraw = 0;			/* "FNUDOFF" */
            } else
            #endif
            /* Check if "standard" InChI generation is requested via STDINCHI shortcut */
            #ifndef ONLY_STDINCHI_STDKEY
            if ( !stricmp( pArg, "STDINCHI" ) ) 
            {
                /* Set defaults according to "standard" InChI generation options */
                
                /* "SABS" (it is default but STDINCHI may me placed after e.g. SRAC) */
                bVer1DefaultMode |= REQ_MODE_STEREO;  
                nMode &= ~(REQ_MODE_RACEMIC_STEREO | REQ_MODE_RELATIVE_STEREO | REQ_MODE_CHIR_FLG_STEREO);
                bVer1DefaultMode |= REQ_MODE_BASIC;	/* "FIXEDH" */
                bReconnectCoord = 1;                /* "RECMET" */
                bPointedEdgeStereo = 1;             /* "NEWPS" */
                bAddPhosphineStereo = 1;            /* "SPXYZ" */
                bAddArsineStereo = 1;               /* "SASXYZ" */ 
                bFixSp3bug = 1;                     /* "FB" */				
                bFixFB2 = 1;						/* "FB2" */
                ip->bFixNonUniformDraw = 1;         /* "FNUD" */
            } else
            #endif

            /* Check if turning OFF "standard" InChI generation options is requested */
            if ( !stricmp( pArg, "NEWPSOFF" ) ) {
                 bPointedEdgeStereo = 0;
            } else
            #ifndef ONLY_STDINCHI_STDKEY
            if ( !stricmp( pArg, "FIXEDHOFF" ) ) {
                 bVer1DefaultMode &= ~REQ_MODE_BASIC; /* only taut. */
            } else
            if ( !stricmp( pArg, "RECMETOFF" ) ) {
                 bReconnectCoord = 0;    
            } else
            if ( !stricmp( pArg, "SPXYZOFF" ) ) {
                 bAddPhosphineStereo = 0;
            } else
            if ( !stricmp( pArg, "SASXYZOFF" ) ) {
                 bAddArsineStereo = 0;
            } else
            if ( !stricmp( pArg, "FixSp3bugOFF" ) ) {
                 bFixSp3bug = 0;
            } else
            if ( !stricmp( pArg, "FBOFF" ) ) {
                 bFixSp3bug = 0;
            } else
            if ( !stricmp( pArg, "FB2OFF" ) ) {
                 bFixFB2 = 0;
            } else
            if ( !stricmp( pArg, "FNUDOFF" ) ) {
                 ip->bFixNonUniformDraw = 0;
            } else
            #endif
                
#endif /* INCLUDE_V102_FEATURES */


            #ifndef ONLY_STDINCHI_STDKEY
            if ( !stricmp( pArg, "RECMET" ) ) {    /* do reconnect metals */
                bReconnectCoord = 1;
            } else
            #endif

            if ( !stricmp( pArg, "AUXNONE" ) ) {    /* no aux. info */
                bINChIOutputOptions |= INCHI_OUT_NO_AUX_INFO;  /* no aux info */
                bINChIOutputOptions &= ~INCHI_OUT_SHORT_AUX_INFO;
            } else
            if ( !stricmp( pArg, "AUXFULL" ) || !stricmp( pArg, "AUXMAX" ) ) {     /* full aux info */
                bINChIOutputOptions &= ~(INCHI_OUT_NO_AUX_INFO | INCHI_OUT_SHORT_AUX_INFO); /* include short aux info */
            } else
            if ( !stricmp( pArg, "AUXMIN" ) ) {     /* minimal aux info */
                bINChIOutputOptions &= ~INCHI_OUT_NO_AUX_INFO; /* include short aux info */
                bINChIOutputOptions |= INCHI_OUT_SHORT_AUX_INFO;
            } else
            if ( !stricmp( pArg, "DONOTADDH" ) ) {
                bDoNotAddH = 1;
            } else
            if ( !stricmp( pArg, "D" ) ) {          /* display the structures */
                bDisplay |= 1;
            } else
#if( READ_INCHI_STRING == 1 )
            if ( !stricmp( pArg, "DDSRC" ) ) {
                bDisplayIfRestoreWarnings = 1;  /* InChI->Structure debugging: Display Debug Structure Restore Components */
            }
            else
#endif
            if ( !stricmp( pArg, "NOLABELS" ) ) {
                 bNoStructLabels = 1;
            } else
            if ( !stricmp( pArg, "WarnOnEmptyStructure" ) ) {
                 ip->bAllowEmptyStructure = 1;
            } else
            if ( !stricmp( pArg, "NoVarH" ) ) {
                 bTgFlagVariableProtons = 0;
            } else
            #ifndef ONLY_STDINCHI_STDKEY
            if ( !stricmp( pArg, "NoADP" ) ) {
                 bTgFlagHardAddRenProtons = 0;
            } else
            if ( !stricmp( pArg, "COMPRESS" ) ) {
                ip->bAbcNumbers = 1;
                ip->bCtPredecessors = 1;             /* compressed output */
            } else
            #endif
            if ( !stricmp( pArg, "FULL" ) ) {
                bVer1DefaultMode       = VER100_DEFAULT_MODE;
                nMode                = 0;
                bReconnectCoord      = 1;            /* full output */
                bINChIOutputOptions   = ((EMBED_REC_METALS_INCHI==1)? INCHI_OUT_EMBED_REC   : 0) | INCHI_OUT_SHORT_AUX_INFO;
                ip->bCtPredecessors  = 0;
                ip->bAbcNumbers      = 0;
                bXml                 |=  INCHI_OUT_PLAIN_TEXT | INCHI_OUT_PLAIN_TEXT_COMMENTS;
                bXml                 &= ~(INCHI_OUT_XML | INCHI_OUT_XML_TEXT_COMMENTS);
            } else
            if ( !stricmp( pArg, "MIN" ) ) {
                bVer1DefaultMode     = VER100_DEFAULT_MODE;
                nMode                = 0;
                bReconnectCoord      = 1;            /* minimal output */
                bINChIOutputOptions   = ((EMBED_REC_METALS_INCHI==1)? INCHI_OUT_EMBED_REC   : 0) | INCHI_OUT_NO_AUX_INFO;            /* minimal compressed output */
                ip->bCtPredecessors  = 1;
                ip->bAbcNumbers      = 1;
                bXml                |= INCHI_OUT_PLAIN_TEXT | INCHI_OUT_PLAIN_TEXT_COMMENTS;
                bXml                &= ~(INCHI_OUT_XML | INCHI_OUT_XML_TEXT_COMMENTS);
            } else
#if( READ_INCHI_STRING == 1 )
            #ifndef ONLY_STDINCHI_STDKEY
            if ( !stricmp( pArg, "InChI2InChI" )  ) {
                 /* Read InChI Identifiers and output InChI Identifiers */
                ip->nInputType = INPUT_INCHI;
                ip->bReadInChIOptions |= READ_INCHI_OUTPUT_INCHI;
                ip->bReadInChIOptions &= ~READ_INCHI_TO_STRUCTURE;
            } else
            if ( !stricmp( pArg, "SplitInChI" )  ) {
                 /* Split InChI Identifiers into components */
                ip->bReadInChIOptions |= READ_INCHI_SPLIT_OUTPUT;
            } else
            #endif /* ifndef ONLY_STDINCHI_STDKEY */
            if ( !stricmp( pArg, "InChI2Struct" )  ) {
                 /* Split InChI Identifiers into components */
                ip->bReadInChIOptions |= READ_INCHI_TO_STRUCTURE;
                ip->bReadInChIOptions &= ~READ_INCHI_OUTPUT_INCHI;
                ip->nInputType = INPUT_INCHI;
            } else
            if ( !stricmp( pArg, "KeepBalanceP" )  ) {
                 /* When spliting InChI Identifiers into components: */
                 /* If MobileH output then add p to each component;  */
                 /* Otherwise add one more component containing balance */
                 /* of protons and exchangeable isotopic H */
                ip->bReadInChIOptions |= READ_INCHI_KEEP_BALANCE_P;
            } else			
#endif
            if ( /* INPUT_NONE == ip->nInputType &&*/
                !memicmp( pArg, "SDF:", 4 )  ) {
                 /* SDfile label */
                k = 0;
                mystrncpy( ip->szSdfDataHeader, pArg+4, MAX_SDF_HEADER+1 );
                LtrimRtrim( ip->szSdfDataHeader, &k );
                if ( k ) {
                    ip->pSdfLabel  = ip->szSdfDataHeader;
                    ip->pSdfValue  = szSdfDataValue;
                    if ( INPUT_NONE == ip->nInputType ) {
                        ip->nInputType = INPUT_SDFILE;
                    }
                } else {
                    ip->pSdfLabel  = NULL;
                    ip->pSdfValue  = NULL;
                    if ( INPUT_NONE == ip->nInputType ) {
                        ip->nInputType = INPUT_MOLFILE;
                    }
                }

/*^^^ post-1.02b */
            } else 
            if ( !stricmp( pArg, "INPAUX" ))
            {
                if (INPUT_NONE == ip->nInputType) 
                    ip->nInputType = INPUT_INCHI_PLAIN;
            } else


#if( ADD_CMLPP == 1 )
            if ( INPUT_NONE == ip->nInputType && !stricmp( pArg, "CML" )  ) {
                 /* CMLfile label */
                ip->nInputType = INPUT_CMLFILE;
            } else
#endif
            /*============= Default options =============*/
            #ifndef ONLY_STDINCHI_STDKEY
            if ( !stricmp( pArg, "RECMET-" ) ) {     /* do not reconnect metals (default) */
                bReconnectCoord = 0;
            } else
            #endif
            
            if ( !stricmp( pArg, "OUTPUTSDF" ) ) {  /* output SDfile */
                bOutputMolfileOnly = 1;
            } else
            if ( !stricmp( pArg, "SdfAtomsDT" ) ) {  /* output isotopes H as D and T in SDfile */
                bOutputMolfileDT = 1;
            } else                                    
            if ( !stricmp( pArg, "SdfSplit" ) ) {  /* Split single Molfiles into disconnected components */
                bOutputMolfileSplit = 1;
            } else
            if ( !stricmp( pArg, "STDIO" ) ) {
                bNameSuffix = 0;
            } else
            if ( !stricmp( pArg, "DCR" ) ) {
                bDisplayCompositeResults = 1;
            } else
#if( FIX_ADJ_RAD == 1 )
            if ( !stricmp( pArg, "FixRad" ) ) {
                bFixAdjacentRad = 1;
            } else
#endif
            /*============= Additional options ==========*/
            /* Tautomer perception off */
            if ( !stricmp( pArg, "EXACT" ) ) {
                bVer1DefaultMode |= REQ_MODE_BASIC;
            } else
            if ( !stricmp( pArg, "MOLFILENUMBER" ) ) {
                ip->bGetMolfileNumber |= 1;
            } else
            if ( !stricmp( pArg, "OutputPLAIN" ) ) {
                bXml |=  INCHI_OUT_PLAIN_TEXT;
                bXml &= ~INCHI_OUT_XML;
            } else
#if( !defined(INCHI_LIBRARY) && !defined(INCHI_LIB) )
            if ( !stricmp( pArg, "Tabbed" ) ) {
                bXml |=  INCHI_OUT_TABBED_OUTPUT;
            } else
#endif
            if ( !stricmp( pArg, "OutputXML" ) ) {
                bXml |= INCHI_OUT_XML;
                bXml &= ~INCHI_OUT_PLAIN_TEXT;
            } else
            if ( !stricmp( pArg, "OutputANNPLAIN" ) ) {
                bXml |=   INCHI_OUT_PLAIN_TEXT_COMMENTS;
                bXml &=  ~INCHI_OUT_XML_TEXT_COMMENTS;
                bXml |=   INCHI_OUT_WINCHI_WINDOW; /* debug */
            } else
            if ( !stricmp( pArg, "OutputANNXML" ) ) {
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
            if ( !memicmp( pArg, "START:", 6 ) ) {
                ip->first_struct_number = strtol(pArg+6, NULL, 10);
            } else
            if ( !memicmp( pArg, "END:", 4 ) ) {
                ip->last_struct_number = strtol(pArg+4, NULL, 10);
            } else /*  RSB: */
            if ( !memicmp( pArg, "RSB:", 4 )) {
                mdbr = (int)strtol(pArg+4, NULL, 10);
            } else
            if ( !stricmp( pArg, "EQU" ) ) {
                bCompareComponents = CMP_COMPONENTS;
            } else
            if ( !stricmp( pArg, "EQUNONISO" ) ) {
                bCompareComponents = CMP_COMPONENTS | CMP_COMPONENTS_NONISO;
            } else
            if ( !memicmp( pArg, "OP:", 3 ) ) {
                bOutputPath = 1;
                strncpy(szOutputPath, pArg+3, sizeof(szOutputPath)-1);
            } else
            /*============== Char+Value options ==============*/
            if ( !memicmp( pArg, "W", 1 ) && (t = strtod( pArg+1, (char**)&q ), q > pArg+1) ) {
                if ( errno == ERANGE || t < 0.0 || t*1000.0 > (double)ULONG_MAX)  {
                    ul = 0;
                } else {
                    ul = (unsigned long)(t*1000.0);  /* max. time per structure */
                }
                ip->msec_MaxTime = ul;
            } else
            if ( !memicmp( pArg, "F", 1 ) && (c =  (int)strtol( pArg+1, (char**)&q, 10 ), q > pArg+1) ) {
                nFontSize = -c;                      /* struct. display font size */
            } else

            /* in-house options */
#ifndef ONLY_STDINCHI_STDKEY    

#if( UNDERIVATIZE == 1 )
            if ( !stricmp( pArg, "DoDRV" ) ) {
                ip->bUnderivatize = 1;
            } else
#endif
#if( RING2CHAIN == 1 )
            if ( !stricmp( pArg, "DoR2C" ) ) {
                ip->bRing2Chain = 1;
            } else
#endif

#if ( RING2CHAIN == 1 || UNDERIVATIZE == 1 )
            if ( !stricmp( pArg, "DoneOnly" ) ) {
                ip->bIngnoreUnchanged = 1;
            } else
#endif

#if ( KETO_ENOL_TAUT == 1 )
            if ( !stricmp( pArg, "KET" ) ) {
                bKetoEnolTaut = 1;
            } else
#endif

#if ( TAUT_15_NON_RING == 1 )
            if ( !stricmp( pArg, "15T" ) ) {
                b15TautNonRing = 1;
            } else
#endif

#endif /* ONLY_STDINCHI_STDKEY     */


            /* Display unrecognized option */
            {
                bRecognizedOption = 0;
#ifdef ONLY_STDINCHI_STDKEY
                inchi_ios_eprint(log_file, "Terminating: there is no option \"%s\" in standard InChI software library\n", pArg);
                return -1;
#else
#ifndef INCHI_LIB
                inchi_ios_eprint(log_file, "Unrecognized option: \"%s\".\n", pArg);
#endif
#endif

            }
            bVer1Options |= bRecognizedOption;
        } else
        if ( ip->num_paths < MAX_NUM_PATHS ) {
            char *sz;
            if ( argv[i] && argv[i][0] ) {
                if ( sz = (char*)inchi_malloc( (strlen(argv[i]) + 1)*sizeof(sz[0])) ) {
                    strcpy( sz, argv[i] );
                }
#ifdef CML_DEBUG
                printf ("1 path %d argv %d %s\n", ip -> num_paths, i, argv [i]);
#endif
                   ip->path[ip->num_paths++] = sz;
            }
        }
    }


    if (ip->bCalcInChIKey!=0)
    /* Suppress InChIKey calculation if: */
    /* compressed output */
    /* Inchi2Struct */
    /* Inchi2Inchi */
    {

        if ( (ip->bAbcNumbers ==1) && (ip->bCtPredecessors == 1) )
        {
#ifndef INCHI_LIB
            inchi_ios_eprint(log_file, "Terminating: generation of InChIKey is not available with 'Compress' option\n");
            return -1;
#endif
            ip->bCalcInChIKey = 0;
        }
        if ( ip->nInputType == INPUT_INCHI )
        {
#ifndef INCHI_LIB
            inchi_ios_eprint(log_file, "Terminating: generation of InChIKey is not available with 'Inchi2Struct' option\n");
            return -1;
#endif
            ip->bCalcInChIKey = 0;
        }
        else
        if ( bOutputMolfileOnly == 1 )
        {
#ifndef INCHI_LIB
            inchi_ios_eprint(log_file, "Terminating: generation of InChIKey is not available with 'OutputSDF' option\n");
            return -1;
#endif
            ip->bCalcInChIKey = 0;
        }

    }


    if ( bNameSuffix || bOutputPath ) {
        const char *p = NULL;
        char       *r = NULL;
        char       *sz;
        int  len;
        const char szNUL[] = "NUL"; /* fix for AMD processor: use const char[] instead of just "NUL" constant 2008-11-5 DT */

        /*  find the 1st path */
        for ( i = 0; i < MAX_NUM_PATHS; i ++ ) {
            if ( !p && ip->path[i] && ip->path[i][0] ) {
                p = ip->path[i];
                break;
            }
        }
        /* fix output path */
        if ( bOutputPath && szOutputPath[0] && p ) {
            /* remove last slash */
            len = strlen(szOutputPath);
            if ( len > 0 && szOutputPath[len-1] != INCHI_PATH_DELIM ) {
                szOutputPath[len++] = INCHI_PATH_DELIM;
                szOutputPath[len]   = '\0';
            }
            if ( len > 0 && (r = (char *)strrchr( p, INCHI_PATH_DELIM ) ) && r[1] ) {
                strcat( szOutputPath, r+1 );
                p = szOutputPath;
            }
        }        /*  add missing paths */
        for ( i = 0; p && i < MAX_NUM_PATHS; i ++ ) { /* fix for AMD processor: changed order 2008-11-5 DT */
            if ( !ip->path[i] || !ip->path[i][0] ) {
                len = strlen( p ) + strlen(szNameSuffix) + strlen( ext[i] );
                if ( sz = (char*)inchi_malloc( (len+1)*sizeof(sz[0]) ) ) {
                    strcpy( sz, p );
                    strcat( sz, szNameSuffix );
                    strcat( sz, ext[i] );
                    ip->num_paths++;
                    ip->path[i] =sz;
                }
             } else
            if ( !stricmp( ip->path[i], szNUL ) ) {
                inchi_free( (char *)ip->path[i] ); /* cast deliberately const qualifier */
                ip->path[i] = NULL;
            }
        }
    }

#endif  /* } NOT ACD_LABS_VERSION */


#if( READ_INCHI_STRING == 1 )
    if ( INPUT_INCHI == ip->nInputType ) {
        bCompareComponents                 = 0;
        /*bDisplayCompositeResults           = 0;*/
#if( I2S_MODIFY_OUTPUT == 1 )
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
#if( I2S_MODIFY_OUTPUT == 1 )
              && !bOutputMolfileOnly
              && !(bXml & (INCHI_OUT_XML | INCHI_OUT_XML_TEXT_COMMENTS | INCHI_OUT_PLAIN_TEXT | INCHI_OUT_XML_TEXT_COMMENTS))
#endif
           ) {
            bINChIOutputOptions |= INCHI_OUT_PLAIN_TEXT;
        }
    }
#endif


    if ( bVer1Options ) {
        nMode |= bVer1DefaultMode;
    } else
    if ( bReleaseVersion ) {
        nMode |= nReleaseMode;
    }

#if( defined(INCHI_ANSI_ONLY) || defined(INCHI_LIB) )
    if ( bCompareComponents && !(bDisplay & 1) ) {
        bCompareComponents = 0;
    }
#endif
    /*  Save original options */
    /* nOrigMode = nMode; */
#ifndef INCHI_ANSI_ONLY
    ip->dp.sdp.nFontSize         = nFontSize;
    ip->dp.sdp.ulDisplTime       = *ulDisplTime;
    ip->bDisplay                 = bDisplay;
#ifdef INCHI_LIB
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
#if( MIN_SB_RING_SIZE > 0 )
    if ( mdbr ) {
        nMinDbRinSize = mdbr;
    }
    nMode |= (nMinDbRinSize << REQ_MODE_MIN_SB_RING_SHFT) & REQ_MODE_MIN_SB_RING_MASK;
#endif
    /*  input file */
    if ( ip->nInputType == INPUT_NONE && ip->num_paths > 0 ) {
        ip->nInputType = INPUT_MOLFILE; /*  default */
#if( ADD_CMLPP == 1 )
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
#if( SDF_OUTPUT_DT == 1 )
        ip->bINChIOutputOptions |= bOutputMolfileDT?    INCHI_OUT_SDFILE_ATOMS_DT : 0;
        ip->bINChIOutputOptions |= bOutputMolfileSplit? INCHI_OUT_SDFILE_SPLIT : 0;
#endif
    }
    if ( bXml & INCHI_OUT_XML ) {
        bXml &= ~(INCHI_OUT_PLAIN_TEXT | INCHI_OUT_XML_TEXT_COMMENTS | INCHI_OUT_TABBED_OUTPUT);
    }
#ifdef INCHI_LIB
    if ( !(bDisplay & 1) ) {
        bXml &= ~(INCHI_OUT_PLAIN_TEXT_COMMENTS | INCHI_OUT_XML_TEXT_COMMENTS); /* do not ouput comments in wINChI text file results */
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
#if( FIX_ADJ_RAD == 1 )
    ip->bTautFlags  |= bFixAdjacentRad?         TG_FLAG_FIX_ADJ_RADICALS : 0;
#endif
#if( ADD_PHOSPHINE_STEREO == 1 )
    ip->bTautFlags  |= bAddPhosphineStereo?     TG_FLAG_PHOSPHINE_STEREO : 0;
#endif
#if( ADD_ARSINE_STEREO == 1 )
    ip->bTautFlags  |= bAddArsineStereo?        TG_FLAG_ARSINE_STEREO : 0;
#endif
    ip->bTautFlags  |= bFixSp3bug?              TG_FLAG_FIX_SP3_BUG   : 0;


    /*^^^ */
    if (bFixFB2)
    {
#if( FIX_ISO_FIXEDH_BUG == 1 )
        ip->bTautFlags  |= TG_FLAG_FIX_ISO_FIXEDH_BUG; /* accomodate FIX_ISO_FIXEDH_BUG */
#endif
#if( FIX_TERM_H_CHRG_BUG == 1 )
        ip->bTautFlags  |= TG_FLAG_FIX_TERM_H_CHRG_BUG; /* accomodate FIX_TERM_H_CHRG_BUG */
#endif
#if( FIX_TRANSPOSITION_CHARGE_BUG == 1 )
        ip->bINChIOutputOptions |= INCHI_OUT_FIX_TRANSPOSITION_CHARGE_BUG;
#endif
    }
#ifdef ONLY_STDINCHI_STDKEY
    ip->bINChIOutputOptions |= INCHI_OUT_STANDARD_FLAG;
#endif	

    if ( !ip->nInputType ) {
        ip->nInputType = INPUT_MOLFILE;
    }



    return 0;
}
/*******************************************************************/
int PrintInputParms( INCHI_IOSTREAM *log_file, INPUT_PARMS *ip )
{
    INCHI_MODE nMode = ip->nMode;
    int i, k;
     
#ifdef INCHI_LIB
    int bInChI2Struct = 0; /* as of December 2008, winchi-1 does not convert InChI to structure */
#else
    int bInChI2Struct = (ip->bReadInChIOptions & READ_INCHI_TO_STRUCTURE) && ip->nInputType == INPUT_INCHI;
#endif
            


    i = 0;
    
    /*  generation/conversion indicator */
    
    #ifdef ONLY_STDINCHI_STDKEY
    if ( !(ip->bINChIOutputOptions & INCHI_OUT_SDFILE_ONLY) && !bInChI2Struct )
        inchi_ios_eprint( log_file, "Generating standard InChI\n" );
    else
#if( !defined(INCHI_LIBRARY) && !defined(INCHI_LIB) && !defined(INCHI_MAIN) )
        if ( bInChI2Struct ) /* effective only in command line program cInChI or stdInChI */
            inchi_ios_eprint( log_file, "Converting InChI(s) to structure(s) in %s\n",
            (ip->bINChIOutputOptions & INCHI_OUT_SDFILE_ONLY)? "MOL format" : "aux. info format" );
#endif
    #else
    inchi_ios_eprint( log_file, "Options: \n" );
    #endif

    /*  output options: line 1 */

    /* SDfile output */
#if 0
    if ( ip->bINChIOutputOptions & INCHI_OUT_SDFILE_ONLY ) {
        inchi_ios_eprint( log_file, "Output SDfile only without stereochemical information and atom coordinates%s\n",
            (ip->bINChIOutputOptions & INCHI_OUT_SDFILE_ATOMS_DT)? "\n(write H isotopes as D, T)":"" );
        ;/*return 0;*/
    }
#endif
    /* tautomerism */
    if( (nMode & (REQ_MODE_BASIC | REQ_MODE_TAUT )) == (REQ_MODE_BASIC | REQ_MODE_TAUT) ) {
        inchi_ios_eprint( log_file, "Mobile H Perception OFF (include FixedH layer)" );
    } else
    if( (nMode & (REQ_MODE_BASIC | REQ_MODE_TAUT )) == (REQ_MODE_TAUT) ) 
    {
#ifndef ONLY_STDINCHI_STDKEY
        inchi_ios_eprint( log_file, "Mobile H Perception ON (omit FixedH layer)" );
#endif
    } else
    if( (nMode & (REQ_MODE_BASIC | REQ_MODE_TAUT )) == (REQ_MODE_BASIC) ) {
        inchi_ios_eprint( log_file, "Mobile H ignored" );
    } else {
        inchi_ios_eprint( log_file, "Undefined Mobile H mode" );
    }
    if ( (ip->bTautFlags & TG_FLAG_VARIABLE_PROTONS) ) { 
         if ( !(ip->bTautFlags & TG_FLAG_HARD_ADD_REM_PROTONS) ) {
            inchi_ios_eprint( log_file, ", Disabled Aggressive (De)protonation" );
         }
    }
    /*inchi_ios_eprint( log_file, "\n");*/
#if( FIX_ADJ_RAD == 1 )
    if ( ip->bTautFlags & TG_FLAG_FIX_ADJ_RADICALS ) {
        inchi_ios_eprint( log_file, "Fix Adjacent Radicals\n" );
    }
#endif
    i = 0;
    /* isotopic */
    if ( nMode & REQ_MODE_ISO ) {
#ifndef ONLY_STDINCHI_STDKEY
        inchi_ios_eprint( log_file, "Isotopic ON\n");
#endif
    } else
    if ( nMode & REQ_MODE_NON_ISO ) {
        inchi_ios_eprint( log_file, "Isotopic OFF\n");
    }
    i ++;
    /*  stereo */
    if ( nMode & REQ_MODE_STEREO ) {
#ifndef ONLY_STDINCHI_STDKEY
        inchi_ios_eprint( log_file,  "%s%s%s%sStereo ON",
                     ( nMode & REQ_MODE_NOEQ_STEREO )?     "Slow ":"",
                     ( nMode & REQ_MODE_REDNDNT_STEREO )?  "Redund. ":"",
                     ( nMode & REQ_MODE_NO_ALT_SBONDS)?    "No AltBond ":"",

                     ( nMode & REQ_MODE_RACEMIC_STEREO)?   "Racemic " :
                     ( nMode &  REQ_MODE_RELATIVE_STEREO)? "Relative " :
                     ( nMode &  REQ_MODE_CHIR_FLG_STEREO)? "Chiral Flag " : "Absolute " );
        if ( 0 == (nMode & (REQ_MODE_SB_IGN_ALL_UU | REQ_MODE_SC_IGN_ALL_UU)) ) {
            inchi_ios_eprint( log_file, "\nInclude undefined/unknown stereogenic centers and bonds");
        } else
        if ( REQ_MODE_SC_IGN_ALL_UU == (nMode & (REQ_MODE_SB_IGN_ALL_UU | REQ_MODE_SC_IGN_ALL_UU)) ) {
            inchi_ios_eprint( log_file, "\nOmit undefined/unknown stereogenic centers");
        } else
        if ( REQ_MODE_SB_IGN_ALL_UU == (nMode & (REQ_MODE_SB_IGN_ALL_UU | REQ_MODE_SC_IGN_ALL_UU)) ) {
            inchi_ios_eprint( log_file, "\nOmit undefined/unknown stereogenic bonds");
        } else {
        /*case REQ_MODE_SB_IGN_ALL_UU | REQ_MODE_SC_IGN_ALL_UU*/
            inchi_ios_eprint( log_file, "\nOmit undefined/unknown stereogenic centers and bonds");
        }
#endif /* ONLY_STDINCHI_STDKEY */

#if( defined(MIN_SB_RING_SIZE) && MIN_SB_RING_SIZE > 0 )
        k = (ip->nMode & REQ_MODE_MIN_SB_RING_MASK) >> REQ_MODE_MIN_SB_RING_SHFT;
        if ( bRELEASE_VERSION != 1 || k != MIN_SB_RING_SIZE ) {
            if ( k >= 3 ) {
                inchi_ios_eprint( log_file, "\nMin. stereobond ring size: %d\n", k );
            } else {
                inchi_ios_eprint( log_file, "\nMin. stereobond ring size: NONE\n" );
            }
            i = 0;
        }
#endif
        if ( TG_FLAG_POINTED_EDGE_STEREO & ip->bTautFlags ) 
        {
#ifndef ONLY_STDINCHI_STDKEY
            inchi_ios_eprint( log_file, "Only narrow end of wedge points to stereocenter\n");
#endif
            i = 0;
        }
        else
        {
            inchi_ios_eprint( log_file, "Both ends of wedge point to stereocenters\n");
        }

#if( ADD_PHOSPHINE_STEREO == 1 )
        if ( TG_FLAG_PHOSPHINE_STEREO & ip->bTautFlags ) 
        {
#ifndef ONLY_STDINCHI_STDKEY
            inchi_ios_eprint( log_file, "Include phosphine stereochemistry\n");
#endif
            i = 0;
        }
#endif
#if( ADD_ARSINE_STEREO == 1 )
        if ( TG_FLAG_ARSINE_STEREO & ip->bTautFlags ) 
        {
#ifndef ONLY_STDINCHI_STDKEY
            inchi_ios_eprint( log_file, "Include arsine stereochemistry\n");
#endif
            i = 0;
        }
#endif
        if ( TG_FLAG_FIX_SP3_BUG & ip->bTautFlags ) {
#ifndef ONLY_STDINCHI_STDKEY
            inchi_ios_eprint( log_file, "Fix bug leading to missing or undefined sp3 parity\n");
#endif
            i = 0;
        }
    } else {
        inchi_ios_eprint( log_file, "Stereo OFF\n");
    }
    if ( ip->bDoNotAddH ) 
    {
        inchi_ios_eprint( log_file, "Do not add H\n");
    }


    /*^^^ FB2, FNUD*/
#ifndef ONLY_STDINCHI_STDKEY
    if ( (TG_FLAG_FIX_ISO_FIXEDH_BUG & ip->bTautFlags) )
    {
        inchi_ios_eprint( log_file, "Fix bugs found after v.1.02b release\n");
    }
    if ( ip->bFixNonUniformDraw ) 
    {
        inchi_ios_eprint( log_file, "Fix non-uniform drawing issues\n");
    }
#endif

    /* metals disconnection */
    if ( ip->bTautFlags & TG_FLAG_DISCONNECT_COORD ) {
        if ( ip->bTautFlags & TG_FLAG_RECONNECT_COORD ) {
            inchi_ios_eprint( log_file, "Include bonds to metals\n");
        }
    } else {
        inchi_ios_eprint( log_file, "Do not disconnect metals\n");
    }




    /*  other options: line 2 */
#if( bRELEASE_VERSION == 1 )
    if ( ip->bCtPredecessors || ip->bAbcNumbers ) {
        if ( ip->bCtPredecessors && ip->bAbcNumbers ) {
            inchi_ios_eprint( log_file, "Representation: Compressed");
            i ++;
        } else {
            inchi_ios_eprint( log_file, "Connection table: %s, %s\n",
                ip->bCtPredecessors? "Predecessor_numbers(closures)":"Canon_numbers(branching, ring closures)",
                ip->bAbcNumbers?     "Shorter alternative":"Numerical");
            i = 0;
        }
    }
#else
    if ( (bRELEASE_VERSION != 1) || ip->bCtPredecessors || ip->bAbcNumbers ) {
        inchi_ios_eprint( log_file, "Connection table: %s, %s\n",
            ip->bCtPredecessors? "Predecessor_numbers(closures)":"Canon_numbers(branching, ring closures)",
            ip->bAbcNumbers?     "Shorter alternative":"Numerical");
        i = 0;
    } else {
        inchi_ios_eprint( log_file, "Representation: Numerical");
        i ++;
    }
#endif

    if ( !(ip->bINChIOutputOptions & INCHI_OUT_SDFILE_ONLY) ) 
    {

        if( ip->bINChIOutputOptions & INCHI_OUT_NO_AUX_INFO ) {
            inchi_ios_eprint( log_file, "Aux. info suppressed\n");
            i ++;
        } else
        if ( ip->bINChIOutputOptions & INCHI_OUT_SHORT_AUX_INFO ) {
            inchi_ios_eprint( log_file, "Minimal Aux. info\n");
            i ++;
        } else {
            inchi_ios_eprint( log_file, "Full Aux. info\n");
            i ++;
        }
    }

    if ( ip->bCalcInChIKey) 
    {
#ifdef ONLY_STDINCHI_STDKEY
        inchi_ios_eprint( log_file, "Generating standard InChIKey\n");
#else
        inchi_ios_eprint( log_file, "Generating InChIKey\n");
#endif
        i ++;
    }
    /*^^^ */



    if ( ip->bAllowEmptyStructure ) {
        inchi_ios_eprint( log_file, "Issue warning on empty structure\n" );
    }

    if ( ip->szSdfDataHeader[0] && ip->nInputType != INPUT_SDFILE ) {
        inchi_ios_eprint( log_file, "SDfile data header: \"%s\"\n", ip->szSdfDataHeader);
    }
    /* input format */
    if ( ip->nInputType ) {
        inchi_ios_eprint( log_file, "Input format: %s",
            ip->nInputType == INPUT_MOLFILE?     "MOLfile"       :
            ip->nInputType == INPUT_SDFILE?      "SDfile"        :
            ip->nInputType == INPUT_CMLFILE?     "CMLfile"       :
#if( READ_INCHI_STRING == 1 )
            ip->nInputType == INPUT_INCHI?       "InChI (plain identifier)" :
#endif
#if( SPECIAL_BUILD == 1 )
            ip->nInputType == INPUT_INCHI_XML?   "MoChI+AuxInfo (xml)"   :
            ip->nInputType == INPUT_INCHI_PLAIN? "MoChI+AuxInfo (plain)" : "Unknown" );
#else
            ip->nInputType == INPUT_INCHI_XML?   "InChI AuxInfo (xml)"   :
            ip->nInputType == INPUT_INCHI_PLAIN? "InChI AuxInfo (plain)" : "Unknown" );
#endif



        if ( (ip->nInputType == INPUT_MOLFILE || ip->nInputType == INPUT_SDFILE) &&
             ip->bGetMolfileNumber ) {
            inchi_ios_eprint( log_file, "  (attempting to read Molfile number)" );
        }
        inchi_ios_eprint( log_file, "\n");
    }

    /*  output format */
    inchi_ios_eprint( log_file, "Output format: %s%s\n", 
        (ip->bINChIOutputOptions & INCHI_OUT_PLAIN_TEXT)?  "Plain text" :
        (ip->bINChIOutputOptions & INCHI_OUT_XML)?         "XML":
        ((ip->bINChIOutputOptions & INCHI_OUT_SDFILE_ONLY) && bInChI2Struct )? "SDfile only (without stereochemical info and atom coordinates)" :
        ((ip->bINChIOutputOptions & INCHI_OUT_SDFILE_ONLY) && !bInChI2Struct)? "SDfile only" : "Unknown",

        ((ip->bINChIOutputOptions & INCHI_OUT_PLAIN_TEXT) &&
        (ip->bINChIOutputOptions & INCHI_OUT_TABBED_OUTPUT))? ", tabbed":"");
    i = 0;
    /*  other options: line 4 */
    if ( ip->msec_MaxTime ) {
        unsigned long seconds = ip->msec_MaxTime/1000;
        unsigned long milliseconds = (ip->msec_MaxTime%1000);
        inchi_ios_eprint( log_file, "Timeout per structure: %lu.%03lu sec", seconds, milliseconds);
        i ++;
    } else {
        inchi_ios_eprint( log_file, "No timeout");
        i ++;
    }
    inchi_ios_eprint( log_file, "%sUp to %d atoms per structure\n", i?"; ":"", MAX_ATOMS);
    i = 0;
    if ( ip->first_struct_number > 1 ) {
        inchi_ios_eprint( log_file, "Skipping %ld structure%s\n", ip->first_struct_number-1, ip->first_struct_number==2? "":"s" );
    }
    if ( ip->last_struct_number > 0 ) {
        inchi_ios_eprint( log_file, "Terminate after structure #%ld\n", ip->last_struct_number );
    }
    if ( ip->bSaveWarningStructsAsProblem && ip->path[3] && ip->path[3][0] ) {
        inchi_ios_eprint( log_file, "Saving warning structures into the problem file\n");
    }
    if ( ip->bSaveAllGoodStructsAsProblem && ip->path[3] && ip->path[3][0] ) {
        inchi_ios_eprint( log_file, "Saving only all good structures into the problem file\n");
    }
    /*  Report debug modes */
    i = 0;
#if( bRELEASE_VERSION != 1 )
    inchi_ios_eprint( log_file, "Release version = NO");
    i ++;
#endif

#if( FIND_RING_SYSTEMS != 1 )
    inchi_ios_eprint( log_file, "%s5-, 6-, 7-memb. ring taut. ignored", i?"; ":"");
    i ++;
#endif


#if( TRACE_MEMORY_LEAKS == 1 && defined(_DEBUG) )
    inchi_ios_eprint( log_file, "%sTracing memory leaks (SLOW)", i?"; ":"");
    i ++;
#endif

    if ( i ) {
        inchi_ios_eprint( log_file, "\n" );
    }



#if( bRELEASE_VERSION != 1 )

#if( FIND_RING_SYSTEMS == 1 )
    inchi_ios_eprint( log_file, "Find ring systems=Y\nTautomers:" );
    inchi_ios_eprint( log_file, " 4-pyridinol=%s", TAUT_4PYRIDINOL_RINGS==1? "Y":"N");
    inchi_ios_eprint( log_file, " pyrazole=%s", TAUT_PYRAZOLE_RINGS==1? "Y":"N");
    inchi_ios_eprint( log_file, " tropolone=%s", TAUT_TROPOLONE_7==1? "Y":"N");
    inchi_ios_eprint( log_file, " tropolone-5=%s", TAUT_TROPOLONE_5==1? "Y":"N");
    inchi_ios_eprint( log_file, "\n" );
    inchi_ios_eprint( log_file, "Only chain attachments to tautomeric rings=%s\n", TAUT_RINGS_ATTACH_CHAIN==1? "Y":"N");
#endif
    if ( ip->bGetSdfileId ) {
        inchi_ios_eprint( log_file, "Extracting SDfile IDs\n");
    }
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
#if( TEST_RENUMB_ATOMS == 1 )
    inchi_ios_eprint( log_file, "\nDbg: TEST_RENUMB_ATOMS=%d; TEST_RENUMB_NEIGH=%d; TEST_RENUMB_SWITCH=%d\n",
                                  TEST_RENUMB_ATOMS,    TEST_RENUMB_NEIGH,    TEST_RENUMB_SWITCH );
    inchi_ios_eprint( log_file, "     TEST_RENUMB_ATOMS_SAVE_LONGEST=%d\n",
                                  TEST_RENUMB_ATOMS_SAVE_LONGEST);
#endif

#endif

    return 0;
}
/************************************************************************************/
void HelpCommandLineParms( INCHI_IOSTREAM *f )
{
    if ( !f )
        return;

#if ( bRELEASE_VERSION == 1 )

    /*^^^ */
     inchi_ios_print_nodisplay( f, 
#ifdef INCHI_MAIN

#ifdef ONLY_STDINCHI_STDKEY
         "%s ver %s%s.\n\nUsage:\nstdinchi_main inFile [outFile [logFile [problemFile]]] [%coption[ %coption...]]\n", 
#else
         "%s ver %s%s.\n\nUsage:\nInChI_MAIN inputFile [outputFile [logFile [problemFile]]] [%coption[ %coption...]]\n", 
#endif
         INCHI_NAME, INCHI_VERSION, SPECIAL_BUILD_STRING, 
         INCHI_OPTION_PREFX, INCHI_OPTION_PREFX); 

#else

#ifdef STDINCHI_AS_DEFAULT
         "%s ver %s%s.\n\nUsage:\nstdinchi-%s inFile [outFile [logFile [problemFile]]] [%coption[ %coption...]]\n", 
         INCHI_NAME, INCHI_VERSION, SPECIAL_BUILD_STRING, 
         INCHI_VERSION, INCHI_OPTION_PREFX, INCHI_OPTION_PREFX); 
#else
         "%s ver %s%s.\n\nUsage:\nc%s-%s inputFile [outputFile [logFile [problemFile]]] [%coption[ %coption...]]\n", 	     
         INCHI_NAME, INCHI_VERSION, SPECIAL_BUILD_STRING, 
         INCHI_NAME, INCHI_VERSION, INCHI_OPTION_PREFX, INCHI_OPTION_PREFX); 
#endif
#endif
    /*^^^ */    
    
     
    inchi_ios_print_nodisplay( f, "\nOptions:\n");
    

    inchi_ios_print_nodisplay( f, "\nInput\n");
    inchi_ios_print_nodisplay( f, "  STDIO       Use standard input/output streams\n");
    /*^^^ post-1.02bv */
    inchi_ios_print_nodisplay( f, "  InpAux      Input structures in %s default aux. info format\n              (for use with STDIO)\n", INCHI_NAME);
    inchi_ios_print_nodisplay( f, "  SDF:DataHeader Read from the input SDfile the ID under this DataHeader\n");
#if( ADD_CMLPP == 1 )
    inchi_ios_print_nodisplay( f, "  CML         Input in CML format (default for input file .CML extension)\n");
#endif


    inchi_ios_print_nodisplay( f, "Output\n");
    
    
    inchi_ios_print_nodisplay( f, "  AuxNone     Omit auxiliary information (default: Include)\n");
    /*inchi_ios_print_nodisplay( f, "  AuxMin      Output minimal auxiliary information\n");*/
    inchi_ios_print_nodisplay( f, "  NoLabels    Omit structure number, DataHeader and ID from %s output\n", INCHI_NAME);
    inchi_ios_print_nodisplay( f, "  Tabbed      Separate structure number, %s, and AuxInfo with tabs\n", INCHI_NAME);
#ifndef ONLY_STDINCHI_STDKEY
    inchi_ios_print_nodisplay( f, "  Compress    Compressed output\n");
#endif
    /*inchi_ios_print_nodisplay( f, "    FULL        Standard set of options for Full Verbose Output\n");*/
    /*inchi_ios_print_nodisplay( f, "    MIN         Standard set of options for Minimal Concise Output\n");*/
#if( defined(_WIN32) && defined(_MSC_VER) && !defined(INCHI_ANSI_ONLY) && !defined(INCHI_LIBRARY) )
    inchi_ios_print_nodisplay( f, "  D           Display the structures\n");
    inchi_ios_print_nodisplay( f, "  EQU         Display sets of identical components\n");
    inchi_ios_print_nodisplay( f, "  Fnumber     Set display Font size in number of points\n");
#endif
    /*inchi_ios_print_nodisplay( f, "    PLAIN       Plain text output (Default: XML format)\n");*/
    inchi_ios_print_nodisplay( f, "  OutputSDF   Convert %s created with default aux. info to SDfile\n", INCHI_NAME);
#if ( SDF_OUTPUT_DT == 1 )
    inchi_ios_print_nodisplay( f, "  SdfAtomsDT  Output Hydrogen Isotopes to SDfile as Atoms D and T\n");
#endif


    inchi_ios_print_nodisplay( f, "Structure perception\n");
#ifndef STDINCHI_AS_DEFAULT
#ifndef ONLY_STDINCHI_STDKEY
    inchi_ios_print_nodisplay( f, "  NEWPS       Narrow end of wedge points to stereocenter (default: both)\n");
#endif
#endif
#ifdef INCLUDE_V102_FEATURES
    inchi_ios_print_nodisplay( f, "  NEWPSOFF    Both ends of wedge point to stereocenters\n");
#endif
    inchi_ios_print_nodisplay( f, "  DoNotAddH   Don't add H according to usual valences: all H are explicit\n");
#ifndef ONLY_STDINCHI_STDKEY
    inchi_ios_print_nodisplay( f, "  NoADP       Disable Aggressive Deprotonation (for testing only)\n");
#endif
#ifdef ONLY_STDINCHI_STDKEY
    inchi_ios_print_nodisplay( f, "  SNon        Exclude stereo\n");
#endif

    inchi_ios_print_nodisplay( f, "Generation\n");
/*^^^ */
    inchi_ios_print_nodisplay( f, "  Wnumber     Set time-out per structure in seconds; W0 means unlimited\n");
    inchi_ios_print_nodisplay( f, "  WarnOnEmptyStructure Warn and produce empty %s for empty structure\n", INCHI_NAME);

#ifdef INCLUDE_V102_FEATURES
#ifdef ONLY_STDINCHI_STDKEY
    inchi_ios_print_nodisplay( f, "  Key         Generate standard InChIKey\n");
#else
    inchi_ios_print_nodisplay( f, "  Key         Generate InChIKey\n");
#endif
#endif

#ifndef ONLY_STDINCHI_STDKEY
    inchi_ios_print_nodisplay( f, "  SNon        Exclude stereo (Default: Include Absolute stereo)\n");
    inchi_ios_print_nodisplay( f, "  SRel        Relative stereo\n");
    inchi_ios_print_nodisplay( f, "  SRac        Racemic stereo\n");
    inchi_ios_print_nodisplay( f, "  SUCF        Use Chiral Flag: On means Absolute stereo, Off - Relative\n"); 
    inchi_ios_print_nodisplay( f, "  SUU         Include omitted unknown/undefined stereo\n");
#endif

#ifndef STDINCHI_AS_DEFAULT
#ifndef ONLY_STDINCHI_STDKEY
    inchi_ios_print_nodisplay( f, "  RecMet      Include reconnected metals results (default: not)\n");
#endif
#endif

#ifdef INCLUDE_V102_FEATURES
#ifndef ONLY_STDINCHI_STDKEY
    inchi_ios_print_nodisplay( f, "  RecMetOFF   Do not include reconnected metals results\n");
#endif
#endif

#ifndef STDINCHI_AS_DEFAULT	
    inchi_ios_print_nodisplay( f, "  FixedH      Include Fixed H layer (default: not)\n");
#endif

#ifdef INCLUDE_V102_FEATURES
#ifndef ONLY_STDINCHI_STDKEY
    inchi_ios_print_nodisplay( f, "  FixedHOFF   Do not include Fixed H layer\n");
#endif
#endif    

#if( ADD_PHOSPHINE_STEREO == 1 )
#ifndef STDINCHI_AS_DEFAULT
    inchi_ios_print_nodisplay( f, "  SPXYZ       Include Phosphines Stereochemistry (default: not)\n");
#endif
#endif

#ifdef INCLUDE_V102_FEATURES
#ifndef ONLY_STDINCHI_STDKEY
    inchi_ios_print_nodisplay( f, "  SPXYZOFF    Do not include Phosphines Stereochemistry\n");
#endif
#endif

#if( ADD_ARSINE_STEREO == 1 )
#ifndef STDINCHI_AS_DEFAULT
    inchi_ios_print_nodisplay( f, "  SAsXYZ      Include Arsines Stereochemistry (default: not)\n");
#endif
#endif

#ifdef INCLUDE_V102_FEATURES
#ifndef ONLY_STDINCHI_STDKEY
    inchi_ios_print_nodisplay( f, "  SAsXYZOFF   Do not include Arsines Stereochemistry\n");
#endif
#endif

#if( FIX_ADJ_RAD == 1 )
    inchi_ios_print_nodisplay( f, "  FixRad      Fix Adjacent Radicals\n");
#endif
/*^^^ */
#ifndef STDINCHI_AS_DEFAULT
    inchi_ios_print_nodisplay( f, "  FB          (or FixSp3Bug) Fix bug leading to missing or undefined sp3 parity\n              (default: not)\n" ); 
#endif

#ifdef INCLUDE_V102_FEATURES
#ifndef ONLY_STDINCHI_STDKEY
    inchi_ios_print_nodisplay( f, "  FBOFF       Do not fix bug leading to missing or undefined sp3 parity\n" );
#endif
#endif

#ifndef STDINCHI_AS_DEFAULT
    inchi_ios_print_nodisplay( f, "  FB2         Fix bugs found after v.1.02b release\n              (default: not)\n" );
#endif

#ifdef INCLUDE_V102_FEATURES
#ifndef ONLY_STDINCHI_STDKEY
    inchi_ios_print_nodisplay( f, "  FB2OFF      Do not fix bugs found after v.1.02b release\n" );
#endif
#endif

#ifndef STDINCHI_AS_DEFAULT
    inchi_ios_print_nodisplay( f, "  FNUD        Fix non-uniform drawing issues (default: not)\n");
#endif

#ifdef INCLUDE_V102_FEATURES
#ifndef ONLY_STDINCHI_STDKEY
    inchi_ios_print_nodisplay( f, "  FNUDOFF     Do not fix non-uniform drawing issues\n" );
#endif
#endif

#ifdef INCLUDE_V102_FEATURES
#ifndef ONLY_STDINCHI_STDKEY
    inchi_ios_print_nodisplay( f, "  STDINCHI    Use standard v. 1.02 options (equivalent to %cFixedH\n",INCHI_OPTION_PREFX);
    inchi_ios_print_nodisplay( f, "              %cRecMet %cSPXYZ %cSAsXYZ %cNEWPS %cFB %cFNUD; Absolute stereo)\n",
                                            INCHI_OPTION_PREFX, INCHI_OPTION_PREFX, INCHI_OPTION_PREFX,
                                            INCHI_OPTION_PREFX, INCHI_OPTION_PREFX, INCHI_OPTION_PREFX);
#ifdef STDINCHI_AS_DEFAULT
    inchi_ios_print_nodisplay( f, "              (default)\n");
#endif
#endif

#ifndef ONLY_STDINCHI_STDKEY

    inchi_ios_print_nodisplay( f, "  STDINCHIOFF Turn off standard v. 1.02 options (equivalent to %cFixedHOFF \n",INCHI_OPTION_PREFX);
    inchi_ios_print_nodisplay( f, "             %cRecMetOFF %cSPXYZOFF %cSAsXYZOFF %cNEWPSOFF %cFBOFF %cFB2OFF %cFNUDOFF)\n",
                                            INCHI_OPTION_PREFX, INCHI_OPTION_PREFX, INCHI_OPTION_PREFX,
                                            INCHI_OPTION_PREFX, INCHI_OPTION_PREFX, INCHI_OPTION_PREFX, INCHI_OPTION_PREFX);
#ifndef STDINCHI_AS_DEFAULT
    inchi_ios_print_nodisplay( f, "              (default)\n");
#endif
#endif

#endif



    inchi_ios_print_nodisplay( f, "Conversion\n");
#ifdef INCHI_MAIN
#ifndef ONLY_STDINCHI_STDKEY
    inchi_ios_print_nodisplay( f, "  InChI2InChI  Test mode: Mol/SDfile->%s->%s\n", INCHI_NAME, INCHI_NAME);
#endif
    inchi_ios_print_nodisplay( f, "  InChI2Struct Test mode: Mol/SDfile->%s->Structure->(%s+AuxInfo)\n", INCHI_NAME, INCHI_NAME);
#else
#ifndef ONLY_STDINCHI_STDKEY
    inchi_ios_print_nodisplay( f, "  InChI2InChI  Convert %s string(s) into %s string(s)\n", INCHI_NAME, INCHI_NAME); 
#endif
    inchi_ios_print_nodisplay( f, "  InChI2Struct Convert InChI string(s) to structure(s) in InChI aux.info format\n");
#endif

       
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
#ifdef INCHI_MAIN

#ifdef ONLY_STDINCHI_STDKEY
         "%s ver %s%s.\n\nUsage:\nstdinchi_main inFile [outFile [logFile [problemFile]]] [%coption[ %coption...]]\n", 
#else
         "%s ver %s%s.\n\nUsage:\nInChI_MAIN inputFile [outputFile [logFile [problemFile]]] [%coption[ %coption...]]\n", 
#endif
         INCHI_NAME, INCHI_VERSION, SPECIAL_BUILD_STRING, 
         INCHI_OPTION_PREFX, INCHI_OPTION_PREFX); 

#else
         "%s ver %s%s.\n\nUsage:\nc%s-%s inputFile [outputFile [logFile [problemFile]]] [%coption[ %coption...]]\n", 
         INCHI_NAME, INCHI_VERSION, SPECIAL_BUILD_STRING, 
         INCHI_NAME, INCHI_VERSION, INCHI_OPTION_PREFX, INCHI_OPTION_PREFX); 		 
#endif
    /*^^^ */    

    inchi_ios_print_nodisplay( f, "\nOptions:\n");


    inchi_ios_print_nodisplay( f, "\nInput\n");
    inchi_ios_print_nodisplay( f, "  STDIO       Use standard input/output streams\n");
    /*^^^ post-1.02bv */
    inchi_ios_print_nodisplay( f, "  InpAux      Input structures in %s default aux. info format\n              (for use with STDIO)\n", INCHI_NAME);
    inchi_ios_print_nodisplay( f, "  SDF:DataHeader Read from the input SDfile the ID under this DataHeader\n");
#if( ADD_CMLPP == 1 )
    inchi_ios_print_nodisplay( f, "  CML         Input in CML format (default for input file .CML extension)\n");
#endif


    inchi_ios_print_nodisplay( f, "Output\n");
    inchi_ios_print_nodisplay( f, "  AuxNone     Omit auxiliary information (default: Include)\n");
    inchi_ios_print_nodisplay( f, "  NoLabels    Omit structure number, DataHeader and ID from %s output\n", INCHI_NAME);
    inchi_ios_print_nodisplay( f, "  Tabbed      Separate structure number, %s, and AuxInfo with tabs\n", INCHI_NAME);
#ifndef ONLY_STDINCHI_STDKEY
    inchi_ios_print_nodisplay( f, "  Compress    Compressed output\n");
#endif
#if( defined(_WIN32) && defined(_MSC_VER) && !defined(INCHI_ANSI_ONLY) && !defined(INCHI_LIBRARY) )
    inchi_ios_print_nodisplay( f, "  D           Display the structures\n");
    inchi_ios_print_nodisplay( f, "  EQU         Display sets of identical components\n");
    inchi_ios_print_nodisplay( f, "  Fnumber     Set display Font size in number of points\n");
#endif
    inchi_ios_print_nodisplay( f, "  OutputSDF   Convert %s created with default aux. info to SDfile\n", INCHI_NAME);
#if ( SDF_OUTPUT_DT == 1 )
    inchi_ios_print_nodisplay( f, "  SdfAtomsDT  Output Hydrogen Isotopes to SDfile as Atoms D and T\n");
#endif

    
    inchi_ios_print_nodisplay( f, "Structure perception\n");
#ifndef STDINCHI_AS_DEFAULT
#ifndef ONLY_STDINCHI_STDKEY
    inchi_ios_print_nodisplay( f, "  NEWPS       Narrow end of wedge points to stereocenter (default: both)\n");
#endif
#endif
#ifdef INCLUDE_V102_FEATURES
    inchi_ios_print_nodisplay( f, "  NEWPSOFF    Both ends of wedge point to stereocenters\n");
#endif
    inchi_ios_print_nodisplay( f, "  DoNotAddH   Don't add H according to usual valences: all H are explicit\n");

#ifndef ONLY_STDINCHI_STDKEY
    inchi_ios_print_nodisplay( f, "  NoADP       Disable Aggressive Deprotonation (for testing only)\n");
#endif
#ifdef ONLY_STDINCHI_STDKEY
    inchi_ios_print_nodisplay( f, "  SNon        Exclude stereo\n");
#endif

    inchi_ios_print_nodisplay( f, "Generation\n");
    inchi_ios_print_nodisplay( f, "  Wnumber     Set time-out per structure in seconds; W0 means unlimited\n");
    inchi_ios_print_nodisplay( f, "  WarnOnEmptyStructure Warn and produce empty %s for empty structure\n", INCHI_NAME);
#ifdef INCLUDE_V102_FEATURES
    inchi_ios_print_nodisplay( f, "  Key         Generate InChIKey\n");
#endif
/*^^^ */
#ifndef ONLY_STDINCHI_STDKEY
    inchi_ios_print_nodisplay( f, "  SNon        Exclude stereo (Default: Include Absolute stereo)\n");
    inchi_ios_print_nodisplay( f, "  SRel        Relative stereo\n");
    inchi_ios_print_nodisplay( f, "  SRac        Racemic stereo\n");
    inchi_ios_print_nodisplay( f, "  SUCF        Use Chiral Flag: On means Absolute stereo, Off - Relative\n"); 
    inchi_ios_print_nodisplay( f, "  SUU         Include omitted unknown/undefined stereo\n");
#endif

#ifndef STDINCHI_AS_DEFAULT
#ifndef ONLY_STDINCHI_STDKEY
    inchi_ios_print_nodisplay( f, "  RecMet      Include reconnected metals results (default: not)\n");
#endif
#endif

#ifdef INCLUDE_V102_FEATURES
#ifndef ONLY_STDINCHI_STDKEY
    inchi_ios_print_nodisplay( f, "  RecMetOFF   Do not include reconnected metals results\n");
#endif
#endif

#ifndef STDINCHI_AS_DEFAULT	
#ifndef ONLY_STDINCHI_STDKEY
    inchi_ios_print_nodisplay( f, "  FixedH      Include Fixed H layer (default: not)\n");
#endif
#endif

#ifdef INCLUDE_V102_FEATURES
#ifndef ONLY_STDINCHI_STDKEY
    inchi_ios_print_nodisplay( f, "  FixedHOFF   Do not include Fixed H layer\n");
#endif
#endif

#if( ADD_PHOSPHINE_STEREO == 1 )
#ifndef STDINCHI_AS_DEFAULT
#ifndef ONLY_STDINCHI_STDKEY
    inchi_ios_print_nodisplay( f, "  SPXYZ       Include Phosphines Stereochemistry (default: not)\n");
#endif
#endif
#endif

#ifdef INCLUDE_V102_FEATURES
#ifndef ONLY_STDINCHI_STDKEY
    inchi_ios_print_nodisplay( f, "  SPXYZOFF    Do not include Phosphines Stereochemistry\n");
#endif
#endif

#if( ADD_ARSINE_STEREO == 1 )
#ifndef STDINCHI_AS_DEFAULT
#ifndef ONLY_STDINCHI_STDKEY
    inchi_ios_print_nodisplay( f, "  SAsXYZ      Include Arsines Stereochemistry (default: not)\n");
#endif
#endif
#endif

#ifdef INCLUDE_V102_FEATURES
#ifndef ONLY_STDINCHI_STDKEY
    inchi_ios_print_nodisplay( f, "  SAsXYZOFF   Do not include Arsines Stereochemistry\n");
#endif
#endif

#if( FIX_ADJ_RAD == 1 )
    inchi_ios_print_nodisplay( f, "  FixRad      Fix Adjacent Radicals\n");
#endif
/*^^^ */
#ifndef STDINCHI_AS_DEFAULT
#ifndef ONLY_STDINCHI_STDKEY
    inchi_ios_print_nodisplay( f, "  FB          (or FixSp3Bug) Fix bug leading to missing or undefined sp3 parity\n              (default: not)\n" );
#endif
#endif

#ifdef INCLUDE_V102_FEATURES
#ifndef ONLY_STDINCHI_STDKEY
    inchi_ios_print_nodisplay( f, "  FBOFF       Do not fix bug leading to missing or undefined sp3 parity\n" );
#endif
#endif

#ifndef STDINCHI_AS_DEFAULT
#ifndef ONLY_STDINCHI_STDKEY
    inchi_ios_print_nodisplay( f, "  FB2         Fix bugs found after v.1.02b release\n              (default: not)\n" );
#endif
#endif

#ifdef INCLUDE_V102_FEATURES
#ifndef ONLY_STDINCHI_STDKEY
    inchi_ios_print_nodisplay( f, "  FB2OFF      Do not fix bugs found after v.1.02b release\n" );
#endif
#endif

#ifndef STDINCHI_AS_DEFAULT
#ifndef ONLY_STDINCHI_STDKEY
    inchi_ios_print_nodisplay( f, "  FNUD        Fix non-uniform drawing issues (default: not)\n");
#endif
#endif

#ifdef INCLUDE_V102_FEATURES
#ifndef ONLY_STDINCHI_STDKEY
    inchi_ios_print_nodisplay( f, "  FNUDOFF     Do not fix non-uniform drawing issues\n" );
#endif
#endif

#ifdef INCLUDE_V102_FEATURES

#ifndef ONLY_STDINCHI_STDKEY
    inchi_ios_print_nodisplay( f, "  STDINCHI    Use standard v. 1.02 options (equivalent to %cFixedH\n",INCHI_OPTION_PREFX);
    inchi_ios_print_nodisplay( f, "              %cRecMet %cSPXYZ %cSAsXYZ %cNEWPS %cFB %cFNUD; Absolute stereo)\n",
                                            INCHI_OPTION_PREFX, INCHI_OPTION_PREFX, INCHI_OPTION_PREFX,
                                            INCHI_OPTION_PREFX, INCHI_OPTION_PREFX, INCHI_OPTION_PREFX);
#ifdef STDINCHI_AS_DEFAULT
    inchi_ios_print_nodisplay( f, "              (default)\n");
#endif
#endif

#ifndef ONLY_STDINCHI_STDKEY
    inchi_ios_print_nodisplay( f, "  STDINCHIOFF Turn off standard v. 1.02 options (equivalent to %cFixedHOFF \n",INCHI_OPTION_PREFX);
    inchi_ios_print_nodisplay( f, "             %cRecMetOFF %cSPXYZOFF %cSAsXYZOFF %cNEWPSOFF %cFBOFF %cFB2OFF %cFNUDOFF)\n",
                                            INCHI_OPTION_PREFX, INCHI_OPTION_PREFX, INCHI_OPTION_PREFX,
                                            INCHI_OPTION_PREFX, INCHI_OPTION_PREFX, INCHI_OPTION_PREFX, INCHI_OPTION_PREFX);
#ifndef STDINCHI_AS_DEFAULT
    inchi_ios_print_nodisplay( f, "              (default)\n");
#endif
#endif

#endif


       
#endif /* #if ( bRELEASE_VERSION == 1 ) */
}




#define fprintf2 inchi_fprintf

#ifndef INCHI_LIBRARY
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
        fprintf2( stderr, "%s version %s%s%s\n", INCHI_NAME, INCHI_VERSION, SPECIAL_BUILD_STRING, bRELEASE_VERSION? "":" (For pre-release testing)" );
        fprintf2( stderr, "Log file not specified. Using standard error output.\n");
        *log_file = stderr;
    } else
    if ( !(*log_file = fopen( ip->path[2], "w" ) ) ) {
        fprintf2( stderr, "%s version %s%s%s\n", INCHI_NAME, INCHI_VERSION, SPECIAL_BUILD_STRING, bRELEASE_VERSION? "":" (For pre-release testing)" );
        fprintf2( stderr, "Cannot open log file '%s'. Using standard error output.\n", ip->path[2] );
        *log_file = stderr;
    } else {
        fprintf2( *log_file, "%s version %s%s%s\n", INCHI_NAME, INCHI_VERSION, SPECIAL_BUILD_STRING, bRELEASE_VERSION? "":" (For pre-release testing)" );
        fprintf2( *log_file, "Opened log file '%s'\n", ip->path[2] );
    }
    /* input file */
    if ( (ip->nInputType == INPUT_MOLFILE || ip->nInputType == INPUT_SDFILE ||
         ip->nInputType == INPUT_CMLFILE  || ip->nInputType == INPUT_INCHI  ||
         ip->nInputType == INPUT_INCHI_PLAIN ) && ip->num_paths > 0 ) 
    {
        const char *fmode = NULL;
#if( defined(_MSC_VER)&&defined(_WIN32) || defined(__BORLANDC__)&&defined(__WIN32__) || defined(__GNUC__)&&defined(__MINGW32__)&&defined(_WIN32) )
        /* compilers that definitely allow fopen "rb" (binary read) mode */
        fmode = "rb";
        if ( !ip->path[0] || !ip->path[0][0] || !(*inp_file = fopen( ip->path[0], "rb" ) ) ) {
            fprintf2( *log_file, "Cannot open input file '%s'. Terminating.\n", ip->path[0]? ip->path[0] : "<No name>" );
            goto exit_function;
        } else
        if ( ip->nInputType == INPUT_CMLFILE ) {
            int c;
#ifdef CML_DEBUG
            printf ("cr %d lf %d ret %d\n", (int) '\r', (int) '\f', (int) '\n');
#endif
            /* read up to the end of the first line */
            while( (c = fgetc( *inp_file )) && c != EOF && c != '\n' && c != '\r' )
                ;
            if ( c == '\r' || c == EOF ) {
                /* text file contains CR; close and reopen as "text" */
                fclose( *inp_file );
                if ( !(*inp_file = fopen( ip->path[0], "r" ) ) ) {
                    fprintf2( *log_file, "Cannot open input file '%s' (2nd attempt). Terminating.\n", ip->path[0] );
                    goto exit_function;
                }
                fprintf2( *log_file, "Opened input file '%s'\n", ip->path[0] );
                fmode = "r";
            } else {
                fclose( *inp_file );
                if ( !(*inp_file = fopen( ip->path[0], "rb" ) ) ) {
                    fprintf2( *log_file, "Cannot open input file '%s' (2nd attempt). Terminating.\n", ip->path[0] );
                    goto exit_function;
                }
                fprintf2( *log_file, "Opened input file '%s': no CR.\n", ip->path[0] );
                fmode = "rb";
            }
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
#if( defined(_MSC_VER)&&defined(_WIN32) || defined(__BORLANDC__)&&defined(__WIN32__) || defined(__GNUC__)&&defined(__MINGW32__)&&defined(_WIN32) )
        if ( ip->nInputType != INPUT_CMLFILE ) {
            fmode = "wb";
        }
#endif
        if (  !(*prb_file = fopen( ip->path[3], fmode ) ) ) {
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

#endif /* INCHI_LIBRARY */

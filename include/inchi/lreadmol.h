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


/* local prototypes */
int bypass_sdf_data_items( FILE* inp, long *cas_reg_no, char* comment, int lcomment, char *name, int lname, int prev_err,
                           const char *pSdfLabel, char *pSdfValue, char *pStrErr );
MOL_DATA* read_mol_file( FILE* inp, MOL_HEADER_BLOCK *OnlyHeaderBlock, MOL_CTAB *OnlyCtab,
                         int bGetOrigCoord, int *err, char *pStrErr );


static int mol_read_hdr(MOL_HEADER_BLOCK *hdr, FILE* inp, char *pStrErr);
static int mol_read_counts_line( MOL_CTAB* ctab, FILE *inp, char *pStrErr );
static int read_atom_block( MOL_CTAB* ctab, FILE *inp, int err, char *pStrErr );
static int read_bonds_block( MOL_CTAB* ctab, FILE *inp, int err, char *pStrErr );
static int read_stext_block( MOL_CTAB* ctab, FILE *inp, int err, char *pStrErr );
static int read_properties_block( MOL_CTAB* ctab, MOL_HEADER_BLOCK *pHdr, FILE *inp, int err, char *pStrErr );

static int identify_sdf_label( char* inp_line, const char *pSdfLabel );
static long extract_cas_rn( char *line );
static int mol_copy_check_empty( char* dest, char* source, int len, char **first_space );
static int mol_read_datum(void* data, int field_len, int  data_type, char** line_ptr);

static int RemoveNonPrintable( char *line );


/******/
#ifndef MOLFILE_ERR_FIN
#define MOLFILE_ERR_FIN(err, new_err, err_fin, msg) \
        if ( !(err) && (new_err) ) { (err) = (new_err);} AddMOLfileError(pStrErr, (msg)); goto err_fin
#endif
#ifndef MOLFILE_ERR_SET
#define MOLFILE_ERR_SET(err, new_err, msg) \
        if ( !(err) && (new_err) ) { (err) = (new_err);} AddMOLfileError(pStrErr, (msg))
#endif

/*************************************************************************/
int AddMOLfileError( char *pStrErr, const char *szMsg )
{
    if ( pStrErr && szMsg && szMsg[0] ) {
        int lenStrErr = strlen( pStrErr );
        int lenMsg    = strlen( szMsg );
        char *p = strstr( pStrErr, szMsg );
        if ( p && (p==pStrErr || *(p-1) == ' ' && (*(p-2) == ';' || *(p-2) == ':' )) &&
                  (p+lenMsg == pStrErr+lenStrErr || 
                  p[lenMsg] == ';' && p[lenMsg+1] == ' ' ||
                  p[lenMsg-1]==':' && p[lenMsg]==' ') ) {
            return 1; /*  reject duplicates */
        }
        if ( lenStrErr + lenMsg + 2*(lenStrErr > 0) < STR_ERR_LEN ) {
            /*  enough room to add */
            if (lenStrErr > 0) {
                if ( pStrErr[lenStrErr-1] != ':' ) {
                    strcat( pStrErr, ";" );
                }
                strcat( pStrErr, " " );
            }
            strcat( pStrErr, szMsg );
            return 1;
        }
        /*  no room */
        if ( strstr( pStrErr, "..." ) ) {
            return 0; /*  no room mark has already been set */
        }
        if ( lenStrErr + 3 < STR_ERR_LEN ) {
            strcat( pStrErr, "..." );
        }
    }
    return 0;
}
/*************** static **********************************************************/
int mol_copy_check_empty( char* dest, char* source, int len, char **first_space )
{
    int i, c;   /* required len >= 0; dest must have at least len+1 bytes */
    if ( len > 0 )
        strncpy( dest, source, len );
    dest[len]='\0'; 
    len = ( len > 0 )? (int)strlen( dest) : 0;
    for ( i = (len-1); i >= 0 && 0 != (c = source[i]) && isspace(UCINT c); i-- )
        ;
    *first_space = dest + (i+1); /* first blank or zero terminating byte in dest */
    return len; /* number of actually processed bytes; zero termination not included */
}
/************* static ************************************************************/
int mol_read_datum(void* data, int field_len, int  data_type, char** line_ptr)
{
/* 1. 'field_len' for MOL_STRING_DATA does not include trailing zero,
 *     that is actual length of the string pointed by 'data'
 *     should be at least field_len+1 bytes.
 *     For numerical data 'field_len' is length of input data field
 *     For numerical integral data field_len <= 0 means read up to first
 *     non-numeric character as strtod() does ("free format")
 * 2.  return value: for MOL_STRING_DATA: number of bytes excluding trailing zero
 *                   for all others:  1=success; 0 = empty; -1= error
 * 3.  on exit *line_ptr points to the next byte after the last entered
 */ 
    char *p = *line_ptr, *q, *p_end;
    int  i, ret=1, c, len;
    long   ldata;
    double ddata;
    
    switch( data_type ) {
    case MOL_STRING_DATA:
        for ( i= 0; i < field_len && 0 != (c = p[i]) && isspace(UCINT c); i++ ) /* pass by all leading spaces */
            ;
        len = mol_copy_check_empty( (char*)data, &p[i], field_len-i, &q );
        ret = ( q - (char*)data );/* actual data length */
        *q = '\0';                /* add zero termination to data if it is not there yet*/
        *line_ptr += (len+i);     /* ptr to the 1st byte of the next input field or to zero termination */
        break;

    case MOL_CHAR_INT_DATA:
    case MOL_SHORT_INT_DATA:
    case MOL_LONG_INT_DATA:
        { /* block start */
            char str[MOL_MAX_VALUE_LEN+1];
            ldata = 0L;
            if ( field_len > MOL_MAX_VALUE_LEN ) {
                ret = -1;
            }else
            if ( field_len > 0 ) { /* fixed length */
                *line_ptr += ( len = mol_copy_check_empty( str, p, field_len, &q ) );
                *q = '\0';
                if ( !len || !(q-str) ) { /* empty string */
                    ret = 0;
                }else
                if ( (ldata=strtol(str,&p_end,10), p_end != q) ){ /* wrong data: incompletely interpreted */
                    ret = -1;
                }
            }else{  /* free format: field_len <= 0 */
                ldata = strtol( p, &p_end, 10 );
                *line_ptr += ( len = p_end - p );
                if ( len == 0 ){
                    ret = 0;
                }
            }    
            
            switch( data_type ) {
            case MOL_CHAR_INT_DATA:
                if ( SCHAR_MIN <= ldata  && ldata <= SCHAR_MAX ){ /* from || to &&: 11-19-96 */
                    *(S_CHAR*)data = (S_CHAR)ldata;
                }else{
                    *(S_CHAR*)data = (S_CHAR)0;
                    ret = -1;
                }
                break;
            case MOL_SHORT_INT_DATA:
                if ( SHRT_MIN <= ldata && ldata <= SHRT_MAX ){
                    *(S_SHORT*)data = (S_SHORT)ldata;
                }else{
                    *(S_SHORT*)data = (S_SHORT)0;
                    ret = -1;
                }
                break;
            case MOL_LONG_INT_DATA:
                if ( LONG_MIN < ldata && ldata < LONG_MAX ){
                    *(long*)data = (long)ldata;
                }else{
                    *(long*)data = 0L;
                    ret = -1;
                }
                break;
            default:
                ret=-1;
            }
            
        } /* block end */
        break;
    case MOL_DOUBLE_DATA:
    case MOL_FLOAT_DATA:
        { /* block start */
            char str[MOL_MAX_VALUE_LEN+1];
            if ( field_len > MOL_MAX_VALUE_LEN ) {
                ret = -1;
                ddata = 0.0;
            }else
            if ( field_len > 0 ) {
                *line_ptr += (len = mol_copy_check_empty( str, p, field_len, &q ));
                *q = '\0';
                if ( !len || !(q-str) ) { /* empty string */
                    ddata = 0.0;
                    ret   = 0;
                }else
                if ( (ddata=strtod(str,&p_end), p_end != q) ){ /* wrong data */
                    ret = -1;
                }
            }else{ /* free format */
                ddata = strtod( p, &p_end );
                *line_ptr += ( len = p_end - p );
                if ( len == 0 ){
                    ret = 0;
                }
            } 
            switch(data_type){
            case MOL_DOUBLE_DATA:
                if ( ddata != HUGE_VAL && /*ldata*/ ddata != -HUGE_VAL ){ /* replaced ldata with ddata 6-30-98 DCh */
                    *(double*)data = ddata;
                }else{
                    *(double*)data = 0.0;
                    ret = -1;
                }
                break;
            case MOL_FLOAT_DATA:
                if ( fabs(ddata) <= (double)FLT_MIN ) {
                    *(float*)data = 0.0;
                }else
                if ( fabs(ddata) >= (double)FLT_MAX ) {
                    *(float*)data = 0.0;
                    ret = -1;
                }else{
                    *(float*)data = (float)ddata;
                }
                break;
            }
        } /* block end */
        break;
    case MOL_JUMP_TO_RIGHT:
        for ( i = 0; i < field_len && p[i]; i++ )
            ;
        *line_ptr += i;
        ret = i;
        break;
    default:
        ret = -1;
    }
    return ret;
}
/************* static ************************************************************/
int mol_read_hdr(MOL_HEADER_BLOCK *hdr, FILE* inp, char *pStrErr)
{   
    /* All input lines can have are up 80 characters */
    /* Header Block */
    char line[MOLFILEINPLINELEN]; /* + cr +lf +zero termination + reserve */
    int  err = 0, len;
    const int  line_len = sizeof(line);
    char *p;
    
    /* memset( &hdr, 0, sizeof( MOL_HEADER_BLOCK ) ); */
    /*------------ header line #1: name ----------------*/
    if ( NULL == ( p = fgets_up_to_lf( line, line_len, inp ) ) ){
        err = 1;             /* can't read the input file line */
        /* AddMOLfileError( pStrErr, "Can't read header block name line" ); */
        goto err_fin;
    }
    remove_one_lf( line );
    /* -- Disabled to relax strictness: allow > 80 chars names.
    if ( line[MOLFILEMAXLINELEN] ){
        err = 2;             // too long line
        goto err_fin;
    }
    */ 
    len = mol_read_datum( hdr->szMoleculeName, sizeof(hdr->szMoleculeName)-1, MOL_STRING_DATA, &p );
    /*----------- header line #2 -----------------------*/
    if ( NULL == ( p = fgets_up_to_lf( line, line_len, inp ) ) ){
        err = 3;             /* can't read the input file line */
        /* AddMOLfileError( pStrErr, "Can't read header block line 2" ); */
        goto err_fin;
    }
    remove_one_lf( line );
    /* -- Disabled to relax strictness: allow > 80 chars names.
    if ( line[MOLFILEMAXLINELEN] ){
        err = 4;             // too long input file line
        goto err_fin;
    } 
    */
    len = mol_read_datum( hdr->szUserInitials, sizeof(hdr->szUserInitials)-1, MOL_STRING_DATA, &p );
    len = mol_read_datum( hdr->szProgramName,  sizeof(hdr->szProgramName)-1,  MOL_STRING_DATA, &p );

    /*------------ Relax strictness -----------------------*/
    len = mol_read_datum( &hdr->cMonth,                  2,  MOL_CHAR_INT_DATA,  &p );
    len = mol_read_datum( &hdr->cDay,                    2,  MOL_CHAR_INT_DATA,  &p );
    len = mol_read_datum( &hdr->cYear,                   2,  MOL_CHAR_INT_DATA,  &p );
    len = mol_read_datum( &hdr->cHour,                   2,  MOL_CHAR_INT_DATA,  &p );
    len = mol_read_datum( &hdr->cMinute,                 2,  MOL_CHAR_INT_DATA,  &p );
    len = mol_read_datum( hdr->szDimCode, sizeof(hdr->szDimCode)-1,  MOL_STRING_DATA, &p );
    len = mol_read_datum( &hdr->nScalingFactor1,         2,  MOL_SHORT_INT_DATA, &p );
    len = mol_read_datum( &hdr->dScalingFactor2,        10,  MOL_DOUBLE_DATA,    &p );
    len = mol_read_datum( &hdr->dEnergy,                12,  MOL_DOUBLE_DATA,    &p );
    len = mol_read_datum( &hdr->lInternalRegistryNumber, 6, MOL_LONG_INT_DATA,   &p );

    /* save the whole line 2 */
    p = line;
    len = mol_read_datum( hdr->szMoleculeLine2, sizeof(hdr->szMoleculeLine2)-1, MOL_STRING_DATA, &p );

    
    /*------------ header line #3: comment ----------------*/
    if ( NULL == ( p = fgets_up_to_lf( line, line_len, inp ) ) ){
        err = 7;             /* can't read the line */
        /* AddMOLfileError( pStrErr, "Can't read header block comment line" ); */
        goto err_fin;
    }
    remove_one_lf( line );
    /* -- Disabled to relax strictness: allow > 80 chars comments.
    if ( line[MOLFILEMAXLINELEN] ){
        err = 8;             // too long line
        goto err_fin;
    }
    */ 
    len = mol_read_datum( hdr->szComment, sizeof(hdr->szComment)-1, MOL_STRING_DATA, &p );

err_fin:

    return err;    
}
/********** static *****************************************************/
int RemoveNonPrintable( char *line )
{
    int i, c, num = 0;
    if ( line ) {
        for ( i = 0; c = UCINT line[i]; i ++ ) {
            /* assuming ASCII charset */
            if ( c < ' ' || c >= 0x7F ) {
                line[i] = '.';
                num ++;
            }
        }
    }
    return num;
}
/************** static *************************************************/
int mol_read_counts_line( MOL_CTAB* ctab, FILE *inp, char *pStrErr )
{
    char *p;
    char line[MOLFILEINPLINELEN];
    const int line_len = sizeof(line);
    int   err = 0, len;
    
    if ( NULL == ( p = fgets_up_to_lf( line, line_len, inp ) ) ){
        MOLFILE_ERR_FIN (err, 1, err_fin, "Cannot read counts line");
        /* can't read the input file line */
    }
    remove_one_lf( line );
    if ( line[MOLFILEMAXLINELEN] ){
        MOLFILE_ERR_SET (err, 0, "Too long counts line");  /* too long input file line */
    } 
    if (    0 > mol_read_datum( &ctab->nNumberOfAtoms,         3,  MOL_SHORT_INT_DATA, &p )
         || 0 > mol_read_datum( &ctab->nNumberOfBonds,         3,  MOL_SHORT_INT_DATA, &p )
#if ( MOL_QUERY == MOL_PRESENT )
         || 0 > mol_read_datum( &ctab->nNumberOfAtomsLists,    3,  MOL_SHORT_INT_DATA, &p )
#else
         || 0 > mol_read_datum( NULL,                          3,  MOL_JUMP_TO_RIGHT,  &p )
#endif
         || 0 > mol_read_datum( NULL, /*obsolete*/             3,  MOL_JUMP_TO_RIGHT,  &p )
         || 0 > mol_read_datum( &ctab->cChiralFlag,            3,  MOL_CHAR_INT_DATA,  &p )
         || 0 > mol_read_datum( &ctab->nNumberOfStextEntries,  3,  MOL_SHORT_INT_DATA, &p )
#if ( MOL_CPSS == MOL_PRESENT )         
         || 0 > mol_read_datum( &ctab->nNumberOfReactionComponentsPlus1, 3, MOL_SHORT_INT_DATA, &p )
         || 0 > mol_read_datum( &ctab->nNumberOfReactants,     3,  MOL_SHORT_INT_DATA, &p )
         || 0 > mol_read_datum( &ctab->nNumberOfProducts,      3,  MOL_SHORT_INT_DATA, &p )
         || 0 > mol_read_datum( &ctab->nNumberOfIntermediates, 3,  MOL_SHORT_INT_DATA, &p )
#else
         || 0 > mol_read_datum( NULL,                          3,  MOL_JUMP_TO_RIGHT,  &p )
         || 0 > mol_read_datum( NULL,                          3,  MOL_JUMP_TO_RIGHT,  &p )
         || 0 > mol_read_datum( NULL,                          3,  MOL_JUMP_TO_RIGHT,  &p )
         || 0 > mol_read_datum( NULL,                          3,  MOL_JUMP_TO_RIGHT,  &p )
#endif
         || 0 > mol_read_datum( &ctab->nNumberOfPropertyLines, 3,  MOL_SHORT_INT_DATA, &p ) ){
        err = 3;  /* can't interpret counts line */
        MOLFILE_ERR_SET (err, 3, "Cannot interpret counts line:");  /* too long input file line */
        RemoveNonPrintable( line );
        AddMOLfileError(pStrErr, line);
        goto err_fin;
    }
    len = mol_read_datum( ctab->csCurrentCtabVersion, sizeof(ctab->csCurrentCtabVersion)-1, MOL_STRING_DATA, &p );
err_fin:
    return err;
}

/************ static *************************************************************/
int read_atom_block( MOL_CTAB* ctab, FILE *inp, int err, char *pStrErr )
{
    char *p;
    char line[MOLFILEINPLINELEN];
    const int line_len = sizeof(line);
    S_SHORT i, chg;
    static S_SHORT charge_val[] = {0, 3, 2, 1, 'R', -1, -2, -3};
    /*                           0  1  2  3   4    5   6   7 */ 
    /*
      if ( NULL == ctab->MolAtom ){
          err = 1;
          goto err_fin; // internal error: MolAtom structure has not been allocated
      }
     */
    
    for ( i = 0; i < ctab->nNumberOfAtoms; i++ ) {
        
        if ( NULL == ( p = fgets_up_to_lf( line, line_len, inp ) ) ){
            if ( !err ) {
                MOLFILE_ERR_SET (err, 2, "Cannot read atom block line");
            }
            break;
        }
        remove_one_lf( line );
        if ( line[MOLFILEMAXLINELEN] ){
            MOLFILE_ERR_SET (err, 0, "Too long atom block line");
        }
        if ( err ) {
            if ( !strcmp( line, SDF_END_OF_DATA ) ) {
                err = -abs(err);
                break;
            }
            continue; /* bypass the rest of the Atom block */
        }
        if ( NULL != ctab->szCoord ) {
            mystrncpy( ctab->szCoord[i], p, 31 ); /* original coordinates */
        }

        if ( NULL != ctab->MolAtom ) {
            if (   0 > mol_read_datum( &ctab->MolAtom[i].fX,   10,  MOL_DOUBLE_DATA, &p )
                || 0 > mol_read_datum( &ctab->MolAtom[i].fY,   10,  MOL_DOUBLE_DATA, &p )
                || 0 > mol_read_datum( &ctab->MolAtom[i].fZ,   10,  MOL_DOUBLE_DATA, &p )
                || 0 > mol_read_datum( NULL, /* undescribed in article*/    1,  MOL_JUMP_TO_RIGHT, &p )
                || 0 == mol_read_datum( &ctab->MolAtom[i].szAtomSymbol,     3,  MOL_STRING_DATA, &p ) /* was sizeof(ctab->MolAtom[0].szAtomSymbol)-1 */
#ifdef INCHI_MAIN
                || 0 > mol_read_datum( &ctab->MolAtom[i].cMassDifference,   2,  MOL_SHORT_INT_DATA, &p )
#else
                || 0 > mol_read_datum( &ctab->MolAtom[i].cMassDifference,   2,  MOL_CHAR_INT_DATA, &p )
#endif
                || 0 > mol_read_datum( &ctab->MolAtom[i].cCharge,           3,  MOL_CHAR_INT_DATA, &p )
                || 0 > mol_read_datum( &ctab->MolAtom[i].cStereoParity,     3,  MOL_CHAR_INT_DATA, &p )
#if ( MOL_QUERY == MOL_PRESENT )
                || 0 > mol_read_datum( &ctab->MolAtom[i].cH_countPlus1,     3,  MOL_CHAR_INT_DATA, &p )
                || 0 > mol_read_datum( &ctab->MolAtom[i].cStereoCare,       3,  MOL_CHAR_INT_DATA, &p )
#else
                || 0 > mol_read_datum( NULL,                                3,  MOL_JUMP_TO_RIGHT,  &p )
                || 0 > mol_read_datum( NULL,                                3,  MOL_JUMP_TO_RIGHT,  &p )
#endif
                || 0 > mol_read_datum( &ctab->MolAtom[i].cValence,          3,  MOL_CHAR_INT_DATA, &p ) ) {
                
                err = 4;
                MOLFILE_ERR_SET (err, 4, "Cannot interpret atom block line:");
                RemoveNonPrintable( line );
                AddMOLfileError(pStrErr, line);
                if ( !strcmp( line, SDF_END_OF_DATA ) ) {
                    err = -abs(err);
                    break;
                }
                continue; /* can't interpret a first half of atom block line */
            }
            if ( 2 == strlen(ctab->MolAtom[i].szAtomSymbol) && isupper(UCINT ctab->MolAtom[i].szAtomSymbol[1]))
                ctab->MolAtom[i].szAtomSymbol[1] = (char)tolower(UCINT ctab->MolAtom[i].szAtomSymbol[1]); /* 5-4-99 DCh*/

            if ( (chg = (S_SHORT) ctab->MolAtom[i].cCharge)< 0 || chg >= (int)(sizeof ( charge_val ) / sizeof( charge_val[0] )) ) {
                /* ctab->MolAtom[i].cCharge = 0; */ /* error; ignore for now */
                ctab->MolAtom[i].cCharge  = (S_CHAR)(4 - chg); /*  allow greater charges to accommodate NCI structures. 8-20-2002 */
                ctab->MolAtom[i].cRadical = 0;
            }else
            if ( 'R' == (chg = charge_val[chg]) ){
                ctab->MolAtom[i].cCharge  = 0;
                ctab->MolAtom[i].cRadical = RADICAL_DOUBLET;
            }else{
                ctab->MolAtom[i].cCharge  = (S_CHAR)chg; /* actual charge value */
                ctab->MolAtom[i].cRadical = 0;
            }
#ifdef INCHI_MAIN
            if ( ctab->MolAtom[i].cMassDifference ) { /* e_ReadMOL.c specific */
                ctab->MolAtom[i].cMassDifference += ISOTOPIC_SHIFT_FLAG;
            }
#endif             

            if (
#if ( MOL_CPSS == MOL_PRESENT )
                   0 > mol_read_datum( &ctab->MolAtom[i].cH0_designator,           3,  MOL_CHAR_INT_DATA, &p )
                || 0 > mol_read_datum( &ctab->MolAtom[i].cReactionComponentType,   3,  MOL_CHAR_INT_DATA, &p )
                || 0 > mol_read_datum( &ctab->MolAtom[i].cReactionComponentNumber, 3,  MOL_CHAR_INT_DATA, &p )
#else
                   0 > mol_read_datum( NULL,                                       3,  MOL_JUMP_TO_RIGHT,  &p )
                || 0 > mol_read_datum( NULL,                                       3,  MOL_JUMP_TO_RIGHT,  &p )
                || 0 > mol_read_datum( NULL,                                       3,  MOL_JUMP_TO_RIGHT,  &p )
#endif        
#if ( MOL_REACT == MOL_PRESENT )
                || 0 > mol_read_datum( &ctab->MolAtom[i].nAtomAtomMappingNumber,   3,  MOL_SHORT_INT_DATA, &p )
                || 0 > mol_read_datum( &ctab->MolAtom[i].cReactionComponentType,   3,  MOL_CHAR_INT_DATA, &p )
#else
                || 0 > mol_read_datum( NULL,                                       3,  MOL_JUMP_TO_RIGHT,  &p )
                || 0 > mol_read_datum( NULL,                                       3,  MOL_JUMP_TO_RIGHT,  &p )
#endif        
#if ( MOL_REACT == MOL_PRESENT || MOL_QUERY == MOL_PRESENT )
                || 0 > mol_read_datum( &ctab->MolAtom[i].cExactChargeFlag,         3,  MOL_CHAR_INT_DATA, &p )
#else
                || 0 > mol_read_datum( NULL,                                       3,  MOL_JUMP_TO_RIGHT,  &p )
#endif
            ){
                err = 5; /* can't interpret a second half of atom block line */
                MOLFILE_ERR_SET (err, 5, "Cannot interpret atom block line:");
                RemoveNonPrintable( line );
                AddMOLfileError(pStrErr, line);
                if ( !strcmp( line, SDF_END_OF_DATA ) ) {
                    err = -abs(err);
                    break;
                }
                continue;
            }
        }
    }
/* err_fin: */
    return err;
}
/************ static *************************************************************/
int read_bonds_block( MOL_CTAB* ctab, FILE *inp, int err, char *pStrErr )
{
    char *p;
    char line[MOLFILEINPLINELEN];
    const int line_len = sizeof(line);
    S_SHORT i;
    /*
      if ( NULL == ctab->MolBond ){
          err = 1;
          goto err_fin;    // internal error: memory has not been allocated for MolBond structure
      }
     */
    for ( i = 0; i < ctab->nNumberOfBonds; i++ ) {
        
        if ( NULL == ( p = fgets_up_to_lf( line, line_len, inp ) ) ){
            if ( !err ) {
                MOLFILE_ERR_SET (err, 2, "Cannot read bond block line");
            }
            break;
        }
        remove_one_lf( line );
        if ( line[MOLFILEMAXLINELEN] ){
            err = err? err : 3;             /* too long input file line */
        }
        if ( err ) {
            if ( !strcmp( line, SDF_END_OF_DATA ) ) {
                err = -abs(err);
                break;
            }
            continue;
        }

        if ( ctab->MolBond ) {
            if (   0 > mol_read_datum( &ctab->MolBond[i].nAtomNo1,      3,  MOL_SHORT_INT_DATA, &p )
                || 0 > mol_read_datum( &ctab->MolBond[i].nAtomNo2,      3,  MOL_SHORT_INT_DATA, &p )
                || 0 > mol_read_datum( &ctab->MolBond[i].cBondType,     3,  MOL_CHAR_INT_DATA,  &p )
                || 0 > mol_read_datum( &ctab->MolBond[i].cBondStereo,   3,  MOL_CHAR_INT_DATA,  &p )
#if ( MOL_QUERY == MOL_PRESENT )
                || 0 > mol_read_datum( &ctab->MolBond[i].cBondTopology, 3,  MOL_CHAR_INT_DATA,  &p ) /* ring/chain */
#else
                || 0 > mol_read_datum( NULL,                            3,  MOL_JUMP_TO_RIGHT,  &p )
#endif
#if ( MOL_REACT == MOL_PRESENT )
                || 0 > mol_read_datum( &ctab->MolBond[i].cReactingCenterStatus, 3,  MOL_CHAR_INT_DATA,  &p )
#else
                || 0 > mol_read_datum( NULL,                            3,  MOL_JUMP_TO_RIGHT,  &p )
#endif
            ){
                if ( !err ) {
                    /* can't interpret bonds block line */
                    MOLFILE_ERR_SET (err, 4, "Cannot interpret bond block line:");
                    RemoveNonPrintable( line );
                    AddMOLfileError(pStrErr, line);
                }
                if ( !strcmp( line, SDF_END_OF_DATA ) ) {
                    err = -abs(err);
                    break;
                }
            }
        }
    }
    /* err_fin: */
    return err;
}
/********** static ***************************************************************/
int read_stext_block( MOL_CTAB* ctab, FILE *inp, int err, char *pStrErr )
{
    /* just pass by all stext enties without attemp to interpret */
    char *p;
    char line[MOLFILEINPLINELEN];
    const int line_len = sizeof(line);
    S_SHORT i;

    for ( i = 0; i < 2*ctab->nNumberOfStextEntries; i++ ) {
        
        if ( NULL == ( p = fgets_up_to_lf( line, line_len, inp ) ) ){
            if ( !err ) {
                MOLFILE_ERR_FIN (err, 2, err_fin, "Cannot read STEXT block line");
            }
            break;
            /* can't read the input file line */
        }
        /*
        remove_one_lf( line );
        if ( line[MOLFILEMAXLINELEN] ){
            MOLFILE_ERR_SET (err, 2, "Warning: Too long STEXT block line");
            // too long input file line
        }
        */
    }
err_fin:
    return err;
} 
/************ static *************************************************************/
int read_properties_block( MOL_CTAB* ctab, MOL_HEADER_BLOCK *pHdr, FILE *inp, int err, char *pStrErr )
{
    enum { MULTI_LINE_MODE_NO_MODE, MULTI_LINE_MODE_ISIS_ALIAS }; 
    char *p;
    char line[MOLFILEINPLINELEN];
    const int line_len = sizeof(line);
    int   nMultiLineMode = MULTI_LINE_MODE_NO_MODE, nAtomNumber=0;
    S_SHORT i, j;
    char  charM[2];
    char  szBlank[3];
    char  szType[4];
    S_SHORT  skip_lines=0;
    S_SHORT  num_entries;
    S_SHORT  num_atoms = ctab->nNumberOfAtoms;

    int  charge_encountered  = 0;
    int  radical_encountered = 0;
    int  isotope_encountered = 0;
    /*
      if ( NULL == ctab->MolAtom ){
          err = 1;
          goto err_fin;    internal error: memory has not been allocated for MolAtom structure
      }
     */
    for ( i = 0; ctab->csCurrentCtabVersion[0]? 1 : (i < ctab->nNumberOfPropertyLines); i++ ) { /* the last line should be M END */
        /* ctab->csCurrentCtabVersion[0] == 0:
              exactly ctab->nNumberOfPropertyLines lines including M END */
        /* ctab->csCurrentCtabVersion[0] != 0:
              read until M END line was encountered */
        if ( NULL == ( p = fgets_up_to_lf( line, line_len, inp ) ) ){
            if ( !err ) {
                MOLFILE_ERR_SET (err, 2, "Cannot read properties block line");
            }
            goto err_fin;
        }
        remove_one_lf( line );
        if ( line[MOLFILEMAXLINELEN] ){
            MOLFILE_ERR_SET (err, 3, "Too long properties block line");
            continue;
        }
        if ( skip_lines > 0 ) {
            skip_lines --;
            continue;
        }
        /* alias. */
        if ( nMultiLineMode == MULTI_LINE_MODE_ISIS_ALIAS && nAtomNumber ) {
            int  len;
            nMultiLineMode = MULTI_LINE_MODE_NO_MODE;
            if ( 0 >= (len=normalize_name( p )) ) {
                nAtomNumber = 0;
                continue;
            }
            if( 0 < len && len < (int)(sizeof(ctab->MolAtom->szAtomSymbol)) ) {
                int  nCharge, nRad;
                MOL_ATOM*  MolAtom = ctab->MolAtom + nAtomNumber-1;
                /* ctab->MolAtom[nAtomNumber-1].cAtomAliasedFlag = 1; */
                /*  extract radicals & charges */
                extract_ChargeRadical( p, &nRad, &nCharge );
                /*  Aliased atom cannot have charge, radical & mass difference */
                /*  in the atom table or "M  CHG", "M  RAD", "M  ISO" */
                /* if ( nCharge ) */
                    MolAtom->cCharge = (S_CHAR)nCharge;
                /* if ( nRad ) */
                    MolAtom->cRadical = (char)nRad;   
                
                if ( 1 == len && 'D' == p[0]    ) {
                    /*  H isotope */
                    p[0] = 'H';
#ifdef INCHI_MAIN
                    MolAtom->cMassDifference=(1 + ISOTOPIC_SHIFT_FLAG);
#else
                    MolAtom->cMassDifference=1;
#endif
                } else
                if ( 1 == len && 'T' == p[0]    ) {
                    /*  H isotope */
                    p[0] = 'H';
#ifdef INCHI_MAIN
                    MolAtom->cMassDifference=(2 + ISOTOPIC_SHIFT_FLAG);
#else
                    MolAtom->cMassDifference=2;
#endif
                } else
                    MolAtom->cMassDifference=0;
                if ( strlen(p) < sizeof(ctab->MolAtom[0].szAtomSymbol) ) {
                    strcpy(MolAtom->szAtomSymbol, p);
                } else {
                    strcpy(MolAtom->szAtomSymbol, "???");
                }
                MolAtom->cAtomAliasedFlag ++;
            }
            skip_lines = 0;
            nAtomNumber = 0;
            continue;
        }
        
        if (   1 != mol_read_datum( charM,     sizeof(charM)   - 1,  MOL_STRING_DATA, &p )
            || 0 != mol_read_datum( szBlank,   sizeof(szBlank) - 1,  MOL_STRING_DATA, &p ) /* must contain 0 bytes */
            || 0 >= mol_read_datum( szType,    sizeof(szType)  - 1,  MOL_STRING_DATA, &p ) /* must contain 3 bytes */
        ) {
            if ( !strcmp( line, SDF_END_OF_DATA ) ) {
                err = err? -abs(err): -4;
                break;
            }
            continue;  /* ignore because cannot recognize */
        }
        if ( charM[0] == 'V' ){
            skip_lines = 0;   /* ISIS/Desktop Atom Value: one-line property */
            continue;
        }
        if ( charM[0] == 'G' ){
            skip_lines = 1;   /* ISIS/Desktop Group abbreviation: two-line property */
            continue;
        }
        if ( charM[0] == 'A' ) {
            if ( NULL != ctab->MolAtom &&
                 0 < ( nAtomNumber = (int)strtol(szType, NULL, 10) ) &&
                 nAtomNumber <= ctab->nNumberOfAtoms  ){
                /* Atom Alias [ISIS/Desktop] two-line property */
                nMultiLineMode = MULTI_LINE_MODE_ISIS_ALIAS;
                continue;
            } else {
                nAtomNumber = 0;
                skip_lines = 1;   
                continue;
            }
        }
        if ( charM[0] == 'S' && !strcmp( szType, "SKP" ) ){  /* skip lines */
            if ( 0 >= mol_read_datum( &skip_lines, 3, MOL_SHORT_INT_DATA, &p ) ) {
                skip_lines = 0;
            }
            continue;
        }
        if ( charM[0] != 'M' ) {/* cannot recognize a line */
            continue;
        }
        if ( !strcmp( szType, "REG" ) ) {
            int len;
            p = p + strspn( p, " " );
            len = strcspn( p, " " );
            len = inchi_min( len, MOL_MAX_VALUE_LEN );
            mol_read_datum( &pHdr->lInternalRegistryNumber, len, MOL_LONG_INT_DATA, &p );
            continue;
        }

        if ( !strcmp( szType, "END" ) ){
            if ( ctab->csCurrentCtabVersion[0] )
                break;  /* end of property lines */
            continue;
        }
        
        if ( NULL == ctab->MolAtom )
            continue; /* ignore because the user requested to bypass all this stuff */
        
        /*----------------------------------- charge: Generic */
        if ( !strcmp( szType, "CHG" ) &&
             0 < mol_read_datum( &num_entries, 3, MOL_SHORT_INT_DATA, &p ) &&
             1 <= num_entries && num_entries <= 8 ) {
            S_SHORT atoms[8];
            S_SHORT charges[8];
            if ( !charge_encountered && !radical_encountered ) {
                /* first charge or radical record clears all Atom Block */
                /* entered charge and radical data to zeroes            */
                charge_encountered = -1;
            }
            for ( j = 0; j < num_entries; j++ ) {
                if ( 0 > mol_read_datum( &atoms[j],    0, MOL_SHORT_INT_DATA, &p ) ||
                     0 > mol_read_datum( &charges[j],  0, MOL_SHORT_INT_DATA, &p ) ||
                     atoms[j]   <=  0 || atoms[j]    > num_atoms ||
                     charges[j] < -15 || charges[j]  > 15 ) {
                    goto charge_error;
                }
            }
            if ( charge_encountered == -1 ) {
                for ( j = 0; j < num_atoms; j++ ) {
                    if ( !ctab->MolAtom[j].cAtomAliasedFlag ) /* do not clear aliased atoms.*/
                        ctab->MolAtom[j].cCharge = ctab->MolAtom[j].cRadical = '\0';
                }
                charge_encountered = 1;
            }
            for ( j = 0; j < num_entries; j++ ) {
                if ( !ctab->MolAtom[atoms[j]-1].cAtomAliasedFlag ) /* do not change aliased atoms.*/
                    ctab->MolAtom[atoms[j]-1].cCharge = (S_CHAR)charges[j];
            }
            continue;
        charge_error:
            MOLFILE_ERR_SET (err, 0, "Charge not recognized:");
            RemoveNonPrintable( line );
            AddMOLfileError(pStrErr, line);
            continue; /* ignore for now */
        }
        /*-------------------------------------- radical: Generic */
        if ( !strcmp( szType, "RAD" ) &&
             0 < mol_read_datum( &num_entries, 3, MOL_SHORT_INT_DATA, &p ) &&
             1 <= num_entries && num_entries <= 8 ) {
            S_SHORT atoms[8];
            S_SHORT radicals[8];
            if ( !charge_encountered && !radical_encountered ) {
                /* first charge or radical record clears all Atom Block */
                /* entered charge and radical data to zeroes            */
                radical_encountered = -1;
            }
            for ( j = 0; j < num_entries; j++ ) {
                if ( 0 > mol_read_datum( &atoms[j],     0, MOL_SHORT_INT_DATA, &p ) ||
                     0 > mol_read_datum( &radicals[j],  0, MOL_SHORT_INT_DATA, &p ) ||
                     atoms[j]    <=  0 || atoms[j]    > num_atoms ||
                     radicals[j] <   0 || radicals[j]  > 3 ) {
                    goto radical_error;
                }
            }
            if ( radical_encountered == -1 ) {
                for ( j = 0; j < num_atoms; j++ ) {
                    if ( !ctab->MolAtom[j].cAtomAliasedFlag )  /* do not clear aliased atoms. 5-3-99 DCh */
                        ctab->MolAtom[j].cCharge = ctab->MolAtom[j].cRadical = '\0';
                }
                radical_encountered = 1;
            }
            for ( j = 0; j < num_entries; j++ ) {
                if ( !ctab->MolAtom[atoms[j]-1].cAtomAliasedFlag ) { /* do not change aliased atoms. 5-3-99 DCh */
                    ctab->MolAtom[atoms[j]-1].cRadical = (S_CHAR)radicals[j];
                }
            }
            continue;
        radical_error:
            MOLFILE_ERR_SET (err, 0, "Radical not recognized:");
            RemoveNonPrintable( line );
            AddMOLfileError(pStrErr, line);
            continue; /* ignore error for now */
        }
        /*-------------------------------------- isotope: Generic */
        if ( !strcmp( szType, "ISO" ) &&
             0 < mol_read_datum( &num_entries, 3, MOL_SHORT_INT_DATA, &p ) &&
             1 <= num_entries && num_entries <= 8 ) {
            S_SHORT atoms[8];
            S_SHORT iso_mass[8]; /*  contains istotope mass number, not difference. 7-14-00 DCh. */
            if ( !isotope_encountered  ) {
                /* first charge or radical record clears all Atom Block */
                /* entered charge and radical data to zeroes            */
                isotope_encountered = -1;
            }
            for ( j = 0; j < num_entries; j++ ) {
                if ( 0 > mol_read_datum( &atoms[j],     0, MOL_SHORT_INT_DATA, &p ) ||
                     0 > mol_read_datum( &iso_mass[j],  0, MOL_SHORT_INT_DATA, &p ) ||
                     atoms[j]    <=  0 || atoms[j]    > num_atoms
                     /*|| iso_mass[j] < -18 || iso_mass[j]  > 12*/ ) {
                    /* goto isotope_error; */
                    atoms[j] = -1; /*  flag error */
                    MOLFILE_ERR_SET (err, 0, "Isotopic data not recognized:");
                    RemoveNonPrintable( line );
                    AddMOLfileError(pStrErr, line);
                    continue; /* ignore isotopic error for now */
                }
            }
            if ( isotope_encountered == -1 ) {
                for ( j = 0; j < num_atoms; j++ ) {
                    /*if ( !ctab->MolAtom[j].cAtomAliasedFlag )*/  /* clear even aliased atoms */
                        ctab->MolAtom[j].cMassDifference = 0;
                }
                isotope_encountered = 1;
            }
            for ( j = 0; j < num_entries; j++ ) {
                if ( atoms[j] <= 0 )
                    continue; /* ignore isotopic error for now */
                if ( 1 /* !ctab->MolAtom[atoms[j]-1].cAtomAliasedFlag */) {
                    char *at =  ctab->MolAtom[atoms[j]-1].szAtomSymbol;
                    if ( at[1] || at[0] != 'D' && at[0] != 'T' ) {  /*  D & T cannot have ISO */
                        /*  need atomic weight to calculate isotope difference. 7-14-00 DCh. */
#ifdef INCHI_MAIN
                        ctab->MolAtom[atoms[j]-1].cMassDifference = iso_mass[j]; /* mass, not difference */
#else
                        int  atw, atw_diff; 
                        if ( (atw = get_atw( at )) && abs( atw_diff = (int)iso_mass[j] - atw ) < 20 ) {
                            ctab->MolAtom[atoms[j]-1].cMassDifference = (char)(atw_diff? atw_diff : ZERO_ATW_DIFF);
                        }
#endif
                    }        
                }    
            }
            continue;
        }
    }
err_fin:
    return err;
}
/************ global *************************************************************/
MOL_DATA* delete_mol_data( MOL_DATA* mol_data )
{
    if ( mol_data ) {
        if ( mol_data->ctab.MolAtom )
            inchi_free( mol_data->ctab.MolAtom );
        if ( mol_data->ctab.MolBond )
            inchi_free( mol_data->ctab.MolBond );
        if ( mol_data->ctab.szCoord )
            inchi_free( mol_data->ctab.szCoord );
        inchi_free( mol_data );
        mol_data = NULL;
    }
    return mol_data;
}
/************* global ************************************************************/
/*  Comletely ingnore STEXT block, queries, and 3D features 
 */
MOL_DATA* read_mol_file( FILE* inp, MOL_HEADER_BLOCK *OnlyHeaderBlock, MOL_CTAB *OnlyCtab,
                         int bGetOrigCoord, int *err, char *pStrErr )
{
    MOL_DATA* mol_data = NULL;
    int       ret      = 0, prev_ret, bEndOfData = 0;
    int       bReadAll = ( OnlyHeaderBlock == NULL );
    MOL_CTAB  ctab,  *pCtab = NULL;
    MOL_HEADER_BLOCK *pHdr  = NULL;
    
    *err = 0;
    if ( bReadAll ) {
        if ( NULL == ( mol_data = ( MOL_DATA* )inchi_calloc( 1, sizeof(MOL_DATA) ) ) ){
            ret = 1; /* can't allocate mol_data structure */
            AddMOLfileError( pStrErr, "Out of RAM" );
            goto err_fin;
        }
        pHdr  = &mol_data->hdr;
        pCtab = &mol_data->ctab;
    } else {
        pHdr  = OnlyHeaderBlock;
        pCtab = OnlyCtab? OnlyCtab : &ctab;
        memset( pHdr,  0, sizeof( MOL_HEADER_BLOCK ) );
        memset( pCtab, 0, sizeof( MOL_CTAB ) );
    }
    pCtab->MolBond = NULL;
    pCtab->MolAtom = NULL;
    pCtab->szCoord = NULL;
    
    if ( 0 != ( ret = mol_read_hdr(pHdr, inp, pStrErr) ) ){
        ret += 10;
        goto err_fin; /*  most probably end of file */
    }
    if ( 0 != ( ret = mol_read_counts_line( pCtab , inp, pStrErr) ) ){
        ret += 20;
        goto err_fin;
    }
    
    if ( bReadAll ) {
        if ( NULL == ( mol_data->ctab.MolAtom = (MOL_ATOM*)inchi_calloc(inchi_max(mol_data->ctab.nNumberOfAtoms,1), sizeof(MOL_ATOM)) ) ){
            ret = 2; /* can't allocate MolAtom structure */
            MOLFILE_ERR_FIN (ret, 2, err_fin, "Out of RAM");
        }
        if ( bGetOrigCoord &&
             NULL == ( mol_data->ctab.szCoord = (MOL_COORD*)inchi_calloc(inchi_max(mol_data->ctab.nNumberOfAtoms,1), sizeof(MOL_COORD)) ) ){
            ret = 2; /* can't allocate MolAtom structure */
            MOLFILE_ERR_FIN (ret, 2, err_fin, "Out of RAM");
        }
    }
    if ( 0 != ( ret = read_atom_block(pCtab, inp, ret, pStrErr) ) ){
        if ( ret < 0 ) {
            ret = -ret;
            bEndOfData = 1;
        }
        ret += 30;
        /* goto err_fin; */
    }
    
    if ( bReadAll && ret < 30 ) {
        if ( !bEndOfData && NULL == ( mol_data->ctab.MolBond = (MOL_BONDS*)inchi_calloc(inchi_max(mol_data->ctab.nNumberOfBonds,1), sizeof(MOL_BONDS)) ) ){
            ret = 3; /* can't allocate MolBond structure */
            MOLFILE_ERR_FIN (ret, 3, err_fin, "Out of RAM");
        }
    }
    prev_ret = ret;
    if ( !bEndOfData && 0 != ( ret = read_bonds_block(pCtab, inp, ret, pStrErr) ) ){
        if ( ret < 0 ) {
            ret = -ret;
            bEndOfData = 1;
        }
        ret = prev_ret? prev_ret : ret + 40;
    }
    prev_ret = ret;
    if ( !bEndOfData && 0 != ( ret = read_stext_block(pCtab, inp, ret, pStrErr) ) ){
        ret = prev_ret? prev_ret : ret + 50;
    }
    prev_ret = ret;
    if ( !bEndOfData && 0 != ( ret = read_properties_block(pCtab, pHdr, inp, ret, pStrErr) ) ){
        if ( ret < 0 ) {
            ret = -ret;
            bEndOfData = 1;
        }
        ret = prev_ret? prev_ret : ret + 60;
    }
    
err_fin:
    *err = bEndOfData? -ret : ret;
    if ( bReadAll ) {
        if ( ret )
            mol_data = delete_mol_data( mol_data ); /* delete all results */
        return mol_data;
    } else {
        if ( ret )
            return NULL;
        else
            return (MOL_DATA*)OnlyHeaderBlock;
    }
}

/******************************************************************/
static const char sdf_data_hdr_name[] = "NAME";
static const char sdf_data_hdr_comm[] = "COMMENT";
enum { SDF_START, SDF_DATA_HEADER, SDF_DATA_HEADER_NAME
     , SDF_DATA_HEADER_COMMENT, SDF_DATA_HEADER_CAS
     , SDF_DATA_HEADER_USER, SDF_DATA_LINE
     , SDF_END_OF_DATA_ITEM, SDF_EMPTY_LINE, SDF_END_OF_DATA_BLOCK };
/********** static ********************************************************/
long extract_cas_rn( char *line )
{
    int i, j;
    i = line[0] == '-'? 1 : 0;
    for ( j = i; line[i]; i ++ ) {
        if ( isdigit( UCINT line[i] ) ) {
            line[j++] = line[i];
        } else
        if ( line[i] != '-' ) {
            break;
        }
    }
    line[j] = '\0';
    return strtol( line, NULL, 10 );
}
/********** static ********************************************************/
int identify_sdf_label( char* inp_line, const char *pSdfLabel )
{
    char line[MOLFILEMAXLINELEN];
    char *p, *q;
    int  i, j, len;
    if ( (p = strchr( inp_line, '<' )) &&
         (q = strchr( p,        '>' )) &&
         (len = q-p-1) > 0 && len < (int)sizeof(line) ) {
        memcpy( line, p+1, len );
        line[len] = '\0';
        for ( i = 0; isspace( UCINT line[i] ); i ++ )
            ;
        for ( j = len-1; j >= i && isspace( UCINT line[i] ); j -- )
            ;
        len = j-i+1;
        p = line+i;
        if ( pSdfLabel && pSdfLabel[0] && len == (int)strlen(pSdfLabel) && !memicmp( p, pSdfLabel, len ) )
            return SDF_DATA_HEADER_USER;
        if ( len == sizeof(sdf_data_hdr_name)-1 && !memicmp( p, sdf_data_hdr_name, len ) )
            return SDF_DATA_HEADER_NAME;
        if ( len == sizeof(sdf_data_hdr_comm)-1 && !memicmp( p, sdf_data_hdr_comm, len ) )
            return SDF_DATA_HEADER_COMMENT;
        if ( !memicmp( p, "CAS", 3 ) )
            return SDF_DATA_HEADER_CAS;
    }
    return SDF_DATA_HEADER;
}
/************* global *****************************************************/
int bypass_sdf_data_items( FILE* inp, long *cas_reg_no, char* comment,
                           int lcomment, char *name, int lname, int prev_err,
                           const char *pSdfLabel, char *pSdfValue, char *pStrErr )
{
    char line[MOLFILEINPLINELEN];
    const int line_len = sizeof(line);
    int   err = 0;
    int   current_state = SDF_START;
    int   n_blank_lines = 0;
    int   n_lines       = 0;
    char* p = NULL;
    int   bNeedsName   = name && lname > 0 && !name[0];
    int   bNeedsComm   = comment && lcomment > 0 && !comment[0];
    int   bNeedsUser   = pSdfLabel && pSdfLabel[0] && pSdfValue;
    int   bNeedsCASrn  = 0;
    int   bCASrnIsUser = 0;

    if ( cas_reg_no != NULL ) {
        bNeedsCASrn = 1;
        *cas_reg_no = 0;
        bCASrnIsUser = (bNeedsUser && !memicmp(pSdfLabel,"CAS", 3));
    }

    while ( err           == 0                    &&
            current_state !=SDF_END_OF_DATA_BLOCK &&
            NULL != ( p = fgets_up_to_lf( line, line_len, inp ) ) ) {
        
        if ( !n_lines && !memcmp(line, "M  END", 6) ) {
            continue; /*  allow subtle errors */
        }
        n_lines++;
        
        remove_trailing_spaces( line );
        if ( line[MOLFILEMAXLINELEN] ){
            if ( current_state != SDF_DATA_HEADER &&
                 current_state != SDF_DATA_LINE   &&
                 current_state != SDF_DATA_HEADER_NAME &&
                 current_state != SDF_DATA_HEADER_USER &&
                 current_state != SDF_DATA_HEADER_COMMENT ) {
                line[MOLFILEMAXLINELEN] = '\0';
                if ( !prev_err ) {
                    MOLFILE_ERR_SET (err, 0, "Too long SData line truncated");
                }
            } else {
                /* allow long lines in SDF data. 9-29-00 DCh */
                line[MOLFILEMAXLINELEN] = '\0';
            }
        }
        
        n_blank_lines += ( *line == '\0' );
        
        switch( current_state ) {
        
        case SDF_START:
        case SDF_END_OF_DATA_ITEM:
        case SDF_EMPTY_LINE:              /* Added 9-25-97 DCh */
        
            if ( 0 == strcmp( line, SDF_END_OF_DATA ) ) {
                current_state = SDF_END_OF_DATA_BLOCK;
            }
            else    
            if ( '>' == *line ) {
                current_state = ( bNeedsName || bNeedsComm || bNeedsCASrn || bNeedsUser )? identify_sdf_label(line, pSdfLabel) : SDF_DATA_HEADER;
            }else
            if ( *line == '\0' ) { /* Added 9-25-97 DCh */
                /* Relax the strictness: Allow more than 1 empty line. */
                current_state=SDF_EMPTY_LINE;
            } else
            if ( !prev_err ) {
                MOLFILE_ERR_SET (err, 3, "Unexpected SData header line:");
                RemoveNonPrintable( line );
                AddMOLfileError(pStrErr, line);
                /* unexpected contents of data header line */
            } else {
                err = 3;
            }
            break;
        
        case SDF_DATA_HEADER_NAME:
             if ( bNeedsName && 0 < normalize_name( line ) ) {
                bNeedsName = 0;
                mystrncpy( name, line, lname );
             }
             goto got_data_line;
                
        case SDF_DATA_HEADER_COMMENT:
             if ( bNeedsComm && 0 < normalize_name( line ) ) {
                bNeedsComm = 0;
                mystrncpy( comment, line, lcomment );
             }
             goto got_data_line;

        case SDF_DATA_HEADER_USER:
             if ( bNeedsUser && 0 < normalize_name( line ) ) {
                 bNeedsUser = 0;
                 mystrncpy( pSdfValue, line, MAX_SDF_VALUE+1 );
                 if ( bCASrnIsUser && bNeedsCASrn ) {
                     *cas_reg_no = extract_cas_rn( line );
                     bNeedsCASrn = (0 == *cas_reg_no);
                 }
             }
             goto got_data_line;

        case SDF_DATA_HEADER_CAS:
             if ( bNeedsCASrn && 0 < normalize_name( line ) ) {
                 *cas_reg_no = extract_cas_rn( line );
                 bNeedsCASrn = (0 == *cas_reg_no);
             }
             goto got_data_line;

        case SDF_DATA_HEADER:
        case SDF_DATA_LINE:
got_data_line:
            current_state = *line? SDF_DATA_LINE : SDF_END_OF_DATA_ITEM;
            break;
        
        }
    }    
    if ( 0 == err && SDF_END_OF_DATA_BLOCK != current_state && NULL == p )
        ; /* err = 4; */ /* unexpected end of file: missing $$$$ */
    else
    if (err && ( n_blank_lines == n_lines && *line == '\0' ) )
        err = 5; /* empty lines -- do not know when this can happen */

    if ( err && err != 5 && current_state != SDF_END_OF_DATA_BLOCK && p ) {
        /*  bypass up to $$$$ */
        while ( ( p = fgets_up_to_lf( line, line_len, inp ) ) && memcmp( line, SDF_END_OF_DATA, 4 ) )
            ;
        if ( p ) {
            err = 9; /*  bypassed to $$$$; non-fatal */
            AddMOLfileError(pStrErr, "Bypassing to next structure");
        }

    }

    return err;
}
/**************** global **************************************************/
MOL_DATA* read_sdfile_segment(FILE* inp, MOL_HEADER_BLOCK *OnlyHeaderBlock, MOL_CTAB *OnlyCtab,
                              int bGetOrigCoord,
                              char *pname, int lname,
                              long *Id, const char *pSdfLabel, char *pSdfValue,
                              int *err, char *pStrErr )
{
    MOL_DATA* mol_data = read_mol_file( inp, OnlyHeaderBlock, OnlyCtab, bGetOrigCoord, err, pStrErr );
    int       err_bypass_sdf = 0;

    if ( pname && lname ) {
        pname[0] = '\0';
    }
    if ( Id ) {
        *Id = 0L;  /* ignore for now */
    }
    /* if ( mol_data && !*err ) { */
    if ( *err < 0 ) {
        *err = -*err; /* end of data encountered */
    } else {
        err_bypass_sdf = bypass_sdf_data_items( inp, Id, NULL, 0, pname, lname, *err, pSdfLabel, pSdfValue, pStrErr );
        if ( err_bypass_sdf ) {
            *err = err_bypass_sdf; /* important to continue to the next good structure */
        }
    }
    /* } */
    return mol_data;
}
/******************* global *********************************************************/
int CopyMOLfile(FILE *inp_file, long fPtrStart, long fPtrEnd, FILE *prb_file, long lNumb)
{
    char line[MOLFILEINPLINELEN], *p;
    long fPtr;
    int  ret = 1;
    char szNumber[32];

    if ( inp_file && prb_file && fPtrStart >= 0L &&
         fPtrEnd > fPtrStart &&
         0 == fseek( inp_file, fPtrStart, SEEK_SET ) ) {

        while ( fPtrEnd > (fPtr = ftell(inp_file)) && fPtr >= 0L &&
                fgets_up_to_lf( line, sizeof(line)-1, inp_file ) ) {
            line[sizeof(line)-1] = '\0'; /*  unnecessary extra precaution */
            if ( fPtr == fPtrStart && lNumb ) {
                int len;
                LtrimRtrim( line, &len );
                len = sprintf( szNumber, "#%ld%s", lNumb, len?"/":"" );
                mystrncpy( line+len, line, sizeof(line)-len-1 );
                memcpy( line, szNumber, len );
            }
            if ( !strchr(line, '\n') ) {
                p = line+strlen(line);
                p[0] = '\n';
                p[1] = '\0';
            }
            fputs( line, prb_file );
        }
        ret = fseek( inp_file, fPtrEnd, SEEK_SET );
    }
    return ret;
}

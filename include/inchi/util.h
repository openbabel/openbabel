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


#ifndef __UTIL_H__
#define __UTIL_H__

#include "inpdef.h"

/* BILLY 8/6/04 */
#ifndef INCHI_ALL_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif

int get_atw(const char *elname);
int get_atw_from_elnum( int nAtNum );
int get_num_H (const char* elname, int inp_num_H, S_CHAR num_iso_H[], int charge, int radical,
              int chem_bonds_valence, int atom_input_valence, int bAliased, int bDoNotAddH, int bHasMetalNeighbor );
int extract_ChargeRadical( char *elname, int *pnRadical, int *pnCharge );
int extract_H_atoms( char *elname, S_CHAR num_iso_H[] );
int normalize_name( char* name );

int mystrncpy(char *target,const char *source,unsigned maxlen);
char* LtrimRtrim( char *p, int* nLen );
char* fgets_up_to_lf( char* line, int line_len, FILE* inp );
void remove_trailing_spaces( char* p );
void remove_one_lf( char* p);
void mystrrev( char *p );
#ifndef INCHI_LIBRARY
int my_fgets( char *szLine, int len, FILE *f, int *bTooLongLine );
int my_fgetsTab( char *szLine, int len, FILE *f, int *bTooLongLine );
int my_fgetsTab1( char *szLine, int len, FILE *f, int *bTooLongLine );
#endif
int my_fgetsUpToLfOrTab( char *szLine, int len, FILE *f );

#define ALPHA_BASE  27
long inchi_strtol( const char *str, const char **p, int base);
double inchi_strtod( const char *str, const char **p );

AT_NUMB *is_in_the_list( AT_NUMB *pathAtom, AT_NUMB nNextAtom, int nPathLen );
int get_periodic_table_number( const char* elname );
int is_el_a_metal( int nPeriodicNum );
int get_el_valence( int nPeriodicNum, int charge, int val_num );
int get_unusual_el_valence( int nPeriodicNum, int charge, int radical, int bonds_valence, int num_H, int num_bonds );
int detect_unusual_el_valence( int nPeriodicNum, int charge, int radical, int bonds_valence, int num_H, int num_bonds );
int needed_unusual_el_valence( int nPeriodicNum, int charge, int radical, int bonds_valence,
                               int actual_bonds_val, int num_H, int num_bonds );
int get_el_type( int nPeriodicNum );
int get_el_number( const char* elname );
int do_not_add_H( int nPeriodicNum );
int GetElementFormulaFromAtNum(int nAtNum, char *szElement );
int MakeRemovedProtonsString( int nNumRemovedProtons, NUM_H *nNumExchgIsotopicH, NUM_H *nNumRemovedProtonsIsotopic,
                              int bIsotopic, char *szRemovedProtons, int *num_removed_iso_H );

/* ion pairs and fixing bonds */
int num_of_H( inp_ATOM *at, int iat );
int has_other_ion_neigh( inp_ATOM *at, int iat, int iat_ion_neigh, const char *el, int el_len );
int has_other_ion_in_sphere_2(inp_ATOM *at, int iat, int iat_ion_neigh, const char *el, int el_len );
int nNoMetalNumBonds( inp_ATOM *at, int at_no );
int nNoMetalBondsValence( inp_ATOM *at, int at_no );
int nNoMetalNeighIndex( inp_ATOM *at, int at_no );
int nNoMetalOtherNeighIndex( inp_ATOM *at, int at_no, int cur_neigh );
int nNoMetalOtherNeighIndex2( inp_ATOM *at, int at_no, int cur_neigh, int cur_neigh2 );

/* mol2atom.c */
int nBondsValToMetal( inp_ATOM* at, int iat );

/* ichi_bns.c */
int nBondsValenceInpAt( const inp_ATOM *at, int *nNumAltBonds, int *nNumWrongBonds );

int bHeteroAtomMayHaveXchgIsoH( inp_ATOM *atom, int iat );

/* IChICan2.c */
int SetBitFree( void );

void WriteCoord( char *str, double x );


extern int ERR_ELEM;
extern int nElDataLen;

/* BILLY 8/6/04 */
#ifndef INCHI_ALL_CPP
#ifdef __cplusplus
}
#endif
#endif

#endif /* __UTIL_H__*/


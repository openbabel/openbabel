/*
 * International Union of Pure and Applied Chemistry (IUPAC)
 * International Chemical Identifier (InChI)
 * Version 1
 * Software version 1.01
 * July 21, 2006
 * Developed at NIST
 */

#ifndef ___INCHISIZE_H__
#define ___INCHISIZE_H__

typedef unsigned short AT_NUMB;
typedef unsigned short AT_RANK;
#define AT_RANK_MASK   ((AT_RANK)~0)

typedef signed short NUM_H;
#define MAX_ATOMS  1024


#define CHAR_MASK  0xFF


typedef AT_RANK  *pAT_RANK;
typedef pAT_RANK *ppAT_RANK;

typedef unsigned long INCHI_MODE;

#define LEN_COORD 10
#define NUM_COORD 3
typedef char MOL_COORD[LEN_COORD*NUM_COORD + NUM_COORD-1]; /*copied 30 bytes from MOLfile */


#endif /* ___INCHISIZE_H__ */


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


#ifndef __INCHISTER_H__
#define __INCHISTER_H__

#ifndef INCHI_ALL_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif
int bCanAtomBeAStereoCenter( char *elname, S_CHAR charge, S_CHAR radical );
int bCanInpAtomBeAStereoCenter( inp_ATOM *at, int cur_at, int bPointedEdgeStereo );
int bCanAtomHaveAStereoBond( char *elname, S_CHAR charge, S_CHAR radical );
int bCanAtomBeTerminalAllene( char *elname, S_CHAR charge, S_CHAR radical );
int bCanAtomBeMiddleAllene( char *elname, S_CHAR charge, S_CHAR radical );
int bAtomHasValence3( char *elname, S_CHAR charge, S_CHAR radical );
int set_stereo_parity( inp_ATOM* at, sp_ATOM* at_output, int num_at, int num_removed_H,
                       int *nMaxNumStereoAtoms, int *nMaxNumStereoBonds, INCHI_MODE nMode,
                       int bPointedEdgeStereo, int vABParityUnknown );
int get_opposite_sb_atom( inp_ATOM *at, int cur_atom, int icur2nxt,
                          int *pnxt_atom, int *pinxt2cur, int *pinxt_sb_parity_ord );

#define PES_BIT_POINT_EDGE_STEREO    1
#define PES_BIT_PHOSPHINE_STEREO     2
#define PES_BIT_ARSINE_STEREO        4
#define PES_BIT_FIX_SP3_BUG          8

#ifndef INCHI_ALL_CPP
#ifdef __cplusplus
}
#endif
#endif

#endif  /* __INCHISTER_H__ */

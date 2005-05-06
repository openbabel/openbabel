/*
 * International Union of Pure and Applied Chemistry (IUPAC)
 * International Chemical Identifier (InChI)
 * Version 1
 * Software version 1.00
 * April 13, 2005
 * Developed at NIST
 */

#ifndef __INCHISTER_H__
#define __INCHISTER_H__

#ifndef INCHI_ALL_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif
int bCanAtomBeAStereoCenter( char *elname, S_CHAR charge, S_CHAR radical );
int bCanInpAtomBeAStereoCenter( inp_ATOM *at, int cur_at );
int bCanAtomHaveAStereoBond( char *elname, S_CHAR charge, S_CHAR radical );
int bCanAtomBeTerminalAllene( char *elname, S_CHAR charge, S_CHAR radical );
int bAtomHasValence3( char *elname, S_CHAR charge, S_CHAR radical );
int set_stereo_parity( inp_ATOM* at, sp_ATOM* at_output, int num_at, int num_removed_H,
                       int *nMaxNumStereoAtoms, int *nMaxNumStereoBonds, INCHI_MODE nMode,
                       int bPointedEdgeStereo );
int get_opposite_sb_atom( inp_ATOM *at, int cur_atom, int icur2nxt,
                          int *pnxt_atom, int *pinxt2cur, int *pinxt_sb_parity_ord );

#ifndef INCHI_ALL_CPP
#ifdef __cplusplus
}
#endif
#endif

#endif  /* __INCHISTER_H__ */

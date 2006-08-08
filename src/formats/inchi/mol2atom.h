/*
 * International Union of Pure and Applied Chemistry (IUPAC)
 * International Chemical Identifier (InChI)
 * Version 1
 * Software version 1.01
 * July 21, 2006
 * Developed at NIST
 */

#ifndef __MOL2ATOM_H__
#define __MOL2ATOM_H__

#include "readmol.h"

#ifndef INCHI_ALL_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif

void calculate_valences (MOL_DATA* mol_data, inp_ATOM* at, int *num_atoms, int bDoNotAddH, int *err, char *pStrErr);
/* void WriteCoord( char *str, double x );*/

#ifndef INCHI_ALL_CPP
#ifdef __cplusplus
}
#endif
#endif


#endif  /* __MOL2ATOM_H__ */

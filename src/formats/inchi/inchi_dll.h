/*
 * International Union of Pure and Applied Chemistry (IUPAC)
 * International Chemical Identifier (InChI)
 * Version 1
 * Software version 1.00
 * April 13, 2005
 * Developed at NIST
 */

#ifndef __INCHI_DLL_H__
#define __INCHI_DLL_H__


#ifndef INCHI_ALL_CPP
#ifdef  __cplusplus
extern "C" {
#endif
#endif


int ExtractOneStructure( STRUCT_DATA *sd, INPUT_PARMS *ip, char *szTitle,
         inchi_Input *inp, INCHI_FILE *log_file, INCHI_FILE *output_file, INCHI_FILE *prb_file,
         ORIG_ATOM_DATA *orig_inp_data, int *num_inp, char *pStr, int nStrLen );



#ifndef INCHI_ALL_CPP
#ifdef  __cplusplus
}
#endif
#endif


#endif /* __INCHI_DLL_H__ */

/*
 * International Chemical Identifier (InChI)
 * Version 1
 * Software version 1.07
 * April 30, 2024
 *
 * MIT License
 *
 * Copyright (c) 2024 IUPAC and InChI Trust
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
*
* The InChI library and programs are free software developed under the
 * auspices of the International Union of Pure and Applied Chemistry (IUPAC).
 * Originally developed at NIST.
 * Modifications and additions by IUPAC and the InChI Trust.
 * Some portions of code were developed/changed by external contributors
 * (either contractor or volunteer) which are listed in the file
 * 'External-contributors' included in this distribution.
 *
 * info@inchi-trust.org
 *
*/

#include "mode.h" /* djb-rwth: necessary header file */
#include "inchi_api.h" /* djb-rwth: necessary header file */
#include "inpdef.h" /* djb-rwth: necessary header file */

#ifndef _READINCH_H_
#define _READINCH_H_

#define INPUT_FILE          INCHI_IOSTREAM

#define MIN_BOND_LENGTH   (1.0e-6)
#define INCHI_LINE_LEN   262144 /*32767*/ /*512*/ /*1024*/ /*256*/
#define INCHI_LINE_ADD   32767 /*384*/  /*128*/  /*64*/


/* Convenience storage for InChI read control data */
typedef struct ReadINCHI_CtlData
{
    unsigned long ulongID;
    int bTooLongLine;
    int bHeaderRead;
    int bErrorMsg;
    int bRestoreInfo;
}
ReadINCHI_CtlData;


/*
Note:
(INCHI_LINE_LEN - INCHI_LINE_ADD) > (length of the longest item: szCoord) = 33
*/

char *FindToken( INCHI_IOSTREAM *inp_molfile,
                 int *bTooLongLine,
                 const char *sToken,
                 int lToken,
                 char *szLine,
                 int nLenLine,
                 char *p,
                 int *res );
char *LoadLine( INPUT_FILE *inp_molfile,
               int *bTooLongLine,
               int *bItemIsOver,
               char **s,
               char *szLine,
               int nLenLine,
               int nMinLen2Load,
               char *p,
               int *res );
inchi_Stereo0D *CreateInchi_Stereo0D( int num_stereo0D );
void FreeInchi_Stereo0D( inchi_Stereo0D **stereo0D );
int Extract0DParities( inp_ATOM *at,
                       int nNumAtoms,
                       inchi_Stereo0D *stereo0D,
                       int num_stereo0D,
                       char *pStrErr,
                       int *err,
                       int vABParityUnknown );
int InchiToInpAtom( INCHI_IOSTREAM *inp_molfile,
                     MOL_COORD **szCoord,
                     int bDoNotAddH,
                     int vABParityUnknown,
                     INPUT_TYPE nInputType,
                     inp_ATOM **at,
                     int max_num_at,
                     int *num_dimensions,
                     int *num_bonds,
                     char *pSdfLabel,
                     char *pSdfValue,
                     unsigned long *Id,
                     INCHI_MODE *pInpAtomFlags,
                     int *err,
                     char *pStrErr );


#endif    /* _READINCH_H_ */

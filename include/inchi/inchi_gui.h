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


#ifndef _INCHI_GUI_H_
#define _INCHI_GUI_H_


#include "strutil.h"
#include "ichicomn.h"


#ifndef COMPILE_ANSI_ONLY

struct tagCANON_GLOBALS;


int DisplayStructure( struct tagCANON_GLOBALS *pCG,
                      inp_ATOM *at,
                      int num_at,
                      OAD_Polymer *polymer,
                      int num_removed_H,
                      int bAdd_DT_to_num_H,
                      int nNumRemovedProtons,
                      NUM_H nNumRemovedProtonsIsotopic[],
                      int bIsotopic,
                      int j /*bTautomeric*/,
                      INChI **cur_INChI,
                      INChI_Aux **cur_INChI_Aux,
                      int bAbcNumbers,
                      DRAW_PARMS *dp,
                      INCHI_MODE nMode,
                      char *szTitle );

int DisplayCompositeStructure( struct tagCANON_GLOBALS *pCG,
                               COMP_ATOM_DATA *composite_norm_data,
                               OAD_Polymer *polymer,
                               int bIsotopic,
                               int bTautomeric,
                               PINChI2 *pINChI2,
                               PINChI_Aux2 *pINChI_Aux2,
                               int bAbcNumbers,
                               DRAW_PARMS *dp,
                               INCHI_MODE nMode,
                               char *szTitle );


int DisplayTheWholeStructure( struct tagCANON_GLOBALS *pCG,
                              struct tagINCHI_CLOCK *ic,
                              struct tagStructData *sd,
                              INPUT_PARMS *ip,
                              char *szTitle,
                              INCHI_IOSTREAM *inp_file,
                              INCHI_IOSTREAM *log_file,
                              ORIG_ATOM_DATA *orig_inp_data,
                              long num_inp,
                              int iINChI,
                              int bShowStruct,
                              int bINCHI_LIB_Flag );

int DisplayTheWholeCompositeStructure( struct tagCANON_GLOBALS *pCG,
                                       struct tagINCHI_CLOCK *ic,
                                       INPUT_PARMS *ip,
                                       struct tagStructData *sd,
                                       long num_inp,
                                       int iINChI,
                                       PINChI2 *pINChI2,
                                       PINChI_Aux2 *pINChI_Aux2,
                                       ORIG_ATOM_DATA *orig_inp_data,
                                       ORIG_ATOM_DATA *prep_inp_data,
                                       COMP_ATOM_DATA composite_norm_data[TAUT_NUM + 1] );


void FillTableParms( SET_DRAW_PARMS *sdp,
                     INChI **cur_INChI,
                     INChI_Aux **cur_INChI_Aux,
                     INCHI_MODE nMode,
                     int bShowIsotopic,
                     int bShowTaut );

void FillCompositeTableParms( SET_DRAW_PARMS *sdp,
                              AT_NUMB StereoFlags,
                              INCHI_MODE nMode,
                              int bShowIsotopic,
                              int bShowTaut );

#endif

#endif /* _INCHI_GUI_H_ */

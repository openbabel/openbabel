/*
 * International Union of Pure and Applied Chemistry (IUPAC)
 * International Chemical Identifier (InChI)
 * Version 1
 * Software version 1.01
 * July 21, 2006
 * Developed at NIST
 */

#ifndef __DISP_STRU_H__
#define __DISP_STRU_H__
#include <windows.h>

/* local types */
/*****************************************************/
typedef struct tagInternalDrawParms {
    double xmin, xmax, ymin, ymax;
    int    max_label_width_char;
    int    max_left_label_width_pix;
    int    max_right_label_width_pix;
    int bInit;
} INT_DRAW_PARMS; /* internal: saved for redisplaying one structure */


/*****************************************************
 *      Window data
 */
typedef struct tagWindowData {

    inp_ATOM  *at0;    /* [MAX_ATOMS]; */
    inp_ATOM  *at1;    /* [MAX_ATOMS]; */
    INF_ATOM_DATA inf_at_data;
    /*inf_ATOM  *inf_at;*/ /* [MAX_ATOMS]; */
    int       num_at;
    int       bOrigAtom;
    int       bHighlight;
    int       bEsc;
    int       bUserIntervened;
    UINT      nTimerId;

    unsigned long  ulDisplTime;
    int            nFontSize;
    RECT           rc;        /* window rectangle size for saving */
    INT_DRAW_PARMS idp;       /* structure geom. parameters for redrawing */
    TBL_DRAW_PARMS tdp;       /* table data for displaying */
    char      *szTitle;       /* for INCHI_LIB printing */

    /* component equivalence info */
    AT_NUMB *nEquLabels;    /* num_at elements or NULL */
    AT_NUMB  nNumEquSets;   /* number of equivalent sets or 0 */
    AT_NUMB  nCurEquLabel;  /* in range 0..nNumEquSets; 0=>do not display equivalent components */

    AT_NUMB  nNewEquLabel;  /* non-zero only if DISPLAY_EQU_COMPONENTS==1 */

} MY_WINDOW_DATA;

#ifndef INCHI_ALL_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif

void FreeWinData( MY_WINDOW_DATA* pWinData );
int CreateInputStructPicture( HDC hDC, MY_WINDOW_DATA *pWinData, RECT *rc, int bPrint, AT_NUMB nNewEquLabel );

#ifndef INCHI_ALL_CPP
#ifdef __cplusplus
}
#endif
#endif


#endif

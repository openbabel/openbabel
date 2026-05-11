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


#ifndef _ICCHIERR_H_
#define _ICCHIERR_H_

#define _IS_OKAY    0
#define _IS_WARNING 1
#define _IS_ERROR   2    /* Microsoft defined its own IS_ERROR() macro */
#define _IS_FATAL   3
#define _IS_UNKNOWN 4    /* unknown error: used in INChI DLL only */
#define _IS_EOF    -1    /* end of file */
#define _IS_SKIP   -2

#define CT_ERR_FIRST         (-30000)
#define CT_OVERFLOW          (CT_ERR_FIRST- 0)  /*(-30000) */
#define CT_LEN_MISMATCH      (CT_ERR_FIRST- 1)  /*(-30001) */
#define CT_OUT_OF_RAM        (CT_ERR_FIRST- 2)  /*(-30002) */
#define CT_RANKING_ERR       (CT_ERR_FIRST- 3)  /*(-30003) */
#define CT_ISOCOUNT_ERR      (CT_ERR_FIRST- 4)  /*(-30004) */
#define CT_TAUCOUNT_ERR      (CT_ERR_FIRST- 5)  /*(-30005) */
#define CT_ISOTAUCOUNT_ERR   (CT_ERR_FIRST- 6)  /*(-30006) */
#define CT_MAPCOUNT_ERR      (CT_ERR_FIRST- 7)  /*(-30007) */
#define CT_TIMEOUT_ERR       (CT_ERR_FIRST- 8)  /*(-30008) */
#define CT_ISO_H_ERR         (CT_ERR_FIRST- 9)  /*(-30009) */
#define CT_STEREOCOUNT_ERR   (CT_ERR_FIRST-10)  /*(-30010) */
#define CT_ATOMCOUNT_ERR     (CT_ERR_FIRST-11)  /*(-30011) */
#define CT_STEREOBOND_ERROR  (CT_ERR_FIRST-12)  /*(-30012) */
#define CT_USER_QUIT_ERR     (CT_ERR_FIRST-13)  /*(-30013) */
#define CT_REMOVE_STEREO_ERR (CT_ERR_FIRST-14)  /*(-30014) */
#define CT_CALC_STEREO_ERR   (CT_ERR_FIRST-15)  /*(-30015) */
#define CT_CANON_ERR         (CT_ERR_FIRST-16)  /*(-30016) */
#define CT_STEREO_CANON_ERR  (CT_ERR_FIRST-17)  /*(-30017) */
#define CT_WRONG_FORMULA     (CT_ERR_FIRST-18)  /*(-30017) */
#define CT_UNKNOWN_ERR       (CT_ERR_FIRST-19)  /*(-30019) */

#define CT_ERR_MIN CT_UNKNOWN_ERR
#define CT_ERR_MAX CT_ERR_FIRST

#define CHECK_OVERFLOW(Len, Maxlen) ( (Len) >= (Maxlen) )
#define RETURNED_ERROR(nVal) (CT_ERR_MIN<=(nVal) && (nVal)<=CT_ERR_MAX)


#define BNS_ERR            -9999
#define BNS_WRONG_PARMS    (BNS_ERR +  0) /*(-9999)*/
#define BNS_OUT_OF_RAM     (BNS_ERR +  1) /*(-9998)*/
#define BNS_PROGRAM_ERR    (BNS_ERR +  2) /*(-9997)*/
#define BNS_ALTPATH_OVFL   (BNS_ERR +  3) /*(-9996)*/
#define BNS_BOND_ERR       (BNS_ERR +  4) /*(-9995)*/
#define BNS_VERT_NUM_ERR   (BNS_ERR +  5) /*(-9994)*/
#define BNS_VERT_EDGE_OVFL (BNS_ERR +  6) /*(-9993)*/
#define BNS_SET_ALTP_ERR   (BNS_ERR +  7) /*(-9992)*/
#define BNS_CPOINT_ERR     (BNS_ERR +  8) /*(-9991)*/
#define BNS_CANT_SET_BOND  (BNS_ERR +  9) /*(-9990)*/
#define BNS_CAP_FLOW_ERR   (BNS_ERR + 10) /*(-9989)*/
#define BNS_RADICAL_ERR    (BNS_ERR + 11) /*(-9988)*/
#define BNS_REINIT_ERR     (BNS_ERR + 12) /*(-9987)*/
#define BNS_ALTBOND_ERR    (BNS_ERR + 13) /*(-9986)*/

#define BNS_TIMEOUT        (BNS_ERR + 14) /*(-9985)*/ /* v. 1.05 */

#define BNS_MAX_ERR_VALUE  (BNS_ERR + 19) /*(-9980)*/

#define IS_BNS_ERROR(X) (BNS_ERR <= (X) && (X) <= BNS_MAX_ERR_VALUE)


#define INCHI_INP_ERROR_ERR   40
#define INCHI_INP_ERROR_RET  (-1)

#define INCHI_INP_FATAL_ERR    1
#define INCHI_INP_FATAL_RET    0

#define INCHI_INP_EOF_ERR     11
#define INCHI_INP_EOF_RET      0

#ifndef COMPILE_ALL_CPP
#ifdef  __cplusplus
extern "C" {
#endif
#endif

    extern int( *UserAction )( void ); /* callback */
    extern int( *ConsoleQuit )( void ); /* Console user issued CTRL+C etc. */

#ifndef COMPILE_ALL_CPP
#ifdef  __cplusplus
}
#endif
#endif

#define LOG_MASK_WARN  1
#define LOG_MASK_ERR   2
#define LOG_MASK_FATAL 4

#define LOG_MASK_ALL     (LOG_MASK_WARN | LOG_MASK_ERR | LOG_MASK_FATAL)
#define LOG_MASK_NO_WARN (LOG_MASK_ERR | LOG_MASK_FATAL)

#ifdef TARGET_LIB_FOR_WINCHI
#include <stdarg.h>

#ifndef COMPILE_ALL_CPP
#ifdef  __cplusplus
extern "C" {
#endif
#endif

    extern void( *FWPRINT ) ( const char * format, va_list argptr );
    extern void( *FWPUSH ) ( const char *s );
    extern void( *DRAWDATA ) ( struct DrawData * pDrawData );
    extern int( *DRAWDATA_EXISTS ) ( int nComponent, int nType, int bReconnected );
    extern struct DrawData * ( *GET_DRAWDATA ) ( int nComponent, int nType, int bReconnected );

#ifndef COMPILE_ALL_CPP
#ifdef  __cplusplus
}
#endif
#endif

#endif

#define USER_ACTION_QUIT   1

#define STR_ERR_LEN 256

const char *ErrMsg( int nErrorCode );

int AddErrorMessage( char *pStrErr, const char *szMsg );

#define WarningMessage AddErrorMessage


#define TREAT_ERR_AND_FIN(err, new_err, err_fin, msg) \
        if ( !(err) && (new_err) ) { (err) = (new_err);} AddErrorMessage(pStrErr, (msg)); goto err_fin
#define TREAT_ERR(err, new_err, msg) \
        if ( !(err) && (new_err) ) { (err) = (new_err);} AddErrorMessage(pStrErr, (msg))


#endif /* _ICCHIERR_H_ */

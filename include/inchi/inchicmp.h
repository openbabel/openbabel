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


#ifndef _INCHICMP_H_
#define _INCHICMP_H_


typedef enum tagInchiCompareDiffBits {
    INCHIDIFF_ZERO = 0x00000000,
    INCHIDIFF_PROBLEM = 0x00000001, /* severe: at least one InChI does not exist */
    INCHIDIFF_NUM_AT = 0x00000001, /* severe: different number of atoms in the skeleton */
    INCHIDIFF_ATOMS = 0x00000001, /* severe: diiferent types of skeleton atoms */
    INCHIDIFF_NUM_EL = 0x00000001, /* severe: formulas differ in another element */
    INCHIDIFF_CON_LEN = 0x00000001, /* severe: different connection table lengths */
    INCHIDIFF_CON_TBL = 0x00000001, /* severe: different connection tables */
    INCHIDIFF_POSITION_H = 0x00000002, /* difference in non-taut (Mobile-H) or all H (Fixed-H) location/number */
    INCHIDIFF_MORE_FH = 0x00000004, /* extra fixed H */
    INCHIDIFF_LESS_FH = 0x00000004, /* missing fixed H */
    INCHIDIFF_MORE_H = 0x00000008, /* formulas differ in number of H */
    INCHIDIFF_LESS_H = 0x00000008, /* formulas differ in number of H */
    INCHIDIFF_NO_TAUT = 0x00000010, /* restored structure has no taut groups while the original InChI has some */
    INCHIDIFF_WRONG_TAUT = 0x00000020, /* restored has tautomerism while the original does not have it */
    INCHIDIFF_SINGLE_TG = 0x00000040, /* restored has 1 taut. group while the original InChI has multiple tg */
    INCHIDIFF_MULTIPLE_TG = 0x00000080, /* restored has multiple tg while the original InChI has only one tg */
    INCHIDIFF_EXTRA_TG_ENDP = 0x00000100, /* extra tautomeric endpoint(s) in restored structure */
    INCHIDIFF_MISS_TG_ENDP = 0x00000100, /* one or more tg endpoint is not in the restored structure */
    INCHIDIFF_DIFF_TG_ENDP = 0x00000100, /* lists of tg endpoints are different */
    INCHIDIFF_NUM_TG = 0x00000200, /* different number of tautomeric groups */
    INCHIDIFF_TG = 0x00000200, /* different tautomeric groups */
    INCHIDIFF_NUM_ISO_AT = 0x00000400, /* ?severe: restored struct. has different number of isotopic atoms */
    INCHIDIFF_ISO_AT = 0x00000400, /* ?severe: restored struct. has different locations/isotopes of isotopic atoms */
    INCHIDIFF_REM_ISO_H = 0x00000800, /* isotopic H removed */
    INCHIDIFF_MOB_ISO_H = 0x00001000, /* different number of mobile exchangeable isotopic H */
    INCHIDIFF_CHARGE = 0x00002000, /* restored structure has different charge */
    INCHIDIFF_REM_PROT = 0x00004000, /* proton(s) removed/added from the restored structure */
    INCHIDIFF_MOBH_PROTONS = 0x00008000, /* different proton balance */
    INCHIDIFF_SC_INV = 0x00010000, /* restores structure has different inversion stereocenter mark */
    INCHIDIFF_SC_PARITY = 0x00020000, /* restored structure has stereoatoms or allenes with different parity */
    INCHIDIFF_SC_EXTRA_UNDF = 0x00040000, /* restored structure has extra undefined stereocenter(s) */
    INCHIDIFF_SC_EXTRA = 0x00080000, /* restored structure has extra stereocenter(s) */
    INCHIDIFF_SC_MISS_UNDF = 0x00100000,  /* restored structure has not some undefined stereocenter(s) */
    INCHIDIFF_SC_MISS = 0x00200000, /* restored structure has not some stereocenters that are not undefined */
    INCHIDIFF_SB_PARITY = 0x00400000, /* restored structure has stereobonds or cumulenes with different parity */
    INCHIDIFF_SB_EXTRA_UNDF = 0x00800000, /* restored structure has extra undefined stereobond(s) */
    INCHIDIFF_SB_EXTRA = 0x01000000, /* restored structure has extra stereobond(s) */
    INCHIDIFF_SB_MISS_UNDF = 0x02000000, /* restored structure has not some undefined stereocenters */
    INCHIDIFF_SB_MISS = 0x04000000, /* restored structure has not some stereobonds that are not undefined */
    INCHIDIFF_COMP_HLAYER = 0x08000000, /* Restored component has Mobile-H layer instead of both Mobile-H & Fixed-H or both instead of one */
    INCHIDIFF_COMP_NUMBER = 0x10000000, /* wrong number of components */
    INCHIDIFF_STR2INCHI_ERR = 0x20000000  /* Restored structure to InChI conversion error */
    /* reserved
                              0x40000000
                              0x80000000
    */
} INCHIDIFF;

typedef enum tagtagCompareInchiMsgGroupID {
    IDGRP_ZERO = 0,
    IDGRP_ERR = 1,
    IDGRP_H = 2,
    IDGRP_MOB_GRP = 3,
    IDGRP_ISO_AT = 4,
    IDGRP_CHARGE = 5,
    IDGRP_PROTONS = 6,
    IDGRP_ISO_H = 7,
    IDGRP_SC = 8,
    IDGRP_SB = 9,
    IDGRP_HLAYER = 10,
    IDGRP_COMP = 11,
    IDGRP_CONV_ERR = 12
} CMP_INCHI_MSG_GROUP_ID;


typedef struct tagCompareInchiMsg {
    INCHIDIFF              nBit;
    CMP_INCHI_MSG_GROUP_ID nGroupID;
    const char            *szMsg;
} CMP_INCHI_MSG;

typedef struct tagCompareInchiMsgGroup {
    CMP_INCHI_MSG_GROUP_ID   nGroupID;
    const char              *szGroupName;
} CMP_INCHI_MSG_GROUP;


#define INCHIDIFF_SB (INCHIDIFF_SB_PARITY | INCHIDIFF_SB_EXTRA_UNDF | INCHIDIFF_SB_EXTRA | INCHIDIFF_SB_MISS_UNDF | INCHIDIFF_SB_MISS)
#define INCHIDIFF_SC (INCHIDIFF_SC_PARITY | INCHIDIFF_SC_EXTRA_UNDF | INCHIDIFF_SC_EXTRA | INCHIDIFF_SC_MISS_UNDF | INCHIDIFF_SC_MISS)

#define INCHIDIFF_CONSTIT (INCHIDIFF_POSITION_H | INCHIDIFF_MORE_FH | INCHIDIFF_LESS_FH | INCHIDIFF_MORE_H | INCHIDIFF_LESS_H |\
                       INCHIDIFF_NO_TAUT | INCHIDIFF_WRONG_TAUT | INCHIDIFF_SINGLE_TG | INCHIDIFF_MULTIPLE_TG | \
                       INCHIDIFF_NUM_TG | INCHIDIFF_EXTRA_TG_ENDP | INCHIDIFF_MISS_TG_ENDP | INCHIDIFF_TG | \
                       INCHIDIFF_NUM_ISO_AT | INCHIDIFF_ISO_AT | INCHIDIFF_CHARGE | INCHIDIFF_REM_PROT | INCHIDIFF_REM_ISO_H |\
                       INCHIDIFF_DIFF_TG_ENDP)
#define INCHIDIFF_STEREO  (INCHIDIFF_SC_INV | INCHIDIFF_SC_PARITY | INCHIDIFF_SC_EXTRA_UNDF | INCHIDIFF_SC_EXTRA | \
                       INCHIDIFF_SC_MISS_UNDF | INCHIDIFF_SC_MISS | INCHIDIFF_SB_PARITY | INCHIDIFF_SB_EXTRA_UNDF |\
                       INCHIDIFF_SB_EXTRA | INCHIDIFF_SB_MISS_UNDF | INCHIDIFF_SB_MISS)


#endif    /* _INCHICMP_H_ */

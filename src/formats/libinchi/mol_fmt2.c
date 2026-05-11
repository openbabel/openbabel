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

#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>

#include "mode.h"
#include "mol_fmt.h"

#include "ichierr.h"
#include "util.h"
#include "ichi_io.h"

#include "bcf_s.h"

/*
    MolFile related procedures - 2

*/

/****************************************************************************
 Read n chars and find where they are terminated with space or trailing 0
****************************************************************************/
int MolfileStrnread(char *dest,
                    char *source,
                    int len,
                    char **first_space)
{
    /* required len >= 0; dest must have at least len+1 bytes */

    int i, c;

    if (len > 0)
    {
        strncpy(dest, source, len);
    }
    dest[len] = '\0';

    len = (len > 0) ? (int)strlen(dest) : 0;

    for (i = ( len - 1 ); i >= 0 && 0 != ( c = source[i] ) && isspace( UCINT c ); i--);

    *first_space = dest + ((long long)i + 1); /* first blank or zero terminating byte in dest */ /* djb-rwth: cast operator added */

    return len; /* number of actually processed bytes excluding zero terminator */
}

/****************************************************************************
 * Extract the 'data' in the mol file field at given text position 'line_ptr'
 *
 *
 * 1. 'field_len' for MOL_FMT_STRING_DATA does not include trailing zero,
 *     that is actual length of the string pointed by 'data'
 *     should be at least field_len+1 bytes.
 *     For numerical data 'field_len' is length of input data field
 *     For numerical integral data field_len <= 0 means read up to first
 *     non-numeric character as strtod() does ("free format")
 * 2.  return value: for MOL_FMT_STRING_DATA: number of bytes excluding trailing zero
 *                   for all others:  1=success; 0 = empty; -1= error
 * 3.  on exit *line_ptr points to the next byte after the last entered
 *
 *
 ****************************************************************************/
int MolfileReadField(void *data,
                     int field_len,
                     int data_type,
                     char **line_ptr)
{
    char *p = *line_ptr, *q, *p_end;
    int  i, c, len, ret = 1;
    long ldata;
    double ddata;

    int DEFINITE_LENGTH_FIELD = 0;
    int FIELD_ENDS_AT_FIRST_NON_DIGIT = 0;
    int TOO_LONG_FIELD = 0;

    if (field_len > MOL_FMT_MAX_VALUE_LEN)
    {
        TOO_LONG_FIELD = 1;
    }
    else if (field_len <= 0)
    {
        FIELD_ENDS_AT_FIRST_NON_DIGIT = 1;
    }
    else
    {
        DEFINITE_LENGTH_FIELD = 1;
    }

    switch (data_type)
    {
    case MOL_FMT_STRING_DATA:
        /* pass by all leading spaces */
        for (i = 0;
             i < field_len && 0 != (c = p[i]) && isspace(UCINT c);
             i++)
        {
            ;
        }

        len = MolfileStrnread((char *)data, &p[i], field_len - i, &q);

            ret = ( q - (char*) data );/* actual data length */
            *q = '\0'; /* add zero termination to data if it is not there yet*/
            *line_ptr += ( (long long)len + (long long)i ); /* ptr to the 1st byte of the next input field or to zero termination */ /* djb-rwth: cast operators added */
        break;

    case MOL_FMT_CHAR_INT_DATA:
    case MOL_FMT_SHORT_INT_DATA:
    case MOL_FMT_LONG_INT_DATA:
    {
        char str[MOL_FMT_MAX_VALUE_LEN + 1];
        ldata = 0L;
        if (TOO_LONG_FIELD)
        {
            ret = -1;
        }
        else if (DEFINITE_LENGTH_FIELD)
        {
            /* fixed length */
            *line_ptr += (len = MolfileStrnread(str, p, field_len, &q));

            *q = '\0';
            if (!len || !(q - str))
            {
                ret = 0; /* empty string */
            }
            else
            {
                if ((ldata = strtol(str, &p_end, 10), p_end != q))
                {
                    ret = -1; /* wrong data: incompletely interpreted */
                }
            }
        }
        else if (FIELD_ENDS_AT_FIRST_NON_DIGIT)
        {
            /* free format: field_len <= 0 */
            ldata = strtol(p, &p_end, 10);
            *line_ptr += (len = p_end - p);
            if (len == 0)
            {
                ret = 0;
            }
        }
        else
        {
            /* should not come here */
            ret = -1;
        }

        switch (data_type)
        {
        case MOL_FMT_CHAR_INT_DATA:
            if (SCHAR_MIN <= ldata && ldata <= SCHAR_MAX)
            {
                /* from || to &&: 11-19-96 */
                *(S_CHAR *)data = (S_CHAR)ldata;
            }
            else
            {
                *(S_CHAR *)data = (S_CHAR)0;
                ret = -1;
            }
            break;
        case MOL_FMT_SHORT_INT_DATA:
            if (SHRT_MIN <= ldata && ldata <= SHRT_MAX)
            {
                *(S_SHORT *)data = (S_SHORT)ldata;
            }
            else
            {
                *(S_SHORT *)data = (S_SHORT)0;
                ret = -1;
            }
            break;
        case MOL_FMT_LONG_INT_DATA:
            if (LONG_MIN < ldata && ldata < LONG_MAX)
            {
                *(long *)data = (long)ldata;
            }
            else
            {
                *(long *)data = 0L;
                ret = -1;
            }
            break;
        default:
            ret = -1;
        }
    } /* MOL_FMT_CHAR_INT_DATA... */
    break;

    case MOL_FMT_DOUBLE_DATA:
    case MOL_FMT_FLOAT_DATA:
    {
        char str[MOL_FMT_MAX_VALUE_LEN + 1];

        if (TOO_LONG_FIELD)
        {
            ret = -1;
            ddata = 0.0;
        }
        else if (DEFINITE_LENGTH_FIELD)
        {
            *line_ptr += (len = MolfileStrnread(str, p, field_len, &q));
            *q = '\0';
            if (!len || !(q - str))
            {
                /* empty string */
                ddata = 0.0;
                ret = 0;
            }
            else if ((ddata = strtod(str, &p_end), p_end != q))
            {
                /* wrong data */
                ret = -1;
            }
        }
        else if (FIELD_ENDS_AT_FIRST_NON_DIGIT)
        {
            /* free format */
            ddata = strtod(p, &p_end);
            *line_ptr += (len = p_end - p);
            if (len == 0)
            {
                ret = 0;
            }
        }
        else
        {
            /* should not come here */
            ret = -1; /* djb-rwth: addressing coverity ID #499478 -- see the original comment above */
        }

        switch (data_type)
        {

        case MOL_FMT_DOUBLE_DATA:
            if (ddata != HUGE_VAL && /*ldata*/ ddata != -HUGE_VAL)
            { /* replaced ldata with ddata 6-30-98 DCh */
                *(double *)data = ddata;
            }
            else
            {
                *(double *)data = 0.0;
                ret = -1;
            }
            break;

        case MOL_FMT_FLOAT_DATA:
            if (fabs(ddata) <= (double)FLT_MIN)
            {
                *(float *)data = 0.0;
            }
            else if (fabs(ddata) >= (double)FLT_MAX)
            {
                *(float *)data = 0.0;
                ret = -1;
            }
            else
            {
                *(float *)data = (float)ddata;
            }
            break;
        }
    } /* MOL_FMT_DOUBLE_DATA... */
    break;

    case MOL_FMT_JUMP_TO_RIGHT:
    {

        for (i = 0; i < field_len && p[i]; i++)
            ;

        *line_ptr += i;
        ret = i;
    }
    break;

    default:
        ret = -1;
    }

    return ret;
}

/****************************************************************************
 Read molfile number from the name line like "Structure #22"
****************************************************************************/
long MolfileExtractStrucNum(MOL_FMT_HEADER_BLOCK *pHdr)
{
    static char sStruct[] = "Structure #";
    static char sINCHI[] = INCHI_NAME;
    long   lMolfileNumber = 0;
    char   *p, *q = NULL;

    if (pHdr)
    {
        if (!inchi_memicmp(pHdr->molname, sStruct, sizeof(sStruct) - 1))
        {
            p = pHdr->molname + sizeof(sStruct) - 1;
            lMolfileNumber = strtol(p, &q, 10);
            p = pHdr->line2;
            if (!q || *q ||
                inchi_memicmp(p, sINCHI, sizeof(sINCHI) - 1) ||
                !strstr(p + sizeof(sINCHI) - 1, "SDfile Output"))
            {
                lMolfileNumber = 0;
            }
        }
    }

    return lMolfileNumber;
}

/****************************************************************************
 Check if MOL file contains no structure
****************************************************************************/
int MolfileHasNoChemStruc(MOL_FMT_DATA *mfdata)
{
    if (!mfdata || !mfdata->ctab.atoms)
    {
        return 1;
    }

    if (mfdata->ctab.n_atoms <= 0)
    {
        return 1;
    }

    if (0 < mfdata->ctab.n_bonds && !mfdata->ctab.bonds)
    {
        return 1;
    }

    return 0;
}

/****************************************************************************
 Copy MOL-formatted data of SDF record or Molfile to another file
****************************************************************************/
int MolfileSaveCopy(INCHI_IOSTREAM *inp_file,
                    long fPtrStart,
                    long fPtrEnd,
                    FILE *outfile,
                    long num)
{
    char line[MOL_FMT_INPLINELEN], *p;
    long fPtr;
    int ret = 1;
    char szNumber[32];

    if (inp_file->type == INCHI_IOS_TYPE_FILE)
    {

        FILE *infile = inp_file->f;

        if (!infile)
        {
            return 1;
        }

        if (!outfile)
        {
            return 1;
        }

        if (fPtrStart < 0L && fPtrEnd <= fPtrStart)
        {
            return 1;
        }

        if (0 != fseek(infile, fPtrStart, SEEK_SET))
        {
            return 1;
        }

        while (fPtrEnd > (fPtr = ftell(infile)) && fPtr >= 0L
                && inchi_fgetsLf(line, sizeof(line) - 1, inp_file))
        {

            line[sizeof(line) - 1] = '\0'; /*  unnecessary extra precaution */

            if (fPtr == fPtrStart && num)
            {
                int len;
                lrtrim(line, &len);
                len = sprintf(szNumber, "#%ld%s", num, len ? "/" : "");
                mystrncpy(line + len, line, sizeof(line) - len - 1);
                memcpy(line, szNumber, len);
            }

            if (!strchr(line, '\n'))
            {
                p = line + strlen(line);
                p[0] = '\n';
                p[1] = '\0';
            }

            fputs(line, outfile);
        }

        ret = fseek(infile, fPtrEnd, SEEK_SET);
    }
    else if (inp_file->type == INCHI_IOS_TYPE_STRING)
    {
        ;
    }
    else
    {
        ;
    }

    return ret;
}


#define MIN_STDATA_X_COORD           0.0
#define MAX_STDATA_X_COORD         256.0
#define MIN_STDATA_Y_COORD           0.0
#define MAX_STDATA_Y_COORD         256.0
#define MIN_STDATA_Z_COORD           0.0
#define MAX_STDATA_Z_COORD         256.0
#define MAX_STDATA_AVE_BOND_LENGTH  20.0
#define MIN_STDATA_AVE_BOND_LENGTH  10.0


/****************************************************************************
 Get xyz dimensionality and normalization factors
****************************************************************************/
int MolfileGetXYZDimAndNormFactors(MOL_FMT_DATA *mfdata,
                                   int find_norm_factors,
                                   double *x0,
                                   double *y0,
                                   double *z0,
                                   double *xmin,
                                   double *ymin,
                                   double *zmin,
                                   double *scaler,
                                   int *err,
                                   char *pStrErr)

{
    int i;
    int num_dimensions = 0, num_atoms, num_bonds;
    double max_x = -1.0e32, max_y = -1.0e32, max_z = -1.0e32;
    double min_x = 1.0e32, min_y = 1.0e32, min_z = 1.0e32;
    double macheps = 1.0e-10, small_coeff = 0.00001;
    double x_coeff, y_coeff, z_coeff, coeff = 1.0, average_bond_length;

    *x0 = MIN_STDATA_X_COORD;
    *y0 = MIN_STDATA_Y_COORD;
    *z0 = MIN_STDATA_Z_COORD;
    *xmin = *ymin = *zmin = 0.0;
    *scaler = coeff;

    if (MolfileHasNoChemStruc(mfdata))
    {
        goto exit_function;
    }

    num_atoms = mfdata->ctab.n_atoms;
    for (i = 0; i < num_atoms; i++)
    {
        max_x = inchi_max(mfdata->ctab.atoms[i].fx, max_x);
        min_x = inchi_min(mfdata->ctab.atoms[i].fx, min_x);
        max_y = inchi_max(mfdata->ctab.atoms[i].fy, max_y);
        min_y = inchi_min(mfdata->ctab.atoms[i].fy, min_y);
        max_z = inchi_max(mfdata->ctab.atoms[i].fz, max_z);
        min_z = inchi_min(mfdata->ctab.atoms[i].fz, min_z);
    }

    num_bonds = 0;
    average_bond_length = 0.0;
    for (i = 0; i < mfdata->ctab.n_bonds; i++)
    {
        double dx, dy, dz;
        int a1 = mfdata->ctab.bonds[i].atnum1 - 1;
        int a2 = mfdata->ctab.bonds[i].atnum2 - 1;

        if (a1 < 0 || a1 >= num_atoms ||
            a2 < 0 || a2 >= num_atoms ||
            a1 == a2)
        {
            *err |= 1; /*  bond for invalid atom number(s); ignored */
            TREAT_ERR(*err, 0, "Bond to nonexistent atom");
            continue;
        }

        dx = mfdata->ctab.atoms[a1].fx - mfdata->ctab.atoms[a2].fx;
        dy = mfdata->ctab.atoms[a1].fy - mfdata->ctab.atoms[a2].fy;
        dz = mfdata->ctab.atoms[a1].fz - mfdata->ctab.atoms[a2].fz;

        average_bond_length += sqrt(dx*dx + dy*dy + dz*dz);
        num_bonds++;
    }

    if (max_x - min_x <= small_coeff * (fabs(max_x) + fabs(min_x)))
    {
        x_coeff = 0.0;
    }
    else
    {
        x_coeff = (MAX_STDATA_X_COORD - MIN_STDATA_X_COORD) / (max_x - min_x);
    }

    if (max_y - min_y <= small_coeff * (fabs(max_y) + fabs(min_y)))
    {
        y_coeff = 0.0;
    }
    else
    {
        y_coeff = (MAX_STDATA_Y_COORD - MIN_STDATA_Y_COORD) / (max_y - min_y);
    }

    if (max_z - min_z <= small_coeff * (fabs(max_z) + fabs(min_z)))
    {
        z_coeff = 0.0;
    }
    else
    {
        z_coeff = (MAX_STDATA_Z_COORD - MIN_STDATA_Z_COORD) / (max_z - min_z);
    }

    num_dimensions = ((x_coeff > macheps || y_coeff > macheps) && fabs(z_coeff) < macheps)
                        ? 2
                        : ( fabs( z_coeff ) > macheps ) ? 3 : 0;


    if (!find_norm_factors)
    {
        goto exit_function;
    }

    /* Find normalization parameters */
    switch (num_dimensions)
    {
    case 0:
        coeff = 0.0;
        break;

    case 2:
        /* choose the smallest stretching coefficient */
        if (x_coeff > macheps && y_coeff > macheps)
        {
            coeff = inchi_min(x_coeff, y_coeff);
        }
        else if (x_coeff > macheps)
        {
            coeff = x_coeff;
        }
        else if (y_coeff > macheps)
        {
            coeff = y_coeff;
        }
        else
        {
            coeff = 1.0;
        }
        break;

    case 3:
        /* choose the smallest stretching coefficient */
        if (x_coeff > macheps && y_coeff > macheps)
        {
            coeff = inchi_min(x_coeff, y_coeff);
            coeff = inchi_min(coeff, z_coeff);
        }
        else if (x_coeff > macheps)
        {
            coeff = inchi_min(x_coeff, z_coeff);
        }
        else if (y_coeff > macheps)
        {
            coeff = inchi_min(y_coeff, z_coeff);
        }
        else
        {
            coeff = z_coeff;
        }
        break;

    default:
        coeff = 0.0;
    }

    if (num_bonds > 0)
    {

        average_bond_length /= (double)num_bonds;
        if (average_bond_length * coeff > MAX_STDATA_AVE_BOND_LENGTH)
        {
            coeff = MAX_STDATA_AVE_BOND_LENGTH / average_bond_length; /* avoid too long bonds */
        }
        else if (average_bond_length * coeff < macheps)
        {
            coeff = 1.0; /* all lengths are of zero length */
        }
        else if (average_bond_length * coeff < MIN_STDATA_AVE_BOND_LENGTH)
        {
            coeff = MIN_STDATA_AVE_BOND_LENGTH / average_bond_length; /* avoid too short bonds */
        }
    }

exit_function:;

    *x0 = min_x;
    *y0 = min_y;
    *z0 = min_z;
    *xmin = MIN_STDATA_X_COORD;
    *ymin = MIN_STDATA_Y_COORD;
    *zmin = MIN_STDATA_Z_COORD;
    *scaler = coeff;

    return num_dimensions;
}

/****************************************************************************
 Clean up MOL-format parser data
****************************************************************************/
MOL_FMT_DATA *FreeMolfileData(MOL_FMT_DATA *mfdata)
{
    if (mfdata)
    {

        if (mfdata->ctab.atoms)
        {
            inchi_free(mfdata->ctab.atoms);
        }

        if (mfdata->ctab.bonds)
        {
            inchi_free(mfdata->ctab.bonds);
        }

        if (mfdata->ctab.coords)
        {
            inchi_free(mfdata->ctab.coords);
        }

        /*if ( 0!=mfdata->ctab.sgroups.used )*/
        MolFmtSgroups_Free(&(mfdata->ctab.sgroups));

        if (mfdata->ctab.v3000)
        {
            DeleteMolfileV3000Info(mfdata->ctab.v3000);
        }

        inchi_free(mfdata);
        mfdata = NULL;
    }

    return mfdata;
}

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
#include "stb_sprintf.h"

/*
    Molfile V3000 related procedures

*/

static int get_actual_atom_number(int index, int n, int *orig, int *fin);

/****************************************************************************
 Init V3000 reader
****************************************************************************/
int MolfileV3000Init(MOL_FMT_CTAB *ctab,
                     char *pStrErr)
{
    int ret = 0;
    int i;

    /* STAR ATOMS */
    ctab->v3000->n_star_atoms = 0;
    ctab->v3000->n_non_star_atoms = 0;

    if (ctab->n_atoms)
    {
        ctab->v3000->atom_index_orig = (int *)inchi_calloc(ctab->n_atoms, sizeof(int));
        ctab->v3000->atom_index_fin = (int *)inchi_calloc(ctab->n_atoms, sizeof(int));
        if (ctab->v3000->atom_index_orig && ctab->v3000->atom_index_fin) /* djb-rwth: fixing a NULL pointer dereference */
        {
            for (i = 0; i < ctab->n_atoms; i++) /* protective */
            {
                ctab->v3000->atom_index_orig[i] = -1;
                ctab->v3000->atom_index_fin[i] = -1;
            }
        }
    }
    else
    {
        ctab->v3000->atom_index_orig = NULL;
        ctab->v3000->atom_index_fin = NULL;
    }

    /* HAPTIC BONDS */
    ctab->v3000->n_haptic_bonds = 0;
    ctab->v3000->haptic_bonds = (NUM_LISTS *)inchi_calloc(1, sizeof(NUM_LISTS));
    if (!ctab->v3000->haptic_bonds)
    {
        AddErrorMessage(pStrErr, "Out of RAM");
        return -1;
    }
    ret = NumLists_Alloc(ctab->v3000->haptic_bonds, 8);
    if (ret < 0)
    {
        ctab->v3000->haptic_bonds = NULL;
        AddErrorMessage(pStrErr, "Out of RAM");
        return -1;
    }

    /* STEABS */
    ctab->v3000->n_steabs = 0;
    ctab->v3000->steabs = (NUM_LISTS *)inchi_calloc(1, sizeof(NUM_LISTS));
    if (!ctab->v3000->steabs)
    {
        AddErrorMessage(pStrErr, "Out of RAM");
        return -1;
    }
    ret = NumLists_Alloc(ctab->v3000->steabs, 1);
    if (ret < 0)
    {
        ctab->v3000->steabs = NULL;
        AddErrorMessage(pStrErr, "Out of RAM");
        return -1;
    }
    /* STEREL */
    ctab->v3000->n_sterel = 0;
    ctab->v3000->sterel = (NUM_LISTS *)inchi_calloc(1, sizeof(NUM_LISTS));
    if (!ctab->v3000->sterel)
    {
        AddErrorMessage(pStrErr, "Out of RAM");
        return -1;
    }
    ret = NumLists_Alloc(ctab->v3000->sterel, 4);
    if (ret < 0)
    {
        ctab->v3000->sterel = NULL;
        AddErrorMessage(pStrErr, "Out of RAM");
        return -1;
    }
    /* STERAC */
    ctab->v3000->n_sterac = 0;
    ctab->v3000->sterac = (NUM_LISTS *)inchi_calloc(1, sizeof(NUM_LISTS));
    if (!ctab->v3000->sterac)
    {
        AddErrorMessage(pStrErr, "Out of RAM");
        return -1;
    }
    ret = NumLists_Alloc(ctab->v3000->sterac, 4);
    if (ret < 0)
    {
        ctab->v3000->sterac = NULL;
        AddErrorMessage(pStrErr, "Out of RAM");
        return -1;
    }

    return ret;
}

/****************************************************************************
 Delete V3000 reader data
****************************************************************************/
int DeleteMolfileV3000Info(MOL_FMT_v3000 *v3000)
{
    if (v3000)
    {

        if (v3000->atom_index_orig)
        {
            inchi_free(v3000->atom_index_orig);
        }

        if (v3000->atom_index_fin)
        {
            inchi_free(v3000->atom_index_fin);
        }

        if (v3000->haptic_bonds)
        {
            NumLists_Free(v3000->haptic_bonds);
            free(v3000->haptic_bonds);
        }

        if (v3000->steabs)
        {
            NumLists_Free(v3000->steabs);
            free(v3000->steabs);
        }

        if (v3000->sterel)
        {
            NumLists_Free(v3000->sterel);
            free(v3000->sterel);
        }

        if (v3000->sterac)
        {
            NumLists_Free(v3000->sterac);
            free(v3000->sterac);
        }

        inchi_free(v3000);
        v3000 = NULL;
    }

    return 0;
}

/****************************************************************************
    Extended version of inchi_fgetsLf which is able of reading
    concatenated lines (ending with '-') of V3000 Molfile.
    Also removes "M  V30 " prefix" and normalizes the rest of string
****************************************************************************/
char *inchi_fgetsLf_V3000(char *line, INCHI_IOSTREAM *inp_stream)
{
    char *p = NULL;
    int len = 0;

    p = inchi_fgetsLf(line, MOL_FMT_V3000_INPLINELEN, inp_stream);
    if (!p)
    {
        return NULL;
    }

    len = (int)strlen(p);
    if (len < 7)
    {
        return NULL;
    }

    if (strncmp(p, "M  V30 ", 7))
    {
        return NULL;
    }

    p += 7;
    len = normalize_string(p); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */

    return p;
}

/****************************************************************************
    Read V3000 field.

    It is MolfileReadField updated for V3000.
    It considers right whitespace as stop sign,
    no predefined len is used.

    Returns -1 on error otherwise number of bytes read.

    NB: ASSUMES THAT STRING HAS BEEN NORMALIZED with normalize_string()

    TODO: treat strings with spaces in double quotes
****************************************************************************/
int MolfileV3000ReadField(void *data,
                          int data_type,
                          char **line_ptr)
{
    int nread = 0;
    char field[MOL_FMT_V3000_MAXFIELDLEN];
    const int max_field_len = sizeof(field);
    long ldata = 0L;
    double ddata = 0.0;
    char *p_end;

    memset(field, 0, max_field_len); /* djb-rwth: memset_s C11/Annex K variant? */

    nread = read_upto_delim(line_ptr, field, max_field_len, " \t\n\v\f\r");

    switch (data_type)
    {
    case MOL_FMT_STRING_DATA:
    {
        if (nread && (nread <= ATOM_EL_LEN)) /* djb-rwth: fixing GHI #133 -- updated 28/09/2025 */
        {
            mystrncpy((char *)data, field, nread + 1);
        }
        else
        {
            ((char *)data)[0] = '\0';
        }
    }
    break;

    case MOL_FMT_CHAR_INT_DATA:
    case MOL_FMT_SHORT_INT_DATA:
    case MOL_FMT_LONG_INT_DATA:
    case MOL_FMT_INT_DATA:
    {
        /* assume that field ends at first non-digit */
        ldata = strtol(field, &p_end, 10);

        if (p_end == field)
        {
            nread = 0;
        }

        if (data_type == MOL_FMT_LONG_INT_DATA)
        {
            if (LONG_MIN < ldata && ldata < LONG_MAX)
            {
                *(long *)data = (long)ldata;
            }
            else
            {
                *(long *)data = 0L;
                nread = -1;
            }
        }
        else if (data_type == MOL_FMT_INT_DATA)
        {
            if (INT_MIN <= ldata && ldata <= INT_MAX) /* djb-rwth: addressing coverity ID #499496/499553 -- ldata check seems to be necessary */
            {
                *(int *)data = (int)ldata;
            }
            else
            {
                *(int *)data = (int)0;
                nread = -1;
            }
        }

        else if (data_type == MOL_FMT_CHAR_INT_DATA)
        {
            if (SCHAR_MIN <= ldata && ldata <= SCHAR_MAX)
            {
                *(S_CHAR *)data = (S_CHAR)ldata;
            }
            else
            {
                *(S_CHAR *)data = (S_CHAR)0;
                nread = -1;
            }
        }

        else if (data_type == MOL_FMT_SHORT_INT_DATA)
        {
            if (SHRT_MIN <= ldata && ldata <= SHRT_MAX)
            {
                *(S_SHORT *)data = (S_SHORT)ldata;
            }
            else
            {
                *(S_SHORT *)data = (S_SHORT)0;
                nread = -1;
            }
        }
        else
        {
            nread = -1;
        }
    }
    break; /* INT's */

    case MOL_FMT_DOUBLE_DATA:
    case MOL_FMT_FLOAT_DATA:
    {
        /* assume that field ends at first non-digit */
        ddata = strtod(field, &p_end);

        if (p_end == field)
        {
            nread = 0;
        }

        if (data_type == MOL_FMT_DOUBLE_DATA)
        {
            if (ddata != HUGE_VAL && /*ldata*/ ddata != -HUGE_VAL)
            {
                *(double *)data = ddata;
            }
            else
            {
                *(double *)data = 0.0;
                nread = -1;
            }
        }
        else if (data_type == MOL_FMT_FLOAT_DATA)
        {
            if (fabs(ddata) <= (double)FLT_MIN)
            {
                *(float *)data = 0.0;
            }
            else if (fabs(ddata) >= (double)FLT_MAX)
            {
                *(float *)data = 0.0;
                nread = -1;
            }
        }
        else
        {
            *(float *)data = (float)ddata; /* djb-rwth: addressing coverity ID #499519 -- probably never reached */
        }
    }
    break; /* REAL's */

    default:
        nread = -1;
    }

    return nread;
}

/****************************************************************************
    Read V3000 keyword.
    TODO: treat strings with spaces in double quotes
****************************************************************************/
int MolfileV3000ReadKeyword(char *key,
                            char **line_ptr)
{
    int nread = 0;
    char field[MOL_FMT_V3000_MAXFIELDLEN];
    const int max_field_len = sizeof(field);

    memset(field, 0, max_field_len); /* djb-rwth: memset_s C11/Annex K variant? */

    nread = read_upto_delim(line_ptr, field, max_field_len, "= \t\n\v\f\r");

    if (nread)
    {
        mystrncpy(key, field, nread + 1);
        /* consume '=' sign if present */
        if (*line_ptr)
        {
            if (*line_ptr[0] == '=')
            {
                *line_ptr = *line_ptr + 1;
            }
        }
    }
    else
    {
        key[0] = '\0';
    }

    return nread;
}

/****************************************************************************
 Read V3000 head of CTab
****************************************************************************/
int MolfileV3000ReadCTABBeginAndCountsLine(MOL_FMT_CTAB *ctab,
                                           INCHI_IOSTREAM *inp_file,
                                           char *pStrErr)
{
    char field[MOL_FMT_V3000_MAXFIELDLEN];
    int err = 0, len; /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    int failed = 0;

    int nc;
    char *p = NULL, *line = NULL;
    INCHI_IOSTREAM tmpin;
    INCHI_IOS_STRING *pin = &tmpin.s;
    inchi_ios_init(&tmpin, INCHI_IOS_TYPE_STRING, NULL);

    field[0] = '\0'; /* djb-rwth: adding zero termination */

    /* Check for proper start */

    /*p = inchi_fgetsLf_V3000( line, inp_file );*/
    inchi_strbuf_reset(pin);
    nc = get_V3000_input_line_to_strbuf(pin, inp_file);
    if (nc < 1)
    {
        p = NULL;
    }
    else
    {
        p = line = pin->pStr;
    }
    if (!p || strcmp(p, "BEGIN CTAB"))
    {
        TREAT_ERR_AND_FIN(err, 1, err_fin, "Error: No V3000 CTab start marker");
    }
    remove_one_lf(line);

    /* Reset all previosly read data from quasi-counts line            */
    /* (which contains only single meaningful value, 'V3000' marker    */
    ctab->n_atoms = -1;
    ctab->n_bonds = -1;
    ctab->chiral_flag = -1;
    ctab->n_stext_entries = -1;
    /* Relax stricthness of V3000 conformance: */
    /* Do not check if '999' supplied, just use this. */
    ctab->n_property_lines = 999;

    /* Read counts line */
    /*p = inchi_fgetsLf_V3000( line, inp_file );*/
    inchi_strbuf_reset(pin);
    nc = get_V3000_input_line_to_strbuf(pin, inp_file);
    if (nc < 1)
    {
        p = NULL;
    }
    else
    {
        p = line = pin->pStr;
    }
    if (!p)
    {
        TREAT_ERR_AND_FIN(err, 1, err_fin, "Cannot read V3000 counts line");
    }
    remove_one_lf(line);

    /* Parse counts line */
    len = MolfileV3000ReadField(field, MOL_FMT_STRING_DATA, &p); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    if (strcmp(field, "COUNTS"))
    {
        TREAT_ERR_AND_FIN(err, 1, err_fin, "Cannot read V3000 counts line");
    }
    failed = 0;
    if (0 > MolfileV3000ReadField(&ctab->n_atoms, MOL_FMT_INT_DATA, &p))
    {
        failed = 2;
    }
    else if (0 > MolfileV3000ReadField(&ctab->n_bonds, MOL_FMT_INT_DATA, &p))
    {
        failed = 1;
    }
    else if (0 > MolfileV3000ReadField(&ctab->v3000->n_sgroups, MOL_FMT_INT_DATA, &p))
    {
        failed = 1;
    }
    else if (0 > MolfileV3000ReadField(&ctab->v3000->n_3d_constraints, MOL_FMT_INT_DATA, &p))
    {
        failed = 1;
    }
    else if (0 > MolfileV3000ReadField(&ctab->chiral_flag, MOL_FMT_CHAR_INT_DATA, &p))
    {
        failed = 1;
    }

    if (failed)
    {
        err = 3;
        if (failed == 2)
        {
            TREAT_ERR(err, 3, "Number of atoms too large. V3000 counts line:");
        }
        else
        {
            /* too long input file line or other value min-max range mismatch */
            TREAT_ERR(err, 3, "Cannot interpret V3000 counts line:");
        }
        dotify_non_printable_chars(line);
        AddErrorMessage(pStrErr, line);
        goto err_fin;
    }

err_fin:
    inchi_strbuf_close(pin);

    return err;
}

/****************************************************************************
 Read V3000 SGroup
****************************************************************************/
int MolfileV3000ReadSGroup(MOL_FMT_CTAB *ctab,
                           INCHI_IOSTREAM *inp_file,
                           int err,
                           char *pStrErr)
{
    int nc;
    char *p = NULL, *line = NULL;
    INCHI_IOSTREAM tmpin;
    INCHI_IOS_STRING *pin = &tmpin.s;

    inchi_ios_init(&tmpin, INCHI_IOS_TYPE_STRING, NULL);
    /*p = inchi_fgetsLf_V3000( line, inp_file );*/

    while (1)
    {
        nc = get_V3000_input_line_to_strbuf(pin, inp_file);

        if (nc < 1)
        {
            p = NULL;
        }
        else
        {
            p = line = pin->pStr;
            remove_one_lf(line);
        }
        if (p && !strcmp(p, "END SGROUP"))
        {
            inchi_ios_close(&tmpin); /* ricrogz: fixing memory leak */
            return 0;
        }
    }

    /* if (  !p  || strcmp(p, "END SGROUP") ) */
    {
        TREAT_ERR_AND_FIN(err, 1, err_fin, "Error: No V3000 SGroup end marker");
    }

err_fin:

    inchi_ios_close(&tmpin); /* ricrogz: fixing memory leak */
    return err;
}

/****************************************************************************
 Read V3000 3DBlock
****************************************************************************/
int MolfileV3000Read3DBlock(MOL_FMT_CTAB *ctab,
                            INCHI_IOSTREAM *inp_file,
                            int err,
                            char *pStrErr)
{
    int nc;
    char *p = NULL, *line = NULL;
    INCHI_IOSTREAM tmpin;
    INCHI_IOS_STRING *pin = &tmpin.s;
    inchi_ios_init(&tmpin, INCHI_IOS_TYPE_STRING, NULL);
    /*p = inchi_fgetsLf_V3000( line, inp_file );*/

    nc = get_V3000_input_line_to_strbuf(pin, inp_file);

    if (nc < 1)
    {
        p = NULL;
    }
    else
    {
        p = line = pin->pStr;
    }
    remove_one_lf(line);

    if (!p || strcmp(p, "END OBJ3D"))
    {
        TREAT_ERR_AND_FIN(err, 1, err_fin, "Error: No V3000 3DBlock end marker");
    }
    goto err_fin;

err_fin:

    inchi_ios_close(&tmpin); /* ricrogz: fixing memory leak */
    return err;
}

/****************************************************************************
 Read V3000 collections
****************************************************************************/
int MolfileV3000ReadCollections(MOL_FMT_CTAB *ctab,
                                INCHI_IOSTREAM *inp_file,
                                int err,
                                char *pStrErr)
{
    char field[MOL_FMT_V3000_MAXFIELDLEN];
    const int max_field_len = sizeof(field);
    int nread, len, n_coll = 0;
    int failed = 0;
    int nc;
    char *p = NULL, *line = NULL;
    INCHI_IOSTREAM tmpin;
    INCHI_IOS_STRING *pin = &tmpin.s;

    inchi_ios_init(&tmpin, INCHI_IOS_TYPE_STRING, NULL);
    /*p = inchi_fgetsLf_V3000( line, inp_file );*/

    nc = get_V3000_input_line_to_strbuf(pin, inp_file);

    if (nc < 1)
    {
        p = NULL;
    }
    else
    {
        p = line = pin->pStr;
    }
    remove_one_lf(line);

    while (p && strcmp(p, "END COLLECTION"))
    {
        int stereo_kind = MOL_FMT_V3000_STENON;
        /* stereo collection of interest */
        NUM_LISTS *ste_coll = NULL;

        nread = read_upto_delim(&p, field, max_field_len, "/");
        if (nread < 6)
        {
            failed = 1;
            break;
        }
        if (strcmp(field, "MDLV30"))
        {
            failed = 1;
            break;
        }

        nread = read_upto_delim(&p, field, max_field_len, "1234567890 \t\n\v\f\r"); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
        if (!strcmp(field, "/STEABS"))
        {
            n_coll = 1;
            stereo_kind = MOL_FMT_V3000_STEABS;
            ste_coll = ctab->v3000->steabs;
        }
        else if (!strcmp(field, "/STEREL"))
        {
            /* get number of collection */
            if (0 > MolfileV3000ReadField(&n_coll, MOL_FMT_CHAR_INT_DATA, &p))
            {
                failed = 1;
                break;
            }
            stereo_kind = MOL_FMT_V3000_STEREL;
            ste_coll = ctab->v3000->sterel;
        }
        else if (!strcmp(field, "/STERAC"))
        {
            /* get number of collection */
            if (0 > MolfileV3000ReadField(&n_coll, MOL_FMT_CHAR_INT_DATA, &p))
            {
                failed = 1;
                break;
            }
            stereo_kind = MOL_FMT_V3000_STERAC;
            ste_coll = ctab->v3000->sterac;
        }
        else
        {
            ;
        }

        if (stereo_kind != MOL_FMT_V3000_STENON)
        /* currently skip non-stereo collections */
        {
            /* consume atoms= */
            if ((len = MolfileV3000ReadKeyword(field, &p) > 0)) /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
            {
                if (!strcmp(field, "ATOMS"))
                {
                    int res, *num_list = NULL;

                    if (0 > MolfileV3000ReadStereoCollection(ctab, &p, &num_list, pStrErr))
                    {
                        failed = 1;
                    }
                    else if (!num_list)
                    {
                        failed = 1;
                    }
                    else
                    {
                        int k, nnum;
                        num_list[0] = n_coll;
                        nnum = num_list[1];
                        for (k = 2; k < nnum; k++)
                        {
                            num_list[k] =
                                get_actual_atom_number(num_list[k],
                                                       ctab->v3000->n_non_star_atoms + ctab->v3000->n_star_atoms,
                                                       ctab->v3000->atom_index_orig,
                                                       ctab->v3000->atom_index_fin);
                        }
                        res = NumLists_Append(ste_coll, num_list);
                        if (res < 0)
                        {
                            failed = 1;
                        }
                        else
                        {
                            if (stereo_kind == MOL_FMT_V3000_STEABS)
                            {
                                ctab->v3000->n_steabs++;
                                ctab->v3000->n_collections++;
                            }
                            else if (stereo_kind == MOL_FMT_V3000_STEREL)
                            {
                                ctab->v3000->n_sterel++;
                                ctab->v3000->n_collections++;
                            }
                            else if (stereo_kind == MOL_FMT_V3000_STERAC)
                            {
                                ctab->v3000->n_sterac++;
                                ctab->v3000->n_collections++;
                            }
                        }
                    }
                }
            }
            else
            {
                failed = 1;
                break;
            }
        }

        /*next_line:*/
        /*p = inchi_fgetsLf_V3000( line, inp_file );*/
        inchi_strbuf_reset(pin);

        nc = get_V3000_input_line_to_strbuf(pin, inp_file);

        if (nc < 1)
        {
            p = NULL;
        }
        else
        {
            p = line = pin->pStr;
            remove_one_lf(line);
        }
    }

    if (failed)
    {
        /*p = inchi_fgetsLf_V3000( line, inp_file );*/
        inchi_strbuf_reset(pin);
        line = NULL; /* Reset line pointer since buffer was freed */

        nc = get_V3000_input_line_to_strbuf(pin, inp_file);

        if (nc < 1)
        {
            p = NULL;
        }
        else
        {
            p = line = pin->pStr;
            remove_one_lf(line);
        }
    }

    if (!p)
    {
        failed = 1;
    }

    if (failed)
    {
        err = 7;
        TREAT_ERR(err, 7, "Cannot interpret V3000 collection line(s)"); /* djb-rwth: addressing coverity ID #499531 -- TREAT_ERR properly used */
        if (line)
        {
            dotify_non_printable_chars(line);
            AddErrorMessage(pStrErr, line);
        }
        goto err_fin;
    }

    /* Error: No V3000 Collection end marker */
    if (ctab->v3000->n_steabs ||
        ctab->v3000->n_sterel ||
        ctab->v3000->n_sterac)
    {
        AddErrorMessage(pStrErr, "V3000 enhanced stereo read/stored but ignored");
    }

err_fin:

    inchi_ios_close(&tmpin); /* ricrogz: fixing memory leak */
    return err;
}

/****************************************************************************
 Read V3000 atoms
****************************************************************************/
int MolfileV3000ReadAtomsBlock(MOL_FMT_CTAB *ctab,
                               INCHI_IOSTREAM *inp_file,
                               int err,
                               char *pStrErr)
{
    int i;
    /* djb-rwth: removing redundant variables */
    char field[MOL_FMT_V3000_MAXFIELDLEN];
    int nc, failed = 0;
    char *p = NULL, *line = NULL;
    INCHI_IOSTREAM tmpin;
    INCHI_IOS_STRING *pin = &tmpin.s;
    inchi_ios_init(&tmpin, INCHI_IOS_TYPE_STRING, NULL);

    /* Check for proper start */
    /*p = inchi_fgetsLf_V3000( line, inp_file );*/

    nc = get_V3000_input_line_to_strbuf(pin, inp_file);

    if (nc < 1)
    {
        p = NULL;
    }
    else
    {
        p = line = pin->pStr;
    }
    if (!p || strcmp(p, "BEGIN ATOM"))
    {
        TREAT_ERR_AND_FIN(err, 1, err_fin, "Error: No V3000 Atom block start marker");
    }
    remove_one_lf(line);

    ctab->v3000->n_non_star_atoms = 0;
    for (i = 0; i < ctab->n_atoms; i++)
    {
        int ii = -1;

        /*p = inchi_fgetsLf_V3000( line, inp_file );*/
        inchi_strbuf_reset(pin);

        nc = get_V3000_input_line_to_strbuf(pin, inp_file);

        if (nc < 1)
        {
            p = NULL;
        }
        else
        {
            p = line = pin->pStr;
        }
        if (!p)
        {
            if (!err)
            {
                TREAT_ERR(err, 2, "Cannot read V3000 atom block line");
            }
            break;
        }
        remove_one_lf(line);

        if (err)
        {
            if (!strcmp(line, SD_FMT_END_OF_DATA))
            {
                err = -abs(err);
                break;
            }
            continue; /* bypass the rest of the Atom block */
        }

        if (ctab->atoms)
        {
            int index, aamap; /* not used actually, just read them */
            int len;
            char symbol[6]; /* TODO: treat possibly long V3000 atom names */
            double fx = 0.0, fy = 0.0, fz = 0.0;
#ifdef GHI100_FIX
#if (SPRINTF_FLAG == 2)
            char *fxs, *fys, *fzs;
            int rfxs, rfys, rfzs;

            fxs = (char *)inchi_malloc((10 + 3) * sizeof(double));
            fys = (char *)inchi_malloc((10 + 3) * sizeof(double));
            fzs = (char *)inchi_malloc((10 + 3) * sizeof(double));

            if (fxs || fys || fzs)
            {
                failed = 1;
            }
#endif
#endif
            symbol[0] = '\0'; /* djb-rwth: adding zero termination */

            /* Read positional parameters */
            failed = 0;

            if (0 > MolfileV3000ReadField(&index, MOL_FMT_INT_DATA, &p))
            {
                failed = 1;
            }
            else if (0 > MolfileV3000ReadField(&symbol, MOL_FMT_STRING_DATA, &p))
            {
                failed = 1;
            }
            else if (0 > MolfileV3000ReadField(&fx, MOL_FMT_DOUBLE_DATA, &p))
            {
                failed = 1;
            }
            else if (0 > MolfileV3000ReadField(&fy, MOL_FMT_DOUBLE_DATA, &p))
            {
                failed = 1;
            }
            else if (0 > MolfileV3000ReadField(&fz, MOL_FMT_DOUBLE_DATA, &p))
            {
                failed = 1;
            }
            else if (0 > MolfileV3000ReadField(&aamap, MOL_FMT_INT_DATA, &p))
            {
                failed = 1;
            }

            if (failed)
            {

                err = 4;
                TREAT_ERR(err, 4, "Cannot interpret V3000 atom block line:"); /* djb-rwth: addressing coverity ID #499547 -- TREAT_ERR properly used */
                dotify_non_printable_chars(line);
                AddErrorMessage(pStrErr, line);

                if (!strcmp(line, SD_FMT_END_OF_DATA))
                {
                    err = -abs(err);
                    break;
                }
                continue; /* can't interpret a first half of atom block line */
            }

            /* emulate V2000 coordinates substring */
            if (ctab->coords)
            {
                char szcoords[40];
#ifdef GHI100_FIX
#if (SPRINTF_FLAG == 2)
                if (fxs)
                {
                    rfxs = dbl2int(fxs, 10, -1, 'g', fx);
                }
                if (fys)
                {
                    rfys = dbl2int(fys, 10, -1, 'g', fy);
                }
                if (fzs)
                {
                    rfzs = dbl2int(fzs, 10, -1, 'g', fz);
                }
                if ((rfxs >= 0) && (rfys >= 0) && (rfzs >= 0))
                {
                    sprintf(szcoords, "%s%s%s", fxs, fys, fzs);
                }
                else
                {
                    failed = 1;
                }
                inchi_free(fxs);
                inchi_free(fys);
                inchi_free(fzs);
#elif (SPRINTF_FLAG == 2)
                stbsp_sprintf(szcoords, "%10g%10g%10g", fx, fy, fz);
#else
                sprintf(szcoords, "%10g%10g%10g", fx, fy, fz);
#endif
#endif
                sprintf(szcoords, "%10g%10g%10g", fx, fy, fz);
                strcpy(ctab->coords[i], szcoords);
            }

            if (!strcmp(symbol, "*"))
            {
                /* ignore star atoms but save index info */
                ctab->v3000->atom_index_orig[i] = index;
                ctab->v3000->atom_index_fin[i] = -1;
                ctab->v3000->n_star_atoms++;
                continue;
            }

            ctab->v3000->n_non_star_atoms++;
            ctab->v3000->atom_index_orig[i] = index;
            ctab->v3000->atom_index_fin[i] = ctab->v3000->n_non_star_atoms;
            ii = ctab->v3000->n_non_star_atoms - 1;

            mystrncpy(ctab->atoms[ii].symbol, symbol, sizeof(ctab->atoms[ii].symbol));
            if (2 == strlen(ctab->atoms[ii].symbol) && isupper(UCINT ctab->atoms[ii].symbol[1]))
            {
                ctab->atoms[ii].symbol[1] = (char)tolower(UCINT ctab->atoms[ii].symbol[1]); /* 5-4-99 DCh*/
            }
            ctab->atoms[ii].fx = fx;
            ctab->atoms[ii].fy = fy;
            ctab->atoms[ii].fz = fz;

            /* Read key-val pairs if any */
            while (p && (len = MolfileV3000ReadKeyword(field, &p)) > 0) /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
            {

                int itmp;
                char ctmp;
                char stmp[MOL_FMT_V3000_MAXFIELDLEN];

                failed = 0;
                if (!strcmp(field, "CHG"))
                {
                    if (0 > MolfileV3000ReadField(&ctab->atoms[ii].charge, MOL_FMT_CHAR_INT_DATA, &p))
                    {
                        failed = 1;
                    }
                }
                else if (!strcmp(field, "RAD"))
                {
                    if (0 > MolfileV3000ReadField(&ctab->atoms[ii].radical, MOL_FMT_CHAR_INT_DATA, &p))
                    {
                        failed = 1;
                    }
                }
                else if (!strcmp(field, "CFG"))
                {
                    if (0 > MolfileV3000ReadField(&ctab->atoms[ii].stereo_parity, MOL_FMT_CHAR_INT_DATA, &p))
                    {
                        failed = 1;
                    }
                }

                else if (!strcmp(field, "MASS"))
                {
                    /*
                        Default = natural abundance
                        A specified value indicates the absolute
                        atomic weight of the designated atom.
                    */
                    S_SHORT iso_mass;
                    if (0 > MolfileV3000ReadField(&iso_mass, MOL_FMT_SHORT_INT_DATA, &p))
                    {
                        failed = 1;
                        TREAT_ERR(err, 0, "Isotopic data not recognized:");
                        AddErrorMessage(pStrErr, line);
                        /* ignore isotopic error for now */
                    }
                    else
                    {
                        /*  What we read is an absolute isotopic mass, by V3000 spec.
                            Adjust this to old convention for further processing:
                            set 'ctab->atoms[ii].mass_difference' to 127
                            if isotopic mass is the same as element mass
                            in Periodic Table (rounded avg by all isotopes), 'atw'
                            delta otherwise, the value of difference 'delta' = ( isotopic mass - 'atw')
                        */
                        int atw, delta;
                        atw = get_atomic_mass(ctab->atoms[ii].symbol);
                        delta = (int)iso_mass - atw;
                        ctab->atoms[ii].mass_difference = (char)(delta ? delta : ZERO_ATW_DIFF);
                    }
                }

                else if (!strcmp(field, "VAL"))
                {
                    if (0 > MolfileV3000ReadField(&itmp, MOL_FMT_INT_DATA, &p))
                    {
                        failed = 1;
                    }
                    else
                    {
                        /* adjust to old convention: was 15 for zero, now -1 for zero */
                        if (itmp == -1)
                        {
                            ctmp = 15;
                        }
                        else
                        {
                            ctmp = (char)itmp;
                        }
                        ctab->atoms[ii].valence = ctmp;
                    }
                }
                else if (!strcmp(field, "HCOUNT"))
                {
                    if (0 > MolfileV3000ReadField(&itmp, MOL_FMT_INT_DATA, &p))
                    {
                        ; /* skip query-related stuff */
                    }
                }
                else if (!strcmp(field, "STBOX"))
                {
                    if (0 > MolfileV3000ReadField(&itmp, MOL_FMT_INT_DATA, &p))
                    {
                        ; /* skip for now */
                    }
                }
                else if (!strcmp(field, "INVRET") || !strcmp(field, "EXACHG"))
                {
                    if (0 > MolfileV3000ReadField(&itmp, MOL_FMT_INT_DATA, &p))
                    {
                        ; /* skip reaction-related stuff */
                    }
                }
                else if (!strcmp(field, "SUBST") || !strcmp(field, "UNSAT") || !strcmp(field, "RBCNT"))
                {
                    if (0 > MolfileV3000ReadField(&itmp, MOL_FMT_INT_DATA, &p))
                    {
                        ; /* skip query-related stuff */
                    }
                }
                else if (!strcmp(field, "ATTCHPT"))
                {
                    if (0 > MolfileV3000ReadField(&itmp, MOL_FMT_INT_DATA, &p))
                    {
                        ;
                    }
                }
                else if (!strcmp(field, "RGROUPS"))
                {
                    if (0 > MolfileV3000ReadField(&stmp, MOL_FMT_STRING_DATA, &p))
                    {
                        ;
                    }
                }
                else if (!strcmp(field, "ATTCHORD"))
                {
                    if (0 > MolfileV3000ReadField(&stmp, MOL_FMT_STRING_DATA, &p))
                    {
                        ;
                    }
                }
                else if (!strcmp(field, "CLASS"))
                {
                    if (0 > MolfileV3000ReadField(&stmp, MOL_FMT_STRING_DATA, &p))
                    {
                        ;
                    }
                }
                else if (!strcmp(field, "SEQID"))
                {
                    if (0 > MolfileV3000ReadField(&itmp, MOL_FMT_INT_DATA, &p))
                    {
                        ;
                    }
                }

                if (failed)
                {
                    err = 4;
                    TREAT_ERR(err, 4, "Cannot interpret V3000 atom block key-value pair");
                    dotify_non_printable_chars(line);
                    AddErrorMessage(pStrErr, line);

                    if (!strcmp(line, SD_FMT_END_OF_DATA))
                    {
                        err = -abs(err);
                        break;
                    }
                    continue;
                }
            }
        } /* if ( NULL != ctab->atoms )  */
    } /* for ( i = 0; i < ctab->n_atoms; i++ )  */

    if (ctab->v3000->n_star_atoms)
    {
        AddErrorMessage(pStrErr, "V3000 star atoms ignored");
        ctab->n_atoms = ctab->v3000->n_non_star_atoms;
    }

    /* Check for proper finish */

    /*p = inchi_fgetsLf_V3000( line, inp_file );*/
    inchi_strbuf_reset(pin);

    nc = get_V3000_input_line_to_strbuf(pin, inp_file);

    if (nc < 1)
    {
        p = NULL;
    }
    else
    {
        p = line = pin->pStr;
    }
    if (!p || strcmp(p, "END ATOM"))
    {
        TREAT_ERR_AND_FIN(err, 1, err_fin, "Error: No V3000 Atom block end marker");
    }
    remove_one_lf(line);

err_fin:
    inchi_strbuf_close(pin);

    return err;
}

/****************************************************************************
 Read V3000 bonds
****************************************************************************/
int MolfileV3000ReadBondsBlock(MOL_FMT_CTAB *ctab,
                               INCHI_IOSTREAM *inp_file,
                               int err,
                               char *pStrErr)
{
    int i;
    char field[MOL_FMT_V3000_MAXFIELDLEN];
    int nc;
    char *p = NULL, *line = NULL;
    INCHI_IOSTREAM tmpin;
    INCHI_IOS_STRING *pin = &tmpin.s;

    if (!ctab->n_bonds)
    {
        return 0;
    }
    inchi_ios_init(&tmpin, INCHI_IOS_TYPE_STRING, NULL);

    /* Check for proper start */
    /*p = inchi_fgetsLf_V3000( line, inp_file );*/

    nc = get_V3000_input_line_to_strbuf(pin, inp_file);

    if (nc < 1)
    {
        p = NULL;
    }
    else
    {
        p = line = pin->pStr;
    }
    if (!p || strcmp(p, "BEGIN BOND"))
    {
        TREAT_ERR_AND_FIN(err, 1, err_fin, "Error: No V3000 Bond block start marker");
    }
    remove_one_lf(line);

    ctab->v3000->n_haptic_bonds = 0;
    ctab->v3000->n_non_haptic_bonds = 0;

    for (i = 0; i < ctab->n_bonds; i++)
    {
        int is_haptic = 0;

        /*p = inchi_fgetsLf_V3000( line, inp_file );*/
        inchi_strbuf_reset(pin);

        nc = get_V3000_input_line_to_strbuf(pin, inp_file);

        if (nc < 1)
        {
            p = NULL;
        }
        else
        {
            p = line = pin->pStr;
        }
        if (!p)
        {
            if (!err)
            {
                TREAT_ERR(err, 2, "Cannot read V3000 bond block line"); /* djb-rwth: addressing coverity ID #499565 -- TREAT_ERR properly used */
            }
            break;
        }
        remove_one_lf(line);

        if (err)
        {
            if (!strcmp(line, SD_FMT_END_OF_DATA))
            {
                err = -abs(err);
                break;
            }
            continue;
        }

        if (ctab->bonds)
        {
            int index, n_orig_at, len;
            short int atnum1 = -1, atnum2 = -1;
            char bond_type = 0, stereo = 0;
            int failed = 0;
            /* djb-rwth: removing redundant variables */

            n_orig_at = ctab->v3000->n_non_star_atoms + ctab->v3000->n_star_atoms;

            /* read positional parameters */
            if (0 > MolfileV3000ReadField(&index, MOL_FMT_INT_DATA, &p))
            {
                failed = 1;
            }
            else if (0 > MolfileV3000ReadField(&bond_type, MOL_FMT_CHAR_INT_DATA, &p))
            {
                failed = 1;
            }
            else if (0 > MolfileV3000ReadField(&atnum1, MOL_FMT_SHORT_INT_DATA, &p))
            {
                failed = 1;
            }
            else if (0 > MolfileV3000ReadField(&atnum2, MOL_FMT_SHORT_INT_DATA, &p))
            {
                failed = 1;
            }

            atnum1 = get_actual_atom_number(atnum1, n_orig_at,
                                            ctab->v3000->atom_index_orig,
                                            ctab->v3000->atom_index_fin);

            atnum2 = get_actual_atom_number(atnum2, n_orig_at,
                                            ctab->v3000->atom_index_orig,
                                            ctab->v3000->atom_index_fin);

            if ((atnum1 < 0) && (atnum2 < 0))
            {
                failed = 1;
            }
            /* djb-rwth: removing redundant code */

            if (failed)
            {

                if (!err)
                {
                    /* can't interpret bonds block line */
                    TREAT_ERR(err, 4, "Cannot interpret V3000 bond block line:");
                    dotify_non_printable_chars(line);
                    AddErrorMessage(pStrErr, line);
                }
                if (!strcmp(line, SD_FMT_END_OF_DATA))
                {
                    err = -abs(err);
                    break;
                }
            }

            /* TODO: treat new bond types  9 10 */
            /* read key-val pairs if any */
            while (p && (len = MolfileV3000ReadKeyword(field, &p)) > 0) /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
            {

                int itmp;
                char stmp[MOL_FMT_V3000_MAXFIELDLEN];
                failed = 0;

                if (!strcmp(field, "CFG"))
                {
                    if (0 > MolfileV3000ReadField(&stereo, MOL_FMT_CHAR_INT_DATA, &p))
                    {
                        failed = 1;
                    }
                    else
                    {
                        /*    adjust stereo to old convention for wedges which was:
                                0 = not stereo, 1 = Up,  4 = Either, 6 = Down
                            now:
                                0 = none (default), 1 = up, 2 = either, 3 = down
                        */
                        if (stereo == 2)
                        {
                            stereo = 4;
                        }
                        else if (stereo == 3)
                        {
                            stereo = 6;
                        }
                    }
                }
                else if (!strcmp(field, "TOPO"))
                {
                    if (0 > MolfileV3000ReadField(&itmp, MOL_FMT_INT_DATA, &p))
                    {
                        ; /* skip query-related stuff */
                    }
                }
                else if (!strcmp(field, "RXCTR"))
                {
                    if (0 > MolfileV3000ReadField(&itmp, MOL_FMT_INT_DATA, &p))
                    {
                        ; /* skip reaction-related stuff */
                    }
                }
                else if (!strcmp(field, "STBOX"))
                {
                    if (0 > MolfileV3000ReadField(&itmp, MOL_FMT_INT_DATA, &p))
                    {
                        ; /* skip for now */
                    }
                }
                else if (!strcmp(field, "ENDPTS"))
                {
                    int res, *num_list = NULL;
                    if (0 > MolfileV3000ReadHapticBond(ctab, &p, &num_list, pStrErr))
                    {
                        failed = 1;
                    }
                    else if (!num_list)
                    {
                        failed = 1;
                    }
                    else
                    {
                        int existent_atom = atnum1;
                        if (existent_atom < 0)
                        {
                            existent_atom = atnum2;
                        }
                        if (existent_atom < 0) /* should not be here */
                        {
                            failed = 1;
                        }
                        else
                        {
                            int k, nnum;
                            nnum = num_list[2];
                            num_list[1] = existent_atom;
                            for (k = 3; k < nnum; k++)
                            {
                                num_list[k] = get_actual_atom_number(num_list[k],
                                                                     n_orig_at,
                                                                     ctab->v3000->atom_index_orig,
                                                                     ctab->v3000->atom_index_fin);
                            }
                            res = NumLists_Append(ctab->v3000->haptic_bonds, num_list);
                            if (res < 0)
                            {
                                failed = 1;
                            }
                            else
                            {
                                is_haptic = 1;
                            }
                        }
                    }
                    /* djb-rwth: addressing coverity ID #499489 -- false positive as num_atoms allocated in MolfileV3000ReadHapticBond and returns a value in this block */
                }
                else if (!strcmp(field, "DISP"))
                {
                    if (0 > MolfileV3000ReadField(&stmp, MOL_FMT_STRING_DATA, &p))
                    {
                        ;
                    }
                }
                else if (!strcmp(field, "ATTACH"))
                {
                    if (0 > MolfileV3000ReadField(&stmp, MOL_FMT_STRING_DATA, &p))
                    {
                        ;
                    }
                }

                if (failed)
                {
                    if (!err)
                    {
                        /* can't interpret bonds block line */
                        TREAT_ERR(err, 4, "Cannot interpret V3000 bond block line:");
                        dotify_non_printable_chars(line);
                        AddErrorMessage(pStrErr, line);
                    }
                    if (!strcmp(line, SD_FMT_END_OF_DATA))
                    {
                        err = -abs(err);
                        break;
                    }
                }
            } /* while ( p && (len=MolfileV3000ReadKeyword(field, &p)) > 0 ) */

            if (is_haptic)
            {
                int ii = ctab->v3000->n_haptic_bonds;
                ctab->v3000->haptic_bonds->lists[ii][0] = bond_type;
                ctab->v3000->n_haptic_bonds++;
                continue;
            }
            else
            {
                int ii = ctab->v3000->n_non_haptic_bonds;
                ctab->bonds[ii].atnum1 = atnum1;
                ctab->bonds[ii].atnum2 = atnum2;
                ctab->bonds[ii].bond_type = bond_type;
                ctab->bonds[ii].bond_stereo = stereo;
                ctab->v3000->n_non_haptic_bonds++;
            }
        } /* if ctab->bonds */
    } /* for ( i = 0; i < ctab->n_bonds; i++ )  */

    if (ctab->v3000->n_haptic_bonds)
    {
        AddErrorMessage(pStrErr, "V3000 haptic bonds read/stored but ignored");
        ctab->n_bonds = ctab->v3000->n_non_haptic_bonds;
    }

    /* Check for proper finish */
    /*p = inchi_fgetsLf_V3000( line, inp_file );*/
    inchi_strbuf_reset(pin);

    nc = get_V3000_input_line_to_strbuf(pin, inp_file);

    if (nc < 1)
    {
        p = NULL;
    }
    else
    {
        p = line = pin->pStr;
    }
    if (!p || strcmp(p, "END BOND"))
    {
        TREAT_ERR_AND_FIN(err, 1, err_fin, "Error: No V3000 Bond block end marker");
    }
    remove_one_lf(line);

err_fin:

    inchi_ios_close(&tmpin); /* ricrogz: fixing memory leak */
    return err;
}

/****************************************************************************
 Convert atom index to the final consequitive atom number (starting from 1)
 Returns -1 for star atom or not found index
****************************************************************************/
int get_actual_atom_number(int index, int n, int *orig, int *fin)
{
    int i;
    for (i = 0; i < n; i++)
    {
        if (orig[i] == index)
        {
            return fin[i];
        }
    }

    return -1;
}

/****************************************************************************
 Read V3000 tail of CTab
****************************************************************************/
int MolfileV3000ReadTailOfCTAB(MOL_FMT_CTAB *ctab,
                               INCHI_IOSTREAM *inp_file,
                               int err,
                               char *pStrErr)
{
    int retcode = err;
    int nc;
    char *p = NULL, *line = NULL;
    INCHI_IOSTREAM tmpin;
    INCHI_IOS_STRING *pin = &tmpin.s;
    inchi_ios_init(&tmpin, INCHI_IOS_TYPE_STRING, NULL);

    /*p = inchi_fgetsLf_V3000( line, inp_file );*/

    nc = get_V3000_input_line_to_strbuf(pin, inp_file);

    if (nc < 1)
    {
        p = NULL;
    }
    else
    {
        p = line = pin->pStr;
    }
    remove_one_lf(line);

    if (p && !strcmp(p, "BEGIN SGROUP"))
    {
        retcode = MolfileV3000ReadSGroup(ctab, inp_file, retcode, pStrErr);
        if (retcode)
        {
            retcode += 70;
            TREAT_ERR_AND_FIN(retcode, 1, err_fin, pStrErr);
        }
        /*p = inchi_fgetsLf_V3000( line, inp_file );*/
        inchi_strbuf_reset(pin);
        nc = get_V3000_input_line_to_strbuf(pin, inp_file);
        if (nc < 1)
        {
            p = NULL;
        }
        else
        {
            p = line = pin->pStr;
            remove_one_lf(line);
        }
    }

    if (p && !strcmp(p, "BEGIN OBJ3D"))
    {
        retcode = MolfileV3000Read3DBlock(ctab, inp_file, retcode, pStrErr);
        if (retcode)
        {
            retcode += 70;
            TREAT_ERR_AND_FIN(retcode, 1, err_fin, pStrErr);
        }
        /*p = inchi_fgetsLf_V3000( line, inp_file );*/
        inchi_strbuf_reset(pin);

        nc = get_V3000_input_line_to_strbuf(pin, inp_file);
        if (nc < 1)
        {
            p = NULL;
        }
        else
        {
            p = line = pin->pStr;
            remove_one_lf(line);
        }
    }

    while (p && !strcmp(p, "LINKNODE"))
    {
        /* skip for now */
        /*p = inchi_fgetsLf_V3000( line, inp_file );*/
        inchi_strbuf_reset(pin);

        nc = get_V3000_input_line_to_strbuf(pin, inp_file);
        if (nc < 1)
        {
            p = NULL;
        }
        else
        {
            p = line = pin->pStr;
            remove_one_lf(line);
        }
    }

    /* Collections */
    while (p && !strcmp(p, "BEGIN COLLECTION"))
    {
        retcode = MolfileV3000ReadCollections(ctab, inp_file, retcode, pStrErr);
        if (retcode)
        {
            retcode += 70;
            TREAT_ERR_AND_FIN(retcode, 1, err_fin, pStrErr);
        }
        /*p = inchi_fgetsLf_V3000( line, inp_file );*/
        inchi_strbuf_reset(pin);

        nc = get_V3000_input_line_to_strbuf(pin, inp_file);

        if (nc < 1)
        {
            p = NULL;
        }
        else
        {
            p = line = pin->pStr;
            remove_one_lf(line);
        }
    }

    if (!p || strcmp(p, "END CTAB"))
    {
        TREAT_ERR_AND_FIN(err, 1, err_fin, "Error: No V3000 CTAB end marker");
    }

    remove_one_lf(line);

err_fin:
    inchi_strbuf_close(pin);

    return err;
}

/****************************************************************************
 Read haptic bond info
****************************************************************************/
int MolfileV3000ReadHapticBond(MOL_FMT_CTAB *ctab,
                               char **line_ptr,
                               int **num_list,
                               char *pStrErr)
{
    int nread = 0;
    char field[MOL_FMT_V3000_MAXFIELDLEN];
    const int max_field_len = sizeof(field);
    char *p_end;
    int i, nnum = 0;

    *num_list = NULL;

    memset(field, 0, max_field_len); /* djb-rwth: memset_s C11/Annex K variant? */

    nread = read_upto_delim(line_ptr, field, max_field_len, "1234567890 \t\n\v\f\r"); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    if (strcmp(field, "("))
    {
        return -1;
    }

    nread = read_upto_delim(line_ptr, field, max_field_len, " \t\n\v\f\r"); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */

    nnum = strtol(field, &p_end, 10);

    if (p_end == field)
    {
        return -1; /* paranoia */
    }
    if (nnum < 0)
    {
        return -1;
    }

    *num_list = (int *)inchi_calloc((long long)nnum + 3, sizeof(int)); /* djb-rwth: cast operator added */

    if (!*num_list)
    {
        nread = -1;
        goto ret;
    }

    (*num_list)[0] = -1; /* will be bond type, to be filled by caller */
    (*num_list)[1] = -1; /* will be atom number, to be filled by caller */
    (*num_list)[2] = nnum;

    for (i = 3; i < nnum + 3; i++)
    {
        if (0 > MolfileV3000ReadField(&((*num_list)[i]), MOL_FMT_INT_DATA, line_ptr))
        {
            nread = -1;
            goto ret;
        }
    }

    /* ')' should have been consumed  by strtol */

    /* check for ATTACH=ALL */

    nread = read_upto_delim(line_ptr, field, max_field_len, " \t\n\v\f\r");
    if (nread > 0)
    {
        if (strcmp(field, "ATTACH=ALL"))
        {
            nread = -1;
            goto ret;
        }
    }

ret:
    if (nread < 0)
    {
        if (*num_list)
        {
            inchi_free(*num_list);
            *num_list = NULL;
        }
    }

    return nread;
}

/****************************************************************************
 Read V3000 stereo collection
****************************************************************************/
int MolfileV3000ReadStereoCollection(MOL_FMT_CTAB *ctab,
                                     char **line_ptr,
                                     int **num_list,
                                     char *pStrErr)
{
    int nread = 0;
    char field[MOL_FMT_V3000_MAXFIELDLEN];
    const int max_field_len = sizeof(field);
    char *p_end;
    int i, nnum = 0;

    *num_list = NULL;

    memset(field, 0, max_field_len); /* djb-rwth: memset_s C11/Annex K variant? */

    nread = read_upto_delim(line_ptr, field, max_field_len, "1234567890 \t\n\v\f\r"); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    if (strcmp(field, "("))
    {
        return -1;
    }

    nread = read_upto_delim(line_ptr, field, max_field_len, " \t\n\v\f\r");

    nnum = strtol(field, &p_end, 10);

    if (p_end == field)
    {
        return -1; /* paranoia */
    }
    if (nnum < 0)
    {
        return -1;
    }

    *num_list = (int *)inchi_calloc((long long)nnum + 3, sizeof(int)); /* djb-rwth: cast operator added */

    if (!*num_list)
    {
        nread = -1;
        goto ret;
    }

    (*num_list)[0] = -1; /* reserved, may be filled by caller */
    (*num_list)[1] = nnum;

    for (i = 2; i < nnum + 2; i++)
    {
        if (0 > MolfileV3000ReadField(&((*num_list)[i]), MOL_FMT_INT_DATA, line_ptr))
        {
            nread = -1;
            goto ret;
        }
    }

    /* ')' should have been consumed  by strtol */

ret:
    if (nread < 0)
    {
        if (*num_list)
        {
            inchi_free(*num_list);
            *num_list = NULL;
        }
    }

    return nread;
}

/****************************************************************************
    Returns -1 @ error
****************************************************************************/
int get_V3000_input_line_to_strbuf(INCHI_IOS_STRING *buf,
                                   INCHI_IOSTREAM *inp_stream)
{
    const int prefix_len = 7; /* "M  V30 " */
    int old_used, crlf2lf = 1, preserve_lf = 0;

    inchi_strbuf_reset(buf);

    old_used = buf->nUsedLength;
    while (1)
    {
        inchi_strbuf_addline(buf, inp_stream, crlf2lf, preserve_lf);

        if (buf->nUsedLength - old_used < 8)
        {
            return -1;
        }
        if (strncmp(buf->pStr + old_used, "M  V30 ", prefix_len))
        {
            return -1;
        }

        memmove((void *)(buf->pStr + old_used), (void *)(buf->pStr + old_used + prefix_len), (long long)buf->nUsedLength - (long long)old_used - (long long)prefix_len + 1); /* djb-rwth: cast operators added */ /* ricrogz: fixing memory overflow error */
        buf->nUsedLength -= prefix_len;

        if (buf->pStr[buf->nUsedLength - 1] != '-')
        {
            break;
        }
        buf->pStr[--buf->nUsedLength] = '\0';

        old_used = buf->nUsedLength;
    }

    remove_trailing_spaces(buf->pStr);
    buf->nUsedLength = strlen(buf->pStr);

    return buf->nUsedLength;
}

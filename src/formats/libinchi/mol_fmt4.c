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
#include "ichimain.h"

#include "bcf_s.h"

/*
    SDFile related procedures

*/

#define ALIASED_AT(i) (0 < NUM_ISO_H(at, i))
#define IS_DEUTERIUM(i) (!strcmp(at[i].elname, "D") || at[i].iso_atw_diff == 2 && !strcmp(at[i].elname, "H"))
#define IS_TRITIUM(i)   (!strcmp(at[i].elname, "T") || at[i].iso_atw_diff == 3 && !strcmp(at[i].elname, "H"))

#define ABNORMAL_ISO(i) (at[i].iso_atw_diff == 1 || at[i].iso_atw_diff < -3 || at[i].iso_atw_diff > 5)
#define ABNORMAL_CHG(i) (abs(at[i].charge) > 3)
#define ABNORMAL_RAD(i) (RADICAL_SINGLET <= at[i].radical && at[i].radical <= RADICAL_TRIPLET)

#define ANY_ISO(i, X)   ((X)? (at[i].iso_atw_diff && !IS_DEUTERIUM(i) && !IS_TRITIUM(i)) :\
                          (at[i].iso_atw_diff ||  IS_DEUTERIUM(i) ||  IS_TRITIUM(i)))
#define ANY_CHG(i)      (0 != at[i].charge)
#define ANY_RAD(i)      (RADICAL_SINGLET <= at[i].radical && at[i].radical <= RADICAL_TRIPLET)

#define NORMAL_ISO(i, X)   (ANY_ISO(i, X) && !ABNORMAL_ISO(i))

/* needs additional M  CHG. M  RAD, M  ISO line */
/* due to ISIS/Draw feature always include M  RAD for any radical */
#define ABNORMAL_AT(i) (at[i].radical || abs(at[i].charge) > 3 || \
                        ABNORMAL_ISO(i))

/* always add M  ISO, M  RAD, M  CHG; Except: (bAtomsDT && D or T) */
#define ADD_LINE_AT(i) (at[i].charge ||  \
                        at[i].radical || \
                        at[i].iso_atw_diff && (bAtomsDT ? (at[i].iso_atw_diff != 1 || strcmp(at[i].elname, "H")) : 1))

/* Local */

static const char sdf_data_hdr_name[] = "NAME";
static const char sdf_data_hdr_comm[] = "COMMENT";

enum
{
    SDF_START,
    SDF_DATA_HEADER,
    SDF_DATA_HEADER_NAME,
    SDF_DATA_HEADER_COMMENT,
    SDF_DATA_HEADER_CAS,
    SDF_DATA_HEADER_USER,
    SDF_DATA_LINE,
    SD_FMT_END_OF_DATA_ITEM,
    SDF_EMPTY_LINE,
    SD_FMT_END_OF_DATA_BLOCK
};

int OrigAtData_WriteToSDfileHeaderAndCountThings(const ORIG_ATOM_DATA *inp_at_data,
                                                 INCHI_IOSTREAM *fcb,
                                                 const char *name,
                                                 const char *comment,
                                                 int bChiralFlag,
                                                 int bAtomsDT,
                                                 const char *szLabel,
                                                 const char *szValue,
                                                 int *nNumAliasLines,
                                                 int *nNumChargeLines,
                                                 int *nNumRadicalLines,
                                                 int *nNumIsoLines,
                                                 int *nNumAddLines,
                                                 int *num_bonds);
int OrigAtData_WriteToSDfileAtomsBlock(const ORIG_ATOM_DATA *inp_at_data,
                                       INCHI_IOSTREAM *fcb,
                                       const char *name,
                                       const char *comment,
                                       int bAtomsDT,
                                       const char *szLabel,
                                       const char *szValue);

int OrigAtData_WriteToSDfileBondsBlock(const ORIG_ATOM_DATA *inp_at_data,
                                       INCHI_IOSTREAM *fcb,
                                       const char *name,
                                       const char *comment,
                                       const char *szLabel,
                                       const char *szValue,
                                       INT_ARRAY *written_bond_ends);

int OrigAtData_WriteToSDfileAdditionalLines(const ORIG_ATOM_DATA *inp_at_data,
                                            INCHI_IOSTREAM *fcb,
                                            const char *name,
                                            const char *comment,
                                            int bAtomsDT,
                                            const char *szLabel,
                                            const char *szValue,
                                            int nNumAliasLines,
                                            int nNumChargeLines,
                                            int nNumRadicalLines,
                                            int nNumIsoLines,
                                            INT_ARRAY *written_bond_ends);

int OrigAtData_WriteToSDfilePolymerData(const ORIG_ATOM_DATA *inp_at_data,
                                        INCHI_IOSTREAM *fcb,
                                        const char *name,
                                        const char *comment,
                                        const char *szLabel,
                                        const char *szValue,
                                        INT_ARRAY *written_bond_ends);

/****************************************************************************
 Skip extra data ( != Molfile) which SDF contains
****************************************************************************/
int SDFileSkipExtraData(INCHI_IOSTREAM *inp_file,
                        unsigned long *CAS_num,
                        char *comment,
                        int lcomment,
                        char *name,
                        int lname,
                        int prev_err,
                        const char *pSdfLabel,
                        char *pSdfValue,
                        char *pStrErr,
                        int bNoWarnings)
{
    char *p = NULL;
    char line[MOL_FMT_INPLINELEN];
    const int line_len = sizeof( line );
    int n_blank_lines = 0, n_lines = 0;
    int current_state = SDF_START;
    int err = 0;
    int wait_for_CAS = 0;
    int CAS_num_is_user = 0;
    int wait_for_name = name && lname > 0 && !name[0];
    int wait_for_comment = comment && lcomment > 0 && !comment[0];
    int wait_for_user = pSdfLabel && pSdfLabel[0] && pSdfValue;

    if (CAS_num != NULL)
    {
        wait_for_CAS = 1;
        *CAS_num = 0LU;
        CAS_num_is_user = (wait_for_user && !inchi_memicmp(pSdfLabel, "CAS", 3));
    }

    while (!err &&
           current_state != SD_FMT_END_OF_DATA_BLOCK &&
           NULL != (p = inchi_fgetsLf(line, line_len, inp_file)))
    {
        if (!n_lines && !memcmp(line, "M  END", 6))
        {
            /*  allow subtle errors */
            continue;
        }

        n_lines++;
        remove_trailing_spaces(line);

        if (line[MOL_FMT_MAXLINELEN])
        {
            if (current_state != SDF_DATA_HEADER &&
                current_state != SDF_DATA_LINE &&
                current_state != SDF_DATA_HEADER_NAME &&
                current_state != SDF_DATA_HEADER_USER &&
                current_state != SDF_DATA_HEADER_COMMENT)
            {
                line[MOL_FMT_MAXLINELEN] = '\0';
                if (!prev_err)
                {
                    TREAT_ERR(err, 0, "Too long SData line truncated");
                }
            }
            else
            {
                /* allow long lines in SDF data. 9-29-00 DCh */
                line[MOL_FMT_MAXLINELEN] = '\0';
            }
        }

        n_blank_lines += (*line == '\0');

        switch (current_state)
        {
            case SDF_START:
            case SD_FMT_END_OF_DATA_ITEM:
            case SDF_EMPTY_LINE: /* Added 9-25-97 DCh */

                if (!strcmp(line, SD_FMT_END_OF_DATA))
                {
                    current_state = SD_FMT_END_OF_DATA_BLOCK;
                }
                else if ('>' == *line)
                {
                    current_state = (wait_for_name || wait_for_comment || wait_for_CAS || wait_for_user) ? SDFileIdentifyLabel( line, pSdfLabel ) : SDF_DATA_HEADER;
                }
                else if (*line == '\0')
                {
                    /* Added 9-25-97 DCh */
                    /* Relax the strictness: Allow more than 1 empty line. */
                    current_state = SDF_EMPTY_LINE;
                }
                else if (!prev_err)
                {
                    TREAT_ERR(err, 3, "Unexpected SData header line:"); /* djb-rwth: addressing coverity ID #499557 -- TREAT_ERR properly used */
                    dotify_non_printable_chars(line);
                    AddErrorMessage(pStrErr, line);
                    /* unexpected contents of data header line */
                }
                else
                {
                    err = 3;
                }
                break;

            case SDF_DATA_HEADER_NAME:

                if (wait_for_name && 0 < normalize_string(line))
                {
                    wait_for_name = 0;
                    mystrncpy(name, line, lname);
                }
                goto got_data_line;

            case SDF_DATA_HEADER_COMMENT:

                if (wait_for_comment && 0 < normalize_string(line))
                {
                    wait_for_comment = 0;
                    mystrncpy(comment, line, lcomment);
                }
                goto got_data_line;

            case SDF_DATA_HEADER_USER:

                if (wait_for_user && 0 < normalize_string(line))
                {
                    wait_for_user = 0;
                    mystrncpy(pSdfValue, line, MAX_SDF_VALUE + 1);

                    if (CAS_num_is_user && wait_for_CAS)
                    {
                        *CAS_num = SDFileExtractCASNo(line);
                        wait_for_CAS = (0LU == *CAS_num);
                    }
                }
                goto got_data_line;

            case SDF_DATA_HEADER_CAS:

                if (wait_for_CAS && 0 < normalize_string(line))
                {
                    *CAS_num = SDFileExtractCASNo(line);
                    wait_for_CAS = (0LU == *CAS_num);
                }
                goto got_data_line;

            case SDF_DATA_HEADER:
            case SDF_DATA_LINE:

            got_data_line:
                current_state = *line ? SDF_DATA_LINE : SD_FMT_END_OF_DATA_ITEM;
                break;
        }
    }

    if (!err && SD_FMT_END_OF_DATA_BLOCK != current_state && NULL == p)
    {
        ; /* err = 4; */ /* unexpected end of file: missing $$$$ */
    }

    else if (err && (n_blank_lines == n_lines && *line == '\0'))
    {
        /* empty lines -- do not know when this can happen */
        err = 5;
    }

    if (err && err != 5 && current_state != SD_FMT_END_OF_DATA_BLOCK && p)
    {
        /*  bypass up to $$$$ */
        while ((p = inchi_fgetsLf(line, line_len, inp_file)) &&
               memcmp(line, SD_FMT_END_OF_DATA, 4))
        {
            ;
        }
        if (p)
        {
            /*  arrived to $$$$; non-fatal */
            err = 9;
            if (!bNoWarnings)
            {
                WarningMessage(pStrErr, "Bypassing to next structure");
            }
        }
    }

    return err;
}

/****************************************************************************/
int SDFileIdentifyLabel(char *inp_line, const char *pSdfLabel)
{
    char line[MOL_FMT_MAXLINELEN];
    char *p, *q;
    int i, j, len, tmp1 = 0, tmp2 = 0, cnd = 0;

    if ((p = strchr(inp_line, '<')) &&
        (q = strchr(p, '>')) &&
        (len = q - p - 1) > 0 && len < (int)sizeof(line))
    {
        memcpy(line, p + 1, len);
        line[len] = '\0';

        for (i = 0; (i < len) && (cnd == 0); i++)
        {
            if (isspace(UCINT line[i]))
            {
                tmp1 = i;
                cnd = 1;
            }
        }

        cnd = 0;

        for (j = len - 1; (j >= tmp1) && (cnd == 0); j--)
        {
            if (isspace(UCINT line[j]))
            {
                tmp2 = j;
                cnd = 1;
            }
        }

        len = tmp2 - tmp1 + 1;
        p = line + tmp1;

        if (pSdfLabel && pSdfLabel[0] && len == (int)strlen(pSdfLabel) && !inchi_memicmp(p, pSdfLabel, len))
        {
            return SDF_DATA_HEADER_USER;
        }

        if (len == sizeof(sdf_data_hdr_name) - 1 && !inchi_memicmp(p, sdf_data_hdr_name, len))
        {
            return SDF_DATA_HEADER_NAME;
        }

        if (len == sizeof(sdf_data_hdr_comm) - 1 && !inchi_memicmp(p, sdf_data_hdr_comm, len))
        {
            return SDF_DATA_HEADER_COMMENT;
        }

        if (!inchi_memicmp(p, "CAS", 3))
        {
            return SDF_DATA_HEADER_CAS;
        }
    }

    return SDF_DATA_HEADER;
}

/****************************************************************************/
unsigned long SDFileExtractCASNo(char *line)
{
    int i, j;

    i = line[0] == '-' ? 1 : 0;

    for (j = i; line[i]; i++)
    {
        if (isdigit(UCINT line[i]))
        {
            line[j++] = line[i];
        }
        else if (line[i] != '-')
        {
            break;
        }
    }

    line[j] = '\0';

    return strtoul(line, NULL, 10);
}

/****************************************************************************
 NUM_LISTS - dynamically growing array of int lists
****************************************************************************/
int NumLists_Alloc(NUM_LISTS *num_lists, int nlists)
{
    if (num_lists)
    {
        if ((num_lists->lists = (int **)inchi_calloc(nlists, sizeof(int *)))) /* djb-rwth: addressing LLVM warning */
        {
            num_lists->increment =
                num_lists->allocated = nlists;
            return 0; /*  ok */
        }
    }

    return -1; /*  error */
}

/****************************************************************************/
int NumLists_ReAlloc(NUM_LISTS *num_lists)
{
    if (num_lists)
    {
        if (num_lists->lists && num_lists->allocated > 0 && num_lists->increment > 0)
        {
            void *p = num_lists->lists;
            if ((num_lists->lists =
                     (int **)inchi_calloc((long long)num_lists->allocated + (long long)num_lists->increment, sizeof(int *)))) /* djb-rwth: cast operators added; addressing LLVM warning */
            {
                memcpy(num_lists->lists, p, num_lists->used * sizeof(num_lists->lists[0]));
                inchi_free(p);
                num_lists->allocated += num_lists->increment;
                return 0; /*  ok */
            }
        }
    }

    return -1; /*  error */
}

/****************************************************************************/
int NumLists_Append(NUM_LISTS *num_lists, int *list)
{
    if (num_lists)
    {
        if (num_lists->used + 1 > num_lists->allocated)
        {
            /* need to expand buffer */
            if (NumLists_ReAlloc(num_lists))
            {
                return -1; /*  error */
            }
        }
        num_lists->lists[num_lists->used++] = list;
        return 0;
    }

    return -1;
}

/****************************************************************************/
void NumLists_Free(NUM_LISTS *num_lists)
{
    if (num_lists)
    {
        int i;
        for (i = 0; i < num_lists->used; i++)
            inchi_free(num_lists->lists[i]); /* djb-rwth: unresolved issue -- revision required? -- false positive as this function just does the clean-up job */
        inchi_free(num_lists->lists);
        memset(num_lists, 0, sizeof(*num_lists)); /* djb-rwth: memset_s C11/Annex K variant? */
    }
}

/****************************************************************************
 INT_ARRAY - dynamically growing array of int
****************************************************************************/

/****************************************************************************
 Allocate new array, return 0 if OK, -1 otherwise
****************************************************************************/
int IntArray_Alloc(INT_ARRAY *items, int nitems)
{
    if ((items->item = (int *)inchi_calloc(nitems, sizeof(int)))) /* djb-rwth: addressing LLVM warning */
    {
        items->increment = items->allocated = nitems;
        items->used = 0;
        return 0;
    }

    return -1;
}

/****************************************************************************
 Expand array, return 0 if OK, -1 otherwise
****************************************************************************/
int IntArray_ReAlloc(INT_ARRAY *items)
{
    if (items)
    {
        if (items->item && items->allocated > 0 && items->increment > 0)
        {
            void *p = items->item;
            if ((items->item =
                     (int *)inchi_calloc((long long)items->allocated + (long long)items->increment, sizeof(items->item[0])))) /* djb-rwth: cast operators added; addressing LLVM warning */
            {
                memcpy(items->item, p, items->used * sizeof(items->item[0]));
                inchi_free(p);
                items->allocated += items->increment;
                return 0;
            }
            inchi_free(p);
        }
    }

    return -1;
}

/****************************************************************************
 Push new item to the end of array
****************************************************************************/
int IntArray_Append(INT_ARRAY *items, int new_item)
{
    if (items)
    {
        if (items->used + 1 > items->allocated)
        {
            /* need to expand buffer */
            if (IntArray_ReAlloc(items))
            {
                return -1;
            }
        }
        items->item[items->used++] = new_item;
        return 0;
    }

    return -1;
}

/****************************************************************************
Push new item to the end of array only if it is absent there
****************************************************************************/
int IntArray_AppendIfAbsent(INT_ARRAY *items, int new_item)
{
    if (!is_in_the_ilist(items->item, new_item, items->used))
    {
        return IntArray_Append(items, new_item);
    }
    return 0;
}

/****************************************************************************/
void IntArray_DebugPrint(INT_ARRAY *items)
{
    if (items)
    {
        int i;
        if (items->used > 0)
        {
            for (i = 0; i < items->used - 1; i++)
            {
                ITRACE_("%-d, ", items->item[i]);
            }
            ITRACE_("%-d\n", items->item[items->used - 1]);
        }
        else
        {
            ; /*ITRACE_( "[None]\n");*/
        }
    }
}

/****************************************************************************/
void IntArray_Reset(INT_ARRAY *items)
{
    items->used = 0;
    return;
}

/****************************************************************************/
void IntArray_Free(INT_ARRAY *items)
{
    if (items)
    {
        if (items->item)
        {
            inchi_free(items->item);
        }
        memset(items, 0, sizeof(*items)); /* djb-rwth: memset_s C11/Annex K variant? */
    }
    return;
}

/****************************************************************************
    MOL_FMT_SGROUPS - dynamically growing array of pointers to SGroups
****************************************************************************/

/*
    SGroup
*/

/****************************************************************************
 Allocate new array Sgroup, return 0 if OK, -1 otherwise
****************************************************************************/
int MolFmtSgroup_Create(MOL_FMT_SGROUP **sgroup, int id, int type)
{
    *sgroup = (MOL_FMT_SGROUP *)inchi_calloc(1, sizeof(MOL_FMT_SGROUP));
    if (*sgroup)
    {
        if (IntArray_Alloc(&((*sgroup)->alist), 8) ||
            IntArray_Alloc(&((*sgroup)->blist), 8))
        {
            MolFmtSgroup_Free(*sgroup);
            return -1;
        }
        (*sgroup)->id = id;
        (*sgroup)->type = type;

        (*sgroup)->subtype = 0;
        (*sgroup)->conn = 0;
        (*sgroup)->label = 0;

        return 0;
    }
    return -1;
}

/****************************************************************************/
void MolFmtSgroup_Free(MOL_FMT_SGROUP *sgroup)
{
    if (sgroup)
    {
        IntArray_Free(&(sgroup->alist));
        IntArray_Free(&(sgroup->blist));
        inchi_free(sgroup);
    }
}

/*
    SGroups
*/

/****************************************************************************
 Allocate new array of Sgroups, return 0 if OK, -1 otherwise
****************************************************************************/
int MolFmtSgroups_Alloc(MOL_FMT_SGROUPS *sgroups, int nsgroups)
{
    if (sgroups)
    {
        if (NULL != (sgroups->group = (MOL_FMT_SGROUP **)inchi_calloc(nsgroups, sizeof(MOL_FMT_SGROUP *))))
        {
            /* ITRACE_( "\nAllocated sgroups->group at %-p \n", sgroups->group ); */
            sgroups->increment = sgroups->allocated = nsgroups;
            return 0;
        }
    }

    return -1;
}

/****************************************************************************
 Expand array of Sgroups, return 0 if OK, -1 otherwise
****************************************************************************/
int MolFmtSgroups_ReAlloc(MOL_FMT_SGROUPS *sgroups)
{
    if (sgroups)
    {
        if (sgroups->group && sgroups->allocated > 0 && sgroups->increment > 0)
        {
            void *p = sgroups->group;
            if ((sgroups->group = (MOL_FMT_SGROUP **)inchi_calloc( (long long)sgroups->allocated + (long long)sgroups->increment,
                sizeof( sgroups->group[0] ) ))) /* djb-rwth: cast operators added; addressing LLVM warning */
            {
                memcpy(sgroups->group, p, sgroups->used * sizeof(sgroups->group[0]));
                inchi_free(p);
                sgroups->allocated += sgroups->increment;
                return 0; /*  ok */
            }
        }
    }

    return -1;
}

/****************************************************************************/
int MolFmtSgroups_Append(MOL_FMT_SGROUPS *sgroups, int id, int type)
{
    if (sgroups)
    {
        /* Make new Sgroup */
        MOL_FMT_SGROUP *sgroup = NULL;
        if (0 != MolFmtSgroup_Create(&sgroup, id, type))
        {
            return -1;
        }
        /* Add new created Sgroup to Sgroups */
        if (sgroups->used + 1 > sgroups->allocated)
        {
            /* expand buffer */
            if (MolFmtSgroups_ReAlloc(sgroups))
            {
                MolFmtSgroup_Free(sgroup); /* djb-rwth: avoiding memory leak */
                return -1; /*  no RAM */
            }
        }
        sgroups->group[sgroups->used++] = sgroup;

        /*
        {
        int num = sgroups->used-1;
        printf("\nCreated/added Sgroup: id=%-d ( num in Sgroups=%-d ) of type=%-d \n", sgroups->group[num]->id, num, sgroups->group[num]->type );
        }
        */

        return 0;
    }

    return -1;
}

/****************************************************************************/
void MolFmtSgroups_Free(MOL_FMT_SGROUPS *sgroups)
{
    if (sgroups)
    {
        int i;
        for (i = 0; i < sgroups->used; i++)
        {
            MolFmtSgroup_Free(sgroups->group[i]);
        }

        /* ITRACE_( "\nAbout to free sgroups->group at %-p\n", sgroups->group ); */
        inchi_free(sgroups->group);

        memset(sgroups, 0, sizeof(MOL_FMT_SGROUPS)); /* djb-rwth: memset_s C11/Annex K variant? */
    }
}

/****************************************************************************/
int MolFmtSgroups_GetIndexBySgroupId(int id, MOL_FMT_SGROUPS *sgroups)
{
    int i;
    for (i = 0; i < sgroups->used; i++)
    {
        if (sgroups->group[i]->id == id)
        {
            return i;
        }
    }
    return -1;
}

/****************************************************************************/
int OrigAtData_WriteToSDfile(const ORIG_ATOM_DATA *inp_at_data,
                             INCHI_IOSTREAM *fcb,
                             const char *name,
                             const char *comment,
                             int bChiralFlag,
                             int bAtomsDT,
                             const char *szLabel,
                             const char *szValue)
{
    int num_bonds = 0, nNumAddLines = 0, nNumIsoLines = 0, nNumChargeLines = 0,
        nNumRadicalLines = 0, nNumAliasLines = 0, ret = 0;
    INT_ARRAY written_bond_ends;

    OrigAtData_WriteToSDfileHeaderAndCountThings((ORIG_ATOM_DATA *)inp_at_data,
                                                 fcb, name, comment,
                                                 bChiralFlag, bAtomsDT,
                                                 szLabel, szValue,
                                                 &nNumAliasLines,
                                                 &nNumChargeLines,
                                                 &nNumRadicalLines,
                                                 &nNumIsoLines,
                                                 &nNumAddLines,
                                                 &num_bonds);

    if (IntArray_Alloc(&written_bond_ends, num_bonds ? num_bonds : 255))
    {
        ret = _IS_ERROR;
        goto exit_function;
    }

    OrigAtData_WriteToSDfileAtomsBlock(inp_at_data, fcb, name, comment,
                                       bAtomsDT, szLabel, szValue);

    OrigAtData_WriteToSDfileBondsBlock(inp_at_data, fcb, name, comment,
                                       szLabel, szValue, &written_bond_ends);

    if (nNumAddLines)
    {
        OrigAtData_WriteToSDfileAdditionalLines(inp_at_data, fcb, name, comment,
                                                bAtomsDT, szLabel, szValue,
                                                nNumAliasLines, nNumChargeLines,
                                                nNumRadicalLines, nNumIsoLines,
                                                &written_bond_ends);
    }

    /* Add field with label/ID if applicable and mark the end of record */
    if (szValue && szValue[0])
    {
        if (szLabel && szLabel[0])
        {
            inchi_ios_print_nodisplay(fcb, "> <%s>\n", szLabel);
        }
        else
        {
            inchi_ios_print_nodisplay(fcb, "> <ID>\n");
        }
        inchi_ios_print_nodisplay(fcb, " %s\n\n", szValue);
    }
    inchi_ios_print_nodisplay(fcb, "$$$$\n");

exit_function:
    IntArray_Free(&written_bond_ends);

    return ret;
}

/****************************************************************************
 OrigAtData : Write To SDfile : Atoms Block
****************************************************************************/
int OrigAtData_WriteToSDfileHeaderAndCountThings(const ORIG_ATOM_DATA *inp_at_data,
                                                 INCHI_IOSTREAM *fcb,
                                                 const char *name,
                                                 const char *comment,
                                                 int bChiralFlag,
                                                 int bAtomsDT,
                                                 const char *szLabel,
                                                 const char *szValue,
                                                 int *nNumAliasLines,
                                                 int *nNumChargeLines,
                                                 int *nNumRadicalLines,
                                                 int *nNumIsoLines,
                                                 int *nNumAddLines,
                                                 int *num_bonds)
{
    int i, ret = 0;
    int bAtomNeedsAlias,
        nNumNecessaryIsoLines = 0,
        nNumNecessaryChgLines = 0,
        nNumNecessaryRadLines = 0;
    int num_atoms = inp_at_data->num_inp_atoms;
    int bV2000 = SDF_OUTPUT_V2000;
    const inp_ATOM *at = inp_at_data->at;

    {
        char strLocName[82];
        memset(strLocName, 0, sizeof(strLocName)); /* djb-rwth: memset_s C11/Annex K variant? */
        if (name && *name)
        {
            strncpy(strLocName, name, 80);
        }
        inchi_ios_print_nodisplay(fcb, "%s\n", strLocName);
    }

    /**********************************************************************/
    /**                                                                  **/
    /** Important: Atoms with alias cannot have charge, radical          **/
    /**            isotope differences are allowed                       **/
    /**                                                                  **/
    /**            Atoms with alias cannot be abnormal.                  **/
    /**                                                                  **/
    /** Abnormal atoms are atoms which need M  CHG, M RAD, M  ISO        **/
    /**                                                                  **/
    /** Output aliased atoms if they have implicit D or T                **/
    /**                                                                  **/
    /**********************************************************************/

    /*                                    F10.5     F12.5       I6
                     IIPPPPPPPPMMDDYYHHmmddSSssssssssssEEEEEEEEEEEERRRRRR
    inchi_ios_eprint( fcb,"NISTTRANHP09089809272D 1   1.0         0.0    %6ld\n", lEpa);*/
    /*^^^
    inchi_ios_print_nodisplay( fcb,"  %s v%s SDfile Output                       \n", INCHI_NAME, INCHI_VERSION);

    Changed 01/10/2009 to conform CTFile specification (by Symyx request)*/

    inchi_ios_print_nodisplay(fcb,
                              /*   IIPPPPPPPPMMDDYYHHmmddSSssssssssssEEEEEEEEEEEERRRRRR*/
                              "  InChIV10                                     \n");
    /*y_fprintf(fcb, "  -CPSS-  1213981200n\n");*/

    {
        char strLocName[82];

        memset(strLocName, 0, sizeof(strLocName)); /* djb-rwth: memset_s C11/Annex K variant? */
        if (comment && *comment)
        {
            strncpy(strLocName, comment, 80);
        }
        inchi_ios_print_nodisplay(fcb, "%s\n", strLocName);
    }

    *num_bonds = 0;
    for (i = 0; i < num_atoms; i++)
    {
        (*num_bonds) += at[i].valence;
    }
    (*num_bonds) /= 2;

    /*find if we need "M  CHG" and "M  RAD"*/
    for (i = 0; i < num_atoms; i++)
    {
        if ((bAtomNeedsAlias = ALIASED_AT(i))) /* djb-rwth: addressing LLVM warning */
        {
            /* has isotopic implicit D or T; ignoring pure 1H */
            (*nNumAliasLines) += 2 * bAtomNeedsAlias;
        }
        else
        {
            /* abnormal means atom needs CHG, RAD, or ISO entry */
            /* nNumAddLines    += ABNORMAL_AT(i); */
            /* nNumIso         += ( 0 == strcmp( at[i].elname, "D" ) || ( 0 == strcmp( at[i].elname, "T" ) || at[i].iso_atw_diff ) ); */
            /* nNumAddIso      += at[i].iso_atw_diff && (at[i].iso_atw_diff == 1 || at[i].iso_atw_diff < -3 || at[i].iso_atw_diff > 5 ); */
            nNumNecessaryIsoLines += ABNORMAL_ISO(i);
            nNumNecessaryChgLines += ABNORMAL_CHG(i);
            nNumNecessaryRadLines += ABNORMAL_RAD(i);
            (*nNumIsoLines) += ANY_ISO(i, bAtomsDT);
            (*nNumChargeLines) += ANY_CHG(i);
            (*nNumRadicalLines) += ANY_RAD(i);
        }
    }

    *nNumChargeLines = (*nNumChargeLines + 7) / 8;
    *nNumRadicalLines = (*nNumRadicalLines + 7) / 8;
    *nNumIsoLines = (*nNumIsoLines + 7) / 8;

    if (!bV2000)
    {
        if (!nNumNecessaryRadLines && !nNumNecessaryChgLines)
        {
            *nNumRadicalLines = 0;
            *nNumChargeLines = 0;
        }
        if (!nNumNecessaryIsoLines)
        {
            *nNumIsoLines = 0;
        }
    }

    /* recalculate number of added lines */
    *nNumAddLines = *nNumChargeLines + *nNumRadicalLines + *nNumIsoLines + *nNumAliasLines; /* 1 for M  END*/

    if (*nNumAddLines || bV2000)
    {
        *nNumAddLines += 1; /* add 1 for "M  END" line*/
    }

    /*                         aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv                                      */

    inchi_ios_print_nodisplay(fcb, "%3d%3d  0  0%3d  0  0  0  0  0%3d%s\n",
                              num_atoms, *num_bonds, bChiralFlag ? 1 : 0, *nNumAddLines, *nNumAddLines ? " V2000" : "");

    return ret;
}

/****************************************************************************
 OrigAtData : Write To SDfile : Atoms Block
****************************************************************************/
int OrigAtData_WriteToSDfileAtomsBlock(const ORIG_ATOM_DATA *inp_at_data,
                                        INCHI_IOSTREAM       *fcb,
                                        const char           *name,
                                        const char           *comment,
                                        int                  bAtomsDT,
                                        const char           *szLabel,
                                        const char           *szValue)
{
    int i, ret = 0;
    int bAtomNeedsAlias;
    /* djb-rwth: removing redundant variables */
    int num_atoms = inp_at_data->num_inp_atoms;
    const inp_ATOM *at = inp_at_data->at;
    double x, y, z;

    for (i = 0; i < num_atoms; i++)
    {
        char elname[ATOM_EL_LEN] = "\0\0\0\0\0";
        int iso = 0;
        int charge = 0;
        int valence = 0;
        int nIsotopeH = IS_DEUTERIUM( i ) ? 1 : IS_TRITIUM( i ) ? 2 : 0;
        int bonds_val;
        bAtomNeedsAlias = ALIASED_AT( i );
        memset( elname, 0, sizeof( elname ) ); /* djb-rwth: memset_s C11/Annex K variant? */

        if (bAtomNeedsAlias)
        {
            /* alias */
            strcpy(elname, "C");
        }
        else
        {
            /* isotope*/
            if (nIsotopeH)
            {
                strcpy(elname, bAtomsDT ? (nIsotopeH == 1 ? "D" : "T") : "H");
            }
            else
            {
                strncpy(elname, at[i].elname, sizeof(elname) - 1);
                elname[sizeof(elname) - 1] = 0; /* adding zero termination after strncpy */
            }
            if (!ABNORMAL_CHG(i) && !ANY_RAD(i))
            {
                /* charge*/
                /* Only atoms without alias can be here*/
                switch (at[i].charge)
                {
                case 3:
                    charge = 1;
                    break;
                case 2:
                    charge = 2;
                    break;
                case 1:
                    charge = 3;
                    break;
                case -1:
                    charge = 5;
                    break;
                case -2:
                    charge = 6;
                    break;
                case -3:
                    charge = 7;
                    break;
                case 0:
                    charge = 0;
                    break;
                default:
                    break; /* djb-rwth: removing redundant code */
                }
            }

            /* radical*/
            if (ANY_RAD(i) && !ANY_CHG(i))
            {
                if (at[i].radical == RADICAL_DOUBLET)
                {
                    charge = 4;
                }
            }
        }

        /* allow isotopic shift for aliased atoms */
        if (NORMAL_ISO(i, bAtomsDT))
        {
            iso = at[i].iso_atw_diff > 0 ? at[i].iso_atw_diff - 1 :
                    at[i].iso_atw_diff < 0 ? at[i].iso_atw_diff :
                        nIsotopeH ? nIsotopeH : 0; /* djb-rwth: removing redundant code */
        }

        x = at[i].x;
        y = at[i].y;
        z = at[i].z;

        /* valence -- set only if needed */
        bonds_val = nBondsValenceInpAt(at + i, NULL, NULL);
        valence = needed_unusual_el_valence(at[i].el_number, at[i].charge, at[i].radical,
                                            at[i].chem_bonds_valence, bonds_val, NUMH(at, i), at[i].valence);

        if (valence < 0)
        {
            valence = 15; /* means no bonds nor H */
        }

        /* Convert "Zz" to "*" element symbol */

        if (!strcmp(elname, "Zz") || !strcmp(elname, "Zy"))
        {
            strcpy(elname, "*");
        }

        /*inchi_ios_eprint(fcb,"%10.4f%10.4f%10.4f %-3.3s%2d%3d  0     0  0  0  0  0  0\n",*/
        /*    (float)at[i].x, (float)(-at[i].y), fzero, at[i].elname, iso, charge);*/
        /*              xxxxxxyyyyyyzzzzzz aaa____ddcccsssnnnbbbvvvrrriiimmmeee  */
        inchi_ios_print_nodisplay(fcb, "%10.4f%10.4f%10.4f %-3.3s%2d%3d  0     0%3d  0  0  0  0\n",
                                  x, y, z, elname, (int)iso, (int)charge, valence /* at[i].special*/);

        /* Reflect image against x-axis;                                    */
        /* when transforming MOLfile back to STDATA in mol_to_stdata(...),  */
        /* make one more reflection to restore original orientation.        */
        /* Reason: in MS Search y-axis is directed from top to bottom,      */
        /*         while in MOLfile y-axis goes from bottom to top.         */
    }

    return ret;
}

/****************************************************************************
 OrigAtData : Write To SDfile : Bonds Block
****************************************************************************/
int OrigAtData_WriteToSDfileBondsBlock( const ORIG_ATOM_DATA *inp_at_data,
                                        INCHI_IOSTREAM       *fcb,
                                        const char           *name,
                                        const char           *comment,
                                        const char           *szLabel,
                                        const char           *szValue,
                                        INT_ARRAY            *written_bond_ends )
{
    int i, j, k, ret = 0;
    int num_atoms = inp_at_data->num_inp_atoms;
    const inp_ATOM *at = inp_at_data->at;

    /* bonds*/
    for (i = 0; i < num_atoms; i++)
    {
        for (j = 0; j < at[i].valence; j++)
        {
            if (i < at[i].neighbor[j])
            {
                unsigned a1, a2;
                if ((k = at[i].bond_stereo[j])) /* djb-rwth: addressing LLVM warning */
                {
                    /* bond stereo */
                    if (k < 0)
                    {
                        /* transposition */
                        a1 = (unsigned)(at[i].neighbor[j] + 1);
                        a2 = (unsigned)(i + 1);
                        inchi_ios_print_nodisplay(fcb, "%3u%3u%3u%3u  0  0  0\n",
                                                  a1, a2, (unsigned)(at[i].bond_type[j]), (unsigned)abs(k));
                    }
                    else
                    {
                        /* no transposition*/
                        a1 = (unsigned)(i + 1);
                        a2 = (unsigned)(at[i].neighbor[j] + 1);
                        inchi_ios_print_nodisplay(fcb, "%3u%3u%3u%3u  0  0  0\n",
                                                  a1, a2, (unsigned)(at[i].bond_type[j]), (unsigned)abs(k));
                    }
                }
                else
                {
                    a1 = (unsigned)(i + 1);
                    a2 = (unsigned)(at[i].neighbor[j] + 1);
                    inchi_ios_print_nodisplay(fcb, "%3u%3u%3u  0  0  0  0\n",
                                              a1, a2, (unsigned)(at[i].bond_type[j]));
                }

                IntArray_Append(written_bond_ends, a1);
                IntArray_Append(written_bond_ends, a2);
            }
        }
    }

    return ret;
}

/****************************************************************************
 OrigAtData : Write To SDfile : Additional Lines
****************************************************************************/
int OrigAtData_WriteToSDfileAdditionalLines( const ORIG_ATOM_DATA *inp_at_data,
                                             INCHI_IOSTREAM       *fcb,
                                             const char           *name,
                                             const char           *comment,
                                             int                  bAtomsDT,
                                             const char           *szLabel,
                                             const char           *szValue,
                                             int                  nNumAliasLines,
                                             int                  nNumChargeLines,
                                             int                  nNumRadicalLines,
                                             int                  nNumIsoLines,
                                             INT_ARRAY            *written_bond_ends )
{
    char str_m[66], entry[25];
    int i, num_m, k, j, ret = 0;

    if (inp_at_data) /* djb-rwth: fixing a NULL pointer dereference */
    {
        int num_atoms = inp_at_data->num_inp_atoms;
        int is_polymer = inp_at_data && inp_at_data->polymer && inp_at_data->polymer->n > 0 && inp_at_data->valid_polymer;
        const inp_ATOM *at = inp_at_data->at;

        /* Aliases. 5-3-99 DCh.*/
        if (nNumAliasLines)
        {
            num_m = 0;
            for (i = 0; i < num_atoms; i++)
            {
                if (ALIASED_AT(i))
                {
                    int len;
                    inchi_ios_print_nodisplay(fcb, "A  %d\n", i + 1);
                    num_m++;
                    len = sprintf(str_m, "%s", at[i].elname);
                    /* add isotopic H to the alias */
                    for (k = 0; k < NUM_H_ISOTOPES; k++)
                    {
                        int num_H = at[i].num_iso_H[k] + (k ? 0 : at[i].num_H);
                        if (num_H)
                        {
                            len += sprintf(str_m + len, "%s", k == 0 ? "H" : k == 1 ? "D" : k == 2 ? "T" : "?");
                            if (num_H != 1)
                            {
                                len += sprintf(str_m + len, "%d", num_H);
                            }
                        }
                    }

                    /* Add charge to the Alias */
                    if (at[i].charge)
                    {
                        len += sprintf(str_m + len, "%s", at[i].charge > 0 ? "+" : "-");
                        if (1 < (j = abs(at[i].charge)))
                        {
                            len += sprintf(str_m + len, "%d", j);
                        }
                    }

                    /* Add radical to the Alias */
                    if (at[i].radical == RADICAL_SINGLET)
                    {
                        len += sprintf(str_m + len, "%s", ":");
                    }
                    else if (at[i].radical == RADICAL_DOUBLET)
                    {
                        len += sprintf(str_m + len, "%s", "^");
                    }
                    else if (at[i].radical == RADICAL_TRIPLET)
                    {
                        len += sprintf(str_m + len, "%s", "^^");
                    }
                    inchi_ios_print_nodisplay(fcb, "%s\n", str_m);
                    num_m++;
                }
            }

            if (num_m != nNumAliasLines)
            {
                /* error in lines counting*/
                ret++;
            }
        }

        /* charges*/
        str_m[0] = 0;
        num_m = 0;
        if (nNumChargeLines)
        {
            for (i = 0; i < num_atoms; i++)
            {
                if (at[i].charge && !ALIASED_AT(i))
                {
                    sprintf(entry, " %3d %3d", i + 1, (int)at[i].charge);
                    strcat(str_m, entry);
                    num_m++;
                }
                if ((i == num_atoms - 1 && num_m) || num_m == 8) /* djb-rwth: addressing LLVM warning */
                {
                    inchi_ios_print_nodisplay(fcb, "M  CHG%3d%s\n", num_m, str_m);
                    str_m[0] = 0;
                    num_m = 0;
                }
            }
        }

        /* radicals*/
        str_m[0] = 0;
        num_m = 0;

        if (nNumRadicalLines)
        {
            for (i = 0; i < num_atoms; i++)
            {
                if (at[i].radical && !ALIASED_AT(i))
                {
                    int radical = (at[i].radical == RADICAL_SINGLET ||
                                   at[i].radical == RADICAL_DOUBLET ||
                                   at[i].radical == RADICAL_TRIPLET) ? at[i].radical : 0;
                    if (radical)
                    {
                        sprintf(entry, " %3d %3d", i + 1, radical);
                        strcat(str_m, entry);
                        num_m++;
                    }
                }
                if ((i == num_atoms - 1 && num_m) || num_m == 8) /* djb-rwth: addressing LLVM warning */
                {
                    inchi_ios_print_nodisplay(fcb, "M  RAD%3d%s\n", num_m, str_m);
                    str_m[0] = 0;
                    num_m = 0;
                }
            }
        }

        /* isotopes*/
        str_m[0] = 0;
        num_m = 0;
        if (nNumIsoLines)
        {
            int el_num, iso;
            for (i = 0; i < num_atoms; i++)
            {
                /*
                if ( 0 == strcmp( at[i].elname, "D" ) ) {
                    sprintf( entry, " %3d %3d", i+1, 2 );
                    strcat( str_m, entry );
                    num_m ++;
                } else
                if ( 0 == strcmp( at[i].elname, "T" ) ) {
                    sprintf( entry, " %3d %3d", i+1, 3 );
                    strcat( str_m, entry );
                    num_m ++;
                } else
                if ( k = at[i].iso_atw_diff ) {
                    int mw = get_atomic_mass_from_elnum( at[i].el_number );
                    mw += (k > 0)? k-1 : k;
                    sprintf( entry, " %3d %3d", i+1, mw );
                    strcat( str_m, entry );
                    num_m ++;
                }
                */

                if (ANY_ISO(i, bAtomsDT) && !ALIASED_AT(i))
                {
                    if (IS_DEUTERIUM(i))
                    {
                        iso = 1;
                        el_num = 1;
                    }
                    else if (IS_TRITIUM(i))
                    {
                        iso = 2;
                        el_num = 1;
                    }
                    else
                    {
                        iso = at[i].iso_atw_diff > 0 ? at[i].iso_atw_diff - 1 : at[i].iso_atw_diff;
                        el_num = at[i].el_number;
                    }
                    iso += get_atomic_mass_from_elnum(el_num);
                    sprintf(entry, " %3d %3d", i + 1, iso);
                    strcat(str_m, entry);
                    num_m++;
                }

                if ((i == num_atoms - 1 && num_m) || num_m == 8) /* djb-rwth: addressing LLVM warning */
                {
                    inchi_ios_print_nodisplay(fcb, "M  ISO%3d%s\n", num_m, str_m);
                    str_m[0] = 0;
                    num_m = 0;
                }
            }
        }

        if (is_polymer)
        {
            OrigAtData_WriteToSDfilePolymerData(inp_at_data, fcb, name, comment,
                                                szLabel, szValue, written_bond_ends);
        }

        inchi_ios_print_nodisplay(fcb, "M  END\n");
    }
    return ret;
}

/****************************************************************************
OrigAtData : Write To SDfile : Polymer Data
****************************************************************************/
int OrigAtData_WriteToSDfilePolymerData( const ORIG_ATOM_DATA *inp_at_data,
                                         INCHI_IOSTREAM       *fcb,
                                         const char           * name,
                                         const char           *comment,
                                         const char           *szLabel,
                                         const char           *szValue,
                                         INT_ARRAY            *written_bond_ends )
{
    int j, k, ju, jj, jprev, ret = 0;
    const char *sty[] = {"NON", "SRU", "MON", "COP", "MOD", "CRO", "MER"};
    const char *sst[] = {"NON", "ALT", "RAN", "BLO"};
    const char *con[] = {"NON", "HT", "HH", "EU"};
    OAD_PolymerUnit *u = NULL;

    /* STY */
    jj = 0;
    jprev = -1;
    for (j = 0; j < inp_at_data->polymer->n; j++)
    {
        u = inp_at_data->polymer->units[j];
        if (u->type > 0 && u->type <= 6)
        {
            jj++;
        }
        if (jj == 8 || j == inp_at_data->polymer->n - 1)
        {
            inchi_ios_print_nodisplay(fcb, "M  STY%3d", jj % 8 ? jj % 8 : 8);
            for (k = jprev + 1; k <= j; k++)
            {
                u = inp_at_data->polymer->units[k];
                if (u->type > 0 && u->type <= 6)
                {
                    inchi_ios_print_nodisplay(fcb, " %3d %3s", u->id, sty[u->type]);
                }
            }
            inchi_ios_print_nodisplay(fcb, "\n");
            jj = 0;
            jprev = j;
        }
    }
    /* SLB */
    /* djb-rwth: removing redundant code */
    jprev = -1;
    for (j = 0; j < inp_at_data->polymer->n; j++)
    {
        /* djb-rwth: removing redundant code */
        if (j == 8 || j == inp_at_data->polymer->n - 1)
        {
            jj = j + 1;
            inchi_ios_print_nodisplay(fcb, "M  SLB%3d", jj % 8 ? jj % 8 : 8);
            for (k = jprev + 1; k < jj; k++)
            {
                u = inp_at_data->polymer->units[k];
                inchi_ios_print_nodisplay(fcb, " %3d %3d", u->id, u->label);
            }
            inchi_ios_print_nodisplay(fcb, "\n");
            /* djb-rwth: removing redundant code */
            jprev = j;
        }
    }

    /* SST */
    jj = 0;
    /* djb-rwth: removing redundant code */
    for (j = 0; j < inp_at_data->polymer->n; j++)
    {
        u = inp_at_data->polymer->units[j];
        if (u->subtype == MOL_FMT_M_SST_ALT || u->subtype == MOL_FMT_M_SST_RAN || u->subtype == MOL_FMT_M_SST_BLK)
        {
            jj++;
        }
    }
    if (jj)
    {
        jj = 0;
        jprev = -1;
        for (j = 0; j < inp_at_data->polymer->n; j++)
        {
            u = inp_at_data->polymer->units[j];
            if (u->subtype == MOL_FMT_M_SST_ALT || u->subtype == MOL_FMT_M_SST_RAN || u->subtype == MOL_FMT_M_SST_BLK)
            {
                jj++;
            }
            if (jj == 8 || j == inp_at_data->polymer->n - 1)
            {
                inchi_ios_print_nodisplay(fcb, "M  SST%3d", jj % 8 ? jj % 8 : 8);
                for (k = jprev + 1; k <= j; k++)
                {
                    u = inp_at_data->polymer->units[k];
                    if (u->subtype == MOL_FMT_M_SST_ALT || u->subtype == MOL_FMT_M_SST_RAN || u->subtype == MOL_FMT_M_SST_BLK)
                    {
                        inchi_ios_print_nodisplay(fcb, " %3d %3s", u->id, sst[u->subtype]);
                    }
                }
                inchi_ios_print_nodisplay(fcb, "\n");
                jj = 0;
                jprev = j;
            }
        }
    }

    /* SCN */
    jj = 0;
    /* djb-rwth: removing redundant code */
    for (j = 0; j < inp_at_data->polymer->n; j++)
    {
        u = inp_at_data->polymer->units[j];
        if (u->conn == MOL_FMT_M_CONN_HT || u->conn == MOL_FMT_M_CONN_HH || u->conn == MOL_FMT_M_CONN_EU)
        {
            jj++;
        }
    }
    if (jj)
    {
        jj = 0;
        jprev = -1;
        for (j = 0; j < inp_at_data->polymer->n; j++)
        {
            u = inp_at_data->polymer->units[j];
            if (u->conn == MOL_FMT_M_CONN_HT || u->conn == MOL_FMT_M_CONN_HH || u->conn == MOL_FMT_M_CONN_EU)
            {
                jj++;
            }
            if (jj == 8 || j == inp_at_data->polymer->n - 1)
            {
                inchi_ios_print_nodisplay(fcb, "M  SCN%3d", jj % 8 ? jj % 8 : 8);
                for (k = jprev + 1; k <= j; k++)
                {
                    u = inp_at_data->polymer->units[k];
                    if (u->conn == MOL_FMT_M_CONN_HT || u->conn == MOL_FMT_M_CONN_HH || u->conn == MOL_FMT_M_CONN_EU)
                    {
                        inchi_ios_print_nodisplay(fcb, " %3d %3s", u->id, con[u->conn]);
                    }
                }
                inchi_ios_print_nodisplay(fcb, "\n");
                jj = 0;
                jprev = j;
            }
        }
    }
    /* SAL */
    for (ju = 0; ju < inp_at_data->polymer->n; ju++)
    {
        u = inp_at_data->polymer->units[ju];
        jj = 0;
        jprev = -1;
        for (j = 0; j < u->na; j++)
        {
            jj++;
            if (jj == 15 || j == u->na - 1)
            {
                inchi_ios_print_nodisplay(fcb, "M  SAL %3d%3d", u->id, jj % 15 ? jj % 15 : 15);
                for (k = jprev + 1; k <= j; k++)
                {
                    inchi_ios_print_nodisplay(fcb, " %3d", u->alist[k]);
                }
                inchi_ios_print_nodisplay(fcb, "\n");
                jj = 0;
                jprev = j;
            }
        }
    }
    /* SBL */
    for (ju = 0; ju < inp_at_data->polymer->n; ju++)
    {
        u = inp_at_data->polymer->units[ju];
        jj = 0;
        jprev = -1;
        for (j = 0; j < u->nb; j++)
        {
            jj++;
            if (jj == 15 || j == u->nb - 1)
            {
                inchi_ios_print_nodisplay(fcb, "M  SBL %3d%3d", u->id, jj % 15 ? jj % 15 : 15);

                for (k = jprev + 1; k <= j; k++)
                {
                    int a1, a2, e1, e2, wb, bond_num = 0;
                    a1 = u->blist[2 * k];
                    a2 = u->blist[2 * k + 1];
                    for (wb = 0; wb < written_bond_ends->used / 2; wb++)
                    {
                        e1 = written_bond_ends->item[2 * wb];
                        e2 = written_bond_ends->item[2 * wb + 1];
                        if ((a1 == e1 && a2 == e2) || (a2 == e1 && a1 == e2))
                        {
                            bond_num = wb + 1;
                            break;
                        }
                    }
                    if (bond_num)
                    {
                        inchi_ios_print_nodisplay(fcb, " %3d", bond_num);
                    }
                }

                inchi_ios_print_nodisplay(fcb, "\n");
                jj = 0;
                jprev = j;
            }
        }
    }

    /* SDI */
    for (j = 0; j < inp_at_data->polymer->n; j++)
    {
        /* better than nothing */
        float xmin, xmax, ymin, ymax;
        xmin = ymin = -1.0 * ((float)j + 1.0); /* djb-rwth: cast operator added, 1->1.0 */
        xmax = ymax = +1.0 * ((float)j + 1.0); /* djb-rwth: cast operator added, 1->1.0 */
        u = inp_at_data->polymer->units[j];
        /* u->xbr1[0], x1, y1, x2, y2 u->xbr1[1], u->xbr1[2], u->xbr1[3] */
        inchi_ios_print_nodisplay(fcb, "M  SDI %3d%3d%10.4f%10.4f%10.4f%10.4f\n", u->id, 4, xmin, ymin, xmin, ymax);
        /* u->xbr2[0], u->xbr2[1], u->xbr2[2], u->xbr2[3] */
        inchi_ios_print_nodisplay(fcb, "M  SDI %3d%3d%10.4f%10.4f%10.4f%10.4f\n", u->id, 4, xmax, ymax, xmax, ymin);
    }

    return ret;
}

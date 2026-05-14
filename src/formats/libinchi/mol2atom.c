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
#include <string.h>
#include <math.h>

#include "mode.h"
#include "mol_fmt.h"
#include "ichierr.h"
#include "util.h"
#include "strutil.h"
#include "inchi_api.h"

#include "bcf_s.h"

#include "logging.h" /*(@nnuk : Nauman Ullah Khan) :: Needed for logging functionality*/

/*
Convert input molecular data to internal representation

*/

/* Local prototypes */

int ReadMolfileToInpAtoms(INCHI_IOSTREAM *inp_file,
                          int bDoNotAddH, inp_ATOM **at,
                          MOL_COORD **szCoord,
                          OAD_Polymer **polymer,
                          OAD_V3000 **v3000,
                          int treat_polymers,
                          int treat_NPZz,
                          int max_num_at,
                          int *num_dimensions,
                          int *num_bonds,
                          const char *pSdfLabel,
                          char *pSdfValue,
                          unsigned long *Id,
                          long *lMolfileNumber,
                          INCHI_MODE *pInpAtomFlags,
                          int *err,
                          char *pStrErr,
                          int bNoWarnings);
inp_ATOM *MakeInpAtomsFromMolfileData(MOL_FMT_DATA *mfdata,
                                      int *num_atoms,
                                      int *num_bonds,
                                      inp_ATOM *at_inp,
                                      int bDoNotAddH,
                                      int *err,
                                      char *pStrErr);
int SetInpAtomsXYZ(MOL_FMT_DATA *mfdata,
                   int num_atoms,
                   inp_ATOM *at,
                   int *err,
                   char *pStrErr);
void calculate_valences(MOL_FMT_DATA *mfdata,
                        inp_ATOM *at,
                        int *num_atoms,
                        int bDoNotAddH,
                        int *err,
                        char *pStrErr);

int SetExtOrigAtDataByMolfileExtInput(MOL_FMT_DATA *mfdata,
                                      OAD_Polymer **polymer,
                                      OAD_V3000 **v3000,
                                      char *pStrErr);

/****************************************************************************
Create OrigInpData From Molfile
****************************************************************************/
int CreateOrigInpDataFromMolfile(INCHI_IOSTREAM *inp_file,
                                 ORIG_ATOM_DATA *orig_at_data,
                                 int bMergeAllInputStructures,
                                 int bGetOrigCoord,
                                 int bDoNotAddH,
                                 int treat_polymers,
                                 int treat_NPZz,
                                 const char *pSdfLabel,
                                 char *pSdfValue,
                                 unsigned long *lSdfId,
                                 long *lMolfileNumber,
                                 INCHI_MODE *pInpAtomFlags,
                                 int *err,
                                 char *pStrErr,
                                 int bNoWarnings)
{
    int i, j, max_num_at;
    int num_dimensions_new = 0; /* djb-rwth: initialisation required to avoid garbage values */
    int num_inp_bonds_new;
    int num_inp_atoms_new;
    int nNumAtoms = 0;
    inp_ATOM *at_new = NULL;
    inp_ATOM *at_old = NULL;
    MOL_COORD *szCoordNew = NULL;
    MOL_COORD *szCoordOld = NULL;
    OAD_Polymer *polymer = NULL;
    OAD_V3000 *v3000 = NULL;

    if (pStrErr)
    {
        pStrErr[0] = '\0';
    }

    /* NB: currently (v. 1.04) legacy CLI option "MERGE" is unsupported
    so the loop below is always a single pass */
    max_num_at = MAX_ATOMS;
    if (!(*pInpAtomFlags & FLAG_SET_INP_LARGE_MOLS))
    {
        max_num_at = NORMALLY_ALLOWED_INP_MAX_ATOMS;
    }

    do /* while ( !*err && bMergeAllInputStructures ) */
    {

        at_old = orig_at_data ? orig_at_data->at : NULL; /*  save pointer to the previous allocation */

        szCoordOld = orig_at_data ? orig_at_data->szCoord : NULL;

        num_inp_atoms_new = ReadMolfileToInpAtoms(inp_file, bDoNotAddH,
                                                  orig_at_data ? &at_new : NULL,
                                                  (bGetOrigCoord && orig_at_data) ? &szCoordNew : NULL,
                                                  &polymer, &v3000,
                                                  treat_polymers, treat_NPZz,
                                                  max_num_at, &num_dimensions_new, &num_inp_bonds_new,
                                                  pSdfLabel, pSdfValue, lSdfId, lMolfileNumber,
                                                  pInpAtomFlags, err, pStrErr, bNoWarnings);

        if (num_inp_atoms_new <= 0 && !*err)
        {
            TREAT_ERR(*err, 0, "Empty structure");
            *err = 98;
        }
        else if (orig_at_data && !num_inp_atoms_new &&
                 10 < *err && *err < 20 &&
                 orig_at_data->num_inp_atoms > 0 &&
                 bMergeAllInputStructures)
        {
            *err = 0; /* end of file */
            break;
        }
        else if (num_inp_atoms_new > 0 && orig_at_data)
        {
            /*  merge pOrigDataTmp + orig_at_data => pOrigDataTmp; */
            nNumAtoms = num_inp_atoms_new + orig_at_data->num_inp_atoms;
            if (nNumAtoms >= max_num_at) /*MAX_ATOMS )  */
            {
                TREAT_ERR(*err, 0, "Too many atoms [did you forget 'LargeMolecules' switch?]");
                *err = 70;
                orig_at_data->num_inp_atoms = -1;
            }
            else if (!at_old)
            {
                /* the first structure */
                orig_at_data->at = at_new;
                orig_at_data->szCoord = szCoordNew;
                at_new = NULL;
                szCoordNew = NULL;
                orig_at_data->num_inp_atoms = num_inp_atoms_new;
                orig_at_data->num_inp_bonds = num_inp_bonds_new;
                orig_at_data->num_dimensions = num_dimensions_new;
                /* v 1.05 */
                orig_at_data->polymer = polymer;
                orig_at_data->v3000 = v3000;
                polymer = NULL;
                v3000 = NULL;
            }
            else if ((orig_at_data->at = (inp_ATOM *)inchi_calloc(nNumAtoms, sizeof(inp_ATOM))) &&
                     (!szCoordNew || (orig_at_data->szCoord = (MOL_COORD *)inchi_calloc(nNumAtoms, sizeof(MOL_COORD)))))
            {
                /*  switch at_new <--> orig_at_data->at; */
                if (orig_at_data->num_inp_atoms)
                {
                    memcpy(orig_at_data->at, at_old, orig_at_data->num_inp_atoms * sizeof(orig_at_data->at[0]));
                    /*  adjust numbering in the newly read structure */
                    for (i = 0; i < num_inp_atoms_new; i++)
                    {
                        if (at_new) /* djb-rwth: fixing a NULL pointer dereference */
                        {
                            for (j = 0; j < at_new[i].valence; j++)
                            {
                                at_new[i].neighbor[j] += orig_at_data->num_inp_atoms;
                            }
                            at_new[i].orig_at_number += orig_at_data->num_inp_atoms; /* 12-19-2003 */
                        }
                    }
                    if (orig_at_data->szCoord && szCoordOld)
                        memcpy(orig_at_data->szCoord, szCoordOld, orig_at_data->num_inp_atoms * sizeof(MOL_COORD));
                }
                if (at_old)
                {
                    inchi_free(at_old);
                    at_old = NULL;
                }
                if (szCoordOld)
                {
                    inchi_free(szCoordOld);
                    szCoordOld = NULL;
                }

                /*  Copy newly read structure */
                if (at_new) /* djb-rwth: fixing a NULL pointer dereference */
                    memcpy(orig_at_data->at + orig_at_data->num_inp_atoms, at_new, num_inp_atoms_new * sizeof(orig_at_data->at[0]));
                if (orig_at_data->szCoord && szCoordNew)
                {
                    memcpy(orig_at_data->szCoord + orig_at_data->num_inp_atoms,
                           szCoordNew,
                           num_inp_atoms_new * sizeof(MOL_COORD));
                }

                /*  Add other things */
                orig_at_data->num_inp_atoms += num_inp_atoms_new;
                orig_at_data->num_inp_bonds += num_inp_bonds_new;
                orig_at_data->num_dimensions = inchi_max(num_dimensions_new, orig_at_data->num_dimensions);
                /* v 1.05 */
                orig_at_data->polymer = polymer;
                orig_at_data->v3000 = v3000;
            }
            else
            {
                TREAT_ERR(*err, 0, "Out of RAM");
                *err = -1;
            }
        }
        else if (num_inp_atoms_new > 0)
        {
            nNumAtoms += num_inp_atoms_new;
        }

        if (at_new)
        {
            inchi_free(at_new);
            at_new = NULL;
        }
        if (polymer)
        {
            inchi_free(polymer);
            polymer = NULL;
        }
        if (v3000)
        {
            inchi_free(v3000);
            v3000 = NULL;
        }

    } while (!*err && bMergeAllInputStructures);

    if (szCoordNew)
    {
        inchi_free(szCoordNew);
    }
    if (at_new)
    {
        inchi_free(at_new);
    }

    if (*err)
    {
        FreeOrigAtData(orig_at_data);
    }

    if (*err && !(10 < *err && *err < 20) && pStrErr && !pStrErr[0])
    {
        TREAT_ERR(*err, 0, "Unknown error"); /*   <BRKPT> */
    }

    return orig_at_data ? orig_at_data->num_inp_atoms : nNumAtoms;
}

/****************************************************************************
ReadMolfileToInpAtoms
****************************************************************************/
int ReadMolfileToInpAtoms(INCHI_IOSTREAM *inp_file,
                          int bDoNotAddH,
                          inp_ATOM **at,
                          MOL_COORD **szCoord,
                          OAD_Polymer **polymer,
                          OAD_V3000 **v3000,
                          int treat_polymers,
                          int treat_NPZz,
                          int max_num_at,
                          int *num_dimensions,
                          int *num_bonds,
                          const char *pSdfLabel,
                          char *pSdfValue,
                          unsigned long *Id,
                          long *lMolfileNumber,
                          INCHI_MODE *pInpAtomFlags,
                          int *err,
                          char *pStrErr,
                          int bNoWarnings)
{
    int num_atoms = 0;
    MOL_FMT_DATA *mfdata = NULL;
    MOL_FMT_HEADER_BLOCK OnlyHeaderBlock,
        *pOnlyHeaderBlock = NULL,
        *pHdr;
    MOL_FMT_CTAB OnlyCTab,
        *pOnlyCTab = NULL;
    char cSdfValueFirstChar = 0;

    if (at)
    {
        pOnlyHeaderBlock = NULL;
        if (*at && max_num_at)
        {
            memset(*at, 0, max_num_at * sizeof(inp_ATOM)); /* djb-rwth: memset_s C11/Annex K variant? */
        }
        if (szCoord && *szCoord)
        {
            inchi_free(*szCoord);
            *szCoord = NULL;
        }
    }
    else
    {
        pOnlyHeaderBlock = &OnlyHeaderBlock;
        pOnlyCTab = &OnlyCTab;
    }

    if (pSdfValue)
    {
        cSdfValueFirstChar = pSdfValue[0];
        pSdfValue[0] = '\0';
    }

    /* Read mol-formatted file (MOL or SD) block  */
    mfdata = ReadMolfile(inp_file, pOnlyHeaderBlock, pOnlyCTab, NULL != szCoord,
                         treat_polymers, treat_NPZz,
                         NULL, 0,
                         Id, pSdfLabel, pSdfValue, err, pStrErr, bNoWarnings);

    pHdr = (mfdata && !pOnlyHeaderBlock)
               ? &mfdata->hdr
           : (!mfdata && pOnlyHeaderBlock) ? pOnlyHeaderBlock
                                           : NULL;

    if (lMolfileNumber && pHdr)
    {
        *lMolfileNumber = MolfileExtractStrucNum(pHdr);
    }

    if (pSdfValue && !pSdfValue[0] &&
        pSdfLabel && pSdfLabel[0] &&
        pHdr)
    {
        if (!inchi_stricmp(pSdfLabel, "MolfileName"))
        {
            mystrncpy(pSdfValue, pHdr->molname, MAX_SDF_VALUE + 1);
            lrtrim(pSdfValue, NULL);
        }
        else if (!inchi_stricmp(pSdfLabel, "MolfileLine2"))
        {
            mystrncpy(pSdfValue, pHdr->line2, MAX_SDF_VALUE + 1);
            lrtrim(pSdfValue, NULL);
        }
        else if (!inchi_stricmp(pSdfLabel, "MolfileComment"))
        {
            mystrncpy(pSdfValue, pHdr->comment, MAX_SDF_VALUE + 1);
            lrtrim(pSdfValue, NULL);
        }
        else if (!inchi_stricmp(pSdfLabel, "MolfileIntRegNo") && pHdr->internal_regno)
        {
            sprintf(pSdfValue, "%ld", pHdr->internal_regno);
        }

        if (!pSdfValue[0])
        {
            pSdfValue[0] = cSdfValueFirstChar;
        }
    }

    if (mfdata && at && !*err)
    {
        /* (*at) either points to already allocated memory or NULL */
        if (mfdata->ctab.n_atoms <= max_num_at)
        {

            *at = MakeInpAtomsFromMolfileData(mfdata, &num_atoms, num_bonds,
                                              *at, bDoNotAddH, err, pStrErr);

            if (*err >= 0)
            {
                *num_dimensions = SetInpAtomsXYZ(mfdata, num_atoms,
                                                 *at, err, pStrErr);

                if (szCoord)
                {
                    *szCoord = mfdata->ctab.coords;
                    mfdata->ctab.coords = NULL;
                }

                /* Check that all atoms nos of Hs are in allowed range ( <=MAXVAL, currently 20 ) */
                if (1)
                {
                    int i;
                    for (i = 0; i < num_atoms; i++)
                    {
                        if ((*at)[i].num_iso_H[0] + (*at)[i].num_iso_H[1] + (*at)[i].num_iso_H[2] > MAXVAL)
                        {
                            *err = 70 + 8;
                            TREAT_ERR(*err, 0, "Too many hydrogens at heavy atom");
                            num_atoms = -1;
                            break;
                        }
                    }
                }

                if (!*err)
                {
                    *err = SetExtOrigAtDataByMolfileExtInput(mfdata, polymer, v3000, pStrErr);
                    if (*err)
                    {
                        /*TREAT_ERR (*err, 0, "Error while getting extended Molfile input");*/
                        *err = 80;
                        num_atoms = -1;
                    }
                }
            }
            else
            {
                if (*err == -2) /* unk element */
                {
                    *err = 90;
                    num_atoms = -1;
                }
            }
        }
        else
        {
            /* non-affordable struct size */
            TREAT_ERR(*err, 0, "Too many atoms [did you forget 'LargeMolecules' switch?]");
            *err = 70;
            num_atoms = -1;
        }

        if (*err > 0)
        {
            *err += 100;
        }

        /* 11-16-2004: use Chiral flag */
        if (num_atoms > 0 && at && *at && mfdata && pInpAtomFlags)
        {
            if (mfdata->ctab.chiral_flag)
            {
                *pInpAtomFlags |= FLAG_INP_AT_CHIRAL;
            }
            else
            {
                *pInpAtomFlags |= FLAG_INP_AT_NONCHIRAL;
            }
        }
    }

    else if (!at)
    {
        num_atoms = pOnlyCTab->n_atoms;
    }

    if (!pOnlyHeaderBlock)
    {
        FreeMolfileData(mfdata);
    }

    return num_atoms;
}

/****************************************************************************
Make Inp Atoms From Molfile Data
****************************************************************************/
inp_ATOM *MakeInpAtomsFromMolfileData(MOL_FMT_DATA *mfdata,
                                      int *num_atoms,
                                      int *num_bonds,
                                      inp_ATOM *at_inp,
                                      int bDoNotAddH,
                                      int *err,
                                      char *pStrErr)
{
    inp_ATOM *at = NULL;
    /* char      *bond_stereo = NULL; */
    AT_NUMB *p1, *p2;
    int i, a1, a2, n1, n2, bonds, iso_atw_diff;
    char bond_stereo, bond_type;

    *err = 0;

    *num_atoms = mfdata->ctab.n_atoms;
    *num_bonds = 0;

    /*(@nnuk : Nauman Ullah Khan) */
    LOG_NO_ARGS("\n############### (L579:mol2atom.c) ################\n");
    LOG_MULT_ARGS("Number of atoms : %d\n", *num_atoms);
    LOG_NO_ARGS("####################################################\n");

    if (MolfileHasNoChemStruc(mfdata))
    {
        goto exit_function;
    }

    /* Allocate memory if necessary */
    if (at_inp)
    {
        at = at_inp;
    }
    else
    {
        at = CreateInpAtom(*num_atoms);
        if (!at)
        {
            *err = -1;
            TREAT_ERR_AND_FIN(*err, -1, exit_function, "Out of RAM");
        }
    }

    /* Copy atoms info */

    /*(@nnuk : Nauman Ullah Khan) */
    LOG_NO_ARGS("\n##################### Atoms Data #########################\n");

    for (i = 0; i < *num_atoms; i++)
    {

        mystrncpy(at[i].elname, mfdata->ctab.atoms[i].symbol, sizeof(at->elname));
        /* at[i].chem_bonds_valence = mfdata->ctab.atoms[i].valence; */
        /*  MOLfile valence; will change */

        at[i].orig_at_number = (AT_NUMB)(i + 1);
        at[i].iso_atw_diff = mfdata->ctab.atoms[i].mass_difference;
        at[i].charge = mfdata->ctab.atoms[i].charge;
        at[i].radical = mfdata->ctab.atoms[i].radical;

        /* see SetInpAtomXYZ()
        at[i].x                  = mfdata->ctab.atoms[i].fx;
        at[i].y                  = mfdata->ctab.atoms[i].fy;
        at[i].z                  = mfdata->ctab.atoms[i].fz;
        */

        iso_atw_diff = mfdata->ctab.atoms[i].mass_difference;

        at[i].iso_atw_diff = iso_atw_diff == ZERO_ATW_DIFF
                                 ? 1
                             : iso_atw_diff > 0 ? iso_atw_diff + 1
                                                : iso_atw_diff;

#if (SINGLET_IS_TRIPLET == 1)
        if (at[i].radical == RADICAL_SINGLET)
        {
            at[i].radical = RADICAL_TRIPLET;
        }
#endif
#if (bRELEASE_VERSION != 1)
        if (isdigit(at[i].elname[0]))
        { /*  for testing */
            mystrncpy(at[i].elname, "C", sizeof(at->elname));
        }
#endif

        if (ERR_ELEM == (n1 = get_periodic_table_number(at[i].elname)))
        {
            /*  Case when elname contains more than 1 element: extract number of H if possible */
            at[i].num_H = extract_H_atoms(at[i].elname, at[i].num_iso_H);

            if (!at[i].elname[0] && NUMH(at, i))
            {
                /* alias contains only H. Added 2004-07-21, fixed 2004-07-22
                 * move the heaviest isotope to the "central atom"
                 * Note: this must be consistent with H-H treatment in remove_terminal_HDT()
                 */
                strcpy(at[i].elname, "H");
                if (NUM_ISO_H(at, i))
                {
                    int j;
                    for (j = NUM_H_ISOTOPES - 1; 0 <= j; j--)
                    {
                        if (at[i].num_iso_H[j])
                        {
                            at[i].num_iso_H[j]--;
                            at[i].iso_atw_diff = 1 + j;
                            break;
                        }
                    }
                }
                else
                {
                    at[i].num_H--;
                }
            }
            if (ERR_ELEM == (n1 = get_periodic_table_number(at[i].elname)))
            {
                n1 = 0;
            }
        } /* if ( ERR_ELEM == */

        at[i].el_number = (U_CHAR)n1;
        if (!n1)
        {
            *err = -2;
            TREAT_ERR(*err, -2, "Unknown element(s):");
            TREAT_ERR_AND_FIN(*err, -2, exit_function, at[i].elname);
        }
        else
        {
            /* replace explicit D or T with isotopic H (added 2003-06-02) */
            if ((int)EL_NUMBER_H == n1 && !at[i].iso_atw_diff)
            {
                switch (at[i].elname[0])
                {
                case 'D':
                    at[i].iso_atw_diff = 2;
                    mystrncpy(at[i].elname, "H", sizeof(at->elname));
                    break;
                case 'T':
                    at[i].iso_atw_diff = 3;
                    mystrncpy(at[i].elname, "H", sizeof(at->elname));
                    break;
                }
            }
        }

        /*(@nnuk : Nauman Ullah Khan) */
        LOG_NO_ARGS("\n############################## (L714:mol2atom.c) #####################################\n");
        LOG_MULT_ARGS("Atom %d: element=%s, x=%f, y=%f, z=%f, chrg=%d, rad=%d, iso=%d\n", i, at[i].elname, at[i].x, at[i].y, at[i].z, at[i].charge, at[i].radical, at[i].iso_atw_diff);
        LOG_NO_ARGS("########################################################################################\n");

    } /* eof copy atom info */

    /*---------------- stereo information notes. ------------------------

    Currently:  1. stereo sign
    =========   --------------
    MOLfile     (atom number = MOLfile atom number - 1, no stdata as an intermediate)
    |        if mfdata->ctab.bonds[i].atnum1 < mfdata->ctab.bonds[i].atnum2
    v        then
    inp_ATOM        stereo > 0
    else
    stereo < 0

    2. neighbor z-coordinate
    ------------------------
    neighbor z-coord > 0 for Up if sign(stdata_bond_no) = sign(at[i].neighbor[j]-i)

    --------------------------------------------------------------------*/

    /* Copy bond info */

    LOG_NO_ARGS("\n######################### Bonds Data ###############################\n");

    for (i = 0, bonds = 0; i < mfdata->ctab.n_bonds; i++)
    {

        bond_stereo = mfdata->ctab.bonds[i].bond_stereo;
        bond_type = mfdata->ctab.bonds[i].bond_type;

        a1 = mfdata->ctab.bonds[i].atnum1 - 1;
        a2 = mfdata->ctab.bonds[i].atnum2 - 1;

        if (a1 < 0 || a1 >= *num_atoms ||
            a2 < 0 || a2 >= *num_atoms ||
            a1 == a2)
        {
            *err |= 1; /*  bond for impossible atom number(s); ignored */
            TREAT_ERR(*err, 0, "Bond to nonexistent atom");
            continue;
        }

        /*  check for multiple bonds between same atoms */
        p1 = is_in_the_list(at[a1].neighbor, (AT_NUMB)a2, at[a1].valence);
        p2 = is_in_the_list(at[a2].neighbor, (AT_NUMB)a1, at[a2].valence);

        /*(@nnuk : Nauman Ullah Khan) */
        LOG_NO_ARGS("\n################## (L771:mol2atom.c) ##################\n");
        LOG_MULT_ARGS("Valence = %d\n", at[i].valence);
        LOG_NO_ARGS("#########################################################\n");

        if ((p1 || p2) && (p1 || at[a1].valence < MAXVAL) && (p2 || at[a2].valence < MAXVAL))
        {
            n1 = p1 ? (p1 - at[a1].neighbor) : at[a1].valence++;
            n2 = p2 ? (p2 - at[a2].neighbor) : at[a2].valence++;
            TREAT_ERR(*err, 0, "Multiple bonds between two atoms");
            *err |= 2; /*  multiple bonds between atoms */
        }
        else if (!p1 && !p2 && at[a1].valence < MAXVAL && at[a2].valence < MAXVAL)
        {
            n1 = at[a1].valence++;
            n2 = at[a2].valence++;
            bonds++;
        }
        else
        {
            char szMsg[64];
            *err |= 4; /*  too large number of bonds. Some bonds ignored. */
            sprintf(szMsg, "Atom '%s' has more than %d bonds",
                    at[a1].valence >= MAXVAL ? at[a1].elname : at[a2].elname, MAXVAL);
            TREAT_ERR(*err, 0, szMsg);
            continue;
        }

        if (bond_type < MIN_INPUT_BOND_TYPE || bond_type > MAX_INPUT_BOND_TYPE)
        {
            char szBondType[16];
            sprintf(szBondType, "%d", bond_type);
            bond_type = 1;
            TREAT_ERR(*err, 0, "Unrecognized bond type:");
            TREAT_ERR(*err, 0, szBondType);
            *err |= 8; /*  Unrecognized Bond type replaced with single bond */
        }

        /* bond type */
        if (n1 < MAXVAL && n2 < MAXVAL) /* djb-rwth: buffer overrun prevention */
        {
            at[a1].bond_type[n1] = at[a2].bond_type[n2] = bond_type;

            /* connection */
            at[a1].neighbor[n1] = (AT_NUMB)a2;
            at[a2].neighbor[n2] = (AT_NUMB)a1;

            /* stereo */
            if (bond_stereo == INPUT_STEREO_DBLE_EITHER /* 3 */)
            {
                at[a1].bond_stereo[n1] =
                    at[a2].bond_stereo[n2] =
                        STEREO_DBLE_EITHER;
            }
            else if (bond_stereo == INPUT_STEREO_SNGL_UP ||     /* 1 */
                     bond_stereo == INPUT_STEREO_SNGL_EITHER || /* 4 */
                     bond_stereo == INPUT_STEREO_SNGL_DOWN /* 6 */)
            {
                char cStereo;
                switch (bond_stereo)
                {
                case INPUT_STEREO_SNGL_UP:
                    cStereo = STEREO_SNGL_UP;
                    break;
                case INPUT_STEREO_SNGL_EITHER:
                    cStereo = STEREO_SNGL_EITHER;
                    break;
                case INPUT_STEREO_SNGL_DOWN:
                    cStereo = STEREO_SNGL_DOWN;
                    break;
                }
                at[a1].bond_stereo[n1] = cStereo;  /*  >0: the wedge (pointed) end is at this atom, a1 */
                at[a2].bond_stereo[n2] = -cStereo; /*  <0: the wedge (pointed) end is at the opposite atom, a1 */
            }
            else if (bond_stereo)
            {
                *err |= 16; /*  Ignored unrecognized Bond stereo */
                TREAT_ERR(*err, 0, "Unrecognized bond stereo");
                continue;
            }
        }

        /*(@nnuk : Nauman Ullah Khan) */
        LOG_NO_ARGS("\n################ (L862:mol2atom.c) ##################\n");
        LOG_MULT_ARGS("Bond %d: atom1=%d, atom2=%d, type=%d, stereo=%d\n", i, a1, a2, bond_type, bond_stereo);
        LOG_NO_ARGS("#######################################################\n");

    } /* eof copy bond info */

    *num_bonds = bonds;

    /* special valences */
    calculate_valences(mfdata, at, num_atoms, bDoNotAddH, err, pStrErr);

exit_function:;

    return at;
}

/****************************************************************************/
void calculate_valences(MOL_FMT_DATA *mfdata,
                        inp_ATOM *at,
                        int *num_atoms,
                        int bDoNotAddH,
                        int *err,
                        char *pStrErr)
{
    int bNonMetal;
    int a1, a2, n1, n2, valence;
    AT_NUMB *p1;

    /* special valences */

    for (bNonMetal = 0; bNonMetal < 2; bNonMetal++)
    {
        for (a1 = 0; a1 < *num_atoms; a1++)
        {
            int num_bond_type[MAX_INPUT_BOND_TYPE - MIN_INPUT_BOND_TYPE + 1],
                bond_type,
                bHasMetalNeighbor;
            /* should the "!=" be replaced with "==" ??? */
            if (bNonMetal == is_el_a_metal(at[a1].el_number))
            {
                /* first process all metals, after that all non-metals */
                continue;
            }

            memset(num_bond_type, 0, sizeof(num_bond_type)); /* djb-rwth: memset_s C11/Annex K variant? */

            /* valence = at[a1].chem_bonds_valence; */ /*  save atom valence if available */
                                                       /* 2006-08-31: fix for uncharged >N(IV)- in an aromatic ring */

            valence =
                (mfdata && mfdata->ctab.atoms) ? mfdata->ctab.atoms[a1].valence
                                               : at[a1].chem_bonds_valence;

            at[a1].chem_bonds_valence = 0;
            bHasMetalNeighbor = 0;
            for (n1 = 0; n1 < at[a1].valence; n1++)
            {
                bond_type = at[a1].bond_type[n1] - MIN_INPUT_BOND_TYPE;
                if (bond_type < 0 || bond_type > MAX_INPUT_BOND_TYPE - MIN_INPUT_BOND_TYPE)
                {
                    bond_type = 0;
                    TREAT_ERR(*err, 0, "Unknown bond type in MOLfile assigned as a single bond");
                }
                num_bond_type[bond_type]++;
                /* -- too a radical solution -- removed from next to ver 1.12B --- */
            }

            for (n1 = 0;
                 MIN_INPUT_BOND_TYPE + n1 <= 3 &&
                 MIN_INPUT_BOND_TYPE + n1 <= MAX_INPUT_BOND_TYPE;
                 n1++)
            {
                /* add all bond orders except for "aromatic" bonds */
                at[a1].chem_bonds_valence += (MIN_INPUT_BOND_TYPE + n1) * num_bond_type[n1];
            }

            n2 = 0; /* djb-rwth: ignoring LLVM warning: value used */
            if (MIN_INPUT_BOND_TYPE <= BOND_TYPE_ALTERN &&
                BOND_TYPE_ALTERN <= MAX_INPUT_BOND_TYPE &&
                (n2 = num_bond_type[BOND_TYPE_ALTERN - MIN_INPUT_BOND_TYPE]))
            {
                /* accept input aromatic bonds for now */
                switch (n2)
                {
                case 2:
                    at[a1].chem_bonds_valence += 3; /* =A- */
                    break;
                case 3:
                    at[a1].chem_bonds_valence += 4; /* =A< */
                    break;
                default:
                    /*  if 1 or >= 4 aromatic bonds then replace   */
                    /* such bonds with single bonds                */
                    /*  and detect an error in the input structure */
                    for (n1 = 0; n1 < at[a1].valence; n1++)
                    {
                        if (at[a1].bond_type[n1] == BOND_TYPE_ALTERN)
                        {
                            a2 = at[a1].neighbor[n1];
                            p1 = is_in_the_list(at[a2].neighbor, (AT_NUMB)a1,
                                                at[a2].valence);
                            if (p1)
                            {
                                at[a1].bond_type[n1] =
                                    at[a2].bond_type[p1 - at[a2].neighbor] = BOND_TYPE_SINGLE;
                            }
                            else
                            {
                                *err = -2; /*  Program error */
                                TREAT_ERR(*err, 0, "Program error interpreting MOLfile");
                                return; /*  no structure */
                            }
                        }
                    }

                    at[a1].chem_bonds_valence += n2;
                    *err |= 32;
                    TREAT_ERR(*err, 0, "Atom has 1 or more than 3 aromatic bonds");
                    n2 = 0;
                    break;
                }
            }

            if (n2 && !valence)
            {
                /* atom has aromatic bonds AND the chemical valence is not known */

                int num_H = NUMH(at, a1);
                /* bug fix 2006-08-25: aliased H result in num_H > 0 */
                /* => wrong call to detect_unusual_el_valence() */
                int chem_valence = at[a1].chem_bonds_valence /*+ num_H*/;

                int bUnusualValenceArom =
                    detect_unusual_el_valence((int)at[a1].el_number, at[a1].charge,
                                              at[a1].radical, chem_valence,
                                              num_H, at[a1].valence);
                int bUnusualValenceNoArom =
                    detect_unusual_el_valence((int)at[a1].el_number, at[a1].charge,
                                              at[a1].radical, chem_valence - 1,
                                              num_H, at[a1].valence);

#if (CHECK_AROMBOND2ALT == 1)
                if (bUnusualValenceArom &&
                    !bUnusualValenceNoArom &&
                    0 == nBondsValToMetal(at, a1))
#else
                if (bUnusualValenceArom && !bUnusualValenceNoArom)
#endif
                {
                    /* typically NH in 5-member aromatic ring */
                    at[a1].chem_bonds_valence--;
                }
            }
            else if (n2 && valence)
            {
                /* atom has aromatic bonds AND the chemical valence is known */
                int num_H = NUMH(at, a1);
                int chem_valence = at[a1].chem_bonds_valence + num_H;
                if (valence == chem_valence - 1)
                {
                    /* typically NH in 5-member aromatic ring */
                    at[a1].chem_bonds_valence--;
                }
            }

            /*  Set number of hydrogen atoms */
            if (mfdata)
            {
                at[a1].num_H = get_num_H(at[a1].elname,
                                         at[a1].num_H,
                                         at[a1].num_iso_H,
                                         at[a1].charge, at[a1].radical,
                                         at[a1].chem_bonds_valence,
                                         mfdata->ctab.atoms[a1].valence, /* instead of valence */
                                         mfdata->ctab.atoms[a1].atom_aliased_flag,
                                         bDoNotAddH,
                                         bHasMetalNeighbor);
            }
        }
    } /* for ( bNonMetal = ... */

    return;
}

/****************************************************************************/
int SetInpAtomsXYZ(MOL_FMT_DATA *mfdata,
                   int num_atoms,
                   inp_ATOM *at,
                   int *err,
                   char *pStrErr)
{
    int i, num_dimensions = 0;

#if (NORMALIZE_INP_COORD == 1)
    int do_scale_xyz = 1;
#else
    int do_scale_xyz = 0;
#endif
    double x0, y0, z0, xmin, ymin, zmin, scaler;

    num_dimensions = MolfileGetXYZDimAndNormFactors(mfdata, do_scale_xyz,
                                                    &x0, &y0, &z0,
                                                    &xmin, &ymin, &zmin,
                                                    &scaler, err, pStrErr);

    if (num_dimensions == 0)
    {
        goto exit_function;
    }

    for (i = 0; i < num_atoms; i++)
    {

        double x = mfdata->ctab.atoms[i].fx;
        double y = mfdata->ctab.atoms[i].fy;
        double z = mfdata->ctab.atoms[i].fz;

        if (!do_scale_xyz)
        {
            at[i].x = x;
            at[i].y = y;
            at[i].z = z;
        }
        else
        {
            x = (x - xmin) * scaler + x0;
            y = (y - ymin) * scaler + y0;
            z = (z - zmin) * scaler + z0;
            /* floor() behavior is not well defined for negative arguments.
             * Use positive arguments only to get nearest integer.
             */
            at[i].x = (x >= 0.0) ? (int)floor(x + 0.5) : -(int)floor(-x + 0.5);
            at[i].y = (y >= 0.0) ? (int)floor(y + 0.5) : -(int)floor(-y + 0.5);
            at[i].z = (z >= 0.0) ? (int)floor(z + 0.5) : -(int)floor(-z + 0.5);
        }
    }

exit_function:;

    return num_dimensions;
}

/****************************************************************************/
inp_ATOM *CreateInpAtom(int num_atoms)
{
    void *p = NULL;
    if (num_atoms >= 0)
    {
        p = inchi_calloc(num_atoms, sizeof(inp_ATOM));
    }

    // void *p = inchi_calloc(num_atoms, sizeof(inp_ATOM));
    return (inp_ATOM *)p;
}

/****************************************************************************/
void FreeInpAtom(inp_ATOM **at)
{
    if (at && *at)
    {
        inchi_free(*at);
        *at = NULL;
    }

    return;
}

/****************************************************************************/
void FreeInpAtomData(INP_ATOM_DATA *inp_at_data)
{
    if (inp_at_data)
    {
        FreeInpAtom(&inp_at_data->at);
        FreeInpAtom(&inp_at_data->at_fixed_bonds);
        memset(inp_at_data, 0, sizeof(*inp_at_data)); /* djb-rwth: memset_s C11/Annex K variant? */
    }

    return;
}

/****************************************************************************/
int CreateInpAtomData(INP_ATOM_DATA *inp_at_data,
                      int num_atoms,
                      int create_at_fixed_bonds)
{
    FreeInpAtomData(inp_at_data);

    if ((inp_at_data->at = CreateInpAtom(num_atoms)) &&
        (!create_at_fixed_bonds || (inp_at_data->at_fixed_bonds = CreateInpAtom(num_atoms))))
    {
        inp_at_data->num_at = num_atoms;
        return 1;
    }

    FreeInpAtomData(inp_at_data);

    return 0;
}

/****************************************************************************/
void FreeCompAtomData(COMP_ATOM_DATA *inp_at_data)
{
    FreeInpAtom(&inp_at_data->at);

    if (inp_at_data->nOffsetAtAndH)
    {
        inchi_free(inp_at_data->nOffsetAtAndH);
    }
    memset(inp_at_data, 0, sizeof(*inp_at_data)); /* djb-rwth: memset_s C11/Annex K variant? */
}

#ifndef TARGET_API_LIB

/****************************************************************************/
int CreateCompAtomData(COMP_ATOM_DATA *inp_at_data,
                       int num_atoms,
                       int num_components,
                       int bIntermediateTaut)
{
    FreeCompAtomData(inp_at_data);

    if ((inp_at_data->at = CreateInpAtom(num_atoms)) &&
        (num_components <= 1 || bIntermediateTaut ||
         (inp_at_data->nOffsetAtAndH = (AT_NUMB *)inchi_calloc(2 * ((long long)num_components + 1), sizeof(inp_at_data->nOffsetAtAndH[0]))))) /* djb-rwth: cast operator added */
    {
        inp_at_data->num_at = num_atoms;
        inp_at_data->num_components = (num_components > 1) ? num_components : 0;
        return 1;
    }

    FreeCompAtomData(inp_at_data);

    return 0;
}
#endif

#ifndef COMPILE_ANSI_ONLY

/****************************************************************************/
void FreeInfAtom(inf_ATOM **at)
{
    if (at && *at)
    {
        inchi_free(*at);
        *at = NULL;
    }
    return;
}

/****************************************************************************/
inf_ATOM *CreateInfAtom(int num_atoms)
{
    return (inf_ATOM *)inchi_calloc(num_atoms, sizeof(inf_ATOM));
}

/****************************************************************************/
void FreeInfoAtomData(INF_ATOM_DATA *inf_at_data)
{
    FreeInfAtom(&inf_at_data->at);
    if (inf_at_data->pStereoFlags)
    {
        inchi_free(inf_at_data->pStereoFlags);
    }
    memset(inf_at_data, 0, sizeof(*inf_at_data)); /* djb-rwth: memset_s C11/Annex K variant? */

    return;
}

/****************************************************************************/
int CreateInfoAtomData(INF_ATOM_DATA *inf_at_data,
                       int num_atoms,
                       int num_components)
{
    FreeInfoAtomData(inf_at_data);

    memset(inf_at_data, 0, sizeof(*inf_at_data)); /* djb-rwth: memset_s C11/Annex K variant? */

    if ((inf_at_data->at = CreateInfAtom(num_atoms)) &&
        (num_components <= 1 ||
         (inf_at_data->pStereoFlags = (AT_NUMB *)inchi_calloc((long long)num_components + 1, sizeof(inf_at_data->pStereoFlags[0]))) /* djb-rwth: cast operator added */
         ))
    {
        inf_at_data->num_at = num_atoms;
        inf_at_data->num_components = num_components;
        return 1;
    }

    FreeInfoAtomData(inf_at_data);

    return 0;
}

/****************************************************************************/
int AllocateInfoAtomData(INF_ATOM_DATA *inf_at_data,
                         int num_atoms,
                         int num_components)
{
    if ((inf_at_data->at = CreateInfAtom(num_atoms))) /* djb-rwth: addressing LLVM warning */
    {
        if (num_components > 1 &&
            !(inf_at_data->pStereoFlags = (AT_NUMB *)inchi_calloc((long long)num_components + 1, sizeof(inf_at_data->pStereoFlags[0])))) /* djb-rwth: cast operator added */
        {
            FreeInfAtom(&inf_at_data->at);
            return 0;
        }
        return 1;
    }

    return 0;
}

/****************************************************************************/
int DuplicateInfoAtomData(INF_ATOM_DATA *inf_at_data_to,
                          const INF_ATOM_DATA *inf_at_data_from)
{
    *inf_at_data_to = *inf_at_data_from;

    if (AllocateInfoAtomData(inf_at_data_to, inf_at_data_from->num_at, inf_at_data_from->num_components))
    {
        memcpy(inf_at_data_to->at, inf_at_data_from->at,
               inf_at_data_from->num_at * sizeof(inf_at_data_to->at[0]));
        if (inf_at_data_to->pStereoFlags && inf_at_data_from->pStereoFlags)
        {
            memcpy(inf_at_data_to->pStereoFlags, inf_at_data_from->pStereoFlags,
                   ((long long)inf_at_data_from->num_components + 1) * sizeof(inf_at_data_to->pStereoFlags[0])); /* djb-rwth: cast operator added */
        }
        return 1;
    }

    return 0;
}
#endif /* COMPILE_ANSI_ONLY */

/****************************************************************************/
void FreeOrigAtData(ORIG_ATOM_DATA *orig_at_data)
{
    if (!orig_at_data)
    {
        return;
    }

    FreeInpAtom(&orig_at_data->at);

    if (orig_at_data->nCurAtLen)
    {
        inchi_free(orig_at_data->nCurAtLen);
    }

    if (orig_at_data->nOldCompNumber)
    {
        inchi_free(orig_at_data->nOldCompNumber);
    }

    if (orig_at_data->szCoord)
    {
        inchi_free(orig_at_data->szCoord);
    }

    if (orig_at_data->nEquLabels)
    {
        inchi_free(orig_at_data->nEquLabels);
    }

    if (orig_at_data->nSortedOrder)
    {
        inchi_free(orig_at_data->nSortedOrder);
    }

    /* v 1.05 */
    FreeExtOrigAtData(orig_at_data->polymer, orig_at_data->v3000);

    memset(orig_at_data, 0, sizeof(*orig_at_data)); /* djb-rwth: memset_s C11/Annex K variant? */

    return;
}

/****************************************************************************
Free v. 1.05 extensions stuff
****************************************************************************/
void FreeExtOrigAtData(OAD_Polymer *pd, OAD_V3000 *v3k)
{
    int k;

    OAD_Polymer_Free(pd);
    pd = NULL;

    if (v3k)
    {
        if (v3k->atom_index_orig)
        {
            inchi_free(v3k->atom_index_orig);
            v3k->atom_index_orig = NULL;
        }
        if (v3k->atom_index_fin)
        {
            inchi_free(v3k->atom_index_fin);
            v3k->atom_index_fin = NULL;
        }
        if (v3k->n_haptic_bonds && v3k->lists_haptic_bonds)
        {
            for (k = 0; k < v3k->n_haptic_bonds; k++)
                if (v3k->lists_haptic_bonds[k])
                {
                    inchi_free(v3k->lists_haptic_bonds[k]);
                    v3k->lists_haptic_bonds[k] = NULL;
                }
            inchi_free(v3k->lists_haptic_bonds);
            v3k->lists_haptic_bonds = NULL;
        }
        if (v3k->n_steabs && v3k->lists_steabs)
        {
            for (k = 0; k < v3k->n_steabs; k++)
                if (v3k->lists_steabs[k])
                {
                    inchi_free(v3k->lists_steabs[k]);
                    v3k->lists_steabs[k] = NULL;
                }
            inchi_free(v3k->lists_steabs);
            v3k->lists_steabs = NULL;
        }
        if (v3k->n_sterel && v3k->lists_sterel)
        {
            for (k = 0; k < v3k->n_sterel; k++)
                if (v3k->lists_sterel[k])
                {
                    inchi_free(v3k->lists_sterel[k]);
                    v3k->lists_sterel[k] = NULL;
                }
            inchi_free(v3k->lists_sterel);
            v3k->lists_sterel = NULL;
        }
        if (v3k->n_sterac && v3k->lists_sterac)
        {
            for (k = 0; k < v3k->n_sterac; k++)
                if (v3k->lists_sterac[k])
                {
                    inchi_free(v3k->lists_sterac[k]);
                    v3k->lists_sterac[k] = NULL;
                }
            inchi_free(v3k->lists_sterac);
            v3k->lists_sterac = NULL;
        }
        memset(v3k, 0, sizeof(*v3k)); /* djb-rwth: memset_s C11/Annex K variant? */
        inchi_free(v3k);
    }

    return;
}

/****************************************************************************/
int SetExtOrigAtDataByMolfileExtInput(MOL_FMT_DATA *mfdata,
                                      OAD_Polymer **ppPolymer,
                                      OAD_V3000 **ppV3000,
                                      char *pStrErr)
{
    int k, m, err = 0;
    OAD_V3000 *pv = NULL;
    int nsgroups = mfdata->ctab.sgroups.used;

    /* djb-rwth: addressing coverity ID #499476 -- TREAT_ERR properly used in all cases */
    /* Polymers */
    if (nsgroups > 0)
    {
        /* Prepare OAD_Polymer container */
        *ppPolymer = (OAD_Polymer *)inchi_calloc(1, sizeof(OAD_Polymer));
        if (!(*ppPolymer))
        {
            TREAT_ERR(err, 9001, "Out of RAM");
            goto exit_function;
        }

        /* Convert Molfile's Sgroup's to OAD_PolymerUnit's */
        (*ppPolymer)->units = (OAD_PolymerUnit **)inchi_calloc(nsgroups, sizeof((*ppPolymer)->units[0]));
        if (!(*ppPolymer)->units)
        {
            TREAT_ERR(err, 9001, "Out of RAM");
            goto exit_function;
        }
        memset((*ppPolymer)->units, 0, sizeof(*(*ppPolymer)->units)); /* djb-rwth: memset_s C11/Annex K variant? */

        (*ppPolymer)->n = nsgroups;
        (*ppPolymer)->is_in_reconn = 0;
        (*ppPolymer)->pzz = NULL;
        (*ppPolymer)->edit_repeats = 0;
        (*ppPolymer)->really_do_frame_shift = 0;
        (*ppPolymer)->frame_shift_scheme = FSS_STARS_CYCLED;
        (*ppPolymer)->treat = POLYMERS_MODERN;

        for (k = 0; k < nsgroups; k++)
        {
            int q = 0;
            MOL_FMT_SGROUP *groupk = mfdata->ctab.sgroups.group[k];

            OAD_PolymerUnit *unitk = (*ppPolymer)->units[k] = (OAD_PolymerUnit *)inchi_calloc(1, sizeof(OAD_PolymerUnit));

            if (!unitk)
            {
                TREAT_ERR(err, 9001, "Out of RAM");
                goto exit_function;
            }

            memset(unitk, 0, sizeof(*unitk)); /* djb-rwth: memset_s C11/Annex K variant? */
            unitk->id = groupk->id;
            unitk->type = groupk->type;
            unitk->subtype = groupk->subtype;
            unitk->conn = groupk->conn;
            unitk->label = groupk->label;

            for (q = 0; q < 4; q++)
            {
                unitk->xbr1[q] = groupk->xbr1[q];
                unitk->xbr2[q] = groupk->xbr2[q];
            }
            strcpy(unitk->smt, groupk->smt);
            unitk->na = groupk->alist.used;
            unitk->alist = (int *)inchi_calloc(unitk->na, sizeof(int));
            if (!unitk->alist)
            {
                TREAT_ERR(err, 9001, "Out of RAM");
                goto exit_function;
            }
            for (m = 0; m < unitk->na; m++)
            {
                unitk->alist[m] = groupk->alist.item[m];
            }
            unitk->nb = groupk->blist.used;
            if (unitk->nb > 0)
            {
                unitk->blist = (int *)inchi_calloc(2 * (long long)unitk->nb, sizeof(int)); /* djb-rwth: cast operator added */
                if (!unitk->blist)
                {
                    TREAT_ERR(err, 9001, "Out of RAM");
                    goto exit_function;
                }
                for (m = 0; m < groupk->blist.used; m++)
                {
                    int ib, ia1, ia2;
                    ib = groupk->blist.item[m];
                    if (ib < 1 || ib > mfdata->ctab.n_bonds)
                    {
                        TREAT_ERR(err, 9004, "Polymer unit in Molfile refers to invalid bond");
                        goto exit_function;
                    }
                    ia1 = mfdata->ctab.bonds[ib - 1].atnum1;
                    ia2 = mfdata->ctab.bonds[ib - 1].atnum2;
                    unitk->blist[2 * m] = ia1;
                    unitk->blist[2 * m + 1] = ia2;
                    if (!strcmp(mfdata->ctab.atoms[ia1 - 1].symbol, "H") ||
                        !strcmp(mfdata->ctab.atoms[ia2 - 1].symbol, "H"))
                    {
                        TREAT_ERR(err, 9002, "Hydrogen as polymer end group is not supported");
                        goto exit_function;
                    }
                }
            }
            else
            {
                unitk->blist = NULL;
            }
        }
    }

    /* V3000 Extensions */
    if (mfdata->ctab.v3000)
    {
        int nn;
        MOL_FMT_v3000 *mpv = mfdata->ctab.v3000;

        *ppV3000 = (OAD_V3000 *)inchi_calloc(1, sizeof(OAD_V3000));
        pv = *ppV3000;
        if (!pv)
        {
            TREAT_ERR(err, 9001, "Out of RAM");
            goto exit_function;
        }
        memset(pv, 0, sizeof(*pv)); /* djb-rwth: memset_s C11/Annex K variant? */

        pv->n_collections = mpv->n_collections;
        pv->n_haptic_bonds = mpv->n_haptic_bonds;
        pv->n_non_haptic_bonds = mpv->n_non_haptic_bonds;
        pv->n_sgroups = mpv->n_sgroups;
        pv->n_non_star_atoms = mpv->n_non_star_atoms;
        pv->n_star_atoms = mpv->n_star_atoms;
        pv->n_steabs = mpv->n_steabs;
        pv->n_sterac = mpv->n_sterac;
        pv->n_sterel = mpv->n_sterel;
        pv->n_3d_constraints = mpv->n_3d_constraints;

        if (mpv->atom_index_orig)
        {
            pv->atom_index_orig = (int *)inchi_calloc(mfdata->ctab.n_atoms, sizeof(int));
            if (NULL == pv->atom_index_orig)
            {
                TREAT_ERR(err, 9001, "Out of RAM");
                goto exit_function;
            }
            memcpy(pv->atom_index_orig, mpv->atom_index_orig, mfdata->ctab.n_atoms);
        }
        if (mpv->atom_index_fin)
        {
            pv->atom_index_fin = (int *)inchi_calloc(mfdata->ctab.n_atoms, sizeof(int));
            if (NULL == pv->atom_index_fin)
            {
                TREAT_ERR(err, 9001, "Out of RAM");
                goto exit_function;
            }
            memcpy(pv->atom_index_fin, mpv->atom_index_fin, mfdata->ctab.n_atoms);
        }
        if (mpv->n_haptic_bonds && mpv->haptic_bonds)
        {
            pv->lists_haptic_bonds = (int **)inchi_calloc(mpv->n_haptic_bonds, sizeof(int *));
            if (NULL == pv->lists_haptic_bonds)
            {
                TREAT_ERR(err, 9001, "Out of RAM");
                goto exit_function;
            }
            for (m = 0; m < mpv->n_haptic_bonds; m++)
            {
                int *lst = NULL;
                int *mol_lst = mpv->haptic_bonds->lists[m];
                nn = mol_lst[2] + 3;
                lst = pv->lists_haptic_bonds[m] = (int *)inchi_calloc(nn, sizeof(int));
                if (NULL == lst)
                {
                    TREAT_ERR(err, 9001, "Out of RAM");
                    goto exit_function;
                }
                for (k = 0; k < nn; k++)
                {
                    lst[k] = mol_lst[k];
                }
            }
        }
        if (mpv->n_steabs && mpv->steabs)
        {
            pv->lists_steabs = (int **)inchi_calloc(mpv->n_steabs, sizeof(int *));
            if (NULL == pv->lists_steabs)
            {
                TREAT_ERR(err, 9001, "Out of RAM");
                goto exit_function;
            }
            for (m = 0; m < mpv->n_steabs; m++)
            {
                int *lst = NULL;
                int *mol_lst = mpv->steabs->lists[m];
                nn = mol_lst[1] + 2;
                lst = pv->lists_steabs[m] = (int *)inchi_calloc(nn, sizeof(int));
                if (NULL == lst)
                {
                    TREAT_ERR(err, 9001, "Out of RAM");
                    goto exit_function;
                }
                for (k = 0; k < nn; k++)
                {
                    lst[k] = mol_lst[k];
                }
            }
        }
        if (mpv->n_sterac && mpv->sterac)
        {
            pv->lists_sterac = (int **)inchi_calloc(mpv->n_sterac, sizeof(int *));
            if (NULL == pv->lists_sterac)
            {
                TREAT_ERR(err, 9001, "Out of RAM");
                goto exit_function;
            }
            for (m = 0; m < mpv->n_sterac; m++)
            {
                int *lst = NULL;
                int *mol_lst = mpv->sterac->lists[m];
                nn = mol_lst[1] + 2;
                lst = pv->lists_sterac[m] = (int *)inchi_calloc(nn, sizeof(int));
                if (NULL == lst)
                {
                    TREAT_ERR(err, 9001, "Out of RAM");
                    goto exit_function;
                }
                for (k = 0; k < nn; k++)
                {
                    lst[k] = mol_lst[k];
                }
            }
        }
        if (mpv->n_sterel && mpv->sterel)
        {
            pv->lists_sterel = (int **)inchi_calloc(mpv->n_sterel, sizeof(int *));
            if (NULL == pv->lists_sterel)
            {
                TREAT_ERR(err, 9001, "Out of RAM");
                goto exit_function;
            }
            for (m = 0; m < mpv->n_sterel; m++)
            {
                int *lst = NULL;
                int *mol_lst = mpv->sterel->lists[m];
                nn = mol_lst[1] + 2;
                lst = pv->lists_sterel[m] = (int *)inchi_calloc(nn, sizeof(int));
                if (NULL == lst)
                {
                    TREAT_ERR(err, 9001, "Out of RAM");
                    goto exit_function;
                }
                for (k = 0; k < nn; k++)
                {
                    lst[k] = mol_lst[k];
                }
            }
        }
    }

exit_function:
    if (err)
    {
        FreeExtOrigAtData((*ppPolymer), pv);
        *ppPolymer = NULL;
    }

    return err;
}

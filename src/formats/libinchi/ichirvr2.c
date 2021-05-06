/*
 * International Chemical Identifier (InChI)
 * Version 1
 * Software version 1.04
 * September 9, 2011
 *
 * The InChI library and programs are free software developed under the
 * auspices of the International Union of Pure and Applied Chemistry (IUPAC).
 * Originally developed at NIST. Modifications and additions by IUPAC 
 * and the InChI Trust.
 *
 * IUPAC/InChI-Trust Licence for the International Chemical Identifier (InChI) 
 * Software version 1.0.
 * Copyright (C) IUPAC and InChI Trust Limited
 * 
 * This library is free software; you can redistribute it and/or modify it under the 
 * terms of the IUPAC/InChI Trust Licence for the International Chemical Identifier 
 * (InChI) Software version 1.0; either version 1.0 of the License, or 
 * (at your option) any later version.
 * 
 * This library is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
 * See the IUPAC/InChI Trust Licence for the International Chemical Identifier (InChI) 
 * Software version 1.0 for more details.
 * 
 * You should have received a copy of the IUPAC/InChI Trust Licence for the 
 * International Chemical Identifier (InChI) Software version 1.0 along with 
 * this library; if not, write to:
 * 
 * The InChI Trust
 * c/o FIZ CHEMIE Berlin
 * Franklinstrasse 11
 * 10587 Berlin
 * GERMANY
 * 
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*^^^ */
/*#define CHECK_WIN32_VC_HEAP*/
#include "mode.h"

#if ( READ_INCHI_STRING == 1 )

#include "ichi.h"
#include "ichitime.h"

#include "inpdef.h"
#include "ichimain.h"
#include "ichierr.h"
#include "incomdef.h" 
#include "ichiring.h"
#include "extr_ct.h"
#include "ichitaut.h"
#include "ichinorm.h"
#include "util.h"

#include "ichicomp.h"
#include "ichister.h"

#include "ichi_bns.h"

#include "strutil.h"

#include "ichirvrs.h"

/******************************************************************************************************/
void CopyAt2St( inp_ATOM *at, inp_ATOM_STEREO * st, int num_atoms )
{
    int i;
    for ( i = 0; i < num_atoms; i ++ ) {
        if ( at[i].p_parity ) {
            memcpy( st[i].p_orig_at_num, at[i].p_orig_at_num, sizeof(st[0].p_orig_at_num) );
            st[i].p_parity = at[i].p_parity;
        }
        if ( at[i].sb_parity[0] ) {
            memcpy( st[i].sb_ord, at[i].sb_ord, sizeof(st[0].sb_ord) );
            memcpy( st[i].sb_parity, at[i].sb_parity, sizeof(st[0].sb_parity) );
            memcpy( st[i].sn_ord, at[i].sn_ord, sizeof(st[0].sn_ord) );
            memcpy( st[i].sn_orig_at_num, at[i].sn_orig_at_num, sizeof(st[0].sn_orig_at_num) );
        }
    }
}
void CopySt2At( inp_ATOM *at, inp_ATOM_STEREO * st, int num_atoms )
{
    int i;
    if ( !st ) {
        return;
    }
    for ( i = 0; i < num_atoms; i ++ ) {
        if ( st[i].p_parity ) {
            memcpy( at[i].p_orig_at_num, st[i].p_orig_at_num, sizeof(at[0].p_orig_at_num) );
            at[i].p_parity = st[i].p_parity;
        }
        if ( st[i].sb_parity[0] ) {
            memcpy( at[i].sb_ord, st[i].sb_ord, sizeof(st[0].sb_ord) );
            memcpy( at[i].sb_parity, st[i].sb_parity, sizeof(at[0].sb_parity) );
            memcpy( at[i].sn_ord, st[i].sn_ord, sizeof(at[0].sn_ord) );
            memcpy( at[i].sn_orig_at_num, st[i].sn_orig_at_num, sizeof(at[0].sn_orig_at_num) );
        }
    }
}

/******************************************************************************************************/
int RestoreAtomConnectionsSetStereo( StrFromINChI *pStruct, int iComponent, int iAtNoOffset, INChI *pInChI, INChI *pInChIMobH)
{
    inp_ATOM     *at  = NULL;
    inp_ATOM_STEREO * st = NULL;
    int           num_atoms, i, jv, jn, n_vertex, n_neigh, num_H, parity;
    int           nNumDeletedH=0, iDeletedH=0, idelH1, idelH2, ret = 0, len;
    int           num_stereo_bonds, num_stereo_centers, num_stereo_bonds2, num_stereo_centers2;
    INChI_Stereo *pStereo = NULL, *pStereo2 = NULL;
    AT_NUMB       nCumulene[MAX_CUMULENE_LEN+2];

    num_atoms = pInChI->nNumberOfAtoms;
    if ( num_atoms <= 0 ) {
        return 0;
    }
    INCHI_HEAPCHK
    /* atoms */
    pStruct->at = at = (inp_ATOM *) inchi_calloc ( num_atoms, sizeof(pStruct->at[0]) );
    if ( !at ) {
        ret = RI_ERR_ALLOC;
        goto exit_function;
    }
    pStruct->num_atoms = num_atoms;
    /* charge */
    pStruct->charge = pInChI->nTotalCharge;
    /* elements, terminal atoms H */
    for ( i = 0; i < num_atoms; i ++ ) {
        at[i].el_number = pInChI->nAtom[i];
        if ( GetElementFormulaFromAtNum(UCINT pInChI->nAtom[i], at[i].elname ) ) {
            ret = RI_ERR_PROGR;
            goto exit_function;
        }
        at[i].orig_at_number = iAtNoOffset + i+1;
        at[i].orig_compt_at_numb = i + 1;
        at[i].component = iComponent + 1;
        num_H = pInChI->nNum_H[i];
        /* --- pInChI->nNum_H_fixed[i] was added to pInChI->nNum_H[i] ---
        if ( pInChI->nNum_H_fixed ) {
            num_H += pInChI->nNum_H_fixed[i];
        }
        */
        at[i].num_H = num_H;
    }
    INCHI_HEAPCHK
    /* connections */
    for ( i = 1, n_vertex = pInChI->nConnTable[0]-1; i < pInChI->lenConnTable; i ++ ) {
        if ( (n_neigh = pInChI->nConnTable[i]-1) < n_vertex ) {
            /*  vertex - neighbor connection */
            jv = at[n_vertex].valence ++;
            at[n_vertex].neighbor[jv] = n_neigh;
            at[n_vertex].bond_type[jv] = BOND_TYPE_SINGLE;
            at[n_vertex].chem_bonds_valence += at[n_vertex].bond_type[jv];
            /*  neighbor - vertex connection */
            jn = at[n_neigh].valence ++;
            at[n_neigh].neighbor[jn] = n_vertex;
            at[n_neigh].bond_type[jn] = BOND_TYPE_SINGLE;
            at[n_neigh].chem_bonds_valence += at[n_neigh].bond_type[jn];
        } else
        if ( (n_vertex = n_neigh) >= num_atoms ) {
            ret = RI_ERR_PROGR;
            goto exit_function;
        }
    }
    INCHI_HEAPCHK
    /* isotopic atoms */
    if ( pInChI->IsotopicAtom && pInChI->nNumberOfIsotopicAtoms ) {
        for ( i = 0; i < pInChI->nNumberOfIsotopicAtoms; i ++ ) {
            n_vertex = pInChI->IsotopicAtom[i].nAtomNumber-1;
            at[n_vertex].iso_atw_diff = (char)pInChI->IsotopicAtom[i].nIsoDifference;
            at[n_vertex].num_iso_H[0] = (char)pInChI->IsotopicAtom[i].nNum_H;
            at[n_vertex].num_iso_H[1] = (char)pInChI->IsotopicAtom[i].nNum_D;
            at[n_vertex].num_iso_H[2] = (char)pInChI->IsotopicAtom[i].nNum_T;
        }
        pStruct->bIsotopic |= 1;
    }
    INCHI_HEAPCHK
    /* tautomeric groups */
    if ( (ret = GetTgroupInfoFromInChI( &pStruct->ti, at, NULL, pInChI )) ) {
        goto exit_function;
    }

    /* coordinates: data from unused members: pInChI->IsotopicTGroup and InChI->nNumberOfIsotopicTGroups */
    if ( pInChI->IsotopicTGroup && !pInChI->nNumberOfIsotopicTGroups ) {
        pStruct->pXYZ = (XYZ_COORD *) pInChI->IsotopicTGroup;
        pInChI->IsotopicTGroup = NULL;
    }
    /* stereo */
    if ( pInChI->StereoIsotopic && 
         (pInChI->StereoIsotopic->nNumberOfStereoBonds +
          pInChI->StereoIsotopic->nNumberOfStereoCenters) ) {
        pStereo = pInChI->StereoIsotopic;
    } else
    if ( pInChI->Stereo && 
         (pInChI->Stereo->nNumberOfStereoBonds +
          pInChI->Stereo->nNumberOfStereoCenters) ) {
        pStereo = pInChI->Stereo;
    } else {
        pStereo = NULL;
    }
    /* stereo2: Mobile-H in addition to Fixed-H*/
    pStereo2 = NULL;
    if ( pInChIMobH && pInChIMobH->nNumberOfAtoms ) {
        if ( pInChIMobH->StereoIsotopic && 
             (pInChIMobH->StereoIsotopic->nNumberOfStereoBonds +
              pInChIMobH->StereoIsotopic->nNumberOfStereoCenters) ) {
            pStereo2 = pInChIMobH->StereoIsotopic;
        } else
        if ( pInChIMobH->Stereo && 
             (pInChIMobH->Stereo->nNumberOfStereoBonds +
              pInChIMobH->Stereo->nNumberOfStereoCenters) ) {
            pStereo2 = pInChIMobH->Stereo;
        }
    }
    INCHI_HEAPCHK

    num_stereo_bonds = num_stereo_bonds2 = 0;
    num_stereo_centers = num_stereo_centers2 = 0;
    /* -- have already been done in the initialization --
       iDeletedH = 0;
       nNumDeletedH = 0;
    */
    if ( pStereo || pStereo2 ) {
        /* count implicit H needed for parities and reallocate at[]; set at[n_vertex].at_type=1 for these atoms */
        int len1  = pStereo? pStereo->nNumberOfStereoCenters : 0;
        int len2 = pStereo2? pStereo2->nNumberOfStereoCenters : 0;
        int i2, diff, diff2;
        for ( i = i2 = 0; i < len1 || i2 < len2; ) {
            if ( i < len1 && i2 < len2 ) {
                diff = (int)pStereo->nNumber[i] - (int)pStereo2->nNumber[i2];
                if ( diff <= 0 ) {
                    n_vertex = pStereo->nNumber[i]-1;
                    i ++;
                    i2 += !diff;
                } else {
                    n_vertex = pStereo2->nNumber[i2]-1;
                    num_stereo_centers2 ++;
                    i2 ++;
                }
            } else
            if ( i < len1 ) {
                n_vertex = pStereo->nNumber[i]-1;
                i ++;
            } else {
                n_vertex = pStereo2->nNumber[i2]-1;
                num_stereo_centers2 ++;
                i2 ++;
            }
            /* find whether it is an allene */
            jv = at[n_vertex].neighbor[0];
            jn = at[n_vertex].neighbor[1];
            if ( at[n_vertex].valence == 2 &&
                 at[n_vertex].num_H   == 0 &&
                 bCanAtomBeMiddleAllene(at[n_vertex].elname, 0, 0) &&
                 at[jv].valence + at[jv].num_H == 3 &&
                 bCanAtomBeTerminalAllene(at[jv].elname, 0, 0)     &&
                 at[jn].valence + at[jn].num_H == 3 &&
                 bCanAtomBeTerminalAllene(at[jn].elname, 0, 0) ) {
                /* allene */
                if ( !at[jv].at_type && at[jv].num_H ) {
                    nNumDeletedH += at[jv].num_H;
                    at[jv].at_type ++;  /* H should be added as an explicit H */
                }
                if ( !at[jn].at_type && at[jn].num_H ) {
                    nNumDeletedH += at[jn].num_H;
                    at[jn].at_type ++;  /* H should be added as an explicit H */
                }
            } else {
                /* stereogenic atom - sp3 */
                if ( !at[n_vertex].at_type && at[n_vertex].num_H ) {
                    nNumDeletedH += at[n_vertex].num_H;
                    at[n_vertex].at_type ++; /* H should be added as an explicit H */
                }
            }
        }
        INCHI_HEAPCHK
        len1  = pStereo? pStereo->nNumberOfStereoBonds : 0;
        len2 = pStereo2? pStereo2->nNumberOfStereoBonds : 0;
        for ( i = i2 = 0; i < len1 || i2 < len2; ) {
            if ( i < len1 && i2 < len2 ) {
                diff  = (int)pStereo->nBondAtom1[i] - (int)pStereo2->nBondAtom1[i2];
                diff2 = (int)pStereo->nBondAtom2[i] - (int)pStereo2->nBondAtom2[i2];
                if ( diff < 0 || (diff == 0 && diff2 <= 0)) {
                    n_vertex = pStereo->nBondAtom1[i]-1;
                    n_neigh  = pStereo->nBondAtom2[i]-1;
                    i ++;
                    i2 += !diff && !diff2;
                } else {
                    n_vertex = pStereo2->nBondAtom1[i2]-1;
                    n_neigh  = pStereo2->nBondAtom2[i2]-1;
                    num_stereo_bonds2 ++;
                    i2 ++;
                }
            } else
            if ( i < len1 ) {
                n_vertex = pStereo->nBondAtom1[i]-1;
                n_neigh  = pStereo->nBondAtom2[i]-1;
                i ++;
            } else {
                n_vertex = pStereo2->nBondAtom1[i2]-1;
                n_neigh  = pStereo2->nBondAtom2[i2]-1;
                num_stereo_bonds2 ++;
                i2 ++;
            }
            if ( !is_in_the_list( at[n_vertex].neighbor, (AT_NUMB)n_neigh, at[n_vertex].valence ) ) {
                /* must be a cumulene */
                if ( !bFindCumuleneChain( at, (AT_NUMB)n_vertex, (AT_NUMB)n_neigh, nCumulene, MAX_CUMULENE_LEN+1 ) ) {
                    ret = RI_ERR_SYNTAX; /* not a cumulene */
                    goto exit_function;
                }
            }
            if ( !at[n_vertex].at_type && at[n_vertex].num_H ) {
                nNumDeletedH += at[n_vertex].num_H;
                at[n_vertex].at_type ++;  /* H should be added as an explicit H */
            }
            if ( !at[n_neigh].at_type && at[n_neigh].num_H ) {
                nNumDeletedH += at[n_neigh].num_H;
                at[n_neigh].at_type ++;   /* H should be added as an explicit H */
            }
        }
        INCHI_HEAPCHK
        if ( nNumDeletedH ) {
            /* add explicit H */
            inp_ATOM *at2 = (inp_ATOM *)inchi_calloc( num_atoms + nNumDeletedH, sizeof(at2[0]) );
            if ( !at2 ) {
                ret = RI_ERR_ALLOC;
                goto exit_function;
            }
            pStruct->num_deleted_H = nNumDeletedH;
            memcpy( at2, at, num_atoms * sizeof(at2[0]) );
            inchi_free( at );
            pStruct->at =  at = at2;
            /* fill out deleted H atom info */
            for ( i = num_atoms; i < num_atoms + nNumDeletedH; i ++ ) {
                strcpy( at[i].elname, "H" );
                at[i].el_number = EL_NUMBER_H;
                at[i].orig_at_number = iAtNoOffset + i+1;
                at[i].orig_compt_at_numb = i + 1;
                at[i].component = iComponent + 1;
            }
            /* connect deleted H */
            for( i = 0; i < num_atoms; i ++ ) {
                if ( at[i].at_type == 1 ) {
                    if ( 0 > (ret = AddExplicitDeletedH( at, i, num_atoms, &iDeletedH, &idelH1, nNumDeletedH, pStereo2 != NULL ))) {
                        goto exit_function;
                    }
                }
            }
        }
        INCHI_HEAPCHK
    }

    if ( pStereo ) {
        /* mark stereo centers, they have already been connected the added explicit H, if any */
        int bInvertedParity   = (pStereo->nCompInv2Abs == -1);
        for ( i = 0; i < pStereo->nNumberOfStereoCenters; i ++ ) {
            n_vertex = pStereo->nNumber[i]-1;
            parity   = pStereo->t_parity[i];
            if ( bInvertedParity ) {
                parity = (parity == AB_PARITY_EVEN)? AB_PARITY_ODD : (parity == AB_PARITY_ODD)? AB_PARITY_EVEN : parity;
            }
            /* find whether it is allene */
            if ( at[n_vertex].valence == 2 &&
                 at[n_vertex].num_H   == 0 &&
                 bCanAtomBeMiddleAllene(at[n_vertex].elname, 0, 0) &&
                 /* allene has exactly 2 double bonds */
                 (jv = at[n_vertex].neighbor[0], at[jv].valence + at[jv].num_H == 3) &&
                 bCanAtomBeTerminalAllene(at[jv].elname, 0, 0)     &&
                 (jn = at[n_vertex].neighbor[1], at[jn].valence + at[jn].num_H == 3) &&
                 bCanAtomBeTerminalAllene(at[jn].elname, 0, 0) ) {
                /* allene: add explicit H if implicit H are present */
                /* iDeletedH = current number of already added explicit H */
                /* idelH1    = index in at[] of the explicit H added to atom jv */
                if ( at[jv].num_H ) {
                    if ( 0 > (ret = AddExplicitDeletedH( at, jv, num_atoms, &iDeletedH, &idelH1, nNumDeletedH, pStereo2 != NULL ))) {
                        goto exit_function;
                    }
                } else {
                    /* index of the stereo atom neighbor */
                    idelH1 = at[jv].neighbor[at[jv].neighbor[0]==n_vertex];
                }
                if ( at[jn].num_H ) {
                    /* iDeletedH = current number of already added explicit H */
                    /* idelH2    = index of the explicit H added to atom jn */
                    if ( 0 > (ret = AddExplicitDeletedH( at, jn, num_atoms, &iDeletedH, &idelH2, nNumDeletedH, pStereo2 != NULL ))) {
                        goto exit_function;
                    }
                } else {
                    idelH2 = at[jn].neighbor[at[jn].neighbor[0]==n_vertex];
                }
                /* allene: set bond types to double */
                /*
                if ( 0 > (ret = set_bond_type( at, (AT_NUMB)n_vertex, (AT_NUMB)jv, BOND_TYPE_DOUBLE ) ) ||
                     0 > (ret = set_bond_type( at, (AT_NUMB)n_vertex, (AT_NUMB)jn, BOND_TYPE_DOUBLE ) ) ) {
                    goto exit_function;
                }
                */
                /* allene: make 0D parity */
                ret = set_cumulene_0D_parity( at, st, num_atoms, idelH1, jv, jn, idelH2, parity, 2 );
                if ( ret < 0 ) {
                    goto exit_function;
                }
            } else {
                /* stereogenic sp3 atom */
                if ( at[n_vertex].num_H ) {
                    if ( 0 > (ret = AddExplicitDeletedH( at, n_vertex, num_atoms, &iDeletedH, &idelH1, nNumDeletedH, pStereo2 != NULL ))) {
                        goto exit_function;
                    }
                }
                ret = set_atom_0D_parity( at, st, num_atoms, nNumDeletedH, n_vertex, parity );
                if ( ret < 0 ) {
                    goto exit_function;
                }
                num_stereo_centers ++;
            }
            if ( ret < 0 ) {
                goto exit_function;
            }
        }
        INCHI_HEAPCHK
        /* mark stereobonds */
        for ( i = 0; i < pStereo->nNumberOfStereoBonds; i ++ ) {
            jv     = pStereo->nBondAtom1[i]-1;
            jn     = pStereo->nBondAtom2[i]-1;
            parity = pStereo->b_parity[i];
            if ( !is_in_the_list( at[jv].neighbor, (AT_NUMB)jn, at[jv].valence ) ) {
                /* must be a cumulene */
                if ( !bFindCumuleneChain( at, (AT_NUMB)jv, (AT_NUMB)jn, nCumulene, MAX_CUMULENE_LEN+1 ) ) {
                    return RI_ERR_SYNTAX; /* not a cumulene */
                }
                len = MAX_CUMULENE_LEN+1;
            } else {
                /* a regular double or alt bond */
                nCumulene[0] = jv;
                nCumulene[1] = jn;
                len = 1; /* cumulene length is number of bonds, not number of atoms */
            }
            /* cumulene or double bond: add explicit H if implicit H are present */
            if ( at[jv].num_H ) {
                if ( 0 > (ret = AddExplicitDeletedH( at, jv, num_atoms, &iDeletedH, &idelH1, nNumDeletedH, pStereo2 != NULL ))) {
                    goto exit_function;
                }
            } else {
                /* double bond neighbor that has the smallest canonical number; it is either 0th or 1st */
                idelH1 = at[jv].neighbor[at[jv].neighbor[0]==nCumulene[1]];
            }
            if ( at[jn].num_H ) {
                if ( 0 > (ret = AddExplicitDeletedH( at, jn, num_atoms, &iDeletedH, &idelH2, nNumDeletedH, pStereo2 != NULL ))) {
                    goto exit_function;
                }
            } else {
                idelH2 = at[jn].neighbor[at[jn].neighbor[0]==nCumulene[len-1]];
            }
            if ( 0 > (ret = set_cumulene_0D_parity( at, st, num_atoms, idelH1, jv, jn, idelH2, parity, len )) ) {
                goto exit_function;
            }
        }
        INCHI_HEAPCHK
    }
    /* allocate memory for Mobile-H-only stereo */
    if ( num_stereo_centers2 + num_stereo_bonds2 ) {
        if ( !(st = (inp_ATOM_STEREO *)inchi_calloc( num_atoms, sizeof(st[0])))) {
            ret = RI_ERR_ALLOC;
            goto exit_function;
        }
        CopyAt2St( at, st, num_atoms );
    }
    pStruct->st = st;
    if ( num_stereo_centers2 ) {
        /* In case of Fixed-H */
        /* mark additional Mobile-H stereo centers, they have already been connected the added explicit H, if any */
        int bInvertedParity   = (pStereo2->nCompInv2Abs == -1);
        for ( i = 0; i < pStereo2->nNumberOfStereoCenters; i ++ ) {
            n_vertex = pStereo2->nNumber[i]-1;
            parity   = pStereo2->t_parity[i];
            if ( at[n_vertex].p_parity ) {
                continue; /* the parity has already been set for Fixed-H */
            }
            if ( bInvertedParity ) {
                parity = (parity == AB_PARITY_EVEN)? AB_PARITY_ODD : (parity == AB_PARITY_ODD)? AB_PARITY_EVEN : parity;
            }
            /* find whether it is allene */
            if ( at[n_vertex].valence == 2 &&
                 at[n_vertex].num_H   == 0 &&
                 bCanAtomBeMiddleAllene(at[n_vertex].elname, 0, 0) &&
                 /* allene has exactly 2 double bonds */
                 (jv = at[n_vertex].neighbor[0], at[jv].valence + at[jv].num_H == 3) &&
                 bCanAtomBeTerminalAllene(at[jv].elname, 0, 0)     &&
                 (jn = at[n_vertex].neighbor[1], at[jn].valence + at[jn].num_H == 3) &&
                 bCanAtomBeTerminalAllene(at[jn].elname, 0, 0) ) {
                /* allene: add explicit H if implicit H are present */
                /* iDeletedH = current number of already added explicit H */
                /* idelH1    = index in at[] of the explicit H added to atom jv */
                if ( at[jv].num_H ) {
                    if ( 0 > (ret = AddExplicitDeletedH( at, jv, num_atoms, &iDeletedH, &idelH1, nNumDeletedH, pStereo2 != NULL ))) {
                        goto exit_function;
                    }
                } else {
                    /* index of the stereo atom neighbor */
                    idelH1 = at[jv].neighbor[at[jv].neighbor[0]==n_vertex];
                }
                if ( at[jn].num_H ) {
                    /* iDeletedH = current number of already added explicit H */
                    /* idelH2    = index of the explicit H added to atom jn */
                    if ( 0 > (ret = AddExplicitDeletedH( at, jn, num_atoms, &iDeletedH, &idelH2, nNumDeletedH, pStereo2 != NULL ))) {
                        goto exit_function;
                    }
                } else {
                    idelH2 = at[jn].neighbor[at[jn].neighbor[0]==n_vertex];
                }
                /* allene: set bond types to double */
                /*
                if ( 0 > (ret = set_bond_type( at, (AT_NUMB)n_vertex, (AT_NUMB)jv, BOND_TYPE_DOUBLE ) ) ||
                     0 > (ret = set_bond_type( at, (AT_NUMB)n_vertex, (AT_NUMB)jn, BOND_TYPE_DOUBLE ) ) ) {
                    goto exit_function;
                }
                */
                /* allene: make 0D parity */
                ret = set_cumulene_0D_parity( at, st, num_atoms, idelH1, jv, jn, idelH2, parity, 2 );
                if ( ret < 0 ) {
                    goto exit_function;
                }
            } else {
                /* stereogenic sp3 atom */
                if ( at[n_vertex].num_H ) {
                    if ( 0 > (ret = AddExplicitDeletedH( at, n_vertex, num_atoms, &iDeletedH, &idelH1, nNumDeletedH, pStereo2 != NULL ))) {
                        goto exit_function;
                    }
                }
                ret = set_atom_0D_parity( at, st, num_atoms, nNumDeletedH, n_vertex, parity );
                if ( ret < 0 ) {
                    goto exit_function;
                }
                num_stereo_centers ++;
            }
            if ( ret < 0 ) {
                goto exit_function;
            }
        }
    }
    if ( num_stereo_bonds2 ) {
        /* In case of Fixed-H */
        /* mark additional Mobile-H stereobonds, they have already been connected the added explicit H, if any */
        for ( i = 0; i < pStereo2->nNumberOfStereoBonds; i ++ ) {
            jv     = pStereo2->nBondAtom1[i]-1;
            jn     = pStereo2->nBondAtom2[i]-1;
            parity = pStereo2->b_parity[i];
            if ( !is_in_the_list( at[jv].neighbor, (AT_NUMB)jn, at[jv].valence ) ) {
                /* must be a cumulene */
                if ( !bFindCumuleneChain( at, (AT_NUMB)jv, (AT_NUMB)jn, nCumulene, MAX_CUMULENE_LEN+1 ) ) {
                    return RI_ERR_SYNTAX; /* not a cumulene */
                }
                len = MAX_CUMULENE_LEN+1;
            } else {
                /* a regular double or alt bond */
                nCumulene[0] = jv;
                nCumulene[1] = jn;
                len = 1; /* cumulene length is number of bonds, not number of atoms */
            }
            /* cumulene or double bond: add explicit H if implicit H are present */
            if ( at[jv].num_H ) {
                if ( 0 > (ret = AddExplicitDeletedH( at, jv, num_atoms, &iDeletedH, &idelH1, nNumDeletedH, pStereo2 != NULL ))) {
                    goto exit_function;
                }
            } else {
                /* double bond neighbor that has the smallest canonical number */
                idelH1 = at[jv].neighbor[at[jv].neighbor[0]==nCumulene[1]];
            }
            if ( at[jn].num_H ) {
                if ( 0 > (ret = AddExplicitDeletedH( at, jn, num_atoms, &iDeletedH, &idelH2, nNumDeletedH, pStereo2 != NULL ))) {
                    goto exit_function;
                }
            } else {
                idelH2 = at[jn].neighbor[at[jn].neighbor[0]==nCumulene[len-1]];
            }
            if ( 0 > (ret = set_cumulene_0D_parity( at, st, num_atoms, idelH1, jv, jn, idelH2, parity, len )) ) {
                goto exit_function;
            }
        }

    }


    ret = num_atoms;

exit_function:
    return ret;
}
/*************************************************************/
int SetStereoBondTypeFor0DParity( inp_ATOM *at, int i1, int m1 )
{
    AT_NUMB nCumulene[MAX_CUMULENE_LEN+2];
    int j, n1, n2, k1, m2, ret, nLenCumulene = 0, bond_type;
    k1 = at[i1].sb_ord[m1];
    n1 = i1;
    nCumulene[nLenCumulene ++] = n1;
    do {
        n2 = at[n1].neighbor[k1]; /* next atom */
        nCumulene[nLenCumulene ++] = n2;
        for (m2 = 0; m2 < MAX_NUM_STEREO_BONDS && at[n2].sb_parity[m2]; m2 ++ ) {
            if ( n1 == at[n2].neighbor[(int)at[n2].sb_ord[m2]] ) {
                /* found the endatom */
                goto found;
            }
        }
        if ( at[n2].num_H || at[n2].valence != 2 || at[n2].endpoint ) {
            break; /* not a middle cumulene */
        }
        k1 = (at[n2].neighbor[0] == n1);
        n1 = n2;
    } while ( at[n1].valence == 2 && !at[n1].num_H && nLenCumulene < MAX_CUMULENE_LEN+2 &&
              bCanAtomBeMiddleAllene( at[n1].elname, at[n1].charge, at[n1].radical ) );
    return RI_ERR_SYNTAX; /* failed */

found:
    if ( nLenCumulene == 2 ) {
        bond_type = BOND_TYPE_STEREO; /* double bond or alternating bond */
    } else {
        bond_type = BOND_TYPE_DOUBLE; /* cumulene or allene */
    }

    for ( j = 1; j < nLenCumulene; j ++ ) {
        /* if bond_type = BOND_TYPE_DOUBLE then increments at->cham_bonds_valence: */
        /* at->cham_bonds_valence += BOND_TYPE_DOUBLE-BOND_TYPE_SINGLE */
        if ( 0 > (ret = set_bond_type( at, (AT_NUMB)nCumulene[j-1], (AT_NUMB)nCumulene[j], bond_type ) ) ) {
            return RI_ERR_PROGR; /* failed */
        }
    }
    return nLenCumulene;
}
/******************************************************************************************************/
int SetStereoBondTypesFrom0DStereo( StrFromINChI *pStruct, INChI *pInChI)
{
    INChI_Stereo *pStereo;
    inp_ATOM     *at        = pStruct->at;
    int           num_atoms = pStruct->num_atoms;
    int           i, j, num_stereo_bonds, ret; 

    if ( pInChI->StereoIsotopic && 
         (pInChI->StereoIsotopic->nNumberOfStereoBonds +
          pInChI->StereoIsotopic->nNumberOfStereoCenters) ) {
        pStereo = pInChI->StereoIsotopic;
    } else
    if ( pInChI->Stereo && 
         (pInChI->Stereo->nNumberOfStereoBonds +
          pInChI->Stereo->nNumberOfStereoCenters) ) {
        pStereo = pInChI->Stereo;
    } else {
        pStereo = NULL;
    }
    
    /************************ set bond types separately from stereo *******************/
    if ( pStereo ) {
        num_stereo_bonds = 0;
        for ( i = 0; i < num_atoms; i ++ ) {
            /* set BOND_TYPE_DOUBLE in allenes and cumulenes */
            /* set BOND_TYPE_STEREO in double bond stereo */
            for ( j = 0; j < MAX_NUM_STEREO_BONDS && at[i].sb_parity[j]; j ++ ) {
                num_stereo_bonds ++;
                if ( 0 > (ret = SetStereoBondTypeFor0DParity( at, i,  j ) ) ) {
                    goto exit_function;
                }
            }
        }
        if ( num_stereo_bonds ) {
            int num_bond_type_stereo;
            int num_bond_type_altern;
            AT_NUMB neigh;
            /* replace adjacent BOND_TYPE_STEREO with BOND_TYPE_ALTERN */
            for ( i = 0; i < num_atoms; i ++ ) {
                num_bond_type_stereo = 0;
                num_bond_type_altern = 0;
                for ( j = 0; j < at[i].valence; j ++ ) {
                    num_bond_type_stereo += ( at[i].bond_type[j] == BOND_TYPE_STEREO );
                    num_bond_type_altern += ( at[i].bond_type[j] == BOND_TYPE_ALTERN );
                }
                if ( num_bond_type_stereo + num_bond_type_altern > 1 && num_bond_type_stereo ) {
                    for ( j = 0; j < at[i].valence; j ++ ) {
                        if  ( at[i].bond_type[j] == BOND_TYPE_STEREO ) {
                            neigh = at[i].neighbor[j];
                            /* does not change at[i].chem_bond_valence in case of BOND_TYPE_ALTERN */
                            if ( 0 > (ret = set_bond_type( at, (AT_NUMB)i, neigh, BOND_TYPE_ALTERN ) ) ) {
                                goto exit_function;
                            }
                        }
                    }
                }
                /* at this point only isolated stereo bonds have type BOND_TYPE_STEREO */
            }
            /* increment at[i].chem_bonds_valence if at[i] has an altern. bond */
            /* replace BOND_TYPE_STEREO with BOND_TYPE_DOUBLE and increment   */
            /* chem_bonds_valence of the adjacent atoms */
            for ( i = 0; i < num_atoms; i ++ ) {
                num_bond_type_stereo = 0;
                num_bond_type_altern = 0;
                for ( j = 0; j < at[i].valence; j ++ ) {
                    num_bond_type_stereo += ( at[i].bond_type[j] == BOND_TYPE_STEREO );
                    num_bond_type_altern += ( at[i].bond_type[j] == BOND_TYPE_ALTERN );
                }
                if ( !num_bond_type_stereo && num_bond_type_altern ) {
                    /* an atom has only BOND_TYPE_ALTERN => adjacent BOND_TYPE_ALTERN case */
                    at[i].chem_bonds_valence += 1;
                } else
                if ( num_bond_type_stereo == 1 ) {
                    /* isolated BOND_TYPE_STEREO => replace with BOND_TYPE_DOUBLE */
                    for ( j = 0; j < at[i].valence; j ++ ) {
                        if  ( at[i].bond_type[j] == BOND_TYPE_STEREO ) {
                            neigh = at[i].neighbor[j];
                            /* replacing BOND_TYPE_STEREO with BOND_TYPE_DOUBLE */
                            /* does not change at->chem_bonds_valence */
                            if ( 0 > (ret = set_bond_type( at, (AT_NUMB)i, neigh, BOND_TYPE_DOUBLE ) ) ) {
                                goto exit_function;
                            }
                            at[i].chem_bonds_valence ++;
                            at[(int)neigh].chem_bonds_valence ++;
                        }
                    }
                } else
                if ( num_bond_type_stereo + num_bond_type_altern ) {
                    /* an atom still has both BOND_TYPE_STEREO and BOND_TYPE_ALTERN */
                    ret = RI_ERR_PROGR;
                    goto exit_function;
                }
            }
            INCHI_HEAPCHK
        }
    }
    ret = 0; /* success */
exit_function:
    return ret;
}

/******************************************************************************************************/
int CopyBnsToAtom( StrFromINChI *pStruct, BN_STRUCT  *pBNS, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups, int bAllowZeroBondOrder )
{
    int i, j, atom_charge, left_charge, charge, ret = 0, v1, nMinorder;
    int          num_at = pStruct->num_atoms;
    inp_ATOM    *at     = pStruct->at;
    ICHICONST SRM *pSrm = pStruct->pSrm;
    BNS_VERTEX  *pv;
    BNS_EDGE    *pe;
    int          chem_bonds_valence, bond_order;

    atom_charge = left_charge = 0;
    for ( i = 0; i < num_at; i ++ ) {
        pv = pBNS->vert + i;
        /* bonds */
        chem_bonds_valence = 0;
        for ( j = 0; j < at[i].valence; j ++ ) {
            pe = pBNS->edge + pv->iedge[j];
            BondFlowMaxcapMinorder( at, pVA, pSrm, i, j, NULL, &nMinorder, NULL );
            bond_order = pe->flow + nMinorder;
            if ( !bAllowZeroBondOrder && !bond_order ) {
                bond_order = 1;
            }
            chem_bonds_valence += bond_order;
            at[i].bond_type[j]  = bond_order;  /* BOND_MARK_HIGHLIGHT */
        }
        at[i].chem_bonds_valence = chem_bonds_valence;
        /* charges (both may be present resulting in zero) */
        at[i].charge = pVA[i].cInitCharge;
        if ( pVA[i].nCMinusGroupEdge ) {
            pe = pBNS->edge + pVA[i].nCMinusGroupEdge - 1;
            if ( (charge = pe->flow) ) {
                at[i].charge  -= charge;
                atom_charge   -= charge;
            }
        }
        if ( pVA[i].nCPlusGroupEdge ) {
            pe = pBNS->edge + pVA[i].nCPlusGroupEdge - 1;
            if ( (charge = pe->cap - pe->flow) ) {
                at[i].charge += charge;
                atom_charge  += charge;
            }
        }
        if ( pv->st_edge.cap > pv->st_edge.flow ) {
            at[i].radical = RADICAL_SINGLET + (pv->st_edge.cap - pv->st_edge.flow);
        }
    }
    /* find charge excess */
    for ( i = num_at; i < pBNS->num_vertices; i ++ ) {
        pv = pBNS->vert + i;
        if ( (charge = pv->st_edge.cap - pv->st_edge.flow) ) {
            if ( IS_BNS_VT_C_OR_CSUPER_GR(pv->type) ) {
                left_charge -= charge;
            } else 
            if ( IS_BNS_VT_YVCONNECTOR(pv->type) ) {
                left_charge += charge;
            }
        }
    }
    /* tautomeric H and (-) */
    for ( i = 0; i < pBNS->num_t_groups; i ++ ) {
        /* tautomeric groups are first non-atom vertices;
           order of them is same as in pTCGroups->pTCG[] */
        int num_H       = pTCGroups->pTCG[i].tg_num_H;
        int num_Minus   = pTCGroups->pTCG[i].tg_num_Minus;
        int bMinusFirst = (pTCGroups->pTCG[i].tg_RestoreFlags & TGRF_MINUS_FIRST);
        int num_at_add;
        Vertex vMinus = NO_VERTEX;
        pv  = pBNS->vert + num_at + i;  /* t-group vertex */
        if ( !(pv->type & BNS_VERT_TYPE_TGROUP) ) {
            return RI_ERR_PROGR;
        }
        if ( pTCGroups->pTCG[i].tg_set_Minus > 0 && num_Minus > 0 ) {
            vMinus = pTCGroups->pTCG[i].tg_set_Minus-1;
            num_Minus --;
        }

        if ( bMinusFirst ) {
            for ( j = 0; j < pv->num_adj_edges; j ++ ) {
                pe = pBNS->edge + pv->iedge[j];
                v1 = pe->neighbor1;
                num_at_add = pe->flow;
                if ( v1 == vMinus ) {
                    if ( num_at_add ) {
                        at[v1].charge = -1;  /* no checking at[v1].charge == 0 for now ??? */
                        num_at_add --;       /* no checking  num_at_add > 0 for now ??? */
                    } else {
                        num_Minus ++;        /* error ??? */
                    }
                    vMinus = NO_VERTEX;
                }
                if ( num_at_add > 0 ) {
                    /* atom has tautomeric attachment; do not allow =N(-) */
                    if ( num_Minus && !at[v1].charge &&
                         at[v1].valence == at[v1].chem_bonds_valence ) {
                        at[v1].charge --;
                        num_at_add --;
                        num_Minus --;
                    }
                    if ( num_at_add > 0 ) {
                        at[v1].num_H += num_at_add;
                        num_H -= num_at_add;
                        num_at_add = 0;
                    }
                }
                at[v1].endpoint = i+1;
            }
            if ( (num_H+num_Minus != pv->st_edge.cap - pv->st_edge.flow) && (num_H || num_Minus || vMinus != NO_VERTEX) ) {
                return RI_ERR_PROGR;
            }
        } else {
            for ( j = pv->num_adj_edges-1; 0 <= j; j -- ) {
                pe = pBNS->edge + pv->iedge[j];
                v1 = pe->neighbor1;
                num_at_add = pe->flow;
                if ( v1 == vMinus ) {
                    if ( num_at_add ) {
                        at[v1].charge = -1;  /* no checking at[v1].charge == 0 for now ??? */
                        num_at_add --;       /* no checking  num_at_add > 0 for now ??? */
                    } else {
                        num_Minus ++;        /* error ??? */
                    }
                    vMinus = NO_VERTEX;
                }
                if ( num_at_add > 0 ) {
                    /* atom has tautomeric attachment; do not allow =N(-) */
                    if ( num_Minus && !at[v1].charge &&
                         at[v1].valence == at[v1].chem_bonds_valence ) {
                        at[v1].charge --;
                        num_at_add --;
                        num_Minus --;
                    }
                    if ( num_at_add > 0 ) {
                        at[v1].num_H += num_at_add;
                        num_H -= num_at_add;
                        num_at_add = 0;
                    }
                }
                at[v1].endpoint = i+1;
            }
            if ( (num_H+num_Minus != pv->st_edge.cap - pv->st_edge.flow) && (num_H || num_Minus || vMinus != NO_VERTEX) ) {
                return RI_ERR_PROGR;
            }
        }
    }

    return ret;
}
/******************************************************************************************************/
int CheckBnsConsistency( StrFromINChI *pStruct, BN_STRUCT  *pBNS, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups, int bNoRad )
{
    int nOutput = 0;
#ifndef TARGET_API_LIB
#if ( bRELEASE_VERSION == 0 )
    char s[128];
    int i, j, atom_charge, left_charge, charge, excess_charge, ret = 0;
    int v1, v2, flow, tot_st_flow, tot_st_cap, num_electrons, nNumMetalAtoms;
    int       num_at = pStruct->num_atoms;
    inp_ATOM *at     = pStruct->at;
    BNS_VERTEX  *pv;
    BNS_EDGE    *pe;
#ifdef _DEBUG
    int       bDebugOutput = 0;
    bNoRad = 1;
    bDebugOutput = 1;
#endif
    /* count electrons and metals */
    num_electrons = -pTCGroups->total_charge;
    nNumMetalAtoms = 0;
    for ( i = 0; i < pTCGroups->num_tgroups; i ++ ) {
        num_electrons += pTCGroups->pTCG[i].tg_num_H;
    }
    for ( i = 0; i < num_at; i ++ ) {
        num_electrons += at[i].el_number + at[i].num_H;
        nNumMetalAtoms += pVA[i].cMetal;
    }
    /* create output string */
    sprintf( s, "%d:%d%sM%02dv%da%de%db%d* ",
                bNoRad, pTCGroups->iComponent+1, num_electrons%2?"O":"E", nNumMetalAtoms,
                pBNS->num_vertices, num_at, pBNS->num_edges, pBNS->num_bonds );


    tot_st_flow = tot_st_cap = 0;
    atom_charge = left_charge = 0;
    if ( pBNS->num_atoms != num_at ) {
        fprintf( stdout, "\n%sNum. atoms discrepancy: %d(BNS) vs. %d(at) ", s, pBNS->num_atoms, num_at);
        nOutput ++;
    }
    /* check edges */
#ifdef _DEBUG
    if ( bDebugOutput && bNoRad ) {
        fprintf( stderr, "\n\n------begin------------------------------------------------\n" );
        fprintf( stderr, "\n\fedge    cap       flow     v1   v2\n\n" );
        /*                  xxxx xxxx/xxxx xxxx/xxxx xxxx xxxx    */
    }
#endif
    for ( i = 0; i < pBNS->num_edges; i ++ ) {
        pe = pBNS->edge + i;
        v1 = pe->neighbor1;
        v2 = v1 ^ pe->neighbor12;
        if ( pe->cap < pe->flow || pe->flow < 0 ) {
            fprintf( stdout, "\n%sedge %d (%d-%d) has cap=%d flow=%d  ", s, i, v1, v2, pe->cap, pe->flow );
            nOutput ++;
        }
#ifdef _DEBUG
        if ( bDebugOutput && bNoRad ) {
            /*                xxxx  xxxx/xxxx xxxx/xxxx xxxx xxxx    */
            fprintf( stderr, "%4d %4d/%-4d %4d/%-4d %4d %4d\n", i, pe->cap, pe->cap0, pe->flow, pe->flow0, v1, v2 ); 
        }
#endif
    }


    /* check vertices */
#ifdef _DEBUG
    if ( bDebugOutput && bNoRad ) {
        fprintf( stderr, "\n\fvert   st-cap   st-flow   type iedge : neigh\n\n" );
        /*                   xxxx xxxx/xxxx xxxx/xxxx 0xXXX xxxx : xxx    */
    }
#endif
    for ( i = 0; i < pBNS->num_vertices; i ++ ) {
        pv = pBNS->vert + i;
#ifdef _DEBUG
        if ( bDebugOutput && bNoRad ) {
            /*                   xxxx xxxx/xxxx xxxx/xxxx 0xXXX xxxx : xxx    */
            int j;
            const char *s;
            char sAtom[6];
            switch( pv->type ) {
            case BNS_VERT_TYPE_ATOM:
                sprintf( sAtom, "At %-2.2s", i < num_at? at[i].elname : "??" );
                s = sAtom;
                break;
            case BNS_VERT_TYPE_ATOM | BNS_VERT_TYPE_ENDPOINT:
                s = "Endpt";
                break;
            case BNS_VT_C_POS:
                s = "(+)  ";
                break;
            case BNS_VT_C_NEG:
                s = "(-)  ";
                break;
            case BNS_VT_C_POS_C:
                s = "(+C) ";
                break;
            case BNS_VT_C_NEG_C:
                s = "(-C) ";
                break;
            case BNS_VT_C_POS_M:
                s = "(+M) ";
                break;
            case BNS_VT_C_NEG_M:
                s = "(-M) ";
                break;
            case BNS_VT_C_POS_ALL:
                s = "(+)Sg";
                break;
            case BNS_VT_C_NEG_ALL:
                s = "(-)Sg";
                break;
            case BNS_VT_M_GROUP:
                s = "M-grp";
                break;
            case BNS_VERT_TYPE__AUX | BNS_VERT_TYPE_TEMP:
                s = "ChStr";
                break;
            case BNS_VERT_TYPE__AUX:
                s = "Yconn";
                break;
            case BNS_VERT_TYPE_TGROUP:
                s = "T-grp";
                break;
            default:
                s = "Unkn.";
                break;
            }
            fprintf( stderr,    "%4d %4d/%-4d %4d/%-4d 0x%03X %5s",
                                i, pv->st_edge.cap, pv->st_edge.cap0, pv->st_edge.flow, pv->st_edge.flow0,
                                pv->type, s );
            for ( j = 0; j < pv->num_adj_edges; j ++ ) {
                fprintf( stderr, " %2d", pv->iedge[j] );
            }
            fprintf( stderr, ":" );
            for ( j = 0; j < pv->num_adj_edges; j ++ ) {
                pe = pBNS->edge + pv->iedge[j];
                fprintf( stderr, " %2d", pe->neighbor12 ^ i );
            }
            fprintf( stderr, "\n" );
        }
#endif
        tot_st_flow += pv->st_edge.flow;
        tot_st_cap  += pv->st_edge.cap;
        if ( pv->num_adj_edges > pv->max_adj_edges ) {
            fprintf( stdout, "\n%s%s %d type 0x%X \"%s\" num_edges=%d > max=%d  ", s,
                i < num_at? "atom":"vertex", i,
                pv->type, at[i].elname, pv->num_adj_edges, pv->max_adj_edges );
            nOutput ++;
        }
        if ( i < num_at ) {
            /* charge on atoms */
            charge = pVA[i].cInitCharge;
            if ( pVA[i].nCMinusGroupEdge ) {
                pe = pBNS->edge + pVA[i].nCMinusGroupEdge - 1;
                if ( pe->flow > 0 ) {
                    charge   -= pe->flow;
                }
            }
            if ( pVA[i].nCPlusGroupEdge ) {
                pe = pBNS->edge + pVA[i].nCPlusGroupEdge - 1;
                if ( pe->cap > pe->flow  ) {
                    charge  += pe->cap - pe->flow;
                }
            }
            if ( bNoRad && pv->st_edge.flow != pv->st_edge.cap ) {
                fprintf( stdout, "\n%s%s %d: type 0x%X \"%s\" unexpected st_cap=%d st_flow=%d ", s,
                    i < num_at? "atom":"vertex", i,
                    pv->type, at[i].elname, pv->st_edge.cap, pv->st_edge.flow );
                nOutput ++;
            } else
            if ( bNoRad && charge && !strcmp(at[i].elname, "C") ) {
                /* ignore carbonyls */
                if ( i == 0 && num_at == 2 && !strcmp(at[1].elname, "O") &&
                     !at[0].num_H && !at[1].num_H && !pTCGroups->total_charge) {
                    ; /* C(-)#O(+) structure */
                } else {
                    fprintf( stdout, "\n%s%s %d: type 0x%X \"%s\" charge=%d ", s,
                        i < num_at? "atom":"vertex", i,
                        pv->type, at[i].elname, charge );
                    nOutput ++;
                }
            }
            atom_charge += charge;
        } else
        if ( (charge = pv->st_edge.cap - pv->st_edge.flow) > 0 ) {
            /* excess charge */
            if ( !bNoRad && IS_BNS_VT_C_OR_CSUPER_GR(pv->type) ) {
                left_charge -= charge;
            } else 
            if ( !bNoRad && IS_BNS_VT_YVCONNECTOR(pv->type) ) {
                left_charge += charge;
            } else
            if ( !bNoRad && IS_BNS_VT_M_GR(pv->type) &&
                 0 <= (j=pTCGroups->nGroup[TCG_MeFlower3]) && 
                 i == pTCGroups->pTCG[j].nVertexNumber ) {
                ; /* additional "radical" on metal flower */
            } else
            if ( !(pv->type & BNS_VERT_TYPE_TGROUP) || bNoRad ) {
                /* t-groups before running BFS should have st_cap > st_flow */
                fprintf( stdout, "\n%s%s %d: type 0x%X unexpected st_cap=%d st_flow=%d ", s,
                    i < num_at? "atom":"vertex", i,
                    pv->type, pv->st_edge.cap, pv->st_edge.flow);
                nOutput ++;
            }
        }
        if ( pv->st_edge.cap < pv->st_edge.flow || pv->st_edge.flow < 0 ) {
            fprintf( stdout, "\n%s%s %d: type 0x%X \"%s\" st_cap=%d st_flow=%d  ", s,
                i < num_at? "atom":"vertex", i,
                pv->type, i < num_at? at[i].elname:"", pv->st_edge.cap, pv->st_edge.flow );
            nOutput ++;
        }
        /* check edge_flow vs. st_flow consistency */
        for( j = 0, flow = 0; j < pv->num_adj_edges; j ++ ) {
            pe = pBNS->edge + pv->iedge[j];
            flow += pe->flow;
        }
        if ( flow != pv->st_edge.flow ) {
            fprintf( stdout, "\n%s%s %d: type 0x%X \"%s\" st_flow=%d edge_flow=%d  ", s,
                i < num_at? "atom":"vertex", i,
                pv->type, i < num_at? at[i].elname:"", pv->st_edge.flow, flow );
            nOutput ++;
        }
    }
#ifdef _DEBUG
    if ( bDebugOutput && bNoRad ) {
        fprintf( stderr, "\n------end--------------------------------------------------\n" );
    }
#endif
    /*
    if ( num_electrons %= 2 ) {
        fprintf( stdout, "\n%d*Odd number of electrons (%d atoms) ", bNoRad, num_at );
        nOutput ++;
    }
    */
    /* tautomeric groups charge */
    for ( i = 0, charge = 0; i < pTCGroups->num_tgroups; i ++ ) {
        charge -= pTCGroups->pTCG[i].tg_num_Minus;
    }
    /* compare */
    if ( charge != pTCGroups->tgroup_charge ) {
        fprintf( stdout, "\n%sCounted t-group charge=%d while %d was saved  ", s,
               charge, pTCGroups->tgroup_charge);
        nOutput ++;
    }
    /* add other charges */
    charge += atom_charge + left_charge;
    excess_charge = pTCGroups->total_charge - pTCGroups->added_charge - pTCGroups->tgroup_charge;
    if ( charge != pTCGroups->total_charge && excess_charge != pTCGroups->total_charge - charge ) {
        fprintf( stdout, "\n%sCounted total charge=%d while %d was saved; excess charge=%d  ", s,
            charge, pTCGroups->total_charge, excess_charge );
        nOutput ++;
    }
    if ( tot_st_cap != pBNS->tot_st_cap || tot_st_flow != pBNS->tot_st_flow ) {
        fprintf( stdout, "\n%sCounted/saved total st_flow=%d/%d st_cap=%d/%d  ", s,
            tot_st_flow, pBNS->tot_st_flow, tot_st_cap, pBNS->tot_st_cap );
        nOutput ++;
    }
    if ( nOutput ) {
        fprintf( stdout, "\n" );
    }
#endif
#endif
    return nOutput;
}

/******************************************************************************************************/
int AddExplicitDeletedH( inp_ATOM *at, int jv, int num_at, int *iDeletedH, int *iH, int nNumDeletedH, int bTwoStereo )
{
    inp_ATOM  *cur_H, *cur_at = at+jv;
    int        tot_num_iso_H = NUM_ISO_H(cur_at, 0);
    int        num_H     = cur_at->num_H;
    int        iso_H     = 0;
    S_CHAR     num_iso_H[NUM_H_ISOTOPES];
    int        i;

    if ( !at[jv].at_type ) {
        return RI_ERR_PROGR;
    }

    if ( at[jv].at_type > 1 ) {
        /* explicit hydrogens have already been added; find them */
        for ( i = 0; i < *iDeletedH; i ++ ) {
            if ( at[num_at + i].neighbor[0] == jv ) {
                *iH = num_at + i; /* return the first found H, it has the smallest canonical pseudo rank */
                return 0;
            }
        }
        return RI_ERR_PROGR;
    }
    /* add all explicit H disconnected from at[jv] in order H, 1H, D, T */
    *iH = *iDeletedH + num_at; /* num_H includes all H, both isotopic and normal */
    for ( i = 0; i < NUM_H_ISOTOPES; i ++ ) {
        num_iso_H[i] = at[jv].num_iso_H[i];
    }
    for ( ; num_H && (*iDeletedH) < nNumDeletedH; (*iDeletedH) ++ ) {
        cur_H = at + num_at + (*iDeletedH); /* first available empty atom will be this explicit H */
        cur_H->neighbor[cur_H->valence] = jv; /* connect this new atom H to the real atom */
        cur_H->bond_type[cur_H->valence] = BOND_TYPE_SINGLE;
        cur_H->valence ++;
        if ( num_H > tot_num_iso_H ) {
            num_H --;
            if ( num_H != tot_num_iso_H ) {
                /* may happen when Mobile-H stereo included in Fixed-H processing */
                if ( bTwoStereo ) {
                    continue;
                } else {
                    return RI_ERR_SYNTAX; /* two identical H neighbors of a stereo atom/bond */
                }
            }
        } else {
            while ( iso_H < NUM_H_ISOTOPES && !num_iso_H[iso_H] )
                iso_H ++;
            if ( iso_H < NUM_H_ISOTOPES ) {
                cur_H->iso_atw_diff = iso_H + 1; /* isotopic shift + 1 */
                num_H --;
                tot_num_iso_H --;
                num_iso_H[iso_H] --;
                if ( num_iso_H[iso_H] ) {
                    return RI_ERR_SYNTAX; /* two identical isotopic H neighbors of a stereo atom/bond */
                }
            } else {
                return RI_ERR_SYNTAX; /* not enough isotopic H */
            }
        }
    }
    if ( num_H ) {
        return RI_ERR_SYNTAX;
    }
    at[jv].at_type ++; /* at[jv].at_type==2 => explicit hydrogens have already been added */
    return 0; /* success */    
}
/******************************************************************************************************/
int bFindCumuleneChain( inp_ATOM *at, AT_NUMB i1, AT_NUMB i2, AT_NUMB nCumulene[], int nMaxLen )
/* nMaxLen = number of bonds in cumulene = 3 = MAX_CUMULENE_LEN+1 */
/* nCumulene[nMaxLen+1] will contain cumulene chain >i1=x=y=i2< in this order */
{
    int i, len, iat, nat;
    nCumulene[0] = i1;
    for ( i = 0; i < at[i1].valence; i ++ ) {
        len = 0;
        iat = i1; /* current */
        nat = at[i1].neighbor[i]; /* next */
        if ( len+1 == nMaxLen ) {
            if ( nat == i2 ) {
                nCumulene[++len] = nat;
                return 1; /* success */
            }
            continue; /* check next at[i1] neighbor */
        }
        while ( at[nat].valence == 2 &&
                at[nat].num_H   == 0 &&
                bCanAtomBeMiddleAllene(at[nat].elname, 0, 0) ) {
            nCumulene[++len] = nat;
            nat = at[nat].neighbor[at[nat].neighbor[0]==iat]; /* new next */
            if ( len+1 == nMaxLen ) {
                if ( nat == i2 ) {
                    nCumulene[++len] = nat;
                    return 1; /* success */
                }
                break; /* check next at[i1] neighbor */
            }
            iat = nCumulene[len]; /* new current */
        }
    }
    return 0; /* failed */
}
/******************************************************************************************************/
int set_bond_type( inp_ATOM *at, AT_NUMB i1, AT_NUMB i2, int bType )
{
    AT_NUMB *p1 = is_in_the_list( at[i1].neighbor, i2, at[i1].valence );
    AT_NUMB *p2 = is_in_the_list( at[i2].neighbor, i1, at[i2].valence );
    if ( p1 && p2 ) {
        int j1 = p1 - at[i1].neighbor;
        int j2 = p2 - at[i2].neighbor;
        int bTypePrev = at[i1].bond_type[j1];
        at[i1].bond_type[j1] = bType;
        at[i2].bond_type[j2] = bType;
        if ( bTypePrev && bTypePrev <= BOND_TYPE_TRIPLE &&
             bType     && bType     <= BOND_TYPE_TRIPLE ) {
            at[i1].chem_bonds_valence += bType - bTypePrev;
            at[i2].chem_bonds_valence += bType - bTypePrev;
        }
        return 0;
    }
    return RI_ERR_SYNTAX;
}
/******************************************************************************************************/
int set_cumulene_0D_parity( inp_ATOM *at, inp_ATOM_STEREO *st, int num_at, int idelH1, int i1, int i2, int idelH2, int parity, int len )
{
    AT_NUMB nCumulene[MAX_CUMULENE_LEN+2];
    AT_NUMB *p1, *p2;
    int     m1, m2, parity1, parity2, sb_ord_m1, sb_ord_m2, k1, k2, num_neigh1, num_neigh2;
    /* the following types must exactly match types in inp_ATOM and inp_ATOM_STEREO */
    S_CHAR  *sb_ord1, *sn_ord1, *sb_parity1;
    S_CHAR  *sb_ord2, *sn_ord2, *sb_parity2;
    AT_NUMB *sn_orig_at_num1;
    AT_NUMB *sn_orig_at_num2;


    if ( !bFindCumuleneChain( at, (AT_NUMB)i1, (AT_NUMB)i2, nCumulene, len ) ) {
        return RI_ERR_SYNTAX; /* not an allene */
    }
    /* stereo bond neighbors: index of a stereo bond in its end-atom adjacency lists */
    if ( (p1 = is_in_the_list( at[i1].neighbor, nCumulene[1], at[i1].valence )) &&
         (p2 = is_in_the_list( at[i2].neighbor, nCumulene[len-1], at[i2].valence )) ) {
        sb_ord_m1 = p1 - at[i1].neighbor; /* indes of stereobond in the atom's adjacency list */
        sb_ord_m2 = p2 - at[i2].neighbor;
    } else {
        return RI_ERR_PROGR;
    }
    num_neigh1 = at[i1].valence + at[i1].num_H;
    num_neigh2 = at[i2].valence + at[i2].num_H;

    if ( num_neigh1 < MIN_NUM_STEREO_BOND_NEIGH || num_neigh1 > MAX_NUM_STEREO_BOND_NEIGH ||
         num_neigh2 < MIN_NUM_STEREO_BOND_NEIGH || num_neigh2 > MAX_NUM_STEREO_BOND_NEIGH ) {
        return RI_ERR_SYNTAX;
    }


    sb_ord1    = st? st[i1].sb_ord : at[i1].sb_ord;
    sb_ord2    = st? st[i2].sb_ord : at[i2].sb_ord;
    sb_parity1 = st? st[i1].sb_parity : at[i1].sb_parity;
    sb_parity2 = st? st[i2].sb_parity : at[i2].sb_parity;

    /* find the first unoccupied locations in the stereobond 0D descriptor lists; check whether the stereo has already been set */
    for( m1 = k1 = 0; m1 < MAX_NUM_STEREO_BONDS && sb_parity1[m1] && !(k1 = sb_ord1[m1] == sb_ord_m1); m1 ++ )
        ;
    for( m2 = k2 = 0; m2 < MAX_NUM_STEREO_BONDS && sb_parity2[m2] && !(k2 = sb_ord2[m2] == sb_ord_m2); m2 ++ )
        ;
    if ( m1 == MAX_NUM_STEREO_BONDS || m2 == MAX_NUM_STEREO_BONDS ) {
        return RI_ERR_SYNTAX;
    }
    if ( k1 && k2 ) {
        return 0; /* the stereo descriptor of this bond/allene/cumulene has already been set */
    }
    if ( k1 || k2 ) {
        return RI_ERR_SYNTAX; /* only half of a bond was set */
    }

    sn_ord1    = st? st[i1].sn_ord : at[i1].sn_ord;
    sn_ord2    = st? st[i2].sn_ord : at[i2].sn_ord;
    sn_orig_at_num1 = st? st[i1].sn_orig_at_num : at[i1].sn_orig_at_num;
    sn_orig_at_num2 = st? st[i2].sn_orig_at_num : at[i2].sn_orig_at_num;

    /* stereo bond neighbors connection index */
    sb_ord1[m1] = sb_ord_m1;
    sb_ord2[m2] = sb_ord_m2;
    /* stereo bond end atom neighbors */
    sn_orig_at_num1[m1] = at[idelH1].orig_at_number;
    if ( idelH1 < num_at ) {
        if ( (p1 = is_in_the_list( at[i1].neighbor, (AT_NUMB)idelH1, at[i1].valence )) ) {
            sn_ord1[m1] = p1 - at[i1].neighbor;
        } else {
            return RI_ERR_PROGR;
        }
    } else {
        sn_ord1[m1] = -1;
    }
    
    sn_orig_at_num2[m2] = at[idelH2].orig_at_number;
    if ( idelH2 < num_at ) {
        if ( (p2 = is_in_the_list( at[i2].neighbor, (AT_NUMB)idelH2, at[i2].valence )) ) {
            sn_ord2[m2] = p2 - at[i2].neighbor;
        } else {
            return RI_ERR_PROGR;
        }
    } else {
        sn_ord2[m2] = -1;
    }
    if ( ATOM_PARITY_WELL_DEF(parity) ) {
        /* special case: 2 bonds to sb atom => inverse parity because */
        /* InChI parity refers to the lone pair as a neighbor */
        int num_inv = (num_neigh1 == MIN_NUM_STEREO_BOND_NEIGH) + (num_neigh2 == MIN_NUM_STEREO_BOND_NEIGH);
        if ( num_inv % 2 ) {
            parity = (parity == AB_PARITY_EVEN)? AB_PARITY_ODD : AB_PARITY_EVEN;
        }
        parity1 = AB_PARITY_EVEN;
        parity2 = (parity == AB_PARITY_EVEN)? AB_PARITY_EVEN : AB_PARITY_ODD;
    } else {
        parity1 = parity2 = parity;
    }
    sb_parity1[m1] = parity1;
    sb_parity2[m2] = parity2;

    return 0;
}
/******************************************************************************************************/
int set_atom_0D_parity( inp_ATOM *at, inp_ATOM_STEREO *st, int num_at, int num_deleted_H, int i1, int parity )
{
    int     m1=0, m2, i, j, tot_num_neigh;
    /* the following types must exactly match types in inp_ATOM and inp_ATOM_STEREO */
    /* Given parity from InChI, the order of stereo center neighbors is: */
    /* 1. The stereocenter itself if the total number of neighbors is 3 (not 4) */
    /* 2. Explicit H: non-isotopic, isotopic in order ofascending  atomic mass */
    /*         Explicit H have already been sorted in this order */
    /* 3. Normal neighboring atoms, atom numbers (=canonical numbers from InChI - 1) in ascending order */
    /*         Normal neighboring atoms have already been sorted in this order */
    S_CHAR  *p_parity;
    AT_NUMB *p_orig_at_num;

    if ( !st || !at[i1].p_parity ) {
        m1            = 0;
        p_parity      = st? &st[i1].p_parity : &at[i1].p_parity;
        p_orig_at_num = st? st[i1].p_orig_at_num : at[i1].p_orig_at_num;

        tot_num_neigh = at[i1].valence + at[i1].num_H;
        if ( tot_num_neigh == MAX_NUM_STEREO_ATOM_NEIGH-1 ) {
            /* only 3 neighbors: the atom itself is the first neighbor */
            p_orig_at_num[m1 ++] = at[i1].orig_at_number;
        } else
        if ( tot_num_neigh != MAX_NUM_STEREO_ATOM_NEIGH ) {
            return RI_ERR_PROGR; /* wrong number of members */
        }
        m2 = m1 + (MAX_NUM_STEREO_ATOM_NEIGH - at[i1].valence);
        /* stereoneighbors: deleted explicit atoms H first, in order of increasing isotopic mass */
        if ( at[i1].num_H ) {
            for ( j = 0; m1 < m2 && j < num_deleted_H; j ++ ) {
                if ( at[j + num_at].neighbor[0] == i1 ) {
                    p_orig_at_num[m1 ++] = at[j + num_at].orig_at_number;
                }
            }
        }
        if ( m1 + at[i1].valence != MAX_NUM_STEREO_ATOM_NEIGH ) {
            return RI_ERR_PROGR; /* wrong number of members */
        }
        /* stereoneighbors: other than explicit H atoms */
        for ( i = 0; i < at[i1].valence; i ++ ) {
            m2 = at[i1].neighbor[i];
            p_orig_at_num[m1 ++] = at[m2].orig_at_number;
        }
        *p_parity = parity;
    }

    return 0;
}
#if ( BNS_RAD_SEARCH == 1 )
/******************************************************************************************************/
int MoveRadToAtomsAddCharges( BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct,
                    inp_ATOM *at, inp_ATOM *at2, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups, int forbidden_mask )
{
    int nNumRad, ret = 0, ret2;
    int i, j, k, num_rad_not_atom, num_moved=0, num_candidates= 0, extra_charge=0, added_charge, delta;
    BNS_EDGE   *pEdge;
    BNS_VERTEX *pv;
    Vertex      v1, v2;
    S_SHORT    *pnRad = NULL, *pnDelta = NULL;
    CC_CAND    *pCand = NULL;
    int         cnBits, bAtomRadRemoved = 0;

    int num_at = pStruct->num_atoms;
    int num_deleted_H = pStruct->num_deleted_H;
    int len_at = num_at + num_deleted_H;

    for ( i = pBNS->num_atoms, num_rad_not_atom=0; i < pBNS->num_vertices; i ++ ) {
        num_rad_not_atom += pBNS->vert[i].st_edge.cap - pBNS->vert[i].st_edge.flow;
    }
    if ( !num_rad_not_atom ) {
        goto exit_function;
    }
    /****************************************************/
    /*                                                  */
    /*    Move radicals from ChargeStruct to atoms      */
    /*                                                  */
    /****************************************************/
    
    /* allocate memory to keep track of moved radicals */
    pnRad   = (S_SHORT *) inchi_malloc(pBNS->num_vertices * sizeof(pnRad[0]));
    pnDelta = (S_SHORT *)inchi_calloc(pBNS->num_atoms, sizeof(pnDelta[0]));
    if ( !pnRad || !pnDelta ) {
        ret = RI_ERR_ALLOC;
        goto exit_function;
    }
    for ( i = 0; i < pBNS->num_vertices; i ++ ) {
        pnRad[i] = pBNS->vert[i].st_edge.cap - pBNS->vert[i].st_edge.flow;
    }
    while( 1 ) {
        /* remove radicals from atoms */
        for ( i = 0; i < pBNS->num_atoms; i ++ ) {
            pnDelta[i] = pBNS->vert[i].st_edge.cap - pBNS->vert[i].st_edge.flow;
            pBNS->vert[i].st_edge.cap -= pnDelta[i];
            bAtomRadRemoved += (0 != pnDelta[i]);
        }
        ret = SetRadEndpoints( pBNS, pBD, RAD_SRCH_FROM_FICT );
        if ( !ret ) {
            break;
        }
        if ( ret < 0 ) {
            goto exit_function;
        }
        nNumRad = ret;
        for ( i = 0; i < nNumRad; i ++ ) {
            pEdge = pBNS->edge + pBD->RadEdges[i];
            v1 = pEdge->neighbor1;
            v2 = pEdge->neighbor12 ^ v1;
            pBNS->vert[v1].st_edge.flow -=   pEdge->flow;
            pBNS->vert[v2].st_edge.flow -=   pEdge->flow;
            pBNS->tot_st_flow           -= 2*pEdge->flow;
            pEdge->flow                  = 0;
            pEdge->forbidden |= forbidden_mask;
            pBNS->edge_forbidden_mask |= forbidden_mask;
        }
        ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
        if ( ret < 0 ) {
            goto exit_function;
        } else {
            num_moved += ret;
        }
        RemoveRadEndpoints( pBNS, pBD, NULL );
        if ( ret == 0 ) {
            break;  /* could not move more radicals */
        }
        if ( bAtomRadRemoved ) {
            /* restore radicals to atoms */
            for ( i = 0; i < pBNS->num_atoms; i ++ ) {
                pBNS->vert[i].st_edge.cap += pnDelta[i];
            }
            bAtomRadRemoved = 0;
        }
    }
    if ( bAtomRadRemoved ) {
        /* restore radicals to atoms */
        for ( i = 0; i < pBNS->num_atoms; i ++ ) {
            pBNS->vert[i].st_edge.cap += pnDelta[i];
        }
        bAtomRadRemoved = 0;
    }
    pBNS->edge_forbidden_mask &= ~forbidden_mask;


    /****************************************************/
    /*                                                  */
    /*    Fix the charges                               */
    /*                                                  */
    /****************************************************/
    if ( num_moved ) {
        /* find reqired charge */
        extra_charge = 0;
        for ( i = pBNS->num_atoms, pv=pBNS->vert+i; i < pBNS->num_vertices; i ++, pv++ ) {
            if ( (delta = pv->st_edge.cap - pv->st_edge.flow) ) {
                if ( IS_BNS_VT_C_OR_CSUPER_GR(pv->type) ) {
                    extra_charge -= delta;
                } else
                if ( BNS_VERT_TYPE__AUX == pv->type ) {
                    extra_charge += delta;
                } else {
                    ret = RI_ERR_PROGR;
                    goto exit_function;
                }
            }
        }
        if ( !extra_charge ) {
            goto exit_function;
        }
        /* find differences */
        num_candidates = 0;
        for ( i = 0; i < pBNS->num_vertices; i ++ ) {
            pnRad[i] = (pBNS->vert[i].st_edge.cap - pBNS->vert[i].st_edge.flow) - pnRad[i];
            if ( i < pBNS->num_atoms && pnRad[i] > 0 && !pVA[i].nTautGroupEdge ) {
                num_candidates ++;
            }
        }
    }
    if ( num_candidates > 0 ) {
        pCand = (CC_CAND *)inchi_calloc( num_candidates, sizeof(pCand[0]) );
        if ( !pCand ) {
            ret = RI_ERR_ALLOC;
            goto exit_function;
        }
        /* create atom */
        memcpy( at2, at, len_at*sizeof(at2[0]));
        pStruct->at = at2;
        ret2 = CopyBnsToAtom( pStruct, pBNS, pVA, pTCGroups, 1 );
        pStruct->at = at;
        if ( ret2 < 0 ) {
            ret = ret2;
            goto exit_function;
        }

        for ( i = 0, j = 0; i < pBNS->num_vertices; i ++ ) {
            if ( i < pBNS->num_atoms && pnRad[i] > 0 && !pVA[i].nTautGroupEdge ) {
                pCand[j].iat                  = i;
                pCand[j].num_bonds            = at2[i].valence;
                pCand[j].chem_valence         = at2[i].chem_bonds_valence;
                pCand[j].cMetal               = pVA[i].cMetal;
                pCand[j].cNumBondsToMetal     = pVA[i].cNumBondsToMetal;
                pCand[j].cNumValenceElectrons = pVA[i].cNumValenceElectrons;
                pCand[j].cPeriodicRowNumber   = pVA[i].cPeriodicRowNumber;
                pCand[j].el_number            = at2[i].el_number;
                cnBits = (pVA[i].cnListIndex > 0)? cnList[pVA[i].cnListIndex-1].bits : 0;
                while ( cnBits > 0 ) {
                    pCand[j].cNumChargeStates ++;
                    cnBits >>= cn_bits_shift;
                }
                j ++;
            }
        }
        if ( j > 1 ) {
            qsort( pCand, j, sizeof(pCand[0]), comp_cc_cand );
        }
        added_charge = 0;
        
        for ( k = 0; k < j; k ++ ) {
            int rest_of_charge = extra_charge - added_charge;
            int charge_per_left_atom = (abs(rest_of_charge) + j-k - 1)/(j-k);
            int this_atom_add_charge = rest_of_charge > 0? charge_per_left_atom : -charge_per_left_atom;
            pVA[pCand[k].iat].cInitCharge += this_atom_add_charge;
            added_charge                  += this_atom_add_charge;
            if ( this_atom_add_charge ) {
                for ( i = pBNS->num_vertices-1, pv = pBNS->vert + i; this_atom_add_charge && pBNS->num_atoms <= i; i --, pv -- ) {
                    if ( (delta = pv->st_edge.cap - pv->st_edge.flow) ) {
                        if ( this_atom_add_charge < 0 && IS_BNS_VT_C_OR_CSUPER_GR(pv->type) ) {
                            if ( delta + this_atom_add_charge > 0 ) {
                                delta = -this_atom_add_charge;
                            }
                            pv->st_edge.cap      -= delta;
                            pBNS->tot_st_cap     -= delta;
                            this_atom_add_charge += delta;
                        } else
                        if ( this_atom_add_charge > 0 && BNS_VERT_TYPE__AUX == pv->type ) {
                            if ( delta > this_atom_add_charge ) {
                                delta = this_atom_add_charge;
                            }
                            pv->st_edge.cap      -= delta;
                            pBNS->tot_st_cap     -= delta;
                            this_atom_add_charge -= delta;
                        }
                    }  
                }
            }
        }
    }

exit_function:
    if ( pnRad ) {
        inchi_free( pnRad );
    }
    if ( pnDelta ) {
        inchi_free( pnDelta );
    }
    if ( pCand ) {
        inchi_free( pCand );
    }
    return ret;
}
#endif
/**************************************************************************************************/
typedef struct tagMobileHGroups {
    AT_NUMB group_number;
    AT_NUMB atom_number;
    AT_NUMB atom_type_pVA;
    S_CHAR  ineigh;
    S_CHAR  bond_type;
    S_CHAR  forbidden;
   /* S_CHAR  el_type;*/
    S_CHAR  endpoint_valence;
    S_CHAR  num_bonds;
    S_CHAR  bonds_valence;
    S_CHAR  num_bonds_non_metal;
    S_CHAR  bonds_valence_non_metal;

} MOBILE_GR;

typedef struct tagMobileGroupList {
    AT_NUMB group_number;
    AT_NUMB num;
} MGROUPS;
/**************************************************************************************************/
int AdjustTgroupsToForbiddenEdges2( BN_STRUCT *pBNS, inp_ATOM *at, VAL_AT *pVA,
                                    int num_atoms, int forbidden_mask )
{
    int i, j, k;
    int centerpoint_type, neigh_type;
    int num_changes;
    int num_donors, num_acceptors, num_donor_endpoints, num_acceptor_endpoints;
    int neigh, tg_number, num_eql_mobile_gr, num_dif_mobile_gr, bond_type, has_mobile_H, has_mobile;
    int num_forbidden, ind_forbidden, forbidden, num_N, num_O, num_P, num_S, num_OSt;
    int val, delta_val, delta_met, num_bonds_non_metal, bonds_valence_non_metal;
    int num_bonds, bonds_valence;
    int inv_forbidden_mask = ~forbidden_mask;
    MOBILE_GR MobileGr[MAXVAL];
    int       num_endpoints;
    MGROUPS   MGroups[MAXVAL];
    int       num_mgroups, num_diff_t_groups;
    BNS_EDGE   *e, *e1, *e2, *ev, *ev1, *ev2;
    BNS_VERTEX *pv1, *pv2;
    num_changes = 0;
    /* search for possible centerpoints */
    for ( i = 0; i < num_atoms; i ++ ) {
        
        if ( at[i].chem_bonds_valence == at[i].valence || at[i].num_H ||
             at[i].endpoint || at[i].charge || at[i].radical ||
             !is_centerpoint_elem(at[i].el_number) ||
             !(centerpoint_type = get_pVA_atom_type( pVA, at, i, 0 )) ||
             2 > (delta_val = at[i].chem_bonds_valence - (val = get_el_valence(at[i].el_number, 0, 0))) ||
             2 > (delta_met = (bonds_valence_non_metal = nNoMetalBondsValence(at, i)) - val )
           ) {
            continue;
        }
        
        num_donors = num_acceptors = num_donor_endpoints = num_acceptor_endpoints = 0;
        num_mgroups = num_endpoints = num_diff_t_groups = 0;
        has_mobile = has_mobile_H = num_eql_mobile_gr = num_dif_mobile_gr = tg_number = 0;
        ind_forbidden = -1;
        num_forbidden = 0;
        num_N = num_O = num_P = num_S = num_OSt = 0;
        num_bonds_non_metal = nNoMetalNumBonds(at, i);
        bonds_valence = at[i].chem_bonds_valence;
        num_bonds     = at[i].valence;

        for ( j = 0; j < at[i].valence; j ++ ) {
            /* collect neighbors info */
            neigh      = at[i].neighbor[j];
            val        = get_endpoint_valence( at[neigh].el_number );
            forbidden  = pBNS->edge[(int)pBNS->vert[i].iedge[j]].forbidden;
            bond_type = (at[i].bond_type[j] & BOND_TYPE_MASK);
            neigh_type = get_pVA_atom_type( pVA, at, neigh, bond_type);
            if ( !forbidden && !at[neigh].endpoint ) {
                /* save forbidden bonds */
                if ( is_el_a_metal(at[neigh].el_number) ) {
                    continue;
                }
                switch( bond_type ) {
                case BOND_TYPE_SINGLE:
                    if ( !at[neigh].num_H && at[neigh].charge != -1 ) {
                        continue; /* not a donor */
                    }
                    break;
                case BOND_TYPE_DOUBLE:
                    if ( !neigh_type ) {
                        continue;
                    }
                    break;
                default:
                    continue;
                }
            }

            MobileGr[num_endpoints].atom_number      = neigh;
            MobileGr[num_endpoints].ineigh           = j;
            MobileGr[num_endpoints].bond_type        = bond_type;
            MobileGr[num_endpoints].group_number     = at[neigh].endpoint;
            MobileGr[num_endpoints].endpoint_valence = val;
            MobileGr[num_endpoints].forbidden        = forbidden;
            MobileGr[num_endpoints].atom_type_pVA    = neigh_type;
            MobileGr[num_endpoints].num_bonds        = at[neigh].valence;
            MobileGr[num_endpoints].bonds_valence    = at[neigh].chem_bonds_valence;
            MobileGr[num_endpoints].num_bonds_non_metal     = nNoMetalNumBonds(at, neigh);
            MobileGr[num_endpoints].bonds_valence_non_metal = nNoMetalBondsValence( at, neigh );

            if ( forbidden & forbidden_mask ) {
                num_forbidden ++;
                ind_forbidden = num_endpoints;
            }
            num_O   += 0 != (neigh_type & EL_TYPE_O) && at[neigh].valence == 1; /* ignore -O- */
            num_N   += 0 != (neigh_type & EL_TYPE_N) &&
                            !(at[neigh].valence == 3 && at[neigh].chem_bonds_valence == 3); /* ignore -N< */
            num_S   += 0 != (neigh_type & EL_TYPE_S) && at[neigh].valence == 1; /* ignore -S- */
            num_P   += 0 != (neigh_type & EL_TYPE_P) &&
                            !(at[neigh].valence == 3 && at[neigh].chem_bonds_valence == 3); /* ignore -P< */
            num_OSt += 0 != (neigh_type & EL_TYPE_OSt);
            num_acceptors += (bond_type == BOND_TYPE_DOUBLE) && (neigh_type & EL_TYPE_PT);
            num_donors    += (bond_type == BOND_TYPE_SINGLE) && (neigh_type & EL_TYPE_PT) &&
                             (at[neigh].num_H || at[neigh].charge==-1 || at[neigh].endpoint);
            if ( at[neigh].endpoint ) {
                num_acceptor_endpoints += (bond_type == BOND_TYPE_DOUBLE);
                num_donor_endpoints    += (bond_type == BOND_TYPE_SINGLE);
                if ( !tg_number ) {
                    tg_number = at[neigh].endpoint;
                    num_eql_mobile_gr = 1;
                } else
                if ( tg_number == at[neigh].endpoint ) {
                    num_eql_mobile_gr ++;
                } else {
                    num_dif_mobile_gr ++;
                }
            } else
            if ( bond_type == BOND_TYPE_SINGLE && val ) {
                if ( at[neigh].endpoint ) {
                    has_mobile_H |= 1; 
                    has_mobile   |= 1; 
                } else {
                    has_mobile_H |= (0 != at[neigh].num_H); 
                    has_mobile   |= (0 != at[neigh].num_H) || (at[neigh].charge == -1); 
                }
            }
            num_endpoints ++;

            if ( at[neigh].endpoint || (neigh_type & EL_TYPE_PT) ) {
                for ( k = 0; k < num_mgroups; k ++ ) {
                    if ( MGroups[k].group_number == at[neigh].endpoint ) {
                        MGroups[k].num ++;
                        break;
                    }
                }
                if ( k == num_mgroups ) {
                    MGroups[k].group_number = at[neigh].endpoint;
                    MGroups[k].num          = 1;
                    num_mgroups ++;
                    num_diff_t_groups += (0 != at[neigh].endpoint);
                }
            }
        }
        if ( !num_acceptors || !num_donors || /* num_acceptors > 2 ||*/
             (num_eql_mobile_gr == num_endpoints && !num_forbidden) ||
             (!tg_number && !has_mobile_H) ) {
            continue; /* nothing to do */
        }
        
/* case_5_1: */
        /***************** determine the case ************************/
        if ( 3 == num_bonds_non_metal &&
             4 == bonds_valence_non_metal &&
             (centerpoint_type == EL_TYPE_C) &&
             2 == num_O && 1 == num_N+num_S && num_OSt &&
             1 == num_forbidden && 3 == num_eql_mobile_gr  ) {
            /******************************************************** 
             ***         InChI Tech. Man., Table 5, case 1        *** 
             ******************************************************** 
                       2                  
                      OH                OH     X  =  N, S, Se, Te
                     /                 /       f  =  fixed bond
                 e  /                 /        tg =  Mobile-H vertex
              HX---C      -->    X===C   
           ev2|| f  \\ ev1         f  \  
              ||     \\                \ 
              tg------O                 OH
                       1
            Problem:
              XH, O, and O belong to the same Mobile-H group.
              Fixed bond prevents the correct structure restoration:
              H cannot migrate from X to O because HX-N bond is fixed
            Solution:
              Move H from X to allow XH-C bond change
              (this unfixes the bond, see SetForbiddenEdges(...) )
            *********************************************************/
            int jXH = -1, jO1 = -1, jO2 = -1, n = 0;
            for ( j = 0; j < num_endpoints; j ++ ) {
                if ( (MobileGr[j].atom_type_pVA & (EL_TYPE_N | EL_TYPE_S)) &&
                     (MobileGr[j].forbidden == forbidden_mask) &&
                     MobileGr[j].bond_type == BOND_TYPE_SINGLE &&
                     jXH < 0 ) {
                    jXH = j;
                    n ++;
                } else
                if ( (MobileGr[j].atom_type_pVA & EL_TYPE_MASK) == EL_TYPE_O &&
                     MobileGr[j].num_bonds_non_metal == 1 &&
                     !MobileGr[j].forbidden ) {
                    if ( MobileGr[j].bond_type == BOND_TYPE_DOUBLE && jO1 < 0 ) {
                        jO1 = j;
                        n ++;
                    } else
                    if ( MobileGr[j].bond_type == BOND_TYPE_SINGLE && jO2 < 0 ) {
                        jO2 = j;
                        n ++;
                    }
                }
            }
            if ( n != 3 ) {
                goto case_5_2;
            }
            /* XH-C edge */
            e   = pBNS->edge + pBNS->vert[i].iedge[MobileGr[jXH].ineigh];
            /* C=O  edge */
            ev1 = pBNS->edge + pBNS->vert[i].iedge[MobileGr[jO1].ineigh];
            /* XH-tg edge */
            ev2 = pBNS->edge + pVA[MobileGr[jXH].atom_number].nTautGroupEdge - 1;

            if ( !ev1->flow || !ev2->flow ) {
                goto case_5_2;
            }

            /* do not remove forbidden edge bit */
            e->flow ++;
            ev1->flow --;
            ev2->flow --;
            pBNS->vert[ev1->neighbor12 ^ i].st_edge.flow --;
            pBNS->vert[ev2->neighbor12 ^ ev2->neighbor1].st_edge.flow --;
            pBNS->tot_st_flow -= 2;
            num_changes ++;
            continue;
        }
case_5_2:
        /*********************************************************************/
        if ( 3 == num_bonds_non_metal &&
             5 == bonds_valence_non_metal &&
             (centerpoint_type == EL_TYPE_N) &&
             2 == num_O && 1 == num_N+num_S &&
             1 == num_forbidden && 3 == num_eql_mobile_gr  ) {
            /******************************************************** 
             ***         InChI Tech. Man., Table 5, case 2        *** 
             ******************************************************** 
                                         
                      O                 OH     X  =  N, S, Se, Te
                     //                /       f  =  fixed bond
                 e  //                /        tg =  Mobile-H vertex
              HX---N      -->    X===N   
           ev2|| f  \\ ev1         f  \\ 
              ||     \\                \\
              tg------O                 O
              
            Problem:
              XH, O, and O belong to the same Mobile-H group.
              Fixed bond prevents the correct structure restoration:
              H cannot migrate from X to O because HX-N bond is fixed
            Solution:
              Move H from X to allow XH-N bond change
              (this unfixes the bond, see SetForbiddenEdges(...) )
            *********************************************************/
            int jXH = -1, jO1 = -1, jO2 = -1, n = 0;
            for ( j = 0; j < num_endpoints; j ++ ) {
                if ( (MobileGr[j].atom_type_pVA & (EL_TYPE_N | EL_TYPE_S)) &&
                     (MobileGr[j].forbidden == forbidden_mask) &&
                     MobileGr[j].bond_type == BOND_TYPE_SINGLE &&
                     jXH < 0 ) {
                    jXH = j;
                    n ++;
                } else
                if ( (MobileGr[j].atom_type_pVA & EL_TYPE_MASK) == EL_TYPE_O &&
                     MobileGr[j].bond_type == BOND_TYPE_DOUBLE &&
                     MobileGr[j].num_bonds_non_metal == 1 &&
                     !MobileGr[j].forbidden ) {
                    if ( jO1 < 0 ) {
                        jO1 = j;
                        n ++;
                    } else
                    if ( jO2 < 0 ) {
                        jO2 = j;
                        n ++;
                    }
                }
            }
            if ( n != 3 ) {
                goto case_5_4;
            }
            /* XH-N edge */
            e   = pBNS->edge + pBNS->vert[i].iedge[MobileGr[jXH].ineigh];
            /* N=O  edge */
            ev1 = pBNS->edge + pBNS->vert[i].iedge[MobileGr[jO1].ineigh];
            /* XH-tg edge */
            ev2 = pBNS->edge + pVA[MobileGr[jXH].atom_number].nTautGroupEdge - 1;

            if ( !ev1->flow || !ev2->flow ) {
                goto case_5_4;
            }
            /* do not remove forbidden edge bit */
            e->flow ++;
            ev1->flow --;
            ev2->flow --;
            pBNS->vert[ev1->neighbor12 ^ i].st_edge.flow --; /* first =O vertex */
            pBNS->vert[ev2->neighbor12 ^ ev2->neighbor1].st_edge.flow --; /* taut group vertex tg */
            pBNS->tot_st_flow -= 2;
            num_changes ++;
            continue;
        }
case_5_4:
        /*********************************************************************/
        if ( 3 == num_bonds_non_metal &&
             5 == bonds_valence_non_metal &&
             (centerpoint_type & (EL_TYPE_N | EL_TYPE_P)) &&
             1 == num_O+num_S && 0 < num_N && 2 == (num_N + num_P) &&
             1 == num_forbidden && num_O+num_S+num_N == num_eql_mobile_gr  ) {
            /******************************************************** 
             ***         InChI Tech. Man., Table 5, case 4        *** 
             ******************************************************** 
                                               O  =  O, S, Se, Te
                       X                 X     X  =  N, P, As
                     //                //      f  =  fixed bond
                 e  // ev2            //       tg =  Mobile-H vertex
               O===N      -->   HO---N   
              || f  \  ev1         f  \\ 
              ||     \                 \\
              tg------NH                N 
                        
            Problem:
              O, NH, and possibly X belong to the same Mobile-H group.
              Fixed bond prevents the correct structure restoration:
              H cannot migrate from NH to O because O=N bond is fixed
            Solution:
              Move H from NH to O to allow O=N bond change
              (this unfixes the bond, see fix_special_bonds(...) )
            *********************************************************/
            int jO = -1, jNH = -1, jX = -1, n = 0;
            for ( j = 0; j < num_endpoints; j ++ ) {
                if ( (MobileGr[j].atom_type_pVA & (EL_TYPE_O | EL_TYPE_S)) &&
                     MobileGr[j].forbidden == forbidden_mask &&
                     MobileGr[j].bond_type == BOND_TYPE_DOUBLE &&
                     MobileGr[j].num_bonds_non_metal == 1 &&
                     jO < 0 ) {
                    jO = j;
                    n ++;
                } else
                if ( (MobileGr[j].atom_type_pVA & (EL_TYPE_N | EL_TYPE_P)) &&
                     !MobileGr[j].forbidden ) {
                    if ( MobileGr[j].bond_type == BOND_TYPE_SINGLE &&
                         (MobileGr[j].atom_type_pVA & EL_TYPE_MASK) == EL_TYPE_N && jNH < 0 ) {
                        jNH = j;
                        n ++;
                    } else
                    if ( MobileGr[j].bond_type == BOND_TYPE_DOUBLE && jX < 0 ) {
                        jX = j;
                        n ++;
                    }
                }
            }
            if ( n != 3 ) {
                goto case_5_6;
            }
            /* O=N edge */
            e   = pBNS->edge + pBNS->vert[i].iedge[MobileGr[jO].ineigh];
            /* N-NH  edge */
            ev1 = pBNS->edge + pBNS->vert[i].iedge[MobileGr[jNH].ineigh];
            /* N=X edge */
            ev2 = pBNS->edge + pBNS->vert[i].iedge[MobileGr[jX].ineigh];

            if ( !e->flow ) {
                goto case_5_6;
            }
            /* do not remove forbidden edge bit */
            e->flow --;
            ev1->flow ++;
            pBNS->vert[e->neighbor12 ^ i].st_edge.flow --;
            pBNS->vert[ev1->neighbor12 ^ i].st_edge.flow --;
            pBNS->tot_st_flow -= 2;
            num_changes ++;
            continue;
        }
case_5_6:
        /********* InChI Tech.Man. Table 5, case 6 **************/
        if ( 2 == delta_met && 4 == num_bonds_non_metal &&
             5 == bonds_valence_non_metal &&
             1 == num_forbidden && 1 < num_eql_mobile_gr &&
             !num_dif_mobile_gr &&
             (centerpoint_type & (EL_TYPE_N | EL_TYPE_P)) &&
             1 <= num_N && 2 <= num_N+num_O+num_S &&
             1 == num_acceptor_endpoints && 0 < num_donor_endpoints ) {
            int jN = -1, njFix = 0, jFix[4], n = 0;
                /* centerpoint is N, P, As, Sb

                  input             output
                  -----             ------
                     end 
                     po- 
                     int 
                     2
                     
                  X  ZH              X  Z     Z=N,O,S,Se,Te [terminal endpoint]
                   \ |                \ ||      
                    \| f   f           \||      
                 Y---N===N---       Y---N---NH---
                       e                  f      
                    cen  end         no bond           
                    ter  po-         fixed  
                    po-  int          
                    int  1        tautomerism O==N--NH is allowed
                    
                  Problem: OH and =N- belong to a Mobile-H group, but
                           forbidden edge e does not allow them to be
                           tautomeric in the restored structure.

                  Solution:
                    
                  1. Decrement flow in edge e
                  2. Decrement st_edge flow in N and N connected by e
                  3. Fix all single order bonds to not terminal tautomeric N around N(centerpoint)
                  4. Run BNS to establist new flow distribution
                */
            /* fixed bond */
            for ( j = 0; j < num_endpoints; j ++ ) {
                neigh = MobileGr[j].atom_number;
                if ( MobileGr[j].bond_type == BOND_TYPE_DOUBLE &&
                     (MobileGr[j].atom_type_pVA & EL_TYPE_MASK) == EL_TYPE_N &&
                     MobileGr[j].num_bonds_non_metal == 2     &&
                     MobileGr[j].bonds_valence_non_metal == 3 &&
                     at[neigh].endpoint &&
                     !at[neigh].num_H && !at[neigh].charge &&
                     !at[neigh].radical &&
                     (MobileGr[j].forbidden & forbidden_mask) && jN < 0 ) {
                    jN = j;
                    n ++;
                } else
                if ( MobileGr[j].bond_type == BOND_TYPE_SINGLE &&
                     at[neigh].endpoint ) {
                    if ( MobileGr[j].num_bonds > 1 ) {
                        jFix[njFix ++] = j;
                    }
                    n ++;
                }
            }

            if ( jN < 0 || n < 2 || 1 + njFix == n ) {
                goto case_5_7;  /* nothing to do */
            }

            e   = pBNS->edge + pBNS->vert[i].iedge[MobileGr[jN].ineigh];  /* fixed edge */
            if ( !e->flow ) {
                goto case_5_7;
            }
            e->flow --;
            pBNS->vert[i].st_edge.flow --;
            pBNS->vert[e->neighbor12 ^ i].st_edge.flow --;
            pBNS->tot_st_flow -= 2;

            for ( j = 0; j < njFix; j ++ ) {
                /* edges to fix */
                ev = pBNS->edge + pBNS->vert[i].iedge[(int)MobileGr[jFix[j]].ineigh];
                ev->forbidden |= forbidden_mask;
            }
            num_changes ++;
            continue;
        }
case_5_7:
        /*********************************************************************/
        if ( 3 == num_bonds_non_metal &&
             4 == bonds_valence_non_metal &&
             (centerpoint_type == EL_TYPE_S) &&
             2 == num_O+num_S && 1 == num_OSt && 1 == num_N &&
             1 == num_forbidden && 3 == num_eql_mobile_gr && 
             MobileGr[ind_forbidden].bond_type == BOND_TYPE_SINGLE ) {
            /******************************************************** 
             ***         InChI Tech. Man., Table 5, case 7        *** 
             ******************************************************** 
                                               O  =  O, S, Se, Te
                      OH                OH     S  =  S, Se, Te
                     /                 /       f  =  fixed bond
                 e  /ev2           f  /        tg =  Mobile-H vertex
              HN---S      -->    N===S         X  =  N or non-endpoint;
              || f  \\                \
           ev2||     \\ev1             \ 
              tg------O                 OH
                                               N, O, O
            Problem:                           =======
              O, NH, OH belong to the same Mobile-H group.
              Fixed bond prevents the correct structure restoration:
              H cannot migrate from NH to O because HN-S bond is fixed
            Solution:
              Move H from NH to =O to allow HN=S bond change by making a 2nd terminal -OH
              (this unfixes the bond, see fix_special_bonds(...) )
            *********************************************************/
            int jO = -1, jNH = -1, jOH = -1, n = 0;
            for ( j = 0; j < num_endpoints; j ++ ) {
                if ( (MobileGr[j].atom_type_pVA & EL_TYPE_MASK) == EL_TYPE_N &&
                     MobileGr[j].forbidden == forbidden_mask &&
                     MobileGr[j].bond_type == BOND_TYPE_SINGLE &&
                     MobileGr[j].num_bonds_non_metal <= 2 &&
                     jNH < 0 ) {
                    jNH = j;
                    n ++;
                } else
                if ( (MobileGr[j].atom_type_pVA & (EL_TYPE_O | EL_TYPE_S)) &&
                     !MobileGr[j].forbidden &&
                     MobileGr[j].num_bonds_non_metal == 1 ) {
                    if ( MobileGr[j].bond_type == BOND_TYPE_DOUBLE &&
                         jO < 0 ) {
                        jO = j;
                        n ++;
                    } else
                    if ( jOH < 0 ) {
                        jOH = j;
                        n ++;
                    }
                }
            }
            if ( n != 3 ) {
                goto case_5_9a;
            }
            /* NH-S edge */
            e   = pBNS->edge + pBNS->vert[i].iedge[MobileGr[jNH].ineigh];
            /* S=O  edge */
            ev1 = pBNS->edge + pBNS->vert[i].iedge[MobileGr[jO].ineigh];
            /* XH-tg edge */
            ev2 = pBNS->edge + pVA[MobileGr[jNH].atom_number].nTautGroupEdge - 1;

            if ( !ev1->flow || !ev2->flow ) {
                goto case_5_9a;
            }

            /* do not remove forbidden edge bit */
            e->flow ++;
            ev1->flow --;
            ev2->flow --;
            pBNS->vert[ev1->neighbor12 ^ i].st_edge.flow --; /* first =O vertex */
            pBNS->vert[ev2->neighbor12 ^ ev2->neighbor1].st_edge.flow --; /* taut group vertex tg */
            pBNS->tot_st_flow -= 2;
            num_changes ++;
            continue;
        }
case_5_9a:
        /*********************************************************************/
        if ( 3 == num_bonds_non_metal &&
             4 == bonds_valence_non_metal &&
             (centerpoint_type == EL_TYPE_S) &&
             1 == num_O+num_S && !num_OSt && 1 <= num_N &&
             1 == num_forbidden &&
             num_O+num_S+num_N == num_eql_mobile_gr && 
             MobileGr[ind_forbidden].bond_type == BOND_TYPE_SINGLE ) {
            /******************************************************** 
             ***         InChI Tech. Man., Table 5, case 9a       *** 
             ******************************************************** 
                                               O  =  O, S, Se, Te
                      X                 X      S  =  S, Se, Te
                     /                 /       f  =  fixed bond
                    /                 /        tg =  Mobile-H vertex
              HN---S      -->    N===S         X  =  N or non-endpoint;
              ||    \\ e              \       -X is not -O(terminal)
              ||    f\\               f\ 
              tg------O                 OH
                                               N, N, O  or N, O
            Problem:                           ================
              O, NH belong to the same Mobile-H group.
              Fixed bond prevents the correct structure restoration:
              H cannot migrate from NH to O because O=S bond is fixed
            Solution:
              Move H from NH to =O to allow O=S bond change by making a terminal -OH
              (this unfixes the bond, see fix_special_bonds(...) )
            *********************************************************/
            int jO = -1, jNH = -1, jX = -1, n = 0;
            for ( j = 0; j < num_endpoints; j ++ ) {
                if ( (MobileGr[j].atom_type_pVA & (EL_TYPE_O | EL_TYPE_S)) &&
                     MobileGr[j].forbidden == forbidden_mask &&
                     MobileGr[j].bond_type == BOND_TYPE_DOUBLE &&
                     MobileGr[j].num_bonds_non_metal == 1 &&
                     jO < 0 ) {
                    jO = j;
                    n ++;
                } else
                if ( (MobileGr[j].atom_type_pVA & EL_TYPE_MASK) == EL_TYPE_N &&
                     !MobileGr[j].forbidden &&
                     MobileGr[j].bond_type == BOND_TYPE_SINGLE &&
                     jNH < 0 ) {
                    jNH = j;
                    n ++;
                } else
                if ( jX < 0 ) {
                    jX = j;
                    n ++;
                }
            }
            if ( jO < 0 || jNH < 0 ) {
                goto case_5_8b_to_9b;
            }

            e = pBNS->edge + pBNS->vert[i].iedge[MobileGr[ind_forbidden].ineigh];
            if ( !e->flow ) {
                goto case_5_8b_to_9b;
            }
            e->flow --;
            pBNS->vert[e->neighbor1].st_edge.flow --;
            pBNS->vert[e->neighbor1 ^ e->neighbor12].st_edge.flow --;
            pBNS->tot_st_flow -= 2;
            num_changes ++;
            continue;
        }
case_5_8b_to_9b: /* #1 */
        /*********************************************************************/
        if ( 3 == num_bonds_non_metal &&
             4 == bonds_valence_non_metal &&
             (centerpoint_type == EL_TYPE_S) &&
             0 == num_O+num_S && 2 == num_N && 0 == num_P && !num_OSt &&
             1 == num_forbidden &&
             1 == num_eql_mobile_gr && 0 == num_dif_mobile_gr &&
             1 == num_donor_endpoints && 0 == num_acceptor_endpoints &&
             1 == num_donors          && 1 == num_acceptors &&
             MobileGr[ind_forbidden].bond_type == BOND_TYPE_SINGLE ) {
            /******************************************************** 
             ***         InChI Tech. Man., Table 5, case 8b->9b   *** 
             ******************************************************** 
                             --->                    O  =  O, S, Se, Te      
                      X                  X           S  =  S, Se, Te         
             \      f/          \      f/            f  =  fixed bond        
              \ ev  /  C=====Z   \     /  C-----ZH   tg =  Mobile-H vertex   
               N===S   |     |    N===S  ||     ||   X  =  is N not an endpoint;
             not    \  |     |         \ ||     ||  -X is not terminal -O,-S,-Se,-Te or
             an      \ |  e  |          \||  e  ||        any N, P, As
           endpoint   NH=====tg          N------tg
                     is an                   f       N, N, X, fixed single
                   endpoint                          =====================

            Problem:
              N is not a Mobile-H endpoint, NH is a Mobile-H endpoint
              Unfixed bond N==S prevents the correct structure restoration:
              H can migrate from NH to N because N=S bond is not fixed
            Solution:
              Move H from NH to =Z to make N=S bond fixed (Table 5, case 9)
              (this unfixes the bond, see fix_special_bonds(...) )
            *********************************************************/
            int jN = -1, jNH = -1, jX = -1, n = 0;
            for ( j = 0; j < num_endpoints; j ++ ) {
                if ( (MobileGr[j].atom_type_pVA & EL_TYPE_MASK) == EL_TYPE_N &&
                     !(MobileGr[j].forbidden & forbidden_mask) ) {
                    if ( MobileGr[j].bond_type == BOND_TYPE_DOUBLE &&
                         !at[MobileGr[j].atom_number].endpoint &&
                         jN < 0 ) {
                        jN = j;
                        n ++;
                    } else
                    if ( MobileGr[j].bond_type == BOND_TYPE_SINGLE &&
                         MobileGr[j].num_bonds == 2 &&
                         MobileGr[j].bonds_valence == 2 &&
                         at[MobileGr[j].atom_number].endpoint &&
                         jNH < 0 ) {
                        jNH = j;
                        n ++;
                    }
                } else
                if ( !((MobileGr[j].atom_type_pVA & (EL_TYPE_N | EL_TYPE_P)) ||
                       ((MobileGr[j].atom_type_pVA & (EL_TYPE_O | EL_TYPE_S)) &&
                        MobileGr[j].num_bonds > 1) ) &&
                     (MobileGr[j].forbidden & forbidden_mask) &&
                     MobileGr[j].bond_type == BOND_TYPE_SINGLE &&
                     jX < 0 ) {
                    jX = j;
                    n ++;
                }
            }
            if ( n != 3 ) {
                goto case_5_8c_to_9c;
            }

            e = pBNS->edge + pVA[MobileGr[jNH].atom_number].nTautGroupEdge - 1;
            if ( !e->flow ) {
                goto case_5_8c_to_9c; /* should not happen ??? */
            }
            e->flow --;
            pBNS->vert[e->neighbor1].st_edge.flow --;
            pBNS->vert[e->neighbor1 ^ e->neighbor12].st_edge.flow --;
            pBNS->tot_st_flow -= 2;
            e->forbidden |= forbidden_mask;
            num_changes ++;
            continue;
        }
case_5_8c_to_9c: /* #2 */
        /*********************************************************************/
        if ( 3 == num_bonds_non_metal &&
             4 == bonds_valence_non_metal &&
             (centerpoint_type == EL_TYPE_S) &&
             0 == num_O+num_S && 3 == num_N && 0 == num_P &&
             1 == num_forbidden &&
             3 == num_eql_mobile_gr && 0 == num_dif_mobile_gr &&
             2 == num_donor_endpoints && 1 == num_acceptor_endpoints &&
             2 == num_donors          && 1 == num_acceptors &&
             MobileGr[ind_forbidden].bond_type == BOND_TYPE_SINGLE ) {
            /******************************************************** 
             ***         InChI Tech. Man., Table 5, case 8c->9c   *** 
             ******************************************************** 
                is an endpoint --->                  O  =  O, S, Se, Te      
                      NH2(X)             NH2(X)      S  =  S, Se, Te         
             \       / pv1  pv2 \       /            f  =  fixed bond        
              \ 1   /  C-----ZH  \     /  C=====Z    tg =  Mobile-H vertex   
               N===S  ||     ||   N===S   |     |    X  =  is N not an endpoint;
            is an   \f||ev1  ||ev2     \f |     |   -X is not terminal -O,-S,-Se,-Te or
          endpoint   \||  e  ||         \ |  e  |         any N, P, As
                    2 N------tg          NH=====tg   C is a centerpoint of a t-group
                     is an                   f       N, N, X, fixed single
                   endpoint                          =====================

            Problem:
              N is not a Mobile-H endpoint, NH is a Mobile-H endpoint
              Unfixed bond N==S prevents the correct structure restoration:
              H can migrate from NH to N because N=S bond is not fixed
            Solution:
              Move H from NH to =Z to make N=S bond fixed (Table 5, case 9)
              (this unfixes the bond, see fix_special_bonds(...) )
            *********************************************************/
            int jN1 = -1, jN2 = -1, jX = -1, n = 0;
            EdgeIndex ie, ie1, ie2;
            for ( j = 0; j < num_endpoints; j ++ ) {
                if ( (MobileGr[j].atom_type_pVA & EL_TYPE_MASK) == EL_TYPE_N &&
                     !(MobileGr[j].forbidden & forbidden_mask) ) {
                    if ( MobileGr[j].bond_type == BOND_TYPE_DOUBLE &&
                         at[MobileGr[j].atom_number].endpoint &&
                         jN1 < 0 ) {
                        jN1 = j;
                        n ++;
                    } else
                    if ( MobileGr[j].bond_type == BOND_TYPE_SINGLE &&
                         MobileGr[j].num_bonds == 2 &&
                         MobileGr[j].bonds_valence == 3 &&
                         MobileGr[j].forbidden == forbidden_mask &&
                         at[MobileGr[j].atom_number].endpoint &&
                         jN2 < 0 ) {
                        jN2 = j;
                        n ++;
                    } else
                    if ( MobileGr[j].bond_type == BOND_TYPE_SINGLE &&
                         MobileGr[j].num_bonds <= 2 &&
                         MobileGr[j].bonds_valence <= 3 &&
                         at[MobileGr[j].atom_number].endpoint &&
                         jX < 0 ) {
                        jX = j;
                        n ++;
                    }
                }
            }
            if ( n != 3 ) {
                goto case_5_9b_to_8b;
            }

            e = pBNS->edge + pVA[MobileGr[jN2].atom_number].nTautGroupEdge - 1;
            if ( e->flow ) {
                goto case_5_9b_to_8b; /* should not happen ??? */
            }
            pv1 = pBNS->vert + e->neighbor1;  /* must be jN2 */
            pv2 = pBNS->vert + (e->neighbor1 ^ e->neighbor12);
            ie = e - pBNS->edge;
            ie1 = ie2 = -1;
            for ( j = 0; j < pv1->num_adj_edges; j ++ ) {
                ev1 = pBNS->edge + pv1->iedge[j];
                if ( ev1->flow && !ev1->forbidden ) {
                    ie1 = ev1 - pBNS->edge;
                    pv1 = pBNS->vert + (ev1->neighbor12 ^ (pv1 - pBNS->vert));
                    break;
                }
            }
            for ( j = 0; j < pv2->num_adj_edges; j ++ ) {
                ev2 = pBNS->edge + pv2->iedge[j];
                if ( ev2->flow && !ev2->forbidden ) {
                    ie2 = ev2 - pBNS->edge;
                    pv2 = pBNS->vert + (ev2->neighbor12 ^ (pv2 - pBNS->vert));
                    break;
                }
            }
            if ( ie1 < 0 || ie2 < 0 ) {
                goto case_5_9b_to_8b;
            }
            e->flow ++;
            e->forbidden |= forbidden_mask;
            ev1->flow --;
            ev2->flow --;
            pv1->st_edge.flow --;
            pv2->st_edge.flow --;
            pBNS->tot_st_flow -= 2;
            num_changes ++;
            continue;
        }
case_5_9b_to_8b: /* #3 */
        /*********************************************************************/
        if ( 3 == num_bonds_non_metal &&
             4 == bonds_valence_non_metal &&
             (centerpoint_type == EL_TYPE_S) &&
             0 == num_O+num_S && 2 == num_N && 0 == num_P &&
             1 == num_forbidden &&
             2 == num_eql_mobile_gr && 0 == num_dif_mobile_gr &&
             1 == num_donor_endpoints && 1 == num_acceptor_endpoints &&
             1 == num_donors          && 1 == num_acceptors &&
             MobileGr[ind_forbidden].bond_type == BOND_TYPE_DOUBLE ) {
            /******************************************************** 
             ***         InChI Tech. Man., Table 5, case 9b->8b   *** 
             ******************************************************** 
                           --->                      O  =  O, S, Se, Te      
                  X    is an              X          S  =  S, Se, Te         
         \       /   endpoint    \       /           f  =  fixed bond        
          \ 1ev /  C-----ZH       \     /  C=====Z   tg =  Mobile-H vertex   
           N===S  ||     ||        N===S   |     |   X  =  is N not an endpoint;
      is an  f  \ ||     ||          f  \  |     |  -X is not terminal -O,-S,-Se,-Te or
   endpoint     2\||  e  ||              \ |  e  |        any N, P, As
                  N------tg               NH=====tg
               is an                          f      N, N, X, fixed double
            endpoint                                 =====================

            Problem:
              N1 and N2 are Mobile-H endpoints and belong to the same Mobile-H group.
              Fixed bond N1==S prevents the correct structure restoration:
              H cannot migrate ZH->N2->N1 because N1=S bond is fixed
            Solution:
              Move H from ZH to N2 to make N1=S bond unfixed and fix S-X bond (Table 5, case 8)
              (see fix_special_bonds(...) for details )
            *********************************************************/
            int jN1 = -1, jN2 = -1, jX = -1, n = 0;
            EdgeIndex ie, ie1, ie2;
            for ( j = 0; j < num_endpoints; j ++ ) {
                if ( (MobileGr[j].atom_type_pVA & EL_TYPE_MASK) == EL_TYPE_N ) {
                    if ( MobileGr[j].bond_type == BOND_TYPE_DOUBLE &&
                         at[MobileGr[j].atom_number].endpoint &&
                         (MobileGr[j].forbidden == forbidden_mask) &&
                         jN1 < 0 ) {
                        jN1 = j;
                        n ++;
                    } else
                    if ( MobileGr[j].bond_type == BOND_TYPE_SINGLE &&
                         MobileGr[j].num_bonds == 2 &&
                         MobileGr[j].bonds_valence == 3 &&
                         at[MobileGr[j].atom_number].endpoint &&
                         !(MobileGr[j].forbidden & forbidden_mask) &&
                         jN2 < 0 ) {
                        jN2 = j;
                        n ++;
                    }
                } else
                if ( !((MobileGr[j].atom_type_pVA & (EL_TYPE_N | EL_TYPE_P)) ||
                       ((MobileGr[j].atom_type_pVA & (EL_TYPE_O | EL_TYPE_S)) &&
                        MobileGr[j].num_bonds > 1) ) &&
                     !(MobileGr[j].forbidden & forbidden_mask) &&
                     MobileGr[j].bond_type == BOND_TYPE_SINGLE &&
                     jX < 0 ) {
                    jX = j;
                    n ++;
                }
            }
            if ( jN1 < 0 || jN2 < 0 ) {
                goto case_5_9c_to_8c;
            }

            ev = pBNS->edge + pBNS->vert[i].iedge[MobileGr[jN1].ineigh];
            e = pBNS->edge + pVA[MobileGr[jN2].atom_number].nTautGroupEdge - 1;
            if ( e->flow ) {
                goto case_5_9c_to_8c; /* should not happen ??? */
            }
            pv1 = pBNS->vert + e->neighbor1;  /* must be jN2 */
            pv2 = pBNS->vert + (e->neighbor1 ^ e->neighbor12);
            ie = e - pBNS->edge;
            ie1 = ie2 = -1;
            ev->forbidden &= inv_forbidden_mask;
            for ( j = 0; j < pv1->num_adj_edges; j ++ ) {
                ev1 = pBNS->edge + pv1->iedge[j];
                if ( ev1->flow && !ev1->forbidden ) {
                    ie1 = ev1 - pBNS->edge;
                    pv1 = pBNS->vert + (ev1->neighbor12 ^ (pv1 - pBNS->vert));
                    break;
                }
            }
            for ( j = 0; j < pv2->num_adj_edges; j ++ ) {
                ev2 = pBNS->edge + pv2->iedge[j];
                if ( ev2->flow && !ev2->forbidden ) {
                    ie2 = ev2 - pBNS->edge;
                    pv2 = pBNS->vert + (ev2->neighbor12 ^ (pv2 - pBNS->vert));
                    break;
                }
            }
            if ( ie1 < 0 || ie2 < 0 ) {
                ev->forbidden |= forbidden_mask; /* failed; restore the forbidden bit */ 
                goto case_5_9c_to_8c;
            }
            e->flow ++;
            e->forbidden |= forbidden_mask;
            ev1->flow --;
            ev2->flow --;
            pv1->st_edge.flow --;
            pv2->st_edge.flow --;
            pBNS->tot_st_flow -= 2;
            num_changes ++;
            continue;
        }
case_5_9c_to_8c: /* #4 */
        /*********************************************************************/
        if ( 3 == num_bonds_non_metal &&
             4 == bonds_valence_non_metal &&
             (centerpoint_type == EL_TYPE_S) &&
             0 == num_O+num_S && 3 == num_N && 0 == num_P &&
             0 == num_forbidden &&
             2 == num_diff_t_groups && 2 == num_mgroups && /* all neighbors belong to 2 t-groups */
             3 == num_eql_mobile_gr + num_dif_mobile_gr && /* all 3 neighbors belong to t-groups */
             2 == num_donor_endpoints && 1 == num_acceptor_endpoints &&
             2 == num_donors          && 1 == num_acceptors ) {
            /******************************************************** 
             ***         InChI Tech. Man., Table 5, case 8c->9c   *** 
             ******************************************************** 
                is an endpoint --->                  O  =  O, S, Se, Te      
                tg1   NH2(X)             NH2(X)      S  =  S, Se, Te         
             \       / pv1  pv2 \       /            f  =  fixed bond        
              \(1)  /  C=====Z   \     /  C-----ZH   tg =  Mobile-H vertex   
               N===S   |     |    N===S  ||     ||   X  =  is N not an endpoint;
            is an   \  |ev1  | ev2     \ ||     ||  -X is not terminal -O,-S,-Se,-Te or
          endpoint   \ |  e  |          \||  e  ||        any N, P, As
          tg1      (2)NH=====tg          N------tg   C is a centerpoint of a t-group
                     is an                   f       N, N, X, fixed single
                   endpoint                          =====================
                   tg2
            Problem:
              N (1) and NH2 are Mobile-H group 1 endpoiints, NH (2) is a Mobile-H group 2 endpoint
              Unfixed bonds N==S--NH(2) allows the two Mobile H groups to merge
              hence prevents the correct structure restoration:
              H can migrate from NH (2) to N (1) because S-NH(1) bond is not fixed
            Solution:
              Move H from NH(2) to =Z to make S-NH(2) bond fixed (Table 5, case 8c)
              (this unfixes the bond, see fix_special_bonds(...) )
            *********************************************************/
            int jN1 = -1, jN2 = -1, jX = -1, n = 0;
            /* find t-group that is represented by only one neighbor */
            for ( j = 0, k = 0; j < num_mgroups; j ++ ) {
                if ( 1 == MGroups[k].num && MGroups[k].group_number ) {
                    k = MGroups[k].group_number;
                    break;
                }
            }
            if ( !k ) {
                goto case_5_9c_to_9d;
            }
            for ( j = 0; j < num_endpoints; j ++ ) {
                if ( (MobileGr[j].atom_type_pVA & EL_TYPE_MASK) == EL_TYPE_N ) {
                    if ( MobileGr[j].bond_type == BOND_TYPE_DOUBLE &&
                         at[MobileGr[j].atom_number].endpoint &&
                         at[MobileGr[j].atom_number].endpoint != k &&
                         jN1 < 0 ) {
                        jN1 = j;
                        n ++;
                    } else
                    if ( MobileGr[j].bond_type == BOND_TYPE_SINGLE &&
                         MobileGr[j].num_bonds == 2 &&
                         MobileGr[j].bonds_valence == 2 &&
                         at[MobileGr[j].atom_number].endpoint == k &&
                         jN2 < 0 ) {
                        jN2 = j;
                        n ++;
                    } else
                    if ( MobileGr[j].bond_type == BOND_TYPE_SINGLE &&
                         MobileGr[j].num_bonds <= 2 &&
                         MobileGr[j].bonds_valence <= 3 &&
                         at[MobileGr[j].atom_number].endpoint &&
                         at[MobileGr[j].atom_number].endpoint != k &&
                         jX < 0 ) {
                        jX = j;
                        n ++;
                    }
                }
            }
            if ( n != 3 ) {
                goto case_5_9c_to_9d;
            }

            e = pBNS->edge + pVA[MobileGr[jN2].atom_number].nTautGroupEdge - 1;
            if ( !e->flow ) {
                goto case_5_9c_to_9d; /* should not happen ??? */
            }
            e->flow --;
            pBNS->vert[e->neighbor1].st_edge.flow --;
            pBNS->vert[e->neighbor1 ^ e->neighbor12].st_edge.flow --;
            pBNS->tot_st_flow -= 2;
            e->forbidden |= forbidden_mask;
            num_changes ++;
            continue;
        }
case_5_9c_to_9d: /* #6 */
        /*********************************************************************/
        if ( 3 == num_bonds_non_metal &&
             4 == bonds_valence_non_metal &&
             (centerpoint_type == EL_TYPE_S) &&
             0 == num_O+num_S && 3 == num_N && 0 == num_P &&
             0 == num_forbidden &&
             3 == num_mgroups && 2 == num_diff_t_groups &&
             2 == num_donor_endpoints && 0 == num_acceptor_endpoints &&
             2 == num_donors          && 1 == num_acceptors ) {
            /******************************************************** 
             ***         InChI Tech. Man., Table 5, case 9b->8b   *** 
             ******************************************************** 
                 3(X)      --->         e2|          O  =  O, S, Se, Te      
                  NH---is an              N======    S  =  S, Se, Te         
         \       /   endpoint    \       /           f  =  fixed bond        
          \ 1   /  C=====Z        \  f  /  C-----Z   tg =  Mobile-H vertex   
           N===S   |     |         N===S  ||     ||  X  =  is N not an endpoint;
     is an  ev  \  |     |           ev \ ||     || -X is not terminal -O,-S,-Se,-Te or
   endpoint     2\ |  e1 |               \||  e1 ||       any N, P, As
                  NH=====tg               N------tg
               is an                          f      N, N, X, fixed double
            endpoint                                 =====================

            Problem:
              N1, N2, and N3 are Mobile-H endpoints and belong to the same Mobile-H group.
              Fixed bond N1==S prevents the correct structure restoration:
              H cannot migrate N3->N2->N1 because N1=S bond is fixed
            Solution:
              Move mobile H to N2 and N3 to make N1=S bond unfixed (Table 5, case 9c)
              (see fix_special_bonds(...) for details )
            *********************************************************/
            int jN1 = -1, jN2 = -1, jX = -1, n = 0;
            for ( j = 0; j < num_endpoints; j ++ ) {
                if ( (MobileGr[j].atom_type_pVA & EL_TYPE_MASK) == EL_TYPE_N ) {
                    if ( MobileGr[j].bond_type == BOND_TYPE_DOUBLE &&
                         !at[MobileGr[j].atom_number].endpoint &&
                         !(MobileGr[j].forbidden & forbidden_mask) &&
                         jN1 < 0 ) {
                        jN1 = j;
                        n ++;
                    } else
                    if ( MobileGr[j].bond_type == BOND_TYPE_SINGLE &&
                         MobileGr[j].num_bonds == 2 &&
                         MobileGr[j].bonds_valence <= 3 &&
                         at[MobileGr[j].atom_number].endpoint &&
                         !(MobileGr[j].forbidden & forbidden_mask) ) {
                        if ( jN2 < 0 ) {
                            jN2 = j;
                            n ++;
                        } else
                        if ( jX < 0 ) {
                            jX = j;
                            n ++;
                        }
                    }
                }
            }
            if ( n != 3 ) {
                goto case_end;
            }
            ev = pBNS->edge + pBNS->vert[i].iedge[MobileGr[jN1].ineigh];
            if ( !e->flow ) {
                goto case_end;
            }
            e1 = pBNS->edge + pVA[MobileGr[jN2].atom_number].nTautGroupEdge - 1;
            if ( !e1->flow ) {
                goto case_end; /* should not happen ??? */
            }
            e2 = pBNS->edge + pVA[MobileGr[jX].atom_number].nTautGroupEdge - 1;
            if ( !e2->flow ) {
                goto case_end; /* should not happen ??? */
            }
            /* take care of edge e1 */
            e = e1;
            e->flow --;
            pBNS->vert[e->neighbor1].st_edge.flow --;
            pBNS->vert[e->neighbor1 ^ e->neighbor12].st_edge.flow --;
            pBNS->tot_st_flow -= 2;
            e->forbidden |= forbidden_mask;
            num_changes ++;
            /* take care of edge e2 */
            e = e2;
            e->flow --;
            pBNS->vert[e->neighbor1].st_edge.flow --;
            pBNS->vert[e->neighbor1 ^ e->neighbor12].st_edge.flow --;
            pBNS->tot_st_flow -= 2;
            e->forbidden |= forbidden_mask;
            num_changes ++;
            /* take care of edge ev: do not let it change */
            ev->forbidden |= forbidden_mask;
            continue;
        }
case_end:;
    }
/*exit_function:*/
    return num_changes;
}
/******************************************************************************************************/
/* Replace ambiguous neutral (+)edge->flow=0, (-)edge->flow=1 with (+)edge->flow=1, (-)edge->flow=0   */
/******************************************************************************************************/
int RearrangePlusMinusEdgesFlow( BN_STRUCT *pBNS, BN_DATA *pBD, VAL_AT *pVA,
                                 ALL_TC_GROUPS *pTCGroups, int forbidden_edge_mask )
{
    int ret, ePlus, eMinus;
    EDGE_LIST NewlyFixedEdges;
    BNS_EDGE *pPlus, *pMinus;
    int i, k1, k2, num_found, num_tot, delta, v1, v2;

    ret = 0;

    AllocEdgeList( &NewlyFixedEdges, EDGE_LIST_CLEAR );
    for ( i = 0, num_found = num_tot = 0; i < pBNS->num_atoms; i ++ ) {
        eMinus = pVA[i].nCMinusGroupEdge - 1;
        ePlus  = pVA[i].nCPlusGroupEdge - 1;
        num_tot += (eMinus >= 0) + (ePlus >= 0);
        if ( eMinus >= 0 && ePlus >= 0 ) {
            pPlus  = pBNS->edge + ePlus;
            pMinus = pBNS->edge + eMinus;
            if ( (k1=pMinus->flow) > 0 &&  (k2=pPlus->cap-pPlus->flow) > 0 ) {
                num_found ++;
            }
        } 
    }
    if ( !num_found ) {
        goto exit_function;
    }
    if ( (ret = AllocEdgeList( &NewlyFixedEdges, num_tot + pBNS->num_bonds )) ) {
        goto exit_function;
    }

    for ( i = 0, num_found = num_tot = 0; i < pBNS->num_atoms; i ++ ) {
        eMinus = pVA[i].nCMinusGroupEdge - 1;
        ePlus  = pVA[i].nCPlusGroupEdge - 1;
        num_tot += (eMinus >= 0) + (ePlus >= 0);
        if ( eMinus >= 0 && ePlus >= 0 ) {
            pPlus  = pBNS->edge + ePlus;
            pMinus = pBNS->edge + eMinus;
            if ( (k1=pMinus->flow) > 0 &&  (k2=pPlus->cap  - pPlus->flow) > 0 ) {
                /* rearrange */
                v1 = pMinus->neighbor1;
                v2 = pMinus->neighbor12 ^ v1;
                delta = inchi_min(k1,k2);
                pMinus->flow -= delta;
                pBNS->vert[v1].st_edge.flow -= delta;
                pBNS->vert[v2].st_edge.flow -= delta;
                pBNS->tot_st_flow -= 2*delta;
            }
            /* fix charges */
            pPlus->forbidden  |= forbidden_edge_mask;
            pMinus->forbidden |= forbidden_edge_mask;
            if ( (ret = AddToEdgeList( &NewlyFixedEdges, eMinus, 0 )) ||
                 (ret = AddToEdgeList( &NewlyFixedEdges, ePlus, 0 ))) {
                goto exit_function;
            }
        } else
        if ( eMinus >= 0 ) {
            /* fix charges */
            pMinus = pBNS->edge + eMinus;
            pMinus->forbidden |= forbidden_edge_mask;
            if ( (ret = AddToEdgeList( &NewlyFixedEdges, eMinus, 0 )) ) {
                goto exit_function;
            }
        } else
        if ( ePlus >= 0 ) {
            /* fix charges */
            pPlus  = pBNS->edge + ePlus;
            pPlus->forbidden  |= forbidden_edge_mask;
            if ( (ret = AddToEdgeList( &NewlyFixedEdges, ePlus, 0 )) ) {
                goto exit_function;
            }
        }
    }
    for ( i = 0; i < pBNS->num_bonds; i ++ ) {
        pBNS->edge[i].forbidden |= forbidden_edge_mask;
        if ( (ret = AddToEdgeList( &NewlyFixedEdges, i, 0 )) ) {
            goto exit_function;
        }
    }
    ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
    RemoveForbiddenEdgeMask( pBNS, &NewlyFixedEdges, forbidden_edge_mask );
    if ( ret < 0 ) {
        goto exit_function;
    }

exit_function:
    AllocEdgeList( &NewlyFixedEdges, EDGE_LIST_FREE );
    return ret;
}
/******************************************************************************************************/
int IncrementZeroOrderBondsToHeteroat( BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct,
                                            inp_ATOM *at, inp_ATOM *at2,
                                            VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups,
                                            int *pnNumRunBNS, int *pnTotalDelta,
                                            int forbidden_edge_mask)
{
#define FIX_BOND_ADD_ALLOC 128    
    Vertex     vPathStart, vPathEnd;
    int        nPathLen, nDeltaH, nDeltaCharge, nNumVisitedAtoms;
    BNS_EDGE   *pe, *peZero, *peNeighMeigh = NULL, *peMeFlower;
    BNS_VERTEX *pMeFlower = NULL, *pNeigh = NULL, *pNeighNeigh=NULL;
    
    int i, j, k, ret2, ret, bFixedCarbonCharges, num_changes, bSuccess;
    int num_at = pStruct->num_atoms;
    int num_deleted_H = pStruct->num_deleted_H;
    int len_at = num_at + num_deleted_H;
    int inv_forbidden_edge_mask = ~forbidden_edge_mask;
    Vertex vMeFlower0, vNeigh, vNeighMeigh = NO_VERTEX;
    
    EDGE_LIST CarbonChargeEdges;
    EDGE_LIST NewlyFixedEdges;

    ret = 0;
    num_changes = 0;

    bFixedCarbonCharges = 0;
    AllocEdgeList( &CarbonChargeEdges, EDGE_LIST_CLEAR );
    AllocEdgeList( &NewlyFixedEdges, EDGE_LIST_CLEAR );

    if ( !pTCGroups->num_metal_atoms || 
         0 > (k=pTCGroups->nGroup[TCG_MeFlower0]) ||
         0 > (vMeFlower0 = pTCGroups->pTCG[k].nVertexNumber)) {
        goto exit_function;
    }

    memcpy( at2, at, len_at*sizeof(at2[0]));
    pStruct->at = at2;
    ret2 = CopyBnsToAtom( pStruct, pBNS, pVA, pTCGroups, 1 );
    pStruct->at = at;
    if ( ret2 < 0 ) {
        ret = ret2;
        goto exit_function;
    }
    for ( i = 0; i < num_at; i ++ ) {
        if ( !pVA[i].cMetal || pVA[i].nMetalGroupEdge <= 0 ) {
            continue;
        }
        peMeFlower = pBNS->edge + pVA[i].nMetalGroupEdge-1;
        if ( vMeFlower0 != (peMeFlower->neighbor12 ^ i) ) {
            ret = RI_ERR_PROGR;
            goto exit_function;
        }
        pMeFlower = pBNS->vert + vMeFlower0;

        for ( j = 0; j < at2[i].valence; j ++ ) {
            if ( !peMeFlower->flow ) {
                break; /* cannot do anything */
            }
            if ( !(at2[i].bond_type[j] & BOND_TYPE_MASK) ) {
                /* found a zero order bond */
                if ( !bFixedCarbonCharges ) {
                    /* do not let carbon atoms get charged */
                    if ( 0 > (ret = ForbidCarbonChargeEdges( pBNS, pTCGroups, &CarbonChargeEdges, forbidden_edge_mask ))) {
                        goto exit_function;
                    }
                    bFixedCarbonCharges ++;
                }
                peZero = pBNS->edge + pBNS->vert[i].iedge[j];
                if ( peZero->flow ) {
                    ret = RI_ERR_PROGR;
                    goto exit_function;
                }
                /* fix other edges */
                for ( k = 0; k < at2[i].valence; k ++ ) {
                    pe = pBNS->edge + pBNS->vert[i].iedge[k];
                    if ( pe->flow == 1 && !(pe->forbidden & forbidden_edge_mask) ) {
                        if ( (ret = AddToEdgeList( &NewlyFixedEdges, pe - pBNS->edge, FIX_BOND_ADD_ALLOC )) ) {
                            goto exit_function;
                        }
                        pe->forbidden |= forbidden_edge_mask;
                    }
                }
                /* do not create =N(+)= in a ring or #O(+) terminal */
                for ( k = 0; k < num_at; k ++ ) {
                    if ( !pVA[k].cMetal && pVA[k].cNumValenceElectrons == 5 &&
                         at2[k].valence == 2 && !at2[k].num_H && pVA[k].cMinRingSize <= 6 &&
                         pVA[k].nCPlusGroupEdge > 0 &&
                         (pe=pBNS->edge + pVA[k].nCPlusGroupEdge-1)->flow==1 &&
                         !(pe->forbidden & forbidden_edge_mask)) {

                        if ( (ret = AddToEdgeList( &NewlyFixedEdges, pe - pBNS->edge, FIX_BOND_ADD_ALLOC )) ) {
                            goto exit_function;
                        }
                        pe->forbidden |= forbidden_edge_mask;
                    }
                }

                /* metal's neighbor connected by a zero-order bond */
                pNeigh = pBNS->vert + (vNeigh = at2[i].neighbor[j]);
                /*for ( k = 0; k < pNeigh->num_adj_edges; k ++ )*/
                for ( k = pNeigh->num_adj_edges-1; 0 <= k; k -- )
                {
                    peNeighMeigh = pBNS->edge + pNeigh->iedge[k];
                    if ( !peNeighMeigh->flow ) {
                        continue;
                    }
                    vNeighMeigh = peNeighMeigh->neighbor12 ^ vNeigh;
                    if ( vNeighMeigh != i && vNeighMeigh != vMeFlower0 ) {
                        /* metal neighbor's neighbor connected by a not-zero-order bond */
                        pNeighNeigh = pBNS->vert + vNeighMeigh;
                        break; /* found */
                    }
                }
                if ( k < 0 ) {
                    continue; /* neighbor not found */
                }
                peZero->flow ++;
                peZero->forbidden |= forbidden_edge_mask;
                peMeFlower->flow --;
                peNeighMeigh->flow --;
                pMeFlower->st_edge.flow --;
                pNeighNeigh->st_edge.flow --;
                pBNS->tot_st_flow -= 2;
                /* test */
                bSuccess = 0;
                ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                      &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

                if ( ret == 1 && ((vPathEnd == vMeFlower0 && vPathStart == vNeighMeigh) ||
                                  (vPathEnd == vNeighMeigh && vPathStart == vMeFlower0)) && abs(nDeltaCharge) <= 2 ) {
                    ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                    if ( ret > 0 ) {
                        (*pnNumRunBNS) ++;
                        *pnTotalDelta += ret;
                        num_changes ++;
                        bSuccess = ret;
                    }
                    if ( (ret = AddToEdgeList( &NewlyFixedEdges, peZero - pBNS->edge, FIX_BOND_ADD_ALLOC )) ) {
                        goto exit_function;
                    }
                } else {
                    peZero->flow --;
                    peZero->forbidden &= inv_forbidden_edge_mask;
                    peMeFlower->flow ++;
                    peNeighMeigh->flow ++;
                    pMeFlower->st_edge.flow ++;
                    pNeighNeigh->st_edge.flow ++;
                    pBNS->tot_st_flow += 2;
                }
                RemoveForbiddenEdgeMask( pBNS, &NewlyFixedEdges, forbidden_edge_mask );
                NewlyFixedEdges.num_edges = 0;
                if ( bSuccess ) {
                    /* update at2[] */
                    memcpy( at2, at, len_at*sizeof(at2[0]));
                    pStruct->at = at2;
                    ret2 = CopyBnsToAtom( pStruct, pBNS, pVA, pTCGroups, 1 );
                    pStruct->at = at;
                    if ( ret2 < 0 ) {
                        ret = ret2;
                        goto exit_function;
                    }
                }
            }
        }
    }
    ret = num_changes;

exit_function:
    RemoveForbiddenEdgeMask( pBNS, &CarbonChargeEdges, forbidden_edge_mask );
    RemoveForbiddenEdgeMask( pBNS, &NewlyFixedEdges, forbidden_edge_mask );
    AllocEdgeList( &CarbonChargeEdges, EDGE_LIST_FREE );
    AllocEdgeList( &NewlyFixedEdges, EDGE_LIST_FREE );
    return ret;
}
/***********************************************************************
        NH2                NH2
           \                  \
            C==S(+)-   =>      C(+)-S-   where NH2 are not tautomeric
           /                  /
        NH2                NH2
************************************************************************/
int MovePlusFromS2DiaminoCarbon( BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct,
                    inp_ATOM *at, inp_ATOM *at2,
                    VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups,
                    int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask)
{
    int i, j, k, ret, ret2, cur_success;
    int delta;
    EdgeIndex        ePlusS, ePlusC, eMinusC, e;
    BNS_VERTEX *pvS, *pvC, *pv1, *pv2;
    BNS_EDGE   *pePlusS, *pePlusC, *pe1, *pe2, *peCN[3], *peSC, *pe;
    Vertex     vC, vN;
    
    int num_at = pStruct->num_atoms;
    int num_deleted_H = pStruct->num_deleted_H;
    int len_at = num_at + num_deleted_H;
    
    Vertex     vPathStart, vPathEnd, v1, v2;
    int        nPathLen, nDeltaH, nDeltaCharge, nNumVisitedAtoms;

    EDGE_LIST AllChargeEdges;
    
    ret = 0;
    cur_success = 0;

    AllocEdgeList( &AllChargeEdges, EDGE_LIST_CLEAR );

    memcpy( at2, at, len_at*sizeof(at2[0]));
    pStruct->at = at2;
    ret2 = CopyBnsToAtom( pStruct, pBNS, pVA, pTCGroups, 1 );
    pStruct->at = at;
    if ( ret2 < 0 ) {
        ret = ret2;
        goto exit_function;
    }
    /* find (NH2)C=S(+) */
    for ( i = 0; i < num_at; i ++ ) {
        pvS = pBNS->vert+i;
        ePlusS = pVA[i].nCPlusGroupEdge-1;
        pePlusS=pBNS->edge+ePlusS;
        pe1=pBNS->edge + pvS->iedge[0];
        pe2=pBNS->edge + pvS->iedge[1];
        peSC=pe1->flow? pe1 : pe2;
        vC = peSC->neighbor12 ^ i;
        ePlusC=pVA[vC].nCPlusGroupEdge-1;
        pePlusC=pBNS->edge+ePlusC;
        eMinusC=pVA[vC].nCMinusGroupEdge-1;
        if ( !pVA[i].cMetal && pVA[i].cNumValenceElectrons == 6 &&
             at2[i].valence == 2 && 
             pvS->st_edge.cap == pvS->st_edge.flow &&
             0 <= ePlusS && !pePlusS->flow && /* S(+) */
             pe1->flow + 
             pe2->flow == 1 /* -S(+)= */ &&
             pVA[vC].cNumValenceElectrons == 4 &&
             at2[vC].valence == 3 &&
             0 <= ePlusC && pePlusC->flow &&
             !(0 <= eMinusC && pBNS->edge[eMinusC].flow ) ) {
            /* found >C=S(+)- */
            pvC = pBNS->vert + vC;
            for ( j = k = 0; j < at[vC].valence; j ++ ) {
                if ( peSC != (peCN[k] = pBNS->edge + pvC->iedge[j]) && !peCN[k]->flow ) {
                    k ++; /* a single bond from C */
                }
            }
            if ( k != 2 ) {
                continue;
            }
            for ( j = 0; j < k; j ++ ) {
                vN = peCN[j]->neighbor12 ^ vC;
                if ( pVA[vN].cNumValenceElectrons != 5 ||
                     pBNS->vert[vN].st_edge.cap != pBNS->vert[vN].st_edge.flow ||
                     at2[vN].num_H != 2 ||
                     at2[vN].endpoint || (pStruct->endpoint && pStruct->endpoint[vN]) ) {
                    break; /* does not fit the pattern */
                }
            }
            if ( j != k ) {
                continue;
            }
            /* fix all charges */
            if ( !AllChargeEdges.num_edges ) {
                for ( j = 0; j < num_at; j ++ ) {
                    if ( 0 <= (e = pVA[j].nCPlusGroupEdge-1) && !pBNS->edge[e].forbidden &&
                         (ret = AddToEdgeList( &AllChargeEdges, e, 2*num_at )) ) {
                        goto exit_function;
                    }
                    if ( 0 <= (e = pVA[j].nCMinusGroupEdge-1) && !pBNS->edge[e].forbidden &&
                         (ret = AddToEdgeList( &AllChargeEdges, e, 2*num_at )) ) {
                        goto exit_function;
                    }
                }
            }
            SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
            pePlusS->forbidden &= ~forbidden_edge_mask;
            pe   = pePlusC;
            if ( !pe->flow )
                continue;
            delta = 1;
            pv1 = pBNS->vert + (v1 = pe->neighbor1);
            pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);

            pe->flow -= delta;
            pv1->st_edge.flow -= delta;
            pv2->st_edge.flow -= delta;
            pBNS->tot_st_flow -= 2*delta;

            ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                  &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

            if ( ret == 1 && ((vPathEnd == v1 && vPathStart == v2) ||
                              (vPathEnd == v2 && vPathStart == v1)) && nDeltaCharge == -1 ) {
                /* Remover (+)charge from S => nDeltaCharge == -1 */
                /* Flow change on pe (+)charge edge (atom S) is not known to RunBnsTestOnce()) */
                ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                if ( ret > 0 ) {
                    (*pnNumRunBNS) ++;
                    cur_success ++;
                }
            } else {
                pe->flow += delta;
                pv1->st_edge.flow += delta;
                pv2->st_edge.flow += delta;
                pBNS->tot_st_flow += 2*delta;
            }
            RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
        }
    }
exit_function:
    AllocEdgeList( &AllChargeEdges, EDGE_LIST_FREE );
    return ret;
}
/******************************************************************************************************/
int EliminateChargeSeparationOnHeteroatoms( BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct,
                                            inp_ATOM *at, inp_ATOM *at2,
                                            VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups,
                                            int *pnNumRunBNS, int *pnTotalDelta,
                                            int forbidden_edge_mask, int forbidden_stereo_edge_mask)
    /********* Avoid charge separation on heteroatoms ******************/
{
    int i, j, k, ret, ret2, num_pos, num_neg, num_min=0, bFixedCarbonCharges;
    int vPlusSuper;  /* (+)super vertex */
    int ePlusSuper;  /* edge from vPlusSuper to (+/-) */
    int vPlusMinus;  /* (+/-) vertex */
    int nDeltaPlus1, nDeltaMinus1, delta;
    BNS_VERTEX *pvPlusSuper, *pvPlusMinus;
    BNS_EDGE   *pEdge;
    
    int num_at = pStruct->num_atoms;
    int num_deleted_H = pStruct->num_deleted_H;
    int len_at = num_at + num_deleted_H;
    int inv_forbidden_edge_mask = ~forbidden_edge_mask;
    
    Vertex     vPathStart, vPathEnd;
    int        nPathLen, nDeltaH, nDeltaCharge, nNumVisitedAtoms;

    EDGE_LIST FixedLargeRingStereoEdges, CarbonChargeEdges;
    
    ret = 0;

    AllocEdgeList( &FixedLargeRingStereoEdges, EDGE_LIST_CLEAR );
    AllocEdgeList( &CarbonChargeEdges, EDGE_LIST_CLEAR );
    bFixedCarbonCharges = 0;

    if ( forbidden_stereo_edge_mask ) {
        for ( i = 0; i < num_at; i ++ ) {
            for ( j = 0; j < at2[i].valence; j ++ ) {
                if ( pBNS->edge[k = pBNS->vert[i].iedge[j]].forbidden == forbidden_stereo_edge_mask ) {
                    int nMinRingSize = is_bond_in_Nmax_memb_ring( at2, i, j, pStruct->pbfsq->q,
                                                             pStruct->pbfsq->nAtomLevel,
                                                             pStruct->pbfsq->cSource, 99 /* max ring size */ );
                    if ( 0 < nMinRingSize && (ret = AddToEdgeList( &FixedLargeRingStereoEdges, k, 64 ))) {
                        goto exit_function;
                    }
                }
            }
        }
        if ( !FixedLargeRingStereoEdges.num_edges ) {
            goto exit_function;
        } else {
            /* allow stereobonds in rings change */
            RemoveForbiddenEdgeMask( pBNS, &FixedLargeRingStereoEdges, forbidden_stereo_edge_mask );
        }
    }

    memcpy( at2, at, len_at*sizeof(at2[0]));
    pStruct->at = at2;
    ret2 = CopyBnsToAtom( pStruct, pBNS, pVA, pTCGroups, 1 );
    pStruct->at = at;
    if ( ret2 < 0 ) {
        ret = ret2;
        goto exit_function;
    }
    /* count charges */
    num_pos = num_neg = 0;
    for ( i = 0; i < num_at; i ++ ) {
        if ( !pVA[i].cMetal && !at2[i].radical ) {
            num_pos += ( at2[i].charge > 0 );
            num_neg += ( at2[i].charge < 0 );
        }
    }
    num_min = inchi_min( num_pos, num_neg );


    if ( num_min && 
         (k          = pTCGroups->nGroup[TCG_Plus]) >= 0 &&
         (ePlusSuper = pTCGroups->pTCG[k].nForwardEdge) > 0 &&
         (vPlusSuper = pTCGroups->pTCG[k].nVertexNumber) >= num_at &&
         !(pEdge=pBNS->edge + ePlusSuper)->forbidden ) {

        vPlusMinus = pEdge->neighbor12 ^ vPlusSuper;
        pvPlusSuper = pBNS->vert + vPlusSuper;
        pvPlusMinus = pBNS->vert + vPlusMinus;
        num_min     = inchi_min( num_min, pEdge->flow );
        nDeltaPlus1  = pvPlusSuper->st_edge.cap - pvPlusSuper->st_edge.flow;
        nDeltaMinus1 = pvPlusMinus->st_edge.cap - pvPlusMinus->st_edge.flow;
        if ( num_min && (!nDeltaPlus1 && !nDeltaMinus1 ) ) {
            if ( !bFixedCarbonCharges ) { /* 02-02-2006 */
                /* do not let carbon atoms get charged */
                if ( 0 > (ret = ForbidCarbonChargeEdges( pBNS, pTCGroups, &CarbonChargeEdges, forbidden_edge_mask ))) {
                    goto exit_function;
                }
                bFixedCarbonCharges ++;
            }
            delta = 1;
            pEdge->forbidden          |= forbidden_edge_mask;
            pBNS->edge_forbidden_mask |= forbidden_edge_mask;
            for ( i = 0; i < num_min; i += delta ) {
                /* cancel 1 pair of charges at a time   */
                /* an attempt to cancel all at once may */
                /* convert a pair of N(IV)(+) into a pair of N(V) neutral with total charge reduced by 2 */
                pEdge->flow               -= delta;
                pvPlusSuper->st_edge.flow -= delta;
                pvPlusMinus->st_edge.flow -= delta;
                pBNS->tot_st_flow         -= 2*delta;
                /* test for charhe cancellation */
                ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                      &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );
                if ( ret < 0 ) {
                    goto exit_function;
                }
                if ( ret == 1 && ((vPathEnd == vPlusSuper && vPathStart == vPlusMinus) ||
                                  (vPathEnd == vPlusMinus && vPathStart == vPlusSuper)) && nDeltaCharge < 0 ) {
                    ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                    (*pnNumRunBNS) ++;
                    if ( ret < 0 ) {
                        goto exit_function;
                    } else
                    if ( ret == 1 ) {
                        *pnTotalDelta += ret;
                    } else {
                        ret = RI_ERR_PROGR;
                        goto exit_function;
                    }
                } else {
                    pEdge->flow               += delta;
                    pvPlusSuper->st_edge.flow += delta;
                    pvPlusMinus->st_edge.flow += delta;
                    pBNS->tot_st_flow         += 2*delta;
                    break;
                }
            }
            num_min -= i; /* how many pairs of charges left */
            pEdge->forbidden &= inv_forbidden_edge_mask;
        }
        nDeltaPlus1  = pvPlusSuper->st_edge.cap - pvPlusSuper->st_edge.flow;
        nDeltaMinus1 = pvPlusMinus->st_edge.cap - pvPlusMinus->st_edge.flow;
        if ( num_min > 1 && (!nDeltaPlus1 && !nDeltaMinus1 ) ) {
            delta = 2;
            pEdge->forbidden          |= forbidden_edge_mask;
            pBNS->edge_forbidden_mask |= forbidden_edge_mask;
            for ( i = 0; i < num_min; i += delta ) {
                /* cancel 2 pairs of opposite charges at a time   */
                /* 1. test cancellation of a pair of (+) charges */
                pvPlusSuper->st_edge.cap += delta;
                pBNS->tot_st_cap         += delta;
                ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                      &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );
                if ( ret < 0 ) {
                    goto exit_function;
                }
                pvPlusSuper->st_edge.cap -= delta;
                pBNS->tot_st_cap         -= delta;
                if ( ret != 1 || (vPathEnd != vPlusSuper || vPathStart != vPlusSuper) || nDeltaCharge >= 0 ) {
                    break;
                }
                /* 2. test cancellation of a pair of (-) charges */
                pvPlusMinus->st_edge.cap += delta;
                pBNS->tot_st_cap         += delta;
                ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                      &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );
                if ( ret < 0 ) {
                    goto exit_function;
                }
                pvPlusMinus->st_edge.cap -= delta;
                pBNS->tot_st_cap         -= delta;
                if ( ret != 1 || (vPathEnd != vPlusMinus || vPathStart != vPlusMinus) || nDeltaCharge >= 0 ) {
                    break;
                }
                /* 3. Actually cancel the pair of charges */
                pEdge->flow               -= delta;
                pvPlusSuper->st_edge.flow -= delta;
                pvPlusMinus->st_edge.flow -= delta;
                pBNS->tot_st_flow         -= 2*delta;
                ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                (*pnNumRunBNS) ++;
                if ( ret < 0 ) {
                    goto exit_function;
                } else
                if ( ret == 2 ) {
                    *pnTotalDelta += ret;
                } else {
                    ret = RI_ERR_PROGR;
                    goto exit_function;
                }
            }
            num_min -= i; /* how many pairs of charges left */
            pEdge->forbidden &= inv_forbidden_edge_mask;
        }
    }
    memcpy( at2, at, len_at*sizeof(at2[0]));
    pStruct->at = at;
exit_function:
    if ( bFixedCarbonCharges ) {
        RemoveForbiddenEdgeMask(pBNS, &CarbonChargeEdges, forbidden_edge_mask ); 
        AllocEdgeList( &CarbonChargeEdges, EDGE_LIST_FREE );
    }
    if ( forbidden_stereo_edge_mask && FixedLargeRingStereoEdges.num_edges ) {
        SetForbiddenEdgeMask( pBNS, &FixedLargeRingStereoEdges, forbidden_stereo_edge_mask );
    }
    AllocEdgeList( &FixedLargeRingStereoEdges, EDGE_LIST_FREE );

    return ret < 0? ret : num_min;
}
#if (MOVE_CHARGES_FROM_HETEREO_TO_METAL == 1 )
/********************** not used *************************************************************************/
int MoveChargeFromHeteroatomsToMetals( BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct,
                                       inp_ATOM *at, inp_ATOM *at2,
                                       VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups,
                                       int *pnNumRunBNS, int *pnTotalDelta,
                                       int forbidden_edge_mask)
    /********* Avoid charge separation on heteroatoms ******************/
{
    int i, k, ret, ret2, num_pos, num_neg, num_min;
    int vPlusSuper, vMinusSuper;  /* (+), (-) super vertices */
    int ePlusSuper, eMinusSuper;  /* edges from vPlusSuper or vMinusSuper to (+/-) */
    int vPlMn;  /* (+/-) vertex */
    int vPlusHeteroat, vMinusHeteroat; /* (+), (-) heteroatom vertices */
    int ePlusHeteroat, eMinusHeteroat; /* edges from (+) or (-) heteroatom vertex to super (+) or (-) */
    int vPlusCarbons,  vMinusCarbons; /* (+), (-) carbons vertices */
    int ePlusCarbons,  eMinusCarbons; /* edges from (+), (-) carbons vertices to super (+) or (-) */
    int vPlusMetals,   vMinusMetals; /* (+), (-) carbons vertices */
    int ePlusMetals,   eMinusMetals; /* edges from (+), (-) carbons vertices to super (+) or (-) */
    int eMinusHeteroToSuper;         /* edge (-)vHetero-[eMinusHeteroat]-Y-[eMinusHeteroToSuper]-(-)vPlusSuper */
    int v1, v2;
    int nDeltaPlus1, nDeltaMinus1, delta;
    BNS_VERTEX *pvPlusSuper, *pvMinusSuper, *pvPlMn;
    BNS_VERTEX *pvPlusHeteroat, *pvMinusHeteroat, *pvPlusCarbons, *pvMinusCarbons;
    BNS_VERTEX *pvPlusMetals,  *pvMinusMetals;
    BNS_EDGE   *pEdgePlusHeteroat, *pEdgeMinusHeteroat, *pEdgeMinusHeteroToSuper;
    BNS_EDGE   *pEdgePlusCarbons, *pEdgeMinusCarbons, *pEdgePlusMetals, *pEdgeMinusMetals;
    BNS_EDGE   *pEdgePlusSuper, *pEdgeMinusSuper;
    
    int num_at = pStruct->num_atoms;
    int num_deleted_H = pStruct->num_deleted_H;
    int len_at = num_at + num_deleted_H;
    int inv_forbidden_edge_mask = ~forbidden_edge_mask;
    
    Vertex     vPathStart, vPathEnd;
    int        nPathLen, nDeltaH, nDeltaCharge, nNumVisitedAtoms;


    ret = 0;

    memcpy( at2, at, len_at*sizeof(at2[0]));
    pStruct->at = at2;
    ret2 = CopyBnsToAtom( pStruct, pBNS, pVA, pTCGroups, 1 );
    pStruct->at = at;
    if ( ret2 < 0 ) {
        ret = ret2;
        goto exit_function;
    }
    /* (+) */
    pEdgePlusSuper = NULL;
    pvPlusSuper    = NULL;
    if ( (k          = pTCGroups->nGroup[TCG_Plus]) >= 0 &&
         (ePlusSuper = pTCGroups->pTCG[k].nForwardEdge) > 0 &&
         (vPlusSuper = pTCGroups->pTCG[k].nVertexNumber) >= num_at ) {
         pEdgePlusSuper = pBNS->edge + ePlusSuper;
         pvPlusSuper    = pBNS->vert + vPlusSuper;
    }
    pEdgePlusCarbons = NULL;
    pvPlusCarbons    = NULL;
    if ( (k            = pTCGroups->nGroup[TCG_Plus_C0] )  > 0 &&
         (ePlusCarbons = pTCGroups->pTCG[k].nForwardEdge)  > 0 &&
         (vPlusCarbons = pTCGroups->pTCG[k].nVertexNumber) >= num_at ) {
        pEdgePlusCarbons = pBNS->edge + ePlusCarbons;
        pvPlusCarbons    = pBNS->vert + vPlusCarbons;
    }
    pEdgePlusHeteroat = NULL;
    pvPlusHeteroat    = NULL;
    if ( (k             = pTCGroups->nGroup[TCG_Plus0] )  > 0 &&
         (ePlusHeteroat = pTCGroups->pTCG[k].nForwardEdge)  > 0 &&
         (vPlusHeteroat = pTCGroups->pTCG[k].nVertexNumber) >= num_at ) {
        pEdgePlusHeteroat = pBNS->edge + ePlusHeteroat;
        pvPlusHeteroat     = pBNS->vert + vPlusHeteroat;
    }
    pEdgePlusMetals = NULL;
    pvPlusMetals    = NULL;
    if ( (k             = pTCGroups->nGroup[TCG_Plus_M0] )  > 0 &&
         (ePlusMetals = pTCGroups->pTCG[k].nForwardEdge)  > 0 &&
         (vPlusMetals = pTCGroups->pTCG[k].nVertexNumber) >= num_at ) {
        pEdgePlusMetals = pBNS->edge + ePlusMetals;
        pvPlusMetals    = pBNS->vert + vPlusMetals;
    }
    /* (-) */
    pEdgeMinusSuper = NULL;
    pvMinusSuper    = NULL;
    if ( (k          = pTCGroups->nGroup[TCG_Minus]) >= 0 &&
         (eMinusSuper = pTCGroups->pTCG[k].nForwardEdge) > 0 &&
         (vMinusSuper = pTCGroups->pTCG[k].nVertexNumber) >= num_at ) {
         pEdgeMinusSuper = pBNS->edge + eMinusSuper;
         pvMinusSuper    = pBNS->vert + vMinusSuper;
    }
    pEdgeMinusCarbons = NULL;
    pvMinusCarbons    = NULL;
    if ( (k            = pTCGroups->nGroup[TCG_Minus_C0] )  > 0 &&
         (eMinusCarbons = pTCGroups->pTCG[k].nForwardEdge)  > 0 &&
         (vMinusCarbons = pTCGroups->pTCG[k].nVertexNumber) >= num_at ) {
        pEdgeMinusCarbons = pBNS->edge + eMinusCarbons;
        pvMinusCarbons    = pBNS->vert + vMinusCarbons;
    }
    pEdgeMinusHeteroat = NULL;
    pvMinusHeteroat    = NULL;
    pEdgeMinusHeteroToSuper = NULL;
    if ( (k             = pTCGroups->nGroup[TCG_Minus0] )  > 0 &&
         (eMinusHeteroat = pTCGroups->pTCG[k].nForwardEdge)  > 0 &&
         (vMinusHeteroat = pTCGroups->pTCG[k].nVertexNumber) >= num_at ) {
        BNS_VERTEX *pvYMinusHetero;
        BNS_EDGE   *pe;
        int         vYMinusHetero;
        pEdgeMinusHeteroat = pBNS->edge + eMinusHeteroat;
        pvMinusHeteroat    = pBNS->vert + vMinusHeteroat;
        /* next edge toward (-)super */
        if ( pvMinusSuper ) {
            vYMinusHetero = pEdgeMinusHeteroat->neighbor12 ^ vMinusHeteroat;
            pvYMinusHetero = pBNS->vert + vYMinusHetero;
            for ( i = 0; i < pvYMinusHetero->num_adj_edges; i ++ ) {
                pe = pBNS->edge + pvYMinusHetero->iedge[i];
                if ( (pe->neighbor12 ^ vYMinusHetero) == vMinusSuper ) {
                    pEdgeMinusHeteroToSuper = pe;
                    eMinusHeteroToSuper     = pe - pBNS->edge;
                    break;
                }
            }
        }
    }
    pEdgeMinusMetals = NULL;
    pvMinusMetals    = NULL;
    if ( (k             = pTCGroups->nGroup[TCG_Minus_M0] )  > 0 &&
         (eMinusMetals = pTCGroups->pTCG[k].nForwardEdge)  > 0 &&
         (vMinusMetals = pTCGroups->pTCG[k].nVertexNumber) >= num_at ) {
        pEdgeMinusMetals = pBNS->edge + eMinusMetals;
        pvMinusMetals    = pBNS->vert + vMinusMetals;
    }
    /* (+/-) */
    pvPlMn = NULL;
    if ( pEdgePlusSuper ) {
        vPlMn = pEdgePlusSuper->neighbor12 ^ vPlusSuper;
        pvPlMn = pBNS->vert + vPlMn;
    } else
    if ( pEdgeMinusSuper ) {
        vPlMn = pEdgeMinusSuper->neighbor12 ^ vMinusSuper;
        pvPlMn = pBNS->vert + vPlMn;
    }
    num_pos = num_neg = 0;
    /***************************************************************/
    /* Positive Charges                                            */
    /***************************************************************/
    if ( pEdgePlusHeteroat && pEdgePlusMetals ) {
        /* count charges */
        for ( i = 0; i < num_at; i ++ ) {
            if ( !at2[i].radical &&
                  at2[i].charge > 0 && 
                  (k = pVA[i].nCPlusGroupEdge-1) >= 0 ) {
                v1 = pBNS->edge[k].neighbor1;
                v2 = pBNS->edge[k].neighbor1 ^ pBNS->edge[k].neighbor12;
                if ( v1 == vPlusHeteroat || v2 == vPlusHeteroat ) {
                    num_pos ++;
                }
            }
        }
        /* attempt to move (+) from heteroatoms to metal atoms */
        num_min = inchi_min( num_pos, pEdgePlusHeteroat->flow );

        nDeltaPlus1  = pvPlusSuper->st_edge.cap - pvPlusSuper->st_edge.flow;
        nDeltaMinus1 = pvPlMn->st_edge.cap - pvPlMn->st_edge.flow;
        if ( num_min && !nDeltaPlus1 && !nDeltaMinus1 ) {
            if ( pEdgePlusSuper ) {
                pEdgePlusSuper->forbidden |= forbidden_edge_mask;
            }
            if ( pEdgeMinusSuper ) {
                pEdgeMinusSuper->forbidden |= forbidden_edge_mask;
            }
            if ( pEdgePlusCarbons ) {
                pEdgePlusCarbons->forbidden |= forbidden_edge_mask;
            }
            if ( pEdgeMinusCarbons ) {
                pEdgeMinusCarbons->forbidden |= forbidden_edge_mask;
            }
            if ( pEdgePlusHeteroat ) {
                pEdgePlusHeteroat->forbidden |= forbidden_edge_mask;
            }
            if ( pEdgeMinusHeteroat ) {
                pEdgeMinusHeteroat->forbidden |= forbidden_edge_mask;
            }
            delta = 1;
            for ( i = 0; i < num_min; i += delta ) {
                v1 = pEdgePlusHeteroat->neighbor1;
                v2 = pEdgePlusHeteroat->neighbor12 ^ v1;
                pEdgePlusHeteroat->flow     -= delta;
                pBNS->vert[v1].st_edge.flow -= delta;
                pBNS->vert[v2].st_edge.flow -= delta;
                pBNS->tot_st_flow           -= 2*delta;
                /* test for charhe cancellation */
                ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                      &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );
                if ( ret < 0 ) {
                    goto exit_function;
                }
                if ( ret == 1 && (vPathEnd == v1 && vPathStart == v2 ||
                                  vPathEnd == v2 && vPathStart == v1) && nDeltaCharge == 0 ) {
                    ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                    (*pnNumRunBNS) ++;
                    if ( ret < 0 ) {
                        goto exit_function;
                    } else
                    if ( ret == 1 ) {
                        *pnTotalDelta += ret;
                    } else {
                        ret = RI_ERR_PROGR;
                        goto exit_function;
                    }
                } else {
                    pEdgePlusHeteroat->flow     += delta;
                    pBNS->vert[v1].st_edge.flow += delta;
                    pBNS->vert[v2].st_edge.flow += delta;
                    pBNS->tot_st_flow           += 2*delta;
                    break;
                }
            }
            if ( pEdgePlusSuper ) {
                pEdgePlusSuper->forbidden &= inv_forbidden_edge_mask;
            }
            if ( pEdgeMinusSuper ) {
                pEdgeMinusSuper->forbidden &= inv_forbidden_edge_mask;
            }
            if ( pEdgePlusCarbons ) {
                pEdgePlusCarbons->forbidden &= inv_forbidden_edge_mask;
            }
            if ( pEdgeMinusCarbons ) {
                pEdgeMinusCarbons->forbidden &= inv_forbidden_edge_mask;
            }
            if ( pEdgePlusHeteroat ) {
                pEdgePlusHeteroat->forbidden &= inv_forbidden_edge_mask;
            }
            if ( pEdgeMinusHeteroat ) {
                pEdgeMinusHeteroat->forbidden &= inv_forbidden_edge_mask;
            }
        }
    }
    /***************************************************************/
    /* Negative Charges                                            */
    /***************************************************************/
    if ( pEdgeMinusHeteroToSuper && pEdgeMinusMetals ) {
        /* count charges */
        for ( i = 0; i < num_at; i ++ ) {
            if ( !at2[i].radical &&
                  at2[i].charge < 0 && 
                  (k = pVA[i].nCMinusGroupEdge-1) >= 0 ) {
                v1 = pBNS->edge[k].neighbor1;
                v2 = pBNS->edge[k].neighbor1 ^ pBNS->edge[k].neighbor12;
                if ( v1 == vMinusHeteroat || v2 == vMinusHeteroat ) {
                    num_neg ++;
                }
            }
        }
        if ( num_neg ) {
            /* attempt to move (+) from heteroatoms to metal atoms */
            num_min = inchi_min( num_neg, pEdgeMinusHeteroToSuper->flow );
        }

        nDeltaPlus1  = pvPlusSuper->st_edge.cap - pvPlusSuper->st_edge.flow;
        nDeltaMinus1 = pvPlMn->st_edge.cap - pvPlMn->st_edge.flow;
        if ( num_min && !nDeltaPlus1 && !nDeltaMinus1 ) {
            if ( pEdgePlusSuper ) {
                pEdgePlusSuper->forbidden |= forbidden_edge_mask;
            }
            if ( pEdgeMinusSuper ) {
                pEdgeMinusSuper->forbidden |= forbidden_edge_mask;
            }
            if ( pEdgePlusCarbons ) {
                pEdgePlusCarbons->forbidden |= forbidden_edge_mask;
            }
            if ( pEdgeMinusCarbons ) {
                pEdgeMinusCarbons->forbidden |= forbidden_edge_mask;
            }
            if ( pEdgePlusHeteroat ) {
                pEdgePlusHeteroat->forbidden |= forbidden_edge_mask;
            }
            if ( pEdgeMinusHeteroToSuper ) {
                pEdgeMinusHeteroToSuper->forbidden |= forbidden_edge_mask;
            }
            delta = 1;
            for ( i = 0; i < num_min; i += delta ) {
                v1 = pEdgeMinusHeteroToSuper->neighbor1;
                v2 = pEdgeMinusHeteroToSuper->neighbor12 ^ v1;
                pEdgeMinusHeteroToSuper->flow -= delta;
                pBNS->vert[v1].st_edge.flow   -= delta;
                pBNS->vert[v2].st_edge.flow   -= delta;
                pBNS->tot_st_flow             -= 2*delta;
                /* test for charhe cancellation */
                ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                      &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );
                if ( ret < 0 ) {
                    goto exit_function;
                }
                if ( ret == 1 && (vPathEnd == v1 && vPathStart == v2 ||
                                  vPathEnd == v2 && vPathStart == v1) && nDeltaCharge == 0 ) {
                    ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                    (*pnNumRunBNS) ++;
                    if ( ret < 0 ) {
                        goto exit_function;
                    } else
                    if ( ret == 1 ) {
                        *pnTotalDelta += ret;
                    } else {
                        ret = RI_ERR_PROGR;
                        goto exit_function;
                    }
                } else {
                    pEdgeMinusHeteroToSuper->flow += delta;
                    pBNS->vert[v1].st_edge.flow   += delta;
                    pBNS->vert[v2].st_edge.flow   += delta;
                    pBNS->tot_st_flow             += 2*delta;
                    break;
                }
            }
            if ( pEdgePlusSuper ) {
                pEdgePlusSuper->forbidden &= inv_forbidden_edge_mask;
            }
            if ( pEdgeMinusSuper ) {
                pEdgeMinusSuper->forbidden &= inv_forbidden_edge_mask;
            }
            if ( pEdgePlusCarbons ) {
                pEdgePlusCarbons->forbidden &= inv_forbidden_edge_mask;
            }
            if ( pEdgeMinusCarbons ) {
                pEdgeMinusCarbons->forbidden &= inv_forbidden_edge_mask;
            }
            if ( pEdgePlusHeteroat ) {
                pEdgePlusHeteroat->forbidden &= inv_forbidden_edge_mask;
            }
            if ( pEdgeMinusHeteroToSuper ) {
                pEdgeMinusHeteroToSuper->forbidden &= inv_forbidden_edge_mask;
            }
        }
    }
exit_function:
    return ret;
}
#endif
/******************************************************************************************************/
int RestoreCyanoGroup( BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct,
                     inp_ATOM *at, inp_ATOM *at2, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups,
                     int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask)
{
    Vertex     vPathStart, vPathEnd;
    int        nPathLen, nDeltaH, nDeltaCharge, nNumVisitedAtoms;
    BNS_EDGE  *pe;
    
    int i, j, ret2, ret;
    int num_at = pStruct->num_atoms;
    int num_deleted_H = pStruct->num_deleted_H;
    int len_at = num_at + num_deleted_H;
    int inv_forbidden_edge_mask = ~forbidden_edge_mask;
    Vertex v1, v2;

    EDGE_LIST CarbonChargeEdges;
    
    ret = 0;
    AllocEdgeList( &CarbonChargeEdges, EDGE_LIST_CLEAR );

    memcpy( at2, at, len_at*sizeof(at2[0]));
    pStruct->at = at2;
    ret2 = CopyBnsToAtom( pStruct, pBNS, pVA, pTCGroups, 1 );
    pStruct->at = at;
    if ( ret2 < 0 ) {
        ret = ret2;
        goto exit_function;
    }

    for ( i = 0; i < num_at && 0 <= ret; i ++ ) {
        if ( at2[i].valence == 1 &&
             at2[i].num_H   == 0 &&
             at2[i].chem_bonds_valence == 2 &&
             at2[i].charge == -1 &&
             at2[i].radical == 0 &&
             pVA[i].cNumValenceElectrons == 5 &&  /* terminal N(-)=, P, As, Sb, Bi */
             pVA[i].nCMinusGroupEdge > 0 &&
             pVA[i].nTautGroupEdge == 0 &&
             at2[j=at2[i].neighbor[0]].valence == 2 &&
             at2[j].num_H == 0 &&
             at2[j].chem_bonds_valence == 4 &&
             at2[j].charge == 0 &&
             at2[j].radical == 0 &&
             pVA[j].cNumValenceElectrons == 4 && /* C or Si or Ge or Sn or Pb */
             pVA[i].cnListIndex > 0 &&
             cnList[pVA[i].cnListIndex-1].bits == cn_bits_MN ) {
            /* found N(-)=C= */
            pe = pBNS->edge + (pVA[i].nCMinusGroupEdge-1); /* N#N(+) triple bond edge */
            
            if ( !pe->flow ) {
                continue; /* wrong atom ??? Strange... */
            }
            v1 = pe->neighbor1;
            v2 = pe->neighbor12 ^ v1;
            pe->flow --;
            pBNS->vert[v1].st_edge.flow --;
            pBNS->vert[v2].st_edge.flow --;
            pBNS->tot_st_flow -= 2;
            pe->forbidden |= forbidden_edge_mask;

            /* do not let carbon atoms get charged */
            if ( 0 > (ret = ForbidCarbonChargeEdges( pBNS, pTCGroups, &CarbonChargeEdges, forbidden_edge_mask ))) {
                goto exit_function;
            }
            
            ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                  &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

            if ( ret == 1 && ((vPathEnd == v1 && vPathStart == v2) ||
                              (vPathEnd == v2 && vPathStart == v1)) && nDeltaCharge == 1 ) {
                ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                (*pnNumRunBNS) ++;
                *pnTotalDelta += ret;
            } else {
                pe->flow ++;
                pBNS->vert[v1].st_edge.flow ++;
                pBNS->vert[v2].st_edge.flow ++;
                pBNS->tot_st_flow += 2;
            }
            RemoveForbiddenEdgeMask( pBNS, &CarbonChargeEdges, forbidden_edge_mask );
            
            pe->forbidden &= inv_forbidden_edge_mask; /* unmask the edges */
        }
    }

exit_function:

    AllocEdgeList( &CarbonChargeEdges, EDGE_LIST_FREE );
    return ret;
}
/******************************************************************************************************/
int RestoreIsoCyanoGroup( BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct,
                     inp_ATOM *at, inp_ATOM *at2, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups,
                     int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask)
{
#define INC_EDGE_LIST 16
    Vertex     vPathStart, vPathEnd;
    int        nPathLen, nDeltaH, nDeltaCharge, nNumVisitedAtoms, num_failed, num_success;
    BNS_EDGE  *pe;
    Vertex    v1, v2;
    
    int i, j, ret2, ret, bIsCarbon;
    int num_at = pStruct->num_atoms;
    int num_deleted_H = pStruct->num_deleted_H;
    int len_at = num_at + num_deleted_H;
    int inv_forbidden_edge_mask = ~forbidden_edge_mask;
    EdgeIndex eNMinusEdge, eNPlusEdge, eNPlusEdge1, eN34Edge;
    EdgeIndex eNFlowerEdge1;

    EDGE_LIST CarbonChargeEdges, AllChargeEdges, IsoCyanoCarbonChargeEdges;

    ret = 0;
    num_failed = num_success = 0;
    AllocEdgeList( &CarbonChargeEdges, EDGE_LIST_CLEAR ); /* carbon charge edges */
    AllocEdgeList( &AllChargeEdges, EDGE_LIST_CLEAR );    /* heteroatom charge edges */
    AllocEdgeList( &IsoCyanoCarbonChargeEdges, EDGE_LIST_CLEAR );   /* C in C(+)#N(+) charge edges */

    memcpy( at2, at, len_at*sizeof(at2[0]));
    pStruct->at = at2;
    ret2 = CopyBnsToAtom( pStruct, pBNS, pVA, pTCGroups, 1 );
    pStruct->at = at;
    if ( ret2 < 0 ) {
        ret = ret2;
        goto exit_function;
    }
    /* 1st attempt: take care of C(+)#N(+)-  => C(-)#N(+)- and remove 2 negative charges */
    /* This would produce nDeltaCharge = 2 */
    AllocEdgeList( &CarbonChargeEdges, 2*num_at );
    for ( i = 0; i < num_at && 0 <= ret; i ++ ) {
        /* accumulate edges for subsequent fixing them */
        bIsCarbon = (pVA[i].cNumValenceElectrons == 4 && pVA[i].cPeriodicRowNumber == 1);
        eNFlowerEdge1 = NO_VERTEX;
        if ( (eNMinusEdge = pVA[i].nCMinusGroupEdge - 1)>= 0 &&
             !pBNS->edge[eNMinusEdge].forbidden ) {
            if ( bIsCarbon ) {
                if ( (ret = AddToEdgeList( &CarbonChargeEdges, eNMinusEdge, INC_EDGE_LIST )) ) {
                    goto exit_function;
                }
            } else
            if ( !pVA[i].cMetal && !at2[i].endpoint && at2[i].charge != -1 ) {
                if ( (ret = AddToEdgeList( &AllChargeEdges, eNMinusEdge, INC_EDGE_LIST )) ) {
                    goto exit_function;
                }
            }
        }
        if ( (eNPlusEdge = pVA[i].nCPlusGroupEdge - 1)>= 0 &&
             !pBNS->edge[eNPlusEdge].forbidden ) {
            if ( bIsCarbon ) {
                if ( (ret = AddToEdgeList( &CarbonChargeEdges, eNPlusEdge, INC_EDGE_LIST )) ) {
                    goto exit_function;
                }
            } else
            if ( !pVA[i].cMetal && !at2[i].endpoint  ) {
                if ( (ret = AddToEdgeList( &AllChargeEdges, eNPlusEdge, INC_EDGE_LIST )) ) {
                    goto exit_function;
                }
                if ( pVA[i].cNumValenceElectrons == 5 &&
                     NO_VERTEX != (eNFlowerEdge1 = GetChargeFlowerUpperEdge( pBNS, pVA, eNPlusEdge )) &&
                     !pBNS->edge[eNFlowerEdge1].flow ) {
                    if ( (ret = AddToEdgeList( &AllChargeEdges, eNFlowerEdge1, INC_EDGE_LIST )) ) {
                        goto exit_function;
                    }
                }
            }
        }
        if ( bIsCarbon &&
             0 <= eNMinusEdge &&
             0 <= eNPlusEdge &&
             at2[i].valence == 1 &&
             at2[i].num_H   == 0 &&
             at2[i].radical   == 0 &&
             !pBNS->edge[eNMinusEdge].forbidden &&
             pBNS->edge[eNMinusEdge].flow == 0 &&
             !pBNS->edge[eNPlusEdge].forbidden &&
             pBNS->edge[eNPlusEdge].flow == 0 &&   /* found terminal C(+) */
             
             at2[j=at2[i].neighbor[0]].valence == 2 &&
             at2[j].num_H == 0 &&
             at2[j].radical == 0 &&
             pVA[j].cNumValenceElectrons == 5 && 
             (eNPlusEdge1  = pVA[j].nCPlusGroupEdge - 1)>= 0 &&
             pBNS->edge[eNPlusEdge].flow == 0  ) {     /* -N(+)- */

#ifdef NEVER /* I have not found a good reason to do this yet */
            /* fix (+) charge on -N(+)- as much as C charges are fixed */
            if ( ret = AddToEdgeList( &CarbonChargeEdges, eNPlusEdge1, INC_EDGE_LIST ) ) {
                goto exit_function;
            }
            /* fix floer edge to prevent N(V) ??? */
            if ( NO_VERTEX != (eNFlowerEdge1 = GetChargeFlowerUpperEdge( pBNS, pVA, eNPlusEdge1 )) &&
                 !pBNS->edge[eNFlowerEdge1].flow ) {
                if ( ret = AddToEdgeList( &CarbonChargeEdges, eNFlowerEdge1, INC_EDGE_LIST ) ) {
                    goto exit_function;
                }
            }
#endif
            /*
               Carbon(+)         Carbon(-) 
               ChargeStruct:     ChargeStruct:
                                 
                         5(+C)             5(+C)
                        /                 //    
               6(-C)  4          6(-C)  4
                   \ //             \\ /        
                    3                 3         
                    |                 |         
                    2                 2         
                    ||                || 
                   -C1-              -C1-       
                    |                 |

             3-6 is (-) Charge Edge; 4-5 is (+) Charge Edge

             To convert the left pattern to the right one:

             We need to release these charge edges and decrement
             edge 3-4 flow to change charge from (+) to (-)
            
            */

            /* find vertices 4 and 5 */
            v1 = pBNS->edge[eNPlusEdge].neighbor1; /* one of two vertices incident with edge 4-5 */
            v2 = pBNS->edge[eNPlusEdge].neighbor12 ^ v1;
            if ( IS_BNS_VT_C_GR(pBNS->vert[v1].type) ) {
                /* v1 is 5(+C) */
                Vertex tmp = v1;
                v1 = v2;
                v2 = tmp;
            }
            /* v1 should be 4, v2 - 5(+C) */
            if ( !IS_BNS_VT_CHRG_STRUCT(pBNS->vert[v1].type) || pBNS->vert[v1].num_adj_edges != 2 ) {
                continue; /* mismatch */
            }
            /* find edge 3-4 */
            eN34Edge = pBNS->vert[v1].iedge[pBNS->vert[v1].iedge[0] == eNPlusEdge];
            if ( pBNS->edge[eN34Edge].forbidden || !pBNS->edge[eN34Edge].flow ) {
                continue;
            }
            /* save 3 edges: 6-3, 4-5, and 3-4 in this order */
            if ( (ret = AddToEdgeList( &IsoCyanoCarbonChargeEdges, eNMinusEdge, INC_EDGE_LIST )) ) {
                goto exit_function;
            }
            if ( (ret = AddToEdgeList( &IsoCyanoCarbonChargeEdges, eNPlusEdge, INC_EDGE_LIST )) ) {
                goto exit_function;
            }
            if ( (ret = AddToEdgeList( &IsoCyanoCarbonChargeEdges, eN34Edge, INC_EDGE_LIST )) ) {
                goto exit_function;
            }
        }
    }
    /* 1st attempt: move (-) charges from heteroatoms to C(+) */
    SetForbiddenEdgeMask( pBNS, &CarbonChargeEdges, forbidden_edge_mask );
    SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
    RemoveForbiddenEdgeMask( pBNS, &IsoCyanoCarbonChargeEdges, forbidden_edge_mask );
    for ( i = IsoCyanoCarbonChargeEdges.num_edges-3; 0 <= i; i -= 3 ) {
        eNMinusEdge = IsoCyanoCarbonChargeEdges.pnEdges[i];
        eNPlusEdge  = IsoCyanoCarbonChargeEdges.pnEdges[i+1];
        eN34Edge    = IsoCyanoCarbonChargeEdges.pnEdges[i+2];
        
        pe = pBNS->edge + eN34Edge;
        pe->forbidden |= forbidden_edge_mask;
        if ( !pe->flow ) {
            continue; /* already done */
        }

        v1 = pe->neighbor1;
        v2 = pe->neighbor12 ^ v1;
        pe->flow --;
        pBNS->vert[v1].st_edge.flow --;
        pBNS->vert[v2].st_edge.flow --;
        pBNS->tot_st_flow -= 2;
        
        ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                              &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

        if ( ret == 1 && ((vPathEnd == v1 && vPathStart == v2) ||
                          (vPathEnd == v2 && vPathStart == v1)) && nDeltaCharge <= -2 ) {
            ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
            (*pnNumRunBNS) ++;
            *pnTotalDelta += ret;
            num_success ++;
        } else {
            pe->flow ++;
            pBNS->vert[v1].st_edge.flow ++;
            pBNS->vert[v2].st_edge.flow ++;
            pBNS->tot_st_flow += 2;
            pe->forbidden &= inv_forbidden_edge_mask;
            num_failed ++;
        }
    }
    if ( num_failed ) {
        /* relax conditions: allow all heteroatoms to change charge */
        RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
        for ( i = IsoCyanoCarbonChargeEdges.num_edges-3; 0 <= i; i -= 3 ) {
            eNMinusEdge = IsoCyanoCarbonChargeEdges.pnEdges[i];
            eNPlusEdge  = IsoCyanoCarbonChargeEdges.pnEdges[i+1];
            eN34Edge    = IsoCyanoCarbonChargeEdges.pnEdges[i+2];
            
            pe = pBNS->edge + eN34Edge;
            pe->forbidden |= forbidden_edge_mask;
            if ( !pe->flow ) {
                continue; /* already done */
            }

            v1 = pe->neighbor1;
            v2 = pe->neighbor12 ^ v1;
            pe->flow --;
            pBNS->vert[v1].st_edge.flow --;
            pBNS->vert[v2].st_edge.flow --;
            pBNS->tot_st_flow -= 2;
            
            ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                  &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

            if ( ret == 1 && ((vPathEnd == v1 && vPathStart == v2) ||
                              (vPathEnd == v2 && vPathStart == v1)) && nDeltaCharge <= 2 ) {
                ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                (*pnNumRunBNS) ++;
                *pnTotalDelta += ret;
                num_success ++;
            } else {
                pe->flow ++;
                pBNS->vert[v1].st_edge.flow ++;
                pBNS->vert[v2].st_edge.flow ++;
                pBNS->tot_st_flow += 2;
                pe->forbidden &= inv_forbidden_edge_mask; /* let it change if it wants */
                num_failed ++;
            }
        }
    }
    RemoveForbiddenEdgeMask( pBNS, &CarbonChargeEdges, forbidden_edge_mask );
    RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
    RemoveForbiddenEdgeMask( pBNS, &IsoCyanoCarbonChargeEdges, forbidden_edge_mask );


exit_function:

    AllocEdgeList( &CarbonChargeEdges, EDGE_LIST_FREE );
    AllocEdgeList( &AllChargeEdges, EDGE_LIST_FREE );
    AllocEdgeList( &IsoCyanoCarbonChargeEdges, EDGE_LIST_FREE );
    return ret;
#undef INC_EDGE_LIST
}

/******************************************************************************************************/
int FixMetal_Nminus_Ominus( BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct,
                     inp_ATOM *at, inp_ATOM *at2, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups,
                     int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask)
{
#define INC_EDGE_LIST 16
    Vertex     vPathStart, vPathEnd;
    int        nPathLen, nDeltaH, nDeltaCharge, nNumVisitedAtoms;
    int        num_failed, num_success, n, nDeltaChargeMax, nMetalCharge;
    BNS_EDGE  *pe;
    Vertex    v1, v2;
    
    int i, j, k, ret2, ret;
    int num_at = pStruct->num_atoms;
    int num_deleted_H = pStruct->num_deleted_H;
    int len_at = num_at + num_deleted_H;
    int inv_forbidden_edge_mask = ~forbidden_edge_mask;
    EdgeIndex e, eNMinusEdge, eNMinusEdge1, eNMinusEdge2, eNPlusEdge2;

    EDGE_LIST AllChargeEdges;

    ret = 0;
    num_failed = num_success = 0;
    AllocEdgeList( &AllChargeEdges, EDGE_LIST_CLEAR );

    memcpy( at2, at, len_at*sizeof(at2[0]));
    pStruct->at = at2;
    ret2 = CopyBnsToAtom( pStruct, pBNS, pVA, pTCGroups, 1 );
    pStruct->at = at;
    if ( ret2 < 0 ) {
        ret = ret2;
        goto exit_function;
    }
    /* attepmt #1 N#N(+)-N => N(-)=N(+)=N */
    for ( i = 0; i < num_at && 0 <= ret; i ++ ) {
        if ( at2[i].valence == 1 &&
             at2[i].num_H   == 0 &&
             at2[i].radical == 0 &&
             pVA[i].cNumValenceElectrons == 6 &&  /* terminal -O */
             (eNMinusEdge = pVA[i].nCMinusGroupEdge - 1)>= 0 && pBNS->edge[eNMinusEdge].flow == 1 &&
             !pBNS->edge[eNMinusEdge].forbidden && /* terminal O(-) */
             
             at2[j=at2[i].neighbor[0]].valence == 2 &&
             at2[j].num_H == 0 &&
             at2[j].radical == 0 &&
             pVA[j].cNumValenceElectrons == 5 &&
             (eNMinusEdge1 = pVA[j].nCMinusGroupEdge - 1)>= 0 && pBNS->edge[eNMinusEdge1].flow == 1 &&
             !pBNS->edge[eNMinusEdge1].forbidden &&
             
             pVA[k=at2[j].neighbor[at2[j].neighbor[0]==i]].cMetal &&

             (eNMinusEdge2 = pVA[k].nCMinusGroupEdge - 1)>= 0 &&
             !pBNS->edge[eNMinusEdge2].forbidden &&
             (eNPlusEdge2 = pVA[k].nCPlusGroupEdge - 1)>= 0 &&
             !pBNS->edge[eNPlusEdge2].forbidden  ) {

            /* found M(q)-N(-)-O(-); convert to M(q-2)-N=O */

            /* find all charge edges to fix */
            if ( 0 == AllChargeEdges.num_edges ) {
                for ( n = 0; n < num_at; n ++ ) {
                    if ( (e = pVA[n].nCMinusGroupEdge - 1)>= 0 &&
                         !pBNS->edge[e].forbidden ) {
                        if ( (ret = AddToEdgeList( &AllChargeEdges, e, num_at )) ) {
                            goto exit_function;
                        }
                    }
                    if ( (e = pVA[n].nCPlusGroupEdge - 1)>= 0 &&
                         !pBNS->edge[e].forbidden ) {
                        if ( (ret = AddToEdgeList( &AllChargeEdges, e, num_at )) ) {
                            goto exit_function;
                        }
                        if ( pVA[n].cNumValenceElectrons == 6 &&
                             NO_VERTEX != (e = GetChargeFlowerUpperEdge( pBNS, pVA, e )) &&
                             pBNS->edge[e].flow == 0 ) {
                            if ( (ret = AddToEdgeList( &AllChargeEdges, e, num_at )) ) {
                                goto exit_function;
                            }
                        }
                    }
                }
            }

            nMetalCharge = (pBNS->edge[eNPlusEdge2].cap - pBNS->edge[eNPlusEdge2].flow)
                           - pBNS->edge[eNMinusEdge2].flow;
            if ( nMetalCharge == 0 ) {
                /* change on O is invisible; charge from N(-) goes, charge comes to Metal */
                nDeltaChargeMax = 0;
            } else
            if ( nMetalCharge == 2 ) {
                /* charges on Metal and N disappear */
                nDeltaChargeMax = -2;
            } else {
                /* charge from N disappears */
                nDeltaChargeMax = -1;
            }

            SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
            pBNS->edge[eNMinusEdge1].forbidden &= inv_forbidden_edge_mask;
            pBNS->edge[eNMinusEdge2].forbidden &= inv_forbidden_edge_mask;
            pBNS->edge[eNPlusEdge2].forbidden &= inv_forbidden_edge_mask;

            pe = pBNS->edge + eNMinusEdge; /* must be already fixed as a charge edge */

            v1 = pe->neighbor1;
            v2 = pe->neighbor12 ^ v1;
            pe->flow --;
            pBNS->vert[v1].st_edge.flow --;
            pBNS->vert[v2].st_edge.flow --;
            pBNS->tot_st_flow -= 2;

            ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                  &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

            if ( ret == 1 && ((vPathEnd == v1 && vPathStart == v2) ||
                              (vPathEnd == v2 && vPathStart == v1)) /*&& nDeltaCharge == nDeltaChargeMax*/ ) {
                ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                (*pnNumRunBNS) ++;
                *pnTotalDelta += ret;
                num_success ++;
            } else {
                pe->flow ++;
                pBNS->vert[v1].st_edge.flow ++;
                pBNS->vert[v2].st_edge.flow ++;
                pBNS->tot_st_flow += 2;
                num_failed ++;
            }
            RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
        }
    }
    ret = num_success;

exit_function:

    AllocEdgeList( &AllChargeEdges, EDGE_LIST_FREE );

    return ret;
#undef INC_EDGE_LIST
}
/******************************************************************************************************/
int RestoreNNNgroup( BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct,
                     inp_ATOM *at, inp_ATOM *at2, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups,
                     int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask)
{
#define INC_EDGE_LIST 16
    Vertex     vPathStart, vPathEnd;
    int        nPathLen, nDeltaH, nDeltaCharge, nNumVisitedAtoms, num_failed, num_success, n, nDeltaChargeMax;
    BNS_EDGE  *pe;
    Vertex    v1, v2;
    
    int i, j, k, ret2, ret;
    int num_at = pStruct->num_atoms;
    int num_deleted_H = pStruct->num_deleted_H;
    int len_at = num_at + num_deleted_H;
    int inv_forbidden_edge_mask = ~forbidden_edge_mask;
    EdgeIndex eNMinusEdge, eNPlusEdge, eNMinusEdge1, eNPlusEdge1, eNMinusEdge2, eNPlusEdge2;
    EdgeIndex eNFlowerEdge1, eNFlowerEdge2;

    EDGE_LIST CarbonChargeEdges, AllChargeEdges, NNNChargeEdges, CurNNNChargeEdges, AllNNNTermAtoms, AllNIIIChargeEdges;

    ret = 0;
    num_failed = num_success = 0;
    AllocEdgeList( &CarbonChargeEdges, EDGE_LIST_CLEAR );
    AllocEdgeList( &AllChargeEdges, EDGE_LIST_CLEAR );
    AllocEdgeList( &NNNChargeEdges, EDGE_LIST_CLEAR );
    AllocEdgeList( &CurNNNChargeEdges, EDGE_LIST_CLEAR );
    AllocEdgeList( &AllNNNTermAtoms, EDGE_LIST_CLEAR );
    AllocEdgeList( &AllNIIIChargeEdges, EDGE_LIST_CLEAR );

    memcpy( at2, at, len_at*sizeof(at2[0]));
    pStruct->at = at2;
    ret2 = CopyBnsToAtom( pStruct, pBNS, pVA, pTCGroups, 1 );
    pStruct->at = at;
    if ( ret2 < 0 ) {
        ret = ret2;
        goto exit_function;
    }
    /* attepmt #1 N#N(+)-N => N(-)=N(+)=N: naive approach; expecting tp move (-) from some other atom */
    for ( i = 0; i < num_at && 0 <= ret; i ++ ) {
        if ( at2[i].valence == 1 &&
             at2[i].num_H   == 0 &&
             at2[i].chem_bonds_valence == 3 &&
             at2[i].charge ==  0 &&
             at2[i].radical == 0 &&
             pVA[i].cNumValenceElectrons == 5 &&  /* terminal N# */
             (eNMinusEdge = pVA[i].nCMinusGroupEdge - 1)>= 0 && pBNS->edge[eNMinusEdge].flow == 0 &&
             !pBNS->edge[eNMinusEdge].forbidden &&
             at2[j=at2[i].neighbor[0]].valence == 2 &&
             at2[j].num_H == 0 &&
             at2[j].chem_bonds_valence == 4 &&
             at2[j].charge == 1 &&
             at2[j].radical == 0 &&
             pVA[j].cNumValenceElectrons == 5 &&
             (eNPlusEdge = pVA[j].nCPlusGroupEdge - 1)>= 0 && pBNS->edge[eNPlusEdge].flow == 0 &&
             !pBNS->edge[eNPlusEdge].forbidden &&
             at2[k=at2[j].neighbor[at2[j].neighbor[0]==i]].valence == 2 &&
             at2[k].num_H == 0 &&
             at2[k].chem_bonds_valence == 3 &&
             pVA[k].cNumValenceElectrons == 5 &&
             (eNPlusEdge2 = pVA[k].nCPlusGroupEdge - 1)>= 0 && pBNS->edge[eNPlusEdge2].flow == 1 &&
             !pBNS->edge[eNPlusEdge2].forbidden &&
             pVA[i].cnListIndex > 0 &&
             cnList[pVA[i].cnListIndex-1].bits == cn_bits_MN ) {
            /* found N#N(+)-N~ where the last N (at2[k]) may be charged */
            pe = pBNS->edge + pBNS->vert[i].iedge[0]; /* N#N(+) triple bond edge */
            
            v1 = pe->neighbor1;
            v2 = pe->neighbor12 ^ v1;
            pe->flow --;
            pBNS->vert[v1].st_edge.flow --;
            pBNS->vert[v2].st_edge.flow --;
            pBNS->tot_st_flow -= 2;

            pe->forbidden                     |= forbidden_edge_mask;
            pBNS->edge[eNPlusEdge].forbidden  |= forbidden_edge_mask;
            pBNS->edge[eNPlusEdge2].forbidden |= forbidden_edge_mask;
            
            
            if ( !CarbonChargeEdges.num_edges ) {
                /* do not let carbon atoms get charged */
                AllocEdgeList( &CarbonChargeEdges, INC_EDGE_LIST ); 
                if ( 0 > (ret = ForbidCarbonChargeEdges( pBNS, pTCGroups, &CarbonChargeEdges, forbidden_edge_mask ))) {
                    goto exit_function;
                }
            } else {
                SetForbiddenEdgeMask( pBNS, &CarbonChargeEdges, forbidden_edge_mask );
            }
            ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                  &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

            if ( ret == 1 && ((vPathEnd == v1 && vPathStart == v2) ||
                              (vPathEnd == v2 && vPathStart == v1)) && nDeltaCharge <= 0 ) {
                ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                (*pnNumRunBNS) ++;
                *pnTotalDelta += ret;
                num_success ++;
                /* fix charges on N(-)=N(+)=N- */
                if ( (ret = AddToEdgeList( &CarbonChargeEdges, eNMinusEdge, INC_EDGE_LIST )) ) {
                    goto exit_function;
                }
                if ( (ret = AddToEdgeList( &CarbonChargeEdges, eNPlusEdge, INC_EDGE_LIST )) ) {
                    goto exit_function;
                }
            } else {
                pe->flow ++;
                pBNS->vert[v1].st_edge.flow ++;
                pBNS->vert[v2].st_edge.flow ++;
                pBNS->tot_st_flow += 2;
                num_failed ++;
            }
            RemoveForbiddenEdgeMask( pBNS, &CarbonChargeEdges, forbidden_edge_mask );

            pe->forbidden                     &= inv_forbidden_edge_mask;
            pBNS->edge[eNPlusEdge].forbidden  &= inv_forbidden_edge_mask;
            pBNS->edge[eNPlusEdge2].forbidden &= inv_forbidden_edge_mask;
        }
    }

    /* 2nd attempt: take care of N#N(+)-N=-...=N-N(-) */
    /* This would produce nDeltaCharge >= 2 */
    
    AllChargeEdges.num_edges = 0;
    AllNNNTermAtoms.num_edges = 0;
    NNNChargeEdges.num_edges = 0;
    AllNIIIChargeEdges.num_edges = 0;

    for ( i = 0; i < num_at && 0 <= ret; i ++ ) {
        if ( (eNMinusEdge = pVA[i].nCMinusGroupEdge - 1)>= 0 &&
             !pBNS->edge[eNMinusEdge].forbidden ) {
            if ( (ret = AddToEdgeList( &AllChargeEdges, eNMinusEdge, INC_EDGE_LIST )) ) {
                goto exit_function;
            }
        } else {
            eNMinusEdge = -1;
        }
        if ( (eNPlusEdge = pVA[i].nCPlusGroupEdge - 1)>= 0 &&
             !pBNS->edge[eNPlusEdge].forbidden ) {
            if ( (ret = AddToEdgeList( &AllChargeEdges, eNPlusEdge, INC_EDGE_LIST )) ) {
                goto exit_function;
            }
            if ( pVA[i].cNumValenceElectrons == 5 && at2[i].valence == 3 && at2[i].chem_bonds_valence == 3) {
                if ( (ret = AddToEdgeList( &AllNIIIChargeEdges, eNPlusEdge, INC_EDGE_LIST )) ) {
                    goto exit_function;
                }
            }
            /* N flower edge */
            if ( pVA[i].cNumValenceElectrons == 5 && pVA[i].cPeriodicRowNumber == 1 &&
                 NO_VERTEX != (eNFlowerEdge1 = GetChargeFlowerUpperEdge( pBNS, pVA, eNPlusEdge )) &&
                 pBNS->edge[eNFlowerEdge1].flow == 0 &&
                 ( ret = AddToEdgeList( &AllChargeEdges, eNFlowerEdge1, INC_EDGE_LIST ) ) ) {
                goto exit_function;
            }
        } else {
            eNPlusEdge = -1;
        }

        if ( 0 <= eNMinusEdge &&
             at2[i].valence == 1 &&
             at2[i].num_H   == 0 &&
             at2[i].radical   == 0 &&
             pVA[i].cNumValenceElectrons == 5 &&  /* terminal N# */
             
             at2[j=at2[i].neighbor[0]].valence == 2 &&
             at2[j].num_H == 0 &&
             at2[j].radical == 0 &&
             pVA[j].cNumValenceElectrons == 5 && 
             (eNMinusEdge1 = pVA[j].nCMinusGroupEdge - 1)>= 0 &&
             (eNPlusEdge1  = pVA[j].nCPlusGroupEdge - 1)>= 0 &&
             !pBNS->edge[eNMinusEdge1].forbidden &&
             !pBNS->edge[eNPlusEdge1].forbidden &&

             at2[k=at2[j].neighbor[at2[j].neighbor[0]==i]].valence == 2 &&
             at2[k].num_H == 0 &&
             at2[k].radical == 0 &&
             pVA[k].cNumValenceElectrons == 5 &&
             (eNMinusEdge2 = pVA[k].nCMinusGroupEdge - 1)>= 0 &&
             (eNPlusEdge2  = pVA[k].nCPlusGroupEdge - 1)>= 0 &&
             !pBNS->edge[eNMinusEdge2].forbidden &&
             !pBNS->edge[eNPlusEdge2].forbidden &&

             pVA[i].cnListIndex > 0 &&
             cnList[pVA[i].cnListIndex-1].bits == cn_bits_MN ) {
            /* found N#N(+)-N~ or N(-)=N-N= where the last N (at2[k]) may be charged */
            
            /* 1. N(-)=N(+)=N- */
            if ( pBNS->edge[eNMinusEdge].flow  == 1 &&  /* N(-) */
                 pBNS->edge[eNMinusEdge1].flow == 0 && pBNS->edge[eNPlusEdge1].flow == 0 && /* N(+) */
                 pBNS->edge[eNMinusEdge2].flow == 0 && pBNS->edge[eNPlusEdge2].flow == 1    /* N */ ) {
                continue; /* already good */
            }
            /* accumulate terminal atoms of all other NNN */
            if ( (ret = AddToEdgeList( &AllNNNTermAtoms, i, INC_EDGE_LIST )) ) {
                goto exit_function;
            }
            /* 2. N#N(+)-N= */
            if ( pBNS->edge[eNMinusEdge].flow  == 0 && /* N */
                 pBNS->edge[eNMinusEdge1].flow == 0 && pBNS->edge[eNPlusEdge1].flow == 0 && /* N(+) */
                 pBNS->edge[eNMinusEdge2].flow == 0 && pBNS->edge[eNPlusEdge2].flow == 1    /* N */ ) {
                /* unfix (-) edge on terminal N# */
                if ( (ret = AddToEdgeList( &NNNChargeEdges, eNMinusEdge, INC_EDGE_LIST )) ) {
                    goto exit_function;
                }
                continue;
            }
            /* 3. N(-)=N-N= */
            if ( pBNS->edge[eNMinusEdge].flow  == 1 && /* N(-) */
                 pBNS->edge[eNMinusEdge1].flow == 0 && pBNS->edge[eNPlusEdge1].flow == 1 && /* N */
                 pBNS->edge[eNMinusEdge2].flow == 0 && pBNS->edge[eNPlusEdge2].flow == 1    /* N */ ) {
                /* unfix (+) edge on middle N */
                if ( (ret = AddToEdgeList( &NNNChargeEdges, eNPlusEdge1, INC_EDGE_LIST )) ) {
                    goto exit_function;
                }
                continue;
            }
            /* 4. N#N(+)-N(-)- */
            if ( pBNS->edge[eNMinusEdge].flow  == 0 &&  /* N */
                 pBNS->edge[eNMinusEdge1].flow == 0 && pBNS->edge[eNPlusEdge1].flow == 0 && /* N(+) */
                 pBNS->edge[eNMinusEdge2].flow == 1 && pBNS->edge[eNPlusEdge2].flow == 1    /* N(-) */ ) {
                /* unfix (-) edge on the 1st and 3rd N */
                if ( (ret = AddToEdgeList( &NNNChargeEdges, eNMinusEdge, INC_EDGE_LIST )) ) {
                    goto exit_function;
                }
                if ( (ret = AddToEdgeList( &NNNChargeEdges, eNMinusEdge2, INC_EDGE_LIST )) ) {
                    goto exit_function;
                }
                continue;
            }
            /* 5. N#N(+)-N(+)# */
            if ( pBNS->edge[eNMinusEdge].flow  == 0 &&  /* N */
                 pBNS->edge[eNMinusEdge1].flow == 0 && pBNS->edge[eNPlusEdge1].flow == 0 && /* N(+) */
                 pBNS->edge[eNMinusEdge2].flow == 0 && pBNS->edge[eNPlusEdge2].flow == 0    /* N(+) */ ) {
                /* unfix (-) edge on the 1st and (+) edge on the 3rd N */
                if ( (ret = AddToEdgeList( &NNNChargeEdges, eNMinusEdge, INC_EDGE_LIST )) ) {
                    goto exit_function;
                }
                if ( (ret = AddToEdgeList( &NNNChargeEdges, eNPlusEdge2, INC_EDGE_LIST )) ) {
                    goto exit_function;
                }
                continue;
            }
        }
    }
    /* try to fix each NNN */
    for ( n = AllNNNTermAtoms.num_edges-1; 0 <= n; n -- ) {
        i = AllNNNTermAtoms.pnEdges[n];
        eNMinusEdge = pVA[i].nCMinusGroupEdge - 1;
        /*eNPlusEdge = pVA[i].nCPlusGroupEdge - 1;*/
        j=at2[i].neighbor[0];
        eNMinusEdge1 = pVA[j].nCMinusGroupEdge - 1;
        eNPlusEdge1  = pVA[j].nCPlusGroupEdge - 1;
        k=at2[j].neighbor[at2[j].neighbor[0]==i];
        eNMinusEdge2 = pVA[k].nCMinusGroupEdge - 1;
        eNPlusEdge2  = pVA[k].nCPlusGroupEdge - 1;
        /*SetForbiddenEdgeMask( pBNS, &CarbonChargeEdges, forbidden_edge_mask );*/
        /* 1. N(-)=N(+)=N- */
        if ( pBNS->edge[eNMinusEdge].flow  == 1 &&  /* N(-) */
             pBNS->edge[eNMinusEdge1].flow == 0 && pBNS->edge[eNPlusEdge1].flow == 0 && /* N(+) */
             pBNS->edge[eNMinusEdge2].flow == 0 && pBNS->edge[eNPlusEdge2].flow == 1    /* N */ ) {

            RemoveFromEdgeListByValue( &NNNChargeEdges, eNMinusEdge );
            RemoveFromEdgeListByValue( &NNNChargeEdges, eNMinusEdge1 );
            RemoveFromEdgeListByValue( &NNNChargeEdges, eNPlusEdge1 );
            RemoveFromEdgeListByValue( &NNNChargeEdges, eNMinusEdge2 );
            RemoveFromEdgeListByValue( &NNNChargeEdges, eNPlusEdge2 );

            pe = NULL;
        } else   /* 2. N#N(+)-N= */
        if ( pBNS->edge[eNMinusEdge].flow  == 0 &&  /* N */
             pBNS->edge[eNMinusEdge1].flow == 0 && pBNS->edge[eNPlusEdge1].flow == 0 && /* N(+) */
             pBNS->edge[eNMinusEdge2].flow == 0 && pBNS->edge[eNPlusEdge2].flow == 1    /* N */ ) {
            /* decrement triple bond on terminal N# */
            pe = pBNS->edge + pBNS->vert[i].iedge[0];
        } else
        /* 3. N(-)=N-N= */
        if ( pBNS->edge[eNMinusEdge].flow  == 1 &&  /* N(-) */
             pBNS->edge[eNMinusEdge1].flow == 0 && pBNS->edge[eNPlusEdge1].flow == 1 && /* N */
             pBNS->edge[eNMinusEdge2].flow == 0 && pBNS->edge[eNPlusEdge2].flow == 1    /* N */ ) {
            /* decrement flow on (+) charge edge of the middle =N- */
            pe = pBNS->edge + eNPlusEdge1;
        } else
        /* 4. N#N(+)-N(-)- */
        if ( pBNS->edge[eNMinusEdge].flow  == 0 &&  /* N */
             pBNS->edge[eNMinusEdge1].flow == 0 && pBNS->edge[eNPlusEdge1].flow == 0 && /* N(+) */
             pBNS->edge[eNMinusEdge2].flow == 1 && pBNS->edge[eNPlusEdge2].flow == 1    /* N(-) */ ) {
            /* decrement triple bond on terminal N# */
            pe = pBNS->edge + pBNS->vert[i].iedge[0];
        } else
        /* 5. N#N(+)-N(+)# */
        if ( pBNS->edge[eNMinusEdge].flow  == 0 &&  /* N */
             pBNS->edge[eNMinusEdge1].flow == 0 && pBNS->edge[eNPlusEdge1].flow == 0 && /* N(+) */
             pBNS->edge[eNMinusEdge2].flow == 0 && pBNS->edge[eNPlusEdge2].flow == 0    /* N(+) */ ) {
            /* decrement triple bond on terminal N# */
            pe = pBNS->edge + pBNS->vert[i].iedge[0];
        } else {
            pe = NULL; /* unknown case */
        }
        if ( pe ) {
            SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
            RemoveForbiddenEdgeMask( pBNS, &NNNChargeEdges, forbidden_edge_mask );
            
            v1 = pe->neighbor1;
            v2 = pe->neighbor12 ^ v1;
            pe->flow --;
            pBNS->vert[v1].st_edge.flow --;
            pBNS->vert[v2].st_edge.flow --;
            pBNS->tot_st_flow -= 2;
            
            pe->forbidden                    |= forbidden_edge_mask;
            
            ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                  &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

            if ( ret == 1 && ((vPathEnd == v1 && vPathStart == v2) ||
                              (vPathEnd == v2 && vPathStart == v1)) /*&& nDeltaCharge <= 2*/ ) {
                ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                (*pnNumRunBNS) ++;
                *pnTotalDelta += ret;
                num_success ++;
                /* fix charges on N(-)=N(+)=N- */
                RemoveFromEdgeListByValue( &NNNChargeEdges, eNMinusEdge );
                RemoveFromEdgeListByValue( &NNNChargeEdges, eNMinusEdge1 );
                RemoveFromEdgeListByValue( &NNNChargeEdges, eNPlusEdge1 );
                RemoveFromEdgeListByValue( &NNNChargeEdges, eNMinusEdge2 );
                RemoveFromEdgeListByValue( &NNNChargeEdges, eNPlusEdge2 );
            } else {
                pe->flow ++;
                pBNS->vert[v1].st_edge.flow ++;
                pBNS->vert[v2].st_edge.flow ++;
                pBNS->tot_st_flow += 2;
                num_failed ++;
            }
            pe->forbidden &= inv_forbidden_edge_mask;
            RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
        }
    }
    
    /* 3rd attempt */

    /*
    AllChargeEdges.num_edges = 0;
    AllNNNTermAtoms.num_edges = 0;
    NNNChargeEdges.num_edges = 0;
    */
    for ( i = 0; i < num_at && 0 <= ret; i ++ ) {

        eNMinusEdge = pVA[i].nCMinusGroupEdge - 1;
        /*eNPlusEdge = pVA[i].nCPlusGroupEdge - 1;*/

        if ( 0 <= eNMinusEdge &&
             at2[i].valence == 1 &&
             at2[i].num_H   == 0 &&
             at2[i].radical   == 0 &&
             pVA[i].cNumValenceElectrons == 5 &&  /* terminal N# */
             
             at2[j=at2[i].neighbor[0]].valence == 2 &&
             at2[j].num_H == 0 &&
             at2[j].radical == 0 &&
             pVA[j].cNumValenceElectrons == 5 && 
             (eNMinusEdge1 = pVA[j].nCMinusGroupEdge - 1)>= 0 &&
             (eNPlusEdge1  = pVA[j].nCPlusGroupEdge - 1)>= 0 &&
             !pBNS->edge[eNMinusEdge1].forbidden &&
             !pBNS->edge[eNPlusEdge1].forbidden &&

             at2[k=at2[j].neighbor[at2[j].neighbor[0]==i]].valence == 2 &&
             at2[k].num_H == 0 &&
             at2[k].radical == 0 &&
             pVA[k].cNumValenceElectrons == 5 &&
             (eNMinusEdge2 = pVA[k].nCMinusGroupEdge - 1)>= 0 &&
             (eNPlusEdge2  = pVA[k].nCPlusGroupEdge - 1)>= 0 &&
             !pBNS->edge[eNMinusEdge2].forbidden &&
             !pBNS->edge[eNPlusEdge2].forbidden &&

             pVA[i].cnListIndex > 0 &&
             cnList[pVA[i].cnListIndex-1].bits == cn_bits_MN ) {

            /* found N#N(+)-N~ or N(-)=N-N= where the last N (at2[k]) may be charged */
            NNNChargeEdges.num_edges = 0;

            eNFlowerEdge1 = GetChargeFlowerUpperEdge( pBNS, pVA, eNPlusEdge1 );
            eNFlowerEdge2 = GetChargeFlowerUpperEdge( pBNS, pVA, eNPlusEdge2 );

            /* 1. N(-)=N(+)=N- */
            if ( pBNS->edge[eNMinusEdge].flow  == 1 &&  /* N(-) */
                 pBNS->edge[eNMinusEdge1].flow == 0 && pBNS->edge[eNPlusEdge1].flow == 0 && /* N(+) */
                 pBNS->edge[eNMinusEdge2].flow == 0 && pBNS->edge[eNPlusEdge2].flow == 1    /* N */ ) {
                /* fix charges on N(-)=N(+)=N- */
                if ( (ret = AddToEdgeList( &CarbonChargeEdges, eNMinusEdge, INC_EDGE_LIST )) ) {
                    goto exit_function;
                }
                if ( (ret = AddToEdgeList( &CarbonChargeEdges, eNPlusEdge1, INC_EDGE_LIST )) ) {
                    goto exit_function;
                }
                if ( (ret = AddToEdgeList( &CarbonChargeEdges, eNMinusEdge1, INC_EDGE_LIST )) ) {
                    goto exit_function;
                }
                if ( (ret = AddToEdgeList( &CarbonChargeEdges, eNPlusEdge2, INC_EDGE_LIST )) ) {
                    goto exit_function;
                }
                if ( (ret = AddToEdgeList( &CarbonChargeEdges, eNMinusEdge2, INC_EDGE_LIST )) ) {
                    goto exit_function;
                }
                continue; /* already good */
            }
            /* 2. N#N(+)-N= */
            if ( pBNS->edge[eNMinusEdge].flow  == 0 && /* N */
                 pBNS->edge[eNMinusEdge1].flow == 0 && pBNS->edge[eNPlusEdge1].flow == 0 && /* N(+) */
                 pBNS->edge[eNMinusEdge2].flow == 0 && pBNS->edge[eNPlusEdge2].flow == 1    /* N */ ) {
                /* unfix (-) edge on terminal N# */
                if ( (ret = AddToEdgeList( &NNNChargeEdges, eNPlusEdge1, INC_EDGE_LIST )) ) {
                    goto exit_function;
                }
                if ( (ret = AddToEdgeList( &NNNChargeEdges, eNMinusEdge1, INC_EDGE_LIST )) ) {
                    goto exit_function;
                }
                if ( (ret = AddToEdgeList( &NNNChargeEdges, eNPlusEdge2, INC_EDGE_LIST )) ) {
                    goto exit_function;
                }
                if ( (ret = AddToEdgeList( &NNNChargeEdges, eNMinusEdge2, INC_EDGE_LIST )) ) {
                    goto exit_function;
                }
                pe = pBNS->edge + pBNS->vert[i].iedge[0];
                nDeltaChargeMax = 0;
                nDeltaChargeMax = (num_failed && !num_success && pStruct->nNumRemovedProtonsMobHInChI > 0)? 2 : 0;
            } else
            /* 3. N(-)=N-N= */
            if ( pBNS->edge[eNMinusEdge].flow  == 1 && /* N(-) */
                 pBNS->edge[eNMinusEdge1].flow == 0 && pBNS->edge[eNPlusEdge1].flow == 1 && /* N */
                 pBNS->edge[eNMinusEdge2].flow == 0 && pBNS->edge[eNPlusEdge2].flow == 1    /* N */ ) {
                /* unfix (+) edge on middle N */
                if ( (ret = AddToEdgeList( &NNNChargeEdges, eNMinusEdge, INC_EDGE_LIST )) ) {
                    goto exit_function;
                }
                if ( (ret = AddToEdgeList( &NNNChargeEdges, eNMinusEdge1, INC_EDGE_LIST )) ) {
                    goto exit_function;
                }
                if ( (ret = AddToEdgeList( &NNNChargeEdges, eNPlusEdge2, INC_EDGE_LIST )) ) {
                    goto exit_function;
                }
                if ( (ret = AddToEdgeList( &NNNChargeEdges, eNMinusEdge2, INC_EDGE_LIST )) ) {
                    goto exit_function;
                }
                if ( NO_VERTEX != eNFlowerEdge1 &&
                     ( ret = AddToEdgeList( &NNNChargeEdges, eNFlowerEdge1, INC_EDGE_LIST ) ) ) {
                    goto exit_function;
                }
                /* decrement flow on (+) charge edge of the middle =N- */
                pe = pBNS->edge + eNPlusEdge1;
                nDeltaChargeMax = 2;
            } else
            /* 4. N#N(+)-N(-)- */
            if ( pBNS->edge[eNMinusEdge].flow  == 0 &&  /* N */
                 pBNS->edge[eNMinusEdge1].flow == 0 && pBNS->edge[eNPlusEdge1].flow == 0 && /* N(+) */
                 pBNS->edge[eNMinusEdge2].flow == 1 && pBNS->edge[eNPlusEdge2].flow == 1    /* N(-) */ ) {
                /* unfix (-) edge on the 1st and 3rd N */
                if ( (ret = AddToEdgeList( &NNNChargeEdges, eNPlusEdge1, INC_EDGE_LIST )) ) {
                    goto exit_function;
                }
                if ( (ret = AddToEdgeList( &NNNChargeEdges, eNMinusEdge1, INC_EDGE_LIST )) ) {
                    goto exit_function;
                }
                if ( (ret = AddToEdgeList( &NNNChargeEdges, eNPlusEdge2, INC_EDGE_LIST )) ) {
                    goto exit_function;
                }
                /* decrement triple bond on terminal N# */
                pe = pBNS->edge + pBNS->vert[i].iedge[0];
                nDeltaChargeMax = 0;
            } else
            /* 5. N#N(+)-N(+)# */
            if ( pBNS->edge[eNMinusEdge].flow  == 0 &&  /* N */
                 pBNS->edge[eNMinusEdge1].flow == 0 && pBNS->edge[eNPlusEdge1].flow == 0 && /* N(+) */
                 pBNS->edge[eNMinusEdge2].flow == 0 && pBNS->edge[eNPlusEdge2].flow == 0    /* N(+) */ ) {
                /* unfix (-) edge on the 1st and (+) edge on the 3rd N */
                if ( (ret = AddToEdgeList( &NNNChargeEdges, eNPlusEdge1, INC_EDGE_LIST )) ) {
                    goto exit_function;
                }
                if ( (ret = AddToEdgeList( &NNNChargeEdges, eNMinusEdge1, INC_EDGE_LIST )) ) {
                    goto exit_function;
                }
                if ( (ret = AddToEdgeList( &NNNChargeEdges, eNMinusEdge2, INC_EDGE_LIST )) ) {
                    goto exit_function;
                }
                /* decrement triple bond on terminal N# */
                pe = pBNS->edge + pBNS->vert[i].iedge[0];
                nDeltaChargeMax = 0;
            } else {
                continue;
            }

            if ( NO_VERTEX != eNFlowerEdge1 && !pBNS->edge[eNFlowerEdge1].flow ) {
                if ( (ret = AddToEdgeList( &NNNChargeEdges, eNFlowerEdge1, INC_EDGE_LIST )) ) {
                    goto exit_function;
                }
            }
            if ( NO_VERTEX != eNFlowerEdge2 && !pBNS->edge[eNFlowerEdge2].flow ) {
                if ( (ret = AddToEdgeList( &NNNChargeEdges, eNFlowerEdge2, INC_EDGE_LIST )) ) {
                    goto exit_function;
                }
            }

            v1 = pe->neighbor1;
            v2 = pe->neighbor12 ^ v1;
            pe->flow --;
            pBNS->vert[v1].st_edge.flow --;
            pBNS->vert[v2].st_edge.flow --;
            pBNS->tot_st_flow -= 2;

            pe->forbidden                    |= forbidden_edge_mask;
            
            if ( !CarbonChargeEdges.num_edges ) {
                /* do not let carbon atoms get charged */
                AllocEdgeList( &CarbonChargeEdges, INC_EDGE_LIST ); 
                if ( 0 > (ret = ForbidCarbonChargeEdges( pBNS, pTCGroups, &CarbonChargeEdges, forbidden_edge_mask ))) {
                    goto exit_function;
                }
            } else {
                SetForbiddenEdgeMask( pBNS, &CarbonChargeEdges, forbidden_edge_mask );
            }
            SetForbiddenEdgeMask( pBNS, &NNNChargeEdges, forbidden_edge_mask );

            ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                  &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

            if ( ret == 1 && ((vPathEnd == v1 && vPathStart == v2) ||
                              (vPathEnd == v2 && vPathStart == v1)) && nDeltaCharge <= nDeltaChargeMax ) {
                ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                (*pnNumRunBNS) ++;
                *pnTotalDelta += ret;
                num_success ++;
                /* fix charges on N(-)=N(+)=N- */
                if ( (ret = AddToEdgeList( &CarbonChargeEdges, eNMinusEdge, INC_EDGE_LIST )) ) {
                    goto exit_function;
                }
                if ( (ret = AddToEdgeList( &CarbonChargeEdges, eNPlusEdge1, INC_EDGE_LIST )) ) {
                    goto exit_function;
                }
                if ( (ret = AddToEdgeList( &CarbonChargeEdges, eNMinusEdge1, INC_EDGE_LIST )) ) {
                    goto exit_function;
                }
                if ( (ret = AddToEdgeList( &CarbonChargeEdges, eNPlusEdge2, INC_EDGE_LIST )) ) {
                    goto exit_function;
                }
                if ( (ret = AddToEdgeList( &CarbonChargeEdges, eNMinusEdge2, INC_EDGE_LIST )) ) {
                    goto exit_function;
                }
            } else {
                pe->flow ++;
                pBNS->vert[v1].st_edge.flow ++;
                pBNS->vert[v2].st_edge.flow ++;
                pBNS->tot_st_flow += 2;
                num_failed ++;
            }
            RemoveForbiddenEdgeMask( pBNS, &CarbonChargeEdges, forbidden_edge_mask );
            RemoveForbiddenEdgeMask( pBNS, &NNNChargeEdges, forbidden_edge_mask );
            pe->forbidden &= inv_forbidden_edge_mask;
            /*pBNS->edge[eNPlusEdge].forbidden &= inv_forbidden_edge_mask;*/ /* BC: array index out of range */
        }
    }
    RemoveForbiddenEdgeMask( pBNS, &CarbonChargeEdges, forbidden_edge_mask );


exit_function:

    AllocEdgeList( &CarbonChargeEdges, EDGE_LIST_FREE );
    AllocEdgeList( &AllChargeEdges, EDGE_LIST_FREE );
    AllocEdgeList( &NNNChargeEdges, EDGE_LIST_FREE );
    AllocEdgeList( &CurNNNChargeEdges, EDGE_LIST_FREE );
    AllocEdgeList( &AllNNNTermAtoms, EDGE_LIST_FREE );
    AllocEdgeList( &AllNIIIChargeEdges, EDGE_LIST_FREE );
    return ret;
#undef INC_EDGE_LIST
}
/******************************************************************************************************/
int EliminateNitrogen5Val3Bonds(BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct,
                     inp_ATOM *at, inp_ATOM *at2, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups,
                     int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask)
{
    int i, j, k, bForbiddenCarbonCharges, ret2, ret;
    int num_at = pStruct->num_atoms;
    int num_deleted_H = pStruct->num_deleted_H;
    int len_at = num_at + num_deleted_H;
    int inv_forbidden_edge_mask = ~forbidden_edge_mask;
    EDGE_LIST CarbonChargeEdges;

    ret = 0;
    bForbiddenCarbonCharges = 0;
    AllocEdgeList( &CarbonChargeEdges, EDGE_LIST_CLEAR );
        
    memcpy( at2, at, len_at*sizeof(at2[0]));
    pStruct->at = at2;
    ret2 = CopyBnsToAtom( pStruct, pBNS, pVA, pTCGroups, 1 );
    if ( ret2 < 0 ) {
        ret = ret2;
        goto exit_function;
    }

    /* forbid creation of other N(V) atoms */
    /* fix single bonds to metals */
    for ( i = 0; i < num_at; i ++ ) {
        if ( pVA[i].cNumValenceElectrons == 5 && 
             0 <= (k = GetChargeFlowerUpperEdge( pBNS, pVA, pVA[i].nCPlusGroupEdge-1 )) &&
             1 == pBNS->edge[k].flow) {
            pBNS->edge[k].forbidden |= forbidden_edge_mask;
        } else
        if ( pVA[i].cMetal ) {
            for ( j = 0; j < at2[i].valence; j ++ ) {
                if ( BOND_TYPE_SINGLE == (at2[i].bond_type[j] & BOND_TYPE_MASK) ) {
                    pBNS->edge[pBNS->vert[i].iedge[j]].forbidden |= forbidden_edge_mask;
                }
            }
        }
    }

    /*------------------------------------------------------------------------------
                 (+)  single line => flow = 0                            (+)-(Y)=(+)super
            01  //    double line => flow = 1                fix-> 01   // <-- fix
         1 --- 0                                                 1 === 0    
          \\ //         edge eij connects vertices i < j:         \   /  02 
        12  2   02 <--- edge number: e02 connects vertices v0   12  2(..) <- double 'radical'  
            |                        v0 and v2                      |       
           =N=      vertex N has number i                          =N=      
            |                                                       |       
    --------------------------------------------------------------------------------*/
    for ( i = 0; i < num_at; i ++ ) {
        if ( pVA[i].cNumValenceElectrons == 5 && at2[i].valence == 3 &&
             at2[i].chem_bonds_valence == 5 && !at2[i].charge && !at2[i].radical &&
             !(at2[i].endpoint || (pStruct->endpoint && pStruct->endpoint[i])) && pVA[i].cnListIndex > 0 &&
             cnList[pVA[i].cnListIndex-1].bits == cn_bits_NPN &&
             pVA[i].nCPlusGroupEdge > 0 ) {
            
            Vertex v, v0 = NO_VERTEX, v1 = NO_VERTEX, v2 = NO_VERTEX;
            EdgeIndex iePlus, ie, ie12 = NO_VERTEX, ie02, ie01;
            BNS_VERTEX *pv0, *pv1, *pv2 = NULL;
            BNS_EDGE   *pePlus, *pe, *pe12 = NULL, *pe02 = NULL, *pe01 = NULL;
            Vertex     vPathStart, vPathEnd;
            int        nPathLen;
            int        nDeltaH, nDeltaCharge, nNumVisitedAtoms;
            
            iePlus = pVA[i].nCPlusGroupEdge - 1;
            pePlus = pBNS->edge + iePlus;

            v0 = IS_BNS_VT_C_GR( pBNS->vert[pePlus->neighbor1].type )?
                 (pePlus->neighbor1 ^ pePlus->neighbor12) : pePlus->neighbor1;
            pv0 = pBNS->vert + v0;
            for ( j = 0; j < pv0->num_adj_edges; j ++ ) {
                ie = pv0->iedge[j];
                if ( ie == iePlus ) {
                    continue;
                }
                pe = pBNS->edge + ie;
                if ( pe->flow == 1 && v2 == NO_VERTEX ) {
                    /* 0 - 2, edge 02 */
                    v2  = pe->neighbor12 ^ v0;
                    pv2 = pBNS->vert + v2;
                    ie02 = ie;
                    pe02 = pe;
                } else
                if ( pe->flow == 0 && v1 == NO_VERTEX ) {
                    /* 0 - 1, edge 01 */
                    v1  = pe->neighbor12 ^ v0;
                    pv1 = pBNS->vert + v2;
                    ie01 = ie;
                    pe01 = pe;
                } else {
                    ret = RI_ERR_PROGR;
                    goto exit_function;
                }
            }
            if ( v1 == NO_VERTEX || v2 == NO_VERTEX ) {
                ret = RI_ERR_PROGR;
                goto exit_function;
            }
            for ( j = 0; j < pv2->num_adj_edges; j ++ ) {
                ie = pv2->iedge[j];
                pe = pBNS->edge + ie;
                v  = pe->neighbor12 ^ v2;
                if ( v == v0 || v == i ) {
                    continue;
                } else
                if ( v == v1 && pe->flow == 1 ) {
                    /* 1 - 2, edge 12 */
                    ie12 = ie;
                    pe12 = pe;
                } else {
                    ret = RI_ERR_PROGR;
                    goto exit_function;
                }
            }
            if ( ie12 == NO_VERTEX ) {
                ret = RI_ERR_PROGR;
                goto exit_function;
            }
            /* rearrange cap and flow, forbid 2 edges */
            pe01->flow = 1;
            pe12->flow = 0;
            pe02->flow = 0;
            pv2->st_edge.flow -= 2;
            pBNS->tot_st_flow -= 2;
            pePlus->forbidden |= forbidden_edge_mask;
            pe01->forbidden   |= forbidden_edge_mask;
            
            if ( !bForbiddenCarbonCharges ) {
                if ( 0 > (ret = ForbidCarbonChargeEdges( pBNS, pTCGroups, &CarbonChargeEdges, forbidden_edge_mask ))) {
                    goto exit_function;
                }
                bForbiddenCarbonCharges = 1;
            }


            ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                  &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );
            if ( ret == 1 && vPathEnd == v2 && vPathStart == v2 && nDeltaCharge <= (pVA[i].cNumBondsToMetal? 2:0) ) {
                ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
            } else {
                pe01->flow = 0;
                pe12->flow = 1;
                pe02->flow = 1;
                pv2->st_edge.flow += 2;
                pBNS->tot_st_flow += 2;
            }
            pePlus->forbidden &= inv_forbidden_edge_mask;
            pe01->forbidden   &= inv_forbidden_edge_mask;

            if ( ret < 0 ) {
                goto exit_function;
            } else
            if ( ret ) {
                memcpy( at2, at, len_at*sizeof(at2[0]));
                pStruct->at = at2;
                ret2 = CopyBnsToAtom( pStruct, pBNS, pVA, pTCGroups, 1 );
                if ( ret2 < 0 ) {
                    ret = ret2;
                    goto exit_function;
                }
            }
        }
    }
exit_function:
    /* allow creation of other N(V) atoms */
    for ( i = 0; i < num_at; i ++ ) {
        if ( pVA[i].cNumValenceElectrons == 5 && 
             0 <= (k = GetChargeFlowerUpperEdge( pBNS, pVA, pVA[i].nCPlusGroupEdge-1 )) &&
             1 == pBNS->edge[k].flow && (pBNS->edge[k].forbidden & forbidden_edge_mask) ) {
            pBNS->edge[k].forbidden &= inv_forbidden_edge_mask;
        } else
        if ( pVA[i].cMetal ) {
            for ( j = 0; j < at2[i].valence; j ++ ) {
                if ( BOND_TYPE_SINGLE == (at2[i].bond_type[j] & BOND_TYPE_MASK) ) {
                    pBNS->edge[pBNS->vert[i].iedge[j]].forbidden &= inv_forbidden_edge_mask;
                }
            }
        }
    }
    RemoveForbiddenEdgeMask( pBNS, &CarbonChargeEdges, forbidden_edge_mask );
    AllocEdgeList( &CarbonChargeEdges, EDGE_LIST_FREE );
    return ret;
}
/******************************************************************************************************/
int Convert_SIV_to_SVI(BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct,
                     inp_ATOM *at, inp_ATOM *at2, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups,
                     int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask)
{
    int i, j, k, neigh, bForbiddenCarbonCharges, nFlowerEdge, delta, ret2, ret;
    int num_at = pStruct->num_atoms;
    int num_deleted_H = pStruct->num_deleted_H;
    int len_at = num_at + num_deleted_H;
    int inv_forbidden_edge_mask = ~forbidden_edge_mask;
    EDGE_LIST CarbonChargeEdges, FlowerEdgesList;

    ret = 0;
    bForbiddenCarbonCharges = 0;
    AllocEdgeList( &CarbonChargeEdges, EDGE_LIST_CLEAR );
    AllocEdgeList( &FlowerEdgesList, EDGE_LIST_CLEAR );
        
    memcpy( at2, at, len_at*sizeof(at2[0]));
    pStruct->at = at2;
    ret2 = CopyBnsToAtom( pStruct, pBNS, pVA, pTCGroups, 1 );
    if ( ret2 < 0 ) {
        ret = ret2;
        goto exit_function;
    }

    /* forbid creation of other S(IV) atoms */
    /* fix single bonds to metals and (N(IV), flow=1), (S(IV), flow=0) */
    for ( i = 0; i < num_at; i ++ ) {
        if ( (pVA[i].cNumValenceElectrons == 5 /* N(IV)*/ || pVA[i].cNumValenceElectrons == 6 /* S(VI)*/) && 
             0 <= (k = GetChargeFlowerUpperEdge( pBNS, pVA, pVA[i].nCPlusGroupEdge-1 )) &&
             !pBNS->edge[k].forbidden &&
             6 == pVA[i].cNumValenceElectrons + pBNS->edge[k].flow ) {

            pBNS->edge[k].forbidden |= forbidden_edge_mask;
            if ( (ret = AddToEdgeList( &FlowerEdgesList, k, 64 )) ) {
                goto exit_function;
            }
        } else
        if ( pVA[i].cMetal ) {
            for ( j = 0; j < at2[i].valence; j ++ ) {
                if ( BOND_TYPE_SINGLE == (at2[i].bond_type[j] & BOND_TYPE_MASK) ) {

                    pBNS->edge[k=pBNS->vert[i].iedge[j]].forbidden |= forbidden_edge_mask;
                    if ( (ret = AddToEdgeList( &FlowerEdgesList, k, 64 )) ) {
                        goto exit_function;
                    }
                }
            }
        } else
        /* fix bonds to neighbors of S(IV) if they are not O,S,Se,Te with 2 or more bonds */
        /* exactly same if(..) as below */
        if ( pVA[i].cNumValenceElectrons == 6 && at2[i].valence == 4 &&
             at2[i].chem_bonds_valence == 4 && !at2[i].charge && !at2[i].radical &&
             !at2[i].endpoint && pVA[i].cnListIndex > 0 &&
             cnList[pVA[i].cnListIndex-1].bits == cn_bits_NPN &&
             0 <= (nFlowerEdge = GetChargeFlowerUpperEdge( pBNS, pVA, pVA[i].nCPlusGroupEdge-1 ) ) &&
             pBNS->edge[nFlowerEdge].flow > 0 ) {

            for ( j = 0; j < at2[i].valence; j ++ ) {
                neigh = at2[i].neighbor[j];
                if ( pVA[neigh].cNumValenceElectrons != 6 && at2[neigh].valence > 1 ) {
                    k = pBNS->vert[i].iedge[j];
                    if ( !pBNS->edge[k].forbidden ) {
                        if ( (ret = AddToEdgeList( &FlowerEdgesList, k, 64 )) ) {
                            goto exit_function;
                        }
                        pBNS->edge[k].forbidden |= forbidden_edge_mask;
                    }
                }
            }
        }
    }
    /*------------------------------------------------------------------------------
                        example: struct #301,
      |          |               disconnected porphyrin with four -SO3(-)
     -S-   =>   =S=
      |          |
    
    -------------------------------------------------------------------------------*/

    /*-------------------------------------------------------------------------------
     found:                                       super(+)=(Y)        super(+)=(Y)   
                                                              \                   \           
             (+)  single line => flow = 0                     (+)                 (+)
        01  //    double line => flow = 1         fix-> 01   //             01   //             
     1 === 0      triple line => flow = 2          (.)1 --- 0(.)  --->    1 --- 0               
      \   /       edge eij connects vertices i<j:      \   /  02  run      \\ //  02            
    12  2   02 <--- edge number: e02 connects        12  2        BFS    12  2                        
       |||              vertices v0 and v2              |||                  |                  
       -S-      vertex S has number i                   -S-                 =S=                 
       / \                                              / \                 / \                 
    --------------------------------------------------------------------------------*/
    for ( i = 0; i < num_at; i ++ ) {
        if ( pVA[i].cNumValenceElectrons == 6 && at2[i].valence == 4 &&
             at2[i].chem_bonds_valence == 4 && !at2[i].charge && !at2[i].radical &&
             !at2[i].endpoint && pVA[i].cnListIndex > 0 &&
             cnList[pVA[i].cnListIndex-1].bits == cn_bits_NPN &&
             /* 01 is nFlowerEdge */
             0 <= (nFlowerEdge = GetChargeFlowerUpperEdge( pBNS, pVA, pVA[i].nCPlusGroupEdge-1 ) ) &&
             pBNS->edge[nFlowerEdge].flow > 0 ) {
            
            Vertex     v1 = NO_VERTEX, v2 = NO_VERTEX;
            BNS_VERTEX *pv1, *pv2;
            BNS_EDGE   *pe;
            Vertex     vPathStart, vPathEnd;
            int        nPathLen;
            int        nDeltaH, nDeltaCharge, nNumVisitedAtoms;

            if ( !bForbiddenCarbonCharges ) {
                if ( 0 > (ret = ForbidCarbonChargeEdges( pBNS, pTCGroups, &CarbonChargeEdges, forbidden_edge_mask ))) {
                    goto exit_function;
                }
                bForbiddenCarbonCharges = 1;
            }

            delta = 1;
            pe = pBNS->edge + nFlowerEdge;                 /* edge  01 */
            pv1 = pBNS->vert + (v1 = pe->neighbor1);       /* vertex 0 */
            pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1); /* vertex 1 */

            pe->forbidden     |= forbidden_edge_mask;
            pe->flow          -= delta;
            pv1->st_edge.flow -= delta;
            pv2->st_edge.flow -= delta;
            pBNS->tot_st_flow -= 2*delta;

            ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                  &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );
            if ( ret == 1 && 
                 ((vPathEnd == v1 && vPathStart == v2) || (vPathEnd == v2 && vPathStart == v1)) &&
                 nDeltaCharge <= (pVA[i].cNumBondsToMetal? 2:0) ) {
                ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
            } else {
                pe->forbidden     &= inv_forbidden_edge_mask;
                pe->flow          += delta;
                pv1->st_edge.flow += delta;
                pv2->st_edge.flow += delta;
                pBNS->tot_st_flow += 2*delta;
            }
            if ( ret < 0 ) {
                goto exit_function;
            } else
            if ( ret ) {
                memcpy( at2, at, len_at*sizeof(at2[0]));
                pStruct->at = at2;
                ret2 = CopyBnsToAtom( pStruct, pBNS, pVA, pTCGroups, 1 );
                if ( ret2 < 0 ) {
                    ret = ret2;
                    goto exit_function;
                }
                /* store the fixed edge to unfix it upon exit */
                if ( (ret = AddToEdgeList( &FlowerEdgesList, nFlowerEdge, 64 )) ) {
                    goto exit_function;
                }
            }
        }
    }
exit_function:
    RemoveForbiddenEdgeMask( pBNS, &CarbonChargeEdges, forbidden_edge_mask );
    AllocEdgeList( &CarbonChargeEdges, EDGE_LIST_FREE );
    RemoveForbiddenEdgeMask( pBNS, &FlowerEdgesList, forbidden_edge_mask );
    AllocEdgeList( &FlowerEdgesList, EDGE_LIST_FREE );
    return ret;
}
/******************************************************************************************************
  

  =N(+)=O       =N-O(-)
            =>
   M(q)         M(q+2)

*******************************************************************************************************/
int PlusFromDB_N_DB_O_to_Metal(BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct,
                     inp_ATOM *at, inp_ATOM *at2, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups,
                     int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask)
{
    int i, j, k, n, bForbiddenCarbonCharges, delta, ret2, ret, num_NO, num_M;
    int num_at = pStruct->num_atoms;
    int num_deleted_H = pStruct->num_deleted_H;
    int len_at = num_at + num_deleted_H;
    int inv_forbidden_edge_mask = ~forbidden_edge_mask;
    EDGE_LIST CarbonChargeEdges, NO_ChargeEdgeList, NO_EdgeList;

    Vertex     v1, v2;
    BNS_VERTEX *pv1, *pv2;
    BNS_EDGE   *pe;
    Vertex     vPathStart, vPathEnd;
    int        nPathLen;
    int        nDeltaH, nDeltaCharge, nNumVisitedAtoms;


    if ( !pTCGroups->num_metal_atoms )
        return 0;

    ret = 0;
    bForbiddenCarbonCharges = 0;
    AllocEdgeList( &CarbonChargeEdges, EDGE_LIST_CLEAR ); /* all charges */
    AllocEdgeList( &NO_ChargeEdgeList, EDGE_LIST_CLEAR ); /* charges to be changed */
    AllocEdgeList( &NO_EdgeList, EDGE_LIST_CLEAR );       /* N(+)=O edges */
        
    memcpy( at2, at, len_at*sizeof(at2[0]));
    pStruct->at = at2;
    ret2 = CopyBnsToAtom( pStruct, pBNS, pVA, pTCGroups, 1 );
    if ( ret2 < 0 ) {
        ret = ret2;
        goto exit_function;
    }
    num_NO = num_M = 0;
    /* forbid creation of other S(IV) atoms */
    /* fix single bonds to metals and (N(IV), flow=1), (S(IV), flow=0) */
    for ( i = 0; i < num_at; i ++ ) {
        if ( !pVA[i].cMetal ) {
            if ( (k = pVA[i].nCPlusGroupEdge-1) >= 0 && !pBNS->edge[k].forbidden ) {
                if ( (ret = AddToEdgeList( &CarbonChargeEdges, k, 64 )) ) {
                    goto exit_function;
                }
            }
            if ( (k = pVA[i].nCMinusGroupEdge-1) >= 0 && !pBNS->edge[k].forbidden ) {
                if ( (ret = AddToEdgeList( &CarbonChargeEdges, k, 64 )) ) {
                    goto exit_function;
                }
            }
        } else {
            num_M ++;
        }
        /*
        if ( pVA[i].cMetal ) {
            if ( (k = pVA[i].nCPlusGroupEdge-1) >= 0 && !pBNS->edge[k].forbidden ) {
                if ( ret = AddToEdgeList( &NO_ChargeEdgeList, k, 64 ) ) {
                    goto exit_function;
                }
            }
            if ( (k = pVA[i].nCMinusGroupEdge-1) >= 0 && !pBNS->edge[k].forbidden ) {
                if ( ret = AddToEdgeList( &NO_ChargeEdgeList, k, 64 ) ) {
                    goto exit_function;
                }
            }
        } else
        */
        if ( !pVA[i].cMetal &&
             pVA[i].cNumValenceElectrons == 6 &&
             at2[i].charge == 0 && !at2[i].num_H &&
             1 == at2[i].valence && 2 == at2[i].chem_bonds_valence &&
             pVA[j=at2[i].neighbor[0]].cNumValenceElectrons == 5 &&
             at2[j].charge == 1 && !at2[j].num_H &&
             2 == at2[j].valence && 4 == at2[j].chem_bonds_valence ) {
            /* found =N(+)=O */
            if ( (k = pVA[i].nCMinusGroupEdge-1) >= 0 && !pBNS->edge[k].forbidden /* O */ &&
                 (n = pVA[j].nCPlusGroupEdge -1) >= 0 && !pBNS->edge[j].forbidden /* N */ ) {
                if ( (ret = AddToEdgeList( &NO_ChargeEdgeList, k, 64 ) ) ||
                     (ret = AddToEdgeList( &NO_ChargeEdgeList, n, 64 ) ) ) {
                    goto exit_function;
                }
                k = pBNS->vert[i].iedge[0];  /* N(+)=O bond */
                if ( !pBNS->edge[k].forbidden ) { 
                    if ( (ret = AddToEdgeList( &NO_EdgeList, k, 64 )) ) {
                        goto exit_function;
                    }
                    num_NO ++;
                }
            }
        }
    }
    if ( num_M && num_NO ) {
        SetForbiddenEdgeMask( pBNS, &CarbonChargeEdges, forbidden_edge_mask );
        SetForbiddenEdgeMask( pBNS, &NO_EdgeList, forbidden_edge_mask );
        RemoveForbiddenEdgeMask( pBNS, &NO_ChargeEdgeList, forbidden_edge_mask  );
        /* now only N(+), O(-) and metal charges are allowed to change */
        for ( i = 0; i < NO_EdgeList.num_edges; i ++ ) {
            k = NO_EdgeList.pnEdges[i];
            delta = 1;
            pe = pBNS->edge + k;                 /* edge  N(+)=O */
            pv1 = pBNS->vert + (v1 = pe->neighbor1);       /* vertex 0 */
            pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1); /* vertex 1 */

            pe->flow          -= delta;
            pv1->st_edge.flow -= delta;
            pv2->st_edge.flow -= delta;
            pBNS->tot_st_flow -= 2*delta;

            ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                  &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );
            if ( ret == 1 && 
                 ((vPathEnd == v1 && vPathStart == v2) || (vPathEnd == v2 && vPathStart == v1)) &&
                 nDeltaCharge == 0 ) {
                ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
            } else {
                pe->forbidden     &= inv_forbidden_edge_mask;
                pe->flow          += delta;
                pv1->st_edge.flow += delta;
                pv2->st_edge.flow += delta;
                pBNS->tot_st_flow += 2*delta;
            }
            if ( ret < 0 ) {
                goto exit_function;
            }
        }
    }
exit_function:
    RemoveForbiddenEdgeMask( pBNS, &CarbonChargeEdges, forbidden_edge_mask );
    RemoveForbiddenEdgeMask( pBNS, &NO_EdgeList, forbidden_edge_mask );
    AllocEdgeList( &CarbonChargeEdges, EDGE_LIST_FREE );
    AllocEdgeList( &NO_EdgeList, EDGE_LIST_FREE );
    AllocEdgeList( &NO_ChargeEdgeList, EDGE_LIST_FREE );
    return ret;
}
/******************************************************************************************************/
int MoveMobileHToAvoidFixedBonds(BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct,
                              inp_ATOM *at, inp_ATOM *at2, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups,
                              int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask)
{
    int ret2, ret;
    int num_at = pStruct->num_atoms;
    int num_deleted_H = pStruct->num_deleted_H;
    int len_at = num_at + num_deleted_H;
    int nNumFixedEdges, nNumAdjEdges;
    
    ret = 0;

    if ( pTCGroups->num_tgroups ) {
        memcpy( at2, at, len_at*sizeof(at2[0]));
        pStruct->at = at2;
        ret2 = CopyBnsToAtom( pStruct, pBNS, pVA, pTCGroups, 1 );
        pStruct->at = at;
        if ( ret2 < 0 ) {
            ret = ret2;
            goto exit_function;
        }
#if ( FIND_RING_SYSTEMS == 1 )
        ret2 = MarkRingSystemsInp( at2, num_at, 0 );
        if ( ret2 < 0 ) {
            ret = ret2;
            goto exit_function;
        }
#endif
        /* --- forbidden edges --- */
        ret2 = SetForbiddenEdges( pBNS, at2, num_at, forbidden_edge_mask );
        if ( ret2 < 0 ) {
            ret2 = -(ret + 1);
        }
        nNumFixedEdges = ret2;
        ret = AdjustTgroupsToForbiddenEdges2( pBNS, at2, pVA, num_at, forbidden_edge_mask );
        nNumAdjEdges = ret;
        if ( ret ) {
            pBNS->edge_forbidden_mask |= forbidden_edge_mask;
            ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
            (*pnNumRunBNS) ++;
            if ( ret < 0 ) {
                goto exit_function;
            } else {
                *pnTotalDelta += ret;
            }
        }
        if ( nNumFixedEdges || nNumAdjEdges ) {
            /* removes this edge mask from ALL edges */
            RemoveForbiddenBondFlowBits( pBNS, forbidden_edge_mask );
        }
    }

exit_function:

    return ret;
}
/******************************************************************************************************/
/* Find and eliminate cases when Mobile H endpoint has radical on it (typical for wrong P(VI)(=O)3OH  */
int RemoveRadFromMobileHEndpoint(BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct,
                              inp_ATOM *at, inp_ATOM *at2, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups,
                              int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask)
{
    int i, num_fixes, tot_num_fixes = 0;
    
    int ret2, ret;
    int num_at = pStruct->num_atoms;
    int num_deleted_H = pStruct->num_deleted_H;
    int len_at = num_at + num_deleted_H;

    int         itg, j, k, n, m;
    Vertex      vtg1, endpoint0=NO_VERTEX, endpoint1, endpoint2, centerpoint;
    Vertex      centerpoint_found=NO_VERTEX;
    BNS_VERTEX *ptg1, *pEndp0=NULL, *pEndp1, *pEndp2, *pCentp, *pCentp_found, *pEndp2_found=NULL;
    BNS_EDGE   *etg0=NULL, *etg1, *etg2, *ecp0, *ecp1, *ecp2;
    BNS_EDGE   *etg1_found=NULL, *ecp0_found=NULL, *ecp1_found=NULL, *ecp2_found=NULL;
    int         tgroup_number, num_endpoints;
    
    ret = 0;

    memcpy( at2, at, len_at*sizeof(at2[0]));
    pStruct->at = at2;
    ret2 = CopyBnsToAtom( pStruct, pBNS, pVA, pTCGroups, 1 );
    if ( ret2 < 0 ) {
        ret = ret2;
        goto exit_function;
    }
    while ( pBNS->tot_st_cap > pBNS->tot_st_flow && pTCGroups->num_tgroups ) {
        num_fixes = 0;
        for ( itg = 0; itg < pTCGroups->num_tgroups; itg ++ ) {
            pCentp_found=NULL;
            tgroup_number = pTCGroups->pTCG[itg].ord_num;
            vtg1 = pTCGroups->pTCG[itg].nVertexNumber;      /* taut group vertex index */
            ptg1 = pBNS->vert + vtg1;                       /* taut group vertex */
            num_endpoints = pTCGroups->pTCG[itg].num_edges;
            for ( i = 0; i < num_endpoints; i ++ ) {
                etg0 = pBNS->edge + ptg1->iedge[i];         /* edge from t-group to endpoint */
                endpoint0 = etg0->neighbor12 ^ vtg1;        /* taut endpoint vertex index */
                pEndp0 = pBNS->vert + endpoint0;            /* taut endpoint vertex (possible location of mobile H */
                if ( pEndp0->st_edge.cap > pEndp0->st_edge.flow ) {
                    /* radical endpoint1 has been detected */
                    /* find a 1-3 centerpoint that has two or more endpoints */
                    /* connected to the t-group vertex by edges with flow>0 and */
                    /* to the centerpoint by edges with flow = 0 */
                    /* after that: (1) increment etg1 flow to eliminate radical */
                    /* (2) increment flow on one of the two other edges to the t-group */
                    /* (3) increment st_cap on the found centerpoint */
                    /* (4) rerun the BNS and re-create the structure */
                    break;
                }
            }
            if ( i == num_endpoints ) {
                continue;
            }
            if ( i < num_endpoints ) {
                /* tautomeric endpoint found; traverse its t-group edges */
                for ( j = 0; j < num_endpoints; j ++ ) {
                    if ( i == j ) {
                        continue; /* avoid the already found radical endpoint */
                    }
                    etg1 = pBNS->edge + ptg1->iedge[j];  /* another edge from t-group to another endpoinr */
                    endpoint1 = etg1->neighbor12 ^ vtg1; /* another endpoint vertex index */
                    pEndp1 = pBNS->vert + endpoint1;     /* another endpoint vertex */
                    if ( pEndp1->st_edge.cap > pEndp1->st_edge.flow ) {
                        continue; /* one more radical-endpoint! What is going on here??? */
                    }
                    if ( !etg1->flow ) {
                        continue; /* avoid enpoints that do not have an attachment */
                    }
                    if ( !(pEndp1->type & BNS_VERT_TYPE_ENDPOINT) ) {
                        continue; /* should not happen */
                    }
                    /* traverse endpoint1 edges to find a single bond connecting it to the centerpoint */
                    for ( k = 0; k < at2[endpoint1].valence; k ++ ) {
                        ecp1 = pBNS->edge + pEndp1->iedge[k];
                        if ( ecp1->flow ) {
                            continue;
                        }
                        centerpoint = ecp1->neighbor12 ^ endpoint1;
                        pCentp = pBNS->vert + centerpoint;
                        /* traverse centerpoint edges to find a single bond to the 2nd endpoint */
                        for ( n = 0; n < at2[centerpoint].valence; n ++ ) {
                            ecp2 = pBNS->edge + pCentp->iedge[n];
                            if ( ecp2->flow ) {
                                continue;
                            }
                            endpoint2 = ecp2->neighbor12 ^ centerpoint;
                            if ( endpoint2 <= endpoint1 || !pVA[endpoint2].nTautGroupEdge ) {
                                continue; /* don't go back: neighbors are in order of ascending ord. numbers */
                            }
                            pEndp2 = pBNS->vert + endpoint2;
                            if ( !(pEndp2->type & BNS_VERT_TYPE_ENDPOINT) ) {
                                continue;
                            }
                            etg2 = pBNS->edge + pVA[endpoint2].nTautGroupEdge - 1;
                            if ( !etg2->flow || (etg2->neighbor12 ^ endpoint2) != vtg1 ) {
                                continue;
                            }
                            /* we have found the path:
                                              Endp1                                Endp1        
                                       etg1 //     \ ecp1                   etg1 /     \\ ecp1  
                                  etg0     //       \                  etg0     /       \\      
                             Endp0-----tg1           Centp  -->   Endp0=====tg1           Centp 
                              ^            \\       /                           \\       /      
                      radical |        etg2 \\     / ecp2                   etg2 \\     / ecp2  
                                              Endp2                                Endp2        
                            */                                                                  

                            /* compare centerpoints */
                            if ( !pCentp_found ||
                                 /* try to avoid carbons */
                                 ((pVA[centerpoint].cNumValenceElectrons != 4 ||
                                  pVA[centerpoint].cPeriodicRowNumber   != 1) &&
                                 pVA[centerpoint_found].cNumValenceElectrons == 4 &&
                                 pVA[centerpoint_found].cPeriodicRowNumber   == 1) ||
                                 /* try a better non-carbon */
                                 ((pVA[centerpoint].cNumValenceElectrons != 4 ||
                                  pVA[centerpoint].cPeriodicRowNumber   != 1  ) &&
                                 (at[centerpoint].valence >  at[centerpoint_found].valence ||
                                  (at[centerpoint].valence == at[centerpoint_found].valence &&
                                  at[centerpoint].el_number > at[centerpoint_found].el_number))) ) {

                                pCentp_found = pCentp;
                                etg1_found   = etg1;
                                ecp1_found   = ecp1;
                                centerpoint_found = centerpoint;
                                break;
                            }                                
                        }
                    }
                }
            }
            if ( pCentp_found ) {
                /* ---- (1) */
                etg0->flow ++;
                pEndp0->st_edge.flow ++;
                /* ---- (2) */
                etg1_found->flow --;
                /* ---- (3) */
                ecp1_found->flow ++;
                /* ---- (4) */
                pCentp_found->st_edge.flow ++;
                pCentp_found->st_edge.cap  ++;

                pBNS->tot_st_flow += 2;
                pBNS->tot_st_cap  += 1;
                pCentp_found = NULL;
                num_fixes ++;
                tot_num_fixes ++;   /* #1 Mob-H */
                continue;
            }

            /* 2nd attempt: increment flow in centerpoint---radical_endpint edge */
            if ( i < num_endpoints ) {
                /* tautomeric endpoint found; traverse its t-group edges */
                for ( j = 0; j < num_endpoints; j ++ ) {
                    if ( i == j ) {
                        continue; /* avoid the found radical endpoint */
                    }
                    etg1 = pBNS->edge + ptg1->iedge[j];
                    endpoint1 = etg1->neighbor12 ^ vtg1;
                    pEndp1 = pBNS->vert + endpoint1;     /* another endpoint */
                    if ( pEndp1->st_edge.cap > pEndp1->st_edge.flow ) {
                        continue; /* one more radical-endpoint! What is going on here??? */
                    }
                    if ( !etg1->flow ) {
                        continue; /* avoid enpoints that do not have an attachment */
                    }
                    if ( !(pEndp1->type & BNS_VERT_TYPE_ENDPOINT) ) {
                        continue; /* should not happen */
                    }
                    /* traverse endpoint1 edges to find the edge connecting it to the centerpoint */
                    for ( k = 0; k < at2[endpoint1].valence; k ++ ) {
                        ecp1 = pBNS->edge + pEndp1->iedge[k];
                        if ( ecp1->flow ) {
                            continue;
                        }
                        centerpoint = ecp1->neighbor12 ^ endpoint1;
                        pCentp = pBNS->vert + centerpoint;
                        if ( pCentp->type & BNS_VERT_TYPE_ENDPOINT ) {
                            continue; /* do not set another endpoint's valence = an unusual value */
                        }

                        /* traverse centerpoint edges to find edge connecting it to the endpoint0 */
                        ecp2 = NULL;
                        pEndp2 = NULL;
                        for ( n = 0; n < at2[centerpoint].valence; n ++ ) {
                            ecp0 = pBNS->edge + pCentp->iedge[n];
                            if ( ecp0->flow ) {
                                endpoint2 = ecp0->neighbor12 ^ centerpoint;
                                if ( (pBNS->vert[endpoint2].type & BNS_VERT_TYPE_ENDPOINT) ) {
                                    continue; /* ignore endpoint2 if it is tautomeric endpoint */
                                }
                                /* check whether ecp0 is stereogenic: if it is then we cannot decrement its flow */
                                for ( m = 0; m < MAX_NUM_STEREO_BONDS && at[centerpoint].sb_parity[m]; m ++ ) {
                                    if ( at[centerpoint].sb_ord[m] == n ) {
                                        endpoint2 = NO_VERTEX;
                                        break;
                                    }
                                }
                                if ( endpoint2 == NO_VERTEX ) {
                                    continue;
                                }
                                pEndp2 = pBNS->vert + endpoint2;  /* found */
                                ecp2   = ecp0;
                                break;
                            }
                        }
                        for ( n = 0; n < at[centerpoint].valence; n ++ ) {
                            ecp0 = pBNS->edge + pCentp->iedge[n];
                            if ( ecp0->flow ) {
                                continue;
                            }
                            if ( endpoint0 != (ecp0->neighbor12 ^ centerpoint) ) {
                                continue;
                            }                    
                            /* Found:
                                              Endp2(not endpoint)           Endp2(not radical)
                                               ||                             |
                                               ||ecp2                         |ecp2
                                          ecp0 ||  ecp1                 ecp0  |  ecp1     
                                     Endp0----Centp----Endp1       Endp0====Centp----Endp1
                                      ^  \             /      -->      \             /    
                              radical |   \           /                 \           /     
                                           \         /                   \         /      
                                       etg0 \       / etg1           etg0 \       / etg1  
                                             \     /                       \     /        
                                               tg1                           tg1          

                             */

                            /* compare centerpoints */
                            if ( !pCentp_found ||
                                 /* try to avoid carbons */
                                 ((pVA[centerpoint].cNumValenceElectrons != 4 ||
                                  pVA[centerpoint].cPeriodicRowNumber   != 1) &&
                                 pVA[centerpoint_found].cNumValenceElectrons == 4 &&
                                 pVA[centerpoint_found].cPeriodicRowNumber   == 1) ||
                                 /* try a better non-carbon */
                                 ((pVA[centerpoint].cNumValenceElectrons != 4 ||
                                  pVA[centerpoint].cPeriodicRowNumber   != 1  ) &&
                                 (at[centerpoint].valence >  at[centerpoint_found].valence ||
                                  (at[centerpoint].valence == at[centerpoint_found].valence &&
                                  at[centerpoint].el_number > at[centerpoint_found].el_number))) ) {

                                pCentp_found = pCentp;
                                etg1_found   = etg1;
                                ecp0_found   = ecp0;
                                centerpoint_found = centerpoint;
                                pEndp2_found = pEndp2;
                                ecp2_found   = ecp2;
                                break;
                            }                                
                        }
                    }
                }
            }
            if ( pCentp_found ) {
                ecp0_found->flow ++;
                if ( ecp0_found->cap < ecp0_found->flow ) {
                    ecp0_found->cap = ecp0_found->flow;
                }
                pEndp0->st_edge.flow ++;
                if ( pEndp2_found && ecp2_found ) {
                    ecp2_found->flow --;
                    pEndp2_found->st_edge.flow --;
                } else {
                    /* Endp2 not found */
                    pCentp_found->st_edge.flow ++;
                    pCentp_found->st_edge.cap  ++;
                    pBNS->tot_st_flow += 2; /* radical elimination */
                    pBNS->tot_st_cap  += 1;
                }

                pCentp_found = NULL;
                num_fixes ++;
                tot_num_fixes ++;      /* #2 Mob-H */
                continue;
            }
            /* 3rd attempt: find =C= and move radical to it */
            if ( i < num_endpoints ) {
                int jj, delta, bNotFixed = 1;
                Vertex     vPathStart, vPathEnd, v1, v2;
                int        nPathLen, nDeltaH, nDeltaCharge, nNumVisitedAtoms;
                for ( jj = 0; jj < num_at && bNotFixed; jj ++ ) {
                    if ( at2[i].endpoint ) {
                        continue;
                    }
                    if ( 2 == at2[jj].valence && pBNS->vert[jj].st_edge.cap == pBNS->vert[jj].st_edge.flow &&
                         4 == pVA[jj].cNumValenceElectrons &&
                         !(ecp0 = pBNS->edge + pBNS->vert[jj].iedge[0])->forbidden &&
                         !(ecp1 = pBNS->edge + pBNS->vert[jj].iedge[1])->forbidden &&
                         1 == ecp0->flow && 1 == ecp1->flow &&
                         !at2[(int)at2[i].neighbor[0]].sb_parity[0] &&
                         !at2[(int)at2[i].neighbor[1]].sb_parity[0]  ) {
                        /* found =C=; make a radical and try to cancel the two radicals */
                        k = ecp0->neighbor12 ^ jj;
                        if ( at2[k].endpoint ) {
                            ecp0 = ecp1;
                            k = ecp0->neighbor12 ^ jj;
                            if ( at2[k].endpoint ) {
                                continue;
                            }
                        }
                        delta = 1;
                        /* decrement C valence */
                        pBNS->vert[jj].st_edge.flow -= delta;
                        pBNS->vert[jj].st_edge.cap  -= delta;
                        /* decrement bond order */
                        ecp0->flow                 -= delta;
                        /* reflect the changes in at2[k] to make it a radical */
                        pBNS->vert[k].st_edge.flow -= delta;
                        pBNS->tot_st_cap           -= delta;
                        pBNS->tot_st_flow          -= 2*delta;

                        v1 = endpoint0;
                        v2 = k;

                        ret2 = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                              &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

                        if ( ret2 == 1 && ((vPathEnd == v1 && vPathStart == v2) ||
                                           (vPathEnd == v2 && vPathStart == v1)) && nDeltaCharge == 0 ) {
                            ret2 = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                            if ( ret2 > 0 ) {
                                num_fixes ++;
                                tot_num_fixes ++;   /* #3 Mob-H */
                                pBNS->vert[jj].st_edge.cap += delta; /* create radical on =C- */
                                pBNS->tot_st_cap          += delta;
                                pCentp_found  = NULL;
                                bNotFixed     = 0; /* exit from the cycle */
                                break;
                            }
                        } else {
                            /* failed */
                            pBNS->vert[jj].st_edge.flow += delta;
                            pBNS->vert[jj].st_edge.cap  += delta;
                            /* decrement bond order */
                            ecp0->flow                 += delta;
                            /* reflect the changes in at2[k] to make it a radical */
                            pBNS->vert[k].st_edge.flow += delta;
                            pBNS->tot_st_cap           += delta;
                            pBNS->tot_st_flow          += 2*delta;
                        }
                        if ( ret2 < 0 ) {
                            ret = ret2;
                            goto exit_function;
                        }
                    }
                }
            }
        }
        if ( !num_fixes ) {
            break;
        }
    }
    ret = tot_num_fixes;
exit_function:
    pStruct->at = at;
    memcpy( at2, at, len_at*sizeof(at2[0]));
    return ret;
}
/******************************************************************************************************/
/* Find and eliminate cases when Mobile H endpoint has radical on it (typical for wrong P(VI)(=O)3OH  */
int RemoveRadFromMobileHEndpointFixH(BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct,
                              inp_ATOM *at, inp_ATOM *at2, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups,
                              int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask)
{
#define IS_C(x) (NO_VERTEX != x && pVA[x].cNumValenceElectrons == 4 && pVA[x].cPeriodicRowNumber == 1)
    int i, num_fixes, tot_num_fixes = 0;
    
    int ret2, ret;
    int num_at = pStruct->num_atoms;
    int num_deleted_H = pStruct->num_deleted_H;
    int len_at = num_at + num_deleted_H;
    int inv_forbidden_edge_mask = ~forbidden_edge_mask;
    EDGE_LIST ChargeEdgeList, BondEdgeList;
    int         itg, j, k, n, m, num_endp;
    Vertex      endpoint0=NO_VERTEX, endpoint1, endpoint2=NO_VERTEX, centerpoint;
    Vertex      centerpoint_found=NO_VERTEX, endpoint2_found=NO_VERTEX;
    BNS_VERTEX *pEndp0=NULL, *pEndp1, *pEndp2, *pCentp, *pCentp_found, *pEndp2_found=NULL;
    BNS_EDGE   *ecp0, *ecp1, *ecp2, *ecp0_found=NULL, *ecp1_found=NULL, *ecp2_found=NULL;
    int          tgroup_number, num_endpoints;
    
    ret = 0;

    if ( pStruct->iMobileH != TAUT_NON )
        return ret;

    AllocEdgeList( &ChargeEdgeList, EDGE_LIST_CLEAR);
    AllocEdgeList( &BondEdgeList, EDGE_LIST_CLEAR);

    memcpy( at2, at, len_at*sizeof(at2[0]));
    pStruct->at = at2;
    ret2 = CopyBnsToAtom( pStruct, pBNS, pVA, pTCGroups, 1 );
    if ( ret2 < 0 ) {
        ret = ret2;
        goto exit_function;
    }
    while ( pBNS->tot_st_cap > pBNS->tot_st_flow && pStruct->ti.num_t_groups ) {
        int iEndpoint = 0;
        num_fixes = 0;
        for ( itg = 0; itg < pStruct->ti.num_t_groups; iEndpoint += num_endpoints, itg ++ ) {
            pCentp_found=NULL;
            tgroup_number = pStruct->ti.t_group[itg].nGroupNumber;
            num_endpoints = pStruct->ti.t_group[itg].nNumEndpoints;
            for ( i = 0; i < num_endpoints; i ++ ) {
                endpoint0 = pStruct->ti.nEndpointAtomNumber[iEndpoint+i];
                pEndp0 = pBNS->vert + endpoint0;  /* taut endpoint vertex (possible location of mobile H */
                if ( pEndp0->st_edge.cap > pEndp0->st_edge.flow ) {
                    /* radical endpoint1 has been detected */
                    /* find a 1-3 centerpoint that has two or more endpoints */
                    /* connected to the t-group vertex by edges with flow>0 and */
                    /* to the centerpoint by edges with flow = 0 */
                    /* after that: (1) increment etg1 flow to eliminate radical */
                    /* (2) increment flow on one of the two other edges to the t-group */
                    /* (3) increment st_cap on the found centerpoint */
                    /* (4) rerun the BNS and re-create the structure */
                    break;
                }
            }

            /* 2nd attempt: increment flow in centerpoint---radical_endpoint edge */
            pCentp_found = NULL;
            if ( i < num_endpoints ) {
                /* tautomeric endpoint found; traverse its t-group edges */
                for ( j = 0; j < num_endpoints; j ++ ) {
                    if ( i == j ) {
                        continue; /* avoid the found radical endpoint */
                    }
                    endpoint1 = pStruct->ti.nEndpointAtomNumber[iEndpoint+j];
                    pEndp1 = pBNS->vert + endpoint1;     /* another endpoint */

                    if ( pEndp1->st_edge.cap > pEndp1->st_edge.flow ) {
                        continue; /* one more radical-endpoint! What is going on here??? */
                    }
                    if ( !at2[endpoint1].num_H && at2[endpoint1].charge != -1 ) {
                        continue; /* avoid enpoints that do not have an attachment */
                    }
                    if ( !pStruct->endpoint[endpoint1] ) {
                        continue; /* should not happen */
                    }
                    /* traverse endpoint1 edges to find the edge connecting it to the centerpoint */
                    for ( k = 0; k < pEndp1->num_adj_edges; k ++ ) {
                        ecp1 = pBNS->edge + pEndp1->iedge[k];
                        if ( ecp1->flow ) {
                            continue;
                        }
                        centerpoint = ecp1->neighbor12 ^ endpoint1;
                        if ( centerpoint >= pBNS->num_atoms ) {
                            break; /* no more edges to atoms */
                        }
                        pCentp = pBNS->vert + centerpoint;
                        if ( pStruct->endpoint[centerpoint] ) {
                            continue; /* do not set another endpoint's valence = an unusual value */
                        }
                        /* traverse centerpoint edges to find edge connecting it to the endpoint0 */
                        /* 1. Find a double bond to an endpoint */
                        ecp2 = NULL;
                        pEndp2 = NULL;
                        for ( n = 0, num_endp = 0; n < at2[centerpoint].valence; n ++ ) {
                            ecp0 = pBNS->edge + pCentp->iedge[n];
                            if ( ecp0->flow ) {
                                endpoint2 = ecp0->neighbor12 ^ centerpoint;
                                if ( pStruct->endpoint[endpoint2] /* ??? */ ) {
                                    continue;
                                }
                                /* check whether ecp0 is stereogenic: if it is then we cannot decrement its flow */
                                for ( m = 0; m < MAX_NUM_STEREO_BONDS && at[centerpoint].sb_parity[m]; m ++ ) {
                                    if ( at[centerpoint].sb_ord[m] == n ) {
                                        endpoint2 = NO_VERTEX;
                                        break;
                                    }
                                }
                                if ( endpoint2 == NO_VERTEX ) {
                                    continue;
                                }
                                pEndp2 = pBNS->vert + endpoint2;
                                ecp2   = ecp0;
                                break;
                            }
                        }
                        if ( !ecp2 ) {
                            continue;
                        }
                        /* 2. Find a single bond to an endpoint0 */
                        for ( n = 0, num_endp = 0; n < at2[centerpoint].valence; n ++ ) {
                            ecp0 = pBNS->edge + pCentp->iedge[n];
                            if ( ecp0->flow ) {
                                continue;
                            }
                            if ( endpoint0 != (ecp0->neighbor12 ^ centerpoint) ) {
                                continue;
                            }                    
                            /* Found:
                                              Endp2                         Endp2(not radical)
                                               ||                             |
                                               ||ecp2                         |ecp2
                                          ecp0 ||  ecp1                 ecp0  |  ecp1     
                                     Endp0----Centp----Endp1       Endp0====Centp----Endp1
                                      ^  \             /      -->      \             /    
                              radical |   \           /                 \           /     
                                           \         /                   \         /      
                                       etg0 \       / etg1           etg0 \       / etg1  
                                             \     /                       \     /        
                                               tg1                           tg1          

                             */

                            /* compare centerpoints */
                            if ( !pCentp_found ||
                                 /* try to avoid carbons */
                                 ((pVA[centerpoint].cNumValenceElectrons != 4 ||
                                  pVA[centerpoint].cPeriodicRowNumber   != 1) &&
                                 pVA[centerpoint_found].cNumValenceElectrons == 4 &&
                                 pVA[centerpoint_found].cPeriodicRowNumber   == 1) ||
                                 /* try a better non-carbon */
                                 ((pVA[centerpoint].cNumValenceElectrons != 4 ||
                                  pVA[centerpoint].cPeriodicRowNumber   != 1  ) &&
                                 (at[centerpoint].valence >  at[centerpoint_found].valence ||
                                  (at[centerpoint].valence == at[centerpoint_found].valence &&
                                  at[centerpoint].el_number > at[centerpoint_found].el_number))) ) {

                                pCentp_found = pCentp;
                                ecp0_found   = ecp0;
                                centerpoint_found = centerpoint;
                                pEndp2_found = pEndp2;
                                ecp2_found   = ecp2;
                                break;
                            }                                
                        }
                    }
                }
            }
            /* decrement st_flow, st_cap on Endp2; decrement flow on ecp2; decrement st_flow on Centp */
            /* result: radicals on Endp0 and Centp => run BNS */
            if ( pCentp_found ) {
                /* make ecp0 a double bond, make ecp2 a single bond, remove radical */
                ecp0_found->flow ++;
                if ( ecp0_found->cap < ecp0_found->flow ) {
                    ecp0_found->cap = ecp0_found->flow;
                }
                pEndp0->st_edge.flow ++;
                if ( pEndp2_found && ecp2_found ) {
                    ecp2_found->flow --;
                    pEndp2_found->st_edge.flow --;
                } else {
                    /* Endp2 not found: only make ecp0 a double bond */
                    pCentp_found->st_edge.flow ++;
                    pCentp_found->st_edge.cap  ++;
                    pBNS->tot_st_flow += 2; /* radical elimination */
                    pBNS->tot_st_cap  += 1;
                }

                pCentp_found = NULL;
                num_fixes ++;    /* #2 */
                tot_num_fixes ++;
                continue;
            }
            /* 1st attempt */
            pCentp_found = NULL;
            if ( i < num_endpoints ) {
                /* tautomeric endpoint found; traverse its t-group edges */
                for ( j = 0; j < num_endpoints; j ++ ) {
                    if ( i == j ) {
                        continue; /* avoid the found radical endpoint */
                    }
                    endpoint1 = pStruct->ti.nEndpointAtomNumber[iEndpoint+j];
                    pEndp1 = pBNS->vert + endpoint1;     /* another endpoint */
                    if ( pEndp1->st_edge.cap > pEndp1->st_edge.flow ) {
                        continue; /* one more radical-endpoint! What is going on here??? */
                    }
                    if ( !at2[endpoint1].num_H && at2[endpoint1].charge != -1 ) {
                        continue; /* avoid enpoints that do not have an attachment */
                    }
                    if ( !pStruct->endpoint[endpoint1] ) {
                        continue; /* should not happen */
                    }
                    /* traverse endpoint1 edges to find the edge connecting it to the centerpoint */
                    for ( k = 0; k < pEndp1->num_adj_edges; k ++ ) {
                        ecp1 = pBNS->edge + pEndp1->iedge[k];
                        if ( ecp1->flow ) {
                            continue;
                        }
                        centerpoint = ecp1->neighbor12 ^ endpoint1;
                        if ( centerpoint >= pBNS->num_atoms ) {
                            break;
                        }
                        pCentp = pBNS->vert + centerpoint;
                        /* traverse centerpoint edges to find the 2nd endpoint */
                        for ( n = 0, num_endp = 0; n < pCentp->num_adj_edges; n ++ ) {
                            ecp2 = pBNS->edge + pCentp->iedge[n];
                            if ( ecp2->flow ) {
                                continue;
                            }
                            endpoint2 = ecp2->neighbor12 ^ centerpoint;
                            if ( endpoint2 >= pBNS->num_atoms ) {
                                break;
                            }
                            if ( !pStruct->endpoint[endpoint2] ) {
                                continue; 
                            }
                            pEndp2 = pBNS->vert + endpoint2;
                            
                            if ( at2[endpoint2].num_H || at2[endpoint1].charge == -1 ) {
                                continue;
                            }

                            /* we have found the path:
                            
                            Endp1 has no attachments, Endp2 has.

                                              Endp1                                Endp1        
                                       etg1 //     \ ecp1                   etg1 /     \\ ecp1  
                                  etg0     //       \                  etg0     /       \\      
                             Endp0-----tg1           Centp  -->   Endp0=====tg1           Centp 
                              ^            \\       /                           \\       /      
                      radical |        etg2 \\     / ecp2                   etg2 \\     / ecp2  
                                              Endp2                                Endp2        
                            */                                                                  

                            /* compare centerpoints */
                            if ( !pCentp_found ||
                                 /* try to avoid carbons */
                                 ((pVA[centerpoint].cNumValenceElectrons != 4 ||
                                  pVA[centerpoint].cPeriodicRowNumber   != 1) &&
                                 pVA[centerpoint_found].cNumValenceElectrons == 4 &&
                                 pVA[centerpoint_found].cPeriodicRowNumber   == 1) ||
                                 /* try a better non-carbon */
                                 ((pVA[centerpoint].cNumValenceElectrons != 4 ||
                                  pVA[centerpoint].cPeriodicRowNumber   != 1  ) &&
                                 (at[centerpoint].valence >  at[centerpoint_found].valence ||
                                  (at[centerpoint].valence == at[centerpoint_found].valence &&
                                  at[centerpoint].el_number > at[centerpoint_found].el_number))) ) {

                                pCentp_found = pCentp;
                                ecp1_found   = ecp1;
                                centerpoint_found = centerpoint;
                                break;
                            }                                
                        }
                    }
                }
            }
            if ( pCentp_found ) {
                /* create a new radical at the centerpoint and try to cancel them */
                int delta = 1, ret3;
                Vertex     vPathStart, vPathEnd, v1, v2;
                int        nPathLen, nDeltaH, nDeltaCharge, nNumVisitedAtoms;
                
                pCentp_found->st_edge.cap += delta;
                pBNS->tot_st_cap          += delta;

                v1 = pCentp_found - pBNS->vert;
                v2 = pEndp0 - pBNS->vert;
                ret3 = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                      &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

                if ( ret3 == 1 && ((vPathEnd == v1 && vPathStart == v2) ||
                                   (vPathEnd == v2 && vPathStart == v1)) && nDeltaCharge % 2 == 0 ) {
                    ret3 = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                    if ( ret3 > 0 ) {
                        num_fixes ++;
                        tot_num_fixes ++;   /* #1 */
                        pCentp_found = NULL;
                        continue;
                    }
                } else {
                    pCentp_found->st_edge.cap -= delta;
                    pBNS->tot_st_cap          -= delta;
                }
                if ( ret3 < 0 ) {
                    ret = ret3;
                    goto exit_function;
                }
            }
            /*----------------------------------------------------------------------------------------
               3rd attempt:                                                         add radical to keep
              ==============           u,f=>unfixed, fixed edges                    (N electrons)%2
                        (-)                       (-)                      (-)          |     
                     e0/  \\ e1          =>     u/  \\u            =>     //  \         v     
                      /    \\ ecp1 ecp2         /    \\   u    f         //    \             
                  C--X*     Y(-)--C==Z      C--X*     Y(-)--C*--Z      --X(-)   Y===C---Z*   
                    rad.    endp     not      rad.   endp  rad not      rad.   endp    not 
                   Endp0   Endp1     endp     endp         to  endp     endp           endp
                                    Endp2                cancel                              
             
               Note: endpoints X and Y may belong to different t-groups
             ----------------------------------------------------------------------------------------*/
            pCentp_found = NULL;
            if ( i < num_endpoints ) {
                int  e0, e1;
                if ( (e0=pVA[endpoint0].nCMinusGroupEdge-1)<0 || pBNS->edge[e0].forbidden ) {
                    continue; /* no negative charge on Endp0 is possible */
                }
                /* a radical-tautomeric endpoint found; traverse all endpoints */
                for ( j = 0; j < pStruct->ti.nNumEndpoints; j ++ ) {
                    if ( iEndpoint+i == j ) {
                        continue; /* avoid the found radical endpoint */
                    }

                    endpoint1 = pStruct->ti.nEndpointAtomNumber[j];
                    pEndp1 = pBNS->vert + endpoint1;     /* another endpoint */

                    if ( pEndp1->st_edge.cap > pEndp1->st_edge.flow ) {
                        continue; /* one more radical-endpoint! What is going on here??? */
                    }
                    if ( ((e1=pVA[endpoint1].nCMinusGroupEdge-1)<0 || !pBNS->edge[e1].flow) || pBNS->edge[e1].forbidden ) {
                        continue; /* no negative charge on Endp1 */
                    }
                    if ( !pStruct->endpoint[endpoint1] ) {
                        continue; /* should not happen */
                    }
                    /* traverse endpoint1 edges to find the edge connecting it to the centerpoint */
                    for ( k = 0; k < pEndp1->num_adj_edges; k ++ ) {
                        ecp1 = pBNS->edge + pEndp1->iedge[k];   /* e1C */
                        if ( ecp1->flow || ecp1->forbidden ) {
                            continue;
                        }
                        centerpoint = ecp1->neighbor12 ^ endpoint1;
                        if ( centerpoint >= pBNS->num_atoms ) {
                            break; /* no more edges to atoms */
                        }
                        pCentp = pBNS->vert + centerpoint;
                        if ( pStruct->endpoint[centerpoint] ) {
                            continue; /* do not set another endpoint's valence = an unusual value */
                        }
                        /* traverse centerpoint edges to find edge connecting it to the endpoint0 */
                        /* 1. Find a double bond to a not endpoint */
                        ecp2 = NULL;
                        pEndp2 = NULL;
                        for ( n = 0, num_endp = 0; n < pCentp->num_adj_edges; n ++ ) {
                            ecp0 = pBNS->edge + pCentp->iedge[n];
                            if ( ecp0->flow && !ecp0->forbidden ) {
                                endpoint2 = ecp0->neighbor12 ^ centerpoint;
                                if ( endpoint2 >= pBNS->num_atoms || pStruct->endpoint[endpoint2] ) {
                                    continue;
                                }
                                /* check whether ecp0 is stereogenic: if it is then we cannot decrement its flow */
                                for ( m = 0; m < MAX_NUM_STEREO_BONDS && at[centerpoint].sb_parity[m]; m ++ ) {
                                    if ( at[centerpoint].sb_ord[m] == n ) {
                                        endpoint2 = NO_VERTEX;
                                        break;
                                    }
                                }
                                if ( endpoint2 == NO_VERTEX ) {
                                    continue;
                                }
                                pEndp2 = pBNS->vert + endpoint2;
                                ecp2   = ecp0;           /* e2C */
                                break;
                            }
                        }
                        if ( !ecp2 )
                            continue;
                        /* compare centerpoints */
                        if ( !pCentp_found ||
                             /* try to find carbons */
                             (!IS_C(endpoint2_found) && IS_C(endpoint2)) ||
                             (IS_C(endpoint2_found) && IS_C(endpoint2) &&
                             !IS_C(centerpoint_found) && IS_C(centerpoint)) ) {

                            pCentp_found = pCentp;
                            centerpoint_found = centerpoint;
                            endpoint2_found   = endpoint2;
                            ecp2_found   = ecp2;
                            ecp1_found   = ecp1;
                            ecp0_found   = pBNS->edge + e0;
                            break;
                        }
                    }
                }
            }
            /* decrement st_flow, st_cap on Endp2; decrement flow on ecp2; decrement st_flow on Centp */
            /* result: radicals on Endp0 and Centp => run BNS */
            if ( pCentp_found ) {
                Vertex     vPathStart, vPathEnd, v1, v2;
                int        nPathLen, nDeltaH, nDeltaCharge, nNumVisitedAtoms;
                int        delta;

                Vertex    vEndp0 = ecp0_found->neighbor1;
                Vertex    vEndp1 = ecp1_found->neighbor12 ^ centerpoint_found;
                Vertex    vEndp2 = ecp2_found->neighbor12 ^ centerpoint_found;
                BNS_EDGE *pe0 = ecp0_found;
                BNS_EDGE *pe1 = pBNS->edge + (pVA[vEndp1].nCMinusGroupEdge - 1);
                pEndp1 = pBNS->vert + vEndp1;
                pEndp2 = pBNS->vert + vEndp2;
                pCentp = pCentp_found;
                if ( !ChargeEdgeList.num_alloc ) {
                    for ( n = 0; n < pStruct->num_atoms; n ++ ) {
                        if ( (k = pVA[n].nCMinusGroupEdge)>= 0 && !pBNS->edge[k].forbidden &&
                             (ret = AddToEdgeList( &ChargeEdgeList, k, pStruct->num_atoms ) ) ) {
                            goto exit_function;
                        }
                        if ( (k = pVA[n].nCPlusGroupEdge)>= 0 && !pBNS->edge[k].forbidden &&
                             (ret = AddToEdgeList( &ChargeEdgeList, k, pStruct->num_atoms ) ) ) {
                            goto exit_function;
                        }
                    }
                }
                if ( !BondEdgeList.num_alloc ) {
                    for ( n = 0; n < pBNS->num_bonds; n ++ ) {
                        if ( (ret = AddToEdgeList( &BondEdgeList, n, pBNS->num_bonds ) ) ) {
                            goto exit_function;
                        }
                    }
                }
                /* fix all bonds and charges */
                SetForbiddenEdgeMask( pBNS, &ChargeEdgeList, forbidden_edge_mask );
                SetForbiddenEdgeMask( pBNS, &BondEdgeList, forbidden_edge_mask );
                /* prepare flow for testing */
                delta = 1;
                ecp2_found->flow     -= delta;
                pCentp->st_edge.flow -= delta;
                pEndp2->st_edge.flow -= delta;
                pBNS->tot_st_flow    -= 2*delta;
                /* unfix edges to be changed */
                pe0->forbidden        &= inv_forbidden_edge_mask;
                pe1->forbidden        &= inv_forbidden_edge_mask;
                ecp1_found->forbidden &= inv_forbidden_edge_mask;

                pBNS->tot_st_cap    += delta;

                v1 = vEndp0;
                v2 = centerpoint_found;
                ret2 = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                      &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

                if ( ret2 == 1 && ((vPathEnd == v1 && vPathStart == v2) ||
                                   (vPathEnd == v2 && vPathStart == v1)) && nDeltaCharge == 0 ) {
                    ret2 = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                    if ( ret2 > 0 ) {
                        num_fixes ++;
                        tot_num_fixes ++; /* #3 */
                        pCentp_found = NULL;
                    }
                } else {
                    /* roll back */
                    ecp2_found->flow     += delta;
                    pCentp->st_edge.flow += delta;
                    pEndp2->st_edge.flow += delta;
                    pBNS->tot_st_flow    += 2*delta;
                }
                RemoveForbiddenEdgeMask( pBNS, &ChargeEdgeList, forbidden_edge_mask );
                RemoveForbiddenEdgeMask( pBNS, &BondEdgeList, forbidden_edge_mask );
                if ( ret2 < 0 ) {
                    ret = ret2;
                    goto exit_function;
                }
                if ( !pCentp_found )
                    continue;
            }

        }
        if ( !num_fixes ) {
            break;
        }
    }

/************ again ***********************************************************/
    while ( pBNS->tot_st_cap > pBNS->tot_st_flow && pStruct->ti.num_t_groups ) {
        int iEndpoint = 0;
        num_fixes = 0;
        for ( itg = 0; itg < pStruct->ti.num_t_groups; iEndpoint += num_endpoints, itg ++ ) {
            pCentp_found=NULL;
            tgroup_number = pStruct->ti.t_group[itg].nGroupNumber;
            num_endpoints = pStruct->ti.t_group[itg].nNumEndpoints;
            for ( i = 0; i < num_endpoints; i ++ ) {
                endpoint0 = pStruct->ti.nEndpointAtomNumber[iEndpoint+i];
                pEndp0 = pBNS->vert + endpoint0;  /* taut endpoint vertex (possible location of mobile H */
                if ( pEndp0->st_edge.cap > pEndp0->st_edge.flow ) {
                    /* radical endpoint1 has been detected */
                    /* find a 1-3 centerpoint that has two or more endpoints */
                    /* connected to the t-group vertex by edges with flow>0 and */
                    /* to the centerpoint by edges with flow = 0 */
                    /* after that: (1) increment etg1 flow to eliminate radical */
                    /* (2) increment flow on one of the two other edges to the t-group */
                    /* (3) increment st_cap on the found centerpoint */
                    /* (4) rerun the BNS and re-create the structure */
                    break;
                }
            }
            /* 4th attempt */
            if ( i < num_endpoints ) {
                /* tautomeric endpoint found; traverse its t-group edges */
                pEndp2_found = NULL;
                for ( j = 0; j < pEndp0->num_adj_edges; j ++ ) {
                    ecp0 = pBNS->edge + pEndp0->iedge[j];
                    centerpoint = ecp0->neighbor12 ^ endpoint0;
                    if ( centerpoint >= pBNS->num_atoms || ecp0->flow || pStruct->endpoint[centerpoint] ) {
                        continue; /* ignore non-single bonds, orig. InChI endpoints, and fictitious atoms */
                    }
                    pCentp = pBNS->vert + centerpoint;
                    for ( k = 0; k < pCentp->num_adj_edges; k ++ ) {
                        ecp1 = pBNS->edge + pCentp->iedge[k];
                        endpoint1 = ecp1->neighbor12 ^ centerpoint;
                        if ( endpoint1 >= pBNS->num_atoms || !ecp1->flow || pStruct->endpoint[endpoint1] ) {
                            continue; /* ignore single bonds, orig. InChI endpoints, and fictitious atoms */
                        }
                        pEndp1 = pBNS->vert + endpoint1;
                        if ( endpoint1 == endpoint0 || pEndp1->st_edge.cap != pEndp1->st_edge.flow ) {
                            continue; /* ignore radicals */
                        }
                        if ( !pEndp2_found ||
                             /* try to find carbons */
                             (!IS_C(endpoint2_found) && IS_C(endpoint1)) ||
                             (IS_C(endpoint2_found) && IS_C(endpoint1) &&
                             !IS_C(centerpoint_found) && IS_C(centerpoint)) ) {
                            pEndp2_found      = pEndp1;
                            pCentp_found      = pCentp;
                            endpoint2_found   = endpoint1;
                            centerpoint_found = centerpoint;
                            ecp1_found        = ecp0;
                            ecp2_found        = ecp1;
                        }
                    }
                }
                if ( pEndp2_found ) {
                    /* move radical from pEndp0 to pEndp2 */
                    pEndp0->st_edge.flow ++;
                    ecp1_found->flow ++;
                    ecp2_found->flow --;
                    pEndp2_found->st_edge.flow --;
                    pEndp2_found = NULL;
                    pCentp_found = NULL;
                    num_fixes ++;    /* #4 */
                    tot_num_fixes ++;
                    continue;
                }
            }
        }
        if ( !num_fixes ) {
            break;
        }
    }


    ret = tot_num_fixes;
exit_function:
    AllocEdgeList( &ChargeEdgeList, EDGE_LIST_FREE);
    AllocEdgeList( &BondEdgeList, EDGE_LIST_FREE);

    pStruct->at = at;
    memcpy( at2, at, len_at*sizeof(at2[0]));
    return ret;
#undef IS_C
}
/************************************************************************************************/
/* move (+) charges to >N- and other centerpoints                                               */
int MoveChargeToMakeCenerpoints(BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct,
                              inp_ATOM *at, inp_ATOM *at2, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups,
                              int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask)
{
    int i, j, neigh, num_endpoints, tg_group=0, num_success;
    int ret2, ret, delta;
    int num_at = pStruct->num_atoms;
    int num_deleted_H = pStruct->num_deleted_H;
    int len_at = num_at + num_deleted_H;
    int inv_forbidden_edge_mask = ~forbidden_edge_mask;

    /* for RunBnsTestOnce */
    Vertex     vPathStart, vPathEnd;
    int        nPathLen, nDeltaH, nDeltaCharge, nNumVisitedAtoms;

    BNS_EDGE   *pEdgePlus, *pEdgeMinus;
    Vertex      v1p, v2p, v1m, v2m;
    BNS_VERTEX *pv1p, *pv2p, *pv1m, *pv2m;

    ret = 0;
    num_success = 0;
    /* to simplify, prepare new at[] from pBNS */
    memcpy( at2, at, len_at*sizeof(at2[0]));
    pStruct->at = at2;
    ret2 = CopyBnsToAtom( pStruct, pBNS, pVA, pTCGroups, 1 );
    pStruct->at = at;
    if ( ret2 < 0 ) {
        ret = ret2;
        goto exit_function;
    }

    for ( i = 0; i < num_at; i ++ ) {
        if ( pVA[i].cNumValenceElectrons != 4 && /* not C, Si, Ge */
             !pVA[i].cMetal && !pVA[i].nTautGroupEdge &&
             !at2[i].num_H && at2[i].valence >= 3 &&
             at2[i].valence == at2[i].chem_bonds_valence &&
             !at2[i].charge && pVA[i].nCPlusGroupEdge > 0 &&
             is_centerpoint_elem( at2[i].el_number ) ) {
            for ( j = 0, num_endpoints = 0; j < at2[i].valence; j ++ ) {
                neigh = at2[i].neighbor[j];
                if ( at2[neigh].endpoint ) {
                    if ( !num_endpoints ) {
                        tg_group = at2[neigh].endpoint;
                    } else
                    if ( tg_group != at2[neigh].endpoint ) {
                        break; /* not a centerpoint */
                    }
                    num_endpoints ++;
                }
            }
            if ( j == at2[i].valence && num_endpoints > 1 ) {
                /* found possible centerpoint */
                pEdgePlus  = pBNS->edge + (pVA[i].nCPlusGroupEdge-1);
                pEdgeMinus = (pVA[i].nCMinusGroupEdge > 0)? pBNS->edge + (pVA[i].nCMinusGroupEdge-1) : NULL;
                if ( pEdgePlus->flow + (pEdgeMinus? pEdgeMinus->flow : 0) != 1 ) {
                    continue;
                }
                v1p = pEdgePlus->neighbor1;
                v2p = pEdgePlus->neighbor12 ^ v1p;
                pv1p = pBNS->vert + v1p;
                pv2p = pBNS->vert + v2p;
                if ( pEdgeMinus ) {
                    v1m  = pEdgeMinus->neighbor1;
                    v2m  = pEdgeMinus->neighbor12 ^ v1m;
                    pv1m = pBNS->vert + v1m;
                    pv2m = pBNS->vert + v2m;
                } else {
                    v1m = NO_VERTEX;
                    v2m = NO_VERTEX;
                    pv1m = NULL;
                    pv2m = NULL;
                }
                ret = 0;
                /* set new flow to run BNS Search */
                if ( (delta = pEdgePlus->flow) ) {
                    /* positive charge <=> flow=0 on (=) edge */
                    pEdgePlus->flow -= delta;
                    pv1p->st_edge.flow -= delta;
                    pv2p->st_edge.flow -= delta;
                    pBNS->tot_st_flow  -= 2*delta;
                    pEdgePlus->forbidden |= forbidden_edge_mask;
                    if ( pEdgeMinus ) {
                        pEdgeMinus->forbidden |= forbidden_edge_mask;
                    }
                    ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                          &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );
                    if ( ret < 0 ) {
                        goto exit_function;
                    }
                    if ( ret == 1 && ((vPathEnd == v1p && vPathStart == v2p) ||
                                      (vPathEnd == v2p && vPathStart == v1p)) &&
                                      nDeltaCharge == -1 /* charge moving to this atom disappers*/ ) {
                        ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                        (*pnNumRunBNS) ++;
                        if ( ret < 0 ) {
                            goto exit_function;
                        } else
                        if ( ret == 1 ) {
                            *pnTotalDelta += ret;
                        } else {
                            ret = RI_ERR_PROGR;
                            goto exit_function;
                        }
                    } else {
                        ret = 0;
                        pEdgePlus->flow    += delta;
                        pv1p->st_edge.flow += delta;
                        pv2p->st_edge.flow += delta;
                        pBNS->tot_st_flow  += 2*delta;
                    }
                    pEdgePlus->forbidden &= inv_forbidden_edge_mask;
                    if ( pEdgeMinus ) {
                        pEdgeMinus->forbidden &= inv_forbidden_edge_mask;
                    }
                } else
                if ( pEdgeMinus && (delta == pEdgeMinus->flow) && pEdgePlus->flow == 0 ) {
                    /* positive charge <=> flow=0 on (=) edge and flow=0 on (-) edge */
                    pEdgeMinus->flow -= delta;
                    pv1m->st_edge.flow -= delta;
                    pv2m->st_edge.flow -= delta;
                    pBNS->tot_st_flow  -= 2*delta;
                    pEdgePlus->forbidden  |= forbidden_edge_mask;
                    pEdgeMinus->forbidden |= forbidden_edge_mask;
                    ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                          &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );
                    if ( ret < 0 ) {
                        goto exit_function;
                    }
                    if ( ret == 1 && ((vPathEnd == v1m && vPathStart == v2m) ||
                                      (vPathEnd == v2m && vPathStart == v1m)) &&
                                      nDeltaCharge == -1  /* charge moving to this atom disappers*/ ) {
                        ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                        (*pnNumRunBNS) ++;
                        if ( ret < 0 ) {
                            goto exit_function;
                        } else
                        if ( ret == 1 ) {
                            *pnTotalDelta += ret;
                        } else {
                            ret = RI_ERR_PROGR;
                            goto exit_function;
                        }
                    } else {
                        ret = 0;
                        pEdgeMinus->flow   += delta;
                        pv1m->st_edge.flow += delta;
                        pv2m->st_edge.flow += delta;
                        pBNS->tot_st_flow  += 2*delta;
                    }
                    pEdgePlus->forbidden  &= inv_forbidden_edge_mask;
                    pEdgeMinus->forbidden &= inv_forbidden_edge_mask;
                }
                if ( ret ) {
                    num_success ++;
                    memcpy( at2, at, len_at*sizeof(at2[0]));
                    pStruct->at = at2;
                    ret2 = CopyBnsToAtom( pStruct, pBNS, pVA, pTCGroups, 1 );
                    pStruct->at = at;
                    if ( ret2 < 0 ) {
                        ret = ret2;
                        goto exit_function;
                    }
                }
            }
        }
    }
    ret = num_success;
exit_function:
    return ret;
}
#endif

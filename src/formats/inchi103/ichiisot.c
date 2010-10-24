/*
 * International Chemical Identifier (InChI)
 * Version 1
 * Software version 1.03
 * May 9, 2010
 *
 * Originally developed at NIST
 * Modifications and additions by IUPAC and the InChI Trust
 *
 * The InChI library and programs are free software developed under the
 * auspices of the International Union of Pure and Applied Chemistry (IUPAC);
 * you can redistribute this software and/or modify it under the terms of 
 * the GNU Lesser General Public License as published by the Free Software 
 * Foundation:
 * http://www.opensource.org/licenses/lgpl-2.1.php
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mode.h"

#include "incomdef.h"
#include "extr_ct.h"
#include "ichitaut.h"
#include "ichicant.h"
#include "ichicomn.h"

#ifdef NEVER
/**********************************************************************************/
int unpack_iso_sort_key( AT_ISO_SORT_KEY iso_sort_key, S_CHAR *num_1H, S_CHAR *num_2H, S_CHAR *num_3H, S_CHAR *iso_atw_diff )
{
    int is_negative;
    AT_ISO_SORT_KEY HOnlyAtwPart;
    static const AT_ISO_SORT_KEY MultAtwDiff = AT_ISO_SORT_KEY_MULT*AT_ISO_SORT_KEY_MULT*AT_ISO_SORT_KEY_MULT;
    if ( !iso_sort_key ) {
        *num_1H = *num_2H = *num_3H = *iso_atw_diff = 0;
        return 0;
    } else
    if ( iso_sort_key < 0 ) {
        is_negative = 1;
        iso_sort_key = -iso_sort_key;
        HOnlyAtwPart = MultAtwDiff - iso_sort_key % MultAtwDiff;
        iso_sort_key += HOnlyAtwPart;
    } else {
        is_negative = 0;
        HOnlyAtwPart  = iso_sort_key % MultAtwDiff;
        iso_sort_key -= HOnlyAtwPart;
    }
    
    iso_sort_key /= MultAtwDiff;
    
    *num_1H = (S_CHAR)(HOnlyAtwPart % AT_ISO_SORT_KEY_MULT);
    HOnlyAtwPart /= AT_ISO_SORT_KEY_MULT;
    *num_2H = (S_CHAR)(HOnlyAtwPart % AT_ISO_SORT_KEY_MULT);
    HOnlyAtwPart /= AT_ISO_SORT_KEY_MULT;
    *num_3H = (S_CHAR)(HOnlyAtwPart % AT_ISO_SORT_KEY_MULT);
    
    *iso_atw_diff = (S_CHAR)(is_negative? -iso_sort_key : iso_sort_key);
    
    return 1;
}
#endif

/**********************************************************************************/
AT_ISO_SORT_KEY make_iso_sort_key( int iso_atw_diff, int num_1H, int num_2H, int num_3H)
{
    AT_ISO_SORT_KEY iso_sort_key = 0, mult=1;
    
    iso_sort_key += mult * num_1H;
    mult *= AT_ISO_SORT_KEY_MULT;
    iso_sort_key += mult * num_2H;
    mult *= AT_ISO_SORT_KEY_MULT;
    iso_sort_key += mult * num_3H;
    mult *= AT_ISO_SORT_KEY_MULT;
    iso_sort_key += mult * iso_atw_diff;
    return iso_sort_key;
}
/**********************************************************************************/
/*  set sp_ATOM isotopic sort keys */
int set_atom_iso_sort_keys(  int num_at, sp_ATOM *at, T_GROUP_INFO* t_group_info, int *bHasIsotopicInTautomerGroups )
{
    int             i, num_isotopic = 0, bMergedTgroup;
    AT_ISO_SORT_KEY iso_sort_key;
    T_GROUP        *t_group =
                       (t_group_info &&
                        t_group_info->t_group &&
                        t_group_info->num_t_groups > 0)? t_group_info->t_group : NULL;

    if ( bHasIsotopicInTautomerGroups )
        *bHasIsotopicInTautomerGroups = 0;
    for ( i = 0; i < num_at; i ++ ) {
        bMergedTgroup = (t_group_info && t_group_info->nIsotopicEndpointAtomNumber && (at[i].cFlags & AT_FLAG_ISO_H_POINT));
        if ( (!at[i].endpoint || !t_group) && !bMergedTgroup ) {
            iso_sort_key = make_iso_sort_key(at[i].iso_atw_diff, at[i].num_iso_H[0], at[i].num_iso_H[1], at[i].num_iso_H[2]);
        } else {
            /*  H isotopes go to the tautomer part of the CT (name) */
            /*  if (at[i].endpoint && t_group) ... */
            iso_sort_key = make_iso_sort_key(at[i].iso_atw_diff, 0, 0, 0);
            if ( bHasIsotopicInTautomerGroups )
                *bHasIsotopicInTautomerGroups += (at[i].num_iso_H[0] || at[i].num_iso_H[1] || at[i].num_iso_H[2] || bMergedTgroup);
        }
        at[i].iso_sort_key = iso_sort_key;
        num_isotopic += (iso_sort_key != 0);
    }
    return num_isotopic;
}


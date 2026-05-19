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
#include <ctype.h>
#include <math.h>

#include "mode.h"

#include "strutil.h"
#include "ichister.h"
#include "ichi_io.h"
#include "ichimain.h"

#include "bcf_s.h"

/* Added fix to remove_ion_pairs() -- 2010-03-17 DT */
#define FIX_P_IV_Plus_O_Minus

#ifndef COMPILE_ALL_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif

    /* Defined in ichisort.c, prototype in ichicomn.h */
int insertions_sort_AT_RANK( AT_RANK *base, int num );

#ifndef COMPILE_ALL_CPP
#ifdef __cplusplus
}
#endif
#endif


typedef struct tagTreeAtom {
    AT_NUMB    neighbor[MAXVAL];        /* positions (from 0) of the neighbors in the inp_ATOM array */
    S_CHAR     valence;                 /* number of bonds = number of neighbors */
    AT_NUMB    nRingSystem;
    AT_NUMB    nBlockSystem;
    S_CHAR     bCutVertex;
} tre_ATOM;


/* Local prototypes */

int cmp_iso_atw_diff_component_no( const void *a1, const void *a2 );
int cmp_components( const void *a1, const void *a2 );
/* int mark_one_struct_component( inp_ATOM* at,
int j,
AT_NUMB *mark,
AT_NUMB num_disconnected_components );
*/
INChI_Stereo *Alloc_INChI_Stereo( int num_at, int num_bonds );
int RemoveInpAtBond( inp_ATOM *at, int iat, int k );
int DisconnectInpAtBond( inp_ATOM *at,
                         AT_NUMB *nOldCompNumber,
                         int iat,
                         int neigh_ord );
int move_explicit_Hcation( inp_ATOM *at,
                           int num_at,
                           int iat,
                           int iat_H,
                           int bInAllComponents );
int DisconnectOneLigand( inp_ATOM *at,
                         AT_NUMB *nOldCompNumber,
                         S_CHAR *bMetal,
                         char *elnumber_Heteroat,
                         int num_halogens,
                         int num_atoms,
                         int iMetal,
                         int jLigand,
                         INCHI_MODE *bTautFlagsDone );
int bIsAmmoniumSalt( inp_ATOM *at,
                     int i,
                     int *piO,
                     int *pk,
                     S_CHAR *num_explicit_H );
int DisconnectAmmoniumSalt( inp_ATOM *at,
                            int i,
                            int iO,
                            int k,
                            S_CHAR *num_explicit_H );
/*int bIsMetalSalt( inp_ATOM *at, int i ); - moved to strutil,h */
int DisconnectMetalSalt( inp_ATOM *at, int i );
int bIsMetalToDisconnect( inp_ATOM *at, int i, int bCheckMetalValence );
int get_iat_number( int el_number );
int tot_unsat( int unsat[] );
int max_unsat( int unsat[] );
double dist3D( inp_ATOM *at1, inp_ATOM *at2 );
double dist2D( inp_ATOM *at1, inp_ATOM *at2 );
double dist_from_segm( double x, double y,
                       double x1, double y1,
                       double x2, double y2 );
int segments_intersect( double x11, double y11,
                        double x12, double y12, /* segment #1 */
                        double x21, double y21,
                        double x22, double y22 );
double GetMinDistDistribution( inp_ATOM *at,
                               int num_at,
                               int iat,
                               int iat_H,
                               int bInAllComponents,
                               double min_dist[],
                               int num_segm );
int nFindOneOM( inp_ATOM *at,
                int at_no,
                int ord_OM[],
                int num_OM );
int the_only_doublet_neigh( inp_ATOM *at, int i1, int *ineigh1, int *ineigh2 );
int fix_non_uniform_drawn_oxoanions( int num_atoms, inp_ATOM *at, int *num_changes );
int fix_non_uniform_drawn_amidiniums( int num_atoms, inp_ATOM *at, int *num_changes );


void add_bond_if_unseen( subgraf_pathfinder *spf,
                         int node0, int node,
                         int *nbonds, int **bonds );

/****************************************************************************/
#ifndef NUMH
#define NUM_ISO_H(AT,N) (AT[N].num_iso_H[0]+AT[N].num_iso_H[1]+AT[N].num_iso_H[2])
#define NUMH(AT,N)     (AT[N].num_H+NUM_ISO_H(AT,N))
#endif
/****************************************************************************/


/****************************************************************************/
int cmp_iso_atw_diff_component_no( const void *a1, const void *a2 )
{
    int ret = (int) ( (const inp_ATOM*) a1 )->iso_atw_diff - (int) ( (const inp_ATOM*) a2 )->iso_atw_diff;
    if (!ret) /*  make the sort stable */
    {
        ret = (int) ( (const inp_ATOM*) a1 )->component - (int) ( (const inp_ATOM*) a2 )->component;
    }

    return ret;
}


/****************************************************************************/
int the_only_doublet_neigh( inp_ATOM *at,
                            int i1,
                            int *ineigh1,
                            int *ineigh2 )
{
    int i, neigh1, num_rad1 = 0, num_rad2 = 0;

    inp_ATOM *a = at + i1, *b;
    if (RADICAL_DOUBLET != a->radical)
    {
        return -1;
    }
    for (i = 0; i < a->valence; i++)
    {
        b = at + ( (int) a->neighbor[i] ); /* djb-rwth: removing redundant code */
        if (RADICAL_DOUBLET == b->radical)
        {
            num_rad1++;
            *ineigh1 = i;
        }
    }

    if (1 == num_rad1)
    {
        a = at + ( neigh1 = (int) a->neighbor[*ineigh1] );
        for (i = 0; i < a->valence; i++)
        {
            b = at + (int) a->neighbor[i];
            if (RADICAL_DOUBLET == b->radical)
            {
                num_rad2++;
                *ineigh2 = i;
            }
        }

        if (1 == num_rad2)
        {
            return neigh1;
        }
    }

    return -1;
}


/****************************************************************************
Correct non-uniformly drawn oxoanions
****************************************************************************/
int fix_non_uniform_drawn_oxoanions( int num_atoms,
                                     inp_ATOM *at,
                                     int *num_changes )
{
    /* For central halogen, apply the following
    correction rules:

    O                             O(-)
    ||                            |
    O=Hal(-)=O         ===>        O=Hal=O
    ||                            ||
    O                             O

    (perchlorate, etc.)


    O                             O(-)
    ||                            |
    Hal(-)=O         ===>         Hal=O
    ||                            ||
    O                             O

    (chlorate, etc.)

    O                             O(-)
    ||                            |
    Hal(-)=O         ===>         Hal=O


    (chlorite, etc.)

    Hal(-)=O         ===>         Hal-O(-)


    (hypochlorite, etc.)



    For halcogenes (S, Se, Te)
    Y                               Y(-)
    ||                              |
    RnX(-)            ===>          RnX

    if:
    1) (X = S, Y = O) || (X = Se, Y = S, O) || (X = Te, Y = O, S, Se)
    2) valence of X exceeds 6, in initially drawn form


    So the following is corrected:

    O                               O(-)
    ||                              |
    O=S(-)-R            ===>        O=S-R
    ||                              ||
    O                               O

    Or

    O                              O-
    ||                             |
    F5Te(-)            ===>        F5Te

    but the following remains unchanged:


    O
    ||
    O=S(-)-R


    The central atom (of IUPAC Group 16-17) is shown as negative but it contains double bond(s)
    to terminal atom of greater electronegativity (of Group 16).
    The central atom is halcogen (S,Se,Te) in highest oxidation state or halogen.

    Fix:
    move negative charge to terminal atom and change double bond to a single one.

    Eligible central atom    Eligible terminal atom at double bond's end
    Cl              O
    Br              O
    I               O
    [At              S,Se,Te]
    S               O
    Se              O,S
    Te              O, S, Se

    Comments:
    1. Central atoms of Groups 13-15 are not considered.
    2. Pauling electronegativities are:
    F(3.98) > O(3.44) > Cl (3.16) > N (3.04) > Br(2.96) > I(2.66) > S(2.58) > Se(2.55) > At (2.2) > Te(2.1)


    */


    /* constants for element numbers */
    enum elems
    {
        dNone, dCl = 17, dBr = 35, dI = 53, dAt = 85, dO = 8,
        dS = 16, dSe = 34, dTe = 52, dP = 15, dC = 6, dN = 7
    };

    static U_CHAR  allowed_elnums_center_halogen[] = { dCl, dBr, dI, dAt };
    static U_CHAR  allowed_elnums_center_halcogen[] = { dS, dSe, dTe };

    int en_center;
    int i, j, k;

    for (i = 0; i < num_atoms; i++)
    {
        /* Find appropriate central atom. This center should be: ...*/

        /* charged exactly (-1) ... */
        if (at[i].charge != -1)
        {
            continue;
        }
        en_center = at[i].el_number;
        /*  from eligible element list ... */
        if (!memchr( allowed_elnums_center_halogen, en_center, sizeof( allowed_elnums_center_halogen ) ))
        {
            /* central atom is not not halogen; check if it is halcogen */
            if (memchr( allowed_elnums_center_halcogen, en_center, sizeof( allowed_elnums_center_halcogen ) ))
            {
                if (at[i].chem_bonds_valence < 7)
                {
                    /* central atom is anionic halcogen, but not in the highest oxidation state */
                    continue;
                }
            }
            else
            {
                continue;
            }
        }

        /* OK, found central halogen or eligible central halcogen. */

        /* non-radical... */
        if (at[i].radical && ( at[i].radical != RADICAL_SINGLET ))
        {
            continue;
        }

        /* Center found, now examine the adjacent terminals... */
        {
            int en_term, kk = 0, jj = 0, min_en = 999, iso = 0, min_iso = 999;
            jj = -1;
            for (k = 0; k < at[i].valence; k++)
            {
                j = at[i].neighbor[k];

                /* Terminal should be: ... */

                /* terminal... */
                if (at[j].valence != 1)
                {
                    continue;
                }
                /* double-bonded ... */
                if (at[i].bond_type[k] != BOND_TYPE_DOUBLE)
                {
                    continue;
                }
                /* zero-charged ... */
                if (at[j].charge != 0)
                {
                    continue;
                }
                /* non-radical */
                if (at[j].radical && ( at[j].radical != RADICAL_SINGLET ))
                {
                    continue;
                }
                /*  of eligible elements list ... */
                en_term = at[j].el_number;
                switch (en_term)
                {
                    case dO:    break;
                    case dS:    if (( en_center == dSe ) || ( en_center == dAt ) || ( en_center == dTe )) break;  continue;
                    case dSe:   if (( en_center == dAt ) || ( en_center == dTe )) break;  continue;
                    case dTe:   if (en_center == dAt) break; continue;
                    default:    continue;
                }

                /* From several candidates, select one with less el. number (==more electronegative). */
                if (en_term < min_en)
                {
                    min_en = en_term; kk = k; jj = j;
                    min_iso = at[j].iso_atw_diff > 0 ? at[i].iso_atw_diff - 1 : at[i].iso_atw_diff;
                    continue;
                }
                /* From same-element candidates, select one with less isotopic mass (arbitrary choice). */
                else if (en_term == min_en)
                {
                    iso = at[j].iso_atw_diff > 0 ? at[i].iso_atw_diff - 1 : at[i].iso_atw_diff;
                    if (iso < min_iso)
                    {
                        min_iso = iso; kk = k; jj = j; continue;
                    }
                }
            } /* end of checking nbrs. */

              /* If OK, apply changes. */
            if (jj >= 0)
            {
                at[i].charge = 0;
                at[jj].charge = -1;
                at[i].bond_type[kk] = BOND_TYPE_SINGLE;
                at[jj].bond_type[0] = BOND_TYPE_SINGLE;
                at[i].bond_stereo[kk] = at[jj].bond_stereo[0] = 0;
                at[i].chem_bonds_valence--;
                at[jj].chem_bonds_valence--;
                ( *num_changes )++;
            }
        }
    }  /* end of search for candidate centers. */

    return 0;
}


/****************************************************************************
Correct non-uniformly drawn amidinium cations.
****************************************************************************/
int fix_non_uniform_drawn_amidiniums( int num_atoms,
                                      inp_ATOM *at,
                                      int *num_changes )

{
    /* Amidines include carboxamidines RC(=NR)NR2,
    sulfinamidines RS(=NR)NR2 and phosphinamidines, R2P(=NR)NR2.


    NR                              NR
    |                               |
    R"-Y-NHR'         ===>          R"-Y=N(+)HR'
    (+)


    Y = C, S, P


    Fix:
    move positive charge to nitrogen and change single bond to a double one.


    Comment:
    Fix is applied only if at least one of R's at N is hydrogen
    (otherwise we have just a '+' delocalization which is already recognized).


    */

    /* constants for element numbers */
    enum elems
    {
        dNone, dCl = 17, dBr = 35, dI = 53, dAt = 85, dO = 8,
        dS = 16, dSe = 34, dTe = 52, dP = 15, dC = 6, dN = 7
    };

    static U_CHAR  allowed_elnums_center[] = { dC, dS, dP };
    int en_center;
    int i, j, k, jj, kk, k1;
    int mismatch = 0, nuH = 0, nuN = 0, nitrogens[MAXVAL];

    for (i = 0; i < num_atoms; i++)
    {
        /* Find appropriate central atom. This center should be: ...*/

        /* charged exactly (+1) ... */
        if (at[i].charge != 1)
        {
            continue;
        }
        en_center = at[i].el_number;
        /*  from eligible element list ... */
        if (!memchr( allowed_elnums_center, en_center, sizeof( allowed_elnums_center ) ))
        {
            continue;
        }
        /* has exactly 3 neighbours connected by single bonds*/
        if (at[i].valence != 3)
        {
            continue;
        }
        if (at[i].chem_bonds_valence != 3)
        {
            continue;
        }

        /* non-radical. */
        if (at[i].radical && ( at[i].radical != RADICAL_SINGLET ))
        {
            continue;
        }

        /* NB: center must have neutral neighbours, two of them are aliphatic N's of which at least one bears H. */
        mismatch = nuH = nuN = kk = 0; /* djb-rwth: removing redundant code */
        memset( nitrogens, 0, sizeof( nitrogens ) ); /* djb-rwth: memset_s C11/Annex K variant? */
        jj = -1;
        for (k = 0; k < at[i].valence; k++)
        {
            j = at[i].neighbor[k];

            if (at[j].charge != 0)
            {
                mismatch = 1;
                break;
            }
            if (at[j].el_number == dN)
            {
                if (( at[j].valence > 3 ) || ( at[j].chem_bonds_valence > 3 ))
                {
                    mismatch = 1;
                    break;
                }
                nuH += NUMH( at, j );
                nuN++;
                if (jj < 0)
                {
                    jj = j;
                    kk = k;
                }
            }
        }

        /* If OK, apply changes. */
        if (mismatch)
        {
            continue;
        }
        if (nuN != 2)
        {
            continue;
        }
        if (nuH < 1)
        {
            continue;
        }
        if (jj >= 0)
        {
            at[i].charge = 0;
            at[jj].charge = 1;
            at[i].bond_type[kk] = BOND_TYPE_DOUBLE;
            for (k1 = 0; k1 < at[jj].valence && i != at[jj].neighbor[k1]; k1++)
                ;
            at[jj].bond_type[k1] = BOND_TYPE_DOUBLE;
            at[i].chem_bonds_valence++;
            at[jj].chem_bonds_valence++;
            /* NB: do nothing with wedge stereo bonds (retain wedge) */

            ( *num_changes )++;
        }
    }  /* end of search for candidate centers. */

    return 0;
}


/****************************************************************************
Not used --
int FixAromaticOxygenAndSulfur( inp_ATOM *atom )
{
if ( !atom->elname[1] && (atom->elname[0]=='O' || atom->elname[0]=='S') &&
atom->valence==2 && !atom->charge && !atom->radical &&
atom->bond_type[0] + atom->bond_type[1] == 3 ) {
atom->charge = 1;
return 1; // fixed
}
return 0;
}
****************************************************************************/


/****************************************************************************/
int fix_odd_things( int num_atoms,
                    inp_ATOM *at,
                    int bFixBug,
                    int bFixNonUniformDraw )
{
    /* N;P;As;Sb;O;S;Se;Te;C;Si */
    static const U_CHAR en[] = {
        EL_NUMBER_N,
        EL_NUMBER_P,
        EL_NUMBER_AS,
        EL_NUMBER_SB,
        EL_NUMBER_O,
        EL_NUMBER_S,
        EL_NUMBER_SE,
        EL_NUMBER_TE
    };
    static int ne = sizeof(en)/sizeof(en[0]);

#define FIRST_NEIGHB2  4
#define FIRST_CENTER2  5

    int i1, i2, k1, k2, c = -1, num_changes = 0;
    /* djb-rwth: removing redundant variables */

    if (bFixNonUniformDraw)
    {
        int ret1; /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
        ret1 = fix_non_uniform_drawn_oxoanions( num_atoms, at, &num_changes ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
        ret1 = fix_non_uniform_drawn_amidiniums( num_atoms, at, &num_changes ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    }

    /* H(-)-X  -> H-X(-);  H(+)-X  -> H-X(+) */
    for (i1 = 0; i1 < num_atoms; i1++)
    {
        if (1 == at[i1].valence &&
             1 == abs( at[i1].charge ) &&
             ( 0 == at[i1].radical || RADICAL_SINGLET == at[i1].radical ) &&
             BOND_TYPE_SINGLE == at[i1].bond_type[0] &&
             EL_NUMBER_H == at[i1].el_number && EL_NUMBER_H != at[i2 = (int) at[i1].neighbor[0]].el_number &&
             !NUMH( at, i1 ) && !NUMH( at, i2 ))
        {
            at[i2].charge += at[i1].charge;
            at[i1].charge = 0;
        }
    }

    /* replace XHm(-)--Y==XHn(+) with XHm==Y--XHn, (n>=0 ,m>=0, X=N,P,As,Sb,O,S,Se,Te) */
    for (i1 = 0; i1 < num_atoms; i1++)
    {
        if (1 != at[i1].charge ||
             (at[i1].radical && RADICAL_SINGLET != at[i1].radical) ||
             at[i1].chem_bonds_valence == at[i1].valence ||
             !memchr( en, at[i1].el_number, ne ) ||
             get_el_valence( at[i1].el_number, at[i1].charge, 0 ) != at[i1].chem_bonds_valence + NUMH( at, i1 )) /* djb-rwth: addressing LLVM warning */
        {
            continue;
        }

        /* found a candidate at[i1] for X in XHn(+) */
        if (1 == at[i1].valence &&
             BOND_TYPE_DOUBLE == at[i1].bond_type[0])
        {
            c = (int) at[i1].neighbor[0];
            for (k2 = 0; k2 < at[c].valence; k2++)
            {
                i2 = at[c].neighbor[k2];
                if (1 == at[i2].valence &&
                     -1 == at[i2].charge  &&
                     at[i2].el_number == at[i1].el_number && /* exact match */
                     ( 0 == at[i2].radical || RADICAL_SINGLET == at[i2].radical ) &&
                     BOND_TYPE_SINGLE == at[i2].bond_type[0] &&
                     /*memchr(en, at[i2].el_number, ne) &&*/
                     get_el_valence( at[i2].el_number, at[i2].charge, 0 ) == at[i2].chem_bonds_valence + NUMH( at, i2 ))
                {
                    /* found both X(-) and X(+); change bonds and remove charges */
                    for (k1 = 0; k1 < at[c].valence && i1 != at[c].neighbor[k1]; k1++)
                        ;
                    at[i1].charge = at[i2].charge = 0;
                    at[i1].bond_type[0] = at[c].bond_type[k1] = BOND_TYPE_SINGLE;
                    at[i1].chem_bonds_valence--;
                    at[i2].bond_type[0] = at[c].bond_type[k2] = BOND_TYPE_DOUBLE;
                    at[i2].chem_bonds_valence++;
                    num_changes++;
                    break;
                }
            }
        }
        else
        {
            /* explicit H case: detect H-neighbors and Y */
            int ineigh, neigh, i1_c, i2_c, num_H_i1, num_H_i2;
            for (ineigh = 0, num_H_i1 = 0, i1_c = -1; ineigh < at[i1].valence; ineigh++)
            {
                neigh = at[i1].neighbor[ineigh];
                if (at[neigh].el_number == EL_NUMBER_H)
                {
                    if (at[neigh].chem_bonds_valence == 1 &&
                        ( 0 == at[neigh].radical || RADICAL_SINGLET == at[neigh].radical ))
                    {
                        num_H_i1++; /* found H-neighbor */
                    }
                    else
                    {
                        break;  /* wrong neighbor */
                    }
                }
                else if (at[i1].bond_type[ineigh] == BOND_TYPE_DOUBLE)
                {
                    /* found a candidate for Y; bond must be double */
                    i1_c = ineigh;
                    c = neigh;
                }
            }
            if (i1_c < 0 || num_H_i1 + 1 != at[i1].valence)
            {
                continue;
            }
            for (k2 = 0; k2 < at[c].valence; k2++)
            {
                i2 = at[c].neighbor[k2];
                if (-1 == at[i2].charge  &&
                     at[i2].el_number == at[i1].el_number && /* exact match */
                     ( 0 == at[i2].radical || RADICAL_SINGLET == at[i2].radical ) &&
                     get_el_valence( at[i2].el_number, at[i2].charge, 0 ) == at[i2].chem_bonds_valence + NUMH( at, i2 ))
                {
                    for (ineigh = 0, num_H_i2 = 0, i2_c = -1; ineigh < at[i2].valence; ineigh++)
                    {
                        neigh = at[i2].neighbor[ineigh];
                        if (at[neigh].el_number == EL_NUMBER_H)
                        {
                            if (at[neigh].chem_bonds_valence == 1 &&
                                ( 0 == at[neigh].radical || RADICAL_SINGLET == at[neigh].radical ))
                            {
                                num_H_i2++;  /* found H-neighbor */
                            }
                            else
                            {
                                break; /* wrong neighbor */
                            }
                        }
                        else
                        {
                            if (c == neigh && at[i2].bond_type[ineigh] == BOND_TYPE_SINGLE)
                            {
                                i2_c = ineigh; /* position of Y neighbor; bond must be single */
                            }
                            else
                            {
                                break;
                            }
                        }
                    }
                    if (num_H_i2 + ( i2_c >= 0 ) != at[i2].valence)
                    {
                        continue;
                    }
                    /* found both X(-) and X(+); change bonds and remove charges */
                    for (k1 = 0; k1 < at[c].valence && i1 != at[c].neighbor[k1]; k1++)
                        ;
                    if ((i1_c >= 0) && (i2_c >= 0)) /* djb-rwth: fixing coverity ID #499537 */
                    {
                        at[i1].charge = at[i2].charge = 0;
                        at[i1].bond_type[i1_c] = at[c].bond_type[k1] = BOND_TYPE_SINGLE;
                        at[i1].chem_bonds_valence--;
                        at[i2].bond_type[i2_c] = at[c].bond_type[k2] = BOND_TYPE_DOUBLE;
                        at[i2].chem_bonds_valence++;
                        num_changes++;
                        break;
                    }
                    else
                    {
                        continue;
                    }
                }
            } /* k2 */
        }
    }

    /* Replace

    X-                X        X=O,S,Se,Te -- terminal atoms  (NEIGHB2)
    \ |               \ ||
    >Y++    with      >Y        Y=S,Se,Te   -- central cation  (CENTER2)
    / |               / ||
    X-                X        Y valence=4, original Y bond valence = 4


    --- the following case of P is processed separately in remove_ion_pairs()
    --- therefire, it has been disabled here, see #ifndef FIX_P_IV_Plus_O_Minus -- 2010-03-17 DT

    X-                X        X=O,S,Se,Te -- terminal atoms  (NEIGHB2)
    \ |               \ ||
    >P+     with      >P
    / |               / |
    X-                X-       Y valence=4, original Y bond valence = 4

    */

    for (i1 = 0; i1 < num_atoms; i1++)
    {
        if (1 == at[i1].valence &&
             -1 == at[i1].charge &&
             ( 0 == at[i1].radical || RADICAL_SINGLET == at[i1].radical ) &&
             !NUMH( at, i1 ) &&
             BOND_TYPE_SINGLE == at[i1].bond_type[0] &&
             memchr( en + FIRST_NEIGHB2, at[i1].el_number, (long long)ne - FIRST_NEIGHB2 )) /* djb-rwth: cast operator added */
        {
            int charge, i;
            /* found a candidate for X */
            c = (int) at[i1].neighbor[0]; /* candidate for Y */
            if (( ( charge = 2 ) == at[c].charge && memchr( en + FIRST_CENTER2, at[c].el_number, (long long)ne - FIRST_CENTER2 ) /* djb-rwth: cast operator added */

#ifndef FIX_P_IV_Plus_O_Minus
                  || ( charge = 1 ) == at[c].charge && EL_NUMBER_P == at[c].el_number
#endif
                  ) &&
                 4 == at[c].valence &&
                 ( 0 == at[c].radical || RADICAL_SINGLET == at[c].radical ) &&
                 at[c].valence == at[c].chem_bonds_valence &&
                 !NUMH( at, c ))
            {
                ;  /* accept */
            }
            else
            {
                continue; /* ignore at[i1] */
            }
            for (k2 = 0; k2 < at[c].valence; k2++)
            {
                i2 = at[c].neighbor[k2];
                if (i2 == i1)
                {
                    continue;
                }
                if (1 == at[i2].valence &&
                     -1 == at[i2].charge  &&
                     memchr( en + FIRST_NEIGHB2, at[i2].el_number, (long long)ne - FIRST_NEIGHB2 ) && /* djb-rwth: cast operator added */
                     /*at[i2].el_number == at[i1].el_number &&*/ /* exact match */
                     ( 0 == at[i2].radical || RADICAL_SINGLET == at[i2].radical ) &&
                     !NUMH( at, i2 ) &&
                     BOND_TYPE_SINGLE == at[i2].bond_type[0])
                {
                    /* found both X(-) and X(-); change bonds and remove charges */
                    for (k1 = 0; k1 < at[c].valence && i1 != at[c].neighbor[k1]; k1++)
                    {
                        ;
                    }
                    for (i = 0; i < charge; i++)
                    {
                        /* in case of P it does not matter which X atom is neutralized
                        because of tautomerism. However, neutral central atom is important
                        for the neutralization of the components */
                        switch (i)
                        {
                            case 0:
                                at[i1].charge++; /* = 0; changed 2010-03-17 DT*/
                                at[i1].bond_type[0] = at[c].bond_type[k1] = BOND_TYPE_DOUBLE;
                                at[i1].bond_stereo[0] = at[c].bond_stereo[k1] = 0;
                                at[i1].chem_bonds_valence++;
                                at[c].chem_bonds_valence++;
                                if (bFixBug) at[c].charge--; /* added 2010-03-17 DT*/
                                num_changes++;
                                break;
                            case 1:
                                at[i2].charge++; /*= 0; changed 2010-03-17 DT*/
                                at[i2].bond_type[0] = at[c].bond_type[k2] = BOND_TYPE_DOUBLE;
                                at[i2].bond_stereo[0] = at[c].bond_stereo[k2] = 0;
                                at[i2].chem_bonds_valence++;
                                at[c].chem_bonds_valence++;
                                if (bFixBug) at[c].charge--; /* added 2010-03-17 DT */
                                num_changes++;
                                break;
                        }
                    }

                    /*   -- removed -- 2010-03-17 DT
                    #if ( FIX_ODD_THINGS_REM_Plus_BUG == 1 )
                    at[c].charge -= charge;
                    #else
                    if ( bFixBug )
                    {
                    at[c].charge -= charge;
                    }
                    #endif
                    */
                    break;
                }
            }
        }
    }

    /* A(doublet)-B(doublet) -> A=B  (A and B have no other doublet neighbors) */
    /* A(doublet)=B(doublet) -> A#B  (A and B have no other doublet neighbors) */
    for (i1 = 0; i1 < num_atoms; i1++)
    {
        if (RADICAL_DOUBLET == at[i1].radical &&
             0 <= ( i2 = the_only_doublet_neigh( at, i1, &k1, &k2 ) ))
        {
            if (at[i1].bond_type[k1] <= BOND_TYPE_DOUBLE)
            {
                at[i1].bond_type[k1] ++;
                at[i1].chem_bonds_valence++;
                at[i2].bond_type[k2] ++;
                at[i2].chem_bonds_valence++;
                at[i1].radical = 0;
                at[i2].radical = 0;
            }
        }
    }

#if ( REMOVE_ION_PAIRS_EARLY == 1 )
    num_changes += remove_ion_pairs( num_atoms, at );
#endif

    return num_changes;
}


/****************************************************************************/
int post_fix_odd_things( int num_atoms, inp_ATOM *at )
{
    int num_changes = 0;
    /* currently does nothing */
    return
        num_changes;
}


/****************************************************************************/
int nFindOneOM( inp_ATOM *at, int at_no, int ord_OM[], int num_OM )
{
    int i, n_OM, best_value, cur_value, diff; /* djb-rwth: removing redundant variables */
    int num_best;

    if (1 == num_OM)
    {
        return ord_OM[0];
    }
    if (1 > num_OM)
    {
        return -1;
    }

    /* select neighbors with min. number of bonds */
    num_best = 1;
    n_OM = (int) at[at_no].neighbor[ord_OM[0]];
    best_value = (int) at[n_OM].valence;
    /* compare number of bonds; move indexes of the best neighbors to the first elements of ord_OM[] */
    for (i = 1; i < num_OM; i++)
    {
        n_OM = at[at_no].neighbor[ord_OM[i]];
        cur_value = (int) at[n_OM].valence;
        diff = cur_value - best_value;
        if (diff < 0)
        {
            /* djb-rwth: removing redundant code */
            best_value = cur_value;
            ord_OM[0] = ord_OM[i];
            num_best = 1;
        }
        else if (diff == 0)
        {    /* was '=', pointed by WDI */
            ord_OM[num_best++] = ord_OM[i];
        }
    }
    num_OM = num_best;
    if (1 == num_OM)
    {
        return ord_OM[0];
    }

    /* select neighbors with min. periodic numbers */
    num_best = 1;
    n_OM = (int) at[at_no].neighbor[ord_OM[0]];
    best_value = (int) at[n_OM].el_number;

    /* compare periodic numbers; move indexes of the best neighbors to the first elements of ord_OM[] */
    for (i = 1; i < num_OM; i++)
    {
        n_OM = at[at_no].neighbor[ord_OM[i]];
        cur_value = (int) at[n_OM].el_number;
        diff = cur_value - best_value;
        if (diff < 0)
        {
            /* djb-rwth: removing redundant code */
            best_value = cur_value;
            ord_OM[0] = ord_OM[i];
            num_best = 1;
        }
        else if (diff == 0)
        {    /* was '=', pointed by WDI */
            ord_OM[num_best++] = ord_OM[i];
        }
    }
    num_OM = num_best;
    if (1 == num_OM)
    {
        return ord_OM[0];
    }

    /* if neighbors are not terminal atoms then reject */
    if (1 < at[n_OM].valence)
    {
        return -1;
    }

    /* if neighbors are terminal atoms then the one without isotope or with lightest isotope */
    num_best = 1;
    n_OM = (int) at[at_no].neighbor[ord_OM[0]];
    best_value = (int) at[n_OM].iso_atw_diff;

    /* compare periodic numbers; move indexes of the best neighbors to the first elements of ord_OM[] */
    for (i = 1; i < num_OM; i++)
    {
        n_OM = at[at_no].neighbor[ord_OM[i]];
        cur_value = (int) at[n_OM].el_number;
        diff = cur_value - best_value;
        if (( !cur_value && best_value ) || diff < 0)
        {
            /* djb-rwth: removing redundant code */
            best_value = cur_value;
            ord_OM[0] = ord_OM[i];
            num_best = 1;
        }
        else if (diff == 0)
        {
            /* was '=', pointed by WDI */
            ord_OM[num_best++] = ord_OM[i];
        }
    }

    num_OM = num_best;
    if (1 == num_OM)
    {
        return ord_OM[0];
    }

    /* return any */
    return ord_OM[0];
}


/****************************************************************************
NB:
the bonds are fixed in fix_special_bonds()
****************************************************************************/
int remove_ion_pairs( int num_atoms, inp_ATOM *at )
{
    int num_changes = 0;
#define MAX_NEIGH 6

    int i, n, n2, i1, i2, i3, i4, type, chrg;
    int num_C_II = 0, num_C_plus = 0, num_C_minus = 0, num_N_plus = 0, num_N_minus = 0, num_O_plus = 0, num_O_minus = 0, num_All;

#ifdef FIX_P_IV_Plus_O_Minus
    int num_P_IV_plus = 0; /* added 2010-03-17 DT */
#endif

    inp_ATOM *a;
    /****** count candidates ********/
    for (i = 0, a = at; i < num_atoms; i++, a++)
    {
        if (1 == ( chrg = a->charge ) || -1 == chrg)
        {
            switch (ion_el_group( a->el_number ))
            {
                case EL_NUMBER_C:
                    if (chrg > 0)
                    {
                        num_C_plus++;
                    }
                    else
                    {
                        num_C_minus++;
                    }
                    break;
                case EL_NUMBER_O:
                    if (chrg > 0)
                    {
                        num_O_plus++;
                    }
                    else
                    {
                        num_O_minus++;
                    }
                    break;
                case EL_NUMBER_N:
                    if (chrg > 0)
                    {
                        num_N_plus++;
                    }
                    else
                    {
                        num_N_minus++;
                    }
#ifdef FIX_P_IV_Plus_O_Minus
                    num_P_IV_plus += a->el_number != EL_NUMBER_N && 
                                     chrg == 1 &&
                                     a->valence == 4 && 
                                     a->chem_bonds_valence == 4; /* added 2010-03-17 DT */
#endif 
                    break;                
            }
        }
        else if (!chrg && a->chem_bonds_valence + NUMH( a, 0 ) == 2 &&
                  get_el_valence( a->el_number, 0, 0 ) == 4 &&
                  ion_el_group( a->el_number ) == EL_NUMBER_C)
        {
            num_C_II++;
        }
    }

    num_All = num_C_II + num_C_plus + num_C_minus + num_N_plus + num_N_minus + num_O_plus + num_O_minus;

    /* do not add num_P_IV_plus ! -- 2010-03-17 DT */
    if (!num_All)
    {
        return 0;
    }

    /**************************************************************************/
    /*************************** Terminal ion pairs ***************************/
    /**************************************************************************/

    /*-------------------------------------------------------------------------
    Pair type 1            N=N,P,As,Sb; O=O,S,Se,Te
    ===========

    X              X     if X is another -O(-) then neutralize O(-)
    |              |     that has the smallest periodic table number
    O=N(+)-O(-) => O=N=O
    i    n
    --------------------------------------------------------------------------*/

    for (type = 1; type <= 18; type++)
    {
        if (( !type || 1 == type ))
        {
            for (i = 0; i < num_atoms && 0 < num_N_plus && 0 < num_O_minus; i++)
            {
                if (1 == at[i].charge && 3 == nNoMetalNumBonds( at, i ) &&
                     4 == nNoMetalBondsValence( at, i ) &&
                     ion_el_group( at[i].el_number ) == EL_NUMBER_N)
                {
                    int num_OM = 0, ord_OM[3]; /* -O(-) */
                    int num_O = 0; /* =O    */
                    int num_O_other = 0;
                    for (i1 = 0; i1 < at[i].valence; i1++)
                    {
                        n = at[i].neighbor[i1];
                        if (1 == nNoMetalNumBonds( at, n ) && 0 == num_of_H( at, n ) &&
                            ion_el_group( at[n].el_number) == EL_NUMBER_O) /* djb-rwth: ignoring LLVM warning: variable used */
                        {
                            if (BOND_TYPE_SINGLE == at[i].bond_type[i1] &&
                                 -1 == at[n].charge)
                            {
                                ord_OM[num_OM++] = i1;
                            }
                            else if (BOND_TYPE_DOUBLE == at[n].bond_type[0] &&
                                      0 == at[n].charge)
                            {
                                num_O++;
                            }
                            else
                            {
                                num_O_other++;
                            }
                        }
                    }
                    if (num_OM > 0 && num_O > 0 && !num_O_other &&
                         0 <= ( i1 = nFindOneOM( at, i, ord_OM, num_OM ) ))
                    {
                        /* remove charges and increase bond order */
                        n = at[i].neighbor[i1];
                        i2 = (int) ( is_in_the_list( at[n].neighbor, (AT_NUMB) i, at[n].valence ) - at[n].neighbor );
                        at[i].bond_type[i1] ++;
                        at[n].bond_type[i2] ++;
                        at[i].chem_bonds_valence++;
                        at[n].chem_bonds_valence++;
                        at[i].charge--;
                        at[n].charge++;
                        at[i].radical = 0;
                        at[n].radical = 0;
                        num_changes++;
                        num_N_plus--;
                        num_O_minus--;
                        num_All -= 2;
                    }
                }
            }

#ifdef FIX_P_IV_Plus_O_Minus

            /*-------------------------------------------------------------------------
            Pair type 1a           P=P,As,Sb; O=O,S,Se,Te  -- added 2010-03-17
            =============

            X              X     if X, Y, or Z is another -O(-) then neutralize O(-)
            |              |     that has the smallest periodic table number
            Y-P(+)-O(-) => Y-P=O
            |i   n         |
            Z              Z

            --------------------------------------------------------------------------*/

            for (i = 0; i < num_atoms && 0 < num_P_IV_plus /*&& 0 < num_N_plus*/ && 0 < num_O_minus; i++)
            {
                if (1 == at[i].charge && 4 == nNoMetalNumBonds( at, i ) &&
                     4 == nNoMetalBondsValence( at, i ) &&
                     at[i].el_number != EL_NUMBER_N && ion_el_group( at[i].el_number ) == EL_NUMBER_N)
                {
                    int num_OM = 0, ord_OM[4]; /* -O(-) */
                                               /*int num_O  = 0;*/ /* =O    */
                    /* djb-rwth: removing redundant variables */
                    for (i1 = 0; i1 < at[i].valence; i1++)
                    {
                        n = at[i].neighbor[i1];
                        if (1 == nNoMetalNumBonds( at, n ) && 0 == num_of_H( at, n ) &&
                            ion_el_group( at[n].el_number) == EL_NUMBER_O) /* djb-rwth: ignoring LLVM warning: variable used */
                        {
                            if (BOND_TYPE_SINGLE == at[i].bond_type[i1] &&
                                 -1 == at[n].charge)
                            {
                                ord_OM[num_OM++] = i1;
                                /*
                                }
                                if ( BOND_TYPE_DOUBLE == at[n].bond_type[0] &&
                                0                == at[n].charge       ) {
                                num_O ++;
                                */
                            }
                            /* djb-rwth: removing redundant code */
                        }
                    }
                    if (num_OM > 0 /*&& num_O > 0 && !num_O_other*/ &&
                         0 <= ( i1 = nFindOneOM( at, i, ord_OM, num_OM ) ))
                    {
                        /* remove charges and increase bond order */
                        n = at[i].neighbor[i1];
                        i2 = is_in_the_list( at[n].neighbor, (AT_NUMB) i, at[n].valence ) - at[n].neighbor;
                        at[i].bond_type[i1] ++;
                        at[n].bond_type[i2] ++;
                        at[i].chem_bonds_valence++;
                        at[n].chem_bonds_valence++;
                        at[i].charge--;
                        at[n].charge++;
                        at[i].radical = 0;
                        at[n].radical = 0;
                        num_changes++;
                        num_N_plus--;
                        num_O_minus--;
                        num_P_IV_plus--;
                        num_All -= 2;
                    }
                }
            }
#endif /* FIX_P_IV_Plus_O_Minus */
        }


        /*-------------------------------------------------------------------------
        Terminal pair types: 2,3,4,5,6,7,8,9   N=N,P,As,Sb; O=O,S,Se,Te; C=C,Si
        ====================================
        type #
        2 2:  O=N-C(II)-     => O=N#C-      N=N,P,As,Sb; O=O,S,Se,Te; C=C,Si
        3 9:  O=O(+)-C(-)(III) => O=O=C(IV)
        4 3:  O(-)-N(+)(IV)  => O=N(V)  (input structure has at least 1 double bond)
        5 4:  O(-)-O(+)(III) => O=O(IV)
        6 8:  O(-)-O-C(+)(III) => O=O=C(IV)
        7 5:  N(-)=N(+)(IV)  => N#N(V)    allow terminal H on N(-)
        8 6:  N(-)=O(+)(III) => N#O-
        9 7:  N(-)=C(+)(III) => N#C-
        --------------------------------------------------------------------------*/

        if (!type || (2 <= type && type <= 9)) /* djb-rwth: addressing LLVM warning */
        {
            for (i = 0; i < num_atoms && 0 < num_All; i++)
            {
                if (0 == at[i].charge && 1 == nNoMetalNumBonds( at, i ) && 2 == nNoMetalBondsValence( at, i ) &&
                     0 == num_of_H( at, i ) &&
                     ion_el_group( at[i].el_number ) == EL_NUMBER_O &&
                     0 <= ( i1 = nNoMetalNeighIndex( at, i ) ) &&
                     at[i].bond_type[i1] <= BOND_TYPE_TRIPLE)
                {
                    /* terminal O= */
                    n = at[i].neighbor[i1];
                    if (( !type || type == 2 ) && 0 < num_C_II)
                    {
                        /* avoid alternating bonds */
                        if (0 == at[n].charge &&
                             2 == nNoMetalNumBonds( at, n ) && 3 == nNoMetalBondsValence( at, n ) &&
                             0 == num_of_H( at, n ) &&
                             ion_el_group( at[n].el_number ) == EL_NUMBER_N &&
                             0 <= ( i2 = nNoMetalOtherNeighIndex( at, n, i ) ) &&
                             at[n].bond_type[i2] <= BOND_TYPE_TRIPLE)
                        {
                            /* i2 = index of opposite to at[i] neighbor of at[n] */
                            /*i2 = (at[n].neighbor[0] == i);*/
                            n2 = at[n].neighbor[i2];
                            if (0 == at[n2].charge &&
                                 2 == at[n2].valence && 2 == at[n2].chem_bonds_valence &&
                                 0 == num_of_H( at, n2 ) &&
                                 ion_el_group( at[n2].el_number ) == EL_NUMBER_C)
                            {
                                /*       i n n2     */
                                /* found O=N-C(II)- */
                                /* convert O=N-C(II)-     => O=N#C- */

                                i3 = ( at[n2].neighbor[0] != n ); /* index of at[n] neighbor of n2 */
                                at[n].chem_bonds_valence = 5; /* N */
                                at[n2].chem_bonds_valence = 4; /* C */
                                at[n].bond_type[i2] = BOND_TYPE_TRIPLE;
                                at[n2].bond_type[i3] = BOND_TYPE_TRIPLE;
                                at[n2].radical = 0;
                                num_changes++;
                                num_C_II--;
                                num_All--;
                                continue;
                            }
                        }
                    }

                    if (( !type || type == 3 ) && 0 < num_O_plus && 0 < num_C_minus)
                    {
                        if (1 == at[n].charge && 2 == nNoMetalNumBonds( at, n ) && 3 == nNoMetalBondsValence( at, n ) &&
                             0 == num_of_H( at, n ) &&
                             ion_el_group( at[n].el_number ) == EL_NUMBER_O &&
                             0 <= ( i2 = nNoMetalOtherNeighIndex( at, n, i ) ) &&
                             at[n].bond_type[i2] <= BOND_TYPE_TRIPLE)
                        {
                            /* found O=O(+)- */
                            /* i2 = index of opposite to at[i] neighbor of at[n] */
                            /*i2 = (at[n].neighbor[0] == i);*/
                            n2 = at[n].neighbor[i2];
                            if (-1 == at[n2].charge && 3 >= nNoMetalNumBonds( at, n2 ) && 3 == nNoMetalBondsValence( at, n2 ) + NUMH( at, n2 ) &&
                                 ion_el_group( at[n2].el_number ) == EL_NUMBER_C)
                            {
                                /*             i n    n2        */
                                /* found found O=O(+)-C(-)(III) */
                                /* convert O=O(+)-C(-)(III)     => O=O=C(IV) */
                                i3 = ( at[n2].neighbor[0] != n ); /* index of at[n] neighbor of n2 */
                                at[n].charge--;
                                at[n2].charge++;
                                at[n].chem_bonds_valence += 1; /* =O- => =O= */
                                at[n2].chem_bonds_valence += 1; /* -C  => =C  */
                                at[n].bond_type[i2] = BOND_TYPE_DOUBLE;
                                at[n2].bond_type[i3] = BOND_TYPE_DOUBLE;
                                num_changes++;
                                num_O_plus--;
                                num_C_minus--;
                                num_All -= 2;
                                continue;
                            }
                        }
                    }
                }

                else if (-1 == at[i].charge &&
                          0 < num_O_minus + num_N_minus &&
                          0 < num_N_plus + num_O_plus + num_C_plus &&
                          1 == nNoMetalNumBonds( at, i ) && 1 == nNoMetalBondsValence( at, i ) &&
                          0 == num_of_H( at, i ) &&
                          ion_el_group( at[i].el_number ) == EL_NUMBER_O &&
                          0 <= ( i1 = nNoMetalNeighIndex( at, i ) ) &&
                          at[i].bond_type[i1] <= BOND_TYPE_TRIPLE)
                {
                    /* terminal O(-)- */
                    n = at[i].neighbor[i1];

                    if (( !type || type == 4 ) && 0 < num_O_minus && 0 < num_N_plus && /* O(-)-N(+)(IV) */
                         1 == at[n].charge && 3 >= nNoMetalNumBonds( at, n ) && 4 == nNoMetalBondsValence( at, n ) &&
                         0 == num_of_H( at, n ) &&
                         ion_el_group( at[n].el_number ) == EL_NUMBER_N /* except >O(+)- */
                         )
                    {
                        /* found O(-)-N(+)(IV) */
                        /* convert O(-)-N(+)(IV)     => O=N(V)  */

                        i2 = (int) ( is_in_the_list( at[n].neighbor, (AT_NUMB) i, at[n].valence ) - at[n].neighbor ); /* index of at[i] neighbor of at[n] */
                        at[i].charge++;
                        at[n].charge--;
                        at[i].chem_bonds_valence++;
                        at[n].chem_bonds_valence++;
                        at[i].bond_type[i1] ++;
                        at[n].bond_type[i2] ++;
                        num_changes++;
                        num_O_minus--;
                        num_N_plus--;
                        num_All -= 2;
                        continue;
                    }

                    if (( !type || type == 5 ) && 0 < num_O_minus && 0 < num_O_plus &&/* O(-)-O(+)(III) */
                         1 == at[n].charge && 3 >= nNoMetalNumBonds( at, n ) && 3 == nNoMetalBondsValence( at, n ) &&
                         0 == num_of_H( at, n ) &&
                         ion_el_group( at[n].el_number ) == EL_NUMBER_O /* except >O(+)- */
                         )
                    {
                        /* found  O(+)(III) */
                        /* convert O(-)-O(+)(III)    => O=O(IV) */

                        i2 = (int) ( is_in_the_list( at[n].neighbor, (AT_NUMB) i, at[n].valence ) - at[n].neighbor ); /* index of at[i] neighbor of at[n] */
                        at[i].charge++;
                        at[n].charge--;
                        at[i].chem_bonds_valence++;
                        at[n].chem_bonds_valence++;
                        at[i].bond_type[i1] ++;
                        at[n].bond_type[i2] ++;
                        num_changes++;
                        num_O_minus--;
                        num_O_plus--;
                        num_All -= 2;
                        continue;
                    }

                    /* i    n n2        */
                    if (( !type || type == 6 ) && /* O(-)-O-C(+)(III) */
                         0 < num_O_minus && 0 < num_C_plus &&
                         0 == at[n].charge && 2 == nNoMetalNumBonds( at, n ) && 2 == nNoMetalBondsValence( at, n ) &&
                         0 == num_of_H( at, n ) &&
                         ion_el_group( at[n].el_number ) == EL_NUMBER_O &&
                         0 <= ( i2 = nNoMetalOtherNeighIndex( at, n, i ) ) &&
                         at[n].bond_type[i2] <= BOND_TYPE_TRIPLE)
                    {
                        /* found O(-)-O- */
                        /* i2 = index of opposite to at[i] neighbor of at[n] */
                        /*i2 = (at[n].neighbor[0] == i);*/
                        n2 = at[n].neighbor[i2];
                        if (1 == at[n2].charge && 3 >= nNoMetalNumBonds( at, n2 ) &&
                             3 == nNoMetalBondsValence( at, n2 ) + NUMH( at, n2 ) &&
                             ion_el_group( at[n2].el_number ) == EL_NUMBER_C)
                        {
                            /*       i    n n2  */
                            /* found O(-)-O-C(+)(III) */
                            /* convert O(-)-O-C(+)(III)     => O=O=C(IV) */
                            /*i3 = (at[n2].neighbor[0] != n);*/ /* i3 = index of at[n] neighbor of at[n2] */
                            i3 = (int) ( is_in_the_list( at[n2].neighbor, (AT_NUMB) n, at[n2].valence ) - at[n2].neighbor );
                            /*i4 = index of at[i] in the adjacency list of at[n] */
                            i4 = (int) ( is_in_the_list( at[n].neighbor, (AT_NUMB) i, at[n].valence ) - at[n].neighbor );
                            at[i].charge++;
                            at[n2].charge--;
                            at[i].chem_bonds_valence += 1; /* O-  => O=  */
                            at[n].chem_bonds_valence += 2; /* -O- => =O= */
                            at[n2].chem_bonds_valence += 1; /* -C  => =C  */
                            at[i].bond_type[i1] = BOND_TYPE_DOUBLE;
                            at[n].bond_type[i4] = BOND_TYPE_DOUBLE;
                            at[n].bond_type[i2] = BOND_TYPE_DOUBLE;
                            at[n2].bond_type[i3] = BOND_TYPE_DOUBLE;
                            num_changes++;
                            num_O_minus--;
                            num_C_plus--;
                            num_All -= 2;
                            continue;
                        }
                    }
                }
                else if (-1 == at[i].charge && 0 < num_N_minus && 0 < num_N_plus + num_O_plus + num_C_plus &&
                          1 == nNoMetalNumBonds( at, i ) && 2 == nNoMetalBondsValence( at, i ) + NUMH( at, i ) &&
                          /*0 == num_of_H( at, i ) &&*/
                          ion_el_group( at[i].el_number ) == EL_NUMBER_N &&
                          0 <= ( i1 = nNoMetalNeighIndex( at, i ) ) &&
                          at[i].bond_type[i1] <= BOND_TYPE_TRIPLE)
                {
                    /* terminal N(-)= */
                    n = at[i].neighbor[i1 = 0];
                    if (( !type || type == 7 ) && 0 < num_N_plus && /* N(-)=N(+)(IV) */
                         1 == at[n].charge && 3 >= nNoMetalNumBonds( at, n ) && 4 == nNoMetalBondsValence( at, n ) &&
                         0 == num_of_H( at, n ) &&
                         ion_el_group( at[n].el_number ) == EL_NUMBER_N)
                    {
                        /* found N(-)-N(+)(IV) */
                        /* convert N(-)=N(+)(IV)     => N#N(V)  */

                        i2 = (int) ( is_in_the_list( at[n].neighbor, (AT_NUMB) i, at[n].valence ) - at[n].neighbor ); /* index of at[i] neighbor of at[n] */
                        at[i].charge++;
                        at[n].charge--;
                        at[i].chem_bonds_valence++;
                        at[n].chem_bonds_valence++;
                        at[i].bond_type[i1] ++;
                        at[n].bond_type[i2] ++;
                        num_changes++;
                        num_N_minus--;
                        num_N_plus--;
                        num_All -= 2;
                        continue;
                    }
                    if (( !type || type == 8 ) && 0 < num_O_plus && /* N(-)=O(+)(III) */
                         1 == at[n].charge && 2 == nNoMetalNumBonds( at, n ) && 3 == nNoMetalBondsValence( at, n ) &&
                         0 == num_of_H( at, n ) &&
                         ion_el_group( at[n].el_number ) == EL_NUMBER_O)
                    {
                        /* found N(-)-O(+)(III) */
                        /* convert N(-)=O(+)(III)    => N#O(IV)- */
                        i2 = (int) ( is_in_the_list( at[n].neighbor, (AT_NUMB) i, at[n].valence ) - at[n].neighbor ); /* index of at[i] neighbor of at[n] */
                        at[i].charge++;
                        at[n].charge--;
                        at[i].chem_bonds_valence++;
                        at[n].chem_bonds_valence++;
                        at[i].bond_type[i1] ++;
                        at[n].bond_type[i2] ++;
                        num_changes++;
                        num_N_minus--;
                        num_O_plus--;
                        num_All -= 2;
                        continue;
                    }
                    if (( !type || type == 9 ) && 0 < num_C_plus && /* N(-)=C(+)(III) */
                         1 == at[n].charge && 2 == at[n].valence && 3 == at[n].chem_bonds_valence &&
                         0 == num_of_H( at, n ) &&
                         ion_el_group( at[n].el_number ) == EL_NUMBER_C)
                    {
                        /* found N(-)=C(+)(III) */
                        /* convert N(-)=C(+)(III)    => N#C(IV)- */

                        i2 = (int) ( is_in_the_list( at[n].neighbor, (AT_NUMB) i, at[n].valence ) - at[n].neighbor ); /* index of at[i] neighbor of at[n] */
                        at[i].charge++;
                        at[n].charge--;
                        at[i].chem_bonds_valence++;
                        at[n].chem_bonds_valence++;
                        at[i].bond_type[i1] ++;
                        at[n].bond_type[i2] ++;
                        num_changes++;
                        num_N_minus--;
                        num_C_plus--;
                        num_All -= 2;
                        continue;
                    }
                }
            }
        }

        /**************************************************************************/
        /*********************** NON-Terminal ion pairs ***************************/
        /**************************************************************************/
        /*-------------------------------------------------------------------------
        Non-Terminal pair types: 10,11,12,13,14   N=N,P,As,Sb; O=O,S,Se,Te; C=C,Si
        ========================================

        10:  N(+)(IV)-C(-)(III)     => N(V)=C(IV)  (N has 3 or 2 bonds)
        11:  N(+)(IV)=C(-)(III)     => N(V)#C(IV)  (N has 3 or 2 bonds)
        12:  N(+)(IV)-N(-)(II)      => N(V)=N(III) (allow terminal H on N(-))
        13: -O(+)-C(-)(III)         => -O=C-
        14: -O(+)=C(-)(III)         => -O#C-
        15:  O(+)(III)-N(-)(II)     => O(IV)=N(III) (allow terminal H on N(-))
        --------------------------------------------------------------------------*/

        if (!type || (10 <= type && type <= 15)) /* djb-rwth: addressing LLVM warning */
        {
            for (i = 0; i < num_atoms && 0 < num_All; i++)
            {
                if (1 == at[i].charge &&
                     0 < num_N_plus + num_O_plus && 0 < num_C_minus + num_N_minus &&
                     4 >= nNoMetalNumBonds( at, i ) && 4 == nNoMetalBondsValence( at, i ) &&
                     0 == num_of_H( at, i ) &&
                     ion_el_group( at[i].el_number ) == EL_NUMBER_N)
                {
                    /* found non-terminal N(+)(IV) */
                    if (( !type || 10 == type ) && 0 < num_N_plus && 0 < num_C_minus)
                    {
                        int num_neigh = 0, pos_neigh = -1;
                        for (i1 = 0; i1 < at[i].valence; i1++)
                        {
                            n = at[i].neighbor[i1];
                            if (-1 == at[n].charge && 3 >= at[n].valence && 3 == at[n].chem_bonds_valence + NUMH( at, n ) &&
                                 /*0 == at[n].num_H &&*/
                                 at[i].bond_type[i1] == BOND_TYPE_SINGLE &&
                                 ion_el_group( at[n].el_number ) == EL_NUMBER_C)
                            {
                                /* found N(+)(IV)-C(-)(III); prepare conversion to N(V)=C(IV) */
                                num_neigh++;
                                pos_neigh = i1;
                            }
                        }
                        i1 = pos_neigh;
                        if (1 == num_neigh &&
                             at[i].bond_type[i1] <= BOND_TYPE_TRIPLE &&
                             !has_other_ion_neigh( at, i, n = at[i].neighbor[i1] ) &&
                             !has_other_ion_neigh( at, n, i ))
                        {
                            /*n = at[i].neighbor[i1=pos_neigh];*/
                            i2 = (int) ( is_in_the_list( at[n].neighbor, (AT_NUMB) i, at[n].valence ) - at[n].neighbor );
                            at[i].charge--;
                            at[n].charge++;
                            at[i].chem_bonds_valence++;
                            at[n].chem_bonds_valence++;
                            at[i].bond_type[i1] ++;
                            at[n].bond_type[i2] ++;
                            num_changes++;
                            num_C_minus--;
                            num_N_plus--;
                            num_All -= 2;
                            continue;
                        }
                    }
                    if (( !type || 11 == type ) && 0 < num_N_plus && 0 < num_C_minus)
                    {
                        int num_neigh = 0, pos_neigh = -1;
                        for (i1 = 0; i1 < at[i].valence; i1++)
                        {
                            n = at[i].neighbor[i1];
                            if (-1 == at[n].charge && 3 >= at[n].valence && 3 == at[n].chem_bonds_valence + NUMH( at, n ) &&
                                 /*0 == at[n].num_H &&*/
                                 at[i].bond_type[i1] == BOND_TYPE_DOUBLE &&
                                 ion_el_group( at[n].el_number ) == EL_NUMBER_C)
                            {
                                /* found N(+)(IV)=C(-)(III); prepare conversion to N(V)#C(IV) */
                                num_neigh++;
                                pos_neigh = i1;
                            }
                        }
                        if (1 == num_neigh &&
                             !has_other_ion_neigh( at, i, n = at[i].neighbor[i1 = pos_neigh]) &&
                             !has_other_ion_neigh( at, n, i))
                        {
                            /*n = at[i].neighbor[i1=pos_neigh];*/
                            i2 = (int) ( is_in_the_list( at[n].neighbor, (AT_NUMB) i, at[n].valence ) - at[n].neighbor );
                            at[i].charge--;
                            at[n].charge++;
                            at[i].chem_bonds_valence++;
                            at[n].chem_bonds_valence++;
                            at[i].bond_type[i1] ++;
                            at[n].bond_type[i2] ++;
                            num_changes++;
                            num_C_minus--;
                            num_N_plus--;
                            num_All -= 2;
                            continue;
                        }
                    }
                    if (!type || (12 == type && 0 < num_N_plus && 0 < num_N_minus)) /* djb-rwth: addressing LLVM warning */
                    {
                        int num_neigh = 0, pos_neigh = -1;
                        for (i1 = 0; i1 < at[i].valence; i1++)
                        {
                            n = at[i].neighbor[i1];
                            if (-1 == at[n].charge && 2 >= nNoMetalNumBonds( at, n ) &&
                                 2 == nNoMetalBondsValence( at, n ) + NUMH( at, n ) &&
                                 /*0 == num_of_H( at, n ) &&*/
                                 at[i].bond_type[i1] == BOND_TYPE_SINGLE &&
                                 ion_el_group( at[n].el_number ) == EL_NUMBER_N)
                            {
                                /* found N(+)(IV)=N(-)(II); prepare conversion to N(V)#N(III) */
                                num_neigh++;
                                pos_neigh = i1;
                            }
                        }
                        if (1 == num_neigh &&
                             !has_other_ion_neigh( at, i, n = at[i].neighbor[i1 = pos_neigh]) &&
                             !has_other_ion_neigh( at, n, i))
                        {
                            /*n = at[i].neighbor[i1=pos_neigh];*/
                            i2 = (int) ( is_in_the_list( at[n].neighbor, (AT_NUMB) i, at[n].valence ) - at[n].neighbor );
                            at[i].charge--;
                            at[n].charge++;
                            at[i].chem_bonds_valence++;
                            at[n].chem_bonds_valence++;
                            at[i].bond_type[i1] ++;
                            at[n].bond_type[i2] ++;
                            num_changes++;
                            num_N_minus--;
                            num_N_plus--;
                            num_All -= 2;
                            continue;
                        }
                    }
                }
                else if (1 == at[i].charge &&
                          0 < num_O_plus && 0 < num_C_minus + num_N_minus &&
                          3 >= nNoMetalNumBonds( at, i ) && 3 == nNoMetalBondsValence( at, i ) &&
                          0 == num_of_H( at, i ) &&
                          ion_el_group( at[i].el_number ) == EL_NUMBER_O)
                {
                    /* found non-terminal O(+)(III) */
                    if (( !type || 13 == type ) && 0 < num_C_minus)
                    {
                        int num_neigh = 0, pos_neigh = -1;
                        for (i1 = 0; i1 < at[i].valence; i1++)
                        {
                            n = at[i].neighbor[i1];
                            if (-1 == at[n].charge && 3 >= at[n].valence && 3 == at[n].chem_bonds_valence + NUMH( at, n ) &&
                                 /*0 == at[n].num_H &&*/
                                 at[i].bond_type[i1] == BOND_TYPE_SINGLE &&
                                 ion_el_group( at[n].el_number ) == EL_NUMBER_C)
                            {
                                /* found O(+)(III)-C(-)(II); prepare conversion to O(IV)=C(IV) */
                                num_neigh++;
                                pos_neigh = i1;
                            }
                        }
                        if (1 == num_neigh &&
                             !has_other_ion_neigh( at, i, n = at[i].neighbor[i1 = pos_neigh]) &&
                             !has_other_ion_neigh( at, n, i))
                        {
                            /*n = at[i].neighbor[i1=pos_neigh];*/
                            i2 = (int) ( is_in_the_list( at[n].neighbor, (AT_NUMB) i, at[n].valence ) - at[n].neighbor );
                            at[i].charge--;
                            at[n].charge++;
                            at[i].chem_bonds_valence++;
                            at[n].chem_bonds_valence++;
                            at[i].bond_type[i1] ++;
                            at[n].bond_type[i2] ++;
                            num_changes++;
                            num_C_minus--;
                            num_O_plus--;
                            num_All -= 2;
                            continue;
                        }
                    }
                    if (( !type || 14 == type ) && 0 < num_C_minus)
                    {
                        int num_neigh = 0, pos_neigh = -1;
                        for (i1 = 0; i1 < at[i].valence; i1++)
                        {
                            n = at[i].neighbor[i1];
                            if (-1 == at[n].charge && 3 >= at[n].valence && 3 == at[n].chem_bonds_valence + NUMH( at, n ) &&
                                 /*0 == at[n].num_H &&*/
                                 at[i].bond_type[i1] == BOND_TYPE_DOUBLE &&
                                 ion_el_group( at[n].el_number ) == EL_NUMBER_C)
                            {
                                /* found O(+)(III)=C(-)(III); prepare conversion to O(IV)#C(IV) */
                                num_neigh++;
                                pos_neigh = i1;
                            }
                        }
                        if (1 == num_neigh &&
                             !has_other_ion_neigh( at, i, n = at[i].neighbor[i1 = pos_neigh]) &&
                             !has_other_ion_neigh( at, n, i))
                        {
                            /*n = at[i].neighbor[i1=pos_neigh];*/
                            i2 = (int) ( is_in_the_list( at[n].neighbor, (AT_NUMB) i, at[n].valence ) - at[n].neighbor );
                            at[i].charge--;
                            at[n].charge++;
                            at[i].chem_bonds_valence++;
                            at[n].chem_bonds_valence++;
                            at[i].bond_type[i1] ++;
                            at[n].bond_type[i2] ++;
                            num_changes++;
                            num_C_minus--;
                            num_O_plus--;
                            num_All -= 2;
                            continue;
                        }
                    }
                    if (( !type || 15 == type ) && 0 < num_N_minus)
                    {
                        int num_neigh = 0, pos_neigh = -1;
                        for (i1 = 0; i1 < at[i].valence; i1++)
                        {
                            n = at[i].neighbor[i1];
                            if (-1 == at[n].charge && 2 >= nNoMetalNumBonds( at, n ) &&
                                 2 == nNoMetalBondsValence( at, n ) + NUMH( at, n ) &&
                                 /*0 == num_of_H( at, n ) &&*/
                                 at[i].bond_type[i1] == BOND_TYPE_SINGLE &&
                                 ion_el_group( at[n].el_number ) == EL_NUMBER_N)
                            {
                                /* found O(+)(III)=N(-)(II); prepare conversion to O(IV)#N(III) */
                                num_neigh++;
                                pos_neigh = i1;
                            }
                        }
                        if (1 == num_neigh &&
                             !has_other_ion_neigh( at, i, n = at[i].neighbor[i1 = pos_neigh]) &&
                             !has_other_ion_neigh( at, n, i))
                        {
                            /*n = at[i].neighbor[i1=pos_neigh];*/
                            i2 = (int) ( is_in_the_list( at[n].neighbor, (AT_NUMB) i, at[n].valence ) - at[n].neighbor );
                            at[i].charge--;
                            at[n].charge++;
                            at[i].chem_bonds_valence++;
                            at[n].chem_bonds_valence++;
                            at[i].bond_type[i1] ++;
                            at[n].bond_type[i2] ++;
                            num_changes++;
                            num_N_minus--;
                            num_O_plus--;
                            num_All -= 2;
                            continue;
                        }
                    }
                }
            }
        }

        /**************************************************************************/
        /*********************** NON-Terminal ion triples *************************/
        /**************************************************************************/
        /*-------------------------------------------------------------------------
        Non-Terminal triple types: 16, 17, 18   N=N,P,As,Sb; O=O,S,Se,Te; C=C,Si
        ========================================
        16: C(+)(III)-O-N(-)(II)  => C(IV)=O=N(III)  (allow terminal H on N(-))
        |                     |
        17: C(+)(III)-N-C(-)(III)  => C(IV)=N=C(IV)

        18: C(-)(III)-N=C(+)(III)  => C(IV)=N#C(IV)   (may have two or no charges)
        C(IV)=N-C(II)          => C(IV)=N#C(IV)

        */

        if (( !type || 16 == type ) && 0 < num_C_plus && 0 < num_N_minus)
        {
            int m[2], j[2], k;
            for (i = 0; i < num_atoms; i++)
            {
                if (0 == at[i].charge && 2 == nNoMetalNumBonds( at, i ) && 2 == nNoMetalBondsValence( at, i ) &&
                     0 == num_of_H( at, i ) &&
                     0 <= ( j[0] = nNoMetalNeighIndex( at, i ) ) &&
                     at[m[0] = at[i].neighbor[j[0]]].charge &&
                     0 <= ( j[1] = nNoMetalOtherNeighIndex( at, i, m[0] ) ) &&
                     0 == at[m[0]].charge + at[m[1] = at[i].neighbor[j[1]]].charge &&
                     5 >= nNoMetalBondsValence( at, m[0] ) + nNoMetalBondsValence( at, m[1] ) &&
                     /*5 >= at[m[0]].chem_bonds_valence + at[m[1]].chem_bonds_valence &&*/
                     ion_el_group( at[i].el_number ) == EL_NUMBER_O)
                {
                    /* found non-terminal A(+)-O-B(-); chem_bond_val of A+B <= 5 */
                    int n_N = -1, n_C = -1, i_C = -1;
                    for (k = 0; k < 2; k++)
                    {
                        n = m[k];
                        if (-1 == at[n].charge && 2 == nNoMetalNumBonds( at, n ) + NUMH( at, n ) &&
                             /*0 == num_of_H( at, n ) &&*/
                             ion_el_group( at[n].el_number ) == EL_NUMBER_N)
                        {
                            n_N = n;
                        }
                        else if (1 == at[n].charge && 3 == at[n].chem_bonds_valence + NUMH( at, n ) &&
                                  ion_el_group( at[n].el_number ) == EL_NUMBER_C)
                        {
                            n_C = n;
                            i_C = k;
                        }
                    }
                    if (n_C < 0 || n_N < 0 ||
                         has_other_ion_in_sphere_2( at, n_C, n_N) ||
                         has_other_ion_in_sphere_2( at, n_N, n_C))
                    {
                        continue;
                    }

                    /* C(+)(III)-O-N(-)(II)  => C(IV)=O=N(III) */
                    for (k = 0; k < 2; k++)
                    {
                        n = k ? n_C : n_N;
                        i1 = k ? j[i_C] : j[1 - i_C];
                        i2 = (int) ( is_in_the_list( at[n].neighbor, (AT_NUMB) i, at[n].valence ) - at[n].neighbor );
                        at[i].bond_type[i1] ++;
                        at[n].bond_type[i2] ++;
                        at[i].chem_bonds_valence++;
                        at[n].chem_bonds_valence++;
                        at[n].charge += ( k ? -1 : 1 );
                    }
                    num_changes++;
                    num_N_minus--;
                    num_C_plus--;
                    num_All -= 2;
                }
            }
        }

        if (( !type || 17 == type ) && 0 < num_C_plus && 0 < num_C_minus)
        {
            int m[3], c[3], j[3], k;
            for (i = 0; i < num_atoms; i++)
            {
                if (0 == at[i].charge && 3 == nNoMetalNumBonds( at, i ) && 3 == nNoMetalBondsValence( at, i ) &&
                     0 == num_of_H( at, i ) &&
                     0 <= ( j[0] = nNoMetalNeighIndex( at, i ) ) &&
                     0 <= ( j[1] = nNoMetalOtherNeighIndex( at, i, m[0] = at[i].neighbor[j[0]] ) ) &&
                     0 <= ( j[2] = nNoMetalOtherNeighIndex2( at, i, m[0], m[1] = at[i].neighbor[j[1]] ) ) &&
                     1 == !( c[0] = at[m[0]].charge )
                     + !( c[1] = at[m[1]].charge )
                     + !( c[2] = at[m[2] = at[i].neighbor[j[2]]].charge ) &&
                     0 == c[0] + c[1] + c[2] &&
                     2 == ( 3 == ( c[0] ? at[m[0]].chem_bonds_valence + NUMH( at, m[0] ) : 0 ) )
                     + ( 3 == ( c[1] ? at[m[1]].chem_bonds_valence + NUMH( at, m[1] ) : 0 ) )
                     + ( 3 == ( c[2] ? at[m[2]].chem_bonds_valence + NUMH( at, m[2] ) : 0 ) ) &&
                     ion_el_group( at[i].el_number ) == EL_NUMBER_N)
                {
                    /* found non-terminal A(+)-O-B(-) */
                    int n_Cp = -1, n_Cm = -1, i_Cp = -1, i_Cm = -1; /* p = positive, m = negatice ion C */
                    for (k = 0; k < 3; k++)
                    {
                        if (c[k])
                        {
                            n = m[k];
                            if (-1 == at[n].charge &&
                                 ion_el_group( at[n].el_number ) == EL_NUMBER_C)
                            {
                                n_Cm = n;
                                i_Cm = k;
                            }
                            else if (1 == at[n].charge &&
                                      ion_el_group( at[n].el_number ) == EL_NUMBER_C)
                            {
                                n_Cp = n;
                                i_Cp = k;
                            }
                        }
                    }
                    if (n_Cp < 0 || n_Cm < 0 ||
                         has_other_ion_in_sphere_2( at, n_Cp, n_Cm) ||
                         has_other_ion_in_sphere_2( at, n_Cm, n_Cp))
                    {
                        continue;
                    }

                    /*           |                     |       */
                    /* C(+)(III)-N-C(-)(III)  => C(IV)=N=C(IV) */
                    for (k = 0; k < 2; k++)
                    {
                        n = k ? n_Cp : n_Cm;
                        i1 = k ? j[i_Cp] : j[i_Cm];
                        i2 = (int) ( is_in_the_list( at[n].neighbor, (AT_NUMB) i, at[n].valence ) - at[n].neighbor );
                        at[i].bond_type[i1] ++;
                        at[n].bond_type[i2] ++;
                        at[i].chem_bonds_valence++;
                        at[n].chem_bonds_valence++;
                        at[n].charge += ( k ? -1 : 1 );
                    }
                    num_changes++;
                    num_C_minus--;
                    num_C_plus--;
                    num_All -= 2;
                }
            }
        }

        if (( !type || 18 == type ) && ( (0 < num_C_plus && 0 < num_C_minus) || 0 < num_C_II )) /* djb-rwth: addressing LLVM warning */
        {
            int m[2], v[2], j[2], k;
            for (i = 0; i < num_atoms; i++)
            {
                if (0 == at[i].charge && 2 == nNoMetalNumBonds( at, i ) && 3 == nNoMetalBondsValence( at, i ) &&
                     0 == num_of_H( at, i ) &&
                     0 <= ( j[0] = nNoMetalNeighIndex( at, i ) ) &&
                     0 <= ( j[1] = nNoMetalOtherNeighIndex( at, i, m[0] = at[i].neighbor[j[0]] ) ) &&
                     0 == at[m[0]].charge
                     + at[m[1] = at[i].neighbor[j[1]]].charge &&
                     6 == ( v[0] = at[m[0]].chem_bonds_valence + NUMH( at, m[0] ) )
                     + ( v[1] = at[m[1]].chem_bonds_valence + NUMH( at, m[1] ) ) &&
                     2 >= abs( v[0] - v[1] ) &&
                     ion_el_group( at[i].el_number ) == EL_NUMBER_N &&
                     ion_el_group( at[m[0]].el_number ) == EL_NUMBER_C &&
                     ion_el_group( at[m[1]].el_number ) == EL_NUMBER_C)
                {
                    /*                    n_Cm      i n_Cp */
                    /* found non-terminal C(-)(III)-N=C(+)(III) or C(IV)=N-C(II): Cm-N-Cp */
                    /* convert to C(IV)=N#C(IV) */
                    int n_Cp = -1, n_Cm = -1, i_Cp = -1, i_Cm = -1; /* p = positive, m = negatice ion C */
                    for (k = 0; k < 2; k++)
                    {
                        n = m[k];
                        if (v[k] == 4 || (v[k] == 3 && at[i].bond_type[j[k]] == BOND_TYPE_SINGLE)) /* djb-rwth: addressing LLVM warning */
                        {
                            n_Cm = n;
                            i_Cm = k;
                        }
                        else if (v[k] == 2 || (v[k] == 3 && at[i].bond_type[j[k]] == BOND_TYPE_DOUBLE)) /* djb-rwth: addressing LLVM warning */
                        {
                            n_Cp = n;
                            i_Cp = k;
                        }
                    }
                    if (n_Cp < 0 || n_Cm < 0 || at[n_Cp].valence + NUMH( at, n_Cp ) != 2)
                    {
                        continue; /* guarantees at[n_Cp].valence <= 2 */
                    }
                    if (v[i_Cp] == 2 || !at[n_Cp].charge)
                    {
                        if (at[n_Cp].valence == 2)
                        {
                            /* neighbor of at[n_Cp] opposite to at[i] */
                            k = at[n_Cp].neighbor[at[n_Cp].neighbor[0] == i];
                            if (ion_el_group( at[k].el_number ) == EL_NUMBER_N)
                            {
                                continue;
                            }
                        }
                    }
                    else if (at[n_Cp].charge)
                    {
                        if (has_other_ion_in_sphere_2( at, n_Cp, n_Cm) ||
                             has_other_ion_in_sphere_2( at, n_Cm, n_Cp))
                        {
                            continue;
                        }
                    }
                    else
                    {
                        continue; /* unknown case */
                    }

                    /*                                         */
                    /* C(-)(III)-N=C(+)(III)  => C(IV)=N#C(IV) */
                    /* C(IV)=N-C(II)          => C(IV)=N#C(IV) */
                    if (at[n_Cp].charge)
                    {
                        num_C_minus--;
                        num_C_plus--;
                        num_All -= 2;
                    }
                    else
                    {
                        num_C_II--;
                        num_All--;
                    }

                    for (k = 0; k < 2; k++)
                    {
                        n = k ? n_Cp : n_Cm;
                        i3 = k ? i_Cp : i_Cm; /* added to fix the bug */
                                              /*i1 = k? j[i_Cp] : j[i_Cm];*/ /* replaced with next line */
                        i1 = j[i3];
                        if (v[i3 /*was i1*/] < 4)
                        {
                            /* WDI found a bug here: bounds violation */
                            int delta = 4 - v[i3 /*was i1*/];
                            i2 = (int) ( is_in_the_list( at[n].neighbor, (AT_NUMB) i, at[n].valence ) - at[n].neighbor );
                            at[i].bond_type[i1] += delta;
                            at[n].bond_type[i2] += delta;
                            at[i].chem_bonds_valence += delta;
                            at[n].chem_bonds_valence += delta;
                            at[n].charge = 0;
                            at[n].radical = 0;
                        }
                    }
                    at[i].charge = 0;
                    at[i].radical = 0;
                    num_changes++;
                }
            }
        }
    }

    return num_changes;
}




/*#if ( DISCONNECT_SALTS == 1 )*/ /* { */


                                  /****************************************************************************/
int RemoveInpAtBond( inp_ATOM *atom, int iat, int k )
{
    int      i, j, m, m2; /* djb-rwth: removing redundant variables */
    inp_ATOM *at = atom + iat;
    inp_ATOM *at2 = NULL;
    int      val = at->valence - 1;

    if (val >= 0)
    {
        int bond = at->bond_type[k];
        if (bond > BOND_TYPE_TRIPLE)
            bond = BOND_TYPE_SINGLE; /* added 08-06-2003 */

                                     /* update CML tetrahedral atom parity. */
        if (at->p_parity)
        {
            for (m = 0; m < MAX_NUM_STEREO_ATOM_NEIGH; m++)
            {
                if (at->p_orig_at_num[m] == at->orig_at_number)
                {
                    at->p_parity = 0;
                    break; /* only 3 bonds are present; removing one bond removes stereo */
                }
            }
            if (at->p_parity /* at->valence == MAX_NUM_STEREO_ATOM_NEIGH*/)
            {
                /* p_orig_at_num is a fixed size array of MAX_NUM_STEREO_ATOM_NEIGH (4) elements */
                for (m = 0; m < at->valence && m < MAX_NUM_STEREO_ATOM_NEIGH; m++) /* djb-rwth: fixing GH PR #72 */
                {
                    if (atom[(int) at->neighbor[k]].orig_at_number == at->p_orig_at_num[m])
                    {
                        break;
                    }
                }
                if (m < at->valence && m < MAX_NUM_STEREO_ATOM_NEIGH) /* djb-rwth: fixing GH PR #72 */
                {
                    at->p_orig_at_num[m] = at->orig_at_number;
                }
                else
                {
                    at->p_parity = 0; /* wrong neighbors: at->neighbor[k] is not in the list of a stereo neighbors */
                }
            }
        }


        /* update CML stereogenic bond parities; at this point no removed explicit H exist yet */
        if (at->sb_parity[0])
        {
            for (m = 0; m < MAX_NUM_STEREO_BONDS && at->sb_parity[m]; )
            {
                if (k == at->sb_ord[m] || (k == at->sn_ord[m] && val < 2 && ATOM_PARITY_WELL_DEF( at->sb_parity[m] ))) /* djb-rwth: addressing LLVM warning */
                {
                    /* !!! FLAW: does take into account removed H !!! */
                    /* stereogenic bond is being removed OR */
                    /* remove stereogenic bond because its only neighbor is being removed */
                    int pnxt_atom, pinxt2cur, pinxt_sb_parity_ord;
                    int len = get_opposite_sb_atom( atom, iat, at->sb_ord[m], &pnxt_atom, &pinxt2cur, &pinxt_sb_parity_ord );
                    if (len)
                    {
                        i = pinxt_sb_parity_ord;
                        at2 = atom + pnxt_atom;
                        /* djb-rwth: removing redundant code */
                    }
                    else
                    {
                        i = MAX_NUM_STEREO_BONDS;
                    }
                    /*
                    at2 = atom + at->neighbor[ (int)at->sb_ord[m] ];
                    for ( i = 0; i < MAX_NUM_STEREO_BONDS && at2->sb_parity[i]; i ++ )
                    {
                    if ( iat == at2->neighbor[ (int)at2->sb_ord[i] ] )
                    break;
                    }
                    */
                    if (i < MAX_NUM_STEREO_BONDS && at2->sb_parity[i])
                    {
                        m2 = i;
                        /* remove bond parity from at */
                        if (m < MAX_NUM_STEREO_BONDS - 1)
                        {
                            memmove(at->sb_parity + m, at->sb_parity + m + 1, (MAX_NUM_STEREO_BONDS - 1 - (long long)m) * sizeof(at->sb_parity[0])); /* djb-rwth: cast operator added */
                            memmove(at->sb_ord + m, at->sb_ord + m + 1, (MAX_NUM_STEREO_BONDS - 1 - (long long)m) * sizeof(at->sb_ord[0])); /* djb-rwth: cast operator added */
                            memmove(at->sn_ord + m, at->sn_ord + m + 1, (MAX_NUM_STEREO_BONDS - 1 - (long long)m) * sizeof(at->sn_ord[0])); /* djb-rwth: cast operator added */
                            memmove(at->sn_orig_at_num + m, at->sn_orig_at_num + m + 1, (MAX_NUM_STEREO_BONDS - 1 - (long long)m) * sizeof(at->sn_orig_at_num[0])); /* djb-rwth: cast operator added */
                        }
                        at->sb_parity[MAX_NUM_STEREO_BONDS - 1] = 0;
                        at->sb_ord[MAX_NUM_STEREO_BONDS - 1] = 0;
                        at->sn_ord[MAX_NUM_STEREO_BONDS - 1] = 0;
                        at->sn_orig_at_num[MAX_NUM_STEREO_BONDS - 1] = 0;
                        /* remove bond parity from at2 */
                        if (m2 < MAX_NUM_STEREO_BONDS - 1)
                        {
                            memmove(at2->sb_parity + m2, at2->sb_parity + m2 + 1, (MAX_NUM_STEREO_BONDS - 1 - (long long)m2) * sizeof(at2->sb_parity[0])); /* djb-rwth: cast operator added */
                            memmove(at2->sb_ord + m2, at2->sb_ord + m2 + 1, (MAX_NUM_STEREO_BONDS - 1 - (long long)m2) * sizeof(at2->sb_ord[0])); /* djb-rwth: cast operator added */
                            memmove(at2->sn_ord + m2, at2->sn_ord + m2 + 1, (MAX_NUM_STEREO_BONDS - 1 - (long long)m2) * sizeof(at2->sn_ord[0])); /* djb-rwth: cast operator added */
                            memmove(at2->sn_orig_at_num + m2, at2->sn_orig_at_num + m2 + 1, (MAX_NUM_STEREO_BONDS - 1 - (long long)m2) * sizeof(at2->sn_orig_at_num[0])); /* djb-rwth: cast operator added */
                        }
                        at2->sb_parity[MAX_NUM_STEREO_BONDS - 1] = 0;
                        at2->sb_ord[MAX_NUM_STEREO_BONDS - 1] = 0;
                        at2->sn_ord[MAX_NUM_STEREO_BONDS - 1] = 0;
                        at2->sn_orig_at_num[MAX_NUM_STEREO_BONDS - 1] = 0;
                        /* do not increment m here because the array elements have been shifted */
                    }
                    else
                    {
                        m++; /* program error: inconsistent stereobond parity */
                    }
                }
                else if (k == at->sn_ord[m])
                {
                    /* stereogenic bond neighbor is being removed; another neighbor remains */
                    /* !!! FLAW: does take into account removed H !!! */
                    for (j = 0, i = -1; j < at->valence; j++)
                    {
                        if (j != k && j != at->sb_ord[m])
                        {
                            i = j;
                            break;
                        }
                    }

                    /* i is the position of the neighbor that will become a new neighbor */
                    /***************************************************************************
                    *  at->sb_parity[m] is the direction (EVEN=clockwise, ODD=counterclockwise)
                    *  from stereobond to the neighbor. If the neighbor is removed then
                    *  the parity should invert, otherwise it should be unchanged.
                    ***************************************************************************/
                    if (i < 0)
                    {
                        /* no alternative neighbor is available */
                        if (ATOM_PARITY_WELL_DEF( at->sb_parity[m] ))
                        {
                            /* parity cannot be not well-defined anymore */
                            int pnxt_atom, pinxt2cur, pinxt_sb_parity_ord;
                            int len = get_opposite_sb_atom( atom, iat, at->sb_ord[m], &pnxt_atom, &pinxt2cur, &pinxt_sb_parity_ord );
                            if (len > 0)
                            {
                                atom[pnxt_atom].sb_parity[pinxt_sb_parity_ord] = at->sb_parity[m] = AB_PARITY_UNDF;
                            }
                        }
                        at->sn_ord[m] = -99; /* sb neighbor has been disconnected */
                        at->sb_ord[m] -= ( at->sb_ord[m] > k ); /* same as above */
                        at->sn_orig_at_num[m] = 0;
                    }
                    else if (i < at->valence)
                    {
                        /* choose another stereogenic bond neighbor, its ord. number is i before bond removal */
                        if (ATOM_PARITY_WELL_DEF( at->sb_parity[m] ))
                        {
                            /* ALL WRONG: 'move' previous stereo bond neighbor to the last position (pos. 2 out of 0,1,2) */
                            /* the parity of the transpositions is (2 - at->sn_ord[m])%2 = at->sn_ord[m] % 2 */
                            /* and replace the neighbor with another; the contribution to the parity is 1 */

                            /*at->sb_parity[m]      =  2 - ( at->sb_parity[m] + at->sn_ord[m] + 1 ) % 2;*/

                            /*at->sb_parity[m]      =  2 - ( at->sb_parity[m] + k + i +
                            (i > k) + (i > at->sb_ord[m]) ) % 2;*/
                            /*=== parity should be INVERTED ===*/
                            at->sb_parity[m] = 3 - at->sb_parity[m];
                        }
                        at->sn_ord[m] = i - ( i > k ); /* ord. number shifted because preceding bond is removed */
                        at->sb_ord[m] -= ( at->sb_ord[m] > k ); /* same as above */
                        at->sn_orig_at_num[m] = atom[(int) at->neighbor[i]].orig_at_number;
                        /*at->sb_parity[m]      =  2 - ( at->sb_parity[m] + 1 ) % 2;*/
                    }
                    else
                    {
                        at->sb_parity[m] = 0; /* program error: inconsistent stereobond parity */
                    }
                    m++;
                }
                else
                {
                    /* removing another neighbor, k: first move it to the last position (pos. 2 out of 0,1,2) */
                    if (k < 2 && ATOM_PARITY_WELL_DEF( at->sb_parity[m] ))
                    {
                        /*at->sb_parity[m] =  2 - ( at->sb_parity[m] + k ) % 2;*/
                        /*at->sb_parity[m] =  2 - ( at->sb_parity[m] + (at->sn_ord[m] > k) + (at->sb_ord[m] > k) ) % 2;*/
                        ;/*==== Parity should remain UNCHANGED ===*/
                    }
                    if (at->sb_ord[m] > k)
                    {
                        at->sb_ord[m] --;
                    }
                    if (at->sn_ord[m] > k)
                    {
                        at->sn_ord[m] --;
                    }
                    m++;
                }
            }
        }

        if (k < val)
        {
            memmove(at->neighbor + k, at->neighbor + k + 1, sizeof(at->neighbor[0])* ((long long)val - (long long)k)); /* djb-rwth: cast operators added */
            memmove(at->bond_stereo + k, at->bond_stereo + k + 1, sizeof(at->bond_stereo[0])* ((long long)val - (long long)k)); /* djb-rwth: cast operators added */
            memmove(at->bond_type + k, at->bond_type + k + 1, sizeof(at->bond_type[0])* ((long long)val - (long long)k)); /* djb-rwth: cast operators added */
        }

        at->neighbor[val] = 0;
        at->bond_stereo[val] = 0;
        at->bond_type[val] = 0;
        at->valence = val;
        at->chem_bonds_valence -= bond;
        return 1;
    }

    return 0;
}


/****************************************************************************/
int DisconnectInpAtBond( inp_ATOM *at,
                         AT_NUMB *nOldCompNumber,
                         int iat,
                         int neigh_ord )
{
    int neigh, i, ret = 0;
    int component;
    neigh = at[iat].neighbor[neigh_ord];

    for (i = 0; i < at[neigh].valence; i++)
    {
        if (iat == (int) at[neigh].neighbor[i])
        {
            break;
        }
    }

    if (i < at[neigh].valence)
    {
        ret += RemoveInpAtBond( at, iat, neigh_ord );
        ret += RemoveInpAtBond( at, neigh, i );
        if (nOldCompNumber && ret)
        {
            if ((component = at[iat].component)) /* djb-rwth: addressing LLVM warning */
            {
                nOldCompNumber[component - 1] = 0;
            }
            if ((component = at[neigh].component)) /* djb-rwth: addressing LLVM warning */
            {
                nOldCompNumber[component - 1] = 0;
            }
        }
    }

    return ( ret == 2 );
}


/****************************************************************************/
int bIsAmmoniumSalt( inp_ATOM *at,
                     int i,
                     int *piO,
                     int *pk,
                     S_CHAR *num_explicit_H )
{
    /* NH4(+charge)-O(-charge)-C -> NH3 + HO-C; any charge including 0, */
    /* any C except charged or radical F, Cl, Br, I                     */

    int num_H, num_non_iso_H, num_impl_iso_H, bDisconnect = 1;
    int j, val, neigh, iO = -1, iC, k = -1;

    if (at[i].el_number != EL_NUMBER_N)
    {
        return 0;
    }

    /* check for NH4-O-C... -> NH3 + HO-C... */
    val = at[i].valence;
    num_impl_iso_H = NUM_ISO_H( at, i );
    num_non_iso_H = at[i].num_H;
    num_H = num_non_iso_H + num_impl_iso_H;
    if (val + num_H == 5)
    {
        int num_O = 0;
        memset( num_explicit_H, 0, ( NUM_H_ISOTOPES + 1 ) * sizeof( num_explicit_H[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
        for (j = 0; j < val; j++)
        { /* looking for O: H4N-O-C... */
            neigh = at[i].neighbor[j];
            if (at[neigh].num_H ||
                 (at[neigh].charge && ( at[neigh].el_number != EL_NUMBER_O || at[neigh].charge + at[i].charge )) ||
                 (at[neigh].radical && at[neigh].radical != RADICAL_SINGLET)) /* djb-rwth: addressing LLVM warnings */
            {
                bDisconnect = 0;
                break; /* reject */
            }
            if (at[neigh].el_number == EL_NUMBER_H && at[neigh].valence == 1 &&
                 !at[neigh].charge && !at[neigh].radical)
            {
                num_H++; /* at this point at[].num_H does not include explicit H count */
                num_non_iso_H += ( 0 == at[neigh].iso_atw_diff );
                num_explicit_H[at[neigh].iso_atw_diff] ++;  /* explicit H on N */
            }
            else if (at[neigh].el_number == EL_NUMBER_O && at[neigh].valence == 2 && !num_O)
            {
                num_O++; /* found O: N-O- */
                iO = neigh;
                k = j;
                iC = at[iO].neighbor[at[iO].neighbor[0] == i];
                if (at[iC].el_number != EL_NUMBER_C || /*
                                                       at[iC].num_H ||
                                                       at[iC].chem_bonds_valence != 4 || */
                     at[iC].charge ||
                     (at[iC].radical && at[iC].radical != RADICAL_SINGLET) /*||
                                                                         at[iC].valence == at[iC].chem_bonds_valence*/) /* djb-rwth: addressing LLVM warning */
                {
                    bDisconnect = 0;
                    break; /* reject */
                }
            }
            else if (( at[neigh].el_number == EL_NUMBER_F ||
                       at[neigh].el_number == EL_NUMBER_CL ||
                       at[neigh].el_number == EL_NUMBER_BR ||
                       at[neigh].el_number == EL_NUMBER_I ) &&
                      at[neigh].valence == 1 && at[neigh].chem_bonds_valence == 1 &&
                      !at[neigh].charge && !NUMH( at, neigh ) && !num_O)
            {
                num_O++; /* found O: N-O- */
                iO = neigh;
                k = j;
                /* djb-rwth: removing redundant code */
            }
            else
            {
                bDisconnect = 0;
                break;  /* reject */
            }
        }
        if (bDisconnect && ( num_O != 1 || num_H != 4 ))
        {
            bDisconnect = 0; /* reject */
        }
    }
    else
    {
        bDisconnect = 0;
    }
    if (bDisconnect)
    {
        *piO = iO;
        *pk = k;
    }

    return bDisconnect;
}


/****************************************************************************/
int DisconnectAmmoniumSalt( inp_ATOM *at,
                            int iN,
                            int iO,
                            int k,
                            S_CHAR *num_explicit_H )
{

    /* disconnect NH4-O from O */
    /* Note: iO = at[iN].neighbor[k], at[iN] is N, at[iO].neighbor[0] is either N=at[iN] or C=at[iC] */

    int nMove_H_iso_diff = -1; /* do not move explicit H */
    int j, neigh, iso_diff, neigh_pos;
    int    val = at[iN].valence;

    if (at[iN].charge && !( at[iN].charge + at[iO].charge ))
    {
        at[iN].charge = at[iO].charge = 0; /* remove charges */
    }

    neigh_pos = ( at[iO].valence == 2 ) ? ( at[iO].neighbor[1] == iN ) : 0; /* position of at[iN] in the neigh list of iO */
                                                                            /* disconnect bond O-N */
    RemoveInpAtBond( at, iO, neigh_pos );
    RemoveInpAtBond( at, iN, k );
    val--;

    /* move 1 H from NH4 to O- or Cl */

    /* find non-isotopic or the lightest isotopic H to move from N to O */
    for (iso_diff = 0; iso_diff <= NUM_H_ISOTOPES; iso_diff++) /* djb-rwth: fixing GH PR #72 */
    {
        if (!iso_diff)
        {
            /* find non-isotopic H */
            if (at[iN].num_H)
            {
                at[iN].num_H--;  /* move non-isotopic implicit H */
                at[iO].num_H++;
                break;
            }
            else if (num_explicit_H[0])
            {
                nMove_H_iso_diff = 0; /* flag: move explicit non-isotopic H */
                break;
            }
        }
        else
        {
            /* find isotopic H */
            /* num_iso_H has length NUM_H_ISOTOPES; do not access out-of-bounds */
            if ((iso_diff < NUM_H_ISOTOPES) && at[iN].num_iso_H[iso_diff]) /* djb-rwth: fixing GH PR #72 */
            {
                at[iN].num_iso_H[iso_diff] --; /* move implicit isotopic H, atw = 1 */
                at[iO].num_iso_H[iso_diff] ++;
                break;
            }
            else
            {
                if (num_explicit_H[iso_diff])
                {
                    nMove_H_iso_diff = iso_diff; /* flag: move explicit isotopic H, atw = 1 */
                    break;
                }
            }
        }
    }

    if (nMove_H_iso_diff >= 0)
    {
        /* move explicit H, it is isotopic if nMove_H_iso_diff > 0 */
        double dist2_H_O, min_dist2_H_O = -1.0;
        int    jH = -1, iH = -1;
        for (j = 0; j < val; j++)
        { /* looking H in N-H such that H-O is shortest */
            neigh = at[iN].neighbor[j];
            if (at[neigh].el_number == EL_NUMBER_H &&
                 at[neigh].iso_atw_diff == nMove_H_iso_diff)
            {
                dist2_H_O = ( at[neigh].x - at[iO].x ) * ( at[neigh].x - at[iO].x ) +
                    ( at[neigh].y - at[iO].y ) * ( at[neigh].y - at[iO].y ) +
                    ( at[neigh].z - at[iO].z ) * ( at[neigh].z - at[iO].z );
                if (min_dist2_H_O < 0.0 || min_dist2_H_O > dist2_H_O)
                {
                    min_dist2_H_O = dist2_H_O;
                    iH = neigh;
                    jH = j;
                }
            }
        }

        /* reconnect; bonds do not need changes except stereo */
        neigh_pos = at[iO].valence;
        at[iO].neighbor[neigh_pos] = iH;
        at[iO].bond_stereo[neigh_pos] = 0;
        at[iO].bond_type[neigh_pos] = at[iH].bond_type[0];
        at[iO].chem_bonds_valence += at[iH].bond_type[0];
        at[iO].valence++;
        at[iH].neighbor[0] = iO;
        at[iH].bond_stereo[0] = 0;

        /* disconnect H from N */
        RemoveInpAtBond( at, iN, jH );
        val--;
        if (k > jH)
        {
            k--;
        }
    }

    return 1;
}


/****************************************************************************/
int bIsMetalSalt( inp_ATOM *at, int i )
{
    int type, val, k, iO, iC, j, neigh;
    int bDisconnect = 1;

    /* check for a metal atom:
    metal atom should be connected and be a metal */
    if (!( val = at[i].valence ) ||
         !( type = get_el_type( at[i].el_number ) ) ||
         !( type & IS_METAL ))
    {
        bDisconnect = 0;  /* reject */
    }
    else if (at[i].num_H)
        /* metal atom should not have adjacent H or multiple bonds or radical */
    {
        bDisconnect = 0; /* reject */
    }
    else
    {
        /* check valence */
        if ((at[i].charge == 0 &&
            ( (( type & 1 ) && val == get_el_valence( at[i].el_number, 0, 0 )) ||
             (( type & 2 ) && val == get_el_valence( at[i].el_number, 0, 1 ) ))) ||
             (at[i].charge > 0 &&
             ( type & 1 ) && val == get_el_valence( at[i].el_number, at[i].charge, 0 ))) /* djb-rwth: addressing LLVM warnings */
        {
            ; /* accept */
        }
        else
        {
            bDisconnect = 0; /* reject */
        }
    }

    if (bDisconnect)
    {
        /*************************************************************************
        *                                                                  |    *
        * check M neighbors. Disconnect if all neighbors are M-O-C# or M-O-C=   *
        *                                                                  |    *
        *************************************************************************/
        for (k = 0; k < at[i].valence; k++)
        {
            iO = at[i].neighbor[k];
            /* halogenide 2004-07-08 */
            if (( at[iO].el_number == EL_NUMBER_F ||
                  at[iO].el_number == EL_NUMBER_CL ||
                  at[iO].el_number == EL_NUMBER_BR ||
                  at[iO].el_number == EL_NUMBER_I ) &&
                 at[iO].valence == 1 && at[iO].chem_bonds_valence == 1 &&
                 !at[iO].charge && !( at[iO].radical && at[iO].radical != RADICAL_SINGLET ) && !NUMH( at, iO ))
            {
                ; /* found */
            }
            else
            {
                /* -O-C= */
                if (at[iO].el_number != EL_NUMBER_O ||
                     NUMH( at, iO ) ||
                     at[iO].valence != 2 ||
                     at[iO].charge ||
                     (at[iO].radical && at[iO].radical != RADICAL_SINGLET) ||
                     at[iO].valence != at[iO].chem_bonds_valence) /* djb-rwth: addressing LLVM warning */
                {
                    bDisconnect = 0; /* reject */
                    break;
                }
                iC = at[iO].neighbor[at[iO].neighbor[0] == i];
                if (at[iC].el_number != EL_NUMBER_C ||
                     at[iC].num_H ||
                     at[iC].chem_bonds_valence != 4 ||
                     at[iC].charge ||
                     (at[iC].radical && at[iC].radical != RADICAL_SINGLET) ||
                     at[iC].valence == at[iC].chem_bonds_valence) /* djb-rwth: addressing LLVM warning */
                {
                    bDisconnect = 0; /* reject */
                    break;
                }
                for (j = 0; j < at[iC].valence; j++)
                {
                    neigh = at[iC].neighbor[j];
                    if (at[neigh].el_number == EL_NUMBER_H)
                    {
                        break;
                    }
                }
                if (j != at[iC].valence)
                {
                    bDisconnect = 0; /* reject */
                    break;
                }
            }
        }
    }

    return bDisconnect;
}


/****************************************************************************/
int DisconnectMetalSalt( inp_ATOM *at, int i )
{
    int k, iO;
    /* disconnect metal atom or ion at[i] */

    for (k = 0; k < at[i].valence; k++)
    {
        iO = at[i].neighbor[k];
        if (at[iO].valence == 2)
        {
            if (at[iO].neighbor[0] == i)
            {
                /* assuming atom O always has 2 bonds */
                /* copy the remaining neighbor to the 0 position */
                at[iO].neighbor[0] = at[iO].neighbor[1];
                at[iO].bond_stereo[0] = at[iO].bond_stereo[1];
                at[iO].bond_type[0] = at[iO].bond_type[1];
            }
            /* clear neighbor at position 1 */
            at[iO].neighbor[1] = 0;
            at[iO].bond_stereo[1] = 0;
            at[iO].bond_type[1] = 0;
        }
        else
        {
            /* clear neighbor at position 1 */
            at[iO].neighbor[0] = 0;
            at[iO].bond_stereo[0] = 0;
            at[iO].bond_type[0] = 0;
        }

        /* make O negatively charged */
        at[iO].charge = -1;

        /* reduce O valence to account for the removed single bond */
        at[iO].valence--;
        at[iO].chem_bonds_valence--;

        /* clear metal neighbor (O) */
        at[i].neighbor[k] = 0;
        at[i].bond_stereo[k] = 0;
        at[i].bond_type[k] = 0;

        /* add a positive charge to the metal */
        at[i].charge++;
    }

    /* set metal valence to zero because it has been disconnected */
    at[i].valence = 0;
    at[i].chem_bonds_valence = 0;

    return k;
}


/****************************************************************************/
int DisconnectSalts( ORIG_ATOM_DATA *orig_inp_data, int bDisconnect )
{
    int i, k, iO, num_changes, val;
    S_CHAR    num_explicit_H[NUM_H_ISOTOPES + 1];
    inp_ATOM *at = orig_inp_data->at;
    int num_at = orig_inp_data->num_inp_atoms;

    /* check each atom */
    for (i = 0, num_changes = 0; i < num_at; i++)
    {

        if (!( val = at[i].valence ) || /* disconnected atom */
             val != at[i].chem_bonds_valence || /* a bond has higher multiplicity than 1 */
             (at[i].radical && at[i].radical != RADICAL_SINGLET) /* radical */) /* djb-rwth: addressing LLVM warning */
        {
            continue;   /* reject */
        }

        if (bIsAmmoniumSalt( at, i, &iO, &k, num_explicit_H ))
        {
            if (bDisconnect)
            {
                DisconnectAmmoniumSalt( at, i, iO, k, num_explicit_H );
                orig_inp_data->num_inp_bonds--;
            }

            /* count disconnected atoms */
            num_changes++;
        }
        else if (bIsMetalSalt( at, i ))
        {
            if (bDisconnect)
            {
                k = DisconnectMetalSalt( at, i );
                orig_inp_data->num_inp_bonds -= k;
            }
            num_changes++;
        }
    }

    return num_changes;
}


/*****************************************************************************/
/* Important: Salt disconnection is independent from coord. disconnection:   */
/* because different atoms are disconnected.                                 */
/* However, sal disconnection may need to be rerun after metal disconnection */
/* because metal disconnection may make certain atoms be eligible for salt   */
/* disconnection                                                             */
/*****************************************************************************/
int bIsMetalToDisconnect( inp_ATOM *at, int i, int bCheckMetalValence )
{
    int type, at_valence, num_H;

    /*
    if ( !at[i].valence )
    */

    if (!( type = get_el_type( at[i].el_number ) ) ||
         !( type & IS_METAL ))
    {
        return 0;
    }

    num_H = NUMH( at, i );
    at_valence = num_H + at[i].chem_bonds_valence;

    if (!at_valence)
    {
        return 0; /* nothing to disconnect */
    }

    if (bCheckMetalValence)
    {
        if (abs( at[i].charge ) > 1)
        {
            return 1; /* multiple charges */
        }

        for (i = 0; i < 2 && ( i & type ); i++)
        {
            if (at_valence == get_el_valence( at[i].el_number, at[i].charge, i )) /* djb-rwth: fixing coverity ID #499532 -- unresolved issue -- revision required */
            {
                return 2; /* atom has normal valence */
            }
        }
    }

    return 1;
}


/****************************************************************************/
int bMayDisconnectMetals( ORIG_ATOM_DATA *orig_inp_data,
                          int bCheckMetalValence,
                          INCHI_MODE *bTautFlagsDone )
{
    int i, j, k, iO, num_changes, val, bRadOrMultBonds, num_impl_H = 0;
    S_CHAR    num_explicit_H[NUM_H_ISOTOPES + 1];
    inp_ATOM *at = orig_inp_data->at;
    int num_at = orig_inp_data->num_inp_atoms;
    int *nNumImplH = &orig_inp_data->bDisconnectCoord;

    /* check each atom */
    for (i = 0, num_changes = 0; i < num_at; i++)
    {

        if (!( val = at[i].valence ) && !NUMH( at, i ))
        {
            continue; /* disconnected atom */
        }

        bRadOrMultBonds = ( val == 0 ) ||
            ( val != at[i].chem_bonds_valence ) || /* a bond has higher multiplicity than 1 */
            ( at[i].radical && at[i].radical != RADICAL_SINGLET ); /* radical */

        if (!bRadOrMultBonds && bIsAmmoniumSalt( at, i, &iO, &k, num_explicit_H ))
        {
            ;
        }
        else if (!bRadOrMultBonds && bIsMetalSalt( at, i ))
        {
            ;
        }
        else if (1 == ( j = bIsMetalToDisconnect( at, i, bCheckMetalValence ) ))
        {
            num_impl_H += NUMH( at, i );
            num_changes++;
        }
        else if (2 == j && bTautFlagsDone)
        {
            *bTautFlagsDone |= TG_FLAG_CHECK_VALENCE_COORD_DONE;
        }
    }

    if (nNumImplH)
    {
        *nNumImplH = num_changes ? num_impl_H + 1 : 0;
    }

    return num_changes;
}


/****************************************************************************/


#if ( bRELEASE_VERSION == 0 && (EXTR_HAS_METAL_ATOM & (EXTR_MASK | EXTR_FLAG) ) )


/****************************************************************************/
int bHasMetalAtom( ORIG_ATOM_DATA *orig_inp_data )
{
    int i;
    inp_ATOM *at;

    if (orig_inp_data && ( at = orig_inp_data->at ))
    {
        int num_at = orig_inp_data->num_inp_atoms;
        /* check each atom */
        for (i = 0; i < num_at; i++)
        {
            if (IS_METAL & get_el_type( at[i].el_number ))
            {
                return 1;
            }
        }
    }

    return 0;
}
#endif



/*****************************************************************************
{ "F",   19,  19,  18.998403220,     0 ,  0, {{0,},       {0,},       {1,},       {2,},       {3,5},      },},
{ "Cl",  35,  35,  34.968852730,     0 ,  0, {{0,},       {0,},       {1,3,5,7},  {2,4,6},    {3,5,},     },},
{ "Br",  80,  79,  78.918336100,     0 ,  0, {{0,},       {0,},       {1,3,5,7,}, {2,4,6,},   {3,5,},     },},
{ "I",  127, 127, 126.904500000,     0 ,  0, {{0,},       {0,},       {1,3,5,7,}, {2,4,6},    {3,5,},     },},
{ "At", 210, 210, 209.987100000,     0 ,  0, {{0,},       {0,},       {1,3,5,7,}, {2,4,6},    {3,5,},     },},
{ "N",   14,  14,  14.003074000,     0 ,  0, {{1,},       {2,},       {3,5},      {4,},       {3,},       },},
{ "P",   31,  31,  30.973762000,     0 ,  0, {{1,3,5,7,}, {2,4,6,},   {3,5,},     {4,},       {3,},       },},
{ "As",  75,  75,  74.921594200,     0 ,  0, {{0,},       {2,4,6,},   {3,5,},     {4,},       {3,},       },},
{ "Sb", 122, 121, 120.903800000,     0 ,  0, {{1,3,5,7,}, {2,4,6,},   {3,5,},     {2,4,},     {3,},       },},
{ "O",   16,  16,  15.994914630,     0 ,  0, {{0,},       {1,},       {2,},       {3,5,},     {4,},       },},
{ "S",   32,  32,  31.972070700,     0 ,  0, {{0,},       {1,3,5,7,}, {2,4,6},    {3,5,},     {4,},       },},
{ "Se",  79,  80,  79.916519600,     0 ,  0, {{0,},       {1,3,5,7,}, {2,4,6,},   {3,5,},     {4,},       },},
{ "Te", 128, 130, 129.906200000,     0 ,  0, {{0,},       {1,3,5,7,}, {2,4,6,},   {3,5,},     {2,4,},     },},
{ "Po", 209, 209, 208.982400000,     0 ,  0, {{0,},       {1,3,5,7,}, {2,4,6,},   {3,5,},     {2,4,},     },},
{ "B",   11,  11,  11.009300000,     0 ,  0, {{3,},       {4,},       {3,},       {2,},       {1,},       },},
*****************************************************************************/



/****************************************************************************/
int DisconnectMetals( ORIG_ATOM_DATA *orig_inp_data,
                      int bCheckMetalValence,
                      INCHI_MODE *bTautFlagsDone )
{
    int i, j, k, n, iO, num_changes, val, bRadOrMultBonds;
    int num_impl_H, num_at, err, num_disconnected;
    S_CHAR num_explicit_H[NUM_H_ISOTOPES + 1];
    static char elnumber_Heteroat[16] = { '\0', };
    static int  num_halogens = 0;
    int num_halogens2;

    inp_ATOM  *at = NULL;
    S_CHAR    *bMetal = NULL;
    inp_ATOM  *atom = orig_inp_data->at;
    int        num_atoms = orig_inp_data->num_inp_atoms;
    int        nNumExplH = ( orig_inp_data->bDisconnectCoord > 0 ) ? orig_inp_data->bDisconnectCoord - 1 : 0;
    AT_NUMB   *nOldCompNumber = orig_inp_data->nOldCompNumber;

    err = 0;
    num_impl_H = 0;
    num_at = num_atoms;
    num_disconnected = 0;

    if (!( at = (inp_ATOM *) inchi_calloc( (long long)num_at + (long long)nNumExplH, sizeof( at[0] ) ) ) || /* djb-rwth: cast operators added */
         !( bMetal = (S_CHAR    *) inchi_calloc( (long long)num_at + (long long)nNumExplH, sizeof( bMetal[0] ) ) )) /* djb-rwth: cast operators added */
    {
        err = 1;
        goto exit_function;
    }

    if (!num_halogens) /* if (!elnumber_Heteroat[0] )  */
    {
        i = 0;
        /* halogens */
        elnumber_Heteroat[i++] = (char) EL_NUMBER_F; /* 0 */
        elnumber_Heteroat[i++] = (char) EL_NUMBER_CL;
        elnumber_Heteroat[i++] = (char) EL_NUMBER_BR;
        elnumber_Heteroat[i++] = (char) EL_NUMBER_I;
        elnumber_Heteroat[i++] = (char) EL_NUMBER_AT; /* 4 */
        num_halogens2 = i;
        /* other non-metal */
        elnumber_Heteroat[i++] = (char) EL_NUMBER_N;
        elnumber_Heteroat[i++] = (char) EL_NUMBER_P;
        elnumber_Heteroat[i++] = (char) EL_NUMBER_AS;
        /*elnumber_Heteroat[i++] = EL_NUMBER_SB;*/ /* metal 10-28-2003 */
        elnumber_Heteroat[i++] = (char) EL_NUMBER_O;
        elnumber_Heteroat[i++] = (char) EL_NUMBER_S;
        elnumber_Heteroat[i++] = (char) EL_NUMBER_SE;
        elnumber_Heteroat[i++] = (char) EL_NUMBER_TE;
        /*elnumber_Heteroat[i++] = EL_NUMBER_PO;*/ /* metal 10-28-2003 */
        elnumber_Heteroat[i++] = (char) EL_NUMBER_B;
        elnumber_Heteroat[i++] = 0;
        num_halogens = num_halogens2;
    }

    memcpy(at, atom, num_atoms * sizeof(at[0]));

    /* check each atom, mark metals */
    for (i = 0, k = 0, num_changes = 0; i < num_atoms; i++)
    {
        if (!( val = at[i].valence ) && !NUMH( at, i ))
        {
            continue; /* disconnected atom */
        }
        bRadOrMultBonds = ( val == 0 ) ||
            ( val != at[i].chem_bonds_valence ) || /* a bond has higher multiplicity than 1 */
            ( at[i].radical && at[i].radical != RADICAL_SINGLET ); /* radical */

        if (!bRadOrMultBonds && bIsAmmoniumSalt( at, i, &iO, &k, num_explicit_H ))
        {
            ;
        }
        else if (!bRadOrMultBonds && bIsMetalSalt( at, i ))
        {
            ;
        }
        else if (1 == ( j = bIsMetalToDisconnect( at, i, bCheckMetalValence ) ))
        {
            num_impl_H += ( k = NUMH( at, i ) );
            bMetal[i] = 1 + k;
            num_changes++;
        }
        else if (2 == j && bTautFlagsDone)
        {
            *bTautFlagsDone |= TG_FLAG_CHECK_VALENCE_COORD_DONE;
        }
    }

    if (num_impl_H != nNumExplH)
    {
        err = 2;
        goto exit_function;
    }

    /* replace implicit H atoms with explicit H atoms */
    for (i = 0; i < num_atoms && 0 < num_impl_H; i++)
    {
        if (bMetal[i] <= 1)
        {
            continue;
        }
        for (k = 0; k < NUM_H_ISOTOPES + 1; k++)
        {
            n = k ? at[i].num_iso_H[k - 1] : at[i].num_H;
            for (j = 0; j < n; j++)
            {
                if (num_at >= num_atoms + nNumExplH)
                {
                    err = 3;
                    goto exit_function;
                }
                at[num_at].elname[0] = 'H';
                at[num_at].el_number = get_periodic_table_number( at[num_at].elname );
                at[num_at].iso_atw_diff = k;
                at[num_at].component = at[i].component;
                move_explicit_Hcation( at, num_at + 1, i, num_at, 1 );
                at[num_at].orig_at_number = num_at + 1;
                num_at++;
                num_impl_H--;
                bMetal[i] --;
                if (k)
                {
                    at[i].num_iso_H[k - 1] --;
                }
                else
                {
                    at[i].num_H--;
                }
            }
        }

        if (bMetal[i] != 1)
        {
            err = 4;
            goto exit_function;
        }
    }

    if (num_at != num_atoms + nNumExplH)
    {
        err = 5;
        goto exit_function;
    }

    /* disconnect metal - ligand bonds */
    for (i = 0; i < num_atoms; i++)
    {
        if (!bMetal[i])
        {
            continue;
        }

        /* disconnect metal atom M

        Note: Defect in case of bridging ligands:

        M     M                          M     M             M     M(+)
        \  /   will be transformed to            , not to
        N(+)                             N(+)                N(-)
        / \                              / \                 / \
        R   R                            R   R               R   R

        Non-bridging are OK:

        M     R           M(+)  R
        \  /                 /
        N(+)    --->      N
        / \               / \
        R   R             R   R

        */

        for (j = at[i].valence - 1; 0 <= j; j--)
        {
            if (j < at[i].valence && !bMetal[(int) at[i].neighbor[j]])
            {
                /* do not break metal-metal bond here */

                num_disconnected += DisconnectOneLigand( at,
                                                         nOldCompNumber,
                                                         bMetal,
                                                         elnumber_Heteroat,
                                                         num_halogens,
                                                         num_atoms,
                                                         i,
                                                         j,
                                                         bTautFlagsDone );
            }
        }
    }

    /* disconnect metal-metal bonds */
    for (i = 0; i < num_atoms; i++)
    {
        if (!bMetal[i])
        {
            continue;
        }
        for (j = at[i].valence - 1; 0 <= j; j--)
        {
            if (j < at[i].valence && bMetal[(int) at[i].neighbor[j]])
            {
                /* break metal-metal bond here */

                num_disconnected += DisconnectOneLigand( at,
                                                         nOldCompNumber,
                                                         bMetal,
                                                         elnumber_Heteroat,
                                                         num_halogens,
                                                         num_atoms,
                                                         i,
                                                         j,
                                                         bTautFlagsDone );
            }
        }
    }

exit_function:

    if (!num_disconnected)
    {
        err = 6;
    }
    if (at && err)
    {
        inchi_free( at );
        at = NULL;
    }
    if (atom && at)
    {    /* changed if ( at ) to if ( atom && at ) 2004-04-03 */
        inchi_free( atom );
        atom = NULL;
    }
    if (bMetal)
        inchi_free( bMetal );

    if (at)
    {
        orig_inp_data->at = at;
        orig_inp_data->num_inp_atoms = num_at;
    }

    return err ? -err : num_disconnected;
}


/****************************************************************************/
int DisconnectOneLigand( inp_ATOM *at,
                         AT_NUMB *nOldCompNumber,
                         S_CHAR *bMetal,
                         char *elnumber_Heteroat,
                         int num_halogens,
                         int num_atoms,
                         int iMetal,
                         int jLigand,
                         INCHI_MODE *bTautFlagsDone )
{
    int i, j, iLigand, neigh, val;
    int metal_neigh_ord[MAXVAL], num_neigh_arom_bonds[MAXVAL];
    int num_metal_neigh, num_disconnections;
    int num_del_arom_bonds, num_tot_arom_bonds, new_charge;
    char *p;

    iLigand = at[iMetal].neighbor[jLigand];
    num_metal_neigh = 0;
    num_disconnections = 0;
    num_del_arom_bonds = num_tot_arom_bonds = 0;

    /* find bonds to disconnect */
    for (i = 0; i < at[iLigand].valence; i++)
    {
        num_neigh_arom_bonds[i] = 0;
        neigh = (int) at[iLigand].neighbor[i];
        if (neigh < num_atoms && bMetal[neigh])
        {
            metal_neigh_ord[num_metal_neigh++] = i;
            if (at[iLigand].bond_type[i] > BOND_TYPE_TRIPLE)
            {
                /* aromatic bond */
                for (j = 0; j < at[neigh].valence; j++)
                {
                    num_neigh_arom_bonds[i] += ( at[neigh].bond_type[j] > BOND_TYPE_TRIPLE );
                }
                num_del_arom_bonds++;
            }
        }
        num_tot_arom_bonds += ( at[iLigand].bond_type[i] > BOND_TYPE_TRIPLE );
    }

    /* Disconnect */
    if (num_del_arom_bonds)
    {
        /* fix chem_valence of the ligand and its neighbors in case of disconnecting arom. bonds */
        /* because in this case special care should be taken of updating at[].chem_bonds_valence */
        for (i = 0; i < num_metal_neigh; i++)
        {
            j = metal_neigh_ord[i];
            if (num_neigh_arom_bonds[j])
            {
                neigh = at[iLigand].neighbor[j];
                at[neigh].chem_bonds_valence -= num_neigh_arom_bonds[j] / 2 - ( num_neigh_arom_bonds[j] - 1 ) / 2;
            }
        }
        at[iLigand].chem_bonds_valence -= num_tot_arom_bonds / 2 - ( num_tot_arom_bonds - num_del_arom_bonds ) / 2;
    }

    /* disconnect in reverse order, otherwise the metal_neigh_ord[i]
    becomes invalid after the first disconnection
    */
    for (i = num_metal_neigh - 1; 0 <= i; i--)
    {
        num_disconnections += DisconnectInpAtBond( at,
                                                   nOldCompNumber,
                                                   iLigand,
                                                   metal_neigh_ord[i] );
    }

    /* attempt to change ligand charge to make its valence 'natural' */
    i = num_tot_arom_bonds - num_del_arom_bonds;
    if ((i && i != 2 && i != 3) ||
         (at[iLigand].radical && at[iLigand].radical != RADICAL_SINGLET) ||
         !( p = strchr( elnumber_Heteroat, at[iLigand].el_number ) )) /* djb-rwth: addressing LLVM warnings */
    {
        goto exit_function;  /* non-standard atom */
    }

    val = at[iLigand].chem_bonds_valence + NUMH( at, iLigand );
    new_charge = MAX_ATOMS; /* impossible value */

    if (!val)
    {
        if (p - elnumber_Heteroat < num_halogens)
        {
            new_charge = -1;
        }
    }
    else
    {
        for (i = -1; i <= 1; i++)
        {
            if (val == get_el_valence( at[iLigand].el_number, i, 0 ))
            {
                new_charge = i; /* found charge that fits chem. valence */
                break;
            }
        }
    }

    if (new_charge != MAX_ATOMS)
    {
        if (( new_charge != at[iLigand].charge ||
            ( at[iLigand].radical && at[iLigand].radical != RADICAL_SINGLET ) ) &&
             1 == num_metal_neigh)
        {
            if (1 == new_charge && 4 == val && 2 == at[iLigand].valence &&
                 4 == at[iLigand].chem_bonds_valence &&
                 at[iLigand].bond_type[0] == at[iLigand].bond_type[1])
            {
                ; /* do not add +1 charge to disconnected =N=, etc. 2004-10-27 */
            }
            else
            {
                if (bTautFlagsDone && new_charge != at[iLigand].charge)
                {
                    *bTautFlagsDone |= TG_FLAG_MOVE_CHARGE_COORD_DONE;
                }
                at[iMetal].charge -= new_charge - at[iLigand].charge;
                at[iLigand].charge = new_charge;
                /*at[iLigand].radical = 0;*/
            }
        }
    }

exit_function:

    return num_disconnections; /* ret;*/
}


/****************************************************************************/
double dist3D( inp_ATOM *at1, inp_ATOM *at2 )
{
    double dx = at1->x - at2->x;
    double dy = at1->y - at2->y;
    double dz = at1->z - at2->z;

    return sqrt( dx*dx + dy*dy + dz*dz );
}


/****************************************************************************/
#define MIN_BOND_LENGTH   (1.0e-6)
#define MIN_COS           (1.0e-6)
#define MIN_BOND_LENGTH2  (MIN_BOND_LENGTH*MIN_BOND_LENGTH)
#define MAX_BOND_LENGTH   (1.0e30)
/****************************************************************************/


/****************************************************************************/
double GetMinDistDistribution( inp_ATOM *at,
                               int num_at,
                               int iat,
                               int iat_H,
                               int bInAllComponents,
                               double min_dist[],
                               int num_segm )
{
    /*    const double one_pi = 2.0*atan2(1.0 , 0.0 ); */
    const double one_pi = 3.14159265358979323846; /* M_PI */
    const double two_pi = 2.0*one_pi;
    const double f_step = two_pi / num_segm;
    const double h_step = f_step / 2.0;
    int i, j, k, kk, ki, kn, n, num_bonds;
    double xi, yi, xn, yn, cross_prod_in, dot_prod_in, xni, yni, rni, tni, rmin;
    double fi, fk, fn, ft = 0, rt = 0, rk, ri, rn, c, ave_bond_len;

    for (i = 0; i < num_segm; i++)
    {
        min_dist[i] = MAX_BOND_LENGTH; /* more than any distance */
    }
    num_bonds = 0;
    ave_bond_len = 0.0;

    for (i = 0; i < num_at; i++)
    {
        if (i != iat && i != iat_H &&
            ( bInAllComponents || at[i].component == at[iat].component ))
        {
            for (j = 0; j < at[i].valence; j++)
            {
                n = at[i].neighbor[j];
                if (( n > i && n != iat ) || n == iat_H)
                {
                    continue;
                }
#if ( bRELEASE_VERSION != 1 && defined(_DEBUG) )
                if (n == iat)
                {
                    int stop = 1;  /* <BRKPT> */
                }
#endif
                xi = at[i].x - at[iat].x;  /* ri; i != iat */
                yi = at[i].y - at[iat].y;
                xn = at[n].x - at[iat].x;  /* rn; possibly n == iat */
                yn = at[n].y - at[iat].y;
                cross_prod_in = xi*yn - xn*yi; /* ((r(i)-r(iat)) x (r(n)-r(iat)) */
                if (cross_prod_in < -0.01*MIN_BOND_LENGTH2)
                {
                    /* make sure the r(i)->r(n) vector is counterclockwise around at[iat] */
                    inchi_swap( (char*) &xi, (char*) &xn, sizeof( xi ) );
                    inchi_swap( (char*) &yi, (char*) &yn, sizeof( yi ) );
                    /* djb-rwth: removing redundant code */
                }

                xni = xn - xi; /* r(n)->r(i) */
                yni = yn - yi;
                rni = xni*xni + yni*yni;
                if (rni > 0.01*MIN_BOND_LENGTH2)
                {
                    /* vector length |ri->rn| is not too small */
                    /* arrowhead of the vector r(t) = ri + (rn-ri)*t; 0 <= t <= 1 points to the bond ri->rn */
                    /* r(tni) is perpendicular to the bond ri->rn so that min|r(t)| = r(tni) = |tni|*rni */
                    tni = -( xni*xi + yni*yi ) / rni;
                    /* find min. distance from n-i bond to at[iat] */
                    if (tni < 0.0)
                    {
                        rmin = sqrt( xi*xi + yi*yi );
                    }
                    else if (tni > 1.0)
                    {
                        rmin = sqrt( xn*xn + yn*yn );
                    }
                    else
                    {
                        rmin = sqrt( tni*tni*rni );
                    }
                    ave_bond_len += sqrt( rni );
                    num_bonds++;
                }
                else
                {
                    /* zero length i-n bond */
                    tni = 0.5; /* fake */
                    rmin = sqrt( xi*xi + yi*yi ); /* arbitrarily choose one */
                }
                if (rmin >= 0.1*MIN_BOND_LENGTH)
                {
                    /* at[iat] does not belong to at[i]-at[n] bond */
                    int    bCalc_rt = 1;
                    fi = atan2( yi, xi );
                    fn = ( n == iat ) ? fi : atan2( yn, xn );
                    if (fi > fn)
                    {
                        /* make sure fn - fi >= 0 */
                        fn += two_pi;
                    }
                    if (fi < 0.0)
                    {
                        fi += two_pi;
                        fn += two_pi;
                    }
                    ki = (int) floor( ( fi + h_step ) / f_step );  /* cast does not match function type */
                    kn = (int) floor( ( fn + h_step ) / f_step );

                    /* the bond may affect several segments */
                    for (k = ki; k <= kn; k++)
                    {
                        kk = k % num_segm;
                        if (min_dist[kk] < rmin)
                        {
                            continue;
                        }
                        if (bCalc_rt)
                        {
                            if (n == iat)
                            {
                                ft = fi;
                                rt = rmin;
                            }
                            else
                            {
                                double xt, yt;
                                xt = xi + xni*tni;
                                yt = yi + yni*tni;
                                ft = atan2( yt, xt );
                                rt = sqrt( xt*xt + yt*yt );
                            }
                            bCalc_rt = 0;
                        }
                        fk = f_step * kk;
                        c = fabs( cos( fk - ft ) );
                        if (c < MIN_COS)
                            c = MIN_COS;
                        rk = rt / c;
                        if (min_dist[kk] > rk)
                        {
                            min_dist[kk] = rk;
                        }
                    }
                }
                else
                {
                    /* rmin < 0.1*MIN_BOND_LENGTH */
                    ri = xi*xi + yi*yi;
                    rn = xn*xn + yn*yn;
                    if (ri > MIN_BOND_LENGTH2 && rn > MIN_BOND_LENGTH2)
                    {
                        dot_prod_in = xn*xi + yn*yi;
                        /* a very short bond */
                        if (dot_prod_in > 0.01*MIN_BOND_LENGTH2)
                        {
                            /* bond does not cross at[iat] */
                            double fyixi = atan2( yi, xi );
                            if (fyixi < 0.0) fyixi += two_pi;
                            kk = (int) floor( ( fyixi + h_step ) / f_step ) % num_segm;
                            if (min_dist[kk] > rmin)
                            {
                                min_dist[kk] = rmin;
                            }
                        }
                        else if (dot_prod_in < -0.01*MIN_BOND_LENGTH2)
                        {
                            /* bond does cross at[iat] */
                            double fyixi = atan2( yi, xi );
                            if (fyixi < 0.0) fyixi += two_pi;
                            kk = (int) floor( ( fyixi + h_step ) / f_step ) % num_segm;
                            if (min_dist[kk] > rmin)
                            {
                                min_dist[kk] = rmin;
                            }
                            fyixi += one_pi;
                            kk = (int) floor( ( fyixi + h_step ) / f_step ) % num_segm;
                            if (min_dist[kk] > rmin)
                            {
                                min_dist[kk] = rmin;
                            }
                        }
                        else
                        {
                            ; /* error, should not happen */
                        }
                    }
                    else if (ri <= MIN_BOND_LENGTH2 && rn <= MIN_BOND_LENGTH2)
                    {
                        /* a very short bond coincides with at[iat]; ignore */
                        ;
                    }
                    else
                    {
                        /* one end of the bond coincides with at[iat] */
                        fi = ri > rn ? atan2( yi, xi ) : atan2( yn, xn );
                        if (fi < 0.0) fi += two_pi;
                        kk = (int) floor( ( fi + h_step ) / f_step ) % num_segm;
                        if (min_dist[kk] > rmin)
                        {
                            min_dist[kk] = rmin;
                        }
                    }
                }
            }
        }
    }

    if (num_bonds)
    {
        return  ave_bond_len / (double) num_bonds;
    }
    else
    {
        return 0.0;
    }
}


/****************************************************************************/
int move_explicit_Hcation( inp_ATOM *at,
                           int num_at,
                           int iat,
                           int iat_H,
                           int bInAllComponents )
{

#define NUM_SEGM 20

    /*    const double one_pi = 2.0*atan2(1.0 , 0.0 ); */
    const double one_pi = 3.14159265358979323846; /* M_PI */
    const double two_pi = 2.0*one_pi;
    const double f_step = two_pi / NUM_SEGM;
    const double h_step = f_step / 2.0;
    double min_dist[NUM_SEGM];
    int nB, i, k, kk, next, val = 0;
    double r, r0, xd, yd, zd, xr, yr, zr, ave_bond_len;
    /*double step = 4.0*atan(1.0)/NUM_SEGM;*/
    /* find at[iat] neighbors coordinates */

    xd = yd = zd = 0.0;

    if (at[iat].valence)
    {
        for (i = 0, nB = 0, r = 0.0; i < at[iat].valence; i++)
        {
            next = at[iat].neighbor[i];
            xd += at[next].x;
            yd += at[next].y;
            zd += at[next].z;
            r += dist3D( at + iat, at + next );
            nB++;
        }
        xd /= (double) nB;
        yd /= (double) nB;
        zd /= (double) nB;
        r /= (double) nB;
        r0 = sqrt( (double) ( xd - at[iat].x )*( xd - at[iat].x )
                   + (double) ( yd - at[iat].y )*( yd - at[iat].y ) );
    }
    else
    {
        if (at[iat_H].valence)
        {
            r = dist3D( at + iat_H, at + (int) at[iat_H].neighbor[0] );
        }
        else
        {
            r = 0.0;
        }
        r0 = 0.0;
    }

    ave_bond_len = GetMinDistDistribution( at, num_at, iat, iat_H,
                                           bInAllComponents, min_dist,
                                           NUM_SEGM );

    if (r < MIN_BOND_LENGTH && ave_bond_len > MIN_BOND_LENGTH)
    {
        r = ave_bond_len; /* ave_bond_len = 0.0 may mean that it is 0D structure */
    }

    if (r > MIN_BOND_LENGTH)
    {
        /* process non-zero bond lengths */
        double f;
        if (10.0*r0 < r)
        {
            xr = -r;     /* arbitrary */
            yr = 0.0;
            zr = 0.0;
        }
        else
        {
            /*
            if ( r0 < MIN_BOND_LENGTH ) {
            r0 = 1.0;
            }
            */
            xr = r * ( at[iat].x - xd ) / r0;
            yr = r * ( at[iat].y - yd ) / r0; /* length = r */
            zr = r * ( at[iat].z - zd ) / r0;

            /*          -- test: opposire direction --
            xr =   -r * ( at[iat].x - xd )/r0;
            yr =   -r * ( at[iat].y - yd )/r0;
            zr =   -r * ( at[iat].z - zd )/r0;
            */
            if (xr*xr + yr*yr < 0.04*r*r)
            {
                xr = -r;
                yr = 0.0;
            }
        }

        r = sqrt( xr*xr + yr*yr );
        f = atan2( yr, xr );

        if (f < 0.0)
        {
            f += two_pi;
        }

        kk = (int) floor( ( f + h_step ) / f_step ) % NUM_SEGM;
        /* cast does not match function type by design */

        if (min_dist[kk] < 1.5* r)
        {
            double dist = 1.5*r;
            int start = -1, len = 0, start_max = -1, len_max = 0;

        again:
            /* look for longest kk interval with min_dist[kk] >= dist */
            for (k = 0, start = 0, len = 0, len_max = 0; k < 2 * NUM_SEGM; k++)
            {
                kk = k % NUM_SEGM;
                if (min_dist[kk] >= dist)
                {
                    if (!len++)
                    {
                        start = k;
                    }
                }
                else
                {
                    if (len > len_max)
                    {
                        len_max = len;
                        start_max = start;
                    }
                    len = 0;
                }
            }
            if (!len_max)
            {
                if (dist > 0.1*r)
                {
                    dist *= 0.75;
                    goto again;
                }
                else
                {
                    goto done; /* do it anyway */
                }
            }
            else
            {
                /* found a good sector */
                f = f_step * ( (double)start_max +  ((double)len_max - 1.0 ) / 2.0 ); /* djb-rwth: cast operators added */
                r0 = dist / 1.5;
                xr = r0 * cos( f );
                yr = r0 * sin( f );
                zr = zr / r*r0;
            }
        }
    }
    else
    {
        xr = yr = zr = 0;
    }

done:
    if (at[iat_H].valence)
    {
        /* disconnect H */
        next = at[iat_H].neighbor[0];
        for (i = 0; i < at[next].valence; i++)
        {
            if (at[next].neighbor[i] == iat_H)
            {
                RemoveInpAtBond( at, next, i );
                i = 0; /* success */
                break;
            }
        }
    }
    else
    {
        /* isolated H+ cation */
        next = iat_H;
        i = 0;
        at[iat_H].valence = 1;
        at[iat_H].chem_bonds_valence = 1;
        at[iat_H].bond_type[0] = BOND_TYPE_SINGLE;
    }

    if (0 == i /*i < at[next].valence*/)
    {
        /* move charge */
        if (at[next].charge > 0 && at[iat].charge < 0)
        {
            at[next].charge--;
            at[iat].charge++;
        }

        /* connect H to at[iat] */
        val = at[iat].valence;
        
#pragma warning (push)
#pragma warning (disable: 6386)
        if (val < MAXVAL)
        {
            at[iat].neighbor[val] = iat_H;
            at[iat].bond_type[val] = at[iat_H].bond_type[0];
            at[iat].bond_stereo[val] = 0;
            at[iat].chem_bonds_valence += at[iat_H].bond_type[0];
            at[iat].valence = val + 1;
        };
#pragma warning (pop)

        at[iat_H].component = at[iat].component;
        at[iat_H].neighbor[0] = iat;
        at[iat_H].bond_stereo[0] = 0; /* possible loss of stereo info */
        at[iat_H].x = at[iat].x + xr;
        at[iat_H].y = at[iat].y + yr;
        at[iat_H].z = at[iat].z + zr;

        return 1; /* success */
    }

    return 0; /* failed */
}


/****************************************************************************/
int add_DT_to_num_H( int num_atoms, inp_ATOM *at )
/*  assume num_1H, num_D and num_T are not included in num_H */
{
    int i, j;
    for (i = 0; i < num_atoms; i++)
    {
        for (j = 0; j < NUM_H_ISOTOPES; j++)
        {
            at[i].num_H += at[i].num_iso_H[j];
        }
    }
    return 0;
}


/****************************************************************************
Return value: new number of atoms > 0 or -1=out of RAM
****************************************************************************/
int remove_terminal_HDT( int num_atoms, inp_ATOM *at, int bFixTermHChrg )
{
    AT_NUMB   *new_ord;
    inp_ATOM  *new_at;
    char *p;
    static const char szHDT[] = "HDT";
    static const int  kMax = sizeof( szHDT ); /*  = 4 */
    int ret = -1;
    int num_hydrogens = 0, num_H = 0;  /*  number of terminal H, D, T */
    int i, j, k, n, m;
    int val;
    AT_RANK new_HydrogenAt_order[NUM_H_ISOTOPES + 1];
    AT_RANK new_OtherNeigh_order[MAXVAL];
    S_CHAR  old_trans[MAX_NUM_STEREO_BONDS];

    int  num_OtherNeigh, num_HydrogenAt;

    new_ord = (AT_NUMB *) inchi_calloc( num_atoms, sizeof( new_ord[0] ) ); /* changed malloc to calloc 9-11-2003 */
    new_at = (inp_ATOM  *) inchi_malloc( sizeof( new_at[0] ) *num_atoms );
    if (!new_ord || !new_at)
    {
        goto exit_function;
    }

    /*  move H. D, T to the end of the list of atoms */
    for (i = 0; i < num_atoms; i++)
    {
        at[i].component = i; /*  temporarily save original numbering */
                             /*  get k = temp. hydrogen isotope/non-hydrogen atom type: */
                             /*  k=0:H, k=2:D, k=3:T, k=4=kMax: not a hydrogen */
        k = at[i].elname[1] ? kMax : ( p = (char*) strchr( szHDT, at[i].elname[0] ) ) ? (int) ( p - szHDT ) : kMax;
        /*  set hydrogen isotope atw differences */
        /*  Notes: k-value of isotopic H is incremented to correct iso_atw_diff value later. */
        /*         1H isotope cannot be detected here. */
        if (k == ATW_H || k == ATW_H + 1)
        {
            /* D or T, k = 1 or 2 */
            at[i].elname[0] = 'H'; /*  hydrogen isotope */
            at[i].iso_atw_diff = ++k; /*  increment k to make k = iso_atw_diff ( 2 for D, 3 for T ) */
        }
        num_H += ( k != kMax && at[i].valence == 1 && at[i].chem_bonds_valence == 1 && !NUMH( at, i ) );
    }

    /* special case: HD, HT, DT, HH: the only non-isotopic H or
    * the lightest isotopic H out of two is removed
    * to become implicit (make the heavier H the "central atom").
    * Note: This must be consistent with MOL_FMT_to_atom()
    * treatment of isotopic Hn aliases.
    */
    if (2 == num_H && 2 == num_atoms && !NUMH( at, 0 ) && !NUMH( at, 1 ))
    {

        if (at[0].iso_atw_diff >= at[1].iso_atw_diff)
        {
            new_ord[0] = 0;
            new_ord[1] = 1;
        }
        else
        {
            new_ord[0] = 1;
            new_ord[1] = 0;
        }
        if (at[new_ord[1]].charge)
        {
            at[new_ord[0]].charge += at[new_ord[1]].charge;
            at[new_ord[1]].charge = 0;
        }
        new_at[new_ord[0]] = at[0];
        new_at[new_ord[1]] = at[1];
        num_hydrogens = 1;
    }
    else
    {
        /* general case except H-H */
        for (i = 0; i < num_atoms; i++)
        {
            k = ( at[i].elname[1] || NUMH( at, i ) ) ? kMax : ( at[i].elname[0] == 'H' ) ? at[i].iso_atw_diff : kMax;
            if (k < kMax && at[i].valence == 1 && at[i].chem_bonds_valence == 1 &&
                 /*  the order of comparison is important */
                ( ( n = (int) at[i].neighbor[0] ) > i               /* at[n] has not been encountered yet*/ ||
                 (int) new_ord[n] < num_atoms - num_hydrogens ) /* at[n] might have been encountered; it has not been moved */)
            {
                /*  found an explicit terminal hydrogen */
                num_hydrogens++;
                if (k == 0 && ATW_H <= at[i].iso_atw_diff && at[i].iso_atw_diff < ATW_H + NUM_H_ISOTOPES)
                {
                    k = at[i].iso_atw_diff; /*  H isotope has already been marked above or elsewhere */ /* djb-rwth: ignoring LLVM warning: variable used */

                }
                if (at[i].charge)
                {
                    /*  transfer charge from the hydrogen */
                    at[n].charge += at[i].charge;
                    at[i].charge = 0;
                    if (bFixTermHChrg)
                    {
                        /*    Fixed bug (July 6, 2008 IPl) :
                        if terminal H was charged (not neutralized before call of remove_terminal_HDT)
                        and had an ordering number > than that of heavy-atom neighbour, then
                        charge on neighbour atom was not adjusted (though charge on H was removed). */
                        if (i > n)
                            /* new_at[new_ord[n]] has been created and filled already */
                            new_at[new_ord[n]].charge = at[n].charge;
                    }
                }
                new_ord[i] = num_atoms - num_hydrogens;  /*  move hydrogens to the end of the list */
            }
            else
            {
                /* atom is not an explicit terminal hydrogen */
                new_ord[i] = i - num_hydrogens;  /*  adjust non-hydrogens positions */
            }

            /*  copy atom to the new position */
            new_at[new_ord[i]] = at[i];
        } /* i */
    } /* general case except H-H */

    if (num_hydrogens)
    {
        int num_others = num_atoms - num_hydrogens; /*  atoms which are not terminal H, D, T */
        if (num_hydrogens > 1)
        {
            /*  sort hydrogen isotopes in ascending order, */
            /*  orig, numbers being the secondary sorting key */
            qsort( new_at + num_others, num_hydrogens, sizeof( new_at[0] ), cmp_iso_atw_diff_component_no );
        }
        /*  save new numbering of hydrogen atoms using temporarily saved orig numbering */
        for (i = num_others; i < num_atoms; i++)
        {
            new_ord[(int) new_at[i].component] = i;
        }

        /*  renumber neighbors according to new_ord[] and detach terminal hydrogens */
        for (i = 0; i < num_others; i++)
        {
            memset( new_HydrogenAt_order, 0, sizeof( new_HydrogenAt_order ) ); /* djb-rwth: memset_s C11/Annex K variant? */
            memset( new_OtherNeigh_order, 0, sizeof( new_OtherNeigh_order ) ); /* djb-rwth: memset_s C11/Annex K variant? */
            num_OtherNeigh = 0;
            num_HydrogenAt = 0;
            num_H = 0;

            for (m = 0; m < MAX_NUM_STEREO_BONDS && new_at[i].sb_parity[m]; m++)
            {
                old_trans[m] = 2 - ( new_at[i].sn_ord[m] + new_at[i].sb_ord[m] + ( new_at[i].sn_ord[m] > new_at[i].sb_ord[m] ) ) % 2;
            }

            for (k = val = 0; k < new_at[i].valence; k++) /* djb-rwth: removing redundant variables/code */
            {
                if (num_others <= ( n = new_ord[new_at[i].neighbor[k]] ))
                {
                    /*  discovered neighbor = disconnected explicit hydrogen
                    *  i = new atom new_at[i] ordering number
                    *  n = new number of the explicit H
                    *  k = ordering number of the explicit H in new_at[i] adjacency list
                    */
                    if (0 < new_at[n].iso_atw_diff && new_at[n].iso_atw_diff < ATW_H + NUM_H_ISOTOPES)
                    {
                        /* make explicit isotopic H implicit */
                        new_at[i].num_iso_H[new_at[n].iso_atw_diff - 1] ++; /*  isotopic H */
                        num_HydrogenAt += !new_HydrogenAt_order[new_at[n].iso_atw_diff];
                        new_HydrogenAt_order[new_at[n].iso_atw_diff] = k + 1;
                    }
                    else
                    {
                        /* make explicit non-isotopic H implicit */
                        new_at[i].num_H++; /*  non-isotopic H */
                        num_HydrogenAt += !num_H;
                        num_H++;
                        new_HydrogenAt_order[0] = k + 1;
                    }
                    /*  decrement chem. bonds valence because one bond is removed */
                    new_at[i].chem_bonds_valence = inchi_max( 0, new_at[i].chem_bonds_valence - 1 );
                    new_at[n].neighbor[0] = i; /*  update removed hydrogen neighbor number */
                    if (new_at[i].sb_parity[0])
                    {
                        /* if the removed H is an SB neighbor then mark it as removed */
                        for (m = 0; m < MAX_NUM_STEREO_BONDS && new_at[i].sb_parity[m]; m++)
                        {
                            if (k == (int) new_at[i].sn_ord[m])
                            {
                                new_at[i].sn_ord[m] = -( new_at[n].iso_atw_diff + 1 );
                                /* means the SB neighbor has been removed; (-4)=H, (-3)=1H, (-2)=D, (-1)=T */
                            }
                        }
                    }
                }
                else
                {
                    /* discovered a regular (not an explicit H) neighbor */
                    if (new_at[i].sb_parity[0])
                    {
                        if (num_OtherNeigh < MAX_NUM_STEREO_BONDS)
                        {
                            new_OtherNeigh_order[num_OtherNeigh] = k + 1;
                        }
                        num_OtherNeigh++; /* increment outside of if() to detect overflow */
                        if (val != k)
                        {
                            /* store new stereobond and sb-neighbor ordering numbers */
                            for (m = 0; m < MAX_NUM_STEREO_BONDS && new_at[i].sb_parity[m]; m++)
                            {
                                if (k == (int) new_at[i].sb_ord[m])
                                    new_at[i].sb_ord[m] = val;
                                else
                                    if (k == (int) new_at[i].sn_ord[m])
                                        new_at[i].sn_ord[m] = val;
                            }
                        }
                    }
                    new_at[i].neighbor[val] = new_ord[new_at[i].neighbor[k]];
                    new_at[i].bond_type[val] = new_at[i].bond_type[k];
                    new_at[i].bond_stereo[val] = new_at[i].bond_stereo[k];
                    val++;
                }
            }
            if (new_at[i].valence > val && new_at[i].sb_parity[0])
            {
                if (num_HydrogenAt == new_at[i].valence - val && num_HydrogenAt + num_OtherNeigh <= MAXVAL)
                {
                    /* recalculate parity so that it would describe neighbor sequence H,1H,D,T,neigh[0],neigh[1]... */
                    memmove(new_OtherNeigh_order + num_HydrogenAt, new_OtherNeigh_order, num_OtherNeigh * sizeof(new_OtherNeigh_order[0]));
                    for (k = 0, j = 1; k <= NUM_H_ISOTOPES; k++)
                    {
                        if (new_HydrogenAt_order[k] && (num_HydrogenAt - j < MAXVAL) && (num_HydrogenAt - j >= 0)) /* djb-rwth: fixing buffer overruns */
                        {
                            new_OtherNeigh_order[num_HydrogenAt - j] = new_HydrogenAt_order[k];
                            for (m = 0; m < MAX_NUM_STEREO_BONDS && new_at[i].sb_parity[m]; m++)
                            {
                                if ((int) new_at[i].sn_ord[m] == -( k + 1 ))
                                {
                                    new_at[i].sn_ord[m] = -j;
                                    /* negative means explicit H isotope ord are
                                    (contiguously) in front of the adjacency list */
                                }
                            }
                            j++;
                        }
                    }
                    /* at this point new_OtherNeigh_order[] contains
                    incremented old ordering numbers in new order */
                    k = insertions_sort_AT_RANK( new_OtherNeigh_order, num_HydrogenAt + num_OtherNeigh ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
                    /* djb-rwth: removing redundant code */
                               /*if ( k ) {*/
                               /*
                               for ( m = 0; m < MAX_NUM_STEREO_BONDS && new_at[i].sb_parity[m]; m ++ ) {
                               if ( PARITY_WELL_DEF(new_at[i].sb_parity[m]) ) {
                               if ( old_trans[m] != 2 - (4 + new_at[i].sn_ord[m] + new_at[i].sb_ord[m] + (new_at[i].sn_ord[m] > new_at[i].sb_ord[m]))%2 ) {
                               new_at[i].sb_parity[m] = 3 - new_at[i].sb_parity[m];
                               }
                               }
                               }
                               */
                               /*}*/
                }
            }
            new_at[i].valence = val;
        }
        memcpy(at, new_at, sizeof(at[0])* num_atoms);
        ret = num_others;
    }
    else
    {
        ret = num_atoms;
    }

exit_function:

    if (new_ord)
    {
        inchi_free( new_ord );
    }
    if (new_at)
    {
        inchi_free( new_at );
    }

    return ret;
}


/*#endif*/ /* } DISCONNECT_SALTS */

typedef enum tagIonAtomType
{
    IAT_H = 0,
    IAT_C,
    IAT_N,
    IAT_P,
    IAT_O,
    IAT_S,
    IAT_Se,
    IAT_Te,
    IAT_F,
    IAT_Cl,
    IAT_Br,
    IAT_I,
    IAT_MAX = 12
} ION_ATOM_TYPE;

/****************************************************************************/
int get_iat_number( int el_number )
{
    switch (el_number) {
        case EL_NUMBER_H:  return IAT_H;
        case EL_NUMBER_C:  return IAT_C;
        case EL_NUMBER_N:  return IAT_N;
        case EL_NUMBER_P:  return IAT_P;
        case EL_NUMBER_O:  return IAT_O;
        case EL_NUMBER_S:  return IAT_S;
        case EL_NUMBER_SE: return IAT_Se;
        case EL_NUMBER_TE: return IAT_Te;
        case EL_NUMBER_F:  return IAT_F;
        case EL_NUMBER_CL: return IAT_Cl;
        case EL_NUMBER_BR: return IAT_Br;
        case EL_NUMBER_I:  return IAT_I;
        default: return -1;
    }
}

#if ( READ_INCHI_STRING == 1 )


/****************************************************************************/
int bHeteroAtomMayHaveXchgIsoH( inp_ATOM *atom, int iat )
{
    inp_ATOM *at = atom + iat, *at2;
    int j, val, is_H = 0, num_H, iat_numb, bAccept; /* djb-rwth: removing redundant variables */

    if (0 > ( iat_numb = get_iat_number( at->el_number ) ))
    {
        return 0;
    }

    if (abs( at->charge ) > 1 || (at->radical && RADICAL_SINGLET != at->radical)) /* djb-rwth: addressing LLVM warning */
    {
        return 0;
    }

    val = -1;
    switch (iat_numb)
    {
        case IAT_N:
        case IAT_P:
            /* djb-rwth: removing redundant code */
            val = 3 + at->charge;
            break;

        case IAT_O:
        case IAT_S:
        case IAT_Se:
        case IAT_Te:
            /* djb-rwth: removing redundant code */
            val = 2 + at->charge;
            break;

        case IAT_F:
        case IAT_Cl:
        case IAT_Br:
        case IAT_I:
            if (at->charge == 0)
            {
                /* djb-rwth: removing redundant code */
                val = 1;
            }
            break;

        case IAT_H:
            if (at->valence == 0 &&
                 at->charge == 1)
            {
                is_H = 1; /* isolated proton */
                val = 0;
            }
    }
    if (val < 0)
    {
        return 0;
    }
    num_H = NUMH( at, 0 );
    if (val != at->chem_bonds_valence + num_H)
    {
        return 0;
    }
    if (is_H)
    {
        return 2; /* H atom */
    }
    else
    {
        /* djb-rwth: removing redundant code */
        for (j = 0, bAccept = 1; j < at->valence && bAccept; j++)
        {
            at2 = atom + (int) at->neighbor[j];
            if ((at2->charge && at->charge) ||
                ( at2->radical && RADICAL_SINGLET != at2->radical )) /* djb-rwth: addressing LLVM warning */
            {
                return 0; /* adjacent charged/radical atoms: do not neutralizate */
            }
        }
    }

    return 1;
}


#endif


/****************************************************************************/
int bNumHeterAtomHasIsotopicH( inp_ATOM *atom, int num_atoms )
{
    int i, j, val, is_H = 0, num_H, iat_numb, bAccept, num_iso_H, cur_num_iso_H, num_iso_atoms; /* djb-rwth: removing redundant variables */
    inp_ATOM *at, *at2;

    num_iso_H = 0;
    num_iso_atoms = 0;

    for (i = 0, at = atom; i < num_atoms; i++, at++)
    {

        num_iso_atoms += ( at->iso_atw_diff != 0 || NUM_ISO_H( at, 0 ) );
        /* isotopic atoms and implicit isotopic H */

        if (0 >( iat_numb = get_iat_number( at->el_number ) ))
        {
            continue;
        }

        if (abs( at->charge ) > 1 || (at->radical && RADICAL_SINGLET != at->radical)) /* djb-rwth: addressing LLVM warning */
        {
            continue;
        }

        val = -1;
        switch (iat_numb)
        {
            case IAT_N:
            case IAT_P:
                /* djb-rwth: removing redundant code */
                val = 3 + at->charge;
                break;

            case IAT_O:
            case IAT_S:
            case IAT_Se:
            case IAT_Te:
                /* djb-rwth: removing redundant code */
                val = 2 + at->charge;
                break;

            case IAT_F:
            case IAT_Cl:
            case IAT_Br:
            case IAT_I:
                if (at->charge == 0)
                {
                    /* djb-rwth: removing redundant code */
                    val = 1;
                }
                break;

            case IAT_H:
                if (at->valence == 0 &&
                     at->charge == 1)
                {
                    is_H = 1; /* isolated proton */
                    val = 0;
                }
        }
        if (val < 0)
        {
            continue;
        }

        num_H = NUMH( at, 0 );
        if (val != at->chem_bonds_valence + num_H)
        {
            continue;
        }

        if (is_H)
        {
            bAccept = 1;
            cur_num_iso_H = ( at->iso_atw_diff != 0 );
        }
        else
        {
            cur_num_iso_H = 0;
            for (j = 0, bAccept = 1; j < at->valence && bAccept; j++)
            {
                at2 = atom + (int) at->neighbor[j];
                if ((at2->charge && at->charge) ||
                    ( at2->radical && RADICAL_SINGLET != at2->radical )) /* djb-rwth: addressing LLVM warning */
                {
                    bAccept = 0; /* adjacent charged/radical atoms: do not neutralizate */
                    break;
                }
                else if (at2->el_number == EL_NUMBER_H &&
                          at2->valence == 1 && at2->iso_atw_diff)
                {
                    cur_num_iso_H++; /* isotopic explicit H */
                }
            }

            if (bAccept)
            {
                num_iso_atoms -= cur_num_iso_H;  /* avoid counting explicit H as isotopic atom */
                cur_num_iso_H += NUM_ISO_H( at, 0 );
            }
        }

        num_iso_H += ( bAccept && cur_num_iso_H ); /* number of acceptable heteroatoms that have isotopic H */
    }

    return
        ( ( num_iso_H ? 1 : 0 ) | ( num_iso_atoms ? 2 : 0 ) );
}


/****************************************************/
/* Mark and count disconnected structure components */
/* by Depth-first searching each component          */
/****************************************************/
int cmp_components( const void *a1, const void *a2 )
{
    int ret;
    AT_NUMB n1;
    AT_NUMB n2;

    n1 = ( (const AT_NUMB *) a1 )[0];
    /* number of atoms in the component -- descending order */
    n2 = ( (const AT_NUMB *) a2 )[0];

    if ((ret = (int) n2 - (int) n1)) /* djb-rwth: addressing LLVM warning */
    {
        return ret;
    }

    /* stable sort */
    n1 = ( (const AT_NUMB *) a1 )[1];
    /* component ordering number -- ascending order */
    n2 = ( (const AT_NUMB *) a2 )[1];

    ret = (int) n1 - (int) n2;

    return ret;
}


/****************************************************************************
Set the (disconnected) component numbers in ORIG_ATOM_DATA 'at[*].component'
NB: components are (stable) sorted by number of heavy atoms
****************************************************************************/
int MarkDisconnectedComponents( ORIG_ATOM_DATA *orig_at_data,
                                int bProcessOldCompNumbers )
{
    typedef AT_NUMB AT_TRIPLE[3];

    inp_ATOM    *at = orig_at_data->at;
    int         num_at = orig_at_data->num_inp_atoms;
    AT_NUMB     *nCurAtLen = NULL;

    AT_NUMB *nNewCompNumber = NULL;
    AT_NUMB *nPrevAtom = NULL;
    S_CHAR  *iNeigh = NULL;

    AT_NUMB *nOldCompNumber = NULL;
    int i, j, num_components, ret;
    int new_comp_no;
    AT_NUMB old_comp_no, another_comp_no, no_component;

    /* component_nbr[i][0] = number of atoms in the component i-1
    * component_nbr[i][1] = original component number (id-1) = i
    * after sorting:
    * component_nbr[j][2] = new number of component #(component_nbr[i][1]+1)
    */
    AT_TRIPLE *component_nbr = NULL;
    int fst_at, nxt_at, cur_at, cur_neq_fst;  /* moved from below 2024-09-01 DT */

    /* initialize */
    if (bProcessOldCompNumbers && !orig_at_data->nOldCompNumber)
    {
        bProcessOldCompNumbers = 0;
    }
    num_components = 0;

    /*
    for ( j = 0; j < num_at; j ++ )
    {
    at[j].component = 0;
    }
    */

    ret = -1;
    if (!num_at)
    {
        return 0;
    }

    nNewCompNumber = (AT_NUMB*)inchi_calloc(num_at, sizeof(nNewCompNumber[0]));
    nPrevAtom = (AT_NUMB*)inchi_calloc(num_at, sizeof(nPrevAtom[0]));
    iNeigh = (S_CHAR*)inchi_calloc(num_at, sizeof(iNeigh[0]));

    if (!nNewCompNumber || !nPrevAtom || !iNeigh) /* nNewCompNumber: for non-recursive DFS only: */ 
    {
        goto exit_function;
    }

    /*printf("\nnum_at = %d\n", num_at);*/

    /* Mark and count; avoid deep DFS recursion: it may make verifying software unhappy */
    /* nNewCompNumber[i] will contain new component number for atoms at[i], i=0..num_at-1 */
    
    for (j = 0; j < num_at; j++)
    {
        if (!nNewCompNumber[j])
        {
            /* mark starting with at[j] */
            fst_at = 0;
            nxt_at = 0;
            cur_at = j;
            cur_neq_fst = 1;
            num_components++;

            /* first time at at[j] */
            fst_at = cur_at;
            nNewCompNumber[fst_at] = (AT_NUMB) num_components;

            /* find next neighbor */
            do
            {
                if (iNeigh[cur_at] < at[cur_at].valence)
                {
                    int ineigh_incr = (int)iNeigh[cur_at];
                    nxt_at = at[cur_at].neighbor[ineigh_incr];
                    iNeigh[cur_at]++;

                    if (!nNewCompNumber[nxt_at])
                    {
                        /* forward edge: found new atom */
                        nNewCompNumber[nxt_at] = (AT_NUMB) num_components;
                        nPrevAtom[nxt_at] = (AT_NUMB) cur_at;
                        cur_at = nxt_at;
                    }
                }
                else if (cur_at == fst_at)
                {
                    cur_neq_fst = 0;
                    /* break;  done */
                }
                else
                {
                    cur_at = nPrevAtom[cur_at]; /* retract */
                }
            } while (cur_neq_fst);
        }
    }

    inchi_free( nPrevAtom );
    nPrevAtom = NULL;
    inchi_free( iNeigh );
    iNeigh = NULL;

    /* Allocate more memory */
    i = inchi_max( num_components, orig_at_data->num_components );

    nCurAtLen = (AT_NUMB*)inchi_calloc((long long)num_components + 1, sizeof(nCurAtLen[0])); /* djb-rwth: cast operator added */
    nOldCompNumber = (AT_NUMB*)inchi_calloc((long long)i + 1, sizeof(nOldCompNumber[0])); /* djb-rwth: cast operator added */
    component_nbr = (AT_TRIPLE*)inchi_calloc((long long)num_components + 1, sizeof(component_nbr[0])); /* djb-rwth: cast operator added */

    if (!nCurAtLen || !nOldCompNumber || !component_nbr)
    {
        goto exit_function;
    }

    /* Count atoms per component and renumber the components */
    for (i = 0; i < num_components; i++)
    {
        component_nbr[i][0] = 0; /* number of atoms in the component */
        component_nbr[i][1] = i; /* component ordering number */
    }

    for (j = 0; j < num_at; j++)
    {
        component_nbr[(int) nNewCompNumber[j] - 1][0] ++; /* count atoms in each component */
    }

    /* Sort settings
    key: number of atoms
    order: descending
    stable sort
    */

    qsort( (void*) component_nbr[0], num_components, sizeof( component_nbr[0] ), cmp_components ); /* djb-rwth: fixed buffer overrun */

    /* Invert the transposition */
    for (i = 0; i < num_components; i++)
    {
        nCurAtLen[i] = component_nbr[i][0];
        component_nbr[component_nbr[i][1]][2] = i + 1;
    }

    /* Renumber the components so that the component with the greatest number of atoms is the first */
    no_component = num_at + 1;

    for (j = 0; j < num_at; j++)
    {
        /* new component number for at[j] */
        new_comp_no = component_nbr[(int) nNewCompNumber[j] - 1][2] - 1; /* starts from 0 */
        if (bProcessOldCompNumbers)
        {
            /* old component number for at[j] */
            old_comp_no = at[j].component;
            /* fill out nOldCompNumber[]; initially it contains zeroes */
            if (!old_comp_no)
            {
                nOldCompNumber[new_comp_no] = no_component; /* atom did not have component number */
            }
            else if (nOldCompNumber[new_comp_no] != old_comp_no)
            {
                if (!nOldCompNumber[new_comp_no])
                {
                    nOldCompNumber[new_comp_no] = old_comp_no;
                }
                else
                {
                    /* at[j] moved from old comp #old_comp_no to old comp #nOldCompNumber[new_comp_no]
                    Both components cannot be equal to any current component */
                    another_comp_no = nOldCompNumber[new_comp_no];
                    for (i = 0; i < num_components; i++)
                    {
                        if (nOldCompNumber[i] == old_comp_no ||
                             nOldCompNumber[i] == another_comp_no)
                        {
                            nOldCompNumber[i] = no_component;
                        }
                    }
                    /* nOldCompNumber[new_comp_no] = num_at+1; */
                }
            }
        }

        /* orig_at_data->nOldCompNumber */

        /* Finally, set the new component number for atom j (NB: starts from 1 ) */
        at[j].component = new_comp_no + 1;
    }

    if (bProcessOldCompNumbers)
    {
        for (j = 0; j < num_components; j++)
        {
            if (nOldCompNumber[j] == no_component)
            {
                /* the component has atom from another component */
                nOldCompNumber[j] = 0;
            }
            else if (nOldCompNumber[j] &&
                      !orig_at_data->nOldCompNumber[nOldCompNumber[j] - 1])
            {
                /* the component has changed in the previous processing  */
                nOldCompNumber[j] = 0;
            }
        }
    }
    else
    {
        for (j = 0; j < num_components; j++)
        {
            nOldCompNumber[j] = j + 1;
        }
    }

    ret = num_components;

exit_function:

    if (nNewCompNumber)
    {
        inchi_free( nNewCompNumber );
    }
    if (component_nbr)
    {
        inchi_free( component_nbr );
    }

    if (ret < 0)
    {
        if (nPrevAtom)
        {
            inchi_free( nPrevAtom );
            nPrevAtom = NULL;
        }
        if (iNeigh)
        {
            inchi_free( iNeigh );
            iNeigh = NULL;
        }
        if (nCurAtLen)
        {
            inchi_free( nCurAtLen );
            nCurAtLen = NULL;
        }
        if (nOldCompNumber)
        {
            inchi_free( nOldCompNumber );
            nOldCompNumber = NULL;
        }
        num_components = ret;
    }

    /* avoid memory leaks */
    if (orig_at_data->nCurAtLen)
    {
        inchi_free( orig_at_data->nCurAtLen );
    }
    if (orig_at_data->nOldCompNumber)
    {
        inchi_free( orig_at_data->nOldCompNumber );
    }

    orig_at_data->nCurAtLen = nCurAtLen;
    orig_at_data->nOldCompNumber = nOldCompNumber;

    orig_at_data->num_components = num_components;

    return ret;  /* number of disconnected components;
              1=>single connected structure        */
}


/****************************************************************************
Extract one (connected) component
****************************************************************************/
int ExtractConnectedComponent( inp_ATOM *at,
                               int num_at,
                               int component_number,
                               inp_ATOM *component_at )
{
    int i, j, num_component_at;
    AT_NUMB *number;

    if (NULL == ( number = (AT_NUMB*) inchi_calloc( num_at, sizeof( AT_NUMB ) ) ))
    {
        return CT_OUT_OF_RAM; /* out of memory */  /*   <BRKPT> */
    }

    /* copy atoms */
    for (i = 0, num_component_at = 0; i < num_at; i++)
    {
        if (at[i].component == component_number)
        {
            number[i] = num_component_at;
            component_at[num_component_at++] = at[i];
        }
    }

    /* renumber neighbors */
    for (i = 0; i < num_component_at; i++)
    {
        component_at[i].orig_compt_at_numb = (AT_NUMB) ( i + 1 );
        for (j = 0; j < component_at[i].valence; j++)
        {
            component_at[i].neighbor[j] = number[(int) component_at[i].neighbor[j]];
        }
    }

    inchi_free( number );

    return num_component_at;
}


/****************************************************************************/
int SetConnectedComponentNumber( inp_ATOM *at, int num_at, int component_number )
{
    int i;
    for (i = 0; i < num_at; i++)
    {
        at[i].component = (AT_NUMB) component_number;
    }

    return 0;
}


/****************************************************************************/
int Free_INChI_Stereo( INChI_Stereo *pINChI_Stereo )
{
    if (pINChI_Stereo)
    {
        qzfree( pINChI_Stereo->nNumber );
        qzfree( pINChI_Stereo->t_parity );
        qzfree( pINChI_Stereo->nNumberInv );
        qzfree( pINChI_Stereo->t_parityInv );
        qzfree( pINChI_Stereo->nBondAtom1 );
        qzfree( pINChI_Stereo->nBondAtom2 );
        qzfree( pINChI_Stereo->b_parity );
    }

    return 0;
}


/****************************************************************************/
INChI_Stereo *Alloc_INChI_Stereo( int num_at, int num_bonds )
{

    INChI_Stereo *pINChI_Stereo = (INChI_Stereo *)
        inchi_calloc( 1, sizeof( INChI_Stereo ) );

    if (pINChI_Stereo)
    {
        if (num_at &&
            ( pINChI_Stereo->nNumber = (AT_NUMB *) inchi_calloc( num_at, sizeof( pINChI_Stereo->nNumber[0] ) ) ) &&
             ( pINChI_Stereo->t_parity = (S_CHAR  *) inchi_calloc( num_at, sizeof( pINChI_Stereo->t_parity[0] ) ) ) &&
             ( pINChI_Stereo->nNumberInv = (AT_NUMB *) inchi_calloc( num_at, sizeof( pINChI_Stereo->nNumberInv[0] ) ) ) &&
             ( pINChI_Stereo->t_parityInv = (S_CHAR  *) inchi_calloc( num_at, sizeof( pINChI_Stereo->t_parityInv[0] ) ) ))
        {
            ;
        }
        else if (num_at)
        {
            goto out_of_RAM;
        }

        if (num_bonds &&
            ( pINChI_Stereo->nBondAtom1 = (AT_NUMB *) inchi_calloc( num_bonds, sizeof( pINChI_Stereo->nBondAtom1[0] ) ) ) &&
             ( pINChI_Stereo->nBondAtom2 = (AT_NUMB *) inchi_calloc( num_bonds, sizeof( pINChI_Stereo->nBondAtom2[0] ) ) ) &&
             ( pINChI_Stereo->b_parity = (S_CHAR  *) inchi_calloc( num_bonds, sizeof( pINChI_Stereo->b_parity[0] ) ) ))
        {
            ;
        }
        else if (num_bonds)
        {
            goto out_of_RAM;
        }

        return pINChI_Stereo;

    out_of_RAM:

        Free_INChI_Stereo( pINChI_Stereo );
        qzfree( pINChI_Stereo );
    } /* if ( pINChI_Stereo )  */

    return NULL;
}


/****************************************************************************/
int Free_INChI( INChI **ppINChI )
{

    INChI *pINChI;

    if ((pINChI = *ppINChI)) /* djb-rwth: addressing LLVM warning */
    {

#if ( bREUSE_INCHI == 1 )
        if (pINChI->nRefCount-- > 0)
            return 1;
#endif

        Free_INChI_Members( pINChI );
        qzfree( pINChI );
        *ppINChI = NULL;
    }

    return 0;
}


/****************************************************************************/
int Free_INChI_Members( INChI *pINChI )
{
    if (pINChI)
    {
        Free_INChI_Stereo(pINChI->Stereo);
        Free_INChI_Stereo(pINChI->StereoIsotopic);
        qzfree(pINChI->nAtom);
        qzfree(pINChI->nConnTable);
        qzfree(pINChI->nTautomer);
        qzfree(pINChI->nNum_H);
        qzfree(pINChI->nNum_H_fixed);
        qzfree(pINChI->IsotopicAtom);
        qzfree(pINChI->IsotopicTGroup);
        qzfree(pINChI->nPossibleLocationsOfIsotopicH);
        qzfree( pINChI->Stereo );       
        qzfree( pINChI->StereoIsotopic );
        qzfree( pINChI->szHillFormula );
    }

    return 0;
}


/****************************************************************************/
INChI *Alloc_INChI( inp_ATOM *at,
                    int num_at,
                    int *found_num_bonds,
                    int *found_num_isotopic,
                    int nAllocMode )
{
    int    i, num_bonds, num_isotopic_atoms;
    INChI  *pINChI;
    int    bIsotopic = ( nAllocMode & REQ_MODE_ISO );
    /* int    bTautomeric = (nAllocMode & REQ_MODE_TAUT); */

    if (num_at <= 0 ||
         NULL == ( pINChI = (INChI *) inchi_calloc( 1, sizeof( INChI ) ) ))
    {
        return NULL;
    }

    for (i = 0, num_bonds = 0, num_isotopic_atoms = 0; i < num_at; i++)
    {
        num_bonds += at[i].valence;
        /* if ( bIsotopic ) { */
        num_isotopic_atoms += ( 0 != at[i].iso_atw_diff ||
                                !strcmp( at[i].elname, "D" ) ||
                                !strcmp( at[i].elname, "T" ) ||
                                at[i].num_iso_H[0] ||
                                at[i].num_iso_H[1] ||
                                at[i].num_iso_H[2] );
        /* } */
    }
    num_bonds /= 2;

    *found_num_bonds = num_bonds;
    *found_num_isotopic = num_isotopic_atoms;

    if (( pINChI->nAtom = (U_CHAR*) inchi_calloc( num_at, sizeof( pINChI->nAtom[0] ) ) ) &&
        ( pINChI->nConnTable = (AT_NUMB*) inchi_calloc( (long long)num_at + (long long)num_bonds, sizeof( pINChI->nConnTable[0] ) ) ) && /* djb-rwth: cast operator added */
         ( pINChI->nTautomer = (AT_NUMB*) inchi_calloc( ( ( 3 + INCHI_T_NUM_MOVABLE )*(long long)num_at ) / 2 + 1, sizeof( pINChI->nTautomer[0] ) ) ) && /* djb-rwth: cast operator added */
         ( pINChI->nNum_H = (S_CHAR*) inchi_calloc( num_at, sizeof( pINChI->nNum_H[0] ) ) ) &&
         ( pINChI->nNum_H_fixed = (S_CHAR*) inchi_calloc( num_at, sizeof( pINChI->nNum_H_fixed[0] ) ) ))
    {
        ;
        /* nTautomer length: max. number of tautomeric groups is num_at/2

        1 word                     -> number of t-groups

        each group has:

        1 word                     -> number of endpoints+INCHI_T_NUM_MOVABLE
        INCHI_T_NUM_MOVABLE words   -> number(s) of moveable attachments
        numbers of endpoints words -> canon. numbers

        max. occurs if each t-group has 2 atoms (num_at/2 t-groups) and all atoms
        belong to t-groups (num_at endpoints)

        Total: 1 + (number of t-groups)*(1+INCHI_T_NUM_MOVABLE) + (number of endpoints) <=
        1 + (num_at/2) * (1+INCHI_T_NUM_MOVABLE) + num_at <=
        1 + (3+INCHI_T_NUM_MOVABLE)*num_at/2 words.
        */
    }
    else
    {
        goto out_of_RAM;
    }

    pINChI->szHillFormula = NULL; /*  the length is unknown */

    if (bIsotopic)
    {
        if (num_isotopic_atoms &&
            ( pINChI->IsotopicAtom = (INChI_IsotopicAtom *) inchi_calloc( num_isotopic_atoms, sizeof( INChI_IsotopicAtom ) ) ) &&
             ( pINChI->IsotopicTGroup = (INChI_IsotopicTGroup *) inchi_calloc( num_isotopic_atoms, sizeof( INChI_IsotopicTGroup ) ) ))
        {
            ;
        }
        else if (num_isotopic_atoms)
        {
            goto out_of_RAM;
        }
        if (!( pINChI->nPossibleLocationsOfIsotopicH = (AT_NUMB *) inchi_calloc( (long long)num_at + 1, sizeof( pINChI->nPossibleLocationsOfIsotopicH[0] ) ) )) /* djb-rwth: cast operator added */
        {
            goto out_of_RAM;
        }
    }

    if (( pINChI->Stereo = Alloc_INChI_Stereo( num_at, num_bonds ) ))
    {
        ;
    }
    else
    {
        goto out_of_RAM;
    }

    if (bIsotopic)
    {
        if (( pINChI->StereoIsotopic = Alloc_INChI_Stereo( num_at, num_bonds ) ))
        {
            ;
        }
        else
        {
            goto out_of_RAM;
        }
    }

    return pINChI;

out_of_RAM:
    if (pINChI)
    {
        Free_INChI( &pINChI );
        /*
        inchi_free(pINChI);
        */
    }

    return NULL;
}


/****************************************************************************/
int Free_INChI_Aux( INChI_Aux **ppINChI_Aux )
{
    INChI_Aux *pINChI_Aux = *ppINChI_Aux;
    if (pINChI_Aux)
    {

#if ( bREUSE_INCHI == 1 )
        if (pINChI_Aux->nRefCount-- > 0)
            return 1;
#endif

        qzfree( pINChI_Aux->nOrigAtNosInCanonOrd );
        qzfree( pINChI_Aux->nIsotopicOrigAtNosInCanonOrd );
        qzfree( pINChI_Aux->nOrigAtNosInCanonOrdInv );
        qzfree( pINChI_Aux->nIsotopicOrigAtNosInCanonOrdInv );
        qzfree( pINChI_Aux->szOrigCoord );
        qzfree( pINChI_Aux->OrigInfo );
        /*
        qzfree( pINChI_Aux->nOriginalAtomNumber          );
        qzfree( pINChI_Aux->nCanonicalTGroupNumbers      );
        qzfree( pINChI_Aux->nIsotopicCanonicalTGroupNumbers);
        qzfree( pINChI_Aux->nTautomer                    );
        qzfree( pINChI_Aux->nNontautomericCanonicalNumbers         );
        qzfree( pINChI_Aux->nIsotopicCanonicalNumbers    );
        qzfree( pINChI_Aux->nNontautomericIsotopicCanonicalNumbers );
        qzfree( pINChI_Aux->nNontautomericEquNumbers               );
        qzfree( pINChI_Aux->nNontautomericIsotopicEquNumbers       );
        */
        qzfree( pINChI_Aux->nConstitEquNumbers );
        qzfree( pINChI_Aux->nConstitEquTGroupNumbers );
        qzfree( pINChI_Aux->nConstitEquIsotopicNumbers );
        qzfree( pINChI_Aux->nConstitEquIsotopicTGroupNumbers );
        qzfree( pINChI_Aux );
        *ppINChI_Aux = NULL;
    }

    return 0;
}


/****************************************************************************/
INChI_Aux *Alloc_INChI_Aux( int num_at,
                            int num_isotopic_atoms,
                            int nAllocMode,
                            int bOrigCoord )
{
    INChI_Aux     *pINChI_Aux;
    int    bIsotopic = ( nAllocMode & REQ_MODE_ISO );
    int    num_at_tg = num_at + num_at / 2;
    /* int    bTautomeric = (nAllocMode & REQ_MODE_TAUT); */

    if (num_at <= 0 ||
         NULL == ( pINChI_Aux = (INChI_Aux *) inchi_calloc( 1, sizeof( INChI_Aux ) ) ))
    {
        return NULL;
    }

    if (( pINChI_Aux->nOrigAtNosInCanonOrd = (AT_NUMB*)
          inchi_calloc( num_at_tg, sizeof( pINChI_Aux->nOrigAtNosInCanonOrd[0] ) ) ) &&
          ( pINChI_Aux->nOrigAtNosInCanonOrdInv = (AT_NUMB*)
            inchi_calloc( num_at_tg, sizeof( pINChI_Aux->nOrigAtNosInCanonOrd[0] ) ) ) &&
            ( pINChI_Aux->nConstitEquNumbers = (AT_NUMB*)
              inchi_calloc( num_at_tg, sizeof( pINChI_Aux->nConstitEquNumbers[0] ) ) ))
    {
        ;
    }
    else
    {
        goto out_of_RAM;
    }

    if (num_at > 1 &&
        ( pINChI_Aux->nConstitEquTGroupNumbers = (AT_NUMB*) inchi_calloc( (long long)num_at / 2 + 1, sizeof( pINChI_Aux->nConstitEquTGroupNumbers[0] ) ) )) /* djb-rwth: cast operator added */
    {
        ;
    }
    else
    {
        if (num_at > 1)
        {
            goto out_of_RAM;
        }
    }

    if (num_at > 0)
    {
        pINChI_Aux->OrigInfo = (ORIG_INFO *) inchi_calloc( num_at, sizeof( pINChI_Aux->OrigInfo[0] ) );
        if (!pINChI_Aux->OrigInfo)
            goto out_of_RAM;
    }

    if (bOrigCoord && num_at > 0)
    {
        pINChI_Aux->szOrigCoord = (MOL_COORD *) inchi_calloc( num_at, sizeof( pINChI_Aux->szOrigCoord[0] ) );
        if (!pINChI_Aux->szOrigCoord)
            goto out_of_RAM;
    }

    if (bIsotopic)
    {
        if ( /*num_isotopic_atoms &&*/
            ( pINChI_Aux->nIsotopicOrigAtNosInCanonOrd = (AT_NUMB*) inchi_calloc( num_at_tg, sizeof( pINChI_Aux->nIsotopicOrigAtNosInCanonOrd[0] ) ) ) &&
             ( pINChI_Aux->nIsotopicOrigAtNosInCanonOrdInv = (AT_NUMB*) inchi_calloc( num_at_tg, sizeof( pINChI_Aux->nIsotopicOrigAtNosInCanonOrd[0] ) ) ) &&
             ( pINChI_Aux->nConstitEquIsotopicNumbers = (AT_NUMB*) inchi_calloc( num_at_tg, sizeof( pINChI_Aux->nConstitEquIsotopicNumbers[0] ) ) ))
        {
            ;
        }
        else if (num_isotopic_atoms)
        {
            goto out_of_RAM;
        }

        if ( /*num_isotopic_atoms && num_at > 1 &&*/
            ( pINChI_Aux->nConstitEquIsotopicTGroupNumbers = (AT_NUMB*) inchi_calloc( (long long)num_at / 2 + 1, sizeof( pINChI_Aux->nConstitEquIsotopicTGroupNumbers[0] ) ) )) /* djb-rwth: cast operator added */
        {
            ;
        }
        else if (num_isotopic_atoms && num_at > 1)
        {
            goto out_of_RAM;
        }
    }

    return pINChI_Aux;


out_of_RAM:

    if (pINChI_Aux)
    {
        Free_INChI_Aux( &pINChI_Aux );
        /*
        inchi_free(pINChI_Aux);
        */
    }

    return NULL;
}


/****************************************************************************
Note that orig_num, curr_num are allocd by caller as (n+1)-long lists
****************************************************************************/
void CompAtomData_GetNumMapping( COMP_ATOM_DATA *adata, int *orig_num, int *curr_num )
{
    int i;
    if (!orig_num || !curr_num)
    {
        return;
    }
    for (i = 0; i < adata->num_at; i++)
    {
        int orig = adata->at[i].orig_at_number;
        orig_num[i] = orig;    /* orig's are from 1 */
        curr_num[orig] = i;    /* curr's are from 0 */
    }
}


/****************************************************************************
 Allocate integer matrix [mxn]
****************************************************************************/
int imat_new( int m, int n, int ***a )
{
    int i;   
    if (m == 0 || n == 0)
    {
        return 0;
    }
    if (*a)
    {
        imat_free(m, *a);
        *a = NULL;
    }
    *a = (int **) inchi_calloc( m, sizeof( int * ) );
    if (NULL == *a)
    {
        return 1;
    }
    for (i = 0; i < m; i++)
    {
        ( *a )[i] = (int *) inchi_calloc( n, sizeof( int ) );
        if (NULL == ( *a )[i])
        {
            return 1;
        }
    }
    return 0;
}

/****************************************************************************
 Free integer matrix [mxn]
****************************************************************************/
void imat_free( int m, int **a )
{
    int i;
    if (NULL != a)
    {
        for (i = 0; i < m; i++)
        {
            if (NULL != a[i]) /* djb-rwth: unresolved issue -- revision required? -- false positive as this function just does the clean-up job */
            {
                inchi_free( a[i] );
            }
        }
        inchi_free( a );
        a = NULL;
    }

    return;
}


/*
Lt-wt subgraph
*/


/****************************************************************************
 Establish light-weight subgraph representing (part of) orig_inp_data
 ****************************************************************************/
subgraf *subgraf_new( ORIG_ATOM_DATA *orig_inp_data,
                      int nnodes,
                      int *nodes )
{
    int i, j, iat, nbr, jat, nj, degree, nat, err = 0;

    subgraf *sg = (subgraf *) inchi_calloc( 1, sizeof( subgraf ) );
    if (!sg)
    {
        return NULL;
    }

    nat = orig_inp_data->num_inp_atoms;

    /* orig2node is mapping of original at numbers --> subgraph node numbers */
    err = 1;
    if (!( sg->orig2node = (int *) inchi_calloc( (long long)nat + 1, sizeof( int ) ) )) /* djb-rwth: cast operator added */
    {
        goto exit_function;
    }
    if (!( sg->nodes = (int *) inchi_calloc( nnodes, sizeof( int ) ) ))
    {
        goto exit_function;
    }
    if (!( sg->degrees = (int *) inchi_calloc( nnodes, sizeof( int ) ) ))
    {
        goto exit_function;
    }
    /* NB:    input list of 'nodes' is assumed to be in 'original_atom_numbering domain' which starts from 1.
    Now it is mapped to current atom numbers which starts from 0/connections using at[j].orig_at_number    */
    sg->nnodes = 0;
    for (i = 0; i < nnodes; i++)
    {
        sg->nodes[sg->nnodes++] = nodes[i];
    }

    for (i = 0; i <= nat; i++)
    {
        sg->orig2node[i] = -1;
    }
    for (i = 0; i < nnodes; i++)
    {
        sg->orig2node[sg->nodes[i]] = i;
    }

    /* Create and fill subgraph adjacency matrix based on nodes/orig atom numbers 
       and connections stored in orig_inp_data */
    sg->adj = (subgraf_edge **) inchi_calloc( nnodes, sizeof( subgraf_edge * ) );
    if (!sg->adj)
    {
        goto exit_function;
    }

    for (i = 0; i < sg->nnodes; i++)
    {
        iat = nodes[i] - 1;    /* current atom  number for this node */
        degree = orig_inp_data->at[iat].valence;
        nj = -1;
        sg->adj[i] = (subgraf_edge *) inchi_calloc( degree, sizeof( subgraf_edge ) );
        if (!sg->adj[i])
        {
            goto exit_function;
        }
        for (j = 0; j < degree; j++)
        {
            jat = orig_inp_data->at[iat].neighbor[j]; /* for curr num jat, a (jat+1) would be an orig num */
            nbr = sg->orig2node[jat + 1];
            if (nbr < 0)
                continue;
            nj++;
            sg->adj[i][nj].nbr = nbr;
            sg->adj[i][nj].etype = orig_inp_data->at[iat].bond_type[j];
        }
        sg->degrees[i] = nj + 1;
    }
    err = 0;

    /* subgraf_debug_trace( sg ); */

exit_function:
    if (err)
    {
        subgraf_free( sg );
        return NULL; /* djb-rwth: avoiding reading from freed memory */
    }

    return sg;
}


/****************************************************************************/
void subgraf_free( subgraf *sg )
{
    int i;
    if (!sg)
    {
        return;
    }
    if (sg->nodes)
    {
        inchi_free( sg->nodes );
    }
    if (sg->degrees)
    {
        inchi_free( sg->degrees );
    }
    if (sg->orig2node)
    {
        inchi_free( sg->orig2node );
    }
    if (sg->adj)
    {
        for (i = 0; i < sg->nnodes; i++)
        {
            if (sg->adj[i]) /* djb-rwth: unresolved issue -- revision required? -- false positive as this function just does the clean-up job */
            {
                inchi_free( sg->adj[i] );
            }
        }
        inchi_free( sg->adj );
    }
    inchi_free( sg );
    sg = NULL;

    return;
}

/****************************************************************************/
void subgraf_debug_trace( subgraf *sg )
{
    int p, q;

    ITRACE_( "\n\n*********************************************************************\n* Subgraf:" );
    ITRACE_( "\n\tNodes: %-d ( ", sg->nnodes );
    for (p = 0; p < sg->nnodes; p++)
    {
        ITRACE_( "%-d ", sg->nodes[p] );
    }
    ITRACE_( ")\n\tAdj lists:\n" );
    for (p = 0; p < sg->nnodes; p++)
    {
        ITRACE_( "\tNode #%-d (orig# %-d) ::: Neighbors (node#, orig#) : ",
                 p, sg->nodes[p] );
        for (q = 0; q < sg->degrees[p]; q++)
        {
            int nbr = sg->adj[p][q].nbr; /* djb-rwth: ignoring LLVM warning: variable used */
            ITRACE_( "(%-d/%-d/%-d)  ", nbr, sg->nodes[nbr] );
        }
        ITRACE_( "\n" );
    }
    ITRACE_( "\n* End Subgraf\n*********************************************************************\n" );

    return;
}


/****************************************************************************/
subgraf_pathfinder * subgraf_pathfinder_new( subgraf *sg,
                                             ORIG_ATOM_DATA *orig_inp_data,
                                             int start,
                                             int end )
{
    subgraf_pathfinder *spf = NULL;

    spf = (subgraf_pathfinder *) inchi_calloc( 1, sizeof( subgraf_pathfinder ) );
    if (!spf)
    {
        goto exit_function;
    }

    spf->sg = sg;
    spf->start = start;
    spf->end = end;
    spf->nbonds = 0;
    spf->nseen = 0;

    spf->seen = (int *) inchi_calloc( spf->sg->nnodes, sizeof( int ) );
    if (!spf->seen)
    {
        inchi_free( spf );
        spf = NULL;
    }

exit_function:
    return spf;
}

/****************************************************************************/
void subgraf_pathfinder_free( subgraf_pathfinder *spf )
{
    if (!spf)
    {
        return;
    }
    if (spf->seen)
    {
        inchi_free( spf->seen );
    }
    inchi_free( spf );
    spf = NULL;
    return;
}


/****************************************************************************
 Find path(s) from subgraf node spf->start to spf->end 
 and fill bonds[nbonds] and atoms[natoms]
 Do not traverse through supplied forbidden edges (if not zero/NULL)
****************************************************************************/
void subgraf_pathfinder_run( subgraf_pathfinder *spf,
                             int nforbidden,		/* number of edges forbidden for traversal	*/
                             int *forbidden,		/* nodes of forbidden edges: [edge1node1,edge1node2, edge2node1, edge2node2, ... ] */
                             int *nbonds,
                             int **bonds,			/* collect subgraf bonds here */
                             int *natoms,
                             int *atoms				/* if not NULL, collect subgraf atoms here	*/
                             )
{
    int j, k, node, node0;
    int f, skip;

    if (spf->nseen < 1)
    {
        /*    Even at very beginning, push start node to seen and set nseen = 1
        and put end node into subgraf_pathfinder's end                            */
        return;
    }

    node0 = spf->seen[spf->nseen - 1];
    for (j = 0; j < spf->sg->degrees[node0]; j++)
    {
        node = spf->sg->adj[node0][j].nbr;
        if (is_in_the_ilist( spf->seen, node, spf->nseen ))
        {
            continue;
        }
        if (nforbidden && forbidden)
        {	
            skip = 0;
            for (f = 0; f < nforbidden; f++)
            {
                if (bIsSameBond(node0, node, forbidden[2 * f], forbidden[2 * f + 1]) )
                {
                    skip = 1;
                    break;
                }
            }
            if (skip)
            {
                continue;
            }
        }
        if (node == spf->end)
        {
            spf->seen[spf->nseen++] = node;

            ITRACE_( "\n\tFound path (in orig atom numbers):\t" );
            for (k = 0; k < spf->nseen; k++)
            {
                int orig_atnum = spf->sg->nodes[spf->seen[k]];
                ITRACE_( "%-d ", orig_atnum);
                if (atoms && !is_in_the_ilist(atoms, orig_atnum, *natoms))
                {
                    atoms[(*natoms)++] = orig_atnum;
                }
            }
            ITRACE_( "\t( In node nums: " );
            for (k = 1; k < spf->nseen; k++)
            {
                int at1 = spf->seen[k - 1];
                int at2 = spf->seen[k];
                add_bond_if_unseen( spf, at1, at2, nbonds, bonds );

                ITRACE_( "%-d ", spf->seen[k] );
            }
            ITRACE_( ")" );

            spf->seen[spf->nseen - 1] = 0;
            spf->nseen--;    /* pop_back        */
            break;
        }
    }
    for (j = 0; j < spf->sg->degrees[node0]; j++)
    {
        node = spf->sg->adj[node0][j].nbr;
        if (node == spf->end || is_in_the_ilist( spf->seen, node, spf->nseen ))
        {
            continue;
        }
        if (nforbidden && forbidden)
        {
            skip = 0;
            for (f = 0; f < nforbidden; f++)
            {
                if (bIsSameBond(node0, node, forbidden[2 * f], forbidden[2 * f + 1]))
                {
                    skip = 1;
                    break;
                }
            }
            if (skip)
            {
                continue;
            }
        }
        spf->seen[spf->nseen++] = node;
        subgraf_pathfinder_run( spf, 0, NULL, nbonds, bonds, natoms, atoms );
        spf->seen[spf->nseen - 1] = 0;
        spf->nseen--;
    }

    return;
}


/****************************************************************************/
void add_bond_if_unseen( subgraf_pathfinder *spf,
                         int node0,
                         int node,
                         int *nbonds,
                         int **bonds )
{
    int seen, p, at1, at2;

    at1 = spf->sg->nodes[node0];
    at2 = spf->sg->nodes[node];
#if 0
    if (at1 > at2)
    {
        int tmp = at1;
        at1 = at2;
        at2 = tmp;
    }
#endif
    seen = 0;
    for (p = 0; p < *nbonds; p++)
    {
        /*if (bonds[p][0] == at1 && bonds[p][1] == at2)*/
        if (bIsSameBond(at1, at2, bonds[p][0], bonds[p][1]))
        {
            seen = 1;
            break;
        }
    }
    if (!seen)
    {
        bonds[*nbonds][0] = at1;
        bonds[*nbonds][1] = at2;
        ( *nbonds )++;
    }

    return;
}


/****************************************************************************
 At the first call, push start node to spf->start and set spf->nseen = 0
****************************************************************************/
int subgraf_pathfinder_collect_all(subgraf_pathfinder *spf,
                                   int nforbidden,		/* number of edges forbidden for traversal	*/
                                   int *forbidden,		/* nodes of forbidden edges: [edge1node1,edge1node2, edge2node1, edge2node2, ... ] */
                                   int *atnums          /* 1-based origs# */
                                    )
{
    int j, f, node, next_node, skip;

    node = spf->start;
    spf->seen[spf->nseen] = node;
    atnums[spf->nseen] = spf->sg->nodes[node];
    spf->nseen++;

    for (j = 0; j < spf->sg->degrees[node]; j++)
    {
        next_node = spf->sg->adj[node][j].nbr;
        if (is_in_the_ilist(spf->seen, next_node, spf->nseen))
        {
            continue;
        }
        if (nforbidden && forbidden)
        {
            skip = 0;
            for (f = 0; f < nforbidden; f++)
            {
                if (bIsSameBond(node, next_node, forbidden[2 * f], forbidden[2 * f + 1]))
                {
                    skip = 1;
                    break;
                }
            }
            if (skip)
            {
                continue;
            }
        }
        spf->start = next_node;
        subgraf_pathfinder_collect_all(spf, nforbidden, forbidden, atnums);
    }

    return spf->nseen;
}


#if ( FIX_ADJ_RAD == 1 )


/****************************************************************************/
int FixNextRadicals( int cur_at, inp_ATOM *at );
int FixNextRadicals( int cur_at, inp_ATOM *at )
{
    int j, neigh, num_found = 0;

    for (j = 0; j < at[cur_at].valence; j++)
    {
        neigh = at[cur_at].neighbor[j];
        if (at[neigh].radical == RADICAL_DOUBLET)
        {
            at[neigh].radical = 0;
            num_found++;
            num_found += FixNextRadicals( neigh, at );
        }
    }

    return num_found;
}


/****************************************************************************/
int FixAdjacentRadicals( int num_inp_atoms, inp_ATOM *at )
{
    int i, j;
    char *bVisited = NULL;
    int  nNumFound = 0, neigh, cur_found;

    for (i = 0; i < num_inp_atoms; i++)
    {
        if (at[i].radical == RADICAL_DOUBLET)
        {
            cur_found = 1;
            for (j = 0; j < at[i].valence; j++)
            {
                neigh = at[i].neighbor[j];
                if (at[neigh].radical == RADICAL_DOUBLET)
                {
                    cur_found++;
                }
            }
            if (cur_found >= 3)
            {
                nNumFound++;
                at[i].radical = 0;
                nNumFound += FixNextRadicals( i, at );
            }
        }
    }

    return nNumFound;
}

#endif

#ifdef COMPILE_ANSI_ONLY

#ifndef TARGET_API_LIB
/*
#include <stdio.h>
#include "inpdef.h"
*/

/****************************************************************************/
void PrintFileName( const char *fmt,
                    FILE *out_file,
                    /* INCHI_IOSTREAM *out_file,  */
                    const char *szFname )
{
    inchi_print_nodisplay( out_file, fmt, szFname );
}
#endif

#endif

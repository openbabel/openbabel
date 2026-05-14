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

#include "mode.h"
#include "ichicant.h"
#include "ichicomn.h"
#include "ichister.h"

#include "bcf_s.h"

typedef struct tagStereoBondNeighbor
{
    /*  *n = sort key */
    AT_RANK       nRank;            /* 1 opposite atom rank; equal ranks mean constit. equivalence  */
    AT_RANK       nNeighRank1;      /* rank of the neighbor in the direction to the opposite atom   */
    AT_RANK       nNeighRank2;      /* rank of the opposite atom neighbor in the direction          */
                                    /* to the current atom                                          */
    AT_RANK       num;              /* number of same type bonds to constitutionally                */
                                    /* equivalent neighbors                                         */
    AT_RANK       num_any_parity;   /*  at least one atom has parity in 1..4 range                  */
    AT_RANK       num_defined_parity; /*  number of neighbors with defined parity <= num            */
    /* AT_RANK       num_undef_parity; */
    /* AT_RANK       num_unkn_parity;  */
    AT_RANK       what2do;
    U_CHAR        cumulene_len;     /*  high nimble bits: (cumulene length - 1)                     */
    U_CHAR        bond_type;        /*  *2 all same, not a real bond type                           */
} STEREO_BOND_NEIGH;



/* Local prototypes */

int SetHalfStereoBondIllDefPariy( sp_ATOM *at,
                                  int     jn, /* atom number*/
                                  int     k1 /* stereo bond number*/,
                                  int     new_parity );

int RemoveHalfStereoBond( sp_ATOM *at,
                          int     jn, /* atom number*/
                          int     k1 /* stereo bond number*/ );

int SetKnownStereoBondParities( CANON_GLOBALS *pCG,
                                sp_ATOM       *at,
                                int           num_atoms,
                                const AT_RANK *nCanonRank,
                                const AT_RANK *nRank,
                                const AT_RANK *nAtomNumber );

int MarkKnownEqualStereoBondParities( sp_ATOM       *at,
                                      int           num_atoms,
                                      const AT_RANK *nRank,
                                      const AT_RANK *nAtomNumber );

int GetNextNeighborAndRank( sp_ATOM        *at,
                            AT_RANK        cur,
                            AT_RANK        prev,
                            AT_RANK        *n,
                            AT_RANK        *cr,
                            const AT_RANK *nCanonRank );

int GetAndCheckNextNeighbors( sp_ATOM       *at,
                              AT_RANK       cur1,
                              AT_RANK       prev1,
                              AT_RANK       cur2,
                              AT_RANK       prev2,
                              AT_RANK       *n1,
                              AT_RANK       *n2,
                              AT_RANK       *nVisited1,
                              AT_RANK       *nVisited2,
                              const AT_RANK *nRank,
                              const AT_RANK *nCanonRank );

AT_RANK PathsHaveIdenticalKnownParities( sp_ATOM       *at,
                                         AT_RANK       prev1,
                                         AT_RANK       cur1,
                                         AT_RANK       prev2,
                                         AT_RANK       cur2,
                                         AT_RANK       *nVisited1,
                                         AT_RANK       *nVisited2,
                                         const AT_RANK *nRank,
                                         const AT_RANK *nCanonRank,
                                         AT_RANK        nLength );

int RemoveKnownNonStereoBondParities( sp_ATOM       *at,
                                      int           num_atoms,
                                      const AT_RANK *nCanonRank,
                                      const AT_RANK *nRank,
                                      CANON_STAT    *pCS );

int SetKnownStereoCenterParities( CANON_GLOBALS  *pCG,
                                  sp_ATOM        *at,
                                  int            num_atoms,
                                  const AT_RANK  *nCanonRank,
                                  const AT_RANK  *nRank,
                                  const AT_RANK  *nAtomNumber );
int RemoveKnownNonStereoCenterParities( CANON_GLOBALS *pCG,
                                        sp_ATOM       *at,
                                        int           num_atoms,
                                        const AT_RANK *nCanonRank,
                                        const AT_RANK *nRank,
                                        CANON_STAT    *pCS );

int MarkKnownEqualStereoCenterParities( sp_ATOM       *at,
                                        int           num_atoms,
                                        const AT_RANK *nRank,
                                        const AT_RANK *nAtomNumber );




/****************************************************************************
  Depth First Search for an atom with parity
****************************************************************************/
int find_atoms_with_parity( sp_ATOM *at,
                            S_CHAR  *visited,
                            int     from_atom,
                            int     cur_atom )
{
    int i, next_atom;

    if (visited[cur_atom])
    {
        return 0;
    }
    if (at[cur_atom].parity)
    {
        return 1;
    }

    visited[cur_atom] = 1;

    for (i = 0; i < at[cur_atom].valence; i++)
    {
        next_atom = at[cur_atom].neighbor[i];

        if (next_atom != from_atom &&
             find_atoms_with_parity( at, visited, cur_atom, next_atom ))
        {
            return 1;
        }
    }

    return 0;
}


/****************************************************************************/
int SetHalfStereoBondIllDefPariy( sp_ATOM *at,
                                  int     jn, /* atom number*/
                                  int     k1 /* stereo bond number*/,
                                  int     new_parity )
{
    int parity;
    if (k1 < MAX_NUM_STEREO_BOND_NEIGH && at[jn].stereo_bond_neighbor[k1])
    {
        parity = at[jn].stereo_bond_parity[k1] ^ PARITY_VAL( at[jn].stereo_bond_parity[k1] );
        at[jn].stereo_bond_parity[k1] = parity | PARITY_VAL( new_parity );
        at[jn].parity = PARITY_VAL( new_parity );
        return 1;  /*  success */
    }

    return 0; /*  failed             */
}


/****************************************************************************/
int RemoveHalfStereoBond( sp_ATOM *at,
                          int     jn, /* atom number*/
                          int     k1 /* stereo bond number*/ )
{
    int k2;
    if (k1 < MAX_NUM_STEREO_BOND_NEIGH && at[jn].stereo_bond_neighbor[k1])
    {
        for (k2 = k1; k2 < MAX_NUM_STEREO_BOND_NEIGH - 1; k2++) /* djb-rwth: loop condition corrected (buffer error) */
        {
            at[jn].stereo_bond_neighbor[k2] = at[jn].stereo_bond_neighbor[k2 + 1];
            at[jn].stereo_bond_ord[k2] = at[jn].stereo_bond_ord[k2 + 1];
            at[jn].stereo_bond_z_prod[k2] = at[jn].stereo_bond_z_prod[k2 + 1];
            at[jn].stereo_bond_parity[k2] = at[jn].stereo_bond_parity[k2 + 1];
        }
        at[jn].stereo_bond_neighbor[k2] = 0;
        at[jn].stereo_bond_ord[k2] = 0;
        at[jn].stereo_bond_z_prod[k2] = 0;
        at[jn].stereo_bond_parity[k2] = 0;

        if (!at[jn].stereo_bond_neighbor[0])
        {   /*  curled braces added 6-6-2002 */
            at[jn].parity = 0;
            at[jn].stereo_atom_parity = 0;
            at[jn].final_parity = 0;
            /* at[jn].bHasStereoOrEquToStereo = 0; */
        }
        return 1; /*  success            */
    }

    return 0; /*  failed             */
}


/****************************************************************************/
int SetOneStereoBondIllDefParity( sp_ATOM *at,
                                  int     jc,       /* atom number              */
                                  int     k,        /* stereo bond ord. number  */
                                  int     new_parity )
{
    int k1, ret = 0, kn, jn = (int) at[jc].stereo_bond_neighbor[k] - 1;

    /*  opposite end */
    for (k1 = ret = 0;
            k1 < MAX_NUM_STEREO_BOND_NEIGH && ( kn = at[jn].stereo_bond_neighbor[k1] );
              k1++) /* djb-rwth: removing redundant code */
    {
        if (kn - 1 == jc)
        {
            ret = SetHalfStereoBondIllDefPariy( at, jn, /* atom number*/ k1 /* stereo bond number*/, new_parity );
            break;
        }
    }

    if (ret)
    {
        ret = SetHalfStereoBondIllDefPariy( at, jc, k, new_parity );
    }

    return ret;
}


/****************************************************************************/
int RemoveOneStereoBond( sp_ATOM *at,
                         int     jc,    /* atom number          */
                         int     k      /* stereo bond number   */
)
{
    int k1, ret = 0, kn, jn = (int) at[jc].stereo_bond_neighbor[k] - 1;

    /*  opposite end */
    for (k1 = ret = 0;
            k1 < MAX_NUM_STEREO_BOND_NEIGH && ( kn = at[jn].stereo_bond_neighbor[k1] );
              k1++) /* djb-rwth: removing redundant code */
    {
        if (kn - 1 == jc)
        {
            ret = RemoveHalfStereoBond( at, jn, k1 );
            break;
        }
    }

    if (ret)
    {
        ret = RemoveHalfStereoBond( at, jc, k );
    }

    return ret;
}


/****************************************************************************/
int RemoveOneStereoCenter( sp_ATOM *at,
                           int     jc /* atom number*/ )
{
    if (at[jc].parity)
    {
        at[jc].parity = 0; /*  remove parity */
        at[jc].stereo_atom_parity = 0;
        at[jc].final_parity = 0;
        /*  at[jc].bHasStereoOrEquToStereo = 0; */
        return 1;
    }

    return 0; /*  failed: not a stereo center */
}


/****************************************************************************
  Remove stereo parity from centers having constitutionally equivalent
  cut-vertex neighbors whose attachments do not have stereogenic elements.
  Currently checks ALL constitutionally equivalent neighbors.
  To optimize, check only one.
****************************************************************************/
int UnmarkNonStereo( CANON_GLOBALS *pCG,
                     sp_ATOM       *at,
                     int           num_atoms,
                     const AT_RANK *nRank,
                     const AT_RANK *nAtomNumber,
                     int           bIsotopic )
{
    int i, i1, i2, j, k, k1, k2, kn /* neigh*/, val, ic/* center*/, jc, num_implicit_H;
    int num_neighbors_with_parity, num_no_parity_atoms, num_removed_parities = -1, num_removed_parities0;
    AT_RANK nNeighborNumber[MAX_NUM_STEREO_ATOM_NEIGH];
    AT_RANK nPrevAtomRank, nPrevNeighRank;
#ifdef FIX_OLEAN_SPIRO_CHIRALITY_DETECTION_BUG
    int num_in_same_ring_system = 1, nRingSystem, num_with_eq_neigh_in_same_ring_system = 0; /* djb-rwth: although unlikely to ever occur, uninitialised num_in_same_ring_system variable can lead to garbage value, including 0 which leads to various errors and inconsistency with 1.06 outputs -- function rewriting and discussion required */
#endif


    S_CHAR *visited = (S_CHAR *) inchi_malloc( num_atoms * sizeof( visited[0] ) );

    if (!visited)
    {
        goto exit_function;
    }

    num_removed_parities = 0;
    num_no_parity_atoms = 0;

    do
    {
        num_removed_parities0 = num_removed_parities;

        for (i = i1 = 0, nPrevAtomRank = 0; i <= num_atoms; i++)
        {
            /*  bounds violation check (i!=num_atoms) added 6-21-2002 */
            if (i == num_atoms || nPrevAtomRank != nRank[j = nAtomNumber[i]]
                 /* at[j].parity && 1 < at[j].valence && at[j].valence < MAX_NUM_STEREO_ATOM_NEIGH*/)
            {
                /*  end of constitutionally equivalent atoms sequence */
                /* nPrevRank = nRank[j]; */
                i2 = i;
                if (i2 - i1 > num_no_parity_atoms /*&& at[jc = nAtomNumber[i1]].parity*/)
                {
                    /*  at[nAtomNumber[i1]]..at[nAtomNumber[i2-1]] are constitutionally equivalent and some of them have parity */
                    jc = nAtomNumber[i1];
                    num_no_parity_atoms = 0;
                    val = at[jc].valence; /*  all equivalent atoms have equal valences, etc. (except parities) */
                    num_implicit_H = at[jc].endpoint ? 0 : at[jc].num_H;
                    /*  Only atoms with valence <= MAX_NUM_STEREO_ATOM_NEIGH may have parity. However, check: */
                    if (val + num_implicit_H > MAX_NUM_STEREO_ATOM_NEIGH)
                    {
                        continue;  /*  program error ??? */ /*   <BRKPT> */
                    }
                    for (k = 0; k < val; k++)
                    {
                        nNeighborNumber[k] = k; /*  initialize an array of indexes for sorting */
                    }
                    /*  check parities */
                    for (ic = i1; ic < i2; ic++)
                    {
                        jc = nAtomNumber[ic];
                        /*  sort neighbors according to their canon. equivalence ranks */
                        pCG->m_pNeighborsForSort = at[jc].neighbor;
                        pCG->m_pn_RankForSort = nRank;
                        insertions_sort( pCG, nNeighborNumber, val, sizeof( nNeighborNumber[0] ), CompNeighborsAT_NUMBER );
                        num_neighbors_with_parity = -1; /*  non-zero */
                        for (k = k1 = 0, nPrevNeighRank = 0; k <= val; k++)
                        {
                            if (k == val || nPrevNeighRank != nRank[at[jc].neighbor[nNeighborNumber[k]]])
                            {
                                k2 = k;
                                if (k2 - k1 > 1)
                                {
                                    /*  found 2 or more constitutionally equivalent neighbors */
                                    /*  Check if they have only non-stereogenic neighbors */
#ifdef FIX_OLEAN_SPIRO_CHIRALITY_DETECTION_BUG
                                    num_in_same_ring_system = nRingSystem = 0;
                                    for (kn = k1; kn < k2; kn++) /* djb-rwth: removing redundant code */
                                    {
                                        int nCurNeighRingSystem = at[(int) at[jc].neighbor[nNeighborNumber[kn]]].nRingSystem;
                                        if (!nRingSystem)
                                        {
                                            nRingSystem = nCurNeighRingSystem;
                                        }
                                        else
                                        {
                                            num_in_same_ring_system += ( nRingSystem == nCurNeighRingSystem );
                                        }
                                    }
#endif

                                    for (kn = k1, num_neighbors_with_parity = 0; kn < k2; kn++)
                                    {
                                        memset( visited, 0, num_atoms * sizeof( visited[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
                                        visited[jc] = 1; /*  starting point; the only atom with parity */
                                        num_neighbors_with_parity +=
                                            find_atoms_with_parity( at, visited, jc, (int) at[jc].neighbor[nNeighborNumber[kn]] );
                                    }
                                }
                                /* if ( !num_neighbors_with_parity ) */
#ifdef FIX_OLEAN_SPIRO_CHIRALITY_DETECTION_BUG
                                if (!num_neighbors_with_parity && !num_in_same_ring_system)
#else
                                if (!num_neighbors_with_parity)
#endif
                                {
                                    break; /*  at[jc] cannot have defined parity */
                                }
                                if (k + 1 < val)
                                {
                                    k1 = k; /*  at least 2 more neighbors left */
                                    nPrevNeighRank = nRank[at[jc].neighbor[nNeighborNumber[k]]];
                                }
                                else
                                {
                                    break;
                                }
                            }
                        }
                        if (num_implicit_H > 1)
                        {
                            if ((bIsotopic && ( at[jc].num_iso_H[0] > 1 ||
                                at[jc].num_iso_H[1] > 1 ||
                                at[jc].num_iso_H[2] > 1 )) ||
                                  num_implicit_H > NUM_H_ISOTOPES ||
                                  !bIsotopic) /* djb-rwth: addressing LLVM warning */
                            {
                                num_neighbors_with_parity = 0;
                            }
                        }
                        /*  increment if: */
                        /*  (a) constitutionally equivalent neighbors do exist, and */
                        /*  (b) all constitutionally equivalent neighbors do not have parity, and */
                        /*  (c) all constitutionally equivalent neighbors are not connected to atoms with parity */
                        num_no_parity_atoms += !num_neighbors_with_parity;
#ifdef FIX_OLEAN_SPIRO_CHIRALITY_DETECTION_BUG
                        num_with_eq_neigh_in_same_ring_system += ( num_in_same_ring_system != 0 ); /* djb-rwth: initialisation of num_in_same_ring_system is required to avoid garbage value */
#endif
                    }
#ifdef FIX_OLEAN_SPIRO_CHIRALITY_DETECTION_BUG
                    if (num_no_parity_atoms == i2 - i1 && num_with_eq_neigh_in_same_ring_system != i2 - i1)
#else
                    if (num_no_parity_atoms == i2 - i1)
#endif

                    {
                        /*  all atoms at[nAtomNumber[i1]]..at[nAtomNumber[i2-1]] cannot be */
                        /*  stereo centers or belong to stereo bonds */
                        for (ic = i1; ic < i2; ic++)
                        {
                            int jn;
                            jc = nAtomNumber[ic];
                            at[jc].parity = 0; /*  remove parity */
                            at[jc].stereo_atom_parity = 0;
                            at[jc].final_parity = 0;
                            at[jc].bHasStereoOrEquToStereo = 0;
                            /*  remove stereo bonds */
                            for (k = 0; k < MAX_NUM_STEREO_BOND_NEIGH && ( jn = at[jc].stereo_bond_neighbor[k] ); k++)
                            {
                                jn--; /*  stereo bond neighbor */
                                /*  opposite end */
                                for (k1 = 0; k1 < MAX_NUM_STEREO_BOND_NEIGH && ( kn = at[jn].stereo_bond_neighbor[k1] ); k1++)
                                {
                                    if (kn - 1 == jc)
                                    {
                                        RemoveHalfStereoBond( at, jn, k1 );
                                        break;
                                    }
                                }
                                /*  at at[jc] stereo bond end; since references to all at[jc] */
                                /*  stereo bond neighbors are to be removed, do not shift them */
                                at[jc].stereo_bond_neighbor[k] = 0;
                                at[jc].stereo_bond_ord[k] = 0;
                                at[jc].stereo_bond_z_prod[k] = 0;
                                at[jc].stereo_bond_parity[k] = 0;
                            }
                        }
                        num_removed_parities += num_no_parity_atoms;
                    }
                }
                if (i < num_atoms)
                {
                    nPrevAtomRank = nRank[j];
                    i1 = i;
                }
                num_no_parity_atoms = 0;
            }
            num_no_parity_atoms += ( i < num_atoms && !at[j].parity );
        }
    }
    while (num_removed_parities != num_removed_parities0);

exit_function:
    if (visited)
    {
        inchi_free( visited );
    }

    return num_removed_parities;
}


/****************************************************************************
  Add stereo descriptor(s) for atom #i
****************************************************************************/
int FillSingleStereoDescriptors( CANON_GLOBALS  *pCG,
                                 sp_ATOM        *at,
                                 int            i,
                                 int            num_trans,
                                 const AT_RANK  *nRank,
                                 AT_STEREO_CARB *LinearCTStereoCarb,
                                 int            *nStereoCarbLen,
                                 int            nMaxStereoCarbLen,
                                 AT_STEREO_DBLE *LinearCTStereoDble,
                                 int            *nStereoDbleLen,
                                 int            nMaxStereoDbleLen,
                                 int            bAllene )
{

    if (!LinearCTStereoDble && !LinearCTStereoCarb)
    {
        return 0; /*  return immediately if no stereo have been requested */
    }

    /***************************************************
    add stereo centers and stereo bonds to the CT
    ***************************************************/
    if (at[i].parity || at[i].stereo_bond_neighbor[0])
    {
        AT_RANK r_neigh, rank = nRank[i];
        AT_NUMB nNeighborNumber2[MAXVAL];
        unsigned parity;
        int      k;
        int num_allene = 0;

        if (ATOM_PARITY_WELL_DEF( at[i].parity ) && num_trans < 0)
        {
            /*  number of neighbors transpositions to the sorted order is unknown. Find it. */
            /*  If parity is not well-defined then doing this is a waste of time */
            int num_neigh = at[i].valence;
            for (k = 0; k < num_neigh; k++)
            {
                nNeighborNumber2[k] = k;
            }

            pCG->m_pNeighborsForSort = at[i].neighbor;
            pCG->m_pn_RankForSort = nRank;
            num_trans = insertions_sort( pCG, nNeighborNumber2, num_neigh, sizeof( nNeighborNumber2[0] ), CompNeighborsAT_NUMBER );

#ifndef CT_NEIGH_INCREASE
            num_trans += ( ( num_neigh*( num_neigh - 1 ) ) / 2 ) % 2;  /*  get correct parity for ascending order */
#endif
        }

        /*  stereo bonds */
        if (LinearCTStereoDble && at[i].stereo_bond_neighbor[0])
        {

            /* HalfStereoBondParity( sp_ATOM *at, int at_no1, int i_sb_neigh, AT_RANK *nRank ) */
            AT_NUMB nStereoNeighNumber[MAX_NUM_STEREO_BONDS], nStereoNeigh[MAX_NUM_STEREO_BONDS], n;
            int       num_stereo, stereo_neigh, stereo_neigh_ord, stereo_bond_parity;
            for (num_stereo = 0;
                      num_stereo < MAX_NUM_STEREO_BONDS &&
                      ( n = at[i].stereo_bond_neighbor[num_stereo] ); num_stereo++)
            {
                nStereoNeighNumber[num_stereo] = num_stereo;
                nStereoNeigh[num_stereo] = n - 1;
                num_allene += IS_ALLENE_CHAIN( at[i].stereo_bond_parity[num_stereo] );
            }
            if ((bAllene > 0 && !num_allene) || (bAllene == 0 && num_allene)) /* djb-rwth: addressing LLVM warning */
            {
                return 0;
            }

            /*  sort stereo bonds according to the ranks of the neighbors */
            pCG->m_pNeighborsForSort = nStereoNeigh;
            pCG->m_pn_RankForSort = nRank;
            insertions_sort( pCG, nStereoNeighNumber, num_stereo, sizeof( nStereoNeighNumber[0] ), CompNeighborsAT_NUMBER );

            /*  process stereo bonds one by one */
            for (k = 0; k < num_stereo; k++)
            {
                stereo_neigh = nStereoNeigh[stereo_neigh_ord = (int) nStereoNeighNumber[k]];

                if (( r_neigh = (AT_NUMB) nRank[stereo_neigh] ) CT_NEIGH_SMALLER_THAN rank)
                {
                    /* accept only neighbors that have smaller ranks */
                    stereo_bond_parity = PARITY_VAL( at[i].stereo_bond_parity[stereo_neigh_ord] );
                    if (stereo_bond_parity == AB_PARITY_NONE)
                    {
                        continue;
                    }

                    /* stereo_neigh      = at[i].stereo_bond_neighbor[nStereoNeighNumber[k]]-1; */
                    if (ATOM_PARITY_KNOWN( stereo_bond_parity ))
                    {
                        parity = stereo_bond_parity;
                    }
                    else if (ATOM_PARITY_WELL_DEF( at[i].parity ) &&
                              ATOM_PARITY_WELL_DEF( at[stereo_neigh].parity ) &&
                              MIN_DOT_PROD <= abs( at[i].stereo_bond_z_prod[stereo_neigh_ord] ))
                    {
                        /*  bond parity can be calculated */
                        int half_parity1, half_parity2, j, nn, stereo_neigh_ord2;
                        stereo_neigh_ord2 = -1;
                        for (j = 0; j < MAX_NUM_STEREO_BONDS &&
                            ( nn = (int) at[stereo_neigh].stereo_bond_neighbor[j] );
                                       j++)
                        {
                            if (i + 1 == nn)
                            {
                                /* found the opposite end of the stereo bond */
                                stereo_neigh_ord2 = j;
                                break;
                            }
                        }
                        if (stereo_neigh_ord2 >= 0)
                        {
                            half_parity1 = HalfStereoBondParity( at, i, stereo_neigh_ord, nRank );
                            half_parity2 = HalfStereoBondParity( at, stereo_neigh, stereo_neigh_ord2, nRank );
                            if (ATOM_PARITY_WELL_DEF( half_parity1 ) &&
                                 ATOM_PARITY_WELL_DEF( half_parity2 ))
                            {
                                parity = 2 - ( half_parity1 + half_parity2
                                         + ( at[i].stereo_bond_z_prod[stereo_neigh_ord] < 0 ) ) % 2;
                            }
                            else
                            {
                                return CT_STEREOBOND_ERROR;  /*   <BRKPT> */
                            }
                        }
                        else
                        {
                            return CT_STEREOBOND_ERROR;  /*   <BRKPT> */
                        }
                    }
                    else
                    {
                        /*  parity cannot be calculated: not enough info or 'unknown' */
                        if (AB_PARITY_NONE == ( parity = inchi_max( at[i].parity, at[stereo_neigh].parity ) ))
                        {
                            continue;
                        }
                        if (ATOM_PARITY_WELL_DEF( parity ))
                        {
                            parity = AB_PARITY_UNDF; /*  should not happen */
                        }
                    }
                    if (CHECK_OVERFLOW( *nStereoDbleLen, nMaxStereoDbleLen ))
                        return CT_OVERFLOW;  /*   <BRKPT> */
                    /*  first stereo bond atom */
                    LinearCTStereoDble[*nStereoDbleLen].at_num1 = rank;
                    /*  second stereo bond atom (opposite end) */
                    LinearCTStereoDble[*nStereoDbleLen].at_num2 = r_neigh;
                    /*  bond parity */
                    LinearCTStereoDble[*nStereoDbleLen].parity = parity;
                    ( *nStereoDbleLen )++;
                }
            }
        }

        /*  stereo carbon */
        if (bAllene > 0)
        {
            return 0;
        }

        if (LinearCTStereoCarb && !at[i].stereo_bond_neighbor[0])
        {
            if (CHECK_OVERFLOW( *nStereoCarbLen, nMaxStereoCarbLen ))
                return CT_OVERFLOW;  /*   <BRKPT> */
            /*  stereo atom rank */
            LinearCTStereoCarb[*nStereoCarbLen].at_num = rank;
            /*  stereo atom parity */
            parity = ATOM_PARITY_WELL_DEF( at[i].parity ) ? ( 2 - ( at[i].parity + num_trans ) % 2 ) : at[i].parity;
            LinearCTStereoCarb[*nStereoCarbLen].parity = parity;
            ( *nStereoCarbLen )++;
        }
    }

    return 0;
}


/****************************************************************************/
void SwitchAtomStereoAndIsotopicStereo( sp_ATOM *at,
                                        int     num_atoms,
                                        int     *bSwitched )
{
    int i;
    /*  switch atom stereo data */
    for (i = 0; i < num_atoms; i++)
    {
        inchi_swap( (char*) &at[i].parity, (char*) &at[i].parity2, sizeof( at[i].parity ) );
        inchi_swap( (char*) &at[i].final_parity, (char*) &at[i].final_parity2, sizeof( at[i].final_parity ) );
        inchi_swap( (char*) &at[i].stereo_atom_parity, (char*) &at[i].stereo_atom_parity2, sizeof( at[i].stereo_atom_parity ) );
        inchi_swap( (char*) &at[i].bHasStereoOrEquToStereo, (char*) &at[i].bHasStereoOrEquToStereo2, sizeof( at[i].bHasStereoOrEquToStereo ) );

        inchi_swap( (char*) at[i].stereo_bond_neighbor, (char*) at[i].stereo_bond_neighbor2, sizeof( at[i].stereo_bond_neighbor ) );
        inchi_swap( (char*) at[i].stereo_bond_ord, (char*) at[i].stereo_bond_ord2, sizeof( at[i].stereo_bond_ord ) );
        inchi_swap( (char*) at[i].stereo_bond_z_prod, (char*) at[i].stereo_bond_z_prod2, sizeof( at[i].stereo_bond_z_prod ) );
        inchi_swap( (char*) at[i].stereo_bond_parity, (char*) at[i].stereo_bond_parity2, sizeof( at[i].stereo_bond_parity ) );
    }

    *bSwitched = !*bSwitched;
}


/****************************************************************************/
void SetCtToIsotopicStereo( CANON_STAT *pCS,
                            CANON_STAT *pCS2 )
{
    pCS->LinearCTStereoDble = pCS2->LinearCTIsotopicStereoDble; /*  enable stereo */
    pCS->LinearCTStereoCarb = pCS2->LinearCTIsotopicStereoCarb;

    pCS->LinearCTStereoDbleInv = pCS2->LinearCTIsotopicStereoDbleInv; /*  enable inv. stereo */
    pCS->LinearCTStereoCarbInv = pCS2->LinearCTIsotopicStereoCarbInv;
    pCS->nMaxLenLinearCTStereoDble = pCS2->nMaxLenLinearCTIsotopicStereoDble;
    pCS->nMaxLenLinearCTStereoCarb = pCS2->nMaxLenLinearCTIsotopicStereoCarb;

    pCS->nLenLinearCTStereoDble = pCS2->nLenLinearCTIsotopicStereoDble;
    pCS->nLenLinearCTStereoCarb = pCS2->nLenLinearCTIsotopicStereoCarb;
}


/****************************************************************************/
void SetCtToNonIsotopicStereo( CANON_STAT *pCS,
                               CANON_STAT *pCS2 )
{
    pCS->LinearCTStereoDble = pCS2->LinearCTStereoDble; /*  enable stereo */
    pCS->LinearCTStereoCarb = pCS2->LinearCTStereoCarb;

    pCS->LinearCTStereoDbleInv = pCS2->LinearCTStereoDbleInv; /*  enable inv. stereo */
    pCS->LinearCTStereoCarbInv = pCS2->LinearCTStereoCarbInv;
    pCS->nMaxLenLinearCTStereoDble = pCS2->nMaxLenLinearCTStereoDble;
    pCS->nMaxLenLinearCTStereoCarb = pCS2->nMaxLenLinearCTStereoCarb;

    pCS->nLenLinearCTStereoDble = pCS2->nLenLinearCTStereoDble;
    pCS->nLenLinearCTStereoCarb = pCS2->nLenLinearCTStereoCarb;

    pCS->nLenLinearCTIsotopicStereoDble = pCS2->nLenLinearCTIsotopicStereoDble;
    pCS->nLenLinearCTIsotopicStereoCarb = pCS2->nLenLinearCTIsotopicStereoCarb;
}


/****************************************************************************/
int FillAllStereoDescriptors( CANON_GLOBALS *pCG,
                              sp_ATOM       *at,
                              int           num_atoms,
                              const         AT_RANK *nCanonRank,
                              const         AT_RANK *nAtomNumberCanon,
                              CANON_STAT    *pCS )
{
    int ret = 0, i;
    /*  initialize zero lengths */
    pCS->nLenLinearCTStereoCarb = 0;
    pCS->nLenLinearCTStereoDble = 0;

    /*  fill atom by atom */
    for (i = 0; !ret && i < num_atoms; i++)
    {
        ret = FillSingleStereoDescriptors( pCG, at, (int) nAtomNumberCanon[i], -1, nCanonRank
                          , pCS->LinearCTStereoCarb, &pCS->nLenLinearCTStereoCarb, pCS->nMaxLenLinearCTStereoCarb
                          , pCS->LinearCTStereoDble, &pCS->nLenLinearCTStereoDble, pCS->nMaxLenLinearCTStereoDble
                          , 0 /* bAllene */ );
    }
    for (i = 0; !ret && i < num_atoms; i++)
    {
        ret = FillSingleStereoDescriptors( pCG, at, (int) nAtomNumberCanon[i], -1, nCanonRank
                          , pCS->LinearCTStereoCarb, &pCS->nLenLinearCTStereoCarb, pCS->nMaxLenLinearCTStereoCarb
                          , pCS->LinearCTStereoDble, &pCS->nLenLinearCTStereoDble, pCS->nMaxLenLinearCTStereoDble
                          , 1 /* bAllene */ );
    }

    return ret;
}


/****************************************************************************
 Find stereo bond parities known in advance
****************************************************************************/
int SetKnownStereoBondParities( CANON_GLOBALS *pCG,
                                sp_ATOM       *at,
                                int           num_atoms,
                                const AT_RANK *nCanonRank,
                                const AT_RANK *nRank,
                                const AT_RANK *nAtomNumber )
{
    int i, j, n, m, j1, k, num_neigh1, num_neigh2, iMax1, parity;
    int trans_i1, trans_i2, trans_k1, trans_k2, prev_trans, trans_k, num_set;
    int i1, i2, k1, k2, n1, n2, m1, m2, /*stereo_bond_parity,*/ cumulene_len;

    AT_RANK nAtomRank1, nAtomRank2, nAtom1NeighRank;
    AT_RANK nNeighRank1[MAX_NUM_STEREO_BONDS],
        nNeighRank2[MAX_NUM_STEREO_BONDS];
    AT_RANK nNeighCanonRank1[MAX_NUM_STEREO_BONDS],
        nNeighCanonRank2[MAX_NUM_STEREO_BONDS];

    for (i1 = 0, num_set = 0; i1 < num_atoms; i1++)
    {
        if (!at[i1].parity || !at[i1].stereo_bond_neighbor[0])
        {
            continue;
        }

        if (!PARITY_WELL_DEF( at[i1].parity ))
        {
            continue;
        }

        nAtomRank1 = nRank[i1];
        iMax1 = (int) nAtomRank1 - 1;
        num_neigh1 = at[i1].valence;

        for (n1 = 0;    n1 < MAX_NUM_STEREO_BONDS &&
            ( i2 = (int) at[i1].stereo_bond_neighbor[n1] );
                           n1++)
        {
            i2--;

            /*  found a stereo bond at[i1]-at[i2] adjacent to at[i1] */
            for (n2 = 0, m = 0;
                    n2 < MAX_NUM_STEREO_BONDS &&
                    ( m = (int) at[i2].stereo_bond_neighbor[n2] ) && m - 1 != i1;
                        n2++)
                ; /* locate stereo bond (#n2) at the opposite atom at[i2] */

            if (m - 1 != i1 || at[i1].stereo_bond_parity[n1] != at[i2].stereo_bond_parity[n2])
            {
                return CT_STEREOCOUNT_ERR; /*  program error */ /*   <BRKPT> */
            }
            if (i1 < i2)
            {
                continue; /* do not process same bond 2 times */
            }
            if (PARITY_KNOWN( at[i1].stereo_bond_parity[n1] ) ||
                !PARITY_VAL( at[i1].stereo_bond_parity[n1] ))
            {
                continue;
            }
            if (!PARITY_WELL_DEF( at[i1].parity ) ||
                 !PARITY_WELL_DEF( at[i2].parity ))
            {
                continue;
            }
            if (PARITY_VAL( at[i1].stereo_bond_parity[n1] ) != AB_PARITY_CALC)
            {
                continue;  /*  ?? program error ?? should not happen */ /*   <BRKPT> */
            }

            /*stereo_bond_parity = PARITY_VAL(at[i1].stereo_bond_parity[n1]);*/
            cumulene_len = BOND_CHAIN_LEN( at[i1].stereo_bond_parity[n1] );
            nAtomRank2 = nRank[i2];
            nAtom1NeighRank = nRank[(int) at[i1].neighbor[(int) at[i1].stereo_bond_ord[n1]]];
            num_neigh2 = at[i2].valence;

            /*  store ranks of at[i1] stereo bond neighbors except one connected by a stereo bond */
            k = (int) at[i1].stereo_bond_ord[n1];
            trans_i1 = 0;
            for (i = j = 0; i < num_neigh1; i++)
            {
                if (i != k)
                {
                    nNeighRank1[j] = nRank[(int) at[i1].neighbor[i]];
                    j++;
                }
            }
            if (j == 2)
            {
                if (nNeighRank1[0] == nNeighRank1[1])
                {
                    /*  neighbors are constitutionally identical, can't find bond parity */
                    continue;
                }
                trans_i1 = insertions_sort( pCG, nNeighRank1, j, sizeof( nNeighRank1[0] ), comp_AT_RANK );
            }

            /*  store ranks of at[i2] stereo bond neighbors except one connected by a stereo bond */
            k = (int) at[i2].stereo_bond_ord[n2];
            trans_i2 = 0;
            for (i = j = 0; i < num_neigh2; i++)
            {
                if (i != k)
                {
                    nNeighRank2[j] = nRank[(int) at[i2].neighbor[i]];
                    j++;
                }
            }

            if (j == 2)
            {
                if (nNeighRank2[0] == nNeighRank2[1])
                {
                    /*  neighbors are constitutionally identical, can't find bond parity */
                    continue;
                }
                trans_i2 = insertions_sort( pCG, nNeighRank2, j, sizeof( nNeighRank2[0] ), comp_AT_RANK );
            }

            prev_trans = -1;
            trans_k1 = -2; /* djb-rwth: ignoring LLVM warning: value used */
            trans_k = -4; /* 2004-04-28 */

            /*  find all pairs of atoms that can be mapped on at[i1], at[i2] pair */
            for (j1 = 0;
                    j1 <= iMax1 && nAtomRank1 == nRank[k1 = (int) nAtomNumber[iMax1 - j1]];
                        j1++)
            {
                /*  at[k1] is constitutionally equivalent to at[i1] */
                /*  find all at[k1] neighbors that have rank nAtomRank2; */
                /*  then find at[k2] constitutionally equivalent at at[i2] */
                if (at[k1].valence != num_neigh1)
                {
                    return CT_STEREOCOUNT_ERR; /*  program error */ /*   <BRKPT> */
                }
                for (m1 = 0; m1 < num_neigh1; m1++)
                {
                    int prev, next, len;
                    if (nAtom1NeighRank != nRank[k2 = (int) at[k1].neighbor[m1]])
                    {
                        continue;
                    }
                    m2 = -1; /*  undefined yet */
                    prev = k1;
                    /* djb-rwth: removing redundant code */
                    if (cumulene_len)
                    {
                        for (len = 0, next = (int) at[k1].neighbor[m1]; len < cumulene_len; len++)
                        {
                            if (at[next].valence == 2 && !at[next].num_H)
                            {
                                j = ( (int) at[next].neighbor[0] == prev );
                                prev = next;
                                next = at[next].neighbor[j];
                            }
                            else
                            {
                                break; /*  cannot continue */
                            }
                        }
                        if (len != cumulene_len || nAtomRank2 != nRank[next])
                        {
                            continue;  /*  not found */
                        }
                        k2 = next;
                    }
                    if (at[k2].valence != num_neigh2)
                    {
                        return CT_STEREOCOUNT_ERR; /*  program error */ /*   <BRKPT> */
                    }

                    /*  store canon. ranks of at[k1] neighbors */ /*  use i,j,k,m,n */
                    for (n = j = 0; n < num_neigh1; n++)
                    {
                        if (n != m1)
                        {
                            i = (int) at[k1].neighbor[n];
                            for (m = 0; m < num_neigh1 - 1; m++)
                            {
                                if (nRank[i] == nNeighRank1[m])
                                {
                                    nNeighCanonRank1[m] = nCanonRank[i];
                                    j++;
                                    break;
                                }
                            }
                        }
                    }
                    if (j != num_neigh1 - 1)
                    {
                        return CT_STEREOCOUNT_ERR;  /*   <BRKPT> */
                    }
                    if (j == 2)
                    {
                        trans_k1 = insertions_sort( pCG, nNeighCanonRank1, j, sizeof( nNeighCanonRank1[0] ), comp_AT_RANK );
                    }
                    else
                    {
                        trans_k1 = 0;
                    }

                    /*  store canon. ranks of at[k2] neighbors */ /*  use i,j,k,m,n */
                    for (n = j = 0; n < num_neigh2; n++)
                    {
                        i = (int) at[k2].neighbor[n];
                        if (i == prev)
                        {
                            /* neighbor belongs to the stereobond */
                            m2 = n;
                        }
                        else
                        {
                            for (m = 0; m < num_neigh2 - 1; m++)
                            {
                                if (nRank[i] == nNeighRank2[m])
                                {
                                    nNeighCanonRank2[m] = nCanonRank[i];
                                    j++;
                                    break;
                                }
                            }
                        }
                    }
                    if (j != num_neigh2 - 1 || m2 < 0)
                    {
                        return CT_STEREOCOUNT_ERR;  /*   <BRKPT> */
                    }
                    if (j == 2)
                    {
                        trans_k2 = insertions_sort( pCG, nNeighCanonRank2, j, sizeof( nNeighCanonRank2[0] ), comp_AT_RANK );
                    }
                    else
                    {
                        trans_k2 = 0;
                    }
                    trans_k = ( trans_k1 + trans_k2 ) % 2;
                    if (prev_trans < 0)
                    {
                        prev_trans = trans_k;
                    }
                    else if (prev_trans != trans_k)
                    {
                        /* was != trans_k1, changed 9-23-2003 */
                        break; /*  different number of transpositions */
                    }
                } /* end of the second atom mapping cycle */
                if (prev_trans >= 0 && prev_trans != trans_k)
                { /* was != trans_k1, changed 9-23-2003 */
                    break;
                }
            } /* end of the first atom mapping cycle */

            if (prev_trans == trans_k)
            {
                /* was == trans_k1, changed 9-23-2003 */
                int z_prod;

                /*  all mappings of canonical numbers on the */
                /*  stereo bond at[i1]-at[i2] produce equivalent numberings. */
                /*  Therefore the stereo bond parity is known at this time. */
                /*  parity_1 = at[i1].parity + (trans_i1 + trans_k1 + num_neigh1 - 1) + (int)at[i1].stereo_bond_ord[n1] */
                /*  expression in parentheses is equivalent to rank[first neigh] > rank[second neigh] */
                /*  same for parity_2. */
                /*  parity_2 = at[i2].parity + (trans_i2 + trans_k2 + num_neigh2 - 1) + (int)at[i2].stereo_bond_ord[n2] */
                /*  Sum of the two parities (without stereo_bond_z_prod) is: */

                parity = ( at[i1].parity + at[i2].parity + prev_trans + trans_i1 + trans_i2
                              + num_neigh1 + num_neigh2
                              + (int) at[i1].stereo_bond_ord[n1] + (int) at[i2].stereo_bond_ord[n2] ) % 2;

                z_prod = at[i1].stereo_bond_z_prod[n1];
                if (MIN_DOT_PROD > abs( z_prod ))
                {
                    parity = AB_PARITY_UNDF; /*  undefined because of geometry */
                }
                else
                {
                    parity = ( z_prod > 0 ) ? 2 - parity : 1 + parity;
                }
                at[i1].stereo_bond_parity[n1] = ALL_BUT_PARITY( at[i1].stereo_bond_parity[n1] ) | parity;
                at[i2].stereo_bond_parity[n2] = ALL_BUT_PARITY( at[i2].stereo_bond_parity[n2] ) | parity;
                num_set++;
            }
        }
    }

    return num_set;
}


/****************************************************************************
 Find stereo center parities known in advance
****************************************************************************/
int MarkKnownEqualStereoBondParities( sp_ATOM       *at,
                                      int           num_atoms,
                                      const AT_RANK *nRank,
                                      const AT_RANK *nAtomNumber )
{
    int j, n, m, j1, num_neigh1, num_neigh2, iMax1;
    int num_set, /*num_sb1, num_sb2,*/ bDifferentParities;
    int i1, i2, k1, k2, n1, n2, m1, m2, s1, s2, stereo_bond_parity, stereo_bond_parity2, cumulene_len;
    AT_RANK nAtomRank1, nAtomRank2, nAtom1NeighRank, nAtom2NeighRank;

    /* djb-rwth: removing redundant code */

    for (i1 = 0, num_set = 0; i1 < num_atoms; i1++)
    {
        if (!at[i1].parity || !at[i1].stereo_bond_neighbor[0])
        {
            continue;
        }

        nAtomRank1 = nRank[i1];
        iMax1 = (int) nAtomRank1 - 1;
        num_neigh1 = at[i1].valence;

        /*  count stereogenic bonds adjacent to at[i1] */
        for (n1 = 0;
                n1 < MAX_NUM_STEREO_BONDS && at[i1].stereo_bond_neighbor[n1];
                    n1++);


        /*num_sb1 = n1;*/
        /*  search for bonds possibly constitutionally equivalent to each of the adjacent bonds */
        /*  and find if all of them have same already known parity */

        for (n1 = 0;
                n1 < MAX_NUM_STEREO_BONDS && ( i2 = (int) at[i1].stereo_bond_neighbor[n1] );
                    n1++) /* djb-rwth: removing redundant code */
        {
            i2--;

            nAtomRank2 = nRank[i2];
            if (nAtomRank2 < nAtomRank1 || (nAtomRank2 == nAtomRank1 && i1 < i2)) /* djb-rwth: addressing LLVM warning */
            {
                /*  An attempt to reduce unnecessary repetitions. */
                /*  We still have repetitions because we do not accumulate a list of */
                /*  processed (nAtomRank2, nAtomRank1) pairs. */
                continue;
            }

            bDifferentParities = -1;   /*  parities have not been compared yet */

            /*  found a stereo bond at[i1]-at[i2] (adjacent to at[i1]) */
            /*
            if ( !PARITY_KNOWN(at[i1].stereo_bond_parity[n1]) || (at[i1].stereo_bond_parity[n1] & KNOWN_PARITIES_EQL) )
            {
                continue;
            }
            */
            if (at[i1].stereo_bond_parity[n1] & KNOWN_PARITIES_EQL)
            {
                continue;
            }


            /*  stereo bond has known or unknown parity; we have not checked it yet */

            for (n2 = 0;
                 n2 < MAX_NUM_STEREO_BONDS && at[i2].stereo_bond_neighbor[n2];
                 n2++)
            {
                ;
            }

            /*num_sb2 = n2;*/
            for (n2 = 0, m = 0;
                 n2 < MAX_NUM_STEREO_BONDS &&
                 ( m = (int) at[i2].stereo_bond_neighbor[n2] ) &&
                    m - 1 != i1;
                 n2++)
            {
                ;
            }

            if (m - 1 != i1 || at[i1].stereo_bond_parity[n1] != at[i2].stereo_bond_parity[n2])
            {
                return CT_STEREOCOUNT_ERR; /*  program error: stereo bonds data in two directions are different */ /*   <BRKPT> */
            }

            stereo_bond_parity = PARITY_VAL( at[i1].stereo_bond_parity[n1] );
            cumulene_len = BOND_CHAIN_LEN( at[i1].stereo_bond_parity[n1] );
            nAtom1NeighRank = nRank[(int) at[i1].neighbor[(int) at[i1].stereo_bond_ord[n1]]];
            nAtom2NeighRank = nRank[(int) at[i2].neighbor[(int) at[i2].stereo_bond_ord[n2]]];
            num_neigh2 = at[i2].valence;

            /*  find all pairs of atoms that possibly can be mapped on at[i1], at[i2] pair */
            /*  (we may also find pairs that cannot be mapped, but we cannot miss any pair */
            /*  that can be mapped) */

            for (j1 = 0; j1 <= iMax1 && nAtomRank1 == nRank[k1 = (int) nAtomNumber[iMax1 - j1]]; j1++)
            {
                /*  at[k1] is constitutionally equivalent to at[i1] */
                /*  find all at[k1] stereo bond neighbors at[k2] that have rank nAtomRank2; */
                /*  then find at[k2] constitutionally equivalent at at[i2] */

                if (at[k1].valence != num_neigh1)
                {
                    return CT_STEREOCOUNT_ERR; /*  program error */ /*   <BRKPT> */
                }

                if (!at[k1].bHasStereoOrEquToStereo)
                {
                    at[k1].bHasStereoOrEquToStereo = 1;
                }

                /* -- do not check number of stereo bonds, check bonds themselves --
                for ( s1 = 0; s1 < MAX_NUM_STEREO_BONDS && at[k1].stereo_bond_neighbor[s1]; s1++ )
                {
                    ;
                }
                if ( num_sb1 != s1 )
                {
                    bDifferentParities = 1;
                }
                */

                for (m1 = 0; m1 < num_neigh1; m1++)
                {
                    /*  Looking for at[k1] neighbor with nRank=nAtom1NeighRank. */
                    /*  This neighbor may be on the bond constit. equivalent to at[i1]-at[i2] stereo bond */
                    /*  (or may be constit. equivalent an adjacent to at[i1] atom in a stereogenic cumulene chain) */
                    int prev, next, len;
                    if (nAtom1NeighRank != nRank[k2 = (int) at[k1].neighbor[m1]])
                        continue;

                    /*  found at[k1] neighbor with nRank=nAtom1NeighRank */

                    m2 = -1; /*  undefined yet */
                    prev = k1;
                    /* djb-rwth: removing redundant code */

                    /*  if cumulene then bypass the cumulene chain */

                    if (cumulene_len)
                    {

                        for (len = 0, next = (int) at[k1].neighbor[m1]; len < cumulene_len; len++)
                        {
                            if (at[next].valence == 2 && !at[next].num_H)
                            {
                                j = ( (int) at[next].neighbor[0] == prev );
                                prev = next;
                                next = at[next].neighbor[j];
                            }
                            else
                            {
                                break; /*  cannot continue: end of cumulene chain */
                            }
                        }

                        if (len != cumulene_len || nAtomRank2 != nRank[next])
                        {
                            continue;  /*  cumulene chain not found at this neighbor */
                        }

                        if (nAtom2NeighRank != nRank[prev])
                        {
                            /* continue; */ /*  ??? program error ??? If not, must be a very rare event */
                            return CT_STEREOCOUNT_ERR;  /*   <BRKPT> */
                        }

                        k2 = next;
                    }

                    /*  a connected pair of constit. equivalent atoms found */

                    if (at[k2].valence != num_neigh2)
                    {
                        return CT_STEREOCOUNT_ERR; /*  program error */ /*   <BRKPT> */
                    }

                    for (n = 0; n < num_neigh2; n++)
                    {
                        if (prev == (int) at[k2].neighbor[n])
                        {
                            m2 = n; /*  found bond from the opposite end of a possibly stereogenic bond */
                            break;
                        }
                    }

                    if (m2 < 0)
                    {
                        return CT_STEREOCOUNT_ERR; /*  program error: opposite direction bond not found */ /*   <BRKPT> */
                    }

                    if (!at[k2].bHasStereoOrEquToStereo)
                    {
                        at[k2].bHasStereoOrEquToStereo = 1;
                    }


                    /*  check if atoms at[k1] and at[k2] are connected by a stereo bond */
                    for (s1 = 0, m = 0;
                         s1 < MAX_NUM_STEREO_BONDS &&
                         ( m = (int) at[k1].stereo_bond_neighbor[s1] ) &&
                                    m - 1 != k2;
                         s1++)
                    {
                        ;
                    }
                    if (m - 1 != k2)
                    {
                        bDifferentParities = 1; /*  cannot find the stereo bond */
                        at[k1].bHasStereoOrEquToStereo =
                            at[k2].bHasStereoOrEquToStereo = 2;
                        continue;
                    }

                    /*  -- do not check number of stereo bonds, check bonds themselves --
                    for ( s2 = 0; s2 < MAX_NUM_STEREO_BONDS && at[k2].stereo_bond_neighbor[s2]; s2++ )
                     {
                        ;
                     }
                    if ( num_sb2 != s2 )
                    {
                        bDifferentParities = 1;
                        continue;
                    }
                    */

                    for (s2 = 0, m = 0;
                         s2 < MAX_NUM_STEREO_BONDS &&
                         ( m = (int) at[k2].stereo_bond_neighbor[s2] ) &&
                                    m - 1 != k1;
                         s2++)
                    {
                        ;
                    }

                    if (m - 1 != k1)
                    {
                        /*
                        bDifferentParities = 1; // cannot find the stereo bond
                        continue;
                        */
                        return CT_STEREOCOUNT_ERR; /*  program error: opposite direction bond not found */ /*   <BRKPT> */
                    }

                    if (at[k1].stereo_bond_parity[s1] != at[k2].stereo_bond_parity[s2])
                    {
                        bDifferentParities = 1;
                        continue;
                    }
                    stereo_bond_parity2 = PARITY_VAL( at[k1].stereo_bond_parity[s1] );
                    if (stereo_bond_parity2 != stereo_bond_parity)
                    {
                        bDifferentParities = 1;
                        continue;
                    }
                    if (stereo_bond_parity2 == stereo_bond_parity && bDifferentParities < 0)
                    {
                        bDifferentParities = 0;
                    }
                }
            }

            /*  mark equal parities */
            if (0 == bDifferentParities && PARITY_KNOWN( stereo_bond_parity ))
            {
                for (j1 = 0;
                        j1 <= iMax1 && nAtomRank1 == nRank[k1 = (int) nAtomNumber[iMax1 - j1]];
                            j1++)
                {
                    /*  at[k1] is constitutionally equivalent to at[i1] */
                    for (s1 = 0; s1 < MAX_NUM_STEREO_BONDS && ( k2 = (int) at[k1].stereo_bond_neighbor[s1] ); s1++) /* djb-rwth: removing redundant code */
                    {
                        k2--;
                        if (nRank[k2] == nAtomRank2)
                        {
                            int b1, b2;
                            for (s2 = 0, m = 0;
                                 s2 < MAX_NUM_STEREO_BONDS &&
                                    ( m = (int) at[k2].stereo_bond_neighbor[s2] ) &&
                                    m - 1 != k1;
                                 s2++)
                            {
                                ;
                            }

                            if (m - 1 != k1)
                            {
                                return CT_STEREOCOUNT_ERR; /*  program error */ /*   <BRKPT> */
                            }
                            /*  mark the stereo bonds */
                            b1 = !( at[k1].stereo_bond_parity[s1] & KNOWN_PARITIES_EQL );
                            b2 = !( at[k2].stereo_bond_parity[s2] & KNOWN_PARITIES_EQL );
                            if (2 == b1 + b2)
                            {
                                at[k1].stereo_bond_parity[s1] |= KNOWN_PARITIES_EQL;
                                at[k2].stereo_bond_parity[s2] |= KNOWN_PARITIES_EQL;
                                num_set++;
                            }
                            else if (b1 || b2)
                            {
                                return CT_STEREOCOUNT_ERR; /*  program error */ /*   <BRKPT> */
                            }
                        }
                    }
                }
            }
        }
    }

    return num_set;
}



#if ( REMOVE_KNOWN_NONSTEREO == 1 ) /* { */


/****************************************************************************
 Return next atom number (and its canon. rank) on the path prev->cur->next
 in order of ascending canonical ranks of the next atoms:
    *cr(output) >  *cr(input)
 To start the sequence let *cr=0
 If no more neighbors available the return value = 0; on success, returns 1.
****************************************************************************/
int GetNextNeighborAndRank( sp_ATOM       *at,
                            AT_RANK       cur,
                            AT_RANK       prev,
                            AT_RANK       *n,
                            AT_RANK       *cr,
                            const AT_RANK *nCanonRank )
{
    int i, val;
    AT_RANK cr1 = MAX_ATOMS + 1, j, j1 = MAX_ATOMS + 1, crj;

    for (i = 0, val = at[(int) cur].valence; i < val; i++)
    {
        if (( j = at[cur].neighbor[i] ) != prev &&
             cr1 > ( crj = nCanonRank[(int) j] ) && crj > *cr)
        {
            cr1 = crj;
            j1 = j;
        }
    }
    if (cr1 <= MAX_ATOMS)
    {
        *cr = cr1;
        *n = (AT_RANK) j1;
        return 1;
    }

    return 0;  /*  program error */ /*   <BRKPT> */
}


/****************************************************************************
 Find next pair of neighbors having the next greater canon. rank
 The neighbors should be constitutionally identical and traversed
 simultaneouslyor not traversed at all
 If a bond cur1-*n1 or cur2-*n2 is a stereo bond then reject if their stereo
 bond parities are different or cannot be calculated without breaking ties.
****************************************************************************/
int GetAndCheckNextNeighbors( sp_ATOM       *at,
                              AT_RANK       cur1,
                              AT_RANK       prev1,
                              AT_RANK       cur2,
                              AT_RANK       prev2,
                              AT_RANK       *n1,
                              AT_RANK       *n2,
                              AT_RANK       *nVisited1,
                              AT_RANK       *nVisited2,
                              const AT_RANK *nRank,
                              const AT_RANK *nCanonRank )
{
    AT_RANK cr1, cr2, s1, s2;
    int     i1, i2, k1, k2;

    cr1 = ( *n1 > MAX_ATOMS ) ? 0 : nCanonRank[(int) *n1];
    cr2 = ( *n2 > MAX_ATOMS ) ? 0 : nCanonRank[(int) *n2];

    if (!GetNextNeighborAndRank( at, cur1, prev1, n1, &cr1, nCanonRank ) ||
         !GetNextNeighborAndRank( at, cur2, prev2, n2, &cr2, nCanonRank ) ||
         nRank[(int) *n1] != nRank[(int) *n2] || nVisited1[(int) *n1] != nVisited2[(int) *n2])
    {
        return 0;  /*  program error; no breakpoint here */ /*   <BRKPT> */
    }

    /*  Even though the bond or cumulene might have already been checked, check it: this is */
    /*  the only place we can check stereo bonds and cumulenes that are not edges of the DFS tree */
    /*  The code works both for a stereo bond and a stereogenic cumulene. */
    for (i1 = 0, k1 = 0;
         i1 < MAX_NUM_STEREO_BONDS &&
         ( s1 = at[cur1].stereo_bond_neighbor[i1] ) &&
            !( k1 = ( at[cur1].neighbor[(int) at[cur1].stereo_bond_ord[i1]] == *n1 ) );
         i1++) /* djb-rwth: ignoring LLVM warning: variable used */
    {
        ;
    }

    for (i2 = 0, k2 = 0;
         i2 < MAX_NUM_STEREO_BONDS &&
         ( s2 = at[cur2].stereo_bond_neighbor[i2] ) &&
            !( k2 = ( at[cur2].neighbor[(int) at[cur2].stereo_bond_ord[i2]] == *n2 ) );
         i2++) /* djb-rwth: ignoring LLVM warning: variable used */
    {
        ;
    }

    if (k1 != k2)
    {
        return 0; /*  possibly not an error: constit. equivalent atoms on a stereo bond and not on a stereo bond */
    }

    if (k1 /* yes, it is a stereo bond */ &&
        ( at[cur1].stereo_bond_parity[i1] != at[cur2].stereo_bond_parity[i2] ||
          /* PARITY_KNOWN (at[cur1].stereo_bond_parity[i1] ) */  /*  replaced 08-13-2002 with the next: */
            !PARITY_WELL_DEF( at[cur1].stereo_bond_parity[i1] ) /*  it suffices to check only one parity */
            ))
    {
        return 0; /*  different or (currently) unknown stereo bond parities */
    }

    return 1; /*  stereo bonds have known parities */
}


/********************************************************************************************/
/*  Simultaneously DFS-traverse 2 paths starting at the bonds prev1->cur1 and prev2->cur2   */
/*  The two paths MUST go through the pairs of constitutionally identical atoms,            */
/*  each atom being on one path.                                                            */
/*  Reject if encountered atoms having currently unknown (without breaking ties)            */
/*      parities or having different known or unknown or undefined parities.                */
/*  Save length of the path into nVisited[cur. atom number].                                */
/*  Only one nVisited[] array is sufficient because the paths from the beginning are        */
/*  in different ring systems.                                                              */
/********************************************************************************************/
AT_RANK PathsHaveIdenticalKnownParities( sp_ATOM       *at,
                                         AT_RANK       prev1,
                                         AT_RANK       cur1,
                                         AT_RANK       prev2,
                                         AT_RANK       cur2,
                                         AT_RANK       *nVisited1,
                                         AT_RANK       *nVisited2,
                                         const AT_RANK *nRank,
                                         const AT_RANK *nCanonRank,
                                         AT_RANK       nLength )
{
    int k;
    AT_RANK n1, n2;

    nLength++;   /*  number of successfully traversed pairs of atoms */
    nVisited1[cur1] = nLength;
    nVisited2[cur2] = nLength;

    /*  the atoms must be either both stereogenic and have well-defined parities or non-stereogenic at all. */
    if (at[cur1].stereo_atom_parity != at[cur2].stereo_atom_parity ||
         (at[cur1].stereo_atom_parity && !PARITY_WELL_DEF( at[cur1].stereo_atom_parity )) ) /* djb-rwth: addressing LLVM warning */
    {
        return 0;  /*  Reject: Different or unknown in advance parities */
    }

    if (at[cur1].valence != at[cur2].valence)
    {
        return 0;  /*  program error */ /*   <BRKPT> */
    }

    if (at[cur1].valence == 1)
    {
        return nLength; /*  so far success */
    }


    for (k = 1, n1 = MAX_ATOMS + 1, n2 = MAX_ATOMS + 1; k < at[cur1].valence; k++)
    {
        /*  start from 1: since we do not go back, we have only (at[cur1].valence-1) bonds to try */

        if (!GetAndCheckNextNeighbors( at, cur1, prev1, cur2, prev2,
            &n1, &n2, nVisited1, nVisited2,
            nRank, nCanonRank ))
        {
            return 0; /*  different neighbors                       */
        }

        /*  In a DFS we do not traverse already visited atoms */
        if (!nVisited1[n1])
        {
            /*  recursion */
            if (!( nLength = PathsHaveIdenticalKnownParities( at, cur1, n1, cur2, n2, nVisited1, nVisited2, nRank, nCanonRank, nLength ) ))
            {
                return 0;
            }
        }
    }

    /*  To be on a safe side, recheck after all nVisited[] have been set */
    for (k = 1, n1 = MAX_ATOMS + 1, n2 = MAX_ATOMS + 1; k < at[cur1].valence; k++)
    {
        /*  start from 1: since we do not go back, we have only (at[cur1].valence-1) bonds to try */
        if (!GetAndCheckNextNeighbors( at, cur1, prev1, cur2, prev2,
            &n1, &n2, nVisited1, nVisited2,
            nRank, nCanonRank ))
        {
            return 0; /*  different neighbors */
        }
    }

    return nLength;
}


/****************************************************************************/
/*  Remove stereo marks from the bonds that are known to be non-stereo      */
/*  (compare neighbors if they are attached by cut-edges)                   */
/****************************************************************************/
int RemoveKnownNonStereoBondParities( sp_ATOM       *at,
                                      int           num_atoms,
                                      const AT_RANK *nCanonRank,
                                      const AT_RANK *nRank,
                                      CANON_STAT    *pCS )
{
    int j, n, m, ret;

    int i1, n1, s2;
    AT_RANK nAtomRank1, nAtomRank2, neigh[3], opposite_atom, *nVisited = NULL;
    ret = 0;
    for (i1 = 0; i1 < num_atoms; i1++)
    {
        if (at[i1].valence != 3 || !at[i1].stereo_bond_neighbor[0])
        {
            continue;
        }
        for (n1 = 0; n1 < MAX_NUM_STEREO_BONDS && ( s2 = at[i1].stereo_bond_neighbor[n1] ); n1++)
        {
            if (!PARITY_CALCULATE( at[i1].stereo_bond_parity[n1] ) && PARITY_WELL_DEF( at[i1].stereo_bond_parity[n1] ))
            {
                continue;
            }
            opposite_atom = (AT_RANK) ( s2 - 1 );
            /* s2 = at[i1].neighbor[m=(int)at[i1].stereo_bond_ord[n1]]; */
            m = (int) at[i1].stereo_bond_ord[n1];
            for (j = 0, n = 0; j < at[i1].valence; j++)
            {
                /* if ( at[i1].neighbor[j] != s2 ) */
                if (j != m)
                {
                    neigh[n++] = at[i1].neighbor[j];
                }
            }
            if (n > 2)
            {
                ret = CT_STEREOBOND_ERROR;  /*   <BRKPT> */
                goto exit_function;
            }
            if (n != 2 || nRank[(int) neigh[0]] != nRank[(int) neigh[1]])
            {
                continue; /*  may happen if another half-bond has not a defined parity */
            }
            if (at[i1].nRingSystem == at[(int) neigh[0]].nRingSystem)
            {
                continue;  /*  no more ring system membership check is necessary because     */
            }              /*  the two neighbors are to be constitutionally equivalent atoms */
            if (!nVisited && !( nVisited = (AT_RANK*) inchi_malloc( sizeof( nVisited[0] )*num_atoms ) ))
            {
                ret = CT_OUT_OF_RAM;  /*   <BRKPT> */
                goto exit_function;
            }
            memset( nVisited, 0, sizeof( nVisited[0] )*num_atoms ); /* djb-rwth: memset_s C11/Annex K variant? */
            nVisited[i1] = 1;
            if (PathsHaveIdenticalKnownParities( at, (AT_RANK) i1, neigh[0], (AT_RANK) i1, neigh[1], nVisited, nVisited, nRank, nCanonRank, 1 ))
            {
                if (!RemoveOneStereoBond( at, i1, /* atom number*/ n1 /* stereo bond number*/ ))
                {
                    ret = CT_STEREOBOND_ERROR;  /*   <BRKPT> */
                    goto exit_function;
                }
                n1--; /*  cycle counter may temporarily become negative */
                /*  Remove from pCS */
                nAtomRank1 = inchi_max( nCanonRank[i1], nCanonRank[opposite_atom] );
                nAtomRank2 = inchi_min( nCanonRank[i1], nCanonRank[opposite_atom] );
                for (n = 0, m = pCS->nLenLinearCTStereoDble - 1; n <= m; n++)
                {
                    if (pCS->LinearCTStereoDble[n].at_num1 == nAtomRank1 &&
                         pCS->LinearCTStereoDble[n].at_num2 == nAtomRank2)
                    {
                        if (n < m)
                        { /*  remove pCS->LinearCTStereoDble[n] */
                            memmove( pCS->LinearCTStereoDble + n, pCS->LinearCTStereoDble + n + 1, ( (long long)m - (long long)n ) * sizeof( pCS->LinearCTStereoDble[0] ) ); /* djb-rwth: cast operators added */
                        }
                        pCS->nLenLinearCTStereoDble--;
#if ( bRELEASE_VERSION == 0 )
                        pCS->bExtract |= EXTR_KNOWN_USED_TO_REMOVE_PARITY;
#endif
                        m = -1;   /*  set flag "found" */
                        break;
                    }
                }
                if (m >= 0)
                {
                    ret = CT_STEREOCOUNT_ERR;  /*  bond not found  <BRKPT> */
                    goto exit_function;
                }
                ret++;  /*  number of removed known in advance non-stereo bonds */
            }
        }
    }

exit_function:
    if (nVisited)
    {
        inchi_free( nVisited );
    }

    return ret;
}
#endif /* } REMOVE_KNOWN_NONSTEREO */


/****************************************************************************
 Find stereo center parities known in advance
****************************************************************************/
int SetKnownStereoCenterParities( CANON_GLOBALS *pCG,
                                  sp_ATOM       *at,
                                  int           num_atoms,
                                  const AT_RANK *nCanonRank,
                                  const AT_RANK *nRank,
                                  const AT_RANK *nAtomNumber )
{
    int i, j, n, m, j1, k, num_neigh, iMax, trans_i, trans_k, prev_trans, num_set;
    AT_RANK nAtomRank;
    AT_RANK nNeighRank[MAX_NUM_STEREO_ATOM_NEIGH];
    AT_RANK nNeighCanonRank[MAX_NUM_STEREO_ATOM_NEIGH];

    for (i = 0, num_set = 0; i < num_atoms; i++)
    {
        if (!at[i].parity || at[i].stereo_bond_neighbor[0])
        {
            continue;
        }
        if (at[i].stereo_atom_parity != AB_PARITY_CALC ||
             !PARITY_WELL_DEF( at[i].parity ))
        {
            continue;
        }
        num_neigh = at[i].valence;
        for (j = 0; j < num_neigh; j++)
        {
            nNeighRank[j] = nRank[(int) at[i].neighbor[j]];
        }
        nAtomRank = nRank[i];
        if (num_neigh == 1)
        {
            /* other neighbors must be implicit H */
            at[i].stereo_atom_parity = at[i].parity;
            trans_i = 0;
        }
        else
        {
            /* sort constitutional equivalence ranks of the neighbors */
            trans_i = insertions_sort( pCG, nNeighRank, num_neigh, sizeof( nNeighRank[0] ), comp_AT_RANK );
            for (j = 1; j < num_neigh; j++)
            {
                if (nNeighRank[j - 1] == nNeighRank[j])
                {
                    break; /* at[i] has consitutionally identical neighbors */
                }
            }
            if (j < num_neigh)
            {
                /*  at least 2 neighbors are const. identical; parity cannot be calculated at this time */
                continue; /*  try next stereo atom */
            }
        }
        prev_trans = -1;
        trans_k = 0;
        /*  find neighbors of constitutionally equivalent stereo centers */
        /*  and at[i] parities in case those centers are mapped on at[i] */
        for (iMax = (int) nAtomRank - 1, j1 = 0;
                j1 <= iMax && nAtomRank == nRank[k = (int) nAtomNumber[iMax - j1]];
                    j1++)
        {
            /*  at[k] is constitutionally equivalent to at[i] */
            if ((int) at[k].valence != num_neigh)
            {
                return CT_STEREOCOUNT_ERR;  /*   <BRKPT> */
            }
            /* -- commented out to accept  non-stereogenic atoms since     --
               -- they may participate in mapping stereocenters 12-16-2003 --
            if ( !PARITY_VAL(at[k].parity) ) {
                continue; // not a stereogenic atom
            }
            */
            for (j = 0, m = 0; m < num_neigh; m++)
            {
                for (n = 0; n < num_neigh; n++)
                {
                    if (nRank[(int) at[k].neighbor[n]] == nNeighRank[m])
                    {
                        /* save canonical numbers (ranks) of the neighbors in
                         * order of increasing constit. equivalence ranks */
                        nNeighCanonRank[m] = nCanonRank[(int) at[k].neighbor[n]];
                        j++;
                        break;
                    }
                }
            }
            if (j != num_neigh)
            {
                return CT_STEREOCOUNT_ERR;  /*   <BRKPT> */
            }
            trans_k = insertions_sort( pCG, nNeighCanonRank, num_neigh, sizeof( nNeighCanonRank[0] ), comp_AT_RANK );
            trans_k %= 2;
            if (prev_trans < 0)
            {
                prev_trans = trans_k;
            }
            else if (trans_k != prev_trans)
            {
                /*  different mappings may produce different parities. Cannot find the parity at this time */
                /*  this may happen when a set of constit. equivalent atoms has non-contiguous canonical numbers */
                break;
            }
        }
        if (trans_k == prev_trans)
        {
            at[i].stereo_atom_parity = 2 - ( at[i].parity + trans_i + prev_trans ) % 2;
            num_set++;
        }
    }

    return num_set;
}


#if ( REMOVE_KNOWN_NONSTEREO == 1 ) /* { */


/****************************************************************************
 DFS along paths starting from the stereocenter through pairs of cut-edges
****************************************************************************/
int RemoveKnownNonStereoCenterParities( CANON_GLOBALS *pCG,
                                        sp_ATOM       *at,
                                        int           num_atoms,
                                        const AT_RANK *nCanonRank,
                                        const AT_RANK *nRank,
                                        CANON_STAT    *pCS )
{
    int i, j, n, m, k, num_neigh, ret = 0;
    /*AT_RANK nAtomRank;*/
    AT_RANK nNeighRank[MAX_NUM_STEREO_ATOM_NEIGH], nNeighOrd[MAX_NUM_STEREO_ATOM_NEIGH];
    AT_RANK *nVisited = NULL;

    for (i = 0; i < num_atoms; i++)
    {
        if (!at[i].parity || at[i].stereo_bond_neighbor[0])
        {
            continue;
        }
        if (!PARITY_CALCULATE( at[i].stereo_atom_parity ) && PARITY_WELL_DEF( at[i].stereo_atom_parity ))
        {
            continue;
        }
        num_neigh = at[i].valence;
        for (j = 0; j < num_neigh; j++)
        {
            nNeighRank[j] = nRank[(int) at[i].neighbor[j]];
            nNeighOrd[j] = j;
        }
        /*nAtomRank = nRank[i];*/
        if (num_neigh == 1)
        {
            continue;
        }
        pCG->m_pn_RankForSort = nNeighRank;
        insertions_sort( pCG, nNeighOrd, num_neigh, sizeof( nNeighRank[0] ), CompRanksOrd );
        for (j = k = 1; k && j < num_neigh; j++)
        {
            if (at[i].nRingSystem != at[(int) at[i].neighbor[(int) nNeighOrd[j]]].nRingSystem &&
                 /*  no more ring system membership check is necessary because */
                 /*  the two neighbors are to be constitutionally equivalent atoms: */
                nNeighRank[nNeighOrd[j - 1]] == nNeighRank[nNeighOrd[j]])
            {
                k = j;
                do
                {
                    if (!nVisited && !( nVisited = (AT_RANK*) inchi_malloc( sizeof( nVisited[0] )*num_atoms ) ))
                    {
                        ret = CT_OUT_OF_RAM;  /*   <BRKPT> */
                        goto exit_function;
                    }
                    memset( nVisited, 0, sizeof( nVisited[0] )*num_atoms ); /* djb-rwth: memset_s C11/Annex K variant? */
                    nVisited[i] = 1;
                    if (PathsHaveIdenticalKnownParities( at, (AT_RANK) i, at[i].neighbor[(int) nNeighOrd[j - 1]],
                        (AT_RANK) i, at[i].neighbor[(int) nNeighOrd[k]],
                        nVisited, nVisited, nRank, nCanonRank, 1 ))
                    {
                        at[i].parity = 0; /*  remove parity */
                        at[i].stereo_atom_parity = 0;
                        at[i].final_parity = 0;
                        /* at[i].bHasStereoOrEquToStereo = 0; */
                        for (n = 0, m = pCS->nLenLinearCTStereoCarb - 1; n <= m; n++)
                        {
                            if (pCS->LinearCTStereoCarb[n].at_num == nCanonRank[i])
                            {
                                if (n < m)
                                {    /*  remove pCS->LinearCTStereoCarb[n] */
                                    memmove( pCS->LinearCTStereoCarb + n, pCS->LinearCTStereoCarb + n + 1, ( (long long)m - (long long)n ) * sizeof( pCS->LinearCTStereoCarb[0] ) ); /* djb-rwth: cast operators added */
                                }
                                pCS->nLenLinearCTStereoCarb--;
                                k = 0;

#if ( bRELEASE_VERSION == 0 )
                                pCS->bExtract |= EXTR_KNOWN_USED_TO_REMOVE_PARITY;
#endif

                                break;
                            }
                        }
                        if (k)
                        {
                            ret = CT_STEREOCOUNT_ERR;  /*   <BRKPT> */
                            goto exit_function;
                        }
                        ret++;
                        break;
                    }
                }
                while (++k < num_neigh && nNeighRank[nNeighOrd[j - 1]] == nNeighRank[nNeighOrd[k]]);
            }
        }
    }

exit_function:
    if (nVisited)
    {
        inchi_free( nVisited );
    }

    return ret;
}



#endif /* } REMOVE_KNOWN_NONSTEREO */


/****************************************************************************
 Find stereo center parities known in advance
****************************************************************************/
int MarkKnownEqualStereoCenterParities( sp_ATOM       *at,
                                        int           num_atoms,
                                        const AT_RANK *nRank,
                                        const AT_RANK *nAtomNumber )
{
    int i, j1, k, num_centers, iMax, bDifferentParities;
    AT_RANK nAtomRank;
    int  parity, parity_k;

    num_centers = 0;
    for (i = 0; i < num_atoms; i++)
    {
        if (!at[i].parity || at[i].stereo_bond_neighbor[0])
        {
            continue;
        }
        if (at[i].bHasStereoOrEquToStereo)
        {
            continue; /* already marked */
        }
        if ( /*!PARITY_KNOWN(at[i].stereo_atom_parity) ||*/ ( at[i].stereo_atom_parity & KNOWN_PARITIES_EQL ))
        {
            continue;
        }
        parity = PARITY_VAL( at[i].stereo_atom_parity );
        if (parity == AB_PARITY_NONE)
        {
            continue;
        }
        nAtomRank = nRank[i];
        bDifferentParities = -1;
        /*  find constitutionally equivalent stereo centers and compare their known at this time parities */
        for (iMax = (int) nAtomRank - 1, j1 = 0; j1 <= iMax && nAtomRank == nRank[k = (int) nAtomNumber[iMax - j1]]; j1++)
        {
            /*  at[k] is constitutionally equivalent to at[i] */
            parity_k = PARITY_VAL( at[k].stereo_atom_parity );
            if (parity_k != parity)
            {
                bDifferentParities = 1;
            }
            else if (parity_k == parity && bDifferentParities < 0)
            {
                bDifferentParities = 0;
            }
            if (!parity_k)
            {
                at[k].bHasStereoOrEquToStereo = 2;
            }
            else if (!at[k].bHasStereoOrEquToStereo)
            {
                at[k].bHasStereoOrEquToStereo = 1;
            }
        }
        if (0 == bDifferentParities && PARITY_KNOWN( parity ))
        {
            for (iMax = (int) nAtomRank - 1, j1 = 0; j1 <= iMax && nAtomRank == nRank[k = (int) nAtomNumber[iMax - j1]]; j1++)
            {
                /*  at[k] is constitutionally equivalent to at[i] */
                at[k].stereo_atom_parity |= KNOWN_PARITIES_EQL;
                num_centers++;
            }
        }
    }

    return num_centers;
}


/*****************************************************************************/
/* invert known parities in at[] and in pCS->LinearCTStereoDble              */
/*                                      pCS->LinearCTStereoCarb              */
/* nCanonRank[] contains canonical ranks used to fill pCS->LinearCTStereo... */
/* nAtomNumberCanon[] will be filled with atom numbers in canonical order    */
/*****************************************************************************/
int InvertStereo( sp_ATOM    *at,
                  int        num_at_tg,
                  AT_RANK    *nCanonRank,
                  AT_RANK    *nAtomNumberCanon,
                  CANON_STAT *pCS,
                  int        bInvertLinearCTStereo )
{
    int i, j, j1, j2, num_changes, parity, cumulene_len;

    num_changes = 0;
    for (i = 0; i < num_at_tg; i++)
    {
        nAtomNumberCanon[(int) nCanonRank[i] - 1] = i;
    }
    for (i = 0; i < pCS->nLenLinearCTStereoCarb; i++)
    {
        parity = pCS->LinearCTStereoCarb[i].parity;
        if (ATOM_PARITY_WELL_DEF( parity ))
        {
            j = nAtomNumberCanon[(int) pCS->LinearCTStereoCarb[i].at_num - 1];
            if (PARITY_WELL_DEF( at[j].parity ))
            {
                at[j].parity ^= AB_INV_PARITY_BITS;
            }
            else
            {
                goto exit_error; /* inconsistency */
            }
            if (bInvertLinearCTStereo)
            {
#ifdef FIX_STEREOCOUNT_ERR
                pCS->LinearCTStereoCarb[i].parity = AB_PARITY_EVEN; /* deliberately worse */
#else
                pCS->LinearCTStereoCarb[i].parity ^= AB_INV_PARITY_BITS;
#endif
            }
            num_changes++;
            if (PARITY_WELL_DEF( at[j].stereo_atom_parity ))
            {
                at[j].stereo_atom_parity ^= AB_INV_PARITY_BITS;
            }
            if (PARITY_WELL_DEF( at[j].final_parity ))
            {
                at[j].final_parity ^= AB_INV_PARITY_BITS;
            }
        }
    }
    for (i = 0; i < pCS->nLenLinearCTStereoDble; i++)
    {
        parity = pCS->LinearCTStereoDble[i].parity;
        if (ATOM_PARITY_WELL_DEF( parity ))
        {
            j1 = nAtomNumberCanon[(int) pCS->LinearCTStereoDble[i].at_num1 - 1];
            cumulene_len = BOND_CHAIN_LEN( at[j1].stereo_bond_parity[0] );
            if (cumulene_len % 2)
            {
                /* invert only in case of allene */
                j2 = nAtomNumberCanon[(int) pCS->LinearCTStereoDble[i].at_num2 - 1];
                /* checks for debug only */
                if (1 < MAX_NUM_STEREO_BONDS)
                {
                    if (at[j1].stereo_bond_neighbor[1] ||
                         at[j2].stereo_bond_neighbor[1])
                    {
                        goto exit_error; /* inconsitency: atom has more than one cumulene bond */
                    }
                }
                if (cumulene_len != BOND_CHAIN_LEN( at[j2].stereo_bond_parity[0] ) ||
                     j1 + 1 != at[j2].stereo_bond_neighbor[0] ||
                     j2 + 1 != at[j1].stereo_bond_neighbor[0])
                {
                    goto exit_error; /* inconsitency: atoms should refer to each other */
                }
                /* invert parities */
                if (PARITY_WELL_DEF( at[j1].parity ) && PARITY_WELL_DEF( at[j2].parity ))
                {
                    j = inchi_min( j1, j2 );
                    at[j].parity ^= AB_INV_PARITY_BITS; /* for reversability always invert only atom with the smaller number */
                }
                else
                {
                    goto exit_error; /* inconsistency */
                }
                if (bInvertLinearCTStereo)
                {
                    pCS->LinearCTStereoDble[i].parity ^= AB_INV_PARITY_BITS;
                }
                num_changes++;
                if (PARITY_WELL_DEF( at[j1].stereo_bond_parity[0] ))
                {
                    at[j1].stereo_bond_parity[0] ^= AB_INV_PARITY_BITS;
                }
                if (PARITY_WELL_DEF( at[j2].stereo_bond_parity[0] ))
                {
                    at[j2].stereo_bond_parity[0] ^= AB_INV_PARITY_BITS;
                }
            }
        }
    }

    return num_changes;

exit_error:
    return CT_STEREOCOUNT_ERR;
}


/********************************************************************************/
/*  Make sure atoms stereo descriptors fit molecular symmetry and remove        */
/*  parity from obviously non-stereo atoms and bonds                            */
/********************************************************************************/
int FillOutStereoParities( sp_ATOM       *at,
                           int           num_atoms,
                           const AT_RANK *nCanonRank,
                           const AT_RANK *nAtomNumberCanon,
                           const AT_RANK *nRank,
                           const AT_RANK *nAtomNumber,
                           CANON_STAT    *pCS,
                           CANON_GLOBALS *pCG,
                           int           bIsotopic )
{
    int ret;
    /*  unmark atoms with 2 or more constitutionally equivalent neighbors */
    /*  such that there is no path through them to an atom with parity */
    ret = UnmarkNonStereo( pCG, at, num_atoms, nRank, nAtomNumber, bIsotopic );
    if (ret < 0)
        return ret;  /*  program error? */ /*   <BRKPT> */
    ret = FillAllStereoDescriptors( pCG, at, num_atoms, nCanonRank, nAtomNumberCanon, pCS ); /*  ret<0: error */
    if (!ret)
    {
        ret = pCS->nLenLinearCTStereoCarb + pCS->nLenLinearCTStereoDble;
    }
    if (ret < 0)
    {
        return ret; /*  program error? */ /*   <BRKPT> */
    }

    if (ret >= 0)
    {
        int ret2;
        ret2 = SetKnownStereoCenterParities( pCG, at, num_atoms, nCanonRank, nRank, nAtomNumber );
        if (ret2 >= 0)
        {
            ret2 = MarkKnownEqualStereoCenterParities( at, num_atoms, nRank, nAtomNumber );
        }
        if (ret2 >= 0)
        {
            ret2 = SetKnownStereoBondParities( pCG, at, num_atoms, nCanonRank, nRank, nAtomNumber );
            if (ret2 >= 0)
            {
                ret2 = MarkKnownEqualStereoBondParities( at, num_atoms, nRank, nAtomNumber );
            }
        }
#if ( REMOVE_KNOWN_NONSTEREO == 1 ) /* { */
        if (ret2 >= 0)
        {
            int ret3;
            do
            {
                ret2 = RemoveKnownNonStereoCenterParities( pCG, at, num_atoms, nCanonRank, nRank, pCS );
                if (ret2 >= 0)
                {
                    ret3 = RemoveKnownNonStereoBondParities( at, num_atoms, nCanonRank, nRank, pCS );
                    ret2 = ret3 >= 0 ? ret2 + ret3 : ret3;
                }
            }
            while (ret2 > 0);
        }
        if (RETURNED_ERROR( ret2 ))
        {
            ret = ret2;
        }
#endif /* } REMOVE_KNOWN_NONSTEREO */
    }

    return ret; /*  non-zero means error */
}



/****************************************************************************/
int GetStereoNeighborPos( sp_ATOM *at,
                          int     iAt1,
                          int     iAt2 )
{
    int k1;
    AT_RANK sNeigh = (AT_RANK) ( iAt2 + 1 );
    AT_RANK s;
    for (k1 = 0; k1 < MAX_NUM_STEREO_BONDS && ( s = at[iAt1].stereo_bond_neighbor[k1] ); k1++)
    {
        if (s == sNeigh)
        {
            return k1;
        }
    }

    return -1; /*  neighbor not found */
}


/****************************************************************************
 Extracted from FillSingleStereoDescriptors(...)
****************************************************************************/
int GetStereoBondParity( sp_ATOM *at,
                         int     i,
                         int     n,
                         AT_RANK *nRank )
{
    int k1, k2, s, parity;

    if (at[i].stereo_bond_neighbor[0])
    {
        for (k1 = 0; k1 < MAX_NUM_STEREO_BONDS && ( s = (int) at[i].stereo_bond_neighbor[k1] ); k1++)
        {
            if (--s == n)
            {
                goto neigh1_found;
            }
        }
        return -1; /*  error: not a stereo neighbor */
    neigh1_found:
        if (PARITY_KNOWN( at[i].stereo_bond_parity[k1] ))
        {
            return PARITY_VAL( at[i].stereo_bond_parity[k1] );
        }
        for (k2 = 0; k2 < MAX_NUM_STEREO_BONDS && ( s = (int) at[n].stereo_bond_neighbor[k2] ); k2++)
        {
            if (--s == i)
            {
                goto neigh2_found;
            }
        }
        return -1; /*  error: not a stereo neighbor */
    neigh2_found:;
    }
    else
    {
        return -1; /*  error: not a stereo bond */
    }

    if (ATOM_PARITY_WELL_DEF( at[i].parity ) &&
         ATOM_PARITY_WELL_DEF( at[n].parity ) &&
         MIN_DOT_PROD <= abs( at[i].stereo_bond_z_prod[k1] ))
    {
        /*  bond parity can be calculated */
        int half_parity1, half_parity2;
        /*  check whether all neighbors are defined */


        half_parity1 = HalfStereoBondParity( at, i, k1, nRank );
        half_parity2 = HalfStereoBondParity( at, n, k2, nRank );
        if (!half_parity1 || !half_parity2)
            return 0; /*  ranks undefined or not a stereo bond */
        if (ATOM_PARITY_WELL_DEF( half_parity1 ) &&
             ATOM_PARITY_WELL_DEF( half_parity2 ))
        {
            parity = 2 - ( half_parity1 + half_parity2
                     + ( at[i].stereo_bond_z_prod[k1] < 0 ) ) % 2;
        }
        else
        {
            return CT_STEREOBOND_ERROR;  /*   <BRKPT> */
        }
    }
    else
    {
        /*  parity cannot be calculated: not enough info or 'unknown' */
        if (AB_PARITY_NONE != ( parity = inchi_max( at[i].parity, at[n].parity ) ))
        {
            parity = AB_PARITY_UNDF; /*  should not happen */
        }
    }

    return parity;
}


/****************************************************************************
 Extracted from FillSingleStereoDescriptors(...)
****************************************************************************/
int GetPermutationParity( CANON_GLOBALS *pCG,
                          sp_ATOM       *at,
                          AT_RANK       nAvoidNeighbor,
                          AT_RANK       *nCanonRank )
{
    AT_RANK nNeighRank[MAX_NUM_STEREO_ATOM_NEIGH];
    int     j, k, parity;
    if (at->valence > MAX_NUM_STEREO_ATOM_NEIGH)
    {
        parity = -1; /*  error */
    }
    else
    {
        for (j = k = 0; j < at->valence; j++)
        {
            if (at->neighbor[j] != nAvoidNeighbor)
            {
                nNeighRank[k++] = nCanonRank[(int) at->neighbor[j]];
            }
        }
        if (k)
        {
            parity = insertions_sort( pCG, nNeighRank, k, sizeof( nNeighRank[0] ), comp_AT_RANK );
            if (nNeighRank[0])
            {
                parity = 2 - parity % 2;
            }
            else
            {
                parity = 0; /*  not all ranks are known */
            }
        }
        else
        {
            /* special case: HX= with implicit H */
            parity = 2;
        }
    }

    return parity;
}


/****************************************************************************/
int GetStereoCenterParity( CANON_GLOBALS *pCG,
                           sp_ATOM       *at,
                           int           i,
                           AT_RANK       *nRank )
{
    AT_NUMB  nNeighborNumber2[MAXVAL];
    int      parity;
    int      k, num_trans;

    if (!at[i].parity)
    {
        return 0;   /*  not a stereo center                     */
    }
    if (at[i].stereo_bond_neighbor[0])
    {
        return -1;  /*  a stereo bond atom, not a stereo center */
    }

    if (ATOM_PARITY_WELL_DEF( at[i].parity ))
    {
        /*  number of neighbors transpositions to the sorted order is unknown. Find it. */
        /*  If parity is not well-defined then doing this is a waste of time            */
        int num_neigh = at[i].valence;
        for (k = 0; k < num_neigh; k++)
        {
            if (!nRank[(int) at[i].neighbor[k]])
                return 0; /*  stereo at[i] does not belong to the traversed part of the structure */
            nNeighborNumber2[k] = k;
        }
        pCG->m_pNeighborsForSort = at[i].neighbor;
        pCG->m_pn_RankForSort = nRank;
        num_trans = insertions_sort( pCG, nNeighborNumber2, num_neigh, sizeof( nNeighborNumber2[0] ), CompNeighborsAT_NUMBER );
#ifndef CT_NEIGH_INCREASE
        num_trans += ( ( num_neigh*( num_neigh - 1 ) ) / 2 ) % 2;  /*  get correct parity for ascending order */
#endif
        parity = 2 - ( at[i].parity + num_trans ) % 2;
    }
    else
    {
        parity = at[i].parity;
    }

    return parity;
}

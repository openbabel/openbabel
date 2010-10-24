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

#include "inpdef.h"
#include "extr_ct.h"
#include "ichitaut.h"
#include "ichinorm.h"
#include "ichierr.h"
#include "util.h"

#include "ichicomp.h"

#ifndef INCHI_ALL_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif

/* defined in ichisort.c, prototype in ichicomn.h */
int insertions_sort_AT_RANK( AT_RANK *base, int num );

#ifndef INCHI_ALL_CPP
#ifdef __cplusplus
}
#endif
#endif

/* local prototypes */
int cmp_iso_atw_diff_component_no( const void *a1, const void *a2 );

/********************************************************************************/
int cmp_iso_atw_diff_component_no( const void *a1, const void *a2 )
{
    int ret = (int)((const inp_ATOM*)a1)->iso_atw_diff - (int)((const inp_ATOM*)a2)->iso_atw_diff;
    if ( !ret ) /*  make the sort stable */
        ret = (int)((const inp_ATOM*)a1)->component - (int)((const inp_ATOM*)a2)->component;
    return ret;
}
typedef struct tagTreeAtom {
    AT_NUMB    neighbor[MAXVAL];        /* positions (from 0) of the neighbors in the inp_ATOM array */
    S_CHAR     valence;                 /* number of bonds = number of neighbors */
    AT_NUMB    nRingSystem;
    AT_NUMB    nBlockSystem;
    S_CHAR     bCutVertex;
} tre_ATOM;

#if( FIND_RING_SYSTEMS == 1 ) /* { */
/********************************************************************************/
int MarkRingSystemsInp( inp_ATOM *at, int num_atoms, int start )
{
    AT_NUMB   *nStackAtom = NULL;
    int        nTopStackAtom=-1;
    AT_NUMB   *nRingStack = NULL;
    int        nTopRingStack=-1; /* was AT_NUMB */
    AT_NUMB   *nDfsNumber = NULL;
    AT_NUMB   *nLowNumber = NULL;
    S_CHAR    *cNeighNumb = NULL;
    AT_NUMB    nDfs;
#if( FIND_RINS_SYSTEMS_DISTANCES == 1 )
    AT_NUMB    nRs, *nRsConnect = NULL;
    int        k;
    AT_NUMB   *tree = NULL;
    int        nNumConnect, nMaxNumConnect, nLenConnect;
#endif
    AT_NUMB    nNumAtInRingSystem;
    int        i, j, u, /*start,*/ nNumRingSystems, nNumStartChildren;

    /*  allocate arrays */
    nStackAtom = (AT_NUMB *)inchi_malloc(num_atoms*sizeof(nStackAtom[0]));
    nRingStack = (AT_NUMB *)inchi_malloc(num_atoms*sizeof(nRingStack[0]));
    nDfsNumber = (AT_NUMB *)inchi_malloc(num_atoms*sizeof(nDfsNumber[0]));
    nLowNumber = (AT_NUMB *)inchi_malloc(num_atoms*sizeof(nLowNumber[0]));
    cNeighNumb = (S_CHAR  *)inchi_malloc(num_atoms*sizeof(cNeighNumb[0]));
#if( FIND_RINS_SYSTEMS_DISTANCES == 1 )
    nRsConnect = (AT_NUMB *)inchi_calloc(3*num_atoms+3,sizeof(nRsConnect[0]));
#endif
    /*  check allocation */
    if ( !nStackAtom || !nRingStack || !nDfsNumber || !nLowNumber || !cNeighNumb 
#if( FIND_RINS_SYSTEMS_DISTANCES == 1 )
        || !nRsConnect
#endif
        ) {
        nNumRingSystems = CT_OUT_OF_RAM;  /*  program error */ /*   <BRKPT> */
        goto exit_function;
    }

    /********************************************
     *
     * Find Cut-vertices & Blocks
     *
     ********************************************/

    /*  initiation */
    /*start           = 0;*/
    nNumRingSystems = 0;
    u               = start; /*  start atom */
    nDfs            = 0;
    nTopStackAtom   =-1;
    nTopRingStack   =-1;
    memset( nDfsNumber, 0, num_atoms*sizeof(nDfsNumber[0]));
    memset( cNeighNumb, 0, num_atoms*sizeof(cNeighNumb[0]));
    /*  push the start atom on the stack */
    nLowNumber[u] = nDfsNumber[u] = ++nDfs;
    nStackAtom[++nTopStackAtom] = (AT_NUMB)u;
    nRingStack[++nTopRingStack] = (AT_NUMB)u;

    nNumStartChildren = 0;

    do {
        /* advance */
advance_block:
        /*if ( (int)at[i=nStackAtom[nTopStackAtom]].valence > (j = (int)cNeighNumb[i]) )*/
        /* replaced due to missing sequence point */
        if ( i=(int)nStackAtom[nTopStackAtom], j = (int)cNeighNumb[i], (int)at[i].valence > j )
        {
            cNeighNumb[i] ++;
            u = (int)at[i].neighbor[j];
            if ( !nDfsNumber[u] ) {
                /* tree edge, 1st visit -- advance */
                nStackAtom[++nTopStackAtom] = (AT_NUMB)u;
                nRingStack[++nTopRingStack] = (AT_NUMB)u;
                nLowNumber[u] = nDfsNumber[u] = ++nDfs;
                nNumStartChildren += (i == start);
            } else
            if ( !nTopStackAtom || u != (int)nStackAtom[nTopStackAtom-1] ) { /*  may comment out ? */
                /* back edge: u is not a predecessor of i */
                if ( nDfsNumber[u] < nDfsNumber[i] ) {
                    /* Back edge, 1st visit: u is an ancestor of i. Compare */
                    if ( nLowNumber[i] > nDfsNumber[u] ) {
                        nLowNumber[i] = nDfsNumber[u];
                    }
                }
            }                                                                /*  may comment out ? */
            goto advance_block;
        } else {
            cNeighNumb[i] = 0;
        }

        /* back up */
        if ( i != start ) {
            u = (int)nStackAtom[nTopStackAtom-1]; /* predecessor of i */
            if ( nLowNumber[i] >= nDfsNumber[u] ) {
                /* output the block */
                nNumRingSystems ++;
                at[u].nBlockSystem = nNumRingSystems;
                if ( u != start || nNumStartChildren > 1 ) {
                    at[u].bCutVertex += 1;
                }
                while ( nTopRingStack >= 0 ) {
                    j = nRingStack[nTopRingStack--];
                    at[j].nBlockSystem = nNumRingSystems; /*  mark the atom */
                    if ( i == j ) {
                        break;
                    }
                }
            } else
            if ( nLowNumber[u] > nLowNumber[i] ) {
                /* inherit */
                nLowNumber[u] = nLowNumber[i];
            }
        }
    } while ( --nTopStackAtom >= 0 );


    /********************************************
     *
     * Find Ring Systems
     * Including chain atoms X: A-X-B, where the bonds (of any kind) are bridges.
     *
     ********************************************/

    /*  initiation */
    /*start           = 0;*/
    nNumRingSystems = 0;
    u               = start; /*  start atom */
    nDfs            = 0;
    nTopStackAtom   =-1;
    nTopRingStack   =-1;
    memset( nDfsNumber, 0, num_atoms*sizeof(nDfsNumber[0]));
    memset( cNeighNumb, 0, num_atoms*sizeof(cNeighNumb[0]));
    /*  push the start atom on the stack */
    nLowNumber[u] = nDfsNumber[u] = ++nDfs;
    nStackAtom[++nTopStackAtom] = (AT_NUMB)u;
    nRingStack[++nTopRingStack] = (AT_NUMB)u;
#if( FIND_RINS_SYSTEMS_DISTANCES == 1 )
    nNumConnect = nLenConnect = nMaxNumConnect = 0;
#endif

    do {
        /* advance */
advance_ring:
        /*if ( (int)at[i=nStackAtom[nTopStackAtom]].valence > (j = (int)cNeighNumb[i]) )*/
        /* replaced due to missing sequence point */
        if ( i=(int)nStackAtom[nTopStackAtom], j = (int)cNeighNumb[i], (int)at[i].valence > j )
        {
            cNeighNumb[i] ++;
            u = (int)at[i].neighbor[j];
            if ( !nDfsNumber[u] ) {
                /* tree edge, 1st visit -- advance */
                nStackAtom[++nTopStackAtom] = (AT_NUMB)u;
                nRingStack[++nTopRingStack] = (AT_NUMB)u;
                nLowNumber[u] = nDfsNumber[u] = ++nDfs;
            } else
            if ( !nTopStackAtom || u != (int)nStackAtom[nTopStackAtom-1] ) {
                /* back edge: u is not a predecessor of i */
                if ( nDfsNumber[u] < nDfsNumber[i] ) {
                    /* Back edge, 1st visit: u is ancestor of i. Compare */
                    if ( nLowNumber[i] > nDfsNumber[u] ) {
                        nLowNumber[i] = nDfsNumber[u];
                    }
                }
            }
            goto advance_ring;
        } else {
            cNeighNumb[i] = 0;
        }

        /* back up */
        if ( nDfsNumber[i] == nLowNumber[i] ) {
            /*  found a ring system */
            nNumRingSystems ++;
            /*  unwind nRingStack[] down to i */
#if( FIND_RINS_SYSTEMS_DISTANCES == 1 )
            nNumConnect = 2;
            /* data structure: for each ring system nRsConnect[] contains:
             * 1) nNumConnect+1 = (number of already discovered neighboring "ring systems" + 1)+1
             * 2) nNumAtInRingSystem
             * 3) (nNumConnect-1) numbers (IDs) of neighboring ring systems.
             * BFS guarantees that each neighboring ring system is encountered only one time
             * Number of all neighboring ring systems = (nNumConnect-1)+1 = nNumConnect
             * (One additional ring system is where the BFS retracts from the vertex #i,
             * except when i=DFS root node. In the latter case there is/are only (nNumConnect-1)
             * neighboring ring system(s).
             */
#endif
            /*  count atoms in a ring system */
            for ( nNumAtInRingSystem = 0, j =  nTopRingStack; 0 <= j; j -- ) {
                nNumAtInRingSystem ++;
                if ( i == (int)nRingStack[j] ) {
                    break;
                }
            }
            while ( nTopRingStack >= 0 ) {
                j = (int)nRingStack[nTopRingStack--];
                at[j].nRingSystem        = (AT_NUMB)nNumRingSystems; /*  ring system id */
                at[j].nNumAtInRingSystem = nNumAtInRingSystem;
#if( FIND_RINS_SYSTEMS_DISTANCES == 1 )
                for ( k = 0; k < at[j].valence; k ++ ) {
                    if ( (nRs = at[at[j].neighbor[k]].nRingSystem) && (int)nRs != nNumRingSystems ) {
                        nRsConnect[nLenConnect + (nNumConnect++)] = nRs; /*  adjacent ring system id */
                    }
                }
#endif
                if ( i == j ) {
                    /*  reached atom on the top of nStackAtom[] stack  */
                    break;
                }
            }
#if( FIND_RINS_SYSTEMS_DISTANCES == 1 )
            nRsConnect[nLenConnect] = nNumConnect;
            nRsConnect[nLenConnect+1] = nNumAtInRingSystem;
            nLenConnect += nNumConnect;
            if ( nMaxNumConnect < nNumConnect ) {
                /*  max number of neighboring ring systems */
                nMaxNumConnect = nNumConnect;
            }
#endif
        } else
        if ( nTopStackAtom > 0 ) {
            j = (int)nStackAtom[nTopStackAtom-1];
            /* inherit nLowNumber */
            if ( nLowNumber[j] > nLowNumber[i] ) {
                nLowNumber[j] = nLowNumber[i];
            }
        }
    } while ( --nTopStackAtom >= 0 );

#if( FIND_RINS_SYSTEMS_DISTANCES == 1 ) /*  normally disabled */
    nMaxNumConnect ++;
    if ( nNumRingSystems > 1 ) {
        int nCol      = nMaxNumConnect+1;
        int nNumInSyst= nMaxNumConnect;
        int nMaxNeigh = nMaxNumConnect-1;
#define T(a,b) tree[(a)*nCol+b]
        if ( tree = (AT_NUMB *)inchi_calloc( nCol * (nNumRingSystems+1), sizeof(tree[0])) ) {
            int len, neigh;
            /*  reuse previous allocations */
            AT_NUMB *nNumVisitedNeighbors  = nStackAtom;
            AT_NUMB *nDistanceFromTerminal = nRingStack;
            AT_NUMB *nCurrActiveRingSystem = nDfsNumber;
            AT_NUMB *nNextActiveRingSystem = nLowNumber;
            int        nNumCurrActiveRingSystems, nNumNextActiveRingSystems, pass;
            /* build a "condensation graph (actually, a tree)" in which
             * each vertex corresponds to a ring system T(row, col) = T(ring syst, neighbors)
             * Number of rows = column length = max. number of ring system neighbors + 2
             * Number of cols = row length    = number of ring systems + 1
             * Neighboring ring systems are contiguously stored in a row
             * T(i,0) = number of neighbors,  1 <= i <= nNumRingSystems;
             * T(i,k) = number of a neighboring ring system, 1 <= k <= T(i,0)
             * T(i,nCol-1) = number of atoms in the system #i
             */
            for ( i = 1, j = 0; len=nRsConnect[j]; i ++ ) {
                T(i, nNumInSyst) = nRsConnect[j+1];
                for ( k = 2; k < len; k ++ ) {
                    neigh = nRsConnect[j+k];
                    if ( T(i,0) < nMaxNeigh && T(neigh,0) < nMaxNeigh ) {
                        T(i,0) ++;
                        T(neigh,0) ++;
                        T(i,T(i,0))         = neigh;
                        T(neigh,T(neigh,0)) = i;
                    } else {
                        nNumRingSystems = CT_OVERFLOW;  /*  program error */ /*   <BRKPT> */
                        goto exit_function;
                    }
                }
                j += len;
            }
            /*  clear memory */
            memset( nNumVisitedNeighbors,  0, nNumRingSystems*sizeof(nNumVisitedNeighbors[0]) );
            memset( nDistanceFromTerminal, 0, nNumRingSystems*sizeof(nDistanceFromTerminal[0]) );
            memset( nCurrActiveRingSystem, 0, nNumRingSystems*sizeof(nCurrActiveRingSystem[0]) );
            memset( nNextActiveRingSystem, 0, nNumRingSystems*sizeof(nNextActiveRingSystem[0]) );
            nNumNextActiveRingSystems = 0;
            for ( i = 0; i < nNumRingSystems; i ++ ) {
                if ( 1 == T(i+1,0) ) {
                    nNextActiveRingSystem[i] = 1; /*  number of traversed neighbors + 1 */
                    nDistanceFromTerminal[i] = 1;
                    nNumNextActiveRingSystems ++;
                } else {
                    nNextActiveRingSystem[i] = 0;
                    nDistanceFromTerminal[i] = 0;
                }
                nNumVisitedNeighbors[i]  = 0;
            }

            /* nCurrActiveRingSystem[i] = a sum of:
             * 1) +1 if it is or was active
             * 2) +(number of neighbors from which it was reached)
             * 3) +1 if it was left and not active anymore
             */
            pass = 0;
            do {
                nNumCurrActiveRingSystems = nNumNextActiveRingSystems;
                nNumNextActiveRingSystems = 0;
                memcpy( nCurrActiveRingSystem, nNextActiveRingSystem,
                        nNumRingSystems*sizeof(nNextActiveRingSystem[0]));
                for ( i = 0; i < nNumRingSystems; i ++ ) {
                    if ( T(i+1,0) == nCurrActiveRingSystem[i] ) {
                        /* on the previous pass currently active ring system i+1 bas been reached
                         * from all neighbors except one;
                         * the neighbors from which it was reached have
                         * T(neigh,0)+1 == nCurrActiveRingSystem[i]
                         * this ring system has not been left yet
                         */
                        for ( k = 1, len=T(i+1,0); k <= len; k ++ ) {
                            neigh = (int)T(i+1,k);
                            if ( T(neigh,0) >= nCurrActiveRingSystem[neigh-1] ) {
                                if ( 0 == pass ) {
                                    nDistanceFromTerminal[i] = 1;
                                }
                                break;
                            }
                        }
                        if ( k <= len ) {
                            /* neigh was not reached from at least 2 neighbors
                             * walk along -R- chain (T(neigh,0)==2) up to
                             * 1)  a terminal system, not including it or
                             * 2)  a branching point.
                             *
                             * pass = 0: started from terminal systems:
                             *     reach the branching point.
                             * If chain system next to a terminal system has already been reached
                             * then walk along it according to Note below
                             *
                             * pass > 0: started from branching points
                             * 2a) If the branching point has not been reached from 2 or more neighbors,
                             *     then include it
                             * 2b) If the branching point has not been reached from 1 neighbor only,
                             *     then do not include it: it will be a starting point later
                             * Note: if a chain atom already has nDistanceFromTerminal[i] > 0, then
                             *     the last atom should be the one such that
                             *     its nDistanceFromTerminal[]+1>= nDistanceFromTerminal[] of the
                             *     next in the chain
                             */
                            int bOk = 0;
                            k = i+1; /*  starting point */
                            if ( 0 == pass && T(k,nNumInSyst) > 1 ) {
                                nNumNextActiveRingSystems ++; /*  request next pass */
                                continue; /*  stop a the terminal ring system */
                            }
                            while( 2 == T(neigh,0) ) {
                                /*  walk along a chain */
                                if ( !nNextActiveRingSystem[neigh-1] ) {
                                    nNextActiveRingSystem[neigh-1] = 1; /*  make neighbor active */
                                } else
                                if ( nDistanceFromTerminal[k-1]+1 <= nDistanceFromTerminal[neigh-1] ) {
                                    /*  walking along the chain; already have had a walk */
                                    /*  in the opposite direction at this pass */
                                } else {
                                     /*  k is the last; neigh (it is a bridge -X-) has not been reached */
                                    bOk = 1;
                                    break;
                                }
                                nNextActiveRingSystem[k-1] ++; /*  leave system k */
                                if ( nNextActiveRingSystem[neigh-1] < T(neigh,0) ) {
                                    nNextActiveRingSystem[neigh-1] ++; /*  add one connection to neigh */
                                }
                                nDistanceFromTerminal[neigh-1] = nDistanceFromTerminal[k-1]+1;
                                j = (T(neigh,1)==k)? 2:1;
                                k = neigh;
                                neigh = T(k,j); /*  next in the chain */
                                nNumNextActiveRingSystems ++;
                                if ( T(k,nNumInSyst) > 1 ) {
                                    bOk = 1;
                                    break; /*  stop on a ring system */
                                }
                            }
                            /*  neigh is a terminal or a bridge or a branching point */
                            if ( 2 > T(neigh,0) ) {
                                /*  neighbor is a terminal atom */
                                if ( 1 < pass ) {
                                    nNumRingSystems = CT_UNKNOWN_ERR; /*  error (debug only) */ /*   <BRKPT> */
                                    goto exit_function;
                                }
                                continue;
                            }
                            if ( 2 == T(neigh,0) ) {
                                /*  neighbor is a bridge */
                                continue;
                            }
                            /*  neighbor is a branching point */
                            if ( T(neigh,0) > nCurrActiveRingSystem[neigh-1] ) {
                                /*  move to the neigh (make neigh active): on previous pass it */
                                /*  has not been reached from 2 or more neighbors */
                                if ( !nNextActiveRingSystem[neigh-1] ) {
                                    nNextActiveRingSystem[neigh-1] = 1;
                                }
                                if ( nDistanceFromTerminal[neigh-1] < nDistanceFromTerminal[k-1]+1 ) {
                                    nDistanceFromTerminal[neigh-1] = nDistanceFromTerminal[k-1]+1;
                                }
                                nNextActiveRingSystem[k-1] ++; /*  leave system k */
                                if ( nNextActiveRingSystem[neigh-1] < T(neigh,0) ) {
                                    nNextActiveRingSystem[neigh-1] ++; /*  add one connection to neigh */
                                }
                                nNumNextActiveRingSystems ++;
                            }
                        }
                    }
                }
                pass ++;
            } while ( nNumNextActiveRingSystems );

            for ( i = 0; i < num_atoms; i ++ ) {
                at[i].nDistanceFromTerminal = nDistanceFromTerminal[(int)at[i].nRingSystem-1];
            }

            inchi_free( tree );
            tree = NULL;
#undef T
        } else {
            nNumRingSystems = CT_OUT_OF_RAM; /*  error */ /*   <BRKPT> */
            goto exit_function;
        }
    }
#endif


exit_function:
    if ( nStackAtom )
        inchi_free( nStackAtom );
    if ( nRingStack )
        inchi_free( nRingStack );
    if ( nDfsNumber )
        inchi_free( nDfsNumber );
    if ( nLowNumber )
        inchi_free( nLowNumber );
    if ( cNeighNumb )
        inchi_free( cNeighNumb );
#if( FIND_RINS_SYSTEMS_DISTANCES == 1 )
    if ( nRsConnect )
        inchi_free( nRsConnect );
    if ( tree )
        inchi_free( tree );
#endif
    return nNumRingSystems;
}


#endif /* } FIND_RING_SYSTEMS */

/********************************************************************************/
/*  Return value: new number of atoms > 0 or -1=out of RAM */
int remove_terminal_HDT( int num_atoms, inp_ATOM *at, int bFixTermHChrg )
{
    AT_NUMB   *new_ord;
    inp_ATOM  *new_at;
    char *p;
    const static char szHDT[]="HDT";
    const static int  kMax = sizeof(szHDT); /*  = 4 */
    int ret = -1;
    int num_hydrogens=0, num_H = 0;  /*  number of terminal H, D, T */
    int i, j, k, n, m;
    int val;
    AT_RANK new_HydrogenAt_order[NUM_H_ISOTOPES+1];
    AT_RANK new_OtherNeigh_order[MAXVAL];
    S_CHAR  old_trans[MAX_NUM_STEREO_BONDS];
    
    int  num_OtherNeigh, num_HydrogenAt;

    new_ord=(AT_NUMB *)inchi_calloc(num_atoms, sizeof(new_ord[0])); /* changed malloc to calloc 9-11-2003 */
    new_at =(inp_ATOM  *)inchi_malloc(sizeof(new_at[0]) *num_atoms);
    if (!new_ord || !new_at)
        goto exit_function;

    /*  move H. D, T to the end of the list of atoms */
    for ( i = 0; i < num_atoms; i ++ ) 
    {
        at[i].component = i; /*  temporarily save original numbering */
        /*  get k = temp. hydrogen isotope/non-hydrogen atom type: */
        /*  k=0:H, k=2:D, k=3:T, k=4=kMax: not a hydrogen */
        k = at[i].elname[1]? kMax : (p=(char*)strchr(szHDT, at[i].elname[0]))? p-szHDT : kMax;
        /*  set hydrogen isotope atw differences */
        /*  Notes: k-value of isotopic H is incremented to correct iso_atw_diff value later. */
        /*         1H isotope cannot be detected here. */
        if ( k == ATW_H || k == ATW_H+1 ) 
        {  
            /* D or T, k = 1 or 2 */
            at[i].elname[0]    = 'H'; /*  hydrogen isotope */
            at[i].iso_atw_diff = ++k; /*  increment k to make k = iso_atw_diff ( 2 for D, 3 for T ) */
        }
        num_H += (k != kMax && at[i].valence == 1 && at[i].chem_bonds_valence == 1 && !NUMH(at,i) );
    }
    
    /* special case: HD, HT, DT, HH: the only non-isotopic H or
     * the lightest isotopic H out of two is removed
     * to become implicit (make the heavier H the "central atom"). 
     * Note: This must be consistent with mol_to_atom()
     * treatment of isotopic Hn aliases.
     */
    if ( 2 == num_H && 2 == num_atoms && !NUMH(at,0) && !NUMH(at,1) ) 
    {

        if ( at[0].iso_atw_diff >= at[1].iso_atw_diff ) {
            new_ord[0] = 0;
            new_ord[1] = 1;
        } else {
            new_ord[0] = 1;
            new_ord[1] = 0;
        }
        if ( at[new_ord[1]].charge ) {
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
        for ( i = 0; i < num_atoms; i ++ ) 
        {
            k = (at[i].elname[1] || NUMH(at,i))? kMax : (at[i].elname[0]=='H')? at[i].iso_atw_diff : kMax;
            if ( k < kMax && at[i].valence == 1 && at[i].chem_bonds_valence == 1 &&
                 /*  the order of comparison is important */
                 ((n=(int)at[i].neighbor[0]) > i               /* at[n] has not been encountered yet*/ ||
                  (int)new_ord[n] < num_atoms - num_hydrogens) /* at[n] might have been encountered; it has not been moved */ ) 
            {
                /*  found an explicit terminal hydrogen */
                num_hydrogens ++;
                if ( k==0 && ATW_H <= at[i].iso_atw_diff && at[i].iso_atw_diff < ATW_H+NUM_H_ISOTOPES ) 
                {
                    k = at[i].iso_atw_diff; /*  H isotope has already been marked above or elsewhere */
                }
                if ( at[i].charge ) 
                { 
                    /*  transfer charge from the hydrogen */
                    at[n].charge += at[i].charge;
                    at[i].charge = 0;
                    if (bFixTermHChrg)
                    {
                        /*^^^^^ Fixed bug (July 6, 2008 IPl) : 
                                if terminal H was charged (not neutralized before call of remove_terminal_HDT) 
                                and had an ordering number > than that of heavy-atom neighbour, then 
                                charge on neighbour atom was not adjusted (though charge on H was removed).
                        ^^^^^ */
                        if ( i > n )
                            /* new_at[new_ord[n]] has been created and filled already */
                            new_at[new_ord[n]].charge = at[n].charge;
                    }
                    /*^^^^^ */
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

    if ( num_hydrogens ) {
        int num_others = num_atoms-num_hydrogens; /*  atoms which are not terminal H, D, T */
        if ( num_hydrogens > 1 ) {
            /*  sort hydrogen isotopes in ascending order, */
            /*  orig, numbers being the secondary sorting key */
            qsort( new_at+num_others, num_hydrogens, sizeof(new_at[0]), cmp_iso_atw_diff_component_no );
        }
        /*  save new numbering of hydrogen atoms using temporarily saved orig numbering */
        for ( i = num_others; i < num_atoms; i ++ ) {
            new_ord[(int)new_at[i].component] = i;
        }

        /*  renumber neighbors according to new_ord[] and detach terminal hydrogens */
        for ( i = 0; i < num_others; i++ ) {
            memset( new_HydrogenAt_order, 0, sizeof(new_HydrogenAt_order) );
            memset( new_OtherNeigh_order, 0, sizeof(new_OtherNeigh_order) );
            num_OtherNeigh = 0;
            num_HydrogenAt = 0;
            num_H          = 0;
            
            for ( m = 0; m < MAX_NUM_STEREO_BONDS && new_at[i].sb_parity[m]; m ++ ) {
                old_trans[m] = 2 - (new_at[i].sn_ord[m] + new_at[i].sb_ord[m] + (new_at[i].sn_ord[m] > new_at[i].sb_ord[m]))%2;
            }
                
            for ( k = j = val = 0; k < new_at[i].valence; k++ ) {
                if ( num_others <= ( n = new_ord[new_at[i].neighbor[k]] ) ) {
                    /*  discovered neighbor = disconnected explicit hydrogen
                     *  i = new atom new_at[i] ordering number
                     *  n = new number of the explicit H
                     *  k = ordering number of the explicit H in new_at[i] adjacency list
                     */
                    if ( 0 < new_at[n].iso_atw_diff && new_at[n].iso_atw_diff < ATW_H+NUM_H_ISOTOPES ) {
                        /* make explicit isotopic H implicit */
                        new_at[i].num_iso_H[new_at[n].iso_atw_diff-1] ++; /*  isotopic H */
                        num_HydrogenAt += !new_HydrogenAt_order[new_at[n].iso_atw_diff];
                        new_HydrogenAt_order[new_at[n].iso_atw_diff] = k+1;
                    } else {
                        /* make explicit non-isotopic H implicit */
                        new_at[i].num_H ++; /*  non-isotopic H */
                        num_HydrogenAt += !num_H;
                        num_H ++;
                        new_HydrogenAt_order[0] = k+1;
                    }
                    /*  decrement chem. bonds valence because one bond is removed */
                    new_at[i].chem_bonds_valence  = inchi_max( 0, new_at[i].chem_bonds_valence-1 );
                    new_at[n].neighbor[0]         = i; /*  update removed hydrogen neighbor number */
                    if ( new_at[i].sb_parity[0] ) {
                        /* if the removed H is an SB neighbor then mark it as removed */
                        for ( m = 0; m < MAX_NUM_STEREO_BONDS && new_at[i].sb_parity[m]; m ++ ) {
                            if ( k == (int)new_at[i].sn_ord[m] ) {
                                new_at[i].sn_ord[m] = -(new_at[n].iso_atw_diff+1); 
                                /* means the SB neighbor has been removed; (-4)=H, (-3)=1H, (-2)=D, (-1)=T */
                            }
                        }
                    }
                } else {
                    /* discovered a regular (not an explicit H) neighbor */
                    if ( new_at[i].sb_parity[0] ) {
                        if ( num_OtherNeigh < MAX_NUM_STEREO_BONDS ) {
                            new_OtherNeigh_order[num_OtherNeigh] = k+1;
                        }
                        num_OtherNeigh ++; /* increment outside of if() to detect overflow */
                        if ( val != k ) {
                            /* store new stereobond and sb-neighbor ordering numbers */
                            for ( m = 0; m < MAX_NUM_STEREO_BONDS && new_at[i].sb_parity[m]; m ++ ) {
                                if ( k == (int)new_at[i].sb_ord[m] )
                                    new_at[i].sb_ord[m] = val;
                                else
                                if ( k == (int)new_at[i].sn_ord[m] )
                                    new_at[i].sn_ord[m] = val;
                            }
                        }
                    }
                    new_at[i].neighbor[val]       = new_ord[new_at[i].neighbor[k]];
                    new_at[i].bond_type[val]      = new_at[i].bond_type[k];
                    new_at[i].bond_stereo[val]    = new_at[i].bond_stereo[k];
                    val ++;
                }
            }
            if ( new_at[i].valence > val && new_at[i].sb_parity[0] ) {
                if ( num_HydrogenAt == new_at[i].valence - val && num_HydrogenAt + num_OtherNeigh <= MAXVAL ) {
                    /* recalculate parity so that it would describe neighbor sequence H,1H,D,T,neigh[0],neigh[1]... */
                    memmove( new_OtherNeigh_order + num_HydrogenAt, new_OtherNeigh_order, num_OtherNeigh*sizeof(new_OtherNeigh_order[0])); 
                    for ( k = 0, j = 1; k <= NUM_H_ISOTOPES; k ++ ) {
                        if ( new_HydrogenAt_order[k] ) {
                            new_OtherNeigh_order[num_HydrogenAt - j] = new_HydrogenAt_order[k];
                            for ( m = 0; m < MAX_NUM_STEREO_BONDS && new_at[i].sb_parity[m]; m ++ ) {
                                if ( (int)new_at[i].sn_ord[m] == -(k+1) ) {
                                    new_at[i].sn_ord[m] = -j; 
                                    /* negative means explicit H isotope ord are
                                       (contiguously) in front of the adjacency list */
                                }
                            }
                            j ++;
                        }
                    }
                    /* at this point new_OtherNeigh_order[] contains
                       incremented old ordering numbers in new order */
                    k = insertions_sort_AT_RANK( new_OtherNeigh_order, num_HydrogenAt + num_OtherNeigh );
                    k = k%2; /* seems to be of no use */
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
#ifdef _DEBUG
                else {
                    /* error */
                    int stop = 1;
                }
#endif
            }
            new_at[i].valence = val;
        }
        memcpy( at, new_at, sizeof(at[0])*num_atoms );
        ret = num_others;
    }  else {
        ret = num_atoms;
    }
exit_function:
    if ( new_ord )
        inchi_free ( new_ord );
    if ( new_at )
        inchi_free ( new_at );
    return ret;
}
/************************************************************************/
int add_DT_to_num_H( int num_atoms, inp_ATOM *at )
/*  assume num_1H, num_D and num_T are not included in num_H */
{
    int i, j;
    for ( i = 0; i < num_atoms; i ++ ) {
        for ( j = 0; j < NUM_H_ISOTOPES; j ++ )
            at[i].num_H += at[i].num_iso_H[j];
    }
    return 0;
}
/***************************************************************/
/* not used ---
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
*/
/********************************************************************

             InChI post-version 1.01 features implementation
             (v. 1.03 : underivatize is still experimental and for engineering mode)

 ********************************************************************/
#if( RING2CHAIN == 1 || UNDERIVATIZE == 1 )

static U_CHAR el_number_O;
static U_CHAR el_number_C;
static U_CHAR el_number_N;
static U_CHAR el_number_P;
static U_CHAR el_number_S;
static U_CHAR el_number_Si;
static U_CHAR el_number_F;
static U_CHAR el_number_Cl;
static U_CHAR el_number_Br;
static U_CHAR el_number_I;
static U_CHAR el_number_B;

typedef struct tagAtPair {
    AT_NUMB at[2];  /* at[0] < at[1] */
} R2C_ATPAIR;

int DisconnectInpAtBond( inp_ATOM *at, AT_NUMB *nOldCompNumber, int iat, int neigh_ord );
int ExtractConnectedComponent(  inp_ATOM *at, int num_at, int component_number, inp_ATOM *component_at );
int mark_arom_bonds( inp_ATOM *at, int num_atoms );

void set_R2C_el_numbers( void );
int UnMarkDisconnectedComponents( ORIG_ATOM_DATA *orig_inp_data );
int UnMarkRingSystemsInp( inp_ATOM *at, int num_atoms );
int UnMarkOtherIndicators( inp_ATOM *at, int num_atoms );
int UnMarkOneComponent( inp_ATOM *at, int num_atoms );
int subtract_DT_from_num_H( int num_atoms, inp_ATOM *at );
int add_inp_ATOM( inp_ATOM *at, int len_at, int len_cur, inp_ATOM *add, int len_add );
int cmp_r2c_atpair( const void *p1, const void *p2 );
int has_atom_pair( R2C_ATPAIR *ap, int num_ap, AT_NUMB at1, AT_NUMB at2 );
int mark_atoms_ap( inp_ATOM *at, AT_NUMB start, R2C_ATPAIR *ap, int num_ap, int num, AT_NUMB cFlags );

/********************************************************************/
void set_R2C_el_numbers( void )
{
    if ( !el_number_O ) {
        el_number_O  = (U_CHAR)get_periodic_table_number( "O" );
        el_number_C  = (U_CHAR)get_periodic_table_number( "C" );
        el_number_N  = (U_CHAR)get_periodic_table_number( "N" );
        el_number_P  = (U_CHAR)get_periodic_table_number( "P" );
        el_number_S  = (U_CHAR)get_periodic_table_number( "S" );
        el_number_Si = (U_CHAR)get_periodic_table_number( "Si" );
        el_number_F  = (U_CHAR)get_periodic_table_number( "F" );
        el_number_Cl = (U_CHAR)get_periodic_table_number( "Cl" );
        el_number_Br = (U_CHAR)get_periodic_table_number( "Br" );
        el_number_I  = (U_CHAR)get_periodic_table_number( "I" );
        el_number_B  = (U_CHAR)get_periodic_table_number( "B" );
    }
}
/***************************************************************/
int UnMarkDisconnectedComponents( ORIG_ATOM_DATA *orig_inp_data )
{
    int i;
    for ( i = 0; i < orig_inp_data->num_inp_atoms; i ++ ) {
        orig_inp_data->at[i].orig_compt_at_numb = 0;
        orig_inp_data->at[i].component          = 0;
    }
    if ( orig_inp_data->nCurAtLen ) {
        inchi_free( orig_inp_data->nCurAtLen );
        orig_inp_data->nCurAtLen = NULL;
    }
    if ( orig_inp_data->nOldCompNumber ) {
        inchi_free( orig_inp_data->nOldCompNumber );
        orig_inp_data->nOldCompNumber = NULL;
    }
    orig_inp_data->num_components = 0;
    return 0;
}
/***************************************************************/
int UnMarkRingSystemsInp( inp_ATOM *at, int num_atoms )
{
    int i;
    for ( i = 0; i < num_atoms; i ++ ) {
        at[i].bCutVertex         = 0;
        at[i].nRingSystem        = 0;
        at[i].nNumAtInRingSystem = 0;
        at[i].nBlockSystem       = 0;
    }
    return 0;
}
/***************************************************************/
int UnMarkOtherIndicators( inp_ATOM *at, int num_atoms )
{
    int i;
    for ( i = 0; i < num_atoms; i ++ ) {
        at[i].at_type         = 0;
        at[i].cFlags          = 0;
    }
    return 0;
}
/***************************************************************/
int UnMarkOneComponent( inp_ATOM *at, int num_atoms )
{
    int i;
    for ( i = 0; i < num_atoms; i ++ ) {
        at[i].orig_compt_at_numb = 0;
        at[i].component          = 0;
    }
    return 0;
}
/***************************************************************/
int subtract_DT_from_num_H( int num_atoms, inp_ATOM *at )
/*  assume num_1H, num_D and num_T are included in num_H */
{
    int i, j;
    for ( i = 0; i < num_atoms; i ++ ) {
        for ( j = 0; j < NUM_H_ISOTOPES; j ++ )
            at[i].num_H -= at[i].num_iso_H[j];
    }
    return 0;
} 
/***************************************************************/
int add_inp_ATOM( inp_ATOM *at, int len_at, int len_cur, inp_ATOM *add, int len_add )
{
    int i, j;
    inp_ATOM *a;
    /* chack correctness */
    if ( len_cur < 0 )
        return len_cur;
    if ( len_add < 0 )
        return len_add;
    if ( len_cur + len_add > len_at )
        return -1;
    /* copy */
    memcpy( at+len_cur, add, len_add*sizeof(at[0]) );
    /* modify */
    if ( len_cur ) {
        a = at + len_cur;
        for ( i = 0; i < len_add; i ++ ) {
            for ( j = 0; j < a[i].valence; j ++ ) {
                a[i].neighbor[j] += len_cur;
            }
        }
    }
    return len_cur + len_add;
}
/****************************************************************************/
int mark_arom_bonds( inp_ATOM *at, int num_atoms )
{
    INCHI_MODE bTautFlags=0, bTautFlagsDone = 0;
    inp_ATOM *at_fixed_bonds_out = NULL;
    T_GROUP_INFO *t_group_info   = NULL;
    int ret;
    ret = mark_alt_bonds_and_taut_groups ( at, at_fixed_bonds_out, num_atoms,
                                           t_group_info, &bTautFlags, &bTautFlagsDone );
    return ret;

}

/********************************************************************/
int cmp_r2c_atpair( const void *p1, const void *p2 )
{
    const R2C_ATPAIR *ap1 = (const R2C_ATPAIR *)p1;
    const R2C_ATPAIR *ap2 = (const R2C_ATPAIR *)p2;
    int diff = (int)ap1->at[0] - (int)ap2->at[0];
    if ( !diff ) {
        diff = (int)ap1->at[1] - (int)ap2->at[1];
    }
    return diff;
}
/***************************************************************/
int has_atom_pair( R2C_ATPAIR *ap, int num_ap, AT_NUMB at1, AT_NUMB at2 )
{
    R2C_ATPAIR ap1;
    int i1, i2, i3, diff;
    int n = at1 > at2;

    ap1.at[n]   = at1;
    ap1.at[1-n] = at2;
    i1 = 0;
    i2 = num_ap-1;
    /* search for ap1 by simple bisections */
    do {
        i3 = (i1 + i2)/2;
        if ( !(diff = cmp_r2c_atpair(&ap1, ap+i3) ) ) {
            return i3+1;  /* found => positive number */
        } else
        if ( diff > 0 ) {
            i1 = i3 + 1;
        } else {
            i2 = i3 - 1;
        }
    }while ( i2 >= i1 );
    return 0; /* not found */
}

/***************************************************************/
/* DFS search for atoms that do not have a flag */ 
int mark_atoms_ap( inp_ATOM *at, AT_NUMB start, R2C_ATPAIR *ap, int num_ap, int num, AT_NUMB cFlags )
{
    if ( !at[start].at_type ) {
        int i;
        AT_NUMB neigh;
        at[start].at_type = cFlags;
        num ++;
        for ( i = 0; i < at[start].valence; i ++ ) {
            neigh = at[start].neighbor[i];
            if ( has_atom_pair( ap, num_ap, start, neigh ) )
                continue;
            num = mark_atoms_ap( at, neigh, ap, num_ap, num, cFlags );
        }
    }
    return num; /* number of atoms traversed forward from at[start] */
}

#endif /* RING2CHAIN || UNDERIVATIZE */

#if( UNDERIVATIZE == 1 )
/***************************************************************/

#ifdef NEVER
typedef struct tagAtTypeBitmap {
    AT_NUMB ord1 : 5; /* up to 2^5-1 = 31 = 0x0037 */
    AT_NUMB ord2 : 5;
    AT_NUMB type : 6; /* up to 2^6-1 = 63 = 0x0077 */
} AtTypeBitmap;
typedef union tagAtTypeUnion {
    AT_NUMB num;
    AtTypeBitmap bits;
} AtTypeUnion;
#endif

#define DERIV_AT_LEN  4
typedef struct tagDerivAttachment {
    char typ[DERIV_AT_LEN];
    char ord[DERIV_AT_LEN];
    char num[DERIV_AT_LEN];
} DERIV_AT;


#define DERIV_BRIDGE_O  0x0001   /* R1-O-R2 => R1-OH + HO-R2 */
#define DERIV_BRIDGE_NH 0x0002   /* R1-NH-R2  amine */ 
#define DERIV_AMINE_tN  0x0004   /* R1-N(-R2)-R3  tertiary amine */ 
#define DERIV_RING_O    0x0008   /* -O- in a ring */
#define DERIV_RING_NH   0x0010   /* -NH- in a ring */
#define DERIV_UNMARK    0x0040   /* unmark the cut */
#define DERIV_NOT       0x1000   /* cannot be a derivatization agent atom */

#define DERIV_DUPLIC    0x0080   /* duplicated disconnection */

#define DERIV_RING     (DERIV_RING_O | DERIV_RING_NH)

#define MAX_AT_DERIV      12
#define NOT_AT_DERIV      99
#define MIN_AT_LEFT_DERIV 3
#define NO_ORD_VAL        0x0037

#define CFLAG_MARK_BRANCH      1          /* for main derivative traversal */
#define CFLAG_MARK_BLOCK       2          /* for block detection */
#define CFLAG_MARK_BLOCK_INV   ((char)~(CFLAG_MARK_BLOCK)) /* for block detection */
#define COUNT_ALL_NOT_DERIV    1      /* 1=> count ALL atoms that are not in deriv. agents */
                                      /* 0=> only atoms that are not in DERIV_RING */

int mark_atoms_cFlags( inp_ATOM *at, int start, int num, char cFlags );
int un_mark_atoms_cFlags( inp_ATOM *at, int start, int num, char cFlags, char cInvFlags );
int is_C_or_S_DB_O( inp_ATOM *at, int i );
int is_C_unsat_not_arom( inp_ATOM *at, int i );
int is_C_Alk( inp_ATOM *at, int i, char cFlags );
int is_Si_IV( inp_ATOM *at, int i );
int is_P_TB_N( inp_ATOM *at, int i );
int is_possibly_deriv_neigh( inp_ATOM *at, int iat, int iord, int type, char cFlags );
int get_traversed_deriv_type( inp_ATOM *at, DERIV_AT *da, int k, DERIV_AT *da1, char cFlags );
int add_to_da( DERIV_AT *da, DERIV_AT *add );
int mark_atoms_deriv( inp_ATOM *at, DERIV_AT *da, int start, int num, char cFlags, int *pbFound );
int count_one_bond_atoms( inp_ATOM *at, DERIV_AT *da, int start, int ord, char cFlags, int *bFound );
int is_silyl( inp_ATOM *at, int start, int ord_prev );
int is_Me_or_Et( inp_ATOM *at, int start, int ord_prev );
int is_CF3_or_linC3F7( inp_ATOM *at, int start, int ord_prev );
int is_phenyl( inp_ATOM *at, int start, int ord_prev );
int is_deriv_ring( inp_ATOM *at, int start, int num_atoms, DERIV_AT *da1, int idrv );
int is_deriv_chain( inp_ATOM *at, int start, int num_atoms, DERIV_AT *da1, int idrv );
int is_deriv_chain_or_ring( inp_ATOM *at, int start, int num_atoms, DERIV_AT *da1, int *idrv );
int remove_deriv( DERIV_AT *da1, int idrv );
int remove_deriv_mark( DERIV_AT *da1, int idrv );
int EliminateDerivNotInList( inp_ATOM *at, DERIV_AT *da, int num_atoms );
int make_single_cut( inp_ATOM *at, DERIV_AT *da, int iat, int icut );

/***************************************************************/
/* DFS search for atoms that do not have a flag */ 
int mark_atoms_cFlags( inp_ATOM *at, int start, int num, char cFlags )
{
    if ( !(at[start].cFlags & cFlags) ) {
        int i;
        at[start].cFlags |= cFlags;
        num ++;
        for ( i = 0; i < at[start].valence; i ++ ) {
            num = mark_atoms_cFlags( at, at[start].neighbor[i], num, cFlags );
        }
    }
    return num; /* number of atoms traversed forward from at[start] */
}
/***************************************************************/
/* DFS search for atoms that do have a flag */ 
int un_mark_atoms_cFlags( inp_ATOM *at, int start, int num, char cFlags, char cInvFlags )
{
    if ( at[start].cFlags & cFlags ) {
        int i;
        at[start].cFlags &= cInvFlags;
        num ++;
        for ( i = 0; i < at[start].valence; i ++ ) {
            num = un_mark_atoms_cFlags( at, at[start].neighbor[i], num, cFlags, cInvFlags );
        }
    }
    return num; /* number of atoms traversed forward from at[start] */
}
/***************************************************************/
int is_C_or_S_DB_O( inp_ATOM *at, int i )
{
    int j, neigh;
    if ( at[i].el_number != el_number_C &&
         at[i].el_number != el_number_S ||
         at[i].charge || at[i].radical )
        return 0;
    for ( j = 0; j < at[i].valence; j ++ ) {
        neigh = at[i].neighbor[j];
        if ( (at[neigh].el_number == el_number_O ||
              at[neigh].el_number == el_number_S ) &&
              !at[neigh].num_H && 1 == at[neigh].valence &&
              2 == at[neigh].chem_bonds_valence ) {
            return 1;
        }
    }
    return 0;
}
/***************************************************************/
int is_C_unsat_not_arom( inp_ATOM *at, int i )
{
    int j, neigh, num_arom, num_DB;
    if ( at[i].el_number != el_number_C ||
         at[i].valence   == at[i].chem_bonds_valence ||
         at[i].valence+1 < at[i].chem_bonds_valence  ||
         at[i].chem_bonds_valence + at[i].num_H != 4 ||
         at[i].charge || at[i].radical )
        return 0;
    num_arom = num_DB = 0;
    for ( j = 0; j < at[i].valence; j ++ ) {
        neigh = at[i].neighbor[j];
        num_arom += at[i].bond_type[j] == BOND_TYPE_ALTERN;
        if ( (at[neigh].el_number == el_number_O ||
              at[neigh].el_number == el_number_S ) &&
              !at[neigh].num_H && 1 == at[neigh].valence &&
              2 == at[neigh].chem_bonds_valence ) {
            continue;
        }
        num_DB += at[i].bond_type[j] == BOND_TYPE_DOUBLE;
    }
    return num_DB && !num_arom;
}
/***************************************************************/
int is_C_Alk( inp_ATOM *at, int i, char cFlags )
{
    if ( at[i].el_number == el_number_C &&
         at[i].valence == at[i].chem_bonds_valence ) {
        int j, k;
        U_CHAR el;
        for ( j = 0; j < at[i].valence; j ++ ) {
            k = at[i].neighbor[j];
            if ( at[k].cFlags & cFlags )
                continue;
            el = at[k].el_number;
            if ( el != el_number_C &&
                 el != el_number_F &&
                 el != el_number_Cl &&
                 el != el_number_Br &&
                 el != el_number_I ) {
                return 0;
            }
        }
        return 1;
    }
    return 0;
}
/***************************************************************/
int is_Si_IV( inp_ATOM *at, int i )
{
    if ( at[i].el_number != el_number_Si ||
         at[i].charge || at[i].radical || at[i].valence != 4 || at[i].chem_bonds_valence != 4 )
        return 0;
    return 1;
}
/***************************************************************/
int is_P_TB_N( inp_ATOM *at, int i )
{
    int j, k;
    if ( at[i].el_number != el_number_P || at[i].chem_bonds_valence - at[i].valence != 2 )
        return 0;
    for ( j = 0; j < at[i].valence; j ++ ) {
        k = at[i].neighbor[j];
        if ( at[k].el_number == el_number_N && at[k].valence == 1 && at[k].chem_bonds_valence == 3 )
            return 1;
    }
    return 0;
}
/***************************************************************/
int is_possibly_deriv_neigh( inp_ATOM *at, int iat, int iord, int type, char cFlags )
{
    int neigh = at[iat].neighbor[iord];
    int neigh_from = -1;
    U_CHAR el = at[neigh].el_number;
    int    bOk = 0;
    switch ( type ) {
    case DERIV_BRIDGE_O:
        neigh_from = at[iat].neighbor[!iord];
        /* -> A--O--B -> traversing from A(neigh_from) to B(neigh); may we cut O--B bond? */
        /* do not cut bond "---" in A=Si(IV), B(=O), B=C: Si(IV)-O---B(=O) */
        if ( !(is_C_or_S_DB_O( at, neigh ) && is_Si_IV( at, neigh_from )) &&
             !is_C_unsat_not_arom( at, neigh ) ) {
            bOk = ( el == el_number_C ||
                    el == el_number_Si ||
                    el == el_number_S ||
                    el == el_number_P );
        }
        break;
    case DERIV_BRIDGE_NH:
        /* -> A--NH--B -> traversing from A(neigh_from) to B(neigh); may we cut NH--B bond? */
        bOk = ( is_C_or_S_DB_O( at, neigh ) ||
                is_C_Alk( at, neigh, cFlags ) ||
                is_Si_IV( at, neigh ) ||
                is_P_TB_N( at, neigh ) ) && !(is_C_unsat_not_arom( at, neigh ));
        break;
    case DERIV_AMINE_tN:
        bOk = ( is_C_or_S_DB_O( at, neigh ) ||
                is_C_Alk( at, neigh, cFlags ) ||
                is_Si_IV( at, neigh ) ||
                is_P_TB_N( at, neigh ) ) && !(is_C_unsat_not_arom( at, neigh ));
        break;
    }
    return bOk;
}
/***************************************************************/
/* determines derivative type on the forward step of the DFS   */
/***************************************************************/
int get_traversed_deriv_type( inp_ATOM *at, DERIV_AT *da, int k, DERIV_AT *da1, char cFlags )
{
    int i, j, m, nBlockSystemFrom, nOrdBack1, nOrdBack2, nOrdBack3, nBackType1, nBackType2;
    memset( da1, 0, sizeof(*da1) );
    if ( at[k].cFlags & cFlags ) {
        return 0;
    }
    for ( m = 0; m < at[k].valence && !(at[(int)at[k].neighbor[m]].cFlags & cFlags); m ++ )
        ;
    if ( m == at[k].valence )
        return -1; /* error: at least one neighbor must have cFlags */
    if ( at[k].valence == 1 && at[k].num_H && (
           at[k].el_number == el_number_O ||
           at[k].el_number == el_number_N ||
           at[k].el_number == el_number_S ||
           at[k].el_number == el_number_P ) ) {
        return DERIV_NOT;
    }
    if ( is_el_a_metal( at[k].el_number ) ) {
        return DERIV_NOT;
    }
#ifdef NEVER
    if ( at[k].el_number == el_number_N && at[k].valence >= 2 && at[k].chem_bonds_valence <= 3 ) {
        return DERIV_NOT; /* prohibit -N-, -N=, allow -N# as in isocyano -N#C or NO2 */
    }
#endif
    /* m is the ord of the bond from which at[k] was reached first time */
    if ( at[k].nNumAtInRingSystem == 1 && (at[k].el_number == el_number_O || at[k].el_number == el_number_S) &&
         at[k].valence == 2 && at[k].chem_bonds_valence == 2 &&
         !at[k].num_H && !at[k].charge && !at[k].radical) {
        /* -> A--O--B -> traversing from A to B; cut O--B */
        /* check for carboxy A(=O)-O-B and A--O--B(=O) */
        /* int has_A_CO   = is_C_or_S_DB_O( at, at[k].neighbor[m] ); */
        int has_B_CO   = is_C_or_S_DB_O( at, at[k].neighbor[!m] );
        int is_A_Si_IV = is_Si_IV( at, at[k].neighbor[m] );
        /* int is_B_Si_IV = is_Si_IV( at, at[k].neighbor[!m] );*/
        if ( is_A_Si_IV && has_B_CO ) {
            ; /* do not cut bond --- in A=Si(IV), B(=O), B=C: Si(IV)-O---B(=O) */
        } else
        if ( is_possibly_deriv_neigh( at, k, !m, DERIV_BRIDGE_O, cFlags ) ) {
            da1->ord[0] =  !m;         /* ord of neighbor B, not traversed yet */
            da1->typ[0] =  DERIV_BRIDGE_O;   /* type */
            return DERIV_BRIDGE_O;   /* R-O-R */
        }
    }
    if ( at[k].bCutVertex && at[k].el_number == el_number_N &&
         at[k].valence == 2 && at[k].chem_bonds_valence == at[k].valence &&
         at[k].valence+at[k].num_H == 3 && !at[k].charge && !at[k].radical) {
        da1->ord[0] =  !m;         /* ord of neighbor B, not traversed yet */
        da1->typ[0] =  DERIV_BRIDGE_NH;   /* type */
        return DERIV_BRIDGE_NH;   /* R1-N(-R2)-R3 or R1-NH-R2  amine */
    }
    if ( at[k].bCutVertex && at[k].el_number == el_number_N &&
         at[k].valence == 3 && at[k].chem_bonds_valence == at[k].valence &&
         at[k].valence+at[k].num_H == 3 && !at[k].charge && !at[k].radical) {
        int rm1 = ( at[at[k].neighbor[m]].nRingSystem == at[at[k].neighbor[(m+1)%3]].nRingSystem );
        int rm2 = ( at[at[k].neighbor[m]].nRingSystem == at[at[k].neighbor[(m+2)%3]].nRingSystem );
        int r12 = ( at[at[k].neighbor[(m+1)%3]].nRingSystem == at[at[k].neighbor[(m+2)%3]].nRingSystem );
        int ord[2]= {-1, -1};
        i = 0; /* type: tertriary amine */
        switch( rm1 + rm2 + r12 ) {
        case 0:
            /* -N< has no ring bonds */
            if ( is_possibly_deriv_neigh( at, k, (m+1)%3, DERIV_AMINE_tN, cFlags ) ) {
                ord[i ++] = (m+1)%3; /* ord of a non-ring neighbor, not traversed yet */
            }
            if ( is_possibly_deriv_neigh( at, k, (m+2)%3, DERIV_AMINE_tN, cFlags ) ) {
                ord[i ++] = (m+2)%3; /* ord of another non-ring neighbor, not traversed yet */
            }
            if ( i == 2 && ord[0] > ord[1] ) {
                int tmp = ord[0];
                ord[0] = ord[1];
                ord[1] = tmp;
            }
            break;

        case 1:
            if ( rm1 && is_possibly_deriv_neigh( at, k, (m+2)%3, DERIV_AMINE_tN, cFlags ) ) {
                ord[i++] = (m+2)%3;   /* ord of a single non-ring neighbor, not traversed yet */
            } else
            if ( rm2 && is_possibly_deriv_neigh( at, k, (m+1)%3, DERIV_AMINE_tN, cFlags ) ) {
                ord[i++] = (m+1)%3; /* ord of a single non-ring neighbor, not traversed yet */
            }
        }
        for ( j = 0; j < i; j ++ ) {
            da1->ord[j] = ord[j];
            da1->typ[j] = DERIV_AMINE_tN;
        }
        if ( i ) {
            return DERIV_AMINE_tN;
        }
        return 0; /* all neighbors or two untraversed edges are in one ring system */
    }
    if ( at[k].bCutVertex && /* DD */
         at[k].valence == at[k].chem_bonds_valence &&
         (!at[k].num_H || at[k].el_number == el_number_C && 1 == at[k].num_H) &&
         !at[k].charge && !at[k].radical &&
         (at[k].el_number == el_number_C  && at[k].valence+at[k].num_H == 4 ||
          at[k].el_number == el_number_Si && at[k].valence == 4 ||
          at[k].el_number == el_number_B  && at[k].valence == 3) ) {

        /*                  entering path: ->X--O--DD
            --X--O          
              |    \   /    DD = C, Si, B
              |     DD
              |    /   \     O = O, NH
            --Y--O
                            X, Y -- ignored for now
         */
        nBlockSystemFrom = 0;
        nBackType1 = nBackType2 = 0;
        nOrdBack1 = nOrdBack2 = nOrdBack3 = -1;
        j=(int)at[k].neighbor[m];
        if ( (at[j].el_number == el_number_O || at[j].el_number == el_number_S) && at[j].valence == 2 &&
             at[j].chem_bonds_valence == at[j].valence &&
             at[j].nNumAtInRingSystem >= 5 &&
             !at[j].num_H && !at[j].charge && !at[j].radical ) {
            nBackType1 = DERIV_RING_O;
            nBlockSystemFrom = at[j].nBlockSystem;
            nOrdBack1 = m; /* predecessor */
        } else
        if ( at[j].el_number == el_number_N && at[j].valence == 2 &&
             at[j].chem_bonds_valence == at[j].valence &&
             at[j].nNumAtInRingSystem >= 5 &&
             1==at[j].num_H && !at[j].charge && !at[j].radical ) {
            nBackType1 = DERIV_RING_NH;
            nBlockSystemFrom = at[j].nBlockSystem;
            nOrdBack1 = m; /* predecessor */
        }
        if ( nBlockSystemFrom ) {
            int num1, num2, bFound;
            at[k].cFlags |= CFLAG_MARK_BLOCK;
            num1 = mark_atoms_cFlags( at, at[k].neighbor[nOrdBack1], 1, CFLAG_MARK_BLOCK );
            for ( i = 0; i < at[k].valence; i ++ ) {
                if ( i == nOrdBack1 )
                    continue;
                j=(int)at[k].neighbor[i];
                bFound = 0;
                if ( at[j].cFlags & CFLAG_MARK_BLOCK ) {
                    if ( (at[j].el_number == el_number_O || at[j].el_number == el_number_S) && at[j].valence == 2 &&
                         at[j].chem_bonds_valence == at[j].valence &&
                         at[j].nNumAtInRingSystem >= 5 &&
                         !at[j].num_H && !at[j].charge && !at[j].radical ) {
                       bFound = 1;
                       if ( nOrdBack2 < 0 ) {
                            nOrdBack2 = i; /* predecessor #2 */
                            nBackType2 = DERIV_RING_O;
                        } else {
                            nOrdBack3 = i; /* predecessor #3 -- should not happen */
                        }
                    }
                    if ( at[j].el_number == el_number_N && at[j].valence == 2 &&
                         at[j].chem_bonds_valence == at[j].valence &&
                         at[j].nNumAtInRingSystem >= 5 &&
                         1==at[j].num_H && !at[j].charge && !at[j].radical ) {
                        bFound = 1;
                        if ( nOrdBack2 < 0 ) {
                            nOrdBack2 = i; /* predecessor #2 */
                            nBackType2 = DERIV_RING_NH;
                        } else {
                            nOrdBack3 = i; /* predecessor #3 -- should not happen */
                        }
                    }
                    if ( !bFound ) {
                        nOrdBack3 = 99; /* reject: wrong neighboring atom in the same block */
                        break;
                    }
                }
            }
            num2 = un_mark_atoms_cFlags( at, k, 0, CFLAG_MARK_BLOCK, CFLAG_MARK_BLOCK_INV );
            if ( num1 != num2 ) {
                return -1; /* mark_atoms_cFlags() program error */ 
            }
            if ( nOrdBack2 >= 0 && nOrdBack3 < 0 ) {
                if ( nOrdBack1 < nOrdBack2 ) {
                    da1->ord[0] = nOrdBack1;  /* ord of a ring neighbor, traversed */
                    da1->typ[0] = nBackType1;
                    da1->ord[1] = nOrdBack2;  /* ord of another ring neighbor, not traversed yet */
                    da1->typ[1] = nBackType2;
                } else {
                    da1->ord[0] = nOrdBack2;  /* ord of a ring neighbor, traversed */
                    da1->typ[0] = nBackType2;
                    da1->ord[1] = nOrdBack1;  /* ord of another ring neighbor, not traversed yet */
                    da1->typ[1] = nBackType1;
                }
                return nBackType1 | nBackType2;
            }
        }
    }
    return 0;
}
/***************************************************************/
int add_to_da( DERIV_AT *da, DERIV_AT *add )
{
    /* if add has more than 1 element the elements are in ascending add.ord[] order */
    int i, j, len;
    for ( len = 0; len < DERIV_AT_LEN && da->typ[len]; len ++ )
        ;
    for ( j = 0; j < DERIV_AT_LEN && add->typ[j]; j ++ ) {
        for ( i = 0; i < len; i ++ ) {
            if ( add->ord[j] == da->ord[i] ) {
                if ( add->typ[j] != da->typ[i] ) {
                    return -1; /* error, should not happen */
                }
                if ( add->num[j] != da->num[i] ) {
                    return -2; /* error, should not happen */
                }
                break;
            }
        }
        if ( i == len ) {
            if ( len < DERIV_AT_LEN ) {
                da->ord[i] = add->ord[j];
                da->typ[i] = add->typ[j];
                da->num[i] = add->num[j];
                len ++;
            } else {
                return -3; /* overflow, should not happen */
            }
        }
    }
    return 0;
}
/***************************************************************/
/* DFS search for atoms that do not have a flag */ 
int mark_atoms_deriv( inp_ATOM *at, DERIV_AT *da, int start, int num, char cFlags, int *pbFound )
{
    int i, nFound=0;
    DERIV_AT da1;
    if ( !(at[start].cFlags & cFlags) ) {
        if ( DERIV_NOT == get_traversed_deriv_type( at, da, start, &da1, cFlags ) ) {
            nFound ++; /* at[start] cannot belong to a derivatizing agent */
        }
        num ++;
        at[start].cFlags |= cFlags;
        if ( da1.typ[0] ) {
            /* possibly a derivatization agent attachment point. */
            /* check neighbors that have not been traversed yet */
            int n1=0, n2=0, i1=-1, i2=-1, nFound1=0, nFound2=0;
            switch( da1.typ[0] ) {
            case DERIV_BRIDGE_O:
            case DERIV_BRIDGE_NH:
                n1 = mark_atoms_deriv( at, da, at[start].neighbor[(int)da1.ord[0]], 0, cFlags, &nFound1 );
                if ( n1 > MAX_AT_DERIV || nFound1 ) {
                    da1.typ[0] = 0;
                } else {
                    da1.num[0] = n1;
                    nFound ++;
                }
                break;
            case DERIV_AMINE_tN:
                n1 = mark_atoms_deriv( at, da, at[start].neighbor[i1=da1.ord[0]], 0, cFlags, &nFound1 );
                if ( da1.typ[1] ) {
                    n2 = mark_atoms_deriv( at, da, at[start].neighbor[i2=da1.ord[1]], 0, cFlags, &nFound2 );
                }
                if ( 0 < n1 && n1 <= MAX_AT_DERIV && !nFound1 ) {
                    da1.num[0] = n1;
                    i = 1;
                    nFound ++;
                } else {
                    da1.ord[0] = da1.ord[1];
                    da1.num[0] = da1.num[1];
                    da1.typ[0] = da1.typ[1];
                    da1.typ[1] = 0;
                    i = 0;
                }
                if ( 0 < n2 && n2 <= MAX_AT_DERIV && !nFound2 ) {
                    da1.num[i] = n2;
                    nFound ++;
                } else {
                    da1.typ[i] = 0;
                }
                break;
            case DERIV_RING_O:
            case DERIV_RING_NH:
                for ( i = 0; i < at[start].valence; i ++ ) {
                    if ( i != da1.ord[0] && i != da1.ord[1] && !(at[at[start].neighbor[i]].cFlags & cFlags) ) {
                        if ( !n1 )
                            n1 = mark_atoms_deriv( at, da, at[start].neighbor[i1=i], 0, cFlags, &nFound1 );
                        else
                            n2 = mark_atoms_deriv( at, da, at[start].neighbor[i2=i], 0, cFlags, &nFound2 );
                    }
                }
                if ( n1+n2+1 > MAX_AT_DERIV || nFound1 || nFound2 ) {
                    da1.typ[1] = da1.typ[0] = 0;
                } else {
                    da1.num[0] = n1;
                    da1.num[1] = n2;
                    nFound ++;
                }
                break;
            }
            if ( n1 < 0 )
                return n1;
            if ( n2 < 0 )
                return n2; /* errors */

            if ( i = add_to_da( da+start, &da1 ) ) {
                return i;  /* error */
            }

            *pbFound += nFound1 + nFound2 + nFound;
            num      += n1 + n2;
        } else {
            *pbFound += nFound;
        }
        for ( i = 0; i < at[start].valence; i ++ ) {
            num = mark_atoms_deriv( at, da, at[start].neighbor[i], num, cFlags, pbFound );
            if ( num < 0 )
                return num;
        }
    }
    /* *pbFound =  number of derivatizer attachment points traversed forward from at[start] */
    return num; /* number of atoms traversed forward from at[start] */
}
/***************************************************************/
int count_one_bond_atoms( inp_ATOM *at, DERIV_AT *da, int start, int ord, char cFlags, int *bFound )
{
    int num = 0;
    if ( !(at[at[start].neighbor[ord]].cFlags & cFlags) ) {
        at[at[start].neighbor[ord]].cFlags |= cFlags;
        num ++;
        num = mark_atoms_deriv( at, da, start, num, cFlags, bFound ); 
    }
    return num;
}
/***************************************************************
  List of allowed derivatives

  Legend:

  ->- marks the bond to be disconnexted: X->-Y => XD + TY
      where TY is a derivatizing agent eventually to be discarded

  Allowed Derivative Types List
  =============================

  DERIV_BRIDGE_O, DERIV_BRIDGE_NH, DERIV_AMINE_tN
  -----------------------------------------------
     CH3                CH3  CH3           CH3   CH3
     |                  |    |             |     |
 X->-Si--CH3        X->-Si---Si--CH3   X->-Si----C--CH3  X= O, NH, N  
     |                  |    |             |     |
     CH3                CH3  CH3           CH3   CH3

     4 atoms           7 atoms             7 atoms        is_silyl()


                                                                      
                                                                      
     O                 O                     O      F         O          
     ||                ||                    ||     |         ||      
  R--C--O->-CH3     R--C--O->-CH2--CH3    R--C--O->-C--F   R--C--O->-CF2-CF2-CF3
                                                    |                 
                                                    F                    
                                                                      
                                                                      

           1 atom             2 atoms            3 atoms          10 atoms
           is_Me_or_Et()      is_Me_or_Et()         is_CF3_or_linC3F7()


 A. DERIV_BRIDGE_NH, DERIV_AMINE_tN 
 -----------------------------------
                                                                      
                                                                      
     O                 O                     O   F             O         
     ||                ||                    ||  |             ||     
 N->-C--CH3        N->-C--CH2--CH3       N->-C---C--F      N->-C--CF2-CF2-CF3   
                                                 |                    
                                                 F                       
                                                                      
                                                                      
                                                       
    3 atoms           5 atoms              8 atoms                12 atoms
    is_Me_or_Et()     is_Me_or_Et()              is_CF3_or_linC3F7()

 DERIV_RING_O (da contains >B- or >C< or >CH- atom)
 ------------

  C----O               R----O                  R----O         
  |     \              |     \     CH3         |     \        
  |      >             |      >   /            |      >       
  |       \            |       \ /             |       \      
  |        B--CH3      |        C              |        CH--Ph
  |       /            |       / \             |       /      
  |      >             |      >   \            |      >       
  |     /              |     /     CH3         |     /        
  C----O               R----O                  R----O         
                                                              
  5-member             5 or 6-member           5 or 6-member


           2 atoms              3 atoms                 7 atoms

 DERIV_RING_NH
 ------------
 
 None in the list

***************************************************************/
int is_silyl( inp_ATOM *at, int start, int ord_prev )
{
    int i, neigh, num_Me=0, iC_IV=-1, iSi_IV=-1, i_C_or_Si_IV;

    if ( at[start].el_number != el_number_Si || at[start].valence != 4 ||
         at[start].valence != at[start].chem_bonds_valence ||
         at[start].charge ||  at[start].radical )
        return 0;
    for ( i = 0; i < at[start].valence; i ++ ) {
        if ( i == ord_prev )
            continue;
        neigh = at[start].neighbor[i];
        if ( at[neigh].charge || at[neigh].radical ||
             at[neigh].valence != at[neigh].chem_bonds_valence)
            return 0;
        if ( at[neigh].valence == 4 ) {
            if ( at[neigh].el_number == el_number_C && iC_IV < 0 && iSi_IV < 0 )
                iC_IV = neigh;
            else
            if ( at[neigh].el_number == el_number_Si && iC_IV < 0 && iSi_IV < 0 )
                iSi_IV = neigh;
            else
                return 0;
        } else
        if ( at[neigh].valence == 1 &&
             at[neigh].valence == at[neigh].chem_bonds_valence &&
             at[neigh].el_number == el_number_C && at[neigh].num_H == 3 ) {
            num_Me ++;
        } else {
            return 0;
        }
    }
    if ( num_Me == 3 && iC_IV < 0 && iSi_IV < 0 )
        return 1; /* Si(CH3)3 */
    
    /* next C(IV) or Si(IV) */
    i_C_or_Si_IV = iC_IV >= 0? iC_IV : iSi_IV >= 0? iSi_IV : -1;
    if ( num_Me != 2 || i_C_or_Si_IV < 0 )
        return 0;

    num_Me = 0;
    for ( i = 0; i < at[i_C_or_Si_IV].valence; i ++ ) {
        neigh = at[i_C_or_Si_IV].neighbor[i];
        if ( neigh == start )
            continue;
        if ( at[neigh].charge || at[neigh].radical ||
             at[neigh].valence != at[neigh].chem_bonds_valence)
            return 0;
        if ( at[neigh].valence == 1 &&
             at[neigh].valence == at[neigh].chem_bonds_valence &&
             at[neigh].el_number == el_number_C && at[neigh].num_H == 3 ) {
            num_Me ++;
        } else {
            return 0;
        }
    }
    if ( num_Me == 3 )
        return 2; /* Si(CH3)2Si/C(CH3)3 */
    return 0;
}
/****************************************************************/
int is_Me_or_Et( inp_ATOM *at, int start, int ord_prev )
{
    int i, neigh = -1;
    if ( at[start].el_number != el_number_C || at[start].valence > 2 ||
         at[start].valence != at[start].chem_bonds_valence ||
         at[start].chem_bonds_valence + at[start].num_H != 4 ||
         at[start].charge ||  at[start].radical )
        return 0;
    for ( i = 0; i < at[start].valence; i ++ ) {
        if ( i == ord_prev )
            continue;
        if ( neigh >= 0 )
            return 0;

        neigh = at[start].neighbor[i];
        if ( at[neigh].el_number != el_number_C || at[neigh].valence > 1 ||
             at[neigh].valence != at[neigh].chem_bonds_valence ||
             at[neigh].chem_bonds_valence + at[neigh].num_H != 4 ||
             at[neigh].charge ||  at[neigh].radical )
            return 0;
    }
    return 1 + (neigh >= 0);
}
#ifdef NEVER
/****************************************************************
              CF3
 -CF3  or -CF<
              CF3
*****************************************************************/
int is_CF3_or_isoC3F7( inp_ATOM *at, int start, int ord_prev )
{
    int i, k, num_C_IV, iC_IV[2], neigh, num_F, iC;
    if ( at[start].el_number != el_number_C || at[start].valence != 4 ||
         at[start].valence != at[start].chem_bonds_valence ||
         at[start].chem_bonds_valence + at[start].num_H != 4 ||
         at[start].charge ||  at[start].radical )
        return 0;

    iC_IV[0] = iC_IV[1] = num_F = 0;
    
    for ( i = num_C_IV = 0; i < at[start].valence; i ++ ) {
        if ( i == ord_prev )
            continue;

        neigh = at[start].neighbor[i];
        if ( at[neigh].valence != at[neigh].chem_bonds_valence ||
             at[neigh].charge ||  at[neigh].radical )
            return 0;
        if ( at[neigh].el_number == el_number_F ) {
            if ( at[neigh].chem_bonds_valence + at[neigh].num_H != 1 )
                return 0;
            num_F ++;
        } else
        if ( at[neigh].el_number == el_number_C &&
             at[neigh].valence == 4 &&
             !at[neigh].num_H && !at[neigh].charge && !at[neigh].radical && num_C_IV < 2 ) {

            if ( num_C_IV > 1 )
                return 0;

            iC_IV[num_C_IV++] = neigh;
        }
    }
    if ( !num_C_IV && 3 == num_F )
        return 1; /* -CF3 */
    if ( 2 != num_C_IV || 1 != num_F )
        return 0;

    /* detect iso-C3F7 */
    for ( k = 0; k < num_C_IV; k ++ ) {
        iC = iC_IV[k];
        num_F = 0;
        for ( i = 0; i < at[iC].valence; i ++ ) {
            neigh = at[iC].neighbor[i];
            if ( neigh == start )
                continue;
            if ( at[neigh].valence != at[neigh].chem_bonds_valence ||
                 at[neigh].charge ||  at[neigh].radical )
                return 0;
            if ( at[neigh].el_number == el_number_F ) {
                if ( at[neigh].chem_bonds_valence + at[neigh].num_H != 1 )
                    return 0;
                num_F ++;
            } else {
                return 0;
            }
        }
        if ( num_F != 3 )
            return 0;
    }
    return 2; /* iso-C3F7 */
}
#endif
/**************************************************************/
int is_CF3_or_linC3F7( inp_ATOM *at, int start, int ord_prev )
{
    int i, num_C_IV, iC_IV, neigh, num_F, num_C=0;
    AT_NUMB *p;

    while( num_C < 4 ) {

        if ( at[start].el_number != el_number_C || at[start].valence != 4 ||
             at[start].valence != at[start].chem_bonds_valence ||
             at[start].chem_bonds_valence + at[start].num_H != 4 ||
             at[start].charge ||  at[start].radical )
            return 0;

        iC_IV = num_F = 0;
        
        for ( i = num_C_IV = 0; i < at[start].valence; i ++ ) {
            if ( i == ord_prev )
                continue;

            neigh = at[start].neighbor[i];
            if ( at[neigh].valence != at[neigh].chem_bonds_valence ||
                 at[neigh].charge ||  at[neigh].radical )
                return 0;
            if ( at[neigh].el_number == el_number_F ) {
                if ( at[neigh].chem_bonds_valence + at[neigh].num_H != 1 )
                    return 0;
                num_F ++;
            } else
            if ( at[neigh].el_number == el_number_C &&
                 at[neigh].valence == 4 &&
                 !at[neigh].num_H && !at[neigh].charge && !at[neigh].radical && num_C_IV < 2 ) {

                if ( num_C_IV )
                    return 0;

                iC_IV = neigh;
                num_C_IV++;
            }
        }
        if ( num_C_IV + num_F != 3 )
            return 0;

        num_C ++; /* found -CF2-C or -CF3 */
        if ( !num_C_IV )
            break; /* -CF3 */

        /* treat next C(IV) as a new start atom */
        if ( p = is_in_the_list( at[iC_IV].neighbor, (AT_NUMB) start, at[iC_IV].valence ) ) {
            start = iC_IV;
            ord_prev = p - at[iC_IV].neighbor;
        } else {
            return -1; /* program error */
        }
    }
    return num_C == 1? 1 : num_C == 3? 2 : 0;
}
/****************************************************************/
int is_phenyl( inp_ATOM *at, int start, int ord_prev )
{
    int k, neigh, cur_at, ord;
    if ( at[start].el_number != el_number_C || at[start].valence != 3 ||
         at[start].valence+1 != at[start].chem_bonds_valence ||
         at[start].chem_bonds_valence + at[start].num_H != 4 ||
         at[start].charge ||  at[start].radical )
        return 0;
    
    ord = (ord_prev + 1) % 3;
    cur_at = start;

    for ( k = 0; k < 5; k ++ ) {
        neigh = at[cur_at].neighbor[ord];
        if ( at[neigh].el_number != el_number_C || at[neigh].valence != 2 ||
             at[neigh].valence+1 != at[neigh].chem_bonds_valence ||
             at[neigh].chem_bonds_valence + at[neigh].num_H != 4 ||
             at[neigh].charge ||  at[neigh].radical )
            return 0;
        ord = (at[neigh].neighbor[0] == cur_at);
        cur_at = neigh;
    }
    return (at[cur_at].neighbor[ord] == start);
}
/****************************************************************/
int is_deriv_ring( inp_ATOM *at, int start, int num_atoms, DERIV_AT *da1, int idrv )
{
    int i, k, neigh_at[2], prev_ord[2], neigh, is_B=0, is_C=0;
    AT_NUMB *p;
    if ( da1->typ[idrv] != DERIV_RING_O || da1->typ[idrv+1] != DERIV_RING_O )
        return 0;
    if ( at[start].charge || at[start].radical || at[start].valence != at[start].chem_bonds_valence )
        return 0;
    if ( at[start].el_number == el_number_B && at[start].valence == 3 )
        is_B = 1;
    else
    if ( at[start].el_number == el_number_C && (at[start].valence == 3 || at[start].valence == 4) &&
         at[start].chem_bonds_valence == at[start].valence &&
         at[start].num_H + at[start].chem_bonds_valence == 4 )
        is_C = 1;
    else
        return 0;
    /* locate bonds connecting >B- or >C< or >C- to the rest of the derivative */     
    for ( i = k = 0; i < at[start].valence; i ++ ) {
        if ( i != da1->ord[idrv] && i != da1->ord[idrv+1] ) {
            neigh = at[start].neighbor[i];
            if ( k < 2 && (p = is_in_the_list( at[neigh].neighbor, (AT_NUMB) start, at[neigh].valence)) ) {
                neigh_at[k] = neigh;
                prev_ord[k] = p - at[neigh].neighbor;
                k ++;
            } else {
                return -1; /* program error */
            }
        }
    }
    if ( is_B && k == 1 && is_Me_or_Et( at, neigh_at[0], prev_ord[0]) )
        return 1;
    if ( is_C && k == 1 && at[start].num_H == 1 && is_phenyl( at, neigh_at[0], prev_ord[0]) )
        return 1;
    if ( is_C && k == 2 && is_Me_or_Et( at, neigh_at[0], prev_ord[0]) &&
                           is_Me_or_Et( at, neigh_at[1], prev_ord[1]) )
        return 1;

    return 0;
}
/****************************************************************/
int is_deriv_chain( inp_ATOM *at, int start, int num_atoms, DERIV_AT *da1, int idrv )
{
    int i, k, prev_ord, neigh, iC, iO, iNxt;
    AT_NUMB *p;
    if ( !da1->typ[idrv] || (da1->typ[idrv] & DERIV_RING) )
        return 0;
    if ( at[start].charge || at[start].radical || at[start].valence != at[start].chem_bonds_valence )
        return 0;

    neigh = at[start].neighbor[(int)da1->ord[idrv]];
    p = is_in_the_list( at[neigh].neighbor, (AT_NUMB) start, at[neigh].valence);
    if ( !p )
        return -1; /* program error */
    prev_ord = p - at[neigh].neighbor;

    /* eliminate silyl possibility */
    if ( is_silyl( at, neigh, prev_ord ) )
        return 1;

    if ( da1->typ[idrv] == DERIV_BRIDGE_O ) {
        /* check acetyl */
        iC = at[start].neighbor[!da1->ord[idrv]];
        if ( at[iC].charge || at[iC].radical || at[iC].num_H ||
             at[iC].el_number != el_number_C || at[iC].valence != 3 ||
             at[iC].valence+1 != at[iC].chem_bonds_valence )
            return 0;
        for ( i = k = 0; i < at[iC].valence; i ++ ) {
            iO = at[iC].neighbor[i];
            if ( at[iO].charge || at[iO].radical || at[iO].num_H ||
                 at[iO].el_number != el_number_O || at[iO].valence != 1 ||
                 at[iO].valence+1 != at[iO].chem_bonds_valence )
                continue;
            k ++; /* number of =O */
        }
        if ( k != 1 )
            return 0;
        /* check derivative */
        return ( is_Me_or_Et( at, neigh, prev_ord ) ||
                 is_CF3_or_linC3F7( at, neigh, prev_ord ) );
    }

    if ( da1->typ[idrv] == DERIV_BRIDGE_NH || da1->typ[idrv] == DERIV_AMINE_tN ) {
        /* check acetyl */
        iNxt = -1;
        iC = at[start].neighbor[(int)da1->ord[idrv]];
        if ( at[iC].charge || at[iC].radical || at[iC].num_H ||
             at[iC].el_number != el_number_C || at[iC].valence != 3 ||
             at[iC].valence+1 != at[iC].chem_bonds_valence )
            return 0;
        for ( i = k = 0; i < at[iC].valence; i ++ ) {
            iO = at[iC].neighbor[i];
            if ( at[iO].charge || at[iO].radical || at[iO].num_H ||
                 at[iO].el_number != el_number_O || at[iO].valence != 1 ||
                 at[iO].valence+1 != at[iO].chem_bonds_valence ) {
                if ( iO != start ) {
                    if ( iNxt < 0 )
                        iNxt = iO;
                    else
                        return 0;
                }
                continue;
            }
            k ++; /* number of =O */
        }
        if ( k != 1 || iNxt < 0 )
            return 0;
        /* find bond from iNxt to iC */
        p = is_in_the_list( at[iNxt].neighbor, (AT_NUMB) iC, at[iNxt].valence);
        if ( !p )
            return -1; /* program error */
        prev_ord = p - at[iNxt].neighbor;
        /* check derivative */
        return ( is_Me_or_Et( at, iNxt, prev_ord ) ||
                 is_CF3_or_linC3F7( at, iNxt, prev_ord ) );
    }
    return 0;
}
/****************************************************************/
int is_deriv_chain_or_ring( inp_ATOM *at, int start, int num_atoms, DERIV_AT *da1, int *idrv )
{
    int i, ret = -1;
    if ( da1->typ[*idrv] & DERIV_RING ) {
        /* find the first ord of this derivative; ord of ring derivatives are in pairs */
        int j = -1;
        for ( i = 0; i < DERIV_AT_LEN && da1->typ[i]; i ++ ) {
            if ( da1->typ[i] & DERIV_RING ) {
                if ( i == *idrv || i+1 == *idrv ) {
                    *idrv = j = i;
                    break;
                }
                i ++; /* bypass the second bond to the same derivatization agent */
            }
        }
        /* check consistency */
        if ( j == -1 || j+1 >= DERIV_AT_LEN ||
             !(da1->typ[j] & DERIV_RING) || !(da1->typ[j+1] & DERIV_RING) ) {
            ret = -1; /* program error */
        } else {
            ret = is_deriv_ring( at, start, num_atoms, da1, j );
        }
    } else
    if ( da1->typ[*idrv] ) {
        ret = is_deriv_ring( at, start, num_atoms, da1, *idrv );
    }
    return ret;
}
/******************************************************/
int remove_deriv( DERIV_AT *da1, int idrv )
{
    int i, j, ret = -1;
    if ( da1->typ[idrv] & DERIV_RING ) {
        /* find the first ord of this derivative; ord of ring derivatives are in pairs */
        j = -1;
        for ( i = 0; i < DERIV_AT_LEN && da1->typ[i]; i ++ ) {
            if ( da1->typ[i] & DERIV_RING ) {
                if ( i == idrv || i+1 == idrv ) {
                    j = i;
                    break;
                }
                i ++; /* bypass the second bond to the same derivatization agent */
            }
        }
        /* delete if data are consistent */
        if ( j >= 0 && j+1 < DERIV_AT_LEN && (da1->typ[j] & DERIV_RING) && (da1->typ[j+1] & DERIV_RING) ) {
            for ( ; j < DERIV_AT_LEN-2 && da1->typ[j+2]; j ++ ) {
                da1->typ[j] = da1->typ[j+2];
                da1->num[j] = da1->num[j+2];
                da1->ord[j] = da1->ord[j+2];
            }
            for ( ; j < DERIV_AT_LEN; j ++ ) {
                da1->typ[j] = 0;
                da1->num[j] = 0;
                da1->ord[j] = 0;
            }
            ret = 0;
        }
    } else {
        j = idrv;

        for ( ; j < DERIV_AT_LEN-1 && da1->typ[j+1]; j ++ ) {
            da1->typ[j] = da1->typ[j+1];
            da1->num[j] = da1->num[j+1];
            da1->ord[j] = da1->ord[j+1];
        }
        for ( ; j < DERIV_AT_LEN; j ++ ) {
            da1->typ[j] = 0;
            da1->num[j] = 0;
            da1->ord[j] = 0;
        }
        ret = 0;
    }
    return ret;
}
/******************************************************/
int remove_deriv_mark( DERIV_AT *da1, int idrv )
{
    int i, j, ret = -1;
    if ( da1->typ[idrv] & DERIV_RING ) {
        /* find the first ord of this derivative; ord of ring derivatives are in pairs */
        j = -1;
        for ( i = 0; i < DERIV_AT_LEN && da1->typ[i]; i ++ ) {
            if ( da1->typ[i] & DERIV_RING ) {
                if ( i == idrv || i+1 == idrv ) {
                    j = i;
                    break;
                }
                i ++; /* bypass the second bond to the same derivatization agent */
            }
        }
        /* delete if data are consistent */
        if ( j >= 0 && j+1 < DERIV_AT_LEN && (da1->typ[j] & DERIV_RING) && (da1->typ[j+1] & DERIV_RING) ) {
            da1->typ[j] |= DERIV_DUPLIC;
            da1->typ[j+1] |= DERIV_DUPLIC;
            ret = 0;
        }
    } else {
        j = idrv;
        da1->typ[j] |= DERIV_DUPLIC;
        ret = 0;
    }
    return ret;
}
/****************************************************************/
int EliminateDerivNotInList( inp_ATOM *at, DERIV_AT *da, int num_atoms )
{
    int i, j, num_da, num_cuts=0, ret=0;
    for ( i = 0; i < num_atoms; i ++ ) {
        if ( !da[i].typ[0] )
            continue;
        /* count deriative attachments */
        for ( num_da = 0; num_da < DERIV_AT_LEN && da[i].typ[num_da]; num_da ++ )
            ;
        if ( num_da > 2 )
            return -1; /* should not happen */
        if ( num_da == 2 && da[i].typ[0] != da[i].typ[1] ) {
            da[i].typ[0] = da[i].typ[1] = 0; /* do not allow */
            continue;
        }
        if ( da[i].typ[0] & DERIV_RING ) {
            ret = 0;
            if ( num_da == 2 && 1 + da[i].num[0] + da[i].num[1] <= MAX_AT_DERIV &&
                 0 < (ret=is_deriv_ring( at, i, num_atoms, da+i, 0 ) ) ) {
                num_cuts += 2;
            } else
            if ( ret < 0 ) {
                return ret;
            } else {
                da[i].typ[0] = da[i].typ[1] = 0; /* not a derivative */
            }
        } else {
            ret = 0;
            if ( da[i].num[0] <= MAX_AT_DERIV && 0 < (ret = is_deriv_chain( at, i, num_atoms, da+i, 0 )) ) {
                num_cuts ++;
                j = 1;
            } else
            if ( ret < 0 ) {
                return ret;
            } else {
                da[i].ord[0] = da[i].ord[1];
                da[i].num[0] = da[i].num[1];
                da[i].typ[0] = da[i].typ[1];
                da[i].typ[1] = 0;
                j = 0;
            }
            if ( da[i].typ[j] && da[i].num[j] <= MAX_AT_DERIV &&
                 0 < (ret = is_deriv_chain( at, i, num_atoms, da+i, j )) ) {
                num_cuts ++;
            } else
            if ( ret < 0 ) {
                return ret;
            } else {
                da[i].typ[j] = 0;
            }
        }
    }
    return num_cuts;
}
/***************************************************************/
int make_single_cut( inp_ATOM *at, DERIV_AT *da, int iat, int icut )
{
    int ret = -1; /* error flag */
    int iord = (int)da[iat].ord[icut]; /* ord of the bond in iat */
    if ( da[iat].typ[icut] & DERIV_DUPLIC ) {
        return 0;
    } else
    if ( iord < 0 ) {
        return -1; /* program error */
    } else {
        /* find other da[] that affect at[iat] */
        int jat  = at[iat].neighbor[iord];  /* opposite atom */
        AT_NUMB *p = is_in_the_list( at[jat].neighbor, (AT_NUMB) iat, at[jat].valence );
        int jord = p? (p - at[jat].neighbor) : -1;
        int i, iD=1, iT=2;
        if ( jord < 0 ) {
            return -1;  /* program error */
        }
        ret = DisconnectInpAtBond( at, NULL, iat, iord );
        if ( ret == 1 ) {
            if ( da[iat].typ[icut] & DERIV_RING ) {
                /* at[jat] belongs to the main structure */
                at[jat].num_H ++;        /* add D to the main structure */
                at[jat].num_iso_H[iD] ++;
                at[iat].num_H ++;        /* add T to the derivatizing fragment */
                at[iat].num_iso_H[iT] ++;
            } else {
                at[iat].num_H ++;        /* add D to the main structure */
                at[iat].num_iso_H[iD] ++;
                at[jat].num_H ++;        /* add T to the derivatizing fragment */
                at[jat].num_iso_H[iT] ++;
            }
            /* adjust ord for other bonds */
            for ( i = 0; i < DERIV_AT_LEN && da[iat].typ[i]; i ++ ) {
                if ( da[iat].ord[i] == iord ) {
                    da[iat].ord[i] = -(1+da[iat].ord[i]); /* mark the bond being disconnected */
                } else
                if ( da[iat].ord[i] > iord ) {
                    da[iat].ord[i] --;
                }
            }
            for ( i = 0; i < DERIV_AT_LEN && da[jat].typ[i]; i ++ ) {
                if ( da[jat].ord[i] == jord ) {
                    /* opposite atom needs the same bond to be disconnected */
                    if ( da[iat].num[icut] == da[jat].num[i] ) {
                        iD=2; /* mark both as fragments */
                    } else
                    if ( da[iat].num[icut] > da[jat].num[i] ) {
                        iD = 2; /* greater as a main structure */
                        iT = 1; /* mark smaller as a derivatizing fragment */
                    }
                    da[jat].ord[i]  = -(1+da[jat].ord[i]);
                    da[jat].typ[i] |= DERIV_DUPLIC;
                } else
                if ( da[jat].ord[i] > jord ) {
                    da[jat].ord[i] --;
                }
            }
        }
    }
    return ret;
}
/***************************************************************/
int underivatize( ORIG_ATOM_DATA *orig_inp_data )
{
#define REMOVE_CUT_DERIV 1  /* remove disconnected derivatizing agents */
    int ret = 0, i, j, k, m, n, num_atoms, num_components, i_component, nFound, num, cur_num_at, len;
    int num_cuts, num_ring_cuts, num_cut_pieces, num_cuts_to_check;
    inp_ATOM *at = orig_inp_data->at;
    INP_ATOM_DATA *inp_cur_data = NULL;
    DERIV_AT      *da           = NULL;
    int  nTotNumCuts = 0;

    set_R2C_el_numbers( );
    /* prepare */
    num_atoms = remove_terminal_HDT( orig_inp_data->num_inp_atoms, at, 1 );
                /*^^^^^ always accomodate accomodate FIX_TERM_H_CHRG_BUG - IPl, July 2008*/
    orig_inp_data->num_inp_atoms = num_atoms;

    /* initialize */
    UnMarkDisconnectedComponents( orig_inp_data );
    num_cuts = 0;
    /* mark */
    num_components = MarkDisconnectedComponents( orig_inp_data, 0 );
    inp_cur_data = (INP_ATOM_DATA *)inchi_calloc( num_components, sizeof(inp_cur_data[0]) );
    for ( i_component = 0; i_component < num_components; i_component ++ ) {
        CreateInpAtomData( inp_cur_data+i_component, orig_inp_data->nCurAtLen[i_component], 0 );
        inp_cur_data[i_component].num_at = ExtractConnectedComponent( orig_inp_data->at, orig_inp_data->num_inp_atoms, i_component+1, inp_cur_data[i_component].at );
        /*  error processing */
        if ( inp_cur_data[i_component].num_at <= 0 || orig_inp_data->nCurAtLen[i_component] != inp_cur_data[i_component].num_at ) {
            ret = -(i_component+1); /* severe error */
            goto exit_function;
        }
        /* initialize */
        num_atoms = inp_cur_data[i_component].num_at;
        at        = inp_cur_data[i_component].at;
        add_DT_to_num_H( num_atoms, at );
        
        UnMarkRingSystemsInp( at, num_atoms );
        UnMarkOtherIndicators( at, num_atoms );
        UnMarkOneComponent( at, num_atoms );
        MarkRingSystemsInp( at, num_atoms, 0 );
        ret = mark_arom_bonds( at, num_atoms );
        if ( ret < 0 ) {
            goto exit_function;
        }
        ret = 0;
        if ( da ) inchi_free( da );
        da = (DERIV_AT *)inchi_calloc( num_atoms, sizeof(da[0]) );
        
        /* detect derivatives */
        nFound = 0;
        for ( i = 0; i < num_atoms; i ++ ) {
            if ( at[i].bCutVertex && !da[i].typ[0] ) {
                for ( k = 0; k < at[i].valence; k ++ ) {
                    num = count_one_bond_atoms( at, da, i, k, CFLAG_MARK_BRANCH, &nFound );
                    UnMarkOtherIndicators( at, num_atoms );
                    if ( num < 0 ) {
                        ret = num; /* severe error */
                        goto exit_function;
                    }
                }
            }
        }
        /* prepare cuts: remove cuts that are not to be done */
        /* in addition, count ring cuts DERIV_RING */
        num_ring_cuts  = 0;
        num_cuts       = 0;
        num_cut_pieces = 0;
        for ( i = num = 0; i < num_atoms; i ++ ) {
            for ( len = 0; len < MAX_AT_DERIV && da[i].typ[len]; len ++ )
                ;
            switch( len ) {
            case 0:
                continue;
            case 1:
                /* single cut: unconditional */
                num_cuts       += 1;
                num_cut_pieces += 1;
                continue;
            case 2:
                if ( (da[i].typ[0] & DERIV_RING) && (da[i].typ[1] & DERIV_RING) ) {
                    /* double cut, unconditional */
                    num_ring_cuts  += 2;
                    num_cuts       += 2;
                    num_cut_pieces += 1;
                    continue;
                } else
                if ( da[i].typ[0] == DERIV_AMINE_tN && da[i].typ[1] == DERIV_AMINE_tN ) {
                    /* double cut, unconditional */
                    num_cuts       += 2;
                    num_cut_pieces += 2;
                    continue;
                }
                if ( da[i].typ[0] == da[i].typ[1] ) {
                    /* DERIV_BRIDGE_O or DERIV_BRIDGE_NH; cut off the smallest */
                    if ( 0 == is_deriv_chain( at, i, num_atoms, da+i, 0 ) ) {
                        da[i].num[0] = NOT_AT_DERIV;
                    }
                    if ( 0 == is_deriv_chain( at, i, num_atoms, da+i, 1 ) ) {
                        da[i].num[1] = NOT_AT_DERIV;
                    }
                    if ( da[i].num[0] > da[i].num[1] ) {
                        da[i].num[0] = da[i].num[1];
                        da[i].ord[0] = da[i].ord[1];
                        da[i].typ[0] = da[i].typ[1];
                        da[i].typ[1] = 0;
                        num_cuts       += 1;
                        num_cut_pieces += 1;
                    } else
                    if ( da[i].num[0] < da[i].num[1] ) {
                        da[i].typ[1] = 0;
                        num_cuts       += 1;
                        num_cut_pieces += 1;
                    } else {
                        /* attachments have same size: ignore both */
                        /* ??? check for standard derivatizations ??? */
                        da[i].typ[0] = 0;
                        da[i].typ[1] = 0;
                    }
                    continue;
                }
                ret = -88;
                goto exit_function; /* unexpected */
            case 3:
                if ( da[i].typ[0] == da[i].typ[1] &&
                     da[i].typ[0] == da[i].typ[2] &&
                     da[i].typ[0] == DERIV_AMINE_tN ) {
                    int x, y, z;
                    if ( 0 == is_deriv_chain( at, i, num_atoms, da+i, 0 ) ) {
                        da[i].num[0] = NOT_AT_DERIV;
                    }
                    if ( 0 == is_deriv_chain( at, i, num_atoms, da+i, 1 ) ) {
                        da[i].num[1] = NOT_AT_DERIV;
                    }
                    if ( 0 == is_deriv_chain( at, i, num_atoms, da+i, 2 ) ) {
                        da[i].num[2] = NOT_AT_DERIV;
                    }

                    x = (da[i].num[0] < da[i].num[1])? 0 : 1;
                    x = (da[i].num[x] < da[i].num[2])? x : 2; /* min */
                    z = (da[i].num[0] < da[i].num[1])? 1 : 0;
                    z = (da[i].num[x] < da[i].num[2])? 2 : z; /* max */
                    y = ((x+1)^(z+1))-1;                      /* median */


                    if (da[i].num[x] == da[i].num[z] ) {
                        /* no cuts */
                        da[i].typ[0] = 0;
                        da[i].typ[1] = 0;
                        da[i].typ[2] = 0;
                        continue; /* all deriv. agents have same size, ignore */
                                  /* ??? check for standard derivatizations ??? */
                    } else
                    if ( da[i].num[x] == da[i].num[y] ) {
                        /* two smallest */
                        switch( z ) {
                        case 0:
                            da[i].num[0] = da[i].num[1];
                            da[i].ord[0] = da[i].ord[1];
                            da[i].typ[0] = da[i].typ[1];
                            
                            da[i].num[1] = da[i].num[2];
                            da[i].ord[1] = da[i].ord[2];
                            da[i].typ[1] = da[i].typ[2];
                            break;
                        case 1:
                            da[i].num[1] = da[i].num[2];
                            da[i].ord[1] = da[i].ord[2];
                            da[i].typ[1] = da[i].typ[2];
                            break;
                        case 2:
                            break;
                        }
                        da[i].typ[2] = 0;
                        num_cuts       += 2;
                        num_cut_pieces += 2;
                    } else {
                        /* one smallest */
                        if ( x ) {
                            da[i].num[0] = da[i].num[x];
                            da[i].ord[0] = da[i].ord[x];
                            da[i].typ[0] = da[i].typ[x];
                        }
                        da[i].typ[1] = 0;
                        da[i].typ[2] = 0;
                        num_cuts       += 1;
                        num_cut_pieces += 1;
                    }
                    continue;
                }
                ret = -88;
                goto exit_function; /* unexpected */
            case 4:
                if ( (da[i].typ[0] & DERIV_RING) && (da[i].typ[1] & DERIV_RING) &&
                     (da[i].typ[2] & DERIV_RING) && (da[i].typ[3] & DERIV_RING) ) {
                    int n01 = inchi_max( da[i].num[0], da[i].num[1] );
                    int n23 = inchi_max( da[i].num[2], da[i].num[3] );
                    if ( n01 < n23 && 0 < is_deriv_ring( at, i, num_atoms, da+i, 0) ) {
                        da[i].typ[2] = 0;
                        da[i].typ[3] = 0;
                        num_cuts       += 2;
                        num_ring_cuts  += 2;
                        num_cut_pieces += 1;
                    } else
                    if ( n01 > n23 && 0 < is_deriv_ring( at, i, num_atoms, da+i, 2) ) {
                        da[i].num[0] = da[i].num[2];
                        da[i].ord[0] = da[i].ord[2];
                        da[i].typ[0] = da[i].typ[2];
                        
                        da[i].num[1] = da[i].num[3];
                        da[i].ord[1] = da[i].ord[3];
                        da[i].typ[1] = da[i].typ[3];

                        da[i].typ[2] = 0;
                        da[i].typ[3] = 0;
                        num_cuts       += 2;
                        num_ring_cuts  += 2;
                        num_cut_pieces += 1;
                    } else {
                        da[i].typ[0] = 0;
                        da[i].typ[1] = 0;
                        da[i].typ[2] = 0;
                        da[i].typ[3] = 0;
                    }
                    continue;
                }
                ret = -88;
                goto exit_function; /* unexpected */
            }
        }
        
        /* eliminate cases when 
             da[i1].typ[j1] && da[i2].typ[j2] &&
             at[i1].neighbor[da[i1].ord[j1]] == i2 && at[i2].neighbor[da[i2].ord[j2]] == i1
           that is, same bond is in the da[] elements of the adjacent atoms */
        nFound = 0; /* number of cuts to remove */
        for ( i = 0; i < num_atoms; i ++ ) {
            for ( j = 0; j < MAX_AT_DERIV && da[i].typ[j]; j ++ ) {
                if ( da[i].typ[j] & DERIV_DUPLIC )
                    continue;
                n = at[i].neighbor[(int)da[i].ord[j]];
                if ( n < i )
                    continue;
                for ( k = 0; k < MAX_AT_DERIV && da[n].typ[k]; k ++ ) {
                    if ( da[n].typ[k] & DERIV_DUPLIC )
                        continue;
                    m = at[n].neighbor[(int)da[n].ord[k]];
                    if ( m == i ) {
                        /* same bond in da[i].typ[j] and da[n].typ[k] */
                        /* check whether both derivatives are acceptable */
                        int k1=k, j1=j;
                        int ret_i = is_deriv_chain_or_ring( at, i, num_atoms, da+i, &j1 );
                        int ret_n = is_deriv_chain_or_ring( at, n, num_atoms, da+n, &k1 );
                        if ( ret_i < 0 ) {
                            ret = ret_i;
                            goto exit_function;
                        }
                        if ( ret_n < 0 ) {
                            ret = ret_n;
                            goto exit_function;
                        }
                        if ( !ret_i || ret_i && ret_n ) {
                            if ( da[i].typ[j1] & DERIV_RING ) {
                                num_cuts       -= 2;
                                num_ring_cuts  -= 2;
                            } else {
                                num_cuts       -= 1;
                            }
                            num_cut_pieces -= 1;
                            if (ret = remove_deriv_mark( da+i, j1 )) {
                                goto exit_function;
                            }
                            nFound ++;
                        }
                        if ( !ret_n || ret_i && ret_n ) {
                            if ( da[n].typ[k1] & DERIV_RING ) {
                                num_cuts       -= 2;
                                num_ring_cuts  -= 2;
                            } else {
                                num_cuts       -= 1;
                            }
                            num_cut_pieces -= 1;
                            if ( ret = remove_deriv_mark( da+n, k1 ) ) {
                                goto exit_function;
                            }
                            nFound ++;
                        }
                    }
                }
            }
        }
        if ( nFound ) {
            for ( i = 0; i < num_atoms; i ++ ) {
                for ( j = 0; j < MAX_AT_DERIV && da[i].typ[j]; j ++ ) {
                   /* attn: j is changed inside the cycle body */
                   if ( da[i].typ[j] & DERIV_DUPLIC ) {
                        if (ret = remove_deriv( da+i, j )) {
                            goto exit_function;
                        }
                        j --;
                    }
                }
            }
        }

        /* make sure DERIV_RING type disconnections increase */
        /* number of components by the number of disconnected derivateves */
        /* Avoid cases like these:

                 O--R--O             DO--R--OD    
                /       \                         
            R--X         Y--R => R--XT2     T2Y--R
                \       /                         
                 O--R--O             DO--R--OD    



                 O--O                 DO--OD
                /    \                      
            R--X--O---Y--R    =>  R--X  OD2 Y--R
                                     T2     T2
                 
        */
        /* count DERIV_RING-type attachments */
#if( COUNT_ALL_NOT_DERIV == 1 )
        num_cuts_to_check = num_cuts;
#else
        num_cuts_to_check = num_ring_cuts;
#endif
        if ( num_cuts_to_check >= 2 )
        {
            /* check */
            R2C_ATPAIR *ap = (R2C_ATPAIR *) inchi_malloc( num_cuts_to_check * sizeof(ap[0]) );
            AT_NUMB    comp_num;
            int        /*n,*/ m_at, m_ord;
            AT_NUMB at1, at2;
            if ( !ap ) {
                ret = -1;
                goto exit_function; /* malloc failure */
            }
repeat_without_deriv_ring:
            comp_num = 0;
            /* fill out the array of bonds to be cut */
            for ( i = j = 0; i < num_atoms; i ++ ) {
                if ( (da[i].typ[0] & DERIV_RING)   && (da[i].typ[1] & DERIV_RING) &&
                      da[i].num[0] <= MAX_AT_DERIV && da[i].num[1] <= MAX_AT_DERIV ) {
                    if ( j+1 >= num_cuts_to_check ) {
                        ret = -2;
                        goto exit_r2c_num; /* wrong number of cuts = num */
                    }
                    for ( k = 0; k < 2; k ++ ) {
                        at1 = i;
                        at2 = at[i].neighbor[(int)da[i].ord[k]];
                        n = ( at1 > at2 );
                        ap[j].at[n] = at1;
                        ap[j].at[1-n] = at2; /* ap[j].at[0] < ap[j].at[1] */
                        j ++;
                    }
                    if ( 0 < cmp_r2c_atpair( ap+j-2, ap+j-1 ) ) {
                        R2C_ATPAIR ap1 = ap[j-2];
                        ap[j-2] = ap[j-1];
                        ap[j-1] = ap1; /* sort each pair */
                    }
                }
#if( COUNT_ALL_NOT_DERIV == 1 )
                else {
                    for ( k = 0; k < DERIV_AT_LEN && da[i].typ[k]; k ++ ) {
                        if ( j >= num_cuts_to_check || (da[i].typ[k] & DERIV_RING) ) {
                            ret = -2;
                            goto exit_r2c_num; /* wrong number of cuts = num or wrong type */
                        }
                        at1 = i;
                        at2 = at[i].neighbor[(int)da[i].ord[k]];
                        n = ( at1 > at2 );
                        ap[j].at[n] = at1;
                        ap[j].at[1-n] = at2; /* ap[j].at[0] < ap[j].at[1] */
                        j ++;
                    }
                }
#endif

            }
            if ( j != num_cuts_to_check ) {
                ret = -3;
                goto exit_r2c_num; /* wrong number of cuts = num */
            }
            /* !!!!!!!! check that there are no derivatives inside a derivative */
            comp_num = 0; /* here it is the number of removed cuts */
            for ( i = 0; i < num_cuts_to_check; i += j ) {
                for ( j = n = 0; j < 2; j ++ ) {
                    int atj = (int)ap[i].at[j];
                    if ( da[atj].typ[0] && at[atj].neighbor[(int)da[atj].ord[0]] == ap[i].at[1-j] ) {
                        k = j;
                        n ++;
                        m_at = atj;
                        m_ord = 0;
                    } else
                    if ( da[atj].typ[1] && at[atj].neighbor[(int)da[atj].ord[1]] == ap[i].at[1-j] ) {
                        k = j;
                        n ++;
                        m_at = atj;
                        m_ord = 1;
                    }
                }
                if ( n != 1 ) {
                    ret = -3;
                    goto exit_r2c_num; /* wrong atom pair */
                }
                if ( (da[m_at].typ[m_ord] & DERIV_RING) ) {
                    n = (int)ap[i].at[k];   /* atom inside the derivation attachment */
                    j = 2;             /* number of bonds to cut */
                    if ( i+j > num_cuts_to_check || (int)ap[i+1].at[0] != n && (int)ap[i+1].at[1] != n ) {
                        ret = -3;
                        goto exit_r2c_num; /* wrong atom pair */
                    }
                } else {
                    n = ap[i].at[1-k]; /* atom inside the derivation attachment */
                    j = 1;             /* number of bonds to cut */
                }

                /* at[n] belongs to the derivation agent */
                cur_num_at = mark_atoms_ap( at, n, ap+i, j, 0, 1 );
                for ( k = 0; k < num_cuts_to_check; k ++ ) {
                    if ( k == i || k==i+j-1 )
                        continue;
                    if ( at[(int)ap[k].at[0]].at_type || at[(int)ap[k].at[1]].at_type ) {
                        /* unmark the cut: found a cut inside the derivatizing agent */
                        da[m_at].typ[m_ord] |= DERIV_UNMARK;
                        num_cuts       -= 1;
                        num_cut_pieces -= 1;
                        if ( j == 2 ) {
                            da[m_at].typ[1-m_ord] |= DERIV_UNMARK;
                            num_cuts       -= 1;
                            num_ring_cuts  -= 2;
                        }
                        comp_num ++;
                        break;
                    }
                }
                UnMarkOtherIndicators( at, num_atoms );
            }
            if ( comp_num ) {
                for ( i = 0; i < num_atoms; i ++ ) {
                    if ( da[i].typ[0] & DERIV_UNMARK ) {
                        da[i].num[0] = da[i].num[1];
                        da[i].ord[0] = da[i].ord[1];
                        da[i].typ[0] = da[i].typ[1];
                        da[i].typ[1] = 0;
                        j = 0;
                    } else {
                        j = 1;
                    }
                    if ( da[i].typ[j] & DERIV_UNMARK ) {
                        da[i].typ[j] = 0;
                    }
                }
#if( COUNT_ALL_NOT_DERIV == 1 )
                num_cuts_to_check = num_cuts;
#else
                num_cuts_to_check = num_ring_cuts;
#endif
                if ( num_cuts < 0 || num_ring_cuts < 0 || num_cut_pieces < 0 ) {
                    ret = -3;
                    goto exit_r2c_num; /* wrong number of cuts = num */
                }
                goto repeat_without_deriv_ring;
            }

            /* sort the bonds for subsequent searching by bisections */
            if ( num_cuts_to_check > 1 ) {
                qsort( ap, num_cuts_to_check, sizeof(ap[0]), cmp_r2c_atpair);
            }
            /* mark components to be disconnected */
            comp_num = 0;   /* number of components */
            cur_num_at = 0; /* number of atoms left after disconnecting the derivatizing agent */
            UnMarkOtherIndicators( at, num_atoms );
            for ( i = 0; i < num_cuts_to_check; i ++ ) {
                n =  0;
                for ( j = 0; j < 2; j ++ ) {
                    if ( da[(int)ap[i].at[j]].typ[0] ) {
                        k = j;
                        n ++;
                    }
                }
                if ( n != 1 ) {
                    ret = -3;
                    goto exit_r2c_num; /* wrong atom pair */
                }
                n = ap[i].at[k]; /* marked atom */
                if ( (da[n].typ[0] & DERIV_RING) ) {
                    n = ap[i].at[1-k];
                }
                /* at[n] belongs to the derivation agent */
                if ( !at[n].at_type ) {
                    comp_num ++;
                    cur_num_at = mark_atoms_ap( at, n, ap, num_cuts_to_check, cur_num_at, comp_num );
                }
            }
            if ( comp_num > 1 ) {
                /* eliminate offending DERIV_RING type derivatives */
                if ( num_ring_cuts <= 2 ) {
                    ret = -99;
                    goto exit_r2c_num;
                }
                n = 0;
                for ( i = j = 0; i < num_atoms; i ++ ) {
                    if ( (da[i].typ[0] & DERIV_RING) && (da[i].typ[1] & DERIV_RING) ) {
                        int at1a = at[i].neighbor[(int)da[i].ord[0]];
                        int at2a = at[i].neighbor[(int)da[i].ord[1]];
                        if ( at[at1a].at_type != at[at2a].at_type ) {
                            da[i].typ[0] = 0; /* eliminate this cut */
                            da[i].typ[1] = 0;
                            n ++;
                            num_cuts_to_check -= 2;
                            num_cuts          -= 2;
                            num_ring_cuts     -= 2;
                            num_cut_pieces    -= 1;
                        }
                    }
                }
                if ( n > 0 && num_cuts_to_check > 2 ) {
                    goto repeat_without_deriv_ring;
                }
            }
            ret = 0;
exit_r2c_num:
            inchi_free( ap );
            UnMarkOtherIndicators( at, num_atoms );
            if ( ret < 0 || num_cuts_to_check >= 2 && cur_num_at < MIN_AT_LEFT_DERIV ) {
                goto exit_function; /* unexpected  error or nothing left */
            }
        }

        if ( !num_cuts ) {
            continue; /*goto exit_function;*/
        }
        /* eliminate derivatives that are not in the list */
        num_cuts = EliminateDerivNotInList( at, da, num_atoms );
        if ( num_cuts < 0 ) {
            ret = num_cuts;
            goto exit_function;
        }


        /* make cuts */
        num_cuts      = 0;
        for ( i = num = 0; i < num_atoms; i ++ ) {
            for ( len = 0; len < MAX_AT_DERIV && da[i].typ[len]; len ++ )
                ;
            switch( len ) {
            case 0:
                continue;
            case 1:
                /* single cut: unconditional */
                make_single_cut( at, da, i, 0 );
                num_cuts += 1;
                continue;
            case 2:
                if ( (da[i].typ[0] & DERIV_RING)    && (da[i].typ[1] & DERIV_RING) ||
                     da[i].typ[0] == DERIV_AMINE_tN && da[i].typ[1] == DERIV_AMINE_tN ) {
                    /* double cut, unconditional */
                    make_single_cut( at, da, i, 1 );
                    make_single_cut( at, da, i, 0 );
                    num_cuts += 1;
                    continue;
                }
                if ( da[i].typ[0] == da[i].typ[1] ) {
                    /* DERIV_BRIDGE_O or DERIV_BRIDGE_NH; cut off the smallest */
                    if ( da[i].num[0] > da[i].num[1] ) {
                        make_single_cut( at, da, i, 1 );
                        num_cuts += 1;
                    } else
                    if ( da[i].num[0] < da[i].num[1] ) {
                        make_single_cut( at, da, i, 0 );
                        num_cuts += 1;
                    }
                    continue;
                }
                ret = -88;
                goto exit_function; /* unexpected */
            case 3:
                if ( da[i].typ[0] == da[i].typ[1] &&
                     da[i].typ[0] == da[i].typ[2] &&
                     da[i].typ[0] == DERIV_AMINE_tN ) {
                    int x, y, z;
                    x = (da[i].num[0] < da[i].num[1])? 0 : 1;
                    x = (da[i].num[x] < da[i].num[2])? x : 2; /* min */
                    z = (da[i].num[0] < da[i].num[1])? 1 : 0;
                    z = (da[i].num[x] < da[i].num[2])? 2 : z; /* max */
                    y = ((x+1)^(z+1))-1;                      /* median */
                    if (da[i].num[x] == da[i].num[z] )
                        continue; /* all deriv. agents have same size */
                    /* two smallest */
                    if ( da[i].num[x] == da[i].num[y] && x < y ) {
                        int t = x; /* first cut x > y */
                        x = y;
                        y = t;
                    }
                    make_single_cut( at, da, i, x );
                    num_cuts += 1;
                    if ( da[i].num[x] == da[i].num[y] ) {
                        /* equally small */
                        make_single_cut( at, da, i, y );
                        num_cuts += 1;
                    }
                    continue;
                }
                ret = -88;
                goto exit_function; /* unexpected */
            case 4:
                if ( (da[i].typ[0] & DERIV_RING) && (da[i].typ[1] & DERIV_RING) &&
                     (da[i].typ[2] & DERIV_RING) && (da[i].typ[3] & DERIV_RING) ) {
                    int n01 = inchi_max( da[i].num[0], da[i].num[1] );
                    int n23 = inchi_max( da[i].num[2], da[i].num[3] );
                    if ( n01 < n23 ) {
                        make_single_cut( at, da, i, 1 );
                        make_single_cut( at, da, i, 0 );
                        num_cuts += 1;
                    } else
                    if ( n01 > n23 ) {
                        make_single_cut( at, da, i, 3 );
                        make_single_cut( at, da, i, 2 );
                        num_cuts += 1;
                    }
                    continue;
                }
            }
        }
        nTotNumCuts += num_cuts;
#if ( REMOVE_CUT_DERIV == 1 )  /* normally YES */
        if ( num_cuts ) {
            /* remove marked with Tritium disconnected derivative attachments */
            ORIG_ATOM_DATA Orig_inp_data1, *orig_inp_data1;
            INP_ATOM_DATA *inp_cur_data1 = NULL;
            int num_components1, i_component1, num_component_left=0;
            orig_inp_data1 = &Orig_inp_data1;
            memset( orig_inp_data1, 0, sizeof(orig_inp_data1[0]) );
            UnMarkRingSystemsInp( at, num_atoms );
            UnMarkOtherIndicators( at, num_atoms );
            UnMarkOneComponent( at, num_atoms );
            for (i = 0; i < num_atoms; i ++ ) {
                orig_inp_data1->num_inp_bonds += at[i].valence;
            }
            orig_inp_data1->num_inp_bonds /= 2;
            orig_inp_data1->num_inp_atoms  = num_atoms;
            orig_inp_data1->at = at; /* = from inp_cur_data[i_component].at */
            num_components1 = MarkDisconnectedComponents( orig_inp_data1, 0 );
            inp_cur_data1 = (INP_ATOM_DATA *)inchi_calloc( num_components1, sizeof(inp_cur_data1[0]) );
            /* extract components and discard disconnected derivatizing agents */
            for ( i_component1 = 0; i_component1 < num_components1; i_component1 ++ ) {
                CreateInpAtomData( inp_cur_data1+i_component1, orig_inp_data1->nCurAtLen[i_component1], 0 );
                inp_cur_data1[i_component1].num_at = ExtractConnectedComponent( orig_inp_data1->at, orig_inp_data1->num_inp_atoms,
                                                                                i_component1+1, inp_cur_data1[i_component1].at );
                /*  error processing */
                if ( inp_cur_data1[i_component1].num_at <= 0 || orig_inp_data1->nCurAtLen[i_component1] != inp_cur_data1[i_component1].num_at ) {
                    ret = -(i_component1+1); /* severe error */
                    break;
                }
                /* if the component has tritium then discard: it is a derivatizing agent */
                for (i = 0; i < inp_cur_data1[i_component1].num_at; i ++ ) {
                    if ( inp_cur_data1[i_component1].at[i].num_iso_H[1] ) {
                        inp_cur_data1[i_component1].at[i].num_iso_H[1] = 0; /* remove deuterium */
                    }
                    if ( inp_cur_data1[i_component1].at[i].num_iso_H[2] ) {
                        FreeInpAtomData( inp_cur_data1+i_component1 );
                        break;
                    }
                }
            }
            /* merge components into one -- must be only one */
            for ( i_component1 = 0, num_atoms = 0; i_component1 < num_components1; i_component1 ++ ) {
                num_atoms += inp_cur_data1[i_component1].num_at;
            }
            at = (inp_ATOM *) inchi_calloc( num_atoms, sizeof(at[0]) );
            cur_num_at = 0;
            for ( i_component1 = 0; i_component1 < num_components1; i_component1 ++ ) {
                /* clean and prepare */
                if ( !inp_cur_data1[i_component1].num_at )
                    continue; /* removed derivatizing object */
                /*UnMarkOneComponent( inp_cur_data1[i_component1].at, inp_cur_data1[i_component1].num_at );*/
                /* merge one by one */
                cur_num_at = add_inp_ATOM( at, num_atoms, cur_num_at, inp_cur_data1[i_component1].at, inp_cur_data1[i_component1].num_at );
                FreeInpAtomData( inp_cur_data1+i_component1 ); /* cleanup */
                num_component_left ++;
            }
            /* replace the component */
            /* order of the following two statements is critically important */
            UnMarkDisconnectedComponents( orig_inp_data1 ); /* orig_inp_data1->at is same as inp_cur_data[i_component].at */
            FreeInpAtomData( inp_cur_data+i_component ); /* cleanup the original component */
            inp_cur_data[i_component].at = at;
            inp_cur_data[i_component].num_at = cur_num_at;
            inchi_free( inp_cur_data1 );
        }
#endif
    }
    if ( nTotNumCuts ) {
        /* merge components into one */
        for ( i = 0, num_atoms = 0; i < num_components; i ++ ) {
            num_atoms += inp_cur_data[i].num_at;
        }
        at = (inp_ATOM *) inchi_calloc( num_atoms, sizeof(at[0]) );
        cur_num_at = 0;
        for ( i = 0; i < num_components; i ++ ) {
            /* clean and prepare */
            UnMarkRingSystemsInp( inp_cur_data[i].at, inp_cur_data[i].num_at );
            UnMarkOtherIndicators( inp_cur_data[i].at, inp_cur_data[i].num_at );
            UnMarkOneComponent( inp_cur_data[i].at, inp_cur_data[i].num_at );
            subtract_DT_from_num_H( inp_cur_data[i].num_at, inp_cur_data[i].at );
            /* merge one by one */
            cur_num_at = add_inp_ATOM( at, num_atoms, cur_num_at, inp_cur_data[i].at, inp_cur_data[i].num_at );
        }
        /* replace orig_inp_data */
        if ( cur_num_at == num_atoms ) {
            inchi_free( orig_inp_data->at );
            orig_inp_data->at = at;
            orig_inp_data->num_inp_atoms = cur_num_at;
            if ( orig_inp_data->szCoord ) {
                inchi_free( orig_inp_data->szCoord );
                orig_inp_data->szCoord = NULL;
            }
            UnMarkDisconnectedComponents( orig_inp_data );
        } else {
            if ( at ) {
                inchi_free( at );
                at = NULL;
            }
            ret = -999; /* num atoms mismatch */
        }
    }
exit_function:
    if ( da ) {
        inchi_free( da );
        da = NULL;
    }
    for ( i_component = 0; i_component < num_components; i_component ++ ) {
        FreeInpAtomData( inp_cur_data+i_component );
    }
    inchi_free( inp_cur_data );
    inp_cur_data = NULL;

    return ret? ret : nTotNumCuts;
}

#endif /* UNDERIVATIZE */
/********************************************************************/
#if( RING2CHAIN == 1 )
/*
    type=1  (incl sugars: W=O, A=C(sat), Z=C(sat), Y=O, B=C(sat)-OH

     A---W               A---WH
    /    |              /     
   |     |        ---> |           
    \    |              \     
     B---Z---YH          B---Z===Y
         |                   |
         |                   |
         C(opt)              C(opt)

    type=2 [not implemented]

     R---W               R---WH
    /     \             /     
   |       Z      ---> |        Z  
    \     /             \     //
     R---YH              R---Y    

*/
#define R2C_EMPTY  127
typedef struct tagRing2Chain {  /* atom Z */
    char type; /* 1 => sugar-like */
    char ordW; /* ordering number of W-neighbor; bond to break; H to add */
    char ordY; /* ordering number of YH-neighbor; bond to increment; H to remove */
    char ordC; /* atom B = C(sat) */
    char ordCopt; /* if exists, saturated C connected by a chain-bond to Z */
} R2C_AT;

int detect_r2c_Zatom( inp_ATOM *at, R2C_AT *da, int iZ );
int cut_ring_to_chain( inp_ATOM *at, R2C_AT *da, int iZ );

/********************************************************************/
int detect_r2c_Zatom( inp_ATOM *at, R2C_AT *da, int iZ )
{
    int i, j, neigh, neighneigh, nRingSystem, num_found;
    R2C_AT da1;
    if ( at[iZ].valence > 4 )
        return 0;
    if ( at[iZ].valence != at[iZ].chem_bonds_valence )
        return 0; /* approach limitation: no double bonds */

    if ( at[iZ].el_number != el_number_C )
        return 0; /* sugar-specific */

    if ( at[iZ].nNumAtInRingSystem < 5 )
        return 0; /* not in a suitable ring */

    if ( !at[iZ].bCutVertex )
        return 0;  /* recognize only type 1 for now */

    nRingSystem = at[iZ].nRingSystem;
    memset ( &da1, R2C_EMPTY, sizeof(da1) );

    for ( i = 0, num_found = 0; i < at[iZ].valence; i ++ ) {
        neigh = at[iZ].neighbor[i];
        if ( at[neigh].charge || at[neigh].radical )
            return 0;
        if ( at[neigh].el_number == el_number_O &&
             at[neigh].valence   == 1 &&
             at[neigh].chem_bonds_valence == 1 &&
             at[neigh].num_H     == 1 ) {
            /* found Z-OH, i.e. Z-YH */
            if ( da1.ordY == R2C_EMPTY ) {
                da1.ordY = i;
                num_found ++;
                continue;
            } else {
                return 0;
            }
        }
        if ( at[neigh].el_number == el_number_O &&
             at[neigh].valence   == 2 &&
             at[neigh].chem_bonds_valence == 2 &&
             at[neigh].num_H     == 0 &&
             at[neigh].nRingSystem == nRingSystem ) {
            /* found Z-O-, i.e. Z-W- */
            if ( da1.ordW == R2C_EMPTY ) {
                 /* j = index of the oppozite to at[iZ] neighbor of at[neigh] */
                j = (at[neigh].neighbor[0] == iZ);
                neighneigh = at[neigh].neighbor[j];
                if ( at[neighneigh].valence != at[neighneigh].chem_bonds_valence ||
                     at[neighneigh].el_number != el_number_C )
                    return 0; /* sugar-specific */
                da1.ordW = i;
                num_found ++;
                continue;
            } else {
                return 0;
            }
        }
        if ( at[neigh].el_number == el_number_C &&
             at[neigh].valence   > 2 &&
             at[neigh].chem_bonds_valence == at[neigh].valence &&
             at[neigh].num_H     <= 1 &&
             at[neigh].nRingSystem == nRingSystem ) {
            /* sugar-specfic: carbon in the ring should have -OH neighbor */
            int iOH;
            for ( j = 0; j < at[neigh].valence; j ++ ) {
                iOH = at[neigh].neighbor[j];
                if ( at[iOH].el_number == el_number_O &&
                     at[iOH].valence == 1 &&
                     at[iOH].chem_bonds_valence == 1 &&
                     at[iOH].num_H == 1 &&
                     !at[iOH].charge && !at[iOH].radical ) {
                    if ( da1.ordC == R2C_EMPTY ) {
                        da1.ordC = i;
                        num_found ++;
                        break;
                    } else {
                        return 0;
                    }
                }
            }
            if ( j < at[neigh].valence )
                continue;
        }
        if ( at[neigh].el_number == el_number_C &&
             at[neigh].chem_bonds_valence == at[neigh].valence &&
             at[neigh].nRingSystem != nRingSystem ) {
            /* extra carbon neighbor of Z */
            if ( da1.ordCopt == R2C_EMPTY ) {
                da1.ordCopt = i;
                continue;
            }
        }
        return 0; /* unexpectd neighbor */
    }
    if (num_found == 3) {
        da1.type = 1;
        da[iZ] = da1;
        return 1; /* disconnection found */
    }
    return 0;
}
/********************************************************************/
int cut_ring_to_chain( inp_ATOM *at, R2C_AT *da, int iZ )
{
    int ret = -1; /* error flag */
    int iordW = (int)da[iZ].ordW; /* ord of the bond in iZ */
    int iordY = (int)da[iZ].ordY; /* ord of the bond in iZ */
    int iordC = (int)da[iZ].ordC;
    int iW, iY, num_iso_H, i, jordZ;
    AT_NUMB *p;

    if ( da[iZ].type != 1 ) {
        return 0;
    }
    if ( 0 > iordW || iordW >= at[iZ].valence ||
         0 > iordY || iordY >= at[iZ].valence ||
         0 > iordC || iordC >= at[iZ].valence /* suger-specific*/) {
        return -1; /* program error */
    }
    /* find other da[] that affect at[iZ] */
    iW  = at[iZ].neighbor[iordW];  /* opposite atom to disconnect and add H */
    iY  = at[iZ].neighbor[iordY];  /* opposite atom to increment the bond and remove H*/
    if ( !at[iY].num_H || at[iZ].bond_type[iordY] != BOND_TYPE_SINGLE ) {
        return -2; /* program error */
    }
    /* increment at[iZ]--at[iY] bond */
    p = is_in_the_list( at[iY].neighbor, (AT_NUMB) iZ, at[iY].valence );
    if ( !p ) {
        return -3; /* program error */
    }
    jordZ = p - at[iY].neighbor;
    at[iZ].bond_type[iordY] ++;
    at[iZ].chem_bonds_valence ++;
    at[iY].bond_type[jordZ] ++;
    at[iY].chem_bonds_valence ++;

    /* disconnect at[iZ]--at[iW] bond */
    ret = DisconnectInpAtBond( at, NULL, iZ, iordW );
    if ( ret != 1 ) {
        return -4; /* program error */
    }
    /* disconnection succeeded */
    /* transfer H from at[iY] to at[iW] */
    num_iso_H = NUM_ISO_H(at, iY);
    if ( at[iY].num_H == num_iso_H ) {
        for ( i = 0; i < NUM_H_ISOTOPES; i ++ ) {
            if ( at[iY].num_iso_H[i] ) {
                at[iY].num_iso_H[i] --;
                at[iW].num_iso_H[i] ++;
            }
        }
    }
    at[iY].num_H --;
    at[iW].num_H ++;
    return 1;
}
/********************************************************************/
int Ring2Chain( ORIG_ATOM_DATA *orig_inp_data )
{
    int ret = 0, i, j, n, num_atoms, num_components, nFound, num, num_cuts, iZ, cur_num_at;
    inp_ATOM *at = orig_inp_data->at;
    INP_ATOM_DATA *inp_cur_data = NULL;
    R2C_AT        *da           = NULL;

    set_R2C_el_numbers( );
    /* prepare */
    num_atoms = remove_terminal_HDT( orig_inp_data->num_inp_atoms, at, 1 );
                /*^^^^^ always accomodate accomodate FIX_TERM_H_CHRG_BUG - IPl, July 2008*/
    orig_inp_data->num_inp_atoms = num_atoms;

    /* initialize */
    UnMarkDisconnectedComponents( orig_inp_data );
    num_cuts = 0;
    /* mark */
    num_components = MarkDisconnectedComponents( orig_inp_data, 0 );
    inp_cur_data = (INP_ATOM_DATA *)inchi_calloc( num_components, sizeof(inp_cur_data[0]) );
    iZ = -1;
    for ( j = 0; j < num_components; j ++ ) {
        CreateInpAtomData( inp_cur_data+j, orig_inp_data->nCurAtLen[j], 0 );
        inp_cur_data[j].num_at = ExtractConnectedComponent( orig_inp_data->at, orig_inp_data->num_inp_atoms, j+1, inp_cur_data[j].at );
        /*  error processing */
        if ( inp_cur_data[j].num_at <= 0 || orig_inp_data->nCurAtLen[j] != inp_cur_data[j].num_at ) {
            ret = -(j+1); /* severe error */
            goto exit_function;
        }
        /* initialize */
        num_atoms = inp_cur_data[j].num_at;
        at        = inp_cur_data[j].at;
        add_DT_to_num_H( num_atoms, at );
        
        UnMarkRingSystemsInp( at, num_atoms );
        UnMarkOtherIndicators( at, num_atoms );
        UnMarkOneComponent( at, num_atoms );
        MarkRingSystemsInp( at, num_atoms, 0 );
        ret = mark_arom_bonds( at, num_atoms );
        if ( ret < 0 ) {
            goto exit_function;
        }
        ret = 0;
        if ( da ) inchi_free( da );
        da = (R2C_AT *)inchi_calloc( num_atoms, sizeof(da[0]) );
        
        /* detect ring-to-chain possibilities */
        nFound = 0;
        for ( i = 0, num=0; i < num_atoms; i ++ ) {
            if ( at[i].bCutVertex /* type 1 specific*/ && !da[i].type ) {
                num += (n=detect_r2c_Zatom( at, da, i ));
                if ( n == 1 )
                    iZ = i;
                UnMarkOtherIndicators( at, num_atoms );
            }
        }
        
        if ( num == 1 ) {
            /* convert ring to chain: make single cut */
            ret = cut_ring_to_chain( at, da, iZ );
            if ( ret < 0 ) {
                goto exit_function;
            }
            num_cuts += (ret == 1);
        } else
        if ( num ) {
            /* allocate an array of bonds to be cut */
            R2C_ATPAIR *ap = (R2C_ATPAIR *)inchi_malloc( sizeof(ap[0]) * num );
            AT_NUMB    comp_num = 0;
            if ( !ap ) {
                ret = -1; /* malloc failure */
                goto exit_function;
            }
            /* fill out the array of bonds to be cut */
            for ( i = j = 0; i < num_atoms; i ++ ) {
                if ( da[i].type == 1 ) {
                    AT_NUMB at1 = i;
                    AT_NUMB at2 = at[i].neighbor[(int)da[i].ordW];
                    if ( j >= num ) {
                        ret = -2;
                        goto exit_r2c_num; /* wrong number of cuts = num */
                    }
                    n = ( at1 > at2 );
                    ap[j].at[n] = at1;
                    ap[j].at[1-n] = at2; /* ap[j].at[0] < ap[j].at[1] */
                    j ++;
                }
            }
            if ( j != num ) {
                ret = -3;
                goto exit_r2c_num; /* wrong number of cuts = num */
            }
            /* sort the bonds for subsequent searching by bisections */
            qsort( ap, num, sizeof(ap[0]), cmp_r2c_atpair);
            /* mark components to be disconnected */
            for ( i = 0; i < num; i ++ ) {
                for ( j = 0; j < 2; j ++ ) {
                    if ( !at[ap[i].at[j]].at_type ) {
                        comp_num ++;
                        mark_atoms_ap( at, (int)ap[i].at[j], ap, num, 0, comp_num );
                    }
                }
            }
            /* convert ring to chain */
            for ( i = 0; i < num; i ++ ) {
                int i1 = ap[i].at[0];
                int i2 = ap[i].at[1];
                iZ = -1;
                if ( at[i1].at_type == at[i2].at_type ) {
                    /* by definition, there are no adjacent iZ atoms; one iZ atom per bond */
                    if ( da[i1].type == 1 && at[i1].neighbor[(int)da[i1].ordW] == i2 ) {
                        iZ = i1;
                    } else
                    if ( da[i2].type == 1 && at[i2].neighbor[(int)da[i2].ordW] == i1 ) {
                        iZ = i2;
                    } else {
                        ret = -3;
                        goto exit_r2c_num;
                    }
                    ret = cut_ring_to_chain( at, da, iZ );
                    if ( ret < 0 ) {
                        goto exit_r2c_num;
                    }
                    num_cuts += (ret == 1);
                }
            }
            ret = 0;
exit_r2c_num:
            inchi_free( ap );
            UnMarkOtherIndicators( at, num_atoms );
            if ( ret < 0 ) {
                goto exit_function;
            }
        }
    }
    if ( num_cuts ) {
        /* merge components into one */
        for ( i = 0, num_atoms = 0; i < num_components; i ++ ) {
            num_atoms += inp_cur_data[i].num_at;
        }
        at = (inp_ATOM *) inchi_calloc( num_atoms, sizeof(at[0]) );
        cur_num_at = 0;
        for ( i = 0; i < num_components; i ++ ) {
            /* clean and prepare */
            UnMarkRingSystemsInp( inp_cur_data[i].at, inp_cur_data[i].num_at );
            UnMarkOtherIndicators( inp_cur_data[i].at, inp_cur_data[i].num_at );
            UnMarkOneComponent( inp_cur_data[i].at, inp_cur_data[i].num_at );
            subtract_DT_from_num_H( inp_cur_data[i].num_at, inp_cur_data[i].at );
            /* merge one by one */
            cur_num_at = add_inp_ATOM( at, num_atoms, cur_num_at, inp_cur_data[i].at, inp_cur_data[i].num_at );
        }
        /* replace orig_inp_data */
        if ( cur_num_at == num_atoms ) {
            inchi_free( orig_inp_data->at );
            orig_inp_data->at = at;
            orig_inp_data->num_inp_atoms = cur_num_at;
            if ( orig_inp_data->szCoord ) {
                inchi_free( orig_inp_data->szCoord );
                orig_inp_data->szCoord = NULL;
            }
            UnMarkDisconnectedComponents( orig_inp_data );
        } else {
            if ( at ) {
                inchi_free( at );
                at = NULL;
            }
            ret = -999; /* num atoms mismatch */
        }
    }
exit_function:
    if ( da ) {
        inchi_free( da );
        da = NULL;
    }
    for ( j = 0; j < num_components; j ++ ) {
        FreeInpAtomData( inp_cur_data+j );
    }
    inchi_free( inp_cur_data );
    inp_cur_data = NULL;

    return ret? ret : num_cuts;
}
#endif /* RING2CHAIN */


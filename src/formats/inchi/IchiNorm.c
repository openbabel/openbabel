/*
 * International Union of Pure and Applied Chemistry (IUPAC)
 * International Chemical Identifier (InChI)
 * Version 1
 * Software version 1.00
 * April 13, 2005
 * Developed at NIST
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
int MarkRingSystemsInp( inp_ATOM *at, int num_atoms )
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
    int        i, j, u, start, nNumRingSystems, nNumStartChildren;

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
    start           = 0;
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
        if ( (int)at[i=nStackAtom[nTopStackAtom]].valence > (j = (int)cNeighNumb[i]) ) {
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
    start           = 0;
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
        if ( (int)at[i=nStackAtom[nTopStackAtom]].valence > (j = (int)cNeighNumb[i]) ) {
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
int remove_terminal_HDT( int num_atoms, inp_ATOM *at )
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
    for ( i = 0; i < num_atoms; i ++ ) {
        at[i].component = i; /*  temporarily save original numbering */
        /*  get k = temp. hydrogen isotope/non-hydrogen atom type: */
        /*  k=0:H, k=2:D, k=3:T, k=4=kMax: not a hydrogen */
        k = at[i].elname[1]? kMax : (p=strchr(szHDT, at[i].elname[0]))? p-szHDT : kMax;
        /*  set hydrogen isotope atw differences */
        /*  Notes: k-value of isotopic H is incremented to correct iso_atw_diff value later. */
        /*         1H isotope cannot be detected here. */
        if ( k == ATW_H || k == ATW_H+1 ) {  /* D or T, k = 1 or 2 */
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
    if ( 2 == num_H && 2 == num_atoms && !NUMH(at,0) && !NUMH(at,1) ) {

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

    } else {
        /* general case except H-H */
        for ( i = 0; i < num_atoms; i ++ ) {
            k = (at[i].elname[1] || NUMH(at,i))? kMax : (at[i].elname[0]=='H')? at[i].iso_atw_diff : kMax;
            if ( k < kMax && at[i].valence == 1 && at[i].chem_bonds_valence == 1 &&
                 /*  the order of comparison is important */
                 ((n=(int)at[i].neighbor[0]) > i               /* at[n] has not been encountered yet*/ ||
                  (int)new_ord[n] < num_atoms - num_hydrogens) /* at[n] might have been encountered; it has not been moved */ ) {
                /*  found an explicit terminal hydrogen */
                num_hydrogens ++;
                if ( k==0 && ATW_H <= at[i].iso_atw_diff && at[i].iso_atw_diff < ATW_H+NUM_H_ISOTOPES ) {
                    k = at[i].iso_atw_diff; /*  H isotope has already been marked above or elsewhere */
                }
                if ( at[i].charge ) { /*  transfer charge from the hydrogen */
                    at[n].charge += at[i].charge;
                    at[i].charge = 0;
                }
                new_ord[i] = num_atoms - num_hydrogens;  /*  move hydrogens to the end of the list */
            } else {
                new_ord[i] = i - num_hydrogens;  /*  adjust non-hydrogens positions */
            }
            new_at[new_ord[i]] = at[i]; /*  copy atoms to their new positions */
        }
    }

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

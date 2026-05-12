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


/*
Underivatization, ring-chain tautomerism, OriGAtData edits, etc.
*/
#include <stdlib.h>
#include <string.h>

#include "mode.h"
#include "ichinorm.h"
#include "ichierr.h"

#include "bcf_s.h"

#if ( FIND_RING_SYSTEMS == 1 ) /* { */



/****************************************************************************/
int MarkRingSystemsInp( inp_ATOM *at, int num_atoms, int start )
{
    AT_NUMB   *nStackAtom = NULL;
    int        nTopStackAtom = -1;
    AT_NUMB   *nRingStack = NULL;
    int        nTopRingStack = -1; /* was AT_NUMB */
    AT_NUMB   *nDfsNumber = NULL;
    AT_NUMB   *nLowNumber = NULL;
    S_CHAR    *cNeighNumb = NULL;
    AT_NUMB    nDfs;
#if ( FIND_RINS_SYSTEMS_DISTANCES == 1 )
    AT_NUMB    nRs, *nRsConnect = NULL;
    int        k;
    AT_NUMB   *tree = NULL;
    int        nNumConnect, nMaxNumConnect, nLenConnect;
#endif
    AT_NUMB    nNumAtInRingSystem;
    int        i, j, u, /*start,*/ nNumRingSystems, nNumStartChildren;

    /*  allocate arrays */
    nStackAtom = (AT_NUMB *) inchi_malloc( num_atoms * sizeof( nStackAtom[0] ) );
    nRingStack = (AT_NUMB *) inchi_malloc( num_atoms * sizeof( nRingStack[0] ) );
    nDfsNumber = (AT_NUMB *) inchi_malloc( num_atoms * sizeof( nDfsNumber[0] ) );
    nLowNumber = (AT_NUMB *) inchi_malloc( num_atoms * sizeof( nLowNumber[0] ) );
    cNeighNumb = (S_CHAR  *) inchi_malloc( num_atoms * sizeof( cNeighNumb[0] ) );
#if ( FIND_RINS_SYSTEMS_DISTANCES == 1 )
    nRsConnect = (AT_NUMB *) inchi_calloc( 3 * num_atoms + 3, sizeof( nRsConnect[0] ) );
#endif
    /*  check allocation */
    if (!nStackAtom || !nRingStack || !nDfsNumber || !nLowNumber || !cNeighNumb
#if ( FIND_RINS_SYSTEMS_DISTANCES == 1 )
         || !nRsConnect
#endif
         )
    {
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
    u = start; /*  start atom */
    nDfs = 0;
    nTopStackAtom = -1;
    nTopRingStack = -1;
    memset( nDfsNumber, 0, num_atoms * sizeof( nDfsNumber[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    memset( cNeighNumb, 0, num_atoms * sizeof( cNeighNumb[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    /*  push the start atom on the stack */
    /* djb-rwth: fixing oss-fuzz issue #66720 */
    if (u <= num_atoms - 1)
    {
        nLowNumber[u] = nDfsNumber[u] = ++nDfs;
        nStackAtom[++nTopStackAtom] = (AT_NUMB)u;
        nRingStack[++nTopRingStack] = (AT_NUMB)u;
    }
    else
    {
        nNumRingSystems = CT_OVERFLOW;  /*  program error */ /*   <BRKPT> */
        goto exit_function;
    }

    nNumStartChildren = 0;

    do
    {

        /* advance */
    advance_block:

        /*if ( (int)at[i=nStackAtom[nTopStackAtom]].valence > (j = (int)cNeighNumb[i]) )*/
        /* replaced due to missing sequence point */
        if (i = (int) nStackAtom[nTopStackAtom], j = (int) cNeighNumb[i], (int) at[i].valence > j)
        {
            cNeighNumb[i] ++;
            u = (int) at[i].neighbor[j];
            if (!nDfsNumber[u])
            {
                /* tree edge, 1st visit -- advance */
                nStackAtom[++nTopStackAtom] = (AT_NUMB) u;
                nRingStack[++nTopRingStack] = (AT_NUMB) u;
                nLowNumber[u] = nDfsNumber[u] = ++nDfs;
                nNumStartChildren += ( i == start );
            }
            else
            {
                if (!nTopStackAtom || u != (int) nStackAtom[nTopStackAtom - 1])
                {
                    /*  may comment out ? */
                    /* back edge: u is not a predecessor of i */
                    if (nDfsNumber[u] < nDfsNumber[i])
                    {
                        /* Back edge, 1st visit: u is an ancestor of i. Compare */
                        if (nLowNumber[i] > nDfsNumber[u])
                        {
                            nLowNumber[i] = nDfsNumber[u];
                        }
                    }
                } /*  may comment out ? */
            }
            goto advance_block;
        }
        else
        {
            cNeighNumb[i] = 0;
        }

        /* back up */
        if (i != start)
        {
            u = (int) nStackAtom[nTopStackAtom - 1]; /* predecessor of i */
            if (nLowNumber[i] >= nDfsNumber[u])
            {
                /* output the block */
                nNumRingSystems++;
                at[u].nBlockSystem = nNumRingSystems;
                if (u != start || nNumStartChildren > 1)
                {
                    at[u].bCutVertex += 1;
                }
                while (nTopRingStack >= 0)
                {
                    j = nRingStack[nTopRingStack--];
                    at[j].nBlockSystem = nNumRingSystems; /*  mark the atom */
                    if (i == j)
                    {
                        break;
                    }
                }
            }
            else
            {
                if (nLowNumber[u] > nLowNumber[i])
                {
                    /* inherit */
                    nLowNumber[u] = nLowNumber[i];
                }
            }
        }
    } while (--nTopStackAtom >= 0);

    /****************************************************************************
    *
    * Find Ring Systems
    * Including chain atoms X: A-X-B, where the bonds (of any kind) are bridges.
    *
    ****************************************************************************/


    /*  initiation */
    /* start           = 0;*/
    nNumRingSystems = 0;
    u = start; /*  start atom */
    nDfs = 0;
    nTopStackAtom = -1;
    nTopRingStack = -1;
    memset( nDfsNumber, 0, num_atoms * sizeof( nDfsNumber[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    memset( cNeighNumb, 0, num_atoms * sizeof( cNeighNumb[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    /*  push the start atom on the stack */
    nLowNumber[u] = nDfsNumber[u] = ++nDfs;
    nStackAtom[++nTopStackAtom] = (AT_NUMB) u;
    nRingStack[++nTopRingStack] = (AT_NUMB) u;

#if ( FIND_RINS_SYSTEMS_DISTANCES == 1 )
    nNumConnect = nLenConnect = nMaxNumConnect = 0;
#endif

    do
    {
        /* advance */
    advance_ring:
        /*if ( (int)at[i=nStackAtom[nTopStackAtom]].valence > (j = (int)cNeighNumb[i]) )*/
        /* replaced due to missing sequence point */
        if (i = (int) nStackAtom[nTopStackAtom], j = (int) cNeighNumb[i], (int) at[i].valence > j)
        {
            cNeighNumb[i] ++;
            u = (int) at[i].neighbor[j];
            if (!nDfsNumber[u])
            {
                /* tree edge, 1st visit -- advance */
                nStackAtom[++nTopStackAtom] = (AT_NUMB) u;
                nRingStack[++nTopRingStack] = (AT_NUMB) u;
                nLowNumber[u] = nDfsNumber[u] = ++nDfs;
            }
            else
            {
                if (!nTopStackAtom || u != (int) nStackAtom[nTopStackAtom - 1])
                {
                    /* back edge: u is not a predecessor of i */
                    if (nDfsNumber[u] < nDfsNumber[i])
                    {
                        /* Back edge, 1st visit: u is ancestor of i. Compare */
                        if (nLowNumber[i] > nDfsNumber[u])
                        {
                            nLowNumber[i] = nDfsNumber[u];
                        }
                    }
                }
            }
            goto advance_ring;
        }
        else
        {
            cNeighNumb[i] = 0;
        }

        /* back up */
        if (nDfsNumber[i] == nLowNumber[i])
        {
            /*  found a ring system */
            nNumRingSystems++;
            /*  unwind nRingStack[] down to i */
#if ( FIND_RINS_SYSTEMS_DISTANCES == 1 )
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
            for (nNumAtInRingSystem = 0, j = nTopRingStack; 0 <= j; j--)
            {
                nNumAtInRingSystem++;
                if (i == (int) nRingStack[j])
                {
                    break;
                }
            }
            while (nTopRingStack >= 0)
            {
                j = (int) nRingStack[nTopRingStack--];
                at[j].nRingSystem = (AT_NUMB) nNumRingSystems; /*  ring system id */
                at[j].nNumAtInRingSystem = nNumAtInRingSystem;
#if ( FIND_RINS_SYSTEMS_DISTANCES == 1 )
                for (k = 0; k < at[j].valence; k++)
                {
                    if (( nRs = at[at[j].neighbor[k]].nRingSystem ) && (int) nRs != nNumRingSystems)
                    {
                        nRsConnect[nLenConnect + ( nNumConnect++ )] = nRs; /*  adjacent ring system id */
                    }
                }
#endif
                if (i == j)
                {
                    /*  reached atom on the top of nStackAtom[] stack  */
                    break;
                }
            }
#if ( FIND_RINS_SYSTEMS_DISTANCES == 1 )
            nRsConnect[nLenConnect] = nNumConnect;
            nRsConnect[nLenConnect + 1] = nNumAtInRingSystem;
            nLenConnect += nNumConnect;
            if (nMaxNumConnect < nNumConnect)
            {
                /*  max number of neighboring ring systems */
                nMaxNumConnect = nNumConnect;
            }
#endif
        }
        else
        {
            if (nTopStackAtom > 0)
            {
                j = (int) nStackAtom[nTopStackAtom - 1];
                /* inherit nLowNumber */
                if (nLowNumber[j] > nLowNumber[i])
                {
                    nLowNumber[j] = nLowNumber[i];
                }
            }
        }
    } while (--nTopStackAtom >= 0);

#if ( FIND_RINS_SYSTEMS_DISTANCES == 1 ) /*  normally disabled */
    nMaxNumConnect++;
    if (nNumRingSystems > 1)
    {
        int nCol = nMaxNumConnect + 1;
        int nNumInSyst = nMaxNumConnect;
        int nMaxNeigh = nMaxNumConnect - 1;
#define T(a,b) tree[(a)*nCol+b]
        if (tree = (AT_NUMB *) inchi_calloc( nCol * ( nNumRingSystems + 1 ), sizeof( tree[0] ) ))
        {
            int len, neigh;
            /*  reuse previous allocations */
            AT_NUMB *nNumVisitedNeighbors = nStackAtom;
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
            for (i = 1, j = 0; len = nRsConnect[j]; i++)
            {
                T( i, nNumInSyst ) = nRsConnect[j + 1];
                for (k = 2; k < len; k++)
                {
                    neigh = nRsConnect[j + k];
                    if (T( i, 0 ) < nMaxNeigh && T( neigh, 0 ) < nMaxNeigh)
                    {
                        T( i, 0 )++;
                        T( neigh, 0 )++;
                        T( i, T( i, 0 ) ) = neigh;
                        T( neigh, T( neigh, 0 ) ) = i;
                    }
                    else
                    {
                        nNumRingSystems = CT_OVERFLOW;  /*  program error */ /*   <BRKPT> */
                        goto exit_function;
                    }
                }
                j += len;
            }
            /*  clear memory */
            memset( nNumVisitedNeighbors, 0, nNumRingSystems * sizeof( nNumVisitedNeighbors[0] ) );
            memset( nDistanceFromTerminal, 0, nNumRingSystems * sizeof( nDistanceFromTerminal[0] ) );
            memset( nCurrActiveRingSystem, 0, nNumRingSystems * sizeof( nCurrActiveRingSystem[0] ) );
            memset( nNextActiveRingSystem, 0, nNumRingSystems * sizeof( nNextActiveRingSystem[0] ) );
            nNumNextActiveRingSystems = 0;
            for (i = 0; i < nNumRingSystems; i++)
            {
                if (1 == T( i + 1, 0 ))
                {
                    nNextActiveRingSystem[i] = 1; /*  number of traversed neighbors + 1 */
                    nDistanceFromTerminal[i] = 1;
                    nNumNextActiveRingSystems++;
                }
                else
                {
                    nNextActiveRingSystem[i] = 0;
                    nDistanceFromTerminal[i] = 0;
                }
                nNumVisitedNeighbors[i] = 0;
            }

            /* nCurrActiveRingSystem[i] = a sum of:
            * 1) +1 if it is or was active
            * 2) +(number of neighbors from which it was reached)
            * 3) +1 if it was left and not active anymore
            */
            pass = 0;
            do
            {
                nNumCurrActiveRingSystems = nNumNextActiveRingSystems;
                nNumNextActiveRingSystems = 0;
                memcpy( nCurrActiveRingSystem, nNextActiveRingSystem,
                        nNumRingSystems * sizeof( nNextActiveRingSystem[0] ) );
                for (i = 0; i < nNumRingSystems; i++)
                {
                    if (T( i + 1, 0 ) == nCurrActiveRingSystem[i])
                    {
                        /* on the previous pass currently active ring system i+1 bas been reached
                        * from all neighbors except one;
                        * the neighbors from which it was reached have
                        * T(neigh,0)+1 == nCurrActiveRingSystem[i]
                        * this ring system has not been left yet
                        */
                        for (k = 1, len = T( i + 1, 0 ); k <= len; k++)
                        {
                            neigh = (int) T( i + 1, k );
                            if (T( neigh, 0 ) >= nCurrActiveRingSystem[neigh - 1])
                            {
                                if (0 == pass)
                                {
                                    nDistanceFromTerminal[i] = 1;
                                }
                                break;
                            }
                        }
                        if (k <= len)
                        {
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
                            k = i + 1; /*  starting point */
                            if (0 == pass && T( k, nNumInSyst ) > 1)
                            {
                                nNumNextActiveRingSystems++; /*  request next pass */
                                continue; /*  stop a the terminal ring system */
                            }
                            while (2 == T( neigh, 0 ))
                            {
                                /*  walk along a chain */
                                if (!nNextActiveRingSystem[neigh - 1])
                                {
                                    nNextActiveRingSystem[neigh - 1] = 1; /*  make neighbor active */
                                }
                                else
                                    if (nDistanceFromTerminal[k - 1] + 1 <= nDistanceFromTerminal[neigh - 1])
                                    {
                                        /*  walking along the chain; already have had a walk */
                                        /*  in the opposite direction at this pass */
                                    }
                                    else
                                    {
                                        /*  k is the last; neigh (it is a bridge -X-) has not been reached */
                                        bOk = 1;
                                        break;
                                    }
                                nNextActiveRingSystem[k - 1] ++; /*  leave system k */
                                if (nNextActiveRingSystem[neigh - 1] < T( neigh, 0 ))
                                {
                                    nNextActiveRingSystem[neigh - 1] ++; /*  add one connection to neigh */
                                }
                                nDistanceFromTerminal[neigh - 1] = nDistanceFromTerminal[k - 1] + 1;
                                j = ( T( neigh, 1 ) == k ) ? 2 : 1;
                                k = neigh;
                                neigh = T( k, j ); /*  next in the chain */
                                nNumNextActiveRingSystems++;
                                if (T( k, nNumInSyst ) > 1)
                                {
                                    bOk = 1;
                                    break; /*  stop on a ring system */
                                }
                            }
                            /*  neigh is a terminal or a bridge or a branching point */
                            if (2 > T( neigh, 0 ))
                            {
                                /*  neighbor is a terminal atom */
                                if (1 < pass)
                                {
                                    nNumRingSystems = CT_UNKNOWN_ERR; /*  error (debug only) */ /*   <BRKPT> */
                                    goto exit_function;
                                }
                                continue;
                            }
                            if (2 == T( neigh, 0 ))
                            {
                                /*  neighbor is a bridge */
                                continue;
                            }
                            /*  neighbor is a branching point */
                            if (T( neigh, 0 ) > nCurrActiveRingSystem[neigh - 1])
                            {
                                /*  move to the neigh (make neigh active): on previous pass it */
                                /*  has not been reached from 2 or more neighbors */
                                if (!nNextActiveRingSystem[neigh - 1])
                                {
                                    nNextActiveRingSystem[neigh - 1] = 1;
                                }
                                if (nDistanceFromTerminal[neigh - 1] < nDistanceFromTerminal[k - 1] + 1)
                                {
                                    nDistanceFromTerminal[neigh - 1] = nDistanceFromTerminal[k - 1] + 1;
                                }
                                nNextActiveRingSystem[k - 1] ++; /*  leave system k */
                                if (nNextActiveRingSystem[neigh - 1] < T( neigh, 0 ))
                                {
                                    nNextActiveRingSystem[neigh - 1] ++; /*  add one connection to neigh */
                                }
                                nNumNextActiveRingSystems++;
                            }
                        }
                    }
                }
                pass++;
            } while (nNumNextActiveRingSystems);

            for (i = 0; i < num_atoms; i++)
            {
                at[i].nDistanceFromTerminal = nDistanceFromTerminal[(int) at[i].nRingSystem - 1];
            }

            inchi_free( tree );
            tree = NULL;
#undef T
        }
        else
        {
            nNumRingSystems = CT_OUT_OF_RAM; /*  error */ /*   <BRKPT> */
            goto exit_function;
        }
    }
#endif

exit_function:

    if (nStackAtom)
    {
        inchi_free( nStackAtom );
    }
    if (nRingStack)
    {
        inchi_free( nRingStack );
    }
    if (nDfsNumber)
    {
        inchi_free( nDfsNumber );
    }
    if (nLowNumber)
    {
        inchi_free( nLowNumber );
    }
    if (cNeighNumb)
    {
        inchi_free( cNeighNumb );
    }

#if ( FIND_RINS_SYSTEMS_DISTANCES == 1 )
    if (nRsConnect)
        inchi_free( nRsConnect );
    if (tree)
        inchi_free( tree );
#endif

    return nNumRingSystems;
}


#endif /* } FIND_RING_SYSTEMS */


/****************************************************************************

InChI post-version 1.01 features implementation
(v. 1.06+ : underivatize is still an experiment available in engineering mode)

****************************************************************************/

#if ( RING2CHAIN == 1 || UNDERIVATIZE == 1 )

typedef struct tagAtPair
{
    AT_NUMB at[2];  /* at[0] < at[1] */
    AT_NUMB atno;   /* atom marked with derivative type */
} R2C_ATPAIR;


/* Local functions */
int mark_arom_bonds( struct tagINCHI_CLOCK *ic,
                     struct tagCANON_GLOBALS *pCG,
                     inp_ATOM *at, int num_atoms );
void set_R2C_el_numbers( void );
int subtract_DT_from_num_H( int num_atoms, inp_ATOM *at );
int add_inp_ATOM( inp_ATOM *at, int len_at, int len_cur,
                  inp_ATOM *add, int len_add );
int cmp_r2c_atpair( const void *p1, const void *p2 );
int has_atom_pair( R2C_ATPAIR *ap, int num_ap, AT_NUMB at1, AT_NUMB at2 );
int mark_atoms_ap( inp_ATOM *at, AT_NUMB start, R2C_ATPAIR *ap,
                   int num_ap, int num, AT_NUMB cFlags );

int UnMarkDisconnectedComponents( ORIG_ATOM_DATA *orig_inp_data );
int UnMarkOtherIndicators( inp_ATOM *at, int num_atoms );
int UnMarkOneComponent( inp_ATOM *at, int num_atoms );

/* Other functions */
int DisconnectInpAtBond( inp_ATOM *at, AT_NUMB *nOldCompNumber,
                         int iat, int neigh_ord );
int ExtractConnectedComponent( inp_ATOM *at, int num_at,
                               int component_number, inp_ATOM *component_at );
int UnMarkRingSystemsInp( inp_ATOM *at, int num_atoms );



/****************************************************************************
Clear the (disconnected) components info in  ORIG_ATOM_DATA
****************************************************************************/
int UnMarkDisconnectedComponents( ORIG_ATOM_DATA *orig_inp_data )
{
    int i;

    for (i = 0; i < orig_inp_data->num_inp_atoms; i++)
    {
        orig_inp_data->at[i].orig_compt_at_numb = 0;
        orig_inp_data->at[i].component = 0;
    }

    if (orig_inp_data->nCurAtLen)
    {
        inchi_free( orig_inp_data->nCurAtLen );
        orig_inp_data->nCurAtLen = NULL;
    }

    if (orig_inp_data->nOldCompNumber)
    {
        inchi_free( orig_inp_data->nOldCompNumber );
        orig_inp_data->nOldCompNumber = NULL;
    }

    orig_inp_data->num_components = 0;

    return 0;
}


/****************************************************************************/
int UnMarkOtherIndicators( inp_ATOM *at, int num_atoms )
{
    int i;
    for (i = 0; i < num_atoms; i++)
    {
        at[i].at_type = 0;
        at[i].cFlags = 0;
    }

    return 0;
}


/****************************************************************************
Clear the (disconnected) component numbers in atoms of inp_ATOM structure
(which typicall came from INP_ATOM_DATA)
****************************************************************************/
int UnMarkOneComponent( inp_ATOM *at, int num_atoms )
{
    int i;
    for (i = 0; i < num_atoms; i++)
    {
        at[i].orig_compt_at_numb = 0;
        at[i].component = 0;
    }

    return 0;
}


/****************************************************************************/
void set_R2C_el_numbers( void )
{
    /*
    if (!el_number_O)
    {
    el_number_O = EL_NUMBER_O;
    el_number_C = EL_NUMBER_C;
    el_number_N = EL_NUMBER_N;
    el_number_P = EL_NUMBER_P;
    el_number_S = EL_NUMBER_S;
    el_number_Si = EL_NUMBER_SI;
    el_number_F = EL_NUMBER_F;
    el_number_Cl = EL_NUMBER_CL;
    el_number_Br = EL_NUMBER_BR;
    el_number_I = EL_NUMBER_I;
    el_number_H = EL_NUMBER_H;
    }
    */
}


/****************************************************************************/
int subtract_DT_from_num_H( int num_atoms, inp_ATOM *at )
/*  assume num_1H, num_D and num_T are included in num_H */
{
    int i, j;
    for (i = 0; i < num_atoms; i++)
    {
        for (j = 0; j < NUM_H_ISOTOPES; j++)
            at[i].num_H -= at[i].num_iso_H[j];
    }

    return 0;
}


/****************************************************************************/
int add_inp_ATOM( inp_ATOM *at,
                  int len_at,
                  int len_cur,
                  inp_ATOM *add,
                  int len_add )
{
    int i, j;
    inp_ATOM *a;
    /* chack correctness */
    if (len_cur < 0)
        return len_cur;
    if (len_add < 0)
        return len_add;
    if (len_cur + len_add > len_at)
        return -1;
    /* copy */
    memcpy(at + len_cur, add, len_add * sizeof(at[0]));
    /* modify */
    if (len_cur)
    {
        a = at + len_cur;
        for (i = 0; i < len_add; i++)
        {
            for (j = 0; j < a[i].valence; j++)
            {
                a[i].neighbor[j] += len_cur;
            }
        }
    }

    return len_cur + len_add;
}


/****************************************************************************/
int mark_arom_bonds( struct tagINCHI_CLOCK *ic, struct tagCANON_GLOBALS *pCG, inp_ATOM *at, int num_atoms )
{
    INCHI_MODE bTautFlags = 0, bTautFlagsDone = 0;
    inp_ATOM *at_fixed_bonds_out = NULL;
    T_GROUP_INFO *t_group_info = NULL;
    int ret;

    ret = mark_alt_bonds_and_taut_groups( ic, pCG, at, at_fixed_bonds_out, num_atoms,
                                          NULL,
                                          t_group_info, &bTautFlags, &bTautFlagsDone, 0, NULL );

    return ret;
}


/****************************************************************************/
int cmp_r2c_atpair( const void *p1, const void *p2 )
{
    const R2C_ATPAIR *ap1 = (const R2C_ATPAIR *) p1;
    const R2C_ATPAIR *ap2 = (const R2C_ATPAIR *) p2;
    int diff = (int) ap1->at[0] - (int) ap2->at[0];
    if (!diff)
    {
        diff = (int) ap1->at[1] - (int) ap2->at[1];
    }

    return diff;
}


/****************************************************************************/
int has_atom_pair_seq( R2C_ATPAIR *ap, int num_ap, AT_NUMB at1, AT_NUMB at2 )
{
    R2C_ATPAIR ap1;
    int i1;
    int n = at1 > at2;

    ap1.at[n] = at1;
    ap1.at[1 - n] = at2;
    for (i1 = 0; i1 < num_ap; i1++)
    {
        if (ap[i1].at[0] == ap1.at[0] &&
             ap[i1].at[1] == ap1.at[1])
            return i1 + 1;
    }

    return 0; /* not found */
}


/****************************************************************************/
int has_atom_pair( R2C_ATPAIR *ap, int num_ap, AT_NUMB at1, AT_NUMB at2 )
{
    R2C_ATPAIR ap1;
    int i1, i2, i3, diff;
    int n = at1 > at2;

    ap1.at[n] = at1;
    ap1.at[1 - n] = at2;
    i1 = 0;
    i2 = num_ap - 1;
    /* search for ap1 by simple bisections */
    do
    {
        i3 = ( i1 + i2 ) / 2;
        if (!( diff = cmp_r2c_atpair( &ap1, ap + i3 ) ))
        {
            return i3 + 1;  /* found => positive number */
        }
        else
            if (diff > 0)
            {
                i1 = i3 + 1;
            }
            else
            {
                i2 = i3 - 1;
            }
    } while (i2 >= i1);

    return 0; /* not found */
}


/****************************************************************************
DFS search for atoms that do not have a flag
****************************************************************************/
int mark_atoms_ap( inp_ATOM *at,
                   AT_NUMB start,
                   R2C_ATPAIR *ap,
                   int num_ap, int num,
                   AT_NUMB cFlags )
{
    if (!at[start].at_type)
    {
        int i;
        AT_NUMB neigh;
        at[start].at_type = cFlags;
        num++;

        for (i = 0; i < at[start].valence; i++)
        {
            neigh = at[start].neighbor[i];
            if (has_atom_pair_seq( ap, num_ap, start, neigh ))
            {
                continue;
            }
            num = mark_atoms_ap( at, neigh, ap, num_ap, num, cFlags );
        }
    }

    return num; /* number of atoms traversed forward from at[start] */
}

#endif /* RING2CHAIN || UNDERIVATIZE */



#if ( UNDERIVATIZE == 1 )

/****************************************************************************/

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


/* Underivatize settings */
/* DERIV_AT::typ begin */
#define DERIV_BRIDGE_O  0x0001   /* R1-O-R2 => R1-OH + HO-R2 */
#define DERIV_BRIDGE_NH 0x0002   /* R1-NH-R2  amine */
#define DERIV_AMINE_tN  0x0004   /* R1-N(-R2)-R3  tertiary amine */
#define DERIV_RING_O_OUTSIDE_PRECURSOR    0x0008   /* -O- in a ring */
#define DERIV_RING_NH_OUTSIDE_PRECURSOR   0x0010   /* -NH- in a ring */
/* MOX_EtOX Underiv: R2-(R3-)C=N-O-R => >C=O (ketone, aldehide only);
-R: -CH3, -CH2-CH3, -Si(CH3)3, -CH2-Phenyl
DERIV_X_OXIME        R3- may be H or any C, R2- is any C
precursor: R2-(R3-)C=O  (note: N replaced with O)-- 2013-08-23 DT */
#define DERIV_X_OXIME   0x0020   /* comment out to disable */
/*  */
#define DERIV_UNMARK    0x0040   /* unmark the cut */
#define DERIV_DUPLIC    0x0080   /* duplicated disconnection */
/* comment out to disable */
#define DERIV_RO_COX    0x0100   /* alcohol derivatives: R-O--C(=O)C[n]F[2n+1] 0<n<4, R-O--C(=O)CH3, R-O--C(=O)-Phenyl */
/* comment out next 2 to disable; DMOX=Dimethyloxazoline DEOX=Diethyloxazoline */
#define DERIV_RING_DMOX_DEOX_N   0x0200   /* =N- in a ring:    /-O--CH2-\      DMOX, DEOX  */
#define DERIV_RING_DMOX_DEOX_O   0x0400   /* -O- in a ring: R-C=N--------C<2Me or 2Et */

#ifdef UNDERIV_PYRROLIDIDES
#define DERIV_RING2_PRRLDD_OUTSIDE_PRECUR 0x0800  /* alcohol derivatives: R(=O)-N<C4H4 5-memb ring Pyrrolidides; replace -N< with -OH */
#endif

#define DERIV_RING2_PPRDN_OUTSIDE_PRECUR  0x1000  /* alcohol derivatives: R(=O)-N<C5H5 6-memb ring Piperidines; replace -N< with -OH */
#define DERIV_DANSYL                      0x2000  /* alcohol derivatives: R-O--SO2-C10H5-N(CH3)2 => R-OH */
/* DERIV_AT::typ end */

/* derivative types, which cannot be absorbed into greater derivatization agents */
static const int DERIV_UNEXPADABLE = 0
#ifdef DERIV_X_OXIME
| DERIV_X_OXIME
#endif
#if( defined(DERIV_RING_DMOX_DEOX_N) && defined(DERIV_RING_DMOX_DEOX_O) )
| ( DERIV_RING_DMOX_DEOX_O | DERIV_RING_DMOX_DEOX_N )
#endif
/*
#ifdef DERIV_RO_COX
|DERIV_RO_COX
#endif
*/
#ifdef DERIV_RING2_PRRLDD_OUTSIDE_PRECUR
| DERIV_RING2_PRRLDD_OUTSIDE_PRECUR
#endif
#ifdef DERIV_RING2_PPRDN_OUTSIDE_PRECUR
| DERIV_RING2_PPRDN_OUTSIDE_PRECUR
#endif
#ifdef DERIV_DANSYL
| DERIV_DANSYL
#endif
;
/* derivative precursor types in which Precur=N- is replaced with Precur=O */
static const int DERIV_REPL_N_WITH_O = 0
#ifdef DERIV_X_OXIME
| DERIV_X_OXIME
#endif
#ifdef DERIV_RING_DMOX_DEOX_N
| DERIV_RING_DMOX_DEOX_N
#endif
;

/* derivative precursor types in which Precur-N< is replaced with Precur-OH */
static const int DERIV_REPL_N_WITH_OH = 0
#ifdef DERIV_RING2_PRRLDD_OUTSIDE_PRECUR
| DERIV_RING2_PRRLDD_OUTSIDE_PRECUR
#endif
#ifdef DERIV_RING2_PPRDN_OUTSIDE_PRECUR
| DERIV_RING2_PPRDN_OUTSIDE_PRECUR
#endif
;

/* combined DERIV_AT::typ */
#if( defined(DERIV_RING_DMOX_DEOX_N) && defined(DERIV_RING_DMOX_DEOX_O) )
#define DERIV_RING_DMOX_DEOX     (DERIV_RING_DMOX_DEOX_O | DERIV_RING_DMOX_DEOX_N)
#endif

#if( defined(DERIV_RING_O_OUTSIDE_PRECURSOR) && defined(DERIV_RING_NH_OUTSIDE_PRECURSOR) )
#define DERIV_RING_OUTSIDE_PRECURSOR     (DERIV_RING_O_OUTSIDE_PRECURSOR | DERIV_RING_NH_OUTSIDE_PRECURSOR)
#elif( defined(DERIV_RING_O_OUTSIDE_PRECURSOR) )
#define DERIV_RING_OUTSIDE_PRECURSOR   DERIV_RING_O_OUTSIDE_PRECURSOR
#elif( defined(DERIV_RING_NH_OUTSIDE_PRECURSOR) )
#define DERIV_RING_OUTSIDE_PRECURSOR   DERIV_RING_NH_OUTSIDE_PRECURSOR
#else
#define DERIV_RING_OUTSIDE_PRECURSOR 0
#endif

#if( defined(DERIV_RING2_PRRLDD_OUTSIDE_PRECUR) && defined(DERIV_RING2_PPRDN_OUTSIDE_PRECUR) )
#define DERIV_RING2_OUTSIDE_PRECUR  (DERIV_RING2_PRRLDD_OUTSIDE_PRECUR | DERIV_RING2_PPRDN_OUTSIDE_PRECUR)
#elif ( defined(DERIV_RING2_PRRLDD_OUTSIDE_PRECUR) )
#define DERIV_RING2_OUTSIDE_PRECUR  DERIV_RING2_PRRLDD_OUTSIDE_PRECUR
#elif ( defined(DERIV_RING2_PPRDN_OUTSIDE_PRECUR) )
#define DERIV_RING2_OUTSIDE_PRECUR DERIV_RING2_PPRDN_OUTSIDE_PRECUR
#endif
/* combined DERIV_AT::typ end */

#define DERIV_AT_LEN  4
typedef struct tagDerivAttachment {
    short typ[DERIV_AT_LEN]; /* changed from char to short on 2013-11-13 DT */
    char  ord[DERIV_AT_LEN]; /* ring: neighbor in precursor; chain: neighbor in derivatizing agent */
    char  num[DERIV_AT_LEN]; /* num. atoms to remove */
#ifdef DERIV_RING_DMOX_DEOX
    AT_NUMB other_atom;      /* other atno+1; for DERIV_RING_DMOX_DEOX */
#endif
} DERIV_AT;

/* return value */
#define DERIV_NOT       0x1000   /* cannot be a derivatization agent atom */

#define MAX_AT_DERIV      13              /* max. num of heavy atoms in removed derivatizing agent: RO-C(O)PheF5 */
#define NOT_AT_DERIV      99              /* DERIV_AT::num; rejected as > MAX_AT_DERIV */
#define MIN_AT_LEFT_DERIV 2               /* was 3 before 2013-11-12; min num heavy atoms in derivative precursor = MIN_AT_LEFT_DERIV-1 */

#define NO_ORD_VAL        0x0037          /* not used */

#define CFLAG_MARK_BRANCH      1          /* for main derivative traversal */
#define CFLAG_MARK_BLOCK       2          /* for block detection */
#define CFLAG_MARK_BLOCK_INV   ((char)~(CFLAG_MARK_BLOCK)) /* for block detection */
#define COUNT_ALL_NOT_DERIV    1      /* 1=> count ALL atoms that are not in deriv. agents */
/* 0=> only atoms that are not in DERIV_RING_OUTSIDE_PRECURSOR */
#define IS_DA_NUM_LE(DA, I, MX)  (((DA)->typ[I] && ((DA)->typ[I] & DERIV_UNEXPADABLE) == (DA)->typ[I]) || (DA)->num[I] <= (MX))


/* derivative classes */
typedef enum tagDerivId {
    DERIV_ID_Acentonate,
    DERIV_ID_Benzlidene,
    DERIV_ID_BenzOX,
    DERIV_ID_BuBorate,
    DERIV_ID_Dansyl,
    DERIV_ID_DEOX,
    DERIV_ID_DMOX,
    DERIV_ID_EtBorate,
    DERIV_ID_EtOX,
    DERIV_ID_HFB,
    DERIV_ID_MeBorate,
    DERIV_ID_MOX,
    DERIV_ID_PFB,
    DERIV_ID_PFP,
    DERIV_ID_Piperidine,
    DERIV_ID_Pyrrolidide,
    DERIV_ID_TBDMS,
    DERIV_ID_TFA,
    DERIV_ID_TMS,
    DERIV_ID_Unknown,


#if defined(UNDERIV_RN_AcMe) || defined(UNDERIV_RNH_AcMe) || defined(UNDERIV_RO_COX_Me)
    DERIV_ID_Acetate,
#endif
#if defined(UNDERIV_RO_COX_BENZOATES)
    DERIV_ID_Benzoate,
#endif


#ifdef UNDERIV_ACETATE_Me
    DERIV_ID_Methylation,
#endif
#ifdef UNDERIV_ACETATE_Et
    DERIV_ID_Ethylation,
#endif
#if defined(UNDERIV_RN_AcEt) || defined(UNDERIV_RNH_AcEt) || defined(UNDERIV_RO_COX_Et)
    DERIV_ID_Propanoate,
#endif


} DerivId;

#ifdef PTR_TO_DERIV_NAMES
extern const char **pszDerivName;
#else
const char *pszDerivName[] = {
    "Acentonate",      /*  DERIV_ID_Acentonate, 00 */
    "Benzlidene",      /*  DERIV_ID_Benzlidene, 01 */
    "BenzOX",          /*  DERIV_ID_BenzOX,     02 */
    "BuBorate",        /*  DERIV_ID_BuBorate,   03 */
    "Dansyl",          /*  DERIV_ID_Dansyl,     04 */
    "DEOX",            /*  DERIV_ID_DEOX,       05 */
    "DMOX",            /*  DERIV_ID_DMOX,       06 */
    "EtBorate",        /*  DERIV_ID_EtBorate,   07 */
    "EtOX",            /*  DERIV_ID_EtOX,       08 */
    "HFB",             /*  DERIV_ID_HFB,        09 */
    "MeBorate",        /*  DERIV_ID_MeBorate,   10 */
    "MOX",             /*  DERIV_ID_MOX,        11 */
    "PFB",             /*  DERIV_ID_PFB,        12 */
    "PFP",             /*  DERIV_ID_PFP,        13 */
    "Piperidine",      /*  DERIV_ID_Piperidine, 14 */
    "Pyrrolidide",     /*  DERIV_ID_Pyrrolidide,15 */
    "TBDMS",           /*  DERIV_ID_TBDMS,      16 */
    "TFA",             /*  DERIV_ID_TFA,        17 */
    "TMS",             /*  DERIV_ID_TMS,        18 */
    "???",             /*  DERIV_ID_Unknown,    19 */


#if defined(UNDERIV_RN_AcMe) || defined(UNDERIV_RNH_AcMe) || defined(UNDERIV_RO_COX_Me)
    "Acetate",         /* DERIV_ID_Acetate,     20 */
#endif
#if defined(UNDERIV_RO_COX_BENZOATES)
    "Benzoate",        /* DERIV_ID_Benzoate,    21 */
#endif

                       /************* disabled ***********************/

#ifdef UNDERIV_ACETATE_Me
                       "Methylation",     /* DERIV_ID_Methylation   */
#endif
#ifdef UNDERIV_ACETATE_Et
                       "Ethylation",      /* DERIV_ID_Ethylation    */
#endif
#if defined(UNDERIV_RN_AcEt) || defined(UNDERIV_RNH_AcEt) || defined(UNDERIV_RO_COX_Et)
                       "Propanoate",      /* DERIV_ID_Propanoate    */
#endif

                       "",  /* for linear searching */
};
#endif /* PTR_TO_DERIV_NAMES */

typedef enum tagDerivBit {
    DERIV_BIT_Acentonate = 1 << DERIV_ID_Acentonate,
    DERIV_BIT_Benzlidene = 1 << DERIV_ID_Benzlidene,
    DERIV_BIT_BenzOX = 1 << DERIV_ID_BenzOX,
    DERIV_BIT_BuBorate = 1 << DERIV_ID_BuBorate,
    DERIV_BIT_Dansyl = 1 << DERIV_ID_Dansyl,
    DERIV_BIT_DEOX = 1 << DERIV_ID_DEOX,
    DERIV_BIT_DMOX = 1 << DERIV_ID_DMOX,
    DERIV_BIT_EtBorate = 1 << DERIV_ID_EtBorate,
    DERIV_BIT_EtOX = 1 << DERIV_ID_EtOX,
    DERIV_BIT_HFB = 1 << DERIV_ID_HFB,
    DERIV_BIT_MeBorate = 1 << DERIV_ID_MeBorate,
    DERIV_BIT_MOX = 1 << DERIV_ID_MOX,
    DERIV_BIT_PFB = 1 << DERIV_ID_PFB,
    DERIV_BIT_PFP = 1 << DERIV_ID_PFP,
    DERIV_BIT_Piperidine = 1 << DERIV_ID_Piperidine,
    DERIV_BIT_Pyrrolidide = 1 << DERIV_ID_Pyrrolidide,
    DERIV_BIT_TBDMS = 1 << DERIV_ID_TBDMS,
    DERIV_BIT_TFA = 1 << DERIV_ID_TFA,
    DERIV_BIT_TMS = 1 << DERIV_ID_TMS,
    DERIV_BIT_Unknown = 1 << DERIV_ID_Unknown,

#if defined(UNDERIV_RN_AcMe) || defined(UNDERIV_RNH_AcMe) || defined(UNDERIV_RO_COX_Me)
    DERIV_BIT_Acetate = 1 << DERIV_ID_Acetate,
#endif
#if defined(UNDERIV_RO_COX_BENZOATES)
    DERIV_BIT_Benzoate = 1 << DERIV_ID_Benzoate,
#endif


#ifdef UNDERIV_ACETATE_Me
    DERIV_BIT_Methylation = 1 << DERIV_ID_Methylation,
#endif
#ifdef UNDERIV_ACETATE_Et
    DERIV_BIT_Ethylation = 1 << DERIV_ID_Ethylation,
#endif
#if defined(UNDERIV_RN_AcEt) || defined(UNDERIV_RNH_AcEt) || defined(UNDERIV_RO_COX_Et)
    DERIV_BIT_Propanoate = 1 << DERIV_ID_Propanoate,
#endif


    DERIV_ID_NUMBER,
} DerivBit;

typedef int BIT_UNDERIV;

int mark_atoms_cFlags( inp_ATOM *at,
                       int start,
                       int num,
                       char cFlags );
int unmark_atoms_cFlags( inp_ATOM *at,
                         int start,
                         int num,
                         char cFlags,
                         char cInvFlags );

int is_C_or_S_DB_O( inp_ATOM *at, int i );
int is_C_DB_O( inp_ATOM *at, int i );
int is_C_unsat_not_arom( inp_ATOM *at, int i );
int is_Aryl( inp_ATOM *at, int outside_point, int attachment_pont );
int is_Saturated_C( inp_ATOM *at, int attachment_pont );
int is_C_Alk( inp_ATOM *at, int i, char cFlags );

#ifdef DERIV_X_OXIME
int is_Phenyl( inp_ATOM *at,
               int outside_point,
               int attachment_point );
int is_PentaFluoroPhenyl( inp_ATOM *at,
                          int outside_point,
                          int attachment_point );
#endif

int is_Methyl( inp_ATOM *at,
               int attachment_point );
int is_Ethyl( inp_ATOM *at,
              int outside_point,
              int attachment_point );
int is_Methyl_or_Etyl( inp_ATOM *at,
                       int outside_point,
                       int attachment_point );
int is_Si_IV( inp_ATOM *at, int i );
int is_P_TB_N( inp_ATOM *at, int i );

#if( defined(DERIV_RING2_PRRLDD_OUTSIDE_PRECUR) || defined(DERIV_RING2_PPRDN_OUTSIDE_PRECUR) )
int is_DERIV_RING2_PRRLDD_PPRDN( inp_ATOM *at,
                                 int cur_atom,
                                 int from_ord,
                                 DERIV_AT *da,
                                 DERIV_AT *da1 );
#endif
#ifdef DERIV_DANSYL
int check_arom_chain( inp_ATOM *at,
                      int first,
                      int first_from,
                      int last,
                      int len );
#endif
int is_Dansyl( inp_ATOM *at,
               int cur_atom,
               int from_ord,
               DERIV_AT *da,
               DERIV_AT *da1 );
int is_possibly_deriv_neigh( inp_ATOM *at,
                             int iat,
                             int iord,
                             int type,
                             char cFlags );
int get_traversed_deriv_type( inp_ATOM *at,
                              DERIV_AT *da,
                              int k, DERIV_AT *da1,
                              char cFlags );
int add_to_da( DERIV_AT *da, DERIV_AT *add );
int mark_atoms_deriv( inp_ATOM *at,
                      DERIV_AT *da,
                      int start,
                      int num,
                      char cFlags,
                      int *pbFound );
int count_one_bond_atoms( inp_ATOM *at,
                          DERIV_AT *da,
                          int start,
                          int ord,
                          char cFlags,
                          int *bFound );
int is_silyl( inp_ATOM *at, int start, int ord_prev );
int is_silyl2( inp_ATOM *at, int start, int from_at );
int is_Me_or_Et( inp_ATOM *at, int start, int ord_prev );
int is_nButyl( inp_ATOM *at, int start, int ord_prev );
int is_CF3_or_linC3F7a( inp_ATOM *at, int start, int iat_prev );
int is_CF3_or_linC3F7( inp_ATOM *at, int start, int ord_prev );
int is_phenyl( inp_ATOM *at, int start, int ord_prev );
int is_DERIV_RING_O_or_NH_OUTSIDE_PRECURSOR( inp_ATOM *at,
                                             int start,
                                             int num_atoms,
                                             DERIV_AT *da1,
                                             int idrv,
                                             char *szUnderiv,
                                             int lenUnderiv,
                                             char *szUnderiv2,
                                             int lenUnderiv2,
                                             BIT_UNDERIV *bitUnderiv );
int is_deriv_chain( inp_ATOM *at,
                    int start,
                    int num_atoms,
                    DERIV_AT *da1,
                    int idrv,
                    char *szUnderiv,
                    int lenUnderiv,
                    char *szUnderiv2,
                    int lenUnderiv2,
                    BIT_UNDERIV *bitUnderiv );
int is_deriv_chain2( inp_ATOM *at,
                     int start,
                     int type,
                     int num,
                     int ord,
                     int idrv,
                     char *szUnderiv,
                     int lenUnderiv,
                     char *szUnderiv2,
                     int lenUnderiv2,
                     BIT_UNDERIV *bitUnderiv );
int is_deriv_chain_or_ring( inp_ATOM *at,
                            int start,
                            int num_atoms,
                            DERIV_AT *da1,
                            int *idrv );

int remove_deriv( DERIV_AT *da1, int idrv );
int remove_deriv_mark( DERIV_AT *da1, int idrv );
int underiv_compare( const void *p1, const void *p2 );
int underiv_list_add_two_cuts( char *szUnderivList,
                               int lenUnderivList,
                               char *szUnderiv,
                               const char cDelim );
int sort_merge_underiv( char *pSdfValue,
                        int bOutputSdf,
                        char *szUnderivList,
                        char cDerivSeparator,
                        const char *pszUnderivPrefix,
                        const char *pszUnderivPostfix );
int eliminate_deriv_not_in_list( inp_ATOM *at,
                                 DERIV_AT *da,
                                 int num_atoms,
                                 char *szUnderivList,
                                 int lenUnderivList,
                                 char *szUnderivList2,
                                 int lenUnderivList2,
                                 BIT_UNDERIV *bitUnderivList );
void underiv_buf_clear( char *szUnderiv );
int underiv_list_add( char *szUnderivList,
                      int lenUnderivList,
                      const char *szUnderiv,
                      char cDelimiter );
const char* underiv_list_get_last( const char *szUnderivList,
                                   char cDelimiter );
int make_single_cut( inp_ATOM *at,
                     DERIV_AT *da,
                     int iat,
                     int icut );
int fill_out_bond_cuts( inp_ATOM *at,
                        DERIV_AT *da,
                        int num_atoms,
                        R2C_ATPAIR *ap,
                        int num_cuts_to_check );
int mark_deriv_agents( inp_ATOM *at,
                       DERIV_AT *da,
                       int num_atoms,
                       R2C_ATPAIR *ap,
                       int num_cuts_to_check,
                       AT_NUMB *pnum_comp,
                       int *pcur_num_at );
int replace_arom_bonds( inp_ATOM *at,
                        int num_atoms,
                        inp_ATOM *at2,
                        int num_atoms2 );
int add_explicit_H( INP_ATOM_DATA *inp_cur_data );
void free_underiv_temp_data( R2C_ATPAIR *ap,
                             DERIV_AT *da,
                             inp_ATOM *at2,
                             INP_ATOM_DATA *inp_cur_data,
                             int num_components );
void remove_cut_derivs( int num_atoms,
                        inp_ATOM *at,
                        INP_ATOM_DATA *inp_cur_data,
                        int i_component,
                        int *errcode );


/****************************************************************************
DFS search for atoms that do not have a flag
****************************************************************************/
int mark_atoms_cFlags( inp_ATOM *at, int start, int num, char cFlags )
{
    if (!( at[start].cFlags & cFlags ))
    {
        int i;
        at[start].cFlags |= cFlags;
        num++;
        for (i = 0; i < at[start].valence; i++)
        {
            num = mark_atoms_cFlags( at, at[start].neighbor[i], num, cFlags );
        }
    }

    return num; /* number of atoms traversed forward from at[start] */
}


/****************************************************************************
DFS search for atoms that do have a flag
****************************************************************************/
int unmark_atoms_cFlags( inp_ATOM *at,
                         int start,
                         int num,
                         char cFlags,
                         char cInvFlags )
{
    if (at[start].cFlags & cFlags)
    {
        int i;
        at[start].cFlags &= cInvFlags;
        num++;
        for (i = 0; i < at[start].valence; i++)
        {
            num = unmark_atoms_cFlags( at, at[start].neighbor[i], num, cFlags, cInvFlags );
        }
    }

    return num; /* number of atoms traversed forward from at[start] */
}


/****************************************************************************/
int is_C_or_S_DB_O( inp_ATOM *at, int i )
{
    int j, neigh;
    if ((at[i].el_number != EL_NUMBER_C &&
         at[i].el_number != EL_NUMBER_S) ||
         at[i].charge || at[i].radical) /* djb-rwth: addressing LLVM warning */
        return 0;
    for (j = 0; j < at[i].valence; j++)
    {
        neigh = at[i].neighbor[j];
        if (( at[neigh].el_number == EL_NUMBER_O ||
              at[neigh].el_number == EL_NUMBER_S ) &&
             !at[neigh].num_H && 1 == at[neigh].valence &&
             2 == at[neigh].chem_bonds_valence)
        {
            return 1;
        }
    }

    return 0;
}


/****************************************************************************/
int is_C_DB_O( inp_ATOM *at, int i )
{
    int j, neigh;
    if (at[i].el_number != EL_NUMBER_C ||
         at[i].charge || at[i].radical ||
         at[i].valence != 3 || at[i].chem_bonds_valence != 4)
        return 0;
    for (j = 0; j < at[i].valence; j++)
    {
        neigh = at[i].neighbor[j];
        if (( at[neigh].el_number == EL_NUMBER_O ) &&
             !at[neigh].num_H && 1 == at[neigh].valence &&
             2 == at[neigh].chem_bonds_valence)
        {
            return j + 1; /* =O ord */
        }
    }

    return 0;
}


/****************************************************************************/
int is_C_unsat_not_arom( inp_ATOM *at, int i )
{
    int j, neigh, num_arom, num_DB;
    if (at[i].el_number != EL_NUMBER_C ||
         at[i].valence == at[i].chem_bonds_valence || /* no double/triple bonds */
         at[i].valence + 1 < at[i].chem_bonds_valence || /* >1 double bond or >=1 triple bond */
         at[i].chem_bonds_valence + at[i].num_H != 4 || /* C has wrong valence */
         at[i].charge || at[i].radical)
        return 0;
    num_arom = num_DB = 0;
    for (j = 0; j < at[i].valence; j++)
    {
        neigh = at[i].neighbor[j];
        num_arom += at[i].bond_type[j] == BOND_TYPE_ALTERN;
        if (( at[neigh].el_number == EL_NUMBER_O ||
              at[neigh].el_number == EL_NUMBER_S ) &&
             !at[neigh].num_H && 1 == at[neigh].valence &&
             2 == at[neigh].chem_bonds_valence)
        {
            continue; /* do not count double bonds to terminal =O or =S */
        }
        num_DB += at[i].bond_type[j] == BOND_TYPE_DOUBLE;
    }

    return num_DB && !num_arom;
}


/****************************************************************************/
int is_Aryl( inp_ATOM *at, int outside_point, int attachment_pont )
{
    int i, num_arom_bonds, neigh;
    if (at[attachment_pont].el_number == EL_NUMBER_C &&
         at[attachment_pont].valence == 3 && at[attachment_pont].chem_bonds_valence == 4 &&
         !at[attachment_pont].num_H && !at[attachment_pont].charge && !at[attachment_pont].radical)
    {
        for (i = 0, num_arom_bonds = 0; i < at[attachment_pont].valence; i++)
        {
            neigh = at[attachment_pont].neighbor[i];
            if (neigh != outside_point)
            {
                num_arom_bonds += ( at[attachment_pont].bond_type[i] == BOND_ALTERN &&
                    ( at[neigh].el_number == EL_NUMBER_C || at[neigh].el_number == EL_NUMBER_N ) );
            }
        }
        return ( num_arom_bonds == 2 );
    }

    return 0;
}

/****************************************************************************/
int is_Saturated_C( inp_ATOM *at, int attachment_pont )
{
    return ( at[attachment_pont].el_number == EL_NUMBER_C &&
             at[attachment_pont].valence == at[attachment_pont].chem_bonds_valence );
}


/****************************************************************************/
int is_C_Alk( inp_ATOM *at, int i, char cFlags )
{
    if (at[i].el_number == EL_NUMBER_C &&
         at[i].valence == at[i].chem_bonds_valence)
    {
        int j, k;
        U_CHAR el;
        for (j = 0; j < at[i].valence; j++)
        {
            k = at[i].neighbor[j];
            if (at[k].cFlags & cFlags)
                continue;
            el = at[k].el_number;
            if (el != EL_NUMBER_C &&
                 el != EL_NUMBER_F &&
                 el != EL_NUMBER_CL &&
                 el != EL_NUMBER_BR &&
                 el != EL_NUMBER_I)
            {
                return 0;
            }
        }
        return 1;
    }

    return 0;
}


#if( defined(DERIV_X_OXIME) || defined(DERIV_RO_COX) || defined(DERIV_DANSYL) )


/****************************************************************************
CH -- CH
/        \
(outside point)-(attachment point, C)          CH
\        /
CH -- CH
note: bond types in the ring are not checked,
we check num_H, valence, charge, radical, RingSystem,
nNumAtInRingSystem, bCutVertex
****************************************************************************/
int is_Phenyl( inp_ATOM *at, int outside_point, int attachment_point )
{
    int iNext, iCur, iNewNext, k;

    if (at[attachment_point].el_number == EL_NUMBER_C &&
         at[attachment_point].valence == 3 && /*at[attachment_point].chem_bonds_valence == 4 &&*/
         !at[attachment_point].num_H && !at[attachment_point].charge && !at[attachment_point].radical &&
         at[attachment_point].nRingSystem != at[outside_point].nRingSystem &&
         at[attachment_point].bCutVertex && at[attachment_point].nNumAtInRingSystem == 6)
    {

        for (iNext = 0; iNext < at[attachment_point].valence; iNext++)
        {
            if (at[attachment_point].neighbor[iNext] != outside_point)
            {
                break;
            }
        }
        if (iNext == at[attachment_point].valence)
        {
            return 0; /* program error*/
        }
        iCur = attachment_point;
        iNext = at[attachment_point].neighbor[iNext];
        for (k = 0; k < 5; k++)
        {
            /* here we do not check bond type in the aromatic ring */
            if (at[iNext].el_number != EL_NUMBER_C || at[iNext].valence != 2 || at[iNext].num_H != 1 || at[iNext].charge || at[iNext].radical)
            {
                return 0;
            }
            iNewNext = at[iNext].neighbor[at[iNext].neighbor[0] == iCur];
            iCur = iNext;
            iNext = iNewNext;
        }
        return ( iNext == attachment_point );
    }

    return 0;
}


/****************************************************************************
F     F
|     |
C --- C
/        \
(outside point)-(attachment point, C)          C---F
\        /
C --- C
|     |
F     F
note: bond types in the ring are not checked,
we check num_H, valence, charge, radical, RingSystem,
nNumAtInRingSystem, bCutVertex
****************************************************************************/
int is_PentaFluoroPhenyl( inp_ATOM *at,
                          int outside_point,
                          int attachment_point )
{
    int iNext, iCur, iNewNext, nF, k, i, neigh;

    if (at[attachment_point].el_number == EL_NUMBER_C &&
         at[attachment_point].valence == 3 && /*at[attachment_point].chem_bonds_valence == 4 &&*/
         !at[attachment_point].num_H && !at[attachment_point].charge && !at[attachment_point].radical &&
         at[attachment_point].nRingSystem != at[outside_point].nRingSystem &&
         at[attachment_point].bCutVertex && at[attachment_point].nNumAtInRingSystem == 6)
    {

        for (iNext = 0; iNext < at[attachment_point].valence; iNext++)
        {
            if (at[attachment_point].neighbor[iNext] != outside_point)
                break;
        }
        if (iNext == at[attachment_point].valence)
            return 0; /* program error*/
        iCur = attachment_point;
        iNext = at[attachment_point].neighbor[iNext];
        for (k = 0; k < 5; k++)
        {
            /* here we do not check bond type in the aromatic ring */
            if (at[iNext].el_number != EL_NUMBER_C || at[iNext].valence != 3 || at[iNext].num_H != 0 || at[iNext].charge || at[iNext].radical)
                return 0;
            for (i = 0, nF = 0, iNewNext = -1; i < at[iNext].valence; i++)
            {
                neigh = at[iNext].neighbor[i];
                if (neigh == iCur)
                {
                    ;
                }
                else
                    if (at[neigh].el_number == EL_NUMBER_F && at[neigh].chem_bonds_valence == 1 && !at[neigh].charge && !at[neigh].radical && !at[neigh].num_H)
                    {
                        nF++; /* terminal flourine */
                    }
                    else
                        if (iNewNext == -1)
                        {
                            iNewNext = neigh; /* Carbon will be checked on the next pass */
                        }
                        else
                        {
                            return 0;
                        }
            }
            if (iNewNext == -1 || nF != 1)
            {
                return 0;
            }
            iCur = iNext;
            iNext = iNewNext;
        }
        return ( iNext == attachment_point );
    }

    return 0;
}


/****************************************************************************/
int is_Methyl( inp_ATOM *at, int attachment_point )
{
    if (at[attachment_point].valence == 1 && at[attachment_point].chem_bonds_valence == 1 &&
         at[attachment_point].el_number == EL_NUMBER_C && at[attachment_point].num_H == 3 &&
         !at[attachment_point].charge && !at[attachment_point].radical)
    {
        /* methyl */
        return 1;
    }

    return 0;
}


/****************************************************************************/
int is_Ethyl( inp_ATOM *at, int outside_point, int attachment_point )
{
    if (at[attachment_point].valence == 2 && at[attachment_point].chem_bonds_valence == 2 &&
         at[attachment_point].el_number == EL_NUMBER_C && at[attachment_point].num_H == 2 &&
         !at[attachment_point].charge && !at[attachment_point].radical)
    {
        /* methanediyl */
        int iat_methyl = at[attachment_point].neighbor[( at[attachment_point].neighbor[0] == outside_point )];
        return is_Methyl( at, iat_methyl );
    }

    return 0;
}

/*****************************************************************************/
int is_Methyl_or_Etyl( inp_ATOM *at, int outside_point, int attachment_point )
{
    if (is_Methyl( at, attachment_point ))
    {
        return 1;
    }
    if (is_Ethyl( at, outside_point, attachment_point ))
    {
        return 2;
    }

    return 0;
}

#endif /* ( defined(DERIV_X_OXIME) || defined(DERIV_RO_COX) || defined(DERIV_DANSYL) ) */


/****************************************************************************/
int is_Si_IV( inp_ATOM *at, int i )
{
    if (at[i].el_number != EL_NUMBER_SI ||
         at[i].charge || at[i].radical || at[i].valence != 4 || at[i].chem_bonds_valence != 4)
    {
        return 0;
    }

    return 1;
}


/****************************************************************************/
int is_P_TB_N( inp_ATOM *at, int i )
{
    int j, k;
    if (at[i].el_number != EL_NUMBER_P || at[i].chem_bonds_valence - at[i].valence != 2)
        return 0;
    for (j = 0; j < at[i].valence; j++)
    {
        k = at[i].neighbor[j];
        if (at[k].el_number == EL_NUMBER_N && at[k].valence == 1 && at[k].chem_bonds_valence == 3)
        {
            return 1;
        }
    }

    return 0;
}


/****************************************************************************
X
||
[iat](iord)--C--(iord_opposite)[iat_oppposite]
****************************************************************************/
int get_CO_opposite( inp_ATOM *at,
                     int iat,
                     int iord,
                     int *iat_opposite,
                     int *iord_opposite )
{
    int i, iOpp, iC = at[iat].neighbor[iord];
    if (at[iat].bond_type[iord] == BOND_SINGLE && 3 == at[iC].valence && 4 == at[iC].chem_bonds_valence)
    {
        /* scan neighbors of iC */
        for (i = 0; i < at[iC].valence; i++)
        {
            if (iat != at[iC].neighbor[i] && BOND_SINGLE == at[iC].bond_type[i])
            {
                iOpp = *iat_opposite = at[iC].neighbor[i];
                goto scan_opposite_atom;
            }
        }
        return 0; /* failed */
    scan_opposite_atom:
        for (i = 0; i < at[iOpp].valence; i++)
        {
            if (iC == at[iOpp].neighbor[i])
            {
                *iord_opposite = i;
                return 1; /* success */
            }
        }
    }

    return 0; /* failed */
}

#if( defined(DERIV_RING_DMOX_DEOX_N) && defined(DERIV_RING_DMOX_DEOX_O) )

#define OX_RING_SIZE 5


/****************************************************************************/
int is_DERIV_RING_DMOX_DEOX_O( inp_ATOM *at,
                               int cur_atom,
                               int from_ord,
                               DERIV_AT *da,
                               DERIV_AT *da1 )
{
    /*-----------------------
    <-
    #4 #3
    R--C==N   Me or Et
    |   \ /
    |    C #2
    |   / \
    at[k]:O--CH2 Me or ET
    #0  #1
    ->
    --------------------------*/
    /*            #0           #1           #2           #3           #4 */
    static const U_CHAR bond_type[OX_RING_SIZE] = { BOND_SINGLE, BOND_SINGLE, BOND_SINGLE, BOND_DOUBLE, BOND_SINGLE };
    static const S_CHAR valence[OX_RING_SIZE] = { 2,           2,           4,           2,           3 };
    static const S_CHAR bonds_valence[OX_RING_SIZE] = { 2,           2,           4,           3,           4 };
    static const S_CHAR num_H[OX_RING_SIZE] = { 0,           2,           0,           0,           0 };

    AT_NUMB from, curr, next, nRingSystem, at_no[OX_RING_SIZE];
    S_CHAR  bond_no[OX_RING_SIZE];
    int     i, n0, n1, attach1, attach2, neigh;

    if (at[cur_atom].el_number == EL_NUMBER_O && at[cur_atom].nNumAtInRingSystem == OX_RING_SIZE)
    {
        AT_NUMB attype[OX_RING_SIZE] = { (AT_NUMB) EL_NUMBER_O, (AT_NUMB) EL_NUMBER_C, (AT_NUMB) EL_NUMBER_C, (AT_NUMB) EL_NUMBER_N, (AT_NUMB) EL_NUMBER_C };

        curr = cur_atom;
        from = at[curr].neighbor[from_ord];
        nRingSystem = at[curr].nRingSystem;
        n0 = 0;

        do
        {
            /* find next atom in a simple ring */
            for (i = 0; i < at[curr].valence &&
                ( from == ( next = at[curr].neighbor[i] ) ||
                  nRingSystem != at[next].nRingSystem ); i++)
            {
                ;
            }
            if (i == at[curr].valence)
            {
                goto check_next_derivative2;
            }
            /* check curr atom */
            if (at[curr].charge || at[curr].radical)
            {
                goto check_next_derivative2;
            }
            if (at[curr].bond_type[i] != bond_type[n0] ||
                 at[curr].valence != valence[n0] ||
                 at[curr].chem_bonds_valence != bonds_valence[n0] ||
                 at[curr].num_H != num_H[n0] ||
                 at[curr].el_number != (U_CHAR) attype[n0])
            {
                goto check_next_derivative2;
            }
            /* save current atom */
            at_no[n0] = curr;
            bond_no[n0] = i;
            /* prepare for the next */
            from = curr;
            curr = next;
            n0++;
        } while (n0 < OX_RING_SIZE && curr != cur_atom);
        /* check completion */
        if (OX_RING_SIZE != n0 || curr != cur_atom)
        {
            goto check_next_derivative2;
        }
        /* check if R is C */
        n1 = at_no[4];
        for (i = 0; i < at[n1].valence; i++)
        {
            neigh = at[n1].neighbor[i];
            if (neigh != at_no[0] && neigh != at_no[3])
            {
                if (at[neigh].el_number != EL_NUMBER_C)
                {
                    goto check_next_derivative2;
                }
                else
                {
                    break; /* checked */
                }
            }
        }
        /* check >C< attachments */
        n1 = at_no[2];
        attach1 = attach2 = 0;
        for (i = 0; i < at[n1].valence; i++)
        {
            if (at[n1].neighbor[i] != at_no[1] &&
                 at[n1].neighbor[i] != at_no[3])
            {
                if (!attach1)
                {
                    attach1 = is_Methyl_or_Etyl( at, n1, at[n1].neighbor[i] );
                }
                else
                {
                    if (!attach2)
                    {
                        attach2 = is_Methyl_or_Etyl( at, n1, at[n1].neighbor[i] );
                    }
                    else
                    {
                        goto check_next_derivative2;
                    }
                }
            }
        }
        if (!attach2 || attach2 != attach1)
        {
            goto check_next_derivative2;
        }
        /* all checks are done */
        if ( /*da &&*/ da1)
        {
            short ord_O = bond_no[0];
            /*short ord_N = !bond_no[3];*/
            AT_NUMB iN = at_no[3];
            /*AT_NUMB iO  = at_no[0];*/
            char num_2remove = 2 + attach1 + attach2;
            da1->typ[0]    /* = da1->typ[1] */ = DERIV_RING_DMOX_DEOX_O;
            da1->ord[0]    /* = da1->ord[1] */ = ord_O;
            da1->num[0]    /* = da1->num[1] */ = num_2remove;
            da1->other_atom = iN + 1;
            /*
            if ( da1->typ[0] ) {
            if ( da1->typ[0]     != DERIV_RING_DMOX_DEOX_O ||
            da1->ord[0]     != ord_O ||
            da1->num[0]     != num_2remove ||
            da1->other_atom != iN ) {
            goto check_next_derivative2;
            }
            } else {
            da1->typ[0]     = DERIV_RING_DMOX_DEOX_O;
            da1->ord[0]     = ord_O;
            da1->num[0]     = num_2remove;
            da1->other_atom = iN;
            }
            */
        }

        return DERIV_RING_DMOX_DEOX_O;
    }

check_next_derivative2:;

    return 0;
}


/****************************************************************************/
int is_DERIV_RING_DMOX_DEOX_N( inp_ATOM *at,
                               int cur_atom,
                               int from_ord,
                               DERIV_AT *da,
                               DERIV_AT *da1 )
{
    /*-----------------------
    ->
    #4 #0
    R--C==N   Me or Et
    |   \ /
    |    C #1
    |   / \
    at[k]:O--CH2 Me or ET
    #3  #2
    <-
    --------------------------*/
    /*            #0           #1           #2           #3           #4 */
    static const U_CHAR bond_type[OX_RING_SIZE] = { BOND_SINGLE, BOND_SINGLE, BOND_SINGLE, BOND_SINGLE, BOND_DOUBLE };
    static const S_CHAR valence[OX_RING_SIZE] = { 2,           4,           2,           2,           3 };
    static const S_CHAR bonds_valence[OX_RING_SIZE] = { 3,           4,           2,           2,           4 };
    static const S_CHAR num_H[OX_RING_SIZE] = { 0,           0,           2,           0,           0 };

    AT_NUMB from, curr, next, nRingSystem, at_no[OX_RING_SIZE];
    S_CHAR  bond_no[OX_RING_SIZE];
    int     i, n0, n1, attach1, attach2, neigh;

    if (at[cur_atom].el_number == EL_NUMBER_N && at[cur_atom].nNumAtInRingSystem == OX_RING_SIZE &&
         at[cur_atom].valence == 2 && at[cur_atom].chem_bonds_valence == 3)
    {
        AT_NUMB attype[OX_RING_SIZE] = { (AT_NUMB) EL_NUMBER_N, (AT_NUMB) EL_NUMBER_C, (AT_NUMB) EL_NUMBER_C, (AT_NUMB) EL_NUMBER_O, (AT_NUMB) EL_NUMBER_C };

        curr = cur_atom;
        from = at[curr].neighbor[from_ord];
        nRingSystem = at[curr].nRingSystem;
        n0 = 0;

        do
        {
            /* find next atom in a simple ring */
            for (i = 0; i < at[curr].valence &&
                ( from == ( next = at[curr].neighbor[i] ) ||
                  nRingSystem != at[next].nRingSystem ); i++)
            {
                ;
            }
            if (i == at[curr].valence)
            {
                goto check_next_derivative2;
            }
            /* check curr atom */
            if (at[curr].charge || at[curr].radical)
            {
                goto check_next_derivative2;
            }
            if (at[curr].bond_type[i] != bond_type[n0] ||
                 at[curr].valence != valence[n0] ||
                 at[curr].chem_bonds_valence != bonds_valence[n0] ||
                 at[curr].num_H != num_H[n0] ||
                 at[curr].el_number != (U_CHAR) attype[n0])
            {
                goto check_next_derivative2;
            }
            /* save current atom */
            at_no[n0] = curr;
            bond_no[n0] = i;
            /* prepare for the next */
            from = curr;
            curr = next;
            n0++;
        } while (n0 < OX_RING_SIZE && curr != cur_atom);

        /* check completion */
        if (OX_RING_SIZE != n0 || curr != cur_atom)
        {
            goto check_next_derivative2;
        }

        /* check if R is C */
        n1 = at_no[4];
        for (i = 0; i < at[n1].valence; i++)
        {
            neigh = at[n1].neighbor[i];
            if (neigh != at_no[0] && neigh != at_no[3])
            {
                if (at[neigh].el_number != EL_NUMBER_C)
                {
                    goto check_next_derivative2;
                }
                else
                {
                    break; /* checked */
                }
            }
        }

        /* check >C< attachments */
        n1 = at_no[1];
        attach1 = attach2 = 0;
        for (i = 0; i < at[n1].valence; i++)
        {
            if (at[n1].neighbor[i] != at_no[0] &&
                 at[n1].neighbor[i] != at_no[2])
            {
                if (!attach1)
                {
                    attach1 = is_Methyl_or_Etyl( at, n1, at[n1].neighbor[i] );
                }
                else
                {
                    if (!attach2)
                    {
                        attach2 = is_Methyl_or_Etyl( at, n1, at[n1].neighbor[i] );
                    }
                    else
                    {
                        goto check_next_derivative2;
                    }
                }
            }
        }
        if (!attach2 || attach2 != attach1)
        {
            goto check_next_derivative2;
        }

        /* all checks are done */
        if ( /*da &&*/ da1)
        {
            /*short ord_O = !bond_no[3];*/
            short ord_N = bond_no[0];
            /*AT_NUMB iN  = at_no[0];*/
            AT_NUMB iO = at_no[3];
            char num_2remove = 2 + attach1 + attach2;
            da1->typ[0]    /* = da1->typ[1] */ = DERIV_RING_DMOX_DEOX_N;
            da1->ord[0]    /* = da1->ord[1] */ = ord_N;
            da1->num[0]    /* = da1->num[1] */ = num_2remove;
            da1->other_atom = iO + 1;
            /*
            if ( da1->typ[0] ) {
            if ( da1->typ[0]     != DERIV_RING_DMOX_DEOX_O ||
            da1->ord[0]     != ord_O ||
            da1->num[0]     != num_2remove ||
            da1->other_atom != iN ) {
            goto check_next_derivative2;
            }
            } else {
            da1->typ[0]     = DERIV_RING_DMOX_DEOX_O;
            da1->ord[0]     = ord_O;
            da1->num[0]     = num_2remove;
            da1->other_atom = iN;
            }
            */
        }

        return DERIV_RING_DMOX_DEOX_N;
    }

check_next_derivative2:;
    return 0;
}


#endif  /* defined(DERIV_RING_DMOX_DEOX_N) && defined(DERIV_RING_DMOX_DEOX_O) */

#if( defined(DERIV_RING2_PRRLDD_OUTSIDE_PRECUR) || defined(DERIV_RING2_PPRDN_OUTSIDE_PRECUR) )

#define PRRLDD_RING_SIZE  5
#define PPRDN_RING_SIZE   6
#if( defined(DERIV_RING2_PRRLDD_OUTSIDE_PRECUR) && defined(DERIV_RING2_PPRDN_OUTSIDE_PRECUR) )
#define MIN_PRRLDD_PPRDN_RING_SIZE PRRLDD_RING_SIZE
#define MAX_PRRLDD_PPRDN_RING_SIZE PPRDN_RING_SIZE
#elif ( defined(DERIV_RING2_PRRLDD_OUTSIDE_PRECUR) )
#define MIN_PRRLDD_PPRDN_RING_SIZE PRRLDD_RING_SIZE
#define MAX_PRRLDD_PPRDN_RING_SIZE PRRLDD_RING_SIZE
#elif ( defined(DERIV_RING2_PPRDN_OUTSIDE_PRECUR) )
#define MIN_PRRLDD_PPRDN_RING_SIZE PPRDN_RING_SIZE
#define MAX_PRRLDD_PPRDN_RING_SIZE PPRDN_RING_SIZE
#else
#define MIN_PRRLDD_PPRDN_RING_SIZE (-1)
#define MAX_PRRLDD_PPRDN_RING_SIZE (-1)
#endif


/****************************************************************************/
int is_DERIV_RING2_PRRLDD_PPRDN( inp_ATOM *at,
                                 int cur_atom,
                                 int from_ord,
                                 DERIV_AT *da,
                                 DERIV_AT *da1 )
{
    /*
    #1   #2                                      #1   #2
    O   CH2--CH2                                 O   CH2--CH2            O
    ||  /     |                                  ||  /     |             ||
    R--C--N #0   |     N is cur_atom             R--C--N #0  CH2 #3  =>  R--C--OH
    from  \     |     C(IV) is from atom         from  \     |
    CH2--CH2                                     CH2--CH2
    #4   #3                                      #5   #4

    Pyrrolidides: Replace -N< with -OH           Piperidines: Replace -N< with -OH
    DERIV_RING2_PRRLDD_OUTSIDE_PRECUR            DERIV_RING2_PPRDN_OUTSIDE_PRECUR
    */
    int iat_from, i, neigh, k;
    char ord[2];

    /* check cur atom */
    if (
        at[cur_atom].el_number == EL_NUMBER_N &&
        MIN_PRRLDD_PPRDN_RING_SIZE <= at[cur_atom].nNumAtInRingSystem &&
        at[cur_atom].nNumAtInRingSystem <= MAX_PRRLDD_PPRDN_RING_SIZE &&
        at[cur_atom].valence == 3 && at[cur_atom].chem_bonds_valence == 3 &&
        !at[cur_atom].charge && !at[cur_atom].radical && !at[cur_atom].num_H &&
        /* check the "from" atom (on the left from cur atom) */
        at[iat_from = at[cur_atom].neighbor[from_ord]].el_number == EL_NUMBER_C &&
        at[iat_from].nNumAtInRingSystem == 1 &&
        at[iat_from].valence == 3 &&
        at[iat_from].chem_bonds_valence == 4 &&
        !at[iat_from].charge && !at[iat_from].radical && !at[iat_from].num_H)
    {
        /* check neighbors of the "from" atom (on the left from cur atom) */
        for (i = 0; i < at[iat_from].valence; i++)
        {
            neigh = at[iat_from].neighbor[i];
            if (neigh == cur_atom)
            {
                ;
            }
            else
            {
                if (at[iat_from].bond_type[i] == BOND_SINGLE)
                {
                    if (at[neigh].el_number != EL_NUMBER_C)
                    {
                        goto check_next_derivative;
                    }
                }
                else
                {
                    if (at[iat_from].bond_type[i] == BOND_DOUBLE)
                    {
                        if (at[neigh].el_number != EL_NUMBER_O ||
                             at[neigh].valence != 1 || at[neigh].chem_bonds_valence != 2 ||
                             at[neigh].charge || at[neigh].radical || at[neigh].num_H)
                        {
                            goto check_next_derivative;
                        }
                    }
                    else
                    {
                        goto check_next_derivative;
                    }
                }
            }
        }

        /* check the ring */
        iat_from = cur_atom;
        neigh = at[cur_atom].neighbor[!from_ord]; /* any except "from" atom */
        i = 1;
        if (at[cur_atom].nRingSystem != at[neigh].nRingSystem)
        {
            goto check_next_derivative;
        }
        do
        {
            if (at[neigh].el_number != EL_NUMBER_C ||
                 at[neigh].valence != 2 ||
                 at[neigh].chem_bonds_valence != 2 ||
                 at[neigh].num_H != 2 ||
                 at[neigh].charge || at[neigh].radical)
            {
                goto check_next_derivative;
            }
            i++; /* at[neigh] satisfied the conditions */
            k = at[neigh].neighbor[at[neigh].neighbor[0] == iat_from];
            iat_from = neigh;
            neigh = k;
        } while (neigh != cur_atom && i <= MAX_PRRLDD_PPRDN_RING_SIZE);

        if (neigh == cur_atom &&
             MIN_PRRLDD_PPRDN_RING_SIZE <= i && i <= MAX_PRRLDD_PPRDN_RING_SIZE)
        {
            if (da1)
            {
#ifdef DERIV_RING2_PRRLDD_OUTSIDE_PRECUR
                if (i == PRRLDD_RING_SIZE)
                {
                    da1->typ[0] = da1->typ[1] = DERIV_RING2_PRRLDD_OUTSIDE_PRECUR;
                }
#endif
#ifdef DERIV_RING2_PPRDN_OUTSIDE_PRECUR
                if (i == PPRDN_RING_SIZE)
                {
                    da1->typ[0] = da1->typ[1] = DERIV_RING2_PPRDN_OUTSIDE_PRECUR;
                }
#endif
                ord[0] = !from_ord;
                ord[1] = ( ( !from_ord + 1 ) ^ ( from_ord + 1 ) ) - 1; /* the 3rd possible index out of (0,1,2) */
                k = ( ord[1] < ord[0] );
                da1->ord[0] = ord[k];  /* smaller */
                da1->ord[1] = ord[!k]; /* greater */
                /*da1->num[0] = */da1->num[0] = i - 1; /* djb-rwth: unresolved issue -- revision required? / da1->num[1] = i-1? */
            }
            return i;
        }
    }

check_next_derivative:

    return 0;
}
#endif /* defined(DERIV_RING2_PRRLDD_OUTSIDE_PRECUR) || defined(DERIV_RING2_PPRDN_OUTSIDE_PRECUR) */


#ifdef DERIV_DANSYL


/****************************************************************************/
int check_arom_chain( inp_ATOM *at,
                      int cur /* first*/,
                      int from,
                      int last,
                      int len )
{
    int i, num;
    num = 0;
    do
    {
        /* check this on all except at[last], which is typically different */
        if (at[cur].el_number != EL_NUMBER_C ||
             at[cur].valence != 2 ||
             at[cur].chem_bonds_valence != 3 ||
             at[cur].num_H != 1)
        {
            goto check_next_derivative;
        }
        /* bond to the next atom - check on all, cur..last, atoms */
        i = ( at[cur].neighbor[0] == from ); /* index of a bond to the next atom */
        if (at[cur].bond_type[i] != BOND_ALTERN)
        {
            goto check_next_derivative;
        }
        num++; /* checks are complete */
               /* prepare for the next atom */
        from = cur;
        cur = at[cur].neighbor[i];
    } while (cur != last && num < len);

    return ( cur == last && ++num == len );

check_next_derivative:

    return 0;
}


/****************************************************************************

O        CH--CH
||       //     \\
R--O[cur]--S[iS]--C[a]     CH                    => R--OH
||       \      /cj bond from c
O     [b]C====C[c]   CH3
/      \    /
CH    [d]C--N[iN]-CH3
\\    //dj bond from d
CH--CH

-O[cur]=  -O- or -S- or -NH-

****************************************************************************/
int is_Dansyl( inp_ATOM *at,
               int cur_atom,
               int to_ord, DERIV_AT *da,
               DERIV_AT *da1 )
{
    int i, a, b, c, d, cj = -1, dj = -1, neigh, k, iS /* S */, iN /* N */;
    if (( (( at[cur_atom].el_number == EL_NUMBER_O || at[cur_atom].el_number == EL_NUMBER_S || at[cur_atom].el_number == EL_NUMBER_N ) &&
          at[cur_atom].valence == 2 && at[cur_atom].num_H == ( at[cur_atom].el_number == EL_NUMBER_N ) &&
          at[cur_atom].nNumAtInRingSystem == 1) ||
          (at[cur_atom].el_number == EL_NUMBER_N && at[cur_atom].valence == 3 && at[cur_atom].num_H == 0) ) &&
         at[cur_atom].valence == at[cur_atom].chem_bonds_valence &&
         at[iS /* S */ = at[cur_atom].neighbor[to_ord]].el_number == EL_NUMBER_S &&
         at[iS].valence == 4 && at[iS].chem_bonds_valence == 6) /* djb-rwth: addressing LLVM warning */
    {
        /* neighbors of S; 6=1+1+2+2, 1+1+1+3 only. Therefore, we do not need to count (=O) neighbors */
        for (i = 0, a = -1; i < at[iS].valence; i++)
        {
            if (cur_atom == ( neigh = at[iS].neighbor[i] ))
            {
                ;
            }
            else
            {
                if (at[iS].bond_type[i] == BOND_DOUBLE)
                {
                    if (at[neigh].el_number != EL_NUMBER_O ||
                         at[neigh].valence != 1 ||
                         at[neigh].chem_bonds_valence != 2 ||
                         at[neigh].num_H || at[neigh].charge || at[neigh].radical)
                    {
                        goto check_next_derivative;
                    }
                }
                else
                {
                    if (at[iS].bond_type[i] == BOND_SINGLE && a == -1)
                    {
                        if (at[neigh].el_number != EL_NUMBER_C ||
                             at[neigh].nNumAtInRingSystem != 10 ||
                             at[neigh].valence != 3 ||
                             at[neigh].chem_bonds_valence != 4 ||
                             at[neigh].num_H || at[neigh].charge || at[neigh].radical ||
                             at[neigh].bond_type[at[neigh].neighbor[0] == iS] != BOND_ALTERN)
                        {
                            goto check_next_derivative;
                        }
                        a = neigh;
                    }
                    else
                    {
                        goto check_next_derivative;
                    }
                }
            }
        }

        if (a < 0)
        {
            goto check_next_derivative;
        }

        /* at[a] - aromatic, j1 - index of its bond to its aromatic neighbor */
        /* find at[b] */
        for (i = k = 0, b = -1; i < at[a].valence; i++)
        {
            neigh = at[a].neighbor[i];
            if (neigh == iS)
            {
                k += 1;
            }
            else
            {
                if (at[a].bond_type[i] != BOND_ALTERN)
                {
                    goto check_next_derivative; /* not Dansyl */
                }
                if (at[neigh].valence == 3 && b == -1)
                {
                    b = neigh;
                    k += 10;
                }
                else
                {
                    if (at[neigh].valence == 2 && at[neigh].chem_bonds_valence == 3)
                    {
                        k += 100;
                    }
                    else
                    {
                        goto check_next_derivative; /* not Dansyl */
                    }
                }
            }
        }

        /* structure check: C[b] */
        if (k != 111)
        {
            goto check_next_derivative; /* not Dansyl */
        }
        if (at[b].el_number != EL_NUMBER_C ||
             at[b].valence != 3 ||
             at[b].chem_bonds_valence != 4 ||
             at[b].charge || at[b].radical || at[b].num_H)
        {
            goto check_next_derivative; /* not Dansyl */
        }

        /* find at[c] */
        for (i = k = 0, c = -1; i < at[b].valence; i++)
        {
            if (at[b].bond_type[i] != BOND_ALTERN)
            {
                goto check_next_derivative; /* not Dansyl */
            }
            neigh = at[b].neighbor[i];
            if (neigh == a)
            {
                k += 1;
            }
            else
            {
                if (at[neigh].valence == 3 && c == -1)
                {
                    c = neigh;
                    k += 10;
                }
                else
                {
                    if (at[neigh].valence == 2 && at[neigh].chem_bonds_valence == 3)
                    {
                        k += 100;
                    }
                    else
                    {
                        goto check_next_derivative; /* not Dansyl */
                    }
                }
            }
        }

        if (k != 111)
        {
            goto check_next_derivative; /* not Dansyl */
        }

        /* structure check: C[c] */
        if (at[c].el_number != EL_NUMBER_C ||
             at[c].valence != 3 ||
             at[c].chem_bonds_valence != 4 ||
             at[c].charge || at[c].radical || at[c].num_H)
        {
            goto check_next_derivative; /* not Dansyl */
        }

        /* find at[d] */
        for (i = k = 0, d = -1; i < at[c].valence; i++)
        {
            if (at[c].bond_type[i] != BOND_ALTERN)
            {
                goto check_next_derivative; /* not Dansyl */
            }
            neigh = at[c].neighbor[i];
            if (neigh == b)
            {
                k += 1;
            }
            else
            {
                if (at[neigh].valence == 3 && d == -1)
                {
                    d = neigh;
                    k += 10;
                }
                else
                {
                    if (at[neigh].valence == 2)
                    {
                        cj = i;
                        k += 100;
                    }
                    else
                    {
                        goto check_next_derivative; /* not Dansyl */
                    }
                }
            }
        }
        if (k != 111)
        {
            goto check_next_derivative; /* not Dansyl */
        }

        /* structure check: C[d] */
        if (at[d].el_number != EL_NUMBER_C ||
             at[d].valence != 3 ||
             at[d].chem_bonds_valence != 4 ||
             at[d].charge || at[d].radical || at[d].num_H)
        {
            goto check_next_derivative; /* not Dansyl */
        }
        /* find at[iN]  */
        for (i = k = 0, iN = -1; i < at[d].valence; i++)
        {
            neigh = at[d].neighbor[i];
            if (neigh == c)
            {
                k += 1;
            }
            else
            {
                if (at[neigh].valence == 3 && iN == -1)
                {
                    iN = neigh;
                    k += 10;
                }
                else
                {
                    if (at[neigh].valence == 2 && at[neigh].chem_bonds_valence == 3 && at[d].bond_type[i] == BOND_ALTERN)
                    {
                        dj = i;
                        k += 100;
                    }
                    else
                    {
                        goto check_next_derivative; /* not Dansyl */
                    }
                }
            }
        }

        /* structure check: at[iN] */
        if (k != 111)
        {
            goto check_next_derivative; /* not Dansyl */
        }
        if (at[iN].el_number != EL_NUMBER_N ||
             at[iN].valence != 3 ||
             at[iN].chem_bonds_valence != 3 ||
             at[iN].charge || at[d].radical || at[d].num_H)
        {
            goto check_next_derivative; /* not Dansyl */
        }

        /* find attached to N 2 methyls */
        for (i = 0; i < at[iN].valence; i++)
        {
            if (( neigh = at[iN].neighbor[i] ) != d &&
                 !is_Methyl( at, neigh ))
            {
                goto check_next_derivative; /* not Dansyl */
            }
        }

        /* check aromatic chain d-dj-...b and c-cj-...a */
        if (check_arom_chain( at, at[d].neighbor[dj] /* first*/, d /*from*/, b /*to*/, 4 ) &&
             check_arom_chain( at, at[c].neighbor[cj] /* first*/, c /*from*/, a /*to*/, 4 ))
        {
            if (da1)
            {
                da1->typ[0] = DERIV_DANSYL;
                da1->ord[0] = to_ord;
                da1->num[0] = 16;
            }
            return DERIV_DANSYL;
        }
    }

check_next_derivative:

    return 0;
}
#endif /* DERIV_DANSYL */


/****************************************************************************/
int is_possibly_deriv_neigh( inp_ATOM *at,
                             int iat,
                             int iord,
                             int type,
                             char cFlags )
{
    int neigh = at[iat].neighbor[iord]; /* inside derivatizing agent */
    int neigh_from = -1;
    U_CHAR el = at[neigh].el_number;
    int    bOk = 0;
    switch (type)
    {
        case DERIV_BRIDGE_O:
            neigh_from = at[iat].neighbor[!iord]; /* inside precursor
                                                  neigh_from  iat
                                                         -> A--O--B -> traversing from A(neigh_from) to B(neigh); may we cut O--B bond? */
                                                  /* do not cut bond "---" in A=Si(IV), B(=O), B=C: Si(IV)-O---B(=O) */
            if (!( is_C_or_S_DB_O( at, neigh ) && is_Si_IV( at, neigh_from ) ) &&
                 !is_C_unsat_not_arom( at, neigh ))
            {
                bOk = ( el == EL_NUMBER_C ||
                        el == EL_NUMBER_SI ||
                        el == EL_NUMBER_S ||
                        el == EL_NUMBER_P ) &&
                    is_deriv_chain2( at, iat, DERIV_BRIDGE_O, -1, iord, 0, NULL, 0, NULL, 0, NULL );
            }
            break;

#ifdef DERIV_RO_COX
        case DERIV_RO_COX:
            /*           iord
            -> R-O--[C(=O)-B]; -B: -CH3, C[n]F[2n+1] 0 < n < 4; may we cut O--C bond? */
            neigh_from = at[iat].neighbor[!iord];
            if (at[neigh_from].el_number == EL_NUMBER_C &&
                 at[iat].el_number == EL_NUMBER_O &&
                 at[neigh].el_number == EL_NUMBER_C &&
                 is_C_or_S_DB_O( at, neigh ))
            {
                bOk = 1; /*is_deriv_chain2( at, iat, DERIV_RO_COX, iord, 0 ); does nothing */
            }
            break;
#endif

        case DERIV_BRIDGE_NH:
            /* -> A--NH--B -> traversing from A(neigh_from) to B(neigh); may we cut NH--B bond? */
            bOk = ( is_C_or_S_DB_O( at, neigh ) ||
                    /*is_C_Alk( at, neigh, cFlags ) ||*/
                    is_Si_IV( at, neigh ) /*||
                                          is_P_TB_N( at, neigh )*/ ) && !( is_C_unsat_not_arom( at, neigh ) ) &&
                is_deriv_chain2( at, iat, DERIV_BRIDGE_NH, -1, iord, 0, NULL, 0, NULL, 0, NULL );
            break;
        case DERIV_AMINE_tN:
            bOk = ( is_C_or_S_DB_O( at, neigh ) ||
                    /*is_C_Alk( at, neigh, cFlags ) ||*/
                    is_Si_IV( at, neigh ) /*||
                                          is_P_TB_N( at, neigh )*/ ) && !( is_C_unsat_not_arom( at, neigh ) ) &&
                is_deriv_chain2( at, iat, DERIV_AMINE_tN, -1, iord, 0, NULL, 0, NULL, 0, NULL );
            break;
    }

    return bOk;
}


/****************************************************************************
Determines derivative type on the forward step of the DFS
****************************************************************************/
int get_traversed_deriv_type( inp_ATOM *at,
                              DERIV_AT *da,
                              int k,
                              DERIV_AT *da1,
                              char cFlags )
{
    /* at[k] is attachment point of the precursor */
    /* at[(int)at[k].neighbor[m]] is inside precursor */
    /* at[(int)at[k].neighbor[!m]] is inside derivatizing agent */
    /* !!! Except DERIV_RING_O_OUTSIDE_PRECURSOR, DERIV_RING_NH_OUTSIDE_PRECURSOR !!! */
    /* when at[k] is B or C attached to two atoms of the precursor */
    int i, j, m, n1, nBlockSystemFrom, nOrdBack1, nOrdBack2, nOrdBack3, nBackType1, nBackType2;

#if( defined(DERIV_X_OXIME) || defined(DERIV_RO_COX) || defined(DERIV_RING_DMOX_DEOX_N) && defined(DERIV_RING_DMOX_DEOX_O) )
    int n0, n2, n3;
#endif

    memset( da1, 0, sizeof( *da1 ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    if (at[k].cFlags & cFlags)
    {
        return 0;
    }

    for (m = 0; m < at[k].valence && !( at[(int) at[k].neighbor[m]].cFlags & cFlags ); m++)
    {
        ;
    }
    if (m == at[k].valence)
    {
        return -1;  /* error: at least one neighbor must have cFlags */
                    /* traversing at[k] from at[(int)at[k].neighbor[m]] */
    }
    if (at[k].valence == 1 && at[k].num_H && (
        at[k].el_number == EL_NUMBER_O ||
        at[k].el_number == EL_NUMBER_N ||
        at[k].el_number == EL_NUMBER_S ||
        at[k].el_number == EL_NUMBER_P ))
    {
        return DERIV_NOT;
    }
    if (is_el_a_metal( at[k].el_number ))
    {
        return DERIV_NOT;
    }
#ifdef NEVER
    if (at[k].el_number == EL_NUMBER_N && at[k].valence >= 2 && at[k].chem_bonds_valence <= 3)
    {
        return DERIV_NOT; /* prohibit -N-, -N=, allow -N# as in isocyano -N#C or NO2 */
    }
#endif
    /* m is the ord of the bond from which at[k] was reached first time */
    if (da[k].typ[0] && ( da[k].typ[0] & DERIV_UNEXPADABLE ) == da[k].typ[0])
    {
        return 0;
    }
#ifdef DERIV_X_OXIME
    if (at[k].nNumAtInRingSystem == 1 && at[k].el_number == EL_NUMBER_N &&
         at[k].valence == 2 && at[k].chem_bonds_valence == 3 &&
         !at[k].num_H && !at[k].charge && !at[k].radical &&
         at[n0 = at[k].neighbor[m]].el_number == EL_NUMBER_C &&  /* inside precursor */
         at[n1 = at[k].neighbor[!m]].el_number == EL_NUMBER_O && /* inside derivatizing agent */
         at[k].bond_type[m] == BOND_DOUBLE && at[k].bond_type[!m] == BOND_SINGLE &&
         at[n0].valence == 3 - at[n0].num_H && at[n0].chem_bonds_valence == 4 - at[n0].num_H && !at[n0].charge && !at[n0].radical && /* C */
         at[n1].valence == 2 && at[n1].chem_bonds_valence == 2 && !at[n1].charge && !at[n1].radical                   /* O */
         )
    {
        /* found C==N--O */
        /* traversing from C to O; C(at[neighbor[m]]) has cFlag; N is at[k]; N-O bond is to be broken */
        /* m     !m
           n0 k  n1 n2       n2    n2  n3   n2  n3      n2  n3...
           C==N--O--R; -R: -CH3, -CH2-CH3, -Si(CH3)3, -CH2-C6H5; cut N-O and replace =N- with =O 2013-08-22 DT */
        /* check other neighbors of C: they should be C,H or C,C */
        for (i = 0; i < at[n0].valence; i++)
        {
            if (at[n0].neighbor[i] != k && at[at[n0].neighbor[i]].el_number != EL_NUMBER_C)
            {
                goto check_next_derivative; /* wrong neighbor */
            }
        }
        /* found    C
        |
        C==N==O
        |
        C
        */
        /* find other neighbor of O */
        n2 = at[n1].neighbor[at[n1].neighbor[0] == k];
        if (is_Si_IV( at, n2 ))
        {
            int n4 = -1;
            for (i = 0; i < at[n2].valence; i++)
            {
                int n3 = at[n2].neighbor[i];
                if (n3 == n1)
                {
                    continue; /* atom O */
                }
                if (at[n3].el_number != EL_NUMBER_C || at[n3].charge ||
                     at[n3].radical || at[n3].chem_bonds_valence != at[n3].valence)
                    goto check_next_derivative; /* wrong neighbor */
                if (n4 == -1 && at[n3].valence == 4 && !at[n3].num_H)
                {
                    n4 = n3; /* possibly tret-butyl */
                }
                else
                {
                    if (at[n3].chem_bonds_valence != 1 || at[n3].num_H != 3)
                    {
                        goto check_next_derivative; /* wrong neighbor */
                                                    /* methyl identified */
                    }
                }
            }
            if (n4 == -1)
            {
#ifdef UNDERIV_X_OXIME_TMS
                /* found    C
                |        n2
                C==N==O--Si(CH3)3
                |
                C
                */
                da1->ord[0] = !m;         /* ord of neighbor O, already checked */
                da1->typ[0] = DERIV_X_OXIME;   /* type */
                da1->num[0] = 5; /* 5 atoms: -O-Si(CH3)3 -O-TMS*/
                return DERIV_X_OXIME;   /* >C=N-O-Si(CH3)3 */
#else
                goto check_next_derivative; /* wrong neighbor */
#endif
            }
#ifndef UNDERIV_X_OXIME_TBDMS
            goto check_next_derivative; /* do not include TBDMS */
#endif

            for (i = 0; i < at[n4].valence; i++)
            {
                int n3 = at[n4].neighbor[i];
                if (n3 == n2)
                {
                    continue; /* atom Si */
                }
                if (at[n3].el_number != EL_NUMBER_C || at[n3].charge ||
                     at[n3].radical || at[n3].chem_bonds_valence != at[n3].valence)
                {
                    goto check_next_derivative; /* wrong neighbor */
                }
                if (at[n3].chem_bonds_valence != 1 || at[n3].num_H != 3)
                {
                    goto check_next_derivative; /* wrong neighbor */
                                                /* methyl identified */
                }
            }
            da1->ord[0] = !m;         /* ord of neighbor O, already checked */
            da1->typ[0] = DERIV_X_OXIME;   /* type */
            da1->num[0] = 8; /* 8 atoms: -O-Si(CH3)2-C(CH3)3 -O-TBDMS  */
            return DERIV_X_OXIME;   /* >C=N-O-Si(CH3)2-C(CH3)3 */

        }
        if (at[n2].el_number != EL_NUMBER_C)
        {
            goto check_next_derivative; /* wrong neighbor */
        }
        if (at[n2].chem_bonds_valence == 1 && at[n2].num_H == 3 && !at[n2].charge && !at[n2].radical)
        {
            da1->ord[0] = !m;         /* ord of neighbor O, already checked */
            da1->typ[0] = DERIV_X_OXIME;   /* type */
            da1->num[0] = 2; /* 2 atoms: -O-CH3 */
            return DERIV_X_OXIME;   /* >C=N-O-CH3 */
        }
        if (at[n2].valence != 2 || at[n2].chem_bonds_valence != 2 || at[n2].num_H != 2 || at[n2].charge || at[n2].radical)
        {
            goto check_next_derivative; /* wrong neighbor */
        }
        n3 = at[n2].neighbor[at[n2].neighbor[0] == n1];
        if (at[n3].chem_bonds_valence == 1 && at[n3].num_H == 3 && !at[n3].charge && !at[n3].radical)
        {
            da1->ord[0] = !m;         /* ord of neighbor O, already checked */
            da1->typ[0] = DERIV_X_OXIME;   /* type */
            da1->num[0] = 3; /* 3 atoms: -O-CH2-CH3 */
            return DERIV_X_OXIME;   /* >C=N-O-CH2-CH3 */
        }
        if (at[n3].valence == 3 && at[n3].chem_bonds_valence == 4 && !at[n3].num_H && !at[n3].charge && !at[n3].radical &&
             at[n3].nRingSystem != at[n2].nRingSystem && at[n3].bCutVertex && at[n3].nNumAtInRingSystem == 6)
        {
            if (is_Phenyl( at, n2, n3 ))
            {
                da1->ord[0] = !m;         /* ord of neighbor O, already checked */
                da1->typ[0] = DERIV_X_OXIME;   /* type */
                da1->num[0] = 8; /* 8 atoms: O--CH2-C6H5 */
                return DERIV_X_OXIME;   /* >C=N-O-CH2-C6H5 */
            }
        }
    }

check_next_derivative:
#endif  /* DERIV_X_OXIME */
#ifdef DERIV_DANSYL
    if (at[k].nNumAtInRingSystem == 1 &&
        ( (( at[k].el_number == EL_NUMBER_O || at[k].el_number == EL_NUMBER_S ) && at[k].valence == 2) ||
          (at[k].el_number == EL_NUMBER_N && at[k].valence == 2 && at[k].num_H == 1) || (at[k].valence == 3 && at[k].num_H == 2) )) /* djb-rwth: addressing LLVM warnings */
    {
        DERIV_AT da2;
        memset( &da2, 0, sizeof( da2 ) ); /* djb-rwth: memset_s C11/Annex K variant? */
        for (j = 0, n1 = 0; j < at[k].valence; j++)
        {
            if (j == m)
            {
                continue;
            }
            if (at[i = at[k].neighbor[j]].el_number == EL_NUMBER_S &&
                 at[i].valence == 4 && at[i].chem_bonds_valence == 6 &&
                 is_Dansyl( at, k, j, da, &da2 ))
            {
                n1++;
            }
        }
        if (n1 == 1)
        {
            *da1 = da2;
            return DERIV_DANSYL;
        }
    }
#endif

    if (at[k].nNumAtInRingSystem == 1 && ( at[k].el_number == EL_NUMBER_O || at[k].el_number == EL_NUMBER_S ) &&
         at[k].valence == 2 && at[k].chem_bonds_valence == 2 &&
         !at[k].num_H && !at[k].charge && !at[k].radical)
    {
        /*   at[k].neighbor[m] k n1==at[k].neighbor[!m]  */
        /*               -> A--O--B -> traversing from A to B; cut O--B */
        /* check for carboxy A(=O)-O-B and A--O--B(=O) */
        /*int has_A_CO   = is_C_or_S_DB_O( at, at[k].neighbor[m] );*/
        int has_B_CO = is_C_or_S_DB_O( at, n1 = at[k].neighbor[!m] );/* B is C(=o) or S(=O) */
        int is_A_Si_IV = is_Si_IV( at, at[k].neighbor[m] ); /* A is >Si< */
                                                            /* int is_B_Si_IV = is_Si_IV( at, at[k].neighbor[!m] );*/

#ifdef DERIV_RO_COX
                                                            /*               n3
                                                            at[k].neighbor[m]  k n1 n2
                                                            R--O--C--X; -X = -CH3,  -Phenyl,  -C[n]F[2n+1] 0 < n < 4
                                                            ||        (acetate)(benzoate)
                                                            O
                                                            */
        n3 = at[k].neighbor[m];  /* R */
        if (has_B_CO && is_C_DB_O( at, n1 ) &&  /* B:n1 is >C=O  */
             !is_silyl2( at, n3, k ) &&
             !is_el_a_metal( at[n3].el_number ))
        {
            for (j = 0; j < at[n1].valence; j++)
            {
                if (at[n1].neighbor[j] != k && at[n1].bond_type[j] == BOND_SINGLE)
                {
                    /* the only suspected neighbor */
                    n2 = at[n1].neighbor[j]; /* X */
                    n0 = is_CF3_or_linC3F7a(at, n2, n1); /* djb-rwth: addressing LLVM warning */
                    if (n0) 
                    {
                        n0 = 2 + 3 * n0 + 1; /* (6,9,12 atoms) -C(=O)C[n]F[2n+1]; is_CF3_or_linC3F7a returns n */
                    }
#ifdef UNDERIV_RO_COX_Me
                    else
                        if (is_Methyl( at, n2 ))
                        {
                            /* methyl */
                            if (!is_Methyl( at, n3 ) && !is_Ethyl( at, k, n3 ))
                            {
                                n0 = 3; /* 3 atoms: -C(=O)-CH3 */
                            }
                        }
#endif /* UNDERIV_RO_COX_Me */
#ifdef UNDERIV_RO_COX_Et
                        else
                            if (is_Ethyl( at, n1, n2 ))
                            { /* ethyl */
                                if (!is_Methyl( at, n3 ) && !is_Ethyl( at, k, n3 ))
                                {
                                    n0 = 4; /* 4 atoms: -C(=O)-CH2-CH3 */
                                }
                            }
#endif /* UNDERIV_RO_COX_Et */
#ifdef UNDERIV_RO_COX_BENZOATES
                            else
                                if (is_Phenyl( at, n1, n2 ))
                                {
                                    if (!is_Methyl( at, n3 ) && !is_Ethyl( at, k, n3 ))
                                    {
                                        n0 = 2 + 6; /* 8  atoms -C(=O)-C6H5 */
                                    }
                                }
#endif /* UNDERIV_RO_COX_BENZOATES */
#ifdef UNDERIV_RO_COX_PENTAFLOUROBENZOATES
                                else
                                    if (is_PentaFluoroPhenyl( at, n1, n2 ))
                                    {
                                        if (!is_Methyl( at, n3 ) && !is_Ethyl( at, k, n3 ))
                                        {
                                            n0 = 13; /* 13  atoms -C(=O)-C6F5 */
                                        }
                                    }
#endif  /* UNDERIV_RO_COX_PENTAFLOUROBENZOATES */
                    if (n0)
                    {
                        da1->ord[0] = !m;         /* ord of at[k]'s, that is, -O-'s neighbor X in C(=O)-X */
                        da1->typ[0] = DERIV_RO_COX;   /* type */
                        da1->num[0] = n0; /* num atoms in derivatizing attachement */
                        return DERIV_RO_COX;   /* R--O--C(=O)--X; -X = -CH3, -C[n]F[2n+1] 0 < n < 4, -Phenyl, -C6F5 */
                    }
                    break;
                }
            }
        }
#endif   /* DERIV_RO_COX */

        if (is_A_Si_IV && has_B_CO)
        {
            /*                                             precursor | deriv.agent */
            ; /* do not cut bond --- in A=>Si<, B(=O), B=C,S: Si(IV)-O---B(=O) */
        }
        else
        {
            /* at[k] is precursor's attachment point;
            at[at[k].neighbor[!m]] belongs to drivatizing agent,
            at[at[k].neighbor[m]]  was marked (from_atom) */
            if (is_possibly_deriv_neigh( at, k, !m, DERIV_BRIDGE_O, cFlags ))
            {
                da1->ord[0] = !m;         /* ord of neighbor B, not traversed yet */
                da1->typ[0] = DERIV_BRIDGE_O;   /* type */
                return DERIV_BRIDGE_O;   /* Representative: R-C(=O)-O(!m)--[D]  */
            }
        }
    }

#if( defined(DERIV_RING_DMOX_DEOX_N) && defined(DERIV_RING_DMOX_DEOX_O) )
    /*
    R--C==N   Me or Et
    |   \ /
    |    C
    |   / \
    O--CH2 Me or ET
    */
    if (at[k].el_number == EL_NUMBER_O && at[k].nNumAtInRingSystem == 5 &&
         is_DERIV_RING_DMOX_DEOX_O( at, k, m, da, da1 ))
    {
        return DERIV_RING_DMOX_DEOX_O;
    }
    if (at[k].el_number == EL_NUMBER_N && at[k].nNumAtInRingSystem == 5 &&
         at[k].valence == 2 && at[k].chem_bonds_valence == 3 &&
         is_DERIV_RING_DMOX_DEOX_N( at, k, m, da, da1 ))
    {
        return DERIV_RING_DMOX_DEOX_N;
    }
#endif

    /* R1--NH--R2 */
    if (at[k].bCutVertex && at[k].el_number == EL_NUMBER_N &&
         at[k].valence == 2 && at[k].chem_bonds_valence == at[k].valence &&
         at[k].valence + at[k].num_H == 3 && !at[k].charge && !at[k].radical)
    {
        da1->ord[0] = !m;         /* ord of neighbor B, not traversed yet */
        da1->typ[0] = DERIV_BRIDGE_NH;   /* type */
        return DERIV_BRIDGE_NH;   /* R1-NH-R2  amine */
    }
    /*
    R2
    /
    R1----N
    \
    R3
    */


    if (at[k].bCutVertex && at[k].el_number == EL_NUMBER_N &&
         at[k].valence == 3 && at[k].chem_bonds_valence == at[k].valence &&
         at[k].valence + at[k].num_H == 3 && !at[k].charge && !at[k].radical)
    {
#if( defined(DERIV_RING2_PRRLDD_OUTSIDE_PRECUR) || defined(DERIV_RING2_PPRDN_OUTSIDE_PRECUR) )
        if (at[j = at[k].neighbor[m]].el_number == EL_NUMBER_C &&
             at[j].valence == 3 && at[j].chem_bonds_valence == 4 &&
             ( i = is_DERIV_RING2_PRRLDD_PPRDN( at, k, m, da, da1 ) ))
        {
            return i;
        }
        else
#endif
        {
            int rm1 = ( at[at[k].neighbor[m]].nRingSystem == at[at[k].neighbor[( m + 1 ) % 3]].nRingSystem );
            int rm2 = ( at[at[k].neighbor[m]].nRingSystem == at[at[k].neighbor[( m + 2 ) % 3]].nRingSystem );
            int r12 = ( at[at[k].neighbor[( m + 1 ) % 3]].nRingSystem == at[at[k].neighbor[( m + 2 ) % 3]].nRingSystem );
            int ord[2] = { -1, -1 };
            i = 0; /* type: tertriary amine: DERIV_AMINE_tN */
            switch (rm1 + rm2 + r12)
            {
                case 0:
                    /* -N< has no ring bonds */
                    if (is_possibly_deriv_neigh( at, k, ( m + 1 ) % 3, DERIV_AMINE_tN, cFlags ))
                    {
                        ord[i++] = ( m + 1 ) % 3; /* ord of a non-ring neighbor, not traversed yet */
                    }
                    if (is_possibly_deriv_neigh( at, k, ( m + 2 ) % 3, DERIV_AMINE_tN, cFlags ))
                    {
                        ord[i++] = ( m + 2 ) % 3; /* ord of another non-ring neighbor, not traversed yet */
                    }
                    if (i == 2 && ord[0] > ord[1])
                    {
                        int tmp = ord[0];
                        ord[0] = ord[1];
                        ord[1] = tmp;
                    }
                    break;

                case 1:
                    /* -N< has one non-ring bond; do not consider [m] because it is "from" bond */
                    if (rm1 && is_possibly_deriv_neigh( at, k, ( m + 2 ) % 3, DERIV_AMINE_tN, cFlags ))
                    {
                        ord[i++] = ( m + 2 ) % 3;   /* ord of a single non-ring neighbor, not traversed yet */
                    }
                    else
                        if (rm2 && is_possibly_deriv_neigh( at, k, ( m + 1 ) % 3, DERIV_AMINE_tN, cFlags ))
                        {
                            ord[i++] = ( m + 1 ) % 3; /* ord of a single non-ring neighbor, not traversed yet */
                        }
                    /* r12 != 0 <=> at[k]neighbor[m] is the only non-ring bond; ignore it because it is "from" bond */
            }
            for (j = 0; j < i; j++)
            {
                da1->ord[j] = ord[j];
                da1->typ[j] = DERIV_AMINE_tN;
            }
            if (i)
            {
                return DERIV_AMINE_tN;
            }
            return 0; /* all neighbors or two untraversed edges are in one ring system */
        }
    }

    /*-----------------------------------------------------------------*/
    /* DERIV_RING_O_OUTSIDE_PRECURSOR, DERIV_RING_NH_OUTSIDE_PRECURSOR */
    if (at[k].bCutVertex && /* DD */
         at[k].valence == at[k].chem_bonds_valence &&
         ( !at[k].num_H || (at[k].el_number == EL_NUMBER_C && 1 == at[k].num_H) ) &&
         !at[k].charge && !at[k].radical &&
         ( (at[k].el_number == EL_NUMBER_C  && at[k].valence + at[k].num_H == 4) ||
           (at[k].el_number == EL_NUMBER_SI && at[k].valence == 4) ||
           (at[k].el_number == EL_NUMBER_B  && at[k].valence == 3) )) /* djb-rwth: addressing LLVM warning */
    {

        /*-->    j \        entering path: ->X--O--DD
        --X--O  v
        |    \ k /    DD = C, CH, Si, B
        |     DD
        |    /   \     O = O, S, NH  = at[j], going from DD
        --Y--O
        X, Y -- must be C
        */
        nBlockSystemFrom = 0;
        nBackType1 = nBackType2 = 0;
        nOrdBack1 = nOrdBack2 = nOrdBack3 = -1;
        j = (int) at[k].neighbor[m];
        ; /* X */
        if (( at[j].el_number == EL_NUMBER_O || at[j].el_number == EL_NUMBER_S ) && at[j].valence == 2 &&
             at[j].chem_bonds_valence == at[j].valence &&
             at[j].nNumAtInRingSystem >= 5 &&
             at[n1 = at[j].neighbor[at[j].neighbor[0] == k]].el_number == EL_NUMBER_C && /* X is C */ 
             !at[j].num_H && !at[j].charge && !at[j].radical) /* djb-rwth: ignoring LLVM warning: variable used */
        {
            nBackType1 = DERIV_RING_O_OUTSIDE_PRECURSOR;
            nBlockSystemFrom = at[j].nBlockSystem;
            nOrdBack1 = m; /* predecessor */
        }
        else
            if (at[j].el_number == EL_NUMBER_N && at[j].valence == 2 &&
                 at[j].chem_bonds_valence == at[j].valence &&
                 at[j].nNumAtInRingSystem >= 5 &&
                 at[n1 = at[j].neighbor[at[j].neighbor[0] == k]].el_number == EL_NUMBER_C && /* X is C */ 
                 1 == at[j].num_H && !at[j].charge && !at[j].radical) /* djb-rwth: ignoring LLVM warning: variable used */
            {
                nBackType1 = DERIV_RING_NH_OUTSIDE_PRECURSOR;
                nBlockSystemFrom = at[j].nBlockSystem;
                nOrdBack1 = m; /* predecessor */
            }
        if (nBlockSystemFrom)
        {
            int num1, num2, bFound;
            at[k].cFlags |= CFLAG_MARK_BLOCK;
            /* mark precursor atoms + at[k] */
            num1 = mark_atoms_cFlags( at, at[k].neighbor[nOrdBack1], 1, CFLAG_MARK_BLOCK );
            for (i = 0; i < at[k].valence; i++)
            {
                if (i == nOrdBack1)
                    continue;
                j = (int) at[k].neighbor[i];
                bFound = 0;
                if (at[j].cFlags & CFLAG_MARK_BLOCK)
                {
                    if (( at[j].el_number == EL_NUMBER_O || at[j].el_number == EL_NUMBER_S ) && at[j].valence == 2 &&
                         at[j].chem_bonds_valence == at[j].valence &&
                         at[j].nNumAtInRingSystem >= 5 &&
                         at[n1 = at[j].neighbor[at[j].neighbor[0] == k]].el_number == EL_NUMBER_C && /* Y is C */ 
                         !at[j].num_H && !at[j].charge && !at[j].radical) /* djb-rwth: ignoring LLVM warning: variable used */
                    {
                        bFound = 1;
                        if (nOrdBack2 < 0)
                        {
                            nOrdBack2 = i; /* predecessor #2 */
                            nBackType2 = DERIV_RING_O_OUTSIDE_PRECURSOR;
                        }
                        else
                        {
                            nOrdBack3 = i; /* predecessor #3 -- should not happen */
                        }
                    }
                    if (at[j].el_number == EL_NUMBER_N && at[j].valence == 2 &&
                         at[j].chem_bonds_valence == at[j].valence &&
                         at[j].nNumAtInRingSystem >= 5 &&
                         at[n1 = at[j].neighbor[at[j].neighbor[0] == k]].el_number == EL_NUMBER_C && /* Y is C */
                         1 == at[j].num_H && !at[j].charge && !at[j].radical) /* djb-rwth: ignoring LLVM warning: variable used */
                    {
                        bFound = 1;
                        if (nOrdBack2 < 0)
                        {
                            nOrdBack2 = i; /* predecessor #2 */
                            nBackType2 = DERIV_RING_NH_OUTSIDE_PRECURSOR;
                        }
                        else
                        {
                            nOrdBack3 = i; /* predecessor #3 -- should not happen */
                        }
                    }
                    if (!bFound)
                    {
                        nOrdBack3 = 99; /* reject: wrong neighboring atom in the same block */
                        break;
                    }
                }
            }
            num2 = unmark_atoms_cFlags( at, k, 0, CFLAG_MARK_BLOCK, CFLAG_MARK_BLOCK_INV );
            if (num1 != num2)
            {
                return -1; /* mark_atoms_cFlags() program error */
            }
            if (nOrdBack2 >= 0 && nOrdBack3 < 0)
            {
                /* note: da1 refers to DD, which is a neighbor of 2 precursor atoms; ord point to precursor attachment points, O */
                if (nOrdBack1 < nOrdBack2)
                {
                    da1->ord[0] = nOrdBack1;  /* ord of a ring neighbor, traversed */
                    da1->typ[0] = nBackType1;
                    da1->ord[1] = nOrdBack2;  /* ord of another ring neighbor, not traversed yet */
                    da1->typ[1] = nBackType2;
                }
                else
                {
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


/****************************************************************************/
int add_to_da( DERIV_AT *da, DERIV_AT *add )
{
    /* if add has more than 1 element the elements are in ascending add.ord[] order */
    int i, j, len_da, len_add, numAddHiPri, numDaHiPri;

    for (len_da = 0, numDaHiPri = 0; len_da < DERIV_AT_LEN && da->typ[len_da]; len_da++)
    {
        numDaHiPri += ( 0 != ( da->typ[len_da] & DERIV_UNEXPADABLE ) );
    }
    for (len_add = 0, numAddHiPri = 0; len_add < DERIV_AT_LEN && da->typ[len_add]; len_add++) /* djb-rwth: addressing coverity ID #499516 -- definitely not a copy-paste error */
    {
        numAddHiPri += ( 0 != ( add->typ[len_add] & DERIV_UNEXPADABLE ) );
    }

    /* HiPri replaces non-HiPri derivatives */
    if (numAddHiPri && !numDaHiPri)
    {
        /* no harm if already len_da=0 */
        memset( da, 0, sizeof( *da ) ); /* djb-rwth: memset_s C11/Annex K variant? */
        len_da = 0;
    }
    else
    {
        /* non-HiPri derivatives cannot be added to HiPri */
        if (!numAddHiPri && numDaHiPri)
        {
            return 0;
        }
    }

    for (j = 0; j < DERIV_AT_LEN && add->typ[j]; j++)
    {
        for (i = 0; i < len_da; i++)
        {
            if (add->ord[j] == da->ord[i])
            {
                if (add->typ[j] != da->typ[i])
                {
                    return -1; /* error, should not happen */
                }
                if (add->num[j] != da->num[i])
                {
                    return -2; /* error, should not happen */
                }
#if( defined(DERIV_RING_DMOX_DEOX_N) && defined(DERIV_RING_DMOX_DEOX_O) )
                if ((( len_da>1 || j ) && ( add->other_atom || da->other_atom )) || (1 == len_da && add->other_atom != da->other_atom)) /* djb-rwth: addressing LLVM warning */
                {
                    return -3; /* other_atom implies single bond to cut */
                }
#endif
                break;
            }
        }

        if (i == len_da)
        {
            /* add->ord[j] has different ord values from all da->ord[]; add or replace */
            if (len_da < DERIV_AT_LEN)
            {
#if( defined(DERIV_RING_DMOX_DEOX_N) && defined(DERIV_RING_DMOX_DEOX_O) )
                if (( i || j ) && add->other_atom)
                {
                    return -3; /* other_atom implies single bond to cut */
                }
#endif
                da->ord[i] = add->ord[j];
                da->typ[i] = add->typ[j];
                da->num[i] = add->num[j];
#if( defined(DERIV_RING_DMOX_DEOX_N) && defined(DERIV_RING_DMOX_DEOX_O) )
                da->other_atom = add->other_atom;
#endif
                len_da++;
            }
            else
            {
                return -4; /* overflow, should not happen */
            }
        }
    }

    return 0;
}


/****************************************************************************
DFS search for atoms that do not have a flag
****************************************************************************/
int mark_atoms_deriv( inp_ATOM *at,
                      DERIV_AT *da,
                      int start,
                      int num,
                      char cFlags,
                      int *pbFound )
{
    int i, nFound = 0, ret; /* djb-rwth: removing redundant variables */
    DERIV_AT da1;
    int      ret2;   /* moved from below 2024-09-01 DT */
    DERIV_AT da2;    /* moved from below 2024-09-01 DT */
    da1.other_atom = 0; /* djb-rwth: initialisation needed for if conditons */
#if( defined(DERIV_RING_DMOX_DEOX_N) && defined(DERIV_RING_DMOX_DEOX_O) )
    /* djb-rwth: initialisation needed to avoid garbage values in add_to_da function call; fixing coverity ID #499492 */
    memset(da2.typ, 0, DERIV_AT_LEN * sizeof(da2.typ[0]));
    memset(da2.ord, '\0', DERIV_AT_LEN * sizeof(da2.ord[0]));
    memset(da2.num, '\0', DERIV_AT_LEN * sizeof(da2.num[0]));
    da2.other_atom = 0; /* djb-rwth: initialisation needed for if conditons */
#endif
    if (!( at[start].cFlags & cFlags ))
    {
        if (DERIV_NOT == ( ret = get_traversed_deriv_type( at, da, start, &da1, cFlags ) ))
        {
            nFound++; /* at[start] cannot belong to a derivatizing agent */
        }
        num++;
        at[start].cFlags |= cFlags;
        if (da1.typ[0])
        {
            /* possibly a derivatization agent attachment point. */
            /* check neighbors that have not been traversed yet */
            int n1 = 0, n2 = 0, i1 = -1, i2 = -1, nFound1 = 0, nFound2 = 0;
            switch (da1.typ[0])
            {
#if( defined(DERIV_RING_DMOX_DEOX_N) && defined(DERIV_RING_DMOX_DEOX_O) )
                case DERIV_RING_DMOX_DEOX_N:
                case DERIV_RING_DMOX_DEOX_O:
                    ret2 = get_traversed_deriv_type( at, da, da1.other_atom - 1, &da2, cFlags );
                    if (ret != ( ret2 ^ DERIV_RING_DMOX_DEOX ))
                    {
                        /* bug */
                        nFound++; /* at[start] cannot belong to a derivatizing agent */
                        goto check_neighbors; /* bypass add_to_da( da+start, &da1 ) */
                    }
                    /*at[da1.other_atom-1].cFlags |= cFlags;*/
                    n1 = da1.num[0]; /* terminal fragment has been identified; don't search subfragments */
                    nFound++;
                    break;
#endif
#ifdef DERIV_X_OXIME
                case DERIV_X_OXIME:
                    n1 = da1.num[0]; /* terminal fragment has been identified; don't search subfragments */
                    nFound++;
                    break;
#endif
#ifdef DERIV_RO_COX
                case DERIV_RO_COX:
                    n1 = da1.num[0]; /* terminal fragment has been identified; don't search subfragments */
                    nFound++;
                    break;
#endif
#if( defined(DERIV_RING2_PRRLDD_OUTSIDE_PRECUR) || defined(DERIV_RING2_PPRDN_OUTSIDE_PRECUR) || defined(DERIV_DANSYL) )
#ifdef DERIV_RING2_PRRLDD_OUTSIDE_PRECUR
                case DERIV_RING2_PRRLDD_OUTSIDE_PRECUR:
#endif
#ifdef DERIV_RING2_PPRDN_OUTSIDE_PRECUR
                case DERIV_RING2_PPRDN_OUTSIDE_PRECUR:
#endif
#ifdef DERIV_DANSYL
                case DERIV_DANSYL:
#endif
                    n1 = da1.num[0]; /* terminal fragment has been identified; don't search subfragments */
                    nFound++;
                    break;
#endif

                case DERIV_BRIDGE_O:
                case DERIV_BRIDGE_NH:
                    n1 = mark_atoms_deriv( at, da, at[start].neighbor[(int) da1.ord[0]], 0, cFlags, &nFound1 );
                    if (n1 > MAX_AT_DERIV || nFound1)
                    {
                        da1.typ[0] = 0;
                    }
                    else
                    {
                        da1.num[0] = n1;
                        nFound++;
                    }
                    break;
                case DERIV_AMINE_tN:
                    n1 = mark_atoms_deriv( at, da, at[start].neighbor[i1 = da1.ord[0]], 0, cFlags, &nFound1 ); /* djb-rwth: ignoring LLVM warning: variable used */
                    if (da1.typ[1])
                    {
                        n2 = mark_atoms_deriv( at, da, at[start].neighbor[i2 = da1.ord[1]], 0, cFlags, &nFound2 ); /* djb-rwth: ignoring LLVM warning: variable used */
                    }
                    if (0 < n1 && n1 <= MAX_AT_DERIV && !nFound1)
                    {
                        da1.num[0] = n1;
                        i = 1;
                        nFound++;
                    }
                    else
                    {
                        da1.ord[0] = da1.ord[1];
                        da1.num[0] = da1.num[1];
                        da1.typ[0] = da1.typ[1];
                        da1.typ[1] = 0;
                        i = 0;
                    }
                    if (0 < n2 && n2 <= MAX_AT_DERIV && !nFound2)
                    {
                        da1.num[i] = n2;
                        nFound++;
                    }
                    else
                    {
                        da1.typ[i] = 0;
                    }
                    break;
                case DERIV_RING_O_OUTSIDE_PRECURSOR:
                case DERIV_RING_NH_OUTSIDE_PRECURSOR:
                    for (i = 0; i < at[start].valence; i++)
                    {
                        if (i != da1.ord[0] && i != da1.ord[1] && !( at[at[start].neighbor[i]].cFlags & cFlags ))
                        {
                            if (!n1)
                            {
                                n1 = mark_atoms_deriv( at, da, at[start].neighbor[i1 = i], 0, cFlags, &nFound1 ); /* djb-rwth: ignoring LLVM warning: variable used */
                            }
                            else
                            {
                                n2 = mark_atoms_deriv( at, da, at[start].neighbor[i2 = i], 0, cFlags, &nFound2 ); /* djb-rwth: ignoring LLVM warning: variable used */
                            }
                        }
                    }
                    if (n1 + n2 + 1 > MAX_AT_DERIV || nFound1 || nFound2)
                    {
                        da1.typ[1] = da1.typ[0] = 0;
                    }
                    else
                    {
                        da1.num[0] = n1;
                        da1.num[1] = n2;
                        nFound++;
                    }
                    break;
            }

            /* --------- end of switch( da1.typ[0] ) ------- */
            if (n1 < 0)
            {
                return n1;
            }
            if (n2 < 0)
            {
                return n2; /* errors */
            }

#if( defined(DERIV_RING_DMOX_DEOX_N) && defined(DERIV_RING_DMOX_DEOX_O) )
            if (da1.other_atom)
            {
                if (da2.other_atom == start + 1)
                {
                    if ((i = add_to_da( da + da1.other_atom - 1, &da2 ))) /* djb-rwth: addressing LLVM warning */
                    {
                        return i;  /* error */
                    }
                }
                else
                {
                    return -4; /* error */
                }
            }
#endif
            if ((i = add_to_da( da + start, &da1 ))) /* djb-rwth: addressing LLVM warning */
            {
                return i;  /* error */
            }
            *pbFound += nFound1 + nFound2 + nFound;
            num += n1 + n2;
        }
        else
        {
            *pbFound += nFound;
        }
#if( defined(DERIV_RING_DMOX_DEOX_N) && defined(DERIV_RING_DMOX_DEOX_O) )
        check_neighbors:
#endif
                       for (i = 0; i < at[start].valence; i++)
                       {
                           num = mark_atoms_deriv( at, da, at[start].neighbor[i], num, cFlags, pbFound );
                           if (num < 0)
                           {
                               return num;
                           }
                       }
    }
    /* *pbFound =  number of derivatizer attachment points traversed forward from at[start] */

    return num; /* number of atoms traversed forward from at[start] */
}


/****************************************************************************/
int count_one_bond_atoms( inp_ATOM *at,
                          DERIV_AT *da,
                          int start,
                          int ord,
                          char cFlags,
                          int *bFound )
{
    int num = 0;
    if (!( at[at[start].neighbor[ord]].cFlags & cFlags ))
    {
        at[at[start].neighbor[ord]].cFlags |= cFlags;
        num++;
        num = mark_atoms_deriv( at, da, start, num, cFlags, bFound );
    }

    return num;
}


/****************************************************************************
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
- eliminated



O                 O                     O      F         O                  O
||                ||                    ||     |         ||                 ||
R--C--O->-CH3     R--C--O->-CH2--CH3    R--C--O->-C--F   R--C--O->-CF2-CF3  R--C--O->-CF2-CF2-CF3
|
F



1 atom             2 atoms            4 atoms     7 atoms          10 atoms
is_Me_or_Et()      is_Me_or_Et()         is_CF3_or_linC3F7()


A. DERIV_BRIDGE_NH, DERIV_AMINE_tN
-----------------------------------


O                 O                     O   F             O                O
||                ||                    ||  |             ||               ||
N->-C--CH3        N->-C--CH2--CH3       N->-C---C--F      N->-C--CF2-CF3   N->-C--CF2-CF2-CF3
|
F



3 atoms           5 atoms              8 atoms           7 atoms         12 atoms
is_Me_or_Et()     is_Me_or_Et()              is_CF3_or_linC3F7()

DERIV_RING_O_OUTSIDE_PRECURSOR (da contains >B- or >C< or >CH- atom)
------------

C----O               R----O                  R----O
|     \              |     \     CH3         |     \
|      >             |      >   /            |      >
|       \            |       \ /             |       \
|        B--X        |        C              |        CH--Ph
|       /            |       / \             |       /
|      >             |      >   \            |      >
|     /              |     /     CH3         |     /
C----O               R----O                  R----O

5-member             5 or 6-member           5 or 6-member
X=Me,Et, n-Butyl

2 atoms              3 atoms                 7 atoms

DERIV_RING_NH_OUTSIDE_PRECURSOR
------------

None in the list

****************************************************************************/


/****************************************************************************/
int is_silyl( inp_ATOM *at, int start, int ord_prev )
{
    int i, neigh, num_Me = 0, iC_IV = -1, iSi_IV = -1, i_C_or_Si_IV;

    if (at[start].el_number != EL_NUMBER_SI || at[start].valence != 4 ||
         at[start].valence != at[start].chem_bonds_valence ||
         at[start].charge || at[start].radical)
    {
        return 0;
    }

    for (i = 0; i < at[start].valence; i++)
    {
        if (i == ord_prev)
        {
            continue;
        }
        neigh = at[start].neighbor[i];
        if (at[neigh].charge || at[neigh].radical ||
             at[neigh].valence != at[neigh].chem_bonds_valence)
        {
            return 0;
        }
        if (at[neigh].valence == 4)
        {
            if (at[neigh].el_number == EL_NUMBER_C && iC_IV < 0 && iSi_IV < 0)
            {
                iC_IV = neigh;
            }
            else
            {
                if (at[neigh].el_number == EL_NUMBER_SI && iC_IV < 0 && iSi_IV < 0)
                {
                    iSi_IV = neigh;
                }
                else
                {
                    return 0;
                }
            }
        }
        else
        {
            if (at[neigh].valence == 1 &&
                 at[neigh].valence == at[neigh].chem_bonds_valence &&
                 at[neigh].el_number == EL_NUMBER_C && at[neigh].num_H == 3)
            {
                num_Me++;
            }
            else
            {
                return 0;
            }
        }
    }

    if (num_Me == 3 && iC_IV < 0 && iSi_IV < 0)
    {
        return 1; /* Si(CH3)3 */
    }

    /* next C(IV) or Si(IV) */
    /* this is a fix requested by Anz. and suggested by Gary 09/21/2011
    it rejects -Si(CH3)2-Si(CH3)3 and allows only -Si(CH3)2-C(CH3)3
    */
    i_C_or_Si_IV = iC_IV >= 0 ? iC_IV : -1;
    if (num_Me != 2 || i_C_or_Si_IV < 0)
    {
        return 0;
    }

    num_Me = 0;
    for (i = 0; i < at[i_C_or_Si_IV].valence; i++)
    {
        neigh = at[i_C_or_Si_IV].neighbor[i];
        if (neigh == start)
        {
            continue;
        }
        if (at[neigh].charge || at[neigh].radical ||
             at[neigh].valence != at[neigh].chem_bonds_valence)
        {
            return 0;
        }
        if (at[neigh].valence == 1 &&
             at[neigh].valence == at[neigh].chem_bonds_valence &&
             at[neigh].el_number == EL_NUMBER_C && at[neigh].num_H == 3)
        {
            num_Me++;
        }
        else
        {
            return 0;
        }
    }
    if (num_Me == 3)
    {
        return 2; /* Si(CH3)2Si/C(CH3)3 */
    }

    return 0;
}


/****************************************************************************/
int is_silyl2( inp_ATOM *at, int start, int from_at )
{
    int i, neigh, num_Me = 0, iC_IV = -1;

    if (at[start].el_number != EL_NUMBER_SI || at[start].valence != 4 ||
         at[start].valence != at[start].chem_bonds_valence ||
         at[start].charge || at[start].radical)
    {
        return 0;
    }

    for (i = 0; i < at[start].valence; i++)
    {
        neigh = at[start].neighbor[i];
        if (neigh == from_at)
        {
            continue;
        }
        if (at[neigh].charge || at[neigh].radical ||
             at[neigh].valence != at[neigh].chem_bonds_valence)
        {
            return 0;
        }
        if (at[neigh].valence == 4)
        {
            if (at[neigh].el_number == EL_NUMBER_C && iC_IV < 0)
            {
                iC_IV = neigh;
            }
            else
            {
                return 0;
            }
        }
        else
        {
            if (at[neigh].valence == 1 &&
                 at[neigh].valence == at[neigh].chem_bonds_valence &&
                 at[neigh].el_number == EL_NUMBER_C && at[neigh].num_H == 3)
            {
                num_Me++;
            }
            else
            {
                return 0;
            }
        }
    }
    if (num_Me == 3 && iC_IV < 0)
    {
        return 1; /* Si(CH3)3 */
    }

    /* next C(IV) */
    if (num_Me != 2 || iC_IV < 0)
    {
        return 0;
    }

    num_Me = 0;
    for (i = 0; i < at[iC_IV].valence; i++)
    {
        neigh = at[iC_IV].neighbor[i];
        if (neigh == start)
        {
            continue;
        }
        if (at[neigh].charge || at[neigh].radical ||
             at[neigh].valence != at[neigh].chem_bonds_valence)
        {
            return 0;
        }
        if (at[neigh].valence == 1 &&
             at[neigh].valence == at[neigh].chem_bonds_valence &&
             at[neigh].el_number == EL_NUMBER_C && at[neigh].num_H == 3)
        {
            num_Me++;
        }
        else
        {
            return 0;
        }
    }
    if (num_Me == 3)
    {
        return 2; /* Si(CH3)2-Si/C(CH3)3, not Si(CH3)2-Si(CH3)3 */
    }

    return 0;
}


/****************************************************************************/
int is_nButyl( inp_ATOM *at, int start, int ord_prev )
{
    int i, next, curr = start;
    int prev = at[curr].neighbor[ord_prev];
    const int len = 4;
    for (i = 0; i < len; i++)
    {
        if (at[curr].el_number != EL_NUMBER_C || at[curr].valence > 2 ||
             at[curr].valence != at[curr].chem_bonds_valence ||
             at[curr].chem_bonds_valence + at[curr].num_H != 4 ||
             at[curr].charge || at[curr].radical)
        {
            return 0;
        }
        if (at[curr].valence == 2)
        {
            next = at[curr].neighbor[prev == at[curr].neighbor[0]];
            prev = curr;
            curr = next;
        }
        else
        {
            return i + 1 == len;
        }
    }

    return 0;
}


/****************************************************************************/
int is_Me_or_Et( inp_ATOM *at, int start, int ord_prev )
{
    int i, neigh = -1;
    if (at[start].el_number != EL_NUMBER_C || at[start].valence > 2 ||
         at[start].valence != at[start].chem_bonds_valence ||
         at[start].chem_bonds_valence + at[start].num_H != 4 ||
         at[start].charge || at[start].radical)
        return 0;
    for (i = 0; i < at[start].valence; i++)
    {
        if (i == ord_prev)
        {
            continue;
        }
        if (neigh >= 0)
        {
            return 0;
        }

        neigh = at[start].neighbor[i];
        if (at[neigh].el_number != EL_NUMBER_C || at[neigh].valence > 1 ||
             at[neigh].valence != at[neigh].chem_bonds_valence ||
             at[neigh].chem_bonds_valence + at[neigh].num_H != 4 ||
             at[neigh].charge || at[neigh].radical)
        {
            return 0;
        }
    }

    return 1 + ( neigh >= 0 );
}


#ifdef NEVER

/****************************************************************************
CF3
-CF3  or -CF<
CF3
****************************************************************************/

/****************************************************************************/
int is_CF3_or_isoC3F7( inp_ATOM *at, int start, int ord_prev )
{
    int i, k, num_C_IV, iC_IV[2], neigh, num_F, iC;
    if (at[start].el_number != EL_NUMBER_C || at[start].valence != 4 ||
         at[start].valence != at[start].chem_bonds_valence ||
         at[start].chem_bonds_valence + at[start].num_H != 4 ||
         at[start].charge || at[start].radical)
        return 0;

    iC_IV[0] = iC_IV[1] = num_F = 0;

    for (i = num_C_IV = 0; i < at[start].valence; i++)
    {
        if (i == ord_prev)
        {
            continue;
        }

        neigh = at[start].neighbor[i];
        if (at[neigh].valence != at[neigh].chem_bonds_valence ||
             at[neigh].charge || at[neigh].radical)
            return 0;
        if (at[neigh].el_number == EL_NUMBER_F)
        {
            if (at[neigh].chem_bonds_valence + at[neigh].num_H != 1)
            {
                return 0;
            }
            num_F++;
        }
        else
        {
            if (at[neigh].el_number == EL_NUMBER_C &&
                 at[neigh].valence == 4 &&
                 !at[neigh].num_H && !at[neigh].charge && !at[neigh].radical && num_C_IV < 2)
            {

                if (num_C_IV > 1)
                    return 0;

                iC_IV[num_C_IV++] = neigh;
            }
        }
    }
    if (!num_C_IV && 3 == num_F)
    {
        return 1; /* -CF3 */
    }
    if (2 != num_C_IV || 1 != num_F)
    {
        return 0;
    }

    /* detect iso-C3F7 */
    for (k = 0; k < num_C_IV; k++)
    {
        iC = iC_IV[k];
        num_F = 0;
        for (i = 0; i < at[iC].valence; i++)
        {
            neigh = at[iC].neighbor[i];
            if (neigh == start)
            {
                continue;
            }
            if (at[neigh].valence != at[neigh].chem_bonds_valence ||
                 at[neigh].charge || at[neigh].radical)
            {
                return 0;
            }
            if (at[neigh].el_number == EL_NUMBER_F)
            {
                if (at[neigh].chem_bonds_valence + at[neigh].num_H != 1)
                {
                    return 0;
                }
                num_F++;
            }
            else
            {
                return 0;
            }
        }
        if (num_F != 3)
        {
            return 0;
        }
    }

    return 2; /* iso-C3F7 */
}

#endif


/****************************************************************************/
int is_CF3_or_linC3F7a( inp_ATOM *at, int start, int iat_prev )
{
    int i;

    for (i = 0; i < at[start].valence; i++)
    {
        if (iat_prev == at[start].neighbor[i])
        {
            return is_CF3_or_linC3F7( at, start, i );
        }
    }

    return 0;
}


/****************************************************************************/
int is_CF3_or_linC3F7( inp_ATOM *at, int start, int ord_prev )
{
    int i, num_C_IV, iC_IV, neigh, num_F, num_C = 0;
    AT_NUMB *p;

    while (num_C < 4)
    {

        if (at[start].el_number != EL_NUMBER_C || at[start].valence != 4 ||
             at[start].valence != at[start].chem_bonds_valence ||
             at[start].chem_bonds_valence + at[start].num_H != 4 ||
             at[start].charge || at[start].radical)
        {
            return 0;
        }

        iC_IV = num_F = 0;

        for (i = num_C_IV = 0; i < at[start].valence; i++)
        {
            if (i == ord_prev)
            {
                continue;
            }

            neigh = at[start].neighbor[i];
            if (at[neigh].valence != at[neigh].chem_bonds_valence ||
                 at[neigh].charge || at[neigh].radical)
            {
                return 0;
            }
            if (at[neigh].el_number == EL_NUMBER_F)
            {
                if (at[neigh].chem_bonds_valence + at[neigh].num_H != 1)
                {
                    return 0;
                }
                num_F++;
            }
            else
            {
                if (at[neigh].el_number == EL_NUMBER_C &&
                     at[neigh].valence == 4 &&
                     !at[neigh].num_H && !at[neigh].charge && !at[neigh].radical && num_C_IV < 2)
                {

                    if (num_C_IV)
                    {
                        return 0;
                    }

                    iC_IV = neigh;
                    num_C_IV++;
                }
            }
        }

        if (num_C_IV + num_F != 3)
        {
            return 0;
        }

        num_C++; /* found -CF2-C or -CF3 */
        if (!num_C_IV)
        {
            break; /* -CF3 */
        }

        /* treat next C(IV) as a new start atom */
        if ((p = is_in_the_list( at[iC_IV].neighbor, (AT_NUMB) start, at[iC_IV].valence ))) /* djb-rwth: addressing LLVM warning */
        {
            start = iC_IV;
            ord_prev = p - at[iC_IV].neighbor;
        }
        else
        {
            return -1; /* program error */
        }
    }

    /* Corrected by DT below - ? was - return num_C == 1 ? 1 : num_C == 3 ? 2 : 0;*/
    return num_C == 1 ? 1 : num_C == 2 ? 2 : num_C == 3 ? 3 : 0;
}


/****************************************************************************/
int is_phenyl( inp_ATOM *at, int start, int ord_prev )
{
    int k, neigh, cur_at, ord;
    if (at[start].el_number != EL_NUMBER_C || at[start].valence != 3 ||
         at[start].valence + 1 != at[start].chem_bonds_valence ||
         at[start].chem_bonds_valence + at[start].num_H != 4 ||
         at[start].charge || at[start].radical)
    {
        return 0;
    }

    ord = ( ord_prev + 1 ) % 3;
    cur_at = start;

    for (k = 0; k < 5; k++)
    {
        neigh = at[cur_at].neighbor[ord];
        if (at[neigh].el_number != EL_NUMBER_C || at[neigh].valence != 2 ||
             at[neigh].valence + 1 != at[neigh].chem_bonds_valence ||
             at[neigh].chem_bonds_valence + at[neigh].num_H != 4 ||
             at[neigh].charge || at[neigh].radical)
        {
            return 0;
        }
        ord = ( at[neigh].neighbor[0] == cur_at );
        cur_at = neigh;
    }

    return ( at[cur_at].neighbor[ord] == start );
}


/****************************************************************************/
int is_DERIV_RING_O_or_NH_OUTSIDE_PRECURSOR( inp_ATOM *at,
                                             int start,
                                             int num_atoms,
                                             DERIV_AT *da1,
                                             int idrv,
                                             char *szUnderiv,
                                             int lenUnderiv,
                                             char *szUnderiv2,
                                             int lenUnderiv2,
                                             BIT_UNDERIV *bitUnderiv )
{
    int i, j, k, neigh_at[2], prev_ord[2], neigh, is_B = 0, is_C = 0, n0, n1, n2, n3, n[4] = {0}, nFound, ind1, ind2; /* djb-rwth: adding variables for char -> int conversion of subscripts */
    AT_NUMB *p;
    const char *pUnk;
    char str[16] = { '\0' };
#if( defined(DERIV_RING_O_OUTSIDE_PRECURSOR) )
    char strO[8] = { '\0' };
#endif
#if( defined(DERIV_RING_NH_OUTSIDE_PRECURSOR) )
    char strN[8] = { '\0' };
#endif

    if (da1->typ[idrv] && ( da1->typ[idrv] & DERIV_RING_OUTSIDE_PRECURSOR ) == da1->typ[idrv] &&
         da1->typ[idrv + 1] && ( da1->typ[idrv + 1] & DERIV_RING_OUTSIDE_PRECURSOR ) == da1->typ[idrv + 1])
    {
        ;
    }
    else
    {
        return 0;
    }

    /*
    if ( (da1->typ[idrv] & DERIV_RING_O_OUTSIDE_PRECURSOR  || da1->typ[idrv+1] != DERIV_RING_O_OUTSIDE_PRECURSOR) &&
    (da1->typ[idrv] != DERIV_RING_NH_OUTSIDE_PRECURSOR || da1->typ[idrv+1] != DERIV_RING_NH_OUTSIDE_PRECURSOR) )
    return 0;
    */
    if (at[start].charge || at[start].radical || at[start].valence != at[start].chem_bonds_valence)
    {
        return 0; /* check bond types start-n0 and start-n3 */
    }

    /* check
    n1    n0
    R1---O
    |     \
    |      B [start]
    |     /
    R2---O
    n2    n3

    All bond are single except n1-n2 (R1-R2), which may be either single or aromatic

    */
    nFound = 0;
    ind1 = da1->ord[0] - '0'; /* djb-rwth: converting char to int for subscript use */
    ind2 = da1->ord[1] - '0'; /* djb-rwth: converting char to int for subscript use */
    n0 = at[start].neighbor[ind1];
    n3 = at[start].neighbor[ind2];
    /* search for i, j, k such that at[at[n1]neighbor[i]].neighbor[k]= at[n2]neighbor[j] */
    for (i = 0; i < at[n0].valence; i++)
    {
        if (( n1 = at[n0].neighbor[i] ) == start)
        {
            continue; /* don't go back */
        }
        if (BOND_SINGLE != at[n0].bond_type[i])
        {
            /* check bond type n0-n1 */
            continue;
        }
        for (j = 0; j < at[n1].valence; j++)
        {
            if (( n2 = at[n1].neighbor[j] ) == n0)
            {
                continue; /* don't go back */
            }
            if ((p = is_in_the_list( at[n3].neighbor, (AT_NUMB) n2, at[n3].valence ))) /* djb-rwth: addressing LLVM warning */
            {
                if (( BOND_SINGLE == at[n1].bond_type[j] || BOND_ALTERN == at[n1].bond_type[j] ) && /* check bond type n1-n2 */
                     BOND_SINGLE == at[n3].bond_type[p - at[n3].neighbor] && /* check bond type n3-n2 */
                     !nFound++)
                {
                    n[0] = n0;
                    n[1] = n1;
                    n[2] = n2;
                    n[3] = n3;
                }
            }
        }
    }

    if (nFound != 1 || at[n[1]].el_number != EL_NUMBER_C || at[n[2]].el_number != EL_NUMBER_C)
    {
        return 0;
    }

    /* n[1] and n[2] cannot have 3 neighboring heteroatoms */
    for (i = 1; i <= 2; i++)
    {
        if (at[n1 = n[i]].valence > 3)
        {
            for (k = 0, j = 0; j < at[n1].valence; j++)
            {
                k += ( at[at[n1].neighbor[j]].el_number != EL_NUMBER_C );
            }
            if (k >= 3)
            {
                return 0;
            }
        }
    }

    n0 = n[0];
    n3 = n[3];

    if (NULL != szUnderiv)
    {
#if( defined(DERIV_RING_O_OUTSIDE_PRECURSOR) )
        if (da1->typ[idrv] == DERIV_RING_O_OUTSIDE_PRECURSOR && da1->typ[idrv + 1] == DERIV_RING_O_OUTSIDE_PRECURSOR)
        {
            if (at[n0].el_number <= at[n3].el_number)
            {
                strcat(strO, at[n0].elname);
                strcat(strO, at[n3].elname);
            }
            else
            {
                strcat(strO, at[n3].elname);
                strcat(strO, at[n0].elname);
            }
        }
        else
        {
            if (da1->typ[idrv] == DERIV_RING_O_OUTSIDE_PRECURSOR)
            {
                strcat(strO, at[n0].elname);
            }
            else
            {
                if (da1->typ[idrv + 1] == DERIV_RING_O_OUTSIDE_PRECURSOR)
                {
                    strcat(strO, at[n3].elname);
                }
            }
        }
#endif

#if( defined(DERIV_RING_NH_OUTSIDE_PRECURSOR) )
        if (da1->typ[idrv] == DERIV_RING_NH_OUTSIDE_PRECURSOR && da1->typ[idrv + 1] == DERIV_RING_NH_OUTSIDE_PRECURSOR)
        {
            if (1 == at[n0].num_H && 1 == at[n3].num_H)
            {
                strcat(strN, "(NH)2");
            }
            else
            {
                if (1 == at[n0].num_H || 1 == at[n3].num_H)
                {
                    strcat(strN, "(NH)N");
                }
                else
                {
                    strcat(strN, "NN");
                }
            }
        }
        else
        {
            if (da1->typ[idrv] == DERIV_RING_NH_OUTSIDE_PRECURSOR)
            {
                strcat(strN, 1 == at[n0].num_H ? "(NH)" : "N");
            }
            else
            {
                if (da1->typ[idrv + 1] == DERIV_RING_NH_OUTSIDE_PRECURSOR)
                {
                    strcat(strN, 1 == at[n3].num_H ? "(NH)" : "N");
                }
            }
        }
#endif


        switch (da1->typ[idrv] | da1->typ[idrv + 1])
        {
#if( defined(DERIV_RING_O_OUTSIDE_PRECURSOR) && defined(DERIV_RING_NH_OUTSIDE_PRECURSOR) )
            case ( DERIV_RING_O_OUTSIDE_PRECURSOR | DERIV_RING_NH_OUTSIDE_PRECURSOR ):
                strcat(str, strN);
                strcat(str, strO); /* "(NH)O" or "(NH)S" */
                break;
#endif
#if( defined(DERIV_RING_O_OUTSIDE_PRECURSOR) )
            case ( DERIV_RING_O_OUTSIDE_PRECURSOR ):
                strcat(str, strO); /* "OO" or "OS" or "SS" */
                break;
#endif
#if( defined(DERIV_RING_NH_OUTSIDE_PRECURSOR) )
            case ( DERIV_RING_NH_OUTSIDE_PRECURSOR ):
                strcat(str, strN);
                break;
#endif
            default:
                strcat(str, "???");
                break;
        }
        strcat(str, "-");
    }
    underiv_list_add( szUnderiv, lenUnderiv, str, 0 );

    /*underiv_list_add( szUnderiv, lenUnderiv, da1->typ[idrv] == DERIV_RING_O_OUTSIDE_PRECURSOR? "OO-" : "(NH)2-", 0 );*/
    if (at[start].el_number == EL_NUMBER_B && at[start].valence == 3)
    {
        is_B = 1;
        underiv_list_add( szUnderiv, lenUnderiv, "B", 0 );
    }
    else
    {
        if (at[start].el_number == EL_NUMBER_C && ( at[start].valence == 3 || at[start].valence == 4 ) &&
             at[start].chem_bonds_valence == at[start].valence &&
             at[start].num_H + at[start].chem_bonds_valence == 4)
        {
            is_C = 1;
            underiv_list_add( szUnderiv, lenUnderiv, at[start].valence == 3 ? "CH" : "C", 0 );
        }
        else
        {
            return 0;
        }
    }

    /* locate bonds connecting >B- or >C< or >C- to the rest of the derivative */
    for (i = k = 0; i < at[start].valence; i++)
    {
        if (i != da1->ord[idrv] && i != da1->ord[idrv + 1])
        {
            neigh = at[start].neighbor[i];
            if (k < 2 && ( p = is_in_the_list( at[neigh].neighbor, (AT_NUMB) start, at[neigh].valence ) ))
            {
                neigh_at[k] = neigh;
                prev_ord[k] = p - at[neigh].neighbor;
                k++;
            }
            else
            {
                return -1; /* program error */
            }
        }
    }

    n1 = n2 = 0; /* djb-rwth: ignoring LLVM warning: variables used */
    if (is_B && k == 1 && da1->typ[idrv] == DERIV_RING_O_OUTSIDE_PRECURSOR)
    {
        if (0 < ( n1 = is_Me_or_Et( at, neigh_at[0], prev_ord[0] ) )
#ifdef UNDERIV_OOB_nButyl
             || 0 < ( n2 = is_nButyl( at, neigh_at[0], prev_ord[0] ) )
#endif
             )
        {
            underiv_list_add( szUnderiv, lenUnderiv, n1 == 1 ? "Me" : n1 == 2 ? "Et" : n2 ? "nButyl" : "???", 0 );
            /* is_B */
            underiv_list_add( szUnderiv2, lenUnderiv2, n1 == 1 ? pszDerivName[DERIV_ID_MeBorate] :
                              n1 == 2 ? pszDerivName[DERIV_ID_EtBorate] :
                              n2 ? pszDerivName[DERIV_ID_BuBorate] : "???", ' ' );
            *bitUnderiv |= n1 == 1 ? DERIV_BIT_MeBorate :
                n1 == 2 ? DERIV_BIT_EtBorate :
                n2 ? DERIV_BIT_BuBorate : DERIV_BIT_Unknown;
            return 1;
        }
    }
    else
    {
        if (is_C)
        {
            if (k == 1 && at[start].num_H == 1 && is_phenyl( at, neigh_at[0], prev_ord[0] ))
            {
                underiv_list_add( szUnderiv, lenUnderiv, "Phe", 0 );
                if (strN[0])
                {
                    if ((pUnk = underiv_list_get_last( szUnderiv, ' ' ))) /* djb-rwth: addressing LLVM warning */
                    {
                        underiv_list_add( szUnderiv2, lenUnderiv2, pUnk, ' ' );
                        *bitUnderiv |= DERIV_BIT_Unknown;
                    }
                }
                else
                {
                    underiv_list_add( szUnderiv2, lenUnderiv2, pszDerivName[DERIV_ID_Benzlidene], ' ' );
                    *bitUnderiv |= DERIV_BIT_Benzlidene;
                }
                return 1;
            }
            if (k == 2 && 0 < ( n1 = is_Me_or_Et( at, neigh_at[0], prev_ord[0] ) ) &&
                 0 < ( n2 = is_Me_or_Et( at, neigh_at[1], prev_ord[1] ) ))
            {
                if (n1 != n2)
                {
                    underiv_list_add( szUnderiv, lenUnderiv, "MeEt", 0 );
                    if ((pUnk = underiv_list_get_last( szUnderiv, ' ' ))) /* djb-rwth: addressing LLVM warning */
                    {
                        underiv_list_add( szUnderiv2, lenUnderiv2, pUnk, ' ' );
                        *bitUnderiv |= DERIV_BIT_Unknown;
                    }
                }
                else
                {
                    if (n1 == 1 || n1 == 2)
                    {
                        underiv_list_add( szUnderiv, lenUnderiv, n1 == 1 ? "Me2" : "Et2", 0 );
                        if (strN[0] || n1 != 1)
                        {
                            if ((pUnk = underiv_list_get_last( szUnderiv, ' ' ))) /* djb-rwth: addressing LLVM warning */
                            {
                                underiv_list_add( szUnderiv2, lenUnderiv2, pUnk, ' ' );
                                *bitUnderiv |= DERIV_BIT_Unknown;
                            }
                        }
                        else
                        {
                            underiv_list_add( szUnderiv2, lenUnderiv2, pszDerivName[n1 == 1 ? DERIV_ID_Acentonate : DERIV_ID_Unknown], ' ' );
                            *bitUnderiv |= ( n1 == 1 ? DERIV_BIT_Acentonate : DERIV_BIT_Unknown );
                        }
                    }
                }

                return 1;
            }
        }
    }

    /*
    if ( is_B && k == 1 && is_Me_or_Et( at, neigh_at[0], prev_ord[0]) )
    return 1;
    if ( is_B && k == 1 && is_nButyl( at, neigh_at[0], prev_ord[0]) )
    return 1;
    if ( is_C && k == 1 && at[start].num_H == 1 && is_phenyl( at, neigh_at[0], prev_ord[0]) )
    return 1;
    if ( is_C && k == 2 && is_Me_or_Et( at, neigh_at[0], prev_ord[0]) &&
    is_Me_or_Et( at, neigh_at[1], prev_ord[1]) )
    return 1;
    */

    return 0;
}


/****************************************************************************/
int is_deriv_chain2( inp_ATOM *at,
                     int start,
                     int type,
                     int num,
                     int ord,
                     int idrv,
                     char *szUnderiv,
                     int lenUnderiv,
                     char *szUnderiv2,
                     int lenUnderiv2,
                     BIT_UNDERIV *bitUnderiv )
{
    int i, k, prev_ord, neigh, iC, iO /* O or N */, iNxt, n1 = 0, n2 = 0;
    AT_NUMB *p;

    if (!type || ( type & DERIV_RING_OUTSIDE_PRECURSOR ))
    {
        return 0;
    }
    /*
    #ifdef DERIV_RING2_OUTSIDE_PRECUR
    if (type & DERIV_RING2_OUTSIDE_PRECUR)
    {
    return 1;
    }
    #endif
    */
    /* reject unexpected unsaturated */
    if (at[start].charge || at[start].radical || (at[start].valence != at[start].chem_bonds_valence
#ifdef DERIV_X_OXIME
         && type != DERIV_X_OXIME
#endif
#ifdef DERIV_RO_COX
         && type != DERIV_RO_COX
#endif
#ifdef DERIV_RING_DMOX_DEOX_N
         && type != DERIV_RING_DMOX_DEOX_N
#endif
         )) /* djb-rwth: addressing LLVM warning */
    {
        return 0;
    }

    neigh = at[start].neighbor[(int) ord];
    p = is_in_the_list( at[neigh].neighbor, (AT_NUMB) start, at[neigh].valence );
    if (!p)
    {
        return -1; /* program error */
    }

    prev_ord = p - at[neigh].neighbor;

    /* eliminate silyl possibility */
    if (type == DERIV_BRIDGE_O || type == DERIV_BRIDGE_NH || type == DERIV_AMINE_tN)
    {
        if ((n1 = is_silyl( at, neigh, prev_ord ))) /* djb-rwth: addressing LLVM warning */
        {
            if (at[start].valence != 2 || /* amine ? */
                                          /*at[start].el_number == EL_NUMBER_O && */
                 at[at[start].neighbor[!ord]].el_number != EL_NUMBER_SI)
            {
                /* Gari's request: disconnect only from C-O->-Si..., not Si-O->-Si... (e.g.  CASr.n.= 141-63-9 ). */
                /* ???? in case of type == DERIV_AMINE_tN why don't we check more neighbors ???? */
                if (NULL != szUnderiv && 0 < lenUnderiv)
                {
                    char szPrecur[16] = "R";
                    switch (type)
                    {
                        case DERIV_BRIDGE_O:
                            strcat(szPrecur, at[start].elname);
                            break;
                        case DERIV_BRIDGE_NH:
                            strcat(szPrecur, "NH");
                            break;
                        case DERIV_AMINE_tN:
                            strcat(szPrecur, "N");
                            break;
                        default:
                            strcat(szPrecur, "??");
                            break;
                    }
                    strcat(szPrecur, "-");
                    switch (n1)
                    {
                        case 1:
                            strcat(szPrecur, "TMS");
                            break;
                        case 2:
                            strcat(szPrecur, "TBDMS");
                            break;
                        default:
                            strcat(szPrecur, "???");
                            break;
                    }
                    underiv_list_add( szUnderiv, lenUnderiv, szPrecur, ' ' );
                    underiv_list_add( szUnderiv2, lenUnderiv2, pszDerivName[n1 == 1 ? DERIV_ID_TMS : n1 == 2 ? DERIV_ID_TBDMS : DERIV_ID_Unknown], ' ' );
                    *bitUnderiv |= n1 == 1 ? DERIV_BIT_TMS : n1 == 2 ? DERIV_BIT_TBDMS : DERIV_BIT_Unknown;

                    /*
                    underiv_list_add( szUnderiv, lenUnderiv, type == DERIV_BRIDGE_O? "RO-" : type == DERIV_BRIDGE_NH? "RNH-" : "RN-", ' ' );
                    underiv_list_add( szUnderiv, lenUnderiv, n1==1?"TMS" : n1==2? "TBDMS" : "???", 0 );
                    */
                }
                return 1;
            }
            else
            {
                return 0;
            }
        }
    }
#ifdef UNDERIV_SYLYL_ONLY
    return 0; /* if it is not Sylyl then it is not a derivative */
#endif

    n1 = n2 = 0; /* djb-rwth: ignoring LLVM warning: variables used */
    if (type == DERIV_BRIDGE_O)
    {
        /* check acetyl */
#if ( !defined(UNDERIV_ACETATE_Me) && !defined(UNDERIV_ACETATE_Et) && !defined(UNDERIV_ACETATE_CnF2np1) )
        return 0;
#else
        iC = at[start].neighbor[!ord];
        if (at[iC].charge || at[iC].radical || at[iC].num_H ||
             at[iC].el_number != EL_NUMBER_C || at[iC].valence != 3 ||
             at[iC].valence + 1 != at[iC].chem_bonds_valence)
            return 0;
        for (i = k = 0; i < at[iC].valence; i++)
        {
            iO = at[iC].neighbor[i];
            if (at[iO].charge || at[iO].radical || at[iO].num_H ||
                 at[iO].el_number != EL_NUMBER_O || at[iO].valence != 1 ||
                 at[iO].valence + 1 != at[iO].chem_bonds_valence)
                continue;
            k++; /* number of =O */
        }
        if (k != 1)
            return 0;
        /* check derivative */
#if defined(UNDERIV_ACETATE_Et) || defined(UNDERIV_ACETATE_Me)
        n1 = is_Me_or_Et( at, neigh, prev_ord );
        if (
#if defined(UNDERIV_ACETATE_Et) && defined(UNDERIV_ACETATE_Me)
            0 !=
#elif defined(UNDERIV_ACETATE_Et)
            2 ==
#elif defined(UNDERIV_ACETATE_Me)
            1 ==
#endif
            n1)
            ;
        else
            n1 = 0;
#endif /* defined(UNDERIV_ACETATE_Et) || defined(UNDERIV_ACETATE_Me) */
#ifdef UNDERIV_ACETATE_CnF2np1
        if (n1 <= 0)
            n2 = is_CF3_or_linC3F7( at, neigh, prev_ord );
#endif
        if (n1 || n2)
        {
            if (szUnderiv)
            {
                underiv_list_add( szUnderiv, lenUnderiv, "RCOO-", ' ' );
                underiv_list_add( szUnderiv, lenUnderiv,
                                  n1 == 1 ? "Me" :
                                  n1 == 2 ? "Et" :
                                  n2 == 1 ? "CF3" :
                                  n2 == 2 ? "C2F5" :
                                  n2 == 3 ? "C3F7" :
                                  "C?F?", 0 );
                underiv_list_add( szUnderiv2, lenUnderiv2, pszDerivName[
#if  defined(UNDERIV_ACETATE_Me)
                    n1 == 1 ? DERIV_ID_Methylation :
#endif
#if  defined(UNDERIV_ACETATE_Et)
                        n1 == 2 ? DERIV_ID_Ethylation :
#endif
                        n2 == 1 ? DERIV_ID_TFA :
                        n2 == 2 ? DERIV_ID_PFP :
                        n2 == 3 ? DERIV_ID_HFB :
                        DERIV_ID_Unknown], ' ' );
                *bitUnderiv |=
#if  defined(UNDERIV_ACETATE_Me)
                    n1 == 1 ? DERIV_BIT_Methylation :
#endif
#if  defined(UNDERIV_ACETATE_Et)
                    n1 == 2 ? DERIV_BIT_Ethylation :
#endif
                    n2 == 1 ? DERIV_BIT_TFA :
                    n2 == 2 ? DERIV_BIT_PFP :
                    n2 == 3 ? DERIV_BIT_HFB :
                    DERIV_BIT_Unknown;
            }
            return 1;
        }
        return 0;
#endif /* !defined(UNDERIV_ACETATE_Me) && !defined(UNDERIV_ACETATE_Et) && !defined(UNDERIV_ACETATE_CnF2np1) */
    }

    n1 = n2 = 0;
    if (type == DERIV_BRIDGE_NH || type == DERIV_AMINE_tN)
    {
        /* check acetyl */
        iNxt = -1;
        iC = at[start].neighbor[(int) ord];
        if (at[iC].charge || at[iC].radical || at[iC].num_H ||
             at[iC].el_number != EL_NUMBER_C || at[iC].valence != 3 ||
             at[iC].valence + 1 != at[iC].chem_bonds_valence)
        {
            return 0;
        }
        for (i = k = 0; i < at[iC].valence; i++)
        {
            iO = at[iC].neighbor[i];
            if (at[iO].charge || at[iO].radical || at[iO].num_H ||
                 at[iO].el_number != EL_NUMBER_O || at[iO].valence != 1 ||
                 at[iO].valence + 1 != at[iO].chem_bonds_valence)
            {
                if (iO != start)
                {
                    if (iNxt < 0)
                    {
                        iNxt = iO;
                    }
                    else
                    {
                        return 0;
                    }
                }
                continue;
            }
            k++; /* number of =O */
        }
        if (k != 1 || iNxt < 0)
        {
            return 0;
        }
        /* find bond from iNxt to iC */
        p = is_in_the_list( at[iNxt].neighbor, (AT_NUMB) iC, at[iNxt].valence );
        if (!p)
        {
            return -1; /* program error */
        }
        prev_ord = p - at[iNxt].neighbor;
        /* check derivative */
#if ( defined( UNDERIV_RN_AcEt ) || defined( UNDERIV_RNH_AcEt )  )
        n1 = is_Me_or_Et( at, iNxt, prev_ord );
#endif
        if (!n1)
        {
            n2 = is_CF3_or_linC3F7( at, iNxt, prev_ord );
        }
        else
        {
            if (n1 == 1)
            {
#if ( !defined(UNDERIV_RN_AcMe) && defined(DERIV_BRIDGE_tN) )
                if (type == DERIV_BRIDGE_tN) n1 = 0;
#elif ( !defined(UNDERIV_RNH_AcMe) && defined(DERIV_BRIDGE_NH) )
                if (type == DERIV_BRIDGE_NH) n1 = 0;
#endif
                ; /* keep C-compiler happy */
            }
            else
            {
                if (n1 == 2)
                {
#if ( !defined(UNDERIV_RN_AcEt) && defined(DERIV_BRIDGE_tN) )
                    if (type == DERIV_BRIDGE_tN) n1 = 0;
#elif ( !defined(UNDERIV_RNH_AcEt) && defined(DERIV_BRIDGE_NH) )
                    if (type == DERIV_BRIDGE_NH) n1 = 0;
#endif
                    ; /* keep C-compiler happy */
                }
                else
                {
                    n1 = 0;
                }
            }
        }

        if (n1 || n2)
        {
            if (szUnderiv)
            {
                underiv_list_add( szUnderiv, lenUnderiv, type == DERIV_AMINE_tN ? "RN-C(O)" : "RNH-C(O)", ' ' );
                underiv_list_add( szUnderiv, lenUnderiv,
                                  n1 == 1 ? "Me" :
                                  n1 == 2 ? "Et" :
                                  n2 == 1 ? "CF3" :
                                  n2 == 2 ? "C2F5" :
                                  n2 == 3 ? "C3F7" :
                                  "C?F?", 0 );
                /* djb-rwth: addressing coverity ID #499506 -- condition is correct for n1 != 1 */
                underiv_list_add( szUnderiv2, lenUnderiv2, pszDerivName[
#if defined(UNDERIV_RN_AcMe) || defined(UNDERIV_RNH_AcMe)
                    n1 == 1 ? DERIV_ID_Acetate :
#endif
#if defined(UNDERIV_RN_AcEt) || defined(UNDERIV_RNH_AcEt)
                        n1 == 2 ? DERIV_ID_Propanoate :
#endif
                        n2 == 1 ? DERIV_ID_TFA :
                        n2 == 2 ? DERIV_ID_PFP :
                        n2 == 3 ? DERIV_ID_HFB :
                        DERIV_ID_Unknown], ' ' );
                *bitUnderiv |=
#if defined(UNDERIV_RN_AcMe) || defined(UNDERIV_RNH_AcMe)
                    n1 == 1 ? DERIV_BIT_Acetate :
#endif
#if defined(UNDERIV_RN_AcEt) || defined(UNDERIV_RNH_AcEt)
                    n1 == 2 ? DERIV_BIT_Propanoate :
#endif
                    n2 == 1 ? DERIV_BIT_TFA :
                    n2 == 2 ? DERIV_BIT_PFP :
                    n2 == 3 ? DERIV_BIT_HFB :
                    DERIV_BIT_Unknown;
            }
            return 1;
        }

        return 0;
    }

    /*underiv_buf_clear( szUnderiv );*/
#ifdef DERIV_X_OXIME
    if (type == DERIV_X_OXIME)
    {
        if (szUnderiv)
        {
            iO = 0;
            if (num == 2)
            {
                underiv_list_add( szUnderiv, lenUnderiv, "MOX", ' ' );
                underiv_list_add( szUnderiv2, lenUnderiv2, pszDerivName[DERIV_ID_MOX], ' ' );
                *bitUnderiv |= DERIV_BIT_MOX;
            }
            else
            {
                if (num == 3)
                {
                    underiv_list_add( szUnderiv, lenUnderiv, "EtOX", ' ' );
                    underiv_list_add( szUnderiv2, lenUnderiv2, pszDerivName[DERIV_ID_EtOX], ' ' );
                    *bitUnderiv |= DERIV_BIT_EtOX;
                }
                else
                {
                    if (num == 8)
                    {
                        neigh = at[start].neighbor[ord];
                        iC = at[neigh].neighbor[start == at[neigh].neighbor[0]];
                        iO = at[iC].el_number == EL_NUMBER_SI;
                    }
                    underiv_list_add( szUnderiv, lenUnderiv, "OX-", ' ' );
                    underiv_list_add( szUnderiv, lenUnderiv, num == 5 ? "TMS" : num == 8 && iO ? "TBDMS" : num == 8 && !iO ? "CH2Phe" : "???", 0 );
                    underiv_list_add( szUnderiv2, lenUnderiv2, pszDerivName[num == 5 ? DERIV_ID_TMS : num == 8 && iO ? DERIV_ID_TBDMS : num == 8 && !iO ? DERIV_ID_BenzOX : DERIV_ID_Unknown], ' ' );
                    *bitUnderiv |= num == 5 ? DERIV_BIT_TMS : num == 8 && iO ? DERIV_BIT_TBDMS : num == 8 && !iO ? DERIV_BIT_BenzOX : DERIV_BIT_Unknown;
                }
            }
        }

        return 1;
    }
#endif

#ifdef DERIV_RO_COX
    if (type == DERIV_RO_COX)
    {
        if (szUnderiv)
        {
            underiv_list_add( szUnderiv, lenUnderiv, at[start].el_number == EL_NUMBER_O ? "RO-C(O)" :
                              at[start].el_number == EL_NUMBER_S ? "RS-C(O)" : "R?-C(O)", ' ' );
            underiv_list_add( szUnderiv, lenUnderiv, num == 3 ? "Me" :
                              num == 4 ? "Et" :
                              num == 6 ? "CF3" :
                              num == 8 ? "Phe" :
                              num == 9 ? "C2F5" :
                              num == 12 ? "C3F7" :
                              num == 13 ? "PheF5" :
                              "???", 0 );
            if (num == 13)
            {
                underiv_list_add( szUnderiv2, lenUnderiv2, underiv_list_get_last( szUnderiv, ' ' ), ' ' );
                *bitUnderiv |= DERIV_BIT_Unknown;
            }
            else
            {
                underiv_list_add( szUnderiv2, lenUnderiv2, pszDerivName[
#if  defined(UNDERIV_RO_COX_Me)
                    num == 3 ? DERIV_ID_Acetate :
#endif
#if  defined(UNDERIV_RO_COX_Et)
                        num == 4 ? DERIV_ID_Propanoate :
#endif
                        num == 6 ? DERIV_ID_TFA :
#if  defined(UNDERIV_RO_COX_BENZOATES)
                        num == 8 ? DERIV_ID_Benzoate :
#endif
                        num == 9 ? DERIV_ID_PFP :
                        num == 12 ? DERIV_ID_HFB :
                        num == 13 ? DERIV_ID_PFB :
                        DERIV_ID_Unknown], ' ' );
                *bitUnderiv |=
#if  defined(UNDERIV_RO_COX_Me)
                    num == 3 ? DERIV_BIT_Acetate :
#endif
#if  defined(UNDERIV_RO_COX_Et)
                    num == 4 ? DERIV_BIT_Propanoate :
#endif
                    num == 6 ? DERIV_BIT_TFA :
#if  defined(UNDERIV_RO_COX_BENZOATES)
                    num == 8 ? DERIV_BIT_Benzoate :
#endif
                    num == 9 ? DERIV_BIT_PFP :
                    num == 12 ? DERIV_BIT_HFB :
                    num == 13 ? DERIV_BIT_PFB :
                    DERIV_ID_Unknown;
            }
        }

        return 1;
    }
#endif

#if( defined(DERIV_RING_DMOX_DEOX_N) && defined(DERIV_RING_DMOX_DEOX_O) )
    if (!idrv && ( type == DERIV_RING_DMOX_DEOX_N ||
                   type == DERIV_RING_DMOX_DEOX_O ))
    {
        if (szUnderiv && type == DERIV_RING_DMOX_DEOX_N)
        {
            /* add only once; do not add upon DERIV_RING_DMOX_DEOX_O */
            underiv_list_add( szUnderiv, lenUnderiv, num == 4 ? "DMOX" : num == 6 ? "DEOX" : num ? "D?OX" : "???", ' ' );
            if (num == 4 || num == 6)
            {
                underiv_list_add( szUnderiv2, lenUnderiv2, pszDerivName[num == 4 ? DERIV_ID_DMOX : num == 6 ? DERIV_ID_DEOX : DERIV_ID_Unknown], ' ' );
                *bitUnderiv |= num == 4 ? DERIV_BIT_DMOX : num == 6 ? DERIV_BIT_DEOX : DERIV_BIT_Unknown;
            }
            else
            {
                underiv_list_add( szUnderiv2, lenUnderiv2, underiv_list_get_last( szUnderiv, ' ' ), ' ' );
                *bitUnderiv |= DERIV_BIT_Unknown;
            }
        }

        return 1;
    }
#endif

#ifdef DERIV_RING2_OUTSIDE_PRECUR
    if (type && ( type & DERIV_RING2_OUTSIDE_PRECUR ) == type)
    {
        if (szUnderiv && !idrv)
        {
            if (num == 4 || num == 5)
            {
                underiv_list_add( szUnderiv, lenUnderiv, num == 4 ? "Pyrrolidide" : num == 5 ? "Piperidine" : "???", ' ' );
                underiv_list_add( szUnderiv2, lenUnderiv2, pszDerivName[num == 4 ? DERIV_ID_Pyrrolidide : num == 5 ? DERIV_ID_Piperidine : DERIV_ID_Unknown], ' ' ); /* djb-rwth: addressing coverity ID #499491 -- working correctly for num == 5 */
                *bitUnderiv |= num == 4 ? DERIV_BIT_Pyrrolidide : num == 5 ? DERIV_BIT_Piperidine : DERIV_BIT_Unknown;
            }
            else
            {
#ifdef _DEBUG
                int stop_DERIV_RING2_OUTSIDE_PRECUR = 1; /* debug only */
#endif
                underiv_list_add( szUnderiv, lenUnderiv, num ? "Pyrrol?Piper?" : "???", ' ' );
                underiv_list_add( szUnderiv2, lenUnderiv2, pszDerivName[DERIV_ID_Unknown], ' ' );
                *bitUnderiv |= DERIV_BIT_Unknown;
            }
        }
        return 1;
    }
#endif

#ifdef DERIV_DANSYL
    if (!idrv && type && ( type & DERIV_DANSYL ) == type)
    {
        if (szUnderiv)
        {
            char szRO[16] = "R";
            strcat(szRO, at[start].elname);
            if (at[start].num_H == 1)
            {
                strcat(szRO, "H");
            }
            strcat(szRO, "-Dansyl");
            underiv_list_add( szUnderiv, lenUnderiv, szRO, ' ' );
            underiv_list_add( szUnderiv2, lenUnderiv2, pszDerivName[DERIV_ID_Dansyl], ' ' );
            *bitUnderiv |= DERIV_BIT_Dansyl;
        }
        return 1;
    }
#endif

    return 0;
}


/****************************************************************************/
int is_deriv_chain( inp_ATOM *at,
                    int start,
                    int num_atoms,
                    DERIV_AT *da1,
                    int idrv,
                    char *szUnderiv,
                    int lenUnderiv,
                    char *szUnderiv2,
                    int lenUnderiv2,
                    BIT_UNDERIV *bitUnderiv )
{
    return is_deriv_chain2( at, start, da1->typ[idrv], da1->num[idrv], da1->ord[idrv], idrv, szUnderiv, lenUnderiv, szUnderiv2, lenUnderiv2, bitUnderiv );
}


/****************************************************************************/
int is_deriv_chain_or_ring( inp_ATOM *at,
                            int start,
                            int num_atoms,
                            DERIV_AT *da1,
                            int *idrv )
{
    int i, ret = -1;
    if (da1->typ[*idrv] & DERIV_RING_OUTSIDE_PRECURSOR)
    {
        /* find the first ord of this derivative; ord of ring derivatives are in pairs */
        int j = -1;
        for (i = 0; i < DERIV_AT_LEN && da1->typ[i]; i++)
        {
            if (da1->typ[i] & DERIV_RING_OUTSIDE_PRECURSOR)
            {
                if (i == *idrv || i + 1 == *idrv)
                {
                    *idrv = j = i;
                    break;
                }
                i++; /* bypass the second bond to the same derivatization agent */
            }
        }
        /* check consistency */
        if (j == -1 || j + 1 >= DERIV_AT_LEN ||
             !( da1->typ[j] & DERIV_RING_OUTSIDE_PRECURSOR ) || !( da1->typ[j + 1] & DERIV_RING_OUTSIDE_PRECURSOR ))
        {
            ret = -1; /* program error */
        }
        else
        {
            ret = is_DERIV_RING_O_or_NH_OUTSIDE_PRECURSOR( at, start, num_atoms, da1, j, NULL, 0, NULL, 0, NULL );
        }
    }
    else
    {
        if (da1->typ[*idrv])
        {
            ret = is_DERIV_RING_O_or_NH_OUTSIDE_PRECURSOR( at, start, num_atoms, da1, *idrv, NULL, 0, NULL, 0, NULL );
        }
    }

    return ret;
}


#define DERIV_RING DERIV_RING_OUTSIDE_PRECURSOR


/****************************************************************************/
int remove_deriv( DERIV_AT *da1, int idrv )
{
    int i, j, ret = -1;

    if (da1->typ[idrv] & DERIV_RING)
    {
        /* find the first ord of this derivative; ord of ring derivatives are in pairs */
        j = -1;
        for (i = 0; i < DERIV_AT_LEN && da1->typ[i]; i++)
        {
            if (da1->typ[i] & DERIV_RING)
            {
                if (i == idrv || i + 1 == idrv)
                {
                    j = i;
                    break;
                }
                i++; /* bypass the second bond to the same derivatization agent */
            }
        }

        /* delete if data are consistent */
        if (j >= 0 && j + 1 < DERIV_AT_LEN && ( da1->typ[j] & DERIV_RING ) && ( da1->typ[j + 1] & DERIV_RING ))
        {
            for (; j < DERIV_AT_LEN - 2 && da1->typ[j + 2]; j++)
            {
                da1->typ[j] = da1->typ[j + 2];
                da1->num[j] = da1->num[j + 2];
                da1->ord[j] = da1->ord[j + 2];
            }
            for (; j < DERIV_AT_LEN; j++)
            {
                da1->typ[j] = 0;
                da1->num[j] = 0;
                da1->ord[j] = 0;
            }
            ret = 0;
        }
    }
    else
    {
        j = idrv;
        for (; j < DERIV_AT_LEN - 1 && da1->typ[j + 1]; j++)
        {
            da1->typ[j] = da1->typ[j + 1];
            da1->num[j] = da1->num[j + 1];
            da1->ord[j] = da1->ord[j + 1];
        }
        for (; j < DERIV_AT_LEN; j++)
        {
            da1->typ[j] = 0;
            da1->num[j] = 0;
            da1->ord[j] = 0;
        }
        ret = 0;
    }

    return ret;
}


/****************************************************************************/
int remove_deriv_mark( DERIV_AT *da1, int idrv )
{
    int i, j, ret = -1;
    if (da1->typ[idrv] & DERIV_RING)
    {
        /* find the first ord of this derivative; ord of ring derivatives are in pairs */
        j = -1;
        for (i = 0; i < DERIV_AT_LEN && da1->typ[i]; i++)
        {
            if (da1->typ[i] & DERIV_RING)
            {
                if (i == idrv || i + 1 == idrv)
                {
                    j = i;
                    break;
                }
                i++; /* bypass the second bond to the same derivatization agent */
            }
        }
        /* delete if data are consistent */
        if (j >= 0 && j + 1 < DERIV_AT_LEN && ( da1->typ[j] & DERIV_RING ) && ( da1->typ[j + 1] & DERIV_RING ))
        {
            da1->typ[j] |= DERIV_DUPLIC;
            da1->typ[j + 1] |= DERIV_DUPLIC;
            ret = 0;
        }
    }
    else
    {
        j = idrv;
        da1->typ[j] |= DERIV_DUPLIC;
        ret = 0;
    }

    return ret;
}


/****************************************************************************/
void underiv_buf_clear( char *szUnderiv )
{
    if (NULL != szUnderiv)
    {
        szUnderiv[0] = '\0';
    }
}

/****************************************************************************/
int underiv_list_add( char *szUnderivList, int lenUnderivList, const char *szUnderiv, char cDelimiter )
{
    if (NULL != szUnderivList && lenUnderivList > 0 && NULL != szUnderiv && *szUnderiv)
    {
        int lenList = strlen( szUnderivList );
        int lenAdd = strlen( szUnderiv );
        int bDelimiter = cDelimiter ? 1 : 0;
        if (lenList + lenAdd + bDelimiter < lenUnderivList)
        {
            if (lenList && bDelimiter)
            {
                szUnderivList[lenList++] = cDelimiter;
            }
            memcpy(szUnderivList + lenList, szUnderiv, (long long)lenAdd + 1); /* +1 adds zero termination */ /* djb-rwth: cast operator added */
            return lenList + lenAdd;
        }
    }

    return 0;
}


/****************************************************************************/
const char* underiv_list_get_last( const char *szUnderivList,
                                   char cDelimiter )
{
    if (szUnderivList && cDelimiter)
    {
        const char *p = strrchr( szUnderivList, cDelimiter );
        if (NULL == p)
        {
            p = szUnderivList;
        }
        return p;
    }

    return NULL;
}


/****************************************************************************/
int underiv_compare( const void *p1, const void *p2 )
{
    return strcmp( *(const char **) p1, *(const char **) p2 );
}


/****************************************************************************/
int underiv_list_add_two_cuts( char *szUnderivList,
                               int lenUnderivList,
                               char *szUnderiv,
                               const char cDelim )
{
    /* may happen only in RN- case, DERIV_AMINE_tN, but not always  */
    const char szDelim[] = { cDelim, 0 };
    char *p1 = strtok( szUnderiv, szDelim );
    char *p2 = p1 ? strtok( NULL, szDelim ) : NULL;
    char *p1m = p1 ? strchr( p1, '-' ) : NULL;
    char *p2m = p2 ? strchr( p2, '-' ) : NULL;
    char *pm;

    if (p1m && p2m)
    {
        if (p1m - p1 == p2m - p2 && !memcmp( p1, p2, p1m - p1 ))
        {
            /* found common prefix */
            int diff;
            *p1m++ = '\0';
            *p2m++ = '\0';
            /* output
            [common prefix]-(suffix1)+(suffix2): (suffix1)<(suffix2)  or
            [common prefix]-2(suffix1):          (suffix1)==(suffix2)
            */
            underiv_list_add( szUnderivList, lenUnderivList, p1, cDelim ); /* [common prefix] */
            underiv_list_add( szUnderivList, lenUnderivList, "-", 0 );  /* - */
            diff = strcmp( p1m, p2m );
            if (diff > 0)
            {
                pm = p1m;
                p1m = p2m; /* (suffix1) - smaller */
                p2m = pm;  /* (suffix2) - greater */
            }
            if (diff)
            {
                underiv_list_add( szUnderivList, lenUnderivList, p1m, 0 );  /* (suffix1) */
                underiv_list_add( szUnderivList, lenUnderivList, "+", 0 );  /* - */
                underiv_list_add( szUnderivList, lenUnderivList, p2m, 0 );  /* (suffix2) */
            }
            else
            {
                underiv_list_add( szUnderivList, lenUnderivList, "2", 0 );  /* 2 */
                underiv_list_add( szUnderivList, lenUnderivList, p1m, 0 );  /* (suffix1) */
            }
        }
        else
        {
            /* should not happen */
            underiv_list_add( szUnderivList, lenUnderivList, p1, cDelim );
            underiv_list_add( szUnderivList, lenUnderivList, p2, cDelim );
        }
    }
    else
    {
        if (p1 && p2)
        {
            /* should not happen */
            underiv_list_add( szUnderivList, lenUnderivList, p1, cDelim );
            underiv_list_add( szUnderivList, lenUnderivList, p2, cDelim );
        }
        else
        {
            if (p1)
            {
                /* happens in case only one num[] is not zero, namely, DERIV_RING2_OUTSIDE_PRECUR */
                underiv_list_add( szUnderivList, lenUnderivList, p1, cDelim );
            }
        }
    }

    return 0;
}


/****************************************************************************/
#if( UNDERIVATIZE_REPORT == 1 )
int sort_merge_underiv( char *pSdfValue,
                        int bOutputSdf,
                        char *szUnderivList,
                        char cDerivSeparator,
                        const char *pszUnderivPrefix,
                        const char *pszUnderivPostfix )
{
#define UNDERIV_MAX_NUM  512   /*max. number of records in szUnderivList */
    int num, numUnderiv = 0, i, j, k, m, n;
    char *q;
    char coeff[16];
    char *pszUnderiv[UNDERIV_MAX_NUM];

    num = strlen( pSdfValue );
    n = num + strlen( pszUnderivPrefix ) + 1;

    if (n + strlen( pszUnderivPostfix ) + 1 + strlen( szUnderivList ) < MAX_SDF_VALUE)
    {
        for (numUnderiv = 0, q = strtok( szUnderivList, " " ); numUnderiv < UNDERIV_MAX_NUM && q; q = strtok( NULL, " " ), numUnderiv++)
        {
            pszUnderiv[numUnderiv] = q;
        }
        /*if ( !bOutputSdf || num ) {*/
        n = underiv_list_add( pSdfValue, MAX_SDF_VALUE, pszUnderivPrefix, 0 );
        /*}*/
        if (numUnderiv > 1)
        {
            qsort( pszUnderiv, numUnderiv, sizeof( pszUnderiv[0] ), underiv_compare );
        }
        for (i = 0; i < numUnderiv; i = j)
        {
            for (j = i + 1; j < numUnderiv && !underiv_compare( pszUnderiv + i, pszUnderiv + j ); j++)
            {
                ; /* find identical derivatives */
            }
            k = strlen( pszUnderiv[i] );
            if (1 < ( num = j - i ))
            {
                k = sprintf(coeff, "%d", num);
            }
            if ((long long)n + (long long)k + sizeof( pszUnderivPostfix ) < MAX_SDF_VALUE) /* djb-rwth: cast operators added */
            {
                m = i ? cDerivSeparator : 0;
                if (num > 1)
                {
                    underiv_list_add( pSdfValue, MAX_SDF_VALUE, coeff, m );
                    m = 0;
                }
                n = underiv_list_add( pSdfValue, MAX_SDF_VALUE, pszUnderiv[i], m );
            }
            else
            {
                underiv_list_add( pSdfValue, MAX_SDF_VALUE, "!", 0 ); /* overflow indicator */
                numUnderiv = -( 1 + numUnderiv );
                break;
            }
        }
        if (!bOutputSdf || num)
        {
            underiv_list_add( pSdfValue, MAX_SDF_VALUE, pszUnderivPostfix, 0 );
        }
    }

    return numUnderiv;
#undef UNDERIV_MAX_NUM
}
#endif  /* UNDERIVATIZE_REPORT == 1 */


/****************************************************************************/
int eliminate_deriv_not_in_list( inp_ATOM *at,
                                 DERIV_AT *da,
                                 int num_atoms,
                                 char *szUnderivList,
                                 int lenUnderivList,
                                 char *szUnderivList2,
                                 int lenUnderivList2,
                                 BIT_UNDERIV *bitUnderivList )
{
    int i, j, num_da, num_cuts = 0, num_cuts_this_atom, ret = 0;
#if( UNDERIVATIZE_REPORT == 1 )
#define UNDERIV_LEN 128
#define UNDERIV_LEN2 128
    char szUnderiv[UNDERIV_LEN];
    char szUnderiv2[UNDERIV_LEN2];
    BIT_UNDERIV  bitUnderiv;
#else
#define UNDERIV_LEN 0
#define UNDERIV_LEN2 0
    char *szUnderiv = NULL;
    char *szUnderiv2 = NULL;
    BIT_UNDERIV  bitUnderiv;
#endif
    for (i = 0; i < num_atoms; i++)
    {
        if (!da[i].typ[0])
        {
            continue;
        }
        /* count deriative attachments */
        for (num_da = 0; num_da < DERIV_AT_LEN && da[i].typ[num_da]; num_da++)
        {
            ;
        }
        if (num_da > 2)
        {
            return -1; /* should not happen */
        }
        if (num_da == 2 && da[i].typ[0] != da[i].typ[1])
        {
            da[i].typ[0] = da[i].typ[1] = 0; /* do not allow */
            continue;
        }

        underiv_buf_clear( szUnderiv );
        underiv_buf_clear( szUnderiv2 );
        bitUnderiv = 0;

        num_cuts_this_atom = 0;
        if (da[i].typ[0] & DERIV_RING_OUTSIDE_PRECURSOR)
        {
            ret = 0;
#ifndef UNDERIV_SYLYL_ONLY
            if (num_da == 2 && 1 + da[i].num[0] + da[i].num[1] <= MAX_AT_DERIV &&
                 0 < ( ret = is_DERIV_RING_O_or_NH_OUTSIDE_PRECURSOR( at, i, num_atoms, da + i, 0, szUnderiv, UNDERIV_LEN, szUnderiv2, UNDERIV_LEN2, &bitUnderiv ) ))
            {
                num_cuts += 2;
                underiv_list_add( szUnderivList, lenUnderivList, szUnderiv, ' ' );
                underiv_list_add( szUnderivList2, lenUnderivList2, szUnderiv2, ' ' );
                *bitUnderivList |= bitUnderiv;
            }
            else
#endif
            {
                if (ret < 0)
                {
                    return ret;
                }
                else
                {
                    da[i].typ[0] = da[i].typ[1] = 0; /* not a derivative */
                }
            }
        }
        else
        {
            ret = 0;
            /*if ( da[i].num[0] <= MAX_AT_DERIV && 0 < (ret = is_deriv_chain( at, i, num_atoms, da+i, 0 )) )*/
            if (IS_DA_NUM_LE( da + i, 0, MAX_AT_DERIV ) && 0 < ( ret = is_deriv_chain( at, i, num_atoms, da + i, 0, szUnderiv, UNDERIV_LEN, szUnderiv2, UNDERIV_LEN2, &bitUnderiv ) ))
            {
                num_cuts++;
                num_cuts_this_atom++;
                j = 1;
                /*underiv_list_add( szUnderivList, lenUnderivList, szUnderiv, ' ' );*/
            }
            else
            {
                if (ret < 0)
                {
                    return ret;
                }
                else
                {
                    da[i].ord[0] = da[i].ord[1];
                    da[i].num[0] = da[i].num[1];
                    da[i].typ[0] = da[i].typ[1];
                    da[i].typ[1] = 0;
                    j = 0;
                }
            }
            /*underiv_buf_clear( szUnderiv );*/
            /*if ( da[i].typ[j] && da[i].num[j] <= MAX_AT_DERIV &&*/
            if (IS_DA_NUM_LE( da + i, j, MAX_AT_DERIV ) &&
                 0 < ( ret = is_deriv_chain( at, i, num_atoms, da + i, j, szUnderiv, UNDERIV_LEN, szUnderiv2, UNDERIV_LEN2, &bitUnderiv ) ))
            {
                num_cuts++;
                num_cuts_this_atom++;
                /*underiv_list_add( szUnderivList, lenUnderivList, szUnderiv, ' ' );*/
            }
            else
            {
                if (ret < 0)
                {
                    return ret;
                }
                else
                {
                    da[i].typ[j] = 0;
                }
            }
            if (num_cuts_this_atom == 2)
            {
                /* may happen only in RN- case, DERIV_AMINE_tN, but not always  */
                underiv_list_add_two_cuts( szUnderivList, lenUnderivList, szUnderiv, ' ' );
                underiv_list_add( szUnderivList2, lenUnderivList2, szUnderiv2, ' ' );
                *bitUnderivList |= bitUnderiv;
            }
            else
            {
                if (num_cuts_this_atom == 1)
                {
                    underiv_list_add( szUnderivList, lenUnderivList, szUnderiv, ' ' );
                    underiv_list_add( szUnderivList2, lenUnderivList2, szUnderiv2, ' ' );
                    *bitUnderivList |= bitUnderiv;
                }
            }
        }
    }

    return num_cuts;
}


/****************************************************************************/
int make_single_cut( inp_ATOM *at, DERIV_AT *da, int iat, int icut )
{
    int ret = -1; /* error flag */
    int iord = (int) da[iat].ord[icut]; /* ord of the bond in iat */

    if (da[iat].typ[icut] & DERIV_DUPLIC)
    {
        return 0;
    }

    if (iord < 0)
    {
        return -1; /* program error */
    }

    else
    {
        /* find other da[] that affect at[iat] */

        int jat = at[iat].neighbor[iord];  /* opposite atom */
        AT_NUMB *p = is_in_the_list( at[jat].neighbor, (AT_NUMB) iat, at[jat].valence );
        int jord = p ? ( p - at[jat].neighbor ) : -1;
        int i;
        const int iT = 2;
#ifdef UNDERIV_ADD_D_TO_PRECURSOR
        const int iD = 1;
#endif

        if (jord < 0)
        {
            return -1;  /* program error */
        }

        ret = DisconnectInpAtBond( at, NULL, iat, iord );

        if (ret == 1)
        {
            if (da[iat].typ[icut] & DERIV_RING_OUTSIDE_PRECURSOR)
            {
                /* at[jat] belongs to the main structure */
                at[jat].num_H++;
#ifdef UNDERIV_ADD_D_TO_PRECURSOR
                at[jat].num_iso_H[iD] ++; /* add D to the main structure */
#endif
                at[iat].num_H++;        /* add T to the derivatizing fragment */
                at[iat].num_iso_H[iT] ++;
            }
            else
            {
                if (da[iat].typ[icut] && ( da[iat].typ[icut] & DERIV_REPL_N_WITH_O ) == da[iat].typ[icut])
                {
                    at[jat].num_H++;
                    at[jat].num_iso_H[iT] ++; /* add T to the derivatizing fragment ??? */
                                              /* replace R=N-DerivAgent with R=O H-DerivAgent */
                    at[iat].elname[0] = 'O';
                    at[iat].el_number = EL_NUMBER_O;   /* since N replaced with O, do not add H */
                }
                else
                {
                    if (!icut && da[iat].typ[icut] && ( da[iat].typ[icut] & DERIV_REPL_N_WITH_OH ) == da[iat].typ[icut])
                    {
                        /* cut #0 is the second; in case of DERIV_RING2_OUTSIDE_PRECUR on the 1st cut H has already been added */
                        at[jat].num_H++;
                        at[jat].num_iso_H[iT] ++; /* add T to the derivatizing fragment ??? */
                                                  /* replace R=N-DerivAgent with R=O H-DerivAgent */
                        at[iat].elname[0] = 'O';
                        at[iat].el_number = EL_NUMBER_O;   /* since N replaced with O, do not add H */
#ifdef DERIV_RING2_OUTSIDE_PRECUR
                        if (!( da[iat].typ[icut] & DERIV_RING2_OUTSIDE_PRECUR ))
                        {
                            at[iat].num_H++;
                        }
#endif
                    }
                    else
                    {
                        at[iat].num_H++;
#ifdef UNDERIV_ADD_D_TO_PRECURSOR
                        at[iat].num_iso_H[iD] ++; /* add D to the main structure */
#endif
                        at[jat].num_H++;        /* add T to the derivatizing fragment */
                        at[jat].num_iso_H[iT] ++;
                    }
                }
            }

            /* adjust ord for other bonds */
            for (i = 0; i < DERIV_AT_LEN && da[iat].typ[i]; i++)
            {
                if (da[iat].ord[i] == iord)
                {
                    da[iat].ord[i] = -( 1 + da[iat].ord[i] ); /* mark the bond being disconnected */
                }
                else if (da[iat].ord[i] > iord)
                {
                    da[iat].ord[i] --;
                }
            }
            for (i = 0; i < DERIV_AT_LEN && da[jat].typ[i]; i++)
            {
                if (da[jat].ord[i] == jord)
                {
                    /* opposite atom needs the same bond to be disconnected */
#ifdef NEVER        /* not needed */
                    if (da[iat].num[icut] == da[jat].num[i])
                    {
                        iD = 2; /* mark both as fragments */
                    }
                    else
                        if (da[iat].num[icut] > da[jat].num[i])
                        {
                            iD = 2; /* greater as a main structure */
                            iT = 1; /* mark smaller as a derivatizing fragment */
                        }
#endif
                    da[jat].ord[i] = -( 1 + da[jat].ord[i] );
                    da[jat].typ[i] |= DERIV_DUPLIC; /* do not cut again */
                }
                else if (da[jat].ord[i] > jord)
                {
                    da[jat].ord[i] --;
                }

            }
        }
    }


    return ret;
}


/****************************************************************************/
int fill_out_bond_cuts( inp_ATOM *at,
                        DERIV_AT *da,
                        int num_atoms,
                        R2C_ATPAIR *ap,
                        int num_cuts_to_check )
{
    int i, j, k, n;
    AT_NUMB at1, at2;
    int ret = 0;
    /* fill out the array of bonds to be cut */
    for (i = j = 0; i < num_atoms; i++)
    {
        if (( da[i].typ[0] & DERIV_RING_OUTSIDE_PRECURSOR ) && ( da[i].typ[1] & DERIV_RING_OUTSIDE_PRECURSOR ) &&
             da[i].num[0] <= MAX_AT_DERIV && da[i].num[1] <= MAX_AT_DERIV)
        {
            if (j + 1 >= num_cuts_to_check)
            {
                ret = -2;
                goto exit_r2c_num; /* wrong number of cuts = num */
            }
            for (k = 0; k < 2; k++)
            {
                at1 = i;
                at2 = at[at1].neighbor[(int) da[at1].ord[k]];
                n = ( at1 > at2 );
                ap[j].at[n] = at1;
                ap[j].at[1 - n] = at2; /* ap[j].at[0] < ap[j].at[1] */
                ap[j].atno = i;
                j++;
            }
            if (0 < cmp_r2c_atpair( ap + j - 2, ap + j - 1 ))
            {
                R2C_ATPAIR ap1 = ap[j - 2];
                ap[j - 2] = ap[j - 1];
                ap[j - 1] = ap1; /* sort each pair */
            }
        }
        else
        {
            /* 2013-12-04 DT */
            if (da[i].typ[0] & DERIV_RING_DMOX_DEOX)
            {
                int other_atom = (int) da[i].other_atom - 1;
                if (da[i].typ[1] || other_atom < 0 || i == other_atom || da[other_atom].other_atom != i + 1 ||
                     !( da[other_atom].typ[0] & DERIV_RING_DMOX_DEOX ) || da[other_atom].typ[1])
                {
                    ret = -3;
                    goto exit_r2c_num;
                    /* no other cut may be at the atom in addition to DERIV_RING_DMOX_DEOX
                    or no other_atom or other_atom has wrong deriv. type */
                }
                /* make sure the ap[] for the two cuts in the same ring are adjacent */
                if (other_atom > i)
                {
                    if (j + 1 >= num_cuts_to_check)
                    {
                        ret = -2;
                        goto exit_r2c_num; /* wrong number of cuts = num */
                    }
                    /* cut #1 */
                    at1 = i;
                    at2 = at[at1].neighbor[(int) da[at1].ord[0]];
                    n = ( at1 > at2 );
                    ap[j].at[n] = at1;
                    ap[j].at[1 - n] = at2; /* ap[j].at[0] < ap[j].at[1] */
                    ap[j].atno = i;
                    j++;
                    /* cut #2 */
                    at1 = other_atom;
                    at2 = at[at1].neighbor[(int) da[at1].ord[0]];
                    n = ( at1 > at2 );
                    ap[j].at[n] = at1;
                    ap[j].at[1 - n] = at2; /* ap[j].at[0] < ap[j].at[1] */
                    ap[j].atno = i;
                    j++;
                }
                else
                {
                    /* add each pair of cuts only once */
                    continue;
                }
            }
#if ( COUNT_ALL_NOT_DERIV == 1 )
            else
            {
                for (k = 0; k < DERIV_AT_LEN && da[i].typ[k]; k++)
                {
                    if (j >= num_cuts_to_check || ( da[i].typ[k] & DERIV_RING_OUTSIDE_PRECURSOR ))
                    {
                        ret = -2;
                        goto exit_r2c_num; /* wrong number of cuts = num or wrong type */
                    }
                    at1 = i;
                    at2 = at[i].neighbor[(int) da[i].ord[k]];
                    n = ( at1 > at2 );
                    /* pair of atoms possibly to be disconnected */
                    ap[j].at[n] = at1;
                    ap[j].at[1 - n] = at2; /* ap[j].at[0] < ap[j].at[1] */
                                           /* precursor's atom */
                    ap[j].atno = i;
                    j++;
                }
            }
#endif
        }
    }
    return j;

exit_r2c_num:

    return ret;
}


/****************************************************************************/
int mark_deriv_agents( inp_ATOM *at,
                       DERIV_AT *da,
                       int num_atoms,
                       R2C_ATPAIR *ap,
                       int num_cuts_to_check,
                       AT_NUMB *pnum_comp,
                       int *pcur_num_at )
{
    /* mark components to be disconnected */
    int comp_num = 0;   /* number of components */
    int cur_num_at = 0; /* number of atoms left after disconnecting the derivatizing agent */
    int ret = 0;
    int i, j, k = -1, n;
    *pnum_comp = 0;
    *pcur_num_at = 0;
    UnMarkOtherIndicators( at, num_atoms );
    for (i = 0; i < num_cuts_to_check; i++)
    {
        n = 0;
        for (j = 0; j < 2; j++)
        {
            if (da[(int) ap[i].at[j]].typ[0])
            {
                k = j;
                n++;
            }
        }
        if (n != 1)
        {
            ret = -3;
            goto exit_r2c_num; /* wrong atom pair */
        }
        n = ap[i].at[k]; /* marked atom */
        if (( da[n].typ[0] & DERIV_RING_OUTSIDE_PRECURSOR ))
        {
            n = ap[i].at[1 - k];
        }
        /* at[n] belongs to the precursor */
        if (!at[n].at_type)
        {
            comp_num++;
            cur_num_at = mark_atoms_ap( at, n, ap, num_cuts_to_check, cur_num_at, comp_num );
        }
    }

    *pnum_comp = comp_num;
    *pcur_num_at = cur_num_at;

exit_r2c_num:

    return ret;
}


#ifdef FIX_UNDERIV_TO_SDF


/****************************************************************************
Input: at2[] has original bonds; at[] has normalized bonds
Description:
this function finds aromatic or other non-single-double-triple
bonds in at[] and replaces them with bonds between atoms, which has
the same orig_at_number, from at2[].
This tolerates permutation of atom locations in at[] because
orig_at_number do not change.
****************************************************************************/
int replace_arom_bonds( inp_ATOM *at,
                        int num_atoms,
                        inp_ATOM *at2,
                        int num_atoms2 )
{
    int i, j, num_err = 0;

    for (i = 0; i < num_atoms; i++)
    {
        for (j = 0; j < at[i].valence; j++)
        {
            if (at[i].bond_type[j] > BOND_TRIPLE)
            {
                /* find pairs of atoms using orig. atom numbers */
                int i1, i2;
                char bSuccess = 0;
                int neigh = at[i].neighbor[j];
                AT_NUMB orig_no1 = at[i].orig_at_number;
                AT_NUMB orig_no2 = at[neigh].orig_at_number;
                for (i1 = 0; i1 < num_atoms2 && at2[i1].orig_at_number != orig_no1; i1++)
                {
                    ;
                }
                for (i2 = 0; i2 < num_atoms2 && at2[i2].orig_at_number != orig_no2; i2++)
                {
                    ;
                }
                if (i1 < num_atoms2 && i2 < num_atoms2)
                {
                    AT_NUMB *p1 = is_in_the_list( at2[i1].neighbor, (AT_NUMB) i2, at[i1].valence );
                    AT_NUMB *pneigh = is_in_the_list( at[neigh].neighbor, (AT_NUMB) i, at[neigh].valence );
                    if (p1 && pneigh)
                    {
                        int n1 = p1 - at2[i1].neighbor;
                        int nneigh = pneigh - at[neigh].neighbor;
                        at[i].bond_type[j] = at[neigh].bond_type[nneigh] = at2[i1].bond_type[n1];
                        bSuccess = 1;
                    }
                }
                if (!bSuccess)
                {
#ifdef _DEBUG
                    int stop_here = 1;
#endif
                    num_err++;
                }
            }
        }
    }

    return num_err;
}
#endif  /* FIX_UNDERIV_TO_SDF */


#ifdef UNDERIV_ADD_EXPLICIT_H


/*****************************************************************************/
int add_explicit_H( INP_ATOM_DATA *inp_cur_data )
{
    /* do not care about stereo parities for now */
    int curRemovedH, num_added_explicit_H, iat, num_removed_H, m, num_H, num_atoms = inp_cur_data->num_at;
    if (( num_removed_H = inp_cur_data->num_removed_H ) > 0)
    {
        inp_ATOM *at_H = (inp_ATOM *) inchi_calloc( num_removed_H, sizeof( inp_ATOM ) );
        inp_ATOM *at = inp_cur_data->at;
        for (curRemovedH = num_atoms, num_added_explicit_H = 0; curRemovedH < num_atoms + num_removed_H; curRemovedH++)
        {
            if (at[curRemovedH].el_number == EL_NUMBER_H &&
                 1 == at[curRemovedH].valence &&
                 ( iat = (int) at[curRemovedH].neighbor[0] ) < num_atoms &&
                 0 <= ( m = at[curRemovedH].iso_atw_diff ) &&
                 m <= NUM_H_ISOTOPES)
            {
                /* num_H is the total number of all implicit H including isotopic H */
                if (at_H && 0 < at[iat].num_H && 0 <= ( num_H = at[iat].num_H - NUM_ISO_H( at, iat ) ) && ( (!m && num_H) || (m && at[iat].num_iso_H[m - 1]) )) /* djb-rwth: addressing LLVM warning; fixing a NULL pointer dereference */
                { /* number of implicit H > 0 */
                    int val = at[iat].valence;
                    /* set hydrogen atom */
                    at_H[num_added_explicit_H] = at[curRemovedH];
                    at_H[num_added_explicit_H].neighbor[0] = iat;
                    at_H[num_added_explicit_H].valence = 1;
                    at_H[num_added_explicit_H].chem_bonds_valence = at_H[num_added_explicit_H].bond_type[0];
                    /* set heavy atom */
                    at[iat].neighbor[val] = num_atoms + num_added_explicit_H;
                    at[iat].bond_type[val] = at_H[num_added_explicit_H].bond_type[0];
                    at[iat].bond_stereo[val] = -at_H[num_added_explicit_H].bond_stereo[0];
                    at[iat].valence++;
                    if (BOND_SINGLE <= at[iat].bond_type[val] && at[iat].bond_type[val] <= BOND_TRIPLE)
                    {
                        at[iat].chem_bonds_valence += at[iat].bond_type[val];
                    }
                    else
                    {
                        /* should not happen */
                        at_H[num_added_explicit_H].bond_type[0] = at[iat].bond_type[val] = BOND_SINGLE;
                        at[iat].chem_bonds_valence += BOND_SINGLE;
                        at_H[num_added_explicit_H].chem_bonds_valence = BOND_SINGLE;
                    }
                    at[iat].num_H--;
                    if (m)
                    {
                        at[iat].num_iso_H[m - 1] --;
                    }
                    num_added_explicit_H++;
                }
            }
        }
        if (0 < num_added_explicit_H && num_added_explicit_H <= num_removed_H)
        {
            memcpy(at + num_atoms, at_H, num_added_explicit_H * sizeof(at));
            inp_cur_data->num_removed_H = 0;
            inp_cur_data->num_at = ( num_atoms += num_added_explicit_H );
        }
        inchi_free( at_H );
    }
    return num_atoms;
}
#endif  /* UNDERIV_ADD_EXPLICIT_H */


/****************************************************************************
Main underivatization procedure
****************************************************************************/
int OAD_Edit_Underivatize( struct tagINCHI_CLOCK *ic,
                           struct tagCANON_GLOBALS *pCG,
                           ORIG_ATOM_DATA *orig_inp_data,
                           int bOutputSdf,
                           int bOutputReport,
                           char *pSdfValue )
{

#define ALLOC_AP \
        if ( 0 < num_cuts_to_check && (lenAllocated_ap < num_cuts_to_check || !ap) ) {\
            if ( ap )\
                inchi_free( ap );\
            ap = (R2C_ATPAIR *) inchi_malloc( num_cuts_to_check * sizeof(ap[0]) );\
            if ( !ap ) {\
                ret = -1;\
                goto exit_function; /* malloc failure */\
            }\
            lenAllocated_ap = num_cuts_to_check;\
        }

    int ret = 0, i, j, k, m, n, num_atoms, num_components, i_component, nFound, num, cur_num_at = 0, len, ind1, ind2, ind3; /* djb-rwth: adding variables for char -> int conversion of subscripts; initialisation added */
    int num_cuts, num_ring_cuts, num_cut_pieces, num_cuts_to_check;
    inp_ATOM *at = orig_inp_data->at; /* djb-rwth: ignoring LLVM warning: value used */
    INP_ATOM_DATA *inp_cur_data = NULL;
    DERIV_AT      *da = NULL;
    R2C_ATPAIR    *ap = NULL;
    int            lenAllocated_ap = 0;
    int  nTotNumCuts = 0;
    int  num_removed_H = 0; /* djb-rwth: ignoring LLVM warning: variable used */
#ifdef FIX_UNDERIV_TO_SDF
    inp_ATOM *at2 = NULL;
#endif
#if ( UNDERIVATIZE_REPORT == 1 )
#define UNDERIV_LIST_LEN 2048
#define UNDERIV_LIST_LEN2 2048
    char szUnderivList[UNDERIV_LIST_LEN] = "";
    char szUnderivList2[UNDERIV_LIST_LEN2] = "";
    char underivPrefix[] = "\tDeriv=";
    char underivPostfix[] = "";
    char underivPrefix2[] = "\tDeriv2=";
    char underivPostfix2[] = "";
    char underivPrefix3[] = "\tDerivBits=";
    char underivPostfix3[] = "";
    char cDerivSeparator = ',';
    int numUnderiv = 0, numUnderiv2 = 0, numUnderiv3 = 0; /* djb-rwth: ignoring LLVM warning: variables used to store function return value */
    BIT_UNDERIV bitUnderivList = 0;
    char szbitUnderivList[16]; /* int32 has at most 32/4=8 hexadecimal digits + 0x prefix + zero termination = 8+2+1=11 */
#else
#define UNDERIV_LIST_LEN 0
#define UNDERIV_LIST_LEN2 0
#define UNDERIV_MAX_NUM  0   /*max. number of records in szUnderivList */
    char *szUnderivList = NULL;
    char *szUnderivList2 = NULL;
    char **pszUnderiv = NULL;
    BIT_UNDERIV bitUnderivList = 0;
#endif


    /* Prepare */

    /*set_R2C_el_numbers( );*/

#ifndef UNDERIV_ADD_EXPLICIT_H
    num_atoms = remove_terminal_HDT( orig_inp_data->num_inp_atoms, at, 1 );
    /*^^^^^ always accomodate accomodate FIX_TERM_H_CHRG_BUG - IPl, July 2008*/
    num_removed_H = orig_inp_data->num_inp_atoms - num_atoms;
    orig_inp_data->num_inp_atoms = num_atoms;
#endif

    /* Initialize */
    UnMarkDisconnectedComponents( orig_inp_data );
    /* djb-rwth: removing redundant code */

    /* Mark */
    num_components = MarkDisconnectedComponents( orig_inp_data, 0 );
    inp_cur_data = (INP_ATOM_DATA *) inchi_calloc( num_components, sizeof( inp_cur_data[0] ) );

    for (i_component = 0; i_component < num_components; i_component++)
    {
        CreateInpAtomData(inp_cur_data + i_component, orig_inp_data->nCurAtLen[i_component], 0);

        inp_cur_data[i_component].num_at = ExtractConnectedComponent(orig_inp_data->at, orig_inp_data->num_inp_atoms, i_component + 1, inp_cur_data[i_component].at);

        /*  error processing */
        if (inp_cur_data[i_component].num_at <= 0 || orig_inp_data->nCurAtLen[i_component] != inp_cur_data[i_component].num_at)
        {
            ret = -(i_component + 1); /* severe error */
            goto exit_function;
        }

#ifdef UNDERIV_ADD_EXPLICIT_H
        num_atoms = remove_terminal_HDT(inp_cur_data[i_component].num_at, inp_cur_data[i_component].at, 1);
        inp_cur_data[i_component].num_removed_H = inp_cur_data[i_component].num_at - num_atoms;
        inp_cur_data[i_component].num_at = num_atoms;
#endif

        /* Initialize */
        num_atoms = inp_cur_data[i_component].num_at;
        at = inp_cur_data[i_component].at;
        add_DT_to_num_H(num_atoms, at);

        UnMarkRingSystemsInp(at, num_atoms);
        UnMarkOtherIndicators(at, num_atoms);
        UnMarkOneComponent(at, num_atoms);
        MarkRingSystemsInp(at, num_atoms, 0);
#ifdef FIX_UNDERIV_TO_SDF
        if (bOutputSdf)
        {
            /* save orig. bond types to restore them after replacing them with aromatic */
            if (at2)
            {
                inchi_free(at2);
            }
            if ((at2 = (inp_ATOM*)inchi_malloc(num_atoms * sizeof(inp_ATOM)))) /* djb-rwth: addressing LLVM warning */
            {
                memcpy(at2, at, num_atoms * sizeof(inp_ATOM));
            }
        }
#endif

        /* Mark aromatic bonds */
        ret = mark_arom_bonds(ic, pCG, at, num_atoms);
        if (ret < 0)
        {
            goto exit_function;
        }
        ret = 0;

        if (da)
        {
            inchi_free(da);
        }
        da = (DERIV_AT*)inchi_calloc(num_atoms, sizeof(da[0]));

        /* Detect derivatives */
        nFound = 0;
        for (i = 0; i < num_atoms; i++)
        {
            if (at[i].bCutVertex && !da[i].typ[inchi_min(at[i].valence, DERIV_AT_LEN) - 1])
            {
                for (k = 0; k < at[i].valence; k++)
                {
                    num = count_one_bond_atoms(at, da, i, k, CFLAG_MARK_BRANCH, &nFound);
                    UnMarkOtherIndicators(at, num_atoms);
                    if (num < 0)
                    {
                        ret = num; /* severe error */
                        goto exit_function;
                    }
                }
            }
        }

        /* Prepare cuts: remove cuts that are not to be done */
        /* in addition, count ring cuts DERIV_RING_OUTSIDE_PRECURSOR */
        num_ring_cuts = 0;
        num_cuts = 0;
        num_cut_pieces = 0;
        if (da) /* djb-rwth: fixing a NULL pointer dereference */
        {
            for (i = num = 0; i < num_atoms; i++) /* djb-rwth: ignoring LLVM warning: variable used */
            {
                /*for ( len = 0; len < MAX_AT_DERIV && da[i].typ[len]; len ++ ) -- bug fixed 2013-11-07 DCh */
                for (len = 0; len < DERIV_AT_LEN && da[i].typ[len]; len++)
                {
                    ;
                }
                switch (len)
                {

                case 0:
                    continue;

                case 1:
#if( defined(DERIV_RING_DMOX_DEOX_N) && defined(DERIV_RING_DMOX_DEOX_O) )
                    if (da[i].typ[0] & DERIV_RING_DMOX_DEOX)
                    {
                        if (!da[i].other_atom || (j = da[i].other_atom - 1) >= num_atoms ||
                            (da[j].typ[0] ^ DERIV_RING_DMOX_DEOX) != da[i].typ[0] ||
                            da[j].other_atom - 1 != i)
                        {
                            /* program error */
                            da[i].num[0] = NOT_AT_DERIV;
                        }
                        else
                        {
                            /* one of 2 ring cuts */
                            num_cuts += 1;
                            num_ring_cuts += 1;
                            num_cut_pieces += (i > j); /* count pieces only once */
                        }
                    }
                    else
#endif
                    {
                        /* single cut: unconditional */
                        num_cuts += 1;
                        num_cut_pieces += 1;
                    }
                    continue;

                case 2:
                    if ((da[i].typ[0] & DERIV_RING_OUTSIDE_PRECURSOR) && (da[i].typ[1] & DERIV_RING_OUTSIDE_PRECURSOR))
                    {
                        /* double cut, unconditional */
                        num_ring_cuts += 2;
                        num_cuts += 2;
                        num_cut_pieces += 1;
                        continue;
                    }
                    else
#ifdef DERIV_RING2_OUTSIDE_PRECUR
                        if (da[i].typ[0] && (da[i].typ[0] & DERIV_RING2_OUTSIDE_PRECUR) == da[i].typ[0] && da[i].typ[1] == da[i].typ[0])
                        {
                            /* double cut, unconditional */
                            num_ring_cuts += 2;
                            num_cuts += 2;
                            num_cut_pieces += 1;
                            continue;
                        }
                        else
#endif
                            if (da[i].typ[0] == DERIV_AMINE_tN && da[i].typ[1] == DERIV_AMINE_tN)
                            {
                                /* double cut, unconditional */
                                num_cuts += 2;
                                num_cut_pieces += 2;
                                continue;
                            }
                            else
#ifdef DERIV_RO_COX
                                if (da[i].typ[0] == DERIV_RO_COX && da[i].typ[1] == DERIV_RO_COX)
                                {
                                    if (da[i].num[0] == da[i].num[1])
                                    {
                                        memset(da + i, 0, sizeof(da[0])); /* don't remove if the two agents are identical */ /* djb-rwth: memset_s C11/Annex K variant? */
                                        continue;
                                    }
                                    else
                                    {
                                        if (da[i].num[0] && da[i].num[1])
                                        {
                                            static char pref_RO_COX[] = {/*likely deriv.agent*/12, 9, 6, 13, 3, 8/*likely precursor*/, 0 };
                                            char* p0 = strchr(pref_RO_COX, da[i].num[0]);
                                            char* p1 = strchr(pref_RO_COX, da[i].num[1]);
                                            if (p0 && p1)
                                            {
                                                j = p1 < p0; /* j=1 => deriv. agent num[1] has higher priority */
                                            }
                                            else
                                            {
                                                j = p0 ? 0 : p1 ? 1 : -1; /* we are here if there is a program error */
                                            }
                                        }
                                        else
                                        {
                                            j = da[i].num[0] ? 0 : da[i].num[1] ? 1 : -1; /* we are here if there is a program error */
                                        }
                                    }

                                    switch (j)
                                    {
                                    case 1:
                                        da[i].num[0] = da[i].num[1];
                                        da[i].ord[0] = da[i].ord[1];
                                        /* fall through */
                                    case 0:
                                        da[i].typ[1] = 0;
                                        da[i].num[1] = 0;
                                        da[i].ord[1] = 0;
                                        num_cuts += 1;
                                        num_cut_pieces += 1;
                                        continue;
                                    case -1:
                                        memset(da + i, 0, sizeof(da[0])); /* djb-rwth: memset_s C11/Annex K variant? */
                                        break; /* will produce error */
                                    }
                                }
#endif  /* DERIV_RO_COX */
                                else
                                    if (da[i].typ[0] == da[i].typ[1])
                                    {
                                        int sy0 = 0, sy1 = 0;
                                        /* DERIV_BRIDGE_O or DERIV_BRIDGE_NH ; cut off the smallest */
                                        if (0 == is_deriv_chain(at, i, num_atoms, da + i, 0, NULL, 0, NULL, 0, NULL))
                                        {
                                            da[i].num[0] = NOT_AT_DERIV;
                                        }
                                        else
                                        {
                                            ind1 = da[i].ord[0] - '0'; /* djb-rwth: converting char to int for subscript use */
                                            sy0 = is_silyl2(at, at[i].neighbor[ind1], i);
                                        }
                                        if (0 == is_deriv_chain(at, i, num_atoms, da + i, 1, NULL, 0, NULL, 0, NULL))
                                        {
                                            da[i].num[1] = NOT_AT_DERIV;
                                        }
                                        else
                                        {
                                            ind2 = da[i].ord[1] - '0'; /* djb-rwth: converting char to int for subscript use */
                                            sy1 = is_silyl2(at, at[i].neighbor[ind2], i);
                                        }
                                        if ((sy1 && (!sy0 || sy1 < sy0)) || (!(sy0 || sy1) && da[i].num[0] > da[i].num[1])) /* djb-rwth: addressing LLVM warnings */
                                        {
                                            da[i].num[0] = da[i].num[1];
                                            da[i].ord[0] = da[i].ord[1];
                                            da[i].typ[0] = da[i].typ[1];
                                            da[i].typ[1] = 0;
                                            num_cuts += 1;
                                            num_cut_pieces += 1;
                                        }
                                        else
                                        {
                                            if ((sy0 && (!sy1 || sy0 < sy1)) || (!(sy0 || sy1) && da[i].num[0] < da[i].num[1])) /* djb-rwth: addressing LLVM warning */
                                            {
                                                da[i].typ[1] = 0;
                                                num_cuts += 1;
                                                num_cut_pieces += 1;
                                            }
                                            else
                                            {
                                                /* attachments have same size: ignore both */
                                                /* ??? check for standard derivatizations ??? */
                                                da[i].typ[0] = 0;
                                                da[i].typ[1] = 0;
                                            }
                                        }
                                        continue;
                                    }
#ifdef DERIV_RO_COX
                                    else
                                        if ((da[i].typ[0] & (DERIV_RO_COX | DERIV_BRIDGE_O)) && (da[i].typ[1] & (DERIV_RO_COX | DERIV_BRIDGE_O)))
                                        {
                                            /*static char pref_RO_COX2[] = {13, 3, 8*/
                                            /*likely precursor*/
                                            /*, 0};*/
                                            j = (da[i].typ[0] == DERIV_BRIDGE_O);  /* da[i].typ[j] == DERIV_RO_COX */
                                            /*------------------------------------------------------------------------*
                                            * discard DERIV_RO_COX only in case [CH3-C(=O)]-O-[CH3]                  *
                                            *                      precursor DERIV_BRIDGE_O   DERIV_RO_COX precursor *
                                            * has already been done in finding DERIV_RO_COX                          *
                                            *------------------------------------------------------------------------*/
#ifdef NEVER
                                            /* methyl/ethyl alcoholes have already been excluded in get_traversed_deriv_type(...)
                                            for etyl/methyl acetate/benzoate. Only acetate/benzoic acids may their precursors */
                                            if (da[i].num[j] == 3 &&  /* -O-[C(=O)-CH3]  : DERIV_RO_COX*/
                                                da[i].num[1 - j] == 1)
                                            { /* [CH3]-O-C(=O)-R : DERIV_BRIDGE_O; derivatizin agent in [] */
                                                j = 1 - j; /* remove da[i].xxx[j] */
                                        }
#endif
                                            /*
                                            O
                                            ||
                                            R----O----C---X

                                            R----O---OC---X
                                            DERIV_BRIDGE_O \___/\____precursor__/

                                            R----O----COX
                                            \___precursor__/ \___/ DERIV_RO_COX

                                            rule: If R is >Si< then select DERIV_BRIDGE_O

                                            */

                                            if (j /* choose DERIV_RO_COX */)
                                            {
                                                /* da[i].typ[1]=DERIV_RO_COX is not likely a precursor */
                                                da[i].typ[0] = da[i].typ[1];
                                                da[i].num[0] = da[i].num[1];
                                                da[i].ord[0] = da[i].ord[1];
                                            }
                                            da[i].typ[1] = 0;
                                            da[i].num[1] = 0;
                                            da[i].ord[1] = 0;

                                            num_cuts += 1;
                                            num_cut_pieces += 1;
                                            continue;
                }
#endif
                    ret = -88;
                    goto exit_function; /* unexpected */

                case 3:
                    if (da[i].typ[0] == da[i].typ[1] &&
                        da[i].typ[0] == da[i].typ[2] &&
                        da[i].typ[0] == DERIV_AMINE_tN)
                    {
                        int x, y, z;
                        int sy[3] = { 0, 0, 0 }; /* silyl */

                        if (0 == is_deriv_chain(at, i, num_atoms, da + i, 0, NULL, 0, NULL, 0, NULL))
                        {
                            da[i].num[0] = NOT_AT_DERIV;
                        }
                        else
                        {
                            ind1 = da[i].ord[0] - '0'; /* djb-rwth: converting char to int for subscript use */
                            sy[0] = is_silyl2(at, at[i].neighbor[ind1], i);
                        }

                        if (0 == is_deriv_chain(at, i, num_atoms, da + i, 1, NULL, 0, NULL, 0, NULL))
                        {
                            da[i].num[1] = NOT_AT_DERIV;
                        }
                        else
                        {
                            ind2 = da[i].ord[1] - '0'; /* djb-rwth: converting char to int for subscript use */
                            sy[1] = is_silyl2(at, at[i].neighbor[ind2], i);
                        }

                        if (0 == is_deriv_chain(at, i, num_atoms, da + i, 2, NULL, 0, NULL, 0, NULL))
                        {
                            da[i].num[2] = NOT_AT_DERIV;
                        }
                        else
                        {
                            ind3 = da[i].ord[2] - '0'; /* djb-rwth: converting char to int for subscript use */
                            sy[2] = is_silyl2(at, at[i].neighbor[ind3], i);
                        }

                        /* djb-rwth: removing redundant code */

                        x = ((sy[0] && (!sy[1] || sy[0] < sy[1])) || (!(sy[0] || sy[1]) && da[i].num[0] < da[i].num[1])) ? 0 : 1; /* djb-rwth: addressing LLVM warning */
                        z = !x;
                        x = ((sy[x] && (!sy[2] || sy[x] < sy[2])) || (!(sy[x] || sy[2]) && da[i].num[x] < da[i].num[2])) ? x : 2; /* min */ /* djb-rwth: addressing LLVM warning */
                        /*z = (da[i].num[0] < da[i].num[1])? 1 : 0;*/
                        /*z = (da[i].num[x] < da[i].num[2])? 2 : z;*/ /* max */
                        z = ((sy[x] && (!sy[2] || sy[x] < sy[2])) || (!(sy[x] || sy[2]) && da[i].num[x] < da[i].num[2])) ? 2 : z; /* max */ /* djb-rwth: addressing LLVM warning */
                        y = ((x + 1) ^ (z + 1)) - 1;                      /* median */

                        if (da[i].num[x] == da[i].num[z] && sy[x] == sy[z])
                        {
                            /* no cuts */
                            da[i].typ[0] = 0;
                            da[i].typ[1] = 0;
                            da[i].typ[2] = 0;
                            continue; /* all deriv. agents have same size, ignore */
                            /* ??? check for standard derivatizations ??? */
                        }
                        else
                        {
                            if ((da[i].num[x] == da[i].num[y] && sy[x] == sy[y]) || (sy[x] && sy[y] && !sy[z])) /* djb-rwth: addressing LLVM warning */
                            {
                                /* two smallest */
                                switch (z)
                                {
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
                                num_cuts += 2;
                                num_cut_pieces += 2;
                            }
                            else
                            {
                                /* one smallest */
                                if (x)
                                {
                                    da[i].num[0] = da[i].num[x];
                                    da[i].ord[0] = da[i].ord[x];
                                    da[i].typ[0] = da[i].typ[x];
                                }
                                da[i].typ[1] = 0;
                                da[i].typ[2] = 0;
                                num_cuts += 1;
                                num_cut_pieces += 1;
                            }
                        }
                        continue;
                    }
                    ret = -88;
                    goto exit_function; /* unexpected */

                case 4:
                    if ((da[i].typ[0] & DERIV_RING_OUTSIDE_PRECURSOR) && (da[i].typ[1] & DERIV_RING_OUTSIDE_PRECURSOR) &&
                        (da[i].typ[2] & DERIV_RING_OUTSIDE_PRECURSOR) && (da[i].typ[3] & DERIV_RING_OUTSIDE_PRECURSOR))
                    {
                        int n01 = inchi_max(da[i].num[0], da[i].num[1]);
                        int n23 = inchi_max(da[i].num[2], da[i].num[3]);
                        if (n01 < n23 && 0 < is_DERIV_RING_O_or_NH_OUTSIDE_PRECURSOR(at, i, num_atoms, da + i, 0, NULL, 0, NULL, 0, NULL))
                        {
                            da[i].typ[2] = 0;
                            da[i].typ[3] = 0;
                            num_cuts += 2;
                            num_ring_cuts += 2;
                            num_cut_pieces += 1;
                        }
                        else
                        {
                            if (n01 > n23 && 0 < is_DERIV_RING_O_or_NH_OUTSIDE_PRECURSOR(at, i, num_atoms, da + i, 2, NULL, 0, NULL, 0, NULL))
                            {
                                da[i].num[0] = da[i].num[2];
                                da[i].ord[0] = da[i].ord[2];
                                da[i].typ[0] = da[i].typ[2];

                                da[i].num[1] = da[i].num[3];
                                da[i].ord[1] = da[i].ord[3];
                                da[i].typ[1] = da[i].typ[3];

                                da[i].typ[2] = 0;
                                da[i].typ[3] = 0;
                                num_cuts += 2;
                                num_ring_cuts += 2;
                                num_cut_pieces += 1;
                            }
                            else
                            {
                                da[i].typ[0] = 0;
                                da[i].typ[1] = 0;
                                da[i].typ[2] = 0;
                                da[i].typ[3] = 0;
                            }
                        }
                        continue;
                    }
                    ret = -88;
                    goto exit_function; /* unexpected */
            }
        }

            /*
            Eliminate cases when
            da[i1].typ[j1] && da[i2].typ[j2] &&
            at[i1].neighbor[da[i1].ord[j1]] == i2 && at[i2].neighbor[da[i2].ord[j2]] == i1
            that is, same bond is in the da[] elements of the adjacent atoms
            */

            nFound = 0; /* number of cuts to remove */
            for (i = 0; i < num_atoms; i++)
            {
                /*for ( j = 0; j < MAX_AT_DERIV && da[i].typ[j]; j ++ ) -- bug fixed 2013-11-07 DCh */
                for (j = 0; j < DERIV_AT_LEN && da[i].typ[j]; j++)
                {
                    if (da[i].typ[j] & DERIV_DUPLIC)
                    {
                        continue;
                    }
                    n = at[i].neighbor[(int)da[i].ord[j]];
                    if (n < i)
                    {
                        continue;
                    }
                    /*for ( k = 0; k < MAX_AT_DERIV && da[n].typ[k]; k ++ ) -- bug fixed 2013-11-07 DCh */
                    for (k = 0; k < DERIV_AT_LEN && n< num_atoms && da[n].typ[k]; k++)
                    {
                        if (da[n].typ[k] & DERIV_DUPLIC)
                        {
                            continue;
                        }
                        m = at[n].neighbor[(int)da[n].ord[k]];
                        if (m == i)
                        {
                            /* same bond in da[i].typ[j] and da[n].typ[k] */
                            /* check whether both derivatives are acceptable */
                            int k1 = k, j1 = j;
                            int ret_i = is_deriv_chain_or_ring(at, i, num_atoms, da + i, &j1);
                            int ret_n = is_deriv_chain_or_ring(at, n, num_atoms, da + n, &k1);
                            if (ret_i < 0)
                            {
                                ret = ret_i;
                                goto exit_function;
                            }
                            if (ret_n < 0)
                            {
                                ret = ret_n;
                                goto exit_function;
                            }
                            if (!ret_i || (ret_i && ret_n)) /* djb-rwth: addressing LLVM warning */
                            {
                                if (da[i].typ[j1] & DERIV_RING_OUTSIDE_PRECURSOR)
                                {
                                    num_cuts -= 2;
                                    num_ring_cuts -= 2;
                                }
                                else
                                {
                                    num_cuts -= 1;
                                }
                                num_cut_pieces -= 1;
                                if ((ret = remove_deriv_mark(da + i, j1))) /* djb-rwth: addressing LLVM warning */
                                {
                                    goto exit_function;
                                }
                                nFound++;
                            }
                            if (!ret_n || (ret_i && ret_n)) /* djb-rwth: addressing LLVM warning */
                            {
                                if (da[n].typ[k1] & DERIV_RING_OUTSIDE_PRECURSOR)
                                {
                                    num_cuts -= 2;
                                    num_ring_cuts -= 2;
                                }
                                else
                                {
                                    num_cuts -= 1;
                                }
                                num_cut_pieces -= 1;
                                if ((ret = remove_deriv_mark(da + n, k1))) /* djb-rwth: addressing LLVM warning */
                                {
                                    goto exit_function;
                                }
                                nFound++;
                            }
                        }
                    }
                }
            }

            if (nFound)
            {
                for (i = 0; i < num_atoms; i++)
                {
                    /*for ( j = 0; j < MAX_AT_DERIV && da[i].typ[j]; j ++ ) -- bug fixed 2013-11-07 DCh */
                    for (j = 0; j < DERIV_AT_LEN && da[i].typ[j]; j++)
                    {
                        /* attn: j is changed inside the cycle body */
                        if (da[i].typ[j] & DERIV_DUPLIC)
                        {
                            if ((ret = remove_deriv(da + i, j))) /* djb-rwth: addressing LLVM warning */
                            {
                                goto exit_function;
                            }
                            j--;
                        }
                    }
                }
            }
        }

        /* make sure DERIV_RING_OUTSIDE_PRECURSOR type disconnections increase */
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
        /* count DERIV_RING_OUTSIDE_PRECURSOR-type attachments */

#if ( COUNT_ALL_NOT_DERIV == 1 )
        num_cuts_to_check = num_cuts; /* STOP HERE */
#else
        num_cuts_to_check = num_ring_cuts;
#endif
        ALLOC_AP

            if (num_cuts_to_check >= 2)
            {
                /* check */
                AT_NUMB    comp_num;
                int        /*n,*/ m_at, m_ord;
                /*AT_NUMB at1, at2;*/

            repeat_without_deriv_ring:

                ALLOC_AP

                    comp_num = 0;
                /* fill out the array of bonds to be cut */
                j = fill_out_bond_cuts( at, da, num_atoms, ap, num_cuts_to_check );
                if (j < 0)
                {
                    ret = j;
                    goto exit_r2c_num; /* wrong number of cuts = num */
                }
                if (j != num_cuts_to_check)
                {
                    ret = -3;
                    goto exit_r2c_num; /* wrong number of cuts = num */
                }
                /* can't sort the bonds for subsequent searching by bisections in
                mark_atoms_ap() -> has_atom_pair() 2013-08-48 DCh */
                /* !!!!!!!! check that there are no derivatives inside a derivative */
                comp_num = 0; /* here it is the number of removed cuts */
                for (i = 0; i < num_cuts_to_check; i += j)
                {
                    for (j = n = 0; j < 2; j++)
                    {
                        int atj = (int) ap[i].at[j];
                        if (da[atj].typ[0] && at[atj].neighbor[(int) da[atj].ord[0]] == ap[i].at[1 - j])
                        {
                            k = j;      /* ap[i].at[k] is precursor atom */
                            n++;
                            m_at = atj; /* precursor atom at[m_at], da[m_at] */
                            m_ord = 0;  /* da[m_at].typ[m_ord] - type of the deriv.bond to break  */
                        }
                        else
                        {
                            if (da[atj].typ[1] && at[atj].neighbor[(int) da[atj].ord[1]] == ap[i].at[1 - j])
                            {
                                k = j;
                                n++;
                                m_at = atj;
                                m_ord = 1;
                            }
                        }
                    }

                    if (n != 1)
                    {
                        ret = -3;
                        goto exit_r2c_num; /* wrong atom pair */
                    }

                    if (( da[m_at].typ[m_ord] & DERIV_RING_OUTSIDE_PRECURSOR ))
                    {
                        n = (int) ap[i].at[k];   /* atom inside the derivation attachment */
                        j = 2;             /* number of bonds to cut */
                        if (i + j > num_cuts_to_check || ((int) ap[i + 1].at[0] != n && (int) ap[i + 1].at[1] != n)) /* djb-rwth: addressing LLVM warning */
                        {
                            ret = -3;
                            goto exit_r2c_num; /* wrong atom pair */
                        }
                    }
                    else
                    {
                        n = ap[i].at[1 - k]; /* atom inside the tentative derivation attachment */
                        j = 1;             /* number of bonds to cut */
#if( defined(DERIV_RING_DMOX_DEOX_N) && defined(DERIV_RING_DMOX_DEOX_O) )
                                           /*j += (0 != (da[m_at].typ[m_ord] & DERIV_RING_DMOX_DEOX)); */
                                           /* these 2 cuts are always adjacent */
                        j += ( da[m_at].typ[m_ord] && da[m_at].typ[m_ord] == ( da[m_at].typ[m_ord] & DERIV_RING_DMOX_DEOX ) ); /* these 2 cuts are always adjacent */
#endif
#ifdef DERIV_RING2_OUTSIDE_PRECUR
                        j += ( da[m_at].typ[m_ord] && da[m_at].typ[m_ord] == ( da[m_at].typ[m_ord] & DERIV_RING2_OUTSIDE_PRECUR ) ); /* these 2 cuts are always adjacent */
#endif
                    }

                    /* at[n] belongs to the derivation agent  */
                    cur_num_at = mark_atoms_ap( at, n, ap + i, j, 0, 1 );
                    for (k = 0; k < num_cuts_to_check; k++)
                    {
                        if (k == i || k == i + j - 1) /* not more than 2 cuts per derivatization agent */
                        {
                            continue; /* skip current 1 or 2 cuts */
                        }
                        if (at[(int) ap[k].at[0]].at_type || at[(int) ap[k].at[1]].at_type)
                        {
                            /*#ifdef DERIV_X_OXIME*/
                            if (( da[m_at].typ[m_ord] & DERIV_UNEXPADABLE ) == da[m_at].typ[m_ord])
                            {
                                /* this derivatization agent cannot be inside another, larger deriv. agent */
                                /* disable cuts of the larger derivatization agent */
                                if (ap[k].atno != m_at)
                                {
                                    for (j = 0; j < 2 && da[ap[k].atno].typ[j]; j++)
                                    {
                                        da[ap[k].atno].typ[j] |= DERIV_UNMARK;
                                    }
                                    num_cuts -= 1;
                                    num_cut_pieces -= 1;
                                    if (j == 2)
                                    {
                                        num_cuts -= 1;
                                        num_ring_cuts -= 2;
                                    }
                                    comp_num++;
                                }
                                /* djb-rwth: removing redundant code */
                            }
                            else
                                /* #endif */
                                /* unmark the cut: found a cut inside the derivatizing agent */
                            {
                                da[m_at].typ[m_ord] |= DERIV_UNMARK;
                                num_cuts -= 1;
                                num_cut_pieces -= 1;
                                if (j == 2)
                                {
                                    da[m_at].typ[1 - m_ord] |= DERIV_UNMARK;
                                    num_cuts -= 1;
                                    num_ring_cuts -= 2;
                                }
                                comp_num++;
                            }
                            break;
                        }
                    }
                    UnMarkOtherIndicators( at, num_atoms );
                }

                if (comp_num)
                {
                    for (i = 0; i < num_atoms; i++)
                    {
                        if (da[i].typ[0] & DERIV_UNMARK)
                        {
                            da[i].num[0] = da[i].num[1];
                            da[i].ord[0] = da[i].ord[1];
                            da[i].typ[0] = da[i].typ[1];
                            da[i].typ[1] = 0;
                            j = 0;
                        }
                        else
                        {
                            j = 1;
                        }
                        if (da[i].typ[j] & DERIV_UNMARK)
                        {
                            da[i].typ[j] = 0;
                        }
                    }

#if ( COUNT_ALL_NOT_DERIV == 1 )
                    num_cuts_to_check = num_cuts;
#else
                    num_cuts_to_check = num_ring_cuts;
#endif
                    if (num_cuts < 0 || num_ring_cuts < 0 || num_cut_pieces < 0)
                    {
                        ret = -3;
                        goto exit_r2c_num; /* wrong number of cuts = num */
                    }
#if( defined(DERIV_X_OXIME) || defined(DERIV_RO_COX) )
                    if (num_cuts_to_check > 0)
#endif
                    {
                        goto repeat_without_deriv_ring;
                    }
                }

                /* Sort the bonds for subsequent searching by bisections
                -- disabled because DERIV_RING_DMOX_DEOX have to be adjacent in ap
                if ( num_cuts_to_check > 1 ) {
                qsort( ap, num_cuts_to_check, sizeof(ap[0]), cmp_r2c_atpair);
                }
                */
                if ((ret = mark_deriv_agents( at, da, num_atoms, ap, num_cuts_to_check, &comp_num, &cur_num_at ))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_r2c_num; /* wrong atom pair */
                }
                if (comp_num > 1)
                {
                    /* eliminate offending DERIV_RING_OUTSIDE_PRECURSOR type derivatives */
                    if (num_ring_cuts <= 2)
                    {
                        ret = -99;
                        goto exit_r2c_num;
                    }
                    n = 0;
                    for (i = j = 0; i < num_atoms; i++) /* djb-rwth: ignoring LLVM warning: variable used */
                    {
                        if (( da[i].typ[0] & DERIV_RING_OUTSIDE_PRECURSOR ) && ( da[i].typ[1] & DERIV_RING_OUTSIDE_PRECURSOR ))
                        {
                            int at1a = at[i].neighbor[(int) da[i].ord[0]];
                            int at2a = at[i].neighbor[(int) da[i].ord[1]];
                            if (at[at1a].at_type != at[at2a].at_type)
                            {
                                da[i].typ[0] = 0; /* eliminate this cut */
                                da[i].typ[1] = 0;
                                n++;
                                num_cuts_to_check -= 2;
                                num_cuts -= 2;
                                num_ring_cuts -= 2;
                                num_cut_pieces -= 1;
                            }
                        }
                    }
                    if (n > 0 && num_cuts_to_check > 2)
                    {
                        UnMarkOtherIndicators( at, num_atoms );
                        goto repeat_without_deriv_ring;
                    }
                }
                ret = 0;

            exit_r2c_num:
                /*inchi_free( ap );*/
                UnMarkOtherIndicators( at, num_atoms );
                /*if ( ret < 0 || num_cuts_to_check >= 2 && cur_num_at < MIN_AT_LEFT_DERIV ) */
                /* -- bug: cur_num_at may include later rejected deriv. agents 2013-11-08 DCh */
                if (ret < 0)
                {
                    goto exit_function; /* unexpected  error or nothing left */
                }
            }

        if (!num_cuts)
        {
            continue; /*goto exit_function;*/
        }

        /* Eliminate derivatives that are not in the list */
        num_cuts = eliminate_deriv_not_in_list( at, da, num_atoms, szUnderivList, UNDERIV_LIST_LEN, szUnderivList2, UNDERIV_LIST_LEN2, &bitUnderivList );
        if (num_cuts < 0)
        {
            ret = num_cuts;
            goto exit_function;
        }

        /* Check how many atoms was left in the precursor - begin - 2013-11-12 DT */
        if (( num_cuts_to_check = num_cuts ) >= 1)
        {
            AT_NUMB comp_num = 0; /* here it is the number of removed cuts */
            int        /*n,*/ m_at, m_ord;

            ALLOC_AP

                /* Fill out the array of bonds to be cut */
                j = fill_out_bond_cuts( at, da, num_atoms, ap, num_cuts_to_check );
            if (j < 0)
            {
                ret = j;
                goto exit_r2c_num2; /* wrong number of cuts = num */
            }
            if (j != num_cuts_to_check)
            {
                ret = -3;
                goto exit_r2c_num2; /* wrong number of cuts = num */
            }

            /* The following code was copied from above to mark removed */
            for (i = 0; i < num_cuts_to_check; i += j)
            {
                for (j = n = 0; j < 2; j++)
                {
                    int atj = (int) ap[i].at[j];
                    if (atj < num_atoms && da[atj].typ[0] && at[atj].neighbor[(int) da[atj].ord[0]] == ap[i].at[1 - j])
                    {
                        k = j;      /* ap[i].at[k] is precursor atom */
                        n++;
                        m_at = atj; /* precursor atom at[m_at], da[m_at] */
                        m_ord = 0;  /* da[m_at].typ[m_ord] - type of the deriv.bond to break  */
                    }
                    else
                    {
                        if (atj < num_atoms && da[atj].typ[1] && at[atj].neighbor[(int) da[atj].ord[1]] == ap[i].at[1 - j])
                        {
                            k = j;
                            n++;
                            m_at = atj;
                            m_ord = 1;
                        }
                    }
                }
                if (n != 1)
                {
                    ret = -3;
                    goto exit_r2c_num2; /* wrong atom pair */
                }
                if (( da[m_at].typ[m_ord] & DERIV_RING_OUTSIDE_PRECURSOR ))
                {
                    n = (int) ap[i].at[k];   /* atom inside the derivation attachment */
                    j = 2;             /* number of bonds to cut */
                    if (i + j > num_cuts_to_check || ((int) ap[i + 1].at[0] != n && (int) ap[i + 1].at[1] != n)) /* djb-rwth: addressing LLVM warning */
                    {
                        ret = -3;
                        goto exit_r2c_num2; /* wrong atom pair */
                    }
                }
                else
                {
                    n = ap[i].at[1 - k];    /* atom inside the tentative derivation attachment */
                    j = 1;                  /* number of bonds to cut */
#if( defined(DERIV_RING_DMOX_DEOX_N) && defined(DERIV_RING_DMOX_DEOX_O) )
                                            /*j += (0 != (da[m_at].typ[m_ord] & DERIV_RING_DMOX_DEOX));*/
                                            /* these 2 cuts are always adjacent */
                    j += ( da[m_at].typ[m_ord] && da[m_at].typ[m_ord] == ( da[m_at].typ[m_ord] & DERIV_RING_DMOX_DEOX ) ); /* these 2 cuts are always adjacent */
#endif
#ifdef DERIV_RING2_OUTSIDE_PRECUR
                    j += ( da[m_at].typ[m_ord] && da[m_at].typ[m_ord] == ( da[m_at].typ[m_ord] & DERIV_RING2_OUTSIDE_PRECUR ) ); /* these 2 cuts are always adjacent */
#endif
                }

                /* at[n] belongs to the derivation agent  */
                cur_num_at = mark_atoms_ap( at, n, ap + i, j, 0, 1 );
                UnMarkOtherIndicators( at, num_atoms );
            }
            ret = mark_deriv_agents( at, da, num_atoms, ap, num_cuts_to_check, &comp_num, &cur_num_at );
            UnMarkOtherIndicators( at, num_atoms );
            if (ap)
            {
                inchi_free( ap );
                ap = NULL;
            }
            if (ret)
            {
                goto exit_r2c_num2; /* wrong atom pair */
            }
        }

    exit_r2c_num2:

        if (ap)
        {
            inchi_free( ap );
            ap = NULL;
        }
        if (ret < 0 || (num_cuts_to_check >= 2 && cur_num_at < MIN_AT_LEFT_DERIV)) /* -- bug: cur_num_at may include later rejected deriv. agents 2013-11-08 DCh */ /* djb-rwth: addressing LLVM warning */
        {
            goto exit_function; /* unexpected  error or nothing left */
        }

        /* Check how many atoms was left in the precursor - end - 2013-11-12 DT */

        /* make cuts */
        num_cuts = 0;
        for (i = num = 0; i < num_atoms; i++) /* djb-rwth: ignoring LLVM warning: variable used */
        {
            /*for ( len = 0; len < MAX_AT_DERIV && da[i].typ[len]; len ++ ) -- bug fixed 2013-11-07 DCh */
            for (len = 0; len < DERIV_AT_LEN && da[i].typ[len]; len++)
            {
                ;
            }
            switch (len)
            {
                case 0:
                    continue;
                case 1:
                    /* single cut: unconditional */
                    make_single_cut( at, da, i, 0 );
                    num_cuts += 1;
                    continue;
                case 2:
                    if (( (da[i].typ[0] & DERIV_RING_OUTSIDE_PRECURSOR ) && ( da[i].typ[1] & DERIV_RING_OUTSIDE_PRECURSOR )) ||
                         (da[i].typ[0] == DERIV_AMINE_tN && da[i].typ[1] == DERIV_AMINE_tN)
#ifdef DERIV_RING2_OUTSIDE_PRECUR
                         || (da[i].typ[0] && da[i].typ[0] == ( da[i].typ[0] & DERIV_RING2_OUTSIDE_PRECUR ) &&
                         da[i].typ[1] == da[i].typ[0])
#endif
                         ) /* djb-rwth: addressing LLVM warning */
                    {
                        /* double cut, unconditional */
                        make_single_cut( at, da, i, 1 );
                        make_single_cut( at, da, i, 0 );
                        num_cuts += 1;
                        continue;
                    }
                    if (da[i].typ[0] == da[i].typ[1])
                    {
                        /* DERIV_BRIDGE_O or DERIV_BRIDGE_NH; cut off the smallest */
                        if (da[i].num[0] > da[i].num[1])
                        {
                            make_single_cut( at, da, i, 1 );
                            num_cuts += 1;
                        }
                        else
                        {
                            if (da[i].num[0] < da[i].num[1])
                            {
                                make_single_cut( at, da, i, 0 );
                                num_cuts += 1;
                            }
                        }
                        continue;
                    }
                    ret = -88;
                    goto exit_function; /* unexpected */
                case 3:
                    if (da[i].typ[0] == da[i].typ[1] &&
                         da[i].typ[0] == da[i].typ[2] &&
                         da[i].typ[0] == DERIV_AMINE_tN)
                    {
                        int x, y, z;
                        x = ( da[i].num[0] < da[i].num[1] ) ? 0 : 1;
                        x = ( da[i].num[x] < da[i].num[2] ) ? x : 2; /* min */
                        z = ( da[i].num[0] < da[i].num[1] ) ? 1 : 0;
                        z = ( da[i].num[x] < da[i].num[2] ) ? 2 : z; /* max */
                        y = ( ( x + 1 ) ^ ( z + 1 ) ) - 1;                      /* median */
                        if (da[i].num[x] == da[i].num[z])
                            continue; /* all deriv. agents have same size */
                                      /* two smallest */
                        if (da[i].num[x] == da[i].num[y] && x < y)
                        {
                            int t = x; /* first cut x > y */
                            x = y;
                            y = t;
                        }
                        make_single_cut( at, da, i, x );
                        num_cuts += 1;
                        if (da[i].num[x] == da[i].num[y])
                        {
                            /* equally small */
                            make_single_cut( at, da, i, y );
                            num_cuts += 1;
                        }
                        continue;
                    }
                    ret = -88;
                    goto exit_function; /* unexpected */
                case 4:
                    if (( da[i].typ[0] & DERIV_RING_OUTSIDE_PRECURSOR ) && ( da[i].typ[1] & DERIV_RING_OUTSIDE_PRECURSOR ) &&
                        ( da[i].typ[2] & DERIV_RING_OUTSIDE_PRECURSOR ) && ( da[i].typ[3] & DERIV_RING_OUTSIDE_PRECURSOR ))
                    {
                        int n01 = inchi_max( da[i].num[0], da[i].num[1] );
                        int n23 = inchi_max( da[i].num[2], da[i].num[3] );
                        if (n01 < n23)
                        {
                            make_single_cut( at, da, i, 1 );
                            make_single_cut( at, da, i, 0 );
                            num_cuts += 1;
                        }
                        else
                        {
                            if (n01 > n23)
                            {
                                make_single_cut( at, da, i, 3 );
                                make_single_cut( at, da, i, 2 );
                                num_cuts += 1;
                            }
                        }
                        continue;
                    }
            }
        }

        nTotNumCuts += num_cuts;

#ifdef FIX_UNDERIV_TO_SDF
        if (bOutputSdf && at2)
        {
            /* replace arom bonds with original */
            replace_arom_bonds( at, num_atoms, at2, num_atoms );
            if (at2)
            {
                inchi_free( at2 );
                at2 = NULL;
            }
        }
#endif

#ifdef UNDERIV_ADD_EXPLICIT_H
        /**********  Add explicit hydrogens ************************/
        num_atoms = add_explicit_H( inp_cur_data + i_component );
#endif

        if (num_cuts)
        {
            remove_cut_derivs( num_atoms, at, inp_cur_data, i_component, &ret );
        }


    } /* for (i_component = 0; i_component < num_components; i_component++) */


    if (nTotNumCuts)
    {
        OAD_Edit_MergeComponentsAndRecreateOAD( orig_inp_data, inp_cur_data, num_components, &ret );
    }

exit_function:

    free_underiv_temp_data( ap, da, at2, inp_cur_data, num_components );

#if( UNDERIVATIZE_REPORT == 1 )
    if (!ret && nTotNumCuts && pSdfValue && bOutputReport)
    {
        numUnderiv = sort_merge_underiv( pSdfValue, bOutputSdf, szUnderivList, cDerivSeparator, underivPrefix, underivPostfix ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
        numUnderiv2 = sort_merge_underiv( pSdfValue, bOutputSdf, szUnderivList2, cDerivSeparator, underivPrefix2, underivPostfix2 ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
        sprintf(szbitUnderivList, "0x%.8X", bitUnderivList); /* djb-rwth: addressing GCC warning about 0 flag being ignored */
        numUnderiv3 = sort_merge_underiv( pSdfValue, bOutputSdf, szbitUnderivList, cDerivSeparator, underivPrefix3, underivPostfix3 ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    }
#endif

    return ret ? ret : nTotNumCuts;
}


#endif /* UNDERIVATIZE */


/****************************************************************************/
#if ( RING2CHAIN == 1 )
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


/****************************************************************************/
int detect_r2c_Zatom( inp_ATOM *at, R2C_AT *da, int iZ )
{
    int i, j, neigh, neighneigh, nRingSystem, num_found;
    R2C_AT da1;

    if (at[iZ].valence > 4)
    {
        return 0;
    }
    if (at[iZ].valence != at[iZ].chem_bonds_valence)
    {
        return 0; /* approach limitation: no double bonds */
    }
    if (at[iZ].el_number != EL_NUMBER_C)
    {
        return 0; /* sugar-specific */
    }
    if (at[iZ].nNumAtInRingSystem < 5)
    {
        return 0; /* not in a suitable ring */
    }
    if (!at[iZ].bCutVertex)
    {
        return 0;  /* recognize only type 1 for now */
    }

    nRingSystem = at[iZ].nRingSystem;
    memset( &da1, R2C_EMPTY, sizeof( da1 ) ); /* djb-rwth: memset_s C11/Annex K variant? */

    for (i = 0, num_found = 0; i < at[iZ].valence; i++)
    {
        neigh = at[iZ].neighbor[i];
        if (at[neigh].charge || at[neigh].radical)
        {
            return 0;
        }
        if (at[neigh].el_number == EL_NUMBER_O &&
             at[neigh].valence == 1 &&
             at[neigh].chem_bonds_valence == 1 &&
             at[neigh].num_H == 1)
        {
            /* found Z-OH, i.e. Z-YH */
            if (da1.ordY == R2C_EMPTY)
            {
                da1.ordY = i;
                num_found++;
                continue;
            }
            else
            {
                return 0;
            }
        }
        if (at[neigh].el_number == EL_NUMBER_O &&
             at[neigh].valence == 2 &&
             at[neigh].chem_bonds_valence == 2 &&
             at[neigh].num_H == 0 &&
             at[neigh].nRingSystem == nRingSystem)
        {
            /* found Z-O-, i.e. Z-W- */
            if (da1.ordW == R2C_EMPTY)
            {
                /* j = index of the oppozite to at[iZ] neighbor of at[neigh] */
                j = ( at[neigh].neighbor[0] == iZ );
                neighneigh = at[neigh].neighbor[j];
                if (at[neighneigh].valence != at[neighneigh].chem_bonds_valence ||
                     at[neighneigh].el_number != EL_NUMBER_C)
                    return 0; /* sugar-specific */
                da1.ordW = i;
                num_found++;
                continue;
            }
            else
            {
                return 0;
            }
        }
        if (at[neigh].el_number == EL_NUMBER_C &&
             at[neigh].valence > 2 &&
             at[neigh].chem_bonds_valence == at[neigh].valence &&
             at[neigh].num_H <= 1 &&
             at[neigh].nRingSystem == nRingSystem)
        {
            /* sugar-specfic: carbon in the ring should have -OH neighbor */
            int iOH;
            for (j = 0; j < at[neigh].valence; j++)
            {
                iOH = at[neigh].neighbor[j];
                if (at[iOH].el_number == EL_NUMBER_O &&
                     at[iOH].valence == 1 &&
                     at[iOH].chem_bonds_valence == 1 &&
                     at[iOH].num_H == 1 &&
                     !at[iOH].charge && !at[iOH].radical)
                {
                    if (da1.ordC == R2C_EMPTY)
                    {
                        da1.ordC = i;
                        num_found++;
                        break;
                    }
                    else
                    {
                        return 0;
                    }
                }
            }
            if (j < at[neigh].valence)
            {
                continue;
            }
        }
        if (at[neigh].el_number == EL_NUMBER_C &&
             at[neigh].chem_bonds_valence == at[neigh].valence &&
             at[neigh].nRingSystem != nRingSystem)
        {
            /* extra carbon neighbor of Z */
            if (da1.ordCopt == R2C_EMPTY)
            {
                da1.ordCopt = i;
                continue;
            }
        }
        return 0; /* unexpectd neighbor */
    }

    if (num_found == 3)
    {
        da1.type = 1;
        da[iZ] = da1;
        return 1; /* disconnection found */
    }

    return 0;
}


/****************************************************************************/
int cut_ring_to_chain( inp_ATOM *at, R2C_AT *da, int iZ )
{
    int ret = -1; /* error flag */
    int iordW = (int) da[iZ].ordW; /* ord of the bond in iZ */
    int iordY = (int) da[iZ].ordY; /* ord of the bond in iZ */
    int iordC = (int) da[iZ].ordC;
    int iW, iY, num_iso_H, i, jordZ;
    AT_NUMB *p;

    if (da[iZ].type != 1)
    {
        return 0;
    }
    if (0 > iordW || iordW >= at[iZ].valence ||
         0 > iordY || iordY >= at[iZ].valence ||
         0 > iordC || iordC >= at[iZ].valence /* suger-specific*/)
    {
        return -1; /* program error */
    }
    /* find other da[] that affect at[iZ] */
    iW = at[iZ].neighbor[iordW];  /* opposite atom to disconnect and add H */
    iY = at[iZ].neighbor[iordY];  /* opposite atom to increment the bond and remove H*/
    if (!at[iY].num_H || at[iZ].bond_type[iordY] != BOND_TYPE_SINGLE)
    {
        return -2; /* program error */
    }
    /* increment at[iZ]--at[iY] bond */
    p = is_in_the_list( at[iY].neighbor, (AT_NUMB) iZ, at[iY].valence );
    if (!p)
    {
        return -3; /* program error */
    }
    jordZ = p - at[iY].neighbor;
    at[iZ].bond_type[iordY] ++;
    at[iZ].chem_bonds_valence++;
    at[iY].bond_type[jordZ] ++;
    at[iY].chem_bonds_valence++;

    /* disconnect at[iZ]--at[iW] bond */
    ret = DisconnectInpAtBond( at, NULL, iZ, iordW );
    if (ret != 1)
    {
        return -4; /* program error */
    }
    /* disconnection succeeded */
    /* transfer H from at[iY] to at[iW] */
    num_iso_H = NUM_ISO_H( at, iY );
    if (at[iY].num_H == num_iso_H)
    {
        for (i = 0; i < NUM_H_ISOTOPES; i++)
        {
            if (at[iY].num_iso_H[i])
            {
                at[iY].num_iso_H[i] --;
                at[iW].num_iso_H[i] ++;
            }
        }
    }
    at[iY].num_H--;
    at[iW].num_H++;

    return 1;
}


/****************************************************************************/
int Ring2Chain( struct tagINCHI_CLOCK *ic,
                struct tagCANON_GLOBALS *pCG,
                ORIG_ATOM_DATA *orig_inp_data )
{
    int ret = 0, i, j, n, num_atoms, num_components, num, num_cuts, iZ; /* djb-rwth: removing redundant variables */
    inp_ATOM *at = orig_inp_data->at;
    INP_ATOM_DATA *inp_cur_data = NULL;
    R2C_AT        *da = NULL;


    /* prepare */

    /*set_R2C_el_numbers( );*/

    num_atoms = remove_terminal_HDT( orig_inp_data->num_inp_atoms, at, 1 );
    /*^^^^^ always accomodate accomodate FIX_TERM_H_CHRG_BUG - IPl, July 2008*/
    orig_inp_data->num_inp_atoms = num_atoms;

    /* initialize */
    UnMarkDisconnectedComponents( orig_inp_data );
    num_cuts = 0;
    /* mark */
    num_components = MarkDisconnectedComponents( orig_inp_data, 0 );
    inp_cur_data = (INP_ATOM_DATA *) inchi_calloc( num_components, sizeof( inp_cur_data[0] ) );
    iZ = -1;
    for (j = 0; j < num_components; j++)
    {
        CreateInpAtomData( inp_cur_data + j, orig_inp_data->nCurAtLen[j], 0 );
        inp_cur_data[j].num_at = ExtractConnectedComponent( orig_inp_data->at, orig_inp_data->num_inp_atoms, j + 1, inp_cur_data[j].at );
        /*  error processing */
        if (inp_cur_data[j].num_at <= 0 || orig_inp_data->nCurAtLen[j] != inp_cur_data[j].num_at)
        {
            ret = -( j + 1 ); /* severe error */
            goto exit_function;
        }
        /* initialize */
        num_atoms = inp_cur_data[j].num_at;
        at = inp_cur_data[j].at;
        add_DT_to_num_H( num_atoms, at );

        UnMarkRingSystemsInp( at, num_atoms );
        UnMarkOtherIndicators( at, num_atoms );
        UnMarkOneComponent( at, num_atoms );
        MarkRingSystemsInp( at, num_atoms, 0 );
        ret = mark_arom_bonds( ic, pCG, at, num_atoms );
        if (ret < 0)
        {
            goto exit_function;
        }
        ret = 0;
        if (da)
        {
            inchi_free( da );
        }
        da = (R2C_AT *) inchi_calloc( num_atoms, sizeof( da[0] ) );

        /* detect ring-to-chain possibilities */
        /* djb-rwth: removing redundant code */
        for (i = 0, num = 0; i < num_atoms; i++)
        {
            if (at[i].bCutVertex /* type 1 specific*/ && !da[i].type)
            {
                num += ( n = detect_r2c_Zatom( at, da, i ) );
                if (n == 1)
                {
                    iZ = i;
                }
                UnMarkOtherIndicators( at, num_atoms );
            }
        }

        if (num == 1)
        {
            /* convert ring to chain: make single cut */
            ret = cut_ring_to_chain( at, da, iZ );
            if (ret < 0)
            {
                goto exit_function;
            }
            num_cuts += ( ret == 1 );
        }
        else
        {
            if (num)
            {
                /* allocate an array of bonds to be cut */
                R2C_ATPAIR *ap = (R2C_ATPAIR *) inchi_malloc( sizeof( ap[0] ) * num );
                AT_NUMB    comp_num = 0;
                if (!ap)
                {
                    ret = -1; /* malloc failure */
                    goto exit_function;
                }
                /* fill out the array of bonds to be cut */
                for (i = j = 0; i < num_atoms; i++)
                {
                    if (da[i].type == 1)
                    {
                        AT_NUMB at1 = i;
                        AT_NUMB at2 = at[i].neighbor[(int) da[i].ordW];
                        if (j >= num)
                        {
                            ret = -2;
                            goto exit_r2c_num; /* wrong number of cuts = num */
                        }
                        n = ( at1 > at2 );
                        ap[j].at[n] = at1;
                        ap[j].at[1 - n] = at2; /* ap[j].at[0] < ap[j].at[1] */
                        j++;
                    }
                }
                if (j != num)
                {
                    ret = -3;
                    goto exit_r2c_num; /* wrong number of cuts = num */
                }
                /* sort the bonds for subsequent searching by bisections */
                qsort( ap, num, sizeof( ap[0] ), cmp_r2c_atpair );
                /* mark components to be disconnected */
                for (i = 0; i < num; i++)
                {
                    for (j = 0; j < 2; j++)
                    {
                        if (!at[ap[i].at[j]].at_type)
                        {
                            comp_num++;
                            mark_atoms_ap( at, (int) ap[i].at[j], ap, num, 0, comp_num );
                        }
                    }
                }
                /* convert ring to chain */
                for (i = 0; i < num; i++)
                {
                    int i1 = ap[i].at[0];
                    int i2 = ap[i].at[1];
                    iZ = -1;
                    if (at[i1].at_type == at[i2].at_type)
                    {
                        /* by definition, there are no adjacent iZ atoms; one iZ atom per bond */
                        if (da[i1].type == 1 && at[i1].neighbor[(int) da[i1].ordW] == i2)
                        {
                            iZ = i1;
                        }
                        else
                            if (da[i2].type == 1 && at[i2].neighbor[(int) da[i2].ordW] == i1)
                            {
                                iZ = i2;
                            }
                            else
                            {
                                ret = -3;
                                goto exit_r2c_num;
                            }
                        ret = cut_ring_to_chain( at, da, iZ );
                        if (ret < 0)
                        {
                            goto exit_r2c_num;
                        }
                        num_cuts += ( ret == 1 );
                    }
                }
                ret = 0;
            exit_r2c_num:
                inchi_free( ap );
                UnMarkOtherIndicators( at, num_atoms );
                if (ret < 0)
                {
                    goto exit_function;
                }
            }
        }
    }

    if (num_cuts)
    {
        OAD_Edit_MergeComponentsAndRecreateOAD( orig_inp_data, inp_cur_data, num_components, &ret );
    }

exit_function:

    if (da)
    {
        inchi_free( da );
        da = NULL;
    }
    for (j = 0; j < num_components; j++)
    {
        FreeInpAtomData( inp_cur_data + j );
    }
    inchi_free( inp_cur_data );
    inp_cur_data = NULL;

    return ret ? ret : num_cuts;
}


#endif /* RING2CHAIN */


#if ( RING2CHAIN == 1 || UNDERIVATIZE == 1 )
/****************************************************************************/
void OAD_Edit_MergeComponentsAndRecreateOAD( ORIG_ATOM_DATA *orig_OrigAtomData,
                                             INP_ATOM_DATA *curr_InpAtomData,
                                             int num_components,
                                             int *errcode )
{
    int i, num_atoms = 0, cur_num_at = 0;
    inp_ATOM *at;

    if (num_components <= 0)
    {
        *errcode = -999; /* num atoms mismatch */
        return;
    }

    /* Merge kept components into 'at' */
    for (i = 0; i < num_components; i++)
    {
        num_atoms += curr_InpAtomData[i].num_at;
    }

    at = (inp_ATOM *) inchi_calloc( num_atoms, sizeof( at[0] ) );
    cur_num_at = 0;

    for (i = 0; i < num_components; i++)
    {
        /* Clean and prepare */
        UnMarkRingSystemsInp( curr_InpAtomData[i].at, curr_InpAtomData[i].num_at );
        UnMarkOtherIndicators( curr_InpAtomData[i].at, curr_InpAtomData[i].num_at );
        UnMarkOneComponent( curr_InpAtomData[i].at, curr_InpAtomData[i].num_at );

        subtract_DT_from_num_H( curr_InpAtomData[i].num_at, curr_InpAtomData[i].at );

        /* Merge one by one */
        cur_num_at = add_inp_ATOM( at, num_atoms, cur_num_at, curr_InpAtomData[i].at, curr_InpAtomData[i].num_at );
    }

    /* Replace original OrigAtomData structure */
    if (cur_num_at == num_atoms)
    {
        inchi_free( orig_OrigAtomData->at );
        orig_OrigAtomData->at = at;

        orig_OrigAtomData->num_inp_atoms = cur_num_at;

        /* Destroy original coordinates as we destroyed a part of original input structure */
        if (orig_OrigAtomData->szCoord)
        {
            inchi_free( orig_OrigAtomData->szCoord );
            orig_OrigAtomData->szCoord = NULL;
        }

        UnMarkDisconnectedComponents( orig_OrigAtomData );
    }
    else
    {
        /* Create copy error! */
        if (at)
        {
            inchi_free( at );
            at = NULL;
        }
        *errcode = -999; /* num atoms mismatch */
    }

    return;
}


/****************************************************************************/
void free_underiv_temp_data( R2C_ATPAIR *ap,
                             DERIV_AT *da,
                             inp_ATOM *at2,
                             INP_ATOM_DATA *inp_cur_data,
                             int num_components )
{
    int i_component;
    if (ap)
    {
        inchi_free( ap );
        ap = NULL;
    }
    if (da)
    {
        inchi_free( da );
        da = NULL;
    }
#ifdef FIX_UNDERIV_TO_SDF
    if (at2)
    {
        inchi_free( at2 );
        at2 = NULL;
    }
#endif
    for (i_component = 0; i_component < num_components; i_component++)
    {
        FreeInpAtomData( inp_cur_data + i_component );
    }
    inchi_free( inp_cur_data );
    inp_cur_data = NULL;
}


/****************************************************************************/
void remove_cut_derivs( int num_atoms,
                        inp_ATOM *at,
                        INP_ATOM_DATA *inp_cur_data,
                        int i_component,
                        int *errcode )
{
    /* Remove marked with Tritium disconnected derivative attachments */
    ORIG_ATOM_DATA Orig_inp_data1, *orig_inp_data1;
    INP_ATOM_DATA *inp_cur_data1 = NULL;
    int i, num_components1, i_component1, cur_num_at = 0; /* djb-rwth: removing redundant variables */

    orig_inp_data1 = &Orig_inp_data1;
    memset( orig_inp_data1, 0, sizeof( orig_inp_data1[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */

    UnMarkRingSystemsInp( at, num_atoms );
    UnMarkOtherIndicators( at, num_atoms );
    UnMarkOneComponent( at, num_atoms );

    for (i = 0; i < num_atoms; i++)
    {
        orig_inp_data1->num_inp_bonds += at[i].valence;
    }
    orig_inp_data1->num_inp_bonds /= 2;
    orig_inp_data1->num_inp_atoms = num_atoms;

    orig_inp_data1->at = at; /* = from inp_cur_data[i_component].at */

    num_components1 = MarkDisconnectedComponents( orig_inp_data1, 0 );

    inp_cur_data1 = (INP_ATOM_DATA *) inchi_calloc( num_components1, sizeof( inp_cur_data1[0] ) );

    if (inp_cur_data1 && (num_components1 > 0)) /* djb-rwth: fixing a NULL pointer dereference */
    {
    /* Extract components and discard disconnected derivatizing agents */
    for (i_component1 = 0; i_component1 < num_components1; i_component1++)
    {
        CreateInpAtomData( inp_cur_data1 + i_component1, orig_inp_data1->nCurAtLen[i_component1], 0 );
        inp_cur_data1[i_component1].num_at = ExtractConnectedComponent( orig_inp_data1->at, orig_inp_data1->num_inp_atoms,
                                                                        i_component1 + 1, inp_cur_data1[i_component1].at );
        /*  error processing */
        if (inp_cur_data1[i_component1].num_at <= 0 || orig_inp_data1->nCurAtLen[i_component1] != inp_cur_data1[i_component1].num_at)
        {
            *errcode = -( i_component1 + 1 ); /* severe error */
            break;
        }
        /* if the component has tritium then discard: it is a derivatizing agent */
        for (i = 0; i < inp_cur_data1[i_component1].num_at; i++)
        {
#ifdef UNDERIV_ADD_D_TO_PRECURSOR
            if (inp_cur_data1[i_component1].at[i].num_iso_H[1])
            {
                inp_cur_data1[i_component1].at[i].num_iso_H[1] = 0; /* remove deuterium */
            }
#endif
            if (inp_cur_data1[i_component1].at[i].num_iso_H[2])
            {
                FreeInpAtomData( inp_cur_data1 + i_component1 );
                break;
            }
        }
    }
    /* Merge components into one -- must be only one */
    for (i_component1 = 0, num_atoms = 0; i_component1 < num_components1; i_component1++)
    {
        num_atoms += inp_cur_data1[i_component1].num_at;
    }
    at = (inp_ATOM *) inchi_calloc( num_atoms, sizeof( at[0] ) );
    cur_num_at = 0;
    for (i_component1 = 0; i_component1 < num_components1; i_component1++)
    {
        /* clean and prepare */
        if (!inp_cur_data1[i_component1].num_at)
        {
            continue; /* removed derivatizing object */
                      /*UnMarkOneComponent( inp_cur_data1[i_component1].at, inp_cur_data1[i_component1].num_at );*/
                      /* merge one by one */
        }
        cur_num_at = add_inp_ATOM( at, num_atoms, cur_num_at, inp_cur_data1[i_component1].at, inp_cur_data1[i_component1].num_at );
        FreeInpAtomData( inp_cur_data1 + i_component1 ); /* cleanup */
        /* djb-rwth: removing redundant code */
    }
    }

    /* Replace the component */
    /* Order of the following two statements is critically important */

    UnMarkDisconnectedComponents( orig_inp_data1 ); /* orig_inp_data1->at is same as inp_cur_data[i_component].at */
    FreeInpAtomData( inp_cur_data + i_component ); /* cleanup the original component */

    inp_cur_data[i_component].at = at;
    inp_cur_data[i_component].num_at = cur_num_at;
    inchi_free( inp_cur_data1 );
}
#endif /* #if ( RING2CHAIN == 1 || UNDERIVATIZE == 1 ) */


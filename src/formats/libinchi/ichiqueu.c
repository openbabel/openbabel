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

#include <string.h>

#include "mode.h"
#include "ichitaut.h"

#include "bcf_s.h"

/****************************************************************************/

#if ( FIND_RING_SYSTEMS == 1 ) /* { */

/* local prototypes */
int are_alt_bonds( U_CHAR *bonds, int len );
int AddBondsPos( inp_ATOM *atom, T_BONDPOS *BondPosTmp, int nNumBondPosTmp, T_BONDPOS *BondPos, int nMaxNumBondPos, int nNumBondPos );
int AddEndPoints( T_ENDPOINT *EndPointTmp, int nNumNewEndPoint, T_ENDPOINT *EndPoint, int nMaxNumEndPoint, int nNumEndPoint );





/******************************************
 *
 *  Tautomerism in 5- and 6-member rings
 *
 ******************************************/

const int NONE = (AT_RANK) ~0;


/*
  1,5 Tautomerism in 6-member alt ring:

   /=\          /==\
 HN   C=O  <-> N    C-OH
   \=/          \\-//


   1,2 Tautomerism in 5-member ring:


  HN--X           N==X
   |   \\         |   \
   |    Z  <->    |    Z
   |   /          |   //
   N==Y          HN--Y


  1,4 tautomerism in 7-member ring

      /C==D             //C-D
   O=B     \        HO-B     \\
     |      E  <->     |      E
  HO-A     //        O=A     /
      \\G-F             \\G-F


  1,4 tautomerism in 5-member ring


   O=B--C            O-B==C
     |   \\            |   \
     |     D  <->      |     D
     |   /             |   //
  HO-A==E           HO=A--E

*/
typedef int CHECK_DFS_RING( struct tagCANON_GLOBALS *pCG,
    inp_ATOM *atom,
    DFS_PATH *DfsPath,
    int nLenDfsPath,
    int nStartAtomNeighbor,
    int nStartAtomNeighbor2,
    int nStartAtomNeighborNeighbor,
    T_ENDPOINT *EndPoint,
    int nMaxNumEndPoint,
    T_BONDPOS  *BondPos,
    int nMaxNumBondPos,
    int *pnNumEndPoint,
    int *pnNumBondPos,
    struct BalancedNetworkStructure *pBNS,
    struct BalancedNetworkData *pBD,
    int num_atoms );

typedef int CHECK_CENTERPOINT( inp_ATOM *atom, int iat );

CHECK_DFS_RING Check7MembTautRing;
CHECK_DFS_RING Check6MembTautRing;
CHECK_DFS_RING Check5MembTautRing;

#if ( TAUT_15_NON_RING == 1 )  /* post v.1 feature */
/* DFS simple alt path for 1,5 tautomerism, post v.1 feature */
typedef int CHECK_DFS_PATH( struct tagCANON_GLOBALS *pCG,
    inp_ATOM *atom,
    DFS_PATH *DfsPath,
    int nLenDfsPath,
    int jNxtNeigh,
    int nStartAtomNeighbor,
    int nStartAtomNeighbor2,
    int nStartAtomNeighborNeighbor,
    T_ENDPOINT *EndPoint,
    int nMaxNumEndPoint,
    T_BONDPOS  *BondPos,
    int nMaxNumBondPos,
    int *pnNumEndPoint,
    int *pnNumBondPos,
    struct BalancedNetworkStructure *pBNS,
    struct BalancedNetworkData *pBD,
    int num_atoms );

typedef int CHECK_DFS_CENTERPOINT( inp_ATOM *atom, DFS_PATH *DfsPath, int nLenDfsPath, int jNxtNeigh,
    struct BalancedNetworkStructure *pBNS,
    struct BalancedNetworkData *pBD, int num_atoms );


CHECK_DFS_PATH        Check15TautPath;
CHECK_DFS_CENTERPOINT Check15TautPathCenterpoint;

int DFS_FindTautAltPath( struct tagCANON_GLOBALS *pCG,
    inp_ATOM *atom,
    int nStartAtom,
    int nStartAtomNeighbor,
    int nStartAtomNeighbor2,
    int nStartAtomNeighborNeighbor,
    int nCycleLen,
    AT_RANK  *nDfsPathPos,
    DFS_PATH *DfsPath,
    CHECK_DFS_PATH *CheckDfsPath,
    CHECK_DFS_CENTERPOINT *CheckCenterPoint,
    T_ENDPOINT *EndPoint,
    int nMaxNumEndPoint,
    T_BONDPOS  *BondPos,
    int nMaxNumBondPos,
    int *pnNumEndPoint,
    int *pnNumBondPos,
    struct BalancedNetworkStructure *pBNS,
    struct BalancedNetworkData *pBD,
    int num_atoms );

#define BOND_WRONG 64
#define IS_ALT_OR_DBLBOND(X) (((X) == BOND_SINGLE || (X) == BOND_DOUBLE)? (X) : \
                              ((X) == BOND_ALTERN || (X) == BOND_TAUTOM || (X) == BOND_ALT12NS)? BOND_ALTERN : \
                              BOND_WRONG);

#endif /* TAUT_15_NON_RING */

int DFS_FindTautInARing( struct tagCANON_GLOBALS *pCG,
    inp_ATOM *atom,
    int nStartAtom,
    int nStartAtomNeighbor,
    int nStartAtomNeighbor2,
    int nStartAtomNeighborNeighbor,
    int nCycleLen,
    AT_RANK  *nDfsPathPos,
    DFS_PATH *DfsPath,
    CHECK_DFS_RING *CheckDfsRing,
    CHECK_CENTERPOINT *CheckCenterPoint,
    T_ENDPOINT *EndPoint,
    int nMaxNumEndPoint,
    T_BONDPOS  *BondPos,
    int nMaxNumBondPos,
    int *pnNumEndPoint,
    int *pnNumBondPos,
    struct BalancedNetworkStructure *pBNS,
    struct BalancedNetworkData *pBD,
    int num_atoms );

#if ( REPLACE_ALT_WITH_TAUT == 1 )
#define REPLACE_THE_BOND(X) ( (X) == BOND_SINGLE || (X) == BOND_DOUBLE || (X) == BOND_ALTERN || (X) == BOND_ALT12NS )
#else
#define REPLACE_THE_BOND(X) ( (X) == BOND_SINGLE || (X) == BOND_DOUBLE )
#endif


/****************************************************************************/
int bIsCenterPointStrict( inp_ATOM *atom, int iat )
{
    if (atom[iat].valence == atom[iat].chem_bonds_valence)
    {
        int endpoint_valence = get_endpoint_valence( atom[iat].el_number );
        if (endpoint_valence && ( (endpoint_valence > atom[iat].valence && /* added a check for negative charge or H 3-31-03 */
            ( atom[iat].num_H || atom[iat].charge == -1 )) ||
            (!atom[iat].charge && atom[iat].c_point) )) /* djb-rwth: addressing LLVM warnings */
        {
            return 1; /*  may appear to be tautomeric or chargable
                          (this increases chem_bonds_valence), should be explored */
        }
        return 0;
    }
    if (atom[iat].valence + 1 == atom[iat].chem_bonds_valence &&
        is_centerpoint_elem_strict( atom[iat].el_number ))
    {
        return 1;
    }

    return 0;
}


/****************************************************************************/
int nGet14TautIn7MembAltRing( struct tagCANON_GLOBALS *pCG,
                              inp_ATOM *atom,
                              int nStartAtom,
                              int nStartAtomNeighbor,
                              int nStartAtomNeighborEndpoint,
                              int nStartAtomNeighborNeighborEndpoint,
                              AT_RANK  *nDfsPathPos,
                              DFS_PATH *DfsPath,
                              int nMaxLenDfsPath,
                              T_ENDPOINT *EndPoint,
                              int nMaxNumEndPoint,
                              T_BONDPOS  *BondPos,
                              int nMaxNumBondPos,
                              int *pnNumEndPoint,
                              int *pnNumBondPos,
                              struct BalancedNetworkStructure *pBNS,
                              struct BalancedNetworkData *pBD,
                              int num_atoms )
{
    int nRet;

    *pnNumEndPoint = 0;
    *pnNumBondPos = 0;

    if (nMaxLenDfsPath <= 7)
    {
        return -1; /*  path is too short */
    }

    nRet = DFS_FindTautInARing( pCG, atom, nStartAtom, nStartAtomNeighbor,
                                nStartAtomNeighborEndpoint,
                                nStartAtomNeighborNeighborEndpoint, 7,
                                nDfsPathPos, DfsPath, Check7MembTautRing,
                                bIsCenterPointStrict, EndPoint, nMaxNumEndPoint,
                                BondPos, nMaxNumBondPos, pnNumEndPoint,
                                pnNumBondPos, pBNS, pBD, num_atoms );


    return nRet;
}


/****************************************************************************/
int nGet14TautIn5MembAltRing( struct tagCANON_GLOBALS *pCG,
                              inp_ATOM *atom,
                              int nStartAtom,
                              int nStartAtomNeighbor,
                              int nStartAtomNeighborEndpoint,
                              int nStartAtomNeighborNeighborEndpoint,
                              AT_RANK  *nDfsPathPos,
                              DFS_PATH *DfsPath,
                              int nMaxLenDfsPath,
                              T_ENDPOINT *EndPoint,
                              int nMaxNumEndPoint,
                              T_BONDPOS  *BondPos,
                              int nMaxNumBondPos,
                              int *pnNumEndPoint,
                              int *pnNumBondPos,
                              struct BalancedNetworkStructure *pBNS,
                              struct BalancedNetworkData *pBD,
                              int num_atoms )
{
    int nRet;

    *pnNumEndPoint = 0;
    *pnNumBondPos = 0;

    if (nMaxLenDfsPath <= 5)
    {
        return -1; /*  path is too short */
    }

    nRet = DFS_FindTautInARing( pCG, atom, nStartAtom, nStartAtomNeighbor,
        nStartAtomNeighborEndpoint, nStartAtomNeighborNeighborEndpoint, 5,
        nDfsPathPos, DfsPath, Check7MembTautRing, bIsCenterPointStrict,
        EndPoint, nMaxNumEndPoint, BondPos, nMaxNumBondPos,
        pnNumEndPoint, pnNumBondPos, pBNS, pBD, num_atoms );

    return nRet;
}

/****************************************************************************/
int nGet12TautIn5MembAltRing( struct tagCANON_GLOBALS *pCG,
                              inp_ATOM *atom,
                              int nStartAtom,
                              int nStartAtomNeighbor,
                              AT_RANK  *nDfsPathPos,
                              DFS_PATH *DfsPath,
                              int nMaxLenDfsPath,
                              T_ENDPOINT *EndPoint,
                              int nMaxNumEndPoint,
                              T_BONDPOS  *BondPos,
                              int nMaxNumBondPos,
                              int *pnNumEndPoint,
                              int *pnNumBondPos,
                              struct BalancedNetworkStructure *pBNS,
                              struct BalancedNetworkData *pBD,
                              int num_atoms )
{
    int nRet;

    *pnNumEndPoint = 0;
    *pnNumBondPos = 0;

    if (nMaxLenDfsPath <= 5)
    {
        return -1; /*  path is too short */
    }

    nRet = DFS_FindTautInARing( pCG, atom, nStartAtom, nStartAtomNeighbor,
                                -1, -1, 5,
                                nDfsPathPos, DfsPath, Check5MembTautRing,
                                bIsCenterPointStrict, EndPoint, nMaxNumEndPoint,
                                BondPos, nMaxNumBondPos, pnNumEndPoint,
                                pnNumBondPos, pBNS, pBD, num_atoms );

    return nRet;
}


/****************************************************************************/
int nGet15TautIn6MembAltRing( struct tagCANON_GLOBALS *pCG,
                              inp_ATOM *atom,
                              int nStartAtom,
                              AT_RANK  *nDfsPathPos,
                              DFS_PATH *DfsPath,
                              int nMaxLenDfsPath,
                              T_ENDPOINT *EndPoint,
                              int nMaxNumEndPoint,
                              T_BONDPOS  *BondPos,
                              int nMaxNumBondPos,
                              int *pnNumEndPoint,
                              int *pnNumBondPos,
                              struct BalancedNetworkStructure *pBNS,
                              struct BalancedNetworkData *pBD,
                              int num_atoms )
{
    int nRet;

    *pnNumEndPoint = 0;
    *pnNumBondPos = 0;

    if (nMaxLenDfsPath <= 7)
    {
        return -1; /*  path is too short */
    }

    nRet = DFS_FindTautInARing( pCG, atom, nStartAtom,
                                -1/*nStartAtomNeighbor*/,
                                -1/*nStartAtomNeighbor2*/,
                                -1/*nStartAtomNeighborNeighbor*/,
                                6 /* nCycleLen*/,
                                nDfsPathPos, DfsPath,
                                Check6MembTautRing, bIsCenterPointStrict,
                                EndPoint, nMaxNumEndPoint,
                                BondPos, nMaxNumBondPos,
                                pnNumEndPoint, pnNumBondPos,
                                pBNS, pBD, num_atoms );

    return nRet;
}


#if ( TAUT_15_NON_RING      == 1 ) /***** post v.1 feature *****/


/****************************************************************************/
int nGet15TautInAltPath( struct tagCANON_GLOBALS *pCG,
                         inp_ATOM *atom,
                         int nStartAtom,
                         AT_RANK  *nDfsPathPos,
                         DFS_PATH *DfsPath,
                         int nMaxLenDfsPath,
                         T_ENDPOINT *EndPoint,
                         int nMaxNumEndPoint,
                         T_BONDPOS  *BondPos,
                         int nMaxNumBondPos,
                         int *pnNumEndPoint,
                         int *pnNumBondPos,
                         struct BalancedNetworkStructure *pBNS,
                         struct BalancedNetworkData *pBD,
                         int num_atoms )
{
    int nRet;

    *pnNumEndPoint = 0;
    *pnNumBondPos = 0;

    if (nMaxLenDfsPath <= 7)
    {
        return -1; /*  path is too short */
    }

    nRet = DFS_FindTautAltPath( pCG, atom, nStartAtom,
                                -1/*nStartAtomNeighbor*/,
                                -1/*nStartAtomNeighbor2*/,
                                -1/*nStartAtomNeighborNeighbor*/,
                                4 /* nCycleLen*/,
                                nDfsPathPos, DfsPath,
                                Check15TautPath, Check15TautPathCenterpoint,
                                EndPoint, nMaxNumEndPoint,
                                BondPos, nMaxNumBondPos,
                                pnNumEndPoint, pnNumBondPos,
                                pBNS, pBD, num_atoms );

    return nRet;
}
#endif


/****************************************************************************
  DFS version
****************************************************************************/
#define MAX_DFS_DEPTH 16



/********************************************************************************/
int DFS_FindTautInARing( struct tagCANON_GLOBALS *pCG,
                         inp_ATOM *atom,
                         int nStartAtom,
                         int nStartAtomNeighbor,
                         int nStartAtomNeighbor2,
                         int nStartAtomNeighborNeighbor,
                         int nCycleLen,
                         AT_RANK  *nDfsPathPos,
                         DFS_PATH *DfsPath,
                         CHECK_DFS_RING *CheckDfsRing,
                         CHECK_CENTERPOINT *CheckCenterPoint,
                         T_ENDPOINT *EndPoint,
                         int nMaxNumEndPoint,
                         T_BONDPOS  *BondPos,
                         int nMaxNumBondPos,
                         int *pnNumEndPoint,
                         int *pnNumBondPos,
                         struct BalancedNetworkStructure *pBNS,
                         struct BalancedNetworkData *pBD,
                         int num_atoms )
{
    /*  Depth First Search */
    /*  Ignore all atoms not belonging to the current ring system (=biconnected component) */
    AT_RANK      nMinLenDfsPath;
    int          j, cur_at, nxt_at, prv_at;
    int          nLenDfsPath, nNumFound, ret;
    /* djb-rwth: removing redundant variables */
    int          nDoNotTouchAtom1 = -1, nDoNotTouchAtom2 = -1;

    nLenDfsPath = 0;
    nNumFound = 0;

    nCycleLen--;

    DfsPath[nLenDfsPath].at_no = cur_at = nStartAtom;
    DfsPath[nLenDfsPath].bond_type = 0;
    DfsPath[nLenDfsPath].bond_pos = -1;
    nDfsPathPos[cur_at] = nLenDfsPath + 1;  /*  mark */
    /* djb-rwth: removing redundant code */
    nMinLenDfsPath = 0;
    if (nStartAtomNeighbor2 >= 0)
    {
        nDoNotTouchAtom1 = (int) atom[cur_at].neighbor[nStartAtomNeighbor2];
    }

    /*  add the first neighbor to the 2nd tree position if required */
    if (nStartAtomNeighbor >= 0)
    {
        j = nStartAtomNeighbor;
        prv_at = cur_at;
        cur_at = atom[prv_at].neighbor[j];
        DfsPath[nLenDfsPath].bond_type = ( atom[prv_at].bond_type[j] & ~BOND_MARK_ALL );
#if ( FIX_BOND23_IN_TAUT == 1 )
        DfsPath[nLenDfsPath].bond_type = ACTUAL_ORDER( pBNS, prv_at, j, DfsPath[nLenDfsPath].bond_type );
#endif
        DfsPath[nLenDfsPath].bond_pos = j;

        nLenDfsPath++;

        DfsPath[nLenDfsPath].at_no = cur_at;
        DfsPath[nLenDfsPath].bond_type = 0;
        DfsPath[nLenDfsPath].bond_pos = -1;
        nDfsPathPos[cur_at] = nLenDfsPath + 1;
        nMinLenDfsPath++;
        if (nStartAtomNeighborNeighbor >= 0)
        {
            nDoNotTouchAtom2 = (int) atom[cur_at].neighbor[nStartAtomNeighborNeighbor];
        }
    }

    /*  MAIN DFS CYCLE: may find one and the same t-group 2 times; saves only one instance */
    /*  traverse *all* paths starting at atom[nStartAtom]; max. path length = (nCycleLen+1)  */
    while (nLenDfsPath >= nMinLenDfsPath)
    {
        j = ++DfsPath[nLenDfsPath].bond_pos;
        if (j < atom[cur_at = (int) DfsPath[nLenDfsPath].at_no].valence)
        {
            DfsPath[nLenDfsPath].bond_type = ( atom[cur_at].bond_type[j] & ~BOND_MARK_ALL );
#if ( FIX_BOND23_IN_TAUT == 1 )
            DfsPath[nLenDfsPath].bond_type = ACTUAL_ORDER( pBNS, cur_at, j, DfsPath[nLenDfsPath].bond_type );
#endif
            nxt_at = (int) atom[cur_at].neighbor[j];
            if (nxt_at == nDoNotTouchAtom1 ||
                nxt_at == nDoNotTouchAtom2)
            {
                ; /*  ignore */
            }
            else
            {
                if (nDfsPathPos[nxt_at])
                {
                    /*  found a ring closure or a step backwards */
                    if (1 == nDfsPathPos[nxt_at] && nLenDfsPath == nCycleLen)
                    {
                        /*  we have found the cycle; check it */
                        ret = ( *CheckDfsRing )( pCG,
                            atom,
                            DfsPath, nLenDfsPath,
                            nStartAtomNeighbor,
                            nStartAtomNeighbor2,
                            nStartAtomNeighborNeighbor,
                            EndPoint, nMaxNumEndPoint,
                            BondPos, nMaxNumBondPos,
                            pnNumEndPoint, pnNumBondPos,
                            pBNS, pBD, num_atoms );
                        if (ret < 0)
                        {
                            nNumFound = ret;
                            goto clear_path;
                        }
                        nNumFound += ret;
                    }
                }
                else
                {
                    if (!( *CheckCenterPoint )( atom, nxt_at ))
                    {
                        ; /*  cannot advance to a non-centerpoint; ignore */
                    }
                    else
                    {
                        if (nLenDfsPath < nCycleLen)
                        {
                            /*  advance */
                            nLenDfsPath++;
                            cur_at = nxt_at;
                            DfsPath[nLenDfsPath].at_no = cur_at;
                            DfsPath[nLenDfsPath].bond_type = 0;
                            DfsPath[nLenDfsPath].bond_pos = -1;
                            nDfsPathPos[cur_at] = nLenDfsPath + 1;  /*  mark */
                        }
                    }
                }
            }
        }
        else
        {
            /*  retract */
            nDfsPathPos[(int) DfsPath[nLenDfsPath].at_no] = 0;
            nLenDfsPath--;
        }
    }

clear_path:
    while (0 <= nLenDfsPath)
    {
        nDfsPathPos[(int) DfsPath[nLenDfsPath].at_no] = 0;
        nLenDfsPath--;
    }

    return nNumFound;
}


#if ( TAUT_15_NON_RING      == 1 ) /***** post v.1 feature *****/


/****************************************************************************/
int DFS_FindTautAltPath( struct tagCANON_GLOBALS *pCG,
                         inp_ATOM *atom,
                         int nStartAtom,
                         int nStartAtomNeighbor,
                         int nStartAtomNeighbor2,
                         int nStartAtomNeighborNeighbor,
                         int nCycleLen,
                         AT_RANK  *nDfsPathPos,
                         DFS_PATH *DfsPath,
                         CHECK_DFS_PATH *CheckDfsPath,
                         CHECK_DFS_CENTERPOINT *CheckCenterPointPath,
                         T_ENDPOINT *EndPoint,
                         int nMaxNumEndPoint,
                         T_BONDPOS  *BondPos,
                         int nMaxNumBondPos,
                         int *pnNumEndPoint,
                         int *pnNumBondPos,
                         struct BalancedNetworkStructure *pBNS,
                         struct BalancedNetworkData *pBD,
                         int num_atoms )
{
    /*  Naive Depth First Search: same atom may be approached along different alt paths */
    /*  Ignore all atoms not belonging to the current ring system (=biconnected component) */
    AT_RANK      nMinLenDfsPath;
    int          j, cur_at, nxt_at, prv_at;
    int          nLenDfsPath, nNumFound, ret;
    /* djb-rwth: removing redundant variables */
    int          nDoNotTouchAtom1 = -1, nDoNotTouchAtom2 = -1;

    nLenDfsPath = 0;
    nNumFound = 0;

    nCycleLen--; /* indef of the last atom in the alt path, statring from 0 */

    DfsPath[nLenDfsPath].at_no = cur_at = nStartAtom;
    DfsPath[nLenDfsPath].bond_type = 0;
    DfsPath[nLenDfsPath].bond_pos = -1;  /* initialize index of the bond to the next atom */
    nDfsPathPos[cur_at] = nLenDfsPath + 1;  /*  mark with distance + 1 */
    /* djb-rwth: removing redundant variables/code */
    nMinLenDfsPath = 0;  /* allow to restart from nStartAtom */
    if (nStartAtomNeighbor2 >= 0)
    {
        nDoNotTouchAtom1 = (int) atom[cur_at].neighbor[nStartAtomNeighbor2];
    }

    /*  add the first neighbor to the 2nd tree position if required */
    if (nStartAtomNeighbor >= 0)
    {
        j = nStartAtomNeighbor;
        prv_at = cur_at;
        cur_at = atom[prv_at].neighbor[j];
        DfsPath[nLenDfsPath].bond_type = ( atom[prv_at].bond_type[j] & ~BOND_MARK_ALL );
#if ( FIX_BOND23_IN_TAUT == 1 )
        DfsPath[nLenDfsPath].bond_type = ACTUAL_ORDER( pBNS, prv_at, j, DfsPath[nLenDfsPath].bond_type );
#endif
        DfsPath[nLenDfsPath].bond_pos = j; /* fix index of the bond to the next atom */

        nLenDfsPath++;

        DfsPath[nLenDfsPath].at_no = cur_at;
        DfsPath[nLenDfsPath].bond_type = 0;
        DfsPath[nLenDfsPath].bond_pos = -1;
        nDfsPathPos[cur_at] = nLenDfsPath + 1; /* mark with distance + 1 */
        nMinLenDfsPath++;                 /* allow to restart from nStartAtom's neighbor */
        if (nStartAtomNeighborNeighbor >= 0)
        {
            nDoNotTouchAtom2 = (int) atom[cur_at].neighbor[nStartAtomNeighborNeighbor];
        }
    }

    /*  MAIN DFS CYCLE: may find one and the same t-group 2 times; saves only one instance */
    /*  traverse *all* paths starting at atom[nStartAtom]; max. path length = (nCycleLen+1)  */
    while (nLenDfsPath >= nMinLenDfsPath)
    {
        j = ++DfsPath[nLenDfsPath].bond_pos;
        if (j < atom[cur_at = (int) DfsPath[nLenDfsPath].at_no].valence)
        {
            DfsPath[nLenDfsPath].bond_type = ( atom[cur_at].bond_type[j] & ~BOND_MARK_ALL );
#if ( FIX_BOND23_IN_TAUT == 1 )
            DfsPath[nLenDfsPath].bond_type = ACTUAL_ORDER( pBNS, cur_at, j, DfsPath[nLenDfsPath].bond_type );
#endif
            nxt_at = (int) atom[cur_at].neighbor[j];
            if (nxt_at == nDoNotTouchAtom1 || /* forbidden */
                nxt_at == nDoNotTouchAtom2 || /* forbidden */
                nDfsPathPos[nxt_at] || /* ring closure */
                (nLenDfsPath && nxt_at == (int) DfsPath[nLenDfsPath - 1].at_no) /* step backwards */
                ) /* djb-rwth: addressing LLVM warning */
            {
                ; /* ignore nxt_at */
            }
            else
            {
                if (nLenDfsPath == nCycleLen &&
                       /* 1,5 and at least one of the endpoints is not in a ring */
                    ( atom[nxt_at].nNumAtInRingSystem == 1 || atom[nStartAtom].nNumAtInRingSystem == 1 ) &&
                       /*  we have found the alt path of the requested length; check it */
                       /* calling Check15TautPath() */
                       ( ret = ( *CheckDfsPath )( pCG,
                                                  atom,
                                                  DfsPath, nLenDfsPath,
                                                  j, nStartAtomNeighbor,
                                                  nStartAtomNeighbor2,
                                                  nStartAtomNeighborNeighbor,
                                                  EndPoint, nMaxNumEndPoint, BondPos, nMaxNumBondPos,
                                                  pnNumEndPoint, pnNumBondPos,
                                                  pBNS, pBD, num_atoms ) ))
                {
                    if (ret < 0)
                    {
                        nNumFound = ret;
                        goto clear_path; /* program error */
                    }
                    nNumFound += ret; /* success */
                }
                else
                {
                    /* calling Check15TautPathCenterpoint() */
                    if (!( *CheckCenterPointPath )( atom, DfsPath, nLenDfsPath, j,
                                                    pBNS, pBD, num_atoms ))
                    {
                        ; /*  cannot advance to a non-centerpoint; ignore */
                    }
                    else
                    {
                        if (nLenDfsPath < nCycleLen)
                        {
                            /*  advance */
                            nLenDfsPath++;
                            cur_at = nxt_at;
                            DfsPath[nLenDfsPath].at_no = cur_at;
                            DfsPath[nLenDfsPath].bond_type = 0;
                            DfsPath[nLenDfsPath].bond_pos = -1;
                            nDfsPathPos[cur_at] = nLenDfsPath + 1;  /*  mark */
                        }
                    }
                }
            }
        }
        else
        {
            /*  retract */
            nDfsPathPos[(int) DfsPath[nLenDfsPath].at_no] = 0;
            nLenDfsPath--;
        }
    }

clear_path:
    while (0 <= nLenDfsPath)
    {
        nDfsPathPos[(int) DfsPath[nLenDfsPath].at_no] = 0;
        nLenDfsPath--;
    }

    return nNumFound;
}

#endif /* TAUT_15_NON_RING */


/****************************************************************************
      Check if bonds are alternating
****************************************************************************/
int are_alt_bonds( U_CHAR *bonds, int len )
{
    U_CHAR next_bond;
    int           i, bAnyBond, bTautBondPresent = BOND_ALTERN;
    if (len < 2 || bonds[0] == BOND_TRIPLE || bonds[0] == BOND_ALT_13)
    {
        return 0;
    }
    next_bond = bonds[0] == BOND_SINGLE ? BOND_DOUBLE : bonds[0] == BOND_DOUBLE ? BOND_SINGLE : 0; /* djb-rwth: removing redundant code; ignoring LLVM warning: possible presence of global variables */
    if (bonds[0] == BOND_TAUTOM)
    {
        bTautBondPresent = BOND_TAUTOM;
        next_bond = 0;
    }
    else
    {
        next_bond = bonds[0] == BOND_SINGLE ? BOND_DOUBLE : bonds[0] == BOND_DOUBLE ? BOND_SINGLE : 0;
    }

    for (i = 1; i < len; i++)
    {
        if (bonds[i] == BOND_TAUTOM)
        {
            bTautBondPresent = BOND_TAUTOM;
            bAnyBond = 1;
        }
        else
        {
            bAnyBond = ( bonds[i] == BOND_ALTERN || bonds[i] == BOND_ALT12NS );
        }
        if (next_bond)
        {
            if (bonds[i] == next_bond || bAnyBond)
            {
                next_bond = ( next_bond == BOND_SINGLE ) ? BOND_DOUBLE : BOND_SINGLE;
                continue;
            }
            return 0;
        }
        else
        {
            if (bonds[i] == BOND_SINGLE)
            {
                next_bond = BOND_DOUBLE;
                continue;
            }
            else
            {
                if (bonds[i] == BOND_DOUBLE)
                {
                    next_bond = BOND_SINGLE;
                    continue;
                }
                else
                {
                    if (!bAnyBond)
                    {
                        return 0;
                    }
                }
            }
        }
    }
    return !next_bond ? bTautBondPresent : ( next_bond == BOND_SINGLE )
                            ? BOND_DOUBLE : BOND_SINGLE; /* bond to the end atom */
}

/****************************************************************************/
int AddBondsPos( inp_ATOM *atom,
                 T_BONDPOS *BondPosTmp,
                 int nNumBondPosTmp,
                 T_BONDPOS *BondPos,
                 int nMaxNumBondPos,
                 int nNumBondPos )
{
    int i, j, k, cur_at, nxt_at;
    /*  add opposite direction bonds to BondPosTmp */
    for (j = 0; j < nNumBondPosTmp; j += 2)
    {
        cur_at = BondPosTmp[j].nAtomNumber;
        nxt_at = atom[cur_at].neighbor[(int) BondPosTmp[j].neighbor_index];
        for (k = 0; k < atom[nxt_at].valence; k++)
        {
            if (cur_at == atom[nxt_at].neighbor[k])
            {
                BondPosTmp[j + 1].nAtomNumber = nxt_at;
                BondPosTmp[j + 1].neighbor_index = k;
                break;
            }
        }
    }
    /*  add new tautomeric bonds */
    for (j = 0; j < nNumBondPosTmp; j += 2)
    {
        for (i = 0; i < nNumBondPos; i++)
        {
            if ((BondPos[i].nAtomNumber == BondPosTmp[j].nAtomNumber &&
                BondPos[i].neighbor_index == BondPosTmp[j].neighbor_index) ||
                (BondPos[i].nAtomNumber == BondPosTmp[j + 1].nAtomNumber &&
                BondPos[i].neighbor_index == BondPosTmp[j + 1].neighbor_index))  /* djb-rwth: addressing LLVM warnings */
            {
                break; /*  bond has already been added */
            }
        }
        if (i == nNumBondPos)
        {
            if (i > nMaxNumBondPos)
            {
                return -1; /*  overflow */
            }
            BondPos[nNumBondPos++] = BondPosTmp[j];
        }
    }

    return nNumBondPos;
}


/****************************************************************************/
int AddEndPoints( T_ENDPOINT *EndPointTmp,
                  int nNumNewEndPoint,
                  T_ENDPOINT *EndPoint,
                  int nMaxNumEndPoint,
                  int nNumEndPoint )
{
    int i, j;
    /*  add new endpoints */
    for (j = 0; j < nNumNewEndPoint; j++)
    {
        for (i = 0; i < nNumEndPoint; i++)
        {
            if (EndPoint[i].nAtomNumber == EndPointTmp[j].nAtomNumber)
            {
                break;
            }
        }
        if (i == nNumEndPoint)
        {
            if (i > nMaxNumEndPoint)
            {
                return -1; /*  overflow */
            }
            EndPoint[nNumEndPoint++] = EndPointTmp[j];
        }
    }

    return nNumEndPoint;
}


/****************************************************************************/
/*

  1,4 tautomerism in 7-member ring

      /C==D             //C-D          A=DfsPath[0].at_no
   O=B     \        HO-B     \\        B=DfsPath[1].at_no
     |      E  <->     |      E        nStartAtomNeighbor2:        from A to HO
  HO-A     //        O=A     /         nStartAtomNeighborNeighbor: from B to O
      \\G-F             \\G-F


  1,4 tautomerism in 5-member ring


   O=B--C            O-B==C
     |   \\            |   \
     |     D  <->      |     D
     |   /             |   //
  HO-A==E           HO=A--E

*/
/****************************************************************************/
int Check7MembTautRing( struct tagCANON_GLOBALS *pCG,
                        inp_ATOM *atom,
                        DFS_PATH *DfsPath,
                        int nLenDfsPath,
                        int nStartAtomNeighbor,
                        int nStartAtomNeighbor2,
                        int nStartAtomNeighborNeighbor,
                        T_ENDPOINT *EndPoint,
                        int nMaxNumEndPoint,
                        T_BONDPOS  *BondPos,
                        int nMaxNumBondPos,
                        int *pnNumEndPoint,
                        int *pnNumBondPos,
                        struct BalancedNetworkStructure *pBNS,
                        struct BalancedNetworkData *pBD,
                        int num_atoms )
{
#define PATH_LEN 8

    int i, j, k, /*m,*/ nNumEndPoint, nNumEndPointTmp, nNumBondPos, nNumBondPosTmp;
    int endpoint, /*nMobile, nMobile1, nMobile2,*/ o1_at, o2_at;
    int ret;
    U_CHAR path_bonds[PATH_LEN + 1], bond_type;
    T_ENDPOINT EndPointTmp[2];
    T_BONDPOS  BondPosTmp[2 * PATH_LEN];
    ENDPOINT_INFO eif1, eif2;
    int nErr = 0;


    if (nLenDfsPath + 2 > PATH_LEN)
    {
        return -1; /*  too long path */
    }
    if (nLenDfsPath != 6 && nLenDfsPath != 4)
    {
        return -1; /*  wrong call */
    }

    nNumBondPos = *pnNumBondPos;
    nNumEndPoint = *pnNumEndPoint;
    nNumBondPosTmp = 0;
    nNumEndPointTmp = 0;
    ret = 0;

    o1_at = atom[(int) DfsPath[1].at_no].neighbor[nStartAtomNeighborNeighbor];
    o2_at = atom[(int) DfsPath[0].at_no].neighbor[nStartAtomNeighbor2];
    /*
    nMobile1 = (atom[o1_at].charge == -1) + atom[o1_at].num_H;
    nMobile2 = (atom[o2_at].charge == -1) + atom[o2_at].num_H;
    */
    if (!nGetEndpointInfo( atom, o1_at, &eif1 ) ||
        !nGetEndpointInfo( atom, o2_at, &eif2 ))
    {
        return 0;
    }

    /*  save endpoints */
    for (j = 0; j < 2; j++)
    {
        endpoint = j ? o2_at : o1_at;
        if (!atom[endpoint].endpoint)
        {
            AddAtom2num( EndPointTmp[nNumEndPointTmp].num, atom, endpoint, 2 ); /* fill out */
            AddAtom2DA( EndPointTmp[nNumEndPointTmp].num_DA, atom, endpoint, 2 );
            /*
                   nMobile  = j? nMobile2 : nMobile1;
               } else {
                   nMobile  = 0;
               }
               if ( nMobile ) {
                   EndPointTmp[nNumEndPointTmp].num[1] = (atom[endpoint].charge == -1);
                   EndPointTmp[nNumEndPointTmp].num[0] = nMobile;
                   for ( m = 0; m < T_NUM_ISOTOPIC; m ++ ) {
                       EndPointTmp[nNumEndPointTmp].num[T_NUM_NO_ISOTOPIC+m] = atom[endpoint].num_iso_H[NUM_H_ISOTOPES-m-1];
                   }
            */
        }
        else
        {
            memset( EndPointTmp + nNumEndPointTmp, 0, sizeof( EndPointTmp[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
        }
        EndPointTmp[nNumEndPointTmp].nAtomNumber = endpoint;
        EndPointTmp[nNumEndPointTmp].nGroupNumber = atom[endpoint].endpoint;
        EndPointTmp[nNumEndPointTmp].nEquNumber = 0;
        nNumEndPointTmp++;
    }

    /*  extract bonds */
    k = (int) DfsPath[1].at_no;
    bond_type = ( atom[k].bond_type[nStartAtomNeighborNeighbor] & ~BOND_MARK_ALL );
#if ( FIX_BOND23_IN_TAUT == 1 )
    bond_type = ACTUAL_ORDER( pBNS, k, nStartAtomNeighborNeighbor, bond_type );
#endif
    path_bonds[0] = bond_type;
    if (REPLACE_THE_BOND( bond_type ))
    {
        BondPosTmp[nNumBondPosTmp].nAtomNumber = k;
        BondPosTmp[nNumBondPosTmp].neighbor_index = nStartAtomNeighborNeighbor;
        nNumBondPosTmp += 2;
    }
    for (i = 1; i <= nLenDfsPath; i++)
    {
        bond_type = DfsPath[i].bond_type;
        path_bonds[i] = bond_type;
        if (REPLACE_THE_BOND( bond_type ))
        {
            BondPosTmp[nNumBondPosTmp].nAtomNumber = DfsPath[i].at_no;
            BondPosTmp[nNumBondPosTmp].neighbor_index = DfsPath[i].bond_pos;
            nNumBondPosTmp += 2;
        }
    }
    bond_type = ( atom[(int) DfsPath[0].at_no].bond_type[nStartAtomNeighbor2] & ~BOND_MARK_ALL );
#if ( FIX_BOND23_IN_TAUT == 1 )
    bond_type = ACTUAL_ORDER( pBNS, (int) DfsPath[0].at_no, nStartAtomNeighbor2, bond_type );
#endif
    path_bonds[i++] = bond_type;
    if (REPLACE_THE_BOND( bond_type ))
    {
        BondPosTmp[nNumBondPosTmp].nAtomNumber = DfsPath[0].at_no;
        BondPosTmp[nNumBondPosTmp].neighbor_index = nStartAtomNeighbor2;
        nNumBondPosTmp += 2;
    }

    if (!are_alt_bonds( path_bonds, i ))
    {
        return 0;
    }

    /* path_bonds is from at_n1 to at_n2 */
    if (!( j = are_alt_bonds( path_bonds, i ) ))
    {
        return 0;
    }
    /* j is a bond type of the last bond to o2_at, the first bond from o1_at is 2-j if j=1 or 2 */

    /* single bond at o2_at: it should have a mobile atom, o1_at should not */
    if ((j == BOND_SINGLE && ( (!atom[o2_at].endpoint && !eif2.cDonor) || (!atom[o1_at].endpoint && !eif1.cAcceptor) )) ||
        /* double bond at o2_at: it should not have a mobile atom, o1_at should */
        (j == BOND_DOUBLE && ( (!atom[o2_at].endpoint && !eif2.cAcceptor) || (!atom[o1_at].endpoint && !eif1.cDonor) ))) /* djb-rwth: addressing LLVM warnings */
    {
        return 0; /* bond pattern does not fit */
    }

    nNumBondPos = AddBondsPos( atom, BondPosTmp, nNumBondPosTmp, BondPos, nMaxNumBondPos, nNumBondPos );
    nNumEndPoint = AddEndPoints( EndPointTmp, nNumEndPointTmp, EndPoint, nMaxNumEndPoint, nNumEndPoint );

    if (nNumBondPos >= 0 && nNumEndPoint >= 0)
    {
        if ((ret = ( nNumBondPos > *pnNumBondPos ) || ( nNumEndPoint > *pnNumEndPoint ))) /* djb-rwth: addressing LLVM warning */
        {
            *pnNumBondPos = nNumBondPos;
            *pnNumEndPoint = nNumEndPoint;
        }
    }

    if (ret)
    {
        /* finally check whether the bonds allow moving the hydrogens */
        if (( atom[o1_at].endpoint != atom[o2_at].endpoint || !atom[o1_at].endpoint ))
        {

            nErr = bExistsAnyAltPath( pCG, pBNS, pBD, atom, num_atoms, o1_at, o2_at, ALT_PATH_MODE_TAUTOM );

            if (nErr <= 0)
            {
                return nErr;
            }
        }
    }

    return ret;


#undef PATH_LEN
}


/****************************************************************************/
/*
  1,5 Tautomerism in 6-member alt ring:

   /=\          /==\          N = DfsPath[0].at_no
 HN   C=O  <-> N    C-OH      C = DfsPath[3].at_no
   \=/          \\-//

*/
/****************************************************************************/

/****************************************************************************
  Check if a tautomeric 6-member ring has been found
*****************************************************************************/
int Check6MembTautRing( struct tagCANON_GLOBALS *pCG,
                        inp_ATOM *atom,
                        DFS_PATH *DfsPath,
                        int nLenDfsPath,
                        int nStartAtomNeighbor,
                        int nStartAtomNeighbor2,
                        int nStartAtomNeighborNeighbor,
                        T_ENDPOINT *EndPoint,
                        int nMaxNumEndPoint,
                        T_BONDPOS  *BondPos,
                        int nMaxNumBondPos,
                        int *pnNumEndPoint,
                        int *pnNumBondPos,
                        struct BalancedNetworkStructure *pBNS,
                        struct BalancedNetworkData *pBD,
                        int num_atoms )
{
#define PATH_LEN 4
    int i, j, k, /*m,*/ nNumBondPos, nNumEndPoint, ept, eptn;
    int nNumEndPointTmp, nNumBondPosTmp, o_at = 0, ret;
    /* int num_taut_endpoints, num_H; */
    int middle_pos;
    int nMobile, endpoint, endpoint_valence, chem_bonds_valence;
    int nMobile1, endpoint_valence1;  /*  o_at */
    int nMobile2, endpoint_valence2;  /*  n_at */
    int nxt_at;
    int n_at;
    U_CHAR path_bonds[2][PATH_LEN + 1], bond_type;
    T_ENDPOINT EndPointTmp[2];
    T_BONDPOS  BondPosTmp[4 * PATH_LEN];
    ENDPOINT_INFO eif1, eif2;

    if (nStartAtomNeighbor >= 0 || nStartAtomNeighbor2 >= 0 || nStartAtomNeighborNeighbor >= 0)
    {
        return -1; /*  wrong call */
    }
    if (nLenDfsPath != 5)
    {
        return -1; /*  wrong call */
    }

    nNumBondPos = *pnNumBondPos;
    nNumEndPoint = *pnNumEndPoint;
    /* djb-rwth: removing redundant code */
    nNumEndPointTmp = 0;
    ret = 0;
    for (ept = 0; ept < 2; ept++) /* djb-rwth: initialisation needed for num array */
        for (eptn = 0; eptn < T_NUM_NO_ISOTOPIC + T_NUM_ISOTOPIC; eptn++)
            EndPointTmp[ept].num[eptn] = 0;

    n_at = (int) DfsPath[0].at_no;   /*  -N= or -NH- atom */
    nxt_at = DfsPath[middle_pos = ( nLenDfsPath + 1 ) / 2].at_no;  /*  must have tautomeric neighbor -OH or =O or -NH2 or =NH */

    if (atom[nxt_at].valence != 3
#if ( TAUT_RINGS_ATTACH_CHAIN == 1 )
        || !atom[nxt_at].bCutVertex
#endif
        )
    {
        return 0;
    }

    for (i = 0; i < atom[nxt_at].valence; i++)
    {
        o_at = atom[nxt_at].neighbor[i];
        if (o_at != DfsPath[middle_pos - 1].at_no && o_at != DfsPath[middle_pos + 1].at_no)
        {
            break; /*  >=O or />-OH has been found */
        }
    }
    if (i == atom[nxt_at].valence)
    {
        return 0; /*  no neighboring atom >=O or />-OH */
    }
    bond_type = ( atom[nxt_at].bond_type[i] & ~BOND_MARK_ALL );
#if ( FIX_BOND23_IN_TAUT == 1 )
    bond_type = ACTUAL_ORDER( pBNS, nxt_at, i, bond_type );
#endif
    if (bond_type != BOND_SINGLE &&
        bond_type != BOND_DOUBLE &&
        bond_type != BOND_TAUTOM &&
        bond_type != BOND_ALT12NS &&
        bond_type != BOND_ALTERN)
    {
        return 0;
    }

    /*  check whether the two atoms already belong to one tautomeric group */
#if ( TAUT_IGNORE_EQL_ENDPOINTS == 1 )
    if (atom[n_at].endpoint && atom[n_at].endpoint == atom[o_at].endpoint)
    {
        return 0;
    }
#endif
    /*  check =O valence; must be 2 for O, S, Se or 3 for N */
    if (!( endpoint_valence1 = nGetEndpointInfo( atom, o_at, &eif1 ) ))
    {
        return 0; /*  n_at has been checked in MarkTautomerGroups(...) */
    }
    /*
        if ( 2 != endpoint_valence1 )
            return 0; // accept only O, S, Se
    */
    /*  check hydrogens/endpoints */
    nMobile1 = atom[o_at].num_H + ( atom[o_at].charge == -1 );
    if (bond_type == BOND_SINGLE && !eif1.cDonor && !atom[o_at].endpoint)
    {
        return 0;
    }
    /* not needed since nGetEndpointInfo returned non-zero
    if ( nMobile1 + atom[o_at].chem_bonds_valence != endpoint_valence1 )
        return 0;
    */

    if (!( endpoint_valence2 = nGetEndpointInfo( atom, n_at, &eif2 ) ))
    {
        return 0; /* should not happen here */
    }
    nMobile2 = atom[n_at].num_H + ( atom[n_at].charge == -1 );

    nMobile = 0;

    /*  can mobile group move from o_at to n_at? */
    nMobile += ( atom[o_at].endpoint || eif1.cDonor ) &&  /*  from o_at */
        bond_type != BOND_DOUBLE &&
        ( atom[n_at].endpoint ||                   /*  to n_at */
            eif2.cNeutralBondsValence > atom[n_at].valence );
    /*  can mobile group move from n_at to o_at? */
    nMobile += ( atom[n_at].endpoint || eif2.cDonor ) && /*  from n_at */
        ( atom[o_at].endpoint ||          /*  to o_at */
            eif1.cNeutralBondsValence > atom[o_at].valence ) &&
        bond_type != BOND_SINGLE;


    if (!nMobile)
    {
        return 0;
    }

    /*
    num_H = atom[n_at].num_H + atom[o_at].num_H;
    num_taut_endpoints = (0!=atom[n_at].endpoint) + (0!=atom[o_at].endpoint); // if O, N already are endpoints
    if ( num_H != 1 && num_taut_endpoints != 2 && !(num_H==2 && num_taut_endpoints >= 1) ) {
        return 0;
    }
    */
    /*  extract -OH bond */
    nNumBondPosTmp = 0;

    path_bonds[0][0] = path_bonds[1][0] = bond_type;
    if (REPLACE_THE_BOND( bond_type ))
    {
        BondPosTmp[nNumBondPosTmp].nAtomNumber = nxt_at; /*  accumulate bonds to be */
        BondPosTmp[nNumBondPosTmp].neighbor_index = i;   /*  marked as tautomeric */
        nNumBondPosTmp += 2; /*  leave room for the same bond in the opposite direction */
    }

    /*  extract other bonds */
    /* path_bonds[] contents:


                    O              OH            OH
                    ||             |             |
                   /  \          //  \          /  \\
                  ||   ||  <-->  |   ||  <-->  ||   |
                   \  /          \\  /          \  //
                    NH             N             N

        path[0]:  O=NH-=-      OH-N...         OH.N...
        path[1]   O=NH-=-      OH-N...         OH.N...
                 bonds are    all bonds       all bonds
                 single and   are either      are either
                 double       alt or taut     alt or taut
    */
    for (j = 0; j < middle_pos; j++)
    {
        for (i = 0; i < 2; i++)
        {
            /* k = i? j : middle_pos-1-j; */
            k = i ? middle_pos + j : middle_pos - 1 - j;
            /*  i=0: from O neighbor i=0: down to N, i=1: up to N */
            bond_type = DfsPath[k].bond_type;

            path_bonds[i][j + 1] = bond_type;
            if (REPLACE_THE_BOND( bond_type ))
            {
                BondPosTmp[nNumBondPosTmp].nAtomNumber = DfsPath[k].at_no;       /*  accumulate bonds to be */
                BondPosTmp[nNumBondPosTmp].neighbor_index = DfsPath[k].bond_pos; /*  marked as tautomeric */
                nNumBondPosTmp += 2;   /*  leave room for the same bond in the opposite direction */
            }
        }
    }

    if (!are_alt_bonds( path_bonds[0], middle_pos + 1 ) || !are_alt_bonds( path_bonds[1], middle_pos + 1 ))
    {
        return 0;
    }

    /* finally check whether the bonds allow moving the hydrogens */
    if (( atom[o_at].endpoint != atom[n_at].endpoint || !atom[o_at].endpoint ))
    {
        int nErr;
        nErr = bExistsAnyAltPath( pCG, pBNS, pBD, atom, num_atoms, n_at, o_at, ALT_PATH_MODE_TAUTOM );
        if (nErr <= 0)
        {
            return nErr;
        }
    }

    /*  save endpoints */
    for (j = 0; j < 2; j++)
    {
        endpoint = j ? n_at :      /*  =N-  2 */
            o_at;       /*  -OH  1 */
        if (!atom[endpoint].endpoint)
        { /* not a known endpoint */
            endpoint_valence = j ? endpoint_valence2 : endpoint_valence1;
            chem_bonds_valence = j ? eif2.cNeutralBondsValence : eif1.cNeutralBondsValence;
            /* endpoint_valence = get_endpoint_valence( atom[endpoint].el_number ); */
            nMobile = j ? nMobile2 : nMobile1;
            /* nMobile  = (atom[endpoint].charge == -1) + atom[endpoint].num_H; */
            /* if ( nMobile + atom[endpoint].chem_bonds_valence != endpoint_valence ) -- fixed 02-06-2003*/
            if (nMobile + chem_bonds_valence != endpoint_valence)
                return 0; /*  abnormal endpoint valence; ignore. */
            AddAtom2num( EndPointTmp[nNumEndPointTmp].num, atom, endpoint, 2 ); /* fill out */
            AddAtom2DA( EndPointTmp[nNumEndPointTmp].num_DA, atom, endpoint, 2 );
            /*
                        EndPointTmp[nNumEndPointTmp].num[1] = (atom[endpoint].charge == -1);
                        EndPointTmp[nNumEndPointTmp].num[0] = nMobile;
                        for ( m = 0; m < T_NUM_ISOTOPIC; m ++ ) {
                            EndPointTmp[nNumEndPointTmp].num[T_NUM_NO_ISOTOPIC+m] = atom[endpoint].num_iso_H[NUM_H_ISOTOPES-m-1];
                        }
            */
        }
        else
        { /* already an endpoint */ /* **now it is wrong:** no mobile atom/charge at this endpoint */
            memset( EndPointTmp + nNumEndPointTmp, 0, sizeof( EndPointTmp[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
        }
        EndPointTmp[nNumEndPointTmp].nAtomNumber = endpoint;
        EndPointTmp[nNumEndPointTmp].nGroupNumber = atom[endpoint].endpoint;
        EndPointTmp[nNumEndPointTmp].nEquNumber = 0;

        nNumEndPointTmp++;
    }
    /*  add collected tautomeric bonds and endpoints to the input/output data */
    nNumBondPos = AddBondsPos( atom, BondPosTmp, nNumBondPosTmp, BondPos, nMaxNumBondPos, nNumBondPos );
    nNumEndPoint = AddEndPoints( EndPointTmp, nNumEndPointTmp, EndPoint, nMaxNumEndPoint, nNumEndPoint );

    if (nNumBondPos >= 0 && nNumEndPoint >= 0)
    {
        if ((ret = ( nNumBondPos > *pnNumBondPos ) || ( nNumEndPoint > *pnNumEndPoint ))) /* djb-rwth: addressing LLVM warning */
        {
            *pnNumBondPos = nNumBondPos;
            *pnNumEndPoint = nNumEndPoint;
        }
    }

    return ret;
#undef PATH_LEN
}


#if ( TAUT_15_NON_RING == 1 )  /* post v.1 feature */

/****************************************************************************
Check (1,5) taut alt path centerpoint (unfinished) [add path checking]
****************************************************************************/
int Check15TautPathCenterpoint( inp_ATOM *atom, DFS_PATH *DfsPath, int nLenDfsPath, int jNxtNeigh,
    struct BalancedNetworkStructure *pBNS,
    struct BalancedNetworkData *pBD, int num_atoms )
{
    int nxt_at = atom[DfsPath[nLenDfsPath].at_no].neighbor[jNxtNeigh];
    /* atom[nxt_at].endpoint below allows for keto-enol -CH< or -CH2- endpoints  */
    return atom[nxt_at].endpoint || bIsCenterPointStrict( atom, nxt_at );
}

/****************************************************************************

  1,5 Tautomerism in general (unfinished) [just a copy from 6-memb case]

  AH--B==C--D==E             C may be carbon exhibiting keto-enol tautomerism
  0   1  2  3  4             as well as A or E may be previously detected such a carbon
            ^   nxt_at
            |
            +-- = nLenDfsPath

****************************************************************************/


/****************************************************************************
  check if 1,5 tautomeric path has been found
****************************************************************************/
int Check15TautPath( struct tagCANON_GLOBALS *pCG,
                     inp_ATOM *atom, DFS_PATH *DfsPath,
                     int nLenDfsPath,
                     int jNxtNeigh,
                     int nStartAtomNeighbor,
                     int nStartAtomNeighbor2,
                     int nStartAtomNeighborNeighbor,
                     T_ENDPOINT *EndPoint,
                     int nMaxNumEndPoint,
                     T_BONDPOS  *BondPos,
                     int nMaxNumBondPos,
                     int *pnNumEndPoint,
                     int *pnNumBondPos,
                     struct BalancedNetworkStructure *pBNS,
                     struct BalancedNetworkData *pBD,
                     int num_atoms )
{
#define PATH_LEN 4
    int i, j, k, /*m,*/ nNumBondPos, nNumEndPoint, cur_at, prv_at, at1, at2 /*, at3, step_at*/;
    int nNumEndPointTmp, nNumBondPosTmp, ret;
    /* int num_taut_endpoints, num_H; */
    int nMobile, endpoint, endpoint_valence, chem_bonds_valence;
    int nMobile1, endpoint_valence1;  /* start atom, at1 */
    int nMobile2, endpoint_valence2;  /* end atom,   at2 */
    /*int nMobile3, endpoint_valence3=-1;*/  /* middle atom, at3 */
    /*int nxt_at;*/
    int alt_bonds[2];
    U_CHAR /*path_bonds[2][PATH_LEN+1],*/ bond_type;
    T_ENDPOINT EndPointTmp[2];
    T_BONDPOS  BondPosTmp[4 * PATH_LEN];
    ENDPOINT_INFO eif1, eif2/*, eif3*/;

    if (nStartAtomNeighbor >= 0 || nStartAtomNeighbor2 >= 0 || nStartAtomNeighborNeighbor >= 0)
    {
        return -1; /*  wrong call */
    }
    if (nLenDfsPath != 3)
    {
        return -1; /*  wrong call */
    }

    nNumBondPos = *pnNumBondPos;
    nNumEndPoint = *pnNumEndPoint;
    /* djb-rwth: removing redundant code */
    nNumEndPointTmp = 0;
    ret = 0;

    /*-------add the last atom, nLenDfsPath=4 --*/
    j = jNxtNeigh;
    prv_at = DfsPath[nLenDfsPath].at_no;
    cur_at = atom[prv_at].neighbor[j];
    DfsPath[nLenDfsPath].bond_type = ( atom[prv_at].bond_type[j] & ~BOND_MARK_ALL );
#if ( FIX_BOND23_IN_TAUT == 1 )
    DfsPath[nLenDfsPath].bond_type = ACTUAL_ORDER( pBNS, prv_at, j, DfsPath[nLenDfsPath].bond_type );
#endif
    DfsPath[nLenDfsPath].bond_pos = j; /* fix index of the bond to the next atom */

    nLenDfsPath++;

    DfsPath[nLenDfsPath].at_no = cur_at;
    DfsPath[nLenDfsPath].bond_type = 0;
    DfsPath[nLenDfsPath].bond_pos = -1;
    /*nDfsPathPos[cur_at]                = nLenDfsPath+1;*/ /* mark with distance + 1 */
/*------------------------------------------*/
    at1 = (int) DfsPath[0].at_no;
    at2 = (int) DfsPath[nLenDfsPath].at_no;
    /*at3 = (int)DfsPath[2].at_no;*/
    if (atom[at1].endpoint && atom[at1].endpoint == atom[at2].endpoint)
    {
        /* start & end already belong to the same taut group */
        goto exit_function;  /* nothing to do */
    }

    /* check bond types along alt path */
    alt_bonds[0] = alt_bonds[1] = 0;
    for (i = 0; i < nLenDfsPath; i++)
    {
        alt_bonds[i % 2] |= IS_ALT_OR_DBLBOND( DfsPath[i].bond_type );
    }
    if (( alt_bonds[0] & alt_bonds[1] & ( BOND_SINGLE | BOND_DOUBLE ) ) ||
        ( alt_bonds[0] & BOND_WRONG ) || ( alt_bonds[1] & BOND_WRONG ))
    {
        goto exit_function; /* incompatible with alt path or wrong bonds */\
    }
    /* check possibly tautomeric endpoints at the ends */
    endpoint_valence1 = nGetEndpointInfo( atom, at1, &eif1 );
    endpoint_valence2 = nGetEndpointInfo( atom, at2, &eif2 );
#ifdef NEVER   /* do not use C-endpoint of keto-enol tautomer to find 1,5 the taut path */
    if (!endpoint_valence1 && !atom[at1].endpoint ||
        !endpoint_valence2 && !atom[at2].endpoint)
        goto exit_function; /* at least one of the end atoms cannot be an endpoint */
#endif
    if (!endpoint_valence1 || !endpoint_valence2)
        goto exit_function;  /* require both endpoints be heteroatoms */
    /*  check hydrogens/endpoints */
    nMobile1 = atom[at1].num_H + ( atom[at1].charge == -1 );
    if (!atom[at1].endpoint)
    {
        if (( alt_bonds[0] & BOND_SINGLE ) && !eif1.cDonor)
        {
            goto exit_function;
        }
        if (( alt_bonds[0] & BOND_DOUBLE ) && !eif1.cAcceptor)
        {
            goto exit_function;
        }
    }
    nMobile2 = atom[at2].num_H + ( atom[at2].charge == -1 );
    if (!atom[at2].endpoint)
    {
        if (( alt_bonds[1] & BOND_SINGLE ) && !eif2.cDonor)
        {
            goto exit_function;
        }
        if (( alt_bonds[1] & BOND_DOUBLE ) && !eif2.cAcceptor)
        {
            goto exit_function;
        }
    }

    nMobile = 0;

    /*  can mobile group move from at1=o_at to at2=n_at? */
    nMobile += ( atom[at1].endpoint || eif1.cDonor ) &&  /*  from o_at */
        !( alt_bonds[0] & BOND_DOUBLE ) &&
        ( atom[at2].endpoint ||                   /*  to n_at */
            eif2.cNeutralBondsValence > atom[at2].valence );
    /*  can mobile group move from at2=n_at to at1=o_at? */
    nMobile += ( atom[at2].endpoint || eif2.cDonor ) && /*  from n_at */
        !( alt_bonds[1] & BOND_DOUBLE ) &&
        ( atom[at1].endpoint ||          /*  to o_at */
            eif1.cNeutralBondsValence > atom[at1].valence );

    if (!nMobile)
    {
        goto exit_function;
    }

    /* check whether the bonds allow moving the hydrogens between at1 and at2 */
    if (( atom[at1].endpoint != atom[at2].endpoint || !atom[at1].endpoint ))
    {
        int nErr;
        nErr = bExistsAnyAltPath( pCG, pBNS, pBD, atom, num_atoms, at1, at2, ALT_PATH_MODE_TAUTOM );
        if (nErr <= 0)
        {
            ret = nErr;
            goto exit_function;
        }
    }

    /* save tautomeric bonds */
    nNumBondPosTmp = 0;
    for (k = 0; k < nLenDfsPath; k++)
    {
        bond_type = DfsPath[k].bond_type;
        if (REPLACE_THE_BOND( bond_type ))
        {
            BondPosTmp[nNumBondPosTmp].nAtomNumber = DfsPath[k].at_no;     /*  accumulate bonds to be */
            BondPosTmp[nNumBondPosTmp].neighbor_index = DfsPath[k].bond_pos;  /*  marked as tautomeric */
            nNumBondPosTmp += 2; /*  leave room for the same bond in opposite direction */
        }
    }
    /*  save endpoints */
    for (j = 0; j < 2; j++)
    {
        endpoint = j ? at2 : at1;
        if (!atom[endpoint].endpoint)
        {
            /* not a known endpoint */
            endpoint_valence = j ? endpoint_valence2 : endpoint_valence1;
            chem_bonds_valence = j ? eif2.cNeutralBondsValence : eif1.cNeutralBondsValence;
            /* endpoint_valence = get_endpoint_valence( atom[endpoint].el_number ); */
            nMobile = j ? nMobile2 : nMobile1;
            /* nMobile  = (atom[endpoint].charge == -1) + atom[endpoint].num_H; */
            /* if ( nMobile + atom[endpoint].chem_bonds_valence != endpoint_valence ) -- fixed 02-06-2003*/
            if (nMobile + chem_bonds_valence != endpoint_valence)
            {
                goto exit_function; /*  abnormal endpoint valence; ignore. */
            }
            AddAtom2num( EndPointTmp[nNumEndPointTmp].num, atom, endpoint, 2 ); /* fill out */
            AddAtom2DA( EndPointTmp[nNumEndPointTmp].num_DA, atom, endpoint, 2 );
        }
        else
        {
            /* already an endpoint */ /* **now it is wrong:** no mobile atom/charge at this endpoint */
            memset( EndPointTmp + nNumEndPointTmp, 0, sizeof( EndPointTmp[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
        }
        EndPointTmp[nNumEndPointTmp].nAtomNumber = endpoint;
        EndPointTmp[nNumEndPointTmp].nGroupNumber = atom[endpoint].endpoint;
        EndPointTmp[nNumEndPointTmp].nEquNumber = 0;

        nNumEndPointTmp++;
    }
    /*  add collected tautomeric bonds and endpoints to the input/output data */
    nNumBondPos = AddBondsPos( atom, BondPosTmp, nNumBondPosTmp, BondPos, nMaxNumBondPos, nNumBondPos );
    nNumEndPoint = AddEndPoints( EndPointTmp, nNumEndPointTmp, EndPoint, nMaxNumEndPoint, nNumEndPoint );

    if (nNumBondPos >= 0 && nNumEndPoint >= 0)
    {
        if ((ret = ( nNumBondPos > *pnNumBondPos ) || ( nNumEndPoint > *pnNumEndPoint ))) /* djb-rwth: addressing LLVM warning */
        {
            *pnNumBondPos = nNumBondPos;
            *pnNumEndPoint = nNumEndPoint;
        }
    }

exit_function:
    /*nDfsPathPos[DfsPath[nLenDfsPath].at_no] = 0;*/

    return ret;

#undef PATH_LEN
}
#endif  /* TAUT_15_NON_RING */


/****************************************************************************

  1,4 tautomerism in 5-member ring


   O=N2-C            O-N2=C          N1 = DfsPath[0].at_no
     |   \\            |   \         N2 = DfsPath[1].at_no
     |     D  <->      |     D
     |   /             |   //
  HO-N1=E           HO=N1-E

****************************************************************************/

/****************************************************************************
  check if a tautomeric 5-member ring (pyrazole derivatives) has been found
****************************************************************************/
int Check5MembTautRing( struct tagCANON_GLOBALS *pCG,
                        inp_ATOM *atom,
                        DFS_PATH *DfsPath,
                        int nLenDfsPath,
                        int nStartAtomNeighbor,
                        int nStartAtomNeighbor2,
                        int nStartAtomNeighborNeighbor,
                        T_ENDPOINT *EndPoint,
                        int nMaxNumEndPoint,
                        T_BONDPOS  *BondPos,
                        int nMaxNumBondPos,
                        int *pnNumEndPoint,
                        int *pnNumBondPos,
                        struct BalancedNetworkStructure *pBNS,
                        struct BalancedNetworkData *pBD,
                        int num_atoms )
{
#define PATH_LEN 4
    int i, j, /*m,*/ nMobile, nMobile1, nMobile2, ept, eptn;
    int num_taut_endpoints, nNumBondPos, nNumBondPosTmp, nNumEndPoint, nNumEndPointTmp, ret;
    int endpoint;
    int n1_at = (int) DfsPath[0].at_no;
    int n2_at = (int) DfsPath[1].at_no;
    U_CHAR path_bonds[PATH_LEN + 1], bond_type;
    T_ENDPOINT EndPointTmp[2];
    T_BONDPOS  BondPosTmp[2 * PATH_LEN];
    ENDPOINT_INFO eif1, eif2;

    /*  the two root atoms (atom[n1_at] and atom[n2_at]) cannot belong */
    /*  to one and the same tautomeric group: it has been verified in MarkTautomerGroups() */

    /*  check hydrogens/endpoints */
    if (nLenDfsPath != 4)
    {
        return 0; /*  program error */
    }
    if (nStartAtomNeighbor2 >= 0 || nStartAtomNeighborNeighbor >= 0)
    {
        return 0; /*  program error: wrong call */
    }

    nNumBondPos = *pnNumBondPos;
    nNumEndPoint = *pnNumEndPoint;
    nNumEndPointTmp = 0;
    /* djb-rwth: removing redundant code */
    ret = 0;
    for (ept = 0; ept < 2; ept++) /* djb-rwth: initialisation needed for num array */
        for (eptn = 0; eptn < T_NUM_NO_ISOTOPIC + T_NUM_ISOTOPIC; eptn++)
            EndPointTmp[ept].num[eptn] = 0;

    if (!nGetEndpointInfo( atom, n1_at, &eif1 ) ||
        !nGetEndpointInfo( atom, n2_at, &eif2 ))
    {
        return 0;
    }

    nMobile1 = atom[n1_at].num_H + ( atom[n1_at].charge == -1 );
    nMobile2 = atom[n2_at].num_H + ( atom[n2_at].charge == -1 );
    nMobile = nMobile1 + nMobile2;
    num_taut_endpoints = ( 0 != atom[n1_at].endpoint ) + ( 0 != atom[n2_at].endpoint ); /*  if both N atoms already are endpoints */
    /*
    if ( !(nMobile == 1 || num_taut_endpoints == 2) && !(nMobile>1 && num_taut_endpoints >= 1) ) {
        return 0;
    }
    */
    if (num_taut_endpoints == 0 && nMobile != 1)
    {
        return 0;
    }

    /* finally check whether the bonds allow moving the hydrogens */
    if (( atom[n1_at].endpoint != atom[n2_at].endpoint || !atom[n1_at].endpoint ))
    {
        int nErr;
        nErr = bExistsAnyAltPath( pCG, pBNS, pBD, atom, num_atoms, n1_at, n2_at, ALT_PATH_MODE_TAUTOM );
        if (nErr <= 0)
        {
            return nErr;
        }
    }

    /*  save endpoints */
    for (j = 0; j < 2; j++)
    {
        endpoint = j ? n1_at : n2_at;
        if (!atom[endpoint].endpoint)
        {
            /* not a known endpoint */
            /*
                        nMobile  = (atom[endpoint].charge == -1) + atom[endpoint].num_H;
                    } else {
                        nMobile  = 0;
                    }
                    if ( nMobile ) {
            */
            AddAtom2num( EndPointTmp[nNumEndPointTmp].num, atom, endpoint, 2 ); /* fill out */
            AddAtom2DA( EndPointTmp[nNumEndPointTmp].num_DA, atom, endpoint, 2 );
            /*
            EndPointTmp[nNumEndPointTmp].num[1] = (atom[endpoint].charge == -1);
            EndPointTmp[nNumEndPointTmp].num[0] = nMobile;
            for ( m = 0; m < T_NUM_ISOTOPIC; m ++ ) {
                EndPointTmp[nNumEndPointTmp].num[T_NUM_NO_ISOTOPIC+m] = atom[endpoint].num_iso_H[NUM_H_ISOTOPES-m-1];
            }
            */
        }
        else
        {
            memset( EndPointTmp + nNumEndPointTmp, 0, sizeof( EndPointTmp[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
        }
        EndPointTmp[nNumEndPointTmp].nAtomNumber = endpoint;
        EndPointTmp[nNumEndPointTmp].nGroupNumber = atom[endpoint].endpoint;
        EndPointTmp[nNumEndPointTmp].nEquNumber = 0;

        nNumEndPointTmp++;
    }

    /*  extract bonds */
    nNumBondPosTmp = 0;
    for (i = 1; i <= nLenDfsPath; i++)
    {
        bond_type = DfsPath[i].bond_type;
        path_bonds[i - 1] = bond_type;
        if (REPLACE_THE_BOND( bond_type ))
        {
            BondPosTmp[nNumBondPosTmp].nAtomNumber = DfsPath[i].at_no;
            BondPosTmp[nNumBondPosTmp].neighbor_index = DfsPath[i].bond_pos;
            nNumBondPosTmp += 2;
        }
    }
    /* path_bonds is from at_n2 to at_n1 */
    if (!( i = are_alt_bonds( path_bonds, nLenDfsPath ) ))
    {
        return 0;
    }
    /* i is a bond type of the last bond to at_n1, the first bond from at_n2 is 2-i if i=1 or 2 */

    /* single bond at n1_at: it should have a mobile atom, n2_at should not */
    if ((i == BOND_SINGLE && ( (!atom[n1_at].endpoint && !eif1.cDonor) || (!atom[n2_at].endpoint && !eif2.cAcceptor) )) ||
        /* double bond at n1_at: it should not have a mobile atom, n2_at should */
        (i == BOND_DOUBLE && ( (!atom[n1_at].endpoint && !eif1.cAcceptor) || (!atom[n2_at].endpoint && !eif2.cDonor) ))) /* djb-rwth: addressing LLVM warnings */
    {
        return 0; /* bond pattern does not fit */
    }

    nNumBondPos = AddBondsPos( atom, BondPosTmp, nNumBondPosTmp, BondPos, nMaxNumBondPos, nNumBondPos );
    nNumEndPoint = AddEndPoints( EndPointTmp, nNumEndPointTmp, EndPoint, nMaxNumEndPoint, nNumEndPoint );

    if (nNumBondPos >= 0 && nNumEndPoint >= 0)
    {
        if ((ret = ( nNumBondPos > *pnNumBondPos ) || ( nNumEndPoint > *pnNumEndPoint ))) /* djb-rwth: addressing LLVM warning */
        {
            *pnNumBondPos = nNumBondPos;
            *pnNumEndPoint = nNumEndPoint;
        }
    }
    return ret;

#undef PATH_LEN
}
#endif /* } */

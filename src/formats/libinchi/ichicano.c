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
#include <ctype.h>

#include "mode.h"

#include "ichitime.h"
#include "ichi.h"
#include "util.h"
#include "extr_ct.h"
#include "ichitaut.h"
#include "inpdef.h"
#include "ichinorm.h"
#include "ichicant.h"
#include "ichicano.h"
#include "ichicomn.h"

#include "ichicomp.h"

/****************************************************************************
 *
 *   Globals for sorting
 */

const NEIGH_LIST      *pNeighList_RankForSort; 
const ATOM_INVARIANT2 *pAtomInvariant2ForSort;
const AT_NUMB         *pNeighborsForSort;
const AT_RANK         *pn_RankForSort;

AT_RANK nMaxAtNeighRankForSort;

int nNumCompNeighborsRanksCountEql;


#define tsort insertions_sort

/* local prototypes */


void FillOutAtomInvariant( sp_ATOM* at, int num_atoms, int num_at_tg, ATOM_INVARIANT* pAtomInvariant, CANON_STAT* pCS );

int Canon_INChI1(  int num_atoms, int num_at_tg, sp_ATOM* at, CANON_STAT* pCS, INCHI_MODE nMode);
int Canon_INChI2(  int num_atoms, int num_at_tg, sp_ATOM* at, CANON_STAT* pCS, INCHI_MODE nMode);
int Canon_INChI3(int num_atoms, int num_at_tg, sp_ATOM* at, 
                 CANON_STAT* pCS, INCHI_MODE nMode, int bTautFtcn);


#ifdef COMPILE_ANSI_ONLY

static clock_t InchiClock(void);

#ifdef INCHI_USETIMES
static clock_t InchiClock(void)
{
    struct tms buf;
    clock_t c = times( &buf );
    if ( c != (clock_t)-1 ) {
        return buf.tms_utime;
    }
    return 0;
}
#else
static clock_t InchiClock(void)
{
    clock_t c = clock();
    if ( c != (clock_t)-1 ) {
        return c;
    }
    return 0;
}
#endif

#define INCHI_MSEC(X)    (long)((1000.0/(double)CLOCKS_PER_SEC)*(X))
#define INCHI_CLOCK_T(X) (clock_t)( (double)(X) / 1000.0 * (double)CLOCKS_PER_SEC )
const clock_t FullMaxClock = (clock_t)(-1);
const clock_t HalfMaxClock = (clock_t)(-1) / 2;
clock_t MaxPositiveClock = 0;
clock_t MinNegativeClock = 0;
clock_t HalfMaxPositiveClock = 0;
clock_t HalfMinNegativeClock = 0;

static void FillMaxMinClock(void); /* keep compiler happy */

static void FillMaxMinClock(void)
{
    if ( !MaxPositiveClock ) {
        clock_t valPos=0, val1 = 1;
        while ( 0 < ((val1 <<= 1), (val1 |= 1)) ) {
            valPos = val1;
        }
        MaxPositiveClock =  valPos;
        MinNegativeClock = -valPos;
        HalfMaxPositiveClock = MaxPositiveClock / 2;
        HalfMinNegativeClock = MinNegativeClock / 2;
    }
}


/******** get current process time ****************************************/
void InchiTimeGet( inchiTime *TickEnd )
{
    TickEnd->clockTime = InchiClock();
}
/******** returns difference TickEnd - TickStart in milliseconds **********/
long InchiTimeMsecDiff( inchiTime *TickEnd, inchiTime *TickStart )
{
    if ( FullMaxClock > 0 ) {
        clock_t delta;
        if ( !TickEnd || !TickStart )
            return 0;
        /* clock_t is unsigned */
        if ( TickEnd->clockTime > TickStart->clockTime ) {
            if ( TickEnd->clockTime > HalfMaxClock &&
                 TickEnd->clockTime - TickStart->clockTime > HalfMaxClock ) {
                /* overflow in TickStart->clockTime, actually TickStart->clockTime was later */
                delta = (FullMaxClock - TickEnd->clockTime) + TickStart->clockTime;
                return -INCHI_MSEC(delta);
            }
            delta = TickEnd->clockTime - TickStart->clockTime;
            return INCHI_MSEC(delta);
        } else
        if ( TickEnd->clockTime < TickStart->clockTime ) {
            if ( TickStart->clockTime > HalfMaxClock &&
                 TickStart->clockTime - TickEnd->clockTime > HalfMaxClock ) {
                /* overflow in TickEnd->clockTime, actually TickEnd->clockTime was later */
                delta = (FullMaxClock - TickStart->clockTime) + TickEnd->clockTime;
                return INCHI_MSEC(delta);
            }
            delta = TickStart->clockTime - TickEnd->clockTime;
            return -INCHI_MSEC(delta);
        }
        return 0; /* TickEnd->clockTime == TickStart->clockTime */
    } else {
        /* may happen under Win32 only where clock_t is SIGNED long */
        clock_t delta;
        FillMaxMinClock( );
        if ( !TickEnd || !TickStart )
            return 0;
        if ( (TickEnd->clockTime >= 0 && TickStart->clockTime >= 0) ||
             (TickEnd->clockTime <= 0 && TickStart->clockTime <= 0)) {
            delta = TickEnd->clockTime - TickStart->clockTime;
        } else
        if ( TickEnd->clockTime >= HalfMaxPositiveClock &&
             TickStart->clockTime <= HalfMinNegativeClock ) {
            /* end is earlier than start */
            delta = (MaxPositiveClock - TickEnd->clockTime) + (TickStart->clockTime - MinNegativeClock);
            delta = -delta;
        } else
        if ( TickEnd->clockTime <= HalfMinNegativeClock &&
             TickStart->clockTime >= HalfMaxPositiveClock ) {
            /* start was earlier than end */
            delta = (MaxPositiveClock - TickStart->clockTime) + (TickEnd->clockTime - MinNegativeClock);
        } else {
            /* there was no overflow, clock passed zero */
            delta = TickEnd->clockTime - TickStart->clockTime;
        }
        return INCHI_MSEC(delta);
    }
}
/******************* get elapsed time from TickStart ************************/
long InchiTimeElapsed( inchiTime *TickStart )
{
    inchiTime TickEnd;
    if ( !TickStart )
        return 0;
    InchiTimeGet( &TickEnd );
    return InchiTimeMsecDiff( &TickEnd, TickStart );
}
/******************* add number of milliseconds to time *********************/
void InchiTimeAddMsec( inchiTime *TickEnd, unsigned long nNumMsec )
{
    clock_t delta;
    if ( !TickEnd )
        return;
    if ( FullMaxClock > 0 ) {
        /* clock_t is unsigned */
        delta = INCHI_CLOCK_T(nNumMsec);
        TickEnd->clockTime += delta;
    } else {
        /* may happen under Win32 only where clock_t is SIGNED long */
        /* clock_t is unsigned */
        FillMaxMinClock( );
        delta = INCHI_CLOCK_T(nNumMsec);
        TickEnd->clockTime += delta;
    }
}
/******************* check whether time has expired *********************/
int bInchiTimeIsOver( inchiTime *TickStart )
{
    if ( FullMaxClock > 0 ) {
        clock_t clockCurrTime;
        if ( !TickStart )
            return 0;
        clockCurrTime = InchiClock();
        /* clock_t is unsigned */
        if ( TickStart->clockTime > clockCurrTime ) {
            if ( TickStart->clockTime > HalfMaxClock &&
                 TickStart->clockTime - clockCurrTime > HalfMaxClock ) {
                /* overflow in clockCurrTime, actually clockCurrTime was later */
                return 1;
            }
            return 0;
        } else
        if ( TickStart->clockTime < clockCurrTime ) {
            if ( clockCurrTime > HalfMaxClock &&
                 clockCurrTime - TickStart->clockTime > HalfMaxClock ) {
                /* overflow in TickStart->clockTime, actually TickStart->clockTime was later */
                return 0;
            }
            return 1;
        }
        return 0; /* TickStart->clockTime == clockCurrTime */
    } else {
        /* may happen under Win32 only where clock_t is SIGNED long */
        clock_t clockCurrTime;
        FillMaxMinClock( );
        if ( !TickStart )
            return 0;
        clockCurrTime = InchiClock();
        if ( (clockCurrTime >= 0 && TickStart->clockTime >= 0) ||
             (clockCurrTime <= 0 && TickStart->clockTime <= 0)) {
            return (clockCurrTime > TickStart->clockTime);
        } else
        if ( clockCurrTime >= HalfMaxPositiveClock &&
             TickStart->clockTime <= HalfMinNegativeClock ) {
            /* curr is earlier than start */
            return 0;
        } else
        if ( clockCurrTime <= HalfMinNegativeClock &&
             TickStart->clockTime >= HalfMaxPositiveClock ) {
            /* start was earlier than curr */
            return 1;
        } else {
            /* there was no overflow, clock passed zero */
            return (clockCurrTime > TickStart->clockTime);
        }
    }
}

#else

/******** get current process time ****************************************/
void InchiTimeGet( inchiTime *TickEnd )
{
    if ( TickEnd ) {
        struct _timeb timeb;
        _ftime( &timeb );
        TickEnd->clockTime = (unsigned long)timeb.time;
        TickEnd->millitime = (long)timeb.millitm;
    }
}
/******** returns difference TickEnd - TickStart in milliseconds **********/
long InchiTimeMsecDiff( inchiTime *TickEnd, inchiTime *TickStart )
{
    long delta;
    if ( !TickEnd || !TickStart ) {
        return 0;
    }
    if ( TickEnd->clockTime >= TickStart->clockTime ) {
        delta = (long)(TickEnd->clockTime - TickStart->clockTime);
        delta *= 1000;
        delta += TickEnd->millitime - TickStart->millitime;
    } else {
        delta = -(long)(TickStart->clockTime - TickEnd->clockTime);
        delta *= 1000;
        delta += TickEnd->millitime - TickStart->millitime;
    }
    return delta;
}
/******************* get elapsed time from TickStart ************************/
long InchiTimeElapsed( inchiTime *TickStart )
{
    inchiTime TickEnd;
    if ( !TickStart )
        return 0;
    InchiTimeGet( &TickEnd );
    return InchiTimeMsecDiff( &TickEnd, TickStart );
}
/******************* add number of milliseconds to time *********************/
void InchiTimeAddMsec( inchiTime *TickEnd, unsigned long nNumMsec )
{
    long delta;
    if ( !TickEnd )
        return;
    TickEnd->clockTime += nNumMsec / 1000;
    delta = nNumMsec % 1000 + TickEnd->millitime;
    TickEnd->clockTime += delta / 1000;
    TickEnd->millitime = delta % 1000;
}
/******************* check whether time has expired *********************/
int bInchiTimeIsOver( inchiTime *TickEnd )
{
    struct _timeb timeb;
    if ( !TickEnd )
        return 0;
    _ftime( &timeb );
    if ( TickEnd->clockTime > (unsigned long)timeb.time )
        return 0;
    if ( TickEnd->clockTime < (unsigned long)timeb.time ||
         TickEnd->millitime < (long)timeb.millitm ) {
        return 1;
    }
    return 0;
}
#endif


/****************************************************************************/
/* length of canonic representation in sizeof(AT_NUMB) units */
int GetCanonLengths( int num_at,  sp_ATOM* at, ATOM_SIZES *s, T_GROUP_INFO *t_group_info )
{ /* include taut. groups as additional "atoms" to the connection table 07-22-2002 */
    int  i, nNumCT, nNumBonds, nNumTBonds=0, nNumDblBondsStereo=0, nNumAsymCarbStereo=0, nNumIsotopic=0;
    T_GROUP *t_group = (s->nLenLinearCTTautomer && t_group_info)? t_group_info->t_group : NULL;
    for (nNumBonds = 0, i = 0; i < num_at; i ++ ) {
        nNumBonds   += at[i].valence;
        if ( at[i].iso_sort_key ) {
            nNumIsotopic ++;  /* not including tautomeric endpoints that are isotopic only due to mobile atoms */
        }

        if ( at[i].parity > 0 ) {  /* ignore hydrogen isotope parities in at[i].parity2 */
            int j = 0, nStereoBondsToAtom=0;  /* number of stereo double bonds at this atom */
            int k;
            for ( ; j < MAX_NUM_STEREO_BONDS && (k=at[i].stereo_bond_neighbor[j]); j ++ ) {
                nStereoBondsToAtom += (at[k-1].parity > 0);
            }
            nNumDblBondsStereo += nStereoBondsToAtom;
            nNumAsymCarbStereo += !j;
        }
    }
    nNumDblBondsStereo /= 2;
    nNumBonds          /= 2;

    s->nLenBonds = inchi_max( s->nLenBonds, nNumBonds );
    nNumCT     = nNumBonds; /* total number of neighbors in the CT */

#if ( CT_ATOMID != CT_ATOMID_DONTINCLUDE )
    nNumCT    += num_at;
#endif

    s->nLenCTAtOnly = inchi_max(s->nLenCTAtOnly, nNumCT);

    if ( t_group ) {
        for ( i = 0; i < t_group_info->num_t_groups; i ++ ) {
            nNumTBonds += t_group[i].nNumEndpoints;
        }
        nNumCT    += nNumTBonds;
#if ( CT_ATOMID != CT_ATOMID_DONTINCLUDE )
        nNumCT    += t_group_info->num_t_groups;
#endif
    }

    nNumCT                    = inchi_max( 1, nNumCT ); /* keep GetBaseCanonRanking() happy */
    s->nLenCT                 = inchi_max(s->nLenCT, nNumCT);
    s->nLenIsotopic           = inchi_max(s->nLenIsotopic, nNumIsotopic);
    s->nLenLinearCTStereoDble = inchi_max(s->nLenLinearCTStereoDble, nNumDblBondsStereo);
    s->nLenLinearCTStereoCarb = inchi_max(s->nLenLinearCTStereoCarb, nNumAsymCarbStereo);
    if ( t_group_info )
        s->nLenIsotopicEndpoints  = inchi_max(s->nLenIsotopicEndpoints, t_group_info->nNumIsotopicEndpoints);

    return 0;
}
/****************************************************************************/
int DeAllocateCS( CANON_STAT *pCS )
{
#define LOCAL_FREE(X) do{if(X){inchi_free(X); X=NULL;}}while(0)    

    /* connection table */
    LOCAL_FREE( pCS->LinearCT );
    LOCAL_FREE( pCS->nCanonOrd );
    LOCAL_FREE( pCS->nSymmRank );
    LOCAL_FREE( pCS->nNum_H    );
    LOCAL_FREE( pCS->nNum_H_fixed );
    LOCAL_FREE( pCS->nExchgIsoH );
    /* isotopic */
    LOCAL_FREE( pCS->LinearCTIsotopic );
    LOCAL_FREE( pCS->nSymmRankIsotopic );
    LOCAL_FREE( pCS->nCanonOrdIsotopic );
    /* isotopic tautomeric */
    LOCAL_FREE( pCS->LinearCTIsotopicTautomer );
    LOCAL_FREE( pCS->nCanonOrdIsotopicTaut );
    LOCAL_FREE( pCS->nSymmRankIsotopicTaut );
    /* stereo */
    LOCAL_FREE( pCS->LinearCTStereoDble );
    LOCAL_FREE( pCS->LinearCTStereoCarb );
    LOCAL_FREE( pCS->LinearCTStereoDbleInv );
    LOCAL_FREE( pCS->LinearCTStereoCarbInv );
    LOCAL_FREE( pCS->nCanonOrdStereo );
    LOCAL_FREE( pCS->nCanonOrdStereoInv );
    LOCAL_FREE( pCS->nCanonOrdStereoTaut );
    /* isotopic stereo */
    LOCAL_FREE( pCS->LinearCTIsotopicStereoDble );
    LOCAL_FREE( pCS->LinearCTIsotopicStereoCarb );
    LOCAL_FREE( pCS->LinearCTIsotopicStereoDbleInv );
    LOCAL_FREE( pCS->LinearCTIsotopicStereoCarbInv );
    LOCAL_FREE( pCS->bRankUsedForStereo );
    LOCAL_FREE( pCS->bAtomUsedForStereo );

    LOCAL_FREE( pCS->nCanonOrdIsotopicStereo );
    LOCAL_FREE( pCS->nCanonOrdIsotopicStereoInv );
    LOCAL_FREE( pCS->nCanonOrdIsotopicStereoTaut );
    /* tautomeric part of the connection table */
    LOCAL_FREE( pCS->LinearCTTautomer );
    LOCAL_FREE( pCS->nCanonOrdTaut );
    LOCAL_FREE( pCS->nSymmRankTaut );

    LOCAL_FREE( pCS->LinearCT2 );

    /* for establishing constitutional equivalence */
    LOCAL_FREE( pCS->nPrevAtomNumber );

    FreeNeighList( pCS->NeighList );
    pCS->NeighList = NULL;

    /* set zero lengths */
    pCS->nMaxLenLinearCTStereoDble = 0;
    pCS->nLenLinearCTStereoDble    = 0;
    pCS->nMaxLenLinearCTStereoCarb = 0;
    pCS->nLenLinearCTStereoCarb    = 0;
    pCS->nMaxLenLinearCTIsotopicStereoDble = 0;
    pCS->nLenLinearCTIsotopicStereoDble = 0;
    pCS->nMaxLenLinearCTIsotopicStereoCarb = 0;
    pCS->nLenLinearCTIsotopicStereoCarb = 0;
    pCS->nMaxLenLinearCTTautomer = 0;
    pCS->nLenLinearCTTautomer = 0;
    pCS->nMaxLenLinearCTIsotopic  = 0;
    pCS->nLenLinearCTIsotopic  = 0;
    pCS->nMaxLenLinearCTIsotopicTautomer = 0;
    pCS->nLenLinearCTIsotopicTautomer = 0;

    /* set canon numbering lengths to zero */
    pCS->nLenCanonOrd = 0;
    pCS->nLenCanonOrdIsotopic = 0;
    pCS->nLenCanonOrdIsotopicTaut = 0;
    pCS->nLenCanonOrdStereo = 0;
    pCS->nLenCanonOrdStereoTaut = 0;
    pCS->nLenCanonOrdIsotopicStereo = 0;
    pCS->nLenCanonOrdIsotopicStereoTaut = 0;
    pCS->nLenCanonOrdTaut       = 0;

    return 0;

#undef LOCAL_FREE
}
/****************************************************************************/
int AllocateCS( CANON_STAT *pCS, int num_at, int num_at_tg, int nLenCT, int nLenCTAtOnly,
                     int nLenLinearCTStereoDble, int nLenLinearCTIsotopicStereoDble,
                     int nLenLinearCTStereoCarb, int nLenLinearCTIsotopicStereoCarb,
                     int nLenLinearCTTautomer, int nLenLinearCTIsotopicTautomer,
                     int nLenIsotopic, INCHI_MODE nMode, BCN *pBCN )
{
#define pCS_CALLOC(PTR,TYPE,LEN) (pCS->PTR=(TYPE*)inchi_calloc((size_t)(LEN),sizeof(*pCS->PTR)))

    int num_err = 0;
    int num_t_groups = num_at_tg - num_at;

    pCS->nMode = nMode;
    /* connection table */
    if ( (nMode & CMODE_CT) && nLenCT > 0 ) {
        num_err += !pCS_CALLOC(LinearCT, AT_NUMB, nLenCT);
        pCS->nMaxLenLinearCT    = 
        pCS->nLenLinearCT       = nLenCT;
        pCS->nLenLinearCTAtOnly = nLenCTAtOnly;
        num_err += !pCS_CALLOC(nCanonOrd, AT_RANK, num_at_tg); 
        num_err += !pCS_CALLOC(nSymmRank, AT_RANK, num_at_tg);
        if ( pBCN ) {
            num_err += !pCS_CALLOC(nNum_H, S_CHAR, num_at);
            num_err += !pCS_CALLOC(nNum_H_fixed, S_CHAR, num_at);
            num_err += !pCS_CALLOC(nExchgIsoH, S_CHAR, num_at);
        }
    }
    /* isotopic */
    if ( (nMode & CMODE_ISO) && nLenIsotopic > 0 ) {
        num_err += !pCS_CALLOC(LinearCTIsotopic, AT_ISOTOPIC, nLenIsotopic);
        pCS->nMaxLenLinearCTIsotopic =
        pCS->nLenLinearCTIsotopic    = nLenIsotopic;
    }
    /* isotopic tautomeric */
    if ( (nMode & CMODE_ISO) && CANON_MODE_TAUT == (nMode & CANON_MODE_TAUT) ) {
        if ( nLenLinearCTIsotopicTautomer > 0 ) {
            num_err += !pCS_CALLOC(LinearCTIsotopicTautomer, AT_ISO_TGROUP, nLenLinearCTIsotopicTautomer);
            pCS->nMaxLenLinearCTIsotopicTautomer =
            pCS->nLenLinearCTIsotopicTautomer    = nLenLinearCTIsotopicTautomer;
        }
        if ( num_t_groups > 0 ) {
            num_err += !pCS_CALLOC(nCanonOrdIsotopicTaut, AT_RANK, num_t_groups);
            num_err += !pCS_CALLOC(nSymmRankIsotopicTaut, AT_RANK, num_t_groups);
        }
    }
    /* isotopic atoms & t-groups */
    if ( (nMode & CMODE_ISO) /*&& nLenIsotopic > 0*/ ||
         ((nMode & CMODE_ISO) && CANON_MODE_TAUT == (nMode & CANON_MODE_TAUT) &&
         nLenLinearCTIsotopicTautomer > 0)
        ) {
        num_err += !pCS_CALLOC(nSymmRankIsotopic, AT_RANK, num_at_tg);
        num_err += !pCS_CALLOC(nCanonOrdIsotopic, AT_RANK, num_at_tg); 
    }
    /* stereo */
    if ( (nMode & CMODE_STEREO) && nLenLinearCTStereoDble > 0 ) {
        num_err += !pCS_CALLOC(LinearCTStereoDble, AT_STEREO_DBLE, nLenLinearCTStereoDble);
        num_err += !pCS_CALLOC(LinearCTStereoDbleInv, AT_STEREO_DBLE, nLenLinearCTStereoDble);
        pCS->nLenLinearCTStereoDbleInv =
        pCS->nMaxLenLinearCTStereoDble =
        pCS->nLenLinearCTStereoDble    = nLenLinearCTStereoDble;
    }
    if ( (nMode & CMODE_STEREO) && nLenLinearCTStereoCarb > 0 ) {
        num_err += !pCS_CALLOC(LinearCTStereoCarb, AT_STEREO_CARB, nLenLinearCTStereoCarb);
        num_err += !pCS_CALLOC(LinearCTStereoCarbInv, AT_STEREO_CARB, nLenLinearCTStereoCarb);
        pCS->nLenLinearCTStereoCarbInv =
        pCS->nMaxLenLinearCTStereoCarb =
        pCS->nLenLinearCTStereoCarb    = nLenLinearCTStereoCarb;
    }
    if ( (nMode & CMODE_STEREO) && (nLenLinearCTStereoDble > 0 || nLenLinearCTStereoCarb > 0 ) ) {
        num_err += !pCS_CALLOC(nCanonOrdStereo, AT_RANK, num_at_tg);
        num_err += !pCS_CALLOC(nCanonOrdStereoInv, AT_RANK, num_at_tg);
        if ( (nMode & CMODE_TAUT) && nLenLinearCTTautomer > 0 && num_t_groups > 0 ) {
            num_err += !pCS_CALLOC(nCanonOrdStereoTaut, AT_RANK, num_t_groups);
        }
    }
    /* isotopic stereo */
    if ( (nMode & CMODE_ISO_STEREO) && nLenLinearCTIsotopicStereoDble > 0 ) {
        num_err += !pCS_CALLOC(LinearCTIsotopicStereoDble, AT_STEREO_DBLE, nLenLinearCTIsotopicStereoDble);
        num_err += !pCS_CALLOC(LinearCTIsotopicStereoDbleInv, AT_STEREO_DBLE, nLenLinearCTIsotopicStereoDble);
        pCS->nLenLinearCTIsotopicStereoDbleInv = 
        pCS->nMaxLenLinearCTIsotopicStereoDble =
        pCS->nLenLinearCTIsotopicStereoDble    = nLenLinearCTIsotopicStereoDble;
    }
    if ( (nMode & CMODE_ISO_STEREO) && nLenLinearCTIsotopicStereoCarb > 0 ) {
        num_err += !pCS_CALLOC(LinearCTIsotopicStereoCarb, AT_STEREO_CARB, nLenLinearCTIsotopicStereoCarb);
        num_err += !pCS_CALLOC(LinearCTIsotopicStereoCarbInv, AT_STEREO_CARB, nLenLinearCTIsotopicStereoCarb);
        pCS->nLenLinearCTIsotopicStereoCarbInv =
        pCS->nMaxLenLinearCTIsotopicStereoCarb =
        pCS->nLenLinearCTIsotopicStereoCarb    = nLenLinearCTIsotopicStereoCarb;
    }
    if ( (nMode & CMODE_ISO_STEREO) && (nLenLinearCTIsotopicStereoDble > 0 || nLenLinearCTIsotopicStereoCarb > 0 ) ) {
        num_err += !pCS_CALLOC(nCanonOrdIsotopicStereo, AT_RANK, num_at_tg);
        num_err += !pCS_CALLOC(nCanonOrdIsotopicStereoInv, AT_RANK, num_at_tg);
        if ( (nMode & CMODE_TAUT) && nLenLinearCTTautomer > 0 && num_t_groups > 0 ) {
            num_err += !pCS_CALLOC(nCanonOrdIsotopicStereoTaut, AT_RANK, num_t_groups);
        }
    }
    if ( ((nMode & CMODE_STEREO) && (nLenLinearCTStereoDble > 0 || nLenLinearCTStereoCarb > 0 )) ||
         ((nMode & CMODE_ISO_STEREO) && (nLenLinearCTIsotopicStereoDble > 0 || nLenLinearCTIsotopicStereoCarb > 0 )) ) {
        num_err += !pCS_CALLOC(bRankUsedForStereo, S_CHAR, num_at);
        num_err += !pCS_CALLOC(bAtomUsedForStereo, S_CHAR, num_at);
    }
    /* tautomeric part of the connection table */
    if ( (nMode & CMODE_CT) && (nMode & CMODE_TAUT) && nLenLinearCTTautomer > 0 ) {
        num_err += !pCS_CALLOC(LinearCTTautomer, AT_TAUTOMER, nLenLinearCTTautomer);
        pCS->nMaxLenLinearCTTautomer =
        pCS->nLenLinearCTTautomer    = nLenLinearCTTautomer;
        if ( num_t_groups > 0 ) {
            num_err += !pCS_CALLOC(nCanonOrdTaut, AT_RANK, num_t_groups);
            num_err += !pCS_CALLOC(nSymmRankTaut, AT_RANK, num_t_groups);
        }
    }
    
    if ( nMode & CMODE_CT )
    num_err += !pCS_CALLOC(LinearCT2, AT_NUMB, nLenCT);

    /* for establishing constitutional equivalence */
    num_err += !pCS_CALLOC(nPrevAtomNumber, AT_RANK, num_at_tg);

    /* set canon numbering lengths to zero */
    pCS->nLenCanonOrd = 0;
    pCS->nLenCanonOrdIsotopic = 0;
    pCS->nLenCanonOrdIsotopicTaut = 0;
    pCS->nLenCanonOrdStereo = 0;
    pCS->nLenCanonOrdStereoTaut = 0;
    pCS->nLenCanonOrdIsotopicStereo = 0;
    pCS->nLenCanonOrdIsotopicStereoTaut = 0;
    pCS->nLenCanonOrdTaut       = 0;


    if ( num_err ) {
        DeAllocateCS( pCS );
        return CT_OUT_OF_RAM;  /*  <BRKPT> */
    }
    return 0;

#undef pCS_CALLOC
}
/****************************************************************************/
#define COMPARE_WITH_CT(CT, CTLEN, VALUE, CONDITION) \
    if ( CONDITION ) { \
        if ( (VALUE) CT_GREATER_THAN (CT)[CTLEN] ) \
            return 1; /* not a minimal CT */ \
        (CONDITION) = (VALUE) == (CT)[CTLEN]; \
    } \
    (CT)[CTLEN] = VALUE; \
    (CTLEN)++

#define COMPARE_WITH_CTVAL(CTVAL, VALUE, CONDITION) \
    if ( CONDITION ) { \
        if ( (VALUE) CT_GREATER_THAN (CTVAL) ) \
            return 1; /* not a minimal CT */ \
        (CONDITION) = (VALUE) == (CTVAL); \
    } \
    (CTVAL) = VALUE

#define COMPARE_WITH_CT2(CT, CTLEN, VALUE, CONDITION, OPER) \
    if ( CONDITION ) { \
        if ( (VALUE) CT_GREATER_THAN (CT)[CTLEN] ) { \
            (OPER); \
            return 1; /* not a minimal CT */ \
        } \
        (CONDITION) = (VALUE) == (CT)[CTLEN]; \
    } \
    (CT)[CTLEN] = VALUE; \
    (CTLEN)++

/****************************************************************************/
int FillIsotopicAtLinearCT(int num_atoms, sp_ATOM* at, 
                           const AT_RANK *nAtomNumber,
                           AT_ISOTOPIC *LinearCTIsotopic, 
                           int nMaxLenLinearCTIsotopic, int *pnLenLinearCTIsotopic)
{
    /* at[i].init_rank = initial ranks before canonizing */
    /* nRank[i]  = new ordering number for atoms: nRank=1,2,.. */
    /* nAtomNumber[r] = orig. atom number= 0,1,...  for r = nRank-1  */
    /* nRank[nAtomNumber[r]] = r; r = 0,1,... */
    /* nAtomNumber[nRank[i]-1] = i; */

    int  i, k, rank;
    int  nLinearCTIsotopicLen=0;

    /* the following parts of the "name" should be compared */
    /* after the connection table comparison is done */
    /* to avoid wrong difference sign. So, these parts */
    /* go to a separate buffers. */
    if ( LinearCTIsotopic && nMaxLenLinearCTIsotopic > 0 ) {
        memset( LinearCTIsotopic, 0, nMaxLenLinearCTIsotopic * sizeof(LinearCTIsotopic[0]) );
    } else {
        return 0;
    }

    /* rank = nRank[nAtomNumber[rank-1]] -- proposed atoms canon. numbers */
    for ( rank = 1; rank <= num_atoms; rank ++ ) {

        i = (int)nAtomNumber[rank-1];  /* current atom */

        /****************************************************
             add isotopic atom info to LinearCTIsotopic
        *****************************************************/

        /* if the atom itself is not isotopic then add it only if */
        /* the atom is not an endpoint AND has attached T or D or 1H.  */
        k = ( !at[i].endpoint && !(at[i].cFlags & AT_FLAG_ISO_H_POINT) && (at[i].num_iso_H[0] || at[i].num_iso_H[1] || at[i].num_iso_H[2]) );
        if ( at[i].iso_atw_diff || k ) {
            if ( CHECK_OVERFLOW(nLinearCTIsotopicLen, nMaxLenLinearCTIsotopic) )
                return CT_OVERFLOW;  /*  <BRKPT> */
            LinearCTIsotopic[nLinearCTIsotopicLen].at_num       = (AT_RANK)rank;
            LinearCTIsotopic[nLinearCTIsotopicLen].iso_atw_diff = at[i].iso_atw_diff;
            LinearCTIsotopic[nLinearCTIsotopicLen].num_1H       = (NUM_H)(k? at[i].num_iso_H[0] : 0);
            LinearCTIsotopic[nLinearCTIsotopicLen].num_D        = (NUM_H)(k? at[i].num_iso_H[1] : 0);
            LinearCTIsotopic[nLinearCTIsotopicLen].num_T        = (NUM_H)(k? at[i].num_iso_H[2] : 0);
            nLinearCTIsotopicLen ++;
        }

    } /* end of cycle over all atoms. */

    if ( LinearCTIsotopic ) {
        if ( *pnLenLinearCTIsotopic ) {
            if ( *pnLenLinearCTIsotopic != nLinearCTIsotopicLen )
                return CT_LEN_MISMATCH;  /*  <BRKPT> */
        }else
            *pnLenLinearCTIsotopic = nLinearCTIsotopicLen;
    }

    /* Return value: >0 => OK */
    return nLinearCTIsotopicLen;
}

/****************************************************************************/
int FillTautLinearCT2(int num_atoms, int num_at_tg, int bIsoTaut,
                      const AT_RANK *nRank, const AT_RANK *nAtomNumber, 
                      const AT_RANK *nSymmRank, const AT_RANK *nRankIso, 
                      const AT_RANK *nAtomNumberIso, const AT_RANK *nSymmRankIso,
                      AT_TAUTOMER   *LinearCTTautomer, 
                      int nMaxLenLinearCTTautomer, int *pnLenLinearCTTautomer,
                      AT_ISO_TGROUP *LinearCTIsotopicTautomer, 
                      int nMaxLenLinearCTIsotopicTautomer, 
                      int *pnLenLinearCTIsotopicTautomer,
                      T_GROUP_INFO *t_group_info)
{
    /* nRank[i]  = Canonical numbers of atoms,.. */
    /* nAtomNumber[r] = orig. atom number= 0,1,...  for r = nRank-1  */
    /* nRank[nAtomNumber[r]] = r; r = 0,1,... */
    /* nAtomNumber[nRank[i]-1] = i; */

    T_GROUP *t_group;

    int      i, j, len=0, g, num_num, offset, max_len = 0, len_iso=0;
    const static int max_num_num = sizeof(t_group->num)/sizeof(t_group->num[0]);
    const static int max_num_iso = sizeof(LinearCTIsotopicTautomer->num)/sizeof(LinearCTIsotopicTautomer->num[0])+T_NUM_NO_ISOTOPIC;

    /****************************************************************************
    
                   Tautomeric groups 07-22-2002, modified 12-2003
    
    ****************************************************************************/

    if ( num_at_tg > num_atoms && t_group_info && t_group_info->num_t_groups ) {
        int        num_t_groups      = t_group_info->num_t_groups;
        AT_NUMB   *tGroupNumber      = t_group_info->tGroupNumber;
        AT_NUMB   *tSymmRank         = tGroupNumber + TGSO_SYMM_RANK*num_t_groups;  /*  equivalence */
        AT_NUMB   *tiSymmRank        = tGroupNumber + TGSO_SYMM_IRANK*num_t_groups;
        AT_NUMB   *tiGroupNumber     = tGroupNumber + TGSO_SYMM_IORDER*num_t_groups;
        AT_RANK    nOffset           = (AT_RANK)num_atoms;
        /*  Fill Canonical ranks and Symmetry Ranks */
        /* memcpy( tPrevGroupNumber, tGroupNumber, num_t_groups*sizeof(tPrevGroupNumber[0])); */
        for ( i = num_atoms, j = 0; i < num_at_tg; i ++, j ++ ) {
            /* tPrevGroupNumber[j] = */
            tGroupNumber[j] = nAtomNumber[i] - nOffset;
            tSymmRank[j]    = nSymmRank[i]   - nOffset;
            if ( bIsoTaut ) {
                tiGroupNumber[j] = nAtomNumberIso[i] - nOffset;
                tiSymmRank[j]    = nSymmRankIso[i]   - nOffset;
            }
        }
        /*  Sort enpoints within each tautomeric group according to the canonical ranks */
        pn_RankForSort = nRank;
        for ( i = 0; i < num_t_groups; i ++ ) {
            qsort( t_group_info->nEndpointAtomNumber + (int)t_group_info->t_group[i].nFirstEndpointAtNoPos,
                   t_group_info->t_group[i].nNumEndpoints,
                   sizeof(t_group_info->nEndpointAtomNumber[0]),
                   CompRank );
        }
        /* fill out LinearCTTautomer */
        if ( nMaxLenLinearCTTautomer ) {
            max_len = T_GROUP_HDR_LEN * t_group_info->num_t_groups + t_group_info->nNumEndpoints+1;
            if ( max_len > nMaxLenLinearCTTautomer )
                return CT_OVERFLOW;  /*   <BRKPT> */
        }
    
        /****************************************************************
         * tautomer group format (#: is an offset)
         ****************************************************************
         *             HEADER (T_GROUP_HDR_LEN=3+3iso)
         * 0:       N = number of endpoints   ( t_group->nNumEndpoints )
         * 1:       number of mobile groups   ( t_group->num[0] )
         * 2:       number of neg. charges    ( t_group->num[1] )  {note: T_NUM_NO_ISOTOPIC=2}
         *             ENDPOINT RANKS
         * 3..N+2:  sorted tautomer group endpoint ranks; the sorting order is in
         *          t_group_info->nEndpointAtomNumber[t_group->nFirstEndpointAtNoPos+j], j=0..N-1
         *
         * End mark : N==0
         ****************************************************************/
        /* num_num = t_group_info->bIgnoreIsotopic? T_NUM_NO_ISOTOPIC : max_num_num; */
        num_num = max_num_num; /*  always include isotopic info; ignore it at the CT comparison step. */
        for ( i = 0; i < t_group_info->num_t_groups; i ++ ) {
            g = tGroupNumber[i]; /*  ith tautomeric group number in canonical order */
            t_group = t_group_info->t_group + g;
            /*******************************************************
             * Tautomer non-isotopic part: LinearCTTautomer
             *******************************************************/
            /*  check length */
            if ( CHECK_OVERFLOW(len + T_GROUP_HDR_LEN + t_group->nNumEndpoints, max_len) )
                return CT_OVERFLOW;  /*   <BRKPT> */

            /*  t_group header: number of endpoints */
            LinearCTTautomer[len++] = t_group->nNumEndpoints;
            /*  t_group header: */
            /*  (a) number of mobile groups in the t_group (number of H + number of (-) ) and */
            /*  (b) number of mobile negative charges (-) in the t_group */
            for ( j = 0; j < T_NUM_NO_ISOTOPIC; j ++ ) {
                LinearCTTautomer[len++] = t_group->num[j];
            }
            /*  t_group endpoint ranks link the group to the tautomeric endpoint atoms in the structure */
            /*  according to their ranks */
            for ( j = 0, offset = t_group->nFirstEndpointAtNoPos; j < t_group->nNumEndpoints; j ++ ) {
                LinearCTTautomer[len++] = nRank[(int)t_group_info->nEndpointAtomNumber[offset+j]];
            }
        }
        if ( nMaxLenLinearCTTautomer ) {
            LinearCTTautomer[len++] = 0; /*  or CT_INITVALUE ??? */
            if ( len != max_len ) {
                len = -len; /*  program error */ /*   <BRKPT> */
            } else
            if ( *pnLenLinearCTTautomer && *pnLenLinearCTTautomer != len ) {
                return CT_LEN_MISMATCH;
            } else {
                *pnLenLinearCTTautomer = len;
            }
        } else {
            *pnLenLinearCTTautomer = 0;
        }
        /******************************************************************
         * Isotopic Tautomeric mobile groups part: LinearCTIsotopicTautomer
         ******************************************************************/
        if ( nMaxLenLinearCTIsotopicTautomer && !t_group_info->nNumIsotopicEndpoints ) {
            for ( i = 0; i < t_group_info->num_t_groups; i ++ ) {
                g = tiGroupNumber[i]; /*  ith tautomeric group number in canonical order */
                t_group = t_group_info->t_group + g;
                /*  find if mobile hydrogens are isotopic */
                if ( !t_group->iWeight ) {
                    continue; /* no isotopic H */
                }
                if ( CHECK_OVERFLOW(len_iso, nMaxLenLinearCTIsotopicTautomer) )
                    return CT_OVERFLOW;  /*   <BRKPT> */
                for ( j = T_NUM_NO_ISOTOPIC; j < max_num_num && j < max_num_iso; j ++ ) {
                     /*  num_T, num_D, num_1H */
                    LinearCTIsotopicTautomer[len_iso].num[j-T_NUM_NO_ISOTOPIC] = t_group->num[j];
                }
                /*  link to tautomer group LinearCTTautomer[i]:  */
                LinearCTIsotopicTautomer[len_iso++].tgroup_num = (AT_NUMB)(i + 1); /*  t_group isotopic rank */
            }
        }
        if ( nMaxLenLinearCTIsotopicTautomer ) {
            if ( *pnLenLinearCTIsotopicTautomer && *pnLenLinearCTIsotopicTautomer != len_iso ) {
                return CT_LEN_MISMATCH;
            }
            *pnLenLinearCTIsotopicTautomer = len_iso;
        } else {
            *pnLenLinearCTIsotopicTautomer = 0;
        }

    }
    return len;
}
/****************************************************************************
 *
 *   Update a linear connection table out of final ranks
 */
int UpdateFullLinearCT( int num_atoms, int num_at_tg, sp_ATOM* at, AT_RANK *nRank, AT_RANK *nAtomNumber,
                        CANON_STAT* pCS, int bFirstTime )
{
    /* at[i].init_rank = initial ranks before canonizing */
    /* nRank[i]  = new ordering number for atoms: nRank=1,2,.. */
    /* nAtomNumber[r] = orig. atom number= 0,1,...  for r = nRank-1  */
    /* nRank[nAtomNumber[r]] = r; r = 0,1,... */
    /* nAtomNumber[nRank[i]-1] = i; */

    AT_NUMB nNeighborNumber[MAXVAL];
    int  i, j, k, num_neigh, rank, bCompare; /*, nRetVal; */

    T_GROUP_INFO *t_group_info        = NULL;
    T_GROUP      *t_group             = NULL;
    AT_NUMB      *nEndpointAtomNumber = NULL;

    int  nCTLen=0, nCTLenAtOnly=0;

    AT_NUMB         r_neigh;
    AT_NUMB        *LinearCT           = pCS->LinearCT;

    /* the following parts of the "name" should be compared */
    /* after the connection table comparison is done */
    /* to avoid wrong difference sign. So, these parts */
    /* go to a separate buffers. */
    /* -- currently not used at all at all -- */

#if CT_ATOMID != CT_ATOMID_DONTINCLUDE
    AT_NUMB          r0_at_type;
#endif

    bCompare = bFirstTime? 0 : 1;
    
    if ( num_at_tg > num_atoms ) {
        t_group_info        = pCS->t_group_info;
        t_group             = t_group_info->t_group;
    } else {
        t_group_info        = NULL;
        t_group             = NULL;
    }

    /**********************************************************************/
    /*                                                                    */
    /*    CYCLE 1: FILL OUT  CONNECTION TABLE(S) FOR ALL ATOMS            */
    /*      ** NOT INCLUDING ISOTOPIC ATOMS AND 1H, 2H(D), 3H(T) **       */
    /*                                                                    */
    /* rank = nRank[nAtomNumber[rank-1]] -- proposed atoms canon. numbers */
    /**********************************************************************/
    for ( rank = 1; rank <= num_atoms; rank ++ ) {
        i = (int)nAtomNumber[rank-1];  /* current atom */
#if ( CT_ATOMID == CT_ATOMID_IS_CURRANK )
        r0_at_type = (AT_NUMB)rank; /* current Rank */
#else
#if ( CT_ATOMID == CT_ATOMID_IS_INITRANK )
        r0_at_type = (AT_NUMB)at[i].init_rank; /* chemical + neighborhood ID */
#else
#if ( CT_ATOMID == CT_ATOMID_DONTINCLUDE )
#else
 #error Undefined or wrong definition of CT_ATOMID
#endif
#endif
#endif
        
        /* add atom to the CT */
#if ( CT_ATOMID != CT_ATOMID_DONTINCLUDE )
        if ( CHECK_OVERFLOW(nCTLen, pCS->nMaxLenLinearCT) )
            return CT_OVERFLOW;  /*  <BRKPT> */
        COMPARE_WITH_CT(LinearCT, nCTLen, r0_at_type, bCompare);
#endif
        /*******************************************************
             add neighbors and (if required) bonds to CT
        ********************************************************/
        
        /* sort neighbors */
        num_neigh = at[i].valence;
        for ( k = 0; k < num_neigh; k ++) {
            nNeighborNumber[k] = (AT_NUMB)k;
        }
        pNeighborsForSort = at[i].neighbor;
        pn_RankForSort    = nRank;
        insertions_sort( nNeighborNumber, (size_t)num_neigh, sizeof(nNeighborNumber[0]), CompNeighborsAT_NUMBER );

        for ( k = 0; k < num_neigh; k ++) {
            /* rank = (new current atom Rank) */
            if ( (int)(r_neigh = (AT_NUMB)nRank[(int)at[i].neighbor[(int)nNeighborNumber[k]]])
                                                              CT_NEIGH_SMALLER_THAN rank ) {
                if ( CHECK_OVERFLOW(nCTLen, pCS->nMaxLenLinearCT) )
                    return CT_OVERFLOW;  /*  <BRKPT> */
                COMPARE_WITH_CT( LinearCT, nCTLen, r_neigh, bCompare);
            }
        }
        
        /* add CT row delimiter */

    } /* end of cycle over all atoms. */

    nCTLenAtOnly = nCTLen;

    /**************************************************************
    
                Tautomeric groups 07-22-2002
    
    ***************************************************************/

    for ( rank = num_atoms + 1; rank <= num_at_tg; rank ++ ) {
        j = (int)nAtomNumber[rank-1];  /* current "atom" */
        i = j - num_atoms;             /* current t-group */
#if ( CT_ATOMID == CT_ATOMID_IS_CURRANK )
        r0_at_type = (AT_NUMB)rank; /* current Rank */
#else
#if ( CT_ATOMID == CT_ATOMID_IS_INITRANK )
        r0_at_type = (AT_NUMB)rank; /* current Rank or  (AT_NUMB)at[i].init_rank; ==> chemical + neighborhood ID */
#else
#if ( CT_ATOMID == CT_ATOMID_DONTINCLUDE )
#else
  #error Undefined or wrong definition of CT_ATOMID
#endif
#endif
#endif

        /* add atom to the CT */
#if ( CT_ATOMID != CT_ATOMID_DONTINCLUDE )
        if ( CHECK_OVERFLOW(nCTLen, pCS->nMaxLenLinearCT) )
            return CT_OVERFLOW;  /*  <BRKPT> */
        COMPARE_WITH_CT(LinearCT, nCTLen, r0_at_type, bCompare);
#endif
        
        /*******************************************************
              add neighbors and (if required) bonds to CT
        ********************************************************/
        
        /* sort endpoints */
        nEndpointAtomNumber = t_group_info->nEndpointAtomNumber+(int)t_group[i].nFirstEndpointAtNoPos;
        pn_RankForSort      = nRank;
        num_neigh           = (int)t_group[i].nNumEndpoints;
        insertions_sort( nEndpointAtomNumber, (size_t)num_neigh, sizeof(nEndpointAtomNumber[0]), CompRank);

        for ( k = 0; k < num_neigh; k ++) {
            /* rank = (new current atom Rank) */
            if ( (int)(r_neigh = (AT_NUMB)nRank[(int)nEndpointAtomNumber[k]])
                                                              CT_NEIGH_SMALLER_THAN rank ) {
                if ( CHECK_OVERFLOW(nCTLen, pCS->nMaxLenLinearCT) )
                    return CT_OVERFLOW;  /*  <BRKPT> */
                COMPARE_WITH_CT( LinearCT, nCTLen, r_neigh, bCompare);
            }
        }
    } /* end of cycle over all tautomeric groups. */

    /* compare bonds types */
    /* compare elements */

    if ( LinearCT ) {
        
        if ( pCS->nLenLinearCT ) {
            if ( pCS->nLenLinearCT != nCTLen )
                return CT_LEN_MISMATCH;  /*  <BRKPT> */
        } else {
            pCS->nLenLinearCT       = nCTLen;
        }
        
        if ( pCS->nLenLinearCT ) {
            if ( pCS->nLenLinearCTAtOnly != nCTLenAtOnly )
                return CT_LEN_MISMATCH;  /*  <BRKPT> */
        } else {
            pCS->nLenLinearCTAtOnly = nCTLenAtOnly;
        }

    }

    /* Return: 0=> identical CT; -1=> new CT is smaller than the previous one */
    return (bCompare-1);
}

/*********************************************************************************************/
/* if (*bChanged & 1) then nSymmRank has been rearranged because for some r
                           min{i: r=nSymmRank[nAtomNumber[i]]}+1 != r
   if (*bChanged & 2) then ranks nTempRank[] from nSymmRank[] differ from input nCurrRank[]
   
     on exit:
    
    nSymmRank[] have been updated if (*bChanged & 1)
    nCurrRank[] have been updated if (*bChanged & 1)
    nTempRank[] is always same as nCurrRank[]
    nAtomNumber[] have been sorted so that
        (i < j) <=> (nSymmRank[nAtomNumber[i]] <= nSymmRank[nAtomNumber[j]])
*/
int FixCanonEquivalenceInfo( int num_at_tg, AT_RANK *nSymmRank, AT_RANK *nCurrRank,
                             AT_RANK *nTempRank, AT_NUMB *nAtomNumber, int *bChanged)
{
        int nNumDiffRanks, bChangeSymmRank, bChangeCurrRank=0;
        /* sort equivalence information */
        /*
        int i;
        for ( i = 0; i < num_at_tg; i ++ ) {
            nAtomNumber[i] = i;
        }
        */
        pn_RankForSort = nSymmRank; /* minimal class representatives: min ranks for equiv. atoms */
        qsort( nAtomNumber, num_at_tg, sizeof(nAtomNumber[0]), CompRanksOrd );
        
        /* convert equivalence information nSymmRank[] into ranks array nTempRank[] */
        /* eq. info contains min. possible ranks for eq. atoms; nCurrRank contains max. possible ranks */
        nNumDiffRanks = SortedEquInfoToRanks( nSymmRank/*inp*/, nTempRank/*out*/, nAtomNumber, num_at_tg, &bChangeSymmRank );
        /* check whether nCurrRank is same as new initial ranks calculated from nSymmRank[] */
        bChangeCurrRank = memcmp( nCurrRank, nTempRank, num_at_tg*sizeof(nTempRank[0]));
        
        /*-----------------------------------------------------------------------
        if ( bChangeSymmRank || bChangeCurrRank ) {
             This is the case when the initial equitable partitioning does not produce
             constitutionally equivalent classes of atoms.
             Rebuild nSymmRank[] according to the new nCurrRank[] := nTempRank[]
             For such structures the found canonical numbers of the constitutionally equivalent atoms
             are not contiguous (see nCanonRank and nSymmRank examples below). Here arrays
             nCurrRank, nAtomNumber, and nSymmRank are changed so that later the
             contiguous canonical numbers for equivalent atoms can be obtained
             (see GetCanonRanking under
             "III. Get final canonical numbering (no stereo, no isotopic)".
            
             Example: for CAS=37520-11-9 (ID=21247: Ethane, 1,2-dicyclopropyl-),
            
                         the numbers are the "final canon. numbers, nCanonRank"
              1
            
              HC   7    5         3
               | \
               |  >CH--CH2        CH
               | /       \      / |
              HC        H2C--CH<  |
                                \ |
              2          6    8   CH
            
                                  4
            
             the arrays (arranged according to ordering in nAtomNumberTemp) are:
                                     before SortedEquInfoToRanks  after SortedRanksToEquInfo
             orig. atom nos.,nAtomNumberTemp:  {4 5 6 7 0 1 2 3}   {4 5 6 7 0 1 2 3}
             order numbers for sorted  ranks:  {0 1 2 3 4 5 6 7}   {0 1 2 3 4 5 6 7}
             canonical numbering, nCanonRank:  {1 2 5 6 3 4 7 8}   {1 2 5 6 3 4 7 8}
             constit. equivalence, nSymmRank:  {1 1 1 1 3 3 7 7}   {1 1 1 1 5 5 7 7} used later
             initial equivalence,  nCurrRank:  {6 6 6 6 6 6 8 8}   {4 4 4 4 6 6 8 8} used later
             initial numbering,  nAtomNumber:  {2 3 4 7 0 1 6 7}   {0 1 2 3 4 5 6 7} used later
             final, no stereo, no isotopic, after  III. GetCanonRanking:
             final canon. numbers, nCanonRank:                     {1 2 3 4 5 6 7 8} final
        }
        ----------------------------------------------------------------------------------*/
        if ( bChangeCurrRank ) {
            memcpy( nCurrRank,   nTempRank,       num_at_tg*sizeof(nCurrRank[0]) );
        }
        if ( bChangeSymmRank ) {
            SortedRanksToEquInfo( nSymmRank/*out*/, nTempRank/*inp*/, nAtomNumber, num_at_tg );
        }
        if ( bChanged ) {
            *bChanged = (0 != bChangeSymmRank) | 2*(0 != bChangeCurrRank);
        }
        return nNumDiffRanks;
}
/* isotopic canonicalization */
/***********************************************************************
 *
 *  Canon_INChI  (former GetCanonRankingUsingEquivInfo)
 *
 */
int Canon_INChI3(int num_atoms, int num_at_tg, sp_ATOM* at, 
                 CANON_STAT* pCS, INCHI_MODE nMode, int bTautFtcn)
{
/****************************************************************

0.    Initiation, Prepare initial ranks for GetCanonRanking()

I.    Find constitutionally equivalent atoms and possibly canonical numbering
I.1      Set tautomer=On, stereo=isotopic=Off
I.2      GetCanonRanking(): Find constitutionally equivalent atoms and possibly canonical numbering
1.3      Fix canonical equivalence info if needed (if the fix is needed then the numbering is not canonical)

II.   Get final non-isotopic canonical numbering. Simultaneously obtain non-minimal isotopic and stereo CTs
         GetCanonRanking() with pCS->bKeepSymmRank = 1
         FillOutStereoParities() (create initial stereo descriptors)
         save non-isotopic canonicalization final results
         hide isotopic and tautomeric results (for historical reasons only)


III.  Find constitutionally equivalent isotopic atoms (for isotopic stereo canonicalization)
III.1    Allocate more memory
III.2    fill allocated memory with the initial data
III.3    duplicate, save old and add isotopic info to the new pCS->t_group_info
III.4    Prepare initial isotopic ranks for GetCanonRanking()
III.5    GetCanonRanking() to Find constitutionally equivalent ISOTOPIC atoms and tautomer groups
III.6    Fix canonical isotopic equivalence information and derive ranks out of it

IV.      Prepare a second Rank/AtomNumber Stack for mapping.

V.    Optimize isotopic part (optimized)
         map_isotopic_atoms2()
         save isotopic canonical numbering

VI.   Optimize stereo descriptors (optimized)
         map_stereo_bonds4()


VII. Optimize isotopic stereo descriptors (optimized)
         SwitchAtomStereoAndIsotopicStereo()
         SetCtToIsotopicStereo()
         FillOutStereoParities()
         SetUseAtomForStereo()
         map_stereo_bonds4()

         SwitchAtomStereoAndIsotopicStereo()
         SetCtToNonIsotopicStereo()




*****************************************************************/
    
    int     nRet = 0, i, n;


    /********************************************************
              input non-stereo canonical info
     ********************************************************/
    BCN            *pBCN             = pCS->pBCN;
    FTCN           *ftcn             = pBCN->ftcn + bTautFtcn;
    
    /********************************************************
              set mode flags
     ********************************************************/
    /* tautomeric structure */
    int bTaut     = (num_at_tg > num_atoms) && pCS->t_group_info && pCS->t_group_info->num_t_groups && pCS->t_group_info->t_group;
    /* special case: induced by exchangable isotopic H inequivalence of atoms in formally non-tautomeric structure */
    int bIsoXchgH = pCS->t_group_info && pCS->t_group_info->nNumIsotopicEndpoints > 1 &&
                    pCS->t_group_info->nIsotopicEndpointAtomNumber && pCS->t_group_info->nIsotopicEndpointAtomNumber[0] &&
                    (pCS->t_group_info->bTautFlagsDone & (TG_FLAG_FOUND_ISOTOPIC_H_DONE|TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE))
                   /* && (ftcn->nCanonFlags & CANON_FLAG_ISO_TAUT_DIFF)*/;
    int bHasIsotopicCanonData = (ftcn->PartitionCtIso.AtNumber && ftcn->PartitionCtIso.Rank && ftcn->nSymmRankCtIso);
    /* bHasIsotopicCanonData==0 means
     *       (1) No isotopic atoms in the component OR 
     *       (2) the component has only exchangable isotopic H that do not change canonical numbering and equivalence.
     */
    T_GROUP_INFO *t_group_info1 = bTaut? pCS->t_group_info : NULL;
    /*int bIsoXchgH = t_group_info1 && t_group_info1->nNumIsotopicEndpoints && t_group_info1->nIsotopicEndpointAtomNumber;*/
    /* isotopic canonicalization */
    int bCanonIsotopic = bHasIsotopicCanonData && ( nMode & CMODE_ISO ) && ( pCS->LinearCTIsotopic || pCS->LinearCTIsotopicTautomer || bIsoXchgH );
    /* stereo canonicalization */
    int bCanonStereo   = ( nMode & CMODE_STEREO ) && ( pCS->LinearCTStereoDble || pCS->LinearCTStereoCarb );
    /* stereo isotopic canonicalization */
    int bCanonIsoStereo = bHasIsotopicCanonData && ( nMode & CMODE_ISO_STEREO ) && (pCS->LinearCTIsotopicStereoDble || pCS->LinearCTIsotopicStereoCarb) && bCanonIsotopic;
    int bIsoTaut = (bTaut && bCanonIsotopic);

    int            bIgnoreIsotopicInputGroups;
    int            bIgnoreIsotopicInputAtoms;
    
    AT_RANK       **pRankStack1      = pBCN->pRankStack;
    int             nRankStackLen    = pBCN->nMaxLenRankStack;
    int             num_max          = pBCN->num_max;     /* allocation lengths in *pRankStack1[] */
    NEIGH_LIST     *NeighList        = ftcn->NeighList;

    int             nNumCurrRanks    = 0;
    AT_RANK        *nTempRank        = NULL;

    AT_RANK        *nSymmRank                  = NULL;
    
    AT_RANK        *nAtomNumber                = NULL;
    AT_RANK        *nRank                      = NULL;

    AT_RANK       **pRankStack2 = NULL;
    AT_RANK        *nCanonRankStereo          = NULL;
    AT_RANK        *nCanonRankStereoInv       = NULL;
    AT_RANK        *nSymmStereo               = NULL;

    AT_RANK        *nCanonRankIsotopicStereo          = NULL;
    AT_RANK        *nCanonRankIsotopicStereoInv       = NULL;

    CUR_TREE *cur_tree = NULL;
    CUR_TREE CurrentTree;


    /*AT_ISO_TGROUP  *LinearCTIsotopicTautomer  = NULL; */

    
    CANON_STAT   CS2;
    CANON_STAT* pCS2 = &CS2;

    inchiTime   ulStartTime, ulEndTime;
    /*=========== Mode Bits (low 8 bits, bit 0 is Least Significant Bit) ===========
    
      Mode      Bits       Description                                
       '0' c    0          Only one connection table canonicalization 
       '1' C    1          Recalculate CT using fixed nSymmRank       
       '2' i    1|2        Isotopic canonicalization (internal)       
       '3' I    1|2|4      Isotopic canonicalization (output)
       '4' s    1|8        Stereo canonicalization                    
       '5' S    1|2|4|16   Stereo isotopic canonicalization           
       '6' A    1|2|4|8|16 Output All                                 

      --- high 8 bits ----
      --- obsolete, only historical interest. ------
      1-2 : 0 => at[i].init_rank from Morgan+NeighList
            1 => at[i].init_rank from Atom Invariants
            2 => at[i].init_rank from nSymmRank[]
                 (at[i].init_rank is included in LinearCT
                       depending on CT_ATOMID definition)
      3   : 1 => Get Stereo canonical info
      4   : 1 => Get Isotopic canonical info
      5   : 1 => Get Charge/Radical canonical info
    ==================================================================*/
    /*int             nOutputMode = 0;*/ /* obsolete */


    int bSwitchedAtomToIsotopic = 0;


    /* vABParityUnknown holds actual value of an internal constant signifying       */
    /* unknown parity: either the same as for undefined parity (default==standard)  */
    /*  or a specific one (non-std; requested by SLUUD switch).                     */
    int vABParityUnknown = AB_PARITY_UNDF;
    if ( 0 != ( nMode & REQ_MODE_DIFF_UU_STEREO) ) 
    {
        /* Make labels for unknown and undefined stereo different */
        vABParityUnknown = AB_PARITY_UNKN;
    }


    InchiTimeGet( &ulStartTime );
    



    *pCS2 = *pCS;  /* save input information and pointers to allocated memory */

    /* set "ignore isotopic differences in tautomer groups" true */
    if ( bTaut ) {
        /* save request for isotopic tautomeric groups */
        bIgnoreIsotopicInputGroups = pCS->t_group_info->bIgnoreIsotopic;
        pCS->t_group_info->bIgnoreIsotopic = 1;
    } else {
        bIgnoreIsotopicInputGroups = 1;
    }
    /* save request for isotopic name */
    bIgnoreIsotopicInputAtoms = pCS->bIgnoreIsotopic;
    /* set "ignore isotopic differences in atoms" true */
    pCS->bIgnoreIsotopic      = 1;


    /* save non-isotopic and isotopic canonicalization results */
    pCS->nCanonFlags = ftcn->nCanonFlags;
    /* 1. non-isotopic */

    /* linear CT, H */
    memcpy( pCS->LinearCT,        ftcn->LinearCt,             ftcn->nLenLinearCt * sizeof(pCS->LinearCT[0]) );
    if ( pCS->nNum_H && ftcn->nNumH ) {
        for ( i = 0; i < num_atoms; i ++ ) {
            pCS->nNum_H[i] = /*(S_CHAR)*/(CHAR_MASK & ftcn->nNumH[i]);
        }
    }
    if ( pCS->nNum_H_fixed && ftcn->nNumHFixH ) {
        for ( i = 0; i < num_atoms; i ++ ) {
            pCS->nNum_H_fixed[i] = /*(S_CHAR)*/(CHAR_MASK & ftcn->nNumHFixH[i]);
        }
    }
    pCS->nLenLinearCT        = ftcn->nLenLinearCt;      
    pCS->nLenLinearCTAtOnly  = ftcn->nLenLinearCtAtOnly;

    /* save non-isotopic atoms equivalence and numbering */
    if ( pCS->nSymmRank ) {
        memcpy( pCS->nSymmRank, ftcn->nSymmRankCt,          num_at_tg * sizeof(pCS->nSymmRank[0]) );
    }
    if ( pCS->nCanonOrd ) {
        memcpy( pCS->nCanonOrd, ftcn->PartitionCt.AtNumber, num_at_tg * sizeof(pCS->nCanonOrd[0]) );
        pCS->nLenCanonOrd = num_atoms;
    }
    if ( ftcn->iso_exchg_atnos && pCS->nExchgIsoH ) {
        for ( i = 0; i < num_atoms; i ++ ) {
            pCS->nExchgIsoH[i] = !ftcn->iso_exchg_atnos[i]; /* (pCS->nExchgIsoH[i]==1) => tautomeric or hetero atoms that may exchange isotopic H */
        }
    }
    /* 2. isotopic */

    if ( bCanonIsotopic ) {
        /* linear CT, num_H are same as non-isotopic */
        /* save atoms equivalence and numbering */
        if ( pCS->nSymmRankIsotopic ) {
            memcpy( pCS->nSymmRankIsotopic, ftcn->nSymmRankCtIso, num_at_tg * sizeof(pCS->nSymmRankIsotopic[0]));
        }
        if ( pCS->nCanonOrdIsotopic ) {
            memcpy( pCS->nCanonOrdIsotopic, ftcn->PartitionCtIso.AtNumber, num_at_tg * sizeof(pCS->nCanonOrdIsotopic[0]) );
            pCS->nLenCanonOrdIsotopic = num_at_tg;
        }
        nRet = FillIsotopicAtLinearCT( num_atoms, at, ftcn->PartitionCtIso.AtNumber,
                                   pCS->LinearCTIsotopic, pCS->nMaxLenLinearCTIsotopic, &pCS->nLenLinearCTIsotopic );
        if ( RETURNED_ERROR(nRet) ) {
            goto exit_function;
        }
        if ( nRet < 0 ) {
            nRet = CT_TAUCOUNT_ERR;
            goto exit_function;
        }
    } else {
        pCS->nMaxLenLinearCTIsotopic = 0;
        pCS->nMaxLenLinearCTIsotopicTautomer = 0;
    }

    /* fill out tautomeric groups, isotopic and non-isotopic tautomeric CT and t_group_info1->tGroupNumber */
    if ( bTaut ) {
        bIsoTaut = bIsoTaut && ftcn->PartitionCtIso.Rank &&
                   ftcn->PartitionCtIso.AtNumber && ftcn->nSymmRankCtIso;
        nRet = FillTautLinearCT2( num_atoms, num_at_tg, bIsoTaut,
                ftcn->PartitionCt.Rank, ftcn->PartitionCt.AtNumber, ftcn->nSymmRankCt,
                ftcn->PartitionCtIso.Rank, ftcn->PartitionCtIso.AtNumber, ftcn->nSymmRankCtIso,
                pCS->LinearCTTautomer, pCS->nMaxLenLinearCTTautomer, &pCS->nLenLinearCTTautomer,
                pCS->LinearCTIsotopicTautomer, pCS->nMaxLenLinearCTIsotopicTautomer, &pCS->nLenLinearCTIsotopicTautomer,
                t_group_info1 );

        if ( RETURNED_ERROR(nRet) ) {
            goto exit_function;
        }
        if ( nRet <= 0 ) {
            nRet = CT_TAUCOUNT_ERR;
            goto exit_function;
        } else {
            /* tautomeric groups: save non-isotopic symmetry & t_group order */
            int num_t_groups = t_group_info1->num_t_groups;
            AT_NUMB *tGroupNumber      = t_group_info1->tGroupNumber;
            AT_NUMB *tSymmRank         = tGroupNumber + TGSO_SYMM_RANK*num_t_groups;
            if ( pCS->nSymmRankTaut ) {
                memcpy( pCS->nSymmRankTaut, tSymmRank, num_t_groups * sizeof(pCS->nSymmRank[0]) ); /* fixed 5-23-02 */
            }
            if ( pCS->nCanonOrdTaut ) {
                memcpy( pCS->nCanonOrdTaut, tGroupNumber, num_t_groups * sizeof(pCS->nCanonOrdTaut[0]) );
                pCS->nLenCanonOrdTaut = num_t_groups;
            }
            if ( bCanonIsotopic /*&& pCS->nLenLinearCTIsotopicTautomer*/ ) {
                /* tautomeric groups: save isotopic symmetry & t_group order */
                /*AT_NUMB ntRankOffset       = (AT_RANK)num_atoms;*/
                AT_NUMB *tiSymmRank        = tGroupNumber + TGSO_SYMM_IRANK*num_t_groups;
                AT_NUMB *tiGroupNumber     = tGroupNumber + TGSO_SYMM_IORDER*num_t_groups;
                if ( pCS->nSymmRankIsotopicTaut ) {
                    memcpy( pCS->nSymmRankIsotopicTaut, tiSymmRank, num_t_groups * sizeof(pCS->nSymmRankIsotopicTaut[0]) );
                }
                memcpy( pCS->nCanonOrdIsotopicTaut, tiGroupNumber,  num_t_groups * sizeof(pCS->nCanonOrdIsotopicTaut[0]) );
                pCS->nLenCanonOrdIsotopicTaut = num_t_groups;
            }
        }
    }
    /* save connection table if requested */
    if ( pCS->LinearCT2 ) {
        memcpy( pCS->LinearCT2, pCS->LinearCT, sizeof(pCS->LinearCT2[0])*pCS->nLenLinearCT );
        pCS->nLenLinearCT2        = pCS->nLenLinearCT;
        pCS->nLenLinearCTAtOnly2  = pCS->nLenLinearCTAtOnly;
    }

    if ( num_atoms <= 1 ) {
        bCanonStereo    = 0;  /* a sinle atom + possibly terminal hydrogen atoms */
        if ( num_atoms < 1 || !at[0].parity2 ) {
            bCanonIsoStereo = 0;  /*  structure; for example Cl- or CH4 */
        }
    }

    if ( !bCanonStereo && !(bCanonIsotopic && bCanonIsoStereo) ) {
        goto exit_function; /* skip stereo canonicalization */
    }



    /**********************************************************
                     Mode
    ***********************************************************/
    nMode = nMode & CANON_MODE_MASK;

    /* memory allocation */

    nAtomNumber      = (AT_RANK *)qmalloc(num_max*sizeof(*nAtomNumber));
    nRank            = (AT_RANK *)qmalloc(num_max*sizeof(*nRank));
    nTempRank        = (AT_RANK *)qmalloc(num_max*sizeof(*nTempRank));
    nSymmRank        = (AT_RANK *)qmalloc(num_max*sizeof(*nSymmRank));
    /***********************************************
                0.1 Initialization
    ************************************************/
    

    if ( !NeighList || !nAtomNumber || !nTempRank ||
         !nRank     || !pCS->LinearCT ) {
        nRet = CT_OUT_OF_RAM;  /* program error */  /*  <BRKPT> */
        goto exit_function;
    }
    
    pCS->NeighList = NeighList;
    
    *pCS2 = *pCS;  /* save input information and pointers to allocated memory */

    if ( !(nMode & CMODE_NOEQ_STEREO) && (bCanonStereo || bCanonIsoStereo ) ) {
        /* will be used to discover vertex equivalences in stereo canonicalization */
        memset( &CurrentTree, 0, sizeof(CurrentTree) );
        cur_tree = &CurrentTree; 
    }


    pCS->bCmpStereo = 0;
    pCS->bCmpIsotopicStereo = 0;


    if ( bCanonStereo || bCanonIsoStereo ) {
        int ii, nn;
        
        /* stereo or isotopic canonicalization: we need a second set of ranks for mapping */
        /* (isotopic atoms or stereo can only increase nNumCurrRanks) */
        pRankStack2 = (AT_RANK **) inchi_calloc( nRankStackLen, sizeof(AT_RANK *) );
        if ( pRankStack2 ) {
            /* prepare for ranks reuse */
            for ( nn = 2; nn < nRankStackLen && pRankStack1[nn]; nn ++ ) {
                pRankStack1[nn][0] = 0; /* means ranks have to be calculated */
            }
            /* reuse memory to reduce number of allocations: */
            /* move last half of pointers from pRankStack1 to pRankStack2 */
            /* The first 2 elements will be assigned separately */
            if ( (nn = (nn-2)/2) > 0 ) {
                for ( ii = 2+nn; ii < nRankStackLen && pRankStack1[ii]; ii ++ ) {
                    pRankStack2[ii-nn] = pRankStack1[ii];
                    pRankStack1[ii] = NULL;
                }
            }
        } else {
            nRet = CT_OUT_OF_RAM;  /*  <BRKPT> */
            goto exit_function; /* program error */
        }
    }
   
    if ( bCanonStereo ) {
        
       /* *pCS2 = *pCS; */ /* save input information and pointers to allocated memory */

        /* initial ranking for non-isotopic mapping */
        memcpy( nAtomNumber,      ftcn->PartitionCt.AtNumber, num_at_tg * sizeof(nAtomNumber[0]) );
        memcpy( nRank,            ftcn->PartitionCt.Rank,     num_at_tg * sizeof(nRank[0]) );
        memcpy( nSymmRank,        ftcn->nSymmRankCt,          num_at_tg * sizeof(nSymmRank[0]) ); 
        
        /* nSymmRank changes if canonical numbers of constitutionally equivalent atoms are not contiguous */
        nNumCurrRanks = FixCanonEquivalenceInfo( num_at_tg, nSymmRank /* in&out*/,
                                           nRank, nTempRank /* out */, nAtomNumber /* in&out */, NULL);
        /* atom numbers in canonical order */
        memcpy( pCS->nPrevAtomNumber, ftcn->PartitionCt.AtNumber, num_at_tg * sizeof(nAtomNumber[0]) );

        /* fill stereo part of the connection table with initial (not optimized) parities */
        /* input
        pCS->LinearCTStereoDble       
        pCS->LinearCTStereoCarb       
        pCS->nMaxLenLinearCTStereoCarb
        pCS->nMaxLenLinearCTStereoDble
        */
        nRet = FillOutStereoParities( at, num_atoms, ftcn->PartitionCt.Rank, ftcn->PartitionCt.AtNumber,
                                      nRank, nAtomNumber, pCS, 0 /* bIsotopic */ );
        /* output
        pCS->LinearCTStereoDble       
        pCS->LinearCTStereoCarb       
        pCS2->nLenLinearCTStereoCarb
        pCS2->nLenLinearCTStereoDble
        */
        if ( RETURNED_ERROR( nRet ) ) {
            goto exit_function;
        }
        if ( nRet < 0 ) {
            nRet = CT_STEREOCOUNT_ERR;
            goto exit_function;
        }
    
        /***************************************************************
         *
         *  VI. Optimize non-isotopic stereo descriptors (optimized)
         *
         ***************************************************************/
    
        /* allocate memory for stereo canonicalization */

        if ( !nCanonRankStereo )
            nCanonRankStereo       = (AT_RANK *)  qmalloc(num_max*sizeof(*nCanonRankStereo));
        if ( !nSymmStereo && !(nMode & CMODE_NOEQ_STEREO) )
            nSymmStereo            = (AT_RANK *)  qmalloc((num_max+1)*sizeof(*nSymmStereo));
        if ( !(nMode & CMODE_NOEQ_STEREO) && 0 > CurTreeAlloc( cur_tree, num_at_tg ) ) {
            nRet = CT_OUT_OF_RAM;  /*  <BRKPT> */
            goto exit_function;
        }
        /* check allocations and assign first 2 elements of pRankStack2 */
        if ( pRankStack1 && pRankStack2 &&
             nCanonRankStereo &&
             /* nCurrRankStereo  && nAtomNumberCurrStereo &&*/
             (nSymmStereo || (nMode & CMODE_NOEQ_STEREO)) ) {
            pRankStack1[0] = pRankStack2[0] = nRank;
            pRankStack1[1] = pRankStack2[1] = nAtomNumber;
        } else {
            nRet = CT_OUT_OF_RAM;  /*  <BRKPT> */
            goto exit_function;
        }

        /****************************************************************
         *
         *  VI-A. Optimize non-isotopic non-inverted stereo descriptors
         *
         ****************************************************************/

        /* set the 1st ranks in the rest of the stack to zero: prepare for ranks reuse */
        for ( n = 2; n < nRankStackLen && pRankStack1[n]; n ++ ) {
            pRankStack1[n][0] = 0; /* means ranks have to be recalculated */
        }
        /* set the 1st ranks to zero: prepare for ranks reuse */
        for ( n = 2; n < nRankStackLen && pRankStack2[n]; n ++ ) {
            pRankStack2[n][0] = 0; /* means ranks have to be recalculated */
        }

        /* for debugging or statistics */
        pCS->lNumBreakTies    =
        pCS->lNumNeighListIter=
        pCS->lNumTotCT        =
        pCS->lNumDecreasedCT  =
        pCS->lNumRejectedCT   =
        pCS->lNumEqualCT      = 0;
        pCS->bKeepSymmRank    = 0;
        pCS->bFirstCT         = 1; /* To fill out nCanonRankStereo[] in map_stero_atoms2() */

        /******************************************************************************
             nCanonRank contains input canonical numbering
             nCanonRankStereo will be filled with a transposition of canonical numbering
               which (1) keeps connection table unchanged and
                     (2) provides minimal stereo descriptors in
                         pCS->LinearCTStereoDble (length=pCS->nLenLinearCTStereoDble)
                         pCS->LinearCTStereoCarb (length=pCS->nLenLinearCTStereoCarb)
         */
        nRet =  map_stereo_bonds4
                        ( at, num_atoms, num_at_tg, num_max, 0, ftcn->PartitionCt.Rank, ftcn->PartitionCt.AtNumber,
                          nCanonRankStereo, nSymmRank,
                          pRankStack1,  pRankStack2, nTempRank, nNumCurrRanks,
                          nSymmStereo, NeighList, pCS, cur_tree, 0 /* nNumMappedBonds */, 
                          vABParityUnknown);
        
        if ( RETURNED_ERROR( nRet ) ) {
            if ( nRet == CT_TIMEOUT_ERR )
                goto exit_function;
            else
                goto exit_function; /* program error */
        } else {
            int bFailed = 0;
            if ( !nRet ) {
                bFailed = 1; /* progrm error */
                pCS2->nLenLinearCTStereoCarb = 
                pCS->nLenLinearCTStereoCarb  = -abs(pCS->nLenLinearCTStereoCarb);
                pCS2->nLenLinearCTStereoDble =
                pCS->nLenLinearCTStereoDble  = -abs(pCS->nLenLinearCTStereoDble);
                nRet = CT_STEREOCOUNT_ERR;  /*  <BRKPT> */
                goto exit_function; /* program error */
            } else {
                /* save non-isotopic lengths */
                pCS2->nLenLinearCTStereoDble = pCS->nLenLinearCTStereoDble;
                pCS2->nLenLinearCTStereoCarb = pCS->nLenLinearCTStereoCarb;
                nRet = 0;
            }
             
            /* save stereo canonical numbering */
            if ( pCS->nCanonOrdStereo ) {
                for ( i = n = 0; i < num_at_tg; i ++ ) {
                    if ( nCanonRankStereo[i] && (int)nCanonRankStereo[i] <= num_at_tg ) {
                        pCS->nCanonOrdStereo[ (int)nCanonRankStereo[i] - 1 ] = (AT_NUMB)i;
                    } else {
                        bFailed ++;
                    }
                }
                pCS->nLenCanonOrdStereo = ( bFailed )? -num_atoms : num_atoms;
            }
            /* save stereo tautomer groups numbering */
            if ( bTaut && pCS->nCanonOrdStereoTaut ) {
                if ( 0 < (nRet = SortTautomerGroupsAndEndpoints( t_group_info1, num_atoms, num_at_tg, nCanonRankStereo ) ) ) {
                    /*non-isotopic contains symmetry ranks */
                    int num_t_groups = t_group_info1->num_t_groups;
                    AT_NUMB *tGroupNumber      = t_group_info1->tGroupNumber;
                    /*AT_NUMB *tiSymmRank        = tGroupNumber + TGSO_SYMM_IRANK*num_t_groups; */
                    memcpy( pCS->nCanonOrdStereoTaut, tGroupNumber,  num_t_groups*sizeof(pCS->nCanonOrdStereoTaut[0]) );
                    pCS->nLenCanonOrdStereoTaut = ( bFailed ) ?
                                                   -num_t_groups : num_t_groups;
                } else
                if ( RETURNED_ERROR( nRet ) ) {
                    goto exit_function;
                } else {
                    nRet = 0;
                }
                /*SortTautomerGroupsAndEndpoints( t_group_info1, nCanonRank ); */ /* ??? return to non-isotopic canonical numbering */
            }
        }
        
        /****************************************************
         *
         *  VI-B. Optimize INVERTED stereo descriptors
         *
         ****************************************************/
        if ( !nCanonRankStereoInv )
            nCanonRankStereoInv  = (AT_RANK *)  qmalloc(num_max*sizeof(*nCanonRankStereoInv));
        if ( !nCanonRankStereoInv ) {
            nRet = CT_OUT_OF_RAM;  /*  <BRKPT> */
            goto exit_function;
        }
        /* copy previous non-isotopic stereo canonicalization results to Inv initial data */
        /* assign pointers */
        pCS->LinearCTStereoDble         = pCS2->LinearCTStereoDbleInv;
        pCS->LinearCTStereoCarb         = pCS2->LinearCTStereoCarbInv;
        
        /* copy the lengths */
        pCS2->nLenLinearCTStereoDbleInv =
        pCS->nLenLinearCTStereoDbleInv  =
        pCS->nLenLinearCTStereoDble     = pCS2->nLenLinearCTStereoDble;
        
        pCS2->nLenLinearCTStereoCarbInv =
        pCS->nLenLinearCTStereoCarbInv  =
        pCS->nLenLinearCTStereoCarb     = pCS2->nLenLinearCTStereoCarb;
        
        if ( pCS->nLenLinearCTStereoDble > 0 || pCS->nLenLinearCTStereoCarb > 0 ) {
            /* copy previous results, the canonical stereo CT */
            memcpy( pCS->LinearCTStereoDble, pCS2->LinearCTStereoDble, pCS->nLenLinearCTStereoDble*sizeof(pCS->LinearCTStereoDble[0]) );
            memcpy( pCS->LinearCTStereoCarb, pCS2->LinearCTStereoCarb, pCS->nLenLinearCTStereoCarb*sizeof(pCS->LinearCTStereoCarb[0]) );
        }
        memcpy( nCanonRankStereoInv, nCanonRankStereo, num_max * sizeof(nCanonRankStereoInv[0]) );
        if ( pCS->nCanonOrdStereoInv && pCS->nCanonOrdStereo ) {
            /* in case there is nothing to invert */
            memcpy( pCS->nCanonOrdStereoInv, pCS->nCanonOrdStereo, num_at_tg*sizeof(pCS->nCanonOrdStereoInv[0]));
        }

        /******************************
         *
         * Invert stereo
         *
         ******************************/

        /*********************************************************************************
         * Create initial approximation for the minimization of the stereo descriptors:
         *  invert stereogenic atom parities, one parity in each allene, all parities in
         *  pCS->LinearCTStereoCarb and allene parities in pCS->nLenLinearCTStereoDble
         */
        nRet = InvertStereo( at, num_at_tg, nCanonRankStereo, nTempRank, pCS, 1 /* bInvertLinearCTStereo */ );
        if ( RETURNED_ERROR( nRet ) ) {
            goto exit_function;
        } else
        if ( nRet > 0 ) {
            /* InvertStereo() has done some changes */
            nRet = 0;
            /* FillOutStereoParities() has already been called to fill out these 2 LinearCTs */

            /* set the 1st ranks in the rest of the stack to zero: prepare for ranks reuse */
            for ( n = 2; n < nRankStackLen && pRankStack1[n]; n ++ ) {
                pRankStack1[n][0] = 0; /* means ranks have to be recalculated */
            }
            /* set the 1st ranks to zero: prepare for ranks reuse */
            for ( n = 2; n < nRankStackLen && pRankStack2[n]; n ++ ) {
                pRankStack2[n][0] = 0; /* means ranks have to be recalculated */
            }

            /* for debugging or statistics */
            pCS->lNumBreakTies    =
            pCS->lNumNeighListIter=
            pCS->lNumTotCT        =
            pCS->lNumDecreasedCT  =
            pCS->lNumRejectedCT   =
            pCS->lNumEqualCT      = 0;
            pCS->bKeepSymmRank    = 0;
            pCS->bFirstCT         = 1; /* To fill out nCanonRankStereo[] in map_stero_atoms2() */

            /******************************************************************************
                 ftcn->PartitionCt.Rank contains input canonical numbering
                 nCanonRankStereoInv will be filled with a transposition of canonical numbering
                   which (1) keeps connection table unchanged and
                         (2) provides minimal stereo descriptors in
                             pCS->LinearCTStereoDble (length=pCS->nLenLinearCTStereoDble)
                             pCS->LinearCTStereoCarb (length=pCS->nLenLinearCTStereoCarb)
             ******************************************************************************/
            nRet = map_stereo_bonds4
                            ( at, num_atoms, num_at_tg, num_max, 0, ftcn->PartitionCt.Rank, ftcn->PartitionCt.AtNumber,
                            nCanonRankStereoInv, nSymmRank,
                            pRankStack1,  pRankStack2, nTempRank, nNumCurrRanks, nSymmStereo,
                            NeighList, pCS, cur_tree, 0,
                            vABParityUnknown);
            if ( RETURNED_ERROR( nRet ) ) {
                if ( nRet == CT_TIMEOUT_ERR )
                    goto exit_function;
                else
                    goto exit_function; /* program error */
            } else {
                int bFailed = 0;
                if ( !nRet ) {
                    bFailed = 1; /* progrm error */
                    pCS2->nLenLinearCTStereoCarb = 
                    pCS->nLenLinearCTStereoCarb  = -abs(pCS->nLenLinearCTStereoCarb);
                    pCS2->nLenLinearCTStereoDble =
                    pCS->nLenLinearCTStereoDble  = -abs(pCS->nLenLinearCTStereoDble);
                    nRet = CT_STEREOCOUNT_ERR;  /*  <BRKPT> */
                    goto exit_function; /* program error */
                }

                /* save non-isotopic pointers & lengths for INVERTED stereo */
                pCS->nLenLinearCTStereoDbleInv  =
                pCS2->nLenLinearCTStereoDbleInv = pCS->nLenLinearCTStereoDble;
                pCS->nLenLinearCTStereoCarbInv  =
                pCS2->nLenLinearCTStereoCarbInv = pCS->nLenLinearCTStereoCarb;
                
                /* restore pointers and lengths to non-inverted stereo    */
                /*  -- this is needed for InvertStereo() back, see below  */
                pCS->LinearCTStereoDble = pCS2->LinearCTStereoDble;
                pCS->LinearCTStereoCarb = pCS2->LinearCTStereoCarb;
                pCS->nLenLinearCTStereoDble = pCS2->nLenLinearCTStereoDble;
                pCS->nLenLinearCTStereoCarb = pCS2->nLenLinearCTStereoCarb;
                /* consistency check */
                if ( pCS->nLenLinearCTStereoDbleInv != pCS->nLenLinearCTStereoDble ||
                     pCS->nLenLinearCTStereoCarbInv != pCS->nLenLinearCTStereoCarb ) {
                    nRet = CT_CALC_STEREO_ERR;
                    goto exit_function; /* program error */
                }
                /******************************
                 *
                 * Invert stereo back
                 *
                 ******************************
                 *  (make sure that pointers
                 *  pCS->LinearCTStereoCarb,
                 *  pCS->LinearCTStereoDble
                 *  and corresponding lengths
                 *  have been restored)
                 ******************************/
            /*********************************************************************************
             *  invert only stereogenic atom parities and one parity in each allene, DO NOT
             *  change parities in pCS->LinearCTStereoCarb and pCS->nLenLinearCTStereoDble
             */
                nRet = InvertStereo( at, num_at_tg, nCanonRankStereo, nTempRank, pCS, 0 );
                if ( RETURNED_ERROR( nRet ) ) {
                    goto exit_function;
                }
                nRet = 0;
                
             
                /* save stereo canonical numbering */
                if ( pCS->nCanonOrdStereoInv ) {
                    for ( i = n = 0; i < num_at_tg; i ++ ) {
                        if ( nCanonRankStereoInv[i] && (int)nCanonRankStereoInv[i] <= num_at_tg ) {
                            pCS->nCanonOrdStereoInv[ (int)nCanonRankStereoInv[i] - 1 ] = (AT_NUMB)i;
                        } else {
                            bFailed ++;
                        }
                    }
                    pCS->nLenCanonOrdStereo = ( bFailed )? -num_atoms : num_atoms;
                }

                /* compare inverted and non-inverted stereo */
                pCS->bCmpStereo = 2 + CompareLinCtStereo(
                                         pCS->LinearCTStereoDbleInv, pCS->nLenLinearCTStereoDbleInv,
                                         pCS->LinearCTStereoCarbInv, pCS->nLenLinearCTStereoCarbInv,
                                         pCS->LinearCTStereoDble,    pCS->nLenLinearCTStereoDble,
                                         pCS->LinearCTStereoCarb,    pCS->nLenLinearCTStereoCarb
                                      );

            }
        } else
        if ( 0 == nRet ) {
            /* nothing has been done, restore pointers and lengths for stereo */
            pCS->LinearCTStereoDble     = pCS2->LinearCTStereoDble;
            pCS->LinearCTStereoCarb     = pCS2->LinearCTStereoCarb;
            pCS->nLenLinearCTStereoDble = pCS2->nLenLinearCTStereoDble;
            pCS->nLenLinearCTStereoCarb = pCS2->nLenLinearCTStereoCarb;
        }
    

    }
    /* restore "ignore isotopic differences in tautomer groups" */
    if ( bTaut ) {
        /* save request for isotopic tautomeric groups */
        pCS->t_group_info->bIgnoreIsotopic = bIgnoreIsotopicInputGroups;
    }
    /* restore request for isotopic name */
    pCS->bIgnoreIsotopic = bIgnoreIsotopicInputAtoms;

    if ( bCanonIsoStereo && bCanonIsotopic ) {

        /****************************************************************
         *
         *   VII. Optimize isotopic stereo descriptors (optimized)
         *
         ****************************************************************/
        /*
        pCS->LinearCTIsotopic     = NULL;
        */

        /* initial ranking for isotopic mapping */
        memcpy( nAtomNumber,      ftcn->PartitionCtIso.AtNumber, num_at_tg * sizeof(nAtomNumber[0]) );
        memcpy( nRank,            ftcn->PartitionCtIso.Rank,     num_at_tg * sizeof(nRank[0]) );
        memcpy( nSymmRank,        ftcn->nSymmRankCtIso,          num_at_tg * sizeof(nSymmRank[0]) ); 
        
        /* nSymmRank will change if canonical numbers of of constitutionally equivalent atoms are not contiguous */
        nNumCurrRanks = FixCanonEquivalenceInfo( num_at_tg, nSymmRank /* in&out*/,
                                           nRank, nTempRank /* out */, nAtomNumber /* in&out */, NULL);

        memcpy( pCS->nPrevAtomNumber, ftcn->PartitionCtIso.AtNumber, num_at_tg * sizeof(nAtomNumber[0]) );

        /* allocate memory for optimized stereo canonicalization */
        /* for stereo canonical numbering to be found. */
        if ( !nCanonRankIsotopicStereo )
            nCanonRankIsotopicStereo       = (AT_RANK *)  qmalloc(num_max*sizeof(*nCanonRankIsotopicStereo));
        if ( !nSymmStereo && !(nMode & CMODE_NOEQ_STEREO) )
            nSymmStereo            = (AT_RANK *)  qmalloc((num_max+1)*sizeof(*nSymmStereo));
        
        if ( !(nMode & CMODE_NOEQ_STEREO) && CurTreeAlloc( cur_tree, num_at_tg ) ) {
            nRet = CT_OUT_OF_RAM;  /*  <BRKPT> */
            goto exit_function;
        }
        /* check allocations and assign first 2 elements of pRankStack2 */
        if ( pRankStack1 && pRankStack2 &&
             nCanonRankIsotopicStereo &&
             (nSymmStereo || (nMode & CMODE_NOEQ_STEREO)) ) {

            pRankStack1[0] = pRankStack2[0] = nRank; /* pRankStack1[0,1] shall be unchanged */
            pRankStack1[1] = pRankStack2[1] = nAtomNumber;
        } else {
            nRet = CT_OUT_OF_RAM;  /*  <BRKPT> */
            goto exit_function;
        }

        /******************************************************************
           Important: fill out a list of stereo atoms and bonds including
           those which are stereo due to isotopic atoms only and create
           LinearCT stereo descriptors for the canonical numbering
         ******************************************************************/

        
        /* at[] has certain members for non-isotopic and isotopic stereo; switch them */
        SwitchAtomStereoAndIsotopicStereo( at, num_atoms, &bSwitchedAtomToIsotopic );
        /* prepare stereo connection tables' pointers */
        SetCtToIsotopicStereo( pCS, pCS2 );

        nRet = FillOutStereoParities( at, num_atoms, ftcn->PartitionCtIso.Rank, ftcn->PartitionCtIso.AtNumber,
                                      nRank, nAtomNumber, pCS, 1 /* bIsotopic */);
        if (RETURNED_ERROR(nRet)) {
            goto exit_function;  /* program error */
        } else
        if ( !nRet ) {
            /* no isotopic stereo */
            pCS2->nLenLinearCTIsotopicStereoDble = pCS->nLenLinearCTIsotopicStereoDble = 0;
            pCS2->nLenLinearCTIsotopicStereoCarb = pCS->nLenLinearCTIsotopicStereoCarb = 0;
            pCS->nLenCanonOrdIsotopicStereo     = 0;   
            pCS->nLenCanonOrdIsotopicStereoTaut = 0;
            pCS2->nLenLinearCTIsotopicStereoDbleInv = pCS->nLenLinearCTIsotopicStereoDbleInv = 0;
            pCS2->nLenLinearCTIsotopicStereoCarbInv = pCS->nLenLinearCTIsotopicStereoCarbInv = 0;
            goto bypass_isotopic_stereo;
        } else {
            nRet = 0; /* not an error */
        }



        /*************************************************************
         *
         *  VII-A. Optimize non-inverted isotopic stereo descriptors
         *
         *************************************************************/

        /* set the 1st ranks in the rest of the stack to zero: prepare for ranks reuse */
        for ( n = 2; n < nRankStackLen && pRankStack1[n]; n ++ ) {
            pRankStack1[n][0] = 0; /* means ranks have to be recalculated */
        }
        /* set the 1st ranks to zero: prepare for ranks reuse */
        for ( n = 2; n < nRankStackLen && pRankStack2[n]; n ++ ) {
            pRankStack2[n][0] = 0; /* means ranks have to be recalculated */
        }

        /* for debugging or statistics */
        pCS->lNumBreakTies    =
        pCS->lNumNeighListIter=
        pCS->lNumTotCT        =
        pCS->lNumDecreasedCT  =
        pCS->lNumRejectedCT   =
        pCS->lNumEqualCT      = 0;
        pCS->bKeepSymmRank    = 0;
        pCS->bFirstCT         = 1; /* To fill out nCanonRankStereo[] in map_stero_atoms2() */

        /**************************************************************************************
          nCanonRankIsotopic contains input canonical numbering
          nCanonRankIsotopicStereo will be filled with a transposition of canonical numbering
            that  (1) keeps connection table unchanged and
                  (2) provides minimal stereo descriptors in
                      pCS->LinearCTStereoDble (length=pCS->nLenLinearCTStereoDble)
                      pCS->LinearCTStereoCarb (length=pCS->nLenLinearCTStereoCarb)
        ***************************************************************************************/
        nRet = map_stereo_bonds4
                       ( at, num_atoms, num_at_tg, num_max, 0, ftcn->PartitionCtIso.Rank, 
                       ftcn->PartitionCtIso.AtNumber,
                        nCanonRankIsotopicStereo, nSymmRank,
                        pRankStack1,  pRankStack2, nTempRank, nNumCurrRanks,
                        nSymmStereo,  NeighList, pCS, cur_tree, 0,
                        vABParityUnknown);
        if ( RETURNED_ERROR( nRet ) ) {
            goto exit_function;
        } else {
            int bFailed = 0;

            if ( !nRet ) {
                bFailed = 1; /* program error */
                pCS2->nLenLinearCTIsotopicStereoDble =
                pCS->nLenLinearCTIsotopicStereoDble  = -abs(pCS->nLenLinearCTStereoDble);
                pCS2->nLenLinearCTIsotopicStereoCarb =
                pCS->nLenLinearCTIsotopicStereoCarb  = -abs(pCS->nLenLinearCTStereoCarb);
                nRet = CT_STEREOCOUNT_ERR;  /*  <BRKPT> */
                goto exit_function; /* program error */
            } else {
                /* save isotopic lengths */
                pCS->nLenLinearCTIsotopicStereoDble  =
                pCS2->nLenLinearCTIsotopicStereoDble = pCS->nLenLinearCTStereoDble;
                pCS->nLenLinearCTIsotopicStereoCarb  =
                pCS2->nLenLinearCTIsotopicStereoCarb = pCS->nLenLinearCTStereoCarb;
             
                /* save stereo canonical numbering */
                if ( pCS->nCanonOrdIsotopicStereo ) {
                    for ( i = n = 0; i < num_at_tg; i ++ ) {
                        if ( nCanonRankIsotopicStereo[i] && (int)nCanonRankIsotopicStereo[i] <= num_at_tg ) {
                            pCS->nCanonOrdIsotopicStereo[ (int)nCanonRankIsotopicStereo[i] - 1 ] = (AT_NUMB)i;
                        } else {
                            bFailed ++;
                        }
                    }
                    pCS->nLenCanonOrdIsotopicStereo = bFailed? -num_atoms : num_atoms;
                }
                /* save stereo tautomer groups numbering */
                if ( pCS->nCanonOrdIsotopicStereoTaut ) {
                    if ( 0 < (nRet=SortTautomerGroupsAndEndpoints( t_group_info1, num_atoms, num_at_tg, nCanonRankIsotopicStereo ) ) ) {
                        /*non-isotopic contains symmetry ranks */
                        int num_t_groups = t_group_info1->num_t_groups;
                        AT_NUMB *tGroupNumber      = t_group_info1->tGroupNumber;
                        /*AT_NUMB *tiSymmRank        = tGroupNumber + TGSO_SYMM_IRANK*num_t_groups; */
                        memcpy( pCS->nCanonOrdIsotopicStereoTaut, tGroupNumber,  num_t_groups*sizeof(pCS->nCanonOrdIsotopicStereoTaut[0]) );
                        pCS->nLenCanonOrdIsotopicStereoTaut = bFailed? -num_t_groups:num_t_groups;

                        /*SortTautomerGroupsAndEndpoints( t_group_info1, nCanonRank ); */ /* ??? return to non-isotopic canonical numbering */
                    } else
                    if ( RETURNED_ERROR( nRet ) ) {
                        goto exit_function;
                    } else {
                        nRet = 0;
                    }
                }
            }
        }

        /**********************************************************
         *
         *  VII-B. Optimize INVERTED isotopic stereo descriptors
         *
         **********************************************************/
        if ( !nCanonRankIsotopicStereoInv )
            nCanonRankIsotopicStereoInv  = (AT_RANK *)  qmalloc(num_max*sizeof(*nCanonRankIsotopicStereoInv));
        if ( !nCanonRankIsotopicStereoInv ) {
            nRet = CT_OUT_OF_RAM;  /*  <BRKPT> */
            goto exit_function;
        }
        /* copy previous isotopic stereo canonicalization results to Inv initial data */
        /* assign pointers */
        pCS->LinearCTStereoDble  = pCS2->LinearCTIsotopicStereoDbleInv; /*  enable stereo */
        pCS->LinearCTStereoCarb  = pCS2->LinearCTIsotopicStereoCarbInv;

        
        /* copy the lengths */
        pCS2->nLenLinearCTIsotopicStereoDbleInv =
        pCS->nLenLinearCTStereoDbleInv  =
        pCS->nLenLinearCTStereoDble     = pCS2->nLenLinearCTIsotopicStereoDble;
        
        pCS2->nLenLinearCTIsotopicStereoCarbInv =
        pCS->nLenLinearCTStereoCarbInv  =
        pCS->nLenLinearCTStereoCarb     = pCS2->nLenLinearCTIsotopicStereoCarb;

        if ( pCS->nLenLinearCTStereoDble > 0 || pCS->nLenLinearCTStereoCarb > 0 ) {
            /* copy previous results, the canonical stereo CT */
            memcpy( pCS->LinearCTStereoDble, pCS2->LinearCTIsotopicStereoDble, pCS->nLenLinearCTStereoDble*sizeof(pCS->LinearCTStereoDble[0]) );
            memcpy( pCS->LinearCTStereoCarb, pCS2->LinearCTIsotopicStereoCarb, pCS->nLenLinearCTStereoCarb*sizeof(pCS->LinearCTStereoCarb[0]) );
        }
        memcpy( nCanonRankIsotopicStereoInv, nCanonRankIsotopicStereo, num_max * sizeof(nCanonRankIsotopicStereoInv[0]) );
        if ( pCS->nCanonOrdIsotopicStereoInv && pCS->nCanonOrdIsotopicStereo ) {
            /* in case there is nothing to invert */
            memcpy( pCS->nCanonOrdIsotopicStereoInv, pCS->nCanonOrdIsotopicStereo, num_at_tg*sizeof(pCS->nCanonOrdIsotopicStereoInv[0]));
        }

        /******************************
         *
         * Invert isotopic stereo
         *
         ******************************/

        /*********************************************************************************
         * Create initial approximation for the minimization of the stereo descriptors:
         *  invert stereogenic atom parities, one parity in each allene, all parities in
         *  pCS->LinearCTStereoCarb and allene parities in pCS->nLenLinearCTStereoDble
         */
        nRet = InvertStereo( at, num_at_tg, nCanonRankIsotopicStereo, nTempRank, pCS, 1 );
        if ( RETURNED_ERROR( nRet ) ) {
            goto exit_function;
        } else
        if ( nRet > 0 ) {
            /* InvertStereo() has done some changes */
            nRet = 0;
            /* FillOutStereoParities() has already been called to fill out these 2 LinearCTs */

            /* set the 1st ranks in the rest of the stack to zero: prepare for ranks reuse */
            for ( n = 2; n < nRankStackLen && pRankStack1[n]; n ++ ) {
                pRankStack1[n][0] = 0; /* means ranks have to be recalculated */
            }
            /* set the 1st ranks to zero: prepare for ranks reuse */
            for ( n = 2; n < nRankStackLen && pRankStack2[n]; n ++ ) {
                pRankStack2[n][0] = 0; /* means ranks have to be recalculated */
            }

            /* for debugging or statistics */
            pCS->lNumBreakTies    =
            pCS->lNumNeighListIter=
            pCS->lNumTotCT        =
            pCS->lNumDecreasedCT  =
            pCS->lNumRejectedCT   =
            pCS->lNumEqualCT      = 0;
            pCS->bKeepSymmRank    = 0;
            pCS->bFirstCT         = 1; /* To fill out nCanonRankStereo[] in map_stero_atoms2() */

            /**************************************************************************************
              nCanonRankIsotopic contains input canonical numbering
              nCanonRankIsotopicStereo will be filled with a transposition of canonical numbering
                that  (1) keeps connection table unchanged and
                      (2) provides minimal stereo descriptors in
                          pCS->LinearCTStereoDble (length=pCS->nLenLinearCTStereoDble)
                          pCS->LinearCTStereoCarb (length=pCS->nLenLinearCTStereoCarb)
            */
            nRet = map_stereo_bonds4
                           ( at, num_atoms, num_at_tg, num_max, 0, ftcn->PartitionCtIso.Rank, ftcn->PartitionCtIso.AtNumber,
                            nCanonRankIsotopicStereoInv, nSymmRank,
                            pRankStack1,  pRankStack2, nTempRank, nNumCurrRanks,
                            nSymmStereo,  NeighList, pCS, cur_tree, 0,
                            vABParityUnknown);
            if ( RETURNED_ERROR( nRet ) ) {
                if ( nRet == CT_TIMEOUT_ERR )
                    goto exit_function;
                else
                    goto exit_function; /* program error */
            } else {
                int bFailed = 0;

                if ( !nRet ) {
                    bFailed = 1; /* program error */
                    pCS2->nLenLinearCTIsotopicStereoDble =
                    pCS->nLenLinearCTIsotopicStereoDble  = -abs(pCS->nLenLinearCTStereoDble);
                    pCS2->nLenLinearCTIsotopicStereoCarb =
                    pCS->nLenLinearCTIsotopicStereoCarb  = -abs(pCS->nLenLinearCTStereoCarb);
                    nRet = CT_STEREOCOUNT_ERR;  /*  <BRKPT> */
                    goto exit_function; /* program error */
                }
                /* save isotopic pointers & lengths for INVERTED stereo */

                /* save isotopic lengths */
                pCS->nLenLinearCTIsotopicStereoDbleInv  =
                pCS2->nLenLinearCTIsotopicStereoDbleInv = pCS->nLenLinearCTStereoDble;
                pCS->nLenLinearCTIsotopicStereoCarbInv  =
                pCS2->nLenLinearCTIsotopicStereoCarbInv = pCS->nLenLinearCTStereoCarb;
         
                /* restore pointers and lengths to non-inverted isotopic stereo */
                /*  -- this is needed for InvertStereo() back, see below        */
                pCS->LinearCTStereoDble = pCS2->LinearCTIsotopicStereoDble;
                pCS->LinearCTStereoCarb = pCS2->LinearCTIsotopicStereoCarb;
                pCS->nLenLinearCTStereoDble = pCS2->nLenLinearCTIsotopicStereoDble;
                pCS->nLenLinearCTStereoCarb = pCS2->nLenLinearCTIsotopicStereoCarb;

                /* consistency check */
                if ( pCS->nLenLinearCTIsotopicStereoDbleInv != pCS->nLenLinearCTIsotopicStereoDble ||
                     pCS->nLenLinearCTIsotopicStereoCarbInv != pCS->nLenLinearCTIsotopicStereoCarb ) {
                    nRet = CT_CALC_STEREO_ERR;
                    goto exit_function; /* program error */
                }
                /******************************
                 *
                 * Invert stereo back
                 *
                 ******************************
                 *  (make sure that pointers
                 *  pCS->LinearCTStereoCarb,
                 *  pCS->LinearCTStereoDble
                 *  and corresponding lengths
                 *  have been restored)
                 ******************************/
                nRet = InvertStereo( at, num_at_tg, nCanonRankIsotopicStereo, nTempRank, pCS, 0 );
                if ( RETURNED_ERROR( nRet ) ) {
                    goto exit_function;
                }
                nRet = 0;

                /* save stereo canonical numbering */
                if ( pCS->nCanonOrdIsotopicStereoInv ) {
                    for ( i = n = 0; i < num_at_tg; i ++ ) {
                        if ( nCanonRankIsotopicStereoInv[i] && (int)nCanonRankIsotopicStereoInv[i] <= num_at_tg ) {
                            pCS->nCanonOrdIsotopicStereoInv[ (int)nCanonRankIsotopicStereoInv[i] - 1 ] = (AT_NUMB)i;
                        } else {
                            bFailed ++;
                        }
                    }
                    pCS->nLenCanonOrdIsotopicStereo = bFailed? -num_atoms : num_atoms;
                }
                /* compare inverted and non-inverted isotopic stereo */
                pCS->bCmpIsotopicStereo = 2 + CompareLinCtStereo(
                                         pCS->LinearCTIsotopicStereoDbleInv, pCS->nLenLinearCTIsotopicStereoDbleInv,
                                         pCS->LinearCTIsotopicStereoCarbInv, pCS->nLenLinearCTIsotopicStereoCarbInv,
                                         pCS->LinearCTIsotopicStereoDble,    pCS->nLenLinearCTIsotopicStereoDble,
                                         pCS->LinearCTIsotopicStereoCarb,    pCS->nLenLinearCTIsotopicStereoCarb
                                      );

            }
        } else
        if ( 0 == nRet ) {
            /* nothing has been done, restore pointers and lengths for stereo */
            pCS->LinearCTStereoDble = pCS2->LinearCTIsotopicStereoDble;
            pCS->LinearCTStereoCarb = pCS2->LinearCTIsotopicStereoCarb;
            pCS->nLenLinearCTStereoDble = pCS2->nLenLinearCTIsotopicStereoDble;
            pCS->nLenLinearCTStereoCarb = pCS2->nLenLinearCTIsotopicStereoCarb;
        }

bypass_isotopic_stereo:;  /* ???       */

        pCS->LinearCTIsotopic = pCS2->LinearCTIsotopic;
    }
        


exit_function:

    if ( bSwitchedAtomToIsotopic ) {
        SwitchAtomStereoAndIsotopicStereo( at, num_atoms, &bSwitchedAtomToIsotopic );
        SetCtToNonIsotopicStereo( pCS, pCS2 ); /* ??? */
    }

    /* restore non-isotopic connection table */
    if ( pCS->LinearCT2 ) {
        inchi_swap( (char*)&pCS->LinearCT, (char*)&pCS->LinearCT2, sizeof(pCS->LinearCT) );
        inchi_swap( (char*)&pCS->nLenLinearCT, (char*)&pCS->nLenLinearCT2, sizeof(pCS->nLenLinearCT) );
        inchi_swap( (char*)&pCS->nLenLinearCTAtOnly, (char*)&pCS->nLenLinearCTAtOnly2, sizeof(pCS->nLenLinearCTAtOnly) );
    }
    
    /* free memory */
    i = 2;
    if ( pRankStack1 ) {
        pRankStack1[0] =
        pRankStack1[1] = NULL; /* deallocated separately */
        for ( ; i < nRankStackLen && pRankStack1[i]; i ++ )
            ;
    }
    if ( pRankStack1 && pRankStack2 ) {
        for ( n = 2; n < nRankStackLen && pRankStack2[n]; n ++ ) {
            if ( i < nRankStackLen - 1 ) {
                pRankStack1[i++]   = pRankStack2[n];
            } else {
                inchi_free( pRankStack2[n] );
            }
        }
        inchi_free( pRankStack2 ); 
    }

    pCS->NeighList = NULL; /* keep the pointer in pBCN->ftcn[bTaut].NeighList for further deallocation */
    qfree ( nAtomNumber );
    qfree ( nTempRank );
    qfree ( nRank );
    qfree ( nSymmRank );

    qfree( nSymmStereo );
    CurTreeFree( cur_tree );
/* memory leak fix */
/*
    qfree ( nCurrRankIsotopicStereo );
    qfree ( nAtomNumberCurrIsotopicStereo);
*/
    qfree ( nCanonRankIsotopicStereo );
    qfree ( nCanonRankIsotopicStereoInv );

    qfree( nCanonRankStereo );
    qfree( nCanonRankStereoInv );

    InchiTimeGet( &ulEndTime );
    pCS->lTotalTime = InchiTimeMsecDiff(&ulEndTime, &ulStartTime);
    return (nRet >= -1)? num_atoms : nRet; /* cannot easily get number of ranks for now */

}

/**************************************************************************************/
int Canon_INChI(int num_atoms, int num_at_tg, sp_ATOM* at, 
                CANON_STAT* pCS, INCHI_MODE nMode, int bTautFtcn)
{
    if ( pCS->pBCN && !pCS->NeighList ) {
        /* new version */
        return Canon_INChI3(  num_atoms, num_at_tg, at, pCS, nMode, bTautFtcn);
    }
    return CT_CANON_ERR;
}

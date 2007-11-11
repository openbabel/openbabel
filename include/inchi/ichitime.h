/*
 * International Chemical Identifier (InChI)
 * Version 1
 * Software version 1.02-beta
 * August 23, 2007
 * Developed at NIST
 *
 * The InChI library and programs are free software developed under the
 * auspices of the International Union of Pure and Applied Chemistry (IUPAC);
 * you can redistribute this software and/or modify it under the terms of 
 * the GNU Lesser General Public License as published by the Free Software 
 * Foundation:
 * http://www.opensource.org/licenses/lgpl-license.php
 */


#ifndef __ICHITIME_H__
#define __ICHITIME_H__


#ifdef INCHI_ANSI_ONLY

#ifdef __FreeBSD__
#include <sys/time.h>
#endif

/* get times() */
#ifdef INCHI_USETIMES
#include <sys/times.h>
#endif

/*#include <sys/timeb.h>*/

#include <time.h>

typedef struct tagInchiTime {
    clock_t clockTime;
} inchiTime;

#else

/* Win32 _ftime(): */
#include <sys/timeb.h>

typedef struct tagInchiTime {
    unsigned long  clockTime; /* Time in seconds since midnight (00:00:00), January 1, 1970;
                                 signed long overflow expected in 2038 */
    long           millitime; /* milliseconds */

} inchiTime;

#endif


#ifdef INCHI_MAIN

#define InchiTimeGet           e_InchiTimeGet
#define InchiTimeMsecDiff      e_InchiTimeMsecDiff
#define InchiTimeAddMsec       e_InchiTimeAddMsec
#define bInchiTimeIsOver       e_bInchiTimeIsOver
#define InchiTimeElapsed       e_InchiTimeElapsed

#define FullMaxClock           e_FullMaxClock
#define HalfMaxClock           e_HalfMaxClock
#define MaxPositiveClock       e_MaxPositiveClock
#define MinNegativeClock       e_MinNegativeClock
#define HalfMaxPositiveClock   e_HalfMaxPositiveClock
#define HalfMinNegativeClock   e_HalfMinNegativeClock



#endif

#ifndef INCHI_ALL_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif


void InchiTimeGet( inchiTime *TickEnd );
long InchiTimeMsecDiff( inchiTime *TickEnd, inchiTime *TickStart );
void InchiTimeAddMsec( inchiTime *TickEnd, unsigned long nNumMsec );
int  bInchiTimeIsOver( inchiTime *TickEnd );
long InchiTimeElapsed( inchiTime *TickStart );

#ifndef INCHI_ALL_CPP
#ifdef __cplusplus
}
#endif
#endif


#endif /* __ICHITIME_H__ */


/**********************************************************************
rand.cpp - Pseudo random number generator.

Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2006 by Geoffrey R. Hutchison

This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

/* derived from sample.c
 * Pseudo-Random Number Generator Generator
 * Roger Sayle, Metaphorics LLC
 * Version 1.2, September 1998
 */
#include <openbabel/babelconfig.h>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "rand.h"

#if TIME_WITH_SYS_TIME
#include <sys/time.h>
#include <time.h>
#else
#if HAVE_SYS_TIME_H
#include <sys/time.h>
#else
#include <time.h>
#endif
#endif

#ifndef True
#define True   1
#define False  0
#endif

#define IsEven(x) (((x)&1)==0)
#define IsOdd(x)  (((x)&1)!=0)

#define BothEven(x,y) IsEven((x)|(y))
#define IsPrime(x)    (!IsEven((x))&&IsOddPrime((x)))

#define HiPart(x)   (x>>16)
#define LoPart(x)   ((x)&0xffff)

/* The maximum number of unique prime factors of a 32 bit number */
/* 2*3*5*7*11*13*17*19*23*29 = 6469693230 > 2^32 = 4294967296    */
#define MAXFACT    10

namespace OpenBabel
{
#define MAXPRIMES  256
  static int primes[MAXPRIMES] = {
    1,    2,    3,    5,    7,   11,   13,   17,   19,   23,
    29,   31,   37,   41,   43,   47,   53,   59,   61,   67,
    71,   73,   79,   83,   89,   97,  101,  103,  107,  109,
    113,  127,  131,  137,  139,  149,  151,  157,  163,  167,
    173,  179,  181,  191,  193,  197,  199,  211,  223,  227,
    229,  233,  239,  241,  251,  257,  263,  269,  271,  277,
    281,  283,  293,  307,  311,  313,  317,  331,  337,  347,
    349,  353,  359,  367,  373,  379,  383,  389,  397,  401,
    409,  419,  421,  431,  433,  439,  443,  449,  457,  461,
    463,  467,  479,  487,  491,  499,  503,  509,  521,  523,
    541,  547,  557,  563,  569,  571,  577,  587,  593,  599,
    601,  607,  613,  617,  619,  631,  641,  643,  647,  653,
    659,  661,  673,  677,  683,  691,  701,  709,  719,  727,
    733,  739,  743,  751,  757,  761,  769,  773,  787,  797,
    809,  811,  821,  823,  827,  829,  839,  853,  857,  859,
    863,  877,  881,  883,  887,  907,  911,  919,  929,  937,
    941,  947,  953,  967,  971,  977,  983,  991,  997, 1009,
    1013, 1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063,
    1069, 1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129,
    1151, 1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217,
    1223, 1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289,
    1291, 1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367,
    1373, 1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447,
    1451, 1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499,
    1511, 1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579,
    1583, 1597, 1601, 1607, 1609, 1613
  };

  static unsigned int isqrt( unsigned int val )
  {
    unsigned int temp;
    unsigned int rem;
    int i,result;

    i = 16;
    while( !(val&((unsigned int)3<<30)) && i )
      {
        val <<= 2;
        i--;
      }

    if( i )
      {
        rem = (val>>30)-1;
        val <<= 2;
        result = 1;
        i--;

        while( i )
          {
            rem = (rem<<2) | (val>>30);
            result <<= 1;
            val <<= 2;

            temp = result<<1;
            if( rem > temp )
              {
                rem -= temp|1;
                result |= 1;
              }
            i--;
          }
        return result;
      }
    else
      return 0;
  }

  static int IsOddPrime( unsigned int x )
  {
    unsigned int root;
    unsigned int i;

    root = isqrt(x);
    for( i=2; i<MAXPRIMES-1; i++ )
      {
        if( (x%primes[i]) == 0 )
          return False;
        if( (unsigned int) primes[i] >= root )
          return True;
      }

    for( i=primes[MAXPRIMES-1]; i<=root; i+=2 )
      if( (x%i) == 0 )
        return False;
    return True;
  }

  static int RelativelyPrime( unsigned int x, unsigned int y )
  {
    if( BothEven(x,y) )
      return False;

    if( IsEven(x) )
      {
        do
          {
            x >>= 1;
          }
        while( IsEven(x) );
      }
    else
      while( IsEven(y) )
        y >>= 1;

    while( x != y )
      {
        if( x > y )
          {
            x -= y;
            do
              {
                x >>= 1;
              }
            while( IsEven(x) );
          }
        else if( x < y )
          {
            y -= x;
            do
              {
                y >>= 1;
              }
            while( IsEven(y) );
          }
      }
    return( x == 1 );
  }

  static void DoubleAdd( DoubleType *x, unsigned int y )
  {
    x->lo += y;
    if( x->lo < y )
      x->hi++;
  }

  static void DoubleMultiply( unsigned int x, unsigned int y, DoubleType *z )
  {
    unsigned int x0, x1, x2, x3;
    unsigned int hx, lx;
    unsigned int hy, ly;

    hx = HiPart(x);
    lx = LoPart(x);
    hy = HiPart(y);
    ly = LoPart(y);

    x0 = lx*ly;
    x1 = lx*hy;
    x2 = hx*ly;
    x3 = hx*hy;

    x1 += HiPart(x0);
    x1 += x2;
    if( x1 < x2 )
      x3 += (1<<16);

    z->hi = HiPart(x1) + x3;
    z->lo = (LoPart(x1)<<16) + LoPart(x0);
  }

  static int LeadingZeros( unsigned int x )
  {
    static int table[256] = {
      0,1,2,2,3,3,3,3,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
      6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
      7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
      7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
      8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,
      8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,
      8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,
      8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8
    };

    if( x >= (1<<16) )
      {
        if( x >= (1<<24) )
          {
            return  8-table[x>>24];
          }
        else
          return 16-table[x>>16];
      }
    else if( x >= (1<<8) )
      {
        return 24-table[x>>8];
      }
    else
      return 32-table[x];
  }

  static unsigned int DoubleModulus( DoubleType *n,  unsigned int d )
  {
    unsigned int d1, d0;
    unsigned int r1, r0;
    unsigned int m,s;

    s = LeadingZeros(d);
    if( s > 0 )
      {
        d = d<<s;
        n->hi = (n->hi<<s) | (n->lo>>(32-s));
        n->lo = n->lo << s;
      }

    d1 = HiPart(d);
    d0 = LoPart(d);

    m = ((unsigned int)(n->hi/d1)) * d0;
    r1 = ((n->hi%d1)<<16) + HiPart(n->lo);
    if( r1 < m )
      {
        r1 += d;
        if( (r1>=d) && (r1<m) )
          r1 += d;
      }
    r1 -= m;

    m = ((unsigned int)(r1/d1)) * d0;
    r0 = ((r1%d1)<<16) + LoPart(n->lo);
    if( r0 < m )
      {
        r0 += d;
        if( (r0>=d) && (r0<m) )
          r0 += d;
      }
    r0 -= m;

    return r0 >> s;
  }

  static int DeterminePotency( unsigned int m, unsigned int a )
  {
    DoubleType d;
    unsigned int k;
    unsigned int b;
    int s;

    b = a-1;
    k = b;
    s = 1;
    while( (k!=0) && (s<100) )
      {
        DoubleMultiply(k,b,&d);
        k = DoubleModulus(&d,m);
        s++;
      }
    return s;
  }

  static int DetermineFactors( unsigned int x, unsigned int *factors )
  {
    unsigned int *ptr;
    unsigned int half;
    unsigned int i;

    half = x/2;
    ptr = factors;
    for( i=1; i<MAXPRIMES; i++ )
      {
        if( (x%primes[i]) == 0 )
          *ptr++ = primes[i];
        if( (unsigned int)(primes[i]) >= half )
          return ptr-factors;
      }

    for( i=primes[MAXPRIMES-1]+2; i<=half; i+=2 )
      if( IsOddPrime(i) && ((x%i)==0) )
        *ptr++ = i;
    return ptr-factors;
  }

  static unsigned int DetermineIncrement( unsigned int m )
  {
    unsigned int hi,lo;
    unsigned int half;
    unsigned int i;

    /* 1/2 + sqrt(3)/6 */
    hi = (int)floor(0.7886751345948*m+0.5);
    if( RelativelyPrime(m,hi) )
      return hi;

    /* 1/2 - sqrt(3)/6 */
    lo = (int)floor(0.2113248654052*m+0.5);
    if( RelativelyPrime(m,lo) )
      return lo;

    half = m/2;
    for( i=1; i<half; i++ )
      {
        if( RelativelyPrime(m,hi+i) )
          return hi+i;
        if( RelativelyPrime(m,hi-i) )
          return hi-i;
        if( RelativelyPrime(m,lo+i) )
          return lo+i;
        if( RelativelyPrime(m,lo-i) )
          return lo-i;
      }
    return 1;
  }

  static int DetermineSequence( unsigned int m, unsigned int *pm,
                         unsigned int *pa,
                         unsigned int *pc )
  {
    unsigned int fact[MAXFACT];
    unsigned int a=0, c;
    unsigned int b;
    int pot,best;
    int count;
    int flag;
    int i;

    do
      {
        best = 0;

        count = DetermineFactors(m,fact);
        if( (m&3) == 0 )
          fact[0] = 4;

        if( count )
          {
            for( b=m-2; b>0; b-- )
              {
                flag = True;
                for( i=0; i<count; i++ )
                  if( b%fact[i] )
                    {
                      flag = False;
                      break;
                    }

                if( flag )
                  {
                    pot = DeterminePotency(m,b+1);
                    if( pot > best )
                      {
                        best = pot;
                        a = b+1;
                      }
                  }
              }

          }

        m++;
      }
    while( best < 3 );
    m--;

    c = DetermineIncrement(m);

    *pm = m;
    *pa = a;
    *pc = c;
    return best;
  }


  static void GenerateSequence( unsigned int p, unsigned int m,
                         unsigned int a, unsigned int c )
  {
    unsigned int i;
    unsigned int x;
    DoubleType d;

    x = 0;  /* seed */
    for( i=0; i<p; i++ )
      {
        //        printf("%u\n",x);

        do
          {
            DoubleMultiply(a,x,&d);
            DoubleAdd(&d,c);
            x = DoubleModulus(&d,m);
          }
        while( x >= p );
      }
  }

  //**********************************************
  //***** Member functions from Random class *****
  //**********************************************

  OBRandom::OBRandom(bool useSysRand)
  {

    this->OBRandomUseSysRand= useSysRand;
    p = 70092;
    DetermineSequence(p,&m,&a,&c);
    x = 0;  /* seed */

    if (useSysRand)
      this->TimeSeed();
  }

  int OBRandom::NextInt()
  {
    if (OBRandomUseSysRand)
      {
        return(rand());
      }
    do
      {
        DoubleMultiply(a,x,&d);
        DoubleAdd(&d,c);
        x = DoubleModulus(&d,m);
      }
    while( x >= p );

    return(x);
  }

  double OBRandom::NextFloat()
  {

    if (OBRandomUseSysRand)
      {
        return(double(rand())/double(RAND_MAX));
      }
    do
      {
        DoubleMultiply(a,x,&d);
        DoubleAdd(&d,c);
        x = DoubleModulus(&d,m);
      }
    while( x >= p );

    return((double)x/p);
  }

  void OBRandom::TimeSeed()
  {
#ifdef WIN32
    // for VC++ do it this way
    time_t ltime;
    time(&ltime);
    const long secs= long(ltime);
    x= secs % p;
    srand( (unsigned)time( NULL ) );
#else

    timeval time;
    gettimeofday(&time,(struct timezone *)NULL);
    x = (time.tv_usec%p);
    srand( x );
#ifdef HAVE_SRANDDEV
    sranddev();
#endif

#endif
  }

} //end namespace OpenBabel

//! \file rand.cpp
//! \brief Pseudo random number generator.

/**********************************************************************
obutil.h - Various utility methods.
 
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2006 by Geoffrey R. Hutchison
 
This file is part of the Open Babel project.
For more information, see <http://openbabel.sourceforge.net/>
 
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#ifndef OB_UTIL_H
#define OB_UTIL_H

#include "babelconfig.h"

#include <string>
#include <iosfwd>

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

namespace OpenBabel
{

  // class introduction in obutil.cpp
  class OBAPI OBStopwatch
    {
#if HAVE_CLOCK_T
      clock_t start, stop;
#else

      timeval start;
      timeval stop;
#endif

    public:
#if HAVE_CLOCK_T

      void  Start()
        {
          start= clock();
        }
      double Lap()
        {
          stop= clock();
          return((double)(stop - start) / CLOCKS_PER_SEC);
        }
#else
      void Start()
        {
          gettimeofday(&start,(struct timezone *)NULL);
        }
      double Lap()
        {
          gettimeofday(&stop,(struct timezone *)NULL);
          return((stop.tv_sec - start.tv_sec)
                 + (double)(stop.tv_usec - start.tv_usec)/1000000.0);
        }
#endif
      double Elapsed()
        {
          return(Lap());
        }
    };


  //! sqrt lookup table - given a distance squared returns distance
  class OBAPI OBSqrtTbl
    {
      double _max,_incr,*_tbl;
    public:
      OBSqrtTbl()
        {
          _tbl=NULL;
          _max = _incr = 0.0;
        }
      OBSqrtTbl(double max,double incr)
        {
          Init(max,incr);
        }
      ~OBSqrtTbl()
        {
          if (_tbl)
            {
              delete [] _tbl;
              _tbl = NULL;
            }
        }
      double Sqrt(double d2) const
        {
          if (_tbl)
            return((d2 < _max) ? _tbl[(int)(d2*_incr)]:sqrt(d2));
          else
            return 0.0;
        }
      void Init(double max,double incr)
        {
          int i;
          double r;
          _max = max*max;
          _incr = incr;
          //array size needs to be large enough to account for fp error
          _tbl = new double [(unsigned int)((_max/_incr)+10)];
          for (r = (_incr/2.0),i=0;r <= _max;r += _incr,i++)
            _tbl[i] = sqrt(r);

          _incr = 1/_incr;
        }
    };



  //******************************************
  //*** Stuff for random number generation ***
  //******************************************

  //! Used for internal random number generation OBRandom (unless the system random generaor is used)
  typedef struct
  {
    unsigned int hi;
    unsigned int lo;
  }
  DoubleType;

  OBAPI void DoubleMultiply( unsigned int,unsigned int,DoubleType*);
  OBAPI void DoubleAdd( DoubleType*,unsigned int);
  OBAPI unsigned int DoubleModulus( DoubleType*,unsigned int);

  //! Random number generator
  class OBAPI OBRandom
    {
      DoubleType d;
      unsigned int m,a,c;
      unsigned int p;
      unsigned int i;
      unsigned int x;
      bool OBRandomUseSysRand;

    public:
      OBRandom(bool useSys= false);
      void Seed(int seed)
        {
          x = seed;
        }
      void TimeSeed();
      int NextInt();
      double NextFloat();
    };

  //***RMS helper methods***/
#ifndef SWIG
  OBAPI void  rotate_coords(double*,double m[3][3],int);
  OBAPI double calc_rms(double*,double*,unsigned int);

  //! \name  String conversion utilities
  //@{
  // Documentation in obutil.cpp
  OBAPI void ToUpper(std::string&);
  OBAPI void ToUpper(char*);
  OBAPI void ToLower(std::string&);
  OBAPI void ToLower(char *);
  //! "Clean" the supplied atom type
  OBAPI void CleanAtomType(char*);
  //@}

  //! Comparison -- returns true if first parameter less than second
  OBAPI bool OBCompareInt(const int &,const int &);
  //! Comparison -- returns true if first parameter less than second
  OBAPI bool OBCompareUnsigned(const unsigned int &,const unsigned int &);
  /*! Safe comparison for floats/doubles: returns fabs(a - b) < epsilon
   * For new code, unless you have a very specific need, you should
   * probably use IsApprox() instead, as it makes more sense w.r.t.
   * floating-point representation.
   */
  OBAPI bool IsNear(const double &, const double &, const double epsilon=2e-6);
  //! Safe comparison for floats/doubles: true if a is less than epsilon
  OBAPI bool IsNearZero(const double &, const double epsilon=2e-6);
  /*! New comparison for doubles: true if
   * fabs(a - b) <= precision * fmin( fabs(a), fabs(b) )
   * This is not a direct replacement for IsNear(). This is a different test.
   * The parameter precision plays the role of 10^-N where N is the number of
   * significant digits to consider.
   */
  OBAPI bool IsApprox(const double &, const double &, const double precision);

  /*! Same thing as IsApprox, but faster. Only works for nonnegative numbers.
   * No check is done.
   */
  OBAPI bool IsApprox_pos(const double &, const double &,
      const double precision);
  /*! Tests whether its argument can be squared without triggering an overflow
   * or underflow.
   */
  OBAPI bool CanBeSquared(const double &);

#endif
  // (end part to be skipped by SWIG)

  //******************triple template*************************
  //! \brief A 3-element templated, based on the design of the STL pair<>
  template <class T1, class T2, class T3>
    struct triple
    {
      //type names for the values
      typedef T1 first_type;
      typedef T2 second_type;
      typedef T3 third_type;

      //member
      T1 first;
      T2 second;
      T3 third;
  
      /** Default constructor
       *	T1() and T2() and T3() force initialization for built in types
       **/
      triple():
        first(T1()),second(T2()),third(T3())
      {}

      //! Constructor for 3 values
      triple(const T1 &a, const T2 &b, const T3 &c):
        first(a), second(b), third(c)
      {}

      //! Copy constructor with implicit conversions
      template<class U, class V, class W>
        triple(const triple<U,V,W> &t):
          first(t.first), second(t.second), third(t.third)
      {}

    };

  //**************quad template********************
  //! \brief A 4-element templated, based on the design of the STL pair<>
  template <class T1, class T2, class T3, class T4>
    struct quad
    {
      //type names for the values
      typedef T1 first_type;
      typedef T2 second_type;
      typedef T3 third_type;
      typedef T4 fourth_type;

      //member
      T1 first;
      T2 second;
      T3 third;
      T4 fourth;

      /*! default constructor
       *	T1() and T2() and T3() force initialization for built in types
       */
      quad():
        first(T1()),second(T2()),third(T3()),fourth(T4())
      {}

      //! constructor for 3 values
      quad(const T1 &a, const T2 &b, const T3 &c, const T4 &d):
        first(a), second(b), third(c), fourth(d)
      {}

      //! copy constructor with implicit conversions
      template<class U, class V, class W, class X>
        quad(const quad<U,V,W,X> &q):
          first(q.first), second(q.second), third(q.third), fourth(q.fourth)
      {}

    };

} // end namespace OpenBabel

#endif // OBUTIL_H

//! \file obutil.h
//! \brief Various utility methods.

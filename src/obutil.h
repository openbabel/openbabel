/**********************************************************************
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (c) 2001-2003 by Geoffrey R. Hutchison

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

#include <string>
#include <iostream>

#if defined(WIN32)
#include <time.h>
#else
#include <sys/time.h>
#endif

namespace OpenBabel {

  // class introduction in obutil.cpp
class OBStopwatch
{
#ifdef WIN32
  clock_t start, stop;
#else
  timeval start;
  timeval stop;
#endif

 public:
#ifdef WIN32
  void  Start() { start= clock();}
  double Lap()   
  { 
	  stop= clock();
	  return((double)(stop - start) / CLOCKS_PER_SEC);
  }
#else
  void Start() {gettimeofday(&start,(struct timezone *)NULL);}
  double Lap() 
    {
      gettimeofday(&stop,(struct timezone *)NULL);
      return((stop.tv_sec - start.tv_sec)
	     + (double)(stop.tv_usec - start.tv_usec)/1000000.0);
    }
#endif
  double Elapsed() {return(Lap());}
};


//! sqrt lookup table - given a distance squared returns distance
class OBSqrtTbl
{
  double _max,_incr,*_tbl;
 public:
  OBSqrtTbl() {_tbl=NULL; _max = _incr = 0.0;}
  OBSqrtTbl(double max,double incr) {Init(max,incr);}
  ~OBSqrtTbl() {if (_tbl) {delete [] _tbl; _tbl = NULL;}}
  double Sqrt(double d2) const
    {
      if (_tbl)
	return((d2 < _max) ? _tbl[(int)(d2*_incr)]:sqrt(d2));
      else
	return 0.0;
    }
  void Init(double max,double incr)
    {
      int i; double r;
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

typedef struct {
        unsigned int hi;
        unsigned int lo;
    } DoubleType;

void DoubleMultiply( unsigned int,unsigned int,DoubleType*);
void DoubleAdd( DoubleType*,unsigned int);
unsigned int DoubleModulus( DoubleType*,unsigned int);

//! Random number generator
class OBRandom
{
  DoubleType d;
  unsigned int m,a,c;
  unsigned int p;
  unsigned int i;
  unsigned int x;
  bool OBRandomUseSysRand;

 public:
    OBRandom(bool useSys= false);
    void Seed(int seed) {x = seed;}
    void TimeSeed();
    int NextInt();
    double NextFloat();
};

//***RMS helper methods***/
void  rotate_coords(double*,double m[3][3],int);
double calc_rms(double*,double*,int);

// String conversion utilities
void ToUpper(std::string&);
void ToUpper(char*);
void ToLower(std::string&);
void ToLower(char *);
void CleanAtomType(char*);

//! Comparison -- returns true if first parameter less than second
bool OBCompareInt(const int &,const int &);
//! Comparison -- returns true if first parameter less than second
bool OBCompareUnsigned(const unsigned int &,const unsigned int &);
//! Safe comparison for floats/doubles: true if a and b are closer than epsilon
bool IsNear(const double &, const double &, const double epsilon=2e-6);
//! Safe comparison for floats/doubles: true if a is less than epsilon
bool IsNearZero(const double &, const double epsilon=2e-6);

//******************triple template*************************
//based on the STL design of the pair<> template
template <class T1, class T2, class T3>
struct triple{
	//type names for the values
	typedef T1 first_type;
	typedef T2 second_type;
	typedef T3 third_type;

	//member
	T1 first;
	T2 second;
	T3 third;

	/*default constructor
	*	T1() and T2() and T3() force initialization for built in types
	*/
	triple()
		:	first(T1()),second(T2()),third(T3())
	{
	}

	//constructor for 3 values
	triple(const T1 &a, const T2 &b, const T3 &c)
		:	first(a), second(b), third(c)
	{
	}

	//copy constructor with implicit conversions
	template<class U, class V, class W>
	triple(const triple<U,V,W> &t)
		:	first(t.first), second(t.second), third(t.third)
	{
	}

};

//**************quad template********************
//based on the design of the STL pair<> template
template <class T1, class T2, class T3, class T4>
struct quad{
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

	/*default constructor
	*	T1() and T2() and T3() force initialization for built in types
	*/
	quad()
		:	first(T1()),second(T2()),third(T3()),fourth(T4())
	{
	}

	//constructor for 3 values
	quad(const T1 &a, const T2 &b, const T3 &c, const T4 &d)
		:	first(a), second(b), third(c), fourth(d)
	{
	}

	//copy constructor with implicit conversions
	template<class U, class V, class W, class X>
	quad(const quad<U,V,W,X> &q)
		:	first(q.first), second(q.second), third(q.third), fourth(q.fourth)
	{
	}

};

} // end namespace OpenBabel


#endif // OBUTIL_H

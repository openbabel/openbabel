/**********************************************************************
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#ifndef __OBUTIL_H__
#define __OBUTIL_H__

#include <string>

#ifdef __sgi
#include <iostream.h>
#else
#include <iostream>
#endif

#if defined(WIN32)
#include <time.h>
#else
#include <sys/time.h>
#endif

namespace OpenBabel {

//*** Stopwatch class used for timing length of execution ***
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
  float Lap()   
  { 
	  stop= clock();
	  return((float)(stop - start) / CLOCKS_PER_SEC);
  }
#else
  void Start() {gettimeofday(&start,(struct timezone *)NULL);}
  float Lap() 
    {
      gettimeofday(&stop,(struct timezone *)NULL);
      return((stop.tv_sec - start.tv_sec) + (float)
	(stop.tv_usec - start.tv_usec)/1000000.0);
    }
#endif
  float Elapsed() {return(Lap());}
};

//
//*** sqrt lookup table - given a distance squared returns distance
//
class OBSqrtTbl
{
  float _max,_incr,*_tbl;
 public:
  OBSqrtTbl() {_tbl=NULL;}
  OBSqrtTbl(float max,float incr) {Init(max,incr);}
  ~OBSqrtTbl() {if (_tbl) delete [] _tbl;}
  float Sqrt(float d2) {return((d2 < _max) ? _tbl[(int)(d2*_incr)]:float(sqrt(d2)));}
  void Init(float max,float incr)
    {
      int i; float r;
      _max = max*max;
      _incr = incr;
      //array size needs to be large enough to account for fp error
      _tbl = new float [(int)((_max/_incr)+10)];
      for (r = (_incr/2.0f),i=0;r <= _max;r += _incr,i++)
	_tbl[i] = float(sqrt(r));

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
    float NextFloat();
};

//***RMS helper methods***/
void rotate_coords(float*,float m[3][3],int);
float calc_rms(float*,float*,int);

//**File stream helper methods
bool SafeOpen(ifstream&,char*);
bool SafeOpen(ofstream&,char*);
bool SafeOpen(ifstream&,string&);
bool SafeOpen(ofstream&,string&);

void ToUpper(string&);
void ToUpper(char*);
void CleanAtomType(char*);

bool OBCompareInt(const int &,const int &);
bool OBCompareUnsigned(const unsigned int &,const unsigned int &);
void SmartsLexReplace(string &,vector<pair<string,string> > &);

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

}


#endif //__OBUTIL_H__

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

#ifndef __BITVEC_H
#define __BITVEC_H

#ifdef WIN32
#pragma warning (disable : 4786)
#endif

#ifdef __sgi 
#include <iostream.h>
#else
#include <iostream>
#endif

#include <algorithm>
#include <vector>
#include <string>

using namespace std;

#ifndef SETWORD
#define SETWORD 32
#endif

#ifndef STARTWORDS
#define STARTWORDS 10
#endif //STARTBITS

namespace OpenBabel {

class OBBitVec
{
  int _size;
  vector<int> _set;
public:
    OBBitVec() {_set.resize(STARTWORDS);_size=_set.size();Clear();}
    OBBitVec(int bits) 
      {_set.resize(bits/SETWORD);_size=_set.size();Clear(); }
    OBBitVec(const OBBitVec&);
    void SetBitOn(int);
    void SetBitOff(int);
    void SetRangeOn(int, int);
    void SetRangeOff(int, int);
    void Fold(int);

    int FirstBit(int) { return (BitIsSet(0) ? 0  : NextBit(0)); }
    int NextBit(int);
    int EndBit() {return(-1);}
    int GetSize() const {return(_size);}
    int CountBits();

    bool Empty() {return(IsEmpty());}
    bool IsEmpty();
    bool Resize(int);
    bool BitIsSet(int bit)
	{return((bit/SETWORD >= GetSize()) ? 
		false : _set[bit/SETWORD]>>(bit%SETWORD)&1);}
    bool BitIsOn(int bit)
	{return((bit/SETWORD >= GetSize()) ? 
		false : _set[bit/SETWORD]>>(bit%SETWORD)&1);}
    
    void FromVecInt(vector<int>&);
    void FromString(string&,int);
    void ToVecInt(vector<int>&);
    void Clear(void);
    void Negate() { for (int i= 0; i != _size; i++) { _set[i] = ~_set[i]; } }

    OBBitVec &operator= (const OBBitVec &);
    OBBitVec &operator&= (OBBitVec &);
    OBBitVec &operator|= (OBBitVec &);
    OBBitVec &operator|= (const int i) {SetBitOn(i);return(*this);}
    OBBitVec &operator^= (OBBitVec &);
    OBBitVec &operator-= (OBBitVec &);
    OBBitVec &operator+= (OBBitVec &bv);
    bool operator[] (int bit) 
	{return((bit/SETWORD >= GetSize()) ? 
		false : _set[bit/SETWORD]>>(bit%SETWORD)&1);}

    friend OBBitVec operator| (OBBitVec &, OBBitVec &);
    friend OBBitVec operator& (OBBitVec &,OBBitVec &);
    friend OBBitVec operator^ (OBBitVec &,OBBitVec &);
    friend OBBitVec operator- (OBBitVec &,OBBitVec &);
    friend bool operator== (const OBBitVec &,const OBBitVec &);

    friend istream& operator>> ( istream&, OBBitVec& );
    friend ostream& operator<< ( ostream&, const OBBitVec& ) ;
};

extern void ThrowError(char *);
float Tanimoto(OBBitVec&,OBBitVec&);

}

#endif //__BITVEC_H

/**********************************************************************
bitvec.cpp - Fast and efficient bitstring class.
 
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2005 by Geoffrey R. Hutchison
 
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

#include "bitvec.h"
#include "oberror.h"

using namespace std;

namespace OpenBabel
{

extern OBMessageHandler obErrorLog;

/*! \class OBBitVec
    \brief Fast and efficient bitstring class

The OBBitVec class is a fast and efficient bitstring class that is
handy to use as a truth table. Truth tables are an easy way to store
whether a list of items has a particular propery. Instances of
OBBitVec can by dynamically resized, and have a number of overloaded
operators that make code simple and readable. The following examples
demonstrate uses of the OBBitVec class:
\code
OBBitVec bv1,bv2,bv3;
bv1.SetBitOn(5);
bv2.SetBitOff(200);
bv1 |= bv2;
bv1 = bv1 & bv2;
if (bv1.Empty()) //Empty() returns true if no bits are set on
{
cout << "bv1 = " << bv1 << endl;
}

int bit;
for (bit = bv1.NextBit(0);bit != bv1.EndBit();bit = bv1.NextBit(bit))
{
cout << "the next bit turned on is " << bit << endl;
}
\endcode
 */

static int bitsoff[SETWORD] =
    {
        0xFFFFFFFF,0xFFFFFFFE,0xFFFFFFFC,0xFFFFFFF8,0xFFFFFFF0,0xFFFFFFE0,0xFFFFFFC0,
        0xFFFFFF80,0xFFFFFF00,0xFFFFFE00,0xFFFFFC00,0xFFFFF800,0xFFFFF000,0xFFFFE000,
        0xFFFFC000,0xFFFF8000,0xFFFF0000,0xFFFE0000,0xFFFC0000,0xFFF80000,0xFFF00000,
        0xFFE00000,0xFFC00000,0xFF800000,0xFF000000,0xFE000000,0xFC000000,0xF8000000,
        0xF0000000,0xE0000000,0xC0000000,0x80000000
    };

#ifndef LowBit
#define LowBit(set, bit)\
  {register int m;\
   if (set != 0)\
   {\
      bit = 31;\
      if (set != 0x80000000) {\
      if ((m = (set & 0x0000ffff))) {set = m; bit -= 16;}\
      if ((m = (set & 0x00ff00ff))) {set = m; bit -= 8;}\
      if ((m = (set & 0x0f0f0f0f))) {set = m; bit -= 4;}\
      if ((m = (set & 0x33333333))) {set = m; bit -= 2;}\
      if ((m = (set & 0x55555555))) {set = m; bit -= 1;}}}\
   else bit = -1;}
#endif

OBBitVec::OBBitVec(const OBBitVec &bv)
{
    _size = 0;
    (*this) = bv;
}

void OBBitVec::SetBitOn(int bit)
{
    int word = bit/SETWORD;
    bit = bit%SETWORD;

    //if (word+1 >= GetSize()) Resize((word+1)*SETWORD);
    if (word >= _size)
        Resize((word+1)*SETWORD);
    _set[word] |= (1<<bit);
}

void OBBitVec::SetRangeOn(int lobit, int hibit)
{
    int i;

    if (lobit > hibit)
        return;
    else if (lobit == hibit)
        SetBitOn(hibit);
    else
    {
        int loword = lobit / SETWORD;
        int hiword = hibit / SETWORD;
        int lobitp = lobit % SETWORD;
        int hibitp = hibit % SETWORD;

        if (hiword >= _size)
            Resize((hiword+1) * SETWORD);

        if (loword == hiword)
        {
            for ( i = lobitp ; i <= hibitp ; i++ )
                _set[loword] |= (1<<i);
        }
        else
        {
            for ( i = lobitp ; i < SETWORD ; i++ )
                _set[loword] |= (1<<i);
            for ( i = loword + 1 ; i < hiword ; i++ )
                _set[i] = 0xFFFFFFFF;
            for ( i = 0 ; i <= hibitp ; i++ )
                _set[hiword] |= (1<<i);
        }
    }
}

void OBBitVec::SetBitOff(int bit)
{
    int word = bit/SETWORD;
    bit = bit%SETWORD;
    _set[word] &= (~(1 << bit));
}

void OBBitVec::SetRangeOff(int lobit, int hibit)
{
    int i;

    if (lobit > hibit)
        return;
    else if (lobit == hibit)
        SetBitOff(hibit);
    else
    {
        int loword = lobit / SETWORD;
        int hiword = hibit / SETWORD;
        int lobitp = lobit % SETWORD;
        int hibitp = hibit % SETWORD;

        if (hiword >= _size)
        {
            hiword = _size - 1;
            hibitp = SETWORD - 1;
        }

        if (loword == hiword)
        {
            for ( i = lobitp ; i <= hibitp ; i++ )
                _set[loword] &= (~(1<<i));
        }
        else
        {
            for ( i = lobitp ; i < SETWORD ; i++ )
                _set[loword] &= (~(1<<i));
            for ( i = loword + 1 ; i < hiword ; i++ )
                _set[i] = 0x00000000;
            for ( i = 0 ; i <= hibitp ; i++ )
                _set[hiword] &= (~(1<<i));
        }
    }
}

int OBBitVec::NextBit(int last)
{
    unsigned s;
    register int bit,wrdcnt;
    last++;

    wrdcnt = last/SETWORD;

    if (wrdcnt >= _size)
        return(-1);

    if (_set[wrdcnt] != 0)
    {
        s = (_set[wrdcnt]) & bitsoff[last - (wrdcnt*SETWORD)];
        if (s)
        {
            LowBit(s,bit);
            if (bit != -1)
                return(bit + (wrdcnt*SETWORD));
        }
    }
    wrdcnt++;

    while(wrdcnt < _size)
    {
        if (this->_set[wrdcnt] != 0)
        {
            s = this->_set[wrdcnt];
            LowBit(s, bit);

            if (bit != -1)
                return(bit+(wrdcnt*SETWORD));
        }
        wrdcnt++;
    }

    return(-1);
}

bool OBBitVec::Resize(int maxbits)
{
    if(!maxbits)
        return(false);
    unsigned int maxword = maxbits/SETWORD;
    if (maxbits%SETWORD)
        maxword++;
    if (maxword >= _set.size())
    {
        _set.resize(maxword);
        _size = _set.size();
    }

    return(true);
}

OBBitVec &OBBitVec::operator= (const OBBitVec &bv)
{
    if (_size != bv._size)
        Resize(bv._size*SETWORD);

    int i;
    for ( i = 0 ; i < bv._size ; i++ )
        _set[i] = bv._set[i];
    for ( ; i < _size ; i++ )
        _set[i] = 0;

    return(*this);
}

OBBitVec &OBBitVec::operator&= ( OBBitVec &bv)
{
    int i;
    int min = (bv._size < _size) ? bv._size : _size;

    for (i = 0;i < min;i++)
        _set[i] &= bv._set[i];
    for (;i < _size;i++)
        _set[i] = 0;

    return(*this);
}

OBBitVec &OBBitVec::operator|= (OBBitVec &bv)
{
    if (_size != bv._size)
    {
        if (_size < bv._size)
            Resize(bv._size*SETWORD);
        else
            bv.Resize(_size*SETWORD);
    }
    for (int i = 0;i < _size;i++)
        _set[i] |= bv._set[i];

    return(*this);
}

OBBitVec &OBBitVec::operator^= (OBBitVec &bv)
{
    int i;
    if (_size != bv._size)
    {
        if (_size < bv._size)
            Resize(bv._size*SETWORD);
        else
            bv.Resize(_size*SETWORD);
    }
    for (i = 0;i < _size;i++)
        _set[i] ^= bv._set[i];

    return(*this);
}

OBBitVec &OBBitVec::operator-= (OBBitVec &bv)
{
    if (GetSize() != bv._size)
      obErrorLog.ThrowError(__FUNCTION__, "Subtracting sets of != size", obDebug);
    else
    {
        OBBitVec tmp;
        tmp = *this ^ bv;
        *this &= tmp;
    }
    return(*this);
}

OBBitVec &OBBitVec::operator+= (OBBitVec &bv)
{
    int old_size = _size;
    Resize(_size*SETWORD+bv._size*SETWORD);
    for (int i = 0;i < bv._size;i++)
        _set[i+old_size] = bv._set[i];
    return(*this);
}

OBBitVec operator| (OBBitVec &bv1,OBBitVec &bv2)
{
    OBBitVec bv;
    bv = bv1;
    bv |= bv2;

    return(bv);
}

OBBitVec operator& (OBBitVec &bv1, OBBitVec &bv2)
{
    OBBitVec bv;
    bv = bv1;
    bv &= bv2;

    return(bv);
}

OBBitVec operator^ (OBBitVec &bv1, OBBitVec &bv2)
{
    OBBitVec bv;
    bv = bv1;
    bv ^= bv2;
    return(bv);
}

bool operator== (const OBBitVec &bv1,const OBBitVec &bv2)
{
    if (bv1._size != bv2._size)
        return(false);

    int i;
    for (i = 0;i < bv1._size;i++)
        if (bv1._set[i] != bv2._set[i])
            return(false);

    return(true);
}

OBBitVec operator- (OBBitVec &bv1,OBBitVec &bv2)
{
    OBBitVec bv;
    bv = bv1 ^ bv2;
    bv &= bv1;
    return(bv);
}

bool operator< (const OBBitVec &bv1, const OBBitVec &bv2) 
{
    if (bv1._size > bv2._size)	return(false);

    int i;
    for (i = 0; i < bv1._size; i++)
      if (bv1._set[i] != (bv1._set[i] & bv2._set[i]))
	return(false);
    return(true);
}  

istream& operator>> ( istream &is, OBBitVec &bv )
{
    size_t startpos = 0, endpos = 0;
    vector<string> tokens;
    string line;

    getline(is,line);

    for (;;)
    {
        startpos = line.find_first_not_of(" \t\n",startpos);
        endpos   = line.find_first_of(" \t\n",startpos);

        if (endpos < line.size() && startpos <= line.size())
            tokens.push_back(line.substr(startpos,endpos-startpos));
        else
            break;

        startpos = endpos + 1;
    }

    for (unsigned int i = 0 ; i < tokens.size() ; i++ )
    {
        if ( tokens[i] == "[" )
            continue;
        else if ( tokens[i] == "]" )
            break;

        int bit = atoi(tokens[i].c_str());

        if (bit >= 0)
            bv.SetBitOn(bit);
        else
	  {
#ifdef HAVE_SSTREAM
		    stringstream errorMsg;
#else
		    strstream errorMsg;
#endif
            errorMsg << "Negative Bit: " << bit << endl;
	    obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obDebug);
	  }
    }

    return is;
}

ostream& operator<< ( ostream &os, const  OBBitVec& bv)
{
    os << "[ " << flush;

    int i,j;
    for (i = 0;i < bv._size;i++)
        for (j = 0;j < SETWORD;j++)
            if (bv._set[i]>>(j%SETWORD)&1)
                os << (j+(i*SETWORD)) << ' ' << flush;

    os << "]" << flush;
    return(os);
}

void OBBitVec::Fold(int nbits)
{
    int nwords = nbits/SETWORD;

    if (_size < nwords)
    {
        _set.resize(nwords);
        _size = nwords;
        return;
    }

    int i,idx = nwords;
    for (i = 0,idx=nwords;idx < _size;idx++)
    {
        _set[i] |= _set[idx];
        if (i+1 < nwords)
            i++;
        else
            i = 0;
    }
    _set.resize(nwords);
    _size = nwords;
}

void OBBitVec::FromVecInt(vector<int> &v)
{
    vector<int>::iterator i;
    int max = 0;

    for (i = v.begin();i != v.end();i++)
        if (*i > max)
            max = *i;

    Resize(max/SETWORD);
    for (i = v.begin();i != v.end();i++)
        SetBitOn(*i);
}

void OBBitVec::FromString(string &line, int bits)
{
    size_t startpos = 0, endpos = 0;
    vector<string> tokens;

    Resize(bits);
    Clear();

    for (;;)
    {
        startpos = line.find_first_not_of(" \t\n",startpos);
        endpos   = line.find_first_of(" \t\n",startpos);

        if (endpos < line.size() && startpos <= line.size())
            tokens.push_back(line.substr(startpos,endpos-startpos));
        else
            break;

        startpos = endpos + 1;
    }

    for (unsigned int i = 0 ; i < tokens.size() ; i++ )
    {
        if ( tokens[i] == "[" )
            continue;
        else if ( tokens[i] == "]" )
            break;

        int bit = atoi(tokens[i].c_str());

        if (bit >= 0)
            SetBitOn(bit);
        else
	  {
#ifdef HAVE_SSTREAM
	    stringstream errorMsg;
#else
	    strstream errorMsg;
#endif
            errorMsg << "Negative Bit: " << bit << endl;
	    obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obDebug);
	  }
    }
}

void OBBitVec::ToVecInt(vector<int> &v)
{
    v.clear();
    v.reserve(CountBits());
    for (int i = NextBit(-1);i != -1;i = NextBit(i))
        v.push_back(i);
}

void OBBitVec::Clear()
{
    vector<int>::iterator i;
    for (i = _set.begin();i != _set.end();i++)
        *i = 0;
}

int OBBitVec::CountBits()
{
    int i,count=0;
    for (i = NextBit(-1);i != -1;i = NextBit(i))
        count++;

    return(count);
}

bool OBBitVec::IsEmpty()
{
    vector<int>::iterator i;
    for (i = _set.begin();i != _set.end();i++)
        if (*i)
            return(false);

    return(true);
}

double Tanimoto(OBBitVec &bv1,OBBitVec &bv2)
{
    OBBitVec bvtmp;
    double andbits,orbits;

    bvtmp = bv1 & bv2;
    andbits = (double)bvtmp.CountBits();
    bvtmp = bv1 | bv2;
    orbits  = (double)bvtmp.CountBits();

    return(andbits/orbits);
}

} // end namespace OpenBabel

//! \file bitvec.cpp
//! \brief Fast and efficient bitstring class


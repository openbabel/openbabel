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

#include "bitvec.h"

using namespace std;

namespace OpenEye {

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

OEBitVec::OEBitVec(const OEBitVec &bv)
{
  _size = 0;
  (*this) = bv;
}

void OEBitVec::SetBitOn(int bit)
{
  int word = bit/SETWORD;
  bit = bit%SETWORD;

  //if (word+1 >= GetSize()) Resize((word+1)*SETWORD);
  if (word >= _size) Resize((word+1)*SETWORD);
  _set[word] |= (1<<bit);
}

void OEBitVec::SetRangeOn(int lobit, int hibit)
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

void OEBitVec::SetBitOff(int bit)
{
    int word = bit/SETWORD;
    bit = bit%SETWORD;
    _set[word] &= (~(1 << bit));
}

void OEBitVec::SetRangeOff(int lobit, int hibit)
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

int OEBitVec::NextBit(int last)
{
  unsigned s;
   register int bit,wrdcnt;
   last++;

   wrdcnt = last/SETWORD;
   
   if (wrdcnt >= _size) return(-1);

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

bool OEBitVec::Resize(int maxbits)
{
  if(!maxbits) return(false);
  unsigned int maxword = maxbits/SETWORD;
  if (maxbits%SETWORD) maxword++;
  if (maxword >= _set.size())
    {
      _set.resize(maxword);
      _size = _set.size();
    }

  return(true);
}

OEBitVec &OEBitVec::operator= (const OEBitVec &bv)
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

OEBitVec &OEBitVec::operator&= ( OEBitVec &bv)
{
  int i;
  int min = (bv._size < _size) ? bv._size : _size;

  for (i = 0;i < min;i++) _set[i] &= bv._set[i];
  for (;i < _size;i++) _set[i] = 0;

    return(*this);
}

OEBitVec &OEBitVec::operator|= (OEBitVec &bv)
{
  if (_size != bv._size) 
    {
      if (_size < bv._size)
	Resize(bv._size*SETWORD);
      else 
	bv.Resize(_size*SETWORD);
    }
  for (int i = 0;i < _size;i++)  _set[i] |= bv._set[i];
  
  return(*this);
}

OEBitVec &OEBitVec::operator^= (OEBitVec &bv)
{
  int i;
  if (_size != bv._size) 
    {
      if (_size < bv._size)
	Resize(bv._size*SETWORD);
      else 
	bv.Resize(_size*SETWORD);
    }
  for (i = 0;i < _size;i++) _set[i] ^= bv._set[i];

  return(*this);
}

OEBitVec &OEBitVec::operator-= (OEBitVec &bv)
{
  if (GetSize() != bv._size)
  ThrowError("Subtracting sets of != size");
  else
  {
    OEBitVec tmp;
    tmp = *this ^ bv;
    *this &= tmp;
  }
  return(*this);
}

OEBitVec &OEBitVec::operator+= (OEBitVec &bv)
{
  int old_size = _size;
  Resize(_size*SETWORD+bv._size*SETWORD);
  for (int i = 0;i < bv._size;i++)  _set[i+old_size] = bv._set[i];
  return(*this);
}

OEBitVec operator| (OEBitVec &bv1,OEBitVec &bv2)
{
    OEBitVec bv;
    bv = bv1;
    bv |= bv2;

    return(bv);
}

OEBitVec operator& (OEBitVec &bv1, OEBitVec &bv2)
{
  OEBitVec bv;
  bv = bv1;
  bv &= bv2;

  return(bv);
}

OEBitVec operator^ (OEBitVec &bv1, OEBitVec &bv2)
{
    OEBitVec bv;
    bv = bv1;
    bv ^= bv2;
    return(bv);
}

bool operator== (const OEBitVec &bv1,const OEBitVec &bv2)
{
    if (bv1._size != bv2._size)	return(false);

    int i;
    for (i = 0;i < bv1._size;i++)
	if (bv1._set[i] != bv2._set[i])
	    return(false);

    return(true);
}

OEBitVec operator- (OEBitVec &bv1,OEBitVec &bv2)
{
    OEBitVec bv;
    bv = bv1 ^ bv2;
    bv &= bv1;
    return(bv);
}

istream& operator>> ( istream &is, OEBitVec &bv )
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
	cerr << "Negative Bit: " << bit << endl;
    }

  return is;
}

ostream& operator<< ( ostream &os, const  OEBitVec& bv) 
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

void OEBitVec::Fold(int nbits)
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
      if (i+1 < nwords) i++;
      else i = 0;
    }
  _set.resize(nwords);
  _size = nwords;
}

void OEBitVec::FromVecInt(vector<int> &v)
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

void OEBitVec::FromString(string &line, int bits)
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
	cerr << "Negative Bit: " << bit << endl;
    }
}

void OEBitVec::ToVecInt(vector<int> &v)
{
    v.clear();
    v.reserve(CountBits());
    for (int i = NextBit(-1);i != -1;i = NextBit(i))
	v.push_back(i);
}

void OEBitVec::Clear()
{
  vector<int>::iterator i;
  for (i = _set.begin();i != _set.end();i++) *i = 0;
}

int OEBitVec::CountBits()
{
    int i,count=0;
    for (i = NextBit(-1);i != -1;i = NextBit(i))
	count++;

    return(count);
}

bool OEBitVec::IsEmpty()
{
    vector<int>::iterator i;
    for (i = _set.begin();i != _set.end();i++)
	if (*i)
	    return(false);
    
    return(true);
}

float Tanimoto(OEBitVec &bv1,OEBitVec &bv2)
{
  OEBitVec bvtmp;
  float andbits,orbits;

  bvtmp = bv1 & bv2;
  andbits = (float)bvtmp.CountBits();
  bvtmp = bv1 | bv2;
  orbits  = (float)bvtmp.CountBits();
  
  return(andbits/orbits);
}


}


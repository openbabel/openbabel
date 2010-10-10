/**********************************************************************
bitvec.cpp - Fast and efficient bitstring class.

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
#include <openbabel/babelconfig.h>

#include <openbabel/bitvec.h>
#include <openbabel/oberror.h>
#include <cstdlib>

namespace OpenBabel
{

  OBERROR extern OBMessageHandler obErrorLog;

  /*! \class OBBitVec bitvec.h <openbabel/bitvec.h>
    \brief Fast and efficient bitstring class

    The OBBitVec class is a fast and efficient bitstring class that is
    handy to use as a truth table. Truth tables are an easy way to store
    whether a list of items has a particular property. Instances of
    OBBitVec can be dynamically resized, and have a number of overloaded
    operators that make code simple and readable. The following examples
    demonstrate uses of the OBBitVec class:
    \code
    OBBitVec bv1,bv2,bv3;
    bv1.SetBitOn(5);
    bv2.SetBitOff(200);
    bv1 |= bv2;
    bv1 = bv1 & bv2;
    if (bv1.IsEmpty()) // IsEmpty() returns true if no bits are set on
    {
       std::cout << "bv1 = " << bv1 << std::endl;
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
#define LowBit(set, bit)                                        \
  {register int m;                                              \
    if (set != 0)                                               \
      {                                                         \
        bit = 31;                                               \
        if (set != 0x80000000) {                                \
          if ((m = (set & 0x0000ffff))!=0) {set = m; bit -= 16;}   \
          if ((m = (set & 0x00ff00ff))!=0) {set = m; bit -= 8;}    \
          if ((m = (set & 0x0f0f0f0f))!=0) {set = m; bit -= 4;}    \
          if ((m = (set & 0x33333333))!=0) {set = m; bit -= 2;}    \
          if ((m = (set & 0x55555555))!=0) {set = m; bit -= 1;}}}  \
    else bit = -1;}
#endif

  /** Set the \p bit_offset 'th bit to 1
    Increases the size of this bit vector if necessary
    \param[in] bit_offset a zero based offset into the bit vector
  */
  void OBBitVec::SetBitOn(unsigned bit_offset)
  {
    unsigned word_offset = bit_offset >> WORDROLL;
    bit_offset &= WORDMASK;

    if (word_offset >= GetSize())
      ResizeWords(word_offset + 1);
    _set[word_offset] |= (1<<bit_offset);
  }

  /** Set the \p bit_offset 'th bit to 0
    \param[in] bit_offset a zero based offset into the bit vector
  */
  void OBBitVec::SetBitOff(unsigned bit_offset)
  {
    unsigned word_offset = bit_offset >> WORDROLL;
    bit_offset &= WORDMASK;

    if (word_offset < GetSize())
      _set[word_offset] &= (~(1 << bit_offset));
  }

  /** Set the range of bits from \p lo_bit_offset to \p hi_bit_offset to 1
    Increases the size of this bit vector if necessary
    \param[in] lo_bit_offset a zero based offset into the bit vector
    \param[in] hi_bit_offset a zero based offset into the bit vector
  */
  void OBBitVec::SetRangeOn(unsigned lo_bit_offset, unsigned hi_bit_offset)
  {
    if (lo_bit_offset > hi_bit_offset)
      return;
    else if (lo_bit_offset == hi_bit_offset)
      SetBitOn(hi_bit_offset);
    else
      {
        unsigned lo_word_offset = lo_bit_offset >> WORDROLL;
        unsigned hi_word_offset = hi_bit_offset >> WORDROLL;
        lo_bit_offset &= WORDMASK;
        hi_bit_offset &= WORDMASK;

        if (hi_word_offset >= GetSize())
          ResizeWords(hi_word_offset + 1);

        if (lo_word_offset == hi_word_offset)
          {
            for ( unsigned i = lo_bit_offset ; i <= hi_bit_offset ; i++ )
              _set[lo_word_offset] |= (1<<i);
          }
        else
          {
            for ( unsigned i = lo_bit_offset ; i < SETWORD ; ++ i )
              _set[lo_word_offset] |= (1<<i);
            for ( unsigned i = lo_word_offset + 1 ; i < hi_word_offset ; ++ i )
              _set[i] = 0xFFFFFFFF;
            for ( unsigned i = 0 ; i <= hi_bit_offset ; ++ i )
              _set[hi_word_offset] |= (1<<i);
          }
      }
  }

  /** Set the range of bits from \p lo_bit_offset to \p hi_bit_offset to 0
    \param[in] lo_bit_offset a zero based offset into the bit vector
    \param[in] hi_bit_offset a zero based offset into the bit vector
  */
  void OBBitVec::SetRangeOff(unsigned lo_bit_offset, unsigned hi_bit_offset)
  {
    if (lo_bit_offset > hi_bit_offset)
      return;
    else if (lo_bit_offset == hi_bit_offset)
      SetBitOff(hi_bit_offset);
    else
      {
        unsigned lo_word_offset = lo_bit_offset >> WORDROLL;
        unsigned hi_word_offset = hi_bit_offset >> WORDROLL;
        lo_bit_offset &= WORDMASK;
        hi_bit_offset &= WORDMASK;

        if (lo_word_offset >= GetSize())
          return;
        if (hi_word_offset >= GetSize())
          {
            hi_word_offset = GetSize() - 1;
            hi_bit_offset = SETWORD - 1;
          }

        if (lo_word_offset == hi_word_offset)
          {
            for ( unsigned i = lo_bit_offset ; i <= hi_bit_offset ; ++ i )
              _set[lo_word_offset] &= (~(1<<i));
          }
        else
          {
            for ( unsigned i = lo_bit_offset ; i < SETWORD ; ++ i )
              _set[lo_word_offset] &= (~(1<<i));
            for ( unsigned i = lo_word_offset + 1 ; i < hi_word_offset ; ++ i )
              _set[i] = 0x00000000;
            for ( unsigned i = 0 ; i <= hi_bit_offset ; ++ i )
              _set[hi_word_offset] &= (~(1<<i));
          }
      }
  }

  /** Reduce the size of the vector to \p new_bit_size
  by or-ing the excess bits over the start of the vector
  \param[in] new_bit_size the size of the resultant vector, in bits
  */
  void OBBitVec::Fold(unsigned new_bit_size)
  {
    unsigned new_word_size = new_bit_size >> WORDROLL;

    if (_size < new_word_size)
      {
        ResizeWords(new_word_size);
        return;
      }

    for (unsigned i = 0, idx = new_word_size; idx < _size; ++idx )
      {
        _set[i] |= _set[idx];
        if (i+1 < new_word_size)
          ++i;
        else
          i = 0;
      }
    ResizeWords(new_word_size);
  }

  /** Searches the vector for the first true value, starting at the \p last_bit_offset 'th bit
      \param[in] last_bit_offset the bit before the first to consider
	    \return the bit offset of the first true bit after \p last_bit_offset, or -1 if there is none
  */
  int OBBitVec::NextBit(int last_bit_offset) const
  {
    unsigned s;
    int bit;
    unsigned wrdcnt;
    ++ last_bit_offset;

    wrdcnt = (unsigned)last_bit_offset >> WORDROLL;

    if (wrdcnt >= GetSize())
      return(-1);

    if (_set[wrdcnt] != 0)
      {
        s = _set[wrdcnt] & bitsoff[last_bit_offset & WORDMASK];
        if (s)
          {
            LowBit(s,bit);
            if (bit != -1)
              return(bit + (wrdcnt << WORDROLL));
          }
      }
    ++ wrdcnt;

    while(wrdcnt < GetSize())
      {
        if (_set[wrdcnt] != 0)
          {
            s = _set[wrdcnt];
            LowBit(s, bit);

            if (bit != -1)
              return(bit + (wrdcnt << WORDROLL));
          }
        ++ wrdcnt;
      }

    return(-1);
  }
  // Used by CountBits
  const unsigned nibble_bit_count[0x10] =
    {
      0, // 0000
      1, // 0001
      1, // 0010
      2, // 0011
      1, // 0100
      2, // 0101
      2, // 0110
      3, // 0111
      1, // 1000
      2, // 1001
      2, // 1010
      3, // 1011
      2, // 1100
      3, // 1101
      3, // 1110
      4  // 1111
    };
  /** Count the number of bits which are set in this vector
      \return the bit count
  */
  unsigned OBBitVec::CountBits() const
  {
    unsigned count = 0;
    for (word_vector::const_iterator sx = _set.begin(), sy = _set.end(); sx != sy; ++ sx)
      {
      unsigned word = * sx;
      while (word)
        {
        count += nibble_bit_count[word & 0xF];
        word >>= 4;
        }
      }
    return count;
  }

	/** Are there no bits set to 1 in this vector?
      \return true for "is empty", false if not empty
  */
  bool OBBitVec::IsEmpty() const
  {
    for (word_vector::const_iterator sx = _set.begin(), sy = _set.end(); sx != sy; ++ sx)
      if (* sx)
        return(false);

    return(true);
  }

  /** Sets bits on, listed as bit offsets
     \param[in] bit_offsets A list of bit offsets
  */
  void OBBitVec::FromVecInt(const std::vector<int> & bit_offsets)
  {
    for (std::vector<int>::const_iterator i = bit_offsets.begin(), j = bit_offsets.end(); i != j; ++i)
      SetBitOn(* i);
  }
  /** Sets bits on, listed as a string of character-represented integers
      This bit vector is first cleared.
      The format is "[ n0 n1 n2 n3 ... ]".
      The square brackets are optional.
      The whitespace can be SPACE, NEWLINE or HTAB
      For example "[ 1 5 6 9 ]"
     \param[in] line A string containing positive integers
     \param[in] new_bit_size The size that the vector should become
  */
  void OBBitVec::FromString(const std::string & line, int new_bit_size)
  {
    size_t startpos = 0, endpos = 0;
    std::vector<std::string> tokens;

    Clear();
    Resize(new_bit_size); // new bits are clear

    for (;;)
      {
        startpos = line.find_first_not_of(" \t\r\n",startpos);
        endpos   = line.find_first_of(" \t\r\n",startpos);

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
            std::stringstream errorMsg;
            errorMsg << "Negative Bit: " << bit << std::endl;
            obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obDebug);
          }
      }
  }

  /** Retrieve a list of bit offsets
     The \p bit_offsets vector is first cleared.
     \param[out] bit_offsets A list of bit offsets, in ascending order
  */
  void OBBitVec::ToVecInt(std::vector<int> & bit_offsets) const
  {
    bit_offsets.clear();
    bit_offsets.reserve(CountBits());
    for (int i = NextBit(-1);i != -1;i = NextBit(i))
      bit_offsets.push_back(i);
  }

  /** Set all the bits in this vector to zero
      Does not currently change the size of the vector.
  */
  void OBBitVec::Clear()
  {
    for (word_vector::iterator wx = _set.begin(), wy = _set.end(); wx != wy; ++wx)
      * wx = 0;
  }

  /** Assign this vector to be a copy of \p bv
      \param[in] bv A bit vector
      \return A reference to this
  */
  OBBitVec & OBBitVec::operator= (const OBBitVec & bv)
  {
    _set = bv._set;
    _size = _set.size();
    return(*this);
  }

  /** Assign this vector to the result of And-ing it with \p bv
      \param[in] bv A bit vector
      \return A reference to this
  */
  OBBitVec & OBBitVec::operator&= (const OBBitVec & bv)
  {
    unsigned min = (bv.GetSize() < _size) ? bv.GetSize() : _size;
    unsigned i;

    for (i = 0;i < min;++i)
      _set[i] &= bv._set[i];
    for (;i < _size;++i)
      _set[i] = 0;

    return(*this);
  }

  /** Assign this vector to the result of Or-ing it with \p bv
      \param[in] bv A bit vector
      \return A reference to this
  */
  OBBitVec & OBBitVec::operator|= (const OBBitVec & bv)
  {
    if (_size < bv.GetSize())
      ResizeWords(bv.GetSize());

    for (unsigned i = 0;i < bv.GetSize(); ++i)
      _set[i] |= bv._set[i];

    return(*this);
  }

  /** Assign this vector to the result of Exclusive-or-ing it with \p bv
      \param[in] bv A bit vector
      \return A reference to this
  */
  OBBitVec & OBBitVec::operator^= (const OBBitVec & bv)
  {
    if (_size < bv.GetSize())
      ResizeWords(bv.GetSize());

    for (unsigned i = 0;i < bv.GetSize(); ++i)
      _set[i] ^= bv._set[i];

    return(*this);
  }

  /** Unset bits in this vector which are set in \p bv
      \param[in] bv A bit vector
      \return A reference to this
  */
  OBBitVec & OBBitVec::operator-= (const OBBitVec & bv)
  {
    if (_size < bv.GetSize())
      ResizeWords(bv.GetSize());

    OBBitVec tmp(*this);
    tmp ^= bv;
    *this &= tmp;
    return(*this);
  }

  /** Append vector \p bv to the end if this vector
      \param[in] bv A bit vector
      \return A reference to this
  */
  OBBitVec & OBBitVec::operator+= (const OBBitVec & bv)
  {
    _set.insert(_set.end(), bv._set.begin(), bv._set.end());
    return(*this);
  }

  /** Return a bit vector of the results of Or-ing each bit in \p bv1 with the corresponding bit in \p bv2
      \param[in] bv1 A bit vector
      \param[in] bv2 Another bit vector
      \return A bit vector
  */
  OBBitVec operator| (const OBBitVec & bv1, const OBBitVec & bv2)
  {
    OBBitVec bv(bv1);
    bv |= bv2;
    return(bv);
  }

  /** Return a bit vector of the results of And-ing each bit in \p bv1 with the corresponding bit in \p bv2
      \param[in] bv1 A bit vector
      \param[in] bv2 Another bit vector
      \return A bit vector
  */
  OBBitVec operator& (const OBBitVec & bv1, const OBBitVec & bv2)
  {
    OBBitVec bv(bv1);
    bv &= bv2;
    return(bv);
  }

  /** Return a bit vector of the results of Exclusive-or-ing each bit in \p bv1 with the corresponding bit in \p bv2
      \param[in] bv1 A bit vector
      \param[in] bv2 Another bit vector
      \return A bit vector
  */
  OBBitVec operator^ (const OBBitVec & bv1, const OBBitVec & bv2)
  {
    OBBitVec bv(bv1);
    bv ^= bv2;
    return(bv);
  }

  /** Return a bit vector of the results of clearing each bit in \p bv1 which is set in \p bv2
      \param[in] bv1 A bit vector
      \param[in] bv2 Another bit vector
      \return A bit vector
  */
  OBBitVec operator- (const OBBitVec & bv1, const OBBitVec & bv2)
  {
    OBBitVec bv;
    bv = bv1 ^ bv2;
    bv &= bv1;
    return(bv);
  }

  /** Return true if \p bv1 and \p bv2 are equivalent
      Not that they may be of different size, and still equivalent provided that the extra bits are all zero.
      \param[in] bv1 A bit vector
      \param[in] bv2 Another bit vector
      \return true if equal, false otherwise
  */
  bool operator== (const OBBitVec & bv1, const OBBitVec & bv2)
  {
    if (bv1.GetSize() < bv2.GetSize())
      { // bv1 smaller than bv2
      unsigned i;
      for (i = 0; i < bv1.GetSize(); ++ i)
        if (bv1._set[i] != bv2._set[i])
          return false;
      for (; i < bv2.GetSize(); ++ i)
        if (bv2._set[i] != 0)
          return false;
      }
    else
      { // bv2 smaller or equal than bv1
      unsigned i;
      for (i = 0; i < bv2.GetSize(); ++ i)
        if (bv1._set[i] != bv2._set[i])
          return false;
      for (; i < bv1.GetSize(); ++ i)
        if (bv1._set[i] != 0)
          return false;
      }
    return true;
  }

  /** Return true if \p bv1 i less than \p bv2
      Lexicographical order, with bit vectors written LSB first.
      \param[in] bv1 A bit vector
      \param[in] bv2 Another bit vector
      \return true if equal, false otherwise
  */
  bool operator< (const OBBitVec & bv1, const OBBitVec & bv2)
  {
    bool rtn = false;
    int next_bit_1 = bv1.NextBit(-1);
    int next_bit_2 = bv2.NextBit(-1);
    bool should_continue = true;
    while (should_continue)
      {
      should_continue = false;
      if (next_bit_1 == -1)
        rtn = (next_bit_2 == -1 ? false : true);
      else if (next_bit_2 == -1)
        rtn = false;
      else if (next_bit_2 < next_bit_1)
        rtn = true;
      else if (next_bit_1 < next_bit_2)
        rtn = false;
      else
        {
        next_bit_1 = bv1.NextBit(next_bit_1);
        next_bit_2 = bv2.NextBit(next_bit_2);
        should_continue = true;
        }
      }
    return rtn;
  }

  /** Sets bits on, listed as a string of character-represented integers in a stream
      Only reads one line of input
      The format is "[ n0 n1 n2 n3 ... ]".
      The square brackets are optional.
      The whitespace can be SPACE or HTAB
      For example "[ 1 5 6 9 ]"
     \param[in,out] is The input stream
     \param[out] bv The bit vector to contain the result
  */
  std::istream & operator>> ( std::istream & is, OBBitVec & bv )
  {
    size_t startpos = 0, endpos = 0;
    std::vector<std::string> tokens;
    std::string line;

    getline(is,line);

    for (;;)
      {
        startpos = line.find_first_not_of(" \t\r\n",startpos);
        endpos   = line.find_first_of(" \t\r\n",startpos);

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
            std::stringstream errorMsg;
            errorMsg << "Negative Bit: " << bit << std::endl;
            obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obDebug);
          }
      }

    return is;
  }

  /** Output this bit vector to a stream
      The format is "[ n0 n1 n2 n3 ... ]".
      The whitespace is SPACE
      For example "[ 1 5 6 9 ]"
     \param[out] os The output stream
     \param[in] bv The bit vector to be output
  */
  std::ostream & operator<< ( std::ostream & os, const OBBitVec & bv)
  {
    os << "[ " << std::flush;

    for (unsigned i = 0;i < bv._size;++i)
      for (unsigned j = 0;j < SETWORD;++j)
        if (bv._set[i]>>(j%SETWORD)&1)
          os << (j+(i*SETWORD)) << ' ' << std::flush;

    os << "]" << std::flush;
    return(os);
  }

  /** The Tanimoto coefficient may be regarded as the proportion of the "on-bits" which are shared.
      \param[in] bv1 the first bit vector
      \param[in] bv2 the second bit vector
      \return the ratio of shared bits to bits which either vector has set.
  */
  double Tanimoto(const OBBitVec & bv1, const OBBitVec & bv2)
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


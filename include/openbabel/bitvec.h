/**********************************************************************
bitvec.h - Vector of bits.

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

#ifndef OB_BITVEC_H
#define OB_BITVEC_H

#include <openbabel/babelconfig.h>

#include <vector>
#include <string>

#if defined(_MSC_VER) && _MSC_VER <= 1600
  // Assuming 32bit integer
  typedef unsigned uint32_t;
#else
  #include <inttypes.h>
#endif

// Use uint32_t
#define SETWORD 32
// SETWORD = 2 ^ WORDROLL
#define WORDROLL 5
// WORDMASK = SETWORD - 1
#define WORDMASK 31

#define WORDSIZE_OF_BITSIZE( bit_size ) ( ( bit_size >> WORDROLL ) + (( bit_size & WORDMASK ) ? 1 : 0) )

#ifndef STARTWORDS
#define STARTWORDS 10
#endif // STARTWORDS

namespace OpenBabel
  {
  /// A speed-optimized vector of bits
  /** This class implements a fast vector of bits
      using internally a vector of fixed size unsigned ints (uint32_t).
      Any bits which are out of reach of the current size
      are considered to be zero.
      Streamlined, corrected and documented by kshepherd1@users.sourceforge.net
  */
  class OBERROR OBBitVec
    {
    public:
      typedef std::vector<uint32_t> word_vector;

	private:
	  /// The number of <b>words</b> currently stored ( NOT bit count )
      size_t _size; //was unsigned
	  /// A vector of words used to store the bit values
      word_vector	_set;

    public:
	  /// Construct a bit vector of the default size
	  /** Construct a bit vector of STARTWORDS size,
	      cleared to all zero bits.
	  */
      OBBitVec()
	  :_set(STARTWORDS, 0)
        { _size = _set.size(); }
	  /// Construct a bit vector of maxbits bits
	  /** Construct a bit vector with a size in bits
	      of \p size_in_bits rounded up to the nearest word
		  and cleared to all zero bits.
		  \param[in]	size_in_bits The number of bits for which to reserve space
	  */
      OBBitVec(unsigned size_in_bits)
	  :_set(WORDSIZE_OF_BITSIZE(size_in_bits), 0)
        { _size = _set.size(); }
      /// Copy constructor (result has same number of bits)
	  /** Construct a bit vector which is an exact
	      duplicate of \p bv.
		  \param[in]	bv The other bit vector to copy to this
	  */
      OBBitVec(const OBBitVec & bv)
	  :_size(0)
	  	{ (*this) = bv; }
	  /// Set the \p bit_offset 'th bit to 1
      void SetBitOn(unsigned bit_offset);
	  /// Set the \p bit_offset 'th bit to 0
      void SetBitOff(unsigned bit_offset);
	  /// Set the range of bits from \p lo_bit_offset to \p hi_bit_offset to 1
      void SetRangeOn(unsigned lo_bit_offset, unsigned hi_bit_offset);
	  /// Set the range of bits from \p lo_bit_offset to \p hi_bit_offset to 0
      void SetRangeOff(unsigned lo_bit_offset, unsigned hi_bit_offset);
	  /// Reduce the size of the vector by or-ing the excess bits over the start
      void Fold(unsigned new_bit_size);
      /// Find the first true bit at or after \p bit_offset
      /** Searches the vector for the first true value, starting at the \p bit_offset 'th bit
          \param[in] bit_offset the first bit to consider
		  \return the bit offset of the first true bit at or after \p bit_offset, or -1 if there is none
	  */
      int FirstBit(unsigned bit_offset = 0) const
        {
          return (BitIsSet(bit_offset) ? 0  : NextBit(bit_offset));
        }
      /// Find the next true bit after \p last_bit_offset
      int NextBit(int last_bit_offset) const;
      /// Return the bit offset of the last bit (for iterating) i.e. -1
      int EndBit() const {  return -1; }
      /// Return the number of words ( NOT the number of bits ).
      size_t GetSize() const    { return(_size);    }
      /// Return the number of bits which are set to 1 in the vector
      unsigned CountBits() const;

	  /// Are there no bits set to 1 in this vector?
      bool IsEmpty() const;
      /// Reserve space for \p size_in_bits bits
	  /** Reserve space for \p size_in_bits bits rounded up
	      \param[in] size_in_bits the number of bits
	      \return true if enlargement was necessary, false otherwise
	  */
      bool Resize(unsigned size_in_bits)
	  	{
		return ResizeWords( WORDSIZE_OF_BITSIZE(size_in_bits) );
		}
      /// Reserve space for \p size_in_words words
	  /** Reserve space for \p size_in_words words
	      \param[in] size_in_words the number of words
	      \return true if enlargement was necessary, false otherwise
	  */
	  bool ResizeWords(unsigned size_in_words)
	  	{
		if (size_in_words <= _size)
		  return false;
		_set.resize(size_in_words, 0); // increase the vector with zeroed bits
		_size = _set.size();
		return true;
		}
      /// Asks if the \p bit_offset 'th bit is set
	  /** Is the \p bit_offset 'th bit set ?
          \param[in] bit_offset a zero based offset into the bit vector
		  \return true if it is set, false otherwise
	  */
      bool BitIsSet(unsigned bit_offset) const
        {
		  bool rtn = false;
		  unsigned word_offset = bit_offset >> WORDROLL;
		  if (word_offset < GetSize())
		  	{
			  bit_offset &= WORDMASK;
			  rtn = (( _set[word_offset] >> bit_offset ) & 1);
			}
          return rtn;
        }
      /// Sets the bits listed as bit offsets
	  void FromVecInt(const std::vector<int> & bit_offsets);
      /// Sets the bits listed as a string of integers
	  void FromString(const std::string & line, int bits);
      /// List the offsets of the bits which are set
	  void ToVecInt(std::vector<int> & bit_offsets) const;
	  /// Set all bits to zero
      void Clear();
      /// Inverts every bit in the vector
	  /** Inverts the entire vector.
	      Note that this may give unexpected results, as the vector
		  can be considered to end in an arbitrary number of zero bits.
	  */
      void Negate()
        {
		  for (word_vector::iterator wx = _set.begin(), wy = _set.end(); wx != wy; ++wx)
		    * wx = ~(* wx);
        }
      /// Return a copy of the internal vector of words, at the end of \p vec
	  /** Copy the internal word vector.
	      The copy is appended to \p vec.
		  \param[out] vec a vector of words to which to append the data
	  */
      void GetWords(word_vector & vec)
        {
		vec.insert(vec.end(), _set.begin(),_set.end());
        }

      /// Assignment operator
      OBBitVec & operator= (const OBBitVec & bv);
      /// And-equals operator
      OBBitVec & operator&= (const OBBitVec & bv);
      /// Or-equals operator
      OBBitVec & operator|= (const OBBitVec & bv);
      /// Or-equals operator for integer
	  /** Or the bit at offset \p bit_offset with 1
	  */
      OBBitVec & operator|= (int bit_offset)
        {
          SetBitOn(bit_offset);
          return(*this);
        }
      /// Exclusive-or-equals operator
      OBBitVec & operator^= (const OBBitVec & bv);
      /// Minus-equals operator
      OBBitVec & operator-= (const OBBitVec & bv);
      /// Plus-equals operator
      OBBitVec & operator+= (const OBBitVec & bv);
      /// Asks if the \p bit_offset 'th bit is set
	  /** Is the \p bit_offset 'th bit set ?
          \param[in] bit_offset a zero based offset into the bit vector
		  \return true if it is set, false otherwise
	  */
      bool operator[] (int bit_offset) const
        { return BitIsSet(bit_offset); }

      /// Or operator
      friend OBERROR OBBitVec operator| (const OBBitVec & bv1, const OBBitVec & bv2);
      /// And operator
      friend OBERROR OBBitVec operator& (const OBBitVec & bv1,const OBBitVec & bv2);
      /// Exclusive-or operator
      friend OBERROR OBBitVec operator^ (const OBBitVec & bv1,const OBBitVec & bv2);
      /// Minus operator
      friend OBERROR OBBitVec operator- (const OBBitVec & bv1,const OBBitVec & bv2);
      /// Equivalency operator
      friend OBERROR bool operator== (const OBBitVec & bv1,const OBBitVec & bv2);
      /// Smaller-than operator
      friend OBERROR bool operator< (const OBBitVec & bv1, const OBBitVec & bv2);

      /// Input from a stream
      friend OBERROR std::istream& operator>> ( std::istream & is, OBBitVec & bv );
      /// Output to a stream
      friend OBERROR std::ostream& operator<< ( std::ostream & os, const OBBitVec & bv ) ;
    };

  /// The Tanimoto coefficient, which may be regarded as the proportion of the "on-bits" which are shared.
  OBERROR double Tanimoto(const OBBitVec & bv1, const OBBitVec & bv2);

  } // end namespace OpenBabel

#endif // OB_BITVEC_H

//! \file bitvec.h
//! \brief Fast and efficient bitstring class

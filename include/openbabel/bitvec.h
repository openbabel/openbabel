/**********************************************************************
bitvec.h - Vector of bits.
 
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

#ifndef OB_BITVEC_H
#define OB_BITVEC_H

#include <openbabel/babelconfig.h>

#ifdef WIN32
#pragma warning (disable : 4786)
#endif

#include <vector>
#include <string>

#ifndef SETWORD
#define SETWORD 32
#endif

#ifndef STARTWORDS
#define STARTWORDS 10
#endif //STARTBITS

namespace OpenBabel
{

  // class introduction in bitvec.cpp
  class OBAPI OBBitVec
    {
      int _size;
      std::vector<int> _set;
    public:
      OBBitVec()
        {
          _set.resize(STARTWORDS);
          _size=_set.size();
          Clear();
        }
      OBBitVec(int bits)
        {
          _set.resize(bits/SETWORD);
          _size=_set.size();
          Clear();
        }
      /// Copy constructor (result has same number of bits)
      OBBitVec(const OBBitVec&);
      void SetBitOn(int);
      void SetBitOff(int);
      void SetRangeOn(int, int);
      void SetRangeOff(int, int);
      void Fold(int);

      //! \return the index of the first bit past @p index that is set to true
      //! \param index the first bit to consider
      int FirstBit(int index = 0)
        {
          return (BitIsSet(index) ? 0  : NextBit(index));
        }
      int NextBit(int);
      //! \return the index of the last bit (for iterating)
      int EndBit()    {        return(-1);    }
      /// \return number of 32 bit words. NOT number of bits.
      int GetSize() const    { return(_size);    }
      /// \return the number of bits
      int CountBits();

      /// \deprecated Use IsEmpty() instead.
      bool Empty()   { return(IsEmpty()); }
      bool IsEmpty();
      ///Number of bits increased if necessary but never decreased
      bool Resize(int maxbits);

      bool BitIsSet(int bit)
        {
          return((bit/SETWORD >= GetSize()) ?
                 false : _set[bit/SETWORD]>>(bit%SETWORD)&1);
        }
      bool BitIsOn(int bit)
        {
          return((bit/SETWORD >= GetSize()) ?
                 false : _set[bit/SETWORD]>>(bit%SETWORD)&1);
        }

      void FromVecInt(std::vector<int>&);
      void FromString(std::string&,int);
      void ToVecInt(std::vector<int>&);
      void Clear(void);
      //! Inverts every bit in the vector
      void Negate()
        {
          for (int i= 0; i != _size; ++i)
            {
              _set[i] = ~_set[i];
            }
        }

      ///Assignment operator but number of bits is not reduced 
      OBBitVec &operator= (const OBBitVec &);
      OBBitVec &operator&= (OBBitVec &);
      OBBitVec &operator|= (OBBitVec &);
      OBBitVec &operator|= (const int i)
        {
          SetBitOn(i);
          return(*this);
        }
      OBBitVec &operator^= (OBBitVec &);
      OBBitVec &operator-= (OBBitVec &);
      OBBitVec &operator+= (OBBitVec &bv);
      bool operator[] (int bit)
        {
          return((bit/SETWORD >= GetSize()) ?
                 false : _set[bit/SETWORD]>>(bit%SETWORD)&1);
        }

      friend OBBitVec operator| (OBBitVec &, OBBitVec &);
      friend OBBitVec operator& (OBBitVec &,OBBitVec &);
      friend OBBitVec operator^ (OBBitVec &,OBBitVec &);
      friend OBBitVec operator- (OBBitVec &,OBBitVec &);
      friend bool operator== (const OBBitVec &,const OBBitVec &);
      friend bool operator< (const OBBitVec &bv1, const OBBitVec &bv2);

      friend std::istream& operator>> ( std::istream&, OBBitVec& );
      friend std::ostream& operator<< ( std::ostream&, const OBBitVec& ) ;
	
      ///Access to data in word size pieces CM
      void GetWords(std::vector<unsigned int>& vec)
        {
          std::vector<int>::iterator itr;
          for(itr=_set.begin();itr!=_set.end();itr++)
            vec.push_back(*itr);
        }
    };

  ///This function can change the size of second parameter. There is an alternative with different parameters.
  OBAPI double Tanimoto(OBBitVec&,OBBitVec&);

}

#endif // OB_BITVEC_H

//! \file bitvec.h
//! \brief Fast and efficient bitstring class

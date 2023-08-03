/**********************************************************************
fingerprint.h - Base class for fingerprints and fast searching

Copyright (C) 2005 by Chris Morley

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

#ifndef OB_FINGERPRINT_H
#define OB_FINGERPRINT_H

#include <list>
#include <map>
#include <set>
#include <vector>
#include <string>

#include <openbabel/plugin.h>

#ifndef OBFPRT
#define OBFPRT
#endif

namespace OpenBabel
{
  class OBBase; //Forward declaration; used only as pointer.

/// \brief The base class for fingerprints
class OBFPRT OBFingerprint : public OBPlugin
{
//see end of cpp file for detailed documentation

MAKE_PLUGIN(OBFingerprint)

const char* TypeID()
	{
		return "fingerprints";
	}

	//Rest of OBFingerprints declarations
public:

  virtual ~OBFingerprint(){}

  /// Sets the nth bit
  void SetBit(std::vector<unsigned int>& vec, const unsigned int n);

  ///return true if the nth bit is set;
  bool GetBit(const std::vector<unsigned int>& vec, const unsigned int n);

    /// Repeatedly ORs the top half with the bottom half until no smaller than nbits
  void Fold(std::vector<unsigned int>& vec, unsigned int nbits);

  /// \return fingerprint in vector, which may be resized, folded to nbits (if nbits!=0)
  virtual bool GetFingerprint(OBBase* pOb, std::vector<unsigned int>& fp, int nbits=0)=0;

  /// Optional flags
  enum FptFlag{FPT_UNIQUEBITS=1, FPT_NOINFO=2};
  virtual unsigned int Flags() { return 0;};
  //// \since version 2.3
  virtual void SetFlags(unsigned int){}

  /// \return a description of each bit that is set (or unset, if bSet=false)
  /// \since version 2.2
  virtual std::string DescribeBits(const std::vector<unsigned int> /* fp */,
                                   bool /* bSet */ =true)
  {
    std::string txt("");
    return txt;
  }

  /// \return the Tanimoto coefficient between two vectors (vector<unsigned int>& SeekPositions)
  static double Tanimoto(const std::vector<unsigned int>& vec1, const std::vector<unsigned int>& vec2);

  /// Inline version of Tanimoto() taking a pointer for the second vector
  static double Tanimoto(const std::vector<unsigned int>& vec1, const unsigned int* p2)
  {
    ///If used for two vectors, vec1 and vec2, call as Tanimoto(vec1, &vec2[0]);
    int andbits=0, orbits=0;
    unsigned int i;
    for (i=0;i<vec1.size();++i)
    {
      int andfp = vec1[i] & p2[i];
      int orfp = vec1[i] | p2[i];
      // Count bits
      /* GCC 3.4 supports a "population count" builtin, which on many targets is
         implemented with a single instruction.  There is a fallback definition
         in libgcc in case a target does not have one, which should be just as
         good as the static function below.  */
#if __GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 4)
      andbits += __builtin_popcount(andfp);
      orbits += __builtin_popcount(orfp);
#else
      for(;andfp;andfp=andfp<<1)
        if(andfp<0) ++andbits;
      for(;orfp;orfp=orfp<<1)
        if(orfp<0) ++orbits;
#endif
    }
      return((double)andbits/(double)orbits);
  };

  static unsigned int Getbitsperint(){ return bitsperint; }

private:
  ///Function object to set bits
  struct bit_or
  {
    unsigned int operator()(const unsigned int a, const unsigned int b)
    {
      return a | b;
    }
  };


public:
/// \return a pointer to a fingerprint (the default if ID is empty), or NULL if not available
  ///For backward compatibility;  a synonym of OBFingerprint::FindType
static OBFingerprint* FindFingerprint(const char* ID){ return FindType(ID);}

private:
  static const unsigned int bitsperint;// = 8 * sizeof(unsigned int);
};

//Fast search routines
/// \struct FptIndexHeader fingerprint.h <openbabel/fingerprint.h>
/// \brief Header for fastsearch index file
struct OBFPRT FptIndexHeader
{
  unsigned int headerlength;///<offset to data: sizeof(FptIndexHeader)
  unsigned int nEntries;    ///<number of fingerprints
  unsigned int words;				///<number 32bit words per fingerprint
  char fpid[15];            ///<ID of the fingerprint type
  char seek64; //if true, seek data consists of 64bit long values (only zero in legacy indices)
  char datafilename[256];   ///<the data that this is an index to
};

/// \struct FptIndex fingerprint.h <openbabel/fingerprint.h>
/// \brief Structure of fastsearch index files
struct OBFPRT FptIndex
{
  FptIndexHeader header;
  std::vector<unsigned int> fptdata;
  std::vector<unsigned long> seekdata;
  bool Read(std::istream* pIndexstream);
  bool ReadIndex(std::istream* pIndexstream);
  bool ReadHeader(std::istream* pIndexstream);

  /// \return A pointer to FP used or NULL and an error message
  OBFingerprint* CheckFP();
};

/// \class FastSearch fingerprint.h <openbabel/fingerprint.h>
/// \brief Class to search fingerprint index files
class OBFPRT FastSearch
{
//see end of cpp file for detailed documentation
public:
  /// \brief Loads an index from a file and returns the name of the datafile
  std::string ReadIndexFile(std::string IndexFilename);
  std::string ReadIndex(std::istream* pIndexstream);

  virtual ~FastSearch(){};

  /// \brief Does substructure search and returns vector of the file positions of matches
  bool    Find(OBBase* pOb, std::vector<unsigned long>& SeekPositions, unsigned int MaxCandidates);

  /// \brief Similar to Find() but all bits of matching fingerprints have to be the same
  /// \since version 2.1
  bool    FindMatch(OBBase* pOb, std::vector<unsigned long>& SeekPositions,
                            unsigned int MaxCandidates);

  /// \return A multimap containing objects whose Tanimoto coefficients with the target
  /// is greater than the value specified.
  bool    FindSimilar(OBBase* pOb, std::multimap<double, unsigned long>& SeekposMap,
    double MinTani, double MaxTani = 1.1 );

  /// \return A multimap containing the nCandidates objects with largest Tanimoto
  ///  coefficients with the target.
  bool    FindSimilar(OBBase* pOb, std::multimap<double, unsigned long>& SeekposMap,
    int nCandidates=0);

  /// \return a pointer to the fingerprint type used to constuct the index
  OBFingerprint* GetFingerprint() const{ return _pFP;};

  /// \return a pointer to the index header containing size info etc.
  const FptIndexHeader& GetIndexHeader() const{ return _index.header;};

private:
  FptIndex   _index;
  OBFingerprint* _pFP;
};

/// \class FastSearchIndexer fingerprint.h <openbabel/fingerprint.h>
/// \brief Class to prepare fingerprint index files See FastSearch class for details
class OBFPRT FastSearchIndexer
{
//see end of cpp file for detailed documentation
public:
  ///\brief Constructor with a new index
  FastSearchIndexer(std::string& datafilename, std::ostream* os, std::string& fpid,
      int FptBits=0, int nmols=0);

  ///\brief Constructor using existing index
  FastSearchIndexer(FptIndex* pindex, std::ostream* os, int nmols=0);

  ~FastSearchIndexer();

  ///\brief Called for each object
  bool Add(OBBase* pOb, std::streampos seekpos);

private:
  std::ostream* _indexstream;
  FptIndex*		_pindex;
  OBFingerprint* _pFP;
  int _nbits;
};

} //namespace OpenBabel
#endif

//! \file fingerprint.h
//! \brief Declaration of OBFingerprint base class and fastsearch classes

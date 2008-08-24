/**********************************************************************
fingerprint.h - Base class for fingerprints and fast searching 
 
Copyright (C) 2005 by Chris Morley
 
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

  /** @brief The base class for fingerprints
   */
  class OBFPRT OBFingerprint : public OBPlugin
  {
    //see end of cpp file for detailed documentation
    MAKE_PLUGIN(OBFingerprint)

    const char* TypeID() { return "fingerprints"; }
    
    //Rest of OBFingerprints declarations
    public:
      /** 
       * @brief Destructor.
       */
      virtual ~OBFingerprint(){}
      /** 
       * @brief Set the n-th bit.
       * 
       * @param vec The fingerprint vector.
       * @param n The bit number to set.
       */
      void SetBit(std::vector<unsigned int>& vec, const unsigned int n);
      /** 
       * @brief Get the n-th bit
       * 
       * @param vec The fingerprint vector.
       * @param n The bit number.
       * 
       * @return The n-th bit.
       */
      bool GetBit(const std::vector<unsigned int>& vec, const unsigned int n);
      /** 
       * @brief Repeatedly ORs the top half with the bottom half until 
       * no smaller than nbits.
       * 
       * @param vec The fingerprint vector.
       * @param nbits Number of bits.
       */
      void Fold(std::vector<unsigned int>& vec, unsigned int nbits); 
      /** 
       * @brief Get the fingerprint @p fp for @p pOb
       * 
       * @param pOb The object to get the fingerprint for (molecule).
       * @param fp The fingerprint vector.
       * @param nbits Fold size.
       * 
       * @return Fingerprint in vector, which may be resized, folded to nbits (if nbits!=0)
       */
      virtual bool GetFingerprint(OBBase* pOb, std::vector<unsigned int>& fp, int nbits=0) = 0;
      //! Optional flags
      enum FptFlag {
        FPT_UNIQUEBITS = 1
      };
      /** 
       * @brief Get the flags.
       */
      virtual unsigned int Flags() { return 0; } 
      /**
       * @since version 2.2
       *
       * @param fp The fingerprint vector..
       * @param bSet True = set bits, false = unset bits.
       * 
       * @return A description of each bit that is set (or unset, if bSet=false).
       */
      virtual std::string DescribeBits(const std::  vector<unsigned int> fp, bool bSet=true)
      {
        std::string txt("");
        return txt;
      }
      /** 
       * @brief Calculate Tanimoto coefficient between two vectors 
       * 
       * @param vec1 First vector.
       * @param vec2 Second vector.
       * 
       * @return The Tanimoto coefficient between two vectors (vector<unsigned int>& SeekPositions)
       */
      static double Tanimoto(const std::vector<unsigned int>& vec1, const std::vector<unsigned int>& vec2);
      /** 
       * @brief Inline version of Tanimoto() taking a pointer for the second vector
       * 
       * @param vec1 First vector.
       * @param p2 Second vector.
       * 
       * @return The Tanimoto coefficient between two vectors (vector<unsigned int>& SeekPositions) 
       */
      static double Tanimoto(const std::vector<unsigned int>& vec1, const unsigned int* p2) 
      {
        /// If used for two vectors, vec1 and vec2, call as Tanimoto(vec1, &vec2[0]);
        int andbits=0, orbits=0;
        unsigned int i;
        for (i=0;i<vec1.size();++i)
        {
          int andfp = vec1[i] & p2[i];
          int orfp = vec1[i] | p2[i];
          // Count bits
          #ifdef __GNUC__
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
      /** 
       * @return The number of bits per int. 
       */
      static unsigned int Getbitsperint(){ return bitsperint; }

    private:
      /** 
       * @brief Function object to set bits
       */
      struct bit_or
      {
        unsigned int operator()(const unsigned int a, const unsigned int b)
        {
          return a | b;	
        }
      };
  
    public:
      /** 
       * @brief For backward compatibility;  a synonym of OBFingerprint::FindType.
       * 
       * @param ID The fingerprint ID.
       * 
       * @return A pointer to a fingerprint (the default if ID is empty), or NULL if not available
       */
      static OBFingerprint* FindFingerprint(const char* ID) { return FindType(ID); }

    private:
      static const unsigned int bitsperint; //!< = 8 * sizeof(unsigned int);
  };

  /** @struct FptIndexHeader fingerprint.h <openbabel/fingerprint.h>
      @brief Header for fastsearch index file
   */
  struct OBFPRT FptIndexHeader
  {
    unsigned int headerlength; //!< offset to data: sizeof(FptIndexHeader)
    unsigned int nEntries;     //!< number of fingerprints
    unsigned int words;	       //!< number 32bit words per fingerprint
    char fpid[16];             //!< ID of the fingerprint type
    char datafilename[256];    //!< the data that this is an index to
  };

  /** @struct FptIndex fingerprint.h <openbabel/fingerprint.h>
      @brief Structure of fastsearch index files
   */
  struct OBFPRT FptIndex
  {
    FptIndexHeader header;
    std::vector<unsigned int> fptdata;
    std::vector<unsigned int> seekdata;
    bool Read(std::istream* pIndexstream);
    /** 
     * @return A pointer to FP used or NULL and an error message.
     */
    OBFingerprint* CheckFP();
  };

  /** @class FastSearch fingerprint.h <openbabel/fingerprint.h>
      @brief Class to search fingerprint index files
   */
  class OBFPRT FastSearch
  {
    //see end of cpp file for detailed documentation
    public:
      /** 
       * @brief Loads an index from a file and returns the name of the datafile
       * 
       * @param IndexFilename The filename containing the index.
       * 
       * @return The name of the datafile.
       */
      std::string ReadIndexFile(std::string IndexFilename);
      /** 
       * @brief Loads an index from a std::istream and returns the name of the datafile
       * 
       * @param pIndexStream The std::istream containing the index.
       * 
       * @return The name of the datafile.
       */
 
      std::string ReadIndex(std::istream* pIndexstream);
      /** 
       * @brief Destructor. 
       */
      virtual ~FastSearch(){};
      /** 
       * @brief Does substructure search and returns vector of the file positions of matches
       * 
       * @param pOb 
       * @param SeekPositions 
       * @param MaxCandidates Maximim number of matches.
       * 
       * @return 
       */
      bool Find(OBBase* pOb, std::vector<unsigned int>& SeekPositions, unsigned int MaxCandidates);
      /** 
       * @brief Similar to Find() but all bits of matching fingerprints have to be the same.
       * @since version 2.1
       * 
       * @param pOb
       * @param SeekPositions
       * @param MaxCandidates Maximim number of matches.
       * 
       * @return 
       */
      bool FindMatch(OBBase* pOb, std::vector<unsigned int>& SeekPositions,
          unsigned int MaxCandidates);
      /** 
       * @param pOb
       * @param SeekposMap
       * @param MinTani Minimum Tanimoto coefficient.
       * 
       * @return A multimap containing objects whose Tanimoto coefficients with 
       * the target is greater than the value specified.
       */
      bool FindSimilar(OBBase* pOb, std::multimap<double, unsigned int>& SeekposMap,
          double MinTani);
      /** 
       * @param pOb
       * @param SeekposMap
       * @param nCandidates
       * 
       * @return  A multimap containing the nCandidates objects with largest 
       * Tanimoto coefficients with the target.
       */
      bool FindSimilar(OBBase* pOb, std::multimap<double, unsigned int>& SeekposMap,
          int nCandidates=0);
      /** 
       * @return A pointer to the fingerprint type used to constuct the index.
       */
      OBFingerprint* GetFingerprint() const { return m_pFP; }
      /** 
       * @return A pointer to the index header containing size info etc.
       */
      const FptIndexHeader& GetIndexHeader() const { return m_index.header; }
    
    private:
      FptIndex       m_index; //!< the index
      OBFingerprint* m_pFP;   //!< the fingerprint
  };

  /** @class FastSearchIndexer fingerprint.h <openbabel/fingerprint.h>
      @brief Class to prepare fingerprint index files See FastSearch class for details
   */
  class OBFPRT FastSearchIndexer
  {
    //see end of cpp file for detailed documentation
    public:
      /** 
       * @brief Constructor with a new index.
       * 
       * @param datafilename The datafile.
       * @param os The index stream.
       * @param fpid The fingerprint ID
       * @param FptBits 
       */
      FastSearchIndexer(std::string& datafilename, std::ostream* os, std::string& fpid,
          int FptBits = 0);
      /** 
       * @brief Constructor using existing index.
       * 
       * @param pindex Pointer to the index.
       * @param os The index stream.
       */
      FastSearchIndexer(FptIndex* pindex, std::ostream* os);
      /** 
       * @brief Constructor.
       */
      ~FastSearchIndexer();
      /** 
       * @brief  Called for each object
       * 
       * @param pOb
       * @param seekpos
       * 
       * @return 
       */
      bool Add(OBBase* pOb, std::streampos seekpos);

    private:
      std::ostream*  m_indexstream; //!< the index stream
      FptIndex*	     m_pindex;      //!< the index
      OBFingerprint* m_pFP;         //!< the fingerprint
      int            m_nbits;       //!< number of bits
  };

} //namespace OpenBabel
#endif

//! @file fingerprint.h
//! @brief Declaration of OBFingerprint base class and fastsearch classes

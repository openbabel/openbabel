/**********************************************************************
fingerprint.cpp - Implementation of fingerpring base class and fastsearching

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

#include <openbabel/babelconfig.h>

#include <vector>
#include <algorithm>
#include <iosfwd>
#include <cstring>
#include <fstream>

#include <openbabel/fingerprint.h>
#include <openbabel/oberror.h>

using namespace std;
namespace OpenBabel
{
#if defined(__CYGWIN__) || defined(__MINGW32__)
  // macro to implement static OBPlugin::PluginMapType& Map()
  PLUGIN_CPP_FILE(OBFingerprint)
#endif

  const unsigned int OBFingerprint::bitsperint = 8 * sizeof(unsigned int);

  void OBFingerprint::SetBit(vector<unsigned int>& vec, const unsigned int n)
  {
    vec[n/Getbitsperint()] |= (1 << (n % Getbitsperint()));
  }

  bool OBFingerprint::GetBit(const vector<unsigned int>& vec, const unsigned int n)
  {
    unsigned int word =vec[n/Getbitsperint()];
    return (word &= (1 << (n % Getbitsperint())))!=0;
  }

  ////////////////////////////////////////
  void OBFingerprint::Fold(vector<unsigned int>& vec, unsigned int nbits)
  {
    if(nbits<Getbitsperint())
    {
      stringstream ss;
      ss << "Can't fold to less than " << Getbitsperint() << "bits";
      obErrorLog.ThrowError(__FUNCTION__, ss.str(), obError);
      return;
    }
    // "folding" to a larger # of bits
    if (nbits > vec.size()*Getbitsperint()) {
      vec.resize(nbits/Getbitsperint(), 0);
    }
    else {
      // normal folding to smaller vector sizes
      while(vec.size()*Getbitsperint()/2 >= nbits)
        vec.erase(transform(vec.begin(),vec.begin()+vec.size()/2,
                            vec.begin()+vec.size()/2, vec.begin(), bit_or()), vec.end());
    }
  }

  ////////////////////////////////////////
/*  bool OBFingerprint::GetNextFPrt(std::string& id, OBFingerprint*& pFPrt)
  {
    Fptpos iter;
    if(id.empty())
      iter=FPtsMap().begin();
    else
      {
        iter=FPtsMap().find(id);
        if(iter!=FPtsMap().end())
          ++iter;
      }
    if(iter==FPtsMap().end())
      return false;
    id    = iter->first;
    pFPrt = iter->second;
    return true;
  }

  OBFingerprint* OBFingerprint::FindFingerprint(const string& ID)
  {
    if(ID.empty())
      return _pDefault;
    Fptpos iter = FPtsMap().find(ID);
    if(iter==FPtsMap().end())
      return NULL;
    else
      return iter->second;
  }
*/
  double OBFingerprint::Tanimoto(const vector<unsigned int>& vec1, const vector<unsigned int>& vec2)
  {
    //Independent of sizeof(unsigned int)
    if(vec1.size()!=vec2.size())
      return -1; //different number of bits
    int andbits=0, orbits=0;
    for (unsigned i=0;i<vec1.size();++i)
      {
        int andfp = vec1[i] & vec2[i];
        int orfp = vec1[i] | vec2[i];
        //Count bits
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
    if(orbits==0)
      return 0.0;
    return((double)andbits/(double)orbits);
  }

  //*****************************************************************
  bool FastSearch::Find(OBBase* pOb, vector<unsigned long>& SeekPositions,
                        unsigned int MaxCandidates)
  {
    ///Finds chemical objects in datafilename (which must previously have been indexed)
    ///that have all the same bits set in their fingerprints as in the fingerprint of
    ///a pattern object. (Usually a substructure search in molecules.)
    ///The type of fingerprint and its degree of folding does not have to be specified
    ///here because the values in the index file are used.
    ///The positions of the candidate matching molecules in the original datafile are returned.

    vector<unsigned int> vecwords;
    _pFP->GetFingerprint(pOb,vecwords, _index.header.words * OBFingerprint::Getbitsperint());

    vector<unsigned int>candidates; //indices of matches from fingerprint screen
    candidates.reserve(MaxCandidates);

    unsigned int dataSize = _index.header.nEntries;
    //	GetFingerprint(mol, vecwords, _index.header.words, _index.header.fptype);

    unsigned int words = _index.header.words;
    unsigned int* nextp = &_index.fptdata[0];
    unsigned int* ppat0 = &vecwords[0];
    unsigned int* p;
    unsigned int* ppat;
    unsigned int i;
    for(i=0;i<dataSize; ++i) //speed critical section
      {
        p=nextp;
        nextp += words;
        ppat=ppat0;
        bool ppat_has_additional_bits = false;
        while(p<nextp)
          {
            if ((*ppat & *p) ^ *ppat) { // any bits in ppat that are not in p?
              ppat_has_additional_bits = true;
              break;
            }
            p++;
            ppat++;
          }
        if(!ppat_has_additional_bits)
          {
            candidates.push_back(i);
            if(candidates.size()>=MaxCandidates)
              break;
          }
      }

    if(i<_index.header.nEntries) //premature end to search
      {
        stringstream errorMsg;
        errorMsg << "Stopped looking after " << i << " molecules." << endl;
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
      }

    vector<unsigned int>::iterator itr;
    for(itr=candidates.begin();itr!=candidates.end();++itr)
      {
        SeekPositions.push_back(_index.seekdata[*itr]);
      }
    return true;
  }

////////////////////////////////////////////////////////////
 bool FastSearch::FindMatch(OBBase* pOb, vector<unsigned long>& SeekPositions,
                            unsigned int MaxCandidates)
{
//Similar to FastSearch::Find() except that successful candidates have all bits the same as the target
  vector<unsigned int> vecwords;
  _pFP->GetFingerprint(pOb,vecwords, _index.header.words * OBFingerprint::Getbitsperint());

  vector<unsigned int>candidates; //indices of matches from fingerprint screen

  unsigned int dataSize = _index.header.nEntries;
  unsigned int words = _index.header.words;
  unsigned int* nextp = &_index.fptdata[0]; // start of next FP in index
  unsigned int* ppat0 = &vecwords[0];       // start of target FP
  unsigned int* p;                          // current position in index
  unsigned int* ppat;                       // current position in target FP
  unsigned int i; // need address of this, can't be register
  for(i=0;i<dataSize; ++i) //speed critical section
  {
    p=nextp;
    nextp += words;
    ppat=ppat0;

    while((*p++ == *ppat++ ) && (p<nextp));

    if(p==nextp)
    {
      candidates.push_back(i);
      if(candidates.size()>=MaxCandidates)
        break;
    }
  }

  vector<unsigned int>::iterator itr;
  for(itr=candidates.begin();itr!=candidates.end();++itr)
    {
      SeekPositions.push_back(_index.seekdata[*itr]);
    }
  return true;
}

  /////////////////////////////////////////////////////////
  bool FastSearch::FindSimilar(OBBase* pOb, multimap<double, unsigned long>& SeekposMap,
                               double MinTani, double MaxTani)
  {
    vector<unsigned int> targetfp;
    _pFP->GetFingerprint(pOb,targetfp, _index.header.words * OBFingerprint::Getbitsperint());

    unsigned int words = _index.header.words;
    unsigned int dataSize = _index.header.nEntries;
    unsigned int* nextp = &_index.fptdata[0];
    unsigned int* p;
    unsigned int i;
    for(i=0;i<dataSize; ++i) //speed critical section
      {
        p=nextp;
        nextp += words;
        double tani = OBFingerprint::Tanimoto(targetfp,p);
        if(tani>MinTani && tani < MaxTani)
          SeekposMap.insert(pair<const double, unsigned long>(tani,_index.seekdata[i]));
      }
    return true;
  }

  /////////////////////////////////////////////////////////
  bool FastSearch::FindSimilar(OBBase* pOb, multimap<double, unsigned long>& SeekposMap,
                               int nCandidates)
  {
    ///If nCandidates is zero or omitted the original size of the multimap is used
    if(nCandidates)
      {
        //initialise the multimap with nCandidate zero entries
        SeekposMap.clear();
        int i;
        for(i=0;i<nCandidates;++i)
          SeekposMap.insert(pair<const double, unsigned long>(0,0));
      }
    else if(SeekposMap.size()==0)
      return false;

    vector<unsigned int> targetfp;
    _pFP->GetFingerprint(pOb,targetfp, _index.header.words * OBFingerprint::Getbitsperint());

    unsigned int words = _index.header.words;
    unsigned int dataSize = _index.header.nEntries;
    unsigned int* nextp = &_index.fptdata[0];
    unsigned int* p;
    unsigned int i;
    for(i=0;i<dataSize; ++i) //speed critical section
      {
        p=nextp;
        nextp += words;
        double tani = OBFingerprint::Tanimoto(targetfp,p);
        if(tani>SeekposMap.begin()->first)
          {
            SeekposMap.insert(pair<const double, unsigned long>(tani,_index.seekdata[i]));
            SeekposMap.erase(SeekposMap.begin());
          }
      }
    return true;
  }

  /////////////////////////////////////////////////////////
  string FastSearch::ReadIndex(istream* pIndexstream)
  {
    //Reads fs index from istream into member variables
    _index.Read(pIndexstream);

    _pFP = _index.CheckFP();
    if(!_pFP)
      *(_index.header.datafilename) = '\0';

    return _index.header.datafilename; //will be empty on error
  }

  //////////////////////////////////////////////////////////
  string FastSearch::ReadIndexFile(string IndexFilename)
  {
    ifstream ifs(IndexFilename.c_str(),ios::binary);
    if(ifs)
      return ReadIndex(&ifs);
    else
    {
      string dum;
      return dum;
    }
  }

  //////////////////////////////////////////////////////////
  bool FptIndex::Read(istream* pIndexstream)
  {
//    pIndexstream->read((char*)&(header), sizeof(FptIndexHeader));
//    pIndexstream->seekg(header.headerlength);//allows header length to be changed

    if(!ReadHeader(pIndexstream))
      {
        *(header.datafilename) = '\0';
        return false;
      }

    unsigned long nwords = header.nEntries * header.words;
    fptdata.resize(nwords);
    seekdata.resize(header.nEntries);

    pIndexstream->read((char*)&(fptdata[0]), sizeof(unsigned int) * nwords);
    if(header.seek64) 
      {
    	pIndexstream->read((char*)&(seekdata[0]), sizeof(unsigned long) * header.nEntries);
      }
    else 
      { //legacy format
	 vector<unsigned int> tmp(header.nEntries);
         pIndexstream->read((char*)&(tmp[0]), sizeof(unsigned int) * header.nEntries);
	 std::copy(tmp.begin(),tmp.end(),seekdata.begin());
      }

    if(pIndexstream->fail())
      {
        *(header.datafilename) = '\0';
        return false;
      }
    return true;
  }

  //////////////////////////////////////////////////////////
  bool FptIndex::ReadHeader(istream* pIndexstream)
  {
    pIndexstream->read( (char*)&header.headerlength, sizeof(unsigned) );
    pIndexstream->read( (char*)&header.nEntries,     sizeof(unsigned) );
    pIndexstream->read( (char*)&header.words,        sizeof(unsigned) );
    pIndexstream->read( (char*)&header.fpid,         sizeof(header.fpid) );
    pIndexstream->read( (char*)&header.seek64,       sizeof(header.seek64) );
    pIndexstream->read( (char*)&header.datafilename, sizeof(header.datafilename) );
    return !pIndexstream->fail();
 }

  //////////////////////////////////////////////////////////
  OBFingerprint* FptIndex::CheckFP()
  {
    //check that fingerprint type is available
    OBFingerprint* pFP = OBFingerprint::FindFingerprint(header.fpid);
    if(!pFP)
      {
        stringstream errorMsg;
        errorMsg << "Index has Fingerprints of type '" << header.fpid
                 << " which is not currently loaded." << endl;
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
      }
    return pFP; //NULL if not available
  }

  //*******************************************************
  FastSearchIndexer::FastSearchIndexer(string& datafilename, ostream* os,
                                       std::string& fpid, int FptBits, int nmols)
  {
    ///Starts indexing process
    _indexstream = os;
    _nbits=FptBits;
    _pindex= new FptIndex;
    _pindex->header.headerlength = 3*sizeof(unsigned)+sizeof(_pindex->header.fpid)
                                    +sizeof(_pindex->header.datafilename);
    strncpy(_pindex->header.fpid,fpid.c_str(),15);
    _pindex->header.fpid[14]='\0'; //ensure fpid is terminated at 14 characters.
    _pindex->header.seek64 = 1;
    strncpy(_pindex->header.datafilename, datafilename.c_str(), 255);

    //just a hint to reserve size of vectors; definitive value set in destructor
    _pindex->header.nEntries = nmols;

    //check that fingerprint type is available
    _pFP = _pindex->CheckFP();
    if(fpid.empty()) // add id of default FP
      strcpy(_pindex->header.fpid, _pFP->GetID());

    //Save a small amount of time by not generating info (FP2 currently)
    _pFP->SetFlags(_pFP->Flags() | OBFingerprint::FPT_NOINFO);
  }

  /////////////////////////////////////////////////////////////
  FastSearchIndexer::FastSearchIndexer(FptIndex* pindex, std::ostream* os, int nmols)
  {
    //nmols is new total number of molecules
    _indexstream = os;
    _pindex = pindex;
    _nbits  = _pindex->header.words * OBFingerprint::Getbitsperint();

    //just a hint to reserve size of vectors; definitive value set in destructor
    _pindex->header.nEntries = nmols;

    //check that fingerprint type is available
    _pFP = _pindex->CheckFP();
  }

  /////////////////////////////////////////////////////////////
  FastSearchIndexer::~FastSearchIndexer()
  {
    ///Saves index file
    FptIndexHeader& hdr = _pindex->header;
    hdr.nEntries = _pindex->seekdata.size();
    //Write header
    //_indexstream->write((const char*)&hdr, sizeof(FptIndexHeader));
    _indexstream->write( (const char*)&hdr.headerlength, sizeof(unsigned) );
    _indexstream->write( (const char*)&hdr.nEntries,     sizeof(unsigned) );
    _indexstream->write( (const char*)&hdr.words,        sizeof(unsigned) );
    _indexstream->write( (const char*)&hdr.fpid,         sizeof(hdr.fpid) );
    _indexstream->write( (const char*)&hdr.seek64,         sizeof(hdr.seek64) );
    _indexstream->write( (const char*)&hdr.datafilename, sizeof(hdr.datafilename) );

    _indexstream->write((const char*)&_pindex->fptdata[0], _pindex->fptdata.size()*sizeof(unsigned int));
    _indexstream->write((const char*)&_pindex->seekdata[0], _pindex->seekdata.size()*sizeof(unsigned long));
    if(!_indexstream)
      obErrorLog.ThrowError(__FUNCTION__,
                            "Difficulty writing index", obWarning);
    delete _pindex;

    _pFP->SetFlags(_pFP->Flags() & ~OBFingerprint::FPT_NOINFO); //Clear
  }

  ///////////////////////////////////////////////////////////////
  bool FastSearchIndexer::Add(OBBase* pOb, std::streampos seekpos)
  {
    ///Adds a fingerprint

    vector<unsigned int> vecwords;
    if(!_pFP)
      return false;
    if(_pFP->GetFingerprint(pOb, vecwords, _nbits))
      {
        _pindex->header.words = vecwords.size(); //Use size as returned from fingerprint
        if(_pindex->fptdata.empty() && _pindex->header.nEntries!=0)
        {
          //Reserve size of vectors at start to avoid multiple realloction and copying later.
          //Done here rather than in constructor because needs the size of the fingerprint.
          _pindex->fptdata.reserve(_pindex->header.nEntries * _pindex->header.words);
          _pindex->seekdata.reserve(_pindex->header.nEntries);
        }
        for(unsigned int i=0;i<_pindex->header.words;++i)
          _pindex->fptdata.push_back(vecwords[i]);
        _pindex->seekdata.push_back(seekpos);
        return true;
      }
    obErrorLog.ThrowError(__FUNCTION__, "Failed to make a fingerprint", obWarning);
    return false;
  }

  /*!
    \class OBFingerprint fingerprint.h <openbabel/fingerprint.h>
    These fingerprints are condensed representation of molecules (or other objects)
    as a list of boolean values (actually bits in a vector<unsigned>) with length
    of a power of 2. The main motivation is for fast searching of data sources
    containing large numbers of molecules (up to several million). Open Babel
    provides some routines which can search text files containing lists of molecules
    in any format. See the documentation on the class FastSearch.

    There are descriptions of molecular fingerprints at <br>
    http://www.daylight.com/dayhtml/doc/theory/theory.finger.html) and <br>
    http://www.mesaac.com/Fingerprint.htm <br>
    Many methods of preparing fingerprints have been described, but the type supported
    currently in OpenBabel has each bit representing a substructure (or other
    molecular property). If a substructure is present in the molecule, then a
    particular bit is set to 1. But because the hashing method may also map other
    substructures to the same bit, a match does not guarantee that a particular
    substructure is present; there may be false positives.  However, with proper design,
    a large fraction of irrelevant molecules in a data set can be eliminated in a
    fast search with boolean methods on the fingerprints.
    It then becomes feasible to make a definitive substructure search by
    conventional methods on this reduced list even if it is slow.

    OpenBabel provides a framework for applying new types of fingerprints without
    changing any existing code. They are derived from OBFingerprint and the
    source file is just compiled with the rest of OpenBabel. Alternatively,
    they can be separately compiled as a DLL or shared library and discovered
    when OpenBabel runs.

    For more on these specific implementations of fingerprints in Open
    Babel, please take a look at the developer's wiki:
    http://openbabel.org/wiki/Fingerprints

    Fingerprints derived from this abstract base class OBFingerprint can be for any
    object derived from OBBase (not just for OBMol).
    Each derived class provides an ID as a string and OBFingerprint keeps a map of
    these to provides a pointer to the class when requested in FindFingerprint.

    <h4>-- To define a fingerprint type --</h4>
    The classes derived form OBFingerprint are required to provide
    a GetFingerprint() routine and a Description() routine
    \code
    class MyFpType : OBFingerprint
    {
       MyFpType(const char* id) : OBFingerprint(id){};

       virtual bool GetFingerprint(OBBase* pOb, vector<unsigned int>& fp, int nbits)
       {
          //Convert pOb to the required type, usually OBMol
          OBMol* pmol = dynamic_cast<OBMol*>(pOb);
          fp.resize(required_number_of_words);
          ...
          use SetBit(fp,n); to set the nth bit

          if(nbits)
             Fold(fp, nbits);
       }

       virtual const char* Description(){ return "Some descriptive text";}
       ...
    };
    \endcode

    Declare a global instance with the ID you will use in -f options to specify
    its use.
    \code
    MyFpType theMyFpType("myfpID");
    \endcode

    <h4>-- To obtain a fingerprint --</h4>
    \code
    OBMol mol;
    ...
    vector<unsigned int> fp;
    OBFingerprint::GetDefault()->GetFingerprint(&mol, fp); //gets default size of fingerprint
    \endcode
    or
    \code
    vector<unsigned int> fp;
    OBFingerPrint* pFP = OBFingerprint::FindFingerprint("myfpID");
    ...and maybe...
    pFP->GetFingerprint(&mol,fp, 128); //fold down to 128bits if was originally larger
    \endcode

    <h4>-- To print a list of available fingerprint types --</h4>
    \code
    std::string id;
    OBFingerPrint* pFPrt=NULL;
    while(OBFingerprint::GetNextFPrt(id, pFPrt))
    {
       cout << id << " -- " << pFPrt->Description() << endl;
    }
    \endcode

    Fingerprints are handled as vector<unsigned int> so that the number of bits
    in this vector and their order will be platform and compiler
    dependent, because of size of int types and endian differences.
    Use fingerprints (and fastsearch indexes containing them) only
    for comparing with other fingerprints prepared on the same machine.

    The FingerprintFormat class is an output format which displays fingerprints
    as hexadecimal. When multiple molecules are supplied it will calculate the
    Tanimoto coefficient from the first molecule to each of the others. It also
    shows whether the first molecule is a possible substructure to all the others,
    i.e. whether all the bits set in the fingerprint for the first molecule are
    set in the fingerprint of the others. To display hexadecimal information when
    multiple molecules are provided it is necessay to use the -xh option.

    To see a list of available format types, type obabel -F on the command line.
    The -xF option of the FingerprintFormat class also provides this output, but due
    to a quirk in the way the program works, it is necessary to have a valid input
    molecule for this option to work.
  */

  /*! \class FastSearch fingerprint.h <openbabel/fingerprint.h>
    The FastSearch class searches an index to a datafile containing a list of molecules
    (or other objects) prepared by FastSearchIndexer.

    OpenBabel can also search files for molecules containing a substructure specified
    by a SMARTS string, using OBSmartsPattern or from the command line:
    \code
    obabel datafile.xxx -O outfile.yyy -sSMARTS
    \endcode
    But when the data file contains more than about 10,000 molecules this becomes
    rather too slow to be used interactively. To do it more quickly, an index
    of the molecules containing their fingerprints (see OBFingerprint) is prepared using
    FastSearchIndexer. The indexing process may take a long time but only has to
    be done once. The index can be searched very quickly with FastSearch. Because
    of the nature of fingerprints a match to a bit does not guarantee
    the presence of a particular substructure or other molecular property, so that
    a definitive answer may require a subsequent search of the (much reduced) list
    of candidates by another method (like OBSmartsPattern).

    Note that the index files are not portable. They need to be prepared on the
    computer that will access them.

    <h4>Using FastSearch and FastSearchIndexer in a program</h4>

    The index has two tables:
    - an array of fingerprints of the molecules
    - an array of the seek positions in the datasource file of all the molecules

    <h4>To prepare an fastsearch index file:</h4>
    - Open an ostream to the index file.
    - Make a FastSearchIndexer object on the heap or the stack, passing in as parameters:
    - the datafile name, the indexfile stream,
    - the id of the fingerprint type to be used,
    -  the number of bits the fingerprint is to be folded down to, If it is to be left
    unfolded, set fpid to 0 or do not specify it.
    .
    - For each molecule, call Add() with its pointer and its position in the datafile.<br>
    Currently the std::streampos value is implicitly cast to unsigned int so that
    for 32bit machines the datafile has to be no longer than about 2Gbyte.
    - The index file is written when the FastSearchIndexer object is deleted or goes out
    of scope.

    <h4>To search in a fastsearch index file</h4>

    - Open a std::istream to the indexfile (in binary mode on some systems)
    - Make a FastSearch object, read the index and open the datafile whose
    name it provides
    \code
    ifstream ifs(indexname,ios::binary);
    FastSearch fs;
    string datafilename = fs.ReadIndex(&ifs);
    if(datafilename.empty()
       return false;

    ifstream datastream(datafilename);
    if(!datastream)
       return false;
    \endcode

    <strong>To do a search for molecules which have all the substructure bits the
    OBMol object, patternMol</strong>
    \code
    vector<unsigned int>& SeekPositions;
    if(!fs.Find(patternMol, SeekPositions, MaxCandidates))
	    for(itr=SeekPositions.begin();itr!=SeekPositions.end();++itr)
      {
         datastream.seekg(*itr);
         ... read the candidate molecule
         and subject to more rigorous test if necessary
      }
    \endcode

    <strong>To do a similarity search based on the Tanimoto coefficient</strong>
    This is defined as: <br>
    <em>Number of bits set in (patternFP & targetFP)/Number of bits in (patternFP | targetFP)</em><br>
    where the boolean operations between the fingerprints are bitwise<br>
    The Tanimoto coefficient has no absolute meaning and depends on
    the design of the fingerprint.
    \code
    multimap<double, unsigned int> SeekposMap;
    // to find n best molecules
    fs.FindSimilar(&patternMol, SeekposMap, n);
    ...or
    // to find molecules with Tanimoto coefficient > MinTani
    fs.FindSimilar(&patternMol, SeekposMap, MinTani);

    multimap<double, unsigned int>::reverse_iterator itr;
    for(itr=SeekposMap.rbegin();itr!=SeekposMap.rend();++itr)
    {
       datastream.seekg(itr->second);
       // ... read the candidate molecule
       double tani = itr->first;
    }
    \endcode

    The FastSearchFormat class facilitates the use of these routine from the
    command line or other front end program. For instance:

    <strong>Prepare an index:</strong>
    \code
    obabel datafile.xxx -O index.fs
    \endcode
    With options you can specify:
    - which type of fingerprint to be used, e.g. -xfFP2,
    -	whether it is folded to a specified number of bits, e.g. -xn128
    (which should be a power of 2)
    - whether to pre-select the molecules which are indexed:
    - by structure e.g only ethers and esters, -sCOC
    - by excluding molecules with bezene rings, -vc1ccccc1
    - by position in the datafile e.g. only mols 10 to 90, -f10 -l90
    .
    <strong>Substructure search</strong> in a previously prepared index file
    \code
    obabel index.fs -O outfile.yyy -sSMILES
    \endcode
    The datafile can also be used as the input file, provided the input format is specified as fs
    \code
    obabel datafile.xxx -O outfile.yyy -sSMILES -ifs
    \endcode
    A "slow" search not using fingerprints would be done on the same data by omitting -ifs.
    A "slow" search can use SMARTS strings, but the fast search is restricted
    to the subset which is valid SMILES.

    With the -S option, the target can be specified as a molecule in a file of any format
    \code
    obabel datafile.xxx -O outfile.yyy -Spattern_mol.zzz -ifs
    \endcode
    These searches have two stages: a fingerprint search which produces
    a number of candidate molecules and a definitive search which selects
    from these using SMARTS pattern matching. The maximum number of candidates
    is 4000 by default but you can change this with an option
    e.g. -al 8000  (Note that you need the space between l and the number.)
    If the fingerprint search reaches the maximum number it will not
    look further and will tell you at what stage it stopped.

    <strong>Similarity search</strong> in a previously prepared index file<br>
    This rather done (rather than a substructure search) if the -at option is used,
    \code
    obabel datafile.xxx -O outfile.yyy -sSMILES -ifs -at0.7
    \endcode
    for instance
    - -at0.7 will recover all molecules with Tanimoto greater than 0.7
    - -at15 (no decimal point) will recover the 15 molecules with largest coefficients.
    - -aa will add the Tanimoto coefficient to the titles of the output molecules.

    All stages, the indexing, the interpretation of the SMILES string in the -s option,
    the file in the -S option and the final SMARTS filter convert to OBMol and apply
    ConvertDativeBonds(). This ensures thatforms such as[N+]([O-])=O  and N(=O)=O
    are recognized as chemically identical.
  */

}//Openbabel

//! \file fingerprint.cpp
//! \brief Definitions for OBFingerprint base class and fastsearch classes

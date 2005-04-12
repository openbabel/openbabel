#ifndef FAST_SEARCH_H
#define FAST_SEARCH_H

#ifdef WIN32
#pragma warning (disable : 4786)
#endif

#include "mol.h"
#include "fingerprint.h"

#ifndef MAXCANDIDATES
#define MAXCANDIDATES 100
#endif

namespace OpenBabel {

struct FptIndexHeader
{
	unsigned int headerlength;//offset to data: sizeof(FptIndexHeader)
	unsigned int nEntries;    //number of fingerprints
	unsigned int words;				//number 32bit words per fingerprint
	int fptype;               //type of fingerprint
	char datafilename[256];   //the data that this is an index to
};
struct FptIndex
{
	FptIndexHeader header;
	std::vector<unsigned int> fptdata;
	std::vector<unsigned int> seekdata;
};

class FastSearch
{
	friend class FastSearchIndexer;
public:
	virtual      ~FastSearch();
	bool         ReadIndex(std::istream* pIndexstream, FptIndex* pindex);
	std::string  Find(OBMol& mol, std::istream* pIndexstream, 
                  	std::vector<unsigned int>& SeekPositions, int MaxCandidates);
};

//**********************************************
///Class to make fingerprint index files
class FastSearchIndexer
{
public:
	FastSearchIndexer(std::string datafilename, std::ostream* os, int FptType,
										int FptWords, unsigned int NfptEstimate=100);
	~FastSearchIndexer();
	bool Add(OBMol& mol, std::streampos seekpos);

private:
	std::ostream* _indexstream;
	FptIndex*		_pindex;
	int _fptype;
};

} //OpenBabel

#endif
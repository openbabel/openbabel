#ifndef FAST_SEARCH_H
#define FAST_SEARCH_H
/**********************************************************************
Copyright (C) 2005 by Chris Morley
 
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

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

class OBAPI FastSearch
{
public:
  std::string ReadIndex(std::istream* pIndexstream);
	virtual    ~FastSearch();
	bool       Find(OBMol& mol, std::vector<unsigned int>& SeekPositions, int MaxCandidates);
private:
	FptIndex   _index;
};

//**********************************************
///Class to make fingerprint index files
class OBAPI FastSearchIndexer
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
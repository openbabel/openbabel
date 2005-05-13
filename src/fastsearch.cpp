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

#if HAVE_CONFIG_H
#include "config.h"
#endif
#include <fstream>
#include <vector>
#include "mol.h"
#include "fingerprint.h"
#include "fastsearch.h"

using namespace std;
namespace OpenBabel {

FastSearch::~FastSearch(){}

//////////////////////////////////////////////////
bool FastSearch::Find(OBMol& mol, vector<unsigned int>& SeekPositions,
												int MaxCandidates)
{
	//Finds molecules in datafilename (which must previously have been indexed)
	//that match a substructure in mol. datafile name can also be the name of the index file.
	//The positions of the candidate matching molecules in the original datafile are returned.
	//When called, datafilename can be either the real data file or the index.
	
	vector<unsigned int>candidates; //indices of matches from fingerprint screen
	candidates.reserve(MaxCandidates);

	unsigned int dataSize = _index.header.nEntries;
	vector<unsigned int> vecwords;
	GetFingerprint(mol, vecwords, _index.header.words, _index.header.fptype); 
	unsigned int words = _index.header.words;
	unsigned int* nextp = &_index.fptdata[0];
	unsigned int* ppat0 = &vecwords[0];
	register unsigned int* p;
	register unsigned int* ppat;
	register unsigned int i;
	register unsigned int a;
	for(i=0;i<dataSize; i++) //speed critical section
	{
		p=nextp;
		nextp += words;
		ppat=ppat0;
		a=0;
		while(p<nextp)
		{
			if(a=((*ppat) & (*p++)) ^ (*ppat++)) break;
		}
		if(!a)
		{
			candidates.push_back(i);
			if(candidates.size()>=MaxCandidates)
				break;
		}
	}

	if(i<_index.header.nEntries) //premature end to search
		cerr << "Stopped looking after " << i << " molecules." << endl;

	vector<unsigned int>::iterator itr;
	for(itr=candidates.begin();itr!=candidates.end();itr++)
	{
		SeekPositions.push_back(_index.seekdata[*itr]);
	}
	return true;
}

/////////////////////////////////////////////////////////
string FastSearch::ReadIndex(istream* pIndexstream)
{
	//Reads fs index from istream into member variables
	// but first checks whether it is already loaded
	FptIndexHeader headercopy = _index.header; 
	pIndexstream->read((char*)&(_index.header), sizeof(FptIndexHeader));
	if(memcmp(&headercopy,&(_index.header),sizeof(FptIndexHeader)))
	{
		pIndexstream->seekg(_index.header.headerlength);//allows header length to be changed

		unsigned int nwords = _index.header.nEntries * _index.header.words;
		_index.fptdata.resize(nwords);
		_index.seekdata.resize(_index.header.nEntries);

		pIndexstream->read((char*)&(_index.fptdata[0]), sizeof(unsigned int) * nwords);
		pIndexstream->read((char*)&(_index.seekdata[0]), sizeof(unsigned int) * _index.header.nEntries);
		if(pIndexstream->fail())
			*(_index.header.datafilename) = '\0';
	}
	return _index.header.datafilename;
}

//*******************************************************
//////////////////////////////////////////////////////////////
FastSearchIndexer::FastSearchIndexer(string datafilename, ostream* os, int FptType,
																		 int FptWords, unsigned int NfptEstimate)
{
	///Starts indexing process
	_indexstream = os;
	_fptype = FptType;

	_pindex= new FptIndex;
	_pindex->fptdata.reserve(NfptEstimate*FptWords);
	_pindex->seekdata.reserve(NfptEstimate);
	_pindex->header.headerlength = sizeof(FptIndexHeader);
	_pindex->header.words = FptWords;
	_pindex->header.fptype = FptType;
	strncpy(_pindex->header.datafilename, datafilename.c_str(), 255);
}

/////////////////////////////////////////////////////////////
FastSearchIndexer::~FastSearchIndexer()
{
	///Saves index file
	_pindex->header.nEntries = _pindex->seekdata.size();
	_indexstream->write((const char*)&_pindex->header, sizeof(FptIndexHeader));
	_indexstream->write((const char*)&_pindex->fptdata[0], _pindex->fptdata.size()*sizeof(unsigned int));
	_indexstream->write((const char*)&_pindex->seekdata[0], _pindex->seekdata.size()*sizeof(unsigned int));
	if(!_indexstream)
		cerr << "Difficulty writing index" << endl;
	delete _pindex;
}

///////////////////////////////////////////////////////////////
bool FastSearchIndexer::Add(OBMol& mol, streampos seekpos)
{
	///Adds a fingerprint if valid 
	
	vector<unsigned int> vecwords;
	if(!GetFingerprint(mol, vecwords, _pindex->header.words, _fptype) || vecwords.empty())
	{
		cerr << "Failed fingerprint for " << mol.GetTitle() << endl;
		return false;
	}
	for(int i=0;i<_pindex->header.words;i++)
	{
		_pindex->fptdata.push_back(vecwords[i]);
	}
	_pindex->seekdata.push_back(seekpos);
	return true;
}

//Fast sub-structure searching using fingerprints
/*
A description of fingerprints is at 
http://www.daylight.com/dayhtml/doc/theory/theory.finger.html)
http://www.mesaac.com/Fingerprint.htm
The data source is a text file with many molecules (several million?) in any format.
A requirement might be to find those molecules which has a substructure
specified by a SMARTS string. By conventional methods this is prohibitively slow.

To do it quickly, an index of the molecules is prepared containing a
their fingerprints - a condensed representation of their structure. 
This process may take a long time but only has to be done once. The index
can be searched very quickly in a process which eliminates almost all the
molecules, leaving only a few to be examined more carefully.

The fingerprints are hash values 32 bits or longer
with each bit representing one or more substructures. If a substructure
is present in the molecule, then a particular bit is 1. It may also be 1 if
some other substructure is present. The fast search prepares a fingerprint
from the requested SMARTS pattern. It finds those	molecules with which have
fingerprints which have 1 in all the same positions as the pattern fingerprint.
These candidate matches will include some false positives because of the use of
each bit for several substructures. But with proper design, the number of
molecules in this list will be a very small fraction of the total. It is
feasible to make a definitive substructure search by conventional means in
this reduced list even if it is slow. It is important, however, that the
checking of the fingerprints against the pattern fingerprint is fast.

The index has two table:
- an array of fingerprints of the molecules
- an array of the seek positions in the datasource file of all the molecules

To search in  fastsearch index file:
 
	ifs.open(indexname,ios_base::binary);
	if(!ifs) return false;

	FastSearch fs;
	vector<unsigned int> SeekPositions;
	string datafilename = fs.ReadIndex(&ifs);
	if(datafilename.empty() || !fs.Find(patternMol, SeekPositions, MaxCandidates))
		return false;
	
	ifstream datastream(datafilename);
	if(!datastream) return false;

	for(itr=SeekPositions.begin();itr!=SeekPositions.end();itr++)
	{
		datastream.seekg(*itr);
		... read the candidate molecule
	}
*/

}//Openbabel



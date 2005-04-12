#include <sstream>
#include <iostream>
#include <fstream>
#include "mol.h"
#include "obconversion.h"
#include "fingerprint.h"
#include "fastsearch.h"

using namespace std;
namespace OpenBabel {

class FastSearchFormat : public OBFormat
{
public:
	//Register this format type ID
	FastSearchFormat() : fsi(NULL)
	{
		OBConversion::RegisterFormat("fs",this);
	}
	
	virtual const char* Description() //required
	{ return
"FastSearching\n \
Uses molecular fingerprints in an index file.\n \
Make an index (slow) by using the fs format for writing:\n \
  obabel datafile.xxx index.fs   or\n \
  obabel datafile.xxx index.??? -ofs\n \
Search an index file by using the fs format for reading:\n \
  obabel index.fs -sSMILES outfile.yyy   or\n \
  obabel datafile.xxx -ifs -sSMILES outfile.yyy\n \
The structure spec can be a molecule from a file: -Spatternfile.zzz\n \
Options (when making index) e.g. -xN25000w2 \n \
-f#  finger print type, default <2>\n \
-w# number of 32bit words in fingerprint default<4>\n \
-N# approx number of molecules to be indexed\n \
";
	};

	virtual unsigned int Flags(){return READONEONLY | WRITEBINARY;};

private:
	FastSearchIndexer* fsi;

public:
	virtual bool ReadChemObject(OBConversion* pConv)
	{
		//Searches index file for structural matches
		//This function is called only once per search

		//Convert the SMARTS string to an OBMol
		OBMol patternMol;
		const char* p = pConv->IsOption('s',true);
		string txt;
		if(p) 
		{
			txt=p;
			txt = txt.substr(0,txt.find('"')); //now has SMARTS string
		}
		if(!txt.empty() && p!=NULL)
		{
			stringstream smarts(txt, stringstream::in);		
			OBConversion Convsm(&smarts);
			if(!Convsm.SetInFormat("smi")) return false;
			Convsm.Read(&patternMol);
		}

		// or Make OBMol from file in -S option
		p = pConv->IsOption('S',true);
		if(p && patternMol.Empty())
		{
			txt=p;
			txt = txt.substr(0,txt.find('"')); //now has filename
			int pos = txt.find_last_of('.');
			if(pos==string::npos)
			{
				cerr << "Filename of pattern molecule in -S option must have an extension" << endl;
				return false;
			}
			ifstream patternstream(txt.c_str());
			if(!patternstream)
			{
				cerr << "Cannot open " << txt << endl;
				return false;
			}
			stringstream smiles(stringstream::out);		

			OBConversion PatternConv(&patternstream,&smiles);
			PatternConv.SetOneObjectOnly();
			if(PatternConv.SetInAndOutFormats(txt.substr(pos+1).c_str(),"smi"))
				PatternConv.Read(&patternMol);

			//Convert to SMILES and generate a -s option for use in the final filtering
			PatternConv.Write(&patternMol);
			string genopts  = pConv->GetGeneralOptions();
			//remove name to leave smiles string
			string smilesstr(smiles.str());
			pos = smilesstr.find(' ');
			if(pos!=string::npos)
				smilesstr = smilesstr.substr(0,pos);
			genopts += "s\"" + smilesstr + "\"";
			pConv->SetGeneralOptions(genopts.c_str());
		}

		if(patternMol.Empty())
		{
			cerr << "Cannot derive a molecule from the -s or -S options" << endl;
			return false;
		}

		//Derive index name
		string indexname = pConv->GetInFilename();
		int pos=indexname.find_last_of('.');
		if(pos!=string::npos)
		{
			indexname.erase(pos);
			indexname += ".fs";
		}

		//Have to open input stream again because needs to be in binary mode
		ifstream ifs(indexname.c_str(),ios_base::binary);
		if(!ifs)
		{
			cerr << "Couldn't open " << indexname << endl;
			return false;
		}

		//TODO Better way of inputting this
		int MaxCandidates = 4000;

		FastSearch fs;
		vector<unsigned int> SeekPositions;
		string datafilename = fs.Find(patternMol, &ifs, SeekPositions, MaxCandidates);

		cerr << SeekPositions.size() << " candidates " << endl;
		if(datafilename.empty())
		{
			cerr << "Difficulty reading from index " << indexname << endl;
			return false;
		}
	
		if(SeekPositions.size()!=0)
		{
			ifstream datastream;
			if(pConv->GetInFilename()==indexname)
			{
				//need to open the datafile and put it in pConv
				//datafile name derived from index file probably won't have a file path
				//but indexname may. Derive a full datafile name
				string path;
				int pos = indexname.find_last_of("/\\");
				if(pos==string::npos)
					path = datafilename;
				else
					path = indexname.substr(0,pos+1) + datafilename;
				
				datastream.open(path.c_str());
				if(!datastream)
				{
					cerr << "Difficulty opening " << path << endl;
					return false;
				}
				pConv->SetInStream(&datastream);
			}
			
			//Input format is currently fs; set it appropriately
			if(!pConv->SetInAndOutFormats(pConv->FormatFromExt(datafilename.c_str()),pConv->GetOutFormat()))
					return false;
			
			//Output the candidate molecules, filtering through s filter
			vector<unsigned int>::iterator itr;
			for(itr=SeekPositions.begin();itr!=SeekPositions.end();itr++)
			{
				datastream.seekg(*itr);
//				datastream.seekg(*itr - datastream.tellg(), ios_base::cur); //Avoid retrieving start
//				string ln;
//				getline(datastream,ln);
//				cerr << ln << endl;
				pConv->SetOneObjectOnly();
				pConv->Convert();
			}
		}

		return false;	//To finish	
	};

/////////////////////////////////////////////////////
	virtual bool WriteChemObject(OBConversion* pConv)
	{
		//Prepares an index file
		if(pConv->GetOutputIndex()==0)
			cerr << "This will prepare an index of " << pConv->GetInFilename()
			<<  " and may take some time..." << flush;
		
		OBStopwatch sw;
		sw.Start();
		
		ostream* pOs = pConv->GetOutStream();
		bool NewOstreamUsed=false;
		if(fsi==NULL)
		{
			//First pass sets up FastSearchIndexer object
			if(pOs==&cout)
			{
				//No index filename specified
				//Derive index name from datfile name
				string indexname=pConv->GetInFilename();
				int pos=indexname.find_last_of('.');
				if(pos!=string::npos)
					indexname.erase(pos);
				indexname += ".fs";

				pOs = new ofstream(indexname.c_str(),ofstream::binary);
				if(!pOs->good())
				{
					cerr << "Cannot open " << indexname << endl;
					return false;
				}
				NewOstreamUsed=true;
			}

			const char* p = pConv->IsOption('N');
			unsigned int nmols = p ? atoi(p) : 200; //estimate of number of molecules
			p = pConv->IsOption('w');
			unsigned int nwords = p ? atoi(p) : 4; //number 32bit words in each fingerprint
			if(nmols==0 || nwords==0)
			{
				cerr << "Bad -x option" << endl;
				return false;
			}
			int fptype=0; //fingerprint type
			p=pConv->IsOption('f');
			if(p)
				fptype=atoi(p);

			//Prepare name without path
			string datafilename = pConv->GetInFilename();
			if(datafilename.empty())
			{
				cerr << "No datafile! " << endl;
				return false;
			}
			int pos = datafilename.find_last_of("/\\");
			if(pos!=string::npos)
				datafilename=datafilename.substr(pos+1);
			fsi = new FastSearchIndexer(datafilename, pOs, fptype, nwords, nmols);
		}

		//All passes provide a molecule for indexing
		OBBase* pOb = pConv->GetChemObject();
		OBMol* pmol = dynamic_cast<OBMol*> (pOb);
		if(pmol)
		{
			streampos seekpos = pConv->GetInPos();
			fsi->Add(*pmol, seekpos );
		}

		if(pmol==NULL || pConv->IsLast())
		{
			//Last pass 
			if(NewOstreamUsed)
				delete pOs;

			delete fsi; //saves index file
			//return to starting conditions
			fsi=NULL;

			cerr << "\n It took " << sw.Elapsed() << endl;
		}
		delete pOb;
		return true;
	};
};

//Make an instance of the format class
FastSearchFormat theFastSearchFormat;

//***************************************************************
//***************************************************************
FastSearch::~FastSearch(){}

//////////////////////////////////////////////////
string FastSearch::Find(OBMol& mol, istream* pIndexstream, vector<unsigned int>& SeekPositions,
												int MaxCandidates)
{
	//Finds molecules in datafilename (which must previously have been indexed)
	//that match a substructure in mol. datafile name can also be the name of the index file.
	//The positions of the candidate matching molecules in the original datafile are returned.
	//When called, datafilename can be either the real data file or the index.
	//On return it will have been altered to the real datafile name
	
	vector<unsigned int>candidates; //indices of matches from fingerprint screen
	candidates.reserve(MaxCandidates);
	FptIndex index; //Large size is in vectors whose content is probably not on the stack
	if(!ReadIndex(pIndexstream, &index))
		return NULL;

	unsigned int dataSize = index.header.nEntries;
	vector<unsigned int> vecwords;
	GetFingerprint(mol, vecwords, index.header.words, index.header.fptype); 
	unsigned int words = index.header.words;
	unsigned int* nextp = &index.fptdata[0];
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

	if(i<index.header.nEntries) //premature end to search
		cerr << "Stopped looking after " << i << " molecules." << endl;

	vector<unsigned int>::iterator itr;
	for(itr=candidates.begin();itr!=candidates.end();itr++)
	{
		SeekPositions.push_back(index.seekdata[*itr]);
	}
	string df(index.header.datafilename);
	return df;
}

/////////////////////////////////////////////////////////
bool FastSearch::ReadIndex(istream* pIndexstream,FptIndex* pindex)
{
	//Reads fs index from istream into memory
	pIndexstream->read((char*)&pindex->header, sizeof(FptIndexHeader));
	pIndexstream->seekg(pindex->header.headerlength);//allows header length to be changed

	unsigned int nwords = pindex->header.nEntries * pindex->header.words;
	pindex->fptdata.resize(nwords);
	pindex->seekdata.resize(pindex->header.nEntries);

	pIndexstream->read((char*)&pindex->fptdata[0], sizeof(unsigned int)*nwords);
	pIndexstream->read((char*)&pindex->seekdata[0], sizeof(unsigned int)*pindex->header.nEntries);
	if(pIndexstream->fail())
		return false;
	return true;
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
	int ndata = _pindex->fptdata.size();
	int nspos = _pindex->seekdata.size();
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
*/

/* Chosen way
Keyboard interface
	To search for SMARTS matches	
obabel databasefile.fs outfile.yyy -s"smarts"
or
obabel databasefile.xxx outfile.yyy -s"smarts" -ifs
or
obabel databasefile.xxx outfile.yyy -sPatternfile.zzz -ifs

	To make index
obabel databasefile.xxx -ofs   (makes databasefile.fs)
or
obabel databasefile.xxx namedindex.fs 

fsformat
does searching in ReadChemObject() and indexing in WriteChemObject()

*/


}//Openbabel



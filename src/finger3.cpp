/**********************************************************************
finger3.cpp: Fingerprint based on list of SMARTS patterns
Copyright (C) 2005 Chris Morley
 
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
#include "mol.h"
#include <fstream>
#include <map>
#include <string>
#include "fingerprint.h"

using namespace std;
namespace OpenBabel
{
/// \brief Fingerprint based on list of SMARTS patterns ID="FP3"
class PatternFP  : public OBFingerprint
{
public:
	PatternFP(const char* ID, const char* filename=NULL, 
			bool IsDefault=false) : OBFingerprint(ID, IsDefault)
	{
		if(filename==NULL)
			_patternsfile="patterns.txt";
		else
			_patternsfile = filename;
	};
	
	virtual string Description()
	{	
		string desc("SMARTS patterns specified in the file ");
		return (desc + _patternsfile);
	};

	//Each bits represents a single substructure; no need for confirmation when substructure searching
	virtual unsigned int Flags() { return FPT_UNIQUEBITS;}; 

	bool GetFingerprint(OBBase* pOb, vector<unsigned int>&fp, int nbits) 
	{
		OBMol* pmol = dynamic_cast<OBMol*>(pOb);
		if(!pmol)
			return false;
		
		//Read patterns file if it has not been done already
		static vector<string> smartsStrings;
		if(smartsStrings.empty())
		{	
			char* datadir = getenv("BABEL_DATADIR");
			if(!datadir)
				datadir = BABEL_DATADIR;
			if(datadir)
			{
				_patternsfile = "/" + _patternsfile;
				_patternsfile = datadir  + _patternsfile;
			}

			ifstream ifpatterns(_patternsfile.c_str());
			if(!ifpatterns)
			{
				cerr << "Cannot open " << _patternsfile << endl;
				return false;
			}
			string smarts, commentline;
			if(!getline(ifpatterns,commentline)) return false;
			
			while(ifpatterns.good())
			{
				if(getline(ifpatterns,smarts))
				{
					int pos = smarts.find(' ');
					if(pos!=string::npos)
						smarts = smarts.substr(0,pos);
					smartsStrings.push_back(smarts);
				}
			}
		}

		//Make fp size the smallest power of two to contain the patterns
		int n=bitsperint;
		while(n<smartsStrings.size())n*=2;
		fp.resize(n/bitsperint);

		for(n=0;n<smartsStrings.size();++n)
		{
			OBSmartsPattern sp;
			sp.Init(smartsStrings[n]);
			if(sp.Match(*pmol))
				SetBit(fp, n);
		}

		if(nbits)
			Fold(fp, nbits);
		return true;
	};
protected:
	string _patternsfile;
};

//***********************************************
//Make a global instance
PatternFP thePatternFP("FP3");
//***********************************************

/*! \class PatternFP
A bit is set when there is a match to one of a list
of SMARTS patterns in the file (default is patterns.txt).
Looks for this file first in the folder in the environment variable
BABEL_DATADIR, then in the folder specified by the macro BABEL_DATADIR
(probably set in babelconfig.h), and then in the current folder. 

The first line of this file is a comment.
On each subsequent line there is a SMARTS string and anything 
after a space is ignored.

Additional fingerprint types using patterns in different files
can be made by just declaring separate instances like:
\code
PatternFP myPatternFP("myFP", "myPatternFile.txt");
\endcode
*/

}//namespace

//! \file finger3.cpp
//! \brief fingerprint3 based on list of SMARTS patterns

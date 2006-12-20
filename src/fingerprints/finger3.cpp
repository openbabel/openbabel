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

#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>

#include <fstream>
#include <map>
#include <string>

#include <openbabel/fingerprint.h>

using namespace std;
namespace OpenBabel
{
/// \brief Fingerprint based on list of SMARTS patterns ID="FP3"
class PatternFP  : public OBFingerprint
{
private:
	vector<string> smartsStrings;
protected:
	string _patternsfile;

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
		if(smartsStrings.empty())
			ReadPatternFile(_patternsfile, smartsStrings);

		//Make fp size the smallest power of two to contain the patterns
		unsigned int n=Getbitsperint();
		while(n<smartsStrings.size())n*=2;
		fp.resize(n/Getbitsperint());

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

	bool ReadPatternFile(const string& filename, vector<string>& lines)
	{	
		//Reads two types of file: SMARTS + comments and vice versa
		//depending on whether the first line is #Comments after SMARTS
		//Output strings in vector are SMARTS + comments
		string file = filename;
		ifstream ifs1, ifs2, *ifsP;
		char* datadir = getenv("BABEL_DATADIR");
		if(!datadir)
			datadir = BABEL_DATADIR;
		if(datadir)
		{
		  file = datadir;
		  file += FILE_SEP_CHAR;
		  file += filename;
		}

		ifs1.open(file.c_str());
		ifsP = &ifs1;
		if(!(*ifsP))
		{
		  file = datadir;
		  file += FILE_SEP_CHAR;
		  file += BABEL_VERSION;
		  file += FILE_SEP_CHAR + filename;
		  ifs2.open(file.c_str());
		  ifsP = &ifs2;
		}

		if(!(*ifsP))
		  {
#ifdef HAVE_SSTREAM 
		    stringstream errorMsg; 
#else 
		    strstream errorMsg; 
#endif 
		    errorMsg << "Cannot open " << filename << endl; 
		    obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError); 
		    return false;
		  }
		string smarts, formatline;

		if(!getline(*ifsP, formatline)) return false;
		if(formatline=="#Comments after SMARTS")
		{
			while(ifsP->good())
			{
				if( getline(*ifsP,smarts) 
						&& smarts.size() > 0 
						&& smarts[0] != '#')
					lines.push_back(smarts); //leave the comments in
			}
		}
		else
		{
			// Christian Laggner's format: SMARTS at end of line
			while(ifsP->good())
			{
				if( getline(*ifsP,smarts) && smarts[0]!='#')
				{
          string::size_type pos = smarts.find(':');
					if(pos!=string::npos)
					{
						pos = smarts.find_first_not_of(" \t", pos+1);
						if(pos!=string::npos)
							lines.push_back(smarts.substr(pos) + ' ' + smarts.substr(0,pos));
					}
				}
			}
		}
		return true;
	}
};

//***********************************************
//Make a global instance
PatternFP thePatternFP("FP3");

PatternFP FP4PatternFP("FP4", "SMARTS_InteLigand.txt");
//***********************************************

/*! \class PatternFP
A bit is set when there is a match to one of a list
of SMARTS patterns in the file (default is patterns.txt).
Looks for this file first in the folder in the environment variable
BABEL_DATADIR, then in the folder specified by the macro BABEL_DATADIR
(probably set during compilation in babelconfig.h), and then in the current folder.

On each line there is a SMARTS string and anything 
after a space is ignored. Lines starting with # are ignored.

Additional fingerprint types using patterns in different files
can be made by just declaring separate instances like:
\code
PatternFP myPatternFP("myFP", "myPatternFile.txt");
\endcode
*/

}//namespace

//! \file finger3.cpp
//! \brief fingerprint3 based on list of SMARTS patterns

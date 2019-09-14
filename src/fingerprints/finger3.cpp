/**********************************************************************
finger3.cpp: Fingerprints based on list of SMARTS patterns
Copyright (C) 2005 Chris Morley

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
#include <openbabel/mol.h>
#include <openbabel/parsmart.h>
#include <openbabel/oberror.h>
#include <sstream>
#include <fstream>
#include <map>
#include <string>

#include <openbabel/fingerprint.h>

using namespace std;
namespace OpenBabel
{
/// \brief Fingerprint based on list of SMARTS patterns
class PatternFP  : public OBFingerprint
{
private:
  struct pattern
  {
    string smartsstring;
    OBSmartsPattern obsmarts;
    string description;
    int numbits;
    int numoccurrences;
    int bitindex;
  };
  vector<pattern> _pats;
  int _bitcount;
  string _version;

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
  }

/////////////////////////////////////////////////////////////////////////////
  virtual const char* Description()
  {
    static string desc;
    //Read patterns file if it has not been done already,
    //because we need _bitcount and _version updated

    // _bitcount and _version are available only after the datafile has been parsed.
    // This is a burden on normal operation (Description() gets called on startup from OBDefine),
    // so the secondline is present only after the fingerprint has been used.
    // the
    string secondline;
    if(!_pats.empty())
      secondline = "\n" + toString(_bitcount) + " bits. Datafile version = " +  _version;
    desc = "SMARTS patterns specified in the file " + _patternsfile
      + secondline
      + "\nPatternFP is definable";
    return (desc.c_str());
  }

//////////////////////////////////////////////////////////////////////////////
  //Each bit represents a single substructure
  virtual unsigned int Flags() { return FPT_UNIQUEBITS;};

///////////////////////////////////////////////////////////////////////////////
  virtual PatternFP* MakeInstance(const std::vector<std::string>& textlines)
  {
    return new PatternFP(textlines[1].c_str(),textlines[2].c_str());
  }

////////////////////////////////////////////////////////////////////////////////
  virtual bool GetFingerprint(OBBase* pOb, vector<unsigned int>&fp, int foldbits)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(!pmol)
      return false;

    //This fingerprint is constructed from a molecule with no explicit hydrogens.
    pmol->DeleteHydrogens();

    unsigned int n;
    //Read patterns file if it has not been done already
    if(_pats.empty())
      ReadPatternFile(_version);

    //Make fp size the smallest power of two to contain the patterns
    n=Getbitsperint();
    while(n < _bitcount)
      n*=2;
    fp.resize(n/Getbitsperint());

    n=0; //bit position
    vector<pattern>::iterator ppat;
    for(ppat=_pats.begin();ppat!=_pats.end();++ppat)
    {
      if(ppat->numbits //ignore pattern if numbits==0
        && ppat->obsmarts.Match(*pmol, ppat->numoccurrences==0))//do single match if all that's needed
      {
        /* Set bits in the fingerprint depending on the number of matches in the molecule
           and the parameters, numbits and numoccurrences, in the pattern.
           The pattern will set or clear numbits bits in the fingerprint.
           They will be in numoccurrences+1 groups, each containing an approximately
           equal number of bits.
           The first group of bits will be set if numMatches > numoccurences;
           The second group will be set if numMatches > numoccurrences - 1;
           and so on.
           So with a pattern with numbits = 4 and numoccurences = 2,
           the groups would be 1, 1, and 2 bits.
           A molecule with
              1 match to the pattern would give 0011
              2 matches to the pattern would give 0111
              3 or more matches to the pattern would give 1111
        */
        int numMatches = ppat->obsmarts.GetUMapList().size();
        int num =  ppat->numbits, div = ppat->numoccurrences+1, ngrp;

        int i = n;
        while(num)
        {
          ngrp = (num -1)/div-- +1; //rounds up
          num -= ngrp;
          while(ngrp--)
            if (numMatches > div) {
              SetBit(fp,i);
            }
          i++;
        }
      }
      n += ppat->numbits;
    }

    if(foldbits)
      Fold(fp, foldbits);
    return true;
  }

  /////////////////////////////////////////////////////////////////////
  bool ReadPatternFile(string& ver)
  {
    //Reads three types of file. See below
    ifstream ifs;
	  stringstream errorMsg;

    if (OpenDatafile(ifs, _patternsfile).length() == 0)
    {
      errorMsg << "Cannot open " << _patternsfile << endl;
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
      return false;
    }

    string line;
    if(!getline(ifs, line)) //first line
      return false;
    bool smartsfirst = (Trim(line)=="#Comments after SMARTS");

    _bitcount=0;
    bool indata=false;
    do
    {
      if(Trim(line).size()>0 && line[0]!='#')
      {
        pattern p;
        p.numbits=1; p.numoccurrences=0; //default values
        p.bitindex = _bitcount;
        istringstream ss(line);
        indata = true;
        if(smartsfirst)
        {
          if(isdigit(line[0]))
          {
            if(!ParseRDKitFormat(ss, p))
              continue;
          }
          else
            //Original format, which looks like:
            //  SMARTS description
            ss >> p.smartsstring >> p.description;
        }
        else
        {
          // Christian Laggner's format:
          //  description: SMARTS [occurrences [numbits]]
          getline(ss, p.description, ':');
          ss >> p.smartsstring;
          ss >> p.numoccurrences >> p.numbits;
        }

        if(!p.obsmarts.Init(p.smartsstring))
        {
          obErrorLog.ThrowError(__FUNCTION__,
            "Faulty SMARTS: " + p.description + ' ' + p.smartsstring, obError);
          continue;
        }
        _pats.push_back(p);
        _bitcount += p.numbits;
      }
      else if(!indata)
      {
        //Find version number
        string::size_type pos = line.find("Version");
        if(pos!=string::npos)
          pos+=8;
        else if(line.find("Extracted from RDKit")!=string::npos)
        {
          pos=20;
          while((pos=line.find('r',pos))!=string::npos)
            if(isdigit(line[++pos]))
              break;
        }
        if(pos!=string::npos)
        {
          ver=line.substr(pos) + ' ';//space fixes bug in while() when number at end of line
          pos=1;
          while(isdigit(ver[++pos]));
          ver.erase(pos);
        }
      }
    }while(getline(ifs,line));

    if (ifs)
      ifs.close();
    return true;
  }

///////////////////////////////////////////////////////////////////////////////
  virtual string DescribeBits(const vector<unsigned int> fp, bool bSet=true)
  {
    //checkmol-type output with tab separated functional group names
    stringstream ss;
    vector<pattern>::iterator ppat;
    for(ppat=_pats.begin();ppat!=_pats.end();++ppat)
    {
      int n = ppat->bitindex;
      int num =  ppat->numbits, div = ppat->numoccurrences+1, ngrp;
      while(num) //for each group of bits
      {
        ngrp = (num + div -1)/div--; //rounds up
        num -= ngrp;
        if(GetBit(fp, n) == bSet)
        {
          ss << ppat->description;
          if(div>0)
            ss << '*' << div+1;
          ss << '\t' ;
          break; //ignore the bits signifying a smaller number of occurrences
        }
        n += ngrp;
      }
    }
    ss << endl;
    return ss.str();
  }

///////////////////////////////////////////////////////////////////////////////////
  bool ParseRDKitFormat(istringstream& ss, pattern& p)
  {
    //rdkit format, e.g.
    //  14:('[S,s]-[S,s]',0), # S-S
    const int dum = 20; //an arbitrary number in case delimiters in ignore statements not found
    string number, comment;
    getline(ss, number, ':');
    ss.ignore(dum, '\'');
    getline(ss, p.smartsstring, '\'');
    if(p.smartsstring[0]=='?') //ignore patterns with SMARTS '?'
      p.smartsstring="[999]";//this seems to match nothing;  was return false;
    ss.ignore(dum,',');
    ss >> p.numoccurrences;
    ss.ignore(dum,'#');
    getline(ss, comment);

    //description is number + edited commment
    Trim(comment);
    string::size_type pos;
    pos = comment.find("FIX");
    if(pos==string::npos)
      pos = comment.find("*NOTE*");
    if(pos!=string::npos)
      comment.erase(pos);
    p.description = number + ": " + comment;
    return true;
  }


}; //class PatternFP

//***********************************************
//Make a global instance
PatternFP FP3PatternFP("FP3");
PatternFP FP4PatternFP("FP4", "SMARTS_InteLigand.txt");
//***********************************************

/*! \class PatternFP
A bit is set when there is a match to one of a list
of SMARTS patterns in the datafile, which is specified in the constructor.
If no filename is given, the default filename is patterns.txt.
Fingerprints can be made by declaring a global variable, as in:

PatternFP FP4PatternFP("FP4", "SMARTS_InteLigand.txt");

Alternatively, an entry in plugindefines.txt like:

PatternFP
MACCS          #ID of this fingerprint type
MACCS.txt      #File containing the SMARTS patterns

defines a fingerprint without the need to recompile.

Three file formats are supported:
 - the preferred format (e.g. SMARTS_InteLigand.txt in FP4)
 - the original format (patterns.txt has an incomplete set of SMARTS patterns)
 - a format made by extracting from an RDKit file (MACCS.txt)
The last two require the first line to be:
#Comments after SMARTS

Lines starting with # are ignored.
For the preferred format each line is of the form:
description: SMARTS [occurrences [numbits]]
A bit is set in the fingerprint for ach SMARTS pattern matched.
The optional integer parameters refine this behaviour; the most obvious uses are:
 - if <occurrences> is present and greater than its default value of 0, the bit
   is set only if the number of matches to the pattern is greater than <occurences>.
 - if <occurences> is 0 and <numbits> is greater than its default value of 1, then
   the fingerprint has <numbits> bits set if there is a match. This gives greater weight
   to the pattern for use in similarity measures like Tanimoto.
 - if the parameters are n-1 and n and the number of matches is n,
   a bit is set for each of the conditions n>=m, n>=m-1, ... , n>=1
   This can be used to distinguish structures with many similar atoms like n-alkanes.
The use of other values for the parameters, which can be any positive integer, can give
other analogous behaviours. If numbits is 0 the pattern is ignored.
*/

}//namespace

//! \file finger3.cpp
//! \brief fingerprints based on list of SMARTS patterns

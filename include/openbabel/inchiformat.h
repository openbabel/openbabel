/**********************************************************************
Copyright (C) 2005,2006,2007 Chris Morley

Based on the IUPAC InChI reference software, which is distributed
under the GNU LGPL:
Copyright (C) 2005 The International Union of Pure and Applied Chemistry
IUPAC International Chemical Identifier (InChI) (contact:secretariat@iupac.org)

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
#include <openbabel/obconversion.h>
#include <openbabel/obmolecformat.h>

#include "inchi_api.h"
#ifdef HAVE_SSTREAM
#include <sstream>
#else
#include <strstream>
#endif
#include <set>
#include <vector>
#include "openbabel/chiral.h"

using namespace std;
namespace OpenBabel
{
extern string GetInChI(istream& is);

class InChIFormat : public OBMoleculeFormat
{
public:
  InChIFormat()
  {
    OBConversion::RegisterFormat("inchi",this);
    OBConversion::RegisterOptionParam("n", this, 0, OBConversion::INOPTIONS);
    OBConversion::RegisterOptionParam("t", this);
    OBConversion::RegisterOptionParam("l", this);
    OBConversion::RegisterOptionParam("X", this, 1, OBConversion::OUTOPTIONS);
    OBConversion::RegisterOptionParam("K", this, 0, OBConversion::OUTOPTIONS);
    OBConversion::RegisterOptionParam("F", this, 0, OBConversion::OUTOPTIONS);
    OBConversion::RegisterOptionParam("X", this, 1, OBConversion::INOPTIONS);
    OBConversion::RegisterOptionParam("T", this, 1, OBConversion::OUTOPTIONS);
  }

  virtual const char* Description()
  {
    return 
    "InChI format\n"
    "IUPAC/NIST molecular identifier\n"
    "Write options, e.g. -xat\n"
    //" n do not use 'recommended' InChI options\n"
    " X <Option string> List of additional InChI options\n"
    //" F include fixed hydrogen layer\n"
    //" M include bonds to metal\n"
    " t add molecule name\n"
    " a output auxilliary information\n"
    " K output InChIKey\n"
    " w don't warn on undef stereo or charge rearrangement\n"
    " l display InChI log\n"
    " u output only unique molecules\n"
    " U output only unique molecules and sort them\n"
    " e compare first molecule to others\n"
    " T <param> truncate InChI, /nostereo etc.\n\n"

    "Input options, e.g. -an\n"
    " X <Option string> List of InChI options\n"
    " n molecule name follows InChI on same line\n"
    " a add InChI string to molecule name\n\n"
    "Currently the output is standard InChI only."
    "InChI options may be reintroduced later."
    " The InChI options should be space delimited in a single quoted string.\n"
    " See InChI documentation for possible options.\n\n"

    " Truncation parameters used with -xT\n"
    "/formula  formula only\n"
    "/connect  formula and connectivity only\n"
    "/nostereo ignore E/Z and sp3 stereochemistry\n"
    "/sp3      ignore sp3 stereochemistry\n"
    "/noEZ     ignore E/Z steroeochemistry\n"
    "/nochg    ignore charge and protonation\n"
    "/noiso    ignore isotopes\n\n"
;
  };

  virtual const char* SpecificationURL()
  { return "http://www.iupac.org/inchi/";};

  virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
  virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);
  virtual int  SkipObjects(int n, OBConversion* pConv);

  static char   CompareInchi(const string& Inchi1, const string& Inchi2);
  static string InChIErrorMessage(const char ch);

  ///Removes layers or truncates InChi, according to\param spec
  ///which can contain any number of:/formula /connect /nostereo /nosp3 /noEZ /nochg /noiso
  static bool EditInchi(std::string& inchi, std::string& spec);

  ///Compare std::strings with embedded numbers so that 
  // "a6b" (or "a06b") is less than "a15b"
  // and "CH4" is less than "C2H6"
  // and "CH4" is less than "ClH" (hydrogen chloride)
  static struct InchiLess
    : public binary_function<const string&, const string&, bool>
  {
    bool operator()(const string& s1, const string& s2) const
    {
      //stop at the first space or the end of the strings
      string::const_iterator p1=s1.begin(), p2=s2.begin(),
        p1end=find(s1.begin(), s1.end(), ' '), p2end=find(s2.begin(), s2.end(), ' ');

      while( p1<p1end && p2<p2end)
      {
        int n1=-1,n2=-1;
        if(isdigit(*p1))
          {
            n1 = atoi(&*p1);
            //skip over number
            while(p1!=s1.end() && isdigit(*p1++)); --p1;
          }
        if(isdigit(*p2))
          {
            n2 = atoi(&*p2);
            while(p2!=s2.end() && isdigit(*p2++)); --p2;
          }
        if(n1<0 && n2 < 0)
          {
            //neither numbers
            if(*p1 != *p2)
        return *p1 < *p2;
          }
        else if(n1>=0 && n2>0)
          {
            //both numbers
            if(n1!=n2)
        return n1 < n2;
          }
        else if(n1>0)
          return islower(*p2)!=0;
        else if(n2>0)
          return !islower(*p1);

        ++p1; ++p2; // iterate
      } // while loop
      return false; //identical
    }
  };

private:
  ///Erases the layer starting with \param str and, if \param all is true, all the subsequent ones 
  static void RemoveLayer (std::string& inchi, const std::string& str, bool all=false);

private:
  OBAtom* GetCommonAtom(OBBond* pb1, OBBond* pb2);
  char* GetInChIOptions(OBConversion* pConv, bool Reading);

  typedef	set<string, InchiLess> nSet;
  nSet allInchi;
  string firstInchi;
  string firstID;
};

//*****************************************************
class InChICompareFormat : public OBMoleculeFormat
{
public:
  InChICompareFormat()
  {
      OBConversion::RegisterFormat("k",this);
  }
  virtual const char* Description() //required
  {
    return 
      "Compare molecules using InChI\n"
      "The first molecule is compared with the rest\n"
      "e.g. babel first.smi second.mol third.cml -ok\n"
      "Same as  -oinchi -xet  and can take the same options as InChI format.\n";
  }
  virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);
  virtual unsigned int Flags() { return NOTREADABLE;};
};

}//namespace OpenBabel

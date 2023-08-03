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
#ifndef OB_INCHIFORMAT_H
#define OB_INCHIFORMAT_H
#include <openbabel/babelconfig.h>
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
#include <cstdlib>
#include <algorithm>

namespace OpenBabel
{
extern std::string GetInChI(std::istream& is);

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
    "IUPAC/NIST molecular identifier\n\n"

    "Write Options, e.g. -xa\n"
    "    Standard InChI is written unless certain InChI options are used\n \n"
    " K output InChIKey only\n"
    " t add molecule name after InChI\n"
    " w ignore less important warnings\n"
    "    These are:\n"
    "    \'Omitted undefined stereo\'\n"
    "    \'Charges were rearranged\'\n"
    "    \'Proton(s) added/removed\'\n"
    "    \'Metal was disconnected\'\n"
    " a output auxiliary information\n"
    " l display InChI log\n"
    " r recalculate InChI; normally an input InChI is reused\n"
    " s recalculate wedge and hash bonds(2D structures only)\n \n"
    "    **Uniqueness options** (see also ``--unique`` and ``--sort`` which are more versatile)\n"
    " u output only unique molecules\n"
    " U output only unique molecules and sort them\n"
    " e compare first molecule to others\n"
    "    This can also be done with :ref:`InChICompare format <Compare_molecules_using_InChI>`::\n \n"
    "      babel first.smi second.mol third.cml -ok\n \n"
    " T <param> truncate InChI according to various parameters\n"
    "    See below for possible truncation parameters.\n"
    
    " X <Option string> Additional InChI options\n"
    "    See InChI documentation.\n"
    "    These options should be space delimited in a single quoted string.\n \n"
    "    - Structure perception (compatible with stdInChI): ``NEWPSOFF``, ``DoNotAddH``, ``SNon``\n"
    "    - Stereo interpretation (produces non-standard InChI): ``SRel``, ``SRac``,\n"
    "      ``SUCF``, ``ChiralFlagON``, ``ChiralFlagOFF``\n"
    "    - InChI creation options (produces non-standard InChI): ``SUU``, ``SLUUD``,\n"
    "      ``FixedH``, ``RecMet``, ``KET``, ``15T``\n \n"
    "    The following options are for convenience, e.g. ``-xF``\n"
    "    but produce non-standard InChI.\n"
    " F include fixed hydrogen layer\n"
    " M include bonds to metal\n\n"

    "Read Options, e.g. -an\n"
    " X <Option string> List of InChI options\n"
    " n molecule name follows InChI on same line\n"
    " a add InChI string to molecule name\n\n"

    "Truncation parameters used with ``-xT``:\n\n"
    "/formula   formula only\n"
    "/connect   formula and connectivity only\n"
    "/nostereo  ignore E/Z and sp3 stereochemistry\n"
    "/nosp3       ignore sp3 stereochemistry\n"
    "/noEZ      ignore E/Z steroeochemistry\n"
    "/nochg     ignore charge and protonation\n"
    "/noiso     ignore isotopes\n\n"
    "Note that these can also be combined, e.g. ``/nochg/noiso``\n"
;
  };

  virtual const char* SpecificationURL()
  { return "http://www.iupac.org/inchi/";};

  virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
  virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);
  virtual int  SkipObjects(int n, OBConversion* pConv);

  static char CompareInchi(const std::string& Inchi1, const std::string& Inchi2);
  static std::string InChIErrorMessage(const char ch);

  ///Removes layers or truncates InChi, according to \param spec
  ///which can contain any number of:/formula /connect /nostereo /nosp3 /noEZ /nochg /noiso
  /// @param inchi The inchi string
  static bool EditInchi(std::string& inchi, std::string& spec);

  ///Compare std::strings with embedded numbers so that
  // "a6b" (or "a06b") is less than "a15b"
  // and "CH4" is less than "C2H6"
  // and "CH4" is less than "ClH" (hydrogen chloride)
  struct InchiLess
    : public std::binary_function<const std::string&, const std::string&, bool>
  {
    bool operator()(const std::string& s1, const std::string& s2) const
    {
      //stop at the first space or the end of the strings
      std::string::const_iterator p1=s1.begin(), p2=s2.begin(),
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
  void SaveInchi(OBMol* pmol, const std::string& s);

  typedef std::set<std::string, InchiLess> nSet;
  nSet allInchi;
  std::string firstInchi;
  std::string firstID;
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
      "A utility format that allows you to compare molecules using their InChIs\n"
      "The first molecule is compared with the rest, e.g.::\n\n"

      "  babel first.smi second.mol third.cml -ok\n\n"

      "This is the same as using ``-oinchi -xet`` and can take the same options as InChI format\n"
      "(see :ref:`InChI_format`).\n";
  }
  virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);
  virtual unsigned int Flags() { return NOTREADABLE;};
};

//*****************************************************
class InChIKeyFormat : public OBMoleculeFormat
{
public:
  InChIKeyFormat()
  {
      OBConversion::RegisterFormat("inchikey",this);
  }
  virtual const char* Description() //required
  {
    return
      "InChIKey\n"
      "A hashed representation of the InChI.\n\n"

      "The InChIKey is a fixed-length (27-character) condensed digital\n"
      "representation of an InChI, developed to make it easy to perform\n"
      "web searches for chemical structures.\n\n"

      "An InChIKey consists of 14 characters (derived from the connectivity\n"
      "layer in the InChI), a hyphen, 9 characters (derived from the\n"
      "remaining layers), a character indicating the InChI version, a hyphen\n"
      "and a final checksum character. Contrast the InChI and InChIKey of the\n"
      "molecule represented by the SMILES string `CC(=O)Cl`::\n\n"

      "  obabel -:CC(=O)Cl -oinchi\n"
      "  InChI=1S/C2H3ClO/c1-2(3)4/h1H3\n\n"

      "  obabel -:CC(=O)Cl -oinchikey\n"
      "  WETWJCDKMRHUPV-UHFFFAOYSA-N\n\n"
      
      "This is the same as using ``-oinchi -xK`` and can take the same options\n"
      "as the InChI format (see :ref:`InChI_format`)::\n\n"

      "  obabel -:CC(=O)Cl -oinchi -xK\n"
      "  WETWJCDKMRHUPV-UHFFFAOYSA-N\n\n"

      "Note that while a molecule with a particular InChI will always give the\n"
      "same InChIKey, the reverse is not true; there may exist more than one\n"
      "molecule which have different InChIs but yield the same InChIKey.\n";
  }
  virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);
  virtual unsigned int Flags() { return NOTREADABLE;};
};

}//namespace OpenBabel

#endif


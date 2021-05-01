/**********************************************************************
unique.cpp - A OBOp for eliminating chemically identical molecules during conversion.

Copyright (C) 2009 by Chris Morley

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
#include <openbabel/op.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/descriptor.h>
#include <openbabel/inchiformat.h>
#if defined(_MSC_VER) || defined(_LIBCPP_VERSION)
  #include <unordered_map>
#elif (__GNUC__ == 4 && __GNUC_MINOR__ >= 1 && !defined(__APPLE_CC__))
  #include <tr1/unordered_map>
#else
  #ifdef USE_BOOST
    #include <boost/tr1/unordered_map.hpp>
  #else
    #define NO_UNORDERED_MAP
    #include <map>
  #endif
#endif

using namespace std;
#ifndef NO_UNORDERED_MAP
  #ifdef _LIBCPP_VERSION
    using std::unordered_map;
  #else
    using std::tr1::unordered_map;
  #endif
#endif
namespace OpenBabel
{

class OpUnique : public OBOp
{
public:
  OpUnique(const char* ID) : OBOp(ID, false){
    OBConversion::RegisterOptionParam("unique", nullptr, 1, OBConversion::GENOPTIONS);
  }

  const char* Description(){ return
    "[param] remove duplicates by descriptor;default inchi\n"
    "param is a descriptor or property, or a truncation spec for InChI\n"
    "(making the comparison less detailed, see below).\n"
    "An OpenBabel warning message is output for each duplicate.\n"
    "Examples: --unique   --unique cansmi   --unique /nostereo\n\n"

    "The duplicates can be output instead by making the first character\n"
    "in the parameter ~  e.g. --unique ~cansmi   --unique ~\n\n"

    "/formula  formula only\n"
    "/connect  formula and connectivity only\n"
    "/nostereo ignore E/Z and sp3 stereochemistry\n"
    "/nosp3    ignore sp3 stereochemistry\n"
    "/noEZ     ignore E/Z steroeochemistry\n"
    "/nochg    ignore charge and protonation\n"
    "/noiso    ignore isotopes\n\n"
; }

  virtual bool WorksWith(OBBase* pOb) const { return dynamic_cast<OBMol*>(pOb) != nullptr; }
  virtual bool Do(OBBase* pOb, const char* OptionText, OpMap* pmap, OBConversion* pConv);

private:

  bool _reportDup;
  std::string _trunc;
  OBDescriptor* _pDesc;
  unsigned _ndups;
  bool _inv;

#ifdef NO_UNORDERED_MAP
  typedef map<std::string, std::string> UMap;
#else
  typedef unordered_map<std::string, std::string> UMap;
#endif

  //key is descriptor text(usually inchi) value is molecule title
  UMap _inchimap;
};

/////////////////////////////////////////////////////////////////
OpUnique theOpUnique("unique"); //Global instance

/////////////////////////////////////////////////////////////////
bool OpUnique::Do(OBBase* pOb, const char* OptionText, OpMap* pmap, OBConversion* pConv)
{
  OBMol* pmol = dynamic_cast<OBMol*>(pOb);
  if(!pmol)
    return false;

  if(pConv->IsFirstInput())
  {
    _ndups=0;
    string descID("inchi"); // the default
    _trunc.clear();
    _inv = OptionText[0]=='~';   //has the parameter a leading ~ ?
    if(_inv)
      clog << "The output has the duplicate structures" << endl;

    if(OptionText[0+_inv]=='/')  //is parameter is /x?
      _trunc = OptionText+_inv;
    else if(OptionText[0+_inv]!='\0') // not empty?
      descID = OptionText+_inv;

    _pDesc = OBDescriptor::FindType(descID.c_str());
    if(!_pDesc)
    {
      obErrorLog.ThrowError(__FUNCTION__,
              "Cannot find descriptor " + descID, obError, onceOnly);
      return false;
    }
    _pDesc->Init();
    _inchimap.clear();

    _reportDup = !_inv; //do not report duplicates when they are the output
  }

  if(!_pDesc)
    return false;
  std::string s;
  _pDesc->GetStringValue(pmol, s);

  if(!_trunc.empty())
    InChIFormat::EditInchi(s, _trunc);
  std::pair<UMap::iterator, bool> result = _inchimap.insert(make_pair(s, pmol->GetTitle()));
  bool ret = true;
  if(!s.empty() && !result.second)
  {
    // InChI is already present in set
    ++_ndups;
    if(_reportDup)
      clog << "Removed " << pmol->GetTitle() << " - a duplicate of " << result.first->second
         << " (#" << _ndups << ")" << endl;
    //delete pOb;
    ret = false; //filtered out
  }
  if(_inv)
    ret = !ret;
  if(!ret)
    delete pOb;
  return ret;
}


}//namespace
/*
Usage: --unique param
During conversion, this option eliminates molecules that are identical by some criterion.
With current babel interface it needs to be last on the command line. With nbabel
it can be anywhere.
If param is missing the criterion is the InChI.
If param starts with / the criterion is a truncated InChI string, see below.
Otherwise param is taken as a descriptor or property ID and the criterion is
its string value. Descriptors which couldbe useful here are cansmi, cansmiNS
(ignores stereo) and possibly title.

OpUnique works by attempting to insert the string value of the descriptor for
each molecule to an internal std::unordered_map. If the string has been seen
previously, the molecule is deleted and OpUnique::Do() returns false, which
causes the molecule not to be output.

InChI trucation values. param can be a concatination of these e.g. /nochg/noiso
/formula  formula only
/connect formula and connectivity only
/nostereo ignore E/Z and sp3 stereochemistry
/nosp3    ignore sp3 stereochemistry
/noEZ     ignore E/Z steroeochemistry
/nochg    ignore charge and protonation
/noiso    ignore isotopes

*/

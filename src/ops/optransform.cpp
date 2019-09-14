/**********************************************************************
optransform.cpp: makes option to transform molecule as specified in a datafile
Copyright (C) 2008 Chris Morley

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
#include <openbabel/optransform.h>
#include <openbabel/mol.h>
#include <openbabel/oberror.h>
#include <openbabel/locale.h>
#include <vector>

using namespace std;
namespace OpenBabel
{
  const char* OpTransform::Description()
  {
    //Adds name of datafile containing SMARTS strings to the description
    static std::string txt;
    txt =  _descr;
    txt += "\n Datafile: ";
    txt += _filename;
    txt += "\nOpTransform is definable";
    return txt.c_str();
  }

bool OpTransform::Initialize()
{
  _dataLoaded=true;
  _transforms.clear();
  ifstream ifs;
  if(ifs.is_open())
    ifs.close();
  char charBuffer[BUFF_SIZE];

  // Set the locale for number parsing to avoid locale issues: PR#1785463
  obLocale.SetLocale();

  if(strcmp(_filename,"*"))
  {
    if(!strncmp(_filename,"TRANSFORM",9))//A single transform can replace the filename
    {
      ParseLine(_filename);
      return true;
    }
    OpenDatafile(ifs, _filename);
    if(!ifs)
    {
      obErrorLog.ThrowError(__FUNCTION__," Could not open " + string(_filename), obError);
      return false;
    }
    while(ifs.getline(charBuffer,BUFF_SIZE))
      ParseLine(charBuffer);
  }
  else //When filename is * use data in lines following
    for(unsigned int i=4; i < _textlines.size(); ++i)
      ParseLine(_textlines[i].c_str());


  // return the locale to the original one
  obLocale.RestoreLocale();

  return true;
}

///////////////////////////////////////////////////
void OpTransform::ParseLine(const char *buffer)
{
  vector<string> vs;

  if (buffer[0] == '#')
    return;

  if (EQn(buffer,"TRANSFORM",7))
  {
    //Split into TRANSFORM reactantSMARTS ProductSMARTS
    tokenize(vs, buffer, " >\t\n");
    OBChemTsfm tr;
    if (vs.empty() || vs.size() < 3 || vs[1].empty() || vs[2].empty())
    {
      string mes("Could not parse line:\n");
      obErrorLog.ThrowError(__FUNCTION__, mes + buffer, obWarning);
    }
    else
    {
      if(!tr.Init(vs[1],vs[2]))
      {
        string mes("Could not make valid transform from the line:\n");
        obErrorLog.ThrowError(__FUNCTION__, mes + buffer, obWarning);
      }
      else
        _transforms.push_back(tr);
    }
  }
}
bool OpTransform::Do(OBBase* pOb, const char* OptionText, OpMap* pOptions, OBConversion* pConv)
{
  OBMol* pmol = dynamic_cast<OBMol*>(pOb);
  if(!pmol)
    return false;

#ifndef _DEBUG //reloads on every call in debug mode
  if(!_dataLoaded)
#endif
    if(!Initialize())
      return false;

  vector<OBChemTsfm>::iterator itr;
  for(itr=_transforms.begin(); itr!=_transforms.end();++itr)
    itr->Apply(*pmol);
  return true;
}
//Dummy instance so that OpTransform will be recognized by define plugin
//and more useful instances created from a textfile
OpTransform dummy("_","","OpTransform Dummy");

/*
////////////////////////////////////////////////////////////////////////////
//A global instance
OpTransform theTautomers("tautomers",    //commandline option is --tautomers
                         "tautomers.txt", //datafile with transforms
                         "Replaces tautomers with their standard forms");

////////////////////////////////////////////////////////////////////////////
*/

}//namespace


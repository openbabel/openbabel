/**********************************************************************
largest.cpp - An OBOp to select N molecules by descriptor value. 

Copyright(C) 2011 by Chris Morley

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
#include <openbabel/descriptor.h>
#include <openbabel/obconversion.h>
#include "deferred.h"
#include <sstream>
#include <algorithm>

using namespace std;

namespace OpenBabel
{

class OpLargest : public OBOp
{
public:
  OpLargest(const char* ID) : OBOp(ID, false){};
  const char* Description()
  {
    //Need to use a member variable so that const char* is valid when it is returned
    description = (strcmp(GetID(),"largest")!=0) ?
     "# <descr> Output # mols with smallest values of descriptor(not displayed in GUI)\n"
     "    obabel infile.xxx -Ooutfile.yyy --smallest 5 MW\n"
     "will convert only the molecules with the 5 smallest molecular weights.\n" :
     "# <descr> Output # mols with largest values\n"
     "of a descriptor <descr>. For example:\n"
     "    obabel infile.xxx -Ooutfile.yyy --largest 5 MW\n"
     "will convert only the molecules with the 5 largest molecular weights.\n";
    description +=
     "A property (OBPairData) can be used instead of a descriptor, but\n"
     "must be present in the first molecule. If the number is omitted,\n"
     "1 is assumed.\n"
     "The parameters can be in either order.\n"
     "Preceding the descriptor by ~ inverts the comparison. (Use this form in the GUI.)\n"
     "If a + follows the descriptor, e.g. MW+ , the value will be added to the title.\n";
    return description.c_str();
  }

  virtual bool WorksWith(OBBase* pOb)const{ return dynamic_cast<OBMol*>(pOb)!=NULL; }
  virtual bool Do(OBBase* pOb, const char* OptionText=NULL, OpMap* pOptions=NULL, OBConversion* pConv=NULL);
  virtual bool ProcessVec(vector<OBBase*>& vec);
  static bool MatchPairData(OBBase* pOb, std::string& s); //Copy of protected OBDescriptor function

private:
  std::string description;
  std::multimap<double, OBBase*> _selmap;
  OBDescriptor* _pDesc;
  std::string _param;
  std::string _prop;
  bool _addDescToTitle;
  bool _rev;
  unsigned _nmols;
  OBConversion* _pConv;
};

/////////////////////////////////////////////////////////////////
OpLargest theOpLargest("largest"); //Global instances
OpLargest theOpSmallest("smallest"); 
/////////////////////////////////////////////////////////////////
bool OpLargest::Do(OBBase* pOb, const char* OptionText, OpMap* pOptions, OBConversion* pConv)
{
  if(!strcmp(OptionText, "inactive"))
    return true;

  OBMol* pmol = dynamic_cast<OBMol*>(pOb);
  if(!pmol)
    return false;
  if(pConv->IsFirstInput())
  {
    _pConv = pConv;
    _selmap.clear();
    _rev = strcmp(GetID(),"largest")!=0; //_rev is initially true for --smallest
    std::vector<string> vec;
    tokenize(vec, OptionText);
    if(vec.size()>1 && isdigit(vec[1][0]))
      swap(vec[0], vec[1]); //ensure n desc or desc 

    unsigned idesc = vec.size()>1 ? 1 : 0;

    if(vec[idesc][0]=='~')
    {
      _rev = !_rev;
      vec[idesc].erase(0,1);
    }
    //If + is last char of descriptor the val is to be added to title on output
    string::iterator it = vec[idesc].end() -1;
    _addDescToTitle = (*it=='+');
    if(_addDescToTitle)
      vec[idesc].erase(it);

    //If the first molecule has matching OBPairData, use that and set _pDesc to NULL
    if(!vec.empty() && MatchPairData(pOb, vec[idesc]))
    {
      _pDesc = NULL;
      _prop = vec[idesc];
    }
    else
    {
      //Read  a descriptor name and parameter 
      istringstream ss;
      ss.str(vec[idesc]);
      pair<string,string> spair = OBDescriptor::GetIdentifier(ss);
      if(vec.empty() || !(_pDesc = OBDescriptor::FindType(spair.first.c_str())))
      {
        obErrorLog.ThrowError(__FUNCTION__,
          "Property or descriptor was not recognized.\n", obError, onceOnly);
        pConv->SetOneObjectOnly(); //stop all conversion
        return false;
      }
      _param = spair.second.empty() ? "" : spair.second;
    }
    
    _nmols=1;
    if(vec.size()>1)
      _nmols = atoi(vec[0].c_str());
    if(_nmols==0)
      _nmols=1;
 
    //It is necessary for this op to act *after* all the other options. They
    //might filter out molecules and the user will expect this op to select
    //from the molecules that are left. By default, ops act on an object
    // *before*, for example, --filter, so that the following trick is used. 
    //This op is removed so that it is not called again normally. A Deferred
    //output format is made, which, because of the true CallDo parameter in the
    //constructor, will immediately call this function,when it receives objects.
    //The handling of the selected molecules is all done here and returning
    //false means that DeferredFormat stores nothing. At the end it calls 
    //this op's ProcessVec() function.

    pConv->AddOption(GetID(),OBConversion::GENOPTIONS,"inactive");//removing messes up DoOps()
    new DeferredFormat(pConv, this, true); //it will delete itself
    //pConv->AddOption("OutputAtEnd",OBConversion::GENOPTIONS);
    return true;
  }

  // All molecules (called from DeferredFormat)
  
  //Save in map if descriptor val is better
  // than the current selection or otherwise delete
  double val;
  if(_pDesc)
    val = _pDesc->Predict(pOb, &_param);
  else
  {
    stringstream ss(pOb->GetData(_prop)->GetValue());
    ss >> val;
  }
 
  if(_selmap.size()<_nmols)
    //populate map of selected mols up to desired number
    _selmap.insert(make_pair(val, pOb));
  else
  {
    //replace mols in selection if new mol is better candidate
    multimap<double, OBBase*>::iterator leastwanted = 
      _rev ? --_selmap.end() : _selmap.begin();
    if((!_rev && val>leastwanted->first) || (_rev && val<leastwanted->first))
    {
      //have a better candidate
      delete leastwanted->second; // delete the worst molecule
      _selmap.erase(leastwanted); // remove it from the map
      _selmap.insert(make_pair(val, pOb)); //add new candidate
    }
    else
      delete pOb; // discard  a mol that did not make selection
  }
  return false; //do not save in DeferredFormat 
}

bool OpLargest::ProcessVec(vector<OBBase*>& vec)
{
  //Called at the end.
  //Add the selected mols to the vec for Deferred format to output
  vec.clear();
  vec.reserve(_selmap.size());
  multimap<double, OBBase*>::reverse_iterator iter;
  for(iter=_selmap.rbegin(); iter!=_selmap.rend(); ++iter)
  {
    if(_addDescToTitle)
    {
      std::stringstream ss;
      ss << iter->second->GetTitle() << ' ' << iter->first;
      iter->second->SetTitle(ss.str().c_str());
    }
    vec.push_back(iter->second);
  }
  if(_rev) //it's difficult to find an elegant and efficient way of doing this
    reverse(vec.begin(),vec.end());
  return true;
}

bool OpLargest::MatchPairData(OBBase* pOb, string& s)
{
  //If s matches a PairData attribute return true
  //else if s with all '_' replaced by spaces matches return true and s is now the form with spaces
  //else return false.
  if(pOb->HasData(s))
    return true;
  if(s.find('_')==string::npos)
    return false;
  string temp(s);
  string::size_type pos = string::npos;
  //replace all underscores by spaces
  while((pos=temp.find('_', ++pos))!=string::npos)
    temp[pos]=' ';
  if(pOb->HasData(temp))
  {
    s = temp;
    return true;
  }
  return false;
}

}//namespace

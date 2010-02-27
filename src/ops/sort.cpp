/**********************************************************************
sort.cpp - A OBOp for sorting molecules during conversion.

Copyright (C) 2009 by Chris Morley
 
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
#include <openbabel/op.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/descriptor.h>
#include <set>

namespace OpenBabel
{

class DeferredFormat : public OBFormat
{
public:
  DeferredFormat(OBConversion* pConv, OBOp* pOp=NULL)
  {
    _pRealOutFormat = pConv->GetOutFormat();
    pConv->SetOutFormat(this);
    _pOp = pOp;
  }
  virtual const char* Description() { return "Read and write an OBBase* array"; }
  virtual bool ReadChemObject(OBConversion* pConv)
  {
    if(_obvec.empty())
    {
      delete this;//self destruction; was made in new in an OBOp
      return false;
    }
    //returned in reverse order, because its easier with a vector
    pConv->AddChemObject(_obvec.back());
    _obvec.pop_back();
    return true;
  }

  virtual bool WriteChemObject(OBConversion* pConv)
  {
    //Store the object pointer.
    //Unlike most formats, no deletion of object here or object constuction in ReadChemObject.
    _obvec.push_back(pConv->GetChemObject());
    if(pConv->IsLast())
    {
      //At the end, sort, or whatever, the vector
      if(_pOp)
      {
        if(_pOp->ProcessVec(_obvec))
        {
          //Now output the processed vector
          pConv->SetInAndOutFormats(this, _pRealOutFormat);

          pConv->SetInStream(&std::ifstream());//Not used, but Convert checks it is ok
          pConv->GetInStream()->clear();

          //clear the options - they have already been applied
          pConv->SetOptions("",OBConversion::GENOPTIONS);
          pConv->SetOutputIndex(0);
          pConv->Convert();
        }
      }
    }
    return true;
  }
private:
  OBFormat* _pRealOutFormat;
  std::vector<OBBase*> _obvec;
  public:
  OBOp* _pOp;
};

//*****************************************************************
template<class T>
struct Order : public std::binary_function<std::pair<OBBase*,T>, std::pair<OBBase*,T>, bool>
{
  Order(OBDescriptor* pDesc, bool rev) : _pDesc(pDesc), _rev(rev){}
  bool operator()(std::pair<OBBase*,T> p1, std::pair<OBBase*,T> p2) const
  {
    return _rev ? 
      _pDesc->Order(p1.second, p2.second) : 
      _pDesc->Order(p2.second, p1.second);
  }
  OBDescriptor* _pDesc;
  bool _rev;
};
//*****************************************************************
class OpSort : public OBOp
{
public:
  OpSort(const char* ID) : OBOp(ID, false)
  {
    OBConversion::RegisterOptionParam(ID, NULL, 1, OBConversion::GENOPTIONS);
  }
  
  const char* Description(){ return "<desc> Sort by descriptor(~desc for reverse)"
    "\n Follow descriptor with + to also add it to the title, e.g. MW+ "
    "\n Custom ordering is possible; see inchi descriptor"; }

  virtual bool WorksWith(OBBase* pOb)const{ return dynamic_cast<OBMol*>(pOb)!=NULL; }
  virtual bool Do(OBBase* pOb, const char* OptionText, OpMap* pmap, OBConversion* pConv);
  virtual bool ProcessVec(std::vector<OBBase*>& vec);
private:
  OBDescriptor* _pDesc;
  bool _rev;
  bool _addDescToTitle;
};

/////////////////////////////////////////////////////////////////
OpSort theOpSort("sort"); //Global instance

/////////////////////////////////////////////////////////////////
bool OpSort::Do(OBBase* pOb, const char* OptionText, OpMap* pmap, OBConversion* pConv)
{
  if(pConv && pConv->IsFirstInput())
  {
    _rev=false;
    if(*OptionText=='~')
    {
      _rev=true;
      ++OptionText;
    }

    const char* pLast = OptionText + strlen(OptionText)-1;
    _addDescToTitle = *(OptionText + strlen(OptionText)-1)=='+';//last char
    if(_addDescToTitle)
      *const_cast<char*>(pLast)='\0';

    _pDesc = OBDescriptor::FindType(OptionText);
    if(!_pDesc)
    {
     obErrorLog.ThrowError(__FUNCTION__, 
              std::string("Unknown descriptor ") + OptionText, obError, onceOnly);
     return false;
    }
    _pDesc->Init();//needed  to clear cache of InChIFilter


    //Make a deferred format and divert the output to it
    new DeferredFormat(pConv, this); //it will delete itself
  }
  return true;
}

//****************************************************************
bool OpSort::ProcessVec(std::vector<OBBase*>& vec)
{
  // Make a vector containing both the OBBase* and the descriptor value and the sort it
  if(!IsNan(_pDesc->Predict(vec[0])))
  {
    //a numerical descriptor
    //Copy into a pair vector
    std::vector<std::pair<OBBase*,double> > valvec;
    valvec.reserve(vec.size());
    std::vector<OBBase*>::iterator iter;
    for(iter=vec.begin();iter!=vec.end();++iter)
      valvec.push_back(std::make_pair<OBBase*,double>(*iter, _pDesc->Predict(*iter)));

    //Sort
    std::sort(valvec.begin(),valvec.end(), Order<double>(_pDesc, _rev));

    //Copy back
    std::vector<std::pair<OBBase*,double> >::iterator valiter;
    iter=vec.begin();
    for(valiter=valvec.begin();valiter!=valvec.end();++valiter, ++iter)
    {
      *iter = valiter->first;
      if(_addDescToTitle)
      {
        std::stringstream ss;
        ss << (*iter)->GetTitle() << ' ' << valiter->second;
        (*iter)->SetTitle(ss.str().c_str());
      }
    }
  }
  else
  {
    //a string descriptor
    //Copy into a pair vector
    std::vector<std::pair<OBBase*,std::string> > valvec;
    valvec.reserve(vec.size());
    std::vector<OBBase*>::iterator iter;
    std::string s;
    for(iter=vec.begin();iter!=vec.end();++iter)
    {
      _pDesc->GetStringValue(*iter, s);
      valvec.push_back(std::make_pair<OBBase*,std::string>(*iter, s));
    }

    //Sort
    std::sort(valvec.begin(),valvec.end(), Order<std::string>(_pDesc, _rev));

    //Copy back
    std::vector<std::pair<OBBase*,std::string> >::iterator valiter;
    iter=vec.begin();
    for(valiter=valvec.begin();valiter!=valvec.end();++valiter, ++iter)
    {
      *iter = valiter->first;
      if(_addDescToTitle)
      {
        std::stringstream ss;
        ss << (*iter)->GetTitle() << ' ' << valiter->second;
        (*iter)->SetTitle(ss.str().c_str());
      }
    }
  }

  return true;
}
/*
This started as a nice compact piece of code! The need to handle descriptors
which return either numbers or strings was originally achieved without testing
the type here by using LessThan() in the descriptor. But this meant that the
descriptor value was recalculated every time it was needed, which is inappropriate
for sorting. A local cache of InChI values was implemented but a more general
solution was needed. The values are now calculated once and stored here in a
vector, which stores numbers or strings and the code is extensively duplicated
because of this. But using templates was not much shorter because four templated
functions were needed, and the code more difficult to understand.
*/
}//namespace

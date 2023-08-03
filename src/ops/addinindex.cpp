/**********************************************************************
addinindex.cpp - Adds input index to title.

Copyright(C) 2007 by Chris Morley
 
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
#include<openbabel/base.h>
#include <sstream>
#include<openbabel/op.h>
#include<openbabel/obconversion.h>

namespace OpenBabel
{

class OpAddInIndex : public OBOp
{
public:
  OpAddInIndex(const char* ID) : OBOp(ID, false){};
  const char* Description(){ return 
    "Append input index to title\n"
    "These are objects before filtering. Use AddOutIndex for objects after filtering\n"; }

  virtual bool WorksWith(OBBase* pOb)const{ return true; } //all objects
  virtual bool Do(OBBase* pOb, const char*, OpMap*, OBConversion* pConv=nullptr);
};

/////////////////////////////////////////////////////////////////
OpAddInIndex theOpAddInIndex("AddInIndex"); //Global instance

/////////////////////////////////////////////////////////////////
bool OpAddInIndex::Do(OBBase* pOb, const char* OptionText, OpMap* pOptions, OBConversion* pConv)
{
  int count = pConv->GetCount();
  if(count>=0) // add nothing unless Convert interface of OBConversion is being used
  {
    std::stringstream ss;
    ss << pOb->GetTitle() << ' ' << count+1;
    pOb->SetTitle(ss.str().c_str());
  }
  return true;
}
}//namespace

/**********************************************************************
opaddfilename.cpp - Add the input file name to the title of the object

Copyright (C) 2011 by Chris Morley

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
#include <string>
#include <openbabel/babelconfig.h>
#include <openbabel/op.h>
#include <openbabel/obconversion.h>
#include <openbabel/base.h>

namespace OpenBabel
{

class OpAddFileName : public OBOp
{
public:
  OpAddFileName(const char* ID) : OBOp(ID, false){};
  const char* Description(){ return 
    "Append input filename to title\n"
    "Any path is removed from the filename\n"
    ; }

  virtual bool WorksWith(OBBase* pOb)const{ return true; } //all OBBase objects
  virtual bool Do(OBBase* pOb, const char* OptionText=NULL, OpMap* pOptions=NULL, OBConversion* pConv=NULL);
};

/////////////////////////////////////////////////////////////////
OpAddFileName theOpAddFileName("addfilename"); //Global instance

/////////////////////////////////////////////////////////////////
bool OpAddFileName::Do(OBBase* pOb, const char* OptionText, OpMap* pOptions, OBConversion* pConv)
{
  if(!pConv)
    return true; //do not stop any conversion but do nothing
  std::string  fname(pConv->GetInFilename());

  //remove path from the filename
  std::string::size_type pos = fname.find_last_of("/\\:");
  if(pos!=std::string::npos)
    fname.erase(0, pos+1);
  fname = " " + fname;
  pOb->SetTitle((pOb->GetTitle() + fname).c_str());
  return true;
}
}//namespace

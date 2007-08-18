/**********************************************************************
Copyright (C) 2005 by Chris Morley
 
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/
#include <openbabel/babelconfig.h>
#include <openbabel/obconversion.h>

using namespace std;
namespace OpenBabel
{

class CopyFormat : public OBFormat
{
public:
  CopyFormat()
  {
    OBConversion::RegisterFormat("copy",this);
  }

  virtual const char* Description() //required
  {
    return
"Copies raw text\n \
Objects can be chemically filtered without the risk\n \
of losing any additional information they contain,\n \
since no format conversion is done.\n \
Note that XML files may be missing non-object elements\n \
at the start or end and so may no longer be well formed.";
  };

  virtual unsigned int Flags()
  {
      return NOTREADABLE;
  };

  /////////////////////////////////////////////////////////////////
  virtual bool WriteChemObject(OBConversion* pConv)
  {
    pConv->GetChemObject();//needed to increment pConv->Index

    istream& ifs = *pConv->GetInStream();
    ostream& ofs = *pConv->GetOutStream();

    streampos startpos = pConv->GetInPos();
    int len = pConv->GetInLen();
    if(len>0)
    {
      streampos curpos = ifs.tellg();
      if(ifs.eof())
        ifs.clear();
      ifs.seekg(startpos);
      
      char* buf = new char[len];
      ifs.read(buf,len);
      ofs.write(buf,len);
      delete[] buf;

      ifs.seekg(curpos);
    }
    else
    {
      //When no length recorded, copy the whole input stream
      //Seem to need to treat stringstreams differently
      stringstream* pss = dynamic_cast<stringstream*>(&ifs);
      if(pss)
        ofs << pss->str() << flush;
      else
        ofs << ifs.rdbuf() << flush;
    }
    return true;
  };

};  

//Make an instance of the format class
CopyFormat theCopyFormat;


} //namespace OpenBabel


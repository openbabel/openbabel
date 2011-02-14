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
"Copy raw text\n"
"A utility format for exactly copying the text of a chemical file format\n"
"This format allows you to filter molecules from multimolecule files\n"
"without the risk of losing any additional information they contain,\n"
"since no format conversion is carried out.\n\n"

".. warning::\n\n"
" Currently not working correctly for files with Windows line endings.\n\n"

"Example:\n\n"

"  Extract only structures that include at least one aromatic carbon\n"
"  (by matching the SMARTS pattern ``[c]``)::\n\n"

"   babel -s '[c]' database.sdf -ocopy new.sd\n\n"

".. note::\n\n"
" XML files may be missing non-object elements\n"
" at the start or end and so may no longer be well formed.\n\n"
;
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

      char* buf = new char[len+1];
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


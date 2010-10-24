/**********************************************************************
obdoc - Automatically generate documentation database for file formats

Copyright (C) 2005 by Geoffrey R. Hutchison
 
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


// used to set import/export for Cygwin DLLs
#ifdef WIN32
#define USING_OBDLL
#endif

#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>

using namespace std;
using namespace OpenBabel;

int main(int argc,char **argv)
{
  char *program_name= argv[0];

  vector< vector<string> > extensionList; // some formats have multiple codes
  vector< OBFormat *> formatList;

  OBConversion conv;
  Formatpos pos;
  OBFormat* pFormat;
  const char* str=NULL;
  bool newFormat;

  while(OBConversion::GetNextFormat(pos,str,pFormat))
    {
      if((pFormat->Flags() & NOTWRITABLE) && (pFormat->Flags() & NOTREADABLE))
        continue;

      if (formatList.size() == 0)
        {
          formatList.push_back(pFormat);
          vector<string> formatExtensions;
          formatExtensions.push_back(pos->first);
          extensionList.push_back(formatExtensions);
        }
      else
        {
          newFormat = true; // until we find otherwise
          for(unsigned int i = 0; i < formatList.size(); ++i)
            if (formatList[i] == pFormat)
              {
                newFormat = false;
                extensionList[i].push_back(pos->first);
                break;
              }
	    
          if (newFormat)
            {
              formatList.push_back(pFormat);
              vector<string> formatExtensions;
              formatExtensions.push_back(pos->first);
              extensionList.push_back(formatExtensions);
            }
        }

    }

  string notes, description, name, code;
  ofstream ofs, indexExt, indexName;

  indexExt.open("code-index.html");

  for(unsigned int i = 0; i < formatList.size(); ++i)
    {
      pFormat = formatList[i];
      description = pFormat->Description();
      name = description.substr(0,description.find('\n'));
      code = extensionList[i][0] + "-data.html";

      ofs.open(code.c_str());

      ofs << "{{Format|" << endl;
      ofs << "|extensions=";

      for (unsigned int j = 0; j < extensionList[i].size(); ++j)
        {
          ofs << extensionList[i][j];

          indexExt << "* " << extensionList[i][j]	 
                   << " - [[" << name << "]]"	 
                   << endl;
 
          if (j != extensionList[i].size() - 1)
            ofs << ", ";
        }
      ofs << endl;

      ofs << "|mime=";
      if (pFormat->GetMIMEType() && strlen(pFormat->GetMIMEType()) != 0)
        ofs << pFormat->GetMIMEType() << endl;
      else
        ofs << endl;
	
      ofs << "|url=";
      if (pFormat->SpecificationURL() && 
          strlen(pFormat->SpecificationURL()) != 0)
        {
          ofs << pFormat->SpecificationURL() << endl;
        }
      else
        ofs << "Unknown" << endl;

      ofs << "|import=";
      if (pFormat->Flags() & NOTREADABLE)
        ofs << "No";
      else
        ofs << "Yes";
      if (pFormat->Flags() & READONEONLY)
        ofs << ", One record per input file";
      if (pFormat->Flags() & READBINARY)
        ofs << ", Input is a binary file";
      ofs << endl;

      ofs << "|export=";
      if (pFormat->Flags() & NOTWRITABLE)
        ofs << "No";
      else
        ofs << "Yes";
      if (pFormat->Flags() & WRITEONEONLY)
        ofs << ", One record per output file";
      if (pFormat->Flags() & WRITEBINARY)
        ofs << ", Output is a binary file";
      ofs << endl;

      ofs << "|version=2.2 and later" << endl; // added in 2.1
      ofs << "|dimensionality=3D" << endl;
      ofs << "|options =\n<pre>";
      ofs << 	description.substr(description.find('\n'), description.size()); 
      ofs << "</pre>" << endl;
      ofs << "}}" << endl << endl;

      ofs << "[[Category:Formats]]" << endl;

      ofs.close();
    }

  indexExt.close();

  return(0);
} // end main



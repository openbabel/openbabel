/**********************************************************************
obdoc - Automatically generate documentation database for file formats

Copyright (C) 2005 by Geoffrey R. Hutchison
 
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
  indexName.open("name-index.html");
  for(unsigned int i = 0; i < formatList.size(); ++i)
    {
      pFormat = formatList[i];
      description = pFormat->Description();
      name = description.substr(0,description.find('\n'));
      code = extensionList[i][0] + "-data.html";

      indexName << "<li>" << name
                << "<a href=\"" << extensionList[i][0] << ".shtml\">" 
                << name << "</a></li>" << endl;

      ofs.open(code.c_str());

      ofs << "<h1>" << name << "<h1>" << endl;
      ofs << "\n\n" << "<table>" << endl;
      ofs << "<tr><th>Filename Extensions</th><td>";
      for (unsigned int j = 0; j < extensionList[i].size(); ++j)
        {
          ofs << extensionList[i][j];
          indexExt << "<li>" << extensionList[i][j] 
                   << "<a href=\"" << extensionList[i][0] << ".shtml\">" 
                   << extensionList[i][j] << " - " 
                   << name << "</a></li>" << endl;

          if (j != extensionList[i].size() - 1)
            ofs << ", ";
        }
      ofs << "</td></tr>" << endl;

      ofs << "<tr><th>Chemical MIME Type</th><td>";
      if (strlen(pFormat->GetMIMEType()) != 0)
        ofs << pFormat->GetMIMEType() << "</td></tr>" << endl;
      else
        ofs << "Undefined</td></tr>" << endl;
	

      ofs << "<tr><th>Specification URL</th><td>";
      if (strlen(pFormat->SpecificationURL()) != 0)
        {
          ofs << "<a href=\"" << pFormat->SpecificationURL() << "\">";
          ofs << pFormat->SpecificationURL() << "</a>";
        }
      else
        ofs << "Unknown";
      ofs << "</td></tr>" << endl;

      ofs << "<tr><th>Notes</th><td>";
      notes.clear();
      if (pFormat->Flags() & NOTREADABLE)
        notes += "Export only";
      else if (pFormat->Flags() & NOTWRITABLE)
        notes += "Import only";

      if (pFormat->Flags() & READONEONLY)
        {
          if (notes.size())
            notes += ", ";
          notes += "One record per input file";
        }
      if (pFormat->Flags() & WRITEONEONLY)
        {
          if (notes.size())
            notes += ", ";
          notes += "One record per output file";
        }
      if (pFormat->Flags() & READBINARY)
        {
          if (notes.size())
            notes += ", ";
          notes += "Input is a binary file";
        }
      if (pFormat->Flags() & WRITEBINARY)
        {
          if (notes.size())
            notes += ", ";
          notes += "Output is a binary file";
        }


      if (notes.size() == 0)
        notes = "None";
      ofs << notes << "</td></tr>" << endl;

      ofs << "</table>" << endl;


      ofs << "<h2>Options</h2>\n<pre>";
      ofs << 	description.substr(description.find('\n'), description.size()); 
      ofs << "\n</pre>" << endl;

      ofs.close();
    }

  indexExt.close();
  indexName.close();
  return(0);
} // end main



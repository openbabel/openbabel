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

#include "mol.h"
#include "obconversion.h"

using namespace std;
using namespace OpenBabel;

int main(int argc,char **argv)
{
    char *program_name= argv[0];

    vector< vector<string> > extensionList; // some formats have multiple codes
    vector< OBFormat *> formatList;

    Formatpos pos;
    OBFormat* pFormat;
    OBConversion conv;
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
	    for(unsigned int i = 0; i < formatList.size(); i++)
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

    string impt, expt, description, name, code;
    ofstream ofs, indexExt, indexName;

    indexExt.open("code-index.html");
    indexName.open("name-index.html");
    for(unsigned int i = 0; i < formatList.size(); i++)
      {
	pFormat = formatList[i];
	description = pFormat->Description();
	name = description.substr(0,description.find('\n'));
	code = extensionList[i][0] + "-data.html";

	indexName << "<li>" << name
		 << "<a href=\"" << extensionList[i][0] << ".shtml\">" 
		 << name << "</a></li>" << endl;

	ofs.open(code.c_str());

	ofs << name << endl << endl;
	ofs << "{{Format|" << endl;
	ofs << " |extensions=";
	for (unsigned int j = 0; j < extensionList[i].size(); j++)
	  {
	    ofs << extensionList[i][j];
	    indexExt << "<li>" << extensionList[i][j] 
		     << "<a href=\"" << extensionList[i][0] << ".shtml\">" 
		     << extensionList[i][j] << " - " 
		     << name << "</a></li>" << endl;

	    if (j != extensionList[i].size() - 1)
	      ofs << ", ";
	  }
	ofs << endl;

	ofs << " |mime=";
	if (strlen(pFormat->GetMIMEType()) != 0)
	  ofs << pFormat->GetMIMEType() << endl;
	else
	  ofs << "Undefined" << endl;

	ofs << " |url=";
	if (strlen(pFormat->SpecificationURL()) != 0)
	  {
	    ofs << pFormat->SpecificationURL();
	  }
	else
	  ofs << "Unknown";
	ofs << endl;

	impt.clear();
	expt.clear();
	if (pFormat->Flags() & NOTREADABLE)
	  impt += "No";
	else
	  impt += "Yes";

	if (pFormat->Flags() & NOTWRITABLE)
	  expt += "No";
	else
	  expt += "Yes";

	if (pFormat->Flags() & READONEONLY)
	  {
	    if (impt.size())
	      impt += ", ";
	    impt += "One record per input file";
	  }
	if (pFormat->Flags() & WRITEONEONLY)
	  {
	    if (expt.size())
	      expt += ", ";
	    expt += "One record per output file";
	  }
	if (pFormat->Flags() & READBINARY)
	  {
	    if (impt.size())
	      impt += ", ";
	    impt += "Input is a binary file";
	  }
	if (pFormat->Flags() & WRITEBINARY)
	  {
	    if (expt.size())
	      expt += ", ";
	    expt += "Output is a binary file";
	  }

	ofs << " |import=";
	ofs << impt << endl;
	ofs << " |export=";
	ofs << expt << endl;

	ofs << " |version=All" << endl;
	ofs << " |dimensionality=3D" << endl;

	ofs << " |options=" << endl;
	ofs << "<pre>";
	ofs << 	description.substr(description.find('\n'), description.size()); 
	ofs << "</pre>" << endl;
	ofs << "}}" << endl << endl;

	ofs << "[[Category:Formats]]" << endl;

	ofs.close();
      }

    indexExt.close();
    indexName.close();
    return(0);
} // end main



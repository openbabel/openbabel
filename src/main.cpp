/**********************************************************************
main.cpp - Main conversion program, command-line handling.

Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (c) 2001-2003 by Geoffrey R. Hutchison
Some portions Copyright (c) 2004 by Chris Morley

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

#include "babelconfig.h"
#if HAVE_IOSTREAM
	#include <iostream>
#elif HAVE_IOSTREAM_H
	#include <iostream.h>
#endif
#if HAVE_FSTREAM
	#include <fstream>
#elif HAVE_FSTREAM_H
	#include <fstream.h>
#endif
#if HAVE_SSTREAM
	#include <sstream>
#elif
	#include <sstream.h>
#endif

#include <string>
#include <map>
#if HAVE_CONIO_H
	#include <conio.h>
#endif

#include "obconversion.h"

using namespace std;
using namespace OpenBabel;

void usage();
void help();

// There isn't a great way to do this -- we need to save argv[0] for usage()
static char *program_name;

int main(int argc,char *argv[])
{

	OBConversion Conv(&cin, &cout); //default input and output are console 

	string GenOptions;
	OBFormat* pInFormat;
	OBFormat* pOutFormat;
	vector<string> FileList, OutputFileList;
	string OutputFileName;

  // Parse commandline
	bool gotInType = false, gotOutType = false;
	bool UseSavedOptions = false;
	bool SplitOrBatch=false;

  char *oext;
  char *iext;

  //Save name of program without its path (and .exe)
	string pn(argv[0]);
	int pos;
	#ifdef _WIN32
		pos = pn.find(".exe");
		if(pos!=string::npos)
			argv[0][pos]='\0';
	#endif
	pos = pn.find_last_of("/\\");
	if(pos==string::npos)
		program_name=argv[0];
	else
		program_name=argv[0]+pos+1;

	int arg;
	for (arg = 1; arg < argc; arg++)
  {
    if (argv[arg])
		{
			if (argv[arg][0] == '-')
	    {	      
	      switch (argv[arg][1])
				{
				case 'V':
					{
						cout << "Open Babel " << BABEL_VERSION << " -- " 
					 << __DATE__ << " -- " << __TIME__ << endl;
						exit(0);
					}
				
				case 'h':
					if (strncmp(argv[arg],"-hpH",4) == 0)
						GenOptions+='p';
					else
						GenOptions+='h';
					break;

				case 'i':
					gotInType = true;
					iext = argv[arg] + 2;

					//The ID provided by the OBFormat class is used as the identifying file extension
					//This is a slight reduction in flexibility (which is not currently used)
					pInFormat = Conv.FindFormat(iext);
					if(pInFormat==NULL)
					{
						cerr << program_name << ": cannot read input format!" << endl;
						usage();
					}
					break;

				case 'o':
					gotOutType = true;
					oext = argv[arg] + 2;
					pOutFormat = Conv.FindFormat(oext);
					if(pOutFormat==NULL)
					{
						cerr << program_name << ": cannot write output format!" << endl;
						usage();
					}
					break;
					
				case 'x':
					Conv.SetOptions(argv[arg]);
					break;

				case 'z':
					UseSavedOptions=true;
					break;

				case 'H':
					if(isalnum(argv[arg][2]))
					{
						if(strncasecmp(argv[arg]+2,"all",3))
						{
							OBFormat* pFormat = Conv.FindFormat(argv[arg]+2);
							if(pFormat)
							{
								cout << argv[arg]+2 << "  " << pFormat->Description() << endl;
								if(strlen(pFormat->SpecificationURL()))
									cout << "Specification at: " << pFormat->SpecificationURL() << endl;
							}
							else
								cout << "Format type: " << argv[arg]+2 << " was not recognized" <<endl;
						}
						else
						{
							Formatpos pos;
							OBFormat* pFormat;
							const char* str=NULL;
							while(OBConversion::GetNextFormat(pos,str,pFormat))
							{
								cout << str << endl;
								char* p = strchr(pFormat->Description(),'\n');
								cout << p+1; //second line of description
								if(strlen(pFormat->SpecificationURL()))
									cout << "Specification at: " << pFormat->SpecificationURL();
								cout << endl << endl;
							}
						}
					}
					else
						help();
					#ifdef _DEBUG
				 //CM keep window open
					cout << "Press any key to finish" <<endl;
					getch();
					#endif
					exit(0);

				case '-': //Do nothing
					/*if (inFileArg == 0)
						inFileArg = -1;
					else
						outFileArg = -1;
					*/
					break;

				case 'm': //multiple output files
					SplitOrBatch=true;
					break;
					
				default:
					char* p = argv[arg]+1;
					GenOptions+=*p++;
					//quotes seem to be stripped out;add them back if there are more chars
					if(*p)
					{
						GenOptions+='\"';
						while(*p)
							GenOptions+=*p++;
						GenOptions+='\"';
					}
					break;
				}
			}
			else 
			{
				//filenames
				if(!gotOutType)
					FileList.push_back(argv[arg]);
				else
					OutputFileName = argv[arg];
			}
		}
	}

	if(!gotOutType) //the last file is the output
	{
		if(FileList.empty())
		{
			cerr << "No output file or format spec!" << endl;
			usage();
		}
		OutputFileName = FileList.back();
		FileList.pop_back();
	}

	#ifdef _WIN32
		//Expand wildcards in input filenames and add to FileList
		vector<string> tempFileList(FileList);
		FileList.clear();
		vector<string>::iterator itr;
		for(itr=tempFileList.begin();itr!=tempFileList.end();itr++)
			DLHandler::findFiles (FileList, *itr);
	#endif

	Conv.SetGeneralOptions(GenOptions.c_str());

	if (!gotInType)
  {
		pInFormat = Conv.FormatFromExt(FileList[0].c_str());
		if(pInFormat==NULL)
		{
			cerr << program_name << ": cannot read input format!" << endl;
			usage();
		}
  }
  if (!gotOutType)
  {
		pOutFormat = Conv.FormatFromExt(OutputFileName.c_str());
		if(pOutFormat==NULL)
		{
			cerr << program_name << ": cannot write output format!" << endl;
			usage();
		}
  }

	Conv.SetInAndOutFormats(pInFormat,pOutFormat);

	if(SplitOrBatch)
	{
		//Put * into output file name
		if(OutputFileName.empty())
		{
			OutputFileName = "*.";
			OutputFileName += oext;
		}
		else
		{
			int pos = OutputFileName.find_last_of('.');
			if(pos==string::npos)
				OutputFileName += '*';
			else
				OutputFileName.insert(pos,"*");
		}
	}

	// Need to make this configurable -- on UNIX should probably be $HOME/.openbabel or similar
	#ifdef _WIN32
	const string OptFile = "options.txt";
	#else
	string OptFile = getenv("HOME");
	OptFile.append("/.openbabel");
	#endif

	if(UseSavedOptions)
		Conv.RestoreOptionsFromFile(OptFile.c_str());
	//Write current options to file
	Conv.SaveOptionsToFile(OptFile.c_str());

		
	int count = Conv.FullConvert(FileList, OutputFileName, OutputFileList);
	cout << count << " molecule(s) converted" <<endl;

	if(OutputFileList.size()>1)
	{
		cout << OutputFileList.size() << " files output. The first is " << OutputFileList[0] <<endl;
	}

#ifdef _DEBUG
 //CM keep window open
  cout << "Press any key to finish" <<endl;
	getch();
#endif

	return 0;
};

void usage()
{
  cout << "Open Babel " << BABEL_VERSION << " -- " << __DATE__ << " -- "
       << __TIME__ << endl;
  cout << "Usage: " << program_name
       << " [-i<input-type>] <name> [-o<output-type>] <name>" << endl;
  cout << "Try  -H option for more information." << endl;

	#ifdef _DEBUG
	 //CM keep window open
		cout << "Press any key to finish" <<endl;
		getch();
	#endif

  exit (0);
}

void help()
{
  cout << "Open Babel converts chemical structures from one file format to another"<< endl;
  cout << "Usage: " << program_name << " <input spec> <output spec> [Options]" << endl;
	cout << "Each spec can be a file whose extension decides the format." << endl;
	cout << "Optionally the format can be specified by preceding the file by" << endl;
	cout << "-i<format-type> e.g. -icml, for input and -o<format-type> for output" << endl;
	cout << "See below for available format-types, which are the same as the " << endl;
	cout << "file extensions and are case independent." << endl; 
	cout << "If no input or output file is given stdin or stdout are used instead." << endl << endl; 
	cout << "More than one input file can be specified and their names can contain" <<endl;
	cout << "wildcard chars(*and?).The molecules are aggregated in the output file.\n" << endl;
	cout << OBConversion::Description(); // Conversion options
  cout << "  -H Outputs this help text" << endl;
	cout << "  -Hxxx (xxx is file format ID e.g. -Hcml) gives format info" <<endl; 
	cout << "  -Hall Outputs details of all formats" <<endl; 
  cout << "  -V Outputs version number" <<endl; 
	cout << "  -m Produces multiple output files, to allow:" <<endl;
  cout << "     Splitting: e.g.        " << program_name << " infile.mol new.smi -m" <<endl;
  cout << "       puts each molecule into new1.smi new2.smi etc" <<endl;
  cout << "     Batch conversion: e.g. " << program_name << " *.mol -osmi -m" <<endl;
  cout << "       converts each input file to a .smi file" << endl;
#ifdef _WIN32
	cout << "     In Windows these can also be done using the forms" <<endl;
  cout << "       " << program_name << " infile.mol new*.smi and " << program_name << " *.mol *.smi respectively.\n" <<endl;
#endif

	cout << OBConversion::GetDefaultFormat()->TargetClassDescription();// some more options probably for OBMol
	cout << "The following file formats are recognized:" << endl;
	
	Formatpos pos;
	OBFormat* pFormat;
	const char* str=NULL;
	while(OBConversion::GetNextFormat(pos,str,pFormat))
	  cout << "  " << str << endl;
	cout << "\nSee further specific info and options using -H<format-type>, e.g. -Hcml" << endl;
}

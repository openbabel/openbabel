/**********************************************************************
Copyright (C) 2004 by Chris Morley

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
// Definition of OBConversion routines

#ifdef _WIN32
	#pragma warning (disable : 4786)
	#ifdef GUI
		#undef DATADIR
		#include "stdafx.h" //(includes<windows.h>
	#endif
#endif

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include "obconversion.h"

using namespace std;
namespace OpenBabel {

const char* OBFormat::TargetClassDescription()
{
	//Provides class of default format unless overridden
	if(OBConversion::GetDefaultFormat())
		return OBConversion::GetDefaultFormat()->TargetClassDescription();
	else
		return "";
};
const type_info& OBFormat::GetType()
{
	//Provides info on class of default format unless overridden
	if(OBConversion::GetDefaultFormat())
		return OBConversion::GetDefaultFormat()->GetType();
	else
		return typeid(this); //rubbish return if DefaultFormat not set
};

const char* OBFormat::SpecificationURL() {return "";};

//***************************************************

/** @class OBConversion
OBConversion maintains a list of the available formats, 
provides information on them, and controls the conversion process.

A conversion is carried out by the calling routine, usually in a
user interface or an application program, making an instance of
OBConversion. It is loaded with the in and out formats, any options
and (usually) the default streams for input and output. Then either
the Convert() function is called, which allows a single input file
to be converted, or the extended functionality of FullConvert()
is used. This allows multiple input and output files, allowing:
 - aggregation      - the contents of many input files converted 
                     and sent to one output file;
 - splitting        - the molecules from one input file sent to 
                     separate output files;
 - batch conversion - each input file converted to an output file.

These procedures constitute the "Convert" interface. OBConversion
and the user interface or application program do not need to be
aware of any other part of OpenBabel - mol.h is not #included. This
allows any chemical object derived from OBBase to be converted;
the type of object is decided by the input format class.
However,currently, almost all the conversions are for molecules of
class OBMol.
///
OBConversion can also be used with an "API" interface
called from programs which manipulate chemical objects. Input/output is
done with the Read() and Write() functions which work with any
chemical object, but need to have its type specified. (The 
ReadMolecule() and WriteMolecule() functions of the format classes
can also be used directly.)


Example code using OBConversion

<b>To read in a molecule, manipulate it and write it out.</b>

Set up an istream and an ostream, to and from files or elsewhere.
(cin and cout are used in the example). Specify the file formats.

@code
	OBConversion conv(&cin,&cout);
	if(conv.SetInAndOutFormats("SMI","MOL"))
	{	
		OBMol mol;
		if(conv.Read(&mol))
		...manipulate molecule 
		
		conv->Write(&mol);
	}
@endcode
	
A two stage construction is used to allow error handling
if the format ID is not recognized. This is necessary now that the
formats are dynamic and errors are not caught at compile time.
OBConversion::Read() is a templated function so that objects derived
from OBBase can also be handled, in addition to OBMol, if the format
routines are written appropriately.

<b>To do a file conversion without manipulating the molecule.</b>

@code
	#include "obconversion.h" //mol.h is not needed
	...set up an istream is and an ostream os 
	OBConversion conv(&is,&os);
	if(conv.SetInAndOutFormats("SMI","MOL"))
	{
		conv.SetOptions("h"); //Optional; (h adds expicit hydrogens)
		conv.Convert();
	}
@endcode

<b>To add automatic format conversion to an existing program.</b>

The existing program inputs from the file identified by the 
const char* filename into the istream is. The file is assumed to have
a format ORIG, but otherformats, identified by their file extensions,
can now be used.

@code
	ifstream ifs(filename); //Original code

	OBConversion conv;
	OBFormat* inFormat = conv.FormatFromExt(filename);
	OBFormat* outFormat = conv.GetFormat("ORIG");
	istream* pIn = &ifs; 
	stringstream newstream;
	if(inFormat && outFormat)
	{
		conv.SetInAndOutFormats(inFormat,outFormat);
		conv.Convert(pIn,&newstream);
		pIn=&newstream;
	}
	//else error; new features not available; fallback to original functionality 

	...Carry on with original code using pIn
@endcode

	In Windows a degree of independence from OpenBabel can be achieved using DLLs.
  This code would be linked with obconv.lib.
	At runtime the following DLLs would be in the executable directory:
	obconv.dll, obdll.dll, one or more *.obf format files.
*/
    
bool OBConversion::FormatFilesLoaded = false;

OBFormat* OBConversion::pDefaultFormat=NULL;

OBConversion::OBConversion(istream* is, ostream* os) : 
	StartNumber(1), EndNumber(0),pInFormat(NULL),pOutFormat(NULL),
	Index(0), pOb1(NULL), Count(-1), m_IsLast(true), MoreFilesToCome(false),
	OneObjectOnly(false)
{
	pInStream=is;
	pOutStream=os;
	strcpy(Dimension,"2D");
	LoadFormatFiles(); 
}

///This static function returns a reference to the FormatsMap
///which, because it is a static local variable is constructed only once.
///This fiddle is to avoid the "static initialization order fiasco"
///See Marshall Cline's C++ FAQ Lite document, www.parashift.com/c++-faq-lite/". 
FMapType& OBConversion::FormatsMap()
{
	static FMapType* fm = new FMapType;
	return *fm;
}

///////////////////////////////////////////////
OBConversion::OBConversion(const OBConversion& O)
{
	//Copy constructor. Needed because of strings
	Options=O.Options;
	GeneralOptions=O.GeneralOptions;
	Title=O.Title;
	pInFormat=O.pInFormat;
	pOutFormat=O.pOutFormat;
	pInStream=O.pInStream;
	pOutStream=O.pOutStream;
	StartNumber=O.StartNumber;
	EndNumber=O.EndNumber;
	strcpy(Dimension,O.Dimension);
	Index=0;
	FormatsMap();//rubbish
}

/////////////////////////////////////////////////
OBConversion::~OBConversion() 
{
}

//////////////////////////////////////////////////////
/// Class information on formats is collected by making an instance of the class
/// derived from OBFormat(only one is ever required).RegisterFormat() is called 
/// from its constructor. 
///
/// If the compiled format were to be stored separately, like in a DLL, the initialization
/// code would make an instance of the imported OBFormat class.
int OBConversion::RegisterFormat(const char* ID, OBFormat* pFormat)
{
	FormatsMap()[ID] = pFormat;
	if(pFormat->Flags() & DEFAULTFORMAT)
		pDefaultFormat=pFormat;
	return FormatsMap().size();
}

//////////////////////////////////////////////////////
int OBConversion::LoadFormatFiles()
{
	int count=0;
	if(FormatFilesLoaded) return 0;
	FormatFilesLoaded=true; //so will load files only once
#ifdef USING_DYNAMIC_LIBS
	//Depending on availablilty, look successively in 
	//FORMATFILE_DIR, executable directory,or current directory
	string TargetDir;
	#ifdef FORMATFILE_DIR
		TargetDir="FORMATFILE_DIR";
	#endif

	DLHandler::getConvDirectory(TargetDir);
	
	vector<string> files;
	if(!DLHandler::findFiles(files,DLHandler::getFormatFilePattern(),TargetDir)) return 0;

	vector<string>::iterator itr;
	for(itr=files.begin();itr!=files.end();itr++)
	{
		if(DLHandler::openLib(*itr))
			count++;
		else
			cerr << *itr << " did not load properly." << endl;
	}
#endif //USING_DYNAMIC_LIBS
	return count;
/*
#ifdef _WIN32
	//Windows specific
	//Load all the format DLLs (*.obf) in current directory
	//Returns the number loaded
	WIN32_FIND_DATA finddata;

	char path[MAX_PATH+1];
	if(GetModuleFileName(NULL,path,MAX_PATH)==0) return FALSE;
	char* p = strrchr(path,'\\');
	if(p && *p) *(p+1)='\0'; //chop off name after last '\'
	strcat(path,"*.obf");
	HANDLE hF = FindFirstFile(path,&finddata);
	if(hF==INVALID_HANDLE_VALUE) return 0;
	int count=0;
	do
	{
		if(LoadLibrary(finddata.cFileName))
			count++;
		else
			cerr << finddata.cFileName << " did not load properly. Error no. " << GetLastError() <<endl;
	}while(FindNextFile(hF,&finddata));
	FindClose(hF);
	return count;
#else
	//Code for other platforms goes here
	return 0;
#endif
*/
}


/**
  *Returns the ID + the first line of the description in str 
  *and a pointer to the format in pFormat.
  *If called with str==NULL the first format is returned;
  *subsequent formats are returned by calling with str!=NULL and the previous value of itr
  *returns false, and str and pFormat NULL, when there are no more formats.
  *Use like:
  *@code
  *	const char* str=NULL;
  *	Formatpos pos;
  *	while(OBConversion::GetNextFormat(pos,str,pFormat))
  *	{
  *		use str and pFormat
  *	}
  *@endcode
  */
/////////////////////////////////////////////////////////
bool OBConversion::GetNextFormat(Formatpos& itr, const char*& str,OBFormat*& pFormat)
{

	pFormat = NULL;
	if(str==NULL) 
		itr = FormatsMap().begin();
	else
		itr++;
	if(itr == FormatsMap().end())
	{
		str=NULL; pFormat=NULL;
		return false;
	}
	static string s; 
	s=itr->first;
	pFormat = itr->second;
	if(pFormat)
	{
		string description(pFormat->Description());
		s += "  :  ";
		s += description.substr(0,description.find('\n'));
	}

	if(pFormat->Flags() & NOTWRITABLE) s+=" [Readonly]";
	if(pFormat->Flags() & NOTREADABLE) s+=" [Writeonly]";

	str = s.c_str();
	return true;
}

//////////////////////////////////////////////////////
/// Sets the formats from their ids, e g CML.
/// If inID is NULL, the input format is left unchanged. Similarly for outID
/// Returns true if both formats have been successfully set at sometime
bool OBConversion::SetInAndOutFormats(const char* inID, const char* outID)
{
	if(inID)
		pInFormat = FindFormat(inID);
  if(outID)
		pOutFormat= FindFormat(outID);
	return pInFormat && pOutFormat
		&& !(pInFormat->Flags() & NOTREADABLE) && !(pOutFormat->Flags() & NOTWRITABLE);
}
//////////////////////////////////////////////////////

bool OBConversion::SetInAndOutFormats(OBFormat* pIn, OBFormat* pOut)
{
	pInFormat=pIn;
	pOutFormat=pOut;
	return !(pInFormat->Flags() & NOTREADABLE) && !(pOutFormat->Flags() & NOTWRITABLE);
}
//////////////////////////////////////////////////////

int OBConversion::Convert(istream* is, ostream* os) 
{
	if(is) pInStream=is;
	if(os) pOutStream=os;
	return Convert();
}

////////////////////////////////////////////////////
/// Actions the "convert" interface.
///	Calls the OBFormat class's ReadMolecule() which 
///	 - makes a new chemical object of its chosen type (e.g. OBMol)
///	 - reads an object from the input file
///	 - subjects the chemical object to 'transformations' as specified by the Options
///	 - calls AddChemObject to add it to a buffer. The previous object is first output 
///	   via the output Format's WriteMolecule(). During the output process calling
/// IsFirst() and GetIndex() (the number of objects including the current one already output.
/// allows more control, for instance writing \<cml\> and \</cml\> tags for multiple molecule outputs only. 
///
///	AddChemObject does not save the object passed to it if it is NULL (as a result of a DoTransformation())
///	or if the number of the object is outside the range defined by
///	StartNumber and EndNumber.This means the start and end counts apply to all chemical objects 
///	found whether or not they	are output.
///	
///	If ReadMolecule returns false the input conversion loop is exited. 
///
int OBConversion::Convert() 
{
	if(pInStream==NULL || pOutStream==NULL)
	{
		cerr << "input or output stream not set" << endl;
		return 0;
	}

	if(!pInFormat) return 0;
	SetStartAndEnd();
	if(OneObjectOnly)
	{
		EndNumber=1;
		OneObjectOnly=false;
	}

//	Index=0;//number objects output
	Count=0;//number objects processed
	ReadyToInput=true;
	m_IsLast=false;
	pOb1=NULL;

	//Input loop
	while(ReadyToInput && pInStream->peek() != EOF && pInStream->good())
	{
		if(pInStream==&cin)
			if(pInStream->peek()=='\n')break;

		if( !pInFormat->ReadChemObject(this)) break; 
		// Objects supplied to AddChemObject() which may output them after a delay
		//ReadyToInput may be made false in AddChemObject()
		// by WriteMolecule() returning false  or by Count==EndNumber		
	}
	
	//Output last object
	if(!MoreFilesToCome)
		m_IsLast=true;
	if(pOutFormat)
		if(!pOutFormat->WriteChemObject(this))
			Index--;
	
	//Put AddChemObject() into non-queue mode
	Count= -1; 
	EndNumber=StartNumber=0; pOb1=NULL;//leave tidy
	MoreFilesToCome=false;

	return Index; //The number actually output
}
//////////////////////////////////////////////////////
void OBConversion::SetStartAndEnd()
{
//Sets starting and ending molecule numbers by parsing GeneralOptions
	const char* p = GetGeneralOptions();
	while(*p)
	{
		switch(*(p++))
		{
		case '\"' : //ignore quoted characters
			p=strchr(p++,'\"');
			if(p++==NULL)
			{
				cerr << "Missing \" in options" <<endl;
			}
			break;
		case 'f':
			StartNumber=atoi(++p);
			p=strchr(p,'\"')+1; //get past ""
			break;
		case 'l':
			EndNumber=atoi(++p);
			p=strchr(p,'\"')+1;
			break;
		}
	}
}

//////////////////////////////////////////////////////
/// Retrieves an object stored by AddChemObject() during output
OBBase* OBConversion::GetChemObject()
{
	Index++;
	return pOb1;
}

//////////////////////////////////////////////////////
///	Called by ReadMolecule() to deliver an object it has read from an input stream.
/// Used in two modes: 
///  - When Count is negative it is left negative and the routine is just a store
///    for an OBBase object.  The negative value returned tells the calling
///    routine that no more objects are required.
///  - When count is >=0, probably set by Convert(), it acts as a queue of 2:
///    writing the currently stored value before accepting the supplied one. This delay
///    allows output routines to respond differently when the written object is the last.
///    Count is incremented with each call, even if pOb=NULL. 
///    Objects are not added to the queue if the count is outside the range  
///    StartNumber to EndNumber. There is no upper limit if EndNumber is zero. 
///    The return value is the number of objects, including this one, which have been
///    input (but not necessarily output).
int OBConversion::AddChemObject(OBBase* pOb)
{
	if(Count<0) 
	{
		pOb1=pOb;
		return Count;
	}
	Count++;
	if(pOb && Count>=StartNumber)//keeps reading objects but does nothing with them
	{	
		if(Count==EndNumber)
			ReadyToInput=false; //stops any more objects being read

		if(pOb1 && pOutFormat) //see if there is an object ready to be output
		{
			//Output object
			if (!pOutFormat->WriteChemObject(this))  
			{
				//faultly write, so finish
				--Index;
				ReadyToInput=false;
				return Count;
			}
		}
		pOb1=pOb;
	}
	return Count;
}
////////////////////////////////////////////////
//////////////////////////////////////////////////////
int OBConversion::GetOutputIndex() const
{
	//The number of objects actually written already from this instance of OBConversion
	return Index;
}
void OBConversion::SetOutputIndex(int indx)
{
	Index=indx;
}
//////////////////////////////////////////////////////
OBFormat* OBConversion::FindFormat(const char* ID)
{
	//Case insensitive
	if(FormatsMap().find(ID) == FormatsMap().end())
		return NULL;
  else
		return FormatsMap()[ID];
}

//////////////////////////////////////////////////
const char* OBConversion::GetOptions() const
{
	return(Options.c_str());
}

void OBConversion::SetOptions(const char* options)
{
	Options=options;
}

const char* OBConversion::GetGeneralOptions() const
{
	return(GeneralOptions.c_str());
}

void OBConversion::SetGeneralOptions(const char* options)
{
	GeneralOptions=options;
}

const char* OBConversion::GetTitle() const
{
	return(Title.c_str());
}
void OBConversion::SetTitle(const char* title)
{
	Title=title;
}
const char* OBConversion::GetDimension() const
{
	return Dimension;
}

void OBConversion::SetDimension(const char* dim)
{
	strcpy(Dimension,dim);
}

void OBConversion::SetMoreFilesToCome()
{
	MoreFilesToCome=true;
}

void OBConversion::SetOneObjectOnly()
{
	OneObjectOnly=true;
}	
/////////////////////////////////////////////////////////
OBFormat* OBConversion::FormatFromExt(const char* filename)
{
	char* p = strrchr(filename,'.');
	if(p)
		return FindFormat(p+1);
	return NULL; //if no extension		
}

//////////////////////////////////////////////////
/// Writes the object pOb but does not delete it afterwards.
/// The output stream is lastingly changed if pout is not NULL
/// Returns true if successful.
bool OBConversion::Write(OBBase* pOb, ostream* pos)
{
	if(pos)
		pOutStream=pos;
	if(!pOutFormat) return false;
	return pOutFormat->WriteMolecule(pOb,this);
}

////////////////////////////////////////////
const char* OBConversion::Description()
{
return "Conversion options\n \
 -f <#> Start import at molecule # specified\n \
 -l <#> End import at molecule # specified\n \
 -z     Use the options used last time\n\n";
}

////////////////////////////////////////////
bool OBConversion::SaveOptionsToFile(const char* filename)
{
	//Returns true if successful
	ofstream os(filename);
	if(os.good())
	{
		os << GetGeneralOptions() << endl;
		os << GetOptions() << ends;
		return true;
	}
	else
		return false;
}

////////////////////////////////////////////
bool OBConversion::RestoreOptionsFromFile(const char* filename)
{
	//Replaces all the conversion options with a string from the specified file
	ifstream is(filename);
	if(is.good())
	{
		const int BUFLEN = 256;
		char buf[BUFLEN];
		is.getline(buf,BUFLEN);
		SetGeneralOptions(buf);
		is.getline(buf,BUFLEN);
		SetOptions(buf);
		return true;
	}
	return false;
}

////////////////////////////////////////////
bool OBConversion::IsLast()
{
	return m_IsLast;
}

////////////////////////////////////////////
/// If ch is not in the option string outside quoted text, returns NULL
/// If it is the return value is not NULL and points to the
/// following quoted text if there is any (but it is not NULL terminated).
const char* OBConversion::IsOption(const char ch) const
{
	const char* p = Options.c_str(); 
	while(*p)
	{
		const char* str=p;
		const char* pp=p;
		if(*p++=='\"') //look for quote
		{
			str=p; //start of quote string
			p=strchr(p,'\"')+1; //...find matching " and move on
			if(p==NULL) 
			{
				cerr << "Unmatched quotes in option string" <<endl;
				return NULL;
			}
		}
		if(*pp==ch) return str;
	}
	return NULL;
}

/////////////////////////////////////////////////
string OBConversion::BatchFileName(string& BaseName, string& InFile)
{
	//Replaces * in BaseName by InFile without extension and path
	string ofname(BaseName);
	int pos = ofname.find('*');
	if(pos>=0)
	{
		//Replace * by input filename
		int posdot=(InFile).rfind('.');
		if(posdot==-1) posdot=(InFile).size();
		int posname=(InFile).find_last_of("\\/");
		ofname.replace(pos,1, (InFile), posname+1, posdot-posname-1);
	}
	return ofname;	
}

////////////////////////////////////////////////
string OBConversion::IncrementedFileName(string& BaseName, const int Count)
{
	//Replaces * in BaseName by Count
	string ofname(BaseName);
	int pos = ofname.find('*');
	if(pos>=0)
	{
		char num[33];
		snprintf(num, 33, "%d", Count);
		//		itoa(Count,num,10);
		ofname.replace(pos,1, num);
	}
	return ofname;		
}
////////////////////////////////////////////////////

/**
Makes input and output streams, and carries out normal,
batch, aggregation, and splitting conversion.

Normal
Done if FileList contains a single file name and OutputFileName
does not contain a *.

Aggregation
Done if FileList has more than one file name and OutputFileName does
not contain * . All the chemical objects are OutputFileName
converted and sent to the single output file.
 
Splitting
Done if FileList contains a single file name and OutputFileName
contains a * . Each chemical object in the input file converted
and sent to a separate file whose name is OutputFileName with the
* replaced by 1, 2, 3, etc.
For example, if OutputFileName is NEW*.smi then the output files are
NEW1.smi, NEW2.smi, etc.

Batch Conversion
Done if FileList has more than one file name and contains a * .
Each input file is converted to an output file whose name is
OutputFileName with the * replaced by the inputfile name without its
path and extension.
So if the input files were inpath/First.cml, inpath/Second.cml
and OutputFileName was NEW.mol, the output files would be
NEWFirst.mol, NEWSecond.mol.
   
If FileList is empty, the input stream already set (usually in the
constructor) is used. If OutputFileName is empty, the output stream
already set is used.

On exit, OutputFileList contains the names of the output files.

Returns the number of Chemical objects converted.
*/
int OBConversion::FullConvert(vector<string>& FileList, string& OutputFileName,
													 vector<string>& OutputFileList)
{
	
	istream* pInStream;
	ostream* pOutStream;
	ifstream is;
	ofstream os;
	bool HasMultipleOutputFiles=false;
	int Count=0;

	try
	{
		ofstream ofs;

		//OUTPUT
		if(OutputFileName.empty())
			pOutStream = NULL; //use existing stream
		else
		{
			if(OutputFileName.find_first_of('*')!=string::npos) HasMultipleOutputFiles = true;
			if(!HasMultipleOutputFiles)
			{
				os.open(OutputFileName.c_str());
				if(!os)
				{
					cerr << "Cannot write to " << OutputFileName <<endl;
					return 0;
				}
				OutputFileList.push_back(OutputFileName);
				pOutStream=&os;
			}
		}

		//INPUT
		if(FileList.empty())
			pInStream = NULL;
		else
		{
			if(FileList.size()>1)
			{
				//multiple input files
				vector<string>::iterator itr, tempitr;
				tempitr = FileList.end();
				tempitr--;
				for(itr=FileList.begin();itr!=FileList.end();itr++)
				{
					ifstream ifs((*itr).c_str());
					if(!ifs)
					{
						cerr << "Cannot open " << *itr <<endl;
						continue;
					}

					if(HasMultipleOutputFiles)
					{
						//Batch conversion
						string batchfile = BatchFileName(OutputFileName,*itr);
						if(ofs.is_open()) ofs.close();
						ofs.open(batchfile.c_str());
						if(!ofs) 
						{
							cerr << "Cannot open " << batchfile << endl;
							return Count;
						}
						OutputFileList.push_back(batchfile);
						SetOutputIndex(0); //reset for new file
						Count += Convert(&ifs,&ofs);					
					}
					else
					{
						//Aggregation
						if(itr!=tempitr) SetMoreFilesToCome();
						Count = Convert(&ifs,pOutStream);					
					}
				}
				return Count;
			}
			else
			{			
				//Single input file
				is.open(FileList[0].c_str());
				if(!is) 
				{
					cerr << "Cannot open " << FileList[0] <<endl;
					return 0;
				}
				pInStream=&is;

				if(HasMultipleOutputFiles)
				{
					//Splitting
					//Output is put in a temporary stream and written to a file
					//with an augmenting name only when it contains a valid object. 
					SetOneObjectOnly();
					int Indx=1;
					for(;;)
					{
						stringstream ss;
						SetOutputIndex(0); //reset for new file
						
						int ThisFileCount = Convert(pInStream,&ss);
						if(ThisFileCount==0) break;
						Count+=ThisFileCount;

						if(ofs.is_open()) ofs.close();
						string incrfile = IncrementedFileName(OutputFileName,Indx++);
						ofs.open(incrfile.c_str());
						if(!ofs)
						{
							cerr << "Cannot write to " << incrfile << endl;
							return Count;
						}
						OutputFileList.push_back(incrfile);
						ofs << ss.rdbuf();
						ofs.close();
						ss.clear();
					}
					return Count;
				}
			}
		}

		//Single input and output files
		Count = Convert(pInStream,pOutStream);
		return Count;
	}
	catch(...)
	{
		cerr << "Conversion failed with an exception" <<endl;
		return Count;
	}
}

}//namespace OpenBabel


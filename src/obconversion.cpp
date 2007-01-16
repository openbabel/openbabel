/**********************************************************************
obconversion.cpp -  Declaration of OBFormat and OBConversion

Copyright (C) 2004 by Chris Morley
Some portions Copyright (C) 2005-2006 by Geoffrey Hutchison

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
#include <openbabel/babelconfig.h>

#ifdef _WIN32
	#pragma warning (disable : 4786)

	//using 'this' in base class initializer
	#pragma warning (disable : 4355)

	#ifdef GUI
		#undef DATADIR
		#include "stdafx.h" //(includes<windows.h>
	#endif
#endif

#include <iosfwd>
#include <fstream>
#include <sstream>
#include <string>
#include <map>

#include <openbabel/obconversion.h>
#ifndef NO_NEWLINEBUF
  #include <openbabel/newlinebuf.h>
#endif

#ifdef HAVE_LIBZ
#include "zipstream.h"
#endif

#if !HAVE_STRNCASECMP
extern "C" int strncasecmp(const char *s1, const char *s2, size_t n);
#endif

#ifndef BUFF_SIZE
#define BUFF_SIZE 32768
#endif

using namespace std;
namespace OpenBabel {

  /** @class OBFormat obconversion.h <openbabel/obconversion.h>
      Two sets of Read and Write functions are specified for each format
      to handle two different requirements.
      The "Convert" interface is for use in file format conversion applications. The
      user interface, a console, a GUI, or another program is kept unaware of the
      details of the chemistry and does not need to \#include mol.h. It is then
      necessary to manipulate only pointers to OBBase in OBConversion and the user
      interface, with all the construction and deletion of OBMol etc objects being
      done in the Format classes or the OB core. The convention  with "Covert"
      interface functions is that chemical objects are made on the heap with new
      in the ReadChemicalObject() functions and and deleted in WriteChemicalObject()
      functions

      The "API" interface is for programatic use of the OB routines in application
      programs where mol.h is \#included. There is generally no creation or
      destruction of objects in ReadMolecule() and WriteMolecule() and no restriction
      on whether the pointers are to the heap or the stack.
  **/
  const char* OBFormat::TargetClassDescription()
  {
    //Provides class of default format unless overridden
    if(OBConversion::GetDefaultFormat())
      return OBConversion::GetDefaultFormat()->TargetClassDescription();
    else
      return "";
  }
  const type_info& OBFormat::GetType()
  {
    //Provides info on class of default format unless overridden
    if(OBConversion::GetDefaultFormat())
      return OBConversion::GetDefaultFormat()->GetType();
    else
      return typeid(this); //rubbish return if DefaultFormat not set
  }

  //***************************************************

  /** @class OBConversion obconversion.h <openbabel/obconversion.h>
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
      aware of any other part of OpenBabel - mol.h is not \#included. This
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
            // ...manipulate molecule 
		
         conv->Write(&mol);
      }
      @endcode
	
      A two stage construction is used to allow error handling
      if the format ID is not recognized. This is necessary now that the
      formats are dynamic and errors are not caught at compile time.
      OBConversion::Read() is a templated function so that objects derived
      from OBBase can also be handled, in addition to OBMol, if the format
      routines are written appropriately.

      <b>To make a molecule from a SMILES string.</b>
      @code
      std::string SmilesString;
      OBMol mol;
      stringstream ss(SmilesString)
      OBConversion conv(&ss);
      if(conv.SetInFormat("smi") && conv.Read(&mol))
         // ...
      @endcode

      <b>To do a file conversion without manipulating the molecule.</b>

      @code
      #include <openbabel/obconversion.h> //mol.h is not needed
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
    
  int OBConversion::FormatFilesLoaded = 0;

  OBFormat* OBConversion::pDefaultFormat=NULL;

  OBConversion::OBConversion(istream* is, ostream* os) : 
    pInFormat(NULL),pOutFormat(NULL), Index(0), StartNumber(1),
    EndNumber(0), Count(-1), m_IsFirstInput(true), m_IsLast(true),
    MoreFilesToCome(false), OneObjectOnly(false),NeedToFreeInStream(false),
    NeedToFreeOutStream(false), pOb1(NULL), pAuxConv(NULL)
  {
    pInStream=is;
    pOutStream=os;
    if (FormatFilesLoaded == 0)
      FormatFilesLoaded = LoadFormatFiles();
	
    //These options take a parameter
    RegisterOptionParam("f", NULL, 1,GENOPTIONS);
    RegisterOptionParam("l", NULL, 1,GENOPTIONS);
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

  ///This static function returns a reference to the FormatsMIMEMap
  ///which, because it is a static local variable is constructed only once.
  ///This fiddle is to avoid the "static initialization order fiasco"
  ///See Marshall Cline's C++ FAQ Lite document, www.parashift.com/c++-faq-lite/". 
  FMapType& OBConversion::FormatsMIMEMap()
  {
    static FMapType* fm = new FMapType;
    return *fm;
  }

  /////////////////////////////////////////////////
  OBConversion::OBConversion(const OBConversion& o)
  {
    Index          = o.Index;
    Count          = o.Count;
    StartNumber    = o.StartNumber;
    EndNumber      = o.EndNumber;
    pInFormat      = o.pInFormat;
    pInStream      = o.pInStream;
    pOutFormat     = o.pOutFormat;
    pOutStream     = o.pOutStream;
    OptionsArray[0]= o.OptionsArray[0];
    OptionsArray[1]= o.OptionsArray[1];
    OptionsArray[2]= o.OptionsArray[2];
    InFilename     = o.InFilename;
    rInpos         = o.rInpos;
    wInpos         = o.wInpos;
    rInlen         = o.rInlen;
    wInlen         = o.wInlen;
    m_IsLast       = o.m_IsLast;
    MoreFilesToCome= o.MoreFilesToCome;
    OneObjectOnly  = o.OneObjectOnly;
    pOb1           = o.pOb1;
    ReadyToInput   = o.ReadyToInput;
    m_IsFirstInput = o.m_IsFirstInput;
    NeedToFreeInStream = o.NeedToFreeInStream;
    NeedToFreeOutStream = o.NeedToFreeOutStream;

    pAuxConv       = NULL;
  }
  ////////////////////////////////////////////////

  OBConversion::~OBConversion() 
  {
    if(pAuxConv!=this)
      delete pAuxConv;
    // Free any remaining streams from convenience functions
    if(pInStream && NeedToFreeInStream) {
      delete pInStream;
      NeedToFreeInStream = false;
    }
    if(pOutStream && NeedToFreeOutStream) {
      delete pOutStream;
      NeedToFreeOutStream = false;
    }
  }
  //////////////////////////////////////////////////////

  /// Class information on formats is collected by making an instance of the class
  /// derived from OBFormat(only one is usually required). RegisterFormat() is called 
  /// from its constructor. 
  ///
  /// If the compiled format is stored separately, like in a DLL or shared library,
  /// the initialization code makes an instance of the imported OBFormat class.
  int OBConversion::RegisterFormat(const char* ID, OBFormat* pFormat, const char* MIME)
  {
    FormatsMap()[ID] = pFormat;
    if (MIME)
      FormatsMIMEMap()[MIME] = pFormat;
    if(pFormat->Flags() & DEFAULTFORMAT)
      pDefaultFormat=pFormat;
    return FormatsMap().size();
  }

  //////////////////////////////////////////////////////
  int OBConversion::LoadFormatFiles()
  {
    int count=0;
    //	if(FormatFilesLoaded) return 0;
    //	FormatFilesLoaded=true; //so will load files only once
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
    for(itr=files.begin();itr!=files.end();++itr)
      {
        if(DLHandler::openLib(*itr))
          count++;
        else
          obErrorLog.ThrowError(__FUNCTION__, *itr + " did not load properly", obError);
      }
#else
    count = 1; //avoid calling this function several times
#endif //USING_DYNAMIC_LIBS
    return count;
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
   *	OBConversion conv; // dummy to make sure static data is available
   *	while(OBConversion::GetNextFormat(pos,str,pFormat))
   *	{
   *		// use str and pFormat
   *	}
   *@endcode
   *
   * NOTE: Because of dynamic loading problems, it is usually necessary to 
   *  declare a "dummy" OBConversion object to access this static method.
   *  (Not elegant, but will hopefully be fixed in the future.)
   */
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
    s =itr->first;
    pFormat = itr->second;
    if(pFormat)
      {
        string description(pFormat->Description());
        s += " -- ";
        s += description.substr(0,description.find('\n'));
      }

    if(pFormat->Flags() & NOTWRITABLE) s+=" [Read-only]";
    if(pFormat->Flags() & NOTREADABLE) s+=" [Write-only]";

    str = s.c_str();
    return true;
  }

  //////////////////////////////////////////////////////
  /// Sets the formats from their ids, e g CML.
  /// If inID is NULL, the input format is left unchanged. Similarly for outID
  /// Returns true if both formats have been successfully set at sometime
  bool OBConversion::SetInAndOutFormats(const char* inID, const char* outID)
  {
    return SetInFormat(inID) && SetOutFormat(outID);
  }
  //////////////////////////////////////////////////////

  bool OBConversion::SetInAndOutFormats(OBFormat* pIn, OBFormat* pOut)
  {
    return SetInFormat(pIn) && SetOutFormat(pOut);
  }
  //////////////////////////////////////////////////////
  bool OBConversion::SetInFormat(OBFormat* pIn)
  {
    if(pIn==NULL)
      return true;
    pInFormat=pIn;
    return !(pInFormat->Flags() & NOTREADABLE);
  }
  //////////////////////////////////////////////////////
  bool OBConversion::SetOutFormat(OBFormat* pOut)
  {
    pOutFormat=pOut;
    return !(pOutFormat->Flags() & NOTWRITABLE);
  }
  //////////////////////////////////////////////////////
  bool OBConversion::SetInFormat(const char* inID)
  {
    if(inID)
      pInFormat = FindFormat(inID);
    return pInFormat && !(pInFormat->Flags() & NOTREADABLE);
  }
  //////////////////////////////////////////////////////

  bool OBConversion::SetOutFormat(const char* outID)
  {
    if(outID)
      pOutFormat= FindFormat(outID);
    return pOutFormat && !(pOutFormat->Flags() & NOTWRITABLE);
  }

  //////////////////////////////////////////////////////
  int OBConversion::Convert(istream* is, ostream* os) 
  {
    if (is) pInStream=is;
    if (os) pOutStream=os;
    ostream* pOrigOutStream = pOutStream;

#ifdef HAVE_LIBZ
    zlib_stream::zip_istream zIn(*pInStream);
    if(zIn.is_gzip())
      pInStream = &zIn;

    zlib_stream::zip_ostream zOut(*pOutStream);
    if(IsOption("z",GENOPTIONS))
      {
        // make sure to output the header
        zOut.make_gzip();
        pOutStream = &zOut;
      }
#endif

#ifndef NO_NEWLINEBUF
    newlinebuf *filter;
    if (!(pInFormat->Flags() & READBINARY) && !(pInFormat->Flags() & READXML))
      {
        filter = new newlinebuf(pInStream->rdbuf());
        pInStream->rdbuf(filter);
      }
#endif
    int count = Convert();
    pOutStream = pOrigOutStream;
    return count;
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
        obErrorLog.ThrowError(__FUNCTION__, "input or output stream not set", obError);
        return 0;
      }

    if(!pInFormat) return 0;
    Count=0;//number objects processed

    if(!SetStartAndEnd())
      return 0;

    ReadyToInput=true;
    m_IsLast=false;
    pOb1=NULL;
    wInlen=0;

    //Input loop
    while(ReadyToInput && pInStream->peek() != EOF && pInStream->good())
      {
        if(pInStream==&cin)
          {
            if(pInStream->peek()=='\n')
              break;
          }
        else
          rInpos = pInStream->tellg();
		
        bool ret=false;
        try
          {
            ret = pInFormat->ReadChemObject(this);
            m_IsFirstInput=false;
          }		
        catch(...)
          {
            if(!IsOption("e", GENOPTIONS) && !OneObjectOnly)
              throw;
          }

        if(!ret)
          {
            //error or termination request: terminate unless
            // -e option requested and sucessfully can skip past current object
            if(!IsOption("e", GENOPTIONS) || pInFormat->SkipObjects(0,this)!=1) 
              break;
          }
        if(OneObjectOnly)
          break;
        // Objects supplied to AddChemObject() which may output them after a delay
        //ReadyToInput may be made false in AddChemObject()
        // by WriteMolecule() returning false  or by Count==EndNumber		
      }
	
    //Output last object
    m_IsLast= !MoreFilesToCome;

    //Output is deferred until the end with -C option
    if(pOutFormat && (!IsOption("C",GENOPTIONS) || m_IsLast))
      if(pOb1 && !pOutFormat->WriteChemObject(this))
        Index--;
	
    //Put AddChemObject() into non-queue mode
    Count= -1; 
    EndNumber=StartNumber=0; pOb1=NULL;//leave tidy
    MoreFilesToCome=false;
    OneObjectOnly=false;

    return Index; //The number actually output
  }
  //////////////////////////////////////////////////////
  bool OBConversion::SetStartAndEnd()
  {
    unsigned int TempStartNumber=0;
    const char* p = IsOption("f",GENOPTIONS);
    if(p)
      {
        StartNumber=atoi(p);
        if(StartNumber>1)
          {
            TempStartNumber=StartNumber;
            //Try to skip objects now
            int ret = pInFormat->SkipObjects(StartNumber-1,this);
            if(ret==-1) //error
              return false; 
            if(ret==1) //success:objects skipped
              {
                Count = StartNumber-1;
                StartNumber=0;
              }
          }
      }

    p = IsOption("l",GENOPTIONS);
    if(p)
      {
        EndNumber=atoi(p);
        if(TempStartNumber && EndNumber<TempStartNumber)
          EndNumber=TempStartNumber;
      }

    return true;
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
  ///    The return value is the return code from WriteChemObject - false on error.
  bool OBConversion::AddChemObject(OBBase* pOb)
  {
    if(Count<0) 
      {
        pOb1=pOb;
        return Count!=0;
      }
    Count++;
    if(Count>=(int)StartNumber)//keeps reading objects but does nothing with them
      {	
        if(Count==(int)EndNumber)
          ReadyToInput=false; //stops any more objects being read

        rInlen = pInStream->tellg() - rInpos;

        if(pOb)
          {
            if(pOb1 && pOutFormat) //see if there is an object ready to be output
              {
                //Output object
                if (!pOutFormat->WriteChemObject(this))  
                  {
                    //faultly write, so finish
                    --Index;
                    //ReadyToInput=false;
                    pOb1=NULL;
                    return false;
                  }
                //Stop after writing with single object output files
                if(pOutFormat->Flags() & WRITEONEONLY)
                  {
                    // if there are more molecules to output, send a warning
                    stringstream errorMsg;
                    errorMsg << "WARNING: You are attempting to convert a file"
                             << " with multiple molecule entries into a format"
                             << " which can only store one molecule. The current"
                             << " output will only contain the first molecule.\n\n";

                    errorMsg << "To convert this input into multiple separate"
                             << " output files, with one molecule per file, try:\n"
                             << "babel [input] [ouptut] -m\n\n";

                    errorMsg << "To pick one particular molecule"
                             << " (e.g., molecule 4), try:\n"
                             << "babel -f 4 -l 4 [input] [output]" << endl;

                    obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);

                    ReadyToInput = false;
                    pOb1 = NULL;
                    return true;
                  }
              }
            pOb1=pOb;
            wInpos = rInpos; //Save the position in the input file to be accessed when writing it
            wInlen = rInlen;
          }
      }
    return true;
  }
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
  const char* OBConversion::GetTitle() const
  {
    return(InFilename.c_str());
  }

  void OBConversion::SetMoreFilesToCome()
  {
    MoreFilesToCome=true;
  }

  void OBConversion::SetOneObjectOnly(bool b)
  {
    OneObjectOnly=b;
    m_IsLast=b;
  }	

  /////////////////////////////////////////////////////////
  OBFormat* OBConversion::FormatFromExt(const char* filename)
  {
    string file = filename;
    string::size_type extPos = file.rfind('.');

    if(extPos!=string::npos)
      {
        // only do this if we actually can read .gz files
#ifdef HAVE_LIBZ
        if (file.substr(extPos,3) == ".gz")
          {
            file.erase(extPos);
            extPos = file.rfind('.');
            if (extPos!=string::npos)
              return FindFormat( (file.substr(extPos + 1, file.size())).c_str() );
          }
        else
#endif
          return FindFormat( (file.substr(extPos + 1, file.size())).c_str() );
      }
    return NULL; //if no extension		
  }

  OBFormat* OBConversion::FormatFromMIME(const char* MIME)
  {
    if(FormatsMIMEMap().find(MIME) == FormatsMIMEMap().end())
      return NULL;
    else
      return FormatsMIMEMap()[MIME];
  }

  bool	OBConversion::Read(OBBase* pOb, std::istream* pin, bool clearFirst)
  {
    if(pin) pInStream=pin;

    if(!pInFormat || !pInStream) return false;

#ifdef HAVE_LIBZ
    zlib_stream::zip_istream zIn(*pInStream);
    if(zIn.is_gzip())
      pInStream = &zIn;
#endif

#ifndef NO_NEWLINEBUF
    newlinebuf *filter;
    if (!(pInFormat->Flags() & READBINARY) && !(pInFormat->Flags() & READXML))
      {
        filter = new newlinebuf(pInStream->rdbuf());
        pInStream->rdbuf(filter);
      }
#endif

    // we really need to do this, but requires including base.h -GH
    //    if (clearFirst)
    //      pOb->Clear();
    return pInFormat->ReadMolecule(pOb, this);
  }
  //////////////////////////////////////////////////
  /// Writes the object pOb but does not delete it afterwards.
  /// The output stream is lastingly changed if pos is not NULL
  /// Returns true if successful.
  bool OBConversion::Write(OBBase* pOb, ostream* pos)
  {
    if(pos) pOutStream=pos;

    if(!pOutFormat || !pOutStream) return false;

    ostream* pOrigOutStream = pOutStream;
#ifdef HAVE_LIBZ
#ifndef _WIN32
    zlib_stream::zip_ostream zOut(*pOutStream);
    if(IsOption("z",GENOPTIONS))
      {
        // make sure to output the header
        zOut.make_gzip();
        pOutStream = &zOut;
      }
#endif
#endif

    bool ret = pOutFormat->WriteMolecule(pOb,this);
    pOutStream = pOrigOutStream;
    return ret;
  }

  //////////////////////////////////////////////////
  /// Writes the object pOb but does not delete it afterwards.
  /// The output stream not changed (since we cannot write to this string later)
  /// Returns true if successful.
  std::string OBConversion::WriteString(OBBase* pOb, bool trimWhitespace)
  {
    ostream *oldStream = pOutStream; // save old output
    stringstream newStream;
    string temp;

    if(pOutFormat)
      {
        Write(pOb, &newStream);
      }
    pOutStream = oldStream;

    temp = newStream.str();
    if (trimWhitespace) // trim the trailing whitespace
      {
        string::size_type notwhite = temp.find_last_not_of(" \t\n\r");
        temp.erase(notwhite+1);
      }
    return temp;
  }

  //////////////////////////////////////////////////
  /// Writes the object pOb but does not delete it afterwards.
  /// The output stream is lastingly changed to point to the file
  /// Returns true if successful.
  bool OBConversion::WriteFile(OBBase* pOb, string filePath)
  {
    if(!pOutFormat) return false;

    // if we have an old stream, free this first before creating a new one
    if (pOutStream && NeedToFreeOutStream) {
      delete pOutStream;
    }

    ofstream *ofs = new ofstream;
    NeedToFreeOutStream = true; // make sure we clean this up later
    ios_base::openmode omode = 
      pOutFormat->Flags() & WRITEBINARY ? ios_base::out|ios_base::binary : ios_base::out;

    ofs->open(filePath.c_str(),omode);
    if(!ofs || !ofs->good())
      {
        obErrorLog.ThrowError(__FUNCTION__,"Cannot write to " + filePath, obError);
        return false;
      }

    return Write(pOb, ofs);
  }

  void OBConversion::CloseOutFile()
  {
    if (pOutStream && NeedToFreeOutStream)
    {
      delete pOutStream;
      NeedToFreeOutStream = false;
      pOutStream = NULL;
    }
  }

  ////////////////////////////////////////////
  bool	OBConversion::ReadString(OBBase* pOb, std::string input, bool clearFirst)
  {
    // if we have an old stream, free this first before creating a new one
    if (pInStream && NeedToFreeInStream) {
      delete pInStream;
    }

    stringstream *pin = new stringstream(input);
    NeedToFreeInStream = true; // make sure we clean this up later
    return Read(pOb, pin, clearFirst);
  }


  ////////////////////////////////////////////
  bool	OBConversion::ReadFile(OBBase* pOb, std::string filePath, bool clearFirst)
  {
    if(!pInFormat) return false;

    // if we have an old stream, free this first before creating a new one
    if (pInStream && NeedToFreeInStream) {
      delete pInStream;
    }

    ifstream *ifs = new ifstream;
    NeedToFreeInStream = true; // make sure we free this
    ios_base::openmode imode = 
      pInFormat->Flags() & READBINARY ? ios_base::in|ios_base::binary : ios_base::in;

    ifs->open(filePath.c_str(),imode);
    if(!ifs || !ifs->good())
      {
        obErrorLog.ThrowError(__FUNCTION__,"Cannot read from " + filePath, obError);
        return false;
      }

    return Read(pOb,ifs, clearFirst);
  }


  ////////////////////////////////////////////
  const char* OBConversion::Description()
  {
#ifdef HAVE_LIBZ
    return "Conversion options\n \
 -f <#> Start import at molecule # specified\n \
 -l <#> End import at molecule # specified\n \
 -e Continue with next object after error, if possible\n \
 -z Compress the output with gzip\n \
 -k Attempt to translate keywords\n";
#else
    return "Conversion options\n \
 -f <#> Start import at molecule # specified\n \
 -l <#> End import at molecule # specified\n \
 -e Continue with next object after error, if possible\n \
 -k Attempt to translate keywords\n";
#endif
    // -t All input files describe a single molecule
  }

  ////////////////////////////////////////////
  bool OBConversion::IsLast()
  {
    return m_IsLast;
  }
  ////////////////////////////////////////////
  bool OBConversion::IsFirstInput()
  {
    return m_IsFirstInput;
  }

  /////////////////////////////////////////////////
  string OBConversion::BatchFileName(string& BaseName, string& InFile)
  {
    //Replaces * in BaseName by InFile without extension and path
    string ofname(BaseName);
    string::size_type pos = ofname.find('*');
    if(pos != string::npos)
      {
        //Replace * by input filename
        string::size_type posdot= InFile.rfind('.');
        if(posdot == string::npos)
          posdot = InFile.size();
        else {
#ifdef HAVE_LIBZ
          if (InFile.substr(posdot,3) == ".gz")
            {
              InFile.erase(posdot);
              posdot = InFile.rfind('.');
              if (posdot == string::npos)
                posdot = InFile.size();
            }
#endif
        }

        string::size_type posname= InFile.find_last_of("\\/");
        ofname.replace(pos,1, InFile, posname+1, posdot-posname-1);
      }
    return ofname;	
  }

  ////////////////////////////////////////////////
  string OBConversion::IncrementedFileName(string& BaseName, const int Count)
  {
    //Replaces * in BaseName by Count
    string ofname(BaseName);
    string::size_type pos = ofname.find('*');
    if(pos!=string::npos)
      {
        char num[33];
        snprintf(num, 33, "%d", Count);
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
     not contain * . All the chemical objects are converted and sent
     to the single output file.
 
     Splitting
     Done if FileList contains a single file name and OutputFileName
     contains a * . Each chemical object in the input file is converted
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
     and OutputFileName was NEW*.mol, the output files would be
     NEWFirst.mol, NEWSecond.mol.

     If FileList is empty, the input stream that has already been set
     (usually in the constructor) is used. If OutputFileName is empty,
     the output stream already set is used.

     On exit, OutputFileList contains the names of the output files.

     Returns the number of Chemical objects converted.
  */
  int OBConversion::FullConvert(std::vector<std::string>& FileList, std::string& OutputFileName,
                                std::vector<std::string>& OutputFileList)
  {
	
    istream* pIs=NULL;
    ostream* pOs=NULL;
    ifstream is;
    ofstream os;
    stringstream ssOut;
    bool HasMultipleOutputFiles=false;
    int Count=0;
    m_IsFirstInput=true;
    bool CommonInFormat = pInFormat ? true:false; //whether set in calling routine
    ios_base::openmode omode = 
      pOutFormat->Flags() & WRITEBINARY ? ios_base::out|ios_base::binary : ios_base::out;
    try
      {
        ofstream ofs;

        //OUTPUT
        if(OutputFileName.empty())
          pOs = NULL; //use existing stream
        else
          {
            if(OutputFileName.find_first_of('*')!=string::npos) HasMultipleOutputFiles = true;
            if(!HasMultipleOutputFiles)
              {
                //If the output file is the same as any of the input
                //files, send the output to a temporary stringstream
                vector<string>::iterator itr;
                for(itr=FileList.begin();itr!=FileList.end();++itr)
                  {
                    if(*itr==OutputFileName)
                      {
                        pOs = &ssOut;
                        break;
                      }
                  }
                if(itr==FileList.end())
                  {
                    os.open(OutputFileName.c_str(),omode);
                    if(!os)
                      {
                        obErrorLog.ThrowError(__FUNCTION__,"Cannot write to " + OutputFileName, obError);
                        return 0;
                      }
                    pOs=&os;
                  }
                OutputFileList.push_back(OutputFileName);
              }
          }

        if(IsOption("t",GENOPTIONS))
          {
            //Concatenate input file option (multiple files, single molecule)
            if(HasMultipleOutputFiles)
              {
                obErrorLog.ThrowError(__FUNCTION__,
                                      "Cannot have multiple output files and also concatenate input files (-t option)",obError);
                return 0;
              }

            stringstream allinput;
            vector<string>::iterator itr;
            for(itr=FileList.begin();itr!=FileList.end();++itr)
              {
                ifstream ifs((*itr).c_str());
                if(!ifs)
                  {
                    obErrorLog.ThrowError(__FUNCTION__,"Cannot open " + *itr, obError);
                    continue;
                  }
                allinput << ifs.rdbuf(); //Copy all file contents
                ifs.close();
              }
            Count = Convert(&allinput,pOs);
            return Count;
          }

        //INPUT
        if(FileList.empty())
          pIs = NULL;
        else
          {
            if(FileList.size()>1)
              {
                //multiple input files
                vector<string>::iterator itr, tempitr;
                tempitr = FileList.end();
                tempitr--;
                for(itr=FileList.begin();itr!=FileList.end();++itr)
                  {
                    InFilename = *itr;
                    ifstream ifs;
                    if(!OpenAndSetFormat(CommonInFormat, &ifs))
                      continue;

                    if(HasMultipleOutputFiles)
                      {
                        //Batch conversion
                        string batchfile = BatchFileName(OutputFileName,*itr);
                        if(ofs.is_open()) ofs.close();
                        ofs.open(batchfile.c_str(), omode);
                        if(!ofs) 
                          {
                            obErrorLog.ThrowError(__FUNCTION__,"Cannot open " + batchfile, obError);
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
                        Count = Convert(&ifs,pOs);					
                      }
                  }

                if(!os.is_open() && !OutputFileName.empty() && !HasMultipleOutputFiles)
                  {
                    //Output was written to temporary string stream. Output it to the file
                    os.open(OutputFileName.c_str(),omode);
                    if(!os)
                      {
                        obErrorLog.ThrowError(__FUNCTION__,"Cannot write to " + OutputFileName, obError);
                        return Count;
                      }
                    os << ssOut.rdbuf();
                  }
                return Count;
              }
            else
              {			
                //Single input file
                InFilename = FileList[0];
                if(!OpenAndSetFormat(CommonInFormat, &is))
                  return 0;
                pIs=&is;

                if(HasMultipleOutputFiles)
                  {
                    //Splitting
                    //Output is put in a temporary stream and written to a file
                    //with an augmenting name only when it contains a valid object. 
                    int Indx=1;
                    SetInStream(&is);
#ifdef HAVE_LIBZ
                    zlib_stream::zip_istream zIn(is);
#endif
                    for(;;)
                      {
                        stringstream ss;
                        SetOutStream(&ss);
                        SetOutputIndex(0); //reset for new file
                        SetOneObjectOnly();

#ifdef HAVE_LIBZ
                        if(Indx==1 && zIn.is_gzip())
                          SetInStream(&zIn);
#endif

                        int ThisFileCount = Convert();
                        if(ThisFileCount==0) break;
                        Count+=ThisFileCount;

                        if(ofs.is_open()) ofs.close();
                        string incrfile = IncrementedFileName(OutputFileName,Indx++);
                        ofs.open(incrfile.c_str(), omode);
                        if(!ofs)
                          {
                            obErrorLog.ThrowError(__FUNCTION__,"Cannot write to " + incrfile, obError);
                            return Count;
                          }
						
                        OutputFileList.push_back(incrfile);
#ifdef HAVE_LIBZ
                        if(IsOption("z",GENOPTIONS))
                          {
                            zlib_stream::zip_ostream zOut(ofs);
                            // make sure to output the header
                            zOut.make_gzip();
                            zOut << ss.rdbuf();
                          }
                        else
#endif
                          ofs << ss.rdbuf();

                        ofs.close();
                        ss.clear();
                      }
                    return Count;
                  }
              }
          }

        //Single input and output files
        Count = Convert(pIs,pOs);

        if(!os.is_open() && !OutputFileName.empty())
          {
            //Output was written to temporary string stream. Output it to the file
            os.open(OutputFileName.c_str(),omode);
            if(!os)
              {
                obErrorLog.ThrowError(__FUNCTION__,"Cannot write to " + OutputFileName, obError);
                return Count;
              }
            os << ssOut.rdbuf();
          }
        return Count;
      }
    catch(...)
      {
        obErrorLog.ThrowError(__FUNCTION__, "Conversion failed with an exception.",obError);
        return Count;
      }
  }

  bool OBConversion::OpenAndSetFormat(bool SetFormat, ifstream* is)
  {
    //Opens file using InFilename and sets pInFormat if requested
    if(!SetFormat)
      {
        pInFormat = FormatFromExt(InFilename.c_str());
        if(pInFormat==NULL)
          {
            string::size_type pos = InFilename.rfind('.');
            string ext;
            if(pos!=string::npos)
              ext = InFilename.substr(pos);
            obErrorLog.ThrowError(__FUNCTION__, "Cannot read input format \"" 
                                  + ext + '\"' + " for file \"" + InFilename + "\"",obError);
            return false;
          }
      }

    ios_base::openmode imode;
#ifdef ALL_READS_BINARY //Makes unix files compatible with VC++6
    imode = ios_base::in|ios_base::binary;
#else
    imode = pInFormat->Flags() & READBINARY ? ios_base::in|ios_base::binary : ios_base::in;
#endif

    is->open(InFilename.c_str(), imode);
    if(!is->good())
      {
        obErrorLog.ThrowError(__FUNCTION__, "Cannot open " + InFilename, obError);
        return false;
      }

    return true;
  }

  ///////////////////////////////////////////////
  void OBConversion::AddOption(const char* opt, Option_type opttyp, const char* txt)
  {
    //Also updates an option
    if(txt==NULL)
      OptionsArray[opttyp][opt]=string();
    else
      OptionsArray[opttyp][opt]=txt;
  }

  const char* OBConversion::IsOption(const char* opt, Option_type opttyp)
  {
    //Returns NULL if option not found or a pointer to the text if it is
    map<string,string>::iterator pos;
    pos = OptionsArray[opttyp].find(opt);
    if(pos==OptionsArray[opttyp].end())
      return NULL;
    return pos->second.c_str();
  }

  bool OBConversion::RemoveOption(const char* opt, Option_type opttyp)
  {
    return OptionsArray[opttyp].erase(opt)!=0;//true if was there
  }

  void OBConversion::SetOptions(const char* options, Option_type opttyp)
  {
    while(*options)
      {
        string ch(1, *options++);
        if(*options=='\"')
          {
            string txt = options+1;
            string::size_type pos = txt.find('\"');
            if(pos==string::npos)
              return; //options is illformed
            txt.erase(pos);
            OptionsArray[opttyp][ch]= txt;
            options += pos+2;
          }
        else
          OptionsArray[opttyp][ch] = string();
      }
  }

  OBConversion::OPAMapType& OBConversion::OptionParamArray(Option_type typ)
  {
    static OPAMapType* opa = new OPAMapType[3];
    return opa[typ];
  }

  void OBConversion::RegisterOptionParam(string name, OBFormat* pFormat,
                                         int numberParams, Option_type typ)
  {
    //Gives error message if the number of parameters conflicts with an existing registration
    map<string,int>::iterator pos;
    pos =	OptionParamArray(typ).find(name);
    if(pos!=OptionParamArray(typ).end())
      {
        if(pos->second!=numberParams)
          {
            string description("API");
            if(pFormat)
              description=pFormat->Description();
            obErrorLog.ThrowError(__FUNCTION__, 
                                  "The number of parameters needed by option \"" + name + "\" in " 
                                  + description.substr(0,description.find('\n'))
                                  + " differs from an earlier registration.", obError);
            return;
          }
      }
    OptionParamArray(typ)[name] = numberParams;
  }

  int OBConversion::GetOptionParams(string name, Option_type typ)
  {
    //returns the number of parameters registered for the option, or 0 if not found
    map<string,int>::iterator pos;
    pos =	OptionParamArray(typ).find(name);
    if(pos==OptionParamArray(typ).end())
      return 0;
    return pos->second;
  }

}//namespace OpenBabel

//! \file obconversion.cpp
//! \brief Implementation of OBFormat and OBConversion classes.

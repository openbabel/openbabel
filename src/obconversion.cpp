/**********************************************************************
obconversion.cpp -  Declaration of OBFormat and OBConversion

Copyright (C) 2004 by Chris Morley
Some portions Copyright (C) 2005-2006 by Geoffrey Hutchison

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
#include <locale>
#include <limits>
#include <typeinfo>

#include <stdlib.h>

#include <openbabel/obconversion.h>
//#include <openbabel/mol.h>
#include <openbabel/locale.h>

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
//using namespace boost::iostreams;

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
      OBConversion::Read() uses a pointer to OBBase, so that, in addition
      to OBMol, other kinds of objects, such as reactions, can also be handled
      if the format routines are written appropriately.

      <b>To make a molecule from a SMILES string.</b>
      @code
      std::string SmilesString;
      OBMol mol;
      stringstream ss(SmilesString)
      OBConversion conv(&ss);
      if(conv.SetInFormat("smi") && conv.Read(&mol))
         // ...
      @endcode
      
      An alternative way is more convenient if using bindings from another language:
      @code
      std::string SmilesString;
      OBMol mol;
      OBConversion conv;
      if(conv.SetInFormat("smi") && conv.ReadString(&mol, SmilesString))
         // ...
      @endcode

      <b>To do a file conversion without manipulating the molecule.</b>

      @code
      #include <openbabel/obconversion.h> //mol.h is not needed
      ...set up an istream is and an ostream os
      OBConversion conv(&is,&os);
      if(conv.SetInAndOutFormats("SMI","MOL"))
      {
         conv.AddOption("h",OBConversion::GENOPTIONS); //Optional; (h adds expicit hydrogens)
         conv.Convert();
      }
      @endcode

      <b>To read a multi-molecule file if using bindings from another language</b>

      The first molecule should be read using ReadFile, and subsequent molecules using Read,
      as follows:
      @code
      #include <openbabel/obconversion.h> //mol.h is not needed
      OBConversion conv;
      OBMol mol;
      bool success = conv.SetInFormat("sdf");
      if(success)
      {
         bool notatend = conv.ReadFile(&mol, "myfile.sdf");
         // Do something with mol
	 while(notatend)
	 {
             notatend = conv.Read(&mol);
	     // Do something with mol
	 }
      }
      @endcode

      <b>To add automatic format conversion to an existing program.</b>

      The existing program inputs from the file identified by the
      const char* filename into the istream is. The file is assumed to have
      a format ORIG, but other formats, identified by their file extensions,
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
  */

  int OBConversion::FormatFilesLoaded = 0;

//  OBFormat* OBConversion::pDefaultFormat=NULL;

  OBConversion::OBConversion(istream* is, ostream* os) :
    pInFormat(NULL),pOutFormat(NULL), Index(0), StartNumber(1),
    EndNumber(0), Count(-1), m_IsFirstInput(true), m_IsLast(true),
    MoreFilesToCome(false), OneObjectOnly(false), CheckedForGzip(false),
    NeedToFreeInStream(false), NeedToFreeOutStream(false),
    pOb1(NULL), pAuxConv(NULL),pLineEndBuf(NULL),wInpos(0),wInlen(0)
  {
    pInStream=is;
    pOutStream=os;
    if (FormatFilesLoaded == 0)
      FormatFilesLoaded = LoadFormatFiles();

    //These options take a parameter
    RegisterOptionParam("f", NULL, 1,GENOPTIONS);
    RegisterOptionParam("l", NULL, 1,GENOPTIONS);
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
    CheckedForGzip = o.CheckedForGzip;
    NeedToFreeInStream = o.NeedToFreeInStream;
    NeedToFreeOutStream = o.NeedToFreeOutStream;
    pLineEndBuf    = o.pLineEndBuf;
    pAuxConv       = NULL;
  }
  ///////////////////////////////////////////////

  OBConversion::~OBConversion()
  {
    if(pAuxConv!=this)
      if(pAuxConv)
      {
        delete pAuxConv;
        //pAuxConv has copies of pInStream, NeedToFreeInStream, pOutStream, NeedToFreeOutStream
        //and may have already deleted the streams. So do not do it again.
        NeedToFreeInStream = NeedToFreeOutStream = false;
      }
    // Free any remaining streams from convenience functions
    if(pInStream && NeedToFreeInStream) {
      delete pInStream;
      pInStream=NULL;
      NeedToFreeInStream = false;
    }
    if(pOutStream && NeedToFreeOutStream) {
      delete pOutStream;
      pOutStream=NULL;
      NeedToFreeOutStream = false;
    }
//    delete pLineEndBuf;
//   pLineEndBuf=NULL;
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
    return pFormat->RegisterFormat(ID, MIME);
  }

  //////////////////////////////////////////////////////
  int OBConversion::LoadFormatFiles()
  {
    int count=0;
    //	if(FormatFilesLoaded) return 0;
    //	FormatFilesLoaded=true; //so will load files only once
#if  defined(USING_DYNAMIC_LIBS)
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
        // Error handling is now handled by DLHandler itself
        //        else
        //          obErrorLog.ThrowError(__FUNCTION__, *itr + " did not load properly", obError);
      }
#else
    count = 1; //avoid calling this function several times
#endif //USING_DYNAMIC_LIBS

    //Make instances for plugin classes defined in the data file
    //This is hook for OBDefine, but does nothing if it is not loaded
    //or if plugindefines.txt is not found.
    OBPlugin* pdef = OBPlugin::GetPlugin("loaders","define");
    if(pdef)
    {
      static vector<string> vec(3);
      vec[1] = string("define");
      vec[2] = string("plugindefines.txt");
      pdef->MakeInstance(vec);
    }
    return count;
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
    return pOut && !(pOutFormat->Flags() & NOTWRITABLE);
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
    if (is) {
      pInStream=is;
      CheckedForGzip = false; // haven't checked this for gzip yet
    }
    if (os) pOutStream=os;
    ostream* pOrigOutStream = pOutStream;

#ifdef HAVE_LIBZ
    zlib_stream::zip_istream *zIn;

    // only try to decode the gzip stream once
    if (!CheckedForGzip) {
      zIn = new zlib_stream::zip_istream(*pInStream);
      if (zIn->is_gzip()) {
        pInStream = zIn;
        CheckedForGzip = true;
      }
      else
        delete zIn;
    }
#ifndef DISABLE_WRITE_COMPRESSION //Unsolved problem with compression under Windows
    zlib_stream::zip_ostream zOut(*pOutStream);
    if(IsOption("z",GENOPTIONS))
      {
        // make sure to output the header
        zOut.make_gzip();
        pOutStream = &zOut;
      }
#endif
#endif

    //The FilteringInputStreambuf delivers characters to the istream, pInStream,
    //and receives characters this stream's original rdbuf.
    //It filters them, converting CRLF and CR line endings to LF.
    //seek and tellg requests to the stream are passed through to the original
    //rdbuf. A FilteringInputStreambuf is installed only for appropriate formats
    //- not for binary or XML formats - if not already present.
    InstallStreamFilter();

    int count = Convert();

    pOutStream = pOrigOutStream;
#ifdef HAVE_LIBZ
    if ( CheckedForGzip ){ // Bug reported by Gert Thijs
		delete zIn;
		pInStream = is;
	}
#endif
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

    if(pInFormat->Flags() & READONEONLY)
      OneObjectOnly=true;

    //Input loop
    while(ReadyToInput && pInStream->good()) //Possible to omit? && pInStream->peek() != EOF
      {
        if(pInStream==&cin)
          {
            if(pInStream->peek()==-1) //Cntl Z Was \n but interfered with piping
              break;
          }
        else
          rInpos = pInStream->tellg();

        bool ret=false;
        try
          {
            ret = pInFormat->ReadChemObject(this);
/*            if (ret && IsOption("readconformers", GENOPTIONS)) {
              std::streampos pos = pInStream->tellg();
              OBMol nextMol;
              OBConversion conv;
              conv.SetOutFormat("smi");
              std::string ref_smiles = conv.WriteString(pOb1);
              while (pInStream->good()) {
                if (!pInFormat->ReadMolecule(&nextMol, this))
                  break;
                std::string smiles = conv.WriteString(&nextMol);
                if (smiles == ref_smiles) {
                  OBMol *pmol = dynamic_cast<OBMol*>(pOb1);
                  if (!pmol)
                    break;
                  unsigned int numCoords = nextMol.NumAtoms() * 3;
                  double *coords = nextMol.GetCoordinates();
                  double *conformer = new double [numCoords];
                  for (unsigned int i = 0; i < numCoords; ++i)
                    conformer[i] = coords[i];
                  pmol->AddConformer(conformer);
                  pos = pInStream->tellg();
                } else {
                  break;
                }
              }
              pInStream->seekg(pos, std::ios::beg);
            }
*/
            SetFirstInput(false);
          }
        catch(...)
          {
            if(!IsOption("e", GENOPTIONS) && !OneObjectOnly)
            {
              obErrorLog.ThrowError(__FUNCTION__, "Convert failed with an exception" , obError);
              return Index; // the number we've actually output so far
            }
          }

        if(!ret)
          {
            //error or termination request: terminate unless
            // -e option requested and successfully can skip past current object
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

    //Output is always occurs at the end with the --OutputAtEnd option
    bool oae = IsOption("OutputAtEnd",GENOPTIONS)!=NULL;
    if(pOutFormat && (!oae || m_IsLast))
      if((oae || pOb1) && !pOutFormat->WriteChemObject(this))
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
  ///    The return value is Count ((>0) or 0 if WriteChemObject returned false.
  int OBConversion::AddChemObject(OBBase* pOb)
  {
    if(Count<0)
      {
        pOb1=pOb;
        return Count; // <0
      }
    Count++;
    if(Count>=(int)StartNumber)//keeps reading objects but does nothing with them
      {
        if(Count==(int)EndNumber)
          ReadyToInput=false; //stops any more objects being read

        rInlen = pInStream ? pInStream->tellg() - rInpos : 0;
         // - (pLineEndBuf ? pLineEndBuf->getCorrection() : 0); //correction for CRLF

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
                    return 0;
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
                    return Count; // >0
                  }
              }
            pOb1=pOb;
            wInpos = rInpos; //Save the position in the input file to be accessed when writing it
            wInlen = rInlen;
          }
      }
    return Count; // >0
  }
  //////////////////////////////////////////////////////
  ///Returns the number of objects which have been output or are currently being output.
  ///The outputindex is incremented when an object for output is fetched by GetChemObject().
  ///So the function will return 1 if called from WriteMolecule() during output of the first object.
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
    return OBFormat::FindType(ID);
  }

  OBFormat* OBConversion::FindFormat(const std::string ID)
  {
    return OBFormat::FindType(ID.c_str());
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

    if(extPos!=string::npos // period found
       && (file.substr(extPos + 1, file.size())).find("/")==string::npos) // and period is after the last "/"
      {
        // only do this if we actually can read .gz files
#ifdef HAVE_LIBZ
        if (file.substr(extPos) == ".gz")
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

    // Check the filename if no extension (e.g. VASP does not use extensions):
    extPos = file.rfind('/');
    if(extPos!=string::npos) {
      return FindFormat( (file.substr(extPos + 1, file.size())).c_str() );
    }
    // If we are just passed the filename with no path, this should catch it:
    return FindFormat( file.c_str() ); //if no format found
  }

  OBFormat* OBConversion::FormatFromExt(const std::string filename)
  {
    return FormatFromExt(filename.c_str());
  }

  OBFormat* OBConversion::FormatFromMIME(const char* MIME)
  {
    return OBFormat::FormatFromMIME(MIME);
  }

  bool	OBConversion::Read(OBBase* pOb, std::istream* pin)
  {
    if(pin) {
      pInStream=pin;
      CheckedForGzip = false; // haven't set this stream to gzip (yet)
    }

    if(!pInFormat || !pInStream) return false;

#ifdef HAVE_LIBZ
    zlib_stream::zip_istream *zIn;

    // only try to decode the gzip stream once
    if (!CheckedForGzip) {
      zIn = new zlib_stream::zip_istream(*pInStream);
      if (zIn->is_gzip()) {
        pInStream = zIn;
        CheckedForGzip = true;
      }
      else
        delete zIn;
    }
#endif

    InstallStreamFilter();

    // Set the locale for number parsing to avoid locale issues: PR#1785463
    obLocale.SetLocale();

    // Also set the C++ stream locale
    locale originalLocale = pInStream->getloc(); // save the original
    locale cNumericLocale(originalLocale, "C", locale::numeric);
    pInStream->imbue(cNumericLocale);

    bool success = pInFormat->ReadMolecule(pOb, this);
    // return the C locale to the original one
    obLocale.RestoreLocale();
    // Restore the original C++ locale as well
    pInStream->imbue(originalLocale);

    // If we failed to read, plus the stream is over, then check if this is a stream from ReadFile
    if (!success && !pInStream->good() && NeedToFreeInStream) {
      ifstream *inFstream = dynamic_cast<ifstream*>(pInStream);
      if (inFstream != 0)
        inFstream->close(); // We will free the stream later, but close the file now
    }

    return success;
  }

    void OBConversion::InstallStreamFilter()
  {
    //Do not install filtering input stream if a binary or XML format
    //or if already installed in the current InStream (which may have changed).
    //Deleting any old LErdbuf before contructing a new one ensures there is
    //only one for each OBConversion object. It is deleted in the destructor.

    if(pInFormat && !(pInFormat->Flags() & (READBINARY | READXML)) && pInStream->rdbuf()!=pLineEndBuf)
    {
      delete pLineEndBuf;
      pLineEndBuf = NULL;
      pLineEndBuf = new LErdbuf(pInStream->rdbuf());
      pInStream->rdbuf(pLineEndBuf);
    }
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
#ifndef DISABLE_WRITE_COMPRESSION
    zlib_stream::zip_ostream zOut(*pOutStream);
    if(IsOption("z",GENOPTIONS))
      {
        // make sure to output the header
        zOut.make_gzip();
        pOutStream = &zOut;
      }
#endif
#endif
    SetOneObjectOnly(); //So that IsLast() returns true, which is important for XML formats

    // Set the locale for number parsing to avoid locale issues: PR#1785463
    obLocale.SetLocale();
    // Also set the C++ stream locale
    locale originalLocale = pOutStream->getloc(); // save the original
    locale cNumericLocale(originalLocale, "C", locale::numeric);
    pOutStream->imbue(cNumericLocale);

    // The actual work is done here
    bool success = pOutFormat->WriteMolecule(pOb,this);

    pOutStream = pOrigOutStream;
    // return the C locale to the original one
    obLocale.RestoreLocale();
    // Restore the C++ stream locale too
    pOutStream->imbue(originalLocale);

    return success;
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
  bool	OBConversion::ReadString(OBBase* pOb, std::string input)
  {
    // if we have an old stream, free this first before creating a new one
    if (pInStream && NeedToFreeInStream) {
      delete pInStream;
    }

    stringstream *pin = new stringstream(input);
    NeedToFreeInStream = true; // make sure we clean this up later
    return Read(pOb, pin);
  }


  ////////////////////////////////////////////
  bool	OBConversion::ReadFile(OBBase* pOb, std::string filePath)
  {
    if(!pInFormat) return false;

    // save the filename
    InFilename = filePath;

    // if we have an old stream, free this first before creating a new one
    if (pInStream && NeedToFreeInStream) {
      delete pInStream;
    }

    ifstream *ifs = new ifstream;
    NeedToFreeInStream = true; // make sure we free this

    ios_base::openmode imode = ios_base::in|ios_base::binary; //now always binary because may be gzipped
//      pInFormat->Flags() & READBINARY ? ios_base::in|ios_base::binary : ios_base::in;
      
    ifs->open(filePath.c_str(),imode);
    if(!ifs || !ifs->good())
      {
        obErrorLog.ThrowError(__FUNCTION__,"Cannot read from " + filePath, obError);
        return false;
      }

    return Read(pOb,ifs);
  }

  ////////////////////////////////////////////
  bool OBConversion::OpenInAndOutFiles(std::string infilepath, std::string outfilepath)
  {
    // if we have an old input stream, free this first before creating a new one
    if (pInStream && NeedToFreeInStream)
      delete pInStream;

    // if we have an old output stream, free this first before creating a new one
    if (pOutStream && NeedToFreeOutStream)
      delete pOutStream;

    ifstream *ifs = new ifstream;
    NeedToFreeInStream = true; // make sure we free this
    ifs->open(infilepath.c_str(),ios_base::in|ios_base::binary); //always open in binary mode
    if(!ifs || !ifs->good())
    {
      obErrorLog.ThrowError(__FUNCTION__,"Cannot read from " + infilepath, obError);
      return false;
    }
    pInStream = ifs;
    InFilename = infilepath;

    if(outfilepath.empty())//Don't open an outfile with an empty name.
      return true;
    ofstream *ofs = new ofstream;
    NeedToFreeOutStream = true; // make sure we clean this up later
    ofs->open(outfilepath.c_str(),ios_base::out|ios_base::binary);//always open in binary mode
    if(!ofs || !ofs->good())
    {
      obErrorLog.ThrowError(__FUNCTION__,"Cannot write to " + outfilepath, obError);
      return false;
    }
    pOutStream = ofs;

    return true;
  }

  ////////////////////////////////////////////
  const char* OBConversion::Description()
  {
    return
      "Conversion options\n"
      "-f <#> Start import at molecule # specified\n"
      "-l <#> End import at molecule # specified\n"
      "-e Continue with next object after error, if possible\n"
      #ifdef HAVE_LIBZ
      #ifndef DISABLE_WRITE_COMPRESSION //Unsolved problem with compression under Windows
      "-z Compress the output with gzip\n"
      #endif
      #endif
      "-k Attempt to translate keywords\n";
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
  void OBConversion::SetFirstInput(bool b)
  {
    m_IsFirstInput = b;
/*    //Also set or clear a general option
    if(b)
      AddOption("firstinput",GENOPTIONS);
    else
      RemoveOption("firstinput",GENOPTIONS);
  */
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
          if (InFile.substr(posdot) == ".gz")
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
  bool OBConversion::CheckForUnintendedBatch(const string& infile, const string& outfile)
  {
    //If infile == outfile issue error message and return false
    //If name without the extensions are the same issue warning and return true;
    //Otherwise return true
    bool ret=true;
    string inname1, inname2;
    string::size_type pos;
    pos = infile.rfind('.');
    if(pos != string::npos)
      inname1 = infile.substr(0,pos);
    pos = outfile.rfind('.');
    if(pos != string::npos)
      inname2 = infile.substr(0,pos);
    if(inname1==inname2)
      obErrorLog.ThrowError(__FUNCTION__,
"This was a batch operation. For splitting, use non-empty base name for the output files", obWarning);

    if(infile==outfile)
      return false;
    return true;
  }
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
     * replaced by 1, 2, 3, etc.  OutputFileName must have at least one
     character other than the * before the extension.
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
    //OBConversion::OutFilename = OutputFileName; ready for 2.4.0

    istream* pIs=NULL;
    ostream* pOs=NULL;
    ifstream is;
    ofstream os;
    stringstream ssOut, ssIn;
    bool HasMultipleOutputFiles=false;
    int Count=0;
    SetFirstInput();
    bool CommonInFormat = pInFormat ? true:false; //whether set in calling routine
    ios_base::openmode omode =
      pOutFormat->Flags() & WRITEBINARY ? ios_base::out|ios_base::binary : ios_base::out;
    obErrorLog.ClearLog();
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
          {
            pIs = NULL;
            if(HasMultipleOutputFiles)
              {
                obErrorLog.ThrowError(__FUNCTION__,"Cannot use multiple output files without an input file", obError);
                return 0;
              }
          }
        else
          {
            if(FileList.size()>1 || OutputFileName.substr(0,2)=="*.")
              {
                //multiple input files
                vector<string>::iterator itr, tempitr;
                tempitr = FileList.end();
                tempitr--;
                for(itr=FileList.begin();itr!=FileList.end();++itr)
                  {
                    InFilename = *itr;
                    ifstream ifs;
                    if(!OpenAndSetFormat(CommonInFormat, &ifs, &ssIn))
                      continue;
                    if(ifs)
                      pIs = &ifs;
                    else
                      pIs = &ssIn;

                    //pIs = ifs ? &ifs : &ssIn;


                    if(HasMultipleOutputFiles)
                      {
                        //Batch conversion
                        string batchfile = BatchFileName(OutputFileName,*itr);

                        //With inputs like babel test.xxx -oyyy -m
                        //the user may have wanted to do a splitting operation
                        //Issue a message and abort if xxx==yyy which would overwrite input file
                        if(FileList.size()==1 && !CheckForUnintendedBatch(batchfile, InFilename))
                          return Count;

                        if(ofs.is_open()) ofs.close();
                        ofs.open(batchfile.c_str(), omode);
                        if(!ofs)
                          {
                            obErrorLog.ThrowError(__FUNCTION__,"Cannot open " + batchfile, obError);
                            return Count;
                          }
                        OutputFileList.push_back(batchfile);
                        SetOutputIndex(0); //reset for new file
                        Count += Convert(pIs,&ofs);
                      }
                    else
                      {
                        //Aggregation
                        if(itr!=tempitr) SetMoreFilesToCome();
                        Count = Convert(pIs,pOs);
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
                if(!OpenAndSetFormat(CommonInFormat, &is, &ssIn))
                  return 0;
                if(is)
                  pIs =&is;
                else
                  pIs = &ssIn;

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
                        if(Indx==1 && zIn.is_gzip()) {
                          SetInStream(&zIn);
                          CheckedForGzip = true; // we know this one is gzip'ed
                        }
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
#ifndef DISABLE_WRITE_COMPRESSION
                        if(IsOption("z",GENOPTIONS))
                          {
                            zlib_stream::zip_ostream zOut(ofs);
                            // make sure to output the header
                            zOut.make_gzip();
                            zOut << ss.rdbuf();
                          }
                        else
#endif
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
    return Count;
  }

  bool OBConversion::OpenAndSetFormat(bool SetFormat, ifstream* is, stringstream* ss)
  {
    //Opens file using InFilename and sets pInFormat if requested
    if(ss && InFilename[0]=='-')
      {
        //InFilename is actually  -:SMILES
        is->setstate(ios::failbit); // do not use the input filestream...
        InFilename.erase(0, 2);
        if(SetFormat || SetInFormat("smi"))
          {
            ss->clear();
            ss->str(InFilename); //...use the stringstream instead
            return true;
          }
      }
      else if(!SetFormat)
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

#ifndef ALL_READS_BINARY
  #define ALL_READS_BINARY //now the default
#endif
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
#ifndef DOXYGEN_SHOULD_SKIP_THIS
/**<pre>
Built-in options for conversion of molecules
Additional options :
-d Delete hydrogens (make implicit)
-h Add hydrogens (make explicit)
-p <pH> Add hydrogens appropriate for this pH
-b Convert dative bonds e.g.[N+]([O-])=O to N(=O)=O
-r Remove all but the largest contiguous fragment
-c Center Coordinates
-C Combine mols in first file with others having same name
--filter <filterstring> Filter: convert only when tests are true:\n
--add <list> Add properties from descriptors\n
--delete <list> Delete properties in list\n
--append <list> Append properties or descriptors in list to title:\n
-s\smarts\ Convert only molecules matching SMARTS:\n
-v\smarts\ Convert only molecules NOT matching SMARTS: (not displayed in GUI)\n
--join Join all input molecules into a single output molecule
--separate Output disconnected fragments separately
--property <attrib> <value> add or replace a property (SDF)
--title <title> Add or replace molecule title
--addtotitle <text> Append to title
--writeconformers Output multiple conformers separately
--addindex Append output index to title
</pre>
**/
#endif
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
    if(!*options) // "" clears all
    {
      OptionsArray[opttyp].clear();
      return;
    }
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
    static OPAMapType opa[3];
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

  /**
   * Returns the list of supported input format
   */
  std::vector<std::string> OBConversion::GetSupportedInputFormat()
  {
    vector<string> vlist;
    OBPlugin::ListAsVector("formats", "in", vlist);
    return vlist;
  }
  /**
   * Returns the list of supported output format
   */
  std::vector<std::string> OBConversion::GetSupportedOutputFormat()
  {
    vector<string> vlist;
    OBPlugin::ListAsVector("formats", "out", vlist);
    return vlist;
  }

  void OBConversion::ReportNumberConverted(int count, OBFormat* pFormat)
  {
    //Send info message to clog. This constructed from the TargetClassDescription
    //of the specified class (or the output format if not specified).
    //Get the last word on the first line of the description which should
    //be "molecules", "reactions", etc and remove the s if only one object converted
    if(!pFormat)
      pFormat = pOutFormat;
    string objectname(pFormat->TargetClassDescription());
    string::size_type pos = objectname.find('\n');
    if(pos==std::string::npos)
      pos=objectname.size();
    if(count==1) --pos;
    objectname.erase(pos);
    pos = objectname.rfind(' ');
    if(pos==std::string::npos)
      pos=0;
    std::clog << count << objectname.substr(pos) << " converted" << endl;
  }

  void OBConversion::CopyOptions(OBConversion* pSourceConv, Option_type typ)
  {
    if(typ==ALL)
    for(int i=0;i<3;++i)
     OptionsArray[i]=pSourceConv->OptionsArray[i];
    else
     OptionsArray[typ]=pSourceConv->OptionsArray[typ];
  }

  int OBConversion::NumInputObjects()
  {
    istream& ifs = *GetInStream();
    ifs.clear(); //it may have been at eof
    //Save position of the input stream
    streampos pos = ifs.tellg();
    if(!ifs)
      return -1;

    //check that the input format supports SkipObjects()
    if(GetInFormat()->SkipObjects(0, this)==0)
    {
      obErrorLog.ThrowError(__FUNCTION__,
        "Input format does not have a SkipObjects function.", obError);
      return -1;
    }

    //counts objects only between the values of -f and -l options
    int nfirst=1, nlast=numeric_limits<int>::max();
    const char* p;
    if( (p=IsOption("f", GENOPTIONS)) ) // extra parens to indicate truth value
      nfirst=atoi(p);
    if( (p=IsOption("l", GENOPTIONS)) ) // extra parens to indicate truth value
      nlast=atoi(p);

    ifs.seekg(0); //rewind
    //Compressed files currently show an error here.***TAKE CHANCE: RESET ifs****
    ifs.clear();

    OBFormat* pFormat = GetInFormat();
    int count=0;
    //skip each object but stop after nlast objects
    while(ifs && pFormat->SkipObjects(1, this)>0  && count<nlast)
      ++count;

    ifs.clear(); //clear eof
    ifs.seekg(pos); //restore old position

    count -= nfirst-1;
    return count;
  }



  //The following function and typedef are deprecated, and are present only
  //for backward compatibility.
  //Use OBConversion::GetSupportedInputFormat(), OBConversion::GetSupportedOutputFormat(),
  //OBPlugin::List(), OBPlugin::OBPlugin::ListAsVector(),OBPlugin::OBPlugin::ListAsString(),
  //or (in extremis) OBPlugin::PluginIterator instead.

  typedef OBPlugin::PluginIterator Formatpos;

  bool OBConversion::GetNextFormat(Formatpos& itr, const char*& str,OBFormat*& pFormat)
  {

    pFormat = NULL;
    if(str==NULL)
      itr = OBPlugin::Begin("formats");
    else
      itr++;
    if(itr == OBPlugin::End("formats"))
      {
        str=NULL; pFormat=NULL;
        return false;
      }
    static string s;
    s =itr->first;
    pFormat = static_cast<OBFormat*>(itr->second);
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

  /**
   * @example obconversion_readstring.cpp
   * Reading a smiles string.
   */

  /**
   * @example obconversion_readstring.py
   * Reading a smiles string in python.
   */


}//namespace OpenBabel

//! \file obconversion.cpp
//! \brief Implementation of OBFormat and OBConversion classes.

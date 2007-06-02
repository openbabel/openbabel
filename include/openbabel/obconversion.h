/**********************************************************************
obconversion.h - Handle file conversions. Declaration of OBFormat, OBConversion

Copyright (C) 2004-2005 by Chris Morley

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

#ifndef OB_CONV_H
#define OB_CONV_H

#include <openbabel/babelconfig.h>

#include <iostream>
#include <fstream>
#include <sstream>

#include <string>
#include <vector>
#include <map>

#include <openbabel/dlhandler.h>
#include <openbabel/oberror.h>
#include <openbabel/lineend.h>

// These macros are used in DLL builds. If they have not
// been set in babelconfig.h, define them as nothing.
#ifndef OBCONV
	#define OBCONV
#endif
#ifndef OBDLL
	#define OBDLL
#endif

//using namespace std;
namespace OpenBabel {


  class OBBase;
  class OBConversion;
  OBERROR extern  OBMessageHandler obErrorLog;

  //*************************************************

  /// @brief Base class for file formats.
  // class introduction in obconversion.cpp
  class OBCONV OBFormat
    {
    public:
      /// @brief The "API" interface Read function.

      /// Reads a single object.
      /// Does not make a new object on the heap; 
      /// can be used with a pointer to an chem object on the heap or the stack.
      virtual bool ReadMolecule(OBBase* /*pOb*/, OBConversion* /*pConv*/)
        { std::cerr << "Not a valid input format"; return false;}

      /// @brief The "Convert" interface Read function.

      /// Possibly reads multiple new objects on the heap and subjects them 
      /// to its DoTransformations() function, which may delete them again. 
      /// Sends result to pConv->AddChemObject()
      virtual bool ReadChemObject(OBConversion* /*pConv*/)
        { std::cerr << "Not a valid input format"; return false;}

      /// @brief The "API" interface Write function.

      /// Writes a single object
      /// Does not delete the object; 
      /// can be used with a pointer to an chem object on the heap or the stack.
      /// \return false on error.
      virtual bool WriteMolecule(OBBase* /*pOb*/, OBConversion* /*pConv*/)
        { std::cerr << "Not a valid output format"; return false;}

      /// @brief The "Convert" interface Write function.

      /// Writes a single object
      /// Deletes the object after writing 
      /// \return false on error
      virtual bool WriteChemObject(OBConversion* /*pConv*/)
        { std::cerr << "Not a valid output format"; return false;}

      /// @brief Information on this format. Printed out in response to -Hxxx option where xxx id the id of the format.

      /// Must be provided by each format class.
      /// Can include a list of command line Options. These may be used to construction
      /// check boxes, radio buttons etc for GUI interface.
      virtual const char* Description()=0;
	
      /// @brief A decription of the chemical object converted by this format.

      /// If not provided, the object type used by the default format is used (usually OBMol). 
      virtual const char* TargetClassDescription();

      /// \return the type of chemical object used by the format.

      /// Defaults to that used by the default format. Useful for checking 
      /// that a format can handle a particular object.
      virtual const std::type_info& GetType();
 	
      /// @brief Web address where the format is defined.
      virtual const char* SpecificationURL() { return ""; }

      /// @brief Chemical MIME type associated with this file type (if any)
      virtual const char* GetMIMEType() { return ""; }

      /// @brief Decribes the capabilities of the format (Read only etc.)
   
      /// Currently, can be a bitwise OR of any of the following
      /// NOTREADABLE READONEONLY NOTWRITABLE WRITEONEONLY DEFAULTFORMAT
      /// READBINARY WRITEBINARY READXML
      virtual unsigned int Flags() { return 0;}; 

      /// @brief Skip past first n objects in input stream (or current one with n=0)

      /// \return 1 on success, -1 on error and 0 if not implemented 
      virtual int SkipObjects(int /*n*/, OBConversion* /*pConv*/)
        {
          return 0; //shows not implemented in the format class
        };

      /// \return a pointer to a new instance of the format, or NULL if fails.

      /// Normally a single global instance is used but this may cause problems
      /// if there are member variables and the format is used in more than one place
      /// in the program.
      virtual OBFormat* MakeNewInstance()
        {
          return NULL; //shows not implemented in the format class
        }

      /// @brief Format classes do not have a destructor
      virtual ~OBFormat(){};
    };

  //*************************************************
  /// @struct CharPtrLess obconversion.h <openbabel/obconversion.h>
  /// @brief Case insensitive string comparison for FormatsMap key.
	struct CharPtrLess : public std::binary_function<const char*,const char*, bool>
	{
		bool operator()(const char* p1,const char* p2) const
		{ return strcasecmp(p1,p2)<0; }
	};

	typedef std::map<const char*,OBFormat*,CharPtrLess > FMapType;
	typedef FMapType::iterator Formatpos;

  //*************************************************
  /// @brief Class to convert from one format to another.
  // Class introduction in obconversion.cpp
  class OBCONV OBConversion
    {
      /// @nosubgrouping
    public:
      /// @name Construction
      //@{
      OBConversion(std::istream* is=NULL, std::ostream* os=NULL);
      /// @brief Copy constructor
      OBConversion(const OBConversion& o);
      virtual     ~OBConversion(); 
      //@}	
      /// @name Collection of formats
      //@{
      /// @brief Called once by each format class
      static int				RegisterFormat(const char* ID, OBFormat* pFormat, const char* MIME = NULL);
      /// @brief Searches registered formats
      static OBFormat*	FindFormat(const char* ID);
      /// @brief Searches registered formats for an ID the same as the file extension
      static OBFormat*	FormatFromExt(const char* filename);
      /// @brief Searches registered formats for a MIME the same as the chemical MIME type passed
      static OBFormat*        FormatFromMIME(const char* MIME);

      ///Repeatedly called to recover available Formats
      static bool	        GetNextFormat(Formatpos& itr, const char*& str,OBFormat*& pFormat);
      //@}
		
      /// @name Information
      //@{
      static const char* Description(); //generic conversion options
      //@}

      /// @name Parameter get and set
      //@{
      std::istream* GetInStream() const {return pInStream;};
      std::ostream* GetOutStream() const {return pOutStream;};
      void          SetInStream(std::istream* pIn)
        { 
          if (pInStream && NeedToFreeInStream) {
            delete pInStream; NeedToFreeInStream = false;
          }
          pInStream=pIn;
          CheckedForGzip = false; // haven't tried to gzip decode this stream
        };
      void          SetOutStream(std::ostream* pOut)
        {
          if (pOutStream && NeedToFreeOutStream) {
            delete pOutStream; NeedToFreeOutStream = false;
          }
          pOutStream=pOut;
        };
      /// Sets the formats from their ids, e g CML
      bool        SetInAndOutFormats(const char* inID, const char* outID);
      bool        SetInAndOutFormats(OBFormat* pIn, OBFormat* pOut);
      /// Sets the input format from an id e.g. CML
      bool	      SetInFormat(const char* inID);
      bool	      SetInFormat(OBFormat* pIn);
      /// Sets the output format from an id e.g. CML
      bool	      SetOutFormat(const char* outID);
      bool	      SetOutFormat(OBFormat* pOut);

      OBFormat*   GetInFormat() const{return pInFormat;};
      OBFormat*   GetOutFormat() const{return pOutFormat;};
      std::string GetInFilename() const{return InFilename;};
	
      ///Get the position in the input stream of the object being read
      std::streampos GetInPos()const{return wInpos;}; 

      ///Get the length in the input stream of the object being read
      size_t GetInLen()const{return wInlen;}; 

      /// \return a default title which is the filename
      const char* GetTitle() const;

      ///@brief Extension method: deleted in ~OBConversion()
      OBConversion* GetAuxConv() const {return pAuxConv;};
      void          SetAuxConv(OBConversion* pConv) {pAuxConv=pConv;};
      //@}
      /// @name Option handling
      //@{
      ///@brief Three types of options set on the the command line by -a? , -x? , or -?
      enum Option_type { INOPTIONS, OUTOPTIONS, GENOPTIONS };

      ///@brief Determine whether an option is set. \return NULL if option not and a pointer to the associated text if it is 
      const char* IsOption(const char* opt,Option_type opttyp=OUTOPTIONS);
	
      ///@brief Access the map with option name as key and any associated text as value
      const std::map<std::string,std::string>* GetOptions(Option_type opttyp)
        { return &OptionsArray[opttyp];};

      ///@brief Set an option of specified type, with optional text
      void AddOption(const char* opt, Option_type opttyp, const char* txt=NULL);
	
      bool RemoveOption(const char* opt, Option_type optype);

      ///@brief Set several single character options of specified type from string like ab"btext"c"ctext"
      void SetOptions(const char* options, Option_type opttyp);

      ///@brief For example -h takes 0 parameters; -f takes 1. Call in a format constructor.
      static void RegisterOptionParam(std::string name, OBFormat* pFormat,
                                      int numberParams=0, Option_type typ=OUTOPTIONS);

      /// \return the number of parameters registered for the option, or 0 if not found
      static int GetOptionParams(std::string name, Option_type typ);
      //@}

      /// @name Supported file format 
      //@{
      // @brief Set and return the list of supported input format
      std::vector<std::string> GetSupportedInputFormat();
      // @brief Set and return the list of supported output format
      std::vector<std::string> GetSupportedOutputFormat();
      //@}

      /// @name Conversion
      //@{
      /// @brief Conversion for single input and output stream
      int         Convert(std::istream* is, std::ostream* os);

      /// @brief Conversion with existing streams
      int         Convert();

      /// @brief Conversion with multiple input/output files:
      /// makes input and output streams, and carries out normal, batch, aggregation, and splitting conversion.
      int					FullConvert(std::vector<std::string>& FileList,
                              std::string& OutputFileName, std::vector<std::string>& OutputFileList);
      //@}

      /// @name Conversion loop control
      //@{
      bool				AddChemObject(OBBase* pOb);///< @brief Adds to internal array during input
      OBBase*			GetChemObject(); ///< @brief Retrieve from internal array during output
      bool				IsLast();///< @brief True if no more objects to be output
      bool				IsFirstInput();///< @brief True if the first input object is being processed
      int         GetOutputIndex() const ;///< @brief Retrieves number of ChemObjects that have been actually output
      void				SetOutputIndex(int indx);///< @brief Sets ouput index (maybe to control whether seen as first object)
      void				SetMoreFilesToCome();///<@brief Used with multiple input files. Off by default.
      void				SetOneObjectOnly(bool b=true);///<@brief Used with multiple input files. Off by default.
      void        SetLast(bool b){SetOneObjectOnly(b);}///@brief.Synonym for SetOneObjectOnly()
      //@}
      /// @name Convenience functions
      //@{
      ///The default format is set in a single OBFormat class (generally it is OBMol) 
      static OBFormat* GetDefaultFormat(){return pDefaultFormat;};

      /// @brief Outputs an object of a class derived from OBBase.
	
      /// Part of "API" interface. 
      /// The output stream can be specified and the change is retained in the OBConversion instance
      bool				Write(OBBase* pOb, std::ostream* pout=NULL);

      /// @brief Outputs an object of a class derived from OBBase as a string
	
      /// Part of "API" interface. 
      /// The output stream is temporarily changed to the string and then restored
      /// This method is primarily intended for scripting languages without "stream" classes
      /// The optional "trimWhitespace" parameter allows trailing whitespace to be removed
      /// (e.g., in a SMILES string or InChI, etc.)
      std::string                     WriteString(OBBase* pOb, bool trimWhitespace = false);

      /// @brief Outputs an object of a class derived from OBBase as a file (with the supplied path)
	
      /// Part of "API" interface. 
      /// The output stream is changed to the supplied file and the change is retained in the
      /// OBConversion instance.
      /// This method is primarily intended for scripting languages without "stream" classes
      bool                            WriteFile(OBBase* pOb, std::string filePath);

      /// @brief Manually closes and deletes the output stream
      /// The file is closed anyway when in the OBConversion destructor or when WriteFile
      /// is called again.
      /// \since version 2.1
      void CloseOutFile();

      /// @brief Reads an object of a class derived from OBBase into pOb.
	
      /// Part of "API" interface. 
      /// The input stream can be specified and the change is retained in the OBConversion instance
      /// \return false and pOb=NULL on error 
      bool	Read(OBBase* pOb, std::istream* pin=NULL);

      /// @brief Reads an object of a class derived from OBBase into pOb from the supplied string
	
      /// Part of "API" interface. 
      /// \return false and pOb=NULL on error
      /// This method is primarily intended for scripting languages without "stream" classes
      bool	ReadString(OBBase* pOb, std::string input);

      /// @brief Reads an object of a class derived from OBBase into pOb from the file specified
	
      /// Part of "API" interface. 
      /// The output stream is changed to the supplied file and the change is retained in the
      /// OBConversion instance.
      /// \return false and pOb=NULL on error 
      /// This method is primarily intended for scripting languages without "stream" classes
      bool	ReadFile(OBBase* pOb, std::string filePath);

protected:
      ///Replaces * in BaseName by InFile without extension and path
      static std::string BatchFileName(std::string& BaseName, std::string& InFile);
      ///Replaces * in BaseName by Count
      static std::string IncrementedFileName(std::string& BaseName, const int Count);
      ///Checks for misunderstandings when using the -m option
      static bool CheckForUnintendedBatch(const std::string& infile, const std::string& outfile);
      ///Adds a filtering rdbuffer to handle line endings if not already installed and not a binary or xml format.
      void InstallStreamFilter();

      //@}

    protected:
      bool             SetStartAndEnd();
      static FMapType& FormatsMap();///<contains ID and pointer to all OBFormat classes
      static FMapType& FormatsMIMEMap();///<contains MIME and pointer to all OBFormat classes
      typedef std::map<std::string,int> OPAMapType;
      static OPAMapType& OptionParamArray(Option_type typ);
      static int       LoadFormatFiles();
      bool             OpenAndSetFormat(bool SetFormat, std::ifstream* is);

      std::string	  InFilename;
      std::istream*     pInStream;
      std::ostream*     pOutStream;
      static OBFormat*  pDefaultFormat;
      OBFormat* 	  pInFormat;
      OBFormat*	  pOutFormat;

      std::map<std::string,std::string> OptionsArray[3];

      int		  Index;
      unsigned int	  StartNumber;
      unsigned int	  EndNumber;
      int	          Count;
      bool			m_IsFirstInput;
      bool		  m_IsLast;
      bool		  MoreFilesToCome;
      bool		  OneObjectOnly;
      bool		  ReadyToInput;
      bool      CheckedForGzip;      ///< input stream is gzip-encoded
      bool      NeedToFreeInStream;
      bool      NeedToFreeOutStream;
      typedef   FilteringInputStreambuf< LineEndingExtractor > LErdbuf;
      LErdbuf*  pLineEndBuf;

      static int FormatFilesLoaded;
      OBBase*		  pOb1;
      std::streampos wInpos; ///<position in the input stream of the object being written
      std::streampos rInpos; ///<position in the input stream of the object being read
      size_t wInlen; ///<length in the input stream of the object being written
      size_t rInlen; ///<length in the input stream of the object being read
	
      OBConversion* pAuxConv;///<Way to extend OBConversion

      std::vector<std::string> SupportedInputFormat; ///< list of supported input format
      std::vector<std::string> SupportedOutputFormat; ///< list of supported output format

    };

  ///For OBFormat::Flags()
#define NOTREADABLE     0x01
#define READONEONLY     0x02
#define READBINARY      0x04
#define ZEROATOMSOK     0x08
#define NOTWRITABLE     0x10
#define WRITEONEONLY    0x20
#define WRITEBINARY     0x40
#define READXML         0x80
#define DEFAULTFORMAT   0x4000

} //namespace OpenBabel
#endif //OB_CONV_H

//! \file
//! \brief Handle file conversions. Declaration of OBFormat, OBConversion.

 

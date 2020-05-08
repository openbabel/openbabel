/**********************************************************************
obconversion.h - Handle file conversions. Declaration of OBFormat, OBConversion

Copyright (C) 2004-2009 by Chris Morley

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

#ifndef OB_CONV_H
#define OB_CONV_H

#include <openbabel/babelconfig.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>

#include <string>
#include <vector>
#include <map>
#ifdef HAVE_STRINGS_H
#include <strings.h>
#endif

#include <openbabel/oberror.h>
#include <openbabel/format.h>
#include <openbabel/lineend.h>

// These macros are used in DLL builds. If they have not
// been set in babelconfig.h, define them as nothing.
#ifndef OBCONV
	#define OBCONV
#endif
#ifndef OBDLL
	#define OBDLL
#endif

namespace OpenBabel {

  // Needed to preserve deprecated API
  typedef OBPlugin::PluginIterator Formatpos;

  OBERROR extern  OBMessageHandler obErrorLog;

  //*************************************************
  /// @brief Class to convert from one format to another.
  // Class introduction in obconversion.cpp
  class OBCONV OBConversion
    {
      /// @nosubgrouping
    public:
      /// @name Construction
      //@{
      OBConversion(std::istream* is=nullptr, std::ostream* os=nullptr);
      OBConversion(std::string inFilename, std::string outFilename="");
      /// @brief Copy constructor.  Stream *ownership* is not copied. Source remains responsible for the memory.
      OBConversion(const OBConversion& o);
      /// @brief Assignment.  Stream *ownership* is not copied.  Source remains responsible for the memory.
      OBConversion& operator=(const OBConversion& rhs);

      virtual     ~OBConversion();
      //@}
      /// @name Collection of formats
      //@{
      /// @brief Called once by each format class
      static int				RegisterFormat(const char* ID, OBFormat* pFormat, const char* MIME = nullptr);
      /// @brief Searches registered formats
      static OBFormat*	FindFormat(const char* ID);
      /// @brief Searches registered formats
      /// \since version 2.3
      static OBFormat*  FindFormat(const std::string ID);
      /// @brief Searches registered formats for an ID the same as the file extension
      static OBFormat*	FormatFromExt(const char* filename);
      static OBFormat*	FormatFromExt(const char* filename, bool& isgzip);
      /// @brief Searches registered formats for an ID the same as the file extension
      /// \since version 2.3
      static OBFormat*	FormatFromExt(const std::string filename);
      static OBFormat*	FormatFromExt(const std::string filename, bool& isgzip);
      /// @brief Searches registered formats for a MIME the same as the chemical MIME type passed
      static OBFormat*        FormatFromMIME(const char* MIME);

      ///Deprecated!.Repeatedly called to recover available Formats
#ifndef SWIG
      static bool	        GetNextFormat(Formatpos& itr, const char*& str,OBFormat*& pFormat);
#endif
      //@}

      /// @name Information
      //@{
      static const char* Description(); //generic conversion options
      //@}

      /// These return a filtered stream for reading/writing (possible filters include compression, decompression, and newline transformation)
      /// @name Parameter get and set
      //@{
      std::istream* GetInStream() const {return pInput;};
      std::ostream* GetOutStream() const {return pOutput;};

      /// @brief Set input stream.  If takeOwnership is true, will deallocate when done.
      /// If isGzipped is true, will treat as a gzipped stream regardless of option settings,
      //  if false, then will be treated as gzipped stream only if z/zin is set.
      void          SetInStream(std::istream* pIn, bool takeOwnership=false);
      void          SetOutStream(std::ostream* pOut, bool takeOwnership=false);

      /// Sets the formats from their ids, e g CML
      bool        SetInAndOutFormats(const char* inID, const char* outID, bool ingzip=false, bool outgzip=false);
      bool        SetInAndOutFormats(OBFormat* pIn, OBFormat* pOut, bool ingzip=false, bool outgzip=false);
      /// Sets the input format from an id e.g. CML
      bool	      SetInFormat(const char* inID, bool isgzip=false);
      bool	      SetInFormat(OBFormat* pIn, bool isgzip=false);
      /// Sets the output format from an id e.g. CML
      bool	      SetOutFormat(const char* outID, bool isgzip=false);
      bool	      SetOutFormat(OBFormat* pOut, bool isgzip=false);

      OBFormat*   GetInFormat() const{return pInFormat;};
      OBFormat*   GetOutFormat() const{return pOutFormat;};
      bool GetInGzipped() const{return inFormatGzip;};
      bool GetOutGzipped() const{return outFormatGzip;};
      std::string GetInFilename() const{return InFilename;};
      std::string GetOutFilename() const{return OutFilename;};

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
      /** @name Option handling
       Three types of Option provide information and control instructions to the
       conversion process, INOPTIONS, OUTOPTIONS, GENOPTIONS, and are stored in each
       OBConversion object in separate maps. Each option has an id and an optional
       text string. They are set individually by AddOption() or (rarely) collectively
       in SetOptions(). Options cannot be altered but can be replaced with AddOption()
       and deleted with RemoveOption(), which, however, should be used in an op derived
       from OBOp (because of iterator invalidation).

       If the "Convert" interface is used, the GENOPTIONS are acted upon in the
       OBBase::DoTransformations() functions (currently only OBMol has one). This
       happens after the object has been input but before it has been output.
       All the options are available to input and output formats, etc. via the IsOption()
       function, and the interpretation of any text string needs to be done subsequently.

       In the commandline interface, options with single character ids are are indicated
       like -s, and those with multiple character ids like --gen3D. An option may have
       one or more parameters which appear, space separated, in the option's text string.
       With babel, unless the option is at the end of the command, it is necessary for
       the number of its parameters to be exactly that specified in RegisterOptionParam().
       The default is 0, but if it is more, and babel is likely to be used, this function
       should be called in the constructor of a format or op.
       With obabel (or the GUI), it is not necessary to call RegisterOptionParam().

       New GENOPTIONS can be defined (as plugins) using the class OBOp.

       It is customary for a format or op to document any INOPTIONS or OUTPTIONS it
       uses in its Description() function. As well as providing documentation during
       use, this is also parsed by the GUI to construct its checkboxes,etc., so it is
       advisable to give new Descriptions the same form as existing ones.

       Some conversion options, such as -f, -l, -m, are unlikely to be used in
       programming, but are listed in OBConversion::Description().  The built-in
       GENOPTIONS for OBMol objects are listed in OBMol::ClassDescription() which
       is in transform.cpp and also in this documentation under AddOption().
       */
      //@{
      ///@brief Three types of options set on the the command line by -a? , -x? , or -?
      enum Option_type { INOPTIONS, OUTOPTIONS, GENOPTIONS, ALL };

      ///@brief Determine whether an option is set. \return NULL if option not and a pointer to the associated text if it is
      const char* IsOption(const char* opt,Option_type opttyp=OUTOPTIONS);

      ///@brief Access the map with option name as key and any associated text as value
      const std::map<std::string,std::string>* GetOptions(Option_type opttyp)
        { return &OptionsArray[opttyp];};

      ///@brief Set an option of specified type, with optional text
      void AddOption(const char* opt, Option_type opttyp=OUTOPTIONS, const char* txt=nullptr);

      bool RemoveOption(const char* opt, Option_type optype);

      ///@brief Set several single character options of specified type from string like ab"btext"c"ctext"
      void SetOptions(const char* options, Option_type opttyp);

      ///@brief For example -h takes 0 parameters; -f takes 1. Call in a format constructor.
      static void RegisterOptionParam(std::string name, OBFormat* pFormat,
                                      int numberParams=0, Option_type typ=OUTOPTIONS);

      /// \return the number of parameters registered for the option, or 0 if not found
      static int GetOptionParams(std::string name, Option_type typ);
      //@}

      ///@brief Copies the options (by default of all types) from one OBConversion Object to another.
      void CopyOptions(OBConversion* pSourceConv, Option_type typ=ALL);

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
      int     AddChemObject(OBBase* pOb);///< @brief Adds to internal array during input
      OBBase*  GetChemObject(); ///< @brief Retrieve from internal array during output
      bool     IsLast();///< @brief True if no more objects to be output
      bool     IsFirstInput();///< @brief True if the first input object is being processed
      void     SetFirstInput(bool b=true);///< @brief Setwhether or not is the first input
      int      GetOutputIndex() const ;///< @brief Retrieves number of ChemObjects that have been actually output
      void     SetOutputIndex(int indx);///< @brief Sets output index (maybe to control whether seen as first object)
      void     SetMoreFilesToCome();///<@brief Used with multiple input files. Off by default.
      void     SetOneObjectOnly(bool b=true);///< @brief Used with multiple input files. Off by default.
      void     SetLast(bool b){SetOneObjectOnly(b);}///< @brief Synonym for SetOneObjectOnly()
      bool     IsLastFile(){ return !MoreFilesToCome;}///< @brief True if no more files to be read
      /// @brief Number of objects read and processed
      /// Incremented after options are processed, so 0 for first object.  Returns -1 if Convert interface not used. 
      int      GetCount()const { return Count; }
      //@}
      /// @name Convenience functions
      //@{
      ///The default format is set in a single OBFormat class (generally it is OBMol)
      static OBFormat* GetDefaultFormat(){return OBFormat::FindType(nullptr);};

      /// @brief Outputs an object of a class derived from OBBase.

      /// Part of "API" interface.
      /// The output stream can be specified and the change is retained in the OBConversion instance
      bool				Write(OBBase* pOb, std::ostream* pout=nullptr);

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
      bool	Read(OBBase* pOb, std::istream* pin=nullptr);

      /// Part of "API" interface.
      /// The input stream can be specified and the change is retained in the OBConversion instance
      /// \return NULL on error
//      OBBase*	ReadObject(std::istream* pin=NULL);

      /// @brief Reads an object of a class derived from OBBase into pOb from the supplied string

      /// Part of "API" interface.
      /// \return false and pOb=NULL on error
      /// This method is primarily intended for scripting languages without "stream" classes
      /// Any existing input stream will be replaced by stringstream.
      bool	ReadString(OBBase* pOb, std::string input);

      /// @brief Reads an object of a class derived from OBBase into pOb from the file specified

      /// Part of "API" interface.
      /// The output stream is changed to the supplied file and the change is retained in the
      /// OBConversion instance. For multi-molecule files, the remaining molecules
      /// can be read by repeatedly calling the Read() method.
      /// \return false and pOb=NULL on error
      /// This method is primarily intended for scripting languages without "stream" classes
      bool	ReadFile(OBBase* pOb, std::string filePath);

      /// Part of the "Convert" interface.
      /// Open the files and update the streams in the OBConversion object.
      /// This method is primarily intended for scripting languages without "stream" classes
      /// and will usually followed by a call to Convert().
      /// Will set format from file extension if format has not already been set.
      /// Files will be opened even if format cannot be determined, but not if file path is empty.
      /// \return false if unsuccessful.
      bool OpenInAndOutFiles(std::string infilepath, std::string outfilepath);

      /// @brief Sends a message like "2 molecules converted" to clog
      /// The type of object is taken from the TargetClassDescription
      /// of the specified class (or the output format if not specified)and
      /// is appropriately singular or plural.
      void ReportNumberConverted(int count, OBFormat* pFormat=nullptr);

      /// \return the number of objects in the inputstream,
      /// or -1 if error or if SkipObjects for the input format is not implemented
      /// Adjusts for the value of -f and -l options (first and last objects).
      int NumInputObjects();


protected:
      ///Replaces * in BaseName by InFile without extension and path
      static std::string BatchFileName(std::string& BaseName, std::string& InFile);
      ///Replaces * in BaseName by Count
      static std::string IncrementedFileName(std::string& BaseName, const int Count);
      ///Checks for misunderstandings when using the -m option
      static bool CheckForUnintendedBatch(const std::string& infile, const std::string& outfile);

      void ClearInStreams();
      //@}

    protected:

      //helper class for saving stream state
      struct StreamState
      {
          std::ios *pStream; //active stream
          std::vector<std::ios *> ownedStreams; //streams we own the memory to

          StreamState(): pStream(nullptr) {}
          ~StreamState()
          {
            assert(ownedStreams.size() == 0); //should be popped
          }

          void pushInput(OBConversion& conv);
          void popInput(OBConversion& conv);

          void pushOutput(OBConversion& conv);
          void popOutput(OBConversion& conv);

          bool isSet() const { return pStream != nullptr; }
      };

      bool             SetStartAndEnd();
//      static FMapType& FormatsMap();///<contains ID and pointer to all OBFormat classes
//      static FMapType& FormatsMIMEMap();///<contains MIME and pointer to all OBFormat classes
      typedef std::map<std::string,int> OPAMapType;
      static OPAMapType& OptionParamArray(Option_type typ);
      bool             OpenAndSetFormat(bool SetFormat, std::ifstream* is, std::stringstream* ss=nullptr);

      std::string	  InFilename, OutFilename; //OutFileName added v2.4.0

      typedef   FilteringInputStream< LineEndingExtractor > LEInStream;

      std::istream *pInput; //input stream, may be filtered
      std::vector<std::istream *> ownedInStreams; //streams we own the memory to

      std::ostream *pOutput; //output stream, may have filters applied
      std::vector<std::ostream *> ownedOutStreams; //streams we own the memory to


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
      bool      SkippedMolecules;    /// skip molecules using -f and -l

      //unlike the z and zin options, these are not sticky - setting formats will reset them
      bool inFormatGzip;
      bool outFormatGzip;

      OBBase*		  pOb1;
      std::streampos wInpos; ///<position in the input stream of the object being written
      std::streampos rInpos; ///<position in the input stream of the object being read
      size_t wInlen; ///<length in the input stream of the object being written
      size_t rInlen; ///<length in the input stream of the object being read

      OBConversion* pAuxConv;///<Way to extend OBConversion

      std::vector<std::string> SupportedInputFormat; ///< list of supported input format
      std::vector<std::string> SupportedOutputFormat; ///< list of supported output format

    };

} //namespace OpenBabel
#endif //OB_CONV_H

//! \file
//! \brief Handle file conversions. Declaration of OBFormat, OBConversion.



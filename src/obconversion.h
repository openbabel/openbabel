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
#ifndef OB_CONV_H
#define OB_CONV_H
//Declaration of OBFormat and OBConversion

#include "babelconfig.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include "dlhandler.h"

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
//*************************************************

/// @brief Base class for file formats.

/// Two sets of Read and Write functions are specified for each format
/// to handle two different requirements.
/// The "Convert" interface is for use in file format conversion applications. The
/// user interface, a console, a GUI, or another program is kept unaware of the
/// details of the chemistry and does not need to #include mol.h. It is then
/// necessary to manipulate only pointers to OBBase in OBConversion and the user
/// interface, with all the construction and deletion of OBMol etc objects being
/// done in the Format classes or the OB core. The convention  with "Covert"
/// interface functions is that chemical objects are made on the heap with new
/// in the ReadChemicalObject() functions and and deleted in WriteChemicalObject()
/// functions
/// 
/// The "API" interface is for programatic use of the OB routines in application
/// programs where mol.h is #included. There is generally no creation or
/// destruction of objects in ReadMolecule() and WriteMolecule() and no restriction
/// on whether the pointers are to the heap or the stack.
/// 
class OBCONV OBFormat
{
public:
	/// @brief The "API" interface Read function.

	/// Reads a single object.
	/// Does not make a new object on the heap; 
	/// can be used with a pointer to an chem object on the heap or the stack.
	virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv)
		{ std::cerr << "Not a valid input format"; return false;}

	/// @brief The "Convert" interface Read function.

	/// Possibly reads multiple new objects on the heap and subjects them 
	/// to its DoTransformations() function, which may delete them again. 
	/// Sends result to pConv->AddChemObject()
	virtual bool ReadChemObject(OBConversion* pConv)
		{ std::cerr << "Not a valid input format"; return false;}

	/// @brief The "API" interface Write function.

	/// Writes a single object
	/// Does not delete the object; 
	/// can be used with a pointer to an chem object on the heap or the stack.
	/// Returns false on error.
	virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv)
		{ std::cerr << "Not a valid output format"; return false;}

	/// @brief The "Convert" interface Write function.

	/// Writes a single object
	/// Deletes the object after writing 
	/// Returns false on error
	virtual bool WriteChemObject(OBConversion* pConv)
		{ std::cerr << "Not a valid output format"; return false;}

	/// @brief Information on this format. Printed out in response to -Hxxx option where xxx id the id of the format.

	/// Must be provided by each format class.
	/// Can include a list of command line Options. These may be used to construction
	/// check boxes, radio buttons etc for GUI interface.
	virtual const char* Description()=0;
	
	/// @brief A decription of the chemical object converted by this format.

	/// If not provided, the object type used by the default format is used (usually OBMol). 
	virtual const char* TargetClassDescription();

	/// @brief Returns the type of chemical object used by the format.

	/// Defaults to that used by the default format. Useful for checking 
	/// that a format can handle a particular object.
	virtual const std::type_info& GetType();
 	
	/// @brief Web address where the format is defined.

	/// (In case somebody wants to do something fancy.)
virtual const char* SpecificationURL();

/// @brief Decribes the capabilities of the format (Read only etc.)
   
/// Currently, can be a bitwise OR of any of the following
/// NOTREADABLE READONEONLY NOTWRITABLE WRITEONEONLY DEFAULTFORMAT
	virtual unsigned int Flags() { return 0;}; 
};

//*************************************************
/// @brief Case insensitive string comparison for FormatsMap key.
	struct CharPtrLess : public std::binary_function<const char*,const char*, bool>
	{
		bool operator()(const char* p1,const char* p2) const
		{ return strcasecmp(p1,p2)<0; }
	};

	typedef std::map<const char*,OBFormat*,CharPtrLess > FMapType;
	typedef FMapType::iterator Formatpos;

//*************************************************

/// Class to convert from one format to another.
class OBCONV OBConversion
{
	/// @nosubgrouping
public:
	/// @name Construction
	//@{
							OBConversion(std::istream* is=NULL, std::ostream* os=NULL);
	/// @brief Copy constructor
							OBConversion(const OBConversion& Conv);
	virtual     ~OBConversion(); 
	//@}	
	/// @name Collection of formats
	//@{
	/// @brief Called once by each format class
	static int				RegisterFormat(const char* ID, OBFormat* pFormat);
	/// @brief Searches registered formats
	static OBFormat*	FindFormat(const char* ID);
	/// @brief Searches registered formats for an ID the same as the file extension
	static OBFormat*	FormatFromExt(const char* filename);

	///Repeatedly called to recover available Formats
	static bool				GetNextFormat(Formatpos& itr, const char*& str,OBFormat*& pFormat);
	//@}		
	/// @name Information
	//@{

	static const char* Description(); //generic conversion options
	//@}
	/// @name Parameter get and set
	//@{
	std::istream*		GetInStream() const {return pInStream;};
	std::ostream*		GetOutStream() const {return pOutStream;};
	bool        SetInAndOutFormats(const char* inID, const char* outID);///< Sets the formats from their ids, e g CML
	bool        SetInAndOutFormats(OBFormat* pIn, OBFormat* pOut);
	OBFormat*		GetInFormat() const{return pInFormat;};
	OBFormat*		GetOutFormat() const{return pOutFormat;};
	
	const char* GetOptions() const;///< @brief Options for output format
	void        SetOptions(const char* options);///< @brief Options for output format

	///@brief Returns NULL if option not set and a valid pointer if it is 
	const char* OBConversion::IsOption(const char ch) const;

	const char* GetGeneralOptions() const;///< @brief Options related to converted chemical object type 
	void        SetGeneralOptions(const char* options);////< @brief Options related to converted chemical object type
	///
	const char* GetTitle() const;
	void        SetTitle(const char* title);
	const char* GetDimension() const;
	void				SetDimension(const char* dim);

	bool				SaveOptionsToFile(const char* filename);
	bool				RestoreOptionsFromFile(const char* filename);
	//@}
	/// @name Conversion
	//@{
	/// @brief Conversion for single input and output stream
	int         Convert(std::istream* is, std::ostream* os);

	/// @brief Conversion for single input and output stream
	int         Convert();

	/// @brief Conversion with multiple input/output files:
	/// makes input and output streams, and carries out normal, batch, aggregation, and splitting conversion.
	int					FullConvert(std::vector<std::string>& FileList,
										std::string& OutputFileName, std::vector<std::string>& OutputFileList);
	//@}
	/// @name Conversion loop control
	//@{
	int					AddChemObject(OBBase* pOb);///< @brief Adds to internal array during input
	OBBase*			GetChemObject(); ///< @brief Retrieve from internal array during output
	bool				IsLast();///< @brief True if no more objects to be output
	int         GetOutputIndex() const ;///< @brief Retrieves number of ChemObjects that have been actually output
	void				SetOutputIndex(int indx);///< @brief Sets ouput index (maybe to control whether seen as first object)
	void				SetMoreFilesToCome();///<@brief Used with multiple input files. Off by default.
	void				SetOneObjectOnly();///<@brief Used with multiple input files. Off by default.
	//@}
	/// @name Convenience functions
	//@{
	///The default format is set in a single OBFormat class (generally it is OBMol) 
	static OBFormat* GetDefaultFormat(){return pDefaultFormat;};

	/// @brief Outputs an object of a class derived from OBBase.
	
	/// Part of "API" interface. 
	/// The output stream can be specified and the change is retained in the OBConversion instance
	bool				Write(OBBase* pOb, std::ostream* pout=NULL);

	/// @brief Reads an object of a class derived from OB base into pOb.
	
	/// Part of "API" interface. 
	/// The input stream can be specified and the change is retained in the OBConversion instance
	/// Returns false and pOb=NULL on error 
	template<class T> 
		bool	Read(T* pOb, std::istream* pin=NULL)
	{
		if(pin)
			pInStream=pin;
		if(!pInFormat) return false;
		if(!pInFormat->ReadMolecule(pOb, this))
		{pOb=NULL; return false;}
		pOb = dynamic_cast<T*>(pOb);
		return (pOb!=NULL);
	};

	///Replaces * in BaseName by InFile without extension and path
	static std::string BatchFileName(std::string& BaseName, std::string& InFile);
	///Replaces * in BaseName by Count
	static std::string IncrementedFileName(std::string& BaseName, const int Count);
	//@}

protected:
	void             SetStartAndEnd();
	static FMapType& FormatsMap();///<contains ID and pointer to all OBFormat classes
	static int       LoadFormatFiles();


	std::istream*     pInStream;
	std::ostream*     pOutStream;
	static OBFormat*  pDefaultFormat;
	OBFormat* 		pInFormat;
	OBFormat*			pOutFormat;
	std::string		Options;
	std::string		GeneralOptions;
	std::string		Title;
	char					Dimension[10];
	int						Index;
	unsigned int	StartNumber;
	unsigned int	EndNumber;
	int	Count;
	bool					m_IsLast;
	bool					MoreFilesToCome;
	bool					OneObjectOnly;
	bool					ReadyToInput;
	static bool		FormatFilesLoaded;
	OBBase*				pOb1;
};

///For OBFormat::Flags()
#define NOTREADABLE     0x01
#define READONEONLY     0x02
#define NOTWRITABLE     0x10
#define WRITEONEONLY    0x20
#define DEFAULTFORMAT 0x4000

} //namespace OpenBabel
#endif //OB_CONV_H

/**********************************************************************
xml.h Declaration of XMLConversion, 
declaration and definition of XMLBaseFormat and XMLMoleculeFormat 
Copyright (C) 2005 by Chris Morley
 
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
#include "libxml/xmlreader.h"
#include "libxml/xmlwriter.h"

namespace OpenBabel
{

//global utility in xml.cpp;
streamsize gettomatch(istream& is, char* buf, streamsize count, const char* match);
//forward declaration
class XMLBaseFormat;

//******************************************************
//XMLConversion class

/// An extended OBConversion class which includes a libxml2 reader for use with XML formats
/// Copies an OBConversion and then extends it with a XML parser. 
/// Instances made on the heap are deleted when the original OBConversion object is.
class XMLConversion : public OBConversion
{
public:
	///Existing OBConversion instance copied
		XMLConversion(OBConversion* pConv);
		
		///Frees reader and writer if necessary
		~XMLConversion();

	bool XMLConversion::SetupReader();///< opens libxml2 reader
	bool XMLConversion::SetupWriter();///< opens libxml2 writer

	///Parses the input xml stream and sends each element to the format's callback routines
	bool ReadXML(XMLBaseFormat* pFormat, OBBase* pOb);

	typedef std::map<std::string, XMLBaseFormat*> NsMapType;

	///This static function returns a reference to the map
	///Avoids "static initialization order fiasco"
	static NsMapType& Namespaces()
	{
		static NsMapType* nsm = new NsMapType;
		return *nsm;
	};

	static void XMLConversion::RegisterXMLFormat(XMLBaseFormat* pFormat,
			bool IsDefault=false, const char* uri=NULL);

	///Returns the extended OBConversion class, making it if necessary
	static XMLConversion* GetDerived(OBConversion* pConv, bool ForReading=true);

	///Because OBConversion::Convert is still using the unextended OBConversion object
	///we need to obtain the conversion paramters from it when requested
	bool XMLConversion::IsLast()
	{ return _pConv->IsLast(); }
	int XMLConversion::GetOutputIndex()
	{ return  _pConv->GetOutputIndex(); }


	xmlTextReaderPtr GetReader() const
	{		return _reader;	};

	xmlTextWriterPtr GetWriter() const
	{		return _writer;	};

	void OutputToStream()
	{
		int ret=xmlOutputBufferFlush(_buf);
	}

	static XMLBaseFormat* GetDefaultXMLClass() //TODO make dependent on object type
	{   return _pDefault;};

	void LookForNamespace()
	{   _LookingForNamespace = true; };

	///Static callback functions for xmlReaderForIO()
	static int ReadStream(void * context, char * buffer, int len);
	static int WriteStream(void * context, const char * buffer, int len);
	//static int CloseStream(void* context);

	string GetAttribute(const char* attrname);

	///Sets value to element content. Returns false if there is no content. 
	string XMLConversion::GetContent();

	///Sets value to element content as an integer. Returns false if there is no content. 
	bool    GetContentInt(int& value);

	///Sets value to element content as an double. Returns false if there is no content. 
	bool GetContentDouble(double& value);

private:
	static XMLBaseFormat* _pDefault;
	OBConversion* _pConv;
	streampos  _requestedpos, _lastpos;  
	xmlTextReaderPtr _reader;
	xmlTextWriterPtr _writer;
	xmlOutputBufferPtr _buf;
	//	xmlBufferPtr _buf;
	bool _LookingForNamespace;
public:	
	bool _SkipNextRead;
};

//*************************************************
/// Abstract class containing common functionality for XML formats.
class XMLBaseFormat : public OBFormat
{
protected:
	XMLConversion* _pxmlConv;
	
	//formating for output
	string _prefix;
	int baseindent, ind;
	string nsdecl;
	int _embedlevel;

public:
	virtual const char* NamespaceURI()const=0;
	virtual bool DoElement(const std::string& ElName){return false;};
	virtual bool EndElement(const std::string& ElName){return false;};
	 /// The tag at the end of the chemical object e.g. "/molecule>"
	virtual const char* EndTag(){return ">";};
	
protected:
	void SetFormatting(OBConversion* pConv, std::ostream& ofs)
	{
		const char* nstxt = pConv->IsOption("N");
		if(nstxt)
		{
			_prefix = nstxt;
			_prefix += ":";
		}
		else
			_prefix.erase();

		const char* p = pConv->IsOption("i");
		baseindent = p ? atoi(p) : 0;

		ind=0;

		if(!pConv->IsOption("x"))
		{
			ofs << "<?xml version=\"1.0\"?>\n";
			nsdecl = " xmlns=\"";
			nsdecl += NamespaceURI(); 
			nsdecl += "\""; 
		}
		else
			nsdecl.erase();
	}

	xmlTextReaderPtr reader() const
	{
		return _pxmlConv->GetReader();
	}

	xmlTextWriterPtr writer() const
	{
		return _pxmlConv->GetWriter();
	}
	
	void OutputToStream()
	{
		_pxmlConv->OutputToStream();
	}
	
};

//*************************************************
///Abstract class containing common functionality for XML formats which represent molecules
class XMLMoleculeFormat : public XMLBaseFormat
{
protected:
	OBMol* _pmol;

public:
	virtual bool ReadChemObject(OBConversion* pConv)
	{
		std::string auditMsg = "OpenBabel::Read molecule ";
		std::string description(Description());
		auditMsg += description.substr(0,description.find('\n'));
		obErrorLog.ThrowError(__FUNCTION__, auditMsg, obAuditMsg);

		//With j option, reuse pmol except for the first mol
		static OBMol* pmol;
		if(!pConv->IsOption("j",OBConversion::GENOPTIONS) || pConv->IsFirstInput())
			pmol = new OBMol;
		
		bool ret=ReadMolecule(pmol,pConv);

		if(ret && pmol->NumAtoms() > 0) //Do transformation and return molecule
			pConv->AddChemObject(pmol->DoTransformations(pConv->GetOptions(OBConversion::GENOPTIONS)));
		else
			pConv->AddChemObject(NULL);

		return ret;
	};

	virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv)
	{
		_pmol = dynamic_cast<OBMol*>(pOb);
		if(!_pmol)
			return false;

		_pxmlConv = XMLConversion::GetDerived(pConv,true);
		if(!_pxmlConv)
			return false;
		_embedlevel = -1;
		return _pxmlConv->ReadXML(this,pOb);
	};


	virtual bool WriteChemObject(OBConversion* pConv)
	{
		//Retrieve the target OBMol
		OBBase* pOb = pConv->GetChemObject();
		OBMol* pmol = dynamic_cast<OBMol*> (pOb);
		bool ret=false;
		if(pmol)
		{	
			if(pmol->NumAtoms()==0)
			{
				std::string auditMsg = "OpenBabel::Molecule ";
				auditMsg += pmol->GetTitle();
				auditMsg += " has 0 atoms";
				obErrorLog.ThrowError(__FUNCTION__,
						auditMsg,
						obInfo);
			}
			ret=true;

			std::string auditMsg = "OpenBabel::Write molecule ";
			std::string description(Description());
			auditMsg += description.substr(0,description.find('\n'));
			obErrorLog.ThrowError(__FUNCTION__,
					      auditMsg,
					      obAuditMsg);

			if(!pConv->IsOption("j",OBConversion::GENOPTIONS) || pConv->IsLast()) //With j option, output only at end
			{
				ret=WriteMolecule(pmol,pConv);
				delete pOb;
			}
		}
		return ret;
	};

	const type_info& GetType()
	{
		return typeid(OBMol*);
	};

};


}//namespace

//! \file
//! \brief Declaration of XMLConversion, 
//!  declaration and definition of XMLBaseFormat and XMLMoleculeFormat 

/**********************************************************************
Copyright (C) 2005 by Chris Morley
 
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include "xml.h"

using namespace std;
namespace OpenBabel
{

//static variable
XMLBaseFormat* XMLConversion::_pDefault=NULL;

XMLConversion::XMLConversion(OBConversion* pConv)
		: OBConversion(*pConv), _reader(NULL), _writer(NULL),
		_LookingForNamespace(false),_SkipNextRead(false),
		  _lastpos(0), _requestedpos(0)
{
	_pConv = pConv;
	pConv->SetAuxConv(this);//marks original OBConversion object as having been extended 
	SetAuxConv(this);//marks this new object as extended (for use with OBConversion pointer)
}

bool XMLConversion::SetupReader()
{
	if(_reader)
		return true; //do not need to make a new reader

	//If the inputstream is not at the start (probably arising in fastsearch),
	//save its position and rewind so that the reader initialization is ok.
	//(Getting the requested object is handled in ReadXML(), when the format is known.) 
	_requestedpos = GetInStream()->tellg();
	if(_requestedpos)
		GetInStream()->seekg(0);

	//Set up a parser from an input stream 
	_reader = xmlReaderForIO(
			ReadStream, //xmlInputReadCallback (static member function)
			NULL,//xmlInputCloseCallback (static member function)
			this,       //context
			"",         //URL
			NULL,       //encoding
			0);         //options

	if (_reader == NULL)
	{
		cerr << "Cannot set up libxml2 reader" << endl;
		return false;
	}
	//A new reader immediately reads 4 bytes (presumably to determine
	//the encoding).
		_lastpos = GetInStream()->tellg();
	return true;
}

bool XMLConversion::SetupWriter()
{
  //Set up XML writer if one does not already exist
	if(_writer)
		return true;
  
	_buf = xmlOutputBufferCreateIO	(
			WriteStream, //xmlOutputWriteCallback 
			NULL,			   //xmlOutputCloseCallback
			this,        //context
			NULL);        //xmlCharEncodingHandlerPtr
	_writer = xmlNewTextWriter(_buf);

/*
	_buf = xmlBufferCreate();
  _writer = xmlNewTextWriterMemory(_buf, 0);
*/ 
	
	if(!_buf || !_writer)
	{
		cerr << "Error setting up xml writer\n" << endl;
    return false;
	}

	int ret = xmlTextWriterSetIndent(_writer,1);
	ret = xmlTextWriterSetIndentString(_writer, BAD_CAST " "); 
	return ret==0;
}

XMLConversion::~XMLConversion()
{
	if(_reader)
		xmlFreeTextReader(_reader);
//	if(_writer)
//		xmlFreeTextWriter(_writer); was crashing
		//xmlBufferFree(_buf);
}

///Called from each XML class during its construction
void XMLConversion::RegisterXMLFormat(XMLBaseFormat* pFormat, bool IsDefault, const char* uri)
{
	if(IsDefault || Namespaces().empty())
		_pDefault=pFormat;
	if(uri)
		Namespaces()[uri] = pFormat;
	else
		Namespaces()[pFormat->NamespaceURI()] = pFormat;
}

XMLConversion* XMLConversion::GetDerived(OBConversion* pConv, bool ForReading)
{
	XMLConversion* pxmlConv;
	if(!pConv->GetAuxConv())
		//Need to make an extended copy. It will be deleted by pConv's destructor
		pxmlConv =  new XMLConversion(pConv);
	else
	{
		//pConv has already had an extended copy made	
		pxmlConv = dynamic_cast<XMLConversion*>(pConv->GetAuxConv());
		if (!pxmlConv)
			return NULL;
	}

	if(ForReading)
	{
		pxmlConv->SetupReader();
		if(pConv->GetInStream()->tellg() < pxmlConv->_lastpos)
		{
			//Probably a new file; copy some member vars and renew the current reader
			pxmlConv->InFilename = pConv->GetInFilename();
			pxmlConv->pInFormat = pConv->GetInFormat();

			if(xmlReaderNewIO( pxmlConv->_reader, ReadStream, NULL, pxmlConv, "", NULL, 0)==-1)
				return false;
		}
	}
	else
		pxmlConv->SetupWriter();

	return pxmlConv;
}


bool XMLConversion::ReadXML(XMLBaseFormat* pFormat, OBBase* pOb)
{	
	if(_requestedpos)
	{
		//The initial stream position was not at the start, probably because of fastsearch
		//Read and discard the first object to synchronize the reader,
		//then continue getting the requested object.
		//Assumes the objects are all at the same level in the DOM tree.
		SetOneObjectOnly(); //probably already set
		streampos SavedReqestedPos = _requestedpos; 
		_requestedpos=0;//don't do this again
		ReadXML(pFormat,pOb);
		GetInStream()->seekg(SavedReqestedPos);
	}

	//**Parse
	int result=1;
	while(_SkipNextRead || (result=xmlTextReaderRead(_reader))==1) //read may not be called
	{
		_SkipNextRead=false;
		if(_LookingForNamespace)
		{
			const xmlChar* puri = xmlTextReaderConstNamespaceUri(_reader);
			if(puri)
			{
				string uri((const char*)puri);
				//Look up appropriate format class from the namespace URI
				NsMapType::iterator nsiter;
				nsiter = Namespaces().find(uri);
				if(nsiter!=Namespaces().end())
				{
					XMLBaseFormat* pNewFormat = nsiter->second;
					//Must have same target, e.g. OBMol, as current format 
					if(pNewFormat->GetType() == pFormat->GetType())
					{
						_LookingForNamespace=false;
						_SkipNextRead=true;
						SetInFormat(pNewFormat);
						pNewFormat->ReadMolecule(pOb,this);
						return true;
					}
				}
			}
		}

		const xmlChar* pname = xmlTextReaderConstLocalName(_reader);
		int typ = xmlTextReaderNodeType(_reader);
		if(typ==XML_READER_TYPE_SIGNIFICANT_WHITESPACE || !pname)
			continue; //Text nodes handled in format class
		string ElName((const char*)pname);

		//Pass the node on to the appropriate format class
		bool ret;
		if(typ==XML_READER_TYPE_ELEMENT)
			ret= pFormat->DoElement(ElName);
		else if(typ==XML_READER_TYPE_END_ELEMENT)
			ret= pFormat->EndElement(ElName);
		
		_lastpos = GetInStream()->tellg();

		if(!ret)
			//derived format callback has stopped processing by returning false;
			//leave reader intact so it can be continued to be used.
			return true;
	}

	if(result==-1)
	{
		cerr << "XML Parser failed in " << GetInFilename() << endl;
		GetInStream()->setstate(ios::eofbit);
	}
	return (result==0);// was result==0;
}

/////////////////////////////////////////////////////////
string XMLConversion::GetAttribute(const char* attrname)
{
	string AttributeValue;
	const xmlChar* pvalue  = xmlTextReaderGetAttribute(_reader, BAD_CAST attrname);
	if(pvalue)
		AttributeValue = (const char*)pvalue;
	return AttributeValue;
}

////////////////////////////////////////////////////////
string XMLConversion::GetContent()
{
	xmlTextReaderRead(_reader);
	const xmlChar* pvalue = xmlTextReaderConstValue(_reader);
	string value((const char*)pvalue);
	return value;
}

////////////////////////////////////////////////////////
bool XMLConversion::GetContentInt(int& value)
{
	xmlTextReaderRead(_reader);
	const xmlChar* pvalue = xmlTextReaderConstValue(_reader);
	if(!pvalue)
		return false;
	value = atoi((const char*)pvalue);
	return true;
}

////////////////////////////////////////////////////////
bool XMLConversion::GetContentDouble(double& value)
{
	xmlTextReaderRead(_reader);
	const xmlChar* pvalue = xmlTextReaderConstValue(_reader);
	if(!pvalue)
		return false;
	value = strtod((const char*)pvalue,NULL);
	return true;
}

//**********************************************
/// Utility function to read an input stream until a specified string is found 
streamsize gettomatch(istream& is, char* buf, streamsize count, const char* match) 
{
	//Reads chars from input stream into a buffer until either: 
	//  count chars have been read or
	//  the string match has been input.
	//The buffer is NOT terminated by a '\0' char.
	//The number of characters stored in buf is returned.  

	int matchlength = 0;
	char lastchar = EOF; //value if no vaild match provided
	if(match)
	{
		matchlength = strlen(match);
		lastchar = match[matchlength-1];
	}
	char* p = buf;
	streambuf* prb = is.rdbuf();
	int i;
	for(i=0;i<count;++i)
	{
		*p = prb->sbumpc();
		if(*p==EOF)
			break;
		if(*p++==lastchar)
		{
			const char* mptr = match + matchlength-2; //last char is already matched
			const char* bptr = p-2;
			while((*mptr-- == *bptr--) && (mptr >= match));

			if(mptr<match)
			{
				i++;
				break;//have found match
			}
		}
	}
	return i;
};
//***********************************************

///Static callback function for xmlReaderForIO()
int XMLConversion::ReadStream(void * context, char * buffer, int len)
{
	//Reads up to the next '>'
	XMLConversion* pConv = static_cast<XMLConversion*>(context);
	istream* ifs = pConv->GetInStream();
	if(ifs->eof())
		return 0;
	const char* endtag = NULL;
	OBFormat* pFormat = pConv->GetInFormat();
	XMLBaseFormat* pxmlFormat = static_cast<XMLBaseFormat*>(pFormat);
	if(pxmlFormat)
		endtag = pxmlFormat->EndTag();

	static char* OrigBuffer;
	if(len==4)
		OrigBuffer = buffer;

	return gettomatch(*ifs, buffer, len , endtag);//was + OrigBuffer - buffer
}

int XMLConversion::WriteStream(void * context, const char * buffer, int len)
{
	XMLConversion* pxmlConv = static_cast<XMLConversion*>(context);
	ostream* ofs = pxmlConv->GetOutStream();
	ofs->write(buffer,len);
	if(!ofs)
		return -1;
	ofs->flush();
	return len;
}

} //namespace OpenBabel
// http://xmlsoft.org/html/libxml-xmlreader.html



/*
Namespaces
This class keeps a map of xml namespaces and the derived classes which implement them.
It is populated on start by the derived classes calling RegisterXMLFormat
from their default constructors.

When ReadChemObject()of a derived class is called, the current format is set by
there and is used for all namespaces. So if CMLFormat is called it will find
all the molecules in a CMLReact file.

When ReadChemObject() of the base class XMLFormat is the one called (e.g for files
with extension .xml), the initial current format is the default format(see below)
providing it handles the same chemical object type as the output format. The 
ReadChemObject() of the default format is called and processing is as if the
default format was called, except that the first explicit namespace declaration
in the xml file that appears in the map can switch the handling to its associated
format.
 
The default format is either the first class to register or one which identifies itself as
the default. */


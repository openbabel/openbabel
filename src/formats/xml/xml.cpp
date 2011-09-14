/**********************************************************************
Copyright (C) 2005-2006 by Chris Morley
Some portions Copyright (C) 2006 by Geoffrey R. Hutchison

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/
#include <openbabel/babelconfig.h>
#include <openbabel/xml.h>

using namespace std;
namespace OpenBabel
{

  //static variable
  XMLBaseFormat* XMLConversion::_pDefault=NULL;

  XMLConversion::XMLConversion(OBConversion* pConv)
    : OBConversion(*pConv),
      _requestedpos(0), _lastpos(0),
      _reader(NULL), _writer(NULL),
      _LookingForNamespace(false), _SkipNextRead(false)
  {
    pLineEndBuf=NULL;
    _pConv = pConv;
    pConv->SetAuxConv(this);//marks original OBConversion object as having been extended
    SetAuxConv(this);//marks this new object as extended (for use with OBConversion pointer)
  }

  bool XMLConversion::SetupReader()
  {
    if(_reader)
      return true; //do not need to make a new reader

	//setup libxml2 for use in a potentially multithreaded
	//environment
	xmlInitParser();

    //If the inputstream is not at the start (probably arising in fastsearch),
    //save its position and rewind so that the reader initialization is ok.
    //(Getting the requested object is handled in ReadXML(), when the format is known.)
    _requestedpos = GetInStream()->tellg();
    if(_requestedpos < 0)
      _requestedpos = 0;
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

    _buf = xmlOutputBufferCreateIO  (
                                     WriteStream, //xmlOutputWriteCallback
                                     NULL,         //xmlOutputCloseCallback
                                     this,        //context
                                     NULL);        //xmlCharEncodingHandlerPtr
    _writer = xmlNewTextWriter(_buf);

    if(!_buf || !_writer)
      {
        cerr << "Error setting up xml writer\n" << endl;
        return false;
      }

    int ret;
    if(IsOption("c"))
      ret = xmlTextWriterSetIndent(_writer,0);
    else
      {
        ret = xmlTextWriterSetIndent(_writer,1);
        ret = xmlTextWriterSetIndentString(_writer, BAD_CAST " ");
      }
    return ret==0;
  }

  XMLConversion::~XMLConversion()
  {
    if(_reader) {
      xmlFreeTextReader(_reader);
      _reader = NULL;
    }
    if(_writer) {
//      xmlTextWriterEndDocument(_writer); //if hasn't been called ealier
        xmlFreeTextWriter(_writer);// was crashing
        _writer = NULL;
    }
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

  ///Returns the extended form of the OBConversion object with an xml reader or writer,
  /// if this has not already been done.
  XMLConversion* XMLConversion::GetDerived(OBConversion* pConv, bool ForReading)
  {
    XMLConversion* pxmlConv;
    if(!pConv->GetAuxConv())
      //Need to make an extended copy. It will be deleted by pConv's destructor
      pxmlConv =  new XMLConversion(pConv);
    else
      {
        //pConv has already had an extended copy made
        *pConv->GetAuxConv() = *pConv; //ensure they have the same OBConversion data
        pxmlConv = dynamic_cast<XMLConversion*>(pConv->GetAuxConv());
        if (!pxmlConv)
          return NULL;
      }

    if(ForReading)
      {
        streampos pos = pConv->GetInStream()->tellg();
        if(pos < pxmlConv->_lastpos || pxmlConv->_lastpos<0)
          {
            //Probably a new file; copy some member vars and renew the current reader
            xmlFreeTextReader(pxmlConv->_reader); //need a new reader to read files with <?xml?>
            pxmlConv->_reader = NULL;
            pxmlConv->InFilename = pConv->GetInFilename();
            pxmlConv->pInFormat = pConv->GetInFormat();

          }
        pxmlConv->SetupReader();
      }
    else
    {
      pxmlConv->SetupWriter();
      pxmlConv->SetLast(pConv->IsLast()); //Copy IsLast flag to the extended object
    }
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
    while(!GetInStream()->bad() && (_SkipNextRead || (result=xmlTextReaderRead(_reader))==1)) //read may not be called
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
      else
        continue;
      _lastpos = GetInStream()->tellg();

      if(!ret)
        //derived format callback has stopped processing by returning false;
        //leave reader intact so it can be continued to be used.
        if(!IsOption("n",OBConversion::INOPTIONS))
          {
            _LookingForNamespace = true;
            return true;
          }
      }

    if(result==-1)
    {
      xmlError* perr = xmlGetLastError();
      if(perr && perr->level!=XML_ERR_NONE)
        {
          obErrorLog.ThrowError("XML Parser " + GetInFilename(),
                                perr->message, obError);
        }
      xmlResetError(perr);
      GetInStream()->setstate(ios::eofbit);
      return false;
    }
    return GetInStream()->good() && result!=0;
  }

  /////////////////////////////////////////////////////////
  ///Read and discard XML text up to the next occurrence of the tag e.g."/molecule>"
  ///This is left as the current node. Returns 1 on success, 0 if not found, -1 if failed.
  int XMLConversion::SkipXML(const char* ctag)
  {
    string tag(ctag);
    tag.erase(--tag.end()); //remove >
    int targettyp = XML_READER_TYPE_ELEMENT;
    if(tag[0]=='/')
      {
        tag.erase(0,1);
        targettyp = XML_READER_TYPE_END_ELEMENT;
      }

    int result;
    while((result = xmlTextReaderRead(_reader))==1)
      {
        if(xmlTextReaderNodeType(_reader)==targettyp
           && !xmlStrcmp(xmlTextReaderConstLocalName(_reader), BAD_CAST	tag.c_str()))
          break;
      }
    return result;
  }
  /////////////////////////////////////////////////////////
  string XMLConversion::GetAttribute(const char* attrname)
  {
    string AttributeValue;
    xmlChar* pvalue  = xmlTextReaderGetAttribute(_reader, BAD_CAST attrname);
    if(pvalue)
    {
      AttributeValue = (const char*)pvalue;
      xmlFree(pvalue);
    }
    return AttributeValue;
  }

  ////////////////////////////////////////////////////////
  string XMLConversion::GetContent()
  {
    xmlTextReaderRead(_reader);
    const xmlChar* pvalue = xmlTextReaderConstValue(_reader);
    string value((const char*)pvalue);
    return Trim(value);
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

  ////////////////////////////////////////////////////////
  ///Static callback function for xmlReaderForIO(). Reads up to the next '>', or len chars.

  int XMLConversion::ReadStream(void * context, char * buffer, int len)
  {
    //@todo worry about non-ascii coding
    XMLConversion* pConv = static_cast<XMLConversion*>(context);
    istream* ifs = pConv->GetInStream();
    if(!ifs->good() || ifs->eof())
      return 0;

    ifs->get(buffer, len+1, '>');
    streamsize count = strlen(buffer);

    if(ifs->peek()=='>')
      {
        ifs->ignore();
        buffer[count] = '>';
        buffer[++count] = '\0';
      }

		if (ifs->peek() == '\n' || ifs->peek() == '\r')
			{
				ifs->get(); // remove any trailing endlines
			}
    return count;
  }

  //////////////////////////////////////////////////////////
  int XMLConversion::WriteStream(void * context, const char * buffer, int len)
  {
    XMLConversion* pxmlConv = static_cast<XMLConversion*>(context);
    ostream* ofs = pxmlConv->GetOutStream();
    if(len>0)                //a call with len=0 coming from xmlFreeTextWriter
    {                        //called from destructor of XMLConversion was causing crash
      ofs->write(buffer,len);
      if(!ofs)
        return -1;
      ofs->flush();
    }
    return len;
  }

} //namespace OpenBabel
// http://xmlsoft.org/html/libxml-xmlreader.html

/*
Programming notes on XML formats

So that there would be no limitation of file sizes, the libxml2
reader was chosen. Rather than build a whole xml tree internally
as a DOM parser does, this provides callbacks when each element,
etc. is encountered (like SAX). Nevertheless it is aware of the
XML structure and will fail if it enounters irregular input. It
is therefore necessary to use a single instance of the reader for
each conversion process, rather than one for each object as would
have been more natural in OB (see below). This input process can
span multiple input files and is associated with the OBConversion
object - in particular the reader object is destroyed at the same
time as the OBConversion object. But it is not as simple as using
an extended OBConversion derived from the base class, because the
base OBConversion object has been constructed before the XML format
has been called. It might have been possible to have the reader as a
member variable in OBConversion, but that would make an undesirable
dependency for obconversion.cpp on the XML formats.

The way it has been done maintains generality and no dependency.
OBConversion is given a member variable pAuxConv which is a pointer
to an OBConversion object. This is deleted in the OBConversion
destructor. By default pAuxConv is NULL.
XMLConversion is a class derived from OBConversion and
contains the interfacing with libxml2 for both reading and writing.
When a conversion involves an XML format, an instance of it is made
and pAuxConv in the original OBConversion is set to point to it.
This process is potentially extendable to allow other, as yet
unwritten, OBConversion extensions by having a chain of pointers
to derived OBConversion objects through their pAuxConv members,
with the last one being NULL.

The design has to make sure that multi-object files are handled
in a way consistent with the rest of OpenBabel. This is based on
formats such as SMILES and MDL mol where the objects are just
concatenated. OpenBabel converts by reading one object at a time
from the input stream. The position in the input stream (obtained
from tellg) is also used to skip objects and as the index in fast
searching. These depend on input file position being left between
objects ready to read the next one.

This causes some difficulty when using libxml2 as the XML parser
because it is a C application and does not have C++ input streams.
xmlReaderForIO is used which requests input data from the callback
routine XMLConversion::ReadStream(). This inputs chunks of characters
up to and including '>'. This ensures that the input stream is between
objects after an object has been parsed, ready for the next one.

Parsing XML
At the start and end of each element the DoElement() and EndElement()
routines respectively in the format are called. The name of the element
is passed as a parameter and up to now it has been considered sufficient
to find the appropriate code using a set of if else statements. Only
those of interest need be handled. The attributes and content of the
element are found by calling libxml2 routines from within the format class.
Parsing is stopped and an object returned to OBConversion when false is
returned from DoElement or (more usually) from EndElement.

Namespaces
XMLConversion class keeps a static map of xml namespaces and the classes
derived from XMLBaseFormat which implement them. It is populated on startup
by the format classes calling RegisterXMLFormat from their default constructors.

When ReadChemObject() of a format class is called, the current format is set by
there and is used for all namespaces. So if CMLFormat is called it will find
all the molecules in a CMLReact file.

When ReadChemObject() of the base class XMLFormat is the one called (e.g for files
with extension .xml), the initial current format is the default format(see below)
providing it handles the same chemical object type as the output format. The
ReadChemObject() of the default format is called and processing is as if the
default format was called, except that the first explicit namespace declaration
in the xml file that appears in the map can switch the handling to its associated
format.

The default format is either the first class to register or one which identifies
itself as the default when calling RegisterXMLFormat().

*/


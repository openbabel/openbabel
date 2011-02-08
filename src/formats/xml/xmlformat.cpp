/**********************************************************************
Definition of XMLFormat
Copyright (C) 2005 by Chris Morley

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

class XMLFormat : XMLBaseFormat
{
public:
  XMLFormat()
  {
      OBConversion::RegisterFormat("xml",this);
  }

  const char* Description()
  {
    return
      "General XML format\n"
      "Calls a particular XML format depending on the XML namespace.\n\n"
      "This is a general XML \"format\" which reads a generic XML file and infers\n"
      "its format from the namespace as given in a xmlns attribute on an element.\n"
      "If a namespace is recognised as associated with one of the XML formats in\n"
      "Open Babel, and the type of the object (e.g. a molecule) is appropriate to\n"
      "the output format then this is used to input a single object. If no namespace\n"
      "declaration is found the default format (currently CML) is used.\n\n"

      "The process is repeated for any subsequent input so that it is possible to\n"
      "input objects written in several different schemas from the same document.\n"
      "The file :file:`CMLandPubChem.xml` illustrates this and contains molecules in\n"
      "both CML and PubChem formats.\n\n"

      "This implementation uses libxml2.\n\n"
      "Read Options, e.g. -an\n"
      " n  Read objects of first namespace only\n\n";
  }

  const char* NamespaceURI()const{return "Undefined";};

  unsigned Flags()
  {
    return READXML|NOTWRITABLE;
  }

  bool ReadChemObject(OBConversion* pConv)
  {
    XMLBaseFormat* pDefault = XMLConversion::GetDefaultXMLClass();
    if(!pDefault || pDefault==this)
    {
      obErrorLog.ThrowError("XML Format", "There is no acceptable default XML Format", obError);
      return false;
    }
    if(pConv->GetOutFormat()->GetType() == pDefault->GetType())
    {
      //Extend the OBConversion (to include XML parsing).
      XMLConversion* pxmlConv = XMLConversion::GetDerived(pConv);
      pxmlConv->LookForNamespace();
      return pDefault->ReadChemObject(pConv);
    }
    else
    {
      //Chemical object type handled by the output format is not
      //the same as that handled by the default format.
      return false;
    }
  };

  bool ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {
    XMLBaseFormat* pDefault = XMLConversion::GetDefaultXMLClass();
    if(pConv->GetOutFormat()->GetType() == pDefault->GetType())
    {
      XMLConversion* pxmlConv = XMLConversion::GetDerived(pConv);
      pxmlConv->LookForNamespace();
      pxmlConv->AddOption("m",OBConversion::INOPTIONS);
      return pDefault->ReadMolecule(pOb, pConv);
    }
    else
    {
      //Chemical object type handled by the output format is not
      //the same as that handled by the default format.
      obErrorLog.ThrowError("XML Format", "Need to specify the input XML format more precisely", obError);
      return false;
    }
  };

};

//Make an instance of the format class
XMLFormat theXMLFormat;

} //namespace OpenBabel

//! \file
//! \brief Definition of XMLFormat,

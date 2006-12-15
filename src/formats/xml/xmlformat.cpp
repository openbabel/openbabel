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
      return " \
General XML format\n \
Calls a particular XML format depending on the XML namespace,\n \
or the default format (which is probably CML).\n \
This implementation uses libxml2.\n \
Read option, e.g. -an\n \
n  Read objects of first namespace only\n \
\n";
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

/**********************************************************************
Copyright (C) 2006 by Geoff Hutchison
 
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

#include <openbabel/babelconfig.h>
#include <openbabel/xml.h>

#ifdef WIN32
#pragma warning (disable : 4800)
#endif

using namespace std;
namespace OpenBabel
{

class ChemDrawXMLFormat : public XMLMoleculeFormat
{
public:
	ChemDrawXMLFormat() 
	{
		OBConversion::RegisterFormat("cdxml", this);
		XMLConversion::RegisterXMLFormat(this);
	}
	virtual const char* NamespaceURI()const{return "http://www.cambridgesoft.com/xml/cdxml.dtd";}
  virtual const char* Description()
  {
      return " \
ChemDraw CDXML format \n \
Minimal extraction of chemical structure information only.\n \
\n";
};

  virtual const char* SpecificationURL()
  {return "http://www.cambridgesoft.com/services/documentation/sdk/chemdraw/cdx/";}


  virtual unsigned int Flags()
  {
      return NOTWRITABLE;
  };

	virtual bool DoElement(const string& name);
	virtual bool EndElement(const string& name);

	// EndTag is used so that the stream buffer is is filled with the XML from
	// complete objects, as far as possible. 
	virtual const char* EndTag(){ return "/fragment"; };

private:
  OBAtom _tempAtom; //!< A temporary atom as the atom tag is read
  OBBond _tempBond; //!< A temporary bond as the bond tag is read
};

////////////////////////////////////////////////////////////////////

ChemDrawXMLFormat theChemDrawXMLFormat;

////////////////////////////////////////////////////////////////////

bool ChemDrawXMLFormat::DoElement(const string& name)
{
	if(name=="fragment") 
	{
		//This is the start of the molecule we are extracting and it will
		//be put into the OBMol* _pmol declared in the parent class.
		//initialise everything
    _tempAtom.Clear();
    _tempBond.Clear();

    _pmol->SetDimension(2);
		_pmol->BeginModify();
	}
	else if(name=="n")
	{
    _tempAtom.SetAtomicNum(6); // default is carbon
    if (_pxmlConv->GetAttribute("Element"))
      _tempAtom.SetAtomicNum(atoi(_pxmlConv->GetAttribute("Element")));

    string coords = _pxmlConv->GetAttribute("p");
	}
	else if(name=="b")
	{

    if (_pxmlConv->GetAttribute("Order"))
      _tempBond.SetBO(atoi(_pxmlConv->GetAttribute("Order")));
    _pxmlConv->GetAttribute("B");
    _pxmlConv->GetAttribute("E");

	}

	return true;
}

bool ChemDrawXMLFormat::EndElement(const string& name)
{
  unsigned int i;
	if(name=="n")
	{
	}
	else if(name=="b")
	{
	}
	else if(name=="fragment") //this is the end of the molecule we are extracting
	{
		_pmol->EndModify();
		return false;//means stop parsing
	}
	return true;
}	


}//namespace

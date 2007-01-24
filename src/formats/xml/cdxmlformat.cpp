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
		OBConversion::RegisterFormat("cdxml", this, "chemical/x-cdxml");
		XMLConversion::RegisterXMLFormat(this, false, "http://www.camsoft.com/xml/cdxml.dtd");
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

  virtual const char* GetMIMEType() 
  { return "chemical/x-cdxml"; };

  virtual const char* SpecificationURL()
  {return "http://www.cambridgesoft.com/services/documentation/sdk/chemdraw/cdx/";}


  virtual unsigned int Flags()
  {
      return READXML | NOTWRITABLE;
  };

	virtual bool DoElement(const string& name);
	virtual bool EndElement(const string& name);

	// EndTag is used so that the stream buffer is is filled with the XML from
	// complete objects, as far as possible. 
	virtual const char* EndTag(){puts("end tag"); return "/fragment>"; };

    //atoms and bonds might have no content, so EndElement is not always called
    // that's why we need to ensure that atoms and bonds are really added.
    void EnsureEndElement(void);

private:
  OBAtom _tempAtom; //!< A temporary atom as the atom tag is read
  OBBond _tempBond; //!< A temporary bond as the bond tag is read
  map<int, int>atoms; //! maps chemdraw atom id to openbabel idx.
};

////////////////////////////////////////////////////////////////////

ChemDrawXMLFormat theChemDrawXMLFormat;

////////////////////////////////////////////////////////////////////

bool ChemDrawXMLFormat::DoElement(const string& name)
{
  string buf;
  if(name=="fragment") 
  {
    //This is the start of the molecule we are extracting and it will
    //be put into the OBMol* _pmol declared in the parent class.
    //initialise everything
    _tempAtom.Clear();
    _tempBond.Clear();
    atoms.clear();

    _pmol->SetDimension(2);
    _pmol->BeginModify();
  }
  else if(name=="n")
  {
    EnsureEndElement();
    _tempAtom.SetAtomicNum(6); // default is carbon
     buf = _pxmlConv->GetAttribute("id");
    if (buf.length())
      _tempAtom.SetIdx(atoi(buf.c_str()));
   buf = _pxmlConv->GetAttribute("Element");
    if (buf.length())
      _tempAtom.SetAtomicNum(atoi(buf.c_str()));

    buf = _pxmlConv->GetAttribute("p"); // coords
    if (buf.length())
    {
      double x = 0., y = 0.;
      sscanf(buf.c_str(), "%lf %lf", &x, &y);
      _tempAtom.SetVector(x, y, 0.);
    }
  }
  else if(name=="b")
  {
    EnsureEndElement();
    _tempBond.SetBO(1); //default value
    buf = _pxmlConv->GetAttribute("Order");
    if (buf.length())
      _tempBond.SetBO(atoi(buf.c_str()));
    buf = _pxmlConv->GetAttribute("B");
    if (buf.length())
      _tempBond.SetBegin(_pmol->GetAtom(atoms[atoi(buf.c_str())]));
    buf = _pxmlConv->GetAttribute("E");
    if (buf.length())
      _tempBond.SetEnd(_pmol->GetAtom(atoms[atoi(buf.c_str())]));
  }

  return true;
}

bool ChemDrawXMLFormat::EndElement(const string& name)
{
  unsigned int i;
  if(name=="n")
  {
    _pmol->AddAtom(_tempAtom);
    atoms[_tempAtom.GetIdx()] = _pmol->NumAtoms();
    _tempAtom.Clear();
  }
  else if(name=="b")
  {
    _pmol->AddBond(_tempBond);
    _tempBond.Clear();
    _tempBond.SetBO(0);
  }
  else if(name=="fragment") //this is the end of the molecule we are extracting
  {
    EnsureEndElement();
    _pmol->EndModify();
    atoms.clear();
    return false;//means stop parsing
  }
  return true;
}	

void ChemDrawXMLFormat::EnsureEndElement(void)
{
  if (_tempAtom.GetAtomicNum() != 0)
  {
    _pmol->AddAtom(_tempAtom);
    atoms[_tempAtom.GetIdx()] = _pmol->NumAtoms();
    _tempAtom.Clear();
  }
  else if (_tempBond.GetBO() != 0)
  {
    _pmol->AddBond(_tempBond);
    _tempBond.Clear();
    _tempBond.SetBO(0);
  }
}


}//namespace

/**********************************************************************
Copyright (C) 2008 by Chris Morley

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/
#include <openbabel/babelconfig.h>
#include <openbabel/obconversion.h>
#include <openbabel/text.h>

using namespace std;
namespace OpenBabel
{

class TextFormat : public OBFormat
{
public:
  TextFormat()
  {
    OBConversion::RegisterFormat("text",this);
  }

  virtual const char* Description() //required
  {
    return
     "Read and write raw text\n"
     "Facilitates the input of boilerplate text with babel commandline" ;
  }

/////////////////////////////////////////////////////////////////
  virtual bool ReadChemObject(OBConversion* pConv)
  {
    //Makes a new OBText
    OBText* pReact = new OBText;
    bool ret=ReadMolecule(pReact,pConv); //call the "API" read function

    std::string auditMsg = "OpenBabel::Read text ";
    std::string description(Description());
    auditMsg += description.substr(0,description.find('\n'));
    obErrorLog.ThrowError(__FUNCTION__,
              auditMsg,
              obAuditMsg);

    if(ret) //Do transformation and return molecule
      return pConv->AddChemObject(pReact->DoTransformations(pConv->GetOptions(OBConversion::GENOPTIONS),pConv))!=0;
    else
        pConv->AddChemObject(nullptr);
    return false;
  }

///////////////////////////////////////////////////////////////////////
  virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {
    //It's really text, not a molecule.
    OBText* pText = dynamic_cast<OBText*>(pOb);
    if (!pText)
      return false;
    string fileText(istreambuf_iterator<char>(*pConv->GetInStream()), istreambuf_iterator<char>());
    pText->SetText(fileText);
    return !fileText.empty();
  }

  virtual bool WriteChemObject(OBConversion* pConv)
  {
    //Output an OBText object and delete any other type.
    OBBase* pOb = pConv->GetChemObject();
    OBText* pText = dynamic_cast<OBText*>(pOb);
    if(!pText)
    {
      delete pOb;
      return false;
    }
    else
    {
      ostream* ofs = pConv->GetOutStream();
      if(ofs)
        *ofs << pText->GetText();
      return (bool)*ofs;
    }
  }
};

  //*********************************************************************
//Make an instance of the format class
TextFormat theTextFormat;


} //namespace OpenBabel


/**********************************************************************
Copyright (C) 2004 by Chris Morley

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
#include "openbabel/babelconfig.h"
#include <string>
#include <iomanip>
#include <typeinfo>
#include "openbabel/mol.h"
#include "openbabel/obconversion.h"
#include "openbabel/reaction.h"

using namespace std;
//using std::tr1::shared_ptr;

namespace OpenBabel
{
class RXNFormat : public OBFormat
{
public:
  //Register this format type ID
  RXNFormat()
  {
      OBConversion::RegisterFormat("rxn",this);
  }

  virtual const char* Description()
  {
      return
        "MDL RXN format\n"
        "The MDL reaction format is used to store information on chemical reactions.\n"
        "Output Options, e.g. -xA\n"
        " A  output in Alias form, e.g. Ph, if present\n\n";
  };

  virtual const char* GetMIMEType()
  { return "chemical/x-mdl-rxn"; };

  virtual const char* TargetClassDescription()
  {
      return OBReaction::ClassDescription();
  };

  const type_info& GetType()
  {
    return typeid(OBReaction*);
  };


  ////////////////////////////////////////////////////
  /// The "API" interface functions
  virtual bool ReadMolecule(OBBase* pReact, OBConversion* pConv);
  virtual bool WriteMolecule(OBBase* pReact, OBConversion* pConv);

  ////////////////////////////////////////////////////
  /// The "Convert" interface functions
  virtual bool ReadChemObject(OBConversion* pConv)
  {
    //Makes a new OBReaction and new associated OBMols
    OBReaction* pReact = new OBReaction;
    bool ret=ReadMolecule(pReact,pConv); //call the "API" read function

    std::string auditMsg = "OpenBabel::Read reaction ";
    std::string description(Description());
    auditMsg += description.substr(0,description.find('\n'));
    obErrorLog.ThrowError(__FUNCTION__,
              auditMsg,
              obAuditMsg);

    if(ret) //Do transformation and return molecule
      return pConv->AddChemObject(pReact->DoTransformations(pConv->GetOptions(OBConversion::GENOPTIONS),pConv))!=0;
    else
    {
      pConv->AddChemObject(NULL);
      delete pReact;
      pReact=NULL;
      return false;
    }
};

  virtual bool WriteChemObject(OBConversion* pConv)
  {
    //WriteChemObject() always deletes the object retrieved by GetChemObject
    //For RXN NO LONGER deletes the associated molecules which are handled by a smart pointer
    //Cast to the class type need, e.g. OBMol
    OBBase* pOb=pConv->GetChemObject();
    OBReaction* pReact = dynamic_cast<OBReaction*>(pOb);
    if(pReact==NULL)
        return false;

    bool ret=false;
    ret=WriteMolecule(pReact,pConv);

    std::string auditMsg = "OpenBabel::Write reaction ";
    std::string description(Description());
          auditMsg += description.substr( 0, description.find('\n') );
          obErrorLog.ThrowError(__FUNCTION__,
                                auditMsg,
                                obAuditMsg);
    delete pOb;
    return ret;
  };
};

//Make an instance of the format class
RXNFormat theRXNFormat;

/////////////////////////////////////////////////////////////////
bool RXNFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
{
    //It's really a reaction, not a molecule.
    //Doesn't make a new OBReactionObject, but does make mew reactant and product OBMols
   OBReaction* pReact = pOb->CastAndClear<OBReaction>();

    OBFormat* pMolFormat = pConv->FindFormat("MOL");
    if(pMolFormat==NULL || !pReact)
        return false;

    //	OBConversion MolConv(*pConv); //new copy to use to read associated MOL

    istream &ifs = *pConv->GetInStream();

    string ln;
    if (!getline(ifs,ln))
      return(false);
    if(Trim(ln).find("$RXN")!=0)
      return false; //Has to start with $RXN
    if (!getline(ifs,ln))
      return(false); //reaction title
    pReact->SetTitle(Trim(ln));

    if (!getline(ifs,ln))
      return false; //creator
    if (!getline(ifs, ln))
      return(false); //comment
    pReact->SetComment(Trim(ln));

    int nReactants, nProducts, i;
    ifs >> setw(3) >> nReactants >> setw(3) >> nProducts >> ws;
    if(!ifs) return false;

    if(nReactants + nProducts)
    {
      //Read the first $MOL. The others are read at the end of the previous MOL
      if(!getline(ifs, ln))
        return false;
      if(Trim(ln).find("$MOL")==string::npos)
        return false;
    }

    OBMol* pmol;

    for(i=0;i<nReactants;i++)
    {
      //Read a MOL file	using the same OBConversion object but with a different format
      pmol=new OBMol;
      if(!pMolFormat->ReadMolecule(pmol,pConv))
        obErrorLog.ThrowError(__FUNCTION__, "Failed to read a reactant", obWarning);
      else
      {
        shared_ptr<OBMol> p(pmol);
        pReact->AddReactant(p);
      }
    }

    for(i=0;i<nProducts;i++)
    {
      //Read a MOL file
      pmol=new OBMol;
      if(!pMolFormat->ReadMolecule(pmol,pConv))
        obErrorLog.ThrowError(__FUNCTION__, "Failed to read a product", obWarning);
      else
      {
        //        pReact->products.push_back(pmol);
        shared_ptr<OBMol> p(pmol);
        pReact->AddProduct(p);
      }
    }

    return(true);
}

/////////////////////////////////////////////////////////////////
bool RXNFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
    //It's really a reaction, not a molecule.
    //Cast output object to the class type need, i.e. OBReaction
    OBReaction* pReact = dynamic_cast<OBReaction*>(pOb);
    if(pReact==NULL)
        return false;

    OBConversion MolConv(*pConv); //new copy to use to write associated MOL
    MolConv.AddOption("no$$$$",OBConversion::OUTOPTIONS);
    MolConv.SetAuxConv(NULL); //temporary until a proper OBConversion copy constructor written

    OBFormat* pMolFormat = pConv->FindFormat("MOL");
    if(pMolFormat==NULL)
    {
      obErrorLog.ThrowError(__FUNCTION__, "MDL MOL format not available", obError);
        return false;
    }

    ostream &ofs = *pConv->GetOutStream();

    ofs << "$RXN" << endl;
    ofs << pReact->GetTitle() << endl;
    ofs << "  OpenBabel" << endl;
    ofs << pReact->GetComment() <<endl;

    ofs << setw(3) << pReact->NumReactants() << setw(3) << pReact->NumProducts() << endl;

    unsigned i;
    for(i=0;i<pReact->NumReactants();i++)
    {
      ofs << "$MOL" << endl;
      //Write reactant in MOL format
      pMolFormat->WriteMolecule(pReact->GetReactant(i).get(), &MolConv);
    }

    for(i=0;i<pReact->NumProducts();i++)
    {
      ofs << "$MOL" << endl;
      //Write reactant in MOL format
      pMolFormat->WriteMolecule(pReact->GetProduct(i).get(), &MolConv);
    }

    return true;
}

} //namespace

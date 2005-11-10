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
#ifdef WIN32
#pragma warning (disable : 4786)
#pragma warning (disable : 4251) //
#endif
#include "mol.h"
#include "obconversion.h"
#include "reaction.h"

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
          "MDL RXN format\n \
          \n";
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
			pConv->AddChemObject(pReact->DoTransformations(pConv->GetOptions(OBConversion::GENOPTIONS)));
		else
				pConv->AddChemObject(NULL);
		return ret;
	};

	virtual bool WriteChemObject(OBConversion* pConv)
	{
		//WriteChemObject() always deletes the object retrieved by GetChemObject
		//For RXN also deletes the associated molecules
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

		vector<OBMol*>::iterator itr;
		for(itr=pReact->reactants.begin();itr!=pReact->reactants.end();itr++)
				delete *itr;
		for(itr=pReact->products.begin();itr!=pReact->products.end();itr++)
				delete *itr;

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
    OBReaction* pReact = dynamic_cast<OBReaction*>(pOb);

    OBFormat* pMolFormat = pConv->FindFormat("MOL");
    if(pMolFormat==NULL)
        return false;

    //	OBConversion MolConv(*pConv); //new copy to use to read associated MOL

    istream &ifs = *pConv->GetInStream();
    char buffer[BUFF_SIZE];

    if (!ifs.getline(buffer,BUFF_SIZE))
        return false;
    if(strncmp(buffer,"$RXN",4))
        return false; //Has to start with $RXN
    if (!ifs.getline(buffer,BUFF_SIZE))
        return(false); //reactiontitle
    pReact->title = buffer;
		Trim(pReact->title);
    if (!ifs.getline(buffer,BUFF_SIZE))
        return(false); //creator
    if (!ifs.getline(buffer,BUFF_SIZE))
        return(false); //comment

    int nReactants, nProducts, i;
    if (!ifs.getline(buffer,BUFF_SIZE))
        return(false); //#reactants,products
    if(sscanf(buffer,"%3i%3i",&nReactants,&nProducts) != 2)
        return false;

    if(nReactants + nProducts)
    {
        //Read the first $MOL. The others are read at the end of the previous MOL
        if (!ifs.getline(buffer,BUFF_SIZE))
            return false;
        if(strncmp(buffer,"$MOL",4))
            return false;
    }

    OBMol* pmol;

    for(i=0;i<nReactants;i++)
    {
        //Read a MOL file	using the same OBConversion object but with a different format
        pmol=new OBMol;
        if(!pMolFormat->ReadMolecule(pmol,pConv))
	  obErrorLog.ThrowError(__FUNCTION__, "Failed to read a reactant", obWarning);
        pReact->reactants.push_back(pmol);
    }

    for(i=0;i<nProducts;i++)
    {
        //Read a MOL file
        pmol=new OBMol;
        if(!pMolFormat->ReadMolecule(pmol,pConv))
	  obErrorLog.ThrowError(__FUNCTION__, "Failed to read a product", obWarning);
        pReact->products.push_back(pmol);
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
		MolConv.SetAuxConv(NULL); //temporary until a proper OBConversion copy constructor written
	   
		OBFormat* pMolFormat = pConv->FindFormat("MOL");
    if(pMolFormat==NULL)
    {
      obErrorLog.ThrowError(__FUNCTION__, "MDL MOL format not available", obError);
        return false;
    }

    ostream &ofs = *pConv->GetOutStream();

    ofs << "$RXN" << endl;
    ofs << pReact->title.c_str() << endl;
    ofs << "  OpenBabel" << endl;
    ofs << "An experimental RXN file" <<endl;

    char buf[10];
    sprintf(buf,"%3u%3u",(unsigned)pReact->reactants.size(),pReact->products.size());
    ofs << buf << endl;

    vector<OBMol*>::iterator itr;
    for(itr=pReact->reactants.begin();itr!=pReact->reactants.end();itr++)
    {
        ofs << "$MOL" << endl;
        //Write reactant in MOL format
        pMolFormat->WriteMolecule(*itr, &MolConv); //does not delete associated molecules
    }

    for(itr=pReact->products.begin();itr!=pReact->products.end();itr++)
    {
        ofs << "$MOL" << endl;
        //Write product in MOL format
        pMolFormat->WriteMolecule(*itr, &MolConv); //does not delete associated molecules
    }
    return true;
}

}

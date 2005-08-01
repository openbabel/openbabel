/**********************************************************************
Copyright (C) 2004 by Chris Morley
 
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

using namespace std;
namespace OpenBabel
{

class ReactionSmilesFormat : public OBFormat
{
public:
    //Register this format type ID
    ReactionSmilesFormat()
    {
        OBConversion::RegisterFormat("rsmi",this);
    }

    virtual const char* Description() //required
    {
        return
            "ReactionSmiles format\n \
            No comments yet\n \
            ";
    };

    virtual const char* SpecificationURL(){return
            "http://www.daylight.com/dayhtml/doc/theory/theory.rxn.html";};

    //Flags() can return be any the following combined by | or be omitted if none apply
    // NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY
    virtual unsigned int Flags()
    {
        return 0;
    };

    //*** This section identical for most OBMol conversions ***
    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

    ////////////////////////////////////////////////////
    /// The "Convert" interface functions
    virtual bool ReadChemObject(OBConversion* pConv)
    {
        //Makes a new OBReaction and new associated OBMols
        OBReaction* pReact = new OBReaction;
        bool ret=ReadMolecule(pReact,pConv); //call the "API" read function

	std::string auditMsg = "OpenBabel::Read molecule ";
	std::string description(Description());
	auditMsg += description.substr(0,description.find('\n'));
	obErrorLog.ThrowError(__FUNCTION__,
			      auditMsg,
			      obAuditMsg);

        if(ret) //Do transformation and return molecule
            pConv->AddChemObject(pReact->DoTransformations(pConv->GetGeneralOptions()));
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

	std::string auditMsg = "OpenBabel::Write molecule ";
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
//***

//Make an instance of the format class
ReactionSmilesFormat theReactionSmilesFormat;

/////////////////////////////////////////////////////////////////
bool ReactionSmilesFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
{

  //It's really a reaction, not a molecule.
  //Doesn't make a new OBReactionObject, but does make mew reactant and product OBMols
  OBReaction* pReact = dynamic_cast<OBReaction*>(pOb);

  OBFormat* pSmiFormat = pConv->FindFormat("smi");
  if(pSmiFormat==NULL)
      return false;

  //Define some references so we can use the old parameter names
  istream &ifs = *pConv->GetInStream();
	
	string ln;
	getline(ifs,ln);
	if(!ifs) return false;
	vector<string> vec;
	tokenize(vec,ln,'>');

	//Read a MOL file	using the same OBConversion object but with a different format
  pmol=new OBMol;
  if(!pSmiFormat->ReadMolecule(pmol,pConv))
      cerr << "Failed to read the reactant" << endl;
  pReact->reactants.push_back(pmol);

	string reactantstring, agentstring, productstring;
	ifs >> reactantstring >> '>' >> agentstring >> '>' >> productstring;


}

////////////////////////////////////////////////////////////////

bool ReactionSmilesFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
        return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    //The old code goes here. Return true, or false on error
}

} //namespace OpenBabel

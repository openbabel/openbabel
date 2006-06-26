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
#include "reaction.h"

#ifdef WIN32
#pragma warning (disable : 4800)
#endif

using namespace std;
namespace OpenBabel
{

class CMLReactFormat : XMLBaseFormat
{
public:
	CMLReactFormat()
  {
      OBConversion::RegisterFormat("cmlr",this);
			XMLConversion::RegisterXMLFormat( this);
			OBConversion::RegisterOptionParam("l", this);
  }
	virtual const char* NamespaceURI()const
	{return "http://www.xml-cml.org/schema/cml2/react";} //guess

  const char* Description()
  {
      return " \
CML Reaction format\n \
Minimal implementation\n \
This implementation uses libxml2.\n \
Write options (e.g. -x1ac)\n \
1  output CML V1.0  or \n \
2  output CML V2.0 (default)\n \
a  output array format for atoms and bonds\n \
l  molecules in list\n \
h  use hydrogenCount for all hydrogens\n \
x  omit XML declaration\n \
N<prefix> add namespace prefix to elements\n \
\n";
	}

  unsigned Flags()
  {
    return 0;
  }
	virtual bool ReadChemObject(OBConversion* pConv);
	virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
  virtual bool WriteChemObject(OBConversion* pConv);
	virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);
	virtual bool DoElement(const string& ElName);
	virtual bool EndElement(const string& ElName);
	virtual const char* EndTag(){ return "/reaction>"; };

	const type_info& GetType()
	{
		return typeid(OBReaction*);
	};

private:
	string AddMolToList(vector<OBMol*>::iterator itr);

private:
	OBReaction* _preact;
	OBMol* pmol;
	map<string,OBMol*> Mols; //used on input
	map<string,OBMol> OMols; //used on output
	int nextmol;
  ostringstream ssout; //temporary output
};

//Make an instance of the format class
CMLReactFormat theCMLReactFormat;

////////////////////////////////////////////////////
bool CMLReactFormat::ReadChemObject(OBConversion* pConv)
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

	 //Do transformation and return reaction, if it has either reactants or products 
	if(ret && (pReact->reactants.size()!=0 || pReact->products.size()!=0)) //Do transformation and return molecule
		pConv->AddChemObject(pReact->DoTransformations(pConv->GetOptions(OBConversion::GENOPTIONS)));
	else
	{
		pConv->AddChemObject(NULL);
		return false;//don't continue after empty reaction
	}
	return ret;
}

bool CMLReactFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
{
	//Not a molecule: the object converted is a reaction
		_preact = dynamic_cast<OBReaction*>(pOb);
		if(!_preact)
			return false;

		_pxmlConv = XMLConversion::GetDerived(pConv);
		if(!_pxmlConv)
			return false;

		return _pxmlConv->ReadXML(this,pOb);
}


bool CMLReactFormat::DoElement(const string& name)
{
	if(name=="reaction")
	{
		pmol=NULL;
		_preact->title = _pxmlConv->GetAttribute("id");
	}

	else if(name=="molecule")
	{
		string reference = _pxmlConv->GetAttribute("ref");
		if(!reference.empty())
		{
			pmol = Mols[reference];
			if(!pmol)
			{
				cerr << " Molecule reference \"" << reference <<"\" not found" << endl;
				return false;
			}
		}
		else
		{
			pmol = new OBMol;	
			OBFormat* pCMLFormat = OBConversion::FindFormat("cml");
			if(!pCMLFormat)
				return false;
			_pxmlConv->_SkipNextRead=true;
			pCMLFormat->ReadMolecule(pmol, _pxmlConv);

			//Store all molecules in map
			string id  = pmol->GetTitle();
			Mols[id] = pmol;
			//TODO worry about deleting unused molecules
		}		
	}
	return true;
}

bool CMLReactFormat::EndElement(const string& name)
{
	if(name=="reactant")
	{
		if(!pmol)
			return false;
    _preact->reactants.push_back(pmol);
	}
	else if(name=="product")
	{
		if(!pmol)
			return false;
    _preact->products.push_back(pmol);
	}
	else if(name=="reaction")
	{
		return false;//means stop parsing
	}
	return true;
}

///////////////////////////////////////////////////////////////////
bool CMLReactFormat::WriteChemObject(OBConversion* pConv)
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
	obErrorLog.ThrowError(__FUNCTION__, auditMsg, obAuditMsg);

	/*Using the list form where several OBReaction objects can contain a pointer
	to an OBMol is a classic case for std::tr1::shared_ptr. Otherwise it
	is difficult to know when to delete the OBMol. Certainly not as here, which
	causes a crash with multiple reactions in list form.
	*/
	vector<OBMol*>::iterator itr;
	for(itr=pReact->reactants.begin();itr!=pReact->reactants.end();itr++)
			delete *itr;
	for(itr=pReact->products.begin();itr!=pReact->products.end();itr++)
			delete *itr;

	delete pOb;
	return ret;
}

////////////////////////////////////////////////////////////////
bool CMLReactFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
  //Badly named function: it's really a reaction, not a molecule.

	_pxmlConv = XMLConversion::GetDerived(pConv,false);
	if(!_pxmlConv)
		return false;

  //Cast output object to the class type need, i.e. OBReaction
  OBReaction* pReact = dynamic_cast<OBReaction*>(pOb);
  if(pReact==NULL)
      return false;

	//Two forms of output are supported:
	//With normal form, if more than one reaction to be output,
	//write <cml> at start and </cml> at end.
	//With list form, use ref attributes when refering to molecules and put
	//actual molecules in a separate list. 
	// <cml><moleculeList>...</moleculeList><reactionList>...</reactionList></cml>

	//The following element and attribute names may be output
	//For strict CMLReact compliance use the alternatives:
	static const xmlChar C_REACTIONLIST[] = "reactionList"; //"list"
	static const xmlChar C_MOLECULELIST[] = "moleculeList"; //"list"
	static const xmlChar C_LISTWRAPPER[]  = "mechanism";     //"cml"
	static const xmlChar C_WRAPPER[]      = "cml";
	static const xmlChar C_MOLECULE[]     = "molecule";
	static const xmlChar C_REACTION[]     = "reaction";
	static const xmlChar C_REACTANT[]     = "reactant";
	static const xmlChar C_PRODUCT[]      = "product";
	static const xmlChar C_REACTANTLIST[] = "reactantList";
	static const xmlChar C_PRODUCTLIST[]  = "productList";
	static const xmlChar C_REF[]          = "ref";
	static const xmlChar C_TITLE[]        = "title";
	
	bool list = _pxmlConv->IsOption("l"); //Output with molecules in a separate list

	ostringstream ssout;
	ostream* pOut = pConv->GetOutStream(); //the original output stream

	xmlChar* prefix = BAD_CAST _pxmlConv->IsOption("N");
	
	xmlChar* uri=NULL;
	
	_pxmlConv->AddOption("MolsNotStandalone",OBConversion::OUTOPTIONS); //inform CMLFormat 

//  OBConversion MolConv(*_pxmlConv); //new copy to use to write associated CML molecules
//	MolConv.SetAuxConv(NULL); //temporary until a proper OBConversion copy constructor written
//	MolConv.SetOneObjectOnly();
	 
	OBFormat* pCMLFormat = _pxmlConv->FindFormat("cml");
  if(pCMLFormat==NULL)
  {
      cerr << "CML format for molecules is not available\n" <<endl;
      return false;
  }
//	MolConv.AddOption("x", OBConversion::OUTOPTIONS); //no xml declaration
//	if(_pxmlConv->IsOption("N"))
//		MolConv.AddOption("N", OBConversion::OUTOPTIONS, _pxmlConv->IsOption("N")); //prefix

	//For first reaction
	if((pConv->GetOutputIndex()==1)) //OBConversion::Convert() is still using original pConv
	{
		if(!_pxmlConv->IsOption("x"))
		{
			xmlTextWriterStartDocument(writer(), NULL, NULL, NULL);
			uri=BAD_CAST NamespaceURI();
		}

		if(list)
		{
			xmlTextWriterStartElementNS(writer(), prefix, C_LISTWRAPPER, uri);
			ssout.clear();
			ssout.seekp(0);
			OMols.clear();
			nextmol=0;
			//With list form, use a temporary output for reactionList so that
			//moleculeList can be output first
			OutputToStream(); //flush what has already been written to pOut 
			_pxmlConv->SetOutStream(&ssout);
		}
		else if(!_pxmlConv->IsLast())
			xmlTextWriterStartElementNS(writer(), prefix, C_WRAPPER, uri);
		uri=NULL; //not needed again
	}

	xmlTextWriterStartElementNS(writer(), prefix, C_REACTION, NULL);
	if(!pReact->title.empty())
		xmlTextWriterWriteFormatAttribute(writer(), C_TITLE,"%s", pReact->title.c_str());

	xmlTextWriterStartElementNS(writer(), prefix, C_REACTANTLIST, NULL);
  vector<OBMol*>::iterator itr;

  for(itr=pReact->reactants.begin();itr!=pReact->reactants.end();itr++)
  {
		xmlTextWriterStartElementNS(writer(), prefix, C_REACTANT, NULL);
		if(list) //put molecules into map and output references
		{
			string id = AddMolToList(itr);
			xmlTextWriterStartElementNS(writer(), prefix, C_MOLECULE, NULL);
			xmlTextWriterWriteFormatAttribute(writer(), C_REF,"%s", id.c_str());
			xmlTextWriterEndElement(writer());//molecule
		}
		else
		{
			//Write reactant in CML format
			pCMLFormat->WriteMolecule(*itr, _pxmlConv);//&MolConv);
		}
		xmlTextWriterEndElement(writer());//reactant
	}

	xmlTextWriterEndElement(writer());//reactantList
	xmlTextWriterStartElementNS(writer(), prefix, C_PRODUCTLIST, NULL);

  for(itr=pReact->products.begin();itr!=pReact->products.end();itr++)
  {
		xmlTextWriterStartElementNS(writer(), prefix, C_PRODUCT, NULL);
		if(list) //put molecules into map and output references
		{
			string id = AddMolToList(itr);
			xmlTextWriterStartElementNS(writer(), prefix, C_MOLECULE, NULL);
			xmlTextWriterWriteFormatAttribute(writer(), C_REF,"%s", id.c_str());
			xmlTextWriterEndElement(writer());//molecule
		}
		else
			//Write product in MOL format
			pCMLFormat->WriteMolecule(*itr, _pxmlConv);//&MolConv);
		xmlTextWriterEndElement(writer());//product
 }

	xmlTextWriterEndElement(writer());//productList
	xmlTextWriterEndElement(writer());//reaction

	if(pConv->IsLast())
	{
		if(list)
		{
			OutputToStream(); //flush the reactionList to ssout 
			//Back to original output stream
			_pxmlConv->SetOutStream(pOut);
			*pOut << ">\n";

			//output moleculeList
			xmlTextWriterStartElementNS(writer(), prefix, C_MOLECULELIST, NULL);
			map<string,OBMol>::iterator mapitr;
			for(mapitr=OMols.begin();mapitr!=OMols.end();++mapitr)
				pCMLFormat->WriteMolecule(&mapitr->second, _pxmlConv);

			xmlTextWriterEndElement(writer());//moleculeList

			xmlTextWriterStartElementNS(writer(), prefix, C_REACTIONLIST, NULL);
			OutputToStream(); //flush the moleculeList to pOut 
			//output delayed reactionList
			*pOut << ssout.str();
			pOut->seekp(-2,ios::cur); //lose last ">"
//			*pOut <<'\n';
			xmlTextWriterFullEndElement(writer());//reactionList

			xmlTextWriterEndElement(writer());//LISTWRAPPER
		}
		else if(_pxmlConv->GetOutputIndex()>1)
			xmlTextWriterEndElement(writer());//WRAPPER

		xmlTextWriterEndDocument(writer());
		OutputToStream();
	}
	return true;
}

string CMLReactFormat::AddMolToList(vector<OBMol*>::iterator itr)
{
	//Adds a molecule to the map
	string id = (*itr)->GetTitle();
	map<string,OBMol>::iterator mapitr;
	if(!id.empty())
		mapitr = OMols.find(id);
	if(id.empty() || mapitr==OMols.end())
	{
		//not in map; need to add
		if(id.empty())
		{
			//no id, so make one
			stringstream ssid;
			ssid << "m" << nextmol++;
			id = ssid.str();
			//TODO Could check map for chemically identical molecules
			//so that diverse sources could be aggregated
			(*itr)->SetTitle(id);
		}
		OMols[id] = **itr;
	}
	return id;
}


} //namespace OpenBabel

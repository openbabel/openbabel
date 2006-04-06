/**********************************************************************
obmolecformat.cpp - Implementation of subclass of OBFormat for conversion of OBMol.

Copyright (C) 2005 Chris Morley

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
#include "babelconfig.h"
#include "obmolecformat.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

using namespace std;
namespace OpenBabel
{

std::map<std::string, OBMol*> OBMoleculeFormat::IMols;
OBMol* OBMoleculeFormat::_jmol;

bool OBMoleculeFormat::ReadChemObjectImpl(OBConversion* pConv, OBFormat* pFormat)
{
  std::istream &ifs = *pConv->GetInStream();
  if (ifs.peek() == EOF || !ifs.good())
    return false;

	OBMol* pmol = new OBMol;

	std::string auditMsg = "OpenBabel::Read molecule ";
	std::string description(pFormat->Description());
	auditMsg += description.substr(0,description.find('\n'));
	obErrorLog.ThrowError(__FUNCTION__,
			auditMsg,
			obAuditMsg);

	if(pConv->IsOption("C",OBConversion::GENOPTIONS))
		return DeferMolOutput(pmol, pConv, pFormat);

	bool ret=pFormat->ReadMolecule(pmol,pConv);

	OBMol* ptmol = NULL; 
	if(ret && (pmol->NumAtoms() > 0 || (pFormat->Flags()&ZEROATOMSOK))) //Do transformation and return molecule
	{
		ptmol = static_cast<OBMol*>(pmol->DoTransformations(pConv->GetOptions(OBConversion::GENOPTIONS)));
		if(ptmol && pConv->IsOption("j",OBConversion::GENOPTIONS))
		{
			//With j option, accumulate all mols in to one stored in this class
			if(pConv->IsFirstInput())
				_jmol = new OBMol;
			*_jmol += *ptmol;
			delete ptmol;
			return true;
		}
	}
	else
		delete pmol;
	ret = ret && pConv->AddChemObject(ptmol); //success of both writing and reading
	return ret;
}

bool OBMoleculeFormat::WriteChemObjectImpl(OBConversion* pConv, OBFormat* pFormat)
{
	if(pConv->IsOption("C",OBConversion::GENOPTIONS))
		return OutputDeferredMols(pConv);
	if(pConv->IsOption("j",OBConversion::GENOPTIONS))
	{
		bool ret=pFormat->WriteMolecule(_jmol,pConv);
		pConv->SetOutputIndex(1);
		delete _jmol;
		return ret;
	}

	//Retrieve the target OBMol
	OBBase* pOb = pConv->GetChemObject();
	OBMol* pmol = dynamic_cast<OBMol*> (pOb);
	bool ret=false;
	if(pmol)
	{	
		if(pmol->NumAtoms()==0)
		{
			std::string auditMsg = "OpenBabel::Molecule ";
			auditMsg += pmol->GetTitle();
			auditMsg += " has 0 atoms";
			obErrorLog.ThrowError(__FUNCTION__,
					auditMsg,
					obInfo);
		}
		ret=true;

		std::string auditMsg = "OpenBabel::Write molecule ";
		std::string description(pFormat->Description());
		auditMsg += description.substr(0,description.find('\n'));
		obErrorLog.ThrowError(__FUNCTION__,
				      auditMsg,
				      obAuditMsg);

		ret=pFormat->WriteMolecule(pmol,pConv);
	}
	delete pOb; //move so that non-OBMol objects are deleted 9March2006
	return ret;
}
bool OBMoleculeFormat::DeferMolOutput(OBMol* pmol, OBConversion* pConv, OBFormat* pF )
{
	/* Instead of sending molecules for output via AddChemObject(), they are
	   saved in here in OBMoleculeFormat or discarded. By default they are 
		 saved only if they are in the first input file. Parts of subsequent
		 molecules, such as chemical structure, coordinates and OBGenericData
		 can replace the parts in molecules with the same title that have already
		 been stored, subject to a set of rules. After all input files have been
		 read, the stored molecules (possibly now having augmented properties)are
		 sent to the output format.

		 Is a static function with <this> as parameter so that it can be called from other
		 format classes like XMLMoleculeFormat which are not derived from OBMoleculeFormat. 
	*/
	static bool IsFirstFile;
	bool OnlyMolsInFirstFile=true;

	if(pConv->IsFirstInput())
	{
		IsFirstFile=true;
		IMols.clear();
	}
	else 
	{
		if((std::streamoff)pConv->GetInStream()->tellg()<=0)
			IsFirstFile=false;//File has changed
	}

	if (!pF->ReadMolecule(pmol,pConv))
	{
		delete pmol;
		return false;
	}
	const char* ptitle = pmol->GetTitle();
	if(*ptitle==0)
		obErrorLog.ThrowError(__FUNCTION__, "Molecule with no title ignored", obWarning);
	else
	{
		string title(ptitle);
		string::size_type pos = title.find_first_of("\t\r\n"); //some title have other data appended
		if(pos!=string::npos)
			title.erase(pos);
		
		map<std::string, OBMol*>::iterator itr;
		itr = IMols.find(title);
		if(itr!=IMols.end())
		{
			//Molecule with the same title has been input previously: update it
			OBMol* pNewMol = MakeCombinedMolecule(itr->second, pmol);
			if(pNewMol)
			{
				delete itr->second;
				IMols[title] = pNewMol;
			}
			else
			{
				//error: cleanup and return false
				delete pmol;
				return DeleteDeferredMols();
			}
		}
		else
		{
			//Molecule not already saved in IMols: save it if in first file
			if(!OnlyMolsInFirstFile || IsFirstFile)
			{
				IMols[title] = pmol;
				return true; //don't delete pmol
			}
		}
	}
	delete pmol;
	return true;
}

OBMol* OBMoleculeFormat::MakeCombinedMolecule(OBMol* pFirst, OBMol* pSecond)
{
	/* Makes a new OBMol on the heap by combining two molecules according to the rule below. 
	  If both have OBGenericData of the same type, or OBPairData with the
		same attribute,	the version from pFirst is used.
		Returns a pointer to a new OBMol which will need deleting by the calling program
		(probably by being sent to an output format). 
		If the molecules cannot be regarded as being the same structure a NULL
	  pointer is returned and an error message logged.

		pFirst and pSecond and the objects they point to are not changed. (const
		modifiers difficult because class OBMol not designed appropriately)

	  Combining molecules: rules for each of the three parts
		Title:
		Use the title of pFirst unless it has none, when use that of pSecond.
		Warning if neither molecule has a title.

		Structure
	- a structure with atoms replaces one with no atoms
	- a structure with bonds replaces one with no bonds,
	  provided the formula is the same, else an error.
	- structures with atoms and bonds are compared by InChI; error if not the same. 
	- a structure with 3D coordinates replaces one with 2D coordinates
	- a structure with 2D coordinates replace one with 0D coordinates

		OBGenericData
		OBPairData
  */
	
	string title("No title");
	if(*pFirst->GetTitle()!=0)
		title = pFirst->GetTitle();
	else
	{
		if(*pSecond->GetTitle()!=0)
			title = pSecond->GetTitle();
		else
			obErrorLog.ThrowError(__FUNCTION__,"Combined molecule has no title", obWarning);
	}

	bool swap=false;
	if(pFirst->NumAtoms()==0 && pSecond->NumAtoms()!=0)
		swap=true;
	else
	{
		if(pFirst->GetSpacedFormula()!=pSecond->GetSpacedFormula())
		{
			obErrorLog.ThrowError(__FUNCTION__, 
				"Molecules with name = " + title + " have different formula",obError);
			return NULL;
		}
		else
		{
			if(pSecond->NumBonds()!=0 && pFirst->NumBonds()==0)
				swap=true;
			else
			{
				//Compare by inchi; error if different NOT YET IMPLEMENTED
				//Use the one with the higher dimension
				if(pSecond->GetDimension() > pFirst->GetDimension())
					swap=true;
			}
		}
	}

	OBMol* pNewMol = new OBMol;
	pNewMol->SetTitle(title);

	OBMol* pMain = swap ? pSecond : pFirst;
	OBMol* pOther = swap ? pFirst : pSecond;
		
	*pNewMol = *pMain; //Now copies all data 

	//Copy some OBGenericData from the OBMol which did not provide the structure
	vector<OBGenericData*>::iterator igd;
	for(igd=pOther->BeginData();igd!=pOther->EndData();++igd)
	{
		//copy only if not already data of the same type from molecule already copied to pNewMol
		unsigned datatype = (*igd)->GetDataType();
		OBGenericData* pData = pNewMol->GetData(datatype);
		if(datatype==OBGenericDataType::PairData)
		{
			if(pData->GetAttribute() == (*igd)->GetAttribute())
				continue;
		}
		else if(pNewMol->GetData(datatype)!=NULL)
			continue;

		OBGenericData* pCopiedData = (*igd)->Clone(pNewMol);
		pNewMol->SetData(pCopiedData);
	}
	return pNewMol;
}

bool OBMoleculeFormat::OutputDeferredMols(OBConversion* pConv)
{
	std::map<std::string, OBMol*>::iterator itr, lastitr;
	bool ret=false;
	int i=1;
	lastitr = IMols.end();
	--lastitr;
	pConv->SetOneObjectOnly(false);
	for(itr=IMols.begin();itr!=IMols.end();++itr,++i)
	{
		if(!(itr->second)->DoTransformations(pConv->GetOptions(OBConversion::GENOPTIONS)))
			continue;
		pConv->SetOutputIndex(i);
		if(itr==lastitr)
			pConv->SetOneObjectOnly(); //to set IsLast

		std::string auditMsg = "OpenBabel::Write molecule ";
		std::string description((pConv->GetOutFormat())->Description());
		auditMsg += description.substr(0,description.find('\n'));
		obErrorLog.ThrowError(__FUNCTION__, auditMsg,  obAuditMsg);

		ret = pConv->GetOutFormat()->WriteMolecule(itr->second, pConv);

		delete itr->second; //always delete OBMol object
		itr->second = NULL; // so can be deleted in DeleteDeferredMols()
		if (!ret) break;
	}
	DeleteDeferredMols();//cleans up in case there have been errors
	return ret;
}

bool OBMoleculeFormat::DeleteDeferredMols()
{
	//Empties IMols, deteting the OBMol objects whose pointers are stored there 
	std::map<std::string, OBMol*>::iterator itr;
	for(itr=IMols.begin();itr!=IMols.end();++itr)
	{
		delete itr->second; //usually NULL
	}
	IMols.clear();
	return false;
}


} //namespace OpenBabel
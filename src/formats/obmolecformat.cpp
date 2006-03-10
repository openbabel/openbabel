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
			return true;
		}
	}
	else
		delete pmol;
	pConv->AddChemObject(ptmol);
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
			OBMol* pNewMol = MakeCombinedMolecule(itr->second,pmol);
			if(pNewMol)
			{
				delete itr->second;
				IMols[title] = pNewMol;
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
	  The OBMol object needs to be deleted at some stage.
		If the molecules cannot be regarded as being the same structure a NULL
	  pointer is returned and an error message logged. 
	
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
	
	OBMol* pNewMol = new OBMol(*pFirst); //Copies the whole molecule (NOT TRUE CURRENTLY)
	string title("No title");
	if(*pFirst->GetTitle()==0)
	{
		pNewMol->SetTitle(pSecond->GetTitle());
		if(*pSecond->GetTitle()==0)
			obErrorLog.ThrowError(__FUNCTION__,"Combined molecule has no title", obWarning);
		else
			title = pNewMol->GetTitle();
	}

	bool swap=false;
	if(pFirst->NumAtoms()==0 && pSecond->NumAtoms()!=0)
		swap=true;
	else
	{
		//Make hydrogens explicit on copies of molecules
		//in order to do comparisons; leave originals as they were.
		OBMol pCopy1(*pSecond);
		OBMol pCopy2(*pFirst);
		if(pCopy1.NumBonds()!=0 || pCopy1.NumAtoms()==1)
			pCopy1.AddHydrogens(false,false);
		if(pCopy2.NumBonds()!=0 || pCopy2.NumAtoms()==1)
			pCopy2.AddHydrogens(false,false);

		if(pCopy1.GetSpacedFormula()!=pCopy2.GetSpacedFormula())
		{
			obErrorLog.ThrowError(__FUNCTION__, 
				title + "cannot be combined with a molecule having a different formula",obError);
			delete pNewMol;
			return NULL;
		}
		else
		{
			if(pCopy2.NumBonds()==0 && pCopy1.NumBonds()!=0)
				swap=true;
			else
			{
				//Compare by inchi; error if different
				//Use the one with the higher dimension
				if(pCopy2.GetDimension() > pCopy1.GetDimension())
					swap=true;
			}
		}
	}
	if(swap)
	{
		*pNewMol = *pSecond; //Doesn't copy all data (yet) 
		pNewMol->SetTitle(title); //may have been overwritten
	}

	vector<OBGenericData*>::iterator igd;
	//Copy OBGenericData from pFirst; may not be necessary when OBMol constructor updated
	for(igd=pFirst->BeginData();igd!=pFirst->EndData();++igd)
	{
		if((*igd)->GetDataType() == OBGenericDataType::RotamerList)//currently copied in OBMol=
			continue;
		OBGenericData* pCopiedData = (*igd)->Clone();
		pNewMol->SetData(pCopiedData);
	}

	//Copy some OBGenericData from pSecond
	for(igd=pSecond->BeginData();igd!=pSecond->EndData();++igd)
	{
		//copy only if not already data of the same type from pFirst
		unsigned datatype = (*igd)->GetDataType();
		OBGenericData* pData = pNewMol->GetData(datatype);
		if(datatype==OBGenericDataType::PairData)
		{
			if(pData->GetAttribute() == (*igd)->GetAttribute())
				continue;
		}
		else if(pNewMol->GetData(datatype)!=NULL)
			continue;

		OBGenericData* pCopiedData = (*igd)->Clone();
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
		ret = pConv->GetOutFormat()->WriteMolecule(itr->second, pConv);
		if (!ret) break;
	}
	IMols.clear();
	return ret;
}
}
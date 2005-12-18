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
#include "obmolecformat.h"

namespace OpenBabel
{

std::map<std::string, OBMol*> OBMoleculeFormat::IMols;

bool OBMoleculeFormat::ReadChemObject(OBConversion* pConv)
{
  std::istream &ifs = *pConv->GetInStream();
  if (ifs.peek() == EOF || !ifs.good())
    return false;

	static OBMol* pmol;

	    std::string auditMsg = "OpenBabel::Read molecule ";
	    std::string description(Description());
	    auditMsg += description.substr(0,description.find('\n'));
	    obErrorLog.ThrowError(__FUNCTION__,
				  auditMsg,
				  obAuditMsg);

	//With j option, reuse pmol except for the first mol
	if(!pConv->IsOption("j",OBConversion::GENOPTIONS) || pConv->IsFirstInput())
		pmol = new OBMol;
	
	if(pConv->IsOption("C",OBConversion::GENOPTIONS))
		return DeferMolOutput(pmol, pConv, this);

	bool ret=ReadMolecule(pmol,pConv);
	
	if(ret && (pmol->NumAtoms() > 0 || (Flags()&ZEROATOMSOK))) //Do transformation and return molecule
		pConv->AddChemObject(pmol->DoTransformations(pConv->GetOptions(OBConversion::GENOPTIONS)));
	else
		pConv->AddChemObject(NULL);

	return ret;
}

bool OBMoleculeFormat::WriteChemObject(OBConversion* pConv)
{
	if(pConv->IsOption("C",OBConversion::GENOPTIONS))
		return OutputDeferredMols(pConv);

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
		std::string description(Description());
		auditMsg += description.substr(0,description.find('\n'));
		obErrorLog.ThrowError(__FUNCTION__,
				      auditMsg,
				      obAuditMsg);

		if(!pConv->IsOption("j",OBConversion::GENOPTIONS) || pConv->IsLast()) //With j option, output only at end
		{
			ret=WriteMolecule(pmol,pConv);
			delete pOb;
		}
	}
	return ret;
}
bool OBMoleculeFormat::DeferMolOutput(OBMol* pmol, OBConversion* pConv, OBFormat* pF )
{
	/* Instead of sending molecules for output vis AddChemObject(), they are
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
	if(!ptitle || !*ptitle)
	{
		obErrorLog.ThrowError(__FUNCTION__, "Molecule with no title ignored", obWarning);
		delete pmol;
		return true;
	}

	std::string title(ptitle);
	std::string::size_type pos = title.find_first_of("\t\r\n"); //some title have other data appended
	if(pos!=std::string::npos)
		title.erase(pos);
	std::map<std::string, OBMol*>::iterator itr;
	itr = IMols.find(title);
	if(itr==IMols.end())
	{
		if((!OnlyMolsInFirstFile || IsFirstFile))
			IMols[title] = pmol; //add molecule not already saved in IMols
		else
			delete pmol;
	}
	else
	{
		OBMol* pold = itr->second;
		/* 
		- ignore molecules with no title
		Rules for handling molecules with the same title
		- a structure with atoms replaces one with no atoms
		- a structure with bonds replaces one with no bonds,
		  provided the formula is the same, else an error.
		- structures with atoms and bonds are compared by InChI; error if not the same. 
		- Any 3D coordinates replace 2D coordinates
		- Any 2D coordinates replace 0D coordinates(0)
		- PairData and other OBGenericData with new attributes are copied; 
		  existing attributes have their values replaced.
    */
		bool swap=false;
		if(pold->NumAtoms()==0 && pmol->NumAtoms()!=0)
			swap=true;
		else
		{
			//Make hydrogens explicit on copies of molecules
			//in order to do comparisons; leave originals as they were.
			OBMol Copymol(*pmol);
			OBMol Copyold(*pold);
			if(Copymol.NumBonds()!=0 || Copymol.NumAtoms()==1)
				Copymol.AddHydrogens(false,false);
			if(Copyold.NumBonds()!=0 || Copyold.NumAtoms()==1)
				Copyold.AddHydrogens(false,false);

			if(Copymol.GetSpacedFormula()!=Copyold.GetSpacedFormula())
			obErrorLog.ThrowError(__FUNCTION__, 
			"Two molecules with id \"" + std::string(ptitle) + "\" have different formula", obWarning);
			else
			{
				if(Copyold.NumBonds()==0 && Copymol.NumBonds()!=0)
					swap=true;
				else
				{
					//Compare by inchi; error if different
					//Use the one with the higher dimension NOT CORRECT!!
					if(Copymol.GetDimension() > Copyold.GetDimension())
						swap=true;
				}
			}
		}
		if(swap)
		{
			IMols[title]=pmol;
			std::swap(pold, pmol);
		}

		//Copy OBGenericData. Currently previous existence not tested
		std::vector<OBGenericData*>::iterator igd;
		for(igd=pmol->BeginData();igd!=pmol->EndData();++igd)
		{
			OBGenericData* pCopiedData = (*igd)->Clone();
			pold->SetData(pCopiedData);
		}

		delete pmol;
	}
	return true;
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
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
#include "mol.h"
#include "obconversion.h"
#include "obmolecformat.h"
#include "fingerprint.h"
#include <vector>
#include <string>
#include <iomanip>

using namespace std;
namespace OpenBabel {

class FingerprintFormat : public OBMoleculeFormat
{
public:
	//Register this format type ID
	FingerprintFormat() {OBConversion::RegisterFormat("fpt",this);}

	virtual const char* Description() //required
	{ return
"Fingerprint format\n \
See Fabien Fontaine's source code\n \
Options e.g. -xn16\n \
 f# finger print type, default <2>\n \
 n# fold to specified number of 32bit words \n \
 h  hex output when multiple molecules\n \
";
	};

	virtual unsigned int Flags(){return NOTREADABLE;};
	virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

private:
	vector<unsigned int> firstfp;
	string firstname;
	bool IsPossibleSubstructure(vector<unsigned int>Mol, vector<unsigned int>Frag);
};

////////////////////////////////////////////////////
//Make an instance of the format class
FingerprintFormat theFingerprintFormat;

//*******************************************************************
bool FingerprintFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
	OBMol* pmol = dynamic_cast<OBMol*>(pOb);
	if(pmol==NULL) return false;

	//Define some references so we can use the old parameter names
	ostream &ofs = *pConv->GetOutStream();
	OBMol &mol = *pmol;

	bool hexoutput=false;
	if(pConv->IsOption('h') || (pConv->GetOutputIndex()==1 && pConv->IsLast()))
		hexoutput=true;

	int type=0, nwords=0;
	const char* p=pConv->IsOption('f');
	if(p)
		type=atoi(p);
	p=pConv->IsOption('n');
	if(p)
		nwords = atoi(p);		

	vector<unsigned int> fptvec;
	if(!GetFingerprint(mol, fptvec, nwords, type))
		return false;
	
  ofs << ">" << mol.GetTitle();
	if(hexoutput)
	{
		int i, bitsset=0;
		for (i=0;i<fptvec.size();++i)
		{
		int wd = fptvec[i];
		for(;wd;wd=wd<<1)//count bits set by shifting into sign bit until word==0
			if(wd<0) ++bitsset;
		}
		ofs  << "   " << bitsset << " bits set. "; 
	}

	if(pConv->GetOutputIndex()==1)
	{
		//store the fingerprint and name of first molecule
		firstfp=fptvec;
		firstname=mol.GetTitle();
		if(firstname.empty())
			firstname = "first mol";
	}
	else
	{
		ofs << "   Tanimoto from " << firstname << " = " << Tanimoto(firstfp, fptvec);
		if(IsPossibleSubstructure(fptvec,firstfp))
			ofs << "\nPossible superstructure of " << firstname;
	}
	ofs << endl;
	
	int i;

	if(hexoutput)
	{
		for(i=fptvec.size()-1;i>=0;i--)
		{
			ofs << hex << setfill('0') << setw(8) << fptvec[i] << " " ;
			if((fptvec.size()-i)%6==0)
				ofs <<endl;
		}
		ofs << dec << endl;
	}
	return true;
}

bool FingerprintFormat::IsPossibleSubstructure(vector<unsigned int>Mol, vector<unsigned int>Frag)
{
	//Returns false if Frag is definitely NOT a substructure of Mol
	int i;
	for (i=0;i<Mol.size();++i)
		if((Mol[i] & Frag[i]) ^ Frag[i]) return false;
	return true;
}

} //namespace OpenBabel

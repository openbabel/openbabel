/**********************************************************************
Copyright (C) 2003 by Dr. Alex M. Clark and Geoffrey Hutchison
Some portions Copyright (C) 2004 by Chris Morley

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
namespace OpenBabel {

class CRK2DFormat : public OBFormat
{
public:
	//Register this format type ID
	CRK2DFormat() {OBConversion::RegisterFormat("CRK2D",this);}

	virtual const char* Description() //required
	{ return
"CRK Diagram Structure (2D)\n \
No comments yet\n \
";
	};

	virtual const char* SpecificationURL(){return
		"http://??";}; //optional

	//Flags() can return be any the following combined by | or be omitted if none apply
	// NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY 
	virtual unsigned int Flags(){return READONEONLY;};

//*** This section identical for most OBMol conversions ***
	////////////////////////////////////////////////////
	/// The "API" interface functions
	virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
	virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

////////////////////////////////////////////////////
	/// The "Convert" interface functions
	virtual bool ReadChemObject(OBConversion* pConv)
	{
		OBMol* pmol = new OBMol;
		bool ret=ReadMolecule(pmol,pConv);
		if(ret) //Do transformation and return molecule
			pConv->AddChemObject(pmol->DoTransformations(pConv->GetGeneralOptions()));
		else
			pConv->AddChemObject(NULL);
		return ret;
	};
	
	virtual bool WriteChemObject(OBConversion* pConv)
	{
		//Retrieve the target OBMol
		OBBase* pOb = pConv->GetChemObject();
		OBMol* pmol = dynamic_cast<OBMol*> (pOb);
		bool ret=false;
		if(pmol)
			ret=WriteMolecule(pmol,pConv);
		delete pOb; 
		return ret;
	};

	static bool CRK2DFormat::ReadCRK(std::istream &ifs,OBMol &mol,const char *classTag);
  static void CRK2DFormat::WriteCRK(std::ostream &ofs,OBMol &mol,bool GroupCharges);

};


//Make an instance of the format class
CRK2DFormat theCRK2DFormat;

/////////////////////////////////////////////////////////////////
bool CRK2DFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
{

	OBMol* pmol = dynamic_cast<OBMol*>(pOb);
	if(pmol==NULL) return false;

	//Define some references so we can use the old parameter names
	istream &ifs = *pConv->GetInStream();
	OBMol &mol = *pmol;
	const char* title = pConv->GetTitle();

	char buffer[BUFF_SIZE];//CM extra buffer

	if (!ifs.getline(buffer,BUFF_SIZE)) {printf("File is empty!\n"); return(false); }
	if (!strstr(buffer,"<Property")) {printf("Not valid CRK XML.\n"); return false;}
	if (!strstr(buffer,"\"DiagramStructure\"")) {printf("Not CRK DiagramStructure (2D).\n"); return false;}

	return ReadCRK(ifs,mol,"<Structure2D>");
}

////////////////////////////////////////////////////////////////

bool CRK2DFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
	OBMol* pmol = dynamic_cast<OBMol*>(pOb);
	if(pmol==NULL) return false;

	//Define some references so we can use the old parameter names
	ostream &ofs = *pConv->GetOutStream();
	OBMol &mol = *pmol;
	const char *dimension = pConv->GetDimension();

	ofs << "<Property Type=\"DiagramStructure\">" <<  endl;
	ofs << " <Structure2D>" << endl;
	
	WriteCRK(ofs,mol,true);
	
	ofs << " </Structure2D>" << endl;
	ofs << "</Property>" << endl;

	return true;
}

//******************************************************
class CRK3DFormat : public OBFormat
{
public:
	//Register this format type ID
	CRK3DFormat() {OBConversion::RegisterFormat("CRK3D",this);}

	virtual const char* Description() //required
	{ return
"CRK(3D) format\n \
No comments yet\n \
";
	};

	virtual const char* SpecificationURL(){return
		"http://??";}; //optional

	//Flags() can return be any the following combined by | or be omitted if none apply
	// NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY 
	virtual unsigned int Flags(){return READONEONLY;};

//*** This section identical for most OBMol conversions ***
	////////////////////////////////////////////////////
	/// The "API" interface functions
	virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
	virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

////////////////////////////////////////////////////
	/// The "Convert" interface functions
	virtual bool ReadChemObject(OBConversion* pConv)
	{
		OBMol* pmol = new OBMol;
		bool ret=ReadMolecule(pmol,pConv);
		if(ret) //Do transformation and return molecule
			pConv->AddChemObject(pmol->DoTransformations(pConv->GetGeneralOptions()));
		else
			pConv->AddChemObject(NULL);
		return ret;
	};
	
	virtual bool WriteChemObject(OBConversion* pConv)
	{
		//Retrieve the target OBMol
		OBBase* pOb = pConv->GetChemObject();
		OBMol* pmol = dynamic_cast<OBMol*> (pOb);
		bool ret=false;
		if(pmol)
			ret=WriteMolecule(pmol,pConv);
		delete pOb; 
		return ret;
	};
};
//***

//Make an instance of the format class
CRK3DFormat theCRK3DFormat;

/////////////////////////////////////////////////////////////////
bool CRK3DFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
{

	OBMol* pmol = dynamic_cast<OBMol*>(pOb);
	if(pmol==NULL) return false;

	//Define some references so we can use the old parameter names
	istream &ifs = *pConv->GetInStream();
	OBMol &mol = *pmol;
	const char* title = pConv->GetTitle();

	char buffer[BUFF_SIZE];//CM extra buffer

	if (!ifs.getline(buffer,BUFF_SIZE)) {printf("File is empty!\n"); return(false); }
	if (!strstr(buffer,"<Property")) {printf("Not valid CRK XML.\n"); return false;}
	if (!strstr(buffer,"\"ModelStructure\"") && !strstr(buffer,"\"XRayStructure\"")) {printf("Not CRK ModelStructure or XRayStructure (3D).\n"); return false;}

	return CRK2DFormat::ReadCRK(ifs,mol,"<Structure3D<");
}

////////////////////////////////////////////////////////////////

bool CRK3DFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
	OBMol* pmol = dynamic_cast<OBMol*>(pOb);
	if(pmol==NULL) return false;

	//Define some references so we can use the old parameter names
	ostream &ofs = *pConv->GetOutStream();
	OBMol &mol = *pmol;
	const char *dimension = pConv->GetDimension();

	ofs << "<Property Type=\"ModelStructure\">" <<  endl;
	ofs << " <Structure3D>" << endl;
	
	CRK2DFormat::WriteCRK(ofs,mol,true);
	
	ofs << " </Structure3D>" << endl;
	ofs << "</Property>" << endl;

	return true;
}

//**************************************************************
bool CRK2DFormat::ReadCRK(std::istream &ifs,OBMol &mol,const char *classTag)
{
	bool foundClass=false;
	
	#define MAX_ATOMS 1000
	int numAtoms=0;
	int statomID[MAX_ATOMS];
	
	#define MAX_BONDS 1000
	int numBonds=0;
	int stbondFrom[MAX_BONDS],stbondTo[MAX_BONDS],stbondStyle[MAX_BONDS];
	double stbondOrder[MAX_BONDS];
	
	bool inAtom=false,inBond=false;
	int atomID,atomNumber;
	double atomX,atomY,atomZ,atomCharge;
	int bondFrom,bondTo,bondStyle;
	double bondOrder;
	char buffer[BUFF_SIZE];//was global
	
	mol.BeginModify();

	while (ifs.getline(buffer,BUFF_SIZE)) 
	{
		if (strstr(buffer,classTag)) foundClass=true;
		else if (strstr(buffer,"<Atom"))
		{
			atomID=0;
			char *tag=strstr(buffer,"ID=\"");
			if (tag) atomID=atoi(tag+4);
			if (atomID>0)
			{
				inAtom=true;
				atomNumber=0;
				atomX=atomY=atomZ=atomCharge=0;
			}
		}
		else if (strstr(buffer,"<Bond"))
		{
			inBond=true;
			bondFrom=bondTo=bondStyle=0;
			bondOrder=0;
		}
		else if (strstr(buffer,"</Atom>"))
		{
			if (inAtom && numAtoms<MAX_ATOMS)
			{
				OBAtom atm; atm.Clear();

				statomID[numAtoms++]=atomID;

				atm.SetAtomicNum(atomNumber);
				atm.SetVector(atomX,atomY,atomZ);
				atm.SetFormalCharge((int)atomCharge);

				if (!mol.AddAtom(atm)) {printf("Unable to add atom.\n"); return false;}
			}
			inAtom=false;
		}
		else if (strstr(buffer,"</Bond>"))
		{
			if (inBond && numBonds<MAX_BONDS)
			{
				stbondFrom[numBonds]=bondFrom;
				stbondTo[numBonds]=bondTo;
				stbondOrder[numBonds]=bondOrder;
				stbondStyle[numBonds]=bondStyle;				
				numBonds++;
			}
			inBond=false;
		}
		else
		{
			char *tag;
			if (inAtom)
			{
				tag=strstr(buffer,"<X>"); if (tag) atomX=atof(tag+3);
				tag=strstr(buffer,"<Y>"); if (tag) atomY=atof(tag+3);
				tag=strstr(buffer,"<Z>"); if (tag) atomZ=atof(tag+3);
				tag=strstr(buffer,"<Element>"); if (tag) 
				{
					char element[3]="\0\0";
					element[0]=tag[9];
					if (tag[10]>='a' && tag[10]<='z') element[1]=tag[10];
					atomNumber=etab.GetAtomicNum(element);
				}
				tag=strstr(buffer,"<Charge>"); if (tag) atomCharge=atof(tag+8);
			}
			if (inBond)
			{
				tag=strstr(buffer,"<From>"); if (tag) bondFrom=atoi(tag+6);
				tag=strstr(buffer,"<To>"); if (tag) bondTo=atoi(tag+4);
				tag=strstr(buffer,"<Order>"); if (tag) bondOrder=atof(tag+7);
				tag=strstr(buffer,"<Style>"); if (tag) bondStyle=atoi(tag+7);
			}
		}
	}
	
	for(int n=0;n<numBonds;n++)
	{
		int fromIdx=0,toIdx=0;
		for(int i=0;i<numAtoms;i++)
		{
			if (stbondFrom[n]==statomID[i]) fromIdx=i+1;
			if (stbondTo[n]==statomID[i]) toIdx=i+1;
		}
		
		if (fromIdx>0 && toIdx>0)
		{
			OBAtom *from=mol.GetAtom(fromIdx),*to=mol.GetAtom(toIdx);

			int order=1;
			if (stbondOrder[n]==2) order=2;
			else if (stbondOrder[n]==3) order=3;
			else if (stbondOrder[n]==1.5) order=5;

			OBBond bnd;
			bnd.Set(n+1,from,to,order,0);
			
			if (stbondStyle[n]==1) bnd.SetUp();
			if (stbondStyle[n]==2) bnd.SetDown();
			if (stbondOrder[n]==1.5) bnd.SetAromatic();
	
			if (!mol.AddBond(bnd)) {printf("Unable to add bond.\n"); return false;}
		}
		else {printf("Unassigned bond ID (%d,%d), source may be invalid.\n",stbondFrom[n],stbondTo[n]); return false;}
	}
	
   mol.EndModify();


	return foundClass;
}

void CRK2DFormat::WriteCRK(std::ostream &ofs,OBMol &mol,bool GroupCharges)
{
	double groupCharge=0;
	if (GroupCharges)
		for(int n=1;n<=mol.NumAtoms();n++) groupCharge+=mol.GetAtom(n)->GetFormalCharge();

	ofs << "  <Group Charge=\"" << groupCharge << "\" Spin=\"0\">" << endl;

	for(int n=1;n<=mol.NumAtoms();n++)
	{
		OBAtom *atm=mol.GetAtom(n);
		
		int id=atm->GetIdx(),atomnum=atm->GetAtomicNum();
		double x=atm->GetX(),y=atm->GetY(),z=atm->GetZ();
		char *element=etab.GetSymbol(atomnum);
		double charge=0;
		if (!GroupCharges) charge=atm->GetFormalCharge();
		
		//ofs << ((int)atm) << endl;
		
		ofs << "   <Atom ID=\"" << id << "\">" << endl;
		ofs << "    <X>" << x << "</X>" << endl;
		ofs << "    <Y>" << y << "</Y>" << endl;
		ofs << "    <Z>" << z << "</Z>" << endl;
		ofs << "    <Element>" << element << "</Element>" << endl;
		if (charge!=0) ofs << "    <Charge>" << charge << "</Charge>" << endl;
		ofs << "   </Atom>" << endl;
	}
	
	for(int m=0;m<mol.NumBonds();m++)
	{
		OBBond *bnd=mol.GetBond(m);
		
		int from=bnd->GetBeginAtom()->GetIdx(),to=bnd->GetEndAtom()->GetIdx();
		double order=bnd->GetBO();
		if (bnd->IsAromatic()) order=1.5;
		int style=0;
		if (bnd->IsUp()) style=1;
		if (bnd->IsDown() || bnd->IsWedge()) style=2;
		
		ofs << "   <Bond>" << endl;
		ofs << "    <From>" << from << "</From>" << endl;
		ofs << "    <To>" << to << "</To>" << endl;
		ofs << "    <Order>" << order << "</Order>" << endl;
		ofs << "    <Style>" << style << "</Style>" << endl;
		ofs << "   </Bond>" << endl;
	}

	ofs << "  </Group>" << endl;
}

} //namespace OpenBabel

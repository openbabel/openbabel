/**********************************************************************
Copyright (C) 2005 by Chris Morley
 
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
#include "mol.h"
#include "obconversion.h"
#include "xml.h"
#include "math/matrix3x3.h"

#ifdef WIN32
#pragma warning (disable : 4800)
#endif

using namespace std;
namespace OpenBabel
{

class CMLFormat : public XMLMoleculeFormat
{
private:
	const char* CML1NamespaceURI()const{return "http://www.xml-cml.org/dtd/cml_1_0_1.dtd";}

public:
  //Constuctor used on startup which registers this format type ID
	CMLFormat() 
  {
      OBConversion::RegisterFormat("cml", this, "chemical/x-cml");
			OBConversion::RegisterOptionParam("1", this);
			OBConversion::RegisterOptionParam("a", this);
			OBConversion::RegisterOptionParam("N", this, 1);
			OBConversion::RegisterOptionParam("m", this);
			OBConversion::RegisterOptionParam("x", this);
			OBConversion::RegisterOptionParam("h", this);

			XMLConversion::RegisterXMLFormat(this, true);	//this is the default XLMformat
			XMLConversion::RegisterXMLFormat(this, false, CML1NamespaceURI());//CML1 also
  }
	virtual const char* NamespaceURI()const{return "http://www.xml-cml.org/schema/cml2/core";}

  virtual const char* Description()
  {
      return " \
Chemical Markup Language\n \
XML format. This implementation uses libxml2.\n \
Write options for CML: -x[flags] (e.g. -x1ac)\n \
1  output CML1 (rather than CML2)\n \
a  output array format for atoms and bonds\n \
h  use hydrogenCount for all hydrogens\n \
m  output metadata\n \
x  omit XML and namespace declarations\n \
N<prefix> add namespace prefix to elements\n \
\n";
};

  virtual const char* SpecificationURL()
  {return "http://wwmm.ch.cam.ac.uk/moin/ChemicalMarkupLanguage";};

  virtual const char* GetMIMEType() 
  { return "chemical/x-cml"; };

  virtual unsigned int Flags()
  {
      return 0;
  };

////////////////////////////////////////////////////////////////////

	virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

protected:
	virtual bool DoElement(const string& name);
	virtual bool EndElement(const string& name);
	virtual const char* EndTag(){ return "/molecule>"; };
private:
	typedef vector< vector< pair<string,string> > > cmlArray;
	bool TransferArray(cmlArray& arr);
	bool TransferElement(cmlArray& arr);
	bool DoAtoms();
	bool DoBonds();
	bool DoMolWideData();
	bool ParseFormula(string& formula, OBMol* pmol);
	
	void WriteFormula(OBMol& mol);
	void WriteMetadataList();
	string getTimestr();
	void WriteBondStereo(OBBond* pbond);
	void WriteCrystal(OBMol& mol);

private:
	map<string,int> AtomMap; //key=atom id, value= ob atom index
	cmlArray AtomArray;
	cmlArray BondArray;
	vector< pair<string,string> > cmlBondOrAtom; //for cml1 only
	vector< pair<string,string> > molWideData;
	bool inBondArray; //for cml1 only
	string RawFormula;
	xmlChar* prefix;
	string CurrentAtomID;
	int CrystalScalarsNeeded;
	vector<double> CrystalVals;
	OBUnitCell* pUnitCell;
};

////////////////////////////////////////////////////////////
//Make an instance of the format class
CMLFormat theCMLFormat;


/*
There are 4 CML styles: CML1, CML2, both with and without array forms.
All styles are converted into the same internal structure in AtomArray
and BondArray which contains pairs of (attribute)name/value pairs for 
each atom or bond. At the end of molecule this is analysed in DoAtoms()
and DoBonds() to construct an OBMol.
*/

//Callback routines
///////////////////////////////////////////////////////
bool CMLFormat::DoElement(const string& name)
{
	//A linear search is good enough for <20 element names; commonest at start.
	string value;
	if(name=="atom")
	{
		cmlBondOrAtom.clear();
		int IsEmpty = xmlTextReaderIsEmptyElement(reader());
		TransferElement(AtomArray);
		if(IsEmpty==1) //have to push here because end atom may not be called
			AtomArray.push_back(cmlBondOrAtom);
	}	
	else if(name=="bond")
	{
		cmlBondOrAtom.clear();
		int IsEmpty = xmlTextReaderIsEmptyElement(reader());
		TransferElement(BondArray);
		if(IsEmpty==1)
			BondArray.push_back(cmlBondOrAtom);
	}
	else if(name=="molecule")
	{
		//Ignore atoms with "ref" attributes
		if(xmlTextReaderGetAttribute(reader(), BAD_CAST "ref"))
			return true;
		AtomArray.clear();
		BondArray.clear();
		inBondArray = false;
		RawFormula.erase();
		molWideData.clear();
		CrystalScalarsNeeded=0;
		CrystalVals.clear();
		pUnitCell = NULL;

		if(++_embedlevel)
			return true; //ignore if already inside a molecule
		AtomMap.clear();

		const xmlChar* ptitle  = xmlTextReaderGetAttribute(reader(), BAD_CAST "title");
		if(!ptitle)
			ptitle  = xmlTextReaderGetAttribute(reader(), BAD_CAST "id");
		if(ptitle)
			_pmol->SetTitle((const char*)ptitle);
	}
	else if(name=="atomArray")
	{
		inBondArray=false;
		TransferArray(AtomArray);
	}
	else if(name=="bondArray")
	{
		inBondArray=true;
		TransferArray(BondArray);
	}
	else if(name=="atomParity" || name=="bondStereo")
	{
		//Save in molWideData:
		//the content,  the atomRefs4 attribute, and (for atomParity only) the centralAtom 
		string atrefs4("atomRefs4");
		value = _pxmlConv->GetAttribute(atrefs4.c_str());
		pair<string,string> atomrefdata(atrefs4,value);					

		xmlTextReaderRead(reader());
		const xmlChar* pvalue = xmlTextReaderConstValue(reader());
		if(pvalue)
		{
			value = (const char*)pvalue;
			pair<string,string> nameAndvalue(name,value);					
			molWideData.push_back(nameAndvalue);			
			molWideData.push_back(atomrefdata);
			
			stringstream ss;
			if(name=="atomParity")
				ss << AtomArray.size()+1; //index of current atom
			else
				ss << BondArray.size(); //index of current bond
			pair<string,string> atdata("centralAtomOrBond",ss.str());
			molWideData.push_back(atdata);
		}
	}
	else if(name=="name")
	{
		if(_pmol)
			_pmol->SetTitle(_pxmlConv->GetContent().c_str());
	}
	else if(name=="formula")
	{
		//Only concise form is currently supported
		const xmlChar* pformula = xmlTextReaderGetAttribute(reader(), BAD_CAST "concise");
		if(pformula)
			RawFormula = (const char*)pformula;
	}
	else if(name=="crystal")
	{
		CrystalScalarsNeeded = 6;	
	}
	else if(name=="scalar")
	{
		if(CrystalScalarsNeeded)
		{
			xmlTextReaderRead(reader());
			const xmlChar* pvalue = xmlTextReaderConstValue(reader());
			if(pvalue)
			{
				CrystalVals.push_back(atof((const char*)pvalue));
				if(--CrystalScalarsNeeded==0)
				{
					pUnitCell = new OBUnitCell;
					pUnitCell->SetData(CrystalVals[0],CrystalVals[1],CrystalVals[2],
						CrystalVals[3],CrystalVals[4],CrystalVals[5]);
					_pmol->SetData(pUnitCell);
				}
			}
		}	
	}

	// CML1 elements
	else	if(name=="string" || name=="float" || name=="integer")
	{
		string name = _pxmlConv->GetAttribute("builtin");
		xmlTextReaderRead(reader());
		const xmlChar* pvalue = xmlTextReaderConstValue(reader());
		if(!pvalue)
			return false;
		string value = (const char*)pvalue;
		pair<string,string> nameAndvalue(name,value);					
		cmlBondOrAtom.push_back(nameAndvalue);			
	}
	else	if(name=="stringArray" || name=="floatArray" || name=="integerArray")
	{
		string name = _pxmlConv->GetAttribute("builtin");
//		cmlArray& arr = (name=="atomRef1" || name=="atomRef2" || name=="order") 
//			? BondArray : AtomArray;
		cmlArray& arr = inBondArray ? BondArray : AtomArray;

		xmlTextReaderRead(reader());
		const xmlChar* pvalue = xmlTextReaderConstValue(reader());
		if(!pvalue)
			return false;
		string value = (const char*)pvalue;

		vector<string> items;
		tokenize(items,value);
		if(arr.size()<items.size())
			arr.resize(items.size());
		int i;
		for(i=0;i<items.size();++i)
		{				
			pair<string,string> nameAndvalue(name,items[i]);					
			arr[i].push_back(nameAndvalue);
		}
	}
	return true;
}

//////////////////////////////////////////////////////
bool CMLFormat::EndElement(const string& name)
{
	if(name=="atom")
	{
		//ok for cml1 but is not called at end of <atom.../>
		AtomArray.push_back(cmlBondOrAtom);
	}
	
	if(name=="bond")
	{
		BondArray.push_back(cmlBondOrAtom);
	}

	if(name=="molecule")
	{
		DoAtoms();
		DoBonds();
		DoMolWideData();

		_pmol->AssignSpinMultiplicity();
		
		//Use formula only if nothing else provided
		if(_pmol->NumAtoms()==0 && !RawFormula.empty())
			if(!ParseFormula(RawFormula, _pmol))
				cerr << "Error in formula" << endl;
		
		_pmol->EndModify();
		return (--_embedlevel>=0); //false to stop parsing if no further embedded mols
//		return false;//means stop parsing
	}
	return true;
}

/////////////////////////////////////////////////////////

///Interprets atoms from AtomArray and writes then to an OBMol
bool CMLFormat::DoAtoms()
{	
	int dim=0; //dimension of molecule
	int nAtoms=_pmol->NumAtoms();//was 0
	cmlArray::iterator AtomIter;
	for(AtomIter=AtomArray.begin();AtomIter!=AtomArray.end();++AtomIter)
	{
//		OBAtom obatom;
		OBAtom* pAtom = _pmol->NewAtom();
		nAtoms++;
		int nhvy = nAtoms;

		double x=0,y=0,z=0;
		
		vector<pair<string,string> >::iterator AttributeIter;
		for(AttributeIter=AtomIter->begin();AttributeIter!=AtomIter->end();++AttributeIter)
		{
			string& attrname = AttributeIter->first;
			string& value    = AttributeIter->second;

			if(attrname=="id" || attrname=="atomId" || attrname=="atomID")//which one correct? 
			{
				Trim(value);
				AtomMap[value] = nhvy;//nAtoms;
			}
			else if(attrname=="elementType")
			{
				int atno, iso=0;
				atno=etab.GetAtomicNum(value.c_str(),iso);
				pAtom->SetAtomicNum(atno);
				if(iso)
					pAtom->SetIsotope(iso);
			}
			
			else if(attrname=="x2" && dim!=3)//ignore 2D dimensions if 3D also provided
				x=strtod(value.c_str(),NULL);

			else if(attrname=="y2" && dim!=3)
			{
				dim=2;
				y=strtod(value.c_str(),NULL);
			}

			else if(attrname=="x3" || (pUnitCell && attrname=="xFract"))
				x=strtod(value.c_str(),NULL);

			else if(attrname=="y3" || (pUnitCell && attrname=="yFract"))
				y=strtod(value.c_str(),NULL);

			else if(attrname=="z3" || (pUnitCell && attrname=="zFract"))
			{
				dim=3;
				z=strtod(value.c_str(),NULL);
			}

			else if(attrname=="xy2" && dim!=3)
			{
				dim=2;
				vector<string> vals;
				tokenize(vals,value);
				if(vals.size()==2)
				{
					x=strtod(vals[0].c_str(),NULL);
					y=strtod(vals[1].c_str(),NULL);
				}
			}
			 else if(attrname=="xyz3" || (pUnitCell && attrname=="xyzFract"))
			{
				dim=3;
				vector<string> vals;
				tokenize(vals,value);
				if(vals.size()==3)
				{
					x=strtod(vals[0].c_str(),NULL);
					y=strtod(vals[1].c_str(),NULL);
					z=strtod(vals[2].c_str(),NULL);
				}
			}

			if(dim)
			{
				if(!pUnitCell)
					pAtom->SetVector(x, y , z);
				else
				{
					//Coordinates are fractional
					vector3 v;
					v.Set(x, y, z);
					v *= pUnitCell->GetOrthoMatrix();
					pAtom->SetVector(v);
				}
			}

			if(attrname=="hydrogenCount")
			{
				int nhvy = nAtoms;
				int i;
				for(i=0;i<atoi(value.c_str());++i)
				{
					OBAtom* hatom = _pmol->NewAtom();
					hatom->SetAtomicNum(1);
					hatom->SetType("H");
					_pmol->AddBond(nhvy,_pmol->NumAtoms(),1);
					++nAtoms;
				}
			}

			else if(attrname=="formalCharge") 
				pAtom->SetFormalCharge(atoi(value.c_str()));

			else if(attrname=="spinMultiplicity")
				pAtom->SetSpinMultiplicity(atoi(value.c_str()));

			else if(attrname=="atomRefs4")//from atomParity element
			{
				vector<string> ids;
				tokenize(ids,value);
				// Have 4 atoms defining the parity
				// but don't currently use them TODO
				//Simply use parity as given to set clockwise/anticlockwise

				attrname = (++AttributeIter)->first;
				if(attrname=="parity")
				{
					value = AttributeIter->second;
					int parity = atoi(value.c_str());
					if(parity>0) pAtom->SetClockwiseStereo();
					if(parity<0) pAtom->SetAntiClockwiseStereo();
				}
			}

			else if(attrname=="radical") //Marvin extension
			{
				int spin=0;
				if(value=="monovalent")
					spin=2;
				else if(value=="divalent")
					spin=3;
				else if(value=="divalent3")
					spin=3;
				else if(value=="divalent1")
					spin=1;
				pAtom->SetSpinMultiplicity(spin);
			}
			else if(attrname=="isotopeNumber" || attrname=="isotope")
				pAtom->SetIsotope(atoi(value.c_str()));

		} //each attribute
		
	}//each atom
	
	_pmol->SetDimension(dim);
	return true;
}
/////////////////////////////////////////////////////////////////////

///Interprets bonds from BondArray and writes then to an OBMol
bool CMLFormat::DoBonds()
{
	vector<pair<string,string> >::iterator AttributeIter;
	cmlArray::iterator BondIter;
	for(BondIter=BondArray.begin();BondIter!=BondArray.end();++BondIter)
	{
		int indx1=0,indx2=0, ord=0;
		string bondstereo, BondStereoRefs;

		for(AttributeIter=BondIter->begin();AttributeIter!=BondIter->end();++AttributeIter)
		{
			string attrname = AttributeIter->first;
			string value    = AttributeIter->second;
		
			if(attrname=="atomRefs2")
			{	
				Trim(value);
				string::size_type pos = value.find(' ');
				indx1 = AtomMap[value.substr(0,pos)];	
				indx2 = AtomMap[value.substr(pos+1)];
			}

			else if(attrname=="atomRef1" || (attrname=="atomRef" && indx1==0))
				indx1 = AtomMap[value];

			else if(attrname=="atomRef2"|| attrname=="atomRef")
				indx2 = AtomMap[value];
			
			else if(attrname=="order")
			{	
				Trim(value);
				const char bo = value[0];
				if(bo=='S')
					ord=1;
				else if(bo=='D')
					ord=2;
				else if(bo=='A')
					ord=5;
				else
					ord=atoi(&bo);
			}

		}
		if(indx1==0 || indx2==0)
		{
			cerr << "Incorrect bond attributes" << endl;
			return false;
		}
		if(ord==0) //Bonds are single if order is not specified
			ord=1;
		_pmol->AddBond(indx1,indx2,ord,0);
	}

	
	return true;
}

/////////////////////////////////////////////////////////////////

bool CMLFormat::DoMolWideData()
{
	//Handle atomParity and bondStereo
	vector<pair<string,string> >::iterator AttributeIter;
	for(AttributeIter=molWideData.begin();AttributeIter!=molWideData.end();++AttributeIter)
	{
		string name  = AttributeIter->first;
		string value = AttributeIter->second;

		if(name=="atomParity" || name=="bondStereo")
		{
			vector<int> AtomRefIdx;
			
			string nextname = (++AttributeIter)->first;
			string atrefsvalue = AttributeIter->second;
			if(nextname=="atomRefs4" && !atrefsvalue.empty())
			{
				vector<string> ids;
				tokenize(ids, atrefsvalue);
				int i;
				for(i=0;i<4;++i)
					AtomRefIdx.push_back(AtomMap[ids[i]]);
			}

			nextname = (++AttributeIter)->first;
			if(!(nextname=="centralAtomOrBond"))
				return false;
			
			int Idx = atoi(AttributeIter->second.c_str());
			if(name=="atomParity")
			{
				int parity =atoi(value.c_str());
				//We now have for the parity for the atom of index AtIdx
				//calculated using the atoms in AtomRefIdx.
				//Need now to adjust the parity to match the standard order
				// ...
				if(parity>0)
					_pmol->GetAtom(Idx)->SetClockwiseStereo();
				else if(parity<0)
					_pmol->GetAtom(Idx)->SetAntiClockwiseStereo();
			}
			else //bondStereo
			{
				OBBond* pbond1;
				OBBond* pbond2;
				if(atrefsvalue.empty())
				{
					OBBond* pDBond = _pmol->GetBond(Idx);
					//With no atomRefs4, the specification is either W, H,
					if(value=="W")
					{
						pDBond->SetWedge();
						return true;
					}
					else if(value=="H")
					{
						pDBond->SetHash();
						return true;
					}
					// ... or ordinary cis/trans
					//which is valid only with one substituent on each C
					
					OBAtom* pAt1 = pDBond->GetBeginAtom();
					OBAtom* pAt2 = pDBond->GetEndAtom();
					FOR_NBORS_OF_ATOM(a1,pAt1)
					  {
					    if(!a1->IsHydrogen() && &*a1!=pAt2)
					      break;
					    pbond1 = _pmol->GetBond(pAt1->GetIdx(),a1->GetIdx());
					  }
					
					FOR_NBORS_OF_ATOM(a2,pAt2)
					  {
					    if(!a2->IsHydrogen() && &*a2!=pAt1)
					      break;
					    pbond2 = _pmol->GetBond(pAt2->GetIdx(),a2->GetIdx());
					  }
				}
				else
				{
					pbond1 = _pmol->GetBond(AtomRefIdx[0],AtomRefIdx[1]);
					pbond2 = _pmol->GetBond(AtomRefIdx[2],AtomRefIdx[3]);
				}

				if(!pbond1 || !pbond2)
					return false;
				if(pbond1->IsUp() || pbond1->IsDown()) //congugated double bonds
				{
					if((pbond1->IsUp() && (value=="C")) || (pbond1->IsDown() && value=="T"))
							pbond2->SetDown();
						else
							pbond2->SetUp();
				}

				else if(pbond2->IsUp() || pbond2->IsDown()) //congugated double bonds
				{
					if((pbond2->IsUp() && (value=="C")) || (pbond2->IsDown() && value=="T"))
							pbond1->SetDown();
						else
							pbond1->SetUp();
				}
				else
				{
					pbond1->SetDown();
					if(value=="C")
						pbond2->SetUp();
					else if(value=="T")
						pbond2->SetDown();				
				}
			}
		}
	}
	//Clear here to aid embedded molecules
	AtomArray.clear();
	BondArray.clear();
	molWideData.clear();

	return true;
}

//////////////////////////////////////////////////////////
bool CMLFormat::TransferArray(cmlArray& arr)
{
	//Reads attributes of the current node, e.g. atomID="a1 a2 a3" 
	//parses each of them into their separate items, e.g. a1, a2, a3
  //and pushes them as a pairs in each of the members of the array
	// e.g. ("atomID", "a1") in AtomArray[0], ("atomID", "a2") in AtomArray[1]

	if(xmlTextReaderHasAttributes(reader()))
	{
		int ret = xmlTextReaderMoveToFirstAttribute(reader());
		while(ret==1)
		{
			const xmlChar* pname = xmlTextReaderConstName(reader());
			string name((const char*)pname);
			const xmlChar* pvalue = xmlTextReaderConstValue(reader());
			string value;
			if(pvalue)
				value = (const char*)pvalue;
			vector<string> items;
			tokenize(items,value);
			if(arr.size()<items.size())
				arr.resize(items.size());
			int i;
			for(i=0;i<items.size();++i)
			{				
				pair<string,string> nameAndvalue(name,items[i]);					
				arr[i].push_back(nameAndvalue);
			}
			ret = xmlTextReaderMoveToNextAttribute(reader());
		}
	}
	return true;
}

bool CMLFormat::TransferElement(cmlArray& arr)
{
	//Reads the attributes of the current node, e.g. <atom id="a1" elementType="C"/> 
	//pushes each of them as a pairs into each of the members of the array
	// e.g. ("id", "a1") and (elementType", "C") will be put into AtomArray[n]
	//where n is the number of times this routine has been called before.

	if(xmlTextReaderHasAttributes(reader()))
	{
		int ret = xmlTextReaderMoveToFirstAttribute(reader());
		while(ret==1)
		{
			const xmlChar* pname = xmlTextReaderConstName(reader());
			string name((const char*)pname);
			const xmlChar* pvalue = xmlTextReaderConstValue(reader());
			string value;
			if(pvalue)
				value = (const char*)pvalue;
			pair<string,string> nameAndvalue(name,value);					
			cmlBondOrAtom.push_back(nameAndvalue);			
			ret = xmlTextReaderMoveToNextAttribute(reader());
		}
	}
	return true;
}

bool CMLFormat::ParseFormula(string& formula, OBMol* pmol)
{
	vector<string> items;
	tokenize(items, formula);
	vector<string>::iterator iSymbol, iNumber;
	for(iSymbol=items.begin();iSymbol!=items.end();++iSymbol)
	{
		iNumber = iSymbol+1;
		if(iNumber==items.end())
			return false;
		int n=atoi(iNumber->c_str());
		int atno, iso=0;
		atno=etab.GetAtomicNum(iSymbol++->c_str(),iso);
		if(atno<=0 || n<=0)
			return false;
		int i;
		for(i=0;i<n;++i)
		{
			OBAtom* pAtom = pmol->NewAtom();
			pAtom->SetAtomicNum(atno);
			if(iso)
				pAtom->SetIsotope(iso);
		}
	}
	return true;
}


void CMLFormat::WriteMetadataList()
{	
	static const xmlChar C_METADATALIST[] = "metadataList";
	static const xmlChar C_METADATA[]     = "metadata";
	static const xmlChar C_TITLE[]        = "title";
	static const xmlChar C_NAME[]         = "name";
	static const xmlChar C_CONTENT[]      = "content";

	xmlTextWriterStartElement(writer(), C_METADATALIST);
	xmlTextWriterWriteAttribute(writer(), C_TITLE, BAD_CAST "generated by OpenBabel");

	xmlTextWriterStartElement(writer(), C_METADATA);
	xmlTextWriterWriteAttribute(writer(), C_NAME, BAD_CAST "dc:creator");
	string version("OpenBabel version ");
	version += BABEL_VERSION;
	xmlTextWriterWriteAttribute(writer(), C_CONTENT, BAD_CAST version.c_str());
	xmlTextWriterEndElement(writer());

	xmlTextWriterStartElement(writer(), C_METADATA);
	xmlTextWriterWriteAttribute(writer(), C_NAME, BAD_CAST "dc:description");
	xmlTextWriterWriteAttribute(writer(), C_CONTENT, BAD_CAST "Conversion of legacy filetype to CML");
	xmlTextWriterEndElement(writer());

	xmlTextWriterStartElement(writer(), C_METADATA);
	xmlTextWriterWriteAttribute(writer(), C_NAME, BAD_CAST "dc:type");
	xmlTextWriterWriteAttribute(writer(), C_CONTENT, BAD_CAST "chemistry");
	xmlTextWriterEndElement(writer());

	xmlTextWriterStartElement(writer(), C_METADATA);
	xmlTextWriterWriteAttribute(writer(), C_NAME, BAD_CAST "dc:contributor");
	xmlTextWriterWriteAttribute(writer(), C_CONTENT, BAD_CAST "unknown");
	xmlTextWriterEndElement(writer());

	xmlTextWriterStartElement(writer(), C_METADATA);
	xmlTextWriterWriteAttribute(writer(), C_NAME, BAD_CAST "dc:date");
	xmlTextWriterWriteAttribute(writer(), C_CONTENT, BAD_CAST getTimestr().c_str());
	xmlTextWriterEndElement(writer());

	xmlTextWriterStartElement(writer(), C_METADATA);
	xmlTextWriterWriteAttribute(writer(), C_NAME, BAD_CAST "cmlm:structure");
	xmlTextWriterWriteAttribute(writer(), C_CONTENT, BAD_CAST "yes");
	xmlTextWriterEndElement(writer());

	xmlTextWriterEndElement(writer());
}

string CMLFormat::getTimestr()
{
  const int TIME_STR_SIZE = 64;
	time_t akttime;                              /* Systemtime                        */
  char timestr[TIME_STR_SIZE + 1] = "";        /* Timestring                        */
  size_t time_res;                             /* Result of strftime                */

  /* ---- Get the system-time ---- */
  akttime = time((time_t *) NULL);
  time_res = strftime(timestr,
                      TIME_STR_SIZE,
                      "%a %b %d %H:%M:%S %Z %Y",
                      localtime((time_t *) &akttime)
                     );
  return timestr;
}

/////////////////////////////////////////////////////////////

bool CMLFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
	static const xmlChar C_MOLECULE[]   = "molecule";
	static const xmlChar C_CML[]        = "cml";
	static const xmlChar C_ATOMARRAY[]  = "atomArray";
	static const xmlChar C_BONDARRAY[]  = "bondArray";
	static const xmlChar C_ATOM[]       = "atom";
	static const xmlChar C_BOND[]       = "bond";
	static const xmlChar C_ID[]         = "id";
	static const xmlChar C_TITLE[]      = "title";
	static const xmlChar C_NAME[]       = "name";
	static const xmlChar C_ATOMPARITY[] = "atomParity";
//	static const xmlChar C_BONDSTEREO[] = "bondStereo";

	static const xmlChar C_X2[]               = "x2";
	static const xmlChar C_Y2[]               = "y2";
	static const xmlChar C_X3[]               = "x3";
	static const xmlChar C_Y3[]               = "y3";
	static const xmlChar C_Z3[]               = "z3";
	static const xmlChar C_XFRACT[]               = "xFract";
	static const xmlChar C_YFRACT[]               = "yFract";
	static const xmlChar C_ZFRACT[]               = "zFract";
	static const xmlChar C_ATOMID[]           = "atomID";
	static const xmlChar C_ELEMENTTYPE[]      = "elementType";
	static const xmlChar C_ISOTOPE[]          = "isotope";
	static const xmlChar C_SPINMULTIPLICITY[] = "spinMultiplicity";
	static const xmlChar C_HYDROGENCOUNT[]    = "hydrogenCount";
	static const xmlChar C_FORMALCHARGE[]     = "formalCharge";
	static const xmlChar C_ATOMREFS2[]        = "atomRefs2";
	static const xmlChar C_ATOMREF1[]         = "atomRef1";
	static const xmlChar C_ATOMREF2[]         = "atomRef2";
	static const xmlChar C_ORDER[]            = "order";
	static const xmlChar C_ATOMREFS4[]        = "atomRefs4";
/* defined in other functions
	static const xmlChar C_FORMULA[] = "formula";
	static const xmlChar C_CONCISE[] = "concise";
*/
	//CML1
	static const xmlChar C_STRING[]       = "string";
	static const xmlChar C_INTEGER[]      = "integer";
	static const xmlChar C_FLOAT[]        = "floatg";
	static const xmlChar C_BUILTIN[]      = "builtin";
	static const xmlChar C_STRINGARRAY[]  = "stringArray";
	static const xmlChar C_INTEGERARRAY[] = "integerArray";
	static const xmlChar C_FLOATARRAY[]   = "floatArray";
/* used as ordinary text
	atomRef  
*/

	const xmlChar* C_X3orFRACT = C_X3; //Non-fraction coordinates are the default
	const xmlChar* C_Y3orFRACT = C_Y3;
	const xmlChar* C_Z3orFRACT = C_Z3;

	_pxmlConv = XMLConversion::GetDerived(pConv,false);
	if(!_pxmlConv)
		return false;

	OBMol* pmol = dynamic_cast<OBMol*>(pOb);
	if(pmol==NULL)
			return false;
	OBMol &mol = *pmol;

	int numbonds = mol.NumBonds(); //Capture this before deleting Hs
	bool UseHydrogenCount=false;
	if(_pxmlConv->IsOption("h"))
	{
		pmol->DeleteHydrogens();
		UseHydrogenCount=true;
	}
	
	bool UseFormulaWithNoBonds=true;

	bool cml1 = _pxmlConv->IsOption("1");
	bool arrayform = _pxmlConv->IsOption("a");
	int dim = mol.GetDimension();
	
	prefix = BAD_CAST _pxmlConv->IsOption("N");
	
	xmlChar* uri=NULL;
	
	if(!_pxmlConv->IsOption("MolsNotStandalone") && _pxmlConv->GetOutputIndex()==1)
	{
		if(!_pxmlConv->IsOption("x"))
		{
			xmlTextWriterStartDocument(writer(), NULL, NULL, NULL);
			if(cml1)
				uri = BAD_CAST CML1NamespaceURI();
			else
				uri=BAD_CAST NamespaceURI();
		}
		//If more than one molecule to be output, write <cml> at start and </cml> at end.
		//Except if Option "0" set, e.g. by CMLReactFormat
		if(!_pxmlConv->IsLast())
		{
			xmlTextWriterStartElementNS(writer(), prefix, C_CML, uri);
			uri=NULL;;
		}
	}
	
	xmlTextWriterStartElementNS(writer(), prefix, C_MOLECULE, uri);
	
	const char* id = mol.GetTitle();
	if(*id)
		xmlTextWriterWriteAttribute(writer(), C_ID, BAD_CAST id);

	if(_pxmlConv->IsOption("m") && _pxmlConv->GetOutputIndex()==1) //only on first molecule
		WriteMetadataList();

	pUnitCell = NULL;
	if (!cml1 && mol.HasData(OBGenericDataType::UnitCell))
	{
		WriteCrystal(mol);//Output will be in crystallographic form
		UseFormulaWithNoBonds = false;
	}

	if(mol.NumAtoms()>0)
	{
		if(numbonds==0 && UseFormulaWithNoBonds)
			WriteFormula(mol);
		else
		{
			xmlTextWriterStartElementNS(writer(), prefix, C_ATOMARRAY, NULL);
			#ifdef HAVE_SSTREAM
					stringstream id, eltyp, iso, chg, spn, hct, x, y, z;
			#else
					strstream id, eltyp, iso, chg, spn, hct, x, y, z;
			#endif
			bool anyChg=false, anySpin=false, anyIsotope=false;
			double X, Y, Z; //atom coordinates

			OBAtom *patom;
			vector<OBNodeBase*>::iterator i;
			for (patom = mol.BeginAtom(i);patom;patom = mol.NextAtom(i))
			{
				string el(etab.GetSymbol(patom->GetAtomicNum()));
				if(el=="Xx")
					el="R";

				int charge = patom->GetFormalCharge();
				int spin = patom->GetSpinMultiplicity();
				int isotope =patom->GetIsotope();

				int hcount=patom->ImplicitHydrogenCount();

				X = patom->GetX();
				Y = patom->GetY();
				Z = patom->GetZ();
				if(pUnitCell)
				{
					//Convert to fractional coordinates
					vector3 v = patom->GetVector();
					v *= pUnitCell->GetFractionalMatrix();
					X = v.x();
					Y = v.y();
					Z = v.z();
					C_X3orFRACT = C_XFRACT;
					C_Y3orFRACT = C_YFRACT;
					C_Z3orFRACT = C_ZFRACT;
					dim=3; //should already be, but make sure
				}

				if(arrayform)
				{
					if(charge)
						anyChg = true;
					if(spin)
						anySpin = true;
					if(isotope)
						anyIsotope = true;
					id << " " << "a" << patom->GetIdx();
					eltyp << " " << el;
					iso << " " << isotope;
					chg << " " << charge;
					spn << " " << spin;
					hct << " " << hcount;

					x << " " << X;
					y << " " << Y;
					z << " " << Z;
				}
				else
				{
					//Non-array form
					xmlTextWriterStartElementNS(writer(), prefix, C_ATOM, NULL);
					xmlTextWriterWriteFormatAttribute(writer(), C_ID,"a%d", patom->GetIdx());
					
					if(!cml1)
					{	
						xmlTextWriterWriteFormatAttribute(writer(), C_ELEMENTTYPE,"%s", el.c_str());
						if(isotope)
							xmlTextWriterWriteFormatAttribute(writer(), C_ISOTOPE,"%d", isotope);

						if(charge)
							xmlTextWriterWriteFormatAttribute(writer(), C_FORMALCHARGE,"%d", charge);

						if(spin)
							xmlTextWriterWriteFormatAttribute(writer(), C_SPINMULTIPLICITY,"%d", spin);

						if(UseHydrogenCount && hcount)
							xmlTextWriterWriteFormatAttribute(writer(), C_HYDROGENCOUNT,"%d", hcount);

						if(dim==2)
						{
							xmlTextWriterWriteFormatAttribute(writer(), C_X2,"%f", X);
							xmlTextWriterWriteFormatAttribute(writer(), C_Y2,"%f", Y);
						}
						if(dim==3)
						{
							xmlTextWriterWriteFormatAttribute(writer(), C_X3orFRACT,"%f", X);
							xmlTextWriterWriteFormatAttribute(writer(), C_Y3orFRACT,"%f", Y);
							xmlTextWriterWriteFormatAttribute(writer(), C_Z3orFRACT,"%f", Z);
						}
						if(patom->HasChiralitySpecified())
						{
							int cfg=0;
							if((patom->IsPositiveStereo() || patom->IsClockwise()))
								cfg=1; //whether +1 or -1 is pure guess. TODO***
							else if(patom->IsNegativeStereo() || patom->IsAntiClockwise())
								cfg=-1;
							if(cfg)
							{
								//in the order they are in OBMol except that any H is put at the end
								vector<int> ref;
								ref.clear();
								int Hidx=0;
								FOR_NBORS_OF_ATOM(a,patom)
								{
									if(a->IsHydrogen())
										Hidx = a->GetIdx();
									else
										ref.push_back(a->GetIdx()); 
								}
								if(Hidx)
									ref.push_back(Hidx);
								if(ref.size()==3)//e.g. using implicit Hs
									ref.push_back(patom->GetIdx());

								xmlTextWriterStartElementNS(writer(), prefix, C_ATOMPARITY, NULL);
								xmlTextWriterWriteFormatAttribute(writer(), C_ATOMREFS4,
									"a%d a%d a%d a%d", ref[0], ref[1], ref[2], ref[3]);														
								xmlTextWriterWriteFormatString(writer(),"%d", cfg);
								xmlTextWriterEndElement(writer());//atomParity
							}
						}
					}
					else
					{
						//CML1
						xmlTextWriterStartElementNS(writer(), prefix, C_STRING, NULL);
						xmlTextWriterWriteFormatAttribute(writer(), C_BUILTIN,"%s", "elementType");
						xmlTextWriterWriteFormatString(writer(),"%s", el.c_str());
						xmlTextWriterEndElement(writer());

						if(charge)
						{
							xmlTextWriterStartElementNS(writer(), prefix, C_INTEGER, NULL);
							xmlTextWriterWriteFormatAttribute(writer(), C_BUILTIN,"%s", "formalCharge");
							xmlTextWriterWriteFormatString(writer(),"%d", charge);
							xmlTextWriterEndElement(writer());
						}

						if(UseHydrogenCount && hcount)
						{
							xmlTextWriterStartElementNS(writer(), prefix, C_INTEGER, NULL);
							xmlTextWriterWriteFormatAttribute(writer(), C_BUILTIN,"%s", "hydrogenCount");
							xmlTextWriterWriteFormatString(writer(),"%d", hcount);
							xmlTextWriterEndElement(writer());
						}

						if(dim==2 || dim==3)
						{
							xmlTextWriterStartElementNS(writer(), prefix, C_FLOAT, NULL);
							xmlTextWriterWriteFormatAttribute(writer(), C_BUILTIN,"%s%d", "x",dim);
							xmlTextWriterWriteFormatString(writer(),"%f", X);
							xmlTextWriterEndElement(writer());

							xmlTextWriterStartElementNS(writer(), prefix, C_FLOAT, NULL);
							xmlTextWriterWriteFormatAttribute(writer(), C_BUILTIN,"%s%d", "y",dim);
							xmlTextWriterWriteFormatString(writer(),"%f", Y);
							xmlTextWriterEndElement(writer());
						}

						if(dim==3)
						{
							xmlTextWriterStartElementNS(writer(), prefix, C_FLOAT, NULL);
							xmlTextWriterWriteFormatAttribute(writer(), C_BUILTIN,"%s%d", "z",dim);
							xmlTextWriterWriteFormatString(writer(),"%f", Z);
							xmlTextWriterEndElement(writer());
						}
						//Stereochemistry currently not written for CML1
					}
					xmlTextWriterEndElement(writer());//atom
				}
			}
		
			if(arrayform)
			{
				if(!cml1)
				{
					xmlTextWriterWriteFormatAttribute(writer(), C_ATOMID,"%s", id.str().c_str());
					xmlTextWriterWriteFormatAttribute(writer(), C_ELEMENTTYPE,"%s", eltyp.str().c_str());
					
					if(anyIsotope)
						xmlTextWriterWriteFormatAttribute(writer(), C_ISOTOPE,"%s", iso.str().c_str());

					if(anyChg)
						xmlTextWriterWriteFormatAttribute(writer(), C_FORMALCHARGE,"%s", chg.str().c_str());

					if(anySpin)
						xmlTextWriterWriteFormatAttribute(writer(), C_SPINMULTIPLICITY,"%s", spn.str().c_str());

					if(UseHydrogenCount)
						xmlTextWriterWriteFormatAttribute(writer(), C_HYDROGENCOUNT,"%s", hct.str().c_str());

					if(dim==2)
					{
						xmlTextWriterWriteFormatAttribute(writer(), C_X2,"%s", x.str().c_str());
						xmlTextWriterWriteFormatAttribute(writer(), C_Y2,"%s", y.str().c_str());
					}
					if(dim==3)
					{
						xmlTextWriterWriteFormatAttribute(writer(), C_X3orFRACT,"%s", x.str().c_str());
						xmlTextWriterWriteFormatAttribute(writer(), C_Y3orFRACT,"%s", y.str().c_str());
						xmlTextWriterWriteFormatAttribute(writer(), C_Z3orFRACT,"%s", z.str().c_str());
					}
				}
				else
				{
					//CML1
					xmlTextWriterStartElementNS(writer(), prefix, C_STRINGARRAY, NULL);
					xmlTextWriterWriteFormatAttribute(writer(), C_BUILTIN,"%s", "atomID");
					xmlTextWriterWriteFormatString(writer(),"%s", id.str().c_str());
					xmlTextWriterEndElement(writer());

					xmlTextWriterStartElementNS(writer(), prefix, C_STRINGARRAY, NULL);
					xmlTextWriterWriteFormatAttribute(writer(), C_BUILTIN,"%s", "elementType");
					xmlTextWriterWriteFormatString(writer(),"%s", eltyp.str().c_str());
					xmlTextWriterEndElement(writer());

					if(anyChg)
					{
						xmlTextWriterStartElementNS(writer(), prefix, C_INTEGERARRAY, NULL);
						xmlTextWriterWriteFormatAttribute(writer(), C_BUILTIN,"%s", "formalCharge");
						xmlTextWriterWriteFormatString(writer(),"%s", chg.str().c_str());
						xmlTextWriterEndElement(writer());
					}

					if(UseHydrogenCount)
					{
						xmlTextWriterStartElementNS(writer(), prefix, C_INTEGERARRAY, NULL);
						xmlTextWriterWriteFormatAttribute(writer(), C_BUILTIN,"%s", "hydrogenCount");
						xmlTextWriterWriteFormatString(writer(),"%s", hct.str().c_str());
						xmlTextWriterEndElement(writer());
					}

					if(dim==2 || dim==3)
					{
						xmlTextWriterStartElementNS(writer(), prefix, C_FLOATARRAY, NULL);
						xmlTextWriterWriteFormatAttribute(writer(), C_BUILTIN,"%s%d", "x",dim);
						xmlTextWriterWriteFormatString(writer(),"%s", x.str().c_str());
						xmlTextWriterEndElement(writer());

						xmlTextWriterStartElementNS(writer(), prefix, C_FLOATARRAY, NULL);
						xmlTextWriterWriteFormatAttribute(writer(), C_BUILTIN,"%s%d", "y",dim);
						xmlTextWriterWriteFormatString(writer(),"%s", y.str().c_str());
						xmlTextWriterEndElement(writer());
					}
					if(dim==3)
					{
						xmlTextWriterStartElementNS(writer(), prefix, C_FLOATARRAY, NULL);
						xmlTextWriterWriteFormatAttribute(writer(), C_BUILTIN,"%s%d", "z",dim);
						xmlTextWriterWriteFormatString(writer(),"%s", z.str().c_str());
						xmlTextWriterEndElement(writer());
					}			
				}
			}
			xmlTextWriterEndElement(writer());//atomArray
		}
	}

	if(mol.NumBonds()>0)
	{
		xmlTextWriterStartElementNS(writer(), prefix, C_BONDARRAY, NULL);

		#ifdef HAVE_SSTREAM
			stringstream ref1, ref2, ord;
		#else
			strstream ref1, ref2, ord;
		#endif
		OBBond *pbond;
		vector<OBEdgeBase*>::iterator ib;
		for (pbond = mol.BeginBond(ib);pbond;pbond = mol.NextBond(ib))
		{
			int bo = pbond->GetBO();
			if(arrayform)
			{
				ref1 << " a" << pbond->GetBeginAtomIdx();
				ref2 << " a" << pbond->GetEndAtomIdx();
				if(bo==5) //aromatic
					ord << " " << 'A';
				else
					ord << " " << bo;
			}
			else
			{
				xmlTextWriterStartElementNS(writer(), prefix, C_BOND, NULL);
//				xmlTextWriterWriteFormatAttribute(writer(), C_ID,"b%d", pbond->GetIdx()); remove bond id
				if(!cml1)
				{
					xmlTextWriterWriteFormatAttribute(writer(), C_ATOMREFS2,"a%d a%d",
						pbond->GetBeginAtomIdx(), pbond->GetEndAtomIdx() );
					if(bo==5) //aromatic
						xmlTextWriterWriteFormatAttribute(writer(), C_ORDER,"%c", 'A');
					else
						xmlTextWriterWriteFormatAttribute(writer(), C_ORDER,"%d", bo);
					
					if(bo==2)
						WriteBondStereo(pbond);
				}
				else
				{
					//CML1
					xmlTextWriterStartElementNS(writer(), prefix, C_STRING, NULL);
					xmlTextWriterWriteFormatAttribute(writer(), C_BUILTIN,"%s", "atomRef");
					xmlTextWriterWriteFormatString(writer(),"a%d", pbond->GetBeginAtomIdx());
					xmlTextWriterEndElement(writer());

					xmlTextWriterStartElementNS(writer(), prefix, C_STRING, NULL);
					xmlTextWriterWriteFormatAttribute(writer(), C_BUILTIN,"%s", "atomRef");
					xmlTextWriterWriteFormatString(writer(),"a%d", pbond->GetEndAtomIdx());
					xmlTextWriterEndElement(writer());

					xmlTextWriterStartElementNS(writer(), prefix, C_STRING, NULL);
					xmlTextWriterWriteFormatAttribute(writer(), C_BUILTIN,"%s", "order");
					xmlTextWriterWriteFormatString(writer(),"%d", bo);
					xmlTextWriterEndElement(writer());
				}
				xmlTextWriterEndElement(writer());//bond
			}
		}
		if(arrayform)
		{
			if(!cml1)
			{
				xmlTextWriterWriteFormatAttribute(writer(), C_ATOMREF1, "%s", ref1.str().c_str());
				xmlTextWriterWriteFormatAttribute(writer(), C_ATOMREF2, "%s", ref2.str().c_str());
				xmlTextWriterWriteFormatAttribute(writer(), C_ORDER, "%s", ord.str().c_str());
			}
			else
			{
				//CML1
				xmlTextWriterStartElementNS(writer(), prefix, C_STRINGARRAY, NULL);
				xmlTextWriterWriteFormatAttribute(writer(), C_BUILTIN,"%s", "atomRef");
				xmlTextWriterWriteFormatString(writer(),"%s", ref1.str().c_str());
				xmlTextWriterEndElement(writer());

				xmlTextWriterStartElementNS(writer(), prefix, C_STRINGARRAY, NULL);
				xmlTextWriterWriteFormatAttribute(writer(), C_BUILTIN,"%s", "atomRef");
				xmlTextWriterWriteFormatString(writer(),"%s", ref2.str().c_str());
				xmlTextWriterEndElement(writer());

				xmlTextWriterStartElementNS(writer(), prefix, C_STRINGARRAY, NULL);
				xmlTextWriterWriteFormatAttribute(writer(), C_BUILTIN,"%s", "order");
				xmlTextWriterWriteFormatString(writer(),"%s", ord.str().c_str());
				xmlTextWriterEndElement(writer());
			}
		}

		xmlTextWriterEndElement(writer());//bondArray
	}

	xmlTextWriterEndElement(writer());//molecule

	if(!_pxmlConv->IsOption("MolsNotStandalone") && _pxmlConv->IsLast())
	{
		xmlTextWriterEndDocument(writer());
		OutputToStream();
	}
	return true;
}

void CMLFormat::WriteFormula(OBMol& mol)
{
	string formula = mol.GetFormula();
	//Insert spaces and add 1s if missing
	string modformula;
	char lastch;
	int i;
	for(i=0; i<formula.size();++i)
	{
		char ch = formula[i];
		if(i>0 && isupper(ch) && !isdigit(lastch))
			modformula += " 1 ";
		else			
		{
			if( (isdigit(ch) && !isdigit(lastch)) || (!isdigit(ch) && isdigit(lastch)) )
				modformula += ' ';
		}
		modformula += ch;
		lastch = ch;
	}
	if(!isdigit(lastch))
		modformula += " 1";

	static const xmlChar C_FORMULA[] = "formula";
	static const xmlChar C_CONCISE[] = "concise";
	xmlTextWriterStartElementNS(writer(), prefix, C_FORMULA, NULL);
	xmlTextWriterWriteFormatAttribute(writer(), C_CONCISE,"%s", modformula.c_str());
	xmlTextWriterEndElement(writer());//formula
}

void CMLFormat::WriteBondStereo(OBBond* pbond)
{
	static const xmlChar C_ATOMREFS4[]  = "atomRefs4";
	static const xmlChar C_BONDSTEREO[] = "bondStereo";

	int ud1=0, ud2=0;
	int idx1, idx2;
	OBAtom* patomA = pbond->GetBeginAtom();
	FOR_BONDS_OF_ATOM(b1,patomA)
	{
		if(b1->IsUp() || b1->IsDown() )
		{
			idx1=(b1->GetNbrAtom(patomA))->GetIdx();
			ud1 = b1->IsDown() ? -1 : 1;
			break;
		}
	}
	OBAtom* patomB = pbond->GetEndAtom();
	FOR_BONDS_OF_ATOM(b2,patomB)
	{
		if(b2->IsUp() || b2->IsDown() )
		{
			idx2=(b2->GetNbrAtom(patomB))->GetIdx();
			ud2 = b2->IsDown() ? -1 : 1;
			break;
		}
	}
	if(!ud1 || !ud2)
		return;
	
	xmlTextWriterStartElementNS(writer(), prefix, C_BONDSTEREO, NULL);
	xmlTextWriterWriteFormatAttribute(writer(), C_ATOMREFS4,
		"a%d a%d a%d a%d", idx1, patomA->GetIdx(), patomB->GetIdx(), idx2);
	char ch = (ud1==ud2) ? 'T' : 'C';
	xmlTextWriterWriteFormatString(writer(),"%c", ch);
	xmlTextWriterEndElement(writer());//bondStereo
}

void CMLFormat::WriteCrystal(OBMol& mol)
{
	static const xmlChar C_CRYSTAL[]  = "crystal";
	static const xmlChar C_SCALAR[] = "scalar";
	static const xmlChar C_Z[] = "z";
	static const xmlChar C_TITLE[] = "title";
	static const xmlChar C_UNITS[] = "units";

	pUnitCell = (OBUnitCell*)mol.GetData(OBGenericDataType::UnitCell);

	xmlTextWriterStartElementNS(writer(), prefix, C_CRYSTAL, NULL);
//	xmlTextWriterWriteFormatAttribute(writer(), C_z,"%d", number of molecules per cell);

	xmlTextWriterStartElementNS(writer(), prefix, C_SCALAR, NULL);
	xmlTextWriterWriteFormatAttribute(writer(), C_TITLE,"%s", "a");
	xmlTextWriterWriteFormatAttribute(writer(), C_UNITS,"%s", "units:angstrom");
	xmlTextWriterWriteFormatString(writer(),"%f", pUnitCell->GetA());
	xmlTextWriterEndElement(writer());//scalar

	xmlTextWriterStartElementNS(writer(), prefix, C_SCALAR, NULL);
	xmlTextWriterWriteFormatAttribute(writer(), C_TITLE,"%s", "b");
	xmlTextWriterWriteFormatAttribute(writer(), C_UNITS,"%s", "units:angstrom");
	xmlTextWriterWriteFormatString(writer(),"%f", pUnitCell->GetB());
	xmlTextWriterEndElement(writer());//scalar

	xmlTextWriterStartElementNS(writer(), prefix, C_SCALAR, NULL);
	xmlTextWriterWriteFormatAttribute(writer(), C_TITLE,"%s", "c");
	xmlTextWriterWriteFormatAttribute(writer(), C_UNITS,"%s", "units:angstrom");
	xmlTextWriterWriteFormatString(writer(),"%f", pUnitCell->GetC());
	xmlTextWriterEndElement(writer());//scalar

	xmlTextWriterStartElementNS(writer(), prefix, C_SCALAR, NULL);
	xmlTextWriterWriteFormatAttribute(writer(), C_TITLE,"%s", "alpha");
	xmlTextWriterWriteFormatAttribute(writer(), C_UNITS,"%s", "units:degree");
	xmlTextWriterWriteFormatString(writer(),"%f", pUnitCell->GetAlpha());
	xmlTextWriterEndElement(writer());//scalar

	xmlTextWriterStartElementNS(writer(), prefix, C_SCALAR, NULL);
	xmlTextWriterWriteFormatAttribute(writer(), C_TITLE,"%s", "beta");
	xmlTextWriterWriteFormatAttribute(writer(), C_UNITS,"%s", "units:degree");
	xmlTextWriterWriteFormatString(writer(),"%f", pUnitCell->GetBeta());
	xmlTextWriterEndElement(writer());//scalar

	xmlTextWriterStartElementNS(writer(), prefix, C_SCALAR, NULL);
	xmlTextWriterWriteFormatAttribute(writer(), C_TITLE,"%s", "gamma");
	xmlTextWriterWriteFormatAttribute(writer(), C_UNITS,"%s", "units:degree");
	xmlTextWriterWriteFormatString(writer(),"%f", pUnitCell->GetGamma());
	xmlTextWriterEndElement(writer());//scalar

	xmlTextWriterEndElement(writer());//crystal	
}

}//namespace

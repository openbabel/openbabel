/**********************************************************************
Copyright (C) 2005 by Chris Morley

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
#include <openbabel/babelconfig.h>
#include <openbabel/xml.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obiter.h>
#include <openbabel/elements.h>

#ifdef WIN32
#pragma warning (disable : 4800)
#endif

using namespace std;
namespace OpenBabel
{

class PubChemFormat : public XMLMoleculeFormat
{
public:
	PubChemFormat()
	{
    OBConversion::RegisterFormat("pc", this, "chemical/x-ncbi-asn1-xml");
		XMLConversion::RegisterXMLFormat(this);
	}
	virtual const char* NamespaceURI()const{return "http://www.ncbi.nlm.nih.gov";}
  virtual const char* Description()
  {
    return
      "PubChem format\n"
      "An XML format containing information on PubChem entries.\n"
      "`PubChem <http://pubchem.ncbi.nlm.nih.gov/>`_ is a freely-available\n"
      "database of chemical compounds and their properties.\n\n"
      "OpenBabel only extracts the chemical structure information, and the\n"
      "potentially large amount of other information is currently ignored.\n"
      "The format seems to handle multiple conformers, but only one is read\n"
      "(this needs testing).\n\n";
  }

  virtual const char* SpecificationURL()
  {return "ftp://ftp.ncbi.nlm.nih.gov/pubchem/data_spec/pubchem.xsd";};

  virtual const char* GetMIMEType()
  { return "chemical/x-ncbi-asn1-xml"; };

  virtual unsigned int Flags()
  {
    return (READXML | NOTWRITABLE);
  };

	virtual bool DoElement(const string& name);
	virtual bool EndElement(const string& name);

	// EndTag is used so that the stream buffer is is filled with the XML from
	// complete objects, as far as possible.
	virtual const char* EndTag(){ return "/PC-Compound>"; };

private:
	int dim;
	vector<int> AtNum;

	vector<int> BondBeginAtIndx;
	vector<int> BondEndAtIndx;
	vector<int> BondOrder;

	vector<int> CoordIndx;
	int ConformerIndx;
	vector<double> Coordx;
	vector<double> Coordy;
	vector<double> Coordz;
};

////////////////////////////////////////////////////////////////////

PubChemFormat thePubChemFormat;

////////////////////////////////////////////////////////////////////

bool PubChemFormat::DoElement(const string& name)
{
	if(name=="PC-Compound")
	{
		//This is the start of the molecule we are extracting and it will
		//be put into the OBMol* _pmol declared in the parent class.
		//initialise everything
		dim=0;
		AtNum.clear();
		BondBeginAtIndx.clear();
		BondEndAtIndx.clear();
		BondOrder.clear();
		CoordIndx.clear();
		ConformerIndx=0;
		Coordx.clear();
		Coordy.clear();
		Coordz.clear();
		_pmol->BeginModify();
	}
	if(name=="PC-Element")
	{
		int AtNumber;
		if(!_pxmlConv->GetContentInt(AtNumber) || !AtNumber)
			return false;
		AtNum.push_back(AtNumber);
	}

	if(name=="PC-CompoundType_id_cid")
	{
		_pmol->SetTitle(_pxmlConv->GetContent().c_str());
	}
	else if(name=="PC-Bonds_aid1_E")
	{
		int indx;
		if(_pxmlConv->GetContentInt(indx))
			BondBeginAtIndx.push_back(indx);
	}
	else if(name=="PC-Bonds_aid2_E")
	{
		int indx;
		if(_pxmlConv->GetContentInt(indx))
			BondEndAtIndx.push_back(indx);
	}
	else if(name=="PC-BondType")
	{
		int order;
		if(_pxmlConv->GetContentInt(order))
			BondOrder.push_back(order);
	}
	else if(name=="PC-CoordinateType")
	{
		if(_pxmlConv->GetAttribute("value")=="twod")
			dim=2;
		else if(_pxmlConv->GetAttribute("value")=="threed")
			dim=3;
		_pmol->SetDimension(dim);
	}

	else if(name=="PC-Coordinates_aid_E")
	{
		int indx;
		if(_pxmlConv->GetContentInt(indx))
			CoordIndx.push_back(indx);
	}
	else if(name=="PC-Conformer_x_E")
	{
		if(ConformerIndx)
			return true; //currently only one conformer is read
		double x;
		if(_pxmlConv->GetContentDouble(x))
			Coordx.push_back(x);
	}
	else if(name=="PC-Conformer_y_E")
	{
		if(ConformerIndx)
			return true; //currently only one conformer is read
		double y;
		if(_pxmlConv->GetContentDouble(y))
			Coordy.push_back(y);
	}
	else if(name=="PC-Conformer_z_E")
	{
		if(ConformerIndx)
			return true; //currently only one conformer is read
		double z;
		if(_pxmlConv->GetContentDouble(z))
			Coordz.push_back(z);
	}
	return true;
}

bool PubChemFormat::EndElement(const string& name)
{
  unsigned int i;
	if(name=="PC-Atoms")
	{
		for(i=0;i<AtNum.size();++i)
		{
			OBAtom* pAtom = _pmol->NewAtom();
			pAtom->SetAtomicNum(AtNum[i]);
		}
	}
	else if(name=="PC-Bonds")
	{
		for(i=0;i<BondBeginAtIndx.size();++i)
			_pmol->AddBond(BondBeginAtIndx[i],BondEndAtIndx[i],BondOrder[i]);
	}
	else if(name=="PC-Conformer")
	{
		++ConformerIndx;
		if(Coordz.size()!=Coordx.size())
			Coordz.resize(Coordx.size());
		for(i=0;i<CoordIndx.size();++i)
		{
			OBAtom* pAtom = _pmol->GetAtom(CoordIndx[i]);
			pAtom->SetVector(Coordx[i],Coordy[i],Coordz[i]);
		}
	}
	else if(name=="PC-Compound") //this is the end of the molecule we are extracting
	{
		_pmol->EndModify();
		return false;//means stop parsing
	}
	return true;
}


}//namespace

//cmlpp version
/**********************************************************************
Copyright (C) 2002-2003 Peter Murray-Rust.
Some portions Copyright (C) 2003-2005 by Geoffrey R. Hutchison
Some Portions Copyright (C) 2004 by Chris Morley
 
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
#include "obmolecformat.h"

// cml includes
#include "cml.hpp"

#include "CMLDocument.hpp"
#include "CMLMolecule.hpp"
#include "CMLAtom.hpp"
#include "CMLBond.hpp"
#include "CMLAtomArray.hpp"
#include "CMLBondArray.hpp"
#include "InputStream.hpp"
#include "InputFile.hpp"

#ifdef WIN32
#pragma warning (disable : 4800)
#endif

using namespace std;
using namespace CML;
namespace OpenBabel
{
//The original routines
extern bool ReadCML(istream &ifs,OBMol &mol, const char *title);
extern bool WriteCML(ostream &ofs,OBMol &mol,const char *dim,const char* xmlOptions);

class CMLFormat : public OBMoleculeFormat
{
public:
	//Register this format type ID
	CMLFormat()
	{
			OBConversion::RegisterFormat("cml", this, "chemical/x-cml");
	}

	virtual const char* Description()
	{
			return " \
Chemical Markup Language\n \
XML format. This implementation uses cmlpp.\n \
Write options for CML: -x[flags] (e.g. -x1ac)\n \
1  output CML V1.0  or \n \
2  output CML V2.0 (default)\n \
a  output array format for atoms and bonds\n \
n<prefix> add namespace prefix to elements\n \
c  use 'cml' as output namespace prefix \n \
x  omit XML declaration\n \
\n";
};

  virtual const char* SpecificationURL()
  {return "http://wwmm.ch.cam.ac.uk/moin/ChemicalMarkupLanguage";};

  virtual const char* GetMIMEType() 
  { return "chemical/x-cml"; };

  ////////////////////////////////////////////////////
  /// The "API" interface functions
  virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
  virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

  ////////////////////////////////////////////////////
  /// The "Convert" interface functions
  virtual bool ReadChemObject(OBConversion* pConv);

};

//Make an instance of the format class
CMLFormat theCMLFormat;

///////////////////////////////////////////////////
bool CMLFormat::ReadChemObject(OBConversion* pConv)
{
	InputStream input(*pConv->GetInStream());

	CMLDocument* documentPtr = new CMLDocument;
	documentPtr -> createDefaultSaxHandler ();

	SaxHandler *sax_p = documentPtr -> getSaxHandler ();

	if (!sax_p)
		return false;
	
	documentPtr->readXML(input);

	Element* pRoot = documentPtr->getRootElementPtr();
  if (!pRoot)
		return false;

	NodeList* pNodes  = pRoot -> getElementsByTagName (C_MOLECULE);
	if(!pNodes)
		return false;
	int iMol;
	for(iMol=0;iMol<pNodes->getLength();++iMol)
	{
		CMLMolecule* pCMLmol = dynamic_cast<CMLMolecule*>(pNodes->item(iMol)); //(*itr);	
		
		//Ignore <molecule ref="xxx"/> elements
		if(!pCMLmol->getAttributeValue("ref").empty())
			continue;

		OBMol* pmol = new OBMol;
		pmol->BeginModify();
		
		pmol->SetTitle(pCMLmol->getAttributeValue(C_TITLE));
		
		std::map<string,int> AtomMap; //key=atom id, value= ob atom index 
		CMLAtomArray* pAtArray = dynamic_cast<CMLAtomArray*>(pCMLmol->getFirstChildWithName(C_ATOMARRAY));
		if(pAtArray)
		{
			int iAtom;
			for(iAtom=0;iAtom<pAtArray->getAtomCount();++iAtom)
			{
				CMLAtom* pCMLAtom = pAtArray->getAtom(iAtom);
				OBAtom obatom;
				
				int atno, iso=0;
				atno=etab.GetAtomicNum(pCMLAtom->getElementType().c_str(),iso);
				obatom.SetAtomicNum(atno);
				if(iso)
					obatom.SetIsotope(iso);
				
				if(pCMLAtom->formalChargeIsSet())
					obatom.SetFormalCharge(pCMLAtom->getFormalCharge());
				if(pCMLAtom->spinMultiplicityIsSet())
					obatom.SetSpinMultiplicity(pCMLAtom->getSpinMultiplicity());
				if(pCMLAtom->isotopeNumberIsSet())				
					obatom.SetIsotope(pCMLAtom->getIsotopeNumber());
				
				double x=0,y=0,z=0;
				int dim=0;
				if(pCMLAtom->x2IsSet())
					dim=2;
				else if(pCMLAtom->x3IsSet())
					dim=3;
				pmol->SetDimension(dim);
				if(dim==2)
				{
					x = pCMLAtom->getX2();
					y = pCMLAtom->getY2();
				}
				if(dim==3)
				{
					x = pCMLAtom->getX3();
					y = pCMLAtom->getY3();
					z = pCMLAtom->getZ3();
				}
				obatom.SetVector(x, y , z);

				pmol->AddAtom(obatom);
				AtomMap[pCMLAtom->getId()] = pmol->NumAtoms(); //GetIdx() return 0 for some reason

				if(pCMLAtom->hydrogenCountIsSet())
				{
					int nprev=pmol->NumAtoms();
					int i;
					for(i=0;i<pCMLAtom->getHydrogenCount();++i)
					{
						OBAtom* hatom = pmol->NewAtom();
						hatom->SetAtomicNum(1);
						hatom->SetType("H");
						pmol->AddBond(nprev,pmol->NumAtoms(),1);
					}
				}

				obatom.Clear();			
			}
		}
		
		CMLBondArray* pBArray = dynamic_cast<CMLBondArray*>(pCMLmol->getFirstChildWithName(C_BONDARRAY));
		if(pBArray)
		{
			int iBond;
			for(iBond=0;iBond<pBArray->getBondCount();++iBond)
			{
				CMLBond* pCMLBond = pBArray->getBond(iBond);
				OBBond obbond;

				string atid1 = pCMLBond->getAtomRefs2().first;
				string atid2 = pCMLBond->getAtomRefs2().second;
				if(atid1.empty() || atid2.empty())
					return false;
				
				std::map<string,int>::iterator pos1, pos2;				
				pos1 = AtomMap.find(atid1);
				pos2 = AtomMap.find(atid2);
				if(pos1==AtomMap.end() || pos1==AtomMap.end())
				{
					cerr << "Incorrect atom references" << endl;
					return false;
				}
				OBAtom* at1, *at2;
				at1=pmol->GetAtom(pos1->second);
				at2=pmol->GetAtom(pos2->second);
				
				int ord=0;
				const char bo = pCMLBond->getOrder()[0];
				if(bo=='S')
					ord=1;
				else if(bo=='D')
					ord=2;
				else if(bo=='A')
					ord=5;
				else
					ord=atoi(&bo);
				if(!ord)
				{
					cerr << "Incorrect bond order" << endl;
					return false;
				}

				obbond.Set(iBond,at1,at2,ord,0);
				pmol->AddBond(obbond);

				//pCMLBond->getBondStereo();
			}
		}
		pmol->EndModify();

		if(!pConv->AddChemObject(pmol->DoTransformations(pConv->GetOptions(OBConversion::GENOPTIONS))))
			return false;
	}
	delete pNodes;
	return true;
}

////////////////////////////////////////////////////////////////////
bool CMLFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
{
	ReadChemObject(pConv); //This makes a new OBMol
	OBMol* pDest = dynamic_cast<OBMol*> (pOb); 
	OBMol* pSrc  = dynamic_cast<OBMol*> (pConv->GetChemObject());
	if((pDest==NULL) || (pSrc==NULL)) return false;

	*pDest = *pSrc; //Copy the read molecule to the desired place
	delete pSrc; //delete the new OBMol
	return true;
}

////////////////////////////////////////////////////////////////
bool CMLFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
	string auditMsg = "OpenBabel::Write molecule ";
	string description(Description());
	auditMsg += description.substr(0,description.find('\n'));
	obErrorLog.ThrowError(__FUNCTION__, auditMsg, obAuditMsg);
	
	OBMol* pmol = dynamic_cast<OBMol*>(pOb);
	if(pmol==NULL)
			return false;

	ostream &ofs = *pConv->GetOutStream();
	OBMol &mol = *pmol;
	int dim = mol.GetDimension();

	bool arrayform = pConv->IsOption("a");
	string ns;
	if(pConv->IsOption("c"))
		ns = "cml:";
	const char* nstxt = pConv->IsOption("n");
	if(nstxt)
	{
		ns = nstxt;
		ns += ":";
	}

	//If more than one molecule to be output, write <cml> at start and </cml> at end.
	if((pConv->GetOutputIndex()==1))
	{
		if(!pConv->IsOption("x"))
			ofs << "<?xml version=\"1.0\"?>\n";
		if(!pConv->IsLast())
			ofs << "<" << ns << "cml>" << endl;
	}

	
	ofs << "<" << ns << "molecule";
	const char* title = mol.GetTitle();
	if(*title)
		ofs << " title=\"" << title << "\"";
	ofs << ">";

	if(mol.NumAtoms()>0)
	{
		ofs << "\n  <" << ns << "atomArray";
		if(!arrayform)
			ofs << ">\n";

		#ifdef HAVE_SSTREAM
				stringstream id, eltyp, chg, spn, x, y, z;
		#else
				strstream id, eltyp, chg, spn, x, y, z;
		#endif

		bool anyChg=false, anySpin=false;

		OBAtom *patom;
		vector<OBNodeBase*>::iterator i;
		for (patom = mol.BeginAtom(i);patom;patom = mol.NextAtom(i))
		{
			string el(etab.GetSymbol(patom->GetAtomicNum()));
			if(el=="Xx")
				el="R";

			int charge = patom->GetFormalCharge();
			if(charge)
				anyChg = true;

			int spin = patom->GetSpinMultiplicity();
			if(spin)
				anySpin = true;
			int isotope =patom->GetIsotope();
			
			if(arrayform)
			{
				id << " " << "a" << patom->GetIdx();
				eltyp << " " << el;
				chg << " " << charge;
				spn << " " << spin;
				x << " " << patom->GetX();
				y << " " << patom->GetY();
				x << " " << patom->GetZ();
			}
			else
			{
				ofs << "    <" << ns << "atom id=\"a" << patom->GetIdx() << "\"";
		//			if(cml1)
		//				ofs << ">\n<string builtin=\"elementType\">" << el << "</string>\n";
		//			else
				ofs << " elementType=\"" << el << "\"";

				if(isotope)
					ofs << " isotope=\"" << isotope << "\"";

				if(charge)
					ofs << " " << " formalCharge=\"" << charge << "\""; 

				if(spin)
					ofs << " " << " spinMultiplicity\"" << spin << "\""; 
				
				if(dim==2 || dim==3)
				{
					ofs << " x" << dim << "=\"" << patom->GetX() << "\"";
					ofs << " y" << dim << "=\"" << patom->GetY() << "\"";
				}
				if(dim==3)
					ofs << " z" << dim << "=\"" << patom->GetZ() << "\"";

				ofs << "/>\n";
			}
		}
		if(arrayform)
		{
			ofs << "\n    atomID=\"" << id.str() << "\"";
			ofs << " elementType=\"" << eltyp.str() << "\"";
			if(anyChg)
				ofs << "\n    formalCharge=\"" << chg.str() << "\"";
			if(anySpin)
				ofs << "\n    spinMultiplicity=\"" << spn.str() << "\"";
			if(dim==2 || dim==3)
			{
				ofs << "\n    x" << dim << "=\"" << x.str() << "\"";
				ofs << "\n    y" << dim << "=\"" << y.str() << "\"";
			}
			if(dim==3)
				ofs << "\n    z" << dim << "=\"" << z.str() << "\"";
			ofs << "/>\n";
			ofs << "  <" << ns << "bondArray";
		}
		else
			ofs << "  </"<< ns << "atomArray>\n  <"<< ns << "bondArray>\n";

		stringstream ref1, ref2, ord;
		OBBond *pbond;
		vector<OBEdgeBase*>::iterator ib;
		for (pbond = mol.BeginBond(ib);pbond;pbond = mol.NextBond(ib))
		{
			stringstream sbo;
			int bo = pbond->GetBO();
			if(bo==5) //aromatic
				sbo << "A";
			else
				sbo << bo;

			if(arrayform)
			{
				ref1 << " a" << pbond->GetBeginAtomIdx();
				ref2 << " a" << pbond->GetEndAtomIdx();
				ord  << " " << sbo.str();
			}
			else
			{
				ofs << "    <" << ns << "bond id=\"b" << pbond->GetIdx() << "\"";
				ofs << " atomRefs2=\"a" << pbond->GetBeginAtomIdx() 
						<< " a" << pbond->GetEndAtomIdx() <<"\"";
				ofs << " " << " order=\"" << sbo.str() << "\"";
				ofs << "/>\n";
			}
		}
		if(arrayform)
		{
			ofs << "\n    atomRef1=\"" << ref1.str() << "\"";
			ofs << "\n    atomRef2=\"" << ref2.str() << "\"";
			ofs << "\n    order=\"" << ord.str() << "\"";
			ofs << "/>\n";
		}
		else
			ofs << "  </"<< ns << "bondArray>\n";
	}
	ofs << "</"<< ns << "molecule>\n";

	if((pConv->GetOutputIndex()>1) && pConv->IsLast())
		ofs << "</" << ns << "cml>" << endl;

	return true;
}


}//namespace

/**********************************************************************
Copyright (C) 2005 Chris Morley
 
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
#include "inchi_api.h"
#include <sstream>
#include <set>
#include <vector>

using namespace std;
namespace OpenBabel
{

//This will eventually taken from NIST code
struct ElData
{
	char sym[3];
	int MW;
	int metal;
	bool AddH;
	int Valence[5][4]; // 4 possible values each for charge -2,-2,0,+1,+2
};

class InChIFormat : public OBMoleculeFormat
{
public:
	InChIFormat()
	{
			OBConversion::RegisterFormat("inchi",this);
	}

	virtual const char* Description() //required
	{
			return 
"INChI format\n \
IUPAC/NIST molecular identifier\n \
Options e.g. -xat \n \
 t add molecule name\n \
 a output auxilliary information\n \
 u output only unique molecules\n \
 U output only unique molecules and sort them\n \
 e compare first molecule to others\n \
 w don't warn on undefined stereo-chemistry alone\n \
";
	};

	virtual const char* SpecificationURL()
	{ return "";}; //optional

  virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
	virtual bool ReadChemObject(OBConversion* pConv);
	virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);
					char CompareInchi(const char* Inchi1, const char* Inchi2);
  virtual unsigned int Flags(){return NOTREADABLE;}; //InChI input disabled pending improvement

private:
	bool ReadInChI(istream& ifs, OBConversion* pConv);
	
	void GetComponent(string& ln, vector<string>& section, vector<char>& sectionchar);

	inline bool DoFormula(string s, OBMol& mol);
	inline bool DoConnections(string s, OBMol& mol);
	inline bool DoHydrogens(string s, OBMol& mol);
	inline bool Dodb(string s, OBMol& mol);
	inline bool Dosp3(string s, OBMol& mol);
	inline bool Dosp3inv(string s, OBMol& mol);
	inline bool DoStereotype(string s, OBMol& mol);

	bool AddHydrogens(OBMol& mol, int AtomIndex, int nH);
	bool AssignBondOrders(OBMol& mol, int Charge);
	bool CalcDeficits(OBMol& mol, vector<int>& Deficits);
	bool FixProtons(OBMol mol, int protons);

	static vector<ElData> Els;
	static bool ElDataRead;	
	bool ReadElData();

	OBAtom* GetCommonAtom(OBBond* pb1, OBBond* pb2);

	///Compare std::strings with embedded numbers so that 
	// "a6b" (or "a06b") is less than "a15b"
	// and "CH4" is less than "C2H6"
	// and "CH4" is less than "ClH" (hydrogen chloride)
	struct InchiLess
		: public binary_function<const string&, const string&, bool>
	{
		bool operator()(const string& s1, const string& s2) const
		{
			string::const_iterator p1, p2;
			p1=s1.begin(); p2=s2.begin();
			while( p1!=s1.end() && p2!=s2.end() )
			{
			  if(iscntrl(*p1) || iscntrl(*p2) || isspace(*p1) || isspace(*p2))
			    return false; //stop comparison at whitespace. Identical up to here 
			  int n1=-1,n2=-1;
			  if(isdigit(*p1))
			    {
			      n1 = atoi(&*p1);
			      //skip over number
			      while(p1!=s1.end() && isdigit(*p1++)); --p1;
			    }
			  if(isdigit(*p2))
			    {
			      n2 = atoi(&*p2);
			      while(p2!=s2.end() && isdigit(*p2++)); --p2;
			    }
			  if(n1<0 && n2 < 0)
			    {
			      //neither numbers
			      if(*p1 != *p2)
				return *p1 < *p2;
			    }
			  else if(n1>=0 && n2>0)
			    {
			      //both numbers
			      if(n1!=n2)
				return n1 < n2;
			    }
			  else if(n1>0)
			    return islower(*p2)!=0;
			  else if(n2>0)
			    return !islower(*p1);

			  ++p1; ++p2; // iterate
			} // while loop
			return false; //identical
		}
	};

	enum charges {MINUS2, MINUS1, ZERO, PLUS1, PLUS2};
  typedef	set<string, InchiLess> nSet;
//  typedef	set<string> nSet; //lexicographical sorting
	nSet allInchi;
	string firstInchi;
	string firstID;
};

bool InChIFormat::ElDataRead=false;
vector<ElData> InChIFormat::Els;

//Make an instance of the format class
InChIFormat theInChIFormat;

/////////////////////////////////////////////////////////////////
bool InChIFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
	OBMol* pmol = dynamic_cast<OBMol*>(pOb);
	if(pmol==NULL) return false;
	OBMol &mol = *pmol;

	stringstream molID;
	if(strlen(mol.GetTitle())==0)
		molID << '#' << pConv->GetOutputIndex() << ' ';
	else
		molID << mol.GetTitle() << ' ';
	if(pConv->GetOutputIndex()==1)
		firstID=molID.str();
	
	mol.AddHydrogens(false,false); //so stereo works

	inchi_Input inp;
	memset(&inp,0,sizeof(inchi_Input));
	bool Is0D=true;
	OBAtom* patom;
  vector<inchi_Atom> inchiAtoms(mol.NumAtoms());
	vector<OBNodeBase*>::iterator itr;
  for (patom = mol.BeginAtom(itr);patom;patom = mol.NextAtom(itr))
	{
		//OB atom index starts at 1; inchi atom index starts at 0
		inchi_Atom& iat = inchiAtoms[patom->GetIdx()-1];
		memset(&iat,0,sizeof(inchi_Atom));
		iat.x = patom->GetX();
		iat.y = patom->GetY();
		iat.z = patom->GetZ();
		if(iat.x!=0 || iat.y!=0 || iat.z!=0)
			Is0D=false;
		
		int nbonds = 0;
		vector<OBEdgeBase*>::iterator itr;
	  OBBond *pbond;
		for (pbond = patom->BeginBond(itr);pbond;pbond = patom->NextBond(itr))
		{
			// do each bond only once. Seems necessary to avoid problems with stereo
			if(pbond->GetNbrAtomIdx(patom)<patom->GetIdx()) continue;
			iat.neighbor[nbonds]      = pbond->GetNbrAtomIdx(patom)-1;
			int bo = pbond->GetBO();
			if(bo==5)
				bo=4;
			iat.bond_type[nbonds]     = bo;

			S_CHAR bs = INChI_BOND_STEREO_NONE;

			if(pbond->IsWedge())
				bs = INChI_BOND_STEREO_SINGLE_1UP;
			if(pbond->IsHash())
				bs = INChI_BOND_STEREO_SINGLE_1DOWN;

			iat.bond_stereo[nbonds++] = bs;
			if(nbonds>MAXVAL)
			{
				cerr << "Too many bonds to " << iat.elname << " atom" << endl;
				return false;
			}
		}
	
		strcpy(iat.elname,etab.GetSymbol(patom->GetAtomicNum()));
		iat.num_bonds = nbonds;
		iat.num_iso_H[0] = -1; //Let inchi add implicit Hs
		if(patom->GetIsotope())
		{
			iat.isotopic_mass = ISOTOPIC_SHIFT_FLAG +
				patom->GetIsotope() - (int)(isotab.GetExactMass(patom->GetAtomicNum())+0.5);
		}
		else
			iat.isotopic_mass = 0 ;
		iat.radical = patom->GetSpinMultiplicity();
		iat.charge  = patom->GetFormalCharge();
	}
	
	inp.atom = &inchiAtoms[0];

	vector<inchi_Stereo0D> stereoVec;
	if(Is0D)
	{
		//Tetrahedral stereo
		mol.FindChiralCenters();
	  OBAtom* patom;
		vector<OBNodeBase*>::iterator itr;
		if(mol.IsChiral())
		{
			for(patom = mol.BeginAtom(itr);patom;patom = mol.NextAtom(itr))
			{
				if(patom->IsChiral())
				{
					inchi_Stereo0D stereo;
					stereo.central_atom = patom->GetIdx()-1;
					stereo.type = INChI_StereoType_Tetrahedral;
					int i=0;
					vector<OBEdgeBase*>::iterator itr;
					OBBond *pbond;
					for (pbond = patom->BeginBond(itr);pbond;pbond = patom->NextBond(itr),++i)
						stereo.neighbor[i] = pbond->GetNbrAtomIdx(patom)-1;
					//if(i!=4)...error

					stereo.parity = INChI_PARITY_UNKNOWN;
					if(patom->IsPositiveStereo() || patom->IsClockwise())
						stereo.parity = INChI_PARITY_ODD;
					if(patom->IsNegativeStereo() || patom->IsAntiClockwise())
						stereo.parity = INChI_PARITY_EVEN;
					stereoVec.push_back(stereo);
				}
			}
		}
		
		//Double bond stereo
		//Currently does not handle cumulenes
		set<OBBond*> UpDown;
		set<OBBond*>::iterator uditr;
		OBBond *pbond;
		vector<OBEdgeBase*>::iterator bitr;
		for (pbond = mol.BeginBond(bitr);pbond;pbond = mol.NextBond(bitr))
		{
		 if(pbond->IsUp() || pbond->IsDown())
			 UpDown.insert(pbond);
		}
		if(!UpDown.empty())
		{
			//For all double bonds look for up/down bonds attached
			for (pbond = mol.BeginBond(bitr);pbond;pbond = mol.NextBond(bitr))
			{
				if(pbond->IsDouble())
				{
					OBAtom* patom, *pA=NULL, *pB=NULL;
					OBBond* pX, *pY;;
					for(uditr=UpDown.begin();uditr!=UpDown.end();++uditr)
					{
						patom = GetCommonAtom(pbond, *uditr);
						if(patom && !pA)
						{
							//first one in pA
							pA = patom;
							pX = *uditr;
						}
						else if(patom && pA)
						{
							//second one in pB
							pB = patom;
							pY = *uditr;
							break;
						}
					}
					if(pA && pB)
					{
						inchi_Stereo0D stereo;
						stereo.central_atom = NO_ATOM;
						stereo.type = INChI_StereoType_DoubleBond;
						stereo.neighbor[0]= pX->GetNbrAtomIdx(pA)-1;
						stereo.neighbor[1] = pA->GetIdx()-1;
						stereo.neighbor[2] = pB->GetIdx()-1;
						stereo.neighbor[3]= pY->GetNbrAtomIdx(pB)-1;
						if((pX->IsUp() && pY->IsUp())||(pX->IsDown() && pY->IsDown()))
							stereo.parity = INChI_PARITY_ODD;
						else
							stereo.parity = INChI_PARITY_EVEN;
						stereoVec.push_back(stereo);
					}
				}
			}
		}
	}

//	*inp.szOptions = '\0';
	inp.num_atoms = mol.NumAtoms();
	inp.stereo0D = &stereoVec[0];
	inp.num_stereo0D = stereoVec.size();

/* CM print out input data
#ifdef _DEBUG
{
	int cmi;
	TRACE("Options'%s'\n",inp.szOptions);
	TRACE("num_atoms=%d\n",inp.num_atoms);
	for(cmi=0;cmi<inp.num_atoms;++cmi)
	{
		inchi_Atom* pat;
		int ib;
		pat=&inp.atom[cmi];
		TRACE("%s num_bonds=%d iso=%d rad=%d charge=%d\n",
			pat->elname, pat->num_bonds, pat->isotopic_mass,pat->radical,pat->charge);
		for(ib=0;ib<pat->num_bonds;++ib)
		{
			int i;
			TRACE("IsoH: %d %d %d %d\n",
				pat->num_iso_H[0],pat->num_iso_H[1],pat->num_iso_H[2],pat->num_iso_H[3]);
			for(i=0;i<MAXVAL;++i)
			{
				if(pat->bond_type[i]==0) break;
				TRACE("%d %d %d\n",pat->neighbor[i],pat->bond_type[i],pat->bond_stereo);
			}
		}
	}
	TRACE("num_stereo0D=%d\n",inp.num_stereo0D);
	for(cmi=0;cmi<inp.num_stereo0D;++cmi)
	{
		inchi_Stereo0D* pst;
		pst=&(inp.stereo0D[cmi]);
		TRACE("neighbors: %d %d %d %d\n",
			pst->neighbor[0],pst->neighbor[1],pst->neighbor[2],pst->neighbor[3]);
		TRACE("type=%d centatom=%d parity=%d\n",pst->type,pst->central_atom,pst->parity);
	}
}
#endif
*/

	inchi_Output inout;
	memset(&inout,0,sizeof(inchi_Output));
	int ret = GetINChI(&inp, &inout);
	if(ret!=inchi_Ret_OKAY)
	{
		string mes(inout.szMessage);
		if(!(pConv->IsOption('w') && mes=="Omitted undefined stereo"))
			cerr << molID.str() << inout.szMessage << endl;
		if(ret!=inchi_Ret_WARNING)
		{
			FreeINChI(&inout);
			return false;
		}
	}
	
	string ostring = inout.szINChI;
	if(pConv->IsOption('t'))
	{
		ostring += ' ';
		ostring +=  mol.GetTitle();
	}
		
	ostream &ofs = *pConv->GetOutStream();

	if(pConv->IsOption('U'))
	{
		if(pConv->GetOutputIndex()==1)
			allInchi.clear();
		//Just add to set and don't output, except at the end
		allInchi.insert(ostring);
		
		if(pConv->IsLast())
		{
			nSet::iterator itr;
			for(itr=allInchi.begin();itr!=allInchi.end();++itr)
				ofs << *itr << endl;
		}
		return true;
	}
	else if(pConv->IsOption('u'))
	{
		if(pConv->GetOutputIndex()==1)
			allInchi.clear();
		if(!allInchi.insert(ostring).second)
			return true; //no output if already in set
	}

	ofs << ostring << endl;

	if (pConv->IsOption('a'))
		ofs << inout.szAuxInfo << endl;
	
	if(pConv->IsOption('e'))
	{
		if(pConv->GetOutputIndex()==1)
			firstInchi = inout.szINChI;
		else
		{
			ofs << "Molecules " << firstID << "and " << molID.str();
			switch (CompareInchi(firstInchi.c_str(), inout.szINChI)) 
			{
			case 0:
				ofs << " are identical" << endl;
				break;
			case '+':
				ofs << " have different formulae" << endl;
				break;
			case 'c':
				ofs << " have different connection tables" << endl;
				break;
			case 'h':
				ofs << " have different bond orders, or radical character" << endl;
				break;
			case 'q':
				ofs << " have different charges" << endl;
				break;
			case 'p':
				ofs << " have different numbers of attached protons" << endl;
				break;
			case 'b':
				ofs << " have different double bond stereochemistry" << endl;
				break;
			case 'm':
			case 't':
				ofs << " have different sp3 stereochemistry" << endl;
				break;
			case 'i':
				ofs << " have different isotopic composition" << endl;
				break;
			default:
				ofs << " are different" << endl;
			}
		}
	}

	FreeINChI(&inout);
	return true;
}

////////////////////////////////////////////////////////////////
bool InChIFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
{
//**This routine needs thinking about
  OBMol* pmol = dynamic_cast<OBMol*>(pOb);
  if(pmol==NULL)
      return false;

  //Define some references so we can use the old parameter names
  istream &ifs = *pConv->GetInStream();
  OBMol &mol = *pmol;

	mol.BeginModify();
	
	bool ret = ReadInChI(ifs, pConv);	
	
	mol.EndModify();
	return ret;
}

bool InChIFormat::ReadChemObject(OBConversion* pConv)
{
	istream &ifs = *pConv->GetInStream();
	bool ret = ReadInChI(ifs, pConv);	
	return ret;
}

/*
From The IUPAC Chemical Identifier - Technical Manual
An algorithm to parse the Identifier may be described in the following way:
1)	Find the first slash. The slash is preceded by the version and followed by a string 
(call it S) that contains all other layers of the identifier.
2)	Search for "/r" in S. If "/r" is found then copy preceding "/r" substring to P[1] 
and the following "/r" string to P[2] else copy S to P[1] 
(P[1] represents the whole identifier or an identifier of a disconnected structure; 
P[2] if not empty represents an identifier of a "reconnected" structure".)
3)	Search for "/f" in each non-empty P[i]. If "/f" was found then copy the 
preceding string to Q[i][1] and the following string to Q[i][2] else copy P[i] to Q[i][1]
 (Q[i][1] represents the Main layer; Q[i][2] represents fixed-H layer)
4)	Search for "/i" in each non-empty Q[i][j]. If "/i" was found then copy the preceding 
string into R[i][j][1] and the following string into R[i][j][2] else copy Q[i][j] to R[i][j][1]
(R[i][j][1] represents the non-isotopic part of the layer;  R[i][j][2] represents the isotopic layer)

At the end, non-empty strings R[i][j][k] (i, j, k = 1 or 2) contain:
i = 1: The identifier or the identifier of a disconnected structure
i = 2: The identifier of the "reconnected" structure
j = 1: The main layer
j = 2: The fixed-H layer
k = 1: The non-isotopic part
k = 2: The isotopic part of the layer

In case of multicomponent compound the parts of the identifier related to 
components are separated by semicolons ";" except for the chemical formula
which is dot-disconnected. The order of the components within the segments 
of the Main layer is same; the fixed-H layer may have a different order of 
the components. In this case a transposition segment (/o) is present. For 
example, transposition (1,2,3) means that component #1 in the Main layer is 
component #2 in Fixed-H layer, component #2 in the Main layer is component 
#3 in Fixed-H layer and component #3 in the Main layer is component #1 in 
the Fixed-H layer. A simple example of a compound that exhibits a transposition 
is on Fig. A3-3 later.

To extract identifiers for individual components 
parse the identifier and obtain array of strings R as explained above
split each of those strings into segments using "/?" as separators and 
identify the separators (see Fig. 2)
split each segment by locating dots or semicolons into parts related to 
individual components and expanded them in case of multipliers and/or 
abbreviations describer in Appendix 2;
transpose components in fixed-H part according to transposition (/o) if 
it is present pick the first entries (corresponding to the first component) and merge 
them together using previously found "/?" separators; add "version/" to 
the beginning of the string. The string is an identifier for the first 
component repeat for all other components

*/
bool InChIFormat::ReadInChI(istream& ifs, OBConversion* pConv)
{
	string ln, version, reconnected, fixedH, title;
	getline(ifs,ln);
	string::size_type pos = ln.find(' ');
	if(pos!=string::npos)
	{
		title = ln.substr(pos+1); 
		ln.erase(pos);
	}

	pos = ln.find('/');
	if(pos!=string::npos)
		version =ln.substr(0,pos);
	else
		return false;
		
	pos = ln.find("/r");
	if(pos!=string::npos)
	{
		reconnected=ln.substr(pos+2);
		ln.erase(pos);	
	}	

	do //First for base structure, then for reconnected 
	{
		pos=ln.find("/f");
		if(pos!=string::npos)
		{
			fixedH=ln.substr(pos+2);
			ln.erase(pos);	
		}	

		//Split the InChI into sections, extract the first component to
		//the section vector and reassemble the rest back into ln
		while(!ln.empty())
		{
			vector<string> section;
			vector<char> sectionchar;
			GetComponent(ln, section, sectionchar);
		
			OBMol* pmol = new OBMol;
			OBMol& mol = *pmol;
			mol.BeginModify();
			mol.SetTitle(title);

			if(!DoFormula(section[1],mol)) return false;

			int n=2;
			//Connections layer
			if(sectionchar[n]=='c')
			{
				if(!DoConnections(section[n++],mol)) 
					return false;
			}
			
			//Hydrogens Layer
			if(sectionchar[n]=='h')
			{
				if(!DoHydrogens(section[n++],mol)) 
					return false;
			}

			//CHARGE Layer
			int charge=0;
			if(sectionchar[n]=='q')
			{
				stringstream ss(section[n++]);
				ss >> charge;
			}
			
			int protons=0;
			if(sectionchar[n]=='p')
			{
				stringstream ss(section[n++]);
				ss >> protons;
			}

			if(!AssignBondOrders(mol,charge)) return false;
			if(protons)
				if(!FixProtons(mol,protons)) return false;
			
			//STEREO Layer
			if(sectionchar[n]=='b')
			{
				if(!Dodb(section[n++],mol)) 
					return false;
			}
			if(sectionchar[n]=='t')
			{
				if(!Dosp3(section[n++],mol)) 
					return false;
			}
			if(sectionchar[n]=='m')
			{
				if(!Dosp3inv(section[n++],mol)) 
					return false;
			}
			if(sectionchar[n]=='s')
			{
				if(!DoStereotype(section[n++],mol)) 
					return false;
			}
			
			//ISOTOPE Layer	
	/*		if(ch=='c')
			{
				if(!DoConnections(section[n++],mol)) 
					return false;
			}
			if(ch=='c')
			{
				if(!DoConnections(section[n++],mol)) 
					return false;
			}

			if(ln==reconnected)
				break;
			else
				ln=reconnected;
	*/
			//Provide the current component
			mol.EndModify();
			if(pConv->AddChemObject(mol.DoTransformations(pConv->GetGeneralOptions())) <0)
				return false;
		}
	}while(!ln.empty());
	return true;
}

///////////////////////////////////////////////////////////
void InChIFormat::GetComponent(string& ln, vector<string>& section, vector<char>& sectionchar)
{
	//Split the InChI into sections, extract the first component to
	//the section vector and reassemble the rest back into ln
	//ln is returned empty when there are no more components
	tokenize(section,ln,"/");
	sectionchar.resize(section.size());
	ln = section[0]; // InChI= xxx/

	//Split each section into subsections for each molecular component
	//separated by . or ; and possibly with multiplier
	int i;
	for(i=1;i<section.size();++i)
	{
		int multiplier;
		bool hasmultiplier=false;
		stringstream ss(section[i]);
		if(i!=1) //no char on first section
			ss >> sectionchar[i];
		if(isdigit(ss.peek()))
		{
			//If there is a possible multiplier
			//Copy whole section with multiplier decreased to ln
			//and multiplied part to section, without subsequent components
			ss >> section[i];
			ss >> multiplier;
			if(ss.peek()=='*')
			{
				ss.ignore();
				hasmultiplier=true;
				ss >> section[i];

				stringstream ssln;
				ssln << '/' ;
				if(i!=1)
					ssln << sectionchar[i];
				if(--multiplier)
					ssln << multiplier << '*';
				ln += ssln.str();
				ln += section[i];
			}
		}

		string::size_type pos = section[i].find_first_of(".;");
		if(pos!=string::npos)
		{
			if(!hasmultiplier)
			{
				ln += '/';
				if(i!=1)
					ln += sectionchar[i];
				ln += section[i].substr(pos+1);
			}
			section[i].erase(pos);
		}
	}
	if(ln==section[0]) ln.erase(); // InChI=xxx/

}
///////////////////////////////////////////////////////////
bool InChIFormat::DoFormula(string s, OBMol& mol)
{
	stringstream ss(s);
	OBAtom atom;
	do
	{
		char sym[3] ={0,0,0};
		ss >> sym[0];
		char ch = ss.peek();
		if(isalpha(ch) && islower(ch))
			ss >> sym[1];
		int n=1;
		if(isdigit(ss.peek()))
			ss >> n;
		int atno = etab.GetAtomicNum(sym);
		if(atno!=1) //ignore hydrogens
		{
			do
			{
				atom.SetAtomicNum(atno);
				atom.SetType(sym);
				if(!mol.AddAtom(atom))
					return false;
				atom.Clear();
			}while(--n > 0);
		}
	}while(ss.good());
	return true;
}

bool InChIFormat::DoConnections(string s, OBMol& mol)
{
	int cur,nxt;
	vector<int> pend;
	char ch;
	stringstream ss(s);
	ss >> cur;
	do
	{
		ss >> ch >> nxt;
		if(ch==',') 
			cur = pend.back();
		else if(ch=='(')
			pend.push_back(cur);
		else if(ch==')')
		{
			cur = pend.back();
			pend.pop_back();
		}
		mol.AddBond(cur, nxt, 1);
		cur = nxt;
	}while(ss.good());
	return true;
}

bool InChIFormat::DoHydrogens(string s, OBMol& mol)
{
	stringstream ss(s);
	int startAt,endAt,nH;
	char ch;
	vector<int> ats;
	do
	{
		if(ss.peek()=='(')
		{
			//shared hydrogen
			//Assign to the first atoms in list for the moment
			while(ss.get()=='(')//all (H...) blocks. Assume all such blocks are together with no comma between
			{			
				ss >> ch;
				if(ch!='H')return false;
				if(isdigit(ss.peek()))
					ss >> nH;
				else
					nH=1;
				ss >> ch;
				if(ch!=',')return false;
				int i;
				for(i=0;i<nH;i++)
				{
					ss >> startAt >> ch;
					AddHydrogens(mol,startAt,1);
				}
				//skip over other atom numbers
				while(ch!=')')
				{
					ss >> ch;
				}
			}
		}
		else
		{
			do
			{
				ss >> startAt >> ch;
				ats.push_back(startAt);
				if(ch=='-')
				{
					ss >> endAt >> ch;
					int i;
					for(i=startAt+1;i<=endAt;++i)
						ats.push_back(i);
				}
			}while(ch!='H');
			if(isdigit(ss.peek()))
				ss >> nH;
			else
				nH=1;
			vector<int>::iterator itr;
			for(itr=ats.begin();itr!=ats.end();++itr)
			{
				AddHydrogens(mol,*itr,nH);
			}
			ats.clear();
			ss >> ch;
		}
	}while(ss.good());
	return true;
}

bool InChIFormat::Dodb(string s, OBMol& mol)
{
	stringstream ss(s);
	return true;
}

bool InChIFormat::Dosp3(string s, OBMol& mol)
{
	stringstream ss(s);
	return true;
}

bool InChIFormat::Dosp3inv(string s, OBMol& mol)
{
	stringstream ss(s);
	return true;
}

bool InChIFormat::DoStereotype(string s, OBMol& mol)
{
	stringstream ss(s);
	return true;
}

//////////////////////////////////////////////////////
bool InChIFormat::AddHydrogens(OBMol& mol, int AtomIndex, int nH)
{
	//Adds nH hydrogens to the atom
	int i;
	for(i=0;i<nH;i++)
	{
		OBAtom* h = mol.NewAtom();
		if(!h) return false;
		h->SetAtomicNum(1);
		h->SetType("H");
    if(!mol.AddBond(AtomIndex,h->GetIdx(),1)) return false;
	}
	return true;
}
/*
Given: connections of all atoms; overall charge or number of extra protons
assign bond orders.
use standard valence
If extra protons, find any atoms with more Hs than appropriate. Assign charge to it.
(May not work e.g. C=O-H not seen until bond)
Calculate valency deficit for all atoms
If deficit is negative and there is a positive charge, assign it to that atom; makes deficit less
For each atom
	if deficit
		choose a bond to another deficit atom
			if there isn't one assign radical centre
		increase bond order
		scan atoms again
if no atoms deficit
	if >1 radical centre: try again with different choice of starting bond
	if 0 or 1 radical centre: success

*/
bool InChIFormat::AssignBondOrders(OBMol& mol, int Charge)
{
	//Find atoms which need more bonds
	vector<int> Deficits(mol.NumAtoms()+1); //atoms index starts at 1
	
	int tries=2; //may need to be more sophisticated
	int skip=0; //to adjust choice of bond on later tries 
	//Make a copy of the molecule to try out bonding solutions on
	OBMol tmol = mol;
	do
	{
		CalcDeficits(mol,Deficits);
		int UnfixedDeficits=0;
		int i;
		for(i=1;i<=mol.NumAtoms();i++)
		{
			if(Deficits[i]>0)
			{
				//For each atom with a deficit
				OBAtom* patom = mol.GetAtom(i);
				OBAtom* pnbr;
				vector<OBEdgeBase*>::iterator itr;
				for (pnbr = patom->BeginNbrAtom(itr);pnbr;pnbr = patom->NextNbrAtom(itr))
				{
					//Find an adjacent atom with a deficit
					if(Deficits[pnbr->GetIdx()]>0)
					{
						if(skip)
						{
							skip--;
							continue;
						}
						//increment bond order
						OBBond* pbond = (OBBond*)(*itr);
						pbond->SetBO(pbond->GetBO() + 1);
						--Deficits[pnbr->GetIdx()];
						--Deficits[patom->GetIdx()];
						break;
					}
				}
				if(pnbr==NULL)
				{
					//Could not find an adjacent deficit atom

					//If there is a charge
					
					//See if an adjacent atom can be given a higher valence
					for (pnbr = patom->BeginNbrAtom(itr);pnbr;pnbr = patom->NextNbrAtom(itr))
					{
						ElData el = Els[patom->GetAtomicNum()];
						int i;
						for (i=0;i<4;i++)
						{
							int newVal = el.Valence[ZERO][i];
							if(newVal==0) break;
							if(newVal > patom->BOSum())
							{
								//Increment bond
								OBBond* pbond = patom->GetBond(pnbr);
								pbond->SetBO(pbond->GetBO() + 1);

								int prevVal = i ? el.Valence[ZERO][i-1] : el.Valence[ZERO][i-1]; 
								Deficits[pnbr->GetIdx()] += newVal - prevVal - 1; 
								--Deficits[patom->GetIdx()];
								break;
							}
						}
					}					
					if(pnbr==NULL)
					{
						patom->SetSpinMultiplicity(Deficits[i]+1); //for rads and (triplet)carbenes
						++UnfixedDeficits;
					}
				}
			}
		}
		if(UnfixedDeficits==0)	return true;
		
		//Need to try again
		skip++;
		mol=tmol;

	}while(tries--);
	//return false;
	string title = mol.GetTitle();
	title += " Bond Assignment failed";
	mol.SetTitle(title);
	return true;
}

////////////////////////////////////////////////
bool InChIFormat::CalcDeficits(OBMol& mol, vector<int>& Deficits)
{
	// Calculates valency deficits for each atom. 

	enum charges {MINUS2, MINUS1, ZERO, PLUS1, PLUS2};
	if(!ElDataRead) ReadElData();

	OBAtom *atom;
	vector<OBNodeBase*>::iterator i;
	for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
	{
		if(atom->IsHydrogen()) continue;
		ElData el = Els[atom->GetAtomicNum()];
		int deficit = el.Valence[ZERO][0] - atom->BOSum();
		Deficits[atom->GetIdx()] = deficit;
	}
	return true;
}

////////////////////////////////////////////////
bool InChIFormat::ReadElData()
{
	const char* filename = "INCHIvalence2.csv";
	ifstream ifs(filename);
	if(!ifs)
	{
		cerr << "Couldn't open " << filename << endl;
		return false;
	}
	string header;
	getline(ifs,header);
	ElData el;
	Els.push_back(el); //dummy atNum=0 (contains rubbish)
	while(ifs.good())
	{
		memset(&el,0,sizeof(ElData));
		
		char ch;
		ifs.getline(el.sym,3,',');
		ifs >> el.MW >> ch >> el.metal >> ch >> el.AddH >> ch;
		
		int nchg,j;
		for(nchg=0;nchg<5;++nchg)
		{
			for(j=0;j<4 && ifs.peek()!=',' && ifs.peek()!='\n';++j)
			{
				if(!(el.Valence[nchg][j] = ifs.get() - '0'))
					break;
			}
			ifs.get(); // ,
		}
		Els.push_back(el);
	}
	return true;	
}

//////////////////////////////////////////////////////
bool InChIFormat::FixProtons(OBMol mol, int protons)
{
	int targ[] = {7,15,8}; //look first for N then P, O atoms to attach protons to
	int i;
	for(i=0;i<3;++i)
	{
		do
		{
			OBAtom* atom;
			vector<OBNodeBase*>::iterator itr;
			for (atom = mol.BeginAtom(itr);atom;atom = mol.NextAtom(itr))
			{
				if((atom->GetAtomicNum()==targ[i]) && (atom->GetValence()<=3))
				{
					AddHydrogens(mol,atom->GetIdx(),1);
					atom->SetFormalCharge(+1);
					--protons;
					break;
				}
				if(!atom) return false;
			}
		}while(protons);
	}
	return true;
}

//////////////////////////////////////////////////////////
char InChIFormat::CompareInchi(const char* Inchi1, const char* Inchi2)
{
	//Returns 0 if identical or an char identifying the layer where they first differed
	string s1(Inchi1), s2(Inchi2);

	//Remove anything after the end of the Inchi
	string::size_type pos;
	pos = s1.find_first_of(" \t\n");
	if(pos!=string::npos)
		s1.erase(pos);
	pos = s2.find_first_of(" \t\n");
	if(pos!=string::npos)
		s2.erase(pos);

	vector<string> layers1, layers2;
	tokenize(layers1,s1,"/\n");
	tokenize(layers2,s2,"/\n");
	int i;
	if(layers1.size()<layers2.size())
		layers1.swap(layers2); //layers1 is the longest
	
	for(i=1;i<layers2.size();++i)
	{
		if(layers1[i]!=layers2[i])
		{
			char ch = '+';
			if(i>1) //not formula layer
				ch=layers1[i][0];
			return ch;
		}
	}
	if(layers1.size()==layers2.size())	
		return 0;
	else
		return layers1[i][0];
}

OBAtom* InChIFormat::GetCommonAtom(OBBond* pb1, OBBond* pb2)
{
	OBAtom* pa1 = pb1->GetBeginAtom();
	if(pa1==pb2->GetBeginAtom() || pa1==pb2->GetEndAtom())
		return pa1;
	pa1 = pb1->GetEndAtom();
	if(pa1==pb2->GetBeginAtom() || pa1==pb2->GetEndAtom())
		return pa1;
	return NULL; //not adjacent bonds
}

//******************************************************
class InChICompareFormat : public OBMoleculeFormat
{
public:
	InChICompareFormat()
	{
			OBConversion::RegisterFormat("k",this);
	}

	virtual const char* Description() //required
	{
			return 
"Compares first molecule to others using INChI.\n \
Same as -oinchi -xet\n \
";
	};
	
	virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv)
	{
		pConv->SetOptions("et");
		return theInChIFormat.WriteMolecule(pOb,pConv);
	};
	
	virtual unsigned int Flags() { return NOTREADABLE;};
};

InChICompareFormat theInChICompareFormat;

}//namespace OpenBabel


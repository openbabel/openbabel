/**********************************************************************
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (c) 2001-2003 by Geoffrey R. Hutchison
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
#include "math/matrix3x3.h"
#include "obconversion.h"

using namespace std;
namespace OpenBabel {

class CacaoFormat : public OBFormat
{
public:
	//Register this format type ID
	CacaoFormat() {OBConversion::RegisterFormat("CACAO",this);}

	virtual const char* Description() //required
	{ return
"Cacao format\n \
No comments yet\n \
";
	};

	virtual const char* SpecificationURL(){return
		"http://www.chembio.uoguelph.ca/oakley/310/cacao/cacao.htm";}; //optional

	//Flags() can return be any the following combined by | or be omitted if none apply
	// NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY 
	virtual unsigned int Flags(){return READONEONLY;};

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

	static void SetHilderbrandt(OBMol&,vector<OBInternalCoord*>&);

};

//Make an instance of the format class
CacaoFormat theCacaoFormat;

/////////////////////////////////////////////////////////////////
bool CacaoFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
{

	OBMol* pmol = dynamic_cast<OBMol*>(pOb);
	if(pmol==NULL) return false;

	//Define some references so we can use the old parameter names
	istream &ifs = *pConv->GetInStream();
	OBMol &mol = *pmol;
	const char* title = pConv->GetTitle();

  char buffer[BUFF_SIZE];
  int natoms;
  double A,B,C,Alpha,Beta,Gamma;
  matrix3x3 m;
  
  ifs.getline(buffer,BUFF_SIZE); mol.SetTitle(buffer);
  ifs.getline(buffer,BUFF_SIZE); sscanf(buffer,"%d",&natoms);
  
  while (ifs.getline(buffer,BUFF_SIZE) &&!EQn(buffer,"CELL",4));

  if (!EQn(buffer,"CELL",4)) return(false);
  vector<string> vs;
  tokenize(vs,buffer," \n\t,");
  if (vs.size() != 7) return(false);

  //parse cell values
  A = atof((char*)vs[1].c_str());
  B = atof((char*)vs[2].c_str());
  C = atof((char*)vs[3].c_str());
  Alpha = atof((char*)vs[4].c_str());
  Beta  = atof((char*)vs[5].c_str());
  Gamma = atof((char*)vs[6].c_str());

  OBUnitCell *uc = new OBUnitCell;
  uc->SetData(A, B, C, Alpha, Beta, Gamma);
  mol.SetData(uc);
  m = uc->GetOrthoMatrix();

  int i;
  double x,y,z;
  char type[10];
  OBAtom *atom;
  vector3 v;

  for (i = 1; i <= natoms;i++)
    {
      if (!ifs.getline(buffer,BUFF_SIZE)) return(false);
      tokenize(vs,buffer," \n\t,");
      if (vs.size() < 4) return(false);
      atom = mol.NewAtom();

      x = atof((char*)vs[1].c_str());
      y = atof((char*)vs[2].c_str());
      z = atof((char*)vs[3].c_str());
      v.Set(x,y,z); v *= m;

      strcpy(type,vs[0].c_str());
      atom->SetAtomicNum(etab.GetAtomicNum(type));
      atom->SetVector(v);
    }

  mol.ConnectTheDots();
  mol.PerceiveBondOrders();

  return(true);
}

////////////////////////////////////////////////////////////////

bool CacaoFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
	OBMol* pmol = dynamic_cast<OBMol*>(pOb);
	if(pmol==NULL) return false;

	//Define some references so we can use the old parameter names
	ostream &ofs = *pConv->GetOutStream();
	OBMol &mol = *pmol;
	const char *dimension = pConv->GetDimension();

  OBAtom *atom;
  char buffer[BUFF_SIZE];
  vector<OBNodeBase*>::iterator i;

  sprintf(buffer,"%s",mol.GetTitle());
  ofs << buffer << endl;
  sprintf(buffer,"%3d   DIST  0  0  0",mol.NumAtoms());
  ofs << buffer << endl;

  if (!mol.HasData(obUnitCell))
    sprintf(buffer,"CELL 1.,1.,1.,90.,90.,90.");
  else
    {
      OBUnitCell *uc = (OBUnitCell*)mol.GetData(obUnitCell);
      sprintf(buffer,"CELL %f,%f,%f,%f,%f,%f",
	      uc->GetA(), uc->GetB(), uc->GetC(),
	      uc->GetAlpha(), uc->GetBeta(), uc->GetGamma());
    }
  ofs << buffer << endl;

  for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
  {
    sprintf(buffer,"%2s %7.4f, %7.4f, %7.4f",
            etab.GetSymbol(atom->GetAtomicNum()),
            atom->x(),
            atom->y(),
            atom->z());
    ofs << buffer << endl;
  }

  return(true);
}

void CacaoFormat::SetHilderbrandt(OBMol &mol,vector<OBInternalCoord*> &vit)
{
  // Roundtrip testing shows that some atoms are NULL
  //  which causes segfaults when dereferencing later
  //   (e.g. in the last "segment" of this routine
  double sum,r;

  OBAtom dummy1,dummy2;
  dummy1.SetVector(0.0,0.0,1.0);
  dummy2.SetVector(1.0,0.0,0.0);

  OBAtom *atom,*a1,*a2,*ref;
  vector<OBNodeBase*>::iterator ai;
  vit.push_back((OBInternalCoord*)NULL);
  for (atom = mol.BeginAtom(ai);atom;atom = mol.NextAtom(ai))
    vit.push_back(new OBInternalCoord (atom));

  vit[1]->_a = &dummy1;
  vit[1]->_b = &dummy2;
  vit[2]->_b = &dummy1;
  vit[2]->_c = &dummy2;
  vit[3]->_c = &dummy1;

  unsigned int i,j;
  for (i = 2;i <= mol.NumAtoms();i++)
    {
      ref = (OBAtom*)NULL;
      a1 = mol.GetAtom(i);
      sum = 100.0;
      for (j = 1;j < i;j++)
	{
	  a2 = mol.GetAtom(j);
	  r = (a1->GetVector()-a2->GetVector()).length_2();
	  if ((r < sum) && (vit[j]->_a != a2) && (vit[j]->_b != a2))
	    {
	      sum = r;
	      ref = a2;
	    }
	}
      vit[i]->_a = ref;
    }

  for (i = 3;i <= mol.NumAtoms();i++)
    vit[i]->_b = vit[vit[i]->_a->GetIdx()]->_a;

  for (i = 4;i <= mol.NumAtoms();i++)
    {
      if (vit[i]->_b && vit[i]->_b->GetIdx())
	vit[i]->_c = vit[vit[i]->_b->GetIdx()]->_b;
      else
	vit[i]->_c = &dummy1;
    }

  OBAtom *a,*b,*c;
  vector3 v1,v2;
  for (i = 2;i <= mol.NumAtoms();i++)
    {
      atom = mol.GetAtom(i);
      a = vit[i]->_a; b = vit[i]->_b; c = vit[i]->_c;
      v1 = atom->GetVector() - a->GetVector();
      v2 = b->GetVector() - a->GetVector();
      vit[i]->_ang = vectorAngle(v1,v2);
      vit[i]->_tor = CalcTorsionAngle(atom->GetVector(),
				      a->GetVector(),
				      b->GetVector(),
				      c->GetVector());
      vit[i]->_dst = (vit[i]->_a->GetVector() - atom->GetVector()).length();
    }

}
//***************************************************************
class CacaoInternalFormat : public OBFormat
{
public:
	//Register this format type ID
	CacaoInternalFormat() {OBConversion::RegisterFormat("CACAOINT",this);}

	virtual const char* Description() //required
	{ return
"CacaoInternal format\n \
No comments yet\n \
";
	};

	virtual const char* SpecificationURL(){return
		"http://www.chembio.uoguelph.ca/oakley/310/cacao/cacao.htm";}; //optional

	//Flags() can return be any the following combined by | or be omitted if none apply
	// NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY 
	virtual unsigned int Flags(){return NOTREADABLE;};

	////////////////////////////////////////////////////
	/// The "API" interface functions
	virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

////////////////////////////////////////////////////
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

//Make an instance of the format class
CacaoInternalFormat theCacaoInternalFormat;

bool CacaoInternalFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
	OBMol* pmol = dynamic_cast<OBMol*>(pOb);
	if(pmol==NULL) return false;

	//Define some references so we can use the old parameter names
	ostream &ofs = *pConv->GetOutStream();
	OBMol &mol = *pmol;
	const char *dimension = pConv->GetDimension();

  unsigned int i;
  vector3 v;
  char tmptype[10],buffer[BUFF_SIZE];

  if (mol.Empty()) return(false);

  //translate first atom to origin
  v = mol.GetAtom(1)->GetVector(); v *= -1.0; mol.Translate(v);

  vector<OBInternalCoord*> vit;
  CacaoFormat::SetHilderbrandt(mol,vit);
  strcpy(tmptype,etab.GetSymbol(mol.GetAtom(1)->GetAtomicNum()));

  sprintf(buffer," # TITLE"); ofs << buffer << endl;
  sprintf(buffer,"%3d  0DIST  0  0  0",mol.NumAtoms()); ofs << buffer << endl;
  sprintf(buffer,"  EL"); ofs << buffer << endl;
  sprintf(buffer,"0.,0.,0., %s",tmptype); ofs << buffer << endl;
  for (i = 2; i <= mol.NumAtoms(); i++)
  {
    strcpy(tmptype,etab.GetSymbol(mol.GetAtom(i)->GetAtomicNum()));

    if (vit[i]->_tor < 0.0) vit[i]->_tor += 360.0;
    sprintf(buffer,"%2d,%d,%2s%7.3f,%7.3f,%7.3f",
	    vit[i]->_a->GetIdx(),i,tmptype,
	    vit[i]->_dst,
	    vit[i]->_ang,
	    vit[i]->_tor);
    ofs << buffer << endl;
  }

  vector<OBInternalCoord*>::iterator j;
  for (j = vit.begin();j != vit.end();j++) 
    if (*j)
      {
	delete *j;
	*j = NULL;
      }

  return(true);
}

} //namespace OpenBabel

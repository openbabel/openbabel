/**********************************************************************
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include "mol.h"
#include "version.h"
#include "typer.h"
#include "resdata.h"

namespace OpenBabel {

extern OBAtomTyper atomtyper;

class OBResidueData : public OBGlobalDataBase
{
    int                                _resnum;
    vector<string>                     _resname;
    vector<vector<string> >            _resatoms;
    vector<vector<pair<string,int> > > _resbonds;

	//variables used only for parsing resdata.txt
    vector<string>                     _vatmtmp;
    vector<pair<string,int> >          _vtmp;
public:
    OBResidueData();
    bool SetResName(string &);
    int  LookupBO(string &);
    int  LookupBO(string &,string&);
    bool LookupType(string &,string&,int&);
    bool AssignBonds(OBMol &,OBBitVec &);
	void ParseLine(char*);
};

static bool ParseAtomRecord(char *, OBMol &,int);
static bool ParseConectRecord(char *,OBMol &);

static OBResidueData resdat;

bool ReadPDB(istream &ifs,OBMol &mol,char *title)
{
  resdat.Init();
  int chainNum = 1;
  char buffer[BUFF_SIZE];
  OBBitVec bs;

  mol.BeginModify();
  while (ifs.getline(buffer,BUFF_SIZE) && !EQn(buffer,"END",3))
  {
    if (EQn(buffer,"TER",3)) chainNum++;
    if (EQn(buffer,"ATOM",4) || EQn(buffer,"HETATM",6))
      {
	ParseAtomRecord(buffer,mol,chainNum);
	if (EQn(buffer,"ATOM",4))     
	  bs.SetBitOn(mol.NumAtoms());
      }

    if (EQn(buffer,"CONECT",6))       
      ParseConectRecord(buffer,mol);
  }
  
  resdat.AssignBonds(mol,bs);
  /*assign hetatm bonds based on distance*/
  mol.ConnectTheDots();

  mol.EndModify();
  mol.SetAtomTypesPerceived();
  atomtyper.AssignImplicitValence(mol);

  if (!mol.NumAtoms()) return(false);
  return(true);
}

bool ReadTerTermPDB(istream &ifs,OBMol &mol,char *title)
{
  resdat.Init();
  int chainNum = 1;
  char buffer[BUFF_SIZE];
  OBBitVec bs;

  mol.BeginModify();
  while (ifs.getline(buffer,BUFF_SIZE) && !EQn(buffer,"END",3) && !EQn(buffer,"TER",3))
  {
    if (EQn(buffer,"ATOM",4) || EQn(buffer,"HETATM",6))
      {
	ParseAtomRecord(buffer,mol,chainNum);
	if (EQn(buffer,"ATOM",4))     
	  bs.SetBitOn(mol.NumAtoms());
      }

    if (EQn(buffer,"CONECT",6))       
      ParseConectRecord(buffer,mol);
  }
  
  resdat.AssignBonds(mol,bs);
  /*assign hetatm bonds based on distance*/
  mol.ConnectTheDots();

  mol.EndModify();
  mol.SetAtomTypesPerceived();
  atomtyper.AssignImplicitValence(mol);

  if (!mol.NumAtoms()) return(false);
  return(true);
}

bool ReadPDB(vector<string> &vpdb,OBMol &mol,char *title)
{
  resdat.Init();
  int chainNum = 1;
  char buffer[BUFF_SIZE];
  vector<string>::iterator i;
  OBBitVec bs;

  mol.BeginModify();
  
  for (i = vpdb.begin();i != vpdb.end();i++)
    {
      strcpy(buffer,i->c_str());
      if (EQn(buffer,"END",3)) break;

      if (EQn(buffer,"TER",3)) chainNum++;
      if (EQn(buffer,"ATOM",4) || EQn(buffer,"HETATM",6))
	{
	  ParseAtomRecord(buffer,mol,chainNum);
	  if (EQn(buffer,"ATOM",4))     
	    bs.SetBitOn(mol.NumAtoms());
	}

      if (EQn(buffer,"CONECT",6))       
	ParseConectRecord(buffer,mol);
    }
  
  resdat.AssignBonds(mol,bs);
  /*assign hetatm bonds based on distance*/

  mol.EndModify();

  if (!mol.NumAtoms()) return(false);
  return(true);
}

bool WriteDelphiPDB(ostream &ofs,OBMol &mol)
{
  OBAtom *atom;
  vector<OBNodeBase*>::iterator i;
  char buffer[BUFF_SIZE];
  
  for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
    {
      sprintf(buffer,"ATOM   %4d %4s %3s    %1d    %8.3f%8.3f%8.3f%6.2f %6.3f",
	      atom->GetIdx(),
	      etab.GetSymbol(atom->GetAtomicNum()),
	      "UNK",0,
	      atom->GetX(),atom->GetY(),atom->GetZ(),
	      etab.GetVdwRad(atom->GetAtomicNum()),
	      atom->GetPartialCharge());
      ofs << buffer << endl;
    }

  int k,bo,count;
  int bond[10];
  OBAtom *nbr;
  vector<OBEdgeBase*>::iterator j;
  for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
    {
      count=0;
      memset(bond,'\0',sizeof(int)*5);
      bond[count++] = atom->GetIdx();
      for (nbr = atom->BeginNbrAtom(j);nbr;nbr = atom->NextNbrAtom(j))
	{
	  bo = ((*j)->GetBO() < 4) ?(*j)->GetBO():1;
	  for (k = 0;k < bo;k++)
	    bond[count++] = nbr->GetIdx();
	}
      sprintf(buffer,"CONECT  %3d  %3d  %3d  %3d  %3d",
	      bond[0],bond[1],bond[2],bond[3],bond[4]);

      ofs << buffer << endl;
    }

  ofs << "TER" << endl;
  return(true);
}

static bool ParseAtomRecord(char *buffer, OBMol &mol,int chainNum)
/* ATOMFORMAT "(i5,1x,a4,a1,a3,1x,a1,i4,a1,3x,3f8.3,2f6.2,1x,i3)" */
{
  string sbuf = &buffer[6];
  if (sbuf.size() < 48) return(false);

  bool hetatm = (EQn(buffer,"HETATM",6)) ? true : false;

  /* serial number */
  string serno = sbuf.substr(0,5);
  //SerialNum(the_atom) = atoi(tmp_str);

  /* atom name */
  string atmid = sbuf.substr(6,4);

  //trim spaces on the right and left sides
  while (!atmid.empty() && atmid[0] == ' ') 
    atmid = atmid.substr(1,atmid.size()-1);

  while (!atmid.empty() && atmid[atmid.size()-1] == ' ') 
    atmid = atmid.substr(0,atmid.size()-1);

  /* residue name */
  
  string resname = sbuf.substr(11,3);
  if (resname == "   ") 
    resname = "UNK";
  else
  {
    while (!resname.empty() && resname[0] == ' ') 
      resname = resname.substr(1,resname.size()-1);

    while (!resname.empty() && resname[resname.size()-1] == ' ') 
      resname = resname.substr(0,resname.size()-1);
  }  

  /* residue sequence number */

  string resnum = sbuf.substr(16,4);

  /* X, Y, Z */
  string xstr = sbuf.substr(24,8);
  string ystr = sbuf.substr(32,8);
  string zstr = sbuf.substr(40,8);

  string type;

  if (EQn(buffer,"ATOM",4))
  {
    type = atmid.substr(0,1);
    if (isdigit(type[0]))
       type = atmid.substr(1,1);

    if (resname.substr(0,2) == "AS" || resname[0] == 'N')
    {
      if (atmid == "AD1") type = "O";
      if (atmid == "AD2") type = "N";
    }
    if (resname.substr(0,3) == "HIS" || resname[0] == 'H')
    {
      if (atmid == "AD1" || atmid == "AE2") type = "N";
      if (atmid == "AE1" || atmid == "AD2") type = "C";
    }
    if (resname.substr(0,2) == "GL" || resname[0] == 'Q')
    {
      if (atmid == "AE1") type = "O";
      if (atmid == "AE2") type = "N";
    }
  }
  else //must be hetatm record
  {
    if (isalpha(atmid[0])) type = atmid.substr(0,1);
    else                   type = atmid.substr(1,1);
    if (atmid == resname)
      {
	type = atmid;
	if (type.size() == 2) type[1] = tolower(type[1]);
      }
    else
    if (resname == "ADR" || resname == "COA" || resname == "FAD" ||
	resname == "GPG" || resname == "NAD" || resname == "NAL" ||
	resname == "NDP")
      {
	if (type.size() > 1)
	  type = type.substr(0,1);
	//type.erase(1,type.size()-1);
      }
    else
      if (isdigit(type[0]))
	{
	  type = type.substr(1,1);
	  //type.erase(0,1);
	  //if (type.size() > 1) type.erase(1,type.size()-1);
	}
      else
	if (type.size() > 1 && isdigit(type[1]))
	  type = type.substr(0,1);
    //type.erase(1,1);
	else
	  if (type.size() > 1 && isalpha(type[1]) && isupper(type[1]))
	    type[1] = tolower(type[1]);

  }

  OBAtom atom;
  Vector v(atof(xstr.c_str()),atof(ystr.c_str()),atof(zstr.c_str()));
  atom.SetVector(v);
  
  atom.SetAtomicNum(etab.GetAtomicNum(type.c_str()));
  atom.SetType(type);

  int        rnum = atoi(resnum.c_str());
  OBResidue *res  = (mol.NumResidues() > 0) ? mol.GetResidue(mol.NumResidues()-1) : NULL;
  if (res == NULL || res->GetName() != resname || res->GetNum() != rnum)
  {
      vector<OBResidue*>::iterator ri;
      for (res = mol.BeginResidue(ri) ; res ; res = mol.NextResidue(ri))
          if (res->GetName() == resname && res->GetNum() == rnum)
              break;

      if (res == NULL)
      {
          res = mol.NewResidue();
          res->SetChainNum(chainNum);
          res->SetName(resname);
          res->SetNum(rnum);
      }
  }

  if (!mol.AddAtom(atom)) 
      return(false);
  else
  {
      OBAtom *atom = mol.GetAtom(mol.NumAtoms());

      res->AddAtom(atom);
      res->SetSerialNum(atom, atoi(serno.c_str()));
      res->SetAtomID(atom, atmid);
      res->SetHetAtom(atom, hetatm);

      return(true);
  }
}

static bool ParseConectRecord(char *buffer,OBMol &mol)
{
  vector<string> vs;
  
  buffer[70] = '\0';
  tokenize(vs,buffer);
  if (vs.empty()) return(false);
  vs.erase(vs.begin());
  int con1,con2,con3,con4;
  con1 = con2 = con3 = con4 = 0;
  int start = atoi(vs[0].c_str());

  if (vs.size() > 1) con1 = atoi(vs[1].c_str());
  if (vs.size() > 2) con2 = atoi(vs[2].c_str());
  if (vs.size() > 3) con3 = atoi(vs[3].c_str());
  if (vs.size() > 4) con4 = atoi(vs[4].c_str());
  if (!con1) return(false);

  OBAtom *a1,*a2;
  OBResidue *r1,*r2;
  vector<OBNodeBase*>::iterator i,j;
  for (a1 = mol.BeginAtom(i);a1;a1 = mol.NextAtom(i))
    {
      r1 = a1->GetResidue();
      if (r1->GetSerialNum(a1) == start)
	  for (a2 = mol.BeginAtom(j);a2;a2 = mol.NextAtom(j))
	    {
	      r2 = a2->GetResidue();
	      if (con1 && r2->GetSerialNum(a2) == con1) 
		mol.AddBond(a1->GetIdx(),a2->GetIdx(),1);
	      if (con2 && r2->GetSerialNum(a2) == con2) 
		mol.AddBond(a1->GetIdx(),a2->GetIdx(),1);
	      if (con3 && r2->GetSerialNum(a2) == con3) 
		mol.AddBond(a1->GetIdx(),a2->GetIdx(),1);
	      if (con4 && r2->GetSerialNum(a2) == con4) 
		mol.AddBond(a1->GetIdx(),a2->GetIdx(),1);
	    }
    }

  return(true);
}

OBResidueData::OBResidueData()
{
  _init = false;
  _dir = DATADIR;
  _envvar = "BABEL_DATADIR";
  _filename = "residue.txt";
  _subdir = "data";
  _dataptr = ResidueData;
}

bool OBResidueData::AssignBonds(OBMol &mol,OBBitVec &bv)
{
  OBAtom *a1,*a2;
  OBResidue *r1,*r2;
  vector<OBNodeBase*>::iterator i,j;
  
  //assign alpha peptide bonds
  for (a1 = mol.BeginAtom(i);a1;a1 = mol.NextAtom(i))
    if (bv.BitIsOn(a1->GetIdx()))
      {
	r1 = a1->GetResidue();
	if (!(r1->GetAtomID(a1) == "C")) continue;
	for (j=i,a2 = mol.NextAtom(j);a2;a2 = mol.NextAtom(j))
	  {
	    r2 = a2->GetResidue();
	    if (!(r2->GetAtomID(a2) == "N")) continue;
	    if (r1->GetNum() < r2->GetNum()-1) break;
	    if (r1->GetChainNum() == r2->GetChainNum())
	      {
		mol.AddBond(a1->GetIdx(),a2->GetIdx(),1);
		break;
	      }
	  }
      }

  Vector v;
  int bo,skipres=0;
  string rname = "";
  //assign other residue bonds
  for (a1 = mol.BeginAtom(i);a1;a1 = mol.NextAtom(i))
    {
      r1 = a1->GetResidue();
      if (skipres && r1->GetNum() == skipres) continue;

      if (r1->GetName() != rname)
	{
	  skipres = SetResName(r1->GetName()) ? 0 : r1->GetNum(); 
	  rname = r1->GetName();
	}
      //assign bonds for each atom
      for (j=i,a2 = mol.NextAtom(j);a2;a2 = mol.NextAtom(j))
	{
	  r2 = a2->GetResidue();
	  if (r1->GetNum() != r2->GetNum()) break;
	  if (r1->GetName() != r2->GetName()) break;

	  if ((bo = LookupBO(r1->GetAtomID(a1),r2->GetAtomID(a2))))
	    {
	      v = a1->GetVector() - a2->GetVector();
	      if (v.length_2() < 3.5) //float check by distance
		  mol.AddBond(a1->GetIdx(),a2->GetIdx(),bo);
	    }
	}
    }

  int hyb;
  string type;

  //types and hybridization
  for (a1 = mol.BeginAtom(i);a1;a1 = mol.NextAtom(i))
    {
      if (a1->IsOxygen() && !a1->GetValence())
	{
	  a1->SetType("O3");
	  continue;
	}

      if (a1->IsHydrogen())
	{
	  a1->SetType("H");
	  continue;
	}

      r1 = a1->GetResidue();
      if (skipres && r1->GetNum() == skipres) continue;

      if (r1->GetName() != rname)
	{
	  skipres = SetResName(r1->GetName()) ? 0 : r1->GetNum(); 
	  rname = r1->GetName();
	}

      //***valence rule for O-
      if (a1->IsOxygen() && a1->GetValence() == 1)
	{
	  OBBond *bond;
	  bond = (OBBond*)*(a1->BeginBonds());
	  if (bond->GetBO() == 2)
	    {
	      a1->SetType("O2"); a1->SetHyb(2);
	    }
	  if (bond->GetBO() == 1)
	    {
	      a1->SetType("O-"); a1->SetHyb(3);
	      a1->SetFormalCharge(-1);
	    }
	}
      else
	if (LookupType(r1->GetAtomID(a1),type,hyb))
	  {
	    a1->SetType(type);
	    a1->SetHyb(hyb);
	  }
	else // try to figure it out by bond order ???
	  {
	  }
    }

  return(true);
}

void OBResidueData::ParseLine(char *buffer)
{
   int bo;
    string s;
    vector<string> vs;

	tokenize(vs,buffer);
	if (!vs.empty())
	  {
	    if (vs[0] == "BOND")
	      {
			s = (vs[1] < vs[2]) ? vs[1] + " " + vs[2] : 
			vs[2] + " " + vs[1];
			bo = atoi(vs[3].c_str());
			_vtmp.push_back(pair<string,int> (s,bo));
	      }

	    if (vs[0] == "ATOM" && vs.size() == 4)
	      {
			_vatmtmp.push_back(vs[1]);
			_vatmtmp.push_back(vs[2]);
			_vatmtmp.push_back(vs[3]);
	      }

	    if (vs[0] == "RES")	_resname.push_back(vs[1]);

	    if (vs[0]== "END")
	      {
			_resatoms.push_back(_vatmtmp);
			_resbonds.push_back(_vtmp);
			_vtmp.clear();
			_vatmtmp.clear();
	      }
	  }
}

bool OBResidueData::SetResName(string &s)
{
  unsigned int i;
    for (i = 0;i < _resname.size();i++)
	if (_resname[i] == s)
	{
	    _resnum = i;
	    return(true);
	}

    _resnum = -1;
    return(false);
}

int OBResidueData::LookupBO(string &s)
{
    if (_resnum == -1) return(0);
    
    unsigned int i;
    for (i = 0;i < _resbonds[_resnum].size();i++)
	if (_resbonds[_resnum][i].first == s)
	    return(_resbonds[_resnum][i].second);

    return(0);
}

int OBResidueData::LookupBO(string &s1,string &s2)
{
    if (_resnum == -1) return(0);
    string s;

    s = (s1 < s2) ? s1 + " " + s2 : s2 + " " + s1;

    unsigned int i;
    for (i = 0;i < _resbonds[_resnum].size();i++)
	if (_resbonds[_resnum][i].first == s)
	    return(_resbonds[_resnum][i].second);

    return(0);
}

bool OBResidueData::LookupType(string &atmid,string &type,int &hyb)
{
    if (_resnum == -1) return(false);

    string s;
    vector<string>::iterator i;

    for (i = _resatoms[_resnum].begin();i != _resatoms[_resnum].end();i+=3)
      if (atmid == *i)
	{
	  i++;	  type = *i;
	  i++;	  hyb = atoi((*i).c_str());
	  return(true);
	}

    return(false);
}

bool ReadBox(vector<string> &vbox, OBMol &mol,char *)
{
  char buffer[BUFF_SIZE];
  vector<string> vs;
  vector<string>::iterator i,j;
  OBAtom atom;

  mol.BeginModify();
  
  for (i = vbox.begin();i != vbox.end();i++)
  {
    strcpy(buffer,i->c_str());

    if (EQn(buffer,"END",3)) break;
    if (EQn(buffer,"ATOM",4))
      {
	string sbuf = &buffer[6]; 
	/* X, Y, Z */
	string x = sbuf.substr(24,8);
	string y = sbuf.substr(32,8);
	string z = sbuf.substr(40,8);
	Vector v(atof(x.c_str()),atof(y.c_str()),atof(z.c_str()));
	atom.SetVector(v);
	if (!mol.AddAtom(atom)) return(false);
      }
    
    if (EQn(buffer,"CONECT",6))
      {
	tokenize(vs,buffer);
	if (!vs.empty() && vs.size() > 2)
	  for (j = vs.begin(),j+=2;j != vs.end();j++)
	    mol.AddBond(atoi(vs[1].c_str()),atoi((*j).c_str()),1);
      }
  }

  mol.EndModify();
  
  return(true);
}

bool WritePDB(ostream &ofs,OBMol &mol)
{
  unsigned int i;
  char buffer[BUFF_SIZE];
  char type_name[10], padded_name[10];
  char the_res[10];
  int res_num;

  sprintf(buffer,"HEADER    PROTEIN");
  ofs << buffer << endl;

  if (strlen(mol.GetTitle()) > 0) 
    sprintf(buffer,"COMPND    %s ",mol.GetTitle());
  else   sprintf(buffer,"COMPND    UNNAMED");
  ofs << buffer << endl;

  sprintf(buffer,"AUTHOR    GENERATED BY BABEL %s",BABEL_VERSION);
  ofs << buffer << endl;

  OBAtom *atom;
  OBResidue *res;
  for (i = 1; i <= mol.NumAtoms(); i++)
  {
    atom = mol.GetAtom(i);
    strcpy(type_name,etab.GetSymbol(atom->GetAtomicNum()));
    if (strlen(type_name) > 1) type_name[1] = toupper(type_name[1]);

    if (atom->HasResidue())
      {
	res = atom->GetResidue();
	strcpy(the_res,(char*)res->GetName().c_str());
	strcpy(type_name,(char*)res->GetAtomID(atom).c_str());
	res_num = res->GetNum();
      }
    else
      {
	strcpy(the_res,"UNK");
	sprintf(padded_name,"%2s",type_name);
	strcpy(type_name,padded_name);
	res_num = 1;
      }
    sprintf(buffer,"ATOM   %4d  %-4s%3s%3s%3d   %9.3f%8.3f%8.3f  1.00  0.00 \n",
	    i,
	    type_name,
	    the_res,
	    "",
	    res_num,
	    atom->GetX(),
	    atom->GetY(),
	    atom->GetZ());
    ofs << buffer;
  }

  OBAtom *nbr;
  vector<OBEdgeBase*>::iterator k;
  for (i = 1; i <= mol.NumAtoms(); i ++)
  {
    atom = mol.GetAtom(i);
    if (atom->GetValence())
      {
	sprintf(buffer,"CONECT %5d",i);
	ofs << buffer;
	for (nbr = atom->BeginNbrAtom(k);nbr;nbr = atom->NextNbrAtom(k))
	  {
	    sprintf(buffer,"%5d",nbr->GetIdx());
	    ofs << buffer;
	  }
	ofs << endl;
      }
  }
  sprintf(buffer,"MASTER        0    0    0    0    0    0    0    0 ");
  ofs << buffer;
  sprintf(buffer,"%4d    0 %4d    0",mol.NumAtoms(),mol.NumAtoms());          
  ofs << buffer << endl;
  sprintf(buffer,"END"); ofs << buffer << endl;
  return(true);
}

}

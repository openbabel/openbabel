/**********************************************************************
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2003 Geoffrey R. Hutchison

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

#if HAVE_CONFIG_H
#include "babelconfig.h"
#endif

#include "mol.h"
#include "typer.h"
#include "resdata.h"

#include <vector>
#include <map>

using namespace std;

namespace OpenBabel {

extern OBAtomTyper atomtyper;

//! \brief Table of common protein residues in PDB files
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
    bool SetResName(const string &);
    int  LookupBO(const string &);
    int  LookupBO(const string &, const string&);
    bool LookupType(const string &,string&,int&);
    bool AssignBonds(OBMol &,OBBitVec &);
    void ParseLine(const char*);
};

class OBSerialNums : public OBGenericData
{
public:

        OBSerialNums()                                           
        { _attr = "obSerialNums"; _type = obSerialNums; }
        OBSerialNums(const OBSerialNums &cp) : OBGenericData(cp) 
        { serialMap = cp.serialMap;                }

        map<int,OBAtom*> &GetData()                     { return serialMap; }
        void              SetData(map<int,OBAtom*> &sm) { serialMap = sm; }

protected:
        map<int, OBAtom*> serialMap;
};

static bool ParseAtomRecord(char *, OBMol &,int);
static bool ParseConectRecord(char *,OBMol &);

static OBResidueData resdat;

bool ReadPDB(istream &ifs,OBMol &mol,const char *title)
{
  resdat.Init();
  int chainNum = 1;
  char buffer[BUFF_SIZE];
  OBBitVec bs;

  mol.SetTitle(title);

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
  mol.EndModify();
  
  mol.ConnectTheDots();

  if (mol.NumAtoms() < 250) // Minimize time required on real proteins
    mol.PerceiveBondOrders();
  mol.SetAtomTypesPerceived();
  atomtyper.AssignImplicitValence(mol);

  if (!mol.NumAtoms()) return(false);
  return(true);
}

bool ReadTerTermPDB(istream &ifs,OBMol &mol,const char *title)
{
  resdat.Init();
  int chainNum = 1;
  char buffer[BUFF_SIZE];
  OBBitVec bs;

  mol.SetTitle(title);
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
  if (mol.NumAtoms() < 250) // Minimize time required on real proteins
    mol.PerceiveBondOrders();
  mol.SetAtomTypesPerceived();
  atomtyper.AssignImplicitValence(mol);

  if (!mol.NumAtoms()) return(false);
  return(true);
}

bool ReadPDB(vector<string> &vpdb,OBMol &mol,const char *title)
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
      sprintf(buffer,"ATOM  %5d %-4s %-3s  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f",
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
      sprintf(buffer,"CONECT%5d%5d%5d%5d%5d",
	      bond[0],bond[1],bond[2],bond[3],bond[4]);
      ofs << buffer << "                                       " << endl;
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

  /* element */
  string element;
  if (sbuf.size() > 71)
    element = sbuf.substr(70,2);
  else
    element = "  ";

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
    type = atmid.substr(0,2);
    if (isdigit(type[0]))
      type = atmid.substr(1,1);
    else if (sbuf[6] == ' ') // one-character element
      type = atmid.substr(0,1);

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
    if (isalpha(element[1]) && (isalpha(element[0]) || (element[0] == ' ')))
      {
        if (isalpha(element[0])) type = element.substr(0,2);
        else type = element.substr(1,1);
	if (type.size() == 2) type[1] = tolower(type[1]);
      }
    else
      {   
	if (isalpha(atmid[0])) type = atmid.substr(0,2);
	else if (atmid[0] == ' ') type = atmid.substr(1,1); // one char element
	else                   type = atmid.substr(1,2);
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

  }

  OBAtom atom;
  vector3 v(atof(xstr.c_str()),atof(ystr.c_str()),atof(zstr.c_str()));
  atom.SetVector(v);
  
  atom.SetAtomicNum(etab.GetAtomicNum(type.c_str()));
  atom.SetType(type);

  int        rnum = atoi(resnum.c_str());
  OBResidue *res  = (mol.NumResidues() > 0) ? mol.GetResidue(mol.NumResidues()-1) : NULL;
  if (res == NULL || res->GetName() != resname || static_cast<int>(res->GetNum()) != rnum)
  {
      vector<OBResidue*>::iterator ri;
      for (res = mol.BeginResidue(ri) ; res ; res = mol.NextResidue(ri))
          if (res->GetName() == resname && static_cast<int>(res->GetNum()) == rnum)
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

//! Utility function to read a 5-digit integer starting from a specified column
/*! This function reads a 5-digit integer, starting from column
  columnAsSpecifiedInPDB from the buffer, converts it to a long
  integer, and returns either false or true, if the conversion was
  successful or not. If the conversion was not successful, the target
  is set to a random value.

  For instance, the PDB Format Description for a CONECT record specifies

  COLUMNS        DATA TYPE        FIELD           DEFINITION
  ---------------------------------------------------------------------------------
  1 -  6         Record name      "CONECT"
  7 - 11         Integer          serial          Atom serial number
  ...

  To read the Atom serial number, you would call

  long int target;
  if ( readIntegerFromRecord(buffer, 7, &target) == false ) {
    cerr << "Could not parse" << endl;
  }
  
  This function does not check the length of the buffer, or
  strlen(buffer). If the buffer is not long enough => SEGFAULT. 
 */
static bool readIntegerFromRecord(char *buffer, unsigned int columnAsSpecifiedInPDB, long int *target)
{
  char integerBuffer[6];
  integerBuffer[5] = 0;

  strncpy(integerBuffer, buffer+columnAsSpecifiedInPDB-1, 5);

  char *errorCheckingEndPtr;
  *target = strtol(integerBuffer, &errorCheckingEndPtr, 10);
  if (integerBuffer == errorCheckingEndPtr)
    return(false);
  return(true);
}

//! Read a CONECT record
/*! This function reads a CONECT record, as specified
  http://www.rcsb.org/pdb/docs/format/pdbguide2.2/guide2.2_frame.html,
  in short:

  COLUMNS         DATA TYPE        FIELD           DEFINITION
  ---------------------------------------------------------------------------------
   1 -  6         Record name      "CONECT"
   7 - 11         Integer          serial          Atom serial number
  12 - 16         Integer          serial          Serial number of bonded atom
  17 - 21         Integer          serial          Serial number of bonded atom
  22 - 26         Integer          serial          Serial number of bonded atom
  27 - 31         Integer          serial          Serial number of bonded atom
  32 - 36         Integer          serial          Serial number of hydrogen bonded atom
  37 - 41         Integer          serial          Serial number of hydrogen bonded atom
  42 - 46         Integer          serial          Serial number of salt bridged atom
  47 - 51         Integer          serial          Serial number of hydrogen bonded atom
  52 - 56         Integer          serial          Serial number of hydrogen bonded atom
  57 - 61         Integer          serial          Serial number of salt bridged atom

  Hydrogen bonds and salt bridges are ignored. --Stefan Kebekus.
*/

static bool ParseConectRecord(char *buffer,OBMol &mol)
{
  // Setup strings and string buffers
  buffer[70] = '\0';
  if (strlen(buffer) < 70) {
    cerr << "WARNING: Problems reading a PDB file, method 'static bool ParseConectRecord(char *, OBMol &)'" << endl
	 << "  Problems reading a CONECT record." << endl
	 << "  OpenBabel found the line '" << buffer << "'" << endl
	 << "  According to the PDB specification (http://www.rcsb.org/pdb/docs/format/pdbguide2.2/guide2.2_frame.html)," << endl
	 << "  the record should have 70 columns, but OpenBabel found " << strlen(buffer) << " columns." << endl
	 << "  THIS CONECT RECORD WILL BE IGNORED." << endl;
    return(false);
  }

  // Serial number of the first atom, read from column 7-11 of the
  // connect record, to which the other atoms connect to.
  long int startAtomSerialNumber;
  if (readIntegerFromRecord(buffer, 7, &startAtomSerialNumber) == false) {
    cerr << "WARNING: Problems reading a PDB file, method 'static bool ParseConectRecord(char *, OBMol &)'" << endl
	 << "  Problems reading a CONECT record." << endl
	 << "  OpenBabel found the line '" << buffer << "'" << endl
	 << "  According to the PDB specification (http://www.rcsb.org/pdb/docs/format/pdbguide2.2/guide2.2_frame.html)," << endl
	 << "  columns 7--11 should contain the serial number of an atom, but OpenBabel was not able" << endl
	 << "  to interpret these columns. " << endl
	 << "  THIS CONECT RECORD WILL BE IGNORED." << endl;
    return(false);
  }

  // Find a pointer to the first atom.
  OBAtom *firstAtom = 0L;
  vector<OBNodeBase*>::iterator i;
  for (OBAtom *a1 = mol.BeginAtom(i);a1;a1 = mol.NextAtom(i)) 
    if (static_cast<long int>(a1->GetResidue()->GetSerialNum(a1)) == startAtomSerialNumber) {
      firstAtom = a1;
      break;
    }
  if (firstAtom == 0L) {
    cerr << "WARNING: Problems reading a PDB file, method 'static bool ParseConectRecord(char *, OBMol &)'" << endl
	 << "  Problems reading a CONECT record." << endl
	 << "  OpenBabel found the line '" << buffer << "'" << endl
	 << "  According to the PDB specification (http://www.rcsb.org/pdb/docs/format/pdbguide2.2/guide2.2_frame.html)," << endl
	 << "  columns 7--11 should contain the serial number of an atom, but OpenBabel was not able" << endl
	 << "  to find an atom with this serial number. " << endl
	 << "  THIS CONECT RECORD WILL BE IGNORED." << endl;
    return(false);
  }
  
  // Serial numbers of the atoms which bind to firstAtom, read from
  // columns 12-16, 17-21, 22-27 and 27-31 of the connect record. Note
  // that we reserve space for 5 integers, but read only four of
  // them. This is to simplify the determination of the bond order;
  // see below.
  long int boundedAtomsSerialNumbers[5]  = {0,0,0,0,0};
  // Bools which tell us which of the serial numbers in
  // boundedAtomsSerialNumbers are read from the file, and which are
  // invalid
  bool boundedAtomsSerialNumbersValid[5] = {false, false, false, false, false};
  
  // Now read the serial numbers. If the first serial number is not
  // present, this connect record probably contains only hydrogen
  // bonds and salt bridges, which we ignore. In that case, we just
  // exit gracefully.
  boundedAtomsSerialNumbersValid[0] = readIntegerFromRecord(buffer, 12, boundedAtomsSerialNumbers+0);
  if (boundedAtomsSerialNumbersValid[0] == false) 
    return(true);
  boundedAtomsSerialNumbersValid[1] = readIntegerFromRecord(buffer, 17, boundedAtomsSerialNumbers+1);
  boundedAtomsSerialNumbersValid[2] = readIntegerFromRecord(buffer, 22, boundedAtomsSerialNumbers+2);
  boundedAtomsSerialNumbersValid[3] = readIntegerFromRecord(buffer, 27, boundedAtomsSerialNumbers+3);
  
  // Now iterate over the VALID boundedAtomsSerialNumbers and connect
  // the atoms.
  for(unsigned int k=0; boundedAtomsSerialNumbersValid[k]; k++) {
    // Find atom that is connected to, write an error message 
    OBAtom *connectedAtom = 0L;
    for (OBAtom *a1 = mol.BeginAtom(i);a1;a1 = mol.NextAtom(i)) 
      if (static_cast<long int>(a1->GetResidue()->GetSerialNum(a1)) == boundedAtomsSerialNumbers[k]) {
	connectedAtom = a1;
	break;
      }
    if (connectedAtom == 0L) {
      cerr << "WARNING: Problems reading a PDB file, method 'static bool ParseConectRecord(char *, OBMol &)'" << endl
	   << "  Problems reading a CONECT record." << endl
	   << "  OpenBabel found the line '" << buffer << "'" << endl
	   << "  According to the PDB specification (http://www.rcsb.org/pdb/docs/format/pdbguide2.2/guide2.2_frame.html)," << endl
	   << "  OpenBabel should connect atoms with serial #" << startAtomSerialNumber << " and #" << boundedAtomsSerialNumbers[k] << endl
	   << "  However, OpenBabel was not able to find an atom with serial #" << boundedAtomsSerialNumbers[k] << "." << endl
	   << "  OpenBabel will proceed, and disregard this particular connection." << endl;
      break;
    }
    
    // Figure the bond order
    unsigned char order = 0;
    while(boundedAtomsSerialNumbersValid[k+order+1] && (boundedAtomsSerialNumbers[k+order] == boundedAtomsSerialNumbers[k+order+1]))
      order++;
    k += order;
    
    // Generate the bond
    mol.AddBond(firstAtom->GetIdx(), connectedAtom->GetIdx(), order+1);
  }
  return(true);
}

OBResidueData::OBResidueData()
{
  _init = false;
  _dir = BABEL_DATADIR;
  _envvar = "BABEL_DATADIR";
  _filename = "resdata.txt";
  _subdir = "data";
  _dataptr = ResidueData;
}

bool OBResidueData::AssignBonds(OBMol &mol,OBBitVec &bv)
{
  OBAtom *a1,*a2;
  OBResidue *r1,*r2;
  vector<OBNodeBase*>::iterator i,j;
  vector3 v;
  
  int bo;
  unsigned int skipres=0;
  string rname = "";
  //assign residue bonds
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
	      if (v.length_2() < 3.5) //check by distance
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

void OBResidueData::ParseLine(const char *buffer)
{
   int bo;
    string s;
    vector<string> vs;

    if (buffer[0] == '#')
      return;

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

bool OBResidueData::SetResName(const string &s)
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

int OBResidueData::LookupBO(const string &s)
{
    if (_resnum == -1) return(0);
    
    unsigned int i;
    for (i = 0;i < _resbonds[_resnum].size();i++)
	if (_resbonds[_resnum][i].first == s)
	    return(_resbonds[_resnum][i].second);

    return(0);
}

int OBResidueData::LookupBO(const string &s1, const string &s2)
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

bool OBResidueData::LookupType(const string &atmid,string &type,int &hyb)
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

bool ReadBox(vector<string> &vbox, OBMol &mol,const char *)
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
	vector3 v(atof(x.c_str()),atof(y.c_str()),atof(z.c_str()));
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
  char *element_name;
  int res_num;
  bool het=true;

  //  sprintf(buffer,"HEADER    PROTEIN");
  //  ofs << buffer << endl;

  if (strlen(mol.GetTitle()) > 0) 
    sprintf(buffer,"COMPND    %s ",mol.GetTitle());
  else   sprintf(buffer,"COMPND    UNNAMED");
  ofs << buffer << endl;

  sprintf(buffer,"AUTHOR    GENERATED BY OPEN BABEL %s",BABEL_VERSION);
  ofs << buffer << endl;

  OBAtom *atom;
  OBResidue *res;
  for (i = 1; i <= mol.NumAtoms(); i++)
  {
    atom = mol.GetAtom(i);
    strcpy(type_name,etab.GetSymbol(atom->GetAtomicNum()));

    //two char. elements are on position 13 and 14 one char. start at 14 
    if (strlen(type_name) > 1)
      type_name[1] = toupper(type_name[1]);
    else
      {
	char tmp[10];
	strcpy(tmp, type_name);
	sprintf(type_name, " %-3s", tmp);
      }

    if (res = atom->GetResidue())
      {
        het = res->IsHetAtom(atom);
        snprintf(the_res,4,"%s",(char*)res->GetName().c_str());
        snprintf(type_name,5,"%s",(char*)res->GetAtomID(atom).c_str());

	//two char. elements are on position 13 and 14 one char. start at 14
	if (strlen(etab.GetSymbol(atom->GetAtomicNum())) == 1)
	{
	  if (strlen(type_name) < 4)
	  {
	    char tmp[10];
	    strcpy(tmp, type_name);
	    sprintf(padded_name," %-3s", tmp);
	    strncpy(type_name,padded_name,4);
	    type_name[4] = '\0';
	  }
	  else
	  {
	    type_name[4] = type_name[3]; type_name[3] = type_name[2];
	    type_name[2] = type_name[1]; type_name[1] = type_name[0];
	    type_name[0] = type_name[4]; type_name[4] = '\0';
	  }
	}
	res_num = res->GetNum();
      }
    else
      {
	strcpy(the_res,"UNK");
	sprintf(padded_name,"%s",type_name);
	strncpy(type_name,padded_name,4);
	type_name[4] = '\0';
	res_num = 1;
      }

    element_name = etab.GetSymbol(atom->GetAtomicNum());
    if (strlen(element_name) == 2)
      element_name[1] = toupper(element_name[1]);
    sprintf(buffer,"%s%5d %-4s %-3s  %4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s  \n",
            het?"HETATM":"ATOM  ",
	    i,
	    type_name,
	    the_res,
	    res_num,
	    atom->GetX(),
	    atom->GetY(),
	    atom->GetZ(),
	    element_name);
    ofs << buffer;
  }

  OBAtom *nbr;
  int count;
  vector<OBEdgeBase*>::iterator k;
  for (i = 1; i <= mol.NumAtoms(); i ++)
  {
    atom = mol.GetAtom(i);
    if (atom->GetValence() <= 4)
      {
	sprintf(buffer,"CONECT%5d", i);
	ofs << buffer;
	for (nbr = atom->BeginNbrAtom(k);nbr;nbr = atom->NextNbrAtom(k))
	  {
	    sprintf(buffer,"%5d", nbr->GetIdx());
	    ofs << buffer;
	  }
	for (count = 0; count < (4 - (int)atom->GetValence()); count++)
	  {
	    sprintf(buffer, "     ");
	    ofs << buffer;
	  }
	ofs << "                                       " << endl;
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

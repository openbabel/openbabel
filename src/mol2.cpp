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
#include "smi.h"

using namespace std;
namespace OpenBabel {

bool ReadMol2(istream &ifs,OBMol &mol,char *title)
{
  bool foundAtomLine = false;
  char buffer[BUFF_SIZE];
  char *comment = NULL;
  string str,str1;
  vector<string> vstr;

  mol.BeginModify();

  for (;;)
    {
      if (!ifs.getline(buffer,BUFF_SIZE)) return(false);
      if (EQn(buffer,"@<TRIPOS>MOLECULE",17)) break;
    }

  int lcount;
  int natoms,nbonds;
  for (lcount=0;;lcount++)
    {
      if (!ifs.getline(buffer,BUFF_SIZE)) return(false);
      if (EQn(buffer,"@<TRIPOS>ATOM",13)) {foundAtomLine = true; break;}
      
      if (lcount == 0)
	{
	  tokenize(vstr,buffer);
	  if (!vstr.empty()) mol.SetTitle(buffer);
	}
      else if (lcount == 1)
	sscanf(buffer,"%d%d",&natoms,&nbonds);
      else if (lcount == 4) //energy
	{
	  tokenize(vstr,buffer);
	  if (!vstr.empty() && vstr[0] == "Energy" && vstr.size() == 3)
	    mol.SetEnergy(atof(vstr[2].c_str()));
	}
      else if (lcount == 5) //comment
	{
	if (strlen(buffer) != 0)
	  {
	    comment = new char [strlen(buffer)+1];
	    strcpy(comment,buffer);
	  }
	}
    }
    
    if (!foundAtomLine)
      {
	mol.EndModify();
	mol.Clear();
	ThrowError("Unable to read Mol2 format file");
	return(false);
      }

    mol.ReserveAtoms(natoms);

    int i;
    Vector v;
    OBAtom atom;
    bool hasPartialCharges=false;
    float x,y,z,pcharge;
    char temp_type[BUFF_SIZE];

    ttab.SetFromType("SYB");
    for (i = 0;i < natoms;i++)
      {
	if (!ifs.getline(buffer,BUFF_SIZE)) return(false);
	sscanf(buffer," %*s %*s %f %f %f %s %*s %*s %f",
	       &x,&y,&z,
	       temp_type,
	       &pcharge);

	v.SetX(x);v.SetY(y);v.SetZ(z);
	atom.SetVector(v);

	str = temp_type;
	ttab.SetToType("ATN"); ttab.Translate(str1,str);
	atom.SetAtomicNum(atoi(str1.c_str()));
	ttab.SetToType("INT");	ttab.Translate(str1,str);
	atom.SetType(str1);
	atom.SetPartialCharge(pcharge);
	if (!mol.AddAtom(atom)) return(false);
	if (pcharge != 0.0f) hasPartialCharges = true;
      }

    for (;;)
      {
	if (!ifs.getline(buffer,BUFF_SIZE)) return(false);
	str = buffer;
	if (EQn(buffer,"@<TRIPOS>BOND",13)) break;
      }
    
    int start,end,order;
    for (i = 0; i < nbonds; i++) 
    {
      if (!ifs.getline(buffer,BUFF_SIZE)) return(false);

      sscanf(buffer,"%*d %d %d %s",&start,&end,temp_type);
      str = temp_type;
      order = 1;
      if (str == "ar" || str == "AR" || str == "Ar") order = 5;
      else
	if (str == "AM" || str == "am" || str == "Am") order = 1;
	else  order = atoi(str.c_str());

      mol.AddBond(start,end,order);
    }

    mol.EndModify();

    //must add generic data after end modify - otherwise it will be blown away
    if (comment)
	{
	  OBCommentData *cd = new OBCommentData;
	  cd->SetData(comment);
	  mol.SetData(cd);
	  delete [] comment;
	  comment = NULL;
	}
    if (hasPartialCharges)
      mol.SetPartialChargesPerceived();

    return(true);
}

bool WriteMol2(ostream &ofs,OBMol &mol,char *dimension)
{
  string str,str1;
  char buffer[BUFF_SIZE],label[BUFF_SIZE];
  char rnum[BUFF_SIZE],rlabel[BUFF_SIZE];

  ofs << "@<TRIPOS>MOLECULE" << endl;
  str = mol.GetTitle();
  if (str.empty()) ofs << "*****" << endl;
  else             ofs << str << endl;

  sprintf(buffer, " %d %d 0 0 0", mol.NumAtoms(),mol.NumBonds());
  ofs << buffer << endl;
  ofs << "SMALL" << endl;

  //if (mol.HasPartialChargesPerceived()) ofs << "GASTEIGER" << endl;
  //  else                                  ofs << "NO_CHARGES" << endl;

  ofs << "GASTEIGER" << endl;
  ofs << "Energy = " << mol.GetEnergy() << endl;
  
  if (mol.HasData("Comment"))
    {
      OBCommentData *cd = (OBCommentData*)mol.GetData(obCommentData);
      ofs << cd->GetData();
    }

  ofs << endl;
  ofs << "@<TRIPOS>ATOM" << endl;

  ttab.SetFromType("INT");ttab.SetToType("SYB");

  OBAtom *atom;
  vector<OBNodeBase*>::iterator i;
  vector<int> labelcount;labelcount.resize(105); //Number of elements
  for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
    {
      sprintf(label,"%s%d",
	      etab.GetSymbol(atom->GetAtomicNum()),
	      ++labelcount[atom->GetAtomicNum()]);

      str = atom->GetType();

      ttab.Translate(str1,str);

      strcpy(rlabel,"<1>");
      strcpy(rnum,"1");

      sprintf(buffer,"%7d%1s%-6s%12.4f%10.4f%10.4f%1s%-5s%4s%1s %-8s%10.4f",
	    atom->GetIdx(),"",label, 
	      atom->GetX(),atom->GetY(),atom->GetZ(),
	      "",str1.c_str(),
	      rnum,"",rlabel,
	      atom->GetPartialCharge());
      ofs << buffer << endl;
    }

  ofs << "@<TRIPOS>BOND" << endl;
  OBBond *bond;
  vector<OBEdgeBase*>::iterator j;
  for (bond = mol.BeginBond(j);bond;bond = mol.NextBond(j))
    {
      if (bond->GetBO() == 5) strcpy(label,"ar");
      else if (bond->IsAmide()) strcpy(label,"am");
      else sprintf(label,"%d",bond->GetBO());
      sprintf(buffer, "%6d%6d%6d%3s%2s", 
	      bond->GetIdx()+1,bond->GetBeginAtomIdx(),bond->GetEndAtomIdx(),
	      "",label);
      ofs << buffer << endl;
    }
  ofs << endl;

  return(true);
}
bool WriteSmiOrderedMol2(ostream &ofs,OBMol &mol,char *dimension)
{
  string str,str1;
  char buffer[BUFF_SIZE],label[BUFF_SIZE];
  char rnum[BUFF_SIZE],rlabel[BUFF_SIZE];

  ofs << "@<TRIPOS>MOLECULE" << endl;
  str = mol.GetTitle();
  if (str.empty()) ofs << "*****" << endl;
  else             ofs << str << endl;

  sprintf(buffer, " %d %d 0 0 0", mol.NumAtoms(),mol.NumBonds());
  ofs << buffer << endl;
  ofs << "SMALL" << endl;

  //if (mol.HasPartialChargesPerceived()) ofs << "GASTEIGER" << endl;
  //  else                                  ofs << "NO_CHARGES" << endl;

  ofs << "GASTEIGER" << endl;
  ofs << "Energy = " << mol.GetEnergy() << endl;
  
  if (mol.HasData("Comment"))
    ofs << (char*)mol.GetData("Comment");

  ofs << endl;
  ofs << "@<TRIPOS>ATOM" << endl;

  ttab.SetFromType("INT");ttab.SetToType("SYB");

	//get smiles order
	OBMol2Smi m2s;
	char smibuffer[BUFF_SIZE];
	m2s.Init();
	m2s.CorrectAromaticAmineCharge(mol);
	m2s.CreateSmiString(mol,smibuffer);
	vector<int>::iterator idx;
	vector<int> smiorder;
	int *backmap = new int[mol.NumAtoms()];
	smiorder = m2s.GetOutputOrder();
	int ct;

  OBAtom *atom;
//  vector<OBNodeBase*>::iterator i;
  vector<int> labelcount;labelcount.resize(105); //Number of elements
//  for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
	for(ct = 1,idx = smiorder.begin();idx != smiorder.end();idx++,ct++) //loop over smiles order
    {
cerr << (*idx) << " ";
			atom = mol.GetAtom(*idx); //set atom for .mol2 files
			backmap[atom->GetIdx()] = ct;

      sprintf(label,"%s%d",
	      etab.GetSymbol(atom->GetAtomicNum()),
	      ++labelcount[atom->GetAtomicNum()]);

      str = atom->GetType();

      ttab.Translate(str1,str);

      strcpy(rlabel,"<1>");
      strcpy(rnum,"1");

      sprintf(buffer,"%7d%1s%-6s%12.4f%10.4f%10.4f%1s%-5s%4s%1s %-8s%10.4f",
	    ct,"",label, 
	      atom->GetX(),atom->GetY(),atom->GetZ(),
	      "",str1.c_str(),
	      rnum,"",rlabel,
	      atom->GetPartialCharge());
      ofs << buffer << endl;
    }
cerr << endl;

  ofs << "@<TRIPOS>BOND" << endl;
  OBBond *bond;
  vector<OBEdgeBase*>::iterator j;
  for (bond = mol.BeginBond(j);bond;bond = mol.NextBond(j))
    {
      if (bond->GetBO() == 5) strcpy(label,"ar");
      else if (bond->IsAmide()) strcpy(label,"am");
      else sprintf(label,"%d",bond->GetBO());
      sprintf(buffer, "%6d%6d%6d%3s%2s", 
	      bond->GetIdx()+1,backmap[bond->GetBeginAtomIdx()],backmap[bond->GetEndAtomIdx()],
	      "",label);
      ofs << buffer << endl;
    }
  ofs << endl;

  return(true);
}

}

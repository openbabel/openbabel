/**********************************************************************
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (c) 2001-2003 by Geoffrey R. Hutchison

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include "mol.h"

using namespace std;
namespace OpenBabel {
  //FF
  //extern OBAtomTyper atomtyper; //CM

int wedgeorhatch(int flag) {
  // convert from wedge/hash bond flags to MDL stereochemistry
  if (flag & OB_WEDGE_BOND)
    return 1;
  else if (flag & OB_HASH_BOND )
    return 6;
  else
    return 0;
}

bool ReadSDFile(istream &ifs,OBMol &mol,const char *title) {
  int i,natoms,nbonds;
  char buffer[BUFF_SIZE];
  char *comment = NULL;
  string r1,r2;

  if (!ifs.getline(buffer,BUFF_SIZE)) return(false);
  mol.SetTitle(buffer);
  if (!ifs.getline(buffer,BUFF_SIZE)) return(false); //creator
  if (!ifs.getline(buffer,BUFF_SIZE)) return(false); //comment
  if (strlen(buffer) > 0) {
    comment = new char [strlen(buffer)+1];
    strcpy(comment,buffer);
  }

  if (!ifs.getline(buffer,BUFF_SIZE)) return(false); //atoms and bonds
  r1 = buffer;
  natoms = atoi((r1.substr(0,3)).c_str());
  nbonds = atoi((r1.substr(3,3)).c_str());

  mol.BeginModify();
  mol.ReserveAtoms(natoms);
  double x,y,z;
  char type[5];
  vector3 v;
  OBAtom atom;
  int charge;

  for (i = 0;i < natoms;i++) {
    if (!ifs.getline(buffer,BUFF_SIZE))
      return(false);

    if (sscanf(buffer,"%lf %lf %lf %s %*d %d",&x,&y,&z,type,&charge) != 5)
      return(false);
    v.SetX(x);v.SetY(y);v.SetZ(z);
    atom.SetVector(x, y, z);
    atom.SetAtomicNum(etab.GetAtomicNum(type));
    atom.SetType(type);

    switch (charge) {
    case 0: break;
    case 3: atom.SetFormalCharge(1); break;
    case 2: atom.SetFormalCharge(2); break;
    case 1: atom.SetFormalCharge(3); break;
    case 5: atom.SetFormalCharge(-1); break;
    case 6: atom.SetFormalCharge(-2); break;
    case 7: atom.SetFormalCharge(-3); break;
    }

    if (!mol.AddAtom(atom))
      return(false);
    atom.Clear();
  }

  int start,end,order,flag,stereo;
  for (i = 0;i < nbonds;i++) {
    flag = 0;
    if (!ifs.getline(buffer,BUFF_SIZE))
      return(false);
    r1 = buffer;
    start = atoi((r1.substr(0,3)).c_str());
    end = atoi((r1.substr(3,3)).c_str());
    order = atoi((r1.substr(6,3)).c_str());
    order = (order == 4) ? 5 : order;
    if (r1.size() >= 12) {  //handle wedge/hash data
      stereo = atoi((r1.substr(9,3)).c_str());
      if (stereo) {
        if (stereo == 1) flag |= OB_WEDGE_BOND;
        if (stereo == 6) flag |= OB_HASH_BOND;
      }
    }

    if (!mol.AddBond(start,end,order,flag)) return(false);
  }

	//CM start 18 Sept 2003
	//Read Properties block, currently only M RAD and M CHG 

  while(ifs.getline(buffer,BUFF_SIZE))
	{
    if(!strchr(buffer,'M')) continue;
    r1 = buffer;
		int n = atoi((r1.substr(6,3)).c_str()); //entries on this line
		if(n==0) break;
		int pos = 10;
		for(;n>0;n--,pos+=8)
		{
			int atomnumber = atoi((r1.substr(pos,3)).c_str());
			if (atomnumber==0) break;
			OBAtom* at;
			at=mol.GetAtom(atomnumber); //atom numbers start at 1
			int value = atoi((r1.substr(pos+4,3)).c_str());
			if(r1.substr(3,3)=="RAD")
				at->SetSpinMultiplicity(value);
			else if(r1.substr(3,3)=="CHG")
				at->SetFormalCharge(value);
			//Although not done here,according to the specification, 
			//previously set formal charges should be reset to zero
			// Lines setting several other properties are not implemented
		}
	}
  //FF spin multiplicity for H-deficient atoms no longer assigned
  //in AssignImplicitValence
  //atomtyper.AssignImplicitValence(mol); //and set _spinmultiplicities for H-deficient atoms
  mol.AssignSpinMultiplicity();
  //FF end
  //
	//CM end

  mol.EndModify();

  if (comment)
  {
	  OBCommentData *cd = new OBCommentData;
	  mol.SetData(cd);
  }

  while (ifs.getline(buffer,BUFF_SIZE)) {
    // RWT 4/7/2001
    // added to get properties
    if (strstr(buffer,"<")) {
      string buff(buffer);
      size_t lt=buff.find("<")+1;
      size_t rt = buff.find_last_of(">");
      string attr = buff.substr(lt,rt-lt);
      ifs.getline(buffer,BUFF_SIZE);

	  OBPairData *dp = new OBPairData;
	  dp->SetAttribute(attr);
	  dp->SetValue(buffer);
      mol.SetData(dp);
    }
    // end RWT    

    if (!strncmp(buffer,"$$$$",4)) break;
  }

  return(true);
}

bool WriteSDFile(ostream &ofs,OBMol &mol,const char *dimension) {
  char buff[BUFF_SIZE];  

  if (mol.NumAtoms() > 999) // Three digits!
    {
      ThrowError("MDL Molfile conversion failed: Molecule is too large to convert.");
      ThrowError("  File format limited to 999 atoms.");
      cerr << "  Molecule size: " << mol.NumAtoms() << " atoms " << endl;
      return(false);
    }

  if (mol.NumBonds() > 999) // Three digits!
    {
      ThrowError("MDL Molfile conversion failed: Molecule is too large to convert.");
      ThrowError("  File format limited to 999 bonds.");
      cerr << "  Molecule size: " << mol.NumBonds() << " atoms " << endl;
      return(false);
    }

  ofs << mol.GetTitle() <<  endl;
//sprintf(buff,"  -ISIS-            %s",dimension); replaced by CM
  sprintf(buff," OpenBabel          %s",dimension); // CM 18 Sept 2003
  ofs << buff << endl;

  if (mol.HasData(obCommentData))
    {
      OBCommentData *cd = (OBCommentData*)mol.GetData(obCommentData);
      ofs << cd->GetData() << endl;
    }
  else
      ofs << endl;

  sprintf(buff,"%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d",
          mol.NumAtoms(),mol.NumBonds(),0,0,0,0,0,0,0,0,1);// CM 18 Sept 2003 1 was 0 (# extra lines)
  ofs << buff << endl;

  OBAtom *atom;
  vector<OBNodeBase*>::iterator i;
  int charge;
  for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i)) {
    switch (atom->GetFormalCharge()) {
    case 1: charge = 3; break;
    case 2: charge = 2; break;
    case 3: charge = 1; break;
    case -1: charge = 5; break;
    case -2: charge = 6; break;
    case -3: charge = 7; break;
    default:
      charge=0; break;
    }

    sprintf(buff,"%10.4f%10.4f%10.4f %-3s%2d%3d%3d%3d%3d",
            atom->GetX(),
            atom->GetY(),
            atom->GetZ(),
            (etab.GetSymbol(atom->GetAtomicNum())),
            0,charge,0,0,0);    
    ofs << buff << endl;
  }

  //so the bonds come out sorted
  OBAtom *nbr;
  OBBond *bond;
  vector<OBEdgeBase*>::iterator j;
  for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
    for (nbr = atom->BeginNbrAtom(j);nbr;nbr = atom->NextNbrAtom(j))
      if (atom->GetIdx() < nbr->GetIdx()) {
        bond = (OBBond*) *j;
        sprintf(buff,"%3d%3d%3d%3d%3d%3d",
                bond->GetBeginAtomIdx(),
                bond->GetEndAtomIdx(),
                (bond->GetBO() == 5) ? 4 : bond->GetBO(),
                wedgeorhatch(bond->GetFlags()),0,0);
        ofs << buff << endl;
      }

  //CM start 18 Sept 2003
  //For radicals
  char txt[50];
  *buff=0;
  int val, radcount=0;
  for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
    {
      if(atom->GetSpinMultiplicity())
	{
	  sprintf(txt,"%3d %3d ",atom->GetIdx(),atom->GetSpinMultiplicity()); //radicals=>2 all carbenes=>3	
	  strcat(buff,txt);
	  radcount++;
	}
    }
  if (radcount)
    {
      sprintf(txt,"M  RAD%3d ",radcount);
      ofs << txt << buff << endl;
    }
  // CM end

  ofs << "M  END" << endl;

  // RWT 4/7/2001
  // now output properties if they exist
  // MTS 4/17/2001
  // changed to use new OBGenericData class
  vector<OBGenericData*>::iterator k;
  vector<OBGenericData*> vdata = mol.GetData();
  for (k = vdata.begin();k != vdata.end();k++)
	  if ((*k)->GetDataType() == obPairData)
	  {
		  ofs << ">  <" << (*k)->GetAttribute() << ">" << endl;
		  ofs << ((OBPairData*)(*k))->GetValue() << endl << endl;
	  }

  // end RWT


  ofs << "$$$$" << endl;

  return(true);
}

}

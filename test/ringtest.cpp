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

bool TestRings(void)
{
  ifstream mifs;
  if (!SafeOpen(mifs,"attype.00.smi")) return(false);

  ifstream rifs;
  if (!SafeOpen(rifs,"ringresults.txt")) return(false);

  unsigned int size;
  OEBond *bond;
  OEAtom *atom;
  int count;
  char buffer[BUFF_SIZE];
  vector<string> vs;
  vector<OERing*> vr;
  vector<bool> vb;
  vector<int> vi;
  OEMol mol(SMI,SMI);
  vector<string>::iterator i;
  vector<OEBond*>::iterator j;
  vector<OEAtom*>::iterator k;
  vector<OERing*>::iterator m;

  for (;mifs;)
    {
      mol.Clear();
      mifs >> mol;
      if (mol.Empty()) continue;
      if (!rifs.getline(buffer,BUFF_SIZE))
	{
	  ThrowError("error reading reference data");
	  return(false);
	}

      vb.clear();
      vb.resize(mol.NumBonds(),false);
      //check ring bonds
      tokenize(vs,buffer);
      for (i = vs.begin();i != vs.end();i++)
	  vb[atoi((char*)i->c_str())] = true;
      
      for (bond = mol.BeginBond(j);bond;bond = mol.NextBond(j))
	if (vb[bond->GetIdx()] != bond->IsInRing())
	  {
	    ThrowError("ring bond data different than reference");
	    ThrowError((char*)mol.GetTitle());
	    return(false);
	  }

      vr = mol.GetSSSR();
      if (!rifs.getline(buffer,BUFF_SIZE))
	{
	  ThrowError("error reading reference data");
	  return(false);
	}
      sscanf(buffer,"%d",&size);
      if (vr.size() != size) //check SSSR size
	{
	  ThrowError("SSSR size different than reference");
	  ThrowError((char*)mol.GetTitle());
	  return(false);
	}

      if (!rifs.getline(buffer,BUFF_SIZE))
	{
	  ThrowError("error reading reference data");
	  return(false);
	}

      tokenize(vs,buffer);
      i = vs.begin();
      for (atom = mol.BeginAtom(k);atom;atom = mol.NextAtom(k))
	{
	  if (i == vs.end())
	    {
	      ThrowError("Error in SSSR count");
	      ThrowError((char*)mol.GetTitle());
	      return(false);
	    }
	  count = 0;
	  for (m = vr.begin();m != vr.end();m++)
	    if ((*m)->_pathset[atom->GetIdx()])
	      count++;

	  if (atoi((char*)i->c_str()) != count)
	    {
	      ThrowError("Ring membership test failed");
	      ThrowError((char*)mol.GetTitle());
	      return(false);
	    }

	  i++;
	}


    }


  return(true);
}

void GenerateRingReference()
{
  ifstream ifs;
  if (!SafeOpen(ifs,"attype.00.smi")) return;

  ofstream ofs;
  if (!SafeOpen(ofs,"ringresults.txt")) return;

  int count;
  OEAtom *atom;
  OEBond *bond;
  char buffer[BUFF_SIZE];
  vector<OERing*> vr;
  vector<OEBond*>::iterator i;
  vector<OEAtom*>::iterator j;
  vector<OERing*>::iterator k;
  OEMol mol(SMI,SMI);

  for (;ifs;)
    {
      mol.Clear();
      ifs >> mol;
      if (mol.Empty()) continue;
      
      //write out ring bonds
      for (bond = mol.BeginBond(i);bond;bond = mol.NextBond(i))
	if (bond->IsInRing())
	  {
	    sprintf(buffer,"%3d",bond->GetIdx());
	    ofs << buffer;
	  }
      ofs << endl;

      vr = mol.GetSSSR();
      //write the total number of rings
      ofs << vr.size() << endl;

      //write the number of rings that each atom is a member of
      for (atom = mol.BeginAtom(j);atom;atom = mol.NextAtom(j))
	{
	  count = 0;
	  for (k = vr.begin();k != vr.end();k++)
	    if ((*k)->_pathset[atom->GetIdx()])
	      count++;

	  sprintf(buffer,"%3d",count);
	  ofs << buffer;
	}
      ofs << endl;

    }

}


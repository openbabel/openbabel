/**********************************************************************
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (c) 2001-2003 Geoffrey R. Hutchison

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include "mol.h"

namespace OpenBabel {
  bool SafeOpen(std::ifstream &fs, char *filename);
  bool SafeOpen(std::ofstream &fs, char *filename);
}

using namespace std;
using namespace OpenBabel;

void GenerateRingReference();

int main(int argc,char *argv[])
{
  if (argc != 1) {
    if (strncmp(argv[1], "-g", 2))
      {
	cout << "Usage: ringtest" << endl;
	cout << "   Tests Open Babel ring perception testing." << endl;
	return 0;
      } else {
	GenerateRingReference();
	return 0;
      }
  }

  cout << endl << "Testing RINGS..." << endl;
  
#ifdef TESTDATADIR
  string testdatadir = TESTDATADIR;
  string results_file = testdatadir + "ringresults.txt";
  string smilestypes_file = testdatadir + "attype.00.smi";
#else
  string results_file = "ringresults.txt";
  string smilestypes_file = "attype.00.smi";
#endif

  std::ifstream mifs;
  if (!SafeOpen(mifs, (char*)smilestypes_file.c_str())) {
    return -1; // test failed
  }

  std::ifstream rifs;
  if (!SafeOpen(rifs, (char*)results_file.c_str())) {
    return -1; // test failed
  }

  unsigned int size;
  OBBond *bond;
  OBAtom *atom;
  int count;
  char buffer[BUFF_SIZE];
  vector<string> vs;
  vector<OBRing*> vr;
  vector<bool> vb;
  vector<int> vi;
  OBMol mol(SMI,SMI);
  vector<string>::iterator i;
  vector<OBEdgeBase*>::iterator j;
  vector<OBNodeBase*>::iterator k;
  vector<OBRing*>::iterator m;
  OBFileFormat ff;

  for (;mifs;)
    {
      mol.Clear();
      ff.ReadMolecule(mifs, mol);
      if (mol.Empty()) continue;
      if (!rifs.getline(buffer,BUFF_SIZE))
	{
	  ThrowError("error reading reference data");
	  return -1; // test failed
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
	    return -1; // test failed
	  }

      vr = mol.GetSSSR();
      if (!rifs.getline(buffer,BUFF_SIZE))
	{
	  ThrowError("error reading reference data");
	  return -1; // test failed
	}
      sscanf(buffer,"%d",&size);
      if (vr.size() != size) //check SSSR size
	{
	  ThrowError("SSSR size different than reference");
	  ThrowError((char*)mol.GetTitle());
	  return -1; // test failed
	}

      if (!rifs.getline(buffer,BUFF_SIZE))
	{
	  ThrowError("error reading reference data");
	  return -1; // test failed
	}

      tokenize(vs,buffer);
      i = vs.begin();
      for (atom = mol.BeginAtom(k);atom;atom = mol.NextAtom(k))
	{
	  if (i == vs.end())
	    {
	      ThrowError("Error in SSSR count");
	      ThrowError((char*)mol.GetTitle());
	      return -1; // test failed
	    }
	  count = 0;
	  for (m = vr.begin();m != vr.end();m++)
	    if ((*m)->_pathset[atom->GetIdx()])
	      count++;

	  if (atoi((char*)i->c_str()) != count)
	    {
	      ThrowError("Ring membership test failed");
	      ThrowError((char*)mol.GetTitle());
	      return -1; // test failed
	    }

	  i++;
	}


    }

  // Passed tests
  return 0;
}

void GenerateRingReference()
{
  std::ifstream ifs;
  if (!SafeOpen(ifs,"attype.00.smi")) return;

  std::ofstream ofs;
  if (!SafeOpen(ofs,"ringresults.txt")) return;

  int count;
  OBAtom *atom;
  OBBond *bond;
  char buffer[BUFF_SIZE];
  vector<OBRing*> vr;
  vector<OBEdgeBase*>::iterator i;
  vector<OBNodeBase*>::iterator j;
  vector<OBRing*>::iterator k;
  OBMol mol(SMI,SMI);
  OBFileFormat ff;

  for (;ifs;)
    {
      mol.Clear();
      ff.ReadMolecule(ifs, mol);
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

  ThrowError("Ring perception test results written successfully");
}


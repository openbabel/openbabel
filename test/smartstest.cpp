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

namespace OpenBabel {
  bool SafeOpen(std::ifstream &fs, char *filename);
  bool SafeOpen(std::ofstream &fs, char *filename);
}

using namespace std;
using namespace OpenBabel;

bool TestSmarts()
{
  std::ifstream ifs;
  if (!SafeOpen(ifs,"smartstest.txt")) return(false);

  //read in the SMARTS test patterns
  char buffer[BUFF_SIZE];
  vector<OBSmartsPattern*> vsp;
  for (;ifs.getline(buffer,BUFF_SIZE);) 
    {
      OBSmartsPattern *sp = new OBSmartsPattern;

      if (sp->Init(buffer)) vsp.push_back(sp);
      else                  delete sp;
    }
  ifs.close();
  
  std::ifstream rifs;
  if (!SafeOpen(rifs,"smartsresults.txt")) return(false);
  unsigned int npats;
  rifs.getline(buffer,BUFF_SIZE);
  sscanf(buffer,"%d %*s",&npats);

  //make sure the number of SMARTS patterns is the same as in the 
  //reference data
  if (npats != vsp.size())
    {
      ThrowError("Correct number of patterns not read in");
      sprintf(buffer,"Read in %d, expected %d", vsp.size(), npats);
      ThrowError(buffer);
      return(false);
    }

  std::ifstream mifs;
  if (!SafeOpen(mifs,"attype.00.smi")) return(false);

  unsigned int k;
  OBMol mol(SMI,SMI);
  vector<string> vs;
  vector<OBSmartsPattern*>::iterator i;
  vector<vector<int> > mlist;
  OBFileFormat ff;

  //read in molecules, match SMARTS, and compare results to reference data
  for (;mifs;)
    {
      mol.Clear();
      ff.ReadMolecule(mifs, mol);
      if (mol.Empty()) continue;
      
      for (i = vsp.begin();i != vsp.end();i++)
	{
	  if (!rifs.getline(buffer,BUFF_SIZE))
	    {
	      ThrowError("error reading reference data");
	      return(false);
	    }

	  tokenize(vs,buffer);
	  (*i)->Match(mol);
	  mlist = (*i)->GetMapList();
	  if (mlist.size() != vs.size())
	    {
	      ThrowError("number of matches different than reference");
	      cerr << "expected " << vs.size() << " got " << mlist.size() << endl;
	      ThrowError((char*)mol.GetTitle());
	      //	      ThrowError((*i)->GetSMARTS());
	      return(false);
	    }
	  
	  if (mlist.size())
	    {
	      for (k = 0;k < vs.size();k++)
		if (atoi((char*)vs[k].c_str()) != mlist[k][0])
		{
		  ThrowError("matching atom numbers different than reference");
		  ThrowError((char*)mol.GetTitle());
		  //		  ThrowError((*i)->GetSMARTS());
		  return(false);
		}
	    }
	}
    }

  return(true);
}

void GenerateSmartsReference()
{
  std::ifstream ifs;
  if (!SafeOpen(ifs,"smartstest.txt")) return;

  char buffer[BUFF_SIZE];
  vector<OBSmartsPattern*> vsp;
  for (;ifs.getline(buffer,BUFF_SIZE);)
    {
      OBSmartsPattern *sp = new OBSmartsPattern;

      if (sp->Init(buffer)) vsp.push_back(sp);
      else                  delete sp;
    }

  std::ofstream ofs;
  if (!SafeOpen(ofs,"smartsresults.txt")) return;
    
  ofs << vsp.size() << " patterns" << endl;
  std::ifstream mifs;
  if (!SafeOpen(mifs,"attype.00.smi")) return;

  vector<int> vm;
  vector<vector<int> > mlist;
  vector<vector<int> >::iterator j;
  vector<OBSmartsPattern*>::iterator i;
  OBMol mol(SMI,SMI);
  OBFileFormat ff;

  for (;mifs;)
    {
      mol.Clear();
      ff.ReadMolecule(mifs, mol);

      if (mol.Empty()) continue;
      for (i = vsp.begin();i != vsp.end();i++)
	{
	  (*i)->Match(mol);
	  mlist = (*i)->GetMapList();
	  for (j = mlist.begin();j != mlist.end();j++)
	    {
	     sprintf(buffer,"%3d",*(j->begin()));
	     ofs << buffer;
	    }
	  ofs << endl; 
	}
    }
  

  ThrowError("SMARTS test results written successfully");
}

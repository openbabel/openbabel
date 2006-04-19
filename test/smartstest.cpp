 /**********************************************************************
 smartstest.cpp - Test SMARTS algorithms and atom typing.

 This file is part of the Open Babel project.
 For more information, see <http://openbabel.sourceforge.net/>

 Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
 Some portions Copyright (C) 2001-2005 Geoffrey R. Hutchison

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation version 2 of the License.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 ***********************************************************************/

#include <fstream>

#include "babelconfig.h"
#include "mol.h"
#include "obconversion.h"

namespace OpenBabel
{
  bool SafeOpen(std::ifstream &fs, char *filename);
  bool SafeOpen(std::ofstream &fs, char *filename);
}

using namespace std;
using namespace OpenBabel;

void GenerateSmartsReference();

int main(int argc,char *argv[])
{
  if (argc != 1)
    {
      if (strncmp(argv[1], "-g", 2))
	{
	  cout << "Usage: smartstest" << endl;
	  cout << "   Tests Open Babel SMILES/SMARTS pattern matching." << endl;
	  return 0;
	}
      else
	{
	  GenerateSmartsReference();
	  return 0;
	}
    }
  
  cout << endl << "# Testing SMARTS...  " << endl;
  
#ifdef TESTDATADIR
  string testdatadir = TESTDATADIR;
  string smarts_file = testdatadir + "smartstest.txt";
  string results_file = testdatadir + "smartsresults.txt";
  string smilestypes_file = testdatadir + "attype.00.smi";
#else

   string smarts_file = "smartstest.txt";
   string results_file = "smartsresults.txt";
   string smilestypes_file = "attype.00.smi";
#endif

   std::ifstream ifs;
   if (!SafeOpen(ifs, (char*)smarts_file.c_str()))
     {
       cout << "Bail out! Cannot read " << smarts_file << endl;
       return -1; // test failed
     }

   //read in the SMARTS test patterns
   char buffer[BUFF_SIZE];
   vector<OBSmartsPattern*> vsp;
   for (;ifs.getline(buffer,BUFF_SIZE);)
     {
       OBSmartsPattern *sp = new OBSmartsPattern;

       if (sp->Init(buffer))
         vsp.push_back(sp);
       else
         delete sp;
     }
   ifs.close();

   std::ifstream rifs;
   if (!SafeOpen(rifs, (char*)results_file.c_str()))
     {
       cout << "Bail out! Cannot read in results file " << results_file << endl;
       return -1; // test failed
     }
   unsigned int npats;
   rifs.getline(buffer,BUFF_SIZE);
   sscanf(buffer,"%d %*s",&npats);

   //make sure the number of SMARTS patterns is the same as in the
   //reference data
   if (npats != vsp.size())
     {
       cout << "Bail out! Correct number of patterns not read in" <<
         "Read in " << vsp.size() << " expected " << npats << endl;
       return -1; // test failed
     }

   std::ifstream mifs;
   if (!SafeOpen(mifs, (char*)smilestypes_file.c_str()))
     {
       cout << "Bail out! Cannot read atom types " << smilestypes_file << endl;
       return -1; // test failed
     }

   unsigned int k;
   unsigned int res_line = 0;
   OBMol mol;
   vector<string> vs;
   vector<OBSmartsPattern*>::iterator i;
   vector<vector<int> > mlist;
   unsigned int currentTest = 1;

   OBConversion conv(&mifs, &cout);
   if (! conv.SetInAndOutFormats("SMI","SMI"))
     {
       cout << "Bail out! SMILES format is not loaded" << endl;
       return -1;
     }

   //read in molecules, match SMARTS, and compare results to reference data
   for (;mifs;)
     {
       mol.Clear();
       conv.Read(&mol);
       if (mol.Empty())
         continue;

       for (i = vsp.begin();i != vsp.end();i++)
         {
           if (!rifs.getline(buffer,BUFF_SIZE))
             {
               cout << "Bail out! Error reading reference data" << endl;
               return -1; // test failed
             }
           res_line++;

           tokenize(vs,buffer);
           (*i)->Match(mol);
           mlist = (*i)->GetMapList();
           if (mlist.size() != vs.size())
             {
               cout << "not ok " << currentTest++ 
                    << " # number of matches different than reference" << endl;
               cout << "# Expected " << vs.size() << " matches, found "
                    << mlist.size() << endl;
               cout << "# Error with molecule " << mol.GetTitle();
               cout << "#  on pattern " << (*i)->GetSMARTS() << endl;
               if (mlist.size())
                 cout << "# First match: atom #" << mlist[0][0] << endl;
             }
           else
             cout << "ok " << currentTest++ 
                  << " # correct number of matches" << endl;

           if (mlist.size())
             {
               for (k = 0;k < vs.size();k++)
		 {
		   if (atoi((char*)vs[k].c_str()) != mlist[k][0])
		     {
		       cout << "not ok " << currentTest++ 
			    << "# matching atom numbers different than reference" << endl;
		       cout << "# Expected " << vs[k] << " but found "
			    << mlist[k][0] << endl;
		       cout << "# Molecule: " << mol.GetTitle() << endl;
		       cout << "# Pattern: " << (*i)->GetSMARTS() << endl;
		     }
		   else
		     cout << "ok " << currentTest++ 
			  << " # correcting matching atom numbers" << endl;
		 }
	     }
	 }
     }

   // output the number of tests run
   cout << "1.." << currentTest-1 << endl;

   // Passed Test
   return 0;
 }

void GenerateSmartsReference()
{
    std::ifstream ifs;
    if (!SafeOpen(ifs,"smartstest.txt"))
        return;

    char buffer[BUFF_SIZE];
    vector<OBSmartsPattern*> vsp;
    for (;ifs.getline(buffer,BUFF_SIZE);)
    {
        OBSmartsPattern *sp = new OBSmartsPattern;

        if (sp->Init(buffer))
            vsp.push_back(sp);
        else
            delete sp;
    }

    std::ofstream ofs;
    if (!SafeOpen(ofs,"smartsresults.txt"))
        return;

    ofs << vsp.size() << " patterns" << endl;
    std::ifstream mifs;
    if (!SafeOpen(mifs,"attype.00.smi"))
        return;

    vector<int> vm;
    vector<vector<int> > mlist;
    vector<vector<int> >::iterator j;
    vector<OBSmartsPattern*>::iterator i;
    OBMol mol;
    OBConversion conv(&mifs, &cout);

    if(! conv.SetInAndOutFormats("SMI","SMI"))
    {
        ThrowError("SMILES format is not loaded");
        return;
    }

    for (;mifs;)
    {
        mol.Clear();
        conv.Read(&mol);

        if (mol.Empty())
            continue;
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


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
#include "oeutil.h"
#include "parsmart.h"
#include "patty.h"

// Simple programmable atom typer
// WPW - 070199
// Usage is in sample main below

namespace OpenEye {

void patty::read_rules(const string &infile)
{
  ifstream ifs, ifs1, *ifsP;
  vector<string> vs;
  char buffer[BUFF_SIZE];
  char tmp_str[BUFF_SIZE];
  char patty_dir[BUFF_SIZE];

  _sp.resize(1000);
  
  ifs.open(infile.c_str());
  ifsP= &ifs;
  if (!ifs)
    {
      if (getenv("OE_DIR") == NULL)
      {
        cerr << "The OE_DIR environment variable is not defined" << endl;
        cerr << "Please define it so the program can find " << infile << endl;
        exit(0);
      }
      else
        strcpy(patty_dir,getenv("OE_DIR"));
      strcat(patty_dir,FILE_SEP_CHAR);  
      strcat(patty_dir,infile.c_str());  
      ifs1.open(patty_dir);  
      ifsP= &ifs1;
 //     if (!ifs1)
  //    {
   //     cerr << "Could not open " << patty_dir << endl;
    //    exit(0);
     // }
    }

  int i = 0;
  if (!ifsP){
        cerr << "Could not open " << patty_dir << endl;
        exit(0);
  }
  while (ifsP->getline(buffer,BUFF_SIZE))
    {
      if (buffer[0] != '#')
    {
      tokenize(vs,buffer," \t\n");
      if (vs.size() >= 2)
        {
          strcpy(tmp_str,vs[0].c_str());
          _sp[i]->Init(tmp_str);
          smarts.push_back(vs[0]);
          typ.push_back(vs[1]);
          i++;
        }
    }
    }
  _sp.resize(i);
}

void patty::assign_rules(vector<string> &rules)
{
	vector<string> vs;
	char buffer[BUFF_SIZE];
	char tmp_str[BUFF_SIZE];
	unsigned int i;

	_sp.resize(1000);
	for ( i = 0 ; i < rules.size() ; i++ )
	{
		strncpy(buffer, rules[i].c_str(), BUFF_SIZE);
		if (buffer[0] != '#')
		{
			tokenize(vs,buffer," \t\n");
			if (vs.size() >= 2)
			{
				strcpy(tmp_str,vs[0].c_str());
				_sp[i]->Init(tmp_str);
				smarts.push_back(vs[0]);
				typ.push_back(vs[1]);
			}
			else
				i--;
		}
		else
			i--;
	}
	_sp.resize(i);  
}

void patty::assign_types(OEMol &mol,vector<string> &atm_typ)
{
  atm_typ.resize(mol.NumAtoms()+1);

  for (unsigned int i = 0; i < _sp.size(); i++)
    {
      _sp[i]->Match(mol);
      vector<vector<int> > match = _sp[i]->GetMapList();
      //vector<vector<int> >& match = _sp[i]->GetMapList();
      if (match.size())
	{
	  if (debug)
	    cout << typ[i] << " " << smarts[i] << " matched " ;
      
	  for (unsigned int j = 0; j < match.size(); j++)
	    {
	      if (debug)
		cout << match[j][0] << " ";
	      atm_typ[match[j][0]] = typ[i];
	    }
	  if (debug)
	    cout << endl;
	}
    }
}

void patty::assign_types(OEMol &mol,vector<int> &atm_typ)
{
  atm_typ.resize(mol.NumAtoms()+1);

  for (unsigned int i = 0; i < _sp.size(); i++)
    {
      _sp[i]->Match(mol);
      vector<vector<int> > match = _sp[i]->GetMapList();
      //vector<vector<int> >& match = _sp[i]->GetMapList();
      if (match.size())
	{
	  if (debug)
	    cout << typ[i] << " " << smarts[i] << " matched " ;
      
	  for (unsigned int j = 0; j < match.size(); j++)
	    {
	      if (debug)
		cout << match[j][0] << " ";
	      atm_typ[match[j][0]] = type_to_int(typ[i]);
	    }
	  if (debug)
	    cout << endl;
	}
    }
}


int patty::type_to_int(const string &type, bool failOnUndefined)
{

  int result;


    switch(toupper(type.c_str()[0]))
      {
      case 'C' : // CAT - CATION
         result = PT_CATION;
         break;
      case 'A' :
         if (toupper(type.c_str()[1]) == 'N') // ANI - ANION
            result = PT_ANION;
         else
            result = PT_ACCEPTOR;
         break;
      case 'P' : // POL - POLAR
         result = PT_POLAR;
         break;
      case 'D' : // DON - DONOR
         result = PT_DONOR;
         break;
      case 'H' : // HYD - HYDROPHOBIC
         result = PT_HYDROPHOBIC;
         break;
			case 'M' : // Metal
				 result = PT_METAL;
				 break;
      case 'O' : // OTH - OTHER
         result = PT_OTHER;
         break;
      default :
          // This was added by Brian,
          // Behavior will fail if type is undefined
          if (failOnUndefined){
              cerr << "Unable to find type of feature passed in " << endl;
              cerr << "Feature passed in is " << type << endl;
              exit(-1);
          }else{
              result = 7;
          }
      }
      return(result);
}

}

#ifdef I_EVER_BUY_A_MICHAEL_BOLTON_ALBUM

int main(int argc, char *argv[])
{
  OEMol mol(SDF,SDF);
  vector<string> types;

  ifstream ifs(argv[1]);
  if (!ifs)
    {
      cerr << "Could not open argv[1] " << endl;
      exit(0);
    }

  patty p("simple.txt");
  for (;;)
    {
      ifs >> mol;
      if (!mol.NumAtoms()) break;
      p.assign_types(mol,types);
      mol.Clear();
    }

  for (int i = 1; i < types.size(); i++)
    {
      cout << i << " " << types[i] << endl;
    }
}

#endif













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

#ifndef OE_PATTY
#define OE_PATTY

namespace OpenEye {
#define PT_CATION      1
#define PT_ANION       2
#define PT_ACCEPTOR    3
#define PT_POLAR       4
#define PT_DONOR       5
#define PT_HYDROPHOBIC 6
#define PT_OTHER       7
#define PT_METAL	   8

class patty
{
  vector<OESmartsPattern*> _sp;
  vector<string> smarts;
  vector<string> typ;
  bool debug;

  public :

  patty() {debug = false;}
  patty(char *s) 
    {
      debug = false;
      read_rules(string(s));
    }
  
  patty(const string &s) 
    {
      debug = false;
      read_rules(s);
    }
  ~patty()
    {
      vector<OESmartsPattern*>::iterator i;
      for (i = _sp.begin();i != _sp.end();i++) delete *i;
    }
  void debug_on() {debug = true;}
  void debug_off() {debug = false;}
  void read_rules(const string &infile);
  void assign_rules(vector<string> &rules);
  void assign_types(OEMol &mol,vector<string> &atm_typ);
  void assign_types(OEMol &mol,vector<int> &atm_typ);
  int type_to_int(const string &type, bool failOnUndefined= false);
};

}

#endif

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

#ifndef MOLVECTOR_H
#define MOLVECTOR_H

#include "mol.h"

namespace OpenBabel {

class OEMolVector
{
  vector<OEMol*> _molvec;

  public :
  
  ~OEMolVector();
  OEMolVector() {};
  OEMolVector(ifstream &ifs, int nToRead= -1)
    {
      Read(ifs,SDF,SDF,nToRead);
    }
  OEMolVector(ifstream &ifs, const io_type in_type, const io_type out_type, int nToRead= -1)
  {
    Read(ifs,in_type,out_type,nToRead);
  }
  void Read(ifstream &ifs, const io_type in_type, const io_type out_type, int nToRead);
  void Write(ofstream &ofs);
  int GetSize() { return(_molvec.size());}
  vector<OEMol *> &GetVector() {return(_molvec);}
  OEMol *GetMol(int i);
  void   SetMol(int i, OEMol *mol) { _molvec[i]= mol; }

  bool ReadConfs(ifstream &ifs, const io_type in_type, const io_type out_type);
  bool ReadConfs(ifstream &ifs)
  {
    bool retval = ReadConfs(ifs,SDF,SDF);
    return(retval);
  }
};

}

#endif


/**********************************************************************
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (c) 2001-2002 by Geoffrey R. Hutchison

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

#ifndef OB_MOLVECTOR_H
#define OB_MOLVECTOR_H

#include "mol.h"

namespace OpenBabel {

class OBMolVector
{
  std::vector<OBMol*> _molvec;

  public :
  
  ~OBMolVector();
  OBMolVector() {};
  OBMolVector(std::ifstream &ifs, int nToRead= -1)
    {
      Read(ifs,SDF,SDF,nToRead);
    }
  OBMolVector(std::ifstream &ifs, const io_type in_type, const io_type out_type, int nToRead= -1)
  {
    Read(ifs,in_type,out_type,nToRead);
  }
  void Read(std::ifstream &ifs, const io_type in_type, const io_type out_type, int nToRead);
  void Write(std::ofstream &ofs);
  int GetSize() { return(_molvec.size());}
  std::vector<OBMol *> &GetVector() {return(_molvec);}
  OBMol *GetMol(int i);
  void   SetMol(int i, OBMol *mol) { _molvec[i]= mol; }

  bool ReadConfs(std::ifstream &ifs, const io_type in_type, const io_type out_type);
  bool ReadConfs(std::ifstream &ifs)
  {
    bool retval = ReadConfs(ifs,SDF,SDF);
    return(retval);
  }
};

} // end namespace OpenBabel

#endif // OB_MOLVECTOR_H


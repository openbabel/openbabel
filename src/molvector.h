/**********************************************************************
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (c) 2001-2003 by Geoffrey R. Hutchison

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

  // Class introduction in molvector.cpp

class OBMolVector
{
  std::vector<OBMol*> _molvec;

  public :
  
  ~OBMolVector();
  OBMolVector() {};
  OBMolVector(std::ifstream &ifs, 
	      const io_type in_type=SDF, const io_type out_type=SDF,
	      int nToRead=-1)
  { Read(ifs,in_type,out_type,nToRead); }

  void Read(std::ifstream &ifs, 
	    const io_type in_type=SDF, const io_type out_type=SDF,
	    int nToRead=-1);
  bool ReadConfs(std::ifstream &ifs,
		 const io_type in_type=SDF, const io_type out_type=SDF);

  void Write(std::ostream &os, const char *options="");

  unsigned int GetSize() { return(_molvec.size());}

  std::vector<OBMol *> &GetVector() {return(_molvec);}
  OBMol *GetMol(unsigned int i);
  void   SetMol(unsigned int i, OBMol *mol) { _molvec[i]= mol; }
  void   PushMol(OBMol *mol);
};

} // end namespace OpenBabel

#endif // OB_MOLVECTOR_H


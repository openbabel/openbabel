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

#ifndef OB_SMI_H
#define OB_SMI_H

#include "mol.h"
#include "typer.h"
#include "bitvec.h"

namespace OpenBabel{

class OBSmiNode
{
  OBAtom *_atom,*_parent;
  std::vector<OBSmiNode*> _nextnode;
  std::vector<OBBond*> _nextbond;
public:
  OBSmiNode(OBAtom *atom);
  ~OBSmiNode();
  int        Size()  {return((_nextnode.empty())?0:_nextnode.size());}
  void       SetParent(OBAtom *a)              {_parent = a;}
  void       SetNextNode(OBSmiNode*,OBBond*);
  OBAtom    *GetAtom()                         {return(_atom);}
  OBAtom    *GetParent()                       {return(_parent);}
  OBAtom    *GetNextAtom(int i)                {return(_nextnode[i]->GetAtom());}
  OBBond    *GetNextBond(int i)                {return(_nextbond[i]);}
  OBSmiNode *GetNextNode(int i)                {return(_nextnode[i]);}
};

class OBMol2Smi
{
  std::vector<int> _atmorder;
  std::vector<int> _storder;
  std::vector<bool> _aromNH;
  OBBitVec _uatoms,_ubonds;
  std::vector<OBEdgeBase*> _vclose;
  std::vector<std::pair<OBAtom*,std::pair<int,int> > > _vopen;
public:
  OBMol2Smi() {_vclose.clear();}
  ~OBMol2Smi() {}
  int          GetUnusedIndex();
  void         Init();
  void         CreateSmiString(OBMol&,char*);
  void         GetClosureAtoms(OBAtom*,std::vector<OBNodeBase*>&);
  void         FindClosureBonds(OBMol&);
  void         ToSmilesString(OBSmiNode *node,char *buffer);
  void         RemoveUsedClosures();
  void         AssignCisTrans(OBSmiNode*);
  bool         BuildTree(OBSmiNode*);
  bool         GetSmilesElement(OBSmiNode*,char*);
  bool         GetChiralStereo(OBSmiNode*,char*);
  void         CorrectAromaticAmineCharge(OBMol&);
  std::vector<std::pair<int,OBBond*> >  GetClosureDigits(OBAtom*);
  std::vector<int> &GetOutputOrder()                   {return(_atmorder);}
};

bool WriteTheSmiles(OBMol & mol,char *out);
bool WriteSmiles(std::ostream &ofs, OBMol &mol,char *title);

} //end namespace OpenBabel

#endif // OB_SMI_H

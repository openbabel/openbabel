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

#ifndef __SMI_H__
#define __SMI_H__

#include "mol.h"
#include "typer.h"
#include "bitvec.h"

namespace OpenEye{

class OESmiNode
{
  OEAtom *_atom,*_parent;
  vector<OESmiNode*> _nextnode;
  vector<OEBond*> _nextbond;
public:
  OESmiNode(OEAtom *atom);
  ~OESmiNode();
  int        Size()  {return((_nextnode.empty())?0:_nextnode.size());}
  void       SetParent(OEAtom *a)              {_parent = a;}
  void       SetNextNode(OESmiNode*,OEBond*);
  OEAtom    *GetAtom()                         {return(_atom);}
  OEAtom    *GetParent()                       {return(_parent);}
  OEAtom    *GetNextAtom(int i)                {return(_nextnode[i]->GetAtom());}
  OEBond    *GetNextBond(int i)                {return(_nextbond[i]);}
  OESmiNode *GetNextNode(int i)                {return(_nextnode[i]);}
};

class OEMol2Smi
{
  vector<int> _atmorder;
  vector<int> _storder;
  vector<bool> _aromNH;
  OEBitVec _uatoms,_ubonds;
  vector<OEBond*> _vclose;
  vector<pair<OEAtom*,pair<int,int> > > _vopen;
public:
  OEMol2Smi() {_vclose.clear();}
  ~OEMol2Smi() {}
  int          GetUnusedIndex();
  void         Init();
  void         CreateSmiString(OEMol&,char*);
  void         GetClosureAtoms(OEAtom*,vector<OEAtom*>&);
  void         FindClosureBonds(OEMol&);
  void         ToSmilesString(OESmiNode *node,char *buffer);
  void         RemoveUsedClosures();
  void         AssignCisTrans(OESmiNode*);
  bool         BuildTree(OESmiNode*);
  bool         GetSmilesElement(OESmiNode*,char*);
  bool         GetChiralStereo(OESmiNode*,char*);
  void         CorrectAromaticAmineCharge(OEMol&);
  vector<pair<int,OEBond*> >  GetClosureDigits(OEAtom*);
  vector<int> &GetOutputOrder()                   {return(_atmorder);}
};

bool WriteTheSmiles(OEMol & mol,char *out);
bool WriteSmiles(ostream &ofs, OEMol &mol,char *title);

}

#endif

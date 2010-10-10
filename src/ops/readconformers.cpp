/**********************************************************************
read-conformers.cpp - A OBOp for combining conformers during conversion.

Copyright (C) 2010 by Chris Morley

This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/
#include <openbabel/babelconfig.h>
#include <openbabel/op.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include "deferred.h"
#include <algorithm>

namespace OpenBabel
{

class OpReadConformers : public OBOp
{
public:
  OpReadConformers(const char* ID) : OBOp(ID, false){};
  const char* Description(){ return
    "Adjacent conformers combined into a single molecule\n"
    "If a molecule has the same structure as the preceding molecule, as determined\n"
    "from its SMILES, it is not output but its coordinates are added to the\n"
    "preceding molecule as an additional conformer. There can be multiple groups\n"
    "of conformers, but the molecules in each group must be adjacent.\n"
    ; }

  virtual bool WorksWith(OBBase* pOb)const{ return dynamic_cast<OBMol*>(pOb)!=NULL; }
  virtual bool Do(OBBase* pOb, const char* OptionText=NULL, OpMap* pOptions=NULL, OBConversion* pConv=NULL);
  virtual bool ProcessVec(std::vector<OBBase*>& vec);
};

/////////////////////////////////////////////////////////////////
OpReadConformers theOpReadConformers("readconformer"); //Global instance

/////////////////////////////////////////////////////////////////
bool OpReadConformers::Do(OBBase* pOb, const char* OptionText, OpMap* pOptions, OBConversion* pConv)
{
  //Make a deferred format and divert the output to it
  if(pConv && pConv->IsFirstInput())
    new DeferredFormat(pConv, this); //it will delete itself

  return true;
}

bool OpReadConformers::ProcessVec(std::vector<OBBase*>& vec)
{
  // DeferredFormat collects all the molecules, they are processed here, and Deferred Format outputs them
  OBConversion smconv;
  smconv.AddOption("n");
  if(!smconv.SetOutFormat("smi"))
  {
    obErrorLog.ThrowError(__FUNCTION__, "SmilesFormat is not loaded" , obError, onceOnly);
    return false;
  }

  std::string smiles, stored_smiles;
  OBMol* stored_pmol=NULL;
  std::vector<OBBase*>::iterator iter;
  for(iter= vec.begin();iter!=vec.end();++iter)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(*iter);
    if(!pmol)
      continue;
    smiles = smconv.WriteString(pmol);
    Trim(smiles);

    if(stored_smiles==smiles)
    {
      //add the coordinates of the current mol to the stored one as a conformer, and delete current mol
      double *confCoord = new double [pmol->NumAtoms() * 3];
      memcpy((char*)confCoord,(char*)pmol->GetCoordinates(),sizeof(double)*3*pmol->NumAtoms());
      stored_pmol->AddConformer(confCoord);
      delete pmol;
      *iter = NULL;
    }
    else
    {
      stored_pmol = pmol;
      stored_smiles = smiles;
    }
  }

  //erase the NULLS
  vec.erase(std::remove(vec.begin(),vec.end(), (void*)NULL), vec.end());
  return true;
}

} //namespace

/*
    To use with OBConversion::Read(), etc.
    OBMol mol;
    vector<OBBase*> vec;
    while(pConv->Read(&mol))
      vec.push_back(&mol);
    OBOp* pOp = OBOp::FindType("readconformers");
    if(!pOp)
      pOp->ProcessVec(vec);
    vec now contains one or more molecules with multiple conformers
*/


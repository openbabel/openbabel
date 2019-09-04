/**********************************************************************
canonical.cpp - A OBOp for generation for renumbering atoms canonical

Copyright (C) 2010 by Tim Vandermeersch

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
#include <openbabel/obiter.h>
#include <openbabel/graphsym.h>
#include <openbabel/canon.h>

namespace OpenBabel
{

class OpCanonical : public OBOp
{
public:
  OpCanonical(const char* ID) : OBOp(ID, false){};
  const char* Description(){ return "Canonicalize the atom order"; }

  virtual bool WorksWith(OBBase* pOb)const{ return dynamic_cast<OBMol*>(pOb)!=NULL; }
  virtual bool Do(OBBase* pOb, const char* OptionText=NULL, OpMap* pOptions=NULL, OBConversion* pConv=NULL);
};

/////////////////////////////////////////////////////////////////
OpCanonical theOpCanonical("canonical"); //Global instance

/////////////////////////////////////////////////////////////////
bool OpCanonical::Do(OBBase* pOb, const char* OptionText, OpMap* pOptions, OBConversion* pConv)
{
  OBMol* pmol = dynamic_cast<OBMol*>(pOb);
  if(!pmol)
    return false;

  std::vector<OBAtom*> atoms;
  FOR_ATOMS_OF_MOL (atom, pmol)
    atoms.push_back(&*atom);

  std::vector<unsigned int> symmetry_classes;
  OBGraphSym gs(pmol);
  gs.GetSymmetry(symmetry_classes);

  std::vector<unsigned int> canon_labels;
  CanonicalLabels(pmol, symmetry_classes, canon_labels);

  std::vector<OBAtom*> newatoms(atoms.size(), 0);
  for (std::size_t i = 0; i < canon_labels.size(); ++i)
    newatoms[canon_labels[i]-1] = atoms[i];

  pmol->RenumberAtoms(newatoms);

  return true;
}
}//namespace

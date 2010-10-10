/**********************************************************************
fillUC.cpp - plugin to fill the unit cell from unique atom positions,
using unit cell & spacegroup operations

Copyright (C) 2009 Vincent Favre-Nicolin

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
#include <openbabel/math/spacegroup.h>
#include <openbabel/generic.h>
#include <openbabel/obconversion.h>
#include <map>
#include <vector>
#include <iostream>

namespace OpenBabel
{

class OpFillUC : public OBOp
{
public:
  OpFillUC(const char* ID) : OBOp(ID, false){
    OBConversion::RegisterOptionParam("fillUC", NULL, 1, OBConversion::GENOPTIONS);
  }
  const char* Description(){ return "<param> Fill the unit cell (strict or keepconnect)\n"
    "using unique positions, unit cell and spacegroup"
    "<param> can be:\n"
    "   strict (keep only atoms inside the UC) => use \"--fillUC strict\"\n"
    "   keepconnect (fill the unit cell but keep the original connectivity => use \"--fillUC keepconnect\""; }

  virtual bool WorksWith(OBBase* pOb)const{ return dynamic_cast<OBMol*>(pOb)!=NULL; }
  virtual bool Do(OBBase* pOb, const char* OptionText=NULL, OpMap* pOptions=NULL, OBConversion* pConv=NULL);
};

/////////////////////////////////////////////////////////////////
OpFillUC theOpFillUC("fillUC"); //Global instance


// Helper function -- transform fractional coordinates to ensure they lie in the unit cell
//
// Copied from generic.cpp - should be defined in some header ? vector.h ?
vector3 transformedFractionalCoordinate2(vector3 originalCoordinate)
{
  // ensure the fractional coordinate is entirely within the unit cell
  vector3 returnValue(originalCoordinate);

  // So if we have -2.08, we take -2.08 - (-2) = -0.08 .... almost what we want
  returnValue.SetX(originalCoordinate.x() - int(originalCoordinate.x()) );
  returnValue.SetY(originalCoordinate.y() - int(originalCoordinate.y()) );
  returnValue.SetZ(originalCoordinate.z() - int(originalCoordinate.z()) );

  if (returnValue.x() < 0.0)
  returnValue.SetX(returnValue.x() + 1.0);
  if (returnValue.y() < 0.0)
  returnValue.SetY(returnValue.y() + 1.0);
  if (returnValue.z() < 0.0)
  returnValue.SetZ(returnValue.z() + 1.0);

  return returnValue;
}

/////////////////////////////////////////////////////////////////
bool OpFillUC::Do(OBBase* pOb, const char* OptionText, OpMap* pOptions, OBConversion* pConv)
{
  OBMol* pmol = dynamic_cast<OBMol*>(pOb);
  if(!pmol)
    return false;

  if (!(pmol->HasData(OBGenericDataType::UnitCell)))
  {
    obErrorLog.ThrowError(__FUNCTION__, "Cannot fill unit cell without a unit cell !" , obWarning);
    return false;
  }
  OBUnitCell *pUC = (OBUnitCell*)pmol->GetData(OBGenericDataType::UnitCell);
  const SpaceGroup* pSG = pUC->GetSpaceGroup();
  if (pSG == NULL)
  {
    obErrorLog.ThrowError(__FUNCTION__, "Cannot fill unit cell without spacegroup information !" , obWarning);
    return false;
  }
  // Now loop over all symmetry operations, and generate symmetric atoms one at a time
  // Avoid creating overlapping atoms (duplicate), and bring back atoms within the unit cell
  // using two options:
  // "--fillUC strict": keep only atoms that are strictly inside the unit cell
  //                    (fractionnal coordinates 0<= <1)
  // "--fillUC keepconnect": generate symmetrics of the molecule, and translate
  //                         it back in the unit cell if necessary

  std::map<OBAtom*,std::vector<vector3> > vatoms;// key: original atoms, value=all generated symmetrics
  FOR_ATOMS_OF_MOL(atom, *pmol)
      vatoms[&(*atom)]=std::vector<vector3>();

  for(std::map<OBAtom*,std::vector<vector3> >:: iterator atom=vatoms.begin();
      atom!=vatoms.end();++atom){
    vector3 orig = atom->first->GetVector();
    orig = pUC->CartesianToFractional(orig);// To fractionnal coordinates

    // Loop over symmetry operators
    transform3dIterator ti;
    const transform3d *t = pSG->BeginTransform(ti);
    while(t){
      atom->second.push_back ( (transform3d)(*t) * orig);
      t = pSG->NextTransform(ti);
    }
  }
  if(0==strncasecmp(OptionText, "keepconnect", 11)){
    // First, bring back all symmetrical molecules back in the UC
    for(unsigned int i=0;i<vatoms.begin()->second.size();++i){
      vector3 ccoord(0,0,0);//geometrical center
      for(std::map<OBAtom*,std::vector<vector3> >:: iterator atom=vatoms.begin();
        atom!=vatoms.end();++atom){
        ccoord+=atom->second[i];
      }
      ccoord/=vatoms.size();
      ccoord=transformedFractionalCoordinate2(ccoord)-ccoord;
      for(std::map<OBAtom*,std::vector<vector3> >:: iterator atom=vatoms.begin();
        atom!=vatoms.end();++atom){
        atom->second[i]+=ccoord;
      }
    }
    // Now add atoms that are not duplicates
    for(std::map<OBAtom*,std::vector<vector3> >:: iterator atom=vatoms.begin();
        atom!=vatoms.end();++atom){
      for(unsigned int i=1;i<atom->second.size();++i){
        bool foundDuplicate = false;
        for(unsigned int j=0;j<i;++j){
          if(atom->second[i].distSq(atom->second[j])<1e-4){
            foundDuplicate=true;
            break;
          }
        }
        if(!foundDuplicate){
          OBAtom *newAtom = pmol->NewAtom();
          newAtom->Duplicate(atom->first);
          newAtom->SetVector( pUC->FractionalToCartesian(atom->second[i]));
        }
      }
    }
  }
  else{
    if(0!=strncasecmp(OptionText, "strict", 6))
      obErrorLog.ThrowError(__FUNCTION__, "fillUC: lacking \"strict\n or \"keepconnect\" option, using strict" , obWarning);
    for(std::map<OBAtom*,std::vector<vector3> >:: iterator atom=vatoms.begin();
        atom!=vatoms.end();++atom){
      // Bring back within unit cell
      for(unsigned int i=0;i<atom->second.size();++i){
        atom->second[i]=transformedFractionalCoordinate2(atom->second[i]);
      }
      for(unsigned int i=1;i<atom->second.size();++i){
        bool foundDuplicate = false;
        for(unsigned int j=0;j<i;++j){
          if(atom->second[i].distSq(atom->second[j])<1e-4){
            foundDuplicate=true;
            break;
          }
        }
        if(!foundDuplicate){
          OBAtom *newAtom = pmol->NewAtom();
          newAtom->Duplicate(atom->first);
          newAtom->SetVector( pUC->FractionalToCartesian(atom->second[i]));
        }
      }
    }
  }

  // Set spacegroup to P1, since we generated all symmetrics
  pUC->SetSpaceGroup("P1");
/*
  list<vector3> transformedVectors; // list of symmetry-defined copies of the atom
  vector3 uniqueV, newV, updatedCoordinate;
    list<vector3> coordinates; // all coordinates to prevent duplicates

    vector3 uniqueV, newV, updatedCoordinate;
    list<vector3> transformedVectors; // list of symmetry-defined copies of the atom
    list<vector3>::iterator transformIterator, duplicateIterator;
    OBAtom *newAtom;
    list<OBAtom*> atoms; // keep the current list of unique atoms -- don't double-create
    list<vector3> coordinates; // all coordinates to prevent duplicates
    bool foundDuplicate;
    FOR_ATOMS_OF_MOL(atom, *mol)
      atoms.push_back(&(*atom));

    list<OBAtom*>::iterator i;
    for (i = atoms.begin(); i != atoms.end(); ++i) {
      uniqueV = (*i)->GetVector();
      uniqueV = CartesianToFractional(uniqueV);
      uniqueV = transformedFractionalCoordinate(uniqueV);
      coordinates.push_back(uniqueV);

      transformedVectors = sg->Transform(uniqueV);
      for (transformIterator = transformedVectors.begin();
           transformIterator != transformedVectors.end(); ++transformIterator) {
        // coordinates are in reciprocal space -- check if it's in the unit cell
        // if not, transform it in place
        updatedCoordinate = transformedFractionalCoordinate(*transformIterator);
        foundDuplicate = false;

        // Check if the transformed coordinate is a duplicate of an atom
        for (duplicateIterator = coordinates.begin();
             duplicateIterator != coordinates.end(); ++duplicateIterator) {
          if (duplicateIterator->distSq(updatedCoordinate) < 1.0e-4) {
            foundDuplicate = true;
            break;
          }
        }
        if (foundDuplicate)
          continue;

        coordinates.push_back(updatedCoordinate); // make sure to check the new atom for dupes
        newAtom = mol->NewAtom();
        newAtom->Duplicate(*i);
        newAtom->SetVector(FractionalToCartesian(updatedCoordinate));
      } // end loop of transformed atoms
      (*i)->SetVector(FractionalToCartesian(uniqueV));
    } // end loop of atoms
*/
  return true;
}
}//namespace

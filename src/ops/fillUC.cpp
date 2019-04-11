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
#include <openbabel/atom.h>
#include <openbabel/obiter.h>
#include <openbabel/math/spacegroup.h>
#include <openbabel/generic.h>
#include <openbabel/obconversion.h>
#include <map>
#include <vector>
#include <iostream>

using namespace std;
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

// Wrap coordinates in the unit cell with some fuzziness, i.e. when one
// coordinate is very very close to 1 (>= 0.999999), we wrap it to
// exactly zero.
vector3 fuzzyWrapFractionalCoordinate (vector3 coord)
{
    double x = fmod(coord.x(), 1);
    double y = fmod(coord.y(), 1);
    double z = fmod(coord.z(), 1);
    if (x < 0) x += 1;
    if (y < 0) y += 1;
    if (z < 0) z += 1;

#define LIMIT 0.999999
    if (x > LIMIT)
      x -= 1;
    if (y > LIMIT)
      y -= 1;
    if (z > LIMIT)
      z -= 1;
#undef LIMIT

    // Fuzzy logic from Francois-Xavier
#define EPSILON 1.0e-6
    if (x > 1 - EPSILON || x < EPSILON)
      x = 0.0;
    if (y > 1 - EPSILON || y < EPSILON)
      y = 0.0;
    if (z > 1 - EPSILON || z < EPSILON)
      z = 0.0;
#undef EPSILON

    return vector3(x, y, z);
}

// Whether two points (given in fractional coordinates) are close enough
// to be considered duplicates.
// This function is duplicate from generic.cpp, these should be merged
bool areDuplicateAtoms2(vector3 v1, vector3 v2)
{
  vector3 dr = fuzzyWrapFractionalCoordinate(v2)
    - fuzzyWrapFractionalCoordinate(v1);

  if (dr.x() < -0.5)
    dr.SetX(dr.x() + 1);
  if (dr.x() > 0.5)
    dr.SetX(dr.x() - 1);
  if (dr.y() < -0.5)
    dr.SetY(dr.y() + 1);
  if (dr.y() > 0.5)
    dr.SetY(dr.y() - 1);
  if (dr.z() < -0.5)
    dr.SetZ(dr.z() + 1);
  if (dr.z() > 0.5)
    dr.SetZ(dr.z() - 1);

  return (dr.length_2() < 1e-3);
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
  SpaceGroup spacegroup;
  const SpaceGroup* pSG;
  map<string,string>::const_iterator itr;

  if(pOptions && pOptions->find("transformations") != pOptions->end())
  {
    itr = pOptions->find("transformations");
    vector<string> vec;
    tokenize(vec, itr->second.c_str());
    for(vector<string>::iterator iter = vec.begin(); iter != vec.end(); ++ iter)
    {
      if (iter == vec.begin()) // Warn user about converting only once
        obErrorLog.ThrowError(__FUNCTION__, "Converting to P 1 cell using available symmetry transformations." , obWarning);
      spacegroup.AddTransform(iter->c_str());
    }

    pSG = &spacegroup;
  }
  else
  {
    pSG = pUC->GetSpaceGroup();
  }

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
    orig = pUC->CartesianToFractional(orig);// To fractional coordinates

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
      ccoord /=vatoms.size();
      ccoord=fuzzyWrapFractionalCoordinate(ccoord)-ccoord;
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
          if(areDuplicateAtoms2(atom->second[i],atom->second[j])){
            foundDuplicate=true;
            break;
          }
        }
        if(!foundDuplicate){
          vector3 transformed = pUC->FractionalToCartesian(atom->second[i]);
          // let's make sure there isn't some *other* atom that's in this spot
          bool foundCartesianDuplicate = false;
          FOR_ATOMS_OF_MOL(a, *pmol) {
            vector3 diff = a->GetVector() - transformed;
            if (diff.length_2() < 1.0e-4) {
              foundCartesianDuplicate = true;
              break;
            }
          }

          if (!foundCartesianDuplicate) {
            OBAtom *newAtom = pmol->NewAtom();
            newAtom->Duplicate(atom->first);
            newAtom->SetVector( transformed );
          }
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
        atom->second[i]=fuzzyWrapFractionalCoordinate(atom->second[i]);
      }
      for(unsigned int i=1;i<atom->second.size();++i){
        bool foundDuplicate = false;
        for(unsigned int j=0;j<i;++j){
          if(areDuplicateAtoms2(atom->second[i],atom->second[j])){
            foundDuplicate=true;
            break;
          }
        }
        if(!foundDuplicate){
          vector3 transformed = pUC->FractionalToCartesian(atom->second[i]);
          // let's make sure there isn't some *other* atom that's in this spot
          bool foundCartesianDuplicate = false;
          FOR_ATOMS_OF_MOL(a, *pmol) {
            vector3 diff = a->GetVector() - transformed;
            if (diff.length_2() < 1.0e-4) {
              foundCartesianDuplicate = true;
              break;
            }
          }

          if (!foundCartesianDuplicate) {
            OBAtom *newAtom = pmol->NewAtom();
            newAtom->Duplicate(atom->first);
            newAtom->SetVector( transformed );
          }
        }
      }
    }
  }

  // Set spacegroup to P1, since we generated all symmetrics
  pUC->SetSpaceGroup("P1");
  return true;
}
}//namespace

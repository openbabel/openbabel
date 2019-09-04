/**********************************************************************
opalign.cpp - Align substructures in multiple molecules

Copyright (C) 2010 by Chris Morley
 
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
#include <openbabel/babelconfig.h>
#include <iostream>
#include<openbabel/op.h>
#include<openbabel/mol.h>
#include<openbabel/math/align.h> // ** This requires Eigen to be installed **
#include<openbabel/obconversion.h>
#include <openbabel/parsmart.h>
#include <openbabel/generic.h>
#include "opisomorph.h"

//#define DEPICTION2D     0x100 // Should have the same value as in format.h!!

namespace OpenBabel
{
using namespace std;

class OpAlign : public OBOp
{
public:
  OpAlign(const char* ID) : OBOp(ID, false), _align(false, false){};
  const char* Description(){ return 
    "Align coordinates to the first molecule\n"
    "Typical use with a -s option:\n"
    "    obabel pattern.www  dataset.xxx  -O outset.yyy  -s SMARTS  --align\n"
    "Only molecules matching SMARTS are converted and are aligned by\n"
    "having all their atom coordinates modified. The atoms that are\n"
    "used in the alignment are those matched by SMARTS in the first\n"
    "output molecule. The subsequent molecules are aligned so that\n"
    "the coordinates of atoms equivalent to these are as nearly as\n"
    "possible the same as those of the pattern atoms.\n"
    "The atoms in the various molecules can be in any order.\n"
    "Tha alignment ignores hydrogen atoms but includes symmetry\n"
    "The standalone program obfit has similar functionality.\n \n"

    "The first input molecule could be part of the data set :\n"
    "    obabel dataset.xxx  -O outset.yyy  -s SMARTS  --align\n"
    "This form also ensures that a particular substructure always\n"
    "has the same orientation in a 2D display of a set of molecules.\n"
    "0D molecules, e.g. from SMILES, are given 2D coordinates before\n"
    "alignment.\n \n"

    "See documentation for the -s option for its other possible\n"
    "parameters. For example, the matching atoms could be those\n"
    "of a molecule in a specified file.\n \n"
    
    "Without an -s option, all the atoms in the first molecule\n"
    "are used as pattern atoms. The order of the atoms must be the same\n"
    "in all the molecules.\n\n"

     "The output molecules have a property (represented internally as\n"
     "OBPairData) called ``rmsd``, which is a measure of the quality\n"
     "of the fit. To attach it to the title of each molecule use\n"
     "--append rmsd.\n"
     "To output the two conformers closest to the first conformer in a dataset:\n"
     "    obabel dataset.xxx  -O outset.yyy  --align  --smallest 2 rmsd\n\n"
    ;

}
  virtual bool WorksWith(OBBase* pOb)const{ return dynamic_cast<OBMol*>(pOb)!=NULL; }
  virtual bool Do(OBBase* pOb, const char* OptionText=NULL, OpMap* pOptions=NULL, OBConversion* pConv=NULL);
private:
  OBAlign _align;
  OBMol _refMol;
  std::vector<vector3> _refvec;
  OpNewS* _pOpIsoM;  //the address of the -s option or NULL if it is not used
  std::string _stext;//the -s option parameters
};

/////////////////////////////////////////////////////////////////
OpAlign theSecondOpAlign("align");

/////////////////////////////////////////////////////////////////
bool OpAlign::Do(OBBase* pOb, const char* OptionText, OpMap* pmap, OBConversion* pConv)
{
  OBMol* pmol = dynamic_cast<OBMol*>(pOb);
  if(!pmol)
    return false;
    
  map<string,string>::const_iterator itr;

  // Is there an -s option?
  if(pConv->IsFirstInput())
  {
    _pOpIsoM = NULL; //assume no -s option
    itr = pmap->find("s");
    if(itr!=pmap->end())
    {
      //There is an -s option; check it is ok
      _pOpIsoM = static_cast<OpNewS*>(OBOp::FindType("s"));
      _stext = itr->second; //get its parameter(s)
      if(!_pOpIsoM || _stext.empty())
      {
        obErrorLog.ThrowError(__FUNCTION__, 
        "No parameter on -s option, or its OBOp version is not loaded", obError);
        pConv->SetOneObjectOnly(); //to finish
        return false;
      }
    }
  }

  // If the output format is a 2D depiction format, then we should align
  // on the 2D coordinates and not the 3D coordinates (if present). This
  //means we need to generate the 2D coordinates at this point.
  if(pmol->GetDimension()==3 && (pConv->GetOutFormat()->Flags() & DEPICTION2D))
  {
    OBOp* pgen = OBOp::FindType("gen2D");
    if(pgen)
      pgen->Do(pmol);
  }

  // All molecules must have coordinates, so add them if 0D
  // They may be added again later when gen2D or gen3D is called, but they will be the same.
  // It would be better if this op was called after them, which would happen
  // if its name was alphabetically after "gen" (and before "s").
  if(pmol->GetDimension()==0)
  {
    //Will the coordinates be 2D or 3D?
    itr = pmap->find("gen3D");
    OBOp* pgen = (itr==pmap->end()) ? OBOp::FindType("gen2D") : OBOp::FindType("gen3D");
    if(pgen)
      pgen->Do(pmol);
  }

  //Do the alignment in 2D if the output format is svg, png etc. and there is no -xn option
  if(pmol->GetDimension()==3 && pConv && !pConv->IsOption("n"))
  {
    OBFormat* pOutFormat = pConv->GetOutFormat();
    if(pOutFormat->Flags() & DEPICTION2D)
    {
      OBOp* pgen = OBOp::FindType("gen2D");
      if(pgen)
        pgen->Do(pmol);
    }
  }

  if(pConv->IsFirstInput() || _refMol.NumAtoms()==0)
  {
    _refvec.clear();
    // Reference molecule is basically the first molecule
    _refMol = *pmol;
    if(!_pOpIsoM)
     //no -s option. Use a molecule reference.
     _align.SetRefMol(_refMol);
    else
    {
      //If there is a -s option, reference molecule has only those atoms that are matched
      //Call the -s option from here
      bool ret = _pOpIsoM->Do(pmol, _stext.c_str(), pmap, pConv);
      // Get the atoms that were matched       
      vector<int> ats = _pOpIsoM->GetMatchAtoms();
      if(!ats.empty())
      {
        // Make a vector of the matching atom coordinates...
        for(vector<int>::iterator iter=ats.begin(); iter!=ats.end(); ++iter)
          _refvec.push_back((pmol->GetAtom(*iter))->GetVector());        
        // ...and use a vector reference
        _align.SetRef(_refvec);
      }
      // Stop -s option being called normally, although it will still be called once
      //  in the DoOps loop already started for the current (first) molecule.
      pConv->RemoveOption("s",OBConversion::GENOPTIONS);
      if(!ret) 
      {
        // the first molecule did not match the -s option so a reference molecule
        // could not be made. Keep trying.
        _refMol.Clear();
        //obErrorLog.ThrowError(__FUNCTION__, "The first molecule did not match the -s option\n"
        //  "so the reference structure was not derived from it", obWarning, onceOnly);
        return false; //not matched
      }
    }
  }

  //All molecules
  if(pmol->GetDimension()!= _refMol.GetDimension())
  {
    stringstream ss;
    ss << "The molecule" << pmol->GetTitle()
       << " does not have the same dimensions as the reference molecule "
       << _refMol.GetTitle() << " and is ignored.";
       obErrorLog.ThrowError(__FUNCTION__, ss.str().c_str(), obError);
    return false;
  }

  if(_pOpIsoM) //Using -s option
  {   
    //Ignore mol if it does not pass -s option
    if(!_pOpIsoM->Do(pmol, "", pmap, pConv)) // "" means will use existing parameters
      return false;

    // Get the atoms equivalent to those in ref molecule        
    vector<int> ats = _pOpIsoM->GetMatchAtoms();

    // Make a vector of their coordinates and get the centroid
    vector<vector3> vec;
    vector3 centroid;
    for(vector<int>::iterator iter=ats.begin(); iter!=ats.end(); ++iter) {
      vector3 v = pmol->GetAtom(*iter)->GetVector();
      centroid += v;
      vec.push_back(v);
    }
    centroid /= vec.size();
    
    // Do the alignment
    _align.SetTarget(vec);
    if(!_align.Align())
      return false;

    // Get the centroid of the reference atoms
    vector3 ref_centroid;
    for(vector<vector3>::iterator iter=_refvec.begin(); iter!=_refvec.end(); ++iter)
      ref_centroid += *iter;
    ref_centroid /= _refvec.size();

    //subtract the centroid, rotate the target molecule, then add the centroid
    matrix3x3 rotmatrix = _align.GetRotMatrix();
    for (unsigned int i = 1; i <= pmol->NumAtoms(); ++i)
    {
      vector3 tmpvec = pmol->GetAtom(i)->GetVector();
      tmpvec -= centroid;
      tmpvec *= rotmatrix; //apply the rotation
      tmpvec += ref_centroid;
      pmol->GetAtom(i)->SetVector(tmpvec);
    }
  }
  else //Not using -s option)
  {
    _align.SetTargetMol(*pmol);
    if(!_align.Align())
      return false;
    _align.UpdateCoords(pmol);
  }

  //Save rmsd as a property
  OBPairData* dp = new OBPairData;
  dp->SetAttribute("rmsd");
  double val = _align.GetRMSD();
  if(val<1e-12)
    val = 0.0;
  dp->SetValue(toString(val));
  dp->SetOrigin(local);
  pmol->SetData(dp);

  return true;
}


}//namespace
/*
With a -s option
Use the atoms from the match of the first molecule to make a vector of their coordinates.
This will be in the order of the match.
For subsequent molecules do the same, and do the alignment with two vectors. The order
of the atoms in the target molecules should not matter.

Without a -s option
Do the alignment with two molecules. The order of the atoms must be the same.
*/

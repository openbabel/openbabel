/**********************************************************************
builder.cpp - Class to create structures.

Copyright (C) 2007-2008 by Tim Vandermeersch 
                           <tim.vandermeersch@gmail.com>
 
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

#include <openbabel/builder.h>

#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obiter.h>
#include <openbabel/math/matrix3x3.h>
#include <openbabel/rotamer.h>
#include <openbabel/rotor.h>
#include <openbabel/obconversion.h>

/* OBBuilder::GetNewBondVector():
 * - is based on OBAtom::GetNewBondVector()
 * - but: when extending a long chain all the bonds are trans
 *
 * fragments.txt:
 * - fragments should be ordered from complex to simple
 *
 * */

using namespace std;

namespace OpenBabel
{
  /** \class OBBuilder builder.h <openbabel/builder.h>
      \brief Class for 3D structure generation

      The OBBuilder class is used for generating 3D structures.

      Below is and example which explain the basics. 
      
      \code
      //
      // code to read molecule from smiles goes here...
      //
      OBBuilder builder;
      builder.Build(mol);
      //
      // code to write molecule to 3D file format goes here...
      //
      \endcode
  **/
  std::vector<std::pair<OBSmartsPattern*, std::vector<vector3> > > OBBuilder::_fragments;

  void OBBuilder::LoadFragments()  {
    // open data/fragments.txt
    ifstream ifs;
    if (OpenDatafile(ifs, "fragments.txt").length() == 0) {
      obErrorLog.ThrowError(__FUNCTION__, "Cannot open fragments.txt", obError);
      return;
    }

    char buffer[BUFF_SIZE];
    vector<string> vs;
    OBSmartsPattern *sp = NULL;
    vector<vector3> coords;
    while (ifs.getline(buffer, BUFF_SIZE)) {
      if (buffer[0] == '#') // skip comment line (at the top)
        continue;
        
      tokenize(vs, buffer);
      
      if (vs.size() == 1) { // SMARTS pattern
        if (sp != NULL)
          _fragments.push_back(pair<OBSmartsPattern*, vector<vector3> > (sp, coords));
        
        coords.clear();
        sp = new OBSmartsPattern;
        if (!sp->Init(vs[0])) {
          delete sp;
          sp = NULL;
          obErrorLog.ThrowError(__FUNCTION__, " Could not parse SMARTS from contribution data file", obInfo);
        }
      } else if (vs.size() == 3) { // XYZ coordinates
        vector3 coord(atof(vs[0].c_str()), atof(vs[1].c_str()), atof(vs[2].c_str()));
        coords.push_back(coord);
      }
 
    }
    
    _fragments.push_back(pair<OBSmartsPattern*, vector<vector3> > (sp, coords));
  }
 
  vector3 OBBuilder::GetNewBondVector(OBAtom *atom){
    vector3 bond1, bond2, bond3, v1, v2, newbond;
    
    bond1 = VZero;
    bond2 = VZero;
    bond3 = VZero;
    
    if (atom == NULL)
      return VZero;

    if (atom->GetValence() == 0) {
      newbond = atom->GetVector() + VX * 1.5;
      return newbond;
    }

    if (atom->GetValence() == 1) {
      FOR_NBORS_OF_ATOM (nbr, atom) {
        bond1 = atom->GetVector() - nbr->GetVector();   
        
        FOR_NBORS_OF_ATOM (nbr2, &*nbr)
          if (&*nbr2 != atom)
            bond2 = nbr->GetVector() - nbr2->GetVector();
      }

      bond1 = bond1.normalize();
      if (bond2 == VZero) {
        if (atom->GetHyb() == 1)
          newbond = bond1;
        if (atom->GetHyb() == 2)
          newbond = bond1 + VY * tan(DEG_TO_RAD*60);
        if (atom->GetHyb() == 3)
          newbond = bond1 + VY * tan(DEG_TO_RAD*70.5);
        
        newbond = newbond.normalize();
        newbond *= 1.5;
        newbond += atom->GetVector();
        return newbond;
      } else {
        v1 = cross(bond1, bond2);
        v2 = cross(bond1, v1);
        
        if (atom->GetHyb() == 1)
          newbond = bond1;
        if (atom->GetHyb() == 2)
          newbond = bond1 - v2 * tan(DEG_TO_RAD*60);
        if (atom->GetHyb() == 3)
          newbond = bond1 - v2 * tan(DEG_TO_RAD*70.5);
        
        newbond = newbond.normalize();
        newbond *= 1.5;
        newbond += atom->GetVector();
        return newbond;
      }
    }
    
    if (atom->GetValence() == 2) {
      FOR_NBORS_OF_ATOM (nbr, atom) {
        if (bond1 == VZero)
          bond1 = atom->GetVector() - nbr->GetVector();
        else
          bond2 = atom->GetVector() - nbr->GetVector();
      }

      v1 = bond1 + bond2;
      v1 = v1.normalize();
     
      if (atom->GetHyb() == 2)
        newbond = v1;
      if (atom->GetHyb() == 3) {
        v2 = cross(bond1, bond2);
        v1 = v1.normalize();
        newbond = v2 + v1 * tan(DEG_TO_RAD*35.25);
      }
      
      newbond = newbond.normalize();
      newbond *= 1.5;
      newbond += atom->GetVector();
      return newbond;
    }
    
    if (atom->GetValence() == 3) {
      FOR_NBORS_OF_ATOM (nbr, atom) {
        if (bond1 == VZero)
          bond1 = atom->GetVector() - nbr->GetVector();
        else if (bond2 == VZero)
          bond2 = atom->GetVector() - nbr->GetVector();
        else
          bond3 = atom->GetVector() - nbr->GetVector();
      }
          
      newbond = bond1 + bond2 + bond3;
      newbond = newbond.normalize();
      newbond *= 1.5;
      newbond += atom->GetVector();
      return newbond;
    }

    return VZero; //previously undefined
  }
  
  bool OBBuilder::Connect(OBMol &mol, int idxA, int idxB, int bondOrder)
  {
    _workMol = mol;
    OBAtom *a = _workMol.GetAtom(idxA);
    OBAtom *b = _workMol.GetAtom(idxB);
    
    if (a == NULL)
      return false;
    if (b == NULL)
      return false;
    
    if (!Connect(a, b))
      return false;

    _workMol.GetBond(idxA, idxB)->SetBondOrder(bondOrder);
    
    // copy back
    mol = _workMol;
    return true;
  }

  // The OBMol mol contains both the molecule to which we want to connect the 
  // fragment and the fragment itself. The fragment containing b will be 
  // rotated and translated. Atom a is the atom from 
  // the main molecule to which we want to connect atom b.
  bool OBBuilder::Connect(OBAtom *a, OBAtom *b)
  {
    OBBitVec fragment = GetFragment(b->GetIdx());
    if (fragment == GetFragment(a->GetIdx()))
      return false; // a and b are in the same fragment

    // Make sure we use a and be from _workMol
    a = _workMol.GetAtom(a->GetIdx());
    b = _workMol.GetAtom(b->GetIdx());
    vector3 posa = a->GetVector();
    vector3 posb = b->GetVector();
    vector3 newpos = GetNewBondVector(a);
    vector3 moldir = newpos - posa;
    // 
    // translate fragment so that atom b is at the origin
    //
    for (unsigned int i = 1; i <= _workMol.NumAtoms(); ++i) {
      if (fragment.BitIsSet(i)) {
        // the atom is part of the fragment, translate it
        vector3 tmpvec = _workMol.GetAtom(i)->GetVector();
        tmpvec -= posb;
        _workMol.GetAtom(i)->SetVector(tmpvec);
      }
    }
    // 
    // rotate the fragment to align the bond directions Mol-a-b and a-b-fragment
    //  
    matrix3x3 xymat, xzmat, yzmat;
    double xyang, yzang, xzang;

    vector3 fragdir = GetNewBondVector(b); // b is at origin 
    xyang = vectorAngle(vector3(moldir.x(), moldir.y(), 0.0), vector3(fragdir.x(), fragdir.y(), 0.0));
    if (cross(vector3(moldir.x(), moldir.y(), 0.0), vector3(fragdir.x(), fragdir.y(), 0.0)).z() > 0) {
      xyang = 180 + xyang;
    } else if (cross(vector3(moldir.x(), moldir.y(), 0.0), vector3(fragdir.x(), fragdir.y(), 0.0)).z() < 0) {
      xyang = 180 - xyang;
    } else {
      xyang = 0.0;
    }
    xymat.SetupRotMat(0.0, 0.0, xyang); 
    for (unsigned int i = 1; i <= _workMol.NumAtoms(); ++i) {
      if (fragment.BitIsSet(i)) {
        vector3 tmpvec = _workMol.GetAtom(i)->GetVector();
        tmpvec *= xymat; //apply the rotation
        _workMol.GetAtom(i)->SetVector(tmpvec);
      }
    }
    
    fragdir = GetNewBondVector(b);
    xzang = vectorAngle(vector3(moldir.x(), moldir.z(), 0.0), vector3(fragdir.x(), fragdir.z(), 0.0));
    if (cross(vector3(moldir.x(), moldir.z(), 0.0), vector3(fragdir.x(), fragdir.z(), 0.0)).z() > 0) {
      xzang = 180 - xzang;
    } else if (cross(vector3(moldir.x(), moldir.z(), 0.0), vector3(fragdir.x(), fragdir.z(), 0.0)).z() < 0) {
      xzang = 180 + xzang;
    } else {
      xzang = 0.0;
    }
    xzmat.SetupRotMat(0.0, xzang, 0.0); 
    for (unsigned int i = 1; i <= _workMol.NumAtoms(); ++i) {
      if (fragment.BitIsSet(i)) {
        vector3 tmpvec = _workMol.GetAtom(i)->GetVector();
        tmpvec *= xzmat; //apply the rotation
        _workMol.GetAtom(i)->SetVector(tmpvec);
      }
    }

    fragdir = GetNewBondVector(b);
    yzang = vectorAngle(vector3(moldir.y(), moldir.z(), 0.0), vector3(fragdir.y(), fragdir.z(), 0.0));
    if (cross(vector3(moldir.y(), moldir.z(), 0.0), vector3(fragdir.y(), fragdir.z(), 0.0)).z() > 0) {
      yzang = 180 + yzang;
    } else if (cross(vector3(moldir.y(), moldir.z(), 0.0), vector3(fragdir.y(), fragdir.z(), 0.0)).z() < 0) {
      yzang = 180 - yzang;
    } else {
      yzang = 0.0;
    }
    yzmat.SetupRotMat(yzang, 0.0, 0.0); 
    for (unsigned int i = 1; i <= _workMol.NumAtoms(); ++i) {
      if (fragment.BitIsSet(i)) {
        vector3 tmpvec = _workMol.GetAtom(i)->GetVector();
        tmpvec *= yzmat; //apply the rotation
        _workMol.GetAtom(i)->SetVector(tmpvec);
      }
    }
    // 
    // translate fragment 
    //
    for (unsigned int i = 1; i <= _workMol.NumAtoms(); ++i) {
      if (fragment.BitIsSet(i)) {
        // translate the fragment
        vector3 tmpvec = _workMol.GetAtom(i)->GetVector();
        tmpvec += newpos;
        _workMol.GetAtom(i)->SetVector(tmpvec);
      }
    }     
    //
    // Create the bond between the two fragments
    //
    OBBond *bond = _workMol.NewBond();
    bond->SetBegin(a);
    bond->SetEnd(b);
    bond->SetBondOrder(1);
    a->AddBond(bond);
    b->AddBond(bond);
  
    return true;
  }

  void OBBuilder::AddNbrs(OBBitVec &fragment, OBAtom *atom)
  {
    FOR_NBORS_OF_ATOM (nbr, atom) {
      if (!fragment.BitIsSet(nbr->GetIdx())) {
        fragment.SetBitOn(nbr->GetIdx());
        AddNbrs(fragment, &*nbr);
      }
    }
  }

  OBBitVec OBBuilder::GetFragment(int index)
  { 
    OBBitVec fragment;
    OBAtom *atom = _workMol.GetAtom(index);

    fragment.SetBitOn(index);
    AddNbrs(fragment, atom);
    
    return fragment; 
  }

  // First we find the most complex fragments in our molecule. Once we have a match,
  // vfrag is set for all the atoms in the fragment. A second match (smaller, more 
  // simple part of the 1st match) is ignored.
  //
  // Next we iterate over the atoms. There are 3 possibilities:
  // 
  // 1) The atom is already added (vdone is set) --> continue
  // 2) The atom belongs to a fragment --> a) This is the first atom to add: leave fragment as it is
  //                                       b) Not the first atom: rotate, translate and connect the fragment
  // 3) The atom doesn't belong to a fragment: a) First atom: place at origin
  //                                           b) Not first atom: Find position and place atom
  bool OBBuilder::Build(OBMol &mol)
  {
    OBBitVec vdone; // Atoms that are done, need no further manipulation.
    OBBitVec vfrag; // Atoms that are part of a fragment found in the database.
                    // These atoms have coordinates, but the fragment still has 
                    // to be rotated and translated.
    vector3 molvec, moldir;
    vector<pair<OBSmartsPattern*, vector<vector3 > > >::iterator i;
    vector<vector<int> >::iterator j;
    vector<int>::iterator k, k2, k3;
    vector<vector3>::iterator l;
    vector<vector<int> > mlist; // match list for fragments
 
    // copy the molecule to private data
    OBMol origMol = mol;
    _workMol = mol;
    
    // delete all bonds in the working molecule
    while (_workMol.NumBonds())
      _workMol.DeleteBond(_workMol.GetBond(0));
    
    _workMol.SetHybridizationPerceived();

    //datafile is read only on first use of Build()
    static bool loaded; //initially false
    if(!loaded) {
      LoadFragments();
      loaded=true;
    }
    // Loop through  the database once and assign the coordinates from
    // the first (most complex) fragment.
    for (i = _fragments.begin();i != _fragments.end();++i) {
      if (i->first != NULL && i->first->Match(origMol)) { 
        mlist = i->first->GetMapList();
          
        for (j = mlist.begin();j != mlist.end();++j) { // for all matches
          if (vfrag.BitIsSet((*j)[0])) // the found match is already added
            continue;

          int index, index2, counter = 0;
          for (k = j->begin(); k != j->end(); ++k) { // for all atoms of the fragment
            index = *k;

            if (vfrag.BitIsSet(index))
              continue;
              
            vfrag.SetBitOn(index); // set vfrag for all atoms of fragment

            // set coordinates for atoms
            OBAtom *atom = _workMol.GetAtom(index);
            atom->SetVector(i->second[counter]);
            counter++;
          }
            
          // add the bonds for the fragment
          for (k = j->begin(); k != j->end(); ++k) {
            index = *k;
            OBAtom *atom1 = origMol.GetAtom(index);
              
            for (k2 = j->begin(); k2 != j->end(); ++k2) {
              index2 = *k2;
              OBAtom *atom2 = origMol.GetAtom(index2);
              OBBond *bond = atom1->GetBond(atom2);

              if (bond != NULL) 
                _workMol.AddBond(*bond);
            }
          }
        }
      }
    }

    // iterate over all atoms to place them in 3D space
    FOR_DFS_OF_MOL (a, origMol) {
      if (vdone.BitIsSet(a->GetIdx())) // continue if the atom is already added
        continue;
      
      // find an atom connected to the current atom that is already added
      OBAtom *prev = NULL;
      FOR_NBORS_OF_ATOM (nbr, &*a) {
        if (vdone.BitIsSet(nbr->GetIdx()))
          prev = &*nbr;
      }
 
      if (vfrag.BitIsSet(a->GetIdx())) { // continue if the atom is already added
        if (prev != NULL) { // if we have a previous atom, translate/rotate the fragment and connect it
          Connect(prev, &*a);
          // set the correct bond order
          int bondOrder = origMol.GetBond(prev->GetIdx(), a->GetIdx())->GetBondOrder();
          _workMol.GetBond(prev->GetIdx(), a->GetIdx())->SetBondOrder(bondOrder);
        }
        
        OBBitVec fragment = GetFragment(a->GetIdx());
        vdone |= fragment; // mark this fragment as done

        continue;
      }

      //
      // below is the code to add non-fragment atoms
      //
   
      // get the position for the new atom, this is done with GetNewBondVector
      if (prev != NULL) {
        molvec = GetNewBondVector(_workMol.GetAtom(prev->GetIdx()));
        moldir = molvec - _workMol.GetAtom(prev->GetIdx())->GetVector();
      } else {
        molvec = VX;
        moldir = VX;
      }
      
      vdone.SetBitOn(a->GetIdx());
      
      // place the atom the atom 
      OBAtom *atom = _workMol.GetAtom(a->GetIdx());
      atom->SetVector(molvec);

      // add bond between previous part and added atom
      if (prev != NULL) {
        OBBond *bond = a->GetBond(prev); // from origMol
        _workMol.AddBond(*bond);
      }

    }

    mol = _workMol;
    mol.SetDimension(3);
    return true;
  }


} // end namespace OpenBabel


//! \file builder.cpp
//! \brief Handle OBBuilder class

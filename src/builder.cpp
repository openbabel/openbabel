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
 * - practically speaking, this means they are ordered by the number of atoms
 * */

using namespace std;

namespace OpenBabel
{
  /** \class OBBuilder builder.h <openbabel/builder.h>
      \brief Class for 3D structure generation
      \since version 2.2
      
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

    // Set the locale for number parsing to avoid locale issues: PR#1785463
    char *old_num_locale = strdup (setlocale (LC_NUMERIC, NULL));
  	setlocale(LC_NUMERIC, "C");

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

    // return the locale to the original one
  	setlocale(LC_NUMERIC, old_num_locale);
  	free (old_num_locale);
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

  bool OBBuilder::Connect(OBAtom *a, OBAtom *b)
  {
    a = _workMol.GetAtom(a->GetIdx());
    b = _workMol.GetAtom(b->GetIdx());
 
    vector3 newpos = GetNewBondVector(a);

    return Connect(a, b, newpos);
  }

  // The OBMol mol contains both the molecule to which we want to connect the 
  // fragment and the fragment itself. The fragment containing b will be 
  // rotated and translated. Atom a is the atom from 
  // the main molecule to which we want to connect atom b.
  bool OBBuilder::Connect(OBAtom *a, OBAtom *b, vector3 &newpos)
  {
    OBBitVec fragment = GetFragment(b->GetIdx());
    if (fragment == GetFragment(a->GetIdx()))
      return false; // a and b are in the same fragment

    // Make sure we use a and be from _workMol
    a = _workMol.GetAtom(a->GetIdx());
    b = _workMol.GetAtom(b->GetIdx());
    vector3 posa = a->GetVector();
    vector3 posb = b->GetVector();
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
    /* 
    matrix3x3 mat;
    vector3 moldir = newpos - posa;
    vector3 fragdir = GetNewBondVector(b); // b is at origin 
    moldir.normalize(); 
    fragdir.normalize(); 
    double angle = acos(dot(moldir, fragdir));
    vector3 axis = cross(moldir, fragdir);
    axis.normalize();

    mat.RotAboutAxisByAngle(axis, angle);
    for (unsigned int i = 1; i <= _workMol.NumAtoms(); ++i) {
      if (fragment.BitIsSet(i)) {
        vector3 tmpvec = _workMol.GetAtom(i)->GetVector();
        tmpvec *= mat; //apply the rotation
        _workMol.GetAtom(i)->SetVector(tmpvec);
      }
    }
    */
    matrix3x3 xymat, xzmat, yzmat;
    vector3 moldir = newpos - posa;
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

  bool OBBuilder::Swap(OBMol &mol, int idxA, int idxB, int idxC, int idxD)
  {
    OBAtom *a = mol.GetAtom(idxA);
    OBAtom *b = mol.GetAtom(idxB);
    OBAtom *c = mol.GetAtom(idxC);
    OBAtom *d = mol.GetAtom(idxD);
    
    // make sure the atoms exist
    if (a == NULL || b == NULL || c == NULL || d == NULL)
      return false;

    OBBond *bond1 = _workMol.GetBond(idxA, idxB);
    OBBond *bond2 = _workMol.GetBond(idxC, idxD);

    // make sure a-b and c-d are connected
    if (bond1 == NULL || bond2 == NULL)
      return false;

    // make sure the bonds are not in a ring 
    if (bond1->IsInRing() || bond2->IsInRing())
      return false;

    // delete the bonds
    mol.DeleteBond(bond1);
    mol.DeleteBond(bond2);

    // save the original positions
    vector3 posB = b->GetVector();
    vector3 posD = d->GetVector();

    // connect the fragments
    if (!Connect(a, d, posB))
      return false;
    if (!Connect(c, b, posD))
      return false;
    
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
    //cerr << "OBBuilder::Build(OBMol &mol)" << endl;
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
    if(_fragments.empty())
      LoadFragments();

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
    
    // correct the chirality
    CorrectStereoBonds(_workMol);
    CorrectStereoAtoms(_workMol);

    mol = _workMol;
    mol.SetDimension(3);

    return true;
  }

  void OBBuilder::CorrectStereoBonds(OBMol &mol)
  {
    /*
    cerr << "bond IsUp IsDown" << endl;
    FOR_BONDS_OF_MOL (bond, mol) {
      cerr << bond->GetIdx() << ": " << bond->IsUp() << " " << bond->IsDown() << " (" 
           << bond->GetBeginAtomIdx() << "-" << bond->GetEndAtomIdx() << ")" << endl;
    }
    */
 
    OBAtom *a, *b, *c, *d;
    OBBond *ab, *bc, *cd;
    vector<int> done;

    mol.DeleteData(OBGenericDataType::TorsionData); // bug #1954233
    FOR_TORSIONS_OF_MOL(t, mol) {
      a = mol.GetAtom((*t)[0] + 1);
      b = mol.GetAtom((*t)[1] + 1);
      c = mol.GetAtom((*t)[2] + 1);
      d = mol.GetAtom((*t)[3] + 1);

      ab = a->GetBond(b);
      bc = b->GetBond(c);
      cd = c->GetBond(d);

      if (bc->IsAromatic() || !bc->IsDouble())
        continue;

      if (!(ab->IsUp() || ab->IsDown()) || !(cd->IsUp() || cd->IsDown())) {
        //cerr << "no double bond stereochemistry set for bond " << bc->GetIdx() << endl;
        continue;
      }

      for (unsigned int i = 0; i < done.size(); ++i)
        if (done[i] == bc->GetIdx())
          continue;

      // Check if atom a is double bonded to an atom.
      // If so, also check if we have visited the bond.
      // This is the case in 
      //      
      //        +---- this bond should be trans, the smiles parser assigns
      //        |     both C/C bonds up. So we have to invert it
      // F/C=C/C=C/C
      //  d   u   u   (down/up)
      //
      // F/C=C/C/C=C/C  
      //  d   u d   u
      //
      // See src/formats/smilesformat.cpp for details on OB's up/down implementation.
      //
      bool invert = false;
      FOR_BONDS_OF_ATOM (bond, a) {
        if (bond->IsDouble())
          for (unsigned int i = 0; i < done.size(); ++i)
            if (done[i] == bond->GetIdx())
              invert = true;
      }

      //cerr << a->GetIdx() << "-" << b->GetIdx() << "-" << c->GetIdx() << "-" << d->GetIdx() << endl;

      double angle;
      if (ab->IsUp() && cd->IsUp())
        angle = 0.0;
      if (ab->IsUp() && cd->IsDown())
        angle = M_PI;
      if (ab->IsDown() && cd->IsUp())
        angle = M_PI;
      if (ab->IsDown() && cd->IsDown())
        angle = 0.0;

      if  (invert)
        mol.SetTorsion(a, b, c, d, angle + M_PI);
      else
        mol.SetTorsion(a, b, c, d, angle);

      done.push_back(bc->GetIdx());
    }

  }

  void OBBuilder::CorrectStereoAtoms(OBMol &mol)
  {
    FOR_ATOMS_OF_MOL (center, mol) {
      if (center->HasData(OBGenericDataType::ChiralData)) {
        OBChiralData *cd = (OBChiralData*) center->GetData(OBGenericDataType::ChiralData);
        vector<unsigned int> refs = cd->GetAtom4Refs(input);
		if (refs.size() < 4)
          continue;
        // We look along the refs[0]-center bond.
        //
        //  1                     1
        //   \        eye        /
        //    0--2    -->   0---C-<2
        //   /                   \
        //  3                     3
        //
        //  fig 1             fig2
        OBAtom *a = mol.GetAtom(refs[0]);
        OBAtom *b = mol.GetAtom(refs[1]);
        OBAtom *c = mol.GetAtom(refs[2]);
        OBAtom *d = mol.GetAtom(refs[3]);

		double angbc = mol.GetTorsion(b, a, &*center, c); // angle 1-0-2 in fig 1
        double angbd = mol.GetTorsion(b, a, &*center, d); // angle 1-0-3 in fig 1

        // this should not never happen if the molecule was build using OBBuilder...
        if ((angbc < 0.0) == (angbd < 0.0))
          continue;

        // invert chirality if needed
        if (angbc > 0.0) {
          if (center->IsAntiClockwise()) {
            Swap(mol, center->GetIdx(), a->GetIdx(), center->GetIdx(), b->GetIdx());
          }
        } else  {
          if (center->IsClockwise()) {
            Swap(mol, center->GetIdx(), a->GetIdx(), center->GetIdx(), b->GetIdx());
          }
        }

      } // HasData
    } // FOR_ATOMS_OF_MOL

  }

} // end namespace OpenBabel


//! \file builder.cpp
//! \brief Handle OBBuilder class

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
#include <openbabel/locale.h>

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
    // TODO: Use OpenDatafile()
    obLocale.SetLocale();

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
    obLocale.RestoreLocale();
  }
 
  vector3 GetCorrectedBondVector(OBAtom *atom1, OBAtom *atom2, int bondOrder = 1)
  {
    double bondLength = 0.0;

    // We create an estimate of the bond length based on the two atoms
    bondLength += etab.CorrectedBondRad(atom1->GetAtomicNum(), atom1->GetHyb());
    bondLength += etab.CorrectedBondRad(atom2->GetAtomicNum(), atom2->GetHyb());

    // These are based on OBBond::GetEquibLength
    if (bondOrder == -1) // aromatic
      bondLength *= 0.93;
    else if (bondOrder == 2)
      bondLength *= 0.91;
    else if (bondOrder == 3)
      bondLength *= 0.87;

    return OBBuilder::GetNewBondVector(atom1, bondLength);
  }

  vector3 OBBuilder::GetNewBondVector(OBAtom *atom)
  {
    return GetNewBondVector(atom, 1.5);
  }

  vector3 OBBuilder::GetNewBondVector(OBAtom *atom, double length) 
  {
    vector3 bond1, bond2, bond3, v1, v2, newbond;
    
    bond1 = VZero;
    bond2 = VZero;
    bond3 = VZero;
    
    if (atom == NULL)
      return VZero;
    
    int dimension = ((OBMol*)atom->GetParent())->GetDimension();

    if (dimension != 2) {
      ////////////
      //   3D   //
      ////////////

      //  
      //  a   --->   a--*
      //
      if (atom->GetValence() == 0) {
        newbond = atom->GetVector() + VX * length;
        return newbond;
      }

      // hyb * = 1
      // ^^^^^^^^^
      //   
      //   (a-1)--a   --->   (a-1)--a--*        angle(a-1, a, *) = 180
      //
      // hyb * = 2
      // ^^^^^^^^^
      // make sure we place the new atom trans to a-2 (if there is an a-2 atom)
      //
      //   (a-2)             (a-2)
      //     \                 \
      //    (a-1)==a   --->   (a-1)==a          angle(a-1, a, *) = 120
      //                              \
      //                               *
      // hyb * = 3
      // ^^^^^^^^^
      // make sure we place the new atom trans to a-2 (if there is an a-2 atom)
      //
      //   (a-2)             (a-2)
      //     \                 \
      //    (a-1)--a   --->   (a-1)--a          angle(a-1, a, *) = 109
      //                              \
      //                               *
      if (atom->GetValence() == 1) {
        bool isCarboxylateO = atom->IsCarboxylOxygen();

        FOR_NBORS_OF_ATOM (nbr, atom) {
          bond1 = atom->GetVector() - nbr->GetVector();

          if (nbr->GetHyb() == 1) // Fix #2119034 & #2013814
            continue;

          FOR_NBORS_OF_ATOM (nbr2, &*nbr) {
            if (&*nbr2 != atom)
              bond2 = nbr->GetVector() - nbr2->GetVector();
            if (isCarboxylateO && nbr2->IsOxygen())
              break; // make sure that the hydrogen is trans to the C=O
          }
        }

        bond1 = bond1.normalize();
        if (bond2 == VZero) {
          // there is no a-2 atom
          v1 = cross(bond1, VY);
          v2 = cross(bond1, v1);
        } else {
          v1 = cross(bond1, bond2);
          v2 = cross(bond1, v1);
          v2 = v2.normalize();
        }      

        if (atom->GetHyb() == 1)
          newbond = bond1;
        else if (atom->GetHyb() == 2)
          newbond = bond1 - v2 * tan(DEG_TO_RAD*60.0);
        else if (atom->GetHyb() == 3)
          newbond = bond1 - v2 * tan(DEG_TO_RAD*70.5);
        
        newbond = newbond.normalize();
        newbond *= length;
        newbond += atom->GetVector();
        return newbond;
      }
    
    //
    //    \         \
    //     X  --->   X--*
    //    /         /
    //
    if (atom->GetValence() == 2) {
      FOR_NBORS_OF_ATOM (nbr, atom) {
        if (bond1 == VZero)
          bond1 = atom->GetVector() - nbr->GetVector();
        else
          bond2 = atom->GetVector() - nbr->GetVector();
      }

      bond1 = bond1.normalize();
      bond2 = bond2.normalize();
      v1 = bond1 + bond2;
      v1 = v1.normalize();
     
      if (atom->GetHyb() == 2)
        newbond = v1;
      if (atom->GetHyb() == 3) {
        v2 = cross(bond1, bond2);
        //v1 = v1.normalize();
        newbond = v2 + v1 * tan(DEG_TO_RAD*35.25);
      }
      
      newbond = newbond.normalize();
      newbond *= length;
      newbond += atom->GetVector();
      return newbond;
    }
    
    
    //
    //    \          \
    //   --X  --->  --X--*
    //    /          /
    //
    if (atom->GetValence() == 3) {
      FOR_NBORS_OF_ATOM (nbr, atom) {
        if (bond1 == VZero)
          bond1 = atom->GetVector() - nbr->GetVector();
        else if (bond2 == VZero)
          bond2 = atom->GetVector() - nbr->GetVector();
        else
          bond3 = atom->GetVector() - nbr->GetVector();
      }
          
      bond1 = bond1.normalize();
      bond2 = bond2.normalize();
      bond3 = bond3.normalize();
      newbond = bond1 + bond2 + bond3;
      newbond = newbond.normalize();
      newbond *= length;
      newbond += atom->GetVector();
      return newbond;
    }

  } else {
    ////////////
    //   2D   //
    ////////////
    OBBondIterator i;
      
    //
    //  a   --->   a---*
    //
    if (atom->GetValence() == 0) {
      newbond = atom->GetVector() + VX * length;
      // Check that the vector is still finite before returning
      if (!isfinite(newbond.x()) || !isfinite(newbond.y()))
        newbond.Set(0.0, 0.0, 0.0);
      return newbond;
    }

    // hyb * = 1                                                                //
    // ^^^^^^^^^                                                                //
    //                                                                          //
    //   (a-1)--a   --->   (a-1)--a--*        angle(a-1, a, *) = 180            //
    //                                                                          //
    // hyb * = 2                                                                //
    // ^^^^^^^^^                                                                //
    // make sure we place the new atom trans to a-2 (if there is an a-2 atom)   //
    //                                                                          //
    //   (a-2)             (a-2)                                                //
    //     \                 \                                                  //
    //    (a-1)==a   --->   (a-1)==a          angle(a-1, a, *) = 120            //
    //                              \                                           //
    //                               *                                          //
    // hyb * = 3                                                                //
    // ^^^^^^^^^                                                                //
    // make sure we place the new atom trans to a-2 (if there is an a-2 atom)   //
    //                                                                          //
    //   (a-2)             (a-2)                                                //
    //     \                 \                                                  //
    //    (a-1)--a   --->   (a-1)--a          angle(a-1, a, *) = 109            //
    //                              \                                           //
    //                               *                                          //
    if (atom->GetValence() == 1) {
      OBAtom *nbr = atom->BeginNbrAtom(i);
      if (!nbr)
        return VZero;
      bond1 = atom->GetVector() - nbr->GetVector(); // bond (a-1)--a

      for (OBAtom *nbr2 = nbr->BeginNbrAtom(i); nbr2; nbr2 = nbr->NextNbrAtom(i))
        if (nbr2 != atom)
          bond2 = nbr->GetVector() - nbr2->GetVector(); // bond (a-2)--(a-1)

      int hyb = atom->GetHyb();
      if (hyb == 1)
        newbond = bond1;
      else if (hyb == 2 || hyb == 3) {
        matrix3x3 m;
        m.RotAboutAxisByAngle(VZ, 60.0);
        newbond = m*bond1;
      }
      newbond.normalize();
      newbond *= length;
      newbond += atom->GetVector();
      return newbond;
    } // GetValence() == 1
  
      //                          //
      //    \         \           //
      //     X  --->   X--*       //
      //    /         /           //
      //                          //
    if (atom->GetValence() == 2) {
      for (OBAtom *nbr = atom->BeginNbrAtom(i); nbr; nbr = atom->NextNbrAtom(i)) {
        if (bond1 == VZero)
          bond1 = atom->GetVector() - nbr->GetVector();
        else
          bond2 = atom->GetVector() - nbr->GetVector();
      }
      bond1.normalize();
      bond2.normalize();
      newbond = bond1 + bond2;
      newbond.normalize();
      newbond *= length;
      newbond += atom->GetVector();
      return newbond;
    }

    //                          //
    //    \          \          //
    //   --X  --->  --X--*      //
    //    /          /          //
    //                          //
    if (atom->GetValence() == 3) {
      if (atom->IsChiral()) {
        OBBond *hash = 0;
        OBBond *wedge = 0;
        vector<OBBond*> plane;
        for (OBAtom *nbr = atom->BeginNbrAtom(i); nbr; nbr = atom->NextNbrAtom(i)) {
          OBBond *bond = atom->GetBond(nbr);

          if (bond->IsWedge()) {
            if (atom == bond->GetBeginAtom())
              wedge = bond;
            else
              hash = bond;
          } else 
            if (bond->IsHash()) {
              if (atom == bond->GetBeginAtom())
                hash = bond;
              else
                wedge = bond;
            } else
              plane.push_back(bond);
        }

        if (wedge && !plane.empty()) {
          bond2 = atom->GetVector() - wedge->GetNbrAtom(atom)->GetVector();
          bond3 = atom->GetVector() - plane[0]->GetNbrAtom(atom)->GetVector();
        } else if (hash && !plane.empty()) {
          bond2 = atom->GetVector() - hash->GetNbrAtom(atom)->GetVector();
          bond3 = atom->GetVector() - plane[0]->GetNbrAtom(atom)->GetVector();
        } else if (plane.size() >= 2) {
          bond2 = atom->GetVector() - plane[0]->GetNbrAtom(atom)->GetVector();
          bond3 = atom->GetVector() - plane[1]->GetNbrAtom(atom)->GetVector();
        } else if (hash && wedge) {
          bond2 = atom->GetVector() - wedge->GetNbrAtom(atom)->GetVector();
          bond3 = atom->GetVector() - hash->GetNbrAtom(atom)->GetVector();
        }
      } else {
        for (OBAtom *nbr = atom->BeginNbrAtom(i); nbr; nbr = atom->NextNbrAtom(i)) {
          if (bond1 == VZero)
            bond1 = atom->GetVector() - nbr->GetVector();
          else if (bond2 == VZero)
            bond2 = atom->GetVector() - nbr->GetVector();
          else
            bond3 = atom->GetVector() - nbr->GetVector();
        }
      }

      bond2.normalize();
      bond3.normalize();
      newbond = -(bond2 + bond3);
      newbond.normalize();
      newbond *= length;
      newbond += atom->GetVector();
      return newbond;
    }
  
  }

  return VZero; //previously undefined
}
  
  // The OBMol mol contains both the molecule to which we want to connect the 
  // fragment and the fragment itself. The fragment containing b will be 
  // rotated and translated. Atom a is the atom from 
  // the main molecule to which we want to connect atom b.
  // NOTE: newpos now uses CorrectedBondVector, so we don't do that below
  bool OBBuilder::Connect(OBMol &mol, int idxA, int idxB, 
                          vector3 &newpos, int bondOrder)
  {
    OBAtom *a = mol.GetAtom(idxA);
    OBAtom *b = mol.GetAtom(idxB);

    if (a == NULL)
      return false;
    if (b == NULL)
      return false;

    OBBitVec fragment = GetFragment(b);
    if (fragment == GetFragment(a))
      return false; // a and b are in the same fragment

    vector3 posa = a->GetVector();
    vector3 posb = b->GetVector();
    // 
    // translate fragment so that atom b is at the origin
    //
    for (unsigned int i = 1; i <= mol.NumAtoms(); ++i) {
      if (fragment.BitIsSet(i)) {
        // the atom is part of the fragment, translate it
        vector3 tmpvec = mol.GetAtom(i)->GetVector();
        tmpvec -= posb;
        mol.GetAtom(i)->SetVector(tmpvec);
      }
    }
    // 
    // rotate the fragment to align the bond directions Mol-a-b and a-b-fragment
    //
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
    for (unsigned int i = 1; i <= mol.NumAtoms(); ++i) {
      if (fragment.BitIsSet(i)) {
        vector3 tmpvec = mol.GetAtom(i)->GetVector();
        tmpvec *= xymat; //apply the rotation
        mol.GetAtom(i)->SetVector(tmpvec);
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
    for (unsigned int i = 1; i <= mol.NumAtoms(); ++i) {
      if (fragment.BitIsSet(i)) {
        vector3 tmpvec = mol.GetAtom(i)->GetVector();
        tmpvec *= xzmat; //apply the rotation
        mol.GetAtom(i)->SetVector(tmpvec);
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
    for (unsigned int i = 1; i <= mol.NumAtoms(); ++i) {
      if (fragment.BitIsSet(i)) {
        vector3 tmpvec = mol.GetAtom(i)->GetVector();
        tmpvec *= yzmat; //apply the rotation
        mol.GetAtom(i)->SetVector(tmpvec);
      }
    }
    // 
    // translate fragment 
    //
    for (unsigned int i = 1; i <= mol.NumAtoms(); ++i) {
      if (fragment.BitIsSet(i)) {
        // translate the fragment
        vector3 tmpvec = mol.GetAtom(i)->GetVector();
        tmpvec += newpos;
        mol.GetAtom(i)->SetVector(tmpvec);
      }
    }     
    //
    // Create the bond between the two fragments
    //
    OBBond *bond = mol.NewBond();
    bond->SetBegin(a);
    bond->SetEnd(b);
    bond->SetBondOrder(bondOrder);
    a->AddBond(bond);
    b->AddBond(bond);
  
    return true;
  }

  bool OBBuilder::Connect(OBMol &mol, int idxA, int idxB, int bondOrder)
  {
    vector3 newpos = GetCorrectedBondVector(mol.GetAtom(idxA), mol.GetAtom(idxB), bondOrder);
    return Connect(mol, idxA, idxB, newpos, bondOrder);
  }


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

  bool OBBuilder::Swap(OBMol &mol, int idxA, int idxB, int idxC, int idxD)
  {
    OBAtom *a = mol.GetAtom(idxA);
    OBAtom *b = mol.GetAtom(idxB);
    OBAtom *c = mol.GetAtom(idxC);
    OBAtom *d = mol.GetAtom(idxD);
    
    // make sure the atoms exist
    if (a == NULL || b == NULL || c == NULL || d == NULL)
      return false;

    OBBond *bond1 = mol.GetBond(idxA, idxB);
    OBBond *bond2 = mol.GetBond(idxC, idxD);

    // save the original bond orders
    int bondOrder1 = bond1->GetBondOrder();
    int bondOrder2 = bond2->GetBondOrder();

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
    if (!Connect(mol, idxA, idxD, posB, bondOrder2))
      return false;
    if (!Connect(mol, idxC, idxB, posD, bondOrder2))
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

  OBBitVec OBBuilder::GetFragment(OBAtom *atom)
  { 
    OBBitVec fragment;

    fragment.SetBitOn(atom->GetIdx());
    AddNbrs(fragment, atom);
    
    return fragment; 
  }

  // Handle the case of 2D molecules
  // Currently just a helper function since we can't add API in 2.2
  void Convert2DCoords(OBMol *obMolecule)
  {
    if (obMolecule == NULL || obMolecule->GetDimension() != 2)
      return;

    double sum = 0.0;
    FOR_BONDS_OF_MOL (bond, obMolecule) {
      sum += bond->GetLength();
    }
    double scale = (1.5 * obMolecule->NumBonds()) / sum;
    FOR_ATOMS_OF_MOL (atom, obMolecule) {
      vector3 vec = atom->GetVector();
      vec.SetX(vec.x() * scale);
      vec.SetY(vec.y() * scale);
      atom->SetVector(vec);
    }
    obMolecule->Center();
    
    // Check for pairs of atoms on top of each other
    // e.g. "OH" or "COOH" labels in 2D files.
    // First check bonding pairs
    FOR_BONDS_OF_MOL(bond, obMolecule)
      {
        if (bond->GetLength() < 1.0e-6) {
          bond->SetLength(bond->GetEquibLength());
        }
      }
    // Now check non-bonded pairs
    FOR_PAIRS_OF_MOL(p, obMolecule)
      {
        OBAtom *a = obMolecule->GetAtom((*p)[0]);
        OBAtom *b = obMolecule->GetAtom((*p)[1]);
        if (fabs(a->GetDistance(b)) < 1.0e-6) {
          vector3 v1;
          v1.randomUnitVector();
          a->SetVector(a->GetVector() + v1);
          b->SetVector(b->GetVector() - v1);
        }
      }
    
    // place end atoms of wedge bonds at +1.0 Z
    // place end atoms of hash bonds at -1.0 Z
    FOR_ATOMS_OF_MOL (atom, obMolecule) {
      FOR_BONDS_OF_ATOM (bond, &*atom) {
        if (bond->IsHash() && (&*atom == bond->GetBeginAtom())) {
          vector3 vec = bond->GetEndAtom()->GetVector();
          vec.SetZ(-1.0);
          bond->GetEndAtom()->SetVector(vec);
        } else if (bond->IsWedge() && (&*atom == bond->GetBeginAtom())) {
          vector3 vec = bond->GetEndAtom()->GetVector();
          vec.SetZ(1.0);
          bond->GetEndAtom()->SetVector(vec);
        }
      }
    }
    // If someone wants actual coordinates from OBBuilder, they should
    // call forcefields themselves
//     OBForceField *ff = OBForceField::FindForceField("UFF");
//     if (ff) {
//       ff->Setup(*obMolecule);
//       ff->ConjugateGradients(250);
//       ff->GetCoordinates(*obMolecule);
//     }
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
    OBMol workMol = mol;
    
    if (workMol.GetDimension() == 2) {
      Convert2DCoords(&workMol);
      mol = workMol;
      mol.SetDimension(3);
      return true;
    }

    // How many ring atoms are there?
    int ratoms = 0;
    FOR_ATOMS_OF_MOL(a, mol)
      if (a->IsInRing())
        ratoms++;

    // delete all bonds in the working molecule
    // we will add them back at the end
    while (workMol.NumBonds())
      workMol.DeleteBond(workMol.GetBond(0));
    
    workMol.SetHybridizationPerceived();
    
    if (ratoms) {
      //datafile is read only on first use of Build()
      if(_fragments.empty())
        LoadFragments();

      // Skip all fragments that are too big to match
      // Note: It would be better to compare to the size of the largest
      //       isolated ring system instead of comparing to ratoms
      for (i = _fragments.begin();i != _fragments.end() && i->first->NumAtoms() > ratoms;++i);
      
      // Loop through  the remaining fragments and assign the coordinates from
      // the first (most complex) fragment.
      // Stop if there are no unassigned ring atoms (ratoms).
      for (;i != _fragments.end() && ratoms;++i) {

        if (i->first != NULL && i->first->Match(mol)) { 
          mlist = i->first->GetUMapList();
            
          for (j = mlist.begin();j != mlist.end();++j) { // for all matches

            // Has any atom of this match already been added?
            bool alreadydone = false;
            for (k = j->begin(); k != j->end(); ++k) { // for all atoms of the fragment
              if (vfrag.BitIsSet(*k)) {
                alreadydone = true;
                break;
              }
            }
            if (alreadydone) // the found match is already added
              continue;

            int index, index2, counter = 0;
            for (k = j->begin(); k != j->end(); ++k) { // for all atoms of the fragment
              index = *k;     
              vfrag.SetBitOn(index); // set vfrag for all atoms of fragment
              if (mol.GetAtom(index)->IsInRing())
                ratoms--;
              
              // set coordinates for atoms
              OBAtom *atom = workMol.GetAtom(index);
              atom->SetVector(i->second[counter]);
              counter++;
            }

            // add the bonds for the fragment
            for (k = j->begin(); k != j->end(); ++k) {
              index = *k;
              OBAtom *atom1 = mol.GetAtom(index);
                
              for (k2 = j->begin(); k2 != j->end(); ++k2) {
                index2 = *k2;
                OBAtom *atom2 = mol.GetAtom(index2);
                OBBond *bond = atom1->GetBond(atom2);

                if (bond != NULL) 
                  workMol.AddBond(*bond);
              }
            }
          }
        }
      }
    } // if (ratoms)

    // iterate over all atoms to place them in 3D space
    FOR_DFS_OF_MOL (a, mol) {
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
          Connect(workMol, prev->GetIdx(), a->GetIdx(), mol.GetBond(prev, &*a)->GetBondOrder());
          // set the correct bond order
          int bondOrder = mol.GetBond(prev->GetIdx(), a->GetIdx())->GetBondOrder();
          workMol.GetBond(prev->GetIdx(), a->GetIdx())->SetBondOrder(bondOrder);
        }
        
        OBBitVec fragment = GetFragment(workMol.GetAtom(a->GetIdx()));
        vdone |= fragment; // mark this fragment as done

        continue;
      }

      //
      // below is the code to add non-fragment atoms
      //
   
      // get the position for the new atom, this is done with GetNewBondVector
      if (prev != NULL) {
        int bondType = a->GetBond(prev)->GetBO();
        if (a->GetBond(prev)->IsAromatic())
          bondType = -1;

        molvec = GetCorrectedBondVector(workMol.GetAtom(prev->GetIdx()),
                                        workMol.GetAtom(a->GetIdx()),
                                        bondType);
        moldir = molvec - workMol.GetAtom(prev->GetIdx())->GetVector();
      } else {
        molvec = VX;
        moldir = VX;
      }
      
      vdone.SetBitOn(a->GetIdx());
      
      // place the atom the atom 
      OBAtom *atom = workMol.GetAtom(a->GetIdx());
      atom->SetVector(molvec);

      // add bond between previous part and added atom
      if (prev != NULL) {
        OBBond *bond = a->GetBond(prev); // from mol
        workMol.AddBond(*bond);
      }

    }

    // Ensure all bonds from the old molecule exist in the new molecule
    int beginIdx, endIdx;
    FOR_BONDS_OF_MOL(b, mol) {
      beginIdx = b->GetBeginAtomIdx();
      endIdx = b->GetEndAtomIdx();
      if (!workMol.GetBond(beginIdx, endIdx)) {
        // We need to duplicate the old bond
        workMol.AddBond(beginIdx, endIdx, b->GetBO(), b->GetFlags());
      }
    }
    
    // correct the chirality
    CorrectStereoBonds(workMol);
    CorrectStereoAtoms(workMol);

    mol = workMol;
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

/**********************************************************************
builder.cpp - Class to create structures.

Copyright (C) 2007-2008 by Tim Vandermeersch
                           <tim.vandermeersch@gmail.com>

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
#include <assert.h>


#include <openbabel/stereo/stereo.h>
#include <openbabel/stereo/cistrans.h>
#include <openbabel/stereo/tetrahedral.h>
#include <openbabel/stereo/squareplanar.h>

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

    if (bondLength < 1.0)
      bondLength = 1.0;

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
    vector3 bond1, bond2, bond3, bond4, bond5, v1, v2, v3, newbond;

    bond1 = VZero;
    bond2 = VZero;
    bond3 = VZero;
    bond4 = VZero;
    bond5 = VZero;

    if (atom == NULL)
      return VZero;

    int dimension = ((OBMol*)atom->GetParent())->GetDimension();

    if (dimension != 2) {
      /////////////////
      //   3D or 0D  //
      /////////////////

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
        v1 = cross(bond1, bond2);
        if (bond2 == VZero || v1 == VZero) {
          vector3 vrand;
          vrand.randomUnitVector();
          double angle = fabs(acos(dot(bond1, vrand)) * RAD_TO_DEG);
          while (angle < 45.0 || angle > 135.0) {
            vrand.randomUnitVector();
            angle = fabs(acos(dot(bond1, vrand)) * RAD_TO_DEG);
          }
          // there is no a-2 atom
          v1 = cross(bond1, vrand); // so find a perpendicular, given the random vector (this doesn't matter here)
          v2 = cross(bond1, v1);
        } else {
          v1 = cross(bond1, bond2);
          v2 = cross(bond1, v1);
        }
        v2 = v2.normalize();

        // check to see if atom is a square planar in disguise
        if (atom->GetHyb() == 3) {
          OBStereoFacade stereoFacade((OBMol*)atom->GetParent(), false); // Don't reperceive
          if (stereoFacade.HasSquarePlanarStereo(atom->GetId()))
            atom->SetHyb(4); // force sq. planar geometry for sq. planar stereo
        }

        if (atom->GetHyb() == 1)
          newbond = bond1; // i.e., in the exact opposite direction
        else if (atom->GetHyb() == 2)
          newbond = bond1 - v2 * tan(DEG_TO_RAD*60.0);
        else if (atom->GetHyb() == 3)
          newbond = bond1 - v2 * tan(DEG_TO_RAD*(180.0 - 109.471));
        else if (atom->GetHyb() == 4)
          newbond = bond1; // like 5-coordinate below, we want a 180-degree bond (trans)
        else if (atom->GetHyb() == 5) {
          /* the first two atoms are the axial ones;  the third, fourth, and fifth atom are equatorial */
          newbond = bond1;
        } else if (atom->GetHyb() == 6) {
          newbond = bond1 - v2 * tan(DEG_TO_RAD*90.0);
        }

        newbond = newbond.normalize();
        newbond *= length;
        newbond += atom->GetVector();
        return newbond;
      }

      //
      //    \	      \
      //     X  --->   X--*
      //    /         /
      //
      if (atom->GetValence() == 2) {
        FOR_NBORS_OF_ATOM (nbr, atom) {
          if (bond1 == VZero)
            bond1 = atom->GetVector() - nbr->GetVector();
          else if (bond2 == VZero)
            bond2 = atom->GetVector() - nbr->GetVector();
        }

        bond1 = bond1.normalize();
        bond2 = bond2.normalize();
        v1 = bond1 + bond2;
        v1 = v1.normalize();

        if (atom->GetHyb() == 2)
          newbond = v1;
        if (atom->GetHyb() == 3) {
          v2 = cross(bond1, bond2); // find the perpendicular
          v2.normalize();
          newbond = bond1 - v2 * tan(DEG_TO_RAD*(180.0 - 109.471));
          newbond = v2 + v1 * (sqrt(2.0) / 2.0); // used to be tan(70.53 degrees/2) which is sqrt(2.0) / 2.0
        }
        if (atom->GetHyb() == 5 || atom->GetHyb() == 4) {
          /* add the first equatorial atom, orthogonally to bond1 (and bond2 = -bond1) */
          /* is atom order correct?  I don't think it matters, but I might have to ask a chemist
           * whether PClF4 would be more likely to have an equatorial or axial Cl-P bond */
          vector3 vrand;
          vrand.randomUnitVector();
          double angle = fabs(acos(dot(bond1, vrand)) * RAD_TO_DEG);
          while (angle < 45.0 || angle > 135.0) {
            vrand.randomUnitVector();
            angle = fabs(acos(dot(bond1, vrand)) * RAD_TO_DEG);
          }
          v1 = cross(bond1, vrand);
          v1 = v1.normalize();
          newbond = v1;
        }
        if (atom->GetHyb() == 6) {
          v2 = cross(bond1, bond2);
          newbond = v2;
        }

        newbond = newbond.normalize();
        newbond *= length;
        newbond += atom->GetVector();
        return newbond;
      }


      /* UFF:
       *    b lg dg  o  y
       *  b - 45 30 45 30
       * lg    - 45  0 45
       * dg       - 45 30
       *  o          - 45
       *  y             -

       * 94s:
       *    b lg dg  o  y
       *  b - 34 34 34 34
       * lg    - 48 21 48
       * dg       - 48 21
       *  o          - 48
       *  y             -

       //
       //    \	       \
       //   --X  --->  --X--*
       //    /          /
       //
       */
      if (atom->GetValence() == 3) {
        if (atom->GetHyb() == 3) {
          FOR_NBORS_OF_ATOM (nbr, atom) {
            if (bond1 == VZero)
              bond1 = atom->GetVector() - nbr->GetVector();
            else if (bond2 == VZero)
              bond2 = atom->GetVector() - nbr->GetVector();
            else if (bond3 == VZero)
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

        if (atom->GetHyb() == 4) { // OK, we want this at -bond3, since bond1 & bond2 are opposite
          FOR_NBORS_OF_ATOM (nbr, atom) {
            if (bond1 == VZero)
              bond1 = atom->GetVector() - nbr->GetVector();
            else if (bond2 == VZero)
              bond2 = atom->GetVector() - nbr->GetVector();
            else if (bond3 == VZero)
              bond3 = atom->GetVector() - nbr->GetVector();
          }

          bond3 = bond3.normalize();

          newbond = bond3;
          newbond = newbond.normalize();
          newbond *= length;
          newbond += atom->GetVector();
          return newbond;
        }

        if (atom->GetHyb() == 5) {
          FOR_NBORS_OF_ATOM (nbr, atom) {
            if (bond1 == VZero)
              bond1 = atom->GetVector() - nbr->GetVector();
            else if (bond2 == VZero)
              bond2 = atom->GetVector() - nbr->GetVector();
            else if (bond3 == VZero)
              bond3 = atom->GetVector() - nbr->GetVector();
          }

          bond1 = bond1.normalize();
          bond2 = bond2.normalize();
          bond3 = bond3.normalize();

          v1 = cross(bond1, bond3);
          v1 = v1.normalize();

          newbond = v1 + tan(DEG_TO_RAD*30.0) * bond3;
          newbond = newbond.normalize();
          newbond *= length;
          newbond += atom->GetVector();
          return newbond;

        }

        if (atom->GetHyb() == 6) {
          FOR_NBORS_OF_ATOM (nbr, atom) {
            if (bond1 == VZero)
              bond1 = atom->GetVector() - nbr->GetVector();
            else if (bond2 == VZero)
              bond2 = atom->GetVector() - nbr->GetVector();
            else if (bond3 == VZero)
              bond3 = atom->GetVector() - nbr->GetVector();
          }

          /* the way things work, newbond is equal to bond1, but will show up at -bond1 next time around */
          newbond = bond1;
          newbond = newbond.normalize();
          newbond *= length;
          newbond += atom->GetVector();
          return newbond;
        }
      }

      if (atom->GetValence() == 4) {
        if (atom->GetHyb() == 6) {
          FOR_NBORS_OF_ATOM (nbr, atom) {
            if (bond1 == VZero)
              bond1 = atom->GetVector() - nbr->GetVector();
            else if (bond2 == VZero)
              bond2 = atom->GetVector() - nbr->GetVector();
            else if (bond3 == VZero)
              bond3 = atom->GetVector() - nbr->GetVector();
            else if (bond4 == VZero)
              bond4 = atom->GetVector() - nbr->GetVector();
          }

          newbond = bond2;
          newbond = newbond.normalize();
          newbond *= length;
          newbond += atom->GetVector();
          return newbond;
        }

        if (atom->GetHyb() == 5) {
          FOR_NBORS_OF_ATOM (nbr, atom) {
            if (bond1 == VZero)
              bond1 = atom->GetVector() - nbr->GetVector();
            else if (bond2 == VZero)
              bond2 = atom->GetVector() - nbr->GetVector();
            else if (bond3 == VZero)
              bond3 = atom->GetVector() - nbr->GetVector();
            else if (bond4 == VZero)
              bond4 = atom->GetVector() - nbr->GetVector();
          }

          bond1 = bond1.normalize();
          bond2 = bond2.normalize();
          bond3 = bond3.normalize();
          bond4 = bond4.normalize();

          v1 = cross(bond1, bond3);
          v1 = v1.normalize();

          newbond = -v1 + tan(DEG_TO_RAD*30.0) * bond3;
          newbond = newbond.normalize();
          newbond *= length;
          newbond += atom->GetVector();
          return newbond;

        }

      }

      if (atom->GetValence() == 5) {
        if (atom->GetHyb() == 6) {
          FOR_NBORS_OF_ATOM (nbr, atom) {
            if (bond1 == VZero)
              bond1 = atom->GetVector() - nbr->GetVector();
            else if (bond2 == VZero)
              bond2 = atom->GetVector() - nbr->GetVector();
            else if (bond3 == VZero)
              bond3 = atom->GetVector() - nbr->GetVector();
            else if (bond4 == VZero)
              bond4 = atom->GetVector() - nbr->GetVector();
            else if (bond5 == VZero)
              bond5 = atom->GetVector() - nbr->GetVector();
          }

          newbond = bond3;
          newbond = newbond.normalize();
          newbond *= length;
          newbond += atom->GetVector();
          return newbond;
        }
      }

      // Undefined case -- return a random vector of length specified
      newbond.randomUnitVector();
      newbond *= length;
      newbond += atom->GetVector();
      return newbond;
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
        else if (hyb == 2 || hyb == 3 || hyb == 0) {
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
        OBStereoFacade stereoFacade((OBMol*)atom->GetParent());
        if (stereoFacade.HasTetrahedralStereo(atom->GetId())) {
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

      // Undefined case -- return a random vector of length specified
      newbond.randomUnitVector();
      newbond.z() = 0.0;
      newbond.normalize();
      newbond *= length;
      newbond += atom->GetVector();
      return newbond;
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

    bool connectedFrag = true; // normal case
    // If we don't have any neighbors, assume that the fragment
    // is inserted at the end of the molecule and set anything after the atom.
    // This lets us place fragments like Cp rings with dummy atoms
    if (b->GetAtomicNum() == 0) {
      connectedFrag = false;
      fragment.SetRangeOn(b->GetIdx(), mol.NumAtoms());
    }

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
    vector3 crossdir;
    if (!connectedFrag) { // nothing bonded to b, like a Cp ring
      vector3 firstDir, secondDir;
      // Try finding the next atom
      OBAtom *nextAtom = mol.GetAtom(b->GetIdx() + 1);
      if (nextAtom) {
        firstDir = nextAtom->GetVector() - b->GetVector();
        // we'll try finding another atom
        OBAtom *secondAtom = mol.GetAtom(b->GetIdx() + 2);
        if (secondAtom)
          secondDir = secondAtom->GetVector() - b->GetVector();
        else
          secondDir.randomUnitVector(); // pick something at random
        // but not too shallow, or the cross product won't work well
        double angle = fabs(acos(dot(firstDir, secondDir)) * RAD_TO_DEG);
        while (angle < 45.0 || angle > 135.0) {
          secondDir.randomUnitVector();
          angle = fabs(acos(dot(firstDir, secondDir)) * RAD_TO_DEG);
        }
        // Now we find a perpendicular vector to the fragment
        crossdir = cross(firstDir, secondDir);
        fragdir = crossdir;
      }
    }
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
    if (!connectedFrag)
      fragdir = crossdir;
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
    if (!connectedFrag)
      fragdir = crossdir;
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
    // Get a neighbor of a and of b for setting the dihedral later
    //
    OBAtom* nbr_a = NULL;
    FOR_NBORS_OF_ATOM(nbr, a) {
      nbr_a = &*nbr;
      break;
    }
    OBAtom* nbr_b = NULL;
    FOR_NBORS_OF_ATOM(nbr, b) {
      if (fragment.BitIsSet(nbr->GetIdx())) {
        nbr_b = &*nbr;
        break;
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
    //
    // Set the dihedral between the two fragments (for the moment, only handle double bonds)
    //
    // For example, if a double bond is coming off a ring, then the dihedral
    // should be 180, e.g. for I/C=C\1/NC1 (don't worry about whether cis or trans
    // at this point - this will be corrected later)
    //
    if(bondOrder==2 && a->GetHyb()==2 && b->GetHyb()==2 && nbr_a && nbr_b)
      mol.SetTorsion(nbr_a, a, b, nbr_b, 180 * DEG_TO_RAD);

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

  // Variation of OBBuilder::Swap that allows swapping with a vector3 rather than
  // an explicit bond. This is useful for correcting stereochemistry at Tet Centers
  // where it is sometimes necessary to swap an existing bond with the location
  // of an implicit hydrogen (or lone pair), in order to correct the stereo.
  bool OBBuilder::SwapWithVector(OBMol &mol, int idxA, int idxB, int idxC, const vector3 &newlocation)
  {
    OBAtom *a = mol.GetAtom(idxA);
    OBAtom *b = mol.GetAtom(idxB);
    OBAtom *c = mol.GetAtom(idxC);

    // make sure the atoms exist
    if (a == NULL || b == NULL || c == NULL)
      return false;

    OBBond *bond1 = mol.GetBond(idxA, idxB);

    // make sure a-b is connected
    if (bond1 == NULL)
      return false;

    // make sure the bond are not in a ring
    if (bond1->IsInRing())
      return false;

    // save the original bond order
    int bondOrder1 = bond1->GetBondOrder();

    // delete the bond
    mol.DeleteBond(bond1);

    // Get the old bond vector
    vector3 bondB = b->GetVector() - a->GetVector();
    vector3 bondD =  newlocation   - c->GetVector();

    // Get the new positions for B and D
    vector3 newB = c->GetVector() + bondB.length() * (bondD/bondD.length());
    vector3 newD = a->GetVector() + bondD.length() * (bondB/bondB.length());

    // connect the fragments
    if (!Connect(mol, idxC, idxB, newB, bondOrder1))
      return false;

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

    OBBond *bond1 = mol.GetBond(idxA, idxB);
    OBBond *bond2 = mol.GetBond(idxC, idxD);

    // make sure a-b and c-d are connected
    if (bond1 == NULL || bond2 == NULL)
      return false;

    // make sure the bonds are not in a ring
    if (bond1->IsInRing() || bond2->IsInRing())
      return false;

    // save the original bond orders
    int bondOrder1 = bond1->GetBondOrder();
    int bondOrder2 = bond2->GetBondOrder();

    // delete the bonds
    mol.DeleteBond(bond1);
    mol.DeleteBond(bond2);

    // Get the bond vectors
    vector3 bondB = b->GetVector() - a->GetVector();
    vector3 bondD = d->GetVector() - c->GetVector();

    // Get the new positions for B and D
    vector3 newB = c->GetVector() + bondB.length() * (bondD/bondD.length());
    vector3 newD = a->GetVector() + bondD.length() * (bondB/bondB.length());

    // connect the fragments
    if (!Connect(mol, idxA, idxD, newD, bondOrder2))
      return false;
    if (!Connect(mol, idxC, idxB, newB, bondOrder1))
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

    // Treat a 2D structure like a 0D one
    if (workMol.GetDimension() == 2)
      workMol.SetDimension(0);

    // Count the number of ring atoms.
    unsigned int ratoms = 0;
    FOR_ATOMS_OF_MOL(a, mol)
      if (a->IsInRing()) {
        ratoms++;
        if (_keeprings) // Mark these as fragments
          vfrag.SetBitOn(a->GetIdx());
      }

    if (_keeprings) {
      // Delete all non-ring bonds
      std::vector<OBBond*> for_deletion;
      FOR_BONDS_OF_MOL(b, workMol)
        if (!b->IsInRing())
          for_deletion.push_back(&(*b));
      for(std::vector<OBBond*>::iterator it=for_deletion.begin(); it!=for_deletion.end(); ++it) {
        workMol.DeleteBond(*it);
      }
    }
    else
      // Delete all bonds in the working molecule
      // (we will add them back at the end)
      while (workMol.NumBonds())
        workMol.DeleteBond(workMol.GetBond(0));

    workMol.SetHybridizationPerceived();

    if (ratoms && !_keeprings) {
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

            // Have any atoms of this match already been added?
            int alreadydone = 0;
            int match_idx;
            for (k = j->begin(); k != j->end(); ++k) // for all atoms of the fragment
              if (vfrag.BitIsSet(*k)) {
                alreadydone += 1;
                if (alreadydone > 1) break;
                match_idx = *k;
              }
            bool spiro = false;
            if (alreadydone > 1) // At least two atoms of the found match have already been added
              continue;
            else if (alreadydone == 1) {
              spiro = IsSpiroAtom(mol.GetAtom(match_idx)->GetId(), mol);
              if (!spiro) continue;
            }

            ratoms += alreadydone - static_cast<int> (j->size()); // Update the number of available ring atoms
            for (k = j->begin(); k != j->end(); ++k)
              vfrag.SetBitOn(*k); // Set vfrag for all atoms of fragment

            if (spiro) {
              vector<int> match_indices(1); // In future, it could be of longer length
              match_indices[0] = match_idx;
              ConnectFrags(mol, workMol, *j, i->second, match_indices);
            }
            else {
              int counter;
              for (k = j->begin(), counter=0; k != j->end(); ++k, ++counter) { // for all atoms of the fragment
                // set coordinates for atoms
                OBAtom *atom = workMol.GetAtom(*k);
                atom->SetVector(i->second[counter]);
              }
            }

            // add the bonds for the fragment
            int index2;
            for (k = j->begin(); k != j->end(); ++k) {
              OBAtom *atom1 = mol.GetAtom(*k);

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

      if (vfrag.BitIsSet(a->GetIdx())) { // Is this atom part of a fragment?
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
        // We don't want to plant all base atoms at exactly the same spot.
        // (or in exactly the same direction)
        // So we'll add a slight tweak -- fixes problem reported by Kasper Thofte
        vector3 randomOffset;
        randomOffset.randomUnitVector();
        molvec = VX + 0.1 * randomOffset;
        moldir = VX + 0.01 * randomOffset;
      }

      vdone.SetBitOn(a->GetIdx());

      // place the atom
      OBAtom *atom = workMol.GetAtom(a->GetIdx());
      atom->SetVector(molvec);

      // add bond between previous part and added atom
      if (prev != NULL) {
        OBBond *bond = a->GetBond(prev); // from mol
        workMol.AddBond(*bond);
      }

    }

    // Make sure we keep the bond indexes the same
    // so we'll delete the bonds again and copy them
    // Fixes PR#3448379 (and likely other topology issues)
    while (workMol.NumBonds())
      workMol.DeleteBond(workMol.GetBond(0));

    int beginIdx, endIdx;
    FOR_BONDS_OF_MOL(b, mol) {
      beginIdx = b->GetBeginAtomIdx();
      endIdx = b->GetEndAtomIdx();
      workMol.AddBond(beginIdx, endIdx, b->GetBO(), b->GetFlags());
    }

    // correct the chirality
    CorrectStereoBonds(workMol);
    CorrectStereoAtoms(workMol);
    workMol.SetChiralityPerceived();
    mol = workMol;
    mol.SetDimension(3);

    return true;
  }

  void OBBuilder::ConnectFrags(OBMol &mol, OBMol &workMol, vector<int> match, vector<vector3> coords,
                               vector<int> pivot)
  {
    if (pivot.size() != 1) // Only handle spiro at the moment
      return;

    OBAtom* p = workMol.GetAtom(pivot[0]);
    OBBitVec frag = GetFragment(p); // The existing fragment
    vector3 posp = p->GetVector();

    // Set coords of new fragment to place the pivot at the origin
    vector3 posp_new;
    vector<int>::iterator match_it;
    int counter;
    for (match_it=match.begin(), counter=0; match_it!=match.end(); ++match_it, ++counter)
      if (*match_it == pivot[0]) {
        posp_new = coords[counter];
        break;
      }
    counter = 0;
    for (match_it=match.begin(), counter=0; match_it!=match.end(); ++match_it, ++counter)
      workMol.GetAtom(*match_it)->SetVector( coords[counter] - posp_new );

    // Find vector that bisects existing angles at the pivot in each fragment
    // and align them
    //                                        \   /
    //  \        \   /    bisect  \             P    align   \                /
    //   P  and    P       --->    P--v1  and   |    --->     P--v1  and v2--P
    //  /                         /             v2           /                \

    // Get v1 (from the existing fragment)
    vector3 bond1, bond2, bond3, bond4, v1;
    bond1 = VZero;
    OBAtom atom1, atom2;
    FOR_NBORS_OF_ATOM(nbr, p) {
      if (bond1 == VZero) {
        atom1 = *nbr;
        bond1 = posp - atom1.GetVector();
      }
      else {
        atom2 = *nbr;
        bond2 = posp - atom2.GetVector();
      }
    }
    bond1 = bond1.normalize();
    bond2 = bond2.normalize();
    v1 = bond1 + bond2;
    v1 = v1.normalize();

    // Get v2 (from the new fragment)
    vector3 v2;
    vector<int> nbrs;
    vector<vector3> nbr_pos;
    FOR_NBORS_OF_ATOM(nbr, mol.GetAtom(pivot[0]))
      if (nbr->GetIdx() != atom1.GetIdx() && nbr->GetIdx() != atom2.GetIdx()) {
        nbrs.push_back(nbr->GetIdx());
        nbr_pos.push_back(workMol.GetAtom(nbr->GetIdx())->GetVector());
      }
    //assert(nbrs.size()==2);
    bond3 = nbr_pos[0] - VZero; // The pivot is at the origin, hence VZero
    bond4 = nbr_pos[1] - VZero;
    bond3 = bond3.normalize();
    bond4 = bond4.normalize();
    v2 = bond3 + bond4;
    v2 = v2.normalize();

    // Set up matrix to rotate around v1 x v2 by the angle between them
    double ang = vectorAngle(v1, v2);
    vector3 cp = cross(v1, v2);
    matrix3x3 mat;
    mat.RotAboutAxisByAngle(cp, ang);

    // Apply rotation
    vector3 tmpvec;
    for (match_it=match.begin(); match_it!=match.end(); ++match_it) {
      tmpvec = workMol.GetAtom(*match_it)->GetVector();
      tmpvec *= mat;
      workMol.GetAtom(*match_it)->SetVector( tmpvec );
    }

    // Rotate the new fragment 90 degrees to make a tetrahedron
    tmpvec = cross(bond1, bond2); // The normal to the ring
    v1 = cross(tmpvec, v1); // In the plane of the ring, orthogonal to tmpvec and the original v1
    v2 = cross(bond3, bond4); // The normal to ring2 - we want to align v2 to v1
    ang = vectorAngle(v1, v2); // Should be 90
    cp = cross(v1, v2);
    mat.RotAboutAxisByAngle(cp, ang);
    for (match_it=match.begin(); match_it!=match.end(); ++match_it) {
      tmpvec = workMol.GetAtom(*match_it)->GetVector();
      tmpvec *= mat;
      workMol.GetAtom(*match_it)->SetVector( tmpvec );
    }

    // Translate to existing pivot location
    for (match_it=match.begin(); match_it!=match.end(); ++match_it)
      workMol.GetAtom(*match_it)->SetVector( workMol.GetAtom(*match_it)->GetVector() + posp );

    // Create the bonds between the two fragments
    for (vector<int>::iterator nbr_id=nbrs.begin(); nbr_id!=nbrs.end(); ++nbr_id)
      workMol.AddBond(p->GetIdx(), *nbr_id, 1, mol.GetBond(p->GetIdx(), *nbr_id)->GetFlags());

    return;
  }

  void OBBuilder::CorrectStereoBonds(OBMol &mol)
  {
    // Get CisTransStereos and make a vector of corresponding OBStereoUnits
    std::vector<OBCisTransStereo*> cistrans, newcistrans;
    OBStereoUnitSet sgunits;
    std::vector<OBGenericData*> vdata = mol.GetAllData(OBGenericDataType::StereoData);
    OBStereo::Ref bond_id;
    for (std::vector<OBGenericData*>::iterator data = vdata.begin(); data != vdata.end(); ++data)
      if (((OBStereoBase*)*data)->GetType() == OBStereo::CisTrans) {
        OBCisTransStereo *ct = dynamic_cast<OBCisTransStereo*>(*data);
        if (ct->GetConfig().specified) {
          cistrans.push_back(ct);
          bond_id = mol.GetBond(mol.GetAtomById(ct->GetConfig().begin),
                                mol.GetAtomById(ct->GetConfig().end))->GetId();
          sgunits.push_back(OBStereoUnit(OBStereo::CisTrans, bond_id));
        }
      }

    // Perceive CisTransStereos
    newcistrans = CisTransFrom3D(&mol, sgunits, false);

    // Compare and correct if necessary
    double newangle, angle;
    OBAtom *a, *b, *c, *d;
    std::vector<OBCisTransStereo*>::iterator origct, newct;
    for (origct=cistrans.begin(), newct=newcistrans.begin(); origct!=cistrans.end(); ++origct, ++newct) {
      OBCisTransStereo::Config config = (*newct)->GetConfig(OBStereo::ShapeU);
      if ((*origct)->GetConfig(OBStereo::ShapeU) != config) { // Wrong cis/trans stereochemistry

        // refs[0]            refs[3]
        //        \          /
        //         begin==end
        //        /          \
        // refs[1]            refs[2]

        a = mol.GetAtomById(config.refs[0]);
        b = mol.GetAtomById(config.begin);
        c = mol.GetAtomById(config.end);
        if (config.refs[3] != OBStereo::ImplicitRef)
          d = mol.GetAtomById(config.refs[3]);
        else
          d = mol.GetAtomById(config.refs[2]);
        angle = mol.GetTorsion(a, b, c, d); // In degrees
        newangle = angle * DEG_TO_RAD + M_PI; // flip the bond by 180 deg (PI radians)
        mol.SetTorsion(a, b, c, d, newangle); // In radians
      }
    }
  }

  bool OBBuilder::IsSpiroAtom(unsigned long atomId, OBMol &mol)
  {
    OBMol workmol = mol; // Make a copy (this invalidates Ids, but not Idxs)
    OBAtom* watom = workmol.GetAtom(mol.GetAtomById(atomId)->GetIdx());
    if (watom->GetHvyValence() != 4) // QUESTION: Do I need to restrict it further?
      return false;

    std::vector <OBBond*> bonds;
    std::vector <OBBond*>::iterator bond_it;
    FOR_BONDS_OF_ATOM(b, watom) {
      if (!b->IsInRing())
        return false;
      else
        bonds.push_back(&(*b));
    }
    // Removing the four bonds should partition the molecule into 3 fragments
    // ASSUMPTION: Molecule is a contiguous fragment to begin with
    for(bond_it=bonds.begin(); bond_it != bonds.end(); ++bond_it)
      workmol.DeleteBond(*bond_it);
    if (workmol.Separate().size() != 3)
      return false;

    return true;
  }

  void OBBuilder::FlipSpiro(OBMol &mol, int idx)
  {
    OBAtom *p = mol.GetAtom(idx); // The pivot atom

    // Find two of the four bonds that are in the same ring
    vector<int> nbrs;
    FOR_NBORS_OF_ATOM(nbr, p)
      nbrs.push_back(nbr->GetIdx());
    //assert(nbrs.size() == 4);

    // Which neighbour is in the same ring as nbrs[0]? The answer is 'ringnbr'.
    vector<int> children;
    mol.FindChildren(children, idx, nbrs[0]);
    int ringnbr = -1;
    for (vector<int>::iterator nbr=nbrs.begin() + 1; nbr!=nbrs.end(); ++nbr)
      if (find(children.begin(), children.end(), *nbr) != children.end()) {
        ringnbr = *nbr;
        break;
      }
    //assert (ringnbr != -1);

    // Split into a fragment to be flipped
    OBMol workMol = mol;
    workMol.DeleteBond(workMol.GetBond(idx, nbrs[0]));
    workMol.DeleteBond(workMol.GetBond(idx, ringnbr));
    OBBitVec fragment = GetFragment(workMol.GetAtom(nbrs[0]));

    // Translate fragment to origin
    vector3 posP = p->GetVector();
    for (unsigned int i = 1; i <= workMol.NumAtoms(); ++i)
      if (fragment.BitIsSet(i))
        workMol.GetAtom(i)->SetVector(workMol.GetAtom(i)->GetVector() - posP);

    // Rotate 180 deg around the bisector of nbrs[0]--p--ringnbr
    vector3 bond1 = posP - mol.GetAtom(nbrs[0])->GetVector();
    vector3 bond2 = posP - mol.GetAtom(ringnbr)->GetVector();
    bond1.normalize();
    bond2.normalize();
    vector3 axis = bond1 + bond2; // The bisector of bond1 and bond2

    matrix3x3 mat;
    mat.RotAboutAxisByAngle(axis, 180);
    vector3 tmpvec;
    for (unsigned int i = 1; i <= workMol.NumAtoms(); ++i)
      if (fragment.BitIsSet(i)) {
        tmpvec = workMol.GetAtom(i)->GetVector();
        tmpvec *= mat;
        workMol.GetAtom(i)->SetVector( tmpvec );
      }

    // Set the coordinates of the original molecule using those of workmol
    for (unsigned int i = 1; i <= workMol.NumAtoms(); ++i)
      if (fragment.BitIsSet(i))
        mol.GetAtom(i)->SetVector(workMol.GetAtom(i)->GetVector() + posP);
  }

  void OBBuilder::CorrectStereoAtoms(OBMol &mol)
  {
    // Get TetrahedralStereos and make a vector of corresponding OBStereoUnits
    std::vector<OBTetrahedralStereo*> tetra, newtetra;
    OBStereoUnitSet sgunits;
    std::vector<OBGenericData*> vdata = mol.GetAllData(OBGenericDataType::StereoData);
    OBStereo::Ref atom_id;
    for (std::vector<OBGenericData*>::iterator data = vdata.begin(); data != vdata.end(); ++data)
      if (((OBStereoBase*)*data)->GetType() == OBStereo::Tetrahedral) {
        OBTetrahedralStereo *th = dynamic_cast<OBTetrahedralStereo*>(*data);
        if (th->GetConfig().specified) {
          tetra.push_back(th);
          atom_id = th->GetConfig().center;
          sgunits.push_back(OBStereoUnit(OBStereo::Tetrahedral, atom_id));
        }
      }

    // Perceive TetrahedralStereos
    newtetra = TetrahedralFrom3D(&mol, sgunits, false);

    // Identify any ring stereochemistry and whether it is right or wrong
    // - ring stereo involves 3 ring bonds, or 4 ring bonds but the
    //   atom must not be spiro
    OBAtom* center;
    bool existswrongstereo = false; // Is there at least one wrong ring stereo?
    typedef std::pair<OBStereo::Ref, bool> IsThisStereoRight;
    std::vector<IsThisStereoRight> ringstereo;
    std::vector<OBTetrahedralStereo*> nonringtetra, nonringnewtetra;
    std::vector<OBTetrahedralStereo*>::iterator origth, newth;
    for (origth=tetra.begin(), newth=newtetra.begin(); origth!=tetra.end(); ++origth, ++newth) {
      OBTetrahedralStereo::Config config = (*newth)->GetConfig(OBStereo::Clockwise, OBStereo::ViewFrom);

      center = mol.GetAtomById(config.center);
      int ringbonds = 0;
      FOR_BONDS_OF_ATOM(b, center)
        if (b->IsInRing())
          ringbonds++;

      if (ringbonds == 3 || (ringbonds==4 && !OBBuilder::IsSpiroAtom(config.center, mol))) {
        bool rightstereo = (*origth)->GetConfig(OBStereo::Clockwise, OBStereo::ViewFrom) == config;
        ringstereo.push_back(IsThisStereoRight(config.center, rightstereo));
        if (!rightstereo)
          existswrongstereo = true;
      }
      else { // A non-ring stereocenter
        nonringtetra.push_back(*origth);
        nonringnewtetra.push_back(*newth);
      }
    }

    if (existswrongstereo) {
      // Fix ring stereo
      OBStereo::Refs unfixed;
      bool inversion = FixRingStereo(ringstereo, mol, unfixed);

      // Output warning message if necessary
      if (unfixed.size() > 0) {
        stringstream errorMsg;
        errorMsg << "Could not correct " << unfixed.size() << " stereocenter(s) in this molecule (" << mol.GetTitle() << ")";
        errorMsg << std::endl << "  with Atom Ids as follows:";
        for (OBStereo::RefIter ref=unfixed.begin(); ref!=unfixed.end(); ++ref)
          errorMsg << " " << *ref;
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
      }

      // Reperceive non-ring TetrahedralStereos if an inversion occured
      if (inversion) {
        sgunits.clear();
        for (origth = nonringtetra.begin(); origth != nonringtetra.end(); ++origth)
          sgunits.push_back(OBStereoUnit(OBStereo::Tetrahedral, (*origth)->GetConfig().center));
        nonringnewtetra = TetrahedralFrom3D(&mol, sgunits, false);
      }
    }

    // Correct the non-ring stereo
    for (origth=nonringtetra.begin(), newth=nonringnewtetra.begin(); origth!=nonringtetra.end(); ++origth, ++newth) {
      OBTetrahedralStereo::Config config = (*newth)->GetConfig(OBStereo::Clockwise, OBStereo::ViewFrom);
      if ((*origth)->GetConfig(OBStereo::Clockwise, OBStereo::ViewFrom) != config) {
        // Wrong tetrahedral stereochemistry

        // Try to find two non-ring bonds
        center = mol.GetAtomById(config.center);
        vector<unsigned int> idxs;
        FOR_BONDS_OF_ATOM(b, center)
          if (!b->IsInRing())
            idxs.push_back(b->GetNbrAtom(center)->GetIdx());

        if (idxs.size() == 0 && OBBuilder::IsSpiroAtom(config.center, mol))
          FlipSpiro(mol, center->GetIdx());
        else if (idxs.size() >= 2)
          Swap(mol, center->GetIdx(), idxs.at(0), center->GetIdx(), idxs.at(1));
        else {
          // It will only reach here if it can only find one non-ring bond
          // -- this is the case if the other non-ring bond is an implicit H
          //    or a lone pair
          // Solution: Find where a new bond vector would be placed, and
          //           replace the atom's coordinates with these
          OBAtom* non_ring_atom = mol.GetAtom(idxs.at(0));
          OBBond* non_ring_bond = mol.GetBond(center, non_ring_atom);
          vector3 newcoords = OBBuilder::GetNewBondVector(center, non_ring_bond->GetLength());
          SwapWithVector(mol, center->GetIdx(), idxs.at(0), center->GetIdx(), newcoords);
        }
      }
    }
  }

  bool OBBuilder::FixRingStereo(std::vector<std::pair<OBStereo::Ref, bool> > atomIds, OBMol &mol,
                                OBStereo::Refs &unfixedcenters)
  {
    bool inversion = false;
    if (atomIds.size() == 0) return inversion;

    // Have we dealt with a particular ring stereo? (Indexed by Id)
    OBBitVec seen;

    for(unsigned int n=0; n<atomIds.size(); ++n) {
      // Keep looping until you come to an unseen wrong stereo
      if (seen.BitIsSet(atomIds[n].first) || atomIds[n].second) continue;

      OBBitVec fragment; // Indexed by Id
      AddRingNbrs(fragment, mol.GetAtomById(atomIds[n].first), mol);

      // Which ring stereos does this fragment contain, and
      // are the majority of them right or wrong?
      OBStereo::Refs wrong, right;
      for (unsigned int i=0; i<atomIds.size(); ++i)
        if (fragment.BitIsSet(atomIds[i].first)) {
          if (atomIds[i].second)
            right.push_back(atomIds[i].first);
          else
            wrong.push_back(atomIds[i].first);
          seen.SetBitOn(atomIds[i].first);
        }

      if (right > wrong) { // Inverting would make things worse!
        unfixedcenters.insert(unfixedcenters.end(), wrong.begin(), wrong.end());
        continue;
      }
      unfixedcenters.insert(unfixedcenters.end(), right.begin(), right.end());

      // Invert the coordinates (QUESTION: should I invert relative to the centroid?)
      inversion = true;
      FOR_ATOMS_OF_MOL(a, mol)
        if (fragment.BitIsSet(a->GetId()))
          a->SetVector( - a->GetVector());

      // Add neighbouring bonds back onto the fragment
      // TODO: Handle spiro
      std::vector<OBBond*> reconnect;
      FOR_ATOMS_OF_MOL(a, mol)
        if (fragment.BitIsSet(a->GetId()))
          FOR_BONDS_OF_ATOM(b, &*a)
            if (!b->IsInRing())
              reconnect.push_back(&*b);

      for (std::vector<OBBond*>::iterator bi=reconnect.begin(); bi!=reconnect.end(); ++bi) {
        OBBond* b = *bi;
        int bo = b->GetBondOrder();
        int begin = b->GetBeginAtomIdx();
        int end = b->GetEndAtomIdx();
        mol.DeleteBond(b);
        OBBuilder::Connect(mol, begin, end, bo);
      }
    }

    return inversion;
  }

  void OBBuilder::AddRingNbrs(OBBitVec &fragment, OBAtom *atom, OBMol &mol)
  {
    // Add the nbrs to the fragment, but don't add the neighbours of a spiro atom.
    FOR_NBORS_OF_ATOM (nbr, atom) {
      if (mol.GetBond(&*nbr, atom)->IsInRing() && !fragment.BitIsSet(nbr->GetId())
          && !OBBuilder::IsSpiroAtom(atom->GetId(), mol)) {
        fragment.SetBitOn(nbr->GetId());
        AddRingNbrs(fragment, &*nbr, mol);
      }
    }
  }

} // end namespace OpenBabel


//! \file builder.cpp
//! \brief Handle OBBuilder class

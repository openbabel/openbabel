/**********************************************************************
distgeom.cpp - Distance Geometry generation and sampling

  Copyright (C) 2011 by Tim Vandermeersch
  Copyright (C) 2012 by Geoffrey Hutchison

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

#ifdef HAVE_EIGEN

#include <openbabel/distgeom.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/ring.h>
#include <openbabel/obiter.h>
#include <openbabel/builder.h>
#include <openbabel/elements.h>
#include <openbabel/generic.h>
#include "rand.h"
#include <LBFGS.h>
#include "stereo/gen3dstereohelper.h"

#include <openbabel/stereo/stereo.h>
#include <openbabel/stereo/cistrans.h>
#include <openbabel/stereo/tetrahedral.h>
#include <openbabel/obconversion.h>

#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>

using namespace std;

#define DIST12_TOL   0.01f
#define DIST13_TOL   0.03f
#define DIST14_TOL   0.05f
#define DIST15_TOL   0.07f

#pragma warning(disable : 4244) // warning C4244: '=' : conversion from 'double' to 'float', possible loss of data
#pragma warning(disable : 4305) // warning C4305: '*=' : truncation from 'double' to 'float'

namespace OpenBabel {

  class DistanceGeometryPrivate {
  public:
    DistanceGeometryPrivate(const unsigned int N)
    {
      bounds = Eigen::MatrixXf(static_cast<int>(N), static_cast<int>(N));
      preMet = Eigen::MatrixXf(bounds);
      debug = false;
      success = false;
    }
    ~DistanceGeometryPrivate()
    { }

    // TODO: Check for overflow on i, j
    void SetLowerBounds(int i, int j, float value)
    {
      if (i > j)
        bounds(i, j) = value;
      else
        bounds(j, i) = value;
    }
    void SetUpperBounds(int i, int j, float value)
    {
      if (i < j)
        bounds(i, j) = value;
      else
        bounds(j, i) = value;
    }
    float GetLowerBounds(int i, int j)
    {
      if (i > j)
        return bounds(i, j);
      else
        return bounds(j, i);
    }
    float GetUpperBounds(int i, int j)
    {
      if (i < j)
        return bounds(i, j);
      else
        return bounds(j, i);
    }
    float GetAvgBounds(int i, int j)
    {
      float lb = GetLowerBounds(i, j);
      return (GetUpperBounds(i, j) - lb) / 2.0 + lb;
    }

    Eigen::MatrixXf bounds, preMet;
    bool debug; double maxBoxSize;
    OBGen3DStereoHelper stereoHelper;
    bool success;
  };


  OBDistanceGeometry::OBDistanceGeometry(): _d(nullptr) {}

  OBDistanceGeometry::OBDistanceGeometry(const OBMol &mol, bool useCurrentGeometry): _d(nullptr)
  {
    Setup(mol, useCurrentGeometry);
  }

  OBDistanceGeometry::~OBDistanceGeometry()
  {
    if (_d != nullptr)
      delete _d;
  }

  float OBDistanceGeometry::GetUpperBounds(int i, int j) {
      return _d->GetUpperBounds(i, j);
  }
  float OBDistanceGeometry::GetLowerBounds(int i, int j) {
        return _d->GetLowerBounds(i, j);
  }
  bool OBDistanceGeometry::Setup(const OBMol &mol, bool useCurrentGeometry)
  {
    if (_d != nullptr)
      delete _d;
    // TODO: add IsSetupNeeded() like OBForceField to prevent duplication of work

    dim = 4;
    _mol = mol;

    _d = new DistanceGeometryPrivate(mol.NumAtoms());

    _d->stereoHelper.Setup(&_mol);

    _mol.SetDimension(3);
    _vdata = _mol.GetAllData(OBGenericDataType::StereoData);

    SetUpperBounds();
    // Do we use the current geometry for default 1-2 and 1-3 bounds?
    Set12Bounds(useCurrentGeometry);
    if (_d->debug) {
      cerr << endl << " 1-2 Matrix\n";
      cerr << _d->bounds << endl;
    }
    Set13Bounds(useCurrentGeometry);
    if (_d->debug) {
      cerr << endl << " 1-3 Matrix\n";
      cerr << _d->bounds << endl;
    }

    vector<OBRing*> rlist = _mol.GetSSSR();
    if (rlist.size() > 0)
      SetAromaticRingBounds();

    Set14Bounds();
    if (_d->debug) {
      cerr << endl << " 1-4 Matrix\n";
      cerr << _d->bounds << endl;
    }
    Set15Bounds();
    SetLowerBounds();
    TriangleSmooth();
    _d->preMet = _d->bounds; // make a copy before metrization
    if (_d->debug) {
      cerr << endl << " Smoothed Matrix\n";
      cerr << _d->bounds << endl;
    }

    // RDKit: Code/DistGeom/ChiralViolationContrib.cpp
    OBStereoFacade facade(&_mol);
    FOR_ATOMS_OF_MOL(atom, _mol) {
      if (facade.HasTetrahedralStereo(atom->GetId())) {
        OBTetrahedralStereo *ts = facade.GetTetrahedralStereo(atom->GetId());
        OBTetrahedralStereo::Config config = ts->GetConfig();
        vector<unsigned long> nbrs;

        unsigned long centerIdx = _mol.GetAtomById(config.center)->GetIdx()-1;

        if (config.from == OBStereo::ImplicitRef)
          nbrs.push_back(centerIdx);
        else
          nbrs.push_back(_mol.GetAtomById(config.from)->GetIdx()-1);

        for(size_t i=0; i<config.refs.size(); i++) {
          if (config.refs[i] == OBStereo::ImplicitRef)
            nbrs.push_back(centerIdx);
          else
            nbrs.push_back(_mol.GetAtomById(config.refs[i])->GetIdx()-1);
        }

        if(config.winding == OBStereo::Clockwise) {
          TetrahedralInfo ti(centerIdx, nbrs, -100.0, -5.0);
          _stereo.push_back(ti);
        } else if(config.winding == OBStereo::AntiClockwise) {
          TetrahedralInfo ti(centerIdx, nbrs, 5.0, 100.0);
          _stereo.push_back(ti);
        }
      }
    }

    // Planar constraint for aromatic ring
    if (rlist.size() > 0) {
      FOR_RINGS_OF_MOL(r, _mol) {
        int size = r->Size();
        if (!r->IsAromatic()) continue;

        std::vector<int> path = r->_path;
        for (int a = 0; a < size; ++a) {
          unsigned long idx = path[a];
          OBAtom* atom = _mol.GetAtom(idx);

          vector<unsigned long> nbrs(1, idx-1);
          FOR_NBORS_OF_ATOM(b, &*atom) {
            nbrs.push_back(b->GetIdx()-1);
          }
          if(nbrs.size() == 4) {
            TetrahedralInfo ti(idx, nbrs, 0.0, 0.0);
            _stereo.push_back(ti);
          }
        }
      }
    }
    return true;
  }

  // Set the default bounds to a maximum distance
  void OBDistanceGeometry::SetUpperBounds()
  {
    if (!_d)
      return;

    unsigned int N = _mol.NumAtoms();
    float maxDist = N*1.5f; // if, somehow all atoms are in a linear chain

    // If we're in a unit cell, the maximum distance is 1/2 the longest body diagonal
    //   (remember that the unit cell wraps around)
    OBUnitCell* pUC = (OBUnitCell*)_mol.GetData(OBGenericDataType::UnitCell);
    if (pUC != nullptr) {
      vector<vector3> cellVectors = pUC->GetCellVectors();

      if (cellVectors.size() == 3) {
        vector3 diagonal = cellVectors[0] + cellVectors[1] + cellVectors[2];
        maxDist = diagonal.length() / 2.0;
      }
    }

    for (unsigned int i = 0; i < N; ++i) {
      // set diagonal to zero
      _d->bounds(i, i) = 0.0f;
      for (unsigned int j = i + 1; j < N; ++j)
        {
          _d->SetLowerBounds(i, j, 0.0); // for now -- allows us to check for visits
          _d->SetUpperBounds(i, j, maxDist);
        }
    }
  }

  void OBDistanceGeometry::Set12Bounds(bool useGeom)
  {
    float length;
    FOR_BONDS_OF_MOL(b, _mol) {
      unsigned int i = b->GetBeginAtomIdx() - 1;
      unsigned int j = b->GetEndAtomIdx() - 1;
      if (useGeom) {
        length = b->GetLength();
        // Allow a tiny amount of slop
        _d->SetLowerBounds(i, j, length - DIST12_TOL );
        _d->SetUpperBounds(i, j, length + DIST12_TOL );
      } else {
        length = b->GetEquibLength(); // ideal length
        // Allow slightly more slop, since that's empirical
        _d->SetLowerBounds(i, j, length - DIST12_TOL*1.5 );
        _d->SetUpperBounds(i, j, length + DIST12_TOL*1.5 );
      }
    }
  }

  // Helper for calculating 13 distances by cosine rule
  //  useful for 14 and 15 relations too
  inline double Calculate13Distance(double ab, double bc, double angle)
  {
    return sqrt(SQUARE(ab) + SQUARE(bc) - 2.0*ab*bc*cos(angle));
  }

  inline double Calculate13Angle(double a, double b, double c)
  {
    double cosine = (SQUARE(a) + SQUARE(b) - SQUARE(c)) / (2.0*a*b);

    // Handle cases where cosine is outside [-1, 1] range.
    if (cosine > 1.0)
      return 0.0;
    if (cosine < -1.0)
      return M_PI;

    return acos(cosine);
  }

  // When atoms i and j are in a 1-3 relationship, the distance
  //   is calculated using the cosine rule: (upper limit ~ lower limit)
  //
  //          b_         ab: bond length
  //         /  \_       ac: bond length
  //        /A    \_     bc = sqrt(ab^2 + ac^2 - 2*ab*ac*cos(A))
  //       a--------c
  //
  void OBDistanceGeometry::Set13Bounds(bool useGeom)
  {
    float dist, rAB, rAC;
    OBAtom *a, *b, *c;
    unsigned int i, j;
    // Angle is    b
    //            /
    //           a----c
    // with a as the vertex
    FOR_ANGLES_OF_MOL(angle, _mol) {
      a = _mol.GetAtom((*angle)[0] + 1);
      b = _mol.GetAtom((*angle)[1] + 1);
      c = _mol.GetAtom((*angle)[2] + 1);
      if (b->GetBond(c) != nullptr)
        continue;
      i = (*angle)[1];
      j = (*angle)[2];

      // Easy case -- use current geometry
      if (useGeom) {
        dist = b->GetDistance(c);
        _d->SetLowerBounds(i, j, dist - DIST13_TOL);
        _d->SetUpperBounds(i, j, dist + DIST13_TOL);
      } else {
        // Guess angle based on central atom
        // TODO: refine if this angle is in a ring
        float theta = 109.5f * DEG_TO_RAD; // in radians

        // If the two endpoints are in the same ring
        //  AND the vertex is in some ring, they're all in it
        int ringSize = AreInSameRing(b, c);
        if (a->IsInRing() && ringSize != 0)
          {
            // Atom is sp2 hybrid, so assume planar
            if (a->IsAromatic() || a->GetHyb() == 2 || ringSize <= 4) {
              theta = 180.0f - (360.0f/float(ringSize));
              theta *= DEG_TO_RAD;
            }
            // Atom is sp3, so approximate
            else if (a->GetHyb() == 3) {
              switch(ringSize) {
              case 3:
                theta = 60.0f * DEG_TO_RAD;
                break;
              case 4:
                theta = 90.0f * DEG_TO_RAD;
                break;
              case 5:
                theta = 104.0f * DEG_TO_RAD;
                break;
              default:
                theta = 109.5f * DEG_TO_RAD;
              }
            } // end sp3
          }
        else { // not all in the same ring
          switch (a->GetHyb()) {
          case 1:
            theta = 180.0f * DEG_TO_RAD;
            break;
          case 2:
            theta = 120.0f * DEG_TO_RAD;
            break;
          case 3:
          default:
            theta = 109.5f * DEG_TO_RAD;
          } // end switch
        }

        // cosine rule
        // Get the 12 distances, since we don't have geometry yet
        // (remember "A" is the vertex in Open Babel
        rAB = _d->GetLowerBounds((*angle)[0], (*angle)[1]) + DIST12_TOL;
        rAC = _d->GetLowerBounds((*angle)[0], (*angle)[2]) + DIST12_TOL;

        dist = Calculate13Distance(rAB, rAC, theta);
        _d->SetLowerBounds(i, j, dist - DIST13_TOL);
        _d->SetUpperBounds(i, j, dist + DIST13_TOL);
      } //end unknown geometry
    }
  }

  // Helper for calculating 14 distances when in the cis conformation
  //      a       d      ab, bc, cd: bond lengths
  //       \ B C /
  //        b---c        ad = bc + ab*cos(180-B) + cd*cos(180-C)
  // ANGLES are in RADIANS!
  inline double Calculate14DistCis(double ab, double bc, double cd,
                                   double B, double C) {
    double lB = M_PI - B;
    double lC = M_PI - C;
    return bc + ab*cos(lB) + cd*cos(lC);
  }

  // Helper for calculating 14 distances when in the trans conformation
  //      a
  //       \ B           delta_x = bc + ab*cos(180-B) + cd*cos(180-C)
  //        b---c        delta_y = ab*sin(180-B) + cd*sin(180-C)
  //           C \       .
  //              d      ad = sqrt(delta_x^2 + delta_y^2)
  // ANGLES are in RADIANS!
  inline double Calculate14DistTrans(double ab, double bc, double cd,
                                     double B, double C) {
    double lB = M_PI - B;
    double lC = M_PI - C;
    double dx = bc + ab*cos(lB) + cd*cos(lC);
    double dy =      ab*sin(lB) + cd*sin(lC);
    return sqrt(SQUARE(dx) + SQUARE(dy));
  }

  void OBDistanceGeometry::SetAromaticRingBounds()
  {
    // Set all pairwise interactions around aromatic rings
    FOR_RINGS_OF_MOL(r, _mol) {
      int size = r->Size();
      if (size != 6 || !r->IsAromatic())
        continue;

      // Aromatic rings must be planar, so atoms in the ring = regular polygon
      float bondDist, radius;

      // We should have 1-2 and 1-3 distances set exactly.
      // So we need to set 1-4 (e.g., para)
      std::vector<int> path = r->_path;
      int a, b, c, d; // entries into path vector
      for (a = 0; a < size; ++a) {
        b = (a + 1) % size;
        c = (a + 2) % size;
        d = (a + 3) % size;

        // get an average distance (e.g., heteroatoms)
        bondDist = _d->GetAvgBounds(path[a] - 1, path[b] - 1)
          + _d->GetAvgBounds(path[b] - 1, path[c] - 1)
          + _d->GetAvgBounds(path[c] - 1, path[d] - 1);
        bondDist = bondDist / 3.0f;

        // so the bonds are the sides of the regular polygon
        // and the circumradius will be half the distance across
        // http://en.wikipedia.org/wiki/Regular_polygon#Circumradius
        radius = bondDist / (2.0f * sin(M_PI / size));

        float lBounds = 2.0f*radius - DIST14_TOL;
        float uBounds = 2.0f*radius + DIST14_TOL;

        _d->SetLowerBounds(path[a] - 1, path[d] - 1, lBounds);
        _d->SetUpperBounds(path[a] - 1, path[d] - 1, uBounds);

      } // done with path
    } // done with rings
  }

  // - when atoms i and j are in a 1-4 relationship, the lower distance
  //   limit is calculated using a torsional angle of 0.0. The upper limit
  //   is calculated using a torsion angle of 180.0.
  void OBDistanceGeometry::Set14Bounds()
  {
    float rAB, rBC, rCD;
    float rAC, rBD, B, C;
    float lBounds, uBounds;
    OBAtom *a, *b, *c, *d;
    OBBond *bc;

    // Loop through all torsions first
    FOR_TORSIONS_OF_MOL(t, _mol) {
      a = _mol.GetAtom((*t)[0] + 1);
      b = _mol.GetAtom((*t)[1] + 1);
      c = _mol.GetAtom((*t)[2] + 1);
      d = _mol.GetAtom((*t)[3] + 1);

      if (a->GetBond(d) != nullptr)
        continue; // these are bonded
      if (_d->GetLowerBounds((*t)[0], (*t)[3]) > 0.01) // we visited this
        continue;

      // We want to know the a-d distance here
      // So ...
      // Get the 12 distances
      rAB = _d->GetLowerBounds((*t)[0], (*t)[1]) + DIST12_TOL;
      rBC = _d->GetLowerBounds((*t)[1], (*t)[2]) + DIST12_TOL;
      rCD = _d->GetLowerBounds((*t)[2], (*t)[3]) + DIST12_TOL;

      // Get the 13 angles
      rAC = _d->GetLowerBounds((*t)[0], (*t)[2]) + DIST13_TOL;
      rBD = _d->GetLowerBounds((*t)[1], (*t)[3]) + DIST13_TOL;
      B = Calculate13Angle(rAB, rBC, rAC);
      C = Calculate13Angle(rBC, rCD, rBD);

      // default bounds
      lBounds = Calculate14DistCis(rAB, rBC, rCD, B, C);
      uBounds = Calculate14DistTrans(rAB, rBC, rCD, B, C);

      // TODO: special cases
      _d->SetLowerBounds((*t)[0], (*t)[3], lBounds - DIST14_TOL);
      _d->SetUpperBounds((*t)[0], (*t)[3], uBounds + DIST14_TOL);
    }

    // OK, check and correct double bond cis/trans stereochemistry
    // Get CisTransStereos and make a vector of corresponding OBStereoUnits
    OBStereoUnitSet sgunits;
    std::vector<OBGenericData*> vdata = _mol.GetAllData(OBGenericDataType::StereoData);
    for (std::vector<OBGenericData*>::iterator data = vdata.begin(); data != vdata.end(); ++data)
      if (((OBStereoBase*)*data)->GetType() == OBStereo::CisTrans) {
        OBCisTransStereo *ct = dynamic_cast<OBCisTransStereo*>(*data);
        if (ct->GetConfig().specified) {
          // OK, get the central bond (bc) and check all the bonded atoms for proper stereo
          b = _mol.GetAtomById(ct->GetConfig().begin);
          c = _mol.GetAtomById(ct->GetConfig().end);
          FOR_NBORS_OF_ATOM(a, b) {
            if (a->GetIdx() == c->GetIdx())
              continue;
            FOR_NBORS_OF_ATOM(d, c) {
              if (d->GetIdx() == b->GetIdx())
                continue;

              float lBounds = _d->GetLowerBounds(a->GetIdx() - 1, d->GetIdx() - 1) + DIST14_TOL;
              float uBounds = _d->GetUpperBounds(a->GetIdx() - 1, d->GetIdx() - 1) - DIST14_TOL;
              if (ct->IsTrans(a, d)) {
                // lower bounds should be trans (current upper bounds)
                _d->SetLowerBounds(a->GetIdx() - 1, d->GetIdx() - 1, uBounds - DIST14_TOL);
              } else if (ct->IsCis(a, d)) {
                // upper bounds should be cis (current lower bounds)
                _d->SetUpperBounds(a->GetIdx() - 1, d->GetIdx() - 1, lBounds + DIST14_TOL);
              }
            } // neighbors of cis/trans c
          } // neighbors of cis/trans b
        }
      } // iterate through cis/trans


    // Now correct ring bonds -- if bc is a ring bond, and a and d are in the same ring, then torsion should be "cis-oid"
    //  If b=c is a double bond, set to exactly cis
    //  Otherwise, give a bit of slack for upper bound
    FOR_RINGS_OF_MOL(r, _mol) {
      int size = r->Size();
      if (size < 4)
        continue;

      std::vector<int> path = r->_path;
      int a, b, c, d; // entries into path vector
      for (a = 0; a < size; ++a) {
        b = (a + 1) % size;
        c = (a + 2) % size;
        d = (a + 3) % size;

        float lBounds = _d->GetLowerBounds(path[a] - 1, path[d] - 1) + DIST14_TOL;
        float uBounds = _d->GetUpperBounds(path[a] - 1, path[d] - 1) - DIST14_TOL;

        bc = _mol.GetBond(path[b], path[c]);
        if (bc->IsAromatic() || bc->GetBondOrder() == 2) {
          uBounds = lBounds + DIST14_TOL;
          // Correct non-ring neighbors too -- these should be out of the ring
          FOR_NBORS_OF_ATOM(nbr, _mol.GetAtom(path[b])) {
            if (nbr->GetIdx() == static_cast<unsigned>(path[a])
            || nbr->GetIdx() == static_cast<unsigned>(path[c]))
              continue;
            // This atom should be trans to atom D
            _d->SetLowerBounds(nbr->GetIdx() - 1, path[d] - 1,
                               _d->GetUpperBounds(nbr->GetIdx() - 1, path[d] - 1) - DIST14_TOL);
          }
          FOR_NBORS_OF_ATOM(nbr, _mol.GetAtom(path[c])) {
            if (nbr->GetIdx() == static_cast<unsigned>(path[d])
            || nbr->GetIdx() == static_cast<unsigned>(path[b]))
              continue;
            // This atom should be trans to atom A
            _d->SetLowerBounds(nbr->GetIdx() - 1, path[a] - 1,
                               _d->GetUpperBounds(nbr->GetIdx() - 1, path[a] - 1) - DIST14_TOL);
          }

        } else if (bc->GetBondOrder() == 1) {
          // Could be anywhere from pure-cis to halfway to trans
          uBounds = _d->GetAvgBounds(path[a] - 1, path[d] - 1);

          // Adjust the non-ring neighbors too -- these should be out of the ring (i.e., trans-oid)
          // Correct non-ring neighbors too -- these should be out of the ring
          FOR_NBORS_OF_ATOM(nbr, _mol.GetAtom(path[b])) {
            if (nbr->GetIdx() == static_cast<unsigned>(path[a])
            || nbr->GetIdx() == static_cast<unsigned>(path[c]))
              continue;
            // This atom should be quasi-trans to atom D
            _d->SetLowerBounds(nbr->GetIdx() - 1, path[d] - 1, _d->GetAvgBounds(nbr->GetIdx() - 1, path[d] - 1));
          }
          FOR_NBORS_OF_ATOM(nbr, _mol.GetAtom(path[c])) {
            if (nbr->GetIdx() == static_cast<unsigned>(path[d])
             || nbr->GetIdx() == static_cast<unsigned>(path[b]))
              continue;
            // This atom should be quasi-trans to atom A
            _d->SetLowerBounds(nbr->GetIdx() - 1, path[a] - 1, _d->GetAvgBounds(nbr->GetIdx() - 1, path[a] - 1));
          }
        }
        // New upper bounds for a-b-c-d
        _d->SetUpperBounds(path[a] - 1, path[d] - 1, uBounds);

      } // done with path
    } // done with rings

    // TODO: More special cases
  }

  // Helper for calculating 15 distances when in the all-cis conformation
  //   also works for trans-cis
  //             e
  //              \
  //      a     D d      ab, bc, cd: bond lengths
  //       \ B C /
  //        b---c        ad = bc + ab*cos(180-B) + cd*cos(180-C)
  // ANGLES are in RADIANS!
  inline double Calculate15DistAnyCis(double ab, double bc, double cd, double de,
                                      double B, double C, double D) {
    // Stolen and adapted from RDKit (http://rdkit.org/)
    // Covered under the BSD license
    double xad = bc - cd*cos(C) - ab*cos(B);
    double yad =      cd*sin(C) - ab*sin(B);
    double ad = sqrt(xad*xad + yad*yad);
    double cval = (cd - bc*cos(C) + ab*cos(B + C))/ad;
    if (cval > 1.0) {
      cval = 1.0;
    } else if (cval < -1.0) {
      cval = -1.0;
    }

    double angADC = acos(cval);
    double angADE = D - angADC;
    return Calculate13Distance(ad, de, angADE);
  }


  // Helper for calculating 15 distances when in the all-trans conformation
  //   also works for cis-trans
  //      a
  //       \ B           delta_x = bc + ab*cos(180-B) + cd*cos(180-C)
  //        b---c        delta_y = ab*sin(180-B) + cd*sin(180-C)
  //           C \  D    .
  //              d--e   ad = sqrt(delta_x^2 + delta_y^2)
  // ANGLES are in RADIANS!
  inline double Calculate15DistAnyTrans(double ab, double bc, double cd, double de,
                                        double B, double C, double D) {
    // Stolen and adapted from RDKit (http://rdkit.org/)
    // Covered under the BSD license
    double xad = bc - cd*cos(C) - ab*cos(B);
    double yad =      cd*sin(C) - ab*sin(B);
    double ad = sqrt(xad*xad + yad*yad);
    double cval = (cd - bc*cos(C) + ab*cos(B + C))/ad;
    if (cval > 1.0) {
      cval = 1.0;
    } else if (cval < -1.0) {
      cval = -1.0;
    }

    double angADC = acos(cval);
    double angADE = D + angADC;
    return Calculate13Distance(ad, de, angADE);
  }


  OBCisTransStereo * OBDistanceGeometry::GetCisTransStereo(OBBond *bond)
  {
    OBAtom *b, *c;
    OBBond *bc;

    OBStereoUnitSet sgunits;
    std::vector<OBGenericData*> vdata = _mol.GetAllData(OBGenericDataType::StereoData);
    for (std::vector<OBGenericData*>::iterator data = vdata.begin(); data != vdata.end(); ++data)
      if (((OBStereoBase*)*data)->GetType() == OBStereo::CisTrans) {
        OBCisTransStereo *ct = dynamic_cast<OBCisTransStereo*>(*data);
        if (ct->GetConfig().specified) {
          // OK, get the central bond (bc) and check all the bonded atoms for proper stereo
          b = _mol.GetAtomById(ct->GetConfig().begin);
          c = _mol.GetAtomById(ct->GetConfig().end);

          bc = _mol.GetBond(b, c);
          if (bc && bc->GetIdx() == bond->GetIdx())
            return ct;
        }
      }
    // didn't find anything, return NULL
    return nullptr;
  }

  // - when atoms i and j are in a 1-5 relationship, the lower distance
  //   cannot be closer than a cis-cis relationship (i.e., curling around)
  //    or farther than a fully extended trans relationship
  void OBDistanceGeometry::Set15Bounds()
  {
    float rAB, rBC, rCD, rDE;
    float rAC, rBD, rCE, rAE, rBE;
    float A, B, C, D;
    float lBounds, uBounds;
    OBAtom *a, *b, *c, *d;
    OBBond *ab, *cd;

    FOR_TORSIONS_OF_MOL(t, _mol) {
      a = _mol.GetAtom((*t)[0] + 1);
      b = _mol.GetAtom((*t)[1] + 1);
      c = _mol.GetAtom((*t)[2] + 1);
      d = _mol.GetAtom((*t)[3] + 1);

      // We want to know the a-d distance here
      // So ...
      // Get the 12 distances
      rAB = _d->GetLowerBounds((*t)[0], (*t)[1]) + DIST12_TOL;
      rBC = _d->GetLowerBounds((*t)[1], (*t)[2]) + DIST12_TOL;
      rCD = _d->GetLowerBounds((*t)[2], (*t)[3]) + DIST12_TOL;

      // Get the 13 angles
      rAC = _d->GetLowerBounds((*t)[0], (*t)[2]) + DIST13_TOL;
      rBD = _d->GetLowerBounds((*t)[1], (*t)[3]) + DIST13_TOL;
      B = Calculate13Angle(rAB, rBC, rAC);
      C = Calculate13Angle(rBC, rCD, rBD);

      OBCisTransStereo *stereo = nullptr;

      // For neighbors of d
      //  Actually depends on stereo of bond C-D
      cd = _mol.GetBond(c, d);
      if (cd && cd->GetBondOrder() == 2 && !cd->IsAromatic()) {
        stereo = GetCisTransStereo(cd);
      }
      FOR_NBORS_OF_ATOM(e, d) {
        if (_mol.GetBond(a, &*e) != nullptr)
          continue; // Already handled by 1,2 interaction
        if (_d->GetLowerBounds((*t)[0], e->GetIdx() - 1) > 0.01) // we visited this
          continue;

        rDE = _d->GetLowerBounds((*t)[3], e->GetIdx() - 1) + DIST12_TOL;
        rCE = _d->GetLowerBounds((*t)[2], e->GetIdx() - 1) + DIST12_TOL;
        D = Calculate13Angle(rCD, rDE, rCE);

        // default bounds
        lBounds = Calculate15DistAnyCis(rAB, rBC, rCD, rDE, B, C, D);
        uBounds = Calculate15DistAnyTrans(rAB, rBC, rCD, rDE, B, C, D);

        // Check stereochemistry
        if (stereo && stereo->IsCis(b->GetId(), e->GetId()))
          uBounds = lBounds; // Must be cis
        if (stereo && stereo->IsTrans(b->GetId(), e->GetId()))
          lBounds = uBounds; // Must be trans

        // Correcting ring shapes -- should be mostly cisoid
        if (AreInSameRing(a, &*e))
          uBounds = (lBounds + uBounds) / 2.0;
        // TODO.. set exo and endo bonds (if needed)

        if (_d->GetLowerBounds((*t)[0], e->GetIdx() - 1) < lBounds)
          _d->SetLowerBounds((*t)[0], e->GetIdx() - 1, lBounds - DIST15_TOL);
        _d->SetUpperBounds((*t)[0], e->GetIdx() - 1, uBounds + DIST15_TOL);
      }

      // OK now for neighbors of a (i.e., z-a-b-c-d)
      //  Now depends on stereo of bond a-b
      stereo = nullptr; // reset
      ab = _mol.GetBond(a, b);
      if (ab && ab->GetBondOrder() == 2 && !ab->IsAromatic()) {
        stereo = GetCisTransStereo(ab);
      }
      FOR_NBORS_OF_ATOM(z, a) {
        if (_mol.GetBond(d, &*z) != nullptr)
          continue; // Already handled by 1,2 interaction
        if (_d->GetLowerBounds((*t)[0], z->GetIdx() - 1) > 0.01) // we visited this
          continue;

        rAE = _d->GetLowerBounds((*t)[0], z->GetIdx() - 1) + DIST12_TOL;
        rBE = _d->GetLowerBounds((*t)[1], z->GetIdx() - 1) + DIST12_TOL;
        A = Calculate13Angle(rAB, rAE, rBE);

        // default bounds
        lBounds = Calculate15DistAnyCis(rAE, rAB, rBC, rCD, A, B, C);
        uBounds = Calculate15DistAnyTrans(rAE, rAB, rBC, rCD, A, B, C);

        // Check stereochemistry
        if (stereo && stereo->IsCis(z->GetId(), c->GetId()))
          uBounds = lBounds; // Must be cis
        if (stereo && stereo->IsTrans(z->GetId(), c->GetId()))
          lBounds = uBounds; // Must be trans

        // Correcting ring shapes -- should be mostly cisoid
        if (AreInSameRing(d, &*z))
          uBounds = (lBounds + uBounds) / 2.0;

        if (_d->GetLowerBounds(z->GetIdx() - 1, (*t)[3]) < lBounds)
          _d->SetLowerBounds(z->GetIdx() - 1, (*t)[3], lBounds - DIST15_TOL);
        _d->SetUpperBounds(z->GetIdx() - 1, (*t)[3], uBounds + DIST15_TOL);
      }
    }
  }

  int OBDistanceGeometry::AreInSameRing(OBAtom *a, OBAtom *b)
  {
    bool a_in, b_in;
    vector<OBRing*> vr;
    vr = _mol.GetSSSR();

    vector<OBRing*>::iterator i;
    vector<int>::iterator j;

    for (i = vr.begin();i != vr.end();++i) {
      a_in = false;
      b_in = false;
      // Go through the path of the ring and see if a and/or b match
      // each node in the path
      for(j = (*i)->_path.begin();j != (*i)->_path.end();++j) {
        if ((unsigned)(*j) == a->GetIdx())
          a_in = true;
        if ((unsigned)(*j) == b->GetIdx())
          b_in = true;
      }

      if (a_in && b_in)
        return (*i)->Size();
    }

    return 0;
  }

  //! Implements the smoothing described by
  //! Dress, AWM, Havel TF; Discrete Applied Mathematics (1988) v. 19 pp. 129-144
  //! "Shortest Path Problems and Molecular Conformation"
  //! https://doi.org/10.1016/0166-218X(88)90009-1
  void OBDistanceGeometry::TriangleSmooth()
  {
    int a, b, c;

    _d->maxBoxSize = 0.0; // size of surrounding space

    float u_ab, u_bc, u_ac; // upper limits
    float l_ab, l_bc, l_ac; // lower limits
    FOR_ATOMS_OF_MOL (_a, _mol) {
      a = _a->GetIdx() - 1;
      FOR_ATOMS_OF_MOL (_b, _mol) {
        if (&*_b == &*_a)
          continue;
        b = _b->GetIdx() - 1;

        // Get upper and lower bounds for ab
        u_ab = _d->GetUpperBounds(a, b);
        l_ab = _d->GetLowerBounds(a, b);
        FOR_ATOMS_OF_MOL (_c, _mol) {
          if (_c->GetIdx() <= _b->GetIdx())
            continue;
          if (&*_c == &*_a)
            continue;

          c = _c->GetIdx() - 1;

          // get the upper and lower limits for bc and ac
          u_bc = _d->GetUpperBounds(b, c);
          l_bc = _d->GetLowerBounds(b, c);
          u_ac = _d->GetUpperBounds(a, c);
          l_ac = _d->GetLowerBounds(a, c);

          // Triangle rule: length can't be longer than the sum of the two other legs
          //   here "a" is the vertex
          if (u_bc > (u_ab + u_ac)) { // u_bc <= u_ab + u_bc
            u_bc = u_ab + u_ac;
            _d->SetUpperBounds(b, c, u_bc);
          }

          // Triangle rule: length can't be shorter than the difference between the legs
          if (l_bc < (l_ab - l_ac)) {
            l_bc = l_ab - l_ac;
            _d->SetLowerBounds(b, c, l_bc);
          } else if (l_bc < (l_ac - l_ab)) {
            l_bc = l_ac - l_ab;
            _d->SetLowerBounds(b, c, l_bc);
          }

          if (u_bc < l_bc) {
            _d->SetUpperBounds(b, c, l_bc);
            //obErrorLog.ThrowError(__FUNCTION__, "Triagle Smoothing: Erroneous Bounds.", obWarning);
          }

        } // loop(c)

        // Update boxSize after all "c" updates for this pair
        if (_d->GetUpperBounds(a, b) > _d->maxBoxSize)
          _d->maxBoxSize = _d->GetUpperBounds(a, b);
      } // loop(b)
    } // loop(a)
  }

  void OBDistanceGeometry::SetLowerBounds()
  {
    // Ensure atoms aren't closer than VDW contacts
    OBAtom *a, *b;
    unsigned int N = _mol.NumAtoms();
    float aRad, bRad, minDist;

    for (unsigned int i = 0; i < N; ++i) {
      a = _mol.GetAtom(i+1);
      aRad = OBElements::GetVdwRad(a->GetAtomicNum());

      for (unsigned int j = i + 1; j < N; ++j)
        {
          b = _mol.GetAtom(j + 1);
          bRad = OBElements::GetVdwRad(b->GetAtomicNum());
          minDist = aRad + bRad;
          if (minDist < 1.0f)
            minDist = 1.0f;

          if (!AreInSameRing(a, b))
            minDist += 0.1; // prevents bonds going through rings

          if (!_mol.GetBond(a, b)
              && _d->GetLowerBounds(i, j) < 0.4f) { // only check for nonobonded contacts
              _d->SetLowerBounds(i, j, minDist);
          }
        }
    }
  }

  bool OBDistanceGeometry::CheckStereoConstraints()
  {
    return _d->stereoHelper.Check(&_mol);

    /*
    // Check stereo by canonical SMILES
    StereoFrom3D(&_mol, true);
    OBConversion conv;
    conv.SetOutFormat("can");
    std::string predicted_smiles = conv.WriteString(&_mol, true);
    return input_smiles == predicted_smiles;
    */

    // Check all stereo constraints
    // First, gather the known, specified stereochemistry
    // Get TetrahedralStereos and make a vector of corresponding OBStereoUnits
    // Get CisTrans and make a vector of those too
    /*std::vector<OBTetrahedralStereo*> tetra, newtetra;
    std::vector<OBCisTransStereo*> cistrans, newcistrans;
    OBStereoUnitSet ctSunits, tetSunits;
    std::vector<OBGenericData*> vdata = _mol.GetAllData(OBGenericDataType::StereoData);
    OBStereo::Ref atom_id;
    OBStereo::Ref bond_id;
    for (std::vector<OBGenericData*>::iterator data = vdata.begin(); data != vdata.end(); ++data) {
      // If it's cis-trans and specified
      if (((OBStereoBase*)*data)->GetType() == OBStereo::CisTrans) {
        OBCisTransStereo *ct = dynamic_cast<OBCisTransStereo*>(*data);
        if (ct->GetConfig().specified) {
          cistrans.push_back(ct);
          bond_id = _mol.GetBond(_mol.GetAtomById(ct->GetConfig().begin),
                                 _mol.GetAtomById(ct->GetConfig().end))->GetId();
          ctSunits.push_back(OBStereoUnit(OBStereo::CisTrans, bond_id));
        }
      }

      if (((OBStereoBase*)*data)->GetType() == OBStereo::Tetrahedral) {
        OBTetrahedralStereo *th = dynamic_cast<OBTetrahedralStereo*>(*data);
        if (th->GetConfig().specified) {
          tetra.push_back(th);
          atom_id = th->GetConfig().center;
          tetSunits.push_back(OBStereoUnit(OBStereo::Tetrahedral, atom_id));
        }
      } // end tetrahedral
    } // end for (i.e., saving the known, specified stereochemistry

    // We'll check cis/trans first
    newcistrans = CisTransFrom3D(&_mol, ctSunits, false);
    std::vector<OBCisTransStereo*>::iterator origct, newct;
    for (origct=cistrans.begin(), newct=newcistrans.begin(); origct!=cistrans.end(); ++origct, ++newct) {
      if ((*origct)->GetConfig(OBStereo::ShapeU)
          !=  (*newct)->GetConfig(OBStereo::ShapeU)) {
        // Wrong cis/trans stereochemistry
        return false;
      }
    } // end checking cis-trans

    // Perceive TetrahedralStereos from current geometry
    newtetra = TetrahedralFrom3D(&_mol, tetSunits, false);
    // Iterate through original and new stereo and validate
    std::vector<OBTetrahedralStereo*>::iterator origth, newth;
    for (origth=tetra.begin(), newth=newtetra.begin(); origth!=tetra.end(); ++origth, ++newth) {
      if ( (*origth)->GetConfig(OBStereo::Clockwise, OBStereo::ViewFrom)
           != (*newth)->GetConfig(OBStereo::Clockwise, OBStereo::ViewFrom) )
        return false; // found an invalid center
    }

    // everything validated
    return true;*/
  }

  Eigen::MatrixXf OBDistanceGeometry::GetBoundsMatrix()
  {
    Eigen::MatrixXf returnValue;
    if (_d != nullptr)
      returnValue = _d->bounds;
    return returnValue;
  }

  bool OBDistanceGeometry::SetBoundsMatrix(const Eigen::MatrixXf bounds)
  {
    if (_d != nullptr) {
      // Check size of bounds matrix
      _d->bounds = bounds;
      return true;
    } else
      return false;
  }

  bool OBDistanceGeometry::generateInitialCoords(void) {
    // place atoms randomly
    unsigned int N = _mol.NumAtoms();
    // random distance matrix
    Eigen::MatrixXd distMat = Eigen::MatrixXd::Zero(N, N);
    OBRandom generator;
    generator.TimeSeed();
    for (size_t i=0; i<N; ++i) {
      for(size_t j=0; j<i; ++j) {
        double lb = _d->GetLowerBounds(i, j);
        double ub = _d->GetUpperBounds(i, j);
        double v = generator.NextFloat() * (ub - lb) + lb;
        distMat(i, j) = v;
        distMat(j, i) = v;
      }
    }
    // metrix matrix
    // https://github.com/rdkit/rdkit/blob/master/Code/DistGeom/DistGeomUtils.cpp
    Eigen::MatrixXd sqMat(N, N);
    double sumSqD2 = 0.0;
    for (size_t i=0; i<N; i++) {
      for (size_t j=0; j<=i; j++) {
        double d2 = distMat(i, j) * distMat(i, j);
        sqMat(i, j) = d2;
        sqMat(j, i) = d2;
        sumSqD2 += d2;
      }
    }
    sumSqD2 /= (N * N);

    Eigen::VectorXd sqD0i = Eigen::VectorXd::Zero(N);
    for (size_t i = 0; i < N; i++) {
      for (size_t j = 0; j < N; j++) {
        sqD0i(i) += sqMat(i, j);
      }
      sqD0i(i) /= N;
      sqD0i(i) -= sumSqD2;
    }

    Eigen::MatrixXd T(N, N);
    for (size_t i = 0; i < N; i++) {
      for (size_t j = 0; j <= i; j++) {
        double v = 0.5 * (sqD0i(i) + sqD0i(j) - sqMat(i, j));
        T(i, j) = v;
        T(j, i) = v;
      }
    }
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(T);
    Eigen::VectorXd eigVals = es.eigenvalues();
    Eigen::MatrixXd eigVecs = es.eigenvectors();

    for (size_t i = 0; i < N; i++) {
      if(eigVals(i) > 0) eigVals(i) = sqrt(eigVals(i));
      else eigVals(i) *= -1;
    }

    _coord.resize(N * dim);
    for (size_t i = 0; i < N; i++) {
      for (size_t j = 0; j < dim; j++) {
        if(N-1-j >= 0) _coord(i*dim + j) = eigVals(N-1-j) * eigVecs(i, N-1-j);
        else _coord(i*dim + j) = 0;
      }
    }
    Eigen::MatrixXd distMat2(N, N);
    for (size_t i = 0; i < N; i++) {
      for (size_t j = 0; j < N; j++) {
        for(size_t k = 0; k < dim; k++)
          distMat2(i, j) += pow(_coord(i*dim + k)-_coord(j*dim + k), 2.0);
        distMat2(i, j) = sqrt(distMat2(i, j));
      }
    }

//    for(size_t i=0; i<N; i++) {
//      for (size_t j=0; j<N; j++) {
//        double lb = _d->GetLowerBounds(i, j);
//        double ub = _d->GetUpperBounds(i, j);
//        cout << "(" << i << ", " << j << ")" << lb << " < " << distMat2(i, j) << " (" << distMat(i, j) << ")" << " < " << ub << endl;
//      }
//    }
    return true;
  }

  bool OBDistanceGeometry::firstMinimization(void) {
    unsigned int N = _mol.NumAtoms();
    for(size_t i=0; i<N; ++i) {
      vector3 v(_coord(i*dim), _coord(i*dim+1), _coord(i*dim+2));
      OBAtom* a = _mol.GetAtom(i+1);
      a->SetVector(v);
    }
    DistGeomFunc f(this);

    LBFGSpp::LBFGSParam<double> param;
    param.epsilon = 1e-6;
    param.max_iterations = 1000;

    // Create solver and function object
    LBFGSpp::LBFGSSolver<double> solver(param);
    DistGeomFunc fun(this);

//    double fx;
//    int niter = solver.minimize(fun, _coord, fx);
//    std::cout << niter << " iterations" << std::endl;
//    std::cout << "f(x) = " << fx << std::endl;

    for(size_t i=0; i<N; ++i) {
      vector3 v(_coord(i*dim), _coord(i*dim+1), _coord(i*dim+2));
      OBAtom* a = _mol.GetAtom(i+1);
      a->SetVector(v);
    }

    Eigen::MatrixXd distMat2(N, N);
    for (size_t i = 0; i < N; i++) {
      for (size_t j = 0; j < N; j++) {
        for(size_t k = 0; k < dim; k++)
          distMat2(i, j) += pow(_coord(i*dim + k)-_coord(j*dim + k), 2.0);
        distMat2(i, j) = sqrt(distMat2(i, j));
      }
    }

//    for(size_t i=0; i<N; i++) {
//      for (size_t j=0; j<N; j++) {
//        double lb = _d->GetLowerBounds(i, j);
//        double ub = _d->GetUpperBounds(i, j);
//        cout << "(" << i << ", " << j << ")" << lb << " < " << distMat2(i, j) << " < " << ub << endl;
//      }
//    }
    return true;
  }

  bool OBDistanceGeometry::minimizeFourthDimension(void) {
    unsigned int N = _mol.NumAtoms();
    for(size_t i=0; i<N; ++i) {
      vector3 v(_coord(i*dim), _coord(i*dim+1), _coord(i*dim+2));
      OBAtom* a = _mol.GetAtom(i+1);
      a->SetVector(v);
    }

    LBFGSpp::LBFGSParam<double> param;
    param.epsilon = 1e-6;
    param.max_iterations = 2000;

    LBFGSpp::LBFGSSolver<double> solver(param);
    DistGeomFunc4D fun(this);

//    double fx;
//    int niter = solver.minimize(fun, _coord, fx);
//    std::cout << niter << " iterations" << std::endl;
//    std::cout << "f(x) = " << fx << std::endl;

    for(size_t i=0; i<N; ++i) {
      vector3 v(_coord(i*dim), _coord(i*dim+1), _coord(i*dim+2));
      OBAtom* a = _mol.GetAtom(i+1);
      a->SetVector(v);
    }
    return true;
  }

  void OBDistanceGeometry::AddConformer()
  {
    // We should use Eigen here, and cast to double*
    double *confCoord = new double [_mol.NumAtoms() * 3]; // initial state (random)
    _mol.AddConformer(confCoord);
    _mol.SetConformer(_mol.NumConformers());

    OBRandom generator(true); // Use system rand() functions
    generator.TimeSeed();

    if (_d->debug) {
      cerr << " max box size: " << _d->maxBoxSize << endl;
    }

    _d->success = false;
    unsigned int maxIter = 10 * _mol.NumAtoms();
    for (unsigned int trial = 0; trial < maxIter; trial++) {
      generateInitialCoords();
      firstMinimization();
      if (dim == 4) minimizeFourthDimension();
      if (CheckStereoConstraints() && CheckBounds()) {
        _d->success = true;
        break;
      }
      if (_d->debug && !_d->success)
        cerr << "Stereo unsatisfied, trying again" << endl;
    }
  }

  bool OBDistanceGeometry::WasSuccessful() const
  {
    return _d->success;
  }

  bool OBDistanceGeometry::CheckBounds()
  {
    // remember atom indexes from 1
    OBAtom *a, *b;
    double dist, aRad, bRad, minDist, uBounds;

    for (unsigned int i = 1; i <= _mol.NumAtoms(); ++i) {
      a = _mol.GetAtom(i);
      aRad = OBElements::GetVdwRad(a->GetAtomicNum());
      for (unsigned int j = i + 1; j <= _mol.NumAtoms(); ++j) {
          b = _mol.GetAtom(j);

          // Compare the current distance to the lower and upper bounds
          dist = a->GetDistance(b);
          // upper first
          uBounds = _d->GetUpperBounds(i - 1, j - 1);
          if (dist - uBounds > 2.5) {
                if (_d->debug) {
                  cerr << " upper violation " << dist << " " << uBounds << endl;
                }
            return false;
          }
          // now lower.. if the two atoms aren't bonded
          if (_mol.GetBond(a, b))
            continue;

          bRad = OBElements::GetVdwRad(b->GetAtomicNum());
          minDist = aRad + bRad - 2.5;
          if (minDist < 0.8)
            minDist = 0.8;

          // Compare the current distance to the lower bounds
          dist = a->GetDistance(b);
          if (dist < minDist) {
            if (_d->debug) {
                  cerr << " lower violation " << dist << " " << minDist << endl;
            }
            return false;
          }
      }
    }

    return true;
  }

  void OBDistanceGeometry::GetConformers(OBMol &mol)
  {

    // Sanity Check
    if (_mol.NumAtoms() != mol.NumAtoms()) {
      obErrorLog.ThrowError(__FUNCTION__, "The number of atoms did not match.", obWarning);
      return;
    }

    mol.SetDimension(3);

    //Copy conformer information
    if (_mol.NumConformers() > 0) {
      int k,l;
      vector<double*> conf;
      double* xyz = nullptr;
      for (k=0 ; k<_mol.NumConformers() ; ++k) {
        xyz = new double [3*_mol.NumAtoms()];
        for (l=0 ; l<(int) (3*_mol.NumAtoms()) ; ++l)
          xyz[l] = _mol.GetConformer(k)[l];
        conf.push_back(xyz);
      }
      mol.SetConformers(conf);
    }
  }

  bool OBDistanceGeometry::GetGeometry(OBMol &mol, bool useCurrentGeom)
  {
    mol.AddHydrogens();
    if (!Setup(mol, useCurrentGeom))
      return false;

    AddConformer();
    GetConformers(mol);

    return _d->success;
  }

  double DistGeomFunc::operator() (const Eigen::VectorXd& x, Eigen::VectorXd& grad) {
    unsigned int dim = owner->GetDimension();
    const size_t size = x.size()/dim;
    double ret = 0.0;
    // calculate distance error
    for(size_t i=0; i<size; ++i) {
        for(size_t j=0; j<size; ++j) {
            double v = 0.0;
            double d2 = 0;
            for(size_t k=0; k<dim; k++)
              d2 += pow(x[i*dim+k]-x[j*dim+k], 2.0);
            double d = sqrt(d2);
            double ub = owner->GetUpperBounds(i, j);
            double lb = owner->GetLowerBounds(i, j);
            double u2 = ub * ub;
            double l2 = lb * lb;
            if (d > ub) v = d2/u2-1.0; 
            else if (d < lb) v = (2.0*l2 / (l2+d2))-1.0;
            if(v > 0.0) ret += v*v;
        }
    }
    // calculate distance error
    for(size_t i = 0; i < owner->_stereo.size(); i++) {
      TetrahedralInfo tetra = owner->_stereo[i];
      vector<unsigned long> nbrs = tetra.GetNeighbors();
      Eigen::Vector3d v1(x(nbrs[0] * dim), x(nbrs[0]*dim+1), x(nbrs[0]*dim+2));
      Eigen::Vector3d v2(x(nbrs[1] * dim), x(nbrs[1]*dim+1), x(nbrs[1]*dim+2));
      Eigen::Vector3d v3(x(nbrs[2] * dim), x(nbrs[2]*dim+1), x(nbrs[2]*dim+2));
      Eigen::Vector3d v4(x(nbrs[3] * dim), x(nbrs[3]*dim+1), x(nbrs[3]*dim+2));
      double vol = (v1-v4).dot((v2-v4).cross(v3-v4));
      double lb = tetra.GetLowerBound();
      double ub = tetra.GetUpperBound();
      if(vol < lb) ret += (vol - lb) * (vol - lb);
      else if(vol > ub) ret += (vol - ub) * (vol - ub);
    }
    // clear gradient
    for (size_t i=0; i<grad.rows(); i++) {
      grad[i] = 0;
    }
    // gradient for distance error
    unsigned int N = x.size() / dim;
    for(size_t i=0; i<N; ++i) {
      for(size_t j=0; j<N; ++j) {
        double preFactor = 0.0;
        double ub = owner->GetUpperBounds(i, j);
        double lb = owner->GetLowerBounds(i, j);
        double d2 = 0;
        for(size_t k=0; k<dim; k++)
          d2 += pow(x[i*dim+k]-x[j*dim+k], 2.0);
        double d = sqrt(d2);
        if (d > ub) {
          double u2 = ub * ub;
          preFactor = 4.0 * (((d2)/u2) - 1.0) * (d/u2);
        } else if(d < lb) {
          double l2 = lb * lb;
          double l2d2 = d2 + l2;
          preFactor = 8.0 * l2 * d * (1.0 - 2.0 * l2 / l2d2) / (l2d2 * l2d2);
        }
        for (size_t k=0; k<dim; ++k) {
          double g = 0;
          if(d > 0) g = preFactor * (x[i*dim+k] - x[j*dim+k]) / d;
          grad[i * dim + k] += g;
          grad[j * dim + k] += -g;
        }
      }
    }
    // gradient for chiral error
    for(size_t i = 0; i < owner->_stereo.size(); i++) {
      TetrahedralInfo tetra = owner->_stereo[i];
      vector<unsigned long> nbrs = tetra.GetNeighbors();
      unsigned long idx1, idx2, idx3, idx4;
      idx1 = nbrs[0];
      idx2 = nbrs[1];
      idx3 = nbrs[2];
      idx4 = nbrs[3];

      Eigen::Vector3d v1(x(idx1 * dim), x(idx1 * dim+1), x(idx1 * dim+2));
      Eigen::Vector3d v2(x(idx2 * dim), x(idx2 * dim+1), x(idx2 * dim+2));
      Eigen::Vector3d v3(x(idx3 * dim), x(idx3 * dim+1), x(idx3 * dim+2));
      Eigen::Vector3d v4(x(idx4 * dim), x(idx4 * dim+1), x(idx4 * dim+2));
      v1 -= v4;
      v2 -= v4;
      v3 -= v4;

      double vol = v1.dot(v2.cross(v3));
      double lb = tetra.GetLowerBound();
      double ub = tetra.GetUpperBound();

      double preFactor;
      if (vol < lb) preFactor = vol - lb;
      else if (vol > ub) preFactor = vol - ub;
      else continue;

      // RDKit: Code/DistGeom/ChiralViolationContrib.cpp
      grad[dim * idx1] += preFactor * (v2.y() * v3.z() - v3.y() * v2.z());
      grad[dim * idx1 + 1] += preFactor * (v3.x() * v2.z() - v2.x() * v3.z());
      grad[dim * idx1 + 2] += preFactor * (v2.x() * v3.y() - v3.x() * v2.y());

      grad[dim * idx2] += preFactor * (v3.y() * v1.z() - v3.z() * v1.y());
      grad[dim * idx2 + 1] += preFactor * (v3.z() * v1.x() - v3.z() * v1.z());
      grad[dim * idx2 + 2] += preFactor * (v3.x() * v1.y() - v3.y() * v1.x());

      grad[dim * idx3] += preFactor * (v2.z() * v1.y() - v2.y() * v1.z());
      grad[dim * idx3 + 1] += preFactor * (v2.x() * v1.z() - v2.z() * v1.x());
      grad[dim * idx3 + 2] += preFactor * (v2.y() * v1.x() - v2.x() * v1.y());

      grad[dim * idx4] += 
        preFactor *
        (x[idx1 * dim + 2] * (x[idx2 * dim + 1] - x[idx3 * dim + 1]) + 
         x[idx2 * dim + 2] * (x[idx3 * dim + 1] - x[idx1 * dim + 1]) + 
         x[idx3 * dim + 2] * (x[idx1 * dim + 1] - x[idx2 * dim + 1]));

      grad[dim * idx4 + 1] += 
        preFactor *
        (x[idx1 * dim] * (x[idx2 * dim + 2] - x[idx3 * dim + 2]) + 
         x[idx2 * dim] * (x[idx3 * dim + 2] - x[idx1 * dim + 2]) + 
         x[idx3 * dim] * (x[idx1 * dim + 2] - x[idx2 * dim + 2]));

      grad[dim * idx4 + 2] += 
        preFactor *
        (x[idx1 * dim + 1] * (x[idx2 * dim] - x[idx3 * dim]) + 
         x[idx2 * dim + 1] * (x[idx3 * dim] - x[idx1 * dim]) + 
         x[idx3 * dim + 1] * (x[idx1 * dim] - x[idx2 * dim]));
    }
    return ret;
  }

  double DistGeomFunc4D::operator() (const Eigen::VectorXd& x, Eigen::VectorXd& grad) {
    unsigned int dim = owner->GetDimension();
    const size_t size = x.size()/dim;
    double ret = 0.0;
    // calculate distance error
    for(size_t i=0; i<size; ++i) {
        for(size_t j=0; j<size; ++j) {
            double v = 0.0;
            double d2 = 0;
            for(size_t k=0; k<dim; k++)
              d2 += pow(x[i*dim+k]-x[j*dim+k], 2.0);
            double d = sqrt(d2);
            double ub = owner->GetUpperBounds(i, j);
            double lb = owner->GetLowerBounds(i, j);
            double u2 = ub * ub;
            double l2 = lb * lb;
            if (d > ub) v = d2/u2-1.0; 
            else if (d < lb) v = (2.0*l2 / (l2+d2))-1.0;
            if(v > 0.0) ret += v*v;
        }
    }
    // calculate distance error
    for(size_t i = 0; i < owner->_stereo.size(); i++) {
      TetrahedralInfo tetra = owner->_stereo[i];
      vector<unsigned long> nbrs = tetra.GetNeighbors();
      Eigen::Vector3d v1(x(nbrs[0] * dim), x(nbrs[0]*dim+1), x(nbrs[0]*dim+2));
      Eigen::Vector3d v2(x(nbrs[1] * dim), x(nbrs[1]*dim+1), x(nbrs[1]*dim+2));
      Eigen::Vector3d v3(x(nbrs[2] * dim), x(nbrs[2]*dim+1), x(nbrs[2]*dim+2));
      Eigen::Vector3d v4(x(nbrs[3] * dim), x(nbrs[3]*dim+1), x(nbrs[3]*dim+2));
      double vol = (v1-v4).dot((v2-v4).cross(v3-v4));
      double lb = tetra.GetLowerBound();
      double ub = tetra.GetUpperBound();
      if(vol < lb) ret += (vol - lb) * (vol - lb);
      else if(vol > ub) ret += (vol - ub) * (vol - ub);
    }
    // clear gradient
    for (size_t i=0; i<grad.rows(); i++) {
      grad[i] = 0;
    }
    // gradient for distance error
    unsigned int N = x.size() / dim;
    for(size_t i=0; i<N; ++i) {
      for(size_t j=0; j<N; ++j) {
        double preFactor = 0.0;
        double ub = owner->GetUpperBounds(i, j);
        double lb = owner->GetLowerBounds(i, j);
        double d2 = 0;
        for(size_t k=0; k<dim; k++)
          d2 += pow(x[i*dim+k]-x[j*dim+k], 2.0);
        double d = sqrt(d2);
        if (d > ub) {
          double u2 = ub * ub;
          preFactor = 4.0 * (((d2)/u2) - 1.0) * (d/u2);
        } else if(d < lb) {
          double l2 = lb * lb;
          double l2d2 = d2 + l2;
          preFactor = 8.0 * l2 * d * (1.0 - 2.0 * l2 / l2d2) / (l2d2 * l2d2);
        }
        for (size_t k=0; k<dim; ++k) {
          double g = 0;
          if(d > 0) g = preFactor * (x[i*dim+k] - x[j*dim+k]) / d;
          grad[i * dim + k] += g;
          grad[j * dim + k] += -g;
        }
      }
    }
    // gradient for chiral error
    for(size_t i = 0; i < owner->_stereo.size(); i++) {
      TetrahedralInfo tetra = owner->_stereo[i];
      vector<unsigned long> nbrs = tetra.GetNeighbors();
      unsigned long idx1, idx2, idx3, idx4;
      idx1 = nbrs[0];
      idx2 = nbrs[1];
      idx3 = nbrs[2];
      idx4 = nbrs[3];

      Eigen::Vector3d v1(x(idx1 * dim), x(idx1 * dim+1), x(idx1 * dim+2));
      Eigen::Vector3d v2(x(idx2 * dim), x(idx2 * dim+1), x(idx2 * dim+2));
      Eigen::Vector3d v3(x(idx3 * dim), x(idx3 * dim+1), x(idx3 * dim+2));
      Eigen::Vector3d v4(x(idx4 * dim), x(idx4 * dim+1), x(idx4 * dim+2));
      v1 -= v4;
      v2 -= v4;
      v3 -= v4;

      double vol = v1.dot(v2.cross(v3));
      double lb = tetra.GetLowerBound();
      double ub = tetra.GetUpperBound();

      double preFactor;
      if (vol < lb) preFactor = vol - lb;
      else if (vol > ub) preFactor = vol - ub;
      else continue;

      // RDKit: Code/DistGeom/ChiralViolationContrib.cpp
      grad[dim * idx1] += preFactor * (v2.y() * v3.z() - v3.y() * v2.z());
      grad[dim * idx1 + 1] += preFactor * (v3.x() * v2.z() - v2.x() * v3.z());
      grad[dim * idx1 + 2] += preFactor * (v2.x() * v3.y() - v3.x() * v2.y());

      grad[dim * idx2] += preFactor * (v3.y() * v1.z() - v3.z() * v1.y());
      grad[dim * idx2 + 1] += preFactor * (v3.z() * v1.x() - v3.z() * v1.z());
      grad[dim * idx2 + 2] += preFactor * (v3.x() * v1.y() - v3.y() * v1.x());

      grad[dim * idx3] += preFactor * (v2.z() * v1.y() - v2.y() * v1.z());
      grad[dim * idx3 + 1] += preFactor * (v2.x() * v1.z() - v2.z() * v1.x());
      grad[dim * idx3 + 2] += preFactor * (v2.y() * v1.x() - v2.x() * v1.y());

      grad[dim * idx4] += 
        preFactor *
        (x[idx1 * dim + 2] * (x[idx2 * dim + 1] - x[idx3 * dim + 1]) + 
         x[idx2 * dim + 2] * (x[idx3 * dim + 1] - x[idx1 * dim + 1]) + 
         x[idx3 * dim + 2] * (x[idx1 * dim + 1] - x[idx2 * dim + 1]));

      grad[dim * idx4 + 1] += 
        preFactor *
        (x[idx1 * dim] * (x[idx2 * dim + 2] - x[idx3 * dim + 2]) + 
         x[idx2 * dim] * (x[idx3 * dim + 2] - x[idx1 * dim + 2]) + 
         x[idx3 * dim] * (x[idx1 * dim + 2] - x[idx2 * dim + 2]));

      grad[dim * idx4 + 2] += 
        preFactor *
        (x[idx1 * dim + 1] * (x[idx2 * dim] - x[idx3 * dim]) + 
         x[idx2 * dim + 1] * (x[idx3 * dim] - x[idx1 * dim]) + 
         x[idx3 * dim + 1] * (x[idx1 * dim] - x[idx2 * dim]));
    }

    for(size_t i=0; i<size; ++i) {
      ret += pow(x[i*dim+3], 2.0);
    }

    for(size_t i=0; i<N; ++i) {
      grad[i * dim + 3] += 2.0 * x[i * dim + 3];
    }
    return ret;
  }

} // end namespace

#endif

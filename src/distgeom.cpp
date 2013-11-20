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
#include <openbabel/builder.h>

using namespace std;

#define DIST12_TOL   0.01f
#define DIST13_TOL   0.1f
#define DIST14_TOL   0.4f

#pragma warning(disable : 4244) // warning C4244: '=' : conversion from 'double' to 'float', possible loss of data
#pragma warning(disable : 4305) // warning C4305: '*=' : truncation from 'double' to 'float'

namespace OpenBabel {

  class DistanceGeometryPrivate {
  public:
    DistanceGeometryPrivate(const unsigned int N)
    {
      bounds = Eigen::MatrixXf(static_cast<int>(N), static_cast<int>(N));
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

    Eigen::MatrixXf bounds;
  };

  OBDistanceGeometry::OBDistanceGeometry(): _d(NULL) {}

  OBDistanceGeometry::OBDistanceGeometry(const OBMol &mol, bool useCurrentGeometry)
  {
    Setup(mol, useCurrentGeometry);
  }

  OBDistanceGeometry::~OBDistanceGeometry()
  {
    if (_d != NULL)
      delete _d;
  }

  bool OBDistanceGeometry::Setup(const OBMol &mol, bool useCurrentGeometry)
  {
    if (_d != NULL)
      delete _d;
    // TODO: add IsSetupNeeded() like OBForceField to prevent duplication of work

    _mol = mol;
    _mol.SetDimension(3);
    _d = new DistanceGeometryPrivate(mol.NumAtoms());

    SetDefaultBounds();
    // Do we use the current geometry for default 1-2 and 1-3 bounds?
    Set12Bounds(useCurrentGeometry);
    cerr << _d->bounds << endl;
    Set13Bounds(useCurrentGeometry);
    cerr << endl << " 13 Matrix\n";
    cerr << _d->bounds << endl;
    Set14Bounds();
    cerr << endl << " 14 Matrix\n";
    cerr << _d->bounds << endl;
    Set15Bounds();
    TriangleSmooth();
    cerr << endl << " Smoothed Matrix\n";
    cerr << _d->bounds << endl;

    return true;
  }

  // Set the default bounds -- vdW contact to maximum distance
  void OBDistanceGeometry::SetDefaultBounds()
  {
    if (!_d)
      return;

    unsigned int N = _mol.NumAtoms();
    float aRad, bRad, minDist, maxDist = N*1.5f; // if, somehow all atoms are in a linear chain
    // TODO: unit cell -- max dist is 1/2 longest body diagonal
    OBAtom *a, *b;
    for (unsigned int i = 0; i < N; ++i) {
      a = _mol.GetAtom(i+1);
      aRad = etab.GetVdwRad(a->GetAtomicNum());

      // set diagonal to zero
      _d->bounds(i, i) = 0.0f;
      for (unsigned int j = i + 1; j < N; ++j)
        {
          b = _mol.GetAtom(j + 1);
          bRad = etab.GetVdwRad(b->GetAtomicNum());
          minDist = aRad + bRad;
          if (minDist < 1.0f)
            minDist = 1.0f;

          _d->SetLowerBounds(i, j, minDist);
          _d->SetUpperBounds(i, j, maxDist);
        }
    }
  }

  void OBDistanceGeometry::Set12Bounds(bool useGeom)
  {
    double length;
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
        _d->SetLowerBounds(i, j, length - DIST12_TOL*4 );
        _d->SetUpperBounds(i, j, length + DIST12_TOL*4 );
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
    return acos((SQUARE(a) + SQUARE(b) - SQUARE(c)) / (2.0*a*b));
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
      if (b->GetBond(c) != NULL)
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
            if (a->GetHyb() == 2 || ringSize <= 4) {
              theta = 180.0f - (360.0f/float(ringSize));
              theta *= DEG_TO_RAD;
            }
            else if (a->GetHyb() == 3 && ringSize == 5)
              theta = 104.0f * DEG_TO_RAD;
            else
              theta = 109.5f * DEG_TO_RAD;
          }
        else { // not all in the same ring
          switch (a->GetHyb()) {
          case (1):
            theta = 180.0f * DEG_TO_RAD;
            break;
          case (2):
            theta = 120.0f * DEG_TO_RAD;
            break;
          case (3):
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
        _d->SetLowerBounds(i, j, dist - DIST13_TOL*2);
        _d->SetUpperBounds(i, j, dist + DIST13_TOL*2);
      } //end unknown geometry
    }
  }

  // - when atoms i and j are in a 1-4 relationship, the lower distance
  //   limit is calculated using a torsional angle of 0.0. The upper limit
  //   is calculated using a torsion angle of 180.0.
  void OBDistanceGeometry::Set14Bounds()
  {
    float rAB, rBC, rCD;
    float rAC, rBD, B, C;
    OBAtom *a, *b, *c, *d;
    unsigned int i, j;

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

      // We really need some special cases (e.g., for rings, amides, etc.)
      _d->SetLowerBounds((*t)[0], (*t)[3], Calculate14DistCis(rAB, rBC, rCD, B, C)   - DIST14_TOL);
      _d->SetUpperBounds((*t)[0], (*t)[3], Calculate14DistTrans(rAB, rBC, rCD, B, C) + DIST14_TOL);
    }
  }

  void OBDistanceGeometry::Set15Bounds()
  {
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

  void OBDistanceGeometry::TriangleSmooth(int iterations)
  {
    int a, b, c;
    int loopCount = 0;
    bool self_consistent = false;
    while (!self_consistent && loopCount < iterations) {
      self_consistent = true;

      double u_ab, u_bc, u_ac; // upper limits
      double l_ab, l_bc, l_ac; // lower limits
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
            if ((&*_c == &*_b) || (&*_c == &*_a))
              continue;
            c = _c->GetIdx() - 1;

            // get the upper and lower limits for bc and ac
            u_bc = _d->GetUpperBounds(b, c);
            //            l_bc = _d->GetLowerBounds(b, c);
            u_ac = _d->GetUpperBounds(a, c);
            l_ac = _d->GetLowerBounds(a, c);

            if (u_ac > (u_ab + u_bc)) { // u_ac <= u_ab + u_bc
              u_ac = u_ab + u_bc;
              self_consistent = false;
            }

            if (l_ac < (l_ab - u_bc)) {// l_ac >= l_ab - u_bc
              l_ac = l_ab - u_bc;
      	      self_consistent = false;
            }

            // store smoothed l_ac and u_ac
            _d->SetUpperBounds(a, c, u_ac);
            _d->SetLowerBounds(a, c, l_ac);
          } // loop(c)
        } // loop(b)
      } // loop(a)
      loopCount++; // Make sure we don't enter an infinite loop
    } // self-consistency
  }

  bool OBDistanceGeometry::CheckChiralConstraints()
  {
    OBBuilder::CorrectStereoBonds(_mol);
    OBBuilder::CorrectStereoAtoms(_mol);

    return true;
  }

  Eigen::MatrixXf OBDistanceGeometry::GetBoundsMatrix()
  {
    Eigen::MatrixXf returnValue;
    if (_d != NULL)
      returnValue = _d->bounds;
    return returnValue;
  }

  bool OBDistanceGeometry::SetBoundsMatrix(const Eigen::MatrixXf bounds)
  {
    if (_d != NULL) {
      // Check size of bounds matrix
      _d->bounds = bounds;
      return true;
    } else
      return false;
  }

  void OBDistanceGeometry::AddConformer()
  {
    // We should use Eigen here, and cast to double*
    double *confCoord = new double [_mol.NumAtoms() * 3]; // initial state (random)
    _mol.AddConformer(confCoord);
    _mol.SetConformer(_mol.NumConformers());

    OBRandom generator(true); // Use system rand() functions
    generator.TimeSeed();

    // Generate initial positions for the atoms
    // We generate a totally random position for atom 1, within a unit sphere
    vector3 start = _mol.GetAtom(1)->GetVector();
    start.randomUnitVector();
    _mol.GetAtom(1)->SetVector(start);
    // Then we place atoms at random positions along their bounds to atom 1
    unsigned int i,j;
    double lBounds, uBounds, dist;
    FOR_ATOMS_OF_MOL(a, _mol) {
      if (a->GetIdx() == 1)
        continue;
      j = a->GetIdx() - 1;

      lBounds = _d->GetLowerBounds(0, j);
      uBounds = _d->GetUpperBounds(0, j);
      dist = lBounds + (uBounds - lBounds) * generator.NextFloat();
      vector3 newPos;
      newPos.randomUnitVector();
      newPos = newPos*dist;
      a->SetVector(newPos + start);
    }

    // Iterate to ensure all atoms satisfy the bounds matrix
    OBAtom *a, *b;
    double lambda;
    for (unsigned int count = 0; count < 10; ++count) {
      lambda = 1.0 - (0.09)*count; // damp the oscillations each cycle
      // remember atom indexes from 1
      for (i = 1; i <= _mol.NumAtoms(); ++i) {
        a = _mol.GetAtom(i);
        for (j = i + 1; j <= _mol.NumAtoms(); ++j) {
          b = _mol.GetAtom(j);

          // Compare the current distance to the lower and upper bounds
          dist = a->GetDistance(b);
          lBounds = _d->GetLowerBounds(i - 1, j - 1);
          uBounds = _d->GetUpperBounds(i - 1, j - 1);

          if (dist < lBounds) {
            vector3 delta = a->GetVector() - b->GetVector();
            double scale = lambda * 0.5 * (lBounds - dist)/dist;
            delta *= scale;
            a->SetVector(a->GetVector() + delta);
            b->SetVector(b->GetVector() - delta);
            cerr << dist << " " << lBounds << " " << a->GetDistance(b) << endl;
          } else if (dist > uBounds) {
            vector3 delta = a->GetVector() - b->GetVector();
            double scale = lambda * 0.5 * (uBounds - dist)/dist;
            delta *= scale;
            a->SetVector(a->GetVector() + delta);
            b->SetVector(b->GetVector() - delta);
            cerr << dist << " " << lBounds << " " << a->GetDistance(b) << endl;
          }
        }
      }
      CheckChiralConstraints();
    }

    // OK, we'll probably need to correct the chirality
    _mol.Center();
  }

  void OBDistanceGeometry::GetConformers(OBMol &mol)
  {
    // Sanity Check
    if (_mol.NumAtoms() != mol.NumAtoms())
      return;

    //Copy conformer information
    if (_mol.NumConformers() > 0) {
      int k,l;
      vector<double*> conf;
      double* xyz = NULL;
      for (k=0 ; k<_mol.NumConformers() ; ++k) {
        xyz = new double [3*_mol.NumAtoms()];
        for (l=0 ; l<(int) (3*_mol.NumAtoms()) ; ++l)
          xyz[l] = _mol.GetConformer(k)[l];
        conf.push_back(xyz);
      }
      mol.SetConformers(conf);
    }
  }

} // end namespace

#endif

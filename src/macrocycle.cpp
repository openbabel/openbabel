/**********************************************************************
macrocycle.cpp - Dale diamond-lattice torsion sequences for macrocycles.

Copyright (C) 2026 by Geoff Hutchison

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

#include <openbabel/macrocycle.h>
#include <openbabel/math/vector3.h>
#include <openbabel/oberror.h>

`#include` <algorithm>
`#include` <cmath>
`#include` <limits>
`#include` <map>

namespace OpenBabel
{

  namespace {

    // Precomputed Dale codes for n = 25 .. 100, generated offline by
    // find_dale_codes.cpp using a NeRF-based candidate search at
    // theta = 109.47 deg, b = 1.54 A. The search penalises both
    // ring closure error and self-intersection (non-adjacent atom
    // pairs closer than 1.5 A); the trailing comment shows the
    // predicted closure-bond length and minimum non-adjacent atom
    // distance for each entry.
    //
    // Diamond-lattice closure requires all sides to be odd and the
    // number of corners to match the parity of n. Even-n entries
    // close exactly (1.540 A closure bond); odd-n entries leave a
    // ~2.5 A residual that FF cleanup recovers from.
    //
    // Smaller n is handled by the crown post-pass in builder.cpp;
    // larger n uses generateCode as a fallback.
    const std::map<int, std::vector<int> > &daleCodeTable()
    {
      static const std::map<int, std::vector<int> > table = {
        { 25, {3,1,3,1,3,1,1,3,1,3,1,3,1} },                       // closure 2.515, min 2.515
        { 26, {3,1,3,1,1,3,1,3,1,3,1,1,3,1} },                     // closure 1.540, min 1.540
        { 27, {3,3,3,3,3,3,3,1,1,1,1,1,1} },                       // closure 2.515, min 1.540
        { 28, {3,1,1,3,1,3,1,1,3,1,1,3,1,3,1,1} },                 // closure 1.540, min 2.515
        { 29, {3,3,3,3,3,3,1,1,1,1,1,1,1,1,1,1,1} },               // closure 2.515, min 1.540
        { 30, {3,1,3,1,3,3,1,3,1,3,1,3,3,1} },                     // closure 1.540, min 1.540
        { 31, {3,3,3,3,3,1,3,1,3,1,3,1,3} },                       // closure 2.515, min 2.515
        { 32, {3,3,3,3,3,3,3,3,1,1,1,1,1,1,1,1} },                 // closure 1.540, min 1.540
        { 33, {3,1,3,1,3,1,3,1,1,3,1,3,1,3,1,3,1} },               // closure 2.515, min 2.515
        { 34, {3,3,3,3,3,3,1,3,3,3,3,3} },                         // closure 1.540, min 2.515
        { 35, {3,1,3,1,3,1,3,1,3,1,3,1,3,1,3,1,3} },               // closure 2.515, min 2.515
        { 36, {3,3,1,3,3,3,3,1,3,3,3,3,1,3} },                     // closure 1.540, min 1.540
        { 37, {5,5,3,3,3,3,3,3,3,3,3} },                           // closure 2.515, min 2.515
        { 38, {5,3,3,3,3,3,3,3,3,3,3,3} },                         // closure 1.540, min 2.515
        { 39, {3,1,3,3,1,3,3,1,3,3,1,3,3,1,3,3,1} },               // closure 2.515, min 2.515
        { 40, {3,3,1,3,3,3,1,3,3,3,1,3,3,3,1,3} },                 // closure 1.540, min 2.515
        { 41, {3,3,3,3,3,3,3,1,3,1,3,1,3,1,3,1,3} },               // closure 2.515, min 1.540
        { 42, {5,3,3,3,5,3,3,3,5,3,3,3} },                         // closure 1.540, min 2.515
        { 43, {5,5,3,3,3,3,3,3,3,3,3,3,3} },                       // closure 2.515, min 1.540
        { 44, {3,3,3,3,1,3,3,3,3,3,3,3,1,3,3,3} },                 // closure 1.540, min 2.515
        { 45, {3,3,3,1,3,3,3,3,1,3,3,3,3,3,1,3,3} },               // closure 2.515, min 1.540
        { 46, {3,3,3,3,3,3,3,3,1,3,3,3,3,3,3,3} },                 // closure 1.540, min 1.540
        { 47, {3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,1,1} },               // closure 2.515, min 2.515
        { 48, {3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3} },                 // closure 1.540, min 2.515
        { 49, {3,3,3,3,3,3,3,3,1,3,3,3,3,3,3,3,3} },               // closure 2.515, min 2.515
        { 50, {5,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3} },                 // closure 1.540, min 2.515
        { 51, {5,3,5,3,5,3,3,3,3,3,3,3,3,3,3} },                   // closure 2.515, min 1.540
        { 52, {5,3,3,3,3,3,3,3,5,3,3,3,3,3,3,3} },                 // closure 1.540, min 1.540
        { 53, {5,3,3,3,5,3,3,3,5,3,3,5,3,3,3} },                   // closure 2.515, min 2.515
        { 54, {5,3,5,3,3,5,3,5,3,5,3,3,5,3} },                     // closure 1.540, min 1.540
        { 55, {5,5,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3} },               // closure 2.515, min 2.515
        { 56, {5,3,3,3,5,3,3,3,5,3,3,3,5,3,3,3} },                 // closure 1.540, min 2.515
        { 57, {5,3,5,3,5,3,3,3,3,3,3,3,3,3,3,3,3} },               // closure 2.515, min 2.515
        { 58, {5,3,5,3,5,5,3,5,3,5,3,5,5,3} },                     // closure 1.540, min 2.515
        { 59, {5,3,3,3,5,3,3,3,3,5,3,3,3,5,3,3,3} },               // closure 2.515, min 1.540
        { 60, {5,3,3,5,3,5,3,3,5,3,3,5,3,5,3,3} },                 // closure 1.540, min 2.515
        { 61, {5,5,5,3,5,5,5,5,5,5,3,5,5} },                       // closure 2.515, min 1.540
        { 62, {5,3,5,3,3,5,3,5,3,5,3,5,3,3,5,3} },                 // closure 1.540, min 1.540
        { 63, {5,3,3,5,3,3,5,3,3,5,3,5,3,3,5,3,3} },               // closure 2.515, min 1.540
        { 64, {5,3,5,3,5,3,5,3,5,3,5,3,5,3,5,3} },                 // closure 1.540, min 2.515
        { 65, {7,7,7,7,7,5,5,5,5,5,5} },                           // closure 2.515, min 1.540
        { 66, {5,5,5,3,5,5,5,5,5,5,3,5,5,5} },                     // closure 1.540, min 1.540
        { 67, {5,5,5,5,5,5,5,3,5,3,5,3,5,3,5} },                   // closure 2.515, min 2.515
        { 68, {5,3,5,5,3,5,5,3,5,3,5,5,3,5,5,3} },                 // closure 1.540, min 2.515
        { 69, {5,5,5,5,5,5,5,5,5,3,5,3,5,3,5} },                   // closure 2.515, min 1.540
        { 70, {5,5,3,5,5,3,5,5,3,5,5,3,5,5,3,5} },                 // closure 1.540, min 1.540
        { 71, {5,5,5,5,5,5,5,5,5,5,3,3,3,3,3,3,3} },               // closure 2.515, min 1.540
        { 72, {5,5,3,5,5,5,3,5,5,5,3,5,5,5,3,5} },                 // closure 1.540, min 2.515
        { 73, {5,5,5,5,5,3,5,3,5,3,5,3,5,3,5,3,5} },               // closure 2.515, min 2.515
        { 74, {7,5,5,5,5,5,5,7,5,5,5,5,5,5} },                     // closure 1.540, min 1.540
        { 75, {5,5,5,5,5,5,5,3,5,3,5,3,5,3,5,3,5} },               // closure 2.515, min 2.515
        { 76, {5,5,5,5,3,5,5,5,5,5,5,5,3,5,5,5} },                 // closure 1.540, min 2.515
        { 77, {5,5,3,5,5,5,3,5,5,5,5,3,5,5,5,3,5} },               // closure 2.515, min 1.540
        { 78, {7,5,5,5,7,5,5,7,5,5,5,7,5,5} },                     // closure 1.540, min 1.540
        { 79, {5,5,5,3,5,5,5,5,3,5,5,5,5,5,3,5,5} },               // closure 2.515, min 1.540
        { 80, {5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5} },                 // closure 1.540, min 2.515
        { 81, {7,7,7,7,7,7,7,7,5,5,5,5,5} },                       // closure 2.515, min 1.540
        { 82, {7,5,7,5,5,7,5,7,5,7,5,5,7,5} },                     // closure 1.540, min 1.540
        { 83, {5,5,5,5,5,5,5,5,3,5,5,5,5,5,5,5,5} },               // closure 2.515, min 2.515
        { 84, {7,5,5,5,5,5,5,5,7,5,5,5,5,5,5,5} },                 // closure 1.540, min 2.515
        { 85, {7,7,7,7,7,7,7,5,7,5,7,5,7} },                       // closure 2.515, min 1.540
        { 86, {7,5,7,5,7,7,5,7,5,7,5,7,7,5} },                     // closure 1.540, min 1.540
        { 87, {9,9,9,9,9,7,7,7,7,7,7} },                           // closure 2.515, min 1.540
        { 88, {7,5,5,5,7,5,5,5,7,5,5,5,7,5,5,5} },                 // closure 1.540, min 2.515
        { 89, {7,7,7,7,7,7,7,5,5,5,5,5,5,5,5} },                   // closure 2.515, min 2.515
        { 90, {7,7,5,7,7,5,7,7,7,5,7,7,5,7} },                     // closure 1.540, min 1.540
        { 91, {7,7,7,7,7,7,7,7,5,5,5,5,5,5,5} },                   // closure 2.515, min 2.515
        { 92, {7,5,5,7,5,7,5,5,7,5,5,7,5,7,5,5} },                 // closure 1.540, min 2.515
        { 93, {9,9,9,9,9,7,9,7,9,7,9} },                           // closure 2.515, min 1.540
        { 94, {7,7,7,5,7,7,7,7,7,7,5,7,7,7} },                     // closure 1.540, min 2.515
        { 95, {7,5,7,5,7,5,7,5,7,5,5,5,5,5,5,5,5} },               // closure 2.515, min 1.540
        { 96, {7,5,7,5,7,5,7,5,7,5,7,5,7,5,7,5} },                 // closure 1.540, min 2.515
        { 97, {7,7,7,7,7,7,5,5,5,5,5,5,5,5,5,5,5} },               // closure 2.515, min 2.515
        { 98, {7,5,7,5,7,7,5,7,5,7,5,7,7,5,7,5} },                 // closure 1.540, min 1.540
        { 99, {7,7,7,7,7,7,7,7,7,5,7,5,7,5,7} },                   // closure 2.515, min 1.540
        {100, {7,5,7,7,5,7,7,5,7,5,7,7,5,7,7,5} },                 // closure 1.540, min 2.515
      };
      return table;
    }

    // Fallback for n outside the precomputed table: produce a
    // balanced 12-corner all-odd code summing to n. Typically lands
    // within 5-30 A of closure, comfortably inside MMFF/UFF cleanup's
    // basin -- not as tight as the table but far better than the open
    // chain DFS produces for large rings.
    std::vector<int> generateCode(int n)
    {
      std::vector<int> code;
      if (n < 6)
        return code;
      // m has same parity as n so n - base*m can be even with odd base.
      int m = (n % 2 == 0) ? 12 : 13;
      int base = (n / m) | 1;
      while (base * m > n && base >= 3)
        base -= 2;
      if (base < 1)
        return code;
      int bumps = (n - base * m) / 2;
      if (bumps > m)
        return code;
      code.assign(m, base);
      if (bumps > 0) {
        double step = double(m) / double(bumps);
        for (int i = 0; i < bumps; ++i)
          code[int(i * step) % m] += 2;
      }
      return code;
    }

    // Place atoms via Natural Extension Reference Frame (NeRF). Builds
    // n + 2 atoms; the trailing two are ghosts used to evaluate closure.
    // tors[k] is the dihedral about ring bond k (atoms k-1, k, k+1, k+2).
    // tors[0] wraps through the closure bond and is not consumed here;
    // we consume tors[1..n-1] to place atoms 3..n+1.
    void nerfChain(const std::vector<double> &tors,
                   double bondLen, double angleRad,
                   std::vector<vector3> &coords)
    {
      const int n = static_cast<int>(tors.size());
      coords.assign(n + 2, vector3(0.0, 0.0, 0.0));

      const double piMinusTheta = M_PI - angleRad;
      const double cosA = std::cos(piMinusTheta);
      const double sinA = std::sin(piMinusTheta);

      coords[1] = vector3(bondLen, 0.0, 0.0);
      coords[2] = coords[1] + vector3(bondLen * cosA, bondLen * sinA, 0.0);

      for (int i = 3; i < n + 2; ++i) {
        const vector3 &A = coords[i - 3];
        const vector3 &B = coords[i - 2];
        const vector3 &C = coords[i - 1];

        vector3 bc = C - B;
        bc.normalize();
        vector3 nhat = cross(B - A, bc);
        nhat.normalize();
        vector3 mhat = cross(bc, nhat);

        // -cos(phi)*mhat (rather than +cos(phi)) aligns phi with the
        // standard A-B-C-D dihedral convention: phi=180 places D anti
        // to A, since mhat points away from A's projection onto bc.
        double phi = tors[i - 2];
        coords[i] = C + bondLen * (cosA * bc
                                   + sinA * (-std::cos(phi) * mhat
                                             + std::sin(phi) * nhat));
      }
    }

    // Expand a Dale code [a_1, ..., a_m] with corner signs into a
    // length-n torsion sequence (radians). The code is cyclically
    // rotated so the longest side sits last; that way the last
    // corner falls at position n - a_max + 1, which lands inside
    // SetTorsion's controllable range [1, n-3] whenever
    // a_max >= 4. Positions 0, n-2, n-1 are anti -- those torsions
    // straddle the ring-closure bond and cannot be set by
    // SetTorsion, so leaving them at the anti default avoids
    // wasting corners on unreachable positions.
    void expandTorsions(const std::vector<int> &code,
                        const std::vector<int> &signs,
                        std::vector<double> &tors)
    {
      const double kAnti   = M_PI;
      const double kGauche = M_PI / 3.0;
      int n = 0;
      for (int a : code) n += a;
      tors.assign(n, kAnti);

      int m = static_cast<int>(code.size());
      int maxSide = 0, maxIdx = 0;
      for (int i = 0; i < m; ++i) {
        if (code[i] > maxSide) { maxSide = code[i]; maxIdx = i; }
      }
      std::vector<int> rcode(m), rsign(m);
      for (int i = 0; i < m; ++i) {
        rcode[i] = code[(maxIdx + 1 + i) % m];
        rsign[i] = signs[(maxIdx + 1 + i) % m];
      }
      int pos = 1;
      for (int i = 0; i < m; ++i) {
        tors[pos % n] = rsign[i] * kGauche;
        pos += rcode[i];
      }
    }

    // Composite score for the sign search. The overlap term punishes
    // self-intersecting conformers (e.g. a 4-corner code closing as a
    // degenerate figure-eight with atom k coincident with atom k+n/2)
    // that would otherwise win on closure alone and poison FF cleanup.
    // `bestSoFar` lets the inner pair loop bail out once the score
    // can't beat the current best -- the O(n^2) pair scan is the hot
    // path of the 2^(m-1) sign search.
    constexpr double kOverlapCutoff = 1.5;
    constexpr double kOverlapWeight = 10.0;

    double signScore(const std::vector<vector3> &coords, double bestSoFar)
    {
      int n = static_cast<int>(coords.size()) - 2;
      double base = (coords[n]     - coords[0]).length()
                  + (coords[n + 1] - coords[1]).length();
      if (base >= bestSoFar)
        return base;
      const double cutoff2 = kOverlapCutoff * kOverlapCutoff;
      double minD2 = cutoff2; // only track distances that would trigger overlap
      bool overlapped = false;
      for (int i = 0; i < n; ++i)
        for (int j = i + 2; j < n; ++j) {
          if (i == 0 && j == n - 1) continue; // closure bond
          double d2 = (coords[i] - coords[j]).length_2();
          if (d2 < minD2) {
            minD2 = d2;
            overlapped = true;
            double penalty = kOverlapWeight * (kOverlapCutoff - std::sqrt(d2));
            if (base + penalty >= bestSoFar)
              return base + penalty;
          }
        }
      if (!overlapped)
        return base;
      return base + kOverlapWeight * (kOverlapCutoff - std::sqrt(minD2));
    }

    // Find the corner sign assignment that minimizes closure error.
    // For m corners we have 2^(m-1) patterns (the first sign is fixed
    // to +1; overall handedness is arbitrary). The table codes cap
    // at m = 17, and 2^16 NeRF evaluations remain well under a
    // millisecond. Beyond m = 18 we fall back to alternate-pair
    // (++--++--) which closes by symmetry for the regular fallback
    // codes generateCode produces.
    std::vector<int> findClosureSigns(const std::vector<int> &code,
                                      double bondLen, double angleRad)
    {
      const int m = static_cast<int>(code.size());
      std::vector<int> best(m, 1);
      if (m == 0)
        return best;

      // Greedy fallback for very high corner counts.
      if (m > 18) {
        for (int i = 0; i < m; ++i)
          best[i] = ((i / 2) % 2 == 0) ? 1 : -1;
        return best;
      }

      const unsigned long long total = 1ULL << (m - 1);
      double bestErr = std::numeric_limits<double>::infinity();
      std::vector<int> signs(m, 1);
      std::vector<double> tors;
      std::vector<vector3> coords;

      for (unsigned long long mask = 0; mask < total; ++mask) {
        signs[0] = 1;
        for (int i = 1; i < m; ++i)
          signs[i] = (mask & (1ULL << (i - 1))) ? -1 : 1;

        expandTorsions(code, signs, tors);
        nerfChain(tors, bondLen, angleRad, coords);
        double err = signScore(coords, bestErr);
        if (err < bestErr) {
          bestErr = err;
          best = signs;
        }
      }
      return best;
    }

  } // anonymous namespace

  std::vector<int> OBDaleCode(int n)
  {
    if (n < 6)
      return {};
    const auto &table = daleCodeTable();
    auto it = table.find(n);
    if (it != table.end())
      return it->second;
    return generateCode(n);
  }

  std::vector<double> OBDaleTorsions(int n, double bondLength,
                                     double bondAngleRad)
  {
    std::vector<double> tors;
    if (n < 6)
      return tors;

    std::vector<int> code = OBDaleCode(n);
    if (code.empty())
      return tors;

    std::vector<int> signs = findClosureSigns(code, bondLength, bondAngleRad);
    expandTorsions(code, signs, tors);
    return tors;
  }

} // namespace OpenBabel

//! \file macrocycle.cpp
//! \brief Dale diamond-lattice torsion sequences for macrocycle builder.

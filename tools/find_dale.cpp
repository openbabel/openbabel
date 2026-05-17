/* find_dale.cpp - Generator for the Dale macrocycle code table.
 *
 * Developer-only utility (not part of the standard build). Scans a
 * range of ring sizes, enumerates candidate Dale codes with all-odd
 * sides, evaluates each one via the NeRF chain construction at ideal
 * tetrahedral geometry, and emits a const-table snippet that can be
 * pasted into src/macrocycle.cpp's daleCodeTable().
 *
 * Build:
 *   clang++ -std=c++17 \
 *       -I include -I build/include \
 *       -L build/lib -lopenbabel \
 *       -Wl,-rpath,$(pwd)/build/lib \
 *       tools/find_dale.cpp -o find_dale
 *
 * Run:
 *   ./find_dale [nMin] [nMax]      # defaults: 25 100
 *
 * Background: Dale codes [a_1, ..., a_m] describe a cycloalkane
 * conformer as anti runs of length a_i separated by +/- gauche
 * corners. Closure on the diamond lattice (at theta = 109.4712 deg)
 * requires all sides to be odd lengths so each anti run switches
 * sublattice, and the corner count m must match the parity of n.
 *
 * Scoring penalises three failure modes:
 *   1. Closure-bond stretch (|c[n-1] - c[0]| - 1.54 A).
 *   2. Tangent mismatch at closure (|c[n+1] - c[1]|).
 *   3. Self-intersection (non-adjacent atoms < 1.5 A apart, which
 *      would cause atom k to fold onto atom k+n/2 in degenerate
 *      figure-eight conformers).
 */

#include <openbabel/math/vector3.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

using namespace OpenBabel;

static void nerfChain(const std::vector<double> &tors, double b, double theta,
                      std::vector<vector3> &out)
{
  const int n = static_cast<int>(tors.size());
  out.assign(n + 2, vector3(0, 0, 0));
  const double piMinusTheta = M_PI - theta;
  const double cA = std::cos(piMinusTheta);
  const double sA = std::sin(piMinusTheta);
  out[1] = vector3(b, 0, 0);
  out[2] = out[1] + vector3(b * cA, b * sA, 0);
  for (int i = 3; i < n + 2; ++i) {
    vector3 bc = out[i - 1] - out[i - 2]; bc.normalize();
    vector3 nh = cross(out[i - 2] - out[i - 3], bc); nh.normalize();
    vector3 mh = cross(bc, nh);
    double phi = tors[i - 2];
    out[i] = out[i - 1] + b * (cA * bc
                               + sA * (-std::cos(phi) * mh
                                       + std::sin(phi) * nh));
  }
}

// Corners at offsets 1, 1+a_0, 1+a_0+a_1, ... mod n; all other
// positions are anti (180 deg).
static std::vector<double> expandTorsions(const std::vector<int> &code,
                                          const std::vector<int> &signs,
                                          int n)
{
  std::vector<double> t(n, M_PI);
  int pos = 1;
  for (size_t i = 0; i < code.size(); ++i) {
    t[pos % n] = signs[i] * M_PI / 3.0;
    pos += code[i];
  }
  return t;
}

struct Score {
  double longestBond;  // |c[n-1] - c[0]|   (closure bond length)
  double tangentErr;   // |c[n+1] - c[1]|   (tangent mismatch)
  double minDist;      // min non-adjacent atom distance
};

static double compositeScore(const Score &s)
{
  double overlap = (s.minDist < 1.5) ? 10.0 * (1.5 - s.minDist) : 0.0;
  return std::fabs(s.longestBond - 1.54) + s.tangentErr + overlap;
}

// Exhaustive sign-pattern search over 2^(m-1) patterns (s_0 fixed to
// +1; overall chirality is arbitrary). For m > 17 we fall back to the
// alternate-pair pattern (++--++--), which closes by symmetry for the
// regular candidates generated below.
static Score evaluate(const std::vector<int> &code, int n, double b,
                      double theta)
{
  int m = static_cast<int>(code.size());
  Score best{1e9, 1e9, 1e9};
  int maxMask = (m <= 17) ? (1 << (m - 1)) : 1;
  for (int mask = 0; mask < maxMask; ++mask) {
    std::vector<int> s(m, 1);
    if (m <= 17) {
      for (int i = 1; i < m; ++i)
        s[i] = (mask & (1 << (i - 1))) ? -1 : 1;
    } else {
      for (int i = 0; i < m; ++i)
        s[i] = ((i / 2) % 2 == 0) ? 1 : -1;
    }
    auto t = expandTorsions(code, s, n);
    std::vector<vector3> c;
    nerfChain(t, b, theta, c);
    Score sc;
    sc.longestBond = (c[n - 1] - c[0]).length();
    sc.tangentErr  = (c[n + 1] - c[1]).length();
    sc.minDist = 1e9;
    for (int i = 0; i < n; ++i)
      for (int j = i + 2; j < n; ++j) {
        if (i == 0 && j == n - 1) continue; // intentional closure bond
        double d = (c[i] - c[j]).length();
        if (d < sc.minDist) sc.minDist = d;
      }
    if (compositeScore(sc) < compositeScore(best))
      best = sc;
  }
  return best;
}

// Balanced all-odd code candidates for given n and m. Sides are
// either `base` or `base+2` (both odd); bump placements explored:
// contiguous, evenly spaced, and alternated. C(m, bumps) full
// permutation is overkill -- the symmetric arrangements dominate.
static std::vector<std::vector<int>> candidateCodes(int n, int m)
{
  std::vector<std::vector<int>> result;
  if (m < 3 || m > n) return result;
  int base = (n / m) | 1;
  if (base * m > n) base -= 2;
  if (base < 1) base = 1;
  int extra = n - base * m;
  if (extra < 0 || extra % 2 != 0) return result;
  int bumps = extra / 2;
  if (bumps > m) return result;
  if (bumps == 0) {
    result.push_back(std::vector<int>(m, base));
    return result;
  }
  // Bumps clustered at the front.
  {
    std::vector<int> code(m, base);
    for (int i = 0; i < bumps; ++i) code[i] += 2;
    result.push_back(code);
  }
  // Bumps evenly spaced (round-robin).
  {
    std::vector<int> code(m, base);
    double step = static_cast<double>(m) / static_cast<double>(bumps);
    for (int i = 0; i < bumps; ++i)
      code[static_cast<int>(std::round(i * step)) % m] += 2;
    result.push_back(code);
  }
  // Alternated longest+shortest.
  if (bumps < m) {
    std::vector<int> code(m, base);
    int idx = 0;
    for (int i = 0; i < bumps && idx < m; ++i) {
      code[idx] += 2;
      idx += 2;
      if (idx >= m) idx = 1;
    }
    result.push_back(code);
  }
  return result;
}

int main(int argc, char **argv)
{
  int nMin = 25, nMax = 100;
  if (argc > 1) nMin = std::atoi(argv[1]);
  if (argc > 2) nMax = std::atoi(argv[2]);

  const double b = 1.54;
  const double theta = std::acos(-1.0 / 3.0); // 109.4712 deg

  printf("// Auto-generated by tools/find_dale.cpp.\n");
  printf("// theta = 109.4712 deg, b = 1.54 A. The trailing comment\n");
  printf("// shows the NeRF-predicted longest bond (closure bond) and\n");
  printf("// the minimum non-adjacent atom distance; min < 1.5 A would\n");
  printf("// indicate a self-intersecting chain and is penalised by\n");
  printf("// the search.\n\n");

  for (int n = nMin; n <= nMax; ++n) {
    std::vector<int> bestCode;
    Score bestScore{1e9, 1e9, 0.0};
    // Even n requires even m; odd n requires odd m.
    int mLo = (n % 2 == 0) ? 4 : 3;
    int mStep = 2;
    for (int m = mLo; m <= std::min(n - 1, 20); m += mStep) {
      auto cands = candidateCodes(n, m);
      for (auto &c : cands) {
        Score sc = evaluate(c, n, b, theta);
        if (compositeScore(sc) < compositeScore(bestScore)) {
          bestScore = sc;
          bestCode = c;
        }
      }
    }
    printf("  { %3d, {", n);
    for (size_t i = 0; i < bestCode.size(); ++i)
      printf("%s%d", i ? "," : "", bestCode[i]);
    printf("} },  // closure bond %.3f A, min atom dist %.3f A\n",
           bestScore.longestBond, bestScore.minDist);
  }
  return 0;
}

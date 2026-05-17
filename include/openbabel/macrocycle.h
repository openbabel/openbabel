/**********************************************************************
macrocycle.h - Dale diamond-lattice torsion sequences for macrocycles.

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

#ifndef OB_MACROCYCLE_H
#define OB_MACROCYCLE_H

#include <openbabel/babelconfig.h>

#include <vector>

namespace OpenBabel
{

  /**
   * Return a Dale code [a_1, ..., a_m] (sum a_i = n) describing the
   * conformer of a cycloalkane of size n. For 6 <= n <= 30 the table
   * is sourced from Dale (1973, 1976) and follow-up work; for n > 30
   * a low-corner-count code is generated. Returns an empty vector
   * when n < 6.
   *
   * A Dale code lists the number of bonds between gauche corners
   * (\f$\pm 60^\circ\f$) around the ring. Bonds between corners are
   * trans (180\f$^\circ\f$). Atoms sit on the diamond lattice for
   * ideal tetrahedral geometry.
   */
  OBAPI std::vector<int> OBDaleCode(int n);

  /**
   * Build a length-n torsion sequence (radians) for a cycloalkane
   * backbone using Dale's diamond-lattice scheme. The corner sign
   * pattern is chosen exhaustively to minimize the ring-closure
   * residual at the supplied bond length and bond angle.
   *
   * The returned array is indexed around the ring: torsions[k] is
   * the dihedral about the bond between ring atoms k and k+1
   * (i.e., the dihedral atom (k-1) -- atom k -- atom (k+1) --
   * atom (k+2), indices modulo n). Returns empty when n < 6.
   *
   * Residual closure error after applying these torsions is small
   * enough that a short force-field minimization closes the ring.
   */
  // Default angle is ideal tetrahedral (acos(-1/3) ~ 109.4712 deg),
  // matching the geometry OBBuilder's DFS placement uses for sp3
  // carbons. The precomputed Dale code table is calibrated for this
  // angle; passing a substantially different angle voids the closure
  // guarantees.
  OBAPI std::vector<double> OBDaleTorsions(int n,
                                           double bondLength = 1.54,
                                           double bondAngleRad = 1.9106332);

} // namespace OpenBabel

#endif // OB_MACROCYCLE_H

//! \file macrocycle.h
//! \brief Dale diamond-lattice torsion sequences for macrocycle builder.

/**********************************************************************
canon.h - Canonical labeling.

  Copyright (C) 2009-2010 by Tim Vandermeersch
  Copyright (C) 2005-2006, eMolecules, Inc. (www.emolecules.com)
                           Craig A. James

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

#include <openbabel/bitvec.h>

namespace OpenBabel {

  class OBMol;

  /**
   * Calculate the canonical labels for the molecule. Stereochemistry is
   * included in the algorithm and the canonical labels. The result will be
   * stored in @p canonical_labels.
   *
   * @param mol The molecule.
   * @param symmetry_classes The symmetry_classes for the molecule. These can
   * be obtained using the OBGraphSym class.
   * @param canonical_labels Reference to the object to store the results in.
   * @param mask The fragment to label. When the bit for an atom is set, it is
   * included in the fragment. If no bits are set, all atoms will be included.
   * Atoms are indexed from 1 (i.e. OBAtom::GetIdx()).
   * @param maxSeconds Timeout in seconds.
   * @param onlyOne If true, the first found labels are returned. These are
   * canonical labels without considering stereochemistry and other attributes
   * not included in the symmetry classes.
   *
   * @return The canonical labels for the molecule in @p canonical_labels.
   *
   * @see @ref canonical_code_algorithm
   * @since 2.3
   */
  void OBAPI CanonicalLabels(OBMol *mol, const std::vector<unsigned int> &symmetry_classes,
      std::vector<unsigned int> &canonical_labels, const OBBitVec &mask = OBBitVec(),
      int maxSeconds = 5, bool onlyOne = false);

} // namespace OpenBabel

//! \file canon.h
//! \brief Canonical labeling.

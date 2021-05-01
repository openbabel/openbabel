/**********************************************************************
graphsym.h - Class for handling graph symmetry.

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

#ifndef OB_GRAPHSYM_H
#define OB_GRAPHSYM_H

#include <openbabel/babelconfig.h>
#include <openbabel/stereo/stereo.h>
#include <vector>

#ifndef OB_EXTERN
#  define OB_EXTERN extern
#endif

namespace OpenBabel {

  class OBBitVec;
  class OBMol;
  class OBAtom;
  class OBBond;
  class OBMol;
  class OBGraphSymPrivate;

  /**
   * @class OBGraphSym graphsym.h <openbabel/graphsym.h>
   * @brief Handle and perceive graph symmtery for canonical numbering .
   * @since version 2.3
   */
  class OBAPI OBGraphSym {

    public:
      //! Constructor
      OBGraphSym(OBMol* pmol, const OBBitVec* frag_atoms = nullptr);
      //! Destructor
      virtual ~OBGraphSym();

      static const unsigned int NoSymmetryClass;

      /**
       * Calculate the symmetry classes for the molecule. The result will be
       * stored in @p symmetry_classes.
       *
       * The results in @p symmetry_classes will be ordered by symmetry
       * classes. Use the OBAtom* pointer in the std::pair to match the atoms
       * with the right symmetry classes.
       *
       * @return The number of symmetry classes.
       */
      int GetSymmetry(std::vector<unsigned int> &symmetry_classes);
      /**
       * Clear the symmetry classes data stored in the molecule specified when
       * construting the OBGraphSym object.
       */
      void ClearSymmetry();

    private:
      OBGraphSymPrivate * const d;
  };

} // namespace OpenBabel

#endif // OB_GRAPHSYM_H

//! \file graphsym.h
//! \brief Handle and perceive graph symmtery for canonical numbering

/**********************************************************************
pointgroup.h - Brute force point group symmetry detection

Copyright  (C) 1996, 2003 S. Patchkovskii, Serguei.Patchkovskii@sympatico.ca
Some portions Copyright (C) 2007 by Geoffrey R. Hutchison

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

#ifndef OB_POINTGROUP_H
#define OB_POINTGROUP_H

#include <openbabel/babelconfig.h>

namespace OpenBabel
{
  class OBMol;
  class PointGroupPrivate;

  /** \class OBPointGroup pointgroup.h <openbabel/pointgroup.h>
      \brief Brute-force point group symmetry perception
      \since version 2.2

      This class performs brute-force point group symmetry perception
      to yield symmetry identifiers. In the future, this should provide
      optimized symmetry-derived coordinates as well.
  */
  class OBAPI OBPointGroup
    {
    public:
      OBPointGroup();
      ~OBPointGroup();

      enum Symbol
      { // The most likely 60 cases
        C1 = 0, Cs, Ci, // 0 to 2
        C2, C3, C4, C5, C6, C7, C8, // 3 to 9
        D2, D3, D4, D5, D6, D7, D8, // 10 to 16
        C2v, C3v, C4v, C5v, C6v, C7v, C8v, // 17 to 23
        C2h, C3h, C4h, C5h, C6h, C7h, C8h, // 24 to 30
        D2d, D3d, D4d, D5d, D6d, D7d, D8d, // 31 to 37
        D2h, D3h, D4h, D5h, D6h, D7h, D8h, // 38 to 44
        S4, S6, S8, // 45 to 47
        T, Th, Td, // 48 to 50
        O, Oh, // 51, 52
        Cinfv, Dinfh, // 53, 54
        I, Ih, // 55, 56
        K, Kh, // 57, 58
        Unknown // 59
      };

      //! Set the point group detection for this molecule
      void Setup(OBMol *);

      /** \return the 3D point group symbol for this molecule. A tolerance of
        *  0.01 is used.
        *
        * \todo Compatibility function; remove at next ABI break, update
        *  default arg in overload.
        */
      const char * IdentifyPointGroup();

      //! \return the 3D point group symbol for this molecule
      const char * IdentifyPointGroup(double tolerance /* = 0.01*/ );

      //! \return the 3D point group symbol for this molecule
      //! \arg tolerance The initial tolerance for determining possibly symmetric atoms
      Symbol IdentifyPointGroupSymbol(double tolerance = 0.01);

      void Symmetrize(OBMol *);

    protected:
      PointGroupPrivate *d;

    }; // class OBPointGroup

}// namespace OpenBabel

#endif   // OB_POINT_GROUP_H

//! \file pointgroup.h
//! \brief Brute-force point group detection

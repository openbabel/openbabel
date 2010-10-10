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

#include <string>

#include <openbabel/mol.h>

namespace OpenBabel
{

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

      //! Set the point group detection for this molecule
      void Setup(OBMol *);

      //! \return the 3D point group symbol for this molecule
      const char *IdentifyPointGroup();

    protected:
      PointGroupPrivate *d;

    }; // class OBPointGroup

}// namespace OpenBabel

#endif   // OB_POINT_GROUP_H

//! \file pointgroup.h
//! \brief Brute-force point group detection

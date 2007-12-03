/**********************************************************************
pointgroup.h - Brute force point group symmetry detection
 
Copyright  (C) 1996, 2003 S. Patchkovskii, Serguei.Patchkovskii@sympatico.ca
Some portions Copyright (C) 2007 by Geoffrey R. Hutchison
 
This file is part of the Open Babel project.
For more information, see <http://openbabel.sourceforge.net/>
 
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

  class OBAPI OBPointGroup
    {
    public:
      OBPointGroup();
      ~OBPointGroup();

      void Setup(OBMol *);
      
      char *IdentifyPointGroup();

      /*! Update coordinates based on perceived symmetry
       *  \param mol the OBMol object to copy the coordinates to
       *  \return true if succesful
       */
      bool UpdateCoordinates(OBMol &mol);
      

    protected:
      PointGroupPrivate *d;

    }; // class OBPointGroup

}// namespace OpenBabel

#endif   // OB_POINT_GROUP_H

//! \file pointgroup.h
//! \brief Brute-force point group detection

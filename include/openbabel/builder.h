/**********************************************************************
builder.h - OBBuilder class.
 
Copyright (C) 2007 by Tim Vandermeersch <tim.vandermeersch@gmail.com>
 
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

#ifndef OB_BUILDER_H
#define OB_BUILDER_H

#include <vector>
#include <string>
#include <map>

#include <list>
#include <set>
#include <openbabel/base.h>
#include <openbabel/mol.h>

namespace OpenBabel
{
  //! \class OBBuilder
  //! \brief Class to build 3D structures
  class OBBuilder {
    public:
      //! constructor
      OBBuilder();
      //! destructor
      ~OBBuilder();
      //! add a fragment to mol, a bond is formed between mol atom with index mi and fragment atom with index sf
      bool AddFragment(OBMol &mol, int mi, OBMol &fragment, int fi);
      vector3 GetNewBondVector(OBMol &mol, OBAtom *atom);
      int GetQuadrant(double x, double y);
      double GetRotationAngle(double x, double y, double u, double v);
      bool Build(OBMol &mol);
    private:
      //! used to hold the fragments loaded in the constructor
      std::vector<std::pair<OBSmartsPattern*, std::vector<vector3> > > _fragments;
  }; // class OBBuilder

}// namespace OpenBabel

#endif   // OB_BUILDER_H

//! \file builder.h
//! \brief Class to build 3D structures

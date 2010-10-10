/**********************************************************************
  tetraplanar.cpp - OBTetraPlanarStereo

  Copyright (C) 2009 by Tim Vandermeersch

  This file is part of the Open Babel project.
  For more information, see <http://openbabel.org/>

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
  02110-1301, USA.
 **********************************************************************/
#include <openbabel/stereo/tetraplanar.h>
#include <cassert>

namespace OpenBabel {

  OBTetraPlanarStereo::OBTetraPlanarStereo(OBMol *mol) : OBStereoBase(mol)
  {
  }

  OBTetraPlanarStereo::~OBTetraPlanarStereo()
  {
  }

  /*
  std::vector<unsigned long> OBTetraPlanarStereo::ToInternal(const std::vector<unsigned long> &refs,
      OBStereo::Shape shape)
  {
    assert( refs.size() == 4 );
    std::vector<unsigned long> result(refs);

    switch (shape) {
      case OBStereo::ShapeU:
        // same as internal, just copy
        return result;
      case OBStereo::ShapeZ:
        // normalize to U shape
        result[1] = refs.at(2);
        result[2] = refs.at(3);
        result[3] = refs.at(1);
        return result;
      case OBStereo::Shape4:
        // normalize to U shape
        result[1] = refs.at(2);
        result[2] = refs.at(1);
        return result;
    }

  }

  std::vector<unsigned long> OBTetraPlanarStereo::ToShape(const std::vector<unsigned long> &refs,
      OBStereo::Shape shape)
  {
    assert( refs.size() == 4 );
    std::vector<unsigned long> result(refs);

    switch (shape) {
      case OBStereo::ShapeU:
        // same as internal, just copy
        return result;
      case OBStereo::ShapeZ:
        // convert to U shape
        result[1] = refs.at(3);
        result[2] = refs.at(1);
        result[3] = refs.at(2);
        return result;
      case OBStereo::Shape4:
        // normalize to U shape
        result[1] = refs.at(2);
        result[2] = refs.at(1);
        return result;
    }

  }
  */

}


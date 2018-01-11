/**********************************************************************
  stereo.cpp - OBStereo

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
#include <openbabel/stereo/stereo.h>

namespace OpenBabel {

  bool OBStereo::ContainsSameRefs(const OBStereo::Refs &refs1, const OBStereo::Refs &refs2)
  {
    if (refs1.size() != refs2.size())
      return false;

    unsigned int count = 0;
    for (ConstRefIter i = refs1.begin(); i != refs1.end(); ++i)
      for (ConstRefIter j = refs2.begin(); j != refs2.end(); ++j)
        if (*i == *j) {
          count++;
          break;
        }

    return (count == refs1.size());
  }

  bool OBStereo::ContainsRef(const OBStereo::Refs &refs, unsigned long id)
  {
    for (ConstRefIter i = refs.begin(); i != refs.end(); ++i)
      if (*i == id)
        return true;

    return false;
  }

  int OBStereo::NumInversions(const OBStereo::Refs &refs)
  {
    OBStereo::Refs invVec; // the inversion vector
    OBStereo::ConstRefIter i, j;
    for (i = refs.begin(); i != refs.end(); ++i) {
      int e = 0; // ith element
      // loop over elements to the right
      for (j = i; j != refs.end(); ++j)
        // increment e if element to the right is lower
        if (*j < *i)
          e++;

      invVec.push_back(e);
    }

    int sum = 0;
    for (OBStereo::RefIter k = invVec.begin(); k != invVec.end(); ++k)
      sum += *k;

    return sum;
  }

  void OBStereo::Permutate(OBStereo::Refs &refs, unsigned int i, unsigned int j)
  {
    if (i >= refs.size())
      return;
    if (j >= refs.size())
      return;
    unsigned long id = refs.at(i);
    refs[i] = refs.at(j);
    refs[j] = id;
  }

  OBStereo::Refs OBStereo::Permutated(const OBStereo::Refs &refs, unsigned int i, unsigned int j)
  {
    if (i >= refs.size())
      return refs;
    if (j >= refs.size())
      return refs;
    OBStereo::Refs result(refs);
    result[i] = refs.at(j);
    result[j] = refs.at(i);
    return result;
  }

}


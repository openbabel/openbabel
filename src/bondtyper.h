/**********************************************************************
Copyright (c) 2003 by Geoffrey R. Hutchison

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

#ifndef OB_BONDTYPER_H
#define OB_BONDTYPER_H

#include "parsmart.h"
#include "data.h"

namespace OpenBabel {

  // class introduction in bondtyper.cpp
  // Used for "perceiving" bonds, e.g. in XYZ or QM files with no bond info.
class OBBondTyper : public OBGlobalDataBase
{
  //! SMARTS patterns for functional group typing
  std::vector<std::pair<OBSmartsPattern*, vector<int> >		_fgbonds;
 public:
  OBBondTyper();
  ~OBBondTyper();

  void ParseLine(const char*);
  void ConnectTheDots(OBMol&);
  void PerceiveBondOrders(OBMol&);
};

}

#endif // OB_BONDTYPER_H

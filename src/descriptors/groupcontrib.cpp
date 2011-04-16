/**********************************************************************
groupcontrib.cpp - Handle logP, PSA, MR, and other group-based predictions

Copyright (C) 2007      by Tim Vandermeersch
              2001-2007 by Stephen Jelfs
              2001-2007 by Joerg Kurt Wegner, me@cheminformatics.eu
              2007      by Chris Morley

Original version: JOELib2, http://joelib.sf.net

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

#include <openbabel/babelconfig.h>
#include <openbabel/groupcontrib.h>

using namespace std;

namespace OpenBabel
{

  //******************************************************
  // Make global instances for descriptors which are all calculated
  // from group contibutions in the same way but with different data.

  // LogP (octanol/water partition coefficient)
  OBGroupContrib thelogP("logP", "logp.txt",
    "octanol/water partition coefficient");

  // TPSA (topological polar surface area)
  OBGroupContrib theTPSA("TPSA", "psa.txt",
    "topological polar surface area");

  // MR (molar refractivity)
  OBGroupContrib theMR("MR", "mr.txt",
    "molar refractivity");

//! \file groupcontrib.cpp
//! \brief Handle logP, PSA and other group-based prediction algorithms.

}//namespace

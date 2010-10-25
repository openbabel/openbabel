/**********************************************************************
obmmff94validate.cpp - Validate the MMFF94 implementation

  This program assumes that the following files are in the current
  working directory:

    MMFF94_dative.mol2
    MMFF94_opti.log

  Both these files are from the MMFF94 validation suite wich can be
  found here: http://server.ccl.net/cca/data/MMFF94/

  For Each molecule in the .mol2 file the energy is calculated using
  the openbabel MMFF94 implementation. The energy values of indivudual 
  energy terms are also calculated (total bond energy, totral angle 
  energy, ...). These resulults are compared with values from the
  MMFF94_opti.log file. The difference is geven and the relative error.
  The atomtypes are also checked, when succesfull, PASSED will be printed.
  When an atom has a wrong type, XXXX FAILED XXXX will be printed next
  to its type.

Copyright (C) 2006 Tim Vandermeersch
 
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

// used to set import/export for Cygwin DLLs
#ifdef WIN32
#define USING_OBDLL
#endif

#include <openbabel/babelconfig.h>
#include <openbabel/base.h>
#include <openbabel/mol.h>
#include <openbabel/forcefield.h>
#include <unistd.h>

using namespace std;
using namespace OpenBabel;

int main(int argc,char **argv)
{
  OBForceField* pFF = OBForceField::FindForceField("MMFF94");
  pFF->Validate();

  // test passed
  return 0;
}

/**********************************************************************
chargemodel.cpp - Base class for partial charge models

Copyright (C) 2010 Geoffrey Hutchison
Some portions Copyright (C) 2009 by Frank Peters

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
#include <openbabel/chargemodel.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/obiter.h>

using namespace std;
namespace OpenBabel
{
#if defined(__CYGWIN__) || defined(__MINGW32__)
  // macro to implement static OBPlugin::PluginMapType& Map()
  PLUGIN_CPP_FILE(OBChargeModel)
#endif

  void OBChargeModel::FillChargeVectors(OBMol &mol)
  {
    OBAtom *atom;
    vector<OBAtom*>::iterator itr;
    m_partialCharges.clear();
    m_partialCharges.reserve(mol.NumAtoms());
    m_formalCharges.clear();
    m_formalCharges.reserve(mol.NumAtoms());

    for (atom = mol.BeginAtom(itr);atom;atom = mol.NextAtom(itr))
      {
        m_partialCharges.push_back(atom->GetPartialCharge());
        m_formalCharges.push_back((double)(atom->GetFormalCharge()));
      }
  }

  vector3 OBChargeModel::GetDipoleMoment(OBMol &mol)
  {
    vector3 dipoleMoment = VZero;

    if (this->ComputeCharges(mol)) {
      FOR_ATOMS_OF_MOL(a, mol) {
        dipoleMoment += a->GetVector() * a->GetPartialCharge();
      }
    }

    // adjust by scaling factor
    return dipoleMoment * this->DipoleScalingFactor();
  }

}

//! \file chargemodel.cpp
//! \brief Base class for molecular partial charge models

/**********************************************************************
chargemodel.h - Base class for partial charge models

Copyright (C) 2010 by Geoffrey Hutchison
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

#ifndef OB_CHARGEMODEL_H
#define OB_CHARGEMODEL_H

#include <openbabel/babelconfig.h>
#include <openbabel/plugin.h>
#include <openbabel/math/vector3.h>

namespace OpenBabel
{
class OBMol; //Forward declaration; used only as pointer.

// Documentation is down below
class OBAPI OBChargeModel : public OBPlugin
{
  MAKE_PLUGIN(OBChargeModel)

  public:
    const char* TypeID(){return "charges";};

    /// \return whether partial charges were successfully assigned to this molecule
    /// \note The method should fill m_partialCharges and m_formalCharges as well
    virtual bool ComputeCharges(OBMol &m ) { return false; };
    virtual bool ComputeCharges(OBMol &m, const char *args) { return ComputeCharges( m ); } 

    /// \return a vector of the formal charges on each atom, indexed from 0
    /// This method returns floating point formal charges since some
    /// charge models consider fractional charges (e.g., 0.5 for an
    /// oxygen in a carboxylate CO2- group).
    /// \note If ComputeCharges() has not been called, this will return an empty vector
    const std::vector<double> & GetFormalCharges() const
    { return m_formalCharges; }

    /// \return a vector of the partial charges on each atom, indexed from 0
    /// \note If ComputeCharges() has not been called, this will return an empty vector
    const std::vector<double> & GetPartialCharges() const
    { return m_partialCharges; }

    /// \return a vector of the dipole moment from this molecule
    vector3 GetDipoleMoment(OBMol &);

 protected:
      std::vector<double> m_partialCharges;
      std::vector<double> m_formalCharges;

      /// Fill the internal partial and formal charge vectors (convenience function)
      void FillChargeVectors(OBMol &mol);

      /// Provide a scaling factor for the dipole moment -- ideally calibrated from many molecules
      virtual double DipoleScalingFactor() { return 1.0; }
};

/** \class OBChargeModel chargemodel.h <openbabel/chargemodel.h>
      \brief Atomic partial charge models
      \since version 2.3

Classes derived from OBChargeModel implement different atomic partial
charge models. It is intended to allow assinging partial charges
beyond the traditional Gasteiger-Marsili sigma charges previously used
in Open Babel. A --partialcharge method is provided for the obabel
command-line, allowing you to override the Gasteiger charge assignment
and use other charge models.

The advantage of plugin classes is that no existing code has to be modified
when a new class is added. You can list those that are present by
obabel -L charges
or from a menu item in the GUI.

Any OBChargeModel derived class works like other plugins and needs to
to have a constructor, a function returning a short description, and a
ComputeCharges() function which does the work. A single global
instance of the class needs to be instantiated to define the ID, by
which the class is subsequently accessed.

Once ComputeCharges() has been called, the atoms of the molecule can
be queried for partial or formal charges using
OBAtom::GetPartialCharge() or in vector form from the model itself:

\code
  OBMol inputMolecule;
  OBChargeModel *mmffCharges = OBChargeModel::FindType("mmff94");
  const std::vector<double> partialCharges;
  if (mmffCharges && mmffCharges->ComputeCharges(inputMolecule)) {
    partialCharges = mmffCharges->GetPartialCharges();
  }
\endcode

Note: Formal charges are also returned as floating point values, since
some charge models consider delocalized charges (e.g., 0.5 for an O in
a carboxylate CO2- group).

\code
  OBChargeModel *gasteiger = OBChargeModel::FindType("gasteiger");
  if (gasteiger) {
    cout << " gasteiger: " << dipoleMagnitude(gasteiger->GetDipoleMoment(mol));
  }
\endcode

By default, Open Babel 2.3 includes Gasteiger and MMFF94 partial
charges. If the Eigen matrix library is found when compiling, the QEq
and QTPIE methods will be added. Future releases will likely add
additional charge models, including the EEM method.

  */

}//namespace
#endif

//! \file chargemodel.h
//! \brief Base class for molecular partial charge models

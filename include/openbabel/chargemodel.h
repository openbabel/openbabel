/**********************************************************************
chargemodel.h - Base class for partial charge models
 
Copyright (C) 2010 by Geoffrey Hutchison
Some portions Copyright (C) 2009 by Frank Peters

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

#ifndef OB_CHARGEMODEL_H
#define OB_CHARGEMODEL_H

#include <openbabel/babelconfig.h>
#include <openbabel/plugin.h>

namespace OpenBabel
{
class OBMol; //Forward declaration; used only as pointer.

class OBAPI OBChargeModel : public OBPlugin
{
  MAKE_PLUGIN(OBChargeModel)

  public:
    const char* TypeID(){return "charges";};

    /// \return whether partial charges were successfully assigned to this molecule
    /// \note The method should fill m_partialCharges and m_formalCharges as well
    virtual bool ComputeCharges(OBMol &) { return false; }

    /// \return a vector of the formal charges on each atom, indexed from 0
    /// \note If ComputeCharges() has not been called, this will return an empty vector
    const std::vector<double> & GetFormalCharges() const
    { return m_formalCharges; }

    /// \return a vector of the partial charges on each atom, indexed from 0
    /// \note If ComputeCharges() has not been called, this will return an empty vector
    const std::vector<double> & GetPartialCharges() const
    { return m_partialCharges; }

 protected:
      std::vector<double> m_partialCharges;
      std::vector<double> m_formalCharges;

      /// Fill the internal partial and formal charge vectors (convenience function)
      void FillChargeVectors(OBMol &mol);
};

}//namespace
#endif

//! \file chargemodel.h
//! \brief Base class for molecular partial charge models

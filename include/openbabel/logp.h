/**********************************************************************
logp.h - Handle logP prediction algorithms.
 
Copyright (C) 2007 by Tim Vandermeersch
 
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

#ifndef OB_LOGP_H
#define OB_LOGP_H


#include <openbabel/mol.h>
#include <openbabel/parsmart.h>

namespace OpenBabel
{

  // class introduction in logp.cpp
  class OBAPI OBLogP
  {
    public:
      //! constructor
      OBLogP();
      //! destructor
      ~OBLogP();
      /*! Predict the logP for molecule mol using the group contributions
       *  algorithm from JOELib2.
       *
       *  \param mol OBMol object for which to predict the logP
       *  \return predicted logP
       */
      double GroupContributions(OBMol &mol);
    private:
      std::vector<std::pair<OBSmartsPattern*, double> > _logPcontribsHeavy; //! logP contributions
      std::vector<std::pair<OBSmartsPattern*, double> > _logPcontribsHydrogen; //! logP contributions
  };

} // end namespace OpenBabel

#endif // OB_LOGP_H

//! \file logp.h
//! \brief Handle logP prediction algorithms.

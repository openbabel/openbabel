/**********************************************************************
groupcontrib.h - Handle group contribution algorithms.

Copyright (C) 2007      by Tim Vandermeersch
              2001-2007 by Stephen Jelfs
              2001-2007 by Joerg Kurt Wegner, me@cheminformatics.eu

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

#ifndef OB_GROUPCONTRIB_H
#define OB_GROUPCONTRIB_H

#include <openbabel/parsmart.h>
#include <openbabel/descriptor.h>

// This macro is used in DLL builds. If it has not
// been set in babelconfig.h, define it as nothing.
#ifndef OBDESC
	#define OBDESC
#endif

namespace OpenBabel
{

  /** \class OBGroupContrib groupcontrib.h <openbabel/groupcontrib.h>
      \brief Handle group contribution algorithms.

      This is the base class for calculations that use the JOELib2 contribution
      algorithm.
    */
class OBDESC OBGroupContrib : public OBDescriptor
{
public:

  /*! Predict the logP, MR, TPSA (each instance of OBGroupContrib
   *  uses different parameters loaded from its own datafile) for
   *  molecule mol using the group contributions algorithm from JOELib2.
   */

  //! constructor. Each instance provides an ID and a datafile.
  OBGroupContrib(const char* ID, const char* filename, const char* descr)
    : OBDescriptor(ID, false), _filename(filename), _descr(descr), _debug(false){}

  virtual const char* Description();

  virtual OBGroupContrib* MakeInstance(const std::vector<std::string>& textlines)
  {
    return new OBGroupContrib(textlines[1].c_str(),textlines[2].c_str(),textlines[3].c_str());
  }


  virtual double Predict(OBBase* pOb, std::string* param=NULL);

 private:
  bool ParseFile();

  const char* _filename;
  const char* _descr;
  std::vector<std::pair<OBSmartsPattern*, double> > _contribsHeavy; //! heavy atom contributions
  std::vector<std::pair<OBSmartsPattern*, double> > _contribsHydrogen; //!  hydrogen contributions
  bool _debug;
};

/* The classes OBLogp, OBPSA and OBMR have been replaced by instances of
OBGroupContrib with different IDs.
So instead of:
      OBLogp logP;
      cout << "logP  " << logP.Predict(mol) << endl;
use:
      OBDescriptor* pDesc = OBDescriptor::FindType("logP")
      if(pDesc)
        cout << "logP  " << pDesc->Predict(&mol) << endl;
*/
} // end namespace OpenBabel

#endif // OB_GROUPCONTRIB_H

//! \file groupcontrib.h
//! \brief Handle group contribution algorithms.

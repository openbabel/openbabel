/**********************************************************************
groupcontrib.h - Handle group contribution algorithms.
 
Copyright (C) 2007-2008 by Tim Vandermeersch
              2001-2007 by Stephen Jelfs
              2001-2007 by Joerg Kurt Wegner, me@cheminformatics.eu

Original version: JOELib2, http://joelib.sf.net
 
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

#ifndef OB_GROUPCONTRIB_H
#define OB_GROUPCONTRIB_H

#include <openbabel/mol.h>
#include <openbabel/parsmart.h>
#include <openbabel/descriptor.h>

// This macro is used in DLL builds. If it has not
// been set in babelconfig.h, define it as nothing.
#ifndef OBDESC
	#define OBDESC
#endif

namespace OpenBabel
{

  /** @class OBGroupContrib groupcontrib.h <openbabel/groupcontrib.h>
      @brief Handle group contribution algorithms.
 
      Predict the logP, MR, TPSA (each instance of OBGroupContrib 
      uses different parameters loaded from its own datafile) for 
      molecule mol using the group contributions algorithm from JOELib2.

      Article for logP and MR: \n
      <b>S. A. Wildman and G. M. Crippen \n
      Prediction of Physicochemical Parameters by Atomic Contributions \n
      J. Chem. Inf. Comput. Sci.,1999, 39,868-873 \n
      abstract (doi = 10.1021/ci990307l): </b>\n
      <i>We present a new atom type classification system for use in atom-based
      calculation of partition coefficient (log P) and molar refractivity
      (MR) designed in part to address published concerns of previous
      atomic methods. The 68 atomic contributions to log P have been determined
      by fitting an extensive training set of 9920 molecules, with r 2) 0.918 and
      0.677. A separate set of 3412 molecules was used
      for the determination of contributions to MR with r 2 ) 0.997 and
      1.43. Both calculations are showin to have high predictive ability.</i>

      Article for TPSA: \n
      <b>P. Ertl and B. Rohde and P. Selzer \n
      Fast Calculation of Molecular Polar Surface Area as a
      Sum of Fragment-Based Contributions and Its Application
      to the Prediction of Drug Transport Properties. \n
      J. Med. Chem., 2000, 43, 3714-3717 \n
      abstract: (doi = 10.1021/jm000942e): </b>\n
      <i>Molecular polar surface area (PSA), i.e., surface belonging to polar
      atoms, is a descriptor that was shown to correlate well with passive
      molecular transport through membranes and, therefore, allows prediction
      of transport properties of drugs. The calculation of PSA, however,
      is rather time-consuming because of the necessity to generate a
      reasonable 3D molecular geometry and the calculation of the surface
      itself. A new approach for the calculation of the PSA is presented
      here, based on the summation of tabulated surface contributions
      of polar fragments. The method, termed topological PSA (TPSA), provides
      results which are practically identical with the 3D PSA (the correlation
      coefficient between 3D PSA and fragment-based TPSA for 34 810 molecules
      from the World Drug Index is 0.99), while the computation speed
      is 2-3 orders of magnitude faster. The new methodology may, therefore,
      be used for fast bioavailability screening of virtual libraries
      having millions of molecules. This article describes the new methodology
      and shows the results of validation studies based on sets of published
      absorption data, including intestinal absorption, Caco-2 monolayer
      penetration, and blood- brain barrier penetration.</i>

      Octanol/water partition coefficient (logP):
      @code
      #include <openbabel/groupcontrib.h>
      OBMol mol; // see OBConversion on how to fill a mol object (from file, smiles, ...)
      OBDescriptor* pDesc = OBDescriptor::FindType("logP");
      if(pDesc)
        cout << "logP  " << pDesc->Predict(&mol) << endl;
      @endcode

      Topological polar surface area (TPSA):
      @code
      #include <openbabel/groupcontrib.h>
      OBMol mol; // see OBConversion on how to fill a mol object (from file, smiles, ...)
      OBDescriptor* pDesc = OBDescriptor::FindType("TPSA");
      if(pDesc)
        cout << "TPSA  " << pDesc->Predict(&mol) << endl;
      @endcode

      Molar refractivity (MR):
      @code
      #include <openbabel/groupcontrib.h>
      OBMol mol; // see OBConversion on how to fill a mol object (from file, smiles, ...)
      OBDescriptor* pDesc = OBDescriptor::FindType("MR");
      if(pDesc)
        cout << "MR  " << pDesc->Predict(&mol) << endl;
      @endcode
   */
  class OBDESC OBGroupContrib : public OBDescriptor
  {
    public:
      /** 
       * Constructor. Each instance provides an ID and a datafile.
       * @param ID The descriptor ID (logP, TPSA, MR)
       * @param filename The filename ro read the contributions from.
       * @param descr Description for this descriptor.
       */
      OBGroupContrib(const char* ID, const char* filename, const char* descr)
      : OBDescriptor(ID, false), m_filename(filename), m_descr(descr){}
      /**
       * @return Description.
       */
      virtual const char* Description();
      /**
       * Create a new instance of this groupcontrib descriptor.
       * @param textlines Vector of 3 strings. textlines[1] is the ID, textlines[2] is 
       * the filename and textlines[3] is the description.
       */
      virtual OBGroupContrib* MakeInstance(const std::vector<std::string>& textlines)
      {
        return new OBGroupContrib(textlines[1].c_str(), textlines[2].c_str(), textlines[3].c_str());
      }
      /**
       * Calculate the descriptor for molecule @p pOb.
       * @param pOb Pointer to OBMol object.
       * @return The predicted descriptor value for molecule @p pOb.
       */
      virtual double Predict(OBBase* pOb); 

    private:
      /**
       * Parse the contribution file.
       */
      bool ParseFile();

      const char* m_filename; //!< filename for the file with contributions
      const char* m_descr; //!< description for the descriptor
      std::vector<std::pair<OBSmartsPattern*, double> > m_contribsHeavy; //!< heavy atom contributions
      std::vector<std::pair<OBSmartsPattern*, double> > m_contribsHydrogen; //!<  hydrogen contributions
  };

} // end namespace OpenBabel

#endif // OB_GROUPCONTRIB_H

//! @file groupcontrib.h
//! @brief Handle group contribution algorithms.

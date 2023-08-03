/**********************************************************************
forcefieldmm2.h - MM2 force field.

Copyright (C) 2006 by Tim Vandermeersch <tim.vandermeersch@gmail.com>

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

#include <vector>
#include <string>
#include <map>

#include <openbabel/forcefield.h>
#include <openbabel/base.h>
#include <openbabel/mol.h>

namespace OpenBabel
{
  // Class OBForceFieldMM2
  // class introduction in forcefield.cpp
  class OBForceFieldMM2: public OBForceField
  {
    protected:

      bool _init; //!< Used to only initialize and read parameters once

      //! \return Parses the parameter file
      bool ParseParamFile();
      //! \return Sets atomtypes to MM2 in _mol
      bool SetMM2Types();

      double bondunit, bond_cubic, bond_quartic;
      double angleunit, angle_sextic;
      double stretchbendunit;
      double torsionunit;
      double outplanebendunit;
      double a_expterm, b_expterm, c_expterm;
      double dielectric;
      std::vector<OBFFParameter> _ffbondparams; // a = atom 1 of bond
                                                // b = atom 2 of bond
						// dpar1 = length
						// dpar2 = force
      std::vector<OBFFParameter> _ffangleparams; // a = atom 1 of angle abc
                                                 // b = atom 2 of angle abc
						 // c = atom 3 of angle abc
						 // dpar1 = angle
						 // dpar2 = force
      std::vector<OBFFParameter> _ffstretchbendparams; // a = atom
                                                        // dpar1 = force
      std::vector<OBFFParameter> _fftorsionparams; // a = atom 1 of torsion
                                                   // b = atom 2 of torsion
                                                   // c = atom 3 of torsion
                                                   // d = atom 4 of torsion
						   // dpar1 = v1
						   // dpar2 = v2
						   // dpar3 = v3
      std::vector<OBFFParameter> _ffoutplanebendparams; // a = atom b        a
                                                        // b = atom d         \
							// dpar1 = force       b---d
							//                    /
							//                   c
      std::vector<OBFFParameter> _ffvdwprparams;  // a = atom 1 of pair
                                                  // b = atom 2 of pair
					    	  // dpar1 = sum of vdw radii
						  // dpar2 = energy parameter
      std::vector<OBFFParameter> _ffvdwparams;    // a = atom 1
					    	  // dpar1 = vdw radii
						  // dpar2 = energy parameter
						  // dpar3 = reduction
      std::vector<OBFFParameter> _ffdipoleparams;    // a = atom 1
                                                     // b = atom 2
						     // dpar1 = dipole
						     // dpar2 = position

      std::vector<vector3> forces; // used to hold forces on each atom

    public:
      //! Setup
      bool Setup(OBMol &mol);
      //! Constructor
      explicit OBForceFieldMM2(const char* ID, bool IsDefault=true) : OBForceField(ID, IsDefault), _init(false)
      {
        // ParseParamFile only called when needed
      }

      //!Clone the current instance. May be desirable in multithreaded environments
      virtual OBForceFieldMM2* MakeNewInstance(){ return new OBForceFieldMM2(*this); }


      virtual const char* Description()
      { return "MM2 force field.";};


      //! Destructor
      virtual ~OBForceFieldMM2();
      //! Assignment
      OBForceFieldMM2 &operator = (OBForceFieldMM2 &);
      //! Returns total energy
      double Energy();
     //! Returns the bond stretching energy
      double E_Bond();
      //! Returns the angle bending energy
      double E_Angle();
      //! Returns the stretch-bend energy
      double E_StrBnd();
      //! Returns the torsional energy
      double E_Torsion();
      //! Returns the out-of-plane bending energy
      double E_OOP();
      //! Returns the Van der Waals energy (Buckingham potential)
      double E_VDW();
      //! Returns the dipole-dipole interaction energy
      double E_Electrostatic();


  }; // class OBForceFieldMM2

}// namespace OpenBabel

//! \file forcefieldmm2.h
//! \brief MM2 force field

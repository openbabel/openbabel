/**********************************************************************
forcefield.h - Handle OBForceField class. Molecular Mechanics force fields.
 
Copyright (C) 2006 by Tim Vandermeersch <tim.vandermeersch@gmail.com>
Some portions Copyright (C) 2003 by Michael Banck
 
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

#ifndef OB_FORCEFIELD_H
#define OB_FORCEFIELD_H

#include <vector>
#include <string>
#include <map>

#include <openbabel/base.h>
#include <openbabel/mol.h>

namespace OpenBabel
{

  class OBFFParameter {
  public:
    int       a, b, c, d; // used to store integer atom types
    char      _a[5], _b[5], _c[5], _d[5]; // used to store ascii atom types

    int       ipar1, ipar2, ipar3, ipar4, ipar5;
    double    dpar1, dpar2, dpar3, dpar4, dpar5;

    void clear () {
      a = 0;
      b = 0;
      c = 0;
      d = 0;
      ipar1 = 0;
      ipar2 = 0;
      ipar3 = 0;
      ipar4 = 0;
      ipar5 = 0;
      dpar1 = 0.0f;
      dpar2 = 0.0f;
      dpar3 = 0.0f;
      dpar4 = 0.0f;
      dpar5 = 0.0f;
    }
  };

  // Class OBForceField
  // class introduction in forcefield.cpp
  class OBAPI OBForceField
  {
  protected:
    std::vector<OBFFParameter> _ffbondparams; // a = atom 1 of bond
    // b = atom 2 of bond
    // dpar1 = lenght
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

    OBMol _mol;

    //! Constructor
    OBForceField();
    //! Destructor
    virtual ~OBForceField();
    //! Assignment
    OBForceField &operator = (OBForceField &);
    //! Get index for std::vector<OBFFParameter> ...
    int GetParameterIdx(int a, int b, int c, int d, std::vector<OBFFParameter> &parameter);
    //! Calculate angle between vector b-d and plane a-b-c
    double PointPlaneAngle(const vector3 &a, const vector3 &b, const vector3 &c, const vector3 &d);
 
  public:
    //! \return the bond stretching energy
    double E_Bond_Harmonic(double bondunit, double bond_cubic, double bond_quartic);
    //! \return the angle bending energy
    double E_Angle_Harmonic(double angleunit, double angle_cubic, double angle_quatric, double angle_pentic, double angle_sextic);
    //! \return the stretch-bend energy
    double E_StretchBend(double stretchbendunit);
    //! \return the torsional energy
    double E_Torsion(double torsionunit);
    //! \return the out-of-plane bending energy
    double E_OutPlaneBend(double outplanebendunit);
    //! \return the Van der Waals energy (Buckingham potential)
    double E_VDWBuckingham(double a_expterm, double b_expterm, double c_expterm);
    //! \return the Van der Waals energy (Buckingham potential)
    double E_DipoleDipole(double dielectric);
  }; // class OBForceField
    
  // Class OBForceFieldMM2
  // class introduction in forcefield.cpp
 class OBAPI OBForceFieldMM2: public OBForceField
  {
  protected:
    double bondunit, bond_cubic, bond_quartic;
    double angleunit, angle_sextic;
    double stretchbendunit;
    double torsionunit;
    double outplanebendunit;
    double a_expterm, b_expterm, c_expterm;
    double dielectric;

    //! \return Parses the parameter file
    bool ParseParamFile();
    //! \return Sets atomtypes to MM2 in _mol
    bool SetMM2Types();

  public:
    //! Setup
    bool Setup(OBMol &mol);
    //! Constructor
    OBForceFieldMM2();
    //! Destructor
    virtual ~OBForceFieldMM2();
    //! Assignment
    OBForceFieldMM2 &operator = (OBForceFieldMM2 &);
    //! \return total energy
    double Energy();
    double E_Bond() { return E_Bond_Harmonic(bondunit, bond_cubic, bond_quartic); }
    double E_Angle() { return E_Angle_Harmonic(angleunit, 0.0f, 0.0f, 0.0f, angle_sextic); }
    double E_StrBnd() { return E_StretchBend(stretchbendunit); }
    double E_Tor() { return E_Torsion(torsionunit); }
    double E_OPBend() { return E_OutPlaneBend(outplanebendunit); }
    double E_VDW() { return E_VDWBuckingham(a_expterm, b_expterm, c_expterm); }
    double E_Dipole() { return E_DipoleDipole(dielectric); }


  }; // class OBForceFieldMM2

}// namespace OpenBabel

#endif   // OB_FORCEFIELD_H

//! \file forcefield.h
//! \brief Handle molecular mechanics force fields

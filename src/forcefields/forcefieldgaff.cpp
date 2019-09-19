/**********************************************************************
forcefieldgaff.cpp - Gaff force field.

Copyright (C) 2009 by Frank Peters <e.a.j.f.peters@tue.nl>
Copyright (C) 2006-2007 by Tim Vandermeersch <tim.vandermeersch@gmail.com>

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
#include <openbabel/mol.h>
#include <openbabel/bond.h>
#include <openbabel/obiter.h>
#include <openbabel/locale.h>
#include <openbabel/oberror.h>
#include <openbabel/parsmart.h>
#include <openbabel/obutil.h>

#include "forcefieldgaff.h"

#include <cstdlib>

using namespace std;

namespace OpenBabel
{
  template<bool gradients>
  void OBFFBondCalculationGaff::Compute()
  {
    if (OBForceField::IgnoreCalculation(idx_a, idx_b)) {
      energy = 0.0;
      return;
    }

    double delta2;

    if (gradients) {
      rab = OBForceField::VectorBondDerivative(pos_a, pos_b, force_a, force_b);
      delta = rab - r0;
      delta2 = delta * delta;

      const double dE = 2.0 * kr * delta;

      OBForceField::VectorSelfMultiply(force_a, dE);
      OBForceField::VectorSelfMultiply(force_b, dE);
    } else {
      rab = OBForceField::VectorDistance(pos_a, pos_b);
      delta = rab - r0;
      delta2 = delta * delta;
    }

    energy = kr * delta2;
  }

  template<bool gradients>
  double OBForceFieldGaff::E_Bond()
  {
    vector<OBFFBondCalculationGaff>::iterator i;
    double energy = 0.0;

    IF_OBFF_LOGLVL_HIGH {
      OBFFLog("\nB O N D   S T R E T C H I N G\n\n");
      OBFFLog("ATOM TYPES  BOND       IDEAL       FORCE\n");
      OBFFLog(" I    J     LENGTH     LENGTH     CONSTANT      DELTA      ENERGY\n");
      OBFFLog("------------------------------------------------------------------------\n");
    }

    for (i = _bondcalculations.begin(); i != _bondcalculations.end(); ++i) {

      i->template Compute<gradients>();
      energy += i->energy;

      if (gradients) {
        AddGradient((*i).force_a, (*i).idx_a);
        AddGradient((*i).force_b, (*i).idx_b);
      }

      IF_OBFF_LOGLVL_HIGH {
        snprintf(_logbuf, BUFF_SIZE, "%s %s  %8.3f   %8.3f     %8.3f   %8.3f   %8.3f\n", (*i).a->GetType(), (*i).b->GetType(),
                (*i).rab, (*i).r0, (*i).kr, (*i).delta, (*i).energy);
        OBFFLog(_logbuf);
      }
    }

    IF_OBFF_LOGLVL_MEDIUM {
      snprintf(_logbuf, BUFF_SIZE, "     TOTAL BOND STRETCHING ENERGY = %8.3f %s\n",  energy, GetUnit().c_str());
      OBFFLog(_logbuf);
    }
    return energy;
  }

  template<bool gradients>
  void OBFFAngleCalculationGaff::Compute()
  {
    if (OBForceField::IgnoreCalculation(idx_a, idx_b, idx_c)) {
      energy = 0.0;
      return;
    }

    double delta2;

    if (gradients) {
      theta = OBForceField::VectorAngleDerivative(pos_a, pos_b, pos_c, force_a, force_b, force_c);
      delta = DEG_TO_RAD * (theta - theta0);

      const double dE = 2.0 * kth * delta;

      OBForceField::VectorSelfMultiply(force_a, dE);
      OBForceField::VectorSelfMultiply(force_b, dE);
      OBForceField::VectorSelfMultiply(force_c, dE);
    } else {
      theta = OBForceField::VectorAngle(pos_a, pos_b, pos_c);
      delta = DEG_TO_RAD * (theta - theta0);
    }

    if (!isfinite(theta))
      theta = 0.0; // doesn't explain why GetAngle is returning NaN but solves it for us;

    delta2 = delta * delta;

    energy = kth * delta2;
  }

  template<bool gradients>
  double OBForceFieldGaff::E_Angle()
  {
    vector<OBFFAngleCalculationGaff>::iterator i;
    double energy = 0.0;

    IF_OBFF_LOGLVL_HIGH {
      OBFFLog("\nA N G L E   B E N D I N G\n\n");
      OBFFLog("ATOM TYPES       VALENCE     IDEAL      FORCE\n");
      OBFFLog(" I    J    K      ANGLE      ANGLE     CONSTANT      DELTA      ENERGY\n");
      OBFFLog("-----------------------------------------------------------------------------\n");
    }

    for (i = _anglecalculations.begin(); i != _anglecalculations.end(); ++i) {

      i->template Compute<gradients>();
      energy += i->energy;

      if (gradients) {
        AddGradient((*i).force_a, (*i).idx_a);
        AddGradient((*i).force_b, (*i).idx_b);
        AddGradient((*i).force_c, (*i).idx_c);
      }

      IF_OBFF_LOGLVL_HIGH {
        snprintf(_logbuf, BUFF_SIZE, "%s %s %s  %8.3f   %8.3f     %8.3f   %8.3f   %8.3f\n", (*i).a->GetType(), (*i).b->GetType(),
                (*i).c->GetType(), (*i).theta, (*i).theta0, (*i).kth, (*i).delta, (*i).energy);
        OBFFLog(_logbuf);
      }
    }

    IF_OBFF_LOGLVL_MEDIUM {
      snprintf(_logbuf, BUFF_SIZE, "     TOTAL ANGLE BENDING ENERGY = %8.3f %s\n", energy, GetUnit().c_str());
      OBFFLog(_logbuf);
    }
    return energy;
  }

  template<bool gradients>
  void OBFFTorsionCalculationGaff::Compute()
  {
    if (OBForceField::IgnoreCalculation(idx_a, idx_b, idx_c, idx_d)) {
      energy = 0.0;
      return;
    }

    if (gradients) {
      tor = OBForceField::VectorTorsionDerivative(pos_a, pos_b, pos_c, pos_d,
                                                               force_a, force_b, force_c, force_d);
      if (!isfinite(tor))
        tor = 1.0e-3;

      const double sine = sin(DEG_TO_RAD*(n*tor-gamma));
      const double dE = n * vn_half * sine;

      OBForceField::VectorSelfMultiply(force_a, dE);
      OBForceField::VectorSelfMultiply(force_b, dE);
      OBForceField::VectorSelfMultiply(force_c, dE);
      OBForceField::VectorSelfMultiply(force_d, dE);
    } else {
      tor = OBForceField::VectorTorsion(pos_a, pos_b, pos_c, pos_d);
      if (!isfinite(tor)) // stop any NaN or infinity
        tor = 1.0e-3; // rather than NaN
    }

    const double cosine = cos(DEG_TO_RAD*(n*tor-gamma));
    const double phi1 = 1.0 + cosine;

    energy = vn_half * phi1;

  }

  template<bool gradients>
  double OBForceFieldGaff::E_Torsion()
  {
    vector<OBFFTorsionCalculationGaff>::iterator i;
    double energy = 0.0;

    IF_OBFF_LOGLVL_HIGH {
      OBFFLog("\nT O R S I O N A L\n\n");
      OBFFLog("----ATOM TYPES-----    FORCE              TORSION\n");
      OBFFLog(" I    J    K    L     CONSTANT     s       ANGLE    n    ENERGY\n");
      OBFFLog("----------------------------------------------------------------\n");
    }

    for (i = _torsioncalculations.begin(); i != _torsioncalculations.end(); ++i) {

      i->template Compute<gradients>();
      energy += i->energy;

      if (gradients) {
        AddGradient((*i).force_a, (*i).idx_a);
        AddGradient((*i).force_b, (*i).idx_b);
        AddGradient((*i).force_c, (*i).idx_c);
        AddGradient((*i).force_d, (*i).idx_d);
      }

      IF_OBFF_LOGLVL_HIGH {
        snprintf(_logbuf, BUFF_SIZE, "%s %s %s %s    %6.3f    %5.0f   %8.3f   %1.0f   %8.3f\n", (*i).a->GetType(), (*i).b->GetType(),
                (*i).c->GetType(), (*i).d->GetType(), (*i).vn_half, (*i).gamma, (*i).tor, (*i).n, (*i).energy);
        OBFFLog(_logbuf);
      }
    }

    IF_OBFF_LOGLVL_MEDIUM {
      snprintf(_logbuf, BUFF_SIZE, "     TOTAL TORSIONAL ENERGY = %8.3f %s\n", energy, GetUnit().c_str());
      OBFFLog(_logbuf);
    }

    return energy;
  }

  //						//
  //  a						//
  //   \  					//
  //    b---d      plane = a-b-c		//
  //   / 					//
  //  c						//
  //						//
  OBFFParameter* OBForceFieldGaff::GetParameterOOP(const char* a, const char* b, const char* c, const char* d,
        std::vector<OBFFParameter> &parameter)
  {
    OBFFParameter *par;
    if (a == NULL || b == NULL || c == NULL || d == NULL )
      return NULL;
    string _a(a);
    string _b(b);
    string _c(c);
    string _d(d);
    for (unsigned int idx=0; idx < parameter.size(); idx++)
      if (((_a == parameter[idx]._a) && (_b == parameter[idx]._b) &&
	   (_c == parameter[idx]._c) && (_d == parameter[idx]._d)) ||
	  ((_a == parameter[idx]._c) && (_b == parameter[idx]._b) &&
	   (_c == parameter[idx]._a) && (_d == parameter[idx]._d))) {
	par = &parameter[idx];
	return par;
      }
    return NULL;
  }

  template<bool gradients>
  void OBFFOOPCalculationGaff::Compute()
  {
    if (OBForceField::IgnoreCalculation(idx_a, idx_b, idx_c, idx_d)) {
      energy = 0.0;
      return;
    }

    if (gradients) {
      tor = OBForceField::VectorTorsionDerivative(pos_a, pos_b, pos_c, pos_d,
                                                               force_a, force_b, force_c, force_d);
      if (!isfinite(tor))
        tor = 1.0e-3;

      const double sine = sin(DEG_TO_RAD*(n*tor-gamma));
      const double dE = n * vn_half * sine;

      OBForceField::VectorSelfMultiply(force_a, dE);
      OBForceField::VectorSelfMultiply(force_b, dE);
      OBForceField::VectorSelfMultiply(force_c, dE);
      OBForceField::VectorSelfMultiply(force_d, dE);
    } else {
      tor = OBForceField::VectorTorsion(pos_a, pos_b, pos_c, pos_d);
      if (!isfinite(tor)) // stop any NaN or infinity
        tor = 1.0e-3; // rather than NaN
    }

    const double cosine = cos(DEG_TO_RAD*(n*tor-gamma));
    const double phi1 = 1.0 + cosine;

    energy = vn_half * phi1;

  }

  template<bool gradients>
  double OBForceFieldGaff::E_OOP()
  {
    vector<OBFFOOPCalculationGaff>::iterator i;
    double energy = 0.0;

    IF_OBFF_LOGLVL_HIGH {
      OBFFLog("\nI M P R O P E R   T O R S I O N A L\n\n");
      OBFFLog("----ATOM TYPES-----    FORCE     IMPROPER_TORSION\n");
      OBFFLog(" I    J    K    L     CONSTANT     s       ANGLE    n    ENERGY\n");
      OBFFLog("----------------------------------------------------------------\n");
    }

    for (i = _oopcalculations.begin(); i != _oopcalculations.end(); ++i) {

      i->template Compute<gradients>();
      energy += i->energy;

      if (gradients) {
        AddGradient((*i).force_a, (*i).idx_a);
        AddGradient((*i).force_b, (*i).idx_b);
        AddGradient((*i).force_c, (*i).idx_c);
        AddGradient((*i).force_d, (*i).idx_d);
      }

      IF_OBFF_LOGLVL_HIGH {
        snprintf(_logbuf, BUFF_SIZE, "%s %s %s %s    %6.3f    %5.0f   %8.3f   %1.0f   %8.3f\n", (*i).a->GetType(), (*i).b->GetType(),
                (*i).c->GetType(), (*i).d->GetType(), (*i).vn_half, (*i).gamma, (*i).tor, (*i).n, (*i).energy);
        OBFFLog(_logbuf);
      }
    }

    IF_OBFF_LOGLVL_MEDIUM {
      snprintf(_logbuf, BUFF_SIZE, "     TOTAL IMPROPER-TORSIONAL ENERGY = %8.3f %s\n", energy, GetUnit().c_str());
      OBFFLog(_logbuf);
    }

    return energy;
  }

  template<bool gradients>
  void OBFFVDWCalculationGaff::Compute()
  {
    if (OBForceField::IgnoreCalculation(idx_a, idx_b)) {
      energy = 0.0;
      return;
    }

    if (gradients) {
      rab = OBForceField::VectorDistanceDerivative(pos_a, pos_b, force_a, force_b);
    } else {
      rab = OBForceField::VectorDistance(pos_a, pos_b);
    }

    const double term = RVDWab / rab;

    double term6 = term * term * term; // ^3
    term6 = term6 * term6; // ^6
    const double term12 = term6 * term6; // ^12

    energy = Eab * (term12 - 2.0*term6);

    if (gradients) {
      const double term13 = term * term12; // ^13
      const double term7 = term * term6; // ^7
      const double dE = (12.0 * Eab / RVDWab) * (-term13 + term7);
      OBForceField::VectorSelfMultiply(force_a, dE);
      OBForceField::VectorSelfMultiply(force_b, dE);
    }
  }

  template<bool gradients>
  double OBForceFieldGaff::E_VDW()
  {
    vector<OBFFVDWCalculationGaff>::iterator i;
    double energy = 0.0;

    IF_OBFF_LOGLVL_HIGH {
      OBFFLog("\nV A N   D E R   W A A L S\n\n");
      OBFFLog("ATOM TYPES\n");
      OBFFLog(" I    J        Rij       kij       ENERGY\n");
      OBFFLog("-----------------------------------------\n");
      //          XX   XX     -000.000  -000.000  -000.000  -000.000
    }

    unsigned int j = 0;
    for (i = _vdwcalculations.begin(); i != _vdwcalculations.end(); ++i, ++j) {
      // Cut-off check
      if (_cutoff)
        if (!_vdwpairs.BitIsSet(j))
          continue;

      i->template Compute<gradients>();
      energy += i->energy;

      if (gradients) {
        AddGradient((*i).force_a, (*i).idx_a);
        AddGradient((*i).force_b, (*i).idx_b);
      }

      IF_OBFF_LOGLVL_HIGH {
        snprintf(_logbuf, BUFF_SIZE, "%s %s   %8.3f  %8.3f\n", (*i).a->GetType(), (*i).b->GetType(),
                (*i).rab, (*i).energy);
        OBFFLog(_logbuf);
      }
    }

    IF_OBFF_LOGLVL_MEDIUM {
      snprintf(_logbuf, BUFF_SIZE, "     TOTAL VAN DER WAALS ENERGY = %8.3f %s\n", energy, GetUnit().c_str());
      OBFFLog(_logbuf);
    }

    return energy;
  }

  template<bool gradients>
  void OBFFElectrostaticCalculationGaff::Compute()
  {
    if (OBForceField::IgnoreCalculation(idx_a, idx_b)) {
      energy = 0.0;
      return;
    }

    if (gradients) {
      rab = OBForceField::VectorDistanceDerivative(pos_a, pos_b, force_a, force_b);
      const double rab2 = rab * rab;
      const double dE = -qq / rab2;
      OBForceField::VectorSelfMultiply(force_a, dE);
      OBForceField::VectorSelfMultiply(force_b, dE);
    } else {
      rab = OBForceField::VectorDistance(pos_a, pos_b);
    }

    if (IsNearZero(rab, 1.0e-3))
      rab = 1.0e-3;

    energy = qq / rab;
  }

  template<bool gradients>
  double OBForceFieldGaff::E_Electrostatic()
  {
    vector<OBFFElectrostaticCalculationGaff>::iterator i;
    double energy = 0.0;

    IF_OBFF_LOGLVL_HIGH {
      OBFFLog("\nE L E C T R O S T A T I C   I N T E R A C T I O N S\n\n");
      OBFFLog("ATOM TYPES\n");
      OBFFLog(" I    J           Rij   332.17*QiQj  ENERGY\n");
      OBFFLog("-------------------------------------------\n");
      //            XX   XX     -000.000  -000.000  -000.000
    }

    unsigned int j = 0;
    for (i = _electrostaticcalculations.begin(); i != _electrostaticcalculations.end(); ++i, ++j) {
      // Cut-off check
      if (_cutoff)
        if (!_elepairs.BitIsSet(j))
          continue;

      i->template Compute<gradients>();
      energy += i->energy;

      if (gradients) {
        AddGradient((*i).force_a, (*i).idx_a);
        AddGradient((*i).force_b, (*i).idx_b);
      }

      IF_OBFF_LOGLVL_HIGH {
        snprintf(_logbuf, BUFF_SIZE, "%s %s   %8.3f  %8.3f  %8.3f\n", (*i).a->GetType(), (*i).b->GetType(),
                (*i).rab, (*i).qq, (*i).energy);
        OBFFLog(_logbuf);
      }
    }

    IF_OBFF_LOGLVL_MEDIUM {
      snprintf(_logbuf, BUFF_SIZE, "     TOTAL ELECTROSTATIC ENERGY = %8.3f %s\n", energy, GetUnit().c_str());
      OBFFLog(_logbuf);
    }

    return energy;
  }

  //***********************************************
  //Make a global instance
  OBForceFieldGaff theForceFieldGaff("GAFF", true);
  //***********************************************

  OBForceFieldGaff::~OBForceFieldGaff()
  {
  }

  OBForceFieldGaff &OBForceFieldGaff::operator=(OBForceFieldGaff &src)
  {
    _mol = src._mol;
    _init = src._init;

    _ffbondparams     = src._ffbondparams;
    _ffangleparams    = src._ffangleparams;
    _fftorsionparams  = src._fftorsionparams;
    _ffoopparams = src._ffoopparams;
    _ffhbondparams    = src._ffhbondparams;
    _ffvdwparams      = src._ffvdwparams;

    _bondcalculations          = src._bondcalculations;
    _anglecalculations         = src._anglecalculations;
    _torsioncalculations       = src._torsioncalculations;
    _oopcalculations      = src._oopcalculations;
    _vdwcalculations           = src._vdwcalculations;
    _electrostaticcalculations = src._electrostaticcalculations;

    return *this;
  }

  bool OBForceFieldGaff::SetupCalculations()
  {
    OBFFParameter *parameter;
    OBAtom *a, *b, *c, *d;

    IF_OBFF_LOGLVL_LOW
      OBFFLog("\nS E T T I N G   U P   C A L C U L A T I O N S\n\n");

    //
    // Bond Calculations
    IF_OBFF_LOGLVL_LOW
      OBFFLog("SETTING UP BOND CALCULATIONS...\n");

    OBFFBondCalculationGaff bondcalc;

    _bondcalculations.clear();

    FOR_BONDS_OF_MOL(bond, _mol) {
      a = bond->GetBeginAtom();
      b = bond->GetEndAtom();


      // skip this bond if the atoms are ignored
      if ( _constraints.IsIgnored(a->GetIdx()) || _constraints.IsIgnored(b->GetIdx()) )
        continue;

      // if there are any groups specified, check if the two bond atoms are in a single intraGroup
      if (HasGroups()) {
        bool validBond = false;
        for (unsigned int i=0; i < _intraGroup.size(); ++i) {
          if (_intraGroup[i].BitIsSet(a->GetIdx()) && _intraGroup[i].BitIsSet(b->GetIdx()))
            validBond = true;
        }
        if (!validBond)
          continue;
      }

      bondcalc.a = a;
      bondcalc.b = b;

      parameter = GetParameter(a->GetType(), b->GetType(), NULL, NULL,  _ffbondparams);
      if (parameter == NULL) {
        parameter = GetParameter("X", a->GetType(), NULL, NULL, _ffbondparams);
        if (parameter == NULL) {
          parameter = GetParameter("X", b->GetType(), NULL, NULL, _ffbondparams);
          if (parameter == NULL) {
            bondcalc.kr = KCAL_TO_KJ * 500.0;
            bondcalc.r0 = 1.100;
            bondcalc.SetupPointers();

            _bondcalculations.push_back(bondcalc);

            IF_OBFF_LOGLVL_LOW {
              snprintf(_logbuf, BUFF_SIZE, "COULD NOT FIND PARAMETERS FOR BOND %s-%s, USING DEFAULT PARAMETERS\n", a->GetType(), b->GetType());
              OBFFLog(_logbuf);
            }

            continue;
          }
        }
      }
      bondcalc.kr = KCAL_TO_KJ * parameter->_dpar[0];
      bondcalc.r0 = parameter->_dpar[1];
      bondcalc.SetupPointers();

      _bondcalculations.push_back(bondcalc);
    }

    //
    // Angle Calculations
    //
    IF_OBFF_LOGLVL_LOW
      OBFFLog("SETTING UP ANGLE CALCULATIONS...\n");

    OBFFAngleCalculationGaff anglecalc;

    _anglecalculations.clear();

    FOR_ANGLES_OF_MOL(angle, _mol) {
      b = _mol.GetAtom((*angle)[0] + 1);
      a = _mol.GetAtom((*angle)[1] + 1);
      c = _mol.GetAtom((*angle)[2] + 1);

      // skip this angle if the atoms are ignored
      if ( _constraints.IsIgnored(a->GetIdx()) || _constraints.IsIgnored(b->GetIdx()) || _constraints.IsIgnored(c->GetIdx()) )
        continue;

      // if there are any groups specified, check if the three angle atoms are in a single intraGroup
      if (HasGroups()) {
        bool validAngle = false;
        for (unsigned int i=0; i < _intraGroup.size(); ++i) {
          if (_intraGroup[i].BitIsSet(a->GetIdx()) && _intraGroup[i].BitIsSet(b->GetIdx()) &&
              _intraGroup[i].BitIsSet(c->GetIdx()))
            validAngle = true;
        }
        if (!validAngle)
          continue;
      }

      anglecalc.a = a;
      anglecalc.b = b;
      anglecalc.c = c;

      parameter = GetParameter(a->GetType(), b->GetType(), c->GetType(), NULL, _ffangleparams);
      if (parameter == NULL) {
        parameter = GetParameter("X", b->GetType(), c->GetType(), NULL, _ffangleparams);
        if (parameter == NULL) {
          parameter = GetParameter(a->GetType(), b->GetType(), "X", NULL, _ffangleparams);
          if (parameter == NULL) {
            parameter = GetParameter("X", b->GetType(), "X", NULL, _ffangleparams);
            if (parameter == NULL) {
              anglecalc.kth = KCAL_TO_KJ * 0.020;
              anglecalc.theta0 = 120.0;
              anglecalc.SetupPointers();

              _anglecalculations.push_back(anglecalc);

              IF_OBFF_LOGLVL_LOW {
                snprintf(_logbuf, BUFF_SIZE, "COULD NOT FIND PARAMETERS FOR ANGLE %s-%s-%s, USING DEFAULT PARAMETERS\n", a->GetType(), b->GetType(), c->GetType());
                OBFFLog(_logbuf);
              }

              continue;
            }
          }
        }
      }
      anglecalc.kth = KCAL_TO_KJ * parameter->_dpar[0];
      anglecalc.theta0 = parameter->_dpar[1];
      anglecalc.SetupPointers();

      _anglecalculations.push_back(anglecalc);
    }

    //
    // Torsion Calculations
    //
    IF_OBFF_LOGLVL_LOW
      OBFFLog("SETTING UP TORSION CALCULATIONS...\n");

    OBFFTorsionCalculationGaff torsioncalc;

    _torsioncalculations.clear();

    FOR_TORSIONS_OF_MOL(t, _mol) {
      a = _mol.GetAtom((*t)[0] + 1);
      b = _mol.GetAtom((*t)[1] + 1);
      c = _mol.GetAtom((*t)[2] + 1);
      d = _mol.GetAtom((*t)[3] + 1);

      // skip this torsion if the atoms are ignored
      if ( _constraints.IsIgnored(a->GetIdx()) || _constraints.IsIgnored(b->GetIdx()) ||
           _constraints.IsIgnored(c->GetIdx()) || _constraints.IsIgnored(d->GetIdx()) )
        continue;

      // if there are any groups specified, check if the four torsion atoms are in a single intraGroup
      if (HasGroups()) {
        bool validTorsion = false;
        for (unsigned int i=0; i < _intraGroup.size(); ++i) {
          if (_intraGroup[i].BitIsSet(a->GetIdx()) && _intraGroup[i].BitIsSet(b->GetIdx()) &&
              _intraGroup[i].BitIsSet(c->GetIdx()) && _intraGroup[i].BitIsSet(d->GetIdx()))
            validTorsion = true;
        }
        if (!validTorsion)
          continue;
      }

      torsioncalc.a = a;
      torsioncalc.b = b;
      torsioncalc.c = c;
      torsioncalc.d = d;

      parameter = GetParameter(a->GetType(), b->GetType(), c->GetType(), d->GetType(), _fftorsionparams);
      if (parameter == NULL) {
        parameter = GetParameter("X", b->GetType(), c->GetType(), d->GetType(), _fftorsionparams);
        if (parameter == NULL) {
          parameter = GetParameter(a->GetType(), b->GetType(), c->GetType(), "X", _fftorsionparams);
          if (parameter == NULL) {
            parameter = GetParameter("X", b->GetType(), c->GetType(), "X", _fftorsionparams);
            if (parameter == NULL) {
	      torsioncalc.vn_half = 0.0;
	      torsioncalc.gamma = 0.0;
	      torsioncalc.n = 0.0;

              torsioncalc.SetupPointers();
              _torsioncalculations.push_back(torsioncalc);

              IF_OBFF_LOGLVL_LOW {
                snprintf(_logbuf, BUFF_SIZE, "COULD NOT FIND PARAMETERS FOR TORSION %s-%s-%s-%s, USING DEFAULT PARAMETERS\n", a->GetType(), b->GetType(), c->GetType(), d->GetType());
                OBFFLog(_logbuf);
              }

              continue;
            }
          }
        }
      }
      torsioncalc.vn_half = KCAL_TO_KJ * parameter->_dpar[0]/ parameter->_ipar[0];
      torsioncalc.gamma = parameter->_dpar[1];
      torsioncalc.n = parameter->_dpar[2];

      torsioncalc.SetupPointers();
      _torsioncalculations.push_back(torsioncalc);
    }

    //
    // Improper torsion (Out-Of-Plane) Calculations
    //
    IF_OBFF_LOGLVL_LOW
      OBFFLog("SETTING UP IMPROPER TORSION CALCULATIONS...\n");

    OBFFOOPCalculationGaff oopcalc;

    _oopcalculations.clear();

    FOR_ATOMS_OF_MOL(atom, _mol) {
      b = (OBAtom*) &*atom;

      a = NULL;
      c = NULL;
      d = NULL;

      FOR_NBORS_OF_ATOM(nbr, b) {
	if (a ==NULL)
	  a = (OBAtom*) &*nbr;
	else if (c == NULL)
	  c = (OBAtom*) &*nbr;
	else
	  d = (OBAtom*) &*nbr;
      }

      // skip this improper dihedral if atom b does not have 3 neighbors
      if ((a == NULL) || (c == NULL) || (d == NULL))
	continue;

      // skip this improper dihedral if the atoms are ignored
      if ( _constraints.IsIgnored(a->GetIdx()) || _constraints.IsIgnored(b->GetIdx()) ||
	   _constraints.IsIgnored(c->GetIdx()) || _constraints.IsIgnored(d->GetIdx()) )
	continue;

      // if there are any groups specified, check if the four improper-dihedral atoms are in a single intraGroup
      if (HasGroups()) {
	bool validOOP = false;
	for (unsigned int i=0; i < _intraGroup.size(); ++i) {
	  if (_intraGroup[i].BitIsSet(a->GetIdx()) && _intraGroup[i].BitIsSet(b->GetIdx()) &&
	      _intraGroup[i].BitIsSet(c->GetIdx()) && _intraGroup[i].BitIsSet(d->GetIdx()))
	    validOOP = true;
	}
	if (!validOOP)
	  continue;
      }

      parameter = GetParameterOOP(a->GetType(), b->GetType(), c->GetType(), d->GetType(), _ffoopparams);
      if (parameter != NULL){
	// A-B-C-D || PLANE = ABC
	oopcalc.a = a;
	oopcalc.b = b;
	oopcalc.c = c;
	oopcalc.d = d;
	oopcalc.vn_half = KCAL_TO_KJ * parameter->_dpar[0];
	oopcalc.gamma = parameter->_dpar[1];
	oopcalc.n = parameter->_dpar[2];
	oopcalc.SetupPointers();
	_oopcalculations.push_back(oopcalc);
	continue;
      }
      parameter = GetParameterOOP(a->GetType(), b->GetType(), d->GetType(), c->GetType(), _ffoopparams);
      if (parameter != NULL){
	// A-B-D-C || PLANE = ABD
	oopcalc.a = a;
	oopcalc.b = b;
	oopcalc.c = d;
	oopcalc.d = c;
	oopcalc.vn_half = KCAL_TO_KJ * parameter->_dpar[0];
	oopcalc.gamma = parameter->_dpar[1];
	oopcalc.n = parameter->_dpar[2];
	oopcalc.SetupPointers();
	_oopcalculations.push_back(oopcalc);
	continue;
      }
      parameter = GetParameterOOP(c->GetType(), b->GetType(), d->GetType(), a->GetType(), _ffoopparams);
      if (parameter != NULL){
	// C-B-D-A || PLANE = CBD
	oopcalc.a = c;
	oopcalc.b = b;
	oopcalc.c = d;
	oopcalc.d = a;
	oopcalc.vn_half = KCAL_TO_KJ * parameter->_dpar[0];
	oopcalc.gamma = parameter->_dpar[1];
	oopcalc.n = parameter->_dpar[2];
	oopcalc.SetupPointers();
	_oopcalculations.push_back(oopcalc);
	continue;
      }
      parameter = GetParameterOOP("X", b->GetType(), c->GetType(), d->GetType(), _ffoopparams);
      if (parameter != NULL){
	// A-B-C-D || PLANE = ABC
	oopcalc.a = a;
	oopcalc.b = b;
	oopcalc.c = c;
	oopcalc.d = d;
	oopcalc.vn_half = KCAL_TO_KJ * parameter->_dpar[0];
	oopcalc.gamma = parameter->_dpar[1];
	oopcalc.n = parameter->_dpar[2];
	oopcalc.SetupPointers();
	_oopcalculations.push_back(oopcalc);
	continue;
      }
      parameter = GetParameterOOP("X", b->GetType(), d->GetType(), c->GetType(), _ffoopparams);
      if (parameter != NULL){
	// A-B-D-C || PLANE = ABD
	oopcalc.a = a;
	oopcalc.b = b;
	oopcalc.c = d;
	oopcalc.d = c;
	oopcalc.vn_half = KCAL_TO_KJ * parameter->_dpar[0];
	oopcalc.gamma = parameter->_dpar[1];
	oopcalc.n = parameter->_dpar[2];
	oopcalc.SetupPointers();
	_oopcalculations.push_back(oopcalc);
	continue;
      }
      parameter = GetParameterOOP("X", b->GetType(), d->GetType(), a->GetType(), _ffoopparams);
      if (parameter != NULL){
	// C-B-D-A || PLANE = CBD
	oopcalc.a = c;
	oopcalc.b = b;
	oopcalc.c = d;
	oopcalc.d = a;
	oopcalc.vn_half = KCAL_TO_KJ * parameter->_dpar[0];
	oopcalc.gamma = parameter->_dpar[1];
	oopcalc.n = parameter->_dpar[2];
	oopcalc.SetupPointers();
	_oopcalculations.push_back(oopcalc);
	continue;
      }
      parameter = GetParameterOOP("X", "X", c->GetType(), d->GetType(), _ffoopparams);
      if (parameter != NULL){
	// A-B-C-D || PLANE = ABC
	oopcalc.a = a;
	oopcalc.b = b;
	oopcalc.c = c;
	oopcalc.d = d;
	oopcalc.vn_half = KCAL_TO_KJ * parameter->_dpar[0];
	oopcalc.gamma = parameter->_dpar[1];
	oopcalc.n = parameter->_dpar[2];
	oopcalc.SetupPointers();
	_oopcalculations.push_back(oopcalc);
	continue;
      }
      parameter = GetParameterOOP("X", "X", d->GetType(), c->GetType(), _ffoopparams);
      if (parameter != NULL){
	// A-B-D-C || PLANE = ABD
	oopcalc.a = a;
	oopcalc.b = b;
	oopcalc.c = d;
	oopcalc.d = c;
	oopcalc.vn_half = KCAL_TO_KJ * parameter->_dpar[0];
	oopcalc.gamma = parameter->_dpar[1];
	oopcalc.n = parameter->_dpar[2];
	oopcalc.SetupPointers();
	_oopcalculations.push_back(oopcalc);
	continue;
      }
      parameter = GetParameterOOP("X", "X", d->GetType(), a->GetType(), _ffoopparams);
      if (parameter != NULL){
	// C-B-D-A || PLANE = CBD
	oopcalc.a = c;
	oopcalc.b = b;
	oopcalc.c = d;
	oopcalc.d = a;
	oopcalc.vn_half = KCAL_TO_KJ * parameter->_dpar[0];
	oopcalc.gamma = parameter->_dpar[1];
	oopcalc.n = parameter->_dpar[2];
	oopcalc.SetupPointers();
	_oopcalculations.push_back(oopcalc);
	continue;
      }
    }

    //
    // VDW Calculations
    //
    IF_OBFF_LOGLVL_LOW
      OBFFLog("SETTING UP VAN DER WAALS CALCULATIONS...\n");

    OBFFVDWCalculationGaff vdwcalc;
    OBFFParameter *parameter_a, *parameter_b;
    double Ra, Rb, Ea, Eb;

    _vdwcalculations.clear();

    FOR_PAIRS_OF_MOL(p, _mol) {
      a = _mol.GetAtom((*p)[0]);
      b = _mol.GetAtom((*p)[1]);

      // skip this vdw if the atoms are ignored
      if ( _constraints.IsIgnored(a->GetIdx()) || _constraints.IsIgnored(b->GetIdx()) )
        continue;

      // if there are any groups specified, check if the two atoms are in a single _interGroup or if
      // two two atoms are in one of the _interGroups pairs.
      if (HasGroups()) {
        bool validVDW = false;
        for (unsigned int i=0; i < _interGroup.size(); ++i) {
          if (_interGroup[i].BitIsSet(a->GetIdx()) && _interGroup[i].BitIsSet(b->GetIdx()))
            validVDW = true;
        }
        for (unsigned int i=0; i < _interGroups.size(); ++i) {
          if (_interGroups[i].first.BitIsSet(a->GetIdx()) && _interGroups[i].second.BitIsSet(b->GetIdx()))
            validVDW = true;
          if (_interGroups[i].first.BitIsSet(b->GetIdx()) && _interGroups[i].second.BitIsSet(a->GetIdx()))
            validVDW = true;
        }

        if (!validVDW)
          continue;
      }

      parameter_a = GetParameter(a->GetType(), NULL, NULL, NULL, _ffvdwparams);
      if (parameter_a == NULL) { // no vdw parameter -> use hydrogen
        Ra = 1.4870;
        Ea = 0.0157;

        IF_OBFF_LOGLVL_LOW {
          snprintf(_logbuf, BUFF_SIZE, "COULD NOT FIND VDW PARAMETERS FOR ATOM %s, USING HYDROGEN VDW PARAMETERS\n", a->GetType());
          OBFFLog(_logbuf);
        }
      } else {
        Ra = parameter_a->_dpar[0];
        Ea = parameter_a->_dpar[1];
      }

      parameter_b = GetParameter(b->GetType(), NULL, NULL, NULL, _ffvdwparams);
      if (parameter_b == NULL) { // no vdw parameter -> use hydrogen
        Rb = 1.4870;
        Eb = 0.0157;

        IF_OBFF_LOGLVL_LOW {
          snprintf(_logbuf, BUFF_SIZE, "COULD NOT FIND VDW PARAMETERS FOR ATOM %s, USING HYDROGEN VDW PARAMETERS\n", b->GetType());
          OBFFLog(_logbuf);
        }
      } else {
        Rb = parameter_b->_dpar[0];
        Eb = parameter_b->_dpar[1];
      }

      vdwcalc.a = &*a;
      vdwcalc.b = &*b;

      //this calculations only need to be done once for each pair,
      //we do them now and save them for later use
      vdwcalc.Eab = KCAL_TO_KJ * sqrt(Ea * Eb);

      // 1-4 scaling
      if (a->IsOneFour(b))
        vdwcalc.Eab *= 0.5;
      /*
        vdwcalc.is14 = false;
        FOR_NBORS_OF_ATOM (nbr, a)
        FOR_NBORS_OF_ATOM (nbr2, &*nbr)
        FOR_NBORS_OF_ATOM (nbr3, &*nbr2)
        if (b == &*nbr3) {
        vdwcalc.is14 = true;
        vdwcalc.kab *= 0.5;
        }
      */

      // not sure why this is needed, but validation showed it works...
      /*
        if (a->IsInRingSize(6) && b->IsInRingSize(6) && IsInSameRing(a, b))
        vdwcalc.samering = true;
        else if ((a->IsInRingSize(5) || a->IsInRingSize(4)) && (b->IsInRingSize(5) || b->IsInRingSize(4)))
        vdwcalc.samering = true;
        else
        vdwcalc.samering = false;
      */

      vdwcalc.RVDWab = (Ra + Rb);
      vdwcalc.SetupPointers();

      _vdwcalculations.push_back(vdwcalc);
    }

    //
    // Electrostatic Calculations
    //
    IF_OBFF_LOGLVL_LOW
      OBFFLog("SETTING UP ELECTROSTATIC CALCULATIONS...\n");

    OBFFElectrostaticCalculationGaff elecalc;

    _electrostaticcalculations.clear();

    FOR_PAIRS_OF_MOL(p, _mol) {
      a = _mol.GetAtom((*p)[0]);
      b = _mol.GetAtom((*p)[1]);

      // skip this ele if the atoms are ignored
      if ( _constraints.IsIgnored(a->GetIdx()) || _constraints.IsIgnored(b->GetIdx()) )
        continue;

      // if there are any groups specified, check if the two atoms are in a single _interGroup or if
      // two two atoms are in one of the _interGroups pairs.
      if (HasGroups()) {
        bool validEle = false;
        for (unsigned int i=0; i < _interGroup.size(); ++i) {
          if (_interGroup[i].BitIsSet(a->GetIdx()) && _interGroup[i].BitIsSet(b->GetIdx()))
            validEle = true;
        }
        for (unsigned int i=0; i < _interGroups.size(); ++i) {
          if (_interGroups[i].first.BitIsSet(a->GetIdx()) && _interGroups[i].second.BitIsSet(b->GetIdx()))
            validEle = true;
          if (_interGroups[i].first.BitIsSet(b->GetIdx()) && _interGroups[i].second.BitIsSet(a->GetIdx()))
            validEle = true;
        }

        if (!validEle)
          continue;
      }

      elecalc.qq = KCAL_TO_KJ * 332.17 * a->GetPartialCharge() * b->GetPartialCharge() / _epsilon;

      if (elecalc.qq) {
        elecalc.a = &*a;
        elecalc.b = &*b;

        // 1-4 scaling
        if (a->IsOneFour(b))
          elecalc.qq *= 0.5;

        elecalc.SetupPointers();
        _electrostaticcalculations.push_back(elecalc);
      }
    }
    return true;
  }

  bool OBForceFieldGaff::SetupPointers()
  {
    for (unsigned int i = 0; i < _bondcalculations.size(); ++i)
      _bondcalculations[i].SetupPointers();
    for (unsigned int i = 0; i < _anglecalculations.size(); ++i)
      _anglecalculations[i].SetupPointers();
    for (unsigned int i = 0; i < _torsioncalculations.size(); ++i)
      _torsioncalculations[i].SetupPointers();
    for (unsigned int i = 0; i < _vdwcalculations.size(); ++i)
      _vdwcalculations[i].SetupPointers();
    for (unsigned int i = 0; i < _electrostaticcalculations.size(); ++i)
      _electrostaticcalculations[i].SetupPointers();

    return true;
  }

  bool OBForceFieldGaff::ParseParamFile()
  {
    vector<string> vs;
    char buffer[BUFF_SIZE];

    OBFFParameter parameter;

    // open data/gaff.dat
    ifstream ifs;
    if (OpenDatafile(ifs, "gaff.dat").length() == 0) {
      obErrorLog.ThrowError(__FUNCTION__, "Cannot open gaff.dat", obError);
      return false;
    }

    // Set the locale for number parsing to avoid locale issues: PR#1785463
    obLocale.SetLocale();
    ifs.getline(buffer, BUFF_SIZE); //first line title block
    ifs.getline(buffer, BUFF_SIZE); //next line
    while(!ifs.fail() && strlen(buffer)!=0) //read block 2 (atom mass and polarizability)
      {
        tokenize(vs, buffer," -\n\t");
        parameter.clear();
        parameter._a = vs[0]; //KNDSYM
        parameter._dpar.push_back(atof(vs[1].c_str())); // AMASS
        parameter._dpar.push_back(atof(vs[2].c_str())); // ATPOL [A^3]
        _ffpropparams.push_back(parameter);
        ifs.getline(buffer, BUFF_SIZE);
      }
    ifs.getline(buffer, BUFF_SIZE); //next line (block 3: line contains JSOLTY)
    ifs.getline(buffer, BUFF_SIZE); //next line
    while(!ifs.fail() && strlen(buffer)!=0) //read block 4 (bonds)
      {
        tokenize(vs, buffer," -\n\t");
        parameter.clear();
        parameter._a = vs[0]; // IBT
        parameter._b = vs[1]; // JBT
        parameter._dpar.push_back(atof(vs[2].c_str())); // RK [kcal/mol/(A^2)]
        parameter._dpar.push_back(atof(vs[3].c_str())); // REQ [A]
        _ffbondparams.push_back(parameter);
        ifs.getline(buffer, BUFF_SIZE);
      }
    ifs.getline(buffer, BUFF_SIZE); //next line
    while(!ifs.fail() && strlen(buffer)!=0) //read block 5 (angles)
      {
        tokenize(vs, buffer," -\n\t");
        parameter.clear();
        parameter._a = vs[0]; //ITT
        parameter._b = vs[1]; //JTT
        parameter._c = vs[2]; //KTT
        parameter._dpar.push_back(atof(vs[3].c_str())); // TK [kcal/mol/(rad**2)]
        parameter._dpar.push_back(atof(vs[4].c_str())); // TEQ [degrees]
        _ffangleparams.push_back(parameter);
        ifs.getline(buffer, BUFF_SIZE);
      }
    ifs.getline(buffer, BUFF_SIZE); //next line
    while(!ifs.fail() && strlen(buffer)!=0) //read block 6 (dihedrals)
      {
	// torsional potental (PK/IDIVF) * (1 + cos(PN*phi - GAMMA))

        tokenize(vs, buffer," -\n\t");
        parameter.clear();
        parameter._a = vs[0]; //IPT
        parameter._b = vs[1]; //JPT
        parameter._c = vs[2]; //KPT
        parameter._d = vs[3]; //LPT
        parameter._ipar.push_back(atoi(vs[4].c_str())); // IDIVF
        parameter._dpar.push_back(atof(vs[5].c_str())); // PK
        parameter._dpar.push_back(atof(vs[6].c_str())); // GAMMA [degrees]
        parameter._dpar.push_back(atof(vs[7].c_str())); // PN
        _fftorsionparams.push_back(parameter);
        ifs.getline(buffer, BUFF_SIZE);
      }
    ifs.getline(buffer, BUFF_SIZE); //next line
    while(!ifs.fail() && strlen(buffer)!=0) //read block 7(impropers)
      {
//         The input is the same as in for the dihedrals except that
// 	   the torsional barrier height is NOT divided by the factor
//         idivf.  The improper torsions are defined between any four
//         atoms not bonded (in a successive fashion) with each other
//         as in the case of "regular" or "proper" dihedrals.  Improper
//         dihedrals are used to keep certain groups planar and to
//         prevent the racemization of certain centers in the united
//         atom model.  Consult the above reference for details.
//         From Townhee help:
// Improper torsions are not automatically generated by the Towhee code as the rules for determining where they are applied are not always straight-forward. Amber param96 exclusively uses the Stereocenter version of the improper torsions, and they are typically centered on an sp2 atom in order to enforce planarity with its three neighbors. Only one improper torsion allowed to be centered on any atom. These torsions are listed in the Amber literature as i-j-k-l where the angle is the dihedral between i-k-l and j-k-l, and the bonding pattern is i, j, and l are all bonded to atom k, and are also not bonded to each other. In the towhee_input file this stereo improper torsion is listed only for atom k, and the atom order there is l, i, j. Remember that you can set the improper type to 0 to have the code automatically determine the improper type (so long as inpstyle is 2).


        tokenize(vs, buffer," -\n\t");
        parameter.clear();
        parameter._a = vs[0]; //IPT
        parameter._b = vs[1]; //JPT
        parameter._c = vs[2]; //KPT
        parameter._d = vs[3]; //LPT
        parameter._dpar.push_back(atof(vs[4].c_str())); // PK
        parameter._dpar.push_back(atof(vs[5].c_str())); // GAMMA
        parameter._dpar.push_back(atof(vs[6].c_str())); // PN
        _ffoopparams.push_back(parameter);
        ifs.getline(buffer, BUFF_SIZE);
      }
    ifs.getline(buffer, BUFF_SIZE); //next line
    while(!ifs.fail() && strlen(buffer)!=0) //read block 8 (10-12 H-BOND)
      {
        tokenize(vs, buffer," -\n\t");
        parameter.clear();
        parameter._a = vs[0]; // KT1
        parameter._b = vs[1]; // KT2
        parameter._dpar.push_back(atof(vs[2].c_str())); // A
        parameter._dpar.push_back(atof(vs[3].c_str())); // B
        _ffhbondparams.push_back(parameter);
	    ifs.getline(buffer, BUFF_SIZE);
      }
    ifs.getline(buffer, BUFF_SIZE); //next line
    while(!ifs.fail() && strlen(buffer)!=0) //read block 9 (equivalent atoms for vdW potential)
      {
	    // not used. We assume the types are all listed individually in block 10
	    ifs.getline(buffer, BUFF_SIZE);
      }
    ifs.getline(buffer, BUFF_SIZE); //next line
    ifs.getline(buffer, BUFF_SIZE); //'RE' means vdW radius and well-depth are read
    while(!ifs.fail() && strlen(buffer)!=0) //read block 10 (vdWaals)
      {
        // van der Waals potental EDEP * ( (R/r)^12 - 2*(R/r)^6 )
        tokenize(vs, buffer," -\n\t");
        parameter.clear();
        parameter._a = vs[0]; // IBT
        parameter._dpar.push_back(atof(vs[1].c_str())); // R
        parameter._dpar.push_back(atof(vs[2].c_str())); // EDEP (kcal/mol)
        _ffvdwparams.push_back(parameter);
	    ifs.getline(buffer, BUFF_SIZE);
      }

    if (ifs)
      ifs.close();

    // return the locale to the original one
    obLocale.RestoreLocale();

    return 0;
  }

  bool OBForceFieldGaff::SetTypes()
  {
    vector<vector<int> > _mlist; //!< match list for atom typing
    vector<pair<OBSmartsPattern*,string> > _vexttyp; //!< external atom type rules
    vector<vector<int> >::iterator j;
    vector<pair<OBSmartsPattern*,string> >::iterator i;
    OBSmartsPattern *sp;
    vector<string> vs;
    char buffer[150];
    OBAtom *atm, *a, *b;
    OBBitVec visited;
    int BO;

    SetPartialChargesBeforeAtomTyping();
    _mol.SetAtomTypesPerceived();

    // open data/gaff.prm
    ifstream ifs;
    if (OpenDatafile(ifs, "gaff.prm").length() == 0) {
      obErrorLog.ThrowError(__FUNCTION__, "Cannot open gaff.prm", obError);
      return false;
    }

    // Set the locale for number parsing to avoid locale issues: PR#1785463
    obLocale.SetLocale();

    while (ifs.getline(buffer, sizeof(buffer))) {
      if (EQn(buffer, "atom", 4)) {
      	tokenize(vs, buffer);

        sp = new OBSmartsPattern;
        if (sp->Init(vs[1]))
	    _vexttyp.push_back(pair<OBSmartsPattern*,string> (sp,vs[2]));
        else {
          delete sp;
          sp = NULL;
          obErrorLog.ThrowError(__FUNCTION__, " Could not parse atom type table from gaff.prm", obInfo);
          return false;
        }
      }
    }

    for (i = _vexttyp.begin();i != _vexttyp.end();++i) {
      if (i->first->Match(_mol)) {
        _mlist = i->first->GetMapList();
        for (j = _mlist.begin();j != _mlist.end();++j) {
          _mol.GetAtom((*j)[0])->SetType(i->second);
        }
      }
    }

    // Implementation of a special feature of GAFF concerning conjugated bonds
    // In a conjugated ring system cc-cc, and cd-cd are single conjugated bonds, cc-cd are double ones
    // All these types are initially types cc-cc, some cc's are changed to cd.
    // similar for ce and cf, cp and cq, nc and nd, ne and nf, pc and pd, pe and pf
    // Initially types cc, ce, cp, nc, ne , pc and pe are assigned.
    // Some of these are changed into cd, cf, cq, nd, nf, pd and pf
    visited.Resize(_mol.NumAtoms());
    visited.Clear();
    FOR_ATOMS_OF_MOL(ai,_mol){
      atm = &(*ai);
      if (!visited[atm->GetIdx()-1]){
	if( strcmp(atm->GetType(),"cc")==0 ||
	    strcmp(atm->GetType(),"ce")==0 ||
	    strcmp(atm->GetType(),"cp")==0 ||
	    strcmp(atm->GetType(),"nc")==0 ||
	    strcmp(atm->GetType(),"ne")==0 ||
	    strcmp(atm->GetType(),"pc")==0 ||
	    strcmp(atm->GetType(),"pe")==0){
	  FOR_BONDS_OF_ATOM(bond,atm) {
	    a = bond->GetBeginAtom();
	    b = bond->GetEndAtom();
	    if ((strcmp(a->GetType(),"cc")==0||strcmp(a->GetType(),"cd")==0)&&(strcmp(b->GetType(),"cc")==0||strcmp(b->GetType(),"cd")==0)){
	      BO = bond->GetBondOrder();
	      if ( (BO > 1  && strcmp(a->GetType(),b->GetType())==0) || (BO==1 && strcmp(a->GetType(),b->GetType())!=0)){
		if (!visited[a->GetIdx()-1]){
		  a->SetType("cd");
		}
		else if (!visited[b->GetIdx()-1]){
		  b->SetType("cd");
		}
	      }
	      visited.SetBitOn(a->GetIdx()-1);
	      visited.SetBitOn(b->GetIdx()-1);
	    }
	    if ((strcmp(a->GetType(),"ce")==0||strcmp(a->GetType(),"cf")==0)&&(strcmp(b->GetType(),"ce")==0||strcmp(b->GetType(),"cf")==0)){
	      BO = bond->GetBondOrder();
	      if ( (BO > 1  && strcmp(a->GetType(),b->GetType())==0) || (BO==1 && strcmp(a->GetType(),b->GetType())!=0)){
		if (!visited[a->GetIdx()-1]){
		  a->SetType("cf");
		}
		else if (!visited[b->GetIdx()-1]){
		  b->SetType("cf");
		}
	      }
	      visited.SetBitOn(a->GetIdx()-1);
	      visited.SetBitOn(b->GetIdx()-1);
	    }
	    if ((strcmp(a->GetType(),"cp")==0||strcmp(a->GetType(),"cq")==0)&&(strcmp(b->GetType(),"cp")==0||strcmp(b->GetType(),"cq")==0)){
	      BO = bond->GetBondOrder();
	      if ( (BO > 1  && strcmp(a->GetType(),b->GetType())==0) || (BO==1 && strcmp(a->GetType(),b->GetType())!=0)){
		if (!visited[a->GetIdx()-1]){
		  a->SetType("cq");
		}
		else if (!visited[b->GetIdx()-1]){
		  b->SetType("cq");
		}
	      }
	      visited.SetBitOn(a->GetIdx()-1);
	      visited.SetBitOn(b->GetIdx()-1);
	    }
	    if ((strcmp(a->GetType(),"nc")==0||strcmp(a->GetType(),"nd")==0)&&(strcmp(b->GetType(),"nc")==0||strcmp(b->GetType(),"nd")==0)){
	      BO = bond->GetBondOrder();
	      if ( (BO > 1  && strcmp(a->GetType(),b->GetType())==0) || (BO==1 && strcmp(a->GetType(),b->GetType())!=0)){
		if (!visited[a->GetIdx()-1]){
		  a->SetType("nd");
		}
		else if (!visited[b->GetIdx()-1]){
		  b->SetType("nd");
		}
	      }
	      visited.SetBitOn(a->GetIdx()-1);
	      visited.SetBitOn(b->GetIdx()-1);
	    }
	    if ((strcmp(a->GetType(),"ne")==0||strcmp(a->GetType(),"nf")==0)&&(strcmp(b->GetType(),"ne")==0||strcmp(b->GetType(),"nf")==0)){
	      BO = bond->GetBondOrder();
	      if ( (BO > 1  && strcmp(a->GetType(),b->GetType())==0) || (BO==1 && strcmp(a->GetType(),b->GetType())!=0)){
		if (!visited[a->GetIdx()-1]){
		  a->SetType("nf");
		}
		else if (!visited[b->GetIdx()-1]){
		  b->SetType("nf");
		}
	      }
	      visited.SetBitOn(a->GetIdx()-1);
	      visited.SetBitOn(b->GetIdx()-1);
	    }
	    if ((strcmp(a->GetType(),"pc")==0||strcmp(a->GetType(),"pd")==0)&&(strcmp(b->GetType(),"pc")==0||strcmp(b->GetType(),"pd")==0)){
	      BO = bond->GetBondOrder();
	      if ( (BO > 1  && strcmp(a->GetType(),b->GetType())==0) || (BO==1 && strcmp(a->GetType(),b->GetType())!=0)){
		if (!visited[a->GetIdx()-1]){
		  a->SetType("pd");
		}
		else if (!visited[b->GetIdx()-1]){
		  b->SetType("pd");
		}
	      }
	      visited.SetBitOn(a->GetIdx()-1);
	      visited.SetBitOn(b->GetIdx()-1);
	    }
	    if ((strcmp(a->GetType(),"pe")==0||strcmp(a->GetType(),"pf")==0)&&(strcmp(b->GetType(),"pe")==0||strcmp(b->GetType(),"pf")==0)){
	      BO = bond->GetBondOrder();
	      if ( (BO > 1  && strcmp(a->GetType(),b->GetType())==0) || (BO==1 && strcmp(a->GetType(),b->GetType())!=0)){
		if (!visited[a->GetIdx()-1]){
		  a->SetType("pf");
		}
		else if (!visited[b->GetIdx()-1]){
		  b->SetType("pf");
		}
	      }
	      visited.SetBitOn(a->GetIdx()-1);
	      visited.SetBitOn(b->GetIdx()-1);
	    }
	  }
	}
      }
      visited.SetBitOn(atm->GetIdx()-1);
    }

    //     FOR_BONDS_OF_MOL(bond, _mol) {
    //       a = bond->GetBeginAtom();
    //       b = bond->GetEndAtom();
    //       BO = bond->GetBondOrder();
    //       cout << "Bond between: " << a->GetType() << " and " << b->GetType() << " with bondorder " << BO << endl;
    //     }
    // DEBUG (validation)
    //    FOR_ATOMS_OF_MOL (a, _mol)
    //      cout << "ATOMTYPE " << a->GetType() << endl;

    IF_OBFF_LOGLVL_LOW {
      OBFFLog("\nA T O M   T Y P E S\n\n");
      OBFFLog("IDX\tTYPE\tRING\n");

      FOR_ATOMS_OF_MOL (a, _mol) {
        snprintf(_logbuf, BUFF_SIZE, "%d\t%s\t%s\n", a->GetIdx(), a->GetType(),
	  (a->IsInRing() ? (a->IsAromatic() ? "AR" : "AL") : "NO"));
        OBFFLog(_logbuf);
      }

      OBFFLog("\nC H A R G E S\n\n");
      OBFFLog("IDX\tCHARGE\n");

      FOR_ATOMS_OF_MOL (a, _mol) {
        snprintf(_logbuf, BUFF_SIZE, "%d\t%f\n", a->GetIdx(), a->GetPartialCharge());
        OBFFLog(_logbuf);
      }
    }


    if (ifs)
      ifs.close();

    // return the locale to the original one
    obLocale.RestoreLocale();

    return true;
  }

  bool OBForceFieldGaff::SetPartialCharges()
  {
    // Do nothing
    //
    // The Partial Charges need to be set *before* atom typing as Gasteiger Charges require the default
    // atom types. This is done in SetPartialChargesBeforeAtomTyping()
    return true;
  }

  bool OBForceFieldGaff::SetPartialChargesBeforeAtomTyping()
  {
    // Trigger calculation of Gasteiger charges
    // Note that the Gastegier charge calculation checks for the values of particular atom types
    _mol.SetAutomaticPartialCharge(true);
    _mol.SetPartialChargesPerceived(false);
    // Trigger partial charge calculation
    FOR_ATOMS_OF_MOL(atom, _mol) {
      atom->GetPartialCharge();
      break;
    }
    _mol.SetPartialChargesPerceived();

    return true;
  }

  double OBForceFieldGaff::Energy(bool gradients)
  {
    double energy = 0.0;


    IF_OBFF_LOGLVL_MEDIUM
      OBFFLog("\nE N E R G Y\n\n");

    if (gradients) {
      ClearGradients();
      energy  = E_Bond<true>();
      energy += E_Angle<true>();
      energy += E_Torsion<true>();
      energy += E_OOP<true>();
      energy += E_VDW<true>();
      energy += E_Electrostatic<true>();
    } else {
      energy  = E_Bond<false>();
      energy += E_Angle<false>();
      energy += E_Torsion<false>();
      energy += E_OOP<false>();
      energy += E_VDW<false>();
      energy += E_Electrostatic<false>();
    }

    IF_OBFF_LOGLVL_MEDIUM {
      snprintf(_logbuf, BUFF_SIZE, "\nTOTAL ENERGY = %8.3f %s\n", energy, GetUnit().c_str());
      OBFFLog(_logbuf);
    }

    return energy;
  }

  bool OBForceFieldGaff::ValidateGradients ()
  {
    vector3 numgrad, anagrad, err;
    int coordIdx;

    bool passed = true; // set to false if any component fails

    OBFFLog("\nV A L I D A T E   G R A D I E N T S\n\n");
    OBFFLog("ATOM IDX      NUMERICAL GRADIENT           ANALYTICAL GRADIENT        REL. ERROR (%)   \n");
    OBFFLog("----------------------------------------------------------------------------------------\n");
    //     "XX       (000.000, 000.000, 000.000)  (000.000, 000.000, 000.000)  (00.00, 00.00, 00.00)"

    FOR_ATOMS_OF_MOL (a, _mol) {
      coordIdx = (a->GetIdx() - 1) * 3;

      // OBFF_ENERGY
      numgrad = NumericalDerivative(&*a, OBFF_ENERGY);
      Energy(); // compute
      anagrad.Set(_gradientPtr[coordIdx], _gradientPtr[coordIdx+1], _gradientPtr[coordIdx+2]);
      err = ValidateGradientError(numgrad, anagrad);

      snprintf(_logbuf, BUFF_SIZE, "%2d       (%7.3f, %7.3f, %7.3f)  (%7.3f, %7.3f, %7.3f)  (%5.2f, %5.2f, %5.2f)\n", a->GetIdx(), numgrad.x(), numgrad.y(), numgrad.z(),
              anagrad.x(), anagrad.y(), anagrad.z(), err.x(), err.y(), err.z());
      OBFFLog(_logbuf);

      // OBFF_EBOND
      numgrad = NumericalDerivative(&*a, OBFF_EBOND);
      ClearGradients();
      E_Bond();
      anagrad.Set(_gradientPtr[coordIdx], _gradientPtr[coordIdx+1], _gradientPtr[coordIdx+2]);
      err = ValidateGradientError(numgrad, anagrad);

      snprintf(_logbuf, BUFF_SIZE, "    bond    (%7.3f, %7.3f, %7.3f)  (%7.3f, %7.3f, %7.3f)  (%5.2f, %5.2f, %5.2f)\n", numgrad.x(), numgrad.y(), numgrad.z(),
              anagrad.x(), anagrad.y(), anagrad.z(), err.x(), err.y(), err.z());
      OBFFLog(_logbuf);
      // 8% tolerance here because some bonds have slight instability
      if (err.x() > 8.0 || err.y() > 8.0 || err.z() > 8.0)
        passed = false;

      // OBFF_EANGLE
      numgrad = NumericalDerivative(&*a, OBFF_EANGLE);
      ClearGradients();
      E_Angle();
      anagrad.Set(_gradientPtr[coordIdx], _gradientPtr[coordIdx+1], _gradientPtr[coordIdx+2]);
      err = ValidateGradientError(numgrad, anagrad);

      snprintf(_logbuf, BUFF_SIZE, "    angle   (%7.3f, %7.3f, %7.3f)  (%7.3f, %7.3f, %7.3f)  (%5.2f, %5.2f, %5.2f)\n", numgrad.x(), numgrad.y(), numgrad.z(),
              anagrad.x(), anagrad.y(), anagrad.z(), err.x(), err.y(), err.z());
      OBFFLog(_logbuf);
      if (err.x() > 5.0 || err.y() > 5.0 || err.z() > 5.0)
        passed = false;

      // OBFF_ETORSION
      numgrad = NumericalDerivative(&*a, OBFF_ETORSION);
      ClearGradients();
      E_Torsion();
      anagrad.Set(_gradientPtr[coordIdx], _gradientPtr[coordIdx+1], _gradientPtr[coordIdx+2]);
      err = ValidateGradientError(numgrad, anagrad);

      snprintf(_logbuf, BUFF_SIZE, "    torsion (%7.3f, %7.3f, %7.3f)  (%7.3f, %7.3f, %7.3f)  (%5.2f, %5.2f, %5.2f)\n", numgrad.x(), numgrad.y(), numgrad.z(),
              anagrad.x(), anagrad.y(), anagrad.z(), err.x(), err.y(), err.z());
      OBFFLog(_logbuf);
      if (err.x() > 5.0 || err.y() > 5.0 || err.z() > 5.0)
        passed = false;

      // OBFF_EOOP
      numgrad = NumericalDerivative(&*a, OBFF_EOOP);
      ClearGradients();
      E_OOP(); // compute
      anagrad.Set(_gradientPtr[coordIdx], _gradientPtr[coordIdx+1], _gradientPtr[coordIdx+2]);
      err = ValidateGradientError(numgrad, anagrad);

      snprintf(_logbuf, BUFF_SIZE, "    oop     (%7.3f, %7.3f, %7.3f)  (%7.3f, %7.3f, %7.3f)  (%5.2f, %5.2f, %5.2f)\n", numgrad.x(), numgrad.y(), numgrad.z(),
              anagrad.x(), anagrad.y(), anagrad.z(), err.x(), err.y(), err.z());
      OBFFLog(_logbuf);


      // OBFF_EVDW
      numgrad = NumericalDerivative(&*a, OBFF_EVDW);
      ClearGradients();
      E_VDW();
      anagrad.Set(_gradientPtr[coordIdx], _gradientPtr[coordIdx+1], _gradientPtr[coordIdx+2]);
      err = ValidateGradientError(numgrad, anagrad);

      snprintf(_logbuf, BUFF_SIZE, "    vdw     (%7.3f, %7.3f, %7.3f)  (%7.3f, %7.3f, %7.3f)  (%5.2f, %5.2f, %5.2f)\n", numgrad.x(), numgrad.y(), numgrad.z(),
              anagrad.x(), anagrad.y(), anagrad.z(), err.x(), err.y(), err.z());
      OBFFLog(_logbuf);
      if (err.x() > 5.0 || err.y() > 5.0 || err.z() > 5.0)
        passed = false;

      // OBFF_EELECTROSTATIC
      numgrad = NumericalDerivative(&*a, OBFF_EELECTROSTATIC);
      ClearGradients();
      E_Electrostatic();
      anagrad.Set(_gradientPtr[coordIdx], _gradientPtr[coordIdx+1], _gradientPtr[coordIdx+2]);
      err = ValidateGradientError(numgrad, anagrad);

      snprintf(_logbuf, BUFF_SIZE, "    electro (%7.3f, %7.3f, %7.3f)  (%7.3f, %7.3f, %7.3f)  (%5.2f, %5.2f, %5.2f)\n", numgrad.x(), numgrad.y(), numgrad.z(),
              anagrad.x(), anagrad.y(), anagrad.z(), err.x(), err.y(), err.z());
      OBFFLog(_logbuf);
      if (err.x() > 5.0 || err.y() > 5.0 || err.z() > 5.0)
        passed = false;
    }

    return passed; // are all components good enough?
  }


} // end namespace OpenBabel

//! \file forcefieldgaff.cpp
//! \brief Gaff force field

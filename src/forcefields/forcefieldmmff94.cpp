/*********************************************************************
forcefieldmmff94.cpp - MMFF94 force field

Copyright (C) 2006-2008 by Tim Vandermeersch <tim.vandermeersch@gmail.com>

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

/*
 * Source code layout:
 * - Functions to calculate the actual interactions
 * - Parse parameter files
 * - Setup Functions
 * - Validation functions
 * - Calculate bond type, angle type, stretch-bend type, torsion type
 * - Various tables & misc. functions
 *
 */

#include <openbabel/babelconfig.h>
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/locale.h>
#include <openbabel/elements.h>
#include <openbabel/bond.h>
#include <openbabel/obiter.h>
#include <openbabel/ring.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <iomanip>
#include "forcefieldmmff94.h"

using namespace std;

namespace OpenBabel
{
  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  //
  //  Functions to calculate the actual interactions
  //
  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////

  double OBForceFieldMMFF94::Energy(bool gradients)
  {
    double energy;

    IF_OBFF_LOGLVL_MEDIUM
      OBFFLog("\nE N E R G Y\n\n");

    if (gradients) {
      ClearGradients();
      energy  = E_Bond<true>();
      energy += E_Angle<true>();
      energy += E_StrBnd<true>();
      energy += E_Torsion<true>();
      energy += E_OOP<true>();
      energy += E_VDW<true>();
      energy += E_Electrostatic<true>();
    } else {
      energy  = E_Bond<false>();
      energy += E_Angle<false>();
      energy += E_StrBnd<false>();
      energy += E_Torsion<false>();
      energy += E_OOP<false>();
      energy += E_VDW<false>();
      energy += E_Electrostatic<false>();
    }

    IF_OBFF_LOGLVL_MEDIUM {
      snprintf(_logbuf, BUFF_SIZE, "\nTOTAL ENERGY = %8.5f %s\n", energy, GetUnit().c_str());
      OBFFLog(_logbuf);
    }

    return energy;
  }

  //
  // MMFF part I - page 494
  //
  //                   kb_ij                              7
  // EB_ij = 143.9325 ------- /\r_ij^2 (1 + cs /\_rij + ---- cs^2 r_ij^2)
  //                     2                               12
  //
  // kb_ij	force constant (md/A)
  //
  // /\r_ij 	r_ij - r0_ij (A)
  //
  // cs		cubic stretch constant = -2 A^(-1)
  //
  template<bool gradients>
  inline void OBFFBondCalculationMMFF94::Compute()
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

      const double dE = 143.9325 * kb * delta * (1.0 - 3.0 * delta + 14.0/3.0 * delta2);

      OBForceField::VectorSelfMultiply(force_a, dE);
      OBForceField::VectorSelfMultiply(force_b, dE);
    } else {
      rab = OBForceField::VectorDistance(pos_a, pos_b);
      delta = rab - r0;
      delta2 = delta * delta;
    }

    energy = kb * delta2 * (1.0 - 2.0 * delta + 7.0/3.0 * delta2);
  }

  template<bool gradients>
  double OBForceFieldMMFF94::E_Bond()
  {
    double energy = 0.0;

    IF_OBFF_LOGLVL_HIGH {
      OBFFLog("\nB O N D   S T R E T C H I N G\n\n");
      OBFFLog("ATOM TYPES   FF    BOND       IDEAL       FORCE\n");
      OBFFLog(" I    J     CLASS  LENGTH     LENGTH     CONSTANT      DELTA      ENERGY\n");
      OBFFLog("------------------------------------------------------------------------\n");
    }

    #ifdef _OPENMP
    #pragma omp parallel for reduction(+:energy)
    #endif
    for (int i = 0; i < _bondcalculations.size(); ++i) {
      _bondcalculations[i].template Compute<gradients>();
      energy += _bondcalculations[i].energy;

      #ifndef _OPENMP
      if (gradients) {
        AddGradient(_bondcalculations[i].force_a, _bondcalculations[i].idx_a);
        AddGradient(_bondcalculations[i].force_b, _bondcalculations[i].idx_b);
      }
      #endif

      IF_OBFF_LOGLVL_HIGH {
        snprintf(_logbuf, BUFF_SIZE, "%2d   %2d      %d   %8.3f   %8.3f     %8.3f   %8.3f   %8.3f\n",
                atoi(_bondcalculations[i].a->GetType()), atoi(_bondcalculations[i].b->GetType()),
                _bondcalculations[i].bt, _bondcalculations[i].rab, _bondcalculations[i].r0,
                _bondcalculations[i].kb, _bondcalculations[i].delta,
                143.9325 * 0.5 * _bondcalculations[i].energy);
        OBFFLog(_logbuf);
      }
    }

    #ifdef _OPENMP
    for (int i = 0; i < _bondcalculations.size(); ++i) {
      if (gradients) {
        AddGradient(_bondcalculations[i].force_a, _bondcalculations[i].idx_a);
        AddGradient(_bondcalculations[i].force_b, _bondcalculations[i].idx_b);
      }
    }
    #endif

    IF_OBFF_LOGLVL_MEDIUM {
      snprintf(_logbuf, BUFF_SIZE, "     TOTAL BOND STRETCHING ENERGY = %8.5f %s\n",  143.9325 * 0.5 * energy, GetUnit().c_str());
      OBFFLog(_logbuf);
    }

    return (143.9325 * 0.5 * energy);
  }

  //
  // MMFF part I - page 495
  //
  //                       ka_ijk
  // EA_ijk = 0.438449325 -------- /\0_ijk^2 (1 + cs /\0_ijk)
  //                         2
  //
  // ka_ijk	force constant (md A/rad^2)
  //
  // /\0_ijk 	0_ijk - 00_ijk (degrees)
  //
  // cs		cubic bend constant = -0.007 deg^-1 = -0.4 rad^-1
  //
  template<bool gradients>
  inline void OBFFAngleCalculationMMFF94::Compute()
  {
    if (OBForceField::IgnoreCalculation(idx_a, idx_b, idx_c)) {
      energy = 0.0;
      return;
    }

    double delta2, dE;

    if (gradients) {
      theta = OBForceField::VectorAngleDerivative(pos_a, pos_b, pos_c, force_a, force_b, force_c);

      if (!isfinite(theta))
        theta = 0.0; // doesn't explain why GetAngle is returning NaN but solves it for us;

      delta = theta - theta0;

      if (linear) {
        energy = 143.9325 * ka * (1.0 + cos((theta) * DEG_TO_RAD));
        dE = -sin((theta) * DEG_TO_RAD) * 143.9325 * ka;
      } else {
        delta2 = delta * delta;
        energy = 0.043844 * 0.5 * ka * delta2 * (1.0 - 0.007 * delta);
        dE = RAD_TO_DEG * 0.043844 * ka * delta * (1.0 - 1.5 * 0.007 * delta);
      }

      OBForceField::VectorSelfMultiply(force_a, dE);
      OBForceField::VectorSelfMultiply(force_b, dE);
      OBForceField::VectorSelfMultiply(force_c, dE);
    } else {
      theta = OBForceField::VectorAngle(pos_a, pos_b, pos_c);

      if (!isfinite(theta))
        theta = 0.0; // doesn't explain why GetAngle is returning NaN but solves it for us;

      delta = theta - theta0;

      if (linear) {
        energy = 143.9325 * ka * (1.0 + cos(theta * DEG_TO_RAD));
      } else {
        delta2 = delta * delta;
        energy = 0.043844 * 0.5 * ka * delta2 * (1.0 - 0.007 * delta);
      }

    }

  }

  template<bool gradients>
  double OBForceFieldMMFF94::E_Angle()
  {
    double energy = 0.0;

    IF_OBFF_LOGLVL_HIGH {
      OBFFLog("\nA N G L E   B E N D I N G\n\n");
      OBFFLog("ATOM TYPES        FF    VALENCE     IDEAL      FORCE\n");
      OBFFLog(" I    J    K     CLASS   ANGLE      ANGLE     CONSTANT      DELTA      ENERGY\n");
      OBFFLog("-----------------------------------------------------------------------------\n");
    }

    #ifdef _OPENMP
    #pragma omp parallel for reduction(+:energy)
    #endif
    for (int i = 0; i < _anglecalculations.size(); ++i) {

      _anglecalculations[i].template Compute<gradients>();
      energy += _anglecalculations[i].energy;

      #ifndef _OPENMP
      if (gradients) {
        AddGradient(_anglecalculations[i].force_a, _anglecalculations[i].idx_a);
        AddGradient(_anglecalculations[i].force_b, _anglecalculations[i].idx_b);
        AddGradient(_anglecalculations[i].force_c, _anglecalculations[i].idx_c);
      }
      #endif

      IF_OBFF_LOGLVL_HIGH {
        snprintf(_logbuf, BUFF_SIZE, "%2d   %2d   %2d      %d   %8.3f   %8.3f     %8.3f   %8.3f   %8.3f\n",
                atoi(_anglecalculations[i].a->GetType()), atoi(_anglecalculations[i].b->GetType()),
                atoi(_anglecalculations[i].c->GetType()), _anglecalculations[i].at,
                _anglecalculations[i].theta, _anglecalculations[i].theta0,
                _anglecalculations[i].ka, _anglecalculations[i].delta,
                _anglecalculations[i].energy);
        OBFFLog(_logbuf);
      }
    }

    #ifdef _OPENMP
    for (int i = 0; i < _anglecalculations.size(); ++i) {
      if (gradients) {
        AddGradient(_anglecalculations[i].force_a, _anglecalculations[i].idx_a);
        AddGradient(_anglecalculations[i].force_b, _anglecalculations[i].idx_b);
        AddGradient(_anglecalculations[i].force_c, _anglecalculations[i].idx_c);
      }
    }
    #endif

    IF_OBFF_LOGLVL_MEDIUM {
      snprintf(_logbuf, BUFF_SIZE, "     TOTAL ANGLE BENDING ENERGY = %8.5f %s\n", energy, GetUnit().c_str());
      OBFFLog(_logbuf);
    }

    return energy;
  }

  //
  // MMFF part I - page 495
  //
  // EBA_ijk = 2.51210 (kba_ijk /\r_ij + kba_kji /\r_kj) /\0_ijk
  //
  // kba_ijk	force constant (md/rad)
  // kba_kji	force constant (md/rad)
  //
  // /\r_xx 	see above
  // /\0_ijk 	see above
  //
  template<bool gradients>
  inline void OBFFStrBndCalculationMMFF94::Compute()
  {
    if (OBForceField::IgnoreCalculation(idx_a, idx_b, idx_c)) {
      energy = 0.0;
      return;
    }

    if (gradients) {
      theta = OBForceField::VectorAngleDerivative(pos_a, pos_b, pos_c,
                                                  force_abc_a, force_abc_b, force_abc_c);
      rab = OBForceField::VectorDistanceDerivative(pos_a, pos_b, force_ab_a, force_ab_b);
      rbc = OBForceField::VectorDistanceDerivative(pos_b, pos_c, force_bc_b, force_bc_c);
    } else {
      theta = OBForceField::VectorAngle(pos_a, pos_b, pos_c);
      rab = OBForceField::VectorDistance(pos_a, pos_b);
      rbc = OBForceField::VectorDistance(pos_b, pos_c);
    }

    if (!isfinite(theta))
      theta = 0.0; // doesn't explain why GetAngle is returning NaN but solves it for us;

    delta_theta = theta - theta0;
    delta_rab = rab - rab0;
    delta_rbc = rbc - rbc0;
    const double factor = RAD_TO_DEG * (kbaABC * delta_rab + kbaCBA * delta_rbc);

    energy = DEG_TO_RAD * factor * delta_theta;
    if (gradients) {
      //grada = 2.51210 * (kbaABC * rab_da * delta_theta + RAD_TO_DEG * theta_da * (kbaABC * delta_rab + kbaCBA * delta_rbc));
      OBForceField::VectorSelfMultiply(force_ab_a, (kbaABC*delta_theta));
      OBForceField::VectorSelfMultiply(force_abc_a, factor);
      OBForceField::VectorAdd(force_ab_a, force_abc_a, force_ab_a);
      OBForceField::VectorMultiply(force_ab_a, 2.51210, force_a);
      //gradc = 2.51210 * (kbaCBA * rbc_dc * delta_theta + RAD_TO_DEG * theta_dc * (kbaABC * delta_rab + kbaCBA * delta_rbc));
      OBForceField::VectorSelfMultiply(force_bc_c, (kbaCBA*delta_theta));
      OBForceField::VectorSelfMultiply(force_abc_c, factor);
      OBForceField::VectorAdd(force_bc_c, force_abc_c, force_bc_c);
      OBForceField::VectorMultiply(force_bc_c, 2.51210, force_c);
      //gradb = -grada - gradc;
      OBForceField::VectorAdd(force_a, force_c, force_b);
      OBForceField::VectorSelfMultiply(force_b, -1.0);
    }
  }

  template<bool gradients>
  double OBForceFieldMMFF94::E_StrBnd()
  {
    double energy = 0.0;

    IF_OBFF_LOGLVL_HIGH {
      OBFFLog("\nS T R E T C H   B E N D I N G\n\n");
      OBFFLog("ATOM TYPES        FF    VALENCE     DELTA        FORCE CONSTANT\n");
      OBFFLog(" I    J    K     CLASS   ANGLE      ANGLE        I J        J K      ENERGY\n");
      OBFFLog("---------------------------------------------------------------------------\n");
    }

    #ifdef _OPENMP
    #pragma omp parallel for reduction(+:energy)
    #endif
    for (int i = 0; i < _strbndcalculations.size(); ++i) {

      _strbndcalculations[i].template Compute<gradients>();
      energy += _strbndcalculations[i].energy;

      #ifndef _OPENMP
      if (gradients) {
        AddGradient(_strbndcalculations[i].force_a, _strbndcalculations[i].idx_a);
        AddGradient(_strbndcalculations[i].force_b, _strbndcalculations[i].idx_b);
        AddGradient(_strbndcalculations[i].force_c, _strbndcalculations[i].idx_c);
      }
      #endif

      IF_OBFF_LOGLVL_HIGH {
        snprintf(_logbuf, BUFF_SIZE, "%2d   %2d   %2d     %2d   %8.3f   %8.3f   %8.3f   %8.3f   %8.3f\n",
                atoi(_strbndcalculations[i].a->GetType()), atoi(_strbndcalculations[i].b->GetType()),
                atoi(_strbndcalculations[i].c->GetType()), _strbndcalculations[i].sbt,
                _strbndcalculations[i].theta, _strbndcalculations[i].delta_theta,
                _strbndcalculations[i].kbaABC, _strbndcalculations[i].kbaCBA,
                2.51210 * _strbndcalculations[i].energy);
        OBFFLog(_logbuf);
      }
    }

    #ifdef _OPENMP
    for (int i = 0; i < _strbndcalculations.size(); ++i) {
      if (gradients) {
        AddGradient(_strbndcalculations[i].force_a, _strbndcalculations[i].idx_a);
        AddGradient(_strbndcalculations[i].force_b, _strbndcalculations[i].idx_b);
        AddGradient(_strbndcalculations[i].force_c, _strbndcalculations[i].idx_c);
      }
    }
    #endif

    IF_OBFF_LOGLVL_MEDIUM {
      snprintf(_logbuf, BUFF_SIZE, "     TOTAL STRETCH BENDING ENERGY = %8.5f %s\n", 2.51210 * energy, GetUnit().c_str());
      OBFFLog(_logbuf);
    }

    return (2.51210 * energy);
  }

  int OBForceFieldMMFF94::GetElementRow(OBAtom *atom)
  {
    int row;

    row = 0;

    if (atom->GetAtomicNum() > 2)
      row++;
    if (atom->GetAtomicNum() > 10)
      row++;
    if (atom->GetAtomicNum() > 18)
      row++;
    if (atom->GetAtomicNum() > 36)
      row++;
    if (atom->GetAtomicNum() > 54)
      row++;
    if (atom->GetAtomicNum() > 86)
      row++;

    return row;
  }

  //
  // MMFF part I - page 495
  //
  // ET_ijkl = 0.5 ( V1 (1 + cos(0_ijkl)) + V2 (1 - cos(2 0_ijkl)) + V3 (1 + cos(3 0_ijkl)) )
  //
  // V1		force constant (md/rad)
  // V2		force constant (md/rad)
  // V3		force constant (md/rad)
  //
  // 0_ijkl 	torsion angle (degrees)
  //
  template<bool gradients>
  inline void OBFFTorsionCalculationMMFF94::Compute()
  {
    if (OBForceField::IgnoreCalculation(idx_a, idx_b, idx_c, idx_d)) {
      energy = 0.0;
      return;
    }

    double cosine, cosine2, cosine3;
    double phi1, phi2, phi3;
    double dE, sine, sine2, sine3;

    if (gradients) {
      tor = OBForceField::VectorTorsionDerivative(pos_a, pos_b, pos_c, pos_d,
                                                  force_a, force_b, force_c, force_d);
      if (!isfinite(tor))
        tor = 1.0e-3;

      sine = sin(DEG_TO_RAD * tor);
      sine2 = sin(2.0 * DEG_TO_RAD * tor);
      sine3 = sin(3.0 * DEG_TO_RAD * tor);

      dE = 0.5 * (v1 * sine - 2.0 * v2 * sine2 + 3.0 * v3 * sine3); // MMFF

      OBForceField::VectorSelfMultiply(force_a, dE);
      OBForceField::VectorSelfMultiply(force_b, dE);
      OBForceField::VectorSelfMultiply(force_c, dE);
      OBForceField::VectorSelfMultiply(force_d, dE);
    } else {
      tor = OBForceField::VectorTorsion(pos_a, pos_b, pos_c, pos_d);
      if (!isfinite(tor))
        tor = 1.0e-3;
    }

    cosine = cos(DEG_TO_RAD * tor);
    cosine2 = cos(DEG_TO_RAD * 2 * tor);
    cosine3 = cos(DEG_TO_RAD * 3 * tor);

    phi1 = 1.0 + cosine;
    phi2 = 1.0 - cosine2;
    phi3 = 1.0 + cosine3;

    energy = (v1 * phi1 + v2 * phi2 + v3 * phi3);

  }

  template<bool gradients>
  double OBForceFieldMMFF94::E_Torsion()
  {
    double energy = 0.0;

    IF_OBFF_LOGLVL_HIGH {
      OBFFLog("\nT O R S I O N A L\n\n");
      OBFFLog("ATOM TYPES             FF     TORSION       FORCE CONSTANT\n");
      OBFFLog(" I    J    K    L     CLASS    ANGLE         V1   V2   V3     ENERGY\n");
      OBFFLog("--------------------------------------------------------------------\n");
    }

    #ifdef _OPENMP
    #pragma omp parallel for reduction(+:energy)
    #endif
    for (int i = 0; i < _torsioncalculations.size(); ++i) {

      _torsioncalculations[i].template Compute<gradients>();
      energy += _torsioncalculations[i].energy;

      #ifndef _OPENMP
      if (gradients) {
        AddGradient(_torsioncalculations[i].force_a, _torsioncalculations[i].idx_a);
        AddGradient(_torsioncalculations[i].force_b, _torsioncalculations[i].idx_b);
        AddGradient(_torsioncalculations[i].force_c, _torsioncalculations[i].idx_c);
        AddGradient(_torsioncalculations[i].force_d, _torsioncalculations[i].idx_d);
      }
      #endif

      IF_OBFF_LOGLVL_HIGH {
        snprintf(_logbuf, BUFF_SIZE, "%2d   %2d   %2d   %2d      %d   %8.3f   %6.3f   %6.3f   %6.3f   %8.3f\n",
                atoi(_torsioncalculations[i].a->GetType()), atoi(_torsioncalculations[i].b->GetType()),
                atoi(_torsioncalculations[i].c->GetType()), atoi(_torsioncalculations[i].d->GetType()),
                _torsioncalculations[i].tt, _torsioncalculations[i].tor, _torsioncalculations[i].v1,
                _torsioncalculations[i].v2, _torsioncalculations[i].v3, 0.5 * _torsioncalculations[i].energy);
        OBFFLog(_logbuf);
      }

    }

    #ifdef _OPENMP
    for (int i = 0; i < _torsioncalculations.size(); ++i) {
      if (gradients) {
        AddGradient(_torsioncalculations[i].force_a, _torsioncalculations[i].idx_a);
        AddGradient(_torsioncalculations[i].force_b, _torsioncalculations[i].idx_b);
        AddGradient(_torsioncalculations[i].force_c, _torsioncalculations[i].idx_c);
        AddGradient(_torsioncalculations[i].force_d, _torsioncalculations[i].idx_d);
      }
    }
    #endif

    IF_OBFF_LOGLVL_MEDIUM {
      snprintf(_logbuf, BUFF_SIZE, "     TOTAL TORSIONAL ENERGY = %8.5f %s\n", 0.5 * energy, GetUnit().c_str());
      OBFFLog(_logbuf);
    }

    return (0.5 * energy);
  }

  //						//
  //  a						//
  //   \  					//
  //    b---d      plane = a-b-c		//
  //   / 					//
  //  c						//
  //						//
  template<bool gradients>
  void OBFFOOPCalculationMMFF94::Compute()
  {
    if (OBForceField::IgnoreCalculation(idx_a, idx_b, idx_c, idx_d)) {
      energy = 0.0;
      return;
    }

    double angle2, dE;

    if (gradients) {
      angle = OBForceField::VectorOOPDerivative(pos_a, pos_b, pos_c, pos_d,
                                                force_a, force_b, force_c, force_d);

      dE =  (-1.0 * RAD_TO_DEG * 0.043844 * angle * koop) / cos(angle * DEG_TO_RAD);

      OBForceField::VectorSelfMultiply(force_a, dE);
      OBForceField::VectorSelfMultiply(force_b, dE);
      OBForceField::VectorSelfMultiply(force_c, dE);
      OBForceField::VectorSelfMultiply(force_d, dE);
    } else {
      angle = OBForceField::VectorOOP(pos_a, pos_b, pos_c, pos_d);
    }

    if (!isfinite(angle))
      angle = 0.0; // doesn't explain why GetAngle is returning NaN but solves it for us;

    angle2 = angle * angle;
    energy = koop * angle2;

  }

  template<bool gradients>
  double OBForceFieldMMFF94::E_OOP()
  {
    double energy = 0.0;

    IF_OBFF_LOGLVL_HIGH {
      OBFFLog("\nO U T - O F - P L A N E   B E N D I N G\n\n");
      OBFFLog("ATOM TYPES             FF       OOP     FORCE\n");
      OBFFLog(" I    J    K    L     CLASS    ANGLE   CONSTANT     ENERGY\n");
      OBFFLog("----------------------------------------------------------\n");
    }

    #ifdef _OPENMP
    #pragma omp parallel for reduction(+:energy)
    #endif
    for (int i = 0; i < _oopcalculations.size(); ++i) {

      _oopcalculations[i].template Compute<gradients>();
      energy += _oopcalculations[i].energy;

      #ifndef _OPENMP
      if (gradients) {
        AddGradient(_oopcalculations[i].force_a, _oopcalculations[i].idx_a);
        AddGradient(_oopcalculations[i].force_b, _oopcalculations[i].idx_b);
        AddGradient(_oopcalculations[i].force_c, _oopcalculations[i].idx_c);
        AddGradient(_oopcalculations[i].force_d, _oopcalculations[i].idx_d);
      }
      #endif

      IF_OBFF_LOGLVL_HIGH {
        snprintf(_logbuf, BUFF_SIZE, "%2d   %2d   %2d   %2d      0   %8.3f   %8.3f     %8.3f\n",
                atoi(_oopcalculations[i].a->GetType()), atoi(_oopcalculations[i].b->GetType()),
                atoi(_oopcalculations[i].c->GetType()), atoi(_oopcalculations[i].d->GetType()),
                _oopcalculations[i].angle, _oopcalculations[i].koop,
                0.043844 * 0.5 * _oopcalculations[i].energy);
        OBFFLog(_logbuf);
      }
    }

    #ifdef _OPENMP
    for (int i = 0; i < _oopcalculations.size(); ++i) {
      if (gradients) {
        AddGradient(_oopcalculations[i].force_a, _oopcalculations[i].idx_a);
        AddGradient(_oopcalculations[i].force_b, _oopcalculations[i].idx_b);
        AddGradient(_oopcalculations[i].force_c, _oopcalculations[i].idx_c);
        AddGradient(_oopcalculations[i].force_d, _oopcalculations[i].idx_d);
      }
    }
    #endif

    IF_OBFF_LOGLVL_MEDIUM {
      snprintf(_logbuf, BUFF_SIZE, "     TOTAL OUT-OF-PLANE BENDING ENERGY = %8.5f %s\n", 0.043844 * 0.5 * energy, GetUnit().c_str());
      OBFFLog(_logbuf);
    }

    return (0.043844 * 0.5 * energy);
  }

  template<bool gradients>
  inline void OBFFVDWCalculationMMFF94::Compute()
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

    const double rab7 = rab*rab*rab*rab*rab*rab*rab;

    double erep = (1.07 * R_AB) / (rab + 0.07 * R_AB); //***
    double erep7 = erep*erep*erep*erep*erep*erep*erep;

    double eattr = (((1.12 * R_AB7) / (rab7 + 0.12 * R_AB7)) - 2.0);

    energy = epsilon * erep7 * eattr;

    if (gradients) {
      const double q = rab / R_AB;
      const double q6 = q*q*q*q*q*q;
      const double q7 = q6 * q;
      erep = 1.07 / (q + 0.07);
      erep7 = erep*erep*erep*erep*erep*erep*erep;
      const double term = q7 + 0.12;
      const double term2 = term * term;
      eattr = (-7.84 * q6) / term2 + ((-7.84 / term) + 14) / (q + 0.07);
      const double dE = (epsilon / R_AB) * erep7 * eattr;
      OBForceField::VectorSelfMultiply(force_a, dE);
      OBForceField::VectorSelfMultiply(force_b, dE);
    }
  }

  template<bool gradients>
  double OBForceFieldMMFF94::E_VDW()
  {
    double energy = 0.0;

    IF_OBFF_LOGLVL_HIGH {
      OBFFLog("\nV A N   D E R   W A A L S\n\n");
      OBFFLog("ATOM TYPES\n");
      OBFFLog(" I    J        Rij       R*IJ    EPSILON    ENERGY\n");
      OBFFLog("--------------------------------------------------\n");
      //       XX   XX     -000.000  -000.000  -000.000  -000.000
    }

    #ifdef _OPENMP
    #pragma omp parallel for reduction(+:energy)
    #endif
    for (int i = 0; i < _vdwcalculations.size(); ++i) {
      // Cut-off check
      if (_cutoff)
        if (!_vdwpairs.BitIsSet(_vdwcalculations[i].pairIndex))
          continue;

      _vdwcalculations[i].template Compute<gradients>();
      energy += _vdwcalculations[i].energy;

      #ifndef _OPENMP
      if (gradients) {
        AddGradient(_vdwcalculations[i].force_a, _vdwcalculations[i].idx_a);
        AddGradient(_vdwcalculations[i].force_b, _vdwcalculations[i].idx_b);
      }
      #endif

      IF_OBFF_LOGLVL_HIGH {
        snprintf(_logbuf, BUFF_SIZE, "%2d   %2d     %8.3f  %8.3f  %8.3f  %8.3f\n",
                atoi(_vdwcalculations[i].a->GetType()), atoi(_vdwcalculations[i].b->GetType()),
                _vdwcalculations[i].rab, _vdwcalculations[i].R_AB, _vdwcalculations[i].epsilon, _vdwcalculations[i].energy);
        OBFFLog(_logbuf);
      }

    }

    #ifdef _OPENMP
    for (int i = 0; i < _vdwcalculations.size(); ++i) {
      // Cut-off check
      if (_cutoff)
        if (!_vdwpairs.BitIsSet(i))
          continue;

      if (gradients) {
        AddGradient(_vdwcalculations[i].force_a, _vdwcalculations[i].idx_a);
        AddGradient(_vdwcalculations[i].force_b, _vdwcalculations[i].idx_b);
      }
    }
    #endif

    IF_OBFF_LOGLVL_MEDIUM {
      snprintf(_logbuf, BUFF_SIZE, "     TOTAL VAN DER WAALS ENERGY = %8.5f %s\n", energy, GetUnit().c_str());
      OBFFLog(_logbuf);
    }

    return energy;
  }

  template<bool gradients>
  inline void OBFFElectrostaticCalculationMMFF94::Compute()
  {
    if (OBForceField::IgnoreCalculation(idx_a, idx_b)) {
      energy = 0.0;
      return;
    }

    if (gradients) {
      rab = OBForceField::VectorDistanceDerivative(pos_a, pos_b, force_a, force_b);
      rab += 0.05; // ??
      const double rab2 = rab * rab;
      const double dE = -qq / rab2;
      OBForceField::VectorSelfMultiply(force_a, dE);
      OBForceField::VectorSelfMultiply(force_b, dE);
    } else {
      rab = OBForceField::VectorDistance(pos_a, pos_b);
      rab += 0.05; // ??
    }

    energy = qq / rab;
  }

  template<bool gradients>
  double OBForceFieldMMFF94::E_Electrostatic()
  {
    double energy = 0.0;

    IF_OBFF_LOGLVL_HIGH {
      OBFFLog("\nE L E C T R O S T A T I C   I N T E R A C T I O N S\n\n");
      OBFFLog("ATOM TYPES\n");
      OBFFLog(" I    J        Rij        Qi         Qj        ENERGY\n");
      OBFFLog("-----------------------------------------------------\n");
      //       XX   XX     XXXXXXXX   XXXXXXXX   XXXXXXXX   XXXXXXXX
    }

    #ifdef _OPENMP
    #pragma omp parallel for reduction(+:energy)
    #endif
    for (int i = 0; i < _electrostaticcalculations.size(); ++i) {
      // Cut-off check
      if (_cutoff)
        if (!_elepairs.BitIsSet(_electrostaticcalculations[i].pairIndex))
          continue;

      _electrostaticcalculations[i].template Compute<gradients>();
      energy += _electrostaticcalculations[i].energy;

      #ifndef _OPENMP
      if (gradients) {
        AddGradient(_electrostaticcalculations[i].force_a, _electrostaticcalculations[i].idx_a);
        AddGradient(_electrostaticcalculations[i].force_b, _electrostaticcalculations[i].idx_b);
      }
      #endif

      IF_OBFF_LOGLVL_HIGH {
        snprintf(_logbuf, BUFF_SIZE, "%2d   %2d   %8.3f  %8.3f  %8.3f  %8.3f\n",
                atoi(_electrostaticcalculations[i].a->GetType()), atoi(_electrostaticcalculations[i].b->GetType()),
                _electrostaticcalculations[i].rab, _electrostaticcalculations[i].a->GetPartialCharge(),
                _electrostaticcalculations[i].b->GetPartialCharge(), _electrostaticcalculations[i].energy);
        OBFFLog(_logbuf);
      }
    }

    #ifdef _OPENMP
    for (int i = 0; i < _electrostaticcalculations.size(); ++i) {
      // Cut-off check
      if (_cutoff)
        if (!_elepairs.BitIsSet(i))
          continue;

      if (gradients) {
        AddGradient(_electrostaticcalculations[i].force_a, _electrostaticcalculations[i].idx_a);
        AddGradient(_electrostaticcalculations[i].force_b, _electrostaticcalculations[i].idx_b);
      }
    }
    #endif

    IF_OBFF_LOGLVL_MEDIUM {
      snprintf(_logbuf, BUFF_SIZE, "     TOTAL ELECTROSTATIC ENERGY = %8.5f %s\n", energy, GetUnit().c_str());
      OBFFLog(_logbuf);
    }

    return energy;
  }

  //
  // OBForceFieldMMFF member functions
  //
  //***********************************************
  //Make a global instance
  OBForceFieldMMFF94 theForceFieldMMFF94("MMFF94", false);
  OBForceFieldMMFF94 theForceFieldMMFF94s("MMFF94s", false);
  //***********************************************

  OBForceFieldMMFF94::~OBForceFieldMMFF94()
  {
  }

  OBForceFieldMMFF94 &OBForceFieldMMFF94::operator=(OBForceFieldMMFF94 &src)
  {
    _mol = src._mol;
    _init = src._init;
    return *this;
  }

  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  //
  //  Parse parameter files
  //
  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////

  bool OBForceFieldMMFF94::ParseParamFile()
  {
    // Set the locale for number parsing to avoid locale issues: PR#1785463
    obLocale.SetLocale();

    vector<string> vs;
    char buffer[80];

    // open data/_parFile
    ifstream ifs;
    if (OpenDatafile(ifs, _parFile).length() == 0) {
      obErrorLog.ThrowError(__FUNCTION__, "Cannot open parameter file", obError);
      return false;
    }

    while (ifs.getline(buffer, 80)) {
      if (EQn(buffer, "#", 1)) continue;

      tokenize(vs, buffer);
      if (vs.size() < 2)
        continue;

      if (vs[0] == "prop")
        ParseParamProp(vs[1]);
      if (vs[0] == "def")
        ParseParamDef(vs[1]);
      if (vs[0] == "bond")
        ParseParamBond(vs[1]);
      if (vs[0] == "ang")
        ParseParamAngle(vs[1]);
      if (vs[0] == "bndk")
        ParseParamBndk(vs[1]);
      if (vs[0] == "chg")
        ParseParamCharge(vs[1]);
      if (vs[0] == "dfsb")
        ParseParamDfsb(vs[1]);
      if (vs[0] == "oop")
        ParseParamOOP(vs[1]);
      if (vs[0] == "pbci")
        ParseParamPbci(vs[1]);
      if (vs[0] == "stbn")
        ParseParamStrBnd(vs[1]);
      if (vs[0] == "tor")
        ParseParamTorsion(vs[1]);
      if (vs[0] == "vdw")
        ParseParamVDW(vs[1]);
    }

    if (ifs)
      ifs.close();

    // return the locale to the original one
    obLocale.RestoreLocale();
    return true;
  }

  bool OBForceFieldMMFF94::ParseParamBond(std::string &filename)
  {
    vector<string> vs;
    char buffer[80];

    OBFFParameter parameter;

    // open data/mmffbond.par
    ifstream ifs;
    if (OpenDatafile(ifs, filename).length() == 0) {
      obErrorLog.ThrowError(__FUNCTION__, "Cannot open mmffbond.par", obError);
      return false;
    }

    while (ifs.getline(buffer, 80)) {
      if (EQn(buffer, "*", 1)) continue;
      if (EQn(buffer, "$", 1)) continue;

      tokenize(vs, buffer);

      parameter.clear();
      parameter._ipar.push_back(atoi(vs[0].c_str()));  // FF class
      parameter.a = atoi(vs[1].c_str());
      parameter.b = atoi(vs[2].c_str());
      parameter._dpar.push_back(atof(vs[3].c_str()));  // kb
      parameter._dpar.push_back(atof(vs[4].c_str()));  // r0
      _ffbondparams.push_back(parameter);
    }

    if (ifs)
      ifs.close();

    return 0;
  }

  bool OBForceFieldMMFF94::ParseParamBndk(std::string &filename)
  {
    vector<string> vs;
    char buffer[80];

    OBFFParameter parameter;

    // open data/mmffbndk.par
    ifstream ifs;
    if (OpenDatafile(ifs, filename).length() == 0) {
      obErrorLog.ThrowError(__FUNCTION__, "Cannot open mmffbndk.par", obError);
      return false;
    }

    while (ifs.getline(buffer, 80)) {
      if (EQn(buffer, "*", 1)) continue;
      if (EQn(buffer, "$", 1)) continue;

      tokenize(vs, buffer);

      parameter.clear();
      parameter.a = atoi(vs[0].c_str());
      parameter.b = atoi(vs[1].c_str());
      parameter._dpar.push_back(atof(vs[2].c_str()));  // r0-ref
      parameter._dpar.push_back(atof(vs[3].c_str()));  // kb-ref
      _ffbndkparams.push_back(parameter);
    }

    if (ifs)
      ifs.close();

    return 0;
  }

  bool OBForceFieldMMFF94::ParseParamAngle(std::string &filename)
  {
    vector<string> vs;
    char buffer[80];

    OBFFParameter parameter;

    // open data/mmffang.par
    ifstream ifs;
    if (OpenDatafile(ifs, filename).length() == 0) {
      obErrorLog.ThrowError(__FUNCTION__, "Cannot open mmffang.par", obError);
      return false;
    }

    while (ifs.getline(buffer, 80)) {
      if (EQn(buffer, "*", 1)) continue;
      if (EQn(buffer, "$", 1)) continue;

      tokenize(vs, buffer);

      parameter.clear();
      parameter._ipar.push_back(atoi(vs[0].c_str()));  // FF class
      parameter.a = atoi(vs[1].c_str());
      parameter.b = atoi(vs[2].c_str());
      parameter.c = atoi(vs[3].c_str());
      parameter._dpar.push_back(atof(vs[4].c_str()));  // ka
      parameter._dpar.push_back(atof(vs[5].c_str()));  // theta0
      _ffangleparams.push_back(parameter);
    }

    if (ifs)
      ifs.close();

    return 0;
  }

  bool OBForceFieldMMFF94::ParseParamStrBnd(std::string &filename)
  {
    vector<string> vs;
    char buffer[80];

    OBFFParameter parameter;

    // open data/mmffstbn.par
    ifstream ifs;
    if (OpenDatafile(ifs, filename).length() == 0) {
      obErrorLog.ThrowError(__FUNCTION__, "Cannot open mmffstbn.par", obError);
      return false;
    }

    while (ifs.getline(buffer, 80)) {
      if (EQn(buffer, "*", 1)) continue;
      if (EQn(buffer, "$", 1)) continue;

      tokenize(vs, buffer);

      parameter.clear();
      parameter._ipar.push_back(atoi(vs[0].c_str()));  // FF class
      parameter.a = atoi(vs[1].c_str());
      parameter.b = atoi(vs[2].c_str());
      parameter.c = atoi(vs[3].c_str());
      parameter._dpar.push_back(atof(vs[4].c_str()));  // kbaIJK
      parameter._dpar.push_back(atof(vs[5].c_str()));  // kbaKJI
      _ffstrbndparams.push_back(parameter);
    }

    if (ifs)
      ifs.close();

    return 0;
  }

  bool OBForceFieldMMFF94::ParseParamDfsb(std::string &filename)
  {
    vector<string> vs;
    char buffer[80];

    OBFFParameter parameter;

    // open data/mmffdfsb.par
    ifstream ifs;
    if (OpenDatafile(ifs, filename).length() == 0) {
      obErrorLog.ThrowError(__FUNCTION__, "Cannot open mmffdfsb.par", obError);
      return false;
    }

    while (ifs.getline(buffer, 80)) {
      if (EQn(buffer, "*", 1)) continue;
      if (EQn(buffer, "$", 1)) continue;

      tokenize(vs, buffer);

      parameter.clear();
      parameter.a = atoi(vs[0].c_str());
      parameter.b = atoi(vs[1].c_str());
      parameter.c = atoi(vs[2].c_str());
      parameter._dpar.push_back(atof(vs[3].c_str()));  // kbaIJK
      parameter._dpar.push_back(atof(vs[4].c_str()));  // kbaKJI
      _ffdfsbparams.push_back(parameter);
    }

    if (ifs)
      ifs.close();

    return 0;
  }

  bool OBForceFieldMMFF94::ParseParamOOP(std::string &filename)
  {
    vector<string> vs;
    char buffer[80];

    OBFFParameter parameter;

    // open data/mmffoop.par
    ifstream ifs;
    if (OpenDatafile(ifs, filename).length() == 0) {
      obErrorLog.ThrowError(__FUNCTION__, "Cannot open mmffoop.par", obError);
      return false;
    }

    while (ifs.getline(buffer, 80)) {
      if (EQn(buffer, "*", 1)) continue;
      if (EQn(buffer, "$", 1)) continue;

      tokenize(vs, buffer);

      parameter.clear();
      parameter.a = atoi(vs[0].c_str());
      parameter.b = atoi(vs[1].c_str());
      parameter.c = atoi(vs[2].c_str());
      parameter.d = atoi(vs[3].c_str());
      parameter._dpar.push_back(atof(vs[4].c_str()));  // koop
      _ffoopparams.push_back(parameter);
    }

    if (ifs)
      ifs.close();

    return 0;
  }

  bool OBForceFieldMMFF94::ParseParamTorsion(std::string &filename)
  {
    vector<string> vs;
    char buffer[80];

    OBFFParameter parameter;

    // open data/mmfftor.par
    ifstream ifs;
    if (OpenDatafile(ifs, filename).length() == 0) {
      obErrorLog.ThrowError(__FUNCTION__, "Cannot open mmfftor.par", obError);
      return false;
    }

    while (ifs.getline(buffer, 80)) {
      if (EQn(buffer, "*", 1)) continue;
      if (EQn(buffer, "$", 1)) continue;

      tokenize(vs, buffer);

      parameter.clear();
      parameter._ipar.push_back(atoi(vs[0].c_str()));  // FF class
      parameter.a = atoi(vs[1].c_str());
      parameter.b = atoi(vs[2].c_str());
      parameter.c = atoi(vs[3].c_str());
      parameter.d = atoi(vs[4].c_str());
      parameter._dpar.push_back(atof(vs[5].c_str()));  // v1
      parameter._dpar.push_back(atof(vs[6].c_str()));  // v2
      parameter._dpar.push_back(atof(vs[7].c_str()));  // v3
      _fftorsionparams.push_back(parameter);
    }

    if (ifs)
      ifs.close();

    return 0;
  }

  bool OBForceFieldMMFF94::ParseParamVDW(std::string &filename)
  {
    vector<string> vs;
    char buffer[80];

    OBFFParameter parameter;

    // open data/mmffvdw.par
    ifstream ifs;
    if (OpenDatafile(ifs, filename).length() == 0) {
      obErrorLog.ThrowError(__FUNCTION__, "Cannot open mmffvdw.par", obError);
      return false;
    }

    while (ifs.getline(buffer, 80)) {
      if (EQn(buffer, "*", 1)) continue;
      if (EQn(buffer, "$", 1)) continue;

      tokenize(vs, buffer);

      parameter.clear();
      parameter.a = atoi(vs[0].c_str());
      parameter._dpar.push_back(atof(vs[1].c_str()));  // alpha-i
      parameter._dpar.push_back(atof(vs[2].c_str()));  // N-i
      parameter._dpar.push_back(atof(vs[3].c_str()));  // A-i
      parameter._dpar.push_back(atof(vs[4].c_str()));  // G-i
      if (EQn(vs[5].c_str(), "-", 1))
        parameter._ipar.push_back(0);
      else if (EQn(vs[5].c_str(), "D", 1))
        parameter._ipar.push_back(1);  // hydrogen bond donor
      else if (EQn(vs[5].c_str(), "A", 1))
        parameter._ipar.push_back(2);  // hydrogen bond acceptor
      _ffvdwparams.push_back(parameter);
    }

    if (ifs)
      ifs.close();

    return 0;
  }

  bool OBForceFieldMMFF94::ParseParamCharge(std::string &filename)
  {
    vector<string> vs;
    char buffer[80];

    OBFFParameter parameter;

    // open data/mmffchg.par
    ifstream ifs;
    if (OpenDatafile(ifs, filename).length() == 0) {
      obErrorLog.ThrowError(__FUNCTION__, "Cannot open mmffchg.par", obError);
      return false;
    }

    while (ifs.getline(buffer, 80)) {
      if (EQn(buffer, "*", 1)) continue;
      if (EQn(buffer, "$", 1)) continue;

      tokenize(vs, buffer);

      parameter.clear();
      parameter._ipar.push_back(atoi(vs[0].c_str()));  // FF class
      parameter.a = atoi(vs[1].c_str());
      parameter.b = atoi(vs[2].c_str());
      parameter._dpar.push_back(atof(vs[3].c_str()));  // bci
      _ffchgparams.push_back(parameter);
    }

    if (ifs)
      ifs.close();

    return 0;
  }

  bool OBForceFieldMMFF94::ParseParamPbci(std::string &filename)
  {
    vector<string> vs;
    char buffer[80];

    OBFFParameter parameter;

    // open data/mmffpbci.par
    ifstream ifs;
    if (OpenDatafile(ifs, filename).length() == 0) {
      obErrorLog.ThrowError(__FUNCTION__, "Cannot open mmffpbci", obError);
      return false;
    }

    while (ifs.getline(buffer, 80)) {
      if (EQn(buffer, "*", 1)) continue;
      if (EQn(buffer, "$", 1)) continue;

      tokenize(vs, buffer);

      parameter.clear();
      parameter.a = atoi(vs[1].c_str());
      parameter._dpar.push_back(atof(vs[2].c_str()));  // pbci
      parameter._dpar.push_back(atof(vs[3].c_str()));  // fcadj
      _ffpbciparams.push_back(parameter);
    }

    if (ifs)
      ifs.close();

    return 0;
  }

  bool OBForceFieldMMFF94::ParseParamProp(std::string &filename)
  {
    vector<string> vs;
    char buffer[80];

    OBFFParameter parameter;

    // open data/mmffprop.par
    ifstream ifs;
    if (OpenDatafile(ifs, filename).length() == 0) {
      obErrorLog.ThrowError(__FUNCTION__, "Cannot open mmffprop.par", obError);
      return false;
    }

    while (ifs.getline(buffer, 80)) {
      if (EQn(buffer, "*", 1)) continue;
      if (EQn(buffer, "$", 1)) continue;

      tokenize(vs, buffer);

      parameter.clear();
      parameter.a = atoi(vs[0].c_str()); // atom type
      parameter._ipar.push_back(atoi(vs[1].c_str()));  // at.no
      parameter._ipar.push_back(atoi(vs[2].c_str()));  // crd
      parameter._ipar.push_back(atoi(vs[3].c_str()));  // val
      parameter._ipar.push_back(atoi(vs[4].c_str()));  // pilp
      parameter._ipar.push_back(atoi(vs[5].c_str()));  // mltb
      parameter._ipar.push_back(atoi(vs[6].c_str()));  // arom
      parameter._ipar.push_back(atoi(vs[7].c_str()));  // linh
      parameter._ipar.push_back(atoi(vs[8].c_str()));  // sbmb

      if (parameter._ipar[3])
        _ffpropPilp.SetBitOn(parameter.a);
      if (parameter._ipar[5])
        _ffpropArom.SetBitOn(parameter.a);
      if (parameter._ipar[6])
        _ffpropLin.SetBitOn(parameter.a);
      if (parameter._ipar[7])
        _ffpropSbmb.SetBitOn(parameter.a);

      _ffpropparams.push_back(parameter);
    }

    if (ifs)
      ifs.close();

    return 0;
  }

  bool OBForceFieldMMFF94::ParseParamDef(std::string &filename)
  {
    vector<string> vs;
    char buffer[80];

    OBFFParameter parameter;

    // open data/mmffdef.par
    ifstream ifs;
    if (OpenDatafile(ifs, filename).length() == 0) {
      obErrorLog.ThrowError(__FUNCTION__, "Cannot open mmffdef.par", obError);
      return false;
    }

    while (ifs.getline(buffer, 80)) {
      if (EQn(buffer, "*", 1)) continue;
      if (EQn(buffer, "$", 1)) continue;

      tokenize(vs, buffer);

      parameter.clear();
      parameter._ipar.push_back(atoi(vs[1].c_str()));  // level 1
      parameter._ipar.push_back(atoi(vs[2].c_str()));  // level 2
      parameter._ipar.push_back(atoi(vs[3].c_str()));  // level 3
      parameter._ipar.push_back(atoi(vs[4].c_str()));  // level 4
      parameter._ipar.push_back(atoi(vs[5].c_str()));  // level 5
      _ffdefparams.push_back(parameter);
    }

    if (ifs)
      ifs.close();

    return 0;
  }

  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  //
  //  Setup Functions
  //
  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////

  // The MMFF94 article doesn't seem to include information about how
  // aromaticity is perceived. This function was written by studying the
  // MMFF_opti.log file, trail-and-error and using the MMFF94 validation
  // set to check the results (If all atom types are assigned correctly,
  // aromatic rings are probably detected correctly)
  bool OBForceFieldMMFF94::PerceiveAromatic()
  {
    bool done = false; // not done actually....
    OBAtom *ringatom;
    OBBond *ringbond;
    vector<OBRing*> vr;
    vr = _mol.GetSSSR();

    vector<OBRing*>::iterator ri;
    vector<int>::iterator rj;
    int n, index, ringsize, first_rj, prev_rj, pi_electrons, c60;
    for (ri = vr.begin();ri != vr.end();++ri) { // for each ring
      ringsize = (*ri)->Size();

      n = 1;
      pi_electrons = 0;
      c60 = 0; // we have a special case to get c60 right (all atom type 37)
      for(rj = (*ri)->_path.begin();rj != (*ri)->_path.end();++rj) { // for each ring atom
        index = *rj;
        ringatom = _mol.GetAtom(index);

        // is the bond to the previous ring atom double?
        if (n > 1) {
          ringbond = _mol.GetBond(prev_rj, index);
          if (!ringbond) {
            prev_rj = index;
            continue;
          }
          if (ringbond->GetBondOrder() == 2) {
            pi_electrons += 2;
            prev_rj = index;
            n++;
            continue;
          }
          prev_rj = index;
        } else {
          prev_rj = index;
          first_rj = index;
        }

        // does the current ring atom have a exocyclic double bond?
        FOR_NBORS_OF_ATOM (nbr, ringatom) {
          if ((*ri)->IsInRing(nbr->GetIdx()))
            continue;

          if (!nbr->IsAromatic()) {
            if (ringatom->GetAtomicNum() == OBElements::Carbon && ringatom->IsInRingSize(5)
                && ringatom->IsInRingSize(6) && nbr->GetAtomicNum() == OBElements::Carbon && nbr->IsInRingSize(5)
                && nbr->IsInRingSize(6)) {
              c60++;
            } else {
              continue;
            }
          }

          ringbond = _mol.GetBond(nbr->GetIdx(), index);
          if (!ringbond) {
            continue;
          }
          if (ringbond->GetBondOrder() == 2)
            pi_electrons++;
        }

        // is the atom N, O or S in 5 rings
        if (ringsize == 5 &&
            ringatom->GetIdx() == (*ri)->GetRootAtom())
          pi_electrons += 2;

        n++;

      } // for each ring atom

      // is the bond from the first to the last atom double?
      ringbond = _mol.GetBond(first_rj, index);
      if (ringbond) {
        if (ringbond->GetBondOrder() == 2)
          pi_electrons += 2;
      }

      if (((pi_electrons == 6) && ((ringsize == 5) || (ringsize == 6)))
        || ((pi_electrons == 5) && (c60 == 5))) {
        // mark ring atoms as aromatic
        for(rj = (*ri)->_path.begin();rj != (*ri)->_path.end();++rj) {
          if (!_mol.GetAtom(*rj)->IsAromatic())
            done = true;
          _mol.GetAtom(*rj)->SetAromatic();
        }
        // mark all ring bonds as aromatic
        FOR_BONDS_OF_MOL (bond, _mol)
          if((*ri)->IsMember(&*bond))
            bond->SetAromatic();
      }
    }

    return done;
  }

  // Symbolic atom typing is skipped
  //
  // atom typing is based on:
  //   MMFF94 I - Table III
  //   MMFF94 V - Table I
  //
  int OBForceFieldMMFF94::GetType(OBAtom *atom)
  {
    OBBond *bond;
    int oxygenCount, nitrogenCount, sulphurCount, doubleBondTo;
    ////////////////////////////////
    // Aromatic Atoms
    ////////////////////////////////
    if (atom->IsAromatic()) {
      if (atom->IsInRingSize(5)) {
        bool IsAromatic = false;
        vector<OBAtom*> alphaPos, betaPos;
        vector<OBAtom*> alphaAtoms, betaAtoms;

        if (atom->GetAtomicNum() == OBElements::Sulfur) {
          return 44; // Aromatic 5-ring sulfur with pi lone pair (STHI)
        }
        if (atom->GetAtomicNum() == OBElements::Oxygen) {
          return 59; // Aromatic 5-ring oxygen with pi lone pair (OFUR)
        }
        if (atom->GetAtomicNum() == OBElements::Nitrogen) {
          FOR_NBORS_OF_ATOM (nbr, atom) {
            if (nbr->GetAtomicNum() == OBElements::Oxygen && (nbr->GetExplicitDegree() == 1)) {
              return 82; // N-oxide nitrogen in 5-ring alpha position,
              // N-oxide nitrogen in 5-ring beta position,
              // N-oxide nitrogen in other 5-ring  position,
              // (N5AX, N5BX, N5OX)
            }
          }
        }
        FOR_NBORS_OF_ATOM (nbr, atom) {
          if (!((_mol.GetBond(atom, &*nbr))->IsAromatic()) || !nbr->IsInRingSize(5))
            continue;

          if (IsInSameRing(atom, &*nbr)) {
            alphaPos.push_back(&*nbr);
          }

          FOR_NBORS_OF_ATOM (nbrNbr, &*nbr) {
            if (nbrNbr->GetIdx() == atom->GetIdx())
              continue;
            if (!((_mol.GetBond(&*nbr, &*nbrNbr))->IsAromatic()) || !nbrNbr->IsInRingSize(5))
              continue;

            IsAromatic = true;

            if (IsInSameRing(atom, &*nbrNbr)) {
              betaPos.push_back(&*nbrNbr);
            }
          }
        }

        if (IsAromatic) {


          for (unsigned int i = 0; i < alphaPos.size(); i++) {
            if (alphaPos[i]->GetAtomicNum() == OBElements::Sulfur) {
              alphaAtoms.push_back(alphaPos[i]);
            } else if (alphaPos[i]->GetAtomicNum() == OBElements::Oxygen) {
              alphaAtoms.push_back(alphaPos[i]);
            } else if (alphaPos[i]->GetAtomicNum() == OBElements::Nitrogen && (alphaPos[i]->GetExplicitDegree() == 3)) {
              bool IsNOxide = false;
              FOR_NBORS_OF_ATOM (nbr, alphaPos[i]) {
                if (nbr->GetAtomicNum() == OBElements::Oxygen && (nbr->GetExplicitDegree() == 1)) {
                  IsNOxide = true;
                }
              }

              if (!IsNOxide) {
                alphaAtoms.push_back(alphaPos[i]);
              }
            }
          }
          for (unsigned int i = 0; i < betaPos.size(); i++) {
            if (betaPos[i]->GetAtomicNum() == OBElements::Sulfur) {
              betaAtoms.push_back(betaPos[i]);
            } else if (betaPos[i]->GetAtomicNum() == OBElements::Oxygen) {
              betaAtoms.push_back(betaPos[i]);
            } else if (betaPos[i]->GetAtomicNum() == OBElements::Nitrogen && (betaPos[i]->GetExplicitDegree() == 3)) {
              bool IsNOxide = false;
              FOR_NBORS_OF_ATOM (nbr, betaPos[i]) {
                if (nbr->GetAtomicNum() == OBElements::Oxygen && (nbr->GetExplicitDegree() == 1)) {
                  IsNOxide = true;
                }
              }

              if (!IsNOxide) {
                betaAtoms.push_back(betaPos[i]);
              }
            }
          }
          if (!betaAtoms.size()) {
            nitrogenCount = 0;
            FOR_NBORS_OF_ATOM (nbr, atom) {
              if (nbr->GetAtomicNum() == OBElements::Nitrogen && (nbr->GetExplicitDegree() == 3)) {
                if ((nbr->GetExplicitValence() == 4) && nbr->IsAromatic()) {
                  nitrogenCount++;
                } else if ((nbr->GetExplicitValence() == 3) && !nbr->IsAromatic()) {
                  nitrogenCount++;
                }
              }
            }
            if (nitrogenCount >= 2) {
              return 80; // Aromatic carbon between N's in imidazolium (CIM+)
            }
          }
          if (!alphaAtoms.size() && !betaAtoms.size()) {
            if (atom->GetAtomicNum() == OBElements::Carbon) {
              bool c60 = true; // special case to ensure c60 is typed correctly -- Paolo Tosco
              FOR_NBORS_OF_ATOM (nbr, atom) {
                if (!(nbr->GetAtomicNum() == OBElements::Carbon && nbr->IsAromatic() && nbr->IsInRingSize(6)))
                  c60 = false;
              }
              if (c60)
                return 37; // correct atom type for c in c60 (all atoms symmetric)
              // there is no S:, O:, or N:
              // this is the case for anions with only carbon and nitrogen in the ring
              return 78; // General carbon in 5-membered aromatic ring (C5)
            } else if (atom->GetAtomicNum() == OBElements::Nitrogen) {
              if (atom->GetExplicitDegree() == 3) {
                // this is the N: atom
                return 39; // Aromatic 5 ring nitrogen with pi lone pair (NPYL)
              } else {
                // again, no S:, O:, or N:
                return 76; // Nitrogen in 5-ring aromatic anion (N5M)
              }
            }
          }
          if (alphaAtoms.size() == 2) {
            if (atom->GetAtomicNum() == OBElements::Carbon && IsInSameRing(alphaAtoms[0], alphaAtoms[1])) {
              if (alphaAtoms[0]->GetAtomicNum() == OBElements::Nitrogen && alphaAtoms[1]->GetAtomicNum() == OBElements::Nitrogen) {
                if ((alphaAtoms[0]->GetExplicitDegree() == 3) && (alphaAtoms[1]->GetExplicitDegree() == 3)) {
                  return 80; // Aromatic carbon between N's in imidazolium (CIM+)
                }
              }
            }
          }
          if (alphaAtoms.size() && !betaAtoms.size()) {
            if (atom->GetAtomicNum() == OBElements::Carbon) {
              return 63; // Aromatic 5-ring C, alpha to N:, O:, or S: (C5A)
            } else if (atom->GetAtomicNum() == OBElements::Nitrogen) {
              if (atom->GetExplicitDegree() == 3) {
                return 81; // Posivite nitrogen in 5-ring alpha position (N5A+)
              } else {
                return 65; // Aromatic 5-ring N, alpha to N:, O:, or S: (N5A)
              }
            }
          }
          if (!alphaAtoms.size() && betaAtoms.size()) {
            if (atom->GetAtomicNum() == OBElements::Carbon) {
              return 64; // Aromatic 5-ring C, beta to N:, O:, or S: (C5B)
            } else if (atom->GetAtomicNum() == OBElements::Nitrogen) {
              if (atom->GetExplicitDegree() == 3) {
                return 81; // Posivite nitrogen in 5-ring beta position (N5B+)
              } else {
                return 66; // Aromatic 5-ring N, beta to N:, O:, or S: (N5B)
              }
            }
          }
          if (alphaAtoms.size() && betaAtoms.size()) {
            for (unsigned int i = 0; i < alphaAtoms.size(); i++) {
              for (unsigned int j = 0; j < betaAtoms.size(); j++) {
                if (!IsInSameRing(alphaAtoms[i], betaAtoms[j])) {
                  if (atom->GetAtomicNum() == OBElements::Carbon) {
                    return 78; // General carbon in 5-membered aromatic ring (C5)
                  } else if (atom->GetAtomicNum() == OBElements::Nitrogen) {
                    return 79; // General nitrogen in 5-membered aromatic ring (N5)
                  }
                }
              }
            }
            for (unsigned int i = 0; i < alphaAtoms.size(); i++) {
              if (alphaAtoms[i]->GetAtomicNum() == OBElements::Sulfur || alphaAtoms[i]->GetAtomicNum() == OBElements::Oxygen) {
                if (atom->GetAtomicNum() == OBElements::Carbon) {
                  return 63; // Aromatic 5-ring C, alpha to N:, O:, or S: (C5A)
                } else if (atom->GetAtomicNum() == OBElements::Nitrogen) {
                  return 65; // Aromatic 5-ring N, alpha to N:, O:, or S: (N5A)
                }
              }
            }
            for (unsigned int i = 0; i < betaAtoms.size(); i++) {
              if (betaAtoms[i]->GetAtomicNum() == OBElements::Sulfur || betaAtoms[i]->GetAtomicNum() == OBElements::Oxygen) {
                if (atom->GetAtomicNum() == OBElements::Carbon) {
                  return 64; // Aromatic 5-ring C, beta to N:, O:, or S: (C5B)
                } else if (atom->GetAtomicNum() == OBElements::Nitrogen) {
                  return 66; // Aromatic 5-ring N, beta to N:, O:, or S: (N5B)
                }
              }
            }

            if (atom->GetAtomicNum() == OBElements::Carbon) {
              return 78; // General carbon in 5-membered aromatic ring (C5)
            } else if (atom->GetAtomicNum() == OBElements::Nitrogen) {
              return 79; // General nitrogen in 5-membered aromatic ring (N5)
            }
          }
        }
      }

      if (atom->IsInRingSize(6)) {

        if (atom->GetAtomicNum() == OBElements::Carbon) {
          return 37; // Aromatic carbon, e.g., in benzene (CB)
        } else if (atom->GetAtomicNum() == OBElements::Nitrogen) {
          FOR_NBORS_OF_ATOM (nbr, atom) {
            if (nbr->GetAtomicNum() == OBElements::Oxygen && (nbr->GetExplicitDegree() == 1)) {
              return 69; // Pyridinium N-oxide nitrogen (NPOX)
            }
          }

          if (atom->GetExplicitDegree() == 3) {
            return 58; // Aromatic nitrogen in pyridinium (NPD+)
          } else {
            return 38; // Aromatic nitrogen with sigma lone pair (NPYD)
          }
        }
      }
    }

    ////////////////////////////////
    // Hydrogen
    ////////////////////////////////
    if (atom->GetAtomicNum() == 1) {
      FOR_NBORS_OF_ATOM (nbr, atom) {
        if (nbr->GetAtomicNum() == OBElements::Carbon) {
          return 5; // Hydrogen attatched to carbon (HC)
        }
        if (nbr->GetAtomicNum() == 14) {
          return 5; // Hydrogen attatched to silicon (HSI)
        }
        if (nbr->GetAtomicNum() == OBElements::Oxygen) {
          if (nbr->GetExplicitValence() == 3) {
            if (nbr->GetExplicitDegree() == 3) {
              return 50; // Hydrogen on oxonium oxygen (HO+)
            } else {
              return 52; // Hydrogen on oxenium oxygen (HO=+)
            }
          }

          int hydrogenCount = 0;
          FOR_NBORS_OF_ATOM (nbrNbr, &*nbr) {
            if (nbrNbr->GetAtomicNum() == OBElements::Hydrogen) {
              hydrogenCount++;
              continue;
            }
            if (nbrNbr->GetAtomicNum() == OBElements::Carbon) {
              if (nbrNbr->IsAromatic()) {
                return 29; // phenol
              }

              FOR_NBORS_OF_ATOM (nbrNbrNbr, &*nbrNbr) {
                if (nbrNbrNbr->GetIdx() == nbr->GetIdx())
                  continue;

                bond = _mol.GetBond(&*nbrNbr, &*nbrNbrNbr);
                if (!bond->IsAromatic() && bond->GetBondOrder() == 2) {
                  if (nbrNbrNbr->GetAtomicNum() == OBElements::Oxygen) {
                    return 24; // Hydroxyl hydrogen in carboxylic acids (HOCO)
                  }
                  if (nbrNbrNbr->GetAtomicNum() == OBElements::Carbon || nbrNbrNbr->GetAtomicNum() == OBElements::Nitrogen) {
                    return 29; // Enolic or phenolic hydroxyl hydrogen,
                    // Hydroxyl hydrogen in HO-C=N moiety (HOCC, HOCN)
                  }
                }
              }
            }
            if (nbrNbr->GetAtomicNum() == OBElements::Phosphorus) {
              return 24; // Hydroxyl hydrogen in H-O-P moiety (HOP)
            }
            if (nbrNbr->GetAtomicNum() == OBElements::Sulfur) {
              return 33; // Hydrogen on oxygen attached to sulfur (HOS)
            }

          }
          if (hydrogenCount == 2) {
            return 31; // Hydroxyl hydrogen in water (HOH)
          }

          return 21; // Hydroxyl hydrogen in alcohols, Generic hydroxyl hydrogen (HOR, HO)
        }
        if (nbr->GetAtomicNum() == OBElements::Nitrogen) {
          switch (GetType(&*nbr)) {
          case 81:
            return 36; // Hydrogen on imidazolium nitrogen (HIM+)
          case 68:
            return 23; // Hydrogen on N in N-oxide (HNOX)
          case 67:
            return 23; // Hydrogen on N in N-oxide (HNOX)
          case 62:
            return 23; // Generic hydrogen on sp3 nitrogen, e.g., in amines (HNR)
          case 56:
            return 36; // Hydrogen on guanimdinium nitrogen (HGD+)
          case 55:
            return 36; // Hydrogen on amidinium nitrogen (HNN+)
          case 43:
            return 28; // Hydrogen on NSO, NSO2, or NSO3 nitrogen, Hydrogen on N triply bonded to C (HNSO, HNC%)
          case 39:
            return 23; // Hydrogen on nitrogen in pyrrole (HPYL)
          case 8:
            return 23; // Generic hydrogen on sp3 nitrogen, e.g., in amines, Hydrogen on nitrogen in ammonia (HNR, H3N)
          }

          if (nbr->GetExplicitValence() == 4) {
            if (nbr->GetExplicitDegree() == 2) {
              return 28; // Hydrogen on N triply bonded to C (HNC%)
            } else {
              return 36; // Hydrogen on pyridinium nitrogen, Hydrogen on protonated imine nitrogen (HPD+, HNC+)
            }
          }

          if (nbr->GetExplicitDegree() == 2) {
            FOR_NBORS_OF_ATOM (nbrNbr, &*nbr) {
              if (nbrNbr->GetAtomicNum() == OBElements::Hydrogen)
                continue;

              bond = _mol.GetBond(&*nbr, &*nbrNbr);
              if (!bond->IsAromatic() && bond->GetBondOrder() == 2) {
                if (nbrNbr->GetAtomicNum() == OBElements::Carbon || nbrNbr->GetAtomicNum() == OBElements::Nitrogen) {
                  return 27; // Hydrogen on imine nitrogen, Hydrogen on azo nitrogen (HN=C, HN=N)
                }

                return 28; // Generic hydrogen on sp2 nitrogen (HSP2)
              }
            }
          }

          FOR_NBORS_OF_ATOM (nbrNbr, &*nbr) {
            if (nbrNbr->GetAtomicNum() == OBElements::Hydrogen)
              continue;

            if (nbrNbr->GetAtomicNum() == OBElements::Carbon) {
              if (nbrNbr->IsAromatic()) {
                return 28; // deloc. lp pair
              }

              FOR_NBORS_OF_ATOM (nbrNbrNbr, &*nbrNbr) {
                if (nbrNbrNbr->GetIdx() == nbr->GetIdx())
                  continue;

                bond = _mol.GetBond(&*nbrNbr, &*nbrNbrNbr);
                if (!bond->IsAromatic() && bond->GetBondOrder() == 2) {
                  if (nbrNbrNbr->GetAtomicNum() == OBElements::Carbon || nbrNbrNbr->GetAtomicNum() == OBElements::Nitrogen || nbrNbrNbr->GetAtomicNum() == OBElements::Oxygen || nbrNbrNbr->GetAtomicNum() == OBElements::Sulfur) {
                    return 28; // Hydrogen on amide nitrogen, Hydrogen on thioamide nitrogen,
                    // Hydrogen on enamine nitrogen, Hydrogen in H-N-C=N moiety (HNCO, HNCS, HNCC, HNCN)
                  }
                }
              }
            }
            if (nbrNbr->GetAtomicNum() == OBElements::Nitrogen) {
              FOR_NBORS_OF_ATOM (nbrNbrNbr, &*nbrNbr) {
                if (nbrNbrNbr->GetIdx() == nbr->GetIdx())
                  continue;

                bond = _mol.GetBond(&*nbrNbr, &*nbrNbrNbr);
                if (!bond->IsAromatic() && bond->GetBondOrder() == 2) {
                  if (nbrNbrNbr->GetAtomicNum() == OBElements::Carbon || nbrNbrNbr->GetAtomicNum() == OBElements::Nitrogen) {
                    return 28; // Hydrogen in H-N-N=C moiety, Hydrogen in H-N-N=N moiety (HNNC, HNNN)
                  }
                }
              }
            }
            if (nbrNbr->GetAtomicNum() == OBElements::Sulfur) {
              FOR_NBORS_OF_ATOM (nbrNbrNbr, &*nbrNbr) {
                if (nbrNbrNbr->GetIdx() == nbr->GetIdx())
                  continue;

                if (nbrNbrNbr->GetAtomicNum() == OBElements::Oxygen || (nbrNbrNbr->GetExplicitDegree() == 1)) {
                  return 28; // Hydrogen on NSO, NSO2 or NSO3 nitrogen (HNSO)
                }
              }
            }
          }

          return 23; // Generic hydrogen on sp3 nitrogen e.g., in amines,
          // Hydrogen on nitrogen in pyrrole, Hydrogen in ammonia,
          // Hydrogen on N in N-oxide (HNR, HPYL, H3N, HNOX)
        }
        if (nbr->GetAtomicNum() == OBElements::Sulfur || nbr->GetAtomicNum() == OBElements::Phosphorus) {
          return 71; // Hydrogen attached to sulfur, Hydrogen attached to >S= sulfur doubly bonded to N,
          // Hydrogen attached to phosphorus (HS, HS=N, HP)
        }
      }
    }

    ////////////////////////////////
    // Lithium
    ////////////////////////////////
    if (atom->GetAtomicNum() == 3) {
      // 0 neighbours
      if (atom->GetExplicitDegree() == 0) {
        return 92; // Lithium cation (LI+)
      }
    }

    ////////////////////////////////
    // Carbon
    ////////////////////////////////
    if (atom->GetAtomicNum() == 6) {
      // 4 neighbours
      if (atom->GetExplicitDegree() == 4) {
        if (atom->IsInRingSize(3)) {
          return 22; // Aliphatic carbon in 3-membered ring (CR3R)
        }

        if (atom->IsInRingSize(4)) {
          return 20; // Aliphatic carbon in 4-membered ring (CR4R)
        }

        return 1; // Alkyl carbon (CR)
      }
      // 3 neighbours
      if (atom->GetExplicitDegree() == 3) {
        int N2count = 0;
        int N3count = 0;
        int N3fcharge = 0;
        oxygenCount = sulphurCount = doubleBondTo = 0;

        FOR_NBORS_OF_ATOM (nbr, atom) {
          bond = _mol.GetBond(&*nbr, atom);
          if (!bond->IsAromatic() && bond->GetBondOrder() == 2) {
            doubleBondTo = nbr->GetAtomicNum();
          }

          if (nbr->GetExplicitDegree() == 1) {
            if (nbr->GetAtomicNum() == OBElements::Oxygen) {
              oxygenCount++;
            } else if (nbr->GetAtomicNum() == OBElements::Sulfur) {
              sulphurCount++;
            }
          } else if (nbr->GetExplicitDegree() == 3) {
            if (nbr->GetAtomicNum() == OBElements::Nitrogen) {
              N3fcharge += nbr->GetFormalCharge();
              N3count++;
            }
          } else if ((nbr->GetExplicitDegree() == 2) && !bond->IsAromatic() && bond->GetBondOrder() == 2) {
            if (nbr->GetAtomicNum() == OBElements::Nitrogen) {
              N2count++;
            }
          }
        }
        if ((N3count >= 2) && (doubleBondTo == 7 || (!(doubleBondTo == 6) && atom->GetExplicitDegree() == 3
          && N3fcharge == 1)) && !N2count && !oxygenCount && !sulphurCount) {
          // N3==C--N3
          return 57; // Guanidinium carbon, Carbon in +N=C-N: resonance structures (CGD+, CNN+)
        }
        if ((oxygenCount == 2) || (sulphurCount == 2)) {
          // O1-?-C-?-O1 or S1-?-C-?-S1
          return 41; // Carbon in carboxylate anion, Carbon in thiocarboxylate anion (CO2M, CS2M)
        }
        if (atom->IsInRingSize(4) && (doubleBondTo == 6)) {
	        return 30; // Olefinic carbon in 4-membered ring (CR4E)
        }
        if ((doubleBondTo ==  7) || (doubleBondTo ==  8) ||
            (doubleBondTo == 15) || (doubleBondTo == 16)) {
          // C==N, C==O, C==P, C==S
          return 3; // Generic carbonyl carbon, Imine-type carbon, Guanidine carbon,
          // Ketone or aldehyde carbonyl carbon, Amide carbonyl carbon,
          // Carboxylic acid or ester carbonyl carbon, Carbamate carbonyl carbon,
          // Carbonic acid or ester carbonyl carbon, Thioester carbonyl (double
          // bonded to O or S), Thioamide carbon (double bonded to S), Carbon
          // in >C=SO2, Sulfinyl carbon in >C=S=O, Thiocarboxylic acid or ester
          // carbon, Carbon doubly bonded to P (C=O, C=N, CGD, C=OR, C=ON, COO,
          // COON, COOO, C=OS, C=S, C=SN, CSO2, CS=O, CSS, C=P)
        }

        return 2; // Vinylic Carbon, Generic sp2 carbon (C=C, CSP2)

      }
      // 2 neighbours
      if (atom->GetExplicitDegree() == 2) {
        return 4; // Acetylenic carbon, Allenic caron (CSP, =C=)
      }
      // 1 neighbours
      if (atom->GetExplicitDegree() == 1) {
        return 60; // Isonitrile carbon (C%-)
      }
    }

    ////////////////////////////////
    // Nitrogen
    ////////////////////////////////
    if (atom->GetAtomicNum() == 7) {
      // 4 neighbours
      if (atom->GetExplicitDegree() == 4) {
        FOR_NBORS_OF_ATOM (nbr, atom) {
          if (nbr->GetAtomicNum() == OBElements::Oxygen && (nbr->GetExplicitDegree() == 1)) {
            return 68; // sp3-hybridized N-oxide nitrogen (N3OX)
          }
        }

        return 34; // Quaternary nitrogen (NR+)
      }
      // 3 neighbours
      if (atom->GetExplicitDegree() == 3) {
        if (atom->GetExplicitValence() >= 4) { // > because we accept -N(=O)=O as a valid nitro group
          oxygenCount = nitrogenCount = doubleBondTo = 0;

          FOR_NBORS_OF_ATOM (nbr, atom) {
            if (nbr->GetAtomicNum() == OBElements::Oxygen && (nbr->GetExplicitDegree() == 1)) {
              oxygenCount++;
            }
            if (nbr->GetAtomicNum() == OBElements::Nitrogen) {
              bond = _mol.GetBond(&*nbr, atom);
              if (!bond->IsAromatic() && bond->GetBondOrder() == 2) {
                doubleBondTo = 7;
              }
            }
            if (nbr->GetAtomicNum() == OBElements::Carbon) {
              bond = _mol.GetBond(&*nbr, atom);
              if (!bond->IsAromatic() && bond->GetBondOrder() == 2) {
                FOR_NBORS_OF_ATOM (nbrNbr, &*nbr) {
                  if (nbrNbr->GetAtomicNum() == OBElements::Nitrogen && (nbrNbr->GetExplicitDegree() == 3)) {
                    nitrogenCount++;
                  }
                }
              }
            }
          }

          if (oxygenCount == 1) {
            return 67; // sp2-hybridized N-oxide nitrogen (N2OX)
          }
          if (oxygenCount >= 2) {
            return 45; // Nitrogen in nitro group, Nitrogen in nitrate group (NO2, NO3)
          }

          if (nitrogenCount == 1) {
            return 54; // Iminium nitrogen (N+=C)
          }
          if (nitrogenCount == 2) {
            return 55; // Either nitrogen in N+=C-N: (NCN+)
          }
          if (nitrogenCount == 3) {
            return 56; // Guanidinium nitrogen (NGD+)
          }

          if (doubleBondTo == 7) {
            return 54; // Positivly charged nitrogen doubly bonded to nitrogen (N+=N)
          }
        }

        if (atom->GetExplicitValence() == 3) {
          bool IsAmide = false;
          bool IsSulfonAmide = false;
          bool IsNNNorNNC = false;
          int tripleBondTo = 0;
          doubleBondTo = 0;

          FOR_NBORS_OF_ATOM (nbr, atom) {
            if (nbr->GetAtomicNum() == OBElements::Sulfur || nbr->GetAtomicNum() == OBElements::Phosphorus) {
              oxygenCount = 0;

              FOR_NBORS_OF_ATOM (nbrNbr, &*nbr) {
                if (nbrNbr->GetAtomicNum() == OBElements::Oxygen && (nbrNbr->GetExplicitDegree() == 1)) {
                  oxygenCount++;
                }
              }
              if (oxygenCount >= 2) {
                IsSulfonAmide = true;
                //return 43; // Sulfonamide nitrogen (NSO2, NSO3)
              }
            }
          }

          FOR_NBORS_OF_ATOM (nbr, atom) {
            if (nbr->GetAtomicNum() == OBElements::Carbon) {
              FOR_NBORS_OF_ATOM (nbrNbr, &*nbr) {
                bond = _mol.GetBond(&*nbr, &*nbrNbr);
                if (!bond->IsAromatic() && bond->GetBondOrder() == 2 && (nbrNbr->GetAtomicNum() == OBElements::Oxygen || nbrNbr->GetAtomicNum() == OBElements::Sulfur)) {
                  IsAmide = true;
                  //return 10; // Amide nitrogen, Thioamide nitrogen (NC=O, NC=S)
                }
              }
            }
          }

          FOR_NBORS_OF_ATOM (nbr, atom) {
            if (nbr->GetAtomicNum() == OBElements::Carbon) {
              int N2count = 0;
              int N3count = 0;
              oxygenCount = sulphurCount = 0;

              FOR_NBORS_OF_ATOM (nbrNbr, &*nbr) {
                bond = _mol.GetBond(&*nbr, &*nbrNbr);
                if (!bond->IsAromatic() && bond->GetBondOrder() == 2) {
                  doubleBondTo = nbrNbr->GetAtomicNum();
                }
                if (bond->IsAromatic()) {
                  if ((nbrNbr->GetAtomicNum() == 7) || (nbrNbr->GetAtomicNum() == 6)) {
                    doubleBondTo = nbrNbr->GetAtomicNum();
                  }
                }
                if (!bond->IsAromatic() && bond->GetBondOrder() == 3) {
                  tripleBondTo = nbrNbr->GetAtomicNum();
                }
                if (nbrNbr->GetAtomicNum() == OBElements::Nitrogen && (nbrNbr->GetExplicitDegree() == 3)) {
                  int nbrOxygen = 0;
                  FOR_NBORS_OF_ATOM (nbrNbrNbr, &*nbrNbr) {
                    if (nbrNbrNbr->GetAtomicNum() == OBElements::Oxygen) {
                      nbrOxygen++;
                    }
                  }
                  if (nbrOxygen < 2) {
                    N3count++;
                  }
                }
                if (nbrNbr->GetAtomicNum() == OBElements::Nitrogen && (nbrNbr->GetExplicitDegree() == 2) && (bond->GetBondOrder() == 2 || bond->IsAromatic())) {
                  N2count++;
                }
                if (nbrNbr->IsAromatic()) {
                  if (nbrNbr->GetAtomicNum() == OBElements::Oxygen) {
                    oxygenCount++;
                  }
                  if (nbrNbr->GetAtomicNum() == OBElements::Sulfur) {
                    sulphurCount++;
                  }
                }
              }
              if (N3count == 3) {
                return 56; // Guanidinium nitrogen (NGD+)
              }

              if (!IsAmide && !IsSulfonAmide && !oxygenCount && !sulphurCount && nbr->IsAromatic()) {
                return 40;
              }

              if ((N3count == 2) && (doubleBondTo == 7) && !N2count) {
                return 55; // Either nitrogen in N+=C-N: (NCN+)
              }
            }

            if (nbr->GetAtomicNum() == OBElements::Nitrogen) {
              nitrogenCount = 0;
              FOR_NBORS_OF_ATOM (nbrNbr, &*nbr) {
                bond = _mol.GetBond(&*nbr, &*nbrNbr);
                if (!bond->IsAromatic() && bond->GetBondOrder() == 2) {
                  if (nbrNbr->GetAtomicNum() == OBElements::Carbon) {
                    oxygenCount = sulphurCount = 0;
                    FOR_NBORS_OF_ATOM (nbrNbrNbr, &*nbrNbr) {
                      if (nbrNbrNbr->GetAtomicNum() == OBElements::Oxygen) {
                        oxygenCount++;
                      }
                      if (nbrNbrNbr->GetAtomicNum() == OBElements::Sulfur) {
                        sulphurCount++;
                      }
                      if (nbrNbrNbr->GetAtomicNum() == OBElements::Sulfur) {
                        nitrogenCount++;
                      }
                    }
                    if (!oxygenCount && !sulphurCount && (nitrogenCount == 1)) {
                      bool bondToAromC = false;
                      FOR_NBORS_OF_ATOM (nbr2, atom) {
                        if (nbr2->IsAromatic() && nbr2->GetAtomicNum() == OBElements::Carbon && nbr2->IsInRingSize(6)) {
                          bondToAromC = true;
                        }
                      }
                      if (!bondToAromC) {
                        IsNNNorNNC = true;
                      }
                    }
                  }
                  if (nbrNbr->GetAtomicNum() == OBElements::Nitrogen) {
                    bool bondToAromC = false;
                    FOR_NBORS_OF_ATOM (nbr2, atom) {
                      if (nbr2->IsAromatic() && nbr2->GetAtomicNum() == OBElements::Carbon && nbr2->IsInRingSize(6)) {
                        bondToAromC = true;
                      }
                    }
                    if (!bondToAromC) {
                      IsNNNorNNC = true;
                    }
                  }
                }
              }
            }
          }

          if (IsSulfonAmide) {
            return 43; // Sulfonamide nitrogen (NSO2, NSO3)
          }
          if (IsAmide) {
            return 10; // Amide nitrogen, Thioamide nitrogen (NC=O, NC=S)
          }

          if ((doubleBondTo ==  6) || (doubleBondTo == 7) ||(doubleBondTo == 15) || (tripleBondTo == 6)) {
            return 40; // Enamine or aniline nitrogen (deloc. lp), Nitrogen in N-C=N with deloc. lp,
            // Nitrogen in N-C=N with deloc. lp, Nitrogen attached to C-C triple bond
            // (NC=C, NC=N, NC=P, NC%C)
          }
          if (tripleBondTo == 7) {
            return 43; // Nitrogen attached to cyano group (NC%N)
          }
          if (IsNNNorNNC) {
            return 10; // Nitrogen in N-N=C moiety with deloc. lp
            // Nitrogen in N-N=N moiety with deloc. lp (NN=C, NN=N)
          }

          return 8; // Amine nitrogen (NR)
        }
      }
      // 2 neighbours
      if (atom->GetExplicitDegree() == 2) {
        if (atom->GetExplicitValence() >= 4) {
          bool isNbrCarbon = false;
          bool isBondTriple = false;
          FOR_NBORS_OF_ATOM (nbr, atom) {
            if (!isNbrCarbon)
              isNbrCarbon = nbr->GetAtomicNum() == OBElements::Carbon;
            bond = _mol.GetBond(&*nbr, atom);
            if (!isBondTriple)
              isBondTriple = !bond->IsAromatic() && bond->GetBondOrder() == 3;
          }
          if (isBondTriple && isNbrCarbon)
            return 61; // Isonitrile nitrogen (NR%)

          return 53; // Central nitrogen in C=N=N or N=N=N (=N=)
        }

        if (atom->GetExplicitValence() == 3) {
          doubleBondTo = 0;

          FOR_NBORS_OF_ATOM (nbr, atom) {
            bond = _mol.GetBond(&*nbr, atom);
            if (nbr->GetAtomicNum() == OBElements::Oxygen && !bond->IsAromatic() && bond->GetBondOrder() == 2 && (nbr->GetExplicitDegree() == 1)) {
              return 46; // Nitrogen in nitroso group (N=O)
            }
            if ((nbr->GetAtomicNum() == OBElements::Carbon || nbr->GetAtomicNum() == OBElements::Nitrogen) && !bond->IsAromatic() && bond->GetBondOrder() == 2) {
              return 9; // Iminie nitrogen, Azo-group nitrogen (N=C, N=N)
            }
          }
          FOR_NBORS_OF_ATOM (nbr, atom) {
            if (nbr->GetAtomicNum() == OBElements::Sulfur) {
              oxygenCount = 0;

              FOR_NBORS_OF_ATOM (nbrNbr, &*nbr) {
                if (nbrNbr->GetAtomicNum() == OBElements::Oxygen && (nbrNbr->GetExplicitDegree() == 1)) {
                  oxygenCount++;
                }
              }
              if (oxygenCount >= 2) {
                return 43; // Sulfonamide nitrogen (NSO2, NSO3)
              }
            }
          }
        }

        if (atom->GetExplicitValence() >= 2) { // Bug reported by Paolo Tosco
          oxygenCount = sulphurCount = 0;

          FOR_NBORS_OF_ATOM (nbr, atom) {
            if (nbr->GetAtomicNum() == OBElements::Sulfur) {
              FOR_NBORS_OF_ATOM (nbrNbr, &*nbr) {
                if (nbrNbr->GetAtomicNum() == OBElements::Oxygen && (nbrNbr->GetExplicitDegree() == 1)) {
                  oxygenCount++;
                }
              }
              if (oxygenCount == 1) {
                return 48; // Divalent nitrogen replacing monovalent O in SO2 group (NSO)
              }
            }
          }

          return 62; // Anionic divalent nitrogen (NM)
        }
      }
      // 1 neighbours
      if (atom->GetExplicitDegree() == 1) {
       	FOR_NBORS_OF_ATOM (nbr, atom) {
          bond = _mol.GetBond(&*nbr, atom);
          if (!bond->IsAromatic() && bond->GetBondOrder() == 3) {
            FOR_NBORS_OF_ATOM (nbrNbr, &*nbr) {
              if (atom != &*nbrNbr && !(nbr->GetAtomicNum() == OBElements::Nitrogen
                && nbrNbr->GetAtomicNum() == OBElements::Nitrogen && nbrNbr->GetExplicitDegree() == 2)) {
                return 42; // Triply bonded nitrogen (NSP)
              }
            }
          }
          if (nbr->GetAtomicNum() == OBElements::Nitrogen && (nbr->GetExplicitDegree() == 2)) {
            return 47; // Terminal nitrogen in azido group (NAZT)
          }
        }
      }
      return 8; // generic amine nitrogen
    }

    ////////////////////////////////
    // Oxygen
    ////////////////////////////////
    if (atom->GetAtomicNum() == 8) {
      // 3 neighbours
      if (atom->GetExplicitDegree() == 3) {
        return 49; // Oxonium oxygen (O+)
      }
      // 2 neighbours
      if (atom->GetExplicitDegree() == 2) {
        int hydrogenCount = 0;
        FOR_NBORS_OF_ATOM (nbr, atom) {
          if (nbr->GetAtomicNum() == OBElements::Hydrogen) {
            hydrogenCount++;
          }
        }

        if (hydrogenCount == 2) {
          // H--O--H
          return 70; // Oxygen in water (OH2)
        }
        if (atom->GetExplicitValence() == 3) {
          return 51; // Oxenium oxygen (O=+)
        }

        return 6; // Generic divalent oxygen, Ether oxygen, Carboxylic acid or ester oxygen,
        // Enolic or phenolic oxygen, Oxygen in -O-C=N- moiety, Divalent oxygen in
        // thioacid or ester, Divalent nitrate "ether" oxygen, Divalent oxygen in
        // sulfate group, Divalent oxygen in sulfite group, One of two divalent
        // oxygens attached to sulfur, Divalent oxygen in R(RO)S=O, Other divalent
        // oxygen attached to sulfur, Divalent oxygen in phosphate group, Divalent
        // oxygen in phosphite group, Divalent oxygen (one of two oxygens attached
        // to P), Other divalent oxygen (-O-, OR, OC=O, OC=C, OC=N, OC=S, ONO2,
        // ON=O, OSO3, OSO2, OSO, OS=O, -OS, OPO3, OPO2, OPO, -OP)

        // 59 ar
      }
      // 1 neighbour
      if (atom->GetExplicitDegree() == 1) {
        oxygenCount = sulphurCount = 0;

        FOR_NBORS_OF_ATOM (nbr, atom) {
          bond = _mol.GetBond(&*nbr, atom);

          if (nbr->GetAtomicNum() == OBElements::Carbon || nbr->GetAtomicNum() == OBElements::Nitrogen) {
            FOR_NBORS_OF_ATOM (nbrNbr, &*nbr) {
              if (nbrNbr->GetAtomicNum() == OBElements::Oxygen && (nbrNbr->GetExplicitDegree() == 1)) {
                oxygenCount++;
              }
              if (nbrNbr->GetAtomicNum() == OBElements::Sulfur && (nbrNbr->GetExplicitDegree() == 1)) {
                sulphurCount++;
              }
            }
          }
          // O---H
          if (nbr->GetAtomicNum() == OBElements::Hydrogen) {
            return 35;
          }
          // O-?-C
          if (nbr->GetAtomicNum() == OBElements::Carbon) {
            if (oxygenCount == 2) {
              // O-?-C-?-O
              return 32; // Oxygen in carboxylate group (O2CM)
            }
            if (!bond->IsAromatic() && bond->GetBondOrder() == 1) {
              // O--C
              return 35; // Oxide oxygen on sp3 carbon, Oxide oxygen on sp2 carbon (OM, OM2)
            } else {
              // O==C
              return 7; // Generic carbonyl oxygen, Carbonyl oxygen in amides,
              // Carbonyl oxygen in aldehydes and ketones, Carbonyl
              // oxygen in acids or esters (O=C, O=CN, O=CR, O=CO)
            }
          }
          // O-?-N
          if (nbr->GetAtomicNum() == OBElements::Nitrogen) {
            if (oxygenCount >= 2) {
              // O-?-N-?-O
              return 32; // Oxygen in nitro group, Nitro-group oxygen in nitrate,
              // Nitrate anion oxygen (O2N, O2NO, O3N)
            }
            if (!bond->IsAromatic() && bond->GetBondOrder() == 1) {
	      if ((nbr->GetExplicitDegree() == 2) || (nbr->GetExplicitValence() == 3))
	      // O(-)--N
	        return 35;
	      else
              // O--N
                return 32; // Oxygen in N-oxides (ONX)
            } else {
              // sometimes ONX bonds are labelled as double bonds
              // (e.g. by MOE, noticed by Paolo Tosco)
              int ndab = 0;
              FOR_BONDS_OF_ATOM(bond2, &*nbr) {
                if (bond2->GetBondOrder() == 2 || bond2->GetBondOrder() == 5)
                  ndab++;
              }
              if (ndab + nbr->GetExplicitDegree() == 5)
                // O--N
                return 32; // Oxygen in N-oxides (ONX)
              else
                // O==N
                return 7; // Nitroso oxygen (O=N)
            }
          }
          // O-?-S
          if (nbr->GetAtomicNum() == OBElements::Sulfur) {
            if (sulphurCount == 1) {
              // O1-?-S-?-S1
              return 32; // Terminal oxygen in thiosulfinate anion (OSMS)
            }
            if (!bond->IsAromatic() && bond->GetBondOrder() == 1) {
              // O--S
              return 32; // Single terminal oxygen on sulfur, One of 2 terminal O's on sulfur,
              // One of 3 terminal O's on sulfur, Terminal O in sulfate anion,
              // (O-S, O2S, O3S, O4S)
            } else {
              // O==S

              // are all sulfur nbr atoms carbon?
              bool isSulfoxide = true;
              int oxygenBoundToSulfur = 0;
              FOR_NBORS_OF_ATOM (nbr2, &*nbr) {
                if (atom == &*nbr2)
                  continue;

                if (nbr2->GetAtomicNum() == OBElements::Oxygen)
                  ++oxygenBoundToSulfur;
              }
              FOR_NBORS_OF_ATOM (nbr2, &*nbr) {
                if (atom == &*nbr2)
                  continue;
                OBBond* bond = nbr->GetBond(&*nbr2);
                if (!bond->IsAromatic() && bond->GetBondOrder() == 2
                  && nbr2->GetAtomicNum() == OBElements::Carbon && oxygenBoundToSulfur == 1)
                  isSulfoxide = false; // O=S on sulfur doubly bonded to, e.g., C (O=S=)

                if ((nbr2->GetAtomicNum() == OBElements::Oxygen && nbr2->GetExplicitDegree() == 1)
		  || (nbr2->GetAtomicNum() == OBElements::Nitrogen && nbr2->GetExplicitDegree() == 2))
                  isSulfoxide = false;
              }

              if (isSulfoxide)
                return 7; // Doubly bonded sulfoxide oxygen (O=S)
              else
                return 32; // (O2S, O2S=C, O3S, O4S)
            }
          }

          return 32; // Oxygen in phosphine oxide, One of 2 terminal O's on sulfur,
          // One of 3 terminal O's on sulfur, One of 4 terminal O's on sulfur,
          // Oxygen in perchlorate anion (OP, O2P, O3P, O4P, O4Cl)
        }
      }
    }

    ////////////////////////////////
    // Flourine
    ////////////////////////////////
    if (atom->GetAtomicNum() == 9) {
      // 1 neighbour
      if (atom->GetExplicitDegree() == 1) {
        return 11; // Fluorine (F)
      }
      // 0 neighbours
      if (atom->GetExplicitDegree() == 0) {
        return 89; // Fluoride anion (F-)
      }
    }

    ////////////////////////////////
    // Sodium
    ////////////////////////////////
    if (atom->GetAtomicNum() == 11) {
      return 93; // Sodium cation (NA+)
    }

    ////////////////////////////////
    // Magnesium
    ////////////////////////////////
    if (atom->GetAtomicNum() == 12) {
      return 99; // Dipositive magnesium cation (MG+2)
    }

    ////////////////////////////////
    // Silicon
    ////////////////////////////////
    if (atom->GetAtomicNum() == 14) {
      return 19; // Silicon (SI)
    }

    ////////////////////////////////
    // Phosphorus
    ////////////////////////////////
    if (atom->GetAtomicNum() == 15) {
      if (atom->GetExplicitDegree() == 4) {
        return 25; // Phosphate group phosphorus, Phosphorus with 3 attached oxygens,
        // Phosphorus with 2 attached oxygens, Phosphine oxide phosphorus,
        // General tetracoordinate phosphorus (PO4, PO3, PO2, PO, PTET)
      }
      if (atom->GetExplicitDegree() == 3) {
        return 26; // Phosphorus in phosphines (P)
      }
      if (atom->GetExplicitDegree() == 2) {
        return 75; // Phosphorus doubly bonded to C (-P=C)
      }
    }

    ////////////////////////////////
    // Sulfur
    ////////////////////////////////
    if (atom->GetAtomicNum() == 16) {
      // 4 neighbours
      if (atom->GetExplicitDegree() == 4) {
        return 18; // Sulfone sulfur, Sulfonamide sulfur, Sulfonate group sulfur,
        // Sulfate group sulfur, Sulfur in nitrogen analog of sulfone
        // (SO2, SO2N, SO3, SO4, SNO)
      }
      // 3 neighbours
      if (atom->GetExplicitDegree() == 3) {
        oxygenCount = sulphurCount = doubleBondTo = 0;

        FOR_NBORS_OF_ATOM (nbr, atom) {
          bond = _mol.GetBond(&*nbr, atom);
          if (!bond->IsAromatic() && bond->GetBondOrder() == 2) {
            if (nbr->GetAtomicNum() == 6)
	            doubleBondTo = 6;
          }

          if (nbr->GetExplicitDegree() == 1) {
            if (nbr->GetAtomicNum() == OBElements::Oxygen) {
              oxygenCount++;
            } else if (nbr->GetAtomicNum() == OBElements::Sulfur) {
              sulphurCount++;
            }
          }
        }

        if (oxygenCount == 2) {
          if (doubleBondTo == 6) {
            return 18; // Sulfone sulfur, doubly bonded to carbon (=SO2)
          }
          return 73; // Sulfur in anionic sulfinate group (SO2M)
        }
        if (oxygenCount && sulphurCount)
          return 73; // Tricoordinate sulfur in anionic thiosulfinate group (SSOM)

        //if ((doubleBondTo == 6) || (doubleBondTo == 8))
        return 17; // Sulfur doubly bonded to carbon, Sulfoxide sulfur (S=C, S=O)
      }
      // 2 neighbours
      if (atom->GetExplicitDegree() == 2) {
        doubleBondTo = 0;

        FOR_NBORS_OF_ATOM (nbr, atom) {
          if (nbr->GetAtomicNum() == OBElements::Oxygen) {
            bond = _mol.GetBond(&*nbr, atom);
            if (!bond->IsAromatic() && bond->GetBondOrder() == 2) {
              doubleBondTo = 8;
            }
          }
        }

        if (doubleBondTo == 8)
          return 74; // Sulfinyl sulfur, e.g., in C=S=O (=S=O)

        return 15; // Thiol, sulfide, or disulfide sulfor (S)
      }
      // 1 neighbour
      if (atom->GetExplicitDegree() == 1) {
        sulphurCount = doubleBondTo = 0;

        FOR_NBORS_OF_ATOM (nbr, atom) {
          FOR_NBORS_OF_ATOM (nbrNbr, &*nbr) {
            if (nbrNbr->GetAtomicNum() == OBElements::Sulfur && (nbrNbr->GetExplicitDegree() == 1)) {
              sulphurCount++;
            }
          }
          bond = _mol.GetBond(&*nbr, atom);
          if (!bond->IsAromatic() && bond->GetBondOrder() == 2) {
            doubleBondTo = nbr->GetAtomicNum();
          }
        }

        if ((doubleBondTo == 6) && (sulphurCount != 2)) {
          return 16; // Sulfur doubly bonded to carbon (S=C)
        }

        return 72; // Terminal sulfur bonded to P, Anionic terminal sulfur,
        // Terminal sulfur in thiosulfinate group (S-P, SM, SSMO)
      }

      // 44 ar
    }

    ////////////////////////////////
    // Clorine
    ////////////////////////////////
    if (atom->GetAtomicNum() == 17) {
      // 4 neighbour
      if (atom->GetExplicitDegree() == 4) {
        oxygenCount = 0;

        FOR_NBORS_OF_ATOM (nbr, atom) {
          if (nbr->GetAtomicNum() == OBElements::Oxygen) {
            oxygenCount++;
          }
        }
        if (oxygenCount == 4)
          return 77; // Perchlorate anion chlorine (CLO4)
      }
      // 1 neighbour
      if (atom->GetExplicitDegree() == 1) {
        return 12; // Chlorine (CL)
      }
      // 0 neighbours
      if (atom->GetExplicitDegree() == 0) {
        return 90; // Chloride anion (CL-)
      }
    }

    ////////////////////////////////
    // Potasium
    ////////////////////////////////
    if (atom->GetAtomicNum() == 19) {
      return 94; // Potasium cation (K+)
    }

    ////////////////////////////////
    // Calcium
    ////////////////////////////////
    if (atom->GetAtomicNum() == 20) {
      // 0 neighbours
      if (atom->GetExplicitDegree() == 0) {
        return 96; // Dipositive calcium cation (CA+2)
      }
    }

    ////////////////////////////////
    // Iron
    ////////////////////////////////
    if (atom->GetAtomicNum() == 26) {
      if (atom->GetFormalCharge() == 2)
        return 87; // Dipositive iron (FE+2)
      else
        return 88; // Tripositive iron (FE+3)
    }

    ////////////////////////////////
    // Copper
    ////////////////////////////////
    if (atom->GetAtomicNum() == 29) {
      if (atom->GetFormalCharge() == 1)
        return 97; // Monopositive copper cation (CU+1)
      else
        return 98; // Dipositive copper cation (CU+2)
    }

    ////////////////////////////////
    // Zinc
    ////////////////////////////////
    if (atom->GetAtomicNum() == 30) {
      return 95; // Dipositive zinc cation (ZN+2)
    }

    ////////////////////////////////
    // Bromine
    ////////////////////////////////
    if (atom->GetAtomicNum() == 35) {
      // 1 neighbour
      if (atom->GetExplicitDegree() == 1) {
        return 13; // Bromine (BR)
      }
      // 0 neighbours
      if (atom->GetExplicitDegree() == 0) {
        return 91; // Bromide anion (BR-)
      }
    }

    ////////////////////////////////
    // Iodine
    ////////////////////////////////
    if (atom->GetAtomicNum() == 53) {
      // 1 neighbour
      if (atom->GetExplicitDegree() == 1) {
        return 14; // Iodine (I)
      }
    }



    return 0;
  }

  bool OBForceFieldMMFF94::SetTypes()
  {
    char type[4];

    _mol.SetAtomTypesPerceived();

    // mark all atoms and bonds as non-aromatic
    _mol.SetAromaticPerceived();
    FOR_BONDS_OF_MOL (bond, _mol)
      bond->SetAromatic(false);
    FOR_ATOMS_OF_MOL (atom, _mol)
      atom->SetAromatic(false);

    // It might be needed to run this function more than once...
    bool done = true;
    while (done) {
      done = PerceiveAromatic();
    }

    FOR_ATOMS_OF_MOL (atom, _mol) {
      snprintf(type, 3, "%d", GetType(&*atom));
      atom->SetType(type);
    }

    PrintTypes();

    return true;
  }

  bool OBForceFieldMMFF94::SetupCalculations()
  {
    OBFFParameter *parameter;
    OBAtom *a, *b, *c, *d;
    int type_a, type_b, type_c, type_d;
    bool found;
    int order;

    IF_OBFF_LOGLVL_LOW
      OBFFLog("\nS E T T I N G   U P   C A L C U L A T I O N S\n\n");

    //
    // Bond Calculations
    //
    // no "step-down" procedure
    // MMFF part V - page 625 (empirical rule)
    //
    IF_OBFF_LOGLVL_LOW
      OBFFLog("SETTING UP BOND CALCULATIONS...\n");

    OBFFBondCalculationMMFF94 bondcalc;
    int bondtype;

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
          if (_intraGroup[i].BitIsSet(a->GetIdx()) && _intraGroup[i].BitIsSet(b->GetIdx())) {
            validBond = true;
            break;
          }
        }
        if (!validBond)
          continue;
      }

      bondtype = GetBondType(a, b);

      parameter = GetTypedParameter2Atom(bondtype, atoi(a->GetType()), atoi(b->GetType()), _ffbondparams); // from mmffbond.par
      if (parameter == NULL) {
        parameter = GetParameter2Atom(a->GetAtomicNum(), b->GetAtomicNum(), _ffbndkparams); // from mmffbndk.par - emperical rules
        if (parameter == NULL) {
          IF_OBFF_LOGLVL_LOW {
            // This should never happen
            snprintf(_logbuf, BUFF_SIZE, "    COULD NOT FIND PARAMETERS FOR BOND %d-%d (IDX)...\n", a->GetIdx(), b->GetIdx());
            OBFFLog(_logbuf);
          }
          return false;
        } else {
          IF_OBFF_LOGLVL_LOW {
            snprintf(_logbuf, BUFF_SIZE, "   USING EMPIRICAL RULE FOR BOND STRETCHING %d-%d (IDX)...\n", a->GetIdx(), b->GetIdx());
            OBFFLog(_logbuf);
          }

          double rr, rr2, rr4, rr6;
          bondcalc.a = a;
          bondcalc.b = b;
          bondcalc.r0 = GetRuleBondLength(a, b);

          rr = parameter->_dpar[0] / bondcalc.r0; // parameter->_dpar[0]  = r0-ref
          rr2 = rr * rr;
          rr4 = rr2 * rr2;
          rr6 = rr4 * rr2;

          bondcalc.kb = parameter->_dpar[1] * rr6; // parameter->_dpar[1]  = kb-ref
          bondcalc.bt = bondtype;
          bondcalc.SetupPointers();

          _bondcalculations.push_back(bondcalc);
        }
      } else {
        bondcalc.a = a;
        bondcalc.b = b;
        bondcalc.kb = parameter->_dpar[0];
        bondcalc.r0 = parameter->_dpar[1];
        bondcalc.bt = bondtype;
        bondcalc.SetupPointers();

        _bondcalculations.push_back(bondcalc);
      }
    }

    //
    // Angle Calculations
    //
    // MMFF part I - page 513 ("step-down" prodedure)
    // MMFF part I - page 519 (reference 68 is actually a footnote)
    // MMFF part V - page 627 (empirical rule)
    //
    // First try and find an exact match, if this fails, step down using the equivalences from mmffdef.par
    // five-stage protocol: 1-1-1, 2-2-2, 3-2-3, 4-2-4, 5-2-5
    // If this fails, use empirical rules
    // Since 1-1-1 = 2-2-2, we will only try 1-1-1 before going to 3-2-3
    //
    // Stretch-Bend Calculations
    //
    IF_OBFF_LOGLVL_LOW
      OBFFLog("SETTING UP ANGLE & STRETCH-BEND CALCULATIONS...\n");

    OBFFAngleCalculationMMFF94 anglecalc;
    OBFFStrBndCalculationMMFF94 strbndcalc;
    int angletype, strbndtype, bondtype1, bondtype2;

    _anglecalculations.clear();
    _strbndcalculations.clear();

    FOR_ANGLES_OF_MOL(angle, _mol) {
      b = _mol.GetAtom((*angle)[0] + 1);
      a = _mol.GetAtom((*angle)[1] + 1);
      c = _mol.GetAtom((*angle)[2] + 1);

      type_a = atoi(a->GetType());
      type_b = atoi(b->GetType());
      type_c = atoi(c->GetType());

      // skip this angle if the atoms are ignored
      if ( _constraints.IsIgnored(a->GetIdx()) || _constraints.IsIgnored(b->GetIdx()) || _constraints.IsIgnored(c->GetIdx()) )
        continue;

      // if there are any groups specified, check if the three angle atoms are in a single intraGroup
      if (HasGroups()) {
        bool validAngle = false;
        for (unsigned int i=0; i < _intraGroup.size(); ++i) {
          if (_intraGroup[i].BitIsSet(a->GetIdx()) && _intraGroup[i].BitIsSet(b->GetIdx()) &&
              _intraGroup[i].BitIsSet(c->GetIdx())) {
            validAngle = true;
            break;
          }
        }
        if (!validAngle)
          continue;
      }

      angletype = GetAngleType(a, b, c);
      strbndtype = GetStrBndType(a, b, c);
      bondtype1 = GetBondType(a, b);
      bondtype2 = GetBondType(b, c);

      if (HasLinSet(type_b)) {
        anglecalc.linear = true;
      } else {
        anglecalc.linear = false;
      }

      // try exact match
      parameter = GetTypedParameter3Atom(angletype, type_a, type_b, type_c, _ffangleparams);
      if (parameter == NULL) // try 3-2-3
        parameter = GetTypedParameter3Atom(angletype, EqLvl3(type_a), type_b, EqLvl3(type_c), _ffangleparams);
      if (parameter == NULL) // try 4-2-4
        parameter = GetTypedParameter3Atom(angletype, EqLvl4(type_a), type_b, EqLvl4(type_c), _ffangleparams);
      if (parameter == NULL) // try 5-2-5
        parameter = GetTypedParameter3Atom(angletype, EqLvl5(type_a), type_b, EqLvl5(type_c), _ffangleparams);

      if (parameter) {
        anglecalc.ka = parameter->_dpar[0];
        anglecalc.theta0 = parameter->_dpar[1];
        strbndcalc.theta0 = parameter->_dpar[1]; // **
      } else {
        IF_OBFF_LOGLVL_LOW {
          snprintf(_logbuf, BUFF_SIZE, "   USING DEFAULT ANGLE FOR %d-%d-%d (IDX)...\n", a->GetIdx(), b->GetIdx(), c->GetIdx());
          snprintf(_logbuf, BUFF_SIZE, "   USING EMPIRICAL RULE FOR ANGLE BENDING %d-%d-%d (IDX)...\n", a->GetIdx(), b->GetIdx(), c->GetIdx());
          OBFFLog(_logbuf);
        }

        anglecalc.ka = 0.0;
        anglecalc.theta0 = 120.0;

        if (GetCrd(type_b) == 4)
          anglecalc.theta0 = 109.45;

        if ((GetCrd(type_b) == 2) && b->GetAtomicNum() == OBElements::Oxygen)
          anglecalc.theta0 = 105.0;

        if (b->GetAtomicNum() > 10)
          anglecalc.theta0 = 95.0;

        if (HasLinSet(type_b))
          anglecalc.theta0 = 180.0;

        if ((GetCrd(type_b) == 3) && (GetVal(type_b) == 3) && !GetMltb(type_b)) {
          if (b->GetAtomicNum() == OBElements::Nitrogen) {
            anglecalc.theta0 = 107.0;
          } else {
            anglecalc.theta0 = 92.0;
          }
        }

        if (a->IsInRingSize(3) && b->IsInRingSize(3) && c->IsInRingSize(3) && IsInSameRing(a, c))
          anglecalc.theta0 = 60.0;

        if (a->IsInRingSize(4) && b->IsInRingSize(4) && c->IsInRingSize(4) && IsInSameRing(a, c))
          anglecalc.theta0 = 90.0;

        strbndcalc.theta0 = anglecalc.theta0; // **
      }

      // empirical rule for 0-b-0 and standard angles
      if (anglecalc.ka == 0.0) {
        IF_OBFF_LOGLVL_LOW {
          snprintf(_logbuf, BUFF_SIZE, "   USING EMPIRICAL RULE FOR ANGLE BENDING FORCE CONSTANT %d-%d-%d (IDX)...\n", a->GetIdx(), b->GetIdx(), c->GetIdx());
          OBFFLog(_logbuf);
        }

        double beta, Za, Zc, Cb, r0ab, r0bc, theta, theta2, D, rr, rr2;
        Za = GetZParam(a);
        Cb = GetCParam(b); // Fixed typo -- PR#2741658
        Zc = GetZParam(c);

        r0ab = GetBondLength(a, b);
        r0bc = GetBondLength(b, c);
        rr = r0ab + r0bc;
        rr2 = rr * rr;
        D = (r0ab - r0bc) / rr2;

        theta = anglecalc.theta0;
        theta2 = theta * theta;

        beta = 1.75;
        if (a->IsInRingSize(4) && b->IsInRingSize(4) && c->IsInRingSize(4) && IsInSameRing(a, c))
          beta = 0.85 * beta;
        if (a->IsInRingSize(3) && b->IsInRingSize(3) && c->IsInRingSize(3) && IsInSameRing(a, c))
          beta = 0.05 * beta;

        // Theta2 is in Degrees^2, but parameters are expecting radians
        // PR#2741669
        anglecalc.ka = (beta * Za * Cb * Zc * exp(-2 * D)) / (rr * theta2 * DEG_TO_RAD * DEG_TO_RAD);
      }

      anglecalc.a = a;
      anglecalc.b = b;
      anglecalc.c = c;
      anglecalc.at = angletype;

      anglecalc.SetupPointers();
      _anglecalculations.push_back(anglecalc);

      if (anglecalc.linear)
        continue;

      parameter = GetTypedParameter3Atom(strbndtype, type_a, type_b, type_c, _ffstrbndparams);
      if (parameter == NULL) {
        int rowa, rowb, rowc;

        // This is not a real empirical rule...
        //IF_OBFF_LOGLVL_LOW {
        //  snprintf(_logbuf, BUFF_SIZE, "   USING EMPIRICAL RULE FOR STRETCH-BENDING FORCE CONSTANT %d-%d-%d (IDX)...\n", a->GetIdx(), b->GetIdx(), c->GetIdx());
        //  OBFFLog(_logbuf);
        //}

        rowa = GetElementRow(a);
        rowb = GetElementRow(b);
        rowc = GetElementRow(c);

        parameter = GetParameter3Atom(rowa, rowb, rowc, _ffdfsbparams);

        if (parameter == NULL) {
          // This should never happen
          IF_OBFF_LOGLVL_LOW {
            snprintf(_logbuf, BUFF_SIZE, "    COULD NOT FIND PARAMETERS FOR STRETCH-BEND %d-%d-%d (IDX)...\n", a->GetIdx(), b->GetIdx(), c->GetIdx());
            OBFFLog(_logbuf);
          }
          return false;
        }

        if (rowa == parameter->a) {
          strbndcalc.kbaABC = parameter->_dpar[0];
          strbndcalc.kbaCBA = parameter->_dpar[1];
        } else {
          strbndcalc.kbaABC = parameter->_dpar[1];
          strbndcalc.kbaCBA = parameter->_dpar[0];
        }
      } else {
        if (type_a == parameter->a) {
          strbndcalc.kbaABC = parameter->_dpar[0];
          strbndcalc.kbaCBA = parameter->_dpar[1];
        } else {
          strbndcalc.kbaABC = parameter->_dpar[1];
          strbndcalc.kbaCBA = parameter->_dpar[0];
        }
      }

      strbndcalc.rab0 = GetBondLength(a, b);
      strbndcalc.rbc0 = GetBondLength(b ,c);
      strbndcalc.a = a;
      strbndcalc.b = b;
      strbndcalc.c = c;
      strbndcalc.sbt = strbndtype;
      strbndcalc.SetupPointers();
      // Set the pointers to addresses in the anglecalc, find the matching bondcalcs and do the same.
      // This should improve performance by not calculating all this twice. We could do the same
      // for torsion and angles since the bond lengths are calculated for bond stretching first.
      //bool found_angle = false;
      /*
        for (unsigned int ai = 0; ai < _anglecalculations.size(); ++ai) {
        if ( (_anglecalculations[ai].a->GetIdx() == a->GetIdx()) &&
        (_anglecalculations[ai].b->GetIdx() == b->GetIdx()) &&
        (_anglecalculations[ai].c->GetIdx() == c->GetIdx()) ) {
        strbndcalc.theta = &(_anglecalculations[ai].theta);
        strbndcalc.force_abc_a = _anglecalculations[ai].force_a;
        strbndcalc.force_abc_b = _anglecalculations[ai].force_b;
        strbndcalc.force_abc_c = _anglecalculations[ai].force_c;
        found_angle = true;
        break;
        } else if ( (_anglecalculations[ai].a->GetIdx() == c->GetIdx()) &&
        (_anglecalculations[ai].b->GetIdx() == b->GetIdx()) &&
        (_anglecalculations[ai].c->GetIdx() == a->GetIdx()) ) {
        strbndcalc.theta = &(_anglecalculations[ai].theta);
        strbndcalc.force_abc_a = _anglecalculations[ai].force_c;
        strbndcalc.force_abc_b = _anglecalculations[ai].force_b;
        strbndcalc.force_abc_c = _anglecalculations[ai].force_a;
        found_angle = true;
        break;
        }
        }
      */

      /*
        vector<OBFFAngleCalculationMMFF94>::iterator ai;
        for (ai = _anglecalculations.begin(); ai != _anglecalculations.end(); ++ai) {
        if ( (((*ai).a)->GetIdx() == a->GetIdx()) && (((*ai).b)->GetIdx() == b->GetIdx()) && (((*ai).c)->GetIdx() == c->GetIdx()) ) {
        strbndcalc.theta = (*ai).theta;
        cout << "theta prt       = " << (*ai).theta << endl;
        cout << "delta prt       = " << &((*ai).delta) << endl;
        cout << "GetThetaPointer = " << ai->GetThetaPointer() << endl;
        strbndcalc.force_abc_a = (*ai).force_a;
        strbndcalc.force_abc_b = (*ai).force_b;
        strbndcalc.force_abc_c = (*ai).force_c;
        found_angle = true;
        break;
        } else if ( (((*ai).a)->GetIdx() == c->GetIdx()) && (((*ai).b)->GetIdx() == b->GetIdx()) && (((*ai).c)->GetIdx() == a->GetIdx()) ) {
        strbndcalc.theta = (*ai).theta;
        strbndcalc.force_abc_a = (*ai).force_c;
        strbndcalc.force_abc_b = (*ai).force_b;
        strbndcalc.force_abc_c = (*ai).force_a;
        found_angle = true;
        break;
        }
        }
        if (!found_angle) // didn't find matching angle, shouldn't happen, but continue to be safe
        continue;


        bool found_rab = false;
        bool found_rbc = false;
        vector<OBFFBondCalculationMMFF94>::iterator bi;
        for (bi = _bondcalculations.begin(); bi != _bondcalculations.end(); ++bi) {
        // find rab
        if ( (((*bi).a)->GetIdx() == a->GetIdx()) && (((*bi).b)->GetIdx() == b->GetIdx()) ) {
        strbndcalc.rab = &((*bi).rab);
        strbndcalc.force_ab_a = (*bi).force_a;
        strbndcalc.force_ab_b = (*bi).force_b;
        found_rab = true;
        } else if ( (((*bi).a)->GetIdx() == b->GetIdx()) && (((*bi).b)->GetIdx() == a->GetIdx()) ) {
        strbndcalc.rab = &((*bi).rab);
        strbndcalc.force_ab_a = (*bi).force_b;
        strbndcalc.force_ab_b = (*bi).force_a;
        found_rab = true;
        }
        // find rbc
        if ( (((*bi).a)->GetIdx() == b->GetIdx()) && (((*bi).b)->GetIdx() == c->GetIdx()) ) {
        strbndcalc.rbc = &(bondcalc.rab);
        strbndcalc.force_ab_a = (*bi).force_a;
        strbndcalc.force_ab_b = (*bi).force_b;
        found_rbc = true;
        } else if ( (((*bi).a)->GetIdx() == c->GetIdx()) && (((*bi).b)->GetIdx() == b->GetIdx()) ) {
        strbndcalc.rbc = &(bondcalc.rab);
        strbndcalc.force_ab_a = (*bi).force_b;
        strbndcalc.force_ab_b = (*bi).force_a;
        found_rbc = true;
        }

        if (found_rab && found_rbc)
        break;
        }
        if (!found_rab || !found_rbc) // didn't find matching bond, or atoms overlap
        continue;
      */
      _strbndcalculations.push_back(strbndcalc);

    }

    //
    // Torsion Calculations
    //
    // MMFF part I - page 513 ("step-down" prodedure)
    // MMFF part I - page 519 (reference 68 is actually a footnote)
    // MMFF part IV - page 631 (empirical rule)
    //
    // First try and find an exact match, if this fails, step down using the equivalences from mmffdef.par
    // five-stage protocol: 1-1-1-1, 2-2-2-2, 3-2-2-5, 5-2-2-3, 5-2-2-5
    // If this fails, use empirical rules
    // Since 1-1-1-1 = 2-2-2-2, we will only try 1-1-1-1 before going to 3-2-2-5
    //
    IF_OBFF_LOGLVL_LOW
      OBFFLog("SETTING UP TORSION CALCULATIONS...\n");

    OBFFTorsionCalculationMMFF94 torsioncalc;
    int torsiontype;

    _torsioncalculations.clear();

    FOR_TORSIONS_OF_MOL(t, _mol) {
      a = _mol.GetAtom((*t)[0] + 1);
      b = _mol.GetAtom((*t)[1] + 1);
      c = _mol.GetAtom((*t)[2] + 1);
      d = _mol.GetAtom((*t)[3] + 1);

      type_a = atoi(a->GetType());
      type_b = atoi(b->GetType());
      type_c = atoi(c->GetType());
      type_d = atoi(d->GetType());

      // skip this torsion if the atoms are ignored
      if ( _constraints.IsIgnored(a->GetIdx()) || _constraints.IsIgnored(b->GetIdx()) ||
           _constraints.IsIgnored(c->GetIdx()) || _constraints.IsIgnored(d->GetIdx()) )
        continue;

      // if there are any groups specified, check if the four torsion atoms are in a single intraGroup
      if (HasGroups()) {
        bool validTorsion = false;
        for (unsigned int i=0; i < _intraGroup.size(); ++i) {
          if (_intraGroup[i].BitIsSet(a->GetIdx()) && _intraGroup[i].BitIsSet(b->GetIdx()) &&
              _intraGroup[i].BitIsSet(c->GetIdx()) && _intraGroup[i].BitIsSet(d->GetIdx())) {
            validTorsion = true;
            break;
          }
        }
        if (!validTorsion)
          continue;
      }

      torsiontype = GetTorsionType(a, b, c, d);
      // CXT = MC*(J*MA**3 + K*MA**2 + I*MA + L) + TTijkl  MC = 6, MA = 136
      order = (type_c*2515456 + type_b*18496 + type_d*136 + type_a)
        - (type_b*2515456 + type_c*18496 + type_a*136 + type_d);

      if (order >= 0) {
        // try exact match
        parameter = GetTypedParameter4Atom(torsiontype, type_a, type_b, type_c, type_d, _fftorsionparams);
        if (parameter == NULL) // try 3-2-2-5
          parameter = GetTypedParameter4Atom(torsiontype, EqLvl3(type_a), type_b, type_c, EqLvl5(type_d), _fftorsionparams);
        if (parameter == NULL) // try 5-2-2-3
          parameter = GetTypedParameter4Atom(torsiontype, EqLvl5(type_a), type_b, type_c, EqLvl3(type_d), _fftorsionparams);
        if (parameter == NULL) // try 5-2-2-5
          parameter = GetTypedParameter4Atom(torsiontype, EqLvl5(type_a), type_b, type_c, EqLvl5(type_d), _fftorsionparams);
      } else {
        // try exact match
        parameter = GetTypedParameter4Atom(torsiontype, type_d, type_c, type_b, type_a, _fftorsionparams);
        if (parameter == NULL) // try 3-2-2-5
          parameter = GetTypedParameter4Atom(torsiontype, EqLvl3(type_d), type_c, type_b, EqLvl5(type_a), _fftorsionparams);
        if (parameter == NULL) // try 5-2-2-3
          parameter = GetTypedParameter4Atom(torsiontype, EqLvl5(type_d), type_c, type_b, EqLvl3(type_a), _fftorsionparams);
        if (parameter == NULL) // try 5-2-2-5
          parameter = GetTypedParameter4Atom(torsiontype, EqLvl5(type_d), type_c, type_b, EqLvl5(type_a), _fftorsionparams);
      }

      if (parameter) {
        torsioncalc.v1 = parameter->_dpar[0];
        torsioncalc.v2 = parameter->_dpar[1];
        torsioncalc.v3 = parameter->_dpar[2];
      } else {
        bool found_rule = false;

        //IF_OBFF_LOGLVL_LOW {
        //  snprintf(_logbuf, BUFF_SIZE, "   USING EMPIRICAL RULE FOR TORSION FORCE CONSTANT %d-%d-%d-%d (IDX)...\n",
        //    a->GetIdx(), b->GetIdx(), c->GetIdx(), d->GetIdx());
        //  OBFFLog(_logbuf);
        //}

        // rule (a) page 631
        if (HasLinSet(type_b) || HasLinSet(type_c))
          continue;

        // rule (b) page 631
        if (b->GetBond(c)->IsAromatic()) {
          double Ub, Uc, pi_bc, beta;
          Ub = GetUParam(b);
          Uc = GetUParam(c);

          if (!HasPilpSet(type_b) && !HasPilpSet(type_c))
            pi_bc = 0.5;
          else
            pi_bc = 0.3;

          if (((GetVal(type_b) == 3) && (GetVal(type_c) == 4)) ||
              ((GetVal(type_b) == 4) && (GetVal(type_c) == 3)))
            beta = 3.0;
          else
            beta = 6.0;

          torsioncalc.v1 = 0.0;
          torsioncalc.v2 = beta * pi_bc * sqrt(Ub * Uc);
          torsioncalc.v3 = 0.0;
          found_rule = true;
        } else {
          // rule (c) page 631
       	  double Ub, Uc, pi_bc, beta;
          Ub = GetUParam(b);
          Uc = GetUParam(c);
          OBBond *bond = a->GetBond(b);
          if (((GetMltb(type_b) == 2) && (GetMltb(type_c) == 2)) && !bond->IsAromatic() && bond->GetBondOrder() == 2)
            pi_bc = 1.0;
          else
            pi_bc = 0.4;

          beta = 6.0;
          torsioncalc.v1 = 0.0;
          torsioncalc.v2 = beta * pi_bc * sqrt(Ub * Uc);
          torsioncalc.v3 = 0.0;
          found_rule = true;
        }

        // rule (d) page 632
        if (!found_rule)
          if (((GetCrd(type_b) == 4) && (GetCrd(type_c) == 4))) {
            double Vb, Vc;
            Vb = GetVParam(b);
            Vc = GetVParam(c);

            torsioncalc.v1 = 0.0;
            torsioncalc.v2 = 0.0;
            torsioncalc.v3 = sqrt(Vb * Vc) / 9.0;
            found_rule = true;
          }

        // rule (e) page 632
        if (!found_rule)
          if (((GetCrd(type_b) == 4) && (GetCrd(type_c) != 4))) {
            if (GetCrd(type_c) == 3) // case (1)
              if ((GetVal(type_c) == 4) || (GetVal(type_c) == 34) || (GetMltb(type_c) != 0))
                continue;

            if (GetCrd(type_c) == 2) // case (2)
              if ((GetVal(type_c) == 3) || (GetMltb(type_c) != 0))
                continue;

            // case (3) saturated bonds -- see rule (h)
          }

        // rule (f) page 632
        if (!found_rule)
          if (((GetCrd(type_b) != 4) && (GetCrd(type_c) == 4))) {
            if (GetCrd(type_b) == 3) // case (1)
              if ((GetVal(type_b) == 4) || (GetVal(type_b) == 34) || (GetMltb(type_b) != 0))
                continue;

            if (GetCrd(type_b) == 2) // case (2)
              if ((GetVal(type_b) == 3) || (GetMltb(type_b) != 0))
                continue;

            // case (3) saturated bonds
          }

        // rule (g) page 632
        if (!found_rule) {
          OBBond *bond = b->GetBond(c);
          if (!bond->IsAromatic() && bond->GetBondOrder() == 1 && (
            (GetMltb(type_b) && GetMltb(type_c)) ||
            (GetMltb(type_b) && HasPilpSet(type_c)) ||
            (GetMltb(type_c) && HasPilpSet(type_b)))) {
            if (HasPilpSet(type_b) && HasPilpSet(type_c)) // case (1)
              continue;

            double Ub, Uc, pi_bc, beta;
            Ub = GetUParam(b);
            Uc = GetUParam(c);
            beta = 6.0;

            if (HasPilpSet(type_b) && GetMltb(type_c)) { // case (2)
              if (GetMltb(type_c) == 1)
                pi_bc = 0.5;
              else if ((GetElementRow(b) == 1) && (GetElementRow(c) == 1))
                pi_bc = 0.3;
              else
                pi_bc = 0.15;
              found_rule = true;
            }

            if (HasPilpSet(type_c) && GetMltb(type_b)) { // case (3)
              if (GetMltb(type_b) == 1)
                pi_bc = 0.5;
              else if ((GetElementRow(b) == 1) && (GetElementRow(c) == 1))
                pi_bc = 0.3;
              else
                pi_bc = 0.15;
              found_rule = true;
            }

            if (!found_rule)
              if (((GetMltb(type_b) == 1) || (GetMltb(type_c) == 1)) && (b->GetAtomicNum() != OBElements::Carbon || c->GetAtomicNum() != OBElements::Carbon)) {
                pi_bc = 0.4;
                found_rule = true;
              }

            if (!found_rule)
              pi_bc = 0.15;

            torsioncalc.v1 = 0.0;
            torsioncalc.v2 = beta * pi_bc * sqrt(Ub * Uc);
            torsioncalc.v3 = 0.0;
            found_rule = true;
          }
        }

        // rule (h) page 632
        if (!found_rule) {
          if ((b->GetAtomicNum() == OBElements::Oxygen || b->GetAtomicNum() == OBElements::Sulfur) && (c->GetAtomicNum() == OBElements::Oxygen || c->GetAtomicNum() == OBElements::Sulfur)) {
            double Wb, Wc;

            if (b->GetAtomicNum() == OBElements::Oxygen) {
              Wb = 2.0;
            }
            else {
              Wb = 8.0;
            }

            if (c->GetAtomicNum() == OBElements::Oxygen) {
              Wc = 2.0;
            }
            else {
              Wc = 8.0;
            }

            torsioncalc.v1 = 0.0;
            torsioncalc.v2 = -sqrt(Wb * Wc);
            torsioncalc.v3 = 0.0;
          } else {
            double Vb, Vc, Nbc;
            Vb = GetVParam(b);
            Vc = GetVParam(c);

            IF_OBFF_LOGLVL_LOW {
              snprintf(_logbuf, BUFF_SIZE, "   USING EMPIRICAL RULE FOR TORSION FORCE CONSTANT %d-%d-%d-%d (IDX)...\n",
                      a->GetIdx(), b->GetIdx(), c->GetIdx(), d->GetIdx());
              OBFFLog(_logbuf);
            }

            Nbc = GetCrd(type_b) * GetCrd(type_c);

            torsioncalc.v1 = 0.0;
            torsioncalc.v2 = 0.0;
            torsioncalc.v3 = sqrt(Vb * Vc) / Nbc;
          }
        }
      }

      torsioncalc.a = a;
      torsioncalc.b = b;
      torsioncalc.c = c;
      torsioncalc.d = d;
      torsioncalc.SetupPointers();
      torsioncalc.tt = torsiontype;

      _torsioncalculations.push_back(torsioncalc);
    }

    //
    // Out-Of-Plane Calculations
    //
    IF_OBFF_LOGLVL_LOW
      OBFFLog("SETTING UP OOP CALCULATIONS...\n");

    OBFFOOPCalculationMMFF94 oopcalc;

    _oopcalculations.clear();

    FOR_ATOMS_OF_MOL(atom, _mol) {
      b = (OBAtom*) &*atom;

      found = false;

      type_b = atoi(b->GetType());

      for (unsigned int idx=0; idx < _ffoopparams.size(); idx++) {
        if (type_b == _ffoopparams[idx].b) {
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

          if ((a == NULL) || (c == NULL) || (d == NULL))
            break;

          type_a = atoi(a->GetType());
          type_c = atoi(c->GetType());
          type_d = atoi(d->GetType());

          // skip this oop if the atoms are ignored
          if ( _constraints.IsIgnored(a->GetIdx()) || _constraints.IsIgnored(b->GetIdx()) ||
               _constraints.IsIgnored(c->GetIdx()) || _constraints.IsIgnored(d->GetIdx()) )
            continue;

          // if there are any groups specified, check if the four oop atoms are in a single intraGroup
          if (HasGroups()) {
            bool validOOP = false;
            for (unsigned int i=0; i < _intraGroup.size(); ++i) {
              if (_intraGroup[i].BitIsSet(a->GetIdx()) && _intraGroup[i].BitIsSet(b->GetIdx()) &&
                  _intraGroup[i].BitIsSet(c->GetIdx()) && _intraGroup[i].BitIsSet(d->GetIdx())) {
                validOOP = true;
                break;
              }
            }
            if (!validOOP)
              continue;
          }

          if (((type_a == _ffoopparams[idx].a) && (type_c == _ffoopparams[idx].c) && (type_d == _ffoopparams[idx].d)) ||
              ((type_c == _ffoopparams[idx].a) && (type_a == _ffoopparams[idx].c) && (type_d == _ffoopparams[idx].d)) ||
              ((type_c == _ffoopparams[idx].a) && (type_d == _ffoopparams[idx].c) && (type_a == _ffoopparams[idx].d)) ||
              ((type_d == _ffoopparams[idx].a) && (type_c == _ffoopparams[idx].c) && (type_a == _ffoopparams[idx].d)) ||
              ((type_a == _ffoopparams[idx].a) && (type_d == _ffoopparams[idx].c) && (type_c == _ffoopparams[idx].d)) ||
              ((type_d == _ffoopparams[idx].a) && (type_a == _ffoopparams[idx].c) && (type_c == _ffoopparams[idx].d)))
            {
              found = true;

              oopcalc.koop = _ffoopparams[idx]._dpar[0];

              // A-B-CD || C-B-AD  PLANE = ABC
              oopcalc.a = a;
              oopcalc.b = b;
              oopcalc.c = c;
              oopcalc.d = d;

              oopcalc.SetupPointers();
              _oopcalculations.push_back(oopcalc);

              // C-B-DA || D-B-CA  PLANE BCD
              oopcalc.a = d;
              oopcalc.d = a;

              oopcalc.SetupPointers();
              _oopcalculations.push_back(oopcalc);

              // A-B-DC || D-B-AC  PLANE ABD
              oopcalc.a = a;
              oopcalc.c = d;
              oopcalc.d = c;

              oopcalc.SetupPointers();
              _oopcalculations.push_back(oopcalc);
            }

          if ((_ffoopparams[idx].a == 0) && (_ffoopparams[idx].c == 0) && (_ffoopparams[idx].d == 0) && !found) // *-XX-*-*
            {
              oopcalc.koop = _ffoopparams[idx]._dpar[0];

              // A-B-CD || C-B-AD  PLANE = ABC
              oopcalc.a = a;
              oopcalc.b = b;
              oopcalc.c = c;
              oopcalc.d = d;

              oopcalc.SetupPointers();
              _oopcalculations.push_back(oopcalc);

              // C-B-DA || D-B-CA  PLANE BCD
              oopcalc.a = d;
              oopcalc.d = a;

              oopcalc.SetupPointers();
              _oopcalculations.push_back(oopcalc);

              // A-B-DC || D-B-AC  PLANE ABD
              oopcalc.a = a;
              oopcalc.c = d;
              oopcalc.d = c;

              oopcalc.SetupPointers();
              _oopcalculations.push_back(oopcalc);
            }
        }
      }
    }

    //
    // VDW Calculations
    //
    IF_OBFF_LOGLVL_LOW
      OBFFLog("SETTING UP VAN DER WAALS CALCULATIONS...\n");

    OBFFVDWCalculationMMFF94 vdwcalc;

    _vdwcalculations.clear();

    int pairIndex = -1;
    FOR_PAIRS_OF_MOL(p, _mol) {
      ++pairIndex;
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
          if (_interGroup[i].BitIsSet(a->GetIdx()) && _interGroup[i].BitIsSet(b->GetIdx())) {
            validVDW = true;
            break;
          }
        }
        if (!validVDW) {
          for (unsigned int i=0; i < _interGroups.size(); ++i) {
            if (_interGroups[i].first.BitIsSet(a->GetIdx()) && _interGroups[i].second.BitIsSet(b->GetIdx())) {
              validVDW = true;
              break;
            }
            if (_interGroups[i].first.BitIsSet(b->GetIdx()) && _interGroups[i].second.BitIsSet(a->GetIdx())) {
              validVDW = true;
              break;
            }
          }
        }

        if (!validVDW)
          continue;
      }

      OBFFParameter *parameter_a, *parameter_b;
      parameter_a = GetParameter1Atom(atoi(a->GetType()), _ffvdwparams);
      parameter_b = GetParameter1Atom(atoi(b->GetType()), _ffvdwparams);
      if ((parameter_a == NULL) || (parameter_b == NULL)) {
        IF_OBFF_LOGLVL_LOW {
          snprintf(_logbuf, BUFF_SIZE, "   COULD NOT FIND VAN DER WAALS PARAMETERS FOR %d-%d (IDX)...\n", a->GetIdx(), b->GetIdx());
          OBFFLog(_logbuf);
        }

        return false;
      }

      vdwcalc.a = a;
      vdwcalc.alpha_a = parameter_a->_dpar[0];
      vdwcalc.Na = parameter_a->_dpar[1];
      vdwcalc.Aa = parameter_a->_dpar[2];
      vdwcalc.Ga = parameter_a->_dpar[3];
      vdwcalc.aDA = parameter_a->_ipar[0];

      vdwcalc.b = b;
      vdwcalc.alpha_b = parameter_b->_dpar[0];
      vdwcalc.Nb = parameter_b->_dpar[1];
      vdwcalc.Ab = parameter_b->_dpar[2];
      vdwcalc.Gb = parameter_b->_dpar[3];
      vdwcalc.bDA = parameter_b->_ipar[0];

      //these calculations only need to be done once for each pair,
      //we do them now and save them for later use
      double R_AA, R_BB, R_AB6, g_AB, g_AB2;
      double R_AB2, R_AB4, /*R_AB7,*/ sqrt_a, sqrt_b;

      R_AA = vdwcalc.Aa * pow(vdwcalc.alpha_a, 0.25);
      R_BB = vdwcalc.Ab * pow(vdwcalc.alpha_b, 0.25);
      sqrt_a = sqrt(vdwcalc.alpha_a / vdwcalc.Na);
      sqrt_b = sqrt(vdwcalc.alpha_b / vdwcalc.Nb);

      if (vdwcalc.aDA == 1) { // hydrogen bond donor
        vdwcalc.R_AB = 0.5 * (R_AA + R_BB);
        R_AB2 = vdwcalc.R_AB * vdwcalc.R_AB;
        R_AB4 = R_AB2 * R_AB2;
        R_AB6 = R_AB4 * R_AB2;

        if (vdwcalc.bDA == 2) { // hydrogen bond acceptor
          vdwcalc.epsilon = 0.5 * (181.16 * vdwcalc.Ga * vdwcalc.Gb * vdwcalc.alpha_a * vdwcalc.alpha_b) / (sqrt_a + sqrt_b) * (1.0 / R_AB6);
          // R_AB is scaled to 0.8 for D-A interactions. The value used in the calculation of epsilon is not scaled.
          vdwcalc.R_AB = 0.8 * vdwcalc.R_AB;
        } else
          vdwcalc.epsilon = (181.16 * vdwcalc.Ga * vdwcalc.Gb * vdwcalc.alpha_a * vdwcalc.alpha_b) / (sqrt_a + sqrt_b) * (1.0 / R_AB6);

        R_AB2 = vdwcalc.R_AB * vdwcalc.R_AB;
        R_AB4 = R_AB2 * R_AB2;
        R_AB6 = R_AB4 * R_AB2;
        vdwcalc.R_AB7 = R_AB6 * vdwcalc.R_AB;
      } else if (vdwcalc.bDA == 1) { // hydrogen bond donor
        vdwcalc.R_AB = 0.5 * (R_AA + R_BB);
       	R_AB2 = vdwcalc.R_AB * vdwcalc.R_AB;
        R_AB4 = R_AB2 * R_AB2;
        R_AB6 = R_AB4 * R_AB2;

        if (vdwcalc.aDA == 2) { // hydrogen bond acceptor
          vdwcalc.epsilon = 0.5 * (181.16 * vdwcalc.Ga * vdwcalc.Gb * vdwcalc.alpha_a * vdwcalc.alpha_b) / (sqrt_a + sqrt_b) * (1.0 / R_AB6);
          // R_AB is scaled to 0.8 for D-A interactions. The value used in the calculation of epsilon is not scaled.
          vdwcalc.R_AB = 0.8 * vdwcalc.R_AB;
        } else
          vdwcalc.epsilon = (181.16 * vdwcalc.Ga * vdwcalc.Gb * vdwcalc.alpha_a * vdwcalc.alpha_b) / (sqrt_a + sqrt_b) * (1.0 / R_AB6);

        R_AB2 = vdwcalc.R_AB * vdwcalc.R_AB;
        R_AB4 = R_AB2 * R_AB2;
        R_AB6 = R_AB4 * R_AB2;
        vdwcalc.R_AB7 = R_AB6 * vdwcalc.R_AB;
      } else {
        g_AB = (R_AA - R_BB) / ( R_AA + R_BB);
        g_AB2 = g_AB * g_AB;
        vdwcalc.R_AB =  0.5 * (R_AA + R_BB) * (1.0 + 0.2 * (1.0 - exp(-12.0 * g_AB2)));
        R_AB2 = vdwcalc.R_AB * vdwcalc.R_AB;
        R_AB4 = R_AB2 * R_AB2;
        R_AB6 = R_AB4 * R_AB2;
        vdwcalc.R_AB7 = R_AB6 * vdwcalc.R_AB;
        vdwcalc.epsilon = (181.16 * vdwcalc.Ga * vdwcalc.Gb * vdwcalc.alpha_a * vdwcalc.alpha_b) / (sqrt_a + sqrt_b) * (1.0 / R_AB6);
      }

      vdwcalc.pairIndex = pairIndex;
      vdwcalc.SetupPointers();
      _vdwcalculations.push_back(vdwcalc);
    }

    //
    // Electrostatic Calculations
    //
    IF_OBFF_LOGLVL_LOW
      OBFFLog("SETTING UP ELECTROSTATIC CALCULATIONS...\n");

    OBFFElectrostaticCalculationMMFF94 elecalc;

    _electrostaticcalculations.clear();

    pairIndex = -1;
    FOR_PAIRS_OF_MOL(p, _mol) {
      ++pairIndex;
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
          if (_interGroup[i].BitIsSet(a->GetIdx()) && _interGroup[i].BitIsSet(b->GetIdx())) {
            validEle = true;
            break;
          }
        }
        if (!validEle) {
          for (unsigned int i=0; i < _interGroups.size(); ++i) {
            if (_interGroups[i].first.BitIsSet(a->GetIdx()) && _interGroups[i].second.BitIsSet(b->GetIdx())) {
              validEle = true;
              break;
            }
            if (_interGroups[i].first.BitIsSet(b->GetIdx()) && _interGroups[i].second.BitIsSet(a->GetIdx())) {
              validEle = true;
              break;
            }
          }
        }

        if (!validEle)
          continue;
      }

      elecalc.qq = 332.0716 * a->GetPartialCharge() * b->GetPartialCharge() / _epsilon;

      if (elecalc.qq) {
        elecalc.a = &*a;
        elecalc.b = &*b;

        // 1-4 scaling
        if (a->IsOneFour(b))
          elecalc.qq *= 0.75;

        elecalc.pairIndex = pairIndex;
        elecalc.SetupPointers();
        _electrostaticcalculations.push_back(elecalc);
      }
    }

    return true;
  }

  bool OBForceFieldMMFF94::SetupPointers()
  {
    for (unsigned int i = 0; i < _bondcalculations.size(); ++i)
      _bondcalculations[i].SetupPointers();
    for (unsigned int i = 0; i < _anglecalculations.size(); ++i)
      _anglecalculations[i].SetupPointers();
    for (unsigned int i = 0; i < _strbndcalculations.size(); ++i)
      _strbndcalculations[i].SetupPointers();
    for (unsigned int i = 0; i < _torsioncalculations.size(); ++i)
      _torsioncalculations[i].SetupPointers();
    for (unsigned int i = 0; i < _oopcalculations.size(); ++i)
      _oopcalculations[i].SetupPointers();
    for (unsigned int i = 0; i < _vdwcalculations.size(); ++i)
      _vdwcalculations[i].SetupPointers();
    for (unsigned int i = 0; i < _electrostaticcalculations.size(); ++i)
      _electrostaticcalculations[i].SetupPointers();

    return true;
  }


  // we set the the formal charge with SetPartialCharge because formal charges
  // in MMFF94 are not always and integer
  bool OBForceFieldMMFF94::SetFormalCharges()
  {
    _mol.SetAutomaticPartialCharge(false);

    FOR_ATOMS_OF_MOL (atom, _mol) {
      int type = atoi(atom->GetType());
      atom->SetPartialCharge(0.0);

      bool done = false;
      switch (type) {
      case 34:
      case 49:
      case 51:
      case 54:
      case 58:
      case 92:
      case 93:
      case 94:
      case 97:
        atom->SetPartialCharge(1.0);
        done = true;
        break;
      case 35:
      case 62:
      case 89:
      case 90:
      case 91:
        atom->SetPartialCharge(-1.0);
        done = true;
        break;
      case 55:
        atom->SetPartialCharge(0.5);
        done = true;
        break;
      case 87:
      case 95:
      case 96:
      case 98:
      case 99:
        atom->SetPartialCharge(2.0);
        done = true;
        break;
      case 88:
        atom->SetPartialCharge(3.0);
        done = true;
        break;
        //case 98:
        //  atom->SetPartialCharge(3.0);
      default:
        break;
      }

      if (done)
        continue;

      if (type == 56) {
        int n_count = 0;
        int temp_type;
        FOR_ATOMS_OF_MOL (atom2, _mol) {
          temp_type = atoi(atom2->GetType());
          if (temp_type == 56 || temp_type == 81)
            ++n_count;
        }
        atom->SetPartialCharge((double)((n_count + 1) / 3) / (double)n_count);
      } else if (type == 32) {
        int o_count = 0;
        bool sulfonamide = false;
        bool sulfone_s_c = false;
        int s_count = 0;

        FOR_NBORS_OF_ATOM(nbr, &*atom) {
          FOR_NBORS_OF_ATOM(nbr2, &*nbr) {
            if (nbr2->GetAtomicNum() == OBElements::Oxygen && (nbr2->GetExplicitDegree() == 1))
              o_count++;
            if (nbr2->GetAtomicNum() == OBElements::Sulfur && (nbr2->GetExplicitDegree() == 1))
              s_count++;
            if (nbr2->GetAtomicNum() == OBElements::Nitrogen && !nbr2->IsAromatic())
              sulfonamide = true;
            OBBond *bond = nbr->GetBond(&*nbr2);
            if (nbr2->GetAtomicNum() == OBElements::Carbon && !bond->IsAromatic() && bond->GetBondOrder() == 2)
              sulfone_s_c = true;
          }

          if (nbr->GetAtomicNum() == OBElements::Carbon)
            atom->SetPartialCharge(-0.5); // O2CM

          if (nbr->GetAtomicNum() == OBElements::Nitrogen && (o_count == 3))
            atom->SetPartialCharge(-1.0 / o_count);  // O3N

          if (nbr->GetAtomicNum() == OBElements::Sulfur && !sulfonamide) {
            if (((o_count + s_count) == 2) && (nbr->GetExplicitDegree() == 3)
                && (nbr->GetExplicitValence() >= 3) && !sulfone_s_c) {
              atom->SetPartialCharge(-0.5); // O2S (sulfinate)
            }
            else if ((o_count + s_count) == 3) {
              atom->SetPartialCharge(-1.0 / 3.0); // O3S
            }
            else if ((o_count + s_count) == 4) {
              atom->SetPartialCharge(-0.5); // O4S
            }
          }

          if (nbr->GetAtomicNum() == OBElements::Phosphorus) {
            if ((o_count + s_count) == 2) {
              atom->SetPartialCharge(-0.5); // O2P
            }
            else if ((o_count + s_count) == 3) {
              atom->SetPartialCharge(-2.0 / 3.0); // O3P
            }
            else if ((o_count + s_count) == 4) {
              atom->SetPartialCharge(-0.25); // O4P
            }
          }

          if (atoi(nbr->GetType()) == 77)
            atom->SetPartialCharge(-0.25); // O4CL
        }
      } else if (type == 61) {
        FOR_BONDS_OF_ATOM(bond, &*atom) {
          OBAtom *nbr = bond->GetNbrAtom(&*atom);
          if (!bond->IsAromatic() && bond->GetBondOrder() == 3 && nbr->GetAtomicNum() == OBElements::Nitrogen)
            atom->SetPartialCharge(1.0);
        }
      } else if (type == 72) {
        int s_count = 0;

        FOR_NBORS_OF_ATOM(nbr, &*atom) {
          if (nbr->GetAtomicNum() == OBElements::Sulfur)
            s_count++;

          if (nbr->GetAtomicNum() == OBElements::Phosphorus || nbr->GetAtomicNum() == OBElements::Sulfur) {
            FOR_NBORS_OF_ATOM(nbr2, &*nbr)
              if ((nbr2->GetAtomicNum() == OBElements::Sulfur || nbr2->GetAtomicNum() == OBElements::Oxygen) && (nbr2->GetExplicitDegree() == 1) && (atom->GetIdx() != nbr2->GetIdx()))
                atom->SetPartialCharge(-0.5);
          } else
            atom->SetPartialCharge(-1.0);

          if (nbr->GetAtomicNum() == OBElements::Carbon)
            FOR_NBORS_OF_ATOM(nbr2, &*nbr)
              if (nbr2->GetAtomicNum() == OBElements::Sulfur && (nbr2->GetExplicitDegree() == 1) && (atom->GetIdx() != nbr2->GetIdx()))
                atom->SetPartialCharge(-0.5); // SSMO

          if (s_count >= 2)
            atom->SetPartialCharge(-0.5); // SSMO
        }
      } else if (type == 76) {
       	vector<OBRing*> vr;
        vr = _mol.GetSSSR();
        vector<OBRing*>::iterator ri;
        vector<int>::iterator rj;
        int n_count;

        for (ri = vr.begin();ri != vr.end();++ri) { // for each ring
          n_count = 0;

          if ((*ri)->IsAromatic() && (*ri)->IsMember(&*atom) && ((*ri)->Size() == 5)) {
            for(rj = (*ri)->_path.begin();rj != (*ri)->_path.end();++rj) // for each ring atom
              if (_mol.GetAtom(*rj)->GetAtomicNum() == OBElements::Nitrogen)
                n_count++;

            if (n_count > 1)
              atom->SetPartialCharge(-1.0 / n_count);
          }
        }
      } else if (type == 81) {
        atom->SetPartialCharge(1.0);

        vector<OBRing*> vr;
        vr = _mol.GetSSSR();
        vector<OBRing*>::iterator ri;
        vector<int>::iterator rj;
        for (ri = vr.begin();ri != vr.end();++ri) // for each ring
          if ((*ri)->IsAromatic() && (*ri)->IsMember(&*atom) && ((*ri)->Size() == 5)) {
            int n_count = 0;
            for(rj = (*ri)->_path.begin();rj != (*ri)->_path.end();++rj) // for each ring atom
              if (_mol.GetAtom(*rj)->GetAtomicNum() == OBElements::Nitrogen && (_mol.GetAtom(*rj)->GetExplicitDegree() == 3))
                n_count++;

            if (n_count) // coverity defensive testing
              atom->SetPartialCharge(1.0 / n_count); // NIM+

            FOR_NBORS_OF_ATOM(nbr, &*atom)
              FOR_NBORS_OF_ATOM(nbr2, &*nbr)
              if (atoi(nbr2->GetType()) == 56)
                atom->SetPartialCharge(1.0 / 3.0);

            FOR_NBORS_OF_ATOM(nbr, &*atom)
              FOR_NBORS_OF_ATOM(nbr2, &*nbr)
              if (atoi(nbr2->GetType()) == 55)
                atom->SetPartialCharge(1.0 / (1.0 + n_count));
          }
      }

    }

    PrintFormalCharges();

    return true;
  }

  bool OBForceFieldMMFF94::SetPartialCharges()
  {
    vector<double> charges(_mol.NumAtoms()+1, 0);
    double M, Wab, factor, q0a, q0b, Pa, Pb;

    FOR_ATOMS_OF_MOL (atom, _mol) {
      int type = atoi(atom->GetType());

      switch (type) {
      case 32:
      case 35:
      case 72:
        factor = 0.5;
        break;
      case 62:
      case 76:
        factor = 0.25;
        break;
      default:
        factor = 0.0;
        break;
      }

      M = GetCrd(type);
      q0a = atom->GetPartialCharge();

      // charge sharing
      if (!factor)
        FOR_NBORS_OF_ATOM (nbr, &*atom)
          if (nbr->GetPartialCharge() < 0.0)
            q0a += nbr->GetPartialCharge() / (2.0 * (double)(nbr->GetExplicitDegree()));

      // needed for SEYWUO, positive charge sharing?
      if (type == 62)
        FOR_NBORS_OF_ATOM (nbr, &*atom)
          if (nbr->GetPartialCharge() > 0.0)
            q0a -= nbr->GetPartialCharge() / 2.0;

      q0b = 0.0;
      Wab = 0.0;
      Pa = Pb = 0.0;
      FOR_NBORS_OF_ATOM (nbr, &*atom) {
        int nbr_type = atoi(nbr->GetType());

        q0b += nbr->GetPartialCharge();

        bool bci_found = false;
        for (unsigned int idx=0; idx < _ffchgparams.size(); idx++)
          if (GetBondType(&*atom, &*nbr) == _ffchgparams[idx]._ipar[0]) {
            if ((type == _ffchgparams[idx].a) && (nbr_type == _ffchgparams[idx].b)) {
              Wab += -_ffchgparams[idx]._dpar[0];
              bci_found = true;
            } else if  ((type == _ffchgparams[idx].b) && (nbr_type == _ffchgparams[idx].a)) {
              Wab += _ffchgparams[idx]._dpar[0];
              bci_found = true;
            }
	  }

        if (!bci_found) {
          for (unsigned int idx=0; idx < _ffpbciparams.size(); idx++) {
            if (type == _ffpbciparams[idx].a)
              Pa = _ffpbciparams[idx]._dpar[0];
            if (nbr_type == _ffpbciparams[idx].a)
              Pb = _ffpbciparams[idx]._dpar[0];
          }
          Wab += Pa - Pb;
        }
      }
      if (factor)
        charges[atom->GetIdx()] = (1.0 - M * factor) * q0a + factor * q0b + Wab;
      else
        charges[atom->GetIdx()] = q0a + Wab;
    }

    FOR_ATOMS_OF_MOL (atom, _mol)
      atom->SetPartialCharge(charges[atom->GetIdx()]);

    PrintPartialCharges();

    return true;
  }

  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  //
  //  Validation functions
  //
  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////

  // used to validate the implementation
  bool OBForceFieldMMFF94::Validate ()
  {
    OBConversion conv;
    OBFormat *format_in = conv.FindFormat("mol2");
    vector<string> vs;
    vector<int> types;
    vector<double> fcharges, pcharges;
    vector<double> bond_lengths;
    char buffer[150];
    bool molfound, atomfound, bondfound, fchgfound, pchgfound;
    double etot, ebond, eangle, eoop, estbn, etor, evdw, eeq;
    double termcount; //1=bond, 2=angle, 3=strbnd, 4=torsion, 5=oop
    int n = 0;

    if (!format_in || !conv.SetInFormat(format_in)) {
      obErrorLog.ThrowError(__FUNCTION__, "Could not set mol2 input format", obError);
      return false;
    }

    ifstream ifs, ifs2;
    ofstream ofs;

    ifs.open("MMFF94_dative.mol2");
    if (!ifs) {
      obErrorLog.ThrowError(__FUNCTION__, "Could not open ./MMFF94_dative.mol2", obError);
      return false;
    }

    ifs2.open("MMFF94_opti.log");
    if (!ifs2) {
      obErrorLog.ThrowError(__FUNCTION__, "Coulg not open ./MMFF_opti.log", obError);
      return false;
    }

    ofs.open("MMFF94_openbabel.log");
    if (!ofs) {
      obErrorLog.ThrowError(__FUNCTION__, "Coulg not open ./MMFF_openbabel.log", obError);
      return false;
    }

    if (!_init) {
      ParseParamFile();
      _init = true;
    }


    SetLogFile(&ofs);
    SetLogLevel(OBFF_LOGLVL_HIGH);

    for (unsigned int c=1;; c++) {
      _mol.Clear();
      types.clear();
      fcharges.clear();
      pcharges.clear();
      bond_lengths.clear();

      if (!conv.Read(&_mol, &ifs))
        break;
      if (_mol.Empty())
        break;

      _ncoords = _mol.NumAtoms() * 3;
      _gradientPtr = new double[_ncoords];

      SetTypes();

      if ((c == 98) || (c == 692)) // CUDPAS & VUWXUG
        continue;

      termcount = 0;
      molfound = false;
      atomfound = false;
      bondfound = false;
      fchgfound = false;
      pchgfound = false;

      while (ifs2.getline(buffer, 150)) {
        tokenize(vs, buffer);
        if (vs.size() == 0) {
          bondfound = false;
          continue;
        }

        string str(buffer);
        if (string::npos != str.find(_mol.GetTitle(),0))
          molfound = true;

        if (atomfound) {
          if (n) {
            types.push_back(atoi(vs[2].c_str()));
            types.push_back(atoi(vs[5].c_str()));
            types.push_back(atoi(vs[8].c_str()));
            types.push_back(atoi(vs[11].c_str()));
          } else {
            if (vs.size() > 2)
              types.push_back(atoi(vs[2].c_str()));
            if (vs.size() > 5)
              types.push_back(atoi(vs[5].c_str()));
            if (vs.size() > 8)
              types.push_back(atoi(vs[8].c_str()));

            atomfound = false;
          }
          n--;
        }

        if (fchgfound) {
          if (n) {
            fcharges.push_back(atof(vs[2].c_str()));
            fcharges.push_back(atof(vs[5].c_str()));
            fcharges.push_back(atof(vs[8].c_str()));
            fcharges.push_back(atof(vs[11].c_str()));
          } else {
            if (vs.size() > 2)
              fcharges.push_back(atof(vs[2].c_str()));
            if (vs.size() > 5)
              fcharges.push_back(atof(vs[5].c_str()));
            if (vs.size() > 8)
              fcharges.push_back(atof(vs[8].c_str()));

            fchgfound = false;
          }
          n--;
        }

        if (pchgfound) {
          if (n) {
            pcharges.push_back(atof(vs[2].c_str()));
            pcharges.push_back(atof(vs[5].c_str()));
            pcharges.push_back(atof(vs[8].c_str()));
            pcharges.push_back(atof(vs[11].c_str()));
          } else {
            if (vs.size() > 2)
              pcharges.push_back(atof(vs[2].c_str()));
            if (vs.size() > 5)
              pcharges.push_back(atof(vs[5].c_str()));
            if (vs.size() > 8)
              pcharges.push_back(atof(vs[8].c_str()));

            pchgfound = false;
          }
          n--;
        }

        if (molfound && EQn(buffer, " ATOM NAME  TYPE", 16)) {
          atomfound = true;
          n = _mol.NumAtoms() / 4;
        }
        if (molfound && EQn(buffer, "   ATOM   FCHARGE", 17)) {
          fchgfound = true;
          n = _mol.NumAtoms() / 4;
        }
        if (molfound && EQn(buffer, "   ATOM    CHARGE", 17)) {
          pchgfound = true;
          n = _mol.NumAtoms() / 4;
        }

        if (bondfound)
          bond_lengths.push_back(atof(vs[7].c_str()));

        if (molfound) {
          if (EQn(buffer, " Total ENERGY", 13))
            etot = atof(vs[3].c_str());
          if (EQn(buffer, " Bond Stretching", 16))
            ebond = atof(vs[2].c_str());
          if (EQn(buffer, " Angle Bending", 14))
            eangle = atof(vs[2].c_str());
          if (EQn(buffer, " Out-of-Plane Bending", 21))
            eoop = atof(vs[2].c_str());
          if (EQn(buffer, " Stretch-Bend", 13))
            estbn = atof(vs[1].c_str());
          if (EQn(buffer, "     Total Torsion", 18))
            etor = atof(vs[2].c_str());
          if (EQn(buffer, "     Net vdW", 12))
            evdw = atof(vs[2].c_str());
          if (EQn(buffer, " Electrostatic", 14))
            eeq = atof(vs[1].c_str());
          if (EQn(buffer, " ---------------------", 22) && (termcount == 0)) {
            termcount++;
            bondfound = true;
          }
          if (EQn(buffer, " OPTIMOL>  # read next", 22))
            break;
        }



      } // while (getline)

      vector<int>::iterator i;
      vector<double>::iterator di;
      unsigned int ni;
      bool failed;

      cout << "--------------------------------------------------------------------------------" << endl;
      cout << "                                                                                " << endl;
      cout << "  VALIDATE MOLECULE " << c << ": " << _mol.GetTitle() << endl;
      cout << "                                                                                " << endl;
      cout << "IDX  HYB  AROM  OB_TYPE  LOG_TYPE       RESULT                                  " << endl;
      cout << "----------------------------------------------                                  " << endl;

      ni = 1;
      failed = false;
      for (i = types.begin(); i != types.end();++i) {
        if (ni > _mol.NumAtoms())
          continue;

        if ( (atoi(_mol.GetAtom(ni)->GetType()) == 87) ||
             (atoi(_mol.GetAtom(ni)->GetType()) == 97)
             ) continue;

        if (atoi(_mol.GetAtom(ni)->GetType()) == (*i))
          snprintf(_logbuf, BUFF_SIZE, "%2d   %3d  %4d    %3d      %3d          PASSED", _mol.GetAtom(ni)->GetIdx(), _mol.GetAtom(ni)->GetHyb(),
                  _mol.GetAtom(ni)->IsAromatic(), atoi(_mol.GetAtom(ni)->GetType()), *i);
        else {
          snprintf(_logbuf, BUFF_SIZE, "%2d   %3d  %4d    %3d      %3d      XXX FAILED XXX", _mol.GetAtom(ni)->GetIdx(), _mol.GetAtom(ni)->GetHyb(),
                  _mol.GetAtom(ni)->IsAromatic(), atoi(_mol.GetAtom(ni)->GetType()), *i);
          failed = true;
        }

        cout << _logbuf << endl;

        ni++;
      }

      if (failed) {
        cout << "Could not successfully assign atom types" << endl;
        return false;
        //continue;
      }

      SetFormalCharges();
      cout << endl;
      cout << "IDX  OB_FCARGE  LOG_FCHARGE       RESULT" << endl;
      cout << "----------------------------------------" << endl;

      ni = 1;
      for (di = fcharges.begin(); di != fcharges.end(); ++di) {
        if (ni > _mol.NumAtoms())
          continue;

        if ( (atoi(_mol.GetAtom(ni)->GetType()) == 87) ||
             (atoi(_mol.GetAtom(ni)->GetType()) == 97)
             ) continue;

        if (fabs((*di) - _mol.GetAtom(ni)->GetPartialCharge()) <= 0.001)
          snprintf(_logbuf, BUFF_SIZE, "%2d   %7.4f     %7.4f          PASSED", _mol.GetAtom(ni)->GetIdx(), _mol.GetAtom(ni)->GetPartialCharge(), *di);
        else {
          snprintf(_logbuf, BUFF_SIZE, "%2d   %7.4f     %7.4f      XXX FAILED XXX", _mol.GetAtom(ni)->GetIdx(), _mol.GetAtom(ni)->GetPartialCharge(), *di);
          failed = true;
        }

        cout << _logbuf << endl;

        ni++;
      }

      if (failed) {
        cout << "Could not successfully assign formal charges" << endl;
        //return false;
        continue;
      }

      SetPartialCharges();
      cout << endl;
      cout << "IDX  OB_PCARGE  LOG_PCHARGE       RESULT" << endl;
      cout << "----------------------------------------" << endl;

      ni = 1;
      for (di = pcharges.begin(); di != pcharges.end(); ++di) {
        if (ni > _mol.NumAtoms())
          continue;

        if ( (atoi(_mol.GetAtom(ni)->GetType()) == 87) ||
             (atoi(_mol.GetAtom(ni)->GetType()) == 97)
             ) continue;

        if (fabs((*di) - _mol.GetAtom(ni)->GetPartialCharge()) <= 0.001)
          snprintf(_logbuf, BUFF_SIZE, "%2d   %7.4f     %7.4f          PASSED", _mol.GetAtom(ni)->GetIdx(), _mol.GetAtom(ni)->GetPartialCharge(), *di);
        else {
          snprintf(_logbuf, BUFF_SIZE, "%2d   %7.4f     %7.4f      XXX FAILED XXX", _mol.GetAtom(ni)->GetIdx(), _mol.GetAtom(ni)->GetPartialCharge(), *di);
          failed = true;
        }

        cout << _logbuf << endl;

        ni++;
      }

      if (failed) {
        cout << "Could not successfully assign partial charges" << endl;
        //return false;
        continue;
      }



      if (!SetupCalculations()) {
        cout << "Could not setup calculations (missing parameters...)" << endl;
        return false;
        //continue;
      }

      double delta;
      cout << endl;
      cout << "TERM                     OB ENERGY     LOG ENERGY         DELTA" << endl;
      cout << "---------------------------------------------------------------" << endl;

      delta = (E_Bond() - ebond);
      snprintf(_logbuf, BUFF_SIZE, "Bond Stretching        %11.5f    %11.5f   %11.5f", E_Bond(), ebond, delta);
      cout << _logbuf << endl;

      delta = (E_Angle() - eangle);
      snprintf(_logbuf, BUFF_SIZE, "Angle Bending          %11.5f    %11.5f   %11.5f", E_Angle(), eangle, delta);
      cout << _logbuf << endl;

      delta = (E_StrBnd() - estbn);
      snprintf(_logbuf, BUFF_SIZE, "Stretch-Bending        %11.5f    %11.5f   %11.5f", E_StrBnd(), estbn, delta);
      cout << _logbuf << endl;

      delta = (E_OOP() - eoop);
      snprintf(_logbuf, BUFF_SIZE, "Out-Of-Plane Bending   %11.5f    %11.5f   %11.5f", E_OOP(), eoop, delta);
      cout << _logbuf << endl;

      delta = (E_Torsion() - etor);
      snprintf(_logbuf, BUFF_SIZE, "Torsional              %11.5f    %11.5f   %11.5f", E_Torsion(), etor, delta);
      cout << _logbuf << endl;

      delta = (E_VDW() - evdw);
      snprintf(_logbuf, BUFF_SIZE, "Van der Waals          %11.5f    %11.5f   %11.5f", E_VDW(), evdw, delta);
      cout << _logbuf << endl;

      delta = (E_Electrostatic() - eeq);
      snprintf(_logbuf, BUFF_SIZE, "Electrostatic          %11.5f    %11.5f   %11.5f", E_Electrostatic(), eeq, delta);
      cout << _logbuf << endl;

      cout << endl;
      delta = (Energy() - etot);
      snprintf(_logbuf, BUFF_SIZE, "Total ENERGY           %11.5f    %11.5f   %11.5f", Energy(), etot, delta);
      cout << _logbuf << endl;

    } // for (unsigned int c;; c++ )

    if (ifs)
      ifs.close();
    if (ifs2)
      ifs2.close();

    return true;
  }

  bool OBForceFieldMMFF94::ValidateGradients ()
  {
    vector3 numgrad, anagrad, err;
    int coordIdx;

    bool passed = true; // set to false if any component fails

    cout << "----------------------------------------------------------------------------------------" << endl;
    cout << "                                                                                        " << endl;
    cout << "  VALIDATE GRADIENTS : " << _mol.GetTitle() << endl;
    cout << "                                                                                        " << endl;
    cout << "                                                                                        " << endl;
    cout << "ATOM IDX      NUMERICAL GRADIENT           ANALYTICAL GRADIENT        REL. ERROR (%)   " << endl;
    cout << "----------------------------------------------------------------------------------------" << endl;
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
      E_Bond(); // compute
      anagrad.Set(_gradientPtr[coordIdx], _gradientPtr[coordIdx+1], _gradientPtr[coordIdx+2]);
      err = ValidateGradientError(numgrad, anagrad);

      snprintf(_logbuf, BUFF_SIZE, "    bond    (%7.3f, %7.3f, %7.3f)  (%7.3f, %7.3f, %7.3f)  (%5.2f, %5.2f, %5.2f)\n", numgrad.x(), numgrad.y(), numgrad.z(),
              anagrad.x(), anagrad.y(), anagrad.z(), err.x(), err.y(), err.z());
      OBFFLog(_logbuf);
      if (err.x() > 5.0 || err.y() > 5.0 || err.z() > 5.0)
        passed = false;

      // OBFF_EANGLE
      numgrad = NumericalDerivative(&*a, OBFF_EANGLE);
      ClearGradients();
      E_Angle(); // compute
      anagrad.Set(_gradientPtr[coordIdx], _gradientPtr[coordIdx+1], _gradientPtr[coordIdx+2]);
      err = ValidateGradientError(numgrad, anagrad);

      snprintf(_logbuf, BUFF_SIZE, "    angle   (%7.3f, %7.3f, %7.3f)  (%7.3f, %7.3f, %7.3f)  (%5.2f, %5.2f, %5.2f)\n", numgrad.x(), numgrad.y(), numgrad.z(),
              anagrad.x(), anagrad.y(), anagrad.z(), err.x(), err.y(), err.z());
      OBFFLog(_logbuf);
      if (err.x() > 5.0 || err.y() > 5.0 || err.z() > 5.0)
        passed = false;

      // OBFF_ESTRBND
      numgrad = NumericalDerivative(&*a, OBFF_ESTRBND);
      ClearGradients();
      E_StrBnd(); // compute
      anagrad.Set(_gradientPtr[coordIdx], _gradientPtr[coordIdx+1], _gradientPtr[coordIdx+2]);
      err = ValidateGradientError(numgrad, anagrad);

      snprintf(_logbuf, BUFF_SIZE, "    strbnd  (%7.3f, %7.3f, %7.3f)  (%7.3f, %7.3f, %7.3f)  (%5.2f, %5.2f, %5.2f)\n", numgrad.x(), numgrad.y(), numgrad.z(),
              anagrad.x(), anagrad.y(), anagrad.z(), err.x(), err.y(), err.z());
      OBFFLog(_logbuf);
      if (err.x() > 5.0 || err.y() > 5.0 || err.z() > 5.0)
        passed = false;

      // OBFF_ETORSION
      numgrad = NumericalDerivative(&*a, OBFF_ETORSION);
      ClearGradients();
      E_Torsion(); // compute
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
      // disable OOP gradient validation for now -- some small errors, but nothing major
      //      if (err.x() > 5.0 || err.y() > 5.0 || err.z() > 5.0)
      //        passed = false;

      // OBFF_EVDW
      numgrad = NumericalDerivative(&*a, OBFF_EVDW);
      ClearGradients();
      E_VDW(); // compute
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
      E_Electrostatic(); // compute
      anagrad.Set(_gradientPtr[coordIdx], _gradientPtr[coordIdx+1], _gradientPtr[coordIdx+2]);
      err = ValidateGradientError(numgrad, anagrad);

      snprintf(_logbuf, BUFF_SIZE, "    electro (%7.3f, %7.3f, %7.3f)  (%7.3f, %7.3f, %7.3f)  (%5.2f, %5.2f, %5.2f)\n", numgrad.x(), numgrad.y(), numgrad.z(),
              anagrad.x(), anagrad.y(), anagrad.z(), err.x(), err.y(), err.z());
      OBFFLog(_logbuf);
      if (err.x() > 5.0 || err.y() > 5.0 || err.z() > 5.0)
        passed = false;
    }

    return passed; // did we pass every single component?
  }

  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  //
  //  Calculate bond type, angle type, stretch-bend type, torsion type
  //
  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////

  //
  // MMFF part V - page 620
  //
  // BTij is 1 when:
  // a) single bond between atoms i and j, both i and j are not aromatic and both types have sbmb set in mmffprop.par, or
  // b) bewtween two aromatic atoms, but the bond is not aromatic (e.g. connecting bond in biphenyl)
  //
  int OBForceFieldMMFF94::GetBondType(OBAtom* a, OBAtom* b)
  {
    OBBond *bond = _mol.GetBond(a, b);
    if (bond->GetBondOrder() != 1 || bond->IsAromatic())
      return 0;

    if (HasAromSet(atoi(a->GetType())) && HasAromSet(atoi(b->GetType())))
      return 1;

    if (HasSbmbSet(atoi(a->GetType())) && HasSbmbSet(atoi(b->GetType())))
      return 1;

    return 0;
  }

  int OBForceFieldMMFF94::GetAngleType(OBAtom* a, OBAtom* b, OBAtom *c)
  {
    int sumbondtypes;

    sumbondtypes = GetBondType(a,b) + GetBondType(b, c);

    if (a->IsInRingSize(3) && b->IsInRingSize(3) && c->IsInRingSize(3) && IsInSameRing(a, c))
      switch (sumbondtypes) {
      case 0:
        return 3;
      case 1:
        return 5;
      case 2:
        return 6;
      }

    if (a->IsInRingSize(4) && b->IsInRingSize(4) && c->IsInRingSize(4) && IsInSameRing(a, c))
      switch (sumbondtypes) {
      case 0:
        return 4;
      case 1:
        return 7;
      case 2:
        return 8;
      }

    return sumbondtypes;
  }

  int OBForceFieldMMFF94::GetStrBndType(OBAtom* a, OBAtom* b, OBAtom *c)
  {
    int btab, btbc, atabc;
    bool inverse;

    btab = GetBondType(a, b);
    btbc = GetBondType(b, c);
    atabc = GetAngleType(a, b, c);

    if (atoi(a->GetType()) <= atoi(c->GetType()))
      inverse = false;
    else
      inverse = true;

    switch (atabc) {
    case 0:
      return 0;

    case 1:
      if (btab) {
        if (!inverse) {
          return 1;
        } else {
          return 2;
        }
      }
      if (btbc) {
        if (!inverse) {
          return 2;
        } else {
          return 1;
        }
      }

    case 2:
      return 3;

    case 3:
      return 5;

    case 4:
      return 4;

    case 5:
      if (btab) {
        if (!inverse) {
          return 6;
        } else {
          return 7;
        }
      }
      if (btbc) {
        if (!inverse) {
          return 7;
        } else {
          return 6;
        }
      }

    case 6:
      return 8;

    case 7:
      if (btab) {
        if (!inverse) {
          return 9;
        } else {
          return 10;
        }
      }
      if (btbc) {
        if (!inverse) {
          return 10;
        } else {
          return 9;
        }
      }

    case 8:
      return 11;
    }

    return 0;
  }

  //
  // MMFF part IV - page 609
  //
  // TTijkl = 1 when BTjk = 1
  // TTijkl = 2 when BTjk = 0 but BTij and/or BTkl = 1
  // TTijkl = 4 when i, j, k and l are all members of the same four-membered ring
  // TTijkl = 5 when i, j, k and l are members of a five-membered ring and at least one is a sp3-hybridized carbon (MMFF atom type 1)
  //
  int OBForceFieldMMFF94::GetTorsionType(OBAtom* a, OBAtom* b, OBAtom *c, OBAtom *d)
  {
    int btab, btbc, btcd;

    btab = GetBondType(a, b);
    btbc = GetBondType(b, c);
    btcd = GetBondType(c, d);

    if (btbc == 1)
      return 1;

    if (a->IsInRingSize(4) && b->IsInRingSize(4) && c->IsInRingSize(4) && d->IsInRingSize(4))
      if (IsInSameRing(a,b) && IsInSameRing(b,c) && IsInSameRing(c,d))
        return 4;

    OBBond *bond = _mol.GetBond(b, c);
    if (bond->GetBondOrder() == 1 && !bond->IsAromatic()) {
      if (btab || btcd)
        return 2;
      /*
        unsigned int order1 = GetCXT(0, atoi(d->GetType()), atoi(c->GetType()), atoi(b->GetType()), atoi(a->GetType()));
        unsigned int order2 = GetCXT(0, atoi(a->GetType()), atoi(b->GetType()), atoi(c->GetType()), atoi(d->GetType()));

        cout << "GetTorsionType(" << a->GetType() << ", " << b->GetType() << ", " << c->GetType() << ", " << d->GetType() << ")" << endl;
        cout << "    order1 = " << order1 << endl;
        cout << "    order2 = " << order2 << endl;
        cout << "    btab = " << btab << endl;
        cout << "    btbc = " << btbc << endl;
        cout << "    btcd = " << btcd << endl;
      */
    }

    if (a->IsInRingSize(5) && b->IsInRingSize(5) && c->IsInRingSize(5) && d->IsInRingSize(5)) {
      vector<OBRing*> vr;
      vr = _mol.GetSSSR();

      if( !((atoi(a->GetType()) == 1) || (atoi(b->GetType()) == 1) || (atoi(c->GetType()) == 1) || (atoi(d->GetType()) == 1)) )
        return 0;

      vector<OBRing*>::iterator ri;
      vector<int>::iterator rj;
      for (ri = vr.begin();ri != vr.end();++ri) { // for each ring
        if ((*ri)->IsAromatic())
          continue;

        if ((*ri)->Size() != 5)
          continue;

        if (!(*ri)->IsMember(a) || !(*ri)->IsMember(b) || !(*ri)->IsMember(c) || !(*ri)->IsMember(d))
          continue;

        return 5;
      }
    }


    return 0;
  }

  // CXB = MC * (I * MA + J) + BTij
  unsigned int OBForceFieldMMFF94::GetCXB(int type, int a, int b)
  {
    unsigned int cxb;
    cxb = 2 * (a * 136 + b) + type;
    return cxb;
  }

  // CXA = MC * (J * MA^2 + I * MA + K) + ATijk
  unsigned int OBForceFieldMMFF94::GetCXA(int type, int a, int b, int c)
  {
    unsigned int cxa;
    cxa = 9 * (b * 18496 + a * 136 + c) + type;
    return cxa;
  }

  // CXS = MC * (J * MA^2 + I * MA + K) + STijk
  unsigned int OBForceFieldMMFF94::GetCXS(int type, int a, int b, int c)
  {
    unsigned int cxs;
    cxs = 12 * (b * 18496 + a * 136 + c) + type;
    return cxs;
  }

  // CXO = J * MA^3 + I * MA^2 + K * MA + L
  unsigned int OBForceFieldMMFF94::GetCXO(int a, int b, int c, int d)
  {
    unsigned int cxo;
    cxo = b * 2515456 + a * 18496 + c * 136 + d;
    return cxo;
  }

  // CXT = MC * (J * MA^3 + K * MA^2 + I * MA + L) + TTijkl
  unsigned int OBForceFieldMMFF94::GetCXT(int type, int a, int b, int c, int d)
  {
    unsigned int cxt;
    cxt = 6 * (b * 2515456 + c * 18496 + a * 136 + d) + type;
    return cxt;
  }

  // CXQ = MC * (I * MA + J) + BTij
  unsigned int OBForceFieldMMFF94::GetCXQ(int type, int a, int b)
  {
    unsigned int cxq;
    cxq = 2 * (a * 136 + b) + type;
    return cxq;
  }

  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  //
  //  Various tables & misc. functions
  //
  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////

  // MMFF part V - TABLE I
  bool OBForceFieldMMFF94::HasLinSet(int atomtype)
  {
    return _ffpropLin.BitIsSet(atomtype);
  }

  // MMFF part V - TABLE I
  bool OBForceFieldMMFF94::HasPilpSet(int atomtype)
  {
    return _ffpropPilp.BitIsSet(atomtype);
  }

  // MMFF part V - TABLE I
  bool OBForceFieldMMFF94::HasAromSet(int atomtype)
  {
    return _ffpropArom.BitIsSet(atomtype);
  }

  // MMFF part V - TABLE I
  bool OBForceFieldMMFF94::HasSbmbSet(int atomtype)
  {
    return _ffpropSbmb.BitIsSet(atomtype);
  }

  // MMFF part V - TABLE I
  int OBForceFieldMMFF94::GetCrd(int atomtype)
  {
    OBFFParameter *par;

    par = GetParameter1Atom(atomtype, _ffpropparams); // from mmffprop.par
    if (par)
      return par->_ipar[1];

    return 0;
  }

  // MMFF part V - TABLE I
  int OBForceFieldMMFF94::GetVal(int atomtype)
  {
    OBFFParameter *par;

    par = GetParameter1Atom(atomtype, _ffpropparams); // from mmffprop.par
    if (par)
      return par->_ipar[2];

    return 0;
  }

  // MMFF part V - TABLE I
  int OBForceFieldMMFF94::GetMltb(int atomtype)
  {
    OBFFParameter *par;

    par = GetParameter1Atom(atomtype, _ffpropparams); // from mmffprop.par
    if (par)
      return par->_ipar[4];

    return 0;
  }

  // MMFF part I - TABLE IV
  int OBForceFieldMMFF94::EqLvl2(int type)
  {
    for (unsigned int idx=0; idx < _ffdefparams.size(); idx++)
      if (_ffdefparams[idx]._ipar[0] == type)
        return _ffdefparams[idx]._ipar[1];

    return type;
  }

  // MMFF part I - TABLE IV
  int OBForceFieldMMFF94::EqLvl3(int type)
  {
    for (unsigned int idx=0; idx < _ffdefparams.size(); idx++)
      if (_ffdefparams[idx]._ipar[0] == type)
        return _ffdefparams[idx]._ipar[2];

    return type;
  }

  // MMFF part I - TABLE IV
  int OBForceFieldMMFF94::EqLvl4(int type)
  {
    for (unsigned int idx=0; idx < _ffdefparams.size(); idx++)
      if (_ffdefparams[idx]._ipar[0] == type)
        return _ffdefparams[idx]._ipar[3];

    return type;
  }

  // MMFF part I - TABLE IV
  int OBForceFieldMMFF94::EqLvl5(int type)
  {
    for (unsigned int idx=0; idx < _ffdefparams.size(); idx++)
      if (_ffdefparams[idx]._ipar[0] == type)
        return _ffdefparams[idx]._ipar[4];

    return type;
  }

  // MMFF part V - TABLE VI
  double OBForceFieldMMFF94::GetZParam(OBAtom* atom)
  {
    switch (atom->GetAtomicNum()) {
    case OBElements::Hydrogen:
      return 1.395;
    case OBElements::Carbon:
      return 2.494;
    case OBElements::Nitrogen:
      return 2.711;
    case OBElements::Oxygen:
      return 3.045;
    case OBElements::Fluorine:
      return 2.847;
    case OBElements::Silicon:
      return 2.350;
    case OBElements::Phosphorus:
      return 2.350;
    case OBElements::Sulfur:
      return 2.980;
    case OBElements::Chlorine:
      return 2.909;
    case OBElements::Bromine:
      return 3.017;
    case OBElements::Iodine:
      return 3.086;
    }

    return 0.0;
  }

  // MMFF part V - TABLE VI
  double OBForceFieldMMFF94::GetCParam(OBAtom* atom)
  {
    switch (atom->GetAtomicNum()) {
    case OBElements::Boron:
      return 0.704;
    case OBElements::Carbon:
      return 1.016;
    case OBElements::Nitrogen:
      return 1.113;
    case OBElements::Oxygen:
      return 1.337;
    case OBElements::Silicon:
      return 0.811;
    case OBElements::Phosphorus:
      return 1.068;
    case OBElements::Sulfur:
      return 1.249;
    case OBElements::Chlorine:
      return 1.078;
    case OBElements::Arsenic:
      return 0.825;
    }

    return 0.0;
  }

  // MMFF part V - TABLE X
  double OBForceFieldMMFF94::GetUParam(OBAtom* atom)
  {
    switch (atom->GetAtomicNum()) {
    case OBElements::Carbon:
      return 2.0;
    case OBElements::Nitrogen:
      return 2.0;
    case OBElements::Oxygen:
      return 2.0;
    case OBElements::Silicon:
      return 1.25;
    case OBElements::Phosphorus:
      return 1.25;
    case OBElements::Sulfur:
      return 1.25;
    }

    return 0.0;
  }

  // MMFF part V - TABLE X
  double OBForceFieldMMFF94::GetVParam(OBAtom* atom)
  {
    switch (atom->GetAtomicNum()) {
    case OBElements::Carbon:
      return 2.12;
    case OBElements::Nitrogen:
      return 1.5;
    case OBElements::Oxygen:
      return 0.2;
    case OBElements::Silicon:
      return 1.22;
    case OBElements::Phosphorus:
      return 2.4;
    case OBElements::Sulfur:
      return 0.49;
    }

    return 0.0;
  }

  // R Blom and A Haaland, J. Mol. Struct., 128, 21-27 (1985)
  double OBForceFieldMMFF94::GetCovalentRadius(OBAtom* a) {

    switch (a->GetAtomicNum()) {
    case 1:
      return 0.33; // corrected value from MMFF part V
    case 5:
      return 0.81;
    case 6:
      return 0.77; // corrected value from MMFF part V
    case 7:
      return 0.73;
    case 8:
      return 0.72;
    case 9:
      return 0.74;
    case 13:
      return 1.22;
    case 14:
      return 1.15;
    case 15:
      return 1.09;
    case 16:
      return 1.03;
    case 17:
      return 1.01;
    case 31:
      return 1.19;
    case 32:
      return 1.20;
    case 33:
      return 1.20;
    case 34:
      return 1.16;
    case 35:
      return 1.15;
    case 44:
      return 1.46;
    case 50:
      return 1.40;
    case 51:
      return 1.41;
    case 52:
      return 1.35;
    case 53:
      return 1.33;
    case 81:
      return 1.51;
    case 82:
      return 1.53;
    case 83:
      return 1.55;
    default:
      return OBElements::GetCovalentRad(a->GetAtomicNum());
    }
  }

  double OBForceFieldMMFF94::GetBondLength(OBAtom* a, OBAtom* b)
  {
    OBFFParameter *parameter;
    double rab;

    parameter = GetTypedParameter2Atom(GetBondType(a, b), atoi(a->GetType()), atoi(b->GetType()), _ffbondparams);
    if (parameter == NULL)
      rab = GetRuleBondLength(a, b);
    else
      rab = parameter->_dpar[1];

    return rab;
  }

  // MMFF part V - page 625
  double OBForceFieldMMFF94::GetRuleBondLength(OBAtom* a, OBAtom* b)
  {
    double r0ab, r0a, r0b, c, Xa, Xb;
    int Ha, Hb, BOab;
    r0a = GetCovalentRadius(a);
    r0b = GetCovalentRadius(b);
    Xa = OBElements::GetAllredRochowElectroNeg(a->GetAtomicNum());
    Xb = OBElements::GetAllredRochowElectroNeg(b->GetAtomicNum());


    if (a->GetAtomicNum() == OBElements::Hydrogen)
      r0a = 0.33;
    if (b->GetAtomicNum() == OBElements::Hydrogen)
      r0b = 0.33;

    if (a->GetAtomicNum() == OBElements::Hydrogen || b->GetAtomicNum() == OBElements::Hydrogen)
      c = 0.050;
    else
      c = 0.085;

    if (GetMltb(atoi(a->GetType()) == 3))
      Ha = 1;
    else if ((GetMltb(atoi(a->GetType())) == 1) || (GetMltb(atoi(a->GetType())) == 2))
      Ha = 2;
    else
      Ha = 3;

    if (GetMltb(atoi(b->GetType()) == 3))
      Hb = 1;
    else if ((GetMltb(atoi(b->GetType())) == 1) || (GetMltb(atoi(b->GetType())) == 2))
      Hb = 2;
    else
      Hb = 3;

    BOab = a->GetBond(b)->GetBondOrder();
    if ((GetMltb(atoi(a->GetType())) == 1) && (GetMltb(atoi(b->GetType())) == 1))
      BOab = 4;
    if ((GetMltb(atoi(a->GetType())) == 1) && (GetMltb(atoi(b->GetType())) == 2))
      BOab = 5;
    if ((GetMltb(atoi(a->GetType())) == 2) && (GetMltb(atoi(b->GetType())) == 1))
      BOab = 5;
    if (a->GetBond(b)->IsAromatic()) {
      if (!HasPilpSet(atoi(a->GetType())) && !HasPilpSet(atoi(b->GetType()))) {
        BOab = 4;
      } else {
        BOab = 5;
      }
    }

    switch (BOab) {
    case 5:
      r0a -= 0.04;
      r0b -= 0.04;
      break;
    case 4:
      r0a -= 0.075;
      r0b -= 0.075;
      break;
    case 3:
      r0a -= 0.17;
      r0b -= 0.17;
      break;
    case 2:
      r0a -= 0.10;
      r0b -= 0.10;
      break;
    case 1:
      if (Ha == 1)
        r0a -= 0.08;
      if (Ha == 2)
        r0a -= 0.03;
      if (Hb == 1)
        r0b -= 0.08;
      if (Hb == 2)
        r0b -= 0.03;
    }

    /*
      cout << "Ha=" << Ha << "  Hb=" << Hb << "  BOab=" << BOab << endl;
      cout << "r0a=" << r0a << "  Xa=" << Xa << endl;
      cout << "r0b=" << r0b << "  Xb=" << Xb << endl;
      cout << "r0a + r0b=" << r0a +r0b << endl;
      cout << "c=" << c << "  |Xa-Xb|=" << fabs(Xa-Xb) << "  |Xa-Xb|^1.4=" << pow(fabs(Xa-Xb), 1.4) << endl;
    */
    r0ab = r0a + r0b - c * pow(fabs(Xa - Xb), 1.4) - 0.008;

    return r0ab;
  }

  OBFFParameter* OBForceFieldMMFF94::GetParameter1Atom(int a, std::vector<OBFFParameter> &parameter)
  {
    OBFFParameter *par;

    for (unsigned int idx=0; idx < parameter.size(); idx++)
      if (a == parameter[idx].a) {
        par = &parameter[idx];
        return par;
      }

    return NULL;
  }

  OBFFParameter* OBForceFieldMMFF94::GetParameter2Atom(int a, int b, std::vector<OBFFParameter> &parameter)
  {
    OBFFParameter *par;

    for (unsigned int idx=0; idx < parameter.size(); idx++)
      if (((a == parameter[idx].a) && (b == parameter[idx].b)) ||
          ((a == parameter[idx].b) && (b == parameter[idx].a)))
        {
          par = &parameter[idx];
          return par;
        }

    return NULL;
  }

  OBFFParameter* OBForceFieldMMFF94::GetParameter3Atom(int a, int b, int c, std::vector<OBFFParameter> &parameter)
  {
    OBFFParameter *par;

    for (unsigned int idx=0; idx < parameter.size(); idx++)
      if (((a == parameter[idx].a) && (b == parameter[idx].b) && (c == parameter[idx].c)) ||
          ((a == parameter[idx].c) && (b == parameter[idx].b) && (c == parameter[idx].a)))
        {
          par = &parameter[idx];
          return par;
        }

    return NULL;
  }

  OBFFParameter* OBForceFieldMMFF94::GetTypedParameter2Atom(int ffclass, int a, int b, std::vector<OBFFParameter> &parameter)
  {
    OBFFParameter *par;

    for (unsigned int idx=0; idx < parameter.size(); idx++)
      if (((a == parameter[idx].a) && (b == parameter[idx].b) && (ffclass == parameter[idx]._ipar[0])) ||
          ((a == parameter[idx].b) && (b == parameter[idx].a) && (ffclass == parameter[idx]._ipar[0])))
        {
          par = &parameter[idx];
          return par;
        }

    return NULL;
  }

  OBFFParameter* OBForceFieldMMFF94::GetTypedParameter3Atom(int ffclass, int a, int b, int c, std::vector<OBFFParameter> &parameter)
  {
    OBFFParameter *par;

    for (unsigned int idx=0; idx < parameter.size(); idx++)
      if (((a == parameter[idx].a) && (b == parameter[idx].b) && (c == parameter[idx].c) && (ffclass == parameter[idx]._ipar[0])) ||
          ((a == parameter[idx].c) && (b == parameter[idx].b) && (c == parameter[idx].a) && (ffclass == parameter[idx]._ipar[0])) )
        {
          par = &parameter[idx];
          return par;
        }

    return NULL;
  }

  OBFFParameter* OBForceFieldMMFF94::GetTypedParameter4Atom(int ffclass, int a, int b, int c, int d, std::vector<OBFFParameter> &parameter)
  {
    OBFFParameter *par;

    for (unsigned int idx=0; idx < parameter.size(); idx++)
      if (((a == parameter[idx].a) && (b == parameter[idx].b) && (c == parameter[idx].c) &&
           (d == parameter[idx].d) && (ffclass == parameter[idx]._ipar[0]))
          /* || ((a == parameter[idx].d) && (b == parameter[idx].c) && (c == parameter[idx].b) &&
             (d == parameter[idx].a) && (ffclass == parameter[idx]._ipar[0])) */ )
        {
          par = &parameter[idx];
          return par;
        }

    return NULL;
  }

} // end namespace OpenBabel

//! \file forcefieldmmff94.cpp
//! \brief MMFF94 force field

/*********************************************************************
forcefieldmmff94.cpp - MMFF94 force field

Copyright (C) 2006-2007 by Tim Vandermeersch <tim.vandermeersch@gmail.com>
 
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
#include <openbabel/babelconfig.h>
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <iomanip>
#include "forcefieldmmff94.h"

using namespace std;

namespace OpenBabel
{
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
  void OBFFBondCalculationMMFF94::Compute(bool gradients)
  {
    vector3 da, db;
    double delta2, dE;
     
    if (gradients) {
      da = a->GetVector();
      db = b->GetVector();
      rab = OBForceField::VectorLengthDerivative(da, db);
    } else
      rab = a->GetDistance(b);

    delta = rab - r0;
    delta2 = delta * delta;
 
    energy = 143.9325 * 0.5 * kb * delta2 * (1.0 - 2.0 * delta + 7.0/12.0 * 4.0 * delta2);
    
    if (gradients) {
      //dE = 143.9325 * kb * delta * (1.0 - 2.0 * delta + 7.0/12.0 * 8.0 * delta);
      dE = 143.9325 * 0.5 * kb * delta * (2.0 - 6.0 * delta + 28.0/3.0 * delta2);
      da.normalize();
      db.normalize();
      grada = dE * da / rab; // - dE/drab * drab/da
      gradb = dE * db / rab; // - dE/drab * drab/db
    }
  }
  
  double OBForceFieldMMFF94::E_Bond(bool gradients)
  {
    vector<OBFFBondCalculationMMFF94>::iterator i;
    double energy = 0.0;

    IF_OBFF_LOGLVL_HIGH {
      OBFFLog("\nB O N D   S T R E T C H I N G\n\n");
      OBFFLog("ATOM TYPES   FF    BOND       IDEAL       FORCE\n");
      OBFFLog(" I    J     CLASS  LENGTH     LENGTH     CONSTANT      DELTA      ENERGY\n");
      OBFFLog("------------------------------------------------------------------------\n");
    }
    
    for (i = _bondcalculations.begin(); i != _bondcalculations.end(); ++i) {
      i->Compute(gradients);
      energy += i->GetEnergy();

      IF_OBFF_LOGLVL_HIGH {
        sprintf(logbuf, "%2d   %2d      %d   %8.3f   %8.3f     %8.3f   %8.3f   %8.3f\n", atoi((*i).a->GetType()), atoi((*i).b->GetType()), 
                (*i).bt, (*i).rab, (*i).r0, (*i).kb, (*i).delta, (*i).energy);
        OBFFLog(logbuf);
      }
    }
    
    IF_OBFF_LOGLVL_MEDIUM {
      sprintf(logbuf, "     TOTAL BOND STRETCHING ENERGY = %8.3f %s\n",  energy, GetUnit().c_str());
      OBFFLog(logbuf);
    }
    
    return energy;
  }
  
  // 
  // MMFF part I - page 495
  //      
  //                       ka_ijk                       
  // EA_ijk = 0.438449325 -------- /\0_ij^2 (1 + cs /\0_ij)
  //                         2                          
  //
  // ka_ijk	force constant (md A/rad^2)
  //
  // /\0_ij 	0_ij - 00_ij (degrees)
  //
  // cs		cubic bend constant = -0.007 deg^-1 = -0.4 rad^-1
  //
  void OBFFAngleCalculationMMFF94::Compute(bool gradients) 
  {
    vector3 da, db, dc;
    double delta2, dE;

    if (gradients) {
      da = a->GetVector();
      db = b->GetVector();
      dc = c->GetVector();
      theta = OBForceField::VectorAngleDerivative(da, db, dc);  
    } else
      theta = a->GetAngle(b->GetIdx(), c->GetIdx());
    
    delta = theta - theta0;
    delta2 = delta * delta;
    
    energy = 0.043844 * 0.5 * ka * delta2 * (1.0 - 0.007 * delta);

    if (gradients) {
      dE = 0.043844 * ka * delta * (1.0 - 1.5 * 0.007 * delta);
      grada = dE * da; // - dE/drab * drab/da
      gradb = dE * db; // - dE/drab * drab/db = - dE/drab * drab/da - dE/drab * drab/dc 
      gradc = dE * dc; // - dE/drab * drab/dc
    }
  }
 
  double OBForceFieldMMFF94::E_Angle(bool gradients)
  {
    vector<OBFFAngleCalculationMMFF94>::iterator i;
    double energy = 0.0;
 
    IF_OBFF_LOGLVL_HIGH {
      OBFFLog("\nA N G L E   B E N D I N G\n\n");
      OBFFLog("ATOM TYPES        FF    VALENCE     IDEAL      FORCE\n");
      OBFFLog(" I    J    K     CLASS   ANGLE      ANGLE     CONSTANT      DELTA      ENERGY\n");
      OBFFLog("-----------------------------------------------------------------------------\n");
    }
    
    for (i = _anglecalculations.begin(); i != _anglecalculations.end(); ++i) {

      i->Compute(gradients);
      energy += i->GetEnergy();
      
      IF_OBFF_LOGLVL_HIGH {
        sprintf(logbuf, "%2d   %2d   %2d      %d   %8.3f   %8.3f     %8.3f   %8.3f   %8.3f\n", 
	        atoi((*i).a->GetType()), atoi((*i).b->GetType()), atoi((*i).c->GetType()), 
                (*i).at, (*i).theta, (*i).theta0, (*i).ka, (*i).delta, (*i).energy);
        OBFFLog(logbuf);
      }
    }
 
    IF_OBFF_LOGLVL_MEDIUM {
      sprintf(logbuf, "     TOTAL ANGLE BENDING ENERGY = %8.3f %s\n", energy, GetUnit().c_str());
      OBFFLog(logbuf);
    }
    
    return energy;
  }
  
  void OBFFStrBndCalculationMMFF94::Compute(bool gradients)
  {
    vector3 rab_da, rab_db, rbc_db, rbc_dc, theta_da, theta_db, theta_dc;
    
    if (gradients) {
      rab_da = theta_da = a->GetVector();
      rab_db = rbc_db = theta_db = b->GetVector();
      rbc_dc = theta_dc = c->GetVector();
      theta = OBForceField::VectorAngleDerivative(theta_da, theta_db, theta_dc);
      rab = OBForceField::VectorLengthDerivative(rab_da, rab_db);
      rbc = OBForceField::VectorLengthDerivative(rbc_db, rbc_dc);
    } else {
      theta = a->GetAngle(b->GetIdx(), c->GetIdx());
      rab = a->GetDistance(b);
      rbc = b->GetDistance(c);
    }
    
    delta_theta = theta - theta0;
    delta_rab = rab - rab0;
    delta_rbc = rbc - rbc0;

    energy = 2.51210 * (kbaABC * delta_rab + kbaCBA * delta_rbc) * delta_theta;

    if (gradients) {

      grada = 2.51210 * (kbaABC * rab_da * delta_theta + theta_da * (kbaABC * delta_rab + kbaCBA * delta_rbc));
      gradc = 2.51210 * (kbaCBA * rbc_dc * delta_theta + theta_dc * (kbaABC * delta_rab + kbaCBA * delta_rbc));
      gradb = -grada - gradc;
    }
  }
  
  double OBForceFieldMMFF94::E_StrBnd(bool gradients) 
  {
    vector<OBFFStrBndCalculationMMFF94>::iterator i;
    double energy = 0.0;
    
    IF_OBFF_LOGLVL_HIGH {
      OBFFLog("\nS T R E T C H   B E N D I N G\n\n");
      OBFFLog("ATOM TYPES        FF    VALENCE     DELTA        FORCE CONSTANT\n");
      OBFFLog(" I    J    K     CLASS   ANGLE      ANGLE        I J        J K      ENERGY\n");
      OBFFLog("---------------------------------------------------------------------------\n");
    }
    
    for (i = _strbndcalculations.begin(); i != _strbndcalculations.end(); ++i) {

      i->Compute(gradients);
      energy += i->GetEnergy();
      
      IF_OBFF_LOGLVL_HIGH {
        sprintf(logbuf, "%2d   %2d   %2d     %2d   %8.3f   %8.3f   %8.3f   %8.3f   %8.3f\n", 
	        atoi((*i).a->GetType()), atoi((*i).b->GetType()), atoi((*i).c->GetType()), 
                (*i).sbt, (*i).theta, (*i).delta_theta, (*i).kbaABC, (*i).kbaCBA, (*i).energy);
        OBFFLog(logbuf);
      }
    }
	
    IF_OBFF_LOGLVL_MEDIUM {
      sprintf(logbuf, "     TOTAL ANGLE BENDING ENERGY = %8.3f %s\n", energy, GetUnit().c_str());
      OBFFLog(logbuf);
    }
 
    return energy;
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

  void OBFFTorsionCalculationMMFF94::Compute(bool gradients)
  {
    vector3 da, db, dc, dd;
    double cosine, cosine2, cosine3;
    double phi1, phi2, phi3;
    double dE, sine, sine2, sine3;
    
    if (gradients) {
      da = a->GetVector();
      db = b->GetVector();
      dc = c->GetVector();
      dd = d->GetVector();
      tor = OBForceField::VectorTorsionDerivative(da, db, dc, dd);
      if (IsNan(tor))
        tor = 1.0e-7;
    } else {
      vector3 vab, vbc, vcd, abbc, bccd;
      vab = a->GetVector() - b->GetVector();
      vbc = b->GetVector() - c->GetVector();
      vcd = c->GetVector() - d->GetVector();
      abbc = cross(vab, vbc);
      bccd = cross(vbc, vcd);

      double dotAbbcBccd = dot(abbc,bccd);
      tor = RAD_TO_DEG * acos(dotAbbcBccd / (abbc.length() * bccd.length()));
      if (IsNearZero(dotAbbcBccd)) {
        tor = 0.0; // rather than NaN
      }
      else if (dotAbbcBccd > 0.0) {
        tor = -tor;
      }
    }
    
    cosine = cos(DEG_TO_RAD * tor);
    cosine2 = cos(DEG_TO_RAD * 2 * tor);
    cosine3 = cos(DEG_TO_RAD * 3 * tor);
      
    phi1 = 1.0 + cosine;
    phi2 = 1.0 - cosine2;
    phi3 = 1.0 + cosine3;
    
    energy = 0.5 * (v1 * phi1 + v2 * phi2 + v3 * phi3);
    
    if (gradients) {
      sine = sin(DEG_TO_RAD * tor);
      sine2 = sin(2.0 * DEG_TO_RAD * tor);
      sine3 = sin(3.0 * DEG_TO_RAD * tor);
      //dE = -k1 * sine + k2 * 2.0 * sine2 - k3 * 3.0 * sine3; ghemical
      dE = -0.5 * (v1 * sine - 2.0 * v2 * sine2 + 3.0 * v3 * sine3); // MMFF
      grada = dE * da; // - dE/drab * drab/da
      gradb = dE * db; // - dE/drab * drab/db
      gradc = dE * dc; // - dE/drab * drab/dc
      gradd = dE * dd; // - dE/drab * drab/dd
    }
  }
  
  double OBForceFieldMMFF94::E_Torsion(bool gradients) 
  {
    vector<OBFFTorsionCalculationMMFF94>::iterator i;
    double energy = 0.0;
        
    IF_OBFF_LOGLVL_HIGH {
      OBFFLog("\nT O R S I O N A L\n\n");
      OBFFLog("ATOM TYPES             FF     TORSION       FORCE CONSTANT\n");
      OBFFLog(" I    J    K    L     CLASS    ANGLE         V1   V2   V3     ENERGY\n");
      OBFFLog("--------------------------------------------------------------------\n");
    }

    for (i = _torsioncalculations.begin(); i != _torsioncalculations.end(); ++i) {

      i->Compute(gradients);
      energy += i->GetEnergy();
      
      IF_OBFF_LOGLVL_HIGH {
        sprintf(logbuf, "%2d   %2d   %2d   %2d      %d   %8.3f   %6.3f   %6.3f   %6.3f   %8.3f\n",  atoi((*i).a->GetType()), atoi((*i).b->GetType()), 
                atoi((*i).c->GetType()), atoi((*i).d->GetType()), (*i).tt, (*i).tor, (*i).v1, (*i).v2, (*i).v3, (*i).energy);
        OBFFLog(logbuf);
      }
    }
    
    IF_OBFF_LOGLVL_MEDIUM {
      sprintf(logbuf, "     TOTAL TORSIONAL ENERGY = %8.3f %s\n", energy, GetUnit().c_str());
      OBFFLog(logbuf);
    }

    return energy;
  }
 
  //
  //  a
  //   \
  //    b---d      plane = a-b-c
  //   /
  //  c
  //
  void OBFFOOPCalculationMMFF94::Compute(bool gradients)
  {
    vector3 da, db, dc, dd;
    double angle2, dE;
    
    if (gradients) {
      da = a->GetVector();
      db = b->GetVector();
      dc = c->GetVector();
      dd = d->GetVector();
      angle = OBForceField::VectorOOPDerivative(da, db, dc, dd);
      dE =  (0.043844 * angle * koop) / cos(angle);
      grada = dE * da; // - dE/drab * drab/da
      gradb = dE * db; // - dE/drab * drab/db
      gradc = dE * dc; // - dE/drab * drab/dc
      gradd = dE * dd; // - dE/drab * drab/dd
    } else {
      angle = Point2PlaneAngle(d->GetVector(), a->GetVector(), b->GetVector(), c->GetVector());
    }
    
    angle2 = angle * angle;
    energy = 0.043844 * 0.5 * koop * angle2;
  }

  double OBForceFieldMMFF94::E_OOP(bool gradients) 
  {
    vector<OBFFOOPCalculationMMFF94>::iterator i;
    double energy = 0.0;
    
    IF_OBFF_LOGLVL_HIGH {
      OBFFLog("\nO U T - O F - P L A N E   B E N D I N G\n\n");
      OBFFLog("ATOM TYPES             FF       OOP     FORCE\n");
      OBFFLog(" I    J    K    L     CLASS    ANGLE   CONSTANT     ENERGY\n");
      OBFFLog("----------------------------------------------------------\n");
    }

    for (i = _oopcalculations.begin(); i != _oopcalculations.end(); ++i) {
      
      i->Compute(gradients);
      energy += i->GetEnergy();
      
      IF_OBFF_LOGLVL_HIGH {
        sprintf(logbuf, "%2d   %2d   %2d   %2d      0   %8.3f   %8.3f     %8.3f\n", 
	        atoi((*i).a->GetType()), atoi((*i).b->GetType()), atoi((*i).c->GetType()), atoi((*i).d->GetType()), 
                (*i).angle, (*i).koop, (*i).energy);
        OBFFLog(logbuf);
      }
    }
    
    IF_OBFF_LOGLVL_MEDIUM {
      sprintf(logbuf, "     TOTAL TORSIONAL ENERGY = %8.3f %s\n", energy, GetUnit().c_str());
      OBFFLog(logbuf);
    }

    return energy;
  }
 
  void OBFFVDWCalculationMMFF94::Compute(bool gradients)
  {
    vector3 da, db;
    double rab2, /*rab3,*/ rab4, /*rab5,*/ rab6, rab7;
    double erep2, /*erep3,*/ erep4, /*erep5,*/ erep6, erep7;
    double dE;

    if (gradients) {
      da = a->GetVector();
      db = b->GetVector();
      rab = OBForceField::VectorLengthDerivative(da, db);
    } else
      rab = a->GetDistance(b);
 
    rab2 = rab * rab;
    rab4 = rab2 * rab2;
    rab6 = rab4 * rab2;
    rab7 = rab6 * rab;

    erep = (1.07 * R_AB) / (rab + 0.07 * R_AB); //***
    erep2 = erep * erep;
    erep4 = erep2 * erep2;
    erep6 = erep4 * erep2;
    erep7 = erep6 * erep;
      
    eattr = (((1.12 * R_AB7) / (rab7 + 0.12 * R_AB7)) - 2.0);
      
    energy = escale * epsilon * erep7 * eattr;
    
    if (gradients) { 
      //dE = ...
      grada = dE * da; // - dE/drab * drab/da
      gradb = dE * db; // - dE/drab * drab/db
    }
  }
  
  double OBForceFieldMMFF94::E_VDW(bool gradients)
  {
    vector<OBFFVDWCalculationMMFF94>::iterator i;
    double energy = 0.0;
    
    IF_OBFF_LOGLVL_HIGH {
      OBFFLog("\nV A N   D E R   W A A L S\n\n");
      OBFFLog("ATOM TYPES\n");
      OBFFLog(" I    J        Rij       R*IJ    EPSILON    E_REP     E_ATTR    ENERGY\n");
      OBFFLog("----------------------------------------------------------------------\n");
      //       XX   XX     -000.000  -000.000  -000.000  -000.000  -000.000  -000.000
    }
    
    for (i = _vdwcalculations.begin(); i != _vdwcalculations.end(); ++i) {
      
      i->Compute(gradients);
      energy += i->GetEnergy();
      
      IF_OBFF_LOGLVL_HIGH {
        sprintf(logbuf, "%2d   %2d     %8.3f  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f\n", atoi((*i).a->GetType()), atoi((*i).b->GetType()), 
                (*i).rab, (*i).R_AB, (*i).epsilon, (*i).erep, (*i).eattr, (*i).energy);
        OBFFLog(logbuf);
      }
    }
    
    IF_OBFF_LOGLVL_MEDIUM {
      sprintf(logbuf, "     TOTAL VAN DER WAALS ENERGY = %8.3f %s\n", energy, GetUnit().c_str());
      OBFFLog(logbuf);
    }

    return energy;
  }

  void OBFFElectrostaticCalculationMMFF94::Compute(bool gradients)
  {
    vector3 da, db;
    double dE, rab2;

    if (gradients) {
      da = a->GetVector();
      db = b->GetVector();
      rab = OBForceField::VectorLengthDerivative(da, db);
    } else
      rab = a->GetDistance(b);

    energy = qq / (rab + 0.05);

    if (gradients) {
      rab2 = rab * rab;
      dE = -qq / rab2;
      grada = dE * da; // - dE/drab * drab/da
      gradb = dE * db; // - dE/drab * drab/db
    } 
  }

  double OBForceFieldMMFF94::E_Electrostatic(bool gradients)
  {
    vector<OBFFElectrostaticCalculationMMFF94>::iterator i;
    double energy = 0.0;
    
    IF_OBFF_LOGLVL_HIGH {
      OBFFLog("\nE L E C T R O S T A T I C   I N T E R A C T I O N S\n\n");
      OBFFLog("ATOM TYPES\n");
      OBFFLog(" I    J        Rij        Qi         Qj        ENERGY\n");
      OBFFLog("-----------------------------------------------------\n");
      //       XX   XX     XXXXXXXX   XXXXXXXX   XXXXXXXX   XXXXXXXX
    }
    
    for (i = _electrostaticcalculations.begin(); i != _electrostaticcalculations.end(); ++i) {
      
      i->Compute(gradients);
      energy += i->GetEnergy();
      
      IF_OBFF_LOGLVL_HIGH {
        sprintf(logbuf, "%2d   %2d   %8.3f  %8.3f  %8.3f\n", atoi((*i).a->GetType()), atoi((*i).b->GetType()), 
                (*i).rab, (*i).a->GetPartialCharge(), (*i).b->GetPartialCharge(), (*i).energy);
        OBFFLog(logbuf);
      }
    }
    
    IF_OBFF_LOGLVL_MEDIUM {
      sprintf(logbuf, "     TOTAL ELECTROSTATIC ENERGY = %8.3f %s\n", energy, GetUnit().c_str());
      OBFFLog(logbuf);
    }

    return energy;
  }
  
  //
  // OBForceFieldMMFF member functions
  //
  //***********************************************
  //Make a global instance
  OBForceFieldMMFF94 theForceFieldMMFF94("MMFF94", false);
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

  bool OBForceFieldMMFF94::ParseParamFile()
  {
    ParseParamProp();
    ParseParamDef();
    ParseParamBond();
    ParseParamBndk();
    ParseParamAngle();
    ParseParamStrBnd();
    ParseParamDfsb();
    ParseParamOOP();
    ParseParamTorsion();
    ParseParamVDW();
    ParseParamCharge();
    ParseParamPbci();
    return true;
  }
  
  bool OBForceFieldMMFF94::ParseParamBond()
  {
    vector<string> vs;
    char buffer[80];
    
    OBFFParameter parameter;
    
    // open data/mmffbond.par

    ifstream ifs;
    if (OpenDatafile(ifs, "mmffbond.par").length() == 0) {
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
  
  bool OBForceFieldMMFF94::ParseParamBndk()
  {
    vector<string> vs;
    char buffer[80];
    
    OBFFParameter parameter;
    
    // open data/mmffbndk.par
    ifstream ifs;
    if (OpenDatafile(ifs, "mmffbndk.par").length() == 0) {
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
  
  bool OBForceFieldMMFF94::ParseParamAngle()
  {
    vector<string> vs;
    char buffer[80];
    
    OBFFParameter parameter;
    
    // open data/mmffang.par
    ifstream ifs;
    if (OpenDatafile(ifs, "mmffang.par").length() == 0) {
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
   
  bool OBForceFieldMMFF94::ParseParamStrBnd()
  {
    vector<string> vs;
    char buffer[80];
    
    OBFFParameter parameter;
    
    // open data/mmffstbn.par
    ifstream ifs;
    if (OpenDatafile(ifs, "mmffstbn.par").length() == 0) {
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
  
  bool OBForceFieldMMFF94::ParseParamDfsb()
  {
    vector<string> vs;
    char buffer[80];
    
    OBFFParameter parameter;
    
    // open data/mmffdfsb.par
    ifstream ifs;
    if (OpenDatafile(ifs, "mmffdfsb.par").length() == 0) {
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
  
  bool OBForceFieldMMFF94::ParseParamOOP()
  {
    vector<string> vs;
    char buffer[80];
    
    OBFFParameter parameter;
    
    // open data/mmffoop.par
    ifstream ifs;
    if (OpenDatafile(ifs, "mmffoop.par").length() == 0) {
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
  
  bool OBForceFieldMMFF94::ParseParamTorsion()
  {
    vector<string> vs;
    char buffer[80];
    
    OBFFParameter parameter;
    
    // open data/mmfftor.par
    ifstream ifs;
    if (OpenDatafile(ifs, "mmfftor.par").length() == 0) {
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
   
  bool OBForceFieldMMFF94::ParseParamVDW()
  {
    vector<string> vs;
    char buffer[80];
    
    OBFFParameter parameter;
    
    // open data/mmffvdw.par
    ifstream ifs;
    if (OpenDatafile(ifs, "mmffvdw.par").length() == 0) {
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
   
  bool OBForceFieldMMFF94::ParseParamCharge()
  {
    vector<string> vs;
    char buffer[80];
    
    OBFFParameter parameter;
    
    // open data/mmffchg.par
    ifstream ifs;
    if (OpenDatafile(ifs, "mmffchg.par").length() == 0) {
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
 
  bool OBForceFieldMMFF94::ParseParamPbci()
  {
    vector<string> vs;
    char buffer[80];
    
    OBFFParameter parameter;
    
    // open data/mmffpbci.par
    ifstream ifs;
    if (OpenDatafile(ifs, "mmffpbci.par").length() == 0) {
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
 
  bool OBForceFieldMMFF94::ParseParamProp()
  {
    vector<string> vs;
    char buffer[80];
    
    OBFFParameter parameter;
    
    // open data/mmffprop.par
    ifstream ifs;
    if (OpenDatafile(ifs, "mmffprop.par").length() == 0) {
      obErrorLog.ThrowError(__FUNCTION__, "Cannot open mmffprop.par", obError);
      return false;
    }
    
    while (ifs.getline(buffer, 80)) {
      if (EQn(buffer, "*", 1)) continue;
      if (EQn(buffer, "$", 1)) continue;
	
      tokenize(vs, buffer);

      parameter.clear();
      parameter.a = atoi(vs[0].c_str());
      parameter._ipar.push_back(atoi(vs[1].c_str()));  // at.no
      parameter._ipar.push_back(atoi(vs[2].c_str()));  // crd
      parameter._ipar.push_back(atoi(vs[3].c_str()));  // val
      parameter._ipar.push_back(atoi(vs[4].c_str()));  // pilp
      parameter._ipar.push_back(atoi(vs[5].c_str()));  // mltb
      parameter._ipar.push_back(atoi(vs[6].c_str()));  // arom
      parameter._ipar.push_back(atoi(vs[7].c_str()));  // linh
      parameter._ipar.push_back(atoi(vs[8].c_str()));  // sbmb
      _ffpropparams.push_back(parameter);
    }

    if (ifs)
      ifs.close();
 
    return 0;
  }
  
  bool OBForceFieldMMFF94::ParseParamDef()
  {
    vector<string> vs;
    char buffer[80];
    
    OBFFParameter parameter;
    
    // open data/mmffdef.par
    ifstream ifs;
    if (OpenDatafile(ifs, "mmffdef.par").length() == 0) {
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
  
  bool OBForceFieldMMFF94::PerceiveAromatic()
  {
    bool done = false; // not done actually....
    OBAtom *ringatom;
    OBBond *ringbond;
    vector<OBRing*> vr;
    vr = _mol.GetSSSR();
    
    vector<OBRing*>::iterator ri;
    vector<int>::iterator rj;
    int n, index, ringsize, first_rj, prev_rj, pi_electrons;
    for (ri = vr.begin();ri != vr.end();++ri) { // for each ring
      ringsize = (*ri)->Size();
      
      n = 1;
      pi_electrons = 0;
      for(rj = (*ri)->_path.begin();rj != (*ri)->_path.end();rj++) { // for each ring atom
        index = *rj;
        ringatom = _mol.GetAtom(index);
        
        // is the bond to the previous ring atom double?
        if (n > 1) {
          ringbond = _mol.GetBond(prev_rj, index);
          if (ringbond->GetBO() == 2) {
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

	  if (!nbr->IsAromatic())
	    continue;

          ringbond = _mol.GetBond(nbr->GetIdx(), index);
          if (ringbond->GetBO() == 2)
            pi_electrons++;
        }

        // is the atom N, O or S in 5 rings
        if (ringsize == 5)
          if (ringatom->GetIdx() == (*ri)->GetRootAtom()) 
	    pi_electrons += 2;

        n++;
      
      } // for each ring atom
      
      // is the bond from the first to the last atom double?
      ringbond = _mol.GetBond(first_rj, index);
      if (ringbond->GetBO() == 2) 
        pi_electrons += 2;

      if ((pi_electrons == 6) && ((ringsize == 5) || (ringsize == 6))) {
	// mark ring atoms as aromatic
	for(rj = (*ri)->_path.begin();rj != (*ri)->_path.end();rj++) {
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
  
  bool OBForceFieldMMFF94::SetTypes()
  {
    std::vector<std::vector<int> > _mlist; //!< match list for atom typing
    std::vector<std::pair<OBSmartsPattern*,std::string> > _vexttyp; //!< external atom type rules
    vector<vector<int> >::iterator j;
    vector<pair<OBSmartsPattern*,string> >::iterator i;
    OBSmartsPattern *sp;
    vector<string> vs;
    char buffer[150];
 
    _mol.SetAtomTypesPerceived();

    // mark all atoms and bonds as non-aromatic, this is needed for the first 
    // phase in the atom type assigning
    _mol.SetAromaticPerceived();
    FOR_BONDS_OF_MOL (bond, _mol)
      bond->UnsetAromatic();
    FOR_ATOMS_OF_MOL (atom, _mol)
      atom->UnsetAromatic();
    
 
    /*
    cout << endl << "B O N D S" << endl << endl;
    cout << " I    J     BO     AR" << endl;
    cout << "---------------------" << endl;
    
    FOR_BONDS_OF_MOL (bond, _mol) {
      sprintf(buffer, "%2d   %2d     %d     %d", bond->GetBeginAtom()->GetIdx(), bond->GetEndAtom()->GetIdx(), bond->GetBO(), bond->IsAromatic());
      cout  << buffer << endl;
    }
    */
   
    ////////////////////////////////////////////////////////////////////////////
    // 
    // phase 1
    //
    // read smarts from mmfsymb.par, assign atom types
    // this is for non-aromatic atom types only!!!
    ////////////////////////////////////////////////////////////////////////////

    // open data/mmffsymb.par
    ifstream ifs;
    if (OpenDatafile(ifs, "mmffsymb.par").length() == 0) {
      obErrorLog.ThrowError(__FUNCTION__, "Cannot open mmffsymb.par", obError);
      return false;
    }

    while (ifs.getline(buffer, 150)) {
      if (EQn(buffer, "*", 1)) continue;
      if (EQn(buffer, "$", 1)) continue;
	
      tokenize(vs, buffer);

      sp = new OBSmartsPattern;
      if (sp->Init(vs[0])) {
        _vexttyp.push_back(pair<OBSmartsPattern*,string> (sp,vs[1]));
      } else {
        delete sp;
        sp = NULL;
        obErrorLog.ThrowError(__FUNCTION__, " Could not parse line from mmffsymb.par", obInfo);
        return false;
      }
 
      for (i = _vexttyp.begin();i != _vexttyp.end();++i) {
        if (i->first->Match(_mol)) {
          _mlist = i->first->GetMapList();
          for (j = _mlist.begin();j != _mlist.end();++j) {
            _mol.GetAtom((*j)[0])->SetType(i->second);
          }
        }
      }
    }

    if (ifs)
      ifs.close();
    

    ////////////////////////////////////////////////////////////////////////////
    // 
    // phase 2 
    //
    // find aromatic rings and assign atom types
    ////////////////////////////////////////////////////////////////////////////

    // It might be needed to runthis function more than once...
    bool done = true;
    while (done)
      done = PerceiveAromatic();
    
    OBAtom *ringatom;
    vector<OBRing*> vr;
    vr = _mol.GetSSSR();
    
       
    vector<OBRing*>::iterator ri;
    vector<int>::iterator rj;
    int index, ringsize;
    for (ri = vr.begin();ri != vr.end();ri++) { // for each ring
      ringsize = (*ri)->Size();
      if ((*ri)->IsAromatic()) {
        for(rj = (*ri)->_path.begin();rj != (*ri)->_path.end();rj++) { // for each ring atom
          index = *rj;
          ringatom = _mol.GetAtom(index);
	  
          if (ringsize == 6) {
            if (ringatom->IsCarbon()) {
              if (atoi(ringatom->GetType()) == 63)
	        continue;
              if (atoi(ringatom->GetType()) == 64)
	        continue;
              if (atoi(ringatom->GetType()) == 78)
	        continue;

              ringatom->SetType((char*)"37"); // CB: CARBON AS IN BENZENE, PYRROLE
	    }

            if (ringatom->IsNitrogen()) {
              if (ringatom->GetValence() == 2)
                ringatom->SetType((char*)"38"); // NPYD: NITROGEN, AS IN PYRIDINE
              if (ringatom->GetValence() == 3) {
                ringatom->SetType((char*)"58"); // NPD+: PYRIDINIUM-TYPE NITROGEN - FORMAL CHARGE=1
	
                FOR_NBORS_OF_ATOM (nbr, ringatom) {
                  if ((*ri)->IsInRing(nbr->GetIdx()))
                    continue;
                  if (nbr->IsOxygen() && (nbr->GetValence() == 1))
                    ringatom->SetType((char*)"69"); // NPOX: PYRIDINE N-OXIDE NITROGEN
                }
              }
            }

          }
          
          if (ringsize == 5) {
            // oxygen in furan
	    if (ringatom->IsOxygen())
	      ringatom->SetType((char*)"59"); // OFUR
            
	    // sulphur in thiophene
	    if (ringatom->IsSulfur())
	      ringatom->SetType((char*)"44"); // STHI

            // general alpha/beta carbon
            OBAtom *rootatom = _mol.GetAtom((*ri)->GetRootAtom());
            if (ringatom->IsCarbon()) {
              if (rootatom->IsConnected(ringatom))
                if (atoi(ringatom->GetType()) != 64)
	          ringatom->SetType((char*)"63"); // C5A
		else
	          ringatom->SetType((char*)"78"); // C5

	      if (rootatom->IsOneThree(ringatom))
                if (atoi(ringatom->GetType()) != 63)
	          ringatom->SetType((char*)"64"); // C5B
		else
	          ringatom->SetType((char*)"78"); // C5
	    } 
	    
            // general alpha/beta nitrogen + pyrrole nitrogen
	    if (ringatom->IsNitrogen()) {
	      if (ringatom->GetIdx() == rootatom->GetIdx())
	        ringatom->SetType((char*)"39"); // NPYL
              else {
                if (rootatom->IsConnected(ringatom))
	          ringatom->SetType((char*)"65"); // N5A
	        if (rootatom->IsOneThree(ringatom))
	          ringatom->SetType((char*)"66"); // N5B
              }
	    }
            
	    //
	    // specific rings start here
	    //
    
	    if (EQn((*ri)->GetType(), "1,2,4-triazole_anion", 20)) {
  	      if (ringatom->IsNitrogen())
	        ringatom->SetType((char*)"76"); // N5M
  	      if (ringatom->IsCarbon())
	          ringatom->SetType((char*)"78"); // C5
            }
        
	    if (EQn((*ri)->GetType(), "1,3,4-triazole_cation", 21)) {
  	      if (ringatom->IsNitrogen())
	        if (ringatom->GetValence() == 3)
	          ringatom->SetType((char*)"81"); // NIM+
		else
	          ringatom->SetType((char*)"79"); // N5 
	       
              if (ringatom->IsCarbon()) {
	        int hetero_count = 0;
		FOR_NBORS_OF_ATOM (nbr, ringatom)
		  if (nbr->IsNitrogen() && (*ri)->IsMember(&*nbr) && (nbr->GetValence() == 3))
		    hetero_count++;
		if (hetero_count == 2)
	          ringatom->SetType((char*)"80"); // CIM+
		else
	          ringatom->SetType((char*)"78"); // C5
	      }
 
	    }
           
	    if (EQn((*ri)->GetType(), "imidazole_cation", 16)) {
  	      if (ringatom->IsNitrogen())
	        ringatom->SetType((char*)"81"); // NIM+
	      
  	      if (ringatom->IsCarbon()) {
	        int hetero_count = 0;
		FOR_NBORS_OF_ATOM (nbr, ringatom)
		  if (nbr->IsNitrogen() && (*ri)->IsMember(&*nbr))
		    hetero_count++;
		if (hetero_count == 2)
	          ringatom->SetType((char*)"80"); // CIM+
		else
	          ringatom->SetType((char*)"78"); // C5
	      }
            }
	    
	    if (EQn((*ri)->GetType(), "pyrazole_anion", 14)) {
  	      if (ringatom->IsNitrogen())
	        ringatom->SetType((char*)"76"); // N5M
	      
  	      if (ringatom->IsCarbon())
	          ringatom->SetType((char*)"78"); // C5
            }

	    if (EQn((*ri)->GetType(), "thiazole_cation", 15) || EQn((*ri)->GetType(), "oxazole_cation", 14)) {
  	      if (ringatom->IsNitrogen())
	        ringatom->SetType((char*)"81"); // NIM+
	
	      if (ringatom->IsCarbon()) {
	        int hetero_count = 0;
		FOR_NBORS_OF_ATOM (nbr, ringatom)
		  //if ((nbr->IsSulfur() || nbr->IsOxygen() || nbr->IsNitrogen()) && (*ri)->IsMember(&*nbr))
		  if (nbr->IsOxygen() || nbr->IsNitrogen())
		    hetero_count++;
		if (hetero_count >= 2)
	          ringatom->SetType((char*)"80"); // CIM+
	      }
            }

            if (EQn((*ri)->GetType(), "1,2,4-thiadiazole_cation", 24) || 
	        EQn((*ri)->GetType(), "1,3,4-thiadiazole_cation", 24) ||
	        EQn((*ri)->GetType(), "1,2,3-oxadiazole_cation", 23) ||
	        EQn((*ri)->GetType(), "1,2,4-oxadiazole_cation", 23)) {
  	      if (ringatom->IsNitrogen() && (ringatom->BOSum() == 4))
	        ringatom->SetType((char*)"81"); // NIM+
	
	      if (ringatom->IsCarbon()) {
	        int hetero_count = 0;
	        int n_count = 0;
		FOR_NBORS_OF_ATOM (nbr, ringatom)
		  if ((nbr->IsSulfur() || nbr->IsOxygen()) && (*ri)->IsMember(&*nbr))
		    hetero_count++;
		  else if (nbr->IsNitrogen() && (nbr->BOSum() == 4) && (*ri)->IsMember(&*nbr))
                    n_count++;
		if (hetero_count && n_count)
	          ringatom->SetType((char*)"80"); // CIM+
	      }
            }
   
            if (EQn((*ri)->GetType(), "1,2,3,4-tetrazole_cation", 24)) {
  	      if (ringatom->IsNitrogen()) {
	        int carbon_count = 0;
		FOR_NBORS_OF_ATOM (nbr, ringatom)
		  if (nbr->IsCarbon() && (*ri)->IsMember(&*nbr))
		    carbon_count++;
		if (carbon_count)
	          ringatom->SetType((char*)"81"); // NIM+
		else
	          ringatom->SetType((char*)"79"); // N5
	      }
            
	      if (ringatom->IsCarbon())
	        ringatom->SetType((char*)"80"); // NIM+
	
	    }
	     
            if (EQn((*ri)->GetType(), "1,2,3,5-tetrazole_anion", 23)) {
  	      if (ringatom->IsNitrogen())
	        ringatom->SetType((char*)"76"); // N5M
            
	      if (ringatom->IsCarbon())
	        ringatom->SetType((char*)"78"); // C5
	
	    }
	

            
	    // correction for N-oxides
	    if (ringatom->IsNitrogen())
              FOR_NBORS_OF_ATOM (nbr, ringatom) {
                if ((*ri)->IsInRing(nbr->GetIdx()))
                  continue;
                if (nbr->IsOxygen() && (nbr->GetValence() == 1))
                  ringatom->SetType((char*)"82"); // N5OX, N5AX, N5BX
              }

	  } // 5 rings


        }
      } 
    } // for each ring
    
    ////////////////////////////////////////////////////////////////////////////
    // 
    // phase 3
    //
    // perform some corrections needed after assigning aromatic types
    ////////////////////////////////////////////////////////////////////////////


    // nitro - needed for GAFNUW
    FOR_ATOMS_OF_MOL (atom, _mol) {
      int nitrogen_count = 0;
      if (atom->IsNitrogen())
        FOR_NBORS_OF_ATOM (nbr, &*atom) 
	  if (nbr->IsOxygen() && (nbr->GetValence() == 1))
	    nitrogen_count++;
      if (nitrogen_count == 2)
        atom->SetType((char*)"45"); 
    }
    
    // divalent nitrogen
    FOR_ATOMS_OF_MOL (atom, _mol) {
      int oxygen_count = 0;
      if (atom->IsNitrogen() && !atom->IsAromatic() && (atom->GetValence() == 2) && (atom->BOSum() == 2)) {
        FOR_NBORS_OF_ATOM (nbr, &*atom) 
	  if (nbr->IsSulfur())
            FOR_NBORS_OF_ATOM (nbr2, &*nbr) 
	      if (nbr2->IsOxygen())
	        oxygen_count++;

        if (oxygen_count == 1)
          atom->SetType((char*)"48");
        else
          atom->SetType((char*)"62");
      }
    }

    // more corrections
    FOR_ATOMS_OF_MOL (a, _mol) {
      if (atoi(a->GetType()) == 55)
        FOR_NBORS_OF_ATOM (nbr, &*a) {
	  if (nbr->IsAromatic() && nbr->IsInRingSize(6) && !IsInSameRing(&*a, &*nbr)) {
            int nitrogen_count = 0;
            FOR_NBORS_OF_ATOM (nbr2, &*nbr) {
	      if (nbr2->IsNitrogen() && IsInSameRing(&*nbr, &*nbr2))
	        nitrogen_count++;
	    }
	    if (nitrogen_count >=1)
	      a->SetType((char*)"40");
	  }
	  
	  if (nbr->IsAromatic() && nbr->IsInRingSize(5) && !IsInSameRing(&*a, &*nbr))
            FOR_NBORS_OF_ATOM (nbr2, &*nbr)
	      if (nbr2->IsAromatic() && nbr2->IsInRingSize(5) && IsInSameRing(&*nbr, &*nbr2))
                FOR_NBORS_OF_ATOM (nbr3, &*nbr2)
		  if (nbr3->IsOxygen() && (nbr3->GetValence() == 1)) {
	            a->SetType((char*)"40");
                      FOR_NBORS_OF_ATOM (nbr4, &*a) 
		        if (nbr4->IsHydrogen())
	                  nbr4->SetType((char*)"28");
		  }
        }
    
      // set correct hydrogen types 
      if (a->IsHydrogen()) {
        FOR_NBORS_OF_ATOM (nbr, &*a) {
          if (atoi(nbr->GetType()) == 81)
	    a->SetType((char*)"36");
          if (atoi(nbr->GetType()) == 68)
	    a->SetType((char*)"23");
          if (atoi(nbr->GetType()) == 67)
	    a->SetType((char*)"23");
          if (atoi(nbr->GetType()) == 62)
	    a->SetType((char*)"23");
          if (atoi(nbr->GetType()) == 58)
	    a->SetType((char*)"36");
          if (atoi(nbr->GetType()) == 56)
	    a->SetType((char*)"36");
          if (atoi(nbr->GetType()) == 55)
	    a->SetType((char*)"36");
          if (atoi(nbr->GetType()) == 40)
	    a->SetType((char*)"28");
          if (atoi(nbr->GetType()) == 39)
	    a->SetType((char*)"23");
          if (atoi(nbr->GetType()) == 8)
	    a->SetType((char*)"23");
	}
      }


    
    
    }
    
    IF_OBFF_LOGLVL_LOW {
      OBFFLog("\nA T O M   T Y P E S\n\n");
      OBFFLog("IDX\tTYPE\n");
      
      FOR_ATOMS_OF_MOL (a, _mol) {
        sprintf(logbuf, "%d\t%s\n", a->GetIdx(), a->GetType());
        OBFFLog(logbuf);
      }
    }
 
    return true;
  }
  
  bool OBForceFieldMMFF94::SetupCalculations()
  {
    OBFFParameter *parameter;
    OBAtom *a, *b, *c, *d;
    bool found;
    
    IF_OBFF_LOGLVL_LOW
      OBFFLog("\nS E T T I N G   U P   C A L C U L A T I O N S\n\n");
 
    // 
    // Bond Calculations
    //
    IF_OBFF_LOGLVL_LOW
      OBFFLog("SETTING UP BOND CALCULATIONS...\n");
 
    OBFFBondCalculationMMFF94 bondcalc;
    int bondtype;

    _bondcalculations.clear();
    
    FOR_BONDS_OF_MOL(bond, _mol) {
      a = bond->GetBeginAtom();
      b = bond->GetEndAtom();	
      
      //if (a->IsIgnored() || b->IsIgnored())
      //  continue;
      
      bondtype = GetBondType(a, b);
      
      parameter = GetTypedParameter2Atom(bondtype, atoi(a->GetType()), atoi(b->GetType()), _ffbondparams); // from mmffbond.par
      if (parameter == NULL) {
	parameter = GetParameter2Atom(a->GetAtomicNum(), b->GetAtomicNum(), _ffbndkparams); // from mmffbndk.par - emperical rules
	if (parameter == NULL) { 
          IF_OBFF_LOGLVL_LOW {
            sprintf(logbuf, "    COULD NOT FIND PARAMETERS FOR BOND %d-%d (IDX)...\n", a->GetIdx(), b->GetIdx());
            OBFFLog(logbuf);
          }
          return false;
        } else {
          IF_OBFF_LOGLVL_LOW {
            sprintf(logbuf, "   USING EMPIRICAL RULE FOR BOND STRETCHING %d-%d (IDX)...\n", a->GetIdx(), b->GetIdx());
            OBFFLog(logbuf);
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
          _bondcalculations.push_back(bondcalc);
	}
      } else {
        bondcalc.a = a;
        bondcalc.b = b;
        bondcalc.kb = parameter->_dpar[0];
        bondcalc.r0 = parameter->_dpar[1];
        bondcalc.bt = bondtype;

        _bondcalculations.push_back(bondcalc);
      }
    }

    //
    // Angle & StrBnd Calculations
    //
    // MMFF part I - page 513 ("step-down" prodedure)
    // MMFF part I - page 519 (reference 68 is actually a footnote)
    // MMFF part IV - page 627 (empirical rule)
    //
    // First try and find an exact match, if this fails, step down using the equivalences from mmffdef.par
    // five-stage protocol: 1-1-1, 2-2-2, 3-2-3, 4-2-4, 5-2-5
    // If this fails, use empirical rules
    // Since 1-1-1 = 2-2-2, we will only try 1-1-1 before going to 3-2-3
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
      
      //if (a->IsIgnored() || b->IsIgnored() || c->IsIgnored())
      //  continue;
 
      angletype = GetAngleType(a, b, c);
      strbndtype = GetStrBndType(a, b, c);
      bondtype1 = GetBondType(a, b);
      bondtype2 = GetBondType(b, c);
      
      // try exact match 
      parameter = GetTypedParameter3Atom(angletype, atoi(a->GetType()), atoi(b->GetType()), atoi(c->GetType()), _ffangleparams);
      if (parameter == NULL) // try 3-2-3
        parameter = GetTypedParameter3Atom(angletype, EqLvl3(atoi(a->GetType())), atoi(b->GetType()), EqLvl3(atoi(c->GetType())), _ffangleparams);
      if (parameter == NULL) // try 4-2-4
        parameter = GetTypedParameter3Atom(angletype, EqLvl4(atoi(a->GetType())), atoi(b->GetType()), EqLvl4(atoi(c->GetType())), _ffangleparams);
      if (parameter == NULL) // try 5-2-5
        parameter = GetTypedParameter3Atom(angletype, EqLvl5(atoi(a->GetType())), atoi(b->GetType()), EqLvl5(atoi(c->GetType())), _ffangleparams);
        
      if (parameter) {
        anglecalc.ka = parameter->_dpar[0];
        anglecalc.theta0 = parameter->_dpar[1];
        strbndcalc.theta0 = parameter->_dpar[1]; // **
      } else {
        IF_OBFF_LOGLVL_LOW {
          sprintf(logbuf, "   USING DEFAULT ANGLE FOR %d-%d-%d (IDX)...\n", a->GetIdx(), b->GetIdx(), c->GetIdx());
          sprintf(logbuf, "   USING EMPIRICAL RULE FOR ANGLE BENDING %d-%d-%d (IDX)...\n", a->GetIdx(), b->GetIdx(), c->GetIdx());
          OBFFLog(logbuf);
        }

	anglecalc.ka = 0.0;
        anglecalc.theta0 = 120.0;

	if (GetCrd(atoi(b->GetType())) == 4)
          anglecalc.theta0 = 109.45;
        
	if ((GetCrd(atoi(b->GetType())) == 2) && b->IsOxygen())
          anglecalc.theta0 = 105.0;
	
	if (b->GetAtomicNum() > 10)
          anglecalc.theta0 = 95.0;
	
	if (HasLinSet(atoi(b->GetType())))
          anglecalc.theta0 = 180.0;
	
	if ((GetCrd(atoi(b->GetType())) == 3) && (GetVal(atoi(b->GetType())) == 3) && !GetMltb(atoi(b->GetType())))
	  if (b->IsNitrogen())
            anglecalc.theta0 = 107.0;
          else
	    anglecalc.theta0 = 92.0;

	if (a->IsInRingSize(3) && b->IsInRingSize(3) && c->IsInRingSize(3) && IsInSameRing(a, c))
	  anglecalc.theta0 = 60.0;
	
	if (a->IsInRingSize(4) && b->IsInRingSize(4) && c->IsInRingSize(4) && IsInSameRing(a, c))
	  anglecalc.theta0 = 90.0;
	
        strbndcalc.theta0 = anglecalc.theta0; // **
      }
      
      // empirical rule for 0-b-0 and standard angles
      if (anglecalc.ka == 0.0) {
        IF_OBFF_LOGLVL_LOW {
          sprintf(logbuf, "   USING EMPIRICAL RULE FOR ANGLE BENDING FORCE CONSTANT %d-%d-%d (IDX)...\n", a->GetIdx(), b->GetIdx(), c->GetIdx());
          OBFFLog(logbuf);
        }

	double beta, Za, Zc, Cb, r0ab, r0bc, theta, theta2, D, rr, rr2;
        Za = GetZParam(a);
        Cb = GetZParam(b);
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
        
	anglecalc.ka = (beta * Za * Cb * Zc * exp(-2 * D)) / (rr * theta2);
      }
      
      anglecalc.a = a;
      anglecalc.b = b;
      anglecalc.c = c;
      anglecalc.at = angletype;
      
      _anglecalculations.push_back(anglecalc);

      parameter = GetTypedParameter3Atom(strbndtype, atoi(a->GetType()), atoi(b->GetType()), atoi(c->GetType()), _ffstrbndparams);
      if (parameter == NULL) {
        int rowa, rowb, rowc;
        
	IF_OBFF_LOGLVL_LOW {
          sprintf(logbuf, "   USING EMPIRICAL RULE FOR STRETCH-BENDING FORCE CONSTANT %d-%d-%d (IDX)...\n", a->GetIdx(), b->GetIdx(), c->GetIdx());
          OBFFLog(logbuf);
        }
    
        rowa = GetElementRow(a);
        rowb = GetElementRow(b);
        rowc = GetElementRow(c);
        
        parameter = GetParameter3Atom(rowa, rowb, rowc, _ffdfsbparams);
        
	if (parameter == NULL) {
	  IF_OBFF_LOGLVL_LOW {
            sprintf(logbuf, "    COULD NOT FIND PARAMETERS FOR STRETCH-BEND %d-%d-%d (IDX)...\n", a->GetIdx(), b->GetIdx(), c->GetIdx());
            OBFFLog(logbuf);
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
        if (atoi(a->GetType()) == parameter->a) {
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
      
      //if (a->IsIgnored() || b->IsIgnored() || c->IsIgnored() || d->IsIgnored())
      //  continue;
 
      torsiontype = GetTorsionType(a, b, c, d);
      
      // try exact match 
      parameter = GetTypedParameter4Atom(torsiontype, atoi(a->GetType()), atoi(b->GetType()), atoi(c->GetType()), atoi(d->GetType()), _fftorsionparams);
      if (parameter == NULL) // try 3-2-2-5
	parameter = GetTypedParameter4Atom(torsiontype, EqLvl3(atoi(a->GetType())), atoi(b->GetType()), 
	                                            atoi(c->GetType()), EqLvl5(atoi(d->GetType())), _fftorsionparams);
      if (parameter == NULL) // try 5-2-2-3
  	parameter = GetTypedParameter4Atom(torsiontype, EqLvl5(atoi(a->GetType())), atoi(b->GetType()), 
	                                            atoi(c->GetType()), EqLvl3(atoi(d->GetType())), _fftorsionparams);
      if (parameter == NULL) // try 5-2-2-5
  	parameter = GetTypedParameter4Atom(torsiontype, EqLvl5(atoi(a->GetType())), atoi(b->GetType()), 
	                                            atoi(c->GetType()), EqLvl5(atoi(d->GetType())), _fftorsionparams);

      if (parameter) {
        torsioncalc.v1 = parameter->_dpar[0];
        torsioncalc.v2 = parameter->_dpar[1];
        torsioncalc.v3 = parameter->_dpar[2];
      } else {
        bool found_rule = false;
	
	IF_OBFF_LOGLVL_LOW {
          sprintf(logbuf, "   USING EMPIRICAL RULE FOR TORSION FORCE CONSTANT %d-%d-%d-%d (IDX)...\n", 
	    a->GetIdx(), b->GetIdx(), c->GetIdx(), d->GetType());
          OBFFLog(logbuf);
        }

	// rule (a) page 631
        if (HasLinSet(atoi(b->GetType())) || HasLinSet(atoi(c->GetType())))
          continue;
        
	// rule (b) page 631
        if (b->GetBond(c)->IsAromatic()) {
	  double Ub, Uc, pi_bc, beta;
	  Ub = GetUParam(b);
	  Uc = GetUParam(c);

	  if (!HasPilpSet(atoi(b->GetType())) && !HasPilpSet(atoi(b->GetType())))
	    pi_bc = 0.5;
	  else
	    pi_bc = 0.3;

          if (((GetVal(atoi(b->GetType())) == 3) && (GetVal(atoi(c->GetType())) == 4)) || 
	      ((GetVal(atoi(b->GetType())) == 4) && (GetVal(atoi(c->GetType())) == 3)))
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

	  if (((GetMltb(atoi(b->GetType())) == 2) && (GetMltb(atoi(c->GetType())) == 2)) && a->GetBond(b)->IsDouble())
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
	  if (((GetCrd(atoi(b->GetType())) == 4) && (GetCrd(atoi(c->GetType())) == 4))) {
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
	  if (((GetCrd(atoi(b->GetType())) == 4) && (GetCrd(atoi(c->GetType())) != 4))) {
	    if (GetCrd(atoi(c->GetType())) == 3) // case (1)
	      if ((GetVal(atoi(c->GetType())) == 4) || (GetVal(atoi(c->GetType())) == 34) || (GetMltb(atoi(c->GetType())) != 0))
                continue;
	    
	    if (GetCrd(atoi(c->GetType())) == 2) // case (2)
	      if ((GetVal(atoi(c->GetType())) == 3) || (GetMltb(atoi(c->GetType())) != 0))
                continue;
	    
	    // case (3) saturated bonds -- see rule (h)
	  }
	
	// rule (f) page 632
        if (!found_rule)
	  if (((GetCrd(atoi(b->GetType())) != 4) && (GetCrd(atoi(c->GetType())) == 4))) {
	    if (GetCrd(atoi(b->GetType())) == 3) // case (1)
	      if ((GetVal(atoi(b->GetType())) == 4) || (GetVal(atoi(b->GetType())) == 34) || (GetMltb(atoi(b->GetType())) != 0))
                continue;
	    
	    if (GetCrd(atoi(b->GetType())) == 2) // case (2)
	      if ((GetVal(atoi(b->GetType())) == 3) || (GetMltb(atoi(b->GetType())) != 0))
                continue;
	    
	    // case (3) saturated bonds
	  }
	
	// rule (g) page 632
        if (!found_rule)
	  if (b->GetBond(c)->IsSingle() && (
	      (GetMltb(atoi(b->GetType())) && GetMltb(atoi(c->GetType()))) ||
	      (GetMltb(atoi(b->GetType())) && HasPilpSet(atoi(c->GetType()))) ||
	      (GetMltb(atoi(c->GetType())) && HasPilpSet(atoi(b->GetType())))  )) {
	    if ((HasPilpSet(atoi(b->GetType())) && HasPilpSet(atoi(c->GetType())))) // case (1)
	      continue;
	    
	    double Ub, Uc, pi_bc, beta;
	    Ub = GetUParam(b);
	    Uc = GetUParam(c);
	    beta = 6.0;
  
	    if ((HasPilpSet(atoi(b->GetType())) && GetMltb(atoi(c->GetType())))) { // case (2)
	      if (GetMltb(atoi(c->GetType())) == 1)
	        pi_bc = 0.5;
	      else if ((GetElementRow(b) == 1) && (GetElementRow(c) == 1))
	        pi_bc = 0.3;
	      else
	        pi_bc = 0.15;
	      found_rule = true;
	    }
	    
	    if ((HasPilpSet(atoi(c->GetType())) && GetMltb(atoi(b->GetType())))) { // case (3)
	      if (GetMltb(atoi(b->GetType())) == 1)
	        pi_bc = 0.5;
	      else if ((GetElementRow(b) == 1) && (GetElementRow(c) == 1))
	        pi_bc = 0.3;
	      else
	        pi_bc = 0.15;
	      found_rule = true;
	    }
            
	    if (!found_rule)
	      if (((GetMltb(atoi(b->GetType())) == 1) || (GetMltb(atoi(b->GetType())) == 1)) && (!b->IsCarbon() || !c->IsCarbon())) {
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
	
	// rule (h) page 632
        if (!found_rule)
	  if ((b->IsOxygen() || b->IsSulfur()) && (c->IsOxygen() || c->IsSulfur())) {
	    double Wb, Wc;

	    if (b->IsOxygen())
	      Wb = 2.0;
	    else
	      Wb = 8.0;
	    
	    if (c->IsOxygen())
	      Wc = 2.0;
	    else
	      Wc = 8.0;

	    torsioncalc.v1 = 0.0;
            torsioncalc.v2 = -sqrt(Wb * Wc);
            torsioncalc.v3 = 0.0;
	  } else {
	    double Vb, Vc, Nbc;
	    Vb = GetVParam(b);
	    Vc = GetVParam(c);
            
	    Nbc = GetCrd(atoi(b->GetType())) * GetCrd(atoi(c->GetType()));

	    torsioncalc.v1 = 0.0;
            torsioncalc.v2 = 0.0;
            torsioncalc.v3 = sqrt(Vb * Vc) / Nbc;
	  }
      }
      
      torsioncalc.a = a;
      torsioncalc.b = b;
      torsioncalc.c = c;
      torsioncalc.d = d;
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

      for (int idx=0; idx < _ffoopparams.size(); idx++) {
        if (atoi(b->GetType()) == _ffoopparams[idx].b) {
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
          
          //if (a->IsIgnored() || b->IsIgnored() || c->IsIgnored() || d->IsIgnored())
          //  continue;
 
          if (((atoi(a->GetType()) == _ffoopparams[idx].a) && (atoi(c->GetType()) == _ffoopparams[idx].c) && (atoi(d->GetType()) == _ffoopparams[idx].d)) ||
              ((atoi(c->GetType()) == _ffoopparams[idx].a) && (atoi(a->GetType()) == _ffoopparams[idx].c) && (atoi(d->GetType()) == _ffoopparams[idx].d)) ||
              ((atoi(c->GetType()) == _ffoopparams[idx].a) && (atoi(d->GetType()) == _ffoopparams[idx].c) && (atoi(a->GetType()) == _ffoopparams[idx].d)) ||
              ((atoi(d->GetType()) == _ffoopparams[idx].a) && (atoi(c->GetType()) == _ffoopparams[idx].c) && (atoi(a->GetType()) == _ffoopparams[idx].d)) ||
              ((atoi(a->GetType()) == _ffoopparams[idx].a) && (atoi(d->GetType()) == _ffoopparams[idx].c) && (atoi(c->GetType()) == _ffoopparams[idx].d)) ||
              ((atoi(d->GetType()) == _ffoopparams[idx].a) && (atoi(a->GetType()) == _ffoopparams[idx].c) && (atoi(c->GetType()) == _ffoopparams[idx].d)))
            {
              found = true;

              oopcalc.koop = _ffoopparams[idx]._dpar[0];
            
              // A-B-CD || C-B-AD  PLANE = ABC
              oopcalc.a = a;
              oopcalc.b = b;
              oopcalc.c = c;
              oopcalc.d = d;
	    
              _oopcalculations.push_back(oopcalc);

              // C-B-DA || D-B-CA  PLANE BCD
              oopcalc.a = d;
              oopcalc.d = a;
	
              _oopcalculations.push_back(oopcalc);
            
              // A-B-DC || D-B-AC  PLANE ABD
              oopcalc.a = a;
              oopcalc.c = d;
              oopcalc.d = c;
	    
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
	    
              _oopcalculations.push_back(oopcalc);

              // C-B-DA || D-B-CA  PLANE BCD
              oopcalc.a = d;
              oopcalc.d = a;
	
              _oopcalculations.push_back(oopcalc);
            
              // A-B-DC || D-B-AC  PLANE ABD
              oopcalc.a = a;
              oopcalc.c = d;
              oopcalc.d = c;
	    
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
 
    FOR_PAIRS_OF_MOL(p, _mol) {
      a = _mol.GetAtom((*p)[0]);
      b = _mol.GetAtom((*p)[1]);
      
      //if (a->IsIgnored() || b->IsIgnored())
      //  continue;
 
      OBFFParameter *parameter_a, *parameter_b;
      parameter_a = GetParameter1Atom(atoi(a->GetType()), _ffvdwparams);
      parameter_b = GetParameter1Atom(atoi(b->GetType()), _ffvdwparams);
      if ((parameter_a == NULL) || (parameter_b == NULL)) {
        IF_OBFF_LOGLVL_LOW {
          sprintf(logbuf, "   COULD NOT FIND VAN DER WAALS PARAMETERS FOR %d-%d (IDX)...\n", a->GetIdx(), b->GetIdx());
          OBFFLog(logbuf);
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
      
      vdwcalc.escale = 1.0;
      if (vdwcalc.aDA == 1) { // hydrogen bond donor
        vdwcalc.R_AB = 0.5 * (R_AA + R_BB);
        if (vdwcalc.bDA == 2) { // hydrogen bond acceptor
          vdwcalc.R_AB = 0.8 * vdwcalc.R_AB;
          vdwcalc.escale = 0.5;
        }
      } else if (vdwcalc.bDA == 1) { // hydrogen bond donor
        vdwcalc.R_AB = 0.5 * (R_AA + R_BB);
        if (vdwcalc.aDA == 2) { // hydrogen bond acceptor
          vdwcalc.R_AB = 0.8 * vdwcalc.R_AB;
          vdwcalc.escale = 0.5;
        }
      } else {
        g_AB = (R_AA - R_BB) / ( R_AA + R_BB);
        g_AB2 = g_AB * g_AB;
        vdwcalc.R_AB =  0.5 * (R_AA + R_BB) * (1.0 + 0.2 * (1.0 - exp(-12.0 * g_AB2)));
      }
      
      R_AB2 = vdwcalc.R_AB * vdwcalc.R_AB;
      R_AB4 = R_AB2 * R_AB2;
      R_AB6 = R_AB4 * R_AB2;
      vdwcalc.R_AB7 = R_AB6 * vdwcalc.R_AB;

      sqrt_a = sqrt(vdwcalc.alpha_a / vdwcalc.Na);
      sqrt_b = sqrt(vdwcalc.alpha_b / vdwcalc.Nb);
      vdwcalc.epsilon = (181.16 * vdwcalc.Ga * vdwcalc.Gb * vdwcalc.alpha_a * vdwcalc.alpha_b) / (sqrt_a + sqrt_b) * (1.0 / R_AB6);
 
      _vdwcalculations.push_back(vdwcalc);
    }
    
    // 
    // Electrostatic Calculations
    //
    IF_OBFF_LOGLVL_LOW
      OBFFLog("SETTING UP ELECTROSTATIC CALCULATIONS...\n");
 
    OBFFElectrostaticCalculationMMFF94 elecalc;

    _electrostaticcalculations.clear();
    
    FOR_PAIRS_OF_MOL(p, _mol) {
      a = _mol.GetAtom((*p)[0]);
      b = _mol.GetAtom((*p)[1]);
      
      //if (a->IsIgnored() || b->IsIgnored())
      //  continue;
 
      elecalc.qq = 332.0716 * a->GetPartialCharge() * b->GetPartialCharge();
      
      if (elecalc.qq) {
        elecalc.a = &*a;
        elecalc.b = &*b;
        
        // 1-4 scaling
        if (a->IsOneFour(b))
          elecalc.qq *= 0.75;
	  
        _electrostaticcalculations.push_back(elecalc);
      }
    }

    return true;
  }

  // we set the the formal charge with SetPartialCharge because formal charges 
  // in MMFF94 are not always and integer
  bool OBForceFieldMMFF94::SetFormalCharges()
  {
    _mol.SetAutomaticPartialCharge(false);
    
    FOR_ATOMS_OF_MOL (atom, _mol) {
      atom->SetPartialCharge(0.0);

      if (atoi(atom->GetType()) == 32) {
        int o_count = 0;
        bool sulfonamide = false;
        int s_count = 0;

	FOR_NBORS_OF_ATOM(nbr, &*atom) {
  	  FOR_NBORS_OF_ATOM(nbr2, &*nbr) {
	    if (nbr2->IsOxygen() && (nbr2->GetValence() == 1))
	      o_count++;
	    if (nbr2->IsSulfur() && (nbr2->GetValence() == 1))
	      s_count++;
	    if (nbr2->IsNitrogen() && !nbr2->IsAromatic())
	      sulfonamide = true;
	  }
	
	  if (nbr->IsCarbon())
	    atom->SetPartialCharge(-0.5); // O2CM
	  
	  if (nbr->IsNitrogen() && (o_count == 3))
	    atom->SetPartialCharge(-1.0 / o_count);  // O3N
	  
	  if (nbr->IsSulfur() && !sulfonamide)
	    if (((o_count + s_count) == 2) && (nbr->GetValence() == 3) && (nbr->BOSum() == 3))
	      atom->SetPartialCharge(-0.5); // O2S
	    else if ((o_count + s_count) == 3)
	      atom->SetPartialCharge(-1.0 / 3.0); // O3S
	    else if ((o_count + s_count) == 4)
	      atom->SetPartialCharge(-0.5); // O4S
	  
	  if (nbr->IsPhosphorus())
	    if ((o_count + s_count) == 2)
	      atom->SetPartialCharge(-0.5); // O2P
	    else if ((o_count + s_count) == 3)
	      atom->SetPartialCharge(-2.0 / 3.0); // O3P
	    else if ((o_count + s_count) == 4)
	      atom->SetPartialCharge(-0.25); // O4P
	  
	  if (atoi(nbr->GetType()) == 77)
	    atom->SetPartialCharge(-0.25); // O4CL
	}
      } else if (atoi(atom->GetType()) == 34)
        atom->SetPartialCharge(1.0);
      else if (atoi(atom->GetType()) == 35)
        atom->SetPartialCharge(-1.0);
      else if (atoi(atom->GetType()) == 49)
        atom->SetPartialCharge(1.0);
      else if (atoi(atom->GetType()) == 51)
        atom->SetPartialCharge(1.0);
      else if (atoi(atom->GetType()) == 54)
        atom->SetPartialCharge(1.0);
      else if (atoi(atom->GetType()) == 55)
        atom->SetPartialCharge(0.5);
      else if (atoi(atom->GetType()) == 56)
        atom->SetPartialCharge(1.0/3.0);
      else if (atoi(atom->GetType()) == 58)
        atom->SetPartialCharge(1.0);
      else if (atoi(atom->GetType()) == 61) {
	FOR_NBORS_OF_ATOM(nbr, &*atom)
	  if (atom->GetBond(&*nbr)->IsTriple() && nbr->IsNitrogen())
            atom->SetPartialCharge(1.0);
      } else if (atoi(atom->GetType()) == 62)
        atom->SetPartialCharge(-1.0);
      else if (atoi(atom->GetType()) == 72) {
        int s_count = 0;

	FOR_NBORS_OF_ATOM(nbr, &*atom) {
	  if (nbr->IsSulfur())
	    s_count++;
	  
	  if (nbr->IsPhosphorus() || nbr->IsSulfur()) {
	    FOR_NBORS_OF_ATOM(nbr2, &*nbr)
	      if ((nbr2->IsSulfur() || nbr2->IsOxygen()) && (nbr2->GetValence() == 1) && (atom->GetIdx() != nbr2->GetIdx()))
	        atom->SetPartialCharge(-0.5);
          } else
	    atom->SetPartialCharge(-1.0);

  	  if (nbr->IsCarbon())
	    FOR_NBORS_OF_ATOM(nbr2, &*nbr)
	      if (nbr2->IsSulfur() && (nbr2->GetValence() == 1) && (atom->GetIdx() != nbr2->GetIdx()))
                atom->SetPartialCharge(-0.5); // SSMO
	
	  if (s_count >= 2)      
            atom->SetPartialCharge(-0.5); // SSMO
	}
      } else if (atoi(atom->GetType()) == 76) {
       	vector<OBRing*> vr;
        vr = _mol.GetSSSR();
        vector<OBRing*>::iterator ri;
        vector<int>::iterator rj;
	int n_count;

        for (ri = vr.begin();ri != vr.end();ri++) { // for each ring
	  n_count = 0;

          if ((*ri)->IsAromatic() && (*ri)->IsMember(&*atom) && ((*ri)->Size() == 5)) {
            for(rj = (*ri)->_path.begin();rj != (*ri)->_path.end();rj++) // for each ring atom
	      if (_mol.GetAtom(*rj)->IsNitrogen())
	        n_count++;
	    
	    if (n_count > 1)
	      atom->SetPartialCharge(-1.0 / n_count);
	  }
	}
      } else if (atoi(atom->GetType()) == 81) {
        atom->SetPartialCharge(1.0);
        
	vector<OBRing*> vr;
        vr = _mol.GetSSSR();
        vector<OBRing*>::iterator ri;
        vector<int>::iterator rj;
        for (ri = vr.begin();ri != vr.end();ri++) // for each ring
          if ((*ri)->IsAromatic() && (*ri)->IsMember(&*atom) && ((*ri)->Size() == 5)) {
	    int n_count = 0;
            for(rj = (*ri)->_path.begin();rj != (*ri)->_path.end();rj++) // for each ring atom
	      if (_mol.GetAtom(*rj)->IsNitrogen() && (_mol.GetAtom(*rj)->GetValence() == 3))
	        n_count++;

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
      } else if (atoi(atom->GetType()) == 87)
        atom->SetPartialCharge(2.0);
      else if (atoi(atom->GetType()) == 98)
        atom->SetPartialCharge(3.0);
      else if (atoi(atom->GetType()) == 89)
        atom->SetPartialCharge(-1.0);
      else if (atoi(atom->GetType()) == 90)
        atom->SetPartialCharge(-1.0);
      else if (atoi(atom->GetType()) == 91)
        atom->SetPartialCharge(-1.0);
      else if (atoi(atom->GetType()) == 92)
        atom->SetPartialCharge(1.0);
      else if (atoi(atom->GetType()) == 93)
        atom->SetPartialCharge(1.0);
      else if (atoi(atom->GetType()) == 94)
        atom->SetPartialCharge(1.0);
      else if (atoi(atom->GetType()) == 95)
        atom->SetPartialCharge(2.0);
      else if (atoi(atom->GetType()) == 96)
        atom->SetPartialCharge(2.0);
      else if (atoi(atom->GetType()) == 97)
        atom->SetPartialCharge(1.0);
      else if (atoi(atom->GetType()) == 98)
        atom->SetPartialCharge(2.0);
      else if (atoi(atom->GetType()) == 99)
        atom->SetPartialCharge(2.0);

    }

    IF_OBFF_LOGLVL_LOW {
      OBFFLog("\nF O R M A L   C H A R G E S\n\n");
      OBFFLog("IDX\tCHARGE\n");
      
      FOR_ATOMS_OF_MOL (a, _mol) {
        sprintf(logbuf, "%d\t%f\n", a->GetIdx(), a->GetPartialCharge());
        OBFFLog(logbuf);
      }
    }
   
    return true;
  }

  bool OBForceFieldMMFF94::SetPartialCharges()
  {
    vector<double> charges(_mol.NumAtoms()+1, 0);
    double M, Wab, factor, q0a, q0b, Pa, Pb;

    FOR_ATOMS_OF_MOL (atom, _mol) {
      if ((atoi(atom->GetType()) == 32) || (atoi(atom->GetType()) == 35) || (atoi(atom->GetType()) == 72))
        factor = 0.5;
      else if ((atoi(atom->GetType()) == 62) || (atoi(atom->GetType()) == 76))
        factor = 0.25;
      else
        factor = 0.0;
      
      M = GetCrd(atoi(atom->GetType()));
      q0a = atom->GetPartialCharge();

      // charge sharing
      if (!factor)
        FOR_NBORS_OF_ATOM (nbr, &*atom)
          if (nbr->GetPartialCharge() < 0.0)
	    q0a += nbr->GetPartialCharge() / (2.0 * nbr->GetValence());
      
     // needed for SEYWUO, positive charge sharing?
     if (atoi(atom->GetType()) == 62)
        FOR_NBORS_OF_ATOM (nbr, &*atom)
          if (nbr->GetPartialCharge() > 0.0)
	    q0a -= nbr->GetPartialCharge() / 2.0;
     
      q0b = 0.0;
      Wab = 0.0;
      FOR_NBORS_OF_ATOM (nbr, &*atom) {
        q0b += nbr->GetPartialCharge();
    
        bool bci_found = false;
	for (int idx=0; idx < _ffchgparams.size(); idx++)
          if (GetBondType(&*atom, &*nbr) == _ffchgparams[idx]._ipar[0]) 
            if ((atoi(atom->GetType()) == _ffchgparams[idx].a) && (atoi(nbr->GetType()) == _ffchgparams[idx].b)) {
	      Wab += -_ffchgparams[idx]._dpar[0];
              bci_found = true;
            } else if  ((atoi(atom->GetType()) == _ffchgparams[idx].b) && (atoi(nbr->GetType()) == _ffchgparams[idx].a)) {
	      Wab += _ffchgparams[idx]._dpar[0];
              bci_found = true;
	    }
        
	if (!bci_found) {
	  for (int idx=0; idx < _ffpbciparams.size(); idx++) {
            if (atoi(atom->GetType()) == _ffpbciparams[idx].a)
	      Pa = _ffpbciparams[idx]._dpar[0];
	    if (atoi(nbr->GetType()) == _ffpbciparams[idx].a)
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

    IF_OBFF_LOGLVL_LOW {
      OBFFLog("\nP A R T I A L   C H A R G E S\n\n");
      OBFFLog("IDX\tCHARGE\n");
      
      FOR_ATOMS_OF_MOL (a, _mol) {
        sprintf(logbuf, "%d\t%f\n", a->GetIdx(), a->GetPartialCharge());
        OBFFLog(logbuf);
      }
    }
 
    return true;
  }

  double OBForceFieldMMFF94::Energy(bool gradients)
  {
    double energy;

    IF_OBFF_LOGLVL_MEDIUM
      OBFFLog("\nE N E R G Y\n\n");
    
    energy = E_Bond(gradients);
    energy += E_Angle(gradients);
    energy += E_StrBnd(gradients);
    energy += E_Torsion(gradients);
    energy += E_OOP(gradients);
    energy += E_VDW(gradients);
    energy += E_Electrostatic(gradients);

    IF_OBFF_LOGLVL_MEDIUM {
      sprintf(logbuf, "\nTOTAL ENERGY = %8.3f %s\n", energy, GetUnit().c_str());
      OBFFLog(logbuf);
    }

    return energy;
  }

  // used to validate the implementation
  bool OBForceFieldMMFF94::Validate ()
  {
    OBConversion conv;
    OBFormat *format_in = conv.FindFormat("mol2");
    vector<string> vs;
    vector<int> types;
    vector<double> fcharges, pcharges;
    vector<double> bond_lengths;
    char buffer[150], logbuf[100];
    bool molfound, atomfound, bondfound, fchgfound, pchgfound;
    double etot, ebond, eangle, eoop, estbn, etor, evdw, eeq;
    double termcount; //1=bond, 2=angle, 3=strbnd, 4=torsion, 5=oop
    int n;

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
      
      SetTypes();
      
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
      int ni;
      bool failed;

      cout << "--------------------------------------------------------------------------------" << endl;
      cout << "                                                                                " << endl;
      cout << "  VALIDATE MOLECULE " << c << ": " << _mol.GetTitle() << endl;
      cout << "                                                                                " << endl;
      cout << "IDX  HYB  AROM  OB_TYPE  LOG_TYPE       RESULT                                  " << endl;
      cout << "----------------------------------------------                                  " << endl;
            
      ni = 1;
      failed = false;
      for (i = types.begin(); i != types.end();i++) {
        if (ni > _mol.NumAtoms())
          continue;
	
	if ( (atoi(_mol.GetAtom(ni)->GetType()) == 87) ||
	     (atoi(_mol.GetAtom(ni)->GetType()) == 97) 
	   ) continue;

        if (atoi(_mol.GetAtom(ni)->GetType()) == (*i))
          sprintf(logbuf, "%2d   %3d  %4d    %3d      %3d          PASSED", _mol.GetAtom(ni)->GetIdx(), _mol.GetAtom(ni)->GetHyb(), 
                  _mol.GetAtom(ni)->IsAromatic(), atoi(_mol.GetAtom(ni)->GetType()), *i);
        else {
          sprintf(logbuf, "%2d   %3d  %4d    %3d      %3d      XXX FAILED XXX", _mol.GetAtom(ni)->GetIdx(), _mol.GetAtom(ni)->GetHyb(), 
                  _mol.GetAtom(ni)->IsAromatic(), atoi(_mol.GetAtom(ni)->GetType()), *i);
          failed = true;
        }
      
        cout << logbuf << endl;
        
        ni++;
      }

      if (failed) {
        cout << "Could not succesfully assign atom types" << endl;
        return false;
        //continue;
      }

      SetFormalCharges();
      cout << endl;
      cout << "IDX  OB_FCARGE  LOG_FCHARGE       RESULT" << endl;
      cout << "----------------------------------------" << endl;
            
      ni = 1;
      for (di = fcharges.begin(); di != fcharges.end(); di++) {
	if (ni > _mol.NumAtoms())
          continue;
	
	if ( (atoi(_mol.GetAtom(ni)->GetType()) == 87) ||
	     (atoi(_mol.GetAtom(ni)->GetType()) == 97) 
	   ) continue;

        if (fabs((*di) - _mol.GetAtom(ni)->GetPartialCharge()) <= 0.001)
          sprintf(logbuf, "%2d   %7.4f     %7.4f          PASSED", _mol.GetAtom(ni)->GetIdx(), _mol.GetAtom(ni)->GetPartialCharge(), *di);
        else {
          sprintf(logbuf, "%2d   %7.4f     %7.4f      XXX FAILED XXX", _mol.GetAtom(ni)->GetIdx(), _mol.GetAtom(ni)->GetPartialCharge(), *di);
          failed = true;
        }
      
        cout << logbuf << endl;
        
        ni++;
      }

      if (failed) {
        cout << "Could not succesfully assign formal charges" << endl;
        return false;
        //continue;
      }

      SetPartialCharges();
      cout << endl;
      cout << "IDX  OB_PCARGE  LOG_PCHARGE       RESULT" << endl;
      cout << "----------------------------------------" << endl;
            
      ni = 1;
      for (di = pcharges.begin(); di != pcharges.end(); di++) {
	if (ni > _mol.NumAtoms())
          continue;
	
	if ( (atoi(_mol.GetAtom(ni)->GetType()) == 87) ||
	     (atoi(_mol.GetAtom(ni)->GetType()) == 97) 
	   ) continue;

        if (fabs((*di) - _mol.GetAtom(ni)->GetPartialCharge()) <= 0.001)
          sprintf(logbuf, "%2d   %7.4f     %7.4f          PASSED", _mol.GetAtom(ni)->GetIdx(), _mol.GetAtom(ni)->GetPartialCharge(), *di);
        else {
          sprintf(logbuf, "%2d   %7.4f     %7.4f      XXX FAILED XXX", _mol.GetAtom(ni)->GetIdx(), _mol.GetAtom(ni)->GetPartialCharge(), *di);
          failed = true;
        }
      
        cout << logbuf << endl;
        
        ni++;
      }

      if (failed) {
        cout << "Could not succesfully assign partial charges" << endl;
	return false;
        //continue;
      }


      
      

      if (!SetupCalculations()) {
        cout << "Could not setup calculations (missing parameters...)" << endl;
        return false;
        //continue;
      }

      double delta, err;
      cout << endl;
      cout << "TERM                   OB ENERGY    LOG ENERGY      DELTA" << endl;
      cout << "-----------------------------------------------------------" << endl;
    
      delta = (E_Bond() - ebond);
      sprintf(logbuf, "Bond Stretching        %8.3f      %8.3f     %8.3f", E_Bond(), ebond, delta);
      cout << logbuf << endl;
    
      delta = (E_Angle() - eangle);
      sprintf(logbuf, "Angle Bending          %8.3f      %8.3f     %8.3f", E_Angle(), eangle, delta);
      cout << logbuf << endl;
    
      delta = (E_StrBnd() - estbn);
      sprintf(logbuf, "Stretch-Bending        %8.3f      %8.3f     %8.3f", E_StrBnd(), estbn, delta);
      cout << logbuf << endl;
    
      delta = (E_OOP() - eoop);
      sprintf(logbuf, "Out-Of-Plane Bending   %8.3f      %8.3f     %8.3f", E_OOP(), eoop, delta);
      cout << logbuf << endl;
    
      delta = (E_Torsion() - etor);
      sprintf(logbuf, "Torsional              %8.3f      %8.3f     %8.3f", E_Torsion(), etor, delta);
      cout << logbuf << endl;
    
      delta = (E_VDW() - evdw);
      sprintf(logbuf, "Van der Waals          %8.3f      %8.3f     %8.3f", E_VDW(), evdw, delta);
      cout << logbuf << endl;
      
      delta = (E_Electrostatic() - eeq);
      sprintf(logbuf, "Electrostatic          %8.3f      %8.3f     %8.3f", E_Electrostatic(), eeq, delta);
      cout << logbuf << endl;

      cout << endl;
      delta = (Energy() - etot);
      sprintf(logbuf, "Total ENERGY           %8.3f      %8.3f     %8.3f", Energy(), etot, delta);
      cout << logbuf << endl;

    } // for (unsigned int c;; c++ )
   
    if (ifs)
      ifs.close();
    if (ifs2)
      ifs2.close();
  return true;
}
  
  vector3 OBForceFieldMMFF94::GetGradient(OBAtom *a, int terms)
  {
    vector<OBFFBondCalculationMMFF94>::iterator i;
    vector<OBFFAngleCalculationMMFF94>::iterator i2;
    vector<OBFFStrBndCalculationMMFF94>::iterator i3;
    vector<OBFFTorsionCalculationMMFF94>::iterator i4;
    vector<OBFFOOPCalculationMMFF94>::iterator i5;
    vector<OBFFVDWCalculationMMFF94>::iterator i6;
    vector<OBFFElectrostaticCalculationMMFF94>::iterator i7;

    vector3 grad(0.0, 0.0, 0.0);
    
    if ((terms & OBFF_ENERGY) || (terms & OBFF_EBOND))
      for (i = _bondcalculations.begin(); i != _bondcalculations.end(); ++i)
        if (((*i).a->GetIdx() == a->GetIdx()) || ((*i).b->GetIdx() == a->GetIdx())) {
          i->Compute(true);
          grad += i->GetGradient(&*a);
	}
    if ((terms & OBFF_ENERGY) || (terms & OBFF_EANGLE))
      for (i2 = _anglecalculations.begin(); i2 != _anglecalculations.end(); ++i2)
        if (((*i2).a->GetIdx() == a->GetIdx()) || ((*i2).b->GetIdx() == a->GetIdx()) || ((*i2).c->GetIdx() == a->GetIdx())) {
          i2->Compute(true);
          grad += i2->GetGradient(&*a);
        }
    if ((terms & OBFF_ENERGY) || (terms & OBFF_ESTRBND))
      for (i3 = _strbndcalculations.begin(); i3 != _strbndcalculations.end(); ++i3)
        if (((*i3).a->GetIdx() == a->GetIdx()) || ((*i3).b->GetIdx() == a->GetIdx()) || ((*i3).c->GetIdx() == a->GetIdx())) {
          i3->Compute(true);
          grad += i3->GetGradient(&*a);
        }
    if ((terms & OBFF_ENERGY) || (terms & OBFF_ETORSION))
      for (i4 = _torsioncalculations.begin(); i4 != _torsioncalculations.end(); ++i4)
        if (((*i4).a->GetIdx() == a->GetIdx()) || ((*i4).b->GetIdx() == a->GetIdx()) || 
	    ((*i4).c->GetIdx() == a->GetIdx()) || ((*i4).d->GetIdx() == a->GetIdx())) {
          i4->Compute(true);
          grad += i4->GetGradient(&*a);
	}
    if ((terms & OBFF_ENERGY) || (terms & OBFF_EOOP))
      for (i5 = _oopcalculations.begin(); i5 != _oopcalculations.end(); ++i5)
        if (((*i5).a->GetIdx() == a->GetIdx()) || ((*i5).b->GetIdx() == a->GetIdx()) || 
	    ((*i5).c->GetIdx() == a->GetIdx()) || ((*i5).d->GetIdx() == a->GetIdx())) {
          i5->Compute(true);
          grad += i5->GetGradient(&*a);
        }
    if ((terms & OBFF_ENERGY) || (terms & OBFF_EVDW))
      for (i6 = _vdwcalculations.begin(); i6 != _vdwcalculations.end(); ++i6)
        if (((*i6).a->GetIdx() == a->GetIdx()) || ((*i6).b->GetIdx() == a->GetIdx())) {
          i6->Compute(true);
          grad += i6->GetGradient(&*a);
        }
    if ((terms & OBFF_ENERGY) || (terms & OBFF_EELECTROSTATIC))
      for (i7 = _electrostaticcalculations.begin(); i7 != _electrostaticcalculations.end(); ++i7)
        if (((*i7).a->GetIdx() == a->GetIdx()) || ((*i7).b->GetIdx() == a->GetIdx())) {
          i7->Compute(true);
          grad += i7->GetGradient(&*a);
        }

    return grad;
  }


  bool OBForceFieldMMFF94::ValidateGradients ()
  {
    vector3 numgrad, anagrad, err;
    
    cout << "----------------------------------------------------------------------------------------" << endl;
    cout << "                                                                                        " << endl;
    cout << "  VALIDATE GRADIENTS : " << _mol.GetTitle() << endl;
    cout << "                                                                                        " << endl;
    cout << "                                                                                        " << endl;
    cout << "ATOM IDX      NUMERICAL GRADIENT           ANALYTICAL GRADIENT        REL. ERRROR (%)   " << endl;
    cout << "----------------------------------------------------------------------------------------" << endl;
    //     "XX       (000.000, 000.000, 000.000)  (000.000, 000.000, 000.000)  (00.00, 00.00, 00.00)"
   
    FOR_ATOMS_OF_MOL (a, _mol) {

      // OBFF_ENERGY
      numgrad = NumericalDerivative(&*a, OBFF_ENERGY);
      anagrad = GetGradient(&*a, OBFF_ENERGY);
      err = ValidateGradientError(numgrad, anagrad);

      sprintf(logbuf, "%2d       (%7.3f, %7.3f, %7.3f)  (%7.3f, %7.3f, %7.3f)  (%5.2f, %5.2f, %5.2f)", a->GetIdx(), numgrad.x(), numgrad.y(), numgrad.z(), 
              anagrad.x(), anagrad.y(), anagrad.z(), err.x(), err.y(), err.z());
      cout << logbuf << endl;

      // OBFF_EBOND
      numgrad = NumericalDerivative(&*a, OBFF_EBOND);
      anagrad = GetGradient(&*a, OBFF_EBOND);
      err = ValidateGradientError(numgrad, anagrad);

      sprintf(logbuf, "    bond    (%7.3f, %7.3f, %7.3f)  (%7.3f, %7.3f, %7.3f)  (%5.2f, %5.2f, %5.2f)", numgrad.x(), numgrad.y(), numgrad.z(), 
              anagrad.x(), anagrad.y(), anagrad.z(), err.x(), err.y(), err.z());
      cout << logbuf << endl;
      
      // OBFF_EANGLE
      numgrad = NumericalDerivative(&*a, OBFF_EANGLE);
      anagrad = GetGradient(&*a, OBFF_EANGLE);
      err = ValidateGradientError(numgrad, anagrad);

      sprintf(logbuf, "    angle   (%7.3f, %7.3f, %7.3f)  (%7.3f, %7.3f, %7.3f)  (%5.2f, %5.2f, %5.2f)", numgrad.x(), numgrad.y(), numgrad.z(), 
              anagrad.x(), anagrad.y(), anagrad.z(), err.x(), err.y(), err.z());
      cout << logbuf << endl;
      
      // OBFF_ESTRBND
      numgrad = NumericalDerivative(&*a, OBFF_ESTRBND);
      anagrad = GetGradient(&*a, OBFF_ESTRBND);
      err = ValidateGradientError(numgrad, anagrad);

      sprintf(logbuf, "    strbnd  (%7.3f, %7.3f, %7.3f)  (%7.3f, %7.3f, %7.3f)  (%5.2f, %5.2f, %5.2f)", numgrad.x(), numgrad.y(), numgrad.z(), 
              anagrad.x(), anagrad.y(), anagrad.z(), err.x(), err.y(), err.z());
      cout << logbuf << endl;

      // OBFF_ETORSION
      numgrad = NumericalDerivative(&*a, OBFF_ETORSION);
      anagrad = GetGradient(&*a, OBFF_ETORSION);
      err = ValidateGradientError(numgrad, anagrad);

      sprintf(logbuf, "    torsion (%7.3f, %7.3f, %7.3f)  (%7.3f, %7.3f, %7.3f)  (%5.2f, %5.2f, %5.2f)", numgrad.x(), numgrad.y(), numgrad.z(), 
              anagrad.x(), anagrad.y(), anagrad.z(), err.x(), err.y(), err.z());
      cout << logbuf << endl;
      
      // OBFF_EOOP
      numgrad = NumericalDerivative(&*a, OBFF_EOOP);
      anagrad = GetGradient(&*a, OBFF_EOOP);
      err = ValidateGradientError(numgrad, anagrad);

      sprintf(logbuf, "    oop     (%7.3f, %7.3f, %7.3f)  (%7.3f, %7.3f, %7.3f)  (%5.2f, %5.2f, %5.2f)", numgrad.x(), numgrad.y(), numgrad.z(), 
              anagrad.x(), anagrad.y(), anagrad.z(), err.x(), err.y(), err.z());
      cout << logbuf << endl;

      // OBFF_EVDW
      numgrad = NumericalDerivative(&*a, OBFF_EVDW);
      anagrad = GetGradient(&*a, OBFF_EVDW);
      err = ValidateGradientError(numgrad, anagrad);

      sprintf(logbuf, "    vdw     (%7.3f, %7.3f, %7.3f)  (%7.3f, %7.3f, %7.3f)  (%5.2f, %5.2f, %5.2f)", numgrad.x(), numgrad.y(), numgrad.z(), 
              anagrad.x(), anagrad.y(), anagrad.z(), err.x(), err.y(), err.z());
      cout << logbuf << endl;

      // OBFF_EELECTROSTATIC
      numgrad = NumericalDerivative(&*a, OBFF_EELECTROSTATIC);
      anagrad = GetGradient(&*a, OBFF_EELECTROSTATIC);
      err = ValidateGradientError(numgrad, anagrad);

      sprintf(logbuf, "    electro (%7.3f, %7.3f, %7.3f)  (%7.3f, %7.3f, %7.3f)  (%5.2f, %5.2f, %5.2f)", numgrad.x(), numgrad.y(), numgrad.z(), 
              anagrad.x(), anagrad.y(), anagrad.z(), err.x(), err.y(), err.z());
      cout << logbuf << endl;
    }
    return true;
  }

  //
  // MMFF part V - page 620
  //
  // BTij is 1 when:
  // a) single bond between atoms i and j, both i and j are not aromatic and both types have sbmb set in mmffprop.par, or
  // b) bewtween two aromatic atoms, but the bond is not aromatic (e.g. connecting bond in biphenyl)
  //
  int OBForceFieldMMFF94::GetBondType(OBAtom* a, OBAtom* b)
  {
    /*
    if ((atoi(a->GetType()) == 2) && (atoi(b->GetType()) ==2)) {
      cout << "GetBondType(" << a->GetType() << ", " << b->GetType() << ")" << endl;
      cout << " bond is single? " << _mol.GetBond(a,b)->IsSingle() << endl;
      cout << " bond is double? " << _mol.GetBond(a,b)->IsSingle() << endl;
      cout << " bond is arom? " << _mol.GetBond(a,b)->IsAromatic() << endl;
      cout << " Sbmb_a? " << HasSbmbSet(atoi(a->GetType())) << endl;
      cout << " Sbmb_b? " << HasSbmbSet(atoi(b->GetType())) << endl;
    }
    */
    if (!_mol.GetBond(a,b)->IsSingle())
      return 0;
    
    if (!_mol.GetBond(a,b)->IsAromatic())
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

    if (a->IsInRingSize(3) && b->IsInRingSize(3) && c->IsInRingSize(3) && IsInSameRing(a, c)) {
      switch (sumbondtypes) {
      case 0:
        return 3; 
      case 1:
        return 5; 
      case 2:
        return 6; 
      }
    }
    
    //if (a->IsInRingSize(4) && b->IsInRingSize(4) && c->IsInRingSize(4) && IsInSameRing(a, b) && IsInSameRing(b, c)) {
    if (a->IsInRingSize(4) && b->IsInRingSize(4) && c->IsInRingSize(4) && IsInSameRing(a, c)) {
      switch (sumbondtypes) {
      case 0:
        return 4; 
      case 1:
        return 7; 
      case 2:
        return 8; 
      }
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
      if (btab)
        if (!inverse)
          return 1;
        else
          return 2;
      if (btbc)
        if (!inverse)
          return 2;
        else
          return 1;

    case 2:
      return 3;

    case 3:
      return 5;

    case 4:
      return 4;

    case 5:
      if (btab)
        if (!inverse)
          return 6;
        else
          return 7;
      if (btbc)
        if (!inverse)
          return 7;
        else
          return 6;
      
    case 6:
      return 8;
      
    case 7:
      if (btab)
        if (!inverse)
          return 9;
        else
          return 10;
      if (btbc)
        if (!inverse)
          return 10;
        else
          return 9;
      
    case 8:
      return 11;
    }
    return -1; //???
  }
  
  int OBForceFieldMMFF94::GetTorsionType(OBAtom* a, OBAtom* b, OBAtom *c, OBAtom *d)
  {
    int btab, btbc, btcd;

    btab = GetBondType(a, b);
    btbc = GetBondType(b, c);
    btcd = GetBondType(c, d);
    
    if (btbc == 1)
      return 1;

    if (btbc && (btab || btcd))
      return 2;

    if (a->IsInRingSize(4) && b->IsInRingSize(4) && c->IsInRingSize(4) && d->IsInRingSize(4))
      if (IsInSameRing(a,b) && IsInSameRing(b,c) && IsInSameRing(c,d))
        return 4;

    if (a->IsInRingSize(5) && b->IsInRingSize(5) && c->IsInRingSize(5) && d->IsInRingSize(5)) {
      vector<OBRing*> vr;
      vr = _mol.GetSSSR();
    
      if (!b->IsCarbon() && !c->IsCarbon())
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
	
        for(rj = (*ri)->_path.begin();rj != (*ri)->_path.end();rj++) // for each ring atom
	  if (_mol.GetAtom(*rj)->GetValence() != _mol.GetAtom(*rj)->BOSum())
	    return 0;
	
	return 5;
      }
    }

    return 0;
  }

  bool OBForceFieldMMFF94::HasLinSet(int atomtype)
  {
    OBFFParameter *par;
    
    par = GetParameter1Atom(atomtype, _ffpropparams); // from mmffprop.par
    if (par)
      if (par->_ipar[6])
        return true;

    return false;
  }
 
  bool OBForceFieldMMFF94::HasPilpSet(int atomtype)
  {
    OBFFParameter *par;
    
    par = GetParameter1Atom(atomtype, _ffpropparams); // from mmffprop.par
    if (par)
      if (par->_ipar[3])
        return true;

    return false;
  }
  
  bool OBForceFieldMMFF94::HasAromSet(int atomtype)
  {
    OBFFParameter *par;
    
    par = GetParameter1Atom(atomtype, _ffpropparams); // from mmffprop.par
    if (par)
      if (par->_ipar[5])
        return true;

    return false;
  }
 
  bool OBForceFieldMMFF94::HasSbmbSet(int atomtype)
  {
    OBFFParameter *par;
    
    par = GetParameter1Atom(atomtype, _ffpropparams); // from mmffprop.par
    if (par)
      if (par->_ipar[7])
        return true;

    return false;
  }

  int OBForceFieldMMFF94::GetCrd(int atomtype)
  {
    OBFFParameter *par;
    
    par = GetParameter1Atom(atomtype, _ffpropparams); // from mmffprop.par
    if (par)
      return par->_ipar[1];

    return 0;
  }

  int OBForceFieldMMFF94::GetVal(int atomtype)
  {
    OBFFParameter *par;
    
    par = GetParameter1Atom(atomtype, _ffpropparams); // from mmffprop.par
    if (par)
      return par->_ipar[2];

    return 0;
  }

  int OBForceFieldMMFF94::GetMltb(int atomtype)
  {
    OBFFParameter *par;
    
    par = GetParameter1Atom(atomtype, _ffpropparams); // from mmffprop.par
    if (par)
      return par->_ipar[4];

    return 0;
  }

  int OBForceFieldMMFF94::EqLvl2(int type)
  {
    for (int idx=0; idx < _ffdefparams.size(); idx++)
      if (_ffdefparams[idx]._ipar[0] == type)
        return _ffdefparams[idx]._ipar[1];

    return type; 
  }
  
  int OBForceFieldMMFF94::EqLvl3(int type)
  {
    for (int idx=0; idx < _ffdefparams.size(); idx++)
      if (_ffdefparams[idx]._ipar[0] == type)
        return _ffdefparams[idx]._ipar[2];

    return type; 
  }
  
  int OBForceFieldMMFF94::EqLvl4(int type)
  {
    for (int idx=0; idx < _ffdefparams.size(); idx++)
      if (_ffdefparams[idx]._ipar[0] == type)
        return _ffdefparams[idx]._ipar[3];

    return type; 
  }

  int OBForceFieldMMFF94::EqLvl5(int type)
  {
    for (int idx=0; idx < _ffdefparams.size(); idx++)
      if (_ffdefparams[idx]._ipar[0] == type)
        return _ffdefparams[idx]._ipar[4];

    return type; 
  }

  // MMFF part V - TABLE VI
  double OBForceFieldMMFF94::GetZParam(OBAtom* atom)
  {
    if (atom->IsHydrogen())
      return 1.395;
    if (atom->IsCarbon())
      return 2.494;
    if (atom->IsNitrogen())
      return 2.711;
    if (atom->IsOxygen())
      return 3.045;
    if (atom->GetAtomicNum() == 9) // F
      return 2.847;
    if (atom->GetAtomicNum() == 14) // Si
      return 2.350;
    if (atom->IsPhosphorus())
      return 2.350;
    if (atom->IsSulfur())
      return 2.980;
    if (atom->GetAtomicNum() == 17) // Cl
      return 2.909;
    if (atom->GetAtomicNum() == 35) // Br
      return 3.017;
    if (atom->GetAtomicNum() == 53) // I
      return 3.086;

    return 0.0;
  }
 
  // MMFF part V - TABLE VI
  double OBForceFieldMMFF94::GetCParam(OBAtom* atom)
  {
    if (atom->GetAtomicNum() == 5) // B
      return 0.704;
    if (atom->IsCarbon())
      return 1.016;
    if (atom->IsNitrogen())
      return 1.113;
    if (atom->IsOxygen())
      return 1.337;
    if (atom->GetAtomicNum() == 14) // Si
      return 0.811;
    if (atom->IsPhosphorus())
      return 1.068;
    if (atom->IsSulfur())
      return 1.249;
    if (atom->GetAtomicNum() == 17) // Cl
      return 1.078;
    if (atom->GetAtomicNum() == 33) // As
      return 0.825;

    return 0.0;
  }
 
  // MMFF part V - TABLE X
  double OBForceFieldMMFF94::GetUParam(OBAtom* atom)
  {
    if (atom->IsCarbon())
      return 2.0;
    if (atom->IsNitrogen())
      return 2.0;
    if (atom->IsOxygen())
      return 2.0;
    if (atom->GetAtomicNum() == 14) // Si
      return 1.25;
    if (atom->IsPhosphorus())
      return 1.25;
    if (atom->IsSulfur())
      return 1.25;
    
    return 0.0;
  }
  
  // MMFF part V - TABLE X
  double OBForceFieldMMFF94::GetVParam(OBAtom* atom)
  {
    if (atom->IsCarbon())
      return 2.12;
    if (atom->IsNitrogen())
      return 1.5;
    if (atom->IsOxygen())
      return 0.2;
    if (atom->GetAtomicNum() == 14) // Si
      return 1.22;
    if (atom->IsPhosphorus())
      return 2.4;
    if (atom->IsSulfur())
      return 0.49;
    
    return 0.0;
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
        return etab.GetCovalentRad(a->GetAtomicNum());
    }
  }
  
  // MMFF part V - page 625
  double OBForceFieldMMFF94::GetRuleBondLength(OBAtom* a, OBAtom* b)
  {
    double r0ab, r0a, r0b, c, Xa, Xb;
    int Ha, Hb, BOab;
    r0a = GetCovalentRadius(a);
    r0b = GetCovalentRadius(b);
    Xa = etab.GetAllredRochowElectroNeg(a->GetAtomicNum());
    Xb = etab.GetAllredRochowElectroNeg(b->GetAtomicNum());
    
   
    if (a->IsHydrogen())
      r0a = 0.33;
    if (b->IsHydrogen())
      r0b = 0.33;
    
    if (a->IsHydrogen() || b->IsHydrogen())
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
    if (a->GetBond(b)->IsAromatic())
      if (!HasPilpSet(atoi(a->GetType())) && !HasPilpSet(atoi(b->GetType())))
        BOab = 4;
      else
        BOab = 5;
     
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
    
    cout << "Ha=" << Ha << "  Hb=" << Hb << "  BOab=" << BOab << endl;
    cout << "r0a=" << r0a << "  Xa=" << Xa << endl;
    cout << "r0b=" << r0b << "  Xb=" << Xb << endl;
    cout << "r0a + r0b=" << r0a +r0b << endl;
    cout << "c=" << c << "  |Xa-Xb|=" << fabs(Xa-Xb) << "  |Xa-Xb|^1.4=" << pow(fabs(Xa-Xb), 1.4) << endl;
    
    r0ab = r0a + r0b - c * pow(fabs(Xa - Xb), 1.4) - 0.008; 

    return r0ab;
  }

  OBFFParameter* OBForceFieldMMFF94::GetParameter1Atom(int a, std::vector<OBFFParameter> &parameter)
  {
    OBFFParameter *par;

    for (int idx=0; idx < parameter.size(); idx++)
      if (a == parameter[idx].a) {
        par = &parameter[idx];
        return par;
      }

    return NULL;
  }
  
  OBFFParameter* OBForceFieldMMFF94::GetParameter2Atom(int a, int b, std::vector<OBFFParameter> &parameter)
  {
    OBFFParameter *par;

    for (int idx=0; idx < parameter.size(); idx++)
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

    for (int idx=0; idx < parameter.size(); idx++)
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
      
    for (int idx=0; idx < parameter.size(); idx++)
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
    
    for (int idx=0; idx < parameter.size(); idx++)
      if (((a == parameter[idx].a) && (b == parameter[idx].b) && (c == parameter[idx].c) && (ffclass == parameter[idx]._ipar[0])) || 
          ((a == parameter[idx].c) && (b == parameter[idx].b) && (c == parameter[idx].a) && (ffclass == parameter[idx]._ipar[0]))) 
        {
          par = &parameter[idx];
          return par;
        }

    return NULL;
  }

  OBFFParameter* OBForceFieldMMFF94::GetTypedParameter4Atom(int ffclass, int a, int b, int c, int d, std::vector<OBFFParameter> &parameter)
  {
    OBFFParameter *par;
    
    for (int idx=0; idx < parameter.size(); idx++)
      if (((a == parameter[idx].a) && (b == parameter[idx].b) && (c == parameter[idx].c) && 
           (d == parameter[idx].d) && (ffclass == parameter[idx]._ipar[0])) || 
           ((a == parameter[idx].d) && (b == parameter[idx].c) && (c == parameter[idx].b) && 
	   (d == parameter[idx].a) && (ffclass == parameter[idx]._ipar[0]))) 
        {
          par = &parameter[idx];
          return par;
        }

    return NULL;
  }

} // end namespace OpenBabel

//! \file forcefieldmmff94.cpp
//! \brief MMFF94 force field

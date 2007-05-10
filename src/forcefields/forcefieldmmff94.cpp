/*********************************************************************
forcefieldmmff94.cpp - Merck Molecular Force Field.

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
  double OBFFBondCalculationMMFF94::GetEnergy()
  {
    double delta2;
    vector3 vab;
     
    vab = a->GetVector() - b->GetVector();
    rab = vab.length();
    delta = rab - r0;
    delta2 = delta * delta;
 
    energy = 143.9325 * 0.5 * kb * delta2 * (1.0 - 2.0 * delta + 7.0/12.0 * 4.0 * delta2);

    return energy;
  }
  
  vector3 OBFFBondCalculationMMFF94::GetGradient(OBAtom *atom) 
  {
    vector3 da, db, gradient;
    double dE;

    if ((atom != a) && (atom != b))
      return  VZero;
     
    da = a->GetVector();
    db = b->GetVector();
    rab = OBForceField::VectorLengthDerivative(da, db);
    delta = rab - r0;
 
    dE = 143.9325 * kb * delta * (1.0 - 2.0 * delta + 7.0/12.0 * 8.0 * delta);

    if (atom == a) {
      gradient = dE * da; // - dE/drab * drab/da
      return gradient;
    } else {
      gradient = dE * db; // - dE/drab * drab/db
      return gradient;
    }
  }

  double OBForceFieldMMFF94::E_Bond()
  {
    vector<OBFFBondCalculationMMFF94>::iterator i;
    double energy = 0.0;

    IF_OBFF_LOGLVL_HIGH {
      *logos << std::endl << "B O N D   S T R E T C H I N G" << endl << endl;
      *logos << "ATOM TYPES   FF    BOND       IDEAL       FORCE" << endl;
      *logos << " I    J     CLASS  LENGTH     LENGTH     CONSTANT      DELTA      ENERGY" << endl;
      *logos << "------------------------------------------------------------------------" << endl;
    }
    
    for (i = _bondcalculations.begin(); i != _bondcalculations.end(); ++i) {

      energy += i->GetEnergy();

      IF_OBFF_LOGLVL_HIGH {
        sprintf(logbuf, "%2d   %2d      %d   %8.3f   %8.3f     %8.3f   %8.3f   %8.3f", atoi((*i).a->GetType()), atoi((*i).b->GetType()), 
                (*i).bt, (*i).rab, (*i).r0, (*i).kb, (*i).delta, (*i).energy);
        *logos  << logbuf << std::endl;
      }
    }
    
    IF_OBFF_LOGLVL_MEDIUM
      *logos << endl << "     TOTAL BOND STRETCHING ENERGY = " << energy << endl << endl;
    return energy;
  }
 
  double OBFFAngleCalculationMMFF94::GetEnergy()
  {
    double delta2;

    theta = a->GetAngle(b->GetIdx(), c->GetIdx());
    
    delta = theta - theta0;
    delta2 = delta * delta;
    
    energy = 0.043844 * 0.5 * ka * delta2 * (1.0 - 0.007 * delta);
     
    return energy;
  }
  
  vector3 OBFFAngleCalculationMMFF94::GetGradient(OBAtom *atom) 
  {
    vector3 da, db, dc, gradient;
    double dE;

    if ((atom->GetIdx() != a->GetIdx()) && (atom->GetIdx() != b->GetIdx()) && (atom->GetIdx() != c->GetIdx()))
      return VZero;
    
    da = a->GetVector();
    db = b->GetVector();
    dc = c->GetVector();
    theta = OBForceField::VectorAngleDerivative(da, db, dc);
    delta = theta - theta0;
 
    dE = 0.043844 * ka * delta * (1.0 - 1.5 * 0.007 * delta);

    if (atom == a) {
      gradient = dE * da; // - dE/drab * drab/da
      return gradient;
    } else if (atom == b) {
      gradient = dE * db; // - dE/drab * drab/db = - dE/drab * drab/da - dE/drab * drab/dc 
      return gradient;
    } else {
      gradient = dE * dc; // - dE/drab * drab/dc
      return gradient;
    }
  }
 
  double OBForceFieldMMFF94::E_Angle()
  {
    vector<OBFFAngleCalculationMMFF94>::iterator i;
    double energy = 0.0;
 
    IF_OBFF_LOGLVL_HIGH {
      *logos << endl << "A N G L E   B E N D I N G" << endl << endl;
      *logos << "ATOM TYPES        FF    VALENCE     IDEAL      FORCE" << endl;
      *logos << " I    J    K     CLASS   ANGLE      ANGLE     CONSTANT      DELTA      ENERGY" << endl;
      *logos << "-----------------------------------------------------------------------------" << endl;
    }
    
    for (i = _anglecalculations.begin(); i != _anglecalculations.end(); ++i) {

      energy += i->GetEnergy();
      
      IF_OBFF_LOGLVL_HIGH {
        sprintf(logbuf, "%2d   %2d   %2d      %d   %8.3f   %8.3f     %8.3f   %8.3f   %8.3f", atoi((*i).a->GetType()), atoi((*i).b->GetType()), atoi((*i).c->GetType()), 
                (*i).at, (*i).theta, (*i).theta0, (*i).ka, (*i).delta, (*i).energy);
        *logos  << logbuf << endl;
      }
    }
 
    IF_OBFF_LOGLVL_MEDIUM
      *logos << endl << "     TOTAL ANGLE BENDING ENERGY = " << energy << endl << endl;
    return energy;
  }
  
  double OBFFStrBndCalculationMMFF94::GetEnergy()
  {
    vector3 vab, vbc;

    vab = a->GetVector() - b->GetVector();
    vbc = b->GetVector() - c->GetVector();
    
    theta = a->GetAngle(b->GetIdx(), c->GetIdx());
    rab = vab.length();
    rbc = vbc.length();
    
    delta_theta = theta - theta0;
    delta_rab = rab - rab0;
    delta_rbc = rbc - rbc0;

    energy = 2.51210 * (kbaABC * delta_rab + kbaCBA * delta_rbc) * delta_theta;

    return energy;
  }
  
  vector3 OBFFStrBndCalculationMMFF94::GetGradient(OBAtom *atom) 
  {
    vector3 da, db, dc, gradient;
    vector3 rab_da, rab_db, rbc_db, rbc_dc, theta_da, theta_db, theta_dc;
    double dE;

    if ((atom->GetIdx() != a->GetIdx()) && (atom->GetIdx() != b->GetIdx()) && (atom->GetIdx() != c->GetIdx()))
      return VZero;
    
    rab_da = theta_da = a->GetVector();
    rab_db = rbc_db = theta_db = b->GetVector();
    rbc_dc = theta_dc = c->GetVector();
    
    rab = OBForceField::VectorLengthDerivative(rab_da, rab_db);
    rbc = OBForceField::VectorLengthDerivative(rbc_db, rbc_dc);
    theta = OBForceField::VectorAngleDerivative(theta_da, theta_db, theta_dc);
    delta_rab = rab - rab0;
    delta_rbc = rbc - rbc0;
    delta_theta = theta - theta0;

    da = 2.51210 * (kbaABC * rab_da * delta_theta + theta_da * (kbaABC * delta_rab + kbaCBA * delta_rbc));
    dc = 2.51210 * (kbaCBA * rbc_dc * delta_theta + theta_dc * (kbaABC * delta_rab + kbaCBA * delta_rbc));
    db = -da - dc;

    if (atom == a) {
      return da;
    } else if (atom == b) {
      return db;
    } else {
      return dc;
    }
  }
 
  double OBForceFieldMMFF94::E_StrBnd() 
  {
    vector<OBFFStrBndCalculationMMFF94>::iterator i;
    double energy = 0.0;
    
    IF_OBFF_LOGLVL_HIGH {
      *logos  << std::endl << "S T R E T C H   B E N D I N G" << std::endl << std::endl;
      *logos << "ATOM TYPES        FF    VALENCE     DELTA        FORCE CONSTANT " << std::endl;
      *logos << " I    J    K     CLASS   ANGLE      ANGLE        I J        J K      ENERGY" << std::endl;
      *logos << "---------------------------------------------------------------------------" << std::endl;
    }
    
    for (i = _strbndcalculations.begin(); i != _strbndcalculations.end(); ++i) {

      energy += i->GetEnergy();
      
      IF_OBFF_LOGLVL_HIGH {
        sprintf(logbuf, "%2d   %2d   %2d     %2d   %8.3f   %8.3f   %8.3f   %8.3f   %8.3f", atoi((*i).a->GetType()), atoi((*i).b->GetType()), atoi((*i).c->GetType()), 
                (*i).sbt, (*i).theta, (*i).delta_theta, (*i).kbaABC, (*i).kbaCBA, (*i).energy);
        *logos  << logbuf << endl;
      }
    }
	
    IF_OBFF_LOGLVL_MEDIUM
      *logos << endl << "     TOTAL STRETCH BENDING ENERGY = " << energy << endl << endl;
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

  double OBFFTorsionCalculationMMFF94::GetEnergy()
  {
    double cosine, cosine2, cosine3;
    double phi1, phi2, phi3;

    tor = CalcTorsionAngle(a->GetVector(), b->GetVector(), c->GetVector(), d->GetVector());

    cosine = cos(DEG_TO_RAD * tor);
    cosine2 = cos(DEG_TO_RAD * 2 * tor);
    cosine3 = cos(DEG_TO_RAD * 3 * tor);
      
    phi1 = 1.0 + cosine;
    phi2 = 1.0 - cosine2;
    phi3 = 1.0 + cosine3;
    
    energy = 0.5 * (v1 * phi1 + v2 * phi2 + v3 * phi3);
 
    return energy;
  }
  
  vector3 OBFFTorsionCalculationMMFF94::GetGradient(OBAtom *atom) 
  {
    vector3 da, db, dc, dd, gradient;
    double dE, sine, sine2, sine3;

    if ((atom->GetIdx() != a->GetIdx()) && (atom->GetIdx() != b->GetIdx()) && (atom->GetIdx() != c->GetIdx()) &&  (atom->GetIdx() != d->GetIdx()))
      return  VZero;
    
    da = a->GetVector();
    db = b->GetVector();
    dc = c->GetVector();
    dd = d->GetVector();
    tor = OBForceField::VectorTorsionDerivative(da, db, dc, dd);
 
    sine = sin(DEG_TO_RAD * (tor));
    sine2 = sin(DEG_TO_RAD * (2.0 * tor));
    sine3 = sin(DEG_TO_RAD * (3.0 * tor));
    dE = -0.5 * (v1 * sine - 2.0 * v2 * sine2 + 3.0 * v3 * sine3);
    
    //if (IsNearZero(dE))
    //  return VZero;

    if (atom == a) {
      gradient = dE * da; // - dE/drab * drab/da
      return gradient;
    } else if (atom == b) {
      gradient = dE * db; // - dE/drab * drab/db
      return gradient;
    } else if (atom == c) {
      gradient = dE * dc; // - dE/drab * drab/dc
      return gradient;
    } else {
      gradient = dE * dd; // - dE/drab * drab/dd
      return gradient;
    }
  }
 
  double OBForceFieldMMFF94::E_Torsion() 
  {
    vector<OBFFTorsionCalculationMMFF94>::iterator i;
    double energy = 0.0;
        
    IF_OBFF_LOGLVL_HIGH {
      *logos << endl << "T O R S I O N A L" << endl << endl;
      *logos << "ATOM TYPES             FF     TORSION       FORCE CONSTANT" << endl;
      *logos << " I    J    K    L     CLASS    ANGLE         V1   V2   V3     ENERGY" << endl;
      *logos << "--------------------------------------------------------------------" << endl;
    }

    for (i = _torsioncalculations.begin(); i != _torsioncalculations.end(); ++i) {

      energy += i->GetEnergy();
      
      IF_OBFF_LOGLVL_HIGH {
        sprintf(logbuf, "%2d   %2d   %2d   %2d      %d   %8.3f   %6.3f   %6.3f   %6.3f   %8.3f",  atoi((*i).a->GetType()), atoi((*i).b->GetType()), 
                atoi((*i).c->GetType()), atoi((*i).d->GetType()), (*i).tt, (*i).tor, (*i).v1, (*i).v2, (*i).v3, (*i).energy);
        *logos  << logbuf << endl;
      }
    }

    IF_OBFF_LOGLVL_MEDIUM
      *logos << std::endl << "     TOTAL TORSIONAL ENERGY = " << energy << std::endl << std::endl;
    return energy;
  }
 
  //
  //  a
  //   \
  //    b---d      plane = a-b-c
  //   /
  //  c
  //
  double OBFFOOPCalculationMMFF94::GetEnergy()
  {
    double angle2;

    angle = Point2PlaneAngle(d->GetVector(), a->GetVector(), b->GetVector(), c->GetVector());
    angle2 = angle * angle;
    
    energy = 0.043844 * 0.5 * koop * angle2;

    return energy;
  }

  double OBForceFieldMMFF94::E_OOP() 
  {
    vector<OBFFOOPCalculationMMFF94>::iterator i;
    double energy = 0.0;
    
    IF_OBFF_LOGLVL_HIGH {
      *logos << std::endl << "O U T - O F - P L A N E   B E N D I N G" << std::endl << std::endl;
      *logos << "ATOM TYPES             FF       OOP     FORCE " << std::endl;
      *logos << " I    J    K    L     CLASS    ANGLE   CONSTANT     ENERGY" << std::endl;
      *logos << "----------------------------------------------------------" << std::endl;
    }

    for (i = _oopcalculations.begin(); i != _oopcalculations.end(); ++i) {
      energy += i->GetEnergy();
      
      IF_OBFF_LOGLVL_HIGH {
        sprintf(logbuf, "%2d   %2d   %2d   %2d      0   %8.3f   %8.3f     %8.3f", atoi((*i).a->GetType()), atoi((*i).b->GetType()), atoi((*i).c->GetType()), atoi((*i).d->GetType()), 
                (*i).angle, (*i).koop, (*i).energy);
        *logos  << logbuf << endl;
      }
    }

    IF_OBFF_LOGLVL_HIGH
      *logos << std::endl << "     TOTAL OUT-OF-PLANE BENDING ENERGY = " << energy << std::endl << std::endl;
    return energy;
  }
 
  double OBFFVDWCalculationMMFF94::GetEnergy()
  {
    double rab2, rab3, rab4, rab5, rab6, rab7;
    double erep2, erep3, erep4, erep5, erep6, erep7;

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
      
    return energy;
  }
  
  double OBForceFieldMMFF94::E_VDW()
  {
    vector<OBFFVDWCalculationMMFF94>::iterator i;
    double energy = 0.0;
    
    IF_OBFF_LOGLVL_HIGH {
      *logos << std::endl << "V A N   D E R   W A A L S" << std::endl << std::endl;
      *logos << "ATOM TYPES          " << std::endl;
      *logos << " I    J        Rij       R*IJ    EPSILON    E_REP     E_ATTR    ENERGY" << std::endl;
      *logos << "----------------------------------------------------------------------" << std::endl;
      //          XX   XX     -000.000  -000.000  -000.000  -000.000  -000.000  -000.000
    }
    
    for (i = _vdwcalculations.begin(); i != _vdwcalculations.end(); ++i) {
      
      energy += i->GetEnergy();
      
      IF_OBFF_LOGLVL_HIGH {
        sprintf(logbuf, "%2d   %2d     %8.3f  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f", atoi((*i).a->GetType()), atoi((*i).b->GetType()), 
                (*i).rab, (*i).R_AB, (*i).epsilon, (*i).erep, (*i).eattr, (*i).energy);
        *logos  << logbuf << endl;
      }
    }

    IF_OBFF_LOGLVL_MEDIUM
      *logos << std::endl << "     TOTAL VAN DER WAALS ENERGY = " << energy << std::endl << std::endl;
    return energy;
  }

  /// \todo Currently unimplemented but needs to be here for linking
  double OBFFElectrostaticCalculationMMFF94::GetEnergy()
  {
    return 0.0;
  }

  /// \todo Currently unimplemented but needs to be here for linking
  double OBForceFieldMMFF94::E_Electrostatic()
  {
    OBAtom *a, *b, *c, *d;
    double e, energy = 0.0;
    
    IF_OBFF_LOGLVL_HIGH {
      *logos << std::endl << "E L E C T R O S T A T I C   I N T E R A C T I O N S" << std::endl << std::endl;
      *logos << "ATOM TYPES          " << std::endl;
      *logos << " I    J        Rij       R*IJ    EPSILON    E_REP     E_ATTR    ENERGY" << std::endl;
      *logos << "----------------------------------------------------------------------" << std::endl;
      //            XX   XX     -000.000  -000.000  -000.000  -000.000  -000.000  -000.000
    }
    
    FOR_PAIRS_OF_MOL(p, _mol) {
    }

    IF_OBFF_LOGLVL_MEDIUM
      *logos << std::endl << "     TOTAL ELECTROSTATIC ENERGY = " << energy << std::endl << std::endl;
    return energy;
  }
  /*
    vector<vector3> OBForceFieldMMFF94::GetForces()
    {
    vector<OBFFBondCalculationMMFF94>::iterator i;
    vector<OBFFAngleCalculationMMFF94>::iterator i2;
    vector<vector3> forces;
    vector3 anagrad;
    
    FOR_ATOMS_OF_MOL (a, _mol) {
    anagrad.Set(0,0,0);

    for (i = _bondcalculations.begin(); i != _bondcalculations.end(); ++i) {
    if (((*i).a->GetIdx() == a->GetIdx()) || ((*i).b->GetIdx() == a->GetIdx())) {
    anagrad += i->GetGradient(&*a);
    }
    }

    for (i2 = _anglecalculations.begin(); i2 != _anglecalculations.end(); i2++) {
    if (((*i2).a->GetIdx() == a->GetIdx()) || ((*i2).b->GetIdx() == a->GetIdx()) || ((*i2).c->GetIdx() == a->GetIdx())) {
    anagrad += i2->GetGradient(&*a);
    }
    }

    forces.push_back(anagrad);
    }

    return forces; 
    }
  */
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
    return *this;
  }

  bool OBForceFieldMMFF94::Setup(OBMol &mol)
  {
    _mol = mol;
    SetMMFFTypes();
    if (!SetupCalculations())
      return false;

    //SetMMFF94Charges();
    
    return true;
  }
 
  bool OBForceFieldMMFF94::ParseParamFile()
  {
    ParseParamProp();
    ParseParamBond();
    ParseParamBndk();
    ParseParamAngle();
    ParseParamStrBnd();
    ParseParamDfsb();
    ParseParamOOP();
    ParseParamTorsion();
    ParseParamVDW();
    ParseParamCharge();
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
      parameter.ipar5 = atoi(vs[0].c_str());  // FF class
      parameter.a = atoi(vs[1].c_str());
      parameter.b = atoi(vs[2].c_str());
      parameter.dpar1 = atof(vs[3].c_str());  // kb
      parameter.dpar2 = atof(vs[4].c_str());  // r0
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
      parameter.a = atoi(vs[1].c_str());
      parameter.b = atoi(vs[2].c_str());
      parameter.dpar1 = atof(vs[3].c_str());  // r0-ref
      parameter.dpar2 = atof(vs[4].c_str());  // kb-ref
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
      parameter.ipar5 = atoi(vs[0].c_str());  // FF class
      parameter.a = atoi(vs[1].c_str());
      parameter.b = atoi(vs[2].c_str());
      parameter.c = atoi(vs[3].c_str());
      parameter.dpar1 = atof(vs[4].c_str());  // ka
      parameter.dpar2 = atof(vs[5].c_str());  // theta0
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
      parameter.ipar5 = atoi(vs[0].c_str());  // FF class
      parameter.a = atoi(vs[1].c_str());
      parameter.b = atoi(vs[2].c_str());
      parameter.c = atoi(vs[3].c_str());
      parameter.dpar1 = atof(vs[4].c_str());  // kbaIJK
      parameter.dpar2 = atof(vs[5].c_str());  // kbaKJI
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
      parameter.dpar1 = atof(vs[3].c_str());  // kbaIJK
      parameter.dpar2 = atof(vs[4].c_str());  // kbaKJI
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
      parameter.dpar1 = atof(vs[4].c_str());  // koop
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
      parameter.ipar5 = atoi(vs[0].c_str());  // FF class
      parameter.a = atoi(vs[1].c_str());
      parameter.b = atoi(vs[2].c_str());
      parameter.c = atoi(vs[3].c_str());
      parameter.d = atoi(vs[4].c_str());
      parameter.dpar1 = atof(vs[5].c_str());  // v1
      parameter.dpar2 = atof(vs[6].c_str());  // v2
      parameter.dpar3 = atof(vs[7].c_str());  // v3
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
      parameter.dpar1 = atof(vs[1].c_str());  // alpha-i
      parameter.dpar2 = atof(vs[2].c_str());  // N-i
      parameter.dpar3 = atof(vs[3].c_str());  // A-i
      parameter.dpar4 = atof(vs[4].c_str());  // G-i
      if (EQn(vs[5].c_str(), "-", 1))
        parameter.ipar1 = 0;
      else if (EQn(vs[5].c_str(), "D", 1))
        parameter.ipar1 = 1;  // hydrogen bond donor
      else if (EQn(vs[5].c_str(), "A", 1))
        parameter.ipar1 = 2;  // hydrogen bond acceptor
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
      parameter.a = atoi(vs[1].c_str());
      parameter.b = atoi(vs[2].c_str());
      parameter.dpar1 = atof(vs[3].c_str());  // bci
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
      parameter.dpar1 = atof(vs[2].c_str());  // pbci
      parameter.dpar2 = atof(vs[3].c_str());  // fcadj
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
    
    // open data/mmffprop.par
    ifstream ifs;
    if (OpenDatafile(ifs, "mmffprop.par").length() == 0) {
      obErrorLog.ThrowError(__FUNCTION__, "Cannot open mmffprop.par", obError);
      return false;
    }
    
    _arom.clear();
    _lin.clear();
    _sbmb.clear();
    while (ifs.getline(buffer, 80)) {
      if (EQn(buffer, "*", 1)) continue;
      if (EQn(buffer, "$", 1)) continue;
	
      tokenize(vs, buffer);

      if (atoi(vs[6].c_str()) == 1)
        _arom.push_back(atoi(vs[0].c_str()));

      if (atoi(vs[7].c_str()) == 1)
        _lin.push_back(atoi(vs[0].c_str()));

      if (atoi(vs[8].c_str()) == 1)
        _sbmb.push_back(atoi(vs[0].c_str()));
    }

    if (ifs)
      ifs.close();
 
    return 0;
  }

  bool OBForceFieldMMFF94::SetMMFFTypes()
  {
    std::vector<std::vector<int> > _mlist; //!< match list for atom typing
    std::vector<std::pair<OBSmartsPattern*,std::string> > _vexttyp; //!< external atom type rules
    vector<vector<int> >::iterator j;
    vector<pair<OBSmartsPattern*,string> >::iterator i;
    OBSmartsPattern *sp;
    vector<string> vs;
    char buffer[150];
 
    *logos << std::endl << "A T O M   T Y P E S" << std::endl << std::endl;
    _mol.SetAtomTypesPerceived();
    
    // DEBUG
    FOR_BONDS_OF_MOL (bond, _mol)
      if (bond->IsDouble())
        cout << bond->GetBeginAtom()->GetIdx() << "=" << bond->GetEndAtom()->GetIdx() << endl;

    cout << "kekulize..." << endl;
    _mol.Kekulize();
    _mol.UnsetAromaticPerceived();
    
    cout << endl << "B O N D S" << endl << endl;
    cout << " I    J     BO     AR" << endl;
    cout << "---------------------" << endl;
    //      "XX   XX     XX     XX
    
    FOR_BONDS_OF_MOL (bond, _mol) {
      sprintf(buffer, "%2d   %2d     %d     %d", bond->GetBeginAtom()->GetIdx(), bond->GetEndAtom()->GetIdx(), bond->GetBO(), bond->IsAromatic());
      cout  << buffer << endl;
    }
   
    ////////////////////////////////////////////////////////////////////////////
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
        obErrorLog.ThrowError(__FUNCTION__, " Could not parse EXTTYP line in atom type table from atomtyp.txt", obInfo);
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
    // find aromatic rings and assign atom types
    ////////////////////////////////////////////////////////////////////////////

    cout << endl << "R I N G S" << endl << endl;
    //      "XX   XX     XX     XX
    
    OBAtom *ringatom;
    OBBond *ringbond;
    vector<OBRing*> vr;
    vr = _mol.GetSSSR();
    
    vector<OBRing*>::iterator ri;
    vector<int>::iterator rj;
    int n, index, ringsize, first_rj, prev_rj, pi_electrons;
    for (ri = vr.begin();ri != vr.end();++ri) { // for each ring
      ringsize = (*ri)->Size();
      cout << "RING: ";
      
      n = 1;
      pi_electrons = 0;
      for(rj = (*ri)->_path.begin();rj != (*ri)->_path.end();rj++) { // for each ring atom
        index = *rj;
        ringatom = _mol.GetAtom(index);
        //cout << index << ": " << ringatom << endl; 
        
        // is the bond to the previous ring atom double?
        if (n > 1) {
          ringbond = _mol.GetBond(prev_rj, index);
          if (ringbond->GetBO() == 2) {
            pi_electrons += 2;
            cout << "=" << index;
            prev_rj = index;
            n++;
            continue;
          } else 
            cout << "-" << index;
          prev_rj = index;
        } else {
          cout << index;
          prev_rj = index;
          first_rj = index;
        }
	
        // does the current ring atom have a exocyclic double bond?
        FOR_NBORS_OF_ATOM (nbr, ringatom) {
          if ((*ri)->IsInRing(nbr->GetIdx()))
            continue;
          ringbond = _mol.GetBond(nbr->GetIdx(), index);
          if (ringbond->GetBO() == 2) {
            pi_electrons++;
            cout << "(=*)";
          }
        }

        // is the atom N, O or S in 5 rings
        if (ringsize == 5) {
          if (ringatom->IsNitrogen() && !(ringatom->BOSum() - ringatom->GetValence())) {
            pi_electrons += 2;
            cout << "(:)";
          }
          if ((ringatom->IsOxygen() || ringatom->IsSulfur()) && ringatom->IsAromatic()) {
            pi_electrons += 2;
            cout << "(:)";
          }
        }

        n++;
      
      } // for each ring atom
      
      // is the bond from the first to the last atom double?
      ringbond = _mol.GetBond(first_rj, index);
      if (ringbond->GetBO() == 2) {
        pi_electrons += 2;
        cout << "=";
      }

      // check if the ring is aromatic
      if (pi_electrons == 6) {
        cout << " size=" << ringsize << " pi=" << pi_electrons << "    AROMATIC!" << endl;
        for(rj = (*ri)->_path.begin();rj != (*ri)->_path.end();rj++) { // for each ring atom
          index = *rj;
          ringatom = _mol.GetAtom(index);
	  
          if (ringsize == 6) {
            if (ringatom->IsCarbon()) 
              ringatom->SetType("37"); // CB: CARBON AS IN BENZENE, PYRROLE

            if (ringatom->IsNitrogen()) {
              if (ringatom->GetValence() == 2)
                ringatom->SetType("38"); // NPYD: NITROGEN, AS IN PYRIDINE
              if (ringatom->GetValence() == 3) {
                ringatom->SetType("58"); // NPD+: PYRIDINIUM-TYPE NITROGEN - FORMAL CHARGE=1
	
                FOR_NBORS_OF_ATOM (nbr, ringatom) {
                  if ((*ri)->IsInRing(nbr->GetIdx()))
                    continue;
                  if (nbr->IsOxygen())
                    ringatom->SetType("69"); // NPOX: PYRIDINE N-OXIDE NITROGEN
                }
              }
            }

          }

        }
      } else {
        cout << " size=" << ringsize << " pi=" << pi_electrons << endl;
      }
    } // for each ring

    IF_OBFF_LOGLVL_MEDIUM {
      *logos << "IDX\tTYPE" << std::endl;
      FOR_ATOMS_OF_MOL (a, _mol)
        *logos << a->GetIdx() << "\t" << a->GetType() << std::endl;
    }
    
    return true;
  }
  
  bool OBForceFieldMMFF94::SetupCalculations()
  {
    OBFFParameter *parameter;
    OBAtom *a, *b, *c, *d;
    bool found;
    
    // 
    // Bond Calculations
    //
    OBFFBondCalculationMMFF94 bondcalc;
    int bondtype;

    _bondcalculations.clear();
    
    FOR_BONDS_OF_MOL(bond, _mol) {
      a = bond->GetBeginAtom();
      b = bond->GetEndAtom();	
      bondtype = GetBondType(a, b);
      
      cout << bondtype << " idx: " << a->GetIdx() << "-" << b->GetIdx() << " types: " << a->GetType() << "-" << b->GetType() << endl;

      parameter = GetParameterMMFF94(bondtype, atoi(a->GetType()), atoi(b->GetType()), 0, 0, _ffbondparams);
      if (parameter == NULL) {
        
        //parameter = GetParameterMMFF94(1, atoi(a->GetType()), atoi(b->GetType()), 0, 0, _ffbondparams);
        //parameter = GetParameter(atoi(a->GetType()), atoi(b->GetType()), 0, 0, _ffbndkparams);
        if (parameter == NULL) {
          obErrorLog.ThrowError(__FUNCTION__, "Could not find all bond parameters", obError);
          return false;
        }
      } else {
        bondcalc.a = a;
        bondcalc.b = b;
        bondcalc.kb = parameter->dpar1;
        bondcalc.r0 = parameter->dpar2;
        bondcalc.bt = bondtype;

        _bondcalculations.push_back(bondcalc);
      }
    }

    //
    // Angle & StrBnd Calculations
    //
    OBFFAngleCalculationMMFF94 anglecalc;
    OBFFStrBndCalculationMMFF94 strbndcalc;
    int angletype, strbndtype, bondtype1, bondtype2;
 
    _anglecalculations.clear();
    _strbndcalculations.clear();
    
    FOR_ANGLES_OF_MOL(angle, _mol) {
      b = _mol.GetAtom((*angle)[0] + 1);
      a = _mol.GetAtom((*angle)[1] + 1);
      c = _mol.GetAtom((*angle)[2] + 1);
      angletype = GetAngleType(a, b, c);
      strbndtype = GetStrBndType(a, b, c);
      bondtype1 = GetBondType(a, b);
      bondtype2 = GetBondType(b, c);
      
      if (HasLinSet(atoi(b->GetType())))
        continue;


      parameter = GetParameterMMFF94(angletype, atoi(a->GetType()), atoi(b->GetType()), atoi(c->GetType()), 0, _ffangleparams);
      if (parameter == NULL) {
        parameter = GetParameterMMFF94(angletype, 0, atoi(b->GetType()), 0, 0, _ffangleparams);
        if (parameter == NULL) {
          cout << angletype << " idx: " << a->GetIdx() << "-" << b->GetIdx() << "-" << c->GetIdx() << " types: " << a->GetType() << "-" << b->GetType() << "-" << c->GetType() << endl;
          obErrorLog.ThrowError(__FUNCTION__, "Could not find all angle parameters", obError);
          return false;
        }
      }
      
      anglecalc.a = a;
      anglecalc.b = b;
      anglecalc.c = c;
      anglecalc.ka   = parameter->dpar1;
      anglecalc.theta0 = parameter->dpar2;
      strbndcalc.theta0 = parameter->dpar2; // **
      anglecalc.at = angletype;
      
      _anglecalculations.push_back(anglecalc);

      parameter = GetParameterMMFF94(strbndtype, atoi(a->GetType()), atoi(b->GetType()), atoi(c->GetType()), 0, _ffstrbndparams);
      if (parameter == NULL) {
        int rowa, rowb, rowc;

        rowa = GetElementRow(a);
        rowb = GetElementRow(b);
        rowc = GetElementRow(c);
        
        parameter = GetParameter(rowa, rowb, rowc, 0, _ffdfsbparams);
        if (parameter == NULL) {
          obErrorLog.ThrowError(__FUNCTION__, "Could not find all stretch-bend parameters", obError);
          return false;
        }

        if (rowa == parameter->a) {
          strbndcalc.kbaABC = parameter->dpar1;
          strbndcalc.kbaCBA = parameter->dpar2;
        } else {
          strbndcalc.kbaABC = parameter->dpar2;
          strbndcalc.kbaCBA = parameter->dpar1;
        }
      } else {
        if (atoi(a->GetType()) == parameter->a) {
          strbndcalc.kbaABC = parameter->dpar1;
          strbndcalc.kbaCBA = parameter->dpar2;
        } else {
          strbndcalc.kbaABC = parameter->dpar2;
          strbndcalc.kbaCBA = parameter->dpar1;
        }
      }
 
      parameter = GetParameterMMFF94(bondtype1, atoi(a->GetType()), atoi(b->GetType()), 0, 0, _ffbondparams);
      strbndcalc.rab0 = parameter->dpar2;
      
      parameter = GetParameterMMFF94(bondtype2, atoi(b->GetType()), atoi(c->GetType()), 0, 0, _ffbondparams);
      strbndcalc.rbc0 = parameter->dpar2;

      strbndcalc.a = a;
      strbndcalc.b = b;
      strbndcalc.c = c;
      strbndcalc.sbt = strbndtype;

      _strbndcalculations.push_back(strbndcalc);
 
    }

    //
    // Torsion Calculations
    //
    OBFFTorsionCalculationMMFF94 torsioncalc;
    int torsiontype;

    _torsioncalculations.clear();
 
    FOR_TORSIONS_OF_MOL(t, _mol) {
      a = _mol.GetAtom((*t)[0] + 1);
      b = _mol.GetAtom((*t)[1] + 1);
      c = _mol.GetAtom((*t)[2] + 1);
      d = _mol.GetAtom((*t)[3] + 1);
      torsiontype = GetTorsionType(a, b, c, d);
      
      if (HasLinSet(atoi(b->GetType())) || HasLinSet(atoi(c->GetType())))
        continue;

      parameter = GetParameterMMFF94(torsiontype, atoi(a->GetType()), atoi(b->GetType()), atoi(c->GetType()), atoi(d->GetType()), _fftorsionparams);

      found = false;
      parameter = GetParameterMMFF94(torsiontype, atoi(a->GetType()), atoi(b->GetType()), atoi(c->GetType()), atoi(d->GetType()), _fftorsionparams);
      if (parameter == NULL) {
        for (int idx=0; idx < _fftorsionparams.size(); idx++) {  // *-XX-XX-XX 
          if (((_fftorsionparams[idx].a == 0) && (atoi(b->GetType()) == _fftorsionparams[idx].b) && (torsiontype == _fftorsionparams[idx].ipar5) &&
               (atoi(c->GetType()) == _fftorsionparams[idx].c) && (atoi(d->GetType()) == _fftorsionparams[idx].d)) ||
              ((_fftorsionparams[idx].a == 0) && (atoi(b->GetType()) == _fftorsionparams[idx].c) && (torsiontype == _fftorsionparams[idx].ipar5) &&
               (atoi(c->GetType()) == _fftorsionparams[idx].b) && (atoi(d->GetType()) == _fftorsionparams[idx].d)) ||
              ((_fftorsionparams[idx].a == 0) && (atoi(b->GetType()) == _fftorsionparams[idx].b) && (torsiontype == _fftorsionparams[idx].ipar5) &&
               (atoi(c->GetType()) == _fftorsionparams[idx].c) && (atoi(d->GetType()) == _fftorsionparams[idx].a)) ||
              ((_fftorsionparams[idx].a == 0) && (atoi(b->GetType()) == _fftorsionparams[idx].c) && (torsiontype == _fftorsionparams[idx].ipar5) &&
               (atoi(c->GetType()) == _fftorsionparams[idx].b) && (atoi(d->GetType()) == _fftorsionparams[idx].a)))
            {
              torsioncalc.v1 = _fftorsionparams[idx].dpar1;
              torsioncalc.v2 = _fftorsionparams[idx].dpar2;
              torsioncalc.v3 = _fftorsionparams[idx].dpar3;
              found = true;  
            }
        }

        if (!found)
          for (int idx=0; idx < _fftorsionparams.size(); idx++) {  // *-XX-XX-*
            if (((_fftorsionparams[idx].a == 0) && (atoi(b->GetType()) == _fftorsionparams[idx].b) &&  (torsiontype == _fftorsionparams[idx].ipar5) &&
                 (atoi(c->GetType()) == _fftorsionparams[idx].c) && (_fftorsionparams[idx].d == 0)) ||
                ((_fftorsionparams[idx].a == 0) && (atoi(c->GetType()) == _fftorsionparams[idx].b) &&  (torsiontype == _fftorsionparams[idx].ipar5) &&
                 (atoi(b->GetType()) == _fftorsionparams[idx].c) && (_fftorsionparams[idx].d == 0)))
              {
                torsioncalc.v1 = _fftorsionparams[idx].dpar1;
                torsioncalc.v2 = _fftorsionparams[idx].dpar2;
                torsioncalc.v3 = _fftorsionparams[idx].dpar3;
                found = true;  
              }
          }

        if (!found) {
          obErrorLog.ThrowError(__FUNCTION__, "Could not find all torsion parameters", obError);
          cout << torsiontype << " idx: " << a->GetIdx() << "-" << b->GetIdx() << "-" << c->GetIdx() << "-" << d->GetIdx() << " types: " << a->GetType() << "-" << b->GetType() << "-" << c->GetType() << "-" << d->GetType() << endl;
          return false;
        }
      } else {
        torsioncalc.v1 = parameter->dpar1;
        torsioncalc.v2 = parameter->dpar2;
        torsioncalc.v3 = parameter->dpar3;
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

          if (((atoi(a->GetType()) == _ffoopparams[idx].a) && (atoi(c->GetType()) == _ffoopparams[idx].c) && (atoi(d->GetType()) == _ffoopparams[idx].d)) ||
              ((atoi(c->GetType()) == _ffoopparams[idx].a) && (atoi(a->GetType()) == _ffoopparams[idx].c) && (atoi(d->GetType()) == _ffoopparams[idx].d)) ||
              ((atoi(c->GetType()) == _ffoopparams[idx].a) && (atoi(d->GetType()) == _ffoopparams[idx].c) && (atoi(a->GetType()) == _ffoopparams[idx].d)) ||
              ((atoi(d->GetType()) == _ffoopparams[idx].a) && (atoi(c->GetType()) == _ffoopparams[idx].c) && (atoi(a->GetType()) == _ffoopparams[idx].d)) ||
              ((atoi(a->GetType()) == _ffoopparams[idx].a) && (atoi(d->GetType()) == _ffoopparams[idx].c) && (atoi(c->GetType()) == _ffoopparams[idx].d)) ||
              ((atoi(d->GetType()) == _ffoopparams[idx].a) && (atoi(a->GetType()) == _ffoopparams[idx].c) && (atoi(c->GetType()) == _ffoopparams[idx].d)))
            {
              found = true;

              oopcalc.koop = _ffoopparams[idx].dpar1;
            
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
              oopcalc.koop = _ffoopparams[idx].dpar1;
	    
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
    OBFFVDWCalculationMMFF94 vdwcalc;

    _vdwcalculations.clear();
 
    FOR_PAIRS_OF_MOL(p, _mol) {
      a = _mol.GetAtom((*p)[0]);
      b = _mol.GetAtom((*p)[1]);

      OBFFParameter *parameter_a, *parameter_b;
      parameter_a = GetParameter(atoi(a->GetType()), 0, 0, 0, _ffvdwparams);
      parameter_b = GetParameter(atoi(b->GetType()), 0, 0, 0, _ffvdwparams);
      if ((parameter_a == NULL) || (parameter_b == NULL)) {
        obErrorLog.ThrowError(__FUNCTION__, "Could not find all Van der Waals parameters", obError);
        return false;
      }
      
      vdwcalc.a = a;
      vdwcalc.alpha_a = parameter_a->dpar1;
      vdwcalc.Na = parameter_a->dpar2;
      vdwcalc.Aa = parameter_a->dpar3;
      vdwcalc.Ga = parameter_a->dpar4;
      vdwcalc.aDA = parameter_a->ipar1;
      
      vdwcalc.b = b;
      vdwcalc.alpha_b = parameter_b->dpar1;
      vdwcalc.Nb = parameter_b->dpar2;
      vdwcalc.Ab = parameter_b->dpar3;
      vdwcalc.Gb = parameter_b->dpar4;
      vdwcalc.bDA = parameter_b->ipar1;
      
      //these calculations only need to be done once for each pair, 
      //we do them now and save them for later use
      double R_AA, R_BB, R_AB6, g_AB, g_AB2;
      double R_AB2, R_AB4, R_AB7, sqrt_a, sqrt_b;
 
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

    return true;
  }

  bool OBForceFieldMMFF94::SetMMFF94Charges()
  {
    OBFFParameter *parameter_a, *parameter_b, *parameter_ab;
    OBAtom *a;
    double q0i;
    
    FOR_ATOMS_OF_MOL (atom, _mol) {
      a = (OBAtom*) &*atom;
      parameter_a = GetParameter(atoi(atom->GetType()), 0, 0, 0, _ffpbciparams);

      if (atom->IsCarboxylOxygen())
        q0i = -0.500;
      else if (atom->IsPhosphateOxygen() && atom->GetHvyValence() == 1)
        q0i = -0.666;
      else if (atom->IsSulfateOxygen())
        q0i = -0.500;
 
      *logos << "fcharge = " << q0i << std::endl;
      //FOR_NBORS_OF_ATOM (nbr, a) {
      // q =  
      //}
    }
    return true;
  }

  double OBForceFieldMMFF94::Energy()
  {
    double energy;

    //energy = E_Bond();
    energy += E_Angle();
    //energy += E_StrBnd();
    //energy += E_OOP();
    //energy += E_Torsion();
    //energy += E_VDW();
    //energy += E_Electrostatic();

    return energy;
  }

  // used to validate the implementation
  bool OBForceFieldMMFF94::Validate ()
  {
    OBConversion conv;
    OBFormat *format_in = conv.FindFormat("mol2");
    vector<string> vs;
    vector<int> types;
    vector<double> bond_lengths;
    char buffer[150], logbuf[100];
    bool molfound, atomfound, bondfound;
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
    
    SetLogFile(&ofs);
    SetLogLevel(OBFF_LOGLVL_HIGH);
   
    for (unsigned int c=1;; c++) {
      _mol.Clear();
      types.clear();
      bond_lengths.clear();

      if (!conv.Read(&_mol, &ifs))
        break;
      if (_mol.Empty())
        break;
      
      SetMMFFTypes();
      
      n = _mol.NumAtoms() / 4;
      termcount = 0;
      molfound = false;
      atomfound = false;
      bondfound = false;

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
        
        if (molfound && EQn(buffer, " ATOM NAME  TYPE", 16))
          atomfound = true;
	
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

      for (di = bond_lengths.begin(); di != bond_lengths.end(); di++)
        cout << "rab = " << *di << endl;

      cout << "--------------------------------------------------------------------------------" << endl;
      cout << "                                                                                " << endl;
      cout << "  VALIDATE MOLECULE : " << _mol.GetTitle() << endl;
      cout << "                                                                                " << endl;
      cout << "IDX  HYB  AROM  OB_TYPE  LOG_TYPE      RESULT                                              " << endl;
      cout << "---------------------------------------                                              " << endl;
            
      ni = 1;
      failed = false;
      for (i = types.begin(); i != types.end();i++) {
        if (ni > _mol.NumAtoms())
          continue;

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
        //return false;
        continue;
      }

      if (!SetupCalculations()) {
        cout << "Could not setup calculations (missing parameters...)" << endl;
        //return false;
        continue;
      }

      // the calculated bond length has an error compared to the one listed in the validation suite logfile. 
      //
      // rlog    (bond length from log file)
      // rab     (bond length calculated by openbabel)
      // error = fabs(rab - rlog)
      //
      //          nbonds                           If this deltaE is equal to the difference in listed
      //           ----                            energy and calculated, we can assume that the 
      // deltaE =  \    E_bond(error)              implementation is correct. The observed errors are a
      //           /                               result of calculations, not the implementation.
      //           ----
      //            n=1
      vector<OBFFBondCalculationMMFF94>::iterator k;
      double energy;
      int idx;

      idx = 0;
      energy = 0.0;
      for (k = _bondcalculations.begin(); k != _bondcalculations.end(); ++k) {
        vector3 vab;
        double delta2;
     
        vab = (*k).a->GetVector() - (*k).b->GetVector();
        (*k).rab = vab.length();
        //(*k).delta = bond_lengths[idx] - (*k).rab;
        (*k).delta = bond_lengths[idx] - (*k).r0;
        delta2 = (*k).delta * (*k).delta;
 
        energy += 143.9325 * 0.5 * (*k).kb * delta2 * (1.0 - 2.0 * (*k).delta + 7.0/12.0 * 4.0 * delta2);

        idx++;
      }
 

      double delta, err;
      cout << endl;
      cout << "TERM                   OB ENERGY    LOG ENERGY      DELTA     REL.ERR (%)     E-ERR        " << endl;
      cout << "-----------------------------------------------------------------------------------------" << endl;
    
      delta = (E_Bond() - ebond);
      err = delta / ebond * 100;
      sprintf(logbuf, "Bond Stretching        %8.3f      %8.3f     %8.3f     %5.2f   %8.3f", E_Bond(), ebond, delta, err, energy);
      cout << logbuf << endl;
    
      delta = (E_Angle() - eangle);
      err = delta / eangle * 100;
      sprintf(logbuf, "Angle Bending          %8.3f      %8.3f     %8.3f     %5.2f", E_Angle(), eangle, delta, err);
      cout << logbuf << endl;
    
      delta = (E_StrBnd() - estbn);
      err = delta / estbn * 100;
      sprintf(logbuf, "Stretch-Bending        %8.3f      %8.3f     %8.3f     %5.2f", E_StrBnd(), estbn, delta, err);
      cout << logbuf << endl;
    
      delta = (E_OOP() - eoop);
      err = delta / eoop * 100;
      sprintf(logbuf, "Out-Of-Plane Bending   %8.3f      %8.3f     %8.3f     %5.2f", E_OOP(), eoop, delta, err);
      cout << logbuf << endl;
    
      delta = (E_Torsion() - etor);
      err = delta / etor * 100;
      sprintf(logbuf, "Torsional              %8.3f      %8.3f     %8.3f     %5.2f", E_Torsion(), etor, delta, err);
      cout << logbuf << endl;
    
      delta = (E_VDW() - evdw);
      err = delta / evdw * 100;
      sprintf(logbuf, "Van der Waals          %8.3f      %8.3f     %8.3f     %5.2f", E_VDW(), evdw, delta, err);
      cout << logbuf << endl;
      
      delta = (E_Electrostatic() - eeq);
      err = delta / ebond * 100;
      sprintf(logbuf, "Electrostatic          %8.3f      %8.3f     %8.3f     %5.2f", E_Electrostatic(), eeq, delta, err);
      cout << logbuf << endl;

      cout << endl;
      delta = (Energy() - etot);
      err = delta / etot * 100;
      sprintf(logbuf, "Total ENERGY           %8.3f      %8.3f     %8.3f     %5.2f", Energy(), etot, delta, err);
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
        if (((*i).a->GetIdx() == a->GetIdx()) || ((*i).b->GetIdx() == a->GetIdx()))
          grad += i->GetGradient(&*a);

    if ((terms & OBFF_ENERGY) || (terms & OBFF_EANGLE))
      for (i2 = _anglecalculations.begin(); i2 != _anglecalculations.end(); ++i2)
        if (((*i2).a->GetIdx() == a->GetIdx()) || ((*i2).b->GetIdx() == a->GetIdx()) || ((*i2).c->GetIdx() == a->GetIdx()))
          grad += i2->GetGradient(&*a);
    
    if ((terms & OBFF_ENERGY) || (terms & OBFF_ESTRBND))
      for (i3 = _strbndcalculations.begin(); i3 != _strbndcalculations.end(); ++i3)
        if (((*i3).a->GetIdx() == a->GetIdx()) || ((*i3).b->GetIdx() == a->GetIdx()) || ((*i3).c->GetIdx() == a->GetIdx()))
          grad += i3->GetGradient(&*a);
      
    if ((terms & OBFF_ENERGY) || (terms & OBFF_ETORSION))
      for (i4 = _torsioncalculations.begin(); i4 != _torsioncalculations.end(); ++i4)
        if (((*i4).a->GetIdx() == a->GetIdx()) || ((*i4).b->GetIdx() == a->GetIdx()) || ((*i4).c->GetIdx() == a->GetIdx()) || ((*i4).d->GetIdx() == a->GetIdx()))
          grad += i4->GetGradient(&*a);
    
    if ((terms & OBFF_ENERGY) || (terms & OBFF_EOOP))
      for (i5 = _oopcalculations.begin(); i5 != _oopcalculations.end(); ++i5)
        if (((*i5).a->GetIdx() == a->GetIdx()) || ((*i5).b->GetIdx() == a->GetIdx()) || ((*i5).c->GetIdx() == a->GetIdx()) || ((*i5).d->GetIdx() == a->GetIdx()))
          grad += i5->GetGradient(&*a);
      
    if ((terms & OBFF_ENERGY) || (terms & OBFF_EVDW))
      for (i6 = _vdwcalculations.begin(); i6 != _vdwcalculations.end(); ++i6)
        if (((*i6).a->GetIdx() == a->GetIdx()) || ((*i6).b->GetIdx() == a->GetIdx()))
          grad += i6->GetGradient(&*a);
    
    if ((terms & OBFF_ENERGY) || (terms & OBFF_EELECTROSTATIC))
      for (i7 = _electrostaticcalculations.begin(); i7 != _electrostaticcalculations.end(); ++i7)
        if (((*i7).a->GetIdx() == a->GetIdx()) || ((*i7).b->GetIdx() == a->GetIdx()))
          grad += i7->GetGradient(&*a);

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

  int OBForceFieldMMFF94::GetBondType(OBAtom* a, OBAtom* b)
  {
    /*
      cout << "IsInSameRing = " << IsInSameRing(a,b) << endl;
      cout << "a: idx = " << a->GetIdx() << "  type = " << a->GetType() << endl;
      cout << "  SB-MB = " << HasSbmbSet(atoi(a->GetType())) << endl;
      cout << "  Arom = " << HasAromSet(atoi(a->GetType())) << endl;
      cout << "  IsAromatic = " << a->IsAromatic() << endl;
      cout << "b: idx = " << b->GetIdx() << "  type = " << b->GetType() << endl;
      cout << "  SB-MB = " << HasSbmbSet(atoi(b->GetType())) << endl;
      cout << "  Arom = " << HasAromSet(atoi(b->GetType())) << endl;
      cout << "  IsAromatic = " << b->IsAromatic() << endl;
    */
    if (_mol.GetBond(a,b)->IsSingle()) {
      if (HasSbmbSet(atoi(a->GetType())) && HasSbmbSet(atoi(b->GetType())) && !HasAromSet(atoi(a->GetType())) && !HasAromSet(atoi(b->GetType())))
        return 1;
        
      //	if (!IsInSameRing(a, b) && HasAromSet(atoi(a->GetType())) && HasAromSet(atoi(b->GetType())))
      //int ia = a->BOSum() - a->GetValence();
      //int ib = b->BOSum() - b->GetValence();
      bool ia, ib;
      ia = ib = false;
      FOR_NBORS_OF_ATOM(nbr, a) {
        OBBond* bond = _mol.GetBond(a, &*nbr);
        if (bond->IsDouble() && !IsInSameRing(a, &*nbr))
          ia = true;
      }
      FOR_NBORS_OF_ATOM(nbr, b) {
        OBBond* bond = _mol.GetBond(b, &*nbr);
        if (bond->IsDouble() && !IsInSameRing(b, &*nbr))
          ib = true;
      }

      if (!IsInSameRing(a, b) && HasSbmbSet(atoi(a->GetType())) && HasSbmbSet(atoi(b->GetType())) && ia && ib)
        return 1;

      if (HasSbmbSet(atoi(a->GetType())) && HasSbmbSet(atoi(b->GetType()))) {
        if ((atoi(a->GetType()) == 80) || (atoi(b->GetType()) == 80))
          return 0;

        if (a->IsAromatic() && b->IsAromatic() &&  
            !HasAromSet(atoi(a->GetType())) && !HasAromSet(atoi(b->GetType())) &&
            IsInSameRing(a,b))
          return 1; // phenone
        //if (a->IsAromatic() && !b->IsAromatic() &&  
        if (a->IsAromatic() && HasAromSet(atoi(a->GetType())) && !HasAromSet(atoi(b->GetType())))
          // && !b->IsInRing())
          return 1; // benzoic acid
        //if (!a->IsAromatic() && b->IsAromatic() &&  
        if (b->IsAromatic() && !HasAromSet(atoi(a->GetType())) && HasAromSet(atoi(b->GetType())))
          // && !a->IsInRing())
          return 1; // benzoic acid
      }
    }

    return 0;
  }
  
  int OBForceFieldMMFF94::GetAngleType(OBAtom* a, OBAtom* b, OBAtom *c)
  {
    int sumbondtypes;

    sumbondtypes = GetBondType(a,b) + GetBondType(b, c);

    if (a->IsInRingSize(3) && b->IsInRingSize(3) && c->IsInRingSize(3)) {
      switch (sumbondtypes) {
      case 0:
        return 3; 
      case 1:
        return 5; 
      case 2:
        return 6; 
      }
    }
    
    if (a->IsInRingSize(4) && b->IsInRingSize(4) && c->IsInRingSize(4)) {
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

    if (a->IsInRingSize(5) && b->IsInRingSize(5) && c->IsInRingSize(5) && d->IsInRingSize(5))
      if (IsInSameRing(a,b) && IsInSameRing(b,c) && IsInSameRing(c,d) && !a->IsAromatic()) {
        OBBond *bond;
        bond = _mol.GetBond(a, b);
        if (!bond->IsSingle())
          return 0;
        bond = _mol.GetBond(b, c);
        if (!bond->IsSingle())
          return 0;
        bond = _mol.GetBond(c, d);
        if (!bond->IsSingle())
          return 0;
        
        return 5; // needs checking for saturations
      }

    return 0;
  }

  // used for bond typing
  // note: when two atoms are in the same ring, but one of them
  // is also in another ring, we return false
  /*
    bool OBForceFieldMMFF94::IsInSameRing(OBAtom* a, OBAtom* b)
    {
    int a_in, b_in, a_in_cur, b_in_cur;
    vector<OBRing*> vr;
    vr = _mol.GetSSSR();
    a_in = 0;
    b_in = 0;
    
    vector<OBRing*>::iterator i;
    vector<int>::iterator j;
    for (i = vr.begin();i != vr.end();i++) {
    for(j = (*i)->_path.begin();j != (*i)->_path.end();j++) {
    if (*j == a->GetIdx())
    a_in++;
    if (*j == b->GetIdx())
    b_in++;
    }
    }

    for (i = vr.begin();i != vr.end();i++) {
    a_in_cur = 0;
    b_in_cur = 0;
    for(j = (*i)->_path.begin();j != (*i)->_path.end();j++) {
    if (*j == a->GetIdx())
    a_in_cur++;
    if (*j == b->GetIdx())
    b_in_cur++;
    }
      
    if (a_in_cur && b_in_cur) {
    if (((a_in >= 2) && (b_in == 1)) || ((b_in >= 2) && (a_in == 1))) // see note above
    return false;

    return true;
    }
    }

    return false;
    }
  */
  bool OBForceFieldMMFF94::HasLinSet(int atomtype)
  {
    vector<int>::iterator i;

    for (i = _lin.begin(); i != _lin.end(); i++)
      if ((*i) == atomtype)
        return true;

    return false;
  }
 
  
  bool OBForceFieldMMFF94::HasAromSet(int atomtype)
  {
    vector<int>::iterator i;

    for (i = _arom.begin(); i != _arom.end(); i++)
      if (*i == atomtype)
        return true;

    return false;
  }
 
  bool OBForceFieldMMFF94::HasSbmbSet(int atomtype)
  {
    vector<int>::iterator i;

    for (i = _sbmb.begin(); i != _sbmb.end(); i++)
      if ((*i) == atomtype)
        return true;

    return false;
  }

  OBFFParameter* OBForceFieldMMFF94::GetParameterMMFF94(int ffclass, int a, int b, int c, int d, std::vector<OBFFParameter> &parameter)
  {
    OBFFParameter *par;

    if (!b)
      for (int idx=0; idx < parameter.size(); idx++)
        if ((a == parameter[idx].a) && (ffclass == parameter[idx].ipar5)) {
          par = &parameter[idx];
          return par;
        }

    if (!c)
      for (int idx=0; idx < parameter.size(); idx++)
        if (((a == parameter[idx].a) && (b == parameter[idx].b) && (ffclass == parameter[idx].ipar5)) || 
            ((a == parameter[idx].b) && (b == parameter[idx].a) && (ffclass == parameter[idx].ipar5))) 
          {
            par = &parameter[idx];
            return par;
          }

    if (!d)
      for (int idx=0; idx < parameter.size(); idx++)
        if (((a == parameter[idx].a) && (b == parameter[idx].b) && (c == parameter[idx].c) && (ffclass == parameter[idx].ipar5)) || 
            ((a == parameter[idx].c) && (b == parameter[idx].b) && (c == parameter[idx].a) && (ffclass == parameter[idx].ipar5))) 
          {
            par = &parameter[idx];
            return par;
          }

    for (int idx=0; idx < parameter.size(); idx++)
      if (((a == parameter[idx].a) && (b == parameter[idx].b) && (c == parameter[idx].c) && (d == parameter[idx].d) && (ffclass == parameter[idx].ipar5)) || 
          ((a == parameter[idx].d) && (b == parameter[idx].c) && (c == parameter[idx].b) && (d == parameter[idx].a) && (ffclass == parameter[idx].ipar5))) 
        {
          par = &parameter[idx];
          return par;
        }

    return NULL;
  }
 
} // end namespace OpenBabel

//! \file forcefieldmmff94.cpp
//! \brief MMFF94 force field

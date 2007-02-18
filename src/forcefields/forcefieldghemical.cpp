/**********************************************************************
forcefieldghemical.cpp - Ghemical force field.
 
Copyright (C) 2006 by Tim Vandermeersch <tim.vandermeersch@gmail.com>
 
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
#include <openbabel/mol.h>
#include "forcefieldghemical.h"

using namespace std;

namespace OpenBabel
{
  double OBFFBondCalculationGhemical::GetEnergy()
  {
    double delta2;
    vector3 vab;
     
    vab = a->GetVector() - b->GetVector();
    rab = vab.length();
    delta = rab - r0;
    delta2 = delta * delta;

    energy = KCAL_TO_KJ * kb * delta2;

    return energy;
  }
   
  vector3 OBFFBondCalculationGhemical::GetGradient(OBAtom *atom) 
  {
    vector3 vab, gradient;

    if ((atom != a) && (atom != b))
      return  vector3(0,0,0);
     
    vab = a->GetVector() - b->GetVector();
    rab = vab.length();
    delta = rab - r0;
 
    gradient = KCAL_TO_KJ * 2.0f * kb * delta * (vab / rab);

    if (atom == a)
      return -gradient;
    else 
      return gradient;
  }

  double OBForceFieldGhemical::E_Bond()
  {
    vector<OBFFBondCalculationGhemical>::iterator i;
    double energy;
    char logbuf[150];
    
    energy = 0.0f;
    
    IF_OBFF_LOGLVL_HIGH {
      *logos << endl << "B O N D   S T R E T C H I N G" << endl << endl;
      *logos << "ATOM TYPES  BOND    BOND       IDEAL       FORCE" << endl;
      *logos << " I    J     TYPE   LENGTH     LENGTH     CONSTANT      DELTA      ENERGY" << endl;
      *logos << "------------------------------------------------------------------------" << endl;
    }
 
    for (i = _bondcalculations.begin(); i != _bondcalculations.end(); i++) {

      energy += i->GetEnergy();

      IF_OBFF_LOGLVL_HIGH {
        sprintf(logbuf, "%s %s    %d   %8.3f   %8.3f     %8.3f   %8.3f   %8.3f", (*i).a->GetType(), (*i).b->GetType(), 
                (*i).bt, (*i).rab, (*i).r0, (*i).kb, (*i).delta, (*i).energy);
        *logos  << logbuf << endl;
      }
    }
    
    IF_OBFF_LOGLVL_MEDIUM
      *logos << "     TOTAL BOND STRETCHING ENERGY = " << energy << endl;
    return energy;
  }
  
  double OBFFAngleCalculationGhemical::GetEnergy()
  {
    double delta2;

    theta = a->GetAngle(b->GetIdx(), c->GetIdx());

    if (theta > 170) {
      delta2 = 1.0f + cos(theta);

      energy = KCAL_TO_KJ * ka * delta2;
    } else {
      delta = theta - theta0;
      delta2 = delta * delta;
    
      energy = KCAL_TO_KJ * ka * delta2;
    }
     
    return energy;
  }
  
  vector3 OBFFAngleCalculationGhemical::GetGradient(OBAtom *atom) 
  {
    vector3 vab, vcb, grada, gradb, gradc;
    double f, rab, rab2, rcb, rcb2, abcb, abcb2;


    if ((atom->GetIdx() != a->GetIdx()) && (atom->GetIdx() != b->GetIdx()) && (atom->GetIdx() != c->GetIdx()))
      return  vector3(0,0,0);
     
    vab = a->GetVector() - b->GetVector();
    vcb = c->GetVector() - b->GetVector();
    rab = vab.length();
    rcb = vcb.length();
    rab2 = rab * rab;
    rcb2 = rcb * rcb;
    
    delta = theta - theta0;
    abcb = dot(vab, vcb) / (rab * rcb);
    abcb2 = 1.0f - abcb * abcb;
    theta = acos(abcb) * RAD_TO_DEG;
    //theta = a->GetAngle(b->GetIdx(), c->GetIdx());

    f = KCAL_TO_KJ * 2.0f * ka * delta;
    
    if (IsNearZero(abcb2))
      return vector3(0.0f, 0.0f, 0.0f);

    grada = f * (vcb * rab * rcb - (vab / rab) * dot(vab, vcb) * rcb) / (sqrt(abcb2) * rab2 * rcb2);
    gradc = -f * ((vcb / rcb) * dot(vab, vcb) * rab - vab * rab * rcb) / (sqrt(abcb2) * rab2 * rcb2);
    gradb = -grada - gradc;

    grada *= (1.0f / DEG_TO_RAD);
    gradb *= (1.0f / DEG_TO_RAD);
    gradc *= (1.0f / DEG_TO_RAD);

    if (atom->GetIdx() == a->GetIdx())
      return grada;
    else if (atom->GetIdx() == c->GetIdx())
      return gradc;
    else
      return gradb;
  }
 
  double OBForceFieldGhemical::E_Angle()
  {
    vector<OBFFAngleCalculationGhemical>::iterator i;
    double energy;
    char logbuf[150];
    
    energy = 0.0f;
    
    IF_OBFF_LOGLVL_HIGH {
      *logos << endl << "A N G L E   B E N D I N G" << endl << endl;
      *logos << "ATOM TYPES       VALENCE     IDEAL      FORCE" << endl;
      *logos << " I    J    K      ANGLE      ANGLE     CONSTANT      DELTA      ENERGY" << endl;
      *logos << "-----------------------------------------------------------------------------" << endl;
    }
    
    for (i = _anglecalculations.begin(); i != _anglecalculations.end(); i++) {

      energy += i->GetEnergy();
      
      IF_OBFF_LOGLVL_HIGH {
        sprintf(logbuf, "%s %s %s  %8.3f   %8.3f     %8.3f   %8.3f   %8.3f", (*i).a->GetType(), (*i).b->GetType(), 
                (*i).c->GetType(), (*i).theta, (*i).theta0, (*i).ka, (*i).delta, (*i).energy);
        *logos  << logbuf << endl;
      }
    }
 
    IF_OBFF_LOGLVL_MEDIUM
      *logos << "     TOTAL ANGLE BENDING ENERGY = " << energy << endl;
    return energy;
  }
  
  double OBFFTorsionCalculationGhemical::GetEnergy()
  {
    double cosine, cosine2, cosine3;
    double phi1, phi2, phi3;
    vector3 vab, vbc, vcd, abbc, bccd;

    //tor = CalcTorsionAngle(a->GetVector(), b->GetVector(), c->GetVector(), d->GetVector());
    vab = a->GetVector() - b->GetVector();
    vbc = b->GetVector() - c->GetVector();
    vcd = c->GetVector() - d->GetVector();
    abbc = cross(vab, vbc);
    bccd = cross(vbc, vcd);
    tor = RAD_TO_DEG * acos(dot(abbc, bccd) / (abbc.length() * bccd.length()));
    if (dot(abbc, bccd) > 0.0f)
      tor = -tor;
    
    cosine = cos(DEG_TO_RAD * (n * tor));
    energy = KCAL_TO_KJ * V * (1 + s * cosine);
    
    return energy;
  }
  
  vector3 OBFFTorsionCalculationGhemical::GetGradient(OBAtom *atom) 
  {
    vector3 vab, vbc, vcd, vac, vbd, grada, gradb, gradc, gradd;
    vector3 abbc, bccd, bcabbc, bcbccd, cdabbc, cdbccd, acabbc, acbccd, ababbc, abbccd, bdabbc, bdbccd;
    double f, rab, rbc, rcd, rabbc, rbccd, rabbc2, rbccd2, rabbc3, rbccd3, abbc_bccd, abbc_bccd2;
    double sine, sine2, sine3;


    if ((atom->GetIdx() != a->GetIdx()) && (atom->GetIdx() != b->GetIdx()) && (atom->GetIdx() != c->GetIdx()))
      return  vector3(0,0,0);
     
    vab = a->GetVector() - b->GetVector();
    vbc = b->GetVector() - c->GetVector();
    vcd = c->GetVector() - d->GetVector();
    vac = a->GetVector() - c->GetVector();
    vbd = b->GetVector() - d->GetVector();
    abbc = cross(vab, vbc);
    bccd = cross(vbc, vcd);
    
    //tor = CalcTorsionAngle(a->GetVector(), b->GetVector(), c->GetVector(), d->GetVector());
    tor = RAD_TO_DEG * acos(dot(abbc, bccd) / (abbc.length() * bccd.length()));
    if (dot(abbc, bccd) > 0.0f)
      tor = -tor;
 
    bcabbc = cross(vbc, abbc);
    bcbccd = cross(vbc, bccd);
    cdabbc = cross(vcd, abbc);
    cdbccd = cross(vcd, bccd);
    acabbc = cross(vac, abbc);
    acbccd = cross(vac, bccd);
    ababbc = cross(vab, abbc);
    abbccd = cross(vab, bccd);
    bdabbc = cross(vbd, abbc);
    bdbccd = cross(vbd, bccd);
    rabbc = abbc.length();
    rbccd = bccd.length();
    rabbc2 = rabbc * rabbc;
    rbccd2 = rbccd * rbccd;
    rabbc3 = rabbc2 * rabbc;
    rbccd3 = rbccd2 * rbccd;
    abbc_bccd = dot(abbc, bccd) / (rabbc * rbccd);
    abbc_bccd2 = 1.0f - abbc_bccd * abbc_bccd;
    
    sine = sin(DEG_TO_RAD * (n * tor));
    f = - KCAL_TO_KJ * V * s * n * sine;
    if (dot(abbc, bccd) > 0.0f)
      f = -f;
 
    if (IsNearZero(f))
      return vector3(0.0f, 0.0f, 0.0f);

    grada = f * (bcbccd / (rabbc*rbccd) - (bcabbc*dot(abbc,bccd)) / (rabbc3*rbccd)) / sqrt(abbc_bccd2);
    gradd = f * (bcabbc / (rabbc*rbccd) - (bcbccd*dot(abbc,bccd)) / (rabbc*rbccd3)) / sqrt(abbc_bccd2);

    gradb = f * ( -(cdbccd*dot(abbc,bccd)) / (rabbc*rbccd3) + 
                  (cdabbc - acbccd) / (rabbc*rbccd) + 
                  (acabbc*dot(abbc,bccd)) / (rabbc3*rbccd)    ) / sqrt(abbc_bccd2);
    
    gradc = f * (  (bdbccd*dot(abbc,bccd)) / (rabbc*rbccd3) + 
                   (abbccd - bdabbc) / (rabbc*rbccd) +
                   -(ababbc*dot(abbc,bccd)) / (rabbc3*rbccd)    ) / sqrt(abbc_bccd2);
    
    if (atom->GetIdx() == a->GetIdx()) 
      return grada;
    else if (atom->GetIdx() == d->GetIdx())
      return gradd;
    else if (atom->GetIdx() == b->GetIdx())
      return gradb;
    else
      return gradc;
  }
 

  double OBForceFieldGhemical::E_Torsion() 
  {
    vector<OBFFTorsionCalculationGhemical>::iterator i;
    double energy;
    char logbuf[150];
 
    energy = 0.0f;

    IF_OBFF_LOGLVL_HIGH {
      *logos << endl << "T O R S I O N A L" << endl << endl;
      *logos << "----ATOM TYPES-----    FORCE              TORSION" << endl;
      *logos << " I    J    K    L     CONSTANT     s       ANGLE    n    ENERGY" << endl;
      *logos << "----------------------------------------------------------------" << endl;
    }
    
    for (i = _torsioncalculations.begin(); i != _torsioncalculations.end(); i++) {

      energy += i->GetEnergy();
      
      IF_OBFF_LOGLVL_HIGH {
        sprintf(logbuf, "%s %s %s %s    %6.3f    %5.0f   %8.3f   %1.0f   %8.3f", (*i).a->GetType(), (*i).b->GetType(), 
                (*i).c->GetType(), (*i).d->GetType(), (*i).V, (*i).s, (*i).tor, (*i).n, (*i).energy);
        *logos << logbuf << endl;
      }
    }

    IF_OBFF_LOGLVL_MEDIUM
      *logos <<  "     TOTAL TORSIONAL ENERGY = " << energy << endl;
    return energy;
  }

  double OBFFVDWCalculationGhemical::GetEnergy()
  {
    double sigma2, sigma4;

    rab = a->GetDistance(b);

    sigma = rab / (Ra + Rb);
    sigma2 = sigma * sigma;
    sigma4 = sigma2 * sigma2;
    sigma6 = sigma4 * sigma2;
    sigma12 = sigma6 * sigma6;
    energy = KCAL_TO_KJ * kab * ((1.0f / sigma12) - (2.0f / sigma6));

    return energy;
  }
  
  vector3 OBFFVDWCalculationGhemical::GetGradient(OBAtom *atom) 
  {
    vector3 vab, gradient;
    double sigma2, sigma4;

    if ((atom != a) && (atom != b))
      return  vector3(0,0,0);
     
    vab = a->GetVector() - b->GetVector();
    rab = vab.length();
    
    sigma = rab / (Ra + Rb);
    sigma2 = sigma * sigma;
    sigma4 = sigma2 * sigma2;
    sigma6 = sigma4 * sigma2;
    sigma12 = sigma6 * sigma6;
 
    gradient = KCAL_TO_KJ * 12.0f * (kab / rab) * ((1.0f / sigma6) - (1.0f / sigma12)) * (vab / rab);
    
    if (atom == a)
      return -gradient;
    else 
      return gradient;
  }


  double OBForceFieldGhemical::E_VDW()
  {
    vector<OBFFVDWCalculationGhemical>::iterator i;
    double energy;
    char logbuf[150];
 
    energy = 0.0f;
    
    IF_OBFF_LOGLVL_HIGH {
      *logos << endl << "V A N   D E R   W A A L S" << endl << endl;
      *logos << "ATOM TYPES          " << endl;
      *logos << " I    J        Rij       kij      SIGMA     ENERGY" << endl;
      *logos << "---------------------------------------------------" << endl;
      //          XX   XX     -000.000  -000.000  -000.000  -000.000
    }
    
    for (i = _vdwcalculations.begin(); i != _vdwcalculations.end(); i++) {
      
      energy += i->GetEnergy();
      
      IF_OBFF_LOGLVL_HIGH {
        sprintf(logbuf, "%s %s   %8.3f  %8.3f  %8.3f  %8.3f", (*i).a->GetType(), (*i).b->GetType(), 
                (*i).rab, (*i).kab, (*i).sigma, (*i).energy);
        *logos  << logbuf << endl;
      }
    }

    IF_OBFF_LOGLVL_MEDIUM
      *logos  << "     TOTAL VAN DER WAALS ENERGY = " << energy << endl;
    return energy;
  }

  double OBFFElectrostaticCalculationGhemical::GetEnergy()
  {
    vector3 vab;

    vab = a->GetVector() - b->GetVector();
    rab = vab.length();

    energy = KCAL_TO_KJ * qq / rab;

    return energy;
  }
  
  vector3 OBFFElectrostaticCalculationGhemical::GetGradient(OBAtom *atom) 
  {
    vector3 vab, gradient;
    double rab2, rab3;

    if ((atom != a) && (atom != b))
      return  vector3(0,0,0);
     
    vab = a->GetVector() - b->GetVector();
    rab = vab.length();
    rab2 = rab * rab;
    rab3 = rab2 * rab;
    
    gradient = KCAL_TO_KJ * qq * (vab / rab3);
    
    if (atom == a)
      return gradient;
    else 
      return -gradient;
  }
  
  double OBForceFieldGhemical::E_Electrostatic()
  {
    vector<OBFFElectrostaticCalculationGhemical>::iterator i;
    double energy;
    char logbuf[150];
 
    energy = 0.0f;
    
    IF_OBFF_LOGLVL_HIGH {
      *logos << std::endl << "E L E C T R O S T A T I C   I N T E R A C T I O N S" << std::endl << std::endl;
      *logos << "ATOM TYPES          " << std::endl;
      *logos << " I    J           Rij   332.17*QiQj  ENERGY" << std::endl;
      *logos << "-------------------------------------------" << std::endl;
      //            XX   XX     -000.000  -000.000  -000.000  
    }

    for (i = _electrostaticcalculations.begin(); i != _electrostaticcalculations.end(); i++) {
      
      energy += i->GetEnergy();
      
      IF_OBFF_LOGLVL_HIGH {
        sprintf(logbuf, "%s %s   %8.3f  %8.3f  %8.3f", (*i).a->GetType(), (*i).b->GetType(), 
                (*i).rab, (*i).qq, (*i).energy);
        *logos  << logbuf << endl;
      }
    }

    IF_OBFF_LOGLVL_MEDIUM
      *logos << "     TOTAL ELECTROSTATIC ENERGY = " << energy << endl << endl;
    return energy;
  }

  //***********************************************
  //Make a global instance
  OBForceFieldGhemical theForceFieldGhemical("Ghemical",false);
  //***********************************************

  OBForceFieldGhemical::~OBForceFieldGhemical()
  {
  }

  OBForceFieldGhemical &OBForceFieldGhemical::operator=(OBForceFieldGhemical &src)
  {
    _mol = src._mol;

    _ffbondparams    = src._ffbondparams;
    _ffangleparams   = src._ffangleparams;
    _fftorsionparams = src._fftorsionparams;
    _ffvdwparams     = src._ffvdwparams;

    _bondcalculations          = src._bondcalculations;
    _anglecalculations         = src._anglecalculations;
    _torsioncalculations       = src._torsioncalculations;
    _vdwcalculations           = src._vdwcalculations;
    _electrostaticcalculations = src._electrostaticcalculations;

    return *this;
  }

  bool OBForceFieldGhemical::Setup(OBMol &mol)
  {
    _mol = mol;
    SetGhemicalTypes();
    
    if (!SetupCalculations())
      return false;
    
    return true;
  }
  
  bool OBForceFieldGhemical::SetupCalculations()
  {
    OBFFParameter *parameter;
    OBAtom *a, *b, *c, *d;
    bool found;
    
    // 
    // Bond Calculations
    //
    OBFFBondCalculationGhemical bondcalc;
    int bondtype;

    _bondcalculations.clear();
    
    FOR_BONDS_OF_MOL(bond, _mol) {
      a = bond->GetBeginAtom();
      b = bond->GetEndAtom();	
      bondtype = bond->GetBondOrder(); 
      if (bond->IsAromatic())
        bondtype = 5;

      bondcalc.a = a;
      bondcalc.b = b;
      bondcalc.bt = bondtype;

      parameter = GetParameterGhemical(bondtype, a->GetType(), b->GetType(), NULL, NULL,  _ffbondparams);
      if (parameter == NULL) {
        parameter = GetParameterGhemical(bondtype, "FFFF", a->GetType(), NULL, NULL, _ffbondparams);

        if (parameter == NULL) {
          parameter = GetParameterGhemical(bondtype, "FFFF", b->GetType(), NULL, NULL, _ffbondparams);

          if (parameter == NULL) {
            //obErrorLog.ThrowError(__FUNCTION__, "Could not find all bond parameters ", obError);
            //return false;
            bondcalc.kb = 500.0f;
            bondcalc.r0 = 1.100f;

            _bondcalculations.push_back(bondcalc);
            continue;
          }
        }
      }
      bondcalc.kb = parameter->dpar2;
      bondcalc.r0 = parameter->dpar1;

      _bondcalculations.push_back(bondcalc);
    }

    //
    // Angle Calculations
    //
    OBFFAngleCalculationGhemical anglecalc;
 
    _anglecalculations.clear();
    
    FOR_ANGLES_OF_MOL(angle, _mol) {
      b = _mol.GetAtom((*angle)[0] + 1);
      a = _mol.GetAtom((*angle)[1] + 1);
      c = _mol.GetAtom((*angle)[2] + 1);
      
      anglecalc.a = a;
      anglecalc.b = b;
      anglecalc.c = c;

      parameter = GetParameter(a->GetType(), b->GetType(), c->GetType(), _ffangleparams);
      if (parameter == NULL) {
        parameter = GetParameter("FFFF", b->GetType(), c->GetType(), _ffangleparams);
        if (parameter == NULL) {
          parameter = GetParameter("FFFF", b->GetType(), "FFFF", _ffangleparams);
          if (parameter == NULL) {
            anglecalc.ka   = 0.020f;
            anglecalc.theta0 = 120.0f;
            
            _anglecalculations.push_back(anglecalc);
            continue;
            //obErrorLog.ThrowError(__FUNCTION__, "Could not find all angle parameters", obError);
            //return false;
          }
        }
      }
      anglecalc.ka   = parameter->dpar2;
      anglecalc.theta0 = parameter->dpar1;
      
      _anglecalculations.push_back(anglecalc);
    }
    
    //
    // Torsion Calculations
    //
    OBFFTorsionCalculationGhemical torsioncalc;
    int torsiontype;

    _torsioncalculations.clear();
 
    FOR_TORSIONS_OF_MOL(t, _mol) {
      a = _mol.GetAtom((*t)[0] + 1);
      b = _mol.GetAtom((*t)[1] + 1);
      c = _mol.GetAtom((*t)[2] + 1);
      d = _mol.GetAtom((*t)[3] + 1);
      OBBond *bc = _mol.GetBond(b, c);
      torsiontype = bc->GetBondOrder(); 
      
      torsioncalc.a = a;
      torsioncalc.b = b;
      torsioncalc.c = c;
      torsioncalc.d = d;
      torsioncalc.tt = torsiontype;

      parameter = GetParameterGhemical(torsiontype, a->GetType(), b->GetType(), c->GetType(), d->GetType(), _fftorsionparams);
      if (parameter == NULL) {
        parameter = GetParameterGhemical(torsiontype, "FFFF", b->GetType(), c->GetType(), d->GetType(), _fftorsionparams);
        if (parameter == NULL) {
          parameter = GetParameterGhemical(torsiontype, a->GetType(), b->GetType(), c->GetType(), "FFFF", _fftorsionparams);
          if (parameter == NULL) {
            parameter = GetParameterGhemical(torsiontype, "FFFF", b->GetType(), c->GetType(), "FFFF", _fftorsionparams);
            if (parameter == NULL) {
              torsioncalc.V = 0.0f;
              torsioncalc.s = 1.0f;
              torsioncalc.n = 1.0f;
              _torsioncalculations.push_back(torsioncalc);
              continue;
              //obErrorLog.ThrowError(__FUNCTION__, "Could not find all torsion parameters", obError);
              //return false;
            }
          }
        }
      }
      torsioncalc.V = parameter->dpar1;
      torsioncalc.s = parameter->dpar2;
      torsioncalc.n = parameter->dpar3;
      _torsioncalculations.push_back(torsioncalc);     
    }
    
    // 
    // VDW Calculations
    //
    OBFFVDWCalculationGhemical vdwcalc;
    OBFFParameter *parameter_a, *parameter_b;

    _vdwcalculations.clear();
    
    FOR_PAIRS_OF_MOL(p, _mol) {
      a = _mol.GetAtom((*p)[0]);
      b = _mol.GetAtom((*p)[1]);

      parameter_a = GetParameter(a->GetType(), _ffvdwparams);
      if (parameter_a == NULL) { // no vdw parameter -> use hydrogen
        vdwcalc.Ra = 1.5f;
        vdwcalc.ka = 0.042f;
      } else {
        vdwcalc.Ra = parameter_a->dpar1;
        vdwcalc.ka = parameter_a->dpar2;
      }

      parameter_b = GetParameter(b->GetType(), _ffvdwparams);
      if (parameter_b == NULL) { // no vdw parameter -> use hydrogen
        vdwcalc.Rb = 1.5f;
        vdwcalc.kb = 0.042;;
      } else {
        vdwcalc.Rb = parameter_b->dpar1;
        vdwcalc.kb = parameter_b->dpar2;
      }

      vdwcalc.a = &*a;
      vdwcalc.b = &*b;
     
      //this calculations only need to be done once for each pair, 
      //we do them now and save them for later use
      vdwcalc.kab = sqrt(vdwcalc.ka * vdwcalc.kb);
      
      // 1-4 scaling
      FOR_NBORS_OF_ATOM (nbr, a)
        FOR_NBORS_OF_ATOM (nbr2, &*nbr)
          FOR_NBORS_OF_ATOM (nbr3, &*nbr2)
            if (b == &*nbr3)
              vdwcalc.kab *= 0.5f;

      _vdwcalculations.push_back(vdwcalc);
    }
    
    // 
    // Electrostatic Calculations
    //
    OBFFElectrostaticCalculationGhemical elecalc;

    _electrostaticcalculations.clear();
    
    FOR_PAIRS_OF_MOL(p, _mol) {
      a = _mol.GetAtom((*p)[0]);
      b = _mol.GetAtom((*p)[1]);
      
      elecalc.qq = 332.17f * a->GetPartialCharge() * b->GetPartialCharge();
      
      if (elecalc.qq) {
        elecalc.a = &*a;
        elecalc.b = &*b;
        
	// 1-4 scaling
        FOR_NBORS_OF_ATOM (nbr, a)
          FOR_NBORS_OF_ATOM (nbr2, &*nbr)
            FOR_NBORS_OF_ATOM (nbr3, &*nbr2)
              if (b == &*nbr3)
                elecalc.qq *= 0.5f;

        _electrostaticcalculations.push_back(elecalc);
      }
    }

    return true;
  }

  bool OBForceFieldGhemical::ParseParamFile()
  {
    vector<string> vs;
    char buffer[80];
    
    OBFFParameter parameter;

    // open data/ghemical.prm
    string buffer2, subbuffer;
    ifstream ifs1, ifs2, *ifsP;
    buffer2 = BABEL_DATADIR;
    buffer2 += FILE_SEP_CHAR;
    subbuffer = buffer2;
    subbuffer += BABEL_VERSION;
    subbuffer += FILE_SEP_CHAR;
    subbuffer += "ghemical.prm";
    buffer2 += "ghemical.prm";

    ifs1.open(subbuffer.c_str());
    ifsP= &ifs1;
    if (!(*ifsP))
      {
        ifs2.open(buffer2.c_str());
        ifsP = &ifs2;
      }

    while (ifsP->getline(buffer, 80)) {
      tokenize(vs, buffer);

      if (EQn(buffer, "bond", 4)) {
        parameter.clear();
        parameter._a = vs[1];
        parameter._b = vs[2];
        parameter.dpar1 = atof(vs[4].c_str()); // length
        parameter.dpar2 = atof(vs[5].c_str()); // force cte
        if (EQn(vs[3].c_str(), "S", 1))
          parameter.ipar5 = 1;
        if (EQn(vs[3].c_str(), "D", 1))
          parameter.ipar5 = 2;
        if (EQn(vs[3].c_str(), "T", 1))
          parameter.ipar5 = 3;
        if (EQn(vs[3].c_str(), "C", 1))
          parameter.ipar5 = 5;
        _ffbondparams.push_back(parameter);
      }
      if (EQn(buffer, "angle", 5)) {
        parameter.clear();
        parameter._a = vs[1];
        parameter._b = vs[2];
        parameter._c = vs[3];
        parameter.dpar1 = atof(vs[5].c_str()); // angle
        parameter.dpar2 = atof(vs[6].c_str()); // force cte
        _ffangleparams.push_back(parameter);
      }
      if (EQn(buffer, "torsion", 7)) {
        parameter.clear();
        parameter._a = vs[1];
        parameter._b = vs[2];
        parameter._c = vs[3];
        parameter._d = vs[4];
        parameter.dpar1 = atof(vs[6].c_str()); // force cte
        parameter.dpar3 = atof(vs[8].c_str()); // n
        if (EQn(vs[7].c_str(), "+", 1))
          parameter.dpar2 = +1; // s
        if (EQn(vs[7].c_str(), "-", 1))
          parameter.dpar2 = -1; // s
        if (EQn(vs[5].c_str(), "?S?", 3))
          parameter.ipar5 = 1;
        if (EQn(vs[5].c_str(), "?D?", 3))
          parameter.ipar5 = 2;
        if (EQn(vs[5].c_str(), "?T?", 3))
          parameter.ipar5 = 3;
        if (EQn(vs[5].c_str(), "?C?", 3)) {
          parameter.ipar5 = 1;
          _fftorsionparams.push_back(parameter);
          parameter.ipar5 = 2;
        }
        _fftorsionparams.push_back(parameter);
      }
      if (EQn(buffer, "vdw", 3)) {
        parameter.clear();
        parameter._a = vs[1];
        parameter.dpar1 = atof(vs[2].c_str()); // r
        parameter.dpar2 = atof(vs[3].c_str()); // force cte
        _ffvdwparams.push_back(parameter);
      }
      if (EQn(buffer, "charge", 6)) {
        parameter.clear();
        parameter._a = vs[1];
        parameter._b = vs[2];
        if (EQn(vs[3].c_str(), "S", 1))
          parameter.ipar5 = 1;
        if (EQn(vs[3].c_str(), "D", 1))
          parameter.ipar5 = 2;
        parameter.dpar1 = atof(vs[4].c_str()); // charge
        _ffchargeparams.push_back(parameter);
      }
    }
	
    if (ifs1)
      ifs1.close();
    if (ifs2)
      ifs2.close();
 
    return 0;
  }
  
  bool OBForceFieldGhemical::SetGhemicalTypes()
  {
    std::vector<std::vector<int> > _mlist; //!< match list for atom typing
    std::vector<std::pair<OBSmartsPattern*,std::string> > _vexttyp; //!< external atom type rules
    vector<vector<int> >::iterator j;
    vector<pair<OBSmartsPattern*,string> >::iterator i;
    OBSmartsPattern *sp;
    vector<string> vs;
    char buffer[80];
 
    _mol.SetAtomTypesPerceived();
    
    // open data/ghemical.prm
    string buffer2, subbuffer;
    ifstream ifs1, ifs2, *ifsP;
    buffer2 = BABEL_DATADIR;
    buffer2 += FILE_SEP_CHAR;
    subbuffer = buffer2;
    subbuffer += BABEL_VERSION;
    subbuffer += FILE_SEP_CHAR;
    subbuffer += "ghemical.prm";
    buffer2 += "ghemical.prm";

    ifs1.open(subbuffer.c_str());
    ifsP= &ifs1;
    if (!(*ifsP))
      {
        ifs2.open(buffer2.c_str());
        ifsP = &ifs2;
      }

    *logos  << std::endl << "A T O M   T Y P E S" << std::endl << std::endl;
    
    while (ifsP->getline(buffer, 80)) {
      if (EQn(buffer, "atom", 4)) {
      	tokenize(vs, buffer);

        sp = new OBSmartsPattern;
        if (sp->Init(vs[1]))
          _vexttyp.push_back(pair<OBSmartsPattern*,string> (sp,vs[2]));
        else {
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
    }

    SetGhemicalCharges();

    *logos << "IDX\tTYPE\tCHARGE" << std::endl;
    FOR_ATOMS_OF_MOL (a, _mol)
      *logos << a->GetIdx() << "\t" << a->GetType() << "\t" << a->GetPartialCharge() << std::endl;
 
    if (ifs1)
      ifs1.close();
    if (ifs2)
      ifs2.close();

    return true;
  }
  
  bool OBForceFieldGhemical::SetGhemicalCharges()
  {
    OBAtom *a, *b;
    int bondtype;

    _mol.SetAutomaticPartialCharge(false);
    _mol.SetPartialChargesPerceived();

    // set all partial charges to 0.0
    FOR_ATOMS_OF_MOL (atom, _mol)
      atom->SetPartialCharge(0.0f);

    FOR_BONDS_OF_MOL (bond, _mol) {
      a = bond->GetBeginAtom();
      b = bond->GetEndAtom();	
      bondtype = bond->GetBondOrder(); 

      string _a(a->GetType());
      string _b(b->GetType());

      for (unsigned int idx=0; idx < _ffchargeparams.size(); idx++) {
        if (((_a == _ffchargeparams[idx]._a) && (_b == _ffchargeparams[idx]._b)) && (bondtype == _ffchargeparams[idx].ipar5)) {
          a->SetPartialCharge( - _ffchargeparams[idx].dpar1);
	  b->SetPartialCharge( + _ffchargeparams[idx].dpar1);
	} else if (((_a == _ffchargeparams[idx]._b) && (_b == _ffchargeparams[idx]._a)) && (bondtype == _ffchargeparams[idx].ipar5)) {
          a->SetPartialCharge( + _ffchargeparams[idx].dpar1);
	  b->SetPartialCharge( - _ffchargeparams[idx].dpar1);
	}
      }
    }

    return true;
  }

  double OBForceFieldGhemical::Energy()
  {
    double energy;
   
    energy = E_Bond();
    energy += E_Angle();
    energy += E_Torsion();
    energy += E_VDW();
    energy += E_Electrostatic();

    return energy;
  }
  
  OBFFParameter* OBForceFieldGhemical::GetParameterGhemical(int type, const char* a, const char* b, const char* c, const char* d, 
                                                            vector<OBFFParameter> &parameter)
  {
    OBFFParameter *par;
    if (a == NULL)
      return NULL;

    if (b == NULL) {
      string _a(a);
      for (unsigned int idx=0; idx < parameter.size(); idx++) 
        if ((_a == parameter[idx]._a) && (type == parameter[idx].ipar5)) {
          par = &parameter[idx];
          return par;
        }
      return NULL;
    }
    if (c == NULL) {
      string _a(a);
      string _b(b);
      for (unsigned int idx=0; idx < parameter.size(); idx++) {
        if (((_a == parameter[idx]._a) && (_b == parameter[idx]._b)) && (type == parameter[idx].ipar5) || 
            ((_a == parameter[idx]._b) && (_b == parameter[idx]._a)) && (type == parameter[idx].ipar5)) {
          par = &parameter[idx];
          return par;
        }
      }
      return NULL;
    }
    if (d == NULL) {
      string _a(a);
      string _b(b);
      string _c(c);
      for (unsigned int idx=0; idx < parameter.size(); idx++) {
        if (((_a == parameter[idx]._a) && (_b == parameter[idx]._b) && (_c == parameter[idx]._c)) && (type == parameter[idx].ipar5)|| 
            ((_a == parameter[idx]._c) && (_b == parameter[idx]._b) && (_c == parameter[idx]._a)) && (type == parameter[idx].ipar5)) {
          par = &parameter[idx];
          return par;
        }
      }
      return NULL;
    }
    string _a(a);
    string _b(b);
    string _c(c);
    string _d(d);
    for (unsigned int idx=0; idx < parameter.size(); idx++)
      if (((_a == parameter[idx]._a) && (_b == parameter[idx]._b) && (_c == parameter[idx]._c) && 
           (_d == parameter[idx]._d)) && (type == parameter[idx].ipar5) || 
          ((_a == parameter[idx]._d) && (_b == parameter[idx]._c) && (_c == parameter[idx]._b) && 
           (_d == parameter[idx]._a)) && (type == parameter[idx].ipar5)) {
        par = &parameter[idx];
        return par;
      }

    return NULL;
  }
  
  vector3 OBForceFieldGhemical::GetGradient(OBAtom *a, int terms)
  {
    vector<OBFFBondCalculationGhemical>::iterator i;
    vector<OBFFAngleCalculationGhemical>::iterator i2;
    vector<OBFFTorsionCalculationGhemical>::iterator i3;
    vector<OBFFVDWCalculationGhemical>::iterator i4;
    vector<OBFFElectrostaticCalculationGhemical>::iterator i5;

    vector3 grad(0.0f, 0.0f, 0.0f);
    
    if ((terms & OBFF_ENERGY) || (terms & OBFF_EBOND))
      for (i = _bondcalculations.begin(); i != _bondcalculations.end(); i++)
        if (((*i).a->GetIdx() == a->GetIdx()) || ((*i).b->GetIdx() == a->GetIdx()))
          grad += i->GetGradient(&*a);
    

    if ((terms & OBFF_ENERGY) || (terms & OBFF_EANGLE))
      for (i2 = _anglecalculations.begin(); i2 != _anglecalculations.end(); i2++)
        if (((*i2).a->GetIdx() == a->GetIdx()) || ((*i2).b->GetIdx() == a->GetIdx()) || ((*i2).c->GetIdx() == a->GetIdx()))
          grad += i2->GetGradient(&*a);
      
    if ((terms & OBFF_ENERGY) || (terms & OBFF_ETORSION))
      for (i3 = _torsioncalculations.begin(); i3 != _torsioncalculations.end(); i3++)
        if (((*i3).a->GetIdx() == a->GetIdx()) || ((*i3).b->GetIdx() == a->GetIdx()) || ((*i3).c->GetIdx() == a->GetIdx()) || ((*i3).d->GetIdx() == a->GetIdx()))
          grad += i3->GetGradient(&*a);
      
    if ((terms & OBFF_ENERGY) || (terms & OBFF_EVDW))
      for (i4 = _vdwcalculations.begin(); i4 != _vdwcalculations.end(); i4++)
        if (((*i4).a->GetIdx() == a->GetIdx()) || ((*i4).b->GetIdx() == a->GetIdx()))
          grad += i4->GetGradient(&*a);
    
    if ((terms & OBFF_ENERGY) || (terms & OBFF_EELECTROSTATIC))
      for (i5 = _electrostaticcalculations.begin(); i5 != _electrostaticcalculations.end(); i5++)
        if (((*i5).a->GetIdx() == a->GetIdx()) || ((*i5).b->GetIdx() == a->GetIdx()))
          grad += i5->GetGradient(&*a);

    return grad;
  }

  bool OBForceFieldGhemical::ValidateGradients ()
  {
    vector<OBFFBondCalculationGhemical>::iterator i1;
    vector<OBFFAngleCalculationGhemical>::iterator i2;
    vector<OBFFTorsionCalculationGhemical>::iterator i3;
    vector<OBFFVDWCalculationGhemical>::iterator i4;
    vector3 numgrad, anagrad, err;
    char logbuf[150];
    
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

      // OBFF_ETORSION
      numgrad = NumericalDerivative(&*a, OBFF_ETORSION);
      anagrad = GetGradient(&*a, OBFF_ETORSION);
      err = ValidateGradientError(numgrad, anagrad);

      sprintf(logbuf, "    torsion (%7.3f, %7.3f, %7.3f)  (%7.3f, %7.3f, %7.3f)  (%5.2f, %5.2f, %5.2f)", numgrad.x(), numgrad.y(), numgrad.z(), 
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
  }
} // end namespace OpenBabel

//! \file forcefieldghemical.cpp
//! \brief Ghemical force field

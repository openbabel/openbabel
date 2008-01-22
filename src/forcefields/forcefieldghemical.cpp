/**********************************************************************
forcefieldghemical.cpp - Ghemical force field.
 
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
#include <openbabel/mol.h>
#include "forcefieldghemical.h"

using namespace std;

namespace OpenBabel
{
  void OBFFBondCalculationGhemical::Compute(bool gradients)
  {
    vector3 da, db;
    double delta2, dE;

    //cout << "compute(" << gradients << ")" << endl;
    
    if (gradients) {
      da = a->GetVector();
      db = b->GetVector();
      rab = OBForceField::VectorLengthDerivative(da, db);
    } else
      rab = a->GetDistance(b);

    delta = rab - r0;
    delta2 = delta * delta;
    energy = kb * delta2;

    if (gradients) {
      dE = 2.0 * kb * delta;
      grada = dE * da; // - dE/drab * drab/da
      gradb = dE * db; // - dE/drab * drab/db
    }
  }
  
  double OBForceFieldGhemical::E_Bond(bool gradients)
  {
    vector<OBFFBondCalculationGhemical>::iterator i;
    double energy = 0.0;
        
    IF_OBFF_LOGLVL_HIGH {
      OBFFLog("\nB O N D   S T R E T C H I N G\n\n");
      OBFFLog("ATOM TYPES  BOND    BOND       IDEAL       FORCE\n");
      OBFFLog(" I    J     TYPE   LENGTH     LENGTH     CONSTANT      DELTA      ENERGY\n");
      OBFFLog("------------------------------------------------------------------------\n");
    }
 
    for (i = _bondcalculations.begin(); i != _bondcalculations.end(); ++i) {

      i->Compute(gradients);
      energy += i->GetEnergy();

      IF_OBFF_LOGLVL_HIGH {
        sprintf(_logbuf, "%s %s    %d   %8.3f   %8.3f     %8.3f   %8.3f   %8.3f\n", (*i).a->GetType(), (*i).b->GetType(), 
                (*i).bt, (*i).rab, (*i).r0, (*i).kb, (*i).delta, (*i).energy);
        OBFFLog(_logbuf);
      }
    }
    
    IF_OBFF_LOGLVL_MEDIUM {
      sprintf(_logbuf, "     TOTAL BOND STRETCHING ENERGY = %8.3f %s\n",  energy, GetUnit().c_str());
      OBFFLog(_logbuf);
    }
    return energy;
  }
  
  void OBFFAngleCalculationGhemical::Compute(bool gradients)
  {
    vector3 da, db, dc;
    double delta2, dE;

    if (gradients) {
      da = a->GetVector();
      db = b->GetVector();
      dc = c->GetVector();
      theta = OBForceField::VectorAngleDerivative(da, db, dc);  
    } else {
      theta = a->GetAngle(b->GetIdx(), c->GetIdx());
    }

    /*
    if (theta0 > 170.0) {
      delta = 1.0 + cos(theta * DEG_TO_RAD);

      energy = ka * delta * RAD_TO_DEG * RAD_TO_DEG;

      if (gradients) {
        grada = ka * da;
        gradb = ka * db;
        gradc = ka * dc;
      } 
    } else {
    */
      delta = theta - theta0;
      delta2 = delta * delta;
    
      energy = ka * delta2;

      if (gradients) {
        dE = 2.0 * ka * delta;
        grada = dE * da; // - dE/drab * drab/da
        gradb = dE * db; // - dE/drab * drab/db = - dE/drab * drab/da - dE/drab * drab/dc 
        gradc = dE * dc; // - dE/drab * drab/dc
      }
    /*
    }
    */
  }
  
  double OBForceFieldGhemical::E_Angle(bool gradients)
  {
    vector<OBFFAngleCalculationGhemical>::iterator i;
    double energy = 0.0;
        
    IF_OBFF_LOGLVL_HIGH {
      OBFFLog("\nA N G L E   B E N D I N G\n\n");
      OBFFLog("ATOM TYPES       VALENCE     IDEAL      FORCE\n");
      OBFFLog(" I    J    K      ANGLE      ANGLE     CONSTANT      DELTA      ENERGY\n");
      OBFFLog("-----------------------------------------------------------------------------\n");
    }
    
    for (i = _anglecalculations.begin(); i != _anglecalculations.end(); ++i) {

      i->Compute(gradients);
      energy += i->GetEnergy();
      
      IF_OBFF_LOGLVL_HIGH {
        sprintf(_logbuf, "%s %s %s  %8.3f   %8.3f     %8.3f   %8.3f   %8.3f\n", (*i).a->GetType(), (*i).b->GetType(), 
                (*i).c->GetType(), (*i).theta, (*i).theta0, (*i).ka, (*i).delta, (*i).energy);
        OBFFLog(_logbuf);
      }
    }
 
    IF_OBFF_LOGLVL_MEDIUM {
      sprintf(_logbuf, "     TOTAL ANGLE BENDING ENERGY = %8.3f %s\n", energy, GetUnit().c_str());
      OBFFLog(_logbuf);
    }
    return energy;
  }
  
  void OBFFTorsionCalculationGhemical::Compute(bool gradients)
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
    cosine2 = cos(2.0 * DEG_TO_RAD * tor);
    cosine3 = cos(3.0 * DEG_TO_RAD * tor);

    phi1 = 1.0 + cosine;
    phi2 = 1.0 - cosine2;
    phi3 = 1.0 + cosine3;

    energy = k1 * phi1 + k2 * phi2 + k3 * phi3;
    
    if (gradients) {
      sine = sin(DEG_TO_RAD * tor);
      sine2 = sin(2.0 * DEG_TO_RAD * tor);
      sine3 = sin(3.0 * DEG_TO_RAD * tor);
      dE = -k1 * sine + k2 * 2.0 * sine2 - k3 * 3.0 * sine3;
      grada = dE * da; // - dE/drab * drab/da
      gradb = dE * db; // - dE/drab * drab/db
      gradc = dE * dc; // - dE/drab * drab/dc
      gradd = dE * dd; // - dE/drab * drab/dd
    }
  }
  
  double OBForceFieldGhemical::E_Torsion(bool gradients) 
  {
    vector<OBFFTorsionCalculationGhemical>::iterator i;
    double energy = 0.0;
 
    IF_OBFF_LOGLVL_HIGH {
      OBFFLog("\nT O R S I O N A L\n\n");
      OBFFLog("----ATOM TYPES-----    FORCE              TORSION\n");
      OBFFLog(" I    J    K    L     CONSTANT     s       ANGLE    n    ENERGY\n");
      OBFFLog("----------------------------------------------------------------\n");
    }
    
    for (i = _torsioncalculations.begin(); i != _torsioncalculations.end(); ++i) {

      i->Compute(gradients);
      energy += i->GetEnergy();
      
      IF_OBFF_LOGLVL_HIGH {
        sprintf(_logbuf, "%s %s %s %s    %6.3f    %5.0f   %8.3f   %1.0f   %8.3f\n", (*i).a->GetType(), (*i).b->GetType(), 
                (*i).c->GetType(), (*i).d->GetType(), (*i).V, (*i).s, (*i).tor, (*i).n, (*i).energy);
        OBFFLog(_logbuf);
      }
    }

    IF_OBFF_LOGLVL_MEDIUM {
      sprintf(_logbuf, "     TOTAL TORSIONAL ENERGY = %8.3f %s\n", energy, GetUnit().c_str());
      OBFFLog(_logbuf);
    }

    return energy;
  }

  void OBFFVDWCalculationGhemical::Compute(bool gradients)
  {
    vector3 da, db;
    double term6, term12, dE, term7, term13;

    if (gradients) {
      da = a->GetVector();
      db = b->GetVector();
      rab = OBForceField::VectorLengthDerivative(da, db);
    } else
      rab = a->GetDistance(b);
    
    term7 = term13 = term6 = rab / ka;

    term6 = term6 * term6 * term6; // ^3
    term6 = term6 * term6; // ^6
    term12 = term6 * term6; // ^12
   
    energy = (1.0 / term12) - (1.0 / term6);
    
    if (gradients) { 
      term13 = term13 * term12; // ^13
      term7 = term7 * term6; // ^7
      dE = - (12.0 / ka) * (1.0 / term13) + (6.0 / kb) * (1.0 / term7);
      grada = dE * da; // - dE/drab * drab/da
      gradb = dE * db; // - dE/drab * drab/db
    }
  }
  
  double OBForceFieldGhemical::E_VDW(bool gradients)
  {
    vector<OBFFVDWCalculationGhemical>::iterator i;
    double energy = 0.0;
     
    IF_OBFF_LOGLVL_HIGH {
      OBFFLog("\nV A N   D E R   W A A L S\n\n");
      OBFFLog("ATOM TYPES\n");
      OBFFLog(" I    J        Rij       kij       ENERGY\n");
      OBFFLog("-----------------------------------------\n");
      //          XX   XX     -000.000  -000.000  -000.000  -000.000
    }
    
    for (i = _vdwcalculations.begin(); i != _vdwcalculations.end(); ++i) {
      
      i->Compute(gradients);
      energy += i->GetEnergy();
      
      IF_OBFF_LOGLVL_HIGH {
        sprintf(_logbuf, "%s %s   %8.3f  %8.3f  %8.3f\n", (*i).a->GetType(), (*i).b->GetType(), 
                (*i).rab, (*i).kab, (*i).energy);
        OBFFLog(_logbuf);
      }
    }

    IF_OBFF_LOGLVL_MEDIUM {
      sprintf(_logbuf, "     TOTAL VAN DER WAALS ENERGY = %8.3f %s\n", energy, GetUnit().c_str());
      OBFFLog(_logbuf);
    }

    return energy;
  }

  void OBFFElectrostaticCalculationGhemical::Compute(bool gradients)
  {
    vector3 da, db;
    double dE, rab2;

    if (gradients) {
      da = a->GetVector();
      db = b->GetVector();
      rab = OBForceField::VectorLengthDerivative(da, db);
    } else
      rab = a->GetDistance(b);

    energy = qq / rab;

    if (gradients) {
      rab2 = rab * rab;
      dE = -qq / rab2;
      grada = dE * da; // - dE/drab * drab/da
      gradb = dE * db; // - dE/drab * drab/db
    } 
  }
  
  double OBForceFieldGhemical::E_Electrostatic(bool gradients)
  {
    vector<OBFFElectrostaticCalculationGhemical>::iterator i;
    double energy = 0.0;
     
    IF_OBFF_LOGLVL_HIGH {
      OBFFLog("\nE L E C T R O S T A T I C   I N T E R A C T I O N S\n\n");
      OBFFLog("ATOM TYPES\n");
      OBFFLog(" I    J           Rij   332.17*QiQj  ENERGY\n");
      OBFFLog("-------------------------------------------\n");
      //            XX   XX     -000.000  -000.000  -000.000  
    }

    for (i = _electrostaticcalculations.begin(); i != _electrostaticcalculations.end(); ++i) {
      
      i->Compute(gradients);
      energy += i->GetEnergy();
      
      IF_OBFF_LOGLVL_HIGH {
        sprintf(_logbuf, "%s %s   %8.3f  %8.3f  %8.3f\n", (*i).a->GetType(), (*i).b->GetType(), 
                (*i).rab, (*i).qq, (*i).energy);
        OBFFLog(_logbuf);
      }
    }

    IF_OBFF_LOGLVL_MEDIUM {
      sprintf(_logbuf, "     TOTAL ELECTROSTATIC ENERGY = %8.3f %s\n", energy, GetUnit().c_str());
      OBFFLog(_logbuf);
    }

    return energy;
  }

  //***********************************************
  //Make a global instance
  OBForceFieldGhemical theForceFieldGhemical("Ghemical", true);
  //***********************************************

  OBForceFieldGhemical::~OBForceFieldGhemical()
  {
  }

  OBForceFieldGhemical &OBForceFieldGhemical::operator=(OBForceFieldGhemical &src)
  {
    _mol = src._mol;
    _init = src._init;

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

  bool OBForceFieldGhemical::SetupCalculations()
  {
    OBFFParameter *parameter;
    OBAtom *a, *b, *c, *d;
    
    IF_OBFF_LOGLVL_LOW
      OBFFLog("\nS E T T I N G   U P   C A L C U L A T I O N S\n\n");
 
    // 
    // Bond Calculations
    IF_OBFF_LOGLVL_LOW
      OBFFLog("SETTING UP BOND CALCULATIONS...\n");
    
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
            bondcalc.kb = KCAL_TO_KJ * 500.0;
            bondcalc.r0 = 1.100;

            _bondcalculations.push_back(bondcalc);

            IF_OBFF_LOGLVL_LOW {
              sprintf(_logbuf, "COULD NOT FIND PARAMETERS FOR BOND %s-%s, USING DEFAULT PARAMETERS\n", a->GetType(), b->GetType());
              OBFFLog(_logbuf);
            }

            continue;
          }
        }
      }
      bondcalc.kb = KCAL_TO_KJ * parameter->_dpar[1];
      bondcalc.r0 = parameter->_dpar[0];

      _bondcalculations.push_back(bondcalc);
    }

    //
    // Angle Calculations
    //
    IF_OBFF_LOGLVL_LOW
      OBFFLog("SETTING UP ANGLE CALCULATIONS...\n");
 
    OBFFAngleCalculationGhemical anglecalc;
 
    _anglecalculations.clear();
    
    FOR_ANGLES_OF_MOL(angle, _mol) {
      b = _mol.GetAtom((*angle)[0] + 1);
      a = _mol.GetAtom((*angle)[1] + 1);
      c = _mol.GetAtom((*angle)[2] + 1);
      
      anglecalc.a = a;
      anglecalc.b = b;
      anglecalc.c = c;

      parameter = GetParameter(a->GetType(), b->GetType(), c->GetType(), NULL, _ffangleparams);
      if (parameter == NULL) {
        parameter = GetParameter("FFFF", b->GetType(), c->GetType(), NULL, _ffangleparams);
        if (parameter == NULL) {
          parameter = GetParameter(a->GetType(), b->GetType(), "FFFF", NULL, _ffangleparams);
          if (parameter == NULL) {
            parameter = GetParameter("FFFF", b->GetType(), "FFFF", NULL, _ffangleparams);
            if (parameter == NULL) {
              anglecalc.ka = KCAL_TO_KJ * 0.020;
              anglecalc.theta0 = 120.0;
            
              _anglecalculations.push_back(anglecalc);
            
              IF_OBFF_LOGLVL_LOW {
                sprintf(_logbuf, "COULD NOT FIND PARAMETERS FOR ANGLE %s-%s-%s, USING DEFAULT PARAMETERS\n", a->GetType(), b->GetType(), c->GetType());
                OBFFLog(_logbuf);
              }

              continue;
            }
          }
        }
      }
      anglecalc.ka = KCAL_TO_KJ * parameter->_dpar[1];
      anglecalc.theta0 = parameter->_dpar[0];
      
      _anglecalculations.push_back(anglecalc);
    }
    
    //
    // Torsion Calculations
    //
    IF_OBFF_LOGLVL_LOW
      OBFFLog("SETTING UP TORSION CALCULATIONS...\n");
 
    OBFFTorsionCalculationGhemical torsioncalc;
    int torsiontype;
    int s;

    _torsioncalculations.clear();
 
    FOR_TORSIONS_OF_MOL(t, _mol) {
      a = _mol.GetAtom((*t)[0] + 1);
      b = _mol.GetAtom((*t)[1] + 1);
      c = _mol.GetAtom((*t)[2] + 1);
      d = _mol.GetAtom((*t)[3] + 1);
      OBBond *bc = _mol.GetBond(b, c);
      torsiontype = bc->GetBondOrder(); 
      if (bc->IsAromatic())
        torsiontype = 5;
      
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
              torsioncalc.V = 0.0;
              torsioncalc.s = 1.0;
              torsioncalc.n = 1.0;
              
              torsioncalc.k1 = 0.0;
              torsioncalc.k2 = 0.0;
              torsioncalc.k3 = 0.0;
              _torsioncalculations.push_back(torsioncalc);

              IF_OBFF_LOGLVL_LOW {
                sprintf(_logbuf, "COULD NOT FIND PARAMETERS FOR TORSION %s-%s-%s-%s, USING DEFAULT PARAMETERS\n", a->GetType(), b->GetType(), c->GetType(), d->GetType());
                OBFFLog(_logbuf);
              }

              continue;
            }
          }
        }
      }
      torsioncalc.V = KCAL_TO_KJ * parameter->_dpar[0];
      torsioncalc.s = parameter->_dpar[1];
      torsioncalc.n = parameter->_dpar[2];

      s = (int) (torsioncalc.s * torsioncalc.n);
      switch(s) {
      case +3:
        torsioncalc.k1 = 0.0;
        torsioncalc.k2 = 0.0;
        torsioncalc.k3 = torsioncalc.V;
        break;
      case +2:
        torsioncalc.k1 = 0.0;
        torsioncalc.k2 = -torsioncalc.V;
        torsioncalc.k3 = 0.0;
        break;
      case +1:
        torsioncalc.k1 = torsioncalc.V;
        torsioncalc.k2 = 0.0;
        torsioncalc.k3 = 0.0;
        break;
      case -1:
        torsioncalc.k1 = -torsioncalc.V;
        torsioncalc.k2 = 0.0;
        torsioncalc.k3 = 0.0;
        break;
      case -2:
        torsioncalc.k1 = 0.0;
        torsioncalc.k2 = torsioncalc.V;
        torsioncalc.k3 = 0.0;
        break;
      case -3:
        torsioncalc.k1 = 0.0;
        torsioncalc.k2 = 0.0;
        torsioncalc.k3 = -torsioncalc.V;
        break;
      }

      _torsioncalculations.push_back(torsioncalc);     
    }
    
    // 
    // VDW Calculations
    //
    IF_OBFF_LOGLVL_LOW
      OBFFLog("SETTING UP VAN DER WAALS CALCULATIONS...\n");
    
    OBFFVDWCalculationGhemical vdwcalc;
    OBFFParameter *parameter_a, *parameter_b;

    _vdwcalculations.clear();
    
    FOR_PAIRS_OF_MOL(p, _mol) {
      a = _mol.GetAtom((*p)[0]);
      b = _mol.GetAtom((*p)[1]);

      parameter_a = GetParameter(a->GetType(), NULL, NULL, NULL, _ffvdwparams);
      if (parameter_a == NULL) { // no vdw parameter -> use hydrogen
        vdwcalc.Ra = 1.5;
        vdwcalc.ka = 0.042;
	
        IF_OBFF_LOGLVL_LOW {
          sprintf(_logbuf, "COULD NOT FIND VDW PARAMETERS FOR ATOM %s, USING HYDROGEN VDW PARAMETERS\n", a->GetType());
          OBFFLog(_logbuf);
        }
      } else {
        vdwcalc.Ra = parameter_a->_dpar[0];
        vdwcalc.ka = parameter_a->_dpar[1];
      }

      parameter_b = GetParameter(b->GetType(), NULL, NULL, NULL, _ffvdwparams);
      if (parameter_b == NULL) { // no vdw parameter -> use hydrogen
        vdwcalc.Rb = 1.5;
        vdwcalc.kb = 0.042;
        
        IF_OBFF_LOGLVL_LOW {
          sprintf(_logbuf, "COULD NOT FIND VDW PARAMETERS FOR ATOM %s, USING HYDROGEN VDW PARAMETERS\n", b->GetType());
          OBFFLog(_logbuf);
        }
      } else {
        vdwcalc.Rb = parameter_b->_dpar[0];
        vdwcalc.kb = parameter_b->_dpar[1];
      }

      vdwcalc.a = &*a;
      vdwcalc.b = &*b;
     
      //this calculations only need to be done once for each pair, 
      //we do them now and save them for later use
      vdwcalc.kab = KCAL_TO_KJ * sqrt(vdwcalc.ka * vdwcalc.kb);
      
      // 1-4 scaling
      if (a->IsOneFour(b))
        vdwcalc.kab *= 0.5;
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

      vdwcalc.ka = (vdwcalc.Ra + vdwcalc.Rb) * pow(1.0 * vdwcalc.kab , 1.0 / 12.0);
      vdwcalc.kb = (vdwcalc.Ra + vdwcalc.Rb) * pow(2.0 * vdwcalc.kab , 1.0 / 6.0);
      
      _vdwcalculations.push_back(vdwcalc);
    }
    
    // 
    // Electrostatic Calculations
    //
    IF_OBFF_LOGLVL_LOW
      OBFFLog("SETTING UP ELECTROSTATIC CALCULATIONS...\n");
 
    OBFFElectrostaticCalculationGhemical elecalc;

    _electrostaticcalculations.clear();
    
    FOR_PAIRS_OF_MOL(p, _mol) {
      a = _mol.GetAtom((*p)[0]);
      b = _mol.GetAtom((*p)[1]);
      
      elecalc.qq = KCAL_TO_KJ * 332.17 * a->GetPartialCharge() * b->GetPartialCharge();
      
      if (elecalc.qq) {
        elecalc.a = &*a;
        elecalc.b = &*b;
        
        // 1-4 scaling
        if (a->IsOneFour(b))
          elecalc.qq *= 0.5;
	  
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
    ifstream ifs;
    if (OpenDatafile(ifs, "ghemical.prm").length() == 0) {
      obErrorLog.ThrowError(__FUNCTION__, "Cannot open ghemical.prm", obError);
      return false;
    }

    while (ifs.getline(buffer, 80)) {
      tokenize(vs, buffer);

      if (EQn(buffer, "bond", 4)) {
        parameter.clear();
        parameter._a = vs[1];
        parameter._b = vs[2];
        parameter._dpar.push_back(atof(vs[4].c_str())); // length
        parameter._dpar.push_back(atof(vs[5].c_str())); // force cte
        parameter._ipar.resize(1);
        if (EQn(vs[3].c_str(), "S", 1))
          parameter._ipar[0] = 1;
        if (EQn(vs[3].c_str(), "D", 1))
          parameter._ipar[0] = 2;
        if (EQn(vs[3].c_str(), "T", 1))
          parameter._ipar[0] = 3;
        if (EQn(vs[3].c_str(), "C", 1))
          parameter._ipar[0] = 5;
        _ffbondparams.push_back(parameter);
      }
      if (EQn(buffer, "angle", 5)) {
        parameter.clear();
        parameter._a = vs[1];
        parameter._b = vs[2];
        parameter._c = vs[3];
        parameter._dpar.push_back(atof(vs[5].c_str())); // angle
        parameter._dpar.push_back(atof(vs[6].c_str())); // force cte
        _ffangleparams.push_back(parameter);
      }
      if (EQn(buffer, "torsion", 7)) {
        parameter.clear();
        parameter._a = vs[1];
        parameter._b = vs[2];
        parameter._c = vs[3];
        parameter._d = vs[4];
        parameter._dpar.resize(3);
        parameter._dpar[0] = atof(vs[6].c_str()); // force cte
        parameter._dpar[2] = atof(vs[8].c_str()); // n
        if (EQn(vs[7].c_str(), "+", 1))
          parameter._dpar[1] = +1; // s
        else if (EQn(vs[7].c_str(), "-", 1))
          parameter._dpar[1] = -1; // s

        parameter._ipar.resize(1);
        if (EQn(vs[5].c_str(), "?S?", 3))
          parameter._ipar[0] = 1;
        else if (EQn(vs[5].c_str(), "?D?", 3))
          parameter._ipar[0] = 2;
        else if (EQn(vs[5].c_str(), "?T?", 3))
          parameter._ipar[0] = 3;
        else if (EQn(vs[5].c_str(), "?C?", 3))
          parameter._ipar[0] = 5;
        _fftorsionparams.push_back(parameter);
      }
      if (EQn(buffer, "vdw", 3)) {
        parameter.clear();
        parameter._a = vs[1];
        parameter._dpar.push_back(atof(vs[2].c_str())); // r
        parameter._dpar.push_back(atof(vs[3].c_str())); // force cte
        _ffvdwparams.push_back(parameter);
      }
      if (EQn(buffer, "charge", 6)) {
        parameter.clear();
        parameter._a = vs[1];
        parameter._b = vs[2];
        parameter._ipar.resize(1);
        if (EQn(vs[3].c_str(), "S", 1))
          parameter._ipar[0] = 1;
        else if (EQn(vs[3].c_str(), "D", 1))
          parameter._ipar[0] = 2;
        parameter._dpar.push_back(atof(vs[4].c_str())); // charge
        _ffchargeparams.push_back(parameter);
      }
    }
	
    if (ifs)
      ifs.close();
 
    return 0;
  }
  
  bool OBForceFieldGhemical::SetTypes()
  {
    vector<vector<int> > _mlist; //!< match list for atom typing
    vector<pair<OBSmartsPattern*,string> > _vexttyp; //!< external atom type rules
    vector<vector<int> >::iterator j;
    vector<pair<OBSmartsPattern*,string> >::iterator i;
    OBSmartsPattern *sp;
    vector<string> vs;
    char buffer[80];
 
    _mol.SetAtomTypesPerceived();
    
    // open data/ghemical.prm
    ifstream ifs;
    if (OpenDatafile(ifs, "ghemical.prm").length() == 0) {
      obErrorLog.ThrowError(__FUNCTION__, "Cannot open ghemical.prm", obError);
      return false;
    }

    while (ifs.getline(buffer, 80)) {
      if (EQn(buffer, "atom", 4)) {
      	tokenize(vs, buffer);

        sp = new OBSmartsPattern;
        if (sp->Init(vs[1]))
          _vexttyp.push_back(pair<OBSmartsPattern*,string> (sp,vs[2]));
        else {
          delete sp;
          sp = NULL;
          obErrorLog.ThrowError(__FUNCTION__, " Could not parse atom type table from ghemical.prm", obInfo);
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

    SetPartialCharges();
 
    IF_OBFF_LOGLVL_LOW {
      OBFFLog("\nA T O M   T Y P E S\n\n");
      OBFFLog("IDX\tTYPE\n");
      
      FOR_ATOMS_OF_MOL (a, _mol) {
        sprintf(_logbuf, "%d\t%s\n", a->GetIdx(), a->GetType());
        OBFFLog(_logbuf);
      }

      OBFFLog("\nC H A R G E S\n\n");
      OBFFLog("IDX\tCHARGE\n");
      
      FOR_ATOMS_OF_MOL (a, _mol) {
        sprintf(_logbuf, "%d\t%f\n", a->GetIdx(), a->GetPartialCharge());
        OBFFLog(_logbuf);
      }
    }
    
    // DEBUG (validation)
    //FOR_ATOMS_OF_MOL (a, _mol)
    //  if (atoi(a->GetType()) != 0)
    //    cout << "ATOMTYPE " << atoi(a->GetType()) << endl;
    //  else
    //    cout << "ATOMTYPE " << a->GetType() << endl;
 
    if (ifs)
      ifs.close();

    return true;
  }
  
  bool OBForceFieldGhemical::SetPartialCharges()
  {
    OBAtom *a, *b;
    int bondtype;

    _mol.SetAutomaticPartialCharge(false);
    _mol.SetPartialChargesPerceived();

    // set all partial charges to 0.0
    FOR_ATOMS_OF_MOL (atom, _mol)
      atom->SetPartialCharge(0.0);

    FOR_BONDS_OF_MOL (bond, _mol) {
      a = bond->GetBeginAtom();
      b = bond->GetEndAtom();	
      bondtype = bond->GetBondOrder(); 

      string _a(a->GetType());
      string _b(b->GetType());

      for (unsigned int idx=0; idx < _ffchargeparams.size(); ++idx) {
        if (((_a == _ffchargeparams[idx]._a) && (_b == _ffchargeparams[idx]._b)) && (bondtype == _ffchargeparams[idx]._ipar[0])) {
          a->SetPartialCharge(a->GetPartialCharge() - _ffchargeparams[idx]._dpar[0]);
          b->SetPartialCharge(b->GetPartialCharge() + _ffchargeparams[idx]._dpar[0]);
        } else if (((_a == _ffchargeparams[idx]._b) && (_b == _ffchargeparams[idx]._a)) && (bondtype == _ffchargeparams[idx]._ipar[0])) {
          a->SetPartialCharge(a->GetPartialCharge() + _ffchargeparams[idx]._dpar[0]);
          b->SetPartialCharge(b->GetPartialCharge() - _ffchargeparams[idx]._dpar[0]);
        }
      }
    }

    return true;
  }

  double OBForceFieldGhemical::Energy(bool gradients)
  {
    double energy;
    
    IF_OBFF_LOGLVL_MEDIUM
      OBFFLog("\nE N E R G Y\n\n");
 
    energy = E_Bond(gradients);
    energy += E_Angle(gradients);
    energy += E_Torsion(gradients);
    energy += E_VDW(gradients);
    energy += E_Electrostatic(gradients);

    IF_OBFF_LOGLVL_MEDIUM {
      sprintf(_logbuf, "\nTOTAL ENERGY = %8.3f %s\n", energy, GetUnit().c_str());
      OBFFLog(_logbuf);
    }

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
      for (unsigned int idx=0; idx < parameter.size(); ++idx) 
        if ((_a == parameter[idx]._a) && (type == parameter[idx]._ipar[0])) {
          par = &parameter[idx];
          return par;
        }
      return NULL;
    }
    if (c == NULL) {
      string _a(a);
      string _b(b);
      for (unsigned int idx=0; idx < parameter.size(); ++idx) {
        if ( ((_a == parameter[idx]._a) && (_b == parameter[idx]._b) &&
              (type == parameter[idx]._ipar[0])) || 
             ((_a == parameter[idx]._b) && (_b == parameter[idx]._a)) &&
             (type == parameter[idx]._ipar[0]) ) {
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
      for (unsigned int idx=0; idx < parameter.size(); ++idx) {
        if ( ((_a == parameter[idx]._a) && (_b == parameter[idx]._b) &&
              (_c == parameter[idx]._c) && 
              (type == parameter[idx]._ipar[0]))|| 
             ((_a == parameter[idx]._c) && (_b == parameter[idx]._b) &&
              (_c == parameter[idx]._a) && (type == parameter[idx]._ipar[0])) ) {
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

    for (unsigned int idx=0; idx < parameter.size(); ++idx) {
      if ( ((_a == parameter[idx]._a) && (_b == parameter[idx]._b) &&
             (_c == parameter[idx]._c) && (_d == parameter[idx]._d) &&
             (type == parameter[idx]._ipar[0])) || 
            ((_a == parameter[idx]._d) && (_b == parameter[idx]._c) && 
             (_c == parameter[idx]._b) && (_d == parameter[idx]._a) &&
             (type == parameter[idx]._ipar[0])) ) {
        par = &parameter[idx];
        return par;
      }
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
    if ((terms & OBFF_ENERGY) || (terms & OBFF_ETORSION))
      for (i3 = _torsioncalculations.begin(); i3 != _torsioncalculations.end(); ++i3)
        if (((*i3).a->GetIdx() == a->GetIdx()) || ((*i3).b->GetIdx() == a->GetIdx()) || 
	    ((*i3).c->GetIdx() == a->GetIdx()) || ((*i3).d->GetIdx() == a->GetIdx())) {
          i3->Compute(true);
          grad += i3->GetGradient(&*a);
        }
    if ((terms & OBFF_ENERGY) || (terms & OBFF_EVDW))
      for (i4 = _vdwcalculations.begin(); i4 != _vdwcalculations.end(); ++i4) 
        if (((*i4).a->GetIdx() == a->GetIdx()) || ((*i4).b->GetIdx() == a->GetIdx())) {
          i4->Compute(true);
          grad += i4->GetGradient(&*a);
    
        }
    if ((terms & OBFF_ENERGY) || (terms & OBFF_EELECTROSTATIC))
      for (i5 = _electrostaticcalculations.begin(); i5 != _electrostaticcalculations.end(); ++i5)
        if (((*i5).a->GetIdx() == a->GetIdx()) || ((*i5).b->GetIdx() == a->GetIdx())) {
          i5->Compute(true);
          grad += i5->GetGradient(&*a);
        }

    return grad;
  }

  bool OBForceFieldGhemical::ValidateGradients ()
  {
    vector3 numgrad, anagrad, err;
    
    OBFFLog("\nV A L I D A T E   G R A D I E N T S\n\n");
    OBFFLog("ATOM IDX      NUMERICAL GRADIENT           ANALYTICAL GRADIENT        REL. ERRROR (%)   \n");
    OBFFLog("----------------------------------------------------------------------------------------\n");
    //     "XX       (000.000, 000.000, 000.000)  (000.000, 000.000, 000.000)  (00.00, 00.00, 00.00)"
   
    FOR_ATOMS_OF_MOL (a, _mol) {

      // OBFF_ENERGY
      numgrad = NumericalDerivative(&*a, OBFF_ENERGY);
      anagrad = GetGradient(&*a, OBFF_ENERGY);
      err = ValidateGradientError(numgrad, anagrad);

      sprintf(_logbuf, "%2d       (%7.3f, %7.3f, %7.3f)  (%7.3f, %7.3f, %7.3f)  (%5.2f, %5.2f, %5.2f)\n", a->GetIdx(), numgrad.x(), numgrad.y(), numgrad.z(), 
              anagrad.x(), anagrad.y(), anagrad.z(), err.x(), err.y(), err.z());
      OBFFLog(_logbuf);

      // OBFF_EBOND
      numgrad = NumericalDerivative(&*a, OBFF_EBOND);
      anagrad = GetGradient(&*a, OBFF_EBOND);
      err = ValidateGradientError(numgrad, anagrad);

      sprintf(_logbuf, "    bond    (%7.3f, %7.3f, %7.3f)  (%7.3f, %7.3f, %7.3f)  (%5.2f, %5.2f, %5.2f)\n", numgrad.x(), numgrad.y(), numgrad.z(), 
              anagrad.x(), anagrad.y(), anagrad.z(), err.x(), err.y(), err.z());
      OBFFLog(_logbuf);
      
      // OBFF_EANGLE
      numgrad = NumericalDerivative(&*a, OBFF_EANGLE);
      anagrad = GetGradient(&*a, OBFF_EANGLE);
      err = ValidateGradientError(numgrad, anagrad);

      sprintf(_logbuf, "    angle   (%7.3f, %7.3f, %7.3f)  (%7.3f, %7.3f, %7.3f)  (%5.2f, %5.2f, %5.2f)\n", numgrad.x(), numgrad.y(), numgrad.z(), 
              anagrad.x(), anagrad.y(), anagrad.z(), err.x(), err.y(), err.z());
      OBFFLog(_logbuf);

      // OBFF_ETORSION
      numgrad = NumericalDerivative(&*a, OBFF_ETORSION);
      anagrad = GetGradient(&*a, OBFF_ETORSION);
      err = ValidateGradientError(numgrad, anagrad);

      sprintf(_logbuf, "    torsion (%7.3f, %7.3f, %7.3f)  (%7.3f, %7.3f, %7.3f)  (%5.2f, %5.2f, %5.2f)\n", numgrad.x(), numgrad.y(), numgrad.z(), 
              anagrad.x(), anagrad.y(), anagrad.z(), err.x(), err.y(), err.z());
      OBFFLog(_logbuf);

      // OBFF_EVDW
      numgrad = NumericalDerivative(&*a, OBFF_EVDW);
      anagrad = GetGradient(&*a, OBFF_EVDW);
      err = ValidateGradientError(numgrad, anagrad);

      sprintf(_logbuf, "    vdw     (%7.3f, %7.3f, %7.3f)  (%7.3f, %7.3f, %7.3f)  (%5.2f, %5.2f, %5.2f)\n", numgrad.x(), numgrad.y(), numgrad.z(), 
              anagrad.x(), anagrad.y(), anagrad.z(), err.x(), err.y(), err.z());
      OBFFLog(_logbuf);

      // OBFF_EELECTROSTATIC
      numgrad = NumericalDerivative(&*a, OBFF_EELECTROSTATIC);
      anagrad = GetGradient(&*a, OBFF_EELECTROSTATIC);
      err = ValidateGradientError(numgrad, anagrad);

      sprintf(_logbuf, "    electro (%7.3f, %7.3f, %7.3f)  (%7.3f, %7.3f, %7.3f)  (%5.2f, %5.2f, %5.2f)\n", numgrad.x(), numgrad.y(), numgrad.z(), 
              anagrad.x(), anagrad.y(), anagrad.z(), err.x(), err.y(), err.z());
      OBFFLog(_logbuf);
    }
    
    // For now, just return true. Should return false if validation fails.
    return true;
  }

} // end namespace OpenBabel

//! \file forcefieldghemical.cpp
//! \brief Ghemical force field

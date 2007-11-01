/**********************************************************************
forcefielduff.cpp - UFF force field.
 
Copyright (C) 2007 by Geoffrey Hutchison
Some portions Copyright (C) 2006-2007 by Tim Vandermeersch

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
#include "forcefielduff.h"

using namespace std;

// This implementation was created based on open code and reference websites:
// http://towhee.sourceforge.net/forcefields/uff.html
// http://rdkit.org/
// http://franklin.chm.colostate.edu/mmac/uff.html
// (for the last, use the Wayback Machine: http://www.archive.org/

// As well, the main UFF paper:
// Rappe, A. K., et. al.; J. Am. Chem. Soc. (1992) 114(25) p. 10024-10035.

namespace OpenBabel
{
  void OBFFBondCalculationUFF::Compute(bool gradients)
  {
    vector3 vab, da, db;
    double delta2, dE;

    if (gradients) {
      da = a->GetVector();
      db = b->GetVector();
      rab = OBForceField::VectorLengthDerivative(da, db);
    } else
      rab = a->GetDistance(b);

    delta = rab - r0; // we pre-compute the r0 below
    delta2 = delta * delta;
    energy = kb * delta2; // we fold the 1/2 into kb below

    if (gradients) {
      dE = 2.0 * kb * delta;
      grada = dE * da; // - dE/drab * drab/da
      gradb = dE * db; // - dE/drab * drab/db
    }
  }
  
  double OBForceFieldUFF::E_Bond(bool gradients)
  {
    vector<OBFFBondCalculationUFF>::iterator i;
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
        sprintf(logbuf, "%s %s    %d   %8.3f   %8.3f     %8.3f   %8.3f   %8.3f\n", (*i).a->GetType(), (*i).b->GetType(), 
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
  
  void OBFFAngleCalculationUFF::Compute(bool gradients)
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

    if (theta0 > 170.0) {
      delta = 1.0 + cos(theta * DEG_TO_RAD);

      energy = ka * delta * RAD_TO_DEG * RAD_TO_DEG;

      if (gradients) {
        grada = ka * da;
        gradb = ka * db;
        gradc = ka * dc;
      } 
    } else {
      delta = theta - theta0;
      delta2 = delta * delta;
    
      energy = ka * delta2;

      if (gradients) {
        dE = 2.0 * ka * delta;
        grada = dE * da; // - dE/drab * drab/da
        gradb = dE * db; // - dE/drab * drab/db = - dE/drab * drab/da - dE/drab * drab/dc 
        gradc = dE * dc; // - dE/drab * drab/dc
      }
    }
    energy = 0.0;
    grada = 0.0;
    gradb = 0.0;
    gradc = 0.0;
  }
  
  double OBForceFieldUFF::E_Angle(bool gradients)
  {
    vector<OBFFAngleCalculationUFF>::iterator i;
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
        sprintf(logbuf, "%s %s %s  %8.3f   %8.3f     %8.3f   %8.3f   %8.3f\n", (*i).a->GetType(), (*i).b->GetType(), 
                (*i).c->GetType(), (*i).theta, (*i).theta0, (*i).ka, (*i).delta, (*i).energy);
        OBFFLog(logbuf);
      }
    }
 
    IF_OBFF_LOGLVL_MEDIUM {
      sprintf(logbuf, "     TOTAL ANGLE BENDING ENERGY = %8.3f %s\n", energy, GetUnit().c_str());
      OBFFLog(logbuf);
    }
    return energy;
  }
  
  void OBFFTorsionCalculationUFF::Compute(bool gradients)
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
    energy = 0.0;
    grada = 0.0;
    gradb = 0.0;
    gradc = 0.0;
    gradd = 0.0;
  }
  
  double OBForceFieldUFF::E_Torsion(bool gradients) 
  {
    vector<OBFFTorsionCalculationUFF>::iterator i;
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
        sprintf(logbuf, "%s %s %s %s    %6.3f    %5.0f   %8.3f   %1.0f   %8.3f\n", (*i).a->GetType(), (*i).b->GetType(), 
                (*i).c->GetType(), (*i).d->GetType(), (*i).V, (*i).s, (*i).tor, (*i).n, (*i).energy);
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
  void OBFFOOPCalculationUFF::Compute(bool gradients)
  {
    double angle2;

    angle = Point2PlaneAngle(d->GetVector(), a->GetVector(), b->GetVector(), c->GetVector());
    angle2 = angle * angle;
    
    //    energy = 0.043844 * 0.5 * koop * angle2;
    energy = 0.0;
  }

  double OBForceFieldUFF::E_OOP(bool gradients) 
  {
     vector<OBFFOOPCalculationUFF>::iterator i;
     double energy = 0.0;
    
     IF_OBFF_LOGLVL_HIGH {
       OBFFLog("O U T - O F - P L A N E   B E N D I N G");
       OBFFLog("ATOM TYPES             FF       OOP     FORCE ");
       OBFFLog(" I    J    K    L     CLASS    ANGLE   CONSTANT     ENERGY");
       OBFFLog("----------------------------------------------------------");
     }

     for (i = _oopcalculations.begin(); i != _oopcalculations.end(); ++i) {
      i->Compute(gradients);
      energy += i->GetEnergy();
      
       IF_OBFF_LOGLVL_HIGH {
         sprintf(logbuf, "%2d   %2d   %2d   %2d      0   %8.3f   %8.3f     %8.3f\n", atoi((*i).a->GetType()), atoi((*i).b->GetType()), atoi((*i).c->GetType()), atoi((*i).d->GetType()), 
                 (*i).angle, (*i).koop, (*i).energy);
         OBFFLog(logbuf);
       }
     }

     IF_OBFF_LOGLVL_HIGH {
       sprintf(logbuf, "     TOTAL OUT-OF-PLANE BENDING ENERGY = %8.3f\n", energy);
       OBFFLog(logbuf);
     }
     return energy;
  }

  void OBFFVDWCalculationUFF::Compute(bool gradients)
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
   
    energy = (1.0 / term12) - (2.0 / term6);
    
    if (gradients) { 
      term13 = term13 * term12; // ^13
      term7 = term7 * term6; // ^7
      dE = - (12.0 / ka) * (1.0 / term13) + (12.0 / kb) * (1.0 / term7);
      grada = dE * da; // - dE/drab * drab/da
      gradb = dE * db; // - dE/drab * drab/db
    }
  }
  
  double OBForceFieldUFF::E_VDW(bool gradients)
  {
    vector<OBFFVDWCalculationUFF>::iterator i;
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
        sprintf(logbuf, "%s %s   %8.3f  %8.3f  %8.3f\n", (*i).a->GetType(), (*i).b->GetType(), 
                (*i).rab, (*i).kab, (*i).energy);
        OBFFLog(logbuf);
      }
    }

    IF_OBFF_LOGLVL_MEDIUM {
      sprintf(logbuf, "     TOTAL VAN DER WAALS ENERGY = %8.3f %s\n", energy, GetUnit().c_str());
      OBFFLog(logbuf);
    }

    return energy;
  }

  void OBFFElectrostaticCalculationUFF::Compute(bool gradients)
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
  
  double OBForceFieldUFF::E_Electrostatic(bool gradients)
  {
    vector<OBFFElectrostaticCalculationUFF>::iterator i;
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
        sprintf(logbuf, "%s %s   %8.3f  %8.3f  %8.3f\n", (*i).a->GetType(), (*i).b->GetType(), 
                (*i).rab, (*i).qq, (*i).energy);
        OBFFLog(logbuf);
      }
    }

    IF_OBFF_LOGLVL_MEDIUM {
      sprintf(logbuf, "     TOTAL ELECTROSTATIC ENERGY = %8.3f %s\n", energy, GetUnit().c_str());
      OBFFLog(logbuf);
    }

    return energy;
  }

  //***********************************************
  //Make a global instance
  OBForceFieldUFF theForceFieldUFF("UFF", true);
  //***********************************************

  OBForceFieldUFF::~OBForceFieldUFF()
  {
  }

  OBForceFieldUFF &OBForceFieldUFF::operator=(OBForceFieldUFF &src)
  {
    _mol = src._mol;

    _ffparams    = src._ffparams;

    _bondcalculations          = src._bondcalculations;
    _anglecalculations         = src._anglecalculations;
    _torsioncalculations       = src._torsioncalculations;
    _oopcalculations           = src._oopcalculations;
    _vdwcalculations           = src._vdwcalculations;
    _electrostaticcalculations = src._electrostaticcalculations;

    return *this;
  }

  bool OBForceFieldUFF::Setup(OBMol &mol)
  {
    _mol = mol;
    
    SetUFFTypes();

    if (!SetupCalculations())
      return false;
    
    return true;
  }
  
  bool OBForceFieldUFF::SetupCalculations()
  {
    OBFFParameter *parameterA, *parameterB, *parameterC, *parameterD;
    OBAtom *a, *b, *c, *d;
    
    IF_OBFF_LOGLVL_LOW
      OBFFLog("\nS E T T I N G   U P   C A L C U L A T I O N S\n\n");
 
    // 
    // Bond Calculations
    IF_OBFF_LOGLVL_LOW
      OBFFLog("SETTING UP BOND CALCULATIONS...\n");
    
    OBFFBondCalculationUFF bondcalc;
    double bondorder;
    double ri, rj, rbo, ren;
    double chiI, chiJ;

    _bondcalculations.clear();
    
    FOR_BONDS_OF_MOL(bond, _mol) {
      a = bond->GetBeginAtom();
      b = bond->GetEndAtom();
      bondorder = bond->GetBondOrder(); 
      if (bond->IsAromatic())
        bondorder = 1.5;
      if (bond->IsAmide())
        bondorder = 1.41;
      
      bondcalc.a = a;
      bondcalc.b = b;

      parameterA = GetParameterUFF(a->GetType(), _ffparams);
      parameterB = GetParameterUFF(b->GetType(), _ffparams);
      ri = parameterA->_dpar[0];
      rj = parameterB->_dpar[0];
      chiI = parameterA->_dpar[8];
      chiJ = parameterB->_dpar[8];

      // precompute the equilibrium geometry
      // From equation 3
      rbo = -0.1332*(ri+rj)*log(bondorder);
      // From equation 4
      ren = ri*rj*(pow((sqrt(chiI) - sqrt(chiJ)),2.0)) / (chiI*ri + chiJ*rj);
      // From equation 2
      // NOTE: See http://towhee.sourceforge.net/forcefields/uff.html
      // There is a typo in the published paper
      bondcalc.r0 = ri + rj + rbo - ren;

      // here we fold the 1/2 into the kij from equation 1a
      // Otherwise, this is equation 6 from the UFF paper.
      bondcalc.kb = (0.5 * KCAL_TO_KJ * 644.12 
        * parameterA->_dpar[5] * parameterB->_dpar[5])
        / (bondcalc.r0 * bondcalc.r0 * bondcalc.r0);

      _bondcalculations.push_back(bondcalc);
    }

    //
    // Angle Calculations
    //
    IF_OBFF_LOGLVL_LOW
      OBFFLog("SETTING UP ANGLE CALCULATIONS...\n");
 
    OBFFAngleCalculationUFF anglecalc;
 
    _anglecalculations.clear();
    
    FOR_ANGLES_OF_MOL(angle, _mol) {
      b = _mol.GetAtom((*angle)[0] + 1);
      a = _mol.GetAtom((*angle)[1] + 1);
      c = _mol.GetAtom((*angle)[2] + 1);
      
      anglecalc.a = a;
      anglecalc.b = b;
      anglecalc.c = c;

      parameterA = GetParameterUFF(a->GetType(), _ffparams);
      parameterB = GetParameterUFF(b->GetType(), _ffparams);
      parameterC = GetParameterUFF(c->GetType(), _ffparams);

      anglecalc.ka = KCAL_TO_KJ * parameterA->_dpar[1];
      anglecalc.theta0 = parameterA->_dpar[0];
      
      _anglecalculations.push_back(anglecalc);
    }
    
    //
    // Torsion Calculations
    //
    IF_OBFF_LOGLVL_LOW
      OBFFLog("SETTING UP TORSION CALCULATIONS...\n");
 
    OBFFTorsionCalculationUFF torsioncalc;
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

      parameterA = GetParameterUFF(a->GetType(), _ffparams);
      parameterB = GetParameterUFF(b->GetType(), _ffparams);
      parameterC = GetParameterUFF(c->GetType(), _ffparams);
      parameterD = GetParameterUFF(d->GetType(), _ffparams);

      torsioncalc.V = KCAL_TO_KJ * parameterA->_dpar[0];
      torsioncalc.s = parameterA->_dpar[1];
      torsioncalc.n = parameterA->_dpar[2];

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
    
    OBFFVDWCalculationUFF vdwcalc;

    _vdwcalculations.clear();
    
    FOR_PAIRS_OF_MOL(p, _mol) {
      a = _mol.GetAtom((*p)[0]);
      b = _mol.GetAtom((*p)[1]);

      parameterA = GetParameterUFF(a->GetType(), _ffparams);
      parameterB = GetParameterUFF(b->GetType(), _ffparams);

      vdwcalc.Ra = parameterA->_dpar[2];
      vdwcalc.ka = parameterA->_dpar[3];
      vdwcalc.Rb = parameterB->_dpar[2];
      vdwcalc.kb = parameterB->_dpar[3];

      vdwcalc.a = &*a;
      vdwcalc.b = &*b;
     
      //this calculations only need to be done once for each pair, 
      //we do them now and save them for later use
      vdwcalc.kab = KCAL_TO_KJ * sqrt(vdwcalc.ka * vdwcalc.kb);
      
      // 1-4 scaling
//       if (a->IsOneFour(b))
//         vdwcalc.kab *= 0.5;

      vdwcalc.ka = (vdwcalc.Ra + vdwcalc.Rb) * pow(1.0 * vdwcalc.kab , 1.0 / 12.0);
      vdwcalc.kb = (vdwcalc.Ra + vdwcalc.Rb) * pow(2.0 * vdwcalc.kab , 1.0 / 6.0);
      
      _vdwcalculations.push_back(vdwcalc);
    }
    
    // 
    // Electrostatic Calculations
    //
    IF_OBFF_LOGLVL_LOW
      OBFFLog("SETTING UP ELECTROSTATIC CALCULATIONS...\n");
 
    OBFFElectrostaticCalculationUFF elecalc;

    _electrostaticcalculations.clear();
    
    // Note that while the UFF paper mentions an electrostatic term,
    // it does not actually use it. Both Towhee and the UFF FAQ
    // discourage the use of electrostatics with UFF.
    // If you wanted to use this term, uncomment the lines below

//     FOR_PAIRS_OF_MOL(p, _mol) {
//       a = _mol.GetAtom((*p)[0]);
//       b = _mol.GetAtom((*p)[1]);
      
//       elecalc.qq = KCAL_TO_KJ * 332.0637 * a->GetPartialCharge() * b->GetPartialCharge();
      
//       if (elecalc.qq) {
//         elecalc.a = &*a;
//         elecalc.b = &*b;
        
//         _electrostaticcalculations.push_back(elecalc);
//       }
//     }

    return true;
  }

  bool OBForceFieldUFF::ParseParamFile()
  {
    vector<string> vs;
    char buffer[80];
    
    OBFFParameter parameter;

    // open data/UFF.prm
    ifstream ifs;
    if (OpenDatafile(ifs, "UFF.prm").length() == 0) {
      obErrorLog.ThrowError(__FUNCTION__, "Cannot open UFF.prm", obError);
      return false;
    }

    while (ifs.getline(buffer, 256)) {
      tokenize(vs, buffer);
      if (vs.size() < 13)
        continue;
      
      if (EQn(buffer, "param", 5)) {
        // set up all parameters from this
        parameter.clear();
        parameter._a = vs[1]; // atom type
        parameter._dpar.push_back(atof(vs[2].c_str())); // r1
        parameter._dpar.push_back(atof(vs[3].c_str())); // theta0
        parameter._dpar.push_back(atof(vs[4].c_str())); // x1
        parameter._dpar.push_back(atof(vs[5].c_str())); // D1
        parameter._dpar.push_back(atof(vs[6].c_str())); // zeta
        parameter._dpar.push_back(atof(vs[7].c_str())); // Z1
        parameter._dpar.push_back(atof(vs[8].c_str())); // Vi
        parameter._dpar.push_back(atof(vs[9].c_str())); // Uj
        parameter._dpar.push_back(atof(vs[10].c_str())); // Xi
        parameter._dpar.push_back(atof(vs[11].c_str())); // Hard
        parameter._dpar.push_back(atof(vs[12].c_str())); // Radius

        char coord = vs[1][2]; // 3rd character of atom type
        cerr << " coordination: " << coord << endl;
        switch (coord) {
        case '1': // linear
          parameter._ipar.push_back(1);
          break;
        case '2': // trigonal planar (sp2)
          parameter._ipar.push_back(2);
          break;
        case '3': // tetrahedral (sp3)
          parameter._ipar.push_back(3);
          break;
        case '4': // square planar
          parameter._ipar.push_back(4);
          break;
        case '5': // trigonal bipyramidal -- not actually in parameterization
          parameter._ipar.push_back(5);
          break;
        case '6': // octahedral
          parameter._ipar.push_back(6);
          break;
        default: // general case (unknown coordination)
          // These atoms appear to generally be linear coordination like Cl
          parameter._ipar.push_back(1);
        }
        
        _ffparams.push_back(parameter);
      }
    }
	
    if (ifs)
      ifs.close();
 
    return 0;
  }
  
  bool OBForceFieldUFF::SetUFFTypes()
  {
    vector<vector<int> > _mlist; //!< match list for atom typing
    vector<pair<OBSmartsPattern*,string> > _vexttyp; //!< external atom type rules
    vector<vector<int> >::iterator j;
    vector<pair<OBSmartsPattern*,string> >::iterator i;
    OBSmartsPattern *sp;
    vector<string> vs;
    char buffer[80];
 
    _mol.SetAtomTypesPerceived();
    
    // open data/UFF.prm
    ifstream ifs;
    if (OpenDatafile(ifs, "UFF.prm").length() == 0) {
      obErrorLog.ThrowError(__FUNCTION__, "Cannot open UFF.prm", obError);
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
          obErrorLog.ThrowError(__FUNCTION__, " Could not parse atom type table from UFF.prm", obInfo);
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
 
    IF_OBFF_LOGLVL_LOW {
      OBFFLog("\nA T O M   T Y P E S\n\n");
      OBFFLog("IDX\tTYPE\n");
      
      FOR_ATOMS_OF_MOL (a, _mol) {
        sprintf(logbuf, "%d\t%s\n", a->GetIdx(), a->GetType());
        OBFFLog(logbuf);
      }

    }
     
    if (ifs)
      ifs.close();

    return true;
  }
  
  double OBForceFieldUFF::Energy(bool gradients)
  {
    double energy;
    
    IF_OBFF_LOGLVL_MEDIUM
      OBFFLog("\nE N E R G Y\n\n");
 
    energy = E_Bond(gradients);
    energy += E_Angle(gradients);
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
  
  OBFFParameter* OBForceFieldUFF::GetParameterUFF(std::string a, vector<OBFFParameter> &parameter)
  {
    OBFFParameter *par;

    for (unsigned int idx=0; idx < parameter.size(); ++idx) {
      if (a == parameter[idx]._a) {
        return &parameter[idx];
      }
    }
    return NULL;
  }
  
  vector3 OBForceFieldUFF::GetGradient(OBAtom *a, int terms)
  {
    vector<OBFFBondCalculationUFF>::iterator i;
    vector<OBFFAngleCalculationUFF>::iterator i2;
    vector<OBFFTorsionCalculationUFF>::iterator i3;
    vector<OBFFVDWCalculationUFF>::iterator i4;
    vector<OBFFElectrostaticCalculationUFF>::iterator i5;
    vector<OBFFOOPCalculationUFF>::iterator i6;

    vector3 grad(0.0, 0.0, 0.0);
    
    if ((terms & OBFF_ENERGY) || (terms & OBFF_EBOND))
      for (i = _bondcalculations.begin(); i != _bondcalculations.end(); ++i)
        if (((*i).a->GetIdx() == a->GetIdx()) || ((*i).b->GetIdx() == a->GetIdx()))
          grad += i->GetGradient(&*a);
    
    if ((terms & OBFF_ENERGY) || (terms & OBFF_EANGLE))
      for (i2 = _anglecalculations.begin(); i2 != _anglecalculations.end(); ++i2)
        if (((*i2).a->GetIdx() == a->GetIdx()) || ((*i2).b->GetIdx() == a->GetIdx()) || ((*i2).c->GetIdx() == a->GetIdx()))
          grad += i2->GetGradient(&*a);
      
    if ((terms & OBFF_ENERGY) || (terms & OBFF_ETORSION))
      for (i3 = _torsioncalculations.begin(); i3 != _torsioncalculations.end(); ++i3)
        if (((*i3).a->GetIdx() == a->GetIdx()) || ((*i3).b->GetIdx() == a->GetIdx()) || ((*i3).c->GetIdx() == a->GetIdx()) || ((*i3).d->GetIdx() == a->GetIdx()))
          grad += i3->GetGradient(&*a);
      
    if ((terms & OBFF_ENERGY) || (terms & OBFF_EVDW))
      for (i4 = _vdwcalculations.begin(); i4 != _vdwcalculations.end(); ++i4)
        if (((*i4).a->GetIdx() == a->GetIdx()) || ((*i4).b->GetIdx() == a->GetIdx()))
          grad += i4->GetGradient(&*a);
    
    if ((terms & OBFF_ENERGY) || (terms & OBFF_EELECTROSTATIC))
      for (i5 = _electrostaticcalculations.begin(); i5 != _electrostaticcalculations.end(); ++i5)
        if (((*i5).a->GetIdx() == a->GetIdx()) || ((*i5).b->GetIdx() == a->GetIdx()))
          grad += i5->GetGradient(&*a);
    
    if ((terms & OBFF_ENERGY) || (terms & OBFF_EOOP))
      for (i6 = _oopcalculations.begin(); i6 != _oopcalculations.end(); ++i6)
        if (((*i6).a->GetIdx() == a->GetIdx()) || ((*i6).b->GetIdx() == a->GetIdx()) || ((*i6).c->GetIdx() == a->GetIdx()) || ((*i6).d->GetIdx() == a->GetIdx()))
          grad += i6->GetGradient(&*a);    
    
    return grad;
  }

  bool OBForceFieldUFF::ValidateGradients ()
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

      sprintf(logbuf, "%2d       (%7.3f, %7.3f, %7.3f)  (%7.3f, %7.3f, %7.3f)  (%5.2f, %5.2f, %5.2f)\n", a->GetIdx(), numgrad.x(), numgrad.y(), numgrad.z(), 
              anagrad.x(), anagrad.y(), anagrad.z(), err.x(), err.y(), err.z());
      OBFFLog(logbuf);

      // OBFF_EBOND
      numgrad = NumericalDerivative(&*a, OBFF_EBOND);
      anagrad = GetGradient(&*a, OBFF_EBOND);
      err = ValidateGradientError(numgrad, anagrad);

      sprintf(logbuf, "    bond    (%7.3f, %7.3f, %7.3f)  (%7.3f, %7.3f, %7.3f)  (%5.2f, %5.2f, %5.2f)\n", numgrad.x(), numgrad.y(), numgrad.z(), 
              anagrad.x(), anagrad.y(), anagrad.z(), err.x(), err.y(), err.z());
      OBFFLog(logbuf);
      
      // OBFF_EANGLE
      numgrad = NumericalDerivative(&*a, OBFF_EANGLE);
      anagrad = GetGradient(&*a, OBFF_EANGLE);
      err = ValidateGradientError(numgrad, anagrad);

      sprintf(logbuf, "    angle   (%7.3f, %7.3f, %7.3f)  (%7.3f, %7.3f, %7.3f)  (%5.2f, %5.2f, %5.2f)\n", numgrad.x(), numgrad.y(), numgrad.z(), 
              anagrad.x(), anagrad.y(), anagrad.z(), err.x(), err.y(), err.z());
      OBFFLog(logbuf);

      // OBFF_ETORSION
      numgrad = NumericalDerivative(&*a, OBFF_ETORSION);
      anagrad = GetGradient(&*a, OBFF_ETORSION);
      err = ValidateGradientError(numgrad, anagrad);

      sprintf(logbuf, "    torsion (%7.3f, %7.3f, %7.3f)  (%7.3f, %7.3f, %7.3f)  (%5.2f, %5.2f, %5.2f)\n", numgrad.x(), numgrad.y(), numgrad.z(), 
              anagrad.x(), anagrad.y(), anagrad.z(), err.x(), err.y(), err.z());
      OBFFLog(logbuf);

      // OBFF_EVDW
      numgrad = NumericalDerivative(&*a, OBFF_EVDW);
      anagrad = GetGradient(&*a, OBFF_EVDW);
      err = ValidateGradientError(numgrad, anagrad);

      sprintf(logbuf, "    vdw     (%7.3f, %7.3f, %7.3f)  (%7.3f, %7.3f, %7.3f)  (%5.2f, %5.2f, %5.2f)\n", numgrad.x(), numgrad.y(), numgrad.z(), 
              anagrad.x(), anagrad.y(), anagrad.z(), err.x(), err.y(), err.z());
      OBFFLog(logbuf);

      // OBFF_EELECTROSTATIC
      numgrad = NumericalDerivative(&*a, OBFF_EELECTROSTATIC);
      anagrad = GetGradient(&*a, OBFF_EELECTROSTATIC);
      err = ValidateGradientError(numgrad, anagrad);

      sprintf(logbuf, "    electro (%7.3f, %7.3f, %7.3f)  (%7.3f, %7.3f, %7.3f)  (%5.2f, %5.2f, %5.2f)\n", numgrad.x(), numgrad.y(), numgrad.z(), 
              anagrad.x(), anagrad.y(), anagrad.z(), err.x(), err.y(), err.z());
      OBFFLog(logbuf);

      // OBFF_EOOP
      numgrad = NumericalDerivative(&*a, OBFF_EOOP);
      anagrad = GetGradient(&*a, OBFF_EOOP);
      err = ValidateGradientError(numgrad, anagrad);
      
      sprintf(logbuf, "    oop (%7.3f, %7.3f, %7.3f)  (%7.3f, %7.3f, %7.3f)  (%5.2f, %5.2f, %5.2f)\n", numgrad.x(), numgrad.y(), numgrad.z(), 
              anagrad.x(), anagrad.y(), anagrad.z(), err.x(), err.y(), err.z());
      OBFFLog(logbuf);      
    }
    
    // For now, just return true. Should return false if validation fails.
    return true;
  }

} // end namespace OpenBabel

//! \file forcefieldUFF.cpp
//! \brief UFF force field

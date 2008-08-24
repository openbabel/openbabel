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
#include <openbabel/locale.h>

#include "forcefieldghemical.h"

using namespace std;

namespace OpenBabel
{
  template<bool gradients>
  double OBForceFieldGhemical::E_Bond()
  {
    vector<OBFFBondCalculationGhemical>::iterator i;
    double energy = 0.0;
        
    IF_OBFF_LOGLVL_HIGH {
      OBFFLog("\nB O N D   S T R E T C H I N G\n\n");
      OBFFLog("ATOM TYPES  BOND    BOND       IDEAL       FORCE\n");
      OBFFLog(" I    J     TYPE   LENGTH     LENGTH     CONSTANT      DELTA      ENERGY\n");
      OBFFLog("------------------------------------------------------------------------\n");
    }

    unsigned int idxA, idxB;
    double rab, delta, delta2, e;
    Eigen::Vector3d Fa, Fb;
    for (i = _bondcalculations.begin(); i != _bondcalculations.end(); ++i) {
      idxA = i->a->GetIdx() - 1;
      idxB = i->b->GetIdx() - 1;
      if (OBForceField::IgnoreCalculation(idxA+1, idxB+1)) {
        continue;
      }

      if (gradients) {
        rab = VectorBondDerivative(GetPositions()[idxA], GetPositions()[idxB], Fa, Fb);
        delta = rab - i->r0;
      
        const double dE = 2.0 * i->kb * delta;
      
        Fa *= dE;
        Fb *= dE;
        GetGradients()[idxA] += Fa;
        GetGradients()[idxB] += Fb;
      } else {
        const Eigen::Vector3d ab = GetPositions()[idxA] - GetPositions()[idxB];
        rab = ab.norm();
        delta = rab - i->r0;
      }

      delta2 = delta * delta;
      e = i->kb * delta2;
      energy += e;

      IF_OBFF_LOGLVL_HIGH {
        snprintf(_logbuf, BUFF_SIZE, "%s %s    %d   %8.3f   %8.3f     %8.3f   %8.3f   %8.3f\n", 
            i->a->GetType(), i->b->GetType(), 
            i->bt, rab, i->r0, i->kb, delta, e);
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
  double OBForceFieldGhemical::E_Angle()
  {
    vector<OBFFAngleCalculationGhemical>::iterator i;
    double energy = 0.0;
        
    IF_OBFF_LOGLVL_HIGH {
      OBFFLog("\nA N G L E   B E N D I N G\n\n");
      OBFFLog("ATOM TYPES       VALENCE     IDEAL      FORCE\n");
      OBFFLog(" I    J    K      ANGLE      ANGLE     CONSTANT      DELTA      ENERGY\n");
      OBFFLog("-----------------------------------------------------------------------------\n");
    }
  
    unsigned int idxA, idxB, idxC;
    double theta, delta, delta2, e;
    Eigen::Vector3d Fa, Fb, Fc;
    for (i = _anglecalculations.begin(); i != _anglecalculations.end(); ++i) {
      idxA = i->a->GetIdx() - 1;
      idxB = i->b->GetIdx() - 1;
      idxC = i->c->GetIdx() - 1;
      if (OBForceField::IgnoreCalculation(idxA+1, idxB+1, idxC+1)) {
        continue;
      }
 
      if (gradients) {
        theta = VectorAngleDerivative(GetPositions()[idxA], GetPositions()[idxB], GetPositions()[idxC], Fa, Fb, Fc);
        delta = theta - i->theta0;
      
        const double dE = RAD_TO_DEG * 2.0 * i->ka * delta;
      
        Fa *= dE;
        Fb *= dE;
        Fc *= dE;
        GetGradients()[idxA] += Fa;
        GetGradients()[idxB] += Fb;
        GetGradients()[idxC] += Fc;
      } else {
        const Eigen::Vector3d ab = GetPositions()[idxA] - GetPositions()[idxB];
        const Eigen::Vector3d bc = GetPositions()[idxC] - GetPositions()[idxB];
        theta = VectorAngle(ab, bc);
        delta = theta - i->theta0;
      }

      if (!isfinite(theta))
        theta = 0.0; // doesn't explain why GetAngle is returning NaN but solves it for us;

      delta2 = delta * delta;
      e = i->ka * delta2;
      energy += e;
      
      IF_OBFF_LOGLVL_HIGH {
        snprintf(_logbuf, BUFF_SIZE, "%s %s %s  %8.3f   %8.3f     %8.3f   %8.3f   %8.3f\n", 
            i->a->GetType(), i->b->GetType(), i->c->GetType(), 
            theta, i->theta0, i->ka, delta, e);
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
  double OBForceFieldGhemical::E_Torsion() 
  {
    vector<OBFFTorsionCalculationGhemical>::iterator i;
    double energy = 0.0;
 
    IF_OBFF_LOGLVL_HIGH {
      OBFFLog("\nT O R S I O N A L\n\n");
      OBFFLog("----ATOM TYPES-----    FORCE              TORSION\n");
      OBFFLog(" I    J    K    L     CONSTANT     s       ANGLE    n    ENERGY\n");
      OBFFLog("----------------------------------------------------------------\n");
    }
  
    unsigned int idxA, idxB, idxC, idxD;
    double tor, e;
    Eigen::Vector3d Fa, Fb, Fc, Fd;
    for (i = _torsioncalculations.begin(); i != _torsioncalculations.end(); ++i) {
      idxA = i->a->GetIdx() - 1;
      idxB = i->b->GetIdx() - 1;
      idxC = i->c->GetIdx() - 1;
      idxD = i->d->GetIdx() - 1;
      if (OBForceField::IgnoreCalculation(idxA+1, idxB+1, idxC+1, idxD+1)) {
        continue;
      }
 
      if (gradients) {
        tor = DEG_TO_RAD * VectorTorsionDerivative(GetPositions()[idxA], GetPositions()[idxB],
            GetPositions()[idxC], GetPositions()[idxD], Fa, Fb, Fc, Fd);
        if (!isfinite(tor))
          tor = 1.0e-3;
      
        const double sine = sin(tor);
        const double sine2 = sin(2.0 * tor);
        const double sine3 = sin(3.0 * tor);
 
        const double dE = i->k1 * sine - i->k2 * 2.0 * sine2 + i->k3 * 3.0 * sine3;
     
        Fa *= dE;
        Fb *= dE;
        Fc *= dE;
        Fd *= dE;
        GetGradients()[idxA] += Fa;
        GetGradients()[idxB] += Fb;
        GetGradients()[idxC] += Fc;
        GetGradients()[idxD] += Fd;
      } else {
        tor = DEG_TO_RAD * VectorTorsion(GetPositions()[idxA], GetPositions()[idxB],
            GetPositions()[idxC], GetPositions()[idxD]);
        if (!isfinite(tor)) // stop any NaN or infinity
          tor = 1.0e-3; // rather than NaN
      }

      const double cosine = cos(tor);
      const double cosine2 = cos(2.0 * tor);
      const double cosine3 = cos(3.0 * tor);
      const double phi1 = 1.0 + cosine;
      const double phi2 = 1.0 - cosine2;
      const double phi3 = 1.0 + cosine3;

      e = i->k1 * phi1 + i->k2 * phi2 + i->k3 * phi3;
      energy += e;
      
      IF_OBFF_LOGLVL_HIGH {
        snprintf(_logbuf, BUFF_SIZE, "%s %s %s %s    %6.3f    %5.0f   %8.3f   %1.0f   %8.3f\n", 
            i->a->GetType(), i->b->GetType(), i->c->GetType(), i->d->GetType(), 
            i->V, i->s, tor, i->n, e);
        OBFFLog(_logbuf);
      }
    }

    IF_OBFF_LOGLVL_MEDIUM {
      snprintf(_logbuf, BUFF_SIZE, "     TOTAL TORSIONAL ENERGY = %8.3f %s\n", energy, GetUnit().c_str());
      OBFFLog(_logbuf);
    }

    return energy;
  }

  template<bool gradients>
  double OBForceFieldGhemical::E_VDW()
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
    
    unsigned int idxA, idxB;
    unsigned int j = 0;
    double rab, term, term6, term7, term12, term13, e, dE;
    Eigen::Vector3d Fa, Fb;
    for (i = _vdwcalculations.begin(); i != _vdwcalculations.end(); ++i, ++j) {
      // Cut-off check
      if (IsCutOffEnabled())
        if (!GetVDWPairs().BitIsSet(j)) 
          continue;

      idxA = i->a->GetIdx() - 1;
      idxB = i->b->GetIdx() - 1;
      if (OBForceField::IgnoreCalculation(idxA+1, idxB+1)) {
        continue;
      }
 
      if (gradients) {
        rab = VectorDistanceDerivative(GetPositions()[idxA], GetPositions()[idxB], Fa, Fb);
        term = rab / i->ka;

        term6 = term * term * term; // ^3
        term6 = term6 * term6; // ^6
        term12 = term6 * term6; // ^12
        term13 = term * term12; // ^13
        term7 = term * term6; // ^7
        dE = - (12.0 / i->ka) * (1.0 / term13) + (6.0 / i->ka) * (1.0 / term7);
        Fa *= dE;
        Fb *= dE;
        GetGradients()[idxA] += Fa;
        GetGradients()[idxB] += Fb;
      } else {
        const Eigen::Vector3d ab = GetPositions()[idxA] - GetPositions()[idxB];
        rab = ab.norm();
        term = rab / i->ka;

        term6 = term * term * term; // ^3
        term6 = term6 * term6; // ^6
        term12 = term6 * term6; // ^12
      } 
      
      e = (1.0 / term12) - (1.0 / term6);
      energy += e;
      
      IF_OBFF_LOGLVL_HIGH {
        snprintf(_logbuf, BUFF_SIZE, "%s %s   %8.3f  %8.3f  %8.3f\n", i->a->GetType(), i->b->GetType(), 
            rab, i->kab, e);
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
  double OBForceFieldGhemical::E_Electrostatic()
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

    unsigned int idxA, idxB;
    unsigned int j = 0;
    double rab, e;
    double rab2, dE;
    Eigen::Vector3d Fa, Fb;
    for (i = _electrostaticcalculations.begin(); i != _electrostaticcalculations.end(); ++i, ++j) {
      // Cut-off check
      if (IsCutOffEnabled())
        if (!GetElePairs().BitIsSet(j)) 
          continue;
 
      idxA = i->a->GetIdx() - 1;
      idxB = i->b->GetIdx() - 1;
      if (OBForceField::IgnoreCalculation(idxA+1, idxB+1)) {
        continue;
      }
 
      if (gradients) {
        rab = VectorDistanceDerivative(GetPositions()[idxA], GetPositions()[idxB], Fa, Fb);
        //const double rab2 = rab * rab;
        rab2 = rab * rab;
        //const double dE = -(i->qq) / rab2;
        dE = -(i->qq) / rab2;
        Fa *= dE;
        Fb *= dE;
        GetGradients()[idxA] += Fa;
        GetGradients()[idxB] += Fb;
      } else {
        const Eigen::Vector3d ab = GetPositions()[idxA] - GetPositions()[idxB];
        rab = ab.norm();
      }
    
      if (IsNearZero(rab, 1.0e-3))
        rab = 1.0e-3;

      e = i->qq / rab;
      energy += e;
      
      IF_OBFF_LOGLVL_HIGH {
        snprintf(_logbuf, BUFF_SIZE, "%s %s   %8.3f  %8.3f  %8.3f\n", i->a->GetType(), i->b->GetType(), 
                rab, i->qq, e);
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
  OBForceFieldGhemical theForceFieldGhemical("Ghemical", true);
  //***********************************************

  OBForceFieldGhemical::~OBForceFieldGhemical()
  {
  }

  bool OBForceFieldGhemical::SetupCalculations()
  {
    IF_OBFF_LOGLVL_LOW
      OBFFLog("\nS E T T I N G   U P   C A L C U L A T I O N S\n\n");

    if (!SetupBondedCalculations())
      return false;
    if (!SetupNonBondedCalculations())
      return false;
   
    return true; 
  }

  bool OBForceFieldGhemical::SetupBondedCalculations()
  {
    OBFFParameter *parameter;
    OBAtom *a, *b, *c, *d;
    
    OBMol *mol = GetMolecule();
    OBFFConstraints constraints = GetConstraints();
    vector<OBBitVec> intraGroup = GetIntraGroup();
    vector<OBBitVec> interGroup = GetInterGroup();
    vector<pair<OBBitVec, OBBitVec> > interGroups = GetInterGroups();
 
    // 
    // Bond Calculations
    IF_OBFF_LOGLVL_LOW
      OBFFLog("SETTING UP BOND CALCULATIONS...\n");
    
    OBFFBondCalculationGhemical bondcalc;
    int bondtype;

    _bondcalculations.clear();
    
    FOR_BONDS_OF_MOL(bond, mol) {
      a = bond->GetBeginAtom();
      b = bond->GetEndAtom();

      // skip this bond if the atoms are ignored 
      if ( constraints.IsIgnored(a->GetIdx()) || constraints.IsIgnored(b->GetIdx()) )
        continue;

      // if there are any groups specified, check if the two bond atoms are in a single intraGroup
      if (HasGroups()) {
        bool validBond = false;
        for (unsigned int i=0; i < intraGroup.size(); ++i) {
          if (intraGroup[i].BitIsOn(a->GetIdx()) && intraGroup[i].BitIsOn(b->GetIdx()))
            validBond = true;
        }
        if (!validBond)
          continue;
      }
      
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
              snprintf(_logbuf, BUFF_SIZE, "COULD NOT FIND PARAMETERS FOR BOND %s-%s, USING DEFAULT PARAMETERS\n", 
                  a->GetType(), b->GetType());
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
    
    FOR_ANGLES_OF_MOL(angle, mol) {
      b = mol->GetAtom((*angle)[0] + 1);
      a = mol->GetAtom((*angle)[1] + 1);
      c = mol->GetAtom((*angle)[2] + 1);
      
      // skip this angle if the atoms are ignored 
      if ( constraints.IsIgnored(a->GetIdx()) || constraints.IsIgnored(b->GetIdx()) || constraints.IsIgnored(c->GetIdx()) ) 
        continue;
 
      // if there are any groups specified, check if the three angle atoms are in a single intraGroup
      if (HasGroups()) {
        bool validAngle = false;
        for (unsigned int i=0; i < intraGroup.size(); ++i) {
          if (intraGroup[i].BitIsOn(a->GetIdx()) && intraGroup[i].BitIsOn(b->GetIdx()) && 
              intraGroup[i].BitIsOn(c->GetIdx()))
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
                snprintf(_logbuf, BUFF_SIZE, "COULD NOT FIND PARAMETERS FOR ANGLE %s-%s-%s, USING DEFAULT PARAMETERS\n", 
                    a->GetType(), b->GetType(), c->GetType());
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
 
    FOR_TORSIONS_OF_MOL(t, mol) {
      a = mol->GetAtom((*t)[0] + 1);
      b = mol->GetAtom((*t)[1] + 1);
      c = mol->GetAtom((*t)[2] + 1);
      d = mol->GetAtom((*t)[3] + 1);

      // skip this torsion if the atoms are ignored 
      if ( constraints.IsIgnored(a->GetIdx()) || constraints.IsIgnored(b->GetIdx()) ||
           constraints.IsIgnored(c->GetIdx()) || constraints.IsIgnored(d->GetIdx()) ) 
        continue;
 
      // if there are any groups specified, check if the four torsion atoms are in a single intraGroup
      if (HasGroups()) {
        bool validTorsion = false;
        for (unsigned int i=0; i < intraGroup.size(); ++i) {
          if (intraGroup[i].BitIsOn(a->GetIdx()) && intraGroup[i].BitIsOn(b->GetIdx()) && 
              intraGroup[i].BitIsOn(c->GetIdx()) && intraGroup[i].BitIsOn(d->GetIdx()))
            validTorsion = true;
        }
        if (!validTorsion)
          continue;
      }
 
      OBBond *bc = mol->GetBond(b, c);
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
                snprintf(_logbuf, BUFF_SIZE, "COULD NOT FIND PARAMETERS FOR TORSION %s-%s-%s-%s, "
                    "USING DEFAULT PARAMETERS\n", a->GetType(), b->GetType(), c->GetType(), d->GetType());
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
  
    return true;
  } 

  bool OBForceFieldGhemical::SetupNonBondedCalculations()
  {
    OBAtom *a, *b;
    
    OBMol *mol = GetMolecule();
    OBFFConstraints constraints = GetConstraints();
    vector<OBBitVec> interGroup = GetInterGroup();
    vector<pair<OBBitVec, OBBitVec> > interGroups = GetInterGroups();
 
 
    // 
    // VDW Calculations
    //
    IF_OBFF_LOGLVL_LOW
      OBFFLog("SETTING UP VAN DER WAALS CALCULATIONS...\n");
    
    OBFFVDWCalculationGhemical vdwcalc;
    OBFFParameter *parameter_a, *parameter_b;

    _vdwcalculations.clear();
    
    FOR_PAIRS_OF_MOL(p, mol) {
      a = mol->GetAtom((*p)[0]);
      b = mol->GetAtom((*p)[1]);

      // skip this vdw if the atoms are ignored 
      if ( constraints.IsIgnored(a->GetIdx()) || constraints.IsIgnored(b->GetIdx()) )
        continue;
  
      // if there are any groups specified, check if the two atoms are in a single interGroup or if
      // two two atoms are in one of the interGroups pairs.
      if (HasGroups()) {
        bool validVDW = false;
        for (unsigned int i=0; i < interGroup.size(); ++i) {
          if (interGroup[i].BitIsOn(a->GetIdx()) && interGroup[i].BitIsOn(b->GetIdx())) 
            validVDW = true;
        }
        for (unsigned int i=0; i < interGroups.size(); ++i) {
          if (interGroups[i].first.BitIsOn(a->GetIdx()) && interGroups[i].second.BitIsOn(b->GetIdx())) 
            validVDW = true;
          if (interGroups[i].first.BitIsOn(b->GetIdx()) && interGroups[i].second.BitIsOn(a->GetIdx())) 
            validVDW = true;
        }
 
        if (!validVDW)
          continue;
      }
 
      parameter_a = GetParameter(a->GetType(), NULL, NULL, NULL, _ffvdwparams);
      if (parameter_a == NULL) { // no vdw parameter -> use hydrogen
        vdwcalc.Ra = 1.5;
        vdwcalc.ka = 0.042;
	
        IF_OBFF_LOGLVL_LOW {
          snprintf(_logbuf, BUFF_SIZE, "COULD NOT FIND VDW PARAMETERS FOR ATOM %s, USING HYDROGEN VDW PARAMETERS\n", 
              a->GetType());
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
          snprintf(_logbuf, BUFF_SIZE, "COULD NOT FIND VDW PARAMETERS FOR ATOM %s, USING HYDROGEN VDW PARAMETERS\n", 
              b->GetType());
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
    
    FOR_PAIRS_OF_MOL(p, mol) {
      a = mol->GetAtom((*p)[0]);
      b = mol->GetAtom((*p)[1]);
 
      // skip this ele if the atoms are ignored 
      if ( constraints.IsIgnored(a->GetIdx()) || constraints.IsIgnored(b->GetIdx()) )
        continue;
  
      // if there are any groups specified, check if the two atoms are in a single interGroup or if
      // two two atoms are in one of the interGroups pairs.
      if (HasGroups()) {
        bool validEle = false;
        for (unsigned int i=0; i < interGroup.size(); ++i) {
          if (interGroup[i].BitIsOn(a->GetIdx()) && interGroup[i].BitIsOn(b->GetIdx())) 
            validEle = true;
        }
        for (unsigned int i=0; i < interGroups.size(); ++i) {
          if (interGroups[i].first.BitIsOn(a->GetIdx()) && interGroups[i].second.BitIsOn(b->GetIdx())) 
            validEle = true;
          if (interGroups[i].first.BitIsOn(b->GetIdx()) && interGroups[i].second.BitIsOn(a->GetIdx())) 
            validEle = true;
        }
 
        if (!validEle)
          continue;
      }
 
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

    // Set the locale for number parsing to avoid locale issues: PR#1785463
    obLocale.SetLocale();

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

    // return the locale to the original one
    obLocale.RestoreLocale();
 
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
 
    OBMol *mol = GetMolecule();

    mol->SetAtomTypesPerceived();
    
    // open data/ghemical.prm
    ifstream ifs;
    if (OpenDatafile(ifs, "ghemical.prm").length() == 0) {
      obErrorLog.ThrowError(__FUNCTION__, "Cannot open ghemical.prm", obError);
      return false;
    }

    // Set the locale for number parsing to avoid locale issues: PR#1785463
    obLocale.SetLocale();

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
      }
    }
    
    for (i = _vexttyp.begin();i != _vexttyp.end();++i) {
      if (i->first->Match(*mol)) {
        _mlist = i->first->GetMapList();
        for (j = _mlist.begin();j != _mlist.end();++j) {
          mol->GetAtom((*j)[0])->SetType(i->second);
        }
      }
    }
 
    SetPartialCharges();
 
    IF_OBFF_LOGLVL_LOW {
      OBFFLog("\nA T O M   T Y P E S\n\n");
      OBFFLog("IDX\tTYPE\n");
      
      FOR_ATOMS_OF_MOL (a, mol) {
        snprintf(_logbuf, BUFF_SIZE, "%d\t%s\n", a->GetIdx(), a->GetType());
        OBFFLog(_logbuf);
      }

      OBFFLog("\nC H A R G E S\n\n");
      OBFFLog("IDX\tCHARGE\n");
      
      FOR_ATOMS_OF_MOL (a, mol) {
        snprintf(_logbuf, BUFF_SIZE, "%d\t%f\n", a->GetIdx(), a->GetPartialCharge());
        OBFFLog(_logbuf);
      }
    }
    
    // DEBUG (validation)
    //FOR_ATOMS_OF_MOL (a, mol)
    //  if (atoi(a->GetType()) != 0)
    //    cout << "ATOMTYPE " << atoi(a->GetType()) << endl;
    //  else
    //    cout << "ATOMTYPE " << a->GetType() << endl;
 
    if (ifs)
      ifs.close();

    // return the locale to the original one
    obLocale.RestoreLocale();

    return true;
  }
  
  bool OBForceFieldGhemical::SetPartialCharges()
  {
    OBMol *mol = GetMolecule();
    OBAtom *a, *b;
    int bondtype;

    mol->SetAutomaticPartialCharge(false);
    mol->SetPartialChargesPerceived();

    // set all partial charges to 0.0
    FOR_ATOMS_OF_MOL (atom, mol)
      atom->SetPartialCharge(0.0);

    FOR_BONDS_OF_MOL (bond, mol) {
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
    double energy = 0.0;
    
    IF_OBFF_LOGLVL_MEDIUM
      OBFFLog("\nE N E R G Y\n\n");
 
    if (gradients) {
      ClearGradients();
      energy  = E_Bond<true>();
      energy += E_Angle<true>();
      energy += E_Torsion<true>();
      energy += E_VDW<true>();
      energy += E_Electrostatic<true>();
    } else {
      energy  = E_Bond<false>();
      energy += E_Angle<false>();
      energy += E_Torsion<false>();
      energy += E_VDW<false>();
      energy += E_Electrostatic<false>();
    }

    IF_OBFF_LOGLVL_MEDIUM {
      snprintf(_logbuf, BUFF_SIZE, "\nTOTAL ENERGY = %8.3f %s\n", energy, GetUnit().c_str());
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
  
  bool OBForceFieldGhemical::ValidateGradients ()
  {
    OBMol *mol = GetMolecule();
    Eigen::Vector3d numgrad, anagrad, err;
    int idx;
    
    bool passed = true; // set to false if any component fails
    
    OBFFLog("\nV A L I D A T E   G R A D I E N T S\n\n");
    OBFFLog("ATOM IDX      NUMERICAL GRADIENT           ANALYTICAL GRADIENT        REL. ERROR (%)   \n");
    OBFFLog("----------------------------------------------------------------------------------------\n");
    //     "XX       (000.000, 000.000, 000.000)  (000.000, 000.000, 000.000)  (00.00, 00.00, 00.00)"


    FOR_ATOMS_OF_MOL (a, mol) {
      idx = (a->GetIdx() - 1);

      // OBFF_ENERGY
      //numgrad = NumericalDerivative(&*a, OBFF_ENERGY);
      Energy(); // compute
      anagrad = GetGradients()[idx];
      err = ValidateGradientError(numgrad, anagrad);

      snprintf(_logbuf, BUFF_SIZE, "%2d       (%7.3f, %7.3f, %7.3f)  (%7.3f, %7.3f, %7.3f)  (%5.2f, %5.2f, %5.2f)\n", a->GetIdx(), numgrad.x(), numgrad.y(), numgrad.z(), 
              anagrad.x(), anagrad.y(), anagrad.z(), err.x(), err.y(), err.z());
      OBFFLog(_logbuf);
      
      // OBFF_EBOND
      //numgrad = NumericalDerivative(&*a, OBFF_EBOND);
      ClearGradients();
      E_Bond();
      anagrad = GetGradients()[idx];
      err = ValidateGradientError(numgrad, anagrad);

      snprintf(_logbuf, BUFF_SIZE, "    bond    (%7.3f, %7.3f, %7.3f)  (%7.3f, %7.3f, %7.3f)  (%5.2f, %5.2f, %5.2f)\n", numgrad.x(), numgrad.y(), numgrad.z(), 
              anagrad.x(), anagrad.y(), anagrad.z(), err.x(), err.y(), err.z());
      OBFFLog(_logbuf);
      // 8% tolerance here because some bonds have slight instability
      if (err.x() > 8.0 || err.y() > 8.0 || err.z() > 8.0)
        passed = false;
      
      // OBFF_EANGLE
      //numgrad = NumericalDerivative(&*a, OBFF_EANGLE);
      ClearGradients();
      E_Angle();
      anagrad = GetGradients()[idx];
      err = ValidateGradientError(numgrad, anagrad);

      snprintf(_logbuf, BUFF_SIZE, "    angle   (%7.3f, %7.3f, %7.3f)  (%7.3f, %7.3f, %7.3f)  (%5.2f, %5.2f, %5.2f)\n", numgrad.x(), numgrad.y(), numgrad.z(), 
              anagrad.x(), anagrad.y(), anagrad.z(), err.x(), err.y(), err.z());
      OBFFLog(_logbuf);
      if (err.x() > 5.0 || err.y() > 5.0 || err.z() > 5.0)
        passed = false;

      // OBFF_ETORSION
      //numgrad = NumericalDerivative(&*a, OBFF_ETORSION);
      ClearGradients();
      E_Torsion();
      anagrad = GetGradients()[idx];
      err = ValidateGradientError(numgrad, anagrad);

      snprintf(_logbuf, BUFF_SIZE, "    torsion (%7.3f, %7.3f, %7.3f)  (%7.3f, %7.3f, %7.3f)  (%5.2f, %5.2f, %5.2f)\n", numgrad.x(), numgrad.y(), numgrad.z(), 
              anagrad.x(), anagrad.y(), anagrad.z(), err.x(), err.y(), err.z());
      OBFFLog(_logbuf);
      if (err.x() > 5.0 || err.y() > 5.0 || err.z() > 5.0)
        passed = false;

      // OBFF_EVDW
      //numgrad = NumericalDerivative(&*a, OBFF_EVDW);
      ClearGradients();
      E_VDW();
      anagrad = GetGradients()[idx];
      err = ValidateGradientError(numgrad, anagrad);

      snprintf(_logbuf, BUFF_SIZE, "    vdw     (%7.3f, %7.3f, %7.3f)  (%7.3f, %7.3f, %7.3f)  (%5.2f, %5.2f, %5.2f)\n", numgrad.x(), numgrad.y(), numgrad.z(), 
              anagrad.x(), anagrad.y(), anagrad.z(), err.x(), err.y(), err.z());
      OBFFLog(_logbuf);
      if (err.x() > 5.0 || err.y() > 5.0 || err.z() > 5.0)
        passed = false;

      // OBFF_EELECTROSTATIC
      //numgrad = NumericalDerivative(&*a, OBFF_EELECTROSTATIC);
      ClearGradients();
      E_Electrostatic();
      anagrad = GetGradients()[idx];
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

//! \file forcefieldghemical.cpp
//! \brief Ghemical force field

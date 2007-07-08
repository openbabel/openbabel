/**********************************************************************
forcefield.cpp - Handle OBForceField class.

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

      development status:
      - src/forcefield.cpp
      - LineSearch(): works, but better algoritmh would be better for performance
      - SteepestDescent(): finished
      - ConjugateGradients(): finished
      - SystematicRotorSearch(): not all combinations are tested but works

      - src/forcefields/forcefieldghemical.cpp
      - Atom typing: finished
      - Charges: finished
      - Energy terms: finished
      - Analytical gradients: finished
      - Validation: finished

***********************************************************************/
#include <openbabel/babelconfig.h>

#include <openbabel/forcefield.h>

#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obiter.h>
#include <openbabel/math/matrix3x3.h>
#include <openbabel/rotamer.h>
#include <openbabel/rotor.h>

using namespace std;

namespace OpenBabel
{
  /** \class OBForceField forcefield.h <openbabel/forcefield.h>
      \brief Base class for molecular mechanics force fields

      The OBForceField class is the base class for molecular mechanics in Open Babel.
      Classes derived from the OBForceField implement specific force fields (Ghemical, MMFF94, ...).
      Other classes such as OBFFParameter, OBFFCalculation and its derived classes are 
      only for internal use. As a user interested in using the available force fields 
      in Open Babel, you don't need these classes. The rest of this short introduction
      is aimed at these users. For information on how to implement additional force
      fields, see the wiki pages or post your questions to the openbabel-devel mailing
      list.

      Before we can start using a force field, we must first select it and set it up.
      This is illustrated in the first example below. The Setup procedure assigns atom
      types, charges and parameters. There are several reasons why this may fail, a log
      message will be written to the logfile before Setup() returns false.

      The force field classes use their own logging functions. You can set the logfile 
      using SetLogFile() and set the log level using SetLogLevel(). If needed you can
      also write to the logfile using OBFFLog(). There are four log levels: 
      BFF_LOGLVL_NONE, OBFF_LOGLVL_LOW, OBFF_LOGLVL_MEDIUM, OBFF_LOGLVL_HIGH.
      See the API documentation to know what kind of output each function writes to
      the logfile for the different log levels.

      Below are two examples which explain the basics. 
      
      This piece of code can always be used to find available forcefields:
      \code
      FOR_EACH(OBForceField, iter) {
      cout << "forcefield ID: " << iter.ID() << endl;
      }
      \endcode

      Calculate the energy for the structure in mol using the Ghemical forcefield.
      \code
      #include <openbabel/forcefield.h>
      #include <openbabel/mol.h>

      // See OBConversion class to fill the mol object.
      OBMol mol;
      // Select the forcefield, this returns a pointer that we 
      // will later use to access the forcefield functions.
      OBForceField* pFF = OBForceField::FindForceField("Ghemical");

      // Make sure we have a valid pointer
      if (!pFF)
      // exit...
     
      // Set the logfile (can also be &cout or &cerr)
      pFF->SetLogFile(&cerr);
      // Set the log level. See indivual functions to know
      // what kind of output each function produces for the
      // different log levels.
      pFF->SetLogLevel(OBFF_LOGLVL_HIGH);

      // We need to setup the forcefield before we can use it. Setup()
      // returns false if it failes to find the atom types, parameters, ...
      if (!pFF->Setup(mol)) {
      cerr << "ERROR: could not setup force field." << endl;
      }
      
      // Calculate the energy. The output will be written to the
      // logfile specified by SetLogFile()
      pFF->Energy();
      \endcode
 
      Minimize the structure in mol using conjugate gradients.
      \code
      #include <openbabel/forcefield.h>
      #include <openbabel/mol.h>

      OBMol mol;
      OBForceField* pFF = OBForceField::FindForceField("Ghemical");

      // Make sure we have a valid pointer
      if (!pFF)
      // exit...
      
      pFF->SetLogFile(&cerr);
      pFF->SetLogLevel(OBFF_LOGLVL_LOW);
      if (!pFF->Setup(mol)) {
      cerr << "ERROR: could not setup force field." << endl;
      }
      
      // Perform the actual minimization, maximum 1000 steps 
      pFF->ConjugateGradients(1000);
      \endcode
  **/

  OBFFParameter* OBForceField::GetParameter(int a, int b, int c, int d, vector<OBFFParameter> &parameter)
  {
    OBFFParameter *par;

    if (!b)
      for (unsigned int idx=0; idx < parameter.size(); idx++)
        if (a == parameter[idx].a) {
          par = &parameter[idx];
          return par;
        }

    if (!c)
      for (unsigned int idx=0; idx < parameter.size(); idx++)
        if (((a == parameter[idx].a) && (b == parameter[idx].b)) || ((a == parameter[idx].b) && (b == parameter[idx].a))) {
          par = &parameter[idx];
          return par;
        }

    if (!d)
      for (unsigned int idx=0; idx < parameter.size(); idx++)
        if (((a == parameter[idx].a) && (b == parameter[idx].b) && (c == parameter[idx].c)) || 
            ((a == parameter[idx].c) && (b == parameter[idx].b) && (c == parameter[idx].a))) {
          par = &parameter[idx];
          return par;
        }

    for (unsigned int idx=0; idx < parameter.size(); idx++)
      if (((a == parameter[idx].a) && (b == parameter[idx].b) && (c == parameter[idx].c) && (d == parameter[idx].d)) || 
          ((a == parameter[idx].d) && (b == parameter[idx].c) && (c == parameter[idx].b) && (d == parameter[idx].a))) {
        par = &parameter[idx];
        return par;
      }

    return NULL;
  }
  
  OBFFParameter* OBForceField::GetParameter(const char* a, const char* b, const char* c, const char* d, vector<OBFFParameter> &parameter)
  {
    OBFFParameter *par;
    if (a == NULL)
      return NULL;

    if (b == NULL) {
      string _a(a);
      for (unsigned int idx=0; idx < parameter.size(); idx++) 
        if (_a == parameter[idx]._a) {
          par = &parameter[idx];
          return par;
        }
      return NULL;
    }
    if (c == NULL) {
      string _a(a);
      string _b(b);
      for (unsigned int idx=0; idx < parameter.size(); idx++) {
        if (((_a == parameter[idx]._a) && (_b == parameter[idx]._b)) || ((_a == parameter[idx]._b) && (_b == parameter[idx]._a))) {
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
        if (((_a == parameter[idx]._a) && (_b == parameter[idx]._b) && (_c == parameter[idx]._c)) || 
            ((_a == parameter[idx]._c) && (_b == parameter[idx]._b) && (_c == parameter[idx]._a))) {
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
      if (((_a == parameter[idx]._a) && (_b == parameter[idx]._b) && (_c == parameter[idx]._c) && (_d == parameter[idx]._d)) || 
          ((_a == parameter[idx]._d) && (_b == parameter[idx]._c) && (_c == parameter[idx]._b) && (_d == parameter[idx]._a))) {
        par = &parameter[idx];
        return par;
      }

    return NULL;
  }
  
  bool OBForceField::SetLogFile(ostream* pos)
  {
    if(pos)
      logos = pos;
    else
      logos = &cout;
    
    return true;
  }
  
  bool OBForceField::SetLogLevel(int level)
  {
    loglvl = level; 

    return true;
  }
 
  void OBForceField::SystematicRotorSearch() 
  {
    OBRotorList rl;
    OBRotamerList rotamers;

    rl.Setup(_mol);
    rotamers.SetBaseCoordinateSets(_mol);
    rotamers.Setup(_mol, rl);
    
    IF_OBFF_LOGLVL_LOW {
      OBFFLog("\nS Y S T E M A T I C   R O T O R   S E A R C H\n\n");
      sprintf(logbuf, "  NUMBER OF ROTATABLE BONDS: %d\n", rl.Size());
      OBFFLog(logbuf);
    }

    if (!rl.Size()) { // only one conformer
      IF_OBFF_LOGLVL_LOW
        OBFFLog("  GENERATED ONLY ONE CONFORMER\n\n");
 
      ConjugateGradients(2500); // final energy minimizatin for best conformation
      
      return;
    }

    std::vector<int> rotorKey(rl.Size() + 1, 0); // indexed from 1

    OBRotorIterator ri;
    OBRotor *rotor = rl.BeginRotor(ri);
    for (int i = 1; i < rl.Size() + 1; ++i, rotor = rl.NextRotor(ri)) { // foreach rotor
      for (unsigned int j = 0; j < rotor->GetResolution().size(); j++) { // foreach torsion
        rotorKey[i] = j;
        rotamers.AddRotamer(rotorKey);
      }
    }

    rotamers.ExpandConformerList(_mol, _mol.GetConformers());
      
    IF_OBFF_LOGLVL_LOW {
      sprintf(logbuf, "  GENERATED %d CONFORMERS\n\n", _mol.NumConformers());
      OBFFLog(logbuf);
      OBFFLog("CONFORMER     ENERGY\n");
      OBFFLog("--------------------\n");
    }
    
    // Calculate energy for all conformers
    char logbuf[100];
    std::vector<double> energies(_mol.NumConformers(), 0.0);
    for (int i = 0; i < _mol.NumConformers(); i++) {
      _mol.SetConformer(i); // select conformer

      energies[i] = Energy(); // calculate and store energy
      
      IF_OBFF_LOGLVL_LOW {
        sprintf(logbuf, "   %3d      %8.3f\n", (i + 1), energies[i]);
        OBFFLog(logbuf);
      }
    }

    // Select conformer with lowest energy
    int best_conformer = 0;
    for (int i = 1; i < _mol.NumConformers(); i++) {
      if (energies[i] < energies[best_conformer])
        best_conformer = i;
    }
  
    IF_OBFF_LOGLVL_LOW {
      sprintf(logbuf, "\n  CONFORMER %d HAS THE LOWEST ENERGY\n\n",  best_conformer + 1);
      OBFFLog(logbuf);
    }

    _mol.SetConformer(best_conformer);
    current_conformer = best_conformer;
  }

  // LineSearch 
  //
  // atom: coordinates of atom at iteration k (x_k)
  // direction: search direction ( d = -grad(x_0) )
  //
  // ALGORITHM:
  // 
  // step = 1
  // for (i = 1 to 100) {                max steps = 100
  //   e_k = energy(x_k)                 energy of current iteration
  //   x_k = x_k + step * d              update coordinates
  //   e_k+1 = energy(x_k+1)             energy of next iteration
  //   
  //   if (e_k+1 < e_k)
  //     step = step * 1.2               increase step size
  //   if (e_k+1 > e_k) {
  //     x_k = x_k - step * d            reset coordinates to previous iteration
  //     step = step * 0.5               reduce step size
  //   }
  //   if (e_k+1 == e_k)
  //     end                             convergence criteria reached, stop
  // }
  vector3 OBForceField::LineSearch(OBAtom *atom, vector3 &direction)
  {
    double e_n1, e_n2, step;
    vector3 old_xyz, orig_xyz, xyz_k, dir(0.0, 0.0, 0.0);

    step = 0.2;
    direction.normalize();
    orig_xyz = atom->GetVector();
    
    e_n1 = Energy(false); // calculate e_k
    
    unsigned int i;
    for (i=0; i<100; i++) {
      old_xyz = atom->GetVector();
      
      xyz_k = atom->GetVector() + direction*step;
      atom->SetVector(xyz_k);  // update coordinates
    
      e_n2 = Energy(false); // calculate e_k+1
      
      // convergence criteria, this is 10 times the default convergence
      // of SteepestDescent or ConjugateGradients. A higher precision here 
      // only takes longer with the same result.
      if (IsNear(e_n2, e_n1, 1.0e-7))
        break;

      if (e_n2 > e_n1) { // decrease stepsize
        step *= 0.5;
        atom->SetVector(old_xyz);
      }
      if (e_n2 < e_n1) {  // increase stepsize
        e_n1 = e_n2;
        step *= 1.2;
        if (step > 1.0)
          step = 1.0;
      }
      
    }
    //cout << "LineSearch steps: " << i << endl;

    dir = atom->GetVector() - orig_xyz;
    atom->SetVector(orig_xyz);     

    // cutoff accuracy
    if (dir.length() < 1.0e-8)
      return VZero;

    return dir;
  }
  
  void OBForceField::SteepestDescentInitialize(int steps, double econv, int method) 
  {
    _nsteps = steps;
    _cstep = 0;
    _econv = econv;
    _method = method;

    _e_n1 = Energy();
    
    IF_OBFF_LOGLVL_LOW {
      OBFFLog("\nS T E E P E S T   D E S C E N T\n\n");
      sprintf(logbuf, "STEPS = %d\n\n",  steps);
      OBFFLog(logbuf);
      OBFFLog("STEP n     E(n)       E(n-1)    \n");
      OBFFLog("--------------------------------\n");
    }
  }
 
  bool OBForceField::SteepestDescentTakeNSteps(int n) 
  {
    double e_n2;
    vector3 grad;

    for (int i = 1; i <= n; i++) {
      _cstep++;
      
      FOR_ATOMS_OF_MOL (a, _mol) {
        if (_method & OBFF_ANALYTICAL_GRADIENT)
          grad = GetGradient(&*a);
        else
          grad = NumericalDerivative(&*a);
        grad = LineSearch(&*a, grad);
        a->SetVector(a->x() + grad.x(), a->y() + grad.y(), a->z() + grad.z());
      }
      e_n2 = Energy();
      
      IF_OBFF_LOGLVL_LOW {
        sprintf(logbuf, " %4d    %8.3f    %8.3f\n", i, e_n2, _e_n1);
        OBFFLog(logbuf);
      }

      if (IsNear(e_n2, _e_n1, _econv)) {
        IF_OBFF_LOGLVL_LOW
          OBFFLog("    STEEPEST DESCENT HAS CONVERGED\n");
        return false;
      }
      
      if (_nsteps == _cstep)
        return false;

      _e_n1 = e_n2;
    }

    return true;  // no convergence reached
  }
 
  void OBForceField::SteepestDescent(int steps, double econv, int method) 
  {
    SteepestDescentInitialize(steps, econv, method);
    SteepestDescentTakeNSteps(steps);
  }
  
  void OBForceField::ConjugateGradientsInitialize(int steps, double econv, 
                                                  int method)
  {
    double e_n2;
    vector3 grad2, dir2;

    _cstep = 1;
    _nsteps = steps;
    _econv = econv;
    _method = method;

    ValidateGradients();

    _e_n1 = Energy();
    
    IF_OBFF_LOGLVL_LOW {
      OBFFLog("\nC O N J U G A T E   G R A D I E N T S\n\n");
      sprintf(logbuf, "STEPS = %d\n\n",  steps);
      OBFFLog(logbuf);
      OBFFLog("STEP n     E(n)       E(n-1)    \n");
      OBFFLog("--------------------------------\n");
    }

    _grad1.resize(_mol.NumAtoms() + 1);
    _dir1.resize(_mol.NumAtoms() + 1);

    // Take the first step (same as steepest descent because there is no 
    // gradient from the previous step.
    FOR_ATOMS_OF_MOL (a, _mol) {
      if (_method & OBFF_ANALYTICAL_GRADIENT)
        grad2 = GetGradient(&*a);
      else
        grad2 = NumericalDerivative(&*a);
      dir2 = grad2;
      dir2 = LineSearch(&*a, dir2);
      a->SetVector(a->x() + dir2.x(), a->y() + dir2.y(), a->z() + dir2.z());
      _grad1[a->GetIdx()] = grad2;
      _dir1[a->GetIdx()] = grad2;
    }
    e_n2 = Energy();
      
    IF_OBFF_LOGLVL_LOW {
      sprintf(logbuf, " %4d    %8.3f    %8.3f\n", _cstep, e_n2, _e_n1);
      OBFFLog(logbuf);
    }
 
    _e_n1 = e_n2;
  }
  
  bool OBForceField::ConjugateGradientsTakeNSteps(int n)
  {
    double e_n2;
    double g2g2, g1g1, g2g1;
    vector3 grad2, dir2;
    
    if (_grad1.size() != (_mol.NumAtoms()+1))
      return false;

    e_n2 = 0.0;
    
    for (int i = 1; i <= n; i++) {
      _cstep++;
      
      FOR_ATOMS_OF_MOL (a, _mol) {
        if (_method & OBFF_ANALYTICAL_GRADIENT)
          grad2 = GetGradient(&*a);
        else
          grad2 = NumericalDerivative(&*a);
        
        g2g2 = dot(grad2, grad2);
        g1g1 = dot(_grad1[a->GetIdx()], _grad1[a->GetIdx()]);
        g2g1 = g2g2 / g1g1;
        dir2 = grad2 + g2g1 * _dir1[a->GetIdx()];

        dir2 = LineSearch(&*a, dir2);
        a->SetVector(a->x() + dir2.x(), a->y() + dir2.y(), a->z() + dir2.z());
	  
        _grad1[a->GetIdx()] = grad2;
        _dir1[a->GetIdx()] = dir2;
        if (e_n2)
          _e_n1 = e_n2;
      }
      e_n2 = Energy();
	
      IF_OBFF_LOGLVL_LOW {
        sprintf(logbuf, " %4d    %8.3f    %8.3f\n", _cstep, e_n2, _e_n1);
        OBFFLog(logbuf);
      }
 
      if (IsNear(e_n2, _e_n1, _econv)) {
        IF_OBFF_LOGLVL_LOW
          OBFFLog("    CONJUGATE GRADIENTS HAS CONVERGED\n");
        return false;
      }

      if (_nsteps == _cstep)
        return false;

      _e_n1 = e_n2;
    }

    return true; // no convergence reached
  }
 
  void OBForceField::ConjugateGradients(int steps, double econv, int method)
  {
    ConjugateGradientsInitialize(steps, econv, method);
    ConjugateGradientsTakeNSteps(steps); // ConjugateGradientsInitialize takes the fisrt step
  }
  
  //  
  //         f(1) - f(0)
  // f'(0) = -----------      f(1) = f(0+h)
  //              h
  //
  vector3 OBForceField::NumericalDerivative(OBAtom *atom, int terms)
  {
    vector3 va, grad;
    double e_orig, e_plus_delta, delta, dx, dy, dz;

    delta = 1.0e-5;

    va = atom->GetVector();

    if (terms & OBFF_ENERGY)
      e_orig = Energy(false);
    else {
      e_orig = 0.0;
      if (terms & OBFF_EBOND)
        e_orig += E_Bond();
      if (terms & OBFF_EANGLE)
        e_orig += E_Angle();
      if (terms & OBFF_ESTRBND)
        e_orig += E_StrBnd();
      if (terms & OBFF_ETORSION)
        e_orig += E_Torsion();
      if (terms & OBFF_EOOP)
        e_orig += E_OOP();
      if (terms & OBFF_EVDW)
        e_orig += E_VDW();
      if (terms & OBFF_EELECTROSTATIC)
        e_orig += E_Electrostatic();
    }
    
    // X direction
    atom->SetVector(va.x() + delta, va.y(), va.z());

    if (terms & OBFF_ENERGY)
      e_plus_delta = Energy(false);
    else {
      e_plus_delta = 0.0;
      if (terms & OBFF_EBOND)
        e_plus_delta += E_Bond();
      if (terms & OBFF_EANGLE)
        e_plus_delta += E_Angle();
      if (terms & OBFF_ESTRBND)
        e_plus_delta += E_StrBnd();
      if (terms & OBFF_ETORSION)
        e_plus_delta += E_Torsion();
      if (terms & OBFF_EOOP)
        e_plus_delta += E_OOP();
      if (terms & OBFF_EVDW)
        e_plus_delta += E_VDW();
      if (terms & OBFF_EELECTROSTATIC)
        e_plus_delta += E_Electrostatic();
    }
    
    dx = (e_plus_delta - e_orig) / delta;
    
    // Y direction
    atom->SetVector(va.x(), va.y() + delta, va.z());
    
    if (terms & OBFF_ENERGY)
      e_plus_delta = Energy(false);
    else {
      e_plus_delta = 0.0;
      if (terms & OBFF_EBOND)
        e_plus_delta += E_Bond();
      if (terms & OBFF_EANGLE)
        e_plus_delta += E_Angle();
      if (terms & OBFF_ESTRBND)
        e_plus_delta += E_StrBnd();
      if (terms & OBFF_ETORSION)
        e_plus_delta += E_Torsion();
      if (terms & OBFF_EOOP)
        e_plus_delta += E_OOP();
      if (terms & OBFF_EVDW)
        e_plus_delta += E_VDW();
      if (terms & OBFF_EELECTROSTATIC)
        e_plus_delta += E_Electrostatic();
    }
 
    dy = (e_plus_delta - e_orig) / delta;
    
    // Z direction
    atom->SetVector(va.x(), va.y(), va.z() + delta);
    
    if (terms & OBFF_ENERGY)
      e_plus_delta = Energy(false);
    else {
      e_plus_delta = 0.0;
      if (terms & OBFF_EBOND)
        e_plus_delta += E_Bond();
      if (terms & OBFF_EANGLE)
        e_plus_delta += E_Angle();
      if (terms & OBFF_ESTRBND)
        e_plus_delta += E_StrBnd();
      if (terms & OBFF_ETORSION)
        e_plus_delta += E_Torsion();
      if (terms & OBFF_EOOP)
        e_plus_delta += E_OOP();
      if (terms & OBFF_EVDW)
        e_plus_delta += E_VDW();
      if (terms & OBFF_EELECTROSTATIC)
        e_plus_delta += E_Electrostatic();
    }
 
    dz = (e_plus_delta - e_orig) / delta;

    // reset coordinates to original
    atom->SetVector(va.x(), va.y(), va.z());

    grad.Set(-dx, -dy, -dz);
    return (grad);
  }
  
  //  
  //         f(2) - 2f(1) + f(0)
  // f'(0) = -------------------      f(1) = f(0+h)
  //                 h^2              f(1) = f(0+2h)
  //
  vector3 OBForceField::NumericalSecondDerivative(OBAtom *atom, int terms)
  {
    vector3 va, grad;
    double e_0, e_1, e_2, delta, dx, dy, dz;

    delta = 1.0e-5;

    va = atom->GetVector();

    // calculate f(0)
    if (terms & OBFF_ENERGY)
      e_0 = Energy(false);
    else {
      e_0 = 0.0;
      if (terms & OBFF_EBOND)
        e_0 += E_Bond(false);
      if (terms & OBFF_EANGLE)
        e_0 += E_Angle(false);
      if (terms & OBFF_ESTRBND)
        e_0 += E_StrBnd(false);
      if (terms & OBFF_ETORSION)
        e_0 += E_Torsion(false);
      if (terms & OBFF_EOOP)
        e_0 += E_OOP(false);
      if (terms & OBFF_EVDW)
        e_0 += E_VDW(false);
      if (terms & OBFF_EELECTROSTATIC)
        e_0 += E_Electrostatic(false);
    }
    
    // 
    // X direction
    //
    
    // calculate f(1)
    atom->SetVector(va.x() + delta, va.y(), va.z());

    if (terms & OBFF_ENERGY)
      e_1 = Energy(false);
    else {
      e_1 = 0.0;
      if (terms & OBFF_EBOND)
        e_1 += E_Bond(false);
      if (terms & OBFF_EANGLE)
        e_1 += E_Angle(false);
      if (terms & OBFF_ESTRBND)
        e_1 += E_StrBnd(false);
      if (terms & OBFF_ETORSION)
        e_1 += E_Torsion(false);
      if (terms & OBFF_EOOP)
        e_1 += E_OOP(false);
      if (terms & OBFF_EVDW)
        e_1 += E_VDW(false);
      if (terms & OBFF_EELECTROSTATIC)
        e_1 += E_Electrostatic(false);
    }
    
    // calculate f(2)
    atom->SetVector(va.x() + 2 * delta, va.y(), va.z());

    if (terms & OBFF_ENERGY)
      e_2 = Energy(false);
    else {
      e_2 = 0.0;
      if (terms & OBFF_EBOND)
        e_2 += E_Bond(false);
      if (terms & OBFF_EANGLE)
        e_2 += E_Angle(false);
      if (terms & OBFF_ESTRBND)
        e_2 += E_StrBnd(false);
      if (terms & OBFF_ETORSION)
        e_2 += E_Torsion(false);
      if (terms & OBFF_EOOP)
        e_2 += E_OOP(false);
      if (terms & OBFF_EVDW)
        e_2 += E_VDW(false);
      if (terms & OBFF_EELECTROSTATIC)
        e_2 += E_Electrostatic(false);
    }
    
    dx = (e_2 - 2 * e_1 + e_0) / (delta * delta);
    
    // 
    // Y direction
    //
    
    // calculate f(1)
    atom->SetVector(va.x(), va.y() + delta, va.z());

    if (terms & OBFF_ENERGY)
      e_1 = Energy(false);
    else {
      e_1 = 0.0;
      if (terms & OBFF_EBOND)
        e_1 += E_Bond(false);
      if (terms & OBFF_EANGLE)
        e_1 += E_Angle(false);
      if (terms & OBFF_ESTRBND)
        e_1 += E_StrBnd(false);
      if (terms & OBFF_ETORSION)
        e_1 += E_Torsion(false);
      if (terms & OBFF_EOOP)
        e_1 += E_OOP(false);
      if (terms & OBFF_EVDW)
        e_1 += E_VDW(false);
      if (terms & OBFF_EELECTROSTATIC)
        e_1 += E_Electrostatic(false);
    }
    
    // calculate f(2)
    atom->SetVector(va.x(), va.y() + 2 * delta, va.z());

    if (terms & OBFF_ENERGY)
      e_2 = Energy(false);
    else {
      e_2 = 0.0;
      if (terms & OBFF_EBOND)
        e_2 += E_Bond(false);
      if (terms & OBFF_EANGLE)
        e_2 += E_Angle(false);
      if (terms & OBFF_ESTRBND)
        e_2 += E_StrBnd(false);
      if (terms & OBFF_ETORSION)
        e_2 += E_Torsion(false);
      if (terms & OBFF_EOOP)
        e_2 += E_OOP(false);
      if (terms & OBFF_EVDW)
        e_2 += E_VDW(false);
      if (terms & OBFF_EELECTROSTATIC)
        e_2 += E_Electrostatic(false);
    }
    
    dy = (e_2 - 2 * e_1 + e_0) / (delta * delta);

    // 
    // Z direction
    //
    
    // calculate f(1)
    atom->SetVector(va.x(), va.y(), va.z() + delta);

    if (terms & OBFF_ENERGY)
      e_1 = Energy(false);
    else {
      e_1 = 0.0;
      if (terms & OBFF_EBOND)
        e_1 += E_Bond(false);
      if (terms & OBFF_EANGLE)
        e_1 += E_Angle(false);
      if (terms & OBFF_ESTRBND)
        e_1 += E_StrBnd(false);
      if (terms & OBFF_ETORSION)
        e_1 += E_Torsion(false);
      if (terms & OBFF_EOOP)
        e_1 += E_OOP(false);
      if (terms & OBFF_EVDW)
        e_1 += E_VDW(false);
      if (terms & OBFF_EELECTROSTATIC)
        e_1 += E_Electrostatic(false);
    }
    
    // calculate f(2)
    atom->SetVector(va.x(), va.y(), va.z() + 2 * delta);

    if (terms & OBFF_ENERGY)
      e_2 = Energy(false);
    else {
      e_2 = 0.0;
      if (terms & OBFF_EBOND)
        e_2 += E_Bond(false);
      if (terms & OBFF_EANGLE)
        e_2 += E_Angle(false);
      if (terms & OBFF_ESTRBND)
        e_2 += E_StrBnd(false);
      if (terms & OBFF_ETORSION)
        e_2 += E_Torsion(false);
      if (terms & OBFF_EOOP)
        e_2 += E_OOP(false);
      if (terms & OBFF_EVDW)
        e_2 += E_VDW(false);
      if (terms & OBFF_EELECTROSTATIC)
        e_2 += E_Electrostatic(false);
    }
    
    dz = (e_2 - 2 * e_1 + e_0) / (delta * delta);

    // reset coordinates to original
    atom->SetVector(va.x(), va.y(), va.z());

    grad.Set(-dx, -dy, -dz);
    return (grad);
  }

  bool OBForceField::UpdateCoordinates(OBMol &mol)
  { 
    OBAtom *atom;
    
    if (_mol.NumAtoms() != mol.NumAtoms())
      return false;
    
    // Copy coordinates for current conformer only
    FOR_ATOMS_OF_MOL (a, _mol) {
      atom = mol.GetAtom(a->GetIdx());
      atom->SetVector(a->GetVector());
    }

    return true;
  }
  
  bool OBForceField::UpdateConformers(OBMol &mol)
  { 
    OBAtom *atom;

    if (_mol.NumAtoms() != mol.NumAtoms())
      return false;
    
    FOR_ATOMS_OF_MOL (a, _mol) {
      atom = mol.GetAtom(a->GetIdx());
      atom->SetVector(a->GetVector());
    }

    //Copy conformer information
    if (_mol.NumConformers() > 1) {
      int k,l;
      vector<double*> conf;
      double* xyz = NULL;
      for (k=0 ; k<_mol.NumConformers() ; ++k) {
        xyz = new double [3*_mol.NumAtoms()];
        for (l=0 ; l<(int) (3*_mol.NumAtoms()) ; ++l)
          xyz[l] = _mol.GetConformer(k)[l];
        conf.push_back(xyz);
      }
      mol.SetConformers(conf);
      mol.SetConformer(current_conformer);
    }
    
    return true;
  }

  vector3 OBForceField::ValidateGradientError(vector3 &numgrad, vector3 &anagrad)
  {
    double errx, erry, errz;

    if (fabs(numgrad.x()) < 1.0)
      errx = numgrad.x() * fabs(numgrad.x() - anagrad.x()) * 100;
    else
      errx = fabs((numgrad.x() - anagrad.x()) / numgrad.x()) * 100;

    if (fabs(numgrad.y()) < 1.0)
      erry = numgrad.y() * fabs(numgrad.y() - anagrad.y()) * 100;
    else
      erry = fabs((numgrad.y() - anagrad.y()) / numgrad.y()) * 100;
    
    if (fabs(numgrad.z()) < 1.0)
      errz = numgrad.z() * fabs(numgrad.z() - anagrad.z()) * 100;
    else
      errz = fabs((numgrad.z() - anagrad.z()) / numgrad.z()) * 100;
    
    errx = fabs(errx);
    erry = fabs(erry);
    errz = fabs(errz);

    return vector3(errx, erry, errz);
  }
  
  double OBForceField::VectorLengthDerivative(vector3 &a, vector3 &b)
  {
    vector3 vab, drab;
    double rab;
    
    vab = a - b;
    rab = vab.length();
    if (rab < 0.1) // atoms are too close to each other
      {
        vab.randomUnitVector();
        vab *= 0.1; // move the atoms a small distance apart
        rab = 0.1;
      }
    drab = vab / rab;

    a = -drab; // -drab/da
    b =  drab; // -drab/db

    return rab;
  }
  
  double OBForceField::VectorAngleDerivative(vector3 &a, vector3 &b, vector3 &c)
  {
    vector3 vab, vcb;
    double theta, rab, rab2, rcb, rcb2, abcb, abcb2;
     
    vab = a - b;
    vcb = c - b;
    rab = vab.length();
    rcb = vcb.length();
    rab2 = rab * rab;
    rcb2 = rcb * rcb;
    
    abcb = dot(vab, vcb) / (rab * rcb);
    abcb2 = 1.0 - abcb * abcb;
    theta = acos(abcb) * RAD_TO_DEG;
    
    if (IsNearZero(abcb2)) {
      a = VZero;
      b = VZero;
      c = VZero;
      return 0.0;
    }

    a = (vcb * rab * rcb - (vab / rab) * dot(vab, vcb) * rcb) / (sqrt(abcb2) * rab2 * rcb2);
    c = -((vcb / rcb) * dot(vab, vcb) * rab - vab * rab * rcb) / (sqrt(abcb2) * rab2 * rcb2);
    b = -a - c;

    a *= (1.0 / DEG_TO_RAD);
    b *= (1.0 / DEG_TO_RAD);
    c *= (1.0 / DEG_TO_RAD);

    return theta;
  }
 
  double OBForceField::VectorTorsionDerivative(vector3 &a, vector3 &b, vector3 &c, vector3 &d)
  {
    vector3 vab, vbc, vcd, vac, vbd, grada, gradb, gradc, gradd;
    vector3 abbc, bccd, bcabbc, bcbccd, cdabbc, cdbccd, acabbc, acbccd, ababbc, abbccd, bdabbc, bdbccd;
    double tor, rabbc, rbccd, rabbc2, rbccd2, rabbc3, rbccd3, abbc_bccd, rabbc_rbccd, abbc_bccd2;

    vab = a - b;
    vbc = b - c;
    vcd = c - d;
    vac = a - c;
    vbd = b - d;
    abbc = cross(vab, vbc);
    bccd = cross(vbc, vcd);
    
    tor = RAD_TO_DEG * acos(dot(abbc, bccd) / (abbc.length() * bccd.length()));
    if (dot(abbc, bccd) > 0.0)
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
    abbc_bccd = dot(abbc, bccd);
    rabbc_rbccd = abbc_bccd / (rabbc * rbccd);
    abbc_bccd2 = sqrt(1.0 - rabbc_rbccd * rabbc_rbccd);

    if (IsNearZero(abbc_bccd2))
      abbc_bccd2 = 1e-7;
    
    /*
      cout << "sqrt(abbc_bccd2) = " << sqrt(abbc_bccd2) << endl;
      cout << "(rabbc*rbccd) = " << (rabbc*rbccd) << endl;
      cout << "dot(abbc,bccd) = " << dot(abbc,bccd) << endl;
      cout << "(rabbc*rbccd3) = " << (rabbc*rbccd3) << endl;
      cout << "(rabbc3*rbccd) = " << (rabbc3*rbccd) << endl;
    */
    
    a = (bcbccd / (rabbc*rbccd) - (bcabbc*abbc_bccd) / (rabbc3*rbccd)) / abbc_bccd2;
    d = (bcabbc / (rabbc*rbccd) - (bcbccd*abbc_bccd) / (rabbc*rbccd3)) / abbc_bccd2;

    b = ( -(cdbccd*abbc_bccd) / (rabbc*rbccd3) + 
          (cdabbc - acbccd) / (rabbc*rbccd) + 
          (acabbc*abbc_bccd) / (rabbc3*rbccd)    ) / abbc_bccd2;
    
    c = (  (bdbccd*abbc_bccd) / (rabbc*rbccd3) + 
           (abbccd - bdabbc) / (rabbc*rbccd) +
           -(ababbc*abbc_bccd) / (rabbc3*rbccd)    ) / abbc_bccd2;

    if (IsNan(a.x()) || IsNan(a.y()) || IsNan(a.z()))
      a = VZero;
    if (IsNan(b.x()) || IsNan(b.y()) || IsNan(b.z()))
      b = VZero;
    if (IsNan(c.x()) || IsNan(c.y()) || IsNan(c.z()))
      c = VZero;
    if (IsNan(d.x()) || IsNan(d.y()) || IsNan(d.z()))
      d = VZero;

    if (dot(abbc, bccd) > 0.0) {
      a = -a;
      b = -b;
      c = -c;
      d = -d;
    }
   
    return tor;  
  }
  
  bool OBForceField::IsInSameRing(OBAtom* a, OBAtom* b)
  {
    bool a_in, b_in;
    vector<OBRing*> vr;
    vr = _mol.GetSSSR();
    
    vector<OBRing*>::iterator i;
    vector<int>::iterator j;
    
    for (i = vr.begin();i != vr.end();i++) {
      a_in = false;
      b_in = false;
      for(j = (*i)->_path.begin();j != (*i)->_path.end();j++) {
        if ((unsigned)(*j) == a->GetIdx())
          a_in = true;
        if ((unsigned)(*j) == b->GetIdx())
          b_in = true;
      }
      
      if (a_in && b_in)
        return true;
    }

    return false;
  }
 
} // end namespace OpenBabel


//! \file forcefield.cpp
//! \brief Handle OBForceField class

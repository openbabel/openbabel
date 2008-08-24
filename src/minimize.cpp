/**********************************************************************
minimize.cpp - Handle OBMinimize class.

Copyright (C) 2006-2007 by Tim Vandermeersch <tim.vandermeersch@gmail.com>
Some portions Copyright (C) 2007-2008 by Geoffrey Hutchison
 
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

#include <openbabel/minimize.h>

using namespace std;

namespace OpenBabel
{
  class OBMinimizePrivate 
  {
    public:
    OBForceField* forcefield;
    // general variables
    bool 	init; //!< Used to make sure we only parse the parameter file once, when needed
    // minimization variables
    double 	econv, e_n1; //!< Used for conjugate gradients and steepest descent(Initialize and TakeNSteps)
    int 	cstep, nsteps; //!< Used for conjugate gradients and steepest descent(Initialize and TakeNSteps)
    std::vector<Eigen::Vector3d> grad1; //!< Used for conjugate gradients and steepest descent(Initialize and TakeNSteps)
    unsigned int nAtoms; //!< Number of atoms
    int         linesearch; //!< LineSearch type

    char        logbuf[BUFF_SIZE];
  };

  OBMinimize::OBMinimize(OBForceField *forcefield) : d(new OBMinimizePrivate)
  {
    d->forcefield = forcefield;
  }
  
  OBMinimize::~OBMinimize()
  {
    delete d;
  }

  //////////////////////////////////////////////////////////////////////////////////
  //
  // Energy Minimization  (XXXX07)
  //
  //////////////////////////////////////////////////////////////////////////////////
  
  void OBMinimize::SetLineSearchType(int type)
  {
    d->linesearch = type;
  }
    
  int OBMinimize::GetLineSearchType()
  {
    return d->linesearch;
  }
 
  // LineSearch 
  //
  // Based on the ghemical code (conjgrad.cpp)
  //
  // Implements several enhancements:
  // 1) Switch to smarter line search method (e.g., Newton's method in 1D)
  //  x(n+1) = x(n) - F(x) / F'(x)   -- can be done numerically
  // 2) Switch to line search for one step of the entire molecule!
  //   This dramatically cuts down on the number of Energy() calls.
  //  (and is more correct anyway)
  //double OBMinimize::Newton2NumLineSearch(double *direction)
  /*
  double OBMinimize::Newton2NumLineSearch(std::vector<Eigen::Vector3d> &direction)
  {
    double e_n1, e_n2, e_n3;
    //double *origCoords = new double [d->ncoords];
    vector<Eigen::Vector3d> origCoords;

    double opt_step = 0.0;
    double opt_e = d->e_n1; // get energy calculated by sd or cg
    const double def_step = 0.025; // default step
    const double max_step = 5.0; // don't move further than 0.3 Angstroms
    
    double sum = 0.0;
    for (unsigned int c = 0; c < d->mol.NumAtoms(); ++c) {
      if (isfinite( direction[c].norm2() )) { 
        sum += direction[c].norm2();
      } else {
        // make sure we don't have NaN or infinity
        direction[c] = Eigen::Vector3d::Zero();
      }
    }

    double scale = sqrt(sum);
    if (IsNearZero(scale)) {
      cout << "WARNING: too small \"scale\" at Newton2NumLineSearch" << endl;
      scale = 1.0e-70; // try to avoid "division by zero" conditions
    }

    double step = def_step / scale;
    double max_scl = max_step / scale;
    
    // Save the current position, before we take a step
    memcpy((char*)origCoords,(char*)d->mol.GetCoordinates(),sizeof(double)*d->ncoords);
    
    int newton = 0;
    while (true) {
      // Take step X(n) + step
      LineSearchTakeStep(origCoords, direction, step);
      e_n1 = Energy(false) + d->constraints.GetConstraintEnergy();

      if (e_n1 < opt_e) {
        opt_step = step;
        opt_e = e_n1;
      }

      if (newton++ > 3) 
        break;
      double delta = step * 0.001;
      
      // Take step X(n) + step + delta
      LineSearchTakeStep(origCoords, direction, step+delta);
      e_n2 = Energy(false) + d->constraints.GetConstraintEnergy();
 
      // Take step X(n) + step + delta * 2.0
      LineSearchTakeStep(origCoords, direction, step+delta*2.0);
      e_n3 = Energy(false) + d->constraints.GetConstraintEnergy();
      
      double denom = e_n3 - 2.0 * e_n2 + e_n1; // f'(x)
      if (denom != 0.0) {
        step = fabs(step - delta * (e_n2 - e_n1) / denom);
        if (step > max_scl) {
          cout << "WARNING: damped steplength " << step << " to " << max_scl << endl;
          step = max_scl;
        }
      } else {
        break;
      }
    }
    
    if (opt_step == 0.0) { // if we still don't have any valid steplength, try a very small step
      step = 0.001 * def_step / scale;
      
      // Take step X(n) + step
      LineSearchTakeStep(origCoords, direction, step);
      e_n1 = Energy(false) + d->constraints.GetConstraintEnergy();

      if (e_n1 < opt_e) {
        opt_step = step;
        opt_e = e_n1;
      }
      
    }

    // Take optimal step 
    LineSearchTakeStep(origCoords, direction, opt_step);

    return opt_step * scale;
  }
  
  void OBMinimize::LineSearchTakeStep(double* origCoords, double *direction, double step)
  {
    double *currentCoords = d->mol.GetCoordinates();

    for (unsigned int c = 0; c < d->ncoords; ++c) {
      if (isfinite(direction[c])) { 
        currentCoords[c] = origCoords[c] + direction[c] * step;
      }
    }
  }
  */

  double OBMinimize::LineSearch(std::vector<Eigen::Vector3d> &currentCoords, std::vector<Eigen::Vector3d> &direction)
  {
    double e_n1, e_n2, step, alpha;//, tempStep;
    Eigen::Vector3d tempStep;
    vector<Eigen::Vector3d> lastStep;

    alpha = 0.0; // Scale factor along direction vector
    step = 0.2;
    double trustRadius = 0.3; // don't move further than 0.3 Angstroms
    double trustRadius2 = 0.9; // use norm2() instead of norm() to avoid sqrt() calls
    
    e_n1 = d->forcefield->Energy(false);
    
    unsigned int i;
    for (i=0; i < 10; ++i) {
      // Save the current position, before we take a step
      //memcpy((char*)lastStep,(char*)currentCoords,sizeof(double)*numCoords);
      lastStep = currentCoords;
      
      // Vectorizing this would be a big benefit
      // Need to look up using BLAS or Eigen or whatever
      for (unsigned int c = 0; c < currentCoords.size(); ++c) {
        if (isfinite(direction[c].norm2())) { 
          // make sure we don't have NaN or infinity
          tempStep = direction[c] * step;

          if (tempStep.norm2() > trustRadius2) { // big step
            tempStep.normalize();
            tempStep *= trustRadius;
          }
          
          currentCoords[c] += tempStep;
        }
      }
    
      e_n2 = d->forcefield->Energy(false);
      
      // convergence criteria: A higher precision here 
      // only takes longer with the same result.
      if (IsNear(e_n2, e_n1, 1.0e-3))
        break;

      if (e_n2 > e_n1) { // decrease stepsize
        step *= 0.1;
        // move back to the last step
        currentCoords = lastStep;
      } else if (e_n2 < e_n1) {  // increase stepsize
        e_n1 = e_n2;
        alpha += step; // we've moved some distance
        step *= 2.15;
        if (step > 1.0)
          step = 1.0;
      }
      
    }
    //cout << "LineSearch steps: " << i << endl;

    //delete [] lastStep;

    return alpha;
  }
  
  void OBMinimize::SteepestDescentInitialize(int steps, double econv) 
  {
    d->nsteps = steps;
    d->cstep = 0;
    d->econv = econv;

    d->e_n1 = d->forcefield->Energy();
    
    if (d->forcefield->GetLogLevel() >= OBFF_LOGLVL_LOW) {
      d->forcefield->OBFFLog("\nS T E E P E S T   D E S C E N T\n\n");
      snprintf(d->logbuf, BUFF_SIZE, "STEPS = %d\n\n",  steps);
      d->forcefield->OBFFLog(d->logbuf);
      d->forcefield->OBFFLog("STEP n       E(n)         E(n-1)    \n");
      d->forcefield->OBFFLog("------------------------------------\n");
      snprintf(d->logbuf, BUFF_SIZE, " %4d    %8.3f      ----\n", d->cstep, d->e_n1);
      d->forcefield->OBFFLog(d->logbuf);
    }
 
  }
 
  bool OBMinimize::SteepestDescentTakeNSteps(int n) 
  {
    double e_n2, alpha;
    for (int i = 1; i <= n; i++) {
      d->cstep++;

      if (!d->forcefield->HasAnalyticalGradients()) {
        for (unsigned int idx = 0; idx < d->forcefield->GetPositions().size(); ++idx) {
          // use numerical gradients
          d->forcefield->GetGradients()[idx] = NumericalDerivative(idx);
        } 
      }
      
      // perform a linesearch
      switch (d->linesearch) {
        case LineSearchType::Newton2Num:
          //alpha = Newton2NumLineSearch(GetGradients());
          break;
        default:
        case LineSearchType::Simple:
          alpha = LineSearch(d->forcefield->GetPositions(), d->forcefield->GetGradients());
          break;
      }
      e_n2 = d->forcefield->Energy();
      
      if (d->forcefield->GetLogLevel() >= OBFF_LOGLVL_LOW) {
        if (d->cstep % 10 == 0) {
          snprintf(d->logbuf, BUFF_SIZE, " %4d    %8.5f    %8.5f\n", d->cstep, e_n2, d->e_n1);
          d->forcefield->OBFFLog(d->logbuf);
        }
      }

      if (IsNear(e_n2, d->e_n1, d->econv)) {
        if (d->forcefield->GetLogLevel() >= OBFF_LOGLVL_LOW)
          d->forcefield->OBFFLog("    STEEPEST DESCENT HAS CONVERGED\n");
        return false;
      }
      
      if (d->nsteps == d->cstep) {
        return false;
      }

      d->e_n1 = e_n2;
    }

    return true;  // no convergence reached
  }
 
  void OBMinimize::SteepestDescent(int steps, double econv) 
  {
    SteepestDescentInitialize(steps, econv);
    SteepestDescentTakeNSteps(steps);
  }

  void OBMinimize::ConjugateGradientsInitialize(int steps, double econv)
  {
    double e_n2, alpha;

    d->cstep = 0;
    d->nsteps = steps;
    d->econv = econv;

    d->e_n1 = d->forcefield->Energy();
    
    if (d->forcefield->GetLogLevel() >= OBFF_LOGLVL_LOW) {
      d->forcefield->OBFFLog("\nC O N J U G A T E   G R A D I E N T S\n\n");
      snprintf(d->logbuf, BUFF_SIZE, "STEPS = %d\n\n",  steps);
      d->forcefield->OBFFLog(d->logbuf);
      d->forcefield->OBFFLog("STEP n     E(n)       E(n-1)    \n");
      d->forcefield->OBFFLog("--------------------------------\n");
    }

    d->grad1.clear();
    d->grad1.resize(d->forcefield->GetPositions().size());
    for (unsigned int idx = 0; idx < d->forcefield->GetPositions().size(); ++idx)
      d->grad1[idx] = Eigen::Vector3d::Zero();
 
    // Take the first step (same as steepest descent because there is no 
    // gradient from the previous step.
    if (!d->forcefield->HasAnalyticalGradients()) {
      for (unsigned int idx = 0; idx < d->forcefield->GetPositions().size(); ++idx) {
        // use numerical gradients
        d->forcefield->GetGradients()[idx] = NumericalDerivative(idx);
      }
    }
    
    // perform a linesearch
    switch (d->linesearch) {
      case LineSearchType::Newton2Num:
        //alpha = Newton2NumLineSearch(GetGradients());
        break;
      default:
      case LineSearchType::Simple:
        alpha = LineSearch(d->forcefield->GetPositions(), d->forcefield->GetGradients());
        break;
    }
    e_n2 = d->forcefield->Energy();
      
    if (d->forcefield->GetLogLevel() >= OBFF_LOGLVL_LOW) {
      snprintf(d->logbuf, BUFF_SIZE, " %4d    %8.3f    %8.3f\n", 1, e_n2, d->e_n1);
      d->forcefield->OBFFLog(d->logbuf);
    }
 
    // save the direction and energy
    d->grad1 = d->forcefield->GetGradients();
    d->e_n1 = e_n2;
  }
  
  bool OBMinimize::ConjugateGradientsTakeNSteps(int n)
  {
    double e_n2;
    double g2g2, g1g1, beta, alpha;
    Eigen::Vector3d grad2, dir2;
    Eigen::Vector3d grad1, dir1; // temporaries to perform dot product, etc.
    
    e_n2 = 0.0;
    
    for (int i = 1; i <= n; i++) {
      d->cstep++;
     
      for (unsigned int idx = 0; idx < d->forcefield->GetPositions().size(); ++idx) {
          if (!d->forcefield->HasAnalyticalGradients()) {
            // use numerical gradients
            grad2 = NumericalDerivative(idx);
          } else {
            // use analytical gradients
            grad2 = d->forcefield->GetGradients()[idx];
          }
 
          // Fletcher-Reeves formula for Beta
          // http://en.wikipedia.org/wiki/Nonlinear_conjugate_gradient_method
          // NOTE: We make sure to reset and use the steepest descent direction
          //   after NumAtoms steps
          if (d->cstep % d->forcefield->GetPositions().size() != 0) {
            g2g2 = grad2.dot(grad2);
            //grad1 = Eigen::Vector3d(d->grad1[coordIdx], d->grad1[coordIdx+1], d->grad1[coordIdx+2]);
            grad1 = d->grad1[idx];
            g1g1 = grad1.dot(grad1);
            beta = g2g2 / g1g1;
            grad2 += beta * grad1;
          } 
 
          d->grad1[idx] = grad2;
      }
      // perform a linesearch
      switch (d->linesearch) {
        case LineSearchType::Newton2Num:
          //alpha = Newton2NumLineSearch(d->grad1);
          break;
        default:
        case LineSearchType::Simple:
          alpha = LineSearch(d->forcefield->GetPositions(), d->grad1);
          break;
      }
      // save the direction
      d->grad1 = d->forcefield->GetGradients();
 
      e_n2 = d->forcefield->Energy();
	
      if (IsNear(e_n2, d->e_n1, d->econv)) {
        if (d->forcefield->GetLogLevel() >= OBFF_LOGLVL_LOW) {
          snprintf(d->logbuf, BUFF_SIZE, " %4d    %8.3f    %8.3f\n", d->cstep, e_n2, d->e_n1);
          d->forcefield->OBFFLog(d->logbuf);
          d->forcefield->OBFFLog("    CONJUGATE GRADIENTS HAS CONVERGED\n");
        }
        return false;
      }

      if (d->forcefield->GetLogLevel() >= OBFF_LOGLVL_LOW) {
        if (d->cstep % 10 == 0) {
          snprintf(d->logbuf, BUFF_SIZE, " %4d    %8.3f    %8.3f\n", d->cstep, e_n2, d->e_n1);
          d->forcefield->OBFFLog(d->logbuf);
        }
      }
 
      if (d->nsteps == d->cstep)
        return false;

      d->e_n1 = e_n2;
    }

    return true; // no convergence reached
  }
  
  void OBMinimize::ConjugateGradients(int steps, double econv)
  {
    ConjugateGradientsInitialize(steps, econv);
    ConjugateGradientsTakeNSteps(steps); // ConjugateGradientsInitialize takes the first step
  }
  
  //  
  //         f(1) - f(0)
  // f'(0) = -----------      f(1) = f(0+h)
  //              h
  //
  Eigen::Vector3d OBMinimize::NumericalDerivative(unsigned int idx) // start at 0!!
  {
    Eigen::Vector3d va, grad;
    double e_orig, e_plus_delta, delta, dx, dy, dz;

    delta = 1.0e-5;
    va = d->forcefield->GetPositions()[idx];

    e_orig = d->forcefield->Energy(false);
    
    // X direction
    d->forcefield->GetPositions()[idx] = Eigen::Vector3d(va.x() + delta, va.y(), va.z());
    e_plus_delta = d->forcefield->Energy(false);
    dx = (e_plus_delta - e_orig) / delta;
    
    // Y direction
    d->forcefield->GetPositions()[idx] = Eigen::Vector3d(va.x(), va.y() + delta, va.z());
    e_plus_delta = d->forcefield->Energy(false);
    dy = (e_plus_delta - e_orig) / delta;
    
    // Z direction
    d->forcefield->GetPositions()[idx] = Eigen::Vector3d(va.x(), va.y(), va.z() + delta);
    e_plus_delta = d->forcefield->Energy(false);
    dz = (e_plus_delta - e_orig) / delta;

    // reset coordinates to original
    d->forcefield->GetPositions()[idx] = Eigen::Vector3d(va.x(), va.y(), va.z());

    grad = Eigen::Vector3d(-dx, -dy, -dz);
    return (grad);
  }
  
  //  
  //         f(2) - 2f(1) + f(0)
  // f'(0) = -------------------      f(1) = f(0+h)
  //                 h^2              f(1) = f(0+2h)
  //
  Eigen::Vector3d OBMinimize::NumericalSecondDerivative(unsigned int idx)
  {
    Eigen::Vector3d va, grad;
    double e_0, e_1, e_2, delta, dx, dy, dz;

    delta = 1.0e-5;

    va = d->forcefield->GetPositions()[idx];

    // calculate f(0)
    e_0 = d->forcefield->Energy(false);
    
    // 
    // X direction
    //
    
    // calculate f(1)
    d->forcefield->GetPositions()[idx] = Eigen::Vector3d(va.x() + delta, va.y(), va.z());
    e_1 = d->forcefield->Energy(false);
    
    // calculate f(2)
    d->forcefield->GetPositions()[idx] = Eigen::Vector3d(va.x() + 2 * delta, va.y(), va.z());
    e_2 = d->forcefield->Energy(false);
    
    dx = (e_2 - 2 * e_1 + e_0) / (delta * delta);
    
    // 
    // Y direction
    //
    
    // calculate f(1)
    d->forcefield->GetPositions()[idx] = Eigen::Vector3d(va.x(), va.y() + delta, va.z());
    e_1 = d->forcefield->Energy(false);
    
    // calculate f(2)
    d->forcefield->GetPositions()[idx] = Eigen::Vector3d(va.x(), va.y() + 2 * delta, va.z());
    e_2 = d->forcefield->Energy(false);
    
    dy = (e_2 - 2 * e_1 + e_0) / (delta * delta);

    // 
    // Z direction
    //
    
    // calculate f(1)
    d->forcefield->GetPositions()[idx] = Eigen::Vector3d(va.x(), va.y(), va.z() + delta);
    e_1 = d->forcefield->Energy(false);
    
    // calculate f(2)
    d->forcefield->GetPositions()[idx] = Eigen::Vector3d(va.x(), va.y(), va.z() + 2 * delta);
    e_2 = d->forcefield->Energy(false);
    
    dz = (e_2 - 2 * e_1 + e_0) / (delta * delta);

    // reset coordinates to original
    d->forcefield->GetPositions()[idx] = Eigen::Vector3d(va.x(), va.y(), va.z());

    grad = Eigen::Vector3d(-dx, -dy, -dz);
    return (grad);
  }
  
} // end namespace OpenBabel


//! \file forcefield.cpp
//! \brief Handle OBMinimize class

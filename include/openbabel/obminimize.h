/**********************************************************************
obminimize.h - Class to minimize OBFunction.
 
Copyright (C) 2008 by Tim Vandermeersch
 
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

#ifndef OB_MINIMIZE_H
#define OB_MINIMIZE_H

#include <vector>
#include <string>

#include <openbabel/babelconfig.h>
#include <openbabel/obforcefield.h>
#include <openbabel/obfunction.h>
#include <float.h>

#ifndef BUFF_SIZE
#define BUFF_SIZE 32768
#endif


namespace OpenBabel
{
  namespace LineSearchType 
  {
    enum {
      Simple, Newton2Num 
    };
  };
  
  // Class OBForceField
  // class introduction in minimize.cpp
  class OBMinimizePrivate;
  class OBFPRT OBMinimize
  {
  protected:
    OBMinimizePrivate * const d; //!< the d-pointer

  public:
    explicit OBMinimize(OBForceField *forcefield);
    /** 
     * @brief Destructor.
     */
    virtual ~OBMinimize();
    
    /////////////////////////////////////////////////////////////////////////
    // Energy Minimization                                                 //
    /////////////////////////////////////////////////////////////////////////
      
    //! \name Methods for energy minimization
    //@{
    /** 
     * @brief Set the LineSearchType. The default type is LineSearchType::Simple.
     * 
     * @param type The LineSearchType to be used in SteepestDescent and ConjugateGradients.
     */ 
    void SetLineSearchType(int type);
    /** 
     * @brief Get the LineSearchType.
     * 
     * @return The current LineSearchType.
     */ 
    int GetLineSearchType();
    /** 
     * @brief Perform a linesearch for the entire molecule in direction p direction. 
     * This function is called when using LineSearchType::Simple.
     * 
     * @param currentCoords Start coordinates.
     * @param direction The search direction.
     * 
     * @return alpha, The scale of the step we moved along the direction vector.
     *
     * @par Output to log:
     *  OBFF_LOGLVL_NONE:   none \n
     *  OBFF_LOGLVL_LOW:    none \n
     *  OBFF_LOGLVL_MEDIUM: none \n
     *  OBFF_LOGLVL_HIGH:   none \n
     */
    double LineSearch(std::vector<Eigen::Vector3d> &currentCoords, std::vector<Eigen::Vector3d> &direction);
    /** 
     * @brief Perform a linesearch for the entire molecule.
     * This function is called when using LineSearchType::Newton2Num.
     *
     * @param direction The search direction.
     * 
     * @return alpha, The scale of the step we moved along the direction vector.
     *
     * @par Output to log:
     *  OBFF_LOGLVL_NONE:   none \n
     *  OBFF_LOGLVL_LOW:    none \n
     *  OBFF_LOGLVL_MEDIUM: none \n
     *  OBFF_LOGLVL_HIGH:   none \n
     */
    double Newton2NumLineSearch(std::vector<Eigen::Vector3d> &direction);
    /** 
     * @brief Set the coordinates of the atoms to origCoord + step.
     * 
     * @param origCoords Start coordinates.
     * @param direction The search direction.
     * @param step The step to take.
     */
    void   LineSearchTakeStep(std::vector<Eigen::Vector3d> &origCoords, 
        std::vector<Eigen::Vector3d> &direction, double step);
    /** 
     * @brief Perform steepest descent optimalization for steps steps or until convergence criteria is reached.
     * 
     * @param steps The number of steps.
     * @param econv Energy convergence criteria. (defualt is 1e-6)
     *
     * @par Output to log:
     *  This function should only be called with the log level set to OBFF_LOGLVL_NONE or OBFF_LOGLVL_LOW. Otherwise
     *  too much information about the energy calculations needed for the minimization will interfere with the list 
     *  of energies for succesive steps. \n\n
     *  OBFF_LOGLVL_NONE:   none \n
     *  OBFF_LOGLVL_LOW:    header including number of steps and first step \n
     *  OBFF_LOGLVL_MEDIUM: see note above \n
     *  OBFF_LOGLVL_HIGH:   see note above \n 
     */
    void SteepestDescent(int steps, double econv = 1e-6f);
    /**
     * @brief Initialize steepest descent optimalization, to be used in combination with SteepestDescentTakeNSteps().
     * 
     * example:
     * @code
     * // pFF is a pointer to a OBForceField class 
     * pFF->SteepestDescentInitialize(100, 1e-5f);
     * while (pFF->SteepestDescentTakeNSteps(5)) {
     *   // do some updating in your program (redraw structure, ...)
     * }
     * @endcode
     *
     * If you don't need any updating in your program, SteepestDescent() is recommended.
     *
     * @param steps The number of steps.
     * @param econv Energy convergence criteria. (defualt is 1e-6)
     *
     * @par Output to log:
     *  This function should only be called with the log level set to OBFF_LOGLVL_NONE or OBFF_LOGLVL_LOW. Otherwise
     *  too much information about the energy calculations needed for the minimization will interfere with the list 
     *  of energies for succesive steps. \n\n
     *  OBFF_LOGLVL_NONE:   none \n
     *  OBFF_LOGLVL_LOW:    header including number of steps \n
     *  OBFF_LOGLVL_MEDIUM: see note above \n
     *  OBFF_LOGLVL_HIGH:   see note above \n
     */
    void SteepestDescentInitialize(int steps = 1000, double econv = 1e-6f);
    /** 
     * @brief Take n steps in a steepestdescent optimalization that was 
     * previously initialized with SteepestDescentInitialize().
     *  
     * @param n The number of steps to take.
     * 
     * @return False if convergence or the number of steps given by SteepestDescentInitialize() has been reached.
     * 
     * @par Output to log:
     *  This function should only be called with the log level set to OBFF_LOGLVL_NONE or OBFF_LOGLVL_LOW. Otherwise
     *  too much information about the energy calculations needed for the minimization will interfere with the list 
     *  of energies for succesive steps. \n\n
     *  OBFF_LOGLVL_NONE:   none \n
     *  OBFF_LOGLVL_LOW:    step number, energy and energy for the previous step \n
     *  OBFF_LOGLVL_MEDIUM: see note above \n
     *  OBFF_LOGLVL_HIGH:   see note above \n
     */
    bool SteepestDescentTakeNSteps(int n);
    /** 
     * @brief Perform conjugate gradient optimalization for steps steps or until convergence criteria is reached.
     * 
     * @param steps The number of steps. 
     * @param econv Energy convergence criteria. (defualt is 1e-6)
     *
     * @par Output to log:
     *  This function should only be called with the log level set to OBFF_LOGLVL_NONE or OBFF_LOGLVL_LOW. Otherwise
     *  too much information about the energy calculations needed for the minimization will interfere with the list 
     *  of energies for succesive steps. \n\n
     *  OBFF_LOGLVL_NONE:   none \n
     *  OBFF_LOGLVL_LOW:    information about the progress of the minimization \n
     *  OBFF_LOGLVL_MEDIUM: see note above \n
     *  OBFF_LOGLVL_HIGH:   see note above \n
     */
    void ConjugateGradients(int steps, double econv = 1e-6f);
    /** 
     * @brief Initialize conjugate gradient optimalization and take the first step, to be 
     * used in combination with ConjugateGradientsTakeNSteps().
     * 
     * example:
     * @code
     * // pFF is a pointer to a OBForceField class 
     * pFF->ConjugateGradientsInitialize(100, 1e-5f);
     * while (pFF->ConjugateGradientsTakeNSteps(5)) {
     *   // do some updating in your program (redraw structure, ...)
     * }
     * @endcode
     * 
     * If you don't need any updating in your program, ConjugateGradients() is recommended.
     *
     * @param steps The number of steps.
     * @param econv Energy convergence criteria. (defualt is 1e-6)
     * 
     * @par Output to log:
     *  This function should only be called with the log level set to OBFF_LOGLVL_NONE or OBFF_LOGLVL_LOW. Otherwise
     *  too much information about the energy calculations needed for the minimization will interfere with the list 
     *  of energies for succesive steps. \n\n
     *  OBFF_LOGLVL_NONE:   none \n
     *  OBFF_LOGLVL_LOW:    header including number of steps and first step \n
     *  OBFF_LOGLVL_MEDIUM: see note above \n
     *  OBFF_LOGLVL_HIGH:   see note above \n
     */
    void ConjugateGradientsInitialize(int steps = 1000, double econv = 1e-6f);
    /** 
     * @brief Take n steps in a conjugate gradient optimalization that was previously 
     * initialized with ConjugateGradientsInitialize().
     *
     * @param n The number of steps to take.
     * @return False if convergence or the number of steps given by ConjugateGradientsInitialize() has been reached.
     *  
     * @par Output to log:
     *  This function should only be called with the log level set to OBFF_LOGLVL_NONE or OBFF_LOGLVL_LOW. Otherwise
     *  too much information about the energy calculations needed for the minimization will interfere with the list 
     *  of energies for succesive steps. \n\n
     *  OBFF_LOGLVL_NONE:   none \n
     *  OBFF_LOGLVL_LOW:    step number, energy and energy for the previous step \n
     *  OBFF_LOGLVL_MEDIUM: see note above \n
     *  OBFF_LOGLVL_HIGH:   see note above \n
    */
    bool ConjugateGradientsTakeNSteps(int n);
    //@}
    
  }; // class OBForceField

}// namespace OpenBabel

#endif   // OB_MINIMIZE_H

//! @file minimize.h
//! @brief Handle minimization

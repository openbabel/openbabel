/**********************************************************************
forcefield.h - Handle OBForceField class.
 
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

#ifndef OB_FORCEFIELD_H
#define OB_FORCEFIELD_H

#include <vector>
#include <string>
#include <map>

#include <list>
#include <set>
#include <openbabel/base.h>
#include <openbabel/mol.h>
#include <openbabel/pluginiter.h>

namespace OpenBabel
{
  // log levels
#define OBFF_LOGLVL_NONE	0   //!< no output
#define OBFF_LOGLVL_LOW	1     //!< SteepestDescent progress... (no output from Energy())
#define OBFF_LOGLVL_MEDIUM	2 //!< individual energy terms
#define OBFF_LOGLVL_HIGH	3   //!< individual calculations and parameters

  // terms
#define OBFF_ENERGY		(1 << 0)   //!< all terms
#define OBFF_EBOND		(1 << 1)   //!< bond term
#define OBFF_EANGLE		(1 << 2)   //!< angle term
#define OBFF_ESTRBND		(1 << 3) //!< strbnd term
#define OBFF_ETORSION		(1 << 4) //!< torsion term
#define OBFF_EOOP		(1 << 5)     //!< oop term
#define OBFF_EVDW		(1 << 6)     //!< vdw term
#define OBFF_EELECTROSTATIC	(1 << 7) //!< electrostatic term

  // mode arguments for SteepestDescent, ConjugateGradients, ...
#define OBFF_NUMERICAL_GRADIENT   (1 << 0)  //!< use numerical gradients
#define OBFF_ANALYTICAL_GRADIENT	(1 << 1)  //!< use analytical gradients

#define KCAL_TO_KJ	4.1868

  // inline if statements for logging.
#define IF_OBFF_LOGLVL_LOW    if(loglvl >= OBFF_LOGLVL_LOW)
#define IF_OBFF_LOGLVL_MEDIUM if(loglvl >= OBFF_LOGLVL_MEDIUM)
#define IF_OBFF_LOGLVL_HIGH   if(loglvl >= OBFF_LOGLVL_HIGH)
  
  //! \class OBFFParameter forcefield.h <openbabel/forcefield.h>
  //! \brief Internal class for OBForceField to hold forcefield parameters
  class OBFFParameter {
  public:
    //! Used to store integer atom types
    int         a, b, c, d;
    //! used to store string atom types
    std::string _a, _b, _c, _d; 

    //! Used to store integer type parameters (bondtypes, multiplicity, ...)
    int       ipar1, ipar2, ipar3, ipar4, ipar5;
    //! Used to store double type parameters (force constants, bond lengths, angles, ...)
    double    dpar1, dpar2, dpar3, dpar4, dpar5;

    //! Assignment 
    OBFFParameter& operator=(const OBFFParameter &ai) 
    {
      if (this != &ai) {
        a = ai.a;
        b = ai.b;
        c = ai.c;
        d = ai.d;
        _a = ai._a;
        _b = ai._b;
        _c = ai._c;
        _d = ai._d;
        ipar1 = ai.ipar1;
        ipar2 = ai.ipar2;
        ipar3 = ai.ipar3;
        ipar4 = ai.ipar4;
        ipar5 = ai.ipar5;
        dpar1 = ai.dpar1;
        dpar2 = ai.dpar2;
        dpar3 = ai.dpar3;
        dpar4 = ai.dpar4;
        dpar5 = ai.dpar5;
      }
        
      return *this;
    }

    //! Reset the atom types and set all parameters to zero
    void clear () 
    {
      a = 0;
      b = 0;
      c = 0;
      d = 0;
      ipar1 = 0;
      ipar2 = 0;
      ipar3 = 0;
      ipar4 = 0;
      ipar5 = 0;
      dpar1 = 0.0f;
      dpar2 = 0.0f;
      dpar3 = 0.0f;
      dpar4 = 0.0f;
      dpar5 = 0.0f;
    }
  }; // class OBFFParameter
  
  // specific class introductions in forcefieldYYYY.cpp (for YYYY calculations)
  //! \class OBFFCalculation forcefield.h <openbabel/forcefield.h>
  //! \brief Internal class for OBForceField to hold energy and gradient calculations on specific force fields
  class OBFFCalculation
  {
    public:
      //! Used to store the energy for this OBFFCalculation
      double energy;
      //! Used to store the gradients for this OBFFCalculation
      vector3 grada, gradb, gradc, gradd;
      //! Used to store the atoms for this OBFFCalculation
      OBAtom *a, *b, *c, *d;

      //! Constructor
      OBFFCalculation() 
        {
	  a = NULL;
	  b = NULL;
	  c = NULL;
	  d = NULL;
	  energy = 0.0f;
          grada = VZero;
          gradb = VZero;
          gradc = VZero;
          gradd = VZero;
        }
      //! Destructor
      virtual ~OBFFCalculation()
        {
        }
      
      //! Compute the energy and gradients for this OBFFCalculation
      virtual void Compute(bool gradients = true) 
        {
        }
      //! \return Energy for this OBFFCalculation (call Compute() first)
      virtual double GetEnergy() 
      {
        if (!energy)
	  Compute(false);

        return energy; 
      }
      //! \return Gradient for this OBFFCalculation with respect to coordinates of atom (call Compute() first)
      virtual vector3 GetGradient(OBAtom *atom) 
      {
        if (atom == a)
          return grada;
        else if (atom == b)
          return gradb;
        else if (atom == c)
          return gradc;
        else if (atom == d)
          return gradd;
        else 
          return  VZero;
      }
  };

  // Class OBForceField
  // class introduction in forcefield.cpp
  class OBAPI OBForceField
  {
  
  MAKE_PLUGIN(OBForceField)
  
    protected:
    /*! 
      Get the correct OBFFParameter from a OBFFParameter vector.
       
      \code vector<OBFFParameter> parameters; \endcode
      
      this vector is filled with entries (as OBFFParameter) from 
      a parameter file. This happens in the Setup() function.
      
      \code GetParameter(a, 0, 0, 0, parameters); \endcode
        
      returns the first OBFFParameter from vector<OBFFParameter> 
      parameters where: pa = a (pa = parameter.a)
      
      use: vdw parameters, ...
      
      \code GetParameter(a, b, 0, 0, parameters); \endcode
      
      returns the first OBFFParameter from vector<OBFFParameter> 
      parameters where: pa = a & pb = b      (ab)
      or: pa = b & pb = a      (ba)
            
      use: bond parameters, vdw parameters (pairs), ...
      
      \code GetParameter(a, b, c, 0, parameters); \endcode
      
      returns the first OBFFParameter from vector<OBFFParameter> 
      parameters where: pa = a & pb = b & pc = c     (abc)
      or: pa = c & pb = b & pc = a     (cba)
      
      use: angle parameters, ...
      
      \code GetParameter(a, b, c, d, parameters); \endcode
      
      returns the first OBFFParameter from vector<OBFFParameter> 
      parameters where: pa = a & pb = b & pc = c & pd = d    (abcd)
      or: pa = d & pb = b & pc = c & pd = a    (dbca)
      or: pa = a & pb = c & pc = b & pd = d    (acbd)
      or: pa = d & pb = c & pc = b & pd = a    (dcba)
      
      use: torsion parameters, ...
    */
    OBFFParameter* GetParameter(int a, int b, int c, int d, std::vector<OBFFParameter> &parameter);
    //! see GetParameter(int a, int b, int c, int d, std::vector<OBFFParameter> &parameter)
    OBFFParameter* GetParameter(const char* a, const char* b, const char* c, const char* d, std::vector<OBFFParameter> &parameter);
           
    //! Calculate the potential energy function derivative numerically with repect to the coordinates of atom with index a (this vector is the gradient)
    /*!
     * \param a  provides coordinates
     * \param terms OBFF_ENERGY, OBFF_EBOND, OBFF_EANGLE, OBFF_ESTRBND, OBFF_ETORSION, OBFF_EOOP, OBFF_EVDW, OBFF_ELECTROSTATIC
     * \return the negative gradient of atom a
     */
    vector3 NumericalDerivative(OBAtom *a, int terms = OBFF_ENERGY);
    //! OB 3.0
    vector3 NumericalSecondDerivative(OBAtom *a, int terms = OBFF_ENERGY);
    /*! Calculate the potential energy function derivative analyticaly with repect to the coordinates of atom with index a (this vector is the gradient)
     *
     *  If the currently selected forcefield doesn't have analytical gradients, 
     *  we can still call this function which will return the result of 
     *  NumericalDerivative()
     *  \param a  provides coordinates
     *  \param terms OBFF_ENERGY, OBFF_EBOND, OBFF_EANGLE, OBFF_ESTRBND, OBFF_ETORSION, OBFF_EOOP, OBFF_EVDW, OBFF_ELECTROSTATIC
     *  \return the negative gradient of atom a
     */
    virtual vector3 GetGradient(OBAtom *a, int terms = OBFF_ENERGY) 
    { 
      return -NumericalDerivative(a, terms); 
    }
    /*! Check if two atoms are in the same ring
     *  \param a atom a
     *  \param b atom b
     *  \return true if atom a and b are in the same ring
     */
    bool IsInSameRing(OBAtom* a, OBAtom* b);
 
    OBMol _mol; //!< Molecule to be evaluated or minimized

    //! Output for logfile
    std::ostream* logos;
    char logbuf[200];
    int loglvl; //!< Log level for output

    //! used to hold i for current conformer (needed by UpdateConformers)
    int current_conformer;

    //! Used for conjugate gradients and steepest descent(Initialize and TakeNSteps)
    double _econv, _e_n1;
    int _method, _cstep, _nsteps;
    std::vector<vector3> _grad1, _dir1;

  public:
    //! Destructor
    virtual ~OBForceField()
      {
      }
    //! short description of the force field type.
    //virtual std::string Description()=0;
    /*! \param ID forcefield id (Ghemical, ...)
     *  \return A pointer to a forcefield (the default if ID is empty), or NULL if not available
     */
    static OBForceField* FindForceField(const std::string& ID)
    { 
      return Iter().FindType(ID);
    } 
    /*! \param ID forcefield id (Ghemical, ...)
     *  \return A pointer to a forcefield (the default if ID is empty), or NULL if not available
     */
    static OBForceField* FindForceField(const char *ID)
    {
      std::string ffname(ID);
      return FindForceField(ffname);
    }
    //! \return The unit (kcal/mol, kJ/mol, ...) in which the energy is expressed as std::string
    virtual std::string GetUnit() { return std::string("au"); }
    /*! Setup the forcefield for mol (assigns atom types, charges, etc.) 
     *  \param mol the OBMol object that contains the atoms and bonds
     *  \return True if succesfull
     */
    virtual bool Setup(OBMol &mol) { return false; }
    /*! Update coordinates for current conformer
     *  \param mol the OBMol object to copy the coordinates to
     *  \return true if succesfull
     */
    bool UpdateCoordinates(OBMol &mol);
    /*! Update coordinates for all conformers
     *  \param mol the OBMol object to copy the coordinates to
     *  \return true if succesfull
     */
    bool UpdateConformers(OBMol &mol);
    /*! Print msg to the logfile
     *  \param msg the message
     */
    void OBFFLog(std::string msg)
    {
      if (!logos)
        return;
      
      *logos << msg;
    }
    /*! Print msg to the logfile
     *  \param msg the message
     */
    void OBFFLog(const char *msg)
    {
      if (!logos)
        return;
      
      *logos << msg;
    }

    /////////////////////////////////////////////////////////////////////////
    // Energy Evaluation                                                   //
    /////////////////////////////////////////////////////////////////////////
      
    //! \name Methods for energy evaluation
    //@{
    /*! \param gradients Set to true when the gradients need to be calculated (needs to be done before calling GetGradient())
     *  \return Total energy
     *   \par Output to log:
     *    OBFF_LOGLVL_NONE:   none \n
     *    OBFF_LOGLVL_LOW:    none \n
     *    OBFF_LOGLVL_MEDIUM: energy for indivudual energy terms \n
     *    OBFF_LOGLVL_HIGH:   energy for individual energy interactions \n
     */
    virtual double Energy(bool gradients = true) { return 0.0f; }
    /*! \param gradients Set to true when the gradients need to be calculated (needs to be done before calling GetGradient())
     *  \return Bond stretching energy
     *   \par Output to log:
     *    see Energy()
     */
    virtual double E_Bond(bool gradients = true) { return 0.0f; }
    /*! \param gradients Set to true when the gradients need to be calculated (needs to be done before calling GetGradient())
     *  \return Angle bending energy
     *  \par Output to log:
     *   see Energy()
     */
    virtual double E_Angle(bool gradients = true) { return 0.0f; }
    /*! \param gradients Set to true when the gradients need to be calculated (needs to be done before calling GetGradient())
     *  \return Stretch bending energy
     *   \par Output to log:
     *    see Energy()
     */ 
    virtual double E_StrBnd(bool gradients = true) { return 0.0f; }
    /*! \param gradients Set to true when the gradients need to be calculated (needs to be done before calling GetGradient())
     *  \return Torsional energy
     *    \par Output to log:
     *	  see Energy()
     */ 
    virtual double E_Torsion(bool gradients = true) { return 0.0f; }
    /*! \param gradients Set to true when the gradients need to be calculated (needs to be done before calling GetGradient())
     *  \return Out-Of-Plane bending energy
     *   \par Output to log:
     *	  see Energy()
     */ 
    virtual double E_OOP(bool gradients = true) { return 0.0f; }
    /*! \param gradients Set to true when the gradients need to be calculated (needs to be done before calling GetGradient())
     *  \return Van der Waals energy
     *   \par Output to log:
     *	  see Energy()
     */ 
    virtual double E_VDW(bool gradients = true) { return 0.0f; }
    /*! \param gradients Set to true when the gradients need to be calculated (needs to be done before calling GetGradient())
     *  \return Electrostatic energy
     *   \par Output to log:
     *	  see Energy()
     */ 
    virtual double E_Electrostatic(bool gradients = true) { return 0.0f; }
    //@}
     
    /////////////////////////////////////////////////////////////////////////
    // Logging                                                             //
    /////////////////////////////////////////////////////////////////////////
      
    //! \name Methods for logging
    //@{
    /*! Set the stream for logging (can also be &cout for logging to screen)
     *  \param pos stream
     *  \return True if succesfull
     */
    bool SetLogFile(std::ostream *pos);
    /*!
      Set the log level (OBFF_LOGLVL_NONE, OBFF_LOGLVL_LOW, OBFF_LOGLVL_MEDIUM, OBFF_LOGLVL_HIGH)

      Inline if statements for logging are available: 
	
      \code
      #define IF_OBFF_LOGLVL_LOW    if(loglvl >= OBFF_LOGLVL_LOW)
      #define IF_OBFF_LOGLVL_MEDIUM if(loglvl >= OBFF_LOGLVL_MEDIUM)
      #define IF_OBFF_LOGLVL_HIGH   if(loglvl >= OBFF_LOGLVL_HIGH)
      \endcode

      example:
      \code
      SetLogLevel(OBFF_LOGLVL_MEDIUM);
  
      IF_OBFF_LOGLVL_HIGH {
      *logos << "this text will NOT be logged..." << endl
      }
   
      IF_OBFF_LOGLVL_LOW {
      *logos << "this text will be logged..." << endl
      }
  
      IF_OBFF_LOGLVL_MEDIUM {
      *logos << "this text will also be logged..." << endl
      }
      \endcode
    */
    bool SetLogLevel(int level);
    //! \return log level
    int GetLogLevel() { return loglvl; }
    //@}
     
    /////////////////////////////////////////////////////////////////////////
    // Structure Generation                                                //
    /////////////////////////////////////////////////////////////////////////
      
    //! \name Methods for structure generation
    //@{
    /*! Generate conformers for the molecule (systematicaly rotating torsions).
     *  
     *  The initial starting structure here is important, this structure should be
     *  minimized for the best results. SystematicRotorSearch works by rotating around
     *  the rotatable bond in a molecule (see OBRotamerList class). This rotating generates 
     *  multiple conformers. The energy for all these conformers is then evaluated and the 
     *  lowest energy conformer is selected.
     *        
     *	\par Output to log:
     *  This function should only be called with the log level set to OBFF_LOGLVL_NONE or OBFF_LOGLVL_LOW. Otherwise
     *	too much information about the energy calculations needed for this function will interfere with the output for
     *	this function. \n\n
     *  OBFF_LOGLVL_NONE:   none \n
     *  OBFF_LOGLVL_LOW:    number of rotatable bonds, energies for the conformers, which one is the lowest, ... \n
     *  OBFF_LOGLVL_MEDIUM: see note above \n
     *  OBFF_LOGLVL_HIGH:   see note above \n 
     */
    void SystematicRotorSearch();
      
    /////////////////////////////////////////////////////////////////////////
    // Energy Minimization                                                 //
    /////////////////////////////////////////////////////////////////////////
      
    //! \name Methods for energy minimization
    //@{
    /*! Perform a linesearch starting at atom in direction direction
      \param atom start coordinates
      \param direction the search direction

      \return vector which starts at atom and stops at the minimum (the length is the ideal stepsize)

      \par Output to log:
        OBFF_LOGLVL_NONE:   none \n
        OBFF_LOGLVL_LOW:    none \n
        OBFF_LOGLVL_MEDIUM: none \n
        OBFF_LOGLVL_HIGH:   none \n
    */
    vector3 LineSearch(OBAtom *atom, vector3 &direction);
    /*! Perform steepest descent optimalization for steps steps or until convergence criteria is reached.
      \param steps the number of steps 
      \param econv energy convergence criteria (defualt is 1e-6)
      \param method OBFF_ANALYTICAL_GRADIENTS (default) or OBFF_NUMERICAL_GRADIENTS

      \par Output to log:
        This function should only be called with the log level set to OBFF_LOGLVL_NONE or OBFF_LOGLVL_LOW. Otherwise
	too much information about the energy calculations needed for the minimization will interfere with the list 
	of energies for succesive steps. \n\n
        OBFF_LOGLVL_NONE:   none \n
        OBFF_LOGLVL_LOW:    header including number of steps and first step \n
        OBFF_LOGLVL_MEDIUM: see note above \n
        OBFF_LOGLVL_HIGH:   see note above \n 
    */
    void SteepestDescent(int steps, double econv = 1e-6f, int method = OBFF_ANALYTICAL_GRADIENT);
    /*! Initialize steepest descent optimalization, to be used in combination with SteepestDescentTakeNSteps().
      
      example:
      \code
      // pFF is a pointer to a OBForceField class 
      pFF->SteepestDescentInitialize(100, 1e-5f);
      while (pFF->SteepestDescentTakeNSteps(5)) {
        // do some updating in your program (redraw structure, ...)
      }
      \endcode
      
      If you don't need any updating in your program, SteepestDescent() is recommended.

      \param steps the number of steps 
      \param econv energy convergence criteria (defualt is 1e-6)
      \param method OBFF_ANALYTICAL_GRADIENTS (default) or OBFF_NUMERICAL_GRADIENTS

      \par Output to log:
        This function should only be called with the log level set to OBFF_LOGLVL_NONE or OBFF_LOGLVL_LOW. Otherwise
	too much information about the energy calculations needed for the minimization will interfere with the list 
	of energies for succesive steps. \n\n
	OBFF_LOGLVL_NONE:   none \n
        OBFF_LOGLVL_LOW:    header including number of steps \n
        OBFF_LOGLVL_MEDIUM: see note above \n
        OBFF_LOGLVL_HIGH:   see note above \n
 
    */
    void SteepestDescentInitialize(int steps = 1000, double econv = 1e-6f, int method = OBFF_ANALYTICAL_GRADIENT);
    /*! Take n steps in a steepestdescent optimalization that was previously initialized with SteepestDescentInitialize().
      \param n the number of steps to take
      
      \return false if convergence or the number of steps given by SteepestDescentInitialize() has been reached 
      
      \par Output to log:
        This function should only be called with the log level set to OBFF_LOGLVL_NONE or OBFF_LOGLVL_LOW. Otherwise
	too much information about the energy calculations needed for the minimization will interfere with the list 
	of energies for succesive steps. \n\n
	OBFF_LOGLVL_NONE:   none \n
        OBFF_LOGLVL_LOW:    step number, energy and energy for the previous step \n
        OBFF_LOGLVL_MEDIUM: see note above \n
        OBFF_LOGLVL_HIGH:   see note above \n
 
    */
    bool SteepestDescentTakeNSteps(int n);
    /*! Perform conjugate gradient optimalization for steps steps or until convergence criteria is reached.
      \param steps the number of steps 
      \param econv energy convergence criteria (defualt is 1e-6)
      \param method OBFF_ANALYTICAL_GRADIENTS (default) or OBFF_NUMERICAL_GRADIENTS

      \par Output to log:
        This function should only be called with the log level set to OBFF_LOGLVL_NONE or OBFF_LOGLVL_LOW. Otherwise
	too much information about the energy calculations needed for the minimization will interfere with the list 
	of energies for succesive steps. \n\n
	OBFF_LOGLVL_NONE:   none \n
        OBFF_LOGLVL_LOW:    information about the progress of the minimization \n
        OBFF_LOGLVL_MEDIUM: see note above \n
        OBFF_LOGLVL_HIGH:   see note above \n
 
    */
    void ConjugateGradients(int steps, double econv = 1e-6f, int method = OBFF_ANALYTICAL_GRADIENT);
    /*! Initialize conjugate gradient optimalization and take the first step, to be used in combination with ConjugateGradientsTakeNSteps().
      
      example:
      \code
      // pFF is a pointer to a OBForceField class 
      pFF->ConjugateGradientsInitialize(100, 1e-5f);
      while (pFF->ConjugateGradientsTakeNSteps(5)) {
        // do some updating in your program (redraw structure, ...)
      }
      \endcode
      
      If you don't need any updating in your program, ConjugateGradients() is recommended.

      \param steps the number of steps 
      \param econv energy convergence criteria (defualt is 1e-6)
      \param method OBFF_ANALYTICAL_GRADIENTS (default) or OBFF_NUMERICAL_GRADIENTS
      
      \par Output to log:
        This function should only be called with the log level set to OBFF_LOGLVL_NONE or OBFF_LOGLVL_LOW. Otherwise
	too much information about the energy calculations needed for the minimization will interfere with the list 
	of energies for succesive steps. \n\n
	OBFF_LOGLVL_NONE:   none \n
        OBFF_LOGLVL_LOW:    header including number of steps and first step \n
        OBFF_LOGLVL_MEDIUM: see note above \n
        OBFF_LOGLVL_HIGH:   see note above \n
    */
    void ConjugateGradientsInitialize(int steps = 1000, double econv = 1e-6f, int method = OBFF_ANALYTICAL_GRADIENT);
    /*! Take n steps in a conjugate gradient optimalization that was previously initialized with ConjugateGradientsInitialize().
      \param n the number of steps to take
      
      \return false if convergence or the number of steps given by ConjugateGradientsInitialize() has been reached 
      
      \par Output to log:
        This function should only be called with the log level set to OBFF_LOGLVL_NONE or OBFF_LOGLVL_LOW. Otherwise
	too much information about the energy calculations needed for the minimization will interfere with the list 
	of energies for succesive steps. \n\n
	OBFF_LOGLVL_NONE:   none \n
        OBFF_LOGLVL_LOW:    step number, energy and energy for the previous step \n
        OBFF_LOGLVL_MEDIUM: see note above \n
        OBFF_LOGLVL_HIGH:   see note above \n
    */
    bool ConjugateGradientsTakeNSteps(int n);
    //@}
      
    /////////////////////////////////////////////////////////////////////////
    // Validation                                                          //
    /////////////////////////////////////////////////////////////////////////
      
    //! \name Methods for forcefield validation
    //@{
    //! Validate the force field implementation (debugging)
    virtual bool Validate() { return false; }
    /*! 
      Validate the analytical gradients by comparing them to numerical ones. This function has to
      be implemented force field specific. (debugging)
    */
    virtual bool ValidateGradients() { return false; }
    /*! 
      Calculate the error of the analytical gradient (debugging)
      \return  error = fabs(numgrad - anagrad) / anagrad * 100% 
   */
    vector3 ValidateGradientError(vector3 &numgrad, vector3 &anagrad);
    //@}
     
    /////////////////////////////////////////////////////////////////////////
    // Vector Analysis                                                     //
    /////////////////////////////////////////////////////////////////////////
      
    //! \name Methods for vector analysis (used by OBFFXXXXCalculationYYYY)
    //@{
    /*! Calculate the derivative of a vector length. The vector is given by a - b, 
     * the length of this vector rab = sqrt(ab.x^2 + ab.y^2 + ab.z^2).
     * \param a atom a (coordinates), will be changed to -drab/da
     * \param b atom b (coordinates), will be changed to -drab/db
     * \return The distance between a and b (bondlength for bond stretching, separation for vdw, electrostatic)
     */
    static double VectorLengthDerivative(vector3 &a, vector3 &b);
    /*! Calculate the derivative of a angle a-b-c. The angle is given by dot(ab,cb)/rab*rcb.
     * \param a atom a (coordinates), will be changed to -dtheta/da
     * \param b atom b (coordinates), will be changed to -dtheta/db
     * \param c atom c (coordinates), will be changed to -dtheta/dc
     * \return The angle between a-b-c
     */
    static double VectorAngleDerivative(vector3 &a, vector3 &b, vector3 &c);
    /*! Calculate the derivative of a torsion angle a-b-c-d. The torsion angle is given by dot(corss(ab,bc),cross(bc,cd)/rabbc*rbccd.
     * \param a atom a (coordinates), will be changed to -dtheta/da
     * \param b atom b (coordinates), will be changed to -dtheta/db
     * \param c atom c (coordinates), will be changed to -dtheta/dc
     * \param d atom d (coordinates), will be changed to -dtheta/dd
     * \return The tosion angle for a-b-c-d
     */
    static double VectorTorsionDerivative(vector3 &a, vector3 &b, vector3 &c, vector3 &d);
    //@}

    ///Do not use.
    ///This function contains rubbish merely to ensure the compiler instantiates
    ///some templated functions which are needed for the Windows Python build.
    ///TODO Find the proper way of doing this.
    void kludge()
    {
      PluginIter<OBForceField> it;
      (++it)->SetLogLevel(1);
      if(it)it.ID();
    }

  }; // class OBForceField

}// namespace OpenBabel

#endif   // OB_FORCEFIELD_H

//! \file forcefield.h
//! \brief Handle forcefields

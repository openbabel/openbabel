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
#include <openbabel/plugin.h>

#ifndef OBFPRT
#define OBFPRT
#endif

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

  // constraint types
#define OBFF_CONST_IGNORE	(1 << 0)   //!< ignore the atom while setting up calculations
#define OBFF_CONST_ATOM		(1 << 1)   //!< fix the atom position
#define OBFF_CONST_ATOM_X	(1 << 2)   //!< fix the x coordinate of the atom position
#define OBFF_CONST_ATOM_Y	(1 << 3)   //!< fix the y coordinate of the atom position
#define OBFF_CONST_ATOM_Z	(1 << 4)   //!< fix the z coordinate of the atom position
#define OBFF_CONST_BOND		(1 << 5)   //!< constrain bond length
#define OBFF_CONST_ANGLE	(1 << 6)   //!< constrain angle
#define OBFF_CONST_TORSION	(1 << 7)   //!< constrain torsion

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
  class OBFPRT OBFFParameter {
  public:
    //! Used to store integer atom types
    int         a, b, c, d;
    //! used to store string atom types
    std::string _a, _b, _c, _d; 

    //! Used to store integer type parameters (bondtypes, multiplicity, ...)
    std::vector<int>    _ipar;
    //! Used to store double type parameters (force constants, bond lengths, angles, ...)
    std::vector<double> _dpar;

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
        _ipar = ai._ipar;
        _dpar = ai._dpar;
      }
        
      return *this;
    }

    //! Reset the atom types and set all parameters to zero
    void clear () 
    {
      a = b = c = d = 0;
      _ipar.clear();
      _dpar.clear();
    }
  }; // class OBFFParameter
  
  // specific class introductions in forcefieldYYYY.cpp (for YYYY calculations)
  //! \class OBFFCalculation forcefield.h <openbabel/forcefield.h>
  //! \brief Internal class for OBForceField to hold energy and gradient calculations on specific force fields
  class OBFPRT OBFFCalculation
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
  
 //! \class OBFFConstraint forcefield.h <openbabel/forcefield.h>
  //! \brief Internal class for OBForceField to hold constraints
  class OBFPRT OBFFConstraint
  {
    public:
      //! Used to store the contraint energy for this OBFFConstraint
      double constraint_energy, constraint_value;
      //! Used to store the contraint type for this OBFFConstraint
      int type;
      //! Used to store the atoms for this OBFFCostraint
      OBAtom *a, *b, *c, *d;
      //! Used to store the gradients for this OBFFCalculation
      vector3 grada, gradb, gradc, gradd;

      //! Constructor
      OBFFConstraint() 
      {
	a = NULL;
	b = NULL;
	c = NULL;
        d = NULL;
        constraint_energy = 0.0f;
      }
      //! Destructor
      ~OBFFConstraint()
      {
      }
      
      //! Compute the constraint energy for this OBFFConstraint
      void Compute();
      //! \return Constraint energy for this OBFFConstraint (call Compute() first)
      double GetConstraintEnergy() 
      {
        if (!constraint_energy)
	  Compute();

        return constraint_energy; 
      }
  };

  //! \class OBFFConstraints forcefield.h <openbabel/forcefield.h>
  //! \brief Internal class for OBForceField to handle constraints
  class OBFPRT OBFFConstraints
  {
    public:
      //! Constructor
      OBFFConstraints()
      {
      }
      //! Destructor
      ~OBFFConstraints()
      {
        _constraints.clear();
      }
      //! Clear all constraints
      void Clear();
      //! Get the constrain energy
      double GetConstraintEnergy();
      //! Get the constrain gradient for the atom
      //GetConstraintGradient(); isn't need I think??? only testing will tell :s
      OBFFConstraints& operator=(const OBFFConstraints &ai) 
      {
        if (this != &ai) {
          _constraints = ai._constraints;
        }
        return *this;
      }

      /////////////////////////////////////////////////////////////////////////
      // Set Constraints                                                     //
      /////////////////////////////////////////////////////////////////////////
      //! \name Methods to set constraints
      //@{
      //! Fix the position of an atom
      void AddAtomConstraint(OBAtom *atom);
      //! Fix the x coordinate of the atom position
      void AddAtomXConstraint(OBAtom *atom);
      //! Fix the y coordinate of the atom position
      void AddAtomYConstraint(OBAtom *atom);
      //! Fix the z coordinate of the atom position
      void AddAtomZConstraint(OBAtom *atom);
      //! Constrain the bond length a-b
      void AddBondConstraint(OBAtom *a, OBAtom *b, double length);
      //! Constrain the angle a-b-c
      void AddAngleConstraint(OBAtom *a, OBAtom *b, OBAtom *c, double angle);
      //! Constrain the torsion angle a-b-c-d
      void AddTorsionConstraint(OBAtom *a, OBAtom *b, OBAtom *c, OBAtom *d, double torsion);
      //@}
      /////////////////////////////////////////////////////////////////////////
      // Get Constraints                                                     //
      /////////////////////////////////////////////////////////////////////////
      //! \name Methods to get information about set constraints
      //@{
      //! \returns the number of set constraints
      int Size() const;
      /*! The following constraint types are known: OBFF_CONST_IGNORE (ignore 
       *  the atom while setting up calculations, forcefield implementations 
       *  need to check this value in their setup function), OBFF_CONST_ATOM
       *  (fix atom position), OBFF_CONST_ATOM_X (fix x coordinate), 
       *  OBFF_CONST_ATOM_Y (fix y coordinate), OBFF_CONST_ATOM_Z (fix z 
       *  coordinate), OBFF_CONST_BOND (constrain bond length), OBFF_CONST_ANGLE
       *  (constrain angle), OBFF_CONST_TORSION (constrain torsion angle)
       *  \return the constraint type
       */
      int GetConstraintType(int index) const;
      /*! \return The constraint value, this can be a bond length, angle or 
       *   torsion angle depending on the constraint type.
       */
      double GetConstraintValue(int index);
      //! \return The constraint atom a (or fixed atom)
      //! \par index constraint index
      OBAtom* GetConstraintAtomA(int index);
      //! \return The constraint atom b
      //! \par index constraint index
      OBAtom* GetConstraintAtomB(int index);
      //! \return The constraint atom c
      //! \par index constraint index
      OBAtom* GetConstraintAtomC(int index);
      //! \return The constraint atom d
      //! \par index constraint index
      OBAtom* GetConstraintAtomD(int index);
      //! \return true if this atom is ignored
      //! \par index atom index
      bool IsIgnored(int index);
      //! \return true if this atom is fixed
      //! \par index atom index
      bool IsFixed(int index);
      //! \return true if the x coordinate for this atom is fixed
      //! \par index atom index
      bool IsXFixed(int index);
      //! \return true if the y coordinate for this atom is fixed
      //! \par index atom index
      bool IsYFixed(int index);
      //! \return true if the z coordinate for this atom is fixed
      //! \par index atom index
      bool IsZFixed(int index);
      //@}
 
    private:
      std::vector<OBFFConstraint> _constraints;
  };
 
  // Class OBForceField
  // class introduction in forcefield.cpp
  class OBFPRT  OBForceField : public OBPlugin
  {
  
    MAKE_PLUGIN(OBForceField)
  
    public:
    //!Clone the current instance. May be desirable in multithreaded environments,
    //!Should be deleted after use
    virtual OBForceField* MakeNewInstance()=0;

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
    //! Get index for vector<OBFFParameter> ...
    int GetParameterIdx(int a, int b, int c, int d, std::vector<OBFFParameter> &parameter);
           
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
 
    /*! Find the first arom in a 1-(level+1) relationship
     *  \param atom atom 1
     *  \param level the 1-(level+1) relationship (1, 2 or 3)
     *  \return index for the atom with 1-(level+1) relationship
     */
    int get_nbr (OBAtom* atom, int level);
    //    bool is14(OBAtom *a, OBAtom *b);
    // use OBAtom::IsOneFour(b)
      
    std::vector<int> _ignore; //!< List of atoms that are ignored while setting up calculations
    std::vector<int> _fix; //!< List of atoms that are fixed while minimizing
    OBMol _mol; //!< Molecule to be evaluated or minimized

    OBFFConstraints _constraints; //!< Constraints

    //! Output for logfile
    std::ostream* logos;
    char logbuf[BUFF_SIZE]; //!< Temporary buffer for logfile output
    int loglvl; //!< Log level for output
    int _origLogLevel;
    
    //! used to hold i for current conformer (needed by UpdateConformers)
    int _current_conformer;
    //! used to hold the energies for all conformers
    std::vector<double> _energies;

    //! Used for conjugate gradients and steepest descent(Initialize and TakeNSteps)
    double _econv, _e_n1;
    int _method, _cstep, _nsteps;
    double *_grad1, *_dir1;
    int _ncoords; //!< Number of coordinates for conjugate gradients

  public:
    //! Destructor
    virtual ~OBForceField()
      {
        if (_grad1 != NULL) {
          delete [] _grad1;
          _grad1 = NULL;
        }
        if (_dir1 != NULL) {
          delete [] _dir1;
          _dir1 = NULL;
        }
      }
    const char* TypeID()
      {
        return "forcefields";
      }

    /*! \param ID forcefield id (Ghemical, ...)
     *  \return A pointer to a forcefield (the default if ID is empty), or NULL if not available
     */
    static OBForceField* FindForceField(const std::string& ID)
    { 
      return FindType(ID.c_str());
    } 
    /*! \param ID forcefield id (Ghemical, ...)
     *  \return A pointer to a forcefield (the default if ID is empty), or NULL if not available
     */
    static OBForceField* FindForceField(const char *ID)
    {
      return FindType(ID);
    }
    //! \return The unit (kcal/mol, kJ/mol, ...) in which the energy is expressed as std::string
    virtual std::string GetUnit() { return std::string("au"); }
    /*! Setup the forcefield for mol (assigns atom types, charges, etc.) 
     *  \param mol the OBMol object that contains the atoms and bonds
     *  \return True if succesfull
     */
    virtual bool Setup(OBMol &mol) { return false; }
    /*! Compare the internal forcefield OBMol object to mol. If the two have the
     *  same number of atoms and bonds, and all atomic numbers are the same, 
     *  this function returns false, and no call to Setup is needed.
     *  \return true if Setup needs to be called
     */
    bool IsSetupNeeded(OBMol &mol);
    /*! Get coordinates for current conformer
     *  \param mol the OBMol object to copy the coordinates to (from OBForceField::_mol)
     *  \return true if succesfull
     */
    bool GetCoordinates(OBMol &mol);
    bool UpdateCoordinates(OBMol &mol) {GetCoordinates(mol); } // = GetCoordinates, depricated
    /*! Get coordinates for all conformers
     *  \param mol the OBMol object to copy the coordinates to (from OBForceField::_mol)
     *  \return true if succesfull
     */
    bool GetConformers(OBMol &mol);
    bool UpdateConformers(OBMol &mol) { GetConformers(mol); } // = GetConformers, depricated
    /*! Set coordinates for current conformer
     *  \param mol the OBMol object to copy the coordinates from (to OBForceField::_mol)
     *  \return true if succesfull
     */
    bool SetCoordinates(OBMol &mol);
    /*! Set coordinates for all conformers
     *  \param mol the OBMol object to copy the coordinates from (to OBForceField::_mol)
     *  \return true if succesfull
     */
    bool SetConformers(OBMol &mol);
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
    //! Generate coordinates for the molecule (distance geometry). (OB 3.0)
    void DistanceGeometry();
    //! Generate coordinates for the molecule (knowledge based, energy minimization). (OB 3.0)
    void GenerateCoordinates();
    /*! Generate conformers for the molecule (systematicaly rotating torsions).
     *  
     *  The initial starting structure here is important, this structure should be
     *  minimized for the best results. SystematicRotorSearch works by rotating around
     *  the rotatable bond in a molecule (see OBRotamerList class). This rotating generates 
     *  multiple conformers. The energy for all these conformers is then evaluated and the 
     *  lowest energy conformer is selected.
     *
     * \param geomSteps the number of steps to take during geometry optimization
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
    void SystematicRotorSearch(unsigned int geomSteps = 2500);
    //! \return the number of conformers
    int SystematicRotorSearchInitialize(unsigned int geomSteps = 2500);
    //! \return true if there are more conformers
    bool SystematicRotorSearchNextConformer(unsigned int geomSteps = 2500);

    /*! Generate conformers for the molecule (randomly rotating torsions).
     *  
     *  The initial starting structure here is important, this structure should be
     *  minimized for the best results. RandomRotorSearch works by randomly rotating around
     *  the rotatable bonds in a molecule (see OBRotamerList class). This rotating generates 
     *  multiple conformers. The energy for all these conformers is then evaluated and the 
     *  lowest energy conformer is selected.
     *
     * \param conformers the number of random conformers to consider during the search
     * \param geomSteps the number of steps to take during geometry optimization for each conformer
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
    void RandomRotorSearch(unsigned int conformers, unsigned int geomSteps = 2500);
    //! Initialize Random Rotor Search
    void RandomRotorSearchInitialize(unsigned int conformers, unsigned int geomSteps = 2500);
    //! \return true if there are more conformers
    bool RandomRotorSearchNextConformer(unsigned int geomSteps = 2500);
     
    /*! Generate conformers for the molecule (randomly rotating torsions).
     *  
     *  The initial starting structure here is important, this structure should be
     *  minimized for the best results. WeightedRotorSearch works by randomly rotating around
     *  the rotatable bonds in a molecule (see OBRotamerList class). Unlike RandomRotorSearch()
     *  the random choice of torsions is reweighted based on the energy of the generated conformer.
     *  Over time, the generated conformers for each step should become increasingly better.
     *  The lowest energy conformer is selected.
     *
     * \param conformers the number of random conformers to consider during the search
     * \param geomSteps the number of steps to take during geometry optimization for each conformer
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
    void WeightedRotorSearch(unsigned int conformers, unsigned int geomSteps);

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

    /*! Perform a linesearch for the entire molecule in direction direction
      \param currentCoords start coordinates
      \param direction the search direction

      \return alpha, the scale of the step we moved along the direction vector

      \par Output to log:
        OBFF_LOGLVL_NONE:   none \n
        OBFF_LOGLVL_LOW:    none \n
        OBFF_LOGLVL_MEDIUM: none \n
        OBFF_LOGLVL_HIGH:   none \n
    */
    double LineSearch(double *currentCoords, double *direction);

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
    //! (debugging)
    bool DetectExplosion();
    //! (debugging)
    vector3 ValidateLineSearch(OBAtom *atom, vector3 &direction);
    //! (debugging)
    void ValidateSteepestDescent(int steps);
    //! (debugging)
    void ValidateConjugateGradients(int steps);
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
    /*! Calculate the derivative of a OOP angle a-b-c-d. b is the central atom, and a-b-c is the plane. 
     * The OOP angle is given by 90° - arccos(dot(corss(ab,cb),db)/rabbc*rdb).
     * \param a atom a (coordinates), will be changed to -dtheta/da
     * \param b atom b (coordinates), will be changed to -dtheta/db
     * \param c atom c (coordinates), will be changed to -dtheta/dc
     * \param d atom d (coordinates), will be changed to -dtheta/dd
     * \return The OOP angle for a-b-c-d
     */
    static double VectorOOPDerivative(vector3 &a, vector3 &b, vector3 &c, vector3 &d);
 
    /*! Calculate the derivative of a torsion angle a-b-c-d. The torsion angle is given by arccos(dot(corss(ab,bc),cross(bc,cd))/rabbc*rbccd).
     * \param a atom a (coordinates), will be changed to -dtheta/da
     * \param b atom b (coordinates), will be changed to -dtheta/db
     * \param c atom c (coordinates), will be changed to -dtheta/dc
     * \param d atom d (coordinates), will be changed to -dtheta/dd
     * \return The tosion angle for a-b-c-d
     */
    static double VectorTorsionDerivative(vector3 &a, vector3 &b, vector3 &c, vector3 &d);
    //@}

  }; // class OBForceField

}// namespace OpenBabel

#endif   // OB_FORCEFIELD_H

//! \file forcefield.h
//! \brief Handle forcefields

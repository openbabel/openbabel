/**********************************************************************
forcefield.h - Handle OBForceField class.

Copyright (C) 2006-2007 by Tim Vandermeersch <tim.vandermeersch@gmail.com>

This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

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

#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>  // TODO: Move OBMol code out of the header (use OBMol*)
#include <openbabel/atom.h> // TODO: Move OBAtom code out of the header
#include <openbabel/plugin.h>
#include <openbabel/bitvec.h>
#include <float.h>

namespace OpenBabel
{
  class OBGridData;

  // log levels
#define OBFF_LOGLVL_NONE	0   //!< no output
#define OBFF_LOGLVL_LOW		1   //!< SteepestDescent progress... (no output from Energy())
#define OBFF_LOGLVL_MEDIUM	2   //!< individual energy terms
#define OBFF_LOGLVL_HIGH	3   //!< individual calculations and parameters

  // terms
#define OBFF_ENERGY		(1 << 0)   //!< all terms
#define OBFF_EBOND		(1 << 1)   //!< bond term
#define OBFF_EANGLE		(1 << 2)   //!< angle term
#define OBFF_ESTRBND		(1 << 3)   //!< strbnd term
#define OBFF_ETORSION		(1 << 4)   //!< torsion term
#define OBFF_EOOP		(1 << 5)   //!< oop term
#define OBFF_EVDW		(1 << 6)   //!< vdw term
#define OBFF_EELECTROSTATIC	(1 << 7)   //!< electrostatic term

  // constraint types
#define OBFF_CONST_IGNORE	(1 << 0)   //!< ignore the atom while setting up calculations
#define OBFF_CONST_ATOM		(1 << 1)   //!< fix the atom position
#define OBFF_CONST_ATOM_X	(1 << 2)   //!< fix the x coordinate of the atom position
#define OBFF_CONST_ATOM_Y	(1 << 3)   //!< fix the y coordinate of the atom position
#define OBFF_CONST_ATOM_Z	(1 << 4)   //!< fix the z coordinate of the atom position
#define OBFF_CONST_DISTANCE	(1 << 5)   //!< constrain distance length
#define OBFF_CONST_ANGLE	(1 << 6)   //!< constrain angle
#define OBFF_CONST_TORSION	(1 << 7)   //!< constrain torsion
#define OBFF_CONST_CHIRAL	(1 << 8)   //!< constrain chiral volume

  // mode arguments for SteepestDescent, ConjugateGradients, ...
#define OBFF_NUMERICAL_GRADIENT  	(1 << 0)  //!< use numerical gradients
#define OBFF_ANALYTICAL_GRADIENT	(1 << 1)  //!< use analytical gradients

const double KCAL_TO_KJ = 4.1868;
const double GAS_CONSTANT = 8.31446261815324e-3 / KCAL_TO_KJ;  //!< kcal mol^-1 K^-1 (2018 CODATA recommended value)

  // inline if statements for logging.
#define IF_OBFF_LOGLVL_LOW    if(_loglvl >= OBFF_LOGLVL_LOW)
#define IF_OBFF_LOGLVL_MEDIUM if(_loglvl >= OBFF_LOGLVL_MEDIUM)
#define IF_OBFF_LOGLVL_HIGH   if(_loglvl >= OBFF_LOGLVL_HIGH)

  //! The type of line search to be used for optimization -- simple or Newton numeric
  struct LineSearchType
  {
    enum {
      Simple, Newton2Num
    };
  };
  /*
  struct ConstraintType
  {
    enum {
      Ignore, Atom, AtomX, AtomY, AtomZ, Distance, Angle, Torsion, Chiral
    };
  };
  */

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

  //! \class OBFFCalculation2 forcefield.h <openbabel/forcefield.h>
  //! \brief Internal class for OBForceField to hold energy and gradient calculations on specific force fields
  class OBFPRT OBFFCalculation2
  {
  public:
    //! Used to store the energy for this OBFFCalculation
    double energy;
    //! Used to store the atoms for this OBFFCalculation
    OBAtom *a, *b;
    //! Used to store the index of atoms for this OBFFCalculation
    int idx_a, idx_b;
    //! Pointer to atom coordinates as double[3]
    double *pos_a, *pos_b;
    //! Pointer to atom forces
    double force_a[3], force_b[3];
    //! Destructor
    virtual ~OBFFCalculation2()
    {
    }
    //! \return Setup pointers to atom positions and forces (To be called
    //!  while setting up calculations). Sets optimized to true.
    virtual void SetupPointers()
    {
      if (!a || !b) return;
      pos_a = a->GetCoordinate();
      idx_a = a->GetIdx();
      pos_b = b->GetCoordinate();
      idx_b = b->GetIdx();
    }
  };

  //! \class OBFFCalculation3 forcefield.h <openbabel/forcefield.h>
  //! \brief Internal class for OBForceField to hold energy and gradient calculations on specific force fields
  class OBFPRT OBFFCalculation3: public OBFFCalculation2
  {
  public:
    //! Used to store the atoms for this OBFFCalculation
    OBAtom *c;
    //! Used to store the index of atoms for this OBFFCalculation
    int idx_c;
    //! Pointer to atom coordinates as double[3]
    double *pos_c;
    //! Pointer to atom forces
    double force_c[3];
    //! Destructor
    virtual ~OBFFCalculation3()
    {
    }
    //! \return Setup pointers to atom positions and forces (To be called
    //!  while setting up calculations). Sets optimized to true.
    virtual void SetupPointers()
    {
      if (!a || !b || !c) return;
      pos_a = a->GetCoordinate();
      idx_a = a->GetIdx();
      pos_b = b->GetCoordinate();
      idx_b = b->GetIdx();
      pos_c = c->GetCoordinate();
      idx_c = c->GetIdx();
    }
  };

  //! \class OBFFCalculation4 forcefield.h <openbabel/forcefield.h>
  //! \brief Internal class for OBForceField to hold energy and gradient calculations on specific force fields
  class OBFPRT OBFFCalculation4: public OBFFCalculation3
  {
  public:
    //! Used to store the atoms for this OBFFCalculation
    OBAtom *d;
    //! Used to store the index of atoms for this OBFFCalculation
    int idx_d;
    //! Pointer to atom coordinates as double[3]
    double *pos_d;
    //! Pointer to atom forces
    double force_d[3];
    //! Destructor
    virtual ~OBFFCalculation4()
    {
    }
    //! \return Setup pointers to atom positions and forces (To be called
    //!  while setting up calculations). Sets optimized to true.
    void SetupPointers()
    {
      if (!a || !b || !c || !d) return;
      pos_a = a->GetCoordinate();
      idx_a = a->GetIdx();
      pos_b = b->GetCoordinate();
      idx_b = b->GetIdx();
      pos_c = c->GetCoordinate();
      idx_c = c->GetIdx();
      pos_d = d->GetCoordinate();
      idx_d = d->GetIdx();
    }
  };

  //! \class OBFFConstraint forcefield.h <openbabel/forcefield.h>
  //! \brief Internal class for OBForceField to hold constraints
  //! \since version 2.2
  class OBFPRT OBFFConstraint
  {
  public:
    //! Used to store the contraint energy for this OBFFConstraint
    double factor, constraint_value;
    double rab0, rbc0;
    //! Used to store the contraint type for this OBFFConstraint
    int type, ia, ib, ic, id;
    //! Used to store the atoms for this OBFFCostraint
    OBAtom *a, *b, *c, *d;
    //! Used to store the gradients for this OBFFCalculation
    vector3 grada, gradb, gradc, gradd;

    //! Constructor
    OBFFConstraint()
      {
        a = b = c = d = nullptr;
        ia = ib = ic = id = 0;
        constraint_value = 0.0;
        factor = 0.0;
      }
    //! Destructor
    ~OBFFConstraint()
      {
      }

    vector3 GetGradient(int a)
    {
      if (a == ia)
        return grada;
      else if (a == ib)
        return gradb;
      else if (a == ic)
        return gradc;
      else if (a == id)
        return gradd;
      else
        return  VZero;
    }
  };

  //! \class OBFFConstraints forcefield.h <openbabel/forcefield.h>
  //! \brief Internal class for OBForceField to handle constraints
  //! \since version 2.2
  class OBFPRT OBFFConstraints
  {
  public:
    //! Constructor
    OBFFConstraints();
    //! Destructor
    ~OBFFConstraints()
      {
        _constraints.clear();
        _ignored.Clear();
        _fixed.Clear();
        _Xfixed.Clear();
        _Yfixed.Clear();
        _Zfixed.Clear();
      }
    //! Clear all constraints
    void Clear();
    //! Get the constraint energy
    double GetConstraintEnergy();
    //! Get the constraint gradient for atom with index a
    vector3 GetGradient(int a);
    //! Get the constrain gradient for the atom
    OBFFConstraints& operator=(const OBFFConstraints &ai)
      {
        if (this != &ai) {
          _constraints = ai._constraints;
          _ignored = ai._ignored;
          _fixed = ai._fixed;
          _Xfixed = ai._Xfixed;
          _Yfixed = ai._Yfixed;
          _Zfixed = ai._Zfixed;
        }
        return *this;
      }

    /*! Translate indices to OBAtom* objects, this function is called from OBForceField::Setup,
     *  this function doesn't have to be called from anywhere else.
     */
    void Setup(OBMol &mol);

    /////////////////////////////////////////////////////////////////////////
    // Set Constraints                                                     //
    /////////////////////////////////////////////////////////////////////////
    //! \name Methods to set constraints
    //@{
    //! Set Constraint factor
    void SetFactor(double factor);
    //! Ignore the atom while setting up calculations
    void AddIgnore(int a);
    //! Fix the position of an atom
    void AddAtomConstraint(int a);
    //! Fix the x coordinate of the atom position
    void AddAtomXConstraint(int a);
    //! Fix the y coordinate of the atom position
    void AddAtomYConstraint(int a);
    //! Fix the z coordinate of the atom position
    void AddAtomZConstraint(int a);
    //! Constrain the bond length a-b
    void AddDistanceConstraint(int a, int b, double length);
    //! Constrain the angle a-b-c
    void AddAngleConstraint(int a, int b, int c, double angle);
    //! Constrain the torsion angle a-b-c-d
    void AddTorsionConstraint(int a, int b, int c, int d, double torsion);
    //! Delete a constraint
    //! \par index constraint index
    void DeleteConstraint(int index);
    //@}
    /////////////////////////////////////////////////////////////////////////
    // Get Constraints                                                     //
    /////////////////////////////////////////////////////////////////////////
    //! \name Methods to get information about set constraints
    //@{
    //! Get Constraint factor
    double GetFactor();
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
    double GetConstraintValue(int index) const;
    //! \return The constraint atom a (or fixed atom)
    //! \par index constraint index
    int GetConstraintAtomA(int index) const;
    //! \return The constraint atom b
    //! \par index constraint index
    int GetConstraintAtomB(int index) const;
    //! \return The constraint atom c
    //! \par index constraint index
    int GetConstraintAtomC(int index) const;
    //! \return The constraint atom d
    //! \par index constraint index
    int GetConstraintAtomD(int index) const;
    //! \return true if this atom is ignored
    //! \par a atom index
    bool IsIgnored(int a);
    //! \return true if this atom is fixed
    //! \par a atom index
    bool IsFixed(int a);
    //! \return true if the x coordinate for this atom is fixed
    //! \par a atom index
    bool IsXFixed(int a);
    //! \return true if the y coordinate for this atom is fixed
    //! \par a atom index
    bool IsYFixed(int a);
    //! \return true if the z coordinate for this atom is fixed
    //! \par a atom index
    bool IsZFixed(int a);
    //! \return the ignored atom indexes as bitvec. (used in
    //! OBForceField::Setup() to determine if a call to
    //! OBForceField::SetupCalculations() is needed)
    OBBitVec GetIgnoredBitVec() { return _ignored; }
    //! \return the fixed atom indexes as bitvec. (used in
    //! OBForceField::SystematicRotorSearch() and similar)
    OBBitVec GetFixedBitVec() { return _fixed; }
    //@}

  private:
    std::vector<OBFFConstraint> _constraints;
    OBBitVec	_ignored;
    OBBitVec	_fixed;
    OBBitVec	_Xfixed;
    OBBitVec	_Yfixed;
    OBBitVec	_Zfixed;
    double _factor;
  };

  // Class OBForceField
  // class introduction in forcefield.cpp
  class OBFPRT OBForceField : public OBPlugin
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
    OBFFParameter* GetParameter(const char* a, const char* b, const char* c, const char* d,
        std::vector<OBFFParameter> &parameter);
    //! Get index for vector<OBFFParameter> ...
    int GetParameterIdx(int a, int b, int c, int d, std::vector<OBFFParameter> &parameter);

    /*! Calculate the potential energy function derivative numerically with
     *  repect to the coordinates of atom with index a (this vector is the gradient)
     *
     * \param a  provides coordinates
     * \param terms OBFF_ENERGY, OBFF_EBOND, OBFF_EANGLE, OBFF_ESTRBND, OBFF_ETORSION,
     * OBFF_EOOP, OBFF_EVDW, OBFF_ELECTROSTATIC
     * \return the negative gradient of atom a
     */
    vector3 NumericalDerivative(OBAtom *a, int terms = OBFF_ENERGY);
    //! OB 3.0
    vector3 NumericalSecondDerivative(OBAtom *a, int terms = OBFF_ENERGY);

    /*
     *   NEW gradients functions
     */

    /*! Set the gradient for atom with index idx to grad
     */
    void SetGradient(double *grad, int idx)
    {
      const int coordIdx = (idx - 1) * 3;
      for (unsigned int i = 0; i < 3; ++i) {
        _gradientPtr[coordIdx + i] = grad[i];
      }
    }

    /*! Add grad to the gradient for atom with index idx
     */
    void AddGradient(double *grad, int idx)
    {
      const int coordIdx = (idx - 1) * 3;
      for (unsigned int i = 0; i < 3; ++i) {
        _gradientPtr[coordIdx + i] += grad[i];
      }
    }

    /*! Set all gradients to zero
     */
    virtual void ClearGradients()
    {
      // We cannot use memset because IEEE floating point representations
      // are not guaranteed by C/C++ standard, but this loop can be
      // unrolled or vectorized by compilers
      for (unsigned int i = 0; i < _ncoords; ++i)
        _gradientPtr[i] = 0.0;
      //      memset(_gradientPtr, '\0', sizeof(double)*_ncoords);
    }

    /*! Check if two atoms are in the same ring. [NOTE: this function uses SSSR,
     *  this means that not all rings are found for bridged rings. This causes
     *  some problems with the MMFF94 validation.]
     *  \param a atom a
     *  \param b atom b
     *  \return true if atom a and b are in the same ring
     */
    bool IsInSameRing(OBAtom* a, OBAtom* b);

    // general variables
    OBMol 	_mol; //!< Molecule to be evaluated or minimized
    bool 	_init; //!< Used to make sure we only parse the parameter file once, when needed
    std::string	_parFile; //! < parameter file name
    bool 	_validSetup; //!< was the last call to Setup successful
    double	*_gradientPtr; //!< pointer to the gradients (used by AddGradient(), minimization functions, ...)
    // logging variables
    std::ostream* _logos; //!< Output for logfile
    char 	_logbuf[BUFF_SIZE+1]; //!< Temporary buffer for logfile output
    int 	_loglvl; //!< Log level for output
    int 	_origLogLevel;
    // conformer genereation (rotor search) variables
    int 	_current_conformer; //!< used to hold i for current conformer (needed by UpdateConformers)
    std::vector<double> _energies; //!< used to hold the energies for all conformers
    // minimization variables
    double 	_econv, _gconv, _e_n1; //!< Used for conjugate gradients and steepest descent(Initialize and TakeNSteps)
    int 	_cstep, _nsteps; //!< Used for conjugate gradients and steepest descent(Initialize and TakeNSteps)
    double 	*_grad1; //!< Used for conjugate gradients and steepest descent(Initialize and TakeNSteps)
    unsigned int _ncoords; //!< Number of coordinates for conjugate gradients
    int         _linesearch; //!< LineSearch type
    // molecular dynamics variables
    double 	_timestep; //!< Molecular dynamics time step in picoseconds
    double 	_temp; //!< Molecular dynamics temperature in Kelvin
    double 	*_velocityPtr; //!< pointer to the velocities
    // contraint varibles
    static OBFFConstraints _constraints; //!< Constraints
    static unsigned int _fixAtom; //!< SetFixAtom()/UnsetFixAtom()
    static unsigned int _ignoreAtom; //!< SetIgnoreAtom()/UnsetIgnoreAtom()
    // cut-off variables
    bool 	_cutoff; //!< true = cut-off enabled
    double 	_rvdw; //!< VDW cut-off distance
    double 	_rele; //!< Electrostatic cut-off distance
    double _epsilon; //!< Dielectric constant for electrostatics
    OBBitVec	_vdwpairs; //!< VDW pairs that should be calculated
    OBBitVec	_elepairs; //!< Electrostatic pairs that should be calculated
    int 	_pairfreq; //!< The frequence to update non-bonded pairs
    // group variables
    std::vector<OBBitVec> _intraGroup; //!< groups for which intra-molecular interactions should be calculated
    std::vector<OBBitVec> _interGroup; //!< groups for which intra-molecular interactions should be calculated
    std::vector<std::pair<OBBitVec, OBBitVec> > _interGroups; //!< groups for which intra-molecular
                                                              //!< interactions should be calculated
  public:
    /*! Clone the current instance. May be desirable in multithreaded environments,
     *  Should be deleted after use
     */
    virtual OBForceField* MakeNewInstance()=0;

    //! Destructor
    virtual ~OBForceField()
    {
      if (_grad1 != nullptr) {
        delete [] _grad1;
        _grad1 = nullptr;
      }
      if (_gradientPtr != nullptr) {
        delete [] _gradientPtr;
	_gradientPtr = nullptr;
      }
    }

    //! \return Plugin type ("forcefields")
    const char* TypeID()
    {
      return "forcefields";
    }

    /*! \param ID forcefield id (Ghemical, MMFF94, UFF, ...).
     *  \return A pointer to a forcefield (the default if ID is empty), or NULL if not available.
     */
    static OBForceField* FindForceField(const std::string& ID)
    {
      return FindType(ID.c_str());
    }
    /*! \param ID forcefield id (Ghemical, MMFF94, UFF, ...).
     *  \return A pointer to a forcefield (the default if ID is empty), or NULL if not available.
     */
    static OBForceField* FindForceField(const char *ID)
    {
      return FindType(ID);
    }
    /*
     *
     */
    void SetParameterFile(const std::string &filename)
    {
      _parFile = filename;
      _init = false;
    }
    /*! \return The unit (kcal/mol, kJ/mol, ...) in which the energy is expressed as std::string.
     */
    virtual std::string GetUnit() { return std::string("au"); }
    /* Does this force field have analytical gradients defined for all
     * calculation components (bonds, angles, non-bonded, etc.)
     * If this is true, code should default to using OBFF_ANALYTICAL_GRADIENT
     * for SteepestDescent() or ConjugateGradients().
     * \return True if all analytical gradients are implemented.
     */
    virtual bool HasAnalyticalGradients() { return false; }
    /*! Setup the forcefield for mol (assigns atom types, charges, etc.). Keep current constraints.
     *  \param mol The OBMol object that contains the atoms and bonds.
     *  \return True if successful.
     */
    bool Setup(OBMol &mol);
    /*! Setup the forcefield for mol (assigns atom types, charges, etc.). Use new constraints.
     *  \param mol The OBMol object that contains the atoms and bonds.
     *  \param constraints The OBFFConstraints object that contains the constraints.
     *  \return True if successful.
     */
    bool Setup(OBMol &mol, OBFFConstraints &constraints);
    /*! Load the parameters (this function is overloaded by the individual forcefields,
     *  and is called autoamically from OBForceField::Setup()).
     */
    // move to protected in future version
    virtual bool ParseParamFile() { return false; }
    /*! Set the atom types (this function is overloaded by the individual forcefields,
     *  and is called autoamically from OBForceField::Setup()).
     */
    // move to protected in future version
    virtual bool SetTypes() { return false; }
    /*! Set the formal charges (this function is overloaded by the individual forcefields,
     *  and is called autoamically from OBForceField::Setup()).
     */
    // move to protected in future version
    virtual bool SetFormalCharges() { return false; }
    /*! Set the partial charges (this function is overloaded by the individual forcefields,
     *  and is called autoamically from OBForceField::Setup()).
     */
    // move to protected in future version
    virtual bool SetPartialCharges() { return false; }
    /*! Setup the calculations (this function is overloaded by the individual forcefields,
     *  and is called autoamically from OBForceField::Setup()).
     */
    // move to protected in future version
    virtual bool SetupCalculations() { return false; }
    /*! Setup the pointers to the atom positions in the OBFFCalculation objects. This method
     *  will iterate over all the calculations and call SetupPointers for each one. (This
     *  function should be implemented by the individual force field implementations).
     */
    // move to protected in future version
    virtual bool SetupPointers() { return false; }
    /*! Compare the internal forcefield OBMol object to mol. If the two have the
     *  same number of atoms and bonds, and all atomic numbers are the same,
     *  this function returns false, and no call to Setup is needed.
     *  \return True if Setup needs to be called.
     */
    bool IsSetupNeeded(OBMol &mol);
    /*! Get the force atom types. The atom types will be added to
     *  the atoms of mol as OBPairData. The attribute will be "FFAtomType".
     *
     *  \code
     *  ...
     *  pFF->Setup(&mol);
     *  pFF->GetAtomTypes(&mol);
     *  FOR_ATOMS_OF_MOL (atom, mol) {
     *    OBPairData *type = (OBPairData*) atom->GetData("FFAtomType");
     *    if (type)
     *      cout << "atom " << atom->GetIdx() << " : " << type->GetValue() << endl;
     *  }
     *  ...
     *  \endcode
     */
    bool GetAtomTypes(OBMol &mol);
    /*! Get the force field formal charges. The formal charges will be added to
     *  the atoms of mol as OBPairData. The attribute will be "FFPartialCharge".
     *
     *  \code
     *  ...
     *  pFF->Setup(&mol);
     *  pFF->GetPartialCharges(&mol);
     *  FOR_ATOMS_OF_MOL (atom, mol) {
     *    OBPairData *chg = (OBPairData*) atom->GetData("FFPartialCharge");
     *    if (chg)
     *      cout << "atom " << atom->GetIdx() << " : " << chg->GetValue() << endl;
     *  }
     *  ...
     *  \endcode
     */
    bool GetPartialCharges(OBMol &mol);



    /*! Get coordinates for current conformer and attach OBConformerData with energies, forces, ... to mol.
     *  \param mol The OBMol object to copy the coordinates to (from OBForceField::_mol).
     *  \return True if successful.
     */
    bool GetCoordinates(OBMol &mol);
    //! \deprecated Use GetCooordinates instead.
    bool UpdateCoordinates(OBMol &mol) {return GetCoordinates(mol); }
    /*! Get coordinates for all conformers and attach OBConformerData with energies, forces, ... to mol.
     *  \param mol The OBMol object to copy the coordinates to (from OBForceField::_mol).
     *  \return True if successful.
     */
    bool GetConformers(OBMol &mol);
    //! \deprecated Use GetConformers instead.
    bool UpdateConformers(OBMol &mol) { return GetConformers(mol); }
    /*! Set coordinates for current conformer.
     *  \param mol the OBMol object to copy the coordinates from (to OBForceField::_mol).
     *  \return true if successful.
     */
    bool SetCoordinates(OBMol &mol);
    /*! Set coordinates for all conformers.
     *  \param mol The OBMol object to copy the coordinates from (to OBForceField::_mol).
     *  \return True if successful.
     */
    bool SetConformers(OBMol &mol);
    /*! Create a grid with spacing @p step and @p padding. Place a probe atom of type probe at every grid point,
     *  calculate the energy and store it in the grid. These grids can then be used to create isosurfaces to
     *  identify locations where the probe atom has favourable interactions with the molecule.
     *  \param step The grid step size in A..
     *  \param padding The padding for the grid in A.
     *  \param type The force field atom type for the probe.
     *  \param pchg The partial charge for the probe atom.
     *  \return Pointer to the grid constaining the results.
     */
    OBGridData *GetGrid(double step, double padding, const char *type, double pchg);

    /////////////////////////////////////////////////////////////////////////
    // Interacting groups                                                  //
    /////////////////////////////////////////////////////////////////////////

    //! \name Methods for specifying interaction groups
    //@{
    /*! Enable intra-molecular interactions for group (bonds, angles, strbnd, torsions, oop).
     *  This function should be called before Setup().
     *  \param group OBBitVec with bits set for the indexes of the atoms which make up the group.
     */
    void AddIntraGroup(OBBitVec &group);
    /*! Enable inter-molecular interactions for group (non-bonded: vdw & ele).
     *  This function should be called before Setup().
     *  \param group OBBitVec with bits set for the indexes of the atoms which make up the group.
     */
    void AddInterGroup(OBBitVec &group);
    /*! Enable inter-molecular interactions between group1 and group2 (non-bonded: vdw & ele).
     *  Note that this function doesn't enable bonded interactions in either group. Non-bonded
     *  interactions in the groups itself are also not enabled.
     *  This function should be called before Setup().
     *  \param group1 OBBitVec with bits set for the indexes of the atoms which make up the first group.
     *  \param group2 OBBitVec with bits set for the indexes of the atoms which make up the second group.
     */
    void AddInterGroups(OBBitVec &group1, OBBitVec &group2);
    /*! Clear all previously specified groups.
     */
    void ClearGroups();
    /*! \return true if there are groups.
     */
    bool HasGroups();
    //@}

    /////////////////////////////////////////////////////////////////////////
    // Cut-off                                                             //
    /////////////////////////////////////////////////////////////////////////

    //! \name Methods for Cut-off distances
    //@{
    /*! Enable or disable Cut-offs. Cut-offs are disabled by default.
     *  \param enable Enable when true, disable when false.
     */
    void EnableCutOff(bool enable)
    {
      _cutoff = enable;
    }
    /*! \return True if Cut-off distances are used.
     */
    bool IsCutOffEnabled()
    {
      return _cutoff;
    }
    /*! Set the VDW cut-off distance to r. Note that this does not enable cut-off distances.
     *  \param r The VDW cut-off distance to be used in A.
     */
    void SetVDWCutOff(double r)
    {
      _rvdw = r;
    }
    /*! Get the VDW cut-off distance.
     *  \return The VDW cut-off distance in A.
     */
    double GetVDWCutOff()
    {
      return _rvdw;
    }
    /*! Set the Electrostatic cut-off distance to r. Note that this does not
     *  enable cut-off distances.
     *  \param r The electrostatic cut-off distance to be used in A.
     */
    void SetElectrostaticCutOff(double r)
    {
      _rele = r;
    }
    /*! Get the Electrostatic cut-off distance.
     *  \return The electrostatic cut-off distance in A.
     */
    double GetElectrostaticCutOff()
    {
      return _rele;
    }
    /*! Set the dielectric constant for electrostatic SetupCalculations
     * \param epsilon The relative permittivity to use (default = 1.0)
     */
     void SetDielectricConstant(double epsilon)
     {
       _epsilon = epsilon;
     }
     /* Get the dielectric permittivity used for electrostatic calculations
     * \rreturn The current relative permittivity
     */
     double GetDielectricConstant()
     {
       return _epsilon;
     }
    /*! Set the frequency by which non-bonded pairs are updated. Values from 10 to 20
     *  are recommended. Too low will decrease performance, too high will cause
     *  non-bonded interactions within cut-off not to be calculated.
     *  \param f The pair list update frequency.
     */
    void SetUpdateFrequency(int f)
    {
      _pairfreq = f;
    }
    /*! Get the frequency by which non-bonded pairs are updated.
     *  \return The pair list update frequency.
     */
    int GetUpdateFrequency()
    {
      return _pairfreq;
    }
    /*! Set the bits in _vdwpairs and _elepairs to 1 for interactions that
     *  are within cut-off distance. This function is called in minimizing
     *  algorithms such as SteepestDescent and ConjugateGradients.
     */
    void UpdatePairsSimple();

    //void UpdatePairsGroup(); TODO

    /*! Get the number of non-bonded pairs in _mol.
     *  \return The number of atom pairs (ignores cutoff)
     */
    unsigned int GetNumPairs();
    /*! Get the number of enabled electrostatic pairs in _mol.
     *  \return The number of pairs currently enabled (within cut-off distance)
     */
    unsigned int GetNumElectrostaticPairs();
    /*! Get the number of enabled VDW pairs in _mol.
     *  \return The number of pairs currently enabled (within cut-off distance)
     */
    unsigned int GetNumVDWPairs();
    /*! Set bits in range 0..._numpairs-1 to 1. Using this means there will
     *  be no cut-off. (not-working: see code for more information.
     */
    void EnableAllPairs()
    {
      // TODO: OBBitVec doesn't seem to be allocating it's memory correctly
      //_vdwpairs.SetRangeOn(0, _numpairs-1);
      //_elepairs.SetRangeOn(0, _numpairs-1);
    }
    //@}

    /*! Get the pointer to the gradients
     */
    virtual vector3 GetGradient(OBAtom *a, int /*terms*/ = OBFF_ENERGY)
    {
      const int coordIdx = (a->GetIdx() - 1) * 3;
      return _gradientPtr + coordIdx;
    }

    /*! Get the pointer to the gradients
     */
    double* GetGradientPtr()
    {
      return _gradientPtr;
    }

    /////////////////////////////////////////////////////////////////////////
    // Energy Evaluation                                                   //
    /////////////////////////////////////////////////////////////////////////

    //! \name Methods for energy evaluation
    //@{
    /*! \param gradients Set to true when the gradients need to be calculated
     *  (needs to be done before calling GetGradient()).
     *  \return Total energy.
     *   \par Output to log:
     *    OBFF_LOGLVL_NONE:   none \n
     *    OBFF_LOGLVL_LOW:    none \n
     *    OBFF_LOGLVL_MEDIUM: energy for individual energy terms \n
     *    OBFF_LOGLVL_HIGH:   energy for individual energy interactions \n
     */
    virtual double Energy(bool UNUSED(gradients) = true) { return 0.0f; }
    /*! \param gradients Set to true when the gradients need to be calculated
     *  (needs to be done before calling GetGradient()).
     *  \return Bond stretching energy.
     *   \par Output to log:
     *    see Energy()
     */
    virtual double E_Bond(bool UNUSED(gradients) = true) { return 0.0f; }
    /*! \param gradients Set to true when the gradients need to be calculated
     *  (needs to be done before calling GetGradient()).
     *  \return Angle bending energy.
     *  \par Output to log:
     *   see Energy()
     */
    virtual double E_Angle(bool UNUSED(gradients) = true) { return 0.0f; }
    /*! \param gradients Set to true when the gradients need to be calculated
     *  (needs to be done before calling GetGradient()).
     *  \return Stretch bending energy.
     *   \par Output to log:
     *    see Energy()
     */
    virtual double E_StrBnd(bool UNUSED(gradients) = true) { return 0.0f; }
    /*! \param gradients Set to true when the gradients need to be calculated
     *  (needs to be done before calling GetGradient()).
     *  \return Torsional energy.
     *    \par Output to log:
     *	  see Energy()
     */
    virtual double E_Torsion(bool UNUSED(gradients) = true) { return 0.0f; }
    /*! \param gradients Set to true when the gradients need to be calculated
     *  (needs to be done before calling GetGradient()).
     *  \return Out-Of-Plane bending energy.
     *   \par Output to log:
     *	  see Energy()
     */
    virtual double E_OOP(bool UNUSED(gradients) = true) { return 0.0f; }
    /*! \param gradients Set to true when the gradients need to be calculated
     *  (needs to be done before calling GetGradient()).
     *  \return Van der Waals energy.
     *   \par Output to log:
     *	  see Energy()
     */
    virtual double E_VDW(bool UNUSED(gradients) = true) { return 0.0f; }
    /*! \param gradients Set to true when the gradients need to be calculated
     *  (needs to be done before calling GetGradient()).
     *  \return Electrostatic energy.
     *   \par Output to log:
     *	  see Energy()
     */
    virtual double E_Electrostatic(bool UNUSED(gradients) = true) { return 0.0f; }
    //@}

    /////////////////////////////////////////////////////////////////////////
    // Logging                                                             //
    /////////////////////////////////////////////////////////////////////////

    //! \name Methods for logging
    //@{
    /*! Print the atom types to the log.
     */
    void PrintTypes();
    /*! Print the formal charges to the log (atom.GetPartialCharge(),
     *  MMFF94 FC's are not always int).
     */
    void PrintFormalCharges();
    /*! Print the partial charges to the log.
     */
    void PrintPartialCharges();
    /*! Print the velocities to the log.
     */
    void PrintVelocities();
    /*! Set the stream for logging (can also be &cout for logging to screen).
     *  \param pos Stream (when pos is 0, std::cout wil be used).
     *  \return True if successful.
     */
    bool SetLogFile(std::ostream *pos);
    /*! Set the log level (OBFF_LOGLVL_NONE, OBFF_LOGLVL_LOW, OBFF_LOGLVL_MEDIUM, OBFF_LOGLVL_HIGH).
     *  Inline if statements for logging are available:
     *  \code
     *  #define IF_OBFF_LOGLVL_LOW    if(_loglvl >= OBFF_LOGLVL_LOW)
     *  #define IF_OBFF_LOGLVL_MEDIUM if(_loglvl >= OBFF_LOGLVL_MEDIUM)
     *  #define IF_OBFF_LOGLVL_HIGH   if(_loglvl >= OBFF_LOGLVL_HIGH)
     *  \endcode
     *
     *  example:
     *  \code
     *  SetLogLevel(OBFF_LOGLVL_MEDIUM);
     *  IF_OBFF_LOGLVL_HIGH {
     *    OBFFLog("this text will NOT be logged...\n");
     *  }
     *
     *  IF_OBFF_LOGLVL_LOW {
     *    OBFFLog"this text will be logged...\n");
     *  }
     *
     *  IF_OBFF_LOGLVL_MEDIUM {
     *    OBFFLog("this text will also be logged...\n");
     *  }
     *  \endcode
     */
    bool SetLogLevel(int level);
    /*! \return The log level.
     */
    int GetLogLevel() { return _loglvl; }
    /*! Print msg to the logfile.
     *  \param msg The message to print.
     */
    void OBFFLog(std::string msg)
    {
      if (!_logos)
        return;

      *_logos << msg;
    }
    /*! Print msg to the logfile.
     *  \param msg The message to print.
     */
    void OBFFLog(const char *msg)
    {
      if (!_logos)
        return;

      *_logos << msg;
    }
    //@}

    /////////////////////////////////////////////////////////////////////////
    // Structure Generation                                                //
    /////////////////////////////////////////////////////////////////////////

    //! \name Methods for structure generation
    //@{
    //! Generate coordinates for the molecule (distance geometry)
    //! \deprecated Use OBDistanceGeometry class instead
    void DistanceGeometry();
    /*! Generate conformers for the molecule (systematicaly rotating torsions).
     *
     *  The initial starting structure here is important, this structure should be
     *  minimized for the best results. SystematicRotorSearch works by rotating around
     *  the rotatable bond in a molecule (see OBRotamerList class). This rotating generates
     *  multiple conformers. The energy for all these conformers is then evaluated and the
     *  lowest energy conformer is selected.
     *
     *  \param geomSteps The number of steps to take during geometry optimization.
     *  \param sampleRingBonds Whether to sample ring torsions.
     *
     *	\par Output to log:
     *  This function should only be called with the log level set to OBFF_LOGLVL_NONE or OBFF_LOGLVL_LOW. Otherwise
     *	too much information about the energy calculations needed for this function will interfere with the output for
     *	this function. \n\n
     *  OBFF_LOGLVL_NONE:   None. \n
     *  OBFF_LOGLVL_LOW:    Number of rotatable bonds, energies for the conformers, which one is the lowest, ... \n
     *  OBFF_LOGLVL_MEDIUM: See note above. \n
     *  OBFF_LOGLVL_HIGH:   See note above. \n
     */
    void SystematicRotorSearch(unsigned int geomSteps = 2500, bool sampleRingBonds = false);
    /*! Generate conformers for the molecule by systematicaly rotating torsions. To be used in combination with
     *  SystematicRotorSearchNexConformer().
     *
     *  example:
     *  \code
     *  // pFF is a pointer to a OBForceField class
     *  pFF->SystematicRotorSearchInitialize(300);
     *  while (pFF->SystematicRotorSearchNextConformer(300)) {
     *    // do some updating in your program (show last generated conformer, ...)
     *  }
     *  \endcode
     *
     *  If you don't need any updating in your program, SystematicRotorSearch() is recommended.
     *
     *  \param geomSteps The number of steps to take during geometry optimization.
     *  \param sampleRingBonds Whether to sample ring torsions.
     *  \return The number of conformers.
     */
    int SystematicRotorSearchInitialize(unsigned int geomSteps = 2500, bool sampleRingBonds = false);
    /*! Evaluate the next conformer.
     *  \param geomSteps The number of steps to take during geometry optimization.
     *  \return True if there are more conformers.
     */
    bool SystematicRotorSearchNextConformer(unsigned int geomSteps = 2500);
    /*! Generate conformers for the molecule (randomly rotating torsions).
     *
     *  The initial starting structure here is important, this structure should be
     *  minimized for the best results. RandomRotorSearch works by randomly rotating around
     *  the rotatable bonds in a molecule (see OBRotamerList class). This rotating generates
     *  multiple conformers. The energy for all these conformers is then evaluated and the
     *  lowest energy conformer is selected.
     *
     *  \param conformers The number of random conformers to consider during the search.
     *  \param geomSteps The number of steps to take during geometry optimization for each conformer.
     *  \param sampleRingBonds Whether to sample ring torsions.
     *
     *	\par Output to log:
     *  This function should only be called with the log level set to OBFF_LOGLVL_NONE or OBFF_LOGLVL_LOW. Otherwise
     *	too much information about the energy calculations needed for this function will interfere with the output for
     *	this function. \n\n
     *  OBFF_LOGLVL_NONE:   None. \n
     *  OBFF_LOGLVL_LOW:    Number of rotatable bonds, energies for the conformers, which one is the lowest, ... \n
     *  OBFF_LOGLVL_MEDIUM: See note above. \n
     *  OBFF_LOGLVL_HIGH:   See note above. \n
     */
    void RandomRotorSearch(unsigned int conformers, unsigned int geomSteps = 2500,
                           bool sampleRingBonds = false);
    /*! Generate conformers for the molecule by randomly rotating torsions. To be used in combination with
     *  RandomRotorSearchNexConformer().
     *
     *  example:
     *  \code
     *  // pFF is a pointer to a OBForceField class
     *  pFF->RandomRotorSearchInitialize(300);
     *  while (pFF->RandomRotorSearchNextConformer(300)) {
     *    // do some updating in your program (show last generated conformer, ...)
     *  }
     *  \endcode
     *
     *  If you don't need any updating in your program, RandomRotorSearch() is recommended.
     *
     *  \param conformers The number of random conformers to consider during the search
     *  \param geomSteps The number of steps to take during geometry optimization
     *  \param sampleRingBonds Whether to sample ring torsions.
     */
    void RandomRotorSearchInitialize(unsigned int conformers, unsigned int geomSteps = 2500,
                                     bool sampleRingBonds = false);
    /*! Evaluate the next conformer.
     *  \param geomSteps The number of steps to take during geometry optimization.
     *  \return True if there are more conformers.
     */
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
     * \param conformers The number of random conformers to consider during the search.
     * \param geomSteps The number of steps to take during geometry optimization for each conformer.
     *  \param sampleRingBonds Whether to sample ring torsions.
     *
     *	\par Output to log:
     *  This function should only be called with the log level set to OBFF_LOGLVL_NONE or OBFF_LOGLVL_LOW. Otherwise
     *	too much information about the energy calculations needed for this function will interfere with the output for
     *	this function. \n\n
     *  OBFF_LOGLVL_NONE:   None. \n
     *  OBFF_LOGLVL_LOW:    Number of rotatable bonds, energies for the conformers, which one is the lowest, ... \n
     *  OBFF_LOGLVL_MEDIUM: See note above. \n
     *  OBFF_LOGLVL_HIGH:   See note above. \n
     */
    void WeightedRotorSearch(unsigned int conformers, unsigned int geomSteps,
                             bool sampleRingBonds = false);
    /**
     * @brief A fast rotor search to find low energy conformations
     *
     * Iterate over each of the rotors, and set the
     * torsion angle to that which minimizes the energy (while keeping the rest of the molecule
     * fixed). In general (for molecules with more than
     * one rotatable bond), this procedure will not find
     * the global minimum, but it will at least get rid of any bad
     * clashes, and it do so quickly.
     *
     * Torsions closer to the center
     * of the molecule will be optimized first as these most likely
     * to generate large clashes.
     *
     * One possible use of this procedure is to prepare a reasonable 3D structure
     * of a molecule for viewing. Another is to prepare the starting structure
     * for a more systematic rotor search (in which case you should geometry
     * optimize the final structure).
     *
     * @param permute Whether or not to permute the order of the 4 most central rotors.
     *                Default is true. This does a more thorough search, but takes 4! = 24 times
     *                as long.
     * @since version 2.4
     */
    int FastRotorSearch(bool permute = true);

#ifdef HAVE_EIGEN
    //! \since version 2.4
    int DiverseConfGen(double rmsd, unsigned int nconfs = 0, double energy_gap = 50, bool verbose = false);
#endif

    /////////////////////////////////////////////////////////////////////////
    // Energy Minimization                                                 //
    /////////////////////////////////////////////////////////////////////////

    //! \name Methods for energy minimization
    //@{
    /*! Set the LineSearchType. The default type is LineSearchType::Newton2Num.
     *  \param type The LineSearchType to be used in SteepestDescent and ConjugateGradients.
     */
    void SetLineSearchType(int type)
    {
      _linesearch = type;
    }
    /*! Get the LineSearchType.
     *  \return The current LineSearchType.
     */
    int GetLineSearchType()
    {
      return _linesearch;
    }
    /*! Perform a linesearch starting at atom in direction direction.
     * \deprecated Current code should use LineSearch(double *, double*) instead.
     */
    vector3 LineSearch(OBAtom *atom, vector3 &direction);
    /*! Perform a linesearch for the entire molecule in direction @p direction.
     *  This function is called when using LineSearchType::Simple.
     *
     *  \param currentCoords Start coordinates.
     *  \param direction The search direction.
     *  \return alpha, The scale of the step we moved along the direction vector.
     *
     *  \par Output to log:
     *  OBFF_LOGLVL_NONE:   none \n
     *  OBFF_LOGLVL_LOW:    none \n
     *  OBFF_LOGLVL_MEDIUM: none \n
     *  OBFF_LOGLVL_HIGH:   none \n
     */
    double LineSearch(double *currentCoords, double *direction);
    /*! Perform a linesearch for the entire molecule.
     *  This function is called when using LineSearchType::Newton2Num.
     *
     *  \param direction The search direction.
     *  \return alpha, The scale of the step we moved along the direction vector.
     *
     *  \par Output to log:
     *  OBFF_LOGLVL_NONE:   none \n
     *  OBFF_LOGLVL_LOW:    none \n
     *  OBFF_LOGLVL_MEDIUM: none \n
     *  OBFF_LOGLVL_HIGH:   none \n
     */
    double Newton2NumLineSearch(double *direction);
    /*! Set the coordinates of the atoms to origCoord + step.
     *  \param origCoords Start coordinates.
     *  \param direction The search direction.
     *  \param step The step to take.
     */
    void   LineSearchTakeStep(double *origCoords, double *direction, double step);
    /*! Perform steepest descent optimalization for steps steps or until convergence criteria is reached.
     *
     *  \param steps The number of steps.
     *  \param econv Energy convergence criteria. (defualt is 1e-6)
     *  \param method Deprecated. (see HasAnalyticalGradients())
     *
     *  \par Output to log:
     *  This function should only be called with the log level set to OBFF_LOGLVL_NONE or OBFF_LOGLVL_LOW. Otherwise
     *  too much information about the energy calculations needed for the minimization will interfere with the list
     *  of energies for succesive steps. \n\n
     *  OBFF_LOGLVL_NONE:   none \n
     *  OBFF_LOGLVL_LOW:    header including number of steps and first step \n
     *  OBFF_LOGLVL_MEDIUM: see note above \n
     *  OBFF_LOGLVL_HIGH:   see note above \n
     */
    void SteepestDescent(int steps, double econv = 1e-6f, int method = OBFF_ANALYTICAL_GRADIENT);
    /*! Initialize steepest descent optimalization, to be used in combination with SteepestDescentTakeNSteps().
     *
     *  example:
     *  \code
     *  // pFF is a pointer to a OBForceField class
     *  pFF->SteepestDescentInitialize(100, 1e-5f);
     *  while (pFF->SteepestDescentTakeNSteps(5)) {
     *    // do some updating in your program (redraw structure, ...)
     *  }
     *  \endcode
     *
     *  If you don't need any updating in your program, SteepestDescent() is recommended.
     *
     *  \param steps The number of steps.
     *  \param econv Energy convergence criteria. (defualt is 1e-6)
     *  \param method Deprecated. (see HasAnalyticalGradients())
     *
     *  \par Output to log:
     *  This function should only be called with the log level set to OBFF_LOGLVL_NONE or OBFF_LOGLVL_LOW. Otherwise
     *  too much information about the energy calculations needed for the minimization will interfere with the list
     *  of energies for succesive steps. \n\n
     *  OBFF_LOGLVL_NONE:   none \n
     *  OBFF_LOGLVL_LOW:    header including number of steps \n
     *  OBFF_LOGLVL_MEDIUM: see note above \n
     *  OBFF_LOGLVL_HIGH:   see note above \n
     */
    void SteepestDescentInitialize(int steps = 1000, double econv = 1e-6f, int method = OBFF_ANALYTICAL_GRADIENT);
    /*! Take n steps in a steepestdescent optimalization that was previously initialized with SteepestDescentInitialize().
     *
     *  \param n The number of steps to take.
     *  \return False if convergence or the number of steps given by SteepestDescentInitialize() has been reached.
     *
     *  \par Output to log:
     *  This function should only be called with the log level set to OBFF_LOGLVL_NONE or OBFF_LOGLVL_LOW. Otherwise
     *  too much information about the energy calculations needed for the minimization will interfere with the list
     *  of energies for succesive steps. \n\n
     *  OBFF_LOGLVL_NONE:   none \n
     *  OBFF_LOGLVL_LOW:    step number, energy and energy for the previous step \n
     *  OBFF_LOGLVL_MEDIUM: see note above \n
     *  OBFF_LOGLVL_HIGH:   see note above \n
     */
    bool SteepestDescentTakeNSteps(int n);
    /*! Perform conjugate gradient optimalization for steps steps or until convergence criteria is reached.
     *
     *  \param steps The number of steps.
     *  \param econv Energy convergence criteria. (defualt is 1e-6)
     *  \param method Deprecated. (see HasAnalyticalGradients())
     *
     *  \par Output to log:
     *  This function should only be called with the log level set to OBFF_LOGLVL_NONE or OBFF_LOGLVL_LOW. Otherwise
     *  too much information about the energy calculations needed for the minimization will interfere with the list
     *  of energies for succesive steps. \n\n
     *  OBFF_LOGLVL_NONE:   none \n
     *  OBFF_LOGLVL_LOW:    information about the progress of the minimization \n
     *  OBFF_LOGLVL_MEDIUM: see note above \n
     *  OBFF_LOGLVL_HIGH:   see note above \n
     */
    void ConjugateGradients(int steps, double econv = 1e-6f, int method = OBFF_ANALYTICAL_GRADIENT);
    /*! Initialize conjugate gradient optimalization and take the first step, to be
     *  used in combination with ConjugateGradientsTakeNSteps().
     *
     *  example:
     *  \code
     *  // pFF is a pointer to a OBForceField class
     *  pFF->ConjugateGradientsInitialize(100, 1e-5f);
     *  while (pFF->ConjugateGradientsTakeNSteps(5)) {
     *    // do some updating in your program (redraw structure, ...)
     *  }
     *  \endcode
     *
     *  If you don't need any updating in your program, ConjugateGradients() is recommended.
     *
     *  \param steps The number of steps.
     *  \param econv Energy convergence criteria. (defualt is 1e-6)
     *  \param method Deprecated. (see HasAnalyticalGradients())
     *
     *  \par Output to log:
     *  This function should only be called with the log level set to OBFF_LOGLVL_NONE or OBFF_LOGLVL_LOW. Otherwise
     *  too much information about the energy calculations needed for the minimization will interfere with the list
     *  of energies for succesive steps. \n\n
     *  OBFF_LOGLVL_NONE:   none \n
     *  OBFF_LOGLVL_LOW:    header including number of steps and first step \n
     *  OBFF_LOGLVL_MEDIUM: see note above \n
     *  OBFF_LOGLVL_HIGH:   see note above \n
     */
    void ConjugateGradientsInitialize(int steps = 1000, double econv = 1e-6f, int method = OBFF_ANALYTICAL_GRADIENT);
    /*! Take n steps in a conjugate gradient optimalization that was previously
     *  initialized with ConjugateGradientsInitialize().
     *
     *  \param n The number of steps to take.
     *  \return False if convergence or the number of steps given by ConjugateGradientsInitialize() has been reached.
     *
     *  \par Output to log:
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

    /////////////////////////////////////////////////////////////////////////
    // Molecular Dynamics                                                  //
    /////////////////////////////////////////////////////////////////////////

    //! \name Methods for molecular dynamics
    //@{
    /*! Generate starting velocities with a Maxwellian distribution.
     */
    void GenerateVelocities();
    /*! Correct the velocities so that the following is true:
     *
     *  \code
     *        3N
     *       ----
     *  0.5  \    m_i * v_i^2 = 0.5 * Ndf * R * T = E_kin
     *       /
     *       ----
     *       i=1
     *
     *  E_kin : kinetic energy
     *  m_i : mass of atom i
     *  v_i : velocity of atom i
     *  Ndf : number of degrees of freedom (3 * number of atoms)
     *  R : gas constant
     *  T : temperature
     *  \endcode
     *
     */
    void CorrectVelocities();
    /*! Take n steps at temperature T. If no velocities are set, they will be generated.
     *
     *  example:
     *  \code
     *  // pFF is a pointer to a OBForceField class
     *  while (pFF->MolecularDynamicsTakeNSteps(5, 300)) {
     *    // do some updating in your program (redraw structure, ...)
     *  }
     * \endcode
     *
     *  \param n The number of steps to take.
     *  \param T Absolute temperature in Kelvin.
     *  \param timestep The time step in picoseconds. (10e-12 s)
     *  \param method OBFF_ANALYTICAL_GRADIENTS (default) or OBFF_NUMERICAL_GRADIENTS
     */
    void MolecularDynamicsTakeNSteps(int n, double T, double timestep = 0.001, int method = OBFF_ANALYTICAL_GRADIENT);
    //@}

    /////////////////////////////////////////////////////////////////////////
    // Constraints                                                         //
    /////////////////////////////////////////////////////////////////////////

    //! \name Methods for constraints
    //@{
    /*! Get the current constraints.
     *  \return The current constrains stored in the force field.
     */
    OBFFConstraints& GetConstraints();
    /*! Set the constraints.
     *  \param constraints The new constraints to be used.
     */
    void SetConstraints(OBFFConstraints& constraints);
    /*! Fix the atom position until UnsetFixAtom() is called. This function
     *  can be used in programs that allow the user to interact with a molecule
     *  that is being minimized without having to check if the atom is already
     *  fixed in the constraints set by Setup() or SetConstraints(). Using this
     *  makes sure the selected atom follows the mouse cursur.
     *  \param index The index for the atom to fix.
     */
    void SetFixAtom(int index);
    /*! Undo last SetFixAtom. This function will not remove the fix atom
     *  constraint for this atom if set by Setup() or SetConstraints().
     */
    void UnsetFixAtom();
    /*! Ignore the atom until UnsetIgnoreAtom() is called. This function
     *  can be used in programs that allow the user to interact with a molecule
     *  that is being minimized without having to check if the atom is already
     *  ignored in the constraints set by Setup() or SetConstraints(). Using this
     *  makes sure, in drawing mode, you can close rings without your newly
     *  created puching the other atoms away.
     *  \param index The index for the atom to ignore.
     */
    void SetIgnoreAtom(int index);
    /*! Undo last SetIgnoreAtom. This function will not remove the ignore atom
     *  constraint for this atom if set by Setup() or SetConstraints().
     */
    void UnsetIgnoreAtom();

    //! internal function
    static bool IgnoreCalculation(int a, int b);
    //! internal function
    static bool IgnoreCalculation(int a, int b, int c);
    //! internal function
    static bool IgnoreCalculation(int a, int b, int c, int d);
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
     * \param pos_a atom a (coordinates)
     * \param pos_b atom b (coordinates)
     * \param force_a - return value for the force on atom a
     * \param force_b - return value for the force on atom b
     * \return The distance between a and b (bondlength for bond stretching, separation for vdw, electrostatic)
     */
    static double VectorBondDerivative(double *pos_a, double *pos_b,
                                       double *force_a, double *force_b);
    /*! To be used for VDW or Electrostatic interactions. This
     *  is faster than VectorBondDerivative, but does no error checking.
     */
    static double VectorDistanceDerivative(const double* const pos_i, const double* const pos_j,
                                           double *force_i, double *force_j);
    //! \deprecated
    static double VectorLengthDerivative(vector3 &a, vector3 &b);

    /*! Calculate the derivative of a angle a-b-c. The angle is given by dot(ab,cb)/rab*rcb.
     *  Used for harmonic (cubic) angle potentials.
     * \param pos_a atom a (coordinates)
     * \param pos_b atom b (coordinates)
     * \param pos_c atom c (coordinates)
     * \param force_a - return value for the force on atom a
     * \param force_b - return value for the force on atom b
     * \param force_c - return value for the force on atom c
     * \return The angle between a-b-c
     */
    static double VectorAngleDerivative(double *pos_a, double *pos_b, double *pos_c,
                                        double *force_a, double *force_b, double *force_c);
    //! \deprecated
    static double VectorAngleDerivative(vector3 &a, vector3 &b, vector3 &c);
    /*! Calculate the derivative of a OOP angle a-b-c-d. b is the central atom, and a-b-c is the plane.
     * The OOP angle is given by 90° - arccos(dot(corss(ab,cb),db)/rabbc*rdb).
     * \param pos_a atom a (coordinates)
     * \param pos_b atom b (coordinates)
     * \param pos_c atom c (coordinates)
     * \param pos_d atom d (coordinates)
     * \param force_a - return value for the force on atom a
     * \param force_b - return value for the force on atom b
     * \param force_c - return value for the force on atom c
     * \param force_d - return value for the force on atom d
     * \return The OOP angle for a-b-c-d
     */
    static double VectorOOPDerivative(double *pos_a, double *pos_b, double *pos_c, double *pos_d,
                                      double *force_a, double *force_b, double *force_c, double *force_d);
    //! \deprecated
    static double VectorOOPDerivative(vector3 &a, vector3 &b, vector3 &c, vector3 &d);
    /*! Calculate the derivative of a torsion angle a-b-c-d. The torsion angle is given by arccos(dot(corss(ab,bc),cross(bc,cd))/rabbc*rbccd).
     * \param pos_a atom a (coordinates)
     * \param pos_b atom b (coordinates)
     * \param pos_c atom c (coordinates)
     * \param pos_d atom d (coordinates)
     * \param force_a - return value for the force on atom a
     * \param force_b - return value for the force on atom b
     * \param force_c - return value for the force on atom c
     * \param force_d - return value for the force on atom d
     * \return The tosion angle for a-b-c-d
     */
    static double VectorTorsionDerivative(double *pos_a, double *pos_b, double *pos_c, double *pos_d,
                                          double *force_a, double *force_b, double *force_c, double *force_d);
    //! \deprecated
    static double VectorTorsionDerivative(vector3 &a, vector3 &b, vector3 &c, vector3 &d);

    /*! inline fuction to speed up minimization speed
     * \param i pointer to i[3]
     * \param j pointer to j[3]
     * \param result pointer to result[3], will be set to i - j
     */
    static void VectorSubtract(double *i, double *j, double *result)
    {
      for (unsigned int c = 0; c < 3; ++c)
        result[c] = i[c] - j[c];
    }

    static void VectorSubtract(const double* const i, const double* const j, double *result)
    {
      for (unsigned int c = 0; c < 3; ++c)
        result[c] = i[c] - j[c];
    }

    /*! inline fuction to speed up minimization speed
     * \param i pointer to i[3]
     * \param j pointer to j[3]
     * \param result pointer to result[3], will be set to i + j
     */
    static void VectorAdd(double *i, double *j, double *result)
    {
      for (unsigned int c = 0; c < 3; ++c)
        result[c] = i[c] + j[c];
    }

    /*! inline fuction to speed up minimization speed
     * \param i pointer to i[3]
     * \param n divide x,y,z with n
     * \param result pointer to result[3]
     */
    static void VectorDivide(double *i, double n, double *result)
    {
      for (unsigned int c = 0; c < 3; ++c)
        result[c] = i[c] / n;
    }

    /*! inline fuction to speed up minimization speed
     * \param i pointer to i[3]
     * \param n multiply x,y,z with n
     * \param result pointer to result[3]
     */
    static void VectorMultiply(double *i, double n, double *result)
    {
      for (unsigned int c = 0; c < 3; ++c)
        result[c] = i[c] * n;
    }

    static void VectorMultiply(const double* const i, const double n, double *result)
    {
      for (unsigned int c = 0; c < 3; ++c)
        result[c] = i[c] * n;
    }

    /*! inline fuction to speed up minimization speed
     * \param i pointer to i[3], multiply this vector by n and set this vector to the result.
     * \param n the scalar value to be multipled
     */
    static void VectorSelfMultiply(double *i, double n)
    {
      for (unsigned int c = 0; c < 3; ++c)
        i[c] *= n;
    }

    /*! inline fuction to speed up minimization speed
     * \param i pointer to i[3] to be normalized
     */
    static void VectorNormalize(double *i)
    {
      double length = VectorLength(i);
      for (unsigned int c = 0; c < 3; ++c)
        i[c] /= length;
    }

    /*! inline fuction to speed up minimization speed
     * \param from pointer to i[3] to be copied from
     * \param to pointer to j[3] to be copied to
     */
    static void VectorCopy(double *from, double *to)
    {
      for (unsigned int c = 0; c < 3; ++c)
        to[c] = from[c];
    }

    /*! inline fuction to speed up minimization speed
     * \param i pointer to i[3]
     * \return the vector length
     */
    static double VectorLength(double *i)
    {
      return sqrt( i[0]*i[0] + i[1]*i[1] + i[2]*i[2] );
    }

    static double VectorDistance(double *pos_i, double *pos_j)
    {
      double ij[3];
      VectorSubtract(pos_i, pos_j, ij);
      const double rij = VectorLength(ij);
      return rij;
    }

    /*! inline fuction to speed up minimization speed
     * \param i pointer to i[3]
     * \param j pointer to j[3]
     * \param k pointer to k[3]
     * \return the vector angle ijk (deg)
     */
    static double VectorAngle(double *i, double *j, double *k);

    /*! inline fuction to speed up minimization speed
     * \param i pointer to i[3]
     * \param j pointer to j[3]
     * \param k pointer to k[3]
     * \param l pointer to l[3]
     * \return the vector torson ijkl (deg)
     */
    static double VectorTorsion(double *i, double *j, double *k, double *l);

    /*! inline fuction to speed up minimization speed
     * \param i pointer to i[3]
     * \param j pointer to j[3]
     * \param k pointer to k[3]
     * \param l pointer to l[3]
     * \return the vector torson ijkl (deg)
     */
    static double VectorOOP(double *i, double *j, double *k, double *l);

    /*! inline fuction to speed up minimization speed
     * \param i pointer to i[3], will set x,y,z to 0,0,0
     */
    static void VectorClear(double *i)
    {
      for (unsigned int c = 0; c < 3; ++c)
        i[c] = 0.0;
    }

    /*! inline fuction to speed up minimization speed
     * \param i pointer to i[3]
     * \param j pointer to j[3]
     * \return the dot product
     */
    static double VectorDot(double *i, double *j)
    {
      double result = 0.0;
      // Written as a loop for vectorization
      // Loop will be unrolled by compiler otherwise
      for (unsigned int c = 0; c < 3; ++c)
        result += i[c]*j[c];
      return result;
    }

    /*! inline fuction to speed up minimization speed
     * \param i pointer to i[3]
     * \param j pointer to j[3]
     * \param result the dot product (as a return value double[3])
     */
    static void VectorCross(double *i, double *j, double *result)
    {
      result[0] =   i[1]*j[2] - i[2]*j[1];
      result[1] = - i[0]*j[2] + i[2]*j[0];
      result[2] =   i[0]*j[1] - i[1]*j[0];
    }

    static void PrintVector(double *i)
    {
      std::cout << "<" << i[0] << ", " << i[1] << ", " << i[2] << ">" << std::endl;
    }
    //@}

  }; // class OBForceField

}// namespace OpenBabel

#endif   // OB_FORCEFIELD_H

//! \file forcefield.h
//! \brief Handle forcefields

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
#include <openbabel/babelconfig.h>
#include <openbabel/base.h>
#include <openbabel/mol.h>
#include <openbabel/plugin.h>
#include <openbabel/grid.h>
#include <openbabel/griddata.h>
#include <openbabel/obfunction.h>
#include <float.h>

namespace OpenBabel
{
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

#define KCAL_TO_KJ	4.1868

  // inline if statements for logging.
#define IF_OBFF_LOGLVL_LOW    if(GetLogLevel() >= OBFF_LOGLVL_LOW)
#define IF_OBFF_LOGLVL_MEDIUM if(GetLogLevel() >= OBFF_LOGLVL_MEDIUM)
#define IF_OBFF_LOGLVL_HIGH   if(GetLogLevel() >= OBFF_LOGLVL_HIGH)

  //! The type of line search to be used for optimization -- simple or Newton numeric
  namespace LineSearchType 
  {
    enum {
      Simple, Newton2Num 
    };
  };
  
  class OBFFConstraint;
  
  /** @class OBFFConstraints forcefield.h <openbabel/forcefield.h>
      @brief Internal class for OBForceField to handle constraints
      @since version 2.2
   */
  class OBFPRT OBFFConstraints
  {
    friend class OBForceField;

    public:
      /** 
       * @brief Constructor.
       */
      OBFFConstraints();
      /** 
       * @brief Destructor.
       */
      ~OBFFConstraints();
      /**
       * @brief Assignment operator.
       */
      OBFFConstraints& operator=(const OBFFConstraints &ai); 
      /**
       *  @brief Clear all constraints
       */
      void Clear();
      /** 
       * @return The constraint energy.
       */
      double GetConstraintEnergy();
  
      /////////////////////////////////////////////////////////////////////////
      // Set Constraints                                                     //
      /////////////////////////////////////////////////////////////////////////
      //! \name Methods to set constraints
      //@{
      /** 
       * @brief Set Constraint factor.
       */
      void SetFactor(double factor);
      /** 
       * @brief Ignore the atom while setting up calculations.
       */
      void AddIgnore(int a);
      /** 
       * @brief Fix the position of an atom.
       */
      void AddAtomConstraint(int a);
      /** 
       * @brief Fix the x coordinate of the atom position.
       */
      void AddAtomXConstraint(int a);
      /** 
       * @brief Fix the y coordinate of the atom position.
       */
      void AddAtomYConstraint(int a);
      /** 
       * @brief Fix the z coordinate of the atom position.
       */
      void AddAtomZConstraint(int a);
      /** 
       * @brief Constrain the bond length a-b.
       */
      void AddDistanceConstraint(int a, int b, double length);
      /** 
       * @brief Constrain the angle a-b-c.
       */
      void AddAngleConstraint(int a, int b, int c, double angle);
      /** 
       * @brief Constrain the torsion angle a-b-c-d.
       */
      void AddTorsionConstraint(int a, int b, int c, int d, double torsion);
      /** 
       * @brief Delete a constraint.
       *
       * @param index constraint index.
       */
      void DeleteConstraint(int index);
      //@}
      /////////////////////////////////////////////////////////////////////////
      // Get Constraints                                                     //
      /////////////////////////////////////////////////////////////////////////
      //! \name Methods to get information about set constraints
      //@{
      /** 
       * @return The constraint factor.
       */
      double GetFactor();
      /** 
       * @returns the number of set constraints.
       */
      int Size() const;
      /** 
       * The following constraint types are known: OBFF_CONST_IGNORE (ignore 
       * the atom while setting up calculations, forcefield implementations 
       * need to check this value in their setup function), OBFF_CONST_ATOM
       * (fix atom position), OBFF_CONST_ATOM_X (fix x coordinate), 
       * OBFF_CONST_ATOM_Y (fix y coordinate), OBFF_CONST_ATOM_Z (fix z 
       * coordinate), OBFF_CONST_BOND (constrain bond length), OBFF_CONST_ANGLE
       * (constrain angle), OBFF_CONST_TORSION (constrain torsion angle)
       *  
       *  @return the constraint type
       */
      int GetConstraintType(unsigned int index) const;
      /** 
       * @return The constraint value, this can be a bond length, angle or 
       * torsion angle depending on the constraint type.
       */
      double GetConstraintValue(unsigned int index) const;
      /** 
       * @param index Constraint index.
       * 
       * @return The constraint atom a (or fixed atom)
       */
      int GetConstraintAtomA(unsigned int index) const;
      /** 
       * @param index Constraint index.
       * 
       * @return The constraint atom b.
       */
      int GetConstraintAtomB(unsigned int index) const;
      /** 
       * @param index Constraint index.
       * 
       * @return The constraint atom c.
       */
      int GetConstraintAtomC(unsigned int index) const;
      /** 
       * @param index Constraint index.
       * 
       * @return The constraint atom d.
       */
      int GetConstraintAtomD(unsigned int index) const;
      /** 
       * @param a Atom index.
       * 
       * @return True if this atom is ignored.
       */
      bool IsIgnored(unsigned int a);
      /** 
       * @param a Atom index.
       *
       * @return True if this atom is fixed.
       */
      bool IsFixed(unsigned int a);
      /** 
       * @param a Atom index.
       *
       * @return True if the x coordinate for this atom is fixed.
       */
      bool IsXFixed(unsigned int a);
      /** 
       * @param a Atom index.
       *
       * @return True if the y coordinate for this atom is fixed.
       */
      bool IsYFixed(unsigned int a);
      /** 
       * @param a Atom index.
       *
       * @return True if the z coordinate for this atom is fixed.
       */
      bool IsZFixed(unsigned int a);
      //@}

    protected: 
      /**
       * @return The constraint gradient for atom with index @p a.
       */
      Eigen::Vector3d GetGradient(int a);
      /** 
       * @brief Translate indices to OBAtom* objects, this function is called 
       * from OBForceField::Setup,this function doesn't have to be called 
       * from anywhere else.
       */
      void Setup(OBMol &mol);
      /** 
       * @return the ignored atom indexes as bitvec. (used in 
       * OBForceField::Setup() to determine if a call to 
       * OBForceField::SetupCalculations() is needed).
       */
      OBBitVec GetIgnoredBitVec() { return m_ignored; }
 
      std::vector<OBFFConstraint> m_constraints;
      OBBitVec m_ignored;
      OBBitVec m_fixed;
      OBBitVec m_Xfixed;
      OBBitVec m_Yfixed;
      OBBitVec m_Zfixed;
      double   m_factor;
  };
 
  // Class OBForceField
  // class introduction in forcefield.cpp
  class OBFFParameter;
  class OBForceFieldPrivate;
  class OBFPRT OBForceField : public OBPlugin, public OBFunction
  {
  
    // Plugin stuff, replaces MAKE_PLUGIN(OBForceField)
    /// @cond DEV
    protected:
      virtual PluginMapType& GetMap() const { return Map(); }
      static PluginMapType& Map() { static PluginMapType m; return m; }
    public:
      const char* TypeID() { return "forcefields"; }
      static OBForceField*& Default() { static OBForceField* d; return d; }
      OBForceField(const char* ID, bool IsDefault=false);
      static OBForceField* FindType(const char* ID);
      /*! Clone the current instance. May be desirable in multithreaded environments,
       *  Should be deleted after use
       */
      virtual OBForceField* MakeNewInstance() = 0;
    /// @endcond

    protected:
      bool m_cutoff; //!< true = cut-off enabled
      int m_loglvl; //!< Log level for output
      //std::vector<Eigen::Vector3d> m_positions;
      //std::vector<Eigen::Vector3d> m_gradients;
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
    /** 
     * see GetParameter(int a, int b, int c, int d, std::vector<OBFFParameter> &parameter).
     */
    OBFFParameter* GetParameter(const char* a, const char* b, const char* c, const char* d, 
        std::vector<OBFFParameter> &parameter);
    /** 
     * Get index for vector<OBFFParameter> ...
     */
    int GetParameterIdx(int a, int b, int c, int d, std::vector<OBFFParameter> &parameter);
    /** 
     * @brief Calculate the potential energy function derivative numerically with 
     * repect to the coordinates of atom with index a (this vector is the gradient)
     *
     * @param a  provides coordinates
     * @param terms OBFF_ENERGY, OBFF_EBOND, OBFF_EANGLE, OBFF_ESTRBND, OBFF_ETORSION, 
     * OBFF_EOOP, OBFF_EVDW, OBFF_ELECTROSTATIC
     * 
     * @return the negative gradient of atom a
     */
    Eigen::Vector3d NumericalDerivative(OBAtom *a, int terms = OBFF_ENERGY);
    //! OB 3.0
    Eigen::Vector3d NumericalSecondDerivative(OBAtom *a, int terms = OBFF_ENERGY);
    /*! 
     * @brief Add @p p grad to the gradient for atom with index @p idx
     */
    void AddGradient(const Eigen::Vector3d &grad, int idx);
    /** 
     * @brief Get the gradient for atom @p a.
     */
    Eigen::Vector3d GetGradient(OBAtom *a); 
    /**
     * @brief Set all gradients to zero
     */
    void ClearGradients(); 
    /** 
     * @brief Check if two atoms are in the same ring. [NOTE: this function uses SSSR, 
     * this means that not all rings are found for bridged rings. This causes 
     * some problems with the MMFF94 validation.]
     *  
     * @param a Atom a.
     * @param b Atom b.
     * 
     * @return True if atom a and b are in the same ring
     */
    bool IsInSameRing(OBAtom* a, OBAtom* b);
    /**
     * @brief Get the intra groups bitvec vector.
     */
    std::vector<OBBitVec>& GetIntraGroup();
    /**
     * @brief Get the inter groups bitvec vector.
     */ 
    std::vector<OBBitVec>& GetInterGroup();
    /**
     * @brief Get the inter groups bitvec vector.
     */
    std::vector<std::pair<OBBitVec, OBBitVec> >& GetInterGroups();
    /**
     * @brief Get the VDW pairs bitvec.
     */
    OBBitVec& GetVDWPairs();
    /**
     * @brief Get the Electrostatic pairs bitvec.
     */
    OBBitVec& GetElePairs();
    
    OBForceFieldPrivate * const d; //!< the d-pointer

  public:
    /** 
     * @brief Destructor.
     */
    virtual ~OBForceField();
    /** 
     * @param ID forcefield id (Ghemical, MMFF94, UFF, ...).
     * 
     * @return A pointer to a forcefield (the default if ID is empty), or 
     * NULL if not available.
     */
    static OBForceField* FindForceField(const std::string& ID)
    { 
      return FindType(ID.c_str());
    } 
    /** 
     * @param ID forcefield id (Ghemical, MMFF94, UFF, ...).
     * 
     * @return A pointer to a forcefield (the default if ID is empty), or 
     * NULL if not available.
     */
    static OBForceField* FindForceField(const char *ID)
    {
      return FindType(ID);
    }
    /**
     * @brief Set the parameter file
     */
    void SetParameterFile(const std::string &filename);
    /**
     * @brief Get the parameter file
     */
    std::string& GetParameterFile();
    /** 
     * @return The unit (kcal/mol, kJ/mol, ...) in which the energy is 
     * expressed as std::string.
     */  
    virtual std::string GetUnit() { return std::string("au"); }
    /** 
     * @brief Does this force field have analytical gradients defined for all
     * calculation components (bonds, angles, non-bonded, etc.)
     * If this is true, code should default to using OBFF_ANALYTICAL_GRADIENT
     * for SteepestDescent() or ConjugateGradients().
     * \return True if all analytical gradients are implemented.
     */
    virtual bool HasAnalyticalGradients() { return false; }
    /** 
     * @brief Setup the forcefield for mol (assigns atom types, charges, etc.). 
     * Keep current constraints.
     * 
     * @param mol The OBMol object that contains the atoms and bonds.
     * 
     * @return True if succesfull.
     */
    bool Setup(OBMol &mol); 
    /** 
     * @brief Setup the forcefield for mol (assigns atom types, charges, etc.). 
     * Use new constraints. 
     * 
     * @param mol The OBMol object that contains the atoms and bonds.
     * @param constraints The OBFFConstraints object that contains the constraints.
     * 
     * @return True if succesfull.
     */
    bool Setup(OBMol &mol, OBFFConstraints &constraints);
    /** 
     * @brief Load the parameters (this function is overloaded by the individual 
     * forcefields, and is called autoamically from OBForceField::Setup()).
     */
    // move to protected in future version
    virtual bool ParseParamFile() { return false; } 
    /** 
     * @brief Set the atom types (this function is overloaded by the individual 
     * forcefields, and is called autoamically from OBForceField::Setup()).
     */
    // move to protected in future version
    virtual bool SetTypes() { return false; }
    /** 
     * @brief Set the formal charges (this function is overloaded by the individual 
     * forcefields, and is called autoamically from OBForceField::Setup()).
     */
    // move to protected in future version
    virtual bool SetFormalCharges() { return false; }
    /** 
     * @brief Set the partial charges (this function is overloaded by the individual 
     * forcefields, and is called autoamically from OBForceField::Setup()).
     */
    // move to protected in future version
    virtual bool SetPartialCharges() { return false; }
    /** 
     * @brief Setup the calculations (this function is overloaded by the individual 
     * forcefields, and is called autoamically from OBForceField::Setup()).
     */
    // move to protected in future version
    virtual bool SetupCalculations() { return false; }
    /** 
     * @brief Compare the internal forcefield OBMol object to mol. If the two have the
     * same number of atoms and bonds, and all atomic numbers are the same, 
     * this function returns false, and no call to Setup is needed.
     * 
     * @return True if Setup needs to be called.
     */
    bool IsSetupNeeded(OBMol &mol);
    /**
     * @brief Get the molecule stored in the force field.
     */
    OBMol* GetMolecule();
    /** 
     * @brief Get the force atom types. The atom types will be added to 
     * the atoms of mol as OBPairData. The attribute will be "FFAtomType".
     *
     * @code
     * ...
     * pFF->Setup(&mol);
     * pFF->GetAtomTypes(&mol);
     * FOR_ATOMS_OF_MOL (atom, mol) {
     *   OBPairData *type = (OBPairData*) atom->GetData("FFAtomType");
     *   if (type)
     *     cout << "atom " << atom->GetIdx() << " : " << type->GetValue() << endl;
     * }
     * ...
     * @endcode
     */
    bool GetAtomTypes(OBMol &mol);
    /** 
     * @brief Get the force field formal charges. The formal charges will be 
     * added to the atoms of mol as OBPairData. The attribute will be 
     * "FFPartialCharge".
     *
     * @code
     * ...
     * pFF->Setup(&mol);
     * pFF->GetPartialCharges(&mol);
     * FOR_ATOMS_OF_MOL (atom, mol) {
     *   OBPairData *chg = (OBPairData*) atom->GetData("FFPartialCharge");
     *   if (chg)
     *     cout << "atom " << atom->GetIdx() << " : " << chg->GetValue() << endl;
     * }
     * ...
     * @endcode
     */
    bool GetPartialCharges(OBMol &mol);
    /** 
     * @brief Get coordinates for current conformer and attach OBConformerData 
     * with energies, forces, ... to mol.
     * 
     * @param mol The OBMol object to copy the coordinates to (from OBForceField::_mol).
     * 
     * @return True if succesfull.
     */
    bool GetCoordinates(OBMol &mol);
    /** 
     * @brief Get coordinates for all conformers and attach OBConformerData with 
     * energies, forces, ... to mol.
     * 
     * @param mol The OBMol object to copy the coordinates to (from OBForceField::_mol).
     * 
     * @return True if succesfull.
     */
    bool GetConformers(OBMol &mol);
    /** 
     * @brief Set coordinates for current conformer.
     * 
     * @param mol the OBMol object to copy the coordinates from (to OBForceField::_mol).
     * 
     * @return true if succesfull.
     */
    bool SetCoordinates(OBMol &mol);
    /** 
     * @brief Set coordinates for all conformers.
     *  
     * @param mol The OBMol object to copy the coordinates from (to OBForceField::_mol).
     * 
     * @return True if succesfull.
     */
    bool SetConformers(OBMol &mol);
    /** 
     * @brief Create a grid with spacing p step and p padding. Place a probe 
     * atom of type probe at every grid point, calculate the energy and store 
     * it in the grid. These grids can then be used to create isosurfaces to
     * identify locations where the probe atom has favourable interactions 
     * with the molecule.
     * 
     * @param step The grid step size in A..
     * @param padding The padding for the grid in A.
     * @param type The force field atom type for the probe.
     * @param pchg The partial charge for the probe atom.
     * 
     * @return Pointer to the grid constaining the results.
     */
    OBGridData *GetGrid(double step, double padding, const char *type, double pchg);
    /** 
     * @brief Get the pointer to the gradients.
     */
    //double* GetGradientPtr(); 
    //std::vector<Eigen::Vector3d>&  GetPositions() { return m_positions; } 
    //std::vector<Eigen::Vector3d>&  GetGradients() { return m_gradients; } 
 
    /////////////////////////////////////////////////////////////////////////
    // Interacting groups                                                  //
    /////////////////////////////////////////////////////////////////////////
      
    //! \name Methods for specifying interaction groups
    //@{
    /** 
     * @brief Enable intra-molecular interactions for group (bonds, angles, 
     * strbnd, torsions, oop). This function should be called before Setup().
     *
     * @param group OBBitVec with bits set for the indexes of the atoms which 
     * make up the group.
     */
    void AddIntraGroup(OBBitVec &group);
    /** 
     * @brief Enable inter-molecular interactions for group (non-bonded: vdw & 
     * ele). This function should be called before Setup().
     *
     * @param group OBBitVec with bits set for the indexes of the atoms which 
     * make up the group.
     */
    void AddInterGroup(OBBitVec &group);
    /** 
     * @brief Enable inter-molecular interactions between group1 and group2 
     * (non-bonded: vdw & ele). Note that this function doesn't enable bonded 
     * interactions in either group. Non-bonded interactions in the groups 
     * itself are also not enabled. This function should be called before Setup().
     *
     * @param group1 OBBitVec with bits set for the indexes of the atoms which 
     * make up the first group.
     * @param group2 OBBitVec with bits set for the indexes of the atoms which 
     * make up the second group.
     */
    void AddInterGroups(OBBitVec &group1, OBBitVec &group2);
    /** 
     * @brief Clear all previously specified groups.
     */
    void ClearGroups(); 
    /** 
     * @return True if there are groups.
     */ 
    bool HasGroups(); 
    //@}
 
    /////////////////////////////////////////////////////////////////////////
    // Cut-off                                                             //
    /////////////////////////////////////////////////////////////////////////
      
    //! \name Methods for Cut-off distances
    //@{
    /** 
     * @brief Enable or disable Cut-offs. Cut-offs are disabled by default.
     *  
     * @param enable Enable when true, disable when false.
     */
    void EnableCutOff(bool enable);
    /** 
     * @return True if Cut-off distances are used.
     */
    bool IsCutOffEnabled() { return m_cutoff; }
    /*! 
     * @brief Set the VDW cut-off distance to r. Note that this does not 
     * enable cut-off distances.
     * 
     * @param r The VDW cut-off distance to be used in A.
     */
    void SetVDWCutOff(double r);
    /** 
     * @brief Get the VDW cut-off distance.
     * 
     * @return The VDW cut-off distance in A.
     */
    double GetVDWCutOff();
    /** 
     * @brief Set the Electrostatic cut-off distance to r. Note that this 
     * does not enable cut-off distances.
     *
     * @param r The electrostatic cut-off distance to be used in A.
     */
    void SetElectrostaticCutOff(double r);
    /** 
     * @brief Get the Electrostatic cut-off distance.
     * 
     * @return The electrostatic cut-off distance in A.
     */
    double GetElectrostaticCutOff();
    /** 
     * @brief Set the frequency by which non-bonded pairs are updated. Values 
     * from 10 to 20 are recommended. Too low will decrease performance, too 
     * high will cause non-bonded interactions within cut-off not to be 
     * calculated.
     *
     * @param f The pair list update frequency.
     */ 
    void SetUpdateFrequency(int f);
    /** 
     * @brief Get the frequency by which non-bonded pairs are updated.
     *
     * @return The pair list update frequency.
     */ 
    int GetUpdateFrequency();
    /** 
     * @brief Set the bits in _vdwpairs and _elepairs to 1 for interactions that 
     * are within cut-off distance. This function is called in minimizing
     * algorithms such as SteepestDescent and ConjugateGradients. 
     */ 
    void UpdatePairsSimple();
    //void UpdatePairsGroup(); TODO
    /** 
     * @brief Get the number of non-bonded pairs in _mol.
     *
     * @return The number of pairs currently enabled (within cut-off distance)
     */ 
    unsigned int GetNumPairs();
    /** 
     * @brief Set bits in range 0..._numpairs-1 to 1. Using this means there will
     * be no cut-off. (not-working: see code for more information.
     */ 
    void EnableAllPairs()
    {
      // TODO: OBBitVec doesn't seem to be allocating it's memory correctly
      //_vdwpairs.SetRangeOn(0, _numpairs-1);
      //_elepairs.SetRangeOn(0, _numpairs-1);
    }
    //@}
 
    /////////////////////////////////////////////////////////////////////////
    // Energy Evaluation                                                   //
    /////////////////////////////////////////////////////////////////////////
      
    double Eval(bool gradients = true) 
    {
      return Energy(gradients);
    }

    //! \name Methods for energy evaluation
    //@{
    /** 
     * @param gradients Set to true when the gradients need to be calculated 
     * (needs to be done before calling GetGradient()).
     * 
     * @return Total energy.
     * 
     * @par Output to log:
     *    OBFF_LOGLVL_NONE:   none \n
     *    OBFF_LOGLVL_LOW:    none \n
     *    OBFF_LOGLVL_MEDIUM: energy for indivudual energy terms \n
     *    OBFF_LOGLVL_HIGH:   energy for individual energy interactions \n
     */
    virtual double Energy(bool gradients = true) { return 0.0f; }
    /** 
     * @param gradients Set to true when the gradients need to be calculated 
     * (needs to be done before calling GetGradient()).
     * 
     * @return Bond stretching energy.
     * 
     * @par Output to log:
     *    see Energy()
     */
    virtual double E_Bond(bool gradients = true) { return 0.0f; }
    /** 
     * @param gradients Set to true when the gradients need to be calculated 
     *  (needs to be done before calling GetGradient()).
     *  
     * @return Angle bending energy.
     *  
     * @par Output to log:
     *    see Energy()
     */
    virtual double E_Angle(bool gradients = true) { return 0.0f; }
    /** 
     * @param gradients Set to true when the gradients need to be calculated 
     * (needs to be done before calling GetGradient()).
     * 
     * @return Stretch bending energy.
     * 
     * @par Output to log:
     *    see Energy()
     */ 
    virtual double E_StrBnd(bool gradients = true) { return 0.0f; }
    /** 
     * @param gradients Set to true when the gradients need to be calculated 
     * (needs to be done before calling GetGradient()).
     * 
     * @return Torsional energy.
     * 
     * @par Output to log:
     *	  see Energy()
     */ 
    virtual double E_Torsion(bool gradients = true) { return 0.0f; }
    /** 
     * @param gradients Set to true when the gradients need to be calculated 
     * (needs to be done before calling GetGradient()).
     * 
     * @return Out-Of-Plane bending energy.
     * 
     * @par Output to log:
     *	  see Energy()
     */ 
    virtual double E_OOP(bool gradients = true) { return 0.0f; }
    /** 
     * @param gradients Set to true when the gradients need to be calculated 
     * (needs to be done before calling GetGradient()).
     * 
     * @return Van der Waals energy.
     * 
     * @par Output to log:
     *	  see Energy()
     */ 
    virtual double E_VDW(bool gradients = true) { return 0.0f; }
    /** 
     * @param gradients Set to true when the gradients need to be calculated 
     * (needs to be done before calling GetGradient()).
     * 
     * @return Electrostatic energy.
     * 
     * @par Output to log:
     *	  see Energy()
     */ 
    virtual double E_Electrostatic(bool gradients = true) { return 0.0f; }
    //@}
     
    /////////////////////////////////////////////////////////////////////////
    // Logging                                                             //
    /////////////////////////////////////////////////////////////////////////
      
    //! \name Methods for logging
    //@{
    /** 
     * @brief Print the atom types to the log.
     */
    void PrintTypes();
    /** 
     * @brief Print the formal charges to the log (atom.GetPartialCharge(), 
     *  MMFF94 FC's are not always int).
     */ 
    void PrintFormalCharges();
    /** 
     * @brief Print the partial charges to the log.
     */ 
    void PrintPartialCharges();
    /** 
     * @brief Print the velocities to the log.
     */ 
    void PrintVelocities();
    /** 
     * @brief Set the stream for logging (can also be &cout for logging to screen).
     * 
     * @param pos Stream (when pos is 0, std::cout wil be used).
     * 
     * @return True if succesfull.
     */
    bool SetLogFile(std::ostream *pos);
    /** 
     * @brief Set the log level (OBFF_LOGLVL_NONE, OBFF_LOGLVL_LOW, OBFF_LOGLVL_MEDIUM, OBFF_LOGLVL_HIGH).
     * Inline if statements for logging are available: 
     * @code
     * #define IF_OBFF_LOGLVL_LOW    if(GetLogLevel() >= OBFF_LOGLVL_LOW)
     * #define IF_OBFF_LOGLVL_MEDIUM if(GetLogLevel() >= OBFF_LOGLVL_MEDIUM)
     * #define IF_OBFF_LOGLVL_HIGH   if(GetLogLevel() >= OBFF_LOGLVL_HIGH)
     * @endcode
     *
     * example:
     * @code
     * SetLogLevel(OBFF_LOGLVL_MEDIUM);
     * IF_OBFF_LOGLVL_HIGH {
     *   OBFFLog("this text will NOT be logged...\n");
     * }
     *
     * IF_OBFF_LOGLVL_LOW {
     *   OBFFLog"this text will be logged...\n");
     * }
     *
     * IF_OBFF_LOGLVL_MEDIUM {
     *   OBFFLog("this text will also be logged...\n");
     * }
     * @endcode
     */
    bool SetLogLevel(int level);
    /** 
     * @return The log level.
     */ 
    int GetLogLevel() { return m_loglvl; }
    /** 
     * @brief Print msg to the logfile.
     * 
     * @param msg The message to print.
     */
    void OBFFLog(std::string msg);
    /** 
     * @brief Print msg to the logfile.
     * 
     * @param msg The message to print.
     */
    void OBFFLog(const char *msg);
    //@}
     
    /////////////////////////////////////////////////////////////////////////
    // Structure Generation                                                //
    /////////////////////////////////////////////////////////////////////////
      
    //! \name Methods for structure generation
    //@{
    /**
     * @brief Generate coordinates for the molecule (distance geometry). (OB 3.0).
     */
    void DistanceGeometry();
    /** 
     * @brief Generate conformers for the molecule (systematicaly rotating torsions).
     *  
     * The initial starting structure here is important, this structure should be
     * minimized for the best results. SystematicRotorSearch works by rotating around
     * the rotatable bond in a molecule (see OBRotamerList class). This rotating generates 
     * multiple conformers. The energy for all these conformers is then evaluated and the 
     * lowest energy conformer is selected.
     *
     * @param geomSteps The number of steps to take during geometry optimization.
     *        
     * @par Output to log:
     *  This function should only be called with the log level set to OBFF_LOGLVL_NONE or OBFF_LOGLVL_LOW. Otherwise
     *	too much information about the energy calculations needed for this function will interfere with the output for
     *	this function. \n\n
     *  OBFF_LOGLVL_NONE:   None. \n
     *  OBFF_LOGLVL_LOW:    Number of rotatable bonds, energies for the conformers, which one is the lowest, ... \n
     *  OBFF_LOGLVL_MEDIUM: See note above. \n
     *  OBFF_LOGLVL_HIGH:   See note above. \n 
     */
    void SystematicRotorSearch(unsigned int geomSteps = 2500);
    /** 
     * @brief Generate conformers for the molecule by systematicaly rotating torsions. To be used in combination with 
     * SystematicRotorSearchNexConformer().
     *
     * example:
     * @code
     * // pFF is a pointer to a OBForceField class 
     * pFF->SystematicRotorSearchInitialize(300);
     * while (pFF->SystematicRotorSearchNextConformer(300)) {
     *   // do some updating in your program (show last generated conformer, ...)
     * }
     * @endcode
     * 
     * If you don't need any updating in your program, SystematicRotorSearch() is recommended.
     *
     * @param geomSteps The number of steps to take during geometry optimization.
     * 
     * @return The number of conformers.
     */
    int SystematicRotorSearchInitialize(unsigned int geomSteps = 2500);
    /** 
     * @brief Evaluate the next conformer. 
     * 
     * @param geomSteps The number of steps to take during geometry optimization.
     * 
     * @return True if there are more conformers.
     */ 
    bool SystematicRotorSearchNextConformer(unsigned int geomSteps = 2500);
    /** 
     * @brief Generate conformers for the molecule (randomly rotating torsions).
     *  
     * The initial starting structure here is important, this structure should be
     * minimized for the best results. RandomRotorSearch works by randomly rotating around
     * the rotatable bonds in a molecule (see OBRotamerList class). This rotating generates 
     * multiple conformers. The energy for all these conformers is then evaluated and the 
     * lowest energy conformer is selected.
     *
     * @param conformers The number of random conformers to consider during the search.
     * @param geomSteps The number of steps to take during geometry optimization for each conformer.
     *        
     * @par Output to log:
     *  This function should only be called with the log level set to OBFF_LOGLVL_NONE or OBFF_LOGLVL_LOW. Otherwise
     *	too much information about the energy calculations needed for this function will interfere with the output for
     *	this function. \n\n
     *  OBFF_LOGLVL_NONE:   None. \n
     *  OBFF_LOGLVL_LOW:    Number of rotatable bonds, energies for the conformers, which one is the lowest, ... \n
     *  OBFF_LOGLVL_MEDIUM: See note above. \n
     *  OBFF_LOGLVL_HIGH:   See note above. \n 
     */
    void RandomRotorSearch(unsigned int conformers, unsigned int geomSteps = 2500);
    /** 
     * @brief Generate conformers for the molecule by randomly rotating torsions. To be used in combination with 
     * RandomRotorSearchNexConformer().
     *
     * example:
     * @code
     * // pFF is a pointer to a OBForceField class 
     * pFF->RandomRotorSearchInitialize(300);
     * while (pFF->RandomRotorSearchNextConformer(300)) {
     *   // do some updating in your program (show last generated conformer, ...)
     * }
     * @endcode
     * 
     * If you don't need any updating in your program, RandomRotorSearch() is recommended.
     *
     * @param conformers The number of random conformers to consider during the search
     * @param geomSteps The number of steps to take during geometry optimization
     */
    void RandomRotorSearchInitialize(unsigned int conformers, unsigned int geomSteps = 2500);
    /** 
     * @brief Evaluate the next conformer. 
     * 
     * @param geomSteps The number of steps to take during geometry optimization.
     * 
     * @return True if there are more conformers.
     */ 
    bool RandomRotorSearchNextConformer(unsigned int geomSteps = 2500);
    /** 
     * @brief Generate conformers for the molecule (randomly rotating torsions).
     *  
     * The initial starting structure here is important, this structure should be
     * minimized for the best results. WeightedRotorSearch works by randomly rotating around
     * the rotatable bonds in a molecule (see OBRotamerList class). Unlike RandomRotorSearch()
     * the random choice of torsions is reweighted based on the energy of the generated conformer.
     * Over time, the generated conformers for each step should become increasingly better.
     * The lowest energy conformer is selected.
     *
     * @param conformers The number of random conformers to consider during the search.
     * @param geomSteps The number of steps to take during geometry optimization for each conformer.
     *        
     * @par Output to log:
     *  This function should only be called with the log level set to OBFF_LOGLVL_NONE or OBFF_LOGLVL_LOW. Otherwise
     *	too much information about the energy calculations needed for this function will interfere with the output for
     *	this function. \n\n
     *  OBFF_LOGLVL_NONE:   None. \n
     *  OBFF_LOGLVL_LOW:    Number of rotatable bonds, energies for the conformers, which one is the lowest, ... \n
     *  OBFF_LOGLVL_MEDIUM: See note above. \n
     *  OBFF_LOGLVL_HIGH:   See note above. \n 
     */
    void WeightedRotorSearch(unsigned int conformers, unsigned int geomSteps);

    /////////////////////////////////////////////////////////////////////////
    // Molecular Dynamics                                                  //
    /////////////////////////////////////////////////////////////////////////
      
    //! \name Methods for molecular dynamics
    //@{
    /** 
     * @brief Generate starting velocities with a Maxwellian distribution.
     */
    void GenerateVelocities();
    /** 
     * @brief Correct the velocities so that the following is true:
     *  
     * @code
     *        3N
     *       ----
     *  0.5  \    m_i * v_i^2 = 0.5 * Ndf * kB * T = E_kin
     *       /
     *       ----
     *       i=1
     *  
     *  E_kin : kinetic energy
     *  m_i : mass of atom i
     *  v_i : velocity of atom i
     *  Ndf : number of degrees of freedom (3 * number of atoms)
     *  kB : Boltzmann's constant
     *  T : temperature
     * @endcode
     *  
     */
    void CorrectVelocities();
    /** 
     * @brief Take n steps at temperature T. If no velocities are set, they will be generated.
     *
     * example:
     * @code
     *  // pFF is a pointer to a OBForceField class 
     *  while (pFF->MolecularDynamicsTakeNSteps(5, 300)) {
     *    // do some updating in your program (redraw structure, ...)
     *  }
     * @endcode
     *
     * @param n The number of steps to take.
     * @param T Absolute temperature in Kelvin.
     * @param timestep The time step in picoseconds. (10e-12 s)
     */
    void MolecularDynamicsTakeNSteps(int n, double T, double timestep = 0.001);
    //@}

    /////////////////////////////////////////////////////////////////////////
    // Constraints                                                         //
    /////////////////////////////////////////////////////////////////////////
      
    //! \name Methods for constraints
    //@{
    /** 
     * @brief Get the current constraints.
     * 
     * @return The current constrains stored in the force field.
     */ 
    OBFFConstraints& GetConstraints();
    /** 
     * @brief Set the constraints.
     * 
     * @param constraints The new constraints to be used.
     */
    void SetConstraints(OBFFConstraints& constraints);
    /** 
     * @brief Fix the atom position until UnsetFixAtom() is called. This function 
     * can be used in programs that allow the user to interact with a molecule
     * that is being minimized without having to check if the atom is already 
     * fixed in the constraints set by Setup() or SetConstraints(). Using this 
     * makes sure the selected atom follows the mouse cursur.
     *  
     * @param index The index for the atom to fix.
     */
    void SetFixAtom(int index);
    /** 
     * @brief Undo last SetFixAtom. This function will not remove the fix atom 
     * constraint for this atom if set by Setup() or SetConstraints().
     */
    void UnsetFixAtom();
    /** 
     * @brief Ignore the atom until UnsetIgnoreAtom() is called. This function 
     * can be used in programs that allow the user to interact with a molecule
     * that is being minimized without having to check if the atom is already 
     * ignored in the constraints set by Setup() or SetConstraints(). Using this 
     * makes sure, in drawing mode, you can close rings without your newly 
     * created puching the other atoms away.
     * 
     * @param index The index for the atom to ignore.
     */
    void SetIgnoreAtom(int index);
    /** 
     * @brief Undo last SetIgnoreAtom. This function will not remove the ignore atom 
     * constraint for this atom if set by Setup() or SetConstraints().
     */
    void UnsetIgnoreAtom();
    /** 
     * @brief 
     * 
     * @param a
     * @param b
     * 
     * @return 
     */
    static bool IgnoreCalculation(unsigned int a, unsigned int b);
    /** 
     * @brief 
     * 
     * @param a
     * @param b
     * @param c
     * 
     * @return 
     */
    static bool IgnoreCalculation(unsigned int a, unsigned int b, unsigned int c);
    /** 
     * @brief 
     * 
     * @param a
     * @param b
     * @param c
     * @param d
     * 
     * @return 
     */
    static bool IgnoreCalculation(unsigned int a, unsigned int b, unsigned int c, unsigned int d);
    //@}

 
    /////////////////////////////////////////////////////////////////////////
    // Validation                                                          //
    /////////////////////////////////////////////////////////////////////////
      
    //! \name Methods for forcefield validation
    //@{
    //! (debugging)
    bool DetectExplosion();
    /** 
     * @brief Validate the analytical gradients by comparing them to 
     * numerical ones. This function has to be implemented force field 
     * specific. (debugging)
     */
    virtual bool ValidateGradients() { return false; }
    /** 
     * @brief Calculate the error of the analytical gradient (debugging)
     * 
     * @return  error = fabs(numgrad - anagrad) / anagrad * 100% 
     */
    Eigen::Vector3d ValidateGradientError(Eigen::Vector3d &numgrad, Eigen::Vector3d &anagrad);
    //@}
     
  }; // class OBForceField

}// namespace OpenBabel

#endif   // OB_FORCEFIELD_H

//! @file forcefield.h
//! @brief Handle forcefields

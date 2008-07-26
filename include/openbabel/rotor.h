/**********************************************************************
rotor.h - Rotate torsional according to rotor rules.
 
Copyright (C) 1998-2000 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2005 by Geoffrey R. Hutchison
Some portions Copyright (C) 2008 by Tim Vandermeersch
 
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

#ifndef OB_ROTOR_H
#define OB_ROTOR_H

#include <openbabel/parsmart.h>
#include <openbabel/typer.h>

namespace OpenBabel
{

#ifndef SQUARE
#define SQUARE(x) ((x)*(x))
#endif

  /** @class OBRotorRule rotor.h <openbabel/rotor.h>
   *  @brief A rule for torsional conformer searching, defined by a SMARTS pattern
   * 
   *  Rules define a SMARTS pattern to match a set of 4 reference atoms
   *  defining the dihedral angle. The rule can either define a set of possible
   *  dihedral angles in degrees and/or a "delta" (i.e., all multiples of delta will
   *  be considered)
   */
  class OBAPI OBRotorRule
  {
    private:
      //! Reference atoms specifying the dihedral angle (as integers), numbered from 1 inside the SMARTS pattern
      unsigned int        m_ref[4]; 
      std::string         m_smarts;    //!< Text of the SMARTS pattern
      OBSmartsPattern*    m_smartspat; //!< The SMARTS pattern for the rotation rule
      std::vector<double> m_vals;      //!< At least one torsion angle (in radians) to evaluate
    public:
      /**
       * Constructor with torsion values.
       * @param buffer The smarts pattern for this rotor rule
       * @param ref Reference atoms specifying the dihedral angle (as integers), numbered from 1 inside the SMARTS pattern
       * @param vals The torsion values in radians.
       */
      OBRotorRule(char *buffer, unsigned int ref[4], std::vector<double> &vals) : m_smarts(buffer), m_vals(vals)
      {
        m_smartspat = new OBSmartsPattern;
        m_smartspat->Init(buffer);
        for (int i = 0; i < 4; ++i)
          m_ref[i] = ref[i];
      }
      /**
       * Constructor with torsion increment.
       * @param buffer The smarts pattern for this rotor rule
       * @param ref Reference atoms specifying the dihedral angle (as integers), numbered from 1 inside the SMARTS pattern
       * @param delta The torsion increment in degrees.
       */
      OBRotorRule(char *buffer, unsigned int ref[4], double delta) : m_smarts(buffer)
      {
        std::vector<double> vals;

        // create torsion values from the increment
        double tor = 0.0;
        while (tor < 360.0) {
          vals.push_back(tor);
          tor += delta;
        }

        OBRotorRule(buffer, ref, vals);
      }
      /**
       * Destructor.
       */
      ~OBRotorRule()
      {
        if (m_smartspat) {
          delete m_smartspat;
          m_smartspat = NULL;
        }
      }
      /**
       * @return whether this rotor rule is valid (i.e., is the SMARTS pattern valid)
       */
      bool IsValid() const { return m_smartspat->IsValid(); }
      /**
       * @return a copy of the reference atom indexes inside the SMARTS pattern
       */
      void GetReferenceAtoms(unsigned int ref[4]) const
      { 
        for (int i = 0; i < 4; ++i)
          ref[i] = m_ref[i];
      }
      /** 
       * @return a reference to the dihedral angles to evaluate (in radians)
       */
      std::vector<double>& GetTorsionVals() { return m_vals; }
      /**
       * @return the text of the SMARTS pattern for this rule
       */
      std::string& GetSmartsString() { return m_smarts; }
      /** 
       * @return the exact OBSmartsPattern object for this rule
       */
      OBSmartsPattern* GetSmartsPattern() const { return m_smartspat; }
  };

  /** @class OBRotorRules rotor.h <openbabel/rotor.h>
   *  @brief Database of default hybridization torsional rules and SMARTS-defined OBRotorRule objects
   * 
   *  Use to automatically evaluate potentially rotatable bonds to generate
   *  lists of dihedral angles to consider.
   *  e.g., rotamer/conformer energy calculations
   */
  class OBAPI OBRotorRules : public OBGlobalDataBase
  {
    private:
      bool                       m_quiet;  //!< Control debugging output from GetRotorIncrements()
      std::vector<OBRotorRule*>  m_vr;     //!< Database of specific OBRotorRules defined by SMARTS patterns
      //! Default dihedral angles to check for generic sp3 - sp3 hybridized rotatable bonds (in radians)
      std::vector<double>        m_sp3sp3; 
      //! Default dihedral angles to check for generic sp3 - sp2 hybridized rotatable bonds (in radians)
      std::vector<double>        m_sp3sp2; 
      //! Default dihedral angles to check for generic sp2 - sp2 hybridized rotatable bonds (in radians)
      std::vector<double>        m_sp2sp2; 
    public:
      /**
       * Constructor.
       */
      OBRotorRules();
      /**
       * Destructor.
       */
      ~OBRotorRules();
      /**
       * Parse a single line from the torsion database file.
       */ 
      void ParseLine(const char*);
      /**
       * @return the number of rotor rules
       */
      //unsigned int Size() const { return m_vr.size();}
      unsigned int GetSize() { return m_vr.size();}
      /** 
       * Set the filename to be used for the database. Default = torlib.txt
       */
      void SetFilename(std::string &s) { _filename = s; }
      /** 
       * Determine the torsional angles to evaluate based on the database
       * @param mol molecule to evaluate
       * @param bond rotatable bond to evaluate
       * @param refs set to be the atom indexes (in mol) of the dihedral angle
       * @param vals set to be the list of angles to evaluate (in radians)
       */
      void GetRotorIncrements(OBMol& mol, OBBond* bond, unsigned int refs[4], std::vector<double> &vals);
      /**
       *  Turn off debugging output from GetRotorIncrements()
       */
      void Quiet() { m_quiet = true; }
  };

  //! \class OBRotor rotor.h <openbabel/rotor.h>
  //! \brief A single rotatable OBBond as part of rotamer searching
  class OBAPI OBRotor
  {
    private:
      int                 m_idx;       //!< the index for this OBRotor
      unsigned int        m_ref[4];    //!< the four torsion atom indexes
      std::vector<int>    m_cidx;      //!< the four torsion atom coordinate indexes
      std::vector<int>    m_rotatoms;  //<! the atoms that should be rotated
      int                 m_numcoords; //!< the number of coordinates
      OBBond             *m_bond;      //!< the bc bond in torsion abcd
      std::vector<double> m_torsions;  //!< torsion values
    public:
      /**
       * Constructor.
       */
      OBRotor();
      /**
       * Destructor.
       */
      ~OBRotor();
      
      //! @name OBRotor modification methods
      //@{
      /**
       * Set the index for this OBRotor.
       */
      void SetIdx(int idx) { m_idx = idx; }
      /**
       * Set the number of coordinates.
       */
      void SetNumCoords(int nc) { m_numcoords = nc; }
      /**
       * Set the bond for this OBRotor.
       */
      void SetBond(OBBond *bond) { m_bond = bond; }
      /**
       * Set the indexes for the four torsion atoms abcd. These indexes will
       * be used when accessing the coordinates in the *c pointers.
       */
      void SetDihedralAtoms(unsigned int ref[4]);
      /**
       * Set the torsion values for this rotor.
       * @param torsions The torsion values in radians.
       */
      void SetTorsionValues(std::vector<double> &torsions) { m_torsions = torsions; }
      /**
       * Specify the atoms to which the rotation will be applied.
       */
      void SetRotAtoms(std::vector<int> &atoms) { m_rotatoms = atoms; }
      /**
       * Remove symmetry. For 2-fold symmetry, all torsion values >180 degrees will 
       * be removed. For 3-fold symmetry, all torsion values >120 degrees will be 
       * removed.
       * @param fold The symmetry (2 or 3).
       */
      void RemoveSymTorsionValues(int fold);
      //@}

      //! \name OBRotor data request methods
      //@{
      /**
       * Get the index for this OBRotor
       */
      int GetIdx() const { return m_idx; }
      /**
       * @return the atoms to which the rotation will be applied.
       */
      std::vector<int>& GetRotAtoms() { return m_rotatoms; }
      /**
       * @return the bond for this OBRotor.
       */
      OBBond* GetBond() const { return m_bond; }
      /**
       * @return the atom indexes for the four torsion atoms abcd.
       */
      void GetDihedralAtoms(unsigned int ref[4]) const
      {
        for (int i=0;i<4;++i)
          ref[i] = m_ref[i];
      }
      /**
       * @return the torsion values for this rotor.
       */ 
      std::vector<double>& GetTorsionValues() { return m_torsions; }
      /**
       * @return the number of torsion values.
       */
      int Size() const { return m_torsions.size(); }
      //@}

      //! \name Rotation methods
      //@{
      /**
       * Rotate the atoms for this rotor so the torsion angle is equal to setang.
       * This function does not take any rotor rules into account.
       * @param c Pointer to all coordinates.
       * @param setang The torsion angle. (in radians)
       */
      void SetToAngle(double *c, double setang);
      //@}
    
      //! @name Iteration methods
      //@{
      std::vector<double>::iterator BeginTorIncrement() { return m_torsions.begin(); }
      std::vector<double>::iterator EndTorIncrement() { return m_torsions.end(); }
      //@}
  };

  //! A standard iterator over a vector of rotors
  typedef std::vector<OBRotor*>::iterator OBRotorIterator;

  /** @class OBRotorList
   *  @brief Given an OBMol, set up a list of possibly rotatable torsions.
   */
  class OBAPI OBRotorList
  {
    bool _quiet;                    //!< Control debugging output
    OBBitVec _fix;                  //!< Bit vector of fixed (i.e., invariant) atoms
    OBRotorRules _rr;               //!< Database of rotatable bonds and dihedral angles to test
    std::vector<int> _dffv;         //!< Distance from fixed
    std::vector<OBRotor*> m_rotors;   //!< List of individual OBRotor torsions
    //!
    std::vector<std::pair<OBSmartsPattern*,std::pair<int,int> > > _vsym2;
    //! 
    std::vector<std::pair<OBSmartsPattern*,std::pair<int,int> > > _vsym3;
  public:
    /**
     * Constructor.
     */
    OBRotorList();
    /**
     * Destructor.
     */
    ~OBRotorList();

    /** 
     * Clear the internal list of rotors and reset.
     */
    void   Clear();

    /**
     *  @return The number of rotors in this list.
     */
    unsigned int Size() const { return m_rotors.size(); }
    //! Intialize the private OBRotorRules database from a specific file
    void Init(std::string &fname)
    {
      _rr.SetFilename(fname);
      _rr.Init();
    }
    /**
     * Turn off debugging output
     */
    void SetQuiet() { _quiet=true; _rr.Quiet(); }

    //! Set the list of fixed (invariant) atoms to the supplied OBBitVec
    void   SetFixAtoms(OBBitVec &fix) { _fix = fix;        }

    //! \brief Return whether this bond is fixed and thus not rotatable
    //! \return true if the bond and at least one neighboring bond has fixed atoms
    bool   IsFixedBond(OBBond*);
    //! \return whether this rotor list has any fixed (invariant) atoms
    bool   HasFixedAtoms()
    {
      return(!_fix.Empty());
    }

    //! \brief Set the atoms to rotate from the dihedral atoms for each rotor
    //! Insures the fixed atoms are respected, but otherwise functions like
    //! SetRotAtoms()
    void   SetRotAtomsByFix(OBMol&);

    //! \brief Set the atoms to rotate from the dihedral atoms for each rotor
    //! Uses OBRotor->GetDihedralAtoms() to call OBRotor->SetRotAtoms()
    //! and standarizes the dihedral angles via OBRotor->SetDihedralAtoms()
    //! \return True
    bool   SetRotAtoms(OBMol&);

    //! Setup this rotor list for the supplied molecule
    //! Calls FindRotors() and AssignTorVals()
    //! \return True if rotatable bonds were found
    bool   Setup(OBMol &);
    //! Find all potentially rotatable bonds in the molecule
    //! Uses OBBond::IsRotor() for initial evaluation
    //! \return True
    bool   FindRotors(OBMol &);
    //! Using the OBRotorRules database, set the torsion values (and delta)
    //! to be evaluated and tested
    bool   AssignTorVals(OBMol &);

    //! Rotates each bond to zero and 180 degrees and tests
    //! if the 2 conformers are duplicates.  if so - the symmetric torsion
    //! values are removed from consideration during a search
    void   RemoveSymVals(OBMol&);

    
    
    //! \name Iterator methods
    //@{
    OBRotor *BeginRotor(OBRotorIterator &i)
    { i = m_rotors.begin(); return ((i == m_rotors.end()) ? 0: *i); }
    OBRotor *NextRotor(OBRotorIterator &i)
    { ++i; return ((i == m_rotors.end()) ? 0: *i); }
    OBRotorIterator BeginRotors() { return m_rotors.begin(); }
    OBRotorIterator EndRotors()   { return m_rotors.end();   }
    //@}
  };


} // end namespace OpenBabel

#endif // OB_ROTOR_H

//! \file rotor.h
//! \brief Rotate torsional according to rotor rules.

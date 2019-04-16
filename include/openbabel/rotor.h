/**********************************************************************
rotor.h - Rotate torsional according to rotor rules.

Copyright (C) 1998-2000 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2005 by Geoffrey R. Hutchison

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

#ifndef OB_ROTOR_H
#define OB_ROTOR_H

#include <openbabel/parsmart.h>
#include <openbabel/typer.h>
#include <openbabel/bitvec.h>

#ifdef UNUSED
#elif (__GNUC__ == 4)
# define UNUSED(x) UNUSED_ ## x __attribute__((unused))
#elif defined(__LCLINT__)
# define UNUSED(x) /*@unused@*/ x
#else
# define UNUSED(x) x
#endif

namespace OpenBabel
{
  class OBRing;

#ifndef SQUARE
#define SQUARE(x) ((x)*(x))
#endif

  //! \class OBRotorRule rotor.h <openbabel/rotor.h>
  //! \brief A rule for torsional conformer searching, defined by a SMARTS pattern
  //!
  //! Rules define a SMARTS pattern to match and a set of 4 reference atoms
  //! defining the dihedral angle. The rule can either define a set of possible
  //! dihedral angles in degrees and/or a "delta" (i.e., all multiples of delta will
  //! be considered)
  class OBAPI OBRotorRule
  {
    int                 _ref[4]; //!< Reference atoms specifying the dihedral angle (as integers), numbered from 1 inside the SMARTS pattern
    double              _delta;  //!< (optional) the resolution of a dihedral step in degrees
    std::string         _s;      //!< Text of the SMARTS pattern
    OBSmartsPattern*    _sp;     //!< The SMARTS pattern for the rotation rule
    std::vector<double> _vals;   //!< At least one torsion angle (in radians) to evaluate
  public:

  OBRotorRule(char *buffer,int ref[4],std::vector<double> &vals,double d):
    _delta(d), _s(buffer), _vals(vals)
    {
      _sp = new OBSmartsPattern;
      _sp->Init(buffer);
      memcpy(_ref,ref,sizeof(int)*4);
    }

    ~OBRotorRule()
      {
        if (_sp)
          {
            delete _sp;
            _sp = NULL;
          }
      }

    //! \return whether this rotor rule is valid (i.e., is the SMARTS pattern valid)
    bool    IsValid()    {        return(_sp->IsValid());       }
    //! \return a copy of the reference atom indexes inside the SMARTS pattern
    //!
    //!  These should be freed after use.
    void    GetReferenceAtoms(int ref[4]) { memcpy(ref,_ref,sizeof(int)*4); }
    //! Set the resolution (delta) of a torsional step in degrees
    void    SetDelta(double d)    {       _delta = d;           }
    //! \return the resolution (delta) of a torsional step in degrees
    double  GetDelta()            {       return(_delta);       }
    //! \return a reference to the dihedral angles to evaluate (in radians)
    std::vector<double>   &GetTorsionVals()    { return(_vals); }
    //! \return the text of the SMARTS pattern for this rule
    std::string  &GetSmartsString(){      return(_s);           }
    //! \return the exact OBSmartsPattern object for this rule
    OBSmartsPattern *GetSmartsPattern() {  return(_sp);         }
  };

  //! \class OBRotorRules rotor.h <openbabel/rotor.h>
  //! \brief Database of default hybridization torsional rules and SMARTS-defined OBRotorRule objects
  //!
  //! Use to automatically evaluate potentially rotatable bonds to generate
  //! lists of dihedral angles to consider.
  //! e.g., rotamer/conformer energy calculations
  class OBAPI OBRotorRules : public OBGlobalDataBase
  {
    bool                       _quiet;  //!< Control debugging output from GetRotorIncrements()
    std::vector<OBRotorRule*>  _vr;     //!< Database of specific OBRotorRules defined by SMARTS patterns
    std::vector<double>        _sp3sp3; //!< Default dihedral angles to check for generic sp3 - sp3 hybridized rotatable bonds (in radians)
    std::vector<double>        _sp3sp2; //!< Default dihedral angles to check for generic sp3 - sp2 hybridized rotatable bonds (in radians)
    std::vector<double>        _sp2sp2; //!< Default dihedral angles to check for generic sp2 - sp2 hybridized rotatable bonds (in radians)
  public:
    OBRotorRules();
    ~OBRotorRules();

    void ParseLine(const char*);
    //! \return the number of rotor rules
    size_t GetSize()                 { return _vr.size();}

    //! Set the filename to be used for the database. Default = torlib.txt
    void SetFilename(std::string &s)       { _filename = s;    }

    //! Determine the torsional angles to evaluate based on the database
    //! \param mol molecule to evaluate
    //! \param bond rotatable bond to evaluate
    //! \param refs set to be the atom indexes (in mol) of the dihedral angle
    //! \param vals set to be the list of angles to evaluate (in radians)
    //! \param delta potential dihedral angle steps (in degrees)
    void GetRotorIncrements(OBMol& mol,OBBond* bond,int refs[4],
                            std::vector<double> &vals,double &delta);
    //! Turn off debugging output from GetRotorIncrements()
    void Quiet()                           { _quiet=true;      }
  };

  /**
   * @class OBRotor rotor.h <openbabel/rotor.h>
   * @brief A single rotatable OBBond as part of rotamer searching
   */
  class OBAPI OBRotor
  {
    int _idx; //!< the index in an OBRotorList
    std::vector<int> _rotatoms; //!< the atoms to rotate
    double _imag, _refang; //!< inverse magnitude and reference angle (see Precompute())
    OBBond *_bond; //!< the bond associated with this rotor
    std::vector<int> _ref, _torsion; //!< indexes for atom coordinates (from 0, multiplied by 3)
    OBBitVec _fixedatoms,_fixedbonds, _evalatoms; //!< fixed atoms/bonds
    std::vector<double> _torsionAngles;  //!< torsion resolution
    std::vector<double> _invmag; //!< the inverse magnitudes (see Precalc)
    std::vector<std::vector<double> > _sn,_cs,_t; //!< the rotation matrix (see Precalc())
    std::vector<OBRing *> _rings; //!< the parent ring (if this is a rotor in a ring)
  public:
    /**
     * Constructor.
     */
    OBRotor();
    /**
     * Destructor.
     */
    ~OBRotor()
      {
      }

    ///@name Setup
    ///@{
    /**
     * Set the OBBond associated with this OBRotor.
     */
    void SetBond(OBBond *bond)
    {
      _bond = bond;
      SetRings();
    }
    /**
     * Set the rings associated with this bond (if it's a ring bond)
     * \since Version 2.4
     */
    void SetRings();
    /**
     * Set the index for this rotor. Used by OBRotorList
     */
    void SetIdx(int idx)
    {
      _idx = idx;
    }
    /**
     * Set the dihedral atoms.
     * @param ref The dihedral atom indexes. These indexes start from 1.
     */
    void SetDihedralAtoms(std::vector<int> &ref);
    /**
     * Set the dihedral atoms.
     * @param ref The dihedral atom indexes. These indexes start from 1.
     */
    void SetDihedralAtoms(int ref[4]);
    /**
     * Set the atom indexes that will be displaced when this rotor
     * changes torsion angle. These indexes start from 0 and are multiplied
     * by 3 for easy coordinate access.
     */
    void SetRotAtoms(std::vector<int> &atoms);
    /**
     * Set the possible torsion values or angles.
     */
    void SetTorsionValues(std::vector<double> &angles)
    {
      _torsionAngles = angles;
    }
    /**
     * Set the bonds that will be fixed.
     */
    void SetFixedBonds(OBBitVec &bv)
    {
      _fixedbonds = bv;
    }
    ///@}

    ///@name Performing rotations
    ///@{
    /**
     * Rotate the atoms in the specified @p coordinates to the specified angle.
     * @param coordinates The coordinates to rotate.
     * @param setang The new torsion angle in radians.
     */
    inline void SetToAngle(double *coordinates, double setang)
    {
      double /*dx,dy,dz,*/ sn,cs,t,ang,mag;
      // compute the angle to rotate (radians)
      ang = setang - CalcTorsion(coordinates);
      // if the angle to rotate is too small, we're done
      if (fabs(ang) < 1e-5)
        return;

      // compute the bond length
      mag = CalcBondLength(coordinates);
      // compute some rotation matrix elements
      sn = sin(ang);
      cs = cos(ang);
      t = 1 - cs;

      // perform rotation
      Set(coordinates, sn, cs, t, 1.0 / mag);
    }
    /**
     * Rotate the atoms in the specified @p coordinates. This function does not
     * require any precomputation and will compute all needed information when
     * needed.
     * @param coordinates The coordinates for the molecules as pointer to double.
     * @param next The index of the new rotor angle. This is an index for the
     * GetTorsionValues() list.
     * @param prev If specified, the torsion current torsion angle can be
     * looked up and does not have to be calculated again.
     */
    void SetRotor(double *coordinates, int next, int prev = -1);
    /**
     * Rotate the specified @p coordinates by using the specified rotation matrix.
     *
     */
    void Set(double *coordinates, double sine, double cosine, double translation, double invmag);
    /**
     * Precompute the reference angle and inverse bond length of this rotor for
     * a single conformer. This function should be used in combination with
     * Set(double *coordinates, int idx).
     * @param coordinates The coordinates to use in the computation.
     *
     * @code
     * OBMol mol;
     * ...
     *
     * unsigned int numCoords = mol.NumAtoms() * 3;
     * double *coords = mol.GetCoordinates();
     * OBRotor rotor;
     * rotor.SetBond(mol.GetBond(3));
     *
     * // set the possible torsion values
     * std::vector<double> angles;
     * angles.push_back(0.0);
     * angles.push_back(3.1415);
     * rotor.SetTorsionValues(angles);
     *
     * // precompute inverse bond length (i.e. the bond length of bond with index 3
     * // using the specified coordinates) and reference angle (i.e. torsion angle
     * //in coords)
     * rotor.Precompute(coords);
     *
     * // copy coordinates to coords_1
     * double *coords_1 = new double[numCoords];
     * for (unsigned int i = 0; i < numCoords; ++i)
     *   coords_1[i] = coords[i];
     * // rotate the atoms in coords_1 to angle with index 0 (i.e. 0.0 degrees)
     * // note: on input, the coordinates should be the same as the coordinates used
     * //       to precompute the inverse bond length and reference angle (in other
     * //       words, the inverse magnitude and reference angle in the specfied
     * //       coordinates should be the same as the one used for Precompute)
     * rotor.Set(coords_1, 0)
     *
     * // copy coordinates to coords_2
     * double *coords_2 = new double[numCoords];
     * for (unsigned int i = 0; i < numCoords; ++i)
     *   coords_2[i] = coords[i];
     * // rotate the atoms in coords_2 to angle with index 1 (i.e. 180.0 degrees)
     * rotor.Set(coords_2, 1)
     *
     * delete coords_1;
     * delete coords_2;
     * @endcode
     */
    void Precompute(double *coordinates);
    /**
     * Rotate the @p coordinates to set the torsion angle of this rotor to the angle
     * specified by the index @p idx. Make sure to call Precompute before calling
     * this function.
     * @param coordinates The coordinates to rotate.
     * @param idx The index of the torsion angle in the GetTorsionValues() list.
     */
    void Set(double *coordinates, int idx);
    /**
     * Precompute the inverse bond lengths, rotation matrices for all
     * specified conformers and all possible torsion values. This method is
     * used in combination with Set(double *coordinates, int conformer, int idx).
     * @param conformers The pointers to the conformer coordinates
     */
    void Precalc(std::vector<double*> &conformers);
    /**
     * Rotate the @p coordinates to set the torsion to the torsion value with the
     * specified @p index. The coordinates should be the same as the conformer used
     * for calling Precalc (i.e. conformers[conformer] == coordinates). Make sure
     * to call Precalc before calling this method.
     * @param coordinates The conformer coordinates.
     * @param conformer The conformer index in the conformer list given to Precalc().
     * @param idx The torsion value index in the GetTorsionValues() list.
     */
    void Set(double *coordinates, int conformer, int idx)
    {
      Set(coordinates, _sn[conformer][idx], _cs[conformer][idx], _t[conformer][idx], _invmag[conformer]);
    }
    ///@}


    ///@name Methods to retrieve information
    ///@{
    /**
     * Get the OBBond object associated with this OBRotor.
     */
    OBBond *GetBond()
    {
      return(_bond);
    }
    /**
     * Get the number of possible torsion angles for this OBRotor. This
     * is the length of the GetTorsionValues() list.
     */
    size_t Size()
    {
      return _torsionAngles.size();
    }
    /**
     * Get the index for this rotor (index in an OBRotorList).
     */
    int GetIdx() const
    {
      return _idx;
    }
    /**
     * Get the dihedral atom indexes. These indexes start from 1.
     */
    void GetDihedralAtoms(int ref[4])
    {
      for (int i = 0; i < 4; ++i)
        ref[i] = _ref[i];
    }
    /**
     * Get the dihedral atom indexes. These indexes start from 1.
     */
    std::vector<int> &GetDihedralAtoms()
      {
        return _ref;
      }
    /**
     * Get the atom indexes that will be displaced when this rotor changes
     * torsion angle. These indexes start from 1.
     */
    const std::vector<int>& GetRotAtoms() const
    {
      return _rotatoms;
    }
    /**
     * Get the possible torsion angles for this OBRotor.
     */
    const std::vector<double> &GetTorsionValues() const
    {
      return _torsionAngles;
    }
    /**
     * Get an OBBitVec objects with bits set for all bonds that are fixed.
     * Bonds are indexed from 0.
     */
    OBBitVec &GetFixedBonds()
      {
        return _fixedbonds;
      }
    /**
     * Calculate the torsion for this OBRotor using the specified coordinates.
     * @param coordinates The coordinates (e.g. OBMol::GetCoordinates()).
     * @return The torsion angle in radians.
     */
    double CalcTorsion(double *coordinates);
    /**
     * Calculate the bond length for this OBRotor using the specified coordinates.
     * @param coordinates The coordinates (e.g. OBMol::GetCoordinates()).
     */
    double CalcBondLength(double *coordinates);
    ///@}


    ///@name Iterator methods
    ///@{
    std::vector<double>::iterator BeginTorIncrement()
      {
        return _torsionAngles.begin();
      }
    std::vector<double>::iterator EndTorIncrement()
      {
        return _torsionAngles.end();
      }
    ///@}

    void RemoveSymTorsionValues(int);

    ///@name Deprecated
    ///@{
    /** @deprecated Has no effect. */
    void SetDelta(double UNUSED(d)) {}
    /** @deprecated Has no effect. */
    double GetDelta() { return 10.0; }
    /** @deprecated */
    OBBitVec &GetFixedAtoms() { return _fixedatoms; }
    /** @deprecated See SetFixedBonds */
    void SetFixedAtoms(OBBitVec &bv) { _fixedatoms = bv; }
    /** @deprecated */
    OBBitVec &GetEvalAtoms() { return _evalatoms; }
    /** @deprecated */
    void SetEvalAtoms(OBBitVec &bv) { _evalatoms = bv; }
    /** @deprecated */
    void* GetRotAtoms() { return &_rotatoms; }
    /** @deprecated Bad name, see GetTorsionValues() */
    std::vector<double> &GetResolution() { return _torsionAngles; }
    /** @deprecated */
    void SetNumCoords(int UNUSED(nc)) {}
    ///@}

  };


  //! A standard iterator over a vector of rotors
  typedef std::vector<OBRotor*>::iterator OBRotorIterator;

  /**
   * @class OBRotorList rotor.h <openbabel/rotor.h>
   * @brief Given an OBMol, set up a list of possibly rotatable torsions,
   */
  class OBAPI OBRotorList
  {
    bool _quiet;                    //!< Control debugging output
    bool _removesym;                //!< Control removal of symmetric rotations
    bool _ringRotors;               //!< Are there ring rotors
    OBBitVec _fixedatoms, _fixedbonds; //!< Bit vector of fixed (i.e., invariant) atoms
    OBRotorRules _rr;               //!< Database of rotatable bonds and dihedral angles to test
    std::vector<int> _dffv;         //!< Distance from fixed
    std::vector<OBRotor*> _rotor;   //!< List of individual OBRotor torsions
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
    void Clear();
    /**
     * @return the number of rotors in this list
     */
    size_t Size()
    {
      return _rotor.size();
    }
    /**
     * When no atoms/bonds are fixed or when bonds are fixed, this function will
     * return true if the bond is fixed. When using the deprecated fixed atoms,
     * this function will return true if the bond and at least one neighboring
     * bond has fixed atoms.
     */
    bool IsFixedBond(OBBond*);
    /**
     * @return True if this rotor list has any fixed bonds.
     */
    bool HasFixedBonds()
    {
      return !_fixedbonds.IsEmpty();
    }
    //! Rotates each bond to zero and 180 degrees and tests
    //! if the 2 conformers are duplicates.  if so - the symmetric torsion
    //! values are removed from consideration during a search
    void RemoveSymVals(OBMol&);

    /**
     * @return True if this rotor list has any ring bonds.
     * @since version 2.4
     */
    bool HasRingRotors()
    { return _ringRotors; }

    ///@name Setup
    /**
     * Setup this rotor list for the supplied molecule. This method calls
     * FindRotors(), SetEvalAtoms(), and AssignTorVals().
     * @param mol The molecule.
     * @param sampleRings Whether to sample ring conformers - default = false
     * @return True if rotatable bonds were found.
     */
    bool Setup(OBMol &mol, bool sampleRings = false);
    /**
     * Set the bonds that will be fixed.
     */
    void SetFixedBonds(OBBitVec &fix)
    {
      _fixedbonds = fix;
      _fixedatoms.Clear();
    }
    /**
     * Intialize the private OBRotorRules database from a specific file.
     */
    void Init(std::string &fname)
    {
      _rr.SetFilename(fname);
      _rr.Init();
    }
    /**
     * Turn off debugging output.
     */
    void SetQuiet() {
      _quiet=true;
      _rr.Quiet();
    }
    //! \brief Set the atoms to rotate from the dihedral atoms for each rotor
    //! Uses OBRotor->GetDihedralAtoms() to call OBRotor->SetRotAtoms()
    //! and standarizes the dihedral angles via OBRotor->SetDihedralAtoms()
    //! \return True
    bool SetRotAtoms(OBMol&);
    /**
     * Find all potentially rotatable bonds in the molecule. This method uses
     * OBBond::IsRotor() for initial evaluation which depends on ring perception
     * (i.e. ring bonds are considered rotatable). Fixed bonds, specified using the
     * deprecated fixed atoms or the new fixed bonds methods are not added to
     * the list. All rotatable bonds will be sorted by their graph theoretical
     * distance (GTD) score (see OBMol::GetGTDVector()). This results in the
     * the rotors going from the inside to the outside of the mol.
     * @param mol The molecule.
     * @param sampleRingBonds whether to sample ring bonds from analysis (default = false)
     * @return True.
     */
    bool FindRotors(OBMol &mol, bool sampleRingBonds = false);
    //! Determines which atoms should be used to calculate the internal energy
    //! if the dihedral angle of the rotor is modified
    //! \return True
    bool SetEvalAtoms(OBMol&);
    /**
     * Using the OBRotorRules database, set the torsion values (and delta)
     * to be evaluated and tested. This method also sets the rotatable atoms
     * for the rotors to the smallest possible set. For each bond there are
     * two candidate sets, one on either side. The smallest of these is used
     * and the torsion indexes for the rotor are inverted if needed
     * (i.e. a-b-c-d -> d-c-b-a).
     */
    bool AssignTorVals(OBMol &);
    ///@}

    //! \name Iterator methods
    //@{
    /**
     * Initialize the @p i iterator and get a pointer to the first OBRotor.
     * @param i OBRotorIterator object.
     */
    OBRotor *BeginRotor(OBRotorIterator &i)
    {
      i = _rotor.begin();
      return((i ==_rotor.end()) ? NULL:*i);
    }
    /**
     * Get a pointer to the next iterator.
     * @param i OBRotorIterator object.
     */
    OBRotor *NextRotor(OBRotorIterator &i)
    {
      ++i;
      return((i ==_rotor.end()) ? NULL:*i);
    }
    /**
     * Get the rotor list begin iterator.
     */
    OBRotorIterator BeginRotors()   { return(_rotor.begin()); }
    /**
     * Get the rotor list end iterator.
     */
    OBRotorIterator EndRotors()     { return(_rotor.end());   }
    //@}

    ///@name Deprecated
    ///@{
    // Not declared
    //! \deprecated Not declared. Use Setup() for top-level functionality
    bool   IdentifyEvalAtoms(OBMol &mol) { return SetEvalAtoms(mol); }
    /**
     * Set the list of fixed (invariant) atoms to the supplied OBBitVec
     * @deprecated See SetFixedBonds()
     */
    void SetFixAtoms(OBBitVec &fix)
    {
      _fixedatoms = fix;
      _fixedbonds.Clear();
    }
    /**
     * @return whether this rotor list has any fixed (invariant) atoms
     * @deprecated See HasFixedBonds()
     */
    bool HasFixedAtoms()
    {
      return(!_fixedatoms.IsEmpty());
    }
    //! Has no effect
    //! \deprecated Currently has no effect
    void IgnoreSymmetryRemoval()    { _removesym = false;}
    //! \brief Set the atoms to rotate from the dihedral atoms for each rotor
    //! Insures the fixed atoms are respected, but otherwise functions like
    //! SetRotAtoms()
    void SetRotAtomsByFix(OBMol&);
    ///@}

  };

  /// @cond DEV
  class rotor_digit {
  public:
    rotor_digit(unsigned int rs)
      {
        resolution_size = rs;
        state = 0;
      }

    rotor_digit()
      {
        resolution_size = 0;
        state = 0;
      }

    void set_size(unsigned int rs)
    {
      resolution_size = rs;
      state = 0;
    }

    void set_state(int st)
    {
      state = st;
    }

    int get_state()
    {
      return state;
    }

    unsigned int size()
    {
      return resolution_size;
    }

    bool next()
    {
      if (state < static_cast<int>(resolution_size - 1)) {
        ++state;
        return false;
      } else
        state = 0;

      return true;
    }
  private:
    unsigned int resolution_size;
    int state;
#ifndef SWIG
  } typedef rotor_digit;
#else
};
#endif
  /// @endcond

  //! \class OBRotorKeys rotor.h <openbabel/rotor.h>
  //! \brief A class to generate all possible rotorKeys
  class OBAPI OBRotorKeys
  {
      /**
      \brief A class to generate all possible rotorKeys

      This class can generate all possible rotor keys for a set of OBRotors
      which can all have their own resolution. Thanks to Yongjin Xu for this
      patch.

      the code blow is taken from  OBForceField::SystematicRotorSearch():
      \code
      #include <openbabel/rotor.h>
      #include <openbabel/mol.h>

      // See OBConversion class to fill the mol object.
      OBMol mol;
      OBRotorList rl;
      OBRotamerList rotamers;

      rl.Setup(_mol);
      rotamers.SetBaseCoordinateSets(_mol);
      rotamers.Setup(_mol, rl);

      cout << "number of rotatable bonds: " <<  rl.Size() << endl;

      if (!rl.Size()) { // only one conformer
        cout << "generated only one conformer" << endl;
        // exit here
      }

      OBRotorKeys rotorKeys;
      OBRotorIterator ri;
      OBRotor *rotor = rl.BeginRotor(ri);
      for (int i = 1; i < rl.Size() + 1; ++i, rotor = rl.NextRotor(ri)) { // foreach rotor
        rotorKeys.AddRotor(rotor->GetResolution().size());
      }

      while (rotorKeys.Next()) {
        std::vector<int> rotorKey = rotorKeys.GetKey();
        cout << "rotorKey = " << rotorKey[1] << " " << rotorKey[2] << endl;
        rotamers.AddRotamer(rotorKey);
      }

      rotamers.ExpandConformerList(_mol, _mol.GetConformers());
      \endcode
      **/

    public:
      //! Constructor
      OBRotorKeys()
      {
        _vr.clear();
      }

      //! Clear all rotors
      void Clear(){
        _vr.clear();
      }

      //! Number of rotor keys (= number of possible conformers)
      unsigned int NumKeys()
      {
        unsigned int numKeys = 0;

        while (Next())
          numKeys++;

        return numKeys;
      }

      //! Add a rotor
      //! \param size the rotor resolution
      void AddRotor(unsigned int size)
      {
        rotor_digit rd(size);
        _vr.push_back(rd);
      }

      //! Select the next rotor key
      //! \return true if there are more rotor keys
      bool Next()
      {
        if(_vr.size() == 0)
          return false;

        bool carry = _vr[0].next();
        unsigned int i = 1;
        while (carry) {
          if(i == _vr.size())
            return false;

          carry = _vr[i].next();
          i++;
        }
        return true;
      }

      //! Get the currently selected rotor key
      //! \return current rotor key
      std::vector<int> GetKey()
      {
        std::vector<int> rt;
        rt.clear();
        rt.push_back(0);
        for(unsigned int i = 0; i < _vr.size(); i++){
          rt.push_back(_vr[i].get_state());
        }

        return rt;
      }

    private:
      std::vector<rotor_digit> _vr;
  };


} // end namespace OpenBabel

#endif // OB_ROTOR_H

//! \file rotor.h
//! \brief Rotate torsional according to rotor rules.

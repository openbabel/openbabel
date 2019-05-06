/**********************************************************************
mol.h - Handle molecules. Declarations of OBMol, OBAtom, OBBond, OBResidue.
        (the main header for Open Babel)

Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2006 by Geoffrey R. Hutchison
Some portions Copyright (C) 2003 by Michael Banck

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

#ifndef OB_MOL_H
#define OB_MOL_H

#include <openbabel/babelconfig.h>

#ifndef EXTERN
#  define EXTERN extern
#endif
#ifndef THREAD_LOCAL
#ifdef SWIG
# define THREAD_LOCAL
# elif (__cplusplus >= 201103L) 
//this is required for correct multi-threading
#  define THREAD_LOCAL thread_local
# else
#  define THREAD_LOCAL
# endif
#endif

#include <math.h>
#include <float.h>

#include <vector>
#include <string>
#include <map>

#include <openbabel/base.h>


namespace OpenBabel
{
  class OBAtom;
  class OBBond;
  class OBResidue;
  class OBRing;
  class OBInternalCoord;
  class OBConversion; //used only as a pointer

  class vector3;
  class OBBitVec;
  class OBMolAtomDFSIter;
  class OBChainsParser;

  typedef std::vector<OBAtom*>::iterator OBAtomIterator;
  typedef std::vector<OBBond*>::iterator OBBondIterator;
  typedef std::vector<OBResidue*>::iterator OBResidueIterator;

  // Class OBMol
  //MOL Property Macros (flags) -- 32+ bits
  //! Smallest Set of Smallest Rings (SSSR) done. See OBRing and OBMol::FindSSSR
#define OB_SSSR_MOL              (1<<1)
  //! Ring flags have been set: See OBRing::FindRingAtomsAndBonds
#define OB_RINGFLAGS_MOL         (1<<2)
  //! Aromatic flags have been set for atoms and bonds
#define OB_AROMATIC_MOL          (1<<3)
  //! Atom typing has been performed. See OBAtomTyper
#define OB_ATOMTYPES_MOL         (1<<4)
  //! Chirality detection has been performed.
#define OB_CHIRALITY_MOL         (1<<5)
  //! Partial charges have been set or percieved
#define OB_PCHARGE_MOL           (1<<6)
  //! Atom hybridizations have been set. See OBAtomTyper
#define OB_HYBRID_MOL            (1<<8)
  //! Ring "closure" bonds have been set. See OBBond::IsClosure
#define OB_CLOSURE_MOL           (1<<11)
  //! Hyrdogen atoms have been added where needed. See OBMol::AddHydrogens
#define OB_H_ADDED_MOL           (1<<12)
  //! pH correction for hydrogen addition has been performed.
#define OB_PH_CORRECTED_MOL      (1<<13)
  //! Biomolecular chains and residues have been set. See OBChainsParser
#define OB_CHAINS_MOL            (1<<15)
  //! Total charge on this molecule has been set. See OBMol::SetTotalCharge
#define OB_TCHARGE_MOL		       (1<<16)
  //! Total spin on this molecule has been set. See OBMol::SetTotalSpinMultiplicity
#define OB_TSPIN_MOL             (1<<17)
  //! Ring typing has been performed. See OBRingTyper
#define OB_RINGTYPES_MOL         (1<<18)
  //! A pattern, not a complete molecule.
#define OB_PATTERN_STRUCTURE     (1<<19)
  //! Largest Set of Smallest Rings (LSSR) done. See OBRing and OBMol::FindLSSR
#define OB_LSSR_MOL              (1<<20)
  //! SpinMultiplicities on atoms have been set in OBMol::AssignSpinMultiplicity()
#define OB_ATOMSPIN_MOL          (1<<21)
  //! Treat as reaction
#define OB_REACTION_MOL          (1<<22)
  // flags 22-32 unspecified

#define SET_OR_UNSET_FLAG(X) \
  if (value) SetFlag(X); \
  else     UnsetFlag(X);

#define OB_CURRENT_CONFORMER	 -1

enum HydrogenType { AllHydrogen, PolarHydrogen, NonPolarHydrogen };

  // class introduction in mol.cpp
 class OBAPI OBMol: public OBBase
  {
  protected:
    int                           _flags;	//!< bitfield of flags
    bool                          _autoPartialCharge;//!< Assign partial charges automatically
    bool                          _autoFormalCharge;//!< Assign formal charges automatically
    std::string                   _title;     	//!< Molecule title
    std::vector<OBAtom*>          _vatom;      	//!< vector of atoms
    std::vector<OBAtom*>          _atomIds;    	//!< vector of atoms indexed by id
    std::vector<OBBond*>          _vbond;      	//!< vector of bonds
    std::vector<OBBond*>          _bondIds;     //!< vector of bonds
    unsigned short int            _dimension;   //!< Dimensionality of coordinates
    int				  _totalCharge; //!< Total charge on the molecule
    unsigned int                  _totalSpin;   //!< Total spin on the molecule (if not specified, assumes lowest possible spin)
    double                        *_c;	        //!< coordinate array
    std::vector<double*>          _vconf;       //!< vector of conformers
    double                        _energy;      //!< heat of formation
    unsigned int                  _natoms;      //!< Number of atoms
    unsigned int                  _nbonds;      //!< Number of bonds
    std::vector<OBResidue*>       _residue;     //!< Residue information (if applicable)
    std::vector<OBInternalCoord*> _internals;   //!< Internal Coordinates (if applicable)
    unsigned short int            _mod;	        //!< Number of nested calls to BeginModify()

  public:

    //! \name Initialization and data (re)size methods
    //@{
    //! Constructor
    OBMol();
    //! Copy constructor, copies atoms,bonds and OBGenericData
    OBMol(const OBMol &);
    //! Destructor
    virtual ~OBMol();
    //! Assignment, copies atoms,bonds and OBGenericData
    OBMol &operator=(const OBMol &mol);
    //! Copies atoms and bonds but not OBGenericData
    OBMol &operator+=(const OBMol &mol);

    //! Reserve a minimum number of atoms for internal storage
    //! This improves performance since the internal atom vector does not grow.
    void ReserveAtoms(int natoms)
    {
      if (natoms > 0 && _mod) {
        _vatom.reserve(natoms);
        _atomIds.reserve(natoms);
      }
    }

    //! Free an OBAtom pointer if defined. Does no bookkeeping
    //! \see DeleteAtom which ensures internal connections
    virtual void DestroyAtom(OBAtom*);
    //! Free an OBBond pointer if defined. Does no bookkeeping
    //! \see DeleteBond which ensures internal connections
    virtual void DestroyBond(OBBond*);
    //! Free an OBResidue pointer if defined. Does no bookkeeping
    //! \see DeleteResidue which ensures internal connections
    virtual void DestroyResidue(OBResidue*);

    //! Add the specified atom to this molecule
    //! \param atom        the atom to add
    //! \param forceNewId  whether to make a new atom Id even if the atom already has one (default is false)
    //! \return Whether the method was successful
    bool AddAtom(OBAtom& atom, bool forceNewId = false);
    //! Add a new atom to this molecule (like AddAtom)
    //! Calls BeginModify() before insertion and EndModify() after insertion
    bool InsertAtom(OBAtom &);
    //! Add a new bond to the molecule with the specified parameters
    //! \param beginIdx  the atom index of the "start" atom
    //! \param endIdx    the atom index of the "end" atom
    //! \param order     the bond order (see OBBond::GetBondOrder())
    //! \param flags     any bond flags such as stereochemistry (default = none)
    //! \param insertpos the position index to insert the bond (default = none)
    //! \return Whether the new bond creation was successful
    bool AddBond(int beginIdx, int endIdx, int order,
                 int flags=0,int insertpos=-1);
    //! Add the specified residue to this molecule and update connections
    //! \return Whether the method was successful
    bool AddBond(OBBond&);
    //! Add the specified residue to this molecule and update connections
    //! \return Whether the method was successful
    bool AddResidue(OBResidue&);

    //! Create a new OBAtom in this molecule and ensure connections
    //! (e.g. OBAtom::GetParent(). A new unique id will be assigned
    //! to this atom.
    OBAtom    *NewAtom();
    //! Create a new OBAtom in this molecule and ensure connections.
    //! (e.g. OBAtom::GetParent(). The @p id will be assigned to this
    //! atom.
    OBAtom    *NewAtom(unsigned long id);
    //! Create a new OBBond in this molecule and ensure connections
    //! (e.g. OBBond::GetParent(). A new unique id will be assigned
    //! to this bond.
    OBBond    *NewBond();
    //! Create a new OBBond in this molecule and ensure connections
    //! (e.g. OBBond::GetParent(). The @p id will be assigned to this
    //! bond.
    OBBond    *NewBond(unsigned long id);
    //! Create a new OBResidue in this molecule and ensure connections.
    OBResidue *NewResidue();
    //! Deletes an atom from this molecule and all appropriate bonds.
    //! Updates the molecule and atom and bond indexes accordingly.
    //! \warning Does not update any residues which may contain this atom
    //! \return Whether deletion was successful
    bool DeleteAtom(OBAtom*, bool destroyAtom = true);
    //! Deletes an bond from this molecule and updates accordingly
    //! \return Whether deletion was successful
    bool DeleteBond(OBBond*, bool destroyBond = true);
    //! Deletes a residue from this molecule and updates accordingly.
    //! \return Whether deletion was successful
    bool DeleteResidue(OBResidue*, bool destroyResidue = true);
    //@}

    //! \name Molecule modification methods
    //@{
    //! Call when making many modifications -- clears conformer/rotomer data.
    //! The method "turns off" perception routines, improving performance.
    //! Changes in molecular structure will be re-considered after modifications.
    virtual void BeginModify(void);
    //! Call when done with modificaions -- re-perceive data as needed.
    //! This method "turns on" perception routines and re-evaluates molecular
    //! structure.
    virtual void EndModify(bool nukePerceivedData=true);
    //! \return The number of nested BeginModify() calls. Used internally.
    int GetMod()           {      return(_mod);    }
    //! Increase the number of nested BeginModify calls. Dangerous!
    //! Instead, properly use BeginModify as needed.
    void IncrementMod()    {      _mod++;          }
    //! Decrease the number of nested BeginModify calls. Dangerous!
    //! Instead, properly use EndModify as needed.
    void DecrementMod()    {      _mod--;          }
    //@}

    //! \name Data retrieval methods
    //@{
    //! \return the entire set of flags. (Internal use, mainly.)
    int          GetFlags()               { return(_flags); }
    //! \return the title of this molecule (often the filename)
    //! \param replaceNewlines whether to replace any newline characters with spaces
    const char  *GetTitle(bool replaceNewlines = true) const;
    //! \return the number of atoms (i.e. OBAtom children)
    unsigned int NumAtoms() const         {  return(_natoms); }
    //! \return the number of bonds (i.e. OBBond children)
    unsigned int NumBonds() const         {  return(_nbonds); }
    //! \return the number of non-hydrogen atoms
    unsigned int NumHvyAtoms();
    //! \return the number of residues (i.e. OBResidue substituents)
    unsigned int NumResidues() const      { return(static_cast<unsigned int> (_residue.size())); }
    //! \return the number of rotatable bonds. If sampleRingBonds is true, will include rotors within rings (see OBBond::IsRotor() for details)
    unsigned int NumRotors(bool sampleRingBonds=false);

    //! \return the atom at index @p idx or NULL if it does not exist.
    //! \warning Atom indexing will change. Use iterator methods instead.
    OBAtom      *GetAtom(int idx) const;
    //! \return the atom with @p id or NULL if it does not exist.
    OBAtom      *GetAtomById(unsigned long id) const;
    //! \return the first atom in this molecule, or NULL if none exist.
    //! \deprecated Will be removed in favor of more standard iterator methods
    OBAtom      *GetFirstAtom() const;
    //! \return the bond at index @p idx or NULL if it does not exist.
    //! \warning Bond indexing may change. Use iterator methods instead.
    OBBond      *GetBond(int idx) const;
    //! \return the bond with @p id or NULL if it does not exist.
    OBBond      *GetBondById(unsigned long id) const;
    //! \return the bond connecting the atom indexed by @p a and @p b or NULL if none exists.
    //! \warning Atom indexing will change. Use atom objects and iterators instead.
    OBBond      *GetBond(int a, int b) const;
    // The safer version of the above method
    //! \return the bond between the atoms @p bgn and @p end or NULL if none exists
    OBBond      *GetBond(OBAtom* bgn, OBAtom* end) const;
    //! \return the residue indexed by @p idx, or NULL if none exists
    //! \warning Residue indexing may change. Use iterator methods instead.
    OBResidue   *GetResidue(int idx) const;
    std::vector<OBInternalCoord*> GetInternalCoord();
    /*! \return the dihedral angle (in degrees) between the four atoms supplied a1-a2-a3-a4)
     *  WARNING: SetTorsion takes an angle in radians while GetTorsion returns it
     *  in degrees
     */
    double       GetTorsion(int,int,int,int);
    /*! \return the dihedral angle (in degrees) between the four atoms @p a, @p b, @p c, and @p d)
     *  WARNING: SetTorsion takes an angle in radians while GetTorsion returns it
     *  in degrees
     */
    double       GetTorsion(OBAtom* a,OBAtom* b,OBAtom* c,OBAtom* d);
    //! \return the angle (in degrees) between the three atoms @p a, @p b and @p c
    //!  (where  a-> b (vertex) -> c )
    double GetAngle(OBAtom* a, OBAtom* b, OBAtom* c);
    //! \return the size of the smallest ring if a and b are in the same ring, 0 otherwise
    //! \since version 2.4
    int AreInSameRing(OBAtom *a, OBAtom *b);
    //! \return the stochoimetric formula (e.g., C4H6O)
    std::string  GetFormula();
    //! \return the stochoimetric formula in spaced format e.g. C 4 H 6 O 1
    std::string  GetSpacedFormula(int ones=0, const char* sp=" ", bool implicitH = true);
    //! \return the heat of formation for this molecule (in kcal/mol)
    double       GetEnergy() const { return _energy; }
    //! \return the standard molar mass given by IUPAC atomic masses (amu)
    double       GetMolWt(bool implicitH = true);
    //! \return the mass given by isotopes (or most abundant isotope, if not specified)
    double	 GetExactMass(bool implicitH = true);
    //! \return the total charge on this molecule (i.e., 0 = neutral, +1, -1...)
    int		 GetTotalCharge();
    //! \return the total spin on this molecule (i.e., 1 = singlet, 2 = doublet...)
    unsigned int GetTotalSpinMultiplicity();
    //! \return the dimensionality of coordinates (i.e., 0 = unknown or no coord, 2=2D, 3=3D)
    unsigned short int GetDimension() const { return _dimension; }
    //! \return the set of all atomic coordinates. See OBAtom::GetCoordPtr for more
    double      *GetCoordinates() { return(_c); }
    //! \return the Smallest Set of Smallest Rings has been run (see OBRing class)
    std::vector<OBRing*> &GetSSSR();
    //! \return the Largest Set of Smallest Rings has been run (see OBRing class)
    std::vector<OBRing*> &GetLSSR();
    //! Get the current flag for whether formal charges are set with pH correction
    bool AutomaticFormalCharge()   { return(_autoFormalCharge);  }
    //! Get the current flag for whether partial charges are auto-determined
    bool AutomaticPartialCharge()  { return(_autoPartialCharge); }
    //@}


    //! \name Data modification methods
    //@{
    //! Set the title of this molecule to @p title
    void   SetTitle(const char *title);
    //! Set the title of this molecule to @p title
    void   SetTitle(std::string &title);
    //! Set the stochiometric formula for this molecule
    void   SetFormula(std::string molFormula);
    //! Set the heat of formation for this molecule (in kcal/mol)
    void   SetEnergy(double energy) { _energy = energy; }
    //! Set the dimension of this molecule (i.e., 0, 1 , 2, 3)
    void   SetDimension(unsigned short int d) { _dimension = d; }
    //! Set the total charge of this molecule to @p charge
    void   SetTotalCharge(int charge);
    //! Set the total spin multiplicity of this molecule to @p spinMultiplicity
    //! Overrides the calculation from spin multiplicity of OBAtoms
    void   SetTotalSpinMultiplicity(unsigned int spinMultiplicity);
    //! Set the internal coordinates to @p int_coord
    //! (Does not call InternalToCartesian to update the 3D cartesian
    //! coordinates).
    //! The size of the @p int_coord has to be the same as the number of atoms
    //! in molecule (+ NULL at the beginning).
    void SetInternalCoord(std::vector<OBInternalCoord*> int_coord);
    //! Set the flag for determining automatic formal charges with pH (default=true)
    void SetAutomaticFormalCharge(bool val)
    { _autoFormalCharge=val;  }
    //! Set the flag for determining partial charges automatically (default=true)
    void SetAutomaticPartialCharge(bool val)
    { _autoPartialCharge=val; }

    //! Mark that aromaticity has been perceived for this molecule (see OBAromaticTyper)
    void   SetAromaticPerceived(bool value = true)    { SET_OR_UNSET_FLAG(OB_AROMATIC_MOL);    }
    //! Mark that Smallest Set of Smallest Rings has been run (see OBRing class)
    void   SetSSSRPerceived(bool value = true)        { SET_OR_UNSET_FLAG(OB_SSSR_MOL);        }
    //! Mark that Largest Set of Smallest Rings has been run (see OBRing class)
    void   SetLSSRPerceived(bool value = true)        { SET_OR_UNSET_FLAG(OB_LSSR_MOL);        }
    //! Mark that rings have been perceived (see OBRing class for details)
    void   SetRingAtomsAndBondsPerceived(bool value = true) { SET_OR_UNSET_FLAG(OB_RINGFLAGS_MOL); }
    //! Mark that atom types have been perceived (see OBAtomTyper for details)
    void   SetAtomTypesPerceived(bool value = true)   { SET_OR_UNSET_FLAG(OB_ATOMTYPES_MOL);   }
    //! Mark that ring types have been perceived (see OBRingTyper for details)
    void   SetRingTypesPerceived(bool value = true)   { SET_OR_UNSET_FLAG(OB_RINGTYPES_MOL);   }
    //! Mark that chains and residues have been perceived (see OBChainsParser)
    void   SetChainsPerceived(bool value = true)      { SET_OR_UNSET_FLAG(OB_CHAINS_MOL);      }
    //! Mark that chirality has been perceived
    void   SetChiralityPerceived(bool value = true)   { SET_OR_UNSET_FLAG(OB_CHIRALITY_MOL);   }
    //! Mark that partial charges have been assigned
    void   SetPartialChargesPerceived(bool value = true) { SET_OR_UNSET_FLAG(OB_PCHARGE_MOL);  }
    //! Mark that hybridization of all atoms has been assigned
    void   SetHybridizationPerceived(bool value = true)  { SET_OR_UNSET_FLAG(OB_HYBRID_MOL);   }
    //! Mark that ring closure bonds have been assigned by graph traversal
    void   SetClosureBondsPerceived(bool value = true)   { SET_OR_UNSET_FLAG(OB_CLOSURE_MOL);  }
    //! Mark that explicit hydrogen atoms have been added
    void   SetHydrogensAdded(bool value = true) { SET_OR_UNSET_FLAG(OB_H_ADDED_MOL); }
    void   SetCorrectedForPH(bool value = true) { SET_OR_UNSET_FLAG(OB_PH_CORRECTED_MOL); }
    void   SetSpinMultiplicityAssigned(bool value = true) { SET_OR_UNSET_FLAG(OB_ATOMSPIN_MOL); }
    //! The OBMol is a pattern, not a complete molecule. Left unchanged by Clear().
    void   SetIsPatternStructure(bool value = true) { SET_OR_UNSET_FLAG(OB_PATTERN_STRUCTURE); }
    void   SetIsReaction(bool value = true)               { SET_OR_UNSET_FLAG(OB_REACTION_MOL) };
    bool   HasFlag(int flag)   { return (_flags & flag) ? true : false; }
    void   SetFlag(int flag)   { _flags |= flag; }
    void   UnsetFlag(int flag) { _flags &= (~(flag)); }
    void   SetFlags(int flags) { _flags = flags; }
    //@}

    //! \name Molecule modification methods
    //@{
    // Description in transform.cpp (command-line transformations to this molecule)
    virtual OBBase*    DoTransformations(const std::map<std::string,std::string>* pOptions,OBConversion* pConv);
    // Ditto (documentation on transformation options)
    static const char* ClassDescription();
    //! Clear all information from a molecule except OB_PATTERN_STRUCTURE left unchanged
    bool Clear();
    //! Renumber the atoms of this molecule according to the order in the supplied vector
    void RenumberAtoms(std::vector<OBAtom*>&);
    //! Renumber the atoms of this molecule using the initial indexes in the supplied vector
    void RenumberAtoms(std::vector<int>);
    //! Set the coordinates for all atoms in this conformer.
    //! \sa OBMol::GetCoordinates()
    void SetCoordinates(double *c);
    //! Translate one conformer and rotate by a rotation matrix (which is returned) to the inertial frame-of-reference
    void ToInertialFrame(int conf, double *rmat);
    //! Translate all conformers to the inertial frame-of-reference
    void ToInertialFrame();
    //! Translates all conformers in the molecule by the supplied vector
    void Translate(const vector3 &v);
    //! Translates one conformer in the molecule by the supplied vector
    void Translate(const vector3 &v, int conf);
    //! Rotate all conformers using the supplied matrix @p u (a 3x3 array of double)
    void Rotate(const double u[3][3]);
    //! Rotate all conformers using the supplied matrix @p m (a linear 3x3 row-major array of double)
    void Rotate(const double m[9]);
    //! Rotate a specific conformer @p nconf using the supplied rotation matrix @p m
    void Rotate(const double m[9],int nconf);
    //! Translate to the center of all coordinates (for this conformer)
    void Center();
    //! Delete all hydrogens from the molecule
    //! \return Success
    bool DeleteHydrogens();
    //! Delete all hydrogens from the supplied atom
    //! \return Success
    bool DeleteHydrogens(OBAtom*);
    //! Delete all hydrogen atoms connected to a polar atom
    //! \see OBAtom::IsPolarHydrogen
    //! \since version 2.4
    bool DeletePolarHydrogens();
    //! Delete all hydrogen atoms connected to a non-polar atom
    //! \see OBAtom::IsNonPolarHydrogen
    bool DeleteNonPolarHydrogens();
    //! Delete the supplied atom if it is a hydrogen
    //! (Helper function for DeleteHydrogens)
    bool DeleteHydrogen(OBAtom*);
    //! Add hydrogens to the entire molecule to fill out implicit valence spots
    //! \param polaronly    Whether to add hydrogens only to polar atoms
    //! (i.e., not to C atoms)
    //! \param correctForPH Whether to call CorrectForPH() first
    //! \param pH The pH to use for CorrectForPH() modification
    //! \return Whether any hydrogens were added
    bool AddHydrogens(bool polaronly=false,bool correctForPH=false, double pH=7.4);
    //! Add hydrogens only to the supplied atom to fill out implicit valence
    bool AddHydrogens(OBAtom*);
    //! Add only polar hydrogens (i.e., attached to polar atoms, not C)
    bool AddPolarHydrogens();
    //! Add only nonpolar hydrogens (i.e., attached to C)
    //! \since version 2.4
    bool AddNonPolarHydrogens();
    //! Add polar and/or nonpolar hydrogens
    //! \since verison 2.4
    bool AddNewHydrogens(HydrogenType whichHydrogen, bool correctForPH=false, double pH=7.4);

    //! If @p threshold is not specified or is zero, remove all but the largest
    //! contiguous fragment. If @p threshold is non-zero, remove any fragments with fewer
    //! than @p threshold atoms.
    bool StripSalts(unsigned int threshold=0);
    //! Copies each disconnected fragment as a separate OBMol
    std::vector<OBMol> Separate(int StartIndex=1);
    //! Iterative component of Separate to copy one fragment at a time
    bool GetNextFragment( OpenBabel::OBMolAtomDFSIter& iter, OBMol& newMol );
    // docs in mol.cpp
    bool CopySubstructure(OBMol& newmol, OBBitVec *includeatoms, OBBitVec *excludebonds = (OBBitVec*)0,
      unsigned int correctvalence=1,
      std::vector<unsigned int> *atomorder=(std::vector<unsigned int>*)0,
      std::vector<unsigned int> *bondorder=(std::vector<unsigned int>*)0);
    //! Converts the charged form of coordinate bonds, e.g.[N+]([O-])=O to N(=O)=O
    bool ConvertDativeBonds();
    //! Converts 5-valent N and P only. Return true if conversion occurred.
    //! \return has charged form of dative bonds(e.g.[N+]([O-])=O from N(=O)=O).
    //! \since version 2.4
    bool MakeDativeBonds();
    /** Convert zero-order bonds to single or double bonds and adjust adjacent atom
     *  charges in an attempt to achieve the correct valence state.
     *  @return Whether any modifications were made
     *  @since version 2.4
     */
    bool ConvertZeroBonds();

    //! Correct for pH by applying the OBPhModel transformations
    bool CorrectForPH(double pH=7.4);
    // docs in mol.cpp
    bool AssignSpinMultiplicity(bool NoImplicitH=false);

    //! Put the specified molecular charge on appropriate atoms.
    //! Assumes all the hydrogen is explicitly included in the molecule.
    //! \since version 2.4
    bool AssignTotalChargeToAtoms(int charge);

    //! \return the center of the supplied conformer @p nconf
    //! \see Center() to actually center all conformers at the origin
    vector3 Center(int nconf);
    /*! Set the torsion defined by these atoms, rotating bonded neighbors
     *  \par ang The torsion angle in radians
     *  WARNING: SetTorsion takes an angle in radians while GetTorsion returns it
     *  in degrees
     */
    void SetTorsion(OBAtom*,OBAtom*,OBAtom*,OBAtom*,double ang);
    //@}

    //! \name Molecule utilities and perception methods
    //@{
    //! Find Smallest Set of Smallest Rings (see OBRing class for more details)
    void FindSSSR();
    //! Find Largest Set of Smallest Rings
    void FindLSSR();
    //! Find all ring atoms and bonds. Does not need to call FindSSSR().
    void FindRingAtomsAndBonds();
    // documented in mol.cpp -- locates all atom indexes which can reach 'end'
    void FindChildren(std::vector<int> & children,int bgnIdx,int endIdx);
    // documented in mol.cpp -- locates all atoms which can reach 'end'
    void FindChildren(std::vector<OBAtom*>& children,OBAtom* bgn,OBAtom* end);
    //! Find the largest fragment in OBMol
    //! (which may include multiple non-connected fragments)
    //! \param frag   Return (by reference) a bit vector indicating the atoms
    //! in the largest fragment
    void FindLargestFragment(OBBitVec &frag);
    //! Sort a list of contig fragments by size from largest to smallest
    //! Each vector<int> contains the atom numbers of a contig fragment
    void ContigFragList(std::vector<std::vector<int> >&);
    //! Aligns atom a on p1 and atom b along p1->p2 vector
    void Align(OBAtom*,OBAtom*,vector3&,vector3&);
    //! Adds single bonds based on atom proximity
    void ConnectTheDots();
    //! Attempts to perceive multiple bonds based on geometries
    void PerceiveBondOrders();
    //! Fills out an OBAngleData with angles from the molecule
    void FindAngles();
    //! Fills out an OBTorsionData with angles from the molecule
    void FindTorsions();
    // documented in mol.cpp: graph-theoretical distance for each atom
    bool         GetGTDVector(std::vector<int> &);
    // documented in mol.cpp: graph-invariant index for each atom
    void         GetGIVector(std::vector<unsigned int> &);
    // documented in mol.cpp: calculate symmetry-unique identifiers
    void         GetGIDVector(std::vector<unsigned int> &);
    //@}

    //! \name Methods to check for existence of properties
    //@{
    //! Are there non-zero coordinates in two dimensions (i.e. X and Y)- and, if Not3D is true, no Z coordinates?
    bool Has2D(bool Not3D=false);
    //! Are there non-zero coordinates in all three dimensions (i.e. X, Y, Z)?
    bool Has3D();
    //! Are there any non-zero coordinates?
    bool HasNonZeroCoords();
    //! Has aromatic perception been performed?
    bool HasAromaticPerceived()     { return(HasFlag(OB_AROMATIC_MOL)); }
    //! Has the smallest set of smallest rings (FindSSSR) been performed?
    bool HasSSSRPerceived()         { return(HasFlag(OB_SSSR_MOL));     }
    //! Has the largest set of smallest rings (FindLSSR) been performed?
    bool HasLSSRPerceived()         { return(HasFlag(OB_LSSR_MOL));     }
    //! Have ring atoms and bonds been assigned?
    bool HasRingAtomsAndBondsPerceived(){return(HasFlag(OB_RINGFLAGS_MOL));}
    //! Have atom types been assigned by OBAtomTyper?
    bool HasAtomTypesPerceived()    { return(HasFlag(OB_ATOMTYPES_MOL));}
    //! Have ring types been assigned by OBRingTyper?
    bool HasRingTypesPerceived()    { return(HasFlag(OB_RINGTYPES_MOL));}
    //! Has atom chirality been assigned?
    bool HasChiralityPerceived()    { return(HasFlag(OB_CHIRALITY_MOL));}
    //! Have atomic Gasteiger partial charges been assigned by OBGastChrg?
    bool HasPartialChargesPerceived() { return(HasFlag(OB_PCHARGE_MOL));}
    //! Has atomic hybridization been assigned by OBAtomTyper?
    bool HasHybridizationPerceived() { return(HasFlag(OB_HYBRID_MOL));  }
    //! Have ring "closure" bonds been assigned? (e.g., OBBond::IsClosure())
    bool HasClosureBondsPerceived() { return(HasFlag(OB_CLOSURE_MOL));  }
    //! Have biomolecule chains and residues been assigned by OBChainsParser?
    bool HasChainsPerceived() { return(HasFlag(OB_CHAINS_MOL));         }
    //! Have hydrogens been added to the molecule?
    bool HasHydrogensAdded() { return(HasFlag(OB_H_ADDED_MOL));         }
    //! Has the molecule been corrected for pH by CorrectForPH?
    bool IsCorrectedForPH() { return(HasFlag(OB_PH_CORRECTED_MOL));     }
    //! Has total spin multiplicity been assigned?
    bool HasSpinMultiplicityAssigned() { return(HasFlag(OB_ATOMSPIN_MOL)); }
    //! Does this OBMol represent a reaction?
    bool IsReaction()                  { return HasFlag(OB_REACTION_MOL); }
    //! Are there any atoms in this molecule?
    bool Empty()                       { return(_natoms == 0);          }
    //@}

    //! \name Multiple conformer member functions
    //@{
    //! \return the number of conformers in this molecule
    int     NumConformers()    { return((_vconf.empty())?0:static_cast<int> (_vconf.size())); }
    //! Set the entire set of conformers for this molecule to @p v
    void    SetConformers(std::vector<double*> &v);
    //! Add a new set of coordinates @p f as a new conformer
    void    AddConformer(double *f)    {  _vconf.push_back(f);    }
    //! Set the molecule's current conformer to @p i
    //! Does nothing if @p i is larger than NumConformers()
    void    SetConformer(unsigned int i);
    //! Copy the conformer @p nconf into the array @p c
    //! \warning Does no checking to see if @p c is large enough
    void    CopyConformer(double* c,int nconf);
    //! Delete the conformer @p nconf
    void    DeleteConformer(int nconf);
    //! \return the coordinates to conformer @p i
    double  *GetConformer(int i)       {  return(_vconf[i]);      }
    //! Set the entire set of conformer energies
    void    SetEnergies(std::vector<double> &energies);
    //! Set the entire set of conformer energies
    std::vector<double> GetEnergies();
    //! Get the energy for conformer ci
    //! \par ci conformer index
    double  GetEnergy(int ci);
    //! Set the iterator to the beginning of the conformer list
    //! \return the array of coordinates for the first conformer
    double  *BeginConformer(std::vector<double*>::iterator&i)
    { i = _vconf.begin();
      return((i == _vconf.end()) ? NULL:*i); }
    //! Advance the iterator to the next confomer, if possible
    //! \return The array of coordinates for the next conformer, or NULL if none exist
    double  *NextConformer(std::vector<double*>::iterator&i)
    { ++i;
      return((i == _vconf.end()) ? NULL:*i); }
    //! \return the entire set of conformers for this molecule as a vector of floating point arrays
    std::vector<double*> &GetConformers() {   return(_vconf);     }
    //@}

    //! \name Iterator methods
    //@{
    //! \return An atom iterator pointing to the beginning of the atom list
    OBAtomIterator BeginAtoms()   { return _vatom.begin(); }
    //! \return An atom iterator pointing to the end of the atom list
    OBAtomIterator EndAtoms()    { return _vatom.begin() + NumAtoms() ; }
    //! \return A bond iterator pointing to the beginning of the bond list
    OBBondIterator BeginBonds()   { return _vbond.begin(); }
    //! \return A bond iterator pointing to the end of the bond list
    OBBondIterator EndBonds()     { return _vbond.begin() + NumBonds() ; }
    //! \return A residue iterator pointing to the beginning of the residue list
    OBResidueIterator BeginResidues() { return _residue.begin(); }
    //! \return A residue iterator pointing to the end of the residue list
    OBResidueIterator EndResidues()   { return _residue.end();   }

    //! Set the iterator @p i to the beginning of the atom list
    //! \return the first atom (or NULL if none exist)
    OBAtom *BeginAtom(OBAtomIterator &i);
    //! Advance the iterator @p i to the next atom in the molecule
    //! \return the next atom (if any, or NULL if none exist)
    OBAtom *NextAtom(OBAtomIterator &i);
    //! Set the iterator @p i to the beginning of the bond list
    //! \return the first bond (or NULL if none exist)
    OBBond *BeginBond(OBBondIterator &i);
    //! Advance the iterator @p i to the next bond in the molecule
    //! \return the next bond (if any, or NULL if none exist)
    OBBond *NextBond(OBBondIterator &i);
    //! Set the iterator @p i to the beginning of the resdiue list
    //! \return the first residue (or NULL if none exist)
    OBResidue *BeginResidue(OBResidueIterator &i)
    {
      i = _residue.begin();
      return((i == _residue.end()) ? NULL:*i);
    }
    //! Advance the iterator @p i to the next residue in the molecule
    //! \return the next residue (if any, or NULL if not possible)
    OBResidue *NextResidue(OBResidueIterator &i)
    {
      ++i;
      return((i == _residue.end()) ? NULL:*i);
    }
    //! Set the iterator to the beginning of the internal coordinate list
    //! \return the first internal coordinate record, or NULL if none exist
    //! \see SetInternalCoord
    OBInternalCoord *BeginInternalCoord(std::vector<OBInternalCoord*>::iterator &i)
    {
      i = _internals.begin();
      return((i == _internals.end()) ? NULL:*i);
    }
    //! Advance the iterator to the next internal coordinate record
    //! \return the next first internal coordinate record, or NULL if none exist
    //! \see SetInternalCoord
    OBInternalCoord *NextInternalCoord(std::vector<OBInternalCoord*>::iterator &i)
    {
      ++i;
      return((i == _internals.end()) ? NULL:*i);
    }
    //@}

  };

  // Utility function prototypes
  //tokenize and Trim declarations moved to base.h
  // Deprecated -- use OBMessageHandler class instead (docs in obutil.cpp)
  OBAPI void ThrowError(char *str);
  // Deprecated -- use OBMessageHandler class instead (docs in obutil.cpp)
  OBAPI void ThrowError(std::string &str);
  //! Convert Cartesian XYZ to a set of OBInternalCoord coordinates
  OBAPI void CartesianToInternal(std::vector<OBInternalCoord*>&,OBMol&);
  //! Convert set of OBInternalCoord coordinates into Cartesian XYZ
  OBAPI void InternalToCartesian(std::vector<OBInternalCoord*>&,OBMol&);
  // Replace the last extension in str with a new one (docs in obutil.cpp)
  OBAPI std::string NewExtension(std::string&,char*);

  //! \brief Nested namespace for max_value templates
  namespace detail {
    //! \struct max_value mol.h <openbabel/mol.h>
    //! \brief a C++ template to return the maximum value of a type (e.g., int)
    template<typename T, int size = sizeof(T)>
    struct max_value
    {
      static const T result = (static_cast<T>(0xFF) << (size-1)*8) + max_value<T, size-1>::result;
    };

    //! \brief a C++ template to return the maximum value of a type (e.g., int)
    template<typename T>
    struct max_value<T, 0>
    {
      static const T result = 0;
    };
  }

  // No unique id
  static const unsigned long NoId = detail::max_value<unsigned long>::result;

  //Utility Macros

#ifndef BUFF_SIZE
#define BUFF_SIZE 32768
#endif

#ifndef EQ
#define EQ(a,b) (!strcmp((a), (b)))
#endif

#ifndef EQn
#define EQn(a,b,n) (!strncmp((a), (b), (n)))
#endif

#ifndef SQUARE
#define SQUARE(x) ((x)*(x))
#endif

#ifndef IsUnsatType
#define IsUnsatType(x)  (EQ(x,"Car") || EQ(x,"C2") || EQ(x,"Sox") || EQ(x,"Sac") || EQ(x,"Pac") || EQ(x,"So2"))
#endif

#ifndef __KCC
  extern "C"
  {
    OBAPI void  get_rmat(double*,double*,double*,int);
    OBAPI void  ob_make_rmat(double mat[3][3],double rmat[9]);
    OBAPI void  qtrfit (double *r,double *f,int size,double u[3][3]);
    OBAPI double superimpose(double*,double*,int);
  }
#else
  OBAPI void get_rmat(double*,double*,double*,int);
  OBAPI void ob_make_rmat(double mat[3][3],double rmat[9]);
  OBAPI void qtrfit (double *r,double *f,int size,double u[3][3]);
  OBAPI double superimpose(double*,double*,int);
#endif // __KCC

//  extern OBMol* (*CreateMolecule) (void);

} // end namespace OpenBabel

#endif // OB_MOL_H

//! \file mol.h
//! \brief Handle molecules. Declarations of OBMol, OBAtom, OBBond, OBResidue.
//!        (the main header for Open Babel)

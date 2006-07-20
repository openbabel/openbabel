/**********************************************************************
mol.h - Handle molecules. Declarations of OBMol, OBAtom, OBBond, OBResidue.
        (the main header for Open Babel)
 
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2006 by Geoffrey R. Hutchison
Some portions Copyright (C) 2003 by Michael Banck
 
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

#ifndef OB_MOL_H
#define OB_MOL_H

#include "babelconfig.h"

#ifndef EXTERN
#  define EXTERN extern
#endif

#include <math.h>

#include <vector>
#include <string>
#include <map>

// Currently includes many headers for 2.x backwards compatibility
// \deprecated -- this will be cleaned up in 3.0 efforts
#include "atom.h"
#include "bond.h"
#include "base.h"
#include "data.h"
#include "chains.h"
#include "math/vector3.h"
#include "bitvec.h"
#include "residue.h"
#include "ring.h"
#include "generic.h"
#include "typer.h"
#include "oberror.h"
#include "obiter.h"
#include "internalcoord.h"

namespace OpenBabel
{

  class OBAtom;
  class OBBond;
  class OBMol;
  class OBInternalCoord;

  // Class OBMol

  //MOL Property Macros (flags) -- 32+ bits
#define OB_SSSR_MOL              (1<<1)
#define OB_RINGFLAGS_MOL         (1<<2)
#define OB_AROMATIC_MOL          (1<<3)
#define OB_ATOMTYPES_MOL         (1<<4)
#define OB_CHIRALITY_MOL         (1<<5)
#define OB_PCHARGE_MOL           (1<<6)
#define OB_HYBRID_MOL            (1<<8)
#define OB_IMPVAL_MOL            (1<<9)
#define OB_KEKULE_MOL            (1<<10)
#define OB_CLOSURE_MOL           (1<<11)
#define OB_H_ADDED_MOL           (1<<12)
#define OB_PH_CORRECTED_MOL      (1<<13)
#define OB_AROM_CORRECTED_MOL    (1<<14)
#define OB_CHAINS_MOL            (1<<15)
#define OB_TCHARGE_MOL		 (1<<16)
#define OB_TSPIN_MOL             (1<<17)
  // flags 18-32 unspecified
#define OB_CURRENT_CONFORMER	 -1

  // class introduction in mol.cpp
  class OBAPI OBMol : public OBGraphBase
    {
    protected:
      int                           _flags;	//!< bitfield of flags
      bool                          _autoPartialCharge; //!< Assign partial charges automatically
      bool                          _autoFormalCharge; //!< Assign formal charges automatically
      std::string                   _title;	//!< Molecule title
      //vector<OBAtom*>             _atom;	//!< not needed (inherited)
      //vector<OBBond*>             _bond;	//!< not needed (inherited)
      unsigned short int            _dimension;   //!< Dimensionality of coordinates
      double                        _energy;      //!< Molecular heat of formation (if applicable)
      int				  _totalCharge; //!< Total charge on the molecule
      unsigned int                  _totalSpin;   //!< Total spin on the molecule (if not specified, assumes lowest possible spin)
      double                       *_c;	        //!< coordinate array
      std::vector<double*>          _vconf;       //!< vector of conformers
      unsigned short int            _natoms;      //!< Number of atoms
      unsigned short int            _nbonds;      //!< Number of bonds
      std::vector<OBResidue*>       _residue;     //!< Residue information (if applicable)
      std::vector<OBInternalCoord*> _internals;   //!< Internal Coordinates (if applicable)
      //    std::vector<OBGenericData*>   _vdata;       //!< Custom data -- see OBGenericData class for more
      unsigned short int            _mod;	        //!< Number of nested calls to BeginModify()

      bool  HasFlag(int flag)    { return((_flags & flag) ? true : false); }
      void  SetFlag(int flag)    { _flags |= flag; }

      //! \name Internal Kekulization routines -- see kekulize.cpp and NewPerceiveKekuleBonds()
      //@{
      void start_kekulize(std::vector <OBAtom*> &cycle, std::vector<int> &electron);
      int expand_kekulize(OBAtom *atom1, OBAtom *atom2, std::vector<int> &currentState, std::vector<int> &initState, std::vector<int> &bcurrentState, std::vector<int> &binitState, std::vector<bool> &mark);
      int getorden(OBAtom *atom);
      void expandcycle(OBAtom *atom, OBBitVec &avisit);
      //@}

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

      void ReserveAtoms(int natoms)
        {
          if (natoms && _mod)
            _vatom.reserve(natoms);
        }
      virtual OBAtom *CreateAtom(void);
      virtual OBBond *CreateBond(void);
      virtual void DestroyAtom(OBNodeBase*);
      virtual void DestroyBond(OBEdgeBase*);
      bool AddAtom(OBAtom&);
      bool AddBond(int,int,int,int flags=0,int insertpos=-1);
      bool AddBond(OBBond&);
      bool AddResidue(OBResidue&);
      bool InsertAtom(OBAtom &);
      bool DeleteAtom(OBAtom*);
      bool DeleteBond(OBBond*);
      bool DeleteResidue(OBResidue*);
      OBAtom    *NewAtom();
      OBResidue *NewResidue();
      //@}

      //! \name Molecule modification methods
      //@{
      //! Call when making many modifications -- clears conformer/rotomer data.
      virtual void BeginModify(void);
      //! Call when done with modificaions -- re-perceive data as needed.
      virtual void EndModify(bool nukePerceivedData=true);
      int GetMod()
        {
          return(_mod);
        }
      void IncrementMod()
        {
          _mod++;
        }
      void DecrementMod()
        {
          _mod--;
        }
      //@}

      //! \name Data retrieval methods
      //@{
      int          GetFlags()               { return(_flags); }
      //! \return the title of this molecule (often the filename)
      const char  *GetTitle() const         { return(_title.c_str()); }
      //! \return the number of atoms (i.e. OBAtom children)
      unsigned int NumAtoms() const         {  return(_natoms); }
      //! \return the number of bonds (i.e. OBBond children)
      unsigned int NumBonds() const         {  return(_nbonds); }
      //! \return the number of non-hydrogen atoms
      unsigned int NumHvyAtoms();
      //! \return the number of residues (i.e. OBResidue substituents)
      unsigned int NumResidues() const      { return(_residue.size()); }
      //! \return the number of rotatble bonds. See OBBond::IsRotor() for details
      unsigned int NumRotors();
    
      OBAtom      *GetAtom(int);
      OBAtom      *GetFirstAtom();
      OBBond      *GetBond(int);
      OBBond      *GetBond(int, int);
      //! \return the bond between the atoms @p bgn and @p end
      OBBond      *GetBond(OBAtom* bgn, OBAtom* end);
      OBResidue   *GetResidue(int);
      std::vector<OBInternalCoord*> GetInternalCoord();
      //! \return the dihedral angle between the four atoms supplied a1-a2-a3-a4)
      double       GetTorsion(int,int,int,int);
      //! \return the dihedral angle between the four atoms @p a, @p b, @p c, and @p d)
      double       GetTorsion(OBAtom* a,OBAtom* b,OBAtom* c,OBAtom* d);
      //! \return the angle between the three atoms @p a, @p b and @p c
      //!  (where  a-> b (vertex) -> c )
      double GetAngle(OBAtom* a, OBAtom* b, OBAtom* c);
      //! \return the stochoimetric formula (e.g., C4H6O)
      std::string  GetFormula();
      //! \return the stochoimetric formula in spaced format e.g. C 4 H 6 O 1
      std::string  GetSpacedFormula(int ones=0, const char* sp=" ");
      //! \return the heat of formation for this molecule (in kcal/mol)
      double       GetEnergy() const { return(_energy); }
      //! \return the standard molar mass given by IUPAC atomic masses (amu)
      double       GetMolWt();
      //! \return the mass given by isotopes (or most abundant isotope, if not specified)
      double	 GetExactMass();
      //! \return the total charge on this molecule (i.e., 0 = neutral, +1, -1...)
      int		 GetTotalCharge();
      //! \return the total spin on this molecule (i.e., 1 = singlet, 2 = doublet...)
      unsigned int GetTotalSpinMultiplicity();
      //! \return the dimensionality of coordinates (i.e., 0 = unknown or no coord, 2=2D, 3=3D)
      unsigned short int GetDimension() const { return _dimension; }
      double      *GetCoordinates() { return(_c); }
      //! \return the Smallest Set of Smallest Rings has been run (see OBRing class
      std::vector<OBRing*> &GetSSSR();
      //! Get the current flag for whether formal charges are set with pH correction
      bool AutomaticFormalCharge()   { return(_autoFormalCharge);  }
      //! Get the current flag for whether partial charges are auto-determined
      bool AutomaticPartialCharge()  { return(_autoPartialCharge); }
      //@}


      //! \name Data modification methods
      //@{
      void   SetTitle(const char *title);
      void   SetTitle(std::string &title);
      //! Set the stochiometric formula for this molecule
      void   SetFormula(std::string molFormula);
      //! Set the heat of formation for this molecule (in kcal/mol)
      void   SetEnergy(double energy) { _energy = energy; }
      //! Set the dimension of this molecule (i.e., 0, 1 , 2, 3)
      void   SetDimension(unsigned short int d) { _dimension = d; }
      void   SetTotalCharge(int charge);
      void   SetTotalSpinMultiplicity(unsigned int spin);
      void   SetInternalCoord(std::vector<OBInternalCoord*> int_coord)
        { _internals = int_coord; }
      //! Set the flag for determining automatic formal charges with pH (default=true)
      void SetAutomaticFormalCharge(bool val)
        { _autoFormalCharge=val;  }
      //! Set the flag for determining partial charges automatically (default=true)
      void SetAutomaticPartialCharge(bool val)
        { _autoPartialCharge=val; }

      //! Mark that aromaticity has been perceived for this molecule (see OBAromaticTyper)
      void   SetAromaticPerceived()    { SetFlag(OB_AROMATIC_MOL);    }
      //! Mark that Smallest Set of Smallest Rings has been run (see OBRing class)
      void   SetSSSRPerceived()        { SetFlag(OB_SSSR_MOL);        }
      //! Mark that rings have been perceived (see OBRing class for details)
      void   SetRingAtomsAndBondsPerceived(){SetFlag(OB_RINGFLAGS_MOL);}
      //! Mark that atom types have been perceived (see OBAtomTyper for details)
      void   SetAtomTypesPerceived()   { SetFlag(OB_ATOMTYPES_MOL);   }
      //! Mark that chains and residues have been perceived (see OBChainsParser)
      void   SetChainsPerceived()      { SetFlag(OB_CHAINS_MOL);      }
      //! Mark that chirality has been perceived
      void   SetChiralityPerceived()   { SetFlag(OB_CHIRALITY_MOL);   }
      //! Mark that partial charges have been assigned
      void   SetPartialChargesPerceived(){ SetFlag(OB_PCHARGE_MOL);   }
      void   SetHybridizationPerceived() { SetFlag(OB_HYBRID_MOL);    }
      void   SetImplicitValencePerceived(){ SetFlag(OB_IMPVAL_MOL);   }
      void   SetKekulePerceived()      { SetFlag(OB_KEKULE_MOL);      }
      void   SetClosureBondsPerceived(){ SetFlag(OB_CLOSURE_MOL);     }
      void   SetHydrogensAdded()       { SetFlag(OB_H_ADDED_MOL);     }
      void   SetCorrectedForPH()       { SetFlag(OB_PH_CORRECTED_MOL);}
      void   SetAromaticCorrected()    { SetFlag(OB_AROM_CORRECTED_MOL);}
      void   SetSpinMultiplicityAssigned(){ SetFlag(OB_TSPIN_MOL);    }
      void   SetFlags(int flags)       { _flags = flags;              }

      void   UnsetAromaticPerceived()  { _flags &= (~(OB_AROMATIC_MOL));   }
      void   UnsetPartialChargesPerceived(){ _flags &= (~(OB_PCHARGE_MOL));}
      void   UnsetImplicitValencePerceived(){_flags &= (~(OB_IMPVAL_MOL)); }
      void   UnsetFlag(int flag)       { _flags &= (~(flag));              }

      //! \name Molecule modification methods
      //@{
      // Description in transform.cpp
      virtual OBBase*    DoTransformations(const std::map<std::string,std::string>* pOptions);
      static const char* ClassDescription();
      //! Clear all information from a molecule
      bool Clear();
      //! Renumber the atoms of this molecule according to the order in the supplied vector
      void RenumberAtoms(std::vector<OBNodeBase*>&);
      //! Translate one conformer and rotate by a rotation matrix (which is returned) to the inertial frame-of-reference
      void ToInertialFrame(int conf, double *rmat);
      //! Translate all conformers to the inertial frame-of-reference
      void ToInertialFrame();
      //! Translates all conformers in the molecule by the supplied vector
      void Translate(const vector3 &v);
      //! Translates one conformer in the molecule by the supplied vector
      void Translate(const vector3 &v, int conf);
      void Rotate(const double u[3][3]);
      void Rotate(const double m[9]);
      void Rotate(const double m[9],int nconf);
      //! Translate to the center of all coordinates (for this conformer)
      void Center();
      //! Transform to standard Kekule bond structure (presumably from an aromatic form)
      bool Kekulize();
      bool PerceiveKekuleBonds();

      void NewPerceiveKekuleBonds();

      bool DeleteHydrogen(OBAtom*);
      bool DeleteHydrogens();
      bool DeleteHydrogens(OBAtom*);
      bool DeleteNonPolarHydrogens();
      bool AddHydrogens(bool polaronly=false,bool correctForPH=true);
      bool AddHydrogens(OBAtom*);
      bool AddPolarHydrogens();

      //! Deletes all atoms except for the largest contiguous fragment
      bool StripSalts();
      //! Converts the charged form of coordinate bonds, e.g.[N+]([O-])=O to N(=O)=O 
      bool ConvertDativeBonds();

      bool CorrectForPH();
      bool AssignSpinMultiplicity();
      vector3 Center(int nconf);
      //! Set the torsion defined by these atoms, rotating bonded neighbors
      void SetTorsion(OBAtom*,OBAtom*,OBAtom*,OBAtom*,double);
      //@}

      //! \name Molecule utilities and perception methods
      //@{
      //! Find Smallest Set of Smallest Rings (see OBRing class for more details)
      void FindSSSR();
      void FindRingAtomsAndBonds();
      void FindChiralCenters();
      void FindChildren(std::vector<int> &,int,int);
      void FindChildren(std::vector<OBAtom*>&,OBAtom*,OBAtom*);
      void FindLargestFragment(OBBitVec &);
      //! Sort a list of contig fragments by size from largest to smallest
      //! Each vector<int> contains the atom numbers of a contig fragment
      void ContigFragList(std::vector<std::vector<int> >&);
      //! Aligns atom a on p1 and atom b along p1->p2 vector
      void Align(OBAtom*,OBAtom*,vector3&,vector3&);
      //! Adds single bonds based on atom proximity
      void ConnectTheDots();
      //! Attempts to perceive multiple bonds based on geometries
      void PerceiveBondOrders();
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
      //! Are there non-zero coordinates in two dimensions (i.e. X and Y)?
      bool Has2D();
      //! Are there non-zero coordinates in all three dimensions (i.e. X, Y, Z)?
      bool Has3D();
      //! Are there any non-zero coordinates?
      bool HasNonZeroCoords();
      bool HasAromaticPerceived()     { return(HasFlag(OB_AROMATIC_MOL)); }
      bool HasSSSRPerceived()         { return(HasFlag(OB_SSSR_MOL));     }
      bool HasRingAtomsAndBondsPerceived(){return(HasFlag(OB_RINGFLAGS_MOL));}
      bool HasAtomTypesPerceived()    { return(HasFlag(OB_ATOMTYPES_MOL));}
      bool HasChiralityPerceived()    { return(HasFlag(OB_CHIRALITY_MOL));}
      bool HasPartialChargesPerceived() { return(HasFlag(OB_PCHARGE_MOL));}
      bool HasHybridizationPerceived() { return(HasFlag(OB_HYBRID_MOL));  }
      bool HasImplicitValencePerceived() { return(HasFlag(OB_IMPVAL_MOL));}
      bool HasKekulePerceived() { return(HasFlag(OB_KEKULE_MOL));         }
      bool HasClosureBondsPerceived() { return(HasFlag(OB_CLOSURE_MOL));  }
      bool HasChainsPerceived() { return(HasFlag(OB_CHAINS_MOL));         }
      bool HasHydrogensAdded() { return(HasFlag(OB_H_ADDED_MOL));         }
      bool HasAromaticCorrected() { return(HasFlag(OB_AROM_CORRECTED_MOL));}
      bool IsCorrectedForPH() { return(HasFlag(OB_PH_CORRECTED_MOL));     }
      bool HasSpinMultiplicityAssigned() { return(HasFlag(OB_TSPIN_MOL)); }
      //! Is this molecule chiral?
      bool IsChiral();
      //! Are there any atoms in this molecule?
      bool Empty()                       { return(_natoms == 0);          }
      //@}

      //! \name Multiple conformer member functions
      //@{
      int     NumConformers()    { return((_vconf.empty())?0:_vconf.size()); }
      void    SetConformers(std::vector<double*> &v);
      void    AddConformer(double *f)    {  _vconf.push_back(f);    }
      void    SetConformer(int i)        {  _c = _vconf[i];         }
      void    CopyConformer(double*,int);
      void    DeleteConformer(int);
      double  *GetConformer(int i)       {  return(_vconf[i]);      }
      double  *BeginConformer(std::vector<double*>::iterator&i)
        { i = _vconf.begin();
        return((i == _vconf.end()) ? NULL:*i); }
      double  *NextConformer(std::vector<double*>::iterator&i)
        { i++;
        return((i == _vconf.end()) ? NULL:*i); }
      std::vector<double*> &GetConformers() {   return(_vconf);     }
      //@}

      //! \name Iterator methods
      //@{
      //! \deprecated Use FOR_ATOMS_OF_MOL and OBMolAtomIter instead
      OBAtom *BeginAtom(std::vector<OBNodeBase*>::iterator &i);
      //! \deprecated Use FOR_ATOMS_OF_MOL and OBMolAtomIter instead
      OBAtom *NextAtom(std::vector<OBNodeBase*>::iterator &i);
      //! \deprecated Use FOR_BONDS_OF_MOL and OBMolBondIter instead
      OBBond *BeginBond(std::vector<OBEdgeBase*>::iterator &i);
      //! \deprecated Use FOR_BONDS_OF_MOL and OBMolBondIter instead
      OBBond *NextBond(std::vector<OBEdgeBase*>::iterator &i);
      //! \deprecated Use FOR_RESIDUES_OF_MOL and OBResidueIter instead
      OBResidue *BeginResidue(std::vector<OBResidue*>::iterator &i)
        {
          i = _residue.begin();
          return((i == _residue.end()) ? NULL:*i);
        }
      //! \deprecated Use FOR_RESIDUES_OF_MOL and OBResidueIter instead
      OBResidue *NextResidue(std::vector<OBResidue*>::iterator &i)
        {
          i++;
          return((i == _residue.end()) ? NULL:*i);
        }
      OBInternalCoord *BeginInternalCoord(std::vector<OBInternalCoord*>::iterator &i)
        {
          i = _internals.begin();
          return((i == _internals.end()) ? NULL:*i);
        }
      OBInternalCoord *NextInternalCoord(std::vector<OBInternalCoord*>::iterator &i)
        {
          i++;
          return((i == _internals.end()) ? NULL:*i);
        }
      //@}

    };

  //function prototypes

  OBAPI bool tokenize(std::vector<std::string>&, const char *buf, const char *delimstr=" \t\n");
  OBAPI bool tokenize(std::vector<std::string>&, std::string&, const char *delimstr=" \t\n", int limit=-1);
  //! remove leading and trailing whitespace from a string
  OBAPI std::string& Trim(std::string& txt);
  //! \deprecated -- use OBMessageHandler class instead
  OBAPI void ThrowError(char *str);
  //! \deprecated -- use OBMessageHandler class instead
  OBAPI void ThrowError(std::string &str);
  OBAPI void CartesianToInternal(std::vector<OBInternalCoord*>&,OBMol&);
  OBAPI void InternalToCartesian(std::vector<OBInternalCoord*>&,OBMol&);
  OBAPI std::string NewExtension(std::string&,char*);
  // Now handled by OBConversion class
  // OBAPI bool SetInputType(OBMol&,std::string&);
  // OBAPI bool SetOutputType(OBMol&,std::string&);

  //global definitions
  //! Global OBElementTable for element properties
  EXTERN  OBElementTable   etab;
  //! Global OBTypeTable for translating between different atom types
  //! (e.g., Sybyl <-> MM2)
  EXTERN  OBTypeTable      ttab;
  //! Global OBIsotopeTable for isotope properties
  EXTERN  OBIsotopeTable   isotab;
  //! Global OBAromaticTyper for detecting aromatic atoms and bonds
  EXTERN  OBAromaticTyper  aromtyper;
  //! Global OBAtomTyper for marking internal valence, hybridization,
  //!  and atom types (for internal and external use)
  EXTERN  OBAtomTyper      atomtyper;
  //! Global OBChainsParser for detecting macromolecular chains and residues
  EXTERN  OBChainsParser   chainsparser;
  //! Global OBMessageHandler error handler
  OBERROR extern  OBMessageHandler obErrorLog;
  //! Global OBResidueData biomolecule residue database
  EXTERN  OBResidueData    resdat;

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

} // end namespace OpenBabel

#endif // OB_MOL_H

//! \file mol.h
//! \brief Handle molecules. Declarations of OBMol, OBAtom, OBBond, OBResidue.
//!        (the main header for Open Babel)

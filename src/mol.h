/**********************************************************************
mol.h - Handle molecules.

Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (c) 2001-2004 by Geoffrey R. Hutchison
Some portions Copyright (c) 2003 by Michael Banck

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

#include "babelconfig.h"

#ifndef OB_MOL_H
#define OB_MOL_H

#include <math.h>

#include <algorithm>
#include <vector>
#include <string>

#if HAVE_IOSTREAM
#include <iostream>
#elif HAVE_IOSTREAM_H
#include <iostream.h>
#endif

#if HAVE_FSTREAM
#include <fstream>
#elif HAVE_FSTREAM_H
#include <fstream.h>
#endif

#include "base.h"
#include "data.h"
#include "chains.h"
#include "math/vector3.h"
#include "bitvec.h"
#include "ring.h"
#include "generic.h"
#include "typer.h"
#include "fileformat.h"

namespace OpenBabel {

class OBAtom;
class OBBond;
class OBMol;
class OBInternalCoord;

// Class OBResidue

// class introduction in residue.cpp
class OBResidue
{
public:

  //! Constructor
  OBResidue(void);
  //! Copy constructor
  OBResidue(const OBResidue &);
  //! Destructor
  virtual ~OBResidue(void);

  OBResidue &operator=(const OBResidue &);

  void AddAtom(OBAtom *atom); 
  void InsertAtom(OBAtom *atom);
  void RemoveAtom(OBAtom *atom);
  void Clear(void);

  void    SetName(const std::string &resname);
  void    SetNum(unsigned int resnum);
  void    SetChain(char chain);
  void    SetChainNum(unsigned int chainnum);
  void    SetIdx(unsigned int idx);

  void    SetAtomID(OBAtom *atom, const std::string &id);
  void    SetHetAtom(OBAtom *atom, bool hetatm);
  void    SetSerialNum(OBAtom *atom, unsigned sernum);

  std::string    GetName(void)			const;
  unsigned int   GetNum(void)			const;
  unsigned int	 GetNumAtoms()			const;
  char           GetChain(void)			const;
  unsigned int   GetChainNum(void)		const;
  unsigned int   GetIdx(void)			const;
  unsigned int	 GetResKey(void)		const;

  std::vector<OBAtom*> GetAtoms(void)		const;
  std::vector<OBBond*> GetBonds(bool = true)	const;

  std::string    GetAtomID(OBAtom *atom)	const;
  unsigned       GetSerialNum(OBAtom *atom)	const;

  bool           GetAminoAcidProperty(int)      const;
  bool           GetAtomProperty(OBAtom *, int) const;
  bool           GetResidueProperty(int)        const;

  bool           IsHetAtom(OBAtom *atom)	const;
  bool		 IsResidueType(int)		const;

  OBAtom *BeginAtom(std::vector<OBAtom*>::iterator &i);
  OBAtom *NextAtom(std::vector<OBAtom*>::iterator &i);

  //! \name Methods for handling generic data
  //@{
  bool                              HasData(std::string &);
  bool                              HasData(const char *);
  bool                              HasData(obDataType);
  void                              DeleteData(obDataType);
  void                              DeleteData(OBGenericData*);
  void                              DeleteData(std::vector<OBGenericData*>&);
  void                              SetData(OBGenericData *d) {_vdata.push_back(d);     }
  //! Return the number of OBGenericData items attached to this molecule.
  unsigned int                      DataSize()                {return(_vdata.size());   }
  OBGenericData                    *GetData(obDataType);
  OBGenericData                    *GetData(std::string&);
  OBGenericData                    *GetData(const char *);
  std::vector<OBGenericData*>           &GetData()                 {return(_vdata);         }
  std::vector<OBGenericData*>::iterator  BeginData()               {return(_vdata.begin()); }
  std::vector<OBGenericData*>::iterator  EndData()                 {return(_vdata.end());   }
  //@}

protected: // members

  unsigned int	 	      _idx;
  char      	              _chain;
  unsigned int		      _aakey;
  unsigned int	              _reskey;
  unsigned int		      _resnum;
  std::string                 _resname;

  std::vector<bool>           _hetatm;
  std::vector<std::string>    _atomid;
  std::vector<OBAtom*>        _atoms;
  std::vector<unsigned int>   _sernum;
  std::vector<OBGenericData*> _vdata;
};


// Class OBAtom

//ATOM Property Macros (flags)
#define OB_4RING_ATOM     (1<<1)
#define OB_3RING_ATOM     (1<<2)
#define OB_AROMATIC_ATOM  (1<<3)
#define OB_RING_ATOM      (1<<4)
#define OB_CSTEREO_ATOM   (1<<5)
#define OB_ACSTEREO_ATOM  (1<<6)
#define OB_DONOR_ATOM     (1<<7)
#define OB_ACCEPTOR_ATOM  (1<<8)
#define OB_CHIRAL_ATOM    (1<<9)
// 10-16 currently unused

// class introduction in atom.cpp
class OBAtom : public OBNodeBase
{
protected:
  char                          _ele;		//!< atomic number
  char                          _impval;	//!< implicit valence
  char                          _type[6];	//!< atomic type
  short int                     _fcharge;	//!< formal charge
  int                           _isotope;	//!< isotope (0 = most abundant)
  int                           _spinmultiplicity;// 2 for radical  1 or 3 for carbene

  //unsigned short int          _idx;		//!< index in parent (inherited)
  unsigned short int            _cidx;		//!< index into coordinate array
  unsigned short int            _hyb;		//!< hybridization
  unsigned short int            _flags;		//!< bitwise flags (e.g. aromaticity)
  double                         _pcharge;	//!< partial charge
  double                       **_c;		//!< coordinate array in double*
  vector3                       _v;		//!< coordinate vector
  OBResidue                    *_residue;	//!< parent residue (if applicable)
  //OBMol                      *_parent;        //!< parent molecule (inherited)
  //vector<OBBond*>             _bond;		//!< connections (inherited)
  std::vector<OBGenericData*>   _vdata;		//!< custom data

  int  GetFlag()                const {return(_flags);}
  void SetFlag(int flag)       {_flags |= flag;}
  bool HasFlag(int flag)       {return((_flags & flag) ? true : false);}

public:
    //! Constructor
    OBAtom();
    //! Destructor
    virtual ~OBAtom();
    //! Assignment
    OBAtom &operator=(OBAtom &);
    //! Clear all data
    void Clear();

    //! \name Methods to set atomic information
    //@{    
    void SetIdx(int idx)                     {_idx = idx;_cidx = (idx-1)*3;}
    void SetHyb(int hyb)                     {_hyb = hyb;}
    void SetAtomicNum(int atomicnum)         {_ele = (char)atomicnum;}
    void SetIsotope(unsigned int iso);
    void SetImplicitValence(int val)         {_impval = (char)val;}
    void IncrementImplicitValence()          {_impval++;}
    void DecrementImplicitValence()          {_impval--;}
    void SetFormalCharge(int fcharge)        {_fcharge = fcharge;}
    void SetSpinMultiplicity(int spin)	     {_spinmultiplicity = spin;} //CM 18 Sept 2003
    void SetType(char *type);
    void SetType(std::string &type);
    void SetPartialCharge(double pcharge)     {_pcharge = pcharge;}
    void SetVector();
    void SetVector(vector3 &v);                
    void SetVector(const double x,const double y,const double z);
    void SetResidue(OBResidue *res)          {_residue=res;}
//  void SetParent(OBMol *ptr)               {_parent=ptr;}
    void SetCoordPtr(double **c)             {_c = c;_cidx = (GetIdx()-1)*3;}
    void SetAromatic()                       {SetFlag(OB_AROMATIC_ATOM);}
    void UnsetAromatic()    		     {_flags &= (~(OB_AROMATIC_ATOM));}
    void SetClockwiseStereo()          {SetFlag(OB_CSTEREO_ATOM|OB_CHIRAL_ATOM);}
    void SetAntiClockwiseStereo()      {SetFlag(OB_ACSTEREO_ATOM|OB_CHIRAL_ATOM);}
    void UnsetStereo()                 
      {
        _flags &= ~(OB_ACSTEREO_ATOM);
        _flags &= ~(OB_CSTEREO_ATOM);
        _flags &= ~(OB_CHIRAL_ATOM);
      }
    void SetInRing()                         {SetFlag(OB_RING_ATOM);}
    void SetChiral()                         {SetFlag(OB_CHIRAL_ATOM);}
    void ClearCoordPtr()                     {_c = NULL;_cidx=0;}
    //@}

    //! \name Methods to retrieve atomic information
    //@{
    //int        GetStereo()        const {return((int)_stereo);}
    int          GetFormalCharge()  const {return(_fcharge);}
    unsigned int GetAtomicNum()     const {return((unsigned int)_ele);}
    unsigned short int GetIsotope() const {return(_isotope);}
    int          GetSpinMultiplicity() const {return(_spinmultiplicity);} //CM 18 Sept 2003   
    //! The atomic mass of this atom given by standard IUPAC average molar mass
    double	 GetAtomicMass()    const;
    //! The atomic mass of given by the isotope (default is most abundant isotope)
    double	 GetExactMass()	    const;
    unsigned int GetIdx()           const {return((int)_idx);}
    unsigned int GetCoordinateIdx() const {return((int)_cidx);}
    unsigned int GetCIdx()          const {return((int)_cidx);}
    //! The current number of explicit connections
    unsigned int GetValence()       const {return((_vbond.empty()) ? 0 : _vbond.size());}
    //! The hybridization of this atom (i.e. 1 for sp, 2 for sp2, 3 for sp3)
    unsigned int GetHyb()             const;
    //! The implicit valence of this atom type (i.e. maximum number of connections expected)
    unsigned int GetImplicitValence() const;
    //! The number of non-hydrogens connected to this atom
    unsigned int GetHvyValence()      const;
    //! The number of heteroatoms connected to an atom
    unsigned int GetHeteroValence()   const;
    char        *GetType();

    //! The x coordinate
    double      GetX()             {return(x());}
    double      x()
      {if (_c) return((*_c)[_cidx]); else return _v.x();}
    //! The y coordinate
    double      GetY()             {return(y());}
    double      y()
      {if (_c) return((*_c)[_cidx+1]); else return _v.y();}
    //! The z coordinate
    double      GetZ()             {return(z());}
    double      z()
      {if (_c) return((*_c)[_cidx+2]); else return _v.z();}
    //! Return the coordinates as a double*
    double     *GetCoordinate()    {return(&(*_c)[_cidx]);}
    //! Return the coordinates as a vector3 object
    vector3   &GetVector();
    double      GetPartialCharge();
    OBResidue *GetResidue();
    //OBMol   *GetParent()        {return((OBMol*)_parent);}
    //! Create a vector for a new bond from this atom, with length supplied by the double
    bool       GetNewBondVector(vector3 &v,double length);
    OBBond    *GetBond(OBAtom *);
    OBAtom    *GetNextAtom();
    //@}

    //! \name Iterator methods
    //@{
    std::vector<OBEdgeBase*>::iterator BeginBonds() {return(_vbond.begin());}
    std::vector<OBEdgeBase*>::iterator EndBonds()   {return(_vbond.end());}
    OBBond *BeginBond(std::vector<OBEdgeBase*>::iterator &i);
    OBBond *NextBond(std::vector<OBEdgeBase*>::iterator &i); 
    OBAtom *BeginNbrAtom(std::vector<OBEdgeBase*>::iterator &); 
    OBAtom *NextNbrAtom(std::vector<OBEdgeBase*>::iterator &); 
    //@}

    double GetDistance(int);
    double GetDistance(OBAtom*);
    double GetAngle(int, int);
    double GetAngle(OBAtom*, OBAtom*);

    //! \name Addition of residue/bond info. for an atom
    //@{
    void NewResidue()                  {if (!_residue) _residue = new OBResidue;}
    void DeleteResidue()               {if (_residue) delete _residue;}
    void AddBond(OBBond *bond)         {_vbond.push_back((OBEdgeBase*)bond);}
    void InsertBond(std::vector<OBEdgeBase*>::iterator &i, OBBond *bond)
      {_vbond.insert(i, (OBEdgeBase*)bond);}
    bool DeleteBond(OBBond*);
    //@}

    //! \name Requests for atomic property information
    //@{
    //! The number of oxygen atoms connected that only have one heavy valence
    unsigned int  CountFreeOxygens()      const;
    //! The number of hydrogens needed to fill the implicit valence of this atom
    unsigned int  ImplicitHydrogenCount() const;
    //! The number of hydrogens explicitly bound to this atom currently
    unsigned int  ExplicitHydrogenCount() const;
    //! The number of rings that contain this atom
    unsigned int  MemberOfRingCount()     const;
    //! The size of the smallest ring that contains this atom (0 if not in a ring)
    unsigned int  MemberOfRingSize()	  const;
    //! The smallest angle of bonds to this atom
    double	  SmallestBondAngle();
    //! The sum of the bond orders of the bonds to the atom (i.e. double bond = 2...)
    unsigned int  BOSum()                 const;
    //! The sum of the bond orders of bonds to the atom, considering only KDouble, KTriple bonds
    unsigned int  KBOSum()                const;
    //@}

    //! \name Builder utilities
    //@{
    //! If this is a hydrogen atom, transform into a methyl group
    bool HtoMethyl();
    //! Change the hybridization of this atom and modify the geometry accordingly
    bool SetHybAndGeom(int);
    //@}

    //! \name Property information
    //@{
    //! Is there any residue information?
    bool HasResidue()             {return(_residue != NULL);}
    bool IsHydrogen()             {return(GetAtomicNum() == 1);}
    bool IsCarbon()               {return(GetAtomicNum() == 6);}
    bool IsNitrogen()             {return(GetAtomicNum() == 7);}
    bool IsOxygen()               {return(GetAtomicNum() == 8);}
    bool IsSulfur()               {return(GetAtomicNum() == 16);}
    bool IsPhosphorus()           {return(GetAtomicNum() == 15);}
    bool IsAromatic()      const;
    bool IsInRing()        const;
    bool IsInRingSize(int) const;
    bool IsHeteroatom();
    //! Is this atom connected to the supplied OBAtom?
    bool IsConnected(OBAtom*);
    //! Is this atom related to the supplied OBAtom in a 1,3 bonding pattern?
    bool IsOneThree(OBAtom*);
    //! Is this atom related to the supplied OBAtom in a 1,4 bonding pattern?
    bool IsOneFour(OBAtom*);
    bool IsCarboxylOxygen();
    bool IsPhosphateOxygen();
    bool IsSulfateOxygen();
    bool IsNitroOxygen();
    bool IsAmideNitrogen();
    bool IsPolarHydrogen();
    bool IsNonPolarHydrogen();
    bool IsAromaticNOxide();
    bool IsChiral();
    bool IsAxial();
    bool IsClockwise()        {return(HasFlag(OB_CSTEREO_ATOM));}
    bool IsAntiClockwise()    {return(HasFlag(OB_ACSTEREO_ATOM));}
    bool HasChiralitySpecified() 
                              {return(HasFlag(OB_CSTEREO_ATOM|OB_ACSTEREO_ATOM));}
    bool HasAlphaBetaUnsat(bool includePandS=true);
    bool HasBondOfOrder(unsigned int);
    int  CountBondsOfOrder(unsigned int);
    bool HasNonSingleBond();
    bool HasSingleBond()      {return(HasBondOfOrder(1));}
    bool HasDoubleBond()      {return(HasBondOfOrder(2));}
    bool HasAromaticBond()    {return(HasBondOfOrder(5));}
    //! Determines if this atom matches the first atom in a given SMARTS pattern
    bool MatchesSMARTS(const char *);
    //@}

    //! \name Methods for handling generic data
    //@{
    bool                              HasData(std::string &);
    bool                              HasData(const char *);
    bool                              HasData(obDataType);
    void                              DeleteData(obDataType);
    void                              DeleteData(OBGenericData*);
    void                              DeleteData(std::vector<OBGenericData*>&);
    void                              SetData(OBGenericData *d) {_vdata.push_back(d);     }
    //! Return the number of OBGenericData items attached to this molecule.
    unsigned int                      DataSize()                {return(_vdata.size());   }
    OBGenericData                    *GetData(obDataType);
    OBGenericData                    *GetData(std::string&);
    OBGenericData                    *GetData(const char *);
    std::vector<OBGenericData*>           &GetData()                 {return(_vdata);         }
    std::vector<OBGenericData*>::iterator  BeginData()               {return(_vdata.begin()); }
    std::vector<OBGenericData*>::iterator  EndData()                 {return(_vdata.end());   }
    //@}
}; // class OBAtom


// Class OBBond

//BOND Property Macros (flags)
#define OB_AROMATIC_BOND  (1<<1)
#define OB_WEDGE_BOND     (1<<2)
#define OB_HASH_BOND      (1<<3)
#define OB_RING_BOND      (1<<4)
#define OB_TORUP_BOND     (1<<5)
#define OB_TORDOWN_BOND   (1<<6)
#define OB_KSINGLE_BOND   (1<<7) 
#define OB_KDOUBLE_BOND   (1<<8)
#define OB_KTRIPLE_BOND   (1<<9)
#define OB_CLOSURE_BOND   (1<<10)
// 11-16 currently unused

// class introduction in bond.cpp
class OBBond : public OBEdgeBase
{
protected:
  char                          _order; //!< Bond order (1, 2, 3, 5=aromatic)
  unsigned short int            _flags; //!< Any flags for this bond
  //OBAtom                     *_bgn;   //!< Not needed, inherited from OBEdgeBase
  //OBAtom                     *_end;   //!< Not needed, inherited from OBEdgeBase
  //OBMol                      *_parent;//!< Not needed, inherited from OBEdgeBase
  //unsigned short int          _idx;   //!< Not needed, inherited from OBEdgeBase
  std::vector<OBGenericData*>   _vdata; //!< Generic data for custom information

  bool HasFlag(int flag)       {return((_flags & flag) != 0);}
  void SetFlag(int flag)       {_flags |= flag;}

public:
   //! Constructor
   OBBond();
   //! Destructor
   virtual ~OBBond();

   //! \name Bond modification methods
   //@{
   void SetIdx(int idx)                     {_idx = idx;}
   void SetBO(int order);
   void SetBegin(OBAtom *begin)             {_bgn = begin;}
   void SetEnd(OBAtom *end)                 {_end = end;}
// void SetParent(OBMol *ptr)               {_parent=ptr;} // (inherited)
   void SetLength(OBAtom*,double);
   void Set(int,OBAtom*,OBAtom*,int,int);
   void SetKSingle();
   void SetKDouble();
   void SetKTriple();
   void SetAromatic()                       {SetFlag(OB_AROMATIC_BOND);}
   void SetUp()                             {SetFlag(OB_TORUP_BOND);}
   void SetDown()                           {SetFlag(OB_TORDOWN_BOND);}
   void SetInRing()                         {SetFlag(OB_RING_BOND);}
   void SetClosure()                        {SetFlag(OB_CLOSURE_BOND);}
   void UnsetAromatic()                     {_flags &= (~(OB_AROMATIC_BOND));}
   void UnsetKekule() {_flags &= (~(OB_KSINGLE_BOND|OB_KDOUBLE_BOND|OB_KTRIPLE_BOND));}
   //@}
   
   //! \name bond data request methods
   //@{
   unsigned int     GetBO()            const    {return((int)_order);}
   unsigned int     GetBondOrder()     const    {return((int)_order);}
   unsigned int     GetFlags()         const    {return(_flags);}
   unsigned int     GetBeginAtomIdx()  const    {return(_bgn->GetIdx());}
   unsigned int     GetEndAtomIdx()    const    {return(_end->GetIdx());}
   OBAtom *GetBeginAtom()              {return((OBAtom*)_bgn);}
   OBAtom *GetEndAtom()                {return((OBAtom*)_end);}
   OBAtom *GetNbrAtom(OBAtom *ptr)     {return((ptr != _bgn)? (OBAtom*)_bgn : (OBAtom*)_end);}
// OBMol  *GetParent()                 {return(_parent);}  // (inherited)
   double   GetEquibLength();
   double   GetLength();
   int     GetNbrAtomIdx(OBAtom *ptr) 
     {return((ptr!=_bgn)?_bgn->GetIdx():_end->GetIdx());}
   //@}

   //! \name property request methods
   //@{
   bool IsAromatic() const;
   bool IsInRing() const;
   bool IsRotor();
   bool IsAmide();
   bool IsPrimaryAmide();
   bool IsSecondaryAmide();
   bool IsEster();
   bool IsCarbonyl();
   bool IsSingle();
   bool IsDouble();
   bool IsTriple();
   bool IsKSingle();
   bool IsKDouble();
   bool IsKTriple();
   bool IsClosure();
   bool IsUp()                        {return(HasFlag(OB_TORUP_BOND));}
   bool IsDown()                      {return(HasFlag(OB_TORDOWN_BOND));}
   bool IsWedge()                     {return(HasFlag(OB_WEDGE_BOND));}
   bool IsHash()                      {return(HasFlag(OB_HASH_BOND));}
   //@}

    //! \name Methods for handling generic data
    //@{
    bool                              HasData(std::string &);
    bool                              HasData(const char *);
    bool                              HasData(obDataType);
    void                              DeleteData(obDataType);
    void                              DeleteData(OBGenericData*);
    void                              DeleteData(std::vector<OBGenericData*>&);
    void                              SetData(OBGenericData *d) {_vdata.push_back(d);     }
    //! Return the number of OBGenericData items attached to this molecule.
    unsigned int                      DataSize()                {return(_vdata.size());   }
    OBGenericData                    *GetData(obDataType);
    OBGenericData                    *GetData(std::string&);
    OBGenericData                    *GetData(const char *);
    std::vector<OBGenericData*>           &GetData()                 {return(_vdata);         }
    std::vector<OBGenericData*>::iterator  BeginData()               {return(_vdata.begin()); }
    std::vector<OBGenericData*>::iterator  EndData()                 {return(_vdata.end());   }
    //@}
}; // class OBBond


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
#define OB_CURRENT_CONFORMER	 -1

// class introduction in mol.cpp
class OBMol : public OBGraphBase
{
protected:
  int                           _flags;	//!< bitfield of flags
  bool                          _autoPartialCharge; //!< Assign partial charges automatically
  bool                          _autoFormalCharge; //!< Assign formal charges automatically
  bool                          _compressed; //!< Is the molecule currently compressed?
  io_type                       _itype;	//!< Input I/O Type
  io_type                       _otype;	//!< Output I/O Type
  std::string                   _title;	//!< Molecule title
  //vector<OBAtom*>             _atom;	//!< not needed (inherited)
  //vector<OBBond*>             _bond;	//!< not needed (inherited)
  std::vector<OBResidue*>       _residue;//!< Residue information (if applicable)
  std::vector<OBInternalCoord*> _internals;//!< Internal Coordinates (if applicable)
  std::vector<OBGenericData*>   _vdata;	//!< Custom data -- see OBGenericData class for more
  double                        _energy;//!< Molecular heat of formation (if applicable)
  int				_totalCharge; //!< Total charge on the molecule
  unsigned int                  _totalSpin; //!< Total spin on the molecule (if not specified, assumes lowest possible spin)
  double                        *_c;	//!< coordinate array
  std::vector<double*>           _vconf;//!< vector of conformers
  unsigned short int            _natoms;//!< Number of atoms
  unsigned short int            _nbonds;//!< Number of bonds
  unsigned short int            _mod;	//!< Number of nested calls to BeginModify()
  unsigned short int            _access;//!< Number of nested calls to BeginAccess()

  bool  HasFlag(int flag)  {return((_flags & flag) ? true : false);}
  void  SetFlag(int flag)  {_flags |= flag;}

public:

    //! \name Initialization and data (re)size methods
    //@{
    //! Constructor (input format) (output format)
    OBMol(io_type=UNDEFINED,io_type=UNDEFINED);
    //! Copy constructor
    OBMol(const OBMol &);
    //! Destructor
    virtual ~OBMol();
    //! Assignment
    OBMol &operator=(const OBMol &mol);
    OBMol &operator+=(const OBMol &mol);
    void ReserveAtoms(int natoms) {if (natoms && _mod) _vatom.reserve(natoms);}
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
    int GetMod() {return(_mod);}
    void IncrementMod() {_mod++;}
    void DecrementMod() {_mod--;}
    //@}

    //! \name Compression methods (via OBCompressData)
    //@{
    //! Save memory by rolling into the OEBinary format as a buffer
    virtual bool Compress(void);
    //! Roll data out from the obCompressData buffer
    virtual bool UnCompress(void);
    //! Call when accessing a compressed molecule to uncompress if needed.
    virtual void BeginAccess(void);
    //! Call when finished accessing a compressed molecule to recompress.
    virtual void EndAccess(void);
    //@}

    //! \name Generic data handling methods (via OBGenericData)
    //@{
     //! Returns true if the generic attribute/value pair exists
    bool                              HasData(std::string &);
     //! Returns true if the generic attribute/value pair exists
    bool                              HasData(const char *);
     //! Returns true if the generic attribute/value pair exists
    bool                              HasData(obDataType);
    void                              DeleteData(obDataType);
    void                              DeleteData(OBGenericData*);
    void                              DeleteData(std::vector<OBGenericData*>&);
    void                              SetData(OBGenericData *d) {_vdata.push_back(d); }
    //! Return the number of OBGenericData items attached to this molecule.
    unsigned int                      DataSize()                {return(_vdata.size()); }
    OBGenericData                    *GetData(obDataType);
    OBGenericData                    *GetData(std::string&);
    OBGenericData                    *GetData(const char *);
    std::vector<OBGenericData*>           &GetData()                 {return(_vdata); }
    std::vector<OBGenericData*>::iterator  BeginData()               {return(_vdata.begin()); }
    std::vector<OBGenericData*>::iterator  EndData()                 {return(_vdata.end()); }
    //@}

    //! \name Data retrieval methods
    //@{
    int          GetFlags()                           {return(_flags);}
    const char  *GetTitle() const                     {return(_title.c_str());}
    io_type      GetInputType() const                 {return(_itype);}
    io_type      GetOutputType() const                {return(_otype);}
    //! The number of atoms (i.e. OBAtom children)
    unsigned int NumAtoms() const                     {return(_natoms);}
    //! The number of bonds (i.e. OBBond children)
    unsigned int NumBonds() const                     {return(_nbonds);}
    //! The number of non-hydrogen atoms
    unsigned int NumHvyAtoms();
    //! The number of residues (i.e. OBResidue substituents)
    unsigned int NumResidues() const                  {return(_residue.size());}
    //! The number of rotatble bonds
    unsigned int NumRotors();
    OBAtom      *GetAtom(int);
    OBAtom      *GetFirstAtom();
    OBBond      *GetBond(int); 
    OBBond      *GetBond(int, int);
    OBBond      *GetBond(OBAtom*,OBAtom*);
    OBResidue   *GetResidue(int);
    std::vector<OBInternalCoord*> GetInternalCoord();
    double       GetTorsion(int,int,int,int);
    double       GetTorsion(OBAtom*,OBAtom*,OBAtom*,OBAtom*);
    //! Heat of formation for this molecule (in kcal/mol)
    double       GetEnergy() const                   {return(_energy);}
    //! Standard molar mass given by IUPAC atomic masses (amu)
    double       GetMolWt();
    //! Mass given by isotopes (or most abundant isotope, if not specified)
    double	 GetExactMass();
    //! Total charge on this molecule (i.e., 0 = neutral, +1, -1...)
    int		 GetTotalCharge();
    //! Total spin on this molecule (i.e., 1 = singlet, 2 = doublet...)
    unsigned int GetTotalSpinMultiplicity();
    double      *GetCoordinates()                    {return(_c);}
    std::vector<OBRing*> &GetSSSR();
    bool         IsCompressed() const                {return _compressed;}
    //! Get the current flag for whether formal charges are set with pH correction
    bool AutomaticFormalCharge()             {return(_autoFormalCharge);}
    //! Get the current flag for whether partial charges are auto-determined
    bool AutomaticPartialCharge()            {return(_autoPartialCharge);}
    //@}


    //! \name Data modification methods
    //@{
    void   SetTitle(const char *title)     {_title = title;}
    void   SetTitle(std::string &title)    {_title = title;}
    //! Set the heat of formation for this molecule (in kcal/mol)
    void   SetEnergy(double energy)        {_energy = energy;}
    void   SetTotalCharge(int charge);
    void   SetTotalSpinMultiplicity(unsigned int spin);
    void   SetInputType(io_type type)      {_itype = type;}
    void   SetOutputType(io_type type)     {_otype = type;}
    void   SetInternalCoord(std::vector<OBInternalCoord*> int_coord) {_internals = int_coord;}
    //! Set the flag for determining automatic formal charges with pH
    void SetAutomaticFormalCharge(bool val)     {_autoFormalCharge=val;}
    //! Set the flag for determining partial charges automatically
    void SetAutomaticPartialCharge(bool val)    {_autoPartialCharge=val;}

    void   SetAromaticPerceived()          {SetFlag(OB_AROMATIC_MOL);}
    void   SetSSSRPerceived()              {SetFlag(OB_SSSR_MOL);}
    void   SetRingAtomsAndBondsPerceived() {SetFlag(OB_RINGFLAGS_MOL);}
    void   SetAtomTypesPerceived()         {SetFlag(OB_ATOMTYPES_MOL);}
    void   SetChainsPerceived()            {SetFlag(OB_CHAINS_MOL);}
    void   SetChiralityPerceived()         {SetFlag(OB_CHIRALITY_MOL);}
    void   SetPartialChargesPerceived()    {SetFlag(OB_PCHARGE_MOL);}
    void   SetHybridizationPerceived()     {SetFlag(OB_HYBRID_MOL);}
    void   SetImplicitValencePerceived()   {SetFlag(OB_IMPVAL_MOL);}
    void   SetKekulePerceived()            {SetFlag(OB_KEKULE_MOL);}
    void   SetClosureBondsPerceived()      {SetFlag(OB_CLOSURE_MOL);}
    void   SetHydrogensAdded()             {SetFlag(OB_H_ADDED_MOL);}
    void   SetCorrectedForPH()             {SetFlag(OB_PH_CORRECTED_MOL);}
    void   SetAromaticCorrected()          {SetFlag(OB_AROM_CORRECTED_MOL);}
    void   UnsetAromaticPerceived()        {_flags &= (~(OB_AROMATIC_MOL));}
    void   UnsetPartialChargesPerceived()  {_flags &= (~(OB_PCHARGE_MOL));}
    void   UnsetImplicitValencePerceived() {_flags &= (~(OB_IMPVAL_MOL));}
    void   UnsetFlag(int flag)		   {_flags &= (~(flag));}
    void   SetFlags(int flags)             {_flags = flags;}

    //! \name Molecule modification methods
    //@{
    //! Clear all information from a molecule
    bool Clear();
    void RenumberAtoms(std::vector<OBNodeBase*>&);
    void ToInertialFrame(int,double*);
    void ToInertialFrame();
    //! Translates all conformers in the molecule
    void Translate(const vector3 &v);
    //! Translates one conformer in the molecule
    void Translate(const vector3 &v, int conf);
    void Rotate(const double u[3][3]);
    void Rotate(const double m[9]);
    void Rotate(const double m[9],int nconf);
    // Translate to the center of all coordinates (for this conformer)
    void Center();
    //! Transform to standard Kekule bond structure (presumably from an aromatic form)
    bool Kekulize();
    bool PerceiveKekuleBonds();
    bool DeleteHydrogen(OBAtom*);
    bool DeleteHydrogens();
    bool DeleteHydrogens(OBAtom*);
    bool DeleteNonPolarHydrogens();
    bool AddHydrogens(bool polaronly=false,bool correctForPH=true);
    bool AddHydrogens(OBAtom*);
    bool AddPolarHydrogens();
    bool StripSalts();
    bool CorrectForPH();
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
//  void ConnectTheDotsSort();
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
    bool HasNonZeroCoords();
    bool HasAromaticPerceived()          {return(HasFlag(OB_AROMATIC_MOL));       }
    bool HasSSSRPerceived()              {return(HasFlag(OB_SSSR_MOL));           }
    bool HasRingAtomsAndBondsPerceived() {return(HasFlag(OB_RINGFLAGS_MOL));      }
    bool HasAtomTypesPerceived()         {return(HasFlag(OB_ATOMTYPES_MOL));      }
    bool HasChiralityPerceived()         {return(HasFlag(OB_CHIRALITY_MOL));      }
    bool HasPartialChargesPerceived()    {return(HasFlag(OB_PCHARGE_MOL));        }
    bool HasHybridizationPerceived()     {return(HasFlag(OB_HYBRID_MOL));         }
    bool HasImplicitValencePerceived()   {return(HasFlag(OB_IMPVAL_MOL));         }
    bool HasKekulePerceived()            {return(HasFlag(OB_KEKULE_MOL));         }
    bool HasClosureBondsPerceived()      {return(HasFlag(OB_CLOSURE_MOL));        }
    bool HasChainsPerceived()            {return(HasFlag(OB_CHAINS_MOL));         }
    bool HasHydrogensAdded()             {return(HasFlag(OB_H_ADDED_MOL));        }
    bool HasAromaticCorrected()          {return(HasFlag(OB_AROM_CORRECTED_MOL)); }
    bool IsCorrectedForPH()              {return(HasFlag(OB_PH_CORRECTED_MOL));   }
    bool IsChiral();
    //! Are there any atoms in this molecule?
    bool Empty() {return(_natoms == 0);}
    //@}

    //! \name Iterator methods
    //@{
    OBAtom *BeginAtom(std::vector<OBNodeBase*>::iterator &i);
    OBAtom *NextAtom(std::vector<OBNodeBase*>::iterator &i);
    OBBond *BeginBond(std::vector<OBEdgeBase*>::iterator &i); 
    OBBond *NextBond(std::vector<OBEdgeBase*>::iterator &i); 
    OBResidue *BeginResidue(std::vector<OBResidue*>::iterator &i)
      {i = _residue.begin();return((i == _residue.end()) ? NULL:*i);}
    OBResidue *NextResidue(std::vector<OBResidue*>::iterator &i)
      {i++;return((i == _residue.end()) ? NULL:*i);}
    OBInternalCoord *BeginInternalCoord(std::vector<OBInternalCoord*>::iterator &i)
      {i = _internals.begin();return((i == _internals.end()) ? NULL:*i);}
    OBInternalCoord *NextInternalCoord(std::vector<OBInternalCoord*>::iterator &i)
      {i++;return((i == _internals.end()) ? NULL:*i);}
    //@}

    //! \name Multiple conformer member functions
    //@{
    int     NumConformers()         {return((_vconf.empty())?0:_vconf.size());}
    void    SetConformers(std::vector<double*> &v);
    void    AddConformer(double *f)                 {_vconf.push_back(f);}
    void    SetConformer(int i)                     {_c = _vconf[i];}
    void    CopyConformer(double*,int);
    void    DeleteConformer(int);
    double  *GetConformer(int i)                    {return(_vconf[i]);}
    double  *BeginConformer(std::vector<double*>::iterator&);
    double  *NextConformer(std::vector<double*>::iterator&);
    std::vector<double*> &GetConformers()           {return(_vconf);}
    //@}

    //! \name Misc bond functions
    //@{
    void AssignResidueBonds(OBBitVec &);
    //void SortBonds() {sort(_vbond.begin(),_vbond.end(),CompareBonds);}
    //@}

    //! \name Convenience functions for I/O
    //@{
    friend std::ostream&       operator<< ( std::ostream&, OBMol& ) ;
    friend std::istream&       operator>> ( std::istream&, OBMol& ) ;
    //@}
};

//! \brief Used to transform from z-matrix to cartesian coordinates.
class OBInternalCoord 
{
public:
  //class members
  OBAtom *_a,*_b,*_c;
  double   _dst,_ang,_tor;
  //! Constructor
  OBInternalCoord(OBAtom *a=(OBAtom*)NULL,
                  OBAtom *b=(OBAtom*)NULL,
                  OBAtom *c=(OBAtom*)NULL)
    {
      _a = a; _b = b; _c = c;
      _dst = _ang = _tor = 0.0;
    }
};

//function prototypes
// (prototypes for Read/Write are in fileformat.h

bool tokenize(std::vector<std::string>&, const char *buf, const char *delimstr=" \t\n");
bool tokenize(std::vector<std::string> &,std::string&, const char*,int limit=-1);
void ThrowError(char *str);
void ThrowError(std::string &str);
void CartesianToInternal(std::vector<OBInternalCoord*>&,OBMol&);
void InternalToCartesian(std::vector<OBInternalCoord*>&,OBMol&);
std::string NewExtension(std::string&,char*);
bool SetInputType(OBMol&,std::string&);
bool SetOutputType(OBMol&,std::string&);

//global definitions
extern  OBExtensionTable extab;
extern  OBElementTable   etab;
extern  OBTypeTable      ttab;
extern  OBIsotopeTable   isotab;
extern  OBChainsParser   chainsparser;

//Utility Macros

#ifndef BUFF_SIZE
#define BUFF_SIZE 1024
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


#if !HAVE_RINT
inline double rint(double x) { return ( (x < 0.0) ? ceil(x-0.5) : floor(x+0.5));}
#endif

#if !HAVE_SNPRINTF
extern "C" int snprintf( char *, size_t, const char *, /* args */ ...);
#endif

#if !HAVE_STRNCASECMP
extern "C" int strncasecmp(const char *s1, const char *s2, size_t n);
#endif

#ifndef __KCC
extern "C"{
void  get_rmat(double*,double*,double*,int);
void  ob_make_rmat(double mat[3][3],double rmat[9]);
void  qtrfit (double *r,double *f,int size,double u[3][3]);
double superimpose(double*,double*,int);
}
#else
void get_rmat(double*,double*,double*,int);
void ob_make_rmat(double mat[3][3],double rmat[9]);
void qtrfit (double *r,double *f,int size,double u[3][3]);
double superimpose(double*,double*,int);
#endif // __KCC


#define obAssert(__b__) \
    if (!(__b__)) {   \
        cerr << "Assert at File " << __FILE__ << " Line " << __LINE__ << endl; \
        int *p= NULL; *p= 10; \
        exit(-1); \
    }

} // end namespace OpenBabel

#endif // OB_MOL_H

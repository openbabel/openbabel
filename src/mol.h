/**********************************************************************
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (c) 2001-2002 by Geoffrey R. Hutchison

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

#include <math.h>

#include <algorithm>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include "base.h"
#include "data.h"
#include "chains.h"
#include "Vector.h"
#include "bitvec.h"
#include "ring.h"
#include "generic.h"
#include "typer.h"
#include "fileformat.h"

namespace OpenBabel {

class OBAtom;
class OBBond;
class OBMol;

namespace OBAminoAcidProperty {
static const unsigned int ACIDIC      =  0;
static const unsigned int ACYCLIC     =  1;
static const unsigned int ALIPHATIC   =  2;
static const unsigned int AROMATIC    =  3;
static const unsigned int BASIC       =  4;
static const unsigned int BURIED      =  5;
static const unsigned int CHARGED     =  6;
static const unsigned int CYCLIC      =  7;
static const unsigned int HYDROPHOBIC =  8;
static const unsigned int LARGE       =  9;
static const unsigned int MEDIUM      = 10;
static const unsigned int NEGATIVE    = 11;
static const unsigned int NEUTRAL     = 12;
static const unsigned int POLAR       = 13;
static const unsigned int POSITIVE    = 14;
static const unsigned int SMALL       = 15;
static const unsigned int SURFACE     = 16;
}

namespace OBResidueAtomProperty {
static const unsigned int ALPHA_CARBON     = 0;
static const unsigned int AMINO_BACKBONE   = 1;
static const unsigned int BACKBONE         = 2;
static const unsigned int CYSTEINE_SULPHUR = 3;
static const unsigned int LIGAND           = 4;
static const unsigned int NUCLEIC_BACKBONE = 5;
static const unsigned int SHAPELY_BACKBONE = 6;
static const unsigned int SHAPELY_SPECIAL  = 7;
static const unsigned int SIDECHAIN        = 8;
static const unsigned int SUGAR_PHOSPHATE  = 9;
}
/// Residue names
namespace OBResidueIndex {
static const unsigned int ALA   =  0;
static const unsigned int GLY   =  1;
static const unsigned int LEU   =  2;
static const unsigned int SER   =  3;
static const unsigned int VAL   =  4;
static const unsigned int THR   =  5;
static const unsigned int LYS   =  6;
static const unsigned int ASP   =  7;
static const unsigned int ILE   =  8;
static const unsigned int ASN   =  9;
static const unsigned int GLU   = 10;
static const unsigned int PRO   = 11;
static const unsigned int ARG   = 12;
static const unsigned int PHE   = 13;
static const unsigned int GLN   = 14;
static const unsigned int TYR   = 15;
static const unsigned int HIS   = 16;
static const unsigned int CYS   = 17;
static const unsigned int MET   = 18;
static const unsigned int TRP   = 19;
static const unsigned int ASX   = 20;
static const unsigned int GLX   = 21;
static const unsigned int PCA   = 22;
static const unsigned int HYP   = 23;
static const unsigned int A     = 24;
static const unsigned int C     = 25;
static const unsigned int G     = 26;
static const unsigned int T     = 27;
static const unsigned int U     = 28;
static const unsigned int UPLUS = 29;
static const unsigned int I     = 30;
//static const unsigned int _1MA  = 31;
//static const unsigned int _5MC  = 32;
static const unsigned int OMC   = 33;
//static const unsigned int _1MG  = 34;
//static const unsigned int _2MG  = 35;
static const unsigned int M2G   = 36;
//static const unsigned int _7MG  = 37;
static const unsigned int OMG   = 38;
static const unsigned int YG    = 39;
static const unsigned int H2U   = 40;
//static const unsigned int _5MU  = 41;
static const unsigned int PSU   = 42;
static const unsigned int UNK   = 43;
static const unsigned int ACE   = 44;
static const unsigned int FOR   = 45;
static const unsigned int HOH   = 46;
static const unsigned int DOD   = 47;
static const unsigned int SO4   = 48;
static const unsigned int PO4   = 49;
static const unsigned int NAD   = 50;
static const unsigned int COA   = 51;
static const unsigned int NAP   = 52;
static const unsigned int NDP   = 53;
}
/// Residue types.
namespace OBResidueProperty {
static const unsigned int AMINO        = 0;
static const unsigned int AMINO_NUCLEO = 1;
static const unsigned int COENZYME     = 2;
static const unsigned int ION          = 3;
static const unsigned int NUCLEO       = 4;
static const unsigned int PROTEIN      = 5;
static const unsigned int PURINE       = 6;
static const unsigned int PYRIMIDINE   = 7;
static const unsigned int SOLVENT      = 8;
static const unsigned int WATER        = 9;
}

/** Residue information.
*
* The residue information is drawn from PDB or MOL2 files,
* and are stored in the OBResidue class. OBResidues are stored inside the 
* OBAtom class. The residue information for an atom can be requested in 
* the following way:
\code
 OBAtom *atom;
 OBResidue *r;
 atom = mol.GetAtom(1);
 r = atom->GetResidue();
\endcode
*
*/
class OBResidue
{
public:

  OBResidue(void);
  OBResidue(const OBResidue &);
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

  vector<OBAtom*> GetAtoms(void)		const;
  vector<OBBond*> GetBonds(bool = true)		const;

  std::string    GetAtomID(OBAtom *atom)	const;
  unsigned       GetSerialNum(OBAtom *atom)	const;

  bool           GetAminoAcidProperty(int)      const;
  bool           GetAtomProperty(OBAtom *, int) const;
  bool           GetResidueProperty(int)        const;

  bool           IsHetAtom(OBAtom *atom)	const;
  bool		 IsResidueType(int)		const;

  OBAtom *BeginAtom(std::vector<OBAtom*>::iterator &i);
  OBAtom *NextAtom(std::vector<OBAtom*>::iterator &i);

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

};

//ATOM Property Macros
#define OB_4RING_ATOM     (1<<1)
#define OB_3RING_ATOM     (1<<2)
#define OB_AROMATIC_ATOM  (1<<3)
#define OB_RING_ATOM      (1<<4)
#define OB_CSTEREO_ATOM   (1<<5)
#define OB_ACSTEREO_ATOM  (1<<6)
#define OB_DONOR_ATOM     (1<<7)
#define OB_ACCEPTOR_ATOM  (1<<8)
#define OB_CHIRAL_ATOM    (1<<9)
/** Atom class
*
* To understand the OBAtom class it is important to state a key
* decision on which the design was based. In OBabel the atom class
* existed, but it was only a data container. All data access and
* modification of atoms was done through the molecule. The result
* was a molecule class that was very large an unwieldy. So the OBAtom
* class was made smarter, and out many of the atom-specific routines
* were separated into the OBAtom thereby decentralizing and
* shrinking the OBMol class. As a result the OBAtom class not only
* holds data, but facilitates extraction of data perceived from both
* the atom and the molecule.
*
* A number of data extraction methods perform what is called
* `Lazy Evaluation,' which is essentially on-the-fly evaluation.
* For example, when an atom is queried as to whether it is cyclic
* or what it's hybridization state is the information is perceived
* automatically. The perception of a particular trait is actually
* performed on the entire molecule the first time it is requested of
* an atom or bond, and stored for subsequent requests for the same
* trait of additional atoms or bonds.The OBAtom class is similar to
* OBMol and the whole of Open Babel in that data access and modification
* is done through Get and Set methods.
*
* The following code demonstrates how to print out the atom numbers,
* element numbers, and coordinates of a molecule:
\code
   OBMol mol(SDF,SDF);
   OBFileFormat::ReadMolecule(cin, mol);
   OBAtom *atom;
   vector<OBNodeBase*>::iterator i;
   for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
   {
       cout << atom->GetIdx() << ` `;
       cout << atom->GetAtomicNum() << ` `;
       cout << atom->GetVector() << endl;
   }
\endcode
* A number of the property member functions indicate that atoms
* have some knowlege of their covalently attached neighbor atoms.
* Bonding information is partly redundant within a molecule in
* that an OBMol has a complete list of bonds in a molecule, and
* an OBAtom has a list bonds of which it is a member. The following
* code demonstrates how an OBAtom uses its bond information to loop
* over atoms attached to itself:
\code
   OBMol mol(SDF,SDF);
   OBFileFormat::ReadMolecule(cin, mol);
   OBAtom *atom,*nbr;
   vector<OBEdgeBase*>::iterator i;
   atom = mol.GetAtom(1);
   for (nbr = atom->BeginNbrAtom(i);nbr;nbr = atom->NextNbrAtom(i))
   {
   cout << "atom #" << atom->GetIdx() << " is attached to atom #" << nbr->GetIdx() << endl;
   }
\endcode
* should produce an output like
\code
   atom #1 is attached to atom #2
\endcode
*/
class OBAtom : public OBNodeBase
{
protected:
  char                          _ele;
  char                          _impval; //implicit valence
  char                          _type[6];
  short int                     _fcharge;
//  unsigned short int            _idx;
  unsigned short int            _cidx; //index into coordinate array
  unsigned short int            _hyb;
  unsigned short int            _flags;
  float                         _pcharge;
  float                       **_c;
  Vector                        _v;
  OBResidue                    *_residue;
//  OBMol                        *_parent;
  //vector<OBBond*>               _bond;

    int  GetFlag()          const {return(_flags);}
    void SetFlag(int flag) {_flags |= flag;}
    bool HasFlag(int flag) {return((_flags & flag) ? true : false);}

public:
    OBAtom();
    virtual ~OBAtom();
    OBAtom &operator=(OBAtom &);
    void Clear();
    
    //***methods to set atomic information***
    void SetIdx(int idx)                     {_idx = idx;_cidx = (idx-1)*3;}
    void SetHyb(int hyb)                     {_hyb = hyb;}
    void SetAtomicNum(int atomicnum)         {_ele = (char)atomicnum;}
    void SetImplicitValence(int val)         {_impval = (char)val;}
    void IncrementImplicitValence()          {_impval++;}
    void DecrementImplicitValence()          {_impval--;}
    void SetFormalCharge(int fcharge)        {_fcharge = fcharge;}
    void SetType(char *type)                 {strcpy(_type,type);}
    void SetType(std::string &type)               {strcpy(_type,type.c_str());}
    void SetPartialCharge(float pcharge)     {_pcharge = pcharge;}
    void SetVector();
    void SetVector(Vector &v);                
    void SetVector(const float,const float,const float);
    void SetResidue(OBResidue *res)          {_residue=res;}
//    void SetParent(OBMol *ptr)               {_parent=ptr;}
    void SetCoordPtr(float **c)              {_c = c;_cidx = (GetIdx()-1)*3;}
    void SetAromatic()                       {SetFlag(OB_AROMATIC_ATOM);}
    void UnsetAromatic()              {_flags &= (~(OB_AROMATIC_ATOM));}
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

    //***methods to retrieve atomic information***
    //int        GetStereo()        const {return((int)_stereo);}
    int          GetFormalCharge()  const {return(_fcharge);}
    unsigned int GetAtomicNum()   const {return((unsigned int)_ele);}
    unsigned int GetIdx()         const {return((int)_idx);}
    unsigned int GetCoordinateIdx() const {return((int)_cidx);}
    unsigned int GetCIdx()          const {return((int)_cidx);}
    unsigned int GetValence()       const {return((_vbond.empty()) ? 0 : _vbond.size());}
    unsigned int GetHyb()             const;
    unsigned int GetImplicitValence() const;
    unsigned int GetHvyValence()      const;
    unsigned int GetHeteroValence()   const;
    char      *GetType();

    float      GetX()             {return(x());}
    float      GetY()             {return(y());}
    float      GetZ()             {return(z());}
      
    float      x()
      {if (_c) return((*_c)[_cidx]); else return _v.x();}
    float      y()
      {if (_c) return((*_c)[_cidx+1]); else return _v.y();}
    float      z()
      {if (_c) return((*_c)[_cidx+2]); else return _v.z();}
    float     *GetCoordinate()    {return(&(*_c)[_cidx]);}
    float      GetPartialCharge();
    Vector     &GetVector();
    OBResidue *GetResidue();
//    OBMol     *GetParent()        {return((OBMol*)_parent);}
    bool       GetNewBondVector(Vector&,float);
    OBBond    *GetBond(OBAtom *);
    OBAtom    *GetNextAtom();

    //***iterator methods***
    std::vector<OBEdgeBase*>::iterator BeginBonds() {return(_vbond.begin());}
    std::vector<OBEdgeBase*>::iterator EndBonds()   {return(_vbond.end());}
    OBBond *BeginBond(std::vector<OBEdgeBase*>::iterator &i);
    OBBond *NextBond(std::vector<OBEdgeBase*>::iterator &i); 
    OBAtom *BeginNbrAtom(std::vector<OBEdgeBase*>::iterator &); 
    OBAtom *NextNbrAtom(std::vector<OBEdgeBase*>::iterator &); 

    //***addition of residue/bond info. for an atom***
    void NewResidue()                                          {if (!_residue) _residue = new OBResidue;}
    void DeleteResidue()                                       {if (_residue) delete _residue;}
    void AddBond(OBBond *bond)                                 {_vbond.push_back((OBEdgeBase*)bond);}
    void InsertBond(std::vector<OBEdgeBase*>::iterator &i, OBBond *bond)
      {_vbond.insert(i, (OBEdgeBase*)bond);}
    bool DeleteBond(OBBond*);

    //***requests for atomic property information***
    unsigned int  CountFreeOxygens()      const;
    unsigned int  ImplicitHydrogenCount() const;
    unsigned int  ExplicitHydrogenCount() const;
    unsigned int  MemberOfRingCount()     const;
    unsigned int  BOSum()                 const; //sum of the bond orders of the bonds to the atom
    unsigned int  KBOSum()                const;

    //builder utilities
    bool HtoMethyl();
    bool SetHybAndGeom(int);

    //property information
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
    bool IsConnected(OBAtom*);
    bool IsOneThree(OBAtom*);
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
    bool HasAromaticBond()    {return(HasBondOfOrder(4));}
};

//BOND Property Macros
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
/** Bond class
* 
*
*/
class OBBond : public OBEdgeBase
{
protected:
  char                    _order;
  unsigned short int      _flags;
//  OBAtom                 *_bgn;
//  OBAtom                 *_end;
//  OBMol                  *_parent;
//  unsigned short int      _idx;

  bool HasFlag(int flag) {return((_flags & flag) != 0);}
  void SetFlag(int flag) {_flags |= flag;}

public:
   OBBond();
   virtual ~OBBond() {}
   //***bond modification methods***
   void SetIdx(int idx)                     {_idx = idx;}
   void SetBO(int order);
   void SetBegin(OBAtom *begin)             {_bgn = begin;}
   void SetEnd(OBAtom *end)                 {_end = end;}
//   void SetParent(OBMol *ptr)               {_parent=ptr;}
   void SetLength(OBAtom*,float);
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

   
   //***bond data request methods***
   unsigned int     GetBO()            const    {return((int)_order);}
   unsigned int     GetBondOrder()     const    {return((int)_order);}
   unsigned int     GetFlags()         const    {return(_flags);}
   unsigned int     GetBeginAtomIdx()  const    {return(_bgn->GetIdx());}
   unsigned int     GetEndAtomIdx()    const    {return(_end->GetIdx());}
   OBAtom *GetBeginAtom()              {return((OBAtom*)_bgn);}
   OBAtom *GetEndAtom()                {return((OBAtom*)_end);}
   OBAtom *GetNbrAtom(OBAtom *ptr)     {return((ptr != _bgn)? (OBAtom*)_bgn : (OBAtom*)_end);}
//   OBMol  *GetParent()                 {return(_parent);}
   float   GetEquibLength();
   float   GetLength();
   int     GetNbrAtomIdx(OBAtom *ptr) 
     {return((ptr!=_bgn)?_bgn->GetIdx():_end->GetIdx());}

   //***property request methods***
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
   bool IsKSingle();
   bool IsKDouble();
   bool IsKTriple();
   bool IsClosure();
   bool IsUp()                        {return(HasFlag(OB_TORUP_BOND));}
   bool IsDown()                      {return(HasFlag(OB_TORDOWN_BOND));}
   bool IsWedge()                     {return(HasFlag(OB_WEDGE_BOND));}
   bool IsHash()                      {return(HasFlag(OB_HASH_BOND));}
};

//prototype necessary for OBMol
//bool CompareBonds(const OBBond &*,const OBBond &*);


//MOL Property Macros
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
#define OB_CURRENT_CONFORMER -1
/** Molecule Class
*
*
* The most important class in Open Babel is OBMol, or the molecule class.
* The OBMol class is designed to store all the basic information
* associated with a molecule, to make manipulations on the connection
* table of a molecule facile, and to provide member functions which
* automatically perceive information about a molecule. A guided tour
* of the OBMol class is a good place to start.
*
* An OBMol class can be declared in either of the following ways:
\code
   OBMol mol;
   //or
   OBMol mol(SDF,MOL2);
\endcode
* The second declaration type sets the input and output formats for a molecule.
* For example:
\code
   #include <iostream.h>
   #include "mol.h"
   int main(int argc,char **argv)
   {
   OBMol mol(SDF,MOL2);
   OBFileFormat::ReadMolecule(cin, mol);
   OBFileFormat::WriteMolecule(cout, mol);
   return(1);
   }
\endcode
*
* will read in a molecule in SD file format from stdin 
* (or the C++ equivalent cin) and write a MOL2 format file out
* to standard out. Additionally, The input and output formats can
* be altered after declaring an OBMol by the member functions
* OBMol::SetInputType(enum io_type type) and
* OBMol::SetOutputType(enum io_type type),
* where the current values of enum io_type are
* \code {
*              ALCHEMY, BALLSTICK, BGF, BIOSYM, BMIN, BOX, CACAO,
*              CACAOINT, CACHE, CADPAC, CCC, CDX, CHARMM, CHEM3D1,
*              CHEM3D2, CHEMDRAW, CIF, CML, CSR, CSSR, DELPDB, DMOL, DOCK,
*              FDAT, FEATURE, FH, FIX, FRACT, GAMESSIN, GAMESSOUT,
*              GAUSSIAN92, GAUSSIAN94, GAUSSIANCART, GAUSSIANZMAT,
*              GHEMICAL, GROMOS96A, GROMOS96N, GSTAT, HIN, ICON8,
*              IDATM, JAGUARIN, JAGUAROUT, M3D, MACCS, MACMOL,
*              MICROWORLD, MM2IN, MM2OUT, MM3, MMADS, MMCIF, MMD,
*              MOL2, MOLDEN, MOLIN, MOLINVENT, MOPACCART, MOPACINT,
*              MOPACOUT, MPQC, MSF, NWCHEMIN, NWCHEMOUT, OEBINARY,
*              PCMODEL, PDB, PREP, QCHEMIN, QCHEMOUT, REPORT,
*              SCHAKAL, SDF, SHELX, SKC, SMI, SPARTAN, SPARTANMM,
*              SPARTANSEMI, TGF, TINKER, TITLE, UNICHEM, VIEWMOL,
*              XED, XYZ
* }\endcode
*
* The following lines of code show how to set the input and output
* types of an OBMol through the member functions:
\code
   OBMol mol;
   mol.SetInputType(SDF);
   mol.SetOutputType(MOL2);
\endcode
* Once a molecule has been read into an OBMol the atoms and bonds
* can be accessed by the following methods:
\code
 OBAtom *atom;
 atom = mol.GetAtom(5); //random access of an atom
\endcode
*  or
\code
  OBBond *bond;
bond = mol.GetBond(14); //random access of a bond
\endcode
* or
\code
   OBAtom *atom;
   vector<OBNodeBase*>::iterator i;
   for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i)) //iterator access
\endcode
*or
\code
   OBBond *bond;
   vector<OBEdgeBase*>::iterator i;
   for (bond = mol.BeginBond(i);bond;bond = mol.NextBond(i)) //iterator access
\endcode
*It is important to note that atom arrays begin at one and bond arrays
* begin at zero. Requesting atom zero (\code
OEAtom *atom = mol.GetAtom(0); \endcode
*will result in an error, but
\code
   OEBond *bond = mol.GetBond(0);
\endcode
* is perfectly valid.
* Note that this is expected to change in the near future to simplify coding
* and improve efficiency.

* The ambiguity of numbering issues and off-by-one errors led to the use
* of iterators in Open Babel. An iterator is essentially just a pointer, but
* when used in conjunction with Standard Template Library (STL) vectors
* it provides an unambiguous way to loop over arrays. OBMols store their
* atom and bond information in STL vectors. Since vectors are template
* based, a vector of any user defined type can be declared. OBMols declare
* vector<OBAtom*> and vector<OBBond*> to store atom and bond information.
* Iterators are then a natural way to loop over the vectors of atoms and bonds.
*
*/
class OBMol : public OBGraphBase
{
protected:
  int                           _flags;
  bool                          _autoPartialCharge;
  bool                          _autoFormalCharge;
  bool                          _compressed;
  io_type                       _itype;
  io_type                       _otype;
  std::string                   _title;
  //vector<OBAtom*>             _atom;
  //vector<OBBond*>             _bond;
  std::vector<OBResidue*>       _residue;
  std::vector<OBGenericData*>   _vdata;
  float                         _energy;
  float                        *_c;
  std::vector<float*>           _vconf;
  unsigned short int            _natoms;
  unsigned short int            _nbonds;
  unsigned short int            _mod;
  unsigned short int            _access;

  bool  HasFlag(int flag)  {return((_flags & flag) ? true : false);}
  void  SetFlag(int flag)  {_flags |= flag;}

public:

    int GetMod() {return(_mod);}

    //***initialization and data (re)size methods***
    virtual ~OBMol();
    OBMol(io_type=UNDEFINED,io_type=UNDEFINED);
    OBMol(const OBMol &);
    OBMol &operator=(const OBMol &mol);
    OBMol &operator+=(const OBMol &mol);
    void IncrementMod() {_mod++;}
    void DecrementMod() {_mod--;}
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
    void SetAutomaticFormalCharge(bool val)     {_autoFormalCharge=val;}
    void SetAutomaticPartialCharge(bool val)    {_autoPartialCharge=val;}
    bool AutomaticFormalCharge()             {return(_autoFormalCharge);}
    bool AutomaticPartialCharge()            {return(_autoPartialCharge);}
    virtual bool Compress(void);
    virtual bool UnCompress(void);
    
    //***methods for handling generic data***
    bool                              HasData(std::string &);
    bool                              HasData(const char *);
    bool                              HasData(obDataType);
    void                              DeleteData(obDataType);
    void                              DeleteData(OBGenericData*);
    void                              DeleteData(std::vector<OBGenericData*>&);
    void                              SetData(OBGenericData *d) {_vdata.push_back(d);     }
    unsigned int                      DataSize()                {return(_vdata.size());   }
    OBGenericData                    *GetData(obDataType);
    OBGenericData                    *GetData(std::string&);
    OBGenericData                    *GetData(const char *);
    std::vector<OBGenericData*>           &GetData()                 {return(_vdata);         }
    std::vector<OBGenericData*>::iterator  BeginData()               {return(_vdata.begin()); }
    std::vector<OBGenericData*>::iterator  EndData()                 {return(_vdata.end());   }

    //OBMol data retrieval methods
    int          GetFlags()                           {return(_flags);}
    bool         GetGTDVector(std::vector<int> &);
    void         GetGIVector(std::vector<unsigned int> &);
    void         GetGIDVector(std::vector<unsigned int> &);
    const char  *GetTitle()                           {return(_title.c_str());}
    io_type      GetInputType()                       {return(_itype);}
    io_type      GetOutputType()                      {return(_otype);}
    unsigned int NumAtoms()                           {return(_natoms);}
    unsigned int NumBonds()                           {return(_nbonds);}
    unsigned int NumHvyAtoms();
    unsigned int NumResidues()                        {return(_residue.size());}
    unsigned int NumRotors();
    OBAtom      *GetAtom(int);
    OBAtom      *GetFirstAtom();
    OBBond      *GetBond(int); 
    OBBond      *GetBond(int, int);
    OBBond      *GetBond(OBAtom*,OBAtom*);
    OBResidue   *GetResidue(int);
    float        GetTorsion(int,int,int,int);
    float        GetTorsion(OBAtom*,OBAtom*,OBAtom*,OBAtom*);
    void         SetTorsion(OBAtom*,OBAtom*,OBAtom*,OBAtom*,float);
    float        GetEnergy()                          {return(_energy);}
    float        GetMolWt();
    float       *GetCoordinates()                     {return(_c);}
    std::vector<OBRing*> &GetSSSR();
    bool         IsCompressed()                       {return _compressed;}

    //***data modification methods***
    void   SetTitle(const char *title)     {_title = title;}
    void   SetTitle(std::string &title)    {_title = title;}
    void   SetEnergy(float energy)         {_energy = energy;}
    void   SetInputType(io_type type)      {_itype = type;}
    void   SetOutputType(io_type type)     {_otype = type;}
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

    //***molecule modification methods***
    bool Clear();
    void RenumberAtoms(std::vector<OBNodeBase*>&);
    void ToInertialFrame(int,float*);
    void ToInertialFrame();
    void Translate(const Vector &v);
    void Translate(const Vector &v,int);
    void Rotate(const float u[3][3]);
    void Rotate(const float m[9]);
    void Rotate(const float m[9],int nconf);
    void Center();
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
    Vector Center(int nconf);

    //***molecule utilities and perception methods***
    void FindSSSR();
    void FindRingAtomsAndBonds();
    void FindChiralCenters();
    void FindChildren(std::vector<int> &,int,int);
    void FindChildren(std::vector<OBAtom*>&,OBAtom*,OBAtom*);
    void FindLargestFragment(OBBitVec &);
    void ContigFragList(std::vector<std::vector<int> >&);
    void Align(OBAtom*,OBAtom*,Vector&,Vector&);
    void ConnectTheDots();
    void PerceiveBondOrders();
    void FindTorsions();
//  void ConnectTheDotsSort();

    virtual void BeginAccess(void);
    virtual void EndAccess(void);
    virtual void BeginModify(void);
    virtual void EndModify(bool nukePerceivedData=true);

    //*** methods to check for existence of properties
    bool Has2D();
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


    //***property request methods***
    bool Empty() {return(_natoms == 0);}

    //***iterator methods***
    OBAtom *BeginAtom(std::vector<OBNodeBase*>::iterator &i);
    OBAtom *NextAtom(std::vector<OBNodeBase*>::iterator &i);
    OBBond *BeginBond(std::vector<OBEdgeBase*>::iterator &i); 
    OBBond *NextBond(std::vector<OBEdgeBase*>::iterator &i); 
    OBResidue *BeginResidue(std::vector<OBResidue*>::iterator &i)
      {i = _residue.begin();return((i == _residue.end()) ? NULL:*i);}
    OBResidue *NextResidue(std::vector<OBResidue*>::iterator &i)
      {i++;return((i == _residue.end()) ? NULL:*i);}

    //*** MultiConf member functions
    int     NumConformers()         {return((_vconf.empty())?0:_vconf.size());}
    void    SetConformers(std::vector<float*> &v);
    void    AddConformer(float *f)                  {_vconf.push_back(f);}
    void    SetConformer(int i)                     {_c = _vconf[i];}
    void    CopyConformer(float*,int);
    void    CopyConformer(double*,int);
    void    DeleteConformer(int);
    float  *GetConformer(int i)                     {return(_vconf[i]);}
    float  *BeginConformer(std::vector<float*>::iterator&);
    float  *NextConformer(std::vector<float*>::iterator&);
    std::vector<float*> &GetConformers()                   {return(_vconf);}

    //misc bond functions
    void AssignResidueBonds(OBBitVec &);
    //void SortBonds() {sort(_vbond.begin(),_vbond.end(),CompareBonds);}

};

/// Class to store internal coordinate values
class OBInternalCoord 
{
public:
  //class members
  OBAtom *_a,*_b,*_c;
  float   _dst,_ang,_tor;
  //constructor
  OBInternalCoord(OBAtom *a=(OBAtom*)NULL,
                  OBAtom *b=(OBAtom*)NULL,
                  OBAtom *c=(OBAtom*)NULL)
    {
      _a = a; _b = b; _c = c;
      _dst = _ang = _tor = 0.0f;
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


#ifdef WIN32
inline float rint(float x) { return ( (x < 0.0f) ? ceil(x-0.5) : floor(x+0.5f));}
#endif

#ifndef __KCC
extern "C"{
void  get_rmat(float*,float*,float*,int);
void  ob_make_rmat(float mat[3][3],float rmat[9]);
void  qtrfit (float *r,float *f,int size,float u[3][3]);
float superimpose(float*,float*,int);
}
#else
void get_rmat(float*,float*,float*,int);
void ob_make_rmat(float mat[3][3],float rmat[9]);
void qtrfit (float *r,float *f,int size,float u[3][3]);
float superimpose(float*,float*,int);
#endif // __KCC


#define obAssert(__b__) \
    if (!(__b__)) {   \
        cerr << "Assert at File " << __FILE__ << " Line " << __LINE__ << endl; \
        int *p= NULL; *p= 10; \
        exit(-1); \
    }

}

#endif // OB_MOL_H

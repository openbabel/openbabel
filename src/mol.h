/**********************************************************************
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#ifdef WIN32
#pragma warning (disable: 4786) // warning: long & complicated stl warning
#pragma warning (disable: 4305) // warning: truncation from 'const double' to 'const float'
#endif

#ifndef MOL_H
#define MOL_H

#include <math.h>

#ifdef __sgi
#include <iostream.h>
#include <fstream.h>
#else
#include <iostream>
#include <fstream>
#endif

#include <algorithm>
#include <vector>
#include <string>

#include "base.h"
#include "data.h"
#include "chains.h"
#include "Vector.h"
#include "bitvec.h"
#include "ring.h"
#include "generic.h"
#include "typer.h"
#include "fileformat.h"
#include "binary_io.h"
#include "ctransform.h"

namespace OpenBabel {


class OBBond;
class OBAtom;
class OBResidue;
class OBMol;
class OBInternalCoord;


#define obAssert(__b__) \
   if (!(__b__)) {   \
       cerr << "Assert at File " << __FILE__ << " Line " << __LINE__ << endl; \
       int *p= NULL; *p= 10;                                                  \
       exit(-1);                                                              \
   }


/*!
**\brief An object to hold pose information.
*/
class OBPose {
    protected:
        unsigned int _conf;
        OBCoordTrans _ct;
    public:
        //Constructor, destructor and copy constructor
        OBPose() { _conf = 0; }
        OBPose(const OBPose& cp) {*this=cp;}
        virtual ~OBPose() { _ct.Clear();}

        //Operator overloads
        OBPose& operator=(const OBPose& cp) {
            _ct = cp._ct;
            _conf = cp._conf;
            return *this;
          }

        //Clear
        void Clear() {_ct.Clear();  _conf = 0;}

        //Set functions
        void SetConformer(unsigned int conf) {_conf = conf;}
        void SetTransform(OBCoordTrans ct) {_ct = ct;}

        //Get functions
        unsigned int ConformerNum() {return _conf;}
        OBCoordTrans& CoordTrans() {return _ct;}

        //Binary read and write
        unsigned int WriteBinary(char *ccc) {
            unsigned int idx=0;
            idx += _ct.WriteBinary(&ccc[idx]);
            idx += OB_io_write_binary(&ccc[idx],(char*)&_conf,sizeof(float),1);
            return idx;
          }
        unsigned int ReadBinary(char *ccc) {
            unsigned int idx=0;
            idx += _ct.ReadBinary(&ccc[idx]);
            idx += OB_io_read_binary(&ccc[idx],(char*)&_conf,sizeof(float),1);
            return idx;
          }
        void WriteBinary(std::ostream& ostr) {
            _ct.WriteBinary(ostr);
            OB_io_write_binary(ostr,(char*)&_conf,sizeof(float),1);
          } 
        void ReadBinary(std::istream& istr) {
            _ct.ReadBinary(istr);
            OB_io_read_binary(istr,(char*)&_conf,sizeof(float),1);
          }
  };


class OBResidue
{
public:

  OBResidue(void);
  virtual ~OBResidue(void);

  OBResidue &operator=(const OBResidue &);

  void AddAtom(OBAtom *atom); 
  void InsertAtom(OBAtom *atom);
  void RemoveAtom(OBAtom *atom);
  void Clear(void);

  void    SetName(std::string &resname)               { _resname  = resname;  }
  void    SetNum(unsigned short resnum)          { _resnum   = resnum;   }
  void    SetChain(char chain)                   { _chain    = chain;    }
  void    SetChainNum(unsigned short chainnum)   { _chainnum = chainnum; }
  void    SetIdx(unsigned short idx)             { _idx      = idx;      }

  void    SetAtomID(OBAtom *atom, std::string &id)    { _atomid[GetIndex(atom)] = id;     }
  void    SetHetAtom(OBAtom *atom, bool hetatm)  { _hetatm[GetIndex(atom)] = hetatm; }
  void    SetSerialNum(OBAtom *atom, unsigned short sernum) { _sernum[GetIndex(atom)] = sernum; }

  std::string&        GetName(void)                   { return(_resname);  }
  unsigned short GetNum(void)                    { return(_resnum);   }
  char           GetChain(void)                  { return(_chain);    }
  unsigned short GetChainNum(void)               { return(_chainnum); }
  unsigned short GetIdx(void)                    { return(_idx);      }

  std::string&        GetAtomID(OBAtom *atom)         { return(_atomid[GetIndex(atom)]); }
  unsigned short GetSerialNum(OBAtom *atom)      { return(_sernum[GetIndex(atom)]); }
  bool           IsHetAtom(OBAtom *atom)         { return(_hetatm[GetIndex(atom)]); }

  OBAtom *BeginAtom(std::vector<OBAtom*>::iterator &i) 
    {i = _atoms.begin();return((i ==_atoms.end()) ? NULL:*i);}
  OBAtom *NextAtom(std::vector<OBAtom*>::iterator &i) 
    {i++;return((i ==_atoms.end()) ? NULL:*i);}

protected: // methods

  int GetIndex(OBAtom *atom) 
  {
    for(unsigned int i=0;i<_atoms.size();i++) if (_atoms[i] == atom) return i;
    return -1; 
  }

protected: // members

  unsigned short         _idx;
  char                   _chain;
  unsigned short         _chainnum;
  unsigned short         _resnum;
  std::string                 _resname;

  std::vector<bool>           _hetatm;
  std::vector<std::string>    _atomid;
  std::vector<OBAtom*>        _atoms;
  std::vector<unsigned short> _sernum;

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
      
    float      x()                {return((*_c)[_cidx]);}
    float      y()                {return((*_c)[_cidx+1]);}
    float      z()                {return((*_c)[_cidx+2]);}
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

class OBMol : public OBGraphBase
{
protected:
  int                           _flags;
  bool                          _autoPartialCharge;
  bool                          _autoFormalCharge;
  bool                          _compressed;
  io_type                       _itype;
  io_type                       _otype;
  std::string                        _title;
  //vector<OBAtom*>               _atom;
  //vector<OBBond*>               _bond;
  std::vector<OBResidue*>            _residue;
  std::vector<OBGenericData*>        _vdata;
  float                         _energy;
  float                        *_c;
  std::vector<float*>                _vconf;
  std::vector<OBPose>                _pose;
  float                        *_xyz_pose;
  unsigned int                  _cur_pose_idx;
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
    void   SetTitle(char *title)           {_title = title;}
    void   SetTitle(std::string &title)         {_title = title;}
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

    //Pose Member functions
    unsigned int NumPoses() const {return ((_pose.empty())?0:_pose.size());}
    unsigned int CurrentPoseIndex();
    void DeletePoses();
    void DeletePose(unsigned int i);
    void AddPose(OBPose& pose); 
    void SetPoses(std::vector<OBPose>& poses);
    void SetPose(unsigned int i);
    void GetPoseCoordinates(unsigned int i, float *xyz);
    OBPose& GetPose(unsigned int i);
    void ChangePosesToConformers();


    //misc bond functions
    void AssignResidueBonds(OBBitVec &);
    //void SortBonds() {sort(_vbond.begin(),_vbond.end(),CompareBonds);}

};

class OBInternalCoord //class to store internal coordinate values
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

bool tokenize(std::vector<std::string>&, char *buf,char *delimstr=" \t\n");
bool tokenize(std::vector<std::string> &,std::string&,char*,int limit=-1);
void ThrowError(char *str);
void ThrowError(std::string &str);
void CartesianToInternal(std::vector<OBInternalCoord*>&,OBMol&);
void InternalToCartesian(std::vector<OBInternalCoord*>&,OBMol&);
std::string NewExtension(std::string&,char*);
bool SetInputType(OBMol&,std::string&);
bool SetOutputType(OBMol&,std::string&);

void qtrfit (float *r,float *f,int size,float u[3][3]);
float superimpose(float*,float*,int);
void get_rmat(float*,float*,float*,int);
void ob_make_rmat(float mat[3][3],float rmat[9]);

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

}

#endif //MOL_H

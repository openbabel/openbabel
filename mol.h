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

using namespace std;

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


class OEBond;
class OEAtom;
class OEResidue;
class OEMol;
class OEInternalCoord;


#define oeAssert(__b__) \
   if (!(__b__)) {   \
       cerr << "Assert at File " << __FILE__ << " Line " << __LINE__ << endl; \
       int *p= NULL; *p= 10;                                                  \
       exit(-1);                                                              \
   }


/*!
**\brief An object to hold pose information.
*/
class OEPose {
    protected:
        unsigned int _conf;
        OECoordTrans _ct;
    public:
        //Constructor, destructor and copy constructor
        OEPose() { _conf = 0; }
        OEPose(const OEPose& cp) {*this=cp;}
        virtual ~OEPose() { _ct.Clear();}

        //Operator overloads
        OEPose& operator=(const OEPose& cp) {
            _ct = cp._ct;
            _conf = cp._conf;
            return *this;
          }

        //Clear
        void Clear() {_ct.Clear();  _conf = 0;}

        //Set functions
        void SetConformer(unsigned int conf) {_conf = conf;}
        void SetTransform(OECoordTrans ct) {_ct = ct;}

        //Get functions
        unsigned int ConformerNum() {return _conf;}
        OECoordTrans& CoordTrans() {return _ct;}

        //Binary read and write
        unsigned int WriteBinary(char *ccc) {
            unsigned int idx=0;
            idx += _ct.WriteBinary(&ccc[idx]);
            idx += OE_io_write_binary(&ccc[idx],(char*)&_conf,sizeof(float),1);
            return idx;
          }
        unsigned int ReadBinary(char *ccc) {
            unsigned int idx=0;
            idx += _ct.ReadBinary(&ccc[idx]);
            idx += OE_io_read_binary(&ccc[idx],(char*)&_conf,sizeof(float),1);
            return idx;
          }
        void WriteBinary(ostream& ostr) {
            _ct.WriteBinary(ostr);
            OE_io_write_binary(ostr,(char*)&_conf,sizeof(float),1);
          } 
        void ReadBinary(istream& istr) {
            _ct.ReadBinary(istr);
            OE_io_read_binary(istr,(char*)&_conf,sizeof(float),1);
          }
  };


class OEResidue
{
public:

  OEResidue(void);
  virtual ~OEResidue(void);

  OEResidue &operator=(const OEResidue &);

  void AddAtom(OEAtom *atom); 
  void InsertAtom(OEAtom *atom);
  void RemoveAtom(OEAtom *atom);
  void Clear(void);

  void    SetName(string &resname)               { _resname  = resname;  }
  void    SetNum(unsigned short resnum)          { _resnum   = resnum;   }
  void    SetChain(char chain)                   { _chain    = chain;    }
  void    SetChainNum(unsigned short chainnum)   { _chainnum = chainnum; }
  void    SetIdx(unsigned short idx)             { _idx      = idx;      }

  void    SetAtomID(OEAtom *atom, string &id)    { _atomid[GetIndex(atom)] = id;     }
  void    SetHetAtom(OEAtom *atom, bool hetatm)  { _hetatm[GetIndex(atom)] = hetatm; }
  void    SetSerialNum(OEAtom *atom, unsigned short sernum) { _sernum[GetIndex(atom)] = sernum; }

  string&        GetName(void)                   { return(_resname);  }
  unsigned short GetNum(void)                    { return(_resnum);   }
  char           GetChain(void)                  { return(_chain);    }
  unsigned short GetChainNum(void)               { return(_chainnum); }
  unsigned short GetIdx(void)                    { return(_idx);      }

  string&        GetAtomID(OEAtom *atom)         { return(_atomid[GetIndex(atom)]); }
  unsigned short GetSerialNum(OEAtom *atom)      { return(_sernum[GetIndex(atom)]); }
  bool           IsHetAtom(OEAtom *atom)         { return(_hetatm[GetIndex(atom)]); }

  OEAtom *BeginAtom(vector<OEAtom*>::iterator &i) 
    {i = _atoms.begin();return((i ==_atoms.end()) ? NULL:*i);}
  OEAtom *NextAtom(vector<OEAtom*>::iterator &i) 
    {i++;return((i ==_atoms.end()) ? NULL:*i);}

protected: // methods

  int GetIndex(OEAtom *atom) 
  {
    for(unsigned int i=0;i<_atoms.size();i++) if (_atoms[i] == atom) return i;
    return -1; 
  }

protected: // members

  unsigned short         _idx;
  char                   _chain;
  unsigned short         _chainnum;
  unsigned short         _resnum;
  string                 _resname;

  vector<bool>           _hetatm;
  vector<string>         _atomid;
  vector<OEAtom*>        _atoms;
  vector<unsigned short> _sernum;

};

//ATOM Property Macros
#define OE_4RING_ATOM     (1<<1)
#define OE_3RING_ATOM     (1<<2)
#define OE_AROMATIC_ATOM  (1<<3)
#define OE_RING_ATOM      (1<<4)
#define OE_CSTEREO_ATOM   (1<<5)
#define OE_ACSTEREO_ATOM  (1<<6)
#define OE_DONOR_ATOM     (1<<7)
#define OE_ACCEPTOR_ATOM  (1<<8)
#define OE_CHIRAL_ATOM    (1<<9)

class OEAtom : public OENodeBase
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
  OEResidue                    *_residue;
//  OEMol                        *_parent;
  //vector<OEBond*>               _bond;

    int  GetFlag()          const {return(_flags);}
    void SetFlag(int flag) {_flags |= flag;}
    bool HasFlag(int flag) {return((_flags & flag) ? true : false);}

public:
    OEAtom();
    virtual ~OEAtom();
    OEAtom &operator=(OEAtom &);
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
    void SetType(string &type)               {strcpy(_type,type.c_str());}
    void SetPartialCharge(float pcharge)     {_pcharge = pcharge;}
    void SetVector();
    void SetVector(Vector &v);                
    void SetVector(const float,const float,const float);
    void SetResidue(OEResidue *res)          {_residue=res;}
//    void SetParent(OEMol *ptr)               {_parent=ptr;}
    void SetCoordPtr(float **c)              {_c = c;_cidx = (GetIdx()-1)*3;}
    void SetAromatic()                       {SetFlag(OE_AROMATIC_ATOM);}
    void UnsetAromatic()              {_flags &= (~(OE_AROMATIC_ATOM));}
    void SetClockwiseStereo()          {SetFlag(OE_CSTEREO_ATOM|OE_CHIRAL_ATOM);}
    void SetAntiClockwiseStereo()      {SetFlag(OE_ACSTEREO_ATOM|OE_CHIRAL_ATOM);}
    void UnsetStereo()                 
      {
        _flags &= ~(OE_ACSTEREO_ATOM);
        _flags &= ~(OE_CSTEREO_ATOM);
        _flags &= ~(OE_CHIRAL_ATOM);
      }
    void SetInRing()                         {SetFlag(OE_RING_ATOM);}
    void SetChiral()                         {SetFlag(OE_CHIRAL_ATOM);}
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
    OEResidue *GetResidue();
//    OEMol     *GetParent()        {return((OEMol*)_parent);}
    bool       GetNewBondVector(Vector&,float);
    OEBond    *GetBond(OEAtom *);
    OEAtom    *GetNextAtom();

    //***iterator methods***
    vector<OEEdgeBase*>::iterator BeginBonds() {return(_vbond.begin());}
    vector<OEEdgeBase*>::iterator EndBonds()   {return(_vbond.end());}
    OEBond *BeginBond(vector<OEEdgeBase*>::iterator &i);
    OEBond *NextBond(vector<OEEdgeBase*>::iterator &i); 
    OEAtom *BeginNbrAtom(vector<OEEdgeBase*>::iterator &); 
    OEAtom *NextNbrAtom(vector<OEEdgeBase*>::iterator &); 

    //***addition of residue/bond info. for an atom***
    void NewResidue()                                          {if (!_residue) _residue = new OEResidue;}
    void DeleteResidue()                                       {if (_residue) delete _residue;}
    void AddBond(OEBond *bond)                                 {_vbond.push_back((OEEdgeBase*)bond);}
    void InsertBond(vector<OEEdgeBase*>::iterator &i, OEBond *bond)
      {_vbond.insert(i, (OEEdgeBase*)bond);}
    bool DeleteBond(OEBond*);

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
    bool IsConnected(OEAtom*);
    bool IsOneThree(OEAtom*);
    bool IsOneFour(OEAtom*);
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
    bool IsClockwise()        {return(HasFlag(OE_CSTEREO_ATOM));}
    bool IsAntiClockwise()    {return(HasFlag(OE_ACSTEREO_ATOM));}
    bool HasChiralitySpecified() 
                              {return(HasFlag(OE_CSTEREO_ATOM|OE_ACSTEREO_ATOM));}
    bool HasAlphaBetaUnsat(bool includePandS=true);
    bool HasBondOfOrder(unsigned int);
    int  CountBondsOfOrder(unsigned int);
    bool HasSingleBond()      {return(HasBondOfOrder(1));}
    bool HasDoubleBond()      {return(HasBondOfOrder(2));}
    bool HasAromaticBond()    {return(HasBondOfOrder(4));}
};

//BOND Property Macros
#define OE_AROMATIC_BOND  (1<<1)
#define OE_WEDGE_BOND     (1<<2)
#define OE_HASH_BOND      (1<<3)
#define OE_RING_BOND      (1<<4)
#define OE_TORUP_BOND     (1<<5)
#define OE_TORDOWN_BOND   (1<<6)
#define OE_KSINGLE_BOND   (1<<7) 
#define OE_KDOUBLE_BOND   (1<<8)
#define OE_KTRIPLE_BOND   (1<<9)
#define OE_CLOSURE_BOND   (1<<10)

class OEBond : public OEEdgeBase
{
protected:
  char                    _order;
  unsigned short int      _flags;
//  OEAtom                 *_bgn;
//  OEAtom                 *_end;
//  OEMol                  *_parent;
//  unsigned short int      _idx;

  bool HasFlag(int flag) {return((_flags & flag) != 0);}
  void SetFlag(int flag) {_flags |= flag;}

public:
   OEBond();
   virtual ~OEBond() {}
   //***bond modification methods***
   void SetIdx(int idx)                     {_idx = idx;}
   void SetBO(int order);
   void SetBegin(OEAtom *begin)             {_bgn = begin;}
   void SetEnd(OEAtom *end)                 {_end = end;}
//   void SetParent(OEMol *ptr)               {_parent=ptr;}
   void SetLength(OEAtom*,float);
   void Set(int,OEAtom*,OEAtom*,int,int);
   void SetKSingle();
   void SetKDouble();
   void SetKTriple();
   void SetAromatic()                       {SetFlag(OE_AROMATIC_BOND);}
   void SetUp()                             {SetFlag(OE_TORUP_BOND);}
   void SetDown()                           {SetFlag(OE_TORDOWN_BOND);}
   void SetInRing()                         {SetFlag(OE_RING_BOND);}
   void SetClosure()                        {SetFlag(OE_CLOSURE_BOND);}
   void UnsetAromatic()                     {_flags &= (~(OE_AROMATIC_BOND));}

   
   //***bond data request methods***
   unsigned int     GetBO()            const    {return((int)_order);}
   unsigned int     GetBondOrder()     const    {return((int)_order);}
   unsigned int     GetFlags()         const    {return(_flags);}
   unsigned int     GetBeginAtomIdx()  const    {return(_bgn->GetIdx());}
   unsigned int     GetEndAtomIdx()    const    {return(_end->GetIdx());}
   OEAtom *GetBeginAtom()              {return((OEAtom*)_bgn);}
   OEAtom *GetEndAtom()                {return((OEAtom*)_end);}
   OEAtom *GetNbrAtom(OEAtom *ptr)     {return((ptr != _bgn)? (OEAtom*)_bgn : (OEAtom*)_end);}
//   OEMol  *GetParent()                 {return(_parent);}
   float   GetEquibLength();
   float   GetLength();
   int     GetNbrAtomIdx(OEAtom *ptr) 
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
   bool IsUp()                        {return(HasFlag(OE_TORUP_BOND));}
   bool IsDown()                      {return(HasFlag(OE_TORDOWN_BOND));}
   bool IsWedge()                     {return(HasFlag(OE_WEDGE_BOND));}
   bool IsHash()                      {return(HasFlag(OE_HASH_BOND));}
};

//prototype necessary for OEMol
//bool CompareBonds(const OEBond &*,const OEBond &*);


//MOL Property Macros
#define OE_SSSR_MOL              (1<<1)
#define OE_RINGFLAGS_MOL         (1<<2)
#define OE_AROMATIC_MOL          (1<<3)
#define OE_ATOMTYPES_MOL         (1<<4)
#define OE_CHIRALITY_MOL         (1<<5)
#define OE_PCHARGE_MOL           (1<<6)
#define OE_HYBRID_MOL            (1<<8)
#define OE_IMPVAL_MOL            (1<<9)
#define OE_KEKULE_MOL            (1<<10)
#define OE_CLOSURE_MOL           (1<<11)
#define OE_H_ADDED_MOL           (1<<12)
#define OE_PH_CORRECTED_MOL      (1<<13)
#define OE_AROM_CORRECTED_MOL    (1<<14)
#define OE_CHAINS_MOL            (1<<15)
#define OE_CURRENT_CONFORMER -1

class OEMol : public OEGraphBase
{
protected:
  int                           _flags;
  bool                          _autoPartialCharge;
  bool                          _autoFormalCharge;
  bool                          _compressed;
  io_type                       _itype;
  io_type                       _otype;
  string                        _title;
  //vector<OEAtom*>               _atom;
  //vector<OEBond*>               _bond;
  vector<OEResidue*>            _residue;
  vector<OEGenericData*>        _vdata;
  float                         _energy;
  float                        *_c;
  vector<float*>                _vconf;
  vector<OEPose>                _pose;
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
    virtual ~OEMol();
    OEMol(io_type=UNDEFINED,io_type=UNDEFINED);
    OEMol(const OEMol &);
    OEMol &operator=(const OEMol &mol);
    OEMol &operator+=(const OEMol &mol);
    void IncrementMod() {_mod++;}
    void DecrementMod() {_mod--;}
    void ReserveAtoms(int natoms) {if (natoms && _mod) _vatom.reserve(natoms);}
      
    virtual OEAtom *CreateAtom(void);
    virtual OEBond *CreateBond(void);
    virtual void DestroyAtom(OENodeBase*);
    virtual void DestroyBond(OEEdgeBase*);
    bool AddAtom(OEAtom&);
    bool AddBond(int,int,int,int flags=0,int insertpos=-1);
    bool AddBond(OEBond&);
    bool AddResidue(OEResidue&);
    bool InsertAtom(OEAtom &);
    bool DeleteAtom(OEAtom*);
    bool DeleteBond(OEBond*);
    bool DeleteResidue(OEResidue*);
    OEAtom    *NewAtom();
    OEResidue *NewResidue();
    void SetAutomaticFormalCharge(bool val)     {_autoFormalCharge=val;}
    void SetAutomaticPartialCharge(bool val)    {_autoPartialCharge=val;}
    bool AutomaticFormalCharge()             {return(_autoFormalCharge);}
    bool AutomaticPartialCharge()            {return(_autoPartialCharge);}
    virtual bool Compress(void);
    virtual bool UnCompress(void);
    
    //***methods for handling generic data***
    bool                              HasData(string &);
    bool                              HasData(const char *);
    bool                              HasData(oeDataType);
    void                              DeleteData(oeDataType);
    void                              DeleteData(OEGenericData*);
    void                              DeleteData(vector<OEGenericData*>&);
    void                              SetData(OEGenericData *d) {_vdata.push_back(d);     }
    unsigned int                      DataSize()                {return(_vdata.size());   }
    OEGenericData                    *GetData(oeDataType);
    OEGenericData                    *GetData(string&);
    OEGenericData                    *GetData(const char *);
    vector<OEGenericData*>           &GetData()                 {return(_vdata);         }
    vector<OEGenericData*>::iterator  BeginData()               {return(_vdata.begin()); }
    vector<OEGenericData*>::iterator  EndData()                 {return(_vdata.end());   }

    //OEMol data retrieval methods
    int          GetFlags()                           {return(_flags);}
    bool         GetGTDVector(vector<int> &);
	void         GetGIVector(vector<unsigned int> &);
	void         GetGIDVector(vector<unsigned int> &);
    const char  *GetTitle()                           {return(_title.c_str());}
    io_type      GetInputType()                       {return(_itype);}
    io_type      GetOutputType()                      {return(_otype);}
    unsigned int NumAtoms()                           {return(_natoms);}
    unsigned int NumBonds()                           {return(_nbonds);}
    unsigned int NumHvyAtoms();
    unsigned int NumResidues()                        {return(_residue.size());}
    unsigned int NumRotors();
    OEAtom      *GetAtom(int);
    OEAtom      *GetFirstAtom();
    OEBond      *GetBond(int); 
    OEBond      *GetBond(int, int);
    OEBond      *GetBond(OEAtom*,OEAtom*);
    OEResidue   *GetResidue(int);
    float        GetTorsion(int,int,int,int);
    float        GetTorsion(OEAtom*,OEAtom*,OEAtom*,OEAtom*);
    void         SetTorsion(OEAtom*,OEAtom*,OEAtom*,OEAtom*,float);
    float        GetEnergy()                          {return(_energy);}
    float        GetMolWt();
    float       *GetCoordinates()                     {return(_c);}
    vector<OERing*> &GetSSSR();
    bool         IsCompressed()                       {return _compressed;}

    //***data modification methods***
    void   SetTitle(char *title)           {_title = title;}
    void   SetTitle(string &title)         {_title = title;}
    void   SetEnergy(float energy)         {_energy = energy;}
    void   SetInputType(io_type type)      {_itype = type;}
    void   SetOutputType(io_type type)     {_otype = type;}
    void   SetAromaticPerceived()          {SetFlag(OE_AROMATIC_MOL);}
    void   SetSSSRPerceived()              {SetFlag(OE_SSSR_MOL);}
    void   SetRingAtomsAndBondsPerceived() {SetFlag(OE_RINGFLAGS_MOL);}
    void   SetAtomTypesPerceived()         {SetFlag(OE_ATOMTYPES_MOL);}
    void   SetChainsPerceived()            {SetFlag(OE_CHAINS_MOL);}
    void   SetChiralityPerceived()         {SetFlag(OE_CHIRALITY_MOL);}
    void   SetPartialChargesPerceived()    {SetFlag(OE_PCHARGE_MOL);}
    void   SetHybridizationPerceived()     {SetFlag(OE_HYBRID_MOL);}
    void   SetImplicitValencePerceived()   {SetFlag(OE_IMPVAL_MOL);}
    void   SetKekulePerceived()            {SetFlag(OE_KEKULE_MOL);}
    void   SetClosureBondsPerceived()      {SetFlag(OE_CLOSURE_MOL);}
    void   SetHydrogensAdded()             {SetFlag(OE_H_ADDED_MOL);}
    void   SetCorrectedForPH()             {SetFlag(OE_PH_CORRECTED_MOL);}
    void   SetAromaticCorrected()          {SetFlag(OE_AROM_CORRECTED_MOL);}
    void   UnsetAromaticPerceived()        {_flags &= (~(OE_AROMATIC_MOL));}
    void   UnsetPartialChargesPerceived()  {_flags &= (~(OE_PCHARGE_MOL));}
    void   UnsetImplicitValencePerceived() {_flags &= (~(OE_IMPVAL_MOL));}

    //***molecule modification methods***
    bool Clear();
    void RenumberAtoms(vector<OENodeBase*>&);
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
    bool DeleteHydrogen(OEAtom*);
    bool DeleteHydrogens();
    bool DeleteHydrogens(OEAtom*);
    bool DeleteNonPolarHydrogens();
    bool AddHydrogens(bool polaronly=false,bool correctForPH=true);
    bool AddHydrogens(OEAtom*);
    bool AddPolarHydrogens();
    bool StripSalts();
    bool CorrectForPH();
    Vector Center(int nconf);

    //***molecule utilities and perception methods***
    void FindSSSR();
    void FindRingAtomsAndBonds();
    void FindChiralCenters();
    void FindChildren(vector<int> &,int,int);
    void FindChildren(vector<OEAtom*>&,OEAtom*,OEAtom*);
    void FindLargestFragment(OEBitVec &);
    void ContigFragList(vector<vector<int> >&);
    void Align(OEAtom*,OEAtom*,Vector&,Vector&);
    void ConnectTheDots();
//  void ConnectTheDotsSort();

    virtual void BeginAccess(void);
    virtual void EndAccess(void);
    virtual void BeginModify(void);
    virtual void EndModify(bool nukePerceivedData=true);

    //*** methods to check for existence of properties
    bool Has2D();
    bool Has3D();
    bool HasNonZeroCoords();
    bool HasAromaticPerceived()          {return(HasFlag(OE_AROMATIC_MOL));       }
    bool HasSSSRPerceived()              {return(HasFlag(OE_SSSR_MOL));           }
    bool HasRingAtomsAndBondsPerceived() {return(HasFlag(OE_RINGFLAGS_MOL));      }
    bool HasAtomTypesPerceived()         {return(HasFlag(OE_ATOMTYPES_MOL));      }
    bool HasChiralityPerceived()         {return(HasFlag(OE_CHIRALITY_MOL));      }
    bool HasPartialChargesPerceived()    {return(HasFlag(OE_PCHARGE_MOL));        }
    bool HasHybridizationPerceived()     {return(HasFlag(OE_HYBRID_MOL));         }
    bool HasImplicitValencePerceived()   {return(HasFlag(OE_IMPVAL_MOL));         }
    bool HasKekulePerceived()            {return(HasFlag(OE_KEKULE_MOL));         }
    bool HasClosureBondsPerceived()      {return(HasFlag(OE_CLOSURE_MOL));        }
    bool HasChainsPerceived()            {return(HasFlag(OE_CHAINS_MOL));         }
    bool HasHydrogensAdded()             {return(HasFlag(OE_H_ADDED_MOL));        }
    bool HasAromaticCorrected()          {return(HasFlag(OE_AROM_CORRECTED_MOL)); }
    bool IsCorrectedForPH()              {return(HasFlag(OE_PH_CORRECTED_MOL));   }
    bool IsChiral();


    //***property request methods***
    bool Empty() {return(_natoms == 0);}

    //***iterator methods***
    OEAtom *BeginAtom(vector<OENodeBase*>::iterator &i);
    OEAtom *NextAtom(vector<OENodeBase*>::iterator &i);
    OEBond *BeginBond(vector<OEEdgeBase*>::iterator &i); 
    OEBond *NextBond(vector<OEEdgeBase*>::iterator &i); 
    OEResidue *BeginResidue(vector<OEResidue*>::iterator &i)
      {i = _residue.begin();return((i == _residue.end()) ? NULL:*i);}
    OEResidue *NextResidue(vector<OEResidue*>::iterator &i)
      {i++;return((i == _residue.end()) ? NULL:*i);}

    //*** MultiConf member functions
    int     NumConformers()         {return((_vconf.empty())?0:_vconf.size());}
    void    SetConformers(vector<float*> &v);
    void    AddConformer(float *f)                  {_vconf.push_back(f);}
    void    SetConformer(int i)                     {_c = _vconf[i];}
    void    CopyConformer(float*,int);
    void    CopyConformer(double*,int);
    void    DeleteConformer(int);
    float  *GetConformer(int i)                     {return(_vconf[i]);}
    float  *BeginConformer(vector<float*>::iterator&);
    float  *NextConformer(vector<float*>::iterator&);
    vector<float*> &GetConformers()                   {return(_vconf);}

    //Pose Member functions
    unsigned int NumPoses() const {return ((_pose.empty())?0:_pose.size());}
    unsigned int CurrentPoseIndex();
    void DeletePoses();
    void DeletePose(unsigned int i);
    void AddPose(OEPose& pose); 
    void SetPoses(vector<OEPose>& poses);
    void SetPose(unsigned int i);
    void GetPoseCoordinates(unsigned int i, float *xyz);
    OEPose& GetPose(unsigned int i);
    void ChangePosesToConformers();


    //misc bond functions
    void AssignResidueBonds(OEBitVec &);
    //void SortBonds() {sort(_vbond.begin(),_vbond.end(),CompareBonds);}

    //***IO operations***
    friend ostream& operator<< ( ostream&, OEMol& ) ;
    friend istream& operator>> ( istream&, OEMol& ) ;
};

class OEInternalCoord //class to store internal coordinate values
{
public:
  //class members
  OEAtom *_a,*_b,*_c;
  float   _dst,_ang,_tor;
  //constructor
  OEInternalCoord(OEAtom *a=(OEAtom*)NULL,
                  OEAtom *b=(OEAtom*)NULL,
                  OEAtom *c=(OEAtom*)NULL)
    {
      _a = a; _b = b; _c = c;
      _dst = _ang = _tor = 0.0f;
    }
};

//function prototypes
// (prototypes for Read/Write are in fileformat.h

bool tokenize(vector<string>&, char *buf,char *delimstr=" \t\n");
bool tokenize(vector<string> &,string&,char*,int limit=-1);
void ThrowError(char *str);
void ThrowError(string &str);
void CartesianToInternal(vector<OEInternalCoord*>&,OEMol&);
void InternalToCartesian(vector<OEInternalCoord*>&,OEMol&);
string NewExtension(string&,char*);
bool SetInputType(OEMol&,string&);
bool SetOutputType(OEMol&,string&);

//global definitions
extern  OEExtensionTable extab;
extern  OEElementTable   etab;
extern  OETypeTable      ttab;
extern  OEChainsParser   chainsparser;

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

#ifndef __KCC
extern "C"{
void qtrfit (float *r,float *f,int size,float u[3][3]);
float superimpose(float*,float*,int);
void get_rmat(float*,float*,float*,int);
void oe_make_rmat(float mat[3][3],float rmat[9]);
}
#else
void qtrfit (float *r,float *f,int size,float u[3][3]);
float superimpose(float*,float*,int);
void get_rmat(float *,float*,float *,int);
void oe_make_rmat(float mat[3][3],float rmat[9]);
#endif //__KCC

#endif //MOL_H

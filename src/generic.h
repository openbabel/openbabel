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

#ifndef GENERIC_H
#define GENERIC_H

//obData0 through obData9 are data slots that are not used in OpenBabel, and
//are meant for use in derivative programs.  Macro definitions can be used
//to define what each data slot is used for.
 
namespace OpenBabel {

class OBAtom;
class OBBond; 

enum obDataType {obUndefinedData,obPairData,obEnergyData,
		 obCommentData,obCompressData,obExternalBondData,obRotamerList,
		 obVirtualBondData,obRingData,obData0,obData1,obData2,obData3,
		 obData4,obData5,obData6,obData7,obData8,obData9};

//base class for generic data - use obData# slots for undefined data types

class OBGenericData
{
protected:
	std::string     _attr; //attribute tag
	obDataType _type;
public:
	OBGenericData();
	OBGenericData(const OBGenericData&);
	virtual ~OBGenericData() {}
	void                  SetAttribute(std::string &v)       {_attr = v;}
	virtual const std::string &GetAttribute()          const {return(_attr);}
	obDataType      GetDataType()                 const {return(_type);}
};

class OBCommentData : public OBGenericData
{
protected:
	std::string _data;
public:
	OBCommentData();
	OBCommentData(const OBCommentData&);
	void          SetData(std::string &data)        {_data = data; }
	void          SetData(const char *d)       {_data = d;    }
	const std::string &GetData()              const {return(_data);}
};

class OBExternalBond
{
  int     _idx;
  OBAtom *_atom;
  OBBond *_bond;
public:
  OBExternalBond() {}
  OBExternalBond(OBAtom *,OBBond *,int);
  OBExternalBond(const OBExternalBond &);
  ~OBExternalBond(){}

  int     GetIdx()  const {return(_idx);  }
  OBAtom *GetAtom() const {return(_atom); }
  OBBond *GetBond() const {return(_bond); }
  void SetIdx(int idx) {_idx = idx;}
  void SetAtom(OBAtom *atom) {_atom = atom;}
  void SetBond(OBBond *bond) {_bond = bond;}
};

class OBExternalBondData : public OBGenericData
{
protected:
  std::vector<OBExternalBond> _vexbnd;
public:
  OBExternalBondData();
  void SetData(OBAtom*,OBBond*,int);
	std::vector<OBExternalBond> *GetData() {return(&_vexbnd);}
};

class OBCompressData : public OBGenericData
{
protected:
  int _size;
  unsigned char *_data;
public:
  OBCompressData();
  ~OBCompressData();
  void SetData(unsigned char *,int);
  int            GetSize() {return(_size);}
  unsigned char *GetData() {return(_data);}
};

class OBPairData : public OBGenericData //use to store attribute/value relationships
{
 protected:
  std::string _value;
 public:
  OBPairData();
  void    SetValue(const char *v) {_value = v;}
  void    SetValue(std::string &v)     {_value = v;}
  std::string &GetValue()              {return(_value);}
};

class OBVirtualBond : public OBGenericData
{
protected:
	int _bgn;
	int _end;
	int _ord;
	int _stereo;
public:
	OBVirtualBond();           
	OBVirtualBond(int,int,int,int stereo=0);
	int GetBgn()                           {return(_bgn);}
	int GetEnd()                           {return(_end);}
	int GetOrder()                         {return(_ord);}
	int GetStereo()                        {return(_stereo);}
};

class OBRingData : public OBGenericData
{
protected:
	std::vector<OBRing*> _vr;
public:
	OBRingData();
	void SetData(std::vector<OBRing*> &vr) {_vr = vr;}
	void PushBack(OBRing *r)          {_vr.push_back(r);}
	std::vector<OBRing*> &GetData()        {return(_vr);}
};

} //end namespace OpenBabel

#endif //GENERIC_H

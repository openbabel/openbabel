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

//oeData0 through oeData9 are data slots that are not used in OELib, and
//are meant for use in derivative programs.  Macro definitions can be used
//to define what each data slot is used for.
 
namespace OpenBabel {

class OEAtom;
class OEBond; 

enum oeDataType {oeUndefinedData,oePairData,oeEnergyData,
		 oeCommentData,oeCompressData,oeExternalBondData,oeRotamerList,
		 oeVirtualBondData,oeRingData,oeData0,oeData1,oeData2,oeData3,
		 oeData4,oeData5,oeData6,oeData7,oeData8,oeData9};

//base class for generic data - use oeData# slots for undefined data types

class OEGenericData
{
protected:
	string     _attr; //attribute tag
	oeDataType _type;
public:
	OEGenericData();
	OEGenericData(const OEGenericData&);
	virtual ~OEGenericData() {}
	void                  SetAttribute(string &v)       {_attr = v;}
	virtual const string &GetAttribute()          const {return(_attr);}
	oeDataType      GetDataType()                 const {return(_type);}
};

class OECommentData : public OEGenericData
{
protected:
	string _data;
public:
	OECommentData();
	OECommentData(const OECommentData&);
	void          SetData(string &data)        {_data = data; }
	void          SetData(const char *d)       {_data = d;    }
	const string &GetData()              const {return(_data);}
};

class OEExternalBond
{
  int     _idx;
  OEAtom *_atom;
  OEBond *_bond;
public:
  OEExternalBond() {}
  OEExternalBond(OEAtom *,OEBond *,int);
  OEExternalBond(const OEExternalBond &);
  ~OEExternalBond(){}

  int     GetIdx()  const {return(_idx);  }
  OEAtom *GetAtom() const {return(_atom); }
  OEBond *GetBond() const {return(_bond); }
  void SetIdx(int idx) {_idx = idx;}
  void SetAtom(OEAtom *atom) {_atom = atom;}
  void SetBond(OEBond *bond) {_bond = bond;}
};

class OEExternalBondData : public OEGenericData
{
protected:
  vector<OEExternalBond> _vexbnd;
public:
  OEExternalBondData();
  void SetData(OEAtom*,OEBond*,int);
	vector<OEExternalBond> *GetData() {return(&_vexbnd);}
};

class OECompressData : public OEGenericData
{
protected:
  int _size;
  unsigned char *_data;
public:
  OECompressData();
  ~OECompressData();
  void SetData(unsigned char *,int);
  int            GetSize() {return(_size);}
  unsigned char *GetData() {return(_data);}
};

class OEPairData : public OEGenericData //use to store attribute/value relationships
{
 protected:
  string _value;
 public:
  OEPairData();
  void    SetValue(const char *v) {_value = v;}
  void    SetValue(string &v)     {_value = v;}
  string &GetValue()              {return(_value);}
};

class OEVirtualBond : public OEGenericData
{
protected:
	int _bgn;
	int _end;
	int _ord;
	int _stereo;
public:
	OEVirtualBond();           
	OEVirtualBond(int,int,int,int stereo=0);
	int GetBgn()                           {return(_bgn);}
	int GetEnd()                           {return(_end);}
	int GetOrder()                         {return(_ord);}
	int GetStereo()                        {return(_stereo);}
};

class OERingData : public OEGenericData
{
protected:
	vector<OERing*> _vr;
public:
	OERingData();
	void SetData(vector<OERing*> &vr) {_vr = vr;}
	void PushBack(OERing *r)          {_vr.push_back(r);}
	vector<OERing*> &GetData()        {return(_vr);}
};

} //end namespace OpenBabel

#endif //GENERIC_H

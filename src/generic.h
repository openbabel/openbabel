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

#ifndef OB_GENERIC_H
#define OB_GENERIC_H

#include <string>
#include <vector>

//obData0 through obData9 are data slots that are not used in OpenBabel, and
//are meant for use in derivative programs.  Macro definitions can be used
//to define what each data slot is used for.
 
namespace OpenBabel {

class OBAtom;
class OBBond; 
class OBRing;

enum obDataType {obUndefinedData,obPairData,obEnergyData,
		 obCommentData,obCompressData,obExternalBondData,obRotamerList,
		 obVirtualBondData,obRingData,obTorsionData,obAngleData,
		 obData0,obData1,obData2,obData3,
		 obData4,obData5,obData6,obData7,
		 obData8,obData9};

//base class for generic data - use obData# slots for undefined data types

class OBGenericData
{
protected:
	std::string     _attr; //attribute tag
	obDataType 	_type;
public:
	OBGenericData();
	OBGenericData(const OBGenericData&);
	virtual ~OBGenericData() {}
	OBGenericData& operator=(const OBGenericData &src);
	
	void                      SetAttribute(std::string &v)  {_attr = v;}
	virtual const std::string &GetAttribute()  const {return(_attr);}
	obDataType                GetDataType()    const {return(_type);}
};

class OBCommentData : public OBGenericData
{
protected:
	std::string _data;
public:
	OBCommentData();
	OBCommentData(const OBCommentData&);
	OBCommentData& operator=(const OBCommentData &src);

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
	OBRingData(const OBRingData &);
	~OBRingData();

	OBRingData &operator=(const OBRingData &);

	void SetData(std::vector<OBRing*> &vr) {_vr = vr;}
	void PushBack(OBRing *r)          {_vr.push_back(r);}
	std::vector<OBRing*> &GetData()        {return(_vr);}
};

/*!
**\brief A class which holds the torsion data for a single rotatable bond and all the four atom sets around it
*/
class OBTorsion
{
    friend class OBMol;
    friend class OBTorsionData;

protected:
	std::pair<OBAtom*,OBAtom*> _bc;
	std::vector<triple<OBAtom*,OBAtom*,float> > _ads; //float is angle in rads

  OBTorsion() { _bc.first=0; _bc.second=0; };  //protected for use only by friend classes
	OBTorsion(OBAtom *, OBAtom *, OBAtom *, OBAtom *);

	std::vector<quad<OBAtom*,OBAtom*,OBAtom*,OBAtom*> > GetTorsions();

public:
  OBTorsion(const OBTorsion &);
  ~OBTorsion() {}

  OBTorsion& operator=(const OBTorsion &);

  void Clear();
  bool Empty() { return(_bc.first == 0 && _bc.second == 0); }

  bool AddTorsion(OBAtom *a,OBAtom *b, OBAtom *c,OBAtom *d);
  bool AddTorsion(quad<OBAtom*,OBAtom*,OBAtom*,OBAtom*> &atoms);

  bool SetAngle(float radians, unsigned int index = 0);
  bool SetData(OBBond *bond);

  bool GetAngle(float &radians, unsigned int index =0);
  unsigned int GetBondIdx();
  unsigned int GetSize() const {return (unsigned int)_ads.size();}

  std::pair<OBAtom*,OBAtom*>                  GetBC()  { return(_bc); }
  std::vector<triple<OBAtom*,OBAtom*,float> > GetADs() { return(_ads) ;}

  bool IsProtonRotor();
};



/*!
**\brief The Object to hold Torsions as Generic data for OBMol
*/
class OBTorsionData : public OBGenericData 
{
    friend class OBMol;

protected:
    std::vector<OBTorsion> _torsions;

    OBTorsionData();
    OBTorsionData(const OBTorsionData &);

public:
    OBTorsionData &operator=(const OBTorsionData &);

    void Clear();

    std::vector<OBTorsion> GetData() const {return _torsions;                     }
    unsigned int      GetSize() const {return (unsigned int)_torsions.size();}

    void SetData(OBTorsion &torsion);

    bool FillTorsionArray(std::vector<std::vector<unsigned int> > &torsions);
};



/*!
**\brief A class which holds the 3 atoms in an angle and the angle itself
*/
class OBAngle 
{
    friend class OBMol;
    friend class OBAngleData;

protected:

	//member data

	OBAtom                *_vertex;
	std::pair<OBAtom*,OBAtom*>  _termini;
	float                  _radians;

	//protected member functions

    OBAngle();	//protect constructor for use only by friend classes
	OBAngle(OBAtom *vertex,OBAtom *a,OBAtom *b);

	triple<OBAtom*,OBAtom*,OBAtom*> GetAtoms();
	void SortByIndex();

public:

	OBAngle(const OBAngle &);
    ~OBAngle() { _vertex = NULL; }

    OBAngle &operator = (const OBAngle &);
	bool     operator ==(const OBAngle &);

	void  Clear();

    float GetAngle() const { return(_radians); }

    void  SetAngle(float radians) { _radians = radians; }
	void  SetAtoms(OBAtom *vertex,OBAtom *a,OBAtom *b);
	void  SetAtoms(triple<OBAtom*,OBAtom*,OBAtom*> &atoms);

};


/*!
**\brief The Object to hold Angle as Generic data for OBMol
*/
class OBAngleData : public OBGenericData
{
    friend class OBMol;

protected:
    std::vector<OBAngle> _angles;

    OBAngleData();
    OBAngleData(const OBAngleData &);
    std::vector<OBAngle> GetData() const { return(_angles); }

public:
    OBAngleData &operator =(const OBAngleData &);

    void Clear();
    unsigned int FillAngleArray(int **angles, unsigned int &size);
    void         SetData(OBAngle &);
    unsigned int GetSize() const { return (unsigned int)_angles.size(); }
};


//****************doxygen for inline functions***********
/*!
**\fn OBTorsionData::GetSize()
**\brief Gets the number of torsion structs
**\return integer count of the number of torsions
*/

/*!
**\fn OBTorsionData::GetData()
**\brief Gets a vector of the OBTorsions
**\return the vector of torsions
*/

/*!
**\fn OBAngleData::GetSize()
**\brief Gets the number of angles 
**\return integer count of the number of angles
*/

/*!
**\fn OBAngleData::GetData()
**\brief Gets the angle vector data
**\return pointer to vector of OBAngles
*/

/*!
**\fn OBAngle::SetAngle()
**\brief Sets the OBAngle angle value
**\param angle in radians
*/

/*!
**\fn OBAngle::GetAngle()
**\brief Gets the OBAngle angle value
**\return angle in radians
*/

/*!
**\fn OBTorsion::GetBondIdx()
**\brief Gets the bond index of the central bond
**\return int bond index
*/

/*!
**\fn OBTorsion::GetBC()
**\brief Gets the two central atoms of ABCD torsion
**\return pair<OBAtom*,OBAtom*>
*/

/*!
**\fn OBTorsion::GetADs()
**\brief Gets the vector of distal atoms of ABCD torsion
**\return vector of A,D atom pointers and a float
*/

} //end namespace OpenBabel

#endif // OB_GENERIC_H

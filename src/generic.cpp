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

#include "mol.h"

namespace OpenBabel {

//
//member functions for OBGenericData class
//

OBGenericData::OBGenericData()
{
  _type = obUndefinedData;
  _attr = "undefined";
}

OBGenericData::OBGenericData(const OBGenericData &src)
{
  _type = src.GetDataType();
  _attr = src.GetAttribute();
}

//
//member functions for OBCommentData class
//

OBCommentData::OBCommentData() 
{
  _type = obCommentData; 
  _attr = "Comment";
}

OBCommentData::OBCommentData(const OBCommentData &src)
{
  _type = obCommentData;
  _attr = "Comment";
  _data = src.GetData();
}

//
//member functions for OBExternalBond class
//
OBExternalBond::OBExternalBond(OBAtom *atom,OBBond *bond,int idx) 
{
  _idx = idx;
  _atom = atom;
  _bond = bond;
}

OBExternalBond::OBExternalBond(const OBExternalBond &src)
{
  _idx = src.GetIdx();
  _atom = src.GetAtom();
  _bond = src.GetBond();
}

//
//member functions for OBExternalBondData class
//

OBExternalBondData::OBExternalBondData()
{
  _type = obExternalBondData;
  _attr = "ExternalBondData";
}

void OBExternalBondData::SetData(OBAtom *atom,OBBond *bond,int idx)
{
  OBExternalBond xb(atom,bond,idx);
  _vexbnd.push_back(xb);
}

//
//member functions for OBCompressData class
//

OBCompressData::OBCompressData() 
{
  _size = 0;
  _data = (unsigned char*)NULL; 
  _type = obCompressData;
  _attr = "CompressData";
}

OBCompressData::~OBCompressData()
{
  if (_data) delete [] _data;
}

void OBCompressData::SetData(unsigned char *d,int size)
{
  if (size <= 0) return;
  
  if (_data) delete [] _data;
  
  _data = new unsigned char[size];
  memcpy(_data,(char*) d, size);
  _size = size;
}

//
//member functions for OBPairData class
//

OBPairData::OBPairData()
{
  _type = obPairData; 
  _attr = "PairData";
}

//
//member functions for OBVirtualBond class
//OBVirtualBond is used to temporarily store bonds that reference
//an atom that has not yet been added to a molecule 
//

OBVirtualBond::OBVirtualBond()
{
	_type = obVirtualBondData;
	_attr = "VirtualBondData";
	_bgn = _end = _ord = 0;
}

OBVirtualBond::OBVirtualBond(int bgn,int end,int ord,int stereo)
{
	_type = obVirtualBondData;
	_attr = "VirtualBondData";
	_bgn = bgn;
	_end = end;
	_ord = ord;
	_stereo = stereo;
}

//
//member functions for OBRingData class - stores SSSR set
//

OBRingData::OBRingData()
{
	_type = obRingData;
	_attr = "RingData";
	_vr.clear();
}

} //end namespace OpenBabel

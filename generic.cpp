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
//member functions for OEGenericData class
//

OEGenericData::OEGenericData()
{
  _type = oeUndefinedData;
  _attr = "undefined";
}

OEGenericData::OEGenericData(const OEGenericData &src)
{
  _type = src.GetDataType();
  _attr = src.GetAttribute();
}

//
//member functions for OECommentData class
//

OECommentData::OECommentData() 
{
  _type = oeCommentData; 
  _attr = "Comment";
}

OECommentData::OECommentData(const OECommentData &src)
{
  _type = oeCommentData;
  _attr = "Comment";
  _data = src.GetData();
}

//
//member functions for OEExternalBond class
//
OEExternalBond::OEExternalBond(OEAtom *atom,OEBond *bond,int idx) 
{
  _idx = idx;
  _atom = atom;
  _bond = bond;
}

OEExternalBond::OEExternalBond(const OEExternalBond &src)
{
  _idx = src.GetIdx();
  _atom = src.GetAtom();
  _bond = src.GetBond();
}

//
//member functions for OEExternalBondData class
//

OEExternalBondData::OEExternalBondData()
{
  _type = oeExternalBondData;
  _attr = "ExternalBondData";
}

void OEExternalBondData::SetData(OEAtom *atom,OEBond *bond,int idx)
{
  OEExternalBond xb(atom,bond,idx);
  _vexbnd.push_back(xb);
}

//
//member functions for OECompressData class
//

OECompressData::OECompressData() 
{
  _size = 0;
  _data = (unsigned char*)NULL; 
  _type = oeCompressData;
  _attr = "CompressData";
}

OECompressData::~OECompressData()
{
  if (_data) delete [] _data;
}

void OECompressData::SetData(unsigned char *d,int size)
{
  if (size <= 0) return;
  
  if (_data) delete [] _data;
  
  _data = new unsigned char[size];
  memcpy(_data,(char*) d, size);
  _size = size;
}

//
//member functions for OEPairData class
//

OEPairData::OEPairData()
{
  _type = oePairData; 
  _attr = "PairData";
}

//
//member functions for OEVirtualBond class
//OEVirtualBond is used to temporarily store bonds that reference
//an atom that has not yet been added to a molecule 
//

OEVirtualBond::OEVirtualBond()
{
	_type = oeVirtualBondData;
	_attr = "VirtualBondData";
	_bgn = _end = _ord = 0;
}

OEVirtualBond::OEVirtualBond(int bgn,int end,int ord,int stereo)
{
	_type = oeVirtualBondData;
	_attr = "VirtualBondData";
	_bgn = bgn;
	_end = end;
	_ord = ord;
	_stereo = stereo;
}

//
//member functions for OERingData class - stores SSSR set
//

OERingData::OERingData()
{
	_type = oeRingData;
	_attr = "RingData";
	_vr.clear();
}

} //end namespace OpenBabel

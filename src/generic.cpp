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
#include "math/matrix3x3.h"

using namespace std;

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


OBGenericData& OBGenericData::operator = (const OBGenericData &src)
{
    if(this == &src) 
        return(*this);
    
    _type = src._type;
    _attr = src._attr;
    
    return(*this);
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
  if (_data) {delete [] _data; _data = (unsigned char*)NULL;}
}

void OBCompressData::SetData(unsigned char *d,int size)
{
  if (size <= 0) return;
  
  if (_data) {delete [] _data; _data = (unsigned char*)NULL;}
  
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
// member functions for OBUnitCell class
//  used for storing information about periodic boundary conditions
//   with conversion from three-vector to a, b, c, alpha, beta, gamma
//
OBUnitCell::OBUnitCell()
{
  _a = _b = _c = _alpha = _beta = _gamma = 0.0;
  _type = obUnitCell;
  _attr = "UnitCell";
}

OBUnitCell::OBUnitCell(const OBUnitCell &src)
{
  _a = src._a;
  _b = src._b;
  _c = src._c;
  _alpha = src._alpha;
  _beta = src._beta;
  _gamma = src._gamma;
}

OBUnitCell & OBUnitCell::operator=(const OBUnitCell &src)
{
  if(this == &src) 
    return(*this);
    
  _a = src._a;
  _b = src._b;
  _c = src._c;
  _alpha = src._alpha;
  _beta = src._beta;
  _gamma = src._gamma;
    
  return(*this);
}

void OBUnitCell::SetData(const vector3 v1, const vector3 v2, const vector3 v3)
{
  _a = v1.length();
  _b = v2.length();
  _c = v3.length();
  _alpha = vectorAngle(v2, v3);
  _beta = vectorAngle(v1, v3);
  _gamma = vectorAngle(v1, v2);
}

vector<vector3> OBUnitCell::GetCellVectors()
{
  vector<vector3> v;
  vector3 cellVec;

  v.reserve(3);
  cellVec.Set(_a, 0.0, 0.0);
  v.push_back(cellVec);
  cellVec.Set(_b*cos(DEG_TO_RAD*_gamma), _b*sin(DEG_TO_RAD*_gamma), 0.0);
  v.push_back(cellVec);
  cellVec.Set(_c*cos(DEG_TO_RAD*_beta)*sin(DEG_TO_RAD*_alpha),
	      _c*sin(DEG_TO_RAD*_beta)*cos(DEG_TO_RAD*_alpha),
	      _c*sin(DEG_TO_RAD*_beta)*sin(DEG_TO_RAD*_alpha));
  v.push_back(cellVec);
  
  return v;
}

matrix3x3 OBUnitCell::GetCellMatrix()
{
  vector3 v1, v2, v3;

  v1.Set(_a, 0.0, 0.0);
  v2.Set(_b*cos(DEG_TO_RAD*_gamma), _b*sin(DEG_TO_RAD*_gamma), 0.0);
  v3.Set(_c*cos(DEG_TO_RAD*_beta)*sin(DEG_TO_RAD*_alpha),
	      _c*sin(DEG_TO_RAD*_beta)*cos(DEG_TO_RAD*_alpha),
	      _c*sin(DEG_TO_RAD*_beta)*sin(DEG_TO_RAD*_alpha));

  matrix3x3 m(v1,v2,v3);
  return m;
}

matrix3x3 OBUnitCell::GetOrthoMatrix()
{
  matrix3x3 m;
  double alphaRad, betaRad, gammaRad;
  double v;

  alphaRad = _alpha * DEG_TO_RAD;
  betaRad = _beta * DEG_TO_RAD;
  gammaRad = _gamma * DEG_TO_RAD;

  v = 1 - SQUARE(cos(alphaRad)) - SQUARE(cos(betaRad)) - SQUARE(cos(gammaRad))
    + 2 * cos(alphaRad) * cos(betaRad) * cos(gammaRad);

  m.Set(0,0, _a);
  m.Set(0,1, _b * cos(gammaRad));
  m.Set(0,2, _c * cos(betaRad));
  m.Set(1,0, 0.0);
  m.Set(1,1, _b * sin(gammaRad));
  m.Set(1,2, _c * (cos(alphaRad)-cos(betaRad)*cos(gammaRad)) / sin(gammaRad));
  m.Set(2,0, 0.0);
  m.Set(2,1, 0.0);
  m.Set(2,2, _c * v);

  return m;
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

/*!
**\brief OBRingData copy constructor
**\param src reference to original OBRingData object (rhs)
*/
OBRingData::OBRingData(const OBRingData &src)
    :	OBGenericData(src),	//chain to base class
		_vr(src._vr)				//chain to member classes
{
	//no other memeber data
	//memory management

	vector<OBRing*>::iterator ring;

	for(ring = _vr.begin();ring != _vr.end();ring++)
	{
		OBRing *newring = new OBRing;
		(*newring) = (**ring);	//copy data to new object
		(*ring)    = newring;	//repoint new pointer to new copy of data
	}
}

OBRingData::~OBRingData()
{
    vector<OBRing*>::iterator ring;
    for (ring = _vr.begin();ring != _vr.end();ring++)
    {
        delete *ring;
    }
}

/*!
**\brief OBRingdata assignment operator
**\param src reference to original OBRingData object (rhs)
**\return reference to changed OBRingData object (lhs)
*/
OBRingData& OBRingData::operator =(const OBRingData &src)
{
    //on identity, return
    if(this == &src)	return(*this);

    //chain to base class
    OBGenericData::operator =(src);

    //member data

    //memory management
    vector<OBRing*>::iterator ring;
    for(ring = _vr.begin();ring != _vr.end();ring++)
    {
	    delete &*ring;	//deallocate old rings to prevent memory leak
    }

    _vr.clear();
    _vr = src._vr;	//copy vector properties

    for(ring = _vr.begin();ring != _vr.end();ring++)
    {
	    if(*ring == 0)
		    continue;
	    
	    //allocate and copy ring data
	    OBRing *newring = new OBRing;
	    (*newring) = (**ring);
	    (*ring) = newring;	//redirect pointer
    }
    return(*this);
}

//
//member functions for OBAngle class - stores all angles
//

/*!
**\brief Angle default constructor
*/
OBAngle::OBAngle()
{
    _vertex         = 0;
    _termini.first  = 0;
    _termini.second = 0;
    _radians        = 0.0;
}

/*!
**\brief Angle constructor
*/
OBAngle::OBAngle(OBAtom *vertex,OBAtom *a,OBAtom *b)
{
	_vertex         = vertex;
	_termini.first  = a; 
    _termini.second = b ;

	SortByIndex();
}

/*!
**\brief OBAngle copy constructor
*/
OBAngle::OBAngle(const OBAngle &src)
	:	_termini(src._termini)
{
	_vertex  = src._vertex;
	_radians = src._radians;
}

/*!
**\brief OBAngle assignment operator
*/
OBAngle& OBAngle::operator = (const OBAngle &src)
{
	if (this == &src) 
        return(*this);

	_vertex         = src._vertex;
	_termini.first  = src._termini.first;
	_termini.second = src._termini.second;
	_radians        = src._radians;

	return(*this);
}

/*!
**\brief Return OBAngle to its original state
*/
void OBAngle::Clear()
{
	_vertex         = 0;
	_termini.first  = 0;
    _termini.second = 0;
	_radians        = 0.0;
	return;
}

/*!
**\brief Sets the 3 atoms in the angle
**\param pointers to each OBAtom
*/
void OBAngle::SetAtoms(OBAtom *vertex,OBAtom *a,OBAtom *b)
{
	_vertex         = vertex;
	_termini.first  = a; 
    _termini.second = b;
	SortByIndex();
	return;
}

/*!
**\brief Sets the 3 atoms in the angle
**\param atoms a triple of OBAtom pointers, the first must be the vertex
*/
void OBAngle::SetAtoms(triple<OBAtom*,OBAtom*,OBAtom*> &atoms)
{
	_vertex         = atoms.first; 
	_termini.first  = atoms.second;
	_termini.second = atoms.third;
	SortByIndex();
	return;
}

/*!
**\brief Retrieves the 3 atom pointer for the angle (vertex first)
**\return triple of OBAtom pointers 
*/
triple<OBAtom*,OBAtom*,OBAtom*> OBAngle::GetAtoms()
{
	triple<OBAtom*,OBAtom*,OBAtom*> atoms;
	atoms.first  = _vertex;
	atoms.second = _termini.first;
	atoms.third  = _termini.second;
    return(atoms);
}

/*!
**\brief sorts atoms in angle by order of indices
*/
void OBAngle::SortByIndex()
{
	OBAtom *tmp;

	if(_termini.first->GetIdx() > _termini.second->GetIdx())
	{
		tmp             = _termini.first; 
		_termini.first  = _termini.second;
		_termini.second = tmp;
	}
}

/*!
**\brief OBAngle equality operator, is same angle, NOT same value
**\return boolean equality
*/
bool OBAngle::operator ==(const OBAngle &other)
{
	return ((_vertex         == other._vertex)        &&
            (_termini.first  == other._termini.first) &&
            (_termini.second == other._termini.second));
}

//
//member functions for OBAngleData class - stores OBAngle set
//

/*!
**\brief OBAngleData constructor
*/
OBAngleData::OBAngleData()
	:	OBGenericData()
{
    _type = obAngleData;
    _attr = "AngleData";
}

/*!
**\brief OBAngleData copy constructor
*/
OBAngleData::OBAngleData(const OBAngleData &src)
	:	OBGenericData(src), _angles(src._angles)
{
    _type = obAngleData;
    _attr = "AngleData";
}

/*!
**\brief OBAngleData assignment operator
*/
OBAngleData& OBAngleData::operator =(const OBAngleData &src)
{
	if (this == &src)
        return(*this);

    _angles = src._angles;

    return(*this);
}

/*!
**\brief sets OBAngleData to its original state
*/
void OBAngleData::Clear()
{
	_angles.clear();
	return;
}

/*!
**\brief Adds a new angle to OBAngleData
*/
void OBAngleData::SetData(OBAngle &angle)
{
	_angles.push_back(angle);
	return;
}

/*!
**\brief Fills an array with the indices of the atoms in the angle (vertex first)
**\param angles pointer to the pointer to an array of angles atom indices
**\param size the current number of rows in the array
**\return int The number of angles
*/
unsigned int OBAngleData::FillAngleArray(int **angles, unsigned int &size)
{
    if(_angles.size() > size)
    {
        delete [] *angles;
        *angles = new int[_angles.size()*3];
        size    = (unsigned int)_angles.size();
    }

    vector<OBAngle>::iterator angle;
    int angleIdx = 0;
    for( angle=_angles.begin(); angle!=_angles.end(); angle++)
    {
        *angles[angleIdx++] = angle->_vertex->GetIdx();
        *angles[angleIdx++] = angle->_termini.first->GetIdx();
        *angles[angleIdx++] = angle->_termini.second->GetIdx();
    }
    return (unsigned int)_angles.size();
}

//
//member functions for OBAngleData class - stores OBAngle set
//

/*!
**\brief OBTorsion constructor
*/
OBTorsion::OBTorsion(OBAtom *a,OBAtom *b, OBAtom *c,OBAtom *d)
{
    triple<OBAtom*,OBAtom*,double> ad(a,d,0.0);
    _ads.push_back(ad);

    _bc.first  = b;
    _bc.second = c;
}

/*!
**\brief OBTorsion copy constructor
*/
OBTorsion::OBTorsion(const OBTorsion &src)
	:	_bc(src._bc), _ads(src._ads)
{
}

/*!
**\brief Returns all the 4 atom sets in OBTorsion
*/
vector<quad<OBAtom*,OBAtom*,OBAtom*,OBAtom*> > OBTorsion::GetTorsions()
{
    quad<OBAtom*,OBAtom*,OBAtom*,OBAtom*> abcd;

    abcd.second = _bc.first;
    abcd.third  = _bc.second;

    vector<quad<OBAtom*,OBAtom*,OBAtom*,OBAtom*> > torsions;
    vector<triple<OBAtom*,OBAtom*,double> >::iterator ad;

    for(ad = _ads.begin();ad != _ads.end();ad++)
    {
        abcd.first = ad->first;
        abcd.fourth = ad->second;
        torsions.push_back(abcd);
    }

    return(torsions);
}

/*!
**\brief OBTorsion assignment operator
*/
OBTorsion& OBTorsion::operator =(const OBTorsion &src)
{
	if (this == &src) 
        return(*this);

    _bc  = src._bc;
    _ads = src._ads;

    return(*this);
}

/*!
**\brief Returns the OBTorsion to its original state
*/
void OBTorsion::Clear()
{
    _bc.first  = 0;
    _bc.second = 0;
    _ads.erase(_ads.begin(),_ads.end());
}

/*!
**\brief Sets the angle of a torsion in OBTorsion
**\param radians the value to assign to the torsion
**\param index the index into the torsion of the OBTorsion
**\return boolean success
*/
bool OBTorsion::SetAngle(double radians,unsigned int index)
{
    if(index >= _ads.size())
        return(false);

    _ads[index].third = radians;

    return(true);
}

/*!
**\brief Obtains the angle of a torsion in OBTorsion
**\param radians the value of the angle is set here
**\param index the index into the torsion of the OBTorsion
**\return boolean success
*/
bool OBTorsion::GetAngle(double &radians, unsigned int index)
{
    if(index >= _ads.size())
        return false;
    radians = _ads[index].third;
    return true;
}

unsigned int OBTorsion::GetBondIdx()
{
    return(_bc.first->GetBond(_bc.second)->GetIdx());
}

/*!
**\brief determines if torsion has only protons on either the a or d end
**\return boolean 
*/
bool OBTorsion::IsProtonRotor()
{
    bool Aprotor = true;
    bool Dprotor = true;
    vector<triple<OBAtom*,OBAtom*,double> >::iterator ad;
    for(ad = _ads.begin();ad != _ads.end() && (Aprotor || Dprotor);ad++)
    {
        if(!ad->first->IsHydrogen())
            Aprotor = false;
        if(!ad->second->IsHydrogen())
            Dprotor = false;
    }
    return (Aprotor || Dprotor);
}

/*!
**\brief adds a new torsion to the OBTorsion object
*/
bool OBTorsion::AddTorsion(OBAtom *a,OBAtom *b, OBAtom *c,OBAtom *d)
{
    if(!Empty() && (b != _bc.first || c != _bc.second))
        return(false);

    if(Empty())
    {
        _bc.first  = b;
        _bc.second = c;
    }

    triple<OBAtom*,OBAtom*,double> ad(a,d,0.0);
    _ads.push_back(ad);

    return(true);
}

/*!
**\brief adds a new torsion to the OBTorsion object
*/
bool OBTorsion::AddTorsion(quad<OBAtom*,OBAtom*,OBAtom*,OBAtom*> &atoms)
{
    if(!Empty() && (atoms.second != _bc.first || atoms.third != _bc.second))
        return(false);

    if(Empty())
    {
        _bc.first  = atoms.second;
        _bc.second = atoms.third;
    }

    triple<OBAtom*,OBAtom*,double> ad(atoms.first,atoms.fourth,0.0);
    _ads.push_back(ad);

    return(true);
}

//\!brief OBTorsionData ctor
OBTorsionData::OBTorsionData()
{
    _type = obTorsionData;
    _attr = "TorsionData";
}

//
//member functions for OBTorsionData class - stores OBTorsion set
//
OBTorsionData::OBTorsionData(const OBTorsionData &src)
	:	OBGenericData(src), _torsions(src._torsions)
{
    _type = obTorsionData;
    _attr = "TorsionData";
}

OBTorsionData& OBTorsionData::operator =(const OBTorsionData &src)
{
    if (this == &src) 
        return(*this);

    OBGenericData::operator =(src);

    _type     = obTorsionData;
    _attr     = "TorsionData";
    _torsions = src._torsions;

    return(*this);
}

void OBTorsionData::Clear()
{
    _torsions.clear();
}

void OBTorsionData::SetData(OBTorsion &torsion)
{
    _torsions.push_back(torsion);
}

/*!
**\brief Fills a vector with the indices of the atoms in torsions (ordered abcd)
**\param torsions reference to the vector of abcd atom sets
**\return boolean success
*/
bool OBTorsionData::FillTorsionArray(vector<vector<unsigned int> > &torsions)
{
    if(_torsions.size() == 0)
        return(false);

    vector<quad<OBAtom*,OBAtom*,OBAtom*,OBAtom*> > tmpquads,quads;
    vector<quad<OBAtom*,OBAtom*,OBAtom*,OBAtom*> >::iterator thisQuad;
    vector<OBTorsion>::iterator torsion;

    //generate set of all 4 atom abcd's from torsion structure
    for (torsion = _torsions.begin();torsion != _torsions.end();torsion++)
    {
        tmpquads = torsion->GetTorsions();
        for(thisQuad = tmpquads.begin();thisQuad != tmpquads.end();thisQuad++)
	        quads.push_back(*thisQuad);
    }

    //fill array of torsion atoms

    torsions.clear();
    torsions.resize(quads.size());

    unsigned int ct = 0;

    for (thisQuad = quads.begin();thisQuad != quads.end();thisQuad++,ct++)
    {
        torsions[ct].resize(4);
        torsions[ct][0] = thisQuad->first->GetIdx()-1;
        torsions[ct][1] = thisQuad->second->GetIdx()-1;
        torsions[ct][2] = thisQuad->third->GetIdx()-1;
        torsions[ct][3] = thisQuad->fourth->GetIdx()-1;
    }

    return(true);
}

} //end namespace OpenBabel

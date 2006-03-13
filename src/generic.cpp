/**********************************************************************
generic.cpp - Handle OBGenericData classes.
 
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2005 by Geoffrey R. Hutchison
 
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

#include "mol.h"
#include "math/matrix3x3.h"

using namespace std;

namespace OpenBabel
{

  /** \class OBGenericData

OBGenericData is an abstract base class which defines an interface for
storage, retrieval, and indexing of arbitrary generic data.
Subclasses of OBGenericData can be used to store custom data
on a per-atom, per-bond, per-molecule, or per-residue basis.
Open Babel currently supports a small subset of chemical functionality
as OBGenericData types, which will expand over time to support additional
interconversion (e.g., spectroscopy, dynamics, surfaces...)

For your own custom data, either define a custom subclass using 
an id from the OBGenericDataType::CustomData0 to OBGenericDataType::CustomData15 slots,
or store your data as a string and use OBPairData for key/value access.
The latter is <strong>highly</strong> recommended for various text descriptors
e.g., in QSAR, atom or bond labels, or other textual data.

Example code using OBGenericData

@code
    if (mol.HasData(OBGenericDataType::UnitCell))
    {
        uc = (OBUnitCell*)mol.GetData(OBGenericDataType::UnitCell);
        sprintf(buffer,
                "%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f",
                uc->GetA(), uc->GetB(), uc->GetC(),
                uc->GetAlpha() , uc->GetBeta(), uc->GetGamma());
        ofs << buffer << endl;
    }

...

  vector<OBGenericData*>::iterator k;
  vector<OBGenericData*> vdata = mol.GetData();
  for (k = vdata.begin();k != vdata.end();k++)
     if ((*k)->GetDataType() == OBGenericDataType::PairData)
          {
                  ofs << ">  <" << (*k)->GetAttribute() << ">" << endl;
                  ofs << ((OBPairData*)(*k))->GetValue() << endl << endl;
          }
@endcode

Similar code also works for OBGenericData stored in an OBAtom or OBBond (or OBResidue).

@code
    if (!atom.HasData("UserLabel")) // stored textual data as an OBPairData
      {
         OBPairData *label = new OBPairData;
	 label->SetAttribute("UserLabel");
	 label->SetValue(userInput);

	 atom.SetData(label);
      }

...

    if (bond.HasData("DisplayType")) // e.g. in a visualization tool
      {
        OBPairData *display = dynamic_cast<OBPairData *> bond.GetData("DisplayType");
	if (display->GetValue() == "wireframe")
	  {
	   ... // display a wireframe view
	  }
      }
@endcode

When designing a class derived from OBGenericData you must add a Clone() function.
For classes used with OBMol this is used when an OBMol object is copied. If your
class member variables contain pointers to atoms or bonds then it will be necessary
to ensure that these are updated in Clone() to refer to the new molecule. Without
these and similar pointers it is more likely that the very simple clone function
@code
		virtual OBGenericData* Clone(OBBase* parent) const{return new MyNewClass(*this);}
@endcode
and the compiler generated copy constructor would be sufficient. It is recommended
that, if possible, OBGenericData classes do not store atom and bond pointers. Using
atom and bond indices instead would allow the simple version of Clone() above. See 
OBRotameterData::Clone for an example of a more complicated version. For classes
which are not intended to support copying, Clone can return NULL 
@code
		virtual OBGenericData* Clone(OBBase* parent) const{return NULL;}
@endcode
Clone is a pure virtual function so that you need to decide what kind of
function you need and include it explicitly.
**/

//
//member functions for OBGenericData class
//

OBGenericData::OBGenericData()
{
    _type = OBGenericDataType::UndefinedData;
    _attr = "undefined";
}

/* Use default copy constructor and assignment operators
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
*/

//
//member functions for OBCommentData class
//

OBCommentData::OBCommentData()
{
    _type = OBGenericDataType::CommentData;
    _attr = "Comment";
}

OBCommentData::OBCommentData(const OBCommentData &src) :
  OBGenericData()
{
    _type = OBGenericDataType::CommentData;
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
    _type = OBGenericDataType::ExternalBondData;
    _attr = "ExternalBondData";
}

void OBExternalBondData::SetData(OBAtom *atom,OBBond *bond,int idx)
{
    OBExternalBond xb(atom,bond,idx);
    _vexbnd.push_back(xb);
}

//
//member functions for OBPairData class
//

OBPairData::OBPairData()
{
    _type = OBGenericDataType::PairData;
    _attr = "PairData";
}

//
//member functions for OBVirtualBond class
//

OBVirtualBond::OBVirtualBond()
{
    _type = OBGenericDataType::VirtualBondData;
    _attr = "VirtualBondData";
    _bgn = _end = _ord = 0;
}

OBVirtualBond::OBVirtualBond(int bgn,int end,int ord,int stereo)
{
    _type = OBGenericDataType::VirtualBondData;
    _attr = "VirtualBondData";
    _bgn = bgn;
    _end = end;
    _ord = ord;
    _stereo = stereo;
}

//
// member functions for OBUnitCell class
//
OBUnitCell::OBUnitCell()
{
    _a = _b = _c = _alpha = _beta = _gamma = 0.0;
    _type = OBGenericDataType::UnitCell;
    _attr = "UnitCell";
}

OBUnitCell::OBUnitCell(const OBUnitCell &src) :
  OBGenericData()
{
    _a = src._a;
    _b = src._b;
    _c = src._c;
    _alpha = src._alpha;
    _beta = src._beta;
    _gamma = src._gamma;
    _offset = src._offset;

    _v1 = src._v1;
    _v2 = src._v2;
    _v3 = src._v3;

    _spaceGroup = src._spaceGroup;
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
    _offset = src._offset;

    _v1 = src._v1;
    _v2 = src._v2;
    _v3 = src._v3;

    _spaceGroup = src._spaceGroup;

    return(*this);
}

//! Implements <a href="http://qsar.sourceforge.net/dicts/blue-obelisk/index.xhtml#convertCartesianIntoNotionalCoordinates">blue-obelisk:convertCartesianIntoNotionalCoordinates</a>
void OBUnitCell::SetData(const vector3 v1, const vector3 v2, const vector3 v3)
{
    _a = v1.length();
    _b = v2.length();
    _c = v3.length();
    _alpha = vectorAngle(v2, v3);
    _beta = vectorAngle(v1, v3);
    _gamma = vectorAngle(v1, v2);
    _v1 = v1;
    _v2 = v2;
    _v3 = v3;
}

//! Implements <a href="http://qsar.sourceforge.net/dicts/blue-obelisk/index.xhtml#convertNotionalIntoCartesianCoordinates">blue-obelisk:convertNotionalIntoCartesianCoordinates</a>
vector<vector3> OBUnitCell::GetCellVectors()
{
    vector<vector3> v;
    v.reserve(3);

    if (_v1.length() == 0 && _v2.length() == 0 && _v3.length() == 0)
      {
	vector3 temp;
	matrix3x3 m = GetOrthoMatrix();

	temp = vector3(1.0f, 0.0f, 0.0f);
	v.push_back(m * temp);
	temp = vector3(0.0f, 1.0f, 0.0f);
	v.push_back(m * temp);
	temp = vector3(0.0f, 0.0f, 1.0f);
	v.push_back(m * temp);
      }
    else
      {
	v.push_back(_v1);
	v.push_back(_v2);
	v.push_back(_v3);
      }

    return v;
}

matrix3x3 OBUnitCell::GetCellMatrix()
{
    matrix3x3 m;

    if (_v1.length() == 0 && _v2.length() == 0 && _v3.length() == 0)
      {
	m = GetOrthoMatrix();
      }
    else
      {
	vector3 v1, v2, v3;
	v1 = _v1;
	v2 = _v2;
	v3 = _v3;
	m = matrix3x3(v1,v2,v3);
      }
    return m;
}

//! Implements <a href="http://qsar.sourceforge.net/dicts/blue-obelisk/index.xhtml#calculateOrthogonalisationMatrix">blue-obelisk:calculateOrthogonalisationMatrix</a>
matrix3x3 OBUnitCell::GetOrthoMatrix()
{
  matrix3x3 m;
  
  // already here, let's not duplicate the work
  m.FillOrth(_alpha, _beta, _gamma, _a, _b, _c);

  return m;
}

// Based on code in PyMMLib: http://pymmlib.sf.net/
//! Matrix to convert from Cartesian to fractional
//! Implements <a href="http://qsar.sourceforge.net/dicts/blue-obelisk/index.xhtml#convertCartesianIntoFractionalCoordinates">blue-obelisk:convertCartesianIntoFractionalCoordinates</a> 
matrix3x3 OBUnitCell::GetFractionalMatrix()
{
  matrix3x3 m;
  double sinAlpha, sinBeta, sinGamma;
  double cosAlpha, cosBeta, cosGamma;
  double v;

  sinAlpha = sin(_alpha * DEG_TO_RAD);
  sinBeta = sin(_beta * DEG_TO_RAD);
  sinGamma = sin(_gamma * DEG_TO_RAD);
  cosAlpha = cos(_alpha * DEG_TO_RAD);
  cosBeta = cos(_beta * DEG_TO_RAD);
  cosGamma = cos(_gamma * DEG_TO_RAD);

  v = sqrt(1 - SQUARE(cosAlpha) - SQUARE(cosBeta) - SQUARE(cosGamma) +
	   2 * cosAlpha*cosBeta*cosGamma);

  m.Set(0,0, 1.0f / _a);
  m.Set(0,1, -cosGamma / (_a * sinGamma) );
  m.Set(0,2, (cosGamma * cosAlpha - cosBeta) / (_a * v * sinGamma) );
  m.Set(1,0, 0.0);
  m.Set(1,1, 1.0f / (_b * sinGamma) );
  m.Set(1,2, (cosGamma * cosBeta - cosAlpha) / (_b * v * sinGamma) );
  m.Set(2,0, 0.0);
  m.Set(2,1, 0.0);
  m.Set(2,2, sinGamma / (_c * v) );

  return m;
}

//
// member functions for OBSymmetryData class
//
OBSymmetryData::OBSymmetryData()
{
    _type = OBGenericDataType::SymmetryData;
    _attr = "Symmetry";
}

OBSymmetryData::OBSymmetryData(const OBSymmetryData &src) :
  OBGenericData()
{
  _pointGroup = src._pointGroup;
  _spaceGroup = src._spaceGroup;
}

OBSymmetryData & OBSymmetryData::operator=(const OBSymmetryData &src)
{
    if(this == &src)
        return(*this);

    _pointGroup = src._pointGroup;
    _spaceGroup = src._spaceGroup;

    return(*this);
}

OBConformerData::OBConformerData()
{
    _type = OBGenericDataType::ConformerData;
    _attr = "Conformers";
}

OBConformerData::OBConformerData(const OBConformerData &src) :
  OBGenericData()
{
  _vDimension = src._vDimension;
  _vEnergies = src._vEnergies;
  _vForces = src._vForces;
  _vVelocity = src._vVelocity;
  _vDisplace = src._vDisplace;
  _vData = src._vData;
}

OBConformerData & OBConformerData::operator=(const OBConformerData &src)
{
    if(this == &src)
        return(*this);

  _vDimension = src._vDimension;
  _vEnergies = src._vEnergies;
  _vForces = src._vForces;
  _vVelocity = src._vVelocity;
  _vDisplace = src._vDisplace;
  _vData = src._vData;

  return(*this);
}

//
//member functions for OBRingData class
//

OBRingData::OBRingData()
{
    _type = OBGenericDataType::RingData;
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
**\brief OBRingData assignment operator
**\param src reference to original OBRingData object (rhs)
**\return reference to changed OBRingData object (lhs)
*/
OBRingData& OBRingData::operator =(const OBRingData &src)
{
    //on identity, return
    if(this == &src)
        return(*this);

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
** Parameters are pointers to each OBAtom
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
    _type = OBGenericDataType::AngleData;
    _attr = "AngleData";
}

/*!
**\brief OBAngleData copy constructor
*/
OBAngleData::OBAngleData(const OBAngleData &src)
        :	OBGenericData(src), _angles(src._angles)
{
    _type = OBGenericDataType::AngleData;
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
{}

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
    _type = OBGenericDataType::TorsionData;
    _attr = "TorsionData";
}

//
//member functions for OBTorsionData class - stores OBTorsion set
//
OBTorsionData::OBTorsionData(const OBTorsionData &src)
        :	OBGenericData(src), _torsions(src._torsions)
{
    _type = OBGenericDataType::TorsionData;
    _attr = "TorsionData";
}

OBTorsionData& OBTorsionData::operator =(const OBTorsionData &src)
{
    if (this == &src)
        return(*this);

    OBGenericData::operator =(src);

    _type     = OBGenericDataType::TorsionData;
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
bool OBTorsionData::FillTorsionArray(std::vector<std::vector<unsigned int> > &torsions)
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

//
// Member functions for OBChiralDarta
//
bool OBChiralData::SetAtom4Refs(std::vector<unsigned int> atom4refs, atomreftype t)
{
     if (atom4refs.size()>4)
       {
	 obErrorLog.ThrowError(__FUNCTION__, "More than 4 atoms in atom4refs", obDebug);
	 return(false);
       }
     switch(t){
               case input: _atom4refs = atom4refs;break;
               case output:_atom4refo = atom4refs;break;
               case calcvolume:_atom4refc = atom4refs;break;
               default: 
		 obErrorLog.ThrowError(__FUNCTION__, "AtomRefType called is invalid", obDebug);
		 return(false);
               }
     return (true);
}
int OBChiralData::AddAtomRef(unsigned int atomref, atomreftype t)
{
    switch(t){
              case input: _atom4refs.push_back(atomref);break;
              case output: _atom4refo.push_back(atomref);break;
              case calcvolume:_atom4refc.push_back(atomref);break;
              default:
		 obErrorLog.ThrowError(__FUNCTION__, "AtomRefType called is invalid", obDebug);
		 return(false);
              }
              
    return (_atom4refs.size());
}

unsigned int OBChiralData::GetAtomRef(int a, atomreftype t)
{
        switch(t){
             case input: return(_atom4refs[a]);break;
             case output: return(_atom4refo[a]);break;
             case calcvolume: return(_atom4refc[a]);break;
             default:
		 obErrorLog.ThrowError(__FUNCTION__, "AtomRefType called is invalid", obDebug);
		 return(false);
              }  
}
std::vector<unsigned int> OBChiralData::GetAtom4Refs(atomreftype t) const
{
            switch (t){
               case output:
                return(_atom4refo);
                break;
               case input:
                return(_atom4refs);
                break;
               case calcvolume:
                return(_atom4refc);
                break;
               default:
		 obErrorLog.ThrowError(__FUNCTION__, "AtomRefType called is invalid", obDebug);
		 return(_atom4refo);
	    }
}

unsigned int OBChiralData::GetSize(atomreftype t) const
{
  switch (t)
    {
    case output:
      return(unsigned int)_atom4refo.size();
      break;
    case input:
      return(unsigned int)_atom4refs.size();
      break;
    case calcvolume:
      return(unsigned int)_atom4refc.size();
    default:
		 obErrorLog.ThrowError(__FUNCTION__, "AtomRefType called is invalid", obDebug);
		 return(0);
    }
}

OBChiralData::OBChiralData()
{
    _type = OBGenericDataType::ChiralData;
    _attr = "ChiralData";
}

OBChiralData::OBChiralData(const OBChiralData &src)
  : OBGenericData()
{
  _atom4refs = src._atom4refs;
  _atom4refo = src._atom4refo;
  _atom4refc = src._atom4refc;
  parity=src.parity;
}

OBChiralData & OBChiralData::operator=(const OBChiralData &src)
{
    if(this == &src)
        return(*this);

    _atom4refs = src._atom4refs;
    _atom4refo = src._atom4refo;
    _atom4refc = src._atom4refc;
    parity=src.parity;
    return(*this);
}

void OBChiralData::Clear()
{
     _atom4refs.clear();
    parity=0;
    _atom4refo.clear();
    _atom4refc.clear();
}

} //end namespace OpenBabel

//! \file generic.cpp
//! \brief Handle OBGenericData classes. Custom data for atoms, bonds, etc.

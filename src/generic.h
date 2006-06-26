/**********************************************************************
generic.h - Handle generic data classes. Custom data for atoms, bonds, etc.
 
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

#ifndef OB_GENERIC_H
#define OB_GENERIC_H

#include "babelconfig.h"

#include <string>
#include <vector>
#include <map>

#include "math/vector3.h"
#include "obutil.h"

namespace OpenBabel
{

class OBBase;
class OBAtom;
class OBBond;
class OBRing;

OBAPI std::string& Trim(std::string& txt);

//! \brief Classification of data stored via OBGenericData class and subclasses.
//!
//! OBGenericDataType can be used as a faster, direct access to a particular category
//! instead of the slower access via GetData(std::string), which must loop
//! through all data to find a match with the supplied key. It is implemented
//! as a set of unsigned integer constants for maximum flexibility and future
//! expansion.
//! 
//! CustomData0 through CustomData15 are data slots that are not used in 
//! OpenBabel directly and are meant for use in derivative programs.
//! Macro definitions can be used to define what each data slot is used in your code.
namespace OBGenericDataType
{
  //! Unknown data type (default)
  static const unsigned int UndefinedData =      0;

  //! Arbitrary key/value data, i.e., OBPairData
  static const unsigned int PairData      =      1;

  //! Energetics data (e.g., total energy, heat of formation, etc.)
  static const unsigned int EnergyData    =      2;

  //! Storing text comments (one per molecule, atom, bond, etc.) (for other data, e.g., author, keyword, ... use OBPairData)
  static const unsigned int CommentData   =      3;

  //! Arbitrary information about conformers, i.e., OBConformerData
  static const unsigned int ConformerData =      4;

  //! Bond data external to OpenBabel, i.e., OBExternalBond, OBExternalBondData
  static const unsigned int ExternalBondData =   5;

  //! Information for generating & manipulating rotamers, i.e. OBRotamerList
  static const unsigned int RotamerList =        6;

  //! Info. for storing bonds to atoms yet to be added, i.e. OBVirtualBond
  static const unsigned int VirtualBondData =    7;

  //! Information on rings in a molecule, i.e., OBRingData
  static const unsigned int RingData =           8;

  //! Information about torsion/dihedral angles, i.e., OBTorsionData and OBTorsion
  static const unsigned int TorsionData =        9;

  //! Bond angles in a molecule, i.e., OBAngle, OBAngleData
  static const unsigned int AngleData =         10;

  //! Residue serial numbers
  static const unsigned int SerialNums =        11;

  //! Crystallographic unit cell data, i.e., OBUnitCell
  static const unsigned int UnitCell =          12;

  //! Spin data, including NMR, atomic and molecular spin, etc.
  static const unsigned int SpinData =          13;

  //! Arbitrary partial and total charges, dipole moments, etc.
  static const unsigned int ChargeData =        14;

  //! Symmetry data -- point and space groups, transforms, etc. i.e., OBSymmetryData
  static const unsigned int SymmetryData =      15;

  //! Arbitrary chiral information (atom, bond, molecule, etc.) i.e., OBChiralData
  static const unsigned int ChiralData =        16;

  //! Atomic and molecular occupation data
  static const unsigned int OccupationData =    17;

  //! Density (cube) data and surfaces
  static const unsigned int DensityData =       18;

  //! Electronic levels, redox states, orbitals, etc.
  static const unsigned int ElectronicData =    19;

  //! Vibrational modes, frequencies, etc.
  static const unsigned int VibrationData =     20;

  //! Rotational energy information
  static const unsigned int RotationData =      21;

  //! Nuclear transitions (e.g., decay, fission, fusion)
  static const unsigned int NuclearData =       22;

  // space for up to 2^14 more entries...

  //! Custom (user-defined data)
  static const unsigned int CustomData0 = 16384;
  static const unsigned int CustomData1 = 16385;
  static const unsigned int CustomData2 = 16386;
  static const unsigned int CustomData3 = 16387;
  static const unsigned int CustomData4 = 16388;
  static const unsigned int CustomData5 = 16389;
  static const unsigned int CustomData6 = 16390;
  static const unsigned int CustomData7 = 16391;
  static const unsigned int CustomData8 = 16392;
  static const unsigned int CustomData9 = 16393;
  static const unsigned int CustomData10 = 16394;
  static const unsigned int CustomData11 = 16395;
  static const unsigned int CustomData12 = 16396;
  static const unsigned int CustomData13 = 16397;
  static const unsigned int CustomData14 = 16398;
  static const unsigned int CustomData15 = 16399;
} // end namespace

//! \brief Base class for generic data
// class introduction in generic.cpp
class OBAPI OBGenericData
{
protected:
    std::string _attr; //!< attribute tag (e.g., "UnitCell", "Comment" or "Author")
    unsigned    _type; //!< attribute type -- declared for each subclass
public:
    OBGenericData();
    //Use default copy constructor and assignment operators
		//OBGenericData(const OBGenericData&);
		
		//Virtual constructors added see http://www.parashift.com/c++-faq-lite/abcs.html#faq-22.5
		//to allow copying given only a base class OBGenericData pointer.
		//It may be necessary to cast the return pointer to the derived class type since we are doing
		//without Covariant Return Types http://www.parashift.com/c++-faq-lite/virtual-functions.html#faq-20.8
		//A derived class may return NULL if copying inappropriate
		virtual OBGenericData* Clone(OBBase* parent) const
      { return NULL; } 
    virtual ~OBGenericData()    {}
    //Use default copy constructor and assignment operators
		//OBGenericData& operator=(const OBGenericData &src);

    //! Set the attribute (key), which can be used to retrieve this data
    void                      SetAttribute(const std::string &v)
    {        _attr = v;        }
    //! \return the attribute (key), which can be used to retrieve this data
    virtual const std::string &GetAttribute()  const
    {        return(_attr);    }
    //! \return the data type for this object as defined in OBGenericDataType
    unsigned int                GetDataType()    const
    {        return(_type);    }
		//! Base class returns value but should never be called
    virtual const std::string &GetValue()  const
		{			return _attr; }
};

//! Used to store a comment string (can be multiple lines long)
class OBAPI OBCommentData : public OBGenericData
{
protected:
    std::string _data;
public:
    OBCommentData();
    OBCommentData(const OBCommentData&);
		virtual OBGenericData* Clone(OBBase* parent) const{return new OBCommentData(*this);}
		
		OBCommentData& operator=(const OBCommentData &src);

    void          SetData(const std::string &data)
    { _data = data; Trim(_data); }
    void          SetData(const char *d)
    {_data = d; Trim(_data);     }
    const std::string &GetData()              const
    {        return(_data);      }
    virtual const std::string &GetValue()              const  
    {        return(_data);      }
};

//! \brief Used to store information on an external bond 
//! (e.g., SMILES fragments)
class OBAPI OBExternalBond
{
    int     _idx;
    OBAtom *_atom;
    OBBond *_bond;
public:
    OBExternalBond()    {}
    OBExternalBond(OBAtom *,OBBond *,int);
    OBExternalBond(const OBExternalBond &);
    ~OBExternalBond()   {}

    int     GetIdx()  const    {        return(_idx);    }
    OBAtom *GetAtom() const    {        return(_atom);   }
    OBBond *GetBond() const    {        return(_bond);   }
    void SetIdx(int idx)       {        _idx = idx;      }
    void SetAtom(OBAtom *atom) {        _atom = atom;    }
    void SetBond(OBBond *bond) {        _bond = bond;    }
};

//! \brief Used to store information on external bonds (e.g., in SMILES fragments)
class OBAPI OBExternalBondData : public OBGenericData
{
protected:
    std::vector<OBExternalBond> _vexbnd;
public:
    OBExternalBondData();
		
		//Copying is not used and too much work to set up
		virtual OBGenericData* Clone(OBBase* parent) const{return NULL;}
    
		void SetData(OBAtom*,OBBond*,int);
    std::vector<OBExternalBond> *GetData()
    {
        return(&_vexbnd);
    }
};

//! \brief Used to store arbitrary attribute/value relationships.
//!
//! Ideal for arbitrary text descriptors for molecules, atoms, bonds, residues,
//!  e.g. in QSAR.
class OBAPI OBPairData : public OBGenericData
{
protected:
    std::string _value;
public:
    OBPairData();
		virtual OBGenericData* Clone(OBBase* parent) const{return new OBPairData(*this);}
    void    SetValue(const char *v)
    {
        _value = v;
    }
    void    SetValue(const std::string &v)
    {
        _value = v;
    }
    virtual const std::string &GetValue() const //now virtual and const
    {
        return(_value);
    }
};

//! \brief Used to temporarily store bonds that reference
//! an atom that has not yet been added to a molecule
class OBAPI OBVirtualBond : public OBGenericData
{
protected:
    int _bgn;
    int _end;
    int _ord;
    int _stereo;
public:
    OBVirtualBond();
		virtual OBGenericData* Clone(OBBase* parent) const{return new OBVirtualBond(*this);}
    OBVirtualBond(int,int,int,int stereo=0);
    int GetBgn()
    {
        return(_bgn);
    }
    int GetEnd()
    {
        return(_end);
    }
    int GetOrder()
    {
        return(_ord);
    }
    int GetStereo()
    {
        return(_stereo);
    }
};

//! Used to store the SSSR set (filled in by OBMol::GetSSSR())
class OBAPI OBRingData : public OBGenericData
{
protected:
    std::vector<OBRing*> _vr;
public:
    OBRingData();
    OBRingData(const OBRingData &);
		virtual OBGenericData* Clone(OBBase* parent) const{return new OBRingData(*this);}
    ~OBRingData();

    OBRingData &operator=(const OBRingData &);

    void SetData(std::vector<OBRing*> &vr)
    {
        _vr = vr;
    }
    void PushBack(OBRing *r)
    {
        _vr.push_back(r);
    }
    std::vector<OBRing*> &GetData()
    {
        return(_vr);
    }
};

//! \brief Used for storing information about periodic boundary conditions
//!   with conversion to/from translation vectors and
//!  (a, b, c, alpha, beta, gamma)
class OBAPI OBUnitCell: public OBGenericData
{
 public:
    enum LatticeType { Undefined, 
                       Triclinic, Monoclinic, Orthorhombic, Tetragonal, 
                       Rhombohedral, Hexagonal, Cubic};

 protected:
    double _a, _b, _c, _alpha, _beta, _gamma;
    vector3 _offset; //!< offset for origin
    vector3 _v1, _v2, _v3; //!< translation vectors
    std::string _spaceGroup;
    LatticeType _lattice;
public:
	//! public contructor
    OBUnitCell();
    OBUnitCell(const OBUnitCell &);
		virtual OBGenericData* Clone(OBBase* parent) const
      {return new OBUnitCell(*this);}
    ~OBUnitCell()    {}

    OBUnitCell &operator=(const OBUnitCell &);

	/*!
	 **\brief Sets the vectors and angles of the unitcell
	 **\param a The length a
	 **\param b The length b
	 **\param c The length c
	 **\param alpha The angle alpha
	 **\param beta The angle beta
	 **\param gamma The angle gamma
	 */
    void SetData(const double a, const double b, const double c,
                 const double alpha, const double beta, const double gamma)
    {   _a = a; _b = b; _c = c;
        _alpha = alpha; _beta = beta; _gamma = gamma; }
    void SetData(const vector3 v1, const vector3 v2, const vector3 v3);
	//! set the offset to the origin to @p v1
    void SetOffset(const vector3 v1) { _offset = v1; }
    //! Set the space group symbol for this unit cell.
    //! Does not create an OBSymmetryData entry or attempt to convert
    //!  between different symbol notations
    void SetSpaceGroup(const std::string sg) { _spaceGroup = sg; }
    //! Set the Bravais lattice type for this unit cell
    void SetLatticeType(const LatticeType lt) { _lattice = lt; }

	//! \return Return vector a
    double GetA()    { return(_a);    }
	//! \return Return vector b
    double GetB()    { return(_b);    }
	//! \return Return vector c
    double GetC()    { return(_c);    }
	//! \return Return angle alpha
    double GetAlpha(){ return(_alpha);}
	//! \return Return angle beta
    double GetBeta() { return(_beta); }
	//! \return Return angle gamma
    double GetGamma(){ return(_gamma);}
    vector3 GetOffset() { return(_offset); }
    const std::string GetSpaceGroup() { return(_spaceGroup); }
    //! \return lattice type (based on angles and distances)
    LatticeType GetLatticeType();

    //! \return v1, v2, v3 cell vectors
    std::vector<vector3> GetCellVectors();
    //! \return v1, v2, v3 cell vectors as a 3x3 matrix
    matrix3x3	GetCellMatrix();
    //! \return The orthogonalization matrix, used for converting from fractional to Cartesian coords.
    matrix3x3 GetOrthoMatrix();
    //! \return The fractionalization matrix, used for converting from Cartesian to fractional coords.
    matrix3x3 GetFractionalMatrix();
    //! \return The cell volume (in Angstroms^3)
    double GetCellVolume();
};

//! \brief Used to hold data on conformers or geometry optimization steps
class OBAPI OBConformerData: public OBGenericData
{
protected:
  //! Dimensionalities of conformers
  std::vector<unsigned short>              _vDimension;
  //! Relative energies of conformers (preferably in kJ/mol)
  std::vector<double>                      _vEnergies;
  //! Atomic forces for each conformer
  std::vector< std::vector< vector3 > >    _vForces;
  //! Atomic velocities for each conformer (e.g., trajectories)
  std::vector< std::vector< vector3 > >    _vVelocity;
  //! Atomic displacements for each conformer (e.g., RMS distances)
  std::vector< std::vector< vector3 > >    _vDisplace;
  //! Additional data (as strings)
  std::vector<std::string>                 _vData;
    
public:
    OBConformerData();
    OBConformerData(const OBConformerData &);
		virtual OBGenericData* Clone(OBBase* parent) const{return new OBConformerData(*this);}
    ~OBConformerData()    {}

    OBConformerData &operator=(const OBConformerData &);

    void SetDimension(std::vector<unsigned short> vd) { _vDimension = vd; }
    void SetEnergies(std::vector<double> ve) { _vEnergies = ve; }
    void SetForces(std::vector< std::vector< vector3 > > vf) {_vForces = vf;}
    void SetVelocities(std::vector< std::vector< vector3 > > vv)
      { _vVelocity = vv; }
    void SetDisplacements(std::vector< std::vector< vector3 > > vd)
      { _vDisplace = vd; }
    void SetData(std::vector<std::string> vdat) { _vData = vdat; }

    std::vector<unsigned short> GetDimension() { return _vDimension; }
    std::vector<double>         GetEnergies()  { return _vEnergies; }
    std::vector< std::vector< vector3 > > GetForces() {return _vForces; }
    std::vector< std::vector< vector3 > > GetVelocities()
      {return _vVelocity;}
    std::vector< std::vector< vector3 > > GetDisplacements()
      {return _vDisplace;}
    std::vector<std::string>    GetData() { return _vData; }

};

//! \brief Used to hold the point-group and/or space-group symmetry
//! \todo Add support for translation between symbol notations.
//!        Add symmetry perception routines.
class OBAPI OBSymmetryData: public OBGenericData
{
protected:
    std::string _spaceGroup;
    std::string _pointGroup;
public:
    OBSymmetryData();
    OBSymmetryData(const OBSymmetryData &);
		virtual OBGenericData* Clone(OBBase* parent) const{return new OBSymmetryData(*this);}
    ~OBSymmetryData()    {}

    OBSymmetryData &operator=(const OBSymmetryData &);

    void SetData(std::string pg, std::string sg = "")
      { _pointGroup = pg; _spaceGroup = sg; }
    void SetPointGroup(std::string pg) { _pointGroup = pg; }
    void SetSpaceGroup(std::string sg) { _spaceGroup = sg; }

    std::string GetPointGroup() { return _pointGroup; }
    std::string GetSpaceGroup() { return _spaceGroup; }
};

//! \brief Used to hold the torsion data for a single rotatable bond
//! and all four atoms around it
class OBAPI OBTorsion
{
    friend class OBMol;
    friend class OBTorsionData;

protected:
    std::pair<OBAtom*,OBAtom*> _bc;
    //! double is angle in rads
    std::vector<triple<OBAtom*,OBAtom*,double> > _ads;

    OBTorsion()
    {
        _bc.first=0;
        _bc.second=0;
    }
    ;  //protected for use only by friend classes
    OBTorsion(OBAtom *, OBAtom *, OBAtom *, OBAtom *);

    std::vector<quad<OBAtom*,OBAtom*,OBAtom*,OBAtom*> > GetTorsions();

public:
    OBTorsion(const OBTorsion &);
    ~OBTorsion()
    {}

    OBTorsion& operator=(const OBTorsion &);

    void Clear();
    bool Empty()
    {
        return(_bc.first == 0 && _bc.second == 0);
    }

    bool AddTorsion(OBAtom *a,OBAtom *b, OBAtom *c,OBAtom *d);
    bool AddTorsion(quad<OBAtom*,OBAtom*,OBAtom*,OBAtom*> &atoms);

    bool SetAngle(double radians, unsigned int index = 0);
    bool SetData(OBBond *bond) { return false; }

    bool GetAngle(double &radians, unsigned int index =0);
    unsigned int GetBondIdx();
    unsigned int GetSize() const
    {
        return (unsigned int)_ads.size();
    }

    std::pair<OBAtom*,OBAtom*>                  GetBC()
    {
        return(_bc);
    }
    std::vector<triple<OBAtom*,OBAtom*,double> > GetADs()
    {
        return(_ads) ;
    }

    bool IsProtonRotor();
};

//! \brief Used to hold torsions as generic data for OBMol.
//! Filled by OBMol::FindTorsions()
class OBAPI OBTorsionData : public OBGenericData
{
    friend class OBMol;

protected:
    std::vector<OBTorsion> _torsions;

    OBTorsionData();
    OBTorsionData(const OBTorsionData &);

public:
    OBTorsionData &operator=(const OBTorsionData &);

		//TO BE MODIFIED TO rebase atom pointers
		virtual OBGenericData* Clone(OBBase* parent) const{return new OBTorsionData(*this);}

    void Clear();

    std::vector<OBTorsion> GetData() const
    {
        return _torsions;
    }
    unsigned int      GetSize() const
    {
        return (unsigned int)_torsions.size();
    }

    void SetData(OBTorsion &torsion);

    bool FillTorsionArray(std::vector<std::vector<unsigned int> > &torsions);
};

//! Used to hold the 3 atoms in an angle and the angle itself
class OBAPI OBAngle
{
    friend class OBMol;
    friend class OBAngleData;

protected:

    //member data

    OBAtom                *_vertex;
    std::pair<OBAtom*,OBAtom*>  _termini;
    double                  _radians;

    //protected member functions

    OBAngle();	//protect constructor for use only by friend classes
    OBAngle(OBAtom *vertex,OBAtom *a,OBAtom *b);

    triple<OBAtom*,OBAtom*,OBAtom*> GetAtoms();
    void SortByIndex();

public:

    OBAngle(const OBAngle &);
    ~OBAngle()
    {
        _vertex = NULL;
    }

    OBAngle &operator = (const OBAngle &);
    bool     operator ==(const OBAngle &);

    void  Clear();

    double GetAngle() const
    {
        return(_radians);
    }

    void  SetAngle(double radians)
    {
        _radians = radians;
    }
    void  SetAtoms(OBAtom *vertex,OBAtom *a,OBAtom *b);
    void  SetAtoms(triple<OBAtom*,OBAtom*,OBAtom*> &atoms);

};


//! \brief Used to hold all angles in a molecule as generic data for OBMol
class OBAPI OBAngleData : public OBGenericData
{
    friend class OBMol;

protected:
    std::vector<OBAngle> _angles;

    OBAngleData();
    OBAngleData(const OBAngleData &);
    std::vector<OBAngle> GetData() const
    {
        return(_angles);
    }

public:
    OBAngleData &operator =(const OBAngleData &);
		virtual OBGenericData* Clone(OBBase* parent) const{return new OBAngleData(*this);}

    void Clear();
    unsigned int FillAngleArray(int **angles, unsigned int &size);
    void         SetData(OBAngle &);
    unsigned int GetSize() const
    {
        return (unsigned int)_angles.size();
    }
};

enum atomreftype{output,input,calcvolume}; // sets which atom4ref is accessed

//! \brief Used to hold chiral inforamtion about the atom as OBGenericData
class OBAPI OBChiralData : public OBGenericData
{
    friend class OBMol;
    friend class OBAtom;

protected:
    std::vector<unsigned int> _atom4refs;
    int parity;
    std::vector<unsigned int> _atom4refo;
    std::vector<unsigned int> _atom4refc;

public:
    std::vector<unsigned int> GetAtom4Refs(atomreftype t) const;
    unsigned int GetAtomRef(int a,atomreftype t);

    OBChiralData();
    OBChiralData(const OBChiralData &src);
		virtual OBGenericData* Clone(OBBase* parent) const{return new OBChiralData(*this);}
    OBChiralData &operator =(const OBChiralData &);
    ~OBChiralData(){}

    void Clear();
    bool SetAtom4Refs(std::vector<unsigned int> atom4refs, atomreftype t);
    int AddAtomRef(unsigned int atomref, atomreftype t);
    unsigned int GetSize(atomreftype t) const;
};

//! Defines a map between serial numbers (e.g., in a PDB file) and OBAtom objects
class OBSerialNums : public OBGenericData
{
protected:
    std::map<int, OBAtom*> _serialMap;

public:

    OBSerialNums()
    {
        _attr = "obSerialNums";
        _type = OBGenericDataType::SerialNums;
    }
    OBSerialNums(const OBSerialNums &cp) : OBGenericData(cp)
    {
        _serialMap = cp._serialMap;
    }
		//Member varaiables contain OBAtom pointers, so copying only vlid within
		//same molecule, unless code modified, as in OBRotameterList 
		virtual OBGenericData* Clone(OBBase* parent) const{return new OBSerialNums(*this);}

    std::map<int,OBAtom*> &GetData()    { return _serialMap;    }
    void SetData(std::map<int,OBAtom*> &sm) { _serialMap = sm;  }

};

//****************doxygen for inline functions***********
/*!
**\fn OBTorsionData::GetSize() const
**\brief Gets the number of torsion structs
**\return integer count of the number of torsions
*/

/*!
**\fn OBTorsionData::GetData() const
**\brief Gets a vector of the OBTorsions
**\return the vector of torsions
*/

/*!
**\fn OBAngleData::GetSize() const
**\brief Gets the number of angles 
**\return integer count of the number of angles
*/

/*!
**\fn OBAngleData::GetData() const
**\brief Gets the angle vector data
**\return pointer to vector<OBAngle>
*/

/*!
**\fn OBAngle::SetAngle(double angle)
**\brief Sets the OBAngle angle value
**\param angle in radians
*/

/*!
**\fn OBAngle::GetAngle() const
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
**\return vector of A,D atom pointers and a double
*/

} //end namespace OpenBabel

#endif // OB_GENERIC_H

//! \file generic.h
//! \brief Handle generic data classes. Custom data for atoms, bonds, etc.

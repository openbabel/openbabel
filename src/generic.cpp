/**********************************************************************
generic.cpp - Handle OBGenericData classes.

Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2006 by Geoffrey R. Hutchison
Some portions Copyright (C) 2010 by David Lonie

This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/
#include <openbabel/babelconfig.h>

#include <string>
#include <set>

#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/ring.h>
#include <openbabel/obiter.h>
#include <openbabel/generic.h>
#include <openbabel/math/matrix3x3.h>
#include <openbabel/elements.h>

// needed for msvc to have at least one reference to AtomClass, AliasData in openbabel library
#include <openbabel/alias.h>

using namespace std;

namespace OpenBabel
{

  /** \class OBGenericData generic.h <openbabel/generic.h>

      OBGenericData is an abstract base class which defines an interface for
      storage, retrieval, and indexing of arbitrary generic data.
      Subclasses of OBGenericData can be used to store custom data
      on a per-atom, per-bond, per-molecule, or per-residue basis.
      Open Babel currently supports a small subset of chemical functionality
      as OBGenericData types, which will expand over time to support additional
      interconversion (e.g., spectroscopy, dynamics, surfaces...)

      For more information on currently supported types, please see
      the developer wiki:
      http://openbabel.org/wiki/Generic_Data

      For your own custom data, either define a custom subclass using
      an id from the OBGenericDataType::CustomData0 to
      OBGenericDataType::CustomData15 slots,
      or store your data as a string and use OBPairData for key/value access.
      The latter is <strong>highly</strong> recommended for various text
      descriptors
      e.g., in QSAR, atom or bond labels, or other textual data.

      <strong>New in Open Babel, version 2.1</strong>
      is the template-based OBPairTemplate,
      which can be used to store arbitrary data types. There are predefined
      types OBPairInteger and OBPairFloatingPoint for storing integers and
      floating-point values without converting to a string representation.

      Also <strong>new</strong> is the "source" or "origin" of a data
      entry, enumerated by DataOrigin. This can be accessed by
      SetOrigin() and GetOrigin(), as well as via "filtering" methods
      in OBBase, allowing you to separate data read in from a file,
      added by a user, or assigned by Open Babel internally.

      While the library and import routines will set DataOrigin correctly,
      you should try to annotate data added by your code. Typically this would
      either be userAdded or external. The former refers to something the
      user requested as an annotation, while the latter refers to annotations
      your code adds automatically.

      Example code using OBGenericData:

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
      for (k = vdata.begin();k != vdata.end();++k)
         if ((*k)->GetDataType() == OBGenericDataType::PairData)
         {
            ofs << ">  <" << (*k)->GetAttribute() << ">" << endl;
            ofs << ((OBPairData*)(*k))->GetValue() << endl << endl;
         }
      @endcode

      Similar code also works for OBGenericData stored in an OBAtom or
      OBBond or OBResidue. These examples show use of DataOrigin outside
      of the Open Babel library.

      @code
      string atomLabel; // e.g., from the user adding annotation to an atom
      if (!atom.HasData("UserLabel")) // stored textual data as an OBPairData
      {
         OBPairData *label = new OBPairData;
         label->SetAttribute("UserLabel");
         label->SetValue(atomLabel);
         label->SetOrigin(userInput); // set by user, not by Open Babel

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

      When designing a class derived from OBGenericData you must add a
      Clone() function. For classes used with OBMol this is used when
      an OBMol object is copied. If your class member variables contain
      pointers to atoms or bonds then it will be necessary to ensure
      that these are updated in Clone() to refer to the new molecule. Without
      these and similar pointers it is more likely that the very simple
      clone function
      @code
      virtual OBGenericData* Clone(OBBase* parent) const
         {return new MyNewClass(*this);}
      @endcode
      and the compiler generated copy constructor would be sufficient.

      It is recommended that, if possible, OBGenericData classes do not
      store atom and bond pointers. Using atom and bond indices instead
      would allow the simple version of Clone() above. See
      OBRotameterData::Clone for an example of a more complicated version.
      For classes which are not intended to support copying, Clone() can
      return NULL
      @code
      virtual OBGenericData* Clone(OBBase* parent) const
         {return NULL;}
      @endcode
      Clone() is a pure virtual function so that you need to decide what
      kind of function you need and include it explicitly.
  **/

  //
  //member functions for OBGenericData class
  //

  OBGenericData::OBGenericData(const std::string attr, const unsigned int type,
                               const DataOrigin  source):
    _attr(attr), _type(type), _source(source)
  { }

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

  OBCommentData::OBCommentData():
    OBGenericData("Comment", OBGenericDataType::CommentData)
  { }

  OBCommentData::OBCommentData(const OBCommentData &src) :
    OBGenericData(src), _data(src._data)
  {  }

  //
  //member functions for OBExternalBond class
  //
  OBExternalBond::OBExternalBond(OBAtom *atom,OBBond *bond,int idx):
    _idx(idx), _atom(atom), _bond(bond)
  {  }

  OBExternalBond::OBExternalBond(const OBExternalBond &src):
    _idx(src._idx), _atom(src._atom), _bond(src._bond)
  { }

  //
  //member functions for OBExternalBondData class
  //

  OBExternalBondData::OBExternalBondData():
    OBGenericData("ExternalBondData", OBGenericDataType::ExternalBondData,
                  perceived)
  { }

  void OBExternalBondData::SetData(OBAtom *atom,OBBond *bond,int idx)
  {
    OBExternalBond xb(atom,bond,idx);
    _vexbnd.push_back(xb);
  }

  //
  //member functions for OBPairData class
  //

  OBPairData::OBPairData() :
    OBGenericData("PairData", OBGenericDataType::PairData)
  { }

  //
  //member functions for OBVirtualBond class
  //

  OBVirtualBond::OBVirtualBond():
    OBGenericData("VirtualBondData", OBGenericDataType::VirtualBondData, perceived),
    _bgn(0), _end(0), _ord(0), _stereo(0)
  {  }

  OBVirtualBond::OBVirtualBond(unsigned int bgn, unsigned int end, unsigned int ord, int stereo):
    OBGenericData("VirtualBondData", OBGenericDataType::VirtualBondData, perceived),
    _bgn(bgn), _end(end), _ord(ord), _stereo(stereo)
  {  }

  //
  // member functions for OBUnitCell class
  //
  OBUnitCell::OBUnitCell():
    OBGenericData("UnitCell", OBGenericDataType::UnitCell),
    _mOrtho(matrix3x3()),
    _mOrient(matrix3x3()),
    _offset(vector3()),
    _spaceGroupName(""), _spaceGroup( NULL ),
    _lattice(Undefined)
  {  }

  OBUnitCell::OBUnitCell(const OBUnitCell &src) :
    OBGenericData("UnitCell", OBGenericDataType::UnitCell),
    _mOrtho(src._mOrtho),
    _mOrient(src._mOrient),
    _offset(src._offset),
    _spaceGroupName(src._spaceGroupName), _spaceGroup(src._spaceGroup),
    _lattice(src._lattice)
  {  }

  OBUnitCell & OBUnitCell::operator=(const OBUnitCell &src)
  {
    if(this == &src)
      return(*this);

    _mOrtho = src._mOrtho;
    _mOrient = src._mOrient;
    _offset = src._offset;

    _spaceGroup = src._spaceGroup;
    _spaceGroupName = src._spaceGroupName;
    _lattice = src._lattice;

    return(*this);
  }

  //! Implements <a href="http://qsar.sourceforge.net/dicts/blue-obelisk/index.xhtml#calculateOrthogonalisationMatrix">blue-obelisk:calculateOrthogonalisationMatrix</a>
  void OBUnitCell::SetData(const double a, const double b, const double c,
                           const double alpha, const double beta, const double gamma)
  {
    _mOrtho.FillOrth(alpha, beta, gamma, a, b, c);
    _mOrient = matrix3x3(1);
    _spaceGroup = NULL;
    _spaceGroupName = "";
    _lattice = OBUnitCell::Undefined;
  }

  //! Implements <a href="http://qsar.sourceforge.net/dicts/blue-obelisk/index.xhtml#calculateOrthogonalisationMatrix">blue-obelisk:calculateOrthogonalisationMatrix</a>
  void OBUnitCell::SetData(const vector3 v1, const vector3 v2, const vector3 v3)
  {
    matrix3x3 m (v1, v2, v3);
    _mOrtho.FillOrth(vectorAngle(v2,v3), // alpha
                     vectorAngle(v1,v3), // beta
                     vectorAngle(v1,v2), // gamma
                     v1.length(),        // a
                     v2.length(),        // b
                     v3.length());       // c
    _mOrient = m.transpose() * _mOrtho.inverse();
    _spaceGroup = NULL;
    _spaceGroupName = "";
    _lattice = OBUnitCell::Undefined;
  }

  //! Implements <a href="http://qsar.sourceforge.net/dicts/blue-obelisk/index.xhtml#calculateOrthogonalisationMatrix">blue-obelisk:calculateOrthogonalisationMatrix</a>
  void OBUnitCell::SetData(const matrix3x3 m)
  {
    SetData(m.GetRow(0), m.GetRow(1), m.GetRow(2));
  }

  void OBUnitCell::SetOffset(const vector3 v1)
  {
    _offset = v1;
  }

  //! Implements <a href="http://qsar.sourceforge.net/dicts/blue-obelisk/index.xhtml#convertNotionalIntoCartesianCoordinates">blue-obelisk:convertNotionalIntoCartesianCoordinates</a>
  vector<vector3> OBUnitCell::GetCellVectors() const
  {
    vector<vector3> v;
    v.reserve(3);

    matrix3x3 m = GetCellMatrix();

    v.push_back(m.GetRow(0));
    v.push_back(m.GetRow(1));
    v.push_back(m.GetRow(2));

    return v;
  }

  matrix3x3 OBUnitCell::GetCellMatrix() const
  {
    return (_mOrient * _mOrtho).transpose();
  }

  matrix3x3 OBUnitCell::GetOrthoMatrix() const
  {
    return _mOrtho;
  }

  matrix3x3 OBUnitCell::GetOrientationMatrix() const
  {
    return _mOrient;
  }

  matrix3x3 OBUnitCell::GetFractionalMatrix() const
  {
    return _mOrtho.inverse();
  }

  vector3 OBUnitCell::FractionalToCartesian(vector3 frac) const
  {
    return _mOrient * _mOrtho * frac + _offset;
  }

  vector3 OBUnitCell::CartesianToFractional(vector3 cart) const
  {
    return _mOrtho.inverse() * _mOrient.inverse() * (cart - _offset);
  }

  vector3 OBUnitCell::WrapCartesianCoordinate(vector3 cart) const
  {
    vector3 v = CartesianToFractional(cart);
    v = WrapFractionalCoordinate(v);
    return FractionalToCartesian(v);
  }

  vector3 OBUnitCell::WrapFractionalCoordinate(vector3 frac) const
  {
    double x = fmod(frac.x(), 1);
    double y = fmod(frac.y(), 1);
    double z = fmod(frac.z(), 1);
    if (x < 0) x += 1;
    if (y < 0) y += 1;
    if (z < 0) z += 1;

#define LIMIT 0.999999
    if (x > LIMIT)
      x -= 1;
    if (y > LIMIT)
      y -= 1;
    if (z > LIMIT)
      z -= 1;
#undef LIMIT

    // Fuzzy logic from Francois-Xavier
#define EPSILON 1.0e-6
    if (x > 1 - EPSILON || x < EPSILON)
      x = 0.0;
    if (y > 1 - EPSILON || y < EPSILON)
      y = 0.0;
    if (z > 1 - EPSILON || z < EPSILON)
      z = 0.0;
#undef EPSILON

    return vector3(x, y, z);
  }

  OBUnitCell::LatticeType OBUnitCell::GetLatticeType( int spacegroup ) const
  {
    //	1-2 	Triclinic
    //	3-15 	Monoclinic
    //	16-74	Orthorhombic
    //	75-142 	Tetragonal
    //	143-167 Rhombohedral
    //	168-194 Hexagonal
    //	195-230 Cubic

    if ( spacegroup == 0  && _spaceGroup)
      spacegroup = _spaceGroup->GetId();

    if ( spacegroup <= 0 )
      return OBUnitCell::Undefined;

    else if ( spacegroup == 1 ||
              spacegroup == 2 )
      return OBUnitCell::Triclinic;

    else if ( spacegroup >= 3 &&
              spacegroup <= 15 )
      return OBUnitCell::Monoclinic;

    else if ( spacegroup >= 16 &&
              spacegroup <= 74 )
      return OBUnitCell::Orthorhombic;

    else if ( spacegroup >= 75 &&
              spacegroup <= 142 )
      return OBUnitCell::Tetragonal;

    else if ( spacegroup >= 143 &&
              spacegroup <= 167 )
      return OBUnitCell::Rhombohedral;

    else if ( spacegroup >= 168 &&
              spacegroup <= 194 )
      return OBUnitCell::Hexagonal;

    else if ( spacegroup >= 195 &&
              spacegroup <= 230 )
      return OBUnitCell::Cubic;

    //just to be extra sure
    else // ( spacegroup > 230 )
      return OBUnitCell::Undefined;
  }

  OBUnitCell::LatticeType OBUnitCell::GetLatticeType() const
  {
    if (_lattice != Undefined)
      return _lattice;
    else if (_spaceGroup != NULL)
      return GetLatticeType(_spaceGroup->GetId());

    double a = GetA();
    double b = GetB();
    double c = GetC();
    double alpha = GetAlpha();
    double beta  = GetBeta();
    double gamma = GetGamma();

    unsigned int rightAngles = 0;
    if (IsApprox(alpha, 90.0, 1.0e-3)) rightAngles++;
    if (IsApprox(beta,  90.0, 1.0e-3)) rightAngles++;
    if (IsApprox(gamma, 90.0, 1.0e-3)) rightAngles++;

    // recast cache member "_lattice" as mutable
    OBUnitCell::LatticeType *lattice =
      const_cast<OBUnitCell::LatticeType*>(&_lattice);

    switch (rightAngles)
      {
      case 3:
        if (IsApprox(a, b, 1.0e-4) && IsApprox(b, c, 1.0e-4))
          *lattice = Cubic;
        else if (IsApprox(a, b, 1.0e-4) || IsApprox(b, c, 1.0e-4))
          *lattice = Tetragonal;
        else
          *lattice = Orthorhombic;
        break;
      case 2:
        if ( (IsApprox(alpha, 120.0, 1.0e-3)
              || IsApprox(beta, 120.0, 1.0e-3)
              || IsApprox(gamma, 120.0f, 1.0e-3))
             && (IsApprox(a, b, 1.0e-4) || IsApprox(b, c, 1.0e-4)) )
          *lattice = Hexagonal;
        else
          *lattice = Monoclinic;
        break;
      default:
        if (IsApprox(a, b, 1.0e-4) && IsApprox(b, c, 1.0e-4))
          *lattice = Rhombohedral;
        else
          *lattice = Triclinic;
      }

    return *lattice;
  }

  int OBUnitCell::GetSpaceGroupNumber( std::string name) const
  {
    static const char * const spacegroups[] = {
      "P1", "P-1", "P2", "P2(1)", "C2", "Pm", "Pc", "Cm", "Cc", "P2/m",
      "P2(1)/m", "C2/m", "P2/c", "P2(1)/c", "C2/c", "P222", "P222(1)",
      "P2(1)2(1)2", "P2(1)2(1)2(1)", "C222(1)", "C222", "F222", "I222",
      "I2(1)2(1)2(1)", "Pmm2", "Pmc2(1)", "Pcc2", "Pma2", "Pca2(1)", "Pnc2",
      "Pmn2(1)", "Pba2", "Pna2(1)", "Pnn2", "Cmm2", "Cmc2(1)", "Ccc2", "Amm2",
      "Abm2", "Ama2", "Aba2", "Fmm2", "Fdd2", "Imm2", "Iba2", "Ima2", "Pmmm",
      "Pnnn", "Pccm", "Pban", "Pmma", "Pnna", "Pmna", "Pcca", "Pbam", "Pccn",
      "Pbcm", "Pnnm", "Pmmn", "Pbcn", "Pbca", "Pnma", "Cmcm", "Cmca", "Cmmm",
      "Cccm", "Cmma", "Ccca", "Fmmm", "Fddd", "Immm", "Ibam", "Ibca", "Imma",
      "P4", "P4(1)", "P4(2)", "P4(3)", "I4", "I4(1)", "P-4", "I-4", "P4/m",
      "P4(2)/m", "P4/n", "P4(2)/n", "I4/m", "I4(1)/a", "P422", "P42(1)2",
      "P4(1)22", "P4(1)2(1)2", "P4(2)22", "P4(2)2(1)2", "P4(3)22", "P4(3)2(1)2",
      "I422", "I4(1)22", "P4mm", "P4bm", "P4(2)cm", "P4(2)nm", "P4cc", "P4nc",
      "P4(2)mc", "P4(2)bc", "I4mm", "I4cm", "I4(1)md", "I4(1)cd", "P-42m",
      "P-42c", "P-42(1)m", "P-42(1)c", "P-4m2", "P-4c2", "P-4b2", "P-4n2",
      "I-4m2", "I-4c2", "I-42m", "I-42d", "P4/mmm", "P4/mcc", "P4/nbm",
      "P4/nnc", "P4/mbm", "P4/mnc", "P4/nmm", "P4/ncc", "P4(2)/mmc",
      "P4(2)/mcm", "P4(2)/nbc", "P4(2)/nnm", "P4(2)/mbc", "P4(2)/mnm",
      "P4(2)/nmc", "P4(2)/ncm", "I4/mmm", "I4/mcm", "I4(1)/amd", "I4(1)/acd",
      "P3", "P3(1)", "P3(2)", "R3", "P-3", "R-3", "P312", "P321", "P3(1)12",
      "P3(1)21", "P3(2)12", "P3(2)21", "R32", "P3m1", "P31m", "P3c1", "P31c",
      "R3m", "R3c", "P-31m", "P-31c", "P-3m1", "P-3c1", "R-3m", "R-3c", "P6",
      "P6(1)", "P6(5)", "P6(2)", "P6(4)", "P6(3)", "P-6", "P6/m", "P6(3)/m",
      "P622", "P6(1)22", "P6(5)22", "P6(2)22", "P6(4)22", "P6(3)22", "P6mm",
      "P6cc", "P6(3)cm", "P6(3)mc", "P-6m2", "P-6c2", "P-62m", "P-62c",
      "P6/mmm", "P6/mcc", "P6(3)/mcm", "P6(3)/mmc", "P23", "F23", "I23",
      "P2(1)3", "I2(1)3", "Pm-3", "Pn-3", "Fm-3", "Fd-3", "Im-3", "Pa-3",
      "Ia-3", "P432", "P4(2)32", "F432", "F4(1)32", "I432", "P4(3)32",
      "P4(1)32", "I4(1)32", "P-43m", "F4-3m", "I-43m", "P-43n", "F-43c",
      "I-43d", "Pm-3m", "Pn-3n", "Pm-3n", "Pn-3m", "Fm-3m", "Fm-3c",
      "Fd-3m", "Fd-3c", "Im-3m", "Ia-3d"
    };

    if (name.length () == 0)
      {
        if (_spaceGroup != NULL)
          return _spaceGroup->GetId();
        else
          name = _spaceGroupName;
      }
    static const int numStrings = sizeof( spacegroups ) / sizeof( spacegroups[0] );
    for ( int i = 0; i < numStrings; ++i ) {
      if (name == spacegroups[i] ) {
        return i+1;
      }
    }
    return 0; //presumably never reached
  }

  // Whether two points (given in fractional coordinates) are close enough
  // to be considered duplicates.
  bool areDuplicateAtoms (vector3 v1, vector3 v2)
  {
    vector3 dr = v2 - v1;
    if (dr.x() < -0.5)
      dr.SetX(dr.x() + 1);
    if (dr.x() > 0.5)
      dr.SetX(dr.x() - 1);
    if (dr.y() < -0.5)
      dr.SetY(dr.y() + 1);
    if (dr.y() > 0.5)
      dr.SetY(dr.y() - 1);
    if (dr.z() < -0.5)
      dr.SetZ(dr.z() + 1);
    if (dr.z() > 0.5)
      dr.SetZ(dr.z() - 1);

    return (dr.length_2() < 1e-6);
  }

  void OBUnitCell::FillUnitCell(OBMol *mol)
  {
    const SpaceGroup *sg = GetSpaceGroup(); // the actual space group and transformations for this unit cell

    if(sg == NULL)
      return ;

    // For each atom, we loop through: convert the coords back to inverse space, apply the transformations and create new atoms
    vector3 baseV, uniqueV, updatedCoordinate;
    list<vector3> transformedVectors; // list of symmetry-defined copies of the atom
    list<vector3>::iterator transformIter;
    list<OBAtom*>::iterator deleteIter, atomIter;
    OBAtom *newAtom;
    list<OBAtom*> atoms, atomsToDelete;
    char hash[22];
    set<string> coordinateSet;

    // Check original mol for duplicates
    FOR_ATOMS_OF_MOL(atom, *mol) {
      baseV = atom->GetVector();
      baseV = CartesianToFractional(baseV);
      baseV = WrapFractionalCoordinate(baseV);
      snprintf(hash, 22, "%03d,%.3f,%.3f,%.3f", atom->GetAtomicNum(), baseV.x(), baseV.y(), baseV.z());
      if (coordinateSet.insert(hash).second) { // True if new entry
        atoms.push_back(&(*atom));
      } else {
        atomsToDelete.push_back(&(*atom));
      }
    }
    for (deleteIter = atomsToDelete.begin(); deleteIter != atomsToDelete.end(); ++deleteIter) {
      mol->DeleteAtom(*deleteIter);
    }

    // Cross-check all transformations for duplicity
    for (atomIter = atoms.begin(); atomIter != atoms.end(); ++atomIter) {
      uniqueV = (*atomIter)->GetVector();
      uniqueV = CartesianToFractional(uniqueV);
      uniqueV = WrapFractionalCoordinate(uniqueV);

      transformedVectors = sg->Transform(uniqueV);
      for (transformIter = transformedVectors.begin();
        transformIter != transformedVectors.end(); ++transformIter) {
        updatedCoordinate = WrapFractionalCoordinate(*transformIter);

        // Check if the transformed coordinate is a duplicate of an atom
        snprintf(hash, 22, "%03d,%.3f,%.3f,%.3f", (*atomIter)->GetAtomicNum(), updatedCoordinate.x(),
                 updatedCoordinate.y(), updatedCoordinate.z());
        if (coordinateSet.insert(hash).second) {
          newAtom = mol->NewAtom();
          newAtom->Duplicate(*atomIter);
          newAtom->SetVector(FractionalToCartesian(updatedCoordinate));
        }
      } // end loop of transformed atoms
    } // end loop of atoms
    SetSpaceGroup(1); // We've now applied the symmetry, so we should act like a P1 unit cell
  }

  /// @todo Remove nonconst overloads in OBUnitCell on next version bump.
#define OBUNITCELL_CALL_CONST_OVERLOAD(_type, _name) \
  _type OBUnitCell::_name() \
  { \
    return const_cast<const OBUnitCell*>(this)->_name(); \
  }
#define OBUNITCELL_CALL_CONST_OVERLOAD_ARG(_type, _name, _argsig) \
  _type OBUnitCell::_name( _argsig arg1 ) \
  { \
    return const_cast<const OBUnitCell*>(this)->_name(arg1); \
  }
  OBUNITCELL_CALL_CONST_OVERLOAD(double, GetA);
  OBUNITCELL_CALL_CONST_OVERLOAD(double, GetB);
  OBUNITCELL_CALL_CONST_OVERLOAD(double, GetC);
  OBUNITCELL_CALL_CONST_OVERLOAD(double, GetAlpha);
  OBUNITCELL_CALL_CONST_OVERLOAD(double, GetBeta);
  OBUNITCELL_CALL_CONST_OVERLOAD(double, GetGamma);
  OBUNITCELL_CALL_CONST_OVERLOAD(vector3, GetOffset);
  OBUNITCELL_CALL_CONST_OVERLOAD_ARG(OBUnitCell::LatticeType,
                                     GetLatticeType, int);
  OBUNITCELL_CALL_CONST_OVERLOAD(OBUnitCell::LatticeType, GetLatticeType);
  OBUNITCELL_CALL_CONST_OVERLOAD(std::vector<vector3>, GetCellVectors);
  OBUNITCELL_CALL_CONST_OVERLOAD(matrix3x3, GetCellMatrix );
  OBUNITCELL_CALL_CONST_OVERLOAD(matrix3x3, GetOrthoMatrix );
  OBUNITCELL_CALL_CONST_OVERLOAD(matrix3x3, GetOrientationMatrix );
  OBUNITCELL_CALL_CONST_OVERLOAD(matrix3x3, GetFractionalMatrix );
  OBUNITCELL_CALL_CONST_OVERLOAD_ARG(vector3, FractionalToCartesian, vector3);
  OBUNITCELL_CALL_CONST_OVERLOAD_ARG(vector3, CartesianToFractional, vector3);
  OBUNITCELL_CALL_CONST_OVERLOAD_ARG(vector3, WrapCartesianCoordinate, vector3);
  OBUNITCELL_CALL_CONST_OVERLOAD_ARG(vector3, WrapFractionalCoordinate, vector3);
  OBUNITCELL_CALL_CONST_OVERLOAD_ARG(int, GetSpaceGroupNumber, std::string);
  OBUNITCELL_CALL_CONST_OVERLOAD(double, GetCellVolume);

  double OBUnitCell::GetA() const
  {
    return _mOrtho.GetColumn(0).length();
  }

  double OBUnitCell::GetB() const
  {
    return _mOrtho.GetColumn(1).length();
  }

  double OBUnitCell::GetC() const
  {
    return _mOrtho.GetColumn(2).length();
  }

  double OBUnitCell::GetAlpha() const
  {
    return vectorAngle(_mOrtho.GetColumn(1), _mOrtho.GetColumn(2));
  }

  double OBUnitCell::GetBeta() const
  {
    return vectorAngle(_mOrtho.GetColumn(0), _mOrtho.GetColumn(2));
  }

  double OBUnitCell::GetGamma() const
  {
    return vectorAngle(_mOrtho.GetColumn(0), _mOrtho.GetColumn(1));
  }

  vector3 OBUnitCell::GetOffset() const
  {
    return _offset;
  }

  double OBUnitCell::GetCellVolume() const
  {
    return fabs(GetCellMatrix().determinant());
  }

  //
  // member functions for OBSymmetryData class
  //
  OBSymmetryData::OBSymmetryData():
    OBGenericData("Symmetry", OBGenericDataType::SymmetryData)
  { }

  OBSymmetryData::OBSymmetryData(const OBSymmetryData &src) :
    OBGenericData(src._attr, src._type, src._source),
    _spaceGroup(src._spaceGroup),
    _pointGroup(src._pointGroup)
  {  }

  OBSymmetryData & OBSymmetryData::operator=(const OBSymmetryData &src)
  {
    if(this == &src)
      return(*this);

    _pointGroup = src._pointGroup;
    _spaceGroup = src._spaceGroup;
    _source = src._source;

    return(*this);
  }

  OBConformerData::OBConformerData() :
    OBGenericData("Conformers", OBGenericDataType::ConformerData)
  {  }

  OBConformerData::OBConformerData(const OBConformerData &src) :
    OBGenericData("Conformers", OBGenericDataType::ConformerData),
    _vDimension(src._vDimension),
    _vEnergies(src._vEnergies), _vForces(src._vForces),
    _vVelocity(src._vVelocity), _vDisplace(src._vDisplace),
    _vData(src._vData)
  {  }

  OBConformerData & OBConformerData::operator=(const OBConformerData &src)
  {
    if(this == &src)
      return(*this);

    _source = src._source;

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

  OBRingData::OBRingData() :
    OBGenericData("RingData", OBGenericDataType::RingData)
  {
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

    for(ring = _vr.begin();ring != _vr.end();++ring)
      {
        OBRing *newring = new OBRing;
        (*newring) = (**ring);	//copy data to new object
        (*ring)    = newring;	//repoint new pointer to new copy of data
      }
  }

  OBRingData::~OBRingData()
  {
    vector<OBRing*>::iterator ring;
    for (ring = _vr.begin();ring != _vr.end();++ring)
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
    for(ring = _vr.begin();ring != _vr.end();++ring)
      {
        delete &*ring;	//deallocate old rings to prevent memory leak
      }

    _vr.clear();
    _vr = src._vr;	//copy vector properties

    for(ring = _vr.begin();ring != _vr.end();++ring)
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

  OBRing *OBRingData::BeginRing(std::vector<OBRing*>::iterator &i)
  {
    i = _vr.begin();
    return((i == _vr.end()) ? (OBRing*)NULL : (OBRing*)*i);
  }

  OBRing *OBRingData::NextRing(std::vector<OBRing*>::iterator &i)
  {
    ++i;
    return((i == _vr.end()) ? (OBRing*)NULL : (OBRing*)*i);
  }

  //
  //member functions for OBAngle class - stores all angles
  //

  /*!
  **\brief Angle default constructor
  */
  OBAngle::OBAngle():
    _vertex((OBAtom *)NULL), _termini((OBAtom *)NULL, (OBAtom *)NULL), _radians(0.0)
  {  }

  /*!
  **\brief Angle constructor
  */
  OBAngle::OBAngle(OBAtom *vertex,OBAtom *a,OBAtom *b):
    _vertex(vertex), _termini(a, b)
  {
    SortByIndex();
  }

  /*!
  **\brief OBAngle copy constructor
  */
  OBAngle::OBAngle(const OBAngle &src):
    _vertex(src._vertex), _termini(src._termini), _radians(src._radians)
  {  }

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
    :	OBGenericData("AngleData", OBGenericDataType::AngleData)
  {  }

  /*!
  **\brief OBAngleData copy constructor
  */
  OBAngleData::OBAngleData(const OBAngleData &src)
    :	OBGenericData(src), _angles(src._angles)
  {  }

  /*!
  **\brief OBAngleData assignment operator
  */
  OBAngleData& OBAngleData::operator =(const OBAngleData &src)
  {
    if (this == &src)
      return(*this);

    _source = src._source;
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
  **\return True if successful
  */
  bool OBAngleData::FillAngleArray(std::vector<std::vector<unsigned int> > &angles)
  {
    if(_angles.empty())
      return(false);

    vector<OBAngle>::iterator angle;

    angles.clear();
    angles.resize(_angles.size());

    unsigned int ct = 0;

    for( angle=_angles.begin(); angle!=_angles.end(); ++angle,ct++)
      {
        angles[ct].resize(3);
        angles[ct][0] = angle->_vertex->GetIdx() - 1;
        angles[ct][1] = angle->_termini.first->GetIdx() - 1;
        angles[ct][2] = angle->_termini.second->GetIdx() - 1;
      }

    return(true);
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
    for( angle=_angles.begin(); angle!=_angles.end(); ++angle)
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

    for(ad = _ads.begin();ad != _ads.end();++ad)
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
    for(ad = _ads.begin();ad != _ads.end() && (Aprotor || Dprotor);++ad)
      {
        if (ad->first->GetAtomicNum() != OBElements::Hydrogen)
          Aprotor = false;
        if (ad->second->GetAtomicNum() != OBElements::Hydrogen)
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
    : OBGenericData("TorsionData", OBGenericDataType::TorsionData)
  {  }

  //
  //member functions for OBTorsionData class - stores OBTorsion set
  //
  OBTorsionData::OBTorsionData(const OBTorsionData &src)
    :	OBGenericData(src), _torsions(src._torsions)
  {  }

  OBTorsionData& OBTorsionData::operator =(const OBTorsionData &src)
  {
    if (this == &src)
      return(*this);

    OBGenericData::operator =(src);

    _source = src._source;
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
    if(_torsions.empty())
      return(false);

    vector<quad<OBAtom*,OBAtom*,OBAtom*,OBAtom*> > tmpquads,quads;
    vector<quad<OBAtom*,OBAtom*,OBAtom*,OBAtom*> >::iterator thisQuad;
    vector<OBTorsion>::iterator torsion;

    //generate set of all 4 atom abcd's from torsion structure
    for (torsion = _torsions.begin();torsion != _torsions.end();++torsion)
      {
        tmpquads = torsion->GetTorsions();
        for(thisQuad = tmpquads.begin();thisQuad != tmpquads.end();++thisQuad)
          quads.push_back(*thisQuad);
      }

    //fill array of torsion atoms

    torsions.clear();
    torsions.resize(quads.size());

    unsigned int ct = 0;

    for (thisQuad = quads.begin();thisQuad != quads.end();++thisQuad,++ct)
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
//member functions for OBDOSData class
//

/*!
**\brief Assign the data
**\param fermi The Fermi energy in eV
**\param vEnergies Energy levels in eV
**\param vDensities Density of states in (number of states) / (unit cell)
**\param vIntegration Integrated DOS
*/
void OBDOSData::SetData(double fermi,
                        const std::vector<double> & vEnergies,
                        const std::vector<double> & vDensities,
                        const std::vector<double> & vIntegration)
{
  this->_fermi = fermi;
  this->_vEnergies = vEnergies;
  this->_vIntegration = vIntegration;
  this->_vDensities = vDensities;
}

  // member functions for OBOrbitalData

  void OBOrbitalData::LoadClosedShellOrbitals(std::vector<double> energies, std::vector<std::string> symmetries, unsigned int alphaHOMO)
  {
    if (energies.size() < symmetries.size())
      return; // something is very weird -- it's OK to pass no symmetries (we'll assume "A")
    if (energies.size() == 0)
      return;
    if (alphaHOMO > energies.size())
      return;

    _alphaHOMO = alphaHOMO;
    _alphaOrbitals.clear();
    _betaHOMO = 0;
    _betaOrbitals.clear();
    _openShell = false;

    if (symmetries.size() < energies.size()) // pad with "A" symmetry
      for (unsigned int i = symmetries.size(); i < energies.size(); ++i)
        symmetries.push_back("A");

    OBOrbital currentOrbital;
    for (unsigned int i = 0; i < energies.size(); ++i)
      {
        if (i < alphaHOMO)
          currentOrbital.SetData(energies[i], 2.0, symmetries[i]);
        else
          currentOrbital.SetData(energies[i], 0.0, symmetries[i]);

        _alphaOrbitals.push_back(currentOrbital);
      }
  }

  void OBOrbitalData::LoadAlphaOrbitals(std::vector<double> energies, std::vector<std::string> symmetries, unsigned int alphaHOMO)
  {
    if (energies.size() < symmetries.size())
      return; // something is very weird -- it's OK to pass no symmetries (we'll assume "A")
    if (energies.size() == 0)
      return;
    if (alphaHOMO > energies.size())
      return;

    _alphaHOMO = alphaHOMO;
    _alphaOrbitals.clear();
    _openShell = true;

    if (symmetries.size() < energies.size()) // pad with "A" symmetry
      for (unsigned int i = symmetries.size(); i < energies.size(); ++i)
        symmetries.push_back("A");

    OBOrbital currentOrbital;
    for (unsigned int i = 0; i < energies.size(); ++i)
      {
        if (i < alphaHOMO)
          currentOrbital.SetData(energies[i], 2.0, symmetries[i]);
        else
          currentOrbital.SetData(energies[i], 0.0, symmetries[i]);

        _alphaOrbitals.push_back(currentOrbital);
      }
  }

  void OBOrbitalData::LoadBetaOrbitals(std::vector<double> energies, std::vector<std::string> symmetries, unsigned int betaHOMO)
  {
    if (energies.size() < symmetries.size())
      return; // something is very weird -- it's OK to pass no symmetries (we'll assume "A")
    if (energies.size() == 0)
      return;
    if (betaHOMO > energies.size())
      return;

    _betaHOMO = betaHOMO;
    _betaOrbitals.clear();
    _openShell = true;

    if (symmetries.size() < energies.size()) // pad with "A" symmetry
      for (unsigned int i = symmetries.size(); i < energies.size(); ++i)
        symmetries.push_back("A");

    OBOrbital currentOrbital;
    for (unsigned int i = 0; i < energies.size(); ++i)
      {
        if (i < betaHOMO)
          currentOrbital.SetData(energies[i], 2.0, symmetries[i]);
        else
          currentOrbital.SetData(energies[i], 0.0, symmetries[i]);

        _betaOrbitals.push_back(currentOrbital);
      }
  }

//
//member functions for OBElectronicTransitionData class
//

/*!
**\brief Assign the basic excitation data
**\param vWavelengths Wavelengths in nm
**\param vForces Oscillator strengths
*/
void OBElectronicTransitionData::SetData(const std::vector<double> & vWavelengths,
                                  const std::vector<double> & vForces)
{
  this->_vWavelengths = vWavelengths;
  this->_vForces = vForces;
}

/*!
**\brief Assign the electronic dipole strengths
**\param vEDipole Electronic dipole moment strength
*/
void OBElectronicTransitionData::SetEDipole(const std::vector<double> & vEDipole)
{
  this->_vEDipole = vEDipole;
}

/*!
**\brief Assign the rotatory strengths (velocity)
**\param vRotatoryStrengthsVelocity Vector containing the rotatory strengths
*/
void OBElectronicTransitionData::SetRotatoryStrengthsVelocity(const std::vector<double> & vRotatoryStrengthsVelocity)
{
  this->_vRotatoryStrengthsVelocity = vRotatoryStrengthsVelocity;
}

/*!
**\brief Assign the rotatory strengths (length)
**\param vRotatoryStrengthsLength Vector containing the rotatory strengths
*/
void OBElectronicTransitionData::SetRotatoryStrengthsLength(const std::vector<double> & vRotatoryStrengthsLength)
{
  this->_vRotatoryStrengthsLength = vRotatoryStrengthsLength;
}

//
//member functions for OBVibrationData class
//

/*!
**\brief Assign the data
**\param vLx Normal modes in 1/sqrt(a.u.)
**\param vFrequencies Harmonic frequencies in inverse centimeters
**\param vIntensities Infrared absorption intensities in KM/Mole
*/
void OBVibrationData::SetData(const std::vector< std::vector< vector3 > > & vLx,
                              const std::vector<double> & vFrequencies,
                              const std::vector<double> & vIntensities)
{
  this->_vLx = vLx;
  this->_vFrequencies = vFrequencies;
  this->_vIntensities = vIntensities;
}

/*!
**\brief Assign the data
**\param vLx Normal modes in 1/sqrt(a.u.)
**\param vFrequencies Harmonic frequencies in inverse centimeters
**\param vIntensities Infrared absorption intensities in KM/Mole
**\param vRamanActivities Raman activities
*/
void OBVibrationData::SetData(const std::vector< std::vector< vector3 > > & vLx,
                              const std::vector<double> & vFrequencies,
                              const std::vector<double> & vIntensities,
                              const std::vector<double> & vRamanActivities)
{
  this->_vLx = vLx;
  this->_vFrequencies = vFrequencies;
  this->_vIntensities = vIntensities;
  this->_vRamanActivities = vRamanActivities;
}


/*!
**\brief Get the number of frequencies
*/
unsigned int OBVibrationData::GetNumberOfFrequencies() const
{
  return this->_vFrequencies.size();
}

void OBFreeGrid::Clear()
{
  _points.clear();
}

} //end namespace OpenBabel

//! \file generic.cpp
//! \brief Handle OBGenericData classes. Custom data for atoms, bonds, etc.

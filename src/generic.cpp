/**********************************************************************
generic.cpp - Handle OBGenericData classes.
 
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2006 by Geoffrey R. Hutchison
 
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
#include <openbabel/babelconfig.h>

#include <string>

#include <openbabel/mol.h>
#include <openbabel/generic.h>
#include <openbabel/math/matrix3x3.h>

// needed for msvc to have at least one reference to AtomClass, AliasData in openbabel library
#include <openbabel/atomclass.h>
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
      http://openbabel.sourceforge.net/wiki/Generic_Data

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

  OBCommentData::OBCommentData() :
      OBGenericData("Comment", OBGenericDataType::CommentData)
  {  }

  OBCommentData::OBCommentData(const OBCommentData &src) :
      OBGenericData(src), m_data(src.m_data)
  {  }

  //
  //member functions for OBExternalBond class
  //
  OBExternalBond::OBExternalBond(OBAtom *atom, OBBond *bond, int idx):
      m_idx(idx), m_atom(atom), m_bond(bond)
  {  }

  OBExternalBond::OBExternalBond(const OBExternalBond &src):
      m_idx(src.m_idx), m_atom(src.m_atom), m_bond(src.m_bond)
  {  }

  //
  //member functions for OBExternalBondData class
  //

  OBExternalBondData::OBExternalBondData():
      OBGenericData("ExternalBondData", OBGenericDataType::ExternalBondData, perceived)
  {  }

  void OBExternalBondData::SetData(OBAtom *atom, OBBond *bond, int idx)
  {
    OBExternalBond xb(atom,bond,idx);
    m_vexbnd.push_back(xb);
  }

  //
  //member functions for OBPairData class
  //

  OBPairData::OBPairData() :
      OBGenericData("PairData", OBGenericDataType::PairData)
  {  }

  //
  //member functions for OBVirtualBond class
  //

  OBVirtualBond::OBVirtualBond():
    OBGenericData("VirtualBondData", OBGenericDataType::VirtualBondData, perceived),
    m_bgn(0), m_end(0), m_ord(0), m_stereo(0)
  {  }

  OBVirtualBond::OBVirtualBond(int bgn,int end,int ord,int stereo):
    OBGenericData("VirtualBondData", OBGenericDataType::VirtualBondData, perceived),
    m_bgn(bgn), m_end(end), m_ord(ord), m_stereo(stereo)
  {  }

  //
  // member functions for OBUnitCell class
  //
  OBUnitCell::OBUnitCell():
    OBGenericData("UnitCell", OBGenericDataType::UnitCell),
    m_a(0.0), m_b(0.0), m_c(0.0), m_alpha(0.0), m_beta(0.0), m_gamma(0.0),
    m_spaceGroup( NULL ), m_lattice(Undefined)
  {  }

  OBUnitCell::OBUnitCell(const OBUnitCell &src) :
    OBGenericData("UnitCell", OBGenericDataType::UnitCell),
    m_a(src.m_a), m_b(src.m_b), m_c(src.m_c), 
    m_alpha(src.m_alpha), m_beta(src.m_beta), m_gamma(src.m_gamma),
    m_offset(src.m_offset),
    m_v1(src.m_v1), m_v2(src.m_v2), m_v3(src.m_v3),
    m_spaceGroupName(src.m_spaceGroupName),
    m_spaceGroup(src.m_spaceGroup),
    m_lattice(src.m_lattice)
  {  }

  OBUnitCell & OBUnitCell::operator=(const OBUnitCell &src)
  {
    if(this == &src)
      return(*this);

    m_a = src.m_a;
    m_b = src.m_b;
    m_c = src.m_c;
    m_alpha = src.m_alpha;
    m_beta = src.m_beta;
    m_gamma = src.m_gamma;
    m_offset = src.m_offset;

    m_v1 = src.m_v1;
    m_v2 = src.m_v2;
    m_v3 = src.m_v3;

    m_spaceGroup = src.m_spaceGroup;
    m_spaceGroupName = src.m_spaceGroupName;
    m_lattice = src.m_lattice;

    return(*this);
  }

  /*!
  ** The angles and lengths of the unitcell will be calculated from the
  ** vectors @p v1, @p v2 and @p v3. Those vectors will as well be
  ** stored internally.
  **Implements <a href="http://qsar.sourceforge.net/dicts/blue-obelisk/index.xhtml#convertCartesianIntoNotionalCoordinates">blue-obelisk:convertCartesianIntoNotionalCoordinates</a>
  **\brief Sets the vectors, angles and lengths of the unitcell
  **\param v1 The x-vector
  **\param v2 The y-vector
  **\param v3 The z-vector
  **\see OBUnitCell::GetCellVectors
  */
  void OBUnitCell::SetData(const Eigen::Vector3d v1, const Eigen::Vector3d v2, const Eigen::Vector3d v3)
  {
    m_a = v1.norm();
    m_b = v2.norm();
    m_c = v3.norm();
    m_alpha = VectorAngle(v2, v3);
    m_beta = VectorAngle(v1, v3);
    m_gamma = VectorAngle(v1, v2);
    m_v1 = v1;
    m_v2 = v2;
    m_v3 = v3;
  }

  //! Implements <a href="http://qsar.sourceforge.net/dicts/blue-obelisk/index.xhtml#convertNotionalIntoCartesianCoordinates">blue-obelisk:convertNotionalIntoCartesianCoordinates</a>
  vector<Eigen::Vector3d> OBUnitCell::GetCellVectors()
  {
    vector<Eigen::Vector3d> v;
    v.reserve(3);

    if (IsNegligible(m_v1.norm(), 1.0, 1.0e-9) &&
        IsNegligible(m_v2.norm(), 1.0, 1.0e-9) &&
        IsNegligible(m_v3.norm(), 1.0, 1.0e-9))
      {
        Eigen::Vector3d temp;
        Eigen::Matrix3d m = GetOrthoMatrix();

        temp = Eigen::Vector3d(1.0, 0.0, 0.0);
        v.push_back(m * temp);
        temp = Eigen::Vector3d(0.0, 1.0, 0.0);
        v.push_back(m * temp);
        temp = Eigen::Vector3d(0.0, 0.0, 1.0);
        v.push_back(m * temp);
      }
    else
      {
        v.push_back(m_v1);
        v.push_back(m_v2);
        v.push_back(m_v3);
      }

    return v;
  }

  Eigen::Matrix3d OBUnitCell::GetCellMatrix()
  {
    Eigen::Matrix3d m;

    if (IsNegligible(m_v1.norm(), 1.0, 1.0e-9) &&
        IsNegligible(m_v2.norm(), 1.0, 1.0e-9) &&
        IsNegligible(m_v3.norm(), 1.0, 1.0e-9))
      {
        m = GetOrthoMatrix();
      }
    else
      {
        m.row(0) = m_v1;
        m.row(1) = m_v2;
        m.row(2) = m_v3;
      }
    return m;
  }

  //! Implements <a href="http://qsar.sourceforge.net/dicts/blue-obelisk/index.xhtml#calculateOrthogonalisationMatrix">blue-obelisk:calculateOrthogonalisationMatrix</a>
  Eigen::Matrix3d OBUnitCell::GetOrthoMatrix()
  {
    Eigen::Matrix3d m;
    double V;

    m_alpha *= DEG_TO_RAD;
    m_beta  *= DEG_TO_RAD;
    m_gamma *= DEG_TO_RAD;

    // from the PDB specification:
    //  http://www.rcsb.org/pdb/docs/format/pdbguide2.2/part_75.html

    // since we'll ultimately divide by (a * b), we've factored those out here
    V = m_c * sqrt(1 - SQUARE(cos(m_alpha)) - SQUARE(cos(m_beta)) - SQUARE(cos(m_gamma))
                 + 2 * cos(m_alpha) * cos(m_beta) * cos(m_gamma));

    m(0,0) = m_a;
    m(0,1) = m_b * cos(m_gamma);
    m(0,2) = m_c * cos(m_beta);

    m(1,0) = 0.0;
    m(1,1) = m_b * sin(m_gamma);
    m(1,2) = m_c * ( cos(m_alpha) - cos(m_beta) * cos(m_gamma) ) / sin(m_gamma);

    m(2,0) = 0.0;
    m(2,1) = 0.0;
    m(2,2) = V / (sin(m_gamma)); // again, we factored out A * B when defining V

    return m;
  }

  // Based on code in PyMMLib: http://pymmlib.sf.net/
  //! Matrix to convert from Cartesian to fractional
  //! Implements <a href="http://qsar.sourceforge.net/dicts/blue-obelisk/index.xhtml#convertCartesianIntoFractionalCoordinates">blue-obelisk:convertCartesianIntoFractionalCoordinates</a> 
  Eigen::Matrix3d OBUnitCell::GetFractionalMatrix()
  {
    Eigen::Matrix3d m;
    double sinAlpha, sinBeta, sinGamma;
    double cosAlpha, cosBeta, cosGamma;
    double v;

    sinAlpha = sin(m_alpha * DEG_TO_RAD);
    sinBeta = sin(m_beta * DEG_TO_RAD);
    sinGamma = sin(m_gamma * DEG_TO_RAD);
    cosAlpha = cos(m_alpha * DEG_TO_RAD);
    cosBeta = cos(m_beta * DEG_TO_RAD);
    cosGamma = cos(m_gamma * DEG_TO_RAD);

    v = sqrt(1 - SQUARE(cosAlpha) - SQUARE(cosBeta) - SQUARE(cosGamma) +
             2 * cosAlpha*cosBeta*cosGamma);

    m(0,0) = 1.0 / m_a;
    m(0,1) = -cosGamma / (m_a * sinGamma);
    m(0,2) = (cosGamma * cosAlpha - cosBeta) / (m_a * v * sinGamma);
    m(1,0) = 0.0;
    m(1,1) = 1.0 / (m_b * sinGamma);
    m(1,2) = (cosGamma * cosBeta - cosAlpha) / (m_b * v * sinGamma);
    m(2,0) = 0.0;
    m(2,1) = 0.0;
    m(2,2) = sinGamma / (m_c * v) ;

    return m;
  }

  OBUnitCell::LatticeType OBUnitCell::GetLatticeType( int spacegroup )
  {
	  //	1-2 	Triclinic
	  //	3-15 	Monoclinic
	  //	16-74	Orthorhombic
	  //	75-142 	Tetragonal
	  //	143-167 Rhombohedral
	  //	168-194 Hexagonal
	  //	195-230 Cubic

      if ( spacegroup == 0  && m_spaceGroup)
          spacegroup = m_spaceGroup->GetId();
	  
	  if ( spacegroup <= 0 )
		  return OBUnitCell::Undefined;

	  else if ( spacegroup == 1 ||
              spacegroup == 2 )
		  return OBUnitCell::Triclinic;
	  
	  else if ( spacegroup >= 3 ||
              spacegroup <= 15 )
		  return OBUnitCell::Monoclinic;
	  
	  else if ( spacegroup >= 16 ||
              spacegroup <= 74 )
		  return OBUnitCell::Orthorhombic;
	  
	  else if ( spacegroup >= 75 ||
              spacegroup <= 142 )
		  return OBUnitCell::Tetragonal;
	  
	  else if ( spacegroup >= 143 ||
              spacegroup <= 167 )
		  return OBUnitCell::Rhombohedral;
	  
	  else if ( spacegroup >= 168 ||
              spacegroup <= 194 )
		  return OBUnitCell::Hexagonal;
	  
	  else if ( spacegroup >= 195 ||
              spacegroup <= 230 )
		  return OBUnitCell::Cubic;

	  //just to be extra sure
	  else // ( spacegroup > 230 )
		  return OBUnitCell::Undefined;
  }
  
  OBUnitCell::LatticeType OBUnitCell::GetLatticeType()
  {
    if (m_lattice != Undefined)
      return m_lattice;
    else if (m_spaceGroup != NULL)
      return GetLatticeType(m_spaceGroup->GetId());

    unsigned int rightAngles = 0;
    if (IsApprox(m_alpha, 90.0, 1.0e-3)) rightAngles++;
    if (IsApprox(m_beta,  90.0, 1.0e-3)) rightAngles++;
    if (IsApprox(m_gamma, 90.0, 1.0e-3)) rightAngles++;

    switch (rightAngles)
      {
      case 3:
        if (IsApprox(m_a, m_b, 1.0e-4) && IsApprox(m_b, m_c, 1.0e-4))
          m_lattice = Cubic;
        else if (IsApprox(m_a, m_b, 1.0e-4) || IsApprox(m_b, m_c, 1.0e-4))
          m_lattice = Tetragonal;
        else
          m_lattice = Orthorhombic;
        break;
      case 2:
        if ( (IsApprox(m_alpha, 120.0, 1.0e-3) 
              || IsApprox(m_beta, 120.0, 1.0e-3) 
              || IsApprox(m_gamma, 120.0f, 1.0e-3))
             && (IsApprox(m_a, m_b, 1.0e-4) || IsApprox(m_b, m_c, 1.0e-4)) )
          m_lattice = Hexagonal;
        else
          m_lattice = Monoclinic;
        break;
      default:
        if (IsApprox(m_a, m_b, 1.0e-4) && IsApprox(m_b, m_c, 1.0e-4))
          m_lattice = Rhombohedral;
        else
          m_lattice = Triclinic;
      }

    return m_lattice;
  }

  int OBUnitCell::GetSpaceGroupNumber( std::string name)
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
        if (m_spaceGroup != NULL)
          return m_spaceGroup->GetId();
        else
          name = m_spaceGroupName;
      }
    static const int numStrings = sizeof( spacegroups ) / sizeof( spacegroups[0] );
    for ( int i = 0; i < numStrings; ++i ) {
      if (name == spacegroups[i] ) {
        return i+1;
      }
    }
    return 0; //presumably never reached
  }
  
  void OBUnitCell::FillUnitCell(OBMol *mol)
  {
    const SpaceGroup *sg = GetSpaceGroup(); // the actual space group and transformations for this unit cell
    
    // For each atom, we loop through: convert the coords back to inverse space, apply the transformations and create new atoms
    Eigen::Vector3d uniqueV, newV;
    list<Eigen::Vector3d> transformedVectors; // list of symmetry-defined copies of the atom
    list<Eigen::Vector3d>::iterator transformIterator;
    OBAtom *newAtom;
    list<OBAtom*> atoms; // keep the current list of unique atoms -- don't double-create
    FOR_ATOMS_OF_MOL(atom, *mol)
      atoms.push_back(&(*atom));

    list<OBAtom*>::iterator i;
    for (i = atoms.begin(); i != atoms.end(); ++i) {
      uniqueV = (*i)->GetVector();
      uniqueV = GetFractionalMatrix() * uniqueV;
        
      transformedVectors = sg->Transform(uniqueV);
      for (transformIterator = transformedVectors.begin();
           transformIterator != transformedVectors.end(); ++transformIterator) {
        // coordinates are in reciprocal space -- check if it's in the unit cell
        // TODO: transform these into the unit cell and check for duplicates
        if (transformIterator->x() < 0.0 || transformIterator->x() > 1.0)
          continue;
        else if (transformIterator->y() < 0.0 || transformIterator->y() > 1.0)
          continue;
        else if (transformIterator->z() < 0.0 || transformIterator->z() > 1.0)
          continue;
             
        newAtom = mol->NewAtom();
        newAtom->Duplicate(*i);
        newAtom->SetVector(GetOrthoMatrix() * (*transformIterator));
      } // end loop of transformed atoms
    } // end loop of atoms
  }
  
  double OBUnitCell::GetCellVolume()
  {
    double result = 0.0;
    
    switch ( GetLatticeType() )
      {
      case Triclinic:
        result = m_a * m_b * m_c 
          * sqrt(1
                 - SQUARE(cos( m_alpha ))
                 - SQUARE(cos( m_beta ))
                 - SQUARE(cos( m_gamma ))
                 + 2 * cos( m_alpha ) * cos( m_beta ) * cos( m_gamma )
                 );
        break;
      case Monoclinic:
        result = m_a * m_b * m_c * sin( m_beta );
        break;
      case Orthorhombic:
        result = m_a * m_b * m_c;
        break;
      case Tetragonal:
        result = m_a * m_a * m_c;
        break;
      case Rhombohedral:
        result = m_a * m_a * m_a
          * sqrt(1
                 - SQUARE(cos( m_alpha ))
                 - SQUARE(cos( m_beta ))
                 - SQUARE(cos( m_gamma ))
                 + 2 * cos( m_alpha ) * cos( m_beta ) * cos( m_gamma )
                 );
        break;
      case Hexagonal:
        result = pow( 3.0, 0.333333333 ) * m_a * m_a * m_c / 2;
        break;
      case Cubic:
        result = m_a * m_a * m_a;
        break;
      default:
        result = 0.0;
      }
    
    return result;
  }
  
  //
  // member functions for OBSymmetryData class
  //
  OBSymmetryData::OBSymmetryData(): 
    OBGenericData("Symmetry", OBGenericDataType::SymmetryData)
  { }

  OBSymmetryData::OBSymmetryData(const OBSymmetryData &src) :
    OBGenericData(src._attr, src._type, src._source), 
    m_pointGroup(src.m_pointGroup), m_spaceGroup(src.m_spaceGroup)
  {  }

  OBSymmetryData & OBSymmetryData::operator=(const OBSymmetryData &src)
  {
    if(this == &src)
      return(*this);

    m_pointGroup = src.m_pointGroup;
    m_spaceGroup = src.m_spaceGroup;
    _source = src._source;

    return(*this);
  }

  OBConformerData::OBConformerData() :
    OBGenericData("Conformers", OBGenericDataType::ConformerData)
  {  }

  OBConformerData::OBConformerData(const OBConformerData &src) :
    OBGenericData("Conformers", OBGenericDataType::ConformerData),
    m_vDimension(src.m_vDimension),
    m_vEnergies(src.m_vEnergies), m_vForces(src.m_vForces),
    m_vVelocity(src.m_vVelocity), m_vDisplace(src.m_vDisplace),
    m_vData(src.m_vData)
  {  }

  OBConformerData & OBConformerData::operator=(const OBConformerData &src)
  {
    if(this == &src)
      return(*this);
    
    _source = src._source;

    m_vDimension = src.m_vDimension;
    m_vEnergies = src.m_vEnergies;
    m_vForces = src.m_vForces;
    m_vVelocity = src.m_vVelocity;
    m_vDisplace = src.m_vDisplace;
    m_vData = src.m_vData;

    return(*this);
  }

  //
  //member functions for OBRingData class
  //

  OBRingData::OBRingData() :
    OBGenericData("RingData", OBGenericDataType::RingData)
  {
    m_vr.clear();
  }

  /*!
  **\brief OBRingData copy constructor
  **\param src reference to original OBRingData object (rhs)
  */
  OBRingData::OBRingData(const OBRingData &src)
    :	OBGenericData(src),	//chain to base class
      m_vr(src.m_vr)				//chain to member classes
  {
    //no other memeber data
    //memory management

    vector<OBRing*>::iterator ring;

    for(ring = m_vr.begin();ring != m_vr.end();++ring)
      {
        OBRing *newring = new OBRing;
        (*newring) = (**ring);	//copy data to new object
        (*ring)    = newring;	//repoint new pointer to new copy of data
      }
  }

  OBRingData::~OBRingData()
  {
    vector<OBRing*>::iterator ring;
    for (ring = m_vr.begin();ring != m_vr.end();++ring)
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
    for(ring = m_vr.begin();ring != m_vr.end();++ring)
      {
        delete &*ring;	//deallocate old rings to prevent memory leak
      }

    m_vr.clear();
    m_vr = src.m_vr;	//copy vector properties

    for(ring = m_vr.begin();ring != m_vr.end();++ring)
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
    i = m_vr.begin();
    return((i == m_vr.end()) ? (OBRing*)NULL : (OBRing*)*i);
  }

  OBRing *OBRingData::NextRing(std::vector<OBRing*>::iterator &i)
  {
    i = m_vr.begin();
    return((i == m_vr.end()) ? (OBRing*)NULL : (OBRing*)*i);
  }

  //
  //member functions for OBAngle class - stores all angles
  //

  OBAngle::OBAngle():
    m_vertex(NULL), m_termini(NULL, NULL), m_radians(0.0)
  {  }

  OBAngle::OBAngle(OBAtom *vertex,OBAtom *a,OBAtom *b):
    m_vertex(vertex), m_termini(a, b)
  {
    SortByIndex();
  }

  OBAngle::OBAngle(const OBAngle &src):
    m_vertex(src.m_vertex), m_termini(src.m_termini), m_radians(src.m_radians)
  {  }

  OBAngle& OBAngle::operator = (const OBAngle &src)
  {
    if (this == &src)
      return(*this);

    m_vertex         = src.m_vertex;
    m_termini.first  = src.m_termini.first;
    m_termini.second = src.m_termini.second;
    m_radians        = src.m_radians;

    return(*this);
  }

  void OBAngle::Clear()
  {
    m_vertex         = 0;
    m_termini.first  = 0;
    m_termini.second = 0;
    m_radians        = 0.0;
    return;
  }

  void OBAngle::SetAtoms(OBAtom *vertex,OBAtom *a,OBAtom *b)
  {
    m_vertex         = vertex;
    m_termini.first  = a;
    m_termini.second = b;
    SortByIndex();
    return;
  }

  void OBAngle::SetAtoms(triple<OBAtom*,OBAtom*,OBAtom*> &atoms)
  {
    m_vertex         = atoms.first;
    m_termini.first  = atoms.second;
    m_termini.second = atoms.third;
    SortByIndex();
    return;
  }

  triple<OBAtom*,OBAtom*,OBAtom*> OBAngle::GetAtoms()
  {
    triple<OBAtom*,OBAtom*,OBAtom*> atoms;
    atoms.first  = m_vertex;
    atoms.second = m_termini.first;
    atoms.third  = m_termini.second;
    return(atoms);
  }

  void OBAngle::SortByIndex()
  {
    OBAtom *tmp;

    if(m_termini.first->GetIdx() > m_termini.second->GetIdx())
      {
        tmp              = m_termini.first;
        m_termini.first  = m_termini.second;
        m_termini.second = tmp;
      }
  }

  bool OBAngle::operator ==(const OBAngle &other)
  {
    return ((m_vertex         == other.m_vertex)        &&
            (m_termini.first  == other.m_termini.first) &&
            (m_termini.second == other.m_termini.second));
  }

  //
  //member functions for OBAngleData class - stores OBAngle set
  //

  OBAngleData::OBAngleData()
    :	OBGenericData("AngleData", OBGenericDataType::AngleData)
  {  }

  OBAngleData::OBAngleData(const OBAngleData &src)
    :	OBGenericData(src), m_angles(src.m_angles)
  {  }

  OBAngleData& OBAngleData::operator =(const OBAngleData &src)
  {
    if (this == &src)
      return(*this);

    _source = src._source;
    m_angles = src.m_angles;

    return(*this);
  }

  void OBAngleData::Clear()
  {
    m_angles.clear();
    return;
  }

  void OBAngleData::SetData(OBAngle &angle)
  {
    m_angles.push_back(angle);
    return;
  }

  bool OBAngleData::FillAngleArray(std::vector<std::vector<unsigned int> > &angles)
  {
    if(m_angles.empty())
      return(false);

    vector<OBAngle>::iterator angle;
    
    angles.clear();
    angles.resize(m_angles.size());

    unsigned int ct = 0;

    for( angle = m_angles.begin(); angle != m_angles.end(); angle++,ct++)
      {
        angles[ct].resize(3);
        angles[ct][0] = angle->m_vertex->GetIdx() - 1;
        angles[ct][1] = angle->m_termini.first->GetIdx() - 1;
        angles[ct][2] = angle->m_termini.second->GetIdx() - 1;
      }
 
    return(true);
  }
  
 unsigned int OBAngleData::FillAngleArray(int **angles, unsigned int &size)
  {
    if(m_angles.size() > size)
      {
        delete [] *angles;
        *angles = new int[m_angles.size()*3];
        size    = (unsigned int)m_angles.size();
      }

    vector<OBAngle>::iterator angle;
    int angleIdx = 0;
    for( angle=m_angles.begin(); angle!=m_angles.end(); ++angle)
      {
        *angles[angleIdx++] = angle->m_vertex->GetIdx();
        *angles[angleIdx++] = angle->m_termini.first->GetIdx();
        *angles[angleIdx++] = angle->m_termini.second->GetIdx();
      }
    return (unsigned int)m_angles.size();
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
    m_ads.push_back(ad);

    m_bc.first  = b;
    m_bc.second = c;
  }

  /*!
  **\brief OBTorsion copy constructor
  */
  OBTorsion::OBTorsion(const OBTorsion &src)
    :	m_bc(src.m_bc), m_ads(src.m_ads)
  {}

  /*!
  **\brief Returns all the 4 atom sets in OBTorsion
  */
  vector<quad<OBAtom*,OBAtom*,OBAtom*,OBAtom*> > OBTorsion::GetTorsions()
  {
    quad<OBAtom*,OBAtom*,OBAtom*,OBAtom*> abcd;

    abcd.second = m_bc.first;
    abcd.third  = m_bc.second;

    vector<quad<OBAtom*,OBAtom*,OBAtom*,OBAtom*> > torsions;
    vector<triple<OBAtom*,OBAtom*,double> >::iterator ad;

    for(ad = m_ads.begin();ad != m_ads.end();++ad)
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

    m_bc  = src.m_bc;
    m_ads = src.m_ads;

    return(*this);
  }

  void OBTorsion::Clear()
  {
    m_bc.first  = 0;
    m_bc.second = 0;
    m_ads.erase(m_ads.begin(), m_ads.end());
  }

  bool OBTorsion::SetAngle(double radians,unsigned int index)
  {
    if(index >= m_ads.size())
      return(false);

    m_ads[index].third = radians;

    return(true);
  }

  bool OBTorsion::GetAngle(double &radians, unsigned int index)
  {
    if(index >= m_ads.size())
      return false;
    radians = m_ads[index].third;
    return true;
  }

  unsigned int OBTorsion::GetBondIdx()
  {
    return(m_bc.first->GetBond(m_bc.second)->GetIdx());
  }

  bool OBTorsion::IsProtonRotor()
  {
    bool Aprotor = true;
    bool Dprotor = true;
    vector<triple<OBAtom*,OBAtom*,double> >::iterator ad;
    for(ad = m_ads.begin();ad != m_ads.end() && (Aprotor || Dprotor);++ad)
      {
        if(!ad->first->IsHydrogen())
          Aprotor = false;
        if(!ad->second->IsHydrogen())
          Dprotor = false;
      }
    return (Aprotor || Dprotor);
  }

  bool OBTorsion::AddTorsion(OBAtom *a,OBAtom *b, OBAtom *c,OBAtom *d)
  {
    if(!Empty() && (b != m_bc.first || c != m_bc.second))
      return(false);

    if(Empty())
      {
        m_bc.first  = b;
        m_bc.second = c;
      }

    triple<OBAtom*,OBAtom*,double> ad(a,d,0.0);
    m_ads.push_back(ad);

    return(true);
  }

  bool OBTorsion::AddTorsion(quad<OBAtom*,OBAtom*,OBAtom*,OBAtom*> &atoms)
  {
    if(!Empty() && (atoms.second != m_bc.first || atoms.third != m_bc.second))
      return(false);

    if(Empty())
      {
        m_bc.first  = atoms.second;
        m_bc.second = atoms.third;
      }

    triple<OBAtom*,OBAtom*,double> ad(atoms.first,atoms.fourth,0.0);
    m_ads.push_back(ad);

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
    :	OBGenericData(src), m_torsions(src.m_torsions)
  {  }

  OBTorsionData& OBTorsionData::operator =(const OBTorsionData &src)
  {
    if (this == &src)
      return(*this);

    OBGenericData::operator =(src);

    _source = src._source;
    m_torsions = src.m_torsions;

    return(*this);
  }

  void OBTorsionData::Clear()
  {
    m_torsions.clear();
  }

  void OBTorsionData::SetData(OBTorsion &torsion)
  {
    m_torsions.push_back(torsion);
  }

  bool OBTorsionData::FillTorsionArray(std::vector<std::vector<unsigned int> > &torsions)
  {
    if(m_torsions.empty())
      return(false);

    vector<quad<OBAtom*,OBAtom*,OBAtom*,OBAtom*> > tmpquads,quads;
    vector<quad<OBAtom*,OBAtom*,OBAtom*,OBAtom*> >::iterator thisQuad;
    vector<OBTorsion>::iterator torsion;

    //generate set of all 4 atom abcd's from torsion structure
    for (torsion = m_torsions.begin();torsion != m_torsions.end();++torsion)
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
  // Member functions for OBChiralDarta
  //
  bool OBChiralData::SetAtom4Refs(std::vector<unsigned int> atom4refs, atomreftype t)
  {
    if (atom4refs.size() > 4) {
      obErrorLog.ThrowError(__FUNCTION__, "More than 4 atoms in atom4refs", obDebug);
      return(false);
    }
    
    switch (t) {
      case input: 
        m_atom4refs = atom4refs;
        break;
      case output:
        m_atom4refo = atom4refs;
        break;
      case calcvolume:
        m_atom4refc = atom4refs;
        break;
      default: 
        obErrorLog.ThrowError(__FUNCTION__, "AtomRefType called is invalid", obDebug);
        return(false);
    }

    return (true);
  }

  int OBChiralData::AddAtomRef(unsigned int atomref, atomreftype t)
  {
    switch (t) {
      case input: 
        m_atom4refs.push_back(atomref);
        break;
      case output: 
        m_atom4refo.push_back(atomref);
        break;
      case calcvolume:
        m_atom4refc.push_back(atomref);
        break;
      default:
        obErrorLog.ThrowError(__FUNCTION__, "AtomRefType called is invalid", obDebug);
        return(false);
    }
              
    return (m_atom4refs.size());
  }

  unsigned int OBChiralData::GetAtomRef(int a, atomreftype t)
  {
    switch (t) {
    case input: 
      return(m_atom4refs[a]);
      break;
    case output: 
      return(m_atom4refo[a]);
      break;
    case calcvolume: 
      return(m_atom4refc[a]);
      break;
    default:
      obErrorLog.ThrowError(__FUNCTION__, "AtomRefType called is invalid", obDebug);
      return(false);
    }  
  }

  std::vector<unsigned int> OBChiralData::GetAtom4Refs(atomreftype t) const
  {
    switch (t) {
      case output:
        return(m_atom4refo);
        break;
      case input:
        return(m_atom4refs);
        break;
      case calcvolume:
        return(m_atom4refc);
        break;
      default:
        obErrorLog.ThrowError(__FUNCTION__, "AtomRefType called is invalid", obDebug);
        return(m_atom4refo);
    }
  }

  unsigned int OBChiralData::GetSize(atomreftype t) const
  {
    switch (t) {
      case output:
        return(unsigned int)m_atom4refo.size();
        break;
      case input:
        return(unsigned int)m_atom4refs.size();
        break;
      case calcvolume:
        return(unsigned int)m_atom4refc.size();
      default:
        obErrorLog.ThrowError(__FUNCTION__, "AtomRefType called is invalid", obDebug);
        return(0);
      }
  }

  OBChiralData::OBChiralData()
    : OBGenericData("ChiralData", OBGenericDataType::ChiralData, perceived)
  {  }

  OBChiralData::OBChiralData(const OBChiralData &src)
    : OBGenericData(src)
  {
    m_atom4refs = src.m_atom4refs;
    m_atom4refo = src.m_atom4refo;
    m_atom4refc = src.m_atom4refc;
    m_parity    = src.m_parity;
  }

  OBChiralData & OBChiralData::operator=(const OBChiralData &src)
  {
    if(this == &src)
      return(*this);

    _source = src._source;

    m_atom4refs = src.m_atom4refs;
    m_atom4refo = src.m_atom4refo;
    m_atom4refc = src.m_atom4refc;
    m_parity    = src.m_parity;
    return(*this);
  }

  void OBChiralData::Clear()
  {
    m_atom4refs.clear();
    m_parity = 0;
    m_atom4refo.clear();
    m_atom4refc.clear();
  }

//
//member functions for OBVibrationData class
//


void OBVibrationData::SetData(const std::vector< std::vector< Eigen::Vector3d > > & vLx,
                              const std::vector<double> & vFrequencies,
                              const std::vector<double> & vIntensities)
{
  this->m_vLx = vLx;
  this->m_vFrequencies = vFrequencies;
  this->m_vIntensities = vIntensities;
}


unsigned int OBVibrationData::GetNumberOfFrequencies() const
{
  return !this->m_vFrequencies.empty() ? this->m_vFrequencies.size() - 1 : 0;
}

} //end namespace OpenBabel

//! \file generic.cpp
//! \brief Handle OBGenericData classes. Custom data for atoms, bonds, etc.

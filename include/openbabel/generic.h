/**********************************************************************
generic.h - Handle generic data classes. Custom data for atoms, bonds, etc.
 
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2006 by Geoffrey R. Hutchison
Some portions Copyright (C) 2008 Tim Vandermeersch
 
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

#include <openbabel/babelconfig.h>

#include <string>
#include <vector>
#include <map>

#include <openbabel/math/spacegroup.h>
#include <openbabel/obutil.h>
#include <openbabel/base.h>

namespace OpenBabel
{
  // Forward declarations
  class OBBase;
  class OBAtom;
  class OBBond;
  class OBMol;
  class OBRing;

  /// @addtogroup data Data classes 
  //@{

  /** @class OBCommentData generic.h <openbabel/generic.h>
      @brief Used to store a comment string (can be multiple lines long).
   */
  class OBAPI OBCommentData : public OBGenericData
  {
    protected:
      std::string m_data; //!< The comment.
    public:
      /**
       * Constructor.
       */
      OBCommentData();
      /**
       * Copy constructor.
       */
      OBCommentData(const OBCommentData&);
      /**
       * Clone this data.
       */
      virtual OBGenericData* Clone(OBBase* /*parent*/) const
      { return new OBCommentData(*this); }
      /**
       * Assign operator.
       */
      OBCommentData& operator=(const OBCommentData &src);
      /**
       * Set the comment.
       */
      void SetData(const std::string &data) { m_data = data; Trim(m_data); }
      /**
       * Set the comment.
       */
      void SetData(const char *d) { m_data = d; Trim(m_data); }
      /**
       * Get the comment.
       */
      const std::string &GetData() const { return m_data; }
      /**
       * Get the comment.
       */
      const std::string &GetValue() const { return m_data; }
  };

  //@}

  /** @class OBExternalBond generic.h <openbabel/generic.h>
      @brief Used to store information on an external bond 
      (e.g., SMILES fragments).
   */
  class OBAPI OBExternalBond
  {
    protected:
      int     m_idx;  //!< Index for the bond.
      OBAtom *m_atom; //!< Start atom for the bond.
      OBBond *m_bond; //!< The bond.
    public:
      /**
       * Constructor.
       */
      OBExternalBond(): m_idx(0), m_atom(NULL), m_bond(NULL) {}
      /**
       * Constructor.
       */ 
      OBExternalBond(OBAtom *, OBBond *, int);
      /**
       * Copy constructor.
       */
      OBExternalBond(const OBExternalBond &);
      /**
       * Destructor.
       */
      ~OBExternalBond() {}
      /**
       * Get the bond index.
       */
      int GetIdx() const { return m_idx; }
      /**
       * Get the start atom for the bond.
       */
      OBAtom *GetAtom() const { return m_atom; }
      /**
       * Get the bond.
       */
      OBBond *GetBond() const { return m_bond; }
      /**
       * Set the bond index.
       */
      void SetIdx(int idx) { m_idx = idx; }
      /**
       * Set the start atom for the bond.
       */
      void SetAtom(OBAtom *atom) { m_atom = atom; }
      /**
       * Set the bond.
       */
      void SetBond(OBBond *bond) { m_bond = bond; }
  };

  /// @addtogroup data Data classes 
  //@{

  /** @class OBExternalBondData generic.h <openbabel/generic.h>
      @brief Used to store information on external bonds (e.g., in SMILES fragments).

      This class is used by the smiles format code to temporarely store 
      external bonds in a molecule.
   */
  class OBAPI OBExternalBondData : public OBGenericData
  {
    protected:
      std::vector<OBExternalBond> m_vexbnd; //!< vector of OBExternalBond objects
    public:
      /**
       * Constructor.
       */
      OBExternalBondData();
      /**
       * Copying is not used and too much work to set up. 
       * @return NULL
       */
      virtual OBGenericData* Clone(OBBase* /*parent*/) const { return NULL; }
      /**
       * Add a new OBExtenalBond to the list.
       */
      void SetData(OBAtom*,OBBond*,int);
      /**
       * Get the list of OBExternalBond objects.
       */
      std::vector<OBExternalBond> *GetData() { return &m_vexbnd; }
  };

  /** @class OBPairData generic.h <openbabel/generic.h>
      @brief Used to store arbitrary text attribute/value relationships.

      Ideal for arbitrary text descriptors for molecules, atoms, bonds, residues,
      e.g. in QSAR.
   */
  class OBAPI OBPairData : public OBGenericData
  {
    protected:
      std::string m_value; //!< The data for this key/value pair
    public:
      /**
       * Constructor.
       */
      OBPairData();
      /**
       * Clone this data.
       */
      virtual OBGenericData* Clone(OBBase* /*parent*/) const
      { return new OBPairData(*this); }
      /**
       * Set the value.
       */
      void SetValue(const char *v)        { m_value = v; }
      /**
       * Set the value.
       */
      void SetValue(const std::string &v) { m_value = v; }
      /**
       * Get the value.
       */
      const std::string &GetValue() const { return m_value; }
  };

  /** @class OBPairTemplate generic.h <openbabel/generic.h>
      @brief Used to store arbitrary attribute/value relationsips of any type.
      
      More detailed description in generic.cpp
   */
  template <class ValueT>
  class OBAPI OBPairTemplate : public OBGenericData
  {
    protected:
      ValueT m_value; //!< The data for this key/value pair
    public:
      /**
       * Constructor.
       */
    OBPairTemplate() :
        OBGenericData("PairData", OBGenericDataType::PairData) 
    {  }
    /**
     * Set the value.
     */
    void SetValue(const ValueT t) { m_value = t; }
    /**
     * Get the data.
     */
    const ValueT &GetGenericValue() const { return m_value; }
  };

  //@}

  //! Store arbitrary key/value integer data like OBPairData
  typedef OBPairTemplate<int>     OBPairInteger;
  //! Store arbitrary key/value floating point data like OBPairData
  typedef OBPairTemplate<double>  OBPairFloatingPoint;

  /// @addtogroup data Data classes 
  //@{

  /** @class OBSetData generic.h <openbabel/generic.h>
      @brief Used to store arbitrary attribute/set relationships.
      
      Should be used to store a set of OBGenericData based on an attribute.
   */
  class OBAPI OBSetData : public OBGenericData
  {
    protected:
      std::vector<OBGenericData *> m_vdata; //!< all the data
    public:
      /**
       * Constructor.
       */
      OBSetData() : OBGenericData("SetData", OBGenericDataType::SetData) {}
      /**
       * Clone this data.
       */
      virtual OBGenericData* Clone(OBBase* /*parent*/) const
      { return new OBSetData(*this); }
      /**
       * Add an OBGenericData element to the set.
       */
      void AddData(OBGenericData *d) { if(d) m_vdata.push_back(d); }
      /** 
       * Set the array of data to a new vector.
       */
      void SetData(std::vector<OBGenericData *> &vdata) { m_vdata = vdata; }
      /**
       * @return the OBGenericData associate with the attribute name parameter.
       */
      OBGenericData *GetData(const char *s) 
      {
        std::vector<OBGenericData*>::iterator i;

        for (i = m_vdata.begin();i != m_vdata.end();++i)
          if ((*i)->GetAttribute() == s)
            return *i;

        return NULL;
      }
      /**
       * @return the OBGenericData associate with the attribute name parameter.
       */
      OBGenericData *GetData(const std::string &s)
      {
        std::vector<OBGenericData*>::iterator i;

        for (i = m_vdata.begin();i != m_vdata.end();++i)
          if ((*i)->GetAttribute() == s)
            return *i;

        return NULL;
      }
      /**
       * Gets the entire set.
       */
      virtual const std::vector<OBGenericData *> &GetData() const //now virtual and const
      { return m_vdata; }
      /**
       * Get the begin iterator.
       */
      std::vector<OBGenericData*>::iterator GetBegin() { return m_vdata.begin(); }
      /**
       * Get the end iterator.
       */
      std::vector<OBGenericData*>::iterator GetEnd() { return m_vdata.end(); }
      /**
       * Delete the matching OBGenericData element.
       */
      void DeleteData(OBGenericData *gd)
      {
        std::vector<OBGenericData*>::iterator i;
        for (i = m_vdata.begin();i != m_vdata.end();++i)
          if (*i == gd) {
            delete *i;
            m_vdata.erase(i);
          }
      }
  }; // OBSetData

  /** @class OBVirtualBond generic.h <openbabel/generic.h>
      @brief Used to temporarily store bonds that reference
      an atom that has not yet been added to a molecule
   */
  class OBAPI OBVirtualBond : public OBGenericData
  {
    protected:
      int m_bgn;    //!< begin atom index
      int m_end;    //!< end atom index
      int m_ord;    //!< bond order
      int m_stereo; //!< bond stereo chemistry
    public:
      /**
       * Constructor.
       */
      OBVirtualBond();
      /**
       * Clone this data.
       */
      virtual OBGenericData* Clone(OBBase* /*parent*/) const
      { return new OBVirtualBond(*this); }
      /**
       * Constructor.
       */
      OBVirtualBond(int,int,int,int stereo=0);
      /**
       * Get the begin atom index.
       */
      int GetBgn() { return m_bgn; }
      /**
       * Get the end atom index.
       */
      int GetEnd() { return m_end; }
      /**
       * Get the bond order.
       */
      int GetOrder() { return m_ord; }
      /**
       * Get the bond stereo chemistry.
       */
      int GetStereo() { return m_stereo; }
  };

  /** @class OBRingData generic.h <openbabel/generic.h>
      @brief Used to store the SSSR set (filled in by OBMol::GetSSSR())
   */
  class OBAPI OBRingData : public OBGenericData
  {
    protected:
      std::vector<OBRing*> m_vr;
    public:
      /**
       * Constructor.
       */
      OBRingData();
      /**
       * Copy constructor.
       */
      OBRingData(const OBRingData &);
      /**
       * Clone this data.
       */
      virtual OBGenericData* Clone(OBBase* /*parent*/) const
      { return new OBRingData(*this); }
      /**
       * Destructor.
       */
      ~OBRingData();
      /**
       * Assign operator.
       */
      OBRingData &operator=(const OBRingData &);
      /**
       * Set the ring list.
       */
      void SetData(std::vector<OBRing*> &vr) { m_vr = vr; }
      /**
       * Add a ring to the list.
       */
      void PushBack(OBRing *r) { m_vr.push_back(r); }
      /**
       * Get the ring list.
       */
      std::vector<OBRing*>& GetData() { return m_vr; }
      /**
       * Get the begin iterator.
       */
      std::vector<OBRing*>::iterator BeginRings() { return m_vr.begin(); }
      /**
       * Get the end iterator.
       */
      std::vector<OBRing*>::iterator EndRings() { return m_vr.end(); }
      /**
       * Get the first ring.
       */
      OBRing *BeginRing(std::vector<OBRing*>::iterator &i);
      /**
       * Get the next ring.
       */
      OBRing *NextRing(std::vector<OBRing*>::iterator &i);
  };

  /** @class OBUnitCell generic.h <openbabel/generic.h>
      @brief Used for storing information about periodic boundary conditions
      with conversion to/from translation vectors and
      (a, b, c, alpha, beta, gamma)
   */
  class OBAPI OBUnitCell: public OBGenericData
  {
    public:
      enum LatticeType { 
        Undefined, 
        Triclinic, 
        Monoclinic, 
        Orthorhombic, 
        Tetragonal, 
        Rhombohedral /**< also called trigonal*/, 
        Hexagonal, 
        Cubic
      };
    protected:
      double m_a, m_b, m_c, m_alpha, m_beta, m_gamma;
      Eigen::Vector3d m_offset; //!< offset for origin
      Eigen::Vector3d m_v1, m_v2, m_v3; //!< translation vectors
      std::string m_spaceGroupName;
      const SpaceGroup* m_spaceGroup;
      LatticeType m_lattice;
    public:
      /** 
       * Public constructor
       */
      OBUnitCell();
      /**
       * Copy constructor.
       */
      OBUnitCell(const OBUnitCell &);
      /**
       * Clone this data.
       */
      virtual OBGenericData* Clone(OBBase* /*parent*/) const
      { return new OBUnitCell(*this); }
      /**
       * Destructor.
       */
      ~OBUnitCell() {}
      /**
       * Assign operator.
       */
      OBUnitCell &operator=(const OBUnitCell &);
      /**
       * Sets the vectors and angles of the unitcell
       * @param a The length a
       * @param b The length b
       * @param c The length c
       * @param alpha The angle alpha
       * @param beta The angle beta
       * @param gamma The angle gamma
       */
      void SetData(const double a, const double b, const double c,
          const double alpha, const double beta, const double gamma)
      {   
        m_a = a; 
        m_b = b; 
        m_c = c;
        m_alpha = alpha; 
        m_beta = beta; 
        m_gamma = gamma; 
      }
      /**
       * Sets the vectors (and angles) of the unitcell.
       */ 
      void SetData(const Eigen::Vector3d v1, const Eigen::Vector3d v2, const Eigen::Vector3d v3);
      /** 
       * Set the offset to the origin to @p v1.
       */
      void SetOffset(const Eigen::Vector3d v1) { m_offset = v1; }
      /** 
       * Set the space group for this unit cell.
       * Does not create an OBSymmetryData entry.
       */
      void SetSpaceGroup(const SpaceGroup* sg) { m_spaceGroup = sg; }
      /** 
       * Set the space group symbol for this unit cell.
       * Does not create an OBSymmetryData entry or attempt to convert
       * between different symbol notations
       */
      void SetSpaceGroup(const std::string sg) 
      { 
        m_spaceGroup = SpaceGroup::GetSpaceGroup (sg); 
        m_spaceGroupName = sg; 
      }
      /** 
       * Set the space group for this unit cell. Each spacegroup-symbol
       * has a numeric equivalent which is easier to use to convert between
       * notations.
       * Does not create an OBSymmetryData entry or attempt to convert
       * between different symbol notations
       */
      void SetSpaceGroup(const int sg) { m_spaceGroup = SpaceGroup::GetSpaceGroup (sg); }
      /** 
       * Set the Bravais lattice type for this unit cell.
       */
      void SetLatticeType(const LatticeType lt) { m_lattice = lt; }
      /** 
       * Duplicate symmetry-unique atoms to fill out the unit cell
       * of the molecule, based on the known space group.
       */
      void FillUnitCell(OBMol *);
      /**
       * @return vector a.
       */
      double GetA() { return m_a; }
      /** 
       * @return vector b.
       */
      double GetB() { return m_b; }
      /** 
       * @return vector c.
       */
      double GetC() { return m_c; }
      /** 
       * @return angle alpha.
       */
      double GetAlpha() { return m_alpha; }
      /** 
       * @return angle beta.
       */
      double GetBeta() { return m_beta; }
      /** 
     * @return angle gamma.
     */
      double GetGamma(){ return m_gamma; }
      /** 
       * @return any offset in the origin of the periodic boundaries.
       */
      Eigen::Vector3d GetOffset() { return m_offset; }
      /**
       * @return the text representation of the space group for this unit cell.
       */
      const SpaceGroup* GetSpaceGroup() { return m_spaceGroup; }
      /**
       * @return the text representation of the space group for this unit cell.
       */
      const std::string GetSpaceGroupName() { return m_spaceGroupName; }
      /** 
       * @return lattice type (based on the @p spacegroup).
       */
      LatticeType GetLatticeType( int spacegroup );
      /**
       * @return lattice type (based on angles and distances).
       */
      LatticeType GetLatticeType();
      /**
       * @return v1, v2, v3 cell vectors.
       */
      std::vector<Eigen::Vector3d> GetCellVectors();
      /**
       * @return v1, v2, v3 cell vectors as a 3x3 matrix.
       */
      Eigen::Matrix3d GetCellMatrix();
      /**
       * @return The orthogonalization matrix, used for converting from fractional to Cartesian coords.
       */
      Eigen::Matrix3d GetOrthoMatrix();
      /** 
       * @return The fractionalization matrix, used for converting from Cartesian to fractional coords.
       */
      Eigen::Matrix3d GetFractionalMatrix();
      /** 
       * @return The numeric value of the given spacegroup.
       */
      int GetSpaceGroupNumber( std::string name = "" );
      /** 
       * @return The cell volume (in Angstroms^3).
       */
      double GetCellVolume();
  };

  /** @class OBConformerData generic.h <openbabel/generic.h>
      @brief Used to hold data on conformers or geometry optimization steps
  
      Supplements the support for multiple coordinate sets in OBMol, e.g.,
      OBMol::GetConformer()

      @todo add example
   */
  class OBAPI OBConformerData: public OBGenericData
  {
    protected:
      //! Dimensionalities of conformers
      std::vector<unsigned short>                   m_vDimension;
      //! Relative energies of conformers (preferably in kJ/mol)
      std::vector<double>                           m_vEnergies;
      //! Atomic forces for each conformer
      std::vector< std::vector< Eigen::Vector3d > > m_vForces;
      //! Atomic velocities for each conformer (e.g., trajectories)
      std::vector< std::vector< Eigen::Vector3d > > m_vVelocity;
      //! Atomic displacements for each conformer (e.g., RMS distances)
      std::vector< std::vector< Eigen::Vector3d > > m_vDisplace;
      //! Additional data (as strings)
      std::vector<std::string>                      m_vData;
    
    public:
      /**
       * Constructor.
       */
      OBConformerData();
      /**
       * Destructor.
       */
      OBConformerData(const OBConformerData &);
      /**
       * Clone this data.
       */
      virtual OBGenericData* Clone(OBBase* /*parent*/) const
      { return new OBConformerData(*this); }
      /**
       * Destructor.
       */
      ~OBConformerData() {}
      /**
       * Assign operator.
       */
      OBConformerData &operator=(const OBConformerData &);
      /**
       * Set the dimensionalities of conformers.
       */
      void SetDimension(std::vector<unsigned short> vd) { m_vDimension = vd; }
      /**
       * Set the relative energies of conformers (preferably in kJ/mol).
       */
      void SetEnergies(std::vector<double> ve) { m_vEnergies = ve; }
      /** 
       * Set the atomic forces for each conformer.
       */
      void SetForces(std::vector< std::vector< Eigen::Vector3d > > vf) { m_vForces = vf; }
      /**
       * Set the atomic velocities for each conformer (e.g., trajectories).
       */
      void SetVelocities(std::vector< std::vector< Eigen::Vector3d > > vv) { m_vVelocity = vv; }
      /** 
       * Set the atomic displacements for each conformer (e.g., RMS distances).
       */
      void SetDisplacements(std::vector< std::vector< Eigen::Vector3d > > vd) { m_vDisplace = vd; }
      /**
       * Set additional data for each conformer (as strings).
       */
      void SetData(std::vector<std::string> vdat) { m_vData = vdat; }
      /**
       * Get the dimensionalities of conformers.
       */
      std::vector<unsigned short> GetDimension() { return m_vDimension; }
      /**
       * Get the relative energies of conformers (preferably in kJ/mol).
       */
      std::vector<double> GetEnergies() { return m_vEnergies; }
      /** 
       * Get the atomic forces for each conformer.
       */
      std::vector<std::vector< Eigen::Vector3d > > GetForces() { return m_vForces; }
      /**
       * Get the atomic velocities for each conformer (e.g., trajectories).
       */
      std::vector< std::vector< Eigen::Vector3d > > GetVelocities(){return m_vVelocity;}
      /** 
       * Get the atomic displacements for each conformer (e.g., RMS distances).
       */
      std::vector< std::vector< Eigen::Vector3d > > GetDisplacements(){return m_vDisplace;}
      /**
       * Get additional data for each conformer (as strings).
       */
      std::vector<std::string> GetData() { return m_vData; }
  };

  /** @class OBSymmetryData generic.h <openbabel/generic.h>
      @brief Used to hold the point-group and/or space-group symmetry
      @todo Add support for translation between symbol notations.
      @todo Add symmetry perception routines.
   */
  class OBAPI OBSymmetryData: public OBGenericData
  {
    protected:
      std::string m_pointGroup; //!< the point-group
      std::string m_spaceGroup; //!< the space-group
    public:
      /**
       * Constructor.
       */
      OBSymmetryData();
      /**
       * Copy constructor.
       */
      OBSymmetryData(const OBSymmetryData &);
      /**
       * Clone this data.
       */
      virtual OBGenericData* Clone(OBBase* /*parent*/) const
      { return new OBSymmetryData(*this); }
      /**
       * Destructor.
       */
      ~OBSymmetryData() {}
      /**
       * Assign operator.
       */
      OBSymmetryData &operator=(const OBSymmetryData &);
      /**
       * Set the point-group and space-group.
       */
      void SetData(std::string pg, std::string sg = "")
      { m_pointGroup = pg; m_spaceGroup = sg; }
      /**
       * Set the point-group.
       */
      void SetPointGroup(std::string pg) { m_pointGroup = pg; }
      /**
       * Set the space group.
       */
      void SetSpaceGroup(std::string sg) { m_spaceGroup = sg; }
      /**
       * Get the point-group.
       */
      std::string GetPointGroup() { return m_pointGroup; }
      /**
       * Get the space group.
       */
      std::string GetSpaceGroup() { return m_spaceGroup; }
  };

  //@}

  /** @class OBTorsion generic.h <openbabel/generic.h>
      @brief Used to hold the torsion data for a single rotatable bond
      and all four atoms around it
   */
  class OBAPI OBTorsion
  {
    friend class OBMol;
    friend class OBTorsionData;

    protected:
      //! atoms b-c in torsion a-b-c-d
      std::pair<OBAtom*,OBAtom*> m_bc;
      //! atoms a & d in torsion a-b-c-d and angle in radians
      std::vector<triple<OBAtom*,OBAtom*,double> > m_ads;
      /**
       * Constructor. (protected for use only by friend classes)
       */
      OBTorsion(): m_bc(NULL, NULL) {}
      /**
       * Constructor. (protected for use only by friend classes)
       */
      OBTorsion(OBAtom *, OBAtom *, OBAtom *, OBAtom *);
      /**
       * Get all the torsions around b-c.
       */  
      std::vector<quad<OBAtom*,OBAtom*,OBAtom*,OBAtom*> > GetTorsions();
    public:
      /**
       * Copy constructor.
       */
      OBTorsion(const OBTorsion &);
      /**
       * Destructor.
       */
      ~OBTorsion() {}
      /**
       * Assign operator.
       */
      OBTorsion& operator=(const OBTorsion &);
      /**
       * Clear all torsions.
       */
      /**
       * Return the OBTorsion to its original state.
       */
      void Clear();
      /**
       * @return True if bc is not set.
       */
      bool Empty() { return (m_bc.first == 0 && m_bc.second == 0); }
      /**
       * Add a new torsion to the OBTorsion object.
       */
      bool AddTorsion(OBAtom *a,OBAtom *b, OBAtom *c,OBAtom *d);
      /**
       * Add a new torsion to the OBTorsion object.
       */
      bool AddTorsion(quad<OBAtom*,OBAtom*,OBAtom*,OBAtom*> &atoms);
      /**
       * Set the angle of a torsion in OBTorsion.
       * @param radians the value to assign to the torsion.
       * @param index the index into the torsion of the OBTorsion.
       * @return boolean success.
       */
      bool SetAngle(double radians, unsigned int index = 0);
      /// @todo remove
      //bool SetData(OBBond * /*bond*/) { return false; }
      /**
       * Obtain the angle of a torsion in OBTorsion.
       * @param radians the value of the angle is set here.
       * @param index the index into the torsion of the OBTorsion.
       * @return boolean success.
       */
      bool GetAngle(double &radians, unsigned int index = 0);
      /** 
       * Get the bond index of the central bond.
       * @return The bond index
       */
      unsigned int GetBondIdx();
      /**
       * @return The number of torsions (a-b-c-d combinations, or # of ads).
       */
      unsigned int GetSize() const { return m_ads.size(); }
      /**
       * Get the two central atoms of ABCD torsion.
       * @return pair<OBAtom*,OBAtom*>
       */
      std::pair<OBAtom*,OBAtom*> GetBC() { return m_bc; }
      /**
       * Get the vector of distal atoms of ABCD torsion.
       * @return vector of A,D atom pointers and a double.
       */
      std::vector<triple<OBAtom*,OBAtom*,double> > GetADs() { return m_ads; }
      /**
       * Determine if torsion has only protons on either the a or d end.
       * @return boolean 
       */
      bool IsProtonRotor();
  };

  /// @addtogroup data Data classes 
  //@{

  /** @class OBTorsionData generic.h <openbabel/generic.h>
      @brief Used to hold torsions as generic data for OBMol.

      Filled by OBMol::FindTorsions()
   */
  class OBAPI OBTorsionData : public OBGenericData
  {
    friend class OBMol;

    protected:
      std::vector<OBTorsion> m_torsions; //!< the torsions

      /**
       * Constructor.
       */
      OBTorsionData();
      /**
       * Copy constructor.
       */
      OBTorsionData(const OBTorsionData &);
    
    public:
      /**
       * Assign operator.
       */
      OBTorsionData &operator=(const OBTorsionData &);
      /**
       * Clone this data.
       * @todo Needs to be updated to rebase atom pointers (or use indexes)
       */
      virtual OBGenericData* Clone(OBBase* /*parent*/) const
      { return new OBTorsionData(*this); }
      /**
       * Clear all torsions.
       */
      void Clear();
      /**
       * Get a vector of the OBTorsion objects.
       * @return the vector of torsions.
       */
      std::vector<OBTorsion> GetData() const { return m_torsions; }
      /**
       * Get the number of torsion structs.
       * @return integer count of the number of torsions.
       */
      unsigned int GetSize() const { return m_torsions.size(); }
      /**
       * Set the torsions.
       */
      void SetData(OBTorsion &torsion);
      /**
       * @brief Fills a vector with the indices of the atoms in torsions (ordered abcd).
       * @param torsions reference to the vector of abcd atom sets.
       * @return boolean success
       */
      bool FillTorsionArray(std::vector<std::vector<unsigned int> > &torsions);
  };

  //@}

  /** @class OBAngle generic.h <openbabel/generic.h>
      @brief Used to hold the 3 atoms in an angle and the angle itself.
   */
  class OBAPI OBAngle
  {
    friend class OBMol;
    friend class OBAngleData;

    protected:
      OBAtom                     *m_vertex;  //!< the vertex (atom b in angle a-b-c)
      std::pair<OBAtom*,OBAtom*>  m_termini; //!< atoms a & c in angle a-b-c
      double                      m_radians; //!< angle 

      /**
       * Constructor. Protected constructor for use only by friend classes.
       */
      OBAngle();
      /**
       * Constructor. Protected constructor for use only by friend classes.
       */
      OBAngle(OBAtom *vertex,OBAtom *a,OBAtom *b);
      /**
       * Get the angle atoms.
       */
      /**
       * Retrieve the 3 atom pointer for the angle (vertex first).
       * @return triple of OBAtom pointers.
       */
      triple<OBAtom*,OBAtom*,OBAtom*> GetAtoms();
      /**
       * Sort atoms in angle by order of indices.
       */
      void SortByIndex();
    public:
      /**
       * Copy constructor.
       */
      OBAngle(const OBAngle &);
      /**
       * Destructor.
       */
      ~OBAngle() { m_vertex = NULL; }
      /**
       * Assign operator.
       */
      OBAngle &operator=(const OBAngle &);
      /**
       * OBAngle equality operator, is same angle, NOT same value.
       * @return boolean equality
       */
      bool operator==(const OBAngle &);
      /**
       * Return OBAngle to its original state.
       */
      void Clear();
      /**
       * Get the OBAngle angle value.
       * @return angle in radians.
       */
      double GetAngle() const { return m_radians; }
      /**
       * Set the OBAngle to @p radians.
       * @param angle in radians
       */
      void SetAngle(double angle) { m_radians = angle; }
      /**
       * Set the 3 atoms in the angle.
       * Parameters are pointers to each OBAtom.
       */
      void  SetAtoms(OBAtom *vertex,OBAtom *a,OBAtom *b);
      /**
       * Set the 3 atoms in the angle.
       * @param atoms a triple of OBAtom pointers, the first must be the vertex.
       */
      void  SetAtoms(triple<OBAtom*,OBAtom*,OBAtom*> &atoms);
  };

  /// @addtogroup data Data classes 
  //@{

  //! \class OBAngleData generic.h <openbabel/generic.h>
  //! \brief Used to hold all angles in a molecule as generic data for OBMol
  class OBAPI OBAngleData : public OBGenericData
  {
    friend class OBMol;

    protected:
      std::vector<OBAngle> m_angles; //!< the angles

      /**
       * OBAngleData constructor.
       */
      OBAngleData();
      /**
       * OBAngleData copy constructor.
       */
      OBAngleData(const OBAngleData &);
      /**
       * Get the angle vector data.
       * @return a vector<OBAngle>
       */
      std::vector<OBAngle> GetData() const { return m_angles; }
      
    public:
      /**
       * OBAngleData assignment operator.
       */
      OBAngleData &operator =(const OBAngleData &);
      /**
       * Copy this data.
       */
      virtual OBGenericData* Clone(OBBase* /*parent*/) const
      { return new OBAngleData(*this); }
      /**
       * Set OBAngleData to its original state.
       */
      void Clear();
      /**
       * Fill an array with the indices of the atoms in the angle (vertex first).
       * @param angles pointer to the pointer to an array of angles atom indices.
       * @param size the current number of rows in the array.
       * @return int The number of angles.
       */
      unsigned int FillAngleArray(int **angles, unsigned int &size);
      /**
       * Fill an array with the indices of the atoms in the angle (vertex first).
       * @param angles pointer to the pointer to an array of angles atom indices.
       * @return True if successful.
       */
      bool FillAngleArray(std::vector<std::vector<unsigned int> > &angles);
      /**
       * Add a new angle to OBAngleData.
       */
      void         SetData(OBAngle &);
      /**
       * Get the number of angles stored.
       * @return integer count of the number of angles.
       */
      unsigned int GetSize() const { return m_angles.size(); }
  };

  //@}

  enum atomreftype{
    output,     //!< 
    input,      //!<
    calcvolume  //!<
  }; // sets which atom4ref is accessed by OBChiralData

  /// @addtogroup data Data classes 
  //@{

  /** @class OBChiralData generic.h <openbabel/generic.h>
      @brief Used to hold chiral inforamtion about the atom as OBGenericData.
       
      Chiral data is a perceived data type. We might read in some chiral info
      but this class derives and converts from whatever is read.
   */
  class OBAPI OBChiralData : public OBGenericData
  {
    friend class OBMol;
    friend class OBAtom;

    protected:
      std::vector<unsigned int> m_atom4refs; //!< input atom references
      std::vector<unsigned int> m_atom4refo; //!< output atom references
      std::vector<unsigned int> m_atom4refc; //!< calcvolume references

      //! The parity of the vector (of length 4)
      //! 1234 returns 0, 1243 returns 1
      int m_parity;

    public:
      /**
       * Constructor.
       */
      OBChiralData();
      /**
       * Copy constructor.
       */
      OBChiralData(const OBChiralData &src);
      /**
       * Clone this data.
       */
      virtual OBGenericData* Clone(OBBase* /*parent*/) const
      { return new OBChiralData(*this); }
      /**
       * Assign operator.
       */
      OBChiralData &operator =(const OBChiralData &);
      /**
       * Destructor.
       */
      ~OBChiralData() {}
      /**
       * Return OBChiralData data to its original state.
       */
      void Clear();
      /**
       * @return A vector of all 4 atom references of type @p t
       */
      std::vector<unsigned int> GetAtom4Refs(atomreftype t) const;
      /**
       * @return the atom reference specified by @p a of type @p t.
       */
      unsigned int GetAtomRef(int a,atomreftype t);
      /**
       * Set the 4 refs for type @p t.
       */
      bool SetAtom4Refs(std::vector<unsigned int> atom4refs, atomreftype t);
      /**
       * Add an atomref for type @p t. Call this functions for times to add 
       * all four refernece atoms.
       */
      int AddAtomRef(unsigned int atomref, atomreftype t);
      /** 
       * @return the size of the atom references of a given type @p t.
       */
      unsigned int GetSize(atomreftype t) const;
  };

  /** @class OBSerialNums generic.h <openbabel/generic.h>
      @brief Defines a map between serial numbers (e.g., in a PDB file) and 
      OBAtom objects inside a molecule.
   */
  class OBSerialNums : public OBGenericData
  {
    protected:
      std::map<int, OBAtom*> m_serialMap; //!< map between serial num
    public:
      /**
       * Constructor.
       */
      OBSerialNums() :
          OBGenericData("obSerialNums", OBGenericDataType::SerialNums)
      {  }
      /**
       * Copy constructor.
       */
      OBSerialNums(const OBSerialNums &cp) : OBGenericData(cp)
      {
        m_serialMap = cp.m_serialMap;
      }
      /** 
       * Clone this data.
       *
       * Member variables contain OBAtom pointers, so copying only valid within
       * same molecule, unless the code is modified, as in OBRotamerList.
       */
      virtual OBGenericData* Clone(OBBase* /*parent*/) const
      { return new OBSerialNums(*this); }
      /**
       * Get the serial number map.
       */
      std::map<int,OBAtom*> &GetData() { return m_serialMap; }
      /**
       * Set the serial number map.
       */ 
      void SetData(std::map<int,OBAtom*> &sm) { m_serialMap = sm; }
  };

  /** @class OBVibrationData generic.h <openbabel/generic.h>
      @brief Used to hold the normal modes of a molecule, etc.
   */
  class OBAPI OBVibrationData: public OBGenericData
  {
    protected:
      //! Normal modes in 1/sqrt(a.u.)
      std::vector< std::vector< Eigen::Vector3d > > m_vLx;
      //! Harmonic frequencies in inverse centimeters
      std::vector<double> m_vFrequencies;
      //! Infrared absorption intensities in KM/Mole
      std::vector<double> m_vIntensities;

    public:
      /**
       * Constructor.
       */
      OBVibrationData(): OBGenericData("VibrationData", OBGenericDataType::VibrationData){};
      /**
       * Destructor.
       */
      virtual ~OBVibrationData() {}
      /**
       * Clone this data.
       */
      virtual OBGenericData* Clone(OBBase*) const
      { return new OBVibrationData(*this); }
      /**
       * Assignment operator
       */
      OBVibrationData & operator=(const OBVibrationData &);
      /**
       * Assign the data.
       * @param vLx Normal modes in 1/sqrt(a.u.).
       * @param vFrequencies Harmonic frequencies in inverse centimeters.
       * @param vIntensities Infrared absorption intensities in KM/Mole.
       */
      void SetData(const std::vector< std::vector< Eigen::Vector3d > > &,
                   const std::vector<double> &,
                   const std::vector<double> &);
      /** 
       * @return Normal modes in 1/sqrt(a.u.).
       */
      std::vector< std::vector< Eigen::Vector3d > > GetLx() const { return this->m_vLx; }
      /**
       * @return Harmonic frequencies in inverse centimeters.
       */
      std::vector<double> GetFrequencies() const { return this->m_vFrequencies; }
      /**
       * @param vIntensities Infrared absorption intensities in KM/Mole.
       */
      std::vector<double> GetIntensities() const { return this->m_vIntensities; }
      /**
       * Get the number of frequencies.
       */
      unsigned int GetNumberOfFrequencies() const;
  };

  /** @class OBRotationData generic.h <openbabel/generic.h>
      @brief Used to hold the rotational constants and symmetry numbers.
   */
  class OBAPI OBRotationData: public OBGenericData
  {
    public:
     enum RType{
       UNKNOWN, 
       ASYMMETRIC, 
       SYMMETRIC, 
       LINEAR
     };
     /**
      * Constructor.
      */
     OBRotationData() : OBGenericData("RotationData", OBGenericDataType::RotationData) {}
     /**
      * Destructor.
      */
     virtual ~OBRotationData() {}
     /**
      * Clone this data.
      */
     virtual OBGenericData* Clone(OBBase*) const
     { return new OBRotationData(*this); }
     /**
      * Set the data.
      * @param type The type (linear, symmetric or asymmetric top)
      * @param RotationalConstants Rotational constants in GHz.
      * @param SymNum Rotational Symmetry Number.
      */
     void SetData(RType RotorType, std::vector<double> RotationalConstants, int SymmetryNumber)
     {
       m_RotConsts = RotationalConstants;
       m_type = RotorType;
       m_SymNum = SymmetryNumber;
     }
     /**
      * @return The rotational constants in GHz.
      */
     std::vector<double> GetRotConsts()const { return m_RotConsts; }
     /**
      * @return The rotational Symmetry Number.
      */
     int GetSymmetryNumber()const { return m_SymNum; }
     /**
      * @return The type (linear, symmetric or asymmetric top)
      */
     RType GetRotorType()const   { return m_type; }

   protected:
     std::vector<double> m_RotConsts;//!< Rotational constants in GHz
     int                 m_SymNum;   //!< Rotational Symmetry Number
     RType               m_type;     //!< linear, symmetric or asymmetric top
 };
 
  /** @class OBVectorData generic.h <openbabel/generic.h>
      @brief Used to hold a 3D vector item (e.g., a dipole moment)
      @since version 2.2
   */
  class OBAPI OBVectorData: public OBGenericData
  {
    public:
      /**
       * Constructor.
       */
      OBVectorData() : OBGenericData("VectorData", OBGenericDataType::VectorData) {}
      /**
       * Destructor.
       */
      virtual ~OBVectorData() {}
      /**
       * Clone this data.
       */
      virtual OBGenericData* Clone(OBBase*) const
      { return new OBVectorData(*this); }
      /**
       * Set the vector.
       */
      void SetData(double x, double y, double z) { m_vec = Eigen::Vector3d(x, y, z); }
      /**
       * Set the vector.
       */
      void SetData(Eigen::Vector3d data) { m_vec = data; }
      /**
       * Get the vector.
       */
      Eigen::Vector3d GetData() const { return m_vec; }
    protected:
      Eigen::Vector3d m_vec; //!< 3D vector to be stored
  };
 
  /** @class OBMatrixData generic.h <openbabel/generic.h>
       @brief Used to hold a 3x3 matrix item (e.g., a quadrupole moment)
       @since version 2.2
   */
  class OBAPI OBMatrixData: public OBGenericData
  {
    public:
      /**
       * Constructor.
       */
      OBMatrixData() : OBGenericData("MatrixData", OBGenericDataType::MatrixData) {}
      /**
       * Destructor
       */
      virtual ~OBMatrixData(){};
      /**
       * Clone this data
       */
      virtual OBGenericData* Clone(OBBase*) const
      { return new OBMatrixData(*this); }
      /**
       * Set the matrix.
       */
      void SetData(Eigen::Matrix3d data) { m_matrix = data; }
      /**
       * Get the matrix
       */
      Eigen::Matrix3d GetData() const { return m_matrix; }
    protected:
      Eigen::Matrix3d m_matrix; //!< 3x3 matrix to be stored
  };

  //@}

 //! A standard iterator over vectors of OBGenericData (e.g., inherited from OBBase)
  typedef std::vector<OBGenericData*>::iterator OBDataIterator;

} //end namespace OpenBabel

#endif // OB_GENERIC_H

//! \file generic.h
//! \brief Handle generic data classes. Custom data for atoms, bonds, etc.

/**********************************************************************
base.h - Base class for OpenBabel objects

Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2006 by Geoffrey R. Hutchison

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

#ifndef OB_BASE_H
#define OB_BASE_H

#include <openbabel/babelconfig.h>

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <openbabel/tokenst.h>

#ifdef UNUSED
#elif (__GNUC__ == 4)
# define UNUSED(x) UNUSED_ ## x __attribute__((unused))
#elif defined(__LCLINT__)
# define UNUSED(x) /*@unused@*/ x
#else
# define UNUSED(x) x
#endif

namespace OpenBabel
{

  //Forward declaration of the base class for OBMol OBReaction, OBAtom, etc.
  //Declaration later in this file.
class OBBase;
class OBConversion; //used only as pointer

//! \return the version of the Open Babel library for feature-detection (e.g. "2.3.1")
  OBAPI std::string OBReleaseVersion();

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
    enum
    {
      //! Unknown data type (default)
      UndefinedData =      0,

      //! Arbitrary key/value data, i.e., OBPairData
      PairData      =      1,

      //! Energetics data (e.g., total energy, heat of formation, etc.)
      EnergyData    =      2,

      //! Storing text comments (one per molecule, atom, bond, etc.) (for other data, e.g., author, keyword, ... use OBPairData)
      CommentData   =      3,

      //! Arbitrary information about conformers, i.e., OBConformerData
      ConformerData =      4,

      //! Bond data external to OpenBabel, i.e., OBExternalBond, OBExternalBondData
      ExternalBondData =   5,

      //! Information for generating & manipulating rotamers, i.e. OBRotamerList
      RotamerList =        6,

      //! Info. for storing bonds to atoms yet to be added, i.e. OBVirtualBond
      VirtualBondData =    7,

      //! Information on rings in a molecule, i.e., OBRingData
      RingData =           8,

      //! Information about torsion/dihedral angles, i.e., OBTorsionData and OBTorsion
      TorsionData =        9,

      //! Bond angles in a molecule, i.e., OBAngle, OBAngleData
      AngleData =         10,

      //! Residue serial numbers
      SerialNums =        11,

      //! Crystallographic unit cell data, i.e., OBUnitCell
      UnitCell =          12,

      //! Spin data, including NMR, atomic and molecular spin, etc.
      SpinData =          13,

      //! Arbitrary partial and total charges, dipole moments, etc.
      ChargeData =        14,

      //! Symmetry data -- point and space groups, transforms, etc. i.e., OBSymmetryData
      SymmetryData =      15,

      // 16 - Value unused, formerly ChiralData

      //! Atomic and molecular occupation data (e.g., for crystal structures)
      OccupationData =    17,

      //! Density (cube) data and surfaces
      DensityData =       18,

      //! Electronic levels, redox states, orbitals, etc.
      ElectronicData =    19,

      //! Vibrational modes, frequencies, etc.
      VibrationData =     20,

      //! Rotational energy information
      RotationData =      21,

      //! Nuclear transitions (e.g., decay, capture, fission, fusion)
      NuclearData =       22,

      //! Set Data (a set of OBGenericData)
      SetData =           23,

      //! Grid Data (e.g., 3D grids of data a.k.a. cubes)
      GridData =          24,

      //! Vector Data (i.e., one vector like a dipole moment)
      VectorData =        25,

      //! Matrix data (i.e., a 3x3 matrix like a rotation or quadrupole moment)
      MatrixData =        26,

      //! Stereochemistry data (see OBStereoBase)
      StereoData =        27,

      //! Density of States data (fermi energy and energy vs. density data)
      DOSData =           28,

      //! Electronic transition data (e.g., UV/Vis, excitation energies, etc.)
      ElectronicTransitionData = 29,

      // space for up to 2^14 more entries...

      //! Custom (user-defined data)
      CustomData0 = 16384,
      CustomData1 = 16385,
      CustomData2 = 16386,
      CustomData3 = 16387,
      CustomData4 = 16388,
      CustomData5 = 16389,
      CustomData6 = 16390,
      CustomData7 = 16391,
      CustomData8 = 16392,
      CustomData9 = 16393,
      CustomData10 = 16394,
      CustomData11 = 16395,
      CustomData12 = 16396,
      CustomData13 = 16397,
      CustomData14 = 16398,
      CustomData15 = 16399
    };
  } // end namespace
  enum DataOrigin {
    any,                 //!< Undefined or unspecified (default)
    fileformatInput,     //!< Read from an input file
    userInput,           //!< Added by the user
    perceived,           //!< Perceived by Open Babel library methods
    external,            //!< Added by an external program
    local                //!< Not for routine external use (e.g. in sdf or cml properties)
  };

  //! \brief Base class for generic data
  // Class introduction in generic.cpp
  // This base class declaration  has no dependence on mol.h
  class OBAPI OBGenericData
  {
  protected:
    std::string  _attr;  //!< attribute tag (e.g., "UnitCell", "Comment" or "Author")
    unsigned int _type;  //!< attribute type -- declared for each subclass
    DataOrigin   _source;//!< source of data for accounting
  public:
    OBGenericData(const std::string attr = "undefined",
                  const unsigned int type =  OBGenericDataType::UndefinedData,
                  const DataOrigin source = any);
    //Use default copy constructor and assignment operators
    //OBGenericData(const OBGenericData&);

    /* Virtual constructors added. see
       http://www.parashift.com/c++-faq-lite/abcs.html#faq-22.5
       to allow copying given only a base class OBGenericData pointer.
       It may be necessary to cast the return pointer to the derived class
       type, since we are doing without Covariant Return Types
       http://www.parashift.com/c++-faq-lite/virtual-functions.html#faq-20.8

       A derived class may return NULL if copying is inappropriate */
    virtual OBGenericData* Clone(OBBase* /*parent*/) const
    { return NULL; }
    virtual ~OBGenericData()    {}
    //Use default copy constructor and assignment operators
    //OBGenericData& operator=(const OBGenericData &src);

    //! Set the attribute (key), which can be used to retrieve this data
    void                      SetAttribute(const std::string &v)
    {        _attr = v;        }
    //! Set the origin of this data, which can be used to filter the data
    void SetOrigin(const DataOrigin s) { _source = s; }
    //! \return The attribute (key), which can be used to retrieve this data
    virtual const std::string &GetAttribute()  const
    {        return(_attr);    }
    //! \return the data type for this object as defined in OBGenericDataType
    unsigned int                GetDataType()    const
    {        return(_type);    }
    //! \brief Base class returns a default value (the attribute type)
    //! but should never be called
    virtual const std::string &GetValue()  const
    {			return _attr; }
    virtual DataOrigin GetOrigin() const
    {     return _source; }
  };

  //! A standard iterator over vectors of OBGenericData (e.g., inherited from OBBase)
  typedef std::vector<OBGenericData*>::iterator OBDataIterator;

  //! Base Class
  // introduction in base.cpp
  class OBAPI OBBase
    {
    public:
      virtual ~OBBase()
        {
          if (!_vdata.empty())
            {
              std::vector<OBGenericData*>::iterator m;
              for (m = _vdata.begin();m != _vdata.end();m++)
                delete *m;
              _vdata.clear();
            }
        }

      //! \brief Clear any and all data associated with this object
      virtual bool Clear();

      //! Perform a set of transformations specified by the user
      //!
      //! Typically these are program options to filter or modify an object
      //! For example, see OBMol::DoTransformations() and OBMol::ClassDescription()
      //! Base type does nothing
      virtual OBBase* DoTransformations(const std::map<std::string,std::string>* /*pOptions*/,
                                        OBConversion* /*pConv*/)
        {
          return this;
        }

      //! \return A list of descriptions of command-line options for DoTransformations()
      static const char* ClassDescription()
        {
          return "";
        }

      //! \brief By default clears the object. Called from ReadMolecule of most format classes
      template< class T >
      T* CastAndClear(bool clear=true)
        {
          T* pOb = dynamic_cast<T*>(this);
          if(pOb && clear)// Clear only if this is of target class
            Clear();
          return pOb;
        }

      //! \brief  Base type does nothing
      //! Made virtual around r3535 to simplify code which passes around OBBase*.
      //Currently no title data member in base class.
      virtual const char  *GetTitle(bool UNUSED(replaceNewlines) = true) const { return "";}
      virtual void  SetTitle(const char *) {}

      //! \name Generic data handling methods (via OBGenericData)
      //@{
      //! \return whether the generic attribute/value pair exists
      bool                              HasData(const std::string &);
      //! \return whether the generic attribute/value pair exists
      bool                              HasData(const char *);
      //! \return whether the generic attribute/value pair exists, for a given OBGenericDataType
      bool                              HasData(const unsigned int type);
      //! Delete any data matching the given OBGenericDataType
      void                              DeleteData(unsigned int type);
      //! Delete the given generic data from this object
      void                              DeleteData(OBGenericData*);
      //! Delete all of the given generic data from this object
      void                              DeleteData(std::vector<OBGenericData*>&);
      //! Deletes the generic data with the specified attribute, returning false if not found
      bool                              DeleteData(const std::string& s);
      //! Adds a data object; does nothing if d==NULL
      void                              SetData(OBGenericData *d)
        {
          if(d) _vdata.push_back(d);
        }
      //! Adds a copy of a data object; does nothing if d == NULL
      //! \since version 2.2
      void                              CloneData(OBGenericData *d);
      //! \return the number of OBGenericData items attached to this molecule.
      size_t                      DataSize() const
        { return(_vdata.size()); }
      //! \return the first matching data for a given type from OBGenericDataType
      //!    or NULL if nothing matches
      OBGenericData                    *GetData(const unsigned int type);
      //! \return any data matching the given attribute name or NULL if nothing matches
      OBGenericData                    *GetData(const std::string&);
      //! \return any data matching the given attribute name or NULL if nothing matches
      OBGenericData                    *GetData(const char *);
      //! \return the all matching data for a given type from OBGenericDataType
      //!    or an empty vector if nothing matches
      //! \since version 2.2
      std::vector<OBGenericData*>       GetAllData(const unsigned int type);
      //! \return all data, suitable for iterating
      std::vector<OBGenericData*>      &GetData() { return(_vdata); }
      //! \return all data with a specific origin, suitable for iterating
      std::vector<OBGenericData*>      GetData(DataOrigin source);
      //! \return An iterator pointing to the beginning of the data
      OBDataIterator  BeginData()
        { return(_vdata.begin());        }
      //! \return An iterator pointing to the end of the data
      OBDataIterator  EndData()
        { return(_vdata.end());          }
      //@}
    protected:
      std::vector<OBGenericData*> _vdata; //!< Custom data

    };

} //namespace OpenBabel

#endif // OB_BASE_H

//! \file base.h
//! \brief Base classes to build a graph

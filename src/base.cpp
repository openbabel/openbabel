/**********************************************************************
base.cpp - Base classes to build a graph

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

#include <string.h>

#include <openbabel/babelconfig.h>
#include <openbabel/base.h>

using namespace std;

//! Global namespace for all Open Babel code
namespace OpenBabel
{

  std::string OBReleaseVersion()
  {
    return std::string(BABEL_VERSION); // defined in babelconfig.h
  }

  /** \class OBBase base.h <openbabel/base.h>

  The various classes (Atom, Bond, Molecule) inherit from base classes--
  OBBase is largely a placeholder class. It also allows adding, deleting, and
  retrieving OBGenericData objects, which are ways to store arbitrary data
  for any atom, bond, molecule, or residue.

  For example, a graphics program may want to allow users to add labels to
  individual atoms:

  \code
  string atomLabel; // e.g., from the user adding annotation to an atom
  if (!atom.HasData("UserLabel")) // stored textual data as an OBPairData
  {
     OBPairData *label = new OBPairData;
     label->SetAttribute("UserLabel");
     label->SetValue(atomLabel);
     label->SetOrigin(userInput); // set by user, not by Open Babel

     atom.SetData(label);
  }
  \endcode

  This class is also important in the OBConversion class. Any derived class
  of OBBase can be supported in reading or writing data. While most OBFormat
  "translators" are designed around reading molecular data, the OBConversion
  framework can support any base object. For example OBReaction supports
  reading and writing reaction files, OBGrid supports reading and writing
  2D or 3D "grids" of numeric data.

  Therefore if you want to expand the range of input or output via the
  OBConversion and OBFormat classes, you will also need to make sure you define
  an appropriate derived class from OBBase.
  */

  //!
  //! This method can be called by OBConversion::Read() before reading data.
  //! Derived classes should be sure to call OBBase::Clear() to remove
  //! inherited generic data.
  //!
  //! \return Whether the call was successful.
  //! \since version 2.1.
  bool OBBase::Clear()
  {
    if (!_vdata.empty()) //clean up generic data
      {
        OBDataIterator m;
        for (m = _vdata.begin();m != _vdata.end();++m)
          delete *m;
        _vdata.clear();
      }

    return(true);
  }

  bool OBBase::HasData(const string &s)
    //returns true if the generic attribute/value pair exists
  {
    if (_vdata.empty())
      return(false);

    OBDataIterator i;

    for (i = _vdata.begin();i != _vdata.end();++i)
      if ((*i)->GetAttribute() == s)
        return(true);

    return(false);
  }

  bool OBBase::HasData(const char *s)
  {
    string temp(s);
    return(HasData(temp));
  }


  bool OBBase::HasData(const unsigned int dt)
    //returns true if the generic attribute/value pair exists
  {
    if (_vdata.empty())
      return(false);

    OBDataIterator i;

    for (i = _vdata.begin();i != _vdata.end();++i)
      if ((*i)->GetDataType() == dt)
        return(true);

    return(false);
  }

  //! \return the value given an attribute name
  OBGenericData *OBBase::GetData(const string &s)
  {
    OBDataIterator i;

    for (i = _vdata.begin();i != _vdata.end();++i)
      if ((*i)->GetAttribute() == s)
        return *i;

    return (OBGenericData*)0;
  }

  //! \return the value given an attribute name
  OBGenericData *OBBase::GetData(const char *s)
  {
    OBDataIterator i;

    for (i = _vdata.begin(); i != _vdata.end(); ++i)
      if (strcmp((*i)->GetAttribute().c_str(), s)==0)
        return *i;

    return (OBGenericData*)0;
  }

  OBGenericData *OBBase::GetData(const unsigned int dt)
  {
    OBDataIterator i;
    for (i = _vdata.begin();i != _vdata.end();++i)
      if ((*i)->GetDataType() == dt)
        return(*i);
    return(NULL);
  }

  std::vector<OBGenericData *> OBBase::GetAllData(const unsigned int dt)
  {
    std::vector<OBGenericData *> matches;

    // return all values matching this type
    OBDataIterator i;
    for (i = _vdata.begin();i != _vdata.end();++i)
      if ((*i)->GetDataType() == dt)
        matches.push_back(*i);

    return(matches);
  }

  std::vector<OBGenericData*>  OBBase::GetData(DataOrigin source)
  {
    std::vector<OBGenericData*> filtered; // filtered data only from source

    OBDataIterator i;
    for (i = _vdata.begin();i != _vdata.end();++i)
      if ((*i)->GetOrigin() == source)
        filtered.push_back((*i));

    return filtered;
  }

  void OBBase::CloneData(OBGenericData *d)
  {
    if(!d)
      return; // Nothing to do for NULL

    // Clone the data, relative to ourselves
    // This creates a new copy -- useable by scripting languages
    OBGenericData *clone = d->Clone(this);
    if (clone)
      _vdata.push_back(clone);

    return;
  }

  void OBBase::DeleteData(unsigned int dt)
  {
    vector<OBGenericData*> vdata;
    OBDataIterator i;
    for (i = _vdata.begin();i != _vdata.end();++i)
      if ((*i)->GetDataType() == dt)
        delete *i;
      else
        vdata.push_back(*i);
    _vdata = vdata;
  }

  void OBBase::DeleteData(vector<OBGenericData*> &vg)
  {
    vector<OBGenericData*> vdata;
    OBDataIterator i,j;

    bool del;
    for (i = _vdata.begin();i != _vdata.end();++i)
      {
        del = false;
        for (j = vg.begin();j != vg.end();++j)
          if (*i == *j)
            {
              del = true;
              break;
            }
        if (del)
          delete *i;
        else
          vdata.push_back(*i);
      }
    _vdata = vdata;
  }

  void OBBase::DeleteData(OBGenericData *gd)
  {
    OBDataIterator i;
    for (i = _vdata.begin();i != _vdata.end();++i)
      if (*i == gd)
        {
          delete *i;
          _vdata.erase(i);
          return; //Must stop since iterators invalidated by erase
        }
  }

  bool OBBase::DeleteData(const string& s)
  {
    OBDataIterator i;
    for (i = _vdata.begin();i != _vdata.end();++i)
    {
      if ((*i)->GetAttribute() == s)
      {
        delete *i;
          _vdata.erase(i);
          return true;
      }
    }
    return false;//not found
  }


  /// @addtogroup main Getting Started
  ///@{

  /*! \mainpage API Documentation

  \section base Introduction

  Open Babel is a full chemical software toolbox. In addition to converting
  file formats, it offers a complete programming library for developing
  chemistry software. The library is written primarily in C++ and also offers
  interfaces to other languages (e.g., Perl, Python, Ruby, and Java)
  using essentially the same API.

  This documentation outlines the Open Babel programming interface, providing
  information on all public classes, methods, and data. In particular, strives
  to provide as much (or as little) detail as needed. More information
  can also be found on the <a href="http://openbabel.org/">
  main website</a> and through the <a
  href="mailto:openbabel-discuss@lists.sourceforge.net">openbabel-discuss</a>
  mailing list.

  - \ref start "Getting Started" \n
     (where to begin, example code, using Open Babel in real life, ...)
  - \ref mainclasses "Classes Overview" \n
     (overview the most important classes ordered by category)
  - \ref changes23 "What's New in Version 2.3" \n
     (Changes since 2.2 releases)
  - \ref changes22 "What's New in Version 2.2" \n
     (Changes since 2.1 releases)
  - \ref changes21 "What's New in Version 2.1" \n
     (Changes since 2.0 releases)
  - <a href="annotated.shtml" class="el">All Classes</a> \n
     (all classes with brief descriptions)

  \page intro Introduction to Open Babel API

  Open Babel is a full chemical software toolbox. In addition to converting
  file formats, it offers a complete programming library for developing
  chemistry software. The library is written primarily in C++ and also offers
  interfaces to other languages (e.g., Perl and Python) using essentially
  the same API.

  The heart of Open Babel lies in the \link OpenBabel::OBMol OBMol\endlink,
  \link OpenBabel::OBAtom OBAtom\endlink, and
  \link OpenBabel::OBBond OBBond\endlink classes,
  which handle operations on atoms, bonds and molecules. Newcomers should
  start with looking at the \link OpenBabel::OBMol OBMol\endlink class,
  designed to store the basic information
  in a molecule and to perceive information about a molecule.

  \section Lazy Evaluation
  One of the key philosophies in the code is that transformations and
  automatic perception of properties are performed in a <a href="http://en.wikipedia.org/wiki/Lazy_evaluation">"lazy"</a>
  manner. That is, until you call for partial atomic charges, no
  charges are calculated. This ensures faster transformations of
  chemical data -- properties that are not needed for your code will
  typically not be calculated. When such data is needed, appropriate
  routines are called, and a "flag" is set (e.g., via OBMol::SetFlag
  or OBAtom::SetFlag etc.) so that the code is only run once.

  \section Arbitrary Data
  Arbitrary custom data and text descriptors can be stored in any atom,
  bond, molecule, or residue using the \link OpenBabel::OBGenericData
  OBGenericData\endlink or \link OpenBabel::OBPairData
  OBPairData\endlink classes.

  Conversion between various chemical file formats is accomplished through
  the \link OpenBabel::OBConversion OBConversion\endlink and \link
  OpenBabel::OBFormat OBFormat\endlink classes, often through use of the \link
  OpenBabel::OBMoleculeFormat OBMoleculeFormat\endlink subclass which is designed
  for easy read/write access to one or more \link OpenBabel::OBMol OBMol\endlink
  objects. The philosophy of the file format codes is to parse as much
  chemical information from a given file as possible (no data left
  behind) and ideally any perception or transformations will occur when
  writing to some other format later.

  \page start Getting Started

  Not surprisingly, the Open Babel library is a full chemical
  toolbox. So to start out, the first example is to read in molecular
  file data and uses the \link OpenBabel::OBMol OBMol\endlink,
  and \link OpenBabel::OBConversion OBConversion\endlink classes. The
  \link OpenBabel::OBMol former\endlink is designed to store the basic
  information in a molecule and to perceive information and chemical
  properties. The \link OpenBabel::OBConversion latter\endlink is
  designed to handle conversion of a variety of data (i.e., not just
  molecules) and import and export.

  This example program shows how to read in a molecule, check the
  number of atoms, and write a SMILES string.
  \code
      #include <iostream.h>

      // Include Open Babel classes for OBMol and OBConversion
      #include <openbabel/mol.h>
      #include <openbabel/obconversion.h>

      int main(int argc,char **argv)
      {
         // Read from STDIN (cin) and Write to STDOUT (cout)
         OBConversion conv(&cin,&cout);

         // Try to set input format to MDL SD file
         // and output to SMILES
         if(conv.SetInAndOutFormats("SDF","SMI"))
         {
            OBMol mol;
            if(conv.Read(&mol))
            {
               //  ...manipulate molecule
               cerr << " Molecule has: " << mol.NumAtoms()
                    << " atoms." << endl;
            }

            // Write SMILES to the standard output
            conv->Write(&mol);
         }
         return 0; // exit with success
      }
  \endcode

  All of the \ref mainclasses "main classes", including OBMol and
  OBConversion, include example code designed to facilitate using the
  Open Babel code in real-world chemistry. Complete examples are provided
  in the doc/examples directory and <a href="examples.shtml">all these examples are included in
  the documentation</a>.

  For a further list of example code, see the
  <a
  href="http://openbabel.org/wiki/Developer:Tutorial">developer
  tutorials</a>. This section includes examples in C++, Perl, Python,
  and other programming languages.

  Also, the <code>tools</code> directory of the Open Babel source
  releases include a variety of programs which are intended to be more
  advanced examples (although usually still under 300 lines of code).

  Please e-mail the openbabel-discuss@lists.sourceforge.net mailing list
  if you have more questions!

  \page other Further Information

  Open Babel is a community project. In addition to this API documentation,
  the website offers a variety of up-to-date and useful information for
  developing with the library.

  Open Babel homepage:
  - <a href="http://openbabel.org/wiki/Develop">Developing with Open Babel</a>
  - <a href="http://openbabel.org/wiki/Developer:Tutorial">Developer
  tutorials</a>
  - <a href="http://openbabel.org/wiki/Contribute">Contributing
  to the Open Babel project</a>

  SourceForge project pages:
  - <a href="http://sourceforge.net/projects/openbabel/">Open Babel project page</a>
  - <a href="http://sourceforge.net/tracker/?group_id=40728&atid=428740">Bug reporter</a>
  - <a href="http://sourceforge.net/tracker/?group_id=40728&atid=428743">Feature requests</a>
  - <a href="http://sourceforge.net/tracker/?group_id=40728&atid=447448">File format requests</a>
  - <a href="http://sourceforge.net/tracker/?group_id=40728&atid=428743">Mailing lists and archives</a>

  */

  ///@}

  /**
  \page mainclasses Classes Overview

  Chemical Storage and Manipulation:
  - OBBase - General base class of OBAtom, OBMol, OBResidue, etc. for storing generic, arbitrary custom data
    - OBMol - Central molecule class, properties, ...
    - OBAtom - Central atom class, properties, coordinates, ...
    - OBBond - Bond connection between two OBAtom, properties, bond orders, ...
    - OBResidue - Biomolecule residues (amino acids, nucleic acids)
    - OBReaction - Used to store chemical reactions (i.e., reactants -> products)
    - OBGrid - A base grid class (e.g., 3D functions like electrostatic potential, evaluated on a grid)
    - OBText - An object containing just text
  - OBRing - Ring cycle perception, Smallest Set of Smallest Rings (SSSR), Largest Set of Smallest Rings (LSSR)

  Import / Export:
  - OBConversion - Conversion between file formats. Interface for import or export
  - OBFormat - Parent class for file format types for any sort of data

  Chemical Utilities:
  - OBFingerprint - Generating molecular fingerprints for databases and similarity searches
  - OBForceField - Molecular mechanics force fields for energy evaluation, conformational searching, coordinate generation, etc.
  - OBChemTsfm - SMARTS based structural modification (chemical transformation)

  Substructure Searching:
  - OBSmartsPattern - Parsing SMARTS chemical search patterns and matching against OBMol objects
  - OBQuery
  - OBQueryAtom
  - OBQueryBond
  - OBIsomorphismMapper

  Utility Classes:
  - OBBitVec - Efficient bit vector (e.g., for marking visit to atoms)
  - OBMessageHandler - Error and warning handling, filtering, and logging
  - OBRandom - Random number generator
  - OBStopwatch - Stopwatch class used for timing length of execution
  - OBSqrtTbl - Square Root lookup table - given a distance squared returns distance

  Coordinate Manipulation:
  - OBBuilder - 3D structure generation
  - OBAlign - Kabsch alignment
  - OBInternalCoord - Internal coordinates
  - OBConformerSearch - Conformer searching using genetic algorithm
  - OBRotorList - Class for detecting rotatable bonds
  - OBRotamerList - Supports a set of rotamer coordinate sets for some number of potentially rotatable bonds

  Stereochemistry:
  - OBStereo
  - OBStereoBase
    - OBTetraNonPlanarStereo
      - OBTetrahedralStereo
    - OBTetraPlanarStereo
      - OBCisTransStereo
      - OBSquarePlanarStereo
  - OBStereoFacade

  Arbitrary Data:
  - OBGenericData - General base class for arbitrary data types
    - OBPairData - Arbitrary text key/value data
    - OBPairTemplate - Arbitrary T key/value data
    - OBUnitCell - Storage and manipulation of crystal structure / reciprocal unit cells
    - AliasData, OBAngleData, OBCommentData, OBConformerData, OBDOSData, OBElectronicTransitionData,
      OBExternalBondData, OBGridData, OBMatrixData, OBNasaThermoData, OBOrbitalEnergyData, OBPairData, OBRateData,
      OBRingData, OBRotamerList, OBRotationData, OBSerialNums, OBSetData, OBStereoBase, OBSymmetryData, OBTorsionData, OBUnitCell,
      OBVectorData, OBVibrationData, OBVirtualBond

  Plugins:
  - OBPlugin
    - OBChargeModel - assign partial (and formal) charges to a molecule
    - OBDescriptor - molecular properties, descriptors or features for filtering or searching
    - OBFingerprint - condensed representations of molecules as a set of bits
    - OBForceField - implementation of molecular mechanics force field methods
    - OBFormat - file formats for chemical data
    - OBOp - operation for modifying a molecule or file (e.g., adding hydrogens)

  Math Utilities:
  - matrix3x3 - Square 3x3 matrices for 3D transformations and rotations
  - vector3 - 3D vector class for translations and planes
  */

  /**
  \page changes23 What's New in Version 2.3

  Throughout the API documentation, new classes and methods are
  indicated with a disclaimer "Since: version 2.3."

  Several new classes and features have been introduced, notably
  - New @ref stereo classes
  - New @ref substructure classes
  - Improved @ref conformer searching
  - OBChargeModel - Cleaner interface for assigning partial charges using different models
  - OBGraphSym - Cleaner interface for computing the stereochemical symmetry classes
  - OBSpectrophore - A new powerful descriptor for screening bioactive molecules
  - OBDepict - 2D depiction

  In addition, a wide variety of new plugins are available, including
  
  - computation of partial charges using GasteigerCharges, MMFF94Charges, QEqCharges, and QTPIECharges
  - "diverse" conformer searching
  - etc.

  \page changes22 What's New in Version 2.2

  Throughout the API documentation, new classes and methods are
  indicated with a disclaimer "Since: version 2.2."

  Several new classes have been introduced, notably
  - OBPlugin - Cleaner interface to dynamic plugins, including OBFormat
  OBForceField, OBFingerprint, and the new OBDescriptor and OBOp.
  - SpaceGroup - Handle crystallographic space group symmetry.
  - AliasData - Handle atoms which serve as "aliases" for larger fragments
  for example in drawing programs - "COOH" vs. four atoms for the carboxylic
  acid. Currently used for MDL 2D files.
  - OBGridData - Store and process numeric grid data involved with orbitals or
  electrostatic potential.
  - OBBuilder - New class to handle generation of 3D coordinates for SMILES,
  InChI and other formats.
  - OBForceField - Revised, improved class, including use of constraints.
  - OBMoleculeFormat - Installed by default, to facilitate new format development.
  - OBLocale - New class to ensure chemical data is read in the "C" numeric locale.
  - OBPhModel - Now accepts pH-dependent transformations.

  In addition, a wide variety of plugins are available, including MMFF94
  and UFF force fields, SMARTS-based descriptors for filtering data sets,
  compound descriptors such as the Lipinsky Rule of Five, and chemical
  transformations.

  Optional headers include integration with the Boost Graph Library.

  \page changes21 What's New in Version 2.1

  Throughout the API documentation, new classes and methods are
  indicated with a disclaimer "Since: version 2.1." Documentation, in
  general should be more thorough and complete.

  One major change is that the library uses the same header
  <code>\#include</code> paths as external code. This documentation
  attempts to uniformly use the correct style, e.g.:
  <code>\#include \<openbabel/mol.h\></code>

  Another major change was to increase the indexing of atoms and
  bonds. Previously only 2^16 atoms or bonds were possible.

  The OBNodeBase, OBEdgeBase, and OBGraphBase classes were
  deprecated. Typedef declarations provided backwards compatibility by
  mapping to OBAtom, OBBond, and OBMol, respectively.

  In addition, this page gives a general list of additions to the library.

  - obiter.h - Additional iterators, including breadth-first and depth-first
      iterators for atoms in OBMol, convenience iterators for rings, angles,
      torsions, and 1-4 non-bonded atom pair interactions.
  - OBGenericData - Generic data is now inherited through OBBase for
      most types, including OBAtom, OBBond, OBMol, OBResidue,
      OBReaction, and OBGrid.
  - OBForceField - Generic interface to molecular mechanics force
      fields, including Ghemical (Tripos-like) methods, and ultimately
      others such as MM2, MMFF94, etc.
  - pluginiter.h - Generic interface to "plugin" classes, including
      force fields, and fingerprints.
  - vector3 - General cleanup.
  - matrix3x3 - General cleanup.
  - OBMol - Copy and assignment now copy associated OBGenericData objects.
  - OBResidueData - Now available for public use, rather than left as private
      to PDB format.
  - OBSetData - Hierarchical sets of generic data.
  - LineEndingExtractor - New class to correct DOS and Classic Mac line endings
     to UNIX-style newline characters.
  */

} // namespace OpenBabel

//! \file base.cpp
//! \brief Implementation of base classes.

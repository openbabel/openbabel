/**********************************************************************
data.h - Global data and resource file parsers.

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

#ifndef OB_DATA_H
#define OB_DATA_H

#include <openbabel/babelconfig.h>

#include <stdio.h>
#include <cstring>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>

namespace OpenBabel
{

  class OBAtom;
  class OBMol;

  /** \class OBGlobalDataBase data.h <openbabel/data.h>
      \brief Base data table class, handles reading data files

      Base data table class--reads ASCII data files in various formats
      -# Checks for the environment variable _envvar (defaults to "BABEL_DATADIR")
      - Tries the _subdir directory if defined (def. "data") and then the main directory
      -# Checks for the directory _dir (def. determined by the build environment)
      - Tries the subdirectory corresponding to this version, then the main directory
      -# Reverts to the compiled-in default data
  **/
  class OBAPI OBGlobalDataBase
    {
    protected:
      bool         _init;		//!< Whether the data been read already
      const char  *_dataptr;//!< Default data table if file is unreadable
      std::string  _filename;//!< File to search for
      std::string  _dir;		//!< Data directory for file if _envvar fails
      std::string  _subdir;	//!< Subdirectory (if using environment variable)
      std::string  _envvar;	//!< Environment variable to check first

    public:
      //! Constructor
      OBGlobalDataBase(): _init(false), _dataptr(NULL) { }
      //! Destructor
      virtual ~OBGlobalDataBase()                  {}
      //! Read in the data file, falling back as needed
      void  Init();
      //! \return the size of the database (for error checking)
      virtual size_t GetSize()                 { return 0;}
      //! Set the directory before calling Init()
      void  SetReadDirectory(char *dir)            { _dir = dir;    }
      //! Set the environment variable to use before calling Init()
      void  SetEnvironmentVariable(char *var)      { _envvar = var; }
      //! Specified by particular table classes (parses an individual data line)
      virtual void ParseLine(const char*)          {}
    };

  /** \class OBAtomHOF data.h <openbabel/data.h>
      \brief helper class for OBAtomicHeatOfFormationTable

      Stores both theoretical and experimental corrections
      needed to compute the Enthalpy of formation. In order to
      use these you need to perform
      Gaussian G2/G3/G4 or CBS-QB3 calculations.
  **/
  class OBAPI OBAtomHOF
  {
  private:
      std::string _element,_method,_desc,_unit;
    double _T,_value;
    int _charge;
    int _multiplicity;

  public:
    /** \brief Initialize Heat of Formation for atom
        
     @param element The element string
     @param charge  The formal charge of the particle (if an ion)
     @param method  The method used for determining the value
     @param desc    Description of the value
     @param T       Temperature
     @param value   The value of the property (energy)
     @param multiplicity The multiplicity of the atomic system
     @param unit    The (energy) unit
    */
    OBAtomHOF(std::string element,int charge,
              std::string method,std::string desc,
              double T,double value,int multiplicity,
              std::string unit)
      {
        _element      = element;
        _charge       = charge;
        _method       = method;
        _desc         = desc;
        _T            = T;
        _value        = value;
        _multiplicity = multiplicity;
        _unit         = unit;
      }

    /** \brief Destructor */
    ~OBAtomHOF() {}
    /** \brief Return the chemical element */
    std::string Element() { return _element; }
    /** \brief Return the formal charge */
    int Charge()          { return _charge; }
    /** \brief Return the method used for the measurement/calculation */
    std::string Method()  { return _method; }
    /** \brief Return specification of the measurement/calculation type */
    std::string Desc()    { return _desc; }
    /** \brief Return the temperature */
    double T()            { return _T; }
    /** \brief Return the (energy) value */
    double Value()        { return _value; }
    /** \brief Return the multiplicity */
    int Multiplicity()    { return _multiplicity; }
    /** \brief Return the (energy) unit */
    std::string Unit()    { return _unit; }
  };

  /** \class OBAtomicHeatOfFormationTable data.h <openbabel/data.h>
      \brief Atomic Heat of Formation Table

      Contributions of atoms to Enthalpy of Formation calculations performed
      in Gaussian, using the G2/G3/G4 or CBS-QB3 methods respectively.
      The energies produced by Gaussian have to be corrected according to their
      document on Thermochemistry with Gaussian. The data in the file
      BABEL_DATA/atomization_energies.txt supplies this information based on
      single atom calculations with Gaussian and the appropriate method and
      experimental data from Curtiss et al., J. Chem. Phys. 106 (1997) 1063-1079.
  */
  class OBAPI OBAtomicHeatOfFormationTable : public OBGlobalDataBase
  {
    std::vector<OBAtomHOF> _atomhof;

    public:
      /** \brief Constructor */
      OBAtomicHeatOfFormationTable(void);
      /** \brief Destructor */
      ~OBAtomicHeatOfFormationTable() {}

      //! \return the number of elements in the Atomic Heat Of Formation table
      size_t GetSize() { return _atomhof.size(); }

      /** \brief Read one line in the file and parse it 
          @param Unnamed the line to be parsed
      */
      void	ParseLine(const char*);
      /** \brief Extract heat of formation and entropy for an atom
       @param elem         The chemical element we're looking for
       @param charge       At this formal charge
       @param method       The method used for computing/measuring
       @param T            The temperature
       @param dhof0        The output energy at 0K
       @param dhof1        The output energy at T
       @param S0T          The entropy at T (it is 0 at 0K)
       \return 1 if the contribution to the Heat of Formation for this atom
       is known at temperature T. If 1 the values
       including all corrections are returned in the dhof variable.
      */
      int	GetHeatOfFormation(std::string elem, 
                               int charge,
                               std::string method,
                               double T, double *dhof0,
                               double *dhofT,double *S0T);
    };

  // class introduction in data.cpp
  class OBAPI OBTypeTable : public OBGlobalDataBase
    {
      int             _linecount;
      unsigned int    _ncols,_nrows;
      int             _from,_to;
      std::vector<std::string> _colnames;
      std::vector<std::vector<std::string> > _table;

    public:

      OBTypeTable(void);
      ~OBTypeTable() {}

      void ParseLine(const char*);

      //! \return the number of atom types in the translation table
      size_t GetSize() { return _table.size(); }

      //! Set the initial atom type to be translated
      bool SetFromType(const char*);
      //! Set the destination atom type for translation
      bool SetToType(const char*);
      //! Translate atom types
      bool Translate(char *to, const char *from); // to, from
      //! Translate atom types
      //! \return whether the translation was successful
      bool Translate(std::string &to, const std::string &from); // to, from
      //! Translate atom types
      //! \return the translated atom type, or an empty string if not possible
      std::string Translate(const std::string &from);

      //! \return the initial atom type to be translated
      std::string GetFromType();
      //! \return the destination atom type for translation
      std::string GetToType();
    };

  //! Global OBTypeTable for translating between different atom types
  //! (e.g., Sybyl <-> MM2)
  EXTERN  OBTypeTable      ttab;

  /** \class OBResidueData data.h <openbabel/data.h>
      \brief Table of common biomolecule residues (for PDB or other files).

      Can assign atom types and bond orders for arbitrary residues
  **/
  class OBAPI OBResidueData : public OBGlobalDataBase
    {
      int                                               _resnum;
      std::vector<std::string>                          _resname;
      std::vector<std::vector<std::string> >            _resatoms;
      std::vector<std::vector<std::pair<std::string,int> > > _resbonds;

      //variables used only temporarily for parsing resdata.txt
      std::vector<std::string>                          _vatmtmp;
      std::vector<std::pair<std::string,int> >          _vtmp;
    public:

      OBResidueData();
      void ParseLine(const char*);

      //! \return the number of residues in the table
      size_t GetSize() { return _resname.size(); }

      //! Sets the table to access the residue information for a specified
      //!  residue name
      //! \return whether this residue name is in the table
      bool SetResName(const std::string &);
      //! \return the bond order for the bond specified in the current residue
      //! \deprecated Easier to use the two-argument form
      int  LookupBO(const std::string &);
      //! \return the bond order for the bond specified between the two specified
      //! atom labels
      int  LookupBO(const std::string &, const std::string&);
      //! Look up the atom type and hybridization for the atom label specified
      //! in the first argument for the current residue
      //! \return whether the atom label specified is found in the current residue
      bool LookupType(const std::string &,std::string&,int&);
      //! Assign bond orders, atom types and residues for the supplied OBMol
      //! based on the residue information assigned to atoms
      bool AssignBonds(OBMol &);
    };

  //! Global OBResidueData biomolecule residue database
  EXTERN  OBResidueData    resdat;


} // end namespace OpenBabel
  
#endif //DATA_H

//! \file data.h
//! \brief Global data and resource file parsers.

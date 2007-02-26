/**********************************************************************
obmolecformat.h - Subclass of OBFormat for conversion of OBMol.

Copyright (C) 2005 Chris Morley

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

#ifndef OB_MOLECULEFORMAT_H
#define OB_MOLECULEFORMAT_H

#ifdef _WIN32
  #include <hash_map>
#endif

#include <openbabel/mol.h>
#include <openbabel/obconversion.h>

namespace OpenBabel {

// This macro is used in DLL builds. If it has not
// been set in babelconfig.h, define it as nothing.
#ifndef OBCOMMON
  #define OBCOMMON
#endif

/** \class OBMoleculeFormat obmolecformat.h <openbabel/obmolecformat.h>
    \brief An OBFormat convenience subclass for conversion to/from OBMol data

    This class is not intended for direct use outside of Open Babel, unless
    you're writing a new format converting to or from an OBMol molecule.
    (e.g., see http://openbabel.sourceforge.net/wiki/HowTo:Add_A_New_File_Format).

    An OBFormat which converts to and/or from OBMol can derive from this class
    to save duplicating the ReadChemObject() and/or WriteChemObject() methods.
    Derive directly from OBFormat if the object converted is not OBMol or 
    if interaction with the framework is required during the execution 
    of ReadMolecule() or WriteMolecule(), as for example in CMLFormat
**/
class OBCOMMON OBMoleculeFormat : public OBFormat
{
private:
  static std::map<std::string, OBMol*> IMols;
  static OBMol* _jmol; //!< Accumulates molecules with the -j option
  static std::vector<OBMol> MolArray; //!< Used in --separate option
  static bool StoredMolsReady; //!< Used in --separate option

public:

  OBMoleculeFormat()
  {
    OBConversion::RegisterOptionParam("b", this, 0, OBConversion::INOPTIONS);
    OBConversion::RegisterOptionParam("s", this, 0, OBConversion::INOPTIONS);
    OBConversion::RegisterOptionParam("title", this, 1,OBConversion::GENOPTIONS);
    OBConversion::RegisterOptionParam("addtotitle", this, 1,OBConversion::GENOPTIONS);
    OBConversion::RegisterOptionParam("property", this, 2, OBConversion::GENOPTIONS);
    OBConversion::RegisterOptionParam("C",        this, 0,OBConversion::GENOPTIONS);
    OBConversion::RegisterOptionParam("j",        this, 0,OBConversion::GENOPTIONS);
    OBConversion::RegisterOptionParam("join",     this, 0,OBConversion::GENOPTIONS);
    OBConversion::RegisterOptionParam("separate", this, 0,OBConversion::GENOPTIONS);

    //The follow are OBMol options, which should not be in OBConversion.
    //But here isn't entirely appropriate either, since one could have
    //OBMol formats loaded but which don't derived from this class.
    //However, this possibility is remote.
    OBConversion::RegisterOptionParam("s", NULL, 1,OBConversion::GENOPTIONS);
    OBConversion::RegisterOptionParam("v", NULL, 1,OBConversion::GENOPTIONS);
    OBConversion::RegisterOptionParam("h", NULL, 0,OBConversion::GENOPTIONS);
    OBConversion::RegisterOptionParam("d", NULL, 0,OBConversion::GENOPTIONS);
    OBConversion::RegisterOptionParam("b", NULL, 0,OBConversion::GENOPTIONS);
    OBConversion::RegisterOptionParam("c", NULL, 0,OBConversion::GENOPTIONS);
    OBConversion::RegisterOptionParam("p", NULL, 0,OBConversion::GENOPTIONS); 
    OBConversion::RegisterOptionParam("t", NULL, 0,OBConversion::GENOPTIONS);
    OBConversion::RegisterOptionParam("k", NULL, 0,OBConversion::GENOPTIONS);
  };

  //! Static routine,  which can be called from elsewhere
  static bool ReadChemObjectImpl(OBConversion* pConv, OBFormat*);
  //! Static routine,  which can be called from elsewhere
  static bool WriteChemObjectImpl(OBConversion* pConv, OBFormat*);

  /// The "Convert" interface for reading a new molecule
  virtual bool ReadChemObject(OBConversion* pConv)
  { return ReadChemObjectImpl(pConv, this);}
    
  /// The "Convert" interface for writing a new molecule
  virtual bool WriteChemObject(OBConversion* pConv)
  { return WriteChemObjectImpl(pConv, this);}
  
  /// \name Routines to handle the -C option for combining data from several OBMols
  //@{
  //! Defer output of a molecule until later, so it can be combined with others
  //! \return Success, or false if no molecule was read.
  static bool   DeferMolOutput(OBMol* pmol, OBConversion* pConv, OBFormat* pF);
  //! Write out all molecules queued with DeferMolOutput
  static bool   OutputDeferredMols(OBConversion* pConv);
  //! Delete the list of queued molecules from DeferMolOutput
  static bool   DeleteDeferredMols();
  //! \return the OBMol which combines @p pFirst and @p pSecond (i.e.)
  static OBMol* MakeCombinedMolecule(OBMol* pFirst, OBMol* pSecond);
  //@}

#ifdef _WIN32
  typedef stdext::hash_map<std::string, unsigned> NameIndexType;
#else
  typedef std::map<std::string, unsigned> NameIndexType;
#endif
  
  // documentation in obmolecformat.cpp
  static bool   ReadNameIndex(NameIndexType& index, const std::string& datafilename,
                  OBFormat* pInFormat);

  //! \return the type of data converted by this format (here, OBMol)
  const std::type_info& GetType()
  {
    return typeid(OBMol*);
  }
//////////////////////////////////////////////////////////////

};

}
#endif //OB_MOLECULEFORMAT_H

//! \file obmolecformat.h
//! \brief Subclass of OBFormat for conversion of OBMol.

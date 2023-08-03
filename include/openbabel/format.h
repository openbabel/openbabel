/**********************************************************************
obconversion.cpp -  Declarations for OBFormat

Copyright (C) 2004-2007 by Chris Morley
Some portions Copyright (C) 2005-2007 by Geoffrey Hutchison

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
#ifndef OB_FORMAT_H
#define OB_FORMAT_H
#include <openbabel/babelconfig.h>
#include <openbabel/plugin.h>
#include <typeinfo>

namespace OpenBabel
{
  class OBBase;
  class OBConversion;

  ///For OBFormat::Flags()
#define NOTREADABLE     0x01
#define READONEONLY     0x02
#define READBINARY      0x04
#define ZEROATOMSOK     0x08
#define NOTWRITABLE     0x10
#define WRITEONEONLY    0x20
#define WRITEBINARY     0x40
#define READXML         0x80
#define DEPICTION2D     0x100
#define DEFAULTFORMAT   0x4000

  /// @brief Base class for file formats.
  // class introduction in obconversion.cpp
class OBCONV OBFormat : public OBPlugin
  {
    //Macro to include functions to handle plugin operations
    MAKE_PLUGIN(OBFormat);

  public:

    ///Default constructor. Registration via RegisterFormat(), not via constructor as in other plugins.
    OBFormat(){}

    const char* TypeID(){ return "formats"; }

    /// @brief The "API" interface Read function.

    /// Reads a single object.
    /// Does not make a new object on the heap;
    /// can be used with a pointer to an chem object on the heap or the stack.
    virtual bool ReadMolecule(OBBase* /*pOb*/, OBConversion* /*pConv*/)
      { std::cerr << "HIER" << std::endl;
std::cerr << "Not a valid input format"; return false;}

    /// @brief The "Convert" interface Read function.

    /// Possibly reads multiple new objects on the heap and subjects them
    /// to its DoTransformations() function, which may delete them again.
    /// Sends result to pConv->AddChemObject()
    virtual bool ReadChemObject(OBConversion* /*pConv*/)
      { std::cerr << "Not a valid input format"; return false;}

    /// @brief The "API" interface Write function.

    /// Writes a single object
    /// Does not delete the object;
    /// can be used with a pointer to an chem object on the heap or the stack.
    /// \return false on error.
    virtual bool WriteMolecule(OBBase* /*pOb*/, OBConversion* /*pConv*/)
      { std::cerr << "Not a valid output format"; return false;}

    /// @brief The "Convert" interface Write function.

    /// Writes a single object
    /// Deletes the object after writing
    /// \return false on error
    virtual bool WriteChemObject(OBConversion* /*pConv*/)
      { std::cerr << "Not a valid output format"; return false;}

    /// @brief Information on this format. Printed out in response to -Hxxx option where xxx id the id of the format.

    /// Must be provided by each format class.
    /// Can include a list of command line Options. These may be used to construction
    /// check boxes, radio buttons etc for GUI interface.
    virtual const char* Description()=0;

    /// @brief A decription of the chemical object converted by this format.

    /// If not provided, the object type used by the default format is used (usually OBMol).
    virtual const char* TargetClassDescription();

    /// \return the type of chemical object used by the format.

    /// Defaults to that used by the default format. Useful for checking
    /// that a format can handle a particular object.
    virtual const std::type_info& GetType();

    /// @brief Web address where the format is defined.
    virtual const char* SpecificationURL() { return ""; }

    /// @brief Chemical MIME type associated with this file type (if any)
    virtual const char* GetMIMEType() { return pMime; }

    /// @brief Decribes the capabilities of the format (Read only etc.)

    /// Currently, can be a bitwise OR of any of the following
    /// NOTREADABLE READONEONLY NOTWRITABLE WRITEONEONLY DEFAULTFORMAT
    /// READBINARY WRITEBINARY READXML
    virtual unsigned int Flags() { return 0;};

    /// @brief Skip past first n objects in input stream (or current one with n=0)

    /// \return 1 on success, -1 on error and 0 if not implemented
    virtual int SkipObjects(int /*n*/, OBConversion* /*pConv*/)
      {
        return 0; //shows not implemented in the format class
      };

    /// \return a pointer to a new instance of the format, or NULL if fails.

    /// Normally a single global instance is used but this may cause problems
    /// if there are member variables and the format is used in more than one place
    /// in the program.
    virtual OBFormat* MakeNewInstance()
      {
        return nullptr; //shows not implemented in the format class
      }

    //New functions since OBFormat is derived from OBPlugin
    //\brief Called from, and an alternative to, OBConversion::RegisterFormat();
    int RegisterFormat(const char* ID, const char* MIME = nullptr);

    ///\brief Provides a description in txt of the format specified by itr.
    ///If param starts with "in", "read", "out" or "write" only the
    ///appropriate formats are output. The others return false.
    ///If param contains "verbose", the whole description is output.
    virtual bool Display(std::string& txt, const char* param, const char* ID=nullptr);

    static OBFormat* FormatFromMIME(const char* MIME);

private:
    static PluginMapType &FormatsMIMEMap()
    {
      static PluginMapType m;
      return m;
    }

    const char* pMime;
/* Functions provided by the MAKE_PLUGIN macro

  ///Constructor that registers the ID of the format
  Not currently used for formats
  OBFormat(const char* ID, bool IsDefault=false);

  ///Returns the sub-type associated with the ID, or the default subtype if ID NULL or empty.
  static OBFormat* FindType(const char* ID);

*/};

}//namespace
#endif

//! \file format.h
//! \brief Declarations for OBFormat

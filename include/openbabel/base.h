/**********************************************************************
base.h - Base classes to build a graph
 
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

#ifndef OB_BASE_H
#define OB_BASE_H

#include <openbabel/babelconfig.h>

#include <vector>
#include <map>

#include <iostream>

// Needed for the OBGenericData handling in this class
// (can't just forward declare it.)
#include <openbabel/generic.h>

namespace OpenBabel
{

  // Utility function prototypes
  OBAPI bool tokenize(std::vector<std::string>&, const char *buf, const char *delimstr=" \t\n\r");
  OBAPI bool tokenize(std::vector<std::string>&, std::string&, const char *delimstr=" \t\n\r", int limit=-1);
  // Remove leading and trailing whitespace from a string (docs in tokenst.cpp)
  OBAPI std::string& Trim(std::string& txt);

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
      virtual OBBase* DoTransformations(const std::map<std::string,std::string>* /*pOptions*/)
        {
          return this;
        } 

      //Base type does nothing
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
        };


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
      //! \return the number of OBGenericData items attached to this molecule.
      unsigned int                      DataSize() const 
        { return(_vdata.size()); }
      //! \return the first matching data for a given type from OBGenericDataType
      //!    or NULL if nothing matches
      OBGenericData                    *GetData(const unsigned int type);
      //! \return any data matching the given attribute name or NULL if nothing matches
      OBGenericData                    *GetData(const std::string&);
      //! \return any data matching the given attribute name or NULL if nothing matches
      OBGenericData                    *GetData(const char *);
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

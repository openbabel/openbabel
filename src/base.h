/**********************************************************************
base.h - Base classes to build a graph
 
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

#ifndef OB_BASE_H
#define OB_BASE_H

#include "babelconfig.h"

#include <vector>
#include <map>

#include <iostream>

#include "generic.h"

namespace OpenBabel
{

  class OBBase;

  /** \brief Base Class
 
  The various classes (Atom, Bond, Molecule) inherit from base classes--
  OBBase is just a placeholder class
  */
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

      virtual OBBase* DoTransformations(const std::map<std::string,std::string>* /*pOptions*/)
        {
          return this;
        } 
      //Base type does nothing
      static const char* ClassDescription()
        {
          return "";
        } 

      //! \name Generic data handling methods (via OBGenericData)
      //@{
      //! \returns whether the generic attribute/value pair exists
      bool                              HasData(const std::string &);
      //! \returns whether the generic attribute/value pair exists
      bool                              HasData(const char *);
      //! \returns whether the generic attribute/value pair exists, for a given
      //!  OBGenericDataType
      bool                              HasData(const unsigned int type);
      //! Delete any data matching the given OBGenericDataType
      void                              DeleteData(unsigned int type);
      //! delete the given generic data from this object
      void                              DeleteData(OBGenericData*);
      //! delete all of the given generic data from this object
      void                              DeleteData(std::vector<OBGenericData*>&);
      //! adds a data object; does nothing if d==NULL
      void                              SetData(OBGenericData *d)
        {
          if(d)
            _vdata.push_back(d);
        }
      //! \return the number of OBGenericData items attached to this molecule.
      unsigned int                      DataSize() const 
        { return(_vdata.size()); }
      //! \return the first matching data for a given type from OBGenericDataType
      //!    or NULL if nothing matches
      OBGenericData                    *GetData(const unsigned int type);
      //! \return any data matching the given attribute name 
      //!     or NULL if nothing matches
      OBGenericData                    *GetData(const std::string&);
      //! \return any data matching the given attribute name 
      //!     or NULL if nothing matches
      OBGenericData                    *GetData(const char *);
      //! \return all data, suitable for iterating
      std::vector<OBGenericData*>      &GetData() { return(_vdata); }
      std::vector<OBGenericData*>::iterator  BeginData()
        {
          return(_vdata.begin());
        }
      std::vector<OBGenericData*>::iterator  EndData()
        {
          return(_vdata.end());
        }
      //@}
    protected:
      std::vector<OBGenericData*> _vdata; //!< Custom data

    };

} //namespace OpenBabel

#endif // OB_BASE_H

//! \file base.h
//! \brief Base classes to build a graph

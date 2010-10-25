/**********************************************************************
Copyright (C) 2007 by Chris Morley

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#ifndef OB_ATOMCLASS_H
#define OB_ATOMCLASS_H

#include <openbabel/babelconfig.h>

#include <vector>
#include <sstream>
#include <openbabel/base.h>

namespace OpenBabel
{
// This macro is used in DLL builds. If it has not
// been set in babelconfig.h, define it as nothing.
#ifndef OBAPI
  #define OBAPI
#endif

  /** \class OBAtomClassData atomclass.h <openbabel/atomclass.h>
      \brief Handle atom classes in reaction SMILES/SMIRKS
      \since version 2.2

      Class for attaching to OBMol to hold the info from SMILES like [C:2]
      Useful for reaction SMILES (SMIRKS). It influences the atom id attribute in CML.
      Not all atoms need have an atom class.
      The atom class can be any positive or negative integer.
  */
class OBAPI OBAtomClassData : public OBGenericData
{
protected:
  std::map<int,int> _map; //index is atom index; value is class
public:
  OBAtomClassData(): OBGenericData("Atom Class", 0x7882){ }
  virtual OBGenericData* Clone(OBBase* /*parent*/) const{return new OBAtomClassData(*this);}

  ///Erase contents
  void Clear(){ _map.clear(); }

  ///Add an individual value
  void Add(int indx, int cl) { _map[indx] = cl;}

  /// \return true if there is an entry for atom index
  bool HasClass(int indx)const { return _map.find(indx)!=_map.end(); }
  /// \return value of class index (Test with HasClass first)
  int GetClass(int indx)const
  {
    std::map<int,int>::const_iterator pos = _map.find(indx);
    if(pos!=_map.end())
      return pos->second;
    return -9999;
  }

  /// If there is an entry for indx, return ":n" where n is the atomclass value;
  /// otherwise return an empty string
  std::string GetClassString(int indx)
  {
    std::stringstream ss;
    std::map<int,int>::const_iterator pos = _map.find(indx);
    if(pos!=_map.end())
      ss << ':' << pos->second;
    return ss.str();
  }
  int size(){ return _map.size(); }

};
} //namespace

#endif // OB_ATOMCLASS_H

//! \file atomclass.h
//! \brief Handle atom classes in reaction SMILES/SMIRKS

/**********************************************************************
alias.h - OBGenericData class to hold alias information on atoms
Copyright (C) Copyright (C) 2007 by Chris Morley

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#ifndef OB_ALIAS_H
#define OB_ALIAS_H

#include <vector>
#include <openbabel/generic.h>
#include <openbabel/shared_ptr.h>

namespace OpenBabel
{
  class OBSmartsPattern;
  class OBMol;

// This macro is used in DLL builds. If it has not
// been set in babelconfig.h, define it as nothing.
#ifndef OBCOMMON
  #define OBCOMMON
#endif

const unsigned int AliasDataType = 0x7883;

/** \class AliasData alias.h <openbabel/alias.h>
  \brief Indicate atoms as aliases for larger functional groups
  \since version 2.2

An object of this class can be attached to an OBAtom if it is considered
to be a placeholder for an alias, such as Ph, COOH, etc.

If the alias has not been interpreted chemically (expanded), the element type of the
placeholder atom should be set to Xx so that the molecule is not interpreted
incorrectly by formats which do not consider this class.

If the alias has been interpreted chemically (which formats, etc. should
normally do immediately), the alias may remain as extra
information or as a hint for an alternative representation, for example to
a chemical drawing program. The _expandedatoms vector would then contains
the ids of the atoms to which the alias is an alternative.
*/
class OBCOMMON AliasData : public OBGenericData
{
protected:
  std::string _alias;
  std::string _right_form;
  std::vector<unsigned long> _expandedatoms; //atom ids (not idxs)
  std::string _color;
public:

  AliasData(): OBGenericData("Alias", AliasDataType){ }

  virtual OBGenericData* Clone(OBBase* /*parent*/) const{return new AliasData(*this);}

  ///Add an alias
  void SetAlias(const std::string& alias) {_alias = alias;}
  void SetAlias(const char* alias) {_alias = alias;}

  /// /return value of alias or its version intended to be connected at its right hand end.
  std::string GetAlias(bool rightAligned = false) const;

  /// Return the color which has been assigned to this alias
  std::string GetColor() const { return _color; }

  /// Assign a color to this alias
  void SetColor(std::string color){ _color = color; }

  bool IsExpanded()const { return !_expandedatoms.empty(); }

  ///Converts all the expanded aliases in a molecule to the alias form.
  ///Note that this deletes atoms and bonds. Use only as a preparation for display.
  static void RevertToAliasForm(OBMol& mol);

  ///Interprets the alias text and adds atoms as appropriate to mol.
  bool Expand(OBMol& mol, const unsigned int atomindex);

 #ifdef HAVE_SHARED_POINTER
 ///Matches univalent fragments in the molecule and adds AliasData to provide an alternative representation.
  ///The molecule remains a respectable OBMol. Data in superatoms.txt decides what aliases are added.
  static bool AddAliases(OBMol* pmol);
 #endif

private:
  /// Interpret the alias as a formula
  bool FormulaParse(OBMol& mol, const unsigned atomindex);

  /// Add one of the atoms that can replace the alias
  void AddExpandedAtom(int id);

  /// Delete the atoms that replace the alias
  void DeleteExpandedAtoms(OBMol& mol);

  struct AliasItem
  {
    std::string right_form;
    std::string smiles;
    std::string color;
  };
  typedef std::map<std::string, AliasItem> SuperAtomTable; //key=alias left-form

  static bool LoadFile(SuperAtomTable& table);
  static SuperAtomTable& table()
  {
    static SuperAtomTable t;
    if(t.empty())
      LoadFile(t);
    return t;
  }
  bool        FromNameLookup(OBMol& mol, const unsigned int atomindex);
#ifdef HAVE_SHARED_POINTER
  typedef std::vector< std::pair<std::string, obsharedptr<OBSmartsPattern> > > SmartsTable;
  static bool LoadFile(SmartsTable& smtable);
#endif
};
} //namespace

#endif // OB_ALIAS_H

//! \file alias.h
//! \brief OBGenericData class to for atom alias data (e.g., in 2D drawing programs for "COOH")

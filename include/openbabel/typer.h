/**********************************************************************
typer.h - Open Babel atom and aromaticity typer.
 
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

#ifndef OB_TYPER_H
#define OB_TYPER_H

#include <openbabel/parsmart.h>
#include <openbabel/data.h>

namespace OpenBabel
{

  // class introduction in typer.cpp
  class OBAPI OBAtomTyper : public OBGlobalDataBase
  {
    protected:
      std::vector<std::vector<int> >                        m_mlist;   //!< match list for atom typing
      std::vector<std::pair<OBSmartsPattern*,int> >         m_vinthyb; //!< internal hybridization rules
      std::vector<std::pair<OBSmartsPattern*,int> >         m_vimpval; //!< internal implicit valence rules
      std::vector<std::pair<OBSmartsPattern*,std::string> > m_vexttyp; //!< external atom type rules

    public:
      /** 
       * @brief Constructor.
       */
      OBAtomTyper();
      /** 
       * @brief Destructor.
       */
      ~OBAtomTyper();
      /** 
       * @brief Parse a line from the database file.
       */
      void ParseLine(const char*);
      /** 
       * @return The number of implicit valence rules
       */
      unsigned int GetSize() { return m_vimpval.size(); }
      /** 
       * @brief Assign atomic hybridization (1 = sp, 2 = sp2, 3 = sp3...).
       * 
       * @param mol The molecule.
       */
      void AssignHyb(OBMol &mol);
      /** 
       * @brief Assign external atomic types (i.e., EXTTYP lines in atomtyp.txt).
       * 
       * @param mol The molecule.
       */
      void AssignTypes(OBMol &mol);
      /** 
       * @brief Assign implicit valence (i.e., given an atomic type, what is the 
       * expected number of bonds to this atom.
       * 
       * @param mol The molecule.
       */
      void AssignImplicitValence(OBMol &mol);
      /** 
       * @brief Correct typing, valence, and hybridization for aromatic nitrogen atoms.
       * 
       * @param mol The molecule.
       */
      void CorrectAromaticNitrogens(OBMol &mol);
  };

  // class introduction in typer.cpp
  class OBAPI OBAromaticTyper : public OBGlobalDataBase
  {
    protected:
      std::vector<bool>                 m_vpa;     //!< potentially aromatic atoms
      std::vector<bool>                 m_visit;   //!< list of visited atoms
      std::vector<bool>                 m_root;    //!< list of start atoms
      std::vector<std::vector<int> >    m_mlist;   //!< smarts match list
      std::vector<OBSmartsPattern*>     m_vsp;     //!< SMARTS of potentially aromatic atoms
      std::vector<std::pair<int,int> >  m_verange; //!< min and max number of electrons
      std::vector<std::pair<int,int> >  m_velec;   //!< # electrons an atom contributes
    public:
      /** 
       * @brief Constructor.
       */
      OBAromaticTyper();
      /** 
       * @brief Destructor.
       */
      ~OBAromaticTyper();
      /** 
       * @return The number of SMARTS patterns.
       */
      unsigned int GetSize() { return m_vsp.size(); }
      /** 
       * @brief Parse a line from the database file. 
       */
      void ParseLine(const char*);
      /** 
       * @brief Assign aromaticity flag to atoms and bonds.
       * 
       * @param The molecule.
       */
      void AssignAromaticFlags(OBMol &mol);
      /** 
       * @brief "Anti-alias" potentially aromatic flags around a molecule
       * (aromatic atoms need to have >= 2 neighboring ring atoms) 
       * 
       * @param atom The current atom.
       */
      void PropagatePotentialAromatic(OBAtom *atom);
      /** 
       * @brief Select the root atoms for traversing atoms in rings. 
       * 
       * @param mol The molecule.
       * @param avoidInnerRingAtoms Inner closure ring atoms with more than 
       * 2 neighbours will be avoided.
       *
       * Picking only the begin atom of a closure bond can cause
       * difficulties when the selected atom is an inner atom
       * with three neighbour ring atoms. Why ? Because this atom
       * can get trapped by the other atoms when determining aromaticity,
       * because a simple visited flag is used in the
       * OBAromaticTyper::TraverseCycle() method.
       *
       * Ported from JOELib, copyright Joerg Wegner, 2003 under the GPL version 2
       */
      void SelectRootAtoms(OBMol &mol, bool avoidInnerRingAtoms = true);
      /** 
       * @brief Remove 3-member rings from consideration.
       * 
       * @param mol The molecule
       */
      void ExcludeSmallRing(OBMol &mol);
      /** 
       * @brief Check aromaticity starting from the root atom, up to a specified depth.
       * 
       * @param root The start atom.
       * @param searchDepth The dpeth.
       */
      void CheckAromaticity(OBAtom *root, int searchDepth);
      /** 
       * @brief Traverse a potentially aromatic cycle starting at @p root. 
       * 
       * @param root  The initial, "root" atom in traversing this ring.
       * @param atom  The current atom to visit and check.
       * @param prev  The bond traversed in moving to this @p atom
       * @param er    The min and max number of pi electrons for this ring.
       * @param depth The maximum number of atoms to visit in a ring (e.g., 6).
       * 
       * @return True if the cycle is likely aromatic.
       *
       * This method traverses a potentially aromatic ring, adding up the possible
       * pi electrons for each atom. At the end (e.g., when @p atom == @p root)
       * the Huekel 4n+2 rule is checked to see if there is a possible electronic
       * configuration which corresponds to aromaticity.
       */
      bool TraverseCycle(OBAtom *root, OBAtom *atom, OBBond *prev,
          std::pair<int,int> &er,int depth);
  };

  // class introduction in typer.cpp
  class OBAPI OBRingTyper : public OBGlobalDataBase
  {
    protected:
      std::vector<std::vector<int> >                        m_mlist;   //!< match list for atom typing
      std::vector<std::pair<OBSmartsPattern*,std::string> > m_ringtyp; //!< ring type rules
    public:
      /** 
       * @brief Constructor.
       */
      OBRingTyper();
      /** 
       * @brief Destructor.
       */
      ~OBRingTyper();
      /** 
       * @brief Parse a line from the database file. 
       */
      void ParseLine(const char*);
      /** 
       * @return The number of SMARTS patterns.
       */
      unsigned int GetSize() { return m_ringtyp.size();}
      /** 
       * @brief Assign external atomic types (ringtyp.txt).
       * 
       * @param mol The molecule.
       */
      void AssignTypes(OBMol &mol);
  };

} //namespace OpenBabel

#endif // OB_TYPER_H

//! @file typer.h
//! @brief Open Babel atom and aromaticity typer.

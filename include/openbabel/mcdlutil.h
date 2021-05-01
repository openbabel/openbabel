/*-*-C++-*-

**********************************************************************
Copyright (C) 2007,2008 by Sergei V. Trepalin sergey_trepalin@chemical-block.com
Copyright (C) 2007,2008 by Andrei Gakh andrei.gakh@nnsa.doe.gov

This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************
*/
/*
  Diagram is generated using templates, which are stored in SD file templates.sdf
  The SD file is usual SD file, which contain chemical structures and might contain data.
  Only chemical structures are used. Subgraph isomorphisme search is executed and coordinates
  of atoms are determined from templates. See Molecules, 11, 129-141 (2006) for algorithm decription.
  Structures in SD file are converted in next manner:
  1. All atoms, except explicit hydrogens, are replaced with generic ANY_ATOM (matched with any atom in subgraph isomorphisme search)
  2. All bonds are replaces with generic ANY_BOND, which can be matched with any bond in molecule
  3. All hydrogen are removed, but they are used for search-query and structure atom matching is believed fo be
     successful if chemical structure contains more or equal number of hydrogens, than query. Using explicitly-defined hydrogens
	 on query enables ones to remove substitutors attachment for atom, which are sterically hidden on templates
  if the file will not be found, predefined templates will be used
*/


namespace OpenBabel {

  //common constants
  static const int MAXBONDS=300;
  static const int MAXFRAGS=200;
  static const int MAXCHARS=1000;
  static const int MAX_DEPTH=10;
  static const int NELEMMAX=120;
  #define NELEMMCDL 121

  // Return valency by hydrogen for given atomic position in the Periodic Table
  OBMCDL int hydrogenValency(int na);
  OBMCDL int maxValency(int na);

  //Alternate overloaded methods
  OBMCDL int alternate(OBMol * pmol, const int nH[], int bondOrders []);  //This method does not work!
  //Zero-based atomic numeration should be in connection matrix arrays iA1 and iA2-so first atom has indez zero
  OBMCDL int alternate(const std::vector<int> aPosition,const std::vector<int> aCharge,
      const std::vector<int> aRad,const std::vector<int> nHydr, const std::vector<int> iA1,
      const std::vector<int> iA2, std::vector<int> & bondOrders, int nAtoms, int nBonds);
  /**
   * @since version 2.3
   */
  //Diagram generation overloaded methods
  OBMCDL void generateDiagram(OBMol * pmol);
  //Zero-based atomic numeration should be in connection matrix arrays iA1 and iA2-so first atom has indez zero
  OBMCDL void generateDiagram(const std::vector<int> iA1, const std::vector<int> iA2,
      std::vector<double>& rx, std::vector<double>& ry, int nAtoms, int nBonds);
  OBMCDL void generateDiagram(OBMol * pmol, std::ostream & ofs);  //for testing purposes only

  //Fragment search - pure subgraph isomorphisme
  bool fragmentSearch(OBMol * query, OBMol * structure);
  bool fragmentSearch(const std::vector<int> aPositionQuery, const std::vector<int> iA1Query,
      const std::vector<int> iA2Query, const std::vector<int> bondTypesQuery,
      const std::vector<int> aPositionStructure, const std::vector<int> iA1Structure,
      const std::vector<int> iA2Structure,  const std::vector<int> bondTypesStructure,
      int nAtomsQuery, int nBondsQuery, int nAtomsStructure, int nBondsStructure);
  ///Equivalence list generation
  OBMCDL void equivalenceList(OBMol * pmol,  std::vector<int>& eqList);
  void equivalenceList(const std::vector<int> aPosition,const std::vector<int> aCharge,
      const std::vector<int> aRad, const std::vector<int> iA1, const std::vector<int> iA2,
      const std::vector<int> bondTypes,  std::vector<int>& eqList, int nAtoms, int nBonds);
  //Fragment addition
  OBMCDL void addFragment(OBMol * molecule, OBMol * fragment, int molAN, int fragAN, int molBN,
      int fragBN, bool isAddition);

  //routines below have no common meaning, but are necessary to process stereo information
  OBMCDL void createStereoLists(OBMol * pmol, std::vector<int>& bondStereoList,
      std::vector<int>& atomStereoList, std::vector<int>& eqList);
  OBMCDL std::string getAtomMCDL(OBMol * pmol, int ntatoms, const std::vector<int> ix,
      const std::vector<int> aNumber, const std::vector<int> atomStereoList, const std::vector<int> eqList);
  OBMCDL std::string getBondMCDL(OBMol * pmol, int nbStore, int ntatoms, const std::vector<int> ix,
      const std::vector<int> aNumber, int bonds[MAXBONDS][4], const std::vector<int> bondStereoList,
      const std::vector<int> eqList);
  OBMCDL void implementAtomStereo(std::vector<int>& iA1, std::vector<int>& iA2, std::vector<int>& stereoBonds,
      const std::vector<double>rx, const std::vector<double> ry, int acount, int bcount, std::string astereo);
  OBMCDL void implementBondStereo(const std::vector<int> iA1, const std::vector<int> iA2,
      std::vector<double>& rx, std::vector<double>& ry, int acount, int bcount, std::string bstereo);


  OBMCDL int groupRedraw(OBMol * pmol, int bondN, int atomN, bool atomNInGroup);
  //int  groupRedrawFrameAtom(OBMol * pmol, int bondN, int atomInFrame);

  OBMCDL int  canonizeMCDL(const std::string atomBlock, std::vector<std::string> & structureList);
  OBMCDL bool parseFormula(const std::string formulaString, std::vector <int>& enumber, int & valency);

  OBMCDL void prepareTest(OBMol * pmol, std::ostream & ofs);

} // namespace OpenBabel

/// @file mcdlutil.h
/// @brief 2D molecule coordinate generation.

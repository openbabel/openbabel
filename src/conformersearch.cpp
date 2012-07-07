/**********************************************************************
conformersearch.cpp - Conformer searching using genetic algorithm.

Copyright (C) 2010 Tim Vandermeersch

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

#include <openbabel/conformersearch.h>
#include <openbabel/math/align.h>
#include <openbabel/forcefield.h>

namespace OpenBabel {

  //////////////////////////////////////////////////////////
  //
  //  OBConformerFilter(s)
  //
  //////////////////////////////////////////////////////////

  bool OBStericConformerFilter::IsGood(const OBMol &mol, const RotorKey &key, double *conformer)
  {
    unsigned int numAtoms = mol.NumAtoms();
    for (unsigned int a1 = 0; a1 < numAtoms; ++a1) {
      for (unsigned int a2 = 0; a2 < numAtoms; ++a2) {
        // skip the pair if the atoms are the same
        // also, only check each pair once
        if (a1 <= a2)
          continue;
        OBAtom *atom1 = mol.GetAtom(a1+1);
        OBAtom *atom2 = mol.GetAtom(a2+1);
        if (atom1->IsHydrogen())
          continue;
        if (atom2->IsHydrogen())
          continue;
        // skip connected atoms
        if (mol.GetAtom(a1+1)->IsConnected(mol.GetAtom(a2+1)))
          continue;
        // compute the distance
        double dx = conformer[a1*3  ] - conformer[a2*3  ];
        double dy = conformer[a1*3+1] - conformer[a2*3+1];
        double dz = conformer[a1*3+2] - conformer[a2*3+2];
        double distance = sqrt(dx*dx + dy*dy + dz*dz);
        // check distance
        if (distance < m_cutoff)
          return false;
      }
    }

    return true;
  }

  //////////////////////////////////////////////////////////
  //
  //  OBConformerScore(s)
  //
  //////////////////////////////////////////////////////////

  double OBRMSDConformerScore::Score(OBMol &mol, unsigned int index,
      const RotorKeys &keys, const std::vector<double*> &conformers)
  {
    unsigned int numAtoms = mol.NumAtoms();

    // create the reference vector3 for conformer with index
    double *conformer_i = conformers[index];
    std::vector<vector3> vi;
    for (unsigned int a = 0; a < numAtoms; ++a)
      vi.push_back(vector3(conformer_i[a*3], conformer_i[a*3+1], conformer_i[a*3+2]));

    OBAlign align(mol, mol, false, false);
    align.SetRef(vi);

    double score_min = 10e10;
    for (unsigned int j = 0; j < conformers.size(); ++j) {
      if (index == j)
        continue;
      double *conformer_j = conformers[j];
      // create vector3 conformer
      std::vector<vector3> vj;
      for (unsigned int a = 0; a < numAtoms; ++a)
        vj.push_back(vector3(conformer_j[a*3], conformer_j[a*3+1], conformer_j[a*3+2]));

      // perform Kabsch alignment
      align.SetTarget(vj);
      align.Align();

      // get the RMSD
      double rmsd = align.GetRMSD();
      // store the rmsd if it is lower than any of the previous
      if (rmsd < score_min)
        score_min = rmsd;
    }

    // return the lowest RMSD
    return score_min;
  }

  double OBEnergyConformerScore::Score(OBMol &mol, unsigned int index,
      const RotorKeys &keys, const std::vector<double*> &conformers)
  {
    double *origCoords = mol.GetCoordinates();
    // copy the original coordinates to coords
    // copy the conformer coordinates to OBMol object
    std::vector<double> coords(mol.NumAtoms() * 3);
    for (unsigned int i = 0; i < mol.NumAtoms() * 3; ++i) {
      coords[i] = origCoords[i];
      origCoords[i] = conformers[index][i];
    }

    OBForceField *ff = OBForceField::FindType("MMFF94");
    if (!ff->Setup(mol))
      return 10e10;
//    ff->SteepestDescent(500);
    double score = ff->Energy();

    // copy original coordinates back
    for (unsigned int i = 0; i < mol.NumAtoms() * 3; ++i)
      origCoords[i] = coords[i];

    return score;
  }

  //////////////////////////////////////////////////////////
  //
  //  OBConformerSearch
  //
  //////////////////////////////////////////////////////////


  OBConformerSearch::OBConformerSearch()
  {
    m_filter = static_cast<OBConformerFilter*>(new OBStericConformerFilter(1.0));
    m_score = static_cast<OBConformerScore*>(new OBRMSDConformerScore());
    //m_score = static_cast<OBConformerScore*>(new OBEnergyConformerScore());
    m_numConformers = 30;
    m_numChildren = 5;
    m_mutability = 10;
  }

  OBConformerSearch::~OBConformerSearch()
  {
    delete m_filter;
    delete m_score;
  }





  bool OBConformerSearch::Setup(const OBMol &mol, int numConformers, int numChildren, int mutability, int convergence)
  {
    // copy some variables
    m_mol = mol;
    m_numConformers = numConformers;
    m_numChildren = numChildren;
    m_mutability = mutability;
    m_convergence = convergence;

    if (m_mol.GetCoordinates() == NULL)
      return false;

    // Initialize the OBRotorList
    m_rotorList.SetFixedBonds(m_fixedBonds);
    m_rotorList.Setup(m_mol);
    if (!m_rotorList.Size()) { // only one conformer
      return false;
    }

    // create initial population
    OBRandom generator;
    generator.TimeSeed();

    RotorKey rotorKey(m_rotorList.Size() + 1, 0); // indexed from 1
    if (IsGood(rotorKey))
      m_rotorKeys.push_back(rotorKey);
    else {
      cout << "Initial conformer does not pass filter!" << endl;
    }

    int tries = 0;
    while (m_rotorKeys.size() < m_numConformers && tries < numConformers * 1000) {
      tries++;
      // perform random mutation(s)
      OBRotorIterator ri;
      OBRotor *rotor = m_rotorList.BeginRotor(ri);
      for (unsigned int i = 1; i < m_rotorList.Size() + 1; ++i, rotor = m_rotorList.NextRotor(ri)) {
        if (generator.NextInt() % m_mutability == 0)
          rotorKey[i] = generator.NextInt() % rotor->GetResolution().size();
      }
      // duplicates are always rejected
      if (!IsUniqueKey(m_rotorKeys, rotorKey))
        continue;
      // execute filter(s)
      if (!IsGood(rotorKey))
        continue;
      // add the key
      m_rotorKeys.push_back(rotorKey);
    }

    // print out initial conformers
    cout << "Initial conformer count: " << m_rotorKeys.size() << endl;
    for (unsigned int i = 0; i < m_rotorKeys.size(); ++i) {
      for (unsigned int j = 1; j < m_rotorKeys[i].size(); ++j)
        std::cout << m_rotorKeys[i][j] << " ";
      std::cout << std::endl;
    }

    return true;
  }

  void OBConformerSearch::NextGeneration()
  {
    // create next generation population
    OBRandom generator;
    generator.TimeSeed();

    // generate the children
    int numConformers = m_rotorKeys.size();
    for (int c = 0; c < numConformers; ++c) {
      for (int child = 0; child < m_numChildren; ++child) {
        bool foundKey = false;
        int tries = 0;
        while (!foundKey) {
          tries++;
          if (tries > 1000)
            foundKey = true;
          RotorKey rotorKey = m_rotorKeys[c]; // copy parent gene
          // perform random mutation(s)
          OBRotorIterator ri;
          OBRotor *rotor = m_rotorList.BeginRotor(ri);
          for (unsigned int i = 1; i < m_rotorList.Size() + 1; ++i, rotor = m_rotorList.NextRotor(ri)) {
            if (generator.NextInt() % m_mutability == 0)
              rotorKey[i] = generator.NextInt() % rotor->GetResolution().size(); // permutate gene
          }
          // duplicates are always rejected
          if (!IsUniqueKey(m_rotorKeys, rotorKey))
            continue;
          // execute the filter(s)
          if (!IsGood(rotorKey))
            continue;
          // add the key
          m_rotorKeys.push_back(rotorKey); // append child to population
          // set foundKey to generate the next child
          foundKey = true;
        }
      }
    }
  }

  // Helper struct to sort conformers by score
  struct ConformerScore {
    ConformerScore(const RotorKey &key, double _score) : rotorKey(key), score(_score) {}
    RotorKey rotorKey;
    double score;
  };

  // Helper struct to compare ConformerScore objects by score
  struct CompareConformerHighScore {
    bool operator()(const ConformerScore &cs1, const ConformerScore &cs2) { return cs1.score > cs2.score; }
  };
  struct CompareConformerLowScore {
    bool operator()(const ConformerScore &cs1, const ConformerScore &cs2) { return cs1.score < cs2.score; }
  };


  double OBConformerSearch::MakeSelection()
  {
    OBRotamerList rotamers;
    rotamers.SetBaseCoordinateSets(m_mol);
    rotamers.Setup(m_mol, m_rotorList);

    // Add all (parent + children) rotor keys
    for (unsigned int i = 0; i < m_rotorKeys.size(); ++i) {
      rotamers.AddRotamer(m_rotorKeys[i]);
    }

    // Get conformers for the rotor keys
    std::vector<double*> conformers;
    rotamers.ExpandConformerList(m_mol, conformers);

    // Score each conformer
    std::vector<ConformerScore> conformer_scores;
    for (unsigned int i = 0; i < conformers.size(); ++i) {
      double score = m_score->Score(m_mol, i, m_rotorKeys, conformers);
      conformer_scores.push_back(ConformerScore(m_rotorKeys[i], score));
    }

    // delete the conformers
    for (unsigned int i = 0; i < conformers.size(); ++i) {
      delete [] conformers[i];
    }

    if (m_score->GetPreferred() == OBConformerScore::HighScore)
      std::sort(conformer_scores.begin(), conformer_scores.end(), CompareConformerHighScore());
    else
      std::sort(conformer_scores.begin(), conformer_scores.end(), CompareConformerLowScore());

    // Rmove the worst scored conformers until we have the disired number of conformers
    while (conformer_scores.size() > m_numConformers) {
      conformer_scores.pop_back();
    }

    // compute sum of all scores, this is a measure of convergence
    double sum = 0.0, lowest, highest;
    m_rotorKeys.clear();
    for (unsigned int i = 0; i < conformer_scores.size(); ++i) {
      switch (m_score->GetConvergence()) {
        case OBConformerScore::Highest:
          if (!i || conformer_scores[i].score > highest) {
            highest = conformer_scores[i].score;
          }
          break;
        case OBConformerScore::Lowest:
          if (!i || conformer_scores[i].score < lowest) {
            lowest = conformer_scores[i].score;
          }
          break;
        default:
          sum += conformer_scores[i].score;
          break;
      }
      // store the best keys
      m_rotorKeys.push_back(conformer_scores[i].rotorKey);
    }

    switch (m_score->GetConvergence()) {
      case OBConformerScore::Highest:
        return highest;
      case OBConformerScore::Lowest:
        return lowest;
      case OBConformerScore::Sum:
        return sum;
      case OBConformerScore::Average:
      default:
        return sum / m_rotorKeys.size();
    }
  }



  void OBConformerSearch::Search()
  {
    int identicalGenerations = 0;
    double last_score = 0.0;
    for (int i = 0; i < 1000; i++) {
      std::cout << "Generation #" << i + 1 << "  " << last_score << std::endl;
      // keep copy of rotor keys if next generation is less fit
      RotorKeys rotorKeys = m_rotorKeys;
      // create the children
      NextGeneration();
      // make the selection
      double score = MakeSelection();

      if (IsNear(last_score, score)) {
        identicalGenerations++;
        last_score = score;
      } else {
        if (score < last_score) {
          m_rotorKeys = rotorKeys;
          identicalGenerations++;
        } else {
          last_score = score;
          identicalGenerations = 0;
        }
      }

//      if (i)
  //      std::cerr << ";";
    //  std::cerr << last_score;

      if (identicalGenerations > m_convergence)
        break;
    }
//    std::cerr << std::endl;


    for (unsigned int i = 0; i < m_rotorKeys.size(); ++i) {
      for (unsigned int j = 1; j < m_rotorKeys[i].size(); ++j)
        std::cout << m_rotorKeys[i][j] << " ";
      std::cout << std::endl;
    }


  }

  void OBConformerSearch::GetConformers(OBMol &mol)
  {
    OBRotamerList rotamers;
    rotamers.SetBaseCoordinateSets(mol);
    rotamers.Setup(mol, m_rotorList);

    std::cout << "GetConformers:" << std::endl;
    // Add all (parent + children) unique rotor keys
    for (unsigned int i = 0; i < m_rotorKeys.size(); ++i) {
      rotamers.AddRotamer(m_rotorKeys[i]);

      for (unsigned int j = 1; j < m_rotorKeys[i].size(); ++j)
        std::cout << m_rotorKeys[i][j] << " ";
      std::cout << std::endl;
    }

    // Get conformers for the rotor keys
    std::vector<double*> conformers;
    rotamers.ExpandConformerList(mol, conformers);
    if (conformers.size())
      mol.SetConformers(conformers);
  }

  bool OBConformerSearch::IsUniqueKey(const RotorKeys &keys, const RotorKey &key) const
  {
    for (unsigned int i = 0; i < keys.size(); ++i)
      if (keys[i] == key)
        return false;
    return true;
  }

  bool OBConformerSearch::IsGood(const RotorKey &key)
  {
    // Setup OBRotamerList
    OBRotamerList rotamers;
    rotamers.SetBaseCoordinateSets(m_mol);
    rotamers.Setup(m_mol, m_rotorList);
    rotamers.AddRotamer(key);

    // Get conformer for the rotor keys
    std::vector<double*> conformers;
    rotamers.ExpandConformerList(m_mol, conformers);
    double *conformer = conformers[0];

    // Execute filter(s)
    bool result = m_filter->IsGood(m_mol, key, conformer);

    delete [] conformer;

    return result;
  }


  /**
   * @example obconformersearch_default.cpp
   */

};

/// @file conformersearch.cpp
/// @brief Conformer searching using genetic algorithm

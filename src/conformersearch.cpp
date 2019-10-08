/**********************************************************************
conformersearch.cpp - Conformer searching using genetic algorithm.

Copyright (C) 2010 Tim Vandermeersch
Some portions Copyright (C) 2016 Torsten Sachse

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
#include <openbabel/elements.h>
#include <openbabel/bond.h>
#include "rand.h"

#if defined(_MSC_VER) && (_MSC_VER < 1800)
 #define OB_ISNAN _isnan
#else
 #define OB_ISNAN std::isnan
#endif

namespace OpenBabel {

  //////////////////////////////////////////////////////////
  //
  //  OBConformerFilter(s)
  //
  //////////////////////////////////////////////////////////

  OBStericConformerFilter::OBStericConformerFilter ()
  {
    m_cutoff = 0.8;
    m_vdw_factor = 0.5;
    m_check_hydrogens = true;
  }

  OBStericConformerFilter::OBStericConformerFilter (double cutoff, double vdw_factor, bool check_hydrogens)
  {
    m_cutoff = cutoff * cutoff;
    m_vdw_factor = vdw_factor;
    m_check_hydrogens = check_hydrogens;
  }

  OBConformerFilter::~OBConformerFilter() {}

  bool OBStericConformerFilter::IsGood(const OBMol &mol, const RotorKey &key, double *conformer)
  {
    unsigned int a1 = 0, a2 = 0;
    unsigned int numAtoms = mol.NumAtoms();
    OBAtom *atom1 = NULL, *atom2 = NULL;
    double dx = 0.0, dy = 0.0, dz = 0.0;
    double distanceSquared = 0.0, vdwCutoff = 0.0;

    for (a1 = 0; a1 < numAtoms; ++a1) {
      for (a2 = a1 + 1; a2 < numAtoms; ++a2) {
        atom1 = mol.GetAtom(a1+1);
        atom2 = mol.GetAtom(a2+1);
        // Default should be to recognize H clashes too
        if (!m_check_hydrogens  && (atom1->GetAtomicNum() == OBElements::Hydrogen || atom2->GetAtomicNum() == OBElements::Hydrogen ))
          continue;

        // skip connected atoms
        if (atom1->IsConnected(atom2))
          continue;
        // compute the distance
        dx = conformer[a1*3  ] - conformer[a2*3  ];
        dy = conformer[a1*3+1] - conformer[a2*3+1];
        dz = conformer[a1*3+2] - conformer[a2*3+2];
        distanceSquared = dx*dx + dy*dy + dz*dz;
        // As we don't check 1-3 and 1-4 bonded atoms, apply a
        // factor of to the sum of VdW radii
        vdwCutoff = m_vdw_factor * (OBElements::GetVdwRad(atom1->GetAtomicNum())
                                    + OBElements::GetVdwRad(atom2->GetAtomicNum()));
        vdwCutoff *= vdwCutoff; // compare squared distances
        //std::cout << vdwCutoff << " " << m_vdw_factor << " " << m_cutoff << " " << distanceSquared << std::endl;

        // check distance
        if (distanceSquared < m_cutoff || distanceSquared < vdwCutoff)
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

  OBConformerScore::~OBConformerScore() {}

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
    energy_nrequest++;
    RotorKey cur_key = keys[index];
    if (energy_map.size () > 0)
      {
        // Check that we haven't already computed this energy);
        mapRotorEnergy::iterator it = energy_map.find (cur_key);
        if (it != energy_map.end ())
          return it->second;
      }
    energy_ncompute++;

    double *origCoords = mol.GetCoordinates();
    // copy the original coordinates to coords
    // copy the conformer coordinates to OBMol object
    std::vector<double> coords(mol.NumAtoms() * 3);
    for (unsigned int i = 0; i < mol.NumAtoms() * 3; ++i) {
      coords[i] = origCoords[i];
      origCoords[i] = conformers[index][i];
    }

    OBForceField *ff = OBForceField::FindType("MMFF94");
    if (!ff->Setup(mol)) {
      ff = OBForceField::FindType("UFF");
      if (!ff->Setup(mol))
        return 10e10;
    }
    double score = ff->Energy(false); // no gradients

    // copy original coordinates back
    for (unsigned int i = 0; i < mol.NumAtoms() * 3; ++i)
      origCoords[i] = coords[i];

    // Save that in the map
    if (energy_map.size () < 50000)
      energy_map[cur_key] = score;

    return score;
  }

  double OBMinimizingEnergyConformerScore::Score(OBMol &mol, unsigned int index,
                                                 const RotorKeys &keys, const std::vector<double*> &conformers)
  {
    energy_nrequest++;
    RotorKey cur_key = keys[index];
    if (energy_map.size () > 0)
      {
        // Check that we haven't already computed this energy);
        mapRotorEnergy::iterator it = energy_map.find (cur_key);
        if (it != energy_map.end ())
          return it->second;
      }
    energy_ncompute++;

    double *origCoords = mol.GetCoordinates();
    // copy the original coordinates to coords
    // copy the conformer coordinates to OBMol object
    std::vector<double> coords(mol.NumAtoms() * 3);
    for (unsigned int i = 0; i < mol.NumAtoms() * 3; ++i) {
      coords[i] = origCoords[i];
      origCoords[i] = conformers[index][i];
    }

    OBForceField *ff = OBForceField::FindType("MMFF94");
    if (!ff->Setup(mol)) {
      ff = OBForceField::FindType("UFF");
      if (!ff->Setup(mol))
        return 10e10;
    }
    ff->ConjugateGradients(50);
    double score = ff->Energy(false); // no gradients

    // copy original coordinates back
    for (unsigned int i = 0; i < mol.NumAtoms() * 3; ++i)
      origCoords[i] = coords[i];

    // Save that in the map
    if (energy_map.size () < 50000)
      energy_map[cur_key] = score;

    return score;
  }

  double OBMinimizingRMSDConformerScore::Score(OBMol &mol, unsigned int index,
                                               const RotorKeys &keys, const std::vector<double*> &conformers)
  {
    unsigned int numAtoms = mol.NumAtoms();
    double *origCoords = mol.GetCoordinates();
    // copy the original coordinates to coords
    // copy the conformer coordinates to OBMol object
    std::vector<double> coords(mol.NumAtoms() * 3);
    for (unsigned int i = 0; i < mol.NumAtoms() * 3; ++i) {
      coords[i] = origCoords[i];
      origCoords[i] = conformers[index][i];
    }

    OBForceField *ff = OBForceField::FindType("MMFF94");
    if (!ff->Setup(mol)) {
      ff = OBForceField::FindType("UFF");
      if (!ff->Setup(mol))
        return 10e10;
    }
    ff->ConjugateGradients(50);
    double score = ff->Energy(false); // no gradients

    // copy original coordinates back
    for (unsigned int i = 0; i < mol.NumAtoms() * 3; ++i)
      origCoords[i] = coords[i];

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

  //////////////////////////////////////////////////////////
  //
  //  OBConformerSearch
  //
  //////////////////////////////////////////////////////////


  OBConformerSearch::OBConformerSearch()
  {
    m_filter = static_cast<OBConformerFilter*>(new OBStericConformerFilter());
    m_score = static_cast<OBConformerScore*>(new OBRMSDConformerScore());
    //m_score = static_cast<OBConformerScore*>(new OBEnergyConformerScore());
    m_numConformers = 30;
    m_numChildren = 5;
    m_mutability = 10;
    // Set some inital values here, but will fine-tune some in the molecule Setup
    use_sharing = false;
    alpha_share = 2.0;
    sigma_share = 2.0;
    nb_niches = 5;
    niche_radius = 1.0;
    p_crossover = 0.7;
    niche_mating = 0.7;
    local_opt_rate = 3;
    // For the moment 'd' is an opaque pointer to an instance of OBRandom*.
    // In future, it could be a pointer to a structure storing all of the
    // private variables.
    d = (void*)new OBRandom();
    ((OBRandom*)d)->TimeSeed();
    m_logstream = &std::cout; 	// Default logging send to standard output
    // m_logstream = NULL;
    m_printrotors = false;  // By default, do not print rotors but perform the conformer search

  }

  OBConformerSearch::~OBConformerSearch()
  {
    delete m_filter;
    delete m_score;
    delete (OBRandom*)d;
  }


  bool OBConformerSearch::Setup(const OBMol &mol, int numConformers, int numChildren, int mutability, int convergence)
  {
    int nb_rotors = 0;
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

    // Print all available rotors if so desired
    if (m_printrotors){
      OBRotorIterator it;
      OBRotorIterator end_it = m_rotorList.EndRotors();
      OBRotor* r = m_rotorList.BeginRotor(it);
      int rotcount = 1;
      std::cout << "Rotors:" << std::endl;
      while(r){
        OBBond* b = r->GetBond();
        int at1,at2;
        at1 = b->GetBeginAtomIdx();
        at2 = b->GetEndAtomIdx();
        std::cout << at1 << "-" << at2 << "  ";
        r = m_rotorList.NextRotor(it);
        if (rotcount%4==0 && r){std::cout << std::endl;}
        ++rotcount;
      }
      std::cout << std::endl;
      return false;
    }
    // Print those that are fixed
    if (!m_fixedBonds.IsEmpty()){
      std::cout << "Fixed Rotors: " << std::endl;
      int end_it = m_fixedBonds.EndBit();
      int it = m_fixedBonds.FirstBit();
      int rotcount = 1;
      while(it!=end_it){
        OBBond* b = m_mol.GetBond(it);
        int at1,at2;
        at1 = b->GetBeginAtomIdx();
        at2 = b->GetEndAtomIdx();
        std::cout << at1 << "-" << at2 << "  ";
        it = m_fixedBonds.NextBit(it);
        if (rotcount%4==0 && it!=end_it){std::cout << std::endl;}
        ++rotcount;
      }
      std::cout << std::endl;
    }

    nb_rotors = m_rotorList.Size();
    if (!nb_rotors) { // only one conformer
      return false;
    }

    // create initial population
    OBRandom generator;
    generator.TimeSeed();

    RotorKey rotorKey(m_rotorList.Size() + 1, 0); // indexed from 1
    if (IsGood(rotorKey))
      m_rotorKeys.push_back(rotorKey);
    else {
      if (m_logstream != NULL)
        (*m_logstream) << "Initial conformer does not pass filter!" << std::endl;
    }

    int tries = 0, ndup = 0, nbad = 0;
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
        {
          ndup++;
          continue;
        }
      // execute filter(s)
      if (!IsGood(rotorKey))
        {
          nbad++;
          continue;
        }
      // add the key
      m_rotorKeys.push_back(rotorKey);
    }

    // print out initial conformers
    if (m_logstream != NULL)
      {
        (*m_logstream) << "Initial conformer count: " << m_rotorKeys.size() << std::endl;
        (*m_logstream) << tries << " attempts,  " << ndup << " duplicates, " << nbad << " failed filter." << std::endl;
        for (unsigned int i = 0; i < m_rotorKeys.size(); ++i) {
          for (unsigned int j = 1; j < m_rotorKeys[i].size(); ++j)
            (*m_logstream) << m_rotorKeys[i][j] << " ";
          (*m_logstream) << std::endl;
        }
      }

    // Setup some default values for dynamic niche sharing according to the molecule.
    nb_niches = (m_rotorKeys.size()) / 10;
    if (nb_niches < 3)
      nb_niches = 3;
    sigma_share = (double)nb_rotors / 3.0;
    if (sigma_share < 1.0)
      sigma_share = 1.0;
    niche_radius =  (double)nb_rotors / 4.0;
    if (niche_radius < 1.0)
      niche_radius = 1.0;

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
    double last_score = 0.0, score = 0.0;

    if (m_logstream != NULL)
      {
        (*m_logstream) << std::endl << "=====> Starting conformers search with a Genetic Algorithm <=====" << std::endl;
        if (use_sharing)
          {
            (*m_logstream) << "Uses fitness sharing (with dynamic niche identification)" << std::endl;
            (*m_logstream) << "Population size :" << m_rotorKeys.size() << std::endl;
            (*m_logstream) << nb_niches << " niches searched, with a key distance radius of " << niche_radius << std::endl;
            (*m_logstream) << "Fitness sharing parameter alpha: " << alpha_share << " \t sigma:" << sigma_share << std::endl;
            (*m_logstream) << "Uniform crossover probability: " << p_crossover << std::endl;
            (*m_logstream) << "Mutation probability: " << (1.0 / (double) m_mutability) << std::endl;
            (*m_logstream) << "Niche mating probability: " << niche_mating << std::endl;
            if (local_opt_rate)
              {
                (*m_logstream) << "Trying to improve best indivual with local search every ";
                (*m_logstream) << local_opt_rate<< "generations" << std::endl;
              }
          }
        else
          {
            (*m_logstream) << "Perform elitist generation replacement with mutation only" << std::endl;
            (*m_logstream) << "Mutation probability: " << (1.0 / (double) m_mutability) << std::endl;
          }
        (*m_logstream) << "Will stop after " << m_convergence << " generations without improvement." << std::endl << std::endl;
      }
    if (use_sharing)
      score_population ();

    for (int i = 0; i < 1000; i++) {
      // keep copy of rotor keys if next generation is less fit
      RotorKeys rotorKeys = m_rotorKeys;

      if (use_sharing)
        {
          if (local_opt_rate && (i % local_opt_rate) == 0)
            local_opt ();		// Try to locally improve the best from times to times
          share_fitness ();
          score = sharing_generation ();
        }
      else
        {
          // create the children
          NextGeneration();
          // make the selection
          score = MakeSelection();
        }
      if (OB_ISNAN(score)) {
          (*m_logstream) << "The current score is not a number, will not continue." << std::endl << std::endl;
          return;
      }
      if (i == 0)
        last_score = score;

      if (IsNear(last_score, score)) {
        identicalGenerations++;
        last_score = score;
      } else {
        if (m_score->GetPreferred() == OBConformerScore::HighScore)
          {
            // Maximize score
            if (score < last_score) {
              if (!use_sharing)
                m_rotorKeys = rotorKeys;
              identicalGenerations++;
            } else {
              last_score = score;
              identicalGenerations = 0;
            }
          }
        else
          {
            // Minimize score
            if (score > last_score) {
              if (!use_sharing)
                m_rotorKeys = rotorKeys;
              identicalGenerations++;
            } else {
              last_score = score;
              identicalGenerations = 0;
            }
          }
      }
      if (m_logstream != NULL)
        {
          if (vscores.size ())
            (*m_logstream) << "Generation #" << i + 1 << "  " << last_score << "\t best " << vscores[0] << std::endl;
          else
            (*m_logstream) << "Generation #" << i + 1 << "  " << last_score << std::endl;
        }
      if (identicalGenerations > m_convergence)
        break;
    }

    if (m_logstream != NULL)
      {
        for (unsigned int i = 0; i < m_rotorKeys.size(); ++i) {
          for (unsigned int j = 1; j < m_rotorKeys[i].size(); ++j)
            (*m_logstream) << m_rotorKeys[i][j] << " ";
          (*m_logstream) << std::endl;
        }
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


  /*! @brief   Genetic similarity measure, i.e. "distance" between two rotor keys.

    @param key1: compute distance from
    @param key2: compute distance to

    @return: number of different alleles, i.e. number of different values in the key vectors
  */
  int
  OBConformerSearch::key_distance (const RotorKey &key1, const RotorKey &key2)
  {
    int dist = 0;
    //    assert(key1.size() > 1 && key1.size()== key2.size());
    std::vector<int>::const_iterator it1 = key1.begin ();
    std::vector<int>::const_iterator it2 = key2.begin ();
    // Skip first values, since  meaningfull valaues are starting at index 1 (Fortran translation inside ;-))
    for (++it1, ++it2; it1 != key1.end ();++it1, ++it2)
      if (*it1 != *it2)
        dist++;
    return dist;
  }

  /*! @brief Make a local search on the best individual.
    Randomly change each rotor key in the key vector.
    This is not a complete coverage of all "neighbors" vectors,
    since not all possible key values are tested.
  */
  int
  OBConformerSearch::local_opt ()
  {
    bool max_flag = (m_score->GetPreferred() == OBConformerScore::HighScore);
    bool flag_improved = false;
    std::vector<double> backup_scores = vscores;
    RotorKey best = m_rotorKeys[0], neighbor, opt_key;
    RotorKeys backup_population = m_rotorKeys;
    double opt_score = 0.0;
    int i = 0, new_val = 0;
    OBRotorIterator ri;
    OBRotor *rotor = NULL;

    // Change 1 and only value, looping on all positions
    m_rotorKeys.clear ();
    rotor = m_rotorList.BeginRotor(ri);
    for (i = 1; i <= m_rotorList.Size(); ++i, rotor = m_rotorList.NextRotor(ri))
      {
        neighbor = best;
        new_val = ((OBRandom*)d)->NextInt() % rotor->GetResolution().size();
        while (new_val == best[i])
          new_val = ((OBRandom*)d)->NextInt() % rotor->GetResolution().size();
        neighbor[i] = new_val;
        if (IsUniqueKey(backup_population, neighbor) && IsGood(neighbor))
          m_rotorKeys.push_back (neighbor);
      }
    score_population ();
    opt_score = vscores[0];
    if ( (max_flag && (opt_score > backup_scores[0]))  || (!max_flag && (opt_score < backup_scores[0])) )
      {				// Found a better conformer
        opt_key = m_rotorKeys[0];
        flag_improved = true;
        if (m_logstream != NULL)
          {
            (*m_logstream) << "    => Best individual improved with local search: ";
            (*m_logstream) << backup_scores[0] << "  --> " << opt_score << std::endl;
          }
      }
    // Set back population and score vector
    m_rotorKeys.clear ();
    m_rotorKeys = backup_population;
    vscores.clear ();
    vscores = backup_scores;
    if (flag_improved)
      { // Replace best one
        m_rotorKeys[0] = opt_key;
        vscores[0] = opt_score;
      }

    return (int) flag_improved;
  }

  /*! @brief  Produces one or two offsprings

    @param key1: reference to the first offspring key
    @param key2: reference to the second offspring key

    @return: 0 if no valid offspring, 1 if first only, 2 if second only, and 3 if 2 are valid offsprings.
  */
  int
  OBConformerSearch::reproduce (RotorKey &key1, RotorKey &key2)
  {
    unsigned int i = 0, iniche = 0, j = 0, pop_size = vscores.size ();
    unsigned int rnd1 = 0, rnd2 = 0, parent1 = 0, parent2 = 0, nsize = 0;
    int ret_code = 0;
    bool flag_crossover = false;
    OBRotorIterator ri;
    OBRotor *rotor = NULL;

    if (pop_size < 2)
      return 0;

    // Make a 2-tournament selection to choose first parent
    i = ((OBRandom*)d)->NextInt() % pop_size;
    j = ((OBRandom*)d)->NextInt() % pop_size;
    parent1 = vshared_fitnes[i] > vshared_fitnes[j] ? i : j;
    iniche = niche_map[parent1];
    if (iniche > -1)
      nsize = dynamic_niches[iniche].size (); // Belongs to a specific niche

    // Do we apply crossover here?
    flag_crossover = (((OBRandom*)d)->NextFloat () <= p_crossover);
    if (flag_crossover && (((OBRandom*)d)->NextFloat () <= niche_mating)  &&  nsize > 1)
      {
        // Apply niche mating: draw second parent in the same niche, if its has
        // at least 2 members. Make a 2-tournament selection whithin this niche
        rnd1 = ((OBRandom*)d)->NextInt() % nsize;
        i =  dynamic_niches[iniche][rnd1];
        rnd2 = ((OBRandom*)d)->NextInt() % nsize;
        j = dynamic_niches[iniche][rnd2];
        parent2 = vshared_fitnes[i] > vshared_fitnes[j] ? i : j;
      }
    else
      {
        // Draw second in the whole population
        i = ((OBRandom*)d)->NextInt() % pop_size;
        j = ((OBRandom*)d)->NextInt() % pop_size;
        parent2 = vshared_fitnes[i] > vshared_fitnes[j] ? i : j;
      }

    if (flag_crossover)
      {
        // Cross the 2 vectors: toss a coin for each position (i.e. uniform crossover)
        for (i = 1; i < key1.size(); i++)
          {
            if (((OBRandom*)d)->NextInt() % 2)
              { // Copy parent1 to offspring 1
                key1[i] = m_rotorKeys[parent1][i];
                key2[i] = m_rotorKeys[parent2][i];
              }
            else
              { // invert
                key2[i] = m_rotorKeys[parent1][i];
                key1[i] = m_rotorKeys[parent2][i];
              }
          }
      }
    else
      {
        key1 =  m_rotorKeys[parent1];
        key2 =  m_rotorKeys[parent2];
      }

    // perform random mutation(s)
    rotor = m_rotorList.BeginRotor(ri);
    for (i = 1; i <= m_rotorList.Size(); ++i, rotor = m_rotorList.NextRotor(ri))
      {
        if (((OBRandom*)d)->NextInt() % m_mutability == 0)
          key1[i] = ((OBRandom*)d)->NextInt() % rotor->GetResolution().size();
        if (((OBRandom*)d)->NextInt() % m_mutability == 0)
          key2[i] = ((OBRandom*)d)->NextInt() % rotor->GetResolution().size();
      }
    if (IsUniqueKey(m_rotorKeys, key1) && IsGood(key1))
      ret_code += 1;
    if (IsUniqueKey(m_rotorKeys, key2) && IsGood(key2))
      ret_code += 2;
    return ret_code;		// Returns number of new distinct individuals (i.e. rotor keys)
  }

  /*! @brief  Score and sort the current population

    @return: 0 when OK, 1 or more if not
  */
  int
  OBConformerSearch::score_population ()
  {
    bool max_flag = (m_score->GetPreferred() == OBConformerScore::HighScore);
    unsigned int i = 0, pop_size = 0;
    double score = 0.0;
    std::vector<double*> conformers;
    std::vector<double>::iterator dit;
    OBRotamerList rotamers;
    std::vector<ConformerScore> conformer_scores;

    rotamers.SetBaseCoordinateSets(m_mol);
    rotamers.Setup(m_mol, m_rotorList);

    // Add all (parent + children) rotor keys
    for (i = 0; i < m_rotorKeys.size(); ++i)
      rotamers.AddRotamer(m_rotorKeys[i]);

    // Get conformers for the rotor keys
    rotamers.ExpandConformerList(m_mol, conformers);

    // Score each conformer
    for (i = 0; i < conformers.size(); ++i)
      {
        score = m_score->Score(m_mol, i, m_rotorKeys, conformers);
        conformer_scores.push_back(ConformerScore(m_rotorKeys[i], score));
      }

    // delete the conformers
    for (i = 0; i < conformers.size(); ++i)
      delete [] conformers[i];

    // Sort conformers
    if (max_flag)
      std::sort(conformer_scores.begin(), conformer_scores.end(), CompareConformerHighScore());
    else
      std::sort(conformer_scores.begin(), conformer_scores.end(), CompareConformerLowScore());
    pop_size = conformer_scores.size();

    // Save the score vector, reorder population
    vscores.clear ();
    m_rotorKeys.clear();
    for (i = 0; i < pop_size; i++)
      {
        vscores.push_back (conformer_scores[i].score);
        m_rotorKeys.push_back(conformer_scores[i].rotorKey);
      }
    return 0;
  }


  /*! @brief   Compute shared fitness values for a given population.
    Assume the population was order prior to this call.

    @return: 0 when OK, 1 or more if not
  */
  int
  OBConformerSearch::share_fitness ()
  {
    bool max_flag = (m_score->GetPreferred() == OBConformerScore::HighScore);
    bool flag_niched = false;
    unsigned int i = 0, iniche = 0, j = 0, k = 0, pop_size = vscores.size ();
    unsigned int max_dist = 0, min_dist = 0, imax = 0;
    int dist = 0;
    double min_score = 0.0, dtmp = 0.0;
    double sh_count = 0.0, sh_ij = 0.0;
    std::vector<double>::iterator dit;
    std::vector<int>::iterator nit;
    std::vector<char> vshared;

    min_score = max_flag?  vscores[pop_size - 1] :  vscores[0];

    // Compute shared fitness from score: rescale score values so the minimum is 1.0
    vshared_fitnes.clear ();
    dtmp = 1.0 - min_score;
    if (max_flag)
      for (dit =  vscores.begin (); dit != vscores.end (); ++dit)
        vshared_fitnes.push_back (*dit + dtmp);
    else			// Minimization: invert score so best has 1.0, others lower values
      for (dit =  vscores.begin (); dit != vscores.end (); ++dit)
        vshared_fitnes.push_back (1.0 / (*dit + dtmp));

    // Use a vector to keep track of assigned indivivduals
    vshared.resize (pop_size);
    // Identify niches
    dynamic_niches.clear ();
    for (k = 0; k < pop_size; k++)
      {
        imax = -1;
        if ((dynamic_niches.size () < nb_niches) && (dynamic_niches.size () %2))
          {
            // Get the most distant individual to current list heads, wihtin the first 2/3 of the population
            // The idea is to find some explicit niches, with maximum disimilarity form existing ones,
            // still with some reasonable head fitness score.
            max_dist = 0;
            for (i = 0; i < ((pop_size * 2) / 3); i++)
              {
                if (vshared[i])
                  continue;
                min_dist = 1000000;
                for(iniche = 0; iniche < dynamic_niches.size (); iniche++)
                  {
                    j = dynamic_niches[iniche][0];
                    dist = key_distance (m_rotorKeys[j], m_rotorKeys[i]);
                    if (dist < min_dist )
                      min_dist = dist;
                  }
                if (min_dist > max_dist)
                  {
                    imax = i;
                    max_dist = min_dist;
                  }
              }
            i = imax;
          }

        // Get the best unassigned individual
        if (imax == -1)
          for (i = 0; i < pop_size; i++)
            if (!vshared[i])
              break;

        for (iniche = 0; iniche < dynamic_niches.size (); iniche++)
          {
            j = dynamic_niches[iniche][0];
            dist = key_distance(m_rotorKeys[j], m_rotorKeys[i]);
            if ((double)dist <= niche_radius)
              {
                dynamic_niches[iniche].push_back (i);
                break;
              }
          }
        if (iniche == dynamic_niches.size ())
          {  // Could not find a niche for this one
            if (dynamic_niches.size () < nb_niches)
              { // Create a new niche.
                iniche = dynamic_niches.size ();
                dynamic_niches.resize (iniche + 1);
                dynamic_niches[iniche].push_back (i);
              }
            else
              { // Do it the hard way: compute sharing factor with all the population
                sh_count = 0.0;
                for (j = 0; j < pop_size; j++)
                  {
                    dist = key_distance (m_rotorKeys[i], m_rotorKeys[j]);
                    if ((double)dist < sigma_share)
                      {
                        sh_ij = 1.0 - pow (((double) dist) / ((double ) sigma_share), alpha_share);
                        sh_count += sh_ij;
                      }
                  }
                vshared_fitnes[i] /= sh_count; // Classical shared fitness
              }
          }
        vshared[i] = 1;
      }

    // Build the niche map: provided an invididual index, provides the niche index it belongs to.
    // -1 means out of explicit niches.
    niche_map.clear ();
    niche_map.resize (pop_size, -1);
    for (iniche = 0; iniche < dynamic_niches.size (); iniche++)
      {
        // Divide each dynamic niche member score by the niche size: i.e. each niche member has the same
        // penalty, i.e. the order is unchanged whitin this niche.
        dtmp = 1.0 / ((double) dynamic_niches[iniche].size ());
        for (nit = dynamic_niches[iniche].begin (); nit != dynamic_niches[iniche].end (); ++nit)
          {
            vshared_fitnes[*nit] *= dtmp;
            niche_map[*nit] = iniche;
          }
      }

    return 0;
  }

  /*! @brief  Create a new generation with fitness sharing

    @return: generation score (lowest, highest, average or sum of fitness)
  */
  double
  OBConformerSearch::sharing_generation ()
  {
    RotorKeys new_generation, offsprings;
    RotorKey key1, key2;
    unsigned int i = 0, j = 0, imax = 0, iniche = 0, nsize = 0, pop_size = 0, half_pop = 0, key_size = 0;
    int ret_code = 0,  ninsert = 0, iround = 0, last_ninsert;
    int dist = 0, dist_max = 0, dist_min = 0;
    double sum = 0.0, lowest = 0.0, highest = 0.0;
    std::vector<int>::iterator nit;
    std::vector<double>::iterator dit;
    std::vector<unsigned char> vflag;
    std::vector<int> out_niches; // Vector of all individuals out of the niches

    // key1 and key2 will be used to create new keys for the new generation
    key_size = m_rotorKeys[0].size();
    key1.resize (key_size);
    key2.resize (key_size);

    pop_size = vscores.size ();
    vflag.resize (pop_size);
    half_pop = pop_size / 2;
    if (pop_size % 2)
      half_pop++;

    // Build the list of individuals out of the niches
    for (i = 0; i < pop_size; i++)
      if (niche_map[i] == -1)
        out_niches.push_back (i);

    if (m_logstream != NULL)
      {
        (*m_logstream) << "  ==> Number of niches: " << dynamic_niches.size ();
        (*m_logstream) << "   # out of niches :" << out_niches.size () << "\t Best :" << vscores[0] << std::endl;
      }

    // Save each niche 1set element, then 2nd until we have half of the population
    // If we still miss some, get from individuals out of the niches (starting from best ones)
    new_generation.push_back (m_rotorKeys[0]);
    vflag[0] = 1;
    ninsert = 1;
    iround = 0;
    while (ninsert < ((half_pop * 2) / 3))
      {
        last_ninsert = ninsert;
        for (iniche = 0; iniche < dynamic_niches.size (); iniche++)
          {
            nsize = dynamic_niches[iniche].size ();
            imax = nsize / 2;
            if (nsize % 2)
              imax++;
            if (iround < imax)
              {
                i = dynamic_niches[iniche][iround];
                if (!vflag[i])
                  {
                    new_generation.push_back (m_rotorKeys[i]);
                    vflag[i] = 1;
                    ninsert++;
                    if (ninsert >= half_pop)
                      break;
                  }
              }
          }
        if (ninsert == last_ninsert)
          {
            for (i = 0; i < out_niches.size () ; i++)
              {
                j = out_niches[i];
                if (!vflag[j])
                  {
                    new_generation.push_back (m_rotorKeys[j]);
                    vflag[j] = 1;
                    ninsert++;
                    break;

                  }
              }
          }
        if (ninsert == last_ninsert)
          break;
        iround++;
      }
    // Now add most distant individuals: force genetic diversity in the new population at the cost of higher energies.
    while (ninsert < half_pop)
      {
        dist_max = 0;
        imax = 0;
        for (i = 0; i < pop_size; i++)
          {
            if (vflag[i])
              continue;
            dist_min = 1000000;
            for (j = 1; j < pop_size; j++)
              {
                // Get closest neigbhor distance for this key
                if (!vflag[j])
                  continue;
                dist = key_distance (m_rotorKeys[i], m_rotorKeys[j]);
                if (dist < dist_min)
                  dist_min = dist;
              }
            if (dist_min > dist_max)
              {			// Find the most distant to its clostest neighbor
                imax = i;
                dist_max = dist_min;
              }
          }
        new_generation.push_back (m_rotorKeys[imax]);
        vflag[imax] = 1;
        ninsert++;
      }

    // Now generate offsprings
    offsprings.clear ();
    while (offsprings.size () < pop_size)
      {
        ret_code = reproduce (key1, key2);
        // Add the first offspring if valid... and not already in the new generation
        if ((ret_code % 2) && IsUniqueKey(new_generation, key1))
          offsprings.push_back (key1);
        // Add the second if valid, and enough space in the new population
        if ((ret_code > 1) && IsUniqueKey(new_generation, key2))
          offsprings.push_back (key2);
      }

    // Score the offsprings and the best in the new population
    m_rotorKeys.clear ();
    for (i = 0; i < offsprings.size(); i++)
      m_rotorKeys.push_back(offsprings[i]);
    score_population ();
    imax = pop_size - new_generation.size ();
    for (i = 0; i < imax; i++)
      new_generation.push_back(m_rotorKeys[i]);
    m_rotorKeys.clear ();
    for (i = 0; i < new_generation.size(); i++)
      m_rotorKeys.push_back(new_generation[i]);
    score_population ();

    // Convergence energy values: average of niches best element
    sum = 0.0;
    if (0) {
      for (iniche = 0; iniche <  dynamic_niches.size (); iniche++)
        {
          i = dynamic_niches[iniche][0];
          sum += vscores[i];
        }
      return sum / ((double) dynamic_niches.size ());
    }
    imax = pop_size / 2;
    for (i = 0; i < imax; i++)
      sum += vscores[i];

    // Convergence criterion: best half population score
    return sum / ((double) imax);
  }

  /**
   * @example obconformersearch_default.cpp
   */

};

/// @file conformersearch.cpp
/// @brief Conformer searching using genetic algorithm

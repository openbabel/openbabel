/**********************************************************************
confsearch.cpp - Conformer Searching Routines (see also forcefield.cpp)

Copyright (C) 2010 Noel O'Boyle <baoilleach@gmail.com>

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

#include <openbabel/babelconfig.h>

#include <openbabel/forcefield.h>
#include <openbabel/rotamer.h>
#include <openbabel/rotor.h>
#include <openbabel/math/align.h>
#include <openbabel/tree/tree.hh>
#include <openbabel/tree/tree_util.hh>
#include <openbabel/math/vector3.h>

#include <float.h> // For DBL_MAX
#include <algorithm> // For min
#include <limits.h> // For UINTS_MAX with certain old GCC4

#include <iomanip> // For setprecision

namespace OpenBabel
{

  class LFSR
  {
    /**
     * Usage:
     *     LFSR lfsr(N); // where N < 2^31
     *     unsigned int d;
     *     do {
     *       d = lfsr.GetNext();
     *       // .... do something with d ...
     *     } while (d != 1);
     **/
  public:
    LFSR(unsigned int range, unsigned int start);
    unsigned int GetNext(); // Return 1 when finished
  private:
    unsigned int _range, _lfsr, _poly;
    static const unsigned int _polynomials[31];
  };

  const unsigned int LFSR::_polynomials[31] = {0x3, 0x6, 0xc, 0x14, 0x30, 0x60, 0xb8, 0x110, 0x240, 0x500,
       0x829, 0x100d, 0x2015, 0x6000, 0xd008, 0x12000, 0x20400, 0x40023, 0x90000, 0x140000,
       0x300000, 0x420000, 0xe10000, 0x1200000, 0x2000023, 0x4000013, 0x9000000, 0x14000000,
       0x20000029, 0x48000000, 0x80200003};

  LFSR::LFSR(unsigned int range, unsigned int start = 1): _range(range), _lfsr(start)
  {
    assert ( _range < (1 << 31) ); // We can only handle up to 2^31 - 1
    // (Currently don't use a start value)
    // assert ( start < _range ); // Otherwise the _start value will never be returned

    int i = 0;
    unsigned int tot = 4;
    while (tot <= _range) {
      i++;
      tot <<= 1;
    }
    _poly = _polynomials[i];
  }

  inline unsigned int LFSR::GetNext()
  {
    do {
      _lfsr = (_lfsr >> 1) ^ (unsigned int)(0 - (_lfsr & 1u) & _poly);
    } while (_lfsr > _range);

    return _lfsr;
  }

  class OBDiversePoses {
    public:
      OBDiversePoses(const OBMol &ref, double RMSD, bool percise=false);
      ~OBDiversePoses() {
        delete palign; // Allocated with 'new'
      }
      bool AddPose(double* coords, double energy);
      bool AddPose(vector<vector3> coords, double energy);
      typedef std::pair<vector<vector3>, double> PosePair;
      typedef tree<PosePair> Tree;
      Tree* GetTree() { return &poses; }
      typedef tree<PosePair>::iterator Tree_it;
      typedef tree<PosePair>::sibling_iterator Tree_sit;
      size_t GetSize();
      inline int GetNRMSD() {
        return n_rmsd;
      }
      inline double GetCutoff() {
        return cutoff;
      }

    private:
      bool _percise;
      vector<vector3> GetHeavyAtomCoords(const vector<vector3> &all_coords);
      int natoms;
      Tree poses;
      std::vector<double> levels;
      OBAlign* palign;
      const double cutoff;
      int n_rmsd;
      OBBitVec hydrogens;
    };

  OBDiversePoses::OBDiversePoses(const OBMol &ref, double RMSD, bool percise):
          palign(new OBAlign(false, percise)), cutoff(RMSD), _percise(percise)
  {
    natoms = ref.NumAtoms();
    palign->SetRefMol(ref);
    palign->SetMethod(OBAlign::QCP);
    n_rmsd = 0;

    const double arr[] = {3.0, 2.0, 1.5, 1.0, 0.5, 0.25};
    std::vector<double> vec (arr, arr + sizeof(arr) / sizeof(arr[0]) );
    vec.erase(std::remove_if(vec.begin(), vec.end(), std::bind2nd(std::less<double>(), (cutoff + 0.1) )), vec.end());
    vec.push_back(cutoff);

    levels = vec;
    vector<vector3> pdummy;
    poses.insert(poses.begin(), PosePair(pdummy, 0.0)); // Add a dummy top node

    // Remember the hydrogens
    hydrogens.Resize(natoms);
    for (int i=1; i<=natoms; i++)
      if (ref.GetAtom(i)->IsHydrogen())
        hydrogens.SetBitOn(i - 1);
  }

  bool OBDiversePoses::AddPose(double* coords, double energy) {
    vector<vector3> vcoords, vcoords_hvy;
    for (unsigned int a = 0; a < natoms; ++a)
      vcoords.push_back(vector3(coords[a*3], coords[a*3+1], coords[a*3+2]));
    return AddPose(vcoords, energy);
  }

  bool OBDiversePoses::AddPose(vector<vector3> vcoords, double energy) {
    Tree_it node = poses.begin();
    int level = 0;
    bool first_time = true;

    // Convert coords to vector<vector3>

    // Only use the heavy-atom coords for the alignment, but store
    // the full set of coordinates in the tree
    vector<vector3> vcoords_hvy = GetHeavyAtomCoords(vcoords);
    palign->SetRef(vcoords_hvy);

    vector<Tree_it> nodes, min_nodes;
    vector<double> min_nodes_rmsds;
    nodes.push_back(node);
    vector<int> stack_levels;
    stack_levels.push_back(level);

    vector<Tree_it> insert_pt;
    vector<PosePair> insert_data;
    vector<int> insert_level;

    while(nodes.size() > 0) { // Using stack-based recursion
      node = nodes.back();
      nodes.pop_back();
      level = stack_levels.back();
      stack_levels.pop_back();

      // Find whether the molecule is similar to any of the children of this node.
      // - min_node will hold the result of this search

      min_nodes.clear();
      min_nodes_rmsds.clear();
      double rmsd;
      double min_rmsd = DBL_MAX;

      Tree_sit sib = poses.begin(node);
      // Skip the first child after the first time through this loop
      // - it will already have been tested against at the end of the previous loop
      if (!first_time)
        ++sib;
      for (; sib != poses.end(node); ++sib) { // Iterate over children of node
        vector<vector3> tcoords = (*sib).first;
        vector<vector3> tcoords_hvy = GetHeavyAtomCoords(tcoords);
        palign->SetTarget(tcoords_hvy);

        palign->Align();
        rmsd = palign->GetRMSD();
        n_rmsd++;
        if (rmsd < levels.at(level)) {
          if (rmsd < cutoff)
            return false;

          min_nodes.push_back(sib);
          min_nodes_rmsds.push_back(rmsd);

          if (!_percise) // Exit as soon as one is found
            break;

        }
      } // end of for loop

      if (min_nodes.size() == 0) {
        // No similar molecule found, so remember it for later so that we can
        // append it the children and add it as the first child all the way down
        // through the levels. The reason we don't add it now is that the molecule
        // could still be rejected for addition to the tree.
        insert_pt.push_back(node);
        insert_level.push_back(level);
        insert_data.push_back(PosePair(vcoords, energy));
        continue;
      }

      // If we reach here, then a similar molecule was found
      int startlevel = level;
      vector<double>::const_iterator n_rmsd = min_nodes_rmsds.begin();
      for (vector<Tree_it>::iterator n = min_nodes.begin(); n != min_nodes.end(); ++n, ++n_rmsd) {
        node = *n;
        rmsd = *n_rmsd;
        level = startlevel + 1;
        while (rmsd < levels.at(level)) {
          node = poses.child(node, 0); // Get the first child
          level++;
        }
        nodes.push_back(node);
        stack_levels.push_back(level);
      }

      first_time = false;

    } // end of while loop


    // If we get here, then the molecule has been accepted for addition to the tree
    vector<PosePair>::iterator b = insert_data.begin();
    vector<int>::iterator c = insert_level.begin();
    for (vector<Tree_it>::iterator a = insert_pt.begin(); a != insert_pt.end(); ++a, ++b, ++c) {
      node = *a;
      for (int k = *c; k < levels.size(); ++k) {
        node = poses.append_child(node, *b);
      }
    }

    return true;
  }

  size_t OBDiversePoses::GetSize() {
    return poses.size() - 1; // Remove the dummy
  }

  vector<vector3> OBDiversePoses::GetHeavyAtomCoords(const vector<vector3> &all_coords) {
    vector<vector3> v_hvyatoms;
    for (unsigned int a = 0; a < natoms; ++a)
      if (!hydrogens.BitIsSet(a))
        v_hvyatoms.push_back(all_coords[a]);
    return v_hvyatoms;
  }

  //bool sortpred(const OBDiversePoses::PosePair *a, const OBDiversePoses::PosePair *b) {
  //  return (a->second < b->second);
  //}
  bool sortpred_b(const OBDiversePoses::PosePair& a, const OBDiversePoses::PosePair& b) {
    return (a.second < b.second);
  }

  int OBForceField::FastRotorSearch(bool permute)
  {
    if (_mol.NumRotors() == 0)
      return 0;

    int origLogLevel = _loglvl;

    // Remove all conformers (e.g. from previous conformer generators) except for current conformer
    double *initialCoord = new double [_mol.NumAtoms() * 3]; // initial state
    double *store_initial = new double [_mol.NumAtoms() * 3]; // store the initial state
    memcpy((char*)initialCoord,(char*)_mol.GetCoordinates(),sizeof(double)*3*_mol.NumAtoms());
    memcpy((char*)store_initial,(char*)_mol.GetCoordinates(),sizeof(double)*3*_mol.NumAtoms());
    std::vector<double *> newConfs(1, initialCoord);
    _mol.SetConformers(newConfs);

    _energies.clear(); // Wipe any energies from previous conformer generators

    OBRotorList rl;
    OBBitVec fixed = _constraints.GetFixedBitVec();
    rl.SetFixAtoms(fixed);
    rl.SetQuiet();
    rl.Setup(_mol);

    OBRotorIterator ri;
    OBRotamerList rotamerlist;
    rotamerlist.SetBaseCoordinateSets(_mol);
    rotamerlist.Setup(_mol, rl);

    // Start with all of the rotors in their 0 position
    // (perhaps instead I should set them randomly?)
    std::vector<int> init_rotorKey(rl.Size() + 1, 0);
    std::vector<int> rotorKey(init_rotorKey);

    unsigned int j, minj;
    double currentE, minE, best_minE;

    double *verybestconf = new double [_mol.NumAtoms() * 3]; // store the best conformer to date
    double *bestconf = new double [_mol.NumAtoms() * 3]; // store the best conformer to date in the current permutation
    double *minconf = new double [_mol.NumAtoms() * 3];  // store the best conformer for the current rotor
    memcpy((char*)bestconf,(char*)_mol.GetCoordinates(),sizeof(double)*3*_mol.NumAtoms());

    double energy_offset;
    // Can take shortcut later, as 4 components of the energy will be constant
    rotamerlist.SetCurrentCoordinates(_mol, rotorKey);
    SetupPointers();
    energy_offset = E_Bond(false) + E_Angle(false) + E_StrBnd(false) + E_OOP(false);

    // This function relies on the fact that Rotors are ordered from the most
    // central to the most peripheral (due to CompareRotors in rotor.cpp)
    std::vector<OBRotor *> vrotors;
    OBRotor *rotor;
    for (rotor = rl.BeginRotor(ri); rotor; rotor = rl.NextRotor(ri))
      vrotors.push_back(rotor);

    // The permutations are ordered so that the first 2 permutations cover the
    // combinations of 2, and the first 6 permutations cover the combinations of 3
    const char permutations[24*4] = {0,1,2,3, 1,0,2,3, 0,2,1,3, 1,2,0,3, 2,0,1,3, 2,1,0,3,
                                     0,1,3,2, 0,2,3,1, 0,3,1,2, 0,3,2,1, 1,0,3,2, 1,2,3,0,
                                     1,3,0,2, 1,3,2,0, 2,0,3,1, 2,1,3,0, 2,3,0,1, 2,3,1,0,
                                     3,0,1,2, 3,0,2,1, 3,1,0,2, 3,1,2,0, 3,2,0,1, 3,2,1,0};
    const char factorial[5] = {0, 1, 2, 6, 24};

    char num_rotors_to_permute, num_permutations;
    if (permute)
      num_rotors_to_permute = std::min<size_t> (4, vrotors.size());
    else
      num_rotors_to_permute = 1; // i.e. just use the original order
    num_permutations = factorial[num_rotors_to_permute];

    // Initialize reordered_rotors - the order in which to test rotors
    std::vector<unsigned int> reordered_rotors(vrotors.size());
    for (int i=0; i<vrotors.size(); ++i)
      reordered_rotors[i] = i;

    std::set<unsigned int> seen;
    best_minE = DBL_MAX;
    for (int N=0; N<num_permutations; ++N) {
      for (int i=0; i<num_rotors_to_permute; ++i)
        reordered_rotors.at(i) = *(permutations + N*4 + i);

      rotorKey = init_rotorKey;
      _mol.SetCoordinates(store_initial);
      bool quit = false;

      for (int i=0; i<reordered_rotors.size(); ++i) {
        unsigned int idx = reordered_rotors[i];
        rotor = vrotors.at(idx);

        minE = DBL_MAX;

        for (j = 0; j < rotor->GetResolution().size(); j++) { // For each rotor position
          // Note: we could do slightly better by skipping the rotor position we already
          //       tested in the last loop (position 0 at the moment). Note that this
          //       isn't as simple as just changing the loop starting point to j = 1.
          _mol.SetCoordinates(bestconf);
          rotorKey[idx + 1] = j;
          rotamerlist.SetCurrentCoordinates(_mol, rotorKey);
          SetupPointers();

          currentE = E_VDW(false) + E_Torsion(false) + E_Electrostatic(false);

          if (currentE < minE) {
            minE = currentE;
            minj = j;
            memcpy((char*)minconf,(char*)_mol.GetCoordinates(),sizeof(double)*3*_mol.NumAtoms());
          }
        } // Finished testing all positions of this rotor
        rotorKey[idx + 1] = minj;

        if (i==4) { // Check whether this rotorKey has already been chosen
          // Create a hash of the rotorKeys (given that the max value of any rotorKey is 11 from torlib.txt)
          unsigned int hash = rotorKey[1] + rotorKey[2]*12 + rotorKey[3]*12*12 + rotorKey[4]*12*12*12;

          if (seen.find(hash) == seen.end()) // Not seen before
            seen.insert(hash);
          else { // Already seen - no point continuing
            quit = true;
            break;
          }
        }

        memcpy((char*)bestconf,(char*)minconf,sizeof(double)*3*_mol.NumAtoms());
      } // end of this permutation
      if (!quit) {
          if (minE < best_minE) {
            best_minE = minE;
            memcpy((char*)verybestconf,(char*)bestconf,sizeof(double)*3*_mol.NumAtoms());
          }
      }

    } // end of final permutation

    _mol.SetCoordinates(verybestconf);
    SetupPointers();

    delete [] store_initial;
    delete [] bestconf;
    delete [] verybestconf;
    delete [] minconf;

    return true;
  }

vector<vector3> GetHeavyAtomCoords(const OBMol* mol, const vector<vector3> &all_coords) {
  vector<vector3> v_hvyatoms;
  for (unsigned int a = 1; a <= mol->NumAtoms(); ++a)
    if (!mol->GetAtom(a)->IsHydrogen())
      v_hvyatoms.push_back(all_coords[a]);
  return v_hvyatoms;
}

void UpdateConformersFromTree(OBMol* mol, vector<double> &energies, OBDiversePoses* divposes) {

  OBDiversePoses::Tree* poses = divposes->GetTree();
  double cutoff = divposes->GetCutoff();

  vector <OBDiversePoses::PosePair> confs, newconfs;

  // The leaf iterator will (in effect) iterate over the nodes just at the loweset level
  for (OBDiversePoses::Tree::leaf_iterator node = poses->begin(); node != poses->end(); ++node)
    if (node->first.size() > 0) // Don't include the dummy head node
      confs.push_back(*node);

  // Sort the confs by energy (lowest first)
  sort(confs.begin(), confs.end(), sortpred_b);

  cout << "..(tree size = " << divposes->GetSize() <<  " confs = " << confs.size() << ")\n";

  typedef vector<OBDiversePoses::PosePair> vpp;

  // Loop through the confs and filter using a tree
  newconfs.clear();
  OBDiversePoses newtree(*mol, cutoff, true);
  for (vpp::iterator conf = confs.begin(); conf!=confs.end(); ++conf) {
    if (newtree.AddPose(conf->first, conf->second)) {
      newconfs.push_back(*conf);
    }
  }

  cout << "..(new tree size = " << newtree.GetSize() <<  " confs = " << newconfs.size() << ")\n";

  // Add confs to the molecule's conformer data and add the energies to molecules's energies
  for (vpp::iterator chosen = newconfs.begin(); chosen!=newconfs.end(); ++chosen) {
    energies.push_back(chosen->second);

    // To avoid making copies of vectors or vector3s, I am using pointers throughout
    vector<vector3> *tmp = &(chosen->first);
    double *confCoord = new double [mol->NumAtoms() * 3];
    for(unsigned int a = 0; a<mol->NumAtoms(); ++a) {
      vector3* pv3 = &(*tmp)[a];
      confCoord[a*3] = pv3->x();
      confCoord[a*3 + 1] = pv3->y();
      confCoord[a*3 + 2] = pv3->z();
    }
    mol->AddConformer(confCoord);
  }
}

int OBForceField::DiverseConfGen(double rmsd, unsigned int nconfs, double energy_gap)
  {
    _energies.clear(); // Wipe any energies from previous conformer generators

    // Remove all conformers (e.g. from previous conformer generators) even the current conformer
    double *initialCoord = new double [_mol.NumAtoms() * 3]; // initial state
    double *store_initial = new double [_mol.NumAtoms() * 3]; // store the initial state
    memcpy((char*)initialCoord,(char*)_mol.GetCoordinates(),sizeof(double)*3*_mol.NumAtoms());
    memcpy((char*)store_initial,(char*)_mol.GetCoordinates(),sizeof(double)*3*_mol.NumAtoms());
    std::vector<double *> newConfs(1, initialCoord);
    _mol.SetConformers(newConfs);

    if (_mol.NumRotors() == 0) {
      SetupPointers();
      _energies.push_back(Energy(false));
      delete [] store_initial;
      return 0;
    }

    // Get estimate of lowest energy conf using FastRotorSearch
    FastRotorSearch(true);
    double lowest_energy = Energy(false);

    int origLogLevel = _loglvl;

    OBRotorList rl;
    OBBitVec fixed = _constraints.GetFixedBitVec();
    rl.SetFixAtoms(fixed);
    if (_loglvl == 0)
      rl.SetQuiet(); // Don't print info on symmetry removal
    rl.Setup(_mol);

    OBRotorIterator ri;
    OBRotamerList rotamerlist;
    rotamerlist.SetBaseCoordinateSets(_mol);
    rotamerlist.Setup(_mol, rl);

    // Can take shortcut later, as 4 components of the energy will be constant
    SetupPointers();
    double energy_offset = E_Bond(false) + E_Angle(false) + E_StrBnd(false) + E_OOP(false);
    lowest_energy -= energy_offset;
    _energies.push_back(lowest_energy);

    OBRotorKeys rotorKeys;
    OBRotor* rotor = rl.BeginRotor(ri);
    unsigned int combinations = 1;
    vector<size_t> rotor_sizes;
    for (int i = 1; i < rl.Size() + 1; ++i, rotor = rl.NextRotor(ri)) { // foreach rotor
      size_t size = rotor->GetResolution().size();
      rotorKeys.AddRotor(size);
      combinations *= size;
      rotor_sizes.push_back(size);
      IF_OBFF_LOGLVL_LOW {
        stringstream ss;
        ss << "....rotor " << i << " from " << rotor->GetBond()->GetBeginAtomIdx() << " to ";
        ss << rotor->GetBond()->GetEndAtomIdx() << " has " << size << " values" << endl;
        OBFFLog(ss.str());
      }
    }
    if (rotor_sizes.size() > 0 && combinations == 0) { // Overflow!
      combinations = UINT_MAX;
    }
    cout << "..tot conformations = " << combinations << "\n";

    if (nconfs == 0)
      nconfs = 1 << 20;
    unsigned int max_combinations = min<unsigned int>(nconfs , combinations);
    LFSR lfsr(max_combinations); // Systematic random number generator
    if (combinations > max_combinations) {
      ostringstream ss;
      ss << "There are " << combinations << " conformers. Using a cutoff of "
        << nconfs << " we will only explore " << std::fixed << setprecision(1)
        << static_cast<float>(nconfs * 100)/static_cast<float>(combinations) << "% of these.";
      obErrorLog.ThrowError(__FUNCTION__, ss.str(), obInfo);
    }

    unsigned int combination;
    OBDiversePoses divposes(_mol, rmsd, false);
    vector<int> my_rotorkey(rotor_sizes.size() + 1, 0);
    int counter = 0;

    // Main loop over rotamers
    unsigned int N_low_energy = 0;
    do {
      _mol.SetCoordinates(store_initial);

      combination = lfsr.GetNext();
      unsigned int t = combination;
      // Convert the combination number into a rotorkey
      for (int i = 0 ; i < rotor_sizes.size(); ++i) {
        my_rotorkey[i + 1] = t % rotor_sizes[i];
        t /= rotor_sizes[i];
      }

      rotamerlist.SetCurrentCoordinates(_mol, my_rotorkey);
      SetupPointers();
      double currentE = E_VDW(false) + E_Torsion(false) + E_Electrostatic(false);
      if (currentE < lowest_energy + energy_gap) { // Don't retain high energy poses
        divposes.AddPose(_mol.GetCoordinates(), currentE);
        N_low_energy++;
        if (currentE < lowest_energy)
          lowest_energy = currentE;
      }
      counter++;
    } while (combination != 1 && counter < nconfs); // The LFSR always terminates with a 1
    cout << "..tot confs tested = " << counter << "\n..below energy threshold = " << N_low_energy << "\n";

    // Reset the coordinates to those of the initial structure
    _mol.SetCoordinates(store_initial);

    // Get results from the tree
    UpdateConformersFromTree(&_mol, _energies, &divposes);

    // Add back the energy offset
    transform(_energies.begin(), _energies.end(), _energies.begin(), bind2nd(plus<double>(), energy_offset));

    // Clean up
    delete [] store_initial;

    return 0;
 }

} // end of namespace OpenBabel

//! \file confsearch.cpp
//! \brief Conformer searching routines

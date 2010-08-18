#include <openbabel/isomorphism.h>
#include <openbabel/query.h>
#include <openbabel/graphsym.h>
#include <cassert>

#define DEBUG 0
//#define DEBUG_PERFORMANCE

using namespace std;

namespace OpenBabel {

  static const char *red    = "\033[1;31m";
  static const char *green  = "\033[1;32m";
  static const char *yellow = "\033[1;33m";
  static const char *blue   = "\033[1;34m";
  static const char *normal = "\033[0m";


  class VF2Mapper : public OBIsomorphismMapper
  {
#ifdef DEBUG_PERFORMANCE
      unsigned int numMapNextCalls;
#endif

    public:
      VF2Mapper(OBQuery *query) : OBIsomorphismMapper(query)
      {
      }

      enum MapperType {
        MapFirstType, MapUniqueType, MapAllType
      };

      struct Candidate {
        Candidate() : queryAtom(0), queriedAtom(0) {}
        Candidate(OBQueryAtom *_queryAtom, OBAtom *_queriedAtom)
            : queryAtom(_queryAtom), queriedAtom(_queriedAtom) {}

        bool operator==(const Candidate &other)
        {
          if (queryAtom != other.queryAtom)
            return false;
          if (queriedAtom != other.queriedAtom)
            return false;
          return true;
        }

        OBQueryAtom *queryAtom;
        OBAtom *queriedAtom;
      };

      struct State {
        State(MapperType _type, const OBQuery *_query, const OBMol *_queried, const OBBitVec &mask)
        {
          type = _type;
          query = _query;
          queried = _queried;
          queriedMask = mask;

          queryTerminalSet.resize(query->NumAtoms(), 0);
          queriedTerminalSet.resize(queried->NumAtoms(), 0);
        }
        MapperType type; // MapFirst, MapUnique, MapAll
        const OBQuery *query; // the query
        const OBMol *queried; // the queried molecule
        OBBitVec queriedMask; // the queriedMask
        std::vector<unsigned int> queryPath; // the path in the query
        std::vector<unsigned int> queriedPath; // the path in the queried molecule

        OBBitVec queryPathBits, queriedPathBits; // the terminal sets
        std::vector<unsigned int> queryTerminalSet, queriedTerminalSet; // the terminal sets
      };

      struct KeyValuePair
      {
        KeyValuePair(unsigned int _key, unsigned int _value) : key(_key), value(_value) {}
        unsigned int key, value;
      };
      struct KeyValuePairCompare
      {
        bool operator()(const KeyValuePair &p1, const KeyValuePair &p2) { return p1.key < p2.key; }
      };

      /**
       * Sort mapping by query atom index. Used during MapAll to detect duplicates.
       */
      void SortMapping(Mapping &map)
      {
        // sort the mapping
        std::vector<KeyValuePair> pairs;
        for (Mapping::iterator it = map.begin(); it != map.end(); ++it)
          pairs.push_back(KeyValuePair(it->first, it->second));
        std::sort(pairs.begin(), pairs.end(), KeyValuePairCompare());

        // create a new map and add the pairs in order
        Mapping newMap;
        for (unsigned int i = 0; i < pairs.size(); ++i)
          newMap[pairs[i].key] = pairs[i].value;

        // copy the sorted map
        map = newMap;
      }

      /**
       * Check all bonds between mapped atoms.
       */
      bool checkBonds(State &state, Mapping &map)
      {
        for (unsigned int i = 0; i < state.query->NumBonds(); ++i) {
          OBQueryBond *qbond = state.query->GetBonds()[i];
          unsigned int beginIndex = qbond->GetBeginAtom()->GetIndex();
          unsigned int endIndex = qbond->GetEndAtom()->GetIndex();

          Mapping::iterator beginIt = map.find(beginIndex);
          Mapping::iterator endIt = map.find(endIndex);

          if (beginIt != map.end() && endIt != map.end()) {
            OBAtom *begin = state.queried->GetAtom(beginIt->second + 1);
            OBAtom *end = state.queried->GetAtom(endIt->second + 1);
            if (!begin || !end)
              return false;
            OBBond *bond = state.queried->GetBond(begin, end);
            if (!bond)
              return false;
            if (!qbond->Matches(bond))
              return false;
          }
        }
        return true;
      }

      /**
       * Check if the current state is a full mapping of the query.
       */
      bool checkForMap(State &state, Mappings &maps)
      {
        // store the mapping if all atoms are mapped
        if (state.queryPath.size() != state.query->NumAtoms())
          return false;

        // create the map
        Mapping map;
        for (unsigned int k = 0; k < state.queryPath.size(); ++k) {
          map[state.queryPath[k]] = state.queriedPath[k];
        }

        if (!checkBonds(state, map))
          return true;

        // Check if the mapping is unique
        if (state.type == MapUniqueType) {
          // get the values from the map
          std::vector<unsigned int> values;
          for (Mapping::iterator it = map.begin(); it != map.end(); ++it)
            values.push_back(it->second);
          std::sort(values.begin(), values.end());

          bool isUnique = true;
          for (unsigned int k = 0; k < maps.size(); ++k) {
            std::vector<unsigned int> kValues;
            for (Mapping::iterator it = maps[k].begin(); it != maps[k].end(); ++it)
              kValues.push_back(it->second);
            std::sort(kValues.begin(), kValues.end());

            if (values == kValues)
              isUnique = false;
          }
          if (isUnique) {
            maps.push_back(map);
          }
        } else if (state.type == MapAllType) {
          SortMapping(map);
          bool duplicate = false;
          for (unsigned int i = 0; i < maps.size(); ++i) {
            if (map == maps[i]) {
              duplicate = true;
              break;
            }
          }
          if (!duplicate) {
            if (DEBUG)
              cout << green << "found mapping" << normal << endl;
            maps.push_back(map);
          }
        } else {
          if (DEBUG)
            cout << green << "found mapping" << normal << endl;
          maps.push_back(map);
        }

        return true;
      }

      /**
       * Match the candidate atoms and bonds.
       */
      bool matchCandidate(State &state, OBQueryAtom *queryAtom, OBAtom *queriedAtom, Mappings &maps)
      {
        // make sure the neighbor atom isn't in the paths already
//        if (std::find(state.queryPath.begin(), state.queryPath.end(), queryAtom->GetIndex()) != state.queryPath.end()) {
        if (state.queryPathBits.BitIsSet(queryAtom->GetIndex())) {
          if (DEBUG)
            cout << "    query index already mapped..." << endl;
          return false;
        }
        if (state.queriedPathBits.BitIsSet(queriedAtom->GetIndex())) {
//        if (std::find(state.queriedPath.begin(), state.queriedPath.end(), queriedAtom->GetIndex()) != state.queriedPath.end()) {
          if (DEBUG)
            cout << "    queried index already mapped..." << endl;
          return false;
        }

        // check if the atoms match
        if (!queryAtom->Matches(queriedAtom)) {
          if (DEBUG)
            cout << "    atoms do not match..." << endl;
          return false;
        }

        // check if the bonds match
        Mapping map;
        for (unsigned int k = 0; k < state.queryPath.size(); ++k) {
          map[state.queryPath[k]] = state.queriedPath[k];
        }
        map[queryAtom->GetIndex()] = queriedAtom->GetIndex();
        if (!checkBonds(state, map)) {
          if (DEBUG)
            cout << "    bonds do not match..." << endl;
          return false;
        }

        // add the neighbors to the paths
        state.queryPath.push_back(queryAtom->GetIndex());
        state.queriedPath.push_back(queriedAtom->GetIndex());
        state.queryPathBits.SetBitOn(queryAtom->GetIndex());
        state.queriedPathBits.SetBitOn(queriedAtom->GetIndex());

        //
        // update queryTerminalSet
        //
        if (!state.queryTerminalSet[queryAtom->GetIndex()])
          state.queryTerminalSet[queryAtom->GetIndex()] = state.queryPath.size();

        std::vector<OBQueryAtom*> queryNbrs = queryAtom->GetNbrs();
        for (unsigned int i = 0; i < queryNbrs.size(); ++i) {
          unsigned int index = queryNbrs[i]->GetIndex();
          if (!state.queryTerminalSet[index])
            state.queryTerminalSet[index] = state.queryPath.size();
        }

        //
        // update queriedTerminalSet
        //
        if (!state.queriedTerminalSet[queriedAtom->GetIndex()])
          state.queriedTerminalSet[queriedAtom->GetIndex()] = state.queryPath.size();

        FOR_NBORS_OF_ATOM (nbr, queriedAtom) {
          unsigned int index = nbr->GetIndex();
          // skip atoms not in the mask
          if (!state.queriedMask.BitIsSet(index + 1))
            continue;
          if (!state.queriedTerminalSet[index])
            state.queriedTerminalSet[index] = state.queryPath.size();
        }

        if (DEBUG) {
          cout << "FOUND:  " << queryAtom->GetIndex() << " -> " << queriedAtom->GetIndex() << "       " << state.queryPath.size() << endl;
          cout << "queryTerminalSet:   ";
          for (unsigned int i = 0; i < state.query->NumAtoms(); ++i)
            cout << state.queryTerminalSet[i] << " ";
          cout <<endl;
          cout << "queriedTerminalSet: ";
          for (unsigned int i = 0; i < state.queried->NumAtoms(); ++i)
            cout << state.queriedTerminalSet[i] << " ";
          cout <<endl;
        }

        //
        // Feasibility rules for the VF2 algorithm:
        //
        //  size( T1(s) ) <= size( T2(s) )
        //
        //  size( N1 - M1(s) - T1(s) ) <= size( N2 - M2(s) - T2(s) )
        //

        // compute T1(s) size
        unsigned int numT1 = 0;
        for (unsigned int i = 0; i < state.query->NumAtoms(); ++i)
          if (state.queryTerminalSet[i])
//            if (std::find(state.queryPath.begin(), state.queryPath.end(), i) == state.queryPath.end())
            if (!state.queryPathBits.BitIsSet(i))
              numT1++;
        // compute T2(s) size
        unsigned int numT2 = 0;
        for (unsigned int i = 0; i < state.queried->NumAtoms(); ++i)
          if (state.queriedTerminalSet[i])
            if (!state.queriedPathBits.BitIsSet(i))
//            if (std::find(state.queriedPath.begin(), state.queriedPath.end(), i) == state.queriedPath.end())
              numT2++;

        // T1(s) > T()
        if (numT1 > numT2) {
          Backtrack(state);
          return false;
        }
        //  N1 - M1(s) - T1(s) > N2 - M2(s) - T2(s)
        if ((state.query->NumAtoms() - state.queryPath.size() - numT1) > (state.queried->NumAtoms() - state.queriedPath.size() - numT2)) {
          Backtrack(state);
          return false;
        }

        // Check if there is a mapping found
        checkForMap(state, maps);

        return true;
      }

      /**
       * The depth-first isomorphism algorithm.
       */
      void MapNext(State &state, OBQueryAtom *queryAtom, OBAtom *queriedAtom, Mappings &maps)
      {
#ifdef DEBUG_PERFORMANCE
        numMapNextCalls++;
#endif

        // load the possible candidates
        std::vector<Candidate> candidates;

        //
        // The terminology used in these comments is taken from the VF2 paper.
        //
        //  G1 : The query graph (i.e. state.query)
        //  G2 : The queried graph (i.e. state.queried)
        //
        //  N1 : The set of query atoms in G1
        //  N2 : The set of queried atoms in G2
        //
        //  M1(s) : The set of already mapped query atoms (i.e. state.queryPath)
        //  M2(s) : The set of already mapped queried atoms (i.e. state.queriedPath)
        //


        //  P(s) = T2(s) x { min T1(s) }        [note: this formula is incorrect in the paper!]
        //
        //  P(s) : The set of mapping candidates to consider for state s
        //
        //  T1(s) : The set of "terminal" query atoms for state s
        //  T2(s) : The set of "terminal" queried atoms for state s
        //
        //  min T1(s) : The element from T1(s) with the lowest index. The { } brackets make the element a set again.
        //
        // T1(s) is computed using state.queryTerminalSet and the current mapping (i.e. state.queryPath).
        // A query atom with index n is said to be in the T1(s) set if state.queryTerminalSet[n] is non-zero
        // and n is not part of the current mapping (i.e. n is not in state.queryPath). T2(s) is computed in
        // a similar way.
        OBQueryAtom *queryTerminal = 0;
        for (unsigned int j = 0; j < state.queried->NumAtoms(); ++j) {
          if (state.queriedTerminalSet[j])
            if (!state.queriedPathBits.BitIsSet(j)) {
//            if (std::find(state.queriedPath.begin(), state.queriedPath.end(), j) == state.queriedPath.end()) {
              OBAtom *queriedTerminal = state.queried->GetAtom(j + 1);

              if (!queryTerminal) {
                for (unsigned int i = 0; i < state.query->NumAtoms(); ++i)
                  if (state.queryTerminalSet[i])
                    if (!state.queryPathBits.BitIsSet(i)) {
//                    if (std::find(state.queryPath.begin(), state.queryPath.end(), i) == state.queryPath.end()) {
                      OBQueryAtom *n = state.query->GetAtoms()[i];
                      if (n->Matches(queriedTerminal)) {
                        queryTerminal = n;
                        break;
                      }
                    }
              }

              if (queryTerminal && queryTerminal->Matches(queriedTerminal))
                candidates.push_back(Candidate(queryTerminal, queriedTerminal));
            }
        }

        // If the P(s) set from above is empty, use this set:
        //
        //  P(s) = (N2 - M2(s)) x { min (N1 - M1(s)) }
        if (candidates.empty()) {
          queryTerminal = 0;
          for (unsigned int j = 0; j < state.queried->NumAtoms(); ++j)
//            if (std::find(state.queriedPath.begin(), state.queriedPath.end(), j) == state.queriedPath.end()) {
            if (!state.queriedPathBits.BitIsSet(j)) {
              OBAtom *queriedTerminal = state.queried->GetAtom(j + 1);

              if (!queryTerminal) {
                for (unsigned int i = 0; i < state.query->NumAtoms(); ++i)
//                  if (std::find(state.queryPath.begin(), state.queryPath.end(), i) == state.queryPath.end()) {
                  if (!state.queryPathBits.BitIsSet(i)) {
                    OBQueryAtom *n = state.query->GetAtoms()[i];
                    if (n->Matches(queriedTerminal)) {
                      queryTerminal = n;
                      break;
                    }
                  }
              }

              if (queryTerminal && queryTerminal->Matches(queriedTerminal))
                candidates.push_back(Candidate(queryTerminal, queriedTerminal));
            }
        }

        if (DEBUG) {
          cout << "Candidates:" << endl;
          for (unsigned int i = 0; i < candidates.size(); ++i)
            cout << "        " << candidates[i].queryAtom->GetIndex() << " -> "
                               << candidates[i].queriedAtom->GetIndex() << endl;
        }

        // do the mapping by checking the candidates
        while (candidates.size()) {
          Candidate candidate(candidates.back());
          candidates.pop_back();

          if (DEBUG)
            cout << yellow << "candidate: " << candidate.queryAtom->GetIndex() << " -> " << candidate.queriedAtom->GetIndex() << normal << endl;

          if (matchCandidate(state, candidate.queryAtom, candidate.queriedAtom, maps)) {
            MapNext(state, candidate.queryAtom, candidate.queriedAtom, maps);
          }
        }

        // backtrack
        Backtrack(state);
      }

      void Backtrack(State &state)
      {
        if (DEBUG)
          cout << red << "backtrack... " << normal << state.queryPath.size()-1 << endl;
        // remove last atoms from the mapping
        if (state.queryPath.size()) {
          state.queryPathBits.SetBitOff(state.queryPath.back());
          state.queryPath.pop_back();
        }
        if (state.queriedPath.size()) {
          state.queriedPathBits.SetBitOff(state.queriedPath.back());
          state.queriedPath.pop_back();
        }
        // restore queryTerminalSet and queriedTerminalSet
        unsigned int depth = state.queryPath.size() + 1;
        for (unsigned int i = 0; i < state.query->NumAtoms(); ++i)
          if (state.queryTerminalSet[i] == depth)
            state.queryTerminalSet[i] = 0;
        for (unsigned int i = 0; i < state.queried->NumAtoms(); ++i)
          if (state.queriedTerminalSet[i] == depth)
            state.queriedTerminalSet[i] = 0;
      }

      /**
       * Get the first mappings of the query on the queried molecule.
       * @param queried The queried molecule.
       * @return The mapping.
       */
      Mapping MapFirst(const OBMol *queried, const OBBitVec &mask)
      {
        // set all atoms to 1 if the mask is empty
        OBBitVec queriedMask = mask;
        if (!queriedMask.CountBits())
          for (unsigned int i = 0; i < queried->NumAtoms(); ++i)
            queriedMask.SetBitOn(i + 1);

        Mappings maps;
        OBQueryAtom *queryAtom = m_query->GetAtoms()[0];
        for (unsigned int i = 0; i < queried->NumAtoms(); ++i) {
          if (!queriedMask.BitIsSet(i + 1))
            continue;
          State state(MapFirstType, m_query, queried, queriedMask);

          OBAtom *queriedAtom = queried->GetAtom(i+1);
          if (!queryAtom->Matches(queriedAtom)) {
            continue;
          }

          if (m_query->NumAtoms() > 1) {
            if (matchCandidate(state, queryAtom, queriedAtom, maps))
              MapNext(state, queryAtom, queriedAtom, maps);
          } else {
            Mapping map;
            map[queryAtom->GetIndex()] = queriedAtom->GetIndex();
            maps.push_back(map);
          }

          if (maps.size())
            return maps[0];
        }

        return Mapping();
      }

      /**
       * Get all unique mappings of the query on the queried molecule.
       * @param queried The queried molecule.
       * @return The unique mappings
       */
      Mappings MapUnique(const OBMol *queried, const OBBitVec &mask)
      {
#ifdef DEBUG_PERFORMANCE
        numMapNextCalls = 0;
#endif
        // set all atoms to 1 if the mask is empty
        OBBitVec queriedMask = mask;
        if (!queriedMask.CountBits())
          for (unsigned int i = 0; i < queried->NumAtoms(); ++i)
            queriedMask.SetBitOn(i + 1);

        Mappings maps;
        OBQueryAtom *queryAtom = m_query->GetAtoms()[0];
        for (unsigned int i = 0; i < queried->NumAtoms(); ++i) {
          if (!queriedMask.BitIsSet(i + 1))
            continue;
          State state(MapUniqueType, m_query, queried, queriedMask);

          OBAtom *queriedAtom = queried->GetAtom(i+1);
          if (!queryAtom->Matches(queriedAtom)) {
            continue;
          }

          if (m_query->NumAtoms() > 1) {
            if (matchCandidate(state, queryAtom, queriedAtom, maps))
              MapNext(state, queryAtom, queriedAtom, maps);
          } else {
            Mapping map;
            map[queryAtom->GetIndex()] = queriedAtom->GetIndex();
            maps.push_back(map);
          }
        }

        if (DEBUG)
          for (unsigned int i =0; i < maps.size(); ++i) {
            cout << "mapping:" << endl;
            for (Mapping::iterator it = maps[i].begin(); it != maps[i].end(); ++it)
              cout << "    " << it->first << " -> " << it->second << endl;
          }

#ifdef DEBUG_PERFORMANCE
        cout << "# MapNext calls: " << numMapNextCalls << endl;
#endif

        return maps;
      }

      /**
       * Get all mappings of the query on the queried molecule. Duplicates are
       * ignored but unlinke MapUnique, multiple mappings of the query on the same
       * part of the queried structure are allowed. This makes it possible to use
       * MapAll for finding the automorphism group.
       * @param queried The queried molecule.
       * @return The mappings.
       */
      Mappings MapAll(const OBMol *queried, const OBBitVec &mask)
      {
#ifdef DEBUG_PERFORMANCE
        numMapNextCalls = 0;
#endif
        // set all atoms to 1 if the mask is empty
        OBBitVec queriedMask = mask;
        if (!queriedMask.CountBits())
          for (unsigned int i = 0; i < queried->NumAtoms(); ++i)
            queriedMask.SetBitOn(i + 1);

        Mappings maps;
        OBQueryAtom *queryAtom = m_query->GetAtoms()[0];
        for (unsigned int i = 0; i < queried->NumAtoms(); ++i) {
          if (!queriedMask.BitIsSet(i + 1)) {
            continue;
          }
          State state(MapAllType, m_query, queried, queriedMask);
          OBAtom *queriedAtom = queried->GetAtom(i+1);
          if (!queryAtom->Matches(queriedAtom)) {
            continue;
          }

          if (m_query->NumAtoms() > 1) {
            if (matchCandidate(state, queryAtom, queriedAtom, maps))
              MapNext(state, queryAtom, queriedAtom, maps);
          } else {
            Mapping map;
            map[queryAtom->GetIndex()] = queriedAtom->GetIndex();
            maps.push_back(map);
          }
        }

        if (DEBUG)
          for (unsigned int i =0; i < maps.size(); ++i) {
            cout << "mapping:" << endl;
            for (Mapping::iterator it = maps[i].begin(); it != maps[i].end(); ++it)
              cout << "    " << it->first << " -> " << it->second << endl;
          }

#ifdef DEBUG_PERFORMANCE
        cout << "# MapNext calls: " << numMapNextCalls << endl;
#endif

        return maps;
      }

  };

  OBIsomorphismMapper* OBIsomorphismMapper::GetInstance(OBQuery *query, const std::string &algorithm)
  {
    if (algorithm == "VF2")
      return new VF2Mapper(query);
    // return VF2 mapper as default
    return new VF2Mapper(query);
  }

  class OBAutomorphismQueryAtom : public OBQueryAtom
  {
    public:
      OBAutomorphismQueryAtom(unsigned int _symClass, const std::vector<unsigned int> &_symClasses)
          : OBQueryAtom(), symClass(_symClass), symClasses(_symClasses)
      {
      }

      bool Matches(const OBAtom *atom) const
      {
        return (symClasses[atom->GetIndex()] == symClass);
      }
      unsigned int symClass;
      std::vector<unsigned int> symClasses;
  };

  OBQuery* CompileAutomorphismQuery(OBMol *mol, const OBBitVec &mask, const std::vector<unsigned int> &symClasses)
  {
    OBQuery *query = new OBQuery;
    unsigned int offset = 0;
    std::vector<unsigned int> indexes;
    FOR_ATOMS_OF_MOL (obatom, mol) {
      indexes.push_back(obatom->GetIndex() - offset);
      if (!mask.BitIsSet(obatom->GetIndex() + 1)) {
        offset++;
        continue;
      }
      query->AddAtom(new OBAutomorphismQueryAtom(symClasses[obatom->GetIndex()], symClasses));
    }
    FOR_BONDS_OF_MOL (obbond, mol) {
      unsigned int beginIndex = obbond->GetBeginAtom()->GetIndex();
      unsigned int endIndex = obbond->GetEndAtom()->GetIndex();
      if (!mask.BitIsSet(beginIndex + 1) || !mask.BitIsSet(endIndex + 1))
        continue;
      query->AddBond(new OBQueryBond(query->GetAtoms()[indexes[beginIndex]], query->GetAtoms()[indexes[endIndex]],
            obbond->GetBondOrder(), obbond->IsAromatic()));
    }

    return query;
  }

  OBIsomorphismMapper::Mappings FindAutomorphisms(OBMol *mol, const OBBitVec &mask)
  {
    // set all atoms to 1 if the mask is empty
    OBBitVec queriedMask = mask;
    if (!queriedMask.CountBits())
      for (unsigned int i = 0; i < mol->NumAtoms(); ++i)
        queriedMask.SetBitOn(i + 1);

    // get the symmetry classes
    OBGraphSym gs(mol, &queriedMask);
    std::vector<unsigned int> symClasses;
    gs.GetSymmetry(symClasses);

    if (DEBUG)
    for (unsigned int i = 0; i < symClasses.size(); ++i)
      cout << i << ": " << symClasses[i] << endl;

    // count the symmetry classes
    std::vector<int> symClassCounts(symClasses.size() + 1, 0);
    for (unsigned int i = 0; i < symClasses.size(); ++i) {
      if (!queriedMask.BitIsSet(i + 1))
        continue;
      unsigned int symClass = symClasses[i];
      symClassCounts[symClass]++;
    }
    // construct the premap
    OBIsomorphismMapper::Mapping premap;
    unsigned int offset = 0;
    for (unsigned int i = 0; i < symClasses.size(); ++i) {
      if (!queriedMask.BitIsSet(i + 1)) {
        offset++;
        continue;
      }

      if (symClassCounts[symClasses[i]] == 1)
        premap[i-offset] = i;
    }

    OBIsomorphismMapper::Mappings maps;
    if (premap.size() < queriedMask.CountBits()) {
      OBQuery *query = CompileAutomorphismQuery(mol, queriedMask, symClasses);
      OBIsomorphismMapper *mapper = OBIsomorphismMapper::GetInstance(query);

      maps = mapper->MapAll(mol, queriedMask);
      delete mapper;
      delete query;
    } else {
      maps.push_back(premap);
    }

    // translate the query indexes if a mask is used
    if (mask.CountBits()) {
      std::vector<unsigned int> indexes;
      unsigned int offset = 0;
      for (unsigned int i = 0; i < mol->NumAtoms(); ++i)
        if (queriedMask.BitIsSet(i + 1))
          indexes.push_back(i);

      OBIsomorphismMapper::Mappings translatedMaps;
      for (OBIsomorphismMapper::Mappings::iterator map = maps.begin(); map != maps.end(); map++) {
        OBIsomorphismMapper::Mapping m;
        for (OBIsomorphismMapper::Mapping::iterator i = map->begin(); i != map->end(); i++)
          m[indexes[i->first]] = i->second;
        translatedMaps.push_back(m);
      }

      return translatedMaps;
    }

    return maps;
  }




}

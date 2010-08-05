#include <openbabel/isomorphism.h>
#include <openbabel/query.h>
#include <openbabel/graphsym.h>
#include <cassert>

#define DEBUG 0

using namespace std;

namespace OpenBabel {
 

  class BFSMapper : public OBIsomorphismMapper
  {
    public:
      BFSMapper(OBQuery *query) : OBIsomorphismMapper(query) 
      {
      }

      enum MapperType {
        MapFirstType, MapUniqueType, MapAllType
      };

      struct Candidate {
        Candidate(OBQueryAtom *_queryAtom, OBQueryAtom *_queryNbr, OBAtom *_queriedAtom, OBAtom *_queriedNbr)
            : queryAtom(_queryAtom), queryNbr(_queryNbr), queriedAtom(_queriedAtom), queriedNbr(_queriedNbr) {}
        OBQueryAtom *queryAtom, *queryNbr;
        OBAtom *queriedAtom, *queriedNbr;
      };

      struct State {
        State(MapperType _type, const OBQuery *_query, const OBMol *_queried, const OBBitVec &mask)
        {
          type = _type;
          query = _query;
          queried = _queried;
          queriedMask = mask;
        }
        MapperType type; // MapFirst, MapUnique, MapAll
        const OBQuery *query; // the query
        const OBMol *queried; // the queried molecule
        OBBitVec queriedMask; // the queriedMask
        std::vector<unsigned int> queryPath; // the path in the query
        std::vector<unsigned int> queriedPath; // the path in the queried molecule
        std::vector<std::vector<Candidate> > candidates; // the candidates to check
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
       * Sort mapping by query atom index
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
            if (!state.queried->GetBond(begin, end))
              return false;
          }
        }
        return true;
      }

      /**
       * Check if the current state is a full mapping of the query.
       */
      void checkForMap(State &state, Mappings &maps) 
      {

        // store the mapping if all atoms are mapped
        if (state.queryPath.size() == state.query->NumAtoms()) {
          // create the map
          Mapping map;
          for (unsigned int k = 0; k < state.queryPath.size(); ++k) {
            map[state.queryPath[k]] = state.queriedPath[k];
          }

          if (!checkBonds(state, map))
            return;

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
                cout << "found mapping" << endl;
              maps.push_back(map);
            }
          } else {
            if (DEBUG)
              cout << "found mapping" << endl;
            maps.push_back(map);
          }
        }
      }

      /**
       * Match the candidate atoms and bonds.
       */
      bool matchCandidate(State &state, OBQueryAtom *queryAtom, OBAtom *queriedAtom, 
          OBQueryAtom *queryNbr, OBAtom *queriedNbr, Mappings &maps) 
      {
        // make sure the neighbor atom isn't in the paths already
        if (std::find(state.queryPath.begin(), state.queryPath.end(), queryNbr->GetIndex()) != state.queryPath.end())
          return false;
        if (std::find(state.queriedPath.begin(), state.queriedPath.end(), queriedNbr->GetIndex()) != state.queriedPath.end())
          return false;

        // check if the atoms match
        if (!queryNbr->Matches(queriedNbr))
          return false;

        OBQueryBond *queryBond = state.query->GetBond(queryAtom, queryNbr);
        OBBond *queriedBond = state.queried->GetBond(queriedAtom, queriedNbr);

        // check if the bonds match
        if (!queryBond->Matches(queriedBond))
          return false;

        // create the map
        Mapping map;
        for (unsigned int k = 0; k < state.queryPath.size(); ++k) {
          map[state.queryPath[k]] = state.queriedPath[k];
        }
        map[queryNbr->GetIndex()] = queriedNbr->GetIndex();
        if (!checkBonds(state, map))
          return false;


        // add the neighbors to the paths
        state.queryPath.push_back(queryNbr->GetIndex());
        state.queriedPath.push_back(queriedNbr->GetIndex());

        if (DEBUG)
          cout << "FOUND:  " << queryNbr->GetIndex() << " -> " << queriedNbr->GetIndex() << "       " << state.queryPath.size() << endl; 

        // check if this is a full match
        checkForMap(state, maps);

        return true;
      }

      /**
       * The depth-first isomorphism algorithm.
       */
      void BFS(State &state, OBQueryAtom *queryAtom, OBAtom *queriedAtom, Mappings &maps) {
        std::vector<OBQueryAtom*> queryNbrs = queryAtom->GetNbrs();
        std::vector<OBAtom*> queriedNbrs;
        FOR_NBORS_OF_ATOM (nbr, queriedAtom)
          queriedNbrs.push_back(&*nbr);

        // load the possible candidates
        std::vector<Candidate> candidates;
        for (unsigned int i = 0; i < queryNbrs.size(); ++i) {
          OBQueryAtom *queryNbr = queryNbrs[i];
          for (unsigned int j = 0; j < queriedNbrs.size(); ++j) {
            OBAtom *queriedNbr = queriedNbrs[j];
            // skip atoms not in the mask
            if (!state.queriedMask.BitIsSet(queriedNbr->GetIndex() + 1))
              continue;
            if (queryNbr->Matches(queriedNbr))
              candidates.push_back(Candidate(queryAtom, queryNbr, queriedAtom, queriedNbr));
          }
        }
        // save the candidates for later (used to explore branches)
        if (state.candidates.size()) {
          std::vector<Candidate> candidate = state.candidates.back();
          for (unsigned int i = 0; i < candidates.size(); ++i)
            candidate.push_back(candidates[i]);
          state.candidates.push_back(candidate);
        } else {
          state.candidates.push_back(candidates);
        }

        if (DEBUG) {
          cout << "Candidates:" << endl;
          for (unsigned int i = 0; i < state.candidates.back().size(); ++i)
            cout << "    " << state.candidates.back()[i].queryNbr->GetIndex() << " -> " 
                           << state.candidates.back()[i].queriedNbr->GetIndex() << endl;
        }

        // do the mapping by checking the candidates
        while (state.candidates.back().size()) {
          Candidate candidate(state.candidates.back().back());
          state.candidates.back().pop_back();
 

          if (matchCandidate(state, candidate.queryAtom, candidate.queriedAtom, candidate.queryNbr, candidate.queriedNbr, maps)) {
            BFS(state, candidate.queryNbr, candidate.queriedNbr, maps);
        
            if (DEBUG)
              cout << "backtrack... " << state.queryPath.size()-1 << endl;
           // backtrack
            if (state.queryPath.size())
              state.queryPath.pop_back();
            if (state.queriedPath.size())
              state.queriedPath.pop_back();
            if (state.candidates.size())
              state.candidates.pop_back();
 
          }


        }
       
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
            state.queryPath.push_back(queryAtom->GetIndex());
            state.queriedPath.push_back(queriedAtom->GetIndex());
            BFS(state, queryAtom, queriedAtom, maps);
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
            state.queryPath.push_back(queryAtom->GetIndex());
            state.queriedPath.push_back(queriedAtom->GetIndex());
            BFS(state, queryAtom, queriedAtom, maps);
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
        // set all atoms to 1 if the mask is empty
        OBBitVec queriedMask = mask;
        if (!queriedMask.CountBits())
          for (unsigned int i = 0; i < queried->NumAtoms(); ++i)
            queriedMask.SetBitOn(i + 1);

        Mappings maps;
        for (unsigned int j = 0; j < m_query->NumAtoms(); ++j) {
          OBQueryAtom *queryAtom = m_query->GetAtoms()[j];
          for (unsigned int i = 0; i < queried->NumAtoms(); ++i) {
            if (!queriedMask.BitIsSet(i + 1))
              continue;
            State state(MapAllType, m_query, queried, queriedMask);
            OBAtom *queriedAtom = queried->GetAtom(i+1);
            if (!queryAtom->Matches(queriedAtom)) {
              continue;
            }

            if (m_query->NumAtoms() > 1) {
              state.queryPath.push_back(queryAtom->GetIndex());
              state.queriedPath.push_back(queriedAtom->GetIndex());
              BFS(state, queryAtom, queriedAtom, maps);
            } else {
              Mapping map;
              map[queryAtom->GetIndex()] = queriedAtom->GetIndex();
              maps.push_back(map);
            }
          }
        }

        if (DEBUG)
          for (unsigned int i =0; i < maps.size(); ++i) {
            cout << "mapping:" << endl;
            for (Mapping::iterator it = maps[i].begin(); it != maps[i].end(); ++it)
              cout << "    " << it->first << " -> " << it->second << endl;
          }

        return maps;
      }
 
  };

  OBIsomorphismMapper* OBIsomorphismMapper::GetInstance(OBQuery *query, const std::string &algorithm)
  {
    if (algorithm == "bfs")
      return new BFSMapper(query);
    return 0;
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

  OBQuery* CompileAutomorphismQuery(OBMol *mol, const OBBitVec &mask)
  {
    OBGraphSym gs(mol, (OBBitVec*)&mask);
    std::vector<unsigned int> symClasses;
    gs.GetSymmetry(symClasses);

    if (DEBUG)
    for (unsigned int i = 0; i < symClasses.size(); ++i)
      cout << i << ": " << symClasses[i] << endl;

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

    OBQuery *query = CompileAutomorphismQuery(mol, queriedMask);
    OBIsomorphismMapper *mapper = OBIsomorphismMapper::GetInstance(query);
    OBIsomorphismMapper::Mappings maps = mapper->MapAll(mol, queriedMask);
    delete mapper;
    delete query;
    return maps;
  }




}

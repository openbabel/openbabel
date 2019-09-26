#include <openbabel/mol.h>
#include <openbabel/isomorphism.h>
#include <openbabel/obiter.h>
#include <openbabel/oberror.h>
#include <openbabel/query.h>
#include <openbabel/graphsym.h>
#include <ctime>
#include <cassert>
#include <algorithm>

#define DEBUG 0

using namespace std;

namespace OpenBabel {

  template<typename T>
    void print_vector(const std::string &label, const std::vector<T> &v)
    {
      std::cout << label << ": ";
      for (std::size_t i = 0; i < v.size(); ++i)
        if (v[i] < 10)
          std::cout << " " << v[i] << " ";
        else
          std::cout << v[i] << " ";

      std::cout << endl;
    }

  static const char *red    = "\033[1;31m";
  static const char *green  = "\033[1;32m";
  static const char *yellow = "\033[1;33m";
  static const char *blue   = "\033[1;34m";
  static const char *normal = "\033[0m";

  class MapAllFunctor : public OBIsomorphismMapper::Functor
  {
    private:
      OBIsomorphismMapper::Mappings &m_maps;
      std::size_t m_memory, m_maxMemory;
    public:
      MapAllFunctor(OBIsomorphismMapper::Mappings &maps, std::size_t maxMemory)
        : m_maps(maps), m_memory(0), m_maxMemory(maxMemory)
      {
      }
      bool operator()(OBIsomorphismMapper::Mapping &map)
      {
        m_maps.push_back(map);
        m_memory += 2 * sizeof(unsigned int) * map.size();

        if (m_memory > m_maxMemory) {
          obErrorLog.ThrowError(__FUNCTION__, "memory limit exceeded...", obError);
          // stop mapping
          return true;
        }

        // continue mapping
        return false;
      }
  };


  class VF2Mapper : public OBIsomorphismMapper
  {
    time_t m_startTime;

    public:
      VF2Mapper(OBQuery *query) : OBIsomorphismMapper(query)
      {
      }

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
        State(Functor &_functor, const OBQuery *_query, const OBMol *_queried, const OBBitVec &mask)
            : functor(_functor), query(_query), queried(_queried), queriedMask(mask)
        {
          abort = false;
          mapping.resize(query->NumAtoms(), 0);
          queryDepths.resize(query->NumAtoms(), 0);
          queriedDepths.resize(queried->NumAtoms(), 0);
        }
        bool abort;
        Functor &functor;
        const OBQuery *query; // the query
        const OBMol *queried; // the queried molecule
        OBBitVec queriedMask; // the queriedMask
        std::vector<unsigned int> queryPath; // the path in the query
        std::vector<unsigned int> queriedPath; // the path in the queried molecule

        std::vector<OBAtom*> mapping;

        OBBitVec queryPathBits, queriedPathBits; // the terminal sets
        std::vector<unsigned int> queryDepths, queriedDepths; // the terminal sets
      };

      bool isInTerminalSet(const std::vector<unsigned int> &depths,
          const OBBitVec &path, std::size_t i)
      {
        if (!depths[i])
          return false;

        if (path.BitIsSet(i))
          return false;

        return true;
      }

      /**
       * Check bonds around newly mapped atom.
       */
      bool checkBonds(State &state, OBQueryAtom *queryAtom)
      {
        const std::vector<OBQueryBond*> &bonds = queryAtom->GetBonds();
        for (unsigned int i = 0; i < bonds.size(); ++i) {
          OBQueryBond *qbond = bonds[i];
          unsigned int beginIndex = qbond->GetBeginAtom()->GetIndex();
          unsigned int endIndex = qbond->GetEndAtom()->GetIndex();

          OBAtom *begin = state.mapping[beginIndex];
          OBAtom *end = state.mapping[endIndex];
          if (!begin || !end)
            continue;
          OBBond *bond = state.queried->GetBond(begin, end);
          if (!bond)
            return false;
          if (!qbond->Matches(bond))
            return false;
        }
        return true;
      }


      /**
       * Check if the current state is a full mapping of the query.
       */
      bool checkForMap(State &state)
      {
        // store the mapping if all atoms are mapped
        if (state.queryPath.size() != state.query->NumAtoms())
          return false;

        if (DEBUG)
          std::cout << green << "-----------------> MATCH" << normal << std::endl;

        // create the map
        Mapping map;
        map.reserve(state.queryPath.size());
        for (unsigned int k = 0; k < state.queryPath.size(); ++k)
          map.push_back(std::make_pair(state.queryPath[k], state.queriedPath[k]));

        return state.functor(map);
      }

      /**
       * Match the candidate atoms and bonds.
       */
      bool matchCandidate(State &state, OBQueryAtom *queryAtom, OBAtom *queriedAtom)
      {
        if (!queryAtom->Matches(queriedAtom))
          return false;

        // add the neighbors to the paths
        state.queryPath.push_back(queryAtom->GetIndex());
        state.queriedPath.push_back(queriedAtom->GetIndex());
        // update the terminal sets
        state.queryPathBits.SetBitOn(queryAtom->GetIndex());
        state.queriedPathBits.SetBitOn(queriedAtom->GetIndex());
        // update mapping
        state.mapping[queryAtom->GetIndex()] = queriedAtom;

        //
        // update queryDepths
        //
        if (!state.queryDepths[queryAtom->GetIndex()])
          state.queryDepths[queryAtom->GetIndex()] = state.queryPath.size();

        std::vector<OBQueryAtom*> queryNbrs = queryAtom->GetNbrs();
        for (unsigned int i = 0; i < queryNbrs.size(); ++i) {
          unsigned int index = queryNbrs[i]->GetIndex();
          if (!state.queryDepths[index])
            state.queryDepths[index] = state.queryPath.size();
        }

        //
        // update queriedDepths
        //
        if (!state.queriedDepths[queriedAtom->GetIndex()])
          state.queriedDepths[queriedAtom->GetIndex()] = state.queriedPath.size();

        FOR_NBORS_OF_ATOM (nbr, queriedAtom) {
          unsigned int index = nbr->GetIndex();
          // skip atoms not in the mask
          if (!state.queriedMask.BitIsSet(index + 1))
            continue;
          if (!state.queriedDepths[index])
            state.queriedDepths[index] = state.queriedPath.size();
        }

        // check if the bonds match
        if (!checkBonds(state, queryAtom)) {
          if (DEBUG)
            cout << "    bonds do not match..." << endl;
          Backtrack(state);
          return false;
        }

        if (DEBUG) {
          cout << "FOUND:  " << queryAtom->GetIndex() << " -> " << queriedAtom->GetIndex() << "       " << state.queryPath.size() << endl;
          cout << "queryDepths:   ";
          for (unsigned int i = 0; i < state.query->NumAtoms(); ++i)
            cout << state.queryDepths[i] << " ";
          cout <<endl;
          cout << "queriedDepths: ";
          for (unsigned int i = 0; i < state.queried->NumAtoms(); ++i)
            cout << state.queriedDepths[i] << " ";
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
          if (isInTerminalSet(state.queryDepths, state.queryPathBits, i))
              numT1++;
        // compute T2(s) size
        unsigned int numT2 = 0;
        for (unsigned int i = 0; i < state.queried->NumAtoms(); ++i)
          if (isInTerminalSet(state.queriedDepths, state.queriedPathBits, i))
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
        state.abort = checkForMap(state);

        return true;
      }


      Candidate NextCandidate(State &state, const Candidate &lastCandidate)
      {
        std::size_t lastQueryAtom = lastCandidate.queryAtom ? lastCandidate.queryAtom->GetIndex() : 0;
        std::size_t lastQueriedAtom = lastCandidate.queriedAtom ? lastCandidate.queriedAtom->GetIndex() + 1 : 0;

        std::size_t querySize = state.query->NumAtoms();
        std::size_t queriedSize = state.queried->NumAtoms();

        std::size_t queryTerminalSize = state.queryDepths.size() - std::count(state.queryDepths.begin(), state.queryDepths.end(), 0);
        std::size_t queriedTerminalSize = state.queriedDepths.size() - std::count(state.queriedDepths.begin(), state.queriedDepths.end(), 0);

        std::size_t mappingSize = state.queryPath.size();

        if (queryTerminalSize > mappingSize && queriedTerminalSize > mappingSize) {
          while (lastQueryAtom < querySize && (state.queryPathBits.BitIsSet(lastQueryAtom) || !state.queryDepths[lastQueryAtom])) {
            lastQueryAtom++;
            lastQueriedAtom = 0;
          }
        } else {
          while(lastQueryAtom < querySize && state.queryPathBits.BitIsSet(lastQueryAtom)) {
            lastQueryAtom++;
            lastQueriedAtom = 0;
          }
        }

        if (queryTerminalSize > mappingSize && queriedTerminalSize > mappingSize) {
          while (lastQueriedAtom < queriedSize && (state.queriedPathBits.BitIsSet(lastQueriedAtom) || !state.queriedDepths[lastQueriedAtom]))
            lastQueriedAtom++;
        } else {
          while(lastQueriedAtom < queriedSize && state.queriedPathBits[lastQueriedAtom])
            lastQueriedAtom++;
        }

        if (lastQueryAtom < querySize && lastQueriedAtom < queriedSize)
          return Candidate(state.query->GetAtoms()[lastQueryAtom], state.queried->GetAtom(lastQueriedAtom + 1));

        return Candidate();
      }

      /**
       * The depth-first isomorphism algorithm.
       */
      void MapNext(State &state, OBQueryAtom *queryAtom, OBAtom *queriedAtom)
      {
        if (time(NULL) - m_startTime > m_timeout)
          return;
        if (state.abort)
          return;

        Candidate candidate;
        while (!state.abort) {
          candidate = NextCandidate(state, candidate);

          if (!candidate.queryAtom)
            return;

          if (DEBUG)
            cout << yellow << "candidate: " << candidate.queryAtom->GetIndex() << " -> " << candidate.queriedAtom->GetIndex() << normal << endl;


          if (matchCandidate(state, candidate.queryAtom, candidate.queriedAtom)) {
            MapNext(state, candidate.queryAtom, candidate.queriedAtom);
            Backtrack(state);
          }
        }

      }

      void Backtrack(State &state)
      {
        if (DEBUG)
          cout << red << "backtrack... " << normal << state.queryPath.size()-1 << endl;
        // remove last atoms from the mapping
        if (state.queryPath.size()) {
          state.mapping[state.queryPath.back()] = 0;
          state.queryPathBits.SetBitOff(state.queryPath.back());
          state.queryPath.pop_back();
        }
        if (state.queriedPath.size()) {
          state.queriedPathBits.SetBitOff(state.queriedPath.back());
          state.queriedPath.pop_back();
        }
        // restore queryDepths and queriedDepths
        unsigned int depth = state.queryPath.size() + 1;
        std::replace(state.queryDepths.begin(), state.queryDepths.end(), depth, static_cast<unsigned int>(0));
        std::replace(state.queriedDepths.begin(), state.queriedDepths.end(), depth, static_cast<unsigned int>(0)); // O(n)  n = # vertices in the queried
      }

      /**
       * Get the first mappings of the query on the queried molecule.
       * @param queried The queried molecule.
       * @return The mapping.
       */
      void MapFirst(const OBMol *queried, Mapping &map, const OBBitVec &mask)
      {
        class MapFirstFunctor : public Functor
        {
          private:
            Mapping &m_map;
          public:
            MapFirstFunctor(Mapping &map) : m_map(map)
            {
            }
            bool operator()(Mapping &map)
            {
              m_map = map;
              // stop mapping
              return true;
            }
        };

        MapFirstFunctor functor(map);
        MapGeneric(functor, queried, mask);
      }

      /**
       * Get all unique mappings of the query on the queried molecule.
       * @param queried The queried molecule.
       * @return The unique mappings
       */
      void MapUnique(const OBMol *queried, Mappings &maps, const OBBitVec &mask)
      {
        class MapUniqueFunctor : public OBIsomorphismMapper::Functor
        {
          private:
            OBIsomorphismMapper::Mappings &m_maps;
          public:
            MapUniqueFunctor(OBIsomorphismMapper::Mappings &maps) : m_maps(maps)
            {
            }
            bool operator()(OBIsomorphismMapper::Mapping &map)
            {
              // get the values from the map
              std::vector<unsigned int> values;
              for (OBIsomorphismMapper::Mapping::const_iterator it = map.begin(); it != map.end(); ++it)
                values.push_back(it->second);
              std::sort(values.begin(), values.end());
               // print_vector("values ", values);

              bool isUnique = true;
              for (unsigned int k = 0; k < m_maps.size(); ++k) {
                std::vector<unsigned int> kValues;
                for (OBIsomorphismMapper::Mapping::iterator it = m_maps[k].begin(); it != m_maps[k].end(); ++it)
                  kValues.push_back(it->second);
                std::sort(kValues.begin(), kValues.end());

              //  print_vector("kValues", kValues);
                if (values == kValues)
                  isUnique = false;
              }

              if (isUnique)
                m_maps.push_back(map);

              // continue mapping
              return false;
            }
        };


        maps.clear();
        MapUniqueFunctor functor(maps);
        MapGeneric(functor, queried, mask);

        if (DEBUG)
          for (unsigned int i =0; i < maps.size(); ++i) {
            cout << "mapping:" << endl;
            for (Mapping::iterator it = maps[i].begin(); it != maps[i].end(); ++it)
              cout << "    " << it->first << " -> " << it->second << endl;
          }
      }

      /**
       * Get all mappings of the query on the queried molecule. Duplicates are
       * ignored but unlinke MapUnique, multiple mappings of the query on the same
       * part of the queried structure are allowed. This makes it possible to use
       * MapAll for finding the automorphism group.
       * @param queried The queried molecule.
       * @return The mappings.
       */
      void MapAll(const OBMol *queried, Mappings &maps, const OBBitVec &mask, std::size_t maxMemory)
      {
        maps.clear();
        MapAllFunctor functor(maps, maxMemory);
        MapGeneric(functor, queried, mask);

        if (DEBUG)
          for (unsigned int i =0; i < maps.size(); ++i) {
            cout << "mapping:" << endl;
            for (Mapping::iterator it = maps[i].begin(); it != maps[i].end(); ++it)
              cout << "    " << it->first << " -> " << it->second << endl;
          }

      }

      void MapGeneric(Functor &functor, const OBMol *queried, const OBBitVec &mask)
      {
        m_startTime = time(NULL);
        if(m_query->NumAtoms() == 0) return;
        // set all atoms to 1 if the mask is empty
        OBBitVec queriedMask = mask;
        if (!queriedMask.CountBits())
          for (unsigned int i = 0; i < queried->NumAtoms(); ++i)
            queriedMask.SetBitOn(i + 1);

        OBQueryAtom *queryAtom = m_query->GetAtoms()[0];
        for (unsigned int i = 0; i < queried->NumAtoms(); ++i) {
          if (!queriedMask.BitIsSet(i + 1)) {
            continue;
          }
          State state(functor, m_query, queried, queriedMask);
          OBAtom *queriedAtom = queried->GetAtom(i+1);
          if (!queryAtom->Matches(queriedAtom)) {
            continue;
          }
          if (DEBUG)
            std::cout << blue << "START: 0 -> " << queriedAtom->GetIndex() << normal << std::endl;

          if (m_query->NumAtoms() > 1) {
            if (matchCandidate(state, queryAtom, queriedAtom))
              MapNext(state, queryAtom, queriedAtom);
          } else {
            Mapping map;
            map.push_back(std::make_pair(queryAtom->GetIndex(), queriedAtom->GetIndex()));
            functor(map);
          }
        }

        if (time(NULL) - m_startTime > m_timeout)
          obErrorLog.ThrowError(__FUNCTION__, "time limit exceeded...", obError);

      }

  };

  OBIsomorphismMapper::OBIsomorphismMapper(OBQuery *query) : m_query(query), m_timeout(60)
  {
  }

  OBIsomorphismMapper::~OBIsomorphismMapper()
  {
  }

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

  bool isFerroceneBond(OBBond *bond);

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
      if (isFerroceneBond(&*obbond))
        continue;
      unsigned int beginIndex = obbond->GetBeginAtom()->GetIndex();
      unsigned int endIndex = obbond->GetEndAtom()->GetIndex();
      if (!mask.BitIsSet(beginIndex + 1) || !mask.BitIsSet(endIndex + 1))
        continue;
      query->AddBond(new OBQueryBond(query->GetAtoms()[indexes[beginIndex]], query->GetAtoms()[indexes[endIndex]],
            obbond->GetBondOrder(), obbond->IsAromatic()));
    }

    return query;
  }

  bool FindAutomorphisms(OBMol *mol, Automorphisms &maps, const OBBitVec &mask, std::size_t maxMemory)
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

    return FindAutomorphisms(mol, maps, symClasses, mask, maxMemory);
  }


  bool FindAutomorphisms(OBMol *mol, Automorphisms &maps, const std::vector<unsigned int> &symClasses, const OBBitVec &mask, std::size_t maxMemory)
  {
    maps.clear();
    MapAllFunctor functor(maps, maxMemory);
    FindAutomorphisms(functor, mol, symClasses, mask);
    return maps.size();
  }

  OBBitVec getFragment(OBAtom *atom, const OBBitVec &mask, const std::vector<OBBond*> &metalloceneBonds = std::vector<OBBond*>());

  void FindAutomorphisms(OBIsomorphismMapper::Functor &functor, OBMol *mol,
      const std::vector<unsigned int> &symClasses, const OBBitVec &mask)
  {
    class AutomorphismFunctor : public OBIsomorphismMapper::Functor
    {
      private:
        OBIsomorphismMapper::Functor &m_functor;
        const OBBitVec &m_fragment;
        std::vector<unsigned int> m_indexes;
      public:
        AutomorphismFunctor(OBIsomorphismMapper::Functor &functor, const OBBitVec &fragment, unsigned int numAtoms)
            : m_functor(functor), m_fragment(fragment)
        {
          for (unsigned int j = 0; j < numAtoms; ++j)
            if (m_fragment.BitIsSet(j + 1))
              m_indexes.push_back(j);
        }
        bool operator()(Automorphism &map)
        {
          // convert the continuous mapping map to a mapping with gaps (considering key values)
          for (Automorphism::iterator it = map.begin(); it != map.end(); ++it)
            it->first = m_indexes[it->first];
          return m_functor(map);
        }
    };

    // set all atoms to 1 if the mask is empty
    OBBitVec queriedMask = mask;
    if (!queriedMask.CountBits())
      for (unsigned int i = 0; i < mol->NumAtoms(); ++i)
        queriedMask.SetBitOn(i + 1);

    if (DEBUG)
    for (unsigned int i = 0; i < symClasses.size(); ++i)
      cout << i << ": " << symClasses[i] << endl;

    // compute the connected fragments
    OBBitVec visited;
    std::vector<OBBitVec> fragments;
    for (std::size_t i = 0; i < mol->NumAtoms(); ++i) {
      if (!queriedMask.BitIsSet(i+1) || visited.BitIsSet(i+1))
        continue;
      fragments.push_back(getFragment(mol->GetAtom(i+1), queriedMask));
      visited |= fragments.back();
    }

    // count the symmetry classes
    std::vector<int> symClassCounts(symClasses.size() + 1, 0);
    for (unsigned int i = 0; i < symClasses.size(); ++i) {
      if (!queriedMask.BitIsSet(i + 1))
        continue;
      unsigned int symClass = symClasses[i];
      symClassCounts[symClass]++;
    }

    for (std::size_t i = 0; i < fragments.size(); ++i) {
      OBQuery *query = CompileAutomorphismQuery(mol, fragments[i], symClasses);
      OBIsomorphismMapper *mapper = OBIsomorphismMapper::GetInstance(query);

      AutomorphismFunctor autFunctor(functor, fragments[i], mol->NumAtoms());
      mapper->MapGeneric(autFunctor, mol, fragments[i]);
      delete mapper;
      delete query;
    }

  }




}

/// @file isomorphism.cpp
/// @brief OBIsomorphismMapper class for finding isomorphisms.

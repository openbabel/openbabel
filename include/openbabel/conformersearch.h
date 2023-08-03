/**********************************************************************
conformersearch.h - Conformer searching using genetic algorithm.

Copyright (C) 2010 Tim Vandermeersch
Some portions Copyright (C) 2012 Materials Design, Inc.
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

#ifndef OB_CONFORMERSEARCH_H
#define OB_CONFORMERSEARCH_H

#include <openbabel/mol.h>
#include <openbabel/rotor.h>
#include <openbabel/rotamer.h>

namespace OpenBabel {

  typedef std::vector<int> RotorKey;
  typedef std::vector<RotorKey> RotorKeys;
  typedef std::map<std::vector<int>,double> mapRotorEnergy;

  ///@addtogroup conformer Conformer Searching
  ///@{

  //////////////////////////////////////////////////////////
  //
  //  OBConformerFilter(s)
  //
  //////////////////////////////////////////////////////////

  /**
   * @class OBConformerFilter conformersearch.h <openbabel/conformersearch.h>
   * @brief Interface used by OBConformerSearch for filtering conformers.
   *
   * The OBConformerFilter class defines an interface used by the OBConformerSearch
   * class. A conformer filter is used as a first selection method when generating
   * conformers. A typical example is a steric filter to ensure atoms are never
   * closer than some specified distance.
   *
   * @see OBConformerSearch OBConformerFilters
   * @since 2.3
   */
  class OBAPI OBConformerFilter
  {
    public:
      /**
       * Check if a conformer is good (i.e. passes the filter). The molecule, rotor
       * key and coords are passed as parameters to this function and can be used
       * for checking the confoormer.
       * @return True if the conformer passes the filter.
       */
      virtual bool IsGood(const OBMol &mol, const RotorKey &key, double *coords) = 0;
      virtual ~OBConformerFilter() = 0;
  };

  /**
   * @class OBConformerFilters conformersearch.h <openbabel/conformersearch.h>
   * @brief Class for combining OBConformerFilter objects.
   *
   * The OBConformerFilters class makes it easy to combine OBConformerFilter
   * objects. A list of OBConformerFilter objects is specified when constructing
   * this class. The IsGood implementation simply checks all the specified
   * filters and returns false if any of the filters returns false.
   *
   * @see OBConformerSearch OBConformerFilter
   * @since 2.3
   */
  class OBAPI OBConformerFilters : public OBConformerFilter
  {
    public:
      /**
       * Constructor specifiying the filters to combine.
       */
      OBConformerFilters(const std::vector<OBConformerFilter*> &filters) : m_filters(filters)
      {
      }
      /**
       * IsGood reimplementation.
       * @return True if all filters pass.
       */
      bool IsGood(const OBMol &mol, const RotorKey &key, double *coords)
      {
        for (unsigned int i = 0; i < m_filters.size(); ++i)
          if (!m_filters[i]->IsGood(mol, key, coords))
            return false;
        return true;
      }
    protected:
      std::vector<OBConformerFilter*> m_filters;
  };

  /**
   * @class OBStericConformerFilter conformersearch.h <openbabel/conformersearch.h>
   * @brief A steric conformer filter class.
   *
   * The OBStericConformerFilter class is the default filter used by OBConformerSearch.
   * The filter removes all molecules which have at least 2 atoms closer together than
   * the specified distance.
   * @since 2.3
   */
  class OBAPI OBStericConformerFilter : public OBConformerFilter
  {
    public:
      OBStericConformerFilter ();
      OBStericConformerFilter (double cutoff, double vdw_factor = 0.5, bool check_hydrogens = true);
      bool IsGood(const OBMol &mol, const RotorKey &key, double *coords);
    private:
      double m_cutoff; //!< Internal cutoff (used as a squared distance)
      double m_vdw_factor;		//!< Factor applied to Van der Waals distance check
      bool m_check_hydrogens;		//!< Check hydrOgens if true 
  };

  //////////////////////////////////////////////////////////
  //
  //  OBConformerScore(s)
  //
  //////////////////////////////////////////////////////////

  /**
   * @class OBConformerScore conformersearch.h <openbabel/conformersearch.h>
   * @brief Interface used by OBConformerSearch for scoring conformers.
   *
   * The OBConformerScore class defines an interface to assign scores to conformers.
   * This class is used by OBConformerSearch to make the class flexible. A higher
   * score means the conformer with index is more favourable. A typical example is
   * the force field energy to find the lowest energy conformer (note: the energy is
   * not directly usable as score since lower values are better). Another example
   * is to use some measure of diversity (e.g. RMSD) to generate a diverse set of
   * conformers.
   * @since 2.3
   */
  class OBAPI OBConformerScore
  {
    public:
      /**
       * Conformer scores can be preferably high or low.
       */
      enum Preferred { HighScore, LowScore };
      /**
       * Preferred order for subclass scoring function.
       */
      virtual Preferred GetPreferred() = 0;
      /**
       * Convergence criteria used.
       */
      enum Convergence { Highest, Lowest, Sum, Average };
      /**
       * Convergence criteria for subclass scoring function.
       */
      virtual Convergence GetConvergence() = 0;
      /**
       * Score an individual conformer specified by index.
       */
      virtual double Score(OBMol &mol, unsigned int index, const RotorKeys &keys,
          const std::vector<double*> &conformers) = 0;
      virtual ~OBConformerScore() = 0;
  };

  /**
   * @class OBRMSDConformerScore conformersearch.h <openbabel/conformersearch.h>
   * @brief A RMSD conformer scoring class.
   *
   * Score conformers by the RMSD between the conformer with specified index and
   * the closest conformer. This results in a diverse set of conformers.
   * @since 2.3
   */
  class OBAPI OBRMSDConformerScore : public OBConformerScore
  {
    public:
      Preferred GetPreferred() { return HighScore; }
      Convergence GetConvergence() { return Average; }
      double Score(OBMol &mol, unsigned int index, const RotorKeys &keys,
          const std::vector<double*> &conformers);
  };

  /**
   * @class OBEnergyConformerScore conformersearch.h <openbabel/conformersearch.h>
   * @brief A lowest energy conformer scoring class.
   * @since 2.3
   */
  class OBAPI OBEnergyConformerScore : public OBConformerScore
  {
    public:
      OBEnergyConformerScore () {
	energy_map.clear ();
	energy_ncompute = 0;
	energy_nrequest = 0;
      }
      long unsigned int GetNbEnergyCompute () {return energy_ncompute;}
      long unsigned int GetNbEnergyRequest () {return energy_nrequest;}
      Preferred GetPreferred() { return LowScore; }
      Convergence GetConvergence() { return Lowest; }
      double Score(OBMol &mol, unsigned int index, const RotorKeys &keys,
          const std::vector<double*> &conformers);
    private:
      mapRotorEnergy energy_map;
      long unsigned int energy_ncompute;
      long unsigned int energy_nrequest;
  };
  
  /**
   * @class OBMinimizingEnergyConformerScore conformersearch.h <openbabel/conformersearch.h>
   * @brief A lowest energy conformer scoring class (after minimization)
   * @since 2.4
   */
  class OBAPI OBMinimizingEnergyConformerScore : public OBConformerScore
  {
    public:
      OBMinimizingEnergyConformerScore () {
	energy_map.clear ();
	energy_ncompute = 0;
	energy_nrequest = 0;
      }
      long unsigned int GetNbEnergyCompute () {return energy_ncompute;}
      long unsigned int GetNbEnergyRequest () {return energy_nrequest;}
      Preferred GetPreferred() { return LowScore; }
      Convergence GetConvergence() { return Lowest; }
      double Score(OBMol &mol, unsigned int index, const RotorKeys &keys,
          const std::vector<double*> &conformers);
    private:
      mapRotorEnergy energy_map;
      long unsigned int energy_ncompute;
      long unsigned int energy_nrequest;
  };

  /**
   * @class OBMinimizingRMSDConformerScore conformersearch.h <openbabel/conformersearch.h>
   * @brief An RMSD conformer scoring class, after a short minimization
   *
   * This scores conformers by the RMSD between the conformer and the closest, to produce a
   * diverse set of conformers, but after minimization. This ensures each conformer is
   * "reasonable" and avoids steric clashes.
   * @since 2.4
   */
  class OBAPI OBMinimizingRMSDConformerScore : public OBConformerScore
  {
    public:
      Preferred GetPreferred() { return HighScore; }
      Convergence GetConvergence() { return Average; }
      double Score(OBMol &mol, unsigned int index, const RotorKeys &keys,
          const std::vector<double*> &conformers);
  };

  //////////////////////////////////////////////////////////
  //
  //  OBConformerSearch
  //
  //////////////////////////////////////////////////////////

  /**
   * @class OBConformerSearch conformersearch.h <openbabel/conformersearch.h>
   * @brief Conformer searching using genetic algorithm.
   * See @ref ConformerSearching
   * @since 2.3
   */
  class OBAPI OBConformerSearch
  {
    public:
      OBConformerSearch();
      ~OBConformerSearch();
      /**
       * Setup this instance with the specified molecule.
       *
       * @param mol The molecule with initial conformer.
       * @param numConformers The number of conformers that should be generated. This
       *        is also the number of conformers selected for each generation.
       * @param numChildren When a new generation is generated, for each of the
       *        numConformers conformers, numChildren children are created.
       * @param mutability The mutability determines how frequent a permutation occurs
       *        when generating the next generation.
       * @param convergence The number of identical generations before considering
       *        the process converged.
       */
      bool Setup(const OBMol &mol, int numConformers = 30, int numChildren = 5,
          int mutability = 5, int convergence = 25);
      /**
       * Set the number of conformers.
       */
      void SetNumConformers(int numConformers) { m_numConformers = numConformers; }
      /**
       * Set the number of children generated for each parent.
       */
      void SetNumChildren(int numChildren) { m_numChildren = numChildren; }
      /**
       * Set the mutability.
       */
      void SetMutability(int mutability) { m_mutability = mutability; }
      /**
       * Set the convergence (i.e. number of generations that did not change the score
       * before considering the iteration converged).
       */
      void SetConvergence(int convergence) { m_convergence = convergence; }
      /**
       * Set the bonds to be fixed.
       */
      void SetFixedBonds(const OBBitVec &fixedBonds) { m_fixedBonds = fixedBonds; }


      /**
       * Set the filter method used to check if a newly generated is acceptable. Typical
       * examples are a steric filter or electrostatic energy filter. The filters make a binary
       * selection and one group will be scored in the next step. A number of filters
       * are provided by default but implementing a custom filter is easy.
       *
       * There is also a OBConformerFilters class to make it easy to use multiple filters.
       */
      void SetFilter(OBConformerFilter *filter)
      {
        delete m_filter;
        m_filter = filter;
      }
      /**
       * All acceptable conformers are scored to select the fittest conformers from
       * the population. A higher score means the conformer is desired in the final
       * set or next generation. Typical scoring would be force field energy to find
       * the lowest energy conformer. Another example is the highest RMSD between
       * conformer i and the closest conformer to create a diverse set of conformers.
       * Diversity could also be expressed in other ways and implementing a scoring
       * class is easy.
       */
      void SetScore(OBConformerScore *score)
      {
        delete m_score;
        m_score = score;
      }

      /**
      * Set whether or not you want rotors to be printed prior to the conformer search.
      */
      void PrintRotors(bool printrotors) { m_printrotors = printrotors; }

      /**
       * Perform conformer search using a genetic algorithm.
       */
      void Search();

      const RotorKeys& GetRotorKeys() const
      {
        return m_rotorKeys;
      }

      void GetConformers(OBMol &mol);

      /* @brief Set an output stream for logging. If NULL pointer is provided, logging is disabled. */
      void SetLogStream  (std::ostream *sptr) {m_logstream = sptr;}

      /*************************************************/
      /* Methods related to fitness sharing parameters */
      
      /* @brief Set the use of fitness sharing ON (default) or OFF*/
      void SetSharing (bool value = true)  {use_sharing = value;}
      
      /* @brief Get the targeted number of niches, for dynamic niche sharing */
      int GetNbNiches () {return nb_niches;}
      
      /* @brief Set the targeted number of niches, for dynamic niche sharing */
      void SetNbNiches (int value) {nb_niches = value;}
      
      /* @brief Get niches radius, for dynamic niche sharing.*/
      double GetNicheRadius () {return niche_radius;}
      
      /* @brief Set niches radius, for dynamic niche sharing.*/
      void SetNicheRadius (double value) {niche_radius = value;}
      
      /* @brief Get the alpha sharing parameter */
      double GetAlphaSharing () {return alpha_share;}
      
      /* @brief Set the alpha sharing parameter */
      void SetAlphaSharing (double value) {alpha_share = value;}
      
      /* @brief Get the sigma sharing parameter */
      double GetSigmaSharing () {return sigma_share;}
      
      /* @brief Set the sigma sharing parameter */
      void SetSigmaSharing (double value) {sigma_share = value;}
      
      /* @brief Get the (uniform) crossover probability */
      double GetCrossoverProbability () {return p_crossover;}
      
      /* @brief Set the (uniform) crossover probability */
      void SetCrossoverProbability (double value) {p_crossover = value;}
      
      /* @brief Get the niche mating probability, for dynamic niche sharing */
      double GetNicheMating () {return niche_mating;}
      
      /* @brief Set the (uniform) crossover probability */
      void SetNicheMating (double value) {niche_mating = value;}
      
      /* @brief Set the local optimization rate */
      void SetLocalOptRate (int value) {local_opt_rate = value;}
      
      /* @brief Get the local optimization rate*/
      int  SetLocalOptRate() {return local_opt_rate;}
    
    private:
      /**
       * Create the next generation.
       */
      void NextGeneration();
      /**
       * Select the fittest numConformers from the parents and their children.
       */
      double MakeSelection();
      /**
       * Check if a conformer key is unique.
       */
      bool IsUniqueKey(const RotorKeys &keys, const RotorKey &key) const;
      /**
       * Check the specified rotor key using the set filter.
       */
      bool IsGood(const RotorKey &key);

      //! @brief Genetic similarity measure, i.e. "distance" between two rotor keys.
      int key_distance (const RotorKey &key1, const RotorKey &key2);      
      //! @brief Make a local search on the best individual
      int local_opt ();
      //! @brief Produces one or two offsprings
      int reproduce (RotorKey &key1, RotorKey &key2);
      //! @brief Score, sort the current population, compute shared fitness values (and niches)
      int score_population ();
      //! @brief Compute shared fitness values for a given population
      int share_fitness ();
      //! @brief Perform one generation with fitness sharing
      double sharing_generation ();

      unsigned int m_numConformers; //!< The desired number of conformers. This is also the population size.
      int m_numChildren; //!< The number of children generated each generation
      int m_mutability; //!< The mutability for generating the next generation
      int m_convergence; //!< Number of generations that remain unchanged before quiting
      
      std::vector<double> vscores;                    //!< Current population score vector
      std::vector<double> vshared_fitnes;             //!< Current population shared fitness vector
      std::vector<std::vector <int> > dynamic_niches; //!< The dynamically found niches
      std::vector<int> niche_map;                     //!< Procide the sharing niche index, given the key inddex
      
      void *d; // Opaque pointer - currently for storing OBRandom* which may be removed in future
      bool use_sharing;		//!< Whether to use sharing or not.
      double alpha_share;	//!< The alpha parameter in sharing function
      double sigma_share;		//!< The sigma parameter in sharing function
      int nb_niches;		//!< The number of dynamic niches to be found
      double niche_radius;		//!< A pre-determined niche radius, for dynamic niche sharing.
      double p_crossover;	//!< Crossover probability
      double niche_mating;	//!< Probability of forcing the second parent in the first parent
      int local_opt_rate;       //!< Perform a random local optimization every local_opt_rate generations. Disabled if set to 
      OBBitVec      m_fixedBonds; //!< Bonds that are fixed
      OBMol         m_mol; //!< The molecule with starting coordinates
      OBRotorList   m_rotorList; //!< The OBRotorList for the molecule
      RotorKeys     m_rotorKeys; //!< The current population
      bool          m_printrotors; //!< Wheter or not to print all rotors that are found instead of performing the conformer search

      OBConformerFilter *m_filter;
      OBConformerScore  *m_score;

      std::ostream *m_logstream;	//!< A pointer to a log stream (NULL means no loogging)
  };

  /**
   * @page ConformerSearching
   *
   * @section introduction Introduction
   * All conformer searching methods in OpenBabel are based on the concept of rotor
   * keys. A rotor key is simply an array of values specifying the rotations around
   * rotatable bonds. The number of possible values for a rotatable bond depends on
   * the bond type. The classes implementing the rotor key concept are OBRotor,
   * OBRotorList and OBRotamer.
   *
   * Previous OpenBabel releases contained only methods for finding stable (low
   * energy) conformers by using the force fields. The 2.3 release introduces a new
   * flexible class (OBConformerSearch) implementing a genetic algorithm. The scoring
   * or ranking of conformers is done by a separate class derived from the abstract
   * OBConformerScore class. Reimplementing this class allows for all sorts of scoring
   * functions (e.g. RMSD, torson, energy, ... based).
   *
   * @section stable Finding stable conformers
   *
   * @subsection systematic Systematic rotor search
   * This is the simplest of all conformer searching methods and only suitable when
   * there are only a few rotatable bonds. Since the number of conformers grows
   * exponentially with the number of bonds, enumerating all conformers for reasonably
   * sized molecules quickly results in unacceptable performance. A systematic rotor
   * search can be done with the OBForceField::SystematicRotorSearch method.
   *
   * @subsection random Random rotor search
   * The random rotor search generates random rotor keys and evaluates the energy.
   * The lowest energy conformer is kept. All rotor keys are random and no attempt
   * is made to improve good scoring keys. In most cases the weighted rotor search
   * is a better choice. A random rotor search can be performed with the
   * OBForceField::RandomRotorSearch method.
   *
   * @subsection weighted Weighted rotor search
   * The weighted rotor search uses the score (energy) from a generated and tries
   * to optimize good scoring keys. This allows the global minimum energy conformer
   * to be identified much faster compared to the random rotor search. However,
   * finding the global minimum can be time consuming for large molecules and there
   * is no guarantee that the found minimum is the global minimum. A weighted rotor
   * search is performed with the OBForceField::WeightedRotorSearch method.
   *
   * @subsection genetic Genetic algorithm
   * The OBConformerSearch class introduced in OpenBabel 2.3 implements a genetic
   * algorithm for finding conformers. The class is configurable and with the
   * right scoring function (OBEnergyConformerScore), stable conformers can be
   * found. Read the OBConformerSearch section for more details on how the genetic
   * algorithm works.
   *
   * @section diverse Finding a diverse set of conformers
   * A divers set of conformers can be generated using the OBConformerSearch class
   * with a scoring function like OBRMSDConformerScore. A diverse set of conformers
   * is often useful when screening bioactive molecules since the interaction with
   * the target can stabilize a higher conformational energy. See the next section
   * for details.
   *
   * @section OBConformerSearch
   * The genetic algorithm starts by generating the initial population of rotor keys.
   * The initial population contains up to numConformers, duplicate keys are ignored.
   * OBConformerFilter objects are used to filter conformers before including them.
   * A typical filter is a steric filter to ignore all conformers with atoms to close
   * together.
   *
   * For each generation, numChildren children are created by permuting the parent
   * rotor keys. The mutability setting determines how frequent a permutation is made
   * (e.g. 5 means 1/5 bonds are permuted, 10 means 1/10).
   * Again, duplicated and filtered molecules are ignored. The population now contains
   * up to numConformer * (1 + numChildren). All these rotor keys are scored using the
   * specified OBConformerScore class. Next, all keys are ordered by their score and
   * the best numConformers conformers are selected as parents for the next
   * generation.
   *
   * New generations are generated until the specified number of generations (i.e.
   * convergence) don't improve the score.
   *
   * @subsection OBConformerFilter
   * Filters are used to exclude certain generated conformers. For example, the
   * OBStericConformerFilter filters out all conformers with atoms to close together.
   * Custom filters can be used to consider additional properties (e.g. electrostatic
   * energy, distance between pharmacophore groups, ...). Filters can easily be
   * combined using the OBConformerFilters class.
   *
   * @subsection OBConformerScore
   * OBConformerScore derived classes compute the score for conformers. The scoring
   * function is probably the most important part of the algorithm since it selects
   * which conformers will be used as parents for the next generation. Derived
   * classes also specify if the scores are descending or ascending and what
   * convergence criteria should be used (i.e. lowest, highest, sum or average).
   * There is no easy way to combine scoring functions since scores are not always
   * additive. However, a custom scoring function could be used that combines several
   * other scoring functions with weights for each score.
   *
   * The scoring function determines what kind of conformers are selected. For
   * example, the OBEnergyConformerScore class selects stable, low energy conformers.
   * A simple variation to this would be to only use the non-bonded energy. The
   * result for these two scoring functions should be similar.
   *
   * The OBRMSDConformerScore scoring function can be used to generate a diverse set
   * of conformers. A diverse set is often useful when screening bioactive molecules
   * since the interaction with the target can stabilize a higher conformational
   * energy. To calculate the RMSD between conformers, Kabsch alignment (OBAlign)
   * is used. The alignment takes symmetry into account.
   *
   * @image html energyconformerscore.png
   * @image html rmsdconformerscore.png "30 Conformers for methotrexate. top: OBEnergyConformerScore bottom: OBRMSDConformerScore"
   */

  ///@}
};

#endif

/// @file conformersearch.h
/// @brief Conformer searching using genetic algorithm

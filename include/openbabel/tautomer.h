/**********************************************************************
tautomer.h - Tautomer support

  Copyright (C) 2011,2020 by Tim Vandermeersch

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

#include <openbabel/babelconfig.h>
#ifndef OBAPI
  #define OBAPI
#endif

#include <string>
#include <vector>

namespace OpenBabel {

  class OBMol;

  /**
   * Functor for tautomer enumeration algorithm. Subclass this class to write
   * your own functor. This functor will be called for every tautomer that is
   * found. This may include duplicates. Subclass UniqueTautomerFunctor to
   * filter these duplicates.
   *
   * Example to write all tautomer smiles to std::cout:
   * @code
   * class WriteSmilesFunctor : public TautomerFunctor
   * {
   *   public:
   *     WriteSmilesFunctor()
   *     {
   *       m_conv.SetOutFormat("can");
   *     }
   *
   *     void operator()(OBMol *mol)
   *     {
   *       std::cout << m_conv.WriteString(mol, true) << std::endl;
   *     }
   *
   *   private:
   *     OBConversion m_conv;
   * };
   * @endcode
   *
   * This functor can be used in the following way:
   * @code
   * OBMol mol;
   * ...
   * WriteSmilesFunctor functor;
   * // write all tautomers to std::cout
   * EnumerateTautomers(mol, functor);
   * @endcode
   */
  class OBAPI TautomerFunctor
  {
    public:
      virtual ~TautomerFunctor() {}
      /**
       * This function is called every time a tautomer is discovered.
       *
       * @param mol The tautomer.
       */
      virtual void operator()(OBMol *mol) = 0;
  };

  /**
   * Functor for tautomer enumeration algorithm. Subclass this class to write
   * your own functor. This functor will be called for every unique tautomer
   * that is found.
   *
   * Example to write all unique tautomer smiles to std::cout:
   * @code
   * class WriteSmilesFunctor : public UniqueTautomerFunctor
   * {
   *   public:
   *     void operator()(OBMol *mol, const std::string &smiles)
   *     {
   *       std::cout << smiles << std::endl;
   *     }
   * };
   * @endcode
   *
   * This functor can be used in the following way:
   * @code
   * OBMol mol;
   * ...
   * WriteSmilesFunctor functor;
   * // write all unique tautomers to std::cout
   * EnumerateTautomers(mol, functor);
   * @endcode
   */
  class OBAPI UniqueTautomerFunctor : public TautomerFunctor
  {
    public:
      virtual ~UniqueTautomerFunctor() {}
      /**
       * This function is called every time a tautomer is discovered.
       *
       * @param mol The tautomer.
       */
      void operator()(OBMol *mol);
      /**
       * This function is called every time a unique tautomer is discovered.
       *
       * @param mol The tautomer.
       * @param smiles The canonical SMILES for the tautomer.
       */
      virtual void operator()(OBMol *mol, const std::string &smiles) = 0;
    private:
      std::vector<std::string> m_smiles; //!< unique tautomer SMILES
  };

  /**
   * Enumerate all tautomers for @p mol. Every time a tautomer is discovered,
   * the @p functor will be invoked with the bonds changed to the tautomer.
   * When the enumeration is complete, the bonds will be changed back to the
   * original bond orders.
   *
   * The algorithm is based on
   * http://www.daylight.com/meetings/emug99/Delany/taut_html/index.htm
   *
   * @warning This function makes hydrogens implicit.
   * @see TautomerFunctor
   */
  void OBAPI EnumerateTautomers(OBMol *mol, TautomerFunctor &functor);

  /**
   * Compuate the canonical tautomer for @p mol. The canonical tautomer is the
   * first found tautomer in the enumeration algorithm and is a mathematically
   * unique tautomer (i.e. the algorithm does not guarantee that this is the
   * dominant tautomer). The bonds in @p mol will be changed to the canonical
   * tautomer.
   *
   * @warning This function makes hydrogens implicit.
   */
  void OBAPI CanonicalTautomer(OBMol *mol);

}

//! \file tautomer.h
//! \brief Tautomer support.

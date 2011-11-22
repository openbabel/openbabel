/**********************************************************************
tautomer.h - Tautomer support

  Copyright (C) 2011 by Tim Vandermeersch

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

namespace OpenBabel {

  class OBMol;

  /**
   * Functor for tautomer enumeration algorithm. Subclass this class to write
   * your own functor.
   *
   * Example to write all tautomer smiles to std::cout:
   * @code
   * class WriteSmilesFunctor : public TautomerFunctor
   * {
   *   public:
   *     WriteSmilesFunctor() : foundTautomer(false)
   *     {
   *       m_conv.SetOutFormat("can");
   *     }
   *
   *     void operator()(OBMol *mol)
   *     {
   *       foundTautomer = true;
   *       std::cout << m_conv.WriteString(mol);
   *     }
   *
   *     bool foundTautomer;
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
   * // write the original smiles if there are no tautomers found
   * if (!functor.foundTautomers) {
   *   OBConversion conv;
   *   conv.SetOutFormat("can");
   *   std::cout << conv.WriteString(&mol);
   * }
   * @endcode
   */  
  class OBAPI TautomerFunctor
  {
    public:
      virtual ~TautomerFunctor() {}
      /**
       * This function is called every time a tautomer is discovered.
       */
      virtual void operator()(OBMol *mol) = 0;
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

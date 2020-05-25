/**********************************************************************
rand.h - Pseudo random number generator.

Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2006 by Geoffrey R. Hutchison

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

#ifndef RAND_H
#define RAND_H

#include <openbabel/babelconfig.h>
#include <memory>
#include <random>

namespace OpenBabel
{

  //! \class OBRandom rand.h
  //! \brief Random number generator
  /**
     The OBRandom class can be used to facilitate cross-platform random number
     generation. The class can be set to specific seed states, or set to
     use the current time or other arbitrary data as a seed.

     \code
     OBRandom generator;        // use a random seed
     //OBRandom generator{42};  // specify a seed
     //generator.Seed(42);      // reset by a seed

     cerr << "New Random Integer " << generator.UniformInt(0, 6) << endl;
     cerr << "New Random Floating-Point " << generator.UniformReal(0.0, 1.0) << endl;
     \endcode
   **/
  class OBRandom
  {
  public:
    OBRandom() { TimeSeed(); }
    explicit OBRandom(uint_fast64_t seed) { Seed(seed); }
    //! \deprecated
    explicit OBRandom(bool OBRandomUseSysRand) { TimeSeed(); }
    //! Use @p seed for the random number generator seed
    void Seed(uint_fast64_t seed) { prng.reset(new std::mt19937_64{seed}); }
    //! Use the current time for the random number generator seed
    //! If sranddev is available (e.g., Mac OS X, BSD...) use this instead
    //! for more random seeds
    void TimeSeed();
    //! \return a random integer in [0, 2^31 - 1] \deprecated
    int NextInt();
    //! \return a random integer in [a, b]
    int UniformInt(int a, int b);
    //! \return a random floating-point number in [0, 1) \deprecated
    double NextFloat();
    //! \return a random floating-point number in [a, b)
    double UniformReal(double a, double b);

  private:
    std::unique_ptr<std::mt19937_64> prng;
  };

} // end namespace OpenBabel

#endif // RAND_H

//! \file rand.h
//! \brief Pseudo random number generator

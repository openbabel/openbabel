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

namespace OpenBabel
{

  //******************************************
  //*** Stuff for random number generation ***
  //******************************************

  //! \struct DoubleType rand.h <openbabel/rand.h>
  //! \brief Used for internal random number generation OBRandom (unless the system random generator is used)
  typedef struct
  {
    unsigned int hi;
    unsigned int lo;
  }
  DoubleType;

  
  //! \class OBRandom rand.h <openbabel/rand.h>
  //! \brief Random number generator
  /**
     The OBRandom class can be used to facilitate cross-platform random number
     generation. The class can be set to specific seed states, or set to
     use the current time or other arbitrary data as a seed.

     \code
     OBRandom generator; // Don't use system rand() functions
     generator.TimeSeed(); // initialize from the current time

     cerr << " New Random Integer " << generator.NextInt() << endl;
     cerr << " New Random Floating-Point " << generator.NextFloat() << endl;
     \endcode

     Alternatively, OBRandom can be used as an interface to the system random
     number generator.

     \code
     OBRandom generator(true); // Use system rand() functions
     generator.Seed(10246);// Use a specific initial seed value for reproducing sequence
     \endcode
   **/
  class OBRandom
  {
    DoubleType d;
    unsigned int m,a,c;
    unsigned int p;
    unsigned int i;
    unsigned int x;
    bool OBRandomUseSysRand;

  public:
    //! \brief Constructor. @p useSys will use the system rand() function
    OBRandom(bool useSys= false);
    //! Use @p seed for the random number generator seed
    void Seed(int seed)
    {
      x = seed;
    }
    //! Use the current time for the random number generator seed
    //! If sranddev is available (e.g., Mac OS X, BSD...) use this instead
    //! for more random seeds
    void TimeSeed();
    //! \return a random integer
    int NextInt();
    //! \return a random floating-point number between 0.0 and 1.0
    double NextFloat();
  };

} // end namespace OpenBabel

#endif // RAND_H

//! \file rand.h
//! \brief Pseudo random number generator

/**********************************************************************
rand.cpp - Pseudo random number generator.

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

#include "rand.h"

#include <cstdlib>

#if TIME_WITH_SYS_TIME
#include <sys/time.h>
#include <ctime>
#else
#if HAVE_SYS_TIME_H
#include <sys/time.h>
#else
#include <ctime>
#endif
#endif


namespace OpenBabel
{

  void OBRandom::Reset()
  {
    auto ob_random_seed = std::getenv("OB_RANDOM_SEED");
    if (ob_random_seed) {
      this->Seed(std::atoll(ob_random_seed));
      return;
    }
#if defined(WIN32) || defined(__MINGW32__)
    // for VC++ do it this way
    time_t ltime;
    time(&ltime);
    const long secs= long(ltime);
    this->Seed(secs);
#else
    // XXX random_device can be implemented as a PRNG in the above envirnment
    this->Seed(std::random_device{}());
#endif
  }

  int OBRandom::UniformInt(int a, int b)
  {
    std::uniform_int_distribution<int> u{a, b};
    return u(*prng);
  }

  double OBRandom::UniformReal(double a, double b)
  {
    std::uniform_real_distribution<double> u{a, b};
    return u(*prng);
  }

  bool OBRandom::Bernoulli(double p)
  {
    std::bernoulli_distribution b{p};
    return b(*prng);
  }

} //end namespace OpenBabel

//! \file rand.cpp
//! \brief Pseudo random number generator.

/**********************************************************************
rand.h - Pseudo random number generator.
 
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2006 by Geoffrey R. Hutchison
 
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

#ifndef RAND_H
#define RAND_H

#include <openbabel/babelconfig.h>

namespace OpenBabel
{

  //******************************************
  //*** Stuff for random number generation ***
  //******************************************

  //! Used for internal random number generation OBRandom (unless the system random generaor is used)
  typedef struct
  {
    unsigned int hi;
    unsigned int lo;
  }
  DoubleType;

  OBAPI void DoubleMultiply( unsigned int,unsigned int,DoubleType*);
  OBAPI void DoubleAdd( DoubleType*,unsigned int);
  OBAPI unsigned int DoubleModulus( DoubleType*,unsigned int);

  //! Random number generator
  class OBAPI OBRandom
  {
    DoubleType d;
    unsigned int m,a,c;
    unsigned int p;
    unsigned int i;
    unsigned int x;
    bool OBRandomUseSysRand;

  public:
    OBRandom(bool useSys= false);
    void Seed(int seed)
    {
      x = seed;
    }
    void TimeSeed();
    int NextInt();
    double NextFloat();
  };

} // end namespace OpenBabel

#endif // RAND_H

//! \file rand.h
//! \brief Pseudo random number generator

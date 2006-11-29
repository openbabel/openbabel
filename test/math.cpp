/**********************************************************************
math.cpp - Unit tests for the math code in OpenBabel

Copyright (C) 2005-2006 Geoffrey R. Hutchison
Copyright (C) 2006 Benoit Jacob
 
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

// used to set import/export for Cygwin DLLs
#ifdef WIN32
#define USING_OBDLL
#endif

#include "babelconfig.h"
#include "math/matrix3x3.h"
#include "obutil.h"

#include <iostream>
#include <stdlib.h>

#define REPEAT 100

# if defined __cplusplus ? __GNUC_PREREQ (2, 6) : __GNUC_PREREQ (2, 4)
#   define FUNCTION    __PRETTY_FUNCTION__
# else
#  if defined __STDC_VERSION__ && __STDC_VERSION__ >= 199901L
#   define FUNCTION    __func__
#  else
#   define FUNCTION    ""
#  endif
# endif

#define VERIFY(expr) \
  if( ! (expr) ) verify_failed( __STRING(expr), FUNCTION, __FILE__, __LINE__ )

#define TEST(func) \
  for( int i = 0; i < REPEAT; i++ ) testBasics_vector3(); \
  cout << "passed " << __STRING(func) << endl

using namespace std;
using namespace OpenBabel;

OBRandom randomizer;

void verify_failed( const char *expr, const char *func, const char *file, int line )
{
  cout << "Failed test ( " << expr
       << " ) in line " << line
       << " in function " << func
       << " in file " << file << endl;
  exit(1);
}

void pickRandom( double & d )
{
  d = randomizer.NextFloat() * 2.0 - 1.0;
}

void pickRandom( vector3 & v )
{
  pickRandom( v.x() );
  pickRandom( v.y() );
  pickRandom( v.z() );
}

void pickRandom( matrix3x3 & m )
{
  for( int row = 0; row < 3; row++ )
    for( int column = 0; column < 3; column++ )
      pickRandom( m( row, column ) );
}

bool compare( const double & d1, const double & d2 )
{
  return IsApprox( d1, d2, 1e-6 );
}

bool compare( const vector3 & v1, const vector3 & v2 )
{
  return v1.IsApprox( v2, 1e-6 );
}

bool compare( const matrix3x3 & m1, const matrix3x3 & m2 )
{
  bool ret = true;
  for( int i = 0; i < 3; i++ )
    ret &= compare( m1.GetColumn(i), m2.GetColumn(i) );
  return ret;
}

// tests the constructors, accessors, mutators of class vector3.
void testBasics_vector3()
{
  double d[3];
  vector3 v1, v2, v3;

  // test compare() itself, assuming Set() works
  v1.Set( 1.0, 0.0, 0.0 );
  v2.Set( 1.0, 1e-20, 0.0 );
  v3.Set( 1.0, 0.1, 0.0 );
  VERIFY( compare( v1, v2 ) );
  VERIFY( ! compare( v2, v3 ) );

  // test constructors and operator=
  pickRandom( d[0] );
  pickRandom( d[1] );
  pickRandom( d[2] );
  v1 = vector3( d[0], d[1], d[2] );
  v2 = vector3( v1 );
  v3 = v1;
  VERIFY( compare( v1, v2 ) );
  VERIFY( compare( v1, v3 ) );

  // test constant operator()
  VERIFY( compare( static_cast<const vector3>(v1).x(), d[0] ) );
  VERIFY( compare( static_cast<const vector3>(v1).y(), d[1] ) );
  VERIFY( compare( static_cast<const vector3>(v1).z(), d[2] ) );

  // test Get
  v1.Get( d );
  VERIFY( compare( v1.x(), d[0] ) );
  VERIFY( compare( v1.y(), d[1] ) );
  VERIFY( compare( v1.z(), d[2] ) );

  // test non-constant operator(), Set, SetX, SetY, SetZ
  pickRandom( d[0] );
  pickRandom( d[1] );
  pickRandom( d[2] );
  v1.x() = d[0];
  v1.y() = d[1];
  v1.z() = d[2];
  VERIFY( compare( v1.x(), d[0] ) );
  VERIFY( compare( v1.y(), d[1] ) );
  VERIFY( compare( v1.z(), d[2] ) );
  v2.Set( d );
  VERIFY( compare( v1, v2 ) );
  v3.SetX( d[0] );
  v3.SetY( d[1] );
  v3.SetZ( d[2] );
  VERIFY( compare( v1, v3 ) );
}

void testBasics_matrix3x3()
{
  double d[3][3] = { { 1.0, 0.0, 0.0 },
                     { 0.0, 1.0, 0.0 },
                     { 0.0, 0.0, 1.0 } };
  double x;
  matrix3x3 m1, m2, m3;
  vector3 v1, v2;
  int i, j;

  // test constructors and operator=
  m1 = matrix3x3(1.0);
  m2 = matrix3x3(d);
  VERIFY( compare( m1, m2 ) );

  // test GetArray
  static double e[3][3]; // without "static", compiler optimizations can
                         // produce bugs.
  pickRandom( m2 );
  m2.GetArray( & e[0][0] );
  m3 = matrix3x3(e);
  VERIFY( compare( m2, m3 ) );

  // test const operator() and Get()
  for( i = 0; i < 3; i++ )
    for( j = 0; j < 3; j++ )
    {
      VERIFY( compare( static_cast<const matrix3x3>(m2)(i,j), e[i][j] ) );
      x = m2.Get( i, j );
      VERIFY( compare( x, e[i][j] ) );
    }

  // test non-const operator() and Set()
  for( i = 0; i < 3; i++ )
    for( j = 0; j < 3; j++ )
      m1(i,j) = m2(i,j);
  VERIFY( compare( m1, m2 ) );
  for( i = 0; i < 3; i++ )
    for( j = 0; j < 3; j++ )
      m3.Set(i,j, m2(i,j) );
  VERIFY( compare( m3, m2 ) );

  // test SetColumn(), GetColumn(), SetRow(), GetRow()
  for( i = 0; i < 3; i++ )
  {
    pickRandom( v1 );
    m1.SetRow( i, v1 );
    v2 = m1.GetRow( i );
    VERIFY( compare( v1, v2 ) );
    m1.SetColumn( i, v1 );
    v2 = m1.GetColumn( i );
    VERIFY( compare( v1, v2 ) );
  }
}

int main(int argc,char *argv[])
{
  if (argc != 1)
  {
    cout << "Usage: math" << endl;
    cout << "   Tests OpenBabel's math code." << endl;
    return 0;
  }
  
  cout << "math: repeating each test " << REPEAT << " times" << endl;
  
  randomizer.TimeSeed();
  
  TEST(testBasics_vector3);
  TEST(testBasics_matrix3x3);

  cout << "math: all tests are successful" << endl;
  return 0;
}

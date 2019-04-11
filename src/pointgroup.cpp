/**********************************************************************
pointgroup.cpp - Brute force symmetry analyzer.

 (C) 1996, 2003 S. Patchkovskii, Serguei.Patchkovskii@sympatico.ca
 Some portions Copyright (C) 2007 by Geoffrey R. Hutchison
     (Ported to C++, integrated with Open Babel)

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

#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/pointgroup.h>
#include <openbabel/obiter.h>
#include <iostream>

#include <string>
#include <math.h>
#include <cstring>
#include <algorithm>

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795028841971694
#endif

#define	DIMENSION 3
#define MAXPARAM  7

namespace OpenBabel {

  /*
   *  All specific structures should have corresponding elements in the
   *  same position generic structure does.
   *
   *  Planes are characterized by the surface normal direction
   *  (taken in the direction *from* the coordinate origin)
   *  and distance from the coordinate origin to the plane
   *  in the direction of the surface normal.
   *
   *  Inversion is characterized by location of the inversion center.
   *
   *  Rotation is characterized by a vector (distance+direction) from the origin
   *  to the rotation axis, axis direction and rotation order. Rotations
   *  are in the clockwise direction looking opposite to the direction
   *  of the axis. Note that this definition of the rotation axis
   *  is *not* unique, since an arbitrary multiple of the axis direction
   *  can be added to the position vector without changing actual operation.
   *
   *  Mirror rotation is defined by the same parameters as normal rotation,
   *  but the origin is now unambiguous since it defines the position of the
   *  plane associated with the axis.
   *
   */

  typedef struct {
    const char *  group_name ;        /* Canonical group name                        */
    const char *  symmetry_code ;     /* Group symmetry code                         */
    int     (*check)( void ) ;        /* Additional verification routine, not used  */
  } POINT_GROUP ;

  /*
   *    Point groups I know about
   */
  POINT_GROUP            PointGroups[]         = {
    {  "C1",    ""},
    {  "Cs",    "(sigma) "},
    {  "Ci",    "(i) "},
    {  "C2",    "(C2) "},
    {  "C3",    "(C3) "},
    {  "C4",    "(C4) (C2) "},
    {  "C5",    "(C5) "},
    {  "C6",    "(C6) (C3) (C2) "},
    {  "C7",    "(C7) "},
    {  "C8",    "(C8) (C4) (C2) "},
    {  "D2",    "3*(C2) "},
    {  "D3",    "(C3) 3*(C2) "},
    {  "D4",    "(C4) 5*(C2) "},
    {  "D5",    "(C5) 5*(C2) "},
    {  "D6",    "(C6) (C3) 7*(C2) "},
    {  "D7",    "(C7) 7*(C2) "},
    {  "D8",    "(C8) (C4) 9*(C2) "},
    {  "C2v",   "(C2) 2*(sigma) "},
    {  "C3v",   "(C3) 3*(sigma) "},
    {  "C4v",   "(C4) (C2) 4*(sigma) "},
    {  "C5v",   "(C5) 5*(sigma) "},
    {  "C6v",   "(C6) (C3) (C2) 6*(sigma) "},
    {  "C7v",   "(C7) 7*(sigma) "},
    {  "C8v",   "(C8) (C4) (C2) 8*(sigma) "},
    {  "C2h",   "(i) (C2) (sigma) "},
    {  "C3h",   "(C3) (S3) (sigma) "},
    {  "C4h",   "(i) (C4) (C2) (S4) (sigma) "},
    {  "C5h",   "(C5) (S5) (sigma) "},
    {  "C6h",   "(i) (C6) (C3) (C2) (S6) (S3) (sigma) "},
    {  "C7h",   "(C7) (S7) (sigma) "},
    {  "C8h",   "(i) (C8) (C4) (C2) (S8) (S4) (sigma) "},
    {  "D2d",   "3*(C2) (S4) 2*(sigma) "},
    {  "D3d",   "(i) (C3) 3*(C2) (S6) 3*(sigma) "},
    {  "D4d",   "(C4) 5*(C2) (S8) 4*(sigma) "},
    {  "D5d",   "(i) (C5) 5*(C2) (S10) 5*(sigma) "},
    {  "D6d",   "(C6) (C3) 7*(C2) (S12) (S4) 6*(sigma) "},
    {  "D7d",   "(i) (C7) 7*(C2) (S14) 7*(sigma) "},
    {  "D8d",   "(C8) (C4) 9*(C2) (S16) 8*(sigma) "},
    {  "D2h",   "(i) 3*(C2) 3*(sigma) "},
    {  "D3h",   "(C3) 3*(C2) (S3) 4*(sigma) "},
    {  "D4h",   "(i) (C4) 5*(C2) (S4) 5*(sigma) "},
    {  "D5h",   "(C5) 5*(C2) (S5) 6*(sigma) "},
    {  "D6h",   "(i) (C6) (C3) 7*(C2) (S6) (S3) 7*(sigma) "},
    {  "D7h",   "(C7) 7*(C2) (S7) 8*(sigma) "},
    {  "D8h",   "(i) (C8) (C4) 9*(C2) (S8) (S4) 9*(sigma) "},
    {  "S4",    "(C2) (S4) "},
    {  "S6",    "(i) (C3) (S6) "},
    {  "S8",    "(C4) (C2) (S8) "},
    {  "T",     "4*(C3) 3*(C2) "},
    {  "Th",    "(i) 4*(C3) 3*(C2) 4*(S6) 3*(sigma) "},
    {  "Td",    "4*(C3) 3*(C2) 3*(S4) 6*(sigma) "},
    {  "O",     "3*(C4) 4*(C3) 9*(C2) "},
    {  "Oh",    "(i) 3*(C4) 4*(C3) 9*(C2) 4*(S6) 3*(S4) 9*(sigma) "},
    {  "Cinfv", "(Cinf) (sigma) "},
    {  "Dinfh", "(i) (Cinf) (C2) 2*(sigma) "},
    {  "I",     "6*(C5) 10*(C3) 15*(C2) "},
    {  "Ih",    "(i) 6*(C5) 10*(C3) 15*(C2) 6*(S10) 10*(S6) 15*(sigma) "},
    {  "K",     "(Cinf)" },
    {  "Kh",    "(i) (Cinf) (sigma) "},
  } ;
#define PointGroupsCount (sizeof(PointGroups)/sizeof(POINT_GROUP))

  class PointGroupPrivate
  {
  public:
    struct _SYMMETRY_ELEMENT_ {
      void    (*transform_atom)( struct _SYMMETRY_ELEMENT_ *el, OBAtom *from, OBAtom *to ) ;
      int *   transform ;     /*   Correspondence table for the transformation         */
      int     order ;         /*   Applying transformation this many times is identity */
      int     nparam ;        /*   4 for inversion and planes, 7 for axes              */
      double  maxdev ;        /*   Largest error associated with the element            */
      double  distance ;
      double  normal[ DIMENSION ] ;
      double  direction[ DIMENSION ] ;
    };

    typedef _SYMMETRY_ELEMENT_ SYMMETRY_ELEMENT;

    PointGroupPrivate()
    {
      ToleranceSame         = 1e-1;
      TolerancePrimary      = 5e-1;
      ToleranceFinal        = 1e-1;
      MaxOptStep            = 5e-1;
      MinOptStep            = 1e-6;
      GradientStep          = 1e-7;
      OptChangeThreshold    = 1e-10;
      DistanceFromCenter    = NULL;
      verbose               = -1;
      MaxOptCycles          = 500;
      OptChangeHits         = 5;
      MaxAxisOrder          = 20;
      PlanesCount           = 0;
      Planes                = NULL;
      MolecularPlane        = NULL;
      InversionCentersCount = 0;
      InversionCenters      = NULL;
      NormalAxesCount       = 0;
      NormalAxes            = NULL;
      ImproperAxesCount     = 0;
      ImproperAxes          = NULL;
      NormalAxesCounts      = NULL;
      ImproperAxesCounts    = NULL;
      BadOptimization       = 0;
      SymmetryCode          = "";
      PointGroupRejectionReason = NULL ;

      StatTotal             = 0 ;
      StatEarly             = 0 ;
      StatPairs             = 0 ;
      StatDups              = 0 ;
      StatOrder             = 0 ;
      StatOpt               = 0 ;
      StatAccept            = 0 ;

      Setup                 = false;
    }

    OBMol  *               _mol;
    double                 ToleranceSame         ;
    double                 TolerancePrimary      ;
    double                 ToleranceFinal        ;
    double                 MaxOptStep            ;
    double                 MinOptStep            ;
    double                 GradientStep          ;
    double                 OptChangeThreshold    ;
    double                 CenterOfSomething[ DIMENSION ];
    double *               DistanceFromCenter    ;
    int                    verbose               ;
    int                    MaxOptCycles          ;
    int                    OptChangeHits         ;
    int                    MaxAxisOrder          ;
    int                    PlanesCount           ;
    SYMMETRY_ELEMENT **    Planes                ;
    SYMMETRY_ELEMENT *     MolecularPlane        ;
    int                    InversionCentersCount ;
    SYMMETRY_ELEMENT **    InversionCenters      ;
    int                    NormalAxesCount       ;
    SYMMETRY_ELEMENT **    NormalAxes            ;
    int                    ImproperAxesCount     ;
    SYMMETRY_ELEMENT **    ImproperAxes          ;
    int *                  NormalAxesCounts      ;
    int *                  ImproperAxesCounts    ;
    int                    BadOptimization       ;
    const char *           SymmetryCode          ;
    char *                 PointGroupRejectionReason;
    std::vector< std::pair<int, int> > PairedAtoms;
    bool                   Setup;

    /*
     *    Statistics
     */
    long                   StatTotal        ;
    long                   StatEarly        ;
    long                   StatPairs        ;
    long                   StatDups         ;
    long                   StatOrder        ;
    long                   StatOpt          ;
    long                   StatAccept       ;

    bool equivalentAtoms(OBAtom &a1, OBAtom &a2)
    {
      if (a1.GetAtomicNum() != a2.GetAtomicNum())
        return false;
      if (a1.GetIsotope() != a2.GetIsotope())
        return false;
      if (a1.GetFormalCharge() != a2.GetFormalCharge())
        return false;
      if (a1.GetSpinMultiplicity() != a2.GetSpinMultiplicity())
        return false;

      return true;
    }

    int
    establish_pairs( SYMMETRY_ELEMENT *elem )
    {
      unsigned int      i, j, best_j;
      char *            atom_used = (char *)calloc( _mol->NumAtoms(), 1 ) ;
      double            distance, best_distance ;
      OBAtom            symmetric;
      OBAtom            *atom;

      PairedAtoms.clear();

      if( atom_used == NULL ){
        //        fprintf( stderr, "Out of memory for tagging array in establish_pairs()\n" ) ;
        return 0;
      }
      for( i = 0 ; i < _mol->NumAtoms() ; i++ ){
        if( elem->transform[i] >= _mol->NumAtoms() ){ /* No symmetric atom yet          */
          if( verbose > 2 ) printf( "        looking for a pair for %d\n", i ) ;
          elem->transform_atom( elem, _mol->GetAtom(i+1), &symmetric ) ;   // ATOM INDEX ISSUE
          if( verbose > 2 ) printf( "        new coordinates are: (%g,%g,%g)\n",
                                    symmetric.x(), symmetric.y(), symmetric.z() ) ;
          best_j        = i ;
          best_distance = 2*TolerancePrimary ;/* Performance value we'll reject */
          for( j = 0 ; j < _mol->NumAtoms() ; j++ ){

            atom = _mol->GetAtom(j+1);
            // START here
            if( atom_used[j] || !equivalentAtoms(*atom, symmetric) )
              continue ;

            distance = symmetric.GetDistance(atom) ;
            if( verbose > 2 ) printf( "        distance to %d is %g\n", j, distance ) ;
            if( distance < best_distance ){
              best_j        = j ;
              best_distance = distance ;
            }
          }
          if( best_distance > TolerancePrimary ){ /* Too bad, there is no symmetric atom */
            if( verbose > 0 )
              printf( "        no pair for atom %d - best was %d with err = %g\n", i, best_j, best_distance ) ;
            free( atom_used ) ;
            return -1 ;
          }
          elem->transform[i] = best_j ;
          atom_used[best_j]  = 1 ;
          if( verbose > 1 ) printf( "        atom %d transforms to the atom %d, err = %g\n", i, best_j, best_distance ) ;
          std::pair<int, int> atomPair(i, best_j);
          PairedAtoms.push_back( atomPair );
        }
      }
      free( atom_used ) ;
      return 0 ;
    }

    int
    check_transform_order( SYMMETRY_ELEMENT *elem )
    {
      unsigned int i, j, k;

      for( i = 0 ; i < _mol->NumAtoms() ; i++ ){
        if( elem->transform[i] == i )   /* Identity transform is Ok for any order */
          continue ;
        if( elem->transform_atom == rotate_reflect_atom ){
          j = elem->transform[i] ;
          if( elem->transform[j] == i )
            continue ; /* Second-order transform is Ok for improper axis */
        }
        for( j = elem->order - 1, k = elem->transform[i] ; j > 0 ; j--, k = elem->transform[k] ){
          if( k == i ){
            if( verbose > 0 ) printf( "        transform looped %d steps too early from atom %d\n", j, i ) ;
            return -1 ;
          }
        }
        if( k != i && elem->transform_atom == rotate_reflect_atom ){
          /* For improper axes, the complete loop may also take twice the order */
          for( j = elem->order ; j > 0 ; j--, k = elem->transform[k] ){
            if( k == i ){
              if( verbose > 0 ) printf( "        (improper) transform looped %d steps too early from atom %d\n", j, i ) ;
              return -1 ;
            }
          }
        }
        if( k != i ){
          if( verbose > 0 ) printf( "        transform failed to loop after %d steps from atom %d\n", elem->order, i ) ;
          return -1 ;
        }
      }
      return 0 ;
    }

    int
    same_transform( SYMMETRY_ELEMENT *a, SYMMETRY_ELEMENT *b )
    {
      unsigned int      i, j;
      int               code;

      if( ( a->order != b->order ) || ( a->nparam != b->nparam ) || ( a->transform_atom != b->transform_atom ) )
        return 0 ;
      for( i = 0, code = 1 ; i < _mol->NumAtoms() ; i++ ){
        if( a->transform[i] != b->transform[i] ){
          code = 0 ;
          break ;
        }
      }
      if( code == 0 && a->order > 2 ){  /* b can also be a reverse transformation for a */
        for( i = 0 ; i < _mol->NumAtoms() ; i++ ){
          j = a->transform[i] ;
          if( b->transform[j] != i )
            return 0 ;
        }
        return 1 ;
      }
      return code ;
    }

    SYMMETRY_ELEMENT *
    alloc_symmetry_element( void )
    {
      SYMMETRY_ELEMENT * elem = (SYMMETRY_ELEMENT *)calloc( 1, sizeof( SYMMETRY_ELEMENT ) ) ;
      unsigned int i;

      if( elem == NULL ){
        //        fprintf( stderr, "Out of memory allocating symmetry element\n" ) ;
        return NULL;
      }
      elem->transform = (int*)calloc( _mol->NumAtoms(), sizeof( int ) ) ;
      if( elem->transform == NULL ){
        //        fprintf( stderr, "Out of memory allocating transform table for symmetry element\n" ) ;
        free(elem);
        return NULL;
      }
      for( i = 0 ; i < _mol->NumAtoms() ; i++ ){
        elem->transform[i] = _mol->NumAtoms() + 1 ; /* An impossible value */
      }
      return elem ;
    }

    void
    destroy_symmetry_element( SYMMETRY_ELEMENT *elem )
    {
      if( elem != NULL ){
        if( elem->transform != NULL )
          free( elem->transform ) ;
        free( elem ) ;
      }
    }

    int
    check_transform_quality( SYMMETRY_ELEMENT *elem )
    {
      unsigned int      i, j;
      OBAtom            symmetric ;
      double            r, max_r ;

      for( i = 0, max_r = 0 ; i < _mol->NumAtoms() ; i++ ){
        j = elem->transform[i] ;
        elem->transform_atom( elem, _mol->GetAtom(i+1), &symmetric ) ;

        r = symmetric.GetDistance(_mol->GetAtom(j+1));
        if( r > ToleranceFinal ){
          if( verbose > 0 ) printf( "        distance to symmetric atom (%g) is too big for %d\n", r, i ) ;
          return -1 ;
        }
        if( r > max_r ) max_r = r ;
      }
      elem->maxdev = max_r ;
      return 0 ;
    }

    double
    eval_optimization_target_function( SYMMETRY_ELEMENT *elem, int *finish )
    {
      unsigned int      i, j, k;
      OBAtom            symmetric ;
      double            target, r, maxr ;

      if( elem->nparam >= 4 ){
        for( k = 0, r = 0 ; k < DIMENSION ; k++ ){
          r += elem->normal[k]*elem->normal[k] ;
        }
        r = sqrt( r ) ;
        if( r < ToleranceSame ){
          //          fprintf( stderr, "Normal collapsed!\n" ) ;
          return 0.0;
        }
        for( k = 0 ; k < DIMENSION ; k++ ){
          elem->normal[k] /= r ;
        }
        if( elem->distance < 0 ){
          elem->distance = -elem->distance ;
          for( k = 0 ; k < DIMENSION ; k++ ){
            elem->normal[k] = -elem->normal[k] ;
          }
        }
      }
      if( elem->nparam >= 7 ){
        for( k = 0, r = 0 ; k < DIMENSION ; k++ ){
          r += elem->direction[k]*elem->direction[k] ;
        }
        r = sqrt( r ) ;
        if( r < ToleranceSame ){
          //          fprintf( stderr, "Direction collapsed!\n" ) ;
          return 0.0;
        }
        for( k = 0 ; k < DIMENSION ; k++ ){
          elem->direction[k] /= r ;
        }
      }
      for( i = 0, target = maxr = 0 ; i < _mol->NumAtoms() ; i++ ){
        elem->transform_atom( elem, _mol->GetAtom(i+1), &symmetric ) ;
        j = elem->transform[i] ;
        r = symmetric.GetDistance(_mol->GetAtom(j+1));
        if( r > maxr ) maxr = r ;
        target += r ;
      }
      if( finish != NULL ){
        *finish = 0 ;
        if( maxr < ToleranceFinal )
          *finish = 1 ;
      }
      return target ;
    }

    void
    get_params( SYMMETRY_ELEMENT *elem, double values[] )
    {
      memcpy( values, &elem->distance, elem->nparam * sizeof( double ) ) ;
    }

    void
    set_params( SYMMETRY_ELEMENT *elem, double values[] )
    {
      memcpy( &elem->distance, values, elem->nparam * sizeof( double ) ) ;
    }

    void
    optimize_transformation_params( SYMMETRY_ELEMENT *elem )
    {
      double            values[ MAXPARAM ] ;
      double            grad  [ MAXPARAM ] ;
      double            force [ MAXPARAM ] ;
      double            step  [ MAXPARAM ] ;
      double            f, fold, fnew, fnew2, fdn, fup, snorm ;
      double            a, b, x ;
      int               vars  = elem->nparam ;
      int               cycle = 0 ;
      int               i, finish ;
      int               hits = 0 ;

      if( vars > MAXPARAM ){
        //        fprintf( stderr, "Catastrophe in optimize_transformation_params()!\n" ) ;
        return;
      }
      f = 0 ;
      do {
        fold = f ;
        f    = eval_optimization_target_function( elem, &finish ) ;
        /* Evaluate function, gradient and diagonal force constants */
        if( verbose > 1 ) printf( "            function value = %g\n", f ) ;
        if( finish ){
          if( verbose > 1 ) printf( "        function value is small enough\n" ) ;
          break ;
        }
        if( cycle > 0 ){
          if( fabs( f-fold ) > OptChangeThreshold )
            hits = 0 ;
          else hits++ ;
          if( hits >= OptChangeHits ){
            if( verbose > 1 ) printf( "        no progress is made, stop optimization\n" ) ;
            break ;
          }
        }
        get_params( elem, values ) ;
        for( i = 0 ; i < vars ; i++ ){
          values[i] -= GradientStep ;
          set_params( elem, values ) ;
          fdn        = eval_optimization_target_function( elem, NULL ) ;
          values[i] += 2*GradientStep ;
          set_params( elem, values ) ;
          fup        = eval_optimization_target_function( elem, NULL ) ;
          values[i] -= GradientStep ;
          grad[i]    = ( fup - fdn ) / ( 2 * GradientStep ) ;
          force[i]   = ( fup + fdn - 2*f ) / ( GradientStep * GradientStep ) ;
          if( verbose > 1 ) printf( "        i = %d, grad = %12.6e, force = %12.6e\n", i, grad[i], force[i] ) ;
        }
        /* Do a quasi-Newton step */
        for( i = 0, snorm = 0 ; i < vars ; i++ ){
          if( force[i] <  0   ) force[i] = -force[i] ;
          if( force[i] < 1e-3 ) force[i] = 1e-3 ;
          if( force[i] > 1e3  ) force[i] = 1e3 ;
          step[i] = - grad[i]/force[i] ;
          snorm += step[i] * step[i] ;
        }
        snorm = sqrt( snorm ) ;
        if( snorm > MaxOptStep ){ /* Renormalize step */
          for( i = 0 ; i < vars ; i++ )
            step[i] *= MaxOptStep/snorm ;
          snorm = MaxOptStep ;
        }
        do {
          for( i = 0 ; i < vars ; i++ ){
            values[i] += step[i] ;
          }
          set_params( elem, values ) ;
          fnew = eval_optimization_target_function( elem, NULL ) ;
          if( fnew < f )
            break ;
          for( i = 0 ; i < vars ; i++ ){
            values[i] -= step[i] ;
            step  [i] /= 2 ;
          }
          set_params( elem, values ) ;
          snorm /= 2 ;
        } while( snorm > MinOptStep ) ;
        if( (snorm > MinOptStep) && (snorm < MaxOptStep / 2) ){  /* try to do quadratic interpolation */
          for( i = 0 ; i < vars ; i++ )
            values[i] += step[i] ;
          set_params( elem, values ) ;
          fnew2 = eval_optimization_target_function( elem, NULL ) ;
          if( verbose > 1 ) printf( "        interpolation base points: %g, %g, %g\n", f, fnew, fnew2 ) ;
          for( i = 0 ; i < vars ; i++ )
            values[i] -= 2*step[i] ;
          a     = ( 4*f - fnew2 - 3*fnew ) / 2 ;
          b     = ( f + fnew2 - 2*fnew ) / 2 ;
          if( verbose > 1 ) printf( "        linear interpolation coefficients %g, %g\n", a, b ) ;
          if( b > 0 ){
            x = -a/(2*b) ;
            if( x > 0.2 && x < 1.8 ){
              if( verbose > 1 ) printf( "        interpolated: %g\n", x ) ;
              for( i = 0 ; i < vars ; i++ )
                values[i] += x*step[i] ;
            }
            else b = 0 ;
          }
          if( b <= 0 ){
            if( fnew2 < fnew ){
              for( i = 0 ; i < vars ; i++ )
                values[i] += 2*step[i] ;
            }
            else {
              for( i = 0 ; i < vars ; i++ )
                values[i] += step[i] ;
            }
          }
          set_params( elem, values ) ;
        }
      } while( snorm > MinOptStep && ++cycle < MaxOptCycles ) ;
      f = eval_optimization_target_function( elem, NULL ) ;
      if( cycle >= MaxOptCycles ) BadOptimization = 1 ;
      if( verbose > 0 ) {
        if( cycle >= MaxOptCycles )
          printf( "        maximum number of optimization cycles made\n" ) ;
        printf( "        optimization completed after %d cycles with f = %g\n", cycle, f ) ;
      }
    }

    int
    refine_symmetry_element( SYMMETRY_ELEMENT *elem, int build_table )
    {
      int               i ;


      if( build_table && (establish_pairs( elem ) < 0) ){
        StatPairs++ ;
        if( verbose > 0 ) printf( "        no transformation correspondence table can be constructed\n" ) ;
        return -1 ;
      }
      for( i = 0 ; i < PlanesCount ; i++ ){
        if( same_transform( Planes[i], elem ) ){
          StatDups++ ;
          if( verbose > 0 ) printf( "        transformation is identical to plane %d\n", i ) ;
          return -1 ;
        }
      }
      for( i = 0 ; i < InversionCentersCount ; i++ ){
        if( same_transform( InversionCenters[i], elem ) ){
          StatDups++ ;
          if( verbose > 0 ) printf( "        transformation is identical to inversion center %d\n", i ) ;
          return -1 ;
        }
      }
      for( i = 0 ; i < NormalAxesCount ; i++ ){
        if( same_transform( NormalAxes[i], elem ) ){
          StatDups++ ;
          if( verbose > 0 ) printf( "        transformation is identical to normal axis %d\n", i ) ;
          return -1 ;
        }
      }
      for( i = 0 ; i < ImproperAxesCount ; i++ ){
        if( same_transform( ImproperAxes[i], elem ) ){
          StatDups++ ;
          if( verbose > 0 ) printf( "        transformation is identical to improper axis %d\n", i ) ;
          return -1 ;
        }
      }
      if( check_transform_order( elem ) < 0 ){
        StatOrder++ ;
        if( verbose > 0 ) printf( "        incorrect transformation order\n" ) ;
        return -1 ;
      }
      optimize_transformation_params( elem ) ;
      if( check_transform_quality( elem ) < 0 ){
        StatOpt++ ;
        if( verbose > 0 ) printf( "        refined transformation does not pass the numeric threshold\n" ) ;
        return -1 ;
      }
      StatAccept++ ;
      return 0 ;
    }

    /*
     *   Plane-specific functions
     */

    static void mirror_atom( SYMMETRY_ELEMENT *plane, OBAtom *from, OBAtom *to )
    {
      double             r = plane->distance;

      r -= from->x() * plane->normal[0];
      r -= from->y() * plane->normal[1];
      r -= from->z() * plane->normal[2];

      // copy the "type" of from into to
      to->SetAtomicNum(from->GetAtomicNum());
      to->SetIsotope(from->GetIsotope());
      to->SetFormalCharge(from->GetFormalCharge());
      to->SetSpinMultiplicity(from->GetSpinMultiplicity());

      to->SetVector(from->x() + 2*r*plane->normal[0],
                    from->y() + 2*r*plane->normal[1],
                    from->z() + 2*r*plane->normal[2]);
    }

    SYMMETRY_ELEMENT *
    init_mirror_plane( int i, int j )
    {
      SYMMETRY_ELEMENT * plane = alloc_symmetry_element() ;
      double             dx[ DIMENSION ], midpoint[ DIMENSION ], rab, r ;
      int                k ;

      if( verbose > 0 ) printf( "Trying mirror plane for atoms %d,%d\n", i, j ) ;
      StatTotal++ ;
      plane->transform_atom = mirror_atom;
      plane->order          = 2 ;
      plane->nparam         = 4 ;

      dx[0]       = _mol->GetAtom(i+1)->x() - _mol->GetAtom(j+1)->x();
      dx[1]       = _mol->GetAtom(i+1)->y() - _mol->GetAtom(j+1)->y();
      dx[2]       = _mol->GetAtom(i+1)->z() - _mol->GetAtom(j+1)->z();

      midpoint[0] = ( _mol->GetAtom(i+1)->x() + _mol->GetAtom(j+1)->x() ) / 2.0 ;
      midpoint[1] = ( _mol->GetAtom(i+1)->y() + _mol->GetAtom(j+1)->y() ) / 2.0 ;
      midpoint[2] = ( _mol->GetAtom(i+1)->z() + _mol->GetAtom(j+1)->z() ) / 2.0 ;

      rab        = _mol->GetAtom(i+1)->GetDistance(_mol->GetAtom(j+1));

      if( rab < ToleranceSame ){
        //        fprintf( stderr, "Atoms %d and %d coincide (r = %g)\n", i, j, rab ) ;
        destroy_symmetry_element(plane);
        return NULL;
      }
      for( k = 0, r = 0 ; k < DIMENSION ; k++ ){
        plane->normal[k] = dx[k]/rab ;
        r += midpoint[k]*plane->normal[k] ;
      }
      if( r < 0 ){  /* Reverse normal direction, distance is always positive! */
        r = -r ;
        for( k = 0 ; k < DIMENSION ; k++ ){
          plane->normal[k] = -plane->normal[k] ;
        }
      }
      plane->distance = r ;
      if( verbose > 0 ) printf( "    initial plane is at %g from the origin\n", r ) ;
      if( refine_symmetry_element( plane, 1 ) < 0 ){
        if( verbose > 0 ) printf( "    refinement failed for the plane\n" ) ;
        destroy_symmetry_element( plane ) ;
        return NULL ;
      }
      return plane ;
    }

    SYMMETRY_ELEMENT *
    init_ultimate_plane( void )
    {
      SYMMETRY_ELEMENT * plane = alloc_symmetry_element() ;
      double             d0[ DIMENSION ], d1[ DIMENSION ], d2[ DIMENSION ] ;
      double             p[ DIMENSION ] ;
      double             r, s0, s1, s2 ;
      double *           d ;
      unsigned int       i, j, k;

      if( verbose > 0 ) printf( "Trying whole-molecule mirror plane\n" ) ;
      StatTotal++ ;
      plane->transform_atom = mirror_atom ;
      plane->order          = 1 ;
      plane->nparam         = 4 ;
      for( k = 0 ; k < DIMENSION ; k++ )
        d0[k] = d1[k] = d2[k] = 0 ;
      d0[0] = 1 ; d1[1] = 1 ; d2[2] = 1 ;
      for( i = 1 ; i < _mol->NumAtoms() ; i++ ){
        for( j = 0 ; j < i ; j++ ){
          p[0] = _mol->GetAtom(i+1)->x() - _mol->GetAtom(j+1)->x();
          p[1] = _mol->GetAtom(i+1)->y() - _mol->GetAtom(j+1)->y();
          p[2] = _mol->GetAtom(i+1)->z() - _mol->GetAtom(j+1)->z();

          r = sqrt(SQUARE(p[0]) + SQUARE(p[1]) + SQUARE(p[2])); // distance between atoms i and j

          for( k = 0, s0=s1=s2=0 ; k < DIMENSION ; k++ ){
            p[k] /= r ;
            s0   += p[k]*d0[k] ;
            s1   += p[k]*d1[k] ;
            s2   += p[k]*d2[k] ;
          }
          for( k = 0 ; k < DIMENSION ; k++ ){
            d0[k] -= s0*p[k] ;
            d1[k] -= s1*p[k] ;
            d2[k] -= s2*p[k] ;
          }
        }
      }
      for( k = 0, s0=s1=s2=0 ; k < DIMENSION ; k++ ){
        s0 += d0[k] ;
        s1 += d1[k] ;
        s2 += d2[k] ;
      }
      d = NULL ;
      if( s0 >= s1 && s0 >= s2 ) d = d0 ;
      if( s1 >= s0 && s1 >= s2 ) d = d1 ;
      if( s2 >= s0 && s2 >= s1 ) d = d2 ;
      if( d == NULL ){
        fprintf( stderr, "Catastrophe in init_ultimate_plane(): %g, %g and %g have no ordering!\n", s0, s1, s2 ) ;
        destroy_symmetry_element(plane);
        return NULL;
      }
      for( k = 0, r = 0 ; k < DIMENSION ; k++ )
        r += d[k]*d[k] ;
      r = sqrt(r) ;
      if( r > 0 ){
        for( k = 0 ; k < DIMENSION ; k++ )
          plane->normal[k] = d[k]/r ;
      }
      else {
        for( k = 1 ; k < DIMENSION ; k++ )
          plane->normal[k] = 0 ;
        plane->normal[0] = 1 ;
      }
      for( k = 0, r = 0 ; k < DIMENSION ; k++ )
        r += CenterOfSomething[k]*plane->normal[k] ;
      plane->distance = r ;
      for( k = 0 ; k < _mol->NumAtoms() ; k++ )
        plane->transform[k] = k ;
      if( refine_symmetry_element( plane, 0 ) < 0 ){
        if( verbose > 0 ) printf( "    refinement failed for the plane\n" ) ;
        destroy_symmetry_element( plane ) ;
        return NULL ;
      }
      return plane ;
    }

    /*
     *   Inversion-center specific functions
     */
    static void
    invert_atom( SYMMETRY_ELEMENT *center, OBAtom *from, OBAtom *to )
    {

      // copy the "type" of from into to
      to->SetAtomicNum(from->GetAtomicNum());
      to->SetIsotope(from->GetIsotope());
      to->SetFormalCharge(from->GetFormalCharge());
      to->SetSpinMultiplicity(from->GetSpinMultiplicity());

      to->SetVector(2*center->distance*center->normal[0] - from->x(),
                    2*center->distance*center->normal[1] - from->y(),
                    2*center->distance*center->normal[2] - from->z());
    }

    SYMMETRY_ELEMENT *
    init_inversion_center( void )
    {
      SYMMETRY_ELEMENT * center = alloc_symmetry_element() ;
      int                k ;
      double             r ;

      if( verbose > 0 ) printf( "Trying inversion center at the center of something\n" ) ;
      StatTotal++ ;
      center->transform_atom = invert_atom ;
      center->order          = 2 ;
      center->nparam         = 4 ;
      for( k = 0, r = 0 ; k < DIMENSION ; k++ )
        r += CenterOfSomething[k]*CenterOfSomething[k] ;
      r = sqrt(r) ;
      if( r > 0 ){
        for( k = 0 ; k < DIMENSION ; k++ )
          center->normal[k] = CenterOfSomething[k]/r ;
      }
      else {
        center->normal[0] = 1 ;
        for( k = 1 ; k < DIMENSION ; k++ )
          center->normal[k] = 0 ;
      }
      center->distance = r ;
      if( verbose > 0 ) printf( "    initial inversion center is at %g from the origin\n", r ) ;
      if( refine_symmetry_element( center, 1 ) < 0 ){
        if( verbose > 0 ) printf( "    refinement failed for the inversion center\n" ) ;
        destroy_symmetry_element( center ) ;
        return NULL ;
      }
      return center ;
    }

    /*
     *   Normal rotation axis-specific routines.
     */
    static void
    rotate_atom( SYMMETRY_ELEMENT *axis, OBAtom *from, OBAtom *to )
    {
      double             x[3], y[3], a[3], b[3], c[3] ;
      double             angle = axis->order ? 2*M_PI/axis->order : 1.0 ;
      double             a_sin = sin( angle ) ;
      double             a_cos = cos( angle ) ;
      double             dot ;
      int                i ;

      if( DIMENSION != 3 ){
        //        fprintf( stderr, "Catastrophe in rotate_atom!\n" ) ;
        return;
      }

      x[0] = from->x() - axis->distance * axis->normal[0];
      x[1] = from->y() - axis->distance * axis->normal[1];
      x[2] = from->z() - axis->distance * axis->normal[2];

      for( i = 0, dot = 0 ; i < 3 ; i++ )
        dot += x[i] * axis->direction[i] ;
      for( i = 0 ; i < 3 ; i++ )
        a[i] = axis->direction[i] * dot ;
      for( i = 0 ; i < 3 ; i++ )
        b[i] = x[i] - a[i] ;
      c[0] = b[1]*axis->direction[2] - b[2]*axis->direction[1] ;
      c[1] = b[2]*axis->direction[0] - b[0]*axis->direction[2] ;
      c[2] = b[0]*axis->direction[1] - b[1]*axis->direction[0] ;
      for( i = 0 ; i < 3 ; i++ )
        y[i] = a[i] + b[i]*a_cos + c[i]*a_sin ;

      to->SetVector(y[0] + axis->distance * axis->normal[0],
                    y[1] + axis->distance * axis->normal[1],
                    y[2] + axis->distance * axis->normal[2]);

      // copy the "type" of from into to
      to->SetAtomicNum(from->GetAtomicNum());
      to->SetIsotope(from->GetIsotope());
      to->SetFormalCharge(from->GetFormalCharge());
      to->SetSpinMultiplicity(from->GetSpinMultiplicity());
    }

    SYMMETRY_ELEMENT *
    init_ultimate_axis(void)
    {
      SYMMETRY_ELEMENT * axis = alloc_symmetry_element() ;
      double             dir[ DIMENSION ], rel[ DIMENSION ] ;
      double             s ;
      unsigned int       i, k;

      if( verbose > 0 ) printf( "Trying infinity axis\n" ) ;
      StatTotal++ ;
      axis->transform_atom = rotate_atom ;
      axis->order          = 0 ;
      axis->nparam         = 7 ;
      for( k = 0 ; k < DIMENSION ; k++ )
        dir[k] = 0 ;
      for( i = 0 ; i < _mol->NumAtoms() ; i++ ){

        rel[0] = _mol->GetAtom(i+1)->x() - CenterOfSomething[0];
        rel[1] = _mol->GetAtom(i+1)->z() - CenterOfSomething[1];
        rel[2] = _mol->GetAtom(i+1)->y() - CenterOfSomething[2];

        s = rel[0]*dir[0] + rel[1]*dir[1] + rel[2]*dir[2];

        if( s >= 0 )
          for( k = 0 ; k < DIMENSION ; k++ )
            dir[k] += rel[k] ;
        else for( k = 0 ; k < DIMENSION ; k++ )
               dir[k] -= rel[k] ;
      }
      for( k = 0, s = 0 ; k < DIMENSION ; k++ )
        s += SQUARE( dir[k] ) ;
      s = sqrt(s) ;
      if( s > 0 )
        for( k = 0 ; k < DIMENSION ; k++ )
          dir[k] /= s ;
      else dir[0] = 1 ;
      for( k = 0 ; k < DIMENSION ; k++ )
        axis->direction[k] = dir[k] ;
      for( k = 0, s = 0 ; k < DIMENSION ; k++ )
        s += SQUARE( CenterOfSomething[k] ) ;
      s = sqrt(s) ;
      if( s > 0 )
        for( k = 0 ; k < DIMENSION ; k++ )
          axis->normal[k] = CenterOfSomething[k]/s ;
      else {
        for( k = 1 ; k < DIMENSION ; k++ )
          axis->normal[k] = 0 ;
        axis->normal[0] = 1 ;
      }
      axis->distance = s ;
      for( k = 0 ; k < _mol->NumAtoms() ; k++ )
        axis->transform[k] = k ;
      if( refine_symmetry_element( axis, 0 ) < 0 ){
        if( verbose > 0 ) printf( "    refinement failed for the infinity axis\n" ) ;
        destroy_symmetry_element( axis ) ;
        return NULL ;
      }
      return axis ;
    }


    SYMMETRY_ELEMENT *
    init_c2_axis( int i, int j, const double support[ DIMENSION ] )
    {
      SYMMETRY_ELEMENT * axis ;
      int                k ;
      double             ris, rjs ;
      double             r, center[ DIMENSION ] ;

      if( verbose > 0 )
        printf( "Trying c2 axis for the pair (%d,%d) with the support (%g,%g,%g)\n",
                i, j, support[0], support[1], support[2] ) ;
      StatTotal++ ;
      /* First, do a quick sanity check */

      vector3 supportVec(support[0], support[1], support[2]);
      ris = vector3(_mol->GetAtom(i+1)->GetVector() - supportVec).length();
      rjs = vector3(_mol->GetAtom(j+1)->GetVector() - supportVec).length();

      if( fabs( ris - rjs ) > TolerancePrimary ){
        StatEarly++ ;
        if( verbose > 0 ) printf( "    Support can't actually define a rotation axis\n" ) ;
        return NULL ;
      }
      axis                 = alloc_symmetry_element() ;
      axis->transform_atom = rotate_atom ;
      axis->order          = 2 ;
      axis->nparam         = 7 ;
      for( k = 0, r = 0 ; k < DIMENSION ; k++ )
        r += CenterOfSomething[k]*CenterOfSomething[k] ;
      r = sqrt(r) ;
      if( r > 0 ){
        for( k = 0 ; k < DIMENSION ; k++ )
          axis->normal[k] = CenterOfSomething[k]/r ;
      }
      else {
        axis->normal[0] = 1 ;
        for( k = 1 ; k < DIMENSION ; k++ )
          axis->normal[k] = 0 ;
      }
      axis->distance = r ;

      center[0] = ( _mol->GetAtom(i+1)->x() + _mol->GetAtom(j+1)->x() ) / 2 - support[0] ;
      center[1] = ( _mol->GetAtom(i+1)->y() + _mol->GetAtom(j+1)->y() ) / 2 - support[1] ;
      center[2] = ( _mol->GetAtom(i+1)->z() + _mol->GetAtom(j+1)->z() ) / 2 - support[2] ;
      r = sqrt(SQUARE(center[0]) + SQUARE(center[1]) + SQUARE(center[2]));

      if( r <= TolerancePrimary ){ /* c2 is underdefined, let's do something special */
        if( MolecularPlane != NULL ){
          if( verbose > 0 ) printf( "    c2 is underdefined, but there is a molecular plane\n" ) ;
          for( k = 0 ; k < DIMENSION ; k++ )
            axis->direction[k] = MolecularPlane->normal[k] ;
        }
        else {
          if( verbose > 0 ) printf( "    c2 is underdefined, trying random direction\n" ) ;

          center[0] = _mol->GetAtom(i+1)->x() - _mol->GetAtom(j+1)->x();
          center[1] = _mol->GetAtom(i+1)->y() - _mol->GetAtom(j+1)->y();
          center[2] = _mol->GetAtom(i+1)->z() - _mol->GetAtom(j+1)->z();

          if( fabs( center[2] ) + fabs( center[1] ) > ToleranceSame ){
            axis->direction[0] =  0 ;
            axis->direction[1] =  center[2] ;
            axis->direction[2] = -center[1] ;
          }
          else {
            axis->direction[0] = -center[2] ;
            axis->direction[1] =  0 ;
            axis->direction[2] =  center[0] ;
          }
          for( k = 0, r = 0 ; k < DIMENSION ; k++ )
            r += axis->direction[k] * axis->direction[k] ;
          r = sqrt(r) ;
          for( k = 0 ; k < DIMENSION ; k++ )
            axis->direction[k] /= r ;
        }
      }
      else { /* direction is Ok, renormalize it */
        for( k = 0 ; k < DIMENSION ; k++ )
          axis->direction[k] = center[k]/r ;
      }
      if( refine_symmetry_element( axis, 1 ) < 0 ){
        if( verbose > 0 ) printf( "    refinement failed for the c2 axis\n" ) ;
        destroy_symmetry_element( axis ) ;
        return NULL ;
      }
      return axis ;
    }

    SYMMETRY_ELEMENT *
    init_axis_parameters( double a[3], double b[3], double c[3] )
    {
      SYMMETRY_ELEMENT * axis ;
      int                i, order, sign ;
      double             ra, rb, rc, rab, rbc, rac, r ;
      double             angle ;

      ra = rb = rc = rab = rbc = rac = 0 ;
      for( i = 0 ; i < DIMENSION ; i++ ){
        ra  += a[i]*a[i] ;
        rb  += b[i]*b[i] ;
        rc  += c[i]*c[i] ;
      }
      ra = sqrt(ra) ; rb  = sqrt(rb) ; rc  = sqrt(rc) ;
      if( fabs( ra - rb ) > TolerancePrimary || fabs( ra - rc ) > TolerancePrimary || fabs( rb - rc ) > TolerancePrimary ){
        StatEarly++ ;
        if( verbose > 0 ) printf( "    points are not on a sphere\n" ) ;
        return NULL ;
      }
      for( i = 0 ; i < DIMENSION ; i++ ){
        rab += (a[i]-b[i])*(a[i]-b[i]) ;
        rac += (a[i]-c[i])*(a[i]-c[i]) ;
        rbc += (c[i]-b[i])*(c[i]-b[i]) ;
      }
      rab = sqrt(rab) ;
      rac = sqrt(rac) ;
      rbc = sqrt(rbc) ;
      if( fabs( rab - rbc ) > TolerancePrimary ){
        StatEarly++ ;
        if( verbose > 0 ) printf( "    points can't be rotation-equivalent\n" ) ;
        return NULL ;
      }
      if( rab <= ToleranceSame || rbc <= ToleranceSame || rac <= ToleranceSame ){
        StatEarly++ ;
        if( verbose > 0 ) printf( "    rotation is underdefined by these points: %8.3f %8.3f %8.3f\n", rab, rbc, rac ) ;
        return NULL ;
      }
      rab   = (rab+rbc)/2 ;
      angle = M_PI - 2*asin( rac/(2*rab) ) ;
      if( verbose > 1 ) printf( "    rotation angle is %f\n", angle ) ;
      if( fabs(angle) <= M_PI/(MaxAxisOrder+1) ){
        StatEarly++ ;
        if( verbose > 0 ) printf( "    atoms are too close to a straight line\n" ) ;
        return NULL ;
      }
      order = static_cast<int> (floor( (2*M_PI)/angle + 0.5 )) ;
      if( order <= 2 || order > MaxAxisOrder ){
        StatEarly++ ;
        if( verbose > 0 ) printf( "    rotation axis order (%d) is not from 3 to %d\n", order, MaxAxisOrder ) ;
        return NULL ;
      }
      axis = alloc_symmetry_element() ;
      axis->order          = order ;
      axis->nparam         = 7 ;
      for( i = 0, r = 0 ; i < DIMENSION ; i++ )
        r += CenterOfSomething[i]*CenterOfSomething[i] ;
      r = sqrt(r) ;
      if( r > 0 ){
        for( i = 0 ; i < DIMENSION ; i++ )
          axis->normal[i] = CenterOfSomething[i]/r ;
      }
      else {
        axis->normal[0] = 1 ;
        for( i = 1 ; i < DIMENSION ; i++ )
          axis->normal[i] = 0 ;
      }
      axis->distance = r ;
      axis->direction[0] = (b[1]-a[1])*(c[2]-b[2]) - (b[2]-a[2])*(c[1]-b[1]) ;
      axis->direction[1] = (b[2]-a[2])*(c[0]-b[0]) - (b[0]-a[0])*(c[2]-b[2]) ;
      axis->direction[2] = (b[0]-a[0])*(c[1]-b[1]) - (b[1]-a[1])*(c[0]-b[0]) ;
      /*
       *  Arbitrarily select axis direction so that first non-zero component
       *  or the direction is positive.
       */
      sign = 0 ;
      if (axis->direction[0] <= 0) {
        if (axis->direction[0] < 0) {
          sign = 1;
        }
        else if (axis->direction[1] <= 0) {
          if (axis->direction[1] < 0) {
            sign = 1;
          }
          else if (axis->direction[2] < 0) {
            sign = 1;
          }
        }
      }
      if( sign )
        for( i = 0 ; i < DIMENSION ; i++ )
          axis->direction[i] = -axis->direction[i] ;
      for( i = 0, r = 0 ; i < DIMENSION ; i++ )
        r += axis->direction[i]*axis->direction[i] ;
      r = sqrt(r) ;
      for( i = 0 ; i < DIMENSION ; i++ )
        axis->direction[i] /= r ;
      if( verbose > 1 ){
        printf( "    axis origin is at (%g,%g,%g)\n",
                axis->normal[0]*axis->distance, axis->normal[1]*axis->distance, axis->normal[2]*axis->distance ) ;
        printf( "    axis is in the direction (%g,%g,%g)\n", axis->direction[0], axis->direction[1], axis->direction[2] ) ;
      }
      return axis ;
    }

    SYMMETRY_ELEMENT *
    init_higher_axis( int ia, int ib, int ic )
    {
      SYMMETRY_ELEMENT * axis ;
      double             a[ DIMENSION ], b[ DIMENSION ], c[ DIMENSION ] ;

      if( verbose > 0 ) printf( "Trying cn axis for the triplet (%d,%d,%d)\n", ia, ib, ic ) ;
      StatTotal++ ;
      /* Do a quick check of geometry validity */

      a[0] = _mol->GetAtom(ia+1)->x() - CenterOfSomething[0];
      a[1] = _mol->GetAtom(ia+1)->y() - CenterOfSomething[1];
      a[2] = _mol->GetAtom(ia+1)->z() - CenterOfSomething[2];

      b[0] = _mol->GetAtom(ib+1)->x() - CenterOfSomething[0];
      b[1] = _mol->GetAtom(ib+1)->y() - CenterOfSomething[1];
      b[2] = _mol->GetAtom(ib+1)->z() - CenterOfSomething[2];

      c[0] = _mol->GetAtom(ic+1)->x() - CenterOfSomething[0];
      c[1] = _mol->GetAtom(ic+1)->y() - CenterOfSomething[1];
      c[2] = _mol->GetAtom(ic+1)->z() - CenterOfSomething[2];

      if( ( axis = init_axis_parameters( a, b, c ) ) == NULL ){
        if( verbose > 0 ) printf( "    no coherrent axis is defined by the points\n" ) ;
        return NULL ;
      }
      axis->transform_atom = rotate_atom ;
      if( refine_symmetry_element( axis, 1 ) < 0 ){
        if( verbose > 0 ) printf( "    refinement failed for the c%d axis\n", axis->order ) ;
        destroy_symmetry_element( axis ) ;
        return NULL ;
      }
      return axis ;
    }

    /*
     *   Improper axes-specific routines.
     *   These are obtained by slight modifications of normal rotation
     *       routines.
     */
    static void
    rotate_reflect_atom( SYMMETRY_ELEMENT *axis, OBAtom *from, OBAtom *to )
    {
      double             x[3], y[3], a[3], b[3], c[3] ;
      double             angle = 2*M_PI/axis->order ;
      double             a_sin = sin( angle ) ;
      double             a_cos = cos( angle ) ;
      double             dot ;
      int                i ;

      if( DIMENSION != 3 ){
        //        fprintf( stderr, "Catastrophe in rotate_reflect_atom!\n" ) ;
        return;
      }

      x[0] = from->x() - axis->distance * axis->normal[0];
      x[1] = from->y() - axis->distance * axis->normal[1];
      x[2] = from->z() - axis->distance * axis->normal[2];

      for( i = 0, dot = 0 ; i < 3 ; i++ )
        dot += x[i] * axis->direction[i] ;
      for( i = 0 ; i < 3 ; i++ )
        a[i] = axis->direction[i] * dot ;
      for( i = 0 ; i < 3 ; i++ )
        b[i] = x[i] - a[i] ;
      c[0] = b[1]*axis->direction[2] - b[2]*axis->direction[1] ;
      c[1] = b[2]*axis->direction[0] - b[0]*axis->direction[2] ;
      c[2] = b[0]*axis->direction[1] - b[1]*axis->direction[0] ;
      for( i = 0 ; i < 3 ; i++ )
        y[i] = -a[i] + b[i]*a_cos + c[i]*a_sin ;

      to->SetVector(y[0] + axis->distance * axis->normal[0],
                    y[1] + axis->distance * axis->normal[1],
                    y[2] + axis->distance * axis->normal[2]);

      // copy the "type" of from into to
      to->SetAtomicNum(from->GetAtomicNum());
      to->SetIsotope(from->GetIsotope());
      to->SetFormalCharge(from->GetFormalCharge());
      to->SetSpinMultiplicity(from->GetSpinMultiplicity());
    }

    SYMMETRY_ELEMENT *
    init_improper_axis( int ia, int ib, int ic )
    {
      SYMMETRY_ELEMENT * axis ;
      double             a[ DIMENSION ], b[ DIMENSION ], c[ DIMENSION ] ;
      double             centerpoint[ DIMENSION ] ;
      double             r ;
      int                i ;

      if( verbose > 0 ) printf( "Trying an axis for the triplet (%d,%d,%d)\n", ia, ib, ic ) ;
      StatTotal++ ;
      /* First, reduce the problem to Cn case */
      a[0] = _mol->GetAtom(ia+1)->x() - CenterOfSomething[0];
      a[1] = _mol->GetAtom(ia+1)->y() - CenterOfSomething[1];
      a[2] = _mol->GetAtom(ia+1)->z() - CenterOfSomething[2];

      b[0] = _mol->GetAtom(ib+1)->x() - CenterOfSomething[0];
      b[1] = _mol->GetAtom(ib+1)->y() - CenterOfSomething[1];
      b[2] = _mol->GetAtom(ib+1)->z() - CenterOfSomething[2];

      c[0] = _mol->GetAtom(ic+1)->x() - CenterOfSomething[0];
      c[1] = _mol->GetAtom(ic+1)->y() - CenterOfSomething[1];
      c[2] = _mol->GetAtom(ic+1)->z() - CenterOfSomething[2];

      for( i = 0, r = 0 ; i < DIMENSION ; i++ ){
        centerpoint[i] = a[i] + c[i] + 2*b[i] ;
        r             += centerpoint[i]*centerpoint[i] ;
      }
      r = sqrt(r) ;
      if( r <= ToleranceSame ){
        StatEarly++ ;
        if( verbose > 0 ) printf( "    atoms can not define improper axis of the order more than 2\n" ) ;
        return NULL ;
      }
      for( i = 0 ; i < DIMENSION ; i++ )
        centerpoint[i] /= r ;
      for( i = 0, r = 0 ; i < DIMENSION ; i++ )
        r += centerpoint[i] * b[i] ;
      for( i = 0 ; i < DIMENSION ; i++ )
        b[i] = 2*r*centerpoint[i] - b[i] ;
      /* Do a quick check of geometry validity */
      if( ( axis = init_axis_parameters( a, b, c ) ) == NULL ){
        if( verbose > 0 ) printf( "    no coherrent improper axis is defined by the points\n" ) ;
        return NULL ;
      }
      axis->transform_atom = rotate_reflect_atom ;
      if( refine_symmetry_element( axis, 1 ) < 0 ){
        if( verbose > 0 ) printf( "    refinement failed for the s%d axis\n", axis->order ) ;
        destroy_symmetry_element( axis ) ;
        return NULL ;
      }
      return axis ;
    }

    /*
     *   Control routines
     */

    void
    find_center_of_something( void )
    {
      unsigned int       i, j;
      double             coord_sum[ DIMENSION ] ;
      double             r ;
      OBAtom             *atom;

      for( j = 0 ; j < DIMENSION ; j++ )
        coord_sum[j] = 0 ;
      for( i = 0 ; i < _mol->NumAtoms() ; i++ ){
        atom = _mol->GetAtom(i+1);
        coord_sum[0] = atom->x() + atom->y() + atom->z();
      }
      for( j = 0 ; j < DIMENSION ; j++ )
        CenterOfSomething[j] = coord_sum[j]/_mol->NumAtoms() ;
      if( verbose > 0 )
        printf( "Center of something is at %15.10f, %15.10f, %15.10f\n",
                CenterOfSomething[0], CenterOfSomething[1], CenterOfSomething[2] ) ;
      DistanceFromCenter = (double *) calloc( _mol->NumAtoms(), sizeof( double ) ) ;
      if( DistanceFromCenter == NULL ){
        //        fprintf( stderr, "Unable to allocate array for the distances\n" ) ;
        return;
      }
      for( i = 0 ; i < _mol->NumAtoms() ; i++ ){
        atom = _mol->GetAtom(i+1);
        r = SQUARE(atom->x() - CenterOfSomething[0])
          + SQUARE(atom->y() - CenterOfSomething[1])
          + SQUARE(atom->z() - CenterOfSomething[2]);

        DistanceFromCenter[i] = r ;
      }
    }

    void
    find_planes(void)
    {
      unsigned int i, j;
      SYMMETRY_ELEMENT * plane ;

      plane = init_ultimate_plane() ;
      if( plane != NULL ){
        MolecularPlane = plane ;
        PlanesCount++ ;
        Planes = (SYMMETRY_ELEMENT **) realloc( Planes, sizeof( SYMMETRY_ELEMENT* ) * PlanesCount ) ;
        if( Planes == NULL ){
          perror( "Out of memory in find_planes" ) ;
          exit( EXIT_FAILURE ) ;
        }
        Planes[ PlanesCount - 1 ] = plane ;
      }
      for( i = 1 ; i < _mol->NumAtoms() ; i++ ){
        for( j = 0 ; j < i ; j++ ){
          if( !equivalentAtoms(*_mol->GetAtom(i+1), *_mol->GetAtom(j+1)) )
            continue ;
          if( ( plane = init_mirror_plane( i, j ) ) != NULL ){
            PlanesCount++ ;
            Planes = (SYMMETRY_ELEMENT **) realloc( Planes, sizeof( SYMMETRY_ELEMENT* ) * PlanesCount ) ;
            if( Planes == NULL ){
              perror( "Out of memory in find_planes" ) ;
              exit( EXIT_FAILURE ) ;
            }
            Planes[ PlanesCount - 1 ] = plane ;
          }
        }
      }
    }

    void
    find_inversion_centers(void)
    {
      SYMMETRY_ELEMENT * center ;

      if( ( center = init_inversion_center() ) != NULL ){
        InversionCenters = (SYMMETRY_ELEMENT **) calloc( 1, sizeof( SYMMETRY_ELEMENT* ) ) ;
        InversionCenters[0]   = center ;
        InversionCentersCount = 1 ;
      }
    }

    void
    find_infinity_axis(void)
    {
      SYMMETRY_ELEMENT * axis ;

      if( ( axis = init_ultimate_axis() ) != NULL ){
        NormalAxesCount++ ;
        NormalAxes = (SYMMETRY_ELEMENT **) realloc( NormalAxes, sizeof( SYMMETRY_ELEMENT* ) * NormalAxesCount ) ;
        if( NormalAxes == NULL ){
          perror( "Out of memory in find_infinity_axes()" ) ;
          return;
        }
        NormalAxes[ NormalAxesCount - 1 ] = axis ;
      }
    }

    void
    find_c2_axes(void)
    {
      unsigned int       i, j, k, l;
      double             center[ DIMENSION ] ;
      double *           distances = (double*)calloc( _mol->NumAtoms(), sizeof( double ) ) ;
      double             r ;
      SYMMETRY_ELEMENT * axis ;
      OBAtom           *a1, *a2, *a3, *a4;

      if( distances == NULL ){
        //        fprintf( stderr, "Out of memory in find_c2_axes()\n" ) ;
        return;
      }
      for( i = 1 ; i < _mol->NumAtoms() ; i++ ){
        for( j = 0 ; j < i ; j++ ){
          if( !equivalentAtoms(*_mol->GetAtom(i+1), *_mol->GetAtom(j+1)) )
            continue ;
          if( fabs( DistanceFromCenter[i] - DistanceFromCenter[j] ) > TolerancePrimary )
            continue ; /* A very cheap, but quite effective check */
          /*
           *   First, let's try to get it cheap and use CenterOfSomething
           */
          a1 = _mol->GetAtom(i+1);
          a2 = _mol->GetAtom(j+1);
          center[0] = (a1->x() + a2->x()) / 2.0;
          center[1] = (a1->y() + a2->z()) / 2.0;
          center[2] = (a1->z() + a2->y()) / 2.0;

          r = (vector3(center[0], center[1], center[2])
               - vector3(CenterOfSomething[0], CenterOfSomething[1], CenterOfSomething[2])).length();

          if( r > 5*TolerancePrimary ){ /* It's Ok to use CenterOfSomething */
            if( ( axis = init_c2_axis( i, j, CenterOfSomething ) ) != NULL ){
              NormalAxesCount++ ;
              NormalAxes = (SYMMETRY_ELEMENT **) realloc( NormalAxes, sizeof( SYMMETRY_ELEMENT* ) * NormalAxesCount ) ;
              if( NormalAxes == NULL ){
                perror( "Out of memory in find_c2_axes" ) ;
                free(distances);
                return;
              }
              NormalAxes[ NormalAxesCount - 1 ] = axis ;
            }
            continue ;
          }
          /*
           *  Now, C2 axis can either pass through an atom, or through the
           *  middle of the other pair.
           */
          for( k = 0 ; k < _mol->NumAtoms() ; k++ ){
            if( ( axis = init_c2_axis( i, j, _mol->GetAtom(k+1)->GetVector().AsArray() ) ) != NULL ){
              NormalAxesCount++ ;
              NormalAxes = (SYMMETRY_ELEMENT **) realloc( NormalAxes, sizeof( SYMMETRY_ELEMENT* ) * NormalAxesCount ) ;
              if( NormalAxes == NULL ){
                perror( "Out of memory in find_c2_axes" ) ;
                free(distances);
                return;
              }
              NormalAxes[ NormalAxesCount - 1 ] = axis ;
            }
          }
          /*
           *  Prepare data for an additional pre-screening check
           */
          for( k = 0 ; k < _mol->NumAtoms() ; k++ ){
            r = SQUARE(_mol->GetAtom(k+1)->x() - center[0])
              + SQUARE(_mol->GetAtom(k+1)->y() - center[1])
              + SQUARE(_mol->GetAtom(k+1)->z() - center[2]);
            distances[k] = sqrt(r) ;
          }
          for( k = 0 ; k < _mol->NumAtoms() ; k++ ){
            a3 = _mol->GetAtom(k+1);
            for( l = 0 ; l < _mol->NumAtoms() ; l++ ){
              a4 = _mol->GetAtom(l+1);
              if( !equivalentAtoms(*a3, *a4) )
                continue ;
              if( fabs( DistanceFromCenter[k] - DistanceFromCenter[l] ) > TolerancePrimary ||
                  fabs( distances[k] - distances[l] ) > TolerancePrimary )
                continue ; /* We really need this one to run reasonably fast! */

              center[0] = (a3->x() + a4->x()) / 2.0;
              center[1] = (a3->y() + a4->y()) / 2.0;
              center[2] = (a3->z() + a4->z()) / 2.0;

              if( ( axis = init_c2_axis( i, j, center ) ) != NULL ){
                NormalAxesCount++ ;
                NormalAxes = (SYMMETRY_ELEMENT **) realloc( NormalAxes, sizeof( SYMMETRY_ELEMENT* ) * NormalAxesCount ) ;
                if( NormalAxes == NULL ){
                  perror( "Out of memory in find_c2_axes" ) ;
                  free(distances);
                  return;
                }
                NormalAxes[ NormalAxesCount - 1 ] = axis ;
              }
            }
          }
        }
      }
      free( distances ) ;
    }

    void
    find_higher_axes(void)
    {
      unsigned int i, j, k;
      SYMMETRY_ELEMENT * axis ;

      for( i = 0 ; i < _mol->NumAtoms() ; i++ ){
        for( j = i + 1 ; j < _mol->NumAtoms() ; j++ ){
          if( !equivalentAtoms(*_mol->GetAtom(i+1), *_mol->GetAtom(j+1)) )
            continue ;
          if( fabs( DistanceFromCenter[i] - DistanceFromCenter[j] ) > TolerancePrimary )
            continue ; /* A very cheap, but quite effective check */
          for( k = 0 ; k < _mol->NumAtoms() ; k++ ){
            if( !equivalentAtoms(*_mol->GetAtom(i+1), *_mol->GetAtom(k+1)) )
              continue ;
            if( ( fabs( DistanceFromCenter[i] - DistanceFromCenter[k] ) > TolerancePrimary ) ||
                ( fabs( DistanceFromCenter[j] - DistanceFromCenter[k] ) > TolerancePrimary ) )
              continue ;
            if( ( axis = init_higher_axis( i, j, k ) ) != NULL ){
              NormalAxesCount++ ;
              NormalAxes = (SYMMETRY_ELEMENT **) realloc( NormalAxes, sizeof( SYMMETRY_ELEMENT* ) * NormalAxesCount ) ;
              if( NormalAxes == NULL ){
                perror( "Out of memory in find_higher_axes" ) ;
                return;
              }
              NormalAxes[ NormalAxesCount - 1 ] = axis ;
            }
          }
        }
      }
    }

    void
    find_improper_axes(void)
    {
      unsigned int i, j, k;
      SYMMETRY_ELEMENT * axis ;

      for( i = 0 ; i < _mol->NumAtoms() ; i++ ){
        for( j = i + 1 ; j < _mol->NumAtoms() ; j++ ){
          for( k = 0 ; k < _mol->NumAtoms() ; k++ ){
            if( ( axis = init_improper_axis( i, j, k ) ) != NULL ){
              ImproperAxesCount++ ;
              ImproperAxes = (SYMMETRY_ELEMENT **) realloc( ImproperAxes, sizeof( SYMMETRY_ELEMENT* ) * ImproperAxesCount ) ;
              if( ImproperAxes == NULL ){
                perror( "Out of memory in find_improper_axes" ) ;
                return;
              }
              ImproperAxes[ ImproperAxesCount - 1 ] = axis ;
            }
          }
        }
      }
    }

    void
    report_planes( void )
    {
      int           i ;

      if( PlanesCount == 0 )
        printf( "There are no planes of symmetry in the molecule\n" ) ;
      else {
        if( PlanesCount == 1 )
          printf( "There is a plane of symmetry in the molecule\n" ) ;
        else printf( "There are %d planes of symmetry in the molecule\n", PlanesCount ) ;
        printf( "     Residual          Direction of the normal           Distance\n" ) ;
        for( i = 0 ; i < PlanesCount ; i++ ){
          printf( "%3d %8.4e ", i, Planes[i]->maxdev ) ;
          printf( "(%11.8f,%11.8f,%11.8f) ", Planes[i]->normal[0], Planes[i]->normal[1], Planes[i]->normal[2] ) ;
          printf( "%14.8f\n", Planes[i]->distance ) ;
        }
      }
    }

    void
    report_inversion_centers( void )
    {
      if( InversionCentersCount == 0 )
        printf( "There is no inversion center in the molecule\n" ) ;
      else {
        printf( "There in an inversion center in the molecule\n" ) ;
        printf( "     Residual                      Position\n" ) ;
        printf( "   %8.4e ", InversionCenters[0]->maxdev ) ;
        printf( "(%14.8f,%14.8f,%14.8f)\n",
                InversionCenters[0]->distance * InversionCenters[0]->normal[0],
                InversionCenters[0]->distance * InversionCenters[0]->normal[1],
                InversionCenters[0]->distance * InversionCenters[0]->normal[2] ) ;
      }
    }

    void
    report_axes( void )
    {
      int           i ;

      if( NormalAxesCount == 0 )
        printf( "There are no normal axes in the molecule\n" ) ;
      else {
        if( NormalAxesCount == 1 )
          printf( "There is a normal axis in the molecule\n" ) ;
        else printf( "There are %d normal axes in the molecule\n", NormalAxesCount ) ;
        printf( "     Residual  Order         Direction of the axis                         Supporting point\n" ) ;
        for( i = 0 ; i < NormalAxesCount ; i++ ){
          printf( "%3d %8.4e ", i, NormalAxes[i]->maxdev ) ;
          if( NormalAxes[i]->order == 0 )
            printf( "Inf " ) ;
          else printf( "%3d ", NormalAxes[i]->order ) ;
          printf( "(%11.8f,%11.8f,%11.8f) ",
                  NormalAxes[i]->direction[0], NormalAxes[i]->direction[1], NormalAxes[i]->direction[2] ) ;
          printf( "(%14.8f,%14.8f,%14.8f)\n",
                  NormalAxes[0]->distance * NormalAxes[0]->normal[0],
                  NormalAxes[0]->distance * NormalAxes[0]->normal[1],
                  NormalAxes[0]->distance * NormalAxes[0]->normal[2] ) ;
        }
      }
    }

    void
    report_improper_axes( void )
    {
      int           i ;

      if( ImproperAxesCount == 0 )
        printf( "There are no improper axes in the molecule\n" ) ;
      else {
        if( ImproperAxesCount == 1 )
          printf( "There is an improper axis in the molecule\n" ) ;
        else printf( "There are %d improper axes in the molecule\n", ImproperAxesCount ) ;
        printf( "     Residual  Order         Direction of the axis                         Supporting point\n" ) ;
        for( i = 0 ; i < ImproperAxesCount ; i++ ){
          printf( "%3d %8.4e ", i, ImproperAxes[i]->maxdev ) ;
          if( ImproperAxes[i]->order == 0 )
            printf( "Inf " ) ;
          else printf( "%3d ", ImproperAxes[i]->order ) ;
          printf( "(%11.8f,%11.8f,%11.8f) ",
                  ImproperAxes[i]->direction[0], ImproperAxes[i]->direction[1], ImproperAxes[i]->direction[2] ) ;
          printf( "(%14.8f,%14.8f,%14.8f)\n",
                  ImproperAxes[0]->distance * ImproperAxes[0]->normal[0],
                  ImproperAxes[0]->distance * ImproperAxes[0]->normal[1],
                  ImproperAxes[0]->distance * ImproperAxes[0]->normal[2] ) ;
        }
      }
    }

    /*
     *  General symmetry handling
     */
    void
    report_and_reset_counters( void )
    {
      if (verbose > -1)
        printf( "  %10ld candidates examined\n"
                "  %10ld removed early\n"
                "  %10ld removed during initial mating stage\n"
                "  %10ld removed as duplicates\n"
                "  %10ld removed because of the wrong transformation order\n"
                "  %10ld removed after unsuccessful optimization\n"
                "  %10ld accepted\n",
                StatTotal, StatEarly, StatPairs, StatDups, StatOrder, StatOpt, StatAccept ) ;
      StatTotal = StatEarly = StatPairs = StatDups = StatOrder = StatOpt = StatAccept = 0 ;
    }

    void
    find_symmetry_elements( void )
    {
      find_center_of_something() ;
      if( verbose > -1 ){
        printf( "Looking for the inversion center\n" ) ;
      }
      find_inversion_centers() ;
      if( verbose > -1 ){
        report_and_reset_counters() ;
        printf( "Looking for the planes of symmetry\n" ) ;
      }
      find_planes() ;
      if( verbose > -1 ){
        report_and_reset_counters() ;
        printf( "Looking for infinity axis\n" ) ;
      }
      find_infinity_axis() ;
      if( verbose > -1 ){
        report_and_reset_counters() ;
        printf( "Looking for C2 axes\n" ) ;
      }
      find_c2_axes() ;
      if( verbose > -1 ){
        report_and_reset_counters() ;
        printf( "Looking for higher axes\n" ) ;
      }
      find_higher_axes() ;
      if( verbose > -1 ){
        report_and_reset_counters() ;
        printf( "Looking for the improper axes\n" ) ;
      }
      find_improper_axes() ;
      if( verbose > -1 ){
        report_and_reset_counters() ;
      }
    }

    static int
    compare_axes( const void *a, const void *b )
    {
      SYMMETRY_ELEMENT * axis_a = *(SYMMETRY_ELEMENT**) a ;
      SYMMETRY_ELEMENT * axis_b = *(SYMMETRY_ELEMENT**) b ;
      int                i, order_a, order_b ;

      order_a = axis_a->order ; if( order_a == 0 ) order_a = 10000 ;
      order_b = axis_b->order ; if( order_b == 0 ) order_b = 10000 ;
      if( ( i = order_b - order_a ) != 0 ) return i ;
      if( axis_a->maxdev > axis_b->maxdev ) return -1 ;
      if( axis_a->maxdev < axis_b->maxdev ) return  1 ;
      return 0 ;
    }

    void
    sort_symmetry_elements( void )
    {
      if( PlanesCount > 1 ){
        qsort( Planes, PlanesCount, sizeof( SYMMETRY_ELEMENT * ), compare_axes ) ;
      }
      if( NormalAxesCount > 1 ){
        qsort( NormalAxes, NormalAxesCount, sizeof( SYMMETRY_ELEMENT * ), compare_axes ) ;
      }
      if( ImproperAxesCount > 1 ){
        qsort( ImproperAxes, ImproperAxesCount, sizeof( SYMMETRY_ELEMENT * ), compare_axes ) ;
      }
    }

    void
    report_symmetry_elements_verbose( void )
    {
      report_inversion_centers() ;
      report_axes() ;
      report_improper_axes() ;
      report_planes() ;
    }

    void
    summarize_symmetry_elements( void )
    {
      int          i ;

      NormalAxesCounts   = (int*) calloc( MaxAxisOrder+1, sizeof( int ) ) ;
      ImproperAxesCounts = (int*) calloc( MaxAxisOrder+1, sizeof( int ) ) ;
      for( i = 0 ; i < NormalAxesCount ; i++ )
        NormalAxesCounts[ NormalAxes[i]->order ]++ ;
      for( i = 0 ; i < ImproperAxesCount ; i++ )
        ImproperAxesCounts[ ImproperAxes[i]->order ]++ ;
    }

    void
    report_symmetry_elements_brief( void )
    {
      int          i ;
      char *       symmetry_code = (char*)calloc( 1, 10*(PlanesCount+NormalAxesCount+ImproperAxesCount+InversionCentersCount+2) ) ;
      char         buf[ 100 ] ;

      if( symmetry_code == NULL ){
        //        fprintf( stderr, "Unable to allocate memory for symmetry ID code in report_symmetry_elements_brief()\n" ) ;
        return;
      }
      if( PlanesCount + NormalAxesCount + ImproperAxesCount + InversionCentersCount == 0 ) {
        SymmetryCode = symmetry_code ;
        return;
        //        printf( "Molecule has no symmetry elements\n" ) ;
      }
      else {
        //        printf( "Molecule has the following symmetry elements: " ) ;
        if( InversionCentersCount > 0 ) strcat( symmetry_code, "(i) " ) ;
        if( NormalAxesCounts[0] == 1 )
          strcat( symmetry_code, "(Cinf) " ) ;
        if( NormalAxesCounts[0] >  1 ) {
          snprintf( buf, 100, "%d*(Cinf) ", NormalAxesCounts[0] ) ;
          strcat( symmetry_code, buf ) ;
        }
        for( i = MaxAxisOrder ; i >= 2 ; i-- ){
          if( NormalAxesCounts[i] == 1 ){ snprintf( buf, 100, "(C%d) ", i ) ; strcat( symmetry_code, buf ) ; }
          if( NormalAxesCounts[i] >  1 ){ snprintf( buf, 100, "%d*(C%d) ", NormalAxesCounts[i], i ) ; strcat( symmetry_code, buf ) ; }
        }
        for( i = MaxAxisOrder ; i >= 2 ; i-- ){
          if( ImproperAxesCounts[i] == 1 ){ snprintf( buf, 100, "(S%d) ", i ) ; strcat( symmetry_code, buf ) ; }
          if( ImproperAxesCounts[i] >  1 ){ snprintf( buf, 100, "%d*(S%d) ", ImproperAxesCounts[i], i ) ; strcat( symmetry_code, buf ) ; }
        }
        if( PlanesCount == 1 ) strcat( symmetry_code, "(sigma) " ) ;
        if( PlanesCount >  1 ){ snprintf( buf, 100, "%d*(sigma) ", PlanesCount ) ; strcat( symmetry_code, buf ) ; }
        //        printf( "%s\n", symmetry_code ) ;
      }
      SymmetryCode = symmetry_code ;
    }

    const char *identify_point_group( void )
    {
      unsigned int   i;
      int            last_matching = -1;
      int            matching_count = 0;

      for( i = 0 ; i < PointGroupsCount ; i++ ){
        if( strcmp( SymmetryCode, PointGroups[i].symmetry_code ) == 0 ){
          last_matching = i ;
          matching_count++ ;
        }
      }
      if( matching_count == 0 ){
        printf( "These symmetry elements match no point group I know of. Sorry.\n" ) ;
      }
      if( matching_count >  1 ){
        printf( "These symmetry elements match more than one group I know of.\n"
                "SOMETHING IS VERY WRONG\n" ) ;
        printf( "Matching groups are:\n" ) ;
        for( i = 0 ; i < PointGroupsCount ; i++ ){
          if( ( strcmp( SymmetryCode, PointGroups[i].symmetry_code ) == 0 )) {
            printf( "    %s\n", PointGroups[i].group_name ) ;
          }
        }
      }
      if( matching_count == 1 ){
        printf( "It seems to be the %s point group\n", PointGroups[last_matching].group_name ) ;
      }
      return PointGroups[last_matching].group_name;
    }

    void clean_paired_atoms(SYMMETRY_ELEMENT *elem)
    {
      if (PairedAtoms.size() == 0)
        return;

      OBAtom *a, *b, symmetric;
      for (unsigned int idx = 0; idx < PairedAtoms.size(); ++idx) {
          std::pair<int, int> atomPair = PairedAtoms[idx];
          a = _mol->GetAtom(atomPair.first + 1); // ATOM INDEX ISSUE
          b = _mol->GetAtom(atomPair.second + 1);
          elem->transform_atom( elem, a, &symmetric ) ;   // ATOM INDEX ISSUE

          // OK, so symmetric is where b *should* be
          vector3 displacement = b->GetVector() - symmetric.GetVector();
          displacement /= 2.0; // take the average displacement
          a->SetVector(a->GetVector() + displacement);
          b->SetVector(b->GetVector() - displacement);
        }
    }

  }; // end class PointGroupPrivate

  OBPointGroup::OBPointGroup()
  {
    d = new PointGroupPrivate;
  }

  OBPointGroup::~OBPointGroup()
  {
    delete d;
  }

  void OBPointGroup::Setup(OBMol *mol)
  {
    d->_mol = mol;
    d->_mol->Center();
    d->Setup = true;
  }

  //! @todo Remove this on next ABI break
  const char* OBPointGroup::IdentifyPointGroup()
  {
    return this->IdentifyPointGroup(0.01);
  }

  const char* OBPointGroup::IdentifyPointGroup(double tolerance)
  {
    // Don't duplicate work, use the more reliable fallback method
    Symbol pg = IdentifyPointGroupSymbol(tolerance);
    if (pg == Unknown)
      pg = C1; // no known symmetry

    return PointGroups[pg].group_name;
  }

  OBPointGroup::Symbol OBPointGroup::IdentifyPointGroupSymbol(double tolerance)
  {
    d->ToleranceSame = tolerance; // allow for "sloppy" perception
    d->find_symmetry_elements();
    d->sort_symmetry_elements();
    d->summarize_symmetry_elements();
    if ( d->BadOptimization ) {
      // error handling
    }

    d->report_symmetry_elements_brief(); // assign a symmetry code
    //    printf("%s\n", d->SymmetryCode);

    Symbol perceived = Unknown; // assume no symmetry

    if( d->PlanesCount + d->NormalAxesCount + d->ImproperAxesCount + d->InversionCentersCount == 0 )
      return C1; // no symmetry

    // OK, let's use the normal decision tree, falling back to an appropriate sub-group as needed
    // Check for linear or K/Kh groups
    if( d->NormalAxesCounts[0] >= 1 ) // has Cinf axis
      {
        if (d->NormalAxesCounts[2] == 1 && d->PlanesCount > 1 && d->InversionCentersCount == 1) // has C2, so Dinfh
          perceived = Dinfh;
        else if (d->InversionCentersCount == 1 && d->PlanesCount == 1) // no C2, but i = Kh
          perceived = Kh;
        else if (d->PlanesCount >= 1) // fallback
          perceived = Cinfv;
        else // really unlikely
          perceived = K;
      }

    if (d->NormalAxesCounts[5] > 1) { // Possible icosahedral
      if ( strcmp( d->SymmetryCode, PointGroups[Ih].symmetry_code ) == 0 ) {
        perceived = Ih;
      }
      else if ( strcmp( d->SymmetryCode, PointGroups[I].symmetry_code) == 0 ) {
        perceived = I;
      }
      // fall back to a subgroup below
    }
    if (d->NormalAxesCounts[4] > 1) { // Possible octahedral
      if ( strcmp( d->SymmetryCode, PointGroups[Oh].symmetry_code ) == 0 ) {
        perceived = Oh;
      }
      else if ( strcmp( d->SymmetryCode, PointGroups[O].symmetry_code) == 0 ) {
        perceived = O;
      }
      // fall back to a subgroup below
    }
    if (d->NormalAxesCounts[3] > 1) { // Possible tetrahedral
      if ( strcmp( d->SymmetryCode, PointGroups[Th].symmetry_code ) == 0 ) {
        perceived = Th;
      }
      else if ( strcmp( d->SymmetryCode, PointGroups[Td].symmetry_code) == 0 ) {
        perceived = Td;
      }
      else if ( strcmp( d->SymmetryCode, PointGroups[T].symmetry_code) == 0 ) {
        perceived = T;
      }
      // fall back to a subgroup below
    }

    // Find the maximum rotational axis (if present)
    unsigned int maxAxis = 0;
    unsigned int maxImproperAxis = 0;
    bool noAxes = true;
    for (int i = d->MaxAxisOrder; i >= 2; i--) { // Find the maximum axis
      if (d->NormalAxesCounts[i] > 0) {
        maxAxis = i;
        noAxes = false;
        break;
      }
      if (d->ImproperAxesCounts[i] > 0) {
        maxImproperAxis = i;
        noAxes = false;
        // continue to loop so we can find the max regular axis
      }
    }

    if (d->NormalAxesCounts[2] > 1 && maxAxis != 0) { // Probably Dihedral
      if ((maxAxis == 2 && d->NormalAxesCounts[2] == 3)
          || d->NormalAxesCounts[2] >= maxAxis) { // Are there the perpendicular C2 axes?

        // If not, we'll handle it later
        if (d->PlanesCount >= maxAxis + 1) { // Likely Dnh
          switch (maxAxis){
          case 8:
            perceived = D8h;
            break;
          case 7:
            perceived = D7h;
            break;
          case 6:
            perceived = D6h;
            break;
          case 5:
            perceived = D5h;
            break;
          case 4:
            perceived = D4h;
            break;
          case 3:
            perceived = D3h;
            break;
          case 2:
          default:
            perceived = D2h;
          }
        }
        else if (d->PlanesCount == maxAxis) { // Likely Dnd
          switch (maxAxis){
          case 8:
            perceived = D8d;
            break;
          case 7:
            perceived = D7d;
            break;
          case 6:
            perceived = D6d;
            break;
          case 5:
            perceived = D5d;
            break;
          case 4:
            perceived = D4d;
            break;
          case 3:
            perceived = D3d;
            break;
          case 2:
          default:
            perceived = D2d;
          }
        }
        else { // Dn groups
          switch (maxAxis){
          case 8:
            perceived = D8;
            break;
          case 7:
            perceived = D7;
            break;
          case 6:
            perceived = D6;
            break;
          case 5:
            perceived = D5;
            break;
          case 4:
            perceived = D4;
            break;
          case 3:
            perceived = D3;
            break;
          case 2:
          default:
            perceived = D2;
          }
        }
        // Check to see if symmetry code is correct, or we'll try again with a subgroup
        if ( strcmp( d->SymmetryCode, PointGroups[perceived].symmetry_code) != 0 ) {
          perceived = Unknown; // try again
        }
      }
    } // end of dihedral group checks

    if (perceived != Unknown)
      return perceived;

    // OK, we'll try again with subgroups & low symmetry
    if (noAxes) { // low symmetry
      if (d->InversionCentersCount > 0)
        perceived = Ci;
      else if (d->PlanesCount > 0)
        perceived = Cs;
      else
        perceived = C1;

      // Don't look for higher symmetry
      return perceived;
    }

    // Must be Cn? or Sn
    // Check Cnh first, usually largest order
    // Must have sigma-h and either i or Sn axis
    if (d->PlanesCount == 1 && (d->InversionCentersCount > 0 || maxImproperAxis > 0)) { // Cnh
      switch (maxAxis){
      case 8:
        perceived = C8h;
        break;
      case 7:
        perceived = C7h;
        break;
      case 6:
        perceived = C6h;
        break;
      case 5:
        perceived = C5h;
        break;
      case 4:
        perceived = C4h;
        break;
      case 3:
        perceived = C3h;
        break;
      case 2:
      default:
        perceived = C2h;
      }
    }
    // Next Easiest is Cnv => n * sigma
    else if (d->PlanesCount >= maxAxis) { // Cnv
      switch (maxAxis){
      case 8:
        perceived = C8v;
        break;
      case 7:
        perceived = C7v;
        break;
      case 6:
        perceived = C6v;
        break;
      case 5:
        perceived = C5v;
        break;
      case 4:
        perceived = C4v;
        break;
      case 3:
        perceived = C3v;
        break;
      case 2:
      default:
        perceived = C2v;
      }
    }
    else if (maxImproperAxis) {
      switch (maxImproperAxis){
      case 8:
        perceived = S8;
        break;
      case 6:
        perceived = S6;
        break;
      case 4:
        perceived = S4;
        break;
      case 2:
      default:
        perceived = Ci;
      }
    }
    else {// Fallback to Cn
      switch (maxAxis){
      case 8:
        perceived = C8;
        break;
      case 7:
        perceived = C7;
        break;
      case 6:
        perceived = C6;
        break;
      case 5:
        perceived = C5;
        break;
      case 4:
        perceived = C4;
        break;
      case 3:
        perceived = C3;
        break;
      case 2:
      default:
        perceived = C2;
      }
    }

    return perceived;
  }

  void OBPointGroup::Symmetrize(OBMol *mol)
  {
    if (!d->Setup) {
      // TODO: We should also check to see if this mol is different from the original setup molecule
      Setup(mol);
      IdentifyPointGroup(); // make sure we run the symmetry analysis
    }

    // We'll do this in several steps
    // First, inversion centers
    if (d->InversionCentersCount) {
      PointGroupPrivate::SYMMETRY_ELEMENT *center = d->InversionCenters[0];
      d->establish_pairs(center);
      d->clean_paired_atoms(center);
    } // inversion centers

    // Mirror planes
    for (unsigned int i = 0; i < d->PlanesCount; i++)
      {
        d->establish_pairs(d->Planes[i]);
        d->clean_paired_atoms(d->Planes[i]);
      }

    // Proper rotations
    for (unsigned int i = 0; i < d->NormalAxesCount; i++)
      {
        d->establish_pairs(d->NormalAxes[i]);
        d->clean_paired_atoms(d->NormalAxes[i]);
      }

    // Improper rotations
    for (unsigned int i = 0; i < d->ImproperAxesCount; i++)
      {
        d->establish_pairs(d->ImproperAxes[i]);
        d->clean_paired_atoms(d->ImproperAxes[i]);
      }

    // Copy back to the molecule
    OBAtom *atom;
    FOR_ATOMS_OF_MOL(a, d->_mol)
      {
        atom = mol->GetAtom(a->GetIdx());
        atom->SetVector(a->GetVector());
      }
  }

} // end namespace OpenBabel

//! \file pointgroup.cpp
//! \brief Brute-force point group detection

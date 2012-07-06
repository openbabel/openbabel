/*
 * International Chemical Identifier (InChI)
 * Version 1
 * Software version 1.04
 * September 9, 2011
 *
 * The InChI library and programs are free software developed under the
 * auspices of the International Union of Pure and Applied Chemistry (IUPAC).
 * Originally developed at NIST. Modifications and additions by IUPAC 
 * and the InChI Trust.
 *
 * IUPAC/InChI-Trust Licence for the International Chemical Identifier (InChI) 
 * Software version 1.0.
 * Copyright (C) IUPAC and InChI Trust Limited
 * 
 * This library is free software; you can redistribute it and/or modify it under the 
 * terms of the IUPAC/InChI Trust Licence for the International Chemical Identifier 
 * (InChI) Software version 1.0; either version 1.0 of the License, or 
 * (at your option) any later version.
 * 
 * This library is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
 * See the IUPAC/InChI Trust Licence for the International Chemical Identifier (InChI) 
 * Software version 1.0 for more details.
 * 
 * You should have received a copy of the IUPAC/InChI Trust Licence for the 
 * International Chemical Identifier (InChI) Software version 1.0 along with 
 * this library; if not, write to:
 * 
 * The InChI Trust
 * c/o FIZ CHEMIE Berlin
 * Franklinstrasse 11
 * 10587 Berlin
 * GERMANY
 * 
 */


#include <stdio.h>	
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "mode.h"

#include "ichierr.h"
#include "inpdef.h"
#include "extr_ct.h"
#include "ichister.h"
#include "ichiring.h"
#include "ichi.h"

#include "ichicomp.h"
#include "util.h"

#define    ZTYPE_DOWN     (-1)  /*  should be equal to -ZTYPE_UP */
#define    ZTYPE_NONE     0
#define    ZTYPE_UP       1     /*  should be equal to -ZTYPE_DOWN */
#define    ZTYPE_3D       3
#define    ZTYPE_EITHER   9999

/*  criteria for ill-defined */
#define MIN_ANGLE             0.10   /*  5.73 degrees */
#define MIN_SINE              0.03   /*  min edge/plane angle in case the tetrahedra has significantly different edge length */
#define MIN_ANGLE_DBOND       0.087156 /* 5 degrees = max angle considered as too small for unambiguous double bond stereo */
#define MIN_SINE_OUTSIDE      0.06   /*  min edge/plane angle to determine whether the central atom is outside of the tetrahedra */
#define MIN_SINE_SQUARE       0.125  /*  min edge/plane angle in case the tetrahedra is somewhat close to a parallelogram */
#define MIN_SINE_EDGE         0.167  /*  min sine/(min.edge) ratio to avoid undefined in case of long edges */
#define MIN_LEN_STRAIGHT      1.900  /*  min length of two normalized to 1 bonds in a straight line */
#define MAX_SINE              0.70710678118654752440084436210485 /*  1/sqrt(2)=sin(pi/4) */
#define MIN_BOND_LEN          0.000001
#define ZERO_LENGTH           MIN_BOND_LEN
#define ZERO_FLOAT            1.0e-12
#define BOND_PARITY_UNDEFINED 64
#if ( STEREO_CENTER_BONDS_NORM == 1 )
#define MPY_SINE              1.00  /*  was 3.0 */
#define MAX_EDGE_RATIO        2.50   /*  max max/min edge ratio for a tetrahedra close to a parallelogram  */
#else
#define MPY_SINE              3.00
#define MAX_EDGE_RATIO        6.00   /*  max max/min edge ratio for a tetrahedra close to a parallelogram  */
#endif
/*  local prototypes */
static int save_a_stereo_bond( int z_prod, int result_action,
                        int at1, int ord1, AT_NUMB *stereo_bond_neighbor1, S_CHAR *stereo_bond_ord1, S_CHAR *stereo_bond_z_prod1, S_CHAR *stereo_bond_parity1, 
                        int at2, int ord2, AT_NUMB *stereo_bond_neighbor2, S_CHAR *stereo_bond_ord2, S_CHAR *stereo_bond_z_prod2, S_CHAR *stereo_bond_parity2 );
static double get_z_coord( inp_ATOM* at, int cur_atom, int neigh_no,  int *nType,int bPointedEdgeStereo );
static double len3( const double c[] );
static double len2( const double c[] );
static double* diff3( const double a[], const double b[], double result[] );
static double* add3( const double a[], const double b[], double result[] );
static double* mult3( const double a[], double b, double result[] );
static double* copy3( const double a[], double result[] );
static double* change_sign3( const double a[], double result[] );
static double dot_prod3( const double a[], const double b[] );
static int dot_prodchar3( const S_CHAR a[], const S_CHAR b[] );
static double* cross_prod3( const double a[], const double b[], double result[] );
static double triple_prod( double a[], double b[], double c[], double *sine_value );
static double triple_prod_and_min_abs_sine(double at_coord[][3], double *min_sine);
static int are_3_vect_in_one_plane( double at_coord[][3], double min_sine);
static int triple_prod_char( inp_ATOM *at, int at_1, int i_next_at_1, S_CHAR *z_dir1,
                                           int at_2, int i_next_at_2, S_CHAR *z_dir2 );

static int CompDble( const void *a1, const void *a2 );
static int Get2DTetrahedralAmbiguity( double at_coord[][3], int bAddExplicitNeighbor, int bFix2DstereoBorderCase );
static double triple_prod_and_min_abs_sine2(double at_coord[][3], double central_at_coord[], int bAddedExplicitNeighbor, double *min_sine, int *bAmbiguous);
static int are_4at_in_one_plane( double at_coord[][3], double min_sine);
static int bInpAtomHasRequirdNeigh ( inp_ATOM *at, int cur_at, int RequirdNeighType, int NumDbleBonds );
static int bIsSuitableHeteroInpAtom( inp_ATOM  *at );
static int bIsOxide( inp_ATOM  *at, int cur_at );
static int half_stereo_bond_parity( inp_ATOM *at, int cur_at, inp_ATOM *at_removed_H, int num_removed_H, S_CHAR *z_dir, 
                                   int bPointedEdgeStereo, int vABParityUnknown );
static int get_allowed_stereo_bond_type( int bond_type );
static int can_be_a_stereo_bond_with_isotopic_H( inp_ATOM *at, int cur_at, INCHI_MODE nMode );
static int half_stereo_bond_action( int nParity, int bUnknown, int bIsotopic, int vABParityUnknown );
static int set_stereo_bonds_parity( sp_ATOM *out_at, inp_ATOM *at, int at_1, inp_ATOM *at_removed_H, int num_removed_H,
                                   INCHI_MODE nMode, QUEUE *q, AT_RANK *nAtomLevel, 
                                   S_CHAR *cSource, AT_RANK min_sb_ring_size, 
                                   int bPointedEdgeStereo, int vABParityUnknown );
static int can_be_a_stereo_atom_with_isotopic_H( inp_ATOM *at, int cur_at, int bPointedEdgeStereo );
static int set_stereo_atom_parity( sp_ATOM *out_at, inp_ATOM *at, int cur_at, inp_ATOM *at_removed_H, int num_removed_H, 
                                  int bPointedEdgeStereo, int vABParityUnknown );
/*
int set_stereo_parity( inp_ATOM* at, sp_ATOM* at_output, int num_at, int num_removed_H,
                       int *nMaxNumStereoAtoms, int *nMaxNumStereoBonds, INCHI_MODE nMode, int bPointedEdgeStereo, vABParityUnknown );
int get_opposite_sb_atom( inp_ATOM *at, int cur_atom, int icur2nxt, int *pnxt_atom, int *pinxt2cur, int *pinxt_sb_parity_ord );
*/
int ReconcileCmlIncidentBondParities( inp_ATOM *at, int cur_atom, int prev_atom, S_CHAR *visited, int bDisconnected );
int comp_AT_NUMB( const void* a1, const void* a2);
int GetHalfStereobond0DParity( inp_ATOM *at, int cur_at, AT_NUMB nSbNeighOrigAtNumb[], int nNumExplictAttachments, int bond_parity, int nFlag );
int GetStereocenter0DParity( inp_ATOM *at, int cur_at, int j1, AT_NUMB nSbNeighOrigAtNumb[], int nFlag );
int GetSbNeighOrigAtNumb( inp_ATOM *at, int cur_at, inp_ATOM *at_removed_H, int num_removed_H, AT_NUMB nSbNeighOrigAtNumb[]);
int FixSb0DParities( inp_ATOM *at, /* inp_ATOM *at_removed_H, int num_removed_H,*/ int chain_length,
                     int at_1, int i_next_at_1, S_CHAR z_dir1[],
                     int at_2, int i_next_at_2, S_CHAR z_dir2[],
                     int *pparity1, int *pparity2 );

/******************************************************************/


static double         *pDoubleForSort;

/**********************************************************************************/
int comp_AT_NUMB( const void* a1, const void* a2)
{
    return (int)*(const AT_NUMB*)a1 - (int)*(const AT_NUMB*)a2;
}
/******************************************************************/
double get_z_coord( inp_ATOM* at, int cur_atom, int neigh_no,  int *nType, int bPointedEdgeStereo )
{
    int stereo_value = at[cur_atom].bond_stereo[neigh_no];
    int stereo_type  = abs( stereo_value );
    int neigh        = (int)at[cur_atom].neighbor[neigh_no];
    double z         = at[neigh].z - at[cur_atom].z;
    int    bFlat;

    if ( bFlat = (fabs(z) < ZERO_LENGTH) ) {
        int i;
        for ( i = 0; i < at[cur_atom].valence; i ++ ) {
            if ( fabs(at[cur_atom].z - at[(int)at[cur_atom].neighbor[i]].z) > ZERO_LENGTH ) {
                bFlat = 0;
                break;
            }
        }
    }

    if ( bFlat ) {
        if ( !bPointedEdgeStereo || bPointedEdgeStereo * stereo_value >= 0 ) {
            /* bPointedEdgeStereo > 0: define stereo from pointed end of the stereo bond only */
            /* bPointedEdgeStereo < 0: define stereo from wide end of the stereo bond only (case of removed H) */
            switch( stereo_type ) {
                /*  1=Up (solid triangle), 6=Down (Dashed triangle), 4=Either (zigzag triangle) */
            case 0: /*  No stereo */
                *nType = ZTYPE_NONE;
                break;
            case STEREO_SNGL_UP: /*  1= Up */
                *nType = ZTYPE_UP;
                break;
            case STEREO_SNGL_EITHER: /*  4 = Either */
                *nType = ZTYPE_EITHER;
                break;
            case STEREO_SNGL_DOWN: /*  6 = Down */
                *nType = ZTYPE_DOWN;
                break;
            default:
                *nType = ZTYPE_NONE; /*  ignore unexpected values */
            }
            if ( stereo_value < 0 && (*nType == ZTYPE_DOWN || *nType == ZTYPE_UP) )
                *nType = -*nType;
        } else {
            *nType = ZTYPE_NONE; /* no stereo */
        }
    } else
    if ( stereo_type == STEREO_SNGL_EITHER &&
         ( !bPointedEdgeStereo || bPointedEdgeStereo * stereo_value >= 0 ) ) {
        *nType = ZTYPE_EITHER; 
    } else {
        *nType = ZTYPE_3D;
    }
    return z;
}
/******************************************************************/
double len3( const double c[] )
{
    return sqrt( c[0]*c[0]   + c[1]*c[1]   + c[2]*c[2] );
}
/******************************************************************/
double len2( const double c[] )
{
    return sqrt( c[0]*c[0]   + c[1]*c[1] );
}
/******************************************************************/
double* diff3( const double a[], const double b[], double result[] )
{

    result[0] =  a[0] - b[0];
    result[1] =  a[1] - b[1];
    result[2] =  a[2] - b[2];

    return result;
}
/******************************************************************/
double* add3( const double a[], const double b[], double result[] )
{
    result[0] =  a[0] + b[0];
    result[1] =  a[1] + b[1];
    result[2] =  a[2] + b[2];

    return result;
}
/******************************************************************/
double* mult3( const double a[], double b, double result[] )
{
    result[0] = a[0] * b;
    result[1] = a[1] * b;
    result[2] = a[2] * b;

    return result;
}
/*************************************************************/
double* copy3( const double a[], double result[] )
{
    result[0] = a[0];
    result[1] = a[1];
    result[2] = a[2];

    return result;
}
/*************************************************************/
double* change_sign3( const double a[], double result[] )
{
    result[0] = -a[0];
    result[1] = -a[1];
    result[2] = -a[2];

    return result;
}
/*************************************************************/
double dot_prod3( const double a[], const double b[] )
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
/*************************************************************/
int dot_prodchar3( const S_CHAR a[], const S_CHAR b[] )
{
    int prod = ((int)a[0]*(int)b[0] + (int)a[1]*(int)b[1] + (int)a[2]*(int)b[2])/100;
    if ( prod > 100 )
        prod = 100;
    else
    if ( prod < -100 )
        prod = -100;
    return prod;
}
/*************************************************************/
double* cross_prod3( const double a[], const double b[], double result[] )
{
    double tmp[3];
    
    tmp[0] =  (a[1]*b[2]-a[2]*b[1]);
    tmp[1] = -(a[0]*b[2]-a[2]*b[0]);
    tmp[2] =  (a[0]*b[1]-a[1]*b[0]);

    result[0] = tmp[0];
    result[1] = tmp[1];
    result[2] = tmp[2];

    return result;
}
/*************************************************************/
double triple_prod( double a[], double b[], double c[], double *sine_value )
{
    double ab[3], dot_prod_ab_c, abs_c, abs_ab;
    cross_prod3( a, b, ab );
    /* ab[0] =  (a[1]*b[2]-a[2]*b[1]); */
    /* ab[1] = -(a[0]*b[2]-a[2]*b[0]); */
    /* ab[2] =  (a[0]*b[1]-a[1]*b[0]); */
    dot_prod_ab_c   =  dot_prod3( ab, c );
    /* dot_prod_ab_c   =  ab[0]*c[0] + ab[1]*c[1] + ab[2]*c[2]; */
    if ( sine_value ) {
        abs_c  = len3( c );
        /* abs_c  = sqrt( c[0]*c[0]   + c[1]*c[1]   + c[2]*c[2] ); */
        abs_ab = len3( ab );
        /* abs_ab = sqrt( ab[0]*ab[0] + ab[1]*ab[1] + ab[2]*ab[2] ); */

        if ( abs_c > 1.e-7 /* otherwise c has zero length */ && abs_ab > 1.e-7 /* otherwise a is parallel to b*/ ) {
            *sine_value = MPY_SINE * dot_prod_ab_c / ( abs_c * abs_ab);
            /*  *sine_value = dot_prod_ab_c / ( abs_c * abs_ab); */
        } else {
            *sine_value = 0.0;
        }
    }
    return dot_prod_ab_c;
}
/*************************************************************/
int CompDble( const void *a1, const void *a2 )
{
    double diff = pDoubleForSort[*(const int*)a1] - pDoubleForSort[*(const int*)a2];
    if ( diff > 0.0 )
        return 1;
    if ( diff < 0.0 )
        return -1;
    return 0;
}
/*************************************************************/
#define T2D_OKAY  1
#define T2D_WARN  2
#define T2D_UNDF  4
int Get2DTetrahedralAmbiguity( double at_coord[][3], int bAddExplicitNeighbor, int bFix2DstereoBorderCase )
{
/*	const double one_pi = 2.0*atan2(1.0 , 0.0 ); */
const double one_pi = 3.14159265358979323846; /* M_PI */
const double two_pi = 2.0*one_pi;
const double dAngleAndPiMaxDiff = 2.0*atan2(1.0, sqrt(7.0)); /*  min sine between 2 InPlane bonds */
int    nBondType[MAX_NUM_STEREO_ATOM_NEIGH], nBondOrder[MAX_NUM_STEREO_ATOM_NEIGH];
double dBondDirection[MAX_NUM_STEREO_ATOM_NEIGH];
volatile double dAngle, dAlpha, dLimit, dBisector; 
    /* 2010-02-10  added 'volatile': workaround ensuring proper behavior for gcc 32-bit */
    /* cml-enabled compiles at >=O1 for SID484922 and alike (both lin&win had problems) */
int  nNumNeigh = MAX_NUM_STEREO_ATOM_NEIGH - (bAddExplicitNeighbor != 0);
int  i, num_Up, num_Dn, bPrev_Up, cur_len_Up, cur_first_Up, len_Up, first_Up;
int  ret=0;

    for ( i = 0, num_Up = num_Dn = 0; i < nNumNeigh; i ++ ) 
    {
        dAngle = atan2( at_coord[i][1], at_coord[i][0] ); /*  range from -pi to +pi */
        if ( dAngle < 0.0 ) {
            dAngle += two_pi;
        }
        dBondDirection[i] = dAngle;
        nBondType[i] = (at_coord[i][2] > 0.0)? 1 : (at_coord[i][2] < 0.0)? -1 : 0; /* z-coord sign */
        if ( nBondType[i] > 0 ) {
            num_Up ++;
        } else
        if ( nBondType[i] < 0 ) {
            num_Dn ++;
        }
        nBondOrder[i] = i;
    }
    if ( num_Up < num_Dn ) {
        for ( i = 0; i < nNumNeigh; i ++ ) {
            nBondType[i] = -nBondType[i];
        }
        inchi_swap( (char*)&num_Dn, (char*)&num_Up, sizeof(num_Dn) );
    }
    if ( !num_Up ) {
        return T2D_UNDF;
    }

    /*  sort according to the bond orientations */
    pDoubleForSort = dBondDirection;
    insertions_sort( nBondOrder, (unsigned) nNumNeigh, sizeof(nBondOrder[0]), CompDble );
    
    /*  find the longest contiguous sequence of Up bonds */
    if ( num_Up == nNumNeigh ) {
        /*  all bonds are Up */
        len_Up   = cur_len_Up = nNumNeigh; /* added cur_len_Up initialization 1/8/2002 */
        first_Up = 0;
    } else {
        /*  at least one bond is not Up */
        cur_len_Up = len_Up = bPrev_Up = 0;
        /* prev. cycle header version ---
        for ( i = 0; 1; i ++ ) {
            if ( i >= nNumNeigh && !bPrev_Up ) {
                break;
            }
        ----------} */
        /* look at all bonds and continue (circle therough the beginning) as long as the current bond is Up */
        for ( i = 0; i < nNumNeigh || bPrev_Up; i ++ ) {
            if ( nBondType[nBondOrder[i % nNumNeigh]] > 0 ) {
                if ( bPrev_Up ) {
                    cur_len_Up ++; /* uncrement number of Up bonds in current contiguous sequence of them */
                } else {
                    bPrev_Up     = 1; /* start new contiguous sequence of Up bonds */
                    cur_len_Up   = 1;
                    cur_first_Up = i % nNumNeigh;
                }
            } else
            if ( bPrev_Up ) { /* end of contiguous sequence of Up bonds */
                if ( cur_len_Up > len_Up ) {
                    first_Up = cur_first_Up; /* store the sequence because it is longer than the ptrvious one */
                    len_Up   = cur_len_Up;
                }
                bPrev_Up = 0;
            }
        }
    }
#if ( FIX_2D_STEREO_BORDER_CASE == 1 )
    /* check if the bonds with ordering numbers first_Up+len_Up and first_Up+len_Up+1 */
    /* have identical angles. In this case switch their order to enlarge the Up sequence */
#define ZERO_ANGLE  0.000001
    if ( nNumNeigh - len_Up >= 2 ) {
        int next1, next2;
        for ( i = 1; i < nNumNeigh - len_Up; i ++ ) {
            next2 = (first_Up+len_Up + i) % nNumNeigh; /* the 2nd after Up sequence */
            if ( nBondType[nBondOrder[next2]] > 0 ) {
                next1 = (first_Up+len_Up) % nNumNeigh; /* the 1st after Up sequence */
                dAngle = dBondDirection[nBondOrder[next1]] - dBondDirection[nBondOrder[next2]];
                if ( fabs(dAngle) < ZERO_ANGLE ) {
                    inchi_swap( (char*)&nBondOrder[next1], (char*)&nBondOrder[next2], sizeof(nBondOrder[0]) );
                    len_Up ++;
                    break;
                }
            }
        }
    }
    /* check whether the not-Up bond (located before the found first-Up) has */
    /* same angle as the Up bond that precedes this not-Up bond */
    if ( nNumNeigh - len_Up >= 2 ) {
        int next1, next2;
        for ( i = 1; i < nNumNeigh - len_Up; i ++ ) {
            next2 = (first_Up+nNumNeigh - i - 1 ) % nNumNeigh; /* the 2nd before Up sequence */
            if ( nBondType[nBondOrder[next2]] > 0 ) {
                next1 = (first_Up+nNumNeigh-1) % nNumNeigh; /* the 1st before Up sequence */
                dAngle = dBondDirection[nBondOrder[next1]] - dBondDirection[nBondOrder[next2]];
                if ( fabs(dAngle) < ZERO_ANGLE ) {
                    inchi_swap( (char*)&nBondOrder[next1], (char*)&nBondOrder[next2], sizeof(nBondOrder[0]) );
                    first_Up = next1; 
                    len_Up ++;
                    break;
                }
            }
        }
    }
#else
    if ( bFix2DstereoBorderCase ) {
        /* check if the bonds with ordering numbers first_Up+len_Up and first_Up+len_Up+1 */
        /* have identical angles. In this case switch their order to enlarge the Up sequence */
#define ZERO_ANGLE  0.000001
        if ( nNumNeigh - len_Up >= 2 ) {
            int next1, next2;
            for ( i = 1; i < nNumNeigh - len_Up; i ++ ) {
                next2 = (first_Up+len_Up + i) % nNumNeigh; /* the 2nd after Up sequence */
                if ( nBondType[nBondOrder[next2]] > 0 ) {
                    next1 = (first_Up+len_Up) % nNumNeigh; /* the 1st after Up sequence */
                    dAngle = dBondDirection[nBondOrder[next1]] - dBondDirection[nBondOrder[next2]];
                    if ( fabs(dAngle) < ZERO_ANGLE ) {
                        inchi_swap( (char*)&nBondOrder[next1], (char*)&nBondOrder[next2], sizeof(nBondOrder[0]) );
                        len_Up ++;
                        break;
                    }
                }
            }
        }
        /* check whether the not-Up bond (located before the found first-Up) has */
        /* same angle as the Up bond that precedes this not-Up bond */
        if ( nNumNeigh - len_Up >= 2 ) {
            int next1, next2;
            for ( i = 1; i < nNumNeigh - len_Up; i ++ ) {
                next2 = (first_Up+nNumNeigh - i - 1 ) % nNumNeigh; /* the 2nd before Up sequence */
                if ( nBondType[nBondOrder[next2]] > 0 ) {
                    next1 = (first_Up+nNumNeigh-1) % nNumNeigh; /* the 1st before Up sequence */
                    dAngle = dBondDirection[nBondOrder[next1]] - dBondDirection[nBondOrder[next2]];
                    if ( fabs(dAngle) < ZERO_ANGLE ) {
                        inchi_swap( (char*)&nBondOrder[next1], (char*)&nBondOrder[next2], sizeof(nBondOrder[0]) );
                        first_Up = next1; 
                        len_Up ++;
                        break;
                    }
                }
            }
        }
    }
#endif
    /*  Turn all the bonds around the center so that */
    /*  the 1st Up bond has zero radian direction */
    dAlpha = dBondDirection[nBondOrder[first_Up]];
    for ( i = 0; i < nNumNeigh; i ++ ) {
        if ( i == nBondOrder[first_Up] ) {
            dBondDirection[i] = 0.0;
        } else {
            dAngle = dBondDirection[i] - dAlpha;
            if ( dAngle < 0.0 ) {
                dAngle += two_pi;
            }
            dBondDirection[i] = dAngle;
        }
    }

    /********************************************************
     * Process particular cases
     ********************************************************/


    if ( nNumNeigh == 3 ) /************************ 3 bonds ************************/
    {
        
        switch( num_Up ) 
        {
            

            
            case 0:	/* 0 Up */
                    return T2D_UNDF;
                
            

            
            case 1:	/* 1 Up */
                    if ( num_Dn ) 
                    {
#ifdef _DEBUG
                        if ( num_Dn != 1 )  /*  debug only */
                            return -1;
#endif
                        ret = (T2D_UNDF | T2D_WARN);
                    } 
                    else 
                    {
                        dAngle = dBondDirection[nBondOrder[(first_Up + 2) % nNumNeigh]] -
                                 dBondDirection[nBondOrder[(first_Up + 1) % nNumNeigh]];

                        if ( dAngle < 0.0 ) 
                            dAngle += two_pi;
                        if ( dAngle - one_pi < -MIN_ANGLE || dAngle - one_pi > MIN_ANGLE  ) 
                        {
                            ret = T2D_OKAY;
                        }
                        else 
                        {
                            ret = (T2D_UNDF | T2D_WARN);
                        }
                    }
                    break;
        

            
            
            case 2:	/* 2 Up */
                    if ( num_Dn ) 
                    {
                        dAlpha = dBondDirection[nBondOrder[(first_Up + 1) % nNumNeigh]] -
                                 dBondDirection[nBondOrder[(first_Up    ) % nNumNeigh]];

                        if ( dAlpha < 0.0 ) 
                            dAlpha += two_pi;

                        if ( dAlpha > one_pi - MIN_ANGLE ) 
                        {
                            ret = T2D_OKAY;
                        } 
                        else if ( dAlpha < two_pi / 3.0 - MIN_ANGLE ) 
                        {
                            ret = (T2D_UNDF | T2D_WARN);
                        } 
                        else 
                        {
                            /*  angle between 2 Up bonds is between 120 and 180 degrees */
                            /*  direction of the (Alpha angle bisector) + 180 degrees	*/
                            dBisector = dBondDirection[nBondOrder[(first_Up    ) % nNumNeigh]];
                            dBisector+= dBondDirection[nBondOrder[(first_Up + 1 ) % nNumNeigh]];
                            dBisector/= 2.0;
                            dBisector-= one_pi;
                            if ( dBisector < 0.0 ) 
                            {
                                dBisector += two_pi;
                            }
                            if ( dAlpha < two_pi / 3.0 + MIN_ANGLE ) 
                            {
                                /*  dAlpha is inside ( 2pi/3 - eps, 2pi/3 + eps ) interval */
                                dLimit = MIN_ANGLE * 3.0 / 2.0;
                            } 
                            else 
                            {
                                dLimit = dAlpha * 3.0 / 2.0 - one_pi;
                            }

                            dAngle = dBondDirection[nBondOrder[(first_Up + 2 ) % nNumNeigh]];

                            if ( dBisector - dAngle < -dLimit ||
                                  dBisector - dAngle >  dLimit  ) 
                            {
                                ret = (T2D_UNDF | T2D_WARN);
                            } 
                            else 
                            {
                                ret = T2D_OKAY;
                            }
                        }
                    } /* if ( num_Dn )  */
                    else 
                    {
                        ret = T2D_OKAY;
                    }
                    break;
        
                    
                    
            case 3:	/* 3 Up */
                    ret = T2D_OKAY;
                    break;
        
                    
            default:/* other Up */
                    return -1;

        } /* eof switch( num_Up ) at  nNumNeigh == 3 */

    } 


    else if ( nNumNeigh == 4) /******************************* 4 bonds ********************/
    {		
        switch( num_Up ) 
        {
        
            case 0:	/* 0 Up */
                    return T2D_UNDF;

        
            case 1:	/* 1 Up */
                    if ( num_Dn ) 
                    {
                        if ( nBondType[nBondOrder[(first_Up + 2) % nNumNeigh]] < 0 ) 
                        {
                            /*
                            * Up, In Plane, Dn, In Plane. Undefined if angle between
                            * two In Plane bonds is wuthin pi +/- 2*arcsine(1/sqrt(8)) interval
                            * That is, 138.5 to 221.4 degrees; for certainty the interval is
                            * increased by 5.7 degrees at each end to
                            * 134.8 to 227.1 degrees
                            */
                            dAngle = dBondDirection[nBondOrder[(first_Up + 3) % nNumNeigh]] -
                                     dBondDirection[nBondOrder[(first_Up + 1) % nNumNeigh]];
                            if ( dAngle < 0.0 ) {
                                dAngle += two_pi;
                            }
                            if ( fabs( dAngle - one_pi ) < dAngleAndPiMaxDiff + MIN_ANGLE ) {
                                ret = (T2D_UNDF | T2D_WARN); 
                            } 
                            else 
                            {
                                ret = T2D_OKAY;
                            }
                        } 
                        else 
                        {
                            ret = T2D_OKAY;
                        }
#ifdef _DEBUG
                        if ( num_Dn != 1 )  /*  debug only */
                        return -1;
#endif
                    } 
                    else 
                    {
                        ret    = T2D_OKAY;
                        dAngle = dBondDirection[nBondOrder[(first_Up + 3) % nNumNeigh]] -
                                 dBondDirection[nBondOrder[(first_Up + 1) % nNumNeigh]];
                        if ( dAngle < 0.0 ) 
                        {
                            dAngle += two_pi;
                        }
                        if ( dAngle < one_pi - MIN_ANGLE ) 
                        {
                            ret |= T2D_WARN;
                        }
                    }
                    break;

        
            case 2:	/* 2 Up */
#if ( FIX_2D_STEREO_BORDER_CASE == 1 )
                    if ( len_Up == 1 ) 
                    {
                        ret = T2D_OKAY;
                    } 
                    else 
                    {
                        dAngle = dBondDirection[nBondOrder[(first_Up + 3) % nNumNeigh]] -
                                 dBondDirection[nBondOrder[(first_Up + 0) % nNumNeigh]];
                        dAngle = fabs(two_pi - dAngle);
                        dAlpha = dBondDirection[nBondOrder[(first_Up + 2) % nNumNeigh]] -
                                 dBondDirection[nBondOrder[(first_Up + 1) % nNumNeigh]];
                        dAlpha = fabs(dAlpha);
                        if ( dAngle < 2.0 * ZERO_ANGLE && dAlpha > MIN_ANGLE ||
                             dAlpha < 2.0 * ZERO_ANGLE && dAngle > MIN_ANGLE  ) 
                        {
                            ret = (T2D_OKAY | T2D_WARN);
                        } 
                        else 
                        {
                            ret = (T2D_UNDF | T2D_WARN);
                        }
                    }
#else
                    if ( bFix2DstereoBorderCase ) 
                    {
                        /* bug fix */
                        if ( len_Up == 1 ) 
                        {
                            ret = T2D_OKAY;
                        } 
                        else 
                        {
                            dAngle = dBondDirection[nBondOrder[(first_Up + 3) % nNumNeigh]] -
                                     dBondDirection[nBondOrder[(first_Up + 0) % nNumNeigh]];
                            dAngle = fabs(two_pi - dAngle);
                            dAlpha = dBondDirection[nBondOrder[(first_Up + 2) % nNumNeigh]] -
                                     dBondDirection[nBondOrder[(first_Up + 1) % nNumNeigh]];
                            dAlpha = fabs(dAlpha);
                            if ( dAngle < 2.0 * ZERO_ANGLE && dAlpha > MIN_ANGLE ||
                                 dAlpha < 2.0 * ZERO_ANGLE && dAngle > MIN_ANGLE  ) 
                            {
                                ret = (T2D_OKAY | T2D_WARN);
                            } 
                            else 
                            {
                            ret = (T2D_UNDF | T2D_WARN);
                            }
                        }
                    } 
                    else 
                    {
                        /* original InChI v. 1 bug */
                        if ( cur_len_Up == 1 ) 
                        {
                            ret = T2D_OKAY;
                        } 
                        else 
                        {
                            ret = (T2D_UNDF | T2D_WARN);
                        }
                    }
#endif
                    break;
        
        
            case 3:	/* 3 Up */
                    ret    = T2D_OKAY;
                    dAngle = dBondDirection[nBondOrder[(first_Up + 2) % nNumNeigh]] -
                             dBondDirection[nBondOrder[(first_Up + 0) % nNumNeigh]];
                    if ( dAngle < 0.0 ) 
                    {
                        dAngle += two_pi;
                    }
                    if ( dAngle < one_pi - MIN_ANGLE ) 
                    {
                        ret |= T2D_WARN;
                    }
                    break;
        
            case 4:	/* 4 Up */
                    ret = (T2D_UNDF | T2D_WARN);
                    break;
        
            default:/* other Up */
                    return -1; /*  program error */
        
        } /* eof switch( num_Up ) at  nNumNeigh == 4 */

        if ( ret == T2D_OKAY ) 
        {
            /*  check whether all bonds are inside a less than 180 degrees sector */
            for ( i = 0; i < nNumNeigh; i ++ ) 
            {
                dAngle = dBondDirection[nBondOrder[(i + nNumNeigh - 1) % nNumNeigh]] -
                         dBondDirection[nBondOrder[ i % nNumNeigh]];
                if ( dAngle < 0.0 ) 
                {
                    dAngle += two_pi;
                }
                if ( dAngle < one_pi - MIN_ANGLE ) 
                {
                    ret |= T2D_WARN;
                    break;
                }
            }
        }

    } /* eof nNumNeigh == 4 */
    
    else /*************************** number of bonds != 3 or 4 ******************/	
    {
        
            return -1; /*  error */   
    } 


    return ret;

}

/*************************************************************/
double triple_prod_and_min_abs_sine2(double at_coord[][3], double central_at_coord[], int bAddedExplicitNeighbor, double *min_sine, int *bAmbiguous)
{
    double min_sine_value=9999.0, sine_value, min_edge_len, max_edge_len, min_edge_len_NoExplNeigh, max_edge_len_NoExplNeigh;
    double s0, s1, s2, s3, e01, e02, e03, e12, e13, e23, tmp[3], e[3][3];
    double prod, ret, central_prod[4];
    int    bLongEdges;

    if ( !min_sine ) {
        return triple_prod( at_coord[0], at_coord[1], at_coord[2], NULL );
    }
    
    ret = triple_prod( at_coord[0], at_coord[1], at_coord[2], &sine_value );
    sine_value = MPY_SINE * fabs( sine_value );
    
    diff3( at_coord[1], at_coord[0], e[2] );
    diff3( at_coord[0], at_coord[2], e[1] );
    diff3( at_coord[2], at_coord[1], e[0] );
    
    /*  lengths of the 6 edges of the tetrahedra */
    e03 = len3( at_coord[0] ); /* 1 */
    e13 = len3( at_coord[1] );
    e23 = len3( at_coord[2] ); /* includes added neighbor if bAddedExplicitNeighbor*/
    e02 = len3( e[1] );        /* includes added neighbor if bAddedExplicitNeighbor*/
    e12 = len3( e[0] );        /* includes added neighbor if bAddedExplicitNeighbor*/
    e01 = len3( e[2] );        
    
    /*  min & max edge length */
    max_edge_len =
    min_edge_len = e03;

    if ( min_edge_len > e13 )
        min_edge_len = e13;
    if ( min_edge_len > e01 )
        min_edge_len = e01;
    min_edge_len_NoExplNeigh = min_edge_len;

    if ( min_edge_len > e23 )
        min_edge_len = e23;
    if ( min_edge_len > e02 )
        min_edge_len = e02;
    if ( min_edge_len > e12 )
        min_edge_len = e12;

    if ( max_edge_len < e13 )
        max_edge_len = e13;
    if ( max_edge_len < e01 )
        max_edge_len = e01;
    max_edge_len_NoExplNeigh = max_edge_len;

    if ( max_edge_len < e23 )
        max_edge_len = e23;
    if ( max_edge_len < e02 )
        max_edge_len = e02;
    if ( max_edge_len < e12 )
        max_edge_len = e12;
    
    if ( !bAddedExplicitNeighbor ) {
        min_edge_len_NoExplNeigh = min_edge_len;
        max_edge_len_NoExplNeigh = max_edge_len;
    }

    bLongEdges = bAddedExplicitNeighbor? 
                  ( max_edge_len_NoExplNeigh < MAX_EDGE_RATIO * min_edge_len_NoExplNeigh ) :
                  ( max_edge_len             < MAX_EDGE_RATIO * min_edge_len );

    if ( sine_value > MIN_SINE && ( min_sine || bAmbiguous ) ) {
        if ( min_sine ) {
            prod = fabs( ret );
            /*  tetrahedra height = volume(prod) / area of a plane(cross_prod) */
            /*  (instead of a tetrahedra calculate parallelogram/parallelepiped area/volume) */

            /*  4 heights from each of the 4 vertices to the opposite plane */
            s0  = prod / len3( cross_prod3( at_coord[1], at_coord[2], tmp ) );
            s1  = prod / len3( cross_prod3( at_coord[0], at_coord[2], tmp ) );
            s2  = prod / len3( cross_prod3( at_coord[0], at_coord[1], tmp ) );
            s3  = prod / len3( cross_prod3( e[0], e[1], tmp ) );
            /*  abs. value of a sine of an angle between each tetrahedra edge and plane */
            /*  sine = height / edge length */
            if ( (sine_value = s0/e01) < min_sine_value )
                min_sine_value = sine_value;
            if ( (sine_value = s0/e02) < min_sine_value )
                min_sine_value = sine_value;
            if ( (sine_value = s0/e03) < min_sine_value )
                min_sine_value = sine_value;

            if ( (sine_value = s1/e01) < min_sine_value )
                min_sine_value = sine_value;
            if ( (sine_value = s1/e12) < min_sine_value )
                min_sine_value = sine_value;
            if ( (sine_value = s1/e13) < min_sine_value )
                min_sine_value = sine_value;

            if ( (sine_value = s2/e02) < min_sine_value )
                min_sine_value = sine_value;
            if ( (sine_value = s2/e12) < min_sine_value )
                min_sine_value = sine_value;
            if ( (sine_value = s2/e23) < min_sine_value )
                min_sine_value = sine_value;

            if ( (sine_value = s3/e03) < min_sine_value )
                min_sine_value = sine_value;
            if ( (sine_value = s3/e13) < min_sine_value )
                min_sine_value = sine_value;
            if ( (sine_value = s3/e23) < min_sine_value )
                min_sine_value = sine_value;
            /*  actually use triple sine */
            *min_sine = sine_value = MPY_SINE * min_sine_value;
        }

        if ( bAmbiguous && sine_value >= MIN_SINE ) {
            /*  check whether the central atom is outside the tetrahedra (0,0,0), at_coord[0,1,2] */
            /*  compare the tetrahedra volume and the volume of a tetrahedra having central_at_coord[] vertex */
            int i;
            diff3( central_at_coord, at_coord[0], tmp );
            central_prod[0] = triple_prod( at_coord[0], at_coord[1], central_at_coord, NULL );
            central_prod[1] = triple_prod( at_coord[1], at_coord[2], central_at_coord, NULL );
            central_prod[2] = triple_prod( at_coord[2], at_coord[0], central_at_coord, NULL );
            central_prod[3] = triple_prod( e[2], e[1], tmp, NULL );
            for ( i = 0; i <= 3; i ++ ) {
                if ( central_prod[i] / ret < -MIN_SINE_OUTSIDE ) {
                    *bAmbiguous |= AMBIGUOUS_STEREO;
                    break;
                }
            }
        }
#if ( STEREO_CENTER_BONDS_NORM == 1 )        
        
        if ( bLongEdges && !bAddedExplicitNeighbor && max_edge_len >= MIN_LEN_STRAIGHT ) {
            /*  possible planar tetragon */
            if ( sine_value < MIN_SINE_SQUARE ) {
                *min_sine = MIN_SINE / 2.0; /*  force parity to be undefined */
                if ( bAmbiguous && !*bAmbiguous ) {
                    *bAmbiguous |= AMBIGUOUS_STEREO;
                }
            }
        }
        
        if ( bLongEdges && sine_value < MIN_SINE_SQUARE && sine_value < MIN_SINE_EDGE * min_edge_len_NoExplNeigh ) {
            *min_sine = MIN_SINE / 2.0; /*  force parity to be undefined */
            if ( bAmbiguous && !*bAmbiguous ) {
                *bAmbiguous |= AMBIGUOUS_STEREO;
            }
        }
#endif

    } else
    if ( min_sine ) {
        *min_sine = sine_value;
    }
    
    return ret;
}
/*************************************************************/
double triple_prod_and_min_abs_sine(double at_coord[][3], double *min_sine)
{
    double min_sine_value=9999.0, sine_value;
    double prod=0.0;

    if ( !min_sine ) {
        return triple_prod( at_coord[0], at_coord[1], at_coord[2], NULL );
    }
    
    prod = triple_prod( at_coord[0], at_coord[1], at_coord[2], &sine_value );
    sine_value = fabs( sine_value );
    min_sine_value = inchi_min( min_sine_value, sine_value );
    
    prod = triple_prod( at_coord[1], at_coord[2], at_coord[0], &sine_value );
    sine_value = fabs( sine_value );
    min_sine_value = inchi_min( min_sine_value, sine_value );
    
    prod = triple_prod( at_coord[2], at_coord[0], at_coord[1], &sine_value );
    sine_value = fabs( sine_value );
    min_sine_value = inchi_min( min_sine_value, sine_value );

    *min_sine = min_sine_value;
    
    return prod;
}
/*************************************************************/
/*  Find if point (0,0,0)a and 3 atoms are in one plane */
int are_3_vect_in_one_plane( double at_coord[][3], double min_sine)
{
    double actual_min_sine;
    double prod;
    prod = triple_prod_and_min_abs_sine( at_coord, &actual_min_sine);
    return actual_min_sine <= min_sine;
}
/*************************************************************/
/*  Find if 4 atoms are in one plane */
int are_4at_in_one_plane( double at_coord[][3], double min_sine)
{
    double actual_min_sine, min_actual_min_sine;
    double coord[3][3], prod;
    int i, k, j;
    for ( k = 0; k < 4; k ++ ) { /* cycle added 4004-08-15 */
        for ( i = j = 0; i < 4; i ++ ) {
            if ( i != k ) {
                diff3( at_coord[i], at_coord[k], coord[j] );
                j ++;
            }
        }
        prod = triple_prod_and_min_abs_sine( coord, &actual_min_sine);
        if ( !k || actual_min_sine < min_actual_min_sine ) {
            min_actual_min_sine = actual_min_sine;
        }
    }
    return min_actual_min_sine <= min_sine;
}
/*************************************************************/
int triple_prod_char( inp_ATOM *at, int at_1, int i_next_at_1, S_CHAR *z_dir1,
                                    int at_2, int i_next_at_2, S_CHAR *z_dir2 )
{
    inp_ATOM *at1, *at2;
    double    pnt[3][3], len;
    int       i;
    int       ret = 0;

    at1 = at + at_1;
    at2 = at + at[at_1].neighbor[i_next_at_1];

    pnt[0][0] = at2->x - at1->x;
    pnt[0][1] = at2->y - at1->y;
    pnt[0][2] = at2->z - at1->z;

    at2 = at + at_2;
    at1 = at + at[at_2].neighbor[i_next_at_2];

    pnt[1][0] = at2->x - at1->x; 
    pnt[1][1] = at2->y - at1->y; 
    pnt[1][2] = at2->z - at1->z;
/*
 *  resultant pnt vector directions:
 *
 *         pnt[0]              pnt[1]
 *
 *   [at_1]---->[...]    [...]---->[at_2]
 *
 *  
 *  add3 below: (pnt[0] + pnt[1]) -> pnt[1]
 */
    add3( pnt[0], pnt[1], pnt[1] );



    for ( i = 0; i < 3; i ++ ) {
        pnt[0][i] = (double)z_dir1[i];
        pnt[2][i] = (double)z_dir2[i];
    }
    for ( i = 0; i < 3; i ++ ) {
        len = len3( pnt[i] );
        if ( len < MIN_BOND_LEN ) {
            if ( i == 1 && (at[at_1].bUsed0DParity || at[at_2].bUsed0DParity) ) {
                pnt[i][0] = 0.0;
                pnt[i][1] = 1.0;
                pnt[i][2] = 0.0;
                len = 1.0; /* standard at_1-->at_2 vector coordinates in case of 0D allene */
            } else {
                goto exit_function; /*  too short bond */
            }
        }
        mult3( pnt[i], 1.0/len, pnt[i] );
    }
    len = 100.0*triple_prod(pnt[0], pnt[1], pnt[2], NULL );
/*
 *   ^ pnt[0]
 *   |                         The orientation on this diagram
 *   |                         produces len = -100
 *  [at_1]------>[at_2]
 *        pnt[1]    /
 *                 /
 *                / pnt[2]  (up from the plane)
 *               v
 *
 * Note: len is invariant upon at_1 <--> at_2 transposition because
 *       triple product changes sign upon pnt[0]<-->pnt[2] transposition and
 *       triple product changes sign upon pnt[1]--> -pnt[1] change of direction:
 *
 * triple_prod(pnt[0],  pnt[1], pnt[2], NULL ) =
 * triple_prod(pnt[2], -pnt[1], pnt[0], NULL )
 *
 */
    
    ret = len >= 0.0? (int)(floor(len+0.5)) : -(int)(floor(0.5-len));

exit_function:

    return ret;
}


/****************************************************************/

#if ( NEW_STEREOCENTER_CHECK == 1 ) /* { */

/********************************************************************************************/
int bInpAtomHasRequirdNeigh ( inp_ATOM *at, int cur_at, int RequirdNeighType, int NumDbleBonds )
{
    /* RequirdNeighType:
          reqired neighbor types (bitmap):
               0 => any neighbors
               1 => no terminal hydrogen atom neighbors
               2 => no terminal -X and -XH together (don't care about -X, -XH bond type, charge, radical)
                    (X = tautomeric endpoint atom)
       NumDbleBonds:
          if non-zero then allow double, alternating and tautomeric bonds
    */
    int i, j, ni, nj, bond_type, num_1s, num_mult, num_other;

    if ( at[cur_at].endpoint ) {  /*  tautomeric endpoint cannot be a stereo center */
        return 0;
    }

    if ( (1 & RequirdNeighType) && at[cur_at].num_H ) {
        return 0;
    }

    if ( 2 & RequirdNeighType ) {
        for ( i = 0; i < at[cur_at].valence; i ++ ) {
            ni = (int)at[cur_at].neighbor[i];
            if ( at[ni].valence != 1 ||
                 !get_endpoint_valence( at[ni].el_number ) ) {
                continue;
            }
            for ( j = i+1; j < at[cur_at].valence; j ++ ) {
                nj = (int)at[cur_at].neighbor[j];
                if ( at[nj].valence != 1 || 
                     at[ni].el_number != at[nj].el_number || 
                     !get_endpoint_valence( at[nj].el_number ) ) {
                    continue;
                }
                /*
                 * if (at[ni].num_H != at[nj].num_H) then the atoms (neighbors of at[cur_at]
                 * are tautomeric endpoints and are indistinguishable => cur_at is not stereogenic
                 * if (at[ni].num_H == at[nj].num_H) then the neighbors are indistinguishable
                 * and cur_at will be found non-sterogenic later
                 * get_endpoint_valence() check will not allow the neighbors to be carbons
                 * Therefore the following "if" is not needed; we may just return 0.
                 */
                if ( at[ni].num_H != at[nj].num_H && strcmp(at[ni].elname, "C" ) ) {
                    return 0; /*  found -X and -XH neighbors */
                }
            }
        }
    }
    
    num_1s = num_mult = num_other = 0;

    for ( i = 0; i < at[cur_at].valence; i ++ ) {
        bond_type = (at[cur_at].bond_type[i] & ~BOND_MARK_ALL);
        switch( bond_type ) {
        case BOND_SINGLE:
            num_1s ++;
            break;
        case BOND_DOUBLE:
        case BOND_ALTERN:
        case BOND_TAUTOM:
        case BOND_ALT12NS:
            num_mult ++;
            break;
        default:
            num_other ++;
            break;
        }
    }
    
    if ( num_other ) {
        return 0;
    }

    if ( NumDbleBonds && NumDbleBonds > num_mult ||
         !NumDbleBonds && at[cur_at].valence != num_1s ) {
        return 0;
    }
    return 1;
}
/********************************************************************************************/
int bCanInpAtomBeAStereoCenter( inp_ATOM *at, int cur_at, int bPointedEdgeStereo )
{
    
/*************************************************************************************
 * current version
 *************************************************************************************
 *       Use #define to split the stereocenter description table into parts
 *       to make it easier to read
 *
 *                      --------- 4 single bonds stereocenters -------
 *                       0       1       2       3      4       5
 *                                                                    
 *                       |       |       |       |      |       |     
 *                      -C-     -Si-    -Ge-    -Sn-   >As[+]  >B[-]  
 *                       |       |       |       |      |       |     
 */
#define SZELEM1         "C\000","Si",   "Ge",   "Sn",   "As", "B\000",
#define CCHARGE1         0,      0,      0,      0,      1,   -1,    
#define CNUMBONDSANDH1   4,      4,      4,      4,      4,    4,    
#define CCHEMVALENCEH1   4,      4,      4,      4,      4,    4,    
#define CHAS3MEMBRING1   0,      0,      0,      0,      0,    0,    
#define CREQUIRDNEIGH1   0,      0,      0,      0,      3,    0,    
/*
 *                      --------------- S, Se stereocenters ----------
 *                       6       7       8       9      10     11    12    13
 *                                                                                
 *                               |       |      ||             |     |     ||     
 *                      -S=     =S=     -S[+]   >S[+]   -Se=  =Se=  -Se[+] >Se[+] 
 *                       |       |       |       |       |     |     |      |     
 */
#define SZELEM2         "S\000","S\000","S\000","S\000","Se", "Se", "Se",  "Se",  
#define CCHARGE2         0,      0,      1,      1,      0,    0,    1,     1,    
#define CNUMBONDSANDH2   3,      4,      3,      4,      3,    4,    3,     4,    
#define CCHEMVALENCEH2   4,      6,      3,      5,      4,    6,    3,     5,    
#define CHAS3MEMBRING2   0,      0,      0,      0,      0,    0,    0,     0,    
#define CREQUIRDNEIGH2   3,      3,      3,      3,      3,    3,    3,     3,
/*
 *                      ------------------ N, P stereocenters -----------------
 *                        14     15       16     17     18       19       20    
 *                                                                               
 *                                                             Phosphine Arsine  
 *                                      X---Y                                    
 *                        |      |       \ /     |       |       \ /      \ /    
 *                       =N-    >N[+]     N     >P[+]   =P-       P        As    
 *                        |      |        |      |       |        |        |     
 */                                                                              
#define SZELEM3         "N\000","N\000","N\000","P\000","P\000","P\000", "As",
#define CCHARGE3         0,      1,      0,      1,      0,      0,       0,     
#define CNUMBONDSANDH3   4,      4,      3,      4,      4,      3,       3,     
#define CCHEMVALENCEH3   5,      4,      3,      4,      5,      3,       3,     
#define CHAS3MEMBRING3   0,      0,      1,      0,      0,      0,       0,     
#define CREQUIRDNEIGH3   3,      3,      1,      3,      3,      2,       2,     

#define PHOSPHINE_STEREO  19  /* the number must match Phosphine number in the comments, see above */
#define ARSINE_STEREO     20  /* the number must match Arsine number in the comments, see above */

    static char        szElem[][3]={ SZELEM1         SZELEM2         SZELEM3        };
    static S_CHAR        cCharge[]={ CCHARGE1        CCHARGE2        CCHARGE3       };
    static S_CHAR  cNumBondsAndH[]={ CNUMBONDSANDH1  CNUMBONDSANDH2  CNUMBONDSANDH3 };
    static S_CHAR  cChemValenceH[]={ CCHEMVALENCEH1  CCHEMVALENCEH2  CCHEMVALENCEH3 };
    static S_CHAR  cHas3MembRing[]={ CHAS3MEMBRING1  CHAS3MEMBRING2  CHAS3MEMBRING3 };
    static S_CHAR  cRequirdNeigh[]={ CREQUIRDNEIGH1  CREQUIRDNEIGH2  CREQUIRDNEIGH3 };

    static int n = sizeof(szElem)/sizeof(szElem[0]);
    /* reqired neighbor types (bitmap):
       0 => check bonds only
       1 => no terminal hydrogen atom neighbors
       2 => no terminal -X and -XH together (don't care the bond type, charge, radical)
            (X = tautomeric endpoint atom)
       Note: whenever cChemValenceH[] > cNumBondsAndH[]
             the tautomeric and/or alternating bonds
             are permitted

    */
    int i, ret = 0;
    for ( i = 0; i < n; i++ ) {
        if ( !strcmp( at[cur_at].elname, szElem[i]) &&
             at[cur_at].charge == cCharge[i] &&
             (!at[cur_at].radical || at[cur_at].radical == 1) &&
             at[cur_at].valence           +at[cur_at].num_H == cNumBondsAndH[i] &&
             at[cur_at].chem_bonds_valence+at[cur_at].num_H == cChemValenceH[i] &&
             (cHas3MembRing[i]? is_atom_in_3memb_ring( at, cur_at ) : 1) &&
             bInpAtomHasRequirdNeigh ( at, cur_at, cRequirdNeigh[i], cChemValenceH[i]-cNumBondsAndH[i]) ) {
            ret = cNumBondsAndH[i];
            break;
        }
    }

    if ( i == PHOSPHINE_STEREO && !(bPointedEdgeStereo & PES_BIT_PHOSPHINE_STEREO) )
        ret = 0;
    if ( i == ARSINE_STEREO && !(bPointedEdgeStereo & PES_BIT_ARSINE_STEREO) )
        ret = 0;
    return ret;
}

#else /* } NEW_STEREOCENTER_CHECK { */

/********************************************************************************************/
int bCanAtomBeAStereoCenter( char *elname, S_CHAR charge, S_CHAR radical )
{
    static const char   szElem[][3] = { "C\000", "Si", "Ge", "N\000", "P\000", "As", "B\000" };
    static const S_CHAR   cCharge[] = {  0,        0,    0,   1,       1,       1,    -1     };
    int i, ret = 0;
    for ( i = 0; i < sizeof(szElem)/sizeof(szElem[0]); i++ ) {
        if ( !strcmp( elname, szElem[i] )  && (charge == cCharge[i]) ) {
            ret = (!radical || radical == RADICAL_SINGLET);
            break;
        }
    }
    return ret;
}
#endif /* } NEW_STEREOCENTER_CHECK */

/****************************************************************/
/*  used for atoms adjacent to stereogenic bonds only */
int bAtomHasValence3( char *elname, S_CHAR charge, S_CHAR radical )
{
    static const char   szElem[][3] = {  "N\000" };
    static const S_CHAR   cCharge[] = {   0,     };
    int i, ret = 0;
    for ( i = 0; i < (int)(sizeof(szElem)/sizeof(szElem[0])); i++ ) {
        if ( !strcmp( elname, szElem[i] ) && (charge == cCharge[i]) ) {
            ret = ( !radical || radical == RADICAL_SINGLET );
            break;
        }
    }
    return ret;
}

/****************************************************************/
/*  used for atoms adjacent to stereogenic bonds only */
int bCanAtomHaveAStereoBond( char *elname, S_CHAR charge, S_CHAR radical )
{
    static const char   szElem[][3] = { "C\000", "Si", "Ge", "N\000", "N\000" };
    static const S_CHAR   cCharge[] = {  0,        0,    0,   0,       1,     };
    static const int       n = sizeof(szElem)/sizeof(szElem[0]);
    int i, ret = 0;
    for ( i = 0; i < n; i++ ) {
        if ( !strcmp( elname, szElem[i] )  && (charge == cCharge[i]) ) {
            ret = (!radical || radical == RADICAL_SINGLET);
            break;
        }
    }
    return ret;
}
/****************************************************************/
/*  used for atoms adjacent to stereogenic bonds only */
int bCanAtomBeMiddleAllene( char *elname, S_CHAR charge, S_CHAR radical )
{
    static const char   szElem[][3] = { "C\000", "Si", "Ge",  };
    static const S_CHAR   cCharge[] = {  0,        0,    0,   };
    static const int       n = sizeof(szElem)/sizeof(szElem[0]);
    int i, ret = 0;
    for ( i = 0; i < n; i++ ) {
        if ( !strcmp( elname, szElem[i] )  && (charge == cCharge[i]) ) {
            ret = (!radical || radical == RADICAL_SINGLET);
            break;
        }
    }
    return ret;
}
/*****************************************************************/
int bIsSuitableHeteroInpAtom( inp_ATOM  *at )
{
    int val, num_H;
    if ( 0 == at->charge &&
         (!at->radical || RADICAL_SINGLET == at->radical) &&
         0 < (val=get_endpoint_valence( at->el_number ) )) {
        num_H = at->num_H;
        if ( val == at->chem_bonds_valence + num_H ) {
            switch( val ) {
            case 2: /* O */
                if ( !num_H && 1 == at->valence )
                    return 0; /* =O */
                break;        /* not found */
            case 3: /* N */
                if ( 1 == at->valence && 1 == num_H ||
                     2 == at->valence && 0 == num_H  )
                    return 1; /* =N- or =NH */
                break;        /* not found */
            }
        }
    }
    return -1;
}
/****************************************************************/
int bIsOxide( inp_ATOM  *at, int cur_at )
{
    int i, bond_type;
    inp_ATOM  *a = at + cur_at, *an;
    for ( i = 0; i < a->valence; i ++ ) {
        bond_type = (a->bond_type[i] &= ~BOND_MARK_ALL);
        if ( bond_type == BOND_DOUBLE ) {
            an = at + (int)a->neighbor[i];
            if ( 1 == an->valence &&
                 !an->charge && !an->num_H && !an->radical &&
                 2 == get_endpoint_valence( an->el_number ) ) {
                return 1;
            }
        } else
        if ( bond_type == BOND_TAUTOM || bond_type == BOND_ALT12NS ) {
            an = at + (int)a->neighbor[i];
            if ( 1 == an->valence &&
                 2 == get_endpoint_valence( an->el_number ) ) {
                return 1;
            }
        }
    }
    return 0;
}
/****************************************************************/
/*  used for atoms adjacent to stereogenic bonds only */
int bCanAtomBeTerminalAllene( char *elname, S_CHAR charge, S_CHAR radical )
{
    static const char   szElem[][3] = { "C\000", "Si", "Ge",  };
    static const S_CHAR   cCharge[] = {  0,        0,    0,   };
    static const int       n = sizeof(szElem)/sizeof(szElem[0]);
    int i, ret = 0;
    for ( i = 0; i < n; i++ ) {
        if ( !strcmp( elname, szElem[i] ) && (charge == cCharge[i]) ) {
            ret = (!radical || radical == RADICAL_SINGLET);
            break;
        }
    }
    return ret;
}
/************************************************************************/
int GetHalfStereobond0DParity( inp_ATOM *at, int cur_at, AT_NUMB nSbNeighOrigAtNumb[],
                               int nNumExplictAttachments, int bond_parity, int nFlag )
{
    int m, last_parity, cur_parity;
    int i, icur2nxt, icur2neigh, cur_order_parity, nxt_at;
    AT_NUMB nNextSbAtOrigNumb;
    /* find atom parities for all valid streobonds incident to at[cur_at] */
    for ( m = 0, last_parity = 0; m < MAX_NUM_STEREO_BONDS && at[cur_at].sb_parity[m]; m ++ ) {
        icur2nxt = icur2neigh = -1; /* ordering number of neighbors in nSbNeighOrigAtNumb[] */
        cur_parity = 0;             /* parity for mth stereobond incident to the cur_at */
        if ( 0 <= at[cur_at].sb_ord[m] && at[cur_at].sb_ord[m] < at[cur_at].valence &&
             0 <= (nxt_at = at[cur_at].neighbor[(int)at[cur_at].sb_ord[m]]) &&
             at[nxt_at].valence <= MAX_NUM_STEREO_BONDS && /* make sure it is a valid stereobond */
             (nNextSbAtOrigNumb = at[nxt_at].orig_at_number) ) {
            /* since at[cur_at].sn_ord[m] = -1 for explicit H use at[cur_at].sn_orig_at_num[m] */
            for ( i = 0; i < nNumExplictAttachments; i ++ ) {
                if ( at[cur_at].sn_orig_at_num[m] == nSbNeighOrigAtNumb[i] ) {
                    icur2neigh = i; /* neighbor */
                } else
                if ( nNextSbAtOrigNumb == nSbNeighOrigAtNumb[i] ) {
                    icur2nxt = i; /* atom connected by a stereobond */
                }
            }
            if ( icur2neigh >= 0 && icur2nxt >= 0 ) {
                if ( ATOM_PARITY_WELL_DEF(at[cur_at].sb_parity[m]) ) {
                    /* parity of at[cur_atom] neighbor permutation to reach this order: { next_atom, neigh_atom, ...} */
                    cur_order_parity = (icur2nxt + icur2neigh + (icur2nxt > icur2neigh) - 1) % 2;
                    cur_parity = 2 - (cur_order_parity + at[cur_at].sb_parity[m]) % 2;
                } else {
                    /* unknowm/undef parities do not depend on the neighbor order */
                    cur_parity = at[cur_at].sb_parity[m];
                }
            }
        } else {
            continue;
        }
        /* use a well-known parity if available; if not then use preferably the unknown */
        if ( !last_parity ) {
            last_parity = cur_parity;
        } else
        if ( last_parity != cur_parity && cur_parity ) {
            if ( ATOM_PARITY_WELL_DEF(last_parity) ) {
                if ( ATOM_PARITY_WELL_DEF(cur_parity) ) {
                    last_parity = 0; /* error: all well-defined parities should be same */
                    break;
                }
            } else
            if ( ATOM_PARITY_WELL_DEF(cur_parity) ) {
                /* replace unknown/undefined parity with well-known */
                last_parity = cur_parity;
            } else {
                /* select min unknown/undefined parity (out of AB_PARITY_UNKN and AB_PARITY_UNDF) */
                last_parity = inchi_min(cur_parity, last_parity);
            }
        }
    }
    if ( last_parity ) {
        bond_parity = last_parity;
        at[cur_at].bUsed0DParity |= nFlag; /* set flag: used stereobond 0D parity */
    }
    return bond_parity;
}
/*******************************************************************************************/
int FixSb0DParities( inp_ATOM *at, /* inp_ATOM *at_removed_H, int num_removed_H,*/ int chain_length,
                     int at_1, int i_next_at_1, S_CHAR z_dir1[],
                     int at_2, int i_next_at_2, S_CHAR z_dir2[],
                     int *pparity1, int *pparity2 )
{
    int k, parity1, parity2, abs_parity1, abs_parity2;
    int j1, j2, parity_sign;
    /*
    AT_NUMB nSbNeighOrigAtNumb1[MAX_NUM_STEREO_BOND_NEIGH], nSbNeighOrigAtNumb2[MAX_NUM_STEREO_BOND_NEIGH];
    int     nNumExplictAttachments1, nNumExplictAttachments2;
    */
    parity1 = parity2 = AB_PARITY_NONE;
    j1      = j2      = -1;
    parity_sign = ( *pparity1 < 0 || *pparity2 < 0 )? -1 : 1;

    abs_parity1 = abs(*pparity1);
    abs_parity2 = abs(*pparity2);

    for ( k = 0; k < MAX_NUM_STEREO_BONDS && at[at_1].sb_parity[k]; k ++ ) {
        if ( at[at_1].sb_ord[k] == i_next_at_1 ) {
            parity1 = at[at_1].sb_parity[k];
            j1 = k;
        }
    }
    for ( k = 0; k < MAX_NUM_STEREO_BONDS && at[at_2].sb_parity[k]; k ++ ) {
        if ( at[at_2].sb_ord[k] == i_next_at_2 ) {
            parity2 = at[at_2].sb_parity[k];
            j2 = k;
        }
    }
    switch( (j1 >= 0) + 2*(j2 >= 0) ) {
    case 0:
        /* the bond has no 0D parity */
        *pparity1 = *pparity2 = parity_sign * AB_PARITY_UNDF;
        return 0;
    case 1:
    case 2:
        /* 0D parity data error */
        *pparity1 = *pparity2 =  AB_PARITY_NONE;
        return -1;
    case 3:
        /* the bond has 0D parity */
        switch (     !(ATOM_PARITY_WELL_DEF( abs_parity1 ) && ATOM_PARITY_WELL_DEF( parity1 )) +
                 2 * !(ATOM_PARITY_WELL_DEF( abs_parity2 ) && ATOM_PARITY_WELL_DEF( parity2 )) ) {
        case 0:
            /* both parities are well-defined; continue */
            break;
        case 1:
            /* 0D parity not well-defined for at_1 */
            *pparity1 = parity_sign * (ATOM_PARITY_WELL_DEF( parity1     )? abs_parity1 :
                                       ATOM_PARITY_WELL_DEF( abs_parity1 )? parity1 :
                                       inchi_min(abs_parity1, parity1));
            *pparity2 = parity_sign * abs_parity2;
            return -1;
        case 2:
            /* 0D parity not well-defined for at_2 */
            *pparity1 = parity_sign * abs_parity1;
            *pparity2 = parity_sign * (ATOM_PARITY_WELL_DEF( parity2     )? abs_parity2 :
                                       ATOM_PARITY_WELL_DEF( abs_parity2 )? parity2 :
                                       inchi_min(abs_parity2, parity2));
            return -1;
        case 3:
            abs_parity1 =  (ATOM_PARITY_WELL_DEF( parity1     )? abs_parity1 :
                            ATOM_PARITY_WELL_DEF( abs_parity1 )? parity1 :
                            inchi_min(abs_parity1, parity1));
            abs_parity2 =  (ATOM_PARITY_WELL_DEF( parity2     )? abs_parity2 :
                            ATOM_PARITY_WELL_DEF( abs_parity2 )? parity2 :
                            inchi_min(abs_parity2, parity2));
            *pparity1 = *pparity2 = parity_sign * inchi_min(abs_parity1, abs_parity2);
            /*return (parity1 == parity2)? 0 : -1;*/
            return -1;
        }
        break;
    }
    /* we are here if both end-atoms of the bond have well-defined 0D parities */
    /*
    nNumExplictAttachments1 = GetSbNeighOrigAtNumb( at, at_1, at_removed_H, num_removed_H, nSbNeighOrigAtNumb1 );
    nNumExplictAttachments2 = GetSbNeighOrigAtNumb( at, at_2, at_removed_H, num_removed_H, nSbNeighOrigAtNumb2 );
    parity1 = GetHalfStereobond0DParity( at, at_1, nSbNeighOrigAtNumb1, nNumExplictAttachments1, *pparity1, 0 );
    parity2 = GetHalfStereobond0DParity( at, at_2, nSbNeighOrigAtNumb2, nNumExplictAttachments2, *pparity2, 0 );
    */
    *pparity1 = parity_sign * abs_parity1;
    *pparity2 = parity_sign * abs_parity2;

    if ( chain_length % 2 ) {
        /* allene; chain_length = (number of double bonds) - 1 */
        /*
        int zer1 = ( !z_dir1[0] && !z_dir1[1] && !z_dir1[2] );
        int zer2 = ( !z_dir2[0] && !z_dir2[1] && !z_dir2[2] );
        */
        int bWrong_z_dir1 = (0 != (at[at_1].bUsed0DParity & FlagSB_0D));
        int bWrong_z_dir2 = (0 != (at[at_2].bUsed0DParity & FlagSB_0D));

        if ( bWrong_z_dir1 && bWrong_z_dir2 ) {
            goto set_default;
        } else
        if ( bWrong_z_dir1 || bWrong_z_dir2 ) {
            double r12[3], zi1[3], zi2[3], abs_r12, abs_zi2;
            int    at_i1, at_i2, j;
            S_CHAR   z_dir[3];
            r12[0] = at[at_2].x - at[at_1].x;
            r12[1] = at[at_2].y - at[at_1].y;
            r12[2] = at[at_2].z - at[at_1].z;
            abs_r12 = len3( r12 );
            if ( abs_r12 < MIN_BOND_LEN ) {
                goto set_default;
            }
            /* make r12[] point to the atom with 'good' z_dir[] */
            if ( bWrong_z_dir1 ) {
                at_i1 = at_2; /* has good z_dir2[] */
                at_i2 = at_1; /* has bad  z_dir1[] */
                zi1[0] = z_dir2[0];
                zi1[1] = z_dir2[1];
                zi1[2] = z_dir2[2];
                mult3( r12, 1.0/abs_r12, r12 ); /* make length = 1 */
            } else {
                at_i1 = at_1; /* has good z_dir1[] */
                at_i2 = at_2; /* has bad  z_dir2[] */
                zi1[0] = z_dir1[0];
                zi1[1] = z_dir1[1];
                zi1[2] = z_dir1[2];
                mult3( r12, -1.0/abs_r12, r12 ); /* make length = 1 */
            }
            cross_prod3( r12, zi1, zi2 );
            abs_zi2 = len3( zi2 );
            mult3( zi2, 100.0/abs_zi2, zi2 ); /* make length = 100 */
            for ( j = 0; j < 3; j ++ ) {
                z_dir[j] = (S_CHAR) (zi2[j]>= 0.0?  floor(0.5 + zi2[j]) :
                                                   -floor(0.5 - zi2[j])); /*  abs(z_dir) = 100 */
            }
            if ( bWrong_z_dir1 ) {
                memcpy( z_dir1, z_dir, sizeof(z_dir) );
            } else {
                memcpy( z_dir2, z_dir, sizeof(z_dir) );
            }
        }
        return 0;

set_default:
        /* z_dir1[] = x-direction; z_dir2[] = z-direction; r12[] = y-direction */
        z_dir1[0] = 100;
        z_dir1[1] = z_dir1[2] = 0;
        z_dir2[0] = z_dir2[1] = 0;
        z_dir2[2] = 100;
    }
    return 0;
}
/**********************************************************/
/* without this InChI fails on reconstructed  CID=450438  */
/* (isotopic, Unknown SB adjacent to SB with known parity) */
/**********************************************************/
int FixUnkn0DStereoBonds(inp_ATOM *at, int num_at)
{
    int i, m, num=0;

    /* add usual Unknown stereobond descriptors to each Unknown bond */
    for( i = 0; i < num_at; i ++ ) {
        for ( m = 0; m < MAX_NUM_STEREO_BONDS && at[i].sb_parity[m]; m ++ ) {
            if ( AB_PARITY_UNKN == at[i].sb_parity[m] ) {
                at[i].bond_stereo[ (int)at[i].sb_ord[m] ] = STEREO_DBLE_EITHER;
                num ++;
            }
        }
    }
#ifdef NEVER
    if ( num ) {
        int j;
        /* how to remove Unknown stereo bond parities */
        for( i = 0; i < num_at; i ++ ) {
            for ( m = 0; m < MAX_NUM_STEREO_BONDS && at[i].sb_parity[m]; m ++ ) {
                if ( AB_PARITY_UNKN == at[i].sb_parity[m] ) {
                    for ( j = m+1; j < MAX_NUM_STEREO_BONDS; j ++ ) {
                        at[i].sb_parity[j-1]      = at[i].sb_parity[j];
                        at[i].sb_ord[j-1]         = at[i].sb_ord[j];
                        at[i].sn_ord[j-1]         = at[i].sn_ord[j];
                        at[i].sn_orig_at_num[j-1] = at[i].sn_orig_at_num[j];
                    }
                    at[i].sb_parity[j-1]      = 0;
                    at[i].sb_ord[j-1]         = 0;
                    at[i].sn_ord[j-1]         = 0;
                    at[i].sn_orig_at_num[j-1] = 0;
                }
            }
        }
    }
#endif
    return num;
}
/*======================================================================================================
 
half_stereo_bond_parity() General Description:

    A) find projections of 3 bonds on a reasonable plane defined
       by a vector z_dir perpendicular to the plane
    B) calculate parity

half_stereo_bond_parity() Detailed Description:

    1) Find at_coord[] = vectors from the central atoms to its neighbors
    2) If only 2 neighbors are present, then create a reasonable 3rd neighbor
       (an implicit H or a fictitious atom in case of =NX) coordinates
    3) Normalize at_coord[] to unit length
    4) Find unit vector pnt[2] perpendicular to the plane containing
       at_coord[] arrow ends.
       Even though it is not necessary, make z-coordinate of pnt[2] positive.
       ** pnt[2] has the new z-axis direction **
    5) Let pnt[0] = perpendicular to pnt[2] component of at_coord[0];
       Normalize pnt[0] to unit length.
       ** pnt[0] has the new x-axis direction **
    6) Let pnt[1] = pnt[2] x pnt[0] (cross-product);
       ** pnt[1] has the new y-axis direction **
    7) Find at_coord[] in the new xyz-basis and normalize their xy-projections
       to a unit length
    8) In the new xy-plane find (counterclockwise) angles:
       tmp1 = (from at_coord[0] to at_coord[1])
       tmp2 = (from at_coord[0] to at_coord[2])
    9) Calculate the parity: if tmp1 < tmp2 then 1 (odd) else 2 (even)
       (even: looking from the arrow end of the new z-axis, 0, 1, and 2 neighbors
        are in clockwise order)
   10) Calculate z_dir = 100*pnt[2].
   
   Note1. If z_dir vectors of atoms located at the opposite ends of a double bond have approximately
          opposite directions (that is, their dot-product is negative) then the parity of the
          stereogenic bond calculated from half-bond-parities should be inverted

   Note2. In case of a tetrahedral cumulene a triple product (z_dir1, (1->2), z_dir2) is used instead
          of the dot-product. (1->2) is a vector from the atom#1 to the atom #2. This triple product
          is invariant with respect to the atom numbering because it does not change upon (1,2)
          permutation.
  
  Stereo ambiguity in case of 2 neighbors:
  ----------------------------------------
  Undefined: single-double bond angle > pi - arcsin(0.03) = 178.28164199834454285275613218975 degrees
  Ambiguous: single-double bond angle > 175 degrees = pi - 0.087156 Rad
 
   Return values 
   (cases: I=only in case of isotopic H atoms the neighbors are different,
           N=in case of non-isotopic H atoms the neighbors are different)

  -4 = AB_PARITY_UNDF => atom is adjacent to a stereogenic bond, but the geometry is undefined, I
  -3 = AB_PARITY_UNKN => atom is adjacent to a stereogenic bond, but the geometry is not known to the iuser, I
  -2 =-AB_PARITY_EVEN => parity of an atom adjacent to a stereogenic bond, I
  -1 =-AB_PARITY_ODD  => parity of an atom adjacent to a stereogenic bond, I
   0 = AB_PARITY_NONE => the atom is not adjacent to a stereogenic bond
   1 = AB_PARITY_ODD  => parity of an atom adjacent to a stereogenic bond, N&I
   2 = AB_PARITY_EVEN => parity of an atom adjacent to a stereogenic bond, N&I
   3 = AB_PARITY_UNKN => atom is adjacent to a stereogenic bond, but the geometry is not known to the iuser, N&I
   4 = AB_PARITY_UNDF => atom is adjacent to a stereogenic bond, but the geometry is undefined, N&I
   5 = AB_PARITY_IISO => atom constitutionally equivalent to this atom may be adjacent to a stereogenic bond, I


=====================================================================================================*/

int half_stereo_bond_parity( inp_ATOM *at, int cur_at, inp_ATOM *at_removed_H, 
                            int num_removed_H, S_CHAR *z_dir, 
                            int bPointedEdgeStereo, int vABParityUnknown )
{
    double at_coord[MAX_NUM_STEREO_BOND_NEIGH][3], c, s, tmp[3], tmp1, tmp2, min_tmp, max_tmp, z;
    double temp[3], pnt[3][3];
    int j, k, p0, p1, p2, next, bValence3=0, num_z, nType, num_either_single, num_either_double;
    int nNumExplictAttachments;
    int bond_parity  =  AB_PARITY_UNDF;
    int    num_H=0, num_iH, num_eH=0, num_nH=0 /* = num_iso_H[0] */;
    int    num_iso_H[NUM_H_ISOTOPES+1];
    int    index_H[5]; /*  cannot have more than 4 elements: 1 H, 1 1H, 1 D, 1 T atom(s) */
    /*	const double one_pi = 2.0*atan2(1.0 , 0.0 ); */
    const double one_pi = 3.14159265358979323846; /* M_PI */
    const double two_pi = 2.0*one_pi;
    int    bIgnoreIsotopicH = (0 != (at[cur_at].cFlags & AT_FLAG_ISO_H_POINT));
    AT_NUMB nSbNeighOrigAtNumb[MAX_NUM_STEREO_BOND_NEIGH];


    if ( z_dir && !z_dir[0] && !z_dir[1] && !z_dir[2] ) {
        z_dir[2]=100;
    }

    num_H  = at[cur_at].num_H;
    if ( num_H > NUM_H_ISOTOPES )
        return 0; /*  at least 2 H atoms are isotopically identical */

    if ( MAX_NUM_STEREO_BOND_NEIGH < at[cur_at].valence + num_H ||
         MIN_NUM_STEREO_BOND_NEIGH > at[cur_at].valence + num_H )
        return 0;

    if ( !bCanAtomHaveAStereoBond( at[cur_at].elname, at[cur_at].charge, at[cur_at].radical ) )
        return 0;
    if ( !bIgnoreIsotopicH ) {
        for ( j = 0, num_nH = num_H; j < NUM_H_ISOTOPES; j ++ ) {
            if ( (k = (int)at[cur_at].num_iso_H[j]) > 1 ) {
                return AB_PARITY_IISO;  /*  two or more identical isotopic H atoms */
            }
            num_nH -= k;
        }
    }
    /*  at this point num_nH = number of non-isotopic H atoms */
    if ( num_nH > 1 )
        return AB_PARITY_IISO; /*  two or more identical non-isotopic H atoms */
    if ( num_nH < 0 )
        return CT_ISO_H_ERR;  /*  program error */ /*   <BRKPT> */

    /********************************************************************
     * Note. At this point all (implicit and explicit) isotopic
     * terminal H neighbors are either different or not present.
     ********************************************************************/

    /*  locate explicit hydrogen atoms */
    /*  (at_removed_H are sorted in ascending isotopic H mass order, non-isotopic first) */
    memset( num_iso_H, 0, sizeof(num_iso_H) );
    if ( at_removed_H && num_removed_H > 0 ) {
        for ( j = 0; j < num_removed_H; j ++ ) {
            if ( at_removed_H[j].neighbor[0] == cur_at ) {
                k = bIgnoreIsotopicH? 0 : at_removed_H[j].iso_atw_diff;
                if ( 0 <= k && k <= NUM_H_ISOTOPES ) {
                    if ( ++num_iso_H[k] > 1 )   /*  num_iso_H[0] = number of non-isotopic H atoms */
                        return CT_ISO_H_ERR;    /*  program error in counting hydrogens */ /*   <BRKPT> */
                    index_H[num_eH++] = j;
                } else {
                    return CT_ISO_H_ERR;  /*  program error */ /*   <BRKPT> */
                }
            }
        }
        num_iH = num_H - num_eH; /*  number of implicit non-isotopic and isotopic H atoms */
        if ( num_iH > 1 ) {
            /*  more than one implicit H: cannot reconstruct the geometry */
            bond_parity = -AB_PARITY_UNDF;
            goto exit_function;
        }
    } else {
        num_iH = num_H;
    }
    /*  at this point num_iH = number of implicit non-isotopic and isotopic H atoms */
    if ( at[cur_at].valence + num_eH < MIN_NUM_STEREO_BOND_NEIGH ) {
        /*  =NH or =CHD when no explicit H is present */
        return num_H == 1? AB_PARITY_UNDF : -AB_PARITY_UNDF;
    }

    bValence3 = bAtomHasValence3( at[cur_at].elname, at[cur_at].charge, at[cur_at].radical );
    /*
     * Can one explicit hydrogen be added to make asymmetric configuration?
     * For now we can add 1 H atom in case of an appropriate geometry if:
     * (a) one non-isotopic H (even if explicit isotopic H atoms are present), or
     * (b) one isotopic or non-isotopic H if NO explicit isotopic or non-isotopic H atom is present
     * This makes sense only in case chem. valence = 4. In case of chem. valence = 3, do not check.
     */
    if ( at[cur_at].valence + num_eH == MIN_NUM_STEREO_BOND_NEIGH && !bValence3 &&
         !(/*(a)*/ 1 == num_nH && !num_iso_H[0] ||
           /*(b)*/ 1 == num_H  && !num_eH) 
       ) {
        goto exit_function;
        /* return num_H == 1? AB_PARITY_UNDF : -AB_PARITY_UNDF; */
    }

    /*  store neighbors coordinates */
    num_z = num_either_single = num_either_double = 0;
    for ( k = nNumExplictAttachments = 0; k < 2; k ++ ) {
        switch( k ) {
        case 0:
            for ( j = 0; j < num_eH; j ++, nNumExplictAttachments ++ ) {
                next = index_H[j];
                at_coord[nNumExplictAttachments][0] = at_removed_H[next].x - at[cur_at].x;
                at_coord[nNumExplictAttachments][1] = at_removed_H[next].y - at[cur_at].y;
                nSbNeighOrigAtNumb[nNumExplictAttachments] = at_removed_H[next].orig_at_number;
                /* use the fact that (at_removed_H - at) = (number of atoms except removed explicit H) */
                z = -get_z_coord( at, (at_removed_H-at)+next, 0 /*neighbor #*/, &nType, -(bPointedEdgeStereo & PES_BIT_POINT_EDGE_STEREO) );
                switch ( nType ) {
                case ZTYPE_EITHER:
                    num_either_single ++; /*  bond in "Either" direction. */
                    break;
                case ZTYPE_UP:
                case ZTYPE_DOWN:
                    nType = -nType; /*  at_removed_H[] contains bonds TO the center, not from */
                    z = len2( at_coord[nNumExplictAttachments] );
                    /*
                    z = sqrt( at_coord[nNumExplictAttachments][0]*at_coord[nNumExplictAttachments][0]
                            + at_coord[nNumExplictAttachments][1]*at_coord[nNumExplictAttachments][1] );
                    */
                    if ( nType == ZTYPE_DOWN )
                        z = -z;
                    /*  no break; here */
                case ZTYPE_3D:
                    num_z ++;
                }
                at_coord[nNumExplictAttachments][2] = z;
            }
            break;
        case 1:
            for ( j = 0; j < at[cur_at].valence; j ++, nNumExplictAttachments ++ ) {
                next = at[cur_at].neighbor[j];
                at_coord[nNumExplictAttachments][0] = at[next].x - at[cur_at].x;
                at_coord[nNumExplictAttachments][1] = at[next].y - at[cur_at].y;
                nSbNeighOrigAtNumb[nNumExplictAttachments] = at[next].orig_at_number;

                z = get_z_coord( at, cur_at, j /*neighbor #*/, &nType, (bPointedEdgeStereo & PES_BIT_POINT_EDGE_STEREO) );
                switch ( nType ) {
                case ZTYPE_EITHER:
                    num_either_single ++; /*  bond in "Either" direction. */
                    break;
                case ZTYPE_UP:
                case ZTYPE_DOWN:
                    z = len2( at_coord[nNumExplictAttachments] );
                    /*
                    z = sqrt( at_coord[nNumExplictAttachments][0]*at_coord[nNumExplictAttachments][0]
                            + at_coord[nNumExplictAttachments][1]*at_coord[nNumExplictAttachments][1] );
                    */
                    if ( nType == ZTYPE_DOWN )
                        z = -z;
                    /*  no break; here */
                case ZTYPE_3D:
                    num_z ++;
                }
                at_coord[nNumExplictAttachments][2] = z;
            }
            break;
        }
    }

    if ( num_either_single ) {
        bond_parity =  vABParityUnknown /*AB_PARITY_UNKN*/;  /*  single bond is 'unknown' */
        goto exit_function;
    }

    /* nNumExplictAttachments is a total number of attachments, including removed explicit terminal hydrogens */
    if ( nNumExplictAttachments == 2 ) {
        /*  create coordinates of the implicit hydrogen (or a fictitious atom in case of ==N-X ), */
        /*  coord[2][], attached to the cur_at. */
        for ( j = 0; j < 3; j ++ ) {
            at_coord[2][j] = - ( at_coord[0][j] + at_coord[1][j] );
        }
        nSbNeighOrigAtNumb[nNumExplictAttachments] = 0; /* implicit H or lone pair */
    }
    for ( j = 0; j < 3; j ++ ) {
        tmp[j] = len3( at_coord[j] );
    }
    min_tmp = inchi_min( tmp[0], inchi_min(tmp[1], tmp[2]) );
    max_tmp = inchi_max( tmp[0], inchi_max(tmp[1], tmp[2]) );
    if ( min_tmp < MIN_BOND_LEN || min_tmp < MIN_SINE*max_tmp ) {
        /*  all bonds or some of bonds are too short */
        if ( at[cur_at].sb_parity[0] ) {
            /* use bond psrity; the reconciliation in ReconcileAllCmlBondParities()
             * has made all ways to calculate parity produce same result
             */
            bond_parity = GetHalfStereobond0DParity( at, cur_at, nSbNeighOrigAtNumb,
                                                     nNumExplictAttachments, bond_parity, FlagSB_0D );
        }
        
        goto exit_function; 
    }
    /*  normalize lengths to 1 */
    for ( j = 0; j < 3; j ++ ) {
        mult3( at_coord[j], 1.0/tmp[j], at_coord[j] );
    }

    /*  find projections of at_coord vector differences on the plane containing their arrowhead ends */
    for ( j = 0; j < 3; j ++ ) {
        /*  pnt[0..2] = {0-1, 1-2, 2-0} */
        tmp[j] = len3(diff3( at_coord[j], at_coord[(j+1)%3], pnt[j] ));
        if ( tmp[j] < MIN_SINE ) {
            goto exit_function; /*  angle #i-cur_at-#j is too small */
        }
        mult3( pnt[j], 1.0/tmp[j], pnt[j] ); /* 2003-10-06 */
    }
    /*  find pnt[p2], a vector perpendicular to the plane, and its length tmp[p2] */
    /*  replace previous pnt[p2], tmp[p2] with new values; the old values do not have any additional */
    /*  information because pnt[p0]+pnt[p1]+pnt[p2]=0 */
    /*  10-6-2003: a cross-product of one pair pnt[j], pnt[(j+1)%3] can be very small. Find the larges one */
    tmp1 = len3( cross_prod3( pnt[0], pnt[1], temp ) );
    for (j = 1, k = 0; j < 3; j ++ ) {
        tmp2 = len3( cross_prod3( pnt[j], pnt[(j+1)%3], temp ) );
        if ( tmp2 > tmp1 ) {
            tmp1 = tmp2;
            k     = j;
        }
    }
    /* previously p0=0, p1=1, p2=2 */
    p0 = k;
    p1 = (k+1)%3;
    p2 = (k+2)%3;
    tmp[p2] = len3( cross_prod3( pnt[p0], pnt[p1], pnt[p2] ) );
    if ( tmp[p2] < MIN_SINE*tmp[p0]*tmp[p1]  ) {
        goto exit_function; /*  pnt[p0] is almost colinear to pnt[p1] */
    }
    /*  new basis: pnt[p0], pnt[p1], pnt[p2]; set z-coord sign and make abs(pnt[p2]) = 1 */
    mult3( pnt[p2], (pnt[p2][2]>0.0? 1.0:-1.0)/tmp[p2], pnt[p2] ); /*  unit vector in the new z-axis direction */

    min_tmp = dot_prod3( at_coord[0], pnt[p2] ); /*  non-planarity measure (sine): hight of at_coord[] pyramid */
    mult3( pnt[p2], min_tmp, pnt[p0] ); /*  vector height of the pyramid, ideally 0 */
    /*  find new pnt[p0] = projection of at_coord[p0] on plane orthogonal to pnt[p2] */
    tmp[p0] = len3(diff3( at_coord[0], pnt[p0], pnt[p0] ));
    mult3( pnt[p0], 1.0/tmp[p0], pnt[p0] );  /*  new x axis basis vector */
    cross_prod3( pnt[p2], pnt[p0], pnt[p1] ); /*  new y axis basis vector */
    /*  find at_coord in the new basis of {pnt[p0], pnt[p1], pnt[p2]} */
    for ( j = 0; j < 3; j ++ ) {
        copy3( at_coord[j], temp );
        for ( k = 0; k < 3; k ++ ) {
            at_coord[j][k] = dot_prod3( temp, pnt[(k+p0)%3] );
        }
        /*  new xy plane projection length */
        tmp[j] = sqrt(at_coord[j][0]*at_coord[j][0] + at_coord[j][1]*at_coord[j][1]);
        /*  make new xy plane projection length = 1 */
        mult3( at_coord[j], 1.0/tmp[j], at_coord[j] );
    }
   
    s = fabs( at_coord[1][0]*at_coord[2][1] - at_coord[1][1]*at_coord[2][0] ); /*  1-2 sine */
    c =       at_coord[1][0]*at_coord[2][0] + at_coord[1][1]*at_coord[2][1];   /*  1-2 cosine */
    if ( s < MIN_SINE && c > 0.5 ) {
        goto exit_function; /*  bonds to neigh. 1 and 2 have almost same direction; relative angles are undefined */
    }
    c = at_coord[0][0]; /*  cosine of the angle between new Ox axis and a bond to the neighbor 0. Should be 1 */
    s = at_coord[0][1]; /*  sine. Should be 0 */
    /*  turn vectors so that vector #1 (at_coord[0]) becomes {1, 0} */
    for ( j = 0; j < MAX_NUM_STEREO_BOND_NEIGH; j ++ ) {
        tmp1 =  c*at_coord[j][0] + s*at_coord[j][1];
        tmp2 = -s*at_coord[j][0] + c*at_coord[j][1];
        at_coord[j][0] = tmp1;
        at_coord[j][1] = tmp2;
    }
    /*  counterclockwise angles from the direction to neigh 0 to to directions to neighbors 1 and 2: */
    tmp1 = atan2( at_coord[1][1], at_coord[1][0] ); /*  range -pi and +pi */
    tmp2 = atan2( at_coord[2][1], at_coord[2][0] );
    if ( tmp1 < 0.0 )
        tmp1 += two_pi; /*  range 0 to 2*pi */
    if ( tmp2 < 0.0 )
        tmp2 += two_pi;
    /*-----------------------------------
                        Example
      1 \               case tmp1 < tmp2
         \              parity is odd
          \             (counterclockwise)
           A------- 0
          /
         /
      2 /

    ------------------------------------*/
    bond_parity = 2 - ( tmp1 < tmp2 );
    for ( j = 0; j < 3; j ++ ) {
        z_dir[j] = (S_CHAR) (pnt[p2][j]>= 0.0?  floor(0.5 + 100.0 * pnt[p2][j]) :
                                               -floor(0.5 - 100.0 * pnt[p2][j])); /*  abs(z_dir) = 100 */
    }
    /*  check for ambiguity */
    if ( nNumExplictAttachments > 2 ) {
        min_tmp = inchi_min( tmp1, tmp2 );
        max_tmp = inchi_max( tmp1, tmp2 );
        if ( min_tmp > one_pi-MIN_SINE || max_tmp < one_pi+MIN_SINE || max_tmp-min_tmp > one_pi - MIN_SINE ) {
            at[cur_at].bAmbiguousStereo |= AMBIGUOUS_STEREO;
        } else /* 3D ambiguity 8-28-2002 */
        if ( fabs(at_coord[0][2]) > MAX_SINE ) { /*  all fabs(at_coord[j][2] (j=0..2) must be equal */
            at[cur_at].bAmbiguousStereo |= AMBIGUOUS_STEREO;
        }
    } else
    if ( nNumExplictAttachments == 2 ) {  /* 10-6-2003: added */
        min_tmp = fabs(tmp1 - one_pi);
        if ( min_tmp < MIN_SINE ) {
            bond_parity = AB_PARITY_UNDF; /* consider as undefined 10-6-2003 */
        } else
        if ( min_tmp < MIN_ANGLE_DBOND ) {
            at[cur_at].bAmbiguousStereo |= AMBIGUOUS_STEREO;
        }
    }


    /*  for 3 neighbors moving implicit H to the index=0 from index=2 position */
    /*  can be done in 2 transpositions and does not change atom's parity */
exit_function:
    if ( num_H > 1 && bond_parity > 0 && !(bond_parity & AB_PARITY_0D) /*&& PARITY_WELL_DEF(bond_parity)*/ ) {
        /*
         * stereo only if isotopes are counted.             Do not inverse
         * Examples:                                        sign for this:
         *     H                            D               
         *    /                            /                    H
         * ==C                      or  ==CH                   / 
         *    \                                             ==N  (bValence3=1)
         *     D                  
         * two explicit         one explicit H isotope (D),
         * isotopic H atoms     one implicit H
         */
        bond_parity = -bond_parity; /*  refers to isotopically substituted structure only */
    }
    return bond_parity;
}

/*************************************************************/
int save_a_stereo_bond( int z_prod, int result_action,
                        int at1, int ord1, AT_NUMB *stereo_bond_neighbor1, S_CHAR *stereo_bond_ord1, S_CHAR *stereo_bond_z_prod1, S_CHAR *stereo_bond_parity1, 
                        int at2, int ord2, AT_NUMB *stereo_bond_neighbor2, S_CHAR *stereo_bond_ord2, S_CHAR *stereo_bond_z_prod2, S_CHAR *stereo_bond_parity2 )
{
    int i1, i2;
    for ( i1 = 0; i1 < MAX_NUM_STEREO_BONDS && stereo_bond_neighbor1[i1]; i1 ++ )
        ;
    for ( i2 = 0; i2 < MAX_NUM_STEREO_BONDS && stereo_bond_neighbor2[i2]; i2 ++ )
        ;
    if ( i1 == MAX_NUM_STEREO_BONDS || i2 == MAX_NUM_STEREO_BONDS )
        return 0;
    
    stereo_bond_parity1[i1] =
    stereo_bond_parity2[i2] = result_action;

    stereo_bond_neighbor1[i1] = (AT_NUMB) (at2+1);
    stereo_bond_ord1[i1]      = (S_CHAR)ord1;
    stereo_bond_neighbor2[i2] = (AT_NUMB) (at1+1);
    stereo_bond_ord2[i2]      = (S_CHAR)ord2;
    stereo_bond_z_prod1[i1]   =
    stereo_bond_z_prod2[i2]   = (S_CHAR)z_prod;
    return 1;
}
/***************************************************************/
int get_allowed_stereo_bond_type( int bond_type )
{
#if (ALLOW_TAUT_ATTACHMENTS_TO_STEREO_BONDS == 0 )
    if ( (bond_type & ~BOND_MARK_ALL) == BOND_TAUTOM )
        return 0; /*  no tautomer bonds allowed */
    else
#endif
#if ( EXCL_ALL_AROM_BOND_PARITY  == 1 )  /* { */
    /*  a stereo bond cannot belong to an aromatic atom */
    if ( (bond_type &= ~BOND_MARK_ALL) == BOND_ALTERN )
    {
        return 0;
    }
#else  /* } { */
#if ( ADD_6MEMB_AROM_BOND_PARITY == 1 )
    /*  accept any aromatic bond as a stereo bond */
    if ( (bond_type &= ~BOND_MARK_ALL) == BOND_ALTERN )
#else
    /*  accept only aromatic bonds in non-6-member rings */
    if ( (bond_type &= ~BOND_MARK_ALL) == BOND_ALTERN ) )
#endif
    {
        return BOND_ALTERN;
    }
#endif  /* } */
    else
    /*  at this point BOND_MARK_ALL bits have been removed from bond_type */
    if ( bond_type == BOND_DOUBLE || bond_type == BOND_SINGLE ) {
        return bond_type;
    }
#if (ALLOW_TAUT_ATTACHMENTS_TO_STEREO_BONDS == 1 )
    else
    if ( bond_type == BOND_TAUTOM ) {
        return BOND_TAUTOM;
    }
#endif

    return 0; /*  wrong bond type */
}

/*************************************************************/
int can_be_a_stereo_bond_with_isotopic_H( inp_ATOM *at, int cur_at, INCHI_MODE nMode )
{
    int i, j, next_at, num_stereo_bonds, bFound;
    int bond_type, num_2s, num_alt;
    int num_2s_next, num_alt_next, num_wrong_bonds_1, num_wrong_bonds_2;
#if ( N_V_STEREOBONDS == 1 )
    int n2sh, num_2s_hetero[2], num_2s_hetero_next[2], next_next_at, type_N, type_N_next;
#endif
    if ( MAX_NUM_STEREO_BOND_NEIGH < at[cur_at].valence+at[cur_at].num_H ||
         MIN_NUM_STEREO_BOND_NEIGH > at[cur_at].valence+at[cur_at].num_H  )
        return 0;
    if ( !bCanAtomHaveAStereoBond( at[cur_at].elname, at[cur_at].charge, at[cur_at].radical ) )
        return 0;
    /*  count bonds and find the second atom on the stereo bond */
    num_2s = num_alt = num_wrong_bonds_1 = 0;
#if ( N_V_STEREOBONDS == 1 )
    num_2s_hetero[0] = num_2s_hetero[1] = type_N = 0;
    if ( 0 == at[cur_at].num_H && 0 == at[cur_at].charge && 0 == at[cur_at].radical &&
         3 == get_endpoint_valence( at[cur_at].el_number ) ) {
        if ( 2 == at[cur_at].valence && 3 == at[cur_at].chem_bonds_valence ) {
            type_N = 1;
        } else
        if ( 3 == at[cur_at].valence && 5 == at[cur_at].chem_bonds_valence ) {
            type_N = 2; /* unfortunately includes >N# */
        }
    }
#endif
    for ( i = 0, num_stereo_bonds = 0; i < at[cur_at].valence; i ++ ) {
        bFound    = 0;
        next_at   = at[cur_at].neighbor[i];
        bond_type = get_allowed_stereo_bond_type( (int)at[cur_at].bond_type[i] );
        if ( bond_type == BOND_ALTERN ) {
            num_alt ++;
            if ( cur_at > next_at && !(nMode & CMODE_NO_ALT_SBONDS) )
                bFound = 1;
        } else
        if ( bond_type == BOND_DOUBLE ) {
            num_2s ++;
#if ( N_V_STEREOBONDS == 1 )
            if ( 0 <= (n2sh = bIsSuitableHeteroInpAtom( at + next_at )) ) {
                num_2s_hetero[n2sh] ++; /* n2sh=0 -> =N- or =NH; n2sh=1 -> =O */
            }
#endif
            if ( cur_at > next_at )
                bFound = 1;
        } else
        if ( bond_type != BOND_SINGLE && bond_type != BOND_TAUTOM ) {
            num_wrong_bonds_1 ++;
#if ( ONE_BAD_SB_NEIGHBOR == 1 )
            if ( num_wrong_bonds_1 > 1 || num_wrong_bonds_1 && 2 >= at[cur_at].valence ) {
                return 0; /* wrong bond type */
            } else {
                continue;
            }
#else
            return 0; /*  wrong bond type */
#endif
        }

        if ( bFound ) {
            /*  check "next_at" atom on the opposite side of the bond */
            if ( MAX_NUM_STEREO_BOND_NEIGH < at[next_at].valence+at[next_at].num_H ||
                 MIN_NUM_STEREO_BOND_NEIGH > at[next_at].valence+at[next_at].num_H )
                continue;
            if ( !bCanAtomHaveAStereoBond( at[next_at].elname, at[next_at].charge, at[next_at].radical ) )
                continue;
            /*  next atom neighbors */
            num_2s_next = num_alt_next = num_wrong_bonds_2 = 0;
#if ( N_V_STEREOBONDS == 1 )
            num_2s_hetero_next[0] = num_2s_hetero_next[1] = type_N_next = 0;
            if ( 0 == at[next_at].num_H && 0 == at[next_at].charge && 0 == at[next_at].radical &&
                 3 == get_endpoint_valence( at[next_at].el_number ) ) {
                if ( 2 == at[next_at].valence && 3 == at[next_at].chem_bonds_valence ) {
                    type_N_next = 1; /* -N= */
                } else
                if ( 3 == at[next_at].valence && 5 == at[next_at].chem_bonds_valence ) {
                    type_N_next = 2; /* unfortunately includes >N# */
                }
            }
#endif
            for ( j = 0; j < at[next_at].valence; j ++ ) {
                bond_type = get_allowed_stereo_bond_type( (int)at[next_at].bond_type[j] );
                if ( bond_type == BOND_ALTERN )
                    num_alt_next ++;
                else
                if ( bond_type == BOND_DOUBLE ) {
                    num_2s_next ++;
#if ( N_V_STEREOBONDS == 1 )
                    next_next_at = at[next_at].neighbor[j];
                    if ( 0 <= (n2sh = bIsSuitableHeteroInpAtom( at + next_next_at )) ) {
                        num_2s_hetero_next[n2sh] ++; /* n2sh=0 -> =N- or =NH; n2sh=1 -> =O */
                    }
#endif
                } else
                if ( bond_type != BOND_SINGLE && bond_type != BOND_TAUTOM ) {
                    num_wrong_bonds_2 ++;
#if ( ONE_BAD_SB_NEIGHBOR == 1 )
                    if ( num_wrong_bonds_1 > 1 || num_wrong_bonds_1 && 2 >= at[cur_at].valence ) {
                        break; /* wrong bond type */
                    } else {
                        continue;
                    }
#else
                    break; /*  wrong bond type */
#endif
                }
            }
            /* figure out whether the at[cur_at]--at[next_at] bond may not be stereogenic */

#if ( N_V_STEREOBONDS == 1 )
            if ( 3 == (type_N | type_N_next) &&
                 ( 2 == type_N && !bIsOxide( at, cur_at ) ||
                   2 == type_N_next && !bIsOxide( at, next_at ) ) ) {
                bFound = 0;
            } else
#endif
            if ( j < at[next_at].valence ||                  /* at[next_at] has a wrong bond type*/
                 (num_alt_next>0) + (num_2s_next>0) != 1     /* only one type of stereogenic bond permitted */
                ) {
                bFound = 0;
            } else
            if ( 2 < num_2s_next ) {
                bFound = 0;
            } else
            if ( 2 == num_2s_next ) {
                if ( 2 == at[next_at].valence ) {
                    ; /* only one double bond permitted except cumulenes */
#if ( N_V_STEREOBONDS == 1 )
                } else
                if ( 1 == (num_2s_hetero_next[0] | num_2s_hetero_next[1]) &&
                     3 == at[next_at].valence + at[next_at].num_H &&
                     5 == at[next_at].chem_bonds_valence + at[next_at].num_H &&
                     3 == get_endpoint_valence( at[next_at].el_number ) &&
                     (!type_N || bIsOxide( at, next_at )) ) {
                    ; /*
                       *   found:
                       *
                       *    \      /    \      /    \      /     
                       *     \    /      \    /      \    /      
                       *      N==C   or   N==C   or   N==N       
                       *    //    \     //    \     //    \      
                       *   O  ^    \   N  ^    \   O  ^    \     
                       *      |           |           |          
                       *      |           |           |          
                       *      at[next_at] at[next_at] at[next_at]
                       */
#endif
                } else {
                    bFound = 0;
                }
            }

        }
        if ( bFound ) {
           num_stereo_bonds++;
        }
    }
    
    if ( (num_alt>0) + (num_2s>0) != 1 || !num_stereo_bonds )
        return 0;
    if ( num_2s > 1 ) {
#if ( N_V_STEREOBONDS == 1 )
        if ( 2 == num_2s &&
             1 == (num_2s_hetero[0] | num_2s_hetero[1]) &&
             3 == at[cur_at].valence + at[cur_at].num_H &&
             5 == at[cur_at].chem_bonds_valence + at[cur_at].num_H &&
             3 == get_endpoint_valence( at[cur_at].el_number ) ) {
            ;
        } else {
            return 0;
        }
#else
        return 0;
#endif
    }

    return num_stereo_bonds;
}
/*************************************************************/
int half_stereo_bond_action( int nParity, int bUnknown, int bIsotopic, int vABParityUnknown )
{
#define AB_NEGATIVE 0x10
#define AB_UNKNOWN  0x20
    int nAction;

    if ( nParity == AB_PARITY_NONE )
        return AB_PARITY_NONE;
    
    /*  Unknown (type 1) in the parity value may come from the 'Either' single bond only */
    /*  Treat it as a known single bond geometry and unknown (Either) double bond */
    if ( nParity == vABParityUnknown /*AB_PARITY_UNKN*/ )
        nParity = AB_PARITY_ODD  | AB_UNKNOWN;
    if ( nParity == -vABParityUnknown /*AB_PARITY_UNKN*/ )
        nParity = AB_PARITY_ODD  | AB_UNKNOWN | AB_NEGATIVE;

    /*  make positive, replace AB_PARITY_EVEN with AB_PARITY_ODD  */
    if ( nParity < 0 )
        nParity = ((nParity == -AB_PARITY_EVEN)? AB_PARITY_ODD : (-nParity)) | AB_NEGATIVE;
    else
    if (nParity == AB_PARITY_EVEN)
        nParity = AB_PARITY_ODD;

    /*  Unknown (type 2): was detected in the double bond attribute */
    /*  (this 'unknown' came from 'Either' double bond) */
    /*  Treat both unknowns in the same way */
    if ( bUnknown )
        nParity |= AB_UNKNOWN;

    if ( bIsotopic ) {
        switch ( nParity ) {
        case AB_PARITY_ODD:
        case AB_PARITY_ODD | AB_NEGATIVE:
            nAction = AB_PARITY_CALC;
            break;
        case AB_PARITY_ODD  | AB_UNKNOWN:
        case AB_PARITY_UNDF | AB_UNKNOWN:
        case AB_PARITY_ODD  | AB_UNKNOWN | AB_NEGATIVE:
        case AB_PARITY_UNDF | AB_UNKNOWN | AB_NEGATIVE:
            nAction = vABParityUnknown /*AB_PARITY_UNKN*/;
            break;
        case AB_PARITY_IISO:
        case AB_PARITY_IISO | AB_UNKNOWN:
            nAction = AB_PARITY_NONE;
            break;
        case AB_PARITY_UNDF:
        case AB_PARITY_UNDF | AB_NEGATIVE:
            nAction = AB_PARITY_UNDF;
            break;
        default:
            nAction = -1; /*  program error */
        }
    } else {
        /*  Non-isotopic */
        switch ( nParity ) {
        case AB_PARITY_ODD:
            nAction = AB_PARITY_CALC;
            break;
        case AB_PARITY_ODD  | AB_UNKNOWN:
        case AB_PARITY_UNDF | AB_UNKNOWN:
            nAction = vABParityUnknown /*AB_PARITY_UNKN*/;
            break;
        /* case AB_PARITY_ODD  | AB_UNKNOWN | AB_NEGATIVE: */
        case AB_PARITY_UNDF:
            nAction = AB_PARITY_UNDF;
            break;
        case AB_PARITY_ODD  | AB_UNKNOWN | AB_NEGATIVE:
        case AB_PARITY_ODD  | AB_NEGATIVE:
        case AB_PARITY_IISO:
        case AB_PARITY_IISO | AB_UNKNOWN:
        case AB_PARITY_UNDF | AB_NEGATIVE:
        case AB_PARITY_UNDF | AB_UNKNOWN | AB_NEGATIVE:
            nAction = AB_PARITY_NONE;
            break;
        default:
            nAction = -1; /*  program error */
        }
    }
    return nAction;
#undef AB_NEGATIVE
#undef AB_UNKNOWN
}
/*************************************************************/
int set_stereo_bonds_parity( sp_ATOM *out_at, inp_ATOM *at, int at_1, inp_ATOM *at_removed_H, int num_removed_H,
  INCHI_MODE nMode, QUEUE *q, AT_RANK *nAtomLevel, S_CHAR *cSource, 
  AT_RANK min_sb_ring_size, int bPointedEdgeStereo, int vABParityUnknown )
{
    int j, k, next_at_1, i_next_at_1, i_next_at_2, at_2, next_at_2, num_stereo_bonds, bFound, bAllene;
    int bond_type, num_2s_1, num_alt_1;
    int num_2s_2, num_alt_2;
#if ( ONE_BAD_SB_NEIGHBOR == 1 )
    int num_wrong_bonds_1, num_wrong_bonds_2;
#endif
#if ( N_V_STEREOBONDS == 1 )
    int n2sh, num_2s_hetero[2], num_2s_hetero_next[2], next_next_at, type_N, type_N_next;
#endif
    int num_stored_stereo_bonds, num_stored_isotopic_stereo_bonds;
    int chain_length, num_chains, cur_chain_length;
    int all_at_2[MAX_NUM_STEREO_BONDS];
    int all_pos_1[MAX_NUM_STEREO_BONDS], all_pos_2[MAX_NUM_STEREO_BONDS];
    S_CHAR all_unkn[MAX_NUM_STEREO_BONDS];
    int /*at_1_parity, at_2_parity,*/ nUnknown, stop=0;
    
    /* at_1_parity = AB_PARITY_NONE; */ /*  do not know */
    /*  check valence */
    if ( MAX_NUM_STEREO_BOND_NEIGH < at[at_1].valence+at[at_1].num_H ||
         MIN_NUM_STEREO_BOND_NEIGH > at[at_1].valence+at[at_1].num_H  )
        return 0;
    if ( !bCanAtomHaveAStereoBond( at[at_1].elname, at[at_1].charge, at[at_1].radical ) )
        return 0;
    if ( at[at_1].c_point )
        return 0; /* rejects atoms that can lose or gain a (positive) charge. 01-24-2003 */

    /*  middle cumulene atoms, for example, =C=, should be ignored here */
    /*  only atoms at the ends of cumulene chains are considered. */
    if ( !at[at_1].num_H && 2 == at[at_1].valence &&
         BOND_DOUBLE == get_allowed_stereo_bond_type( (int)at[at_1].bond_type[0] ) &&
         BOND_DOUBLE == get_allowed_stereo_bond_type( (int)at[at_1].bond_type[1] ) ) {
        return 0;
    }

    /*  count bonds and find the second atom on the stereo bond */
    num_2s_1 = num_alt_1 = 0;
    chain_length      = 0;
    num_chains        = 0;
#if ( ONE_BAD_SB_NEIGHBOR == 1 )
    num_wrong_bonds_1 = 0;
#endif
#if ( N_V_STEREOBONDS == 1 )
    num_2s_hetero[0] = num_2s_hetero[1] = type_N = 0;
    if ( 0 == at[at_1].num_H && 0 == at[at_1].charge && 0 == at[at_1].radical &&
         3 == get_endpoint_valence( at[at_1].el_number ) ) {
        if ( 2 == at[at_1].valence && 3 == at[at_1].chem_bonds_valence ) {
            type_N = 1;
        } else
        if ( 3 == at[at_1].valence && 5 == at[at_1].chem_bonds_valence ) {
            type_N = 2; /* unfortunately includes >N# */
        }
    }
#endif
    for ( i_next_at_1 = 0, num_stereo_bonds = 0; i_next_at_1 < at[at_1].valence; i_next_at_1 ++ ) {
        nUnknown = (at[at_1].bond_stereo[i_next_at_1] == STEREO_DBLE_EITHER);
        bond_type = get_allowed_stereo_bond_type( (int)at[at_1].bond_type[i_next_at_1] );
        at_2 = -1; /* not found */
        if ( bond_type == BOND_ALTERN ||
             bond_type == BOND_DOUBLE ) {
            next_at_1 = at_2 = at[at_1].neighbor[i_next_at_1];
            next_at_2 = at_1;
        }
        switch ( bond_type ) {
        case BOND_ALTERN:
            num_alt_1 ++;
#if ( FIND_RING_SYSTEMS == 1 )
            if ( at[at_1].nRingSystem != at[at_2].nRingSystem )
                continue; /* reject alt. bond connecting different ring systems */
#endif
            if ( (nMode & CMODE_NO_ALT_SBONDS) ||
                 !bCanAtomHaveAStereoBond( at[at_2].elname, at[at_2].charge, at[at_2].radical ) ) {
                continue; /*  reject non-stereogenic bond to neighbor ord. #i_next_at_1 */
            }
            break;
        case BOND_DOUBLE:
            /*  check for cumulene/allene */
            num_2s_1++;
            cur_chain_length = 0;
            if ( bCanAtomBeTerminalAllene( at[at_1].elname, at[at_1].charge, at[at_1].radical ) ) {
                /*
                 * Example of cumulene
                 * chain length = 2:     >X=C=C=Y<
                 *                        | | | |
                 *  1st cumulene atom= at_1 | | at_2 =last cumlene chain atom
                 *  next to at_1=   next_at_1 next_at_2  =previous to at_2
                 *
                 *  chain length odd:  stereocenter on the middle atom ( 1=> allene )
                 *  chain length even: "long stereogenic bond"
                 */
                while ((bAllene = 
                        !at[at_2].num_H && at[at_2].valence == 2 &&
                        BOND_DOUBLE == get_allowed_stereo_bond_type( (int)at[at_2].bond_type[0] )  &&
                        BOND_DOUBLE == get_allowed_stereo_bond_type( (int)at[at_2].bond_type[1] )) &&          
                        bCanAtomBeMiddleAllene( at[at_2].elname, at[at_2].charge, at[at_2].radical ) ) {
                    k = ((int)at[at_2].neighbor[0]==next_at_2); /*  opposite neighbor position */
                    next_at_2 = at_2;
                    nUnknown += (at[at_2].bond_stereo[k] == STEREO_DBLE_EITHER);
                    at_2 = (int)at[at_2].neighbor[k];
                    cur_chain_length ++;  /*  count =C= atoms */
                }
                if ( cur_chain_length ) {
                    num_chains ++;
                    if ( bAllene /* at the end of the chain atom Y is =Y=, not =Y< or =Y- */ ||
                         !bCanAtomBeTerminalAllene( at[at_2].elname, at[at_2].charge, at[at_2].radical ) ) {
                        cur_chain_length = 0;
                        continue; /*  ignore: does not fit cumulene description; go to check next at_1 neighbor */
                    }
                    chain_length = cur_chain_length; /*  accept a stereogenic cumulele */
                }
            }
#if ( N_V_STEREOBONDS == 1 )
            if ( !cur_chain_length &&
                 0 <= (n2sh = bIsSuitableHeteroInpAtom( at + at_2 )) ) {
                num_2s_hetero[n2sh] ++; /* n2sh=0 -> =N- or =NH; n2sh=1 -> =O */
            }
#endif
            if ( !cur_chain_length &&
                 !bCanAtomHaveAStereoBond( at[at_2].elname, at[at_2].charge, at[at_2].radical ) ) {
                    continue; /*  reject non-stereogenic bond to neighbor #i_next_at_1 */
            }

            break;

        case BOND_SINGLE:
        case BOND_TAUTOM:
            continue; /*  reject non-stereogenic bond to neighbor #i_next_at_1 */
        default:
#if ( ONE_BAD_SB_NEIGHBOR == 1 )
            num_wrong_bonds_1 ++;
            continue;
#else
            return 0; /*  wrong bond type; */
#endif
        }

        /*  check atom at the opposite end of possibly stereogenic bond */

        bFound   = (at_2 >= 0 && at_1 > at_2 ); /*  i_next_at_1 = at_1 stereogenic bond neighbor attachment number */

        if ( bFound ) {
            /*  check "at_2" atom on the opposite side of the bond or cumulene chain */
            if ( MAX_NUM_STEREO_BOND_NEIGH < at[at_2].valence+at[at_2].num_H ||
                 MIN_NUM_STEREO_BOND_NEIGH > at[at_2].valence+at[at_2].num_H )
                continue;
            
            /*  check at_2 neighbors and bonds */
            num_2s_2 = num_alt_2 = 0;
#if ( N_V_STEREOBONDS == 1 )
            num_2s_hetero_next[0] = num_2s_hetero_next[1] = type_N_next = 0;
            if ( 0 == at[at_2].num_H && 0 == at[at_2].charge && 0 == at[at_2].radical &&
                 3 == get_endpoint_valence( at[at_2].el_number ) ) {
                if ( 2 == at[at_2].valence && 3 == at[at_2].chem_bonds_valence ) {
                    type_N_next = 1; /* -N= */
                } else
                if ( 3 == at[at_2].valence && 5 == at[at_2].chem_bonds_valence ) {
                    type_N_next = 2; /* unfortunately includes >N# */
                }
            }
#endif
            i_next_at_2 = -1;  /*  unassigned mark */
#if ( ONE_BAD_SB_NEIGHBOR == 1 )
            num_wrong_bonds_2 = 0;
#endif
            for ( j = 0; j < at[at_2].valence; j ++ ) {
                bond_type = get_allowed_stereo_bond_type( (int)at[at_2].bond_type[j] );
                if ( !bond_type ) {
#if ( ONE_BAD_SB_NEIGHBOR == 1 )
                    num_wrong_bonds_2 ++;
                    continue;  /*  this bond type is not allowed to be adjacent to a stereo bond */
#else
                    break;
#endif
                }
                if ( bond_type == BOND_DOUBLE ) {
                    num_2s_2 ++;
#if ( N_V_STEREOBONDS == 1 )
                    next_next_at = at[at_2].neighbor[j];
                    if ( 0 <= (n2sh = bIsSuitableHeteroInpAtom( at + next_next_at )) ) {
                        num_2s_hetero_next[n2sh] ++; /* n2sh=0 -> =N- or =NH; n2sh=1 -> =O */
                    }
#endif
                } else {
                    num_alt_2  += ( bond_type == BOND_ALTERN );
                }
                if ( (int)at[at_2].neighbor[j] == next_at_2 )
                    i_next_at_2 = j; /*  assigned */
            }
            if ( 
#if ( ONE_BAD_SB_NEIGHBOR == 1 )
                 num_wrong_bonds_2 > 1 || num_wrong_bonds_2 && 2 >= at[at_2].valence ||
#else
                 j < at[at_2].valence /* "next" has a wrong bond type*/ ||
#endif
                 (num_alt_2>0) + (num_2s_2>0) != 1 || /* all double XOR all alt bonds only */
                  /* num_2s_2 > 1  ||*/ /* only one double bond permitted */
                  i_next_at_2 < 0 /* atom next to the opposite atom not found */ ) {
                bFound = 0;
            } else
            if ( at[at_2].c_point ) {
                bFound = 0; /* rejects atoms that can lose or gain a (positive) charge. 01-24-2003 */
            } else
            if ( num_2s_2 > 2 ) {
                bFound = 0;
            } else
#if ( N_V_STEREOBONDS == 1 )
            if ( 3 == (type_N | type_N_next) &&
                 ( 2 == type_N && !bIsOxide( at, at_1 ) ||
                   2 == type_N_next && !bIsOxide( at, at_2 ) ) ) {
                bFound = 0;
            } else
#endif
            if ( 2 == num_2s_2 ) {
#if ( N_V_STEREOBONDS == 1 )
                if ( !chain_length &&
                     1 == (num_2s_hetero_next[0] | num_2s_hetero_next[1]) &&
                     3 == at[at_2].valence + at[at_2].num_H &&
                     5 == at[at_2].chem_bonds_valence + at[at_2].num_H &&
                     3 == get_endpoint_valence( at[at_2].el_number ) &&
                     (!type_N || bIsOxide( at, at_2 )) ) {
                      /*
                       *   found:
                       *
                       *    \      /    \      /    \      /     
                       *     \    /      \    /      \    /      
                       *      N==C   or   N==C   or   N==N       
                       *    //    \     //    \     //    \      
                       *   O  ^    \   N  ^    \   O  ^    \     
                       *      |           |           |          
                       *      |           |           |          
                       *      at[at_2]    at[at_2]    at[at_2]
                       */
                    ;
                } else {
                    bFound = 0;
                }
#else
                bFound = 0;
#endif
            }


            if ( chain_length && num_alt_2 )
                return 0; /*  allow no alt bonds in cumulenes */
        }
        if ( bFound ) {
            all_pos_1[num_stereo_bonds]   = i_next_at_1; /* neighbor to at_1 position */
            all_pos_2[num_stereo_bonds]   = i_next_at_2; /* neighbor to at_2 position */
            all_at_2[num_stereo_bonds]    = at_2;        /* at_2 */
            all_unkn[num_stereo_bonds]    = nUnknown;    /* stereogenic bond has Unknown configuration */
            /*
            if ( (at[at_1].bUsed0DParity & 2) || (at[at_2].bUsed0DParity & 2) ) {
                for ( k = 0; k < MAX_NUM_STEREO_BONDS && at[at_1].sb_parity[k]; k ++ ) {
                    if ( at[at_1].sb_neigh[k] == i_next_at_1 ) {
                        if ( at[at_1].sb_parity[k] == AB_PARITY_UNKN && !nUnknown ) {
                            all_unkn[num_stereo_bonds] = 1;
                        }
                        break;
                    }
                }
            }
            */
            num_stereo_bonds ++;
        }
    }
    if ( num_chains > 1 ) {
        return 0; /*  cannot be more than 1 cumulene chain. */
    }
#if ( ONE_BAD_SB_NEIGHBOR == 1 )
    if ( num_wrong_bonds_1 > 1 || num_wrong_bonds_1 && 2 >= at[at_1].valence ) {
        return 0; /* wrong bond type */
    }
#endif
    /*  accept only short chains for now */
    /*  chain_length=1: >C=C=C<      tetrahedral center, allene */
    /*  chain_length=2: >C=C=C=C<    stereogenic bond, cumulene */
    if ( chain_length && (num_stereo_bonds != 1 || num_alt_1 || chain_length > MAX_CUMULENE_LEN) ) {
        return 0;
    }

    /*  we need 1 double bond/chain XOR up to 3 arom. bonds */
    /*  to have a stereogenic bond */
    if ( (num_alt_1>0) + (num_2s_1>0) != 1 || !num_stereo_bonds /*|| num_2s_1 > 1*/ )
        return 0;

    if ( num_2s_1 > 1 ) {
#if ( N_V_STEREOBONDS == 1 )
        if ( 2 == num_2s_1 &&
             2 == type_N   &&
             1 == (num_2s_hetero[0] | num_2s_hetero[1]) &&
             3 == at[at_1].valence + at[at_1].num_H &&
             5 == at[at_1].chem_bonds_valence + at[at_1].num_H &&
             3 == get_endpoint_valence( at[at_1].el_number ) ) {
            ;
        } else {
            return 0;
        }
#else
        return 0;
#endif
    }

    /* ================== calculate parities ====================== */
   
    
    /*  find possibly stereo bonds and save them */
    num_stored_isotopic_stereo_bonds = 0;
    num_stored_stereo_bonds          = 0;
    for ( k = 0; k < num_stereo_bonds; k ++ ) {
        
        int cur_parity, next_parity, abs_cur_parity, abs_next_parity, dot_prod_z;
        S_CHAR z_dir1[3], z_dir2[3]; /*  3D vectors for half stereo bond parity direction */
        int  chain_len_bits = MAKE_BITS_CUMULENE_LEN(chain_length);
        int  cur_parity_defined, next_parity_defined;
        int  cur_action, next_action, result_action;
        
        at_2 = all_at_2[k];
        i_next_at_1 = all_pos_1[k];

#if ( MIN_SB_RING_SIZE > 0 )
        if( at[at_1].nRingSystem == at[at_2].nRingSystem ) {
            /*  check min. ring size only if both double bond/cumulene */
            /*  ending atoms belong to the same ring system */
            j = is_bond_in_Nmax_memb_ring( at, at_1, i_next_at_1, q, nAtomLevel, cSource, min_sb_ring_size );
            if ( j > 0 ) {
                continue;
            } else
            if ( j < 0 ) {
                return CT_STEREOBOND_ERROR;
            }
        }
#endif

        i_next_at_2 = all_pos_2[k];
        nUnknown    = all_unkn[k];
        memset(z_dir1, 0, sizeof(z_dir1));
        memset(z_dir2, 0, sizeof(z_dir2));
        
        /********************************************************************************
         * find atom parities (negative means parity due to H-isotopes only)
         * and half stereo bond parity directions z_dir1, z_dir2.
         *
         * Bond can have unknown or undefined parity or no parity because of:
         * 1. Geometry (poorly defined, cannot calculate, for example linear =C-F
         *    or =CHD with no geometry) -- Undefined parity
         *                                                              H
         * 2. Identical H atoms (no parity in principle, for example =C<  )
         *    -- No parity                                              H
         *
         * 3. The user said double bond stereo is unknown
         *    or at least one of single bonds is in unknown direction
         *    -- Unknown parity
         *
         * These 3 cases (see above) are referred below as 1, 2, 3.
         * Each of the cases may be present or not (2 possibilities)
         * Total number of combination is 2*2*2=8
         *
         * Since a case when all 3 are not present is a well-defined parity,
         * we do not consider this case here. Then 2*2*2-1=7 cases are left.
         *
         * If several cases are present, list them below separated by "+".
         * For example, 1+2 means (1) undefined geometry and (2) no parity
         * is possible because of identical H atoms.
         *
         * N) Decision table, Non-isotopic, 2*2*2-1=7 cases:
         * =================================================
         * none     : 2+any: 1+2(e.g.=CH2); 1+2+3; 2; 2+3  AB_PARITY_NONE=0
         * undefined: 1                                    AB_PARITY_UNDF
         * unknown  : 1+3; 3                               AB_PARITY_UNKN
         *
         * I) Decision table, Isotopic, 2*2*2-1=7 cases:
         * =============================================
         * none     : none
         * undefined: 1; 1+2; 1+2+3; 2; 2+3
         * unknown  : 1+3; 3
         *
         * Note: When defining identical atoms H atoms in case 2,
         *       Isotopic and Non-isotopic cases are different:
         *  N: do NOT take into account the isotopic composition of H atoms
         *  I: DO take into account the isotopic composition of H atoms
         *     (it is assumed that H isotopes are always different)
         *
         * half_stereo_bond_parity() returns:
         * ==================================
         * Note: half_stereo_bond_parity() is unaware of case 3.
         *
         * can't be a half of a stereo bond    AB_PARITY_NONE
         * 1, isotopic & non-isotopic:         AB_PARITY_UNDF
         * 1, isotopic only                   -AB_PARITY_UNDF
         * 2, no parity: identical H isotopes  AB_PARITY_IISO
         * 3, 'Either' single bond(s)          AB_PARITY_UNKN ???
         * 3, 'Either' single bond(s), iso H  -AB_PARITY_UNKN ???
         * defined parity                      AB_PARITY_ODD,  AB_PARITY_EVEN
         * defined parity for isotopic only:  -AB_PARITY_ODD, -AB_PARITY_EVEN
         *
         * Resultant value for the stereo bond parity
         * ---+-------------------+-------+--------+----------------+
         * 3? | half_stereo_bond_ | N or I| case 1,| bond parity    |
         *    |  parity()=        |       | 2 or 3 |                |
         * ---+-------------------+-------+--------+----------------+
         *   ( AB_PARITY_ODD/EVEN) => N&I: -       => AB_PARITY_CALC (=6, calc.later)
         * 3+( AB_PARITY_ODD/EVEN) => N&I: 3       => AB_PARITY_UNKN (=3)
         *   (-AB_PARITY_ODD/EVEN) => N:   2       => AB_PARITY_NONE (=0)
         *   (-AB_PARITY_ODD/EVEN) => I:   -       => AB_PARITY_CALC
         * 3+(-AB_PARITY_ODD/EVEN) => N:   2+3     => AB_PARITY_UNDF (=4)
         * 3+(-AB_PARITY_ODD/EVEN) => I:   3       => AB_PARITY_UNKN
         *   ( AB_PARITY_IISO )    => N:   1+2, 2  => AB_PARITY_NONE (=0)
         *   ( AB_PARITY_IISO )    => I:   1+2, 2  => AB_PARITY_UNDF
         * 3+( AB_PARITY_IISO )    => N:  1+2+3,2+3=> AB_PARITY_NONE
         * 3+( AB_PARITY_IISO )    => I:  1+2+3,2+3=> AB_PARITY_UNDF
         *   ( AB_PARITY_UNDF )    => N&I: 1       => AB_PARITY_UNDF
         * 3+( AB_PARITY_UNDF )    => N&I: 1+3     => AB_PARITY_UNKN
         *   (-AB_PARITY_UNDF )    => N:   1+2     => AB_PARITY_NONE
         *   (-AB_PARITY_UNDF )    => I:   1       => AB_PARITY_UNDF
         * 3+(-AB_PARITY_UNDF )    => N:   1+2+3   => AB_PARITY_NONE
         * 3+(-AB_PARITY_UNDF )    => I:   1+3     => AB_PARITY_UNKN
         * ---+-------------------+-------+--------+----------------+
        
         * If bond parity is undefined because abs(dot_prod_z) < MIN_DOT_PROD
         * then replace: AB_PARITY_CALC 
         *         with: AB_PARITY_UNDF
         * Joining two half_bond_parity() results:
         * 
         *
         * atom1 \ atom2   | AB_PARITY_NONE  AB_PARITY_UNKN  AB_PARITY_UNDF  AB_PARITY_CALC
         * ----------------+---------------------------------------------------------------
         *0=AB_PARITY_NONE | AB_PARITY_NONE  AB_PARITY_NONE  AB_PARITY_NONE  AB_PARITY_NONE
         *3=AB_PARITY_UNKN |                 AB_PARITY_UNKN  AB_PARITY_UNKN  AB_PARITY_UNKN
         *4=AB_PARITY_UNDF |                                 AB_PARITY_UNDF  AB_PARITY_UNDF
         *6=AB_PARITY_CALC |                                                 AB_PARITY_CALC
         *
         * that is, take min out of the two
         *********************************************************************************/

        cur_parity  =  half_stereo_bond_parity( at, at_1, at_removed_H, num_removed_H, 
                                                z_dir1, bPointedEdgeStereo, vABParityUnknown );
        next_parity =  half_stereo_bond_parity( at, at_2, at_removed_H, num_removed_H, 
                                                z_dir2, bPointedEdgeStereo, vABParityUnknown );

        if ( RETURNED_ERROR(cur_parity) || RETURNED_ERROR(next_parity) ) {
            return CT_CALC_STEREO_ERR;
        }
        if ( (at[at_1].bUsed0DParity & FlagSB_0D) || (at[at_1].bUsed0DParity & FlagSB_0D) ) {
            FixSb0DParities( at, /* at_removed_H, num_removed_H,*/ chain_length,
                             at_1, i_next_at_1, z_dir1,
                             at_2, i_next_at_2, z_dir2, &cur_parity, &next_parity );
        }

        if ( cur_parity == AB_PARITY_NONE || abs(cur_parity) == AB_PARITY_IISO ) {
            continue;
        }
        if ( next_parity == AB_PARITY_NONE || abs(next_parity) == AB_PARITY_IISO ) {
            continue;
        }
            
        cur_action    = half_stereo_bond_action( cur_parity, nUnknown, 0, vABParityUnknown ); /*  -1 => program error */
        next_action   = half_stereo_bond_action( next_parity, nUnknown, 0, vABParityUnknown );
        result_action = inchi_min(cur_action, next_action);

        if ( result_action == -1 ) {
            stop = 1; /*  program error <BRKPT> */
        }

        abs_cur_parity   = abs( cur_parity );
        abs_next_parity  = abs( next_parity );
        cur_parity_defined  = ATOM_PARITY_WELL_DEF(abs_cur_parity);
        next_parity_defined = ATOM_PARITY_WELL_DEF(abs_next_parity);
        
        
        if ( cur_parity_defined && next_parity_defined ) {
            /*  find how the whole bond parity depend on geometry */
            /*  if dot_prod_z < 0 then bond_parity := 3-bond_parity */
            /*  can be done only for a well-defined geometry */
            /*
            dot_prod_z  = (chain_len_bits & BIT_CUMULENE_CHI)? 
                           triple_prod_char( at, at_1, i_next_at_1, z_dir1, at_2, i_next_at_2, z_dir2 ) :
                           dot_prodchar3(z_dir1, z_dir2);
            */
            dot_prod_z  = (chain_len_bits && BOND_CHAIN_LEN(chain_len_bits)%2)? 
                           triple_prod_char( at, at_1, i_next_at_1, z_dir1, at_2, i_next_at_2, z_dir2 ) :
                           dot_prodchar3(z_dir1, z_dir2);
    
            if ( abs(dot_prod_z) < MIN_DOT_PROD ) {
                /*  The geometry is not well-defined. Eliminate AB_PARITY_CALC */
                result_action = inchi_min( result_action, AB_PARITY_UNDF );
            }
        } else {
            dot_prod_z = 0;
        }

        if ( result_action != AB_PARITY_NONE && result_action != -1 ) {
            /*  stereo, no isotopes (only positive) */
            if ( cur_parity > 0 && next_parity > 0 ) {
                if ( save_a_stereo_bond( dot_prod_z, result_action | chain_len_bits,
                                     at_1, i_next_at_1, out_at[at_1].stereo_bond_neighbor,
                                     out_at[at_1].stereo_bond_ord, out_at[at_1].stereo_bond_z_prod,
                                     out_at[at_1].stereo_bond_parity,
                                     at_2, i_next_at_2, out_at[at_2].stereo_bond_neighbor,
                                     out_at[at_2].stereo_bond_ord, out_at[at_2].stereo_bond_z_prod,
                                     out_at[at_2].stereo_bond_parity) ) {
                    if ( !out_at[at_1].parity ||
                         cur_parity_defined && !ATOM_PARITY_WELL_DEF(abs(out_at[at_1].parity)) ) {
                        out_at[at_1].parity = cur_parity;
                        memcpy( out_at[at_1].z_dir, z_dir1, sizeof(out_at[0].z_dir) );
                    }
                    if ( !out_at[at_2].parity ||
                         next_parity_defined && !ATOM_PARITY_WELL_DEF(abs(out_at[at_2].parity)) ) {
                        out_at[at_2].parity = next_parity;
                        memcpy( out_at[at_2].z_dir, z_dir2, sizeof(out_at[0].z_dir) );
                    }
                    out_at[at_1].bAmbiguousStereo |= at[at_1].bAmbiguousStereo;
                    out_at[at_2].bAmbiguousStereo |= at[at_2].bAmbiguousStereo;
                    num_stored_stereo_bonds ++;
                }
            }
        }
        
        /*  stereo + isotopic (all non-zero) */
        cur_action    = half_stereo_bond_action( cur_parity, nUnknown, 1, vABParityUnknown ); /*  -1 => program error */
        next_action   = half_stereo_bond_action( next_parity, nUnknown, 1, vABParityUnknown );
        result_action = inchi_min(cur_action, next_action);
        cur_parity  = abs_cur_parity;
        next_parity = abs_next_parity;
        if ( result_action != AB_PARITY_NONE && result_action != -1 ) {
            /*  stero, isotopic */
            if ( cur_parity > 0 && next_parity > 0 ) {
                if( save_a_stereo_bond( dot_prod_z, result_action | chain_len_bits,
                                     at_1, i_next_at_1, out_at[at_1].stereo_bond_neighbor2,
                                     out_at[at_1].stereo_bond_ord2, out_at[at_1].stereo_bond_z_prod2,
                                     out_at[at_1].stereo_bond_parity2,
                                     at_2, i_next_at_2, out_at[at_2].stereo_bond_neighbor2,
                                     out_at[at_2].stereo_bond_ord2, out_at[at_2].stereo_bond_z_prod2,
                                     out_at[at_2].stereo_bond_parity2) ) {
                    if ( !out_at[at_1].parity2 ||
                         cur_parity_defined && !ATOM_PARITY_WELL_DEF(abs(out_at[at_1].parity2)) ) {
                        out_at[at_1].parity2 = cur_parity /*| chain_len_bits*/;
                        if ( !out_at[at_1].parity ) {
                            memcpy( out_at[at_1].z_dir, z_dir1, sizeof(out_at[0].z_dir) );
                        }
                    }
                    if ( !out_at[at_2].parity2 || /* next line changed from abs(out_at[at_2].parity) 2006-03-05 */
                         next_parity_defined && !ATOM_PARITY_WELL_DEF(abs(out_at[at_2].parity2)) ) {
                        out_at[at_2].parity2 = next_parity /*| chain_len_bits*/;
                        if ( !out_at[at_2].parity ) {
                            memcpy( out_at[at_2].z_dir, z_dir2, sizeof(out_at[0].z_dir) );
                        }
                    }
                    out_at[at_1].bAmbiguousStereo |= at[at_1].bAmbiguousStereo;
                    out_at[at_2].bAmbiguousStereo |= at[at_2].bAmbiguousStereo;
                    num_stored_isotopic_stereo_bonds ++;
                }
            }
        } else
        if ( result_action == -1 ) {
            stop = 1; /*  program error? <BRKPT> */
        }

    }
    if ( stop ) {
        return CT_CALC_STEREO_ERR;
    }
    return /*num_stored_stereo_bonds+*/ num_stored_isotopic_stereo_bonds;
}


/*********************************************************************/
/*  if isotopic H, D, T added, can the atom be a stereo center? */
#if ( NEW_STEREOCENTER_CHECK == 1 )
/* int bCanInpAtomBeAStereoCenter( inp_ATOM *at, int cur_at ) */
int can_be_a_stereo_atom_with_isotopic_H( inp_ATOM *at, int cur_at, int bPointedEdgeStereo )
{
    int nNumNeigh;
    if ( (nNumNeigh = bCanInpAtomBeAStereoCenter( at, cur_at, bPointedEdgeStereo )) &&
         at[cur_at].valence + at[cur_at].num_H == nNumNeigh &&
         at[cur_at].num_H <= NUM_H_ISOTOPES
       ) {
        return 1;
    }
    return 0;
}
#else
int can_be_a_stereo_atom_with_isotopic_H( inp_ATOM *at, int cur_at )
{
    int j, ret = 0;
    if ( bCanAtomBeAStereoCenter( at[cur_at].elname, at[cur_at].charge, at[cur_at].radical ) &&
         at[cur_at].valence + at[cur_at].num_H == MAX_NUM_STEREO_ATOM_NEIGH &&
         at[cur_at].num_H < MAX_NUM_STEREO_ATOM_NEIGH
       ) {
        
        for ( j = 0, ret=1; ret && j < at[cur_at].valence; j ++ ) {
            if ( (at[cur_at].bond_type[j] & ~BOND_MARK_ALL) != BOND_SINGLE ) {
                ret = 0;
            }
        }
    }
    return ret;
}
#endif
/***************************************************************/
int GetStereocenter0DParity( inp_ATOM *at, int cur_at, int j1, AT_NUMB nSbNeighOrigAtNumb[], int nFlag )
{
    int parity = AB_PARITY_NONE;
    if ( at[cur_at].p_parity && (j1 == MAX_NUM_STEREO_ATOM_NEIGH-1 || j1 == MAX_NUM_STEREO_ATOM_NEIGH) ) {
        int i, num_trans_inp, num_trans_neigh;
        AT_NUMB nInpNeighOrigAtNumb[MAX_NUM_STEREO_ATOM_NEIGH];
        for ( i = 0; i < MAX_NUM_STEREO_ATOM_NEIGH; i ++ ) {
            nInpNeighOrigAtNumb[i] = at[cur_at].p_orig_at_num[i];
            if ( nInpNeighOrigAtNumb[i] == at[cur_at].orig_at_number ) {
                nInpNeighOrigAtNumb[i] = 0; /* lone pair or explicit H */
            }
        }
        num_trans_inp   = insertions_sort( nInpNeighOrigAtNumb, MAX_NUM_STEREO_ATOM_NEIGH, sizeof(nInpNeighOrigAtNumb[0]), comp_AT_NUMB );
        num_trans_neigh = insertions_sort( nSbNeighOrigAtNumb, j1, sizeof(nSbNeighOrigAtNumb[0]), comp_AT_NUMB );
        if ( j1 == MAX_NUM_STEREO_ATOM_NEIGH-1 ) {
            ; /*num_trans_neigh += j1;*/ /* the lone pair or implicit H is implicitly at the top of the list */
        }
        if ( !memcmp( nInpNeighOrigAtNumb + MAX_NUM_STEREO_ATOM_NEIGH-j1, nSbNeighOrigAtNumb, j1*sizeof(AT_NUMB) ) ) {
            if ( ATOM_PARITY_WELL_DEF(at[cur_at].p_parity) ) {
                parity = 2 - (num_trans_inp + num_trans_neigh + at[cur_at].p_parity) % 2;
            } else {
                parity = at[cur_at].p_parity;
            }
            at[cur_at].bUsed0DParity |= nFlag; /* 0D parity used for streocenter parity */
        }
    }
    return parity;
}

/***************************************************************
 * Get stereo atom parity for the current order of attachments
 * The result in at[cur_at].parity is valid for previously removed
 * explicit hydrogen atoms, including isotopic ones, that are located in at_removed_H[]
 * The return value is a calculated parity.
 */
#define ADD_EXPLICIT_HYDROGEN_NEIGH  1
#define ADD_EXPLICIT_LONE_PAIR_NEIGH 2
int set_stereo_atom_parity( sp_ATOM *out_at, inp_ATOM *at, int cur_at, inp_ATOM *at_removed_H,
                            int num_removed_H, int bPointedEdgeStereo, int vABParityUnknown)
{
    int    j, k, next_at, num_z, j1, nType, num_explicit_H, tot_num_iso_H, nMustHaveNumNeigh;
    int    num_explicit_iso_H[NUM_H_ISOTOPES+1]; /*  numbers of removed hydrogen atoms */
    int    index_H[MAX_NUM_STEREO_ATOM_NEIGH]; /*  cannot have more than 4 elements: 1 H, 1 D, 1 T atom(s) */
    double z, sum_xyz[3], min_sine, triple_product;
    double at_coord[MAX_NUM_STEREO_ATOM_NEIGH][3];
    double bond_len_xy[4], rmax=0.0, rmin=0.0;
    double at_coord_center[3];
    int    parity, bAmbiguous = 0, bAddExplicitNeighbor = 0, b2D = 0, n2DTetrahedralAmbiguity = 0;
    int    bIgnoreIsotopicH = (0 != (at[cur_at].cFlags & AT_FLAG_ISO_H_POINT));
    AT_NUMB nSbNeighOrigAtNumb[MAX_NUM_STEREO_ATOM_NEIGH];
    
    out_at[cur_at].parity  =
    out_at[cur_at].parity2 =
    out_at[cur_at].stereo_atom_parity  =
    out_at[cur_at].stereo_atom_parity2 = AB_PARITY_NONE;
    parity                             = AB_PARITY_NONE;

    memset(num_explicit_iso_H, 0, sizeof(num_explicit_iso_H));
    num_explicit_H = 0;

#if ( NEW_STEREOCENTER_CHECK == 1 )
    if ( !(nMustHaveNumNeigh = bCanInpAtomBeAStereoCenter( at, cur_at, bPointedEdgeStereo ) ) ||
         at[cur_at].num_H > NUM_H_ISOTOPES
       ) {
        goto exit_function;
    }
#else
    nMustHaveNumNeigh = MAX_NUM_STEREO_ATOM_NEIGH;
    if ( !bCanAtomBeAStereoCenter( at[cur_at].elname, at[cur_at].charge, at[cur_at].radical ) ||
         at[cur_at].valence + at[cur_at].num_H != nMustHaveNumNeigh ||
         at[cur_at].num_H > NUM_H_ISOTOPES
       ) {
        goto exit_function;
    }
    for ( j = 0; j < at[cur_at].valence; j ++ ) {
        if ( (at[cur_at].bond_type[j] & ~BOND_MARK_ALL) != BOND_SINGLE ) {
            goto exit_function;
        }
    }
#endif

    /*  numbers of isotopic H atoms */
    for ( j = 0, tot_num_iso_H = 0; j < NUM_H_ISOTOPES; j ++ ) {
        if ( at[cur_at].num_iso_H[j] > 1 ) {
            goto exit_function; /*  two or more identical hydrogen isotopic neighbors */
        }
        tot_num_iso_H += at[cur_at].num_iso_H[j];
    }
    if ( bIgnoreIsotopicH ) {
        tot_num_iso_H = 0; /* isotopic H considered subject to exchange => ignore isotopic */
    }
    /*  number of non-isotopic H atoms */
    if ( at[cur_at].num_H - tot_num_iso_H > 1 ) {
        goto exit_function; /*  two or more identical hydrogen non-isotopic neighbors */
    }
    
    /*  count removed explicit terminal hydrogens attached to at[cur_at]. */
    /*  the result is num_explicit_H. */
    /*  Removed hydrogens are sorted in increasing isotopic shift order */
    if ( at_removed_H && num_removed_H > 0 ) {
        for ( j = 0; j < num_removed_H; j ++ ) {
            if ( at_removed_H[j].neighbor[0] == cur_at ) {
                k = at_removed_H[j].iso_atw_diff;
                /*  iso_atw_diff values: H=>0, 1H=>1, D=2H=>2, T=3H=>3 */
                if ( k < 0 || k > NUM_H_ISOTOPES || bIgnoreIsotopicH )
                    k = 0; /*  treat wrong H isotopes as non-isotopic H */
                num_explicit_iso_H[k] ++;
                index_H[num_explicit_H++] = j;
            }
        }
    }

    /*  coordinates initialization */
    num_z = 0;
    sum_xyz[0] = sum_xyz[1] = sum_xyz[2] = 0.0;

    at_coord_center[0] =
    at_coord_center[1] =
    at_coord_center[2] = 0.0;

    /*  fill out stereo center neighbors coordinates */
    /*  and obtain the parity from the geometry */

    for ( k = 0, j1 = 0; k < 2; k ++ ) {
        switch( k ) {

        case 0:
            /*   add coordinates of removed hydrogens */
            for ( j = 0; j < num_explicit_H; j ++, j1 ++ ) {
                next_at = index_H[j];
                /*  use bond description located at removed_H atom */
                /*  minus sign at get_z_coord: at_removed_H[] contains bonds TO at[cur_at], not FROM it. */
                /*  Note: &at[(at_removed_H-at)+ next_at] == &at_removed_H[next_at] */
                z = -get_z_coord( at, (at_removed_H-at)+ next_at, 0 /*neighbor #*/, &nType, -(bPointedEdgeStereo & PES_BIT_POINT_EDGE_STEREO) );
                switch ( nType ) {
                case ZTYPE_EITHER:
                    parity = vABParityUnknown /*AB_PARITY_UNKN*/ ; /*  no parity: bond in "Either" direction. */
                    goto exit_function;
                case ZTYPE_UP:
                case ZTYPE_DOWN:
                    nType = -nType; /*  at_removed_H[] contains bonds TO the center, not from */
                    b2D ++;
                    /*  no break; here */
                case ZTYPE_3D:
                    num_z ++;
                }

                nSbNeighOrigAtNumb[j1] = at_removed_H[next_at].orig_at_number;
                at_coord[j1][0] = at_removed_H[next_at].x-at[cur_at].x;
                at_coord[j1][1] = at_removed_H[next_at].y-at[cur_at].y;
                bond_len_xy[j1] = len2(at_coord[j1]);
                /* bond_len_xy[j1] = sqrt(at_coord[j1][0]*at_coord[j1][0]+at_coord[j1][1]*at_coord[j1][1]); */
                at_coord[j1][2] = (nType==ZTYPE_3D?    z :
                                   nType==ZTYPE_UP?    bond_len_xy[j1] :
                                   nType==ZTYPE_DOWN? -bond_len_xy[j1] : 0.0 );
            }
            break;
        case 1:
            /*  add all coordinates of other neighboring atoms */
            for ( j = 0; j < at[cur_at].valence; j ++, j1 ++ ) {
                next_at = at[cur_at].neighbor[j];
                z = get_z_coord( at, cur_at, j, &nType, (bPointedEdgeStereo & PES_BIT_POINT_EDGE_STEREO) );
                switch ( nType ) {
                case ZTYPE_EITHER:
                    parity = vABParityUnknown /*AB_PARITY_UNKN*/; /*  unknown parity: bond in "Either" direction. */
                    goto exit_function;
                case ZTYPE_UP:
                case ZTYPE_DOWN:
                    b2D ++;
                case ZTYPE_3D:
                    num_z ++;
                }

                nSbNeighOrigAtNumb[j1] = at[next_at].orig_at_number;
                at_coord[j1][0] = at[next_at].x-at[cur_at].x;
                at_coord[j1][1] = at[next_at].y-at[cur_at].y;
                bond_len_xy[j1] = len2(at_coord[j1]);
                /* bond_len_xy[j1] = sqrt(at_coord[j1][0]*at_coord[j1][0]+at_coord[j1][1]*at_coord[j1][1]); */
                at_coord[j1][2] = (nType==ZTYPE_3D?    z :
                                   nType==ZTYPE_UP?    bond_len_xy[j1] :
                                   nType==ZTYPE_DOWN? -bond_len_xy[j1] : 0.0 );
            }
            break;
        }
    }
    /* j1 is the number of explicit neighbors (that is, all neighbors except implicit H) */

    b2D = (b2D == num_z && num_z);  /*  1 => two-dimensional */

    if ( MAX_NUM_STEREO_ATOM_NEIGH   != at[cur_at].valence+num_explicit_H &&
         MAX_NUM_STEREO_ATOM_NEIGH-1 != at[cur_at].valence+num_explicit_H ) {
        /*  not enough geometry data to find the central atom parity */
        if ( nMustHaveNumNeigh == at[cur_at].valence+at[cur_at].num_H &&
             at[cur_at].num_H > 1 ) {
            /*  only isotopic parity is possible; no non-isotopic parity */
            if ( parity == vABParityUnknown /*AB_PARITY_UNKN*/ ) {
                parity = -vABParityUnknown /*AB_PARITY_UNKN*/; /*  the user marked the center as "unknown" */
            } else {
                parity = -AB_PARITY_UNDF; /*  not enough geometry; only isotopic parity is possible */
            }
        } else {
            parity = AB_PARITY_NONE;      /*  not a stereocenter at all */
        }
        goto exit_function;
    }
    /*  make all vector lengths equal to 1; exit if too short. 9-10-2002 */
    for ( j = 0; j < j1; j ++ ) {
        z = len3( at_coord[j] );
        if ( z < MIN_BOND_LEN ) {
            /* bond length is too small: use 0D parities */
            if ( AB_PARITY_NONE == (parity = GetStereocenter0DParity( at, cur_at, j1, nSbNeighOrigAtNumb, FlagSC_0D )) ) {
                parity = AB_PARITY_UNDF;
            }
            goto exit_function;
        }
#if ( STEREO_CENTER_BONDS_NORM == 1 )
        else {
            mult3( at_coord[j], 1.0/z, at_coord[j] );
        }
#endif
        rmax = j? inchi_max( rmax, z) : z;
        rmin = j? inchi_min( rmin, z) : z;
    }
    if ( rmin / rmax < MIN_SINE ) {
        /* bond ratio is too small: use 0D parities */
        if ( AB_PARITY_NONE == (parity = GetStereocenter0DParity( at, cur_at, j1, nSbNeighOrigAtNumb, FlagSC_0D )) ) {
            parity = AB_PARITY_UNDF;
        }
        goto exit_function;
    }
    for ( j = 0; j < j1; j ++ ) {
        add3( sum_xyz, at_coord[j], sum_xyz );
    }



    /*  here j1 is a number of neighbors including explicit terminal isotopic H */
    /*  num_explicit_iso_H[0] = number of explicit non-isotopic hydrogen atom neighbors */
    j = j1;
    /*  Add Explicit Neighbor */
    if ( j1 == MAX_NUM_STEREO_ATOM_NEIGH-1 ) {
        /*  add an explicit neighbor if possible */
        if ( nMustHaveNumNeigh == MAX_NUM_STEREO_ATOM_NEIGH-1 ) {
            bAddExplicitNeighbor = ADD_EXPLICIT_LONE_PAIR_NEIGH;
        } else
        if ( nMustHaveNumNeigh == MAX_NUM_STEREO_ATOM_NEIGH ) {
            /*  check whether an explicit non-isotopic hydrogen can be added */
            /*  to an atom that is a stereogenic atom */
            if ( 1 == at[cur_at].num_H - num_explicit_H &&     /*  the atom has only one one implicit hydrogen */
                 1 == at[cur_at].num_H - tot_num_iso_H ) {     /*  this hydrogen is non-isotopic */
                bAddExplicitNeighbor = ADD_EXPLICIT_HYDROGEN_NEIGH;
            }
        }
    }

    if ( bAddExplicitNeighbor ) {
        /***********************************************************
         * May happen only if (j1 == MAX_NUM_STEREO_ATOM_NEIGH-1)
         * 3 neighbors only, no H-neighbors. Create and add coordinates of an implicit H
         * or a fake 4th neighbor, that is, a lone pair 
         */
        if ( parity == vABParityUnknown /*AB_PARITY_UNKN*/ ) {
            goto exit_function;  /*  the user insists the parity is unknown and the isotopic */
                                 /*  composition of the neighbors does not contradict */
        } else
        if ( num_z == 0 || are_3_vect_in_one_plane(at_coord, MIN_SINE) ) {
            /*  "hydrogen down" rule is needed to resolve an ambiguity */
            if ( num_z > 0 ) {
                bAmbiguous |= AMBIGUOUS_STEREO;
            }
#if ( APPLY_IMPLICIT_H_DOWN_RULE == 1 )  /* { */
            /*  Although H should be at the top of the list, add it to the bottom. */
            /*  This will be taken care of later by inverting parity 1<->2 */
            at_coord[j][0] = 0.0;
            at_coord[j][1] = 0.0;
#if ( STEREO_CENTER_BONDS_NORM == 1 )
            at_coord[j][2] = -1.0;
#else
            at_coord[j][2] = -(bond_len_xy[0]+bond_len_xy[1]+bond_len_xy[2])/3.0;
#endif
#else /* } APPLY_IMPLICIT_H_DOWN_RULE { */
#if (ALWAYS_SET_STEREO_PARITY == 1)
            parity =  AB_PARITY_EVEN; /*  suppose atoms are pre-sorted (testing) */
#else
            /* all 3 bonds are in one plain: try to get 0D parities */
            if ( AB_PARITY_NONE == (parity = GetStereocenter0DParity( at, cur_at, j1, nSbNeighOrigAtNumb, FlagSC_0D )) ) {
                parity = AB_PARITY_UNDF;
            }
            /*parity =  AB_PARITY_UNDF;*/ /*  no parity can be calculated found */
#endif
            goto exit_function;
#endif /* } APPLY_IMPLICIT_H_DOWN_RULE */
        } else {
            /*  we have enough information to find implicit hydrogen coordinates */
            /*
            at_coord[j][0] = -sum_x;
            at_coord[j][1] = -sum_y;
            at_coord[j][2] = -sum_z;
            */
            copy3( sum_xyz, at_coord[j] );
            change_sign3( at_coord[j], at_coord[j] );
            z = len3( at_coord[j] );
#if ( FIX_STEREO_SCALING_BUG == 1 )
            if ( z > 1.0 ) {
                rmax *= z;
            } else {
                rmin *= z;
            }
#else
            /* Comparing the original bond lengths to lenghts derived from normalized to 1 */
            /* This bug leads to pronouncing legitimate stereogenic atoms */
            /* connected by 3 bonds "undefined" if in a nicely drawn 2D structure */
            /* bond lengths are about 20 or greater. Reported by Reinhard Dunkel 2005-08-05 */
            if ( bPointedEdgeStereo & PES_BIT_FIX_SP3_BUG ) {
                /* coordinate scaling bug fixed here */
                if ( z > 1.0 ) {
                    rmax *= z;
                } else {
                    rmin *= z;
                }
            } else {
                /* original InChI v.1 bug */
                rmax = inchi_max( rmax, z );
                rmin = inchi_min( rmin, z );
            }
#endif
            if ( z < MIN_BOND_LEN || rmin/rmax < MIN_SINE ) {
                /* the new 4th bond is too short: try to get 0D parities */
                if ( AB_PARITY_NONE == (parity = GetStereocenter0DParity( at, cur_at, j1, nSbNeighOrigAtNumb, FlagSC_0D )) ) {
                    parity = AB_PARITY_UNDF;
                }
                goto exit_function;
            }
#if ( STEREO_CENTER_BOND4_NORM == 1 )
            else {
                mult3( at_coord[j], 1.0/z, at_coord[j] );
            }
#endif
        }
    } else
    if ( j1 != MAX_NUM_STEREO_ATOM_NEIGH ) {
        if ( parity == vABParityUnknown /*AB_PARITY_UNKN*/ ) {
            parity = -AB_PARITY_UNDF; /*  isotopic composition of H-neighbors contradicts 'unknown' */
        }
        goto exit_function;
    } else /*  j1 == MAX_NUM_STEREO_ATOM_NEIGH */
    if ( num_z == 0 || are_4at_in_one_plane(at_coord, MIN_SINE) ) {
        /*  all four neighours in xy plane: undefined geometry. */
        if ( num_z > 0 ) {
            bAmbiguous |= AMBIGUOUS_STEREO;
        }
        if ( parity != vABParityUnknown /*AB_PARITY_UNKN*/ ) {
#if (ALWAYS_SET_STEREO_PARITY == 1)
            parity =  AB_PARITY_EVEN; /*  suppose atoms are pre-sorted (testing) */
#else
            /* all 4 bonds are in one plain: try to get 0D parities */
            if ( AB_PARITY_NONE == (parity = GetStereocenter0DParity( at, cur_at, j1, nSbNeighOrigAtNumb, FlagSC_0D )) ) {
                parity = AB_PARITY_UNDF;
            } else
            if ( ATOM_PARITY_WELL_DEF( parity ) ) {
                bAmbiguous &= ~AMBIGUOUS_STEREO; /* 0D parity has resolved the ambiguity */
            }
#endif
        }
        goto exit_function;
    }
    /***********************************************************
     * At this point we have 4 neighboring atoms.
     * check for tetrahedral ambiguity in 2D case
     */
    if ( b2D ) 
    {

        n2DTetrahedralAmbiguity = Get2DTetrahedralAmbiguity( at_coord, bAddExplicitNeighbor, (bPointedEdgeStereo & PES_BIT_FIX_SP3_BUG ) );

        if ( 0 < n2DTetrahedralAmbiguity ) 
        {
            if ( T2D_WARN & n2DTetrahedralAmbiguity ) {
                bAmbiguous |= AMBIGUOUS_STEREO;
            }
            if ( T2D_UNDF & n2DTetrahedralAmbiguity ) {
                if ( parity != vABParityUnknown /*AB_PARITY_UNKN*/ ) {
#if (ALWAYS_SET_STEREO_PARITY == 1)
                    parity =  AB_PARITY_EVEN; /*  suppose atoms are pre-sorted (testing) */
#else
                    parity =  AB_PARITY_UNDF; /*  no parity */
#endif
                }
                goto exit_function;
            }
        } else
        if ( n2DTetrahedralAmbiguity < 0 ) {
            bAmbiguous |= AMBIGUOUS_STEREO_ERROR; /*  error */
            parity = AB_PARITY_UNDF;
            goto exit_function;
        }
    }

    /************************************************************/
    /*  Move coordinates origin to the neighbor #0 */
    for ( j = 1; j < MAX_NUM_STEREO_ATOM_NEIGH; j ++ ) {
        diff3(at_coord[j], at_coord[0], at_coord[j]);
    }
    diff3(at_coord_center, at_coord[0], at_coord_center);

    /*
    for ( k = 0; k < 3; k++ ) {
        for ( j = 1; j < MAX_NUM_STEREO_ATOM_NEIGH; j ++ ) {
            at_coord[j][k] -= at_coord[0][k];
        }
        at_coord_center[k] -= at_coord[0][k];
    }
    */
    /********************************************************
     * find the central (cur_at) atom's parity
     * (orientation of atoms #1-3 when looking from #0)
     ********************************************************/
    triple_product = triple_prod_and_min_abs_sine2(&at_coord[1], at_coord_center, bAddExplicitNeighbor, &min_sine, &bAmbiguous);
    /*
     * check for tetrahedral ambiguity -- leave it out for now
     */
    if ( fabs(triple_product) > ZERO_FLOAT && (min_sine > MIN_SINE || fabs(min_sine) > ZERO_FLOAT && (n2DTetrahedralAmbiguity & T2D_OKAY ) ) ) {
         /* Even => sorted in correct order, Odd=>transposed */
        parity = triple_product > 0.0? AB_PARITY_EVEN : AB_PARITY_ODD;
        /* if ( num_explicit_H && at[cur_at].removed_H_parity % 2 )  */
              /* odd transposition of the removed implicit H */
        /*     out_at[cur_at].parity = 3 - out_at[cur_at].parity; */
        
        /*  moved; see below */
        /* out_at[cur_at].bAmbiguousStereo |= bAmbiguous; */
        /* at[cur_at].bAmbiguousStereo |= bAmbiguous; */

        /*  for 4 attached atoms, moving the implicit H from index=3 to index=0 */
        /*  can be done in odd number (3) transpositions: (23)(12)(01), which inverts the parity */
        if ( j1 == MAX_NUM_STEREO_ATOM_NEIGH-1 ) {
            parity = 3 - parity;
        }
    } else {
#if (ALWAYS_SET_STEREO_PARITY == 1)
        parity =  AT_PARITY_EVEN; /*  suppose atoms are pre-sorted (testing) */
#else
        if ( num_z > 0 ) {
            bAmbiguous |= AMBIGUOUS_STEREO;
        }
        parity =  AB_PARITY_UNDF; /*  no parity: 4 bonds are in one plane. */
#endif
    }
exit_function:

    if ( parity ) {
        out_at[cur_at].bAmbiguousStereo |= bAmbiguous;
        at[cur_at].bAmbiguousStereo |= bAmbiguous;
    }

    
    /*  non-isotopic parity */
    if ( at[cur_at].num_H > 1 || parity <= 0 )
        ; /*  no non-isotopic parity */
    else
        out_at[cur_at].parity = parity;

    /*  isotopic parity */
    if ( parity == -AB_PARITY_UNDF || parity == -vABParityUnknown /*AB_PARITY_UNKN*/ )
        parity = -parity;
    if ( parity < 0 )
        parity = AB_PARITY_NONE;
    out_at[cur_at].parity2 = parity;


    parity = PARITY_VAL(out_at[cur_at].parity);
    out_at[cur_at].stereo_atom_parity = ATOM_PARITY_WELL_DEF( parity )? AB_PARITY_CALC : parity;
    parity = PARITY_VAL(out_at[cur_at].parity2);
    out_at[cur_at].stereo_atom_parity2 = ATOM_PARITY_WELL_DEF( parity )? AB_PARITY_CALC : parity;
    /*
    out_at[cur_at].parity2 = out_at[cur_at].parity; // save for stereo + isotopic canon.
    if ( out_at[cur_at].parity ) {
        if ( num_explicit_H > 1 || j1 == MAX_NUM_STEREO_ATOM_NEIGH-1 && num_explicit_H ) {
            //              X   H      X
            // for example,  >C<   or   >C-D
            //              Y   D      Y
            // parity exists for stereo + isotopic atoms canonicalization only
            out_at[cur_at].parity  = 0;
        }
    }
    // returning 0 means this can be an adjacent to a stereogenic bond atom
    */
    
    return (int)out_at[cur_at].parity2;

}
#undef ADD_EXPLICIT_HYDROGEN_NEIGH
#undef ADD_EXPLICIT_LONE_PAIR_NEIGH

/*************************************************************/
int set_stereo_parity( inp_ATOM* at, sp_ATOM* at_output, int num_at, int num_removed_H,
                       int *nMaxNumStereoAtoms, int *nMaxNumStereoBonds, INCHI_MODE nMode,
                       int bPointedEdgeStereo, int vABParityUnknown )
{
    int num_3D_stereo_atoms=0;
    int num_stereo_bonds=0; /* added to fix allene stereo bug reported for FClC=C=CFCl by Burt Leland - 2009-02-05 DT */

    int i, is_stereo, num_stereo, max_stereo_atoms=0, max_stereo_bonds=0;
    QUEUE *q = NULL;
    AT_RANK *nAtomLevel = NULL;
    S_CHAR  *cSource    = NULL;
    AT_RANK min_sb_ring_size = 0;

    /**********************************************************
     *
     * Note: this parity reflects only relative positions of
     *       the atoms-neighbors and their ordering in the
     *       lists of neighbors.
     *
     *       To obtain the actual parity, the parity of a number 
     *       of neighbors transpositions (to obtain a sorted
     *       list of numbers assigned to the atoms) should be
     *       added.
     *
     **********************************************************/

    /*********************************************************************************

     An example of parity=1 for stereogenic center, tetrahedral asymmetric atom
                   
                       
                    
                  (1)
                   |
                   |
               [C] |
                   |
         (2)------(0)
                  /
                /
              /
            /
         (3)
         

     Notation: (n) is a tetrahedral atom neighbor; n is an index of a neighbor in
     the central_at->neighbor[] array : neighbor atom number is central_at->neighbor[n].

     (0)-(1), (0)-(2), (0)-(3) are lines connecting atom [C] neighbors to neighbor (0)
     (0), (1) and (2) are in the plane
     (0)-(3) is directed from the plain to the viewer
     [C] is somewhere between (0), (1), (2), (3)
     Since (1)-(2)-(3)  are in a clockwise order when looking from (0), parity is 2, or even;
     otherwise parity would be 1, or odd.

    **********************************************************************************
    
      Examples of a stereogenic bond.
    
      Notation:   [atom number], (index of a neighbor):
                  [1] and [2] are atoms connected by the stereogenic bond
                  numbers in () are indexes of neighbors of [1] or [2].
                  (12 x 16)z = z-component of [1]-[2] and [1]-[6] cross-product

                                     atom [1]                     atom [2]                   
     [8]                    [4]      prod01 = (12 x 16)z < 0      prod01 = (21 x 24)z < 0    
       \                    /        prod02 = (12 x 18)z > 0      prod02 = (21 x 25)z > 0    
        (2)               (1)        0 transpositions because     0 transpositions because   
          \              /           double bond is in 0 posit.   double bond is in 0 position
          [1]==(0)(0)==[2]           0 = (prod01 > prod02)        0 = (prod01 > prod02)
          /              \                                                                   
        (1)               (2)        result: parity = 2, even     result: parity=2, even                                
       /                    \                                                                
     [6]                    [5]                                                              
                                                                                             
                                                                                             
                                                                                             
                                     atom [1]                     atom [2]                   
     [8]                    [5]      prod01 = (12 x 18)z > 0      prod01 = (21 x 24)z > 0    
       \                    /        prod02 = (12 x 16)z < 0      prod02 = (21 x 25)z < 0    
        (0)               (2)        2 transpositions to move     1 transposition to move    
          \              /           at [2] from 2 to 0 pos.      at [1] from 1 to 0 position
          [1]==(2)(1)==[2]           1 = (prod01 > prod02)        1 = (prod01 > prod02)      
          /              \                                                                   
        (1)               (0)        result: parity = (1+2)       result: parity=(1+1)           
       /                    \        2-(1+2)%2 = 1, odd           2-(1+1)%2 = 2, even
     [6]                    [4]           
                                          

    ***********************************************************************************
    Note: atoms' numbers [1], [2], [4],... are not used to calculate parity at this
          point. They will be used for each numbering in the canonicalization.
    Note: parity=3 for a stereo atom means entered undefined bond direction
          parity=4 for an atom means parity cannot be determined from the given geometry
    ***********************************************************************************/
    
    if ( !at_output || !at ) {
        return -1;
    }

    /*  clear stereo descriptors */

    for( i = 0; i < num_at; i ++ ) {
        at_output[i].parity  = 0;
        at_output[i].parity2 = 0;
        memset(&at_output[i].stereo_bond_neighbor[0], 0, sizeof(at_output[0].stereo_bond_neighbor) );
        memset(&at_output[i].stereo_bond_neighbor2[0], 0, sizeof(at_output[0].stereo_bond_neighbor2) );
        memset(&at_output[i].stereo_bond_ord[0], 0, sizeof(at_output[0].stereo_bond_ord) );
        memset(&at_output[i].stereo_bond_ord2[0], 0, sizeof(at_output[0].stereo_bond_ord2) );
        memset(&at_output[i].stereo_bond_z_prod[0], 0, sizeof(at_output[0].stereo_bond_z_prod) );
        memset(&at_output[i].stereo_bond_z_prod2[0], 0, sizeof(at_output[0].stereo_bond_z_prod2) );
        memset(&at_output[i].stereo_bond_parity[0], 0, sizeof(at_output[0].stereo_bond_parity) );
        memset(&at_output[i].stereo_bond_parity2[0], 0, sizeof(at_output[0].stereo_bond_parity2) );
    }
    /*  estimate max numbers of stereo atoms and bonds if isotopic H are added */
    if ( nMaxNumStereoAtoms || nMaxNumStereoBonds ) {
        for( i = 0, num_stereo = 0; i < num_at; i ++ ) {
            int num;
            num = can_be_a_stereo_atom_with_isotopic_H( at, i, bPointedEdgeStereo );
            if ( num ) {
                max_stereo_atoms += num;
            } else
            if ( (num = can_be_a_stereo_bond_with_isotopic_H( at, i, nMode ) ) ) { /*  accept cumulenes */
                max_stereo_bonds += num;
            }
        }
        if ( nMaxNumStereoAtoms )
            *nMaxNumStereoAtoms = max_stereo_atoms;
        if ( nMaxNumStereoBonds )
            *nMaxNumStereoBonds = max_stereo_bonds;
    }
    /*  calculate stereo descriptors */
#if ( MIN_SB_RING_SIZE > 0 )
    min_sb_ring_size = (AT_RANK)(((nMode & REQ_MODE_MIN_SB_RING_MASK) >> REQ_MODE_MIN_SB_RING_SHFT) & AT_RANK_MASK);
    if ( min_sb_ring_size >= 3 ) {
        /* create BFS data structure for finding for each stereo bond its min. ring sizes */
        q          = QueueCreate( num_at+1, sizeof(qInt) );
        nAtomLevel = (AT_RANK*)inchi_calloc(sizeof(nAtomLevel[0]),num_at);
        cSource    = (S_CHAR *)inchi_calloc(sizeof(cSource[0]),num_at);
        if ( !q || !cSource || !nAtomLevel ) {
            num_3D_stereo_atoms = CT_OUT_OF_RAM;
            goto exit_function;
        }
    } else {
        min_sb_ring_size = 2;
    }
#endif
    /* main cycle: set stereo parities */
    for( i = 0, num_stereo = 0; i < num_at; i ++ ) 
    {
        is_stereo = set_stereo_atom_parity( at_output, at, i, at+num_at, num_removed_H,
                                            bPointedEdgeStereo, vABParityUnknown ) ;
        if ( is_stereo ) 
        {
            num_3D_stereo_atoms += ATOM_PARITY_WELL_DEF( is_stereo );
        } 
        else 
        {
            is_stereo = set_stereo_bonds_parity( at_output, at, i, at+num_at, num_removed_H, nMode,
                                       q, nAtomLevel, cSource, min_sb_ring_size, 
                                       bPointedEdgeStereo, vABParityUnknown );
            if ( RETURNED_ERROR( is_stereo ) ) {
                num_3D_stereo_atoms = is_stereo;
                break;
            }
            num_stereo_bonds += (is_stereo != 0); /* added to fix bug reported by Burt Leland - 2009-02-05 DT */
        }
        num_stereo += (is_stereo != 0);
        is_stereo = is_stereo;
    }

    /* added to fix bug reported by Burt Leland - 2009-02-05 DT */
    if ( max_stereo_atoms < num_3D_stereo_atoms && nMaxNumStereoAtoms )
            *nMaxNumStereoAtoms = num_3D_stereo_atoms;
    if ( max_stereo_bonds < num_stereo_bonds && nMaxNumStereoBonds )
            *nMaxNumStereoBonds = num_stereo_bonds;

    /*
    if ( (nMode & REQ_MODE_SC_IGN_ALL_UU ) 
    REQ_MODE_SC_IGN_ALL_UU
    REQ_MODE_SB_IGN_ALL_UU
    */

#if ( MIN_SB_RING_SIZE > 0 )
    if ( q ) {
        q = QueueDelete( q );
    }
    if ( nAtomLevel )
        inchi_free( nAtomLevel );
    if ( cSource )
        inchi_free( cSource );
exit_function:
#endif


    return num_3D_stereo_atoms;
}
/*****************************************************************
 *  Functions that disconnect bonds
 *
 *=== During Preprocessing ===
 *
 *  RemoveInpAtBond
 *  DisconnectMetalSalt (is not aware of bond parities)
 *  DisconnectAmmoniumSalt
 *
 *=== Before Normalization ===
 *
 *  remove_terminal_HDT
 *
 *=== During the Normalization ===
 *
 *  AddOrRemoveExplOrImplH
 *
 *****************************************************************/
int ReconcileAllCmlBondParities( inp_ATOM *at, int num_atoms, int bDisconnected )
{
    int i, ret = 0;
    S_CHAR *visited = (S_CHAR*) inchi_calloc( num_atoms, sizeof(*visited) );
    if ( !visited )
        return -1; /* out of RAM */
    for ( i = 0; i < num_atoms; i ++ ) {
        if ( at[i].sb_parity[0] && !visited[i] && !(bDisconnected && is_el_a_metal(at[i].el_number)) ) {
            if ( ret = ReconcileCmlIncidentBondParities( at, i, -1, visited, bDisconnected ) ) {
                break; /* error */
            }
        }
    }
    inchi_free ( visited );
    return ret;
}
/*****************************************************************/
int ReconcileCmlIncidentBondParities( inp_ATOM *at, int cur_atom, int prev_atom, S_CHAR *visited, int bDisconnected )
{
    /* visited = 0 or parity => atom has not been visited
                 10 + parity => currently is on the stack + its final parity
                 20 + parity => has been visited; is not on the stack anymore + its final parity */
    int i, j, nxt_atom, ret = 0, len;
    int icur2nxt,  icur2neigh;   /* cur atom neighbors */
    int inxt2cur,  inxt2neigh;   /* next atom neighbors */
    int cur_parity, nxt_parity;
    int cur_order_parity, nxt_order_parity, cur_sb_parity, nxt_sb_parity, bCurMask, bNxtMask;
    /* !(bDisconnected && is_el_a_metal(at[i].el_number) */
    
    if ( at[cur_atom].valence > MAX_NUM_STEREO_BONDS )
        return 0; /* ignore */
    
    if ( !at[cur_atom].sb_parity[0] )
        return 1; /* wrong call */
    
    if ( visited[cur_atom] >= 10 )
        return 2; /* program error */

    cur_parity = visited[cur_atom] % 10;

    visited[cur_atom] += 10;

    for ( i = 0; i < MAX_NUM_STEREO_BONDS && at[cur_atom].sb_parity[i]; i ++ ) {
        icur2nxt = (int)at[cur_atom].sb_ord[i];
        len = get_opposite_sb_atom( at, cur_atom, icur2nxt, &nxt_atom, &inxt2cur, &j );
        if ( !len ) {
            return 4; /* could not find the opposite atom: bond parity data error */
        }
        if ( nxt_atom == prev_atom )
            continue;
        if ( visited[nxt_atom] >= 20 )
            continue; /* back edge, second visit: ignore */
        if ( at[nxt_atom].valence > MAX_NUM_STEREO_BONDS )
            continue; /* may be treated only after metal disconnection */

        if ( bDisconnected && (at[cur_atom].sb_parity[i] & SB_PARITY_FLAG) ) {
            cur_sb_parity = (at[cur_atom].sb_parity[i] >> SB_PARITY_SHFT);
            bCurMask = 3 << SB_PARITY_SHFT;
        } else {
            cur_sb_parity = (at[cur_atom].sb_parity[i] & SB_PARITY_MASK);
            bCurMask = 3;
        }
        if ( bDisconnected && (at[nxt_atom].sb_parity[j] & SB_PARITY_FLAG) ) {
            nxt_sb_parity = (at[nxt_atom].sb_parity[j] >> SB_PARITY_SHFT);
            bNxtMask = 3 << SB_PARITY_SHFT;
        } else {
            nxt_sb_parity = (at[nxt_atom].sb_parity[j] & SB_PARITY_MASK);
            bNxtMask = 3;
        }
        


        if ( !ATOM_PARITY_WELL_DEF(cur_sb_parity) ||
             !ATOM_PARITY_WELL_DEF(nxt_sb_parity)  ) {
            if ( cur_sb_parity == nxt_sb_parity ) {
                continue;
                /*goto move_forward;*/ /* bypass unknown/undefined */
            }
            return 3; /* sb parities do not match: bond parity data error */
        }

        icur2neigh  = (int)at[cur_atom].sn_ord[i];
        inxt2neigh  = (int)at[nxt_atom].sn_ord[j];
        /* parity of at[cur_atom].neighbor[] premutation to reach this order: { next_atom, neigh_atom, ...} */

        /* 1. move next_atom  from position=icur2nxt to position=0 =>
         *         icur2nxt permutations
         * 2. move neigh_atom from position=inxt2neigh+(inxt2cur > inxt2neigh) to position=1 =>
         *         inxt2neigh+(inxt2cur > inxt2neigh)-1 permutations.
         * Note if (inxt2cur > inxt2neigh) then move #1 increments neigh_atom position
         * Note add 4 because icur2neigh may be negative due to isotopic H removal
         */
        cur_order_parity = (4+icur2nxt + icur2neigh + (icur2neigh > icur2nxt)) % 2;
        /* same for next atom: */
        /* parity of at[nxt_atom].neighbor[] premutation to reach this order: { cur_atom, neigh_atom, ...} */
        nxt_order_parity = (4+inxt2cur + inxt2neigh + (inxt2neigh > inxt2cur)) % 2;
        
        nxt_parity = visited[nxt_atom] % 10;
        
        if ( !cur_parity ) {
            cur_parity = 2 - (cur_order_parity + cur_sb_parity) % 2;
            visited[cur_atom] += cur_parity;
        } else
        if ( cur_parity != 2 - (cur_order_parity + cur_sb_parity) % 2 ) {

            /***** reconcile bond parities *****/

            /* Each bond parity is split into two values located at the end atoms.
               For T (trans) the values are (1,1) or (2,2)
               For C (cis)   the values are (1,2) or (2,1)
               The fact that one pair = another with inverted parities, namely
               Inv(1,1) = (2,2) and Inv(1,2) = (2,1), allows one to
               simultaneouly invert parities of the current bond end atoms
               (at[cur_atom].sb_parity[i], at[nxt_atom].sb_parity[j])
               so that the final current atom parity cur_parity
               calculated later in stereochemical canonicalization for
               each stereobond incident with the current atomis same.
               Achieving this is called here RECONCILIATION.
               If at the closure of an aromatic circuit the parities of
               next atom cannot be reconciled with already calculated then
               this function returns 5 (error).
            */
            
            at[cur_atom].sb_parity[i] ^= bCurMask;
            at[nxt_atom].sb_parity[j] ^= bNxtMask;
            cur_sb_parity ^= 3;
            nxt_sb_parity ^= 3;
        }
        
        if ( !nxt_parity ) {
            nxt_parity = 2 - (nxt_order_parity + nxt_sb_parity) % 2;
            visited[nxt_atom] += nxt_parity;
        } else
        if ( nxt_parity != 2 - (nxt_order_parity + nxt_sb_parity) % 2 ) {
            return 5; /* algorithm does not work for Mebius-like structures */
        }

/* move_forward: */
        if ( visited[nxt_atom] < 10 ) {
            ret = ReconcileCmlIncidentBondParities( at, nxt_atom, cur_atom, visited, bDisconnected );
            if ( ret ) {
                break;
            }
        }
    }
    visited[cur_atom] += 10; /* all bonds incident to the current atom have
                                been processed or an error occurred. */ 
    return ret;
}
/*****************************************************************/
int get_opposite_sb_atom( inp_ATOM *at, int cur_atom, int icur2nxt, int *pnxt_atom, int *pinxt2cur, int *pinxt_sb_parity_ord )
{
    AT_NUMB nxt_atom;
    int     j, len;
    
    len = 0;
    while ( len ++ < 20 ) { /* arbitrarily set cumulene length limit to avoid infinite loop */
        nxt_atom = at[cur_atom].neighbor[icur2nxt];
        for ( j = 0; j < MAX_NUM_STEREO_BONDS && at[nxt_atom].sb_parity[j]; j ++ ) {
            if ( cur_atom == at[nxt_atom].neighbor[(int)at[nxt_atom].sb_ord[j]] ) {
                /* found the opposite atom */
                *pnxt_atom = nxt_atom;
                *pinxt2cur = at[nxt_atom].sb_ord[j];
                *pinxt_sb_parity_ord = j;
                return len;
            }
        }
        if ( j ) {
            return 0; /* reached atom(s) with stereobond (sb) parity, the opposite atom has not been found */
        }
        if ( at[nxt_atom].valence == 2 && 2*BOND_TYPE_DOUBLE == at[nxt_atom].chem_bonds_valence ) {
            /* follow cumulene =X= path */
            icur2nxt = (at[nxt_atom].neighbor[0] == cur_atom);
            cur_atom = nxt_atom;
        } else {
            return 0; /* neither atom with a sb parity not middle cumulene could be reached */
        }
    }
    return 0; /* too long chain of cumulene was found */
}

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
#include <string.h>

/*^^^ */
/*#define CHECK_WIN32_VC_HEAP*/
#include "mode.h"

#if ( READ_INCHI_STRING == 1 )

#include "ichi.h"
#include "ichitime.h"

#include "inpdef.h"
#include "ichimain.h"
#include "ichierr.h"
#include "incomdef.h" 
#include "ichiring.h"
#include "extr_ct.h"
#include "ichitaut.h"
#include "ichinorm.h"
#include "util.h"

#include "ichicomp.h"
#include "ichister.h"

#include "ichi_bns.h"

#include "strutil.h"

#include "ichirvrs.h"


/**************************************************************************************************

  ChargeStruct fictitios structures MY_CONST CN_LIST cnList[*]
  ============================================================


  bond         flow      (+)  => Positive charge c-group
  -----------------      (-)  => Negative charge c-group
  Single        0        (+C) => Positive charge group for C, Si, Ge, Sn, Pb
  Double        1        (-C) => Negative charge group for C, Si, Ge, Sn, Pb
  Triple        2        (.)  => additional one unit of st_cap

  A) Interpretation:

  X-(-)  or X=(+)  => zero charge
  X=(-)  or X-(+)  => charge = -1 or +1, respectively

  B) Information to keep:
  
  ordering zero-based number of the edge
  to (+) or (-) from the Interpretation (A) section

  vCap  = vertex cap
  vFlow = vertex flow
  val   = number of edges incident to the vertex
  neigh = 1-based ordering number of the adjacent vertex; 0 => no more adjacent vertices
  cap   = cap of the edge to the adjacent vertex
  flow  = flow of the edge to the adjacent vertex

  atom (c-point) always has number 1
  c-group(s) always are the last vertices
  always adjacent_neigh_number > vertex_number, that is, neigh > vertex

  Contribution to the Total Charge:
  ----------------------------------
  edge_cap(+) - edge_flow(+) - edge_flow(-) - Delta(+) - Delta(-)

  where edge_cap(+)  is edge capacity to c-group (+);
        edge_flow(+) is edge capacity to c-group (?), (?)= (+) or (-);
        Delta(?) = st_cap(?) - st_floe(?) of the c-group vertex (?), (?)= (+) or (-);

***************************************************************************************************/
/**************************************************************************************************

 Important:
 vCap and vFlow      Note: since metal charge group (vert. 2-4) MUST be registered before
 marked with empty         the "metal flower" (5-8) all charge group vertex numbers are
 comments for vertices     less than metal flower vertices: (2,3,4) < (5,6,7,8)
 1 and 5(M) should         This MAY be neded for c-group enumeration. The order is: 
 be set separately         t-groups, c-groups, M-flower. All types BNS_VT_M_GROUP allow one to avoid duplications.

              3(+)           
               ||         (Metal)
               ||          \|/       init charge=0; MAX_METAL_CHARGE = 16
   4(-)  5(M)  2           -Fe-                     CAP(BOND_TO_BNS_VT_M_GROUP) = NUM_BONDS*CAP
      \   |   /             |
        \ | /
          1              X(V), V=valence
*/
MY_CONST C_NODE cnMe[5] = {
    /*  vertex type       vCap  vFlow val;  neigh cap flow; neigh cap flow; neigh cap flow    vertex */
    { {BNS_VERT_TYPE_ATOM,0/**/ ,0/**/,3}, {{ 2, 16,0, 0 },{ 4,  16,0, 0 }, { 0,  0,0, 0 }} },  /* 1 */
    { {BNS_VT_CHRG_STRUCT,16,   16,    2}, {{ 3, 16,0,16 },{ 0,   0,0, 0 }, { 0,  0,0, 0 }} },  /* 2 */
    { {BNS_VT_C_POS_M,    16,   16,    1}, {{ 0,  0,0, 0 },{ 0,   0,0, 2 }, { 0,  0,0, 0 }} },  /* 3 */
    { {BNS_VT_C_NEG_M,    0+16,  0,    1}, {{ 0,  0,0, 0 },{ 0,   0,0, 0 }, { 0,  0,0, 0 }} },  /* 4 */
    { {BNS_VT_M_GROUP,    0/**/ ,0/**/,3}, {{ 1,  3,0, 0 },{ 0,   0,0, 0 }, { 0,  0,0, 0 }} },  /* 5 */
};
/*
#define cn_bits_Me (-1)
*/

/**************************************************************************************************
     c=2     5(+.)
     _____  /          (PNPN)
    4=====3             ||||       init charge=0
 c=2 \   / c=1          -N-
       2                 |
      ||||
       1              X+(V), X(V+1), X+(V+2), X(V+3); V=valence
*/
MY_CONST C_NODE cnPNPN[5] = {
    /*  vertex type       vCap  vFlow val;  neigh cap flow; neigh cap flow; neigh cap flow    vertex */
    { {BNS_VERT_TYPE_ATOM, 3,    3,    1}, {{ 2,  3,0, 3 },{ 0,   0,0, 0 }, { 0,  0,0, 0 }} },  /* 1 */
    { {BNS_VT_CHRG_STRUCT, 3,    3,    3}, {{ 3,  1,0, 0 },{ 4,   2,0, 0 }, { 0,  0,0, 0 }} },  /* 2 */
    { {BNS_VT_CHRG_STRUCT, 2,    2,    3}, {{ 5,  1,0, 0 },{ 4,   2,0, 2 }, { 0,  0,0, 0 }} },  /* 3 */
    { {BNS_VT_CHRG_STRUCT, 2,    2,    2}, {{ 0,  0,0, 0 },{ 0,   0,0, 0 }, { 0,  0,0, 0 }} },  /* 4 */
    { {BNS_VT_C_POS,       1+1,  0,    1}, {{ 0,  0,0, 0 },{ 0,   0,0, 0 }, { 0,  0,0, 0 }} },  /* 5 */
};
/*
#define cn_bits_PNPN MAKE_CN_BITS(cn_bits_P, cn_bits_N, cn_bits_P, cn_bits_N)
*/
/**************************************************************************************************
             5(+)
     c=1   //          (NPNP)
    4=====3             ||||       init charge=0
 c=1 \   / c=2          -N-
       2                 |
      ||||
       1              X(V), X+(V+1), X(V+2), X+(V+3); V=valence
*/
MY_CONST C_NODE cnNPNP[5] = {
    /*  vertex type       vCap  vFlow val;  neigh cap flow; neigh cap flow; neigh cap flow    vertex */
    { {BNS_VERT_TYPE_ATOM, 3,    3,    1}, {{ 2,  3,0, 3 },{ 0,   0,0, 0 }, { 0,  0,0, 0 }} },  /* 1 */
    { {BNS_VT_CHRG_STRUCT, 3,    3,    3}, {{ 3,  2,0, 0 },{ 4,   1,0, 0 }, { 0,  0,0, 0 }} },  /* 2 */
    { {BNS_VT_CHRG_STRUCT, 2,    2,    3}, {{ 5,  1,0, 1 },{ 4,   1,0, 1 }, { 0,  0,0, 0 }} },  /* 3 */
    { {BNS_VT_CHRG_STRUCT, 1,    1,    2}, {{ 0,  0,0, 0 },{ 0,   0,0, 0 }, { 0,  0,0, 0 }} },  /* 4 */
    { {BNS_VT_C_POS,       0+1,  1,    1}, {{ 0,  0,0, 0 },{ 0,   0,0, 0 }, { 0,  0,0, 0 }} },  /* 5 */
};
/*
#define cn_bits_NPNP MAKE_CN_BITS(cn_bits_N, cn_bits_P, cn_bits_N, cn_bits_P)
*/
/********************* end new ********************************************************************/


/**************************************************************************************************
             5(+)
           //          (NPN)
    4=====3             |||       init charge=0
     \   /              -N-
       2                 |
      |||
       1              X(V), X+(V+1), X(V+2); V=valence
*/
MY_CONST C_NODE cnNPN[5] = {
    /*  vertex type       vCap  vFlow val;  neigh cap flow; neigh cap flow; neigh cap flow    vertex */
    { {BNS_VERT_TYPE_ATOM, 2,    2,    1}, {{ 2,  2,0, 2 },{ 0,   0,0, 0 }, { 0,  0,0, 0 }} },  /* 1 */
    { {BNS_VT_CHRG_STRUCT, 2,    2,    3}, {{ 3,  1,0, 0 },{ 4,   1,0, 0 }, { 0,  0,0, 0 }} },  /* 2 */
    { {BNS_VT_CHRG_STRUCT, 2,    2,    3}, {{ 5,  1,0, 1 },{ 4,   1,0, 1 }, { 0,  0,0, 0 }} },  /* 3 */
    { {BNS_VT_CHRG_STRUCT, 1,    1,    2}, {{ 0,  0,0, 0 },{ 0,   0,0, 0 }, { 0,  0,0, 0 }} },  /* 4 */
    { {BNS_VT_C_POS,       1,    1,    1}, {{ 0,  0,0, 0 },{ 0,   0,0, 0 }, { 0,  0,0, 0 }} },  /* 5 */
};
/*
#define cn_bits_NPN MAKE_CN_BITS(cn_bits_N, cn_bits_P, cn_bits_N, 0)
*/
/**************************************************************************************************
             5(+.)
            /          (PNP)
    4=====3             |||       init charge=0
     \   /             -Cl-
       2                /\
      |||
       1              X+(V), X(V+1), X+(V+2); V=valence
*/
MY_CONST C_NODE cnPNP[5] = {
    /*  vertex type       vCap  vFlow val;  neigh cap flow; neigh cap flow; neigh cap flow    vertex */
    { {BNS_VERT_TYPE_ATOM, 2,    2,    1}, {{ 2,  2,0, 2 },{ 0,   0,0, 0 }, { 0,  0,0, 0 }} },  /* 1 */
    { {BNS_VT_CHRG_STRUCT, 2,    2,    3}, {{ 3,  1,0, 0 },{ 4,   1,0, 0 }, { 0,  0,0, 0 }} },  /* 2 */
    { {BNS_VT_CHRG_STRUCT, 1,    1,    3}, {{ 5,  1,0, 0 },{ 4,   1,0, 1 }, { 0,  0,0, 0 }} },  /* 3 */
    { {BNS_VT_CHRG_STRUCT, 1,    1,    2}, {{ 0,  0,0, 0 },{ 0,   0,0, 0 }, { 0,  0,0, 0 }} },  /* 4 */
    { {BNS_VT_C_POS,       1+1,  0,    1}, {{ 0,  0,0, 0 },{ 0,   0,0, 0 }, { 0,  0,0, 0 }} },  /* 5 */
};
/*
#define cn_bits_PNP MAKE_CN_BITS(cn_bits_P, cn_bits_N, cn_bits_P, 0)
*/
/**************************************************************************************************
                  
                      (MNP)
                       \ /        init charge=0
                        N(.) 
   3(-)  2(+)          / \
      \ //
      1(.)           X-(V), X(V+1), X+(V+2); V=valence
*/
MY_CONST C_NODE cnMNP[3] = {
    /*  vertex type       vCap  vFlow val;  neigh cap flow; neigh cap flow; neigh cap flow    vertex */
    { {BNS_VERT_TYPE_ATOM, 2,    1,    2}, {{ 2,  1,0, 1 },{ 3,   1,0, 0 }, { 0,  0,0, 0 }} },  /* 1 */
    { {BNS_VT_C_POS,       1,    1,    1}, {{ 0,  0,0, 0 },{ 0,   0,0, 0 }, { 0,  0,0, 0 }} },  /* 2 */
    { {BNS_VT_C_NEG,       0+1,  0,    1}, {{ 0,  0,0, 0 },{ 0,   0,0, 0 }, { 0,  0,0, 0 }} }   /* 3 */
};                                      
/*                                        
#define cn_bits_MNP  MAKE_CN_BITS(cn_bits_M, cn_bits_N, cn_bits_P, 0)
*/
#ifdef NEVER
/************************** not used **************************************************************
                        (PNM)
 5(-)      4(+)         \\ /
   \       //             B(.)    init charge=0
    3     2              / \
     \\  /
       1(.)           X+(V), X(V+1), X+(V+2); V=valence
*/
MY_CONST C_NODE cnPNM[5] = {
    /*  vertex type       vCap  vFlow val;  neigh cap flow; neigh cap flow; neigh cap flow    vertex */
    { {BNS_VERT_TYPE_ATOM, 2,    1,    2}, {{ 2,  1,0, 0 },{ 3,   1,0, 1 }, { 0,  0,0, 0 }} },  /* 1 */
    { {BNS_VT_CHRG_STRUCT, 1,    1,    2}, {{ 4,  1,0, 1 },{ 0,   0,0, 0 }, { 0,  0,0, 0 }} },  /* 2 */
    { {BNS_VT_CHRG_STRUCT, 1,    1,    2}, {{ 5,  1,0, 0 },{ 0,   0,0, 0 }, { 0,  0,0, 0 }} },  /* 3 */
    { {BNS_VT_C_POS,       1,    1,    1}, {{ 0,  0,0, 0 },{ 0,   0,0, 0 }, { 0,  0,0, 0 }} },  /* 4 */
    { {BNS_VT_C_NEG,       0+1,  0,    1}, {{ 0,  0,0, 0 },{ 0,   0,0, 0 }, { 0,  0,0, 0 }} },  /* 5 */
};

#define cn_bits_PNM  MAKE_CN_BITS(cn_bits_P, cn_bits_N, cn_bits_M, 0)

#endif

/**************************************************************************************************
   4(-)   3(+.)        (PNM)
     \   /              |||        init charge=0
       2               --P--
      |||               / \
       1              X-(V), X(V+1), X+(V+2); V=valence
*/
MY_CONST C_NODE cnPNM[4] = {
    /*  vertex type       vCap  vFlow val;  neigh cap flow; neigh cap flow; neigh cap flow    vertex */
    { {BNS_VERT_TYPE_ATOM, 2,    2,    1}, {{ 2,  2,0, 2 },{ 0,   0,0, 0 }, { 0,  0,0, 0 }} },  /* 1 */
    { {BNS_VT_CHRG_STRUCT, 2,    2,    3}, {{ 3,  1,0, 0 },{ 4,   1,0, 0 }, { 0,  0,0, 0 }} },  /* 2 */
    { {BNS_VT_C_POS,       1+1,  0,    1}, {{ 0,  0,0, 0 },{ 0,   0,0, 0 }, { 0,  0,0, 0 }} },  /* 3 */
    { {BNS_VT_C_NEG,       0+1,  0,    1}, {{ 0,  0,0, 0 },{ 0,   0,0, 0 }, { 0,  0,0, 0 }} },  /* 4 */
}; /* explanaton of vCap:  ^ ^  */      
/*         additional dot:/   \ capacity of the edge to (+) or (-) vertex */
/*
#define cn_bits_PNM  MAKE_CN_BITS(cn_bits_P, cn_bits_N, cn_bits_M, 0)
*/
/**************************************************************************************************
           5(+C)
          //                        init charge=0
 6(-C)  4
     \ /               (EN)  E=either +1 or -1
      3                 |
      ||               -C(.)-
      2                 |
      |  
      1(.)            X-(V), X+(V), X(V+1); V=valence
*/
MY_CONST C_NODE cnEN[6] = {
    /*  vertex type       vCap  vFlow val;  neigh cap flow; neigh cap flow; neigh cap flow    vertex */
    { {BNS_VERT_TYPE_ATOM, 1,    0,    1}, {{ 2,  1,0, 0 },{ 0,   0,0, 0 }, { 0,  0,0, 0 }} },  /* 1 */
    { {BNS_VT_CHRG_STRUCT, 1,    1,    2}, {{ 3,  1,0, 1 },{ 0,   0,0, 0 }, { 0,  0,0, 0 }} },  /* 2 */
    { {BNS_VT_CHRG_STRUCT, 1,    1,    3}, {{ 4,  1,0, 0 },{ 6,   1,0, 0 }, { 0,  0,0, 0 }} },  /* 3 */
    { {BNS_VT_CHRG_STRUCT, 1,    1,    2}, {{ 5,  1,0, 1 },{ 0,   0,0, 0 }, { 0,  0,0, 0 }} },  /* 4 */
    { {BNS_VT_C_POS_C,     0+1,  1,    1}, {{ 0,  0,0, 0 },{ 0,   0,0, 0 }, { 0,  0,0, 0 }} },  /* 5 */
    { {BNS_VT_C_NEG_C,     0+1,  0,    1}, {{ 0,  0,0, 0 },{ 0,   0,0, 0 }, { 0,  0,0, 0 }} }   /* 6 */
};
/*
#define cn_bits_EN  MAKE_CN_BITS(cn_bits_P | cn_bits_M, cn_bits_N, 0, 0)
*/
/**************************************************************************************************
             5(-)
            /          (NMN)        init charge=0
    4=====3             |||
     \   /              -X-
       2                /\
      |||
       1              X(V), X-(V+1), X(V+2); V=valence
*/
MY_CONST C_NODE cnNMN[5] = {
    /*  vertex type       vCap  vFlow val;  neigh cap flow; neigh cap flow; neigh cap flow    vertex */
    { {BNS_VERT_TYPE_ATOM, 2,    2,    1}, {{ 2,  2,0, 2 },{ 0,   0,0, 0 }, { 0,  0,0, 0 }} },  /* 1 */
    { {BNS_VT_CHRG_STRUCT, 2,    2,    3}, {{ 3,  1,0, 0 },{ 4,   1,0, 0 }, { 0,  0,0, 0 }} },  /* 2 */
    { {BNS_VT_CHRG_STRUCT, 1,    1,    3}, {{ 5,  1,0, 0 },{ 4,   1,0, 1 }, { 0,  0,0, 0 }} },  /* 3 */
    { {BNS_VT_CHRG_STRUCT, 1,    1,    2}, {{ 0,  0,0, 0 },{ 0,   0,0, 0 }, { 0,  0,0, 0 }} },  /* 4 */
    { {BNS_VT_C_NEG,       0+1,  0,    1}, {{ 0,  0,0, 0 },{ 0,   0,0, 0 }, { 0,  0,0, 0 }} }   /* 5 */
};
/*
#define cn_bits_NMN  MAKE_CN_BITS(cn_bits_N, cn_bits_M, cn_bits_N, 0)
*/
/**************************************************************************************************
           4(+)
          //           (NE)  E=either +1 or -1  
  5(-)  3               ||
     \ /               -X-          init charge=0
      2                 |                       
      ||         
      1               X(V), X+(V+1), X-(V+1); V=valence
*/
MY_CONST C_NODE cnNE[5] = {
    /*  vertex type       vCap  vFlow val;  neigh cap flow; neigh cap flow; neigh cap flow    vertex */
    { {BNS_VERT_TYPE_ATOM, 1,    1,    1}, {{ 2,  1,0, 1 },{ 0,   0,0, 0 }, { 0,  0,0, 0 }} },  /* 1 */
    { {BNS_VT_CHRG_STRUCT, 1,    1,    3}, {{ 3,  1,0, 0 },{ 5,   1,0, 0 }, { 0,  0,0, 0 }} },  /* 2 */
    { {BNS_VT_CHRG_STRUCT, 1,    1,    2}, {{ 4,  1,0, 1 },{ 0,   0,0, 0 }, { 0,  0,0, 0 }} },  /* 3 */
    { {BNS_VT_C_POS,       0+1,  1,    1}, {{ 0,  0,0, 0 },{ 0,   0,0, 0 }, { 0,  0,0, 0 }} },  /* 4 */
    { {BNS_VT_C_NEG,       0+1,  0,    1}, {{ 0,  0,0, 0 },{ 0,   0,0, 0 }, { 0,  0,0, 0 }} }   /* 5 */
};
/*
#define cn_bits_NE  MAKE_CN_BITS(cn_bits_N, cn_bits_P | cn_bits_M, 0, 0)
*/
/**************************************************************************************************
6(-)         5(+)
  \        //          (NEN)
    4=====3             |||         init charge=0
     \   /              -X-
       2                 |
      |||
       1              X(V), X+(V+1), X-(V+1), X(V+2); V=valence
*/
MY_CONST C_NODE cnNEN[6] = {
    /*  vertex type       vCap  vFlow val;  neigh cap flow; neigh cap flow; neigh cap flow    vertex */
    { {BNS_VERT_TYPE_ATOM, 2,    2,    1}, {{ 2,  2,0, 2 },{ 0,   0,0, 0 }, { 0,  0,0, 0 }} },  /* 1 */
    { {BNS_VT_CHRG_STRUCT, 2,    2,    3}, {{ 3,  1,0, 0 },{ 4,   1,0, 0 }, { 0,  0,0, 0 }} },  /* 2 */
    { {BNS_VT_CHRG_STRUCT, 2,    2,    3}, {{ 5,  1,0, 1 },{ 4,   1,0, 1 }, { 0,  0,0, 0 }} },  /* 3 */
    { {BNS_VT_CHRG_STRUCT, 1,    1,    3}, {{ 6,  1,0, 0 },{ 0,   0,0, 0 }, { 0,  0,0, 0 }} },  /* 4 */
    { {BNS_VT_C_POS,       0+1,  1,    1}, {{ 0,  0,0, 0 },{ 0,   0,0, 0 }, { 0,  0,0, 0 }} },  /* 5 */
    { {BNS_VT_C_NEG,       0+1,  0,    1}, {{ 0,  0,0, 0 },{ 0,   0,0, 0 }, { 0,  0,0, 0 }} },  /* 6 */
};
/*
#define cn_bits_NEN  MAKE_CN_BITS(cn_bits_N, cn_bits_M | cn_bits_N, cn_bits_N, 0)
*/
/*=======================================================*/
/**************************************************************************************************
                        (NP)
                         ||
                        -X-         init charge=0
      2(+)               |
      || 
       1              X(V), X+(V+1); V=valence
*/
MY_CONST C_NODE cnNP[2] = {
    /*  vertex type       vCap  vFlow val;  neigh cap flow; neigh cap flow; neigh cap flow    vertex */
    { {BNS_VERT_TYPE_ATOM, 1,    1,    1}, {{ 2,  1,0, 1 },{ 0,   0,0, 0 }, { 0,  0,0, 0 }} },  /* 1 */
    { {BNS_VT_C_POS,       0+1,  1,    1}, {{ 0,  0,0, 0 },{ 0,   0,0, 0 }, { 0,  0,0, 0 }} },  /* 2 */
};
/*
#define cn_bits_NP  MAKE_CN_BITS(cn_bits_N, cn_bits_P, 0, 0)
*/
/**************************************************************************************************
                        (PN)
      3(+.)              ||         init charge=0  [because cap(+)-flow(+)-Delta=1-0-1=0]
       |                -X-
       2                 |
      || 
       1              X+(V), X(V+1); V=valence
*/
MY_CONST C_NODE cnPN[3] = {
    /*  vertex type       vCap  vFlow val;  neigh cap flow; neigh cap flow; neigh cap flow    vertex */
    { {BNS_VERT_TYPE_ATOM, 1,    1,    1}, {{ 2,  1,0, 1 },{ 0,   0,0, 0 }, { 0,  0,0, 0 }} },  /* 1 */
    { {BNS_VT_CHRG_STRUCT, 1,    1,    2}, {{ 3,  1,0, 0 },{ 0,   0,0, 0 }, { 0,  0,0, 0 }} },  /* 2 */
    { {BNS_VT_C_POS,       1+1,  0,    1}, {{ 0,  0,0, 0 },{ 0,   0,0, 0 }, { 0,  0,0, 0 }} },  /* 3 */
};                                      
/*                                        
#define cn_bits_PN  MAKE_CN_BITS(cn_bits_P, cn_bits_N, 0, 0)
*/
/**************************************************************************************************
                        (NM)
      3(-)               ||         init charge=0
       |                -X-
       2                 |
      || 
       1              X(V), X-(V+1); V=valence
*/
MY_CONST C_NODE cnNM[3] = {
    /*  vertex type       vCap  vFlow; val; neigh cap flow; neigh cap flow; neigh cap flow    vertex */
    { {BNS_VERT_TYPE_ATOM, 1,    1,    1}, {{ 2,  1,0, 1 },{ 0,   0,0, 0 }, { 0,  0,0, 0 }} },  /* 1 */
    { {BNS_VT_CHRG_STRUCT, 1,    1,    2}, {{ 3,  1,0, 0 },{ 0,   0,0, 0 }, { 0,  0,0, 0 }} },  /* 2 */
    { {BNS_VT_C_NEG,       0+1,  0,    1}, {{ 0,  0,0, 0 },{ 0,   0,0, 0 }, { 0,  0,0, 0 }} },  /* 3 */
};                                      
/*                                        
#define cn_bits_NM  MAKE_CN_BITS(cn_bits_N, cn_bits_M, 0, 0)
*/
/**************************************************************************************************
                        (MN)
                         | 
                        -X-         init charge=0
      2(-)               |
       | 
       1(.)           X-(V), X(V+1); V=valence
*/
MY_CONST C_NODE cnMN[2] = {
    /*  vertex type       vCap  vFlow val;  neigh cap flow; neigh cap flow; neigh cap flow    vertex */
    { {BNS_VERT_TYPE_ATOM, 1,    0,    1}, {{ 2,  1,0, 0 },{ 0,   0,0, 0 }, { 0,  0,0, 0 }} },  /* 1 */
    { {BNS_VT_C_NEG,       0+1,  0,    1}, {{ 0,  0,0, 0 },{ 0,   0,0, 0 }, { 0,  0,0, 0 }} },  /* 2 */
};
/*
#define cn_bits_MN  MAKE_CN_BITS(cn_bits_M, cn_bits_N, 0, 0)
*/
/**************************************************************************************************
                        (P)
                         | 
                        -X-         init charge=0
      2(+.)              |
       | 
       1              X+(V); V=valence; all chemical (real) bonds to X have cap=0
*/
MY_CONST C_NODE cnP_[2] = {
    /*  vertex type       vCap  vFlow val;  neigh cap flow; neigh cap flow; neigh cap flow    vertex */
    { {BNS_VERT_TYPE_ATOM, 0,    0,    1}, {{ 2,  1,1, 0 },{ 0,   0,0, 0 }, { 0,  0,0, 0 }} },  /* 1 */
    { {BNS_VT_C_POS,       1+1,  0,    1}, {{ 0,  0,0, 0 },{ 0,   0,0, 0 }, { 0,  0,0, 0 }} },  /* 2 */
};
/*
#define cn_bits_P_  MAKE_CN_BITS(cn_bits_P, 0, 0, 0)
*/
#ifdef NEVER
/**************************************************************************************************
                        (M)
                         | 
                        -X-         init charge=-1 on atom
      2(-)               |
       | 
       1(.)           X+(V); V=valence; all chemical (real) bonds to X have cap=0
*/
MY_CONST C_NODE cnM_[2] = {
    /*  vertex type       vCap  vFlow val;  neigh cap flow; neigh cap flow; neigh cap flow    vertex */
    { {BNS_VERT_TYPE_ATOM, 1,    0,    1}, {{ 2,  1,0, 0 },{ 0,   0,0, 0 }, { 0,  0,0, 0 }} },  /* 1 */
    { {BNS_VT_C_NEG,       0+1,  0,    1}, {{ 0,  0,0, 0 },{ 0,   0,0, 0 }, { 0,  0,0, 0 }} },  /* 2 */
};
#endif

MY_CONST C_NODE cnM_[1] = {
    /*  vertex type       vCap  vFlow val;  neigh cap flow; neigh cap flow; neigh cap flow    vertex */
    { {BNS_VERT_TYPE_ATOM, 0,    0,    0}, {{ 0,  0,0, 0 },{ 0,   0,0, 0 }, { 0,  0,0, 0 }} },  /* 1 */
};
/*
#define cn_bits_M_  MAKE_CN_BITS(cn_bits_M, 0, 0, 0)
*/
/**************************************************************************************************
                           
                           
                        -X-         init charge=0
                         |
         
       1              X(V); V=valence;
*/
MY_CONST C_NODE cnN_[1] = {
    /*  vertex type       vCap  vFlow val;  neigh cap flow; neigh cap flow; neigh cap flow    vertex */
    { {BNS_VERT_TYPE_ATOM, 0,    0,    0}, {{ 0,  0,0, 0 },{ 0,   0,0, 0 }, { 0,  0,0, 0 }} },  /* 1 */
};
/*
#define cn_bits_N_  MAKE_CN_BITS(cn_bits_N, 0, 0, 0)
*/
/**************************************************************************************************/


MY_CONST CN_LIST cnList[] = {
    {cnPNPN, cn_bits_PNPN, 0, sizeof(cnPNPN)/sizeof(cnPNPN[0])},    /*  0 */
    {cnNPNP, cn_bits_NPNP, 0, sizeof(cnNPNP)/sizeof(cnNPNP[0])},    /*  1 */
    {cnNPN,  cn_bits_NPN,  0, sizeof(cnNPN)/sizeof(cnNPN[0])},      /*  2 */
    {cnPNP,  cn_bits_PNP,  0, sizeof(cnPNP)/sizeof(cnPNP[0])},      /*  3 */
    {cnMNP,  cn_bits_MNP,  0, sizeof(cnMNP)/sizeof(cnMNP[0])},      /*  4 */
    {cnPNM,  cn_bits_PNM,  0, sizeof(cnPNM)/sizeof(cnPNM[0])},      /*  5 */
    {cnEN,   cn_bits_EN ,  0, sizeof(cnEN)/sizeof(cnEN[0])},        /*  6 */
    {cnNMN,  cn_bits_NMN,  0, sizeof(cnNMN)/sizeof(cnNMN[0])},      /*  7 */
    {cnNE,   cn_bits_NE ,  0, sizeof(cnNE)/sizeof(cnNE[0])},        /*  8 */
    {cnNEN,  cn_bits_NEN,  0, sizeof(cnNEN)/sizeof(cnNEN[0])},      /*  9 */
    {cnNP,   cn_bits_NP,   0, sizeof(cnNP)/sizeof(cnNP[0])},        /* 10 */
    {cnPN,   cn_bits_PN,   0, sizeof(cnPN)/sizeof(cnPN[0])},        /* 11 */
    {cnNM,   cn_bits_NM,   0, sizeof(cnNM)/sizeof(cnNM[0])},        /* 12 */
    {cnMN,   cn_bits_MN,   0, sizeof(cnMN)/sizeof(cnMN[0])},        /* 13 */
    {cnP_,   cn_bits_P_,   0, sizeof(cnP_)/sizeof(cnP_[0])},        /* 14 */
    {cnM_,   cn_bits_M_,  -1, sizeof(cnM_)/sizeof(cnM_[0])},        /* 15 */
    {cnN_,   cn_bits_N_,   0, sizeof(cnN_)/sizeof(cnN_[0])},        /* 16 */
    {cnMe,   cn_bits_Me,   0, sizeof(cnMe)/sizeof(cnMe[0])}         /* 17 */
};

#define cnListIndexMe (17)  /* index of {cnMe, cn_bits_Me,... } element of cnList[] */

int cnListNumEl = (int)(sizeof(cnList)/sizeof(cnList[0]));
/**********************************************************************/
void clear_t_group_info( T_GROUP_INFO *ti )
{
    if ( !ti ) {
        return;
    } else {
        T_GROUP   *t_group                     = ti->t_group;
        int       max_num_t_groups             = ti->max_num_t_groups;
        AT_NUMB   *tGroupNumber                = ti->tGroupNumber;
        int       num_t_groups                 = ti->num_t_groups;
        AT_NUMB   *nEndpointAtomNumber         = ti->nEndpointAtomNumber;
        int       nNumEndpoints                = ti->nNumEndpoints;
        AT_NUMB   *nIsotopicEndpointAtomNumber = ti->nIsotopicEndpointAtomNumber;
        int       nNumIsotopicEndpoints        = ti->nNumIsotopicEndpoints;
        memset( ti, 0, sizeof(*ti) );
        if ( t_group ) {
            memset( t_group, 0, sizeof(t_group[0])*max_num_t_groups );
        } else {
            max_num_t_groups = 0;
        }
        if ( tGroupNumber ) {
            memset( tGroupNumber, 0, sizeof(tGroupNumber[0])*num_t_groups );
        } else {
            num_t_groups = 0;
        }
        if ( nEndpointAtomNumber ) {
            memset( nEndpointAtomNumber, 0, sizeof(nEndpointAtomNumber[0])*nNumEndpoints );
        } else {
            nNumEndpoints = 0;
        }
        if ( nIsotopicEndpointAtomNumber ) {
            memset( nIsotopicEndpointAtomNumber, 0, sizeof(nIsotopicEndpointAtomNumber[0])*nNumIsotopicEndpoints );
        } else {
            nNumIsotopicEndpoints = 0;
        }
        ti->t_group                      = t_group;
        ti->max_num_t_groups             = max_num_t_groups;
        ti->tGroupNumber                 = tGroupNumber;
        ti->num_t_groups                 = num_t_groups;
        ti->nEndpointAtomNumber          = nEndpointAtomNumber;
        ti->nNumEndpoints                = nNumEndpoints;
        ti->nIsotopicEndpointAtomNumber  = nIsotopicEndpointAtomNumber;
        ti->nNumIsotopicEndpoints        = nNumIsotopicEndpoints;
    }
    return;
}
/******************************************************************************************************/
int GetTgroupInfoFromInChI( T_GROUP_INFO *ti, inp_ATOM *at, AT_NUMB *endpoint, INChI *pInChI )
{
    int ret, i, j, k, itg, num_atoms, len_tg, bIso, num_t_groups;
    AT_NUMB   *tGroupNumber      = NULL;
    AT_NUMB   *tSymmRank         = NULL;
    AT_NUMB   *tiSymmRank        = NULL;
    AT_NUMB   *tiGroupNumber     = NULL;

    ret       = 0;

    clear_t_group_info( ti );
    if ( pInChI && pInChI->lenTautomer > 1 && pInChI->nTautomer && pInChI->nTautomer[0] > 0 ) {
        num_atoms = pInChI->nNumberOfAtoms;
        bIso      = pInChI->IsotopicAtom && pInChI->nNumberOfIsotopicAtoms;
        num_t_groups = pInChI->nTautomer[0];
        len_tg = pInChI->lenTautomer - T_GROUP_HDR_LEN*pInChI->nTautomer[0] - 1; /* number of endpoints */

        /* allocation ti->t_group */
        if ( ti->max_num_t_groups != num_atoms/2+1 || !ti->t_group ) {
            ti->max_num_t_groups = num_atoms/2+1;
            if ( ti->t_group )
                inchi_free( ti->t_group );
            ti->t_group = (T_GROUP *)inchi_calloc( ti->max_num_t_groups, sizeof(ti->t_group[0]));
        }
        /* allocation ti->tGroupNumber */
        if ( ti->num_t_groups != num_t_groups || !ti->tGroupNumber ) {
            ti->num_t_groups = num_t_groups;
            if ( ti->tGroupNumber )
                inchi_free( ti->tGroupNumber );
            ti->tGroupNumber = (AT_NUMB *)inchi_calloc((ti->num_t_groups+1)*TGSO_TOTAL_LEN, sizeof(ti->tGroupNumber[0]));
        }
        /* allocation ti->tGroupNumber */
        if ( len_tg != ti->nNumEndpoints || !ti->nEndpointAtomNumber ) {
            ti->nNumEndpoints = len_tg;
            if ( ti->nEndpointAtomNumber )
                inchi_free( ti->nEndpointAtomNumber );
            ti->nEndpointAtomNumber = (AT_NUMB *)inchi_calloc(len_tg+1, sizeof(ti->nEndpointAtomNumber[0]));
        }
        
        
        /* check */
        if ( !ti->t_group || !ti->tGroupNumber || !ti->nEndpointAtomNumber ) {
            ret = RI_ERR_ALLOC;
            goto exit_function;
        }
        tGroupNumber      = ti->tGroupNumber;
        tSymmRank         = tGroupNumber + TGSO_SYMM_RANK   * ti->num_t_groups;  /*  equivalence; cannot restore */
        tiSymmRank        = tGroupNumber + TGSO_SYMM_IRANK  * ti->num_t_groups;
        tiGroupNumber     = tGroupNumber + TGSO_SYMM_IORDER * ti->num_t_groups;


        INCHI_HEAPCHK
        j = 1; /* index in pInChI->nTautomer[] */
        i = 0; /* index in ti->nEndpointAtomNumber[] */
        for ( itg = 0; itg < pInChI->nTautomer[0]; itg ++ ) {
            len_tg = pInChI->nTautomer[j]; /* t-group length not including pInChI->nTautomer[j] */
            ti->t_group[itg].num[0] = pInChI->nTautomer[j+1]+pInChI->nTautomer[j+2]; /* num mobile H & (-) */
            ti->t_group[itg].num[1] = pInChI->nTautomer[j+2]; /* num mobile (-) */
            tGroupNumber[itg] = tiGroupNumber[itg] = itg; /* index */
            ti->t_group[itg].nGroupNumber = /*tSymmRank[itg] = tiSymmRank[itg] =*/ itg+1; /* t-group number */
            j      += T_GROUP_HDR_LEN;   /* skip t-group header */
            len_tg -= T_GROUP_HDR_LEN-1;
            
            ti->t_group[itg].nNumEndpoints = len_tg;
            ti->t_group[itg].nFirstEndpointAtNoPos = i;

            for( ; 0 < len_tg --; j ++, i ++ ) {
                k = ti->nEndpointAtomNumber[i] = pInChI->nTautomer[j]-1;
                if ( at ) {
                    at[k].endpoint = itg+1;
                }
                if ( endpoint ) {
                    endpoint[k] = itg+1;
                }
            }
        }
        if ( i != ti->nNumEndpoints ) {
            ret = RI_ERR_PROGR;
        }
        INCHI_HEAPCHK
    }
exit_function:
    return ret;
}
/******************************************************************************************************/
int FillOutpStructEndpointFromInChI( INChI *pInChI, AT_NUMB **pEndpoint )
{
    int num_at = pInChI->nNumberOfAtoms;
    AT_NUMB *endpoint = *pEndpoint;
    int     itg, i, j, k, len_tg;

    if ( !endpoint && !(endpoint = (AT_NUMB*) inchi_malloc(num_at * sizeof(endpoint[0]) ) ) ) {
        return RI_ERR_ALLOC;
    }
    memset( endpoint, 0, num_at * sizeof(endpoint[0]) );
    if ( pInChI->lenTautomer <= 1 || !pInChI->nTautomer ) {
        goto exit_function;
    }
    j = 1; /* index in pInChI->nTautomer[] */
    i = 0; /* index in ti->nEndpointAtomNumber[] */
    for ( itg = 0; itg < pInChI->nTautomer[0]; itg ++ ) {
        len_tg = pInChI->nTautomer[j]; /* t-group length not including pInChI->nTautomer[j] */
        j      += T_GROUP_HDR_LEN;   /* skip t-group header */
        len_tg -= T_GROUP_HDR_LEN-1;
        /* ti->t_group[itg].nNumEndpoints = len_tg; */
        for( ; 0 < len_tg --; j ++, i ++ ) {
            k = pInChI->nTautomer[j]-1;
            endpoint[k] = itg+1;
        }
    }
exit_function:
    *pEndpoint = endpoint;
    return 0;
}


/************************************************************************************/
int cmp_charge_val( const void *a1, const void *a2 )
{
    const CHARGE_VAL *p1 = (const CHARGE_VAL *) a1;
    const CHARGE_VAL *p2 = (const CHARGE_VAL *) a2;
    int    diff;
    
    if ( (diff = (int)p1->nValence - (int)p2->nValence) )  /* smaller valence first */
        return diff;
    if ( (diff = abs((int)p1->nCharge) - abs((int)p2->nCharge)) ) /* smaller abs charge first */
        return diff;
    if ( (diff = (int)p2->nCharge - (int)p1->nCharge) ) /* (+) first, (-) second */
        return diff;
    return (int)p1->nValenceOrderingNumber - (int)p2->nValenceOrderingNumber;
}
/************************************************************************************/
int bMayBeACationInMobileHLayer( inp_ATOM *at, VAL_AT *pVA, int iat, int bMobileH )
{
    static const char szEl[] = "N;P;O;S;Se;Te;";
    static const char cVal[] = {4,4,3,3, 3, 3, 0};
    static char en[8];
    static int  ne;
    int    i, j, neigh;
    char   *p;
    if ( !bMobileH || !at[iat].num_H ) {
        return 1;
    }
    if ( !ne ) { /* one time initialization */
        const char *b, *e;
        int  len;
        char elname[ATOM_EL_LEN];
        for ( b = szEl; (e = strchr( b, ';')); b = e+1 ) {
            len = e-b;
            memcpy( elname, b, len );
            elname[len] = '\0';
            en[ne++] = get_periodic_table_number( elname );
        }
        en[ne] = '\0';
    }
    if ( (p = (char *)memchr( en, at[iat].el_number, ne )) ) {
        i = p - en;
        /* >B(-)< exception */
        if ( at[iat].valence + at[iat].num_H <= cVal[i] ) {
            for ( j = 0; j < at[iat].valence; j ++ ) {
                neigh = at[iat].neighbor[j];
                if ( at[neigh].valence == 4 && at[neigh].chem_bonds_valence == 4 && !at[neigh].num_H &&
                     pVA[neigh].cNumValenceElectrons == 3 && pVA[neigh].cPeriodicRowNumber == 1 ) {
                    return 1;
                }
            }
            return 0;
        }
    }
    return 1;
}
/************************************************************************************/
int clean_charge_val( CHARGE_VAL *pChargeVal, int len, inp_ATOM *atom, VAL_AT *pVA,
                      int iat, int bIsMetal, int bMobileH, AT_NUMB *endpoint )
{
    inp_ATOM *at      = atom + iat;
    int nPeriodicNum  = at->el_number;
    int num_bonds     = at->valence;
    int min_valence   = at->valence + at->num_H;
    /* in fixed-H case treat tautomeric -O as tautomeric to avoid #O(+) */
    int bTautomeric   = (at->endpoint != 0);
    int bFixedHTautomeric = !bMobileH && (endpoint && endpoint[iat] &&
                            pVA[iat].cNumValenceElectrons == 6 && 1==num_bonds &&
                            !at->num_H && !bIsMetal);
    /* int bIsMetal      = is_el_a_metal( nPeriodicNum );*/
    int bDoNotAddH    = do_not_add_H( nPeriodicNum );
    int nPeriod, nNumEqAbsCharges;
    int nNumValenceEl = get_sp_element_type( nPeriodicNum, &nPeriod ) - 1;

    int i, j;

    if ( !len )
        return len;

    insertions_sort( pChargeVal, len, sizeof(pChargeVal[0]), cmp_charge_val );
    /* metals -- very preliminary code */
    if ( bIsMetal && bDoNotAddH ) {
        /* keep the 1st found */
        return inchi_min( 1, len );
    }
    /* Mobile-H layer cannot have H on positively charged N, P (all IV), O, S, Se, Te (all III) */

    /*
    if ( abs( pChargeVal[0].nCharge ) > 1 && pChargeVal[0].nValence >= min_valence ) {
        return inchi_min( 1, len );
    }
    */
    nNumEqAbsCharges = 0;
    for ( i = j = 0; i < len && j < (nNumEqAbsCharges? 3+nNumEqAbsCharges:4); i ++ ) {
        /* for now accept only charge = 0, -1, +1 */
        if ( abs( pChargeVal[i].nCharge ) > 1 ) {
            continue;
        }
        if ( BOND_TYPE_TRIPLE + BOND_TYPE_DOUBLE * (min_valence - 1) < pChargeVal[i].nValence ) {
            continue; /* not more than one triple and the rest - double bonds per atom */
        }
        if ( (bTautomeric || (j && bFixedHTautomeric)) && pChargeVal[i].nCharge < 0 ) {
            continue; /* negative charge must be included in the tautomeric group */
        }
        if ( (bTautomeric || bFixedHTautomeric) && pChargeVal[i].nCharge > 0 ) {
            continue; /* positive charge for now cannot reach a tautomeric group */
        }
        if ( j && !bMayBeACationInMobileHLayer( atom, pVA, iat, bMobileH ) && pChargeVal[i].nCharge > 0 ) {
            if ( i+1 < len &&
                 pChargeVal[i].nValence == pChargeVal[i+1].nValence &&
                 pChargeVal[i].nCharge == -pChargeVal[i+1].nCharge ) {
                /* (-) if exists is always after (+) */
                i += 1; /* also skip the next element */
            }
            continue; /* in case of Mobile-H, a hydrogen cannot be on a (+)-charged heteroatom */
        }
        /* accept same valence opposite charges only for C and its group in Periodic Table */
        if ( j && !bTautomeric &&
             pChargeVal[i].nValence == pChargeVal[j-1].nValence &&
             pChargeVal[i].nCharge  == -pChargeVal[j-1].nCharge ) {
            if ( nNumValenceEl == VALUE_OCTET/2 && pChargeVal[i].nCharge && !nNumEqAbsCharges ) {
                pChargeVal[j ++] = pChargeVal[i];
                nNumEqAbsCharges ++;
            }
            continue;
        }
        /* do not accept valence=5 for neutral NHn in case of not Mobile-H 2005-01-26 ???? */
        if ( nNumValenceEl == 5 && nPeriod == 1 && at->num_H &&
             j && !bMobileH &&
             pChargeVal[i].nValence == 5 && !pChargeVal[i].nCharge ) {
            continue;
        }
        /* do not accept gaps in allowed valences */
        if ( j && pChargeVal[i].nValence > pChargeVal[j-1].nValence+1 ) {
            break;
        }
        pChargeVal[j ++] = pChargeVal[i];
    }
    len = j;
    if ( !nNumEqAbsCharges && num_bonds < 3 && len == 4 ) {
        len --; /* prohibit =S#  where # is a triple bond */
    }
    return len;


}
/************************************************************************************
int GetAtomRestoreInfo( inp_ATOM *atom, int iat, VAL_AT *pVArray )
 
    pVA->cDoNotAddH
    pVA->cMetal
    pVA->cNumValenceElectrons
    pVA->cPeriodicRowNumber
    pVA->cInitFreeValences
    pVA->cnListIndex        = index+1

return value:
     -1 => error
      0 => do not know what to do; leave the atom unchanged
      1 => success
*************************************************************************************/
int GetAtomRestoreInfo( inp_ATOM *atom, int iat, VAL_AT *pVArray, ICHICONST SRM *pSrm, int bMobileH, AT_NUMB *endpoint )
{
/* #defines from util.c */
#define MIN_ATOM_CHARGE        (-2)
#define MAX_ATOM_CHARGE         2
#define NEUTRAL_STATE          (-MIN_ATOM_CHARGE)  
#define NUM_ATOM_CHARGES       (MAX_ATOM_CHARGE - MIN_ATOM_CHARGE + 1)
#define MAX_NUM_VALENCES        5  /* max. number + 1 to provide zero termination */

    int i, j, j2, k, k2, charge, cur_charge, num_non_bonding_electrons;
    int nNumStates, nNumSelectedStates, num_H, num_bonds;
    int nOctetNeutralValenceExcess, nFirstNeutralValenceExcess;
    int nFoundNeutralValenceExcess, nFoundNeutralValenceOrdNumber;
    int nLastFoundValenceOrdNumber, nLastFoundValenceState;
    int cn_bits, cn_bits_array[5], len_cn_bits_array;
    inp_ATOM    *at  = atom+iat;
    VAL_AT      *pVA = pVArray + iat;
    int nPeriodicNum = at->el_number;
    int cur_chem_valence, cur_chem_valence_fixed, min_chem_valence, known_chem_valence;
    int metal_bonds_chem_valence, not_metal_bonds_chem_valence, alt_bonds_delta_valence, bonds_chem_valence, bond_type;
    CHARGE_VAL  ChargeVal[NUM_ATOM_CHARGES*MAX_NUM_VALENCES];

    memset( ChargeVal, 0, sizeof(ChargeVal) );

    pVA->cDoNotAddH = do_not_add_H( nPeriodicNum );  /* InChI never adds H to this atom */
    /*pVA->cMetal     = is_el_a_metal( nPeriodicNum );*/ /* the atom is a metal */
    
    /* count bonds to metal atoms; metals have already been marked */
    metal_bonds_chem_valence = not_metal_bonds_chem_valence = alt_bonds_delta_valence = 0;
    if ( pVA->cMetal ) {
        j = at->valence; /* all bonds to metal */
        for ( i = k = j2 = k2 = 0; i < at->valence; i ++ ) {
            bond_type = (at->bond_type[i] & BOND_TYPE_MASK);
            if ( bond_type <= BOND_TYPE_TRIPLE ) {
                metal_bonds_chem_valence += inchi_max(BOND_TYPE_SINGLE, bond_type);
            } else {
                metal_bonds_chem_valence += BOND_TYPE_SINGLE;
                k ++; /* count alternating bonds */
            }
        }
    } else {
        for ( i = j = j2 = k = k2 = 0; i < at->valence; i ++ ) {
            bond_type = (at->bond_type[i] & BOND_TYPE_MASK);
            if ( pVArray[ (int)at->neighbor[i] ].cMetal ) {
                j ++;  /* number of bonds to metal atoms */
                if ( bond_type <= BOND_TYPE_TRIPLE ) {
                    metal_bonds_chem_valence += inchi_max(BOND_TYPE_SINGLE, bond_type);
                } else {
                    metal_bonds_chem_valence += BOND_TYPE_SINGLE;
                    k ++; /* count alternating bonds */
                }
            } else {
                j2 ++;
                if ( bond_type <= BOND_TYPE_TRIPLE ) {
                    not_metal_bonds_chem_valence += inchi_max(BOND_TYPE_SINGLE, bond_type);
                } else {
                    not_metal_bonds_chem_valence += BOND_TYPE_SINGLE;
                    k2 ++; /* count alternating bonds */
                }
            }
        }
    }
    bonds_chem_valence = metal_bonds_chem_valence + not_metal_bonds_chem_valence;
    if ( at->chem_bonds_valence > bonds_chem_valence ) {
        if ( at->chem_bonds_valence - bonds_chem_valence > 1 ) {
            at->chem_bonds_valence = bonds_chem_valence + 1;  /* should not happen */
        }
        alt_bonds_delta_valence = at->chem_bonds_valence - bonds_chem_valence;
    }
   
    pVA->cNumBondsToMetal = j;

    if ( nPeriodicNum == EL_NUMBER_H ) { 
        /* ignore bridging H; ??? later add ??? */
        return 0;
    }

    num_H     = at->num_H;
    num_bonds = at->valence;

    if ( !num_bonds && !num_H ) {
        return 0; /* do not know the answer: isolated atom */
    }
    /* at the beginning all bonds are single */
    min_chem_valence = num_bonds + num_H;
    cur_chem_valence = bonds_chem_valence + alt_bonds_delta_valence + num_H; /* includes double & alternating bond contribution */

    /* number of non-bonding electrons in case of all single bonds */
    num_non_bonding_electrons = (int)pVA->cNumValenceElectrons - min_chem_valence;
    /* Octet rule: charge = bonds_valence + NumValenceElectrons - 8 */
    charge = min_chem_valence + (int)pVA->cNumValenceElectrons - VALUE_OCTET; /* wrong */

    /* typical (ad hoc) minimal neutral valence */
    known_chem_valence = ( pVA->cNumValenceElectrons > VALUE_OCTET/2 )? 
                                VALUE_OCTET - pVA->cNumValenceElectrons :
                                pVA->cNumValenceElectrons;
    /* excess of typical valence over all-single-bonds valence */
    nOctetNeutralValenceExcess = known_chem_valence - min_chem_valence;
    /*  (NB=num.bonds, NV=neutral valence, NVX=neutral valence excess, LFVS=last found valence state, val.=valence)

       element NB  knownFst octet Last  octetNVX  firstNVX foundNVX chargeLFVS  LFVS
                   valence  val.  NV>=                                 
                                                                     
       -B       1    3        3    3        2         2   =     2        +2      
       >B       2    3        3    3        1         1   =     1        +1
       >B-      3    3        3    3        0         0   =     0         0
       >B<      4    3        3    3       -1        -1   <>   N/A       -1
                                                                      
       -C       1    4        4    4        3         3   =     3        N/A
       >C       2    4        4    4        2         2   =     2        +2  (-2)
       >C-      3    4        4    4        1         1   =     1        +1  (-1)
       >C<      4    4        4    4        0         0   =     0         0
       C(V)     5    4        4    N/A     -1        -1   <>   N/A       N/A
                                                                     
       -Si      1    4        4    4        3         3   =     3        N/A
       >Si      2    4        4    4        2         2   =     2        +2  (-2)
       >Si-     3    4        4    4        1         1   =     1        +1  (-1)
       >Si-     4    4        4    4        0         0   =     0         0
       Si(V)    5    4        4    N/A     -1        -1   <>   N/A       -1 
                                                                     
       -N       1    3        3    3        2         2   =     2        -2
       >N       2    3        3    3        1         1   =     1        -1
       >N-      3    3        3    3        0         0   =     0         0  (+2)
       >N<      4    3        3    5       -1        -1   <>    1        +1
       N(V)     5    3        3    5       -2        -2   <>    0         0
       N(VI)    6    3        3    N/A     -3        -3   <>   N/A       N/A
       N(VII)   7    3        3    N/A     -4        -4   <>   N/A       N/A

       -P       1    3        3    3        2         2   =     2        -2
       >P       2    3        3    3        1         1   =     1        -1
       >P-      3    3        3    3        0         0   =     0         0  (-2, +2)
       >P<      4    3        3    5       -1        -1   <>    1        +1  (-1)
       P(V)     5    3        3    5       -2        -2   <>    0         0  (-2)
       P(VI)    6    3        3    N/A     -3        -3   <>   N/A       -1
       P(VII)   7    3        3    N/A     -4        -4   <>   N/A       -2
       P(VIII)  8    3        3    N/A     -5        -5   <>   N/A       N/A

       -O       1    2        2    2        1         1   =     1        -1
       >O       2    2        2    2        0         0   =     0         0
       >O-      3    2        2    N/A     -1        -1   <>   N/A       +1
       >O<      4    2        2    N/A     -2        -2   <>   N/A       +2
       O(V)     5    2        2    N/A     -3        -3   <>   N/A       +1
       O(VI)    6    2        2    N/A     -4        -4   <>   N/A       N/A

       -S       1    2        2    2        1         1   =     1        -1
       >S       2    2        2    2        0         0   =     0         0       NPNP - prohibit
       >S-      3    2        2    4       -1        -1   <>    1        +1  (-1) PNPN
       >S<      4    2        2    4       -2        -2   <>    0         0  (+2)
       S(V)     5    2        2    6       -3        -3   <>    1        +1  (-1)
       S(VI)    6    2        2    6       -4        -4   <>    0         0 
       S(VII)   7    2        2    N/A     -5        -5   <>    0        -1 
       S(VIII)  8    2        2    N/A     -6        -6   <>   N/A       N/A

       -F       1    1        1    1        0         0   =     0         0
       >F       2    1        1    1       -1        -1   <>   N/A       +1
       >F-      3    1        1    1       -2        -2   <>   N/A       +2
       >F<      4    1        1    1       -3        -3   <>   N/A       N/A
       F(V)     5    1        1    1       -4        -4   <>   N/A       +2
       F(VI)    6    1        1    1       -5        -5   <>   N/A       N/A

       -Cl      1    1        1    1        0         0   =     0         0      NPNP - prohibit
       >Cl      2    1        1    3       -1        -1   <>    1        +1      PNPN - prohibit
       >Cl-     3    1        1    3       -2        -2   <>    0         0 (+2) NPNP
       >Cl<     4    1        1    5       -3        -3   <>    1        +1      PNPN
       Cl(V)    5    1        1    5       -4        -4   <>    0         0
       Cl(VI)   6    1        1    7       -5        -5   <>    1        +1 
       Cl(VII)  7    1        1    7       -6        -6   <>    0         0
       Cl(VIII) 8    1        1    N/A     -7        -7   <>   N/A       N/A
       

      NB                 = num_bonds+num_H

      knownFst valence   = nFirstNeutralValenceExcess + min_chem_valence
      octet val.         = nOctetNeutralValenceExcess + min_chem_valence
      Last NV>=          = nFoundNeutralValenceExcess + min_chem_valence

      octetNVX           = nOctetNeutralValenceExcess
      firstNVX           = nFirstNeutralValenceExcess
      foundNVX           = nFoundNeutralValenceExcess

      chargeLFVS         = ChargeVal[nLastFoundValenceState].nCharge

    */
    /* minimal known neutral atom valence; different for Sn(2/4), Tl(1/3), Pb(2/4): (known/typical ad hoc) */
    known_chem_valence = get_el_valence( nPeriodicNum, 0, 0 );

    if ( pSrm->bMetalAddFlower ) {
        /* bond orders of bonds to metal may be as they are (pSrm->nMetalInitBondOrder==1)
                                        or decreased by one (pSrm->nMetalInitBondOrder==0)
           nMetalInitBondOrder == nMetalMinBondOrder + nMetalInitEdgeFlow
        */
        cur_chem_valence_fixed = cur_chem_valence - pVA->cNumBondsToMetal * (1-pSrm->nMetalInitBondOrder);
        pVA->cInitOrigValenceToMetal = metal_bonds_chem_valence;
        pVA->cInitValenceToMetal = metal_bonds_chem_valence - pVA->cNumBondsToMetal * (1-pSrm->nMetalInitBondOrder);
        pVA->cInitFlowToMetal  = pVA->cInitValenceToMetal - pVA->cNumBondsToMetal * pSrm->nMetalMinBondOrder;
        if ( pVA->cMetal ) {
            pVA->cInitFreeValences += alt_bonds_delta_valence;
        }
        if ( pSrm->nMetalInitEdgeFlow < pSrm->nMetalInitBondOrder - pSrm->nMetalMinBondOrder ) {
            /* single bond has zero initial flow + 2 radicals at incident atoms */
            if ( pVA->cInitFlowToMetal <= pVA->cNumBondsToMetal ) {
                if ( pVA->cMetal ) {
                    pVA->cInitFreeValences += pVA->cInitFlowToMetal;
                }
                pVA->cInitFlowToMetal = 0;
            } else {
                if ( pVA->cMetal ) {
                    pVA->cInitFreeValences += pVA->cNumBondsToMetal * (1 - pSrm->nMetalInitEdgeFlow);
                }
                pVA->cInitFlowToMetal -= pVA->cNumBondsToMetal * (1 - pSrm->nMetalInitEdgeFlow);
            }
        }
             
    } else {
        /* treat metal atoms as ordinary non-metal atoms */
        cur_chem_valence_fixed = cur_chem_valence;
        pVA->cInitFlowToMetal  = metal_bonds_chem_valence - pVA->cNumBondsToMetal;
        pVA->cInitValenceToMetal = metal_bonds_chem_valence;
        pVA->cInitOrigValenceToMetal = metal_bonds_chem_valence;
    }


    if ( pVA->cMetal && pSrm->bMetalAddFlower ) {
        pVA->cnListIndex = cnListIndexMe + 1;
        /*
        pVA->cInitOrigValenceToMetal += alt_bonds_delta_valence;
        pVA->cInitValenceToMetal     += alt_bonds_delta_valence;
        pVA->cInitFreeValences = (pSrm->nMetalInitBondOrder + alt_bonds_delta_valence
                                 - (pSrm->nMetalMinBondOrder + pSrm->nMetalInitEdgeFlow)) * pVA->cNumBondsToMetal;
        */
        return 0;  /* metal */
    }

    if ( !known_chem_valence ) {
        /* a noble gas like He, Ne, ... */
        pVA->cInitFreeValences = at->chem_bonds_valence - at->valence;
        return TREAT_ATOM_AS_METAL; /* do not know anything about this atom; needs 2nd pass */
    }

    nFirstNeutralValenceExcess = known_chem_valence - min_chem_valence;
    
    nFoundNeutralValenceExcess    = NO_VALUE_INT;
    nFoundNeutralValenceOrdNumber = NO_VALUE_INT;
    nLastFoundValenceOrdNumber    = NO_VALUE_INT;
    nLastFoundValenceState        = NO_VALUE_INT;
    
    /* find the lowest known valence >= all-single-bonds valence */
    for ( cur_charge = MIN_ATOM_CHARGE, nNumStates = 0; cur_charge <= MAX_ATOM_CHARGE; cur_charge ++ ) {
        for ( i = 0; i < MAX_NUM_VALENCES; i ++ ) {
            known_chem_valence = get_el_valence( nPeriodicNum, cur_charge, i );
            if ( cur_chem_valence_fixed > known_chem_valence || !known_chem_valence ) {
                continue; /* known valence < all-single-bonds valence */
            }
            if ( BOND_TYPE_TRIPLE + BOND_TYPE_DOUBLE * (num_bonds - 1) + num_H < known_chem_valence ) {
                continue; /* not more than one triple and the rest - double bonds per atom */
            }
            /* keep all found */
            ChargeVal[nNumStates].nValence               = known_chem_valence;
            ChargeVal[nNumStates].nCharge                = cur_charge;
            ChargeVal[nNumStates].nValenceOrderingNumber = i;
            if ( !cur_charge && nFoundNeutralValenceExcess == NO_VALUE_INT ) {
                /* neutral state; compare to the lowest typical valence */
                nFoundNeutralValenceExcess    = known_chem_valence - min_chem_valence;
                nFoundNeutralValenceOrdNumber = i;
            }
            if ( min_chem_valence == known_chem_valence ) {
                if ( nLastFoundValenceState == NO_VALUE_INT ) {
                    /* accept the first found */
                    nLastFoundValenceState = nNumStates;
                } else
                if ( abs( ChargeVal[nLastFoundValenceState].nCharge )  >= abs( cur_charge ) ) {
                    /* accept smaller abs(charge); if abs(charges) are same, accept (+) */
                    nLastFoundValenceState = nNumStates;
                }
            }
            nNumStates ++;
        }
    }
    /***********************************************************************************/
    /* select only appropriate charge & valence so that a suitable ChargeStruct exists */
    /***********************************************************************************/

    nNumSelectedStates = clean_charge_val( ChargeVal, nNumStates, atom, pVArray, iat, pVA->cMetal, bMobileH, endpoint );

    if ( !nNumSelectedStates ) {
        return TREAT_ATOM_AS_METAL; /* nothing to do */
    }
    /***********************************************************************************/
    /*       Find an appropriate ChargeStruct index for the ChargeVal found            */
    /***********************************************************************************/
    cn_bits = 0;
    memset( cn_bits_array, 0, sizeof(cn_bits_array) );
    /***** set bits identifying a suitable ChargeStruct ******/ 
    for ( i = len_cn_bits_array = 0; i < nNumSelectedStates && len_cn_bits_array < 4; i ++ ) {
        switch( ChargeVal[i].nCharge ) {
        case -1:
            cn_bits_array[len_cn_bits_array] |= cn_bits_M; /* Minus 1 */
            break;
        case 0:
            cn_bits_array[len_cn_bits_array] |= cn_bits_N; /* Neutral */
            break;
        case 1:
            cn_bits_array[len_cn_bits_array] |= cn_bits_P; /* Plus 1 */
            break;
        default:
            return RI_ERR_PROGR; /* program error */
        }
        if ( i+1 < nNumSelectedStates &&
             ChargeVal[i].nValence == ChargeVal[i+1].nValence &&
             ChargeVal[i].nCharge &&
             ChargeVal[i].nCharge == -ChargeVal[i+1].nCharge ) {
            ; /* add opposite charge to the same element of cn_bits_array[] */
        } else {
            len_cn_bits_array ++;
        }
    }
    if ( !len_cn_bits_array || len_cn_bits_array > 4 ) {
        return RI_ERR_PROGR; /* program error */
    }
    /* accommodate added 4-state ChargeStruct: +/- cannot be in case of 4 states */
    if ( len_cn_bits_array + 1 == nNumSelectedStates && nNumSelectedStates == 4 ) {
        len_cn_bits_array --;
        nNumSelectedStates --;
        cn_bits_array[len_cn_bits_array] = 0;
    }
    /* fix for terminal hydrogenless -C as in isocyano or CO: there is no just cnE_[] ChargeStruct */
    
    if ( len_cn_bits_array == 1 &&
         cn_bits_array[0] == (cn_bits_P | cn_bits_M) &&
         ChargeVal[0].nValence + 1 > BOND_TYPE_TRIPLE + BOND_TYPE_DOUBLE * (num_bonds - 1) + num_H ) {
        cn_bits_array[len_cn_bits_array ++] = cn_bits_N;
        ChargeVal[nNumSelectedStates].nValence = ChargeVal[nNumSelectedStates-1].nValence;
        ChargeVal[nNumSelectedStates].nCharge  = 0;
        ChargeVal[nNumSelectedStates].nValenceOrderingNumber = 0;
    }
    
make_cn_bits:
    cn_bits = MAKE_CN_BITS(cn_bits_array[0], cn_bits_array[1], cn_bits_array[2], cn_bits_array[3]);
    /*********** find ChargeStructure **************/
    for ( i = 0, j = -1; i < cnListNumEl; i ++ ) {
        if ( cnList[i].bits == cn_bits ) {
            j = i;
            break; /* found */
        }
    }
    if ( j < 0 ) {
        /* ChargeStructure was not found */
        if ( 1 < len_cn_bits_array && len_cn_bits_array + 1 == nNumSelectedStates ) {
            /* a pair of opposite charges was combined */
            len_cn_bits_array --;
            cn_bits_array[len_cn_bits_array] = 0;
            goto make_cn_bits;
        } else
        if ( nNumSelectedStates == 4 ) {
            /* reduce number of states */
            len_cn_bits_array --;
            cn_bits_array[len_cn_bits_array] = 0;
            nNumSelectedStates --;
            goto make_cn_bits;
        }
        return RI_ERR_PROGR; /* charge structure not found */
    }
    /********** ChargeStructure has been found **********/
    pVA->cnListIndex = j+1; /* charge structure index + 1 */
    pVA->cInitCharge = cnList[j].nInitialCharge;
    /********** Calculate "Free Valence" ****************/
#if ( ALLOW_METAL_BOND_ZERO == 1 )

#if ( INIT_METAL_BOND_ZERO == 1 )
    if ( pVA->cMetal ) {
        j = 0;
    } else {
        j = ChargeVal[0].nValence - cur_chem_valence_fixed;
    }
#else
    j = ChargeVal[0].nValence - cur_chem_valence_fixed;
#endif

#else
    j = ChargeVal[0].nValence - cur_chem_valence_fixed;
#endif
    if ( j < 0 ) {
        return RI_ERR_PROGR; /* program error */
    }
    pVA->cInitFreeValences = j; /* number of initial unsatisfied valences; should be combined with */
                                /* (cap - flow) of vertex=0 in the charge structure[pVA->cnListIndex-1] */ 
    return 1; /* success */

#undef MIN_ATOM_CHARGE 
#undef MAX_ATOM_CHARGE 
#undef NEUTRAL_STATE   
#undef NUM_ATOM_CHARGES
#undef MAX_NUM_VALENCES

}


#ifdef NEVER
/******************************************************************************************************/
int get_bonds_valences( int nPeriodicNum, int bonds_valence, int num_H, VAL_AT *pVA )
{
    int i, j, charge, chem_valence, known_chem_valence;
#define MAX_NUM_VALENCES 5  /* defined in util.c */

    memset( pVA, 0, sizeof( pVA[0] ) );
    
    if ( !bonds_valence && !num_H )
        return 0; /* do not know the answer */
    
    chem_valence = bonds_valence + num_H;
    for ( charge = VAL_MIN_CHARGE; charge <= VAL_MAX_CHARGE; charge ++ ) {
        for ( i = 0, j = 0; i < MAX_NUM_VALENCES, j < VAL_NUMBER; i ++ ) {
            if ( chem_valence <= (known_chem_valence = get_el_valence( nPeriodicNum, charge, i ) ) ) {
                if ( !charge ) {
                    pVA->cValence[j][VAL_NEUTR_ORDER] = i+1;
                }
                pVA->cValence[j++][charge+VAL_BASE] = known_chem_valence - num_H;
            }
        }
    }
    pVA->cDoNotAddH = do_not_add_H( nPeriodicNum );
    pVA->cMetal     = is_el_a_metal( nPeriodicNum );
    return pVA->cValence[0][VAL_BASE];  /* 0 means do not know the answer */
#undef MAX_NUM_VALENCES
}
#endif
/*********** calculate s or p-element type ************/
int get_sp_element_type( int nPeriodicNumber, int *nRow )
/* 
                                                num  el
                                                el   neg
   1 => H                          ATYPE_H   1    1  21
   2 => Li, Na, K,  Rb, Cs, Fr     ATYPE_Na  2    1  10 09 08 08 07 
   3 => Be, Mg, Ca, Sr, Ba, Ra     ATYPE_Mg  3    2  15 12 10 10 09
   4 => B,  Al, Ga, In, Tl         ATYPE_B   4    3  20 15 18 17 18
   5 => C,  Si, Ge, Sn, Pb         ATYPE_C   5    4  25 18 18 18 18 
   6 => N,  P,  As, Sb, Bi         ATYPE_N   6    5  30 21 20 19 19
   7 => O,  S,  Se, Te, Po         ATYPE_O   7    6  35 25 24 21 20
   8 => F,  Cl, Br, I,  At         ATYPE_Cl  8    7  40 30 28 25 22

number of valence electrons = (type>1)? type-1: type

  */
{
    int row = 0, type = 0;
    if ( nPeriodicNumber == 1 ) {
        type = 1; /* H: 1 */
        row  = 0;
    } else
    if ( nPeriodicNumber == 2 ) {
        type = 0; row = 0;
    } else
    if ( nPeriodicNumber <= 10 ) {
        /* Li: 2, Be: 3, B: 4, C: 5, N: 6, O: 7, F: 8, Ne: 9; later subtract 1 */
        type = nPeriodicNumber-1; row = 1;
    } else
    if ( nPeriodicNumber <= 18 ) {
        type = nPeriodicNumber - 9; row = 2;
    } else
    if ( nPeriodicNumber <= 20 ) {
        type = nPeriodicNumber - 17; row = 3;
    } else
    if ( nPeriodicNumber <= 30 ) {
        type = 0; row = 3;
    } else
    if ( nPeriodicNumber <= 36 ) {
        type = nPeriodicNumber - 27; row = 3;
    } else
    if ( nPeriodicNumber <= 38 ) {
        type = nPeriodicNumber - 35; row = 4;
    } else
    if ( nPeriodicNumber <= 48 ) {
        type = 0; row = 4;
    } else
    if ( nPeriodicNumber <= 54 ) {
        type = nPeriodicNumber - 45; row = 4;
    } else
    if ( nPeriodicNumber <= 56 ) {
        type = nPeriodicNumber - 53; row = 5;
    } else
    if ( nPeriodicNumber <= 80 ) {
        type = 0; row = 5;
    } else
    if ( nPeriodicNumber <= 86 ) {
        type = nPeriodicNumber - 77; row = 5;
    } else
    if ( nPeriodicNumber <= 88 ) {
        type = nPeriodicNumber - 85; row = 6;
    } else {
        type = 0; row = 6;
    }
    *nRow = row;
    return type==9? 0 : type;
}
/******************************************************************************************************/
int ReallocTCGroups( ALL_TC_GROUPS *pTCGroups, int nAdd )
{
    TC_GROUP *pTCGroup = (TC_GROUP *) inchi_malloc( sizeof(pTCGroup[0])*(pTCGroups->max_tc_groups + nAdd) );
    if ( pTCGroup ) {
        if ( pTCGroups->num_tc_groups ) {
            memcpy( pTCGroup, pTCGroups->pTCG, sizeof(pTCGroup[0])*pTCGroups->num_tc_groups );
        }
        memset( pTCGroup + pTCGroups->max_tc_groups, 0, sizeof(pTCGroup[0])*nAdd );
        if ( pTCGroups->pTCG ) {
            inchi_free( pTCGroups->pTCG );
        }
        pTCGroups->pTCG = pTCGroup;
        pTCGroups->max_tc_groups += nAdd;
        return 0;
    }
    return RI_ERR_ALLOC;
}
/******************************************************************************************************/
int RegisterTCGroup( ALL_TC_GROUPS *pTCGroups, int nGroupType, int nGroupOrdNum,
                     int nVertexCap, int nVertexFlow, int nEdgeCap, int nEdgeFlow, int nNumEdges)
{
    int i, ret = 0;
    /* search */
    for ( i = 0; i < pTCGroups->num_tc_groups; i ++ ) {
        if ( pTCGroups->pTCG[i].type    == nGroupType &&
             pTCGroups->pTCG[i].ord_num == nGroupOrdNum ) {
            break;
        }
    }
    if ( i == pTCGroups->num_tc_groups ) {
        /* add one more group */
        if ( pTCGroups->num_tc_groups == pTCGroups->max_tc_groups ) {
            ret = ReallocTCGroups( pTCGroups, INC_NUM_TCGROUPS );
            if ( ret ) {
                goto exit_function;
            }
        }
        ret = i+1; /* added new group */
        pTCGroups->num_tc_groups ++;
        pTCGroups->pTCG[i].type    = nGroupType;
        pTCGroups->pTCG[i].ord_num = nGroupOrdNum;
    }
    pTCGroups->pTCG[i].num_edges  += nNumEdges;

    pTCGroups->pTCG[i].st_cap     += nVertexCap;
    pTCGroups->pTCG[i].st_flow    += nVertexFlow;

    pTCGroups->pTCG[i].edges_cap  += nEdgeCap;
    pTCGroups->pTCG[i].edges_flow += nEdgeFlow;

exit_function:
    return ret;         
}
/******************************************************************************************************/
int nTautEndpointEdgeCap( inp_ATOM *at, VAL_AT *pVA, int i )
{
    /* There are 3 sources of cap-flow = number of unsatisfied valences:
       -----------------------------------------------------------------
       1. pVA[i].cInitFreeValences 
       2. pCN[0].v.cap - pCN[0].v.flow
       3. st[i].chem_bonds_valence - SUM(SINGLE, DOUBLE, TRIPLE bond orders)
          Reasons: (a) This sum will not include 'ALTERN' bonds
                   (b) until now at[i].chem_bonds_valence was used as a
                       number of satisfied valences. In case of adjacent
                       stereobonds marked as BOND_TYPE_ALTERN the value of
                       at[i].chem_bonds_valence may be = at[i].valence+1.
       4. Since tautomerism is defined for a neutral atom, do not add
          initial flows from the atom to the ChargeStruct
          CORRECTION: tautomeric endpoints do not have ChargeStruct.
                       
     */
    int j, k, nEdgeCap, bonds_valence, stereo_bond_excess_valence;
    MY_CONST C_NODE *pCN = pVA[i].cnListIndex>0? cnList[pVA[i].cnListIndex-1].pCN:NULL;

    /* 1: free valences to reach the minimum known atom valence */
    nEdgeCap = pVA[i].cInitFreeValences;
    /* 2: atom free valence in the ChargeStruct */
    if ( pCN ) {
        nEdgeCap += pCN[0].v.cap - pCN[0].v.flow; /* normally should not happen */
    }
    /* 3: atom free valence due to known from stereochemistry stereogenic bond types */
    /*
    for ( j = 0, bonds_valence = 0; j < at[i].valence; j ++ ) {
        if ( at[i].bond_type[j] <= BOND_TYPE_TRIPLE ) {
            bonds_valence += at[i].bond_type[j];
        }
    }
    */
    /* bonds > SINGLE are assumed fixed stereobonds; fixed bond cannot increase t-group edge flow */
    for ( stereo_bond_excess_valence=0, j = 0; j < MAX_NUM_STEREO_BONDS && at[i].sb_parity[j]; j ++ ) {
        k = at[i].sb_ord[j];
        if ( at[i].bond_type[k] < BOND_TYPE_TRIPLE ) {
            stereo_bond_excess_valence += at[i].bond_type[k] - BOND_TYPE_SINGLE;
        }
    }
    /*
    bonds_valence = (at[i].chem_bonds_valence - bonds_valence) + (bonds_valence -at[i].valence - stereo_bond_excess_valence);
    */
    bonds_valence = (at[i].chem_bonds_valence - at[i].valence) - stereo_bond_excess_valence;


    /*---- add 1, 2, 3 ----*/
    if ( bonds_valence >= 0 ) {
        nEdgeCap += bonds_valence;
    } else {
        nEdgeCap = RI_ERR_PROGR;
    }
    return nEdgeCap;
}
/******************************************************************************************************/
/* If Metal flowers are allowed ( pSrm->bMetalAddFlower != 0), then:                                  */
/*                                                                                                    */
/*   bond to a metal atom                min_bond_order[i] = pSrm->nMetalMinBondOrder                 */
/*   taut endpoint - metal               min_bond_order[i] = pSrm->nMetal2EndpointMinBondOrder        */
/*   single bond to metal atom:      initial_bond_order[i] = pSrm->nMetalInitBondOrder                */
/*   n-order bond to metal atom      initial_bond_order[i] = pSrm->nMetalInitBondOrder + n-1          */
/*                                 = bond_order[i]-BOND_TYPE_SINGLE+pSrm->nMetalInitBondOrder         */
/*   single t-endpoint--atom bond    initial_bond_order[i] = pSrm->nMetal2EndpointInitBondOrder       */
/*   n-order t-endpoint--metal bond  initial_bond_order[i] = pSrm->nMetal2EndpointInitBondOrder+n-1   */
/*                                 = bond_order[i]-BOND_TYPE_SINGLE+pSrm->nMetal2EndpointInitBondOrder*/
/*                                                                                                    */
/* Exceptions from simple atom-metal conditions:                                                      */
/*   1. Atom is a tautomeric endpoint: use pSrm->nMetal2Endpoint* instead of pSrm->nMetal*            */
/*   2. Atom is sp3-stereogenic and pSrm->bFixStereoBonds != 0: use atom-atom rules                   */
/*   3. Atom has a sp2-stereo   and pSrm->bFixStereoBonds != 0: use atom-atom rules                   */
/*                                                                                                    */
/* Atom-atom rules (applies to all atoms if pSrm->bMetalAddFlower=0)                                  */
/*                                                                                                    */
/*   min_bond_order[i]      = BOND_TYPE_SINGLE  (BOND_TYPE_SINGLE = 1)                                */
/*   initial_bond_order[i]  = bond_type[i]                                                            */
/*                                                                                                    */
/* General rules:                                                                                     */
/*   initial_bond_flow[i]   = initial_bond_order[i]-min_bond_order[i]                                 */
/*   atom[k] initial_st_cap = at[k].chem_bonds_valence - SUM{i; initial_bond_order[i]}                */
/*   bond_cap[i]            =  BOND_TYPE_TRIPLE - min_bond_order[i]                                   */
/*     (reason: quadruple and higher order bonds are not allowed)                                     */
/* Exception: in case of metal-atom bond, if pSrm->nMetal2EndpointInitEdgeFlow = 0 AND                */
/* pSrm->nMetalInitBondOrder - pSrm->nMetalMinBondOrder = 1 then                                      */
/*   reduce bond to metal order by 1 and increase st_cap of both neighbors by 1:                      */
/*   initial_bond_flow[i] --; metal_initial_st_cap += num_bonds;                                      */
/* ==== Note: ONLY the INCREASE is already included in pVA->cInitFreeValences of both atoms           */
/*                                                                                                    */
/* Notes: initial_st_cap does not include:                                                            */
/*   1. atom[k] additional st_cap from ChargeStruct pCN[0].v.cap                                      */
/*   2. pVA[k].cInitFreeValences due to a difference between the smallest known valence and st_cap    */
/*                                                                                                    */
/*  here k=atom at[k] index,                                                                          */
/*       i=bond index; i = 0..at[k].valence;                                                          */
/*       SUM{i; M[i]} is a sum of M[i] over all i                                                     */
/*       bond_order[i] = at[k].bond_type[i] >= BOND_TYPE_SINGLE - input bond order                    */
/******************************************************************************************************/

/***************** new *************************************************************************************/

int BondFlowMaxcapMinorder( inp_ATOM *atom, VAL_AT *pVA, ICHICONST SRM *pSrm, int iat, int ineigh,
                            int *pnMaxcap, int *pnMinorder, int *pbNeedsFlower )
{
    int nFlow, nMaxcap, nMinorder, nInitorder, bNeedsFlower = 0;
    inp_ATOM *at        = atom + iat;
    int       neigh     = at->neighbor[ineigh];
    int       bond_type = at->bond_type[ineigh] & BOND_TYPE_MASK;
    int       nMetal    = (0 != pVA[iat].cMetal) + (0 != pVA[neigh].cMetal);
    int       nEndpoint = (0 != at->endpoint) + (0 != atom[neigh].endpoint);
    int       nStereo   = (at->p_parity || at->sb_parity[0]) + (atom[neigh].p_parity || atom[neigh].sb_parity[0]);

    if ( bond_type > BOND_TYPE_TRIPLE ) {
        bond_type = BOND_TYPE_SINGLE;
    }
    /* M=metal, A=non-metal atom, e=endpoint */
    if ( (nStereo && pSrm->bFixStereoBonds) || !nMetal || !pSrm->bMetalAddFlower ) {
        /* atom-atom rules, no metal atoms involved (1: A-A, A-Ae, Ae-Ae) */
        nMinorder  = BOND_TYPE_SINGLE;
        nInitorder = bond_type;
        nFlow      = nInitorder - nMinorder;
    } else
    if ( nMetal && !nEndpoint ) { 
        /* M-a, M-M */
        /* atom - metal or metal-metal, none of them is an endpoint (2: M-M, M-A) */
        nMinorder  = pSrm->nMetalMinBondOrder;
        nInitorder = pSrm->nMetalInitBondOrder + bond_type - BOND_TYPE_SINGLE;
        nFlow      = nInitorder - nMinorder;
        if ( !pSrm->nMetalInitEdgeFlow &&
             pSrm->nMetalInitBondOrder > pSrm->nMetalMinBondOrder &&
             nFlow > 0 ) {
            /* reduce initial flow by 1 and increase st_cap on metal by 1 */
            nFlow --;
        }
        bNeedsFlower = (0 != pVA[iat].cMetal);
    } else
    if ( (pVA[iat].cMetal   && !at->endpoint         && !pVA[neigh].cMetal && atom[neigh].endpoint) ||
         (pVA[neigh].cMetal && !atom[neigh].endpoint && !pVA[iat].cMetal   && at->endpoint) ) {
        /* M-ae */
        /* metal connected to a non-metal endpoint (3: M-Ae) */
        nMinorder  = pSrm->nMetal2EndpointMinBondOrder;
        nInitorder = pSrm->nMetal2EndpointInitBondOrder + bond_type - BOND_TYPE_SINGLE;
        nFlow      = nInitorder - nMinorder;
        if ( !pSrm->nMetal2EndpointInitEdgeFlow &&
             pSrm->nMetal2EndpointInitBondOrder > pSrm->nMetal2EndpointMinBondOrder &&
             nFlow > 0 ) {
            /* reduce initial flow by 1 and increase st_cap on metal by 1 */
            nFlow --;
        }
        bNeedsFlower = (0 != pVA[iat].cMetal);
    } else {
        /* endpoint is metal => no flower (4: M-Me, Me-Me, Me-A, Me-Ae) */
        nMinorder  = pSrm->nMetal2EndpointMinBondOrder;
        nInitorder = pSrm->nMetal2EndpointInitBondOrder + bond_type - BOND_TYPE_SINGLE;
        nFlow      = nInitorder - nMinorder;
        if ( !pSrm->nMetal2EndpointInitEdgeFlow &&
             pSrm->nMetal2EndpointInitBondOrder > pSrm->nMetal2EndpointMinBondOrder &&
             nFlow > 0 ) {
            /* reduce initial flow by 1 and increase st_cap on metal by 1 */
            nFlow --;
        }
        bNeedsFlower = (pVA[iat].cMetal && !at->endpoint);
    }
    nMaxcap = BOND_TYPE_TRIPLE - nMinorder;
    if ( pnMaxcap ) {
        *pnMaxcap = nMaxcap;
    }
    if ( pnMinorder ) {
        *pnMinorder = nMinorder;
    }
    if ( pbNeedsFlower ) {
        *pbNeedsFlower = bNeedsFlower;
    }
    return nFlow;
}
/*********** new *******************************************************************************************/
int AtomStcapStflow( inp_ATOM *atom, VAL_AT *pVA, ICHICONST SRM *pSrm, int iat, int *pnStcap, int *pnStflow,
                     EdgeFlow *pnMGroupEdgeCap, EdgeFlow *pnMGroupEdgeFlow )
{
    int ineigh, bFlower;
    int nStflow=0, nMaxBondCap, nMinBondOrder, bNeedsFlower = 0;
    int valence = atom[iat].valence;
    int nStcap  = atom[iat].chem_bonds_valence;
    int nMGroupEdgeCap = 0, nMGroupEdgeFlow = 0, nFlow;

    if ( pSrm->bMetalAddFlower ) {
        nStcap  -= pVA[iat].cInitOrigValenceToMetal - pVA[iat].cInitValenceToMetal;
    }

    for ( ineigh = 0; ineigh < valence; ineigh ++ ) {
        nFlow     = BondFlowMaxcapMinorder( atom, pVA, pSrm, iat, ineigh, &nMaxBondCap, &nMinBondOrder, &bFlower );
        nStflow  += nFlow;
        nStcap   -= nMinBondOrder;
        if ( bFlower ) {
            bNeedsFlower ++;
            nMGroupEdgeFlow += nFlow;
            nMGroupEdgeCap  += BOND_TYPE_TRIPLE - nMinBondOrder + pSrm->nMetalMaxCharge_D; 
        }
    }
    if ( pnStcap ) {
        *pnStcap = bNeedsFlower? nStflow : nStcap; /* initially, metal atoms are not radicals */
    }
    if ( pnStflow ) {
        *pnStflow = nStflow;
    }
    if ( pnMGroupEdgeFlow ) {
        *pnMGroupEdgeFlow = nMGroupEdgeCap - nMGroupEdgeFlow;
    }
    if ( pnMGroupEdgeCap ) {
        *pnMGroupEdgeCap = nMGroupEdgeCap;
    }
    return bNeedsFlower; /* number of variable bonds to metal */
}
/**************************************************************************************
int nCountBnsSizes( inp_ATOM *at, int num_at, int nAddEdges2eachAtom, int nAddVertices,
                    T_GROUP_INFO *ti, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups )

fills out totals:

  pTCGroups->num_atoms        = number of atoms
  pTCGroups->num_bonds        = number of bonds between atoms
  pTCGroups->num_tgroups      = number of tautomeric groups
  pTCGroups->num_tgroup_edges = number of edges to tautomeric groups
  pTCGroups->tgroup_charge    = total charge on thautomeric atoms (negative)
  pTCGroups->num_tc_groups    = total number of all groups
  pTCGroups->nVertices        = total number of vertices excluding groups interconnections
  pTCGroups->nEdges           = total number of edges excluding groups interconnections

creates entries for the groups and adds to each group:
    
  TC_GROUP::type       =  BNS_VERT_TYPE_TGROUP, BNS_VT_C_POS, BNS_VT_C_NEG, BNS_VT_C_POS_C, BNS_VT_C_NEG_C 
  TC_GROUP::ord_num    =  ordering number within the type, e.g. t-group number
  TC_GROUP::st_cap     =  all from the atoms in ChargeStruct or tautomeric group info.
  TC_GROUP::st_flow    =  all from the atoms in ChargeStruct (0 for t-groups).
  TC_GROUP::num_edges  =  number of edges to the atoms or ChargeStruct vertices.
  TC_GROUP::edges_cap  =  sum of all incoming edge caps; see also nTautEndpointEdgeCap(..).
  TC_GROUP::edges_flow =  sum of all incoming edge flows; 0 for t-groups.

  TC_GROUP::nVertexNumber - NO FILLED WITH ANYTHING

  Note: the nDelta = st_cap - st_flow needs to be preserved when adding more vertices
  
Return value: =0 => success
              <0 => error
 **************************************************************************************/
int nCountBnsSizes( inp_ATOM *at, int num_at, int nAddEdges2eachAtom, int nAddVertices,
                    T_GROUP_INFO *ti, VAL_AT *pVA, ICHICONST SRM *pSrm, ALL_TC_GROUPS *pTCGroups )
{
    int i, j, n, k, ret = 0, nBonds, nOtherEdges, nVertices, bMetalAtoms, bNeedsFlower;
    int nTgroupEdges, nTgroupEdgesFromTg, nTotNegChargInTgroups, cap, flow;
    MY_CONST C_NODE *pCN = NULL;
    nVertices = nBonds = nOtherEdges = nTgroupEdges = nTgroupEdgesFromTg = nTotNegChargInTgroups = 0;
    
    /* count metal atoms and electrons */
    for ( i = 0; i < num_at; i ++ ) {
        pTCGroups->num_metal_atoms += (pVA[i].cMetal != 0);
        pTCGroups->num_metal_bonds += pVA[i].cNumBondsToMetal;
        pTCGroups->total_electrons += at[i].el_number;
        pTCGroups->total_electrons_metals += pVA[i].cMetal? at[i].el_number : 0;
    }
    pTCGroups->total_electrons -= pTCGroups->total_charge;
    pTCGroups->num_metal_bonds /= 2;

    /* register tautomeric groups */
    for ( i = 0; i < ti->num_t_groups; i ++ ) {
        ret = RegisterTCGroup( pTCGroups, BNS_VERT_TYPE_TGROUP, ti->t_group[i].nGroupNumber,
                               ti->t_group[i].num[0] /* st_cap */, 0 /* st_flow */,
                               0 /* edge cap */, 0 /* edge flow */,  ti->t_group[i].nNumEndpoints /* num Edges */ );
        if ( ret < 0 ) {
            goto exit_function;
        }
        /* edges to tautomeric groups */
        nOtherEdges           += ti->t_group[i].nNumEndpoints;
        nTgroupEdgesFromTg    += ti->t_group[i].nNumEndpoints;
        /* total negative charge in t-groups */
        nTotNegChargInTgroups += ti->t_group[i].num[1];
        if ( ret > 0 ) {
            /* should always happen since this is the first time this t-group is added */
            j = ret-1;
            pTCGroups->pTCG[j].tg_num_H     = ti->t_group[i].num[0] - ti->t_group[i].num[1];
            pTCGroups->pTCG[j].tg_num_Minus = ti->t_group[i].num[1];
        }
    }

    bMetalAtoms = 0;

repeat_for_metals:    

    /* count vertices and register ChargeValence groups */
    /* for now an atom may belong either to a t-group or to a ChargeValence group, but not to both */
    for ( i = 0; i < num_at; i ++ ) {
        /* number of bonds */
        nBonds += at[i].valence;
        /* Process ChargeStruct vertices and edges */
        if ( pVA[i].cnListIndex ) {
            /* count vertices & edges in the ChargeValence Substructure attached to an atom */
            /* Important: unlike inp_ATOM, each edge e appears in pCN[*].e[*] only ONE time */
            int     len     = cnList[j = pVA[i].cnListIndex-1].len;
            int     bits    = cnList[j].bits;
            int     type, neigh_type, metal_group_number;
            pCN     = cnList[j].pCN;
            
            /* first process all non-metals, after that -- all metals */
            if ( (bits != cn_bits_Me) != !bMetalAtoms ) {
                continue;
            }
            metal_group_number = 0;
            for ( j = 0; j < len; j ++) {
                type = pCN[j].v.type; /* ChargeStruct vertex type: atom is the first, c-groups are last */

                /* process all pCN[j] neighbors */
                for ( k = 0; k < MAX_CN_VAL && (n = pCN[j].e[k].neigh); k ++ ) {
                    nOtherEdges ++;  /* edges inside ChargeStruct */
                    n --; /* neighbor vertex position inside cnList[j].pCN */
                    neigh_type = pCN[n].v.type; /* type of the neighboring atom */

                    if ( IS_BNS_VT_C_GR(neigh_type) ) {
                        /* register this edge to a CN-group vertex */
                        cap  = !bMetalAtoms? pCN[j].e[k].cap  : pCN[j].e[k].cap?  pSrm->nMetalMaxCharge_D : 0;
                        flow = !bMetalAtoms? pCN[j].e[k].flow : pCN[j].e[k].flow? pSrm->nMetalMaxCharge_D : 0;

                        ret = RegisterTCGroup( pTCGroups, neigh_type, 0 /* ord_num*/,
                                               0 /* st_cap */, 0 /* st_flow */,
                                               cap /* edge cap*/, flow /* edge flow */, 1 /* nNumEdges*/);
                        if ( ret < 0 ) {
                            goto exit_function;
                        }
                        if ( ret > 0 ) {
                            /* the group has just been created; add one more edge to (+/-) or supergroup */
                            ret = RegisterTCGroup( pTCGroups, neigh_type, 0 /* ord_num*/,
                                                   0 /* st_cap */, 0 /* st_flow */,
                                                   0 /* edge cap*/, 0/* edge flow*/, 1 /* nNumEdges*/);
                            if ( ret < 0 ) {
                                goto exit_function;
                            }
                            nOtherEdges ++;
                        }
                    }

                    if ( IS_BNS_VT_C_GR(type) ) {
                        /* register this edge to a CN-group vertex; normally this does not happen */
                        cap  = !bMetalAtoms? pCN[j].e[k].cap  : pCN[j].e[k].cap?  pSrm->nMetalMaxCharge_D : 0;
                        flow = !bMetalAtoms? pCN[j].e[k].flow : pCN[j].e[k].flow? pSrm->nMetalMaxCharge_D : 0;
                        ret = RegisterTCGroup( pTCGroups, type, 0 /* ord_num*/,
                                               0 /* st_cap */, 0 /* st_flow */,
                                               cap /* edge cap*/, flow /* edge flow */, 1 /* nNumEdges*/);
                        if ( ret < 0 ) {
                            goto exit_function;
                        }
                        if ( ret > 0 ) {
                            /* the group has just been created; add one more edge to (+/-) or supergroup */
                            ret = RegisterTCGroup( pTCGroups, type, 0 /* ord_num*/,
                                                   0 /* st_cap */, 0 /* st_flow */,
                                                   0 /* edge cap*/, 0/* edge flow*/, 1 /* nNumEdges*/);
                            if ( ret < 0 ) {
                                goto exit_function;
                            }
                            nOtherEdges ++;
                        }
                    }
                } /* end of the current vertex pCN[j] neighbors */

                /* process  pCN[j] vertex */

                if ( type & BNS_VERT_TYPE_ATOM ) {
                    continue;  /* do not count regular atoms here */
                }
                if ( IS_BNS_VT_CHRG_STRUCT(type) ) {
                    nVertices ++;
                    continue;
                }

                if ( pSrm->bMetalAddFlower && IS_BNS_VT_M_GR( type ) ) {
                    /* special treatment: flow and cap are known as well as structure */
                    /* initial bond valence to metal is either 0 or 1 */
                    EdgeFlow nEdgeFlow, nEdgeCap;
                    bNeedsFlower = AtomStcapStflow( at, pVA, pSrm, i, NULL /*pnStcap*/, NULL /*pnStflow*/,
                                                    &nEdgeCap, &nEdgeFlow );
                    if ( !bNeedsFlower ) {
                        ret = RI_ERR_PROGR;
                        goto exit_function;
                    }
                    /*
                    GetAtomToMCGroupInitEdgeCapFlow( &nEdgeCap, &nEdgeFlow, pSrm, at,  pVA, i );
                    GetAtomToMCGroupInitEdgeCapFlow( &nEdgeCap, &nEdgeFlow, pSrm );
                    */
                    /* the 1st is the flower base */
                    /* atom - G0 edge and G0 vertex */
                    ret = RegisterTCGroup( pTCGroups, type, 0 /* ord_num*/,
                                           /*pVA[i].cInitFreeValences*/ 0 /* st_cap */, 0 /* st_flow */,
                                           (int)nEdgeCap, (int)nEdgeFlow, 1 /* nNumEdges*/);
                    if ( ret < 0 ) {
                        goto exit_function;
                    }
                    /* count edge atom-G0 */
                    nOtherEdges ++;
                    if ( ret > 0 ) {
                        /* first time registration: add G0-G1 and G0-G2 edges to G0 */
                        ret = RegisterTCGroup( pTCGroups, type, 0 /* ord_num*/,
                                               0 /* st_cap */, 0 /* st_flow */,
                                               0,/* edge cap*/ 0 /*edge flow*/, 2 /* nNumEdges*/);
                        
                        if ( ret < 0 ) {
                            goto exit_function;
                        }
                        /* first time registration: add G1; it has 3 edges */
                        ret = RegisterTCGroup( pTCGroups, type, 1 /* ord_num*/,
                                               0 /* st_cap */, 0 /* st_flow */,
                                               0,/* edge cap*/ 0 /*edge flow*/, 3 /* nNumEdges*/);
                        
                        if ( ret <= 0 ) {
                            ret = !ret? RI_ERR_PROGR : ret;
                            goto exit_function;
                        }
                        /* first time registration: add G2; it has 3 edges */
                        ret = RegisterTCGroup( pTCGroups, type, 2 /* ord_num*/,
                                               0 /* st_cap */, 0 /* st_flow */,
                                               0,/* edge cap*/ 0 /*edge flow*/, 3 /* nNumEdges*/);
                        
                        if ( ret <= 0 ) {
                            ret = !ret? RI_ERR_PROGR : ret;
                            goto exit_function;
                        }
                        /* first time registration: add G3; it has 2 edges */
                        ret = RegisterTCGroup( pTCGroups, type, 3 /* ord_num*/,
                                               0 /* st_cap */, 0 /* st_flow */,
                                               0,/* edge cap*/ 0 /*edge flow*/, 2 /* nNumEdges*/);
                        
                        if ( ret <= 0 ) {
                            ret = !ret? RI_ERR_PROGR : ret;
                            goto exit_function;
                        }
                        /* count added metal flower vertices: G0, G1, G2, G3 */
                        nVertices += 4;
                        /* count added metal flower edges: C0-C1, C0-C2, C1-C2, C1-C3, C2-C3 */
                        nOtherEdges += 5;
                        /* add connections of G0 to G1 and G2 */
                    }
                    continue;
                }

                nVertices ++; /* count BNS_VT_C_POS* types; all contain BNS_VERT_TYPE_C_GROUP bit */
                if ( !IS_BNS_VT_C_GR(type) ) {  /* check */
                    ret = RI_ERR_PROGR;
                    goto exit_function;
                }
                /* add st_cap and st_flow for a charge group */
                cap  = !bMetalAtoms? pCN[j].v.cap  : pCN[j].v.cap?  pSrm->nMetalMaxCharge_D : 0;
                flow = !bMetalAtoms? pCN[j].v.flow : pCN[j].v.flow? pSrm->nMetalMaxCharge_D : 0;
                ret = RegisterTCGroup( pTCGroups, type, 0 /* ord_num*/,
                                       cap /* st-cap*/, flow /* st-flow */,
                                       0 /* edge cap */, 0 /* edge flow */, 0 /* edges already counted */ );
                if ( ret < 0 ) {
                    goto exit_function;
                }
            }
        } else {
            pCN = NULL;
        }
        /* count edge caps to t-groups */
        if ( at[i].endpoint ) {
            int nEdgeCap = nTautEndpointEdgeCap( at, pVA, i );
            nTgroupEdges ++;
            if ( nEdgeCap < 0 ) {
                ret = nEdgeCap;
                goto exit_function;
            }
            /* add number of unsatisfied valences for a t-group; the unknown flow = 0 */
            ret = RegisterTCGroup( pTCGroups, BNS_VERT_TYPE_TGROUP, at[i].endpoint,
                                   0 /* st_cap */, 0 /* st_flow */,
                                   nEdgeCap /* edge cap */, 0 /* edge flow */,
                                   0 /* t-group edges have already been counted */ );
            if ( ret < 0 ) {
                goto exit_function;
            }

        }
    }
    if ( !bMetalAtoms && pTCGroups->num_metal_atoms ) {
        bMetalAtoms = 1;
        nBonds      = 0; /* added 2006-05-15 */
        goto repeat_for_metals;
    }

    /* count real atoms and bonds */
    nBonds /= 2;
    pTCGroups->num_atoms        = num_at;
    pTCGroups->num_bonds        = nBonds;

    pTCGroups->num_tgroups      = ti->num_t_groups;
    pTCGroups->num_tgroup_edges = nTgroupEdges;
    pTCGroups->tgroup_charge    = -nTotNegChargInTgroups;

    if ( 0 <= ret && nTgroupEdgesFromTg != nTgroupEdges ) {
        ret = BNS_PROGRAM_ERR;
    }

    nVertices += num_at;


    /* count other vertices */
    nVertices += ti->num_t_groups;
    nBonds    += nOtherEdges;

    /* return edges and vertices */
    pTCGroups->nVertices     = nVertices;
    pTCGroups->nEdges        = nBonds;

exit_function:
    return ret;
}

/****************************************************************
  int nAddSuperCGroups( ALL_TC_GROUPS *pTCGroups )

  1. adds BNS_VT_C_POS_ALL and BNS_VT_C_NEG_ALL ONLY if both
     {TCG_Plus0  and TCG_Plus_C0} and/or
     {TCG_Minus0 and TCG_Minus_C0} are present, respectively
    
  2. fills pTCGroups->nGroup[]:

  pTCGroups->nGroup[k] < 0  => does not exist
  pTCGroups->nGroup[k] = i  => the group is pTCGroups->pTCG[i]

    where           group           group 
      k =           type            number
    TCG_Plus0       BNS_VT_C_POS      0 
    TCG_Plus1,      BNS_VT_C_POS      1 
    TCG_Minus0,     BNS_VT_C_NEG      0 
    TCG_Minus1,     BNS_VT_C_NEG      1 
    TCG_Plus_C0,    BNS_VT_C_POS_C    0 
    TCG_Plus_C1,    BNS_VT_C_POS_C    1 
    TCG_Minus_C0,   BNS_VT_C_NEG_C    0 
    TCG_Minus_C1,   BNS_VT_C_NEG_C    1 
    TCG_Plus,       BNS_VT_C_POS_ALL  0 
    TCG_Minus,      BNS_VT_C_NEG_ALL  0 

only groups with number 0 are processed

  3. If only one of the groups in pairs mentioned in (1) above
     is present then 
     
     pTCGroups->nGroup[TCG_Plus] := pTCGroups->nGroup[TCG_Plus0] or
     pTCGroups->nGroup[TCG_Plus] := pTCGroups->nGroup[TCG_Plus_C0];
     an additional BNS_VT_C_POS_ALL vertex is not created

     same for pTCGroups->nGroup[TCG_Minus] and BNS_VT_C_NEG_ALL

  4. Adds to these new "supergroups" (TCG_Plus, TCG_Minus)
     descriptions in pTCGroups->pTCG[k]
     st_cap, st_flow, edges cap and flow from the corresponding
     groups {TCG_Plus0  and TCG_Plus_C0}. Same for the Minus groups.
     Stores indexes k in
     pTCGroups->nGroup[TCG_Plus], pTCGroups->nGroup[TCG_Minus]

 ****************************************************************/
int nAddSuperCGroups( ALL_TC_GROUPS *pTCGroups )
{
    int i, k, n, n1, n2, n3, nNumTg = 0, ret = 0, nNumToConnect;

    for ( i = 0; i < pTCGroups->num_tc_groups; i ++ ) {
        if ( pTCGroups->pTCG[i].type & BNS_VERT_TYPE_TGROUP ) {
            nNumTg ++;
            continue; /* t-group */
        }
        if ( IS_BNS_VT_C_GR(pTCGroups->pTCG[i].type) ||
             IS_BNS_VT_M_GR(pTCGroups->pTCG[i].type) ) {
            /* ChargeValence (cn) group */
            switch( pTCGroups->pTCG[i].type ) {
            case BNS_VT_C_POS:
                k = TCG_Plus0;
                break;
            case BNS_VT_C_NEG:
                k = TCG_Minus0;
                break;
            case BNS_VT_C_POS_C:
                k = TCG_Plus_C0;
                break;
            case BNS_VT_C_NEG_C:
                k = TCG_Minus_C0;
                break;
            case BNS_VT_C_POS_M:
                k = TCG_Plus_M0;
                break;
            case BNS_VT_C_NEG_M:
                k = TCG_Minus_M0;
                break;
            case BNS_VT_M_GROUP:
                switch( pTCGroups->pTCG[i].ord_num ) {
                case 0:
                    k = TCG_MeFlower0;
                    break;
                case 1:
                    k = TCG_MeFlower1;
                    break;
                case 2:
                    k = TCG_MeFlower2;
                    break;
                case 3:
                    k = TCG_MeFlower3;
                    break;
                default:
                    ret = RI_ERR_PROGR; /* unexpected group type */
                    goto exit_function;
                }
                break;

            default:
                ret = RI_ERR_PROGR; /* unexpected group type */
                goto exit_function;
            }
            if ( pTCGroups->nGroup[k] >= 0 || (pTCGroups->pTCG[i].ord_num && !IS_BNS_VT_M_GR(pTCGroups->pTCG[i].type)) ) {
                ret = RI_ERR_PROGR;
                goto exit_function;
            }
            pTCGroups->nGroup[k] = i; /* ordering number of the Charge group, starting from 0 */
        }
    }
    /* add (+) supergroup */
    n1 = pTCGroups->nGroup[TCG_Plus0];
    n2 = pTCGroups->nGroup[TCG_Plus_C0];
    n3 = pTCGroups->nGroup[TCG_Plus_M0];
    nNumToConnect = (n1>=0) + (n2>=0) + (n3>=0);
    if ( nNumToConnect ) {
        /* if both groups are present then add a supergroup */
        ret = RegisterTCGroup( pTCGroups, BNS_VT_C_POS_ALL, 0,
                               0 /* st_cap */,
                               0 /* st_flow */,
                               0 /* edge cap */,
                               0 /* edge flow */,
                               1+nNumToConnect /* one more edge to connect to an additional (+/-) vertex */ );

        if ( ret <= 0 ) {
            ret = !ret? RI_ERR_PROGR : ret;
            goto exit_function;
        }
        pTCGroups->nGroup[TCG_Plus] = ret - 1; /* newly added group number */
        pTCGroups->nVertices += 2; /* two vertices including itself */
        pTCGroups->nEdges    += 1 + nNumToConnect; /* one more edge to connect to an additional (+/-) vertex */
    }
    /* add (-) supergroup */
    n1 = pTCGroups->nGroup[TCG_Minus0];
    n2 = pTCGroups->nGroup[TCG_Minus_C0];
    n3 = pTCGroups->nGroup[TCG_Minus_M0];
    nNumToConnect = (n1>=0) + (n2>=0) + (n3>=0);
    if ( nNumToConnect ) {
        /* if both groups are present then add a supergroup */
        ret = RegisterTCGroup( pTCGroups, BNS_VT_C_NEG_ALL, 0,
                               0 /* st_cap */,
                               0 /* st_flow */,
                               0 /* edge cap */,
                               0 /* edge flow */,
                               1+nNumToConnect /* one more edge to connect to an additional (+/-) vertex */ );

        if ( ret < 0 ) {
            goto exit_function;
        }
        pTCGroups->nGroup[TCG_Minus] = ret - 1; /* newly added group number */
        pTCGroups->nVertices += 2; /* needs two vertices including itself */
        pTCGroups->nEdges    += 1 + nNumToConnect; /* one more edge to connect to an additional (+/-) vertex */
    }

    /* add neutralization vertex: (+)-()=(-) connection */
    k = pTCGroups->nGroup[TCG_Minus];
    n = pTCGroups->nGroup[TCG_Plus];
    nNumToConnect = (k>=0) + (n>=0);
    if ( nNumToConnect ) {
        pTCGroups->nVertices += 1;
        pTCGroups->nEdges    += nNumToConnect; /* one edge per super-c-group */
    }

    ret = 0;

exit_function:
    return ret;
}
/*********************************************************************************/
int AddTGroups2TCGBnStruct( BN_STRUCT *pBNS, StrFromINChI *pStruct, VAL_AT *pVA,
                            ALL_TC_GROUPS *pTCGroups, int nMaxAddEdges )
{
    int ret = 0;
    inp_ATOM *at        = pStruct->at;
    int       num_atoms = pStruct->num_atoms;
    int tot_st_cap, tot_st_flow;
    /* ret = ReInitBnStruct( pBNS ); */
    if ( pTCGroups->num_tgroups /* tgi && tgi->num_t_groups && tgi->t_group*/ ) {
        int         i, k, endpoint, /*centerpoint,*/ fictpoint;
        int         num_tg       = pTCGroups->num_tgroups;
        int         num_edges    = pBNS->num_edges;
        int         num_vertices = pBNS->num_vertices;
        BNS_VERTEX *vert_ficpoint, *vert_ficpoint_prev;  /* fictitious vertex describing t-group */
        BNS_VERTEX *vert_endpoint;
        BNS_EDGE   *edge;      /* edge between that vertex and the tautomeric endpoint */
        int        nMaxTGroupNumber = 0;
        /*ENDPOINT_INFO eif;*/

        /* Debug: check overflow */
        if ( num_vertices + num_tg >= pBNS->max_vertices ) {
            return BNS_VERT_EDGE_OVFL;
        }
        if ( num_edges + pTCGroups->num_tgroup_edges >= pBNS->max_edges ) {
            return BNS_VERT_EDGE_OVFL;
        }
        /* find the largest t-group ID */
        for ( i = 0; i < pTCGroups->num_tc_groups; i ++ ) {
            if ( pTCGroups->pTCG[i].type & BNS_VERT_TYPE_TGROUP ) {
                k = pTCGroups->pTCG[i].ord_num;
                if ( k <= 0 ) {
                    return BNS_CPOINT_ERR; /* t-group does not have a number or has a wrong number */
                }
                if ( k > pTCGroups->num_tc_groups ) {
                    return BNS_CPOINT_ERR; /* t-group has a wrong number */
                }
                if ( k != nMaxTGroupNumber + 1 ) {
                    return BNS_CPOINT_ERR; /* t-group numbers are not contiguously ascending */
                }
                nMaxTGroupNumber = k;
            } else {
                break; /* t-groups are contiguous and first in the list */
            }
        }
        if ( i != num_tg ) {
            return BNS_CPOINT_ERR; /* number of t-groups is wrong */
        }
        /* since t-group IDs may be not contiguous, clear all vertices that will be added.
           all-zeroes-vertex will be ignored by the BNS
        */
        memset( pBNS->vert+num_vertices, 0, nMaxTGroupNumber*sizeof(pBNS->vert[0]) );
        /* initialize new fictitious vertices */
        vert_ficpoint_prev = pBNS->vert+num_vertices - 1;

        tot_st_cap = tot_st_flow = 0;

        for ( i = 0; i < num_tg; i ++, vert_ficpoint_prev = vert_ficpoint ) {
            /*
              vert_ficpoint-1 is the last vertex;
              vert_ficpoint   is the vertex that is being added
              Note: nGroupNumber are not contiguous
            */
            vert_ficpoint                = pBNS->vert+num_vertices + pTCGroups->pTCG[i].ord_num - 1;
            vert_ficpoint->iedge         = vert_ficpoint_prev->iedge + vert_ficpoint_prev->max_adj_edges;
            vert_ficpoint->max_adj_edges = pTCGroups->pTCG[i].num_edges+nMaxAddEdges+BNS_ADD_SUPER_TGROUP;
            vert_ficpoint->num_adj_edges = 0;
            vert_ficpoint->st_edge.flow  = vert_ficpoint->st_edge.flow0  = 0;
            vert_ficpoint->st_edge.cap   = vert_ficpoint->st_edge.cap0   = pTCGroups->pTCG[i].st_cap;
            tot_st_cap                   += pTCGroups->pTCG[i].st_cap;
            vert_ficpoint->type          = pTCGroups->pTCG[i].type;
            pTCGroups->pTCG[i].nVertexNumber = vert_ficpoint - pBNS->vert;
        }

        for ( endpoint = 0; endpoint < num_atoms; endpoint ++ ) {
            if ( !at[endpoint].endpoint )
                continue;
            fictpoint = at[endpoint].endpoint + num_vertices - 1;
            vert_ficpoint = pBNS->vert + fictpoint;  /* t-group vertex */
            vert_endpoint = pBNS->vert + endpoint;   /* endpoint vertex */
            /* Debug: check overflow */
            if ( fictpoint >= pBNS->max_vertices ||
                 num_edges >= pBNS->max_edges    ||
                 vert_ficpoint->num_adj_edges >= vert_ficpoint->max_adj_edges ||
                 vert_endpoint->num_adj_edges >= vert_endpoint->max_adj_edges ) {
                ret = BNS_VERT_EDGE_OVFL;
                break;
            }
#ifdef NEVER
            /* obtain donor/acceptor info */
            if ( !nGetEndpointInfo( at, endpoint, &eif ) ) {
                ret = BNS_BOND_ERR;
                break;
            }
#endif
            vert_endpoint->type |= BNS_VERT_TYPE_ENDPOINT;
#ifdef NEVER
            /* set capacity = 1 to the edges from the endpoint to the centerpoint(s) */
            for ( k = 0; k < vert_endpoint->num_adj_edges; k ++ ) {
                int iedge = vert_endpoint->iedge[k];
                if ( !pBNS->edge[iedge].cap ) {
                    /* single bond, possibly between endpoint and centerpoint */
                    centerpoint = (pBNS->edge[iedge].neighbor12 ^ endpoint);
                    if ( centerpoint < pBNS->num_atoms &&
                         pBNS->vert[centerpoint].st_edge.cap >= 1 ) {
                        int bond_type = (at[endpoint].bond_type[k] & BOND_TYPE_MASK);
                        if (bond_type == BOND_TAUTOM  ||
                            bond_type == BOND_ALTERN  ||
                            bond_type == BOND_ALT12NS ||
                            bond_type == BOND_SINGLE ) {
                            pBNS->edge[iedge].cap = 1;
                        }
                    }
                }
            }
#endif
            /* create a new edge connecting endpoint to the new fictitious t-group vertex vert_ficpoint */
            edge = pBNS->edge + num_edges;
            edge->cap       = vert_endpoint->st_edge.cap - vert_endpoint->st_edge.flow;
            edge->cap       = inchi_min( edge->cap, MAX_TGROUP_EDGE_CAP );
            edge->cap       = inchi_max( edge->cap, 0 );
            edge->flow      = 0;
            edge->pass      = 0;
#if ( RESET_EDGE_FORBIDDEN_MASK == 1 )
            edge->forbidden &= pBNS->edge_forbidden_mask;
#endif

#ifdef NEVER
            /* later include case when the charge change allows the endpoint to become tautomeric */
            /* mark endoint having moveable H atom with flow=1 */

            /* -- old "no charges" version -- */
            /* if (at[endpoint].chem_bonds_valence == at[endpoint].valence) */
            /* -- the following line takes charges into account -- */
            if ( eif.cDonor ) /* means the endpoint has an H-atom to donate */
            {
                /* increment edge flow */
                edge->flow ++;
                /* increment one vertex st-flow & cap */
                vert_ficpoint->st_edge.flow ++;
                vert_ficpoint->st_edge.cap ++;
                /* increment another vertex st-flow & cap */
                vert_endpoint->st_edge.flow ++;
                vert_endpoint->st_edge.cap ++;
            }
#endif
            /* connect edge to endpoint and fictpoint and increment the counters of neighbors and edges */
            ret = ConnectTwoVertices( vert_endpoint, vert_ficpoint, edge, pBNS, 0 );
            if ( IS_BNS_ERROR( ret ) ) {
                break;
            }
            num_edges ++;
            edge->cap0  = edge->cap;
            edge->flow0 = edge->flow;
            pVA[endpoint].nTautGroupEdge = num_edges; /* edge index + 1 */
        }

        pBNS->num_edges     = num_edges;
        pBNS->num_vertices += nMaxTGroupNumber;
        pBNS->num_t_groups  = num_tg;
        pBNS->tot_st_cap   += tot_st_cap;
        pBNS->tot_st_flow  += tot_st_flow;

    }
    return ret;
}
/*****************************************************************************************************/
int ConnectTwoVertices( BNS_VERTEX *p1, BNS_VERTEX *p2, BNS_EDGE *e, BN_STRUCT *pBNS, int bClearEdge )
{
    int ip1 = p1 - pBNS->vert;
    int ip2 = p2 - pBNS->vert;
    int ie  = e  - pBNS->edge;
    /* debug: check bounds */
    if ( ip1 >= pBNS->max_vertices || ip1 < 0 ||
         ip2 >= pBNS->max_vertices || ip2 < 0 ||
         ie  >= pBNS->max_edges    || ie  < 0 ||
         (p1->iedge - pBNS->iedge) < 0 ||
         (p1->iedge - pBNS->iedge) + p1->max_adj_edges > pBNS->max_iedges ||
         (p2->iedge - pBNS->iedge) < 0 ||
         (p2->iedge - pBNS->iedge) + p2->max_adj_edges > pBNS->max_iedges ||
         p1->num_adj_edges >= p1->max_adj_edges ||
         p2->num_adj_edges >= p2->max_adj_edges  ) {
        return BNS_VERT_EDGE_OVFL;
    }
    /* clear the edge */
    if ( bClearEdge ) {
        memset( e, 0, sizeof(*e) );
    } else
    if ( e->neighbor1 || e->neighbor12 ) {
        return BNS_PROGRAM_ERR;
    }
    /* connect */
    e->neighbor1  = inchi_min( ip1, ip2 );
    e->neighbor12 = ip1 ^ ip2;
    p1->iedge[p1->num_adj_edges] = ie;
    p2->iedge[p2->num_adj_edges] = ie;
    e->neigh_ord[ip1 > ip2] = p1->num_adj_edges ++;
    e->neigh_ord[ip1 < ip2] = p2->num_adj_edges ++;
    return 0;
}

/***********************************************************************************************************
                     METAL ATOMS' FLOWER - Provides a source/sink of "free valences"
 ***********************************************************************************************************

                    c1+...+cn = 2c+dc  - total cap and flow of edges to the flower base from metal atoms
                    f1+...+fn = 2f+df    they should allow changing bonds to metals from 0-order to triple
                    dc,df = 0 or 1       hence c=3*n, f=0 (initial zero bond order) or n
 Gi=vertex(M-group)
 Ci=its st_cap                  [C3,F3]  C0 = F0 = 2c + 2D + dc    (st_cap & st_flow)
 Fi=its st_flow               G3         C2 = F2 =  c + 2D                           
                              / \        C1 = F1 =  c + 2D + dc-df                   
 ci=cap of edge i       cx,fx/   \cy,fy  C3 = F3 =  0                                 
 fi=edge flow               /     \                                             Constraints
                    [C2,F2]/ cd,fd \[C1,F1]                                     -----------------
                         G2--------G1                                           fa+fb+2f+df=F0=C0
                          \        /          ca =  c + 2D         (edge cap)   fa+fd      =F2=C2
                      ca,fa \    / cb,fb      fa =  c +  D - f     (edge flow)  fb+fd      =C1=F1
                             \  /                                               fi <= ci
                              G0 [C0,F0]      cb =  c + 2D + dc                 -----------------
                              /\              fb =  c +  D + dc - (f + df)    
   ci=3, fi=0 or 1   c1,f1  /... \ cn,fn                            ------------------------------------
                          /       \           cd =  c + 2D          D is an arbitrary integer > 0
   all n Metal atoms:    M1 ...    Mn         fd =  f +  D          it allows one to apply
                                                                    C3++ (add st_flow to cancel radicals)
  For each Mi add cap and flow=cap            cx = cy = D           D times.
  to M-charge group                           fx = fy = 0
       --------------------------------------------------------------------------------------
                 |  f=0                   |  f=c, dc>=df         |  0 <= 2f+df <= 2c+dc
         edge    +------------------------+-----------+----------+-------------+-------------
                 |  flow      |  rescap   |  flow     |  rescap  |  flow       |  rescap  
       ----------+------------+-----------+-----------+----------+-------------+-------------
        f1+..+fn |  df        |  2c+dc-df |  2c+df    |  dc-df   |  2f+df      |  2c-2f+dc-df
              fa |  c+D       |  D        |  D        |  c+D     |  c+D-f      |  c+D     
              fb |  c+D+dc-df |  D+df     |  D+dc-df  |  c+D+df  |  c+D+dc-f-df|  c+D+df  
              fd |  D         |  c+D      |  c+D      |  D       |  f+D        |  D       
       --------------------------------------------------------------------------------------
***********************************************************************************************************/
int AddRadicalToMetal( int *tot_st_cap, int *tot_st_flow, ICHICONST SRM *pSrm, BN_STRUCT *pBNS, ALL_TC_GROUPS *pTCGroups )
{
    int iG0 = pTCGroups->nGroup[TCG_MeFlower0]; /* index in pTCGroups->pTCG[] */
    int iG1 = pTCGroups->nGroup[TCG_MeFlower1];
    int iG2 = pTCGroups->nGroup[TCG_MeFlower2];
    int iG3 = pTCGroups->nGroup[TCG_MeFlower3];
    int n   = (iG0>=0) + (iG1>=0) + (iG2>=0) + (iG3>=0);
    int vG0, vG1, vG2, vG3;  /* M-vertex number */
    BNS_VERTEX *pG0=NULL, *pG1=NULL, *pG2=NULL, *pG3=NULL;

    if ( pTCGroups->num_metal_atoms &&
         pSrm->bMetalAddFlower      &&
         *tot_st_cap % 2            &&
         n == 4 ) {
        vG0 = pTCGroups->pTCG[iG0].nVertexNumber;
        vG1 = pTCGroups->pTCG[iG1].nVertexNumber;
        vG2 = pTCGroups->pTCG[iG2].nVertexNumber;
        vG3 = pTCGroups->pTCG[iG3].nVertexNumber;

        pG0 = pBNS->vert+vG0;
        pG1 = pBNS->vert+vG1;
        pG2 = pBNS->vert+vG2;
        pG3 = pBNS->vert+vG3;

        /* add 1 unit to metal flower st_cap */
        pG3->st_edge.cap  ++;
        pG3->st_edge.cap0 ++;
        (*tot_st_cap)     ++;
        return 1;
    }
    return 0;
}
/***********************************************************************************************************/
int ConnectMetalFlower( int *pcur_num_vertices, int *pcur_num_edges,
                        int *tot_st_cap, int *tot_st_flow, ICHICONST SRM *pSrm,
                        BN_STRUCT *pBNS, ALL_TC_GROUPS *pTCGroups )
{
    int iG0 = pTCGroups->nGroup[TCG_MeFlower0]; /* index in pTCGroups->pTCG[] */
    int iG1 = pTCGroups->nGroup[TCG_MeFlower1];
    int iG2 = pTCGroups->nGroup[TCG_MeFlower2];
    int iG3 = pTCGroups->nGroup[TCG_MeFlower3];
    int n   = (iG0>=0) + (iG1>=0) + (iG2>=0) + (iG3>=0);
    int vG0, vG1, vG2, vG3;  /* M-vertex number */
    int cur_num_edges    = *pcur_num_edges;
    int cur_num_vertices = *pcur_num_vertices;
    BNS_VERTEX *pG0=NULL, *pG1=NULL, *pG2=NULL, *pG3=NULL;
    BNS_EDGE   *ea=NULL, *eb=NULL, *ed=NULL, *ex=NULL, *ey=NULL, *e;
    int         ia, ib, id, ix, iy;
    int         c, f, dc, df, ca, fa, cb, fb, cd, fd, cx, fx, cy, fy;
    int         C0, F0, C1, F1, C2, F2, C3, F3, D;
    int         ret = 0, i;

    if ( 0 == n ) {
        goto exit_function;
    }
    if ( 4 != n ) {
        ret = RI_ERR_PROGR;
        goto exit_function;
    }
    vG0 = pTCGroups->pTCG[iG0].nVertexNumber;
    vG1 = pTCGroups->pTCG[iG1].nVertexNumber;
    vG2 = pTCGroups->pTCG[iG2].nVertexNumber;
    vG3 = pTCGroups->pTCG[iG3].nVertexNumber;

    pG0 = pBNS->vert+vG0;
    pG1 = pBNS->vert+vG1;
    pG2 = pBNS->vert+vG2;
    pG3 = pBNS->vert+vG3;

    /* count G0 edges cap and flow (currently only atoms are connected to G0) */
    for ( i = 0, c = 0, f = 0; i < pG0->num_adj_edges; i ++ ) {
        e = pBNS->edge + pG0->iedge[i];
        c += e->cap;
        f += e->flow;
    }

    /* consistency checks */
    if ( !IS_BNS_VT_M_GR(pTCGroups->pTCG[iG0].type) &&
         (pTCGroups->pTCG[iG0].edges_cap != pG0->st_edge.cap ||
          pTCGroups->pTCG[iG0].edges_flow != pG0->st_edge.flow) ) {
        ret = RI_ERR_PROGR;
        goto exit_function;
    }
    if ( pTCGroups->pTCG[iG0].edges_cap  != c ||
         pTCGroups->pTCG[iG0].edges_flow != f ) {
        ret = RI_ERR_PROGR;
        goto exit_function;
    }
    
    /* get new edges */

    ea = pBNS->edge + (ia=cur_num_edges++);
    eb = pBNS->edge + (ib=cur_num_edges++);
    ed = pBNS->edge + (id=cur_num_edges++);
    ex = pBNS->edge + (ix=cur_num_edges++);
    ey = pBNS->edge + (iy=cur_num_edges++);

    /* connect vertices with edges */
    ret = ConnectTwoVertices( pG0, pG1, eb, pBNS, 1 );
    if ( IS_BNS_ERROR( ret ) ) {
        goto exit_function;
    }
    ret = ConnectTwoVertices( pG0, pG2, ea, pBNS, 1 );
    if ( IS_BNS_ERROR( ret ) ) {
        goto exit_function;
    }
    ret = ConnectTwoVertices( pG1, pG2, ed, pBNS, 1 );
    if ( IS_BNS_ERROR( ret ) ) {
        goto exit_function;
    }
    ret = ConnectTwoVertices( pG1, pG3, ey, pBNS, 1 );
    if ( IS_BNS_ERROR( ret ) ) {
        goto exit_function;
    }
    ret = ConnectTwoVertices( pG2, pG3, ex, pBNS, 1 );
    if ( IS_BNS_ERROR( ret ) ) {
        goto exit_function;
    }

    /* calculate caps and flows */

    dc  = c % 2;
    c  /= 2;
    df  = f % 2;
    f  /= 2;

    D = pSrm->nMetalFlowerParam_D;
    
    C0 = F0 = 2*c + 2*D + dc;
    C1 = F1 =   c + 2*D + dc - df;
    C2 = F2 =   c + 2*D;
    C3 = F3 =   0;
    
    ca =  c + 2*D;
    fa =  c +   D - f;
    
    cb =  c + 2*D + dc;
    fb =  c +   D + dc - ( f + df );
    
    cd =  c + 2*D;
    fd =  f +   D;
    
    cx = cy = D;
    fx = fy = 0;

    /* check overflow */
    if ( C0 >= EDGE_FLOW_ST_MASK || F0 >= EDGE_FLOW_ST_MASK ||
         C1 >= EDGE_FLOW_ST_MASK || F1 >= EDGE_FLOW_ST_MASK ||
         C2 >= EDGE_FLOW_ST_MASK || F2 >= EDGE_FLOW_ST_MASK ||
         C3 >= EDGE_FLOW_ST_MASK || F3 >= EDGE_FLOW_ST_MASK ) {
        return BNS_PROGRAM_ERR; /* cannot handle too large st-cap or st-flow */
    }

    /* set st caps and flows */

    SetStCapFlow( pG0, tot_st_flow, tot_st_cap, C0, F0 );
    SetStCapFlow( pG1, tot_st_flow, tot_st_cap, C1, F1 );
    SetStCapFlow( pG2, tot_st_flow, tot_st_cap, C2, F2 );
    SetStCapFlow( pG3, tot_st_flow, tot_st_cap, C3, F3 );

    SetEdgeCapFlow( ea, ca, fa );
    SetEdgeCapFlow( eb, cb, fb );
    SetEdgeCapFlow( ed, cd, fd );
    SetEdgeCapFlow( ex, cx, fx );
    SetEdgeCapFlow( ey, cy, fy );


    *pcur_num_edges    = cur_num_edges;
    *pcur_num_vertices = cur_num_vertices;

    ret = 0;
    
exit_function:

    return ret;
}
/********************************************************************************/
void SetEdgeCapFlow( BNS_EDGE *e, int edge_cap, int edge_flow )
{
    e->cap  = e->cap0  = edge_cap;
    e->flow = e->flow0 = edge_flow;
}

/*********************************************************************************
  Add cap and flow to an edge
  Add edge flow to the source vertex st_flow
  Add edge cap & flow to the destination vertex cap and flow
 *********************************************************************************/
int AddEdgeFlow( int edge_cap, int edge_flow, BNS_EDGE *e01, BNS_VERTEX *pSrc /*src*/,
                  BNS_VERTEX *pDst/*dest*/, int *tot_st_cap, int *tot_st_flow )
{
    /* overflow chaeck */
    if ( e01->cap < 0 || edge_cap < 0 ||  (int)e01->cap + edge_cap >= EDGE_FLOW_MASK ) {
        return BNS_PROGRAM_ERR;
    }
    if ( pDst->st_edge.cap  < 0 || (int)pDst->st_edge.cap  + edge_cap  >= EDGE_FLOW_ST_MASK ||
         pDst->st_edge.flow < 0 || (int)pDst->st_edge.flow + edge_flow >= EDGE_FLOW_ST_MASK ||
         pSrc->st_edge.cap  < 0 || pSrc->st_edge.flow < 0 ||
                                   (int)pSrc->st_edge.flow + edge_flow >= EDGE_FLOW_ST_MASK ) {
        return BNS_PROGRAM_ERR;
    }
    /* add flow */
    e01->cap    += edge_cap;
    e01->flow   += edge_flow;
    e01->cap0    = e01->cap;
    e01->flow0   = e01->flow;

    pDst->st_edge.cap  += edge_cap;
    pDst->st_edge.cap0  = pDst->st_edge.cap;
    *tot_st_cap        += edge_cap;

    pDst->st_edge.flow += edge_flow;
    pDst->st_edge.flow0 = pDst->st_edge.flow;
    *tot_st_flow       += edge_flow;

    pSrc->st_edge.flow += edge_flow;
    pSrc->st_edge.flow0 = pSrc->st_edge.flow;
    *tot_st_flow       += edge_flow;

/*
    pDst->st_edge.cap  += e01->cap;
    pDst->st_edge.cap0  = pDst->st_edge.cap;
    *tot_st_cap       += e01->cap;

    pDst->st_edge.flow += e01->flow;
    pDst->st_edge.flow0 = pDst->st_edge.flow;
    *tot_st_flow      += e01->flow;

    pSrc->st_edge.flow += e01->flow;
    pSrc->st_edge.flow0 = pSrc->st_edge.flow;
    *tot_st_flow      += e01->flow;
*/
    return 0;
}
/**************************************************************
     (+) and (-) group V - connection
     ================================

  BNS_VERT_TYPE__AUX    (+/-)-connection
                    (v)  st_cap  = 
                    / \  st_flow = (cap0 - Delta0 - flow0) + (cap1 - Delta1 -flow1)
                   /   \   
                  /     \   cap  =  cap1
                 /       \  flow = (cap1 - Delta1 - flow1)
                /         \
              (-)         (+)  st_cap  = cap1
             /   \       /   \ st_flow = cap1 - Delta1
            /cap0 \     /cap1 \
             flow0       flow1
    
***************************************************************

     (+) supergroup Y - connection
     ==============================
     
                           (+) BNS_VT_C_POS_ALL  (+) supergroup
      Delta0                |                    ==============
      not shown             |  cap  = cap0+cap1
                            |  flow = flow0+flow1-Delta0-Delta1
       BNS_VERT_TYPE__AUX  (y) <------------------ additional vertex: st_cap  = cap0+cap1
                           / \                                        st_flow = cap0+cap1
         cap=cap0         /   \  cap  = cap1
         flow=cap0-flow0 /     \ flow = cap1 - flow1 - Delta1
             -Delta0    /       \
                not-C (+)       (+) Carbons         st_cap  = cap1
         BNS_VT_C_POS / \       / \ BNS_VT_C_POS_C  st_flow = cap1 - Delta1
                     /   \     /   \
        totals      cap0        cap1 = sum of all cap  going up into (+) from atoms or ChargeStruct
        to (+):     flow0       flow1= sum of all flow going up into (+)
                   Delta0       Delta1 = st_cap(+)-st_flow(+) before connection
  Observations
  ============
  A. Any Delta > 0 on (+) or (-) group decreases total (signed) charge by Delta

  B. Any alt path from an atom through ChargeStruct to an atom
     does not change the total charge

  C. st_flow(+/-) = cap(+)-flow(+)-Delta(+) + cap(-)-flow(-)-Delta(-) =
                  = charge(+) + |max (-) charge| + charge(-) = const
     (charge conservation)

  D. To decrease total charge: increase st_cap on (+) or (-) group, including supergroup
  E. To increase total charge: increase st_cap on any (y) or (v)-connecting vertex

  F. To cancel charges: 
       1. Forbid (+/-)-(+) or (+/-)-(-) edge
       2. Add delta>0 to (+/-) st_cap
       3. Add same delta to (+) or (-) st_cap


****************************************************************/

/************************************************************************************
     j2,j3 < j1 < j0

                                (+/-) <---- next step; if does not exist then 
                                /   \                  st_cap1' := st_flow1'
                               /     \
                              /       \ st_cap1' := cap01'
                        pv1 (+)super    st_flow1':= cap01'-flow01' = flow02'+flow03'     
                        j1   |
                             |           cap01' := st_cap0'
                             |           flow01':= st_cap0'-flow02'-flow03'
                             |
                             |           st_cap0' :=
                            ( )  pv0,j0  st_flow0':= cap2+st_cap3
                           /   \
                          /     \   cap03'  = cap3
                         /       \  flow03' = cap3 - flow3 - Delta3
                        /         \
   st_cap2, st_flow2  (+)         (+C)        st_cap3'  := cap3
                      pv2,j2      pv3,j3      st_flow3' := cap3-Delta3
                    /    \       /    \       Delta3    := st_cap3 - st_flow3
                cap2, flow2     cap3, flow3 = sums of incoming
 **************************************************************************************/
 
int ConnectSuperCGroup( int nSuperCGroup, int nAddGroups[], int num_add,
                        int *pcur_num_vertices, int *pcur_num_edges,
                        int *tot_st_cap, int *tot_st_flow,
                        BN_STRUCT *pBNS, ALL_TC_GROUPS *pTCGroups )
{
    BNS_EDGE   **e0X = NULL, *e;
    BNS_VERTEX **pvX = NULL, *pv0=NULL, *pv1=NULL, *pv=NULL;
    int         *jX  = NULL, *iX = NULL;
    int        i, j, num_groups, j0, i1, j1, iXX, ret = 0, fst=0;
    int        cur_num_vertices = *pcur_num_vertices;
    int        cur_num_edges    = *pcur_num_edges;

    if ( nSuperCGroup >= 0 ) {
        i1 = pTCGroups->nGroup[nSuperCGroup]; /* the supergroup */
        if ( i1 < 0 )
            return 0;
    } else {
        i1 = -1;
        fst = 1;
    }

    for ( i = num_groups = 0; i < num_add; i ++ ) {
        iXX = pTCGroups->nGroup[nAddGroups[i]];
        num_groups += (iXX >= 0 && iXX != i1);
    }
    if ( num_groups < 1 ) {  /* Y connect only 2 or more groups; V connects even 1 group */
        return 0;
    }

    e0X = (BNS_EDGE   **)inchi_calloc( num_groups + 1, sizeof(e0X[0]) );
    pvX = (BNS_VERTEX **)inchi_calloc( num_groups + 1, sizeof(pvX[0]) );
    jX  = (int         *)inchi_calloc( num_groups + 1, sizeof(jX[0]) );
    iX  = (int         *)inchi_calloc( num_groups + 1, sizeof(iX[0]) );
    if ( !e0X || !pvX || !jX || !iX ) {
        ret = RI_ERR_ALLOC;
        goto exit_function;
    }
    /* create vert_ficpoint -- central Y-connection vertex */
    j0 = cur_num_vertices;
    pv0 = pBNS->vert + j0; /* center of the Y-connection; has number j0 */
    pv0->iedge = (pv0 - 1)->iedge + (pv0 - 1)->max_adj_edges;
    pv0->max_adj_edges = num_groups + 1 + BNS_ADD_EDGES; /* Y-connection num. edges */
    pv0->num_adj_edges = 0; /* nothing connected yet */
    pv0->type          = BNS_VT_YVCONNECTOR;
    cur_num_vertices ++;

    if ( fst == 0 ) {
        /* find super c-group vertex pv1, number j1 */
        jX[0]  = j1 = pTCGroups->pTCG[i1].nVertexNumber;
        iX[0]  = i1;
        pvX[0] = pv1 = pBNS->vert + j1;
    }
    /* find other c-group vertices */
    for( i = 0, j = 1; i < num_add; i ++ ) {
        iXX = pTCGroups->nGroup[nAddGroups[i]];
        if ( (iXX >= 0) && (iXX != i1) ) {
            iX[j] = iXX;
            jX[j] = pTCGroups->pTCG[iXX].nVertexNumber;
            pvX[j] = pBNS->vert + jX[j];
            j ++;
        }
    }

    /* grab (num_groups+1) free edges */
    for ( i = fst; i <= num_groups; i ++ ) {
        e = e0X[i] = pBNS->edge + cur_num_edges;
        pv = pvX[i];
        j  = jX[i];
        iXX = iX[i];
        /* connect all to pv0 */
        ret = ConnectTwoVertices( pv0, pv, e, pBNS, 1 );
        if ( IS_BNS_ERROR( ret ) ) {
            goto exit_function;
        }
        if ( i ) {
            /* from c-group to central Y-connecting vertex of from supergroup to (+/-) vertex */
            pTCGroups->pTCG[iX[i]].nForwardEdge = cur_num_edges;
        } else {
            /* from central Y-connecting vertex to supergroup */
            pTCGroups->pTCG[iX[i]].nBackwardEdge = cur_num_edges;
        }
        cur_num_edges ++;
    }
    /* set flow and cap for incoming into pv0 edges */
    for ( i = 1; i <= num_groups; i ++ ) {
        int nDelta    = pTCGroups->pTCG[iX[i]].st_cap - pTCGroups->pTCG[iX[i]].edges_cap;
        int edge_cap  = pTCGroups->pTCG[iX[i]].edges_cap + nDelta; /* added nDelta */
        int edge_flow = pTCGroups->pTCG[iX[i]].edges_cap-pTCGroups->pTCG[iX[i]].edges_flow /*-nDelta*/;
        ret = AddEdgeFlow( edge_cap, edge_flow,
                     e0X[i], pvX[i]/*src*/, pv0 /* dest*/, tot_st_cap, tot_st_flow );
        if ( IS_BNS_ERROR( ret ) ) {
            goto exit_function;
        }
    }
    if ( fst == 0 ) {
        /* set flow and cap for going out of pv0 and into pv1 edge */
        int edge_cap  = pv0->st_edge.cap;
        int edge_flow = pv0->st_edge.cap - pv0->st_edge.flow;
        ret = AddEdgeFlow( pv0->st_edge.cap, pv0->st_edge.cap - pv0->st_edge.flow,
                     e0X[0], pv0/*src*/, pv1 /* dest*/, tot_st_cap, tot_st_flow );
        if ( IS_BNS_ERROR( ret ) ) {
            goto exit_function;
        }
        pTCGroups->pTCG[iX[0]].edges_cap  += edge_cap;
        pTCGroups->pTCG[iX[0]].edges_flow += edge_flow;
        pTCGroups->pTCG[iX[0]].st_cap     += edge_cap;
        pTCGroups->pTCG[iX[0]].st_flow    += edge_flow;
    } else {
        /* no supergroup => change cap to flow */
        *tot_st_cap      += pv0->st_edge.flow - pv0->st_edge.cap;
        pv0->st_edge.cap += pv0->st_edge.flow - pv0->st_edge.cap;
        pv0->st_edge.cap0 = pv0->st_edge.cap;
    }

    *pcur_num_vertices = cur_num_vertices;
    *pcur_num_edges    = cur_num_edges;
    ret = num_groups;
exit_function:
    if ( e0X ) inchi_free( e0X );
    if ( pvX ) inchi_free( pvX );
    if ( jX  ) inchi_free( jX  );
    if ( iX  ) inchi_free( iX  );
    return ret;
}
/*********************************************************************************/
void AddStCapFlow( BNS_VERTEX *vert_ficpoint, int *tot_st_flow, int *tot_st_cap, int cap, int flow )
{
    vert_ficpoint->st_edge.flow += flow;
    *tot_st_flow                += flow;
    vert_ficpoint->st_edge.cap  += cap;
    *tot_st_cap                 += cap;

    vert_ficpoint->st_edge.flow0 = vert_ficpoint->st_edge.flow;
    vert_ficpoint->st_edge.cap0  = vert_ficpoint->st_edge.cap;
}
/*********************************************************************************/
void SetStCapFlow( BNS_VERTEX *vert_ficpoint, int *tot_st_flow, int *tot_st_cap, int cap, int flow )
{
    *tot_st_flow                += flow - vert_ficpoint->st_edge.flow;
    vert_ficpoint->st_edge.flow  = flow;
    *tot_st_cap                 += cap - vert_ficpoint->st_edge.cap;
    vert_ficpoint->st_edge.cap   = cap;

    vert_ficpoint->st_edge.flow0 = vert_ficpoint->st_edge.flow;
    vert_ficpoint->st_edge.cap0  = vert_ficpoint->st_edge.cap;
}

/*********************************************************************************
int AddCGroups2TCGBnStruct( BN_STRUCT *pBNS, StrFromINChI *pStruct,
                            VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups )


 *********************************************************************************/
int AddCGroups2TCGBnStruct( BN_STRUCT *pBNS, StrFromINChI *pStruct, VAL_AT *pVA,
                            ALL_TC_GROUPS *pTCGroups, int nMaxAddEdges )
{
    int ret = 0, ret1, ret2, ret3, bNeedsFlower;
    inp_ATOM *at               = pStruct->at;
    int       num_atoms        = pStruct->num_atoms;
    /*int       num_tg           = pTCGroups->num_tgroups;*/
    int       num_cg           = pTCGroups->num_tc_groups - pTCGroups->num_tgroups;
    int       fst_cg_vertex    = pBNS->num_vertices;
    int       fst_cg_group     = pTCGroups->num_tgroups;
    int       num_vertices     = pBNS->num_vertices;
    int       num_edges        = pBNS->num_edges;
    int       cg_charge        = 0;
    ICHICONST SRM *pSrm        = pStruct->pSrm;
    /* ret = ReInitBnStruct( pBNS ); */
    if ( num_cg > 0 ) {
    /* if ( cgi && cgi->num_c_groups && cgi->c_group ) */
        int         i, i1, i2, j, j1, j2, k, k1, k2, n, c_point, c_neigh, cap, flow;
        int         cur_num_vertices, cur_num_edges;
        BNS_VERTEX *vert_ficpoint, *vert_ficpoint_prev, *vert_ficpoint_base;  /* fictitious vertex describing charge c-group */
        BNS_VERTEX *pv1, *pv2;
        BNS_EDGE   *edge;      /* edge between that vertex and the tautomeric c_point */
        int        nMaxCGroupNumber = 0;
        MY_CONST C_NODE *pCN;
        int              cn_len, cn_bits, bMetalAtoms;
        int              type;
        int        tot_st_cap, tot_st_flow;
        int        nAddGroups[16];


        /* Debug: check overflow */
        if ( num_vertices >= pBNS->max_vertices ) {
            return BNS_VERT_EDGE_OVFL;
        }
        nMaxCGroupNumber = num_cg;
        /* clear all vertices not used until now */
        memset( pBNS->vert+num_vertices, 0, (pBNS->max_vertices - num_vertices)*sizeof(pBNS->vert[0]) );
        tot_st_cap  = pBNS->tot_st_cap;
        tot_st_flow = pBNS->tot_st_flow;
        /*****************************************/
        /* initialize new fictitious vertices    */
        /* representing c-point groups, c-groups */
        /*****************************************/
        vert_ficpoint_prev = pBNS->vert + fst_cg_vertex - 1;
        
        for ( i = 0; i < num_cg; i ++ ) {
            /*
              vert_ficpoint-1 is the last vertex;
              vert_ficpoint   is the being added vertex
              Note: nGroupNumber are not contiguous
            */
            vert_ficpoint                = vert_ficpoint_prev + 1;
            vert_ficpoint->iedge         = vert_ficpoint_prev->iedge + vert_ficpoint_prev->max_adj_edges;
            vert_ficpoint->max_adj_edges = pTCGroups->pTCG[i+fst_cg_group].num_edges+nMaxAddEdges;
            vert_ficpoint->num_adj_edges = 0;

            vert_ficpoint->st_edge.flow += pTCGroups->pTCG[i+fst_cg_group].st_flow;
            tot_st_flow                 += pTCGroups->pTCG[i+fst_cg_group].st_flow;
            vert_ficpoint->st_edge.cap  += pTCGroups->pTCG[i+fst_cg_group].st_cap;
            tot_st_cap                  += pTCGroups->pTCG[i+fst_cg_group].st_cap;
            
            vert_ficpoint->st_edge.flow0 = vert_ficpoint->st_edge.flow;
            vert_ficpoint->st_edge.cap0  = vert_ficpoint->st_edge.cap;
            
            vert_ficpoint->type          = pTCGroups->pTCG[i+fst_cg_group].type;
            /* save the vertex number */
            pTCGroups->pTCG[i+fst_cg_group].nVertexNumber = vert_ficpoint - pBNS->vert;

            vert_ficpoint_prev             = vert_ficpoint; /* keep track of iedges */
        }
        cur_num_vertices = (vert_ficpoint_prev - pBNS->vert) + 1;
        cur_num_edges    = num_edges;
        
        /*************************************************************/
        /* pass 1:                                                   */
        /* create ChargeStruct for c-points and connect them to      */
        /* the vertices representing c-point groups;                 */
        /* set final atom st_cap, st_flow                            */
        /*************************************************************/
        for ( c_point = 0; c_point < num_atoms; c_point ++ ) {
            if ( !(k=pVA[c_point].cnListIndex) )
                continue;  /* not a c-point */
            k --;
            pCN     = cnList[k].pCN;   /* pointer to the ChargeStruct */
            cn_len  = cnList[k].len;   /* length of the ChargeStruct  */
            cn_bits = cnList[k].bits;  /* bits: for M-recognition */
            /* cn_bits = cnList[k].bits; */ /* ChargeStruct type */
            bMetalAtoms = (cn_bits == cn_bits_Me);
            vert_ficpoint_base = vert_ficpoint_prev; /* add aux vertices after this */
            /* create disconnected auxiliary vertices of the at[c_point] ChargeStruct; add to them st_flow & st_cap */
            for ( i1 = 0; i1 < cn_len; i1 ++ ) {
                if ( !IS_BNS_VT_CHRG_STRUCT(pCN[i1].v.type) ) {
                    continue;
                }
                /* the atom is always the first; the attached c-points are always the last */
                vert_ficpoint                = vert_ficpoint_base + i1; /* i1 = 1, 2,.. less number of attached c-points */
                vert_ficpoint->iedge         = vert_ficpoint_prev->iedge + vert_ficpoint_prev->max_adj_edges;
                vert_ficpoint->max_adj_edges = pCN[i1].v.valence; /* do not add additional edges to aux vertices */
                vert_ficpoint->num_adj_edges = 0;

                cap  = !bMetalAtoms? pCN[i1].v.cap : pCN[i1].v.cap? pSrm->nMetalMaxCharge_D : 0;
                flow = !bMetalAtoms? pCN[i1].v.flow : pCN[i1].v.flow? pSrm->nMetalMaxCharge_D : 0;

                AddStCapFlow( vert_ficpoint, &tot_st_flow, &tot_st_cap, cap, flow );
                vert_ficpoint->type          = pCN[i1].v.type; /* =BNS_VERT_TYPE__AUX */

                vert_ficpoint_prev = vert_ficpoint; /* the last one will be vert_ficpoint for the next c-point */
                cur_num_vertices   = (vert_ficpoint - pBNS->vert) + 1;

                if ( vert_ficpoint->iedge + vert_ficpoint->max_adj_edges - pBNS->iedge >= pBNS->max_iedges ) {
                    return BNS_VERT_EDGE_OVFL;
                }
                if ( cur_num_vertices >= pBNS->max_vertices ) {
                    return BNS_VERT_EDGE_OVFL;
                }
            }
            /* connect the vertices with new edges, add edge flow and cap */
            for ( i1 = 0; i1 < cn_len; i1 ++ ) {
                pv1 = NULL;
                k1  = -1;
                /* find vertex cooresponding to i1 */
                if ( pCN[i1].v.type & BNS_VERT_TYPE_ATOM ) {
                    pv1 = pBNS->vert+c_point; /* may be only one atom -- the current c_point at i1==0 */
                    /* add atom vertex st_cap and st_flow */
                    cap  = !bMetalAtoms? pCN[i1].v.cap : pCN[i1].v.cap? pSrm->nMetalMaxCharge_D : 0;
                    flow = !bMetalAtoms? pCN[i1].v.flow : pCN[i1].v.flow? pSrm->nMetalMaxCharge_D : 0;
                    AddStCapFlow( pv1, &tot_st_flow, &tot_st_cap, cap, flow );
                } else
                if ( IS_BNS_VT_C_GR(pCN[i1].v.type) ) {
                    /* find c-group vertex by looking for its type */
                    for( j = 0; j < num_cg; j ++ ) {
                        if ( pCN[i1].v.type == pBNS->vert[fst_cg_vertex + j].type ) {
                            pv1 = pBNS->vert + fst_cg_vertex + j;
                            break;
                        }
                    }
                    /* index of the pTCGroups->pTCG[] */
                    if ( pv1 ) {
                        k1 = j + fst_cg_group;
                        if ( pTCGroups->pTCG[k1].type != pCN[i1].v.type ||
                             pTCGroups->pTCG[k1].ord_num ) {
                            return RI_ERR_PROGR;
                        }
                    }
                } else
                if ( IS_BNS_VT_M_GR( pCN[i1].v.type ) ) {
                    k1 = pTCGroups->nGroup[TCG_MeFlower0];
                    if ( k1 < 0 ||
                         pTCGroups->pTCG[k1].type != pCN[i1].v.type  ||
                             pTCGroups->pTCG[k1].ord_num ||
                             !pSrm->bMetalAddFlower ) {
                            return RI_ERR_PROGR;
                    }
                    pv1 = pBNS->vert + pTCGroups->pTCG[k1].nVertexNumber;
                } else
                if ( IS_BNS_VT_CHRG_STRUCT(pCN[i1].v.type) ) {
                    /* aux vertex */
                    pv1 = vert_ficpoint_base + i1;
                }
                if ( !pv1 ) {
                    return BNS_BOND_ERR;
                }

                /* connect pairs of vertices with new edges */
                for ( k = 0; k < MAX_CN_VAL && (i2=pCN[i1].e[k].neigh); k ++ ) {
                    pv2 = NULL;
                    k2  = -1;
                    i2 --; /* neighbor */
                    /* find vertex cooresponding to i2 */
                    if ( pCN[i2].v.type & BNS_VERT_TYPE_ATOM ) {
                        pv2 = pBNS->vert+c_point;
                        cap  = !bMetalAtoms? pCN[i2].v.cap : pCN[i2].v.cap? pSrm->nMetalMaxCharge_D : 0;
                        flow = !bMetalAtoms? pCN[i2].v.flow : pCN[i2].v.flow? pSrm->nMetalMaxCharge_D : 0;
                        /* add atom vertex st_cap and st_flow; this normally should not happen */
                        AddStCapFlow( pv2, &tot_st_flow, &tot_st_cap, cap, flow );
                    } else
                    if ( IS_BNS_VT_C_GR(pCN[i2].v.type) ) {
                        /* find c-group vertex by looking for its type */
                        for( j = 0; j < num_cg; j ++ ) {
                            if ( pCN[i2].v.type == pBNS->vert[fst_cg_vertex + j].type ) {
                                pv2 = pBNS->vert + fst_cg_vertex + j;
                                break;
                            }
                        }
                        if ( pv2 ) {
                            k2 = j + fst_cg_group;
                            if ( pTCGroups->pTCG[k2].type != pCN[i2].v.type ||
                                 pTCGroups->pTCG[k2].ord_num ) {
                                return RI_ERR_PROGR;
                            }
                        }
                    } else
                    if ( IS_BNS_VT_M_GR( pCN[i2].v.type ) ) {
                        k2 = pTCGroups->nGroup[TCG_MeFlower0];
                        if ( k2 < 0 ||
                             pTCGroups->pTCG[k2].type != pCN[i2].v.type  ||
                                 pTCGroups->pTCG[k2].ord_num ||
                                 !pSrm->bMetalAddFlower ) {
                                return RI_ERR_PROGR;
                        }
                        pv2 = pBNS->vert + pTCGroups->pTCG[k2].nVertexNumber;
                    } else
                    if ( IS_BNS_VT_CHRG_STRUCT(pCN[i2].v.type) ){
                        pv2 = vert_ficpoint_base + i2;
                    }
                    
                    /* connect pv1 and pv2 */
                    if ( !pv1 || !pv2 || pv1 == pv2 ) {
                        return BNS_BOND_ERR;
                    }
                    j1 = pv1 - pBNS->vert;
                    j2 = pv2 - pBNS->vert;
                    /* create a new edge connecting pv1 and pv2 */
                    edge = pBNS->edge + cur_num_edges;
                    if ( (IS_BNS_VT_M_GR( pCN[i1].v.type ) && IS_BNS_VT_ATOM( pCN[i2].v.type )) ||
                         (IS_BNS_VT_M_GR( pCN[i2].v.type ) && IS_BNS_VT_ATOM( pCN[i1].v.type )) ) {
                        /* at[c_point] is a metal or is treated as a metal; connect it to M-group */
                        /* metal - M-group (i.e. Metal-Flower) edge */
                        int nStCap, nStFlow;
                        bNeedsFlower = AtomStcapStflow( at, pVA, pSrm, c_point, &nStCap, &nStFlow, &edge->cap, &edge->flow );
                        /* GetAtomToMCGroupInitEdgeCapFlow( &edge->cap, &edge->flow, pSrm, at, pVA, c_point ); */
                        if ( !bNeedsFlower ) {
                            return RI_ERR_PROGR;
                        }
                        pVA[c_point].nMetalGroupEdge = cur_num_edges + 1;
                        /* pBNS->vert[c_point].st_edge.cap  += edge->flow;*/ /* where was this done ???*/
                        pBNS->vert[c_point].st_edge.flow += edge->flow;
                        pBNS->vert[c_point].st_edge.cap  += edge->flow + pVA[c_point].cInitFreeValences;
                        pBNS->vert[c_point].st_edge.flow0 = pBNS->vert[c_point].st_edge.flow;
                        pBNS->vert[c_point].st_edge.cap0  = pBNS->vert[c_point].st_edge.cap;
                        tot_st_flow += edge->flow;
                        tot_st_cap  += edge->flow + pVA[c_point].cInitFreeValences;
                    } else {
                        edge->cap     = !bMetalAtoms? pCN[i1].e[k].cap : pCN[i1].e[k].cap? pSrm->nMetalMaxCharge_D : 0;
                        edge->flow    = !bMetalAtoms? pCN[i1].e[k].flow : pCN[i1].e[k].flow? pSrm->nMetalMaxCharge_D : 0;
                    }
                    edge->forbidden = pCN[i1].e[k].bForbiddenEdge? BNS_EDGE_FORBIDDEN_MASK : 0;
                    /* c-group incoming edges cap and flow needed in ConnectSuperCGroup() */
                    /*
                    if ( k1 >= 0 ) {
                        pTCGroups->pTCG[k1].edges_cap  += pCN[i1].e[k].cap;
                        pTCGroups->pTCG[k1].edges_flow += pCN[i1].e[k].flow;
                    }
                    if ( k2 >= 0 ) {
                        pTCGroups->pTCG[k2].edges_cap  += pCN[i1].e[k].cap;
                        pTCGroups->pTCG[k2].edges_flow += pCN[i1].e[k].flow;
                    }
                    */
                    edge->pass      = 0;
#if ( RESET_EDGE_FORBIDDEN_MASK == 1 )
                    edge->forbidden &= pBNS->edge_forbidden_mask;
#endif
                    /* check edge overflow */
                    if ( pv1->num_adj_edges >= pv1->max_adj_edges ||
                         pv2->num_adj_edges >= pv2->max_adj_edges ||
                         cur_num_edges      >= pBNS->max_edges     ) {
                        return BNS_VERT_EDGE_OVFL;
                    }
                    
                    /* connect edge to the incident vertices and increment the counters of neighbors and edges */
                    ret = ConnectTwoVertices( pv1, pv2, edge, pBNS, 0 );
                    if ( IS_BNS_ERROR( ret ) ) {
                        return ret;
                    }
                    edge->cap0  = edge->cap;
                    edge->flow0 = edge->flow;
                    /* save the edge index */
                    type = IS_BNS_VT_C_GR(pv1->type)? pv1->type :
                           IS_BNS_VT_C_GR(pv2->type)? pv2->type : 0;
                    if ( type ) {
                        /* the edge connects to a c-group */
                        if ( type & BNS_VERT_TYPE_C_NEGATIVE ) {
                            pVA[c_point].nCMinusGroupEdge = cur_num_edges+1;
                        } else {
                            pVA[c_point].nCPlusGroupEdge  = cur_num_edges+1;
                        }
                    }
                    cur_num_edges ++; /* end of new edge creation */
                }
            }
        }
        /*************************************************************/
        /* pass 2:                                                   */
        /* adjust bond cap, flow from the final atom st_cap, st_flow */
        /*************************************************************/
        for ( c_point = 0; c_point < num_atoms; c_point ++ ) {
            int st_cap, st_cap2, max_edge_flow;
            pv1 = pBNS->vert + c_point;  /* atom vertex */
            st_cap = pv1->st_edge.cap;
            for ( k = 0; k < pv1->num_adj_edges; k ++ ) {
                edge = pBNS->edge + pv1->iedge[k];      /* incident edge */
                c_neigh = edge->neighbor12 ^ c_point;   /* adjacent vertex */
                pv2 = pBNS->vert + c_neigh;
                if ( c_neigh > c_point || !(pv2->type & BNS_VERT_TYPE_ATOM) ) {
                    continue;
                }
                /* adjacent vertex is an atom; the edge is a bond; process each bond only once */
                st_cap2 = pv2->st_edge.cap;
                /* the edge flow <= min( incident atom st_caps) */
                max_edge_flow = inchi_min( st_cap, st_cap2 );
                /* bond order <= triple bond (flow=2) */
                if ( pSrm->bMetalAddFlower && !pSrm->nMetalMinBondOrder &&
                    ((pVA[c_point].cMetal && pVA[c_point].cNumBondsToMetal) ||
                     (pVA[c_neigh].cMetal && pVA[c_neigh].cNumBondsToMetal)) ) {
                    max_edge_flow = inchi_min( max_edge_flow, MAX_BOND_EDGE_CAP+1 );
                } else {
                    max_edge_flow = inchi_min( max_edge_flow, MAX_BOND_EDGE_CAP );
                }
                if ( at[c_point].bond_type[k] == BOND_TYPE_SINGLE ) {
                    /* the bond has not been changed due to stereo */
                    edge->cap = edge->cap0 = max_edge_flow;
                }
            }
        }
        /***********************************************************/
        /**************                                 ************/
        /************** connect M-flower with new edges ************/
        /**************                                 ************/
        /***********************************************************/
        ret = ConnectMetalFlower(&cur_num_vertices, &cur_num_edges, &tot_st_cap, &tot_st_flow, pSrm, pBNS, pTCGroups);
        if ( ret < 0 ) {
            goto exit_function;
        }
        /***********************************************************/
        /**************                                 ************/
        /************** add additional vertices & edges ************/
        /************** to connect c-groups             ************/
        /**************                                 ************/
        /***********************************************************/
        /* (+) supergroup, Y-connection */
        k = 0;
        nAddGroups[k ++] = TCG_Plus0;
        nAddGroups[k ++] = TCG_Plus_C0;
        nAddGroups[k ++] = TCG_Plus_M0;
        ret1 = ConnectSuperCGroup( TCG_Plus, nAddGroups, k, &cur_num_vertices, &cur_num_edges,
                              &tot_st_cap, &tot_st_flow, pBNS, pTCGroups );
        /* (-) supergroup, Y-connection */
        k = 0;
        nAddGroups[k ++] = TCG_Minus0;
        nAddGroups[k ++] = TCG_Minus_C0;
        nAddGroups[k ++] = TCG_Minus_M0;
        ret2 = ConnectSuperCGroup( TCG_Minus, nAddGroups, k, &cur_num_vertices, &cur_num_edges,
                              &tot_st_cap, &tot_st_flow, pBNS, pTCGroups );
        /******** connect (+) and (-) ***************/
        k = 0;
        nAddGroups[k ++] = TCG_Plus;
        nAddGroups[k ++] = TCG_Minus;
        ret3 = ConnectSuperCGroup( -1, nAddGroups, k, &cur_num_vertices, &cur_num_edges,
                              &tot_st_cap, &tot_st_flow, pBNS, pTCGroups );

        /* Take care of the full charge */
        cg_charge = pTCGroups->total_charge - pTCGroups->tgroup_charge - pTCGroups->charge_on_atoms;
        ret = 1;
        if ( ret3 > 0 ) {
            /* (+) and (-) or at least one of them have been connected */
            int nVertPlusMinus = cur_num_vertices - 1;
            BNS_VERTEX *pVertPlusMinus = pBNS->vert + nVertPlusMinus;
            BNS_VERTEX *pVertPlus = NULL, *pVertMinus = NULL, *pVert=NULL;
            BNS_EDGE   *pEdgePlus = NULL, *pEdgeMinus = NULL, *pEdge=NULL;
            n = pTCGroups->nGroup[TCG_Plus] >= 0;   /* (+)-supergroup exists */
            k = pTCGroups->nGroup[TCG_Minus] >= 0;  /* (-)-supergroup exists */
            if ( pVertPlusMinus->num_adj_edges == 2 && k+n==2 ) {
                pEdgePlus      = pBNS->edge + pVertPlusMinus->iedge[0];  /* TCG_Plus was the 1st */
                pEdgeMinus     = pBNS->edge + pVertPlusMinus->iedge[1];  /* TCG_Minus was the 2nd */
            } else
            if ( pVertPlusMinus->num_adj_edges == 1 && k+n==1 ) {
                if ( pTCGroups->nGroup[TCG_Plus] >= 0 ) {
                    pEdgePlus      = pBNS->edge + pVertPlusMinus->iedge[0];
                } else
                if ( pTCGroups->nGroup[TCG_Minus] >= 0 ) {
                    pEdgeMinus     = pBNS->edge + pVertPlusMinus->iedge[0];
                }
            } else
            if ( k+n ) {
                /* program error check */
                ret = BNS_BOND_ERR;
                goto exit_function;
            }
            if ( pEdgePlus ) {
                pVertPlus  = pBNS->vert + (pEdgePlus->neighbor12 ^ nVertPlusMinus);
            }
            if ( pEdgeMinus ) {
                pVertMinus = pBNS->vert + (pEdgeMinus->neighbor12 ^ nVertPlusMinus);
            }
            pVert = pVertPlus? pVertPlus : pVertMinus? pVertMinus : NULL;
            pEdge = pEdgePlus? pEdgePlus : pEdgeMinus? pEdgeMinus : NULL;
            if ( pEdgeMinus ) {
                pTCGroups->nEdgeMinus = pEdgeMinus - pBNS->edge;
            }
            if ( pEdgePlus ) {
                pTCGroups->nEdgePlus = pEdgePlus - pBNS->edge;
            }
            if ( pEdge ) {
                pTCGroups->nEdge4charge = pEdge - pBNS->edge;
            }
            /* set total charge */
            if ( pVert && pEdge ) {
                /* do not check rescaps for now */
                if ( cg_charge > 0 ) {
                    pVertPlusMinus->st_edge.cap += cg_charge;
                    tot_st_cap                  += cg_charge;
                    pVertPlusMinus->st_edge.cap0 = pVertPlusMinus->st_edge.cap;
                }
                if ( cg_charge < 0 ) {
                    pVert->st_edge.cap -= cg_charge;
                    tot_st_cap         -= cg_charge;
                    pVert->st_edge.cap0 = pVert->st_edge.cap;

                    if ( pEdge->cap - pEdge->flow + cg_charge < 0 ) {
                        /* 2006-02-06: increase edge capacity to avoid clogging */
                        pEdge->cap = pEdge->flow - cg_charge;
                    }
                }
                pTCGroups->added_charge = cg_charge;
                
            }
            if ( !cg_charge || (pVert && pEdge) ) {
                ret = 2;
            }
        }

        AddRadicalToMetal( &tot_st_cap, &tot_st_flow, pSrm, pBNS, pTCGroups );


        pBNS->num_edges     = cur_num_edges;
        pBNS->num_vertices  = cur_num_vertices;
        pBNS->num_c_groups  = num_cg;
        pBNS->tot_st_cap    = tot_st_cap;
        pBNS->tot_st_flow   = tot_st_flow;

    }
exit_function:
    return ret;
}
/********************************************************************************/
int nNumEdgesToCnVertex( MY_CONST C_NODE *pCN, int len, int v )
{
    int i, j, n, num_edges, v1 = v+1;
    for ( i = 0, num_edges = 0; i < len; i ++ ) {
        for ( j = 0; j < MAX_CN_VAL && (n = pCN[i].e[j].neigh); j ++ ) {
            num_edges += ( i == v || n == v1 );
        }
    }
    return num_edges;
}

/*********************************************************************************
BN_STRUCT* AllocateAndInitTCGBnStruct( StrFromINChI *pStruct, VAL_AT *pVA,
                                       ALL_TC_GROUPS *pTCGroups,
                                       int nMaxAddAtoms, int nMaxAddEdges,
                                       int max_altp, int *pNum_changed_bonds )
allocate BN_STRUCT that has:

  pBNS->max_vertices = pTCGroups->nVertices + nMaxAddAtoms   
  pBNS->max_edges    = pTCGroups->nEdges +
                       pBNS->max_vertices * (nMaxAddEdges + NUM_KINDS_OF_GROUPS)
  pBNS->max_iedges   = 2*pBNS->max_edges + pTCGroups->nAddIedges

  pBNS->len_alt_path = pBNS->max_vertices + iALTP_HDR_LEN + 1 +
                       max( pBNS->max_vertices/2, 16 )
  pBNS->max_altp     = max_altp
  
other members:

  pBNS->num_atoms       = num_atoms;
  pBNS->num_bonds       = num_bonds;
  pBNS->num_added_atoms = 0;           
  pBNS->num_t_groups    = 0;           
  pBNS->num_c_groups    = 0;           
  pBNS->nMaxAddAtoms    = nMaxAddAtoms;
  pBNS->nMaxAddEdges    = nMaxAddEdges;

atom vertices and bond edges:
  --- vertex(atom) ---
  st_cap  = (at[].chem_bonds_valence - at[].valence) + pVA[].cInitFreeValences
  st_flow = SUM{bond_orders; ALT_BOND counted as SINGLE} - at[].valence
  --- edge(bond) ---
  flow    = bond_order - 1; for ALT_BOND flow = 0
  cap     = min(min(st_cap of neighbors),2); for ALT_BOND cap = 1
  max number of edges per atom = number of bonds +
                             number of edges to ChargeStruct +
                             1 (if atom is a tautomeric endpoint) +
                             nMaxAddEdges
  --- NOTE ---
  Here are not included nDelta(dots) from ChargeStruct and flow to ChargeStruct

 *********************************************************************************/
BN_STRUCT* AllocateAndInitTCGBnStruct( StrFromINChI *pStruct, VAL_AT *pVA,
                                       ALL_TC_GROUPS *pTCGroups,
                                       int nMaxAddAtoms, int nMaxAddEdges,
                                       int max_altp, int *pNum_changed_bonds )
{
    inp_ATOM    *at        = pStruct->at;
    int          num_atoms = pStruct->num_atoms;
    ICHICONST SRM *pSrm    = pStruct->pSrm;
    
    BN_STRUCT   *pBNS      = NULL;
    BNS_VERTEX  *vert;
    BNS_IEDGE   *iedge;

    int    neigh, num_changed_bonds=0;
    U_CHAR bond_type, bond_mark;
    int bNeedsFlower1, bNeedsFlower2, min_order;

    int i, j, k, m, n_edges, num_bonds, num_edges;
    int f1, f2, c1, c2, edge_cap, edge_flow, st_cap, st_flow, flag_alt_bond;
    int tot_st_cap, tot_st_flow;
    int max_tg, max_edges, max_vertices, len_alt_path, max_iedges, num_iedges, num_altp;

    /* count vertices */
    max_tg = pTCGroups->num_tgroups;
    /* +1 for a super-tautomeric group */
    /* max_vertices = num_atoms + nMaxAddAtoms + max_tg + 1; */
    max_vertices = pTCGroups->nVertices + nMaxAddAtoms;
    
    /* count edges */
    num_changed_bonds = 0;
    num_bonds = pTCGroups->num_bonds;

    /* each atom has enough edges to belong to a tautomeric group + nMaxAddEdges */
    /* number of atoms is large enough to accommodate max. possible number of t-groups + nMaxAddAtoms */
    /* max_altp cannot be larger than BN_MAX_ALTP = 16 */
    num_edges    = pTCGroups->nEdges;
    /* +max_tg for edges between t-groups and super-tautomeric group */
    max_edges    = num_edges + (nMaxAddEdges + NUM_KINDS_OF_GROUPS)*max_vertices;
    max_iedges   = 2*max_edges + pTCGroups->nAddIedges;
    len_alt_path = max_vertices+iALTP_HDR_LEN + 1; /* may overflow if an edge is traversed in 2 directions */
    len_alt_path += inchi_max( max_vertices/2, 16 ); /* to avoid the overflow */

    if ( !( pBNS           = (BN_STRUCT   *)inchi_calloc( 1,           sizeof(BN_STRUCT)) )  ||
         !( pBNS->edge     = (BNS_EDGE    *)inchi_calloc( max_edges,   sizeof(BNS_EDGE)) )   ||
         !( pBNS->vert     = (BNS_VERTEX  *)inchi_calloc( max_vertices,sizeof(BNS_VERTEX)) ) ||
         !( pBNS->iedge    = (BNS_IEDGE   *)inchi_calloc( max_iedges,  sizeof(BNS_IEDGE)) ) ) { 
        return DeAllocateBnStruct( pBNS );
    }
    /* alt path init (standard spell) */
    for ( num_altp = 0; num_altp < max_altp && num_altp < BN_MAX_ALTP; num_altp ++ ) {
        if ( !( pBNS->altp[num_altp] = (BNS_ALT_PATH*)inchi_calloc( len_alt_path,sizeof(BNS_ALT_PATH))) ) {
            return DeAllocateBnStruct( pBNS );
        }
        ALTP_ALLOCATED_LEN(pBNS->altp[num_altp]) = len_alt_path;
        pBNS->len_alt_path                 = len_alt_path;  /* ??? duplication ??? */
        /* re-init */
        ALTP_DELTA(pBNS->altp[num_altp])         = 0;
        ALTP_START_ATOM(pBNS->altp[num_altp])    = NO_VERTEX;
        ALTP_END_ATOM(pBNS->altp[num_altp])      = NO_VERTEX;
        ALTP_PATH_LEN(pBNS->altp[num_altp])      = 0;
    }
    pBNS->alt_path = NULL;
    pBNS->num_altp = 0;
    pBNS->max_altp = num_altp;


    /* fill vertices (no connectivity) */
    iedge = pBNS->iedge;
    num_iedges = 0;
    tot_st_cap = tot_st_flow = 0;
    for ( i = 0; i < num_atoms; i ++ ) {
        /* count edges incident to pBNS->vert[i] */
        k = at[i].valence + (at[i].endpoint != 0) +  (nMaxAddEdges /*+ NUM_KINDS_OF_GROUPS*/);
        if ( (j = pVA[i].cnListIndex-1) >= 0 ) {
            /* add number of neighbors in the ChargeStruct */
            k += nNumEdgesToCnVertex( cnList[j].pCN, cnList[j].len, 0 );
        }
        /* set max number of edges for the vertex */
        pBNS->vert[i].max_adj_edges = k;
        pBNS->vert[i].iedge = iedge;
        iedge += k;
        /* add atom vertex cap */
        st_cap  = 0;
        st_flow = 0;
        bNeedsFlower1 = AtomStcapStflow( at, pVA, pSrm, i, &c1, &f1, NULL, NULL );
        /* pVA[i].cNumBondsToMetal = bNeedsFlower1; */
        /* GetAtomStCapFlow( at, pVA, pSrm, i, &c1, &f1 ); */
        st_cap += c1;
        st_cap += bNeedsFlower1? 0 : pVA[i].cInitFreeValences;
        pBNS->vert[i].st_edge.cap = st_cap; /* the 1st time st_cap is set */
        pBNS->vert[i].st_edge.cap0 = pBNS->vert[i].st_edge.cap;
        tot_st_cap += st_cap;
    }
    num_iedges = iedge - pBNS->iedge;
    if ( max_iedges - num_iedges < (nMaxAddEdges + NUM_KINDS_OF_GROUPS)*max_vertices ) {
        return DeAllocateBnStruct( pBNS );
    }
    
    pBNS->num_atoms       = num_atoms;      /* number of real atoms */
    pBNS->num_added_atoms = 0;
    pBNS->num_t_groups    = 0;              /* number of added t-groups */
    pBNS->num_c_groups    = 0;
    pBNS->nMaxAddAtoms    = nMaxAddAtoms;
    pBNS->nMaxAddEdges    = nMaxAddEdges;

    pBNS->num_vertices    = num_atoms;      /* current number of vertices, in general a sum of
                                               pBNS->num_atoms
                                               pBNS->num_t_groups
                                               number of c-groups
                                               number of auxiliary vertices
                                               pBNS->num_added_atoms
                                            */
    pBNS->max_vertices    = max_vertices;
    

    pBNS->num_bonds       = num_bonds;      /* number of real edges (bonds) */
    pBNS->max_edges       = max_edges;
    pBNS->max_iedges      = max_iedges;


    /* 
       To remove t-groups and added atoms:

        for ( i = 0; i < pBNS->num_atoms; i ++ ) {
            for ( j = pBNS->vert[i].num_adj_edges-1; 0 <= j; j -- ) {
                k = pBNS->edge[pBNS->vert[i].iedge[j]].neighbor12 ^ i;
                if ( pBNS->vert[k].type & BNS_VERT_TYPE_ATOM ) {
                    pBNS->vert[i].num_adj_edges = j+1;
                    break;
                }
            }
        }

       pBNS->num_vertices    = pBNS->num_atoms;
       pBNS->num_edges       = pBNS->num_bonds;
       pBNS->num_added_atoms = 0;
       pBNS->num_t_groups    = 0;
       pBNS->num_added_edges = 0;

       ALTP_DELTA(pBNS->alt_path)      = 0;
       ALTP_START_ATOM(pBNS->alt_path) = NO_VERTEX;
       ALTP_END_ATOM(pBNS->alt_path)   = NO_VERTEX;
       ALTP_PATH_LEN(pBNS->alt_path)   = 0;

    */


    /* add and fill edges and connectivity */
    for ( i = 0, n_edges = 0; i < num_atoms; i ++ ) {
        vert    = pBNS->vert + i; /* pointer to the ith vertex */
        st_cap  = 0;
        st_flow = 0;
        flag_alt_bond = 0;
        for ( j = 0; j < at[i].valence; j ++ ) {
            neigh = at[i].neighbor[j];
            /* find this bond at the neighbor */
            for ( k = 0; k < at[neigh].valence; k ++ ) {
                if ( at[neigh].neighbor[k] == i ) {
                    break;
                }
            }
            bond_type = (at[i].bond_type[j] & BOND_TYPE_MASK);
            bond_mark = (at[i].bond_type[j] & ~BOND_TYPE_MASK);
            if ( bond_type != BOND_SINGLE && bond_type != BOND_DOUBLE && bond_type != BOND_TRIPLE ) {
                /* make unknown bonds single */
                bond_type = BOND_SINGLE;
                at[i].bond_type[j] = bond_mark | bond_type;
                num_changed_bonds ++;
            }
            if ( neigh > i ) {
                /* this is the first time we encounter this bond */
                bNeedsFlower1 = AtomStcapStflow( at, pVA, pSrm, i, &c1, &f1, NULL, NULL );
                /* GetAtomStCapFlow( at, pVA, pSrm, i, &c1, &f1 ); */
                c1 += bNeedsFlower1? 0 : pVA[i].cInitFreeValences;  /* elevate cap to the lowest valence in ChargeStruct */
                bNeedsFlower2 = AtomStcapStflow( at, pVA, pSrm, neigh, &c2, &f2, NULL, NULL );
                /* GetAtomStCapFlow( at, pVA, pSrm, neigh, &c2, &f2 ); */
                c2 += bNeedsFlower2? 0 : pVA[neigh].cInitFreeValences; /* elevate cap to the lowest valence in ChargeStruct */
                /* at this point -O would have st_cap=st_flow=0 because the lowest valence=1 for charge=-1 */
                /* however, if -O belongs to a t-group its cap would be 1, flow = 0 */
                /*f1 = MAX_AT_FLOW(at[i]);*/
                /*f2 = MAX_AT_FLOW(at[neigh]);*/
                edge_flow = BondFlowMaxcapMinorder( at, pVA, pSrm, i, j, &edge_cap, &min_order, NULL);

                pBNS->edge[n_edges].neighbor1    = (AT_NUMB)i;
                pBNS->edge[n_edges].neighbor12   = (AT_NUMB)(i ^ neigh);
                pBNS->edge[n_edges].flow =
                pBNS->edge[n_edges].flow0        = edge_flow;
                pBNS->edge[n_edges].cap  =
                pBNS->edge[n_edges].cap0         = edge_cap;
                pBNS->edge[n_edges].neigh_ord[0] = j;  /* iedge to neigh index at vertex[i],     i < neigh */
                pBNS->edge[n_edges].neigh_ord[1] = k;  /* iedge to i     index at vertex[neigh], i < neigh */
                pBNS->edge[n_edges].pass         = 0;
                pBNS->edge[n_edges].forbidden    = 0; /* may be forbidden if edge_flow = 1: stereogenic fixed double bond */
                if ( bond_type == BOND_TYPE_DOUBLE ) {
                    /* forbid changing stereogenic double bonds */
                    for ( m = 0; m < MAX_NUM_STEREO_BONDS && at[i].sb_parity[m]; m ++ ) {
                        if ( at[i].sb_ord[m] == j ) {
                            pBNS->edge[n_edges].forbidden |= BNS_EDGE_FORBIDDEN_MASK;
                            break;
                        }
                    }
                }
                vert->iedge[j] = pBNS->vert[neigh].iedge[k] = n_edges ++; /* same iedge index as neighbor index in at[] */
            } else {
                /* this is the second time we encounter this bond. It was stored at */
                int  iedge2 = pBNS->vert[neigh].iedge[k];
                edge_cap    = pBNS->edge[iedge2].cap;
                edge_flow   = pBNS->edge[iedge2].flow;
            }
            st_flow += edge_flow;
            /*
            st_cap  += edge_cap;
            */
        }
        vert->num_adj_edges = j;
        /*
        vert->st_edge.cap   =
        vert->st_edge.cap0  = st_cap;
        */
        vert->st_edge.flow  =
        vert->st_edge.flow0 = st_flow;
        vert->type          = BNS_VERT_TYPE_ATOM;
        /*
        tot_st_cap  += vert->st_edge.cap;
        */
        tot_st_flow += vert->st_edge.flow;
    }
    *pNum_changed_bonds = num_changed_bonds/2;

    pBNS->num_edges       = n_edges;   /* number of edges */
    pBNS->num_iedges      = num_iedges;
    pBNS->num_added_edges = 0;

    pBNS->tot_st_cap  = tot_st_cap;
    pBNS->tot_st_flow = tot_st_flow;

/* exit_function: */

    return pBNS;
}
/******************************************************************************************************/
void IncrZeroBondsAndClearEndpts(inp_ATOM *at, int num_at, int iComponent )
{
    int i, j;
    for ( i = 0; i < num_at; i ++ ) {
        at[i].endpoint = 0;
        at[i].component = iComponent;
        for ( j = 0; j < at[i].valence; j ++ ) {
            if ( !at[i].bond_type[j] ) {
                at[i].bond_type[j]        = BOND_TYPE_SINGLE;
                at[i].chem_bonds_valence += BOND_TYPE_SINGLE;
            }
        }
    }
}
void IncrZeroBonds(inp_ATOM *at, int num_at, int iComponent )
{
    int i, j;
    for ( i = 0; i < num_at; i ++ ) {
        at[i].component = iComponent;
        for ( j = 0; j < at[i].valence; j ++ ) {
            if ( !at[i].bond_type[j] ) {
                at[i].bond_type[j]        = BOND_TYPE_SINGLE;
                at[i].chem_bonds_valence += BOND_TYPE_SINGLE;
            }
        }
    }
}
void ClearEndpts(inp_ATOM *at, int num_at )
{
    int i;
    for ( i = 0; i < num_at; i ++ ) {
        at[i].endpoint = 0;
    }
}


/******************************************************************************************************/
#define ANY_VERT_TYPE(X) (((X) & (BNS_VERT_TYPE_ATOM | BNS_VERT_TYPE_TGROUP | BNS_VERT_TYPE_C_GROUP)) && \
                           !((X) & (BNS_VERT_TYPE_SUPER_TGROUP)))
#define GRP_VERT_TYPE(X) (((X) & (BNS_VERT_TYPE_TGROUP | BNS_VERT_TYPE_C_GROUP)) && \
                           !((X) & (BNS_VERT_TYPE_SUPER_TGROUP)))
typedef struct tagVertexFlow {
    int        type;
    Vertex     v;
    EdgeIndex  e_In;
    EdgeIndex  e_Out;
    EdgeFlow   delta_In;
    EdgeFlow   delta_Out;
    Vertex     bUsed; /* indicates the charge edge belongs to already processed atom */
} VF;
#define NUM_VF 3
#define VF_USED_IN  1
#define VF_USED_OUT 2
#define VF_USED_ALL (VF_USED_IN | VF_USED_OUT)

int GetDeltaChargeFromVF( BN_STRUCT *pBNS, VAL_AT *pVA, VF *vf );
/******************************************************************************************************/
int GetDeltaChargeFromVF( BN_STRUCT *pBNS, VAL_AT *pVA, VF *vf )
{
    int i, v = NO_VERTEX;
    int ieIn1  = (!(vf->bUsed & VF_USED_IN)  && vf->e_In  >= 0 && vf->delta_In )? vf->e_In+1  : NO_VERTEX;
    int ieOut1 = (!(vf->bUsed & VF_USED_OUT) && vf->e_Out >= 0 && vf->delta_Out)? vf->e_Out+1 : NO_VERTEX;
    int nInitCharge, nPlusFlow, nMinusFlow, nDeltaCharge, nNumDeltaCharge, eCPlus, eCMinus;

    if ( !(vf->type & BNS_VERT_TYPE_C_GROUP) ||
          (vf->type & BNS_VERT_TYPE_SUPER_TGROUP) ||
          (ieIn1 == NO_VERTEX && ieOut1 == NO_VERTEX ) ) {
        return 0;
    }
    if ( vf->type & BNS_VERT_TYPE_C_NEGATIVE ) {
        /* negative charge edge */
        for ( i = 0; i < pBNS->num_atoms; i ++ ) {
            if ( pVA[i].nCMinusGroupEdge == ieIn1 || pVA[i].nCMinusGroupEdge == ieOut1 ) {
                v = i;
                break;
            }
        }
    } else {
        /* positive charge edge */
        for ( i = 0; i < pBNS->num_atoms; i ++ ) {
            if ( pVA[i].nCPlusGroupEdge == ieIn1 || pVA[i].nCPlusGroupEdge == ieOut1 ) {
                v = i;
                break;
            }
        }
    }
    if ( v == NO_VERTEX )
        return 0;
    
    nInitCharge = pVA[v].cInitCharge;
    nPlusFlow = nMinusFlow = 0;
    nNumDeltaCharge = 0;

    if ( (eCPlus = pVA[v].nCPlusGroupEdge-1) >= 0 ) {
        nPlusFlow = pBNS->edge[eCPlus].cap
                  - pBNS->edge[eCPlus].flow;
    }
    if ( (eCMinus = pVA[v].nCMinusGroupEdge-1) >= 0 ) {
        nMinusFlow = -pBNS->edge[eCMinus].flow;
    }
    nInitCharge += nPlusFlow + nMinusFlow;

    nDeltaCharge = 0;
    
    if ( !(vf[0].bUsed & VF_USED_OUT) ) {
        if ( vf[0].e_Out==eCPlus || vf[0].e_Out==eCMinus ) {
            nDeltaCharge -= vf[0].delta_Out;
            vf[0].bUsed  |= VF_USED_OUT;
        }
    }
    
    if ( !(vf[0].bUsed & VF_USED_IN) ) {
        if ( vf[0].e_In==eCPlus || vf[0].e_In==eCMinus ) {
            nDeltaCharge -= vf[0].delta_In;
            vf[0].bUsed |= VF_USED_IN;
        }
    }
    if ( !nInitCharge && nDeltaCharge ) {
        nNumDeltaCharge ++;
    } else
    if ( nInitCharge && 0 == nInitCharge + nDeltaCharge ) {
        nNumDeltaCharge --;
    }
    return nNumDeltaCharge;
}
/******************************************************************************************************/
int EvaluateChargeChanges( BN_STRUCT *pBNS, VAL_AT *pVA, int *pnDeltaH, int *pnDeltaCharge, int *pnNumVisitedAtoms )
{
    int       pass, i, j, v0, v1, v2, v, ineigh1, /*ineigh2,*/ vLast, n, delta, ret, ie, err = 0;
    BNS_EDGE *edge;
    int       nDeltaH, nDeltaCharge, iPrev, nInitCharge, nPlusFlow, nMinusFlow;
    int       nNumDeltaH = 0;
    int       nNumDeltaCharge = 0;
    int       nNumVisitedAtoms = 0;
    VF   vf[NUM_VF+1];
    
    *pnDeltaH           = 0;
    *pnDeltaCharge      = 0;
    *pnNumVisitedAtoms  = 0;

    for ( pass = pBNS->num_altp-1, ret = 0; 0 <= pass; pass -- ) {
        
        pBNS->alt_path = pBNS->altp[pass];
        v1 = ALTP_START_ATOM(pBNS->alt_path);
        n  = ALTP_PATH_LEN(pBNS->alt_path);
        delta  = ALTP_DELTA(pBNS->alt_path);
        vLast = ALTP_END_ATOM(pBNS->alt_path);
        v0 = v2 = NO_VERTEX;

        memset( vf, 0, sizeof(vf) );
        for ( i = 0; i < (int)(sizeof(vf)/sizeof(vf[0])); i ++ ) {
            vf[i].v      = NO_VERTEX; /* = -2 */
            vf[i].e_In   = NO_VERTEX;
            vf[i].e_Out  = NO_VERTEX;
        }
        iPrev = 0;
        /* add to the queue */
        if ( ANY_VERT_TYPE(pBNS->vert[v1].type) ) {
            if (pBNS->vert[v1].type & BNS_VERT_TYPE_ATOM) {
                nNumVisitedAtoms ++;
            }
            vf[2].type = pBNS->vert[v1].type;
            vf[2].v    = v1;
            iPrev = 2;
        }

        nNumDeltaH = 0;      
        nNumDeltaCharge = 0; 
        nNumVisitedAtoms = 0;

        for ( i = 0; i < n; i ++, delta = -delta, v0 = v1, v1 = v2 ) {
            ineigh1 = ALTP_THIS_ATOM_NEIGHBOR(pBNS->alt_path, i);  /* v1->v2 neighbor */
            /*ineigh2 = ALTP_NEXT_ATOM_NEIGHBOR(pBNS->alt_path, i);*/  /* v2->v1 neighbor */
            edge = pBNS->edge + (ie=pBNS->vert[v1].iedge[ineigh1]);
            /* follow the BN Structure, not the inp_ATOM, to take care of swithching to
               t-groups, c-groups or other fictitious edges/vertices
            */
            
            if ( iPrev ) {
                /* add exit delta and edge */
                vf[2].e_Out     = ie;
                vf[2].delta_Out = delta;
            }
            
            v2 = edge->neighbor12 ^ v1;  /* next vertex */
            if (pBNS->vert[v2].type & BNS_VERT_TYPE_ATOM) {
                nNumVisitedAtoms ++;
            }

            if ( (ANY_VERT_TYPE(pBNS->vert[v2].type) || i == n-1) &&
                 (vf[0].type & BNS_VERT_TYPE_C_GROUP) && vf[0].bUsed != VF_USED_ALL ) {
                /* unused vertex is about to be discarded */
                nNumDeltaCharge += GetDeltaChargeFromVF( pBNS, pVA, &vf[0] );
            }

            if ( ANY_VERT_TYPE(pBNS->vert[v2].type) ) {
                /* shift the queue */
                vf[0] = vf[1];
                vf[1] = vf[2];
                vf[2] = vf[3]; /* make vf[2] empty */
                /* add next vertex */
                vf[2].v        = v2;
                vf[2].type     = pBNS->vert[v2].type;
                vf[2].e_In     = ie;
                vf[2].delta_In = delta;
                iPrev = 2; /* indicates a newly added vertex */
            } else
            if ( i == n-1 ) {
                /* shift the queue */
                vf[0] = vf[1];
                vf[1] = vf[2];
                vf[2] = vf[3]; /* make vf[2] empty */
                iPrev = 1; /* indicates the last vertex */
            } else {
                iPrev = 0; /* no new vertex has been added */
            }
            
            if ( iPrev && (vf[1].type & BNS_VERT_TYPE_ATOM)) {
                /* a new vertex has just been added and  */
                /* an atom is in the middle of the queue */
                EdgeIndex eCPlus, eCMinus;
                v = vf[1].v;
                nInitCharge = pVA[v].cInitCharge;
                nPlusFlow = nMinusFlow = 0;
                if ( (eCPlus = pVA[v].nCPlusGroupEdge-1) >= 0 ) {
                    nPlusFlow = pBNS->edge[eCPlus].cap
                              - pBNS->edge[eCPlus].flow;
                }
                if ( (eCMinus = pVA[v].nCMinusGroupEdge-1) >= 0 ) {
                    nMinusFlow = -pBNS->edge[eCMinus].flow;
                }
                nInitCharge += nPlusFlow + nMinusFlow;

                nDeltaH = nDeltaCharge = 0;
                
                if ( vf[0].type & BNS_VERT_TYPE_TGROUP ) {
                    nDeltaH -= delta;
                } else
                if ( (vf[0].type & BNS_VERT_TYPE_C_GROUP) && !(vf[0].bUsed & VF_USED_OUT) ) {
                    if ( vf[0].e_Out==eCPlus || vf[0].e_Out==eCMinus ) {
                        nDeltaCharge -= vf[0].delta_Out;
                        vf[0].bUsed  |= VF_USED_OUT;
                    }
                }
                
                if ( vf[2].type & BNS_VERT_TYPE_TGROUP ) {
                    nDeltaH += delta;
                } else
                if ( (vf[2].type & BNS_VERT_TYPE_C_GROUP) && !(vf[2].bUsed & VF_USED_IN) ) {
                    if ( vf[2].e_In==eCPlus || vf[2].e_In==eCMinus ) {
                        nDeltaCharge -= vf[2].delta_In;
                        vf[2].bUsed |= VF_USED_IN;
                    }
                }
                if ( !nInitCharge && nDeltaCharge ) {
                    nNumDeltaCharge ++;
                } else
                if ( nInitCharge && 0 == nInitCharge + nDeltaCharge ) {
                    nNumDeltaCharge --;
                }
                
                nNumDeltaH      += abs(nDeltaH);
                /* nNumDeltaCharge += abs(nDeltaCharge); */
                vf[1].bUsed = VF_USED_ALL;
            }
        }
        for ( j = 0; j < 3; j ++ ) {
            nNumDeltaCharge += GetDeltaChargeFromVF( pBNS, pVA, &vf[j] );
        }
       
        *pnDeltaH           += nNumDeltaH;
        *pnDeltaCharge      += nNumDeltaCharge;
        *pnNumVisitedAtoms  += nNumVisitedAtoms;


        if ( v2 != vLast ) {
            err = BNS_PROGRAM_ERR;
        }
    }
    return err? err : ret;
}

/******************************************************************************************************/
int RunBnsTestOnce( BN_STRUCT *pBNS, BN_DATA *pBD, VAL_AT *pVA, Vertex *pvFirst, Vertex *pvLast,
                    int *pPathLen, int *pnDeltaH, int *pnDeltaCharge, int *pnNumVisitedAtoms  )
{
    int bChangeFlow = 0; /* do not change flow */
    int delta, ret, ret2, pass;
    
    ReInitBnStructAltPaths( pBNS );
    pass = 0;
    pBNS->alt_path = pBNS->altp[pass];
    pBNS->num_altp = 0;
    pBNS->bChangeFlow = 0;
    delta=BalancedNetworkSearch ( pBNS, pBD, bChangeFlow );
    if ( delta > 0 ) {
        pBNS->alt_path = pBNS->altp[pass];
        *pvFirst   = ALTP_START_ATOM(pBNS->alt_path);
        *pPathLen  = ALTP_PATH_LEN(pBNS->alt_path);
        *pvLast    = ALTP_END_ATOM(pBNS->alt_path);
        pBNS->num_altp ++;
        ret2 = EvaluateChargeChanges( pBNS, pVA, pnDeltaH, pnDeltaCharge, pnNumVisitedAtoms );
    } else {
        *pvFirst   = NO_VERTEX;
        *pPathLen  = 0;
        *pvLast    = NO_VERTEX;
        ret2       = 0;
    }
    ReInitBnStructAltPaths( pBNS );
    ret = ReInitBnData( pBD );
    return (delta >= 0 && ret > 0 )? -ret : delta;

}
/******************************************************************************************************/
int RunBnsRestoreOnce( BN_STRUCT *pBNS, BN_DATA *pBD, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups )
{
    /* run BNS for the first time */
    int nTotalDelta = 0, ret = 0;
    int nDelta;
    ReInitBnStructAltPaths( pBNS );
    do {
        nDelta = RunBalancedNetworkSearch( pBNS, pBD, BNS_EF_CHNG_FLOW );
        if ( IS_BNS_ERROR(nDelta) ) {
            ret = nDelta;
            goto exit_function;
        }
        nTotalDelta += nDelta;
        ReInitBnStructAltPaths( pBNS );
        ret    = ReInitBnData( pBD );
        if ( ret > 0 ) {
            ret = -ret;
            goto exit_function;
        }
    } while( nDelta > 0 && ret == 0 );
    pBNS->tot_st_flow += 2*nTotalDelta;
    ret = nTotalDelta;
exit_function:
    return ret;
}
/******************************************************************************************************/
int comp_cc_cand( const void *a1, const void *a2 )
{
    const CC_CAND *p1 = (const CC_CAND *) a1;
    const CC_CAND *p2 = (const CC_CAND *) a2;
    int            ret;
    if ( (ret = (int)p2->cMetal - (int)p1->cMetal) )
        return ret; /* metal first */
    if ( (ret = (int)p2->cNumBondsToMetal - (int)p1->cNumBondsToMetal) )
        return ret; /* connected to metal first */
    if ( (ret = (int)p2->cPeriodicRowNumber - (int)p1->cPeriodicRowNumber) )
        return ret; /* heaviest first */
    if ( (ret = (int)p2->num_bonds - (int)p1->num_bonds) )
        return ret; /* more bonds first */
    if ( (ret = (int)p1->chem_valence - (int)p2->chem_valence) )
        return ret; /* less bond order first */
    if ( !p1->cNumValenceElectrons && p2->cNumValenceElectrons )
        return -1; /* no valence electrons first */
    if ( !p2->cNumValenceElectrons && p1->cNumValenceElectrons )
        return -1; /* no valence electrons first */
    if ( (int)p2->cNumValenceElectrons - (int)p1->cNumValenceElectrons )
        return ret; /* more valence electrons first */
    ret = (int)p2->iat - (int)p1->iat; /* greater canon number first */
    return ret;
}
/*****************************************************************************************************
Locate E1=C-E2 where
         e ev  are the edges

       E1 and E2 are atoms that belong to the same t-group
       C         is an atom that does not belong to any t-group
       e         is a forbidden edge
       ev        is not a forbidden edge

  Make changes so that:
       E1(d)-C(d)-E2

       where (d) means doublet radical

*/
/**************************************************************************************************/
int get_pVA_atom_type( VAL_AT *pVA, inp_ATOM *at, int iat, int bond_type )
{
    int type = 0, val;
    if ( pVA[iat].cNumValenceElectrons == 4 ) {
        if ( pVA[iat].cPeriodicRowNumber == 1 ) {
            type |= EL_TYPE_C;
        }
    } else
    if ( pVA[iat].cNumValenceElectrons == 6 ) {
        if ( pVA[iat].cPeriodicRowNumber == 1 ) {
            type |= EL_TYPE_O;
        } else
        if ( pVA[iat].cPeriodicRowNumber < 5 ) {
            type |= EL_TYPE_S;
        }
        if ( bond_type == BOND_TYPE_SINGLE &&
             (type & (EL_TYPE_O | EL_TYPE_S)) &&
             1 == nNoMetalBondsValence(at, iat ) &&
             1 == nNoMetalNumBonds(at, iat) ) {
            type |= EL_TYPE_OSt;
        }
    } else
    if ( pVA[iat].cNumValenceElectrons == 5 ) {
        if ( pVA[iat].cPeriodicRowNumber == 1 ) {
            type |= EL_TYPE_N;
        } else {
            type |= EL_TYPE_P;
        }
    } else
    if ( !is_el_a_metal(pVA[iat].cPeriodicNumber) ) {
        type |= EL_TYPE_X;
    }
    /* check for possibility to be a tautomeric endpoint (that is, be a Mobile H site) */
    val = get_endpoint_valence( at[iat].el_number );
    if ( val && val > at[iat].valence && !at[iat].radical &&
         -1 <= at[iat].charge && at[iat].charge <= 0 &&
         val == at[iat].chem_bonds_valence - at[iat].charge + at[iat].num_H ) {
        type |= EL_TYPE_PT;
    }
    return type;
}

/*************************************************************************************/
int AllocEdgeList( EDGE_LIST *pEdges, int nLen )
{
    switch( nLen ) {
    case EDGE_LIST_FREE:
        if ( NULL != pEdges->pnEdges ) {
            inchi_free( pEdges->pnEdges );
        }
        /* fall through */
    case EDGE_LIST_CLEAR:
        memset( pEdges, 0, sizeof(*pEdges) );
        break;
    default:
        if ( nLen > 0 && nLen != pEdges->num_alloc ) {
            EdgeIndex *tmp_edges = pEdges->pnEdges;
            int        tmp_num   = pEdges->num_edges;
            pEdges->pnEdges = (EdgeIndex *)inchi_calloc( nLen, sizeof(pEdges->pnEdges[0]));
            if ( !pEdges->pnEdges ) {
                return RI_ERR_ALLOC;
            }
            tmp_num = inchi_min( tmp_num, nLen );
            if ( tmp_edges && tmp_num > 0 ) {
                memcpy( pEdges->pnEdges, tmp_edges, tmp_num * sizeof(pEdges->pnEdges[0]) );
                pEdges->num_edges = tmp_num;
            } else {
                pEdges->num_edges = 0;
            }
            if ( tmp_edges ) {
                inchi_free( tmp_edges );
            }
            pEdges->num_alloc = nLen;
            return 0;
        }
        break;
    }
    return 0;
}
/********************************************************************/
int AddToEdgeList( EDGE_LIST *pEdges, int iedge, int nAddLen )
{
    if ( pEdges->num_alloc == pEdges->num_edges ) {
        int ret;
        if ( nAddLen <= 0 ) {
            return RI_ERR_PROGR;
        }
        if ( (ret = AllocEdgeList( pEdges, pEdges->num_alloc + nAddLen )) ) {
            return ret;
        }
    }
    pEdges->pnEdges[pEdges->num_edges ++] = (EdgeIndex)iedge;
    return 0;
}
/********************************************************************/
int RemoveFromEdgeListByIndex( EDGE_LIST *pEdges, int index )
{
    int len;
    if ( 0 <= (len = pEdges->num_edges - index - 1) ) {
        if ( len ) {
            memmove( pEdges->pnEdges+index, pEdges->pnEdges+index+1, len*sizeof(pEdges->pnEdges[0]));
        }
        pEdges->num_edges --;
        pEdges->pnEdges[pEdges->num_edges] = 0;
        return 0;
    }
    return -1;
}
/********************************************************************/
int FindInEdgeList( EDGE_LIST *pEdges, int iedge )
{
    int i;
    EdgeIndex ie = iedge;
    for ( i = pEdges->num_edges-1; 0 <= i; i -- ) {
        if ( ie == pEdges->pnEdges[i] ) {
                return i;
        }
    }
    return -1;
}
/********************************************************************/
int RemoveFromEdgeListByValue( EDGE_LIST *pEdges, int iedge )
{
    int i, ret, n = 0;
    EdgeIndex ie = iedge;
    for ( i = pEdges->num_edges-1; 0 <= i; i -- ) {
        if ( ie == pEdges->pnEdges[i] ) {
            if ( (ret = RemoveFromEdgeListByIndex( pEdges, i )) ) {
                return ret;
            }
            n ++;
        }
    }
    return n;
}
/********************************************************************/
int AllocBfsQueue( BFS_Q *pQ, int num_at, int min_ring_size )
{
    int ret = 0;
    switch( num_at ) {
    case BFS_Q_FREE:
        if ( pQ->q ) {
            pQ->q = QueueDelete( pQ->q );
        }
        if ( pQ->nAtomLevel ) {
            inchi_free( pQ->nAtomLevel );
        }
        if ( pQ->cSource ) {
            inchi_free( pQ->cSource );
        }
        /* fall through */
    case BFS_Q_CLEAR:
        memset( pQ, 0, sizeof( *pQ ) );
        return 0;
    default:
        if ( num_at <= 0 ) {
            ret = RI_ERR_PROGR;
            goto exit_function;
        }
        if ( num_at > pQ->num_at ) {
            if ( pQ->num_at ) {
                AllocBfsQueue( pQ, BFS_Q_FREE, 0 );
            }
            pQ->q          = QueueCreate( num_at+1, sizeof(qInt) );
            pQ->nAtomLevel = (AT_RANK*)inchi_calloc( sizeof(pQ->nAtomLevel[0]), num_at );
            pQ->cSource    = (S_CHAR *)inchi_calloc( sizeof(pQ->cSource[0]), num_at );
            if ( !pQ->q || !pQ->cSource || !pQ->nAtomLevel ) {
                ret = RI_ERR_ALLOC;
                goto exit_function;
            }
            pQ->num_at = num_at;
        }
        pQ->min_ring_size = min_ring_size;
    }
exit_function:
    return ret;
}
/*************************************************************************************/
void RemoveForbiddenEdgeMask( BN_STRUCT *pBNS, EDGE_LIST *pEdges, int forbidden_edge_mask  )
{
    int i, mask = ~forbidden_edge_mask;
    for ( i = 0; i < pEdges->num_edges; i ++ ) {
        pBNS->edge[pEdges->pnEdges[i]].forbidden &= mask;
    }
}
/*************************************************************************************/
void SetForbiddenEdgeMask( BN_STRUCT *pBNS, EDGE_LIST *pEdges, int forbidden_edge_mask  )
{
    int i;
    for ( i = 0; i < pEdges->num_edges; i ++ ) {
        pBNS->edge[pEdges->pnEdges[i]].forbidden |= forbidden_edge_mask;
    }
}
/******************************************************************************************************/
void RemoveForbiddenBondFlowBits( BN_STRUCT *pBNS, int forbidden_edge_mask_int )
{
    BNS_EDGE   *e;
    int         i;
    int         inv_forbidden_edge_mask = ~forbidden_edge_mask_int;
    for ( i = 0, e = pBNS->edge; i < pBNS->num_bonds; i ++, e ++ ) {
        e->forbidden &= inv_forbidden_edge_mask;
    }
}
/******************************************************************************************************
         upper   vc
         edge   /
     v1[i0]---v0
       \     /
        \   /
         \ /
          v1[i1]
          |
          |
          atom
*/
int GetChargeFlowerUpperEdge( BN_STRUCT *pBNS, VAL_AT *pVA, int nChargeEdge )
{
    int ret = NO_VERTEX, i, j, k, i0, i1;
    Vertex  v0, v1[3], vc, v_t, v;
    BNS_EDGE   *pe, *pe1[3], *pe_t;
    BNS_VERTEX *pv0, *pv1[3], *pv_t;

    if ( nChargeEdge < 0 ) {
        goto exit_function;
    }
    pe = pBNS->edge + nChargeEdge;
    vc = pe->neighbor1; /* charge vertex */
    if ( !IS_BNS_VT_C_GR(pBNS->vert[vc].type) ) {
        vc = vc ^ pe->neighbor12;
    }
    v0 = vc ^ pe->neighbor12; /* ChargeStruct vertex ? */
    pv0 = pBNS->vert + v0;
    if ( IS_BNS_VT_ATOM(pv0->type) ) {
        goto exit_function; /* no charge flower exists */
    }
    /* 2 edges from v0 */
    for ( i = j = 0; i < pv0->num_adj_edges && j < 3; i ++ ) {
        pe1[j] = pBNS->edge + pv0->iedge[i];
        if ( vc != ( v1[j] = pe1[j]->neighbor12 ^ v0 ) &&
             (pv1[j] = pBNS->vert + v1[j], 
              !IS_BNS_VT_ATOM(pv1[j]->type) && !IS_BNS_VT_C_GR(pv1[j]->type)) ) {
            j ++;
        }
    }
    if ( j != 2 || i != pv0->num_adj_edges ) {
        goto exit_function;
    }
    
    if ( pv1[1]->num_adj_edges == 2 &&
         pv1[0]->num_adj_edges == 3 ) {
        i0 = 1;
        i1 = 0;
    } else
    if ( pv1[0]->num_adj_edges == 2 &&
         pv1[1]->num_adj_edges == 3 ) {
        i0 = 0;
        i1 = 1;
    } else {
        goto exit_function;
    }
    /* additional check: traverse edges around v1[i1] */
    pv_t = pv1[i1];
    v_t  = v1[i1];
    for ( i = k = 0; i < pv_t->num_adj_edges; i ++ ) {
        pe_t = pBNS->edge + pv_t->iedge[i];
        v  = pe_t->neighbor12 ^ v_t; /* v1[i1] neighbor */
        if ( v == v0 ) {
            k += 1;
        }
        if ( v == v1[i0] ) {
            k += 2;
        }
        if ( IS_BNS_VT_ATOM(pBNS->vert[v].type) ) {
            k += 4;
        }
    }
    if ( k != 7 ) {
        goto exit_function;
    }
    ret = pe1[i0] - pBNS->edge;
        
exit_function:
    return ret;

}
#if (INCLUDE_NORMALIZATION_ENTRY_POINT == 1 )
/********************************************************************************************
input: allocate (num_at+num_deleted_H) atoms in inp_ATOM *at_norm, *at_fixed_bonds_out
       allocate t_group_info
*********************************************************************************************/
int NormalizeStructure( ICHICONST INPUT_PARMS *ip, STRUCT_DATA *sd, BN_STRUCT *pBNS, 
                        StrFromINChI *pStruct, inp_ATOM *at, inp_ATOM *at2,
                        VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups,
                        inp_ATOM *at_norm, inp_ATOM *at_fixed_bonds_out, T_GROUP_INFO *t_group_info )
{
    int i, ret, num_endpoints, nLenTaut;
    int num_at = pStruct->num_atoms;
    int num_deleted_H = pStruct->num_deleted_H;
    int len_at = num_at + num_deleted_H;
    /*
    T_GROUP_INFO tgi;
    T_GROUP_INFO *t_group_info = &tgi;
    inp_ATOM *at_fixed_bonds_out = NULL;
    inp_ATOM *at_norm = NULL;

    at_norm = (inp_ATOM *)inchi_calloc( len_at, sizeof(at_norm[0]) );
    at_fixed_bonds_out = (inp_ATOM *)inchi_calloc( len_at, sizeof(at_fixed_bonds_out[0]) );
    if ( !at_norm || !at_fixed_bonds_out ) {
        if ( at_norm ) inchi_free( at_norm );
        if ( at_fixed_bonds_out ) inchi_free( at_fixed_bonds_out );
        ret = RI_ERR_ALLOC;
        goto exit_function;
    }
    */
/* call normalization only */
    memset( t_group_info, 0, sizeof(t_group_info[0]) );
    t_group_info->tni.nNumRemovedExplicitH = pStruct->num_deleted_H;
    t_group_info->bTautFlags     = ip->bTautFlags;
    t_group_info->bTautFlagsDone = 0; /* (ip->bTautFlagsDone | sd->bTautFlagsDone[INCHI_BAS]);*/
    
    memcpy( at2, at, len_at*sizeof(at2[0]));
    pStruct->at = at2;
    ret = CopyBnsToAtom( pStruct, pBNS, pVA, pTCGroups, 1 );
    pStruct->at = at;
    if ( ret < 0 ) {
        goto exit_function;
    }
#if ( FIND_RING_SYSTEMS == 1 )
    ret = MarkRingSystemsInp( at2, num_at, 0 );
    if ( ret < 0 ) {
        goto exit_function;
    }
#endif
    memcpy( at_norm, at2, len_at * sizeof(at_norm[0]) );
    for ( i = 0, num_endpoints = 0; i < num_at; i ++ ) {
        num_endpoints += (0 != at_norm[i].endpoint);
        at_norm[i].endpoint = 0;
    }
    
    ret = mark_alt_bonds_and_taut_groups ( at_norm, at_fixed_bonds_out, num_at, t_group_info,
                                           NULL /* &inpbTautFlags*/, NULL /*inpbTautFlagsDone*/ );
    if ( ret < 0 ) {
        goto exit_function;/*  out of RAM or other normalization problem */
    }
    /* after normalization, t_group_info->t_group[i].num[0] = number of H + number of (-)   */
    /*                         t_group_info->t_group[i].num[1] = number of (-)              */

    /* --- count t-groups, remove (-)-only t-groups, replace -------------------------------*/
    /* t_group_info->t_group[i].num[0] with                                                 */
    /* t_group_info->t_group[i].num[0]-t_group_info->t_group[i].num[1]                      */
    nLenTaut = CountTautomerGroupsInpAt( at_norm, num_at, t_group_info );
    ret = nLenTaut;
exit_function:
    return ret;
}
#endif
/******************************************************************************************************/
int MakeOneInChIOutOfStrFromINChI2(ICHICONST INPUT_PARMS *ip_inp, STRUCT_DATA *sd_inp, 
                                   BN_STRUCT *pBNS, StrFromINChI *pStruct,
                                   inp_ATOM *at, inp_ATOM *at2, inp_ATOM *at3, 
                                   VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups,
                                   T_GROUP_INFO **t_group_info, 
                                   inp_ATOM **at_norm, inp_ATOM **at_prep )
{ 
    int ret;
    INPUT_PARMS ip_loc, *ip;
    STRUCT_DATA sd_loc, *sd;

    ip_loc = *ip_inp;
    sd_loc = *sd_inp;
    ip = &ip_loc;
    sd = &sd_loc;
    memset( sd, 0, sizeof(*sd) );
    /* create structure out of BNS */
    memcpy( at2, at, (pStruct->num_atoms + pStruct->num_deleted_H)*sizeof(at2[0]));
    pStruct->at = at2;
    ret = CopyBnsToAtom( pStruct, pBNS, pVA, pTCGroups, 1 );
    pStruct->at = at;
    if ( ret < 0 ) {
        goto exit_function;/*  out of RAM or other normalization problem */
    }
    pStruct->at = at;
    ret = MakeOneInChIOutOfStrFromINChI( ip, sd, pStruct, at2, at3, pTCGroups );
    if ( ret < 0 ) {
        goto exit_function;/*  out of RAM or other normalization problem */
    }
    if ( at_norm ) {
        *at_norm        = pStruct->pOne_norm_data[0]->at;
    }
    if ( at_prep ) {
        if ( pStruct->pOne_norm_data[0]->bTautPreprocessed && pStruct->pOne_norm_data[0]->at_fixed_bonds ) {
            *at_prep = pStruct->pOne_norm_data[0]->at_fixed_bonds;
        } else
        /* get preprocessed structure in case of Fixed-H */
        if ( pStruct->iMobileH == TAUT_NON && pStruct->pOne_norm_data[1] && pStruct->pOne_norm_data[1]->bTautPreprocessed ) {
            *at_prep = pStruct->pOne_norm_data[1]->at_fixed_bonds;
        } else {
            *at_prep = NULL;
        }
    }
    if ( t_group_info ) {
        if ( pStruct->iMobileH == TAUT_YES &&
             pStruct->One_ti.num_t_groups &&
             pStruct->One_ti.t_group && pStruct->One_ti.nEndpointAtomNumber ) {
            *t_group_info   = &pStruct->One_ti;
        } else {
            *t_group_info   = NULL;
        }
    }
exit_function:
    return ret;
}

/******************************************************************************************************/
int MakeOneInChIOutOfStrFromINChI( ICHICONST INPUT_PARMS *ip, STRUCT_DATA *sd, StrFromINChI *pStruct,
                                   inp_ATOM *at2, inp_ATOM *at3, ALL_TC_GROUPS *pTCGroups )
{

    INCHI_MODE     bTautFlags     = ip->bTautFlags | TG_FLAG_H_ALREADY_REMOVED;
    INCHI_MODE     bTautFlagsDone = 0; /*(ip->bTautFlagsDone | sd->bTautFlagsDone[INCHI_BAS]);*/
    INChI       *cur_INChI[TAUT_NUM];
    INChI_Aux   *cur_INChI_Aux[TAUT_NUM];
    int           i, j, k;
    int           iComponent  = pTCGroups->iComponent;
    int           len_at = pStruct->num_atoms + pStruct->num_deleted_H;
    int           num_atoms = pStruct->num_atoms;
    long          ulStructTime;
    
    INP_ATOM_DATA InpCurAtData;
    INP_ATOM_DATA *inp_cur_data;

    INP_ATOM_DATA InpNormAtData, InpNormTautData;
    INP_ATOM_DATA *inp_norm_data[TAUT_NUM]; /*  = { &InpNormAtData, &InpNormTautData }; */

    int            bOrigCoord = 0;
    int            num_at, ret = RI_ERR_PROGR;
    struct tagInchiTime ulMaxTime;

    T_GROUP_INFO *t_group_info = NULL;
    /* initialization */
    inp_cur_data     = &InpCurAtData;
    inp_norm_data[TAUT_NON] = &InpNormAtData;
    inp_norm_data[TAUT_YES] = &InpNormTautData;

    memset( inp_cur_data      , 0, sizeof( *inp_cur_data     ) );
    memset( inp_norm_data[TAUT_NON], 0, sizeof( *inp_norm_data[0] ) );
    memset( inp_norm_data[TAUT_YES], 0, sizeof( *inp_norm_data[0] ) );
    ulStructTime = sd->ulStructTime;
    memset( sd, 0, sizeof(*sd) );
    
    /* deallocate old results */
    free_t_group_info( &pStruct->One_ti );
    for ( k = 0; k < TAUT_NUM; k ++ ) {
        Free_INChI( &pStruct->pOneINChI[k] );
        Free_INChI_Aux( &pStruct->pOneINChI_Aux[k] );
        if ( pStruct->pOne_norm_data[k] ) {
            FreeInpAtomData( pStruct->pOne_norm_data[k] );
            inchi_free( pStruct->pOne_norm_data[k] );
            pStruct->pOne_norm_data[k] = NULL;
        }
        cur_INChI[k]      = NULL;
        cur_INChI_Aux[k]  = NULL;
    }
    memcpy( at3, at2, sizeof(at3[0])*len_at );
    /* prepare the structure */
    IncrZeroBondsAndClearEndpts(at3, num_atoms, iComponent+1);
    CopySt2At( at3, pStruct->st, pStruct->num_atoms );
    FixUnkn0DStereoBonds( at3, pStruct->num_atoms);
    ret = ReconcileAllCmlBondParities( at3, pStruct->num_atoms, 0 );
    if ( ret < 0 ) {
        goto exit_function;
    }
    if ( 0 < fix_odd_things( num_atoms, at3, 1, ip->bFixNonUniformDraw ) ) 
    {
        if ( sd->nErrorType < _IS_WARNING ) 
        {
            sd->nErrorType = _IS_WARNING;
        }
        sd->bTautFlagsDone[INCHI_BAS] |= TG_FLAG_FIX_ODD_THINGS_DONE;
    }
    /* allocate and set parameters */
    inp_cur_data->at = at3;
    inp_cur_data->num_at = num_atoms;
    inp_cur_data->num_removed_H = pStruct->num_deleted_H;

    bTautFlagsDone &= ~(TG_FLAG_FOUND_ISOTOPIC_H_DONE | TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE);

    if ( (i = bNumHeterAtomHasIsotopicH( at3, num_atoms )) ) {
        if ( i & 1 ) {
            bTautFlagsDone |= TG_FLAG_FOUND_ISOTOPIC_H_DONE;
        }
        if ( i & 2 ) {
            bTautFlagsDone |= TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE;
        }
    }

    memset( &ulMaxTime, 0, sizeof(ulMaxTime));

    /*  allocate memory for non-tautimeric (k=0) and tautomeric (k=1) results */
    for ( k = 0; k < TAUT_NUM; k ++ ) {

        if ( !pStruct->bMobileH || k == pStruct->bMobileH ) {
            /* pStruct->bMobileH=0: k = 0, 1   => allow allocation of both Fixed-H and Mobile-H InChI
               pStruct->bMobileH=1: k = 1 only => allow allocation of only Mobile-H InChI              */ 
            int nAllocMode = (k==TAUT_YES? REQ_MODE_TAUT:0) |
                             (bTautFlagsDone & ( TG_FLAG_FOUND_ISOTOPIC_H_DONE |
                                                 TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE ))?
                             (ip->nMode & REQ_MODE_ISO):0;

            if ( (k==TAUT_NON && (ip->nMode & REQ_MODE_BASIC )) ||
                 (k==TAUT_YES && (ip->nMode & REQ_MODE_TAUT ))     ) {
                /*  alloc INChI and INChI_Aux only if ip->nMode allows this */
                cur_INChI[k]     = Alloc_INChI( inp_cur_data->at, inp_cur_data->num_at, &inp_cur_data->num_bonds,
                                              &inp_cur_data->num_isotopic, nAllocMode );
                cur_INChI_Aux[k] = Alloc_INChI_Aux( inp_cur_data->num_at,
                                              inp_cur_data->num_isotopic, nAllocMode, bOrigCoord );
                if ( cur_INChI_Aux[k] ) {
                    cur_INChI_Aux[k]->bIsIsotopic = inp_cur_data->num_isotopic;
                }
                /*  alloc memory for the output structure: non-tautomeric and tautomeric (for displaying) */
                CreateInpAtomData( inp_norm_data[k], inp_cur_data->num_at+inp_cur_data->num_removed_H, k );
                inp_norm_data[k]->num_removed_H = inp_cur_data->num_removed_H;
            } else {
                FreeInpAtomData( inp_norm_data[k] );
            }
        } else {
            FreeInpAtomData( inp_norm_data[k] );
        }
    }
    k = pStruct->bMobileH;
    /* In case of Fixed-H we have to create InChI for both Fixed-H and Mobile-H */
    num_at = Create_INChI( cur_INChI, cur_INChI_Aux, NULL/* not used */, inp_cur_data->at,
                          inp_norm_data,
                          inp_cur_data->num_at+inp_cur_data->num_removed_H,
                          ip->nMode, &bTautFlags, &bTautFlagsDone, NULL /* &ulMaxTime*/, &pStruct->One_ti, sd->pStrErrStruct);
    SetConnectedComponentNumber( inp_cur_data->at, inp_cur_data->num_at, iComponent+1 ); /*  normalization alters structure component number */
    
    /* detect InChI errors */

    if ( num_at < 0 ) {
        ret = num_at;
    } else
    if ( cur_INChI[k] && cur_INChI[k]->nErrorCode ) {
        ret = cur_INChI[k]->nErrorCode;
    } else
    if ( cur_INChI_Aux[k] && cur_INChI_Aux[k]->nErrorCode ) {
        ret = cur_INChI_Aux[k]->nErrorCode;
    } else {
        ret = 0;
    }

    /* fill out the output */

    if ( !ret ) {
        int bMobileH = pStruct->bMobileH;
        if ( bMobileH == TAUT_NON &&
             0 == cur_INChI[TAUT_NON]->nNumberOfAtoms &&
             0 <  cur_INChI[TAUT_YES]->nNumberOfAtoms ) {
            /* tautomerism or H(+) removal/addition was not discovered */
            bMobileH = TAUT_YES;
        }
        pStruct->nChargeRevrs = cur_INChI[TAUT_YES]->nTotalCharge;

        pStruct->pOneINChI[0]       = cur_INChI[bMobileH];
        pStruct->pOneINChI_Aux[0]   = cur_INChI_Aux[bMobileH];
        pStruct->nOneINChI_bMobileH = bMobileH;
        cur_INChI[bMobileH]         = NULL;  /* remove pointer to avoid deallocation at exit_function */
        cur_INChI_Aux[bMobileH]     = NULL;  /* remove pointer to avoid deallocation at exit_function */

        pStruct->nNumRemovedProtons = (pStruct->iMobileH == TAUT_YES)? pStruct->One_ti.tni.nNumRemovedProtons : 0;


        /* set correct t-group numbers to endpoints */
        t_group_info       = &pStruct->One_ti;
        if ( t_group_info->num_t_groups && t_group_info->t_group && t_group_info->nEndpointAtomNumber ) {
            inp_ATOM     *at_norm      = inp_norm_data[TAUT_YES]->at;
            int          num_at_norm   = inp_norm_data[TAUT_YES]->num_at;
            for ( i = 0; i < num_at_norm; i ++ ) {
                at_norm[i].endpoint = 0;
            }
            for ( i = 0; i < t_group_info->num_t_groups; i ++ ) {
                k = t_group_info->t_group[i].nFirstEndpointAtNoPos;
                /* add number of mobile (-) to the number of mobile H */
                t_group_info->t_group[i].num[0] += t_group_info->t_group[i].num[1];
                for ( j = 0; j < t_group_info->t_group[i].nNumEndpoints; j ++, k ++ ) {
                    at_norm[t_group_info->nEndpointAtomNumber[k]].endpoint = t_group_info->t_group[i].nGroupNumber;
                }
            }
        }
        pStruct->pOne_norm_data[0] = (INP_ATOM_DATA *) inchi_malloc( sizeof(pStruct->pOne_norm_data[0][0]) );
        if ( pStruct->pOne_norm_data[0] ) {
            memcpy( pStruct->pOne_norm_data[0], inp_norm_data[bMobileH], sizeof(pStruct->pOne_norm_data[0][0]));
            memset( inp_norm_data[bMobileH], 0, sizeof(*inp_norm_data[0]) );
        } else {
            ret = RI_ERR_ALLOC;
        }
        if ( bMobileH == TAUT_NON && cur_INChI[TAUT_YES]->nNumberOfAtoms > 0 ) {
            int bMobileHalt = ALT_TAUT(bMobileH); /* = TAUT_YES */
            pStruct->pOneINChI[1]       = cur_INChI[bMobileHalt];
            pStruct->pOneINChI_Aux[1]   = cur_INChI_Aux[bMobileHalt];
            cur_INChI[bMobileHalt]         = NULL;
            cur_INChI_Aux[bMobileHalt]     = NULL;
            pStruct->pOne_norm_data[1] = (INP_ATOM_DATA *) inchi_malloc( sizeof(pStruct->pOne_norm_data[0][0]) );
            if ( pStruct->pOne_norm_data[1] ) {
                memcpy( pStruct->pOne_norm_data[1], inp_norm_data[bMobileHalt], sizeof(pStruct->pOne_norm_data[0][0]));
                memset( inp_norm_data[bMobileHalt], 0, sizeof(*inp_norm_data[0]) );
            } else {
                ret = RI_ERR_ALLOC;
            }
        }
    } else {
#if ( bRELEASE_VERSION != 1 )
#ifndef TARGET_API_LIB
        fprintf( stdout, "ERROR: Create_INChI returned %d\n", ret );
#endif
#endif
    }

exit_function:
    /* deallocate unused */
    for ( k = 0; k < TAUT_NUM; k ++ ) {
        Free_INChI( &cur_INChI[k] );
        Free_INChI_Aux( &cur_INChI_Aux[k] );
        FreeInpAtomData( inp_norm_data[k] );
    }
    sd->ulStructTime = ulStructTime;

    return ret;
}

/******************************************************************************************************
Input: 
       at[].num_H       = total number of all terminal H connected to the atom
       at[].num_iso_H[] = numbers of isotopic H among at[].num_H
       Explicit H are disconnected
       Calculate InChI with normalization only in MakeOneInChIOutOfStrFromINChI()
       with (TG_FLAG_H_ALREADY_REMOVED & bTautFlags) != 0
Output:
       at[].num_H       = number of implicit non-isotopic H connected to the atom
       at[].num_iso_H[] = numbers of implicit isotopic H (not included in at[].num_H)
       Explicit H are connected
       Calculate InChI with full preprocessing MakeInChIOutOfStrFromINChI2()
       with (TG_FLAG_H_ALREADY_REMOVED & bTautFlags) == 0
*******************************************************************************************************/
int ConnectDisconnectedH( inp_ATOM *at, int num_atoms, int num_deleted_H )
{
    int i, j, k, n, m, num_H;
    int tot_atoms = num_atoms + num_deleted_H;
    
    for ( i = num_atoms; i < tot_atoms; i = j ) {
        k = at[i].neighbor[0]; /* a[k] is the atom connected to the explicit hydrogen at[i] */
        for ( j = i; j < tot_atoms && at[j].neighbor[0] == k; j ++ )
            ;
        num_H = j-i; /* number of explicit H for at[k] */
        if ( num_H > at[k].num_H ) {
            return RI_ERR_PROGR;
        }
        if ( num_H + at[k].valence > MAXVAL ) {
            return RI_ERR_SYNTAX;
        }
        /* insert links to explicit H before all other links in the connection list */
        n = at[k].valence;
        memmove( at[k].neighbor   +num_H, at[k].neighbor,    sizeof(at[k].neighbor[0]) * n );
        memmove( at[k].bond_stereo+num_H, at[k].bond_stereo, sizeof(at[k].bond_stereo[0]) * n );
        memmove( at[k].bond_type  +num_H, at[k].bond_type   , sizeof(at[k].bond_type[0]) * n );
        for ( n = 0; n < num_H; n ++ ) {
            at[k].neighbor[n]    = i + n;
            at[k].bond_stereo[n] = 0;
            at[k].bond_type[n]   = BOND_TYPE_SINGLE;
        }
        for ( m = 0; m < MAX_NUM_STEREO_BONDS && at[k].sb_parity[m]; m ++ ) {
            at[k].sb_ord[m] += num_H;
            if ( at[k].sn_ord[m] < 0 ) {
                for ( n = i; n < j; n ++ ) {
                    if ( at[n].orig_at_number == at[k].sn_orig_at_num[m] ) {
                        at[k].sn_ord[m] = n-i;
                        break;
                    }
                }
                if ( n == j ) {
                    return RI_ERR_PROGR;
                }
            } else {
                at[k].sn_ord[m] += num_H;
            }
        }
        at[k].valence            += num_H;
        at[k].chem_bonds_valence += num_H;
        at[k].num_H -= num_H; /* cannot be negative */
        /*memset( at[k].num_iso_H, 0, sizeof(at[0].num_iso_H) );*/ /* attached H must carry all isotopic shifts */
        for ( n = i; n < j; n ++ ) {
            at[n].chem_bonds_valence = BOND_TYPE_SINGLE;
        }
        /* isotopic H */
        for ( m = j-1; i <= m && at[m].iso_atw_diff > 0 ; m -- ) {
            if ( at[m].iso_atw_diff > NUM_H_ISOTOPES ) {
                return RI_ERR_PROGR;
            }
            if ( 0 >= at[k].num_iso_H[(int)at[m].iso_atw_diff-1] -- ) {
                return RI_ERR_PROGR;
            }
        }

    }
    /* subtract isotopic H */
    for ( i = 0; i < num_atoms; i ++ ) {
        for ( m = 0; m < NUM_H_ISOTOPES; m ++ ) {
            at[i].num_H -= at[i].num_iso_H[m];
        }
        if ( 0 > at[i].num_H ) {
            return RI_ERR_PROGR;
        }
    }

    return tot_atoms;
}
/******************************************************************************************************
Input:
       at[].num_H       = number of implicit non-isotopic H connected to the atom
       at[].num_iso_H[] = numbers of implicit isotopic H (not included in at[].num_H)
       Explicit H are connected
       Calculate InChI with (TG_FLAG_H_ALREADY_REMOVED & bTautFlags) == 0
Output:
       at[].num_H       = total number of all terminal H connected to the atom
       at[].num_iso_H[] = numbers of isotopic H among at[].num_H
       Explicit H are disconnected
       Calculate InChI with (TG_FLAG_H_ALREADY_REMOVED & bTautFlags) != 0
*******************************************************************************************************/
int DisconnectedConnectedH( inp_ATOM *at, int num_atoms, int num_deleted_H )
{
    int i, j, k, n, m, num_H, num_iso_H;
    int tot_atoms = num_atoms + num_deleted_H;
    
    /* add implicit isotopic H to total implicit H */
    for ( i = 0; i < num_atoms; i ++ ) {
        for ( m = 0; m < NUM_H_ISOTOPES; m ++ ) {
            at[i].num_H += at[i].num_iso_H[m];
        }
    }
    for ( i = num_atoms; i < tot_atoms; i = j ) {
        k = at[i].neighbor[0]; /* a[k] is the atom connected to the explicit hydrogen at[i] */
        for ( j = i; j < tot_atoms && at[j].neighbor[0] == k; j ++ ) {
            at[j].chem_bonds_valence = 0;
        }
        num_H = j-i; /* number of explicit H for at[k] */
        /* verify correct number of explicit H */
        for ( n = 0; n < at[k].valence && at[k].neighbor[n] >= num_atoms; n ++ )
            ;
        if ( n != num_H ) {
            return RI_ERR_PROGR;
        }
        /* remove bonds to explicit H located in front of all other bonds in the connection list */
        n = (at[k].valence -= num_H); /* new number of bonds */
        at[k].chem_bonds_valence -= num_H; /* new no-H valence */
        if ( n ) {
            memmove( at[k].neighbor,    at[k].neighbor    + num_H, sizeof(at[k].neighbor[0]) * n );
            memmove( at[k].bond_stereo, at[k].bond_stereo + num_H, sizeof(at[k].bond_stereo[0]) * n );
            memmove( at[k].bond_type,   at[k].bond_type   + num_H, sizeof(at[k].bond_type[0]) * n );
        }
        /* clear the 'tails' */
        memset( at[k].neighbor+n,    0, sizeof(at[k].neighbor[0]) * num_H );
        memset( at[k].bond_stereo+n, 0, sizeof(at[k].bond_stereo[0]) * num_H );
        memset( at[k].bond_type+n,   0, sizeof(at[k].bond_type[0]) * num_H );

        for ( m = 0; m < MAX_NUM_STEREO_BONDS && at[k].sb_parity[m]; m ++ ) {
            at[k].sb_ord[m] -= num_H;
            if ( 0 <= at[k].sn_ord[m] && at[k].sn_ord[m] < num_H ) {
                at[k].sn_ord[m] = -1; /* disconnected explicit H */
            }
        }
        /* add explicit isotopic H (already included in num_H) */
        for ( num_iso_H = 0, m = j-1; i <= m && at[m].iso_atw_diff > 0 ; m -- ) {
            if ( at[m].iso_atw_diff > NUM_H_ISOTOPES ) {
                return RI_ERR_PROGR;
            }
            at[k].num_iso_H[(int)at[m].iso_atw_diff-1] ++;
        }
        at[k].num_H += num_H; /* add all explicit H including isotopic */
    }
    return tot_atoms;
}
/******************************************************************************************************/
int MakeInChIOutOfStrFromINChI2( ICHICONST INPUT_PARMS *ip_inp, STRUCT_DATA *sd_inp, StrFromINChI *pStruct,
                                 int iComponent, int iAtNoOffset, long num_inp )
{
    char szTitle[MAX_SDF_HEADER+MAX_SDF_VALUE+256];

    int   len, ret;
    /*
    PINChI2     *pINChI[INCHI_NUM];
    PINChI_Aux2 *pINChI_Aux[INCHI_NUM];
    */
    char        pStr[256];
    INPUT_PARMS local_ip;
    STRUCT_DATA local_sd;
    INPUT_PARMS *ip = &local_ip;
    STRUCT_DATA *sd = &local_sd;
    
    ORIG_ATOM_DATA OrigAtData; /* 0=> disconnected, 1=> original */
    ORIG_ATOM_DATA *orig_inp_data = &OrigAtData;
    ORIG_ATOM_DATA PrepAtData[2]; /* 0=> disconnected, 1=> original */
    ORIG_ATOM_DATA *prep_inp_data = PrepAtData;
    
    *ip = *ip_inp;
    ip->bDisplay = 0;
    ip->bDisplayCompositeResults = 0;
    ip->bDisplayEachComponentINChI = 0;
    ip->bDisplayIfRestoreWarnings  = 0;
    ip->bINChIOutputOptions = INCHI_OUT_NO_AUX_INFO;
    /*
    if ( pStruct->bMobileH ) {
        ip->nMode &= ~REQ_MODE_BASIC;
        ip->nMode |= REQ_MODE_TAUT;
    } else {
        ip->nMode |= (REQ_MODE_TAUT | REQ_MODE_BASIC);
    }
    */
    memset( sd, 0, sizeof(*sd) );
    sd->fPtrStart = -1;
    sd->fPtrEnd = -1;
    /*
    if ( ip->nMode & REQ_MODE_STEREO ) {
        if ( ip->nMode & (REQ_MODE_RELATIVE_STEREO | REQ_MODE_RACEMIC_STEREO) ) {
            sd->bChiralFlag |= FLAG_INP_AT_NONCHIRAL;
        } else {
            sd->bChiralFlag |= FLAG_INP_AT_CHIRAL;
        }
    }
    */
    memset( orig_inp_data     , 0,   sizeof( *orig_inp_data  ) );
    memset( prep_inp_data     , 0, 2*sizeof( *prep_inp_data  ) );
    memset( pStruct->RevInChI.pINChI,     0, sizeof(pStruct->RevInChI.pINChI    ) );
    memset( pStruct->RevInChI.pINChI_Aux, 0, sizeof(pStruct->RevInChI.pINChI_Aux) );
    memset( pStr, 0, sizeof(pStr) );
    memset( szTitle, 0, sizeof(szTitle) );
    
    len = sizeof(orig_inp_data->at[0])*(pStruct->num_atoms + pStruct->num_deleted_H);
    orig_inp_data->at = (inp_ATOM *) inchi_malloc( len );
    if ( orig_inp_data->at ) {
        /*memcpy( orig_inp_data->at, pStruct->at2, len );*/
        /*ret = ConnectDisconnectedH( orig_inp_data->at, pStruct->num_atoms, pStruct->num_deleted_H );*/
        CopySt2At( pStruct->at2, pStruct->st, pStruct->num_atoms );
        ret = ConnectDisconnectedH( pStruct->at2, pStruct->num_atoms, pStruct->num_deleted_H );
        if ( ret < 0 ) {
            goto exit_error;
        }
        orig_inp_data->num_inp_atoms = ret;
        /* connections changed => reconcile parities even if they were reconciled before */
        /* remove t-group markings and increment zero-order bonds,
           otherwise MakeInChIOutOfStrFromINChI2() woild fail */
        /*
        IncrZeroBondsAndClearEndpts(pStruct->at2, pStruct->num_atoms, iComponent+1);
        */
        IncrZeroBonds(pStruct->at2, pStruct->num_atoms, iComponent+1);
        
        /* CopySt2At() moved to the position before ConnectDisconnectedH() because 
           in case stereo exists only in Mobile-H layer and the processd here
           component is restored in Fixed-H layer the parities needed by
           ConnectDisconnectedH() must be there before calling
           ConnectDisconnectedH()
        */
        /*CopySt2At( pStruct->at2, pStruct->st, pStruct->num_atoms );*/

        ret = ReconcileAllCmlBondParities( pStruct->at2, orig_inp_data->num_inp_atoms, 0 );
        if ( ret < 0 ) {
            goto exit_error;
        }
        memcpy( orig_inp_data->at, pStruct->at2, len );
        ClearEndpts(orig_inp_data->at, pStruct->num_atoms);
        if ( FixUnkn0DStereoBonds(orig_inp_data->at, pStruct->num_atoms) ) {
            ret = ReconcileAllCmlBondParities( pStruct->at2, orig_inp_data->num_inp_atoms, 0 );
            if ( ret < 0 ) {
                goto exit_error;
            }
        }
        /* keep endpoint[] markings in at2[] for subsequent add/remove protons */
    } else {
        ret = RI_ERR_ALLOC;
        goto exit_error;
    }
    memset( sd->num_components, 0, sizeof(sd->num_components) );
    memset( sd->num_taut, 0, sizeof(sd->num_taut) );
    memset( sd->num_non_taut, 0, sizeof(sd->num_non_taut) );
    memset( sd->bTautFlagsDone, 0, sizeof(sd->bTautFlagsDone) );
    memset( sd->bTautFlags, 0, sizeof(sd->bTautFlags) );

    ret = ProcessOneStructure( sd, ip, szTitle, pStruct->RevInChI.pINChI, pStruct->RevInChI.pINChI_Aux,
                             NULL /*inp_file*/, NULL /*log_file*/, NULL /*output_file*/, NULL /*prb_file*/,
                             orig_inp_data, prep_inp_data,
                             num_inp, pStr, sizeof(pStr), 
                             0 /* save_opt_bits */);

    memcpy(pStruct->RevInChI.num_components, sd->num_components, sizeof(pStruct->RevInChI.num_components) );
    memcpy(sd_inp->pStrErrStruct, sd->pStrErrStruct, sizeof(sd_inp->pStrErrStruct) );
    pStruct->RevInChI.nRetVal = ret;
    /* translate returned value */
    if ( ret == _IS_ERROR || ret == _IS_FATAL || ret == _IS_UNKNOWN ) {
        ret = RI_ERR_PROGR;
    } else
    if ( ret == _IS_OKAY ) {
        ret = 0;
    } else
    if ( ret == _IS_WARNING ) {
        ret = 1;
    } else {
        ret = RI_ERR_PROGR;
    }
    /* save total charge from Mobile-H layer */
    pStruct->nChargeRevrs = 0;
    if ( ret >= 0 ) {
        if ( bRevInchiComponentExists( pStruct, INCHI_REC, TAUT_YES, 0 ) ) {
            pStruct->nChargeRevrs = pStruct->RevInChI.pINChI[INCHI_REC][0][TAUT_YES]->nTotalCharge;
        } else
        if ( bRevInchiComponentExists( pStruct, INCHI_BAS, TAUT_YES, 0 ) ) {
            pStruct->nChargeRevrs = pStruct->RevInChI.pINChI[INCHI_BAS][0][TAUT_YES]->nTotalCharge;
        }
    }
    
    /* free structure data */
    FreeOrigAtData( orig_inp_data );
    FreeOrigAtData( prep_inp_data );
    FreeOrigAtData( prep_inp_data+1 );

exit_error:
    return ret;
}
/******************************************************************************************************/
int OutputInChIOutOfStrFromINChI(ICHICONST INPUT_PARMS *ip_inp, STRUCT_DATA *sd_inp, 
                                 long num_inp, int bINChIOutputOptions,
                                 INCHI_IOSTREAM *pout, INCHI_IOSTREAM *plog, 
                                 InpInChI *pOneInput, int bHasSomeFixedH,
                                 unsigned char save_opt_bits)
{
    char szTitle[MAX_SDF_HEADER+MAX_SDF_VALUE+256];

    int   len, ret;
/*    
    PINChI2     *pINChI[INCHI_NUM];
    PINChI_Aux2 *pINChI_Aux[INCHI_NUM];
*/    
    REV_INCHI   RevInChI;
    int nStrLen = INCHI_SEGM_BUFLEN;
    char *pStr  = NULL;

    INPUT_PARMS local_ip;
    STRUCT_DATA local_sd;
    INPUT_PARMS *ip = &local_ip;
    STRUCT_DATA *sd = &local_sd;
    
    ORIG_ATOM_DATA OrigAtData; /* 0=> disconnected, 1=> original */
    ORIG_ATOM_DATA *orig_inp_data = &OrigAtData;
    ORIG_ATOM_DATA PrepAtData[2]; /* 0=> disconnected, 1=> original */
    ORIG_ATOM_DATA *prep_inp_data = PrepAtData;
    
    *ip = *ip_inp;
    ip->bNoStructLabels = 1;
    ip->bDisplay = 0;
    ip->bDisplayCompositeResults = 0;
    ip->bDisplayEachComponentINChI = 0;
    ip->bDisplayIfRestoreWarnings  = 0;
#if ( I2S_MODIFY_OUTPUT == 1 )
    if ( bINChIOutputOptions & INCHI_OUT_SDFILE_ONLY )
        ip->bINChIOutputOptions = bINChIOutputOptions & ~(INCHI_OUT_PLAIN_TEXT | INCHI_OUT_XML | INCHI_OUT_PLAIN_TEXT_COMMENTS | INCHI_OUT_XML_TEXT_COMMENTS);
    else
    if ( bINChIOutputOptions & INCHI_OUT_XML )
        ip->bINChIOutputOptions = bINChIOutputOptions & ~(INCHI_OUT_PLAIN_TEXT | INCHI_OUT_SDFILE_ONLY) | INCHI_OUT_EMBED_REC;
    else
    if ( bINChIOutputOptions & INCHI_OUT_PLAIN_TEXT )
        ip->bINChIOutputOptions = bINChIOutputOptions & ~(INCHI_OUT_XML | INCHI_OUT_SDFILE_ONLY) | INCHI_OUT_EMBED_REC;
    else
    if ( bINChIOutputOptions & (INCHI_OUT_NO_AUX_INFO | INCHI_OUT_SHORT_AUX_INFO | INCHI_OUT_ONLY_AUX_INFO | INCHI_OUT_TABBED_OUTPUT))
        ip->bINChIOutputOptions = (INCHI_OUT_PLAIN_TEXT | INCHI_OUT_EMBED_REC | bINChIOutputOptions);
    else
        ip->bINChIOutputOptions = (INCHI_OUT_PLAIN_TEXT | INCHI_OUT_EMBED_REC);
#else
    ip->bINChIOutputOptions = (INCHI_OUT_PLAIN_TEXT | INCHI_OUT_EMBED_REC );
#endif

    if ( bHasSomeFixedH ) 
    {
        ip->nMode |= (REQ_MODE_TAUT | REQ_MODE_BASIC);
    } else 
    {
        ip->nMode &= ~REQ_MODE_BASIC;
        ip->nMode |= REQ_MODE_TAUT;
    }

    memset( sd, 0, sizeof(*sd) );
    sd->fPtrStart = -1;
    sd->fPtrEnd = -1;
    /*
    if ( ip->nMode & REQ_MODE_STEREO ) {
        if ( ip->nMode & (REQ_MODE_RELATIVE_STEREO | REQ_MODE_RACEMIC_STEREO) ) {
            sd->bChiralFlag |= FLAG_INP_AT_NONCHIRAL;
        } else {
            sd->bChiralFlag |= FLAG_INP_AT_CHIRAL;
        }
    }
    */
    memset( orig_inp_data,       0,   sizeof( *orig_inp_data  ) );
    memset( prep_inp_data,       0, 2*sizeof( *prep_inp_data  ) );
    memset( RevInChI.pINChI,     0, sizeof(RevInChI.pINChI    ) );
    memset( RevInChI.pINChI_Aux, 0, sizeof(RevInChI.pINChI_Aux) );
    
    len = sizeof(orig_inp_data->at[0]) * pOneInput->num_atoms;
    orig_inp_data->at = (inp_ATOM *) inchi_malloc( len );
    orig_inp_data->szCoord = (MOL_COORD *)inchi_calloc( pOneInput->num_atoms, sizeof(orig_inp_data->szCoord[0]));
    pStr  = (char *)inchi_calloc( nStrLen, sizeof(char) );
    if ( orig_inp_data->at && orig_inp_data->szCoord && pStr ) {
        int i, k;
        memcpy( orig_inp_data->at, pOneInput->atom, len );
        orig_inp_data->num_inp_atoms = pOneInput->num_atoms;
        ClearEndpts( orig_inp_data->at, orig_inp_data->num_inp_atoms );
        /* otherwise fails on CID=450438 */
        if ( FixUnkn0DStereoBonds(orig_inp_data->at, orig_inp_data->num_inp_atoms) ) {
            ret = ReconcileAllCmlBondParities( orig_inp_data->at, orig_inp_data->num_inp_atoms, 0 );
            if ( ret < 0 ) {
                goto exit_error;
            }
        }
        /* To obtain rA,rB,rC in AuxInfo we have to emulate input coordinates; make all of them zeroes */
        for ( i = 0; i < pOneInput->num_atoms; i ++ ) {
            for ( k = 0; k < NUM_COORD*LEN_COORD; k += LEN_COORD ) {
                orig_inp_data->szCoord[i][k] = '0';
            }
        }
    } else {
        ret = RI_ERR_ALLOC;
        goto exit_error;
    }
    memset( sd->num_components, 0, sizeof(sd->num_components) );
    memset( sd->num_taut, 0, sizeof(sd->num_taut) );
    memset( sd->num_non_taut, 0, sizeof(sd->num_non_taut) );
    memset( sd->bTautFlagsDone, 0, sizeof(sd->bTautFlagsDone) );
    memset( sd->bTautFlags, 0, sizeof(sd->bTautFlags) );
    memset( szTitle, 0, sizeof(szTitle) );

    ret = ProcessOneStructure(sd, ip, szTitle, RevInChI.pINChI, RevInChI.pINChI_Aux,
                             NULL /*inp_file*/, plog /*log_file*/, pout /*output_file*/, NULL /*prb_file*/,
                             orig_inp_data, prep_inp_data,
                             num_inp, pStr, nStrLen,
                             save_opt_bits);
    memcpy(RevInChI.num_components, sd->num_components, sizeof(RevInChI.num_components) );
    /*
    memcpy(sd_inp->pStrErrStruct, sd->pStrErrStruct, sizeof(sd_inp->pStrErrStruct) );
    */
    RevInChI.nRetVal = ret;
    /* translate returned value */
    if ( ret == _IS_ERROR || ret == _IS_FATAL || ret == _IS_UNKNOWN ) {
        ret = RI_ERR_PROGR;
    } else
    if ( ret == _IS_OKAY ) {
        ret = 0;
    } else
    if ( ret == _IS_WARNING ) {
        ret = 1;
    } else {
        ret = RI_ERR_PROGR;
    }
    
    /* free structure data */
    FreeOrigAtData( orig_inp_data );
    FreeOrigAtData( prep_inp_data );
    FreeOrigAtData( prep_inp_data+1 );
    FreeAllINChIArrays( RevInChI.pINChI,
                        RevInChI.pINChI_Aux,
                        RevInChI.num_components );


exit_error:
    if ( pStr ) inchi_free( pStr );
    return ret;
}
#endif

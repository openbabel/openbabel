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
#include <ctype.h>
#include <math.h>

#include "mode.h"

#include "inpdef.h"
#include "util.h"
#include "ichi.h"
#include "strutil.h"
#include "ichierr.h"

#include "ichicomp.h"
#include "extr_ct.h"
#include "ichister.h"

#include "ichi_io.h"


#define FIX_P_IV_Plus_O_Minus    /* added fix to remove_ion_pairs() -- 2010-03-17 DT */

/* local prototypes */
int cmp_components( const void *a1, const void *a2 );
/*int mark_one_struct_component( inp_ATOM* at, int j, AT_NUMB *mark, AT_NUMB num_disconnected_components );*/
INChI_Stereo *Alloc_INChI_Stereo(int num_at, int num_bonds);
int RemoveInpAtBond( inp_ATOM *at, int iat, int k );
int DisconnectInpAtBond( inp_ATOM *at, AT_NUMB *nOldCompNumber, int iat, int neigh_ord );
int move_explicit_Hcation(inp_ATOM *at, int num_at, int iat, int iat_H, int bInAllComponents);
int DisconnectOneLigand( inp_ATOM *at, AT_NUMB *nOldCompNumber, S_CHAR *bMetal, char *elnumber_Heteroat,
                         int num_halogens, int num_atoms, int iMetal, int jLigand, INCHI_MODE *bTautFlagsDone );
int bIsAmmoniumSalt( inp_ATOM *at, int i, int *piO, int *pk, S_CHAR *num_explicit_H );
int DisconnectAmmoniumSalt ( inp_ATOM *at, int i, int iO, int k, S_CHAR *num_explicit_H );
/*int bIsMetalSalt( inp_ATOM *at, int i ); - moved to strutil,h */
int DisconnectMetalSalt( inp_ATOM *at, int i );
int bIsMetalToDisconnect(inp_ATOM *at, int i, int bCheckMetalValence);

int get_iat_number( int el_number, const int el_num[], int el_num_len );
int tot_unsat( int unsat[] );
int max_unsat( int unsat[] );

double dist3D( inp_ATOM *at1, inp_ATOM *at2 );
double dist2D( inp_ATOM *at1, inp_ATOM *at2 );
double dist_from_segm( double x, double y, double x1, double y1, double x2, double y2);
int segments_intersect( double x11, double y11, double x12, double y12, /* segment #1 */
                        double x21, double y21, double x22, double y22 );
double GetMinDistDistribution( inp_ATOM *at, int num_at, int iat, int iat_H,
                                int bInAllComponents, double min_dist[], int num_segm );

int nFindOneOM(inp_ATOM *at, int at_no, int ord_OM[], int num_OM);
int the_only_doublet_neigh(inp_ATOM *at, int i1, int *ineigh1, int *ineigh2);


#ifndef NUMH
#define NUM_ISO_H(AT,N) (AT[N].num_iso_H[0]+AT[N].num_iso_H[1]+AT[N].num_iso_H[2])
#define NUMH(AT,N)     (AT[N].num_H+NUM_ISO_H(AT,N))
#endif

/************************************************************************/
int the_only_doublet_neigh(inp_ATOM *at, int i1, int *ineigh1, int *ineigh2)
{
    int i, neigh1, num_rad1=0, num_rad2=0;
    inp_ATOM *a = at+i1, *b;
    if ( RADICAL_DOUBLET != a->radical )
        return -1;
    for ( i = 0; i < a->valence; i ++ ) {
        b = at + (neigh1 = (int)a->neighbor[i]);
        if ( RADICAL_DOUBLET == b->radical ) {
            num_rad1 ++;
            *ineigh1 = i;
        }
    }
    if ( 1 == num_rad1 ) {
        a = at + (neigh1 = (int)a->neighbor[*ineigh1]);
        for ( i = 0; i < a->valence; i ++ ) {
            b = at +(int)a->neighbor[i];
            if ( RADICAL_DOUBLET == b->radical ) {
                num_rad2 ++;
                *ineigh2 = i;
            }
        }
        if ( 1 == num_rad2 ) {
            return neigh1;
        }
    }
    return -1;
}

/************************************************************************/
int fix_odd_things( int num_atoms, inp_ATOM *at, int bFixBug, int bFixNonUniformDraw )
{   /*                           0 1 2  3  4 5 6  7                       8  9  */
    static const char    el[] = "N;P;As;Sb;O;S;Se;Te;";   /* 8 elements + C, Si */
    static U_CHAR  en[10];              /* same number: 8 elements */
    static int     ne=0, ne2;           /* will be 8 and 10 */
    static int     el_number_P;
    static int     el_number_H;
    static int     el_number_C;
    static int     el_number_O;
    static int     el_number_Si;

#define FIRST_NEIGHB2  4
#define FIRST_CENTER2  5
#define NUM_CENTERS_N  4

    int i1, i2, k1, k2, c, num_changes = 0;
    char elname[ATOM_EL_LEN];


    /* constants for element numbers */
    enum elems { dNone, dCl=17,dBr=35,dI=53,dAt=85,dO=8,dS=16,dSe=34,dTe=52, dP=15, dC=6, dN=7 } ; 




    if (bFixNonUniformDraw)
    {
    

        /* Correct non-uniformly drawn oxoanions and amidinium cations. */                
        {

            /* For central halogen, apply the following correftion rules:
    
                O                             O(-)
                ||                            |
              O=Hal(-)=O         ===>        O=Hal=O
                ||                            ||
                O                             O
        
            (perchlorate, etc.)


                O                             O(-)
                ||                            |
                Hal(-)=O         ===>         Hal=O
                ||                            ||
                O                             O

            (chlorate, etc.)

                O                             O(-)
                ||                            |
                Hal(-)=O         ===>         Hal=O
            

            (chlorite, etc.)

                Hal(-)=O         ===>         Hal-O(-)
    

            (hypochlorite, etc.)



            For halcogenes (S, Se, Te)
                Y                               Y(-)
                ||                              |
              RnX(-)            ===>          RnX

           if:
           1) (X = S, Y = O) || (X = Se, Y = S, O) || (X = Te, Y = O, S, Se)
           2) valence of X exceeds 6, in initially drawn form 


            So the following is corrected:  

                O                               O(-)
                ||                              |
              O=S(-)-R            ===>        O=S-R
                ||                              ||
                O                               O

            or

                  O                              O-
                  ||                             |
                F5Te(-)            ===>        F5Te

            but the following remains unchanged:


                O                
                ||               
              O=S(-)-R           


            The central atom (of IUPAC Group 16-17) is shown as negative but it contains double bond(s) 
            to terminal atom of greater electronegativity (of Group 16).
            The central atom is halcogen (S,Se,Te) in highest oxidation state or halogen.

            Fix:
            move negative charge to terminal atom and change double bond to a single one.

            Eligible central atom    Eligible terminal atom at double bond's end
                        Cl              O
                        Br              O
                        I               O
                       [At              S,Se,Te]
                        S               O
                        Se              O,S
                        Te              O, S, Se
          
            Comments: 
            1. Central atoms of Groups 13-15 are not considered.
            2. Pauling electronegativities are: 
            F(3.98) > O(3.44) > Cl (3.16) > N (3.04) > Br(2.96) > I(2.66) > S(2.58) > Se(2.55) > At (2.2) > Te(2.1)


            */
    

            static U_CHAR  allowed_elnums_center_halogen[] = {dCl, dBr, dI, dAt} ; 
            static U_CHAR  allowed_elnums_center_halcogen[] = {dS, dSe, dTe} ; 
   
            int en_center;
            int i, j, k; 

            for (i=0; i<num_atoms; i++)
            {
                /* Find appropriate central atom. This center should be: ...*/ 
        
                /* charged exactly (-1) ... */
                if ( at[i].charge != -1 ) 
                    continue;
                en_center = at[i].el_number;
                /*  from eligible element list ... */
                if (!memchr(allowed_elnums_center_halogen, en_center , sizeof(allowed_elnums_center_halogen)) ) 
                {
                    /* central atom is not not halogen; check if it is halcogen */
                    if ( memchr(allowed_elnums_center_halcogen, en_center , sizeof(allowed_elnums_center_halcogen)) ) 
                    {
    
                        if (at[i].chem_bonds_valence<7)
                        /* central atom is anionic halcogen, but not in the highest oxidation state */
                            continue;
                    }
                    else
                        continue;
                } 
        
                /* OK, found central halogen or eligible central halcogen. */            
            
                /* non-radical... */
                if ( at[i].radical && (at[i].radical!=RADICAL_SINGLET) ) 
                    continue;        

                /* Center found, now examine the adjacent terminals... */    
                {            
                    int en_term, kk=0, jj=0, min_en=999, iso=0, min_iso=999;
                    jj = -1;
                    for (k=0; k<at[i].valence; k++)
                    {
                        j = at[i].neighbor[k];
                
                        /* Terminal should be: ... */ 

                        /* terminal... */
                        if ( at[j].valence!=1 ) 
                            continue;
                        /* double-bonded ... */
                        if ( at[i].bond_type[k] != BOND_TYPE_DOUBLE ) 
                            continue;
                        /* zero-charged ... */
                        if ( at[j].charge != 0 ) 
                            continue;
                        /* non-radical */
                        if ( at[j].radical && (at[j].radical!=RADICAL_SINGLET) ) 
                            continue;
                        /*  of eligible elements list ... */
                        en_term = at[j].el_number;
                        switch (en_term)
                        {
                        case dO:    break;
                        case dS:    if ( (en_center==dSe)||(en_center==dAt)||(en_center==dTe) ) break;  continue;
                        case dSe:   if ( (en_center==dAt)||(en_center==dTe) ) break;  continue;
                        case dTe:   if ( en_center==dAt ) break; continue;
                        default:    continue;
                        }
                
                        /* From several candidates, select one with less el. number (==more electronegative). */
                        if ( en_term < min_en)  
                        { 
                            min_en = en_term; kk = k; jj = j; 
                            min_iso = at[j].iso_atw_diff > 0? at[i].iso_atw_diff-1 : at[i].iso_atw_diff;
                            continue; 
                        }   
                        /* From same-element candidates, select one with less isotopic mass (arbitrary choice). */
                        else if (en_term==min_en)
                        {
                            iso = at[j].iso_atw_diff > 0? at[i].iso_atw_diff-1 : at[i].iso_atw_diff;
                            if ( iso <min_iso )
                            { 
                                min_iso = iso; kk = k; jj = j; continue; 
                            }
                        }
                    } /* end of checking nbrs. */
            
                    /* If OK, apply changes. */
                    if (jj>=0)
                    {
                        at[i].charge = 0; 
                        at[jj].charge = -1;
                        at[i].bond_type[kk] = BOND_TYPE_SINGLE; 
                        at[jj].bond_type[0] = BOND_TYPE_SINGLE;
                        at[i].bond_stereo[kk] = at[jj].bond_stereo[0] = 0;
                        at[i].chem_bonds_valence--; 
                        at[jj].chem_bonds_valence--;
                        num_changes ++;
                    }
                }  

            }  /* end of search for candidate centers. */    

        } /* end of correcting oxoanions */




        /* Correct non-uniformly drawn amidinium cations. */
        {   

            /* Amidines include carboxamidines RC(=NR)NR2, 
            sulfinamidines RS(=NR)NR2 and phosphinamidines, R2P(=NR)NR2.
            
            
            NR                              NR
            |                               |
            R"-Y-NHR'         ===>          R"-Y=N(+)HR'        
            (+)                              

    
            Y = C, S, P
                                 

            Fix:            
            move positive charge to nitrogen and change single bond to a double one.

          
            Comment: 
            Fix is applied only if at least one of R's at N is hydrogen 
            (otherwise we have just a '+' delocalization which is already recognized).


            */
 
            static U_CHAR  allowed_elnums_center[] = {dC, dS, dP} ; 
            int en_center;
            int i, j, k, jj, kk; 
            int mismatch = 0, nuH=0, nuN = 0, nitrogens[MAXVAL];

            for (i=0; i<num_atoms; i++)
            {
                /* Find appropriate central atom. This center should be: ...*/ 
        
                /* charged exactly (+1) ... */
                if ( at[i].charge != 1 ) 
                    continue;
                en_center = at[i].el_number;
                /*  from eligible element list ... */
                if ( !memchr(allowed_elnums_center, en_center , sizeof(allowed_elnums_center)) ) 
                    continue;
                /* has exactly 3 neighbours connected by single bonds*/
                if ( at[i].valence != 3 ) 
                    continue;
                if ( at[i].chem_bonds_valence != 3 ) 
                    continue;

                /* non-radical. */
                if ( at[i].radical && (at[i].radical!=RADICAL_SINGLET) ) 
                    continue;        
        
                /* NB: center must have neutral neighbours, two of them are aliphatic N's of which at least one bears H. */
                mismatch = nuH = nuN = jj = kk = 0;
                memset(nitrogens, 0, sizeof(nitrogens));           
                jj = -1;
                for (k=0; k<at[i].valence; k++)
                {            
                    j = at[i].neighbor[k];
                
                    if ( at[j].charge != 0 ) 
                    { 
                        mismatch = 1; 
                        break; 
                    }
                    if ( at[j].el_number == dN) 
                    {
                        if ( ( at[j].valence > 3 ) || ( at[j].chem_bonds_valence > 3 ) )
                        { 
                            mismatch = 1; 
                            break; 
                        }
                        nuH+= NUMH(at,j);
                        nuN++;
                        if (jj<0) 
                        {
                            jj = j;               
                            kk = k;
                        }
                    }
                }

                /* If OK, apply changes. */
                if (mismatch)   continue;
                if (nuN!=2)     continue;
                if (nuH<1)      continue;
                if (jj>=0)
                {
                    at[i].charge = 0; 
                    at[jj].charge = 1;
                    at[i].bond_type[kk] = BOND_TYPE_DOUBLE; 
                    for ( k1 = 0; k1 < at[jj].valence && i != at[jj].neighbor[k1]; k1 ++ ) 
                        ;
                    at[jj].bond_type[k1] = BOND_TYPE_DOUBLE; 
                    at[i].chem_bonds_valence++; 
                    at[jj].chem_bonds_valence++;
                    /* NB: do nothing with wedge stereo bonds (retain wedge) */
                    num_changes ++;
                }

            }  /* end of search for candidate centers. */    
    
        } /* end of correcting amidiniums */



    } /*( if (bFixNonUniformDraw) */




    if ( !ne ) 
    { 
        /* one time initialization */
        const char *b, *e;
        int  len;
        for ( b = el; (e = strchr( b, ';')); b = e+1 )
        {
            len = e-b;
            memcpy( elname, b, len );
            elname[len] = '\0';
            en[ne++] = get_periodic_table_number( elname );
        }
        ne2 = ne;
        el_number_P  = get_periodic_table_number( "P" );
        el_number_H  = get_periodic_table_number( "H" );
        el_number_O  = get_periodic_table_number( "O" );
        en[ne2++] = el_number_C  = get_periodic_table_number( "C" );
        en[ne2++] = el_number_Si = get_periodic_table_number( "Si" );
    }
    
    /* H(-)-X  -> H-X(-);  H(+)-X  -> H-X(+) */
    for ( i1 = 0; i1 < num_atoms; i1 ++ ) 
    {
        if ( 1 == at[i1].valence &&
             1 == abs(at[i1].charge)  &&
             (0 == at[i1].radical || RADICAL_SINGLET == at[i1].radical) &&
             BOND_TYPE_SINGLE == at[i1].bond_type[0] &&
             el_number_H == at[i1].el_number &&
             el_number_H != at[i2=(int)at[i1].neighbor[0]].el_number &&
             !NUMH(at,i1) &&
             !NUMH(at,i2)
             ) 
        {
            at[i2].charge += at[i1].charge;
            at[i1].charge  = 0;
        }
    }

    /* replace XHm(-)--Y==XHn(+) with XHm==Y--XHn, (n>=0 ,m>=0, X=N,P,As,Sb,O,S,Se,Te) */
    for ( i1 = 0; i1 < num_atoms; i1 ++ ) 
    {
        if ( 1 != at[i1].charge ||
             (at[i1].radical && RADICAL_SINGLET != at[i1].radical) ||
             at[i1].chem_bonds_valence == at[i1].valence ||
             !memchr(en, at[i1].el_number, ne) ||
             get_el_valence( at[i1].el_number, at[i1].charge, 0 ) != at[i1].chem_bonds_valence+NUMH(at,i1) ) 
        {
            continue;
        }
        
        /* found a candidate at[i1] for X in XHn(+) */
        if ( 1 == at[i1].valence &&
             BOND_TYPE_DOUBLE == at[i1].bond_type[0] ) 
        {
            c = (int)at[i1].neighbor[0];
            for ( k2 = 0; k2 < at[c].valence; k2 ++ ) 
            {
                i2 = at[c].neighbor[k2];
                if ( 1 == at[i2].valence &&
                    -1 == at[i2].charge  &&
                     at[i2].el_number == at[i1].el_number && /* exact match */
                     (0 == at[i2].radical || RADICAL_SINGLET == at[i2].radical) &&
                     BOND_TYPE_SINGLE == at[i2].bond_type[0] &&
                     /*memchr(en, at[i2].el_number, ne) &&*/
                     get_el_valence( at[i2].el_number, at[i2].charge, 0 ) == at[i2].chem_bonds_valence+NUMH(at,i2) ) {
                    /* found both X(-) and X(+); change bonds and remove charges */
                    for ( k1 = 0; k1 < at[c].valence && i1 != at[c].neighbor[k1]; k1 ++ )
                        ;
                    at[i1].charge = at[i2].charge = 0;
                    at[i1].bond_type[0] = at[c].bond_type[k1] = BOND_TYPE_SINGLE;
                    at[i1].chem_bonds_valence --;
                    at[i2].bond_type[0] = at[c].bond_type[k2] = BOND_TYPE_DOUBLE;
                    at[i2].chem_bonds_valence ++;
                    num_changes ++;
                    break;
                }
            }
        }
        else 
        {
            /* explicit H case: detect H-neighbors and Y */
            int ineigh, neigh, i1_c, i2_c, num_H_i1, num_H_i2;
            for ( ineigh = 0, num_H_i1 = 0, i1_c = -1; ineigh < at[i1].valence; ineigh ++ ) 
            {
                neigh = at[i1].neighbor[ineigh];
                if ( at[neigh].el_number == el_number_H ) 
                {
                    if ( at[neigh].chem_bonds_valence == 1  &&
                         (0 == at[neigh].radical || RADICAL_SINGLET == at[neigh].radical) ) 
                    {
                        num_H_i1 ++; /* found H-neighbor */
                    } 
                    else 
                    {
                        break;  /* wrong neighbor */
                    }
                } 
                else if ( at[i1].bond_type[ineigh] == BOND_TYPE_DOUBLE ) 
                {
                    /* found a candidate for Y; bond must be double */
                    i1_c = ineigh;
                    c    = neigh;
                }
            }
            if ( i1_c < 0 || num_H_i1 + 1 != at[i1].valence ) 
            {
                continue;
            }
            for ( k2 = 0; k2 < at[c].valence; k2 ++ ) 
            {
                i2 = at[c].neighbor[k2];
                if (-1 == at[i2].charge  &&
                     at[i2].el_number == at[i1].el_number && /* exact match */
                     (0 == at[i2].radical || RADICAL_SINGLET == at[i2].radical) &&
                     get_el_valence( at[i2].el_number, at[i2].charge, 0 ) == at[i2].chem_bonds_valence+NUMH(at,i2) ) 
                {
                    for ( ineigh = 0, num_H_i2 = 0, i2_c = -1; ineigh < at[i2].valence; ineigh ++ ) 
                    {
                        neigh = at[i2].neighbor[ineigh];
                        if ( at[neigh].el_number == el_number_H ) 
                        {
                            if ( at[neigh].chem_bonds_valence == 1  &&
                                 (0 == at[neigh].radical || RADICAL_SINGLET == at[neigh].radical) ) 
                            {
                                num_H_i2 ++;  /* found H-neighbor */
                            } 
                            else 
                            {
                                break; /* wrong neighbor */
                            }
                        } 
                        else
                        if ( c == neigh && at[i2].bond_type[ineigh] == BOND_TYPE_SINGLE ) 
                        { 
                            i2_c = ineigh; /* position of Y neighbor; bond must be single */
                        } 
                        else 
                        {
                            break;
                        }
                    }
                    if ( num_H_i2 + (i2_c >= 0) != at[i2].valence ) 
                    {
                        continue;
                    }
                    /* found both X(-) and X(+); change bonds and remove charges */
                    for ( k1 = 0; k1 < at[c].valence && i1 != at[c].neighbor[k1]; k1 ++ )
                        ;
                    at[i1].charge = at[i2].charge = 0;
                    at[i1].bond_type[i1_c] = at[c].bond_type[k1] = BOND_TYPE_SINGLE;
                    at[i1].chem_bonds_valence --;
                    at[i2].bond_type[i2_c] = at[c].bond_type[k2] = BOND_TYPE_DOUBLE;
                    at[i2].chem_bonds_valence ++;
                    num_changes ++;
                    break;
                }
            }
        }
    }

    /* Replace
                                                                             
          X-                X        X=O,S,Se,Te -- terminal atoms  (NEIGHB2)         
        \ |               \ ||                                               
         >Y++    with      >Y        Y=S,Se,Te   -- central cation  (CENTER2)          
        / |               / ||                                               
          X-                X        Y valence=4, original Y bond valence = 4


       --- the following case of P is processed separately in remove_ion_pairs()
       --- therefire, it has been disabled here, see #ifndef FIX_P_IV_Plus_O_Minus -- 2010-03-17 DT

          X-                X        X=O,S,Se,Te -- terminal atoms  (NEIGHB2)         
        \ |               \ ||                                               
         >P+     with      >P                 
        / |               / |                                                
          X-                X-       Y valence=4, original Y bond valence = 4

    */

    for ( i1 = 0; i1 < num_atoms; i1 ++ ) 
    {
        if ( 1 == at[i1].valence &&
             -1 == at[i1].charge  &&
             (0 == at[i1].radical || RADICAL_SINGLET == at[i1].radical) &&
              !NUMH(at,i1) &&
             BOND_TYPE_SINGLE == at[i1].bond_type[0] &&
             memchr( en+FIRST_NEIGHB2, at[i1].el_number, ne-FIRST_NEIGHB2 ) ) 
        {
            int charge, i;
            /* found a candidate for X */
            c = (int)at[i1].neighbor[0]; /* candidate for Y */
            if ( ((charge=2) == at[c].charge && memchr( en+FIRST_CENTER2, at[c].el_number, ne-FIRST_CENTER2) 
#ifndef FIX_P_IV_Plus_O_Minus                 
                 || (charge=1) == at[c].charge && el_number_P==at[c].el_number
#endif
                 ) &&
                 4 == at[c].valence &&
                 (0 == at[c].radical || RADICAL_SINGLET == at[c].radical ) &&
                 at[c].valence == at[c].chem_bonds_valence &&
                 !NUMH(at,c) ) 
            {
                ;  /* accept */
            } 
            else 
            {
                continue; /* ignore at[i1] */
            }
            for ( k2 = 0; k2 < at[c].valence; k2 ++ ) 
            {
                i2 = at[c].neighbor[k2];
                if ( i2 == i1 ) 
                {
                    continue;
                }
                if ( 1 == at[i2].valence &&
                    -1 == at[i2].charge  &&
                     memchr( en+FIRST_NEIGHB2, at[i2].el_number, ne-FIRST_NEIGHB2 ) &&
                     /*at[i2].el_number == at[i1].el_number &&*/ /* exact match */
                     (0 == at[i2].radical || RADICAL_SINGLET == at[i2].radical) &&
                     !NUMH(at,i2) &&
                     BOND_TYPE_SINGLE == at[i2].bond_type[0]  ) 
                {
                    /* found both X(-) and X(-); change bonds and remove charges */
                    for ( k1 = 0; k1 < at[c].valence && i1 != at[c].neighbor[k1]; k1 ++ )
                        ;
                    for ( i = 0; i < charge; i ++ ) 
                    {
                        /* in case of P it does not matter which X atom is neutralized
                           because of tautomerism. However, neutral central atom is important
                           for the neutralization of the components */
                        switch ( i ) 
                        {
                        case 0:
                            at[i1].charge ++; /* = 0; changed 2010-03-17 DT*/
                            at[i1].bond_type[0] = at[c].bond_type[k1] = BOND_TYPE_DOUBLE;
                            at[i1].bond_stereo[0] = at[c].bond_stereo[k1] = 0;
                            at[i1].chem_bonds_valence ++;
                            at[c].chem_bonds_valence ++;
                            if ( bFixBug ) at[c].charge --; /* added 2010-03-17 DT*/
                            num_changes ++;
                            break;
                        case 1:
                            at[i2].charge ++; /*= 0; changed 2010-03-17 DT*/
                            at[i2].bond_type[0] = at[c].bond_type[k2] = BOND_TYPE_DOUBLE;
                            at[i2].bond_stereo[0] = at[c].bond_stereo[k2] = 0;
                            at[i2].chem_bonds_valence ++;
                            at[c].chem_bonds_valence ++;
                            if ( bFixBug ) at[c].charge --; /* added 2010-03-17 DT */
                            num_changes ++;
                            break;
                        }
                    }
/*   -- removed -- 2010-03-17 DT
#if ( FIX_ODD_THINGS_REM_Plus_BUG == 1 )
                    at[c].charge -= charge;
#else
                    if ( bFixBug ) 
                    {
                        at[c].charge -= charge;
                    }
#endif
*/
                    break;
                }
            }
        }
    }




    /* A(doublet)-B(doublet) -> A=B  (A and B have no other doublet neighbors) */
    /* A(doublet)=B(doublet) -> A#B  (A and B have no other doublet neighbors) */
    for( i1 = 0; i1 < num_atoms; i1 ++ ) 
    {
        if ( RADICAL_DOUBLET == at[i1].radical &&
             0 <= (i2=the_only_doublet_neigh(at, i1, &k1, &k2)) ) 
        {
            if ( at[i1].bond_type[k1] <= BOND_TYPE_DOUBLE ) 
            {
                at[i1].bond_type[k1] ++;
                at[i1].chem_bonds_valence ++;
                at[i2].bond_type[k2] ++;
                at[i2].chem_bonds_valence ++;
                at[i1].radical = 0;
                at[i2].radical = 0;
            }
        }
    }

#if ( REMOVE_ION_PAIRS_EARLY == 1 )
    num_changes += remove_ion_pairs( num_atoms, at );
#endif




    return num_changes;
}








/************************************************************************/
int post_fix_odd_things( int num_atoms, inp_ATOM *at )
{   
    int num_changes = 0;
    /* currently does nothing */
    return num_changes;
}



/************************************************************************/
int nFindOneOM(inp_ATOM *at, int at_no, int ord_OM[], int num_OM)
{
    int i, n_OM, n_OM_best, best_value, cur_value, diff;
    int num_best;

    if ( 1 == num_OM ) {
        return ord_OM[0];
    }
    if ( 1 > num_OM ) {
        return -1;
    }

    /* select neighbors with min. number of bonds */
    num_best = 1;
    n_OM       = (int)at[at_no].neighbor[ord_OM[0]];
    best_value = (int)at[n_OM].valence;
    /* compare number of bonds; move indexes of the best neighbors to the first elements of ord_OM[] */
    for ( i = 1; i < num_OM; i ++ ) {
        n_OM = at[at_no].neighbor[ord_OM[i]];
        cur_value = (int)at[n_OM].valence;
        diff = cur_value - best_value;
        if ( diff < 0 ) {
            n_OM_best  = n_OM;
            best_value = cur_value;
            ord_OM[0]  = ord_OM[i];
            num_best   = 1;
        } else
        if ( diff == 0 ) {  /* was '=', pointed by WDI */
            ord_OM[num_best ++] = ord_OM[i];
        }
    }
    num_OM = num_best;
    if ( 1 == num_OM ) {
        return ord_OM[0];
    }
    /* select neighbors with min. periodic numbers */
    num_best = 1;
    n_OM       = (int)at[at_no].neighbor[ord_OM[0]];
    best_value = (int)at[n_OM].el_number;
    /* compare periodic numbers; move indexes of the best neighbors to the first elements of ord_OM[] */
    for ( i = 1; i < num_OM; i ++ ) {
        n_OM = at[at_no].neighbor[ord_OM[i]];
        cur_value = (int)at[n_OM].el_number;
        diff = cur_value - best_value;
        if ( diff < 0 ) {
            n_OM_best  = n_OM;
            best_value = cur_value;
            ord_OM[0]  = ord_OM[i];
            num_best   = 1;
        } else
        if ( diff == 0 ) {  /* was '=', pointed by WDI */
            ord_OM[num_best ++] = ord_OM[i];
        }
    }
    num_OM = num_best;
    if ( 1 == num_OM ) {
        return ord_OM[0];
    }
    /* if neighbors are not terminal atoms then reject */
    if ( 1 < at[n_OM].valence ) {
        return -1;
    }
    /* if neighbors are terminal atoms then the one without isotope or with lightest isotope */
    num_best = 1;
    n_OM       = (int)at[at_no].neighbor[ord_OM[0]];
    best_value = (int)at[n_OM].iso_atw_diff;
    /* compare periodic numbers; move indexes of the best neighbors to the first elements of ord_OM[] */
    for ( i = 1; i < num_OM; i ++ ) {
        n_OM = at[at_no].neighbor[ord_OM[i]];
        cur_value = (int)at[n_OM].el_number;
        diff = cur_value - best_value;
        if ( (!cur_value && best_value) || diff < 0 ) {
            n_OM_best  = n_OM;
            best_value = cur_value;
            ord_OM[0]  = ord_OM[i];
            num_best   = 1;
        } else
        if ( diff == 0 ) {    /* was '=', pointed by WDI */
            ord_OM[num_best ++] = ord_OM[i];
        }
    }
    num_OM = num_best;
    if ( 1 == num_OM ) {
        return ord_OM[0];
    }
    /* return any */
    return ord_OM[0];
}



/************************************************************************/
/* the bonds are fixed in fix_special_bonds() */
int remove_ion_pairs( int num_atoms, inp_ATOM *at )
{   
    int num_changes = 0;

    /*                           0 1 2  3  4 5 6  7  8  9                   8  9  */
#if ( FIX_REM_ION_PAIRS_Si_BUG == 1 )
    static const char    el[] = "N;P;As;Sb;O;S;Se;Te;C;Si;";   /* 8 elements + C, Si */
#else
    static const char    el[] = "N;P;As;Sb;O;S;Se;Te;C;Si";   /* 8 elements + C, Si */
#endif
    static char    en[12];         /* same number: 8 elements */
    static int     ne=0;           /* will be 8 and 10 */

#define ELEM_N_FST  0
#define ELEM_N_LEN  4
#define ELEM_O_FST  4
#define ELEM_O_LEN  4
#define ELEM_C_FST  8
#define ELEM_C_LEN  2

#define MAX_NEIGH 6

    int i, n, n2, i1, i2, i3, i4, type, chrg;
    int num_C_II=0, num_C_plus=0, num_C_minus=0, num_N_plus=0, num_N_minus=0, num_O_plus=0, num_O_minus=0, num_All;
#ifdef FIX_P_IV_Plus_O_Minus
    int num_P_IV_plus=0; /* added 2010-03-17 DT */
#endif
    inp_ATOM *a;
    char elname[ATOM_EL_LEN], *p;
    if ( !ne ) { /* one time initialization */
        const char *b, *e;
        int  len;
        for ( b = el; (e = strchr( b, ';')); b = e+1 ) {
            len = e-b;
            memcpy( elname, b, len );
            elname[len] = '\0';
            en[ne++] = get_periodic_table_number( elname );
        }
        en[ne] = '\0';
    }

    /****** count candidates ********/
    for ( i = 0, a = at; i < num_atoms; i ++, a++ ) {
        if ( 1 == (chrg=a->charge) || -1 == chrg ) {
            if ( (p = (char*)memchr( en, a->el_number, ne)) ) {
                n = p - en;
                if ( n >= ELEM_C_FST ) {
                    if ( chrg > 0 )
                        num_C_plus ++;
                    else
                        num_C_minus ++;
                } else
                if ( n >= ELEM_O_FST ) {
                    if ( chrg > 0 )
                        num_O_plus ++;
                    else
                        num_O_minus ++;
                } else {
                    if ( chrg > 0 )
                        num_N_plus ++;
                    else
                        num_N_minus ++;
#ifdef FIX_P_IV_Plus_O_Minus
                    num_P_IV_plus += n > 0 && chrg == 1 && a->valence == 4 && a->chem_bonds_valence == 4; /* added 2010-03-17 DT */
#endif
                }
            }
        } else
        if ( !chrg && a->chem_bonds_valence + NUMH(a, 0) == 2 &&
             get_el_valence( a->el_number, 0, 0 ) == 4     &&
             NULL != memchr( en+ELEM_C_FST, a->el_number, ELEM_C_LEN) ) {
            num_C_II ++;
        }
    }
    num_All = num_C_II + num_C_plus + num_C_minus + num_N_plus + num_N_minus + num_O_plus + num_O_minus;
    /* do not add num_P_IV_plus ! -- 2010-03-17 DT */
    if ( !num_All ) {
        return 0;
    }

    /**************************************************************************/
    /*************************** Terminal ion pairs ***************************/
    /**************************************************************************/

    /*-------------------------------------------------------------------------
       Pair type 1            N=N,P,As,Sb; O=O,S,Se,Te
       ===========

         X              X     if X is another -O(-) then neutralize O(-) 
         |              |     that has the smallest periodic table number
       O=N(+)-O(-) => O=N=O
         i    n
     --------------------------------------------------------------------------*/
    for ( type = 1; type <= 18; type ++ ) {
        if ( (!type || 1 == type) ) {
            for ( i = 0; i < num_atoms && 0 < num_N_plus && 0 < num_O_minus; i ++ ) {
                if ( 1 == at[i].charge && 3 == nNoMetalNumBonds(at, i) &&
                     4 == nNoMetalBondsValence(at, i) &&
                     NULL != memchr( en+ELEM_N_FST, at[i].el_number, ELEM_N_LEN) ) {
                    int num_OM = 0, ord_OM[3]; /* -O(-) */
                    int num_O  = 0; /* =O    */
                    int num_O_other = 0;
                    for ( i1 = 0; i1 < at[i].valence; i1 ++ ) {
                        n = at[i].neighbor[i1];
                        if ( 1 == nNoMetalNumBonds(at, n) && 0 == num_of_H( at, n ) &&
                             NULL != (p = (char*)memchr( en+ELEM_O_FST, at[n].el_number, ELEM_O_LEN)) ) {
                            if ( BOND_TYPE_SINGLE == at[i].bond_type[i1] &&
                                 -1               == at[n].charge       ) {
                                ord_OM[num_OM ++]  = i1;
                            } else
                            if ( BOND_TYPE_DOUBLE == at[n].bond_type[0] &&
                                 0                == at[n].charge       ) {
                                num_O ++;
                            } else {
                                num_O_other ++;
                            }
                        }
                    }
                    if ( num_OM > 0 && num_O > 0 && !num_O_other &&
                         0 <= (i1=nFindOneOM(at, i, ord_OM, num_OM)) ) {
                        /* remove charges and increase bond order */
                        n = at[i].neighbor[i1];
                        i2 =  is_in_the_list( at[n].neighbor, (AT_NUMB)i, at[n].valence ) - at[n].neighbor;
                        at[i].bond_type[i1] ++;
                        at[n].bond_type[i2] ++;
                        at[i].chem_bonds_valence ++;
                        at[n].chem_bonds_valence ++;
                        at[i].charge --;
                        at[n].charge ++;
                        at[i].radical = 0;
                        at[n].radical = 0;
                        num_changes ++;
                        num_N_plus --;
                        num_O_minus --;
                        num_All -= 2;
                    }
                }
            }
#ifdef FIX_P_IV_Plus_O_Minus
    /*-------------------------------------------------------------------------
       Pair type 1a           P=P,As,Sb; O=O,S,Se,Te  -- added 2010-03-17
       =============

         X              X     if X, Y, or Z is another -O(-) then neutralize O(-) 
         |              |     that has the smallest periodic table number
       Y-P(+)-O(-) => Y-P=O
         |i   n         |
         Z              Z
         
     --------------------------------------------------------------------------*/
            for ( i = 0; i < num_atoms && 0 < num_P_IV_plus /*&& 0 < num_N_plus*/ && 0 < num_O_minus; i ++ ) {
                if ( 1 == at[i].charge && 4 == nNoMetalNumBonds(at, i) &&
                     4 == nNoMetalBondsValence(at, i) &&
                     NULL != memchr( en+ELEM_N_FST+1, at[i].el_number, ELEM_N_LEN-1) ) {
                    int num_OM = 0, ord_OM[4]; /* -O(-) */
                    /*int num_O  = 0;*/ /* =O    */
                    int num_O_other = 0;
                    for ( i1 = 0; i1 < at[i].valence; i1 ++ ) {
                        n = at[i].neighbor[i1];
                        if ( 1 == nNoMetalNumBonds(at, n) && 0 == num_of_H( at, n ) &&
                             NULL != (p = (char*)memchr( en+ELEM_O_FST, at[n].el_number, ELEM_O_LEN)) ) {
                            if ( BOND_TYPE_SINGLE == at[i].bond_type[i1] &&
                                 -1               == at[n].charge       ) {
                                ord_OM[num_OM ++]  = i1;
                            /*
                            }
                            if ( BOND_TYPE_DOUBLE == at[n].bond_type[0] &&
                                 0                == at[n].charge       ) {
                                num_O ++;
                            */
                            } else {
                                num_O_other ++;
                            }
                        }
                    }
                    if ( num_OM > 0 /*&& num_O > 0 && !num_O_other*/ &&
                         0 <= (i1=nFindOneOM(at, i, ord_OM, num_OM)) ) {
                        /* remove charges and increase bond order */
                        n = at[i].neighbor[i1];
                        i2 =  is_in_the_list( at[n].neighbor, (AT_NUMB)i, at[n].valence ) - at[n].neighbor;
                        at[i].bond_type[i1] ++;
                        at[n].bond_type[i2] ++;
                        at[i].chem_bonds_valence ++;
                        at[n].chem_bonds_valence ++;
                        at[i].charge --;
                        at[n].charge ++;
                        at[i].radical = 0;
                        at[n].radical = 0;
                        num_changes ++;
                        num_N_plus --;
                        num_O_minus --;
                        num_P_IV_plus --;
                        num_All -= 2;
                    }
                }
            }
#endif /* FIX_P_IV_Plus_O_Minus */
        }
        /*-------------------------------------------------------------------------
           Terminal pair types: 2,3,4,5,6,7,8,9   N=N,P,As,Sb; O=O,S,Se,Te; C=C,Si
           ====================================
    type #
       2 2:  O=N-C(II)-     => O=N#C-      N=N,P,As,Sb; O=O,S,Se,Te; C=C,Si
       3 9:  O=O(+)-C(-)(III) => O=O=C(IV)
       4 3:  O(-)-N(+)(IV)  => O=N(V)  (input structure has at least 1 double bond)
       5 4:  O(-)-O(+)(III) => O=O(IV)
       6 8:  O(-)-O-C(+)(III) => O=O=C(IV)
       7 5:  N(-)=N(+)(IV)  => N#N(V)    allow terminal H on N(-)
       8 6:  N(-)=O(+)(III) => N#O-
       9 7:  N(-)=C(+)(III) => N#C-
         --------------------------------------------------------------------------*/
        if ( !type || (2 <= type && type <= 9) ) {
            for ( i = 0; i < num_atoms && 0 < num_All; i ++ ) {
                if ( 0 == at[i].charge && 1 == nNoMetalNumBonds(at, i) && 2 == nNoMetalBondsValence(at, i) &&
                     0 == num_of_H( at, i ) &&
                     NULL != memchr( en+ELEM_O_FST, at[i].el_number, ELEM_O_LEN) &&
                     0 <= ( i1 = nNoMetalNeighIndex(at, i)) &&
                     at[i].bond_type[i1] <= BOND_TYPE_TRIPLE ) {
                    /* terminal O= */
                    n = at[i].neighbor[i1];
                    if ( (!type || type == 2) && 0 < num_C_II ) {   /* avoid alternating bonds */
                        if ( 0 == at[n].charge &&
                             2 == nNoMetalNumBonds(at, n) && 3 == nNoMetalBondsValence(at, n) &&
                             0 == num_of_H( at, n ) &&
                             NULL != memchr( en+ELEM_N_FST, at[n].el_number, ELEM_N_LEN) &&
                             0 <= (i2 = nNoMetalOtherNeighIndex( at, n, i ) ) &&
                             at[n].bond_type[i2] <= BOND_TYPE_TRIPLE ) {
                            /* i2 = index of opposite to at[i] neighbor of at[n] */
                            /*i2 = (at[n].neighbor[0] == i);*/
                            n2 = at[n].neighbor[i2];
                            if ( 0 == at[n2].charge &&
                                 2 == at[n2].valence && 2 == at[n2].chem_bonds_valence &&
                                 0 == num_of_H( at, n2 ) &&
                                 NULL != memchr( en+ELEM_C_FST, at[n2].el_number, ELEM_C_LEN) ) {
                                 /*       i n n2     */
                                 /* found O=N-C(II)- */
                                 /* convert O=N-C(II)-     => O=N#C- */
                                i3 = (at[n2].neighbor[0] != n); /* index of at[n] neighbor of n2 */
                                at[ n].chem_bonds_valence = 5; /* N */
                                at[n2].chem_bonds_valence = 4; /* C */
                                at[ n].bond_type[i2] = BOND_TYPE_TRIPLE;
                                at[n2].bond_type[i3] = BOND_TYPE_TRIPLE;
                                at[n2].radical = 0;
                                num_changes ++;
                                num_C_II --;
                                num_All --;
                                continue;
                            }
                        }
                    }
                    if ( (!type || type == 3) && 0 < num_O_plus && 0 < num_C_minus ) {
                        if ( 1 == at[n].charge && 2 == nNoMetalNumBonds(at, n) && 3 == nNoMetalBondsValence(at, n) &&
                             0 == num_of_H( at, n ) &&
                             NULL != memchr( en+ELEM_O_FST, at[n].el_number, ELEM_O_LEN)  &&
                             0 <= (i2 = nNoMetalOtherNeighIndex( at, n, i ) ) &&
                             at[n].bond_type[i2] <= BOND_TYPE_TRIPLE ) {
                            /* found O=O(+)- */
                            /* i2 = index of opposite to at[i] neighbor of at[n] */
                            /*i2 = (at[n].neighbor[0] == i);*/
                            n2 = at[n].neighbor[i2];
                            if ( -1 == at[n2].charge && 3 >= nNoMetalNumBonds(at, n2) && 3 == nNoMetalBondsValence(at, n2)+NUMH(at,n2) &&
                                 NULL != memchr( en+ELEM_C_FST, at[n2].el_number, ELEM_C_LEN) ) {
                                 /*             i n    n2        */
                                 /* found found O=O(+)-C(-)(III) */
                                 /* convert O=O(+)-C(-)(III)     => O=O=C(IV) */
                                i3 = (at[n2].neighbor[0] != n); /* index of at[n] neighbor of n2 */
                                at[ n].charge --;
                                at[n2].charge ++;
                                at[ n].chem_bonds_valence += 1; /* =O- => =O= */
                                at[n2].chem_bonds_valence += 1; /* -C  => =C  */
                                at[ n].bond_type[i2] = BOND_TYPE_DOUBLE;
                                at[n2].bond_type[i3] = BOND_TYPE_DOUBLE;
                                num_changes ++;
                                num_O_plus --;
                                num_C_minus --;
                                num_All -= 2;
                                continue;
                            }
                        }
                    }
                } else
                if ( -1 == at[i].charge &&
                      0 < num_O_minus + num_N_minus &&
                      0 < num_N_plus + num_O_plus + num_C_plus &&
                      1 == nNoMetalNumBonds(at, i) && 1 == nNoMetalBondsValence(at, i) &&
                      0 == num_of_H( at, i ) &&
                     NULL != memchr( en+ELEM_O_FST, at[i].el_number, ELEM_O_LEN) &&
                     0 <= (i1 = nNoMetalNeighIndex( at, i )) &&
                     at[i].bond_type[i1] <= BOND_TYPE_TRIPLE ) {
                    
                    /* terminal O(-)- */
                    
                    n = at[i].neighbor[i1];
                    
                    if ( (!type || type == 4) && 0 < num_O_minus && 0 < num_N_plus && /* O(-)-N(+)(IV) */
                         1 == at[n].charge && 3 >= nNoMetalNumBonds(at, n) && 4 == nNoMetalBondsValence(at, n) &&
                         0 == num_of_H( at, n ) &&
                         NULL != memchr( en+ELEM_N_FST, at[n].el_number, ELEM_N_LEN) /* except >O(+)- */
                       ) {
                         /* found O(-)-N(+)(IV) */
                         /* convert O(-)-N(+)(IV)     => O=N(V)  */

                        i2 =  is_in_the_list( at[n].neighbor, (AT_NUMB)i, at[n].valence ) - at[n].neighbor; /* index of at[i] neighbor of at[n] */
                        at[i].charge ++;
                        at[n].charge --;
                        at[i].chem_bonds_valence ++;
                        at[n].chem_bonds_valence ++;
                        at[i].bond_type[i1] ++;
                        at[n].bond_type[i2] ++;
                        num_changes ++;
                        num_O_minus --;
                        num_N_plus --;
                        num_All -= 2;
                        continue;
                    }
                    
                    if ( (!type || type == 5) && 0 < num_O_minus && 0 < num_O_plus &&/* O(-)-O(+)(III) */
                         1 == at[n].charge && 3 >= nNoMetalNumBonds(at, n) && 3 == nNoMetalBondsValence(at, n) &&
                         0 == num_of_H( at, n ) &&
                         NULL != memchr( en+ELEM_O_FST, at[n].el_number, ELEM_O_LEN) /* except >O(+)- */
                       ) {
                         /* found  O(+)(III) */
                         /* convert O(-)-O(+)(III)    => O=O(IV) */

                        i2 =  is_in_the_list( at[n].neighbor, (AT_NUMB)i, at[n].valence ) - at[n].neighbor; /* index of at[i] neighbor of at[n] */
                        at[i].charge ++;
                        at[n].charge --;
                        at[i].chem_bonds_valence ++;
                        at[n].chem_bonds_valence ++;
                        at[i].bond_type[i1] ++;
                        at[n].bond_type[i2] ++;
                        num_changes ++;
                        num_O_minus --;
                        num_O_plus --;
                        num_All -= 2;
                        continue;
                    }
                                                 /* i    n n2        */
                    if ( (!type || type == 6) && /* O(-)-O-C(+)(III) */
                         0 < num_O_minus && 0 < num_C_plus &&
                         0 == at[n].charge && 2 == nNoMetalNumBonds(at, n) && 2 == nNoMetalBondsValence(at, n) &&
                         0 == num_of_H( at, n ) &&
                         NULL != memchr( en+ELEM_O_FST, at[n].el_number, ELEM_O_LEN) &&
                         0 <= (i2=nNoMetalOtherNeighIndex( at, n, i )) &&
                         at[n].bond_type[i2] <= BOND_TYPE_TRIPLE ) {
                        /* found O(-)-O- */
                        /* i2 = index of opposite to at[i] neighbor of at[n] */
                        /*i2 = (at[n].neighbor[0] == i);*/
                        n2 = at[n].neighbor[i2];
                        if ( 1 == at[n2].charge && 3 >= nNoMetalNumBonds(at, n2) &&
                             3 == nNoMetalBondsValence(at, n2)+NUMH(at,n2) &&
                             NULL != memchr( en+ELEM_C_FST, at[n2].el_number, ELEM_C_LEN) ) {
                            /*       i    n n2  */
                            /* found O(-)-O-C(+)(III) */
                            /* convert O(-)-O-C(+)(III)     => O=O=C(IV) */
                            /*i3 = (at[n2].neighbor[0] != n);*/ /* i3 = index of at[n] neighbor of at[n2] */
                            i3 = is_in_the_list( at[n2].neighbor, (AT_NUMB)n, at[n2].valence ) - at[n2].neighbor;
                            /*i4 = index of at[i] in the adjacency list of at[n] */
                            i4 = is_in_the_list( at[n].neighbor, (AT_NUMB)i, at[n].valence ) - at[n].neighbor;
                            at[ i].charge ++;
                            at[n2].charge --;
                            at[ i].chem_bonds_valence += 1; /* O-  => O=  */
                            at[ n].chem_bonds_valence += 2; /* -O- => =O= */
                            at[n2].chem_bonds_valence += 1; /* -C  => =C  */
                            at[ i].bond_type[i1]   = BOND_TYPE_DOUBLE;
                            at[ n].bond_type[i4]   = BOND_TYPE_DOUBLE;
                            at[ n].bond_type[i2]   = BOND_TYPE_DOUBLE;
                            at[n2].bond_type[i3]   = BOND_TYPE_DOUBLE;
                            num_changes ++;
                            num_O_minus --;
                            num_C_plus --;
                            num_All -= 2;
                            continue;
                        }
                    }
                } else
                if ( -1 == at[i].charge && 0 < num_N_minus && 0 < num_N_plus+num_O_plus+num_C_plus &&
                      1 == nNoMetalNumBonds(at, i) && 2 == nNoMetalBondsValence(at, i)+NUMH(at, i) &&
                      /*0 == num_of_H( at, i ) &&*/
                     NULL != memchr( en+ELEM_N_FST, at[i].el_number, ELEM_N_LEN) &&
                     0 <= (i1 = nNoMetalNeighIndex( at, i )) &&
                     at[i].bond_type[i1] <= BOND_TYPE_TRIPLE ) {
                    /* terminal N(-)= */
                    n = at[i].neighbor[i1 = 0];
                    if ( (!type || type == 7) && 0 < num_N_plus && /* N(-)=N(+)(IV) */
                         1 == at[n].charge && 3 >= nNoMetalNumBonds(at, n) && 4 == nNoMetalBondsValence(at, n) &&
                         0 == num_of_H( at, n ) &&
                         NULL != memchr( en+ELEM_N_FST, at[n].el_number, ELEM_N_LEN)
                       ) {
                         /* found N(-)-N(+)(IV) */
                         /* convert N(-)=N(+)(IV)     => N#N(V)  */

                        i2 =  is_in_the_list( at[n].neighbor, (AT_NUMB)i, at[n].valence ) - at[n].neighbor; /* index of at[i] neighbor of at[n] */
                        at[i].charge ++;
                        at[n].charge --;
                        at[i].chem_bonds_valence ++;
                        at[n].chem_bonds_valence ++;
                        at[i].bond_type[i1] ++;
                        at[n].bond_type[i2] ++;
                        num_changes ++;
                        num_N_minus --;
                        num_N_plus --;
                        num_All -= 2;
                        continue;
                    }
                    if ( (!type || type == 8) && 0 < num_O_plus && /* N(-)=O(+)(III) */
                         1 == at[n].charge && 2 == nNoMetalNumBonds(at, n) && 3 == nNoMetalBondsValence(at, n) &&
                         0 == num_of_H( at, n ) &&
                         NULL != memchr( en+ELEM_O_FST, at[n].el_number, ELEM_O_LEN)
                         ) {
                         /* found N(-)-O(+)(III) */
                         /* convert N(-)=O(+)(III)    => N#O(IV)- */
                        i2 =  is_in_the_list( at[n].neighbor, (AT_NUMB)i, at[n].valence ) - at[n].neighbor; /* index of at[i] neighbor of at[n] */
                        at[i].charge ++;
                        at[n].charge --;
                        at[i].chem_bonds_valence ++;
                        at[n].chem_bonds_valence ++;
                        at[i].bond_type[i1] ++;
                        at[n].bond_type[i2] ++;
                        num_changes ++;
                        num_N_minus --;
                        num_O_plus --;
                        num_All -= 2;
                        continue;
                    }
                    if ( (!type || type == 9) && 0 < num_C_plus && /* N(-)=C(+)(III) */
                         1 == at[n].charge && 2 == at[n].valence && 3 == at[n].chem_bonds_valence &&
                         0 == num_of_H( at, n ) &&
                         NULL != memchr( en+ELEM_C_FST, at[n].el_number, ELEM_C_LEN)
                         ) {
                         /* found N(-)=C(+)(III) */
                         /* convert N(-)=C(+)(III)    => N#C(IV)- */

                        i2 =  is_in_the_list( at[n].neighbor, (AT_NUMB)i, at[n].valence ) - at[n].neighbor; /* index of at[i] neighbor of at[n] */
                        at[i].charge ++;
                        at[n].charge --;
                        at[i].chem_bonds_valence ++;
                        at[n].chem_bonds_valence ++;
                        at[i].bond_type[i1] ++;
                        at[n].bond_type[i2] ++;
                        num_changes ++;
                        num_N_minus --;
                        num_C_plus --;
                        num_All -= 2;
                        continue;
                    }
                }
            }
        }

        /**************************************************************************/
        /*********************** NON-Terminal ion pairs ***************************/
        /**************************************************************************/
        /*-------------------------------------------------------------------------
           Non-Terminal pair types: 10,11,12,13,14   N=N,P,As,Sb; O=O,S,Se,Te; C=C,Si
           ========================================

        10:  N(+)(IV)-C(-)(III)     => N(V)=C(IV)  (N has 3 or 2 bonds)
        11:  N(+)(IV)=C(-)(III)     => N(V)#C(IV)  (N has 3 or 2 bonds)
        12:  N(+)(IV)-N(-)(II)      => N(V)=N(III) (allow terminal H on N(-))
        13: -O(+)-C(-)(III)         => -O=C-
        14: -O(+)=C(-)(III)         => -O#C-
        15:  O(+)(III)-N(-)(II)     => O(IV)=N(III) (allow terminal H on N(-))
         --------------------------------------------------------------------------*/
        if ( !type || (10 <= type && type <= 15) ) {
            for ( i = 0; i < num_atoms && 0 < num_All; i ++ ) {
                if ( 1 == at[i].charge &&
                     0 < num_N_plus + num_O_plus && 0 < num_C_minus + num_N_minus &&
                     4 >= nNoMetalNumBonds(at, i) && 4 == nNoMetalBondsValence(at, i) &&
                     0 == num_of_H( at, i ) &&
                     NULL != memchr( en+ELEM_N_FST, at[i].el_number, ELEM_N_LEN) ) {
                    /* found non-terminal N(+)(IV) */
                    if ( (!type || 10 == type) && 0 < num_N_plus && 0 < num_C_minus ) {
                        int num_neigh = 0, pos_neigh = -1;
                        for ( i1 = 0; i1 < at[i].valence; i1 ++ ) {
                            n = at[i].neighbor[i1];
                            if ( -1 == at[n].charge && 3 >= at[n].valence && 3 == at[n].chem_bonds_valence+NUMH(at,n) &&
                                  /*0 == at[n].num_H &&*/
                                  at[i].bond_type[i1] == BOND_TYPE_SINGLE &&
                                  NULL != memchr( en+ELEM_C_FST, at[n].el_number, ELEM_C_LEN) ) {
                                /* found N(+)(IV)-C(-)(III); prepare conversion to N(V)=C(IV) */
                                num_neigh ++;
                                pos_neigh = i1;
                            }
                        }
                        i1=pos_neigh;
                        if ( 1 == num_neigh &&
                             at[i].bond_type[i1] <= BOND_TYPE_TRIPLE &&
                             !has_other_ion_neigh( at, i, n=at[i].neighbor[i1], en, ne ) &&
                             !has_other_ion_neigh( at, n, i, en, ne )) {
                            /*n = at[i].neighbor[i1=pos_neigh];*/
                            i2 = is_in_the_list( at[n].neighbor, (AT_NUMB)i, at[n].valence ) - at[n].neighbor;
                            at[i].charge --;
                            at[n].charge ++;
                            at[i].chem_bonds_valence ++;
                            at[n].chem_bonds_valence ++;
                            at[i].bond_type[i1] ++;
                            at[n].bond_type[i2] ++;
                            num_changes ++;
                            num_C_minus --;
                            num_N_plus --;
                            num_All -= 2;
                            continue;
                        }
                    }
                    if ( (!type || 11 == type) && 0 < num_N_plus && 0 < num_C_minus ) {
                        int num_neigh = 0, pos_neigh = -1;
                        for ( i1 = 0; i1 < at[i].valence; i1 ++ ) {
                            n = at[i].neighbor[i1];
                            if ( -1 == at[n].charge && 3 >= at[n].valence && 3 == at[n].chem_bonds_valence+NUMH(at,n) &&
                                  /*0 == at[n].num_H &&*/
                                  at[i].bond_type[i1] == BOND_TYPE_DOUBLE &&
                                  NULL != memchr( en+ELEM_C_FST, at[n].el_number, ELEM_C_LEN) ) {
                                /* found N(+)(IV)=C(-)(III); prepare conversion to N(V)#C(IV) */
                                num_neigh ++;
                                pos_neigh = i1;
                            }
                        }
                        if ( 1 == num_neigh  &&
                             !has_other_ion_neigh( at, i, n=at[i].neighbor[i1=pos_neigh], en, ne ) &&
                             !has_other_ion_neigh( at, n, i, en, ne )) {
                            /*n = at[i].neighbor[i1=pos_neigh];*/
                            i2 = is_in_the_list( at[n].neighbor, (AT_NUMB)i, at[n].valence ) - at[n].neighbor;
                            at[i].charge --;
                            at[n].charge ++;
                            at[i].chem_bonds_valence ++;
                            at[n].chem_bonds_valence ++;
                            at[i].bond_type[i1] ++;
                            at[n].bond_type[i2] ++;
                            num_changes ++;
                            num_C_minus --;
                            num_N_plus --;
                            num_All -= 2;
                            continue;
                        }
                    }
                    if ( !type || (12 == type && 0 < num_N_plus && 0 < num_N_minus) ) {
                        int num_neigh = 0, pos_neigh = -1;
                        for ( i1 = 0; i1 < at[i].valence; i1 ++ ) {
                            n = at[i].neighbor[i1];
                            if ( -1 == at[n].charge && 2 >= nNoMetalNumBonds(at, n) &&
                                  2 == nNoMetalBondsValence(at, n)+NUMH(at, n) &&
                                  /*0 == num_of_H( at, n ) &&*/
                                  at[i].bond_type[i1] == BOND_TYPE_SINGLE &&
                                  NULL != memchr( en+ELEM_N_FST, at[n].el_number, ELEM_N_LEN) ) {
                                /* found N(+)(IV)=N(-)(II); prepare conversion to N(V)#N(III) */
                                num_neigh ++;
                                pos_neigh = i1;
                            }
                        }
                        if ( 1 == num_neigh  &&
                             !has_other_ion_neigh( at, i, n=at[i].neighbor[i1=pos_neigh], en, ne ) &&
                             !has_other_ion_neigh( at, n, i, en, ne )) {
                            /*n = at[i].neighbor[i1=pos_neigh];*/
                            i2 = is_in_the_list( at[n].neighbor, (AT_NUMB)i, at[n].valence ) - at[n].neighbor;
                            at[i].charge --;
                            at[n].charge ++;
                            at[i].chem_bonds_valence ++;
                            at[n].chem_bonds_valence ++;
                            at[i].bond_type[i1] ++;
                            at[n].bond_type[i2] ++;
                            num_changes ++;
                            num_N_minus --;
                            num_N_plus --;
                            num_All -= 2;
                            continue;
                        }
                    }
                } else
                if ( 1 == at[i].charge &&
                     0 < num_O_plus && 0 < num_C_minus + num_N_minus &&
                     3 >= nNoMetalNumBonds(at, i) && 3 == nNoMetalBondsValence(at, i) &&
                     0 == num_of_H( at, i ) &&
                     NULL != memchr( en+ELEM_O_FST, at[i].el_number, ELEM_O_LEN) ) {
                    /* found non-terminal O(+)(III) */
                    if ( (!type || 13 == type) && 0 < num_C_minus ) {
                        int num_neigh = 0, pos_neigh = -1;
                        for ( i1 = 0; i1 < at[i].valence; i1 ++ ) {
                            n = at[i].neighbor[i1];
                            if ( -1 == at[n].charge && 3 >= at[n].valence && 3 == at[n].chem_bonds_valence+NUMH(at,n) &&
                                  /*0 == at[n].num_H &&*/
                                  at[i].bond_type[i1] == BOND_TYPE_SINGLE &&
                                  NULL != memchr( en+ELEM_C_FST, at[n].el_number, ELEM_C_LEN) ) {
                                /* found O(+)(III)-C(-)(II); prepare conversion to O(IV)=C(IV) */
                                num_neigh ++;
                                pos_neigh = i1;
                            }
                        }
                        if ( 1 == num_neigh  &&
                             !has_other_ion_neigh( at, i, n=at[i].neighbor[i1=pos_neigh], en, ne ) &&
                             !has_other_ion_neigh( at, n, i, en, ne )) {
                            /*n = at[i].neighbor[i1=pos_neigh];*/
                            i2 = is_in_the_list( at[n].neighbor, (AT_NUMB)i, at[n].valence ) - at[n].neighbor;
                            at[i].charge --;
                            at[n].charge ++;
                            at[i].chem_bonds_valence ++;
                            at[n].chem_bonds_valence ++;
                            at[i].bond_type[i1] ++;
                            at[n].bond_type[i2] ++;
                            num_changes ++;
                            num_C_minus --;
                            num_O_plus --;
                            num_All -= 2;
                            continue;
                        }
                    }
                    if ( (!type || 14 == type) && 0 < num_C_minus ) {
                        int num_neigh = 0, pos_neigh = -1;
                        for ( i1 = 0; i1 < at[i].valence; i1 ++ ) {
                            n = at[i].neighbor[i1];
                            if ( -1 == at[n].charge && 3 >= at[n].valence && 3 == at[n].chem_bonds_valence+NUMH(at,n) &&
                                  /*0 == at[n].num_H &&*/
                                  at[i].bond_type[i1] == BOND_TYPE_DOUBLE &&
                                  NULL != memchr( en+ELEM_C_FST, at[n].el_number, ELEM_C_LEN) ) {
                                /* found O(+)(III)=C(-)(III); prepare conversion to O(IV)#C(IV) */
                                num_neigh ++;
                                pos_neigh = i1;
                            }
                        }
                        if ( 1 == num_neigh  &&
                             !has_other_ion_neigh( at, i, n=at[i].neighbor[i1=pos_neigh], en, ne ) &&
                             !has_other_ion_neigh( at, n, i, en, ne )) {
                            /*n = at[i].neighbor[i1=pos_neigh];*/
                            i2 = is_in_the_list( at[n].neighbor, (AT_NUMB)i, at[n].valence ) - at[n].neighbor;
                            at[i].charge --;
                            at[n].charge ++;
                            at[i].chem_bonds_valence ++;
                            at[n].chem_bonds_valence ++;
                            at[i].bond_type[i1] ++;
                            at[n].bond_type[i2] ++;
                            num_changes ++;
                            num_C_minus --;
                            num_O_plus --;
                            num_All -= 2;
                            continue;
                        }
                    }
                    if ( (!type || 15 == type) && 0 < num_N_minus ) {
                        int num_neigh = 0, pos_neigh = -1;
                        for ( i1 = 0; i1 < at[i].valence; i1 ++ ) {
                            n = at[i].neighbor[i1];
                            if ( -1 == at[n].charge && 2 >= nNoMetalNumBonds(at, n) &&
                                  2 == nNoMetalBondsValence(at, n)+NUMH(at, n) &&
                                  /*0 == num_of_H( at, n ) &&*/
                                  at[i].bond_type[i1] == BOND_TYPE_SINGLE &&
                                  NULL != memchr( en+ELEM_N_FST, at[n].el_number, ELEM_N_LEN) ) {
                                /* found O(+)(III)=N(-)(II); prepare conversion to O(IV)#N(III) */
                                num_neigh ++;
                                pos_neigh = i1;
                            }
                        }
                        if ( 1 == num_neigh  &&
                             !has_other_ion_neigh( at, i, n=at[i].neighbor[i1=pos_neigh], en, ne ) &&
                             !has_other_ion_neigh( at, n, i, en, ne )) {
                            /*n = at[i].neighbor[i1=pos_neigh];*/
                            i2 = is_in_the_list( at[n].neighbor, (AT_NUMB)i, at[n].valence ) - at[n].neighbor;
                            at[i].charge --;
                            at[n].charge ++;
                            at[i].chem_bonds_valence ++;
                            at[n].chem_bonds_valence ++;
                            at[i].bond_type[i1] ++;
                            at[n].bond_type[i2] ++;
                            num_changes ++;
                            num_N_minus --;
                            num_O_plus --;
                            num_All -= 2;
                            continue;
                        }
                    }
                }
            }
        }
        /**************************************************************************/
        /*********************** NON-Terminal ion triples *************************/
        /**************************************************************************/
        /*-------------------------------------------------------------------------
           Non-Terminal triple types: 16, 17, 18   N=N,P,As,Sb; O=O,S,Se,Te; C=C,Si
           ========================================
          16: C(+)(III)-O-N(-)(II)  => C(IV)=O=N(III)  (allow terminal H on N(-))

                        |                     |      
          17: C(+)(III)-N-C(-)(III)  => C(IV)=N=C(IV)

          18: C(-)(III)-N=C(+)(III)  => C(IV)=N#C(IV)   (may have two or no charges)
              C(IV)=N-C(II)          => C(IV)=N#C(IV)

        */
        if ( (!type || 16 == type) && 0 < num_C_plus && 0 < num_N_minus ) {
            int m[2], j[2], k;
            for ( i = 0; i < num_atoms; i ++ ) {
                if ( 0 == at[i].charge && 2 == nNoMetalNumBonds(at, i) && 2 == nNoMetalBondsValence(at, i) &&
                     0 == num_of_H( at, i ) &&
                     0 <= (j[0] = nNoMetalNeighIndex( at, i )) &&
                     at[m[0]=at[i].neighbor[j[0]]].charge &&
                     0 <= (j[1] = nNoMetalOtherNeighIndex( at, i, m[0] )) &&
                     0 == at[m[0]].charge + at[m[1]=at[i].neighbor[j[1]]].charge &&
                     5 >= nNoMetalBondsValence(at, m[0]) + nNoMetalBondsValence(at, m[1]) &&
                     /*5 >= at[m[0]].chem_bonds_valence + at[m[1]].chem_bonds_valence &&*/
                     NULL != memchr( en+ELEM_O_FST, at[i].el_number, ELEM_O_LEN) ) {
                    /* found non-terminal A(+)-O-B(-); chem_bond_val of A+B <= 5 */
                    int n_N=-1, n_C=-1, i_C=-1;
                    for ( k = 0; k < 2; k ++ ) {
                        n = m[k];
                        if ( -1 == at[n].charge && 2 == nNoMetalNumBonds(at, n)+NUMH(at, n) &&
                             /*0 == num_of_H( at, n ) &&*/
                             NULL != memchr( en+ELEM_N_FST, at[n].el_number, ELEM_N_LEN) ) {
                            n_N = n;
                        } else
                        if ( 1 == at[n].charge && 3 == at[n].chem_bonds_valence+NUMH(at,n) &&
                             NULL != memchr( en+ELEM_C_FST, at[n].el_number, ELEM_C_LEN) ) {
                            n_C = n;
                            i_C = k;
                        }
                    }
                    if ( n_C < 0 || n_N < 0 ||
                         has_other_ion_in_sphere_2(at, n_C, n_N, en, ne ) ||
                         has_other_ion_in_sphere_2(at, n_N, n_C, en, ne ) ) {
                        continue;
                    }
                    /* C(+)(III)-O-N(-)(II)  => C(IV)=O=N(III) */
                    for ( k = 0; k < 2; k ++ ) {
                        n  = k? n_C : n_N;
                        i1 = k? j[i_C] : j[1-i_C];
                        i2 = is_in_the_list( at[n].neighbor, (AT_NUMB)i, at[n].valence ) - at[n].neighbor;
                        at[i].bond_type[i1] ++;
                        at[n].bond_type[i2] ++;
                        at[i].chem_bonds_valence ++;
                        at[n].chem_bonds_valence ++;
                        at[n].charge += (k? -1:1);
                    }
                    num_changes ++;
                    num_N_minus --;
                    num_C_plus --;
                    num_All -= 2;
                }
            }
        }
        if ( (!type || 17 == type) && 0 < num_C_plus && 0 < num_C_minus ) {
            int m[3], c[3], j[3], k;
            for ( i = 0; i < num_atoms; i ++ ) {
                if ( 0 == at[i].charge && 3 == nNoMetalNumBonds(at, i) && 3 == nNoMetalBondsValence(at, i) &&
                     0 == num_of_H( at, i ) &&
                     0 <= ( j[0] = nNoMetalNeighIndex(at, i) ) &&
                     0 <= ( j[1] = nNoMetalOtherNeighIndex( at, i, m[0] = at[i].neighbor[j[0]] ) ) &&
                     0 <= ( j[2] = nNoMetalOtherNeighIndex2( at, i, m[0], m[1] = at[i].neighbor[j[1]] ) ) &&
                     1 == !(c[0]=at[m[0]].charge) 
                        + !(c[1]=at[m[1]].charge)
                        + !(c[2]=at[m[2]=at[i].neighbor[j[2]]].charge) &&
                     0 == c[0] + c[1] + c[2] &&
                     2 == (3== (c[0]? at[m[0]].chem_bonds_valence+NUMH(at,m[0]):0))
                        + (3== (c[1]? at[m[1]].chem_bonds_valence+NUMH(at,m[1]):0))
                        + (3== (c[2]? at[m[2]].chem_bonds_valence+NUMH(at,m[2]):0)) &&
                     NULL != memchr( en+ELEM_N_FST, at[i].el_number, ELEM_N_LEN) ) {
                    /* found non-terminal A(+)-O-B(-) */
                    int n_Cp=-1, n_Cm=-1, i_Cp=-1, i_Cm=-1; /* p = positive, m = negatice ion C */
                    for ( k = 0; k < 3; k ++ ) {
                        if ( c[k] ) {
                            n = m[k];
                            if ( -1 == at[n].charge &&
                                 NULL != memchr( en+ELEM_C_FST, at[n].el_number, ELEM_C_LEN) ) {
                                n_Cm = n;
                                i_Cm = k;
                            } else
                            if ( 1 == at[n].charge &&
                                 NULL != memchr( en+ELEM_C_FST, at[n].el_number, ELEM_C_LEN) ) {
                                n_Cp = n;
                                i_Cp = k;
                            }
                        }
                    }
                    if ( n_Cp < 0 || n_Cm < 0 ||
                         has_other_ion_in_sphere_2(at, n_Cp, n_Cm, en, ne ) ||
                         has_other_ion_in_sphere_2(at, n_Cm, n_Cp, en, ne )) {
                        continue;
                    }
                    /*           |                     |       */
                    /* C(+)(III)-N-C(-)(III)  => C(IV)=N=C(IV) */
                    for ( k = 0; k < 2; k ++ ) {
                        n  = k? n_Cp : n_Cm;
                        i1 = k? j[i_Cp] : j[i_Cm];
                        i2 = is_in_the_list( at[n].neighbor, (AT_NUMB)i, at[n].valence ) - at[n].neighbor;
                        at[i].bond_type[i1] ++;
                        at[n].bond_type[i2] ++;
                        at[i].chem_bonds_valence ++;
                        at[n].chem_bonds_valence ++;
                        at[n].charge += (k? -1:1);
                    }
                    num_changes ++;
                    num_C_minus --;
                    num_C_plus --;
                    num_All -= 2;
                }
            }
        }
        if ( (!type || 18 == type) && ((0 < num_C_plus && 0 < num_C_minus) || 0 < num_C_II) ) {
            int m[2], v[2], j[2], k;
            for ( i = 0; i < num_atoms; i ++ ) {
                if ( 0 == at[i].charge && 2 == nNoMetalNumBonds(at, i) && 3 == nNoMetalBondsValence(at, i) &&
                     0 == num_of_H( at, i ) &&
                     0 <= (j[0] = nNoMetalNeighIndex( at, i )) &&
                     0 <= (j[1] = nNoMetalOtherNeighIndex( at, i, m[0] = at[i].neighbor[j[0]] )) &&
                     0 == at[m[0]].charge
                         +at[m[1]=at[i].neighbor[j[1]]].charge &&
                     6 == (v[0]=at[m[0]].chem_bonds_valence+NUMH(at,m[0]))
                         +(v[1]=at[m[1]].chem_bonds_valence+NUMH(at,m[1])) &&
                     2 >= abs(v[0]-v[1]) &&
                     NULL != memchr( en+ELEM_N_FST, at[i].el_number, ELEM_N_LEN) &&
                     NULL != memchr( en+ELEM_C_FST, at[m[0]].el_number, ELEM_C_LEN) &&
                     NULL != memchr( en+ELEM_C_FST, at[m[1]].el_number, ELEM_C_LEN)
                    ) {
                    /*                    n_Cm      i n_Cp */
                    /* found non-terminal C(-)(III)-N=C(+)(III) or C(IV)=N-C(II): Cm-N-Cp */
                    /* convert to C(IV)=N#C(IV) */
                    int n_Cp=-1, n_Cm=-1, i_Cp=-1, i_Cm=-1; /* p = positive, m = negatice ion C */
                    for ( k = 0; k < 2; k ++ ) {
                        n = m[k];
                        if ( v[k] == 4 || (v[k] == 3 && at[i].bond_type[j[k]] == BOND_TYPE_SINGLE) ) {
                            n_Cm = n;
                            i_Cm = k;
                        } else
                        if ( v[k] == 2 || (v[k] == 3 && at[i].bond_type[j[k]] == BOND_TYPE_DOUBLE) ) {
                            n_Cp = n;
                            i_Cp = k;
                        }
                    }
                    if ( n_Cp < 0 || n_Cm < 0 || at[n_Cp].valence+NUMH(at,n_Cp) != 2 ) {
                        continue; /* guarantees at[n_Cp].valence <= 2 */
                    }
                    if ( v[i_Cp] == 2 || !at[n_Cp].charge ) {
                        if ( at[n_Cp].valence == 2 ) {
                            /* neighbor of at[n_Cp] opposite to at[i] */
                            k = at[n_Cp].neighbor[at[n_Cp].neighbor[0]==i];
                            if ( NULL != memchr( en+ELEM_N_FST, at[k].el_number, ELEM_N_LEN) ) {
                                continue;
                            }
                        }
                    } else
                    if ( at[n_Cp].charge ) {
                        if ( has_other_ion_in_sphere_2(at, n_Cp, n_Cm, en, ne ) ||
                             has_other_ion_in_sphere_2(at, n_Cm, n_Cp, en, ne )) {
                            continue;
                        }
                    } else {
                        continue; /* unknown case */
                    }
                    /*                                         */
                    /* C(-)(III)-N=C(+)(III)  => C(IV)=N#C(IV) */
                    /* C(IV)=N-C(II)          => C(IV)=N#C(IV) */
                    if ( at[n_Cp].charge ) {
                        num_C_minus --;
                        num_C_plus --;
                        num_All -= 2;
                    } else {
                        num_C_II --;
                        num_All --;
                    }
                    
                    for ( k = 0; k < 2; k ++ ) {
                        n  = k? n_Cp : n_Cm;
                        i3 = k? i_Cp : i_Cm; /* added to fix the bug */
                        /*i1 = k? j[i_Cp] : j[i_Cm];*/ /* replaced with next line */
                        i1 = j[i3];
                        if ( v[i3 /*was i1*/] < 4 ) { /* WDI found a bug here: bounds violation */
                            int delta = 4 - v[i3 /*was i1*/];
                            i2 = is_in_the_list( at[n].neighbor, (AT_NUMB)i, at[n].valence ) - at[n].neighbor;
                            at[i].bond_type[i1] += delta;
                            at[n].bond_type[i2] += delta;
                            at[i].chem_bonds_valence += delta;
                            at[n].chem_bonds_valence += delta;
                            at[n].charge = 0;
                            at[n].radical = 0;
                        }
                    }
                    at[i].charge = 0;
                    at[i].radical = 0;
                    num_changes ++;
                }
            }
        }
    }

    return num_changes;
}

/*#if ( DISCONNECT_SALTS == 1 )*/ /* { */



/*************************************************************************************************/
int RemoveInpAtBond( inp_ATOM *atom, int iat, int k )
{
    int      i, j, m, m2, k2;
    inp_ATOM *at = atom + iat;
    inp_ATOM *at2 = NULL;
    int      val = at->valence - 1;
    if ( val >= 0 ) {
        int bond = at->bond_type[k];
        if ( bond > BOND_TYPE_TRIPLE )
            bond = BOND_TYPE_SINGLE; /* added 08-06-2003 */
        
        /* update CML tetrahedral atom parity. */
        if ( at->p_parity ) {
            for( m = 0; m < MAX_NUM_STEREO_ATOM_NEIGH; m ++ ) {
                if ( at->p_orig_at_num[m] == at->orig_at_number ) {
                    at->p_parity = 0;
                    break; /* only 3 bonds are present; removing one bond removes stereo */
                }
            }
            if ( at->p_parity /* at->valence == MAX_NUM_STEREO_ATOM_NEIGH*/ ) {
                for ( m = 0; m < at->valence; m ++ ) {
                    if ( atom[(int)at->neighbor[k]].orig_at_number == at->p_orig_at_num[m] ) {
                        break;
                    }
                }
                if ( m < at->valence ) {
                    at->p_orig_at_num[m] = at->orig_at_number;
                } else {
                    at->p_parity = 0; /* wrong neighbors: at->neighbor[k] is not in the list of a stereo neighbors */
                }
            }
        }
        
        /* update CML stereogenic bond parities; at this point no removed explicit H exist yet */
        if ( at->sb_parity[0] ) {
            for ( m = 0; m < MAX_NUM_STEREO_BONDS && at->sb_parity[m]; ) {
                if ( k == at->sb_ord[m] || (k == at->sn_ord[m] && val < 2 && ATOM_PARITY_WELL_DEF(at->sb_parity[m])) ) {
                    /* !!! FLAW: does take into account removed H !!! */
                    /* stereogenic bond is being removed OR */
                    /* remove stereogenic bond because its only neighbor is being removed */
                    int pnxt_atom, pinxt2cur, pinxt_sb_parity_ord;
                    int len= get_opposite_sb_atom( atom, iat, at->sb_ord[m], &pnxt_atom, &pinxt2cur, &pinxt_sb_parity_ord );
                    if ( len ) {
                        i   = pinxt_sb_parity_ord;
                        at2 = atom + pnxt_atom;
                        k2  = pinxt2cur;
                    } else {
                        i = MAX_NUM_STEREO_BONDS;
                    }
                    /*
                    at2 = atom + at->neighbor[ (int)at->sb_ord[m] ];
                    for ( i = 0; i < MAX_NUM_STEREO_BONDS && at2->sb_parity[i]; i ++ ) {
                        if ( iat == at2->neighbor[ (int)at2->sb_ord[i] ] )
                            break;
                    }
                    */
                    if ( i < MAX_NUM_STEREO_BONDS && at2->sb_parity[i] ) {
                        m2 = i;
                        /* remove bond parity from at */
                        if ( m < MAX_NUM_STEREO_BONDS-1 ) {
                            memmove( at->sb_parity+m, at->sb_parity+m+1, (MAX_NUM_STEREO_BONDS-1 - m) * sizeof(at->sb_parity[0]));
                            memmove( at->sb_ord+m, at->sb_ord+m+1, (MAX_NUM_STEREO_BONDS-1 - m) * sizeof(at->sb_ord[0]));
                            memmove( at->sn_ord+m, at->sn_ord+m+1, (MAX_NUM_STEREO_BONDS-1 - m) * sizeof(at->sn_ord[0]));
                            memmove( at->sn_orig_at_num+m, at->sn_orig_at_num+m+1, (MAX_NUM_STEREO_BONDS-1 - m) * sizeof(at->sn_orig_at_num[0]));
                        }
                        at->sb_parity[MAX_NUM_STEREO_BONDS-1] = 0;
                        at->sb_ord[MAX_NUM_STEREO_BONDS-1] = 0;
                        at->sn_ord[MAX_NUM_STEREO_BONDS-1] = 0;
                        at->sn_orig_at_num[MAX_NUM_STEREO_BONDS-1] = 0;
                        /* remove bond parity from at2 */
                        if ( m2 < MAX_NUM_STEREO_BONDS-1 ) {
                            memmove( at2->sb_parity+m2, at2->sb_parity+m2+1, (MAX_NUM_STEREO_BONDS-1 - m2) * sizeof(at2->sb_parity[0]));
                            memmove( at2->sb_ord+m2, at2->sb_ord+m2+1, (MAX_NUM_STEREO_BONDS-1 - m2) * sizeof(at2->sb_ord[0]));
                            memmove( at2->sn_ord+m2, at2->sn_ord+m2+1, (MAX_NUM_STEREO_BONDS-1 - m2) * sizeof(at2->sn_ord[0]));
                            memmove( at2->sn_orig_at_num+m2, at2->sn_orig_at_num+m2+1, (MAX_NUM_STEREO_BONDS-1 - m2) * sizeof(at2->sn_orig_at_num[0]));
                        }
                        at2->sb_parity[MAX_NUM_STEREO_BONDS-1] = 0;
                        at2->sb_ord[MAX_NUM_STEREO_BONDS-1] = 0;
                        at2->sn_ord[MAX_NUM_STEREO_BONDS-1] = 0;
                        at2->sn_orig_at_num[MAX_NUM_STEREO_BONDS-1] = 0;
                        /* do not increment m here because the array elements have been shifted */
                    } else {
                        m ++; /* program error: inconsistent stereobond parity */
                    }
                } else
                if ( k == at->sn_ord[m] ) {
                    /* stereogenic bond neighbor is being removed; another neighbor remains */
                    /* !!! FLAW: does take into account removed H !!! */
                    for ( j = 0, i = -1; j < at->valence; j ++ ) {
                        if ( j != k && j != at->sb_ord[m] ) {
                            i = j;
                            break;
                        }
                    }
                    /* i is the position of the neighbor that will become a new neighbor */
                    /***************************************************************************
                     *  at->sb_parity[m] is the direction (EVEN=clockwise, ODD=counterclockwise)
                     *  from stereobond to the neighbor. If the neighbor is removed then
                     *  the parity should invert, otherwise it should be unchanged.
                     ***************************************************************************/
                    if ( i < 0 ) {
                        /* no alternative neighbor is available */
                        if ( ATOM_PARITY_WELL_DEF(at->sb_parity[m] ) ) {
                            /* parity cannot be not well-defined anymore */
                            int pnxt_atom, pinxt2cur, pinxt_sb_parity_ord;
                            int len= get_opposite_sb_atom( atom, iat, at->sb_ord[m], &pnxt_atom, &pinxt2cur, &pinxt_sb_parity_ord );
                            if ( len > 0 ) {
                                atom[pnxt_atom].sb_parity[pinxt_sb_parity_ord] = at->sb_parity[m] = AB_PARITY_UNDF;
                            }
#ifdef _DEBUG
                            else {
                                int stop = 1; /* sb parities error */
                            }
#endif
                        }
                        at->sn_ord[m]         = -99; /* sb neighbor has been disconnected */
                        at->sb_ord[m]        -= (at->sb_ord[m] > k); /* same as above */
                        at->sn_orig_at_num[m] = 0;
                    } else
                    if ( i < at->valence ) {
                        /* choose another stereogenic bond neighbor, its ord. number is i before bond removal */
                        if ( ATOM_PARITY_WELL_DEF(at->sb_parity[m]) ) {
                            /* ALL WRONG: 'move' previous stereo bond neighbor to the last position (pos. 2 out of 0,1,2) */
                            /* the parity of the transpositions is (2 - at->sn_ord[m])%2 = at->sn_ord[m] % 2 */
                            /* and replace the neighbor with another; the contribution to the parity is 1 */
                            
                            /*at->sb_parity[m]      =  2 - ( at->sb_parity[m] + at->sn_ord[m] + 1 ) % 2;*/
                            
                            /*at->sb_parity[m]      =  2 - ( at->sb_parity[m] + k + i +
                                                           (i > k) + (i > at->sb_ord[m]) ) % 2;*/
                            /*=== parity should be INVERTED ===*/
                            at->sb_parity[m] = 3 - at->sb_parity[m];
                        }
                        at->sn_ord[m]         = i - (i > k); /* ord. number shifted because preceding bond is removed */
                        at->sb_ord[m]        -= (at->sb_ord[m] > k); /* same as above */
                        at->sn_orig_at_num[m] = atom[(int)at->neighbor[i]].orig_at_number;
                        /*at->sb_parity[m]      =  2 - ( at->sb_parity[m] + 1 ) % 2;*/
                    } else {
                        at->sb_parity[m] = 0; /* program error: inconsistent stereobond parity */
                    }
                    m ++;
                } else {
                    /* removing another neighbor, k: first move it to the last position (pos. 2 out of 0,1,2) */
                    if ( k < 2 && ATOM_PARITY_WELL_DEF(at->sb_parity[m]) ) {
                        /*at->sb_parity[m] =  2 - ( at->sb_parity[m] + k ) % 2;*/
                        /*at->sb_parity[m] =  2 - ( at->sb_parity[m] + (at->sn_ord[m] > k) + (at->sb_ord[m] > k) ) % 2;*/
                        ;/*==== Parity should remain UNCHANGED ===*/
                    }
                    if ( at->sb_ord[m] > k ) {
                        at->sb_ord[m] --;
                    }
                    if ( at->sn_ord[m] > k ) {
                        at->sn_ord[m] --;
                    }
                    m ++;
                }
            }
        }

        if ( k < val ) {
            memmove( at->neighbor+k, at->neighbor+k+1, sizeof(at->neighbor[0])*(val-k) );
            memmove( at->bond_stereo+k, at->bond_stereo+k+1, sizeof(at->bond_stereo[0])*(val-k) );
            memmove( at->bond_type+k, at->bond_type+k+1, sizeof(at->bond_type[0])*(val-k) );
        }
        at->neighbor[val]    = 0;
        at->bond_stereo[val] = 0;
        at->bond_type[val]   = 0;
        at->valence = val;
        at->chem_bonds_valence -= bond;
        return 1;
    }
    return 0;
}



/*************************************************************************************************/
int DisconnectInpAtBond( inp_ATOM *at, AT_NUMB *nOldCompNumber, int iat, int neigh_ord )
{
    int neigh, i, ret = 0;
    int component;
    neigh = at[iat].neighbor[neigh_ord];
    for ( i = 0; i < at[neigh].valence; i ++ ) {
        if ( iat == (int)at[neigh].neighbor[i] )
            break;
    }
    if ( i < at[neigh].valence ) {
        ret += RemoveInpAtBond( at, iat, neigh_ord );
        ret += RemoveInpAtBond( at, neigh, i );
        if ( nOldCompNumber && ret ) {
            if ( (component = at[iat].component) ) {
                nOldCompNumber[component-1] = 0;
            }
            if ( (component = at[neigh].component) ) {
                nOldCompNumber[component-1] = 0;
            }
        }
    }
    return (ret == 2);
}



/*************************************************************************************************/
int bIsAmmoniumSalt( inp_ATOM *at, int i, int *piO, int *pk, S_CHAR *num_explicit_H ) 
{
    /* NH4(+charge)-O(-charge)-C -> NH3 + HO-C; any charge including 0, any C except charged or radical */
    /* F, Cl, Br, I */
    static U_CHAR el_number_C=0, el_number_O=0, el_number_H=0, el_number_N=0;
    static U_CHAR el_number_F=0, el_number_Cl=0, el_number_Br=0, el_number_I=0;
    int num_H, num_non_iso_H, num_impl_iso_H, bDisconnect = 1;
    int j, val, neigh, iO=-1, iC, k=-1;
    if ( 0 == el_number_C ) {
        /* one time initialization */
        el_number_C = get_periodic_table_number( "C" );
        el_number_O = get_periodic_table_number( "O" );
        el_number_H = get_periodic_table_number( "H" );
        el_number_N = get_periodic_table_number( "N" );
        el_number_F = get_periodic_table_number( "F" );
        el_number_Cl= get_periodic_table_number( "Cl" );
        el_number_Br= get_periodic_table_number( "Br" );
        el_number_I = get_periodic_table_number( "I" );
    }
    if ( at[i].el_number != el_number_N )
        return 0;

    /* check for NH4-O-C... -> NH3 + HO-C... */
    val            = at[i].valence;
    num_impl_iso_H = NUM_ISO_H(at,i);
    num_non_iso_H  = at[i].num_H;
    num_H = num_non_iso_H + num_impl_iso_H;
    if ( val + num_H == 5 ) {
        int num_O = 0;
        memset( num_explicit_H, 0, (NUM_H_ISOTOPES+1)*sizeof(num_explicit_H[0]) );
        for ( j = 0; j < val; j ++ ) { /* looking for O: H4N-O-C... */
            neigh = at[i].neighbor[j];
            if ( at[neigh].num_H ||
                 (at[neigh].charge && (at[neigh].el_number != el_number_O || at[neigh].charge + at[i].charge)) ||
                 (at[neigh].radical && at[neigh].radical != RADICAL_SINGLET) ) {
                bDisconnect = 0;
                break; /* reject */
            }
            if ( at[neigh].el_number == el_number_H && at[neigh].valence == 1 &&
                 !at[neigh].charge && !at[neigh].radical ) {
                num_H ++; /* at this point at[].num_H does not include explicit H count */
                num_non_iso_H += (0==at[neigh].iso_atw_diff);
                num_explicit_H[at[neigh].iso_atw_diff] ++;  /* explicit H on N */
            } else
            if ( at[neigh].el_number == el_number_O && at[neigh].valence == 2 && !num_O ) {
                num_O ++; /* found O: N-O- */
                iO = neigh;
                k  = j;
                iC = at[iO].neighbor[at[iO].neighbor[0] == i];
                if ( at[iC].el_number != el_number_C || /*
                     at[iC].num_H ||
                     at[iC].chem_bonds_valence != 4 || */
                     at[iC].charge         ||
                     (at[iC].radical && at[iC].radical != RADICAL_SINGLET) /*||
                     at[iC].valence == at[iC].chem_bonds_valence*/ ) {
                    bDisconnect = 0;
                    break; /* reject */
                }
            } else
            if ( (at[neigh].el_number == el_number_F  ||
                  at[neigh].el_number == el_number_Cl ||
                  at[neigh].el_number == el_number_Br ||
                  at[neigh].el_number == el_number_I ) &&
                  at[neigh].valence == 1 && at[neigh].chem_bonds_valence == 1 &&
                  !at[neigh].charge && !NUMH(at,neigh) && !num_O ) {
                num_O ++; /* found O: N-O- */
                iO = neigh;
                k  = j;
                iC = -1;
            } else {
                bDisconnect = 0;
                break;  /* reject */
            }
        }
        if ( bDisconnect && (num_O != 1 || num_H != 4) ) {
            bDisconnect = 0; /* reject */
        }
    } else {
        bDisconnect = 0;
    }
    if ( bDisconnect ) {
        *piO = iO;
        *pk  = k;
    }
    return bDisconnect;
}



/*************************************************************************************************/
int DisconnectAmmoniumSalt ( inp_ATOM *at, int iN, int iO, int k, S_CHAR *num_explicit_H )
{
    /* disconnect NH4-O from O */
    /* Note: iO = at[iN].neighbor[k], at[iN] is N, at[iO].neighbor[0] is either N=at[iN] or C=at[iC] */
    int nMove_H_iso_diff = -1; /* do not move explicit H */
    int j, neigh, iso_diff, neigh_pos;
    static U_CHAR el_number_H = 0;
    int    val = at[iN].valence;

    if ( !el_number_H ) {
        el_number_H = get_periodic_table_number( "H" );
    }
    if ( at[iN].charge && !(at[iN].charge + at[iO].charge) ) {
        at[iN].charge = at[iO].charge = 0; /* remove charges */
    }
    neigh_pos = (at[iO].valence == 2)? (at[iO].neighbor[1] == iN) : 0; /* position of at[iN] in the neigh list of iO */
    /* disconnect bond O-N */
    RemoveInpAtBond( at, iO, neigh_pos );
    RemoveInpAtBond( at, iN, k );
    val --;

    /* move 1 H from NH4 to O- or Cl */

    /* find non-isotopic or the lightest isotopic H to move from N to O */
    for ( iso_diff = 0; iso_diff <= NUM_H_ISOTOPES; iso_diff ++ ) {
        if ( !iso_diff ) {
            /* find non-isotopic H */
            if ( at[iN].num_H ) {
                at[iN].num_H --;  /* move non-isotopic implicit H */
                at[iO].num_H ++;
                break;
            } else
            if ( num_explicit_H[0] ) {
                nMove_H_iso_diff = 0; /* flag: move explicit non-isotopic H */
                break;
            }
        } else {
            /* find isotopic H */
            if ( at[iN].num_iso_H[iso_diff] ) {
                at[iN].num_iso_H[iso_diff] --; /* move implicit isotopic H, atw = 1 */
                at[iO].num_iso_H[iso_diff] ++;
                break;
            } else
            if ( num_explicit_H[iso_diff] ) {
                nMove_H_iso_diff = iso_diff; /* flag: move explicit isotopic H, atw = 1 */
                break;
            }
        }
    }
    if ( nMove_H_iso_diff >= 0 ) {
        /* move explicit H, it is isotopic if nMove_H_iso_diff > 0 */
        double dist2_H_O, min_dist2_H_O = -1.0;
        int    jH = -1, iH = -1;
        for ( j = 0; j < val; j ++ ) { /* looking H in N-H such that H-O is shortest */
            neigh = at[iN].neighbor[j];
            if ( at[neigh].el_number    == el_number_H &&
                at[neigh].iso_atw_diff == nMove_H_iso_diff ) {
                dist2_H_O =  (at[neigh].x - at[iO].x) * (at[neigh].x - at[iO].x) +
                             (at[neigh].y - at[iO].y) * (at[neigh].y - at[iO].y) +
                             (at[neigh].z - at[iO].z) * (at[neigh].z - at[iO].z);
                if ( min_dist2_H_O < 0.0 || min_dist2_H_O > dist2_H_O ) {
                    min_dist2_H_O = dist2_H_O;
                    iH = neigh;
                    jH = j;
                }
            }
        }
        /* reconnect; bonds do not need changes except stereo */
        neigh_pos = at[iO].valence;
        at[iO].neighbor[neigh_pos]         = iH;
        at[iO].bond_stereo[neigh_pos]      = 0;
        at[iO].bond_type[neigh_pos]        = at[iH].bond_type[0];
        at[iO].chem_bonds_valence         += at[iH].bond_type[0];
        at[iO].valence ++;
        at[iH].neighbor[0]                 = iO;
        at[iH].bond_stereo[0]              = 0;
        /* disconnect H from N */
        RemoveInpAtBond( at, iN, jH );
        val --;
        if ( k > jH ) {
            k --;
        }
    }
    return 1;
}



/*************************************************************************************************/
int bIsMetalSalt( inp_ATOM *at, int i ) 
{
    int type, val, k, iO, iC, j, neigh;
    int bDisconnect = 1;
    static U_CHAR el_number_C=0, el_number_O=0, el_number_H=0;
    static U_CHAR el_number_F=0, el_number_Cl=0, el_number_Br=0, el_number_I=0;
    if ( 0 == el_number_C ) {
        /* one time initialization */
        el_number_C = get_periodic_table_number( "C" );
        el_number_O = get_periodic_table_number( "O" );
        el_number_H = get_periodic_table_number( "H" );
        el_number_F = get_periodic_table_number( "F" );
        el_number_Cl= get_periodic_table_number( "Cl" );
        el_number_Br= get_periodic_table_number( "Br" );
        el_number_I = get_periodic_table_number( "I" );
    }
    /* check for a metal atom:
       metal atom should be connected and be a metal */
    if ( !(val = at[i].valence) ||
         !(type = get_el_type( at[i].el_number )) ||
         !(type & IS_METAL) ) {
        bDisconnect = 0;  /* reject */
    } else
    /* metal atom should not have adjacent H or multiple bonds or radical */
    if ( at[i].num_H ) {
        bDisconnect = 0; /* reject */
    } else
    /* check valence */
    if ( (at[i].charge == 0 &&
         ( ((type & 1) && val == get_el_valence( at[i].el_number, 0, 0 ))   ||
           ((type & 2) && val == get_el_valence( at[i].el_number, 0, 1 )) )) ||
         (at[i].charge > 0 &&
         (type & 1) && val == get_el_valence( at[i].el_number, at[i].charge, 0 )) ) {
        ; /* accept */
    } else {
        bDisconnect = 0; /* reject */
    }
    if ( bDisconnect ) {
        /*************************************************************************
         *                                                                  |    *
         * check M neighbors. Disconnect if all neighbors are M-O-C# or M-O-C=   *
         *                                                                  |    *
         *************************************************************************/
        for ( k = 0; k < at[i].valence; k ++ ) {
            iO = at[i].neighbor[k];
            /* halogenide 2004-07-08 */
            if ( (at[iO].el_number == el_number_F  ||
                  at[iO].el_number == el_number_Cl ||
                  at[iO].el_number == el_number_Br ||
                  at[iO].el_number == el_number_I ) &&
                  at[iO].valence == 1 && at[iO].chem_bonds_valence == 1 &&
                  !at[iO].charge && !(at[iO].radical && at[iO].radical != RADICAL_SINGLET) && !NUMH(at,iO) ) {
                    ; /* found */
                  } else {
                /* -O-C= */
                if ( at[iO].el_number != el_number_O ||
                    NUMH(at, iO) ||
                    at[iO].valence   != 2 ||
                    at[iO].charge         ||
                    (at[iO].radical && at[iO].radical != RADICAL_SINGLET) ||
                    at[iO].valence != at[iO].chem_bonds_valence ) {
                    bDisconnect = 0; /* reject */
                    break;
                }
                iC = at[iO].neighbor[at[iO].neighbor[0] == i];
                if ( at[iC].el_number != el_number_C ||
                    at[iC].num_H ||
                    at[iC].chem_bonds_valence != 4 ||
                    at[iC].charge         ||
                    (at[iC].radical && at[iC].radical != RADICAL_SINGLET) ||
                    at[iC].valence == at[iC].chem_bonds_valence ) {
                    bDisconnect = 0; /* reject */
                    break;
                }
                for ( j = 0; j < at[iC].valence; j ++ ) {
                    neigh = at[iC].neighbor[j];
                    if ( at[neigh].el_number == el_number_H ) {
                        break;
                    }
                }
                if ( j != at[iC].valence ) {
                    bDisconnect = 0; /* reject */
                    break;
                }
            }
        }
    }
    return bDisconnect;
}



/*************************************************************************************************/
int DisconnectMetalSalt( inp_ATOM *at, int i )
{
    int k, iO;
    /* disconnect metal atom or ion at[i] */
    for ( k = 0; k < at[i].valence; k ++ ) {
        iO = at[i].neighbor[k];
        if ( at[iO].valence == 2 ) {
            if ( at[iO].neighbor[0] == i ) { /* assuming atom O always has 2 bonds */
                /* copy the remaining neighbor to the 0 position */
                at[iO].neighbor[0]    = at[iO].neighbor[1];
                at[iO].bond_stereo[0] = at[iO].bond_stereo[1];
                at[iO].bond_type[0]   = at[iO].bond_type[1];
            }
            /* clear neighbor at position 1 */
            at[iO].neighbor[1]    = 0;
            at[iO].bond_stereo[1] = 0;
            at[iO].bond_type[1]   = 0;
        } else {
            /* clear neighbor at position 1 */
            at[iO].neighbor[0]    = 0;
            at[iO].bond_stereo[0] = 0;
            at[iO].bond_type[0]   = 0;
        }
        /* make O negatively charged */
        at[iO].charge = -1;
        /* reduce O valence to account for the removed single bond */
        at[iO].valence --;
        at[iO].chem_bonds_valence --;

        /* clear metal neighbor (O) */
        at[i].neighbor[k]    = 0;
        at[i].bond_stereo[k] = 0;
        at[i].bond_type[k]   = 0;
        /* add a positive charge to the metal */
        at[i].charge ++;
    }
    /* set metal valence to zero because it has been disconnected */
    at[i].valence            = 0;
    at[i].chem_bonds_valence = 0;
    return k;
}



/*************************************************************************************************/
int DisconnectSalts( ORIG_ATOM_DATA *orig_inp_data, int bDisconnect )
{
    int i, k, iO, num_changes, val;
    S_CHAR    num_explicit_H[NUM_H_ISOTOPES+1];
    inp_ATOM *at = orig_inp_data->at;
    int num_at   = orig_inp_data->num_inp_atoms;

    /* check each atom */
    for ( i = 0, num_changes = 0; i < num_at; i ++ ) {

        if ( !(val = at[i].valence) || /* disconnected atom */
             val != at[i].chem_bonds_valence || /* a bond has higher multiplicity than 1 */
             (at[i].radical && at[i].radical != RADICAL_SINGLET) /* radical */ ) {
            continue;   /* reject */
        }
        if ( bIsAmmoniumSalt( at, i, &iO, &k, num_explicit_H ) ) {
            if ( bDisconnect ) {
                DisconnectAmmoniumSalt ( at, i, iO, k, num_explicit_H );
                orig_inp_data->num_inp_bonds --;
            }
            /* count disconnected atoms */
            num_changes ++;
        } else
        if ( bIsMetalSalt( at, i ) ) {
            if ( bDisconnect ) {
                k = DisconnectMetalSalt( at, i );
                orig_inp_data->num_inp_bonds -= k;
            }
            num_changes ++;
        }
    }   
    return num_changes;
}



/*****************************************************************************/
/* Important: Salt disconnection is independent from coord. disconnection:   */
/* because different atoms are disconnected.                                 */
/* However, sal disconnection may need to be rerun after metal disconnection */
/* because metal disconnection may make certain atoms be eligible for salt   */
/* disconnection                                                             */
/*****************************************************************************/
int bIsMetalToDisconnect(inp_ATOM *at, int i, int bCheckMetalValence)
{
    int type, at_valence, num_H;
/*
    if ( !at[i].valence )
*/
    if ( !(type = get_el_type( at[i].el_number )) ||
         !(type & IS_METAL ) ) {
        return 0;
    }
    num_H      = NUMH(at,i);
    at_valence = num_H + at[i].chem_bonds_valence;
    if ( !at_valence ) {
        return 0; /* nothing to disconnect */
    }
    if ( bCheckMetalValence ) {
        if ( abs(at[i].charge) > 1 ) {
            return 1; /* multiple charges */
        }
        for ( i = 0; i < 2 && (i & type); i ++ ) {
            if ( at_valence ==  get_el_valence( at[i].el_number, at[i].charge, i ) ) {
                return 2; /* atom has normal valence */
            }
        }
    }
    return 1;

}



/*****************************************************************************/
int bMayDisconnectMetals( ORIG_ATOM_DATA *orig_inp_data, int bCheckMetalValence, INCHI_MODE *bTautFlagsDone )
{
    int i, j, k, iO, num_changes, val, bRadOrMultBonds, num_impl_H = 0;
    S_CHAR    num_explicit_H[NUM_H_ISOTOPES+1];
    inp_ATOM *at   =  orig_inp_data->at;
    int num_at     =  orig_inp_data->num_inp_atoms;
    int *nNumImplH = &orig_inp_data->bDisconnectCoord;
    /* check each atom */
    for ( i = 0, num_changes = 0; i < num_at; i ++ ) {

        if ( !(val = at[i].valence) && !NUMH(at,i) ) {
            continue; /* disconnected atom */
        }
        bRadOrMultBonds = (val == 0) ||
             (val != at[i].chem_bonds_valence) || /* a bond has higher multiplicity than 1 */
             (at[i].radical && at[i].radical != RADICAL_SINGLET); /* radical */
        
        if ( !bRadOrMultBonds && bIsAmmoniumSalt( at, i, &iO, &k, num_explicit_H ) ) {
            ;
        } else
        if ( !bRadOrMultBonds && bIsMetalSalt( at, i ) ) {
            ;
        } else
        if ( 1 == (j = bIsMetalToDisconnect(at, i, bCheckMetalValence)) ) {
            num_impl_H += NUMH(at,i);
            num_changes ++;
        } else
        if ( 2 == j && bTautFlagsDone ) {
            *bTautFlagsDone |= TG_FLAG_CHECK_VALENCE_COORD_DONE;
        }
    }
    if ( nNumImplH )
        *nNumImplH = num_changes? num_impl_H+1 : 0;
    return num_changes;
}



/*****************************************************************************/
#if ( bRELEASE_VERSION == 0 && (EXTR_HAS_METAL_ATOM & (EXTR_MASK | EXTR_FLAG) ) )
int bHasMetalAtom( ORIG_ATOM_DATA *orig_inp_data )
{
    int i;
    inp_ATOM *at;
    if ( orig_inp_data && (at   =  orig_inp_data->at) ) {
        int num_at     =  orig_inp_data->num_inp_atoms;
        /* check each atom */
        for ( i = 0; i < num_at; i ++ ) {
            if ( IS_METAL & get_el_type( at[i].el_number ) ) {
                return 1;
            }
        }
    }
    return 0;
}
#endif
/*****************************************************************************
{ "F",   19,  19,  18.998403220,     0 ,  0, {{0,},       {0,},       {1,},       {2,},       {3,5},      },},
{ "Cl",  35,  35,  34.968852730,     0 ,  0, {{0,},       {0,},       {1,3,5,7},  {2,4,6},    {3,5,},     },},
{ "Br",  80,  79,  78.918336100,     0 ,  0, {{0,},       {0,},       {1,3,5,7,}, {2,4,6,},   {3,5,},     },},
{ "I",  127, 127, 126.904500000,     0 ,  0, {{0,},       {0,},       {1,3,5,7,}, {2,4,6},    {3,5,},     },},
{ "At", 210, 210, 209.987100000,     0 ,  0, {{0,},       {0,},       {1,3,5,7,}, {2,4,6},    {3,5,},     },},
{ "N",   14,  14,  14.003074000,     0 ,  0, {{1,},       {2,},       {3,5},      {4,},       {3,},       },},
{ "P",   31,  31,  30.973762000,     0 ,  0, {{1,3,5,7,}, {2,4,6,},   {3,5,},     {4,},       {3,},       },},
{ "As",  75,  75,  74.921594200,     0 ,  0, {{0,},       {2,4,6,},   {3,5,},     {4,},       {3,},       },},
{ "Sb", 122, 121, 120.903800000,     0 ,  0, {{1,3,5,7,}, {2,4,6,},   {3,5,},     {2,4,},     {3,},       },},
{ "O",   16,  16,  15.994914630,     0 ,  0, {{0,},       {1,},       {2,},       {3,5,},     {4,},       },},
{ "S",   32,  32,  31.972070700,     0 ,  0, {{0,},       {1,3,5,7,}, {2,4,6},    {3,5,},     {4,},       },},
{ "Se",  79,  80,  79.916519600,     0 ,  0, {{0,},       {1,3,5,7,}, {2,4,6,},   {3,5,},     {4,},       },},
{ "Te", 128, 130, 129.906200000,     0 ,  0, {{0,},       {1,3,5,7,}, {2,4,6,},   {3,5,},     {2,4,},     },},
{ "Po", 209, 209, 208.982400000,     0 ,  0, {{0,},       {1,3,5,7,}, {2,4,6,},   {3,5,},     {2,4,},     },},
{ "B",   11,  11,  11.009300000,     0 ,  0, {{3,},       {4,},       {3,},       {2,},       {1,},       },},
*****************************************************************************/



int DisconnectMetals( ORIG_ATOM_DATA *orig_inp_data, int bCheckMetalValence, INCHI_MODE *bTautFlagsDone ) 
                           /*inp_ATOM *atom, int num_atoms, int nNumExplH, int *new_num_atoms */
{
    int i, j, k, n, iO, num_changes, val, bRadOrMultBonds;
    int num_impl_H, num_at, err, num_disconnected;
    S_CHAR num_explicit_H[NUM_H_ISOTOPES+1];
    static char elnumber_Heteroat[16] = {'\0', };
    static int  num_halogens;
    inp_ATOM  *at             = NULL;
    S_CHAR    *bMetal         = NULL;
    inp_ATOM  *atom           = orig_inp_data->at;
    int        num_atoms      = orig_inp_data->num_inp_atoms;
    int        nNumExplH      = (orig_inp_data->bDisconnectCoord > 0)? orig_inp_data->bDisconnectCoord - 1 : 0;
    AT_NUMB   *nOldCompNumber = orig_inp_data->nOldCompNumber;

    err              = 0;
    num_impl_H       = 0;
    num_at           = num_atoms;
    num_disconnected = 0;
    if ( !(at     = (inp_ATOM *)inchi_calloc( num_at + nNumExplH, sizeof(at[0]    ) )) ||
         !(bMetal = ( S_CHAR    *)inchi_calloc( num_at + nNumExplH, sizeof(bMetal[0]) )) ) {
        err = 1;
        goto exit_function;
    }
    if (!elnumber_Heteroat[0] ) {
        i = 0;
        /* halogens */
        elnumber_Heteroat[i++] = (char)get_periodic_table_number( "F"  ); /* 0 */
        elnumber_Heteroat[i++] = (char)get_periodic_table_number( "Cl" );
        elnumber_Heteroat[i++] = (char)get_periodic_table_number( "Br" );
        elnumber_Heteroat[i++] = (char)get_periodic_table_number( "I"  );
        elnumber_Heteroat[i++] = (char)get_periodic_table_number( "At" ); /* 4 */
        num_halogens = i;
        /* other non-metal */
        elnumber_Heteroat[i++] = (char)get_periodic_table_number( "N"  );
        elnumber_Heteroat[i++] = (char)get_periodic_table_number( "P"  );
        elnumber_Heteroat[i++] = (char)get_periodic_table_number( "As" );
        /*elnumber_Heteroat[i++] = get_periodic_table_number( "Sb" );*/ /* metal 10-28-2003 */
        elnumber_Heteroat[i++] = (char)get_periodic_table_number( "O"  );
        elnumber_Heteroat[i++] = (char)get_periodic_table_number( "S"  );
        elnumber_Heteroat[i++] = (char)get_periodic_table_number( "Se" );
        elnumber_Heteroat[i++] = (char)get_periodic_table_number( "Te" );
        /*elnumber_Heteroat[i++] = get_periodic_table_number( "Po" );*/ /* metal 10-28-2003 */
        elnumber_Heteroat[i++] = (char)get_periodic_table_number( "B"  );
        elnumber_Heteroat[i++] = 0;
    }

    memcpy( at, atom, num_atoms * sizeof(at[0]) );
    
    /* check each atom, mark metals */
    for ( i = 0, k = 0, num_changes = 0; i < num_atoms; i ++ ) {

        if ( !(val = at[i].valence) && !NUMH(at,i) ) {
            continue; /* disconnected atom */
        }
        bRadOrMultBonds = (val == 0) ||
             (val != at[i].chem_bonds_valence) || /* a bond has higher multiplicity than 1 */
             (at[i].radical && at[i].radical != RADICAL_SINGLET); /* radical */
        
        if ( !bRadOrMultBonds && bIsAmmoniumSalt( at, i, &iO, &k, num_explicit_H ) ) {
            ;
        } else
        if ( !bRadOrMultBonds && bIsMetalSalt( at, i ) ) {
            ;
        } else
        if ( 1 == (j = bIsMetalToDisconnect(at, i, bCheckMetalValence)) ) {
            num_impl_H += (k = NUMH(at,i));
            bMetal[i] = 1+k;
            num_changes ++;
        } else
        if ( 2 == j && bTautFlagsDone ) {
            *bTautFlagsDone |= TG_FLAG_CHECK_VALENCE_COORD_DONE;
        }
    }
    if ( num_impl_H != nNumExplH ) {
        err = 2;
        goto exit_function;
    }


    /* replace implicit H atoms with explicit H atoms */
    for ( i = 0; i < num_atoms && 0 < num_impl_H; i ++ ) {
        if ( bMetal[i] <= 1 ) {
            continue;
        }
        for ( k = 0; k < NUM_H_ISOTOPES+1; k ++ ) {
            n = k? at[i].num_iso_H[k-1] : at[i].num_H;
            for ( j = 0; j < n; j ++ ) {
                if ( num_at >= num_atoms + nNumExplH ) {
                    err = 3;
                    goto exit_function;
                }
                at[num_at].elname[0] = 'H';
                at[num_at].el_number = get_periodic_table_number(at[num_at].elname);
                at[num_at].iso_atw_diff = k;
                at[num_at].component = at[i].component;
                move_explicit_Hcation(at, num_at+1, i, num_at, 1);
                at[num_at].orig_at_number = num_at+1;
                num_at ++;
                num_impl_H --;
                bMetal[i] --;
                if ( k ) {
                    at[i].num_iso_H[k-1] --;
                } else {
                    at[i].num_H --;
                }
            }
        }
        if ( bMetal[i] != 1 ) {
            err = 4;
            goto exit_function;
        }
    }
    if ( num_at != num_atoms + nNumExplH ) {
        err = 5;
        goto exit_function;
    }

    /* disconnect metal - ligand bonds */
    for ( i = 0; i < num_atoms; i ++ ) {
        if ( !bMetal[i] ) {
            continue;
        }
        /* disconnect metal atom M

        Note: Defect in case of bridging ligands:

        M     M                          M     M             M     M(+)
          \  /   will be transformed to            , not to         
           N(+)                             N(+)                N(-)
          / \                              / \                 / \  
         R   R                            R   R               R   R
         
        Non-bridging are OK:

        M     R           M(+)  R
          \  /                 / 
           N(+)    --->      N   
          / \               / \  
         R   R             R   R
         
        */
        for ( j = at[i].valence-1; 0 <= j; j -- ) {
            if ( j < at[i].valence && !bMetal[ (int)at[i].neighbor[j] ] ) {
                /* do not break metal-metal bond here */
                num_disconnected += DisconnectOneLigand( at, nOldCompNumber, bMetal, elnumber_Heteroat,
                                                         num_halogens, num_atoms, i, j, bTautFlagsDone );
            }
        }
    }
    /* disconnect metal-metal bonds */
    for ( i = 0; i < num_atoms; i ++ ) {
        if ( !bMetal[i] ) {
            continue;
        }
        for ( j = at[i].valence-1; 0 <= j; j --  ) {
            if ( j < at[i].valence && bMetal[ (int)at[i].neighbor[j] ] ) {
                /* break metal-metal bond here */
                num_disconnected += DisconnectOneLigand( at, nOldCompNumber, bMetal, elnumber_Heteroat,
                                                         num_halogens, num_atoms, i, j, bTautFlagsDone );
            }
        }
    }


exit_function:
    if ( !num_disconnected ) {
        err = 6;
    }
    if ( at && err ) {
        inchi_free( at );
        at = NULL;
    }
    if ( atom && at ) { /* changed if ( at ) to if ( atom && at ) 2004-04-03 */
        inchi_free( atom );
        atom = NULL;
    }
    if ( bMetal )
        inchi_free( bMetal );
    
    if ( at ) {
        orig_inp_data->at = at;
        orig_inp_data->num_inp_atoms = num_at;
    }
    return err? -err : num_disconnected;
}



/*****************************************************************************/
int DisconnectOneLigand( inp_ATOM *at, AT_NUMB *nOldCompNumber, S_CHAR *bMetal, char *elnumber_Heteroat,
                         int num_halogens, int num_atoms, int iMetal, int jLigand, INCHI_MODE *bTautFlagsDone )
{
    int i, j, iLigand, neigh, val;
    int metal_neigh_ord[MAXVAL], num_neigh_arom_bonds[MAXVAL];
    int num_metal_neigh, num_disconnections;
    int num_del_arom_bonds, num_tot_arom_bonds, new_charge;
    char *p;

    iLigand = at[iMetal].neighbor[jLigand];
    num_metal_neigh    = 0;
    num_disconnections = 0;
    num_del_arom_bonds  = num_tot_arom_bonds = 0;

    /* find bonds to disconnect */
    for ( i = 0; i < at[iLigand].valence; i ++ ) {
        num_neigh_arom_bonds[i] = 0;
        neigh = (int)at[iLigand].neighbor[i];
        if ( neigh < num_atoms && bMetal[ neigh ] ) {
            metal_neigh_ord[ num_metal_neigh ++ ] = i;
            if ( at[iLigand].bond_type[i] > BOND_TYPE_TRIPLE ) {
                /* aromatic bond */
                for ( j = 0; j < at[neigh].valence; j ++ ) {
                    num_neigh_arom_bonds[i] += ( at[neigh].bond_type[j] > BOND_TYPE_TRIPLE );
                }
                num_del_arom_bonds ++;
            }
        }
        num_tot_arom_bonds += (at[iLigand].bond_type[i] > BOND_TYPE_TRIPLE);
    }
    /* Disconnect */
    if ( num_del_arom_bonds ) {
        /* fix chem_valence of the ligand and its neighbors in case of disconnecting arom. bonds */
        /* because in this case special care should be taken of updating at[].chem_bonds_valence */
        for ( i = 0; i < num_metal_neigh; i ++ ) {
            j = metal_neigh_ord[i];
            if ( num_neigh_arom_bonds[j] ) {
                neigh = at[iLigand].neighbor[j];
                at[neigh].chem_bonds_valence -= num_neigh_arom_bonds[j]/2 - (num_neigh_arom_bonds[j]-1)/2;
            }
        }
        at[iLigand].chem_bonds_valence -= num_tot_arom_bonds/2 - (num_tot_arom_bonds-num_del_arom_bonds)/2;
    }
    /* disconnect in reverse order, otherwise the metal_neigh_ord[i]
       becomes invalid after the first disconnection
    */
    for ( i = num_metal_neigh-1; 0 <= i; i -- ) {
        num_disconnections += DisconnectInpAtBond( at, nOldCompNumber, iLigand, metal_neigh_ord[i] );
    }

    /* attempt to change ligand charge to make its valence 'natural' */
    i = num_tot_arom_bonds - num_del_arom_bonds;
    if ( (i  && i != 2 && i != 3) ||
         (at[iLigand].radical && at[iLigand].radical != RADICAL_SINGLET) ||
         !(p = strchr( elnumber_Heteroat, at[iLigand].el_number ) ) ) {
        goto exit_function;  /* non-standard atom */
    }
    val = at[iLigand].chem_bonds_valence + NUMH(at, iLigand);
    new_charge = MAX_ATOMS; /* impossible value */
    if ( !val ) {
        if ( p - elnumber_Heteroat < num_halogens ) {
            new_charge = -1;
        }
    } else {
        for ( i = -1; i <= 1; i ++ ) {
            if ( val == get_el_valence( at[iLigand].el_number, i, 0 ) ) {
                new_charge = i; /* found charge that fits chem. valence */
                break;
            }
        }
    }
    if ( new_charge != MAX_ATOMS ) {
        if ( (new_charge != at[iLigand].charge ||
              (at[iLigand].radical && at[iLigand].radical != RADICAL_SINGLET)) &&
             1 == num_metal_neigh ) {
            if ( 1 == new_charge && 4 == val && 2 == at[iLigand].valence &&
                 4 == at[iLigand].chem_bonds_valence &&
                at[iLigand].bond_type[0] == at[iLigand].bond_type[1] ) {
                ; /* do not add +1 charge to disconnected =N=, etc. 2004-10-27 */
            } else {
                if ( bTautFlagsDone && new_charge != at[iLigand].charge ) {
                    *bTautFlagsDone |= TG_FLAG_MOVE_CHARGE_COORD_DONE;
                }
                at[iMetal].charge -= new_charge - at[iLigand].charge;
                at[iLigand].charge = new_charge;
                /*at[iLigand].radical = 0;*/
            }
        }
    }
exit_function:
    return num_disconnections; /* ret;*/
}




/****************************************************************************************/
double dist3D( inp_ATOM *at1, inp_ATOM *at2 )
{
    double dx = at1->x - at2->x;
    double dy = at1->y - at2->y;
    double dz = at1->z - at2->z;
    return sqrt( dx*dx+dy*dy+dz*dz );
}
/****************************************************************************************/
#define MIN_BOND_LENGTH   (1.0e-6)
#define MIN_COS           (1.0e-6)
#define MIN_BOND_LENGTH2  (MIN_BOND_LENGTH*MIN_BOND_LENGTH)
#define MAX_BOND_LENGTH   (1.0e30)
/****************************************************************************************/
double GetMinDistDistribution( inp_ATOM *at, int num_at, int iat, int iat_H,
                                int bInAllComponents, double min_dist[], int num_segm )
{
/*	const double one_pi = 2.0*atan2(1.0 , 0.0 ); */
    const double one_pi = 3.14159265358979323846; /* M_PI */
    const double two_pi = 2.0*one_pi;
    const double f_step = two_pi / num_segm;
    const double h_step = f_step/2.0;

    int i, j, k, kk, ki, kn, n, num_bonds;
    double xi, yi, xn, yn, cross_prod_in, dot_prod_in, xni, yni, rni, tni, rmin;
    double fi, fk, fn, ft, rt, rk, ri, rn, c, ave_bond_len;

    for ( i = 0; i < num_segm; i ++ ) {
        min_dist[i] = MAX_BOND_LENGTH; /* more than any distance */
    }
    num_bonds    = 0;
    ave_bond_len = 0.0;
    for ( i = 0; i < num_at; i ++ ) {
        if ( i != iat && i != iat_H && (bInAllComponents || at[i].component == at[iat].component) ) {
            for ( j = 0; j < at[i].valence; j ++ ) {
                n = at[i].neighbor[j];
                if ( (n > i && n != iat) || n == iat_H )
                    continue;
#if ( bRELEASE_VERSION != 1 && defined(_DEBUG) )
                if ( n == iat ) {
                    int stop = 1;  /* <BRKPT> */
                }
#endif
                xi = at[i].x - at[iat].x;  /* ri; i != iat */
                yi = at[i].y - at[iat].y;
                xn = at[n].x - at[iat].x;  /* rn; possibly n == iat */
                yn = at[n].y - at[iat].y;
                cross_prod_in = xi*yn - xn*yi; /* ((r(i)-r(iat)) x (r(n)-r(iat)) */
                if ( cross_prod_in < -0.01*MIN_BOND_LENGTH2 ) {
                    /* make sure the r(i)->r(n) vector is counterclockwise around at[iat] */
                    inchi_swap( (char*)&xi, (char*)&xn, sizeof(xi) );
                    inchi_swap( (char*)&yi, (char*)&yn, sizeof(yi) );
                    cross_prod_in = -cross_prod_in;
                }
                xni = xn - xi; /* r(n)->r(i) */
                yni = yn - yi;
                rni = xni*xni + yni*yni;
                if ( rni > 0.01*MIN_BOND_LENGTH2 ) {
                    /* vector length |ri->rn| is not too small */
                    /* arrowhead of the vector r(t) = ri + (rn-ri)*t; 0 <= t <= 1 points to the bond ri->rn */
                    /* r(tni) is perpendicular to the bond ri->rn so that min|r(t)| = r(tni) = |tni|*rni */
                    tni = -(xni*xi + yni*yi)/rni;
                    /* find min. distance from n-i bond to at[iat] */
                    if ( tni < 0.0 ) {
                        rmin = sqrt( xi*xi + yi*yi );
                    } else
                    if ( tni > 1.0 ) {
                        rmin = sqrt( xn*xn + yn*yn );
                    } else {
                        rmin = sqrt(tni*tni*rni);
                    }
                    ave_bond_len += sqrt( rni );
                    num_bonds ++;
                } else {
                    /* zero length i-n bond */
                    tni  = 0.5; /* fake */
                    rmin = sqrt( xi*xi + yi*yi ); /* arbitrarily choose one */
                }
                if ( rmin >= 0.1*MIN_BOND_LENGTH ) {
                    /* at[iat] does not belong to at[i]-at[n] bond */
                    int    bCalc_rt = 1;
                    fi = atan2( yi, xi );
                    fn = (n == iat)? fi : atan2( yn, xn );
                    if ( fi > fn ) {
                        /* make sure fn - fi >= 0 */
                        fn += two_pi;
                    }
                    if ( fi < 0.0 ) {
                        fi += two_pi;
                        fn += two_pi;
                    }
                    ki = (int)floor((fi+h_step)/f_step);  /* cast does not match function type */
                    kn = (int)floor((fn+h_step)/f_step);
                    /* the bond may affect several segments */
                    for ( k = ki; k <= kn; k ++ ) {
                        kk = k % num_segm;
                        if ( min_dist[kk] < rmin )
                            continue;
                        if ( bCalc_rt ) {
                            if ( n == iat ) {
                                ft = fi;
                                rt = rmin;
                            } else {
                                double xt, yt;
                                xt  = xi + xni*tni;
                                yt  = yi + yni*tni;
                                ft  = atan2( yt, xt );
                                rt  = sqrt(xt*xt + yt*yt);
                            }
                            bCalc_rt = 0;
                        }
                        fk = f_step * kk;
                        c  = fabs(cos( fk - ft ));
                        if ( c < MIN_COS )
                            c = MIN_COS;
                        rk = rt / c;
                        if ( min_dist[kk] > rk ) {
                            min_dist[kk] = rk;
                        }
                    }
                } else {
                    /* rmin < 0.1*MIN_BOND_LENGTH */
                    ri = xi*xi + yi*yi;
                    rn = xn*xn + yn*yn;
                    if ( ri > MIN_BOND_LENGTH2 && rn > MIN_BOND_LENGTH2 ) {
                        dot_prod_in = xn*xi + yn*yi;
                        /* a very short bond */
                        if ( dot_prod_in > 0.01*MIN_BOND_LENGTH2 ) {
                            /* bond does not cross at[iat] */
                            double fyixi = atan2( yi, xi );
                            if ( fyixi < 0.0 ) fyixi += two_pi;
                            kk = (int)floor((fyixi+h_step)/f_step) % num_segm; 
                            if ( min_dist[kk] > rmin ) {
                                min_dist[kk] = rmin;
                            }
                        } else
                        if ( dot_prod_in < -0.01*MIN_BOND_LENGTH2 ) {
                            /* bond does cross at[iat] */
                            double fyixi = atan2( yi, xi );
                            if ( fyixi < 0.0 ) fyixi += two_pi;
                            kk = (int)floor((fyixi+h_step)/f_step) % num_segm; 
                            if ( min_dist[kk] > rmin ) {
                                min_dist[kk] = rmin;
                            }
                            fyixi += one_pi;
                            kk = (int)floor((fyixi+h_step)/f_step) % num_segm; 
                            if ( min_dist[kk] > rmin ) {
                                min_dist[kk] = rmin;
                            }
                        } else {
                            ; /* error, should not happen */
                        }
                    } else
                    if ( ri <= MIN_BOND_LENGTH2 && rn <= MIN_BOND_LENGTH2 ) {
                        /* a very short bond coincides with at[iat]; ignore */
                        ;
                    } else {
                        /* one end of the bond coincides with at[iat] */
                        fi = ri>rn? atan2( yi, xi) : atan2( yn, xn );
                        if ( fi < 0.0 ) fi += two_pi;
                        kk = (int)floor((fi+h_step)/f_step) % num_segm; 
                        if ( min_dist[kk] > rmin ) {
                            min_dist[kk] = rmin;
                        }
                    }
                }
            }
        }
    }
    if ( num_bonds ) {
        return  ave_bond_len / (double)num_bonds;
    } else {
        return 0.0;
    }
}



/****************************************************************************************/
int move_explicit_Hcation(inp_ATOM *at, int num_at, int iat, int iat_H, int bInAllComponents)
{
#define NUM_SEGM 20
    /*	const double one_pi = 2.0*atan2(1.0 , 0.0 ); */
    const double one_pi = 3.14159265358979323846; /* M_PI */
    const double two_pi = 2.0*one_pi;
    const double f_step = two_pi / NUM_SEGM;
    const double h_step = f_step/2.0;
    double min_dist[NUM_SEGM];
    int nB, i, k, kk, next, val;
    double r, r0, xd, yd, zd, xr, yr, zr, ave_bond_len;
    /*double step = 4.0*atan(1.0)/NUM_SEGM;*/
    /* find at[iat] neighbors coordinates */
    xd=yd=zd=0.0;
    if ( at[iat].valence ) {
        for ( i = 0, nB=0, r = 0.0; i < at[iat].valence; i ++ ) {
            next = at[iat].neighbor[i];
            xd += at[next].x;
            yd += at[next].y;
            zd += at[next].z;
            r  += dist3D( at+iat, at+next );
            nB ++;
        }
        xd /= (double)nB;
        yd /= (double)nB;
        zd /= (double)nB;
        r  /= (double)nB;
        r0 = sqrt((double)(xd-at[iat].x)*(xd-at[iat].x)  
                + (double)(yd-at[iat].y)*(yd-at[iat].y));
    } else {
        if ( at[iat_H].valence ) {
            r = dist3D( at+iat_H, at+ (int)at[iat_H].neighbor[0] );
        } else {
            r = 0.0;
        }
        r0 = 0.0;
    }
    ave_bond_len = GetMinDistDistribution( at, num_at, iat, iat_H, bInAllComponents, min_dist, NUM_SEGM );
    if ( r < MIN_BOND_LENGTH && ave_bond_len > MIN_BOND_LENGTH ) {
            r = ave_bond_len; /* ave_bond_len = 0.0 may mean that it is 0D structure */
    }
    if ( r > MIN_BOND_LENGTH ) {
        /* process non-zero bond lengths */
        double f;
        if ( 10.0*r0 < r ) {
            xr =  -r;     /* arbitrary */
            yr =  0.0;
            zr =  0.0;
        } else {
            /*
            if ( r0 < MIN_BOND_LENGTH ) {
                r0 = 1.0;
            }
            */
            xr =   r * ( at[iat].x - xd )/r0;
            yr =   r * ( at[iat].y - yd )/r0; /* length = r */
            zr =   r * ( at[iat].z - zd )/r0;

/*          -- test: opposire direction --
            xr =   -r * ( at[iat].x - xd )/r0;
            yr =   -r * ( at[iat].y - yd )/r0;
            zr =   -r * ( at[iat].z - zd )/r0;
*/            
            if ( xr*xr + yr*yr < 0.04*r*r ) {
                xr = -r;
                yr = 0.0;
            }
        }
        r = sqrt( xr*xr + yr*yr );
        f = atan2( yr, xr );

        if ( f < 0.0 )
            f += two_pi;

        
        
        kk = (int)floor((f+h_step)/f_step) % NUM_SEGM; /* cast does not match function type by design */
        if ( min_dist[kk] < 1.5* r ) {
            double dist = 1.5*r;
            int start=-1, len=0, start_max=-1, len_max=0;
again:
            /* look for longest kk interval with min_dist[kk] >= dist */
            for ( k = 0, start = 0, len = 0, len_max = 0; k < 2*NUM_SEGM; k ++ ) {
                kk = k % NUM_SEGM;
                if ( min_dist[kk] >= dist ) {
                    if ( !len ++) {
                        start = k;
                    }
                } else {
                    if ( len > len_max ) {
                        len_max = len;
                        start_max = start;
                    }
                    len = 0;
                }
            }
            if ( !len_max ) {
                if ( dist > 0.1*r ) {
                    dist *= 0.75;
                    goto again;
                } else {
                    goto done; /* do it anyway */
                }
            } else {
                /* found a good sector */
                f = f_step * (start_max + (double)(len_max - 1)/2.0);
                r0 = dist / 1.5;
                xr = r0 * cos(f);
                yr = r0 * sin(f);
                zr = zr/r*r0;
            }
        }
    } else {
        xr = yr = zr = 0;
    }

done:

    if ( at[iat_H].valence ) {
        /* disconnect H */
        next = at[iat_H].neighbor[0];
        for ( i = 0; i < at[next].valence; i ++ ) {
            if ( at[next].neighbor[i] == iat_H ) {
                RemoveInpAtBond( at, next, i );
                i = 0; /* success */
                break;
            }
        }
    } else {
        /* isolated H+ cation */
        next = iat_H;
        i    = 0;
        at[iat_H].valence = 1;
        at[iat_H].chem_bonds_valence = 1;
        at[iat_H].bond_type[0] = BOND_TYPE_SINGLE;
    }
    if ( 0 == i /*i < at[next].valence*/ ) {
        /* move charge */
        if ( at[next].charge > 0 && at[iat].charge < 0 ) {
            at[next].charge --;
            at[iat].charge ++;
        }
        /* connect H to at[iat] */
        val = at[iat].valence;
        at[iat].neighbor[val] = iat_H;
        at[iat].bond_type[val] = at[iat_H].bond_type[0];
        at[iat].bond_stereo[val] = 0;
        at[iat].chem_bonds_valence += at[iat_H].bond_type[0];
        at[iat].valence = val+1;

        at[iat_H].component = at[iat].component;
        at[iat_H].neighbor[0] = iat;
        at[iat_H].bond_stereo[0] = 0; /* possible loss of stereo info */
        at[iat_H].x = at[iat].x + xr;
        at[iat_H].y = at[iat].y + yr;
        at[iat_H].z = at[iat].z + zr;
        return 1; /* success */
    }
    return 0; /* failed */
}
/****************************************************************************************/
int get_iat_number( int el_number, const int el_num[], int el_num_len )
{
    int i;
    for ( i = 0; i < el_num_len; i ++ ) {
        if ( el_num[i] == el_number )
            return i;
    }
    return -1;
}


/*#endif*/ /* } DISCONNECT_SALTS */
 typedef enum tagIonAtomType {
     IAT_H=0,
     IAT_C,
     IAT_N,
     IAT_P,
     IAT_O,
     IAT_S,
     IAT_Se,
     IAT_Te,
     IAT_F,
     IAT_Cl,
     IAT_Br,
     IAT_I,
     IAT_MAX
 } ION_ATOM_TYPE;

#if ( READ_INCHI_STRING == 1 )
/****************************************************************************************/
int bHeteroAtomMayHaveXchgIsoH( inp_ATOM *atom, int iat )
{
    inp_ATOM *at = atom + iat, *at2;
    static int el_num[IAT_MAX];
    int j, val, is_O=0, is_Cl=0, is_N=0, is_H=0, num_H, iat_numb, bAccept, cur_num_iso_H;
    
    if ( !el_num[IAT_H]) {
        el_num[IAT_H ] = get_periodic_table_number( "H" ); 
        el_num[IAT_C ] = get_periodic_table_number( "C" ); 
        el_num[IAT_N ] = get_periodic_table_number( "N" ); 
        el_num[IAT_P ] = get_periodic_table_number( "P" ); 
        el_num[IAT_O ] = get_periodic_table_number( "O" ); 
        el_num[IAT_S ] = get_periodic_table_number( "S" ); 
        el_num[IAT_Se] = get_periodic_table_number( "Se"); 
        el_num[IAT_Te] = get_periodic_table_number( "Te"); 
        el_num[IAT_F ] = get_periodic_table_number( "F" ); 
        el_num[IAT_Cl] = get_periodic_table_number( "Cl"); 
        el_num[IAT_Br] = get_periodic_table_number( "Br"); 
        el_num[IAT_I ] = get_periodic_table_number( "I" ); 
    }
    if ( 0 > (iat_numb = get_iat_number( at->el_number, el_num, IAT_MAX )) ) {
        return 0;
    }
    if ( abs(at->charge) > 1 || (at->radical && RADICAL_SINGLET != at->radical) ) {
        return 0;
    }
    val = -1;
    switch( iat_numb ) {
    case IAT_N:
    case IAT_P:
        is_N = 1;
        val  = 3+at->charge;
        break;
    case IAT_O:
    case IAT_S:
    case IAT_Se:
    case IAT_Te:
        is_O = 1;
        val  = 2+at->charge;
        break;
    case IAT_F:
    case IAT_Cl:
    case IAT_Br:
    case IAT_I:
        if ( at->charge == 0 ) {
            is_Cl = 1; /* isolated HCl */
            val   = 1;
        }
        break;
    case IAT_H:
        if ( at->valence   == 0 &&
             at->charge    == 1 ) {
            is_H = 1; /* isolated proton */
            val  = 0;
        }
    }
    if ( val < 0 ) {
        return 0;
    }
    num_H = NUMH(at,0);
    if ( val != at->chem_bonds_valence + num_H ) {
        return 0;
    }
    if ( is_H ) {
        return 2; /* H atom */
    } else {
        cur_num_iso_H = 0;
        for ( j = 0, bAccept = 1; j < at->valence && bAccept; j ++ ) {
            at2 = atom + (int)at->neighbor[j];
            if ( (at2->charge && at->charge) || (at2->radical && RADICAL_SINGLET != at2->radical ) ) {
                return 0; /* adjacent charged/radical atoms: do not neutralizate */
            }
        }
    }
    return 1;
}
#endif
/****************************************************************************************/
int bNumHeterAtomHasIsotopicH( inp_ATOM *atom, int num_atoms )
{
    static int el_num[IAT_MAX];
    int i, j, val, is_O=0, is_Cl=0, is_N=0, is_H=0, num_H, iat_numb, bAccept, num_iso_H, cur_num_iso_H, num_iso_atoms;
    inp_ATOM *at, *at2;
    /* one time initialization */
    if ( !el_num[IAT_H]) {
        el_num[IAT_H ] = get_periodic_table_number( "H" ); 
        el_num[IAT_C ] = get_periodic_table_number( "C" ); 
        el_num[IAT_N ] = get_periodic_table_number( "N" ); 
        el_num[IAT_P ] = get_periodic_table_number( "P" ); 
        el_num[IAT_O ] = get_periodic_table_number( "O" ); 
        el_num[IAT_S ] = get_periodic_table_number( "S" ); 
        el_num[IAT_Se] = get_periodic_table_number( "Se"); 
        el_num[IAT_Te] = get_periodic_table_number( "Te"); 
        el_num[IAT_F ] = get_periodic_table_number( "F" ); 
        el_num[IAT_Cl] = get_periodic_table_number( "Cl"); 
        el_num[IAT_Br] = get_periodic_table_number( "Br"); 
        el_num[IAT_I ] = get_periodic_table_number( "I" ); 
    }
    num_iso_H     = 0;
    num_iso_atoms = 0;
    for ( i = 0, at = atom; i < num_atoms; i ++, at ++ ) {

        num_iso_atoms += ( at->iso_atw_diff != 0 || NUM_ISO_H(at,0) ); /* isotopic atoms and implicit isotopic H */

        if ( 0 > (iat_numb = get_iat_number( at->el_number, el_num, IAT_MAX )) ) {
            continue;
        }
        
        if ( abs(at->charge) > 1 || (at->radical && RADICAL_SINGLET != at->radical) ) {
            continue;
        }

        val = -1;
        switch( iat_numb ) {
        case IAT_N:
        case IAT_P:
            is_N = 1;
            val  = 3+at->charge;
            break;
        case IAT_O:
        case IAT_S:
        case IAT_Se:
        case IAT_Te:
            is_O = 1;
            val  = 2+at->charge;
            break;
        case IAT_F:
        case IAT_Cl:
        case IAT_Br:
        case IAT_I:
            if ( at->charge == 0 ) {
                is_Cl = 1; /* isolated HCl */
                val   = 1;
            }
            break;
        case IAT_H:
            if ( at->valence   == 0 &&
                 at->charge    == 1 ) {
                is_H = 1; /* isolated proton */
                val  = 0;
            }
        }
        if ( val < 0 ) {
            continue;
        }
        num_H = NUMH(at,0);
        if ( val != at->chem_bonds_valence + num_H ) {
            continue;
        }
        if ( is_H ) {
            bAccept = 1;
            cur_num_iso_H = (at->iso_atw_diff != 0);
        } else {
            cur_num_iso_H = 0;
            for ( j = 0, bAccept = 1; j < at->valence && bAccept; j ++ ) {
                at2 = atom + (int)at->neighbor[j];
                if ( (at2->charge && at->charge) || (at2->radical && RADICAL_SINGLET != at2->radical ) ) {
                    bAccept = 0; /* adjacent charged/radical atoms: do not neutralizate */
                    break;
                } else
                if ( at2->el_number == el_num[IAT_H ] && at2->valence == 1 && at2->iso_atw_diff ) {
                    cur_num_iso_H ++; /* isotopic explicit H */
                }
            }
            if ( bAccept ) {
                num_iso_atoms -= cur_num_iso_H;  /* avoid counting explicit H as isotopic atom */
                cur_num_iso_H += NUM_ISO_H(at,0);
            }
            
        }
        num_iso_H += (bAccept && cur_num_iso_H); /* number of acceptable heteroatoms that have isotopic H */
    }
    return ((num_iso_H? 1:0) | (num_iso_atoms? 2:0));
}


/****************************************************/
/* Mark and count disconnected structure components */
/* by Depth-first searching each component          */
/****************************************************/
int cmp_components( const void *a1, const void *a2 )
{
    int ret;
    AT_NUMB n1;
    AT_NUMB n2;

    n1 = ((const AT_NUMB *)a1)[0];     /* number of atoms in the component -- descending order */
    n2 = ((const AT_NUMB *)a2)[0];
    if ( (ret = (int)n2 - (int)n1) ) {
        return ret;
    }
    /* stable sort */
    n1 = ((const AT_NUMB *)a1)[1];    /* component ordering number -- ascending order */
    n2 = ((const AT_NUMB *)a2)[1];
    ret = (int)n1 - (int)n2;

    return ret;
    
}
/*************************************************************************************************/
int MarkDisconnectedComponents( ORIG_ATOM_DATA *orig_at_data, int bProcessOldCompNumbers )
{
    typedef AT_NUMB AT_TRIPLE[3];

    inp_ATOM  *at                  = orig_at_data->at;
    int        num_at              = orig_at_data->num_inp_atoms;
    AT_NUMB *nCurAtLen           = NULL;
    
    AT_NUMB *nNewCompNumber      = NULL;
    AT_NUMB *nPrevAtom           = NULL;
    S_CHAR  *iNeigh              = NULL;

    AT_NUMB *nOldCompNumber      = NULL;
    int i, j, num_components, ret;
    int new_comp_no;
    AT_NUMB old_comp_no, another_comp_no, no_component;

    /* component_nbr[i][0] = number of atoms in the component i-1
     * component_nbr[i][1] = original component number (id-1) = i
     * after sorting:
     * component_nbr[j][2] = new number of component #(component_nbr[i][1]+1)
     */
    AT_TRIPLE *component_nbr = NULL;

    /* initialize */
    if ( bProcessOldCompNumbers && !orig_at_data->nOldCompNumber ) {
        bProcessOldCompNumbers = 0; 
    }
    num_components = 0;
    /*
    for ( j = 0; j < num_at; j ++ ) {
        at[j].component = 0;
    }
    */
    ret = -1;
    if ( !num_at ) {
        return 0;
    }
    if ( !( nNewCompNumber = (AT_NUMB *) inchi_calloc( num_at, sizeof(nNewCompNumber[0]) ) ) ||
         /* for non-recursive DFS only: */
         !( nPrevAtom      = (AT_NUMB *) inchi_calloc( num_at, sizeof(nPrevAtom[0]) ) ) ||
         !( iNeigh         = (S_CHAR  *) inchi_calloc( num_at, sizeof(iNeigh[0]) ) )) {
        goto exit_function;
    }
    /* mark and count; avoid deep DFS recursion: it may make verifying software unhappy */
    /* nNewCompNumber[i] will contain new component number for atoms at[i], i=0..num_at-1 */
    for ( j = 0; j < num_at; j++ ) {
        if ( !nNewCompNumber[j] ) {
            /* mark starting with at[j] */
            int fst_at, nxt_at, cur_at = j;
            num_components ++;
            /* first time at at[j] */
            nNewCompNumber[fst_at = cur_at] = (AT_NUMB) num_components;
            /* find next neighbor */
            while ( 1 ) {
                if ( iNeigh[cur_at] < at[cur_at].valence ) {
                    nxt_at = at[cur_at].neighbor[(int)iNeigh[cur_at] ++];
                    if ( !nNewCompNumber[nxt_at] ) {
                        /* forward edge: found new atom */
                        nNewCompNumber[nxt_at] = (AT_NUMB) num_components;
                        nPrevAtom[nxt_at]      = (AT_NUMB) cur_at;
                        cur_at = nxt_at;
                    }
                } else
                if ( cur_at == fst_at ) {
                    break; /* done */
                } else {
                    cur_at = nPrevAtom[cur_at]; /* retract */
                }
            }
        }
    }
    inchi_free( nPrevAtom ); nPrevAtom = NULL;
    inchi_free( iNeigh );    iNeigh    = NULL;

    /* Allocate more memory */
    i = inchi_max( num_components, orig_at_data->num_components ); 
    if ( !(nCurAtLen      = (AT_NUMB *) inchi_calloc( num_components+1, sizeof(nCurAtLen[0]) ) ) ||
         !(nOldCompNumber = (AT_NUMB *) inchi_calloc( i             +1, sizeof(nOldCompNumber[0]) ) ) ||
         !(component_nbr  = (AT_TRIPLE *) inchi_calloc( num_components+1, sizeof(component_nbr[0]) ) ) ) {
        goto exit_function;
    }

    /* count atoms per component and renumber the components */
    for ( i = 0; i < num_components; i ++ ) {
        component_nbr[i][0] = 0; /* number of atoms in the component */
        component_nbr[i][1] = i; /* component ordering number */
    }
    for ( j = 0; j < num_at; j ++ ) {
        component_nbr[(int)nNewCompNumber[j]-1][0] ++; /* count atoms in each component */
    }
    /* sort key: number of atoms; order: descending */
    qsort( (void*)component_nbr[0], num_components,
           sizeof(component_nbr[0]), cmp_components);
    /* invert the transposition */
    for ( i = 0; i < num_components; i ++ ) {
        nCurAtLen[i] = component_nbr[i][0];
        component_nbr[ component_nbr[i][1] ][2] = i+1;
    }
    /* renumber the components so that the component with the greatest number of atoms is the first */
    no_component = num_at+1;
    for ( j = 0; j < num_at; j ++ ) {
        /* new component number for at[j] */
        new_comp_no = component_nbr[(int)nNewCompNumber[j]-1][2]-1; /* starts from 0 */
        if ( bProcessOldCompNumbers ) {
            /* old component number for at[j] */
            old_comp_no = at[j].component;
            /* fill out nOldCompNumber[]; initially it contains zeroes */
            if ( !old_comp_no ) {
                nOldCompNumber[new_comp_no] = no_component; /* atom did not have component number */
            } else
            if ( nOldCompNumber[new_comp_no] !=  old_comp_no ) {
                if ( !nOldCompNumber[new_comp_no] ) {
                    nOldCompNumber[new_comp_no] = old_comp_no;
                } else {
                    /* at[j] moved from old comp #old_comp_no to old comp #nOldCompNumber[new_comp_no]
                       Both components cannot be equal to any current component */
                    another_comp_no = nOldCompNumber[new_comp_no];
                    for ( i = 0; i < num_components; i ++ ) {
                        if ( nOldCompNumber[i] == old_comp_no ||
                             nOldCompNumber[i] == another_comp_no ) {
                            nOldCompNumber[i] = no_component;
                        }
                    }
                    /* nOldCompNumber[new_comp_no] = num_at+1; */
                }
            }
        }
        /* orig_at_data->nOldCompNumber */
        at[j].component = new_comp_no+1;  /* starts from 1 */
    }
    if ( bProcessOldCompNumbers ) { 
        for ( j = 0; j < num_components; j ++ ) {
            if ( nOldCompNumber[j] == no_component ) {
                /* the component has atom from another component */
                nOldCompNumber[j] = 0;
            } else
            if ( nOldCompNumber[j] && 
                 !orig_at_data->nOldCompNumber[nOldCompNumber[j]-1] ) {
                /* the component has changed in the previous processing  */
                nOldCompNumber[j] = 0;
            }
        }
    } else {
        for ( j = 0; j < num_components; j ++ ) {
            nOldCompNumber[j] = j + 1;
        }
    }
    ret = num_components;
exit_function:
    if ( nNewCompNumber )
        inchi_free( nNewCompNumber );
    if ( component_nbr )
        inchi_free( component_nbr );

    if ( ret < 0 ) {
        if ( nPrevAtom ) {
            inchi_free( nPrevAtom );
            nPrevAtom = NULL;
        }
        if ( iNeigh ) {
            inchi_free( iNeigh );
            iNeigh = NULL;
        }
        if ( nCurAtLen ) {
            inchi_free( nCurAtLen );
            nCurAtLen = NULL;
        }
        if ( nOldCompNumber ) {
            inchi_free( nOldCompNumber );
            nOldCompNumber = NULL;
        }
        num_components = ret;
    }
    /* avoid memory leaks */
    if ( orig_at_data->nCurAtLen )
        inchi_free ( orig_at_data->nCurAtLen );
    if ( orig_at_data->nOldCompNumber )
        inchi_free ( orig_at_data->nOldCompNumber );

    orig_at_data->nCurAtLen      = nCurAtLen;
    orig_at_data->nOldCompNumber = nOldCompNumber;

    orig_at_data->num_components = num_components;

    return ret;  /* number of disconnected components; 1=>single connected structure*/
}
/******************************************************************************/
/*                        Extract one (connected) component                   */
/******************************************************************************/
int ExtractConnectedComponent(  inp_ATOM *at, int num_at, int component_number, inp_ATOM *component_at )
{
    int i, j, num_component_at;
    AT_NUMB *number;
    if ( NULL == (number = (AT_NUMB*)inchi_calloc(num_at, sizeof(AT_NUMB)))){
        return CT_OUT_OF_RAM; /* out of memory */  /*   <BRKPT> */
    }
    /* copy atoms */
    for ( i = 0, num_component_at = 0; i < num_at; i ++ ) {
        if ( at[i].component == component_number ) {
            number[i] = num_component_at;
            component_at[num_component_at ++] = at[i];
        }
    }
    /* renumber neighbors */
    for ( i = 0; i < num_component_at; i ++ ) {
        component_at[i].orig_compt_at_numb = (AT_NUMB)(i + 1);
        for ( j = 0; j < component_at[i].valence; j ++ ) {
            component_at[i].neighbor[j] = number[(int)component_at[i].neighbor[j]];
        }
    }
    inchi_free( number );
    return num_component_at;
}
/****************************************************************/
int SetConnectedComponentNumber( inp_ATOM *at, int num_at, int component_number )
{
    int i;
    for ( i = 0; i < num_at; i ++ ) {
        at[i].component = (AT_NUMB)component_number;
    }
    return 0;
}
/****************************************************************/

int Free_INChI_Stereo( INChI_Stereo *pINChI_Stereo )
{
    if ( pINChI_Stereo ) {
        qzfree( pINChI_Stereo->nNumber );
        qzfree( pINChI_Stereo->t_parity );
        qzfree( pINChI_Stereo->nNumberInv );
        qzfree( pINChI_Stereo->t_parityInv );
        qzfree( pINChI_Stereo->nBondAtom1 );
        qzfree( pINChI_Stereo->nBondAtom2 );
        qzfree( pINChI_Stereo->b_parity );
    }
    return 0;
}

/****************************************************************/
INChI_Stereo *Alloc_INChI_Stereo(int num_at, int num_bonds)
{
    INChI_Stereo *pINChI_Stereo = (INChI_Stereo *)inchi_calloc(1, sizeof(INChI_Stereo));
    if ( pINChI_Stereo ) {
        if ( num_at &&
             (pINChI_Stereo->nNumber  = (AT_NUMB *)inchi_calloc(num_at, sizeof(pINChI_Stereo->nNumber[0]))) &&
             (pINChI_Stereo->t_parity = (S_CHAR  *)inchi_calloc(num_at, sizeof(pINChI_Stereo->t_parity[0]))) 
             && (pINChI_Stereo->nNumberInv  = (AT_NUMB *)inchi_calloc(num_at, sizeof(pINChI_Stereo->nNumberInv[0])))
             && (pINChI_Stereo->t_parityInv = (S_CHAR  *)inchi_calloc(num_at, sizeof(pINChI_Stereo->t_parityInv[0]))) 
             ) {
            ;
        } else
        if ( num_at ) {
            goto out_of_RAM;
        }
        if ( num_bonds &&
             (pINChI_Stereo->nBondAtom1 =(AT_NUMB *)inchi_calloc(num_bonds, sizeof(pINChI_Stereo->nBondAtom1[0]))) &&
             (pINChI_Stereo->nBondAtom2 =(AT_NUMB *)inchi_calloc(num_bonds, sizeof(pINChI_Stereo->nBondAtom2[0]))) &&
             (pINChI_Stereo->b_parity   =(S_CHAR  *)inchi_calloc(num_bonds, sizeof(pINChI_Stereo->b_parity[0]))) ) {
            ;
        } else
        if ( num_bonds ) {
            goto out_of_RAM;
        }
        return pINChI_Stereo;

out_of_RAM:
        Free_INChI_Stereo( pINChI_Stereo );
        qzfree( pINChI_Stereo );
    }
    return NULL;
}
/****************************************************************/
int Free_INChI(INChI **ppINChI)
{
    INChI *pINChI;

    if ( (pINChI = *ppINChI) ) {
#if ( bREUSE_INCHI == 1 )
        if ( pINChI->nRefCount -- > 0 )
            return 1;
#endif
        Free_INChI_Members(pINChI);

        qzfree( pINChI );
        *ppINChI = NULL;

    }
    return 0;
}
/****************************************************************/
int Free_INChI_Members(INChI *pINChI)
{
    if ( pINChI ) {
        
        Free_INChI_Stereo(pINChI->Stereo           );
        Free_INChI_Stereo(pINChI->StereoIsotopic   );
        qzfree(pINChI->nAtom                    );
        qzfree(pINChI->nConnTable               );
        qzfree(pINChI->nTautomer                );
        qzfree(pINChI->nNum_H                   );
        qzfree(pINChI->nNum_H_fixed             );
        qzfree(pINChI->IsotopicAtom             );
        qzfree(pINChI->IsotopicTGroup           );
        qzfree(pINChI->nPossibleLocationsOfIsotopicH);
        qzfree(pINChI->Stereo           );
        qzfree(pINChI->StereoIsotopic   );
        qzfree(pINChI->szHillFormula );
    }
    return 0;
}
/****************************************************************/
INChI *Alloc_INChI( inp_ATOM *at, int num_at, int *found_num_bonds, int *found_num_isotopic, int nAllocMode )
{
    int    i, num_bonds, num_isotopic_atoms;
    INChI  *pINChI;
    int    bIsotopic   = (nAllocMode & REQ_MODE_ISO);
    /* int    bTautomeric = (nAllocMode & REQ_MODE_TAUT); */

    if ( num_at <= 0 || NULL == (pINChI = (INChI *)inchi_calloc( 1, sizeof(INChI)))) {
        return NULL;
    }

    for ( i = 0, num_bonds = 0, num_isotopic_atoms = 0; i < num_at; i ++ ) {
        num_bonds += at[i].valence;
        /* if ( bIsotopic ) { */
            num_isotopic_atoms += (0 != at[i].iso_atw_diff ||
                                   !strcmp(at[i].elname, "D") ||
                                   !strcmp(at[i].elname, "T") ||
                                    at[i].num_iso_H[0] ||
                                    at[i].num_iso_H[1] ||
                                    at[i].num_iso_H[2]);
        /* } */

    }
    num_bonds /= 2;

    *found_num_bonds    = num_bonds;
    *found_num_isotopic = num_isotopic_atoms;


    if ( (pINChI->nAtom       = (U_CHAR*) inchi_calloc( num_at, sizeof(pINChI->nAtom[0]))) &&
         (pINChI->nConnTable  = (AT_NUMB*)inchi_calloc( num_at+num_bonds, sizeof(pINChI->nConnTable[0]))) &&
         (pINChI->nTautomer   = (AT_NUMB*)inchi_calloc( ((3+INCHI_T_NUM_MOVABLE)*num_at)/2+1, sizeof(pINChI->nTautomer[0]))) &&
         (pINChI->nNum_H      = (S_CHAR*) inchi_calloc( num_at, sizeof(pINChI->nNum_H[0]))) &&
         (pINChI->nNum_H_fixed= (S_CHAR*) inchi_calloc( num_at, sizeof(pINChI->nNum_H_fixed[0])))
         ) {
        ; /* nTautomer length: max. number of tautomeric groups is num_at/2
           
           1 word                     -> number of t-groups
           
           each group has:

           1 word                     -> number of endpoints+INCHI_T_NUM_MOVABLE
           INCHI_T_NUM_MOVABLE words   -> number(s) of moveable attachments
           numbers of endpoints words -> canon. numbers

           max. occurs if each t-group has 2 atoms (num_at/2 t-groups) and all atoms
                          belong to t-groups (num_at endpoints)

           Total: 1 + (number of t-groups)*(1+INCHI_T_NUM_MOVABLE) + (number of endpoints) <=
                  1 + (num_at/2) * (1+INCHI_T_NUM_MOVABLE) + num_at <=
                  1 + (3+INCHI_T_NUM_MOVABLE)*num_at/2 words.
           */

    } else {
        goto out_of_RAM;
    }
    pINChI->szHillFormula = NULL; /*  the length is unknown */
    if ( bIsotopic ) {
        if ( num_isotopic_atoms &&
             (pINChI->IsotopicAtom              = (INChI_IsotopicAtom *)inchi_calloc(num_isotopic_atoms, sizeof(INChI_IsotopicAtom) )) &&
             (pINChI->IsotopicTGroup            = (INChI_IsotopicTGroup *)inchi_calloc(num_isotopic_atoms, sizeof(INChI_IsotopicTGroup) ))
            ) {
            ;
        } else
        if ( num_isotopic_atoms ) {
            goto out_of_RAM;
        }
        if ( !(pINChI->nPossibleLocationsOfIsotopicH = (AT_NUMB *)inchi_calloc( num_at+1, sizeof(pINChI->nPossibleLocationsOfIsotopicH[0]) ) ) ) {
            goto out_of_RAM;
        }
    }

    if ((pINChI->Stereo            = Alloc_INChI_Stereo(num_at, num_bonds))
       ) {
        ;
    } else {
        goto out_of_RAM;
    }
    if ( bIsotopic ) {
        if ((pINChI->StereoIsotopic    = Alloc_INChI_Stereo(num_at, num_bonds))
           ) {
            ;
        } else {
            goto out_of_RAM;
        }
    }
    return pINChI;

out_of_RAM:
    if ( pINChI ) {
        Free_INChI(&pINChI);
        /*
        inchi_free(pINChI);
        */
    }
    return NULL;
}
/****************************************************************/
int Free_INChI_Aux( INChI_Aux **ppINChI_Aux )
{
    INChI_Aux *pINChI_Aux = *ppINChI_Aux;
    if ( pINChI_Aux ) {
#if ( bREUSE_INCHI == 1 )
        if ( pINChI_Aux->nRefCount -- > 0 )
            return 1;
#endif

        qzfree( pINChI_Aux->nOrigAtNosInCanonOrd            );
        qzfree( pINChI_Aux->nIsotopicOrigAtNosInCanonOrd    );
        qzfree( pINChI_Aux->nOrigAtNosInCanonOrdInv         );
        qzfree( pINChI_Aux->nIsotopicOrigAtNosInCanonOrdInv );
        qzfree( pINChI_Aux->szOrigCoord );
        qzfree( pINChI_Aux->OrigInfo );
/*
        qzfree( pINChI_Aux->nOriginalAtomNumber          );
        qzfree( pINChI_Aux->nCanonicalTGroupNumbers      );
        qzfree( pINChI_Aux->nIsotopicCanonicalTGroupNumbers);
        qzfree( pINChI_Aux->nTautomer                    );
        qzfree( pINChI_Aux->nNontautomericCanonicalNumbers         );
        qzfree( pINChI_Aux->nIsotopicCanonicalNumbers    );
        qzfree( pINChI_Aux->nNontautomericIsotopicCanonicalNumbers );
        qzfree( pINChI_Aux->nNontautomericEquNumbers               );
        qzfree( pINChI_Aux->nNontautomericIsotopicEquNumbers       );
*/
        qzfree( pINChI_Aux->nConstitEquNumbers                  );
        qzfree( pINChI_Aux->nConstitEquTGroupNumbers             );
        qzfree( pINChI_Aux->nConstitEquIsotopicNumbers          );
        qzfree( pINChI_Aux->nConstitEquIsotopicTGroupNumbers     );
        qzfree( pINChI_Aux );
        *ppINChI_Aux = NULL;
    }
    return 0;
}
/****************************************************************/
INChI_Aux *Alloc_INChI_Aux( int num_at, int num_isotopic_atoms, int nAllocMode, int bOrigCoord )
{
    INChI_Aux     *pINChI_Aux;
    int    bIsotopic   = (nAllocMode & REQ_MODE_ISO);
    int    num_at_tg   = num_at + num_at/2;
    /* int    bTautomeric = (nAllocMode & REQ_MODE_TAUT); */

    if ( num_at <= 0 || NULL == (pINChI_Aux = (INChI_Aux *)inchi_calloc(sizeof(INChI_Aux), 1))) {
        return NULL;
    }
    if ( (pINChI_Aux->nOrigAtNosInCanonOrd      = (AT_NUMB*)inchi_calloc(sizeof(pINChI_Aux->nOrigAtNosInCanonOrd[0]), num_at_tg)) &&
         (pINChI_Aux->nOrigAtNosInCanonOrdInv   = (AT_NUMB*)inchi_calloc(sizeof(pINChI_Aux->nOrigAtNosInCanonOrd[0]), num_at_tg)) &&
         (pINChI_Aux->nConstitEquNumbers        = (AT_NUMB*)inchi_calloc(sizeof(pINChI_Aux->nConstitEquNumbers[0]), num_at_tg)) ) {
        ;
    } else {
        goto out_of_RAM;
    }

    if ( num_at > 1 &&
         (pINChI_Aux->nConstitEquTGroupNumbers  = (AT_NUMB*)inchi_calloc(sizeof(pINChI_Aux->nConstitEquTGroupNumbers[0]), num_at/2+1)) ) {
        ;
    } else
    if ( num_at > 1 ) {
        goto out_of_RAM;
    }
    
    if ( num_at > 0 ) {
        pINChI_Aux->OrigInfo    = (ORIG_INFO *)inchi_calloc(sizeof(pINChI_Aux->OrigInfo[0]), num_at);
        if ( !pINChI_Aux->OrigInfo )
            goto out_of_RAM;
    }
    if ( bOrigCoord && num_at > 0 ) {
        pINChI_Aux->szOrigCoord = (MOL_COORD *)inchi_calloc(sizeof(pINChI_Aux->szOrigCoord[0]), num_at);
        if ( !pINChI_Aux->szOrigCoord )
            goto out_of_RAM;
    }
    if ( bIsotopic ) {
        if ( /*num_isotopic_atoms &&*/
             (pINChI_Aux->nIsotopicOrigAtNosInCanonOrd    = (AT_NUMB*)inchi_calloc(sizeof(pINChI_Aux->nIsotopicOrigAtNosInCanonOrd[0]), num_at_tg)) &&
             (pINChI_Aux->nIsotopicOrigAtNosInCanonOrdInv = (AT_NUMB*)inchi_calloc(sizeof(pINChI_Aux->nIsotopicOrigAtNosInCanonOrd[0]), num_at_tg)) &&
             (pINChI_Aux->nConstitEquIsotopicNumbers      = (AT_NUMB*)inchi_calloc(sizeof(pINChI_Aux->nConstitEquIsotopicNumbers[0]), num_at_tg)) ) {
            ;
        } else
        if ( num_isotopic_atoms ) {
            goto out_of_RAM;
        }
        if ( /*num_isotopic_atoms && num_at > 1 &&*/
             (pINChI_Aux->nConstitEquIsotopicTGroupNumbers = (AT_NUMB*)inchi_calloc(sizeof(pINChI_Aux->nConstitEquIsotopicTGroupNumbers[0]), num_at/2+1)) ) {
            ;
        } else
        if ( num_isotopic_atoms && num_at > 1 ) {
            goto out_of_RAM;
        }
    }
    return pINChI_Aux;


out_of_RAM:
    if ( pINChI_Aux ) {
        Free_INChI_Aux(&pINChI_Aux);
        /*
        inchi_free(pINChI_Aux);
        */
    }
    return NULL;
}
/***********************************************************************************/

#define IS_DEUTERIUM(i) (!strcmp( at[i].elname, "D" ) || at[i].iso_atw_diff == 2 && !strcmp( at[i].elname, "H" ))
#define IS_TRITIUM(i)   (!strcmp( at[i].elname, "T" ) || at[i].iso_atw_diff == 3 && !strcmp( at[i].elname, "H" ))

#define ABNORMAL_ISO(i) (at[i].iso_atw_diff == 1 || at[i].iso_atw_diff < -3 || at[i].iso_atw_diff > 5 )
#define ABNORMAL_CHG(i) (abs(at[i].charge) > 3)
#define ABNORMAL_RAD(i) (RADICAL_SINGLET <= at[i].radical && at[i].radical <= RADICAL_TRIPLET )

#define ANY_ISO(i, X)   ((X)? (at[i].iso_atw_diff && !IS_DEUTERIUM(i) && !IS_TRITIUM(i)) :\
                              (at[i].iso_atw_diff ||  IS_DEUTERIUM(i) ||  IS_TRITIUM(i)))
#define ANY_CHG(i)      (0 != at[i].charge)
#define ANY_RAD(i)      (RADICAL_SINGLET <= at[i].radical && at[i].radical <= RADICAL_TRIPLET )

#define NORMAL_ISO(i, X)   (ANY_ISO(i, X) && !ABNORMAL_ISO(i))


/* needs additional M  CHG. M  RAD, M  ISO line */
/* due to ISIS/Draw feature always include M  RAD for any radical */
#define ABNORMAL_AT(i) ( at[i].radical || abs(at[i].charge) > 3 || \
                         ABNORMAL_ISO(i) )

/* always add M  ISO, M  RAD, M  CHG; Except: (bAtomsDT && D or T) */
#define ADD_LINE_AT(i) ( at[i].charge  || \
                         at[i].radical || \
                         at[i].iso_atw_diff && (bAtomsDT? (at[i].iso_atw_diff != 1 || strcmp(at[i].elname, "H")) : 1) )
#define ALIASED_AT(i) (0 < NUM_ISO_H(at, i))
/***********************************************************************************/
#if ( TEST_RENUMB_ATOMS_SAVE_LONGEST == 1 || TEST_RENUMB_SWITCH == 1 )
int WriteToSDfile( const INP_ATOM_DATA *inp_at_data, INCHI_IOSTREAM* fcb, const char* name, const char* comment,
                   const char *szLabel, const char *szValue)
{
    int i, j, k, num_bonds=0, ret=0, bAtomsDT = 1 /* treat D, T as normal atoms */, bV2000 = 0 /*V2000 Molfile */;
    int bAtomNeedsAlias;
    int flag_bad_charge=0, flag_bad_iso=0, nNumAddLines=0, nNumIsoLines=0, nNumChargeLines=0, nNumRadicalLines=0, nNumAliasLines=0;
    int nNumNecessaryIsoLines = 0, nNumNecessaryChgLines = 0, nNumNecessaryRadLines = 0;
    /*sp_ATOM *at; */
    /*float fzero=0.0F;*/
    double x, y, z;
    int bNext /*, s*/;
    const inp_ATOM *at = inp_at_data->at_fixed_bonds? inp_at_data->at_fixed_bonds : inp_at_data->at;
    int num_atoms      = inp_at_data->num_at;
    /*at = species->atom;*/
    
    /*inchi_ios_eprint(fcb,"%ld.MOL\n",species->casno);*/
    {   /* block start */
        char strLocName[82];
        memset(strLocName, 0, sizeof(strLocName) );
        if ( name && *name ) {
            strncpy( strLocName, name, 80 );
        }
        inchi_ios_print_nodisplay( fcb,"%s\n", strLocName );
    }   /* block end */
    /**********************************************************************/
    /**                                                                  **/
    /** Important: Atoms with alias cannot have charge, radical, or      **/
    /**            isotope differences.                                  **/
    /**                                                                  **/
    /**            Atoms with alias cannot be abnormal.                  **/
    /**                                                                  **/
    /** Abnormal atoms are atoms which need M  CHG, M RAD, M  ISO        **/
    /**                                                                  **/
    /**********************************************************************/

    /*                                    F10.5     F12.5       I6
                     IIPPPPPPPPMMDDYYHHmmddSSssssssssssEEEEEEEEEEEERRRRRR
    inchi_ios_eprint( fcb,"NISTTRANHP09089809272D 1   1.0         0.0    %6ld\n", lEpa);*/

    /*^^^
    inchi_ios_print_nodisplay( fcb,"  -%s v%s SDfile Output                         \n", INCHI_NAME, INCHI_VERSION);

    Changed 01/10/2009 to conform CTFile specification (by Symyx request)*/
    inchi_ios_print_nodisplay( fcb,
    /*   IIPPPPPPPPMMDDYYHHmmddSSssssssssssEEEEEEEEEEEERRRRRR*/
        "  InChIV10                                     \n");



    /*y_fprintf(fcb, "  -CPSS-  1213981200n\n");*/

    {   /*block start*/
        char strLocName[82];
        
        memset(strLocName, 0, sizeof(strLocName) );
        if ( comment && *comment ) {
            strncpy( strLocName, comment, 80 );
        }
        inchi_ios_print_nodisplay( fcb,"%s\n", strLocName );
    }   /*block end*/
    for (i=0; i< num_atoms; i++)
        num_bonds += at[i].valence;
    num_bonds /= 2;
    
    /*find if we need "M  CHG", "M  RAD", "M  ISO" */
    for (i=0, nNumAddLines = 0; i < num_atoms; i++) {
        if ( bAtomNeedsAlias = ALIASED_AT(i) ) {
            nNumAliasLines  += 2 * bAtomNeedsAlias;
        } else {
            nNumNecessaryIsoLines += ABNORMAL_ISO(i);
            nNumNecessaryChgLines += ABNORMAL_CHG(i);
            nNumNecessaryRadLines += ABNORMAL_RAD(i);
            nNumIsoLines          += ANY_ISO(i, bAtomsDT); 
            nNumChargeLines       += ANY_CHG(i);
            nNumRadicalLines      += ANY_RAD(i);
        }
    }
    if ( !bV2000 ) {
        if ( !nNumNecessaryRadLines && !nNumNecessaryChgLines ) {
            nNumRadicalLines = 0;
            nNumChargeLines  = 0;
        }
        if ( !nNumNecessaryIsoLines ) {
            nNumIsoLines = 0;
        }
    }


    /* count additional M lines*/
    nNumChargeLines  = ( nNumChargeLines  + 7 ) / 8;
    nNumRadicalLines = ( nNumRadicalLines + 7 ) / 8;
    nNumIsoLines     = ( nNumIsoLines     + 7 ) / 8;
    
    nNumAddLines = nNumChargeLines + nNumRadicalLines + nNumIsoLines + nNumAliasLines; /* 1 for M  END*/
    
    if ( nNumAddLines || bV2000 ) {
        nNumAddLines += 1; /* add 1 for "M  END" line*/
    }
    
    /*                         aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv*/
    inchi_ios_print_nodisplay(fcb,"%3d%3d  0  0  0  0  0  0  0  0%3d%s\n",num_atoms, num_bonds, nNumAddLines,nNumAddLines?" V2000":"");
    /* atoms block*/
    for (i=0; i < num_atoms; i++)  {
        char elname[ATOM_EL_LEN];
        int  iso       = 0;
        int  charge    = 0;
        int  valence   = 0;
        int  nIsotopeH = IS_DEUTERIUM(i)? 1 : IS_TRITIUM(i)? 2 : 0;
        bAtomNeedsAlias = ALIASED_AT(i);   /* Has implicit D and/or T neighbors */
        memset( elname, 0, sizeof(elname) );
        
        if ( bAtomNeedsAlias ) {
            /* alias */
            strcpy ( elname, "C" );
        } else {
            /* isotope*/
            if ( nIsotopeH ) {
                strcpy( elname, bAtomsDT? ( nIsotopeH==1? "D" : "T" ) : "H" );
            } else {
                strncpy ( elname, at[i].elname, sizeof(elname)-1 );
            }
            if ( !ABNORMAL_CHG(i) && !ANY_RAD(i) ) {
                /* charge*/
                /* Only atoms without alias can be here*/
                switch ( at[i].charge ) {
                   case  3: charge = 1; break;
                   case  2: charge = 2; break;
                   case  1: charge = 3; break;
                   case -1: charge = 5; break;
                   case -2: charge = 6; break;
                   case -3: charge = 7; break;
                   case  0: charge = 0; break;
                   default: flag_bad_charge = 1; break;
                };
            }
            /* radical*/
            if ( ANY_RAD(i) && !ANY_CHG(i) ) {
                if ( at[i].radical == RADICAL_DOUBLET ) {
                    charge = 4;
                }
            }
        }
        /* allow isotopic shift for aliased atoms */
        if ( NORMAL_ISO(i, bAtomsDT) ) {
            iso = at[i].iso_atw_diff > 0? at[i].iso_atw_diff-1:
                  at[i].iso_atw_diff < 0? at[i].iso_atw_diff  :
                  nIsotopeH? nIsotopeH : (flag_bad_iso ++, 0);
        }
        
        x = at[i].x;
        y = at[i].y;
        z = at[i].z;

        if( at[i].num_H > 0 ) {
            for ( j = 0, valence = 0; j < at[i].valence; j++ ) {
                switch( k = at[i].bond_type[j] ) { /* fixed valence calculation 12-23-99 DCh.*/
                case 2:
                case 3:
                    valence += 2*k;
                    break;
                case 4:
                    valence += 3;
                    break;
                default:
                    valence += 2;
                }
            }
            valence = valence/2 + at[i].num_H;
        } else
        /* Added 07-09-2003 DCh*/
        if ( at[i].chem_bonds_valence > 0 ) {
            valence = at[i].chem_bonds_valence;
        } else
        /* Added 07-09-2003 DCh*/
        if ( !at[i].valence && !at[i].num_H && !at[i].chem_bonds_valence ) {
            valence = 15;
        }
        /*inchi_ios_eprint(fcb,"%10.4f%10.4f%10.4f %-3.3s%2d%3d  0     0  0  0  0  0  0\n",*/
        /*    (float)at[i].x, (float)(-at[i].y), fzero, at[i].elname, iso, charge);*/
        /*              xxxxxxyyyyyyzzzzzz aaa____ddcccsssnnnbbbvvvrrriiimmmeee  */
        inchi_ios_print_nodisplay(fcb,"%10.4f%10.4f%10.4f %-3.3s%2d%3d  0     0%3d  0  0  0  0\n",
                   x, y, z, elname, (int)iso, (int)charge, valence /* at[i].special*/);
            /* reflect image against x-axis;
               when transforming MOLfile back to STDATA in mol_to_stdata(...),
               make one more reflection to restore original orientation.
               Reason: in MS Search y-axis is directed from top to bottom,
                       while in MOLfile y-axis goes from bottom to top.
             */
    }        
    bNext = 0; /* debug only*/
    
    /* bonds*/
    for (i=0; i< num_atoms; i++) {
        for (j=0; j<at[i].valence; j++) {
            if (i < at[i].neighbor[j]) {
                if ( k=at[i].bond_stereo[j] ) {
                    /* bond stereo */
                    if ( k < 0 ) {
                        /* transposition */
                        inchi_ios_print_nodisplay(fcb,"%3u%3u%3u%3u  0  0  0\n",
                            (unsigned)(at[i].neighbor[j]+1), (unsigned)(i+1), (unsigned)(at[i].bond_type[j]), (unsigned)abs(k));
                    } else {
                        /* no transposition*/
                        inchi_ios_print_nodisplay(fcb,"%3u%3u%3u%3u  0  0  0\n",
                            (unsigned)(i+1), (unsigned)(at[i].neighbor[j]+1), (unsigned)(at[i].bond_type[j]), (unsigned)abs(k));
                    }
                } else {
                    inchi_ios_print_nodisplay(fcb,"%3u%3u%3u  0  0  0  0\n",
                        (unsigned)(i+1), (unsigned)(at[i].neighbor[j]+1), (unsigned)(at[i].bond_type[j]));
                }
            }
        }
    }
    if ( nNumAddLines ) {
        char str_m[66], entry[10];
        int  num_m;
        
        /* Aliases. 5-3-99 DCh.*/
        if ( nNumAliasLines ) {
            num_m = 0;
            for (i=0; i < num_atoms; i++) {
                if ( ALIASED_AT(i) ) {
                    int num_H;
                    inchi_ios_print_nodisplay( fcb, "A  %d\n", i+1 );
                    num_m ++;
                    strcpy( str_m, at[i].elname );
                    /* Add H, D, T */
                    if ( num_H = at[i].num_H + at[i].num_iso_H[0] ) { /* protium is lost here */
                        strcat( str_m, "H" );
                        if ( num_H > 1 ) {
                            sprintf( str_m + strlen(str_m), "%d", num_H );
                        }
                    }
                    if ( num_H = at[i].num_iso_H[1] ) { /* deuterium */
                        strcat( str_m, "D" );
                        if ( num_H > 1 ) {
                            sprintf( str_m + strlen(str_m), "%d", num_H );
                        }
                    }
                    if ( num_H = at[i].num_iso_H[2] ) { /* Tritium */
                        strcat( str_m, "T" );
                        if ( num_H > 1 ) {
                            sprintf( str_m + strlen(str_m), "%d", num_H );
                        }
                    }
                    /* Add charge to the Alias */
                    if ( at[i].charge){
                        strcat(str_m, at[i].charge>0? "+" : "-");
                        if ( 1 < (j=abs(at[i].charge)) )
                            sprintf( str_m+strlen(str_m), "%d", j );
                    }
                    /* Add radical to the Alias */
                    switch( at[i].radical ) {
                    case RADICAL_SINGLET:
                        strcat( str_m, ":" );
                        break;
                    case RADICAL_DOUBLET:
                        strcat( str_m, "^" );
                        break;
                    case RADICAL_TRIPLET:
                        strcat( str_m, "^^" );
                        break;
                    }
                    inchi_ios_print_nodisplay( fcb, "%s\n", str_m );
                    num_m ++;
                }
            }
            if ( num_m != nNumAliasLines ) {
                /* error in lines counting*/
                ret ++;
            }
        }
        /* charges*/
        str_m[0] = 0;
        num_m    = 0;
        if ( nNumChargeLines ) {
            for (i=0; i < num_atoms; i++) {
                if ( ANY_CHG(i) && !ALIASED_AT(i) ) {
                    sprintf( entry, " %3d %3d", i+1, (int)at[i].charge );
                    strcat( str_m, entry );
                    num_m ++;
                }
                if ( i == num_atoms-1 && num_m || num_m == 8 ) {
                    inchi_ios_print_nodisplay( fcb, "M  CHG%3d%s\n", num_m, str_m );
                    str_m[0] = 0;
                    num_m    = 0;
                }
            }
        }
        /* radicals*/
        str_m[0] = 0;
        num_m    = 0;
        if ( nNumRadicalLines ) {
            for (i=0; i < num_atoms; i++) {
                if ( ANY_RAD(i) && !ALIASED_AT(i) ) {
                    int radical = (at[i].radical==RADICAL_SINGLET ||
                                   at[i].radical==RADICAL_DOUBLET ||
                                   at[i].radical==RADICAL_TRIPLET)? at[i].radical : 0;
                    if ( radical ) {
                        sprintf( entry, " %3d %3d", i+1, radical );
                        strcat( str_m, entry );
                        num_m ++;
                    }
                }
                if ( i == num_atoms-1 && num_m || num_m == 8 ) {
                    inchi_ios_print_nodisplay( fcb, "M  RAD%3d%s\n", num_m, str_m );
                    str_m[0] = 0;
                    num_m    = 0;
                }
            }
        }
        /* isotopes*/
        str_m[0] = 0;
        num_m    = 0;
        if ( nNumIsoLines ) {
            int el_num, iso;
            for (i=0; i < num_atoms; i++) {
                if ( ANY_ISO(i,bAtomsDT) && !ALIASED_AT(i) ) {
                    if ( IS_DEUTERIUM(i) ) {
                        iso = 1;
                        el_num = 1;
                    } else
                    if ( IS_TRITIUM(i) ) {
                        iso = 2;
                        el_num = 1;
                    } else {
                        iso = at[i].iso_atw_diff > 0? at[i].iso_atw_diff-1 : at[i].iso_atw_diff;
                        el_num = at[i].el_number;
                    }
                    iso += get_atw_from_elnum( el_num );

                    sprintf( entry, " %3d %3d", i+1, iso );
                    strcat( str_m, entry );
                    num_m ++;
                }

                if ( i == num_atoms-1 && num_m || num_m == 8 ) {
                    inchi_ios_print_nodisplay( fcb, "M  ISO%3d%s\n", num_m, str_m );
                    str_m[0] = 0;
                    num_m    = 0;
                }
            }
        }
        inchi_ios_print_nodisplay( fcb, "M  END\n" );
    }
    if ( szValue && szValue[0] ) {
        if ( szLabel && szLabel[0] ) {
            inchi_ios_print_nodisplay( fcb, ">  <%s>\n", szLabel );
        } else {
            inchi_ios_print_nodisplay( fcb, ">  <ID>\n" );
        }
        inchi_ios_print_nodisplay( fcb, "%s\n\n", szValue );
    }
    inchi_ios_print_nodisplay(fcb, "$$$$\n");
    
    
    return ret;
    
}
#endif
/***************************************************************************************************/
int WriteOrigAtomDataToSDfile(const ORIG_ATOM_DATA *inp_at_data, INCHI_IOSTREAM * fcb, const char* name, 
                              const char* comment, int bChiralFlag, int bAtomsDT, const char *szLabel, const char *szValue)
{
    int i, j, k, num_bonds=0, ret=0;
    int bAtomNeedsAlias;
    int flag_bad_charge=0, flag_bad_iso = 0;
    int nNumAddLines=0, nNumIsoLines=0, nNumChargeLines=0, nNumRadicalLines=0, nNumAliasLines=0;
    int nNumNecessaryIsoLines = 0, nNumNecessaryChgLines = 0, nNumNecessaryRadLines = 0;
    int bV2000 = SDF_OUTPUT_V2000;
    double x, y, z;
    int bNext /*, s*/;
    const inp_ATOM *at = inp_at_data->at;
    int num_atoms      = inp_at_data->num_inp_atoms;
    /*at = species->atom;*/
    
    /*inchi_ios_eprint(fcb,"%ld.MOL\n",species->casno);*/
    {   /* block start */
        char strLocName[82];
        memset(strLocName, 0, sizeof(strLocName) );
        if ( name && *name ) {
            strncpy( strLocName, name, 80 );
             /* --- debug only ---
            if ( strstr( name, "#3959" ) ) {
                int stop = 1;
            }
            */
        }
        inchi_ios_print_nodisplay( fcb,"%s\n", strLocName );
    }   /* block end */
    /**********************************************************************/
    /**                                                                  **/
    /** Important: Atoms with alias cannot have charge, radical          **/
    /**            isotope differences are allowed                       **/
    /**                                                                  **/
    /**            Atoms with alias cannot be abnormal.                  **/
    /**                                                                  **/
    /** Abnormal atoms are atoms which need M  CHG, M RAD, M  ISO        **/
    /**                                                                  **/
    /** Output aliased atoms if they have implicit D or T                **/
    /**                                                                  **/
    /**********************************************************************/

    /*                                    F10.5     F12.5       I6
                     IIPPPPPPPPMMDDYYHHmmddSSssssssssssEEEEEEEEEEEERRRRRR
    inchi_ios_eprint( fcb,"NISTTRANHP09089809272D 1   1.0         0.0    %6ld\n", lEpa);*/
    /*^^^
    inchi_ios_print_nodisplay( fcb,"  %s v%s SDfile Output                       \n", INCHI_NAME, INCHI_VERSION);
    
    Changed 01/10/2009 to conform CTFile specification (by Symyx request)*/
    inchi_ios_print_nodisplay( fcb,
    /*   IIPPPPPPPPMMDDYYHHmmddSSssssssssssEEEEEEEEEEEERRRRRR*/
        "  InChIV10                                     \n");
    /*y_fprintf(fcb, "  -CPSS-  1213981200n\n");*/

    {   /*block start*/
        char strLocName[82];
        
        memset(strLocName, 0, sizeof(strLocName) );
        if ( comment && *comment ) {
            strncpy( strLocName, comment, 80 );
        }
        inchi_ios_print_nodisplay( fcb,"%s\n", strLocName );
    }   /*block end*/
    for (i=0; i< num_atoms; i++)
        num_bonds += at[i].valence;
    num_bonds /= 2;
    
    /*find if we need "M  CHG" and "M  RAD"*/
    for (i=0; i < num_atoms; i++) {
        if ( (bAtomNeedsAlias = ALIASED_AT(i)) ) {     /* has isotopic implicit D or T; ignoring pure 1H */
            nNumAliasLines  += 2 * bAtomNeedsAlias;
        } else {
        /* abnormal means atom needs CHG, RAD, or ISO entry */
        /* nNumAddLines    += ABNORMAL_AT(i); */ 
        /* nNumIso         += ( 0 == strcmp( at[i].elname, "D" ) || ( 0 == strcmp( at[i].elname, "T" ) || at[i].iso_atw_diff ) ); */
        /* nNumAddIso      += at[i].iso_atw_diff && (at[i].iso_atw_diff == 1 || at[i].iso_atw_diff < -3 || at[i].iso_atw_diff > 5 ); */
            nNumNecessaryIsoLines += ABNORMAL_ISO(i);
            nNumNecessaryChgLines += ABNORMAL_CHG(i);
            nNumNecessaryRadLines += ABNORMAL_RAD(i);
            nNumIsoLines          += ANY_ISO(i, bAtomsDT); 
            nNumChargeLines       += ANY_CHG(i);
            nNumRadicalLines      += ANY_RAD(i);
        }
    }
    nNumChargeLines  = ( nNumChargeLines  + 7 ) / 8;
    nNumRadicalLines = ( nNumRadicalLines + 7 ) / 8;
    nNumIsoLines     = ( nNumIsoLines     + 7 ) / 8;

    if ( !bV2000 ) {
        if ( !nNumNecessaryRadLines && !nNumNecessaryChgLines ) {
            nNumRadicalLines = 0;
            nNumChargeLines  = 0;
        }
        if ( !nNumNecessaryIsoLines ) {
            nNumIsoLines = 0;
        }
    }


    /* recalculate number of added lines */
    nNumAddLines = nNumChargeLines + nNumRadicalLines + nNumIsoLines + nNumAliasLines; /* 1 for M  END*/
    
    if ( nNumAddLines || bV2000 ) {
        nNumAddLines += 1; /* add 1 for "M  END" line*/
    }
    
    /*                         aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv*/
    inchi_ios_print_nodisplay(fcb,"%3d%3d  0  0%3d  0  0  0  0  0%3d%s\n",
                          num_atoms, num_bonds, bChiralFlag?1:0, nNumAddLines,nNumAddLines?" V2000":"");
    /* atoms block*/
    for (i=0; i < num_atoms; i++)  {
        char elname[ATOM_EL_LEN] = "\0\0\0\0\0";
        int  iso       = 0;
        int  charge    = 0;
        int  valence   = 0;
        int  nIsotopeH = IS_DEUTERIUM(i)? 1 : IS_TRITIUM(i)? 2 : 0;
        int  bonds_val;
        bAtomNeedsAlias = ALIASED_AT(i);
        memset( elname, 0, sizeof(elname) );
        
        if ( bAtomNeedsAlias ) {
            /* alias */
            strcpy ( elname, "C" );
        } else {
            /* isotope*/
            if ( nIsotopeH ) {
                strcpy( elname, bAtomsDT? ( nIsotopeH==1? "D" : "T" ) : "H" );
            } else {
                strncpy ( elname, at[i].elname, sizeof(elname)-1 );
            }
            if ( !ABNORMAL_CHG(i) && !ANY_RAD(i) ) {
                /* charge*/
                /* Only atoms without alias can be here*/
                switch ( at[i].charge ) {
                   case  3: charge = 1; break;
                   case  2: charge = 2; break;
                   case  1: charge = 3; break;
                   case -1: charge = 5; break;
                   case -2: charge = 6; break;
                   case -3: charge = 7; break;
                   case  0: charge = 0; break;
                   default: flag_bad_charge = 1; break;
                };
            }
            /* radical*/
            if ( ANY_RAD(i) && !ANY_CHG(i) ) {
                if ( at[i].radical == RADICAL_DOUBLET ) {
                    charge = 4;
                }
            }
        }
        /* allow isotopic shift for aliased atoms */
        if ( NORMAL_ISO(i, bAtomsDT) ) {
            iso = at[i].iso_atw_diff > 0? at[i].iso_atw_diff-1:
                  at[i].iso_atw_diff < 0? at[i].iso_atw_diff  :
                  nIsotopeH? nIsotopeH : (flag_bad_iso ++, 0);
        }
        
        x = at[i].x;
        y = at[i].y;
        z = at[i].z;

        /* valence -- set only if needed */
        bonds_val = nBondsValenceInpAt( at+i, NULL, NULL );
        valence=needed_unusual_el_valence( at[i].el_number, at[i].charge, at[i].radical,
                                 at[i].chem_bonds_valence, bonds_val, NUMH(at, i), at[i].valence );
        if ( valence < 0 ) {
            valence = 15;  /* means no bonds nor H */
        }

        /*inchi_ios_eprint(fcb,"%10.4f%10.4f%10.4f %-3.3s%2d%3d  0     0  0  0  0  0  0\n",*/
        /*    (float)at[i].x, (float)(-at[i].y), fzero, at[i].elname, iso, charge);*/
        /*              xxxxxxyyyyyyzzzzzz aaa____ddcccsssnnnbbbvvvrrriiimmmeee  */
        inchi_ios_print_nodisplay(fcb,"%10.4f%10.4f%10.4f %-3.3s%2d%3d  0     0%3d  0  0  0  0\n",
                   x, y, z, elname, (int)iso, (int)charge, valence /* at[i].special*/);
            /* reflect image against x-axis;
               when transforming MOLfile back to STDATA in mol_to_stdata(...),
               make one more reflection to restore original orientation.
               Reason: in MS Search y-axis is directed from top to bottom,
                       while in MOLfile y-axis goes from bottom to top.
             */
    }        
    bNext = 0; /* debug only*/
    
    /* bonds*/
    for (i=0; i< num_atoms; i++) {
        for (j=0; j<at[i].valence; j++) {
            if (i < at[i].neighbor[j]) {
                if ( (k=at[i].bond_stereo[j]) ) {
                    /* bond stereo */
                    if ( k < 0 ) {
                        /* transposition */
                        inchi_ios_print_nodisplay(fcb,"%3u%3u%3u%3u  0  0  0\n",
                            (unsigned)(at[i].neighbor[j]+1), (unsigned)(i+1), (unsigned)(at[i].bond_type[j]), (unsigned)abs(k));
                    } else {
                        /* no transposition*/
                        inchi_ios_print_nodisplay(fcb,"%3u%3u%3u%3u  0  0  0\n",
                            (unsigned)(i+1), (unsigned)(at[i].neighbor[j]+1), (unsigned)(at[i].bond_type[j]), (unsigned)abs(k));
                    }
                } else {
                    inchi_ios_print_nodisplay(fcb,"%3u%3u%3u  0  0  0  0\n",
                        (unsigned)(i+1), (unsigned)(at[i].neighbor[j]+1), (unsigned)(at[i].bond_type[j]));
                }
            }
        }
    }
    if ( nNumAddLines ) {
        char str_m[66], entry[10];
        int  num_m;
        
        /* Aliases. 5-3-99 DCh.*/
        if ( nNumAliasLines ) {
            num_m = 0;
            for (i=0; i < num_atoms; i++) {
                if ( ALIASED_AT(i) ) {
                    int len;
                    inchi_ios_print_nodisplay( fcb, "A  %d\n", i+1 );
                    num_m ++;
                    len = sprintf( str_m, "%s", at[i].elname );
                    /* add isotopic H to the alias */
                    for ( k = 0; k < NUM_H_ISOTOPES; k ++ ) {
                        int num_H = at[i].num_iso_H[k] + (k? 0:at[i].num_H);
                        if ( num_H ) {
                            len += sprintf( str_m+len, "%s", k == 0? "H" : k==1? "D" : k==2? "T" : "?" );
                            if ( num_H != 1 ) {
                                len += sprintf( str_m+len, "%d", num_H );
                            }
                        }
                    }
                    /* Add charge to the Alias */
                    if ( at[i].charge){
                        len += sprintf(str_m+len, "%s", at[i].charge>0? "+" : "-");
                        if ( 1 < (j=abs(at[i].charge)) ) {
                            len += sprintf( str_m+len, "%d", j );
                        }
                    }
                    /* Add radical to the Alias */
                    if ( at[i].radical == RADICAL_SINGLET ) {
                        len += sprintf( str_m+len, "%s", ":" );
                    } else
                    if ( at[i].radical == RADICAL_DOUBLET ) {
                        len += sprintf( str_m+len, "%s", "^" );
                    } else
                    if ( at[i].radical == RADICAL_TRIPLET ) {
                        len += sprintf( str_m+len, "%s", "^^" );
                    }
                    inchi_ios_print_nodisplay( fcb, "%s\n", str_m );
                    num_m ++;
                }
            }
            if ( num_m != nNumAliasLines ) {
                /* error in lines counting*/
                ret ++;
            }
        }
        /* charges*/
        str_m[0] = 0;
        num_m    = 0;
        if ( nNumChargeLines ) {
            for (i=0; i < num_atoms; i++) {
                if ( at[i].charge && !ALIASED_AT(i) ) {
                    sprintf( entry, " %3d %3d", i+1, (int)at[i].charge );
                    strcat( str_m, entry );
                    num_m ++;
                }
                if ( (i == num_atoms-1 && num_m) || num_m == 8 ) {
                    inchi_ios_print_nodisplay( fcb, "M  CHG%3d%s\n", num_m, str_m );
                    str_m[0] = 0;
                    num_m    = 0;
                }
            }
        }
        /* radicals*/
        str_m[0] = 0;
        num_m    = 0;
        if ( nNumRadicalLines ) {
            for (i=0; i < num_atoms; i++) {
                if ( at[i].radical && !ALIASED_AT(i) ) {
                    int radical = (at[i].radical==RADICAL_SINGLET ||
                                   at[i].radical==RADICAL_DOUBLET ||
                                   at[i].radical==RADICAL_TRIPLET)? at[i].radical : 0;
                    if ( radical ) {
                        sprintf( entry, " %3d %3d", i+1, radical );
                        strcat( str_m, entry );
                        num_m ++;
                    }
                }
                if ( (i == num_atoms-1 && num_m) || num_m == 8 ) {
                    inchi_ios_print_nodisplay( fcb, "M  RAD%3d%s\n", num_m, str_m );
                    str_m[0] = 0;
                    num_m    = 0;
                }
            }
        }
        /* isotopes*/
        str_m[0] = 0;
        num_m    = 0;
        if ( nNumIsoLines ) {
            int el_num, iso;
            for (i=0; i < num_atoms; i++) {
                /*
                if ( 0 == strcmp( at[i].elname, "D" ) ) {
                    sprintf( entry, " %3d %3d", i+1, 2 );
                    strcat( str_m, entry );
                    num_m ++;
                } else
                if ( 0 == strcmp( at[i].elname, "T" ) ) {
                    sprintf( entry, " %3d %3d", i+1, 3 );
                    strcat( str_m, entry );
                    num_m ++;
                } else
                if ( k = at[i].iso_atw_diff ) {
                    int mw = get_atw_from_elnum( at[i].el_number );
                    mw += (k > 0)? k-1 : k;
                    sprintf( entry, " %3d %3d", i+1, mw );
                    strcat( str_m, entry );
                    num_m ++;
                }
                */
                if ( ANY_ISO(i, bAtomsDT) && !ALIASED_AT(i) ) {
                    if ( IS_DEUTERIUM(i) ) {
                        iso = 1;
                        el_num = 1;
                    } else
                    if ( IS_TRITIUM(i) ) {
                        iso = 2;
                        el_num = 1;
                    } else {
                        iso = at[i].iso_atw_diff > 0? at[i].iso_atw_diff-1 : at[i].iso_atw_diff;
                        el_num = at[i].el_number;
                    }
                    iso += get_atw_from_elnum( el_num );

                    sprintf( entry, " %3d %3d", i+1, iso );
                    strcat( str_m, entry );
                    num_m ++;
                }


                if ( (i == num_atoms-1 && num_m) || num_m == 8 ) {
                    inchi_ios_print_nodisplay( fcb, "M  ISO%3d%s\n", num_m, str_m );
                    str_m[0] = 0;
                    num_m    = 0;
                }
            }
        }
        inchi_ios_print_nodisplay( fcb, "M  END\n" );
    }
    if ( szValue && szValue[0] ) {
        if ( szLabel && szLabel[0] ) {
            inchi_ios_print_nodisplay( fcb, "> <%s>\n", szLabel );
        } else {
            inchi_ios_print_nodisplay( fcb, "> <ID>\n" );
        }
        inchi_ios_print_nodisplay( fcb, " %s\n\n", szValue );
    }
    inchi_ios_print_nodisplay(fcb, "$$$$\n");
    
    
    return ret;
    
}
#if ( FIX_ADJ_RAD == 1 )
/*************************************************************************/
int FixNextRadicals( int cur_at, inp_ATOM *at );
int FixNextRadicals( int cur_at, inp_ATOM *at )
{
    int j, neigh, num_found = 0;
    for ( j = 0; j < at[cur_at].valence; j ++ ) {
        neigh = at[cur_at].neighbor[j];
        if ( at[neigh].radical == RADICAL_DOUBLET ) {
            at[neigh].radical = 0;
            num_found ++;
            num_found += FixNextRadicals( neigh, at );
        }
    }
    return num_found;
}
/*************************************************************************/
int FixAdjacentRadicals( int num_inp_atoms, inp_ATOM *at )
{
    int i, j;
    char *bVisited = NULL;
    int  nNumFound = 0, neigh, cur_found;
    for ( i = 0; i < num_inp_atoms; i ++ ) {
        if ( at[i].radical == RADICAL_DOUBLET ) {
            cur_found = 1;
            for ( j = 0; j < at[i].valence; j ++ ) {
                neigh = at[i].neighbor[j];
                if ( at[neigh].radical == RADICAL_DOUBLET ) {
                    cur_found ++;
                }
            }
            if ( cur_found >= 3 ) {
                nNumFound ++;
                at[i].radical = 0;
                nNumFound += FixNextRadicals( i, at );
            }
        }
    }
    return nNumFound;
}
#endif

#ifdef COMPILE_ANSI_ONLY
#ifndef TARGET_API_LIB
/*
#include <stdio.h>
#include "inpdef.h"
*/
void PrintFileName( const char *fmt, 
                   FILE *output_file, /* INCHI_IOSTREAM *output_file,  */
                   const char *szFname )
{
    inchi_print_nodisplay( output_file, fmt, szFname );
}
#endif
#endif

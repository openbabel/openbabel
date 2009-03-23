/*
 * International Chemical Identifier (InChI)
 * Version 1
 * Software version 1.02
 * October 31, 2008
 * Developed at NIST
 *
 * The InChI library and programs are free software developed under the
 * auspices of the International Union of Pure and Applied Chemistry (IUPAC);
 * you can redistribute this software and/or modify it under the terms of 
 * the GNU Lesser General Public License as published by the Free Software 
 * Foundation:
 * http://www.opensource.org/licenses/lgpl-license.php
 */


#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <ctype.h>

#include "mode.h"

#include "inpdef.h"
#include "util.h"
#include "readmol.h"

#include "ichicomp.h"
#include "ichierr.h"
#include "extr_ct.h"

#include "ichi_io.h"

#define NO_ATOM          (-1) /* non-existent (central) atom */

#define SB_PARITY_FLAG  0x38 /* disconnected structure has undef. parity */
#define SB_PARITY_SHFT  3
#define SB_PARITY_MASK  0x07
#define SB_PARITY_1(X) (X & SB_PARITY_MASK)  /* refers to connected structure */
#define SB_PARITY_2(X) (((X) >> SB_PARITY_SHFT) & SB_PARITY_MASK) /* refers to connected structure */

/* 0D parity types */
typedef enum tagINCHIStereoType0D {
   INCHI_StereoType_None        = 0,
   INCHI_StereoType_DoubleBond  = 1,
   INCHI_StereoType_Tetrahedral = 2,
   INCHI_StereoType_Allene      = 3
} inchi_StereoType0D;

/* 0D parities */
typedef enum tagINCHIStereoParity0D {
   INCHI_PARITY_NONE      = AB_PARITY_NONE,
   INCHI_PARITY_ODD       = AB_PARITY_ODD ,  /* 'o' */
   INCHI_PARITY_EVEN      = AB_PARITY_EVEN,  /* 'e' */
   INCHI_PARITY_UNKNOWN   = AB_PARITY_UNKN,  /* 'u' */
   INCHI_PARITY_UNDEFINED = AB_PARITY_UNDF   /* '?' -- should not be used; however, see Note above */
} inchi_StereoParity0D;

/* bond type definitions */
typedef enum tagINCHIBondType {
   INCHI_BOND_TYPE_NONE    =  0,
   INCHI_BOND_TYPE_SINGLE  =  BOND_TYPE_SINGLE,
   INCHI_BOND_TYPE_DOUBLE  =  BOND_TYPE_DOUBLE,
   INCHI_BOND_TYPE_TRIPLE  =  BOND_TYPE_TRIPLE,
   INCHI_BOND_TYPE_ALTERN  =  BOND_TYPE_ALTERN
} inchi_BondType;

/* 2D stereo definitions */
typedef enum tagINCHIBondStereo2D {
   /* stereocenter-related; positive: the sharp end points to this atom  */
   INCHI_BOND_STEREO_NONE           =  0,
   INCHI_BOND_STEREO_SINGLE_1UP     =  STEREO_SNGL_UP,
   INCHI_BOND_STEREO_SINGLE_1EITHER =  STEREO_SNGL_EITHER,
   INCHI_BOND_STEREO_SINGLE_1DOWN   =  STEREO_SNGL_DOWN,
   /* stereocenter-related; negative: the sharp end points to the opposite atom  */
   INCHI_BOND_STEREO_SINGLE_2UP     = -STEREO_SNGL_UP,
   INCHI_BOND_STEREO_SINGLE_2EITHER = -STEREO_SNGL_EITHER,
   INCHI_BOND_STEREO_SINGLE_2DOWN   = -STEREO_SNGL_DOWN,
   /* stereobond-related */
   INCHI_BOND_STEREO_DOUBLE_EITHER  =  STEREO_DBLE_EITHER
} inchi_BondStereo2D;

/****************************************************************************/
typedef struct tagINCHIStereo0D {
    S_SHORT  neighbor[4];    /* 4 atoms always */
    S_SHORT  central_atom;   /* central tetrahedral atom or a central */
                            /* atom of allene; otherwise NO_ATOM */
    S_CHAR  type;           /* inchi_StereoType0D */
    S_CHAR  parity;         /* inchi_StereoParity0D: may be a combination of two parities: */
                            /* ParityOfConnected | (ParityOfDisconnected << 3), see Note above */
}inchi_Stereo0D;
/****************************************************************************/

/* This contains executable code. Included in lReadAux.c, e_ReadINCH.c, ReadINCH.c,  */
#include "aux2atom.h"


int INChIToOrigAtom( INCHI_IOSTREAM *inp_molfile, ORIG_ATOM_DATA *orig_at_data, int bMergeAllInputStructures,
                       int bGetOrigCoord, int bDoNotAddH, INPUT_TYPE nInputType,
                       char *pSdfLabel, char *pSdfValue, long *lSdfId,
                       INCHI_MODE *pInpAtomFlags, int *err, char *pStrErr )
{
    /* inp_ATOM       *at = NULL; */
    int            num_dimensions_new;
    int            num_inp_bonds_new;
    int            num_inp_atoms_new;
    inp_ATOM      *at_new     = NULL;
    inp_ATOM      *at_old     = NULL;
    int            nNumAtoms  = 0;
    MOL_COORD     *szCoordNew = NULL;
    MOL_COORD     *szCoordOld = NULL;
    int            i, j;

    if ( pStrErr ) {
        pStrErr[0] = '\0';
    }

    /*FreeOrigAtData( orig_at_data );*/
    if ( lSdfId )
        *lSdfId = 0;
    do {
        
        at_old     = orig_at_data? orig_at_data->at      : NULL; /*  save pointer to the previous allocation */
        szCoordOld = orig_at_data? orig_at_data->szCoord : NULL;
        num_inp_atoms_new =
            INChIToInpAtom( inp_molfile, (bGetOrigCoord && orig_at_data)? &szCoordNew : NULL,
                          bDoNotAddH, nInputType, orig_at_data? &at_new:NULL,
                          MAX_ATOMS,
                          &num_dimensions_new, &num_inp_bonds_new,
                          pSdfLabel, pSdfValue, lSdfId, pInpAtomFlags, err, pStrErr );
        if ( num_inp_atoms_new <= 0 && !*err ) {
            MOLFILE_ERR_SET (*err, 0, "Empty structure");
            *err = 98;
        } else
        if ( orig_at_data && !num_inp_atoms_new && 10 < *err && *err < 20 && orig_at_data->num_inp_atoms > 0 && bMergeAllInputStructures ) {
            *err = 0; /* end of file */
            break;
        } else
        if ( num_inp_atoms_new > 0 && orig_at_data ) {
            /*  merge pOrigDataTmp + orig_at_data => pOrigDataTmp; */
            nNumAtoms = num_inp_atoms_new + orig_at_data->num_inp_atoms;
            if ( nNumAtoms >= MAX_ATOMS ) {
                MOLFILE_ERR_SET (*err, 0, "Too many atoms");
                *err = 70;
                orig_at_data->num_inp_atoms = -1;
            } else
            if ( !at_old ) {
                /* the first structure */
                orig_at_data->at      = at_new;
                orig_at_data->szCoord = szCoordNew; 
                at_new = NULL;
                szCoordNew = NULL;
                orig_at_data->num_inp_atoms  = num_inp_atoms_new;
                orig_at_data->num_inp_bonds  = num_inp_bonds_new;
                orig_at_data->num_dimensions = num_dimensions_new;
            } else
            if ( (orig_at_data->at = ( inp_ATOM* ) inchi_calloc( nNumAtoms, sizeof(inp_ATOM) )) &&
                  (!szCoordNew || (orig_at_data->szCoord = (MOL_COORD *) inchi_calloc( nNumAtoms, sizeof(MOL_COORD) ))) ) {
                /*  switch at_new <--> orig_at_data->at; */
                if ( orig_at_data->num_inp_atoms ) {
                    memcpy( orig_at_data->at, at_old, orig_at_data->num_inp_atoms * sizeof(orig_at_data->at[0]) );
                    /*  adjust numbering in the newly read structure */
                    for ( i = 0; i < num_inp_atoms_new; i ++ ) {
                        for ( j = 0; j < at_new[i].valence; j ++ ) {
                            at_new[i].neighbor[j] += orig_at_data->num_inp_atoms;
                        }
                        at_new[i].orig_at_number += orig_at_data->num_inp_atoms; /* 12-19-2003 */
                    }
                    if ( orig_at_data->szCoord && szCoordOld ) {
                        memcpy( orig_at_data->szCoord, szCoordOld, orig_at_data->num_inp_atoms * sizeof(MOL_COORD) );
                    }
                }
                if ( at_old ) {
                    inchi_free( at_old );
                    at_old = NULL;
                }
                if ( szCoordOld ) {
                    inchi_free( szCoordOld );
                    szCoordOld = NULL;
                }
                /*  copy newly read structure */
                memcpy( orig_at_data->at + orig_at_data->num_inp_atoms,
                        at_new,
                        num_inp_atoms_new * sizeof(orig_at_data->at[0]) );
                if ( orig_at_data->szCoord && szCoordNew ) {
                    memcpy( orig_at_data->szCoord + orig_at_data->num_inp_atoms,
                            szCoordNew,
                            num_inp_atoms_new * sizeof(MOL_COORD) );
                }
                /*  add other things */
                orig_at_data->num_inp_atoms += num_inp_atoms_new;
                orig_at_data->num_inp_bonds += num_inp_bonds_new;
                orig_at_data->num_dimensions = inchi_max(num_dimensions_new, orig_at_data->num_dimensions);
            } else {
                MOLFILE_ERR_SET (*err, 0, "Out of RAM");
                *err = -1;
            }
        } else
        if ( num_inp_atoms_new > 0 ) {
            nNumAtoms += num_inp_atoms_new;
        }
        if ( at_new ) {
            inchi_free( at_new );
            at_new = NULL;
        }

    } while ( !*err && bMergeAllInputStructures );
    /*
    if ( !*err ) {
        orig_at_data->num_components =
            MarkDisconnectedComponents( orig_at_data );
        if ( orig_at_data->num_components == 0 ) {
            MOLFILE_ERR_SET (*err, 0, "No components found");
            *err = 99;
        }
        if ( orig_at_data->num_components < 0 ) {
            MOLFILE_ERR_SET (*err, 0, "Too many components");
            *err = 99;
        }
    }
    */
    if ( szCoordNew ) {
        inchi_free( szCoordNew );
    }
    if ( at_new ) {
        inchi_free( at_new );
    }
    if ( !*err && orig_at_data ) { /* added testing (orig_at_data != NULL) */
        if ( ReconcileAllCmlBondParities( orig_at_data->at, orig_at_data->num_inp_atoms, 0 ) ) {
            MOLFILE_ERR_SET (*err, 0, "Cannot reconcile stereobond parities");  /*   <BRKPT> */
            if (!orig_at_data->num_dimensions) {
                *err = 1;
            }
        }
    }
    if ( *err ) {
        FreeOrigAtData( orig_at_data );
    }
    if ( *err && !(10 < *err && *err < 20) && pStrErr && !pStrErr[0] ) {
        MOLFILE_ERR_SET (*err, 0, "Unknown error");  /*   <BRKPT> */
    }
    return orig_at_data? orig_at_data->num_inp_atoms : nNumAtoms;
}



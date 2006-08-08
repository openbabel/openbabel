/*
 * International Union of Pure and Applied Chemistry (IUPAC)
 * International Chemical Identifier (InChI)
 * Version 1
 * Software version 1.01
 * July 21, 2006
 * Developed at NIST
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "mode.h"

#include "comdef.h"
#include "readmol.h"
#include "inpdef.h"
#include "util.h"

#include "ichicomp.h"

#if( ADD_CMLPP == 1 )
#include "debug.h"
#endif
#include "mol2atom.h"

#define MIN_STDATA_X_COORD           0.0
#define MAX_STDATA_X_COORD         256.0
#define MIN_STDATA_Y_COORD           0.0
#define MAX_STDATA_Y_COORD         256.0
#define MIN_STDATA_Z_COORD           0.0
#define MAX_STDATA_Z_COORD         256.0
#define MAX_STDATA_AVE_BOND_LENGTH  20.0
#define MIN_STDATA_AVE_BOND_LENGTH  10.0


/* local prototypes */
inp_ATOM* mol_to_atom( MOL_DATA* mol_data, int *num_atoms, int *num_bonds, inp_ATOM* at_inp, int bDoNotAddH, int *err, char *pStrErr );
int mol_to_atom_xyz( MOL_DATA* mol_data, int num_atoms, inp_ATOM* at, int *err, char *pStrErr );
long GetMolfileNumber( MOL_HEADER_BLOCK *pHdr );

/******************************************************************************************************/
void FreeInpAtom( inp_ATOM **at )
{
    if ( at && *at ) {
        inchi_free( *at );
        *at = NULL;
    }
}
/******************************************************************************************************/
inp_ATOM *CreateInpAtom( int num_atoms )
{
    /*
    void *p = inchi_calloc(num_atoms, sizeof(inp_ATOM) );
    if ( p == (void*)0x009143A8 ) {
        int stop = 1;
    }
    return (inp_ATOM* )p;
    */
   return (inp_ATOM* ) inchi_calloc(num_atoms, sizeof(inp_ATOM) );
}
/******************************************************************************************************/
void FreeInpAtomData( INP_ATOM_DATA *inp_at_data )
{
    if ( inp_at_data ) {
        FreeInpAtom( &inp_at_data->at );
        FreeInpAtom( &inp_at_data->at_fixed_bonds );
        memset( inp_at_data, 0, sizeof(*inp_at_data) );
    }
}
/******************************************************************************************************/
int CreateInpAtomData( INP_ATOM_DATA *inp_at_data, int num_atoms, int create_at_fixed_bonds )
{
    FreeInpAtomData( inp_at_data );
    if ( (inp_at_data->at = CreateInpAtom( num_atoms )) &&
         (!create_at_fixed_bonds || (inp_at_data->at_fixed_bonds = CreateInpAtom( num_atoms) ) ) ) {
        inp_at_data->num_at = num_atoms;
        return 1;
    }
    FreeInpAtomData( inp_at_data );
    return 0;
}
/******************************************************************************************************/
void FreeCompAtomData( COMP_ATOM_DATA *inp_at_data )
{
    FreeInpAtom( &inp_at_data->at );
    if ( inp_at_data->nOffsetAtAndH )
        inchi_free( inp_at_data->nOffsetAtAndH );
    memset( inp_at_data, 0, sizeof(*inp_at_data) );
}
#ifndef INCHI_ANSI_ONLY
/******************************************************************************************************/
int CreateCompAtomData( COMP_ATOM_DATA *inp_at_data, int num_atoms, int num_components, int bIntermediateTaut )
{
    FreeCompAtomData( inp_at_data );
    if ( (inp_at_data->at = CreateInpAtom( num_atoms )) &&
         (num_components <= 1 || bIntermediateTaut ||
            (inp_at_data->nOffsetAtAndH = (AT_NUMB*)inchi_calloc(sizeof(inp_at_data->nOffsetAtAndH[0]), 2*(num_components+1))))) {

        inp_at_data->num_at = num_atoms;
        inp_at_data->num_components = (num_components>1)? num_components : 0;
        return 1;
    }
    FreeCompAtomData( inp_at_data );
    return 0;
}

/******************************************************************************************************/
void FreeInfAtom( inf_ATOM **at )
{
    if ( at && *at ) {
        inchi_free( *at );
        *at = NULL;
    }
}
/******************************************************************************************************/
inf_ATOM *CreateInfAtom( int num_atoms )
{
    return (inf_ATOM* ) inchi_calloc(num_atoms, sizeof(inf_ATOM) );
}
/******************************************************************************************************/
void FreeInfoAtomData( INF_ATOM_DATA *inf_at_data )
{
    FreeInfAtom( &inf_at_data->at );
    if ( inf_at_data->pStereoFlags )
        inchi_free( inf_at_data->pStereoFlags );
    memset(inf_at_data, 0, sizeof(*inf_at_data));
}
/******************************************************************************************************/
int CreateInfoAtomData( INF_ATOM_DATA *inf_at_data, int num_atoms, int num_components )
{
    FreeInfoAtomData( inf_at_data );
    memset( inf_at_data, 0, sizeof(*inf_at_data) );
    if ( (inf_at_data->at = CreateInfAtom( num_atoms )) &&
         (num_components <= 1 ||
          (inf_at_data->pStereoFlags = (AT_NUMB *)inchi_calloc(num_components+1, sizeof(inf_at_data->pStereoFlags[0])))
         )
       ) {
        inf_at_data->num_at = num_atoms;
        inf_at_data->num_components = num_components;
        return 1;
    }
    FreeInfoAtomData( inf_at_data );
    return 0;
}
/******************************************************************************************************/
int AllocateInfoAtomData( INF_ATOM_DATA *inf_at_data, int num_atoms, int num_components )
{
    if ( inf_at_data->at = CreateInfAtom( num_atoms ) ) {
        if ( num_components > 1 &&
            !(inf_at_data->pStereoFlags = (AT_NUMB *)inchi_calloc(num_components+1, sizeof(inf_at_data->pStereoFlags[0]))) ) {
            FreeInfAtom( &inf_at_data->at );
            return 0;
        }
        return 1;
    }
    return 0;
}
/******************************************************************************************************/
int DuplicateInfoAtomData( INF_ATOM_DATA *inf_at_data_to, const INF_ATOM_DATA *inf_at_data_from)
{
    *inf_at_data_to = *inf_at_data_from;
    if ( AllocateInfoAtomData( inf_at_data_to, inf_at_data_from->num_at, inf_at_data_from->num_components ) ) {
        memcpy( inf_at_data_to->at, inf_at_data_from->at,
                inf_at_data_from->num_at * sizeof(inf_at_data_to->at[0]));
        if ( inf_at_data_to->pStereoFlags && inf_at_data_from->pStereoFlags ) {
            memcpy( inf_at_data_to->pStereoFlags, inf_at_data_from->pStereoFlags,
                    (inf_at_data_from->num_components+1)*sizeof(inf_at_data_to->pStereoFlags[0]));
        }
        return 1;
    }
    return 0;
}
#endif /* ifndef INCHI_ANSI_ONLY */


#if( TEST_RENUMB_ATOMS == 1 )  /*  { */

/******************************************************************************************************/
int CopyInpAtomData( INP_ATOM_DATA *dest_inp_at_data, INP_ATOM_DATA *src_inp_at_data )
{
    int ret = 1;
    if ( !dest_inp_at_data->at  || dest_inp_at_data->num_at != src_inp_at_data->num_at ) {
        ret = CreateInpAtomData( dest_inp_at_data, src_inp_at_data->num_at, (NULL != src_inp_at_data->at_fixed_bonds) );
    } else {
        inp_ATOM *at  = dest_inp_at_data->at;  /*  save ptr to already allocated memory */
        inp_ATOM *at2 = dest_inp_at_data->at_fixed_bonds;
        *dest_inp_at_data = *src_inp_at_data; /*  copy all other (scalar) data */
        dest_inp_at_data->at = at;            /*  restore ptr to already allocated memory */
        dest_inp_at_data->at_fixed_bonds = at2;
    }
    if ( ret ) {
        memcpy( dest_inp_at_data->at, src_inp_at_data->at,
                src_inp_at_data->num_at*sizeof(dest_inp_at_data->at[0]) );
        if ( dest_inp_at_data->at_fixed_bonds && src_inp_at_data->at_fixed_bonds ) {
            memcpy( dest_inp_at_data->at_fixed_bonds, src_inp_at_data->at_fixed_bonds,
                src_inp_at_data->num_at*sizeof(dest_inp_at_data->at_fixed_bonds[0]) );
        }
    }
    return ret;
}
/******************************************************************************************************/
void RenumbInpAtomData( INP_ATOM_DATA *dest_inp_at_data, INP_ATOM_DATA *src_inp_at_data, AT_RANK *new_ord )
{
    int j, n, m, val;
#if( TEST_RENUMB_NEIGH == 1 )
    int i, k;
#endif
    int       num_atoms = src_inp_at_data->num_at;
    inp_ATOM *dest_at   = dest_inp_at_data->at;
    for ( n = 0; n < num_atoms; n ++ ) {
        m = new_ord[n];
        dest_at[m] = src_inp_at_data->at[n];
        dest_at[m].orig_compt_at_numb = (AT_NUMB)(m+1);  /*  new ordering number within the component */
        val = dest_at[m].valence;
        for ( j = 0; j < val; j ++ ) {
            dest_at[m].neighbor[j] = new_ord[dest_at[m].neighbor[j]];
        }
#if( TEST_RENUMB_NEIGH == 1 )
        for ( i = 0; i < val; i ++ ) {
            j = i;
            k = j + (rand() * (val-j)) / (RAND_MAX+1);
            if ( k >= val || j == k ) {
                continue;
            }
            swap( (char*)&dest_at[m].neighbor[j],    (char*)&dest_at[m].neighbor[k],    sizeof(dest_at[0].neighbor[0]) );
            swap( (char*)&dest_at[m].bond_stereo[j], (char*)&dest_at[m].bond_stereo[k], sizeof(dest_at[0].bond_stereo[0]) );
            swap( (char*)&dest_at[m].bond_type[j],   (char*)&dest_at[m].bond_type[k],   sizeof(dest_at[0].bond_type[0]) );
            /* adjust stereo bond links */
            if ( dest_at[m].sb_parity[0] ) {
                int a;
                for ( a = 0; a < MAX_NUM_STEREO_BONDS && dest_at[m].sb_parity[a]; a ++ ) {
                    
                    if ( k == (int)dest_at[m].sb_ord[a] ) {
                        dest_at[m].sb_ord[a] = j;
                    } else
                    if ( j == (int)dest_at[m].sb_ord[a] ) {
                        dest_at[m].sb_ord[a] = k;
                    }

                    if ( k == (int)dest_at[m].sn_ord[a] ) {
                        dest_at[m].sn_ord[a] = j;
                    } else
                    if ( j == (int)dest_at[m].sn_ord[a] ) {
                        dest_at[m].sn_ord[a] = k;
                    }
                }
            }
        }
#endif
    }
}
/******************************************************************************************************/
void MakeNewOrd( int num_atoms, AT_RANK *new_ord )
{
    int i, j, k;
    for ( i = 0; i < num_atoms; i ++ ) {
        j = i;
        k = (rand() * (num_atoms-i)) / (RAND_MAX+1);
        if ( k >= num_atoms || j == k ) {
            continue;
        }
        swap( (char*)&new_ord[j], (char*)&new_ord[k], sizeof(new_ord[0]) );
    }
}
#endif /*  } TEST_RENUMB_ATOMS == 1  */

/******************************************************************************************************/
inp_ATOM* mol_to_atom( MOL_DATA* mol_data, int *num_atoms, int *num_bonds, inp_ATOM* at_inp,
                       int bDoNotAddH, int *err, char *pStrErr )
{
    inp_ATOM *at = NULL;
    /* char      *bond_stereo = NULL; */
    AT_NUMB  *p1, *p2;
    int       i, a1, a2, n1, n2, bonds, iso_atw_diff;
    char      cBondStereo, cBondType;
    static int el_number_H = 0;


    if ( !el_number_H ) {
        el_number_H = get_periodic_table_number( "H" ); /* one-time initialization */
    }

    *err = 0;
    *num_atoms = *num_bonds = 0;
    /* check if MOLfile contains atoms */
    if ( !mol_data || !mol_data->ctab.MolAtom ||
         0 < mol_data->ctab.nNumberOfBonds && !mol_data->ctab.MolBond ||
         0 >= (*num_atoms = mol_data->ctab.nNumberOfAtoms) ) {
        /* MOLFILE_ERR_SET (*err, 0, "Empty structure"); */
        goto exit_function; /*  no structure */
    }
    /* allocate memory if necessary */
    if ( at_inp ) {
        at = at_inp;
    } else
    if ( !(at = CreateInpAtom( *num_atoms ) ) ) {
        *err = -1;
        MOLFILE_ERR_FIN (*err, -1, exit_function, "Out of RAM");
    }

    /* copy atom info */
    for ( i = 0; i < *num_atoms; i ++ ) {
        mystrncpy( at[i].elname, mol_data->ctab.MolAtom[i].szAtomSymbol, sizeof(at->elname) );
        /* at[i].chem_bonds_valence = mol_data->ctab.MolAtom[i].cValence; */ /*  MOLfile valence; will change */
        at[i].orig_at_number     = (AT_NUMB)(i+1);
        at[i].iso_atw_diff       = mol_data->ctab.MolAtom[i].cMassDifference;
        at[i].charge             = mol_data->ctab.MolAtom[i].cCharge;
        at[i].radical            = mol_data->ctab.MolAtom[i].cRadical;
        /* see mol_to_atom_xyz()
        at[i].x                  = mol_data->ctab.MolAtom[i].fX;
        at[i].y                  = mol_data->ctab.MolAtom[i].fY;
        at[i].z                  = mol_data->ctab.MolAtom[i].fZ;
        */
        iso_atw_diff             = mol_data->ctab.MolAtom[i].cMassDifference;
        at[i].iso_atw_diff       = iso_atw_diff==ZERO_ATW_DIFF? 1:
                                   iso_atw_diff>  0?            iso_atw_diff+1:
                                                                iso_atw_diff;
#if( SINGLET_IS_TRIPLET == 1 )
        if ( at[i].radical == RADICAL_SINGLET ) {
            at[i].radical = RADICAL_TRIPLET;
        }
#endif
#if( bRELEASE_VERSION != 1 )
        if ( isdigit( at[i].elname[0] ) ) { /*  for testing */
            mystrncpy( at[i].elname, "C", sizeof(at->elname) );
        }
#endif
        if ( ERR_ELEM == (n1 = get_periodic_table_number( at[i].elname ) ) ) {
            /*  Case when elname contains more than 1 element: extract number of H if possible */
            at[i].num_H = extract_H_atoms( at[i].elname, at[i].num_iso_H );
            if ( !at[i].elname[0] && NUMH(at, i) ) {
                /* alias contains only H. Added 2004-07-21, fixed 2004-07-22
                 * move the heaviest isotope to the "central atom"
                 * Note: this must be consistent with H-H treatment in remove_terminal_HDT()
                 */
                strcpy( at[i].elname, "H" );
                if ( NUM_ISO_H(at,i) ) {
                    int j;
                    for ( j = NUM_H_ISOTOPES-1; 0 <= j; j -- ) {
                        if ( at[i].num_iso_H[j] ) {
                            at[i].num_iso_H[j] --;
                            at[i].iso_atw_diff = 1 + j;
                            break;
                        }
                    }
                } else {
                    at[i].num_H --;
                }
            }
            if ( ERR_ELEM == (n1 = get_periodic_table_number( at[i].elname ) ) ) {
                n1 = 0;
            }
        }

        at[i].el_number = (U_CHAR) n1;
        if ( !n1 ) {
            *err |= 64; /*  Unrecognized aromatic bond(s) replaced with single */
            MOLFILE_ERR_SET (*err, 0, "Unknown element(s):");
            MOLFILE_ERR_SET (*err, 0, at[i].elname);
        } else
        /* replace explicit D or T with isotopic H (added 2003-06-02) */
        if ( el_number_H == n1 && !at[i].iso_atw_diff ) {
            switch( at[i].elname[0] ) {
            case 'D':
                at[i].iso_atw_diff = 2;
                mystrncpy( at[i].elname, "H", sizeof(at->elname) );
                break;
            case 'T':
                at[i].iso_atw_diff = 3;
                mystrncpy( at[i].elname, "H", sizeof(at->elname) );
                break;
            }
        }
    }


    /*---------------- stereo information notes. ------------------------

      Currently:  1. stereo sign
      =========   --------------
      MOLfile     (atom number = MOLfile atom number - 1, no stdata as an intermediate)
         |        if mol_data->ctab.MolBond[i].nAtomNo1 < mol_data->ctab.MolBond[i].nAtomNo2
         v        then
      inp_ATOM        stereo > 0
                  else
                     stereo < 0

                  2. neighbor z-coordinate
                  ------------------------
                  neighbor z-coord > 0 for Up if sign(stdata_bond_no) = sign(at[i].neighbor[j]-i)

    --------------------------------------------------------------------*/

    /* copy bond info */
    for ( i = 0, bonds = 0; i < mol_data->ctab.nNumberOfBonds; i ++ ) {
        cBondStereo = mol_data->ctab.MolBond[i].cBondStereo;
        cBondType   = mol_data->ctab.MolBond[i].cBondType;
        a1 = mol_data->ctab.MolBond[i].nAtomNo1-1;
        a2 = mol_data->ctab.MolBond[i].nAtomNo2-1;

        if ( a1 < 0 || a1 >= *num_atoms ||
             a2 < 0 || a2 >= *num_atoms ||
             a1 == a2 ) {
            *err |= 1; /*  bond for impossible atom number(s); ignored */
            MOLFILE_ERR_SET (*err, 0, "Bond to nonexistent atom");
            continue;
        }
        /*  check for multiple bonds between same atoms */
        p1 = is_in_the_list( at[a1].neighbor, (AT_NUMB)a2, at[a1].valence );
        p2 = is_in_the_list( at[a2].neighbor, (AT_NUMB)a1, at[a2].valence );
        if ( (p1 || p2) && (p1 || at[a1].valence < MAXVAL) && (p2 || at[a2].valence < MAXVAL) ) {
            n1 = p1? (p1 - at[a1].neighbor) : at[a1].valence ++;
            n2 = p2? (p2 - at[a2].neighbor) : at[a2].valence ++;
            MOLFILE_ERR_SET (*err, 0, "Multiple bonds between two atoms");
            *err |= 2; /*  multiple bonds between atoms */
        } else
        if ( !p1 && !p2 && at[a1].valence < MAXVAL && at[a2].valence < MAXVAL ) {
            n1 = at[a1].valence ++;
            n2 = at[a2].valence ++;
            bonds ++;
        } else {
            char szMsg[64];
            *err |= 4; /*  too large number of bonds. Some bonds ignored. */
            sprintf( szMsg, "Atom '%s' has more than %d bonds",
                            at[a1].valence>= MAXVAL? at[a1].elname:at[a2].elname, MAXVAL );
            MOLFILE_ERR_SET (*err, 0, szMsg);
            continue;
        }
        if ( cBondType < MIN_INPUT_BOND_TYPE || cBondType > MAX_INPUT_BOND_TYPE ) {
            char szBondType[16];
            sprintf( szBondType, "%d", cBondType );
            cBondType = 1;
            MOLFILE_ERR_SET (*err, 0, "Unrecognized bond type:");
            MOLFILE_ERR_SET (*err, 0, szBondType);
            *err |= 8; /*  Unrecognized Bond type replaced with single bond */
        }
        /* bond type */
        at[a1].bond_type[n1] =
        at[a2].bond_type[n2] = cBondType;
        /* connection */
        at[a1].neighbor[n1] = (AT_NUMB)a2;
        at[a2].neighbor[n2] = (AT_NUMB)a1;
        /* stereo */
        if ( cBondStereo == INPUT_STEREO_DBLE_EITHER /* 3 */ ) {
            at[a1].bond_stereo[n1] =
            at[a2].bond_stereo[n2] = STEREO_DBLE_EITHER;
        } else
        if ( cBondStereo == INPUT_STEREO_SNGL_UP     ||  /* 1 */
             cBondStereo == INPUT_STEREO_SNGL_EITHER ||  /* 4 */
             cBondStereo == INPUT_STEREO_SNGL_DOWN       /* 6 */ ) {
            char cStereo;
            switch ( cBondStereo ) {
            case INPUT_STEREO_SNGL_UP:
                cStereo = STEREO_SNGL_UP;
                break;
            case INPUT_STEREO_SNGL_EITHER:
                cStereo = STEREO_SNGL_EITHER;
                break;
            case INPUT_STEREO_SNGL_DOWN:
                cStereo = STEREO_SNGL_DOWN;
                break;
            }
            at[a1].bond_stereo[n1] =  cStereo; /*  >0: the wedge (pointed) end is at this atom, a1 */
            at[a2].bond_stereo[n2] = -cStereo; /*  <0: the wedge (pointed) end is at the opposite atom, a1 */
        } else
        if ( cBondStereo ) {
            *err |= 16; /*  Ignored unrecognized Bond stereo */
            MOLFILE_ERR_SET (*err, 0, "Unrecognized bond stereo");
            continue;
        }
    }
    *num_bonds = bonds;


    /* special valences */
    calculate_valences (mol_data, at, num_atoms, bDoNotAddH, err, pStrErr);

exit_function:;
    return at;
}
/******************************************************************************************************/
void calculate_valences (MOL_DATA* mol_data, inp_ATOM* at, int *num_atoms, int bDoNotAddH, int *err, char *pStrErr)
{
    int bNonMetal;
    int a1, a2, n1, n2, valence;
    AT_NUMB  *p1;

    /* special valences */
    for ( bNonMetal = 0; bNonMetal < 2; bNonMetal ++ ) {
        for ( a1 = 0; a1 < *num_atoms; a1 ++ ) {
            int num_bond_type[MAX_INPUT_BOND_TYPE - MIN_INPUT_BOND_TYPE + 1], bond_type, bHasMetalNeighbor;
            /* should the "!=" be replaced with "==" ??? */
            if ( bNonMetal == is_el_a_metal( at[a1].el_number ) ) {
                continue; /* first process all metals, after that all non-metals */
            }
            memset( num_bond_type, 0, sizeof(num_bond_type) );

            valence = at[a1].chem_bonds_valence; /*  save atom valence if available */

            at[a1].chem_bonds_valence = 0;
            bHasMetalNeighbor = 0;
            for ( n1 = 0; n1 < at[a1].valence; n1 ++ ) {
                bond_type = at[a1].bond_type[n1] - MIN_INPUT_BOND_TYPE;
                if (  bond_type < 0 || bond_type > MAX_INPUT_BOND_TYPE - MIN_INPUT_BOND_TYPE ) {
                    bond_type = 0;
                    MOLFILE_ERR_SET (*err, 0, "Unknown bond type in MOLfile assigned as a single bond");
                }
                num_bond_type[ bond_type ] ++;
                /* -- too a radical solution -- removed from next to ver 1.12B --- */
            }
            for ( n1 = 0; MIN_INPUT_BOND_TYPE + n1 <= 3 && MIN_INPUT_BOND_TYPE + n1 <= MAX_INPUT_BOND_TYPE; n1 ++ ) {
                /* add all bond orders except for "aromatic" bonds */
                at[a1].chem_bonds_valence += (MIN_INPUT_BOND_TYPE + n1) * num_bond_type[n1];
            }
            n2 = 0;
            if ( MIN_INPUT_BOND_TYPE <= BOND_TYPE_ALTERN && BOND_TYPE_ALTERN <= MAX_INPUT_BOND_TYPE &&
                 ( n2 = num_bond_type[BOND_TYPE_ALTERN-MIN_INPUT_BOND_TYPE] ) ) {
                /* accept input aromatic bonds for now */
                switch ( n2 ) {
                case 2:
                    at[a1].chem_bonds_valence += 3;  /* =A- */
                    break;
                case 3:
                    at[a1].chem_bonds_valence += 4;  /* =A< */
                    break;
                default:
                    /*  if 1 or >= 4 aromatic bonds then replace such bonds with single bonds */
                    /*  and detect an error in the input structure */
                    for ( n1 = 0; n1 < at[a1].valence; n1 ++ ) {
                        if ( at[a1].bond_type[n1] == BOND_TYPE_ALTERN ) {
                            a2 = at[a1].neighbor[n1];
                            p1 = is_in_the_list( at[a2].neighbor, (AT_NUMB)a1, at[a2].valence );
                            if ( p1 ) {
                                at[a1].bond_type[n1] =
                                at[a2].bond_type[p1-at[a2].neighbor] = BOND_TYPE_SINGLE;
                            } else {
                                *err = -2;  /*  Program error */
                                MOLFILE_ERR_SET (*err, 0, "Program error interpreting MOLfile");
                                return; /*  no structure */
                            }
                        }
                    }
                    at[a1].chem_bonds_valence += n2;
                    *err |= 32;
                    MOLFILE_ERR_SET (*err, 0, "Atom has 1 or more than 3 aromatic bonds");
                    n2 = 0;
                    break;
                }
            }
            if ( n2 && !valence ) {
                /* atom has aromatic bonds AND the chemical valence is not known */
                int num_H = NUMH(at, a1);
                int chem_valence = at[a1].chem_bonds_valence + num_H;
                int bUnusualValenceArom = 
                    detect_unusual_el_valence( (int)at[a1].el_number, at[a1].charge,
                                                at[a1].radical, chem_valence,
                                                num_H, at[a1].valence );
                int bUnusualValenceNoArom = 
                    detect_unusual_el_valence( (int)at[a1].el_number, at[a1].charge,
                                                at[a1].radical, chem_valence-1,
                                                num_H, at[a1].valence );
#if ( CHECK_AROMBOND2ALT == 1 )
                if ( bUnusualValenceArom && !bUnusualValenceNoArom && 0 == nBondsValToMetal( at, a1) )
#else
                if ( bUnusualValenceArom && !bUnusualValenceNoArom )
#endif                     
                {
                    /* typically NH in 5-member aromatic ring */
                    at[a1].chem_bonds_valence --;
                }
            } else
            if ( n2 && valence ) {
                /* atom has aromatic bonds AND the chemical valence is known */
                int num_H = NUMH(at, a1);
                int chem_valence = at[a1].chem_bonds_valence + num_H;
                if ( valence == chem_valence-1 ) {
                    /* typically NH in 5-member aromatic ring */
                    at[a1].chem_bonds_valence --;
                }
            }

            /*************************************************************************************
             *
             *  Set number of hydrogen atoms
             */
            if (mol_data) {
                at[a1].num_H = get_num_H( at[a1].elname, at[a1].num_H, at[a1].num_iso_H,
                                        at[a1].charge, at[a1].radical,
                                        at[a1].chem_bonds_valence,
                                        mol_data->ctab.MolAtom[a1].cValence, /* instead of valence */
                                        mol_data->ctab.MolAtom[a1].cAtomAliasedFlag,
                                        bDoNotAddH, bHasMetalNeighbor );
            }
        }
    }
}
/******************************************************************************************************/
int mol_to_atom_xyz( MOL_DATA* mol_data, int num_atoms, inp_ATOM* at, int *err, char *pStrErr )
{
    int      i, num_dimensions=0;
    int      num_bonds;
    double max_x=-1.0e32, max_y=-1.0e32, max_z=-1.0e32;
    double min_x= 1.0e32, min_y= 1.0e32, min_z= 1.0e32;
    double macheps = 1.0e-10, small_coeff = 0.00001;
    double x_coeff, y_coeff, z_coeff, coeff, average_bond_length;

    /*  *err = 0; */
    /* check if MOLfile contains atoms */
    if ( !mol_data || !mol_data->ctab.MolAtom ||
         0 < mol_data->ctab.nNumberOfBonds && !mol_data->ctab.MolBond ||
         0 >= (num_atoms = mol_data->ctab.nNumberOfAtoms) ) {
        goto exit_function; /*  no structure */
    }
    /* copy atom info */
    for ( i = 0; i < num_atoms; i ++ ) {
        max_x = inchi_max(mol_data->ctab.MolAtom[i].fX, max_x);
        min_x = inchi_min(mol_data->ctab.MolAtom[i].fX, min_x);
        max_y = inchi_max(mol_data->ctab.MolAtom[i].fY, max_y);
        min_y = inchi_min(mol_data->ctab.MolAtom[i].fY, min_y);
        max_z = inchi_max(mol_data->ctab.MolAtom[i].fZ, max_z);
        min_z = inchi_min(mol_data->ctab.MolAtom[i].fZ, min_z);
    }

    /* copy bond info */
    num_bonds = 0;
    average_bond_length = 0.0;
    for ( i = 0; i < mol_data->ctab.nNumberOfBonds; i ++ ) {
        int  a1 = mol_data->ctab.MolBond[i].nAtomNo1-1;
        int  a2 = mol_data->ctab.MolBond[i].nAtomNo2-1;
        double dx = mol_data->ctab.MolAtom[a1].fX-mol_data->ctab.MolAtom[a2].fX;
        double dy = mol_data->ctab.MolAtom[a1].fY-mol_data->ctab.MolAtom[a2].fY;
        double dz = mol_data->ctab.MolAtom[a1].fZ-mol_data->ctab.MolAtom[a2].fZ;

        if ( a1 < 0 || a1 >= num_atoms ||
             a2 < 0 || a2 >= num_atoms ||
             a1 == a2 ) {
            *err |= 1; /*  bond for impossible atom number(s); ignored */
            MOLFILE_ERR_SET (*err, 0, "Bond to nonexistent atom");
            continue;
        }
        average_bond_length += sqrt( dx*dx + dy*dy + dz*dz );
        num_bonds ++;
    }

    /* convert to integral coordinates */

    if ( max_x - min_x <= small_coeff*(fabs(max_x) + fabs(min_x)) )
        x_coeff = 0.0;
    else
        x_coeff = (MAX_STDATA_X_COORD - MIN_STDATA_X_COORD)/(max_x - min_x);

    if ( max_y - min_y <= small_coeff*(fabs(max_y) + fabs(min_y)) )
        y_coeff = 0.0;
    else
        y_coeff = (MAX_STDATA_Y_COORD - MIN_STDATA_Y_COORD)/(max_y - min_y);
    if ( max_z - min_z <= small_coeff*(fabs(max_z) + fabs(min_z)) )
        z_coeff = 0.0;
    else
        z_coeff = (MAX_STDATA_Z_COORD - MIN_STDATA_Z_COORD)/(max_z - min_z);

    num_dimensions = ((x_coeff > macheps || y_coeff >macheps ) && fabs(z_coeff) < macheps)? 2:
                     (fabs(z_coeff) > macheps)? 3: 0;

    switch ( num_dimensions ) {
    case 0:
        coeff = 0.0;
        break;
    case 2:
        /* choose the smallest stretching coefficient */
        if ( x_coeff > macheps && y_coeff > macheps ) {
            coeff = inchi_min( x_coeff, y_coeff );
        }else
        if ( x_coeff > macheps ){
            coeff = x_coeff;
        }else
        if ( y_coeff > macheps ){
            coeff = y_coeff;
        }else{
            coeff = 1.0;
        }
        break;
    case 3:
        /* choose the smallest stretching coefficient */
        if ( x_coeff > macheps && y_coeff > macheps ) {
            coeff = inchi_min( x_coeff, y_coeff );
            coeff = inchi_min( coeff, z_coeff );
        }else
        if ( x_coeff > macheps ){
            coeff = inchi_min( x_coeff, z_coeff );
        }else
        if ( y_coeff > macheps ){
            coeff = inchi_min( y_coeff, z_coeff );
        }else{
            coeff = z_coeff;
        }
        break;
    default:
        coeff = 0.0;
    }

    if ( num_bonds > 0 ) {
        average_bond_length /= (double)num_bonds;
        if ( average_bond_length * coeff > MAX_STDATA_AVE_BOND_LENGTH ) {
            coeff = MAX_STDATA_AVE_BOND_LENGTH / average_bond_length; /* avoid too long bonds */
        } else
        if ( average_bond_length * coeff < macheps ) {
            coeff = 1.0; /* all lengths are of zero length */
        } else
        if ( average_bond_length * coeff < MIN_STDATA_AVE_BOND_LENGTH ) {
            coeff = MIN_STDATA_AVE_BOND_LENGTH / average_bond_length; /* avoid too short bonds */
        }
    }
#if( NORMALIZE_INP_COORD == 1 )
    /* set integral coordinates */
    for ( i = 0; i < num_atoms; i ++ ) {
        double x = mol_data->ctab.MolAtom[i].fX;
        double y = mol_data->ctab.MolAtom[i].fY;
        double z = mol_data->ctab.MolAtom[i].fZ;
        x = (x - min_x)*coeff + MIN_STDATA_X_COORD;
        y = (y - min_y)*coeff + MIN_STDATA_Y_COORD;
        z = (z - min_z)*coeff + MIN_STDATA_Z_COORD;
        /* floor() behavior is not well defined for negative arguments.
         * Use positive arguments only to get nearest integer.
         */
        at[i].x = ( x >= 0.0 )? (int)floor( x + 0.5 ) : -(int)floor( -x + 0.5 );
        at[i].y = ( y >= 0.0 )? (int)floor( y + 0.5 ) : -(int)floor( -y + 0.5 );
        at[i].z = ( z >= 0.0 )? (int)floor( z + 0.5 ) : -(int)floor( -z + 0.5 );
    }
#else
    /* set input coordinates */
    for ( i = 0; i < num_atoms; i ++ ) {
        double x = mol_data->ctab.MolAtom[i].fX;
        double y = mol_data->ctab.MolAtom[i].fY;
        double z = mol_data->ctab.MolAtom[i].fZ;
        at[i].x = x;
        at[i].y = y;
        at[i].z = z;
    }
#endif

exit_function:;
    return num_dimensions;
}
/****************************************************************************/
long GetMolfileNumber( MOL_HEADER_BLOCK *pHdr )
{
    static char sStruct[] = "Structure #";
    static char sINCHI[]  = INCHI_NAME;
    long   lMolfileNumber = 0;
    char   *p, *q = NULL;
    if ( pHdr ) {
        if ( !memicmp( pHdr->szMoleculeName, sStruct, sizeof(sStruct)-1 ) ) {
            p = pHdr->szMoleculeName + sizeof(sStruct)-1;
            lMolfileNumber = strtol( p, &q, 10 );
            p = pHdr->szMoleculeLine2;
            if ( !q || *q ||
                 memicmp( p, sINCHI, sizeof(sINCHI)-1) ||
                 !strstr( p+sizeof(sINCHI)-1, "SDfile Output" ) ) {
                lMolfileNumber = 0;
            }
        }
    }
    return lMolfileNumber;
}
/****************************************************************************/
int MolfileToInpAtom( FILE *inp_molfile, int bDoNotAddH, inp_ATOM **at, MOL_COORD **szCoord, int max_num_at,
                      int *num_dimensions, int *num_bonds, const char *pSdfLabel, char *pSdfValue,
                      long *Id, long *lMolfileNumber, INCHI_MODE *pInpAtomFlags, int *err, char *pStrErr )
{
    int      num_atoms = 0;
    MOL_DATA *mol_data = NULL;
    MOL_HEADER_BLOCK OnlyHeaderBlock, *pOnlyHeaderBlock = NULL, *pHdr;
    MOL_CTAB         OnlyCtab,        *pOnlyCtab = NULL;
    char             cSdfValueFirstChar;
#ifdef CML_DEBUG
        FILE *f_p;
#endif
    if ( at ) {
        pOnlyHeaderBlock = NULL;
        if ( *at && max_num_at ) {
            memset( *at, 0, max_num_at * sizeof(inp_ATOM) );
        }
        if ( szCoord && *szCoord ) {
            inchi_free( *szCoord );
            *szCoord = NULL;
        }
    } else {
        pOnlyHeaderBlock = &OnlyHeaderBlock;
        pOnlyCtab        = &OnlyCtab;
    }
    if ( pSdfValue ) {
        cSdfValueFirstChar = pSdfValue[0];
        pSdfValue[0] = '\0';
    }

    mol_data = read_sdfile_segment(inp_molfile, pOnlyHeaderBlock, pOnlyCtab, NULL != szCoord,
                                   NULL, 0, Id, pSdfLabel, pSdfValue, err, pStrErr );

    pHdr = (  mol_data && !pOnlyHeaderBlock )? &mol_data->hdr :
           (  !mol_data && pOnlyHeaderBlock )? pOnlyHeaderBlock : NULL;
    if ( lMolfileNumber && pHdr ) {
         *lMolfileNumber = GetMolfileNumber( pHdr );
    }
    if ( pSdfValue && !pSdfValue[0] &&
         pSdfLabel && pSdfLabel[0]  && pHdr ) {
         if ( !stricmp(pSdfLabel, "MolfileName") ) {
             mystrncpy( pSdfValue, pHdr->szMoleculeName, MAX_SDF_VALUE+1 );
             LtrimRtrim( pSdfValue, NULL );
         } else
         if ( !stricmp(pSdfLabel, "MolfileLine2") ) {
             mystrncpy( pSdfValue, pHdr->szMoleculeLine2, MAX_SDF_VALUE+1 );
             LtrimRtrim( pSdfValue, NULL );
         } else
         if ( !stricmp(pSdfLabel, "MolfileComment") ) {
             mystrncpy( pSdfValue, pHdr->szComment, MAX_SDF_VALUE+1 );
             LtrimRtrim( pSdfValue, NULL );
         } else
         if ( !stricmp(pSdfLabel, "MolfileIntRegNo") && pHdr->lInternalRegistryNumber ) {
             sprintf( pSdfValue, "%ld", pHdr->lInternalRegistryNumber );
         }
         if ( !pSdfValue[0] ) {
             pSdfValue[0] = cSdfValueFirstChar;
         }
    }

    if ( mol_data && at && !*err ) {
        /* *at points to an allocated memory */
        if ( *at && mol_data->ctab.nNumberOfAtoms <= max_num_at ) {
            *at = mol_to_atom( mol_data, &num_atoms, num_bonds, *at, bDoNotAddH, err, pStrErr );
            if ( *err >= 0 ) {
                *num_dimensions = mol_to_atom_xyz( mol_data, num_atoms, *at, err, pStrErr );

                if ( szCoord ) {
                    *szCoord = mol_data->ctab.szCoord;
                    mol_data->ctab.szCoord = NULL;
                }

            }
        } else
        /* *at points to NULL */
        if ( !*at && mol_data->ctab.nNumberOfAtoms <= max_num_at ) {
            *at = mol_to_atom( mol_data, &num_atoms, num_bonds, *at, bDoNotAddH, err, pStrErr );
            if ( *err >= 0 ) {
                *num_dimensions = mol_to_atom_xyz( mol_data, num_atoms, *at, err, pStrErr );

                if ( szCoord ) {
                    *szCoord = mol_data->ctab.szCoord;
                    mol_data->ctab.szCoord = NULL;
                }

            }
        } else {
            MOLFILE_ERR_SET (*err, 0, "Too many atoms");
            *err = 70;
            num_atoms = -1;
        }
        if ( *err > 0 ) {
            *err += 100;
        }
        /* 11-16-2004: use Chiral flag */
        if ( num_atoms > 0 && at && *at && mol_data && pInpAtomFlags ) {
            if ( mol_data->ctab.cChiralFlag ) {
                *pInpAtomFlags |= FLAG_INP_AT_CHIRAL;
            } else {
                *pInpAtomFlags |= FLAG_INP_AT_NONCHIRAL;
            }
        }
    } else
    if ( !at ) {
        num_atoms = pOnlyCtab->nNumberOfAtoms;
    }

    if ( !pOnlyHeaderBlock ) {
        delete_mol_data( mol_data );
    }
#ifdef CML_DEBUG
        puts ("MOL");
        f_p = fopen ("mol.dbg", "a");
        if (f_p)
            {
                PrintInpAtom (f_p, *at, num_atoms);
                fclose (f_p);
            }
        else
            {
                puts ("Couldn't open file");
            }
#endif

    return num_atoms;
}
/**********************************************************************************/
void FreeOrigAtData( ORIG_ATOM_DATA *orig_at_data )
{
    if ( !orig_at_data )
        return;
    FreeInpAtom( &orig_at_data->at );
    if ( orig_at_data->nCurAtLen ) {
        inchi_free( orig_at_data->nCurAtLen );
    }
    if ( orig_at_data->nOldCompNumber ) {
        inchi_free( orig_at_data->nOldCompNumber );
    }
    if ( orig_at_data->szCoord ) {
        inchi_free( orig_at_data->szCoord );
    }

    if ( orig_at_data->nEquLabels ) {
        inchi_free( orig_at_data->nEquLabels );
    }
    if ( orig_at_data->nSortedOrder ) {
        inchi_free( orig_at_data->nSortedOrder );
    }

    memset( orig_at_data, 0, sizeof(*orig_at_data) );
}
/**********************************************************************************/
int MolfileToOrigAtom( FILE *inp_molfile, ORIG_ATOM_DATA *orig_at_data, int bMergeAllInputStructures,
                       int bGetOrigCoord, int bDoNotAddH,
                       const char *pSdfLabel, char *pSdfValue, long *lSdfId, long *lMolfileNumber,
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

    do {

        at_old     = orig_at_data? orig_at_data->at      : NULL; /*  save pointer to the previous allocation */
        szCoordOld = orig_at_data? orig_at_data->szCoord : NULL;
        num_inp_atoms_new =
            MolfileToInpAtom( inp_molfile, bDoNotAddH, orig_at_data? &at_new:NULL, (bGetOrigCoord && orig_at_data)? &szCoordNew : NULL, MAX_ATOMS,
                          &num_dimensions_new, &num_inp_bonds_new,
                          pSdfLabel, pSdfValue, lSdfId, lMolfileNumber, pInpAtomFlags, err, pStrErr );


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
    if ( *err ) {
        FreeOrigAtData( orig_at_data );
    }
    if ( *err && !(10 < *err && *err < 20) && pStrErr && !pStrErr[0] ) {
        MOLFILE_ERR_SET (*err, 0, "Unknown error");  /*   <BRKPT> */
    }
    return orig_at_data? orig_at_data->num_inp_atoms : nNumAtoms;
}



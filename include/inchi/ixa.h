/*
 * International Chemical Identifier (InChI)
 * Version 1
 * Software version 1.07
 * April 30, 2024
 *
 * MIT License
 *
 * Copyright (c) 2024 IUPAC and InChI Trust
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
*
* The InChI library and programs are free software developed under the
 * auspices of the International Union of Pure and Applied Chemistry (IUPAC).
 * Originally developed at NIST.
 * Modifications and additions by IUPAC and the InChI Trust.
 * Some portions of code were developed/changed by external contributors
 * (either contractor or volunteer) which are listed in the file
 * 'External-contributors' included in this distribution.
 *
 * info@inchi-trust.org
 *
*/


#ifndef __IXA_H__
#define __IXA_H__


/****************************************************************************/
/*                 InChI Extensible API Object Handles                      */
/*
 * These are the "handles" which can be used to refer to "Objects"
 * used in the InChI Extensible API.
 *
 * To ensure that each handle has a different formal type, each is declared
 * as a pointer to a C struct containing a dummy integer field (which is in
 * fact not used). In practice, the values are simply cast from int, but
 * this approach provides greater type security, and prevents handles for
 * different objects getting mixed up.
*/

typedef struct { int dummy; } IXA_STATUS_HANDLE_STRUCT;
typedef IXA_STATUS_HANDLE_STRUCT *IXA_STATUS_HANDLE;

typedef struct { int dummy; } IXA_MOL_HANDLE_STRUCT;
typedef IXA_MOL_HANDLE_STRUCT *IXA_MOL_HANDLE;

typedef struct { int dummy; } IXA_INCHIBUILDER_HANDLE_STRUCT;
typedef IXA_INCHIBUILDER_HANDLE_STRUCT *IXA_INCHIBUILDER_HANDLE;

typedef struct { int dummy; } IXA_INCHIKEYBUILDER_HANDLE_STRUCT;
typedef IXA_INCHIKEYBUILDER_HANDLE_STRUCT *IXA_INCHIKEYBUILDER_HANDLE;

/****************************************************************************/
/*        Types for Atom, Bond and Stereo Descriptor Identifiers            */
/*
 * These types are for the identifiers for individual atoms, bonds and
 * stereodescriptors in an IXA Molecule Object.
 *
 * To ensure that each identifier has a different formal type, each is
 * declared as a pointer to a C struct containing a dummy integer field
 * (which is in fact not used). In practice, the values are simply cast
 * from int, but this approach provides greater type security, and
 * prevents different sorts of identifier from getting mixed up.
*/

typedef struct { int dummy; } IXA_ATOMID_STRUCT;
typedef IXA_ATOMID_STRUCT *IXA_ATOMID;

typedef struct { int dummy; } IXA_BONDID_STRUCT;
typedef IXA_BONDID_STRUCT *IXA_BONDID;

typedef struct { int dummy; } IXA_STEREOID_STRUCT;
typedef IXA_STEREOID_STRUCT *IXA_STEREOID;

/* Extended mol data */

typedef struct { int dummy; } IXA_POLYMERUNITID_STRUCT;
typedef IXA_POLYMERUNITID_STRUCT *IXA_POLYMERUNITID;   /* polymer unit */

/****************************************************************************/
/*                Constants and enumerated types                            */

#define IXA_ATOMID_INVALID ((IXA_ATOMID)0)

#define IXA_ATOMID_IMPLICIT_H ((IXA_ATOMID)-1)

#define IXA_BONDID_INVALID ((IXA_BONDID)0)

#define IXA_STEREOID_INVALID ((IXA_STEREOID)0)

#define IXA_ATOM_NATURAL_MASS 0

#define IXA_POLYMERUNITID_INVALID ((IXA_POLYMERUNITID)0)
#define IXA_EXT_MOLDATA_INVALID (-1)
#define IXA_EXT_POLYMER_INVALID (-1)
#define IXA_EXT_V3000_INVALID (-1)

typedef enum
{
    IXA_STATUS_SUCCESS,
    IXA_STATUS_WARNING,
    IXA_STATUS_ERROR
} IXA_STATUS;

typedef enum
{
    IXA_FALSE = 0,
    IXA_TRUE = 1
} IXA_BOOL;

typedef enum
{
    IXA_ATOM_RADICAL_NONE = 0,
    IXA_ATOM_RADICAL_SINGLET = 1,
    IXA_ATOM_RADICAL_DOUBLET = 2,
    IXA_ATOM_RADICAL_TRIPLET = 3
} IXA_ATOM_RADICAL;

typedef enum
{
    IXA_BOND_TYPE_SINGLE = 1,
    IXA_BOND_TYPE_DOUBLE = 2,
    IXA_BOND_TYPE_TRIPLE = 3,
    IXA_BOND_TYPE_AROMATIC = 4
} IXA_BOND_TYPE;

typedef enum
{
    IXA_BOND_WEDGE_NONE = 0,
    IXA_BOND_WEDGE_UP = 1,
    IXA_BOND_WEDGE_DOWN = 2,
    IXA_BOND_WEDGE_EITHER = 3
} IXA_BOND_WEDGE;

typedef enum
{
    IXA_DBLBOND_CONFIG_PERCEIVE = 0,
    IXA_DBLBOND_CONFIG_EITHER = 1
} IXA_DBLBOND_CONFIG;

typedef enum
{
    IXA_STEREO_TOPOLOGY_INVALID = 0,
    IXA_STEREO_TOPOLOGY_TETRAHEDRON = 2,
    IXA_STEREO_TOPOLOGY_RECTANGLE = 3,
    IXA_STEREO_TOPOLOGY_ANTIRECTANGLE = 4
} IXA_STEREO_TOPOLOGY;

typedef enum
{
    IXA_STEREO_PARITY_NONE = 0,
    IXA_STEREO_PARITY_ODD = 1,
    IXA_STEREO_PARITY_EVEN = 2,
    IXA_STEREO_PARITY_UNKNOWN = 3
} IXA_STEREO_PARITY;

typedef enum
{
    IXA_INCHIBUILDER_OPTION_NewPsOff,
    IXA_INCHIBUILDER_OPTION_DoNotAddH,
    IXA_INCHIBUILDER_OPTION_SUU,
    IXA_INCHIBUILDER_OPTION_SLUUD,
    IXA_INCHIBUILDER_OPTION_FixedH,
    IXA_INCHIBUILDER_OPTION_RecMet,
    IXA_INCHIBUILDER_OPTION_KET,
    IXA_INCHIBUILDER_OPTION_15T,
    IXA_INCHIBUILDER_OPTION_SaveOpt,
    IXA_INCHIBUILDER_OPTION_AuxNone,
    IXA_INCHIBUILDER_OPTION_WarnOnEmptyStructure,
    IXA_INCHIBUILDER_OPTION_LargeMolecules,
    IXA_INCHIBUILDER_OPTION_Polymers,
    IXA_INCHIBUILDER_OPTION_Polymers105,
    IXA_INCHIBUILDER_OPTION_Polymers105Plus,
    IXA_INCHIBUILDER_OPTION_FilterSS,
    IXA_INCHIBUILDER_OPTION_InvFilterSS,
    IXA_INCHIBUILDER_OPTION_NPZZ,
    IXA_INCHIBUILDER_OPTION_SATZZ,
    IXA_INCHIBUILDER_OPTION_NoFrameShift,
    IXA_INCHIBUILDER_OPTION_FoldCRU,
    IXA_INCHIBUILDER_OPTION_NoEdits,
    IXA_INCHIBUILDER_OPTION_LooseTSACheck,
    IXA_INCHIBUILDER_OPTION_OutErrInChI,
    IXA_INCHIBUILDER_OPTION_NoWarnings
/*#if BUILD_WITH_ENG_OPTIONS==1*/
    ,IXA_INCHIBUILDER_OPTION_DoDrv,
    IXA_INCHIBUILDER_OPTION_DoDrvReport,
    IXA_INCHIBUILDER_OPTION_DoR2C,
    IXA_INCHIBUILDER_OPTION_DoneOnly,
    IXA_INCHIBUILDER_OPTION_OnlyRecSalt,
    IXA_INCHIBUILDER_OPTION_OnlyExact,
    IXA_INCHIBUILDER_OPTION_OnlyRecMet
/*#endif*/   
    , IXA_INCHIBUILDER_OPTION_PT_22_00,
    IXA_INCHIBUILDER_OPTION_PT_16_00,
    IXA_INCHIBUILDER_OPTION_PT_06_00,
    IXA_INCHIBUILDER_OPTION_PT_39_00,
    IXA_INCHIBUILDER_OPTION_PT_13_00,
    IXA_INCHIBUILDER_OPTION_PT_18_00
} IXA_INCHIBUILDER_OPTION;

typedef enum
{
    IXA_INCHIBUILDER_STEREOOPTION_SAbs,
    IXA_INCHIBUILDER_STEREOOPTION_SNon,
    IXA_INCHIBUILDER_STEREOOPTION_SRel,
    IXA_INCHIBUILDER_STEREOOPTION_SRac,
    IXA_INCHIBUILDER_STEREOOPTION_SUCF
} IXA_INCHIBUILDER_STEREOOPTION;



/* Uncomment the next line if old-API coverage is intended - instead of GetINCHIEx(), GetStructFromINCHIEx() */
/*#define IXA_USES_NON_EX_CORE_API 1*/

/* Comment the next line to disable improved storage mechanism (pre-allocated growing arrays) for IXA_MOL   */
/* data structs (added after report by Daniel Lowe on performance deterioraion due to numerous alloc calls) */
#define IXA_USES_SMART_ALLOCS 1


/****************************************************************************/

#ifdef __cplusplus
extern "C" {
#endif

/****************************************************************************/
/*              Functions handling IXA Status Objects                      */

    EXPIMP_TEMPLATE INCHI_API IXA_STATUS_HANDLE INCHI_DECL IXA_STATUS_Create( void );

    EXPIMP_TEMPLATE INCHI_API void INCHI_DECL IXA_STATUS_Clear( IXA_STATUS_HANDLE hStatus );

    EXPIMP_TEMPLATE INCHI_API void INCHI_DECL IXA_STATUS_Destroy( IXA_STATUS_HANDLE hStatus );

    EXPIMP_TEMPLATE INCHI_API IXA_BOOL INCHI_DECL IXA_STATUS_HasError( IXA_STATUS_HANDLE hStatus );

    EXPIMP_TEMPLATE INCHI_API IXA_BOOL INCHI_DECL IXA_STATUS_HasWarning( IXA_STATUS_HANDLE hStatus );

    EXPIMP_TEMPLATE INCHI_API int INCHI_DECL IXA_STATUS_GetCount( IXA_STATUS_HANDLE hStatus );

    EXPIMP_TEMPLATE INCHI_API IXA_STATUS INCHI_DECL IXA_STATUS_GetSeverity( IXA_STATUS_HANDLE hStatus,
                                                                           int               vIndex );

    EXPIMP_TEMPLATE INCHI_API const char* INCHI_DECL IXA_STATUS_GetMessage( IXA_STATUS_HANDLE hStatus,
                                                                           int               vIndex );




    /****************************************************************************/
    /*       Functions to Create, Clear and Destroy Molecule Objects            */

    EXPIMP_TEMPLATE INCHI_API IXA_MOL_HANDLE INCHI_DECL IXA_MOL_Create( IXA_STATUS_HANDLE hStatus );

    EXPIMP_TEMPLATE INCHI_API void INCHI_DECL IXA_MOL_Clear( IXA_STATUS_HANDLE hStatus,
                                                            IXA_MOL_HANDLE    hMolecule );

    EXPIMP_TEMPLATE INCHI_API void INCHI_DECL IXA_MOL_Destroy( IXA_STATUS_HANDLE hStatus,
                                                              IXA_MOL_HANDLE    hMolecule );


    /****************************************************************************/
    /*               Functions Operating on Complete Molecules                  */

    EXPIMP_TEMPLATE INCHI_API void INCHI_DECL IXA_MOL_ReadMolfile( IXA_STATUS_HANDLE hStatus,
                                                                   IXA_MOL_HANDLE    hMolecule,
                                                                   const char*       pBytes );

    
    EXPIMP_TEMPLATE INCHI_API void INCHI_DECL IXA_MOL_ReadInChI( IXA_STATUS_HANDLE hStatus,
                                                                 IXA_MOL_HANDLE    hMolecule,
                                                                 const char*       pInChI );

    EXPIMP_TEMPLATE INCHI_API void INCHI_DECL IXA_MOL_ReadAuxInfo(IXA_STATUS_HANDLE hStatus,
                                                                    IXA_MOL_HANDLE    hMolecule,
                                                                    char* pAuxInfo,
                                                                    int               bDoNotAddH,
                                                                    int               bDiffUnkUndfStereo);

    EXPIMP_TEMPLATE INCHI_API void INCHI_DECL IXA_MOL_SetChiral( IXA_STATUS_HANDLE hStatus,
                                                                 IXA_MOL_HANDLE    hMolecule,
                                                                 IXA_BOOL          vChiral );

    EXPIMP_TEMPLATE INCHI_API IXA_BOOL INCHI_DECL IXA_MOL_GetChiral( IXA_STATUS_HANDLE hStatus,
                                                                     IXA_MOL_HANDLE    hMolecule );


    /****************************************************************************/
    /*               Functions to Add and Define Atoms                          */

    EXPIMP_TEMPLATE INCHI_API IXA_ATOMID INCHI_DECL IXA_MOL_CreateAtom( IXA_STATUS_HANDLE hStatus,
                                                                        IXA_MOL_HANDLE    hMolecule );

    EXPIMP_TEMPLATE INCHI_API void INCHI_DECL IXA_MOL_SetAtomElement( IXA_STATUS_HANDLE hStatus,
                                                                      IXA_MOL_HANDLE    hMolecule,
                                                                      IXA_ATOMID        vAtom,
                                                                      const char*       pElement );

    EXPIMP_TEMPLATE INCHI_API void INCHI_DECL IXA_MOL_SetAtomAtomicNumber( IXA_STATUS_HANDLE hStatus,
                                                                           IXA_MOL_HANDLE    hMolecule,
                                                                           IXA_ATOMID        vAtom,
                                                                           int               vAtomicNumber );

    EXPIMP_TEMPLATE INCHI_API void INCHI_DECL IXA_MOL_SetAtomMass( IXA_STATUS_HANDLE hStatus,
                                                                   IXA_MOL_HANDLE    hMolecule,
                                                                   IXA_ATOMID        vAtom,
                                                                   int               vMassNumber );

    EXPIMP_TEMPLATE INCHI_API void INCHI_DECL IXA_MOL_SetAtomCharge( IXA_STATUS_HANDLE hStatus,
                                                                     IXA_MOL_HANDLE    hMolecule,
                                                                     IXA_ATOMID        vAtom,
                                                                     int               vCharge );

    EXPIMP_TEMPLATE INCHI_API void INCHI_DECL IXA_MOL_SetAtomRadical( IXA_STATUS_HANDLE hStatus,
                                                                     IXA_MOL_HANDLE    hMolecule,
                                                                     IXA_ATOMID        vAtom,
                                                                     IXA_ATOM_RADICAL  vRadical );

    EXPIMP_TEMPLATE INCHI_API void INCHI_DECL IXA_MOL_SetAtomHydrogens( IXA_STATUS_HANDLE hStatus,
                                                                       IXA_MOL_HANDLE    hMolecule,
                                                                       IXA_ATOMID        vAtom,
                                                                       int               vHydrogenMassNumber,
                                                                       int               vHydrogenCount );

    EXPIMP_TEMPLATE INCHI_API void INCHI_DECL IXA_MOL_SetAtomX( IXA_STATUS_HANDLE hStatus,
                                                               IXA_MOL_HANDLE    hMolecule,
                                                               IXA_ATOMID        vAtom,
                                                               double            vX );

    EXPIMP_TEMPLATE INCHI_API void INCHI_DECL IXA_MOL_SetAtomY( IXA_STATUS_HANDLE hStatus,
                                                               IXA_MOL_HANDLE    hMolecule,
                                                               IXA_ATOMID        vAtom,
                                                               double            vY );

    EXPIMP_TEMPLATE INCHI_API void INCHI_DECL IXA_MOL_SetAtomZ( IXA_STATUS_HANDLE hStatus,
                                                               IXA_MOL_HANDLE    hMolecule,
                                                               IXA_ATOMID        vAtom,
                                                               double            vZ );


    /****************************************************************************/
    /*               Functions to Add and Define Bonds                          */

    EXPIMP_TEMPLATE INCHI_API IXA_BONDID INCHI_DECL IXA_MOL_CreateBond( IXA_STATUS_HANDLE hStatus,
                                                                       IXA_MOL_HANDLE    hMolecule,
                                                                       IXA_ATOMID        vAtom1,
                                                                       IXA_ATOMID        vAtom2 );

    EXPIMP_TEMPLATE INCHI_API void INCHI_DECL IXA_MOL_SetBondType( IXA_STATUS_HANDLE hStatus,
                                                                  IXA_MOL_HANDLE    hMolecule,
                                                                  IXA_BONDID        vBond,
                                                                  IXA_BOND_TYPE     vType );

    EXPIMP_TEMPLATE INCHI_API void INCHI_DECL IXA_MOL_SetBondWedge( IXA_STATUS_HANDLE hStatus,
                                                                   IXA_MOL_HANDLE    hMolecule,
                                                                   IXA_BONDID        vBond,
                                                                   IXA_ATOMID        vRefAtom,
                                                                   IXA_BOND_WEDGE    vDirection );

    EXPIMP_TEMPLATE INCHI_API void INCHI_DECL IXA_MOL_SetDblBondConfig( IXA_STATUS_HANDLE  hStatus,
                                                                       IXA_MOL_HANDLE     hMolecule,
                                                                       IXA_BONDID         vBond,
                                                                       IXA_DBLBOND_CONFIG vConfig );


    /*****************************************************************************/
    /*              Functions to Add and Define Stereodescriptors                */

    EXPIMP_TEMPLATE INCHI_API IXA_STEREOID INCHI_DECL IXA_MOL_CreateStereoTetrahedron( IXA_STATUS_HANDLE hStatus,
                                                                                      IXA_MOL_HANDLE    hMolecule,
                                                                                      IXA_ATOMID        vCentralAtom,
                                                                                      IXA_ATOMID        vVertex1,
                                                                                      IXA_ATOMID        vVertex2,
                                                                                      IXA_ATOMID        vVertex3,
                                                                                      IXA_ATOMID        vVertex4 );

    EXPIMP_TEMPLATE INCHI_API IXA_STEREOID INCHI_DECL IXA_MOL_CreateStereoRectangle( IXA_STATUS_HANDLE hStatus,
                                                                                    IXA_MOL_HANDLE    hMolecule,
                                                                                    IXA_BONDID        vCentralBond,
                                                                                    IXA_ATOMID        vVertex1,
                                                                                    IXA_ATOMID        vVertex2,
                                                                                    IXA_ATOMID        vVertex3,
                                                                                    IXA_ATOMID        vVertex4 );

    EXPIMP_TEMPLATE INCHI_API IXA_STEREOID INCHI_DECL IXA_MOL_CreateStereoAntiRectangle( IXA_STATUS_HANDLE hStatus,
                                                                                        IXA_MOL_HANDLE    hMolecule,
                                                                                        IXA_ATOMID        vCentralAtom,
                                                                                        IXA_ATOMID        vVertex1,
                                                                                        IXA_ATOMID        vVertex2,
                                                                                        IXA_ATOMID        vVertex3,
                                                                                        IXA_ATOMID        vVertex4 );

    EXPIMP_TEMPLATE INCHI_API void INCHI_DECL IXA_MOL_SetStereoParity( IXA_STATUS_HANDLE hStatus,
                                                                      IXA_MOL_HANDLE    hMolecule,
                                                                      IXA_STEREOID      vStereo,
                                                                      IXA_STEREO_PARITY vParity );


    EXPIMP_TEMPLATE INCHI_API int INCHI_DECL IXA_MOL_ReserveSpace( IXA_STATUS_HANDLE hStatus,
                                                                   IXA_MOL_HANDLE    hMolecule,
                                                                   int num_atoms,
                                                                   int num_bonds,
                                                                   int num_stereos );


    /****************************************************************************/
    /*               Functions to to Treat Extended molecular data              */

    EXPIMP_TEMPLATE INCHI_API IXA_POLYMERUNITID INCHI_DECL IXA_MOL_CreatePolymerUnit( IXA_STATUS_HANDLE hStatus,
                                                                                      IXA_MOL_HANDLE    hMolecule );

    EXPIMP_TEMPLATE INCHI_API void INCHI_DECL IXA_MOL_SetPolymerUnit( IXA_STATUS_HANDLE hStatus,
                                                                      IXA_MOL_HANDLE    hMolecule,
                                                                      IXA_POLYMERUNITID vPunit,
                                                                      int               vid,
                                                                      int               vtype,
                                                                      int               vsubtype,
                                                                      int               vconn,
                                                                      int               vlabel,
                                                                      int               vna,
                                                                      int               vnb,
                                                                      double            vxbr1[4],
                                                                      double            vxbr2[4],
                                                                      char              vsmt[80],
                                                                      int               *valist,
                                                                      int               *vblist );



    /****************************************************************************/
    /*               Functions to Navigate Within a Molecule                    */

    EXPIMP_TEMPLATE INCHI_API int INCHI_DECL IXA_MOL_GetNumAtoms( IXA_STATUS_HANDLE hStatus,
                                                                  IXA_MOL_HANDLE    hMolecule );

    EXPIMP_TEMPLATE INCHI_API int INCHI_DECL IXA_MOL_GetNumBonds( IXA_STATUS_HANDLE hStatus,
                                                                  IXA_MOL_HANDLE    hMolecule );

    EXPIMP_TEMPLATE INCHI_API IXA_ATOMID INCHI_DECL IXA_MOL_GetAtomId( IXA_STATUS_HANDLE  hStatus,
                                                                       IXA_MOL_HANDLE hMolecule,
                                                                       int            vAtomIndex );

    EXPIMP_TEMPLATE INCHI_API IXA_BONDID INCHI_DECL IXA_MOL_GetBondId( IXA_STATUS_HANDLE  hStatus,
                                                                       IXA_MOL_HANDLE hMolecule,
                                                                       int            vBondIndex );

    EXPIMP_TEMPLATE INCHI_API int INCHI_DECL IXA_MOL_GetAtomIndex( IXA_STATUS_HANDLE hStatus,
                                                                   IXA_MOL_HANDLE    hMolecule,
                                                                   IXA_ATOMID        vAtom );

    EXPIMP_TEMPLATE INCHI_API int INCHI_DECL IXA_MOL_GetBondIndex( IXA_STATUS_HANDLE hStatus,
                                                                  IXA_MOL_HANDLE    hMolecule,
                                                                  IXA_BONDID        vBond );

    EXPIMP_TEMPLATE INCHI_API int INCHI_DECL IXA_MOL_GetAtomNumBonds( IXA_STATUS_HANDLE hStatus,
                                                                     IXA_MOL_HANDLE    hMolecule,
                                                                     IXA_ATOMID    vAtom );

    EXPIMP_TEMPLATE INCHI_API IXA_POLYMERUNITID INCHI_DECL IXA_MOL_GetPolymerUnitId( IXA_STATUS_HANDLE  hStatus,
                                                                                     IXA_MOL_HANDLE     hMolecule,
                                                                                     int                vPolymerUnitIndex );
    EXPIMP_TEMPLATE INCHI_API int INCHI_DECL IXA_MOL_GetPolymerUnitIndex( IXA_STATUS_HANDLE hStatus,
                                                                          IXA_MOL_HANDLE    hMolecule,
                                                                          IXA_POLYMERUNITID vPolymerUnit );


    EXPIMP_TEMPLATE INCHI_API IXA_BONDID INCHI_DECL IXA_MOL_GetAtomBond( IXA_STATUS_HANDLE  hStatus,
                                                                            IXA_MOL_HANDLE hMolecule,
                                                                            IXA_ATOMID     vAtom,
                                                                            int            vBondIndex );

    EXPIMP_TEMPLATE INCHI_API IXA_BONDID INCHI_DECL IXA_MOL_GetCommonBond( IXA_STATUS_HANDLE  hStatus,
                                                                              IXA_MOL_HANDLE hMolecule,
                                                                              IXA_ATOMID     vAtom1,
                                                                              IXA_ATOMID     vAtom2 );

    EXPIMP_TEMPLATE INCHI_API IXA_ATOMID INCHI_DECL IXA_MOL_GetBondAtom1( IXA_STATUS_HANDLE  hStatus,
                                                                             IXA_MOL_HANDLE hMolecule,
                                                                             IXA_BONDID     vBond );

    EXPIMP_TEMPLATE INCHI_API IXA_ATOMID INCHI_DECL IXA_MOL_GetBondAtom2( IXA_STATUS_HANDLE  hStatus,
                                                                             IXA_MOL_HANDLE hMolecule,
                                                                             IXA_BONDID     vBond );

    EXPIMP_TEMPLATE INCHI_API IXA_ATOMID INCHI_DECL IXA_MOL_GetBondOtherAtom( IXA_STATUS_HANDLE  hStatus,
                                                                             IXA_MOL_HANDLE     hMolecule,
                                                                             IXA_BONDID         vBond,
                                                                             IXA_ATOMID         vAtom );

    /*****************************************************************************/
    /*             Functions to Return Information About Atoms                   */

    EXPIMP_TEMPLATE INCHI_API const char* INCHI_DECL IXA_MOL_GetAtomElement( IXA_STATUS_HANDLE hStatus,
                                                                            IXA_MOL_HANDLE    hMolecule,
                                                                            IXA_ATOMID        vAtom );

    EXPIMP_TEMPLATE INCHI_API int INCHI_DECL IXA_MOL_GetAtomAtomicNumber( IXA_STATUS_HANDLE hStatus,
                                                                         IXA_MOL_HANDLE    hMolecule,
                                                                         IXA_ATOMID        vAtom );

    EXPIMP_TEMPLATE INCHI_API int INCHI_DECL IXA_MOL_GetAtomMass( IXA_STATUS_HANDLE hStatus,
                                                                 IXA_MOL_HANDLE    hMolecule,
                                                                 IXA_ATOMID        vAtom );

    EXPIMP_TEMPLATE INCHI_API int INCHI_DECL IXA_MOL_GetAtomCharge( IXA_STATUS_HANDLE hStatus,
                                                                   IXA_MOL_HANDLE    hMolecule,
                                                                   IXA_ATOMID        vAtom );

    EXPIMP_TEMPLATE INCHI_API IXA_ATOM_RADICAL INCHI_DECL IXA_MOL_GetAtomRadical( IXA_STATUS_HANDLE hStatus,
                                                                                 IXA_MOL_HANDLE    hMolecule,
                                                                                 IXA_ATOMID        vAtom );

    EXPIMP_TEMPLATE INCHI_API int INCHI_DECL IXA_MOL_GetAtomHydrogens( IXA_STATUS_HANDLE hStatus,
                                                                      IXA_MOL_HANDLE    hMolecule,
                                                                      IXA_ATOMID        vAtom,
                                                                      int               vHydrogenMassNumber );

    EXPIMP_TEMPLATE INCHI_API double INCHI_DECL IXA_MOL_GetAtomX( IXA_STATUS_HANDLE hStatus,
                                                                 IXA_MOL_HANDLE    hMolecule,
                                                                 IXA_ATOMID        vAtom );

    EXPIMP_TEMPLATE INCHI_API double INCHI_DECL IXA_MOL_GetAtomY( IXA_STATUS_HANDLE hStatus,
                                                                 IXA_MOL_HANDLE    hMolecule,
                                                                 IXA_ATOMID        vAtom );

    EXPIMP_TEMPLATE INCHI_API double INCHI_DECL IXA_MOL_GetAtomZ( IXA_STATUS_HANDLE hStatus,
                                                                 IXA_MOL_HANDLE    hMolecule,
                                                                 IXA_ATOMID        vAtom );


    /*****************************************************************************/
    /*             Functions to Return Information About Bonds                   */

    EXPIMP_TEMPLATE INCHI_API IXA_BOND_TYPE INCHI_DECL IXA_MOL_GetBondType( IXA_STATUS_HANDLE hStatus,
                                                                           IXA_MOL_HANDLE    hMolecule,
                                                                           IXA_BONDID        vBond );

    EXPIMP_TEMPLATE INCHI_API IXA_BOND_WEDGE INCHI_DECL IXA_MOL_GetBondWedge( IXA_STATUS_HANDLE hStatus,
                                                                             IXA_MOL_HANDLE    hMolecule,
                                                                             IXA_BONDID        vBond,
                                                                             IXA_ATOMID        vRefAtom );

    EXPIMP_TEMPLATE INCHI_API IXA_DBLBOND_CONFIG INCHI_DECL IXA_MOL_GetDblBondConfig( IXA_STATUS_HANDLE hStatus,
                                                                                        IXA_MOL_HANDLE hMolecule,
                                                                                        IXA_BONDID     vBond );

    /*****************************************************************************/
    /*      Functions to return Information About Stereodescriptors              */


    EXPIMP_TEMPLATE INCHI_API int INCHI_DECL IXA_MOL_GetNumStereos( IXA_STATUS_HANDLE hStatus,
                                                                   IXA_MOL_HANDLE    hMolecule );

    EXPIMP_TEMPLATE INCHI_API IXA_STEREOID INCHI_DECL IXA_MOL_GetStereoId( IXA_STATUS_HANDLE hStatus,
                                                                          IXA_MOL_HANDLE    hMolecule,
                                                                          int               vStereoIndex );

    EXPIMP_TEMPLATE INCHI_API int INCHI_DECL IXA_MOL_GetStereoIndex( IXA_STATUS_HANDLE hStatus,
                                                                    IXA_MOL_HANDLE    hMolecule,
                                                                    IXA_STEREOID      vStereo );

    EXPIMP_TEMPLATE INCHI_API IXA_STEREO_TOPOLOGY INCHI_DECL IXA_MOL_GetStereoTopology( IXA_STATUS_HANDLE hStatus,
                                                                                       IXA_MOL_HANDLE    hMolecule,
                                                                                       IXA_STEREOID      vStereo );

    EXPIMP_TEMPLATE INCHI_API IXA_ATOMID INCHI_DECL IXA_MOL_GetStereoCentralAtom( IXA_STATUS_HANDLE hStatus,
                                                                                 IXA_MOL_HANDLE    hMolecule,
                                                                                 IXA_STEREOID      vStereo );

    EXPIMP_TEMPLATE INCHI_API IXA_BONDID INCHI_DECL IXA_MOL_GetStereoCentralBond( IXA_STATUS_HANDLE hStatus,
                                                                                 IXA_MOL_HANDLE    hMolecule,
                                                                                 IXA_STEREOID      vStereo );

    EXPIMP_TEMPLATE INCHI_API int INCHI_DECL IXA_MOL_GetStereoNumVertices( IXA_STATUS_HANDLE hStatus,
                                                                          IXA_MOL_HANDLE    hMolecule,
                                                                          IXA_STEREOID      vStereo );

    EXPIMP_TEMPLATE INCHI_API IXA_ATOMID INCHI_DECL IXA_MOL_GetStereoVertex( IXA_STATUS_HANDLE hStatus,
                                                                            IXA_MOL_HANDLE    hMolecule,
                                                                            IXA_STEREOID      vStereo,
                                                                            int               vVertexIndex );

    EXPIMP_TEMPLATE INCHI_API IXA_STEREO_PARITY INCHI_DECL IXA_MOL_GetStereoParity( IXA_STATUS_HANDLE hStatus,
                                                                                      IXA_MOL_HANDLE hMolecule,
                                                                                      IXA_STEREOID   vStereo );

    /****************************************************************************/
    /*                    Functions for Generating InChIs                       */

    EXPIMP_TEMPLATE INCHI_API IXA_INCHIBUILDER_HANDLE INCHI_DECL IXA_INCHIBUILDER_Create( IXA_STATUS_HANDLE hStatus );

    EXPIMP_TEMPLATE INCHI_API void INCHI_DECL IXA_INCHIBUILDER_SetMolecule( IXA_STATUS_HANDLE       hStatus,
                                                                           IXA_INCHIBUILDER_HANDLE hInChIBuilder,
                                                                           IXA_MOL_HANDLE          hMolecule );

    EXPIMP_TEMPLATE INCHI_API const char* INCHI_DECL IXA_INCHIBUILDER_GetInChIVersion(IXA_BOOL vFullDescription);

    EXPIMP_TEMPLATE INCHI_API const char* INCHI_DECL IXA_INCHIBUILDER_GetInChI( IXA_STATUS_HANDLE       hStatus,
                                                                               IXA_INCHIBUILDER_HANDLE hInChIBuilder );

    EXPIMP_TEMPLATE INCHI_API const char* INCHI_DECL IXA_INCHIBUILDER_GetInChIEx( IXA_STATUS_HANDLE       hStatus,
                                                                                 IXA_INCHIBUILDER_HANDLE hBuilder );

    EXPIMP_TEMPLATE INCHI_API const char* INCHI_DECL IXA_INCHIBUILDER_GetAuxInfo( IXA_STATUS_HANDLE       hStatus,
                                                                                 IXA_INCHIBUILDER_HANDLE hInChIBuilder );

    EXPIMP_TEMPLATE INCHI_API const char* INCHI_DECL IXA_INCHIBUILDER_GetLog( IXA_STATUS_HANDLE       hStatus,
                                                                             IXA_INCHIBUILDER_HANDLE hInChIBuilder );

    EXPIMP_TEMPLATE INCHI_API void INCHI_DECL IXA_INCHIBUILDER_Destroy( IXA_STATUS_HANDLE       hStatus,
                                                                       IXA_INCHIBUILDER_HANDLE hInChIBuilder );


    /****************************************************************************/
    /*       Functions for Specifying/checking InChI Generation Options         */


    EXPIMP_TEMPLATE INCHI_API void INCHI_DECL IXA_INCHIBUILDER_SetOption( IXA_STATUS_HANDLE       hStatus,
                                                                         IXA_INCHIBUILDER_HANDLE hInChIBuilder,
                                                                         IXA_INCHIBUILDER_OPTION vOption,
                                                                         IXA_BOOL                vValue );

    EXPIMP_TEMPLATE INCHI_API void INCHI_DECL IXA_INCHIBUILDER_SetOption_Stereo( IXA_STATUS_HANDLE             hStatus,
                                                                                IXA_INCHIBUILDER_HANDLE       hInChIBuilder,
                                                                                IXA_INCHIBUILDER_STEREOOPTION vValue );

    EXPIMP_TEMPLATE INCHI_API void INCHI_DECL IXA_INCHIBUILDER_SetOption_Timeout( IXA_STATUS_HANDLE       hStatus,
                                                                                 IXA_INCHIBUILDER_HANDLE hInChIBuilder,
                                                                                 int                     vValue );
    EXPIMP_TEMPLATE INCHI_API void INCHI_DECL IXA_INCHIBUILDER_SetOption_Timeout_MilliSeconds( IXA_STATUS_HANDLE       hStatus,
                                                                                               IXA_INCHIBUILDER_HANDLE hInChIBuilder,
                                                                                               long                     vValue );

    EXPIMP_TEMPLATE INCHI_API IXA_BOOL INCHI_DECL IXA_INCHIBUILDER_CheckOption( IXA_STATUS_HANDLE       hStatus,
                                                                                IXA_INCHIBUILDER_HANDLE hInChIBuilder,
                                                                                IXA_INCHIBUILDER_OPTION vOption);

    EXPIMP_TEMPLATE INCHI_API IXA_BOOL INCHI_DECL IXA_INCHIBUILDER_CheckOption_Stereo(IXA_STATUS_HANDLE       hStatus,
                                                                                      IXA_INCHIBUILDER_HANDLE hInChIBuilder,
                                                                                      IXA_INCHIBUILDER_STEREOOPTION vValue);

    EXPIMP_TEMPLATE INCHI_API long INCHI_DECL IXA_INCHIBUILDER_GetOption_Timeout_MilliSeconds(IXA_STATUS_HANDLE       hStatus,
                                                                                              IXA_INCHIBUILDER_HANDLE hInChIBuilder);



    /****************************************************************************/
    /*                    Functions for Generating InChI Keys                   */


    EXPIMP_TEMPLATE INCHI_API IXA_INCHIKEYBUILDER_HANDLE INCHI_DECL IXA_INCHIKEYBUILDER_Create( IXA_STATUS_HANDLE hStatus );

    EXPIMP_TEMPLATE INCHI_API void INCHI_DECL IXA_INCHIKEYBUILDER_SetInChI( IXA_STATUS_HANDLE          hStatus,
                                                                           IXA_INCHIKEYBUILDER_HANDLE hInChIKeyBuilder,
                                                                           const char*                pInChI );

    EXPIMP_TEMPLATE INCHI_API const char* INCHI_DECL IXA_INCHIKEYBUILDER_GetInChIKey( IXA_STATUS_HANDLE          hStatus,
                                                                                     IXA_INCHIKEYBUILDER_HANDLE hInChIKeyBuilder );

    EXPIMP_TEMPLATE INCHI_API void INCHI_DECL IXA_INCHIKEYBUILDER_Destroy( IXA_STATUS_HANDLE          hStatus,
                                                                          IXA_INCHIKEYBUILDER_HANDLE hInChIKeyBuilder );


#ifdef __cplusplus
}
#endif

#endif

/* chains.cpp:
 *
 * Protein/Nucleic Acid Perception
 * Original Author:  Roger Sayle, Metaphorics
 * Original Version: Version 1.6, March 1998
 *
 * Ported & Modified By: Joe Corkery, OpenEye Scientific Software
 * Current Version:      March 2001
 */

////////////////////////////////////////////////////////////////////////////////
// File Includes
////////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include "mol.h"
#include "chains.h"
#include <map>

using namespace std;

////////////////////////////////////////////////////////////////////////////////
// Preprocessor Definitions
////////////////////////////////////////////////////////////////////////////////

#define RESIDMIN       3
#define RESIDMAX       32
#define ATOMMAX        68
#define ATOMMINAMINO   4
#define ATOMMINNUCLEIC 50
#define ELEMMAX        104
#define MAXPEPTIDE     11
#define MAXNUCLEIC     15
#define AMINOMAX       21
#define NUCLEOMAX      6
#define STACKSIZE      20
#define BUFMAX         8192
#define MAXCOVAL       2.0
#define SLOPFACTOR     0.56
#define THRESHOLD      12

#define AI_N           0
#define AI_CA          1
#define AI_C           2
#define AI_O           3
#define AI_OXT         37

#define AI_P           38
#define AI_O1P         39
#define AI_O2P         40
#define AI_O5          41
#define AI_C5          42
#define AI_C4          43
#define AI_O4          44
#define AI_C3          45
#define AI_O3          46
#define AI_C2          47
#define AI_O2          48
#define AI_C1          49

#define BitN           0x0001
#define BitNTer        0x0002
#define BitNPro        0x0004
#define BitNPT         0x0008
#define BitCA          0x0010
#define BitCAGly       0x0020
#define BitC           0x0100
#define BitCTer        0x0200
#define BitCOXT        0x0400
#define BitO           0x1000
#define BitOXT         0x2000

#define BitNAll        0x000F
#define BitCAAll       0x0030
#define BitCAll        0x0700
#define BitOAll        0x3000

#define BitP           0x0001
#define BitPTer        0x0002
#define BitOP          0x0004
#define BitO5          0x0008
#define BitO5Ter       0x0010
#define BitC5          0x0020
#define BitC4          0x0040
#define BitO4          0x0080
#define BitC3          0x0100
#define BitO3          0x0200
#define BitO3Ter       0x0400
#define BitC2RNA       0x0800
#define BitC2DNA       0x1000
#define BitO2          0x2000
#define BitC1          0x4000

#define BitPAll        0x0003
#define Bit05All       0x0018
#define BitO3All       0x0600
#define BitC2All       0x1800

#define BC_ASSIGN      0x01
#define BC_COUNT       0x02
#define BC_ELEM        0x03
#define BC_EVAL        0x04
#define BC_IDENT       0x05
#define BC_LOCAL       0x06

#define BF_SINGLE      0x01
#define BF_DOUBLE      0x02
#define BF_TRIPLE      0x04
#define BF_AROMATIC    0x08

////////////////////////////////////////////////////////////////////////////////
// Begin OpenEye Namespace
////////////////////////////////////////////////////////////////////////////////

namespace OpenBabel { 

OBChainsParser chainsparser;

////////////////////////////////////////////////////////////////////////////////
// Structure / Type Definitions
////////////////////////////////////////////////////////////////////////////////

typedef struct 
{
    char *name;
    char *data;
} ResidType;

typedef struct 
{
    int atomid,elem;
    int bcount;
    int index;
} MonoAtomType;

typedef struct 
{
    int src,dst;
    int index;
    int flag;
} MonoBondType;

typedef struct 
{
    int type;
    union _ByteCode *next;
} MonOpStruct;

typedef struct 
{
    int type;
    int value;
    union _ByteCode *tcond;
    union _ByteCode *fcond;
} BinOpStruct;

typedef struct 
{
    int type;
    int resid;
    int *atomid;
    int *bflags;
} AssignStruct;

typedef union _ByteCode 
{
    int type;
    MonOpStruct eval;     /* BC_EVAL   */
    BinOpStruct count;    /* BC_COUNT  */
    BinOpStruct elem;     /* BC_ELEM   */
    BinOpStruct ident;    /* BC_IDENT  */
    BinOpStruct local;    /* BC_LOCAL  */
    AssignStruct assign;  /* BC_ASSIGN */
} ByteCode;

typedef struct 
{
    int atom,bond;
    int prev;
} StackType;

////////////////////////////////////////////////////////////////////////////////
// Global Variables / Tables
////////////////////////////////////////////////////////////////////////////////

static char ChainsAtomName[ATOMMAX][4] = {
    /*  0 */  { ' ', 'N', ' ', ' ' },
    /*  1 */  { ' ', 'C', 'A', ' ' },
    /*  2 */  { ' ', 'C', ' ', ' ' },
    /*  3 */  { ' ', 'O', ' ', ' ' },
    /*  4 */  { ' ', 'C', 'B', ' ' },
    /*  5 */  { ' ', 'S', 'G', ' ' },
    /*  6 */  { ' ', 'O', 'G', ' ' },
    /*  7 */  { ' ', 'C', 'G', ' ' },
    /*  8 */  { ' ', 'O', 'G', '1' },
    /*  9 */  { ' ', 'C', 'G', '1' },
    /* 10 */  { ' ', 'C', 'G', '2' },
    /* 11 */  { ' ', 'C', 'D', ' ' },
    /* 12 */  { ' ', 'O', 'D', ' ' },
    /* 13 */  { ' ', 'S', 'D', ' ' },
    /* 14 */  { ' ', 'C', 'D', '1' },
    /* 15 */  { ' ', 'O', 'D', '1' },
    /* 16 */  { ' ', 'N', 'D', '1' },
    /* 17 */  { ' ', 'C', 'D', '2' },
    /* 18 */  { ' ', 'O', 'D', '2' },
    /* 19 */  { ' ', 'N', 'D', '2' },
    /* 20 */  { ' ', 'C', 'E', ' ' },
    /* 21 */  { ' ', 'N', 'E', ' ' },
    /* 22 */  { ' ', 'C', 'E', '1' },
    /* 23 */  { ' ', 'O', 'E', '1' },
    /* 24 */  { ' ', 'N', 'E', '1' },
    /* 25 */  { ' ', 'C', 'E', '2' },
    /* 26 */  { ' ', 'O', 'E', '2' },
    /* 27 */  { ' ', 'N', 'E', '2' },
    /* 28 */  { ' ', 'C', 'E', '3' },
    /* 29 */  { ' ', 'C', 'Z', ' ' },
    /* 30 */  { ' ', 'N', 'Z', ' ' },
    /* 31 */  { ' ', 'C', 'Z', '2' },
    /* 32 */  { ' ', 'C', 'Z', '3' },
    /* 33 */  { ' ', 'O', 'H', ' ' },
    /* 34 */  { ' ', 'N', 'H', '1' },
    /* 35 */  { ' ', 'N', 'H', '2' },
    /* 36 */  { ' ', 'C', 'H', '2' },
    /* 37 */  { ' ', 'O', 'X', 'T' },

    /* 38 */  { ' ', 'P', ' ', ' ' },
    /* 39 */  { ' ', 'O', '1', 'P' },
    /* 40 */  { ' ', 'O', '2', 'P' },
    /* 41 */  { ' ', 'O', '5', '*' },
    /* 42 */  { ' ', 'C', '5', '*' },
    /* 43 */  { ' ', 'C', '4', '*' },
    /* 44 */  { ' ', 'O', '4', '*' },
    /* 45 */  { ' ', 'C', '3', '*' },
    /* 46 */  { ' ', 'O', '3', '*' },
    /* 47 */  { ' ', 'C', '2', '*' },
    /* 48 */  { ' ', 'O', '2', '*' },
    /* 49 */  { ' ', 'C', '1', '*' },
    /* 50 */  { ' ', 'N', '9', ' ' },
    /* 51 */  { ' ', 'C', '8', ' ' },
    /* 52 */  { ' ', 'N', '7', ' ' },
    /* 53 */  { ' ', 'C', '5', ' ' },
    /* 54 */  { ' ', 'C', '6', ' ' },
    /* 55 */  { ' ', 'O', '6', ' ' },
    /* 56 */  { ' ', 'N', '6', ' ' },
    /* 57 */  { ' ', 'N', '1', ' ' },
    /* 58 */  { ' ', 'C', '2', ' ' },
    /* 59 */  { ' ', 'O', '2', ' ' },
    /* 60 */  { ' ', 'N', '2', ' ' },
    /* 61 */  { ' ', 'N', '3', ' ' },
    /* 62 */  { ' ', 'C', '4', ' ' },
    /* 63 */  { ' ', 'O', '4', ' ' },
    /* 64 */  { ' ', 'N', '4', ' ' },
    /* 65 */  { ' ', 'C', '5', ' ' },
    /* 66 */  { ' ', 'C', '5', 'M' },
    /* 67 */  { ' ', 'C', '6', ' ' }
        };

static Template Peptide[MAXPEPTIDE] = {
    /* N     */    {  0x0001, 7, 2, 0x0030, 0x0100,      0, 0 },
    /* NTer  */    {  0x0002, 7, 1, 0x0030,      0,      0, 0 },
    /* NPro  */    {  0x0004, 7, 3, 0x0030, 0x0100,     -6, 0 },
    /* NPT   */    {  0x0008, 7, 2, 0x0030,     -6,      0, 0 },
    /* CA    */    {  0x0010, 6, 3, 0x000F, 0x0700,     -6, 0 },
    /* CAGly */    {  0x0020, 6, 2, 0x0003, 0x0700,      0, 0 },
    /* C     */    {  0x0100, 6, 3, 0x0030, 0x1000, 0x0005, 0 },
    /* CTer  */    {  0x0200, 6, 2, 0x0030, 0x1000,      0, 0 },
    /* COXT  */    {  0x0400, 6, 3, 0x0030, 0x1000, 0x2000, 0 },
    /* O     */    {  0x1000, 8, 1, 0x0700,      0,      0, 0 },
    /* OXT   */    {  0x2000, 8, 1, 0x0400,      0,      0, 0 }
        };

static Template Nucleotide[MAXNUCLEIC] = {
    /* P     */    {  0x0001, 15, 4, 0x0004, 0x0004, 0x0008, 0x0200 },
    /* PTer  */    {  0x0002, 15, 3, 0x0004, 0x0004, 0x0008,      0 },
    /* OP    */    {  0x0004,  8, 1, 0x0003,      0,      0,      0 },
    /* O5    */    {  0x0008,  8, 2, 0x0020, 0x0003,      0,      0 },
    /* O5Ter */    {  0x0010,  8, 1, 0x0020,      0,      0,      0 },
    /* C5    */    {  0x0020,  6, 2, 0x0018, 0x0040,      0,      0 },
    /* C4    */    {  0x0040,  6, 3, 0x0020, 0x0080, 0x0100,      0 },
    /* O4    */    {  0x0080,  8, 2, 0x0040, 0x4000,      0,      0 },
    /* C3    */    {  0x0100,  6, 3, 0x0040, 0x0600, 0x1800,      0 },
    /* O3    */    {  0x0200,  8, 2, 0x0100, 0x0001,      0,      0 },
    /* O3Ter */    {  0x0400,  8, 1, 0x0100,      0,      0,      0 },
    /* C2RNA */    {  0x0800,  6, 3, 0x0100, 0x4000, 0x2000,      0 },
    /* C2DNA */    {  0x1000,  6, 2, 0x0100, 0x4000,      0,      0 },
    /* O2    */    {  0x2000,  8, 1, 0x0800,      0,      0,      0 },
    /* C1    */    {  0x4000,  6, 3, 0x0080, 0x1800,     -7,      0 }
        };

static char ChainsResName[RESIDMAX][4] = {
    /*0*/ "UNK",  /*1*/ "HOH",  /*2*/ "LIG"
};

static ResidType AminoAcids[AMINOMAX] = {
    { "ILE", "1-4(-9-14)-10"                        },
    { "VAL", "1-4(-9)-10"                           },

    { "ALA", "1-4"                                  },
    { "ASN", "1-4-7(=15)-19"                        },
    { "ASP", "1-4-7(=15)-18"                        },
    { "ARG", "1-4-7-11-21-29(=34)-35"               },
    { "CYS", "1-4-5"                                },
    { "GLN", "1-4-7-11(=23)-27"                     },
    { "GLU", "1-4-7-11(=23)-26"                     },
    { "GLY", "1"                                    },
    { "HIS", "1-4-7^16~22^27^17~7"                  },
    { "HYP", "1-4-7(-12)-11-0"                      }, /* ??? */
    { "LEU", "1-4-7(-14)-17"                        },
    { "LYS", "1-4-7-11-20-30"                       },
    { "MET", "1-4-7-13-20"                          },
    { "PHE", "1-4-7~14^22~29^25~17^7"               },
    { "PRO", "1-4-7-11-0"                           },
    { "SER", "1-4-6"                                },
    { "THR", "1-4(-8)-10"                           },
    { "TRP", "1-4-7~14^24^25~17(^7)^28~32^36~31^25" },
    { "TYR", "1-4-7~14^22~29(-33)^25~17^7"          } 
        };

/* Pyroglutamate (PCA):        1-4-7-11(=" OB ")-0  PDB Example: 1CEL */
/* Amino-N-Butyric Acid (ABA): 1-4-7                PDB Example: 1BBO */
/* Selenic Acid (SEC):         1-4-"SEG "(-15)-18   PDB Example: 1GP1 */

static ResidType Nucleotides[NUCLEOMAX] = {
    { "  A", "49-50-51-52-53-54(-56)-57-58-61-62(-53)-50"      },
    { "  C", "49-57-58(-59)-61-62(-64)-65-67-57"               },
    { "  G", "49-50-51-52-53-54(-55)-57-58(-60)-61-62(-53)-50" },
    { "  T", "49-57-58(-59)-61-62(-63)-65(-66)-67-57"          },
    { "  U", "49-57-58(-59)-61-62(-63)-65-67-57"               },
    { "  I", "49-50-51-52-53-54(-55)-57-58-61-62(-53)-50"      }
        };

static MonoAtomType MonoAtom[MaxMonoAtom];
static MonoBondType MonoBond[MaxMonoBond];
static int MonoAtomCount;
static int MonoBondCount;

static StackType Stack[STACKSIZE];
static int StackPtr;

static int  AtomIndex;
static int  BondIndex;
static bool StrictFlag = false;

////////////////////////////////////////////////////////////////////////////////
// Static Functions
////////////////////////////////////////////////////////////////////////////////

static ByteCode *AllocateByteCode(int type)
{
    ByteCode *result;

    result = (ByteCode*)malloc(sizeof(ByteCode));
    if( !result )
    {   
      cerr << "Error: Unable to allocate byte codes in chains.cpp!" << endl;
      exit(1);
    }
    result->type = type;
    
	return (result);
}

static void FatalMemoryError(void)
{
  cerr << "Error: Fatal memory allocation error in chains.cpp!" << endl;
  exit(1);
}

void GenerateByteCodes(ByteCode **node, int resid, int curr, int prev, int bond)
{
	StackType neighbour[4];
	StackType original;
	int count,i,j;
	ByteCode *ptr;
	bool done,found;

	if( curr != prev )
	{   
		if( MonoAtom[curr].atomid < ATOMMINAMINO )
		{   
			found = false;
			while( *node && ((*node)->type==BC_IDENT) )
			{
				if( (*node)->ident.value == MonoAtom[curr].atomid )
				{   
					node  = (ByteCode**)&(*node)->ident.tcond;
					found = true;
					break;
				} 
				else 
					node = (ByteCode**)&(*node)->ident.fcond;
			}

			if (!found)
			{   
				ptr = AllocateByteCode(BC_IDENT);
				ptr->ident.tcond = (ByteCode*)0;
				ptr->ident.fcond = *node; *node = ptr;
				node = (ByteCode**)&ptr->ident.tcond;
				ptr->ident.value = MonoAtom[curr].atomid;
			}
			MonoBond[bond].index = BondIndex++;
			done = true;
		} 
		else if( MonoAtom[curr].index != -1 )
		{   
			while( *node && ((*node)->type==BC_IDENT) )
				node = (ByteCode**)&(*node)->ident.fcond;

			found = false;
			while( *node && ((*node)->type==BC_LOCAL) )
			{
				if( (*node)->local.value == MonoAtom[curr].index )
				{   
					node = (ByteCode**)&(*node)->local.tcond;
					found = true;
					break;
				} 
				else 
					node = (ByteCode**)&(*node)->local.fcond;
			}

			if (!found)
			{   
				ptr = AllocateByteCode(BC_LOCAL);
				ptr->local.tcond = (ByteCode*)0;
				ptr->local.fcond = *node; *node = ptr;
				node = (ByteCode**)&ptr->local.tcond;
				ptr->local.value = MonoAtom[curr].index;
			}

			MonoBond[bond].index = BondIndex++;
			done = true;
		} 
		else 
		{   
			while( *node && ((*node)->type==BC_IDENT) )
				node = (ByteCode**)&(*node)->ident.fcond;
			while( *node && ((*node)->type==BC_LOCAL) )
				node = (ByteCode**)&(*node)->local.fcond;

			found = false;
			while( *node && ((*node)->type==BC_ELEM) )
			{
				if( (*node)->elem.value == MonoAtom[curr].elem )
				{   
					node = (ByteCode**)&(*node)->elem.tcond;
					found = true;
					break;
				} 
				else 
					node = (ByteCode**)&(*node)->elem.fcond;
			}

			if( !found )
			{   
				ptr = AllocateByteCode(BC_ELEM);
				ptr->elem.tcond = (ByteCode*)0;
				ptr->elem.fcond = *node; *node = ptr;
				node = (ByteCode**)&ptr->elem.tcond;
				ptr->elem.value = MonoAtom[curr].elem;
			}

			MonoAtom[curr].index = AtomIndex++;
			MonoBond[bond].index = BondIndex++;
			done = false;
		}
	} 
	else 
	{   
		MonoAtom[curr].index = AtomIndex++;
		done = false;
	}

	count = 0;
	if (!done)
	{   
		for( i=0; i<MonoBondCount; i++ )
		{
			if( MonoBond[i].src == curr )
			{   
				if( MonoBond[i].dst != prev )
				{   
					neighbour[count].atom = MonoBond[i].dst;
					neighbour[count].bond = i;
					count++;
				}
			} 
			else if( MonoBond[i].dst == curr )
			{   
				if( MonoBond[i].src != prev )
				{   
					neighbour[count].atom = MonoBond[i].src;
					neighbour[count].bond = i;
					count++;
				}
			}
		}

		if ( *node && ((*node)->type==BC_EVAL) )
		{   
			found = false;
			node  = (ByteCode**)&(*node)->eval.next;
			while( *node && ((*node)->type==BC_COUNT) )
			{
				if( (*node)->count.value == count )
				{   
					node = (ByteCode**)&(*node)->count.tcond;
					found = true;
					break;
				} 
				else 
					node = (ByteCode**)&(*node)->count.fcond;
			}

			if( !found )
			{   
				ptr = AllocateByteCode(BC_COUNT);
				ptr->count.tcond = (ByteCode*)0;
				ptr->count.fcond = *node; *node = ptr;
				node = (ByteCode**)&ptr->count.tcond;
				ptr->count.value = count;
			}
		} 
		else if( count || StrictFlag || StackPtr )
		{   
			ptr = AllocateByteCode(BC_EVAL);
			ptr->eval.next = *node;  *node = ptr;
			node = (ByteCode**)&ptr->eval.next;

			ptr = AllocateByteCode(BC_COUNT);
			ptr->count.tcond = (ByteCode*)0;
			ptr->count.fcond = *node; *node = ptr;
			node = (ByteCode**)&ptr->count.tcond;
			ptr->count.value = count;
		}
	}

	if( count == 1 )
	{   
		GenerateByteCodes(node,resid,neighbour[0].atom, curr,neighbour[0].bond);
	} 
	else if( count == 2 )
	{   
		original = Stack[StackPtr++];
		Stack[StackPtr-1] = neighbour[0];
		Stack[StackPtr-1].prev = curr;
		GenerateByteCodes(node,resid,neighbour[1].atom,
								curr,neighbour[1].bond);
		Stack[StackPtr-1] = neighbour[1];
		Stack[StackPtr-1].prev = curr;
		GenerateByteCodes(node,resid,neighbour[0].atom,
								curr,neighbour[0].bond);
		Stack[--StackPtr] = original;
	} 
	else if( count )
	{   
	  cerr << "Error: Maximum Monomer Fanout Exceeded!" << endl;
		fprintf(stderr,"Residue %s atom %d\n",ChainsResName[resid],curr);
		fprintf(stderr,"Previous = %d  Fanout = %d\n",prev,count);
		exit(1);
	} 
	else if( StackPtr )
	{   
		StackPtr--;
		GenerateByteCodes(node,resid,Stack[StackPtr].atom,
				Stack[StackPtr].prev,Stack[StackPtr].bond);
		StackPtr++;
	} 
	else if( !(*node) )
	{   
		ptr = AllocateByteCode(BC_ASSIGN);
		ptr->assign.resid = resid;
		ptr->assign.atomid = (int*)malloc(AtomIndex*sizeof(int));
		if( !ptr->assign.atomid ) FatalMemoryError();
		for( i=0; i<MonoAtomCount; i++ )
			if( (j=MonoAtom[i].index) != -1 )
			   ptr->assign.atomid[j] = MonoAtom[i].atomid;
		if( BondIndex )
		{   
			ptr->assign.bflags = (int*)malloc(BondIndex*sizeof(int));
			for( i=0; i<MonoBondCount; i++ )
				if( (j=MonoBond[i].index) != -1 )
					ptr->assign.bflags[j] = MonoBond[i].flag;
		}
		*node = ptr;
	} 
	else if( (*node)->type == BC_ASSIGN )
	{   
		if( (*node)->assign.resid != resid )
		{   
			fputs("Error: Duplicated Monomer Specification!\n",stderr);
			fprintf(stderr,"Residue %s matches resid",ChainsResName[resid]);
			fprintf(stderr,"ue %s!\n",ChainsResName[(*node)->assign.resid]);
		}
	}

	/* Restore State! */
	if( curr != prev )
	{   
		if( !done )
		{   
			MonoAtom[curr].index = -1;
			AtomIndex--;
		}
		MonoBond[bond].index = -1;
		BondIndex--;
	}
}
    
////////////////////////////////////////////////////////////////////////////////
// Constructors / Destructors
////////////////////////////////////////////////////////////////////////////////

// validated
OBChainsParser::OBChainsParser(void)
{
	int i, res = RESIDMIN;

	PDecisionTree = (ByteCode*)0;
	for( i=0 ; i < AMINOMAX ; i++ )
	{   
		strcpy(ChainsResName[res],AminoAcids[i].name);
		DefineMonomer(&PDecisionTree,res,AminoAcids[i].data);
		res++;
	}

	NDecisionTree = (ByteCode*)0;
	for( i=0 ; i< NUCLEOMAX ; i++ )
	{    
		strcpy(ChainsResName[res],Nucleotides[i].name);
		DefineMonomer(&NDecisionTree,res,Nucleotides[i].data);
		res++;
	}

	bitmasks = NULL;
	hetflags = NULL;
	atomids  = NULL;
	resids   = NULL;
	resnos   = NULL;
	sernos   = NULL;
	hcounts  = NULL;
	chains   = NULL;
	flags    = NULL;
}

OBChainsParser::~OBChainsParser(void)
{
}

////////////////////////////////////////////////////////////////////////////////
// Setup / Cleanup Functions
////////////////////////////////////////////////////////////////////////////////

void OBChainsParser::SetupMol(OBMol &mol)
{
	CleanupMol();

	int i;
	int asize = mol.NumAtoms();
	int bsize = mol.NumBonds();

	bitmasks = new unsigned short[asize];
	resids   = new unsigned char[asize];
	flags    = new unsigned char[bsize];
	hetflags = new bool[asize];
	atomids  = new short[asize];
	resnos   = new short[asize];
	sernos   = new short[asize];
	hcounts  = new char[asize];
	chains   = new char[asize];

	for ( i = 0 ; i < asize ; i++ )
	{
		hetflags[i] = false;
		bitmasks[i] = 0;
		atomids[i]  = -1;
		resids[i]   = 0;
		resnos[i]   = 0;
		sernos[i]   = 0;
		hcounts[i]  = 0;
		chains[i]   = ' ';
	}

	for ( i = 0 ; i < bsize ; i++ )
		flags[i] = 0;
}

void OBChainsParser::CleanupMol(void)
{
  if (bitmasks != NULL) {delete bitmasks; bitmasks = NULL;}
  if (hetflags != NULL) {delete hetflags; hetflags = NULL;}
  if (atomids  != NULL) {delete atomids; atomids = NULL;}
  if (resids   != NULL) {delete resids; resids = NULL;}
  if (resnos   != NULL) {delete resnos; resnos = NULL;}
  if (sernos   != NULL) {delete sernos; sernos = NULL;}
  if (hcounts  != NULL) {delete hcounts; hcounts = NULL;}
  if (chains   != NULL) {delete chains; chains = NULL;}
  if (flags    != NULL) {delete flags; flags = NULL;}
}

void OBChainsParser::ClearResidueInformation(OBMol &mol)
{
	OBResidue *residue;
	vector<OBResidue*> residues;
	vector<OBResidue*>::iterator r;

	for (residue = mol.BeginResidue(r) ; residue ; residue = mol.NextResidue(r))
		residues.push_back(residue);

	for ( unsigned int i = 0 ; i < residues.size() ; i++ )
		mol.DeleteResidue(residues[i]);

	residues.clear();
}

void OBChainsParser::SetResidueInformation(OBMol &mol)
{
    char buffer[256];
    string atomid, name;

    OBAtom    *atom;
    OBResidue *residue;
	map<short, OBResidue *> resmap;

    int size = mol.NumAtoms();
    for ( int i = 0 ; i < size ; i++ )
    {
        atom = mol.GetAtom(i+1); // WARNING: ATOM INDEX ISSUE

        if (atomids[i] == -1)
            sprintf(buffer, "%s", etab.GetSymbol(atom->GetAtomicNum()));
        else if (atom->IsHydrogen())
        {
            if (hcounts[i])
                sprintf(buffer, "%cH%.2s", hcounts[i]+'0', ChainsAtomName[atomids[i]]+2);
            else 
                sprintf(buffer, "H%.2s", ChainsAtomName[atomids[i]]+2);
        }
        else
            sprintf(buffer, "%.4s", ChainsAtomName[atomids[i]]);

        if (buffer[3] == ' ')
            buffer[3] = '\0';

        atomid = (buffer[0] == ' ') ? buffer + 1 : buffer;

        if (resmap.find(resnos[i]) != resmap.end())
        {
            residue = resmap[resnos[i]];
            residue->AddAtom(atom);
            residue->SetAtomID(atom, atomid);
            residue->SetHetAtom(atom, hetflags[i]);
            residue->SetSerialNum(atom, sernos[i]);
        }
        else
        {
            name    = ChainsResName[resids[i]];
            residue = mol.NewResidue();

            residue->SetName(name);
            residue->SetNum(resnos[i]);
            residue->SetChain(chains[i]);
            residue->SetChainNum((chains[i] > 'A') ? (int)(chains[i] - 'A') : 1);

            residue->AddAtom(atom);
            residue->SetAtomID(atom, atomid);
            residue->SetHetAtom(atom, hetflags[i]);
            residue->SetSerialNum(atom, sernos[i]);

            resmap[resnos[i]] = residue;
        }
    }

    if (mol.NumResidues() == 1)
        mol.DeleteResidue(mol.GetResidue(0));
}

////////////////////////////////////////////////////////////////////////////////
// Perception Functions
////////////////////////////////////////////////////////////////////////////////

bool OBChainsParser::PerceiveChains(OBMol &mol)
{
	bool result = true;

	SetupMol(mol);
	ClearResidueInformation(mol);

	result = DetermineHetAtoms(mol)          && result;
	result = DetermineConnectedChains(mol)   && result;
	result = DeterminePeptideBackbone(mol)   && result;
	result = DeterminePeptideSidechains(mol) && result;
	result = DetermineNucleicBackbone(mol)   && result;
	result = DetermineNucleicSidechains(mol) && result;
	result = DetermineHydrogens(mol)         && result;

	SetResidueInformation(mol);
	CleanupMol();

	return result;
}

////////////////////////////////////////////////////////////////////////////////
// Hetero Atom Perception
////////////////////////////////////////////////////////////////////////////////

bool OBChainsParser::DetermineHetAtoms(OBMol &mol)
{
	OBAtom *atom;
	vector<OBNodeBase *>::iterator a;
	for (atom = mol.BeginAtom(a) ; atom ; atom = mol.NextAtom(a))
		if (!atom->IsHydrogen() && atom->GetValence() == 0)
		{
			resids[atom->GetIdx()-1]   = (atom->IsOxygen()) ? 1 : 2;
			hetflags[atom->GetIdx()-1] = true;
		}
	return true;
}

////////////////////////////////////////////////////////////////////////////////
// Connected Chain Perception
////////////////////////////////////////////////////////////////////////////////

bool OBChainsParser::DetermineConnectedChains(OBMol &mol)
{
	int resid;
	int resno;
    int count;
    int size;
    int i,idx;
	int numAtoms;

	resno    = 1;
    count    = 0;
	numAtoms = mol.NumAtoms();

	OBAtom *atom;
	vector<OBNodeBase *>::iterator a;
	for (atom = mol.BeginAtom(a) ; atom ; atom = mol.NextAtom(a))
	{
		idx = atom->GetIdx() - 1;
		if (!hetflags[idx] && chains[idx] == ' ' && !atom->IsHydrogen())
		{
			size = RecurseChain(mol, idx, 'A' + count);
			if (size < 10)
			{
				if (size == 1 && atom->IsOxygen())
					resid = 1; /* HOH */
				else
					resid = 2;
				
				for (i = 0 ; i < numAtoms ; i++)
				{
					if (chains[i] == ('A' + count))
					{
						hetflags[i] = true;
						resids[i]   = resid;
						resnos[i]   = resno;
						chains[i]   = ' ';
					}
				}
				resno++;
            } 
			else 
				count++;
		}
	}

    if( count == 1 )
		for ( i = 0 ; i < numAtoms ; i++ )
			chains[i] = ' ';

	return true;
}

int OBChainsParser::RecurseChain(OBMol &mol, int i, int c)
{
	OBAtom *atom, *nbr;
	vector<OBEdgeBase *>::iterator b;
    int result;

	atom      = mol.GetAtom(i+1);
    result    = (atom->IsHydrogen()) ? 0 : 1;
	chains[i] = c;

	for (nbr = atom->BeginNbrAtom(b) ; nbr ; nbr = atom->NextNbrAtom(b))
		if (chains[nbr->GetIdx()-1] == ' ')
			result += RecurseChain(mol,nbr->GetIdx()-1,c);

    return (result);
}

////////////////////////////////////////////////////////////////////////////////
// Peptide Backbone Perception
////////////////////////////////////////////////////////////////////////////////

bool OBChainsParser::DeterminePeptideBackbone(OBMol &mol)
{
    ConstrainBackbone(mol, Peptide, MAXPEPTIDE);
        
	int i, max = mol.NumAtoms(); 

    /*
    int count = 0;
    for ( i = 0 ; i < max ; i++ )
        if ( bitmasks[i]&BitCAAll ) 
			count++;

    fprintf(stderr,"%d alpha carbons\n",count);
    */

    /* Order Peptide Backbone */

	for ( i = 0 ; i < max ; i++ )
		if (atomids[i] == -1)
        {   
			if( bitmasks[i] & BitNTer )
            {   
				atomids[i] = AI_N;
                TracePeptideChain(mol,i,1);
            } 
			else if( (bitmasks[i]&BitNPT) && !(bitmasks[i]&BitN) )
            {   
				atomids[i] = AI_N;
                TracePeptideChain(mol,i,1);
            }
        }

    /* Carbonyl Double Bond */

	OBBond *bond;
	vector<OBEdgeBase*>::iterator b;
	for (bond = mol.BeginBond(b) ; bond ; bond = mol.NextBond(b))
	{
		if ((atomids[bond->GetBeginAtomIdx()-1] == 2 && atomids[bond->GetEndAtomIdx()-1] == 3) ||
			(atomids[bond->GetBeginAtomIdx()-1] == 3 && atomids[bond->GetEndAtomIdx()-1] == 2))
			flags[bond->GetIdx()] |= BF_DOUBLE;
	}	

	return true;
}

void OBChainsParser::ConstrainBackbone(OBMol &mol, Template *templ, int tmax)
{
    static OBAtom *neighbour[6];
    Template *pep;
    OBAtom *na,*nb,*nc,*nd;
    OBAtom *atom, *nbr;
    bool change, result;
    int  count;
    int  i,idx;

	vector<OBNodeBase *>::iterator a;
	vector<OBEdgeBase *>::iterator b;	

    /* First Pass */

	for (atom = mol.BeginAtom(a) ; atom ; atom = mol.NextAtom(a))
	{
		idx = atom->GetIdx() - 1;
		bitmasks[idx] = 0;
		for ( i = 0 ; i < tmax ; i++ )
			if ( (static_cast<unsigned int>(templ[i].elem)  == atom->GetAtomicNum()) &&
				 (static_cast<unsigned int>(templ[i].count) == atom->GetValence()))
				bitmasks[idx] |= templ[i].flag;
	}

    /* Second Pass */

    do 
	{
        change = false;
		for (atom = mol.BeginAtom(a) ; atom ; atom = mol.NextAtom(a))
		{
			idx = atom->GetIdx() - 1;
			if (bitmasks[idx]) // Determine Neighbours
			{
				count = 0;
				for (nbr = atom->BeginNbrAtom(b) ; nbr ; nbr = atom->NextNbrAtom(b))
					if (!nbr->IsHydrogen())
						neighbour[count++] = nbr;

				na = neighbour[0];
				nb = neighbour[1];
				nc = neighbour[2];
				nd = neighbour[3];

				for ( i = 0 ; i < tmax ; i++ )
					if ( templ[i].flag & bitmasks[idx] )
					{
						pep    = &templ[i];
						result = true;

						if (count == 4)
							result = Match4Constraints(pep,na,nb,nc,nd);
						else if (count == 3)
							result = Match3Constraints(pep,na,nb,nc);
						else if (count == 2)
							result = Match2Constraints(pep,na,nb);
						else // count == 1
							result = MatchConstraint(na,pep->n1);
		
						if(result == false)
						{
							bitmasks[idx] &= ~pep->flag;
							change = true;
						}
					}
			}
		}
    } while( change );
}

bool OBChainsParser::MatchConstraint(OBAtom *atom, int mask)
{
    if( mask < 0 )
		return(atom->GetAtomicNum() == static_cast<unsigned int>(-mask));
    else 
		return(((bitmasks[atom->GetIdx()-1]&mask) == 0) ? false : true);
}

bool OBChainsParser::Match2Constraints(Template *tmpl, OBAtom *na, OBAtom *nb)
{
    if( MatchConstraint(na,tmpl->n2) )
        if( MatchConstraint(nb,tmpl->n1) )
            return( true );
    if( MatchConstraint(nb,tmpl->n2) )
        if( MatchConstraint(na,tmpl->n1) )
            return( true );
    return( false );
}

bool OBChainsParser::Match3Constraints(Template *tmpl, OBAtom *na, OBAtom *nb, OBAtom *nc)
{
    if( MatchConstraint(na,tmpl->n3) )
        if( Match2Constraints(tmpl,nb,nc) )
            return( true );
    if( MatchConstraint(nb,tmpl->n3) )
        if( Match2Constraints(tmpl,na,nc) )
            return( true );
    if( MatchConstraint(nc,tmpl->n3) )
        if( Match2Constraints(tmpl,na,nb) )
            return( true );
    return( false );
}

bool OBChainsParser::Match4Constraints(Template *tmpl, OBAtom *na, OBAtom *nb, OBAtom *nc, OBAtom *nd)
{
    if( MatchConstraint(na,tmpl->n4) )
        if( Match3Constraints(tmpl,nb,nc,nd) )
            return( true );
    if( MatchConstraint(nb,tmpl->n4) )
        if( Match3Constraints(tmpl,na,nc,nd) )
            return( true );
    if( MatchConstraint(nc,tmpl->n4) )
        if( Match3Constraints(tmpl,na,nb,nd) )
            return( true );
    if( MatchConstraint(nd,tmpl->n4) )
        if( Match3Constraints(tmpl,na,nb,nc) )
            return( true );
    return( false );
}

void OBChainsParser::TracePeptideChain(OBMol &mol, int i, int r)
{
    int neighbour[4];
    int na,nb,nc;
    OBAtom *atom, *nbr;
    int count;
    int j,k,idx;

	vector<OBEdgeBase *>::iterator b;

    /* Determine Neighbours */

	atom = mol.GetAtom(i+1);
	idx  = atom->GetIdx() - 1;

    count = 0;
	for (nbr = atom->BeginNbrAtom(b) ; nbr ; nbr = atom->NextNbrAtom(b))
		if (!nbr->IsHydrogen())
			neighbour[count++] = nbr->GetIdx()-1;

	resnos[idx] = r;

    na = neighbour[0];
    nb = neighbour[1];
    nc = neighbour[2];

    switch( atomids[i] )
    {   
		case(AI_N):   
			for( j=0; j<count; j++ )
				if( bitmasks[neighbour[j]] & BitCAAll )
				{   
					atomids[neighbour[j]] = AI_CA;
					TracePeptideChain(mol,neighbour[j],r);
				}
			break;

        case(AI_CA):  
			if( count == 3 )
			{   
				if ( bitmasks[na] & BitNAll )
					na = nc;
				else if ( bitmasks[nb] & BitNAll )
					nb = nc;

				if ( bitmasks[na] & BitC )
				{   
					j = na;  
					k = nb;
				} 
				else if ( bitmasks[nb] & BitC )
				{   
					j = nb;
					k = na;
				} 
				else if( bitmasks[na] & BitCAll )
				{   
					j = na;
					k = nb;
				} 
				else /* bitmasks[nb] & BitCAll */
				{   
					j = nb;
					k = na;
				}

				atomids[j]  = AI_C;
				bitmasks[k] = 0;

				TracePeptideChain(mol,j,r);
			} 
			else /* count == 2 */
			{   
				if ( bitmasks[na] & BitCAll )
				{   
					atomids[na] = AI_C;
					TracePeptideChain(mol,na,r);
				} 
				else 
				{   
					atomids[nb] = AI_C;
					TracePeptideChain(mol,nb,r);
				}
			}
			break;

        case(AI_C):  
			k = AI_O;
			for ( j = 0; j < count; j++ )
			{
				if ( bitmasks[neighbour[j]] & BitNAll )
				{   
					atomids[neighbour[j]] = AI_N;
					TracePeptideChain(mol,neighbour[j],r+1);
				} 
				else if( bitmasks[neighbour[j]] & BitOAll )
				{   
					atomids[neighbour[j]] = k;
					resnos[neighbour[j]]  = r;
					k = AI_OXT;  /* OXT */
				}
			}
			break;
    }
}

////////////////////////////////////////////////////////////////////////////////
// Peptide Sidechains Perception
////////////////////////////////////////////////////////////////////////////////

bool OBChainsParser::DeterminePeptideSidechains(OBMol &mol)
{
	int resid;
	int max = mol.NumAtoms();

	for (int i = 0 ; i < max ; i++)
		if (atomids[i] == 1)
		{
			resid = IdentifyResidue(PDecisionTree, mol, i, resnos[i]);
			AssignResidue(mol,resnos[i],chains[i],resid);
		}

	return true;
}

void OBChainsParser::AssignResidue(OBMol &mol, int r, int c, int i)
{
	int max = mol.NumAtoms();
	for (int j = 0 ; j < max ; j++)
		if ((resnos[j] == r) && (chains[j] == c) && !hetflags[j])
			resids[j] = i;
}

int OBChainsParser::IdentifyResidue(void *tree, OBMol &mol, int seed, int resno)
{
    ByteCode *ptr;

    int AtomCount, BondCount;
    int curr,prev,bond;
    int bcount;
    int i,j;

    ptr    = (ByteCode *) tree;
    bcount = 0;

    Stack[0].atom = seed;
    Stack[0].prev = seed;
    StackPtr = 0;

    ResMonoAtom[0] = seed;
    AtomCount = 1;
    BondCount = 0;

    OBAtom *atom, *nbr;
    vector<OBEdgeBase *>::iterator b;

    while( ptr )
        switch(ptr->type)
        {   case(BC_IDENT):  curr = Stack[StackPtr-1].atom;
                             if( atomids[curr] == ptr->ident.value )
                             {   bond = Stack[StackPtr-1].bond;
                                 ResMonoBond[BondCount++] = bond;
                                 ptr = ptr->ident.tcond;
                                 StackPtr--;
                             } else ptr = ptr->ident.fcond;
                             break;

            case(BC_LOCAL):  curr = Stack[StackPtr-1].atom;
                             if( curr == ResMonoAtom[ptr->local.value] )
                             {   bond = Stack[StackPtr-1].bond;
                                 ResMonoBond[BondCount++] = bond;
                                 ptr = ptr->local.tcond;
                                 StackPtr--;
                             } else ptr = ptr->local.fcond;
                             break;

            case(BC_ELEM):   curr = Stack[StackPtr-1].atom;
                             if( mol.GetAtom(curr+1)->GetAtomicNum() == static_cast<unsigned int>(ptr->elem.value) )
                             {   bond = Stack[StackPtr-1].bond;
                                 ResMonoAtom[AtomCount++] = curr;
                                 ResMonoBond[BondCount++] = bond;
                                 resnos[curr] = resno;
                                 ptr = ptr->elem.tcond;
                                 StackPtr--;
                             } else ptr = ptr->elem.fcond;
                             break;

            case(BC_EVAL):   bcount = 0;
                             curr = Stack[StackPtr].atom;
                             prev = Stack[StackPtr].prev;

                             atom = mol.GetAtom(curr+1);
                             for (nbr = atom->BeginNbrAtom(b) ; nbr ; nbr = atom->NextNbrAtom(b))
                             {
                                 j = nbr->GetIdx() - 1;
                                 if (!((curr == prev) && bitmasks[j]) && (j != prev) && !(nbr->IsHydrogen()))
                                 {
									 Stack[StackPtr].prev = curr;
                                     Stack[StackPtr].atom = j;
                                     Stack[StackPtr].bond = (*b)->GetIdx();
                                     StackPtr++;
                                     bcount++;
                                 }
                             }

                             ptr = ptr->eval.next;
                             break;

            case(BC_COUNT):  if( bcount == ptr->count.value )
                             {      ptr = ptr->count.tcond;
                             } else ptr = ptr->count.fcond;
                             break;

            case(BC_ASSIGN): for( i=0; i<AtomCount; i++ )
                                 if( !bitmasks[ResMonoAtom[i]] )
                                 {   j = ptr->assign.atomid[i];
                                     atomids[ResMonoAtom[i]] = j;
                                 }
                             for( i=0; i<BondCount; i++ )
                             {   j = ptr->assign.bflags[i];
                                 flags[ResMonoBond[i]] = j;
                             }
                             return( ptr->assign.resid );

            default:  /* Illegal Instruction! */
                      return( 0 );
        }
    return 0;
}

////////////////////////////////////////////////////////////////////////////////
// Nucleic Backbone Perception
////////////////////////////////////////////////////////////////////////////////

bool OBChainsParser::DetermineNucleicBackbone(OBMol &mol)
{
    ConstrainBackbone(mol, Nucleotide, MAXNUCLEIC);

	int i, max = mol.NumAtoms();    

    /*
    int count = 0;
	for ( i = 0 ; i < max ; i++ )
        if ( bitmasks[i] & BitC5 ) 
			count++;

    fprintf(stderr,"%d sugar phosphates\n",count);
    */

    /* Order Nucleic Backbone */

    for( i = 0 ; i < max ; i++ )
		if( atomids[i] == -1 )
        {   
			if( bitmasks[i] & BitPTer )
            {   
				atomids[i] = AI_P;
                TraceNucleicChain(mol,i,1);
            } 
			else if( bitmasks[i] & BitO5Ter )
            {   
				atomids[i] = AI_O5;
                TraceNucleicChain(mol,i,1);
            }
        }

	return true;
}

void OBChainsParser::TraceNucleicChain(OBMol &mol, int i, int r)
{
    int neighbour[4];
    int na,nb,nc;
    int count;
    int j,k;

	OBAtom *atom, *nbr;
	vector<OBEdgeBase *>::iterator b;

	count = 0;
	atom  = mol.GetAtom(i + 1);	
	for (nbr = atom->BeginNbrAtom(b) ; nbr ; nbr = atom->NextNbrAtom(b))
		if (!nbr->IsHydrogen())
			neighbour[count++] = nbr->GetIdx() - 1;

	resnos[i] = r;

    na = neighbour[0];
    nb = neighbour[1];
    nc = neighbour[2];

    switch( atomids[i] )
    {   
		case(AI_P): 
			k = AI_O1P;  /* O1P */
            for( j=0; j<count; j++ )
			{
				if( bitmasks[neighbour[j]] & BitO5 )
				{   
					atomids[neighbour[j]] = AI_O5;
					TraceNucleicChain(mol,neighbour[j],r);
				} 
				else if( bitmasks[neighbour[j]] & BitOP )
				{   
					atomids[neighbour[j]] = k;
					resnos[neighbour[j]]  = r;
					k = AI_O2P;  /* O2P */
				}
			}

			break;

        case(AI_O5):  
			for( j=0; j<count; j++ )
				if( bitmasks[neighbour[j]] & BitC5 )
				{   
					atomids[neighbour[j]] = AI_C5;
					TraceNucleicChain(mol,neighbour[j],r);
				}
			
			break;

        case(AI_C5):  
			for( j=0 ; j<count; j++ )
				if( bitmasks[neighbour[j]] & BitC4 )
				{   
					atomids[neighbour[j]] = AI_C4;
					TraceNucleicChain(mol,neighbour[j],r);
				}
			
			break;

        case(AI_C4):  
			for( j=0; j<count; j++ )
			{
				if( bitmasks[neighbour[j]] & BitC3 )
				{   
					atomids[neighbour[j]] = AI_C3;
					TraceNucleicChain(mol,neighbour[j],r);		
				} 
				else if( bitmasks[neighbour[j]] & BitO4 )
				{   
					atomids[neighbour[j]] = AI_O4;
					resnos[neighbour[j]]  = r;
				}
			}		

			break;

        case(AI_C3):  
			for( j=0; j<count; j++ )
			{
				if( bitmasks[neighbour[j]] & BitO3All )
				{   
					atomids[neighbour[j]] = AI_O3;
					TraceNucleicChain(mol,neighbour[j],r);
				} 
				else if( bitmasks[neighbour[j]] & BitC2All )
				{   
					atomids[neighbour[j]] = AI_C2;
					TraceNucleicChain(mol,neighbour[j],r);
				}
			}

			break;

        case(AI_O3):  
			for( j=0; j<count; j++ )
				if( bitmasks[neighbour[j]] & BitP )
				{   
					atomids[neighbour[j]] = AI_P;
					TraceNucleicChain(mol,neighbour[j],r+1);
				}

			break;

        case(AI_C2):
			for( j=0; j<count; j++ )
			{
				if( bitmasks[neighbour[j]] & BitC1 )
				{   
					atomids[neighbour[j]] = AI_C1;
					resnos[neighbour[j]]  = r;
				} 
				else if( bitmasks[neighbour[j]] & BitO2 )
				{   
					atomids[neighbour[j]] = AI_O2;
					resnos[neighbour[j]]  = r;
				}
			}

			break;
    }
}

////////////////////////////////////////////////////////////////////////////////
// Nucleic Sidechains Perception
////////////////////////////////////////////////////////////////////////////////

bool OBChainsParser::DetermineNucleicSidechains(OBMol &mol)
{
    for( unsigned int i = 0 ; i < mol.NumAtoms() ; i++ )
		if( atomids[i] == 49 )
        {   
			int resid = IdentifyResidue(NDecisionTree,mol,i,resnos[i]);
            AssignResidue(mol,resnos[i],chains[i],resid);
        } 

	return true;
}

////////////////////////////////////////////////////////////////////////////////
// Hydrogens Perception
////////////////////////////////////////////////////////////////////////////////

bool OBChainsParser::DetermineHydrogens(OBMol &mol)
{
    OBAtom *atom, *nbr;
    int idx,sidx;

	int max = mol.NumAtoms();
	for ( int i = 0 ; i < max ; i++ )
		hcounts[i] = 0;

    /* First Pass */

	vector<OBNodeBase*>::iterator a;
	vector<OBEdgeBase*>::iterator b;

    for(atom = mol.BeginAtom(a); atom ; atom = mol.NextAtom(a))
        if(atom->IsHydrogen())
        {   
			nbr = atom->BeginNbrAtom(b);
			if (nbr != NULL)
			{
				idx  = atom->GetIdx() - 1;
				sidx = nbr->GetIdx() - 1;

				hcounts[idx]  = ++hcounts[sidx];
				hetflags[idx] = hetflags[sidx];
				atomids[idx]  = atomids[sidx];
				resids[idx]   = resids[sidx];
				resnos[idx]   = resnos[sidx];
			}
        }

    /* Second Pass */

	for(atom = mol.BeginAtom(a) ; atom ; atom = mol.NextAtom(a))
		if (atom->IsHydrogen())
		{
			nbr = atom->BeginNbrAtom(b);
			if (nbr != NULL && hcounts[nbr->GetIdx()-1] == 1)
				hcounts[atom->GetIdx()-1] = 0;	
		}

	return true;
}

////////////////////////////////////////////////////////////////////////////////
// Utility Functions
////////////////////////////////////////////////////////////////////////////////

// validated
void OBChainsParser::DefineMonomer(void **tree, int resid, char *smiles)
{
    int i;

    MonoAtomCount = 0;
    MonoBondCount = 0;

    ParseSmiles(smiles,-1);

    for( i=0; i<MonoBondCount; i++ )
        MonoBond[i].index = -1;
    for( i=0; i<MonoAtomCount; i++ )
        MonoAtom[i].index = -1;
    AtomIndex = BondIndex = 0;

    StackPtr = 0;
    GenerateByteCodes((ByteCode**)tree, resid, 0, 0, 0 );
}

int OBChainsParser::IdentifyElement(char *ptr)
{
    int ch;

    ch = toupper(ptr[1]);
    switch( toupper(ptr[0]) )
    {   case(' '):  switch( ch )
                    {   case('B'):  return(  5 );
                        case('C'):  return(  6 );
                        case('D'):  return(  1 );
                        case('F'):  return(  9 );
                        case('H'):  return(  1 );
                        case('I'):  return( 53 );
                        case('K'):  return( 19 );
                        case('L'):  return(  1 );
                        case('N'):  return(  7 );
                        case('O'):  return(  8 );
                        case('P'):  return( 15 );
                        case('S'):  return( 16 );
                        case('U'):  return( 92 );
                        case('V'):  return( 23 );
                        case('W'):  return( 74 );
                        case('Y'):  return( 39 );
                    }
                    break;

        case('A'):  switch( ch )
                    {   case('C'):  return( 89 );
                        case('G'):  return( 47 );
                        case('L'):  return( 13 );
                        case('M'):  return( 95 );
                        case('R'):  return( 18 );
                        case('S'):  return( 33 );
                        case('T'):  return( 85 );
                        case('U'):  return( 79 );
                    }
                    break;

        case('B'):  switch( ch )
                    {   case('A'):  return( 56 );
                        case('E'):  return(  4 );
                        case('I'):  return( 83 );
                        case('K'):  return( 97 );
                        case('R'):  return( 35 );
                        case(' '):  return(  5 );
                    }
                    break;

        case('C'):  switch( ch )
                    {   case('A'):  return( 20 );
                        case('D'):  return( 48 );
                        case('E'):  return( 58 );
                        case('F'):  return( 98 );
                        case('L'):  return( 17 );
                        case('M'):  return( 96 );
                        case('O'):  return( 27 );
                        case('R'):  return( 24 );
                        case('S'):  return( 55 );
                        case('U'):  return( 29 );
                        case(' '):  return(  6 );
                    }
                    break;

        case('D'):  if( ch=='Y' )
                    {   return( 66 );
                    } else if( ch==' ' )
                        return( 1 );
                    break;

        case('E'):  if( ch=='R' )
                    {   return( 68 );
                    } else if( ch=='S' )
                    {   return( 99 );
                    } else if( ch=='U' )
                        return( 63 );
                    break;

        case('F'):  if( ch=='E' )
                    {   return(  26 );
                    } else if( ch=='M' )
                    {   return( 100 );
                    } else if( ch=='R' )
                    {   return(  87 );
                    } else if( ch=='F' )
                        return(   9 );
                    break;

        case('G'):  if( ch=='A' )
                    {   return( 31 );
                    } else if( ch=='D' )
                    {   return( 64 );
                    } else if( ch=='E' )
                        return( 32 );
                    break;

        case('H'):  if( ch=='E' )
                    {   return(  2 );
                    } else if( ch=='F' )
                    {   return( 72 );
                    } else if( ch=='G' )
                    {   return( 80 );
                    } else if( ch=='O' )
                    {   return( 67 );
                    } else if( ch==' ' )
                        return(  1 );
                    break;

        case('I'):  if( ch=='N' )
                    {   return( 49 );
                    } else if( ch=='R' )
                    {   return( 77 );
                    } else if( ch==' ' )
                        return( 53 );
                    break;

        case('K'):  if( ch=='R' )
                    {   return( 36 );
                    } else if( ch==' ' )
                        return( 19 );
                    break;

        case('L'):  if( ch=='A' )
                    {   return(  57 );
                    } else if( ch=='I' )
                    {   return(   3 );
                    } else if( (ch=='R') || (ch=='W') )
                    {   return( 103 );
                    } else if( ch=='U' )
                    {   return(  71 );
                    } else if( ch==' ' )
                        return(   1 );
                    break;

        case('M'):  if( ch=='D' )
                    {   return( 101 );
                    } else if( ch=='G' )
                    {   return(  12 );
                    } else if( ch=='N' )
                    {   return(  25 );
                    } else if( ch=='O' )
                        return(  42 );
                    break;

        case('N'):  switch( ch )
                    {   case('A'):  return(  11 );
                        case('B'):  return(  41 );
                        case('D'):  return(  60 );
                        case('E'):  return(  10 );
                        case('I'):  return(  28 );
                        case('O'):  return( 102 );
                        case('P'):  return(  93 );
                        case(' '):  return(   7 );
                    }
                    break;

        case('O'):  if( ch=='S' )
                    {   return( 76 );
                    } else if( ch==' ' )
                        return( 8 );
                    break;

        case('P'):  switch( ch )
                    {   case('A'):  return( 91 );
                        case('B'):  return( 82 );
                        case('D'):  return( 46 );
                        case('M'):  return( 61 );
                        case('O'):  return( 84 );
                        case('R'):  return( 59 );
                        case('T'):  return( 78 );
                        case('U'):  return( 94 );
                        case(' '):  return( 15 );
                    }
                    break;

        case('R'):  switch( ch )
                    {   case('A'):  return( 88 );
                        case('B'):  return( 37 );
                        case('E'):  return( 75 );
                        case('H'):  return( 45 );
                        case('N'):  return( 86 );
                        case('U'):  return( 44 );
                    }
                    break;

        case('S'):  switch( ch )
                    {   case('B'):  return( 51 );
                        case('C'):  return( 21 );
                        case('E'):  return( 34 );
                        case('I'):  return( 14 );
                        case('M'):  return( 62 );
                        case('N'):  return( 50 );
                        case('R'):  return( 38 );
                        case(' '):  return( 16 );
                    }
                    break;

        case('T'):  switch( ch )
                    {   case('A'):  return( 73 );
                        case('B'):  return( 65 );
                        case('C'):  return( 43 );
                        case('E'):  return( 52 );
                        case('H'):  return( 90 );
                        case('I'):  return( 22 );
                        case('L'):  return( 81 );
                        case('M'):  return( 69 );
                    }
                    break;

        case('U'):  if( ch==' ' )
                        return( 92 );
                    break;

        case('V'):  if( ch==' ' )
                        return( 23 );
                    break;

        case('W'):  if( ch==' ' )
                        return( 74 );
                    break;

        case('X'):  if( ch=='E' )
                        return( 54 );
                    break;

        case('Y'):  if( ch=='B' )
                    {   return( 70 );
                    } else if( ch==' ' )
                        return( 39 );
                    break;

        case('Z'):  if( ch=='N' )
                    {   return( 30 );
                    } else if( ch=='R' )
                        return( 40 );
                    break;
    }

    if( (*ptr>='0') && (*ptr<='9') )
        if( (ch=='H') || (ch=='D') )
            return( 1 ); /* Hydrogen */

    return( 0 );
}

char *OBChainsParser::ParseSmiles(char *ptr, int prev)
{
    char *name;
    int atomid;
    int next;
    int type;
    int ch;

    type = 0;
    while( (ch = *ptr++) )
    {   switch( ch )
        {   case('-'): type = BF_SINGLE;  break;
            case('='): type = BF_DOUBLE;  break;
            case('#'): type = BF_TRIPLE;  break;
            case('^'): type = BF_SINGLE|BF_AROMATIC; break;
            case('~'): type = BF_DOUBLE|BF_AROMATIC; break;

            case(')'): return( ptr );
            case('.'): prev = -1;  break;
            case('('): ptr = ParseSmiles(ptr,prev);
                       break;

            default:   atomid = ch-'0';
                       while( isdigit(*ptr) )
                           atomid = (atomid*10)+(*ptr++)-'0';

                       for( next=0; next<MonoAtomCount; next++ )
                           if( MonoAtom[next].atomid == atomid )
                               break;

                       if( next == MonoAtomCount )
                       {   name = ChainsAtomName[atomid];
                           MonoAtom[next].elem = IdentifyElement(name);
                           MonoAtom[next].atomid = atomid;
                           MonoAtom[next].bcount = 0;
                           MonoAtomCount++;
                       }

                       if( prev != -1 )
                       {   MonoBond[MonoBondCount].flag = type;
                           MonoBond[MonoBondCount].src = prev;
                           MonoBond[MonoBondCount].dst = next;
                           MonoBondCount++;

                           MonoAtom[prev].bcount++;
                           MonoAtom[next].bcount++;
                       }
                       prev = next;
        }
    }
    return( ptr-1 );
}

#ifdef _I_WANT_TO_OUTPUT_PDB_

static ChainsAtom *PDBOrder[MaxChainsAtom];

int PDBSort(ChainsAtom **arg1, ChainsAtom **arg2)
{
    ChainsAtom *atom1;
    ChainsAtom *atom2;

    atom1 = *arg1;
    atom2 = *arg2;

    if( atom1->chain != atom2->chain )
        return( atom1->chain - atom2->chain );

    if( atom1->hetflag != atom2->hetflag )
       return( atom1->hetflag? 1 : -1 );

    if( atom1->resno != atom2->resno )
        return( atom1->resno - atom2->resno );

    if( (atom1->elem==1) && (atom2->elem!=1) )
        return( 1 );
    if( (atom1->elem!=1) && (atom2->elem==1) )
        return( -1 );

    if( atom1->atomid != atom2->atomid )
        return( atom1->atomid - atom2->atomid );

    if( (atom1->elem==1) && (atom2->elem==1) )
        return( atom1->hcount - atom2->hcount );
    return( 0 );
}

static void OutputPDBFile(ChainsMolecule *mol, FILE *fp)
{
    int src,dst;
    ChainsAtom *atom;
    char *ptr;
    int i;

    for( i=0; i<mol->acount; i++ )
        PDBOrder[i] = &mol->atom[i];

#ifdef __STDC__
    qsort(PDBOrder,mol->acount,sizeof(ChainsAtom*),
          (int(*)(const void*,const void*))PDBSort);
#else
    qsort(PDBOrder,mol->acount,sizeof(ChainsAtom*),PDBSort);
#endif

    ptr = mol->name;
    while( *ptr == ' ' )
        ptr++;

    if( *ptr )
    {   fputs("COMPND    ",fp);
        while( *ptr )
            fputc(*ptr++,fp);
        fputc('\n',fp);
    }

    for( i=0; i<mol->acount; i++ )
    {   atom = PDBOrder[i];
        atom->serno = i+1;

        if( atom->hetflag )
        {   fputs("HETATM ",fp);
        } else fputs("ATOM   ",fp);

        fprintf(fp,"%4d ",atom->serno);

        if( atom->atomid == -1 )
        {   fprintf(fp,"%s  ", etab.GetSymbol(atom->elem));
        } else if( atom->elem == 1 )
        {   if( atom->hcount )
            {   fputc(atom->hcount+'0',fp);
            } else fputc(' ',fp);
            fprintf(fp,"H%.2s",ChainsAtomName[atom->atomid]+2);
        } else fprintf(fp,"%.4s",ChainsAtomName[atom->atomid]);

        fprintf(fp," %s ",ChainsResName[atom->resid]);
        fprintf(fp,"%c%4d",atom->chain,atom->resno);
        fprintf(fp,"    %8.3lf%8.3lf%8.3lf",atom->x,atom->y,atom->z);
        fputs("  1.00  0.00\n",fp);
    }

    for( i=0; i<mol->bcount; i++ )
        if( mol->bond[i].flag & BF_DOUBLE )
        {   src = mol->atom[mol->bond[i].src].serno;
            dst = mol->atom[mol->bond[i].dst].serno;
            fprintf(fp,"CONECT%5d%5d%5d\n",src,dst,dst);
            fprintf(fp,"CONECT%5d%5d%5d\n",dst,src,src);
        }
    fputs("END \n",fp);
}

#endif

} // End OpenEye Namespace

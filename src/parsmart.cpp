/**********************************************************************
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (c) 2001-2003 by Geoffrey R. Hutchison

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include <ctype.h>

#include "mol.h"
#include "bitvec.h"
#include "parsmart.h"

#ifndef True
#define True   1
#define False  0
#endif

/* Strict syntax checking! */
//#define STRICT
#define VERBOSE

#ifdef __sgi
#define UnusedArgument(x)   ((x)=(x))
#else
#define UnusedArgument(x)
#endif

using namespace std;

namespace OpenBabel {

  /*! \class OBSmartsPattern
      \brief

Substructure search is an incredibly useful tool in the context of a
small molecule programming library. Having an efficient substructure
search engine reduces the amount of hard code needed for molecule
perception, as well as increases the flexibility of certain
operations. For instance, atom typing can be easily performed based on
hard coded rules of element type and bond orders (or
hybridization). Alternatively, atom typing can also be done by
matching a set of substructure rules read at run time. In the latter
case customization based on application (such as changing the pH)
becomes a facile operation. Fortunately for Open Babel and its users,
Roger Sayle donated a SMARTS parser which became the basis for SMARTS
matching in Open Babel.

The SMARTS matcher, or OBSmartsPattern, is a separate object which can
match patterns in the OBMol class. The following code demonstrates how
to use the OBSmartsPattern class:
\code
OBMol mol(SDF,SDF);
cin >> mol;
OBSmartsPattern sp;
sp.Init("CC");
sp.Match(mol);
vector<vector<int> > maplist;
maplist = sp.GetMapList();
//or maplist = sp.GetUMapList();
//print out the results
vector<vector<int> >::iterator i;
vector<int>::iterator j;
for (i = maplist.begin();i != maplist.end();i++)
{
  for (j = i->begin();j != i->end();j++)
    cout << j << ' `;
  cout << endl;
}
\endcode
The preceding code reads in a molecule, initializes a smarts pattern
 of two single-bonded carbons, and locates all instances of the
 pattern in the molecule. Note that calling the Match() function
 does not return the results of the substructure match. The results
 from a match are stored in the OBSmartsPattern, and a call to
 GetMapList() or GetUMapList() must be made to extract the
 results. The function GetMapList() returns all matches of a
 particular pattern while GetUMapList() returns only the unique
 matches. For instance, the pattern [OD1]~C~[OD1] describes a
 carboxylate group. This pattern will match both atom number
 permutations of the carboxylate, and if GetMapList() is called, both
 matches will be returned. If GetUMapList() is called only unique
 matches of the pattern will be returned. A unique match is defined as
 one which does not cover the identical atoms that a previous match
 has covered.
  */

/*============================*/
/*  Period Table of Elements  */
/*============================*/

std::vector<std::pair<Pattern*,std::vector<bool> > > RSCACHE; //recursive smarts cache

typedef struct 
{
  char *symbol;
  int organic;
  int aromflag;
  double weight;
} Element;

#define ELEMMAX  104

typedef struct 
{
  BondExpr *closord[100];
  int       closure[100];
  int       closindex;
} ParseState;

/*
#define ATOMEXPRPOOL  16
#define BONDEXPRPOOL  16
#define ATOMPOOL      16
#define BONDPOOL      16
*/

#define ATOMEXPRPOOL  1
#define BONDEXPRPOOL  1
#define ATOMPOOL      1
#define BONDPOOL      1

/*=====================*/
/*  BondExpr Bit Sets  */
/*=====================*/

#define BF_NONRINGUNSPEC   0x0001
#define BF_NONRINGDOWN     0x0002
#define BF_NONRINGUP       0x0004
#define BF_NONRINGDOUBLE   0x0008
#define BF_NONRINGTRIPLE   0x0010
#define BF_RINGUNSPEC      0x0020
#define BF_RINGDOWN        0x0040
#define BF_RINGUP          0x0080
#define BF_RINGAROM        0x0100
#define BF_RINGDOUBLE      0x0200
#define BF_RINGTRIPLE      0x0400

#define BS_ALL             0x07FF
#define BS_SINGLE          0x00E7
#define BS_DOUBLE          0x0208
#define BS_TRIPLE          0x0410
#define BS_AROM            0x0100
#define BS_UP              0x0084
#define BS_DOWN            0x0042
#define BS_UPUNSPEC        0x00A5
#define BS_DOWNUNSPEC      0x0063
#define BS_RING            0x07E0
#define BS_DEFAULT         0x01E7

static char *MainPtr;
static char *LexPtr;

#define BUFMAX  1024
static char Buffer[BUFMAX];
static char Descr[BUFMAX];

static bool match(OBMol &mol,Pattern *pat,std::vector<std::vector<int> > &mlist,bool single=false);
static bool EvalAtomExpr(AtomExpr *expr,OBAtom *atom);
static bool EvalBondExpr(BondExpr *expr,OBBond *bond);
static int GetVectorBinding();
static int CreateAtom(Pattern*,AtomExpr*,int,int vb=0);

/*=============================*/
/*  Standard Utility Routines  */
/*=============================*/

static void FatalAllocationError( char *ptr )
{
    printf("Error: Unable to allocate %s!\n",ptr);
    exit(1);
}

/*================================*/
/*  Atom Expression Manipulation  */
/*================================*/

static void FreePattern( Pattern* );
static Pattern *CopyPattern( Pattern* );

static AtomExpr *AllocAtomExpr( void )
{
    register AtomExpr *result;

    result = (AtomExpr*)malloc(sizeof(AtomExpr));
    return result;
}

static AtomExpr *CopyAtomExpr( AtomExpr *expr )
{
    register AtomExpr *result;

    result = AllocAtomExpr();
    result->type = expr->type;
    switch( expr->type )
    {   case(AE_ANDHI):
        case(AE_ANDLO):
        case(AE_OR):    result->bin.lft = CopyAtomExpr(expr->bin.lft);
                        result->bin.rgt = CopyAtomExpr(expr->bin.rgt);
                        break;

        case(AE_NOT):   result->mon.arg = CopyAtomExpr(expr->mon.arg);
                        break;

        case(AE_RECUR): result->recur.recur = CopyPattern(
                                              (Pattern*)expr->recur.recur );
                        break;

        case(AE_LEAF):  result->leaf.prop = expr->leaf.prop;
                        result->leaf.value = expr->leaf.value;
                        break;
    }
    return result;
}

static void FreeAtomExpr( AtomExpr *expr )
{
    if( expr )
    {
		switch( expr->type )
        {   case(AE_ANDHI):
            case(AE_ANDLO):
            case(AE_OR):     FreeAtomExpr(expr->bin.lft);
                             FreeAtomExpr(expr->bin.rgt);
                             break;

            case(AE_NOT):    FreeAtomExpr(expr->mon.arg);
                             break;

            case(AE_RECUR):  FreePattern( (Pattern*)expr->recur.recur );
                             break;
        }
		if (expr)
		{
			free(expr);
			expr = (AtomExpr*)NULL;
		}
    }
}

static AtomExpr *BuildAtomLeaf( int prop, int val )
{
    register AtomExpr *result;

    result = AllocAtomExpr();
    result->leaf.type = AE_LEAF;
    result->leaf.prop = prop;
    result->leaf.value = val;
    return result;
}

static AtomExpr *BuildAtomNot( AtomExpr *expr )
{
    register AtomExpr *result;

    result = AllocAtomExpr();
    result->mon.type = AE_NOT;
    result->mon.arg = expr;
    return result;
}

static AtomExpr *BuildAtomBin( int op, AtomExpr *lft, AtomExpr *rgt )
{
    register AtomExpr *result;

    result = AllocAtomExpr();
    result->bin.type = op;
    result->bin.lft = lft;
    result->bin.rgt = rgt;
    return result;
}

static AtomExpr *BuildAtomRecurs( Pattern *pat )
{
  register AtomExpr *result;

  result = AllocAtomExpr();
  result->recur.type = AE_RECUR;
  result->recur.recur = (void*)pat;
  return result;
}

static AtomExpr *GenerateElement( int elem )
{
  return BuildAtomLeaf(AL_ELEM,elem);
}

static AtomExpr *GenerateAromElem( int elem, int flag )
{
  AtomExpr *expr1;
  AtomExpr *expr2;

  expr1 = BuildAtomLeaf(AL_AROM,flag);
  expr2 = BuildAtomLeaf(AL_ELEM,elem);
  return BuildAtomBin(AE_ANDHI,expr1,expr2);
}

static int IsInvalidAtom( AtomExpr *expr )
{
  if( !expr ) return True;
  return( (expr->type==AE_LEAF) &&
	  (expr->leaf.prop==AL_CONST)
	  && !expr->leaf.value );
}

/*================================*/
/*  Bond Expression Manipulation  */
/*================================*/

static BondExpr *AllocBondExpr( void )
{
    register BondExpr *result;

    result = (BondExpr*)malloc(sizeof(BondExpr));
    return result;
}

static BondExpr *CopyBondExpr( BondExpr *expr )
{
    register BondExpr *result;

    result = AllocBondExpr();
    result->type = expr->type;
    switch( expr->type )
    {   case(AE_ANDHI):
        case(AE_ANDLO):
        case(AE_OR):    result->bin.lft = CopyBondExpr(expr->bin.lft);
                        result->bin.rgt = CopyBondExpr(expr->bin.rgt);
                        break;

        case(AE_NOT):   result->mon.arg = CopyBondExpr(expr->mon.arg);
                        break;

        case(AE_LEAF):  result->leaf.prop = expr->leaf.prop;
                        result->leaf.value = expr->leaf.value;
                        break;
    }
    return result;
}

static void FreeBondExpr( BondExpr *expr )
{
    if( expr )
    {   switch( expr->type )
        {   case(BE_ANDHI):
            case(BE_ANDLO):
            case(BE_OR):     FreeBondExpr(expr->bin.lft);
                             FreeBondExpr(expr->bin.rgt);
                             break;

            case(BE_NOT):    FreeBondExpr(expr->mon.arg);
                             break;
        }

		if (expr)
		{
			free(expr);
			expr = (BondExpr*)NULL;
		}
    }
}

static BondExpr *BuildBondLeaf( int prop, int val )
{
   register BondExpr *result;

   result = AllocBondExpr();
   result->leaf.type = BE_LEAF;
   result->leaf.prop = prop;
   result->leaf.value = val;
   return result;
}

static BondExpr *BuildBondNot( BondExpr *expr )
{
    register BondExpr *result;

    result = AllocBondExpr();
    result->mon.type = BE_NOT;
    result->mon.arg = expr;
    return result;
}

static BondExpr *BuildBondBin( int op, BondExpr *lft, BondExpr *rgt )
{
    register BondExpr *result;

    result = AllocBondExpr();
    result->bin.type = op;
    result->bin.lft = lft;
    result->bin.rgt = rgt;
    return result;
}

static BondExpr *GenerateDefaultBond( void )
{
  register BondExpr *expr1;
  register BondExpr *expr2;

  expr1 = BuildBondLeaf(BL_TYPE,BT_SINGLE);
  expr2 = BuildBondLeaf(BL_TYPE,BT_AROM);
  return(BuildBondBin(BE_OR,expr1,expr2));
}

/*===============================*/
/*  SMARTS Pattern Manipulation  */
/*===============================*/

static Pattern *AllocPattern( void )
{
  Pattern *ptr;

  ptr = (Pattern*)malloc(sizeof(Pattern));
  if( !ptr ) FatalAllocationError("pattern");

  ptr->atom = (AtomSpec*)0;
  ptr->aalloc = 0;
  ptr->acount = 0;

  ptr->bond = (BondSpec*)0;
  ptr->balloc = 0;
  ptr->bcount = 0;

  ptr->parts = 1;
  return ptr;
}

static int CreateAtom( Pattern *pat, AtomExpr *expr, int part,int vb)
{
  int index,size;

  if( pat->acount == pat->aalloc )
    {   pat->aalloc += ATOMPOOL;
    size = (int)(pat->aalloc*sizeof(AtomSpec));
    if( pat->atom )
      {   pat->atom = (AtomSpec*)realloc(pat->atom,size);
      } else pat->atom = (AtomSpec*)malloc(size);
    if( !pat->atom ) FatalAllocationError("atom pool");
    }

  index = pat->acount++;
  pat->atom[index].part = part;
  pat->atom[index].expr = expr;
  pat->atom[index].vb = vb; //std::vector binding

  return index;
}

static int CreateBond( Pattern *pat, BondExpr *expr, int src, int dst )
{
  int index,size;

  if( pat->bcount == pat->balloc )
    {   pat->balloc += BONDPOOL;
    size = (int)(pat->balloc*sizeof(BondSpec));
    if( pat->bond )
      {   pat->bond = (BondSpec*)realloc(pat->bond,size);
      } else pat->bond = (BondSpec*)malloc(size);
    if( !pat->bond ) FatalAllocationError("bond pool");
    }

  index = pat->bcount++;
  pat->bond[index].expr = expr;
  pat->bond[index].src = src;
  pat->bond[index].dst = dst;
  return(index);
}

static Pattern *CopyPattern( Pattern *pat )
{
  Pattern *result;
  AtomExpr *aexpr;
  BondExpr *bexpr;
  int i;

  result = AllocPattern();
  result->parts = pat->parts;
  for( i=0; i<pat->acount; i++ )
    {   aexpr = CopyAtomExpr(pat->atom[i].expr);
    CreateAtom(result,aexpr,pat->atom[i].part);
    }

  for( i=0; i<pat->bcount; i++ )
    {   bexpr = CopyBondExpr(pat->bond[i].expr);
    CreateBond(result,bexpr,pat->bond[i].src,pat->bond[i].dst);
    }

  return result;
}

static void FreePattern( Pattern *pat )
{
  int i;

  if( pat )
    {
	  if( pat->aalloc )
      {   
		for( i=0; i<pat->acount; i++ )
			FreeAtomExpr(pat->atom[i].expr);
		free(pat->atom);
      }

	  if( pat->balloc )
      {
		  for( i=0; i<pat->bcount; i++ )
			  FreeBondExpr(pat->bond[i].expr);
		  free(pat->bond);
      }
	  free(pat);
    }
}

/*=========================*/
/*  SMARTS Syntax Parsing  */
/*=========================*/

static Pattern *ParseSMARTSPattern( void );
static Pattern *ParseSMARTSPart( Pattern*, int );

static Pattern *SMARTSError( Pattern *pat )
{
  char *ptr;

  fprintf(stderr,"SMARTS Error: %s\n",MainPtr);

  fputs("              ",stdout);
  for( ptr=MainPtr; ptr<LexPtr; ptr++ )
    fputc(' ',stdout);
  fputs("^\n",stdout);

  FreePattern(pat);
  return (Pattern*)0;
}

static AtomExpr *ParseSimpleAtomPrimitive( void )
{
    switch( *LexPtr++ )
    {   case '*':  return BuildAtomLeaf(AL_CONST,True);
#ifndef STRICT
        case 'A':  return BuildAtomLeaf(AL_AROM,False);
#endif
        case 'B':  if( *LexPtr == 'r' )
                   {   LexPtr++;
                       return GenerateElement(35);
                   }
                   return GenerateElement(5);

        case 'C':  if( *LexPtr == 'l' )
                   {   LexPtr++;
                       return GenerateElement(17);
                   }
                   return GenerateAromElem(6,False);

        case 'F':  return GenerateElement( 9);
        case 'I':  return GenerateElement(53);
        case 'N':  return GenerateAromElem( 7,False);
        case 'O':  return GenerateAromElem( 8,False);
        case 'P':  return GenerateElement(15);
        case 'S':  return GenerateAromElem(16,False);

#ifndef STRICT
        case 'a':  return BuildAtomLeaf(AL_AROM,True);
#endif
        case 'c':  return GenerateAromElem( 6,True);
        case 'n':  return GenerateAromElem( 7,True);
        case 'o':  return GenerateAromElem( 8,True);
        case 'p':  return GenerateAromElem(15,True);
        case 's':  return GenerateAromElem(16,True);
    }
    LexPtr--;
    return (AtomExpr*)0;
}

static AtomExpr *ParseComplexAtomPrimitive( void )
{
    register Pattern *pat;
    register int index;

    switch( *LexPtr++ )
    {   case('#'):  if( !isdigit(*LexPtr) )
                        return( (AtomExpr*)0 );

                    index = 0;
                    while( isdigit(*LexPtr) )
                        index = index*10 + ((*LexPtr++)-'0');
                    if( index > ELEMMAX )
                    {   LexPtr--;
                        return( (AtomExpr*)0 );
                    } else if( !index )
                        return( (AtomExpr*)0 );
                    return( GenerateElement(index) );

        case('$'):  if( *LexPtr != '(' )
                        return( (AtomExpr*)0 );
                    LexPtr++;
#ifdef STRICT
                    pat = ParseSMARTSPart(AllocPattern(),0);
#else
                    pat = ParseSMARTSPattern();
#endif
                    if( !pat ) return( (AtomExpr*)0 );
                    if( *LexPtr != ')' )
                    {   FreePattern(pat);
                        return( (AtomExpr*)0 );
                    }
                    LexPtr++;
                    return( BuildAtomRecurs(pat) );

        case('*'):  return( BuildAtomLeaf(AL_CONST,True) );

        case('+'):  if( isdigit(*LexPtr) )
                    {   index = 0;
                        while( isdigit(*LexPtr) )
                            index = index*10 + ((*LexPtr++)-'0');
                    } else
                    {   index = 1;
                        while( *LexPtr == '+' )
                        {   LexPtr++;
                            index++;
                        }
                    }
                    return( BuildAtomLeaf(AL_POSITIVE,index) );

        case('-'):  if( isdigit(*LexPtr) )
                    {   index = 0;
                        while( isdigit(*LexPtr) )
                            index = index*10 + ((*LexPtr++)-'0');
                    } else
                    {   index = 1;
                        while( *LexPtr == '-' )
                        {   LexPtr++;
                            index++;
                        }
                    }
                    return BuildAtomLeaf(AL_NEGATIVE,index);

        case '@':
	  if (*LexPtr != '@')
	    return(BuildAtomLeaf(AL_CHIRAL,AL_ANTICLOCKWISE));
	  LexPtr++;
	  return(BuildAtomLeaf(AL_CHIRAL,AL_CLOCKWISE));

        case '^': 
	  if (isdigit(*LexPtr))
	    {
	      index = 0;
	      while( isdigit(*LexPtr) )
		index = index*10 + ((*LexPtr++)-'0');
	      return(BuildAtomLeaf(AL_HYB,index));
	    }
	  else return(BuildAtomLeaf(AL_HYB,1));

        case('0'): case('1'): case('2'): case('3'): case('4'):
        case('5'): case('6'): case('7'): case('8'): case('9'):
                    index = LexPtr[-1]-'0';
                    while( isdigit(*LexPtr) )
                            index = index*10 + ((*LexPtr++)-'0');
                    return BuildAtomLeaf(AL_MASS,index);

        case('A'):  switch( *LexPtr++ )
                    {   case('c'):  return GenerateElement(89);
                        case('g'):  return GenerateElement(47);
                        case('l'):  return GenerateElement(13);
                        case('m'):  return GenerateElement(95);
                        case('r'):  return GenerateElement(18);
                        case('s'):  return GenerateElement(33);
                        case('t'):  return GenerateElement(85);
                        case('u'):  return GenerateElement(79);
                    }
                    LexPtr--;
                    return BuildAtomLeaf(AL_AROM,False);

        case('B'):  switch( *LexPtr++ )
                    {   case('a'):  return GenerateElement(56);
                        case('e'):  return GenerateElement( 4);
                        case('i'):  return GenerateElement(83);
                        case('k'):  return GenerateElement(97);
                        case('r'):  return GenerateElement(35);
                    }
                    LexPtr--;
                    return GenerateElement(5);

        case('C'):  switch( *LexPtr++ )
                    {   case('a'):  return GenerateElement(20);
                        case('d'):  return GenerateElement(48);
                        case('e'):  return GenerateElement(58);
                        case('f'):  return GenerateElement(98);
                        case('l'):  return GenerateElement(17);
                        case('m'):  return GenerateElement(96);
                        case('o'):  return GenerateElement(27);
                        case('r'):  return GenerateElement(24);
                        case('s'):  return GenerateElement(55);
                        case('u'):  return GenerateElement(29);
                    }
                    LexPtr--;
                    return GenerateAromElem(6,False);

        case('D'):  if( *LexPtr == 'y' )
                    {   LexPtr++;  return GenerateElement(66);
                    } else if( isdigit(*LexPtr) )
                    {   index = 0;
                        while( isdigit(*LexPtr) )
                            index = index*10 + ((*LexPtr++)-'0');
                        return BuildAtomLeaf(AL_DEGREE,index);
                    }
                    break;

        case('E'):  if( *LexPtr == 'r' )
                    {   LexPtr++;  return GenerateElement(68);
                    } else if( *LexPtr == 's' )
                    {   LexPtr++;  return GenerateElement(99);
                    } else if( *LexPtr == 'u' )
                    {   LexPtr++;  return GenerateElement(63);
                    }
                    break;

        case('F'):  if( *LexPtr == 'e' )
                    {   LexPtr++;  return GenerateElement(26);
                    } else if( *LexPtr == 'm' )
                    {   LexPtr++;  return GenerateElement(100);
                    } else if( *LexPtr == 'r' )
                    {   LexPtr++;  return GenerateElement(87);
                    }
                    return GenerateElement(9);

        case('G'):  if( *LexPtr == 'a' )
                    {   LexPtr++;  return( GenerateElement(31) );
                    } else if( *LexPtr == 'd' )
                    {   LexPtr++;  return( GenerateElement(64) );
                    } else if( *LexPtr == 'e' )
                    {   LexPtr++;  return( GenerateElement(32) );
                    }
                    break;

    case('H'):      if( *LexPtr == 'e' )
                    {   LexPtr++;  return( GenerateElement( 2) );
                    } else if( *LexPtr == 'f' )
                    {   LexPtr++;  return( GenerateElement(72) );
                    } else if( *LexPtr == 'g' )
                    {   LexPtr++;  return( GenerateElement(80) );
                    } else if( *LexPtr == 'o' )
                    {   LexPtr++;  return( GenerateElement(67) );
                    } else if( isdigit(*LexPtr) )
                    {   index = 0;
                        while( isdigit(*LexPtr) )
                            index = index*10 + ((*LexPtr++)-'0');
                        return( BuildAtomLeaf(AL_HCOUNT,index) );
                    }
                    return( BuildAtomLeaf(AL_HCOUNT,1) );
                    /* BuildAtomLeaf(AL_HCOUNT,1) ??? */
                    /* or else GenerateElement(1) ??? */

        case('I'):  if( *LexPtr == 'n' )
                    {   LexPtr++;  return( GenerateElement(49) );
                    } else if( *LexPtr == 'r' )
                    {   LexPtr++;  return( GenerateElement(77) );
                    }
                    return( GenerateElement(53) );

        case('K'):  if( *LexPtr == 'r' )
                    {   LexPtr++;  return( GenerateElement(36) );
                    }
                    return( GenerateElement(19) );

        case('L'):  if( *LexPtr == 'a' )
                    {   LexPtr++;  return( GenerateElement( 57) );
                    } else if( *LexPtr == 'i' )
                    {   LexPtr++;  return( GenerateElement(  3) );
                    } else if( *LexPtr == 'r' )
                    {   LexPtr++;  return( GenerateElement(103) );
                    } else if( *LexPtr == 'u' )
                    {   LexPtr++;  return( GenerateElement( 71) );
                    }
                    break;

        case('M'):  if( *LexPtr == 'd' )
                    {   LexPtr++;  return( GenerateElement(101) );
                    } else if( *LexPtr == 'g' )
                    {   LexPtr++;  return( GenerateElement( 12) );
                    } else if( *LexPtr == 'n' )
                    {   LexPtr++;  return( GenerateElement( 25) );
                    } else if( *LexPtr == 'o' )
                    {   LexPtr++;  return( GenerateElement( 42) );
                    }
                    break;

        case('N'):  switch( *LexPtr++ )
                    {   case('a'):  return( GenerateElement( 11) );
                        case('b'):  return( GenerateElement( 41) );
                        case('d'):  return( GenerateElement( 60) );
                        case('e'):  return( GenerateElement( 10) );
                        case('i'):  return( GenerateElement( 28) );
                        case('o'):  return( GenerateElement(102) );
                        case('p'):  return( GenerateElement( 93) );
                    }
                    LexPtr--;
                    return( GenerateAromElem(7,False) );

        case('O'):  if( *LexPtr == 's' )
                    {   LexPtr++;  return( GenerateElement(76) );
                    }
                    return( GenerateAromElem(8,False) );

        case('P'):  switch( *LexPtr++ )
                    {   case('a'):  return( GenerateElement(91) );
                        case('b'):  return( GenerateElement(82) );
                        case('d'):  return( GenerateElement(46) );
                        case('m'):  return( GenerateElement(61) );
                        case('o'):  return( GenerateElement(84) );
                        case('r'):  return( GenerateElement(59) );
                        case('t'):  return( GenerateElement(78) );
                        case('u'):  return( GenerateElement(94) );
                    }
                    LexPtr--;
                    return( GenerateElement(15) );

        case('R'):  switch( *LexPtr++ )
                    {   case('a'):  return( GenerateElement(88) );
                        case('b'):  return( GenerateElement(37) );
                        case('e'):  return( GenerateElement(75) );
                        case('h'):  return( GenerateElement(45) );
                        case('n'):  return( GenerateElement(86) );
                        case('u'):  return( GenerateElement(44) );
                    }
                    LexPtr--;
                    if( isdigit(*LexPtr) )
                    {   index = 0;
                        while( isdigit(*LexPtr) )
                            index = index*10 + ((*LexPtr++)-'0');
                    } else index = -1;
                    return( BuildAtomLeaf(AL_RINGS,index) );

        case('S'):  switch( *LexPtr++ )
                    {   case('b'):  return( GenerateElement(51) );
                        case('c'):  return( GenerateElement(21) );
                        case('e'):  return( GenerateElement(34) );
                        case('i'):  return( GenerateElement(14) );
                        case('m'):  return( GenerateElement(62) );
                        case('n'):  return( GenerateElement(50) );
                        case('r'):  return( GenerateElement(38) );
                    }
                    LexPtr--;
                    return( GenerateAromElem(16,False) );

        case('T'):  switch( *LexPtr++ )
                    {   case('a'):  return( GenerateElement(73) );
                        case('b'):  return( GenerateElement(65) );
                        case('c'):  return( GenerateElement(43) );
                        case('e'):  return( GenerateElement(52) );
                        case('h'):  return( GenerateElement(90) );
                        case('i'):  return( GenerateElement(22) );
                        case('l'):  return( GenerateElement(81) );
                        case('m'):  return( GenerateElement(69) );
                    }
                    LexPtr--;
                    break;

        case('U'):  return( GenerateElement(92) );
        case('V'):  return( GenerateElement(23) );
        case('W'):  return( GenerateElement(74) );

        case('X'):  if( *LexPtr == 'e' )
                    {   LexPtr++;  return( GenerateElement(54) );
                    } else if( isdigit(*LexPtr) )
                    {   index = 0;
                        while( isdigit(*LexPtr) )
                            index = index*10 + ((*LexPtr++)-'0');
                        return( BuildAtomLeaf(AL_CONNECT,index) );
                    }
                    break;

        case('Y'):  if( *LexPtr == 'b' )
                    {   LexPtr++;  return( GenerateElement(70) );
                    }
                    return( GenerateElement(39) );

        case('Z'):  if( *LexPtr == 'n' )
                    {   LexPtr++;  return GenerateElement(30);
                    } else if( *LexPtr == 'r' )
                    {   LexPtr++;  return GenerateElement(40);
                    }
                    break;

        case('a'):  if( *LexPtr == 's' )
                    {   LexPtr++;  return GenerateAromElem(33,True);
                    }
                    return BuildAtomLeaf(AL_AROM,True);

        case('c'):  return GenerateAromElem(6,True);

        case('h'):  if( isdigit(*LexPtr) )
                    {   index = 0;
                        while( isdigit(*LexPtr) )
                            index = index*10 + ((*LexPtr++)-'0');
                    } else index = 1;
                    return BuildAtomLeaf(AL_IMPLICIT,index);

        case('n'):  return GenerateAromElem(7,True);
        case('o'):  return GenerateAromElem(8,True);
        case('p'):  return GenerateAromElem(15,True);

        case('r'):  if( isdigit(*LexPtr) )
                    {   index = 0;
                        while( isdigit(*LexPtr) )
                            index = index*10 + ((*LexPtr++)-'0');
                        if( index == 0 )
                            return BuildAtomLeaf(AL_RINGS,0);
                        return BuildAtomLeaf(AL_SIZE,index);
                    }
                    return BuildAtomLeaf(AL_RINGS,-1);

        case('s'):  if( *LexPtr == 'i' )
                    {   LexPtr++;  return GenerateAromElem(14,True);
                    }
                    return GenerateAromElem(16,True);

        case('v'):  if( isdigit(*LexPtr) )
                    {   index = 0;
                        while( isdigit(*LexPtr) )
                            index = index*10 + ((*LexPtr++)-'0');
                        return BuildAtomLeaf(AL_VALENCE,index);
                    }
                    break;
    }
    LexPtr--;
    return (AtomExpr*)0;
}

static AtomExpr *ParseAtomExpr( int level )
{
    register AtomExpr *expr1;
    register AtomExpr *expr2;
    register char *prev;

    switch( level )
    {   case(0): /* Low Precedence Conjunction */
                 if( !(expr1=ParseAtomExpr(1)) )
                     return (AtomExpr*)0;

                 while( *LexPtr == ';' )
                 {   LexPtr++;
                     if( !(expr2=ParseAtomExpr(1)) )
                     {   FreeAtomExpr(expr1);
                         return (AtomExpr*)0;
                     }
                     expr1 = BuildAtomBin(AE_ANDLO,expr1,expr2);
                 }
                 return expr1;

        case(1): /* Disjunction */
                 if( !(expr1=ParseAtomExpr(2)) )
                     return (AtomExpr*)0;

                 while( *LexPtr == ',' )
                 {   LexPtr++;
                     if( !(expr2=ParseAtomExpr(2)) )
                     {   FreeAtomExpr(expr1);
                         return( (AtomExpr*)0 );
                     }
                     expr1 = BuildAtomBin(AE_OR,expr1,expr2);
                 }
                 return( expr1 );

        case(2): /* High Precedence Conjunction */
                 if( !(expr1=ParseAtomExpr(3)) )
                     return( (AtomExpr*)0 );

                 while( (*LexPtr!=']') && (*LexPtr!=';') &&
                        (*LexPtr!=',') && *LexPtr )
                 {   if( *LexPtr=='&' ) LexPtr++;
                     prev = LexPtr;
                     if( !(expr2=ParseAtomExpr(3)) )
                     {   if( prev != LexPtr )
                         {   FreeAtomExpr(expr1);
                             return( (AtomExpr*)0 );
                         } else return( expr1 );
                     }
                     expr1 = BuildAtomBin(AE_ANDHI,expr1,expr2);
                 }
                 return( expr1 );

        case(3): /* Negation or Primitive */
                 if( *LexPtr == '!' )
                 {   LexPtr++;
                     if( !(expr1=ParseAtomExpr(3)) )
                         return( (AtomExpr*)0 );
                     return( BuildAtomNot(expr1) );
                 }
                 return( ParseComplexAtomPrimitive() );
    }
    return (AtomExpr*)0;
}

static BondExpr *ParseBondPrimitive( void )
{
    switch( *LexPtr++ )
    {   case('-'):  return BuildBondLeaf(BL_TYPE,BT_SINGLE);
        case('='):  return BuildBondLeaf(BL_TYPE,BT_DOUBLE);
        case('#'):  return BuildBondLeaf(BL_TYPE,BT_TRIPLE);
        case(':'):  return BuildBondLeaf(BL_TYPE,BT_AROM);
        case('@'):  return BuildBondLeaf(BL_TYPE,BT_RING);
        case('~'):  return BuildBondLeaf(BL_CONST,True);

        case('/'):  if( *LexPtr == '?' )
                    {   LexPtr++;
                        return BuildBondLeaf(BL_TYPE,BT_UPUNSPEC);
                    }
                    return BuildBondLeaf(BL_TYPE,BT_UP);

        case('\\'): if( *LexPtr == '?' )
                    {   LexPtr++;
                        return BuildBondLeaf(BL_TYPE,BT_DOWNUNSPEC);
                    }
                    return BuildBondLeaf(BL_TYPE,BT_DOWN);
    }
    LexPtr--;
    return (BondExpr*)0;
}

static BondExpr *ParseBondExpr( int level )
{
    register BondExpr *expr1;
    register BondExpr *expr2;
    register char *prev;

    switch( level )
    {   case(0): /* Low Precedence Conjunction */
                 if( !(expr1=ParseBondExpr(1)) )
                     return (BondExpr*)0;

                 while( *LexPtr == ';' )
                 {   LexPtr++;
                     if( !(expr2=ParseBondExpr(1)) )
                     {   FreeBondExpr(expr1);
                         return (BondExpr*)0;
                     }
                     expr1 = BuildBondBin(BE_ANDLO,expr1,expr2);
                 }
                 return expr1;

        case(1): /* Disjunction */
                 if( !(expr1=ParseBondExpr(2)) )
                     return (BondExpr*)0;

                 while( *LexPtr == ',' )
                 {   LexPtr++;
                     if( !(expr2=ParseBondExpr(2)) )
                     {   FreeBondExpr(expr1);
                         return (BondExpr*)0;
                     }
                     expr1 = BuildBondBin(BE_OR,expr1,expr2);
                 }
                 return expr1;

        case(2): /* High Precedence Conjunction */
                 if( !(expr1=ParseBondExpr(3)) )
                     return (BondExpr*)0;

                 while( (*LexPtr!=']') && (*LexPtr!=';') &&
                        (*LexPtr!=',') && *LexPtr )
                 {   if( *LexPtr == '&' ) LexPtr++;
                     prev = LexPtr;
                     if( !(expr2=ParseBondExpr(3)) )
                     {   if( prev != LexPtr )
                         {   FreeBondExpr(expr1);
                             return (BondExpr*)0;
                         } else return expr1;
                     }
                     expr1 = BuildBondBin(BE_ANDHI,expr1,expr2);
                 }
                 return expr1;

        case(3): /* Negation or Primitive */
                 if( *LexPtr == '!' )
                 {   LexPtr++;
                     if( !(expr1=ParseBondExpr(3)) )
                         return (BondExpr*)0;
                     return BuildBondNot(expr1);
                 }
                 return ParseBondPrimitive();
    }
    return (BondExpr*)0;
}

static int GetVectorBinding()
{
  int vb=0;

  LexPtr++; //skip colon
  if(isdigit(*LexPtr))
    {
      vb = 0;
      while( isdigit(*LexPtr) )
	vb = vb*10 + ((*LexPtr++)-'0');
    }

  return(vb);
}

static Pattern *ParseSMARTSError( Pattern *pat, BondExpr *expr )
{
    if( expr ) FreeBondExpr(expr);
    return SMARTSError(pat);
}

static Pattern *SMARTSParser( Pattern *pat, ParseState *stat,
                              int prev, int part )
{
  int vb = 0;
  register AtomExpr *aexpr;
  register BondExpr *bexpr;
  register int index;

    bexpr = (BondExpr*)0;

    while( *LexPtr )
    {   switch( *LexPtr++ )
        {   case('.'):  if( bexpr || (prev==-1) )
                            return ParseSMARTSError(pat,bexpr);
                        prev = -1;
                        break;

            case('-'):  case('='):  case('#'):
            case(':'):  case('~'):  case('@'):
            case('/'):  case('\\'): case('!'):
                        LexPtr--;
                        if( (prev==-1) || bexpr )
                            return ParseSMARTSError(pat,bexpr);
                        if( !(bexpr=ParseBondExpr(0)) )
                            return ParseSMARTSError(pat,bexpr);
                        break;

            case('('):
#ifdef STRICT
                        if( (prev==-1) || bexpr )
                        {   LexPtr--;
                            return ParseSMARTSError(pat,bexpr);
                        }
                        pat = SMARTSParser(pat,stat,prev,part);
                        if( !pat ) return (Pattern*)0;
#else /* STRICT */
                        if( bexpr )
                        {   LexPtr--;
                            return ParseSMARTSError(pat,bexpr);
                        }
                        if( prev == -1 )
                        {   index = pat->acount;
                            pat = SMARTSParser(pat,stat,-1,part);
                            if( !pat ) return( (Pattern*)0 );
                            if( index == pat->acount )
                                return ParseSMARTSError(pat,bexpr);
                            prev = index;
                        } else
                        {   pat = SMARTSParser(pat,stat,prev,part);
                            if( !pat ) return (Pattern*)0;
                        }
#endif /* STRICT */

                        if( *LexPtr != ')' )
                            return ParseSMARTSError(pat,bexpr);
                        LexPtr++;
                        break;

            case(')'):  LexPtr--;
                        if( (prev==-1) || bexpr )
                            return ParseSMARTSError(pat,bexpr);
                        return pat;

            case('%'):  if( prev == -1 )
                        {   LexPtr--;
                            return ParseSMARTSError(pat,bexpr);
                        }

                        if( isdigit(LexPtr[0]) && isdigit(LexPtr[1]) )
                        {   index = 10*(LexPtr[0]-'0') + (LexPtr[1]-'0');
                            LexPtr += 2;
                        } else return ParseSMARTSError(pat,bexpr);

                        if( stat->closure[index] == -1 )
                        {   stat->closord[index] = bexpr;
                            stat->closure[index] = prev;
                        } else if( stat->closure[index] != prev )
                        {   FreeBondExpr(stat->closord[index]);
                            if( !bexpr ) bexpr = GenerateDefaultBond();
                            CreateBond(pat,bexpr,prev,stat->closure[index]);
                            stat->closure[index] = -1;
                            bexpr = (BondExpr*)0;
                        } else return ParseSMARTSError(pat,bexpr);
                        break;

            case('0'):  case('1'):  case('2'):
            case('3'):  case('4'):  case('5'):
            case('6'):  case('7'):  case('8'):
            case('9'):  LexPtr--;
                        if( prev == -1 )
                            return ParseSMARTSError(pat,bexpr);
                        index = (*LexPtr++)-'0';

                        if( stat->closure[index] == -1 )
                        {   stat->closord[index] = bexpr;
                            stat->closure[index] = prev;
                            bexpr = (BondExpr*)0;
                        } else if( stat->closure[index] != prev )
                        {   FreeBondExpr(stat->closord[index]);
                            if( !bexpr ) bexpr = GenerateDefaultBond();
                            CreateBond(pat,bexpr,prev,stat->closure[index]);
                            stat->closure[index] = -1;
                            bexpr = (BondExpr*)0;
                        } else return ParseSMARTSError(pat,bexpr);
                        break;

            case('['):  aexpr = ParseAtomExpr(0);
	                vb = (*LexPtr == ':') ? GetVectorBinding():0;
                        if( !aexpr || (*LexPtr!=']') )
                            return ParseSMARTSError(pat,bexpr);
                        index = CreateAtom(pat,aexpr,part,vb);
                        if( prev != -1 )
                        {   if( !bexpr ) bexpr = GenerateDefaultBond();
                            CreateBond(pat,bexpr,prev,index);
                            bexpr = (BondExpr*)0;
                        }
                        prev = index;
                        LexPtr++;
                        break;

            default:    LexPtr--;
                        aexpr = ParseSimpleAtomPrimitive();
                        if( !aexpr ) return ParseSMARTSError(pat,bexpr);
                        index = CreateAtom(pat,aexpr,part);
                        if( prev != -1 )
                        {   if( !bexpr ) bexpr = GenerateDefaultBond();
                            CreateBond(pat,bexpr,prev,index);
                            bexpr = (BondExpr*)0;
                        }
                        prev = index;
        }
    }

    if( (prev==-1) || bexpr )
        return ParseSMARTSError(pat,bexpr);

    return pat;
}

static void MarkGrowBonds(Pattern *pat)
{
  int i;
  OBBitVec bv;

  for (i = 0;i < pat->bcount;i++)
    {
      pat->bond[i].grow = (bv[pat->bond[i].src] && bv[pat->bond[i].dst])?
	false:true;

      bv.SetBitOn(pat->bond[i].src);
      bv.SetBitOn(pat->bond[i].dst);
    }
}

static int GetChiralFlag(AtomExpr *expr)
{
  int size=0;
#define OB_EVAL_STACKSIZE 40
  AtomExpr *stack[OB_EVAL_STACKSIZE];
  memset(stack,'\0',sizeof(AtomExpr*)*OB_EVAL_STACKSIZE);
#undef OB_EVAL_STACKSIZE
  bool lftest=true;

  for (size=0,stack[size] = expr;size >= 0;expr=stack[size])
    {
    switch (expr->type)
      {
         case AE_LEAF:
	   if (expr->leaf.prop == AL_CHIRAL) return(expr->leaf.value);
	   size--;
	   break;

         case AE_ANDHI: 
         case AE_ANDLO: 

	   if (stack[size+1] == expr->bin.rgt)      size--;
	   else if (stack[size+1] == expr->bin.lft)
	     {
	       if (lftest){size++;stack[size] = expr->bin.rgt;}
	       else	  size--;
	     }
	   else     {size++;stack[size] = expr->bin.lft;}
	   break;

         case AE_OR:

	   if (stack[size+1] == expr->bin.rgt)      size--;
	   else if (stack[size+1] == expr->bin.lft)
	     {
	       if (!lftest) {size++;stack[size] = expr->bin.rgt;}
	       else size--;
	     }
	   else  {size++;stack[size] = expr->bin.lft;}
	   break;

         case AE_NOT:
	   if (stack[size+1] != expr->mon.arg) 
	     {size++;stack[size] = expr->mon.arg;}
	   else {lftest = !lftest; size--;}
	   break;

      case AE_RECUR:
	size--;
	break;
      }
    }

  return((int)false);
}

static Pattern *ParseSMARTSPart( Pattern *result, int part )
{
  auto ParseState stat;
  int i,flag;

  for( i=0; i<100; i++ )
    stat.closure[i] = -1;

  result = SMARTSParser(result,&stat,-1,part);

  flag = False;
  for( i=0; i<100; i++ )
    if( stat.closure[i] != -1 )
      {   FreeBondExpr(stat.closord[i]);
      flag = True;
      }

  if( result )
    {
      if( flag )
	return(SMARTSError(result));
      else
	{
	  MarkGrowBonds(result);
	  result->ischiral = false;
	  for (i = 0;i < result->acount;i++)
	    {
	      result->atom[i].chiral_flag = GetChiralFlag(result->atom[i].expr);
	      if (result->atom[i].chiral_flag)
		result->ischiral = true;
	    }
	  return(result);
	}
    } else return (Pattern*)0;
}


static Pattern *ParseSMARTSPattern( void )
{
  Pattern *result;
  result = AllocPattern();

  while( *LexPtr == '(' )
    {   LexPtr++;
    result = ParseSMARTSPart(result,result->parts);
    if( !result ) return (Pattern*)0;
    result->parts++;

    if( *LexPtr != ')' )
      return SMARTSError(result);
    LexPtr++;

    if( !*LexPtr || (*LexPtr==')') )
      return result;

    if( *LexPtr != '.' )
      return SMARTSError(result);
    LexPtr++;
    }

  return ParseSMARTSPart(result,0);
}

static Pattern *ParseSMARTSString( char *ptr )
{
    register Pattern *result;

    if( !ptr || !*ptr )
        return (Pattern*)0;

    LexPtr = MainPtr = ptr;
    result = ParseSMARTSPattern();
    if( result && *LexPtr )
        return SMARTSError(result);
    return result;
}

Pattern *ParseSMARTSRecord( char *ptr )
{
    register char *src,*dst;

    src = ptr;
    while( *src && !isspace(*src) )
        src++;

    if( isspace(*src) )
    {   *src++ = '\0';
        while( isspace(*src) )
            src++;
    }

    dst = Descr;
    while( *src && (dst<Descr+78) )
    {   if( isspace(*src) )
        {   *dst++ = ' ';
            while( isspace(*src) )
                src++;
        } else *dst++ = *src++;
    }
    *dst = '\0';

    return ParseSMARTSString(Buffer);
}

/*==============================*/
/*  SMARTS Component Traversal  */
/*==============================*/

static void TraverseSMARTS( Pattern *pat, int i )
{
    register int j,k;

    pat->atom[i].visit = True;
    for( j=0; j<pat->bcount; j++ )
        if( pat->bond[j].visit == -1 )
        {   if( pat->bond[j].src == i )
            {   pat->bond[j].visit = i;
                k = pat->bond[j].dst;
                if( !pat->atom[k].visit )
                    TraverseSMARTS(pat,k);
            } else if( pat->bond[j].dst == i )
            {   pat->bond[j].visit = i;
                k = pat->bond[j].src;
                if( !pat->atom[k].visit )
                    TraverseSMARTS(pat,k);
            }
        }
}

/*============================*/
/*  Canonical SMARTS Pattern  */
/*============================*/

static AtomExpr *NotAtomExpr( AtomExpr* );
static AtomExpr *AndAtomExpr( AtomExpr*, AtomExpr* );
static AtomExpr *OrAtomExpr( AtomExpr*, AtomExpr* );
//static AtomExpr *TransformAtomExpr( AtomExpr* );
//static Pattern *CanonicaliseSMARTS( Pattern* );

static int IsBooleanAtomLeaf( AtomExpr *expr )
{
    return (expr->leaf.prop==AL_AROM) ||
           (expr->leaf.prop==AL_CONST);
}

static int IsNegatingAtomLeaf( AtomExpr *expr )
{
    return (expr->leaf.prop==AL_RINGS);
}

static int EqualAtomExpr( AtomExpr *lft, AtomExpr *rgt )
{
    if( lft->type != rgt->type )
        return False;

    if( lft->type == AE_LEAF )
    {   return( (lft->leaf.prop==rgt->leaf.prop) &&
                (lft->leaf.value==rgt->leaf.value) );
    } else if( lft->type == AE_NOT )
    {   return EqualAtomExpr(lft->mon.arg,rgt->mon.arg);
    } else if( lft->type == AE_RECUR )
        return False;

    return EqualAtomExpr(lft->bin.lft,rgt->bin.lft) &&
           EqualAtomExpr(lft->bin.rgt,rgt->bin.rgt);
}

static int OrderAtomExpr( AtomExpr *lft, AtomExpr *rgt )
{
    register AtomExpr *larg;
    register AtomExpr *rarg;
    register int stat;

    if( lft->type == AE_NOT )
    {   /* larg->type == AE_LEAF */
        larg = lft->mon.arg;
    } else larg = lft;

    if( rgt->type == AE_NOT )
    {   /* rarg->type == AE_LEAF */
        rarg = rgt->mon.arg;
    } else rarg = rgt;

    if( larg->type > rarg->type )
    {   return  1;
    } else if( larg->type < rarg->type )
        return -1;

    if( larg->type == AE_LEAF )
    {   if( larg->leaf.prop > rarg->leaf.prop )
            return  1;
        if( larg->leaf.prop < rarg->leaf.prop )
            return -1;
        return( larg->leaf.value - rarg->leaf.value );
    }

    stat = OrderAtomExpr(lft->bin.lft,rgt->bin.lft);
    if( stat != 0 ) return stat;
    return OrderAtomExpr(lft->bin.rgt,rgt->bin.rgt);
}

static int AtomLeafConflict( AtomExpr *lft, AtomExpr *rgt )
{
    register AtomExpr *tmp;

    if( (lft->type==AE_LEAF) && (rgt->type==AE_LEAF) )
    {   if( lft->leaf.prop == rgt->leaf.prop )
        {   if( IsNegatingAtomLeaf(lft) )
            {   if( lft->leaf.value == 0 )
                {   return rgt->leaf.value != 0;
                } else if( lft->leaf.value == -1 )
                    return rgt->leaf.value == 0;

                if( rgt->leaf.value == 0 )
                {   return lft->leaf.value != 0;
                } else if( rgt->leaf.value == -1 )
                    return lft->leaf.value == 0;
            }
            return lft->leaf.value != rgt->leaf.value;
        }

        if( lft->leaf.prop > rgt->leaf.prop )
        {   tmp = lft;
            lft = rgt;
            rgt = tmp;
        }

        /* Aromaticity -> Ring */
        if( (lft->leaf.prop==AL_AROM) && (rgt->leaf.prop==AL_RINGS) )
            return( lft->leaf.value && !rgt->leaf.value );

        /* Positive charge ~ Negative charge */
        if( (lft->leaf.prop==AL_NEGATIVE) && (rgt->leaf.prop==AL_POSITIVE) )
            return( (lft->leaf.value!=0) || (rgt->leaf.value!=0) );

        /* Total hcount >= Implicit hcount */
        if( (lft->leaf.prop==AL_HCOUNT) && (rgt->leaf.prop==AL_IMPLICIT) )
            return( lft->leaf.value < rgt->leaf.value );
    }

    if( (lft->type==AE_LEAF) && (rgt->type==AE_NOT) )
    {   rgt = rgt->mon.arg;
        if( (lft->leaf.prop==AL_NEGATIVE) && (rgt->leaf.prop==AL_POSITIVE) )
            return( (lft->leaf.value==0) && (rgt->leaf.value==0) );
        if( (lft->leaf.prop==AL_POSITIVE) && (rgt->leaf.prop==AL_NEGATIVE) )
            return( (lft->leaf.value==0) && (rgt->leaf.value==0) );
        return False;
    }

    if( (lft->type==AE_NOT) && (rgt->type==AE_LEAF) )
    {   lft = lft->mon.arg;
        if( (lft->leaf.prop==AL_NEGATIVE) && (rgt->leaf.prop==AL_POSITIVE) )
            return( (lft->leaf.value==0) && (rgt->leaf.value==0) );
        if( (lft->leaf.prop==AL_POSITIVE) && (rgt->leaf.prop==AL_NEGATIVE) )
            return( (lft->leaf.value==0) && (rgt->leaf.value==0) );
        return False;
    }

    return False;
}

static int AtomExprConflict( AtomExpr *lft, AtomExpr *rgt )
{
    while( rgt->type == AE_ANDHI )
    {   if( AtomLeafConflict(lft,rgt->bin.lft) )
            return True;
        rgt = rgt->bin.rgt;
    }
    return AtomLeafConflict(lft,rgt);
}

/* return LEAF(lft) => LEAF(rgt); */
static int AtomLeafImplies( AtomExpr *lft, AtomExpr *rgt )
{
    if( (lft->type==AE_LEAF) && (rgt->type==AE_LEAF) )
    {   /* Implied Ring Membership */
        if( (rgt->leaf.prop==AL_RINGS) && (rgt->leaf.value==-1) )
        {   if( lft->leaf.prop == AL_AROM )
                return lft->leaf.value;

            if( lft->leaf.prop == AL_RINGS )
                return lft->leaf.value > 0;

            if( lft->leaf.prop == AL_SIZE )
                return lft->leaf.value > 0;
        }

        /* Positive charge ~ Negative charge */
        if( (lft->leaf.prop==AL_POSITIVE) && (rgt->leaf.prop==AL_NEGATIVE) )
            return (lft->leaf.value==0) && (rgt->leaf.value==0);
        return False;
    }

    if( (lft->type==AE_LEAF) && (rgt->type==AE_NOT) )
    {   rgt = rgt->mon.arg;
        if( lft->leaf.prop == rgt->leaf.prop )
            return lft->leaf.value != rgt->leaf.value;

        if( (lft->leaf.prop==AL_POSITIVE) && (rgt->leaf.prop==AL_NEGATIVE) )
            return True;
        if( (lft->leaf.prop==AL_NEGATIVE) && (rgt->leaf.prop==AL_POSITIVE) )
            return True;
        return False;
    }

    return False;
}

/* return EXPR(rgt) => LEAF(lft); */
static int AtomExprImplied( AtomExpr *lft, AtomExpr *rgt )
{
    while( rgt->type == AE_ANDHI )
    {   if( AtomLeafImplies(rgt->bin.lft,lft) )
            return True;
        rgt = rgt->bin.rgt;
    }
    return AtomLeafImplies(rgt,lft);
}

/* remove implied nodes from EXPR(rgt) */
static AtomExpr *AtomExprImplies( AtomExpr *lft, AtomExpr *rgt )
{
    register AtomExpr *tmp;

    if( rgt->type != AE_ANDHI )
    {   if( AtomLeafImplies(lft,rgt) )
        {   FreeAtomExpr(rgt);
            return (AtomExpr*)0;
        }
        return rgt;
    }

    tmp = AtomExprImplies(lft,rgt->bin.rgt);

    if( tmp )
    {   if( AtomLeafImplies(lft,rgt->bin.lft) )
        {   rgt->bin.rgt = (AtomExpr*)0;
            FreeAtomExpr(rgt);
            return tmp;
        }
        rgt->bin.rgt = tmp;
        return rgt;
    } else
    {   rgt->bin.rgt = (AtomExpr*)0;
        if( AtomLeafImplies(lft,rgt->bin.lft) )
        {   FreeAtomExpr(rgt);
            return (AtomExpr*)0;
        }
        tmp = rgt->bin.lft;
        rgt->bin.lft = (AtomExpr*)0;
        FreeAtomExpr(rgt);
        return tmp;
    }
}

static AtomExpr *AndAtomExprLeaf( AtomExpr *lft, AtomExpr *rgt )
{
    if( AtomExprConflict(lft,rgt) )
    {   FreeAtomExpr(lft);
        FreeAtomExpr(rgt);
        return BuildAtomLeaf(AL_CONST,False);
    }

    if( AtomExprImplied(lft,rgt) )
    {   FreeAtomExpr(lft);
        return rgt;
    }

    rgt = AtomExprImplies(lft,rgt);
    if( !rgt ) return lft;

    return BuildAtomBin(AE_ANDHI,lft,rgt);
}

static AtomExpr *ConstrainRecursion( AtomExpr *recur, AtomExpr *expr )
{
    register AtomExpr *head;
    register Pattern *pat;

    pat = (Pattern*)recur->recur.recur;
    head = AndAtomExpr(pat->atom[0].expr,expr);
    pat->atom[0].expr = head;

    if( IsInvalidAtom(head) )
    {   FreePattern(pat);
        return BuildAtomLeaf(AL_CONST,False);
    }
    return recur;
}

static AtomExpr *AndAtomExpr( AtomExpr *lft, AtomExpr *rgt )
{
    register AtomExpr *expr;
    register int order;

    /* Identities */
    if( EqualAtomExpr(lft,rgt) )
    {   FreeAtomExpr(rgt);
        return lft;
    }

    if( (lft->type==AE_LEAF) && (lft->leaf.prop==AL_CONST) )
    {   if( lft->leaf.value )
        {   FreeAtomExpr(lft);
            return rgt;
        } else
        {   FreeAtomExpr(rgt);
            return lft;
        }
    }

    if( (rgt->type==AE_LEAF) && (rgt->leaf.prop==AL_CONST) )
    {   if( rgt->leaf.value )
        {   FreeAtomExpr(rgt);
            return lft;
        } else
        {   FreeAtomExpr(lft);
            return rgt;
        }
    }

    /*  Distributivity  */
    if( lft->type == AE_OR )
    {   expr = CopyAtomExpr(rgt);
        expr = OrAtomExpr(AndAtomExpr(expr,lft->bin.lft),
                          AndAtomExpr(rgt, lft->bin.rgt));
        lft->bin.lft = (AtomExpr*)0;
        lft->bin.rgt = (AtomExpr*)0;
        FreeAtomExpr(lft);
        return( expr );
    }

    if( rgt->type == AE_OR )
    {   expr = CopyAtomExpr(lft);
        expr = OrAtomExpr(AndAtomExpr(expr,rgt->bin.lft),
                          AndAtomExpr(lft, rgt->bin.rgt));
        rgt->bin.lft = (AtomExpr*)0;
        rgt->bin.rgt = (AtomExpr*)0;
        FreeAtomExpr(rgt);
        return( expr );
    }

    /* Recursion */
    if( (rgt->type==AE_RECUR) && (lft->type!=AE_RECUR) )
        return ConstrainRecursion(rgt,lft);

    if( (rgt->type!=AE_RECUR) && (lft->type==AE_RECUR) )
        return ConstrainRecursion(lft,rgt);

    order = OrderAtomExpr(lft,rgt);
    if( order > 0 )
    {   expr = lft;
        lft = rgt;
        rgt = expr;
    }

    if( lft->type == AE_ANDHI )
    {   expr = AndAtomExpr(lft->bin.rgt,rgt);
        expr = AndAtomExpr(lft->bin.lft,expr);
        lft->bin.lft = (AtomExpr*)0;
        lft->bin.rgt = (AtomExpr*)0;
        FreeAtomExpr(lft);
        return expr;
    }

    if( rgt->type == AE_ANDHI )
    {   if( OrderAtomExpr(lft,rgt->bin.lft) > 0 )
        {   expr = AndAtomExpr(lft,rgt->bin.rgt);
            expr = AndAtomExpr(rgt->bin.lft,expr);
            rgt->bin.lft = (AtomExpr*)0;
            rgt->bin.rgt = (AtomExpr*)0;
            FreeAtomExpr(rgt);
            return expr;
        }

        if( EqualAtomExpr(lft,rgt->bin.lft) )
        {   FreeAtomExpr(lft);
            return rgt;
        }
    }

    return AndAtomExprLeaf(lft,rgt);
}

static AtomExpr *OrAtomExprLeaf( AtomExpr *lft, AtomExpr *rgt )
{
    return BuildAtomBin(AE_OR,lft,rgt);
}

static AtomExpr *OrAtomExpr( AtomExpr *lft, AtomExpr *rgt )
{
    register AtomExpr *expr;
    register int order;

    /* Identities */
    if( EqualAtomExpr(lft,rgt) )
    {   FreeAtomExpr(rgt);
        return lft;
    }

    if( (lft->type==AE_LEAF) && (lft->leaf.prop==AL_CONST) )
    {   if( lft->leaf.value )
        {   FreeAtomExpr(rgt);
            return lft;
        } else
        {   FreeAtomExpr(lft);
            return rgt;
        }
    }

    if( (rgt->type==AE_LEAF) && (rgt->leaf.prop==AL_CONST) )
    {   if( rgt->leaf.value )
        {   FreeAtomExpr(lft);
            return rgt;
        } else
        {   FreeAtomExpr(rgt);
            return lft;
        }
    }

    order = OrderAtomExpr(lft,rgt);
    if( order > 0 )
    {   expr = lft;
        lft = rgt;
        rgt = expr;
    }

    if( lft->type == AE_OR )
    {   expr = OrAtomExpr(lft->bin.rgt,rgt);
        expr = OrAtomExpr(lft->bin.lft,expr);
        lft->bin.lft = (AtomExpr*)0;
        lft->bin.rgt = (AtomExpr*)0;
        FreeAtomExpr(lft);
        return expr;
    }

    if( rgt->type == AE_OR )
    {   if( OrderAtomExpr(lft,rgt->bin.lft) > 0 )
        {   expr = OrAtomExpr(lft,rgt->bin.rgt);
            expr = OrAtomExpr(rgt->bin.lft,expr);
            rgt->bin.lft = (AtomExpr*)0;
            rgt->bin.rgt = (AtomExpr*)0;
            FreeAtomExpr(rgt);
            return expr;
        }

        if( EqualAtomExpr(lft,rgt->bin.lft) )
        {   FreeAtomExpr(lft);
            return rgt;
        }
    }

    return OrAtomExprLeaf(lft,rgt);
}

static AtomExpr *NotAtomExpr( AtomExpr *expr )
{
    register AtomExpr *result;
    register AtomExpr *lft;
    register AtomExpr *rgt;

    if( expr->type == AE_LEAF )
    {   if( IsBooleanAtomLeaf(expr) )
        {   expr->leaf.value = !expr->leaf.value;
            return expr;
        } else if( IsNegatingAtomLeaf(expr) )
        {   if( expr->leaf.value == -1 )
            {   expr->leaf.value = 0;
                return expr;
            } else if( expr->leaf.value == 0 )
            {   expr->leaf.value = -1;
                return expr;
            }
        }
    } else if( expr->type == AE_NOT )
    {   result = expr->mon.arg;
        expr->mon.arg = (AtomExpr*)0;
        FreeAtomExpr(expr);
        return result;
    } else if( (expr->type==AE_ANDHI) ||
               (expr->type==AE_ANDLO) )
    {
        lft = NotAtomExpr(expr->bin.lft);
        rgt = NotAtomExpr(expr->bin.rgt);
        expr->bin.lft = (AtomExpr*)0;
        expr->bin.rgt = (AtomExpr*)0;
        FreeAtomExpr(expr);
        return OrAtomExpr(lft,rgt);
    } else if( expr->type == AE_OR )
    {   lft = NotAtomExpr(expr->bin.lft);
        rgt = NotAtomExpr(expr->bin.rgt);
        expr->bin.lft = (AtomExpr*)0;
        expr->bin.rgt = (AtomExpr*)0;
        FreeAtomExpr(expr);
        return AndAtomExpr(lft,rgt);
    }
    return BuildAtomNot(expr);
}

/*==============================*/
/*  Canonical Bond Expressions  */
/*==============================*/

static int GetBondLeafIndex( BondExpr *expr )
{
    if( expr->leaf.prop == BL_CONST )
    {   if( expr->leaf.value )
        {   return( BS_ALL );
        } else return( 0 );
    } else /* expr->leaf.prop == BL_TYPE */
        switch( expr->leaf.value )
        {   case(BT_SINGLE):     return( BS_SINGLE );
            case(BT_DOUBLE):     return( BS_DOUBLE );
            case(BT_TRIPLE):     return( BS_TRIPLE );
            case(BT_AROM):       return( BS_AROM );
            case(BT_UP):         return( BS_UP );
            case(BT_DOWN):       return( BS_DOWN );
            case(BT_UPUNSPEC):   return( BS_UPUNSPEC );
            case(BT_DOWNUNSPEC): return( BS_DOWNUNSPEC );
            case(BT_RING):       return( BS_RING );
        }
    return 0;
}

static int GetBondExprIndex( BondExpr *expr )
{
    register int lft,rgt;
    register int arg;

    switch( expr->type )
    {   case(BE_LEAF):   return GetBondLeafIndex(expr);

        case(BE_NOT):    arg = GetBondExprIndex(expr->mon.arg);
                         return( arg ^ BS_ALL );

        case(BE_ANDHI):
        case(BE_ANDLO):  lft = GetBondExprIndex(expr->bin.lft);
                         rgt = GetBondExprIndex(expr->bin.rgt);
                         return( lft & rgt );

        case(BE_OR):     lft = GetBondExprIndex(expr->bin.lft);
                         rgt = GetBondExprIndex(expr->bin.rgt);
                         return( lft | rgt );
    }
    /* Avoid Compiler Warning */
    return 0;
}

static BondExpr *NotBondExpr( BondExpr *expr )
{
    register BondExpr *result;

    if( expr->type == BE_LEAF )
    {   if( expr->leaf.prop == BL_CONST )
        {   expr->leaf.value = !expr->leaf.value;
            return expr;
        }
    } else if( expr->type == BE_NOT )
    {   result = expr->mon.arg;
        expr->mon.arg = (BondExpr*)0;
        FreeBondExpr(expr);
        return result;
    }
    return BuildBondNot(expr);
}

static BondExpr *TransformBondExpr( BondExpr *expr )
{
    register BondExpr *lft,*rgt;
    register BondExpr *arg;

    if( expr->type == BE_LEAF )
    {   return expr;
    } else if( expr->type == BE_NOT )
    {   arg = expr->mon.arg;
        arg = TransformBondExpr(arg);
        expr->mon.arg = (BondExpr*)0;
        FreeBondExpr(expr);
        return NotBondExpr(arg);
    } else if( expr->type == BE_ANDHI )
    {   lft = expr->bin.lft;
        rgt = expr->bin.rgt;
        lft = TransformBondExpr(lft);
        rgt = TransformBondExpr(rgt);
        expr->bin.lft = lft;
        expr->bin.rgt = rgt;
        return expr;
    } else if( expr->type == BE_ANDLO )
    {   lft = expr->bin.lft;
        rgt = expr->bin.rgt;
        lft = TransformBondExpr(lft);
        rgt = TransformBondExpr(rgt);
        expr->bin.lft = lft;
        expr->bin.rgt = rgt;
        return expr;
    } else if( expr->type == BE_OR )
    {   lft = expr->bin.lft;
        rgt = expr->bin.rgt;
        lft = TransformBondExpr(lft);
        rgt = TransformBondExpr(rgt);
        expr->bin.lft = lft;
        expr->bin.rgt = rgt;
        return expr;
    }
    return expr;
}

#ifdef FOO
static BondExpr *CanonicaliseBond( BondExpr *expr )
{
#ifndef ORIG
    register int index;

    index = GetBondExprIndex(expr);
    FreeBondExpr(expr);

    LexPtr = CanBondExpr[index];
    if( *LexPtr )
    {   expr = ParseBondExpr(0);
    } else expr = GenerateDefaultBond();
#endif
    return TransformBondExpr(expr);
}
#endif


//**********************************
//********Pattern Matching**********
//**********************************

bool OBSmartsPattern::Init(const char *buffer)
{
    strcpy(Buffer,buffer);

    _pat = ParseSMARTSRecord(Buffer);
    _str = buffer;

    return(_pat != (Pattern*)NULL);
}

bool OBSmartsPattern::Init(const std::string &s)
{
    strcpy(Buffer, s.c_str());

    _pat = ParseSMARTSRecord(Buffer);
    _str = s;

    return(_pat != (Pattern*)NULL);
}

OBSmartsPattern::~OBSmartsPattern()
{
  if (_pat) FreePattern(_pat);
}

bool OBSmartsPattern::Match(OBMol &mol,bool single)
{
  RSCACHE.clear();
  return(match(mol,_pat,_mlist,single));
}

bool OBSmartsPattern::RestrictedMatch(OBMol &mol,std::vector<std::pair<int,int> > &pr,bool single)
{
  bool ok;
  std::vector<std::vector<int> > mlist;
  std::vector<std::vector<int> >::iterator i;
  std::vector<std::pair<int,int> >::iterator j;

  RSCACHE.clear();
  match(mol,_pat,mlist);
  _mlist.clear();
  if (mlist.empty()) return(false);

  for (i = mlist.begin();i != mlist.end();i++)
    {
      ok = true;
      for (j = pr.begin();j != pr.end() && ok;j++)
		  if ((*i)[j->first] != j->second)
			  ok = false;
  
      if (ok) _mlist.push_back(*i);
      if (single && !_mlist.empty()) return(true);
    }

  return((_mlist.empty()) ? false:true);
}

bool OBSmartsPattern::RestrictedMatch(OBMol &mol,OBBitVec &vres,bool single)
{
  bool ok;
  std::vector<int>::iterator j;
  std::vector<std::vector<int> > mlist;
  std::vector<std::vector<int> >::iterator i;

  RSCACHE.clear();
  match(mol,_pat,mlist);

  _mlist.clear();
  if (mlist.empty()) return(false);

  for (i = mlist.begin();i != mlist.end();i++)
    {
      ok = true;
      for (j = i->begin();j != i->end();j++)
	if (!vres[*j])
	  {
	    ok = false;
	    break;
	  }
      if (!ok) continue;
  
      _mlist.push_back(*i);
      if (single && !_mlist.empty()) return(true);
    }

  return((_mlist.empty()) ? false:true);
}

void SetupAtomMatchTable(std::vector<std::vector<bool> > &ttab,Pattern *pat,OBMol &mol)
{
  int i;

  ttab.resize(pat->acount);
  for (i = 0;i < pat->acount;i++)
    ttab[i].resize(mol.NumAtoms()+1);

  OBAtom *atom;
  std::vector<OBNodeBase*>::iterator j;
  for (i = 0;i < pat->acount;i++)
    for (atom = mol.BeginAtom(j);atom;atom = mol.NextAtom(j))
      if (EvalAtomExpr(pat->atom[0].expr,atom))
	ttab[i][atom->GetIdx()] = true;
}

static void FastSingleMatch(OBMol &mol,Pattern *pat,std::vector<std::vector<int> > &mlist)
{
  OBAtom *atom,*a1,*nbr;
  std::vector<OBNodeBase*>::iterator i;

  OBBitVec bv(mol.NumAtoms()+1);
  std::vector<int> map; map.resize(pat->acount);
  std::vector<std::vector<OBEdgeBase*>::iterator> vi;
  std::vector<bool> vif;

  if (pat->bcount) {
     vif.resize(pat->bcount);
     vi.resize(pat->bcount);
  }

  int bcount;
  for (atom = mol.BeginAtom(i);atom;atom=mol.NextAtom(i))
    if (EvalAtomExpr(pat->atom[0].expr,atom))
      {
	map[0] = atom->GetIdx();
	if (pat->bcount) vif[0] = false;
	bv.Clear();
        bv.SetBitOn(atom->GetIdx());

	for (bcount=0;bcount >=0;)
	  {
	    //***entire pattern matched***
	    if (bcount == pat->bcount) //save full match here
	      {
		mlist.push_back(map);
		bcount--;
		return; //found a single match
	      }
	    
	    //***match the next bond***
	    if (!pat->bond[bcount].grow) //just check bond here
	      {
		if ( !vif[bcount] )
		  {
		    OBBond *bond = mol.GetBond(map[pat->bond[bcount].src],
                                               map[pat->bond[bcount].dst]);
		    if (bond && EvalBondExpr(pat->bond[bcount].expr,bond))
		      {
			vif[bcount++] = true;
			if (bcount < pat->bcount)
                            vif[bcount] = false;
		      }
		    else
		      bcount--;
		  }
		else //bond must have already been visited - backtrack
		  bcount--;
	      }
	    else //need to map atom and check bond
	      {
		a1 = mol.GetAtom(map[pat->bond[bcount].src]);

		if (!vif[bcount]) //figure out which nbr atom we are mapping
		  {
                    nbr = a1->BeginNbrAtom(vi[bcount]);
                  }
		else
		  {
		    bv.SetBitOff(map[pat->bond[bcount].dst]);
		    nbr = a1->NextNbrAtom(vi[bcount]);
		  }

		for (;nbr;nbr=a1->NextNbrAtom(vi[bcount]))
		  if (!bv[nbr->GetIdx()])
		  if (EvalAtomExpr(pat->atom[pat->bond[bcount].dst].expr,nbr)
		      && EvalBondExpr(pat->bond[bcount].expr,(OBBond *)*(vi[bcount])))
		    {
		      bv.SetBitOn(nbr->GetIdx());
		      map[pat->bond[bcount].dst] = nbr->GetIdx();
                      vif[bcount] = true;
		      bcount++;
		      if (bcount < pat->bcount) 
                         vif[bcount] = false;
		      break;
		    }
		
		if (!nbr)//no match - time to backtrack
		    bcount--;
	      }
	  }
      }
}


static bool match(OBMol &mol,Pattern *pat,std::vector<std::vector<int> > &mlist,bool single)
{
  mlist.clear();
  if (!pat || pat->acount == 0) return(false);//shouldn't ever happen

  if (single && !pat->ischiral)
    FastSingleMatch(mol,pat,mlist);
  else
    {
      OBSSMatch ssm(mol,pat);
      ssm.Match(mlist);
    }

  if (pat->ischiral && mol.Has3D())
    {
      int j,k,r1,r2,r3,r4;
      std::vector<std::vector<int> >::iterator m;
      OBAtom *ra1,*ra2,*ra3,*ra4;
      std::vector<std::vector<int> > tmpmlist;

      for (j = 0;j < pat->acount;j++)
	if (pat->atom[j].chiral_flag)
	  {
	    r1 = r2 = r3 = r4 = -1;
	    r2 = j;
	    for (k = 0;k < pat->bcount;k++)
	      if (pat->bond[k].dst == r2)
		if (r1 == -1) r1 = pat->bond[k].src;
		else if (r3 == -1) r3 = pat->bond[k].src;
		else if (r4 == -1) r4 = pat->bond[k].src;

	    for (k = 0;k < pat->bcount;k++)
	      if (pat->bond[k].src == r2)
		if (r1 == -1) r1 = pat->bond[k].dst;
		else if (r3 == -1) r3 = pat->bond[k].dst;
		else if (r4 == -1) r4 = pat->bond[k].dst;
		
	    if (r1 == -1 || r2 == -1 || r3 == -1 || r4 == -1) continue;

	    tmpmlist.clear();
	    for (m = mlist.begin();m != mlist.end();m++)
	      {
		ra1 = mol.GetAtom((*m)[r1]);
		ra2 = mol.GetAtom((*m)[r2]);
		ra3 = mol.GetAtom((*m)[r3]);
		ra4 = mol.GetAtom((*m)[r4]);
		double sign = CalcTorsionAngle(ra1->GetVector(),
					      ra2->GetVector(),
					      ra3->GetVector(),
					      ra4->GetVector());
		if (sign > 0.0 && pat->atom[j].chiral_flag == AL_ANTICLOCKWISE)
		  continue;
		if (sign < 0.0 && pat->atom[j].chiral_flag == AL_CLOCKWISE)
		  continue;

		//ok - go ahead and save it
		tmpmlist.push_back(*m);
	      }
	    mlist = tmpmlist;
	  }
    }
  
  return(!mlist.empty());
}

#define RECURSIVE

#ifdef RECURSIVE
static bool EvalAtomExpr(AtomExpr *expr,OBAtom *atom)
{
  for (;;)
    switch (expr->type)
      {
      case AE_LEAF:
	switch( expr->leaf.prop )
	  {
	  case AL_ELEM:
            return(expr->leaf.value == (int)atom->GetAtomicNum());
	  case AL_AROM:
            if( !expr->leaf.value )
                return !atom->IsAromatic();
            return atom->IsAromatic(); 
	  case AL_HCOUNT:
	    if (atom->ExplicitHydrogenCount() > atom->ImplicitHydrogenCount())
	      return (expr->leaf.value==(signed int)atom->ExplicitHydrogenCount());
	    else
	      return (expr->leaf.value==(signed int)atom->ImplicitHydrogenCount());
	    /* Roger's broken code --mts
            register int hcount;
            hcount = atom->ExplicitHydrogenCount() +
                     atom->ImplicitHydrogenCount();
		     return expr->leaf.value == hcount;
	    */

	  case AL_DEGREE:  
	    return(expr->leaf.value == (int)atom->GetHvyValence());
	  case AL_VALENCE: 
	    return(expr->leaf.value == (int)atom->KBOSum());
	  case AL_CONNECT: 
	    return(expr->leaf.value == (int)atom->GetImplicitValence());
	  case AL_NEGATIVE: 
	    return(expr->leaf.value == -(atom->GetFormalCharge()));
	  case AL_POSITIVE: 
	    return(expr->leaf.value == atom->GetFormalCharge());
	  case AL_HYB:
            return(expr->leaf.value == (int)atom->GetHyb());

	  case AL_RINGS:     
	    if( expr->leaf.value == -1 ) return atom->IsInRing();
	    else if( expr->leaf.value == 0 ) return !atom->IsInRing();
	    else return expr->leaf.value == (int)atom->MemberOfRingCount();

	  case AL_SIZE:
            if( expr->leaf.value == -1 ) return atom->IsInRing();
	    if (!expr->leaf.value) return !atom->IsInRing();
	    else return atom->IsInRingSize(expr->leaf.value);

	  case AL_IMPLICIT:  
	    return expr->leaf.value == (int)atom->ImplicitHydrogenCount();

	  case AL_CONST:
            if( !expr->leaf.value )
                return false;
            return(!atom->IsHydrogen());  /* ??? FIXME */
	  default:
            return false;
	  }
 
      case AE_NOT:   return(!EvalAtomExpr(expr->mon.arg,atom));
      case AE_ANDHI: /* Same as AE_ANDLO */
      case AE_ANDLO: 
	if( !EvalAtomExpr(expr->bin.lft,atom)) return(false);
	expr = expr->bin.rgt;
	break;
      case AE_OR:
	if(EvalAtomExpr(expr->bin.lft,atom)) return(true);
	expr = expr->bin.rgt;
	break;

      case AE_RECUR:
      {
	//see if pattern has been matched
	std::vector<std::pair<Pattern*,std::vector<bool> > >::iterator i;
	for (i = RSCACHE.begin();i != RSCACHE.end();i++)
	  if (i->first == (Pattern*)expr->recur.recur)
	    return(i->second[atom->GetIdx()]);

	//perceive and match pattern
	std::vector<std::vector<int> >::iterator j;
	std::vector<bool> vb(((OBMol*) atom->GetParent())->NumAtoms()+1);
	std::vector<std::vector<int> > mlist;
	if (match( *((OBMol *) atom->GetParent()),
		   (Pattern*)expr->recur.recur,mlist))
	  for (j = mlist.begin();j != mlist.end();j++)
	    vb[(*j)[0]] = true;
	
	RSCACHE.push_back(std::pair<Pattern*,std::vector<bool> > ((Pattern*)expr->recur.recur,vb));

	return(vb[atom->GetIdx()]);
      }

      default:
          return(false);
      }
}

#else

static bool EvalAtomExpr(AtomExpr *expr,OBAtom *atom)
{
  int size=0;
#define OB_EVAL_STACKSIZE 40
  AtomExpr *stack[OB_EVAL_STACKSIZE];
  memset(stack,'\0',sizeof(AtomExpr*)*OB_EVAL_STACKSIZE);
#undef OB_EVAL_STACKSIZE
  bool lftest=true;

  for (size=0,stack[size] = expr;size >= 0;expr=stack[size])
    {
    switch (expr->type)
      {
         case AE_LEAF:
	   switch( expr->leaf.prop )
	     {
	       //expr->leaf.value
	     case AL_ELEM: 
	       lftest = (expr->leaf.value == atom->GetAtomicNum());
	       break;
	     case AL_AROM:
	       lftest = (expr->leaf.value == (int)atom->IsAromatic()); 
	       break;
	     case AL_HCOUNT:
	       if (atom->ExplicitHydrogenCount() > atom->ImplicitHydrogenCount())
		 lftest=(expr->leaf.value==atom->ExplicitHydrogenCount());
	       else
		 lftest=(expr->leaf.value==atom->ImplicitHydrogenCount());
	       break;
	     case AL_DEGREE:
	       lftest = (expr->leaf.value == atom->GetHvyValence());
	       break;
	     case AL_VALENCE:
	       lftest = (expr->leaf.value == atom->BOSum());
	       break;
	     case AL_CONNECT:   //X
	       lftest = (expr->leaf.value == atom->GetImplicitValence());
	       break;
	     case AL_NEGATIVE:  
	       lftest=(expr->leaf.value == -1*(atom->GetFormalCharge()));
	       break;
	     case AL_POSITIVE: 
	       lftest=(expr->leaf.value == atom->GetFormalCharge());
	       break;
	     case AL_HYB:
	       lftest=(expr->leaf.value == atom->GetHyb());
	       break;
	     case AL_RINGS:     
	       if (expr->leaf.value == -1)  lftest = (atom->IsInRing());
	       else
		 if (expr->leaf.value == 0) lftest = !(atom->IsInRing());
		 else
		   lftest=(atom->MemberOfRingCount()==expr->leaf.value);
	       break;
	     case AL_SIZE:
	       if (!expr->leaf.value) lftest = !atom->IsInRing();
	       else lftest = atom->IsInRingSize(expr->leaf.value);
	       break; 

	     case AL_IMPLICIT:  
	       lftest=(expr->leaf.value==atom->ImplicitHydrogenCount());
	       break;
	     case AL_CONST:     lftest=!atom->IsHydrogen();  break;
	     case AL_MASS:      break;
	     default:           break;
	     }
	   size--;
	   break;

         case AE_ANDHI: 

	   if (stack[size+1] == expr->bin.rgt)      size--;
	   else if (stack[size+1] == expr->bin.lft)
	     {
	       if (lftest){size++;stack[size] = expr->bin.rgt;}
	       else	  size--;
	     }
	   else     {size++;stack[size] = expr->bin.lft;}
	   break;

         case AE_OR:

	   if (stack[size+1] == expr->bin.rgt)      size--;
	   else if (stack[size+1] == expr->bin.lft)
	     {
	       if (!lftest) {size++;stack[size] = expr->bin.rgt;}
	       else size--;
	     }
	   else  {size++;stack[size] = expr->bin.lft;}
	   break;

         case AE_ANDLO: 

	   if (stack[size+1] == expr->bin.rgt)      size--;
	   else if (stack[size+1] == expr->bin.lft)
	     {
	       if (lftest){size++;stack[size] = expr->bin.rgt;}
	       else	  size--;
	     }
	   else     {size++;stack[size] = expr->bin.lft;}
	   break;

         case AE_NOT:
	   if (stack[size+1] != expr->mon.arg) 
	     {size++;stack[size] = expr->mon.arg;}
	   else {lftest = !lftest; size--;}
	   break;

      case AE_RECUR:
	//see if pattern has been matched
	bool matched=false;

	std::vector<std::pair<Pattern*,std::vector<bool> > >::iterator i;
	for (i = RSCACHE.begin();i != RSCACHE.end();i++)
	  if (i->first == (Pattern*)expr->recur.recur)
	    {
	      lftest = i->second[atom->GetIdx()];
	      matched = true;
	      break;
	    }

	if (!matched)
	  {
	    std::vector<bool> vb(atom->GetParent()->NumAtoms()+1);
	    std::vector<std::vector<int> > mlist;
	    lftest = false;
	    if (match((*atom->GetParent()),(Pattern*)expr->recur.recur,mlist))
	      {
		std::vector<std::vector<int> >::iterator i;
		for (i = mlist.begin();i != mlist.end();i++)
		  {
		    if ((*i)[0] == atom->GetIdx()) lftest = true;
		    vb[(*i)[0]] = true;
		  }
	      }
	    RSCACHE.push_back(std::pair<Pattern*,std::vector<bool> > ((Pattern*)expr->recur.recur,vb));
	  }

	size--;
	break;
      }
    }

  return(lftest);
}
#endif

#ifdef RECURSIVE

static bool EvalBondExpr(BondExpr *expr,OBBond *bond)
{
  for (;;)
    switch( expr->type )
      {
      case BE_LEAF: 

		  if( expr->leaf.prop == BL_CONST ) return((expr->leaf.value != 0) ? true : false);
	else 
	  switch( expr->leaf.value )
	    {   
	    case BT_SINGLE: return(bond->GetBO() == 1 && !bond->IsAromatic());
	    case BT_AROM: return(bond->IsAromatic());
	    case BT_DOUBLE: return(bond->GetBO()==2 && !bond->IsAromatic());    
	    case BT_TRIPLE: return(bond->GetBO()==3);
	    case BT_RING: return(bond->IsInRing());
	    default: return(false);
	      /*
		need to handle these cases
		case BT_UP:
		case BT_DOWN:       
		case BT_UPUNSPEC:
		case BT_DOWNUNSPEC:
	      */
	    }


      case BE_NOT: return(!EvalBondExpr(expr->mon.arg,bond));
      case BE_ANDHI:
      case BE_ANDLO:
	if (!EvalBondExpr(expr->bin.lft,bond)) return(false);
	expr = expr->bin.rgt;
	break;

      case BE_OR:
	if (EvalBondExpr(expr->bin.lft,bond)) return(true);
	expr = expr->bin.rgt;
	break;
      default:
        return false;
      }
}

#else

static bool EvalBondExpr(BondExpr *expr,OBBond *bond)
{
  int size=0;
#define OB_EVAL_STACKSIZE 40
  BondExpr *stack[OB_EVAL_STACKSIZE];
  memset(stack,'\0',sizeof(AtomExpr*)*OB_EVAL_STACKSIZE);
#undef OB_EVAL_STACKSIZE
  bool lftest=true;

  for (size=0,stack[size] = expr;size >= 0;expr=stack[size])
    switch( expr->type )
    {   case(BE_LEAF): 

          if( expr->leaf.prop == BL_CONST )
	    lftest = (expr->leaf.value)?true:false;
	  else /* expr->leaf.prop == BL_TYPE */
	    switch( expr->leaf.value )
	      {   
	         case(BT_SINGLE):    
		   lftest = (bond->GetBO() == 1 && !bond->IsAromatic());  break;
	         case(BT_DOUBLE):     
		   lftest = (bond->GetBO()==2 && !bond->IsAromatic());  break;
	         case(BT_TRIPLE):
		   lftest = (bond->GetBO()==3);  break;
	         case(BT_AROM): lftest=bond->IsAromatic(); break;
	         case(BT_RING): lftest=bond->IsInRing();   break;
	         case(BT_UP):         break;
	         case(BT_DOWN):       break;
	         case(BT_UPUNSPEC):   break;
	         case(BT_DOWNUNSPEC): break;
	      }
	      size--;
	      break;

        case(BE_NOT):    
	   if (stack[size+1] != expr->mon.arg) 
	     {size++;stack[size] = expr->mon.arg;}
	   else {lftest = !lftest; size--;}
	   break;

        case(BE_ANDHI):
	  if (stack[size+1] == expr->bin.rgt)      size--;
	  else if (stack[size+1] == expr->bin.lft)
	    {
	      if (lftest){size++;stack[size] = expr->bin.rgt;}
	      else	  size--;
	    }
	  else     {size++;stack[size] = expr->bin.lft;}
	  break;
	   
        case(BE_ANDLO):
	   if (stack[size+1] == expr->bin.rgt)      size--;
	   else if (stack[size+1] == expr->bin.lft)
	     {
	       if (lftest){size++;stack[size] = expr->bin.rgt;}
	       else	  size--;
	     }
	   else     {size++;stack[size] = expr->bin.lft;}
	   break;

        case(BE_OR):
	  if (stack[size+1] == expr->bin.rgt)      size--;
	   else if (stack[size+1] == expr->bin.lft)
	     {
	       if (!lftest) {size++;stack[size] = expr->bin.rgt;}
	       else size--;
	     }
	   else  {size++;stack[size] = expr->bin.lft;}
	   break;
    }
  return(lftest);
}
#endif

std::vector<std::vector<int> > &OBSmartsPattern::GetUMapList()
{
  if (_mlist.empty() || _mlist.size() == 1) return(_mlist);

  bool ok;
  OBBitVec bv;
  std::vector<OBBitVec> vbv;
  std::vector<std::vector<int> > mlist;
  std::vector<std::vector<int> >::iterator i;
  std::vector<OBBitVec>::iterator j;

  for (i = _mlist.begin();i != _mlist.end();i++)
    {
      ok = true;
      bv.Clear(); bv.FromVecInt(*i);
      for (j = vbv.begin();j != vbv.end() && ok;j++)
	if ((*j) == bv)
	  ok = false;

      if (ok)
	{
	  mlist.push_back(*i);
	  vbv.push_back(bv);
	}
    }

  _mlist = mlist;
  return(_mlist);
}

void OBSmartsPattern::WriteMapList(ostream &ofs)
{
    std::vector<std::vector<int> >::iterator i;
    std::vector<int>::iterator j;

    for ( i = _mlist.begin() ; i != _mlist.end() ; i++ )
    {
        for (j = (*i).begin();j != (*i).end();j++)
	  ofs << *j << ' ' << ends;
        ofs << endl;
    }
}

//*******************************************************************
//  The OBSSMatch class performs exhaustive matching using recursion
//  Explicit stack handling is used to find just a single match in 
//  match()
//*******************************************************************

OBSSMatch::OBSSMatch(OBMol &mol,Pattern *pat)
{
  _mol = &mol;
  _pat = pat;
  _map.resize(pat->acount);

  if (!mol.Empty())
    {
      _uatoms = new bool [mol.NumAtoms()+1];
      memset((char*)_uatoms,'\0',sizeof(bool)*(mol.NumAtoms()+1));
    }
  else
    _uatoms = (bool*)NULL;
}

OBSSMatch::~OBSSMatch()
{
  if (_uatoms) delete [] _uatoms;
}

void OBSSMatch::Match(std::vector<std::vector<int> > &mlist,int bidx)
{
  if (bidx == -1)
    {
      OBAtom *atom;
      std::vector<OBNodeBase*>::iterator i;
      for (atom = _mol->BeginAtom(i);atom;atom = _mol->NextAtom(i))
	if (EvalAtomExpr(_pat->atom[0].expr,atom))
	  {
	    _map[0] = atom->GetIdx();
	    _uatoms[atom->GetIdx()] = true;
	    Match(mlist,0);
	    _map[0] = 0;
	    _uatoms[atom->GetIdx()] = false;
	  }
      return;
    }

  if (bidx == _pat->bcount) //save full match here
    {
      mlist.push_back(_map);
      return;
    }

  if (_pat->bond[bidx].grow) //match the next bond
    {
      int src,dst;
      src = _pat->bond[bidx].src;
      dst = _pat->bond[bidx].dst;

      AtomExpr *aexpr = _pat->atom[dst].expr;
      BondExpr *bexpr = _pat->bond[bidx].expr;
      OBAtom *atom,*nbr;
      std::vector<OBEdgeBase*>::iterator i;
      atom = _mol->GetAtom(_map[src]);
      for (nbr = atom->BeginNbrAtom(i);nbr;nbr = atom->NextNbrAtom(i))
	if (!_uatoms[nbr->GetIdx()] && EvalAtomExpr(aexpr,nbr) && 
	    EvalBondExpr(bexpr,((OBBond*) *i)))
	  {
	    _map[dst] = nbr->GetIdx();
	    _uatoms[nbr->GetIdx()] = true;
	    Match(mlist,bidx+1);
	    _uatoms[nbr->GetIdx()] = false;
	    _map[dst] = 0;
	  }
    } 
  else //just check bond here
    {
      OBBond *bond = _mol->GetBond(_map[_pat->bond[bidx].src],
				   _map[_pat->bond[bidx].dst]);
      if (bond && EvalBondExpr(_pat->bond[bidx].expr,bond))
	Match(mlist,bidx+1);
    }
}

static int GetExprOrder(BondExpr *expr)
{
  int size=0;
  BondExpr *stack[15];
  memset(stack,'\0',sizeof(AtomExpr*)*15);
  bool lftest=true;

  for (size=0,stack[size] = expr;size >= 0;expr=stack[size])
    switch( expr->type )
    {   case(BE_LEAF): 

          if( expr->leaf.prop == BL_CONST ) lftest = true;
	  else /* expr->leaf.prop == BL_TYPE */
	    switch( expr->leaf.value )
	      {   
	         case(BT_SINGLE):    return(1);
	         case(BT_DOUBLE):    return(2);
	         case(BT_TRIPLE):    return(3);
	         case(BT_AROM):      return(5);
	      default: lftest = true; 
	      }
	      size--;
	      break;

        case(BE_NOT):    return(0);
        case(BE_ANDHI):
        case(BE_ANDLO):
        case(BE_OR):
	  if (stack[size+1] == expr->bin.rgt)      size--;
	  else if (stack[size+1] == expr->bin.lft)
	    {
	      if (lftest){size++;stack[size] = expr->bin.rgt;}
	      else	  size--;
	    }
	  else     {size++;stack[size] = expr->bin.lft;}
	  break;
    }

  return(0);
}

int OBSmartsPattern::GetCharge(int idx)
{
  AtomExpr *expr = _pat->atom[idx].expr;

  int size=0;
  AtomExpr *stack[15];
  memset(stack,'\0',sizeof(AtomExpr*)*15);
  bool lftest=true;

  for (size=0,stack[size] = expr;size >= 0;expr=stack[size])
    {
    switch (expr->type)
      {
         case AE_LEAF:
	   switch( expr->leaf.prop )
	     {
	     case AL_NEGATIVE:  
	       return(-1*(int)expr->leaf.value);
	     case AL_POSITIVE: 
	       return((int)expr->leaf.value);
	     default: lftest=true;
	     }
	   size--;
	   break;

         case AE_OR:
         case AE_ANDHI: 
         case AE_ANDLO: 

	   if (stack[size+1] == expr->bin.rgt)      size--;
	   else if (stack[size+1] == expr->bin.lft)
	     {
	       if (lftest){size++;stack[size] = expr->bin.rgt;}
	       else	  size--;
	     }
	   else     {size++;stack[size] = expr->bin.lft;}
	   break;

         case AE_NOT: return(0);
         case AE_RECUR: return(0);
      }
    }

  return(0);
}

int OBSmartsPattern::GetAtomicNum(int idx)
{
  AtomExpr *expr = _pat->atom[idx].expr;

  int size=0;
  AtomExpr *stack[15];
  memset(stack,'\0',sizeof(AtomExpr*)*15);
  bool lftest=true;

  for (size=0,stack[size] = expr;size >= 0;expr=stack[size])
    {
    switch (expr->type)
      {
         case AE_LEAF:
	   if ( expr->leaf.prop == AL_ELEM) return(expr->leaf.value);
	   lftest = true;
	   size--;
	   break;

         case AE_OR:
         case AE_ANDHI: 
         case AE_ANDLO: 

	   if (stack[size+1] == expr->bin.rgt)      size--;
	   else if (stack[size+1] == expr->bin.lft)
	     {
	       if (lftest){size++;stack[size] = expr->bin.rgt;}
	       else	  size--;
	     }
	   else     {size++;stack[size] = expr->bin.lft;}
	   break;

         case AE_NOT:   return(0);
         case AE_RECUR: return(0);
      }
    }

  return(0);
}

void OBSmartsPattern::GetBond(int &src,int &dst,int &ord,int idx)
{
  src = _pat->bond[idx].src;
  dst = _pat->bond[idx].dst;
  ord = GetExprOrder(_pat->bond[idx].expr);
}

void SmartsLexReplace(std::string &s,std::vector<std::pair<std::string,std::string> > &vlex)
{
  size_t j,pos;
  std::string token,repstr;
  std::vector<std::pair<std::string,std::string> >::iterator i;
  
  for (pos = 0,pos = s.find("$",pos);pos < s.size();pos = s.find("$",pos))
    //for (pos = 0,pos = s.find("$",pos);pos != std::string::npos;pos = s.find("$",pos))
    {
      pos++;
      for (j = pos;j < s.size();j++)
	if (!isalpha(s[j]) && !isdigit(s[j]) && s[j] != '_')
	  break;
      if (pos == j) continue;

      token = s.substr(pos,j-pos); 
      for (i = vlex.begin();i != vlex.end();i++)
	if (token == i->first)
	  {
	    repstr = "(" + i->second + ")";
	    s.replace(pos,j-pos,repstr);
	    j = 0;
	  }
      pos = j;
    }
}

}

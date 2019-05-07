/**********************************************************************
parsmart.cpp - SMARTS parser.

Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2006 by Geoffrey R. Hutchison

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

#include <ctype.h>
#include <iomanip>
#include <cstring>

#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/parsmart.h>
#include <openbabel/stereo/stereo.h>
#include <openbabel/stereo/tetrahedral.h>

using namespace std;

namespace OpenBabel
{
  /*! \class OBSmartsPattern parsmart.h <openbabel/parsmart.h>

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

    For more information on the SMARTS support in Open Babel, see the wiki page:
    http://openbabel.org/wiki/SMARTS

    The SMARTS matcher, or OBSmartsPattern, is a separate object which can
    match patterns in the OBMol class. The following code demonstrates how
    to use the OBSmartsPattern class:
    \code
    OBMol mol;
    ...
    OBSmartsPattern sp;
    sp.Init("CC");
    sp.Match(mol);
    vector<vector<int> > maplist;
    maplist = sp.GetMapList();
    //or maplist = sp.GetUMapList();
    //print out the results
    vector<vector<int> >::iterator i;
    vector<int>::iterator j;
    for (i = maplist.begin();i != maplist.end();++i)
    {
    for (j = i->begin();j != i->end();++j)
    cout << j << ' `;
    cout << endl;
    }
    \endcode

    The preceding code reads in a molecule, initializes a SMARTS pattern
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

#define ATOMPOOL      1
#define BONDPOOL      1

#define AE_ANDHI        1
#define AE_ANDLO        2
#define AE_OR           3
#define AE_RECUR        4
#define AE_NOT          5
#define AE_TRUE         6
#define AE_FALSE        7
#define AE_AROMATIC     8
#define AE_ALIPHATIC    9
#define AE_CYCLIC       10
#define AE_ACYCLIC      11
#define AE_MASS         12
#define AE_ELEM         13
#define AE_AROMELEM     14
#define AE_ALIPHELEM    15
#define AE_HCOUNT       16
#define AE_CHARGE       17
#define AE_CONNECT      18
#define AE_DEGREE       19
#define AE_IMPLICIT     20
#define AE_RINGS        21
#define AE_SIZE         22
#define AE_VALENCE      23
#define AE_CHIRAL       24
#define AE_HYB          25
#define AE_RINGCONNECT  26

#define AL_CLOCKWISE      1
#define AL_ANTICLOCKWISE  2
#define AL_UNSPECIFIED    0

/* Each BondExpr node is identified by an integer type field
   selected from the list of BE_ codes below.  BE_ANDHI,
   BE_ANDLO and BE_OR are binary nodes of type BondExpr.bin,
   BE_NOT is unary node of type BondExpr.mon, and the remaining
   code are all leaf nodes.  */

#define BE_ANDHI        1
#define BE_ANDLO        2
#define BE_OR           3
#define BE_NOT          4
#define BE_ANY          5
#define BE_DEFAULT      6
#define BE_SINGLE       7
#define BE_DOUBLE       8
#define BE_TRIPLE       9
#define BE_QUAD         10
#define BE_AROM         11
#define BE_RING         12
#define BE_UP           13
#define BE_DOWN         14
#define BE_UPUNSPEC     15
#define BE_DOWNUNSPEC   16


  static int CreateAtom(Pattern*,AtomExpr*,int,int vb=0);

  const int SmartsImplicitRef = -9999; // Used as a placeholder when recording atom nbrs for chiral atoms


  /*=============================*/
  /*  Standard Utility Routines  */
  /*=============================*/

  static void FatalAllocationError( const char *ptr )
  {
    stringstream errorMsg;
    errorMsg << "Error: Unable to allocate" << ptr << endl;
    obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
  }

  /*================================*/
  /*  Atom Expression Manipulation  */
  /*================================*/

  static void FreePattern( Pattern* );
  static Pattern *CopyPattern( Pattern* );

  static AtomExpr *CopyAtomExpr( AtomExpr *expr )
  {
    AtomExpr *result;

    result = new AtomExpr;
    result->type = expr->type;
    switch( expr->type )
      {
      case AE_ANDHI:
      case AE_ANDLO:
      case AE_OR:
        result->bin.lft = CopyAtomExpr(expr->bin.lft);
        result->bin.rgt = CopyAtomExpr(expr->bin.rgt);
        break;

      case AE_NOT:
        result->mon.arg = CopyAtomExpr(expr->mon.arg);
        break;

      case AE_RECUR:
        result->recur.recur = CopyPattern((Pattern*)expr->recur.recur);
        break;

      default:
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
          {
          case AE_ANDHI:
          case AE_ANDLO:
          case AE_OR:
            FreeAtomExpr(expr->bin.lft);
            FreeAtomExpr(expr->bin.rgt);
            break;

          case AE_NOT:
            FreeAtomExpr(expr->mon.arg);
            break;

          case AE_RECUR:
            FreePattern((Pattern*)expr->recur.recur);
            break;
          }

        delete expr;
      }
  }

  static AtomExpr *BuildAtomPred( int type )
  {
    AtomExpr *result;

    result = new AtomExpr;
    result->leaf.type = type;
    result->leaf.value = 0;
    return result;
  }

  static AtomExpr *BuildAtomLeaf( int type, int val )
  {
    AtomExpr *result;

    result = new AtomExpr;
    result->leaf.type = type;
    result->leaf.value = val;
    return result;
  }

  static AtomExpr *BuildAtomNot( AtomExpr *expr )
  {
    AtomExpr *result;

    result = new AtomExpr;
    result->mon.type = AE_NOT;
    result->mon.arg = expr;
    return result;
  }

  static AtomExpr *BuildAtomBin( int op, AtomExpr *lft, AtomExpr *rgt )
  {
    AtomExpr *result;

    result = new AtomExpr;
    result->bin.type = op;
    result->bin.lft = lft;
    result->bin.rgt = rgt;
    return result;
  }

  static AtomExpr *BuildAtomRecurs( Pattern *pat )
  {
    AtomExpr *result;

    result = new AtomExpr;
    result->recur.type = AE_RECUR;
    result->recur.recur = (void*)pat;
    return result;
  }

  static AtomExpr *GenerateElement( int elem )
  {
    return BuildAtomLeaf(AE_ELEM,elem);
  }

  static AtomExpr *GenerateAromElem( int elem, int flag )
  {
    return flag ? BuildAtomLeaf(AE_AROMELEM,elem)
                : BuildAtomLeaf(AE_ALIPHELEM,elem);
  }

  /*================================*/
  /*  Bond Expression Manipulation  */
  /*================================*/

  static BondExpr *CopyBondExpr( BondExpr *expr )
  {
    BondExpr *result;

    result = new BondExpr;
    result->type = expr->type;
    switch( expr->type )
      {
      case BE_ANDHI:
      case BE_ANDLO:
      case BE_OR:
        result->bin.lft = CopyBondExpr(expr->bin.lft);
        result->bin.rgt = CopyBondExpr(expr->bin.rgt);
        break;

      case BE_NOT:
        result->mon.arg = CopyBondExpr(expr->mon.arg);
        break;
      }
    return result;
  }

  /**
   * Check if two BondExpr objects are the same. This is used for ring closures
   * to identify invalid SMARTS like:
   *
   *   C-1CCCCC#1
   *   C=1CCCCC:1
   *
   * However, the SMARTS below are valid and the bond expression next to the the
   * second closure digit is used.
   *
   *   C1CCCCC#1
   *   C1CCCCC=1
   */
  static bool EquivalentBondExpr( BondExpr *expr1, BondExpr *expr2 )
  {
    if (expr1 == NULL && expr2 == NULL)
      return true;
    if (expr1 == NULL && expr2 != NULL)
      return false;
    if (expr1 != NULL && expr2 == NULL)
      return false;

    if (expr1->type != expr2->type)
      return false;

    switch( expr1->type )
      {
      case BE_ANDHI:
      case BE_ANDLO:
      case BE_OR:
        return EquivalentBondExpr(expr1->bin.lft, expr2->bin.lft) &&
               EquivalentBondExpr(expr1->bin.rgt, expr2->bin.rgt);

      case BE_NOT:
        return EquivalentBondExpr(expr1->mon.arg, expr2->mon.arg);
      }
    return true;
  }

  static void FreeBondExpr( BondExpr *expr )
  {
    if( expr )
      {
        switch( expr->type )
          {
          case BE_ANDHI:
          case BE_ANDLO:
          case BE_OR:
            FreeBondExpr(expr->bin.lft);
            FreeBondExpr(expr->bin.rgt);
            break;

          case BE_NOT:
            FreeBondExpr(expr->mon.arg);
            break;
          }

        delete expr;
      }
  }

  static BondExpr *BuildBondLeaf( int type )
  {
    BondExpr *result;

    result = new BondExpr;
    result->type = type;
    return result;
  }

  static BondExpr *BuildBondNot( BondExpr *expr )
  {
    BondExpr *result;

    result = new BondExpr;
    result->mon.type = BE_NOT;
    result->mon.arg = expr;
    return result;
  }

  static BondExpr *BuildBondBin( int op, BondExpr *lft, BondExpr *rgt )
  {
    BondExpr *result;

    result = new BondExpr;
    result->bin.type = op;
    result->bin.lft = lft;
    result->bin.rgt = rgt;
    return result;
  }

  static BondExpr *GenerateDefaultBond( void )
  {
    return BuildBondLeaf(BE_DEFAULT);
  }

  /*===============================*/
  /*  SMARTS Pattern Manipulation  */
  /*===============================*/

  static Pattern *AllocPattern( void )
  {
    Pattern *ptr;

    ptr = new Pattern;
    if( !ptr ) {
      FatalAllocationError("pattern");
      return NULL;
    }

    ptr->atom = (AtomSpec*)0;
    ptr->aalloc = 0;
    ptr->acount = 0;

    ptr->bond = (BondSpec*)0;
    ptr->balloc = 0;
    ptr->bcount = 0;

    ptr->parts = 1;

    ptr->hasExplicitH=false;
    return ptr;
  }

  static int CreateAtom( Pattern *pat, AtomExpr *expr, int part,int vb)
  {
    int index,size;

    if (!pat)
      return -1; // should never happen

    if( pat->acount == pat->aalloc )
      {
        pat->aalloc += ATOMPOOL;
        size = (int)(pat->aalloc*sizeof(AtomSpec));
        if( pat->atom )
          {
            AtomSpec *tmp = new AtomSpec[pat->aalloc];
            copy(pat->atom, pat->atom + pat->aalloc - ATOMPOOL, tmp);
            delete [] pat->atom;
            pat->atom = tmp;
          }
        else
          pat->atom = new AtomSpec[pat->aalloc];
        if( !pat->atom )
          FatalAllocationError("atom pool");
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

    if (!pat)
      return -1; // should never happen

    if( pat->bcount == pat->balloc )
      {
        pat->balloc += BONDPOOL;
        size = (int)(pat->balloc*sizeof(BondSpec));
        if( pat->bond )
          {
            BondSpec *tmp = new BondSpec[pat->balloc];
            copy(pat->bond, pat->bond + pat->balloc - BONDPOOL, tmp);
            delete [] pat->bond;
            pat->bond = tmp;
          }
        else
          pat->bond = new BondSpec[pat->balloc];
        if( !pat->bond )
          FatalAllocationError("bond pool");
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
      {
        aexpr = CopyAtomExpr(pat->atom[i].expr);
        CreateAtom(result,aexpr,pat->atom[i].part);
      }

    for( i=0; i<pat->bcount; i++ )
      {
        bexpr = CopyBondExpr(pat->bond[i].expr);
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
            if (pat->atom != NULL) {
              //free(pat->atom);
              delete [] pat->atom;
              pat->atom = NULL;
            }
          }

        if( pat->balloc )
          {
            for( i=0; i<pat->bcount; i++ )
              FreeBondExpr(pat->bond[i].expr);
            if (pat->bond != NULL) {
              //free(pat->bond);
              delete [] pat->bond;
              pat->bond = NULL;
            }
          }
        delete pat;
        pat = NULL;
      }
  }

  /*=========================*/
  /*  SMARTS Syntax Parsing  */
  /*=========================*/


  Pattern *OBSmartsPattern::SMARTSError( Pattern *pat )
  {
    stringstream errorMsg;
    errorMsg << "SMARTS Error:\n" << MainPtr << endl;
    errorMsg << setw(LexPtr-MainPtr+1) << '^' << endl;
    obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError, onceOnly);

    FreePattern(pat);
    return (Pattern*)0;
  }

  AtomExpr *OBSmartsPattern::ParseSimpleAtomPrimitive( void )
  {
    switch( *LexPtr++ )
      {
      case '*':
        return BuildAtomPred(AE_TRUE);
      case 'A':
        return BuildAtomPred(AE_ALIPHATIC);
      case 'B':
        if( *LexPtr == 'r' )
          {
            LexPtr++;
            return GenerateElement(35);
          }
        return GenerateElement(5);
      case 'C':
        if( *LexPtr == 'l' )
          {
            LexPtr++;
            return GenerateElement(17);
          }
        return GenerateAromElem(6,false);
      case 'F':
        return GenerateElement( 9);
      case 'I':
        return GenerateElement(53);
      case 'N':
        return GenerateAromElem(7,false);
      case 'O':
        return GenerateAromElem(8,false);
      case 'P':
        return GenerateAromElem(15, false);
      case 'S':
        return GenerateAromElem(16,false);
      case 'a':
        if( *LexPtr == 's' )
          {
            LexPtr++;
            return GenerateAromElem(33, true);
          }
        return BuildAtomPred(AE_AROMATIC);
      case 'c':
        return GenerateAromElem( 6,true);
      case 'n':
        return GenerateAromElem( 7,true);
      case 'o':
        return GenerateAromElem( 8,true);
      case 'p':
        return GenerateAromElem(15,true);
      case 's':
        if( *LexPtr == 'e' )
          {
            LexPtr++;
            return GenerateAromElem(34, true);
          }
        return GenerateAromElem(16,true);
      }
    LexPtr--;
    return (AtomExpr*)0;
  }

  AtomExpr *OBSmartsPattern::ParseComplexAtomPrimitive( void )
  {
    Pattern *pat;
    int index;

    switch( *LexPtr++ )
      {
      case('#'):
        if( !isdigit(*LexPtr) )
          return( (AtomExpr*)0 );

        index = 0;
        while( isdigit(*LexPtr) )
          index = index*10 + ((*LexPtr++)-'0');
        if( index > 255 )
          {
            LexPtr--;
            return( (AtomExpr*)0 );
          }
        else if( !index )
          return( (AtomExpr*)0 );
        return( GenerateElement(index) );

      case('$'):
        if( *LexPtr != '(' )
          return( (AtomExpr*)0 );
        LexPtr++;
        pat = ParseSMARTSPattern();

        if( !pat )
          return( (AtomExpr*)0 );
        if( *LexPtr != ')' )
          {
            FreePattern(pat);
            return( (AtomExpr*)0 );
          }
        LexPtr++;
        return( BuildAtomRecurs(pat) );

      case('*'):
        return BuildAtomPred(AE_TRUE);

      case('+'):
        if( isdigit(*LexPtr) )
          {
            index = 0;
            while( isdigit(*LexPtr) )
              index = index*10 + ((*LexPtr++)-'0');
          }
        else
          {
            index = 1;
            while( *LexPtr == '+' )
              {
                LexPtr++;
                index++;
              }
          }
        return BuildAtomLeaf(AE_CHARGE,index);

      case('-'):
        if( isdigit(*LexPtr) )
          {
            index = 0;
            while( isdigit(*LexPtr) )
              index = index*10 + ((*LexPtr++)-'0');
          }
        else
          {
            index = 1;
            while( *LexPtr == '-' )
              {
                LexPtr++;
                index++;
              }
          }
        return BuildAtomLeaf(AE_CHARGE,-index);

      case '@':
        if (*LexPtr == '?')
          {
            LexPtr++;
            return BuildAtomLeaf(AE_CHIRAL,AL_UNSPECIFIED); // unspecified
          }
        else if (*LexPtr != '@')
          return BuildAtomLeaf(AE_CHIRAL,AL_ANTICLOCKWISE);
        else
          {
            LexPtr++;
            return BuildAtomLeaf(AE_CHIRAL,AL_CLOCKWISE);
          }

      case '^':
        if (isdigit(*LexPtr))
          {
            index = 0;
            while( isdigit(*LexPtr) )
              index = index*10 + ((*LexPtr++)-'0');
            return BuildAtomLeaf(AE_HYB,index);
          }
        else
          return BuildAtomLeaf(AE_HYB,1);

      case('0'): case('1'): case('2'): case('3'): case('4'):
      case('5'): case('6'): case('7'): case('8'): case('9'):
        index = LexPtr[-1]-'0';
        while( isdigit(*LexPtr) )
          index = index*10 + ((*LexPtr++)-'0');
        return BuildAtomLeaf(AE_MASS,index);

      case('A'):
        switch( *LexPtr++ )
          {
          case('c'):  return GenerateElement(89);
          case('g'):  return GenerateElement(47);
          case('l'):  return GenerateElement(13);
          case('m'):  return GenerateElement(95);
          case('r'):  return GenerateElement(18);
          case('s'):  return GenerateElement(33);
          case('t'):  return GenerateElement(85);
          case('u'):  return GenerateElement(79);
          }
        LexPtr--;
        return BuildAtomPred(AE_ALIPHATIC);

      case('B'):
        switch( *LexPtr++ )
          {
          case('a'):  return GenerateElement(56);
          case('e'):  return GenerateElement( 4);
          case('i'):  return GenerateElement(83);
          case('k'):  return GenerateElement(97);
          case('r'):  return GenerateElement(35);
          }
        LexPtr--;
        return GenerateElement(5);

      case('C'):
        switch( *LexPtr++ )
          {
          case('a'):  return GenerateElement(20);
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
        return GenerateAromElem(6,false);

      case('D'):
        if( *LexPtr == 'y' )
          {
            LexPtr++;
            return GenerateElement(66);
          }
        else if( isdigit(*LexPtr) )
          {
            index = 0;
            while( isdigit(*LexPtr) )
              index = index*10 + ((*LexPtr++)-'0');
            return BuildAtomLeaf(AE_DEGREE,index);
          }
        return BuildAtomLeaf(AE_DEGREE,1);

      case('E'):
        if( *LexPtr == 'r' )
          {
            LexPtr++;
            return GenerateElement(68);
          }
        else if( *LexPtr == 's' )
          {
            LexPtr++;
            return GenerateElement(99);
          }
        else if( *LexPtr == 'u' )
          {
            LexPtr++;
            return GenerateElement(63);
          }
        break;

      case('F'):
        if( *LexPtr == 'e' )
          {
            LexPtr++;
            return GenerateElement(26);
          }
        else if( *LexPtr == 'm' )
          {
            LexPtr++;
            return GenerateElement(100);
          }
        else if( *LexPtr == 'r' )
          {
            LexPtr++;
            return GenerateElement(87);
          }
        return GenerateElement(9);

      case('G'):
        if( *LexPtr == 'a' )
          {
            LexPtr++;
            return( GenerateElement(31) );
          }
        else if( *LexPtr == 'd' )
          {
            LexPtr++;
            return( GenerateElement(64) );
          }
        else if( *LexPtr == 'e' )
          {
            LexPtr++;
            return( GenerateElement(32) );
          }
        break;

      case('H'):
        if( *LexPtr == 'e' )
          {
            LexPtr++;
            return( GenerateElement( 2) );
          }
        else if( *LexPtr == 'f' )
          {
            LexPtr++;
            return( GenerateElement(72) );
          }
        else if( *LexPtr == 'g' )
          {
            LexPtr++;
            return( GenerateElement(80) );
          }
        else if( *LexPtr == 'o' )
          {
            LexPtr++;
            return( GenerateElement(67) );
          }
        else if( isdigit(*LexPtr) )
          {
            index = 0;
            while( isdigit(*LexPtr) )
              index = index*10 + ((*LexPtr++)-'0');
            return BuildAtomLeaf(AE_HCOUNT,index);
          }
        return BuildAtomLeaf(AE_HCOUNT,1);

      case('I'):
        if( *LexPtr == 'n' )
          {
            LexPtr++;
            return( GenerateElement(49) );
          }
        else if( *LexPtr == 'r' )
          {
            LexPtr++;
            return( GenerateElement(77) );
          }
        return( GenerateElement(53) );

      case('K'):
        if( *LexPtr == 'r' )
          {
            LexPtr++;
            return( GenerateElement(36) );
          }
        return( GenerateElement(19) );

      case('L'):
        if( *LexPtr == 'a' )
          {
            LexPtr++;
            return( GenerateElement( 57) );
          }
        else if( *LexPtr == 'i' )
          {
            LexPtr++;
            return( GenerateElement(  3) );
          }
        else if( *LexPtr == 'r' )
          {
            LexPtr++;
            return( GenerateElement(103) );
          }
        else if( *LexPtr == 'u' )
          {
            LexPtr++;
            return( GenerateElement( 71) );
          }
        break;

      case('M'):
        if( *LexPtr == 'd' )
          {
            LexPtr++;
            return( GenerateElement(101) );
          }
        else if( *LexPtr == 'g' )
          {
            LexPtr++;
            return( GenerateElement( 12) );
          }
        else if( *LexPtr == 'n' )
          {
            LexPtr++;
            return( GenerateElement( 25) );
          }
        else if( *LexPtr == 'o' )
          {
            LexPtr++;
            return( GenerateElement( 42) );
          }
        break;

      case('N'):
        switch( *LexPtr++ )
          {
          case('a'):  return( GenerateElement( 11) );
          case('b'):  return( GenerateElement( 41) );
          case('d'):  return( GenerateElement( 60) );
          case('e'):  return( GenerateElement( 10) );
          case('i'):  return( GenerateElement( 28) );
          case('o'):  return( GenerateElement(102) );
          case('p'):  return( GenerateElement( 93) );
          }
        LexPtr--;
        return( GenerateAromElem(7,false) );

      case('O'):
        if( *LexPtr == 's' )
          {
            LexPtr++;
            return( GenerateElement(76) );
          }
        return( GenerateAromElem(8,false) );

      case('P'):
        switch( *LexPtr++ )
          {
          case('a'):  return( GenerateElement(91) );
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

      case('R'):
        switch( *LexPtr++ )
          {
          case('a'):  return( GenerateElement(88) );
          case('b'):  return( GenerateElement(37) );
          case('e'):  return( GenerateElement(75) );
          case('h'):  return( GenerateElement(45) );
          case('n'):  return( GenerateElement(86) );
          case('u'):  return( GenerateElement(44) );
          }
        LexPtr--;
        if( isdigit(*LexPtr) )
          {
            index = 0;
            while( isdigit(*LexPtr) )
              index = index*10 + ((*LexPtr++)-'0');
            if( index == 0 )
              return BuildAtomPred(AE_ACYCLIC);
            return BuildAtomLeaf(AE_RINGS,index);
          }
        return BuildAtomPred(AE_CYCLIC);

      case('S'):
        switch( *LexPtr++ )
          {
          case('b'):  return( GenerateElement(51) );
          case('c'):  return( GenerateElement(21) );
          case('e'):  return( GenerateElement(34) );
          case('i'):  return( GenerateElement(14) );
          case('m'):  return( GenerateElement(62) );
          case('n'):  return( GenerateElement(50) );
          case('r'):  return( GenerateElement(38) );
          }
        LexPtr--;
        return( GenerateAromElem(16,false) );

      case('T'):
        switch( *LexPtr++ )
          {
          case('a'):  return( GenerateElement(73) );
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

      case('X'):
        if( *LexPtr == 'e' )
          {
            LexPtr++;
            return( GenerateElement(54) );
          }
        else if( isdigit(*LexPtr) )
          {
            index = 0;
            while( isdigit(*LexPtr) )
              index = index*10 + ((*LexPtr++)-'0');
            if (index == 0) // default to 1 (if no number present)
              index = 1;
            return BuildAtomLeaf(AE_CONNECT,index);
          }
        return BuildAtomLeaf(AE_CONNECT,1);

      case('Y'):
        if( *LexPtr == 'b' )
          {
            LexPtr++;
            return( GenerateElement(70) );
          }
        return( GenerateElement(39) );

      case('Z'):
        if( *LexPtr == 'n' )
          {
            LexPtr++;
            return GenerateElement(30);
          }
        else if( *LexPtr == 'r' )
          {
            LexPtr++;
            return GenerateElement(40);
          }
        break;

      case('a'):
        if( *LexPtr == 's' )
          {
            LexPtr++;
            return GenerateAromElem(33,true);
          }
        return BuildAtomPred(AE_AROMATIC);

      case('c'):
        return GenerateAromElem(6,true);

      case('h'):
        if( isdigit(*LexPtr) )
          {
            index = 0;
            while( isdigit(*LexPtr) )
              index = index*10 + ((*LexPtr++)-'0');
          }
        else
          index = 1;
        return BuildAtomLeaf(AE_IMPLICIT,index);

      case('n'):  return GenerateAromElem(7,true);
      case('o'):  return GenerateAromElem(8,true);
      case('p'):  return GenerateAromElem(15,true);

      case('r'):
        if( isdigit(*LexPtr) )
          {
            index = 0;
            while( isdigit(*LexPtr) )
              index = index*10 + ((*LexPtr++)-'0');
            if( index == 0 )
              return BuildAtomPred(AE_ACYCLIC);
            return BuildAtomLeaf(AE_SIZE,index);
          }
        return BuildAtomPred(AE_CYCLIC);

      case('s'):
        if( *LexPtr == 'e' )
          {
            LexPtr++;
            return GenerateAromElem(34,true);
          }
        return GenerateAromElem(16,true);

      case('v'):
        if( isdigit(*LexPtr) )
          {
            index = 0;
            while( isdigit(*LexPtr) )
              index = index*10 + ((*LexPtr++)-'0');
            return BuildAtomLeaf(AE_VALENCE,index);
          }
        return BuildAtomLeaf(AE_VALENCE,1);

      case('x'):
        if( isdigit(*LexPtr) )
          {
            index = 0;
            while( isdigit(*LexPtr) )
              index = index*10 + ((*LexPtr++)-'0');
            return BuildAtomLeaf(AE_RINGCONNECT,index);
          }
        return BuildAtomPred(AE_CYCLIC);
      }
    LexPtr--;
    return (AtomExpr*)0;
  }

  AtomExpr *OBSmartsPattern::ParseAtomExpr( int level )
  {
    AtomExpr *expr1 = NULL;
    AtomExpr *expr2 = NULL;
    char *prev;

    switch( level )
      {
      case(0): /* Low Precedence Conjunction */
        if( !(expr1=ParseAtomExpr(1)) )
          return (AtomExpr*)0;

        while( *LexPtr == ';' )
          {
            LexPtr++;
            if( !(expr2=ParseAtomExpr(1)) )
              {
                FreeAtomExpr(expr1);
                return (AtomExpr*)0;
              }
            expr1 = BuildAtomBin(AE_ANDLO,expr1,expr2);
          }
        return expr1;

      case(1): /* Disjunction */
        if( !(expr1=ParseAtomExpr(2)) )
          return (AtomExpr*)0;

        while( *LexPtr == ',' )
          {
            LexPtr++;
            if( !(expr2=ParseAtomExpr(2)) )
              {
                FreeAtomExpr(expr1);
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
          {
            if( *LexPtr=='&' )
              LexPtr++;
            prev = LexPtr;
            if( !(expr2=ParseAtomExpr(3)) )
              {
                if( prev != LexPtr )
                  {
                    FreeAtomExpr(expr1);
                    return( (AtomExpr*)0 );
                  }
                else
                  return( expr1 );
              }
            expr1 = BuildAtomBin(AE_ANDHI,expr1,expr2);
          }
        return( expr1 );

      case(3): /* Negation or Primitive */
        if( *LexPtr == '!' )
          {
            LexPtr++;
            if( !(expr1=ParseAtomExpr(3)) )
              return( (AtomExpr*)0 );
            return( BuildAtomNot(expr1) );
          }
        return( ParseComplexAtomPrimitive() );
      }
    return (AtomExpr*)0;
  }

  BondExpr *OBSmartsPattern::ParseBondPrimitive( void )
  {
    char bsym = *LexPtr++;

    switch(bsym)
      {
      case '-':  return BuildBondLeaf(BE_SINGLE);
      case '=':  return BuildBondLeaf(BE_DOUBLE);
      case '#':  return BuildBondLeaf(BE_TRIPLE);
      case '$':  return BuildBondLeaf(BE_QUAD);
      case ':':  return BuildBondLeaf(BE_AROM);
      case '@':  return BuildBondLeaf(BE_RING);
      case '~':  return BuildBondLeaf(BE_ANY);

      // return BuildBondLeaf(*LexPtr == '?' ? BE_UPUNSPEC : BE_UP);
      case '/':  return BuildBondLeaf(BE_SINGLE);
      // return BuildBondLeaf(*LexPtr == '?' ? BE_DOWNUNSPEC : BE_DOWN);
      case '\\': return BuildBondLeaf(BE_SINGLE);
      }
    LexPtr--;
    return (BondExpr*)0;
  }

  BondExpr *OBSmartsPattern::ParseBondExpr( int level )
  {
    BondExpr *expr1 = NULL;
    BondExpr *expr2 = NULL;
    char *prev;

    switch( level )
      {
      case(0): /* Low Precedence Conjunction */
        if( !(expr1=ParseBondExpr(1)) )
          return (BondExpr*)0;

        while( *LexPtr == ';' )
          {
            LexPtr++;
            if( !(expr2=ParseBondExpr(1)) )
              {
                FreeBondExpr(expr1);
                return (BondExpr*)0;
              }
            expr1 = BuildBondBin(BE_ANDLO,expr1,expr2);
          }
        return expr1;

      case(1): /* Disjunction */
        if( !(expr1=ParseBondExpr(2)) )
          return (BondExpr*)0;

        while( *LexPtr == ',' )
          {
            LexPtr++;
            if( !(expr2=ParseBondExpr(2)) )
              {
                FreeBondExpr(expr1);
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
          {
            if( *LexPtr == '&' )
              LexPtr++;
            prev = LexPtr;
            if( !(expr2=ParseBondExpr(3)) )
              {
                if( prev != LexPtr )
                  {
                    FreeBondExpr(expr1);
                    return (BondExpr*)0;
                  }
                else
                  return expr1;
              }
            expr1 = BuildBondBin(BE_ANDHI,expr1,expr2);
          }
        return expr1;

      case(3): /* Negation or Primitive */
        if( *LexPtr == '!' )
          {
            LexPtr++;
            if( !(expr1=ParseBondExpr(3)) )
              return (BondExpr*)0;
            return BuildBondNot(expr1);
          }
        return ParseBondPrimitive();
      }
    return (BondExpr*)0;
  }

  int OBSmartsPattern::GetVectorBinding()
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

  Pattern *OBSmartsPattern::ParseSMARTSError( Pattern *pat, BondExpr *expr )
  {
    if( expr )
      FreeBondExpr(expr);
    return SMARTSError(pat);
  }

  Pattern *OBSmartsPattern::SMARTSParser( Pattern *pat, ParseState *stat,
                                int prev, int part )
  {
    int vb = 0;
    AtomExpr *aexpr;
    BondExpr *bexpr;
    int index;

    bexpr = (BondExpr*)0;

    while( *LexPtr )
      {
        switch( *LexPtr++ )
          {
          case('.'):
            return ParseSMARTSError(pat,bexpr);

          case('-'):  case('='):  case('#'): case('$'):
          case(':'):  case('~'):  case('@'):
          case('/'):  case('\\'): case('!'):
            LexPtr--;
            if( (prev==-1) || bexpr )
              return ParseSMARTSError(pat,bexpr);
            if( !(bexpr=ParseBondExpr(0)) )
              return ParseSMARTSError(pat,bexpr);
            break;

          case('('):
            if( bexpr )
              {
                LexPtr--;
                return ParseSMARTSError(pat,bexpr);
              }
            if( prev == -1 )
              {
                index = pat->acount;
                pat = SMARTSParser(pat,stat,-1,part);
                if( !pat )
                  return( (Pattern*)0 );
                if( index == pat->acount )
                  return ParseSMARTSError(pat,bexpr);
                prev = index;
              }
            else
              {
                pat = SMARTSParser(pat,stat,prev,part);
                if( !pat )
                  return (Pattern*)0;
              }

            if( *LexPtr != ')' )
              return ParseSMARTSError(pat,bexpr);
            LexPtr++;
            break;

          case(')'):  LexPtr--;
            if( (prev==-1) || bexpr )
              return ParseSMARTSError(pat,bexpr);
            return pat;

          case('%'):  if( prev == -1 )
              {
                LexPtr--;
                return ParseSMARTSError(pat,bexpr);
              }

            if( isdigit(LexPtr[0]) && isdigit(LexPtr[1]) )
              {
                index = 10*(LexPtr[0]-'0') + (LexPtr[1]-'0');
                LexPtr += 2;
              }
            else
              return ParseSMARTSError(pat,bexpr);

            if( stat->closure[index] == -1 )
              {
                stat->closord[index] = bexpr;
                stat->closure[index] = prev;
              }
            else if( stat->closure[index] != prev )
              {
                if( !bexpr ) {
                  if (!stat->closord[index]) {
                    bexpr = GenerateDefaultBond();
                    FreeBondExpr(stat->closord[index]);
                  } else
                    bexpr = stat->closord[index];
                } else if (stat->closord[index] && !EquivalentBondExpr(bexpr, stat->closord[index]))
                  return ParseSMARTSError(pat,bexpr);

                CreateBond(pat,bexpr,prev,stat->closure[index]);
                stat->closure[index] = -1;
                bexpr = (BondExpr*)0;
              }
            else
              return ParseSMARTSError(pat,bexpr);
            break;

          case('0'):  case('1'):  case('2'):
          case('3'):  case('4'):  case('5'):
          case('6'):  case('7'):  case('8'):
          case('9'):  LexPtr--;
            if( prev == -1 )
              return ParseSMARTSError(pat,bexpr);
            index = (*LexPtr++)-'0';

            if( stat->closure[index] == -1 )
              { // Ring opening
                stat->closord[index] = bexpr;
                stat->closure[index] = prev;
                pat->atom[prev].nbrs.push_back(-index); // Store the BC idx as a -ve
                bexpr = (BondExpr*)0;
              }
            else if( stat->closure[index] != prev )
              { // Ring closure
                if( !bexpr ) {
                  if (!stat->closord[index]) {
                    bexpr = GenerateDefaultBond();
                    FreeBondExpr(stat->closord[index]);
                  } else
                    bexpr = stat->closord[index];
                } else if (stat->closord[index] && !EquivalentBondExpr(bexpr, stat->closord[index]))
                  return ParseSMARTSError(pat,bexpr);

                CreateBond(pat,bexpr,prev,stat->closure[index]);
                pat->atom[prev].nbrs.push_back(stat->closure[index]);
                for (unsigned int nbr_idx=0; nbr_idx < pat->atom[stat->closure[index]].nbrs.size(); ++nbr_idx) {
                  if (pat->atom[stat->closure[index]].nbrs[nbr_idx] == -index)
                    pat->atom[stat->closure[index]].nbrs[nbr_idx] = prev;
                }
                stat->closure[index] = -1;
                bexpr = (BondExpr*)0;
              }
            else
              return ParseSMARTSError(pat,bexpr);
            break;

          case('['):
            // shortcut for '[H]' primitive (PR#1463791)
            if (*LexPtr == 'H' && *(LexPtr+1) == ']')
              {
                aexpr = GenerateElement(1);
                LexPtr++; // skip the 'H'
                pat->hasExplicitH = true;
              }
            else
              aexpr = ParseAtomExpr(0);
            vb = (*LexPtr == ':') ? GetVectorBinding():0;
            if( !aexpr || (*LexPtr!=']') )
              return ParseSMARTSError(pat,bexpr);
            index = CreateAtom(pat,aexpr,part,vb);
            if( prev != -1 )
              {
                if( !bexpr )
                  bexpr = GenerateDefaultBond();
                CreateBond(pat,bexpr,prev,index);
                pat->atom[index].nbrs.push_back(prev);
                pat->atom[prev].nbrs.push_back(index);
                bexpr = (BondExpr*)0;
              }
            if(*(LexPtr-1) == 'H' && *(LexPtr-2) == '@' ) { // i.e. [C@H] or [C@@H]
              pat->atom[index].nbrs.push_back(SmartsImplicitRef);
            }
            prev = index;
            LexPtr++;
            break;

          default:
            LexPtr--;
            aexpr = ParseSimpleAtomPrimitive();
            if( !aexpr )
              return ParseSMARTSError(pat,bexpr);
            index = CreateAtom(pat,aexpr,part);
            if( prev != -1 )
              {
                if( !bexpr )
                  bexpr = GenerateDefaultBond();
                CreateBond(pat,bexpr,prev,index);
                pat->atom[index].nbrs.push_back(prev);
                pat->atom[prev].nbrs.push_back(index);
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

    for (i = 0;i < pat->bcount;++i)
      {
        pat->bond[i].grow = (bv[pat->bond[i].src] && bv[pat->bond[i].dst])?
          false:true;

        bv.SetBitOn(pat->bond[i].src);
        bv.SetBitOn(pat->bond[i].dst);
      }
  }

  static int GetChiralFlag(AtomExpr *expr)
  {
    int tmp1,tmp2;

    switch (expr->type)
      {
        case AE_CHIRAL:
          return expr->leaf.value;

        case AE_ANDHI:
        case AE_ANDLO:
          tmp1 = GetChiralFlag(expr->bin.lft);
          tmp2 = GetChiralFlag(expr->bin.rgt);
          if (tmp1 == 0) return tmp2;
          if (tmp2 == 0) return tmp1;
          if (tmp1 == tmp2) return tmp1;
          break;

        case AE_OR:
          tmp1 = GetChiralFlag(expr->bin.lft);
          tmp2 = GetChiralFlag(expr->bin.rgt);
          if (tmp1 == 0 || tmp2 == 0) return 0;
          if (tmp1 == tmp2) return tmp1;
          break;

        case AE_NOT:
          // Treat [!@] as [@@], and [!@@] as [@]
          tmp1 = GetChiralFlag(expr->mon.arg);
          if (tmp1 == AL_ANTICLOCKWISE) return AL_CLOCKWISE;
          if (tmp1 == AL_CLOCKWISE) return AL_ANTICLOCKWISE;
          break;
      }

    return 0;
  }

  Pattern *OBSmartsPattern::ParseSMARTSPart( Pattern *result, int part )
  {
    ParseState stat;
    int i,flag;

    for( i=0; i<100; i++ )
      stat.closure[i] = -1;

    result = SMARTSParser(result,&stat,-1,part);

    flag = false;
    for( i=0; i<100; i++ )
      if( stat.closure[i] != -1 )
        {
          FreeBondExpr(stat.closord[i]);
          flag = true;
        }

    if( result )
      {
        if( flag )
          return(SMARTSError(result));
        else
          {
            MarkGrowBonds(result);
            result->ischiral = false;
            for (i = 0;i < result->acount;++i)
              {
                result->atom[i].chiral_flag = GetChiralFlag(result->atom[i].expr);
                if (result->atom[i].chiral_flag)
                  result->ischiral = true;
              }
            return(result);
          }
      }
    else
      return (Pattern*)0;
  }


  Pattern *OBSmartsPattern::ParseSMARTSPattern( void )
  {
    Pattern *result;
    result = AllocPattern();

    while( *LexPtr == '(' )
      {
        if (!result)
          return (Pattern*)0; // ensure we don't get a null dereference

        LexPtr++;
        result = ParseSMARTSPart(result,result->parts);
        if( !result )
          return (Pattern*)0;
        result->parts++;

        if( *LexPtr != ')' )
          return SMARTSError(result);
        LexPtr++;

        if( !*LexPtr || (*LexPtr==')') )
          return result;

        if( *LexPtr != '.' )
          return SMARTSError(result);

        // Here's where we'd handle fragments
        //        cerr << " conjunction " << LexPtr[0] << endl;
        LexPtr++;
      }

    return ParseSMARTSPart(result,0);
  }

  Pattern *OBSmartsPattern::ParseSMARTSString( char *ptr )
  {
    Pattern *result;

    if( !ptr || !*ptr )
      return (Pattern*)0;

    LexPtr = MainPtr = ptr;
    result = ParseSMARTSPattern();
    if( result && *LexPtr )
      return SMARTSError(result);
    return result;
  }

  Pattern *OBSmartsPattern::ParseSMARTSRecord( char *ptr )
  {
    char *src;

    src = ptr;
    while( *src && !isspace(*src) )
      src++;

    if( isspace(*src) )
      {
        *src++ = '\0';
        while( isspace(*src) )
          src++;
      }

    return ParseSMARTSString(ptr);
  }

  //**********************************
  //********Pattern Matching**********
  //**********************************

  bool OBSmartsPattern::Init(const char *buffer)
  {
    if (_buffer != NULL)
      delete[] _buffer;
    _buffer = new char[strlen(buffer) + 1];
    strcpy(_buffer,buffer);

    if (_pat != NULL)
      FreePattern(_pat);
    _pat = ParseSMARTSRecord(_buffer);
    _str = _buffer;

    return(_pat != (Pattern*)NULL);
  }

  bool OBSmartsPattern::Init(const std::string &s)
  {
    if (_buffer != NULL)
      delete[] _buffer;
    _buffer = new char[s.length() + 1];
    strcpy(_buffer, s.c_str());

    if (_pat != NULL)
      FreePattern(_pat);
    _pat = ParseSMARTSRecord(_buffer);
    _str = s;

    return (_pat != (Pattern*) NULL);
  }

  OBSmartsPattern::~OBSmartsPattern()
  {
    if (_pat)
      FreePattern(_pat);
    if(_buffer)
    	delete [] _buffer;
  }

  bool OBSmartsPattern::Match(OBMol &mol,bool single)
  {
	OBSmartsMatcher matcher;
	if(_pat == NULL)
      return false;
    if(_pat->hasExplicitH) //The SMARTS pattern contains [H]
      {
        //Do matching on a copy of mol with explicit hydrogens
        OBMol tmol = mol;
        tmol.AddHydrogens(false,false);
        return(matcher.match(tmol,_pat,_mlist,single));
      }
    return(matcher.match(mol,_pat,_mlist,single));
  }

  bool OBSmartsPattern::HasMatch(OBMol &mol) const
  {
	  //a convenience function
	  std::vector<std::vector<int> > dummy;
	  return Match(mol, dummy, Single);
  }

  bool OBSmartsPattern::Match(OBMol &mol, std::vector<std::vector<int> > & mlist,
		  MatchType mtype /*=All*/) const
  {
	OBSmartsMatcher matcher;
	mlist.clear();
	if(_pat == NULL)
      return false;
    if(_pat->hasExplicitH) //The SMARTS pattern contains [H]
      {
        //Do matching on a copy of mol with explicit hydrogens
        OBMol tmol = mol;
        tmol.AddHydrogens(false,false);
        if(!matcher.match(tmol,_pat,mlist,mtype == Single))
        	return false;
      }
    else if(!matcher.match(mol,_pat,mlist,mtype == Single))
    	return false;

    if((mtype == AllUnique) && mlist.size() > 1)
    {
    	//uniquify
         bool ok;
        OBBitVec bv;
        std::vector<OBBitVec> vbv;
        std::vector<std::vector<int> > ulist;
        std::vector<std::vector<int> >::iterator i;
        std::vector<OBBitVec>::iterator j;

        for (i = mlist.begin();i != mlist.end();++i)
          {
            ok = true;
            bv.Clear();
            bv.FromVecInt(*i);
            for (j = vbv.begin();j != vbv.end() && ok;++j)
              if ((*j) == bv)
                ok = false;

            if (ok)
              {
                ulist.push_back(*i);
                vbv.push_back(bv);
              }
          }

        mlist = ulist;
    }
    return true;
  }


  bool OBSmartsPattern::RestrictedMatch(OBMol &mol,
                                        std::vector<std::pair<int,int> > &pr,
                                        bool single)
  {
    bool ok;
    std::vector<std::vector<int> > mlist;
    std::vector<std::vector<int> >::iterator i;
    std::vector<std::pair<int,int> >::iterator j;

    OBSmartsMatcher matcher;
    matcher.match(mol,_pat,mlist);
    _mlist.clear();
    if (mlist.empty())
      return(false);

    for (i = mlist.begin();i != mlist.end();++i)
      {
        ok = true;
        for (j = pr.begin();j != pr.end() && ok;++j)
          if ((*i)[j->first] != j->second)
            ok = false;

        if (ok)
          _mlist.push_back(*i);
        if (single && !_mlist.empty())
          return(true);
      }

    return((_mlist.empty()) ? false:true);
  }

  bool OBSmartsPattern::RestrictedMatch(OBMol &mol,OBBitVec &vres, bool single)
  {
    bool ok;
    std::vector<int>::iterator j;
    std::vector<std::vector<int> > mlist;
    std::vector<std::vector<int> >::iterator i;

    OBSmartsMatcher matcher;
    matcher.match(mol,_pat,mlist);

    _mlist.clear();
    if (mlist.empty())
      return(false);

    for (i = mlist.begin();i != mlist.end();++i)
      {
        ok = true;
        for (j = i->begin();j != i->end();++j)
          if (!vres[*j])
            {
              ok = false;
              break;
            }
        if (!ok)
          continue;

        _mlist.push_back(*i);
        if (single && !_mlist.empty())
          return(true);
      }

    return((_mlist.empty()) ? false:true);
  }

  void OBSmartsMatcher::SetupAtomMatchTable(std::vector<std::vector<bool> > &ttab,
                           const Pattern *pat, OBMol &mol)
  {
    int i;

    ttab.resize(pat->acount);
    for (i = 0;i < pat->acount;++i)
      ttab[i].resize(mol.NumAtoms()+1);

    OBAtom *atom;
    std::vector<OBAtom*>::iterator j;
    for (i = 0;i < pat->acount;++i)
      for (atom = mol.BeginAtom(j);atom;atom = mol.NextAtom(j))
        if (EvalAtomExpr(pat->atom[0].expr,atom))
          ttab[i][atom->GetIdx()] = true;
  }

  void OBSmartsMatcher::FastSingleMatch(OBMol &mol, const Pattern *pat,
                              std::vector<std::vector<int> > &mlist)
  {
    OBAtom *atom,*a1,*nbr;
    std::vector<OBAtom*>::iterator i;

    OBBitVec bv(mol.NumAtoms()+1);
    std::vector<int> map;
    map.resize(pat->acount);
    std::vector<std::vector<OBBond*>::iterator> vi;
    std::vector<bool> vif;

    if (pat->bcount)
      {
        vif.resize(pat->bcount);
        vi.resize(pat->bcount);
      }

    int bcount;
    for (atom = mol.BeginAtom(i);atom;atom=mol.NextAtom(i))
      if (EvalAtomExpr(pat->atom[0].expr,atom))
        {
          map[0] = atom->GetIdx();
          if (pat->bcount)
            vif[0] = false;
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


  bool OBSmartsMatcher::match(OBMol &mol, const Pattern *pat,
                    std::vector<std::vector<int> > &mlist,bool single)
  {
    mlist.clear();
    if (!pat || pat->acount == 0)
      return(false);//shouldn't ever happen

    if (single && !pat->ischiral) {
      // perform a fast single match (only works for non-chiral SMARTS)
      FastSingleMatch(mol,pat,mlist);
    } else {
      // perform normal match (chirality ignored and checked below)
      OBSSMatch ssm(mol,pat);
      ssm.Match(mlist);
    }

    if (pat->ischiral) {
      std::vector<std::vector<int> >::iterator m;
      std::vector<std::vector<int> > tmpmlist;

      tmpmlist.clear();
      // iterate over the atom mappings
      for (m = mlist.begin();m != mlist.end();++m) {

        bool allStereoCentersMatch = true;

        // for each pattern atom
        for (int j = 0; j < pat->acount; ++j) {
          // skip non-chiral pattern atoms
          if (!pat->atom[j].chiral_flag)
            continue;
          // ignore @? in smarts, parse like any other smarts
          if (pat->atom[j].chiral_flag == AL_UNSPECIFIED)
            continue;

          // use the mapping the get the chiral atom in the molecule being queried
          OBAtom *center = mol.GetAtom((*m)[j]);

          // get the OBTetrahedralStereo::Config from the molecule
          OBStereoFacade stereo(&mol);
          OBTetrahedralStereo *ts = stereo.GetTetrahedralStereo(center->GetId());
          if (!ts || !ts->GetConfig().specified) {
            // no stereochemistry specified in molecule for the atom
            // corresponding to the chiral pattern atom using the current
            // mapping --> no match
            allStereoCentersMatch = false;
            break;
          }

          std::vector<int> nbrs = pat->atom[j].nbrs;

          if (nbrs.size() != 4) { // 3 nbrs currently not supported. Other values are errors.
            //stringstream ss;
            //ss << "Ignoring stereochemistry. There are " << nbrs.size() << " connections to this atom instead of 4. Title: " << mol.GetTitle();
            //obErrorLog.ThrowError(__FUNCTION__, ss.str(), obWarning);
            continue;
          }

          // construct a OBTetrahedralStereo::Config using the smarts pattern
          OBTetrahedralStereo::Config smartsConfig;
          smartsConfig.center = center->GetId();
          if (nbrs.at(0) == SmartsImplicitRef)
            smartsConfig.from = OBStereo::ImplicitRef;
          else
            smartsConfig.from = mol.GetAtom( (*m)[nbrs.at(0)] )->GetId();
          OBStereo::Ref firstref;
          if (nbrs.at(1) == SmartsImplicitRef)
            firstref = OBStereo::ImplicitRef;
          else
            firstref = mol.GetAtom( (*m)[nbrs.at(1)] )->GetId();
          OBAtom *ra2 = mol.GetAtom( (*m)[nbrs.at(2)] );
          OBAtom *ra3 = mol.GetAtom( (*m)[nbrs.at(3)] );
          smartsConfig.refs = OBStereo::MakeRefs(firstref, ra2->GetId(), ra3->GetId());

          smartsConfig.view = OBStereo::ViewFrom;
          switch (pat->atom[j].chiral_flag) {
            case AL_CLOCKWISE:
              smartsConfig.winding = OBStereo::Clockwise;
              break;
            case AL_ANTICLOCKWISE:
              smartsConfig.winding = OBStereo::AntiClockwise;
              break;
            default:
              smartsConfig.specified = false;
          }

          // cout << "smarts config = " << smartsConfig << endl;
          // cout << "molecule config = " << ts->GetConfig() << endl;
          // cout << "match = " << (ts->GetConfig() == smartsConfig) << endl;

          // and save the match if the two configurations are the same
          if (ts->GetConfig() != smartsConfig)
            allStereoCentersMatch = false;

          // don't waste time checking more stereocenters using this mapping if one didn't match
          if (!allStereoCentersMatch)
            break;
        }

        // if all the atoms in the molecule match the stereochemistry specified
        // in the smarts pattern, save this mapping as a match
        if (allStereoCentersMatch)
          tmpmlist.push_back(*m);
      }

      mlist = tmpmlist;
    }

    return(!mlist.empty());
  }

  bool OBSmartsMatcher::EvalAtomExpr(AtomExpr *expr,OBAtom *atom)
  {
    for (;;)
      switch (expr->type)
        {
        case AE_TRUE:
          return true;
        case AE_FALSE:
          return false;
        case AE_AROMATIC:
          return atom->IsAromatic();
        case AE_ALIPHATIC:
          return !atom->IsAromatic();
        case AE_CYCLIC:
          return atom->IsInRing();
        case AE_ACYCLIC:
          return !atom->IsInRing();

        case AE_MASS:
          return expr->leaf.value == atom->GetIsotope();
        case AE_ELEM:
          return expr->leaf.value == (int)atom->GetAtomicNum();
        case AE_AROMELEM:
          return expr->leaf.value == (int)atom->GetAtomicNum() &&
                 atom->IsAromatic();
        case AE_ALIPHELEM:
          return expr->leaf.value == (int)atom->GetAtomicNum() &&
                 !atom->IsAromatic();
        case AE_HCOUNT:
          return expr->leaf.value == ((int)atom->ExplicitHydrogenCount() +
                                      (int)atom->GetImplicitHCount());
        case AE_CHARGE:
          return expr->leaf.value == atom->GetFormalCharge();
        case AE_CONNECT:
          return expr->leaf.value == (int)atom->GetTotalDegree();
        case AE_DEGREE:
          return expr->leaf.value == (int)atom->GetExplicitDegree();
        case AE_IMPLICIT:
          return expr->leaf.value == (int)atom->GetImplicitHCount();
        case AE_RINGS:
          return expr->leaf.value == (int)atom->MemberOfRingCount();
        case AE_SIZE:
          return atom->IsInRingSize(expr->leaf.value);
        case AE_VALENCE:
          return expr->leaf.value == (int)atom->GetTotalValence();
        case AE_CHIRAL:
          // always return true (i.e. accept the match) and check later
          return true;
        case AE_HYB:
          return expr->leaf.value == (int)atom->GetHyb();
        case AE_RINGCONNECT:
          return expr->leaf.value == (int)atom->CountRingBonds();

        case AE_NOT:
          return !EvalAtomExpr(expr->mon.arg,atom);
        case AE_ANDHI: /* Same as AE_ANDLO */
        case AE_ANDLO:
          if (!EvalAtomExpr(expr->bin.lft,atom))
            return false;
          expr = expr->bin.rgt;
          break;
        case AE_OR:
          if (EvalAtomExpr(expr->bin.lft,atom))
            return true;
          expr = expr->bin.rgt;
          break;

        case AE_RECUR:
          {
            //see if pattern has been matched
            std::vector<std::pair<const Pattern*,std::vector<bool> > >::iterator i;
            for (i = RSCACHE.begin();i != RSCACHE.end();++i)
              if (i->first == (Pattern*)expr->recur.recur)
                return(i->second[atom->GetIdx()]);

            //perceive and match pattern
            std::vector<std::vector<int> >::iterator j;
            std::vector<bool> vb(((OBMol*) atom->GetParent())->NumAtoms()+1);
            std::vector<std::vector<int> > mlist;
            if (match( *((OBMol *) atom->GetParent()),
                       (Pattern*)expr->recur.recur,mlist))
              for (j = mlist.begin();j != mlist.end();++j)
                vb[(*j)[0]] = true;

            RSCACHE.push_back(std::pair<const Pattern*,
                              std::vector<bool> > ((const Pattern*)expr->recur.recur,
                                                   vb));

            return(vb[atom->GetIdx()]);
          }

        default:
          return false;
        }
  }

  bool OBSmartsMatcher::EvalBondExpr(BondExpr *expr,OBBond *bond)
  {
    for (;;)
      switch( expr->type )
        {
        case BE_ANDHI:
        case BE_ANDLO:
          if (!EvalBondExpr(expr->bin.lft,bond))
            return false;
          expr = expr->bin.rgt;
          break;

        case BE_OR:
          if (EvalBondExpr(expr->bin.lft,bond))
            return true;
          expr = expr->bin.rgt;
          break;

        case BE_NOT:
          return !EvalBondExpr(expr->mon.arg,bond);

        case BE_ANY:
          return true;
        case BE_DEFAULT:
          return bond->GetBondOrder()==1 || bond->IsAromatic();
        case BE_SINGLE:
          return bond->GetBondOrder()==1 && !bond->IsAromatic();
        case BE_DOUBLE:
          return bond->GetBondOrder()==2 && !bond->IsAromatic();
        case BE_TRIPLE:
          return bond->GetBondOrder() == 3;
        case BE_QUAD:
          return bond->GetBondOrder() == 4;
        case BE_AROM:
          return bond->IsAromatic();
        case BE_RING:
          return bond->IsInRing();
        //case BE_UP:
        //  return bond->IsUp();
        //case BE_DOWN:
        //  return bond->IsDown();
        //case BE_UPUNSPEC: // up or unspecified (i.e., not down)
        //  return !bond->IsDown();
        //case BE_DOWNUNSPEC: // down or unspecified (i.e., not up)
        //  return !bond->IsUp();
        default:
          return false;
        }
  }

  std::vector<std::vector<int> > &OBSmartsPattern::GetUMapList()
  {
    if (_mlist.empty() || _mlist.size() == 1)
      return(_mlist);

    bool ok;
    OBBitVec bv;
    std::vector<OBBitVec> vbv;
    std::vector<std::vector<int> > mlist;
    std::vector<std::vector<int> >::iterator i;
    std::vector<OBBitVec>::iterator j;

    for (i = _mlist.begin();i != _mlist.end();++i)
      {
        ok = true;
        bv.Clear();
        bv.FromVecInt(*i);
        for (j = vbv.begin();j != vbv.end() && ok;++j)
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

    for ( i = _mlist.begin() ; i != _mlist.end() ; ++i )
      {
        for (j = (*i).begin();j != (*i).end();++j)
          ofs << *j << ' ' << ends;
        ofs << endl;
      }
  }

  //*******************************************************************
  //  The OBSSMatch class performs exhaustive matching using recursion
  //  Explicit stack handling is used to find just a single match in
  //  match()
  //*******************************************************************

  OBSSMatch::OBSSMatch(OBMol &mol, const Pattern *pat)
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
    if (_uatoms)
      delete [] _uatoms;
  }

  void OBSSMatch::Match(std::vector<std::vector<int> > &mlist,int bidx)
  {
	  OBSmartsMatcher matcher;
    if (bidx == -1)
      {
        OBAtom *atom;
        std::vector<OBAtom*>::iterator i;
        for (atom = _mol->BeginAtom(i);atom;atom = _mol->NextAtom(i))
          if (matcher.EvalAtomExpr(_pat->atom[0].expr,atom))
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
        int src = _pat->bond[bidx].src;
        int dst = _pat->bond[bidx].dst;

        if (_map[src] <= 0 || _map[src] > (signed)_mol->NumAtoms())
          return;

        AtomExpr *aexpr = _pat->atom[dst].expr;
        BondExpr *bexpr = _pat->bond[bidx].expr;
        OBAtom *atom,*nbr;
        std::vector<OBBond*>::iterator i;

        atom = _mol->GetAtom(_map[src]);
        for (nbr = atom->BeginNbrAtom(i);nbr;nbr = atom->NextNbrAtom(i))
          if (!_uatoms[nbr->GetIdx()] && matcher.EvalAtomExpr(aexpr,nbr) &&
        		  matcher.EvalBondExpr(bexpr,((OBBond*) *i)))
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
        if (bond && matcher.EvalBondExpr(_pat->bond[bidx].expr,bond))
          Match(mlist,bidx+1);
      }
  }

  static int GetExprOrder(BondExpr *expr)
  {
    int tmp1,tmp2;

    switch( expr->type )
      {
      case BE_SINGLE:  return 1;
      case BE_DOUBLE:  return 2;
      case BE_TRIPLE:  return 3;
      case BE_QUAD:    return 4;
      case BE_AROM:    return 5;

      case BE_UP:
      case BE_DOWN:
      case BE_UPUNSPEC:
      case BE_DOWNUNSPEC:
        return 1;
      
      case BE_ANDHI:
      case BE_ANDLO:
        tmp1 = GetExprOrder(expr->bin.lft);
        tmp2 = GetExprOrder(expr->bin.rgt);
        if (tmp1 == 0) return tmp2;
        if (tmp2 == 0) return tmp1;
        if (tmp1 == tmp2) return tmp1;
        break;

      case BE_OR:
        tmp1 = GetExprOrder(expr->bin.lft);
        if (tmp1 == 0) return 0;
        tmp2 = GetExprOrder(expr->bin.rgt);
        if (tmp2 == 0) return 0;
        if (tmp1 == tmp2) return tmp1;
        break;
      }

    return 0;
  }

  static int GetExprCharge(AtomExpr *expr)
  {
    int tmp1,tmp2;

    switch( expr->type )
      {
      case AE_CHARGE:
        return expr->leaf.value;

      case AE_ANDHI:
      case AE_ANDLO:
        tmp1 = GetExprCharge(expr->bin.lft);
        tmp2 = GetExprCharge(expr->bin.rgt);
        if (tmp1 == 0) return tmp2;
        if (tmp2 == 0) return tmp1;
        if (tmp1 == tmp2) return tmp1;
        break;

      case AE_OR:
        tmp1 = GetExprCharge(expr->bin.lft);
        if (tmp1 == 0) return 0;
        tmp2 = GetExprCharge(expr->bin.rgt);
        if (tmp2 == 0) return 0;
        if (tmp1 == tmp2) return tmp1;
        break;
      }

    return 0;
  }

  int OBSmartsPattern::GetCharge(int idx)
  {
    return GetExprCharge(_pat->atom[idx].expr);
  }

  static int GetExprAtomicNum(AtomExpr *expr)
  {
    int tmp1,tmp2;

    switch( expr->type )
      {
      case AE_ELEM:
      case AE_AROMELEM:
      case AE_ALIPHELEM:
        return expr->leaf.value;

      case AE_ANDHI:
      case AE_ANDLO:
        tmp1 = GetExprAtomicNum(expr->bin.lft);
        tmp2 = GetExprAtomicNum(expr->bin.rgt);
        if (tmp1 == 0) return tmp2;
        if (tmp2 == 0) return tmp1;
        if (tmp1 == tmp2) return tmp1;
        break;

      case AE_OR:
        tmp1 = GetExprAtomicNum(expr->bin.lft);
        if (tmp1 == 0) return 0;
        tmp2 = GetExprAtomicNum(expr->bin.rgt);
        if (tmp2 == 0) return 0;
        if (tmp1 == tmp2) return tmp1;
        break;
      }

    return 0;
  }

  int OBSmartsPattern::GetAtomicNum(int idx)
  {
    return GetExprAtomicNum(_pat->atom[idx].expr);
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
        for (j = pos;j < s.size();++j)
          if (!isalpha(s[j]) && !isdigit(s[j]) && s[j] != '_')
            break;
        if (pos == j)
          continue;

        token = s.substr(pos,j-pos);
        for (i = vlex.begin();i != vlex.end();++i)
          if (token == i->first)
            {
              repstr = "(" + i->second + ")";
              s.replace(pos,j-pos,repstr);
              j = 0;
            }
        pos = j;
      }
  }

} // end namespace OpenBabel

//! \file parsmart.cpp
//! \brief Implementation of Daylight SMARTS parser.

/**********************************************************************
parsmart.h - Daylight SMARTS parser.
 
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2005 by Geoffrey R. Hutchison
 
This file is part of the Open Babel project.
For more information, see <http://openbabel.sourceforge.net/>
 
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#ifndef OB_PARSMART_H
#define OB_PARSMART_H

#include <string>
#include <vector>

#include <openbabel/mol.h>

/*==========================*/
/*  SMARTS Data Structures  */
/*==========================*/

namespace OpenBabel
{

  // mark this so that SWIG will not attempt to wrap for scripting languages

#ifndef SWIG

  //! \union _AtomExpr parsmart.h <openbabel/parsmart.h>
  //! \brief An internal (SMARTS parser) atomic expression
  typedef union _AtomExpr {
    int type;
    struct
    {
      int type;
      int prop;
      int value;
    }
      leaf;
    struct
    {
      int type;
      void *recur;
    }
      recur;
    struct
    {
      int type;
      union _AtomExpr *arg;
    }
      mon;
    struct
    {
      int type;
      union _AtomExpr *lft;
      union _AtomExpr *rgt;
    }
      bin;
  } AtomExpr;

#define BE_LEAF      0x01
#define BE_ANDHI     0x02
#define BE_ANDLO     0x03
#define BE_NOT       0x04
#define BE_OR        0x05

#define BL_CONST     0x01
#define BL_TYPE      0x02

#define BT_SINGLE     0x01
#define BT_DOUBLE     0x02
#define BT_TRIPLE     0x03
#define BT_AROM       0x04
#define BT_UP         0x05
#define BT_DOWN       0x06
#define BT_UPUNSPEC   0x07
#define BT_DOWNUNSPEC 0x08
#define BT_RING       0x09

  //! \union _BondExpr parsmart.h <openbabel/parsmart.h>
  //! \brief An internal (SMARTS parser) bond expression
  typedef union _BondExpr {
    int type;
    struct
    {
      int type;
      int prop;
      int value;
    }
      leaf;
    struct
    {
      int type;
      union _BondExpr *arg;
    }
      mon;
    struct
    {
      int type;
      union _BondExpr *lft;
      union _BondExpr *rgt;
    }
      bin;
  } BondExpr;

  //! \struct BondSpec parsmart.h <openbabel/parsmart.h>
  //! \brief An internal (SMARTS parser) bond specification
  typedef struct
  {
    BondExpr *expr;
    int src,dst;
    int visit;
    bool grow;
  }
  BondSpec;

  //! \struct AtomSpec parsmart.h <openbabel/parsmart.h>
  //! \brief An internal (SMARTS parser) atom specification
  typedef struct
  {
    AtomExpr *expr;
    int visit;
    int part;
    int chiral_flag;
    int vb;
  }
  AtomSpec;

  //! \struct Pattern parsmart.h <openbabel/parsmart.h>
  //! \brief A SMARTS parser internal pattern
  typedef struct
  {
    int aalloc,acount;
    int balloc,bcount;
    bool ischiral;
    AtomSpec *atom;
    BondSpec *bond;
    int parts;
  }
  Pattern;
#else
  // for SWIG, just forward declare that we have some Pattern struct
  // (but this is private and not wrapped for scripting languages)
  struct Pattern;
#endif

  // class introduction in parsmart.cpp
  //! \brief SMARTS (SMiles ARbitrary Target Specification) substructure searching
  class OBAPI OBSmartsPattern
  {
  protected:
    std::vector<bool>          		  _growbond; //!< \deprecated (Not used)
    std::vector<std::vector<int> >	_mlist;    //!< The list of matches
    Pattern                        *_pat;      //!< The parsed SMARTS pattern
    std::string				              _str;      //!< The string of the SMARTS expression

  public:
    OBSmartsPattern() : _pat(NULL) { }
    virtual ~OBSmartsPattern();

  OBSmartsPattern(const OBSmartsPattern& cp): _pat(NULL)
      {
        *this = cp;
      }
    OBSmartsPattern& operator=(const OBSmartsPattern& cp)
      {
        if (_pat)
          delete [] _pat;
        _pat = NULL;
        std::string s = cp._str;
        Init(s);
        return (*this);
      }
    
    //! \name Initialization Methods
    //@{
    //! Parse the @p pattern SMARTS string.
    //! \return Whether the pattern is a valid SMARTS expression
    bool         Init(const char* pattern);
    //! Parse the @p pattern SMARTS string.
    //! \return Whether the pattern is a valid SMARTS expression
    bool         Init(const std::string& pattern);
    //@}

    //! \name Pattern Properties
    //@{
    //! \return the SMARTS string which is currently used
    const std::string &GetSMARTS() const    {      return _str;    }
    //! \return the SMARTS string which is currently used
    std::string  &GetSMARTS()               {      return _str;    }

    //! \return If the SMARTS pattern is an empty expression (e.g., invalid)
    bool         Empty() const     {      return(_pat == NULL);    }
    //! \return If the SMARTS pattern is a valid expression
    bool         IsValid() const   {      return(_pat != NULL);    }

    //! \return the number of atoms in the SMARTS pattern
    unsigned int NumAtoms()   const
    {
      return _pat ? _pat->acount : 0;
    }
    //! \return the number of bonds in the SMARTS pattern
    unsigned int NumBonds()   const
    {
      return _pat ? _pat->bcount : 0;
    }

    //! Access the bond @p idx in the internal pattern
    //! \param src The index of the beginning atom
    //! \param dst The index of the end atom
    //! \param ord The bond order of this bond
    //! \param idx The index of the bond in the SMARTS pattern
    void         GetBond(int& src,int& dst,int& ord,int idx);
    //! \return the atomic number of the atom @p idx in the internal pattern
    int          GetAtomicNum(int idx);
    //! \return the formal charge of the atom @p idx in the internal pattern
    int          GetCharge(int idx);

    //! \return the vector binding of the atom @p idx in the internal pattern
    int          GetVectorBinding(int idx) const
    {
      return(_pat->atom[idx].vb);
    }
    //@}

    //! \name Matching methods (SMARTS on a specific OBMol)
    //@{
    //! Perform SMARTS matching for the pattern specified using Init().
    //! \param mol The molecule to use for matching
    //! \param single Whether only a single match is required (faster). Default is false.
    //! \return Whether matches occurred
    bool Match(OBMol &mol, bool single=false);

    bool RestrictedMatch(OBMol &mol, std::vector<std::pair<int,int> > &pairs, bool single=false);

    bool RestrictedMatch(OBMol &mol, OBBitVec &bv, bool single=false);
    //! \return the number of non-unique SMARTS matches 
    //! To get the number of unique SMARTS matches, query GetUMapList()->size()
    unsigned int NumMatches() const
    {
      return (unsigned int)_mlist.size();
    }

    //! \return the entire list of non-unique matches for this pattern
    //! \see GetUMapList()
    std::vector<std::vector<int> > &GetMapList()
      {
        return(_mlist);
      }
    //! \return An iterator over the (non-unique) match list, starting at the beginning
    std::vector<std::vector<int> >::iterator BeginMList()
      {
        return(_mlist.begin());
      }
    //! \return An iterator over the non-unique match list, set to the end
    std::vector<std::vector<int> >::iterator EndMList()
      {
        return(_mlist.end());
      }

    //! \return the entire list of unique matches for this pattern
    /**
        A unique match is defined as one which does not cover the 
        identical atoms that a previous match has covered.
        
        For instance, the pattern [OD1]~C~[OD1] describes a
        carboxylate group. This pattern will match both atom number
        permutations of the carboxylate, and if GetMapList() is called, both
        matches will be returned. If GetUMapList() is called only unique
        matches of the pattern will be returned.
    **/
    std::vector<std::vector<int> > &GetUMapList();
    //@}

    //! Debugging -- write a list of matches to the output stream
    void         WriteMapList(std::ostream&);
  };

  //! \class OBSSMatch parsmart.h <openbabel/parsmart.h>
  //! \brief Internal class: performs fast, exhaustive matching used to find 
  //! just a single match in match() using recursion and explicit stack handling.
  class OBAPI OBSSMatch
  {
  protected:
    bool        *_uatoms;
    OBMol       *_mol;
    Pattern     *_pat;
    std::vector<int>  _map;

  public:
    OBSSMatch(OBMol&,Pattern*);
    ~OBSSMatch();
    void Match(std::vector<std::vector<int> > &v, int bidx=-1);
  };

  OBAPI void SmartsLexReplace(std::string &,
                              std::vector<std::pair<std::string,std::string> > &);

} // end namespace OpenBabel

#endif // OB_PARSMART_H

//! \file parsmart.h
//! \brief Daylight SMARTS parser.

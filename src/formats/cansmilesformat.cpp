/**********************************************************************
Some portions Copyright (C) 2005-2006 by Craig A. James, eMolecules Inc.
Some portions Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2005 by Geoffrey R. Hutchison
Some portions Copyright (C) 2004 by Chris Morley

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/


/****************************************************************************
 *
 * This is derived from "smilesformat.cpp", the SMILES reader/writer for
 * OpenBabel 2.0.  This implements the CANSMI datatype, which is a write-only
 * format (the regular SMILES parser does the reading).
 * 
 ****************************************************************************/

/*----------------------------------------------------------------------
 * PLEASE READ THIS
 *
 * I originally intended to make this a standard file-format plugin using
 * the new OpenBabel 2.x dynamic linking system.  However, it proved to 
 * be too difficult to debug the code and get the thing to work, so 
 * this is a "fallback" to static linking.  The macro OPENBABEL_FORMAT
 * is used to zap all of the OpenBabel 2.x stuff and revert to an 
 * ordinary function call to generate the canonical SMILES.  If, in 
 * the future, we figure out how to work with the OpenBabel plug-in
 * system, the code is all here (but untested).  -- CJ.
 ----------------------------------------------------------------------*/

#define DEBUG 0

#define OPENBABEL_FORMAT 1


#include <openbabel/babelconfig.h>
#include <openbabel/base.h>
#include <openbabel/mol.h>

#if OPENBABEL_FORMAT
#include <openbabel/obconversion.h>
#include <openbabel/obmolecformat.h>
#endif

#include <openbabel/canon.h>

#if OPENBABEL_FORMAT
#else
#include "cansmilesformat.h"
#endif

using namespace std;

namespace OpenBabel
{

/*----------------------------------------------------------------------
 * CLASS: OBBondClosureInfo: For recording bond-closure digits as
 * work progresses on canonical SMILES.
 ----------------------------------------------------------------------*/

class OBBondClosureInfo
{
 public:
  OBAtom *toatom;       // second atom in SMILES order
  OBAtom *fromatom;     // first atom in SMILES order
  OBBond *bond;
  int    ringdigit;
  int    is_open;       // TRUE if SMILES processing hasn't reached 'toatom' yet

  OBBondClosureInfo(OBAtom *, OBAtom*, OBBond*, int, bool);
  ~OBBondClosureInfo();
};

OBBondClosureInfo::OBBondClosureInfo(OBAtom *a1, OBAtom *a2, OBBond *b, int rd, bool open)
{
  toatom    = a1;
  fromatom  = a2;
  bond      = b;
  ringdigit = rd;
  is_open   = open;
}

OBBondClosureInfo::~OBBondClosureInfo()
{
}


/*----------------------------------------------------------------------
 * CLASS: OBCanSmiNode: A Tree structure, each node of which is an atom in
 * the tree being built to write out the SMILES.
 ----------------------------------------------------------------------*/

class OBCanSmiNode
{
  OBAtom *_atom,*_parent;
  std::vector<OBCanSmiNode*> _child_nodes;
  std::vector<OBBond*> _child_bonds;

public:
  OBCanSmiNode(OBAtom *atom);
  ~OBCanSmiNode();

  int Size()
  {
    return(_child_nodes.empty() ? 0 : _child_nodes.size());
  }

  void SetParent(OBAtom *a)
  {
    _parent = a;
  }

  void AddChildNode(OBCanSmiNode*,OBBond*);

  OBAtom *GetAtom()
  {
    return(_atom);
  }

  OBAtom *GetParent()
  {
    return(_parent);
  }

  OBAtom *GetChildAtom(int i)
  {
    return(_child_nodes[i]->GetAtom());
  }

  OBBond *GetChildBond(int i)
  {
    return(_child_bonds[i]);
  }

  OBCanSmiNode *GetChildNode(int i)
  {
    return(_child_nodes[i]);
  }
};


OBCanSmiNode::OBCanSmiNode(OBAtom *atom)
{
  _atom = atom;
  _parent = NULL;
  _child_nodes.clear();
  _child_bonds.clear();
}

void OBCanSmiNode::AddChildNode(OBCanSmiNode *node,OBBond *bond)
{
  _child_nodes.push_back(node);
  _child_bonds.push_back(bond);
}

OBCanSmiNode::~OBCanSmiNode()
{
  vector<OBCanSmiNode*>::iterator i;
  for (i = _child_nodes.begin();i != _child_nodes.end();i++)
    delete (*i);
}

/*----------------------------------------------------------------------
* CLASS OBMol2Cansmi - Declarations
----------------------------------------------------------------------*/

class OBMol2Cansmi
{
  std::vector<int> _atmorder;
  std::vector<bool> _aromNH;
  OBBitVec _uatoms,_ubonds;
  std::vector<OBBondClosureInfo> _vopen;
#if OPENBABEL_FORMAT
  OBConversion* _pconv;
#endif

public:
  OBMol2Cansmi()
  {
  }
  ~OBMol2Cansmi() {}

#if OPENBABEL_FORMAT
  void         Init(OBConversion* pconv=NULL);
#else
  void         Init();
#endif

  void         AssignCisTrans(OBMol*);
  void         AddHydrogenToChiralCenters(OBMol &mol, OBBitVec &frag_atoms);
  bool         AtomIsChiral(OBAtom *atom);
  bool         BuildCanonTree(OBMol &mol, OBBitVec &frag_atoms,
                              vector<unsigned int> &canonical_order,
                              OBCanSmiNode *node);
  void         CorrectAromaticAmineCharge(OBMol&);
  void         CreateFragCansmiString(OBMol&, OBBitVec&, char *);
  bool         GetChiralStereo(OBCanSmiNode*,
                               vector<OBAtom*>&chiral_neighbors,
                                vector<unsigned int> &symmetry_classes,
                               char*);
  bool         GetSmilesElement(OBCanSmiNode*,
                                vector<OBAtom*>&chiral_neighbors,
                                vector<unsigned int> &symmetry_classes,
                                char*);
  int          GetSmilesValence(OBAtom *atom);
  int          GetUnusedIndex();
  vector<OBBondClosureInfo>
               GetCanonClosureDigits(OBAtom *atom,
                                     OBBitVec &frag_atoms,
                                     vector<unsigned int> &canonical_order);
  bool         IsSuppressedHydrogen(OBAtom *atom);
  bool         SameChirality(vector<OBAtom*> &v1, vector<OBAtom*> &v2);
  void         ToCansmilesString(OBCanSmiNode *node,
                                 char *buffer,
                                 OBBitVec &frag_atoms,
                                 vector<unsigned int> &symmetry_classes,
                                 vector<unsigned int> &canonical_order);
  
  std::vector<int> &GetOutputOrder()
  {
    return(_atmorder);
  }
};


/*----------------------------------------------------------------------
* CLASS OBMol2Cansmi - implementation
----------------------------------------------------------------------*/

/***************************************************************************
* FUNCTION: Init
*
* DESCRIPTION:
*       Initializes the OBMol2Cansmi writer object.
***************************************************************************/

#if OPENBABEL_FORMAT
void OBMol2Cansmi::Init(OBConversion* pconv)
#else
void OBMol2Cansmi::Init()
#endif
{
  _atmorder.clear();
  _aromNH.clear();
  _uatoms.Clear();
  _ubonds.Clear();
  _vopen.clear();
#if OPENBABEL_FORMAT
  _pconv = pconv;
#endif
}


/***************************************************************************
* FUNCTION: GetUnusedIndex
*
* DESCRIPTION:
*       Returns the next available bond-closure index for a SMILES.
*
*       You could just do this sequentially, not reusing bond-closure
*       digits, thus:
*
*               c1cc2ccccc2cc1          napthalene
*               c1ccccc1c2ccccc2        biphenyl
*
*       But molecules with more than ten rings, this requires the use of
*       two-digit ring closures (like c1ccccc1C...c%11ccccc%11).  To help
*       avoid digit reuse, this finds the lowest digit that's not currently
*       "open", thus
*
*               c1cc2ccccc2cc1          napthalene (same)
*               c1ccccc1c1ccccc1        biphenyl (reuses "1")
*
***************************************************************************/


int OBMol2Cansmi::GetUnusedIndex()
{
  int idx=1;

  vector<OBBondClosureInfo>::iterator j;
  for (j = _vopen.begin();j != _vopen.end();)
    if (j->ringdigit == idx)
      {
        idx++; //increment idx and start over if digit is already used
        j = _vopen.begin();
      }
    else j++;

  return(idx);
}

/***************************************************************************
* FUNCTION: CorrectAromaticAmineCharge
*
* DESCRIPTION:
*       Finds all aromatic nitrogens, and updates the _aromNH vector to
*       note which aromatic nitrogens require an H to balance the charge.
***************************************************************************/

void OBMol2Cansmi::CorrectAromaticAmineCharge(OBMol &mol)
{
  OBAtom *atom;
  vector<OBNodeBase*>::iterator i;

  _aromNH.clear();
  _aromNH.resize(mol.NumAtoms()+1);

  for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
    if (atom->IsNitrogen() && atom->IsAromatic())
      if (atom->GetHvyValence() == 2)
        {
          if (atom->GetValence() == 3 || atom->GetImplicitValence() == 3)
            _aromNH[atom->GetIdx()] = true;
        }
}

/***************************************************************************
* FUNCTION: OBMol2Cansmi
*
* DESCRIPTION:
*       Traverse the tree searching for acyclic olefins. If an olefin is
*       found with at least one heavy atom attachment on each end, assign
*       stereochemistry.
***************************************************************************/

void OBMol2Cansmi::AssignCisTrans(OBMol *pmol)
{
  OBBond *bond;
  vector<OBEdgeBase*>::iterator j, k;

  FOR_BONDS_OF_MOL(dbi, pmol) {

    bond = &(*dbi);

    // Not double, or in a ring?  Skip it.
    if (bond->GetBO() != 2 || bond->IsInRing())
      continue;

    OBAtom *b = bond->GetBeginAtom();
    OBAtom *c = bond->GetEndAtom();

    // skip allenes
    if (b->GetHyb() == 1 || c->GetHyb() == 1)
      continue;

    // Skip if only hydrogen on either end (Note that GetHvyValence()
    // is counting the atom across the double bond, too, so the atom
    // must have at least two heavy atoms, i.e. at most one hydrogen.)
    if (b->GetHvyValence() < 2 || c->GetHvyValence() < 2)
      continue;

    // Ok, looks like a cis/trans double bond.

    OBAtom *a,*d;

    // Look for bond with assigned stereo as in poly-ene
    for (a = b->BeginNbrAtom(j); a; a = b->NextNbrAtom(j))
      if (((OBBond*)*j)->IsUp() || ((OBBond*)*j)->IsDown())
        break;

    if (!a) {
      for (a = b->BeginNbrAtom(j);a;a = b->NextNbrAtom(j))
        if (a != c && !a->IsHydrogen())
          break;
    }
    for (d = c->BeginNbrAtom(k);d;d = c->NextNbrAtom(k)) {
      if (d != b && !d->IsHydrogen())
        break;
    }

    if (((OBBond*)*j)->IsUp() || ((OBBond*)*j)->IsDown()) { //stereo already assigned
      if (fabs(CalcTorsionAngle(a->GetVector(),b->GetVector(), c->GetVector(),d->GetVector())) > 10.0) {
        if (((OBBond*)*j)->IsUp()) {
          ((OBBond*)*k)->SetUp();
        } else {
          ((OBBond*)*k)->SetDown();
        }
      }
      else {
        if (((OBBond*)*j)->IsUp()) {
          ((OBBond*)*k)->SetDown();
        } else {
          ((OBBond*)*k)->SetUp();
        }
      }
    }
    else {      //assign stereo to both ends
      ((OBBond*)*j)->SetUp();
      if (fabs(CalcTorsionAngle(a->GetVector(),b->GetVector(), c->GetVector(),d->GetVector())) > 10.0) {
        ((OBBond*)*k)->SetUp();
      } else {
        ((OBBond*)*k)->SetDown();
      }
    }
  }
}

/***************************************************************************
* FUNCTION: GetSmilesElement
*
* DESCRIPTION:
*       Writes the symbol for an atom, e.g. "C" or "[NH2]" or "[C@@H]".
*
* RETURNS: true (always)
***************************************************************************/


bool OBMol2Cansmi::GetSmilesElement(OBCanSmiNode *node,
                                    vector<OBAtom*>&chiral_neighbors,
                                    vector<unsigned int> &symmetry_classes,
                                    char *buffer)
{
  char symbol[10];
  bool bracketElement = false;
  bool normalValence = true;

  OBAtom *atom = node->GetAtom();

  int bosum = atom->KBOSum();

  switch (atom->GetAtomicNum()) {
    case 0: break;
    case 5: break;
    case 6: break;
    case 7:
      if (atom->IsAromatic() && atom->GetHvyValence() == 2 && atom->GetImplicitValence() == 3) {
        bracketElement = !(normalValence = false);
        break;
      }
      else
        bracketElement = !(normalValence = (bosum == 3 || bosum == 5));
      break;
    case 8: break;
    case 9: break;
    case 15: break;
    case 16:
      bracketElement = !(normalValence = (bosum == 2 || bosum == 4 || bosum == 6));
      break;
    case 17: break;
    case 35: break;
    case 53: break;
    default: bracketElement = true;
  }

  if (atom->GetFormalCharge() != 0) //bracket charged elements
    bracketElement = true;

  if(atom->GetIsotope())
    bracketElement = true;

  char stereo[5] = "";
  if (GetSmilesValence(atom) > 2 && atom->IsChiral()) {
    if (GetChiralStereo(node, chiral_neighbors, symmetry_classes, stereo))
      strcat(buffer,stereo);
  }
  if (stereo[0] != '\0') 
    bracketElement = true;


#if OPENBABEL_FORMAT
  if (atom->GetSpinMultiplicity()) {
    //For radicals output bracket form anyway unless r option specified
    if(!(_pconv && _pconv->IsOption ("r")))
      bracketElement = true;
  }
#endif

  if (!bracketElement) {

    // Ordinary non-bracket element
    if (atom->GetAtomicNum()) {
      strcpy(symbol,etab.GetSymbol(atom->GetAtomicNum()));
      if (atom->IsAromatic())
        symbol[0] = tolower(symbol[0]);

#if OPENBABEL_FORMAT
      //Radical centres lc if r option set
      if(atom->GetSpinMultiplicity() && _pconv && _pconv->IsOption ("r"))
        symbol[0] = tolower(symbol[0]);
#endif
    }

    // Atomic number zero - either '*' or an external atom
    else {
      bool external = false;
      vector<pair<int,pair<OBAtom *,OBBond *> > > *externalBonds =
        (vector<pair<int,pair<OBAtom *,OBBond *> > > *)((OBMol*)atom->GetParent())->GetData("extBonds");
      vector<pair<int,pair<OBAtom *,OBBond *> > >::iterator externalBond;

      if (externalBonds)
        for(externalBond = externalBonds->begin();externalBond != externalBonds->end();externalBond++) {
          if (externalBond->second.first == atom) {
            external = true;
            strcpy(symbol,"&");
            OBBond *bond = externalBond->second.second;
            if (bond->IsUp()) {
              if ( (bond->GetBeginAtom())->HasDoubleBond() ||
                   (bond->GetEndAtom())->HasDoubleBond() )
                strcat(symbol,"\\");
            }
            if (bond->IsDown()) {
              if ( (bond->GetBeginAtom())->HasDoubleBond() ||
                   (bond->GetEndAtom())->HasDoubleBond() )
                strcat(symbol,"/");
            }
            if (bond->GetBO() == 2 && !bond->IsAromatic())
              strcat(symbol,"=");
            if (bond->GetBO() == 2 && bond->IsAromatic())
              strcat(symbol,":");
            if (bond->GetBO() == 3)
              strcat(symbol,"#");
            sprintf(symbol,"%s%d",symbol,externalBond->first);
            break;
          }
        }

      if(!external)
        strcpy(symbol,"*");
    }

    strcpy(buffer, symbol);
    return(true);
  }

  // Bracketed elements, e.g. [Pb], [OH-], [C@]

  strcpy(buffer, "[");
  if (atom->GetIsotope()) {
    char iso[4];
    sprintf(iso,"%d",atom->GetIsotope());
    strcat(buffer,iso);
  }
  if (!atom->GetAtomicNum())
    strcpy(symbol,"*");
  else {
    strcpy(symbol,etab.GetSymbol(atom->GetAtomicNum()));
    if (atom->IsAromatic())
      symbol[0] = tolower(symbol[0]);
  }
  strcat(buffer,symbol);

  // If chiral, append '@' or '@@'
  if (stereo[0] != '\0')
    strcat(buffer, stereo);

  // Add extra hydrogens
  if (!atom->IsHydrogen()) {      
    int hcount = atom->ImplicitHydrogenCount() + atom->ExplicitHydrogenCount();
    if (hcount != 0) {
      strcat(buffer,"H");
      if (hcount > 1) {
        char tcount[10];
        sprintf(tcount,"%d", hcount);
        strcat(buffer,tcount);
      }
    }
  }

  // Append charge to the end
  if (atom->GetFormalCharge() != 0) {
    if (atom->GetFormalCharge() > 0)
      strcat(buffer,"+");
    else
      strcat(buffer,"-");

    if (abs(atom->GetFormalCharge()) > 1)
      sprintf(buffer+strlen(buffer), "%d", abs(atom->GetFormalCharge()));
  }

  strcat(buffer,"]");

  return(true);
}

/***************************************************************************
* FUNCTION: OBMol2Cansmi::SameChirality
*
* DESCRIPTION:
*       Given two atom vectors representing the chiral configuration around
*       an atom, returns true/false indicating that they represent the same
*       or opposite forms.
*
*       This is used when canonicalizing a SMILES.  During canonicalization,
*       the order in which the atoms are printed is often changed, and we
*       need to compare "before and after" to see if we've altered the
*       chirality.
*
*               (NOTE: This should be integrated with OBChiralData, but
*               isn't yet because that is a bigger project that requires
*               rewriting this same section of smilesformat.cpp, as well as
*               other code that uses the OBChiralData object.)
*
*       Throughout this code, we represent chirality as an ordered set of
*       neighbor atoms, as follows.  Call the neighbors A, B, C, and D, and the
*       central (chiral) atom X.  If the SMILES is A[X@](B)(C)D, then the
*       vector could contain (A, B, C, D), in that order.  When "writing"
*       down these vectors, we ALWAYS write them in anti-clockwise order,
*       and we leave out the center, chiral atoms X.
*
*       However, there are many possible ways to write each chiral center 
*       (hence the complexity of this function).  For example, the following
*       all represent the exact same chirality:
*
*               A[X@](B)(C)D            "looking" from A to X
*               B[X@](A)(D)C            "looking" from B to X
*               C[X@](A)(B)D            "looking" from C to X
*               D[X@](A)(C)B            "looking" from D to X
*
*       Furthermore, the choice of the second atom in the vector is arbitrary;
*       you can "rotate" the last three atoms in the SMILES without altering
*       the implied chirality, e.g. the following three represent the same
*       chirality:
*
*               A[X@](B)(C)D
*               A[X@](C)(D)B
*               A[X@](D)(B)C
*       
*       These two sets of equalities (choice of first atom, choice of second
*       atom) mean there are transformations of the vector that don't alter
*       its meaning.  Using the first atom, we see that the following transformations
*       are available:
*       
*               0 1 2 3                 original order
*               1 0 3 2                 B A D C
*               2 0 1 3                 C A B D
*               3 0 2 1                 D A C B
*
*       Since the choice of the second atom is also arbitrary, by "rotatating" the
*       last three atoms, the following transformations are available:
*
*               0 1 2 3                 A B C D         Original order
*               0 2 3 1                 A C D B
*               0 3 1 2                 A D B C
*
*       This function uses these transformations to determine whether two
*       vectors represent the same or opposite chirality for a particular atom.
*       Given two vectors, v1 and v2:
*
*               Transform v2 so that v2[0] == v1[0]
*               Transform v2 so that v2[1] == v1[1]
*
*       After these two transformations, the third and fourth atoms of v1
*       and v2 will either be the same or opposite, indication that the
*       chirality represented by the vectors is the same or opposite.
***************************************************************************/

bool OBMol2Cansmi::SameChirality(vector<OBAtom*> &v1, vector<OBAtom*> &v2)
{
  vector<OBAtom*> vtmp;

  // First transform v2 so that the first atom matches v1
  if (v2[1] == v1[0]) {
    vtmp[0] = v2[1];
    vtmp[1] = v2[0];
    vtmp[2] = v2[3];
    vtmp[3] = v2[2];
    v2 = vtmp;
  }
  else if (v2[2] == v1[0]) {
    vtmp[0] = v2[2];
    vtmp[1] = v2[0];
    vtmp[2] = v2[1];
    vtmp[3] = v2[3];
    v2 = vtmp;
  }
  else if (v2[3] == v1[0]) {
    vtmp[0] = v2[3];
    vtmp[1] = v2[0];
    vtmp[2] = v2[2];
    vtmp[3] = v2[1];
    v2 = vtmp;
  }
  // else -- the first atoms already match.

  // Now rotate the last three atoms of v2 so that the
  // second atom matchs v1's second atom

  if (v1[1] == v2[2]) {
    v2[2] = v2[3];
    v2[3] = v2[1];
    v2[1] = v1[1];      // use v1[1] rather than tmp var since it's got what we need
  }
  else if (v1[1] == v2[3]) {
    v2[3] = v2[2];
    v2[2] = v2[1];
    v2[1] = v1[1];      // ditto re tmp usage
  }

  // Now, the third and fourth atoms of v1/v2 are the same or opposite, indicating
  // the same or opposite chirality.
  return (v1[3] == v2[3]);
}

/***************************************************************************
* FUNCTION: AtomIsChiral
*
* DESCRIPTION:
*       Returns TRUE if the atom is genuinely chiral, that is, it meets
*       the criteria from OBAtom::IsChiral, and additionally it actually
*       has a connected hash or wedge bond.
* 
*       We arbitrarily reject chiral nitrogen because for our purposes there's
*       no need to consider it.
*
*       NOTE: This is a simplistic test.  When the full SMILES canonicalization
*       includes chiral markings, this should check the symmetry classes
*       of the neighbors, not the hash/wedge bonds.
***************************************************************************/

bool OBMol2Cansmi::AtomIsChiral(OBAtom *atom)
{
  if (!atom->IsChiral())
    return false;
  if (atom->IsNitrogen())
    return false;
  vector<int> symclass;
  FOR_BONDS_OF_ATOM(bond, atom) {
    if (bond->IsHash() || bond->IsWedge())
      return true;
  }
  return false;
}

/***************************************************************************
* FUNCTION: GetChiralStereo
*
* DESCRIPTION:
*       If the atom is chiral, fills in the string with either '@' or '@@',
*       and returns true, otherwise returns false.
***************************************************************************/

bool OBMol2Cansmi::GetChiralStereo(OBCanSmiNode *node,
                                   vector<OBAtom*> &chiral_neighbors,
                                   vector<unsigned int> &symmetry_classes,
                                   char *stereo)
{
  bool is2D=false;
  double torsion;
  OBAtom *atom = node->GetAtom();
  OBMol *mol = (OBMol*) atom->GetParent();

  // If the molecule has no coordinates but DOES have chirality specified, it
  // must have come from a SMILES.  In this case, the atoms' GetIdx() values 
  // will be in the same order they appeared in the original SMILES, so we
  // can deduce the meaning of @ or @@ via "IsClockwise()" or "IsAnticlockwise()".
  // For example, if X is the center atom and A,B,C,D are the connected atoms,
  // appearing sequentially in the input SMILES, then A[X@](B)(C)D is
  //              
  //             B
  //            /
  //      A -- X
  //           |\ 
  //           C D (up wedge bond on D)
  //
  // and "@@" would be the opposite (with C and D switched).
  //

  if (!mol->HasNonZeroCoords()) {               // no coordinates?

    // NOTE: THIS SECTION IS WRONG.  IT'S JUST A COPY OF THE ORIGINAL OPENBABEL
    // CODE, AND DOESN'T ACCOUNT FOR THE FACT THAT THE CANONICAL SMILES IS REORDERED.
    // NEEDS TO BE REWRITTEN, BUT IN COORDINATION WITH A REWRITE OF CHIRALITY IN
    // THE smilesformat.cpp FILE.  -- CJ

    if (!atom->HasChiralitySpecified())         //   and no chirality on this atom?
      return(false);                            //   not a chiral atom -- all done.

    // Ok, it's a chiral atom, so we need to get the A, B, C, D atoms, in the order in
    // which they appeared in the original SMILES.  (NYI!!)
    if (atom->IsClockwise())
      strcpy(stereo,"@@");
    else if (atom->IsAntiClockwise())
      strcpy(stereo,"@");
    else
      return(false);
    return(true);
  }

  // If no chiral neighbors were passed in, we're done
  if (chiral_neighbors.size() < 4)
    return false;

  // If any of the neighbors have the same symmetry class, we're done.
  for (int i = 0; i < chiral_neighbors.size(); i++) {
    int idx = chiral_neighbors[i]->GetIdx();
    int symclass = symmetry_classes[idx-1];
    for (int j = i+1; j < chiral_neighbors.size(); j++) {
      int idx = chiral_neighbors[j]->GetIdx();
      if (symclass == symmetry_classes[idx-1])
        return false;
    }
  }

  // We have 3D coordinates for the four neighbor atoms of the chiral
  // center.  Use the "torsion angle" to deduce chirality.  If you're not
  // familiar with this, it helps to draw it on paper.  Imagine you have
  // A[X](B)(C)D.  Draw three vectors: A--X, X--B, B--C.  These three
  // vectors form a "torsion angle": If you imagine looking at X--B "end
  // on", the vectors A--X and B--C would form an angle.  If you flip the
  // chirality of X, that angle stays the same magnitude, but its sign
  // changes; thus, we can tell the chirality by whether the torsion angle
  // is positive or negative.  (Note: GetVector() is a bad name; it should
  // be called GetXYZ()).

  torsion = CalcTorsionAngle(chiral_neighbors[0]->GetVector(),
                             chiral_neighbors[1]->GetVector(),
                             chiral_neighbors[2]->GetVector(),
                             chiral_neighbors[3]->GetVector());

  strcpy(stereo,(torsion < 0.0) ? "@" : "@@");

  return(true);
}


/***************************************************************************
* FUNCTION: BuildCanonTree
*
* DESCRIPTION:
*       Builds the SMILES tree, in canonical order, for the specified
*       molecular fragment.
***************************************************************************/

bool OBMol2Cansmi::BuildCanonTree(OBMol &mol,
                                  OBBitVec &frag_atoms,
                                  vector<unsigned int> &canonical_order,
                                  OBCanSmiNode *node)
{
  vector<OBEdgeBase*>::iterator i;
  OBAtom *nbr, *atom;
  vector<OBAtom *> sort_nbrs;
  vector<OBAtom *>::iterator ai;
  OBBond *bond;
  OBCanSmiNode *next;
  int idx, canorder;

  atom = node->GetAtom();

#if DEBUG
  cout << "BuildCanonTree: " << etab.GetSymbol(atom->GetAtomicNum()) << ", " << atom->GetIdx() << ", canorder " << canonical_order[atom->GetIdx()-1] << "\n";
#endif

  // Create a vector of neighbors sorted by canonical order, but favor
  // double and triple bonds over single and aromatic.  This causes
  // ring-closure digits to avoid double and triple bonds.
  //
  // Since there are typically just one to three neighbors, we just do a
  // ordered insertion rather than sorting.

  for (nbr = atom->BeginNbrAtom(i); nbr; nbr = atom->NextNbrAtom(i)) {

    idx = nbr->GetIdx();
    if (   (nbr->IsHydrogen() && IsSuppressedHydrogen(nbr))
        || _uatoms[idx]
        || !frag_atoms.BitIsOn(idx))
      continue;

    OBBond *nbr_bond = atom->GetBond(nbr);
    int new_needs_bsymbol = nbr_bond->IsDouble() || nbr_bond->IsTriple();

    for (ai = sort_nbrs.begin(); ai != sort_nbrs.end(); ai++) {
      bond = atom->GetBond(*ai);
      int sorted_needs_bsymbol = bond->IsDouble() || bond->IsTriple();
      if (new_needs_bsymbol && !sorted_needs_bsymbol) {
        sort_nbrs.insert(ai, nbr);
        break;
      }
      if (   new_needs_bsymbol == sorted_needs_bsymbol
          && canonical_order[idx-1] < canonical_order[(*ai)->GetIdx()-1]) {
        sort_nbrs.insert(ai, nbr);
        break;
      }
    }
    if (ai == sort_nbrs.end())
      sort_nbrs.push_back(nbr);
  }

  _uatoms.SetBitOn(atom->GetIdx());     //mark the atom as visited
  _atmorder.push_back(atom->GetIdx());  //store the atom ordering

  // Build the next layer of nodes, in canonical order
  for (ai = sort_nbrs.begin(); ai != sort_nbrs.end(); ai++) {
    nbr = *ai;
    idx = nbr->GetIdx();
    if (_uatoms[idx])   // depth-first search may have used this atom since
      continue;         // we sorted the bonds above
    bond = atom->GetBond(nbr);
    _ubonds.SetBitOn(bond->GetIdx());
    next = new OBCanSmiNode(nbr);
    next->SetParent(atom);
    node->AddChildNode(next, bond);
    BuildCanonTree(mol, frag_atoms, canonical_order, next);
  }

  return(true);
}



/***************************************************************************
* FUNCTION: GetCanonClosureDigits
*
* DESCRIPTION:
*       Given an atom, returns the ring-closure digits for that atom, in
*       the form of a vector of digit/OBBond* pair.  Some of the digits may
*       be for newly-opened rings (the matching digit occurs later in the
*       SMILES string), and some may be for closing rings (the matching
*       digit occured earlier in the string).
*
*       Canonicalization requires that atoms with more than one digit
*       have the digits assigned in a canonical fashion.  For example,
*       the SMILES  "CC12(NCCC2)CCC1" and "CC12(NCCC1)CCC2" are the
*       same molecule; we need to assign the digits of the first "C12"
*       such that it always comes out one way or the other.
*
*       This needs to find closing bonds (ring bonds already assigned a 
*       digit) and opening bonds (ring bonds not encountered yet).
*
*    Closing Bonds:
*       This is easy: open bonds are already stored in the _vopen vector, 
*       in canonical order.  Just find open bonds to this atom and copy
*       them from _vopen to our return vector.
*
*    Opening Bonds:
*       This function looks through the bonds for this atoms and finds
*       any that aren't on the _ubonds "used" list, (and also are non-H
*       and are in this fragment).  Any such bonds must be ring-closure
*       bonds.  If there is more than one, they are ordered by the
*       canonical order of the bonds' neighbor atoms; that is, the bond 
*       to the lowest canonical-ordered neighbor is assigned the first
*       available number, and upwards in neighbor-atom canonical order.
***************************************************************************/

vector<OBBondClosureInfo>
OBMol2Cansmi::GetCanonClosureDigits(OBAtom *atom,
                                    OBBitVec &frag_atoms,
                                    vector<unsigned int> &canonical_order)
{
  vector<OBBondClosureInfo> vp_closures;
  vector<OBBond*> vbonds;
  vector<OBBond*>::iterator bi;
  vector<OBEdgeBase*>::iterator i;
  OBBond *bond1, *bond2;
  OBAtom *nbr1, *nbr2;
  int nbr1_canorder, nbr2_canorder;

  vp_closures.clear();
  vbonds.clear();

  // Find new ring-closure bonds for this atom
  for (bond1 = atom->BeginBond(i); bond1; bond1 = atom->NextBond(i)) {

    // Is this a ring-closure neighbor?
    if (_ubonds.BitIsOn(bond1->GetIdx()))
      continue;
    nbr1 = bond1->GetNbrAtom(atom);
    nbr1_canorder = canonical_order[nbr1->GetIdx()-1];
    if (   (nbr1->IsHydrogen() && IsSuppressedHydrogen(nbr1))
        || !frag_atoms.BitIsOn(nbr1->GetIdx()))
      continue;

    // Insert into the bond-vector in canonical order (by neighbor atom order)
    for (bi = vbonds.begin(); bi != vbonds.end(); bi++) {
      bond2 = *bi;
      nbr2 = bond2->GetNbrAtom(atom);
      nbr2_canorder = canonical_order[nbr2->GetIdx()-1];
      if (nbr1_canorder < nbr2_canorder) {
        vbonds.insert(bi, bond1);
        break;
      }
    }
    if (bi == vbonds.end())     // highest one (or first one) - append to end
      vbonds.push_back(bond1);
  }

  // If we found new open bonds, assign a bond-closure digits to each one,
  // add it to _vopen, and add it to the return vector.
  for (bi = vbonds.begin(); bi != vbonds.end(); bi++) {
    bond1 = *bi;
    _ubonds.SetBitOn(bond1->GetIdx());
    int digit = GetUnusedIndex();
    int bo = (bond1->IsAromatic())? 1 : bond1->GetBO();
    _vopen.push_back(OBBondClosureInfo(bond1->GetNbrAtom(atom), atom, bond1, digit, true));
    vp_closures.push_back(OBBondClosureInfo(bond1->GetNbrAtom(atom), atom, bond1, digit, true));
  }


  // Now look through the list of open closure-bonds and find any to this
  // atom (but watch out for the ones we just added).  For each one found,
  // add it to the return vector, and erase it from _vopen.

  if (!_vopen.empty()) {
    vector<OBBondClosureInfo>::iterator j;
    for (j = _vopen.begin(); j != _vopen.end(); ) {
      if (j->toatom == atom) {
        OBBondClosureInfo bci = *j;
        _vopen.erase(j);                // take bond off "open" list
        bci.is_open = false;            // mark it "closed"
        vp_closures.push_back(bci);     // and add it to this atom's list
        j = _vopen.begin();             // reset iterator
      }
      else
        j++;
    }
  }

  return(vp_closures);
}


/***************************************************************************
* FUNCTION: IsSuppressedHydrogen
*
* DESCRIPTION:
*       For a hydrogen atom, returns TRUE if the atom is not [2H] or [3H], only
*       has one bond, and is not bonded to another hydrogen. 
*
*       NOTE: Return value is nonsensical if you pass it a non-hydrogen
*       atom.  Presumably, you're calling this because you've found a 
*       hydrogen and want to know if it goes in the SMILES.
***************************************************************************/

bool OBMol2Cansmi::IsSuppressedHydrogen(OBAtom *atom)
{
  if (atom->GetIsotope() != 0)          // Deuterium or Tritium
    return false;
  if (atom->GetValence() != 1)          // not exactly one bond
    return false;

  FOR_NBORS_OF_ATOM(nbr, atom) {
    if (nbr->GetAtomicNum() == 1)       // neighbor is hydrogen
      return false;
  }

  return true;
}

/***************************************************************************
* FUNCTION: GetSmilesValence
*
* DESCRIPTION:
*       This is like GetHvyValence(), but it returns the "valence" of an
*       atom as it appears in the SMILES string.  In particular, hydrogens
*       count if they will appear explicitely -- see IsSuppressedHydrogen()
*       above.
***************************************************************************/

int OBMol2Cansmi::GetSmilesValence(OBAtom *atom)
{
  int count = 0;

  if (atom->IsHydrogen())
    return atom->GetValence();

  FOR_NBORS_OF_ATOM(nbr, atom) {
    if (  !nbr->IsHydrogen()
        || nbr->GetIsotope() != 0
        || nbr->GetValence() != 1)
      count++;
  }
  return(count);
}


/***************************************************************************
* FUNCTION: ToCansmilesString
*
* DESCRIPTION:
*       Recursively writes the canonical SMILES string to a buffer.  Writes
*       this node, then selects each of the child nodes (in canonical
*       order) and writes them.
*
*       Chirality is the tricky bit here.  Before we can write out a chiral
*       atom, we have to "look ahead" to determine the order in which the
*       neighbor atoms are/will be written.
*
*       The SMILES language defines the order-of-appearance of a ring-closure
*       bond as the position of the digit, in the SMILES, not the actual atom.
*       For example, the fragments N[C@H](C)Br, and N[C@H]1(Br)CCCC1 have
*       the same chiral center, because the "1" in the second one is a "stand
*       in" for the "C" in the first, even though the actual carbon atom appears
*       after the Bromine atom in the second string.
***************************************************************************/

void OBMol2Cansmi::ToCansmilesString(OBCanSmiNode *node,
                                     char *buffer,
                                     OBBitVec &frag_atoms,
                                     vector<unsigned int> &symmetry_classes,
                                     vector<unsigned int> &canonical_order)
{
  OBAtom *atom = node->GetAtom();
  vector<OBAtom *> chiral_neighbors;

  // Get the ring-closure digits in canonical order.  We'll use these in
  // two places: First, for figuring out chirality, then later for writing
  // the actual ring-closure digits to the string.
  vector<OBBondClosureInfo> vclose_bonds = GetCanonClosureDigits(atom, frag_atoms, canonical_order);

  // First thing: Figure out chirality.  We start by creating a vector of the neighbors
  // in the order in which they'll appear in the canonical SMILES string.  This is more
  // complex than you'd guess because of implicit/explicit H and ring-closure digits.

  bool is_chiral = AtomIsChiral(atom);
  if (is_chiral) {

    // If there's a parent node, it's the first atom in the ordered neighbor-vector
    // used for chirality.
    if (node->GetParent()) {
      chiral_neighbors.push_back(node->GetParent());
    }

    // Next for chirality order will be hydrogen -- since it occurs
    // inside the atom's [] brackets, it's always before other neighbors.
    //
    // Note that we check the regular neighbor list, NOT the canonical
    // SMILES tree, because hydrogens normally aren't part of the canonical
    // SMILES, but we still need them to figure out chirality.
    //
    // There are two cases: it's explicit in the OBMol object but should be
    // written inside the brackets, i.e. "[C@H]", or it is explicit and
    // must be outside the brackets, such as for deuterium.  (A hydrogen
    // that will appear explicitely in the SMILES as a separate atom is
    // treated like any other atom when calculating the chirality.)

    FOR_NBORS_OF_ATOM(i_nbr, atom) {
      OBAtom *nbr = &(*i_nbr);
      if (nbr->IsHydrogen() && IsSuppressedHydrogen(nbr) ) {
        chiral_neighbors.push_back(nbr);
        break;        // quit loop: only be one H if atom is chiral
      }
    }

    // Ok, done with H.  Next in the SMILES will be the ring-closure characters.
    // So we need to find the corresponding atoms and add them to the list.
    // (We got the canonical ring-closure list earlier.)
    if (!vclose_bonds.empty()) {
      vector<OBBondClosureInfo>::iterator i;
      for (i = vclose_bonds.begin();i != vclose_bonds.end();i++) {
        OBBond *bond = i->bond;
        OBAtom *nbr = bond->GetNbrAtom(atom);
        chiral_neighbors.push_back(nbr);
      }
    }
    
    // Finally, add the "regular" neighbors, the "child" nodes in the
    // canonical-SMILES tree, to the chiral-neighbors list.
    for (int i = 0; i < node->Size(); i++) {
      OBAtom *nbr = node->GetChildAtom(i);
      chiral_neighbors.push_back(nbr);
    }
  }

  // Write the current atom to the string
  GetSmilesElement(node, chiral_neighbors, symmetry_classes, buffer+strlen(buffer));

  // Write ring-closure digits
  if (!vclose_bonds.empty()) {
    vector<OBBondClosureInfo>::iterator bci;
    for (bci = vclose_bonds.begin();bci != vclose_bonds.end();bci++) {
      if (!bci->is_open) {
        if (bci->bond->IsUp())                                    strcat(buffer,"\\");
        if (bci->bond->IsDown())                                  strcat(buffer,"/");
        if (bci->bond->GetBO() == 2 && !bci->bond->IsAromatic())  strcat(buffer,"=");
        if (bci->bond->GetBO() == 3)                              strcat(buffer,"#");
      }
      if (bci->ringdigit > 9) strcat(buffer,"%");
      sprintf(buffer+strlen(buffer), "%d", bci->ringdigit); 
    }
  }

  // Write child bonds, then recursively follow paths to child nodes
  // to print the SMILES for each child branch.
  //
  // Note: Cis/trans bonds are tricky, for example, C/C=C/C is trans,
  // but C(/C)=C/C is cis.  If a '/' or '\' bond symbol immediately follows
  // a left parenthesis '(', then it needs to be reversed.  This is a bad
  // hack, and should be replaced with an internal mechanism that represents
  // true cis/trans as an ordered set of four atoms around a double bond.

  OBBond *bond;
  for (int i = 0;i < node->Size();i++) {
    bond = node->GetChildBond(i);
    int up   = bond->IsUp();
    int down = bond->IsDown();
    int open_parens = 0;                // see note above
    if (i+1 < node->Size()) {
      strcat(buffer,"(");
      open_parens = 1;
    }
    if (up && !open_parens || down && open_parens) strcat(buffer,"\\");
    if (down && !open_parens || up && open_parens) strcat(buffer,"/");
    if (bond->GetBO() == 2 && !bond->IsAromatic()) strcat(buffer,"=");
    if (bond->GetBO() == 3)                        strcat(buffer,"#");
    
    ToCansmilesString(node->GetChildNode(i),buffer, frag_atoms, symmetry_classes, canonical_order);
    if (i+1 < node->Size()) strcat(buffer,")");
  }
}


/***************************************************************************
* FUNCTION: CreateFragCansmiString
*
* DESCRIPTION:
*       Selects the "root" atom, which will be first in the SMILES, then
*       builds a tree in canonical order, and finally generates the SMILES.
*       If there are then atoms that haven't been visited (i.e. a molecule
*       with disconnected parts), selects a new root from the remaining
*       atoms and repeats the process.
***************************************************************************/

void OBMol2Cansmi::CreateFragCansmiString(OBMol &mol, OBBitVec &frag_atoms, char *buffer)
{
  OBAtom *atom;
  OBCanSmiNode *root;
  buffer[0] = '\0';
  vector<OBNodeBase*>::iterator ai;
  vector<unsigned int> symmetry_classes, canonical_order;

  // First, create a canonical ordering vector for the atoms.  Canonical
  // labels are zero indexed, corresponding to "atom->GetIdx()-1".
  CanonicalLabels(&mol, frag_atoms, symmetry_classes, canonical_order);

  // OUTER LOOP: Handles dot-disconnected structures.  Finds the 
  // lowest unmarked canorder atom, and starts there to generate a SMILES.
  // Repeats until no atoms remain unmarked.

  while (1) {

    // It happens that the lowest canonically-numbered atom is usually 
    // a good place to start the canonical SMILES.
    OBAtom *root_atom;
    int lowest_canorder = 999999;
    root_atom = NULL;
    for (atom = mol.BeginAtom(ai); atom; atom = mol.NextAtom(ai)) {
      int idx = atom->GetIdx();
      if (!atom->IsHydrogen()           // don't start with a hydrogen
          && !_uatoms[idx]              // skip atoms already used (for fragments)
          && frag_atoms.BitIsOn(idx)            // skip atoms not in this fragment
          //&& !atom->IsChiral()                // don't use chiral atoms as root node
          && canonical_order[idx-1] < lowest_canorder) {
        root_atom = atom;
        lowest_canorder = canonical_order[idx-1];
      }
    }
    if (lowest_canorder == 999999)
      break;

    // Clear out closures in case structure is dot disconnected
    _atmorder.clear();
    _vopen.clear();

    // Dot disconnected structure?
    if (strlen(buffer) > 0) strcat(buffer,"."); 
    root = new OBCanSmiNode (root_atom);

    BuildCanonTree(mol, frag_atoms, canonical_order, root);
    ToCansmilesString(root, buffer, frag_atoms, symmetry_classes, canonical_order);
    delete root;
  }
}

/***************************************************************************
* FUNCTION: OBMol2Cansmi::AddHydrogenToChiralCenters
*
* DESCRIPTION:
*       Adds an explicit hydrogen to any chiral center that only has three
*       atoms.  This makes analysis much easier since the algorithms can
*       assume that all tetrahedral carbons have four neighbors.
***************************************************************************/

void OBMol2Cansmi::AddHydrogenToChiralCenters(OBMol &mol, OBBitVec &frag_atoms)
{
  bool is_modified = false;

  FOR_ATOMS_OF_MOL(atom, mol)
    {

      if (!frag_atoms[atom->GetIdx()] || !AtomIsChiral(&*atom))
        continue;
      
      if (GetSmilesValence(&*atom) == 3 && atom->GetValence() == 3) {       // implicit H?
        
        // Get the (x,y,z) coordinates where best to put the H
        vector3 v;
        atom->GetNewBondVector(v, 1.0);   // Returns (x,y,z) of the "empty" area, for a new bond
        
        // If we haven't put the molecule into "modify" mode yet, do so now
        if (!is_modified) {
          is_modified = true;
          mol.BeginModify();
        }
        
#if DEBUG
        cout << "AddHydrogenToChiralCenters: Adding H to atom " << atom->GetIdx() << "\n";
#endif
        
        // Add the H atom
        OBAtom *h = mol.NewAtom();
        h->SetAtomicNum(1);
        h->SetType("H");
        mol.AddBond(atom->GetIdx(), h->GetIdx(), 1, 0, -1);
        
        // Set its (x,y,z) coordinates
        h->SetVector(v);
        
        frag_atoms.SetBitOn(h->GetIdx());
      }
    }
  if (is_modified)
    mol.EndModify();
}

/*----------------------------------------------------------------------
 * END OF CLASS: OBMol2Cansmi
 ----------------------------------------------------------------------*/



/***************************************************************************
* FUNCTION: CreateCansmiString
*
* DESCRIPTION:
*       Writes the canonical SMILES for a molecule or molecular fragment
*       to the given buffer.
*
*       frag_atoms represents atoms in a fragment of the molecule; the
*       SMILES will contain those atoms only.
*
*       (Note: This is an ordinary public C++ function, not a member
*       of any class.)
*
***************************************************************************/

void CreateCansmiString(OBMol &mol, char *buffer, OBBitVec &frag_atoms, bool iso)
{
  char tmp[BUFF_SIZE];
  int chg;
  char *p, *pp;

  // This is a hack to prevent recursion problems.
  //  we still need to fix the underlying problem -GRH
  if (mol.NumAtoms() > 1000) {
#ifdef HAVE_SSTREAM
    stringstream errorMsg;
#else
    strstream errorMsg;
#endif
    errorMsg <<
      "SMILES Conversion failed: Molecule is too large to convert."
      "Open Babel is currently limited to 1000 atoms." << endl;
    errorMsg << "  Molecule size: " << mol.NumAtoms() << " atoms " << endl;
    obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
    return;
  }

  // If we're doing isomeric (stereo), make a copy
  OBMol *pmol;
  if (iso) 
    pmol = new OBMol(mol);
  else
    pmol = &mol;

  OBMol2Cansmi m2s;
  m2s.Init();

  // Figure out Cis/Trans 
  if (mol.Has2D())
    m2s.AssignCisTrans(pmol);

  // If the molecule has 2D coordinate AND has hash/wedge bonds,
  // create pseudo-Z coordinates by "pushing" the up/down bonds
  // to +/-1 in the Z direction.  This will be used in the next
  // section when we deduce chirality from the coordinates.

  if (iso) {
    if (!pmol->Has3D()) {

      FOR_ATOMS_OF_MOL(iatom, *pmol) {
        OBAtom *atom = &(*iatom);

        if (!atom->IsChiral()) continue;
        if (m2s.GetSmilesValence(atom) < 3) continue;

        vector3 v;
        OBAtom *nbr;
        OBBond *bond;

        FOR_BONDS_OF_ATOM(bond, atom) {

          // The bond's "start atom" is the pointy end of the hash or wedge
          // bond, so we need to know whether the pointy end of the bond is
          // toward the center atom (normal case) or toward the neighbor atom
          // (poor drawing style, but it happens).  The amount to push up/down
          // is "z", and is normally 1.0, but is set to 0.5 for non-terminal
          // atoms.  This keeps adjacent chiral centers from screwing each other up.

          nbr = bond->GetNbrAtom(atom);
          double z = (nbr->GetHvyValence() > 1) ? 0.5 : 1.0;
          v = nbr->GetVector();
          if (bond->GetBeginAtom() == atom) {       // The pointy end is at the central atom
            if (bond->IsWedge())
              v.SetZ(z);
            else if (bond->IsHash())
              v.SetZ(-z);
          }
          else {                                    // The pointy end is at the neighbor atom
            if (bond->IsWedge())
              v.SetZ(-z);
            else if (bond->IsHash())
              v.SetZ(z);
          }
          nbr->SetVector(v);
        }
      }
    }

    m2s.AddHydrogenToChiralCenters(*pmol, frag_atoms);
  }

  else {
    // Not isomeric - be sure there are no Z coordinates, clear
    // all stereo-center and cis/trans information.
    OBBond *bond;
    OBAtom *atom;
    vector<OBEdgeBase*>::iterator bi;
    vector<OBNodeBase*>::iterator ai;
    for (bond = pmol->BeginBond(bi); bond; bond = pmol->NextBond(bi)) {
      bond->UnsetUp();
      bond->UnsetDown();
      bond->UnsetHash();
      bond->UnsetWedge();
    }
    for (atom = pmol->BeginAtom(ai); atom; atom = pmol->NextAtom(ai)) {
      atom->UnsetStereo();
      vector3 v = atom->GetVector();
      if (v[2] != 0.0) {
        v.SetZ(0.0);
        atom->SetVector(v);
      }
    }
  }

  // If the fragment includes ordinary hydrogens, get rid of them.
  // They won't appear in the SMILES (unless they're attached to a chiral
  // center) anyway.
  FOR_ATOMS_OF_MOL(iatom, *pmol) {
    OBAtom *atom = &(*iatom);
    if (frag_atoms.BitIsOn(atom->GetIdx()) && atom->IsHydrogen() && (!iso || m2s.IsSuppressedHydrogen(atom))) {
      frag_atoms.SetBitOff(atom->GetIdx());
    }
  }

  m2s.CreateFragCansmiString(*pmol, frag_atoms, buffer);
  if (iso) {
    pmol->Clear();
    delete pmol;
  }
}


/*----------------------------------------------------------------------
* CLASS: CANSMIFormat 
----------------------------------------------------------------------*/

#if OPENBABEL_FORMAT


class CANSMIFormat : public OBMoleculeFormat
{
public:
  //Register this format type ID
  CANSMIFormat()
  {
    OBConversion::RegisterFormat("can", this, "chemical/x-daylight-cansmiles");
    OBConversion::RegisterOptionParam("n", this);
    OBConversion::RegisterOptionParam("t", this);
  }

  virtual const char* GetMIMEType() { return "chemical/x-daylight-smiles"; };

  // CANSMI format is write only
  virtual unsigned int Flags()  {return NOTREADABLE;};

  ////////////////////////////////////////////////////
  /// The "API" interface functions
  virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

  ///////////////////////////////////////////////////////

  virtual const char* Description() {
    return
      "Canonical SMILES format.\n"
      " A linear text format which can describe the connectivity\n"
      " and chirality of a molecule, and has a single 'canonical'\n"
      " form for any particular molecule.\n"
      " Write Options e.g. -xt\n"
      "   -i  Includes isotopic and chiral markings\n"
      "   -n  No molecule name\n"
      "   -t  Molecule name only\n"
      "   -r  Radicals lower case eg ethyl is Cc\n\n";
  };

  virtual const char* SpecificationURL()
  {return "http://www.daylight.com/smiles/f_smiles.html";};

  virtual int SkipObjects(int n, OBConversion* pConv)
  {
    if (n==0) return 1;         // already points after current line
    string temp;
    istream& ifs = *pConv->GetInStream();
    int i;
    for(i=0;i<n && ifs.good();i++)
      getline(ifs, temp);
    return ifs.good() ? 1 : -1; 
  };  
};

// Make an instance of the format class
CANSMIFormat theCANSMIFormat;

//////////////////////////////////////////////////

bool CANSMIFormat::WriteMolecule(OBBase* pOb,OBConversion* pConv)
{
  OBMol* pmol = dynamic_cast<OBMol*>(pOb);

  // Define some references so we can use the old parameter names
  ostream &ofs = *pConv->GetOutStream();
  OBMol &mol = *pmol;

  // Title only option?
  if(pConv->IsOption("t")) {
      ofs << mol.GetTitle() <<endl;
      return true;
    }

  char buffer[BUFF_SIZE];
  *buffer = '\0'; // clear the buffer

  // This is a hack to prevent recursion problems.
  //  we still need to fix the underlying problem (mainly chiral centers) -GRH
  if (mol.NumAtoms() > 1000) {
#ifdef HAVE_SSTREAM
    stringstream errorMsg;
#else
    strstream errorMsg;
#endif
    errorMsg <<
      "SMILES Conversion failed: Molecule is too large to convert."
      "Open Babel is currently limited to 1000 atoms." << endl;
    errorMsg << "  Molecule size: " << mol.NumAtoms() << " atoms " << endl;
    obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
    return(false);
  }

  OBBitVec allbits(mol.NumAtoms());
  FOR_ATOMS_OF_MOL(a, mol)
    {
      allbits.SetBitOn(a->GetIdx());
    }

  if (mol.NumAtoms() != 0) {
      OBMol2Cansmi m2s;
      m2s.Init(pConv);
      m2s.CorrectAromaticAmineCharge(mol);
      CreateCansmiString(mol, buffer, allbits, true);
    }

  ofs << buffer ;
  if(!pConv->IsOption("n"))
    ofs << '\t' <<  mol.GetTitle();
  ofs << endl;

  return true;
}

#endif  // OPENBABEL_FORMAT



}

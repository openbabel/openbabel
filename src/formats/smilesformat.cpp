/**********************************************************************
Copyright (C) 2005-2007 by Craig A. James, eMolecules Inc.
Some portions Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2008 by Geoffrey R. Hutchison
Some portions Copyright (C) 2004 by Chris Morley
Some portions Copyright (C) 2019 by NextMove Software.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

// This code uses the old OpenEye SMILES parser
// but replaces the SMILES export with Craig James canonical smiles
// (For regular SMILES, the canonical order is not computed and ignored)

#include <openbabel/babelconfig.h>
#include <openbabel/obmolecformat.h>

#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obiter.h>
#include <openbabel/elements.h>
#include <openbabel/generic.h>


#include <openbabel/stereo/tetrahedral.h>
#include <openbabel/stereo/cistrans.h>
#include <openbabel/stereo/squareplanar.h>
#include <openbabel/reactionfacade.h>
#include <openbabel/stereo/stereo.h>
#include <openbabel/obfunctions.h>
#include <openbabel/graphsym.h>
#include <openbabel/kekulize.h>
#include <openbabel/canon.h>

#include "smilesvalence.h"

#include <limits>
#include <iostream>
#include <cassert>
#include <string>

//#define DEBUG 1
#define IMPLICIT_CIS_RING_SIZE 8

using namespace std;

namespace OpenBabel {

  // some constant variables
  const char BondUpChar = '\\';
  const char BondDownChar = '/';

  // This function return true for sulfur and nitrogen
  // (I'm not sure that is the right approach, longterm)
  static bool CanHaveLonePair(unsigned int elem)
  {
    switch (elem) {
    case OBElements::Nitrogen:
    case OBElements::Sulfur:
      return true;
    default:
      return false;
    }
  }

  //Base class for SMIFormat and CANSIFormat with most of the functionality
  class SMIBaseFormat : public OBMoleculeFormat
  {
  public:
    virtual const char* GetMIMEType()
    { return "chemical/x-daylight-smiles"; };

    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

    ///////////////////////////////////////////////////////

    virtual const char* TargetClassDescription(){return OBMol::ClassDescription();};

    virtual const char* SpecificationURL()
    {return "http://www.daylight.com/smiles/";};

    virtual int SkipObjects(int n, OBConversion* pConv)
    {
      if(n==0) return 1; //already points after current line
      istream& ifs = *pConv->GetInStream();
      if (ifs.eof())
        return -1;

      int i=0;
      while(i<n && ifs.good())
        {
          if(ifs.peek()!='#')
            i++;
          ifs.ignore(numeric_limits<streamsize>::max(),'\n');
        }
      return ifs ? 1 : -1;
    }
  private:
    bool GetInchifiedSMILESMolecule(OBMol *mol, bool useFixedHRecMet);
  };

  //**************************************************
  class SMIFormat : public SMIBaseFormat
  {
  public:
    //Register this format type ID
    SMIFormat()
    {
      OBConversion::RegisterFormat("smi",this, "chemical/x-daylight-smiles");
      OBConversion::RegisterFormat("smiles",this, "chemical/x-daylight-smiles");
      OBConversion::RegisterOptionParam("n", this);
      OBConversion::RegisterOptionParam("t", this);
      OBConversion::RegisterOptionParam("r", this);
      OBConversion::RegisterOptionParam("a", this);
      OBConversion::RegisterOptionParam("h", this);
      OBConversion::RegisterOptionParam("x", this);
      OBConversion::RegisterOptionParam("C", this);	// "anti-canonical" form (random order)
    }
    virtual const char* Description()
    {
      return
        "SMILES format\n"
        "A linear text format which can describe the connectivity and chirality of a molecule\n"
        "Open Babel implements the `OpenSMILES specification <http://opensmiles.org>`_.\n\n"

        "It also implements an extension to this specification for radicals.\n\n"

        "Note that the ``l <atomno>`` option, used to specify a \"last\" atom, is\n"
        "intended for the generation of SMILES strings to which additional atoms\n"
        "will be concatenated. If the atom specified has an explicit H within a bracket\n"
        "(e.g. ``[nH]`` or ``[C@@H]``) the output will have the H removed along with any\n"
        "associated stereo symbols.\n\n"

        ".. seealso::\n\n"

        "  The :ref:`Canonical_SMILES_format` produces a canonical representation\n"
        "  of the molecule in SMILES format. This is the same as the ``c`` option\n"
        "  below but may be more convenient to use.\n\n"

        "Write Options e.g. -xt\n"
        "  a  Output atomclass like [C:2], if available\n"
        "  c  Output in canonical form\n"
        "  U  Universal SMILES\n"
        "  I  Inchified SMILES\n"
        "  h  Output explicit hydrogens as such\n"
        "  i  Do not include isotopic or chiral markings\n"
        "  k  Create Kekule SMILES instead of aromatic\n"
        "  n  No molecule name\n"
        "  r  Radicals lower case eg ethyl is Cc\n"
        "  t  Molecule name only\n"
        "  x  append X/Y coordinates in canonical-SMILES order\n"
        "  C  'anti-canonical' random order (mostly for testing)\n"
        "  o  <ordering> Output in user-specified order\n"
        "     Ordering should be specified like 4-2-1-3 for a 4-atom molecule.\n"
        "     This gives canonical labels 1,2,3,4 to atoms 4,2,1,3 respectively,\n"
        "     so that atom 4 will be visited first and the remaining atoms\n"
        "     visited in a depth-first manner following the lowest canonical labels.\n"
        "  O  Store the SMILES atom order as a space-separated string\n"
        "     The string is stored as an OBPairData wth the name\n"
        "     'SMILES Atom Order'.\n"
        "  F  <atom numbers> Generate SMILES for a fragment\n"
        "     The atom numbers should be specified like \"1 2 4 7\".\n"
        "  R  Do not reuse bond closure symbols\n"
        "  f  <atomno> Specify the first atom\n"
        "     This atom will be used to begin the SMILES string.\n"
        "  l  <atomno> Specify the last atom\n"
        "     The output will be rearranged so that any additional\n"
        "     SMILES added to the end will be attached to this atom.\n"
        "  T  <max seconds> Specify the canonicalization timeout\n"
        "     Canonicalization can take a while for symmetric molecules and a\n"
        "     timeout is used. The default is 5 seconds.\n\n"

        "Read Options e.g. -aa\n"
        "  a  Preserve aromaticity present in the SMILES\n"
        "     This option should only be used if reading aromatic SMILES\n"
        "     generated by the same version of Open Babel. Any other\n"
        "     use will lead to undefined behavior. The advantage of this\n"
        "     option is that it avoids aromaticity perception, thus speeding\n"
        "     up reading SMILES.\n"
        "  S  Clean stereochemistry\n"
        "     By default, stereochemistry is accepted as given. If you wish\n"
        "     to clean up stereochemistry (e.g. by removing tetrahedral\n"
        "     stereochemistry where two of the substituents are identical)\n"
        "     then specifying this option will reperceive stereocenters.\n"
        "\n\n"
        ;
    }


  };

  //Make an instance of the format class
  SMIFormat theSMIFormat;

  //**************************************************
  class CANSMIFormat : public SMIBaseFormat
  {
  public:
    //Register this format type ID
    CANSMIFormat()
    {
      OBConversion::RegisterFormat("can", this, "chemical/x-daylight-cansmiles");
    }

    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv)
    {
      //The "c" option sets us to use canonical ordering
      pConv->AddOption("c",OBConversion::OUTOPTIONS);
      return SMIBaseFormat::WriteMolecule(pOb, pConv);
    }

    ///////////////////////////////////////////////////////

    virtual const char* Description() {
      return
        "Canonical SMILES format\n"
        "A canonical form of the SMILES linear text format\n"
        "The SMILES format is a linear text format which can describe the\n"
       	"connectivity "
        "and chirality of a molecule. Canonical SMILES gives a single\n"
       	"'canonical' form for any particular molecule.\n\n"

        ".. seealso::\n\n"

        "  The \"regular\" :ref:`SMILES_format` gives faster\n"
        "  output, since no canonical numbering is performed.\n\n"

        "Write Options e.g. -xt\n"
        "  a  Output atomclass like [C:2], if available\n"
        "  h  Output explicit hydrogens as such\n"
        "  i  Do not include isotopic or chiral markings\n"
        "  n  No molecule name\n"
        "  r  Radicals lower case eg ethyl is Cc\n"
        "  t  Molecule name only\n"
        "  F  <atom numbers> Generate Canonical SMILES for a fragment\n"
        "     The atom numbers should be specified like \"1 2 4 7\".\n"
        "  f  <atomno> Specify the first atom\n"
        "     This atom will be used to begin the SMILES string.\n"
        "  l  <atomno> Specify the last atom\n"
        "     The output will be rearranged so that any additional\n"
        "     SMILES added to the end will be attached to this atom.\n"
        "     See the :ref:`SMILES_format` for more information.\n\n";
    };

  };

  // Make an instance of the format class
  CANSMIFormat theCANSMIFormat;

  //************************************************************

  class OBSmilesParser
  {
    // simple structs to make code more readable

    // see _extbond
    struct ExternalBond
    {
      int digit;
      int prev;
      int order;
      char updown;
    };
    // see _rclose
    struct RingClosureBond
    {
      int digit;
      int prev;
      int order;
      char updown;
      int numConnections;
    };


    char _updown;
    int _order;
    int _prev;
    int _rxnrole;
    const char *_ptr;
    bool _preserve_aromaticity;
    vector<int>             _vprev;
    vector<RingClosureBond> _rclose;
    vector<ExternalBond>    _extbond;
    vector<int>             _path;
    vector<bool>            _avisit;
    vector<bool>            _bvisit;
    vector<int>             _hcount;
    vector<int> PosDouble; //for extension: lc atoms as conjugated double bonds

    struct StereoRingBond
    {
      vector<OBAtom*> atoms;
      vector<char> updown;
    };
    map<OBBond*, StereoRingBond> _stereorbond; // Remember info on the stereo ring closure bonds

    // stereochimistry
    bool chiralWatch; // set when a tetrahedral atom is read
    map<OBAtom*, OBTetrahedralStereo::Config*> _tetrahedralMap; // map of tetrahedral atoms and their data
    map<OBBond*, char> _upDownMap; // store the '/' & '\' as they occured in smiles
    map<unsigned int, char> _chiralLonePair; // for atoms with potential chiral lone pairs, remember when the l.p. was encountered
    bool squarePlanarWatch; // set when a square planar atom is read
    map<OBAtom*, OBSquarePlanarStereo::Config*> _squarePlanarMap;

  public:

    OBSmilesParser(bool preserve_aromaticity=false): _preserve_aromaticity(preserve_aromaticity), _rxnrole(1) { }
    ~OBSmilesParser() { }

    bool SmiToMol(OBMol&,const string&);
    bool ParseSmiles(OBMol&, const string&);
    bool ParseSimple(OBMol&);
    bool ParseComplex(OBMol&);
    bool ParseRingBond(OBMol&);
    bool ParseExternalBond(OBMol&);
    bool CapExternalBonds(OBMol &mol);
    int NumConnections(OBAtom *, bool isImplicitRef=false);
    void CreateCisTrans(OBMol &mol);
    char SetRingClosureStereo(StereoRingBond rcstereo, OBBond* dbl_bond);
    void InsertTetrahedralRef(OBMol &mol, unsigned long id);
    void InsertSquarePlanarRef(OBMol &mol, unsigned long id);

    bool IsUp(OBBond*);
    bool IsDown(OBBond*);
  };

  /////////////////////////////////////////////////////////////////
  /* Lines starting with # are ignored. Whitespace at the start (including
     blank lines) terminate the input unless -e option is used.
     Valid SMILES reactions such as [C]=O.O>[Fe]>O=C=O.[H][H] with non-null
     reactant and product are accepted and the reactant, product and
     possibly the agent molecules are output when using the Convert interface
     (babel commandline). With the OBConversion functions Read, ReadString
     and ReadFile all SMILES reactions give an error when read with this format.
  */
  bool SMIBaseFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = pOb->CastAndClear<OBMol>();

    istream &ifs = *pConv->GetInStream();
    string ln, smiles, title;
    string::size_type pos;

    //Ignore lines that start with #
    while(ifs && ifs.peek()=='#')
      if(!getline(ifs, ln))
        return false;

    //Get title
    if(getline(ifs, ln))
    {
      pos = ln.find_first_of(" \t");
      if(pos!=string::npos)
      {
        smiles = ln.substr(0,pos);
        title = ln.substr(pos+1);
        Trim(title);
        pmol->SetTitle(title.c_str());
      }
      else
        smiles = ln;
    }

    pmol->SetDimension(0);
    OBSmilesParser sp(pConv->IsOption("a", OBConversion::INOPTIONS));
    if (!pConv->IsOption("S", OBConversion::INOPTIONS))
      pmol->SetChiralityPerceived();

    return sp.SmiToMol(*pmol, smiles); //normal return
  }

  //////////////////////////////////////////////

  bool OBSmilesParser::SmiToMol(OBMol &mol,const string &s)
  {
    _vprev.clear();
    _rclose.clear();
    _prev=0;
    chiralWatch=false;
    squarePlanarWatch = false;

    // We allow the empty reaction (">>") but not the empty molecule ("")
    if (!ParseSmiles(mol, s) || (!mol.IsReaction() && mol.NumAtoms() == 0))
      {
        mol.Clear();
        return(false);
      }

    // TODO: Is the following a memory leak? - there are return statements above
    map<OBAtom*, OBTetrahedralStereo::Config*>::iterator i;
    for (i = _tetrahedralMap.begin(); i != _tetrahedralMap.end(); ++i)
      delete i->second;
    _tetrahedralMap.clear();

    map<OBAtom*, OBSquarePlanarStereo::Config*>::iterator j;
    for (j = _squarePlanarMap.begin(); j != _squarePlanarMap.end(); ++j)
      delete j->second;
    _squarePlanarMap.clear();

    mol.SetAutomaticFormalCharge(false);

    return(true);
  }

  bool OBSmilesParser::ParseSmiles(OBMol &mol, const std::string &smiles)
  {
    mol.SetAromaticPerceived(); // Turn off perception until the end of this function
    mol.BeginModify();

    for (_ptr=smiles.c_str();*_ptr;_ptr++)
    {
      switch(*_ptr)
      {
      case '\r':
        if (*(_ptr+1) == '\0') // may have a terminating '\r' due to Windows line-endings
          break;
        return false;
      case '0': case '1': case '2': case '3': case '4':
      case '5': case '6': case '7': case '8': case '9':
      case '%':  //ring open/close
        if (_prev == 0)
          return false;
        if (!ParseRingBond(mol))
          return false;
        break;
      case '&': //external bond
        if (_prev == 0)
          return false;
        if (!ParseExternalBond(mol))
          return false;
        break;
      case '.':
        _prev=0;
        break;
      case '>':
        _prev = 0;
        _rxnrole++;
        if (_rxnrole == 2) {
          mol.SetIsReaction();
          // Handle all the reactant atoms
          // - the remaining atoms will be handled on-the-fly
          FOR_ATOMS_OF_MOL(atom, mol) {
            OBPairInteger *pi = new OBPairInteger();
            pi->SetAttribute("rxnrole");
            pi->SetValue(1);
            atom->SetData(pi);
          }
        }
        else if (_rxnrole == 4) {
          stringstream errorMsg;
          errorMsg << "Too many greater-than signs in SMILES string";
          std::string title = mol.GetTitle();
          if (!title.empty())
            errorMsg << " (title is " << title << ")";
          errorMsg << endl;
          obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
          return false;
        }
        break;
      case '(':
        _vprev.push_back(_prev);
        break;
      case ')':
        if(_vprev.empty()) //CM
          return false;
        _prev = _vprev.back();
        _vprev.pop_back();
        break;
      case '[':
        if (!ParseComplex(mol))
        {
          mol.EndModify();
          mol.Clear();
          return false;
        }
        break;
      case '-':
        if (_prev == 0)
          return false;
        _order = 1;
        break;
      case '=':
        if (_prev == 0)
          return false;
        _order = 2;
        break;
      case '#':
        if (_prev == 0)
          return false;
        _order = 3;
        break;
      case '$':
        if (_prev == 0)
          return false;
        _order = 4;
        break;
      case ':':
        if (_prev == 0)
          return false;
        _order = 0; // no-op
        break;
      case '/':
        if (_prev == 0)
          return false;
        _order = 1;
        _updown = BondDownChar;
        break;
      case '\\':
        if (_prev == 0)
          return false;
        _order = 1;
        _updown = BondUpChar;
        break;
      default:
        if (!ParseSimple(mol))
        {
          mol.EndModify();
          mol.Clear();
          return false;
        }
      } // end switch
    } // end for _ptr

    // place dummy atoms for each unfilled external bond
    if(!_extbond.empty())
      CapExternalBonds(mol);

    // Check to see if we've balanced out all ring closures
    // They are removed from _rclose when matched
    if (!_rclose.empty()) {
      mol.EndModify();
      mol.Clear();

      stringstream errorMsg;
      errorMsg << "Invalid SMILES string: " << _rclose.size() << " unmatched "
               << "ring bonds." << endl;
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
      return false; // invalid SMILES since rings aren't properly closed
    }

    // Check to see if we've the right number of '>' for reactions
    if (_rxnrole > 1 && _rxnrole !=3) {
      mol.EndModify();
      stringstream errorMsg;
      errorMsg << "Invalid reaction SMILES string: only a single '>' sign found (two required to be valid).";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
      return false; // invalid SMILES since rings aren't properly closed
    }
    if (mol.IsReaction()) {
      OBReactionFacade facade(&mol);
      facade.AssignComponentIds();
    }

    // Apply the SMILES valence model
    FOR_ATOMS_OF_MOL(atom, mol) {
      unsigned int idx = atom->GetIdx();
      int hcount = _hcount[idx - 1];
      if (hcount == -1) { // Apply SMILES implicit valence model
        unsigned int bosum = 0;
        FOR_BONDS_OF_ATOM(bond, &(*atom)) {
          bosum += bond->GetBondOrder();
        }
        unsigned int impval = SmilesValence(atom->GetAtomicNum(), bosum);
        unsigned int imph = impval - bosum;
        if (imph > 0 && atom->IsAromatic())
          imph--;
        atom->SetImplicitHCount(imph);
      }
      else // valence is explicit e.g. [CH3]
        atom->SetImplicitHCount(hcount);
    }

    mol.EndModify(false);

    // Unset any aromatic bonds that *are not* in rings where the two aromatic atoms *are* in a ring
    // This is rather subtle, but it's correct and reduces the burden of kekulization
    FOR_BONDS_OF_MOL(bond, mol) {
      if (bond->IsAromatic() && !bond->IsInRing()) {
        if (bond->GetBeginAtom()->IsInRing() && bond->GetEndAtom()->IsInRing())
          bond->SetAromatic(false);
      }
    }

    // TODO: Only Kekulize if the molecule has a lower case atom
    bool ok = OBKekulize(&mol);
    if (!ok) {
      stringstream errorMsg;
      errorMsg << "Failed to kekulize aromatic SMILES";
      std::string title = mol.GetTitle();
      if (!title.empty())
        errorMsg << " (title is " << title << ")";
      errorMsg << endl;
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
      // return false; // Should we return false for a kekulization failure?
    }

    // Add the data stored inside the _tetrahedralMap to the atoms now after end
    // modify so they don't get lost.
    if(!_tetrahedralMap.empty()) {
      OBAtom* atom;
      map<OBAtom*, OBTetrahedralStereo::Config*>::iterator ChiralSearch;
      for(ChiralSearch = _tetrahedralMap.begin(); ChiralSearch != _tetrahedralMap.end(); ++ChiralSearch) {
        atom = ChiralSearch->first;
        OBTetrahedralStereo::Config *ts = ChiralSearch->second;
        if (!ts)
          continue;
        if (ts->refs.size() != 3)
          continue;
        if (ts->refs[2] == OBStereo::NoRef) {
          // This happens where there is chiral lone pair or where there simply aren't enough connections
          // around a chiral atom. We handle the case where there is a S with a chiral lone pair.
          // All other cases are ignored, and raise a warning. (Note that S can be chiral even without
          // a lone pair, think of C[S@](=X)(=Y)Cl.

          // We have remembered where to insert the lone pair in the _chiralLonePair map
          map<unsigned int, char>::iterator m_it = _chiralLonePair.find(atom->GetIdx());
          if (CanHaveLonePair(atom->GetAtomicNum()) && m_it != _chiralLonePair.end()) {
            ts->refs[2] = ts->refs[1]; ts->refs[1] = ts->refs[0];
            if (m_it->second == 0) { // Insert in the 'from' position
              ts->refs[0] = ts->from;
              ts->from = OBStereo::ImplicitRef;
            }
            else // Insert in the refs[0] position
              ts->refs[0] = OBStereo::ImplicitRef;
          }
          else { // Ignored by Open Babel
            stringstream ss;
            ss << "Ignoring stereochemistry. Not enough connections to this atom. " << mol.GetTitle();
            obErrorLog.ThrowError(__FUNCTION__, ss.str(), obWarning);
            continue;
          }
        }

        // cout << "*ts = " << *ts << endl;
        OBTetrahedralStereo *obts = new OBTetrahedralStereo(&mol);
        obts->SetConfig(*ts);
        mol.SetData(obts);
      }
    }

    // Add the data stored inside the _squarePlanarMap to the atoms now after end
    // modify so they don't get lost.
    if(!_squarePlanarMap.empty()) {
      OBAtom* atom;
      map<OBAtom*, OBSquarePlanarStereo::Config*>::iterator ChiralSearch;
      for(ChiralSearch = _squarePlanarMap.begin(); ChiralSearch != _squarePlanarMap.end(); ++ChiralSearch) {
        atom = ChiralSearch->first;
        OBSquarePlanarStereo::Config *sp = ChiralSearch->second;
        if (!sp)
          continue;
        if (sp->refs.size() != 4)
          continue;

        // cout << "*ts = " << *ts << endl;
        OBSquarePlanarStereo *obsp = new OBSquarePlanarStereo(&mol);
        obsp->SetConfig(*sp);
        mol.SetData(obsp);
      }
    }

    if (!_preserve_aromaticity)
      mol.SetAromaticPerceived(false);

    CreateCisTrans(mol);

    return(true);
  }

  bool OBSmilesParser::IsUp(OBBond *bond)
  {
    map<OBBond*, char>::iterator UpDownSearch;
    UpDownSearch = _upDownMap.find(bond);
    if (UpDownSearch != _upDownMap.end())
      if (UpDownSearch->second == BondUpChar)
        return true;
    return false;
  }

  bool OBSmilesParser::IsDown(OBBond *bond)
  {
    map<OBBond*, char>::iterator UpDownSearch;
    UpDownSearch = _upDownMap.find(bond);
    if (UpDownSearch != _upDownMap.end())
      if (UpDownSearch->second == BondDownChar)
        return true;
    return false;
  }

  char OBSmilesParser::SetRingClosureStereo(StereoRingBond rcstereo, OBBond* dbl_bond)
  {
    // Ring Closure bonds appear twice (at opening and closure).
    // If involved in cis/trans stereo, then the stereo may be
    // specified at either end or indeed both. Although Open Babel
    // will only write out SMILES with the stereo at one end (the end
    // on the double bond), it must handle all cases when reading.

    // For example:
    //
    //         C
    //        /|
    //   C = C |
    //  /     \|
    // C       N
    //
    // Can be written as:
    // (a) C/C=C/1\NC1 -- preferred
    // (b) C/C=C1\NC\1 
    // (c) C/C=C/1\NC\1
    //  or indeed by replacing the "\N" with "N".

    // If the stereo chemistry for a ring closure is inconsistently specified,
    // it is ignored. In that case, if a stereo symbol does not exist for its
    // partner bond on the double bond (e.g. (b) below), then the stereo is unspecified.

    // (a) C/C=C/1NC\1 -- specified stereo
    // (b) C/C=C/1NC/1  -- ignore ring closure stereo => treated as C/C=C1NC1  => CC=C1NC1
    // (c) C/C=C/1\NC/1 -- ignore ring closure stereo => treated as C/C=C1\NC1 => C/C=C/1\NC1

    // The ring closure bond is either up or down with respect 
    // to the double bond. Our task here is to figure out which it is,
    // based on the contents of _stereorbond.

    bool found = false; // We have found the answer
    bool updown = true; // The answer

    if (rcstereo.updown[0] == BondUpChar || rcstereo.updown[0] == BondDownChar) { // Is there a stereo symbol at the opening?
      bool on_dbl_bond = (rcstereo.atoms[0] == dbl_bond->GetBeginAtom() || rcstereo.atoms[0] == dbl_bond->GetEndAtom());
      updown = (rcstereo.updown[0]==BondUpChar) ^ on_dbl_bond;
      found = true;
    }
    if (rcstereo.updown[1] == BondUpChar || rcstereo.updown[1] == BondDownChar) { // Is there a stereo symbol at the closing?
      bool on_dbl_bond = (rcstereo.atoms[1] == dbl_bond->GetBeginAtom() || rcstereo.atoms[1] == dbl_bond->GetEndAtom());
      bool new_updown = (rcstereo.updown[1]==BondUpChar) ^ on_dbl_bond;
      if (!found) {
        updown = new_updown;
        found = true;
      }
      else if (new_updown != updown) {
        obErrorLog.ThrowError(__FUNCTION__, "Ignoring the cis/trans stereochemistry specified for the ring closure\n  as it is inconsistent.", obWarning);
        found = false;
      }
    }
    
    if (!found)
      return 0;
    else
      return updown ? 1 : 2;
  }

  void OBSmilesParser::CreateCisTrans(OBMol &mol)
  {
    // Create a vector of CisTransStereo objects for the molecule
    FOR_BONDS_OF_MOL(dbi, mol) {

      OBBond *dbl_bond = &(*dbi);

      // Not a double bond?
      if (dbl_bond->GetBondOrder()!=2 || dbl_bond->IsAromatic())
        continue;

      // Find the single bonds around the atoms connected by the double bond.

      OBAtom *a1 = dbl_bond->GetBeginAtom();
      OBAtom *a2 = dbl_bond->GetEndAtom();

      // Check that both atoms on the double bond have at least one
      // other neighbor, but not more than two other neighbors;
      // Note: In theory, we could relax the second requirement but we would
      //       need to change the data structure we use to store cis/trans
      //       stereo to only store 2 refs instead of 4
      int v1 = a1->GetExplicitDegree();
      int v2 = a2->GetExplicitDegree();
      if (v1 < 2 || v1 > 3 || v2 < 2 || v2 > 3) {
        continue;
      }

      vector<OBAtom*> dbl_bond_atoms;
      dbl_bond_atoms.push_back(a1);
      dbl_bond_atoms.push_back(a2);

      vector<bool> bond_stereo(2, true); // Store the stereo of the chosen bonds at each end of the dbl bond
      vector<OBBond*> stereo_bond(2, (OBBond*) NULL); // These are the chosen stereo bonds
      vector<OBBond*> other_bond(2, (OBBond*) NULL);  // These are the 'other' bonds at each end
      
      for (int i = 0; i < 2; ++i) { // Loop over each end of the double bond in turn

        FOR_BONDS_OF_ATOM(bi, dbl_bond_atoms[i]) {
          OBBond *b = &(*bi);
          if (b == dbl_bond) continue;
          if (!(IsUp(b) || IsDown(b))) {
            other_bond[i] = b; // Use this for the 'other' bond
            continue;
          }

          bool found = true;
          bool stereo;
          map<OBBond*, StereoRingBond>::iterator sb_it = _stereorbond.find(b);
          if (sb_it == _stereorbond.end()) // Not a ring closure
            // True/False for "up/down if moved to before the double bond C"
            stereo = !(IsUp(b) ^ (b->GetNbrAtomIdx(dbl_bond_atoms[i]) < dbl_bond_atoms[i]->GetIdx())) ;
          else  {                                                               // Is a ring closure
            char bc_result = SetRingClosureStereo(sb_it->second, dbl_bond);
            if (bc_result)
              stereo = bc_result == 1 ? true : false;
            else
              found = false;
          }

          if (!found) { // This cannot be used as the stereo bond
            other_bond[i] = b; // Use this for the 'other' bond
            continue;
          }
          
          if (stereo_bond[i] == NULL) { // This is a first stereo bond
            stereo_bond[i] = b; // Use this for the 'stereo' bond
            bond_stereo[i] = stereo;
          }
          else {               // This is a second stereo bond
            if (stereo != bond_stereo[i]) { // Verify that the other stereo bond (on the same atom) has opposite stereo
              other_bond[i] = b; // Use this for the 'other' bond
            }
            else  {
              obErrorLog.ThrowError(__FUNCTION__, "Error in cis/trans stereochemistry specified for the double bond\n", obWarning);
              stereo_bond[i] = (OBBond*) NULL;
            }
          }
        }
      }

      if (stereo_bond[0] == NULL || stereo_bond[1] == NULL) continue; // No cis/trans

      // other_bond will contain NULLs if there are bonds to implicit hydrogens
      unsigned int second = (other_bond[0] == NULL) ? OBStereo::ImplicitRef : other_bond[0]->GetNbrAtom(a1)->GetId();
      unsigned int fourth = (other_bond[1] == NULL) ? OBStereo::ImplicitRef : other_bond[1]->GetNbrAtom(a2)->GetId();


      OBCisTransStereo *ct = new OBCisTransStereo(&mol);
      OBCisTransStereo::Config cfg;
      cfg.begin = a1->GetId();
      cfg.end = a2->GetId();

       // If bond_stereo[0]==bond_stereo[1], this means cis for stereo_bond[0] and stereo_bond[1].
      if (bond_stereo[0] == bond_stereo[1])
        cfg.refs = OBStereo::MakeRefs(stereo_bond[0]->GetNbrAtom(a1)->GetId(), second,
                                      fourth, stereo_bond[1]->GetNbrAtom(a2)->GetId());
      else
        cfg.refs = OBStereo::MakeRefs(stereo_bond[0]->GetNbrAtom(a1)->GetId(), second,
                                      stereo_bond[1]->GetNbrAtom(a2)->GetId(), fourth);
      ct->SetConfig(cfg);
      // add the data to the atom
      mol.SetData(ct);
    }
  }

  void OBSmilesParser::InsertTetrahedralRef(OBMol &mol, unsigned long id)
  {
    map<OBAtom*, OBTetrahedralStereo::Config*>::iterator ChiralSearch;
    ChiralSearch = _tetrahedralMap.find(mol.GetAtom(_prev));
    if (ChiralSearch != _tetrahedralMap.end() && ChiralSearch->second != NULL)
    {
      int insertpos = NumConnections(ChiralSearch->first, id == OBStereo::ImplicitRef) - 2; // -1 indicates "from"
      if (insertpos > 2)
        return;
      if (insertpos < 0) {
        if (ChiralSearch->second->from != OBStereo::NoRef)
          obErrorLog.ThrowError(__FUNCTION__, "Warning: Overwriting previous from reference id.", obWarning);

        (ChiralSearch->second)->from = id;
        // cerr << "Adding " << id << " at Config.from to " << ChiralSearch->second << endl;
      } else {
        if (ChiralSearch->second->refs[insertpos] != OBStereo::NoRef)
          obErrorLog.ThrowError(__FUNCTION__, "Warning: Overwriting previously set reference id.", obWarning);

        (ChiralSearch->second)->refs[insertpos] = id;
        // cerr << "Adding " << id << " at " << insertpos << " to " << ChiralSearch->second << endl;
      }
    }
  }

  void OBSmilesParser::InsertSquarePlanarRef(OBMol &mol, unsigned long id)
  {
    map<OBAtom*, OBSquarePlanarStereo::Config*>::iterator ChiralSearch;
    ChiralSearch = _squarePlanarMap.find(mol.GetAtom(_prev));
    if (ChiralSearch != _squarePlanarMap.end() && ChiralSearch->second != NULL)
    {
      int insertpos = NumConnections(ChiralSearch->first) - 1;
      switch(insertpos) {
      case -1:
        if (ChiralSearch->second->refs[0] != OBStereo::NoRef)
          obErrorLog.ThrowError(__FUNCTION__, "Warning: Overwriting previous from reference id.", obWarning);
        (ChiralSearch->second)->refs[0] = id;
        break;
      case 0: case 1: case 2: case 3:
        if (ChiralSearch->second->refs[insertpos] != OBStereo::NoRef)
          obErrorLog.ThrowError(__FUNCTION__, "Warning: Overwriting previously set reference id.", obWarning);
        (ChiralSearch->second)->refs[insertpos] = id;
        break;
      default:
        obErrorLog.ThrowError(__FUNCTION__, "Warning: Square planar stereo specified for atom with more than 4 connections.", obWarning);
        break;
      }
    }
  }

  bool OBSmilesParser::ParseSimple(OBMol &mol)
  {
    int element;
    bool arom=false;

    switch(*_ptr)
      {
      case '*':
        element = 0;
        arom = false;
        break;

      case 'C':
        _ptr++;
        if (*_ptr == 'l')
          {
            element = 17;
          }
        else
          {
            element = 6;
            _ptr--;
          }
        break;

      case 'N':
        element = 7;
        break;
      case 'O':
        element = 8;
        break;
      case 'S':
        element = 16;
        break;
      case 'P':
        element = 15;
        break;
      case 'F':
        element = 9;
        break;
      case 'I':
        element = 53;
        break;

      case 'B':
        _ptr++;
        if (*_ptr == 'r')
          {
            element = 35;
          }
        else
          {
            element = 5;
            _ptr--;
          }
        break;

      // aromatics
      case 'b':
        arom = true;
        element = 5;
        break;
      case 'c':
        arom = true;
        element = 6;
        break;
      case 'n':
        arom = true;
        element = 7;
        break;
      case 'o':
        arom = true;
        element = 8;
        break;
      case 'p':
        arom = true;
        element = 15;
        break;
      case 's':
        arom = true;
        element = 16;
        break;
      default:
        {
        std::string err;
        err += "SMILES string contains a character '";
        err += *_ptr;
        err += "' which is invalid";
        obErrorLog.ThrowError(__FUNCTION__,
          err, obError);
        return false;
        }
      }

    OBAtom *atom = mol.NewAtom();
    atom->SetAtomicNum(element);
    if (_rxnrole > 1) { // Quick test for reaction
      // Set reaction role
      OBPairInteger *pi = new OBPairInteger();
      pi->SetAttribute("rxnrole");
      pi->SetValue(_rxnrole);
      atom->SetData(pi);
    }

    if (arom)
      atom->SetAromatic();

    if (_prev) //need to add bond
      {
        OBAtom* prevatom = mol.GetAtom(_prev);
        assert(prevatom);
        if (arom && prevatom->IsAromatic() && _order == 0)
          mol.AddBond(_prev, mol.NumAtoms(), 1, OB_AROMATIC_BOND); // this will be kekulized later
        else
          mol.AddBond(_prev, mol.NumAtoms(), _order == 0 ? 1 : _order);
        // store up/down
        if (_updown == BondUpChar || _updown == BondDownChar)
          _upDownMap[mol.GetBond(_prev, mol.NumAtoms())] = _updown;

        InsertTetrahedralRef(mol, mol.NumAtoms() - 1);
        InsertSquarePlanarRef(mol, mol.NumAtoms() - 1);
      }

    //set values
    _prev = mol.NumAtoms();
    _order = 0; // the default is that no bond symbol has been seen
    _updown = ' ';

    _hcount.push_back(-1); // implicit hydrogen count

    return(true);
  }

  bool OBSmilesParser::ParseComplex(OBMol &mol)
  {
    int element=0;
    bool arom=false;

    _ptr++;

    // Parse isotope information
    // - we parse anything with 1 to 4 digits
    // - any bigger and we risk overflowing the short int used to
    //   store the isotope information (max 65536)
    int isotope = 0;
    unsigned int size = 0;
    for (;*_ptr && isdigit(*_ptr) && size < 5;_ptr++) {
      isotope *= 10;
      isotope += *_ptr - '0';
      size++;
    }
    if (size == 5)
      return false;

    //parse element data
    switch(*_ptr)
    {
      case '*':
        element = 0;
        break;
    
      case 'C':
        _ptr++;
        switch(*_ptr)
          {
          case 'a':
            element = 20;
            break;
          case 'd':
            element = 48;
            break;
          case 'e':
            element = 58;
            break;
          case 'f':
            element = 98;
            break;
          case 'l':
            element = 17;
            break;
          case 'm':
            element = 96;
            break;
          case 'n':
            element = 112;
            break;
          case 'o':
            element = 27;
            break;
          case 'r':
            element = 24;
            break;
          case 's':
            element = 55;
            break;
          case 'u':
            element = 29;
            break;
          default:
            element =  6;
            _ptr--;
          }
        break;

      case 'N':
        _ptr++;
        switch(*_ptr)
          {
          case 'a':
            element =  11;
            break;
          case 'b':
            element =  41;
            break;
          case 'd':
            element =  60;
            break;
          case 'e':
            element =  10;
            break;
          case 'h':
            element = 113;
            break;
          case 'i':
            element =  28;
            break;
          case 'o':
            element = 102;
            break;
          case 'p':
            element =  93;
            break;
          default:
            element =   7;
            _ptr--;
          }
        break;

      case 'O':
        _ptr++;
        switch(*_ptr) {
        case 'g':
          element = 118;
          break;
        case 's':
          element = 76;
          break;
        default:
          element = 8;
          _ptr--;
        }
        break;

      case 'P':
        _ptr++;
        switch(*_ptr)
          {
          case 'a':
            element = 91;
            break;
          case 'b':
            element = 82;
            break;
          case 'd':
            element = 46;
            break;
          case 'm':
            element = 61;
            break;
          case 'o':
            element = 84;
            break;
          case 'r':
            element = 59;
            break;
          case 't':
            element = 78;
            break;
          case 'u':
            element = 94;
            break;
          default:
            element = 15;
            _ptr--;
          }
        break;

      case('S'):
        _ptr++;
        switch(*_ptr)
          {
          case 'b':
            element = 51;
            break;
          case 'c':
            element = 21;
            break;
          case 'e':
            element = 34;
            break;
          case 'g':
            element = 106;
            break;
          case 'i':
            element = 14;
            break;
          case 'm':
            element = 62;
            break;
          case 'n':
            element = 50;
            break;
          case 'r':
            element = 38;
            break;
          default:
            element = 16;
            _ptr--;
          }
        break;

      case 'B':
        _ptr++;
        switch(*_ptr)
          {
          case 'a':
            element = 56;
            break;
          case 'e':
            element =  4;
            break;
          case 'h':
            element =  107;
            break;
          case 'i':
            element = 83;
            break;
          case 'k':
            element = 97;
            break;
          case 'r':
            element = 35;
            break;
          default:
            element = 5;
            _ptr--;
          }
        break;

      case 'F':
        _ptr++;
        switch(*_ptr)
          {
          case 'e':
            element = 26;
            break;
          case 'l':
            element = 114;
            break;
          case 'm':
            element = 100;
            break;
          case 'r':
            element = 87;
            break;
          default:
            element = 9;
            _ptr--;
          }
        break;

      case 'I':
        _ptr++;
        switch(*_ptr)
          {
          case 'n':
            element = 49;
            break;
          case 'r':
            element = 77;
            break;
          default:
            element = 53;
            _ptr--;
          }
        break;

      case 'A':
        _ptr++;
        switch(*_ptr)
          {
          case 'c':
            element = 89;
            break;
          case 'g':
            element = 47;
            break;
          case 'l':
            element = 13;
            break;
          case 'm':
            element = 95;
            break;
          case 'r':
            element = 18;
            break;
          case 's':
            element = 33;
            break;
          case 't':
            element = 85;
            break;
          case 'u':
            element = 79;
            break;
          default:
            return(false);
          }
        break;

      case 'D':
        _ptr++;
        switch(*_ptr)
          {
          case 'b':
            element = 105;
            break;
          case 's':
            element = 110;
            break;
          case 'y':
            element = 66;
            break;
          default:
            return(false);
          }
        break;

      case 'E':
        _ptr++;
        switch(*_ptr)
          {
          case 'r':
            element = 68;
            break;
          case 's':
            element = 99;
            break;
          case 'u':
            element = 63;
            break;
          default:
            return(false);
          }
        break;

      case 'G':
        _ptr++;
        switch (*_ptr)
          {
          case 'a':
            element = 31;
            break;
          case 'd':
            element = 64;
            break;
          case 'e':
            element = 32;
            break;
          default:
            return(false);
          }
        break;

      case 'H':
        _ptr++;
        switch (*_ptr)
          {
          case 'e':
            element =  2;
            break;
          case 'f':
            element = 72;
            break;
          case 'g':
            element = 80;
            break;
          case 'o':
            element = 67;
            break;
          case 's':
            element = 108;
            break;
          default:
            element = 1;
            _ptr--;
          }
        break;

      case 'K':
        _ptr++;
        if(*_ptr == 'r')
          {
            element = 36;
          }
        else
          {
            element = 19;
            _ptr--;
          }
        break;

      case 'L':
        _ptr++;
        switch(*_ptr)
          {
          case 'a':
            element =  57;
            break;
          case 'i':
            element =   3;
            break;
          case 'r':
            element = 103;
            break;
          case 'u':
            element =  71;
            break;
          case 'v':
            element = 116;
            break;
          default:
            return(false);
          }
        break;

      case 'M':
        _ptr++;
        switch(*_ptr)
          {
          case 'c':
            element = 115;
            break;
          case 'd':
            element = 101;
            break;
          case 'g':
            element =  12;
            break;
          case 'n':
            element =  25;
            break;
          case 'o':
            element =  42;
            break;
          case 't':
            element =  109;
            break;
          default:
            return(false);
          }
        break;

      case 'R':
        _ptr++;
        switch(*_ptr)
          {
          case 'a':
            element = 88;
            break;
          case 'b':
            element = 37;
            break;
          case 'e':
            element = 75;
            break;
          case 'f':
            element = 104;
            break;
          case 'g':
            element = 111;
            break;
          case 'h':
            element = 45;
            break;
          case 'n':
            element = 86;
            break;
          case 'u':
            element = 44;
            break;
          default:
            return(false);
          }
        break;

      case 'T':
        _ptr++;
        switch(*_ptr)
          {
          case 'a':
            element = 73;
            break;
          case 'b':
            element = 65;
            break;
          case 'c':
            element = 43;
            break;
          case 'e':
            element = 52;
            break;
          case 'h':
            element = 90;
            break;
          case 'i':
            element = 22;
            break;
          case 'l':
            element = 81;
            break;
          case 'm':
            element = 69;
            break;
          case 's':
            element = 117;
            break;
          default:
            return(false);
          }
        break;

      case('U'):  element = 92;
        break;
      case('V'):  element = 23;
        break;
      case('W'):  element = 74;
        break;

      case('X'):
        _ptr++;
        if (*_ptr == 'e')
          {
            element = 54;
          }
        else
          {
            return(false);
          }
        break;

      case('Y'):
        _ptr++;
        if (*_ptr == 'b')
          {
            element = 70;
          }
        else
          {
            element = 39;
            _ptr--;
          }
        break;

      case('Z'):
        _ptr++;
        switch(*_ptr)
          {
          case 'n':
            element = 30;
            break;
          case 'r':
            element = 40;
            break;
          default:
            return false;
          }
        break;

      case 'a':
        _ptr++;
        if (*_ptr == 's') {
          arom = true;
          element = 33;
        }
        else
          return false;
        break;

      case 'b':
        _ptr++;
        if (*_ptr == 'i') {
          arom = true;
          element = 83;
        }
        else {
          arom = true;
          element = 5;
          _ptr--;
        }
        break;

      case 'c':
        arom = true;
        element = 6;
        break;

      case 'g':
        _ptr++;
        if (*_ptr == 'e') {
          arom = true;
          element = 32;
        }
        else
          return false;
        break;

      case 'n':
        arom = true;
        element = 7;
        break;

      case 'o':
        arom = true;
        element = 8;
        break;

      case 'p':
        arom = true;
        element = 15;
        break;

      case 's':
        arom = true;
        _ptr++;
        switch (*_ptr) {
        case 'e':
          element = 34;
          break;
        case 'i':
          element = 14;
          break;
        case 'n':
          element = 50;
          break;
        case 'b':
          element = 51;
          break;
        default:
          element = 16;
          _ptr--;
        }
        break;

      case 't':
        _ptr++;
        if (*_ptr == 'e') {
          arom = true;
          element = 52;
        }
        else
          return false;
        break;

      case '#':
        // Only support three digits for this extension
        if ((_ptr[1] == '1' || _ptr[1] == '2') &&
            (_ptr[2] >= '0' && _ptr[2] <= '9') &&
            (_ptr[3] >= '0' && _ptr[3] <= '9')) {
          element = (_ptr[1]-'0')*100 + (_ptr[2]-'0')*10 + (_ptr[3]-'0');
          if (element > 255) {
            std::string err = "Element number must be <= 255)";
            obErrorLog.ThrowError(__FUNCTION__,
              err, obError);
            return false;
          }
          _ptr += 3;
          break;
        }
        /* fall through to default */

      default:
        {
          std::string err;
          err += "SMILES string contains a character '";
          err += *_ptr;
          err += "' which is invalid";
          obErrorLog.ThrowError(__FUNCTION__,
            err, obError);
          return false;
        }
      }

    //handle hydrogen count, stereochemistry, and charge

    OBAtom *atom = mol.NewAtom();
    int hcount = 0;
    int charge=0;
    int rad=0;
    int clval=0;
    char tmpc[2];
    tmpc[1] = '\0';

    stringstream errorMsg;

    for (_ptr++;*_ptr && *_ptr != ']';_ptr++)
      {
        switch(*_ptr)
          {
          case '@':
            _ptr++;
            if (*_ptr == 'S' && _ptr[1] == 'P') { // @SP1/2/3
              // square planar atom found
              squarePlanarWatch = true;
              if (_squarePlanarMap.find(atom)==_squarePlanarMap.end()) // Prevent memory leak for malformed smiles (PR#3428432)
                _squarePlanarMap[atom] = new OBSquarePlanarStereo::Config;
              _squarePlanarMap[atom]->refs = OBStereo::Refs(4, OBStereo::NoRef);
              _squarePlanarMap[atom]->center = atom->GetId();
              _ptr += 2;
              switch(*_ptr) {
              case '1':
                _squarePlanarMap[atom]->shape = OBStereo::ShapeU; break;
              case '2':
                _squarePlanarMap[atom]->shape = OBStereo::Shape4; break;
              case '3':
                _squarePlanarMap[atom]->shape = OBStereo::ShapeZ; break;
              default:
                obErrorLog.ThrowError(__FUNCTION__, "Square planar stereochemistry must be one of SP1, SP2 or SP3", obWarning);
                return false;
              }
            } else {
              // tetrahedral atom found
              chiralWatch=true;
              if (_tetrahedralMap.find(atom)==_tetrahedralMap.end()) // Prevent memory leak for malformed smiles (PR#3428432)
                _tetrahedralMap[atom] = new OBTetrahedralStereo::Config;
              _tetrahedralMap[atom]->refs = OBStereo::Refs(3, OBStereo::NoRef);
              _tetrahedralMap[atom]->center = atom->GetId();
              if (*_ptr == '@') {
                _tetrahedralMap[atom]->winding = OBStereo::Clockwise;
              } else if (*_ptr == '?') {
                _tetrahedralMap[atom]->specified = false;
              } else {
                _tetrahedralMap[atom]->winding = OBStereo::AntiClockwise;
                _ptr--;
              }
            }
            break;
          case '-':
            _ptr++;
            if (!isdigit(*_ptr))
              charge--;
            while( isdigit(*_ptr) ) // go number by number
              charge = charge*10 - ((*_ptr++)-'0');
            _ptr--;
            break;
          case '+':
            _ptr++;
            if (!isdigit(*_ptr))
              charge++;
            while( isdigit(*_ptr) ) // go number by number
              charge = charge*10 + ((*_ptr++)-'0');
            _ptr--;
            break;
          case 'H':
            _ptr++;
            if (isdigit(*_ptr))
              {
                tmpc[0] = *_ptr;
                hcount = atoi(tmpc);
              }
            else
              {
                hcount = 1;
                _ptr--;
              }
            break;
          case '.': //CM Feb05
            rad=2;
            if(*(++_ptr)=='.')
              rad=3;
            else
              _ptr--;
            break;

          case ':':
            if(!isdigit(*(++_ptr)))
              {
                obErrorLog.ThrowError(__FUNCTION__,"The atom class following : must be a number", obWarning);
                return false;
              }
            while( isdigit(*_ptr) && clval < 100000000)
              clval = clval*10 + ((*_ptr++)-'0');
            --_ptr;
            { // a block is needed here to scope the OBPairInteger assignment
              OBPairInteger *atomclass = new OBPairInteger();
              atomclass->SetAttribute("Atom Class"); 
              atomclass->SetValue(clval);
              atomclass->SetOrigin(fileformatInput);
              atom->SetData(atomclass);
            }
            break;

          default:
            return(false);
          }
      }

    if (!*_ptr || *_ptr != ']')
      return(false); // we should have a trailing ']' now

    if (charge) {
      atom->SetFormalCharge(charge);
      if (abs(charge) > 10 || (element && charge > element)) { // if the charge is +/- 10 or more than the number of electrons
        errorMsg << "Atom " << atom->GetIdx() << " had an unrealistic charge of " << charge 
                 << "." << endl;
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
      }
    }
    if (rad)
      atom->SetSpinMultiplicity(rad);
    atom->SetAtomicNum(element);
    atom->SetIsotope(isotope);
    if (arom)
      atom->SetAromatic();
    if (_rxnrole > 1) { // Quick test for reaction
      // Set reaction role
      OBPairInteger *pi = new OBPairInteger();
      pi->SetAttribute("rxnrole");
      pi->SetValue(_rxnrole);
      atom->SetData(pi);
    }

    if (_prev) //need to add bond
      {
        OBAtom* prevatom = mol.GetAtom(_prev);
        if (arom && prevatom->IsAromatic() && _order == 0)
          mol.AddBond(_prev, mol.NumAtoms(), 1, OB_AROMATIC_BOND); // this will be kekulized later
        else
          mol.AddBond(_prev, mol.NumAtoms(), _order == 0 ? 1 : _order);
        // store up/down
        if (_updown == BondUpChar || _updown == BondDownChar)
          _upDownMap[mol.GetBond(_prev, mol.NumAtoms())] = _updown;

        if(chiralWatch) { // if tetrahedral atom, set previous as from atom
          _tetrahedralMap[atom]->from = mol.GetAtom(_prev)->GetId();
          if (CanHaveLonePair(element)) // Handle chiral lone pair as in X[S@@](Y)Z
            _chiralLonePair[mol.NumAtoms()] = 1; // First of the refs

          //cerr <<"NB7: line 1622: Added atom ref "<<_prev<<" at " << 0 << " to "<<_mapcd[atom]<<endl;
        }
        if (squarePlanarWatch) { // if squareplanar atom, set previous atom as first ref
          _squarePlanarMap[atom]->refs[0] = mol.GetAtom(_prev)->GetId();
          //cerr <<"TV7: line 1748: Added atom ref " << mol.GetAtom(_prev)->GetId()
          //     << " at " << 0 << " to " << _squarePlanarMap[atom] << endl;
        }
        InsertTetrahedralRef(mol, atom->GetId());
        InsertSquarePlanarRef(mol, atom->GetId());
      }
    else
      {
        // Handle chiral lone pair as in [S@@](X)(Y)Z
        if (chiralWatch && CanHaveLonePair(element)) // Handle chiral lone pair (only S at the moment)
          _chiralLonePair[mol.NumAtoms()] = 0; // 'from' atom
      }

    //set values
    _prev = mol.NumAtoms();
    _order = 0;
    _updown = ' ';
    
    if (hcount > 0) {
      if (chiralWatch)
        InsertTetrahedralRef(mol, OBStereo::ImplicitRef);
      if (squarePlanarWatch)
        InsertSquarePlanarRef(mol, OBStereo::ImplicitRef);
    }
    _hcount.push_back(hcount);

    chiralWatch=false;
    squarePlanarWatch = false;
    return(true);
  }

  bool OBSmilesParser::CapExternalBonds(OBMol &mol)
  {
    if (_extbond.empty())
      return true;

    OBAtom *atom;
    vector<ExternalBond>::iterator bond;
    for (bond = _extbond.begin(); bond != _extbond.end(); ++bond) {
      // create new dummy atom
      atom = mol.NewAtom();
      atom->SetAtomicNum(0);

      // bond dummy atom to mol via external bond
      mol.AddBond(bond->prev, atom->GetIdx(), bond->order);
      // store up/down
      if (bond->updown == BondUpChar || bond->updown == BondDownChar)
        _upDownMap[mol.GetBond(bond->prev, atom->GetIdx())] = bond->updown;

      OBBond *refbond = atom->GetBond(mol.GetAtom(bond->prev));

      //record external bond information
      OBExternalBondData *xbd;
      if (mol.HasData(OBGenericDataType::ExternalBondData)) {
        xbd = (OBExternalBondData*) mol.GetData(OBGenericDataType::ExternalBondData);
      } else {
        xbd = new OBExternalBondData;
        xbd->SetOrigin(fileformatInput);
        mol.SetData(xbd);
      }
      xbd->SetData(atom,refbond, bond->digit);
      //this data gets cleaned up in mol.Clear.
    }

    return true;
  }

  bool OBSmilesParser::ParseExternalBond(OBMol &mol)
  {
    int digit;
    char str[10];

    //*_ptr should == '&'
    _ptr++;

    switch (*_ptr) // check for bond order indicators CC&=1.C&1
      {
      case '-':
        _order = 1;
        _ptr++;
        break;
      case '=':
        _order = 2;
        _ptr++;
        break;
      case '#':
        _order = 3;
        _ptr++;
        break;
      case '$':
        _order = 4;
        _ptr++;
        break;
      case ';':
        _order = 5;
        _ptr++;
        break;
      case '/': //chiral, but _order still == 1
        _order = 1;
        _updown = BondDownChar;
        _ptr++;
        break;
      case '\\': // chiral, but _order still == 1
        _order = 1;
        _updown = BondUpChar;
        _ptr++;
        break;
      default: // no bond indicator just leave order = 0
        break;
      }

    if (*_ptr == '%') // external bond indicator > 10
      {
        _ptr++;
        str[0] = *_ptr;
        _ptr++;
        str[1] = *_ptr;
        str[2] = '\0';
      }
    else // simple single digit external bond indicator
      {
        str[0] = *_ptr;
        str[1] = '\0';
      }
    digit = atoi(str);  // convert indicator to digit

    //check for dot disconnect closures
    vector<ExternalBond>::iterator bond;
    int upDown, bondOrder;
    for (bond = _extbond.begin(); bond != _extbond.end(); ++bond) {

      if (bond->digit == digit) {
        upDown = (_updown > bond->updown) ? _updown : bond->updown;
        bondOrder = (_order > bond->order) ? _order : bond->order;
        mol.AddBond(bond->prev, _prev, bondOrder);
        // store up/down
        if (upDown == BondUpChar || upDown == BondDownChar)
          _upDownMap[mol.GetBond(bond->prev, _prev)] = upDown;


        // after adding a bond to atom "_prev"
        // search to see if atom is bonded to a chiral atom
        InsertTetrahedralRef(mol, bond->prev - 1);
        InsertSquarePlanarRef(mol, bond->prev - 1);

        _extbond.erase(bond);
        _updown = ' ';
        _order = 0;
        return true;
      }
    }

    //since no closures save another ext bond
    ExternalBond extBond;
    extBond.digit  = digit;
    extBond.prev   = _prev;
    extBond.order  = _order;
    extBond.updown = _updown;

    _extbond.push_back(extBond);
    _order = 0;
    _updown = ' ';

    return(true);

  }

  bool OBSmilesParser::ParseRingBond(OBMol &mol)
  {
    // The ring closure must be associated with a 'prev' atom
    OBAtom* prevatom = mol.GetAtom(_prev);
    if (!prevatom) {
      obErrorLog.ThrowError(__FUNCTION__,"Number not parsed correctly as a ring bond", obWarning);
      return false;
    } 

    // Parse the ring closure digit or digits
    int digit = 0;
    if (*_ptr == '%') {
      _ptr++;
      if (*_ptr == '(') { // %(NNN) extension to OpenSMILES
        _ptr++;
        const char* start = _ptr;
        while (isdigit(*_ptr)) {
          digit *= 10;
          digit += *_ptr - '0';
          _ptr++;
          if (_ptr - start > 5) {
            obErrorLog.ThrowError(__FUNCTION__, "Ring closure numbers with more than 5 digits are not supported", obWarning);
            return false;
          }
        }
        if (*_ptr != ')') {
          obErrorLog.ThrowError(__FUNCTION__, "Matching close parenthesis not found for ring closure number", obWarning);
          return false;
        }
      }
      else { // % followed by two-digit ring closure
        if (!isdigit(*_ptr) || !isdigit(*(_ptr + 1))) {
          obErrorLog.ThrowError(__FUNCTION__, "Two digits expected after %", obWarning);
          return false;
        }
        digit = (*_ptr - '0') * 10 + *(_ptr+1) - '0';
        _ptr++;
      }
    }
    else {
      digit = *_ptr - '0';
    }

    vector<RingClosureBond>::iterator bond;
    int upDown, bondOrder;
    for (bond = _rclose.begin(); bond != _rclose.end(); ++bond) {
      if (bond->digit == digit) {
        // Check for self-bonding, e.g. C11
        if (bond->prev == _prev) {
          obErrorLog.ThrowError(__FUNCTION__, "Invalid SMILES: Ring closures imply atom bonded to itself.", obWarning);
          return false;
        }
        upDown = (_updown > bond->updown) ? _updown : bond->updown;
        bondOrder = (_order > bond->order) ? _order : bond->order;
        // Check if this ring closure bond may be aromatic and set order accordingly
        bool aromatic_bond = false;
        if (bondOrder == 0) {
          OBAtom *a1 = mol.GetAtom(bond->prev);
          OBAtom *a2 = mol.GetAtom(_prev);
          if (a1->IsAromatic() && a2->IsAromatic())
            aromatic_bond = true;
        }
        mol.AddBond(bond->prev, _prev, bondOrder == 0 ? 1 : bondOrder, aromatic_bond ? OB_AROMATIC_BOND : 0, bond->numConnections);
        // store up/down
        if (upDown == BondUpChar || upDown == BondDownChar)
          _upDownMap[mol.GetBond(bond->prev, _prev)] = upDown;

        // For assigning cis/trans in the presence of bond closures, we need to
        // remember info on all bond closure bonds.
        StereoRingBond sb;
        sb.updown.push_back(_updown);
        sb.atoms.push_back(mol.GetAtom(_prev));
        sb.updown.push_back(bond->updown);
        sb.atoms.push_back(mol.GetAtom(bond->prev));
        _stereorbond[mol.GetBond(bond->prev, _prev)] = sb; // Store for later
        
        // after adding a bond to atom "_prev"
        // search to see if atom is bonded to a chiral atom
        // need to check both _prev and bond->prev as closure is direction independent
        InsertTetrahedralRef(mol, bond->prev - 1);
        InsertSquarePlanarRef(mol, bond->prev - 1);

        // FIXME: needed for squreplanar too??
        map<OBAtom*, OBTetrahedralStereo::Config*>::iterator ChiralSearch;
        ChiralSearch = _tetrahedralMap.find(mol.GetAtom(bond->prev));
        if (ChiralSearch != _tetrahedralMap.end() && ChiralSearch->second != NULL) {
          int insertpos = bond->numConnections - 1;
          switch(insertpos) {
          case -1:
            if (ChiralSearch->second->from != OBStereo::NoRef)
              obErrorLog.ThrowError(__FUNCTION__, "Warning: Overwriting previous from reference id.", obWarning);
            (ChiralSearch->second)->from = mol.GetAtom(_prev)->GetId();
            break;
          case 0: case 1: case 2:
            if (ChiralSearch->second->refs[insertpos] != OBStereo::NoRef)
              obErrorLog.ThrowError(__FUNCTION__, "Warning: Overwriting previously set reference id.", obWarning);
            (ChiralSearch->second)->refs[insertpos] = mol.GetAtom(_prev)->GetId();
            break;
          default:
            obErrorLog.ThrowError(__FUNCTION__, "Warning: Tetrahedral stereo specified for atom with more than 4 connections.", obWarning);
            break;
          }
        }

        //CM ensure neither atoms in ring closure is a radical centre
        OBAtom* patom = mol.GetAtom(_prev);
        patom->SetSpinMultiplicity(0);
        patom = mol.GetAtom(bond->prev);
        patom->SetSpinMultiplicity(0);
        //CM end
        _rclose.erase(bond);
        _updown = ' ';
        _order = 0;
        return true;
      }
    }

    //since no closures save another rclose bond
    RingClosureBond ringClosure;
    ringClosure.digit  = digit;
    ringClosure.prev   = _prev;
    ringClosure.order  = _order;
    ringClosure.updown = _updown;

    OBAtom* atom = mol.GetAtom(_prev);
    ringClosure.numConnections = NumConnections(atom); //store position to insert closure bond
    _rclose.push_back(ringClosure);
    _order = 0;
    _updown = ' ';

    return(true);
  }

  // NumConnections finds the number of connections already made to
  // a particular atom. This is used to figure out the correct position
  // to insert an atom ID into atom4refs
  int OBSmilesParser::NumConnections(OBAtom *atom, bool isImplicitRef)
  {
    int val = atom->GetExplicitDegree();
    // The implicit H is not included in "val" so we need to adjust by 1
    if (isImplicitRef)
      return val+1;

    int idx = atom->GetIdx();
    // Need to adjust for any implicit H (e.g. [C@@H]) but only for atoms after the H.
    // The following line controls for this. It uses the fact the _hcount is only set
    // after this function is called to handle inserting the stereo ref for the implicit H.
    if (idx-1 < _hcount.size() && _hcount[idx-1] > 0)
      val += _hcount[idx-1];
    vector<RingClosureBond>::iterator bond;
    //correct for multiple closure bonds to a single atom
    for (bond = _rclose.begin(); bond != _rclose.end(); ++bond)
      if (bond->prev == idx)
        val++;

    return val;
  }


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
    for (i = _child_nodes.begin();i != _child_nodes.end();++i)
      delete (*i);
  }

  struct OutOptions
  {
    bool isomeric;
    bool kekulesmi;
    bool showatomclass;
    bool showexplicitH;
    bool smarts;
    const char* ordering; // This is a pointer to the string in the original map
    OutOptions(bool _isomeric, bool _kekulesmi, bool _showatomclass, bool _showexplicitH, bool _smarts,
               const char* _ordering):
      isomeric(_isomeric), kekulesmi(_kekulesmi), showatomclass(_showatomclass), showexplicitH(_showexplicitH),
      smarts(_smarts),
      ordering(_ordering)
      {}
  };

  /*----------------------------------------------------------------------
   * CLASS OBMol2Cansmi - Declarations
   ----------------------------------------------------------------------*/

  class OBMol2Cansmi
  {
    std::vector<int> _atmorder;
    OBBitVec _uatoms,_ubonds;
    std::vector<OBBondClosureInfo> _vopen;
    unsigned int _bcdigit; // Unused unless option "R" is specified
    std::vector<OBCisTransStereo> _cistrans, _unvisited_cistrans;
    std::map<OBBond *, bool> _isup;

    bool          _canonicalOutput; // regular or canonical SMILES

    OBMol* _pmol;
    OBStereoFacade *_stereoFacade;
    OBConversion* _pconv;

    OBAtom* _endatom;
    OBAtom* _startatom;

    OutOptions &options;

  public:
    OBMol2Cansmi(OutOptions &_options): options(_options)
    {
    }
    ~OBMol2Cansmi()
    {
      delete _stereoFacade;
    }

    void         Init(OBMol* pmol, bool canonicalOutput = true, OBConversion* pconv=NULL);

    void         CreateCisTrans(OBMol&);
    char         GetCisTransBondSymbol(OBBond *, OBCanSmiNode *);
    bool         AtomIsChiral(OBAtom *atom);
    bool         BuildCanonTree(OBMol &mol, OBBitVec &frag_atoms,
                                vector<unsigned int> &canonical_order,
                                OBCanSmiNode *node);
    void         CreateFragCansmiString(OBMol&, OBBitVec&, std::string&);
    const char * GetTetrahedralStereo(OBCanSmiNode*,
                                      vector<OBAtom*>&chiral_neighbors,
                                      vector<unsigned int> &symmetry_classes);
    const char*  GetSquarePlanarStereo(OBCanSmiNode*,
                                       vector<OBAtom*>&chiral_neighbors,
                                       vector<unsigned int> &symmetry_classes);
    bool         GetSmilesElement(OBCanSmiNode*,
                                  vector<OBAtom*>&chiral_neighbors,
                                  vector<unsigned int> &symmetry_classes,
                                  std::string&);
    int          GetSmilesValence(OBAtom *atom);
    int          GetUnusedIndex();
    vector<OBBondClosureInfo>
    GetCanonClosureDigits(OBAtom *atom,
                          OBBitVec &frag_atoms,
                          vector<unsigned int> &canonical_order);
    bool         IsSuppressedHydrogen(OBAtom *atom);
    void         ToCansmilesString(OBCanSmiNode *node,
                                   std::string &buffer,
                                   OBBitVec &frag_atoms,
                                   vector<unsigned int> &symmetry_classes,
                                   vector<unsigned int> &canonical_order);
    bool         HasStereoDblBond(OBBond *, OBAtom *atom);
    void MyFindChildren(OBMol &mol, vector<OBAtom*> &children, OBBitVec &seen, OBAtom *end);
    void GetOutputOrder(std::string &outorder);
    bool         ParseInChI(OBMol &mol, vector<int> &atom_order);
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

  void OBMol2Cansmi::Init(OBMol* pmol, bool canonical, OBConversion* pconv)
  {
    _atmorder.clear();
    _atmorder.reserve(pmol->NumAtoms());
    _uatoms.Clear();
    _ubonds.Clear();
    _vopen.clear();

    _pmol = pmol;
    _stereoFacade = new OBStereoFacade(_pmol); // needs to be destroyed in dtor
    _pconv = pconv;
    _canonicalOutput = canonical;

    _endatom = NULL;
    _startatom = NULL;
  }


  /***************************************************************************
   * FUNCTION: GetUnusedIndex
   *
   * DESCRIPTION:
   *       Returns the next available bond-closure index for a SMILES.
   *
   *       You could just do this sequentially, not reusing bond-closure
   *       digits, thus (chosen by Option("R")):
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
    if (_pconv->IsOption("R")) {
      // Keep incrementing the bond closure digits (for each connected component)
      _bcdigit++;
      return _bcdigit;
    }

    int idx=1;
    vector<OBBondClosureInfo>::iterator j;
    for (j = _vopen.begin();j != _vopen.end();)
      if (j->ringdigit == idx)
        {
          idx++; //increment idx and start over if digit is already used
          j = _vopen.begin();
        }
      else ++j;

    return(idx);
  }

  void OBMol2Cansmi::CreateCisTrans(OBMol &mol)
  {
    std::vector<OBGenericData*> vdata = mol.GetAllData(OBGenericDataType::StereoData);
    for (std::vector<OBGenericData*>::iterator data = vdata.begin(); data != vdata.end(); ++data) {
      if (((OBStereoBase*)*data)->GetType() != OBStereo::CisTrans)
        continue;
      OBCisTransStereo *ct = dynamic_cast<OBCisTransStereo*>(*data);
      if (ct && ct->GetConfig().specified) {
        OBCisTransStereo::Config config = ct->GetConfig();
        OBBond* dbl_bond = mol.GetBond(mol.GetAtomById(config.begin), mol.GetAtomById(config.end));
        if (!dbl_bond)
          continue;
        // Do not output cis/trans bond symbols for double bonds in rings of size IMPLICIT_CIS_RING_SIZE or less
        unsigned int boundedringsize = OBBondGetSmallestRingSize(dbl_bond, IMPLICIT_CIS_RING_SIZE);
        if (boundedringsize == 0) // either not in ring at all, or not in small ring
          _cistrans.push_back(*ct);
      }
    }

    _unvisited_cistrans = _cistrans; // Make a copy of _cistrans
  }

  bool OBMol2Cansmi::HasStereoDblBond(OBBond *bond, OBAtom *atom)
  {
    // This is a helper function for determining whether to
    // consider writing a cis/trans bond symbol for bond closures.
    // Returns TRUE only if the atom is connected to the cis/trans
    // double bond. To handle the case of conjugated bonds, one must
    // remember that the ring opening preceded the closure, so if the
    // ring opening bond was on a stereocenter, it got the symbol already.

    if (!bond || !atom)
      return false;

    std::vector<OBCisTransStereo>::iterator ChiralSearch;
    OBAtom *nbr_atom = bond->GetNbrAtom(atom);

    bool stereo_dbl = false;
    if (atom->HasDoubleBond()) {
      stereo_dbl = true;
      if (nbr_atom->HasDoubleBond())
        // Check whether the nbr_atom is a begin or end in any CisTransStereo. If so,
        // then the ring opening already had the symbol.
        for (ChiralSearch = _cistrans.begin(); ChiralSearch != _cistrans.end(); ++ChiralSearch) {
          OBCisTransStereo::Config cfg = ChiralSearch->GetConfig();
          if (nbr_atom->GetId() == cfg.begin || nbr_atom->GetId() == cfg.end) {
            // I don't think I need to check whether it has a bond with atom
            stereo_dbl = false;
            break;
          }
        }
    }
    return stereo_dbl;
  }

  char OBMol2Cansmi::GetCisTransBondSymbol(OBBond *bond, OBCanSmiNode *node)
  {
    // Given a cis/trans bond and the node in the SMILES tree, figures out
    // whether to write a '/' or '\' symbol.
    // See the comments smilesformat.cpp: FixCisTransBonds().
    //
    // The OBCanSmiNode is the most-recently-written atom in the SMILES string
    // we're creating.  If it is the double-bonded atom, then the substituent
    // follows, so that "up" means '/' and "down" means '\'.  If the OBCanSmiNode
    // atom is the single-bonded atom then the double-bonded atom comes next,
    // in which case "up" means '\' and "down" means '/'.
    //
    // Note that the story is not so simple for conjugated systems where
    // we need to take into account what symbol was already used.

    if (!bond /*|| (!bond->IsUp() && !bond->IsDown())*/)
      return '\0';
    OBAtom *atom = node->GetAtom();
    OBAtom *nbr_atom = bond->GetNbrAtom(atom);
    OBMol *mol = atom->GetParent();

    // If this bond is in two different obcistransstereos (e.g. a conjugated system)
    // choose the one where the dbl bond atom is *atom (i.e. the one which comes first)
    std::vector<OBCisTransStereo>::iterator ChiralSearch;
    std::vector<unsigned long>::iterator lookup;

    bool dbl_bond_first = false;
    if (atom->HasDoubleBond())
    {
      if (nbr_atom->HasDoubleBond())
        // Check whether the atom is a center in any CisTransStereo. If so,#
        // then this CisTransStereo takes precedence over any other
        for (ChiralSearch = _cistrans.begin(); ChiralSearch != _cistrans.end(); ++ChiralSearch)
        {
          OBCisTransStereo::Config cfg = ChiralSearch->GetConfig();
          if (atom->GetId() == cfg.begin || atom->GetId() == cfg.end) {
            // I don't think I need to check whether it has a bond with nbr_atom
            dbl_bond_first = true;
            break;
          }
        }
      else
        dbl_bond_first = true;
    }

    // Has the symbol for this bond already been set?
    if (_isup.find(bond) == _isup.end()) // No it hasn't
    {
      unsigned int endatom, centeratom;
      if (dbl_bond_first) {
        if (atom->IsAromatic())
          FOR_BONDS_OF_ATOM (bond, atom)
            if (bond->IsAromatic() && bond->GetBondOrder()==2)
              return 0;
        endatom = nbr_atom->GetId();
        centeratom = atom->GetId();
      }
      else {
        if (nbr_atom->IsAromatic())
          FOR_BONDS_OF_ATOM (bond, nbr_atom)
            if (bond->IsAromatic() && bond->GetBondOrder() == 2)
              return 0;
        endatom = atom->GetId();
        centeratom = nbr_atom->GetId();
      }

      for (ChiralSearch = _unvisited_cistrans.begin(); ChiralSearch != _unvisited_cistrans.end(); ++ChiralSearch)
      {
        OBCisTransStereo::Config cfg = ChiralSearch->GetConfig(OBStereo::ShapeU);
        lookup = std::find(cfg.refs.begin(), cfg.refs.end(), endatom);
        if (lookup != cfg.refs.end() && (cfg.begin == centeratom || cfg.end == centeratom))
        { // Atoms endatom and centeratom are in this OBCisTransStereo

          std::vector<OBBond *> refbonds(4, (OBBond*)NULL);
          refbonds[0] = mol->GetBond(mol->GetAtomById(cfg.refs[0]), mol->GetAtomById(cfg.begin));

          if (cfg.refs[1] != OBStereo::ImplicitRef) // Could be a hydrogen
            refbonds[1] = mol->GetBond(mol->GetAtomById(cfg.refs[1]), mol->GetAtomById(cfg.begin));

          if (cfg.refs[2] != OBStereo::ImplicitRef) // Could be a hydrogen
            refbonds[2] = mol->GetBond(mol->GetAtomById(cfg.refs[2]), mol->GetAtomById(cfg.end));

          if (cfg.refs[3] != OBStereo::ImplicitRef) // Could be a hydrogen
            refbonds[3] = mol->GetBond(mol->GetAtomById(cfg.refs[3]), mol->GetAtomById(cfg.end));

          // What symbol would the four refs use if before the dbl bond?
          bool config[4] = {true, false, false, true};
          bool use_same_config = true;
          // (The actual config used will be config ^ use_same_config)

          // Make sure that the symbol for this bond is true. This ensures
          // a canonical string, so that it's always C/C=C/C and not C\C=C\C.
          for(int i=0;i<4;i++)
            if (refbonds[i] == bond)
              if (!config[i])
              {
                  use_same_config = false;
                  break;
              }
          // If any of the bonds have been previously set, now set them all
          // in the opposite sense
          for(int i=0;i<4;i++)
            if (_isup.find(refbonds[i]) != _isup.end()) // We have already set this one (conjugated bond)
              if (_isup[refbonds[i]] == (config[i] ^ use_same_config))
              {
                use_same_config = !use_same_config;
                break;
              }
          // Set the configuration
          for(int i=0;i<4;i++)
            if (refbonds[i] != NULL)
              _isup[refbonds[i]] = config[i] ^ use_same_config;
          _unvisited_cistrans.erase(ChiralSearch);
          break; // break out of the ChiralSearch
        }
      }
    }

     // If ChiralSearch didn't find the bond, we can't set this symbol
    if (_isup.find(bond) == _isup.end()) return '\0';

    if (dbl_bond_first) { // double-bonded atom is first in the SMILES
      if (_isup[bond])
        return '/';
      else
        return '\\';
    }
    else { // double-bonded atom is second in the SMILES
      if (_isup[bond])
        return '\\';
      else
        return '/';
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
                                      std::string &buffer)
  {
    char symbol[10];
    symbol[0] = '\0'; // make sure to initialize for all paths below

    bool bracketElement = false;
    bool normalValence = true;
    bool writeExplicitHydrogen = false;

    OBAtom *atom = node->GetAtom();
    int element = atom->GetAtomicNum();

    // Handle SMILES Valence model, and explicit and implicit hydrogens
    if (IsOutsideOrganicSubset(element))
      bracketElement = true;

    unsigned int numExplicitHsToSuppress = 0;
    // Don't suppress any explicit Hs attached if the atom is an H itself (e.g. [H][H]) or -xh was specified
    if (atom->GetAtomicNum() != OBElements::Hydrogen && !options.showexplicitH) {
      FOR_NBORS_OF_ATOM(nbr, atom) {
        if (nbr->GetAtomicNum() == OBElements::Hydrogen && (!options.isomeric || nbr->GetIsotope() == 0) && nbr->GetExplicitDegree() == 1 &&
          nbr->GetFormalCharge() == 0 && (!options.showatomclass || !nbr->GetData("Atom Class")))
          numExplicitHsToSuppress++;
      }
    }

    unsigned int numImplicitHs = 0;
    if (options.smarts) {
      if (numExplicitHsToSuppress > 0) {
        bracketElement = true;
        numImplicitHs = numExplicitHsToSuppress;
      }
    }
    else {
      numImplicitHs = atom->GetImplicitHCount() + numExplicitHsToSuppress;
      if (!bracketElement) {
        if (element == 0) { // asterisk is always hypervalent but we don't bracket it unless has Hs
          if (numImplicitHs > 0)
            bracketElement = true;
        }
        else {
          int bosum = atom->GetExplicitValence() - numExplicitHsToSuppress;
          unsigned int implicitValence = SmilesValence(element, bosum, false);
          unsigned int defaultNumImplicitHs = implicitValence - bosum;
          if (implicitValence == 0 // hypervalent
             ||  numImplicitHs != defaultNumImplicitHs // undervalent
             || (!options.kekulesmi && element != 6 && atom->IsAromatic() && numImplicitHs != 0) ) // aromatic nitrogen/phosphorus
            bracketElement = true;
        }
      }
    }

    if (atom->GetFormalCharge() != 0 // charged elements
      || (options.isomeric && atom->GetIsotope()) // isotopes
      || (options.showatomclass && atom->HasData("Atom Class")) ) // If the molecule has Atom Class data and -xa option set and atom has data
      bracketElement = true;

    const char* stereo = (const char*)0;
    if (GetSmilesValence(atom) > 2 && options.isomeric) {
      stereo = GetTetrahedralStereo(node, chiral_neighbors, symmetry_classes);
      if (stereo == (const char*)0)
        stereo = GetSquarePlanarStereo(node, chiral_neighbors, symmetry_classes);
    }
    if (stereo != (const char*)0)
      bracketElement = true;

    if (!bracketElement) {

      // Ordinary non-bracket element
      if (element) {
        const char* symbol = OBElements::GetSymbol(atom->GetAtomicNum());
        if ((!options.kekulesmi && atom->IsAromatic()) || // aromatic atom
            (atom->GetSpinMultiplicity() && _pconv->IsOption("r"))) //Radical centres lowercase if r option set
        {
          buffer += symbol[0] + ('a' - 'A');
          if (symbol[1])
            buffer += symbol[1];
        }
        else
          buffer += symbol;
      }

      // Atomic number zero - either '*' or an external atom
      else {
        bool external = false;
        vector<pair<int,pair<OBAtom *,OBBond *> > > *externalBonds =
          (vector<pair<int,pair<OBAtom *,OBBond *> > > *)((OBMol*)atom->GetParent())->GetData("extBonds");
        vector<pair<int,pair<OBAtom *,OBBond *> > >::iterator externalBond;

        if (externalBonds) // TODO: This code has bit-rotted and needs some love
          for(externalBond = externalBonds->begin();externalBond != externalBonds->end();++externalBond) {
            if (externalBond->second.first == atom) {
              external = true;
              buffer += '&';
              OBBond *bond = externalBond->second.second;
              //if (bond->IsUp()) {
              //  if ( (bond->GetBeginAtom())->HasDoubleBond() ||
              //       (bond->GetEndAtom())->HasDoubleBond() )
              //    buffer += '\\';
              //}
              //if (bond->IsDown()) {
              //  if ( (bond->GetBeginAtom())->HasDoubleBond() ||
              //       (bond->GetEndAtom())->HasDoubleBond() )
              //    buffer += '/';
              //}
              if (bond->GetBondOrder() == 2 && !bond->IsAromatic()) // TODO: need to check for kekulesmi
                buffer += '=';
              if (bond->GetBondOrder() == 2 && bond->IsAromatic())
                buffer += ':';
              if (bond->GetBondOrder() == 3)
                buffer += '#';
              if (bond->GetBondOrder() == 4)
                buffer += '$';
              char tmp[10];
              snprintf(tmp, 10, "%d", externalBond->first);
              buffer += tmp;
              break;
            }
          }

        if(!external)
          buffer += '*';
      }

      return true;
    }

    // Bracketed atoms, e.g. [Pb], [OH-], [C@]
    buffer += '[';
    unsigned short iso = atom->GetIsotope();
    if (options.isomeric && iso) {
      if (iso >= 10000) // max 4 characters
        obErrorLog.ThrowError(__FUNCTION__, "Isotope value larger than 9999. Ignoring value.", obWarning);
      else {
        char iso[8]; // 7 characters plus null
        snprintf(iso, 8, "%u", atom->GetIsotope());
        buffer += iso;
      }
    }
    if (!atom->GetAtomicNum())
      buffer += '*';
    else {
      if (atom->GetAtomicNum() == OBElements::Hydrogen && options.smarts)
        buffer += "#1";
      else {
        unsigned int elem = atom->GetAtomicNum();
        const char* symbol = OBElements::GetSymbol(elem);
        if (*symbol == '\0') {
          char atomnum[8];  // '#' plus 3 digits plus null
          snprintf(atomnum, 8, "#%u", elem);
          buffer += atomnum;
        } else if (!options.kekulesmi && atom->IsAromatic()) { // aromatic atom
          buffer += symbol[0] + ('a' - 'A');
          if (symbol[1])
            buffer += symbol[1];
        }
        else
          buffer += symbol;
      }
    }

    // If chiral, append '@' or '@@'...unless we're creating a SMARTS ("s") and it's @H or @@H
    if (stereo != (const char*)0 && !(options.smarts && atom->GetImplicitHCount() > 0))
      buffer += stereo;

    // Add extra hydrogens.
    int hcount = numImplicitHs;
    if (hcount > 0 && (atom == _endatom || atom == _startatom)) // Leave a free valence for attachment
      hcount--;
    if (hcount > 0) {
      if (options.smarts && stereo == (const char*)0) {
        char tcount[10];
        for (int i = 0; i < hcount; ++i) {
          buffer += "!H";
          snprintf(tcount, 10, "%d", i);
          buffer += tcount;
        }
      }
      else {
        buffer += 'H';
        if (hcount > 1) {
          char tcount[10];
          snprintf(tcount, 10, "%d", hcount);
          buffer += tcount;
        }
      }
    }

    // Append charge to the end
    int charge = atom->GetFormalCharge();
    if (charge != 0) {
      if (charge > 0)
        buffer += '+';
      else
        buffer += '-';

      if (abs(charge) > 1) {
        char tchar[10];
        snprintf(tchar, 10, "%d", abs(charge));
        buffer += tchar;
      }
    }

    //atom class e.g. [C:2]
    if (options.showatomclass) {
      OBGenericData *data = atom->GetData("Atom Class");
      if (data) {
        OBPairInteger* acdata = dynamic_cast<OBPairInteger*>(data); // Could replace with C-style cast if willing to live dangerously
        if (acdata) {
          int ac = acdata->GetGenericValue();
          if (ac >= 0) { // Allow 0, why not?
            buffer += ':';
            char tchar[12]; // maxint has 10 digits
            snprintf(tchar, 12, "%d", ac);
            buffer += tchar;
          }
        }
      }
    }

    buffer += ']';

    return true;
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
    unsigned long atomid = atom->GetId();
    return _stereoFacade->HasTetrahedralStereo(atomid) || _stereoFacade->HasSquarePlanarStereo(atomid);
  }

  /***************************************************************************
   * FUNCTION: GetTetrahedralStereo
   *
   * DESCRIPTION:
   *       If the atom is chiral, return either "@" or "@@". Otherwise 0.
   ***************************************************************************/

  const char* OBMol2Cansmi::GetTetrahedralStereo(OBCanSmiNode *node,
                                          vector<OBAtom*> &chiral_neighbors,
                                          vector<unsigned int> &symmetry_classes)
  {
    // If not enough chiral neighbors were passed in, we're done
    if (chiral_neighbors.size() < 4)
      return (const char*)0;

    OBAtom *atom = node->GetAtom();
    OBTetrahedralStereo *ts = _stereoFacade->GetTetrahedralStereo(atom->GetId());
    // If atom is not a tetrahedral center, we're done
    if (!ts)
      return (const char*)0;

    // get the Config struct defining the stereochemistry
    OBTetrahedralStereo::Config atomConfig = ts->GetConfig();

    // Unspecified or unknown stereochemistry
    if (!atomConfig.specified || (atomConfig.specified && atomConfig.winding==OBStereo::UnknownWinding))
      return (const char*)0;

    // create a Config struct with the chiral_neighbors in canonical output order
    OBStereo::Refs canonRefs;
    for (vector<OBAtom*>::const_iterator atom_it = chiral_neighbors.begin() + 1; atom_it != chiral_neighbors.end(); ++atom_it) {
      if (*atom_it)
        canonRefs.push_back((*atom_it)->GetId());
      else // Handle a chiral lone pair, represented by a NULL OBAtom* in chiral_neighbors
        canonRefs.push_back(OBStereo::ImplicitRef);
    }
    OBTetrahedralStereo::Config canConfig;
    canConfig.center = atom->GetId();
    if (chiral_neighbors[0])
      canConfig.from = chiral_neighbors[0]->GetId();
    else // Handle a chiral lone pair, represented by a NULL OBAtom* in chiral_neighbors
      canConfig.from = OBStereo::ImplicitRef;
    canConfig.refs = canonRefs;

    // canConfig is clockwise
    if (atomConfig == canConfig)
      return "@@";
    else
      return "@";
  }

  /***************************************************************************
   * FUNCTION: GetSquarePlanarStereo
   *
   * DESCRIPTION:
   *       If the atom is chiral, return either '@SP1', '@SP2' or '@SP'.
   *       Otherwise, return 0.
   ***************************************************************************/

  const char* OBMol2Cansmi::GetSquarePlanarStereo(OBCanSmiNode *node,
                                           vector<OBAtom*> &chiral_neighbors,
                                           vector<unsigned int> &symmetry_classes)
  {
    // If no chiral neighbors were passed in, we're done
    if (chiral_neighbors.size() < 4)
      return (const char*)0;
    
    OBAtom *atom = node->GetAtom();

    OBSquarePlanarStereo *sp = _stereoFacade->GetSquarePlanarStereo(atom->GetId());
    // If atom is not a square-planar center, we're done
    if (!sp)
      return (const char*)0;

    // get the Config struct defining the stereochemistry
    OBSquarePlanarStereo::Config atomConfig = sp->GetConfig();

    if (!atomConfig.specified)
      return (const char*)0;

    // create a Config struct with the chiral_neighbors in canonical output order
    OBStereo::Refs canonRefs = OBStereo::MakeRefs(chiral_neighbors[0]->GetId(),
        chiral_neighbors[1]->GetId(), chiral_neighbors[2]->GetId(), chiral_neighbors[3]->GetId());
    OBSquarePlanarStereo::Config canConfig;
    canConfig.center = atom->GetId();
    canConfig.refs = canonRefs;

    // canConfig is U shape
    if (atomConfig == canConfig)
      return "@SP1";

    canConfig.shape = OBStereo::Shape4;
    if (atomConfig == canConfig)
      return "@SP2";

    canConfig.shape = OBStereo::ShapeZ;
    if (atomConfig == canConfig)
      return "@SP3";

    return (const char*)0;
  }

  //! Adaptation of OBMol::FindChildren to allow a vector of OBAtoms to be passed in
  //  MOVE THIS TO OBMOL FOR 2.4
  void OBMol2Cansmi::MyFindChildren(OBMol &mol, vector<OBAtom*> &children, OBBitVec &seen, OBAtom *end)
  {
    OBBitVec curr,next;

    OBBitVec used(seen);
    used |= end->GetIdx();
    curr |= end->GetIdx();
    children.clear();

    int i;
    OBAtom *atom,*nbr;
    vector<OBBond*>::iterator j;

    for (;;)
      {
        next.Clear();
        for (i = curr.NextBit(-1);i != curr.EndBit();i = curr.NextBit(i))
          {
            atom = mol.GetAtom(i);
            for (nbr = atom->BeginNbrAtom(j);nbr;nbr = atom->NextNbrAtom(j))
              if (!used[nbr->GetIdx()])
                {
                  children.push_back(nbr);
                  next |= nbr->GetIdx();
                  used |= nbr->GetIdx();
                }
          }
        if (next.IsEmpty())
          break;
        curr = next;
      }
  }

  // Do we need to write out a bond symbol for this bond?
  // No - if it's aromatic
  // Otherwise, yes if the bond order is not 1
  // If the bond order *is* 1, then only if the bond is in a ring and between aromatic atoms
  static bool NeedsBondSymbol(OBBond* bond)
  {
    if (bond->IsAromatic())
      return false;
    switch (bond->GetBondOrder()) {
    case 1:
      if (bond->IsInRing() && bond->GetBeginAtom()->IsAromatic() && bond->GetEndAtom()->IsAromatic())
        return true;
      return false;
    default: // bond orders != 1
      return true;
    }
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
    vector<OBBond*>::iterator i;
    OBAtom *nbr, *atom;
    vector<OBAtom *> sort_nbrs;
    vector<OBAtom *>::iterator ai;
    OBBond *bond;
    OBCanSmiNode *next;
    int idx;

    atom = node->GetAtom();

#if DEBUG
    cout << "BuildCanonTree: " << OBElements::GetSymbol(atom->GetAtomicNum()) << ", " << atom->GetIdx() << ", canorder " << canonical_order[atom->GetIdx()-1] << "\n";
#endif

    // Create a vector of neighbors sorted by canonical order, but favor
    // double and triple bonds over single and aromatic.  This causes
    // ring-closure digits to avoid double and triple bonds.
    //
    // Since there are typically just one to three neighbors, we just do a
    // ordered insertion rather than sorting.

    bool favor_multiple = true; // Visit 'multiple' bonds first
    if (options.ordering)
      favor_multiple = false; // Visit in strict canonical order (if using user-specified order)

    for (nbr = atom->BeginNbrAtom(i); nbr; nbr = atom->NextNbrAtom(i)) {

      idx = nbr->GetIdx();
      //if (nbr->GetAtomicNum() == OBElements::Hydrogen && IsSuppressedHydrogen(nbr)) {
      //  _uatoms.SetBitOn(nbr->GetIdx());        // mark suppressed hydrogen, so it won't be considered
      //  continue;                               // later when looking for more fragments.
      //}
      if (_uatoms[idx] || !frag_atoms.BitIsSet(idx))
        continue;

      OBBond *nbr_bond = atom->GetBond(nbr);
      unsigned int nbr_bond_order = nbr_bond->GetBondOrder();
      int new_needs_bsymbol = NeedsBondSymbol(nbr_bond);

      for (ai = sort_nbrs.begin(); ai != sort_nbrs.end(); ++ai) {
        bond = atom->GetBond(*ai);
        unsigned int bond_order = bond->GetBondOrder();
        int sorted_needs_bsymbol = NeedsBondSymbol(bond);
        if (favor_multiple && new_needs_bsymbol && !sorted_needs_bsymbol) {
          sort_nbrs.insert(ai, nbr);
          ai = sort_nbrs.begin();//insert invalidated ai; set it to fail next test
          break;
        }
        if (   (!favor_multiple || new_needs_bsymbol == sorted_needs_bsymbol)
               && canonical_order[idx-1] < canonical_order[(*ai)->GetIdx()-1]) {
          sort_nbrs.insert(ai, nbr);
          ai = sort_nbrs.begin();//insert invalidated ai; set it to fail next test
          break;
        }
      }
      if (ai == sort_nbrs.end())
        sort_nbrs.push_back(nbr);
    }

    _uatoms.SetBitOn(atom->GetIdx());     //mark the atom as visited

    if (_endatom && !_uatoms.BitIsSet(_endatom->GetIdx()) && sort_nbrs.size() > 1) {
      // If you have specified an _endatom, the following section rearranges
      // sort_nbrs as follows:
      //   - if a branch does not lead to the end atom, move it to the front
      //     (i.e. visit it first)
      //   - otherwise move it to the end
      // This section is skipped if sort_nbrs has only a single member, or if
      // we have already visited _endatom.
      
      vector<OBAtom*> children;
      MyFindChildren(mol, children, _uatoms, _endatom);
      
      vector<OBAtom*> front, end;
      for (vector<OBAtom *>::iterator it=sort_nbrs.begin(); it!=sort_nbrs.end(); ++it)
        if (std::find(children.begin(), children.end(), *it) == children.end() && *it != _endatom)
          front.push_back(*it);
        else
          end.push_back(*it);
      sort_nbrs = front;
      sort_nbrs.insert(sort_nbrs.end(), end.begin(), end.end());
    }

    // Build the next layer of nodes, in canonical order
    for (ai = sort_nbrs.begin(); ai != sort_nbrs.end(); ++ai) {
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
    vector<OBBond*>::iterator i;
    OBBond *bond1, *bond2;
    OBAtom *nbr1, *nbr2;
    int nbr1_canorder, nbr2_canorder;

    vp_closures.clear();
    vbonds.clear();

    // Find new ring-closure bonds for this atom
    for (bond1 = atom->BeginBond(i); bond1; bond1 = atom->NextBond(i)) {

      // Is this a ring-closure neighbor?
      if (_ubonds.BitIsSet(bond1->GetIdx()))
        continue;
      nbr1 = bond1->GetNbrAtom(atom);
      // Skip hydrogens before checking canonical_order
      // PR#1999348
      if (   (nbr1->GetAtomicNum() == OBElements::Hydrogen && IsSuppressedHydrogen(nbr1))
             || !frag_atoms.BitIsSet(nbr1->GetIdx()))
        continue;

      nbr1_canorder = canonical_order[nbr1->GetIdx()-1];

      // Insert into the bond-vector in canonical order (by neighbor atom order)
      for (bi = vbonds.begin(); bi != vbonds.end(); ++bi) {
        bond2 = *bi;
        nbr2 = bond2->GetNbrAtom(atom);
        nbr2_canorder = canonical_order[nbr2->GetIdx()-1];
        if (nbr1_canorder < nbr2_canorder) {
          vbonds.insert(bi, bond1);
          bi = vbonds.begin();//insert invalidated bi; set it to fail next test
          break;
        }
      }
      if (bi == vbonds.end())     // highest one (or first one) - append to end
        vbonds.push_back(bond1);
    }

    // If we found new open bonds, assign a bond-closure digits to each one,
    // add it to _vopen, and add it to the return vector.
    for (bi = vbonds.begin(); bi != vbonds.end(); ++bi) {
      bond1 = *bi;
      _ubonds.SetBitOn(bond1->GetIdx());
      int digit = GetUnusedIndex();
      int bo = (bond1->IsAromatic())? 1 : bond1->GetBondOrder();  // CJ: why was this line added?  bo is never used?
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
          ++j;
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
    if (atom->GetExplicitDegree() != 1)          // not exactly one bond
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
   *       This is like GetHvyDegree(), but it returns the "valence" of an
   *       atom as it appears in the SMILES string.  In particular, hydrogens
   *       count if they will appear explicitly -- see IsSuppressedHydrogen()
   *       above.
   ***************************************************************************/

  int OBMol2Cansmi::GetSmilesValence(OBAtom *atom)
  {
    int count = 0;

    if (atom->GetAtomicNum() == OBElements::Hydrogen)
      return atom->GetExplicitDegree();

    if (options.showexplicitH)
      return atom->GetExplicitDegree();

    FOR_NBORS_OF_ATOM(nbr, atom) {
      if (nbr->GetAtomicNum() != OBElements::Hydrogen
            || nbr->GetIsotope() != 0
            || nbr->GetExplicitDegree() != 1)
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
                                       std::string &buffer,
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

    // Don't include chiral symbol on _endatom or _startatom.
    // Otherwise, we end up with C[C@@H](Br)(Cl), where the C has 4 neighbours already
    // and we cannot concatenate another SMILES string without creating a 5-valent C.

    bool is_chiral = AtomIsChiral(atom);
    if (is_chiral && atom!=_endatom && atom!=_startatom) {

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
      // that will appear explicitly in the SMILES as a separate atom is
      // treated like any other atom when calculating the chirality.)

      if (!options.showexplicitH) {
        FOR_NBORS_OF_ATOM(i_nbr, atom) {
          OBAtom *nbr = &(*i_nbr);
          if (nbr->GetAtomicNum() == OBElements::Hydrogen && IsSuppressedHydrogen(nbr)) {
            chiral_neighbors.push_back(nbr);
            break;        // quit loop: only be one H if atom is chiral
          }
        }
      }

      // Handle implict H by adding a NULL OBAtom*
      if(atom->GetImplicitHCount() == 1)
        chiral_neighbors.push_back(static_cast<OBAtom*> (NULL));

      // Ok, done with H. Now we need to consider the case where there is a chiral
      // lone pair. If it exists (and we won't know for sure until we've counted up
      // all the neighbours) it will go in here
      int lonepair_location = chiral_neighbors.size();

      // Ok, done with all that. Next in the SMILES will be the ring-closure characters.
      // So we need to find the corresponding atoms and add them to the list.
      // (We got the canonical ring-closure list earlier.)
      if (!vclose_bonds.empty()) {
        vector<OBBondClosureInfo>::iterator i;
        for (i = vclose_bonds.begin();i != vclose_bonds.end();++i) {
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

      // Handle a chiral lone-pair on a sulfur, by inserting a NULL OBAtom* at the
      // appropriate location
      if (chiral_neighbors.size() == 3 && CanHaveLonePair(atom->GetAtomicNum())) // Handle sulfur
        chiral_neighbors.insert(chiral_neighbors.begin() + lonepair_location, static_cast<OBAtom*> (NULL));

    }

    // Write the current atom to the string
    GetSmilesElement(node, chiral_neighbors, symmetry_classes, buffer);

    _atmorder.push_back(atom->GetIdx());  //store the atom ordering

    // Write ring-closure digits
    if (!vclose_bonds.empty()) {
      vector<OBBondClosureInfo>::iterator bci;
      for (bci = vclose_bonds.begin(); bci != vclose_bonds.end(); ++bci) {
        if (!bci->is_open)
        { // Ring closure
          char bs[2] = { '\0', '\0' };
          // Only get symbol for ring closures on the dbl bond
          if (HasStereoDblBond(bci->bond, node->GetAtom()))
            bs[0] = GetCisTransBondSymbol(bci->bond, node);
          if (bs[0])
            buffer += bs;	// append "/" or "\"
          else
          {
            switch (bci->bond->GetBondOrder())
            {
            case 1:
              if (!bci->bond->IsAromatic() && bci->bond->IsInRing() && bci->bond->GetBeginAtom()->IsAromatic() && bci->bond->GetEndAtom()->IsAromatic())
                buffer += '-';
              break;
            case 2:
              if (options.kekulesmi || !bci->bond->IsAromatic())
                buffer += '=';
              break;
            case 3:
              buffer += '#'; break;
            case 4:
              buffer += '$'; break;
            }
          }
        }
        else
        { // Ring opening
          char bs[2] = { '\0', '\0' };
          // Only get symbol for ring openings on the dbl bond
          if (!HasStereoDblBond(bci->bond, bci->bond->GetNbrAtom(node->GetAtom())))
            bs[0] = GetCisTransBondSymbol(bci->bond, node);
          if (bs[0])
            buffer += bs;	// append "/" or "\"
        }
        if (bci->ringdigit > 9) {
          buffer += '%';
          if (bci->ringdigit > 99)
            buffer += '(';
          char digits[8];
          snprintf(digits, 8, "%d", bci->ringdigit);
          buffer += digits;
          if (bci->ringdigit > 99)
            buffer += ')';
        }
        else
          buffer += (char)(bci->ringdigit + '0');
      }
    }

    // Write child bonds, then recursively follow paths to child nodes
    // to print the SMILES for each child branch.

    OBBond *bond;
    for (int i = 0;i < node->Size();i++) {
      bond = node->GetChildBond(i);
      if (i+1 < node->Size() || node->GetAtom() == _endatom)
        buffer += '(';

      switch (bond->GetBondOrder()) {
      case 1:
        char cc[2];
        cc[0] = GetCisTransBondSymbol(bond, node);
        if (cc[0] != '\0') {
          cc[1] = '\0';
          buffer += cc;
        }
        else {
          // Write a single bond symbol if not aromatic but end atoms are both aromatic
          // This will speed up reading as it will avoid ring perception around line 563 (bond->IsInRing())
          // TODO: Consider making the test for IsInRing() an option
          if (!bond->IsAromatic() && bond->IsInRing() && bond->GetBeginAtom()->IsAromatic() && bond->GetEndAtom()->IsAromatic())
            buffer += '-';
        }
        break;
      case 2:
        if (options.kekulesmi || !bond->IsAromatic())
          buffer += '=';
        break;
      case 3:
        buffer += '#';
        break;
      case 4:
        buffer += '$';
      }

      ToCansmilesString(node->GetChildNode(i),buffer, frag_atoms, symmetry_classes, canonical_order);
      if (i+1 < node->Size() || node->GetAtom() == _endatom)
        buffer += ')';
    }
  }

  /****************************************************************************
   * FUNCTION: StandardLabels
   *
   * DESCRIPTION:
   *        Creates a set of non-canonical labels for the fragment atoms
   * ***************************************************************************/
  void StandardLabels(OBMol *pMol, OBBitVec *frag_atoms,
                      vector<unsigned int> &symmetry_classes,
                      vector<unsigned int> &labels)
  {
    FOR_ATOMS_OF_MOL(atom, *pMol) {
      if (frag_atoms->BitIsSet(atom->GetIdx())) {
        labels.push_back(atom->GetIdx() - 1);
        symmetry_classes.push_back(atom->GetIdx() - 1);
      }
      else{
        labels.push_back(OBStereo::ImplicitRef); //to match situation when canonical ordering. Just a big number?
        symmetry_classes.push_back(OBStereo::ImplicitRef);
      }
    }
  }

  /***************************************************************************
   * FUNCTION: RandomLabels
   *
   * DESCRIPTION:
   *    Creates a set of random labels for the fragment atoms.  Primarily
   *    for testing: you can create a bunch of random SMILES for the same
   *    molecule, and use those to test the canonicalizer.
   ***************************************************************************/

  void RandomLabels(OBMol *pMol, OBBitVec &frag_atoms,
      vector<unsigned int> &symmetry_classes,
      vector<unsigned int> &labels)
  {
    int natoms = pMol->NumAtoms();
    OBBitVec used(natoms);

    FOR_ATOMS_OF_MOL(atom, *pMol) {
      if (frag_atoms.BitIsSet(atom->GetIdx())) {
        int r = rand() % natoms;
        while (used.BitIsSet(r)) {
          r = (r + 1) % natoms;         // find an unused number
        }
        used.SetBitOn(r);
        labels.push_back(r);
        symmetry_classes.push_back(r);
      }
      else{
        labels.push_back(OBStereo::ImplicitRef); //to match situation when canonical ordering. Just a big number?
        symmetry_classes.push_back(OBStereo::ImplicitRef);
      }
    }
  }

  //! Same as tokenize, except in treatment of multiple delimiters. Tokenize
  //! treats multiple delimiters as a single delimiter. It also ignores delimiters
  //! in the first or last position. In contrast, mytokenize treats each instance of
  //! the delimiter as the end/start of a new token.
  bool mytokenize(std::vector<std::string> &vcr, std::string &s,
                      const char *delimstr)
  {
    vcr.clear();
    size_t startpos=0,endpos=0;

    size_t s_size = s.size();
    for (;;)
      {
        //startpos = s.find_first_not_of(delimstr,startpos);
        endpos = s.find_first_of(delimstr,startpos);
        if (endpos <= s_size && startpos <= s_size)
          {
            vcr.push_back(s.substr(startpos,endpos-startpos));
          }
        else
          {
            if (startpos <= s_size)
              vcr.push_back(s.substr(startpos,s_size-startpos));
            break;
          }

        startpos = endpos+1;
      }
    return(true);
  }

  // Returns canonical label order
  bool OBMol2Cansmi::ParseInChI(OBMol &mol, vector<int> &atom_order)
  {
    /*OBConversion MolConv(*_pconv); //new copy to use to write associated MOL
    MolConv.SetAuxConv(NULL); //temporary until a proper OBConversion copy constructor written

    OBFormat* pInChIFormat = _pconv->FindFormat("InChI");
    if(pInChIFormat==NULL) {
      obErrorLog.ThrowError(__FUNCTION__, "InChI format not available", obError);
      return false;
    }*/
    OBConversion MolConv;
    MolConv.SetOutFormat("InChI");
    MolConv.SetAuxConv(NULL); //temporary until a proper OBConversion copy constructor written
    stringstream newstream;
    MolConv.SetOutStream(&newstream);
    // I'm sure there's a better way of preventing InChI warning output
    MolConv.AddOption("w", OBConversion::OUTOPTIONS);
    MolConv.AddOption("a", OBConversion::OUTOPTIONS);
    MolConv.AddOption("X", OBConversion::OUTOPTIONS, "RecMet FixedH");
    //pInChIFormat->WriteMolecule(&mol, &MolConv);
    MolConv.Write(&mol);

    vector<string> splitlines;
    string tmp = newstream.str();
    tokenize(splitlines, tmp,"\n");
    vector<string> split, split_aux;
    string aux_part;

    size_t rm_start = splitlines.at(0).find("/r"); // Correct for reconnected metal if necessary
    if (rm_start == string::npos) {
      tokenize(split, splitlines.at(0),"/");
      aux_part = splitlines.at(1); // Use the normal labels
    }
    else { 
      tmp = splitlines.at(0).substr(rm_start);
      tokenize(split, tmp, "/");
      split.insert(split.begin(), "");
      size_t rm_start_b = splitlines.at(1).find("/R:");
      aux_part = splitlines.at(1).substr(rm_start_b); // Use the reconnected metal labels
    }
    tokenize(split_aux, aux_part, "/");

    // Parse the canonical labels
    vector<vector<int> > canonical_labels;
    vector<string> s_components, s_atoms;

    tmp = split_aux.at(2).substr(2);
    tokenize(s_components, tmp, ";");
    for(vector<string>::iterator it=s_components.begin(); it!=s_components.end(); ++it) {
      tokenize(s_atoms, *it, ",");
      vector<int> atoms;
      for(vector<string>::iterator itb=s_atoms.begin(); itb!=s_atoms.end(); ++itb)
        atoms.push_back(atoi(itb->c_str()));
      canonical_labels.push_back(atoms);
    }

    // Adjust the canonical labels if necessary using a /F section
    size_t f_start = aux_part.find("/F:");
    if (f_start != string::npos) {
      tmp = aux_part.substr(f_start+3);
      tokenize(split_aux, tmp, "/");
      tokenize(s_components, split_aux.at(0), ";");
      vector<vector<int> > new_canonical_labels;
      int total = 0;
      for(vector<string>::iterator it=s_components.begin(); it!=s_components.end(); ++it) {
        // e.g. "1,2,3;2m" means replace the first component by "1,2,3"
        //                       but keep the next two unchanged
        if (*(it->rbegin()) == 'm') {
          int mult;
          if (it->size()==1)
            mult = 1;
          else
            mult = atoi(it->substr(0, it->size()-1).c_str());
          new_canonical_labels.insert(new_canonical_labels.end(), 
            canonical_labels.begin()+total, canonical_labels.begin()+total+mult);
          total += mult;
        }
        else {
          tokenize(s_atoms, *it, ",");
          vector<int> atoms;
          for(vector<string>::iterator itb=s_atoms.begin(); itb!=s_atoms.end(); ++itb)
            atoms.push_back(atoi(itb->c_str()));
          new_canonical_labels.push_back(atoms);
          total++;
        }
      }
      canonical_labels = new_canonical_labels;
    }

    // Flatten the canonical_labels
    for(vector<vector<int> >::iterator it=canonical_labels.begin(); it!=canonical_labels.end(); ++it) {
      atom_order.insert(atom_order.end(), it->begin(), it->end());      
    }

    return true;
  }

  /**
   * Helper function for getFragment below.
   */
  void addNbrs(OBBitVec &fragment, OBAtom *atom, const OBBitVec &mask)
  {
    FOR_NBORS_OF_ATOM (nbr, atom) {
      if (!mask.BitIsSet(nbr->GetIdx()))
        continue;
      // skip visited atoms
      if (fragment.BitIsSet(nbr->GetIdx()))
        continue;
      // add the neighbor atom to the fragment
      fragment.SetBitOn(nbr->GetIdx());
      // recurse...
      addNbrs(fragment, &*nbr, mask);
    }
  }

  /**
   * Create an OBBitVec objects with bets set for the fragment consisting of all
   * atoms for which there is a path to atom without going through skip. These
   * fragment bitvecs are indexed by atom idx (i.e. OBAtom::GetIdx()).
   */
  OBBitVec getFragment(OBAtom *atom, const OBBitVec &mask)
  {
    OBBitVec fragment;
    fragment.SetBitOn(atom->GetIdx());
    // start the recursion
    addNbrs(fragment, atom, mask);
    return fragment;
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

  void OBMol2Cansmi::CreateFragCansmiString(OBMol &mol, OBBitVec &frag_atoms, std::string &buffer)
  {
    //cout << "CreateFragCansmiString()" << endl;
    OBAtom *atom;
    OBCanSmiNode *root;
    buffer[0] = '\0';
    vector<OBNodeBase*>::iterator ai;
    vector<unsigned int> symmetry_classes, canonical_order;
    symmetry_classes.reserve(mol.NumAtoms());
    canonical_order.reserve(mol.NumAtoms());

    // Remember the desired endatom, if specified
    const char* pp = _pconv->IsOption("l");
    unsigned int atom_idx  = pp ? atoi(pp) : 0;
    if (atom_idx >= 1 && atom_idx <= mol.NumAtoms())
      _endatom = mol.GetAtom(atom_idx);
    // Was a start atom specified?
    pp = _pconv->IsOption("f");
    atom_idx  = pp ? atoi(pp) : 0;
    if (atom_idx >= 1 && atom_idx <= mol.NumAtoms())
      _startatom = mol.GetAtom(atom_idx);

    // Was an atom ordering specified?
    const char* ppo = options.ordering;
    vector<string> s_atom_order;
    vector<int> atom_order;
    if (ppo) {
      tokenize(s_atom_order,ppo,"-()");
      if (s_atom_order.size() != mol.NumHvyAtoms())
        ppo = NULL;
      else {
        for (vector<string>::const_iterator cit=s_atom_order.begin(); cit!=s_atom_order.end(); ++cit)
          atom_order.push_back(atoi(cit->c_str()));
        atom_idx = atom_order.at(0);
        if (atom_idx >= 1 && atom_idx <= mol.NumAtoms())
          _startatom = mol.GetAtom(atom_idx);
      }
    }
    // Was Universal SMILES requested?
    bool universal_smiles = _pconv->IsOption("U");
    if (universal_smiles) {
      bool parsedOkay = ParseInChI(mol, atom_order);
      if (!parsedOkay)
        universal_smiles = false;
    }

    // First, create a canonical ordering vector for the atoms.  Canonical
    // labels are zero indexed, corresponding to "atom->GetIdx()-1".
    if (_canonicalOutput) {

      // Find the (dis)connected fragments.
      OBBitVec visited;
      std::vector<OBBitVec> fragments;
      for (std::size_t i = 0; i < mol.NumAtoms(); ++i) {
        if (!frag_atoms.BitIsSet(i+1) || visited.BitIsSet(i+1))
          continue;
        fragments.push_back(getFragment(mol.GetAtom(i+1), frag_atoms));
        visited |= fragments.back();
      }

      // Determine symmetry classes for each disconnected fragment seperatly
      symmetry_classes.resize(mol.NumAtoms());
      for (std::size_t i = 0; i < fragments.size(); ++i) {
        OBGraphSym gs(&mol, &(fragments[i]));
        vector<unsigned int> tmp;
        gs.GetSymmetry(tmp);

        for (std::size_t j = 0; j < mol.NumAtoms(); ++j)
          if (fragments[i].BitIsSet(j+1))
            symmetry_classes[j] = tmp[j];
      }

      /*
      OBGraphSym gs(&mol, &frag_atoms);
      gs.GetSymmetry(symmetry_classes);
      */

      // Was a canonicalization timeout given?
      unsigned int maxSeconds = 5;
      const char *timeoutString = _pconv->IsOption("T");
      if (timeoutString) {
        std::stringstream ss(timeoutString);
        if (!(ss >> maxSeconds)) {
          obErrorLog.ThrowError(__FUNCTION__, "Canonicalization timeout should be a number", obWarning);
          maxSeconds = 5;
        }
      }

      CanonicalLabels(&mol, symmetry_classes, canonical_order, frag_atoms, maxSeconds);
    }
    else {
      if (_pconv->IsOption("C")) {      // "C" == "anti-canonical form"
        RandomLabels(&mol, frag_atoms, symmetry_classes, canonical_order);
      } else if (ppo || universal_smiles) { // user-specified or InChI canonical labels
        canonical_order.resize(mol.NumAtoms());
        symmetry_classes.resize(mol.NumAtoms());
        int idx = 3; // Start the labels at 3 (to leave space for special values 0, 1 and 2)
        for (int i=0; i<atom_order.size(); ++i)
          if (canonical_order[atom_order[i] - 1] == 0) { // Ignore ring closures (for "U")
            canonical_order[atom_order[i] - 1] = idx;
            symmetry_classes[atom_order[i] - 1] = idx;
            ++idx;
          }
        for (int i=0; i<canonical_order.size(); ++i)
          if (canonical_order[i] == 0) { // Explicit hydrogens
            if (mol.GetAtom(i+1)->GetAtomicNum() == OBElements::Hydrogen && mol.GetAtom(i+1)->GetIsotope()!=0) { // [2H] or [3H]
              canonical_order[i] = mol.GetAtom(i+1)->GetIsotope() - 1; // i.e. 1 or 2
              symmetry_classes[i] = canonical_order[i];
            }
          }
      } else {
        StandardLabels(&mol, &frag_atoms, symmetry_classes, canonical_order);
      }
    }

    // OUTER LOOP: Handles dot-disconnected structures and reactions.  Finds the
    // lowest unmarked canorder atom in the current reaction role, and starts there
    // to generate a SMILES.
    // Repeats until no atoms remain unmarked.

    bool new_rxn_role = false; // flag to indicate whether we have started a new reaction role
    bool isrxn = mol.IsReaction();
    OBReactionFacade rxn(&mol);
    unsigned int rxnrole = 1; // reactants
    while (1) {
      if (_pconv->IsOption("R"))
        _bcdigit = 0; // Reset the bond closure index for each disconnected component

      // It happens that the lowest canonically-numbered atom is usually
      // a good place to start the canonical SMILES.
      OBAtom *root_atom;
      unsigned int lowest_canorder = 999999;
      root_atom = NULL;

      // If we specified a startatom_idx & it's in this fragment, use it to start the fragment
      if (_startatom)
        if (!_uatoms[_startatom->GetIdx()] && 
           frag_atoms.BitIsSet(_startatom->GetIdx()) && 
           (!isrxn || rxn.GetRole(_startatom)==rxnrole))
          root_atom = _startatom;

      if (root_atom == NULL) {
        for (atom = mol.BeginAtom(ai); atom; atom = mol.NextAtom(ai)) {
          int idx = atom->GetIdx();
          if (//atom->GetAtomicNum() != OBElements::Hydrogen       // don't start with a hydrogen
              !_uatoms[idx]          // skip atoms already used (for fragments)
              && frag_atoms.BitIsSet(idx)// skip atoms not in this fragment
              && (!isrxn || rxn.GetRole(atom)==rxnrole) // skip atoms not in this rxn role
              //&& !atom->IsChiral()    // don't use chiral atoms as root node
              && canonical_order[idx-1] < lowest_canorder) {
            root_atom = atom;
            lowest_canorder = canonical_order[idx-1];
          }
        }
        // For Inchified or Universal SMILES, if the start atom is an [O-] attached to atom X, choose any =O attached to X instead.
        //          Ditto for [S-] and =S.
        if ((_pconv->IsOption("I") || universal_smiles)
             && root_atom && root_atom->GetFormalCharge()==-1  && root_atom->GetExplicitDegree() == 1
             && root_atom->HasSingleBond() && (root_atom->GetAtomicNum() == OBElements::Oxygen || root_atom->GetAtomicNum() == OBElements::Sulfur)) {
          OBBondIterator bi = root_atom->BeginBonds();
          OBAtom* central = root_atom->BeginNbrAtom(bi);
          FOR_NBORS_OF_ATOM(nbr, central) {
            if (root_atom == &*nbr) continue;
            if (nbr->GetAtomicNum() == root_atom->GetAtomicNum() && nbr->GetExplicitDegree() == 1 && nbr->HasDoubleBond()) {
              root_atom = &*nbr;
              break;
            }
          }
        }
      }

      // No atom found?  If it's not a reaction, then we've done all fragments.
      // If it is, then increment the rxn role and try again.
      if (root_atom == NULL) {
        if (mol.IsReaction()) {
          rxnrole++;
          if (rxnrole==4)
            break;
          buffer += '>';
          new_rxn_role = true;
          continue;
        }
        else
          break;
      }

      // Clear out closures in case structure is dot disconnected
      //      _atmorder.clear();
      _vopen.clear();

      // Dot disconnected structure or new rxn role?
      if (new_rxn_role)
        new_rxn_role = false;
      else if (!buffer.empty())
        buffer += '.';

      root = new OBCanSmiNode (root_atom);

      BuildCanonTree(mol, frag_atoms, canonical_order, root);
      ToCansmilesString(root, buffer, frag_atoms, symmetry_classes, canonical_order);
      delete root;
    }
  }

  void OBMol2Cansmi::GetOutputOrder(std::string &outorder)
  {
    std::vector<int>::iterator it = _atmorder.begin();
    if (it != _atmorder.end()) {
      char tmp[15];
      snprintf(tmp, 15, "%d", *it);
      outorder += tmp;
      ++it;
      for (; it != _atmorder.end(); ++it) {
        snprintf(tmp, 15, "%d", *it);
        outorder += ' ';
        outorder += tmp;
      }
    }
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

  void CreateCansmiString(OBMol &mol, std::string &buffer, OBBitVec &frag_atoms, OBConversion* pConv)
  {
    bool canonical = pConv->IsOption("c") != NULL;

    OutOptions options(!pConv->IsOption("i"), pConv->IsOption("k"),
      pConv->IsOption("a"),
      pConv->IsOption("h"), pConv->IsOption("s"),
      pConv->IsOption("o"));

    OBMol2Cansmi m2s(options);
    m2s.Init(&mol, canonical, pConv);

    if (options.isomeric) {
      PerceiveStereo(&mol);
      m2s.CreateCisTrans(mol); // No need for this if not iso
    }
    else {
      // Not isomeric - be sure there are no Z coordinates, clear
      // all stereo-center and cis/trans information.
      OBBond *bond;
      vector<OBBond*>::iterator bi;
      for (bond = mol.BeginBond(bi); bond; bond = mol.NextBond(bi)) {
        bond->SetHash(false);
        bond->SetWedge(false);
      }
    }

    if (!options.showexplicitH) {
      // If the fragment includes explicit hydrogens, exclude them.
      // They won't appear in the SMILES anyway (unless they're attached to
      // a chiral center, or it's something like [H][H]).
      FOR_ATOMS_OF_MOL(iatom, mol) {
        OBAtom *atom = &(*iatom);
        if (frag_atoms.BitIsSet(atom->GetIdx()) && atom->GetAtomicNum() == OBElements::Hydrogen
          && (!options.isomeric || m2s.IsSuppressedHydrogen(atom))) {
          frag_atoms.SetBitOff(atom->GetIdx());
        }
      }
    }

    m2s.CreateFragCansmiString(mol, frag_atoms, buffer);

    if (pConv->IsOption("O")) { // record smiles atom order info
      // This atom order data is useful not just for canonical SMILES
      // Could also save canonical bond order if anyone desires
      OBPairData *canData;
      if (!mol.HasData("SMILES Atom Order")) {
        // Create new OBPairData
        canData = new OBPairData;
        canData->SetAttribute("SMILES Atom Order");
        canData->SetOrigin(OpenBabel::local);
        mol.SetData(canData);
      }
      else {
        // Recanonicalizing - update existing new OBPairData
        canData = (OBPairData *)mol.GetData("SMILES Atom Order");
      }
      std::string atmorder;
      m2s.GetOutputOrder(atmorder);
      canData->SetValue(atmorder);
    }
  }

  bool SMIBaseFormat::GetInchifiedSMILESMolecule(OBMol *mol, bool useFixedHRecMet)
  {
    OBConversion MolConv;

    OBFormat* pInChIFormat = MolConv.FindFormat("InChI");
    if(pInChIFormat==NULL) {
      obErrorLog.ThrowError(__FUNCTION__, "InChI format not available", obError);
      return false;
    }
    stringstream newstream;
    MolConv.SetOutStream(&newstream);
    if (useFixedHRecMet) {
      MolConv.AddOption("w", OBConversion::OUTOPTIONS);
      MolConv.AddOption("X", OBConversion::OUTOPTIONS, "RecMet FixedH");
    }
    else
      MolConv.AddOption("w", OBConversion::OUTOPTIONS);
    bool success = pInChIFormat->WriteMolecule(mol, &MolConv);
    if (!success) return false;
    string inchi = newstream.str();
    if (inchi.size() == 0) return false;
    vector<string> vs;
    tokenize(vs, inchi);
    MolConv.SetInFormat(pInChIFormat);
    success = MolConv.ReadString(mol, vs.at(0));
    mol->DeleteData("inchi"); // Tidy up this side-effect
    return success;
  }

  static const char* NewLine()
  {
#ifdef _WIN32
    return "\r\n";
#else
    return "\n";
#endif
  }

  //////////////////////////////////////////////////
  bool SMIBaseFormat::WriteMolecule(OBBase* pOb,OBConversion* pConv)
  {
    //cout << "SMIBaseFromat::WriteMolecule()" << endl;
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);

    // Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();

    // Inchified SMILES? If so, then replace mol with the new 'normalised' one
    if (pConv->IsOption("I")) {
      bool success = GetInchifiedSMILESMolecule(pmol, false);
      if (!success) {
        ofs << NewLine();
        obErrorLog.ThrowError(__FUNCTION__, "Cannot generate Universal NSMILES for this molecule", obError);
        return false;
      }
    }

    // Title only option?
    if(pConv->IsOption("t")) {
      ofs << pmol->GetTitle() << NewLine();
      return true;
    }

    // Option 'x' needs "SMILES Atom Order" to be set
    // FIXME: When we support features of CXN extended SMILES
    //        we can rewrite this
    if (pConv->IsOption("x"))
      pConv->AddOption("O");

    std::string buffer;
    buffer.reserve(1000);

    // If there is data attached called "SMILES_Fragment", then it's
    // an ascii OBBitVec, representing the atoms of a fragment.  The
    // SMILES generated will only include these fragment atoms.

    OBBitVec fragatoms(pmol->NumAtoms());

    OBPairData *dp = (OBPairData *) pmol->GetData("SMILES_Fragment");
    const char* ppF = pConv->IsOption("F");
    if (dp) {
      fragatoms.FromString(dp->GetValue(), pmol->NumAtoms());
    }
    else if (ppF) { // Use info from option "F"
      fragatoms.FromString(ppF, pmol->NumAtoms());
    }
    // If no "SMILES_Fragment" data, fill the entire OBBitVec
    // with 1's so that the SMILES will be for the whole molecule.
    else {
      FOR_ATOMS_OF_MOL(a, *pmol)
        {
          fragatoms.SetBitOn(a->GetIdx());
        }
    }

    if (pmol->NumAtoms() > 0 || pmol->IsReaction()) {
      CreateCansmiString(*pmol, buffer, fragatoms, pConv);
    }

    bool writenewline = false;
    if(!pConv->IsOption("smilesonly")) {

      if(!pConv->IsOption("n")) {
        buffer += '\t';
        buffer += pmol->GetTitle();
      }

      if (pConv->IsOption("x") && pmol->HasData("SMILES Atom Order")) {
        vector<string> vs;
        string canorder = pmol->GetData("SMILES Atom Order")->GetValue();
        tokenize(vs, canorder);
        buffer += '\t';
        char tmp[15];
        for (unsigned int i = 0; i < vs.size(); i++) {
          int idx = atoi(vs[i].c_str());
          OBAtom *atom = pmol->GetAtom(idx);
          if (i > 0)
            buffer += ',';
          snprintf(tmp, 15, "%.4f", atom->GetX());
          buffer += tmp;
          buffer += ',';
          snprintf(tmp, 15, "%.4f", atom->GetY());
          buffer += tmp;
        }
      }

      if(!pConv->IsOption("nonewline"))
        writenewline = true;
    }

    ofs << buffer;
    if (writenewline)
      ofs << NewLine();

    return true;
  }

  //********************************************************
  class FIXFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    FIXFormat()
    {
      OBConversion::RegisterFormat("fix",this);
    }

    virtual const char* Description() //required
    {
      return
        "SMILES FIX format\n"
        "  No comments yet\n";
    };

    virtual const char* SpecificationURL()
    {return "";}; //optional

    //Flags() can return be any the following combined by | or be omitted if none apply
    // NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY
    virtual unsigned int Flags()
    {
      return NOTREADABLE;
    };

    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

  };

  //Make an instance of the format class
  FIXFormat theFIXFormat;

  /////////////////////////////////////////////////////////////////

  bool FIXFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    std::string buffer;
    OutOptions options(!pConv->IsOption("i"), pConv->IsOption("k"),
      pConv->IsOption("a"),
      pConv->IsOption("h"), pConv->IsOption("s"),
      pConv->IsOption("o"));
    OBMol2Cansmi m2s(options);

    m2s.Init(pmol, true, pConv);

    // We're outputting a full molecule
    // so we pass a bitvec for all atoms
    OBBitVec allbits(mol.NumAtoms());
    FOR_ATOMS_OF_MOL(a, mol) {
      allbits.SetBitOn(a->GetIdx());
    }

    if (mol.NumAtoms() > 0) {
      CreateCansmiString(mol, buffer, allbits, pConv);
    }
    ofs << buffer << endl;

    OBAtom *atom;
    vector<int>::iterator i;
    // Retrieve the canonical order of the molecule
    std::string orderString;
    m2s.GetOutputOrder(orderString);
    vector<string> canonical_order;
    tokenize(canonical_order, orderString);

    int j;
    int atomIdx;
    char coords[100];
    for (j = 0;j < mol.NumConformers();j++)
      {
        mol.SetConformer(j);
        for (unsigned int index = 0; index < canonical_order.size();
             ++index) {
          atomIdx = atoi(canonical_order[index].c_str());
          atom = mol.GetAtom(atomIdx);
          snprintf(coords, 100, "%9.3f %9.3f %9.3f", atom->GetX(), atom->GetY(), atom->GetZ());
          ofs << coords << endl;
        }
      }
    return(true);
  }

} // end namespace OpenBabel

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

#ifdef WIN32
#pragma warning (disable: 4786) // warning: long & complicated stl warning
#endif

#include "mol.h"
#include "typer.h"

using namespace std;

namespace OpenBabel {

  //FF extern OBAromaticTyper  aromtyper;
  //FF extern OBAtomTyper      atomtyper;

class OBSmilesParser
{
  int _bondflags;
  int _order;
  int _prev;
  char *_ptr;
  vector<int> _vprev;
  vector<vector<int> > _rclose;
  vector<vector<int> > _extbond;
  vector<int>          _path;
  vector<bool>         _avisit;
  vector<bool>         _bvisit;
  char _buffer[BUFF_SIZE];
public:
  bool SmiToMol(OBMol&,string&);
  bool ParseSmiles(OBMol&);
  bool ParseSimple(OBMol&);
  bool ParseComplex(OBMol&);
  bool ParseRingBond(OBMol&);
	bool ParseExternalBond(OBMol&);
	bool CapExternalBonds(OBMol &mol);
  void FindAromaticBonds(OBMol &mol,OBAtom*,int);
  void FindAromaticBonds(OBMol&);
  void FindOrphanAromaticAtoms(OBMol &mol); //CM 18 Sept 2003
};

bool OBSmilesParser::SmiToMol(OBMol &mol,string &s)
{
  strcpy(_buffer,s.c_str());

  _vprev.clear();
  _rclose.clear();
  _prev=0;

  if (!ParseSmiles(mol))
    {
      mol.Clear();
      return(false);
    }

  return(true);
}

bool OBSmilesParser::ParseSmiles(OBMol &mol)
{
  mol.BeginModify();

  for (_ptr=_buffer;*_ptr;_ptr++)
    {
      if (isspace(*_ptr))
	continue;
      else if (isdigit(*_ptr) || *_ptr == '%') //ring open/close
	  {
		  ParseRingBond(mol);
		  continue;
	  }
	  else if(*_ptr == '&') //external bond
	  {
		  ParseExternalBond(mol);
		  continue;
	  }
      else switch(*_ptr)
	  {
        case '.': _prev=0; break;
		case '(': 
			_vprev.push_back(_prev); 
			break;
		case ')': 
			_prev = _vprev.back(); 
			_vprev.pop_back(); 
			break;
		case '[': 
			if (!ParseComplex(mol)) 
			{
				mol.EndModify();mol.Clear();
				return(false); 
			}
			break; 
		case '-':  _order = 1; break;
		case '=':  _order = 2; break;
		case '#':  _order = 3; break;
		case ':':  _order = 5; break;
		case '/':  _bondflags |= OB_TORDOWN_BOND; break;
		case '\\': _bondflags |= OB_TORUP_BOND; break;
		default: 
			if (!ParseSimple(mol))
			{
				mol.EndModify(); mol.Clear();
				return(false);
			}
		} // end switch
		} // end for _ptr

	// place dummy atoms for each unfilled external bond
  if(!_extbond.empty())	
    CapExternalBonds(mol);

  //set aromatic bond orders
  mol.SetAromaticPerceived();
  FindAromaticBonds(mol);
  FindOrphanAromaticAtoms(mol);// CM 18 Sept 2003
  //FF spin multiplicity for H-deficient atoms no longer assigned
  //in AssignImplicitValence
  //atomtyper.AssignImplicitValence(mol); // CM and set _spinmultiplicities for H-deficient atoms
  mol.AssignSpinMultiplicity();
  //FF end
  mol.UnsetAromaticPerceived();
  mol.EndModify();

  return(true);
}

// CM 18 Sept 2003
void OBSmilesParser::FindOrphanAromaticAtoms(OBMol &mol)
{
  //Facilitates the use shorthand for radical entry: CcC = isopropyl radical
  //Atoms which are marked as aromatic but have no aromatic bonds
  //are taken to be radical centres 
  OBAtom *atom;
  vector<OBNodeBase*>::iterator j;
 
  for (atom = mol.BeginAtom(j);atom;atom = mol.NextAtom(j))
    if(atom->IsAromatic() && !atom->HasBondOfOrder(5)) //bonds order 5 set in FindAromaticBonds()
      atom->UnsetAromatic();
}
// CM end
 
void OBSmilesParser::FindAromaticBonds(OBMol &mol)
{
  _path.clear(); _avisit.clear(); _bvisit.clear();
  _avisit.resize(mol.NumAtoms()+1);
  _bvisit.resize(mol.NumBonds());
  _path.resize(mol.NumAtoms()+1);

  OBBond *bond;
  vector<OBEdgeBase*>::iterator i;
  for (bond = mol.BeginBond(i);bond;bond = mol.NextBond(i))
    if (!bond->GetBeginAtom()->IsAromatic() || 
	!bond->GetEndAtom()->IsAromatic())// ||
	//      bond->IsInRing())   // SM 18 Sept 2003
      _bvisit[bond->GetIdx()] = true;

  OBAtom *atom;
  vector<OBNodeBase*>::iterator j;

  for (atom = mol.BeginAtom(j);atom;atom = mol.NextAtom(j))
    if(!_avisit[atom->GetIdx()] && atom->IsAromatic())
      FindAromaticBonds(mol,atom,0);
}

void OBSmilesParser::FindAromaticBonds(OBMol &mol,OBAtom *atom,int depth )
{
  OBBond *bond;
  vector<OBEdgeBase*>::iterator k;

  if (_avisit[atom->GetIdx()])
    {
      int j = depth-1;
      bond=mol.GetBond(_path[j--]); 
      bond->SetBO(5);
      while( j >= 0 )
        { 
	  bond=mol.GetBond(_path[j--]);
	  bond->SetBO(5);
	  if(bond->GetBeginAtom() == atom || bond->GetEndAtom() == atom)
	    break;
        }
    } 
  else
    {
      _avisit[atom->GetIdx()] = true;
      for(bond = atom->BeginBond(k);bond;bond = atom->NextBond(k))
	if( !_bvisit[bond->GetIdx()])
	  {  
	    _path[depth] = bond->GetIdx();
	    _bvisit[bond->GetIdx()] = true;
	    FindAromaticBonds(mol,bond->GetNbrAtom(atom),depth+1);
	  }
    }
}


bool OBSmilesParser::ParseSimple(OBMol &mol)
{
  char symbol[3];
  int element;
  bool arom=false;
  memset(symbol,'\0',sizeof(char)*3);

  if (isupper(*_ptr))
    switch(*_ptr)
    {
    case 'C': 
      _ptr++;
      if (*_ptr == 'l') {strcpy(symbol,"Cl"); element = 17;}
      else              {symbol[0] = 'C';element = 6;_ptr--;}
      break;

    case 'N':  element = 7;  symbol[0] = 'N'; break;
    case 'O':  element = 8;  symbol[0] = 'O'; break;
    case 'S':  element = 16; symbol[0] = 'S'; break;
    case 'P':  element = 15; symbol[0] = 'P'; break;
    case 'F':  element = 9;  symbol[0] = 'F'; break;
    case 'I':  element = 53; symbol[0] = 'I'; break;

    case 'B':  
      _ptr++;
      if (*_ptr == 'r') {element = 35;strcpy(symbol,"Br");}
      else              {element = 5;symbol[0] = 'B';_ptr--;}
      break;
    default: return(false);
    }
  else
    {
      arom = true;
      switch(*_ptr)
	{
	case 'c': element = 6;  symbol[0] = 'C'; break;
	case 'n': element = 7;  symbol[0] = 'N'; break;
	case 'o': element = 8;  symbol[0] = 'O'; break;
	case 'p': element = 15; symbol[0] = 'P'; break;
	case 's': element = 16; symbol[0] = 'S'; break;
	case '*': element = 0;  strcpy(symbol,"Du"); break;
	default: return(false);
	}
    }

  OBAtom *atom = mol.NewAtom();
  atom->SetAtomicNum(element);
  atom->SetType(symbol);
  if (arom)
    {
      atom->SetAromatic();
      //      atom->SetSpinMultiplicity(2); // CM 18 Sept 2003
    }

  if (_prev) //need to add bond 
    {
      /* CM 18 Sept 2003
	 Extension so that lower case c can represent a radical centre
	 and cccc... a carbon chain bonded by conjugated double bonds.
	 Atoms c,n,o, etc initially added as radical centre
	 unless _prev is a radical centre when both are made a normal atoms
	 connected by a double bond. 
	 Since they are still marked as aromatic, FindAromaticBonds() will
	 replace the bonds by aromatic bonds if they are in a ring and remove the
	 aromatic tag from the atoms if they are not.
      */
      if(arom)
	{
	  OBAtom* prevatom = mol.GetAtom(_prev);
	  if (prevatom->GetSpinMultiplicity())
	    {
	      prevatom->SetSpinMultiplicity(0);
	      atom->SetSpinMultiplicity(0);
	      _order=2;
	    }
  
	}
      // CM end
      mol.AddBond(_prev,mol.NumAtoms(),_order,_bondflags);
    }  
  
  //set values
  _prev = mol.NumAtoms();
  _order = 1;
  _bondflags = 0;

  return(true);
}

bool OBSmilesParser::ParseComplex(OBMol &mol)
{
  char symbol[3];
  int element=0;
  int isotope=0;
  int isoPtr=0;
  bool arom=false;
  memset(symbol,'\0',sizeof(char)*3);

  _ptr++;

  //grab isotope information
  for (;*_ptr && isdigit(*_ptr);_ptr++)
    {
      symbol[isoPtr] = *_ptr;
      isoPtr++;
    }
  isotope = atoi(symbol);
  
  //parse element data
  if (isupper(*_ptr))
    switch(*_ptr)
      {
      case 'C':
	_ptr++;
	switch(*_ptr)
	  {
	  case 'a':  element = 20;  strcpy(symbol,"Ca"); break;
	  case 'd':  element = 48;  strcpy(symbol,"Cd"); break;
	  case 'e':  element = 58;  strcpy(symbol,"Ce"); break;
	  case 'f':  element = 98;  strcpy(symbol,"Cf"); break;
	  case 'l':  element = 17;  strcpy(symbol,"Cl"); break;
	  case 'm':  element = 96;  strcpy(symbol,"Cm"); break;
	  case 'o':  element = 27;  strcpy(symbol,"Co"); break;
	  case 'r':  element = 24;  strcpy(symbol,"Cr"); break;
	  case 's':  element = 55;  strcpy(symbol,"Cs"); break;
	  case 'u':  element = 29;  strcpy(symbol,"Cu"); break;
	  default:
	    element =  6;
	    symbol[0] = 'C';
	    _ptr--;
	  }
	break;

      case 'N':
	_ptr++;
	switch(*_ptr)
	  {
	  case 'a':  element =  11; strcpy(symbol,"Na"); break;
	  case 'b':  element =  41; strcpy(symbol,"Nb"); break;
	  case 'd':  element =  60; strcpy(symbol,"Nd"); break;
	  case 'e':  element =  10; strcpy(symbol,"Ne"); break;
	  case 'i':  element =  28; strcpy(symbol,"Ni"); break;
	  case 'o':  element = 102; strcpy(symbol,"No"); break;
	  case 'p':  element =  93; strcpy(symbol,"Np"); break;
	  default:    
	    element =   7;
	    symbol[0] = 'N';
	    _ptr--;
	  }
	break;

      case('O'):
	_ptr++;
	if(*_ptr == 's') {element = 76; strcpy(symbol,"Os");}
	else             {element = 8; symbol[0] = 'O';_ptr--;}
	break;

      case 'P':  
	_ptr++;
	switch(*_ptr)
	  {
	  case 'a':  element = 91; strcpy(symbol,"Pa"); break;
	  case 'b':  element = 82; strcpy(symbol,"Pb"); break;
	  case 'd':  element = 46; strcpy(symbol,"Pd"); break;
	  case 'm':  element = 61; strcpy(symbol,"Pm"); break;
	  case 'o':  element = 84; strcpy(symbol,"Po"); break;
	  case 'r':  element = 59; strcpy(symbol,"Pr"); break;
	  case 't':  element = 78; strcpy(symbol,"Pt"); break;
	  case 'u':  element = 94; strcpy(symbol,"Pu"); break;
	  default:    
	    element = 15;
	    symbol[0] = 'P';
	    _ptr--;
	  }
	break;

      case('S'):
	_ptr++;
	switch(*_ptr)
	  {
	  case 'b':  element = 51; strcpy(symbol,"Sb"); break;
	  case 'c':  element = 21; strcpy(symbol,"Sc"); break;
	  case 'e':  element = 34; strcpy(symbol,"Se"); break;
	  case 'i':  element = 14; strcpy(symbol,"Si"); break;
	  case 'm':  element = 62; strcpy(symbol,"Sm"); break;
	  case 'n':  element = 50; strcpy(symbol,"Sn"); break;
	  case 'r':  element = 38; strcpy(symbol,"Sr"); break;
	  default:    
	    element = 16; 
	    symbol[0] = 'S';
	    _ptr--;
	  }
	break;

      case 'B': 
	_ptr++;
	switch(*_ptr)
	  {
	  case 'a':  element = 56;  strcpy(symbol,"Ba"); break;
	  case 'e':  element =  4;  strcpy(symbol,"Be"); break;
	  case 'i':  element = 83;  strcpy(symbol,"Bi"); break;
	  case 'k':  element = 97;  strcpy(symbol,"Bk"); break;
	  case 'r':  element = 35;  strcpy(symbol,"Br"); break;
	  default:    
	    element = 5;   
	    symbol[0] = 'B'; 
	    _ptr--; 
	  }
	break;
	
      case 'F':
	_ptr++;
	switch(*_ptr)
	  {
	  case 'e': element = 26;  strcpy(symbol,"Fe"); break;
	  case 'm': element = 100; strcpy(symbol,"Fm"); break;
	  case 'r': element = 87;  strcpy(symbol,"Fr"); break;
	  default:
	    element = 9;
	    symbol[0] = 'F';
	    _ptr--;
	  }
	break;

      case 'I':  
	_ptr++;
	switch(*_ptr)
	  {
	  case 'n': element = 49; strcpy(symbol,"In"); break;
	  case 'r': element = 77; strcpy(symbol,"Ir"); break;
	  default:
	    element = 53;
	    symbol[0] = 'I';
	    _ptr--;
	  }
	break;

      case 'A':  
	_ptr++;
	switch(*_ptr)
	  {
	  case 'c':  element = 89;  strcpy(symbol,"Ac"); break;
	  case 'g':  element = 47;  strcpy(symbol,"Ag"); break;
	  case 'l':  element = 13;  strcpy(symbol,"Al"); break;
	  case 'm':  element = 95;  strcpy(symbol,"Am"); break;
	  case 'r':  element = 18;  strcpy(symbol,"Ar"); break;
	  case 's':  element = 33;  strcpy(symbol,"As"); break;
	  case 't':  element = 85;  strcpy(symbol,"At"); break;
	  case 'u':  element = 79;  strcpy(symbol,"Au"); break;
	  default: 
	    _ptr--; 
	    return(false);
	  }
	break;

      case 'D':
	_ptr++;
	if (*_ptr == 'y') {element = 66; strcpy(symbol,"Dy");}
	else              {_ptr--; return(false);}
	break;

      case 'E': 
	_ptr++;
	switch(*_ptr)
	  {
	  case 'r': element = 68; strcpy(symbol,"Er"); break;
	  case 's': element = 99; strcpy(symbol,"Es"); break;
	  case 'u': element = 63; strcpy(symbol,"Eu"); break;
	  default:
	    _ptr--;
	    return(false);
	  }
	break;

      case 'G':
	_ptr++;
	switch (*_ptr)
	  {
	  case 'a': element = 31; strcpy(symbol,"Ga"); break;
	  case 'd': element = 64; strcpy(symbol,"Gd"); break;
	  case 'e': element = 32; strcpy(symbol,"Ge"); break;
	  default:
	    _ptr--;
	    return(false);
	  }
	break;

      case 'H': 
	_ptr++;
	switch (*_ptr)
	  {
	  case 'e': element =  2; strcpy(symbol,"He"); break;
	  case 'f': element = 72; strcpy(symbol,"Hf"); break;
	  case 'g': element = 80; strcpy(symbol,"Hg"); break;
	  case 'o': element = 67; strcpy(symbol,"Ho"); break;
	  default:
	    element = 1;
	    symbol[0] = 'H';
	    _ptr--;
	  }
	break;

      case 'K':
	_ptr++;
	if(*_ptr == 'r') {element = 36; strcpy(symbol,"Kr");}
	else             {element = 19; symbol[0] = 'K'; _ptr--;}
	break;

      case 'L':
	_ptr++;
	switch(*_ptr)
	  {
	  case 'a': element =  57; strcpy(symbol,"La"); break;
	  case 'i': element =   3; strcpy(symbol,"Li"); break;
	  case 'r': element = 103; strcpy(symbol,"Lr"); break;
	  case 'u': element =  71; strcpy(symbol,"Lu"); break;
	  default:
	    _ptr--;
	    return(false);
	  }
	break;

      case 'M':
	_ptr++;
	switch(*_ptr)
	  {
	  case 'd': element = 101; strcpy(symbol,"Md"); break;
	  case 'g': element =  12; strcpy(symbol,"Mg"); break;
	  case 'n': element =  25; strcpy(symbol,"Mn"); break;
	  case 'o': element =  42; strcpy(symbol,"Mo"); break;
	  default:
	    _ptr--;
	    return(false);
	  }
	break;

      case 'R':
	_ptr++;
	switch(*_ptr)
	  {
	  case 'a': element = 88; strcpy(symbol,"Ra"); break;
	  case 'b': element = 37; strcpy(symbol,"Rb"); break;
	  case 'e': element = 75; strcpy(symbol,"Re"); break;
	  case 'h': element = 45; strcpy(symbol,"Rh"); break;
	  case 'n': element = 86; strcpy(symbol,"Rn"); break;
	  case 'u': element = 44; strcpy(symbol,"Ru"); break;
	  default:
	    _ptr--;
	    return(false);
	  }
	break;

      case 'T':
	_ptr++;
	switch(*_ptr)
	  {
	  case 'a': element = 73; strcpy(symbol,"Ta"); break;
	  case 'b': element = 65; strcpy(symbol,"Tb"); break;
	  case 'c': element = 43; strcpy(symbol,"Tc"); break;
	  case 'e': element = 52; strcpy(symbol,"Te"); break;
	  case 'h': element = 90; strcpy(symbol,"Th"); break;
	  case 'i': element = 22; strcpy(symbol,"Ti"); break;
	  case 'l': element = 81; strcpy(symbol,"Tl"); break;
	  case 'm': element = 69; strcpy(symbol,"Tm"); break;
	  default:
	    _ptr--;
	    return(false);
	  }
	break;

      case('U'):  element = 92;  symbol[0] = 'U'; break;
      case('V'):  element = 23;  symbol[0] = 'V'; break;
      case('W'):  element = 74;  symbol[0] = 'W'; break;

      case('X'):  
	_ptr++;
	if (*_ptr == 'e') {element = 54; strcpy(symbol,"Xe");}
	else              {_ptr--; return(false);}
	break;

      case('Y'):
	_ptr++;
	if (*_ptr == 'b') {element = 70; strcpy(symbol,"Yb");}
	else              {element = 39; symbol[0] = 'Y'; _ptr--;}
	break;

      case('Z'):  
	_ptr++;
	switch(*_ptr)
	  {
	  case 'n': element = 30; strcpy(symbol,"Zn"); break;
	  case 'r': element = 40; strcpy(symbol,"Zr"); break;
	  default:
	    _ptr--;
	    return(false);
	  }
	break;
      }
  else
    {
      arom = true;
      switch(*_ptr)
	{
	case 'c': element = 6;  symbol[0] = 'C'; break;
	case 'n': element = 7;  symbol[0] = 'N'; break;
	case 'o': element = 8;  symbol[0] = 'O'; break;
	case 'p': element = 15; symbol[0] = 'P'; break;
	case 's': 
	  _ptr++;
	  if (*_ptr == 'e')
	    {element = 34; strcpy(symbol,"Se");}
	  else
	    {element = 16; symbol[0] = 'S'; _ptr--;}
	  break;
	case 'a': 
	  _ptr++;
	  if (*_ptr == 's')
	    {element = 33; strcpy(symbol,"As");}
	  else
			return(false);
	  break;
	default: return(false);
	}
    }

  //handle hydrogen count, stereochemistry, and charge

  OBAtom *atom = mol.NewAtom();
  int hcount = 0;
  int charge=0;
  char tmpc[2]; tmpc[1] = '\0';
  for (_ptr++;*_ptr && *_ptr != ']';_ptr++)
	{
		switch(*_ptr)
      {
      case '@': 
	_ptr++;
	if (*_ptr == '@') atom->SetClockwiseStereo();
	else             {atom->SetAntiClockwiseStereo(); _ptr--;}
	break;
      case '-':
	_ptr++;
	if (isdigit(*_ptr))
	  {tmpc[0] = *_ptr; charge = -atoi(tmpc);}
	else
	  {charge--; _ptr--;}
	break;
      case '+': 
	_ptr++;
	if (isdigit(*_ptr))
	  {tmpc[0] = *_ptr; charge = atoi(tmpc);}
	else
	  {charge++; _ptr--;}
	break;

  case 'H': 
			_ptr++;
	if (isdigit(*_ptr)) {tmpc[0] = *_ptr; hcount = atoi(tmpc);}
	else                {hcount = 1; _ptr--;}
	break;

      default:
	return(false);
      }
	}

  if (charge) atom->SetFormalCharge(charge);
  atom->SetAtomicNum(element);
  atom->SetIsotope(isotope);
  atom->SetType(symbol);
  if (arom) atom->SetAromatic();

  if (_prev) //need to add bond 
    mol.AddBond(_prev,mol.NumAtoms(),_order,_bondflags);

  //set values
  _prev = mol.NumAtoms();
  _order = 1;
  _bondflags = 0;

  //now add hydrogens

  for (int i = 0;i < hcount;i++)
    {
      atom = mol.NewAtom();
      atom->SetAtomicNum(1);
      atom->SetType("H");
      mol.AddBond(_prev,mol.NumAtoms(),1);
    }

  return(true);
}

bool OBSmilesParser::CapExternalBonds(OBMol &mol)
{

	if(_extbond.empty())
		return(true);

	OBAtom *atom;
	vector<vector<int> >::iterator bond;
	
	for(bond = _extbond.begin();bond != _extbond.end();bond++)
	{
		// create new dummy atom
		atom = mol.NewAtom();
		atom->SetAtomicNum(0);
		atom->SetType("*");

		// bond dummy atom to mol via external bond
		mol.AddBond((*bond)[1],atom->GetIdx(),(*bond)[2],(*bond)[3]);
		OBBond *refbond = atom->GetBond(mol.GetAtom((*bond)[1]));

		//record external bond information
		OBExternalBondData *xbd;
		if(mol.HasData(obExternalBondData))
			xbd = (OBExternalBondData*)mol.GetData(obExternalBondData);
		else
		{
			xbd = new OBExternalBondData;
			mol.SetData(xbd);
		}
		xbd->SetData(atom,refbond,(*bond)[0]);

		/* old code written by AGS -- mts
		{
			externalbonds = (vector<pair<int,pair<OBAtom *,OBBond *> > > *)mol.GetData("extBonds");
		}
		else
		{
			externalbonds = new vector<pair<int,pair<OBAtom *,OBBond *> > >;
		}

		//save data <external bond count, bond index>
		externalbonds->push_back(pair<int,pair<OBAtom *,OBBond *> > ((*bond)[0], pair<OBAtom *,OBBond *> (atom,mol.GetBond((*bond)[1],atom->GetIdx()))));
		mol.SetData("extBonds",externalbonds);
		*/

		//this data gets cleaned up in mol.Clear.
	}
		
	return(true);
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
		case ';':
			_order = 5;
			_ptr++;
			break;
		case '/':	//chiral, but _order still == 1
			_bondflags |= OB_TORDOWN_BOND;
			_ptr++;
			break;
			_ptr++;
		case '\\': // chiral, but _order still == 1
			_bondflags |= OB_TORUP_BOND;
			_ptr++;
			break;
		default: // no bond indicator just leave order = 1
			break;
	}

	if (*_ptr == '%') // external bond indicator > 10
	{
		_ptr++;
		str[0] = *_ptr; _ptr++;
		str[1] = *_ptr; str[2] = '\0';
	}
	else // simple single digit external bond indicator
	{
		str[0] = *_ptr;
		str[1] = '\0';
	}
	digit = atoi(str);	// convert indicator to digit

	//check for dot disconnect closures
	vector<vector<int> >::iterator j;
	int bondFlags,bondOrder;
	for(j = _extbond.begin();j != _extbond.end();j++)
	{
		if((*j)[0] == digit)
		{
			bondFlags = (_bondflags > (*j)[3]) ? _bondflags : (*j)[3];
			bondOrder = (_order > (*j)[2]) ? _order : (*j)[2];
			mol.AddBond((*j)[1],_prev,bondOrder,bondFlags);
			_extbond.erase(j);
			_bondflags = 0;
			_order = 0;
			return(true);
		}
	}

	//since no closures save another ext bond
  vector<int> vtmp(4); 
  vtmp[0] = digit; vtmp[1] = _prev;
  vtmp[2] = _order; vtmp[3] = _bondflags;

  _extbond.push_back(vtmp);
  _order = 1;
  _bondflags = 0;

  return(true);

}

bool OBSmilesParser::ParseRingBond(OBMol &mol)
{
  int digit;
  char str[10];

  if (*_ptr == '%')
    {
      _ptr++;
      str[0] = *_ptr; _ptr++;
      str[1] = *_ptr; str[2] = '\0';
    }
  else
    {str[0] = *_ptr; str[1] = '\0';}
  digit = atoi(str);

  int bf,ord;
  vector<vector<int> >::iterator j; 
  for (j = _rclose.begin();j != _rclose.end();j++)
    if ((*j)[0] == digit)
      {
		bf = (_bondflags > (*j)[3]) ? _bondflags : (*j)[3];
		ord = (_order > (*j)[2]) ? _order : (*j)[2];
		mol.AddBond((*j)[1],_prev,ord,bf,(*j)[4]);
		_rclose.erase(j);
		_bondflags = 0; _order = 1;
		return(true);
      }

  vector<int> vtmp(5); 
  vtmp[0] = digit; vtmp[1] = _prev;
  vtmp[2] = _order; vtmp[3] = _bondflags;
  vtmp[4] = mol.GetAtom(_prev)->GetValence(); //store position to insert closure bond
  for (j = _rclose.begin();j != _rclose.end();j++) //correct for multiple closure bonds to a single atom
  if ((*j)[1] == _prev)
	vtmp[4]++;

  _rclose.push_back(vtmp);
  _order = 1;
  _bondflags = 0;

  return(true);
}

bool SmiToMol(OBMol &mol,string &smi,const char *title)
{
  OBSmilesParser sp;
  mol.SetTitle(title);

  if (!sp.SmiToMol(mol,smi))
    return(false);

  return(true);
}

bool ReadSmiles(istream &ifs,OBMol &mol,const char *title)
{
  char buffer[BUFF_SIZE];

  if (!ifs.getline(buffer,BUFF_SIZE)) return(false);
  vector<string> vs;
  tokenize(vs,buffer);

  // RWT 10/3/2000
  // 
  // added the following to allow spaces in compound names (titles).
  // Essentially everything after the first space on a SMILES file line
  // is treated as the name.
  // Also had to change the condition a few lines below from:
  // if (vs.size() == 2) ...  to
  // if (vs.size() >= 2)  

  if (vs.size()>2) 
  {
    for (unsigned int i=2;i<vs.size(); i++) 
    {
      vs[1]=vs[1]+" "+vs[i];
    }
  }

  if (!vs.empty())
    {
      if (vs.size() >= 2) SmiToMol(mol,vs[0],(char*)vs[1].c_str());
      if (vs.size() == 1) SmiToMol(mol,vs[0],"");
    }

  return(true);
}

}

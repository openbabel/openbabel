/**********************************************************************
atom.cpp - Handle OBAtom class.

Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2008 by Geoffrey R. Hutchison
Some portions Copyright (C) 2003 by Michael Banck

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

#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/mol.h>
#include <openbabel/obiter.h>
#include <openbabel/molchrg.h>
#include <openbabel/ring.h>
#include <openbabel/phmodel.h>
#include <openbabel/builder.h>
#include <openbabel/elements.h>
#include <openbabel/chains.h>
#include <openbabel/obutil.h>
#include <openbabel/residue.h>
#include <openbabel/chains.h>

#include <openbabel/math/matrix3x3.h>

#if !HAVE_STRNCASECMP
extern "C" int strncasecmp(const char *s1, const char *s2, size_t n);
#endif

using namespace std;


namespace OpenBabel
{
  EXTERN OBChainsParser chainsparser;
  /** \class OBAtom atom.h <openbabel/atom.h>
      \brief Atom class

      To understand the OBAtom class it is important to state a key
      decision on which the design was based. The OBAtom class not only
      holds data, but facilitates extraction of data perceived from both
      the atom and the molecule. This prevents the OBMol class from
      becoming overly large and complicated.

      A number of data extraction methods perform what is called
      <a href="http://en.wikipedia.org/wiki/Lazy_evaluation">`Lazy Evaluation,'</a>
      which is essentially on-the-fly evaluation.
      For example, when an atom is queried as to whether it is cyclic
      or what it's hybridization state is the information is perceived
      automatically. The perception of a particular trait is actually
      performed on the entire molecule the first time it is requested of
      an atom or bond, and stored for subsequent requests for the same
      trait of additional atoms or bonds.The OBAtom class is similar to
      OBMol and the whole of Open Babel in that data access and modification
      is done through Get and Set methods.

      The following code demonstrates how to print out the atom numbers,
      element numbers, and coordinates of a molecule:
      \code
      OBMol mol;

      FOR_ATOMS_OF_MOL(atom, mol)
      {
         cout << atom->GetIdx() << ` `;
         cout << atom->GetAtomicNum() << ` `;
         cout << atom->GetVector() << endl;
      }
      \endcode
      A number of the property member functions indicate that atoms
      have some knowledge of their covalently attached neighbor atoms.
      Bonding information is partly redundant within a molecule in
      that an OBMol has a complete list of bonds in a molecule, and
      an OBAtom has a list bonds of which it is a member. The following
      code demonstrates how an OBAtom uses its bond information to loop
      over atoms attached to itself:
      \code
      OBMol mol;
      OBAtom *atom;

      atom = mol.GetAtom(1);
      FOR_NBORS_OF_ATOM(nbr, atom)
      {
         cout << "atom #" << atom->GetIdx() << " is attached to atom #"
              << nbr->GetIdx() << endl;
      }
      \endcode

      should produce an output like
      \code
      atom #1 is attached to atom #2
      \endcode
  */

  extern THREAD_LOCAL OBAromaticTyper  aromtyper;
  extern THREAD_LOCAL OBAtomTyper      atomtyper;
  extern THREAD_LOCAL OBPhModel        phmodel;
  EXTERN OBTypeTable      ttab;
  
  //
  // OBAtom member functions
  //

  OBAtom::OBAtom()
  {
    _parent = (OBMol*)NULL;
    Clear();
  }

  OBAtom::~OBAtom()
  {
    if (_residue != NULL)
      {
        _residue->RemoveAtom(this);
      }
    /*
      if (!_vdata.empty())
      {
      vector<OBGenericData*>::iterator m;
      for (m = _vdata.begin();m != _vdata.end();++m)
      delete *m;
      _vdata.clear();
      }
    */
  }

  bool OBAtom::Clear()
  {
    _c = (double**)NULL;
    _cidx = 0;
    _flags=0;
    _idx = 0;
    _hyb = 0;
    _ele = (char)0;
    _isotope = 0;
    _spinmultiplicity=0; // CM 18 Sept 2003
    _imph = 0;
    _fcharge = 0;
    _type[0] = '\0';
    _pcharge = 0.0;
    _vbond.clear();
    _vbond.reserve(4);
    _residue = (OBResidue*)NULL;
    _id = NoId;

    return(OBBase::Clear());
  }

  OBAtom &OBAtom::operator=(OBAtom &src)
  //copy atom information
  //bond info is not copied here as ptrs may be invalid
  {
    if (this != &src) {
      _idx = src.GetIdx();
      Duplicate(&src);
    }
    return(*this);
  }

  void OBAtom::Duplicate(OBAtom *src)
  {
    if (!src)
      return;

    _hyb = src->_hyb;
    _ele = src->_ele;
    _imph = src->_imph;
    _isotope = src->_isotope;
    _fcharge = src->_fcharge;
    _spinmultiplicity = src->_spinmultiplicity;
    strncpy(_type,src->_type, sizeof(_type) - 1);
    _type[sizeof(_type) - 1] = '\0';
    _pcharge = src->_pcharge;
    _v = src->GetVector();
    _flags = src->_flags;
    _residue = (OBResidue*)NULL;
    _id = src->_id;

    _vdata.clear();
    //Copy all the OBGenericData, providing the new atom
    vector<OBGenericData*>::iterator itr;
    for(itr=src->BeginData();itr!=src->EndData();++itr)
      {
        OBGenericData* pCopiedData = (*itr)->Clone(this);
        SetData(pCopiedData);
      }
  }

  bool OBAtom::IsConnected(OBAtom *a1)
  {
    OBBondIterator i;
    OBBond *bond;

    for (bond = BeginBond(i);bond;bond = NextBond(i))
      if (bond->GetBeginAtom() == a1 || bond->GetEndAtom() == a1)
        return(true);

    return(false);
  }

  bool OBAtom::IsOneThree(OBAtom *a1)
  {
    OBAtom *atom1,*atom2;
    OBBond *bond1,*bond2;
    OBBondIterator i,j;
    atom1 = this;
    atom2 = a1;

    for (bond1 = atom1->BeginBond(i);bond1;bond1 = atom1->NextBond(i))
      for (bond2 = atom2->BeginBond(j);bond2;bond2 = atom2->NextBond(j))
        if (bond1->GetNbrAtom(atom1) == bond2->GetNbrAtom(atom2))
          return(true);

    return(false);
  }

  bool OBAtom::IsOneFour(OBAtom *a1)
  {
    OBAtom *atom1,*atom2;
    OBBond *bond1,*bond2;
    OBBondIterator i,j;
    atom1 = this;
    atom2 = a1;

    for (bond1 = atom1->BeginBond(i);bond1;bond1 = atom1->NextBond(i))
      for (bond2 = atom2->BeginBond(j);bond2;bond2 = atom2->NextBond(j))
        if ((bond1->GetNbrAtom(atom1))->IsConnected(bond2->GetNbrAtom(atom2)))
          return(true);

    return(false);
  }

  bool OBAtom::IsAxial()
  {
    double tor;
    OBAtom *a,*b,*c;
    OBBondIterator i,j,k;

    for (a = BeginNbrAtom(i);a;a = NextNbrAtom(i))
      if (a->GetHyb() == 3 && a->IsInRing() && !(*i)->IsInRing())
        for (b = a->BeginNbrAtom(j);b;b = a->NextNbrAtom(j))
          if (b != this && b->IsInRing() && b->GetHyb() == 3)
            for (c = b->BeginNbrAtom(k);c;c = b->NextNbrAtom(k))
              if (c != a && c->IsInRing())
                {
                  tor = fabs(((OBMol*)GetParent())->GetTorsion(this,a,b,c));
                  return(tor > 55.0 && tor < 75.0);
                }

    return(false);
  }

  /**     This can be sketched as follows
          \code
              '*'
                 \
                  a=b
          \endcode
          where a and b are the 'apha' and 'beta' atoms, respectively and '*'
          indicates the current atom.
    **/
  bool OBAtom::HasAlphaBetaUnsat(bool includePandS)
  {
    OBAtom *a1,*a2;
    OBBondIterator i,j;

    for (a1 = BeginNbrAtom(i);a1;a1 = NextNbrAtom(i))
      if (includePandS || (a1->GetAtomicNum() != OBElements::Phosphorus && a1->GetAtomicNum() != OBElements::Sulfur))
        for (a2 = a1->BeginNbrAtom(j);a2;a2 = a1->NextNbrAtom(j))
          if (a2 != this && ((*j)->GetBondOrder() == 2 || (*j)->GetBondOrder() == 3 || (*j)->GetBondOrder() == 5))
            return(true);

    return(false);
  }

  bool OBAtom::HasBondOfOrder(unsigned int order)
  {
    OBBond *bond;
    OBBondIterator i;
    for (bond = BeginBond(i);bond;bond = NextBond(i))
      if (bond->GetBondOrder() == order)
        return(true);

    return(false);
  }

  int OBAtom::CountBondsOfOrder(unsigned int order)
  {
    int count = 0;
    OBBond *bond;
    OBBondIterator i;
    for (bond = BeginBond(i);bond;bond = NextBond(i))
      if (bond->GetBondOrder() == order)
        count++;

    return(count);
  }

  int OBAtom::HighestBondOrder()
  {
    unsigned int highest = 0;
    OBBond *bond;
    OBBondIterator i;
    for(bond = BeginBond(i); bond; bond = NextBond(i))
      if(bond->GetBondOrder() > highest)
        highest = bond->GetBondOrder();

    return(highest);
  }

  bool OBAtom::HasNonSingleBond()
  {
    OBBond *bond;
    OBBondIterator i;
    for (bond = BeginBond(i);bond;bond = NextBond(i))
      if (bond->GetBondOrder() != 1)
        return(true);

    return(false);
  }

  bool OBAtom::IsPolarHydrogen()
  {
    if (GetAtomicNum() != OBElements::Hydrogen)
      return(false);

    OBAtom *atom;
    OBBond *bond;
    OBBondIterator i;
    for (bond = BeginBond(i);bond;bond = NextBond(i))
      {
        atom = bond->GetNbrAtom(this);
        if (atom->GetAtomicNum() == 7)
          return(true);
        if (atom->GetAtomicNum() == 8)
          return(true);
        if (atom->GetAtomicNum() == 15)
          return(true);
        if (atom->GetAtomicNum() == 16)
          return(true);
      }

    return(false);
  }

  bool OBAtom::IsNonPolarHydrogen()
  {
    if (GetAtomicNum() != OBElements::Hydrogen)
      return(false);

    OBAtom *atom;
    OBBond *bond;
    OBBondIterator i;
    for (bond = BeginBond(i);bond;bond = NextBond(i))
      {
        atom = bond->GetNbrAtom(this);
        if (atom->GetAtomicNum() == 6)
          return(true);
      }

    return(false);
  }

  vector3 &OBAtom::GetVector()
  {
    if (!_c)
      return(_v);

    _v.Set((*_c)[_cidx],(*_c)[_cidx+1],(*_c)[_cidx+2]);
    return(_v);
  }

  const vector3 &OBAtom::GetVector() const
  {
    if (!_c)
      return(_v);

    _v.Set((*_c)[_cidx],(*_c)[_cidx+1],(*_c)[_cidx+2]);
    return(_v);
  }

  void OBAtom::SetVector()
  {
    //    obAssert(_c);
    if (_c)
      _v.Set((*_c)[_cidx],(*_c)[_cidx+1],(*_c)[_cidx+2]);
  }

  void OBAtom::SetVector(const vector3 &v)
  {
    if (!_c)
      _v = v;
    else
      {
        (*_c)[_cidx  ] = v.x();
        (*_c)[_cidx+1] = v.y();
        (*_c)[_cidx+2] = v.z();
      }
  }

  void OBAtom::SetVector(const double v_x,const double v_y,const double v_z)
  {
    if (!_c)
      _v.Set(v_x,v_y,v_z);
    else
      {
        (*_c)[_cidx  ] = v_x;
        (*_c)[_cidx+1] = v_y;
        (*_c)[_cidx+2] = v_z;
      }
  }

  void OBAtom::SetType(const char *type)
  {
    strncpy(_type,type, sizeof(_type) - 1);
    _type[sizeof(_type) - 1] = '\0';
    if (_ele == 1 && type[0] == 'D')
      _isotope = 2;
  }

  void OBAtom::SetType(const string &type)
  {
    strncpy(_type,type.c_str(), sizeof(_type) - 1);
    _type[sizeof(_type) - 1] = '\0';
    if (_ele == 1 && type[0] == 'D')
      _isotope = 2;
  }

  void OBAtom::SetIsotope(unsigned int iso)
  {
    _isotope = iso;
  }

  OBResidue *OBAtom::GetResidue()
  {
    OBMol *mol = this->GetParent();
    if (!mol->HasChainsPerceived())
      chainsparser.PerceiveChains(*mol);

    return _residue;
  }

  double OBAtom::GetAtomicMass() const
  {
    if (_isotope == 0)
      return OBElements::GetMass(_ele);
    else
      return OBElements::GetExactMass(_ele, _isotope);
  }

  double OBAtom::GetExactMass() const
  {
    return OBElements::GetExactMass(_ele, _isotope);
  }

  char *OBAtom::GetType()
  {
    OBMol *mol = (OBMol*)GetParent();
    if (mol)
      if (!mol->HasAtomTypesPerceived())
        atomtyper.AssignTypes(*((OBMol*)GetParent()));

    if (strlen(_type) == 0) // Somehow we still don't have a type!
      {
        char num[6];
        string fromType = ttab.GetFromType(); // save previous types
        string toType = ttab.GetToType();

        ttab.SetFromType("ATN");
        ttab.SetToType("INT");
        snprintf(num, 6, "%d", GetAtomicNum());
        ttab.Translate(_type, num);

        ttab.SetFromType(fromType.c_str());
        ttab.SetToType(toType.c_str());
      }
    if (_ele == 1 && _isotope == 2)
      snprintf(_type, 6, "%s", "D");

    return(_type);
  }

  unsigned int OBAtom::GetHyb() const
  {
    //hybridization is assigned when atoms are typed
    OBMol *mol = (OBMol*)((OBAtom*)this)->GetParent();
    if (mol && !mol->HasHybridizationPerceived())
      atomtyper.AssignHyb(*mol);

    return(_hyb);
  }


  unsigned int OBAtom::GetHvyDegree() const
  {
    unsigned int count=0;

    OBBond *bond;
    OBBondIterator i;
    for (bond = ((OBAtom*)this)->BeginBond(i); bond; bond = ((OBAtom*)this)->NextBond(i))
      if (bond->GetNbrAtom((OBAtom*)this)->GetAtomicNum() != OBElements::Hydrogen)
        count++;

    return(count);
  }

  unsigned int OBAtom::GetHeteroDegree() const
  {
    unsigned int count=0;
    OBBond *bond;
    OBBondIterator i;
    for (bond = ((OBAtom*)this)->BeginBond(i);bond;bond = ((OBAtom*)this)->NextBond(i))
      if (bond->GetNbrAtom((OBAtom*)this)->IsHeteroatom())
        count++;

    return((unsigned int)count);
  }

  double OBAtom::GetPartialCharge()
  {
    if (!GetParent())
      return(_pcharge);
    if (!((OBMol*)GetParent())->AutomaticPartialCharge())
      return(_pcharge);

    if (!((OBMol*)GetParent())->HasPartialChargesPerceived())
      {
        //seed partial charges are set in the atom typing procedure
        OBAtom *atom;
        OBMol *mol = (OBMol*)GetParent();
        vector<OBAtom*>::iterator i;
        for (atom = mol->BeginAtom(i); atom; atom = mol->NextAtom(i))
          atom->SetPartialCharge(0.0);

        phmodel.AssignSeedPartialCharge(*((OBMol*)GetParent()));
        OBGastChrg gc;
        gc.AssignPartialCharges(*((OBMol*)GetParent()));
      }

    return(_pcharge);
  }

  //! Returns true if nitrogen is part of an amide
  bool OBAtom::IsAmideNitrogen()
  {
    if (GetAtomicNum() != OBElements::Nitrogen)
      return(false);

    OBAtom *nbratom,*atom;
    OBBond *abbond,*bond;

    OBBondIterator i,j;
    atom = this;
    for (bond = BeginBond(i);bond;bond = NextBond(i))
      {
        nbratom = bond->GetNbrAtom(atom);
        for (abbond = nbratom->BeginBond(j);abbond;abbond = nbratom->NextBond(j))
          if (abbond->GetBondOrder() == 2 &&
              (((abbond->GetNbrAtom(nbratom))->GetAtomicNum() == 8) ||
               ((abbond->GetNbrAtom(nbratom))->GetAtomicNum() == 16)))
            return(true);
      }

    return(false);
  }

  bool OBAtom::IsAromaticNOxide()
  {
    if (GetAtomicNum() != OBElements::Nitrogen || !IsAromatic())
      return(false);

    OBAtom *atom;
    OBBondIterator i;

    for (atom = BeginNbrAtom(i);atom;atom = NextNbrAtom(i))
      if (atom->GetAtomicNum() == OBElements::Oxygen && !(*i)->IsInRing() && (*i)->GetBondOrder() == 2)
        return(true);

    return(false);
  }

  bool OBAtom::IsCarboxylOxygen()
  {
    if (GetAtomicNum() != OBElements::Oxygen)
      return(false);
    if (GetHvyDegree() != 1)
      return(false);

    OBAtom *atom;
    OBBond *bond;
    OBBondIterator i;

    atom = NULL;
    for (bond = BeginBond(i);bond;bond = NextBond(i))
      if ((bond->GetNbrAtom(this))->GetAtomicNum() == OBElements::Carbon)
        {
          atom = bond->GetNbrAtom(this);
          break;
        }
    if (!atom)
      return(false);
    if (!(atom->CountFreeOxygens() == 2)
      && !(atom->CountFreeOxygens() == 1 && atom->CountFreeSulfurs() == 1))
      return(false);

    //atom is connected to a carbon that has a total
    //of 2 attached free oxygens or 1 free oxygen and 1 free sulfur
    return(true);
  }

  bool OBAtom::IsPhosphateOxygen()
  {
    if (GetAtomicNum() != OBElements::Oxygen)
      return(false);
    if (GetHvyDegree() != 1)
      return(false);
    OBAtom *atom;
    OBBond *bond;
    OBBondIterator i;

    atom = NULL;
    for (bond = BeginBond(i);bond;bond = NextBond(i))
      if ((bond->GetNbrAtom(this))->GetAtomicNum() == OBElements::Phosphorus)
        {
          atom = bond->GetNbrAtom(this);
          break;
        }
    if (!atom)
      return(false);
    if (atom->CountFreeOxygens() > 2)
      return(true);

    //atom is connected to a carbon that has a total
    //of 2 attached free oxygens
    return(false);
  }

  bool OBAtom::IsSulfateOxygen()
  {
    if (GetAtomicNum() != OBElements::Oxygen)
      return(false);
    if (GetHvyDegree() != 1)
      return(false);

    OBAtom *atom;
    OBBond *bond;
    OBBondIterator i;

    atom = NULL;
    for (bond = BeginBond(i);bond;bond = NextBond(i))
      if ((bond->GetNbrAtom(this))->GetAtomicNum() == OBElements::Sulfur)
        {
          atom = bond->GetNbrAtom(this);
          break;
        }
    if (!atom)
      return(false);
    if (atom->CountFreeOxygens() < 3)
      return(false);

    //atom is connected to a carbon that has a total
    //of 2 attached free oxygens
    return(true);
  }

  // Helper function for IsHBondAcceptor
  static bool IsSulfoneOxygen(OBAtom* atm)
  // Stefano Forli 
  //atom is connected to a sulfur that has a total
  //of 2 attached free oxygens, and it's not a sulfonamide
  //e.g. C-SO2-C
  // Is this atom an oxygen in a sulfone(R1 - SO2 - R2) group ?
  {
    if (atm->GetAtomicNum() != OBElements::Oxygen)
      return(false);
    if (atm->GetHvyDegree() != 1){
      //cerr << "sulfone> O valence is not 1\n";
      return(false);
      }

    OBAtom *nbr = NULL;
    OBBond *bond1,*bond2;
    OBBondIterator i,j;

    // searching for attached sulfur
    for (bond1 = atm->BeginBond(i); bond1; bond1 = atm->NextBond(i))
      if ((bond1->GetNbrAtom(atm))->GetAtomicNum() == OBElements::Sulfur)
        { nbr = bond1->GetNbrAtom(atm);
          break; }
    if (!nbr){
      //cerr << "sulfone> atom null\n" ;
      return(false); }

    // check for sulfate
    //cerr << "sulfone> If we're here... " << atom->GetAtomicNum() <<"\n" << atom->GetAtomicNum() == OBElements::Sulfur << "\n";
    //cerr << "sulfone> number of free oxygens:" << atom->CountFreeOxygens() << "\n";
    if (nbr->CountFreeOxygens() != 2){
      //cerr << "sulfone> count of free oxygens not 2" << atom->CountFreeOxygens() << '\n' ;
      return(false); }

    // check for sulfonamide
    for (bond2 = nbr->BeginBond(j);bond2;bond2 = nbr->NextBond(j)){
      //cerr<<"NEIGH: " << (bond2->GetNbrAtom(atom))->GetAtomicNum()<<"\n";
      if ((bond2->GetNbrAtom(nbr))->GetAtomicNum() == OBElements::Nitrogen){
        //cerr << "sulfone> sulfonamide null\n" ;
        return(false);}}
    //cerr << "sulfone> none of the above\n";
    return(true); // true sulfone
  }





  bool OBAtom::IsNitroOxygen()
  {
    if (GetAtomicNum() != OBElements::Oxygen)
      return(false);
    if (GetHvyDegree() != 1)
      return(false);

    OBAtom *atom;
    OBBond *bond;
    OBBondIterator i;

    atom = NULL;
    for (bond = BeginBond(i);bond;bond = NextBond(i))
      if ((bond->GetNbrAtom(this))->GetAtomicNum() == OBElements::Nitrogen)
        {
          atom = bond->GetNbrAtom(this);
          break;
        }
    if (!atom)
      return(false);
    if (atom->CountFreeOxygens() != 2)
      return(false);

    //atom is connected to a nitrogen that has a total
    //of 2 attached free oxygens
    return(true);
  }

  bool OBAtom::IsHeteroatom()
  {
    switch(GetAtomicNum())
      {
      case 7:
      case 8:
      case 15:
      case 16:
      case 33:
      case 34:
      case 51:
      case 52:
      case 83:
      case 84:
        return(true);
      }
    return(false);
  }

  bool OBAtom::IsAromatic() const
  {
    OBMol *mol = ((OBAtom*)this)->GetParent();
    if (!mol->HasAromaticPerceived())
      aromtyper.AssignAromaticFlags(*mol);

    if (((OBAtom*)this)->HasFlag(OB_AROMATIC_ATOM))
      return true;

    return false;
  }

  bool OBAtom::IsInRing() const
  {
    OBMol *mol = ((OBAtom*)this)->GetParent();
    if (!mol->HasRingAtomsAndBondsPerceived())
      mol->FindRingAtomsAndBonds();

    if (((OBAtom*)this)->HasFlag(OB_RING_ATOM))
      return true;

    return false;
  }

  //! @todo
  bool OBAtom::IsChiral()
  {
    OBMol *mol = (OBMol*) GetParent();
    OBStereoFacade stereoFacade(mol);
    return stereoFacade.HasTetrahedralStereo(_id);
  }

  bool OBAtom::IsInRingSize(int size) const
  {
    vector<OBRing*> rlist;
    vector<OBRing*>::iterator i;

    OBMol *mol = (OBMol*)((OBAtom*)this)->GetParent();
    if (!mol->HasSSSRPerceived())
      mol->FindSSSR();

    if (!((OBAtom*)this)->HasFlag(OB_RING_ATOM))
      return(false);

    rlist = mol->GetSSSR();
    for (i = rlist.begin();i != rlist.end();++i)
      if ((*i)->IsInRing(GetIdx()) && static_cast<int>((*i)->PathSize()) == size)
        return(true);

    return(false);
  }

  unsigned int OBAtom::MemberOfRingCount() const
  {
    vector<OBRing*> rlist;
    vector<OBRing*>::iterator i;
    unsigned int count=0;

    OBMol *mol = (OBMol*)((OBAtom*)this)->GetParent();

    if (!mol->HasSSSRPerceived())
      mol->FindSSSR();

    if (!((OBAtom*)this)->IsInRing())
      return(0);

    rlist = mol->GetSSSR();

    for (i = rlist.begin();i != rlist.end();++i)
      if ((*i)->IsInRing(GetIdx()))
        count++;

    return((unsigned int)count);
  }

  unsigned int OBAtom::MemberOfRingSize() const
  {
    vector<OBRing*> rlist;
    vector<OBRing*>::iterator i;

    OBMol *mol = (OBMol*)((OBAtom*)this)->GetParent();

    if (!mol->HasSSSRPerceived())
      mol->FindSSSR();

    if (!((OBAtom*)this)->IsInRing())
      return(0);

    rlist = mol->GetSSSR();

    for (i = rlist.begin();i != rlist.end();++i)
      if ((*i)->IsInRing(GetIdx()))
        return((*i)->Size());

    return(0);
  }

  unsigned int OBAtom::CountRingBonds() const
  {
    unsigned int count = 0;
    OBBond *bond;
    OBBondIterator i;

    for (bond = ((OBAtom*)this)->BeginBond(i);
         bond;
         bond = ((OBAtom*)this)->NextBond(i))
      {
        if (bond->IsInRing())
          count++;
      }
    return count;
  }

  double	  OBAtom::SmallestBondAngle()
  {
    OBAtom *b, *c;
    vector3 v1, v2;
    double degrees, minDegrees;
    //    vector<OBAtom*>::iterator i;
    OBBondIterator j,k;

    minDegrees = 360.0;

    for (b = BeginNbrAtom(j); b; b = NextNbrAtom(j))
      {
        k = j;
        for (c = NextNbrAtom(k); c; c = NextNbrAtom(k))
          {
            v1 = b->GetVector() - GetVector();
            v2 = c->GetVector() - GetVector();
            degrees = vectorAngle(v1, v2);
            if (degrees < minDegrees)
              minDegrees = degrees;
          }
      }
    return minDegrees;
  }

  double	  OBAtom::AverageBondAngle()
  {
    OBAtom *b, *c;
    vector3 v1, v2;
    double degrees, avgDegrees;
    //    vector<OBAtom*>::iterator i;
    OBBondIterator j,k;
    int n=0;

    avgDegrees = 0.0;

    for (b = BeginNbrAtom(j); b; b = NextNbrAtom(j))
      {
        k = j;
        for (c = NextNbrAtom(k); c; c = NextNbrAtom(k))
          {
            v1 = b->GetVector() - GetVector();
            v2 = c->GetVector() - GetVector();
            degrees = vectorAngle(v1, v2);
            avgDegrees += degrees;
            n++;
          }
      }

    if (n >=1)
      avgDegrees /= n;

    return avgDegrees;
  }

  unsigned int OBAtom::CountFreeOxygens() const
  {
    unsigned int count = 0;
    OBAtom *atom;
    OBBond *bond;
    OBBondIterator i;

    for (bond = ((OBAtom*)this)->BeginBond(i);bond;bond = ((OBAtom*)this)->NextBond(i))
      {
        atom = bond->GetNbrAtom((OBAtom*)this);
        if (atom->GetAtomicNum() == OBElements::Oxygen && atom->GetHvyDegree() == 1)
          count++;
      }

    return(count);
  }

  unsigned int OBAtom::CountFreeSulfurs() const
  {
    unsigned int count = 0;
    OBAtom *atom;
    OBBond *bond;
    OBBondIterator i;

    for (bond = ((OBAtom*)this)->BeginBond(i);bond;bond = ((OBAtom*)this)->NextBond(i))
      {
        atom = bond->GetNbrAtom((OBAtom*)this);
        if (atom->GetAtomicNum() == OBElements::Sulfur && atom->GetHvyDegree() == 1)
          count++;
      }

    return(count);
  }

  unsigned int OBAtom::GetExplicitValence() const
  {
    unsigned int bosum = 0;
    
    OBBondIterator i;
    for (OBBond *bond = ((OBAtom*)this)->BeginBond(i); bond; bond = ((OBAtom*)this)->NextBond(i))
      bosum += bond->GetBondOrder();

    return bosum;
  }

  unsigned int OBAtom::GetTotalValence() const
  {
    return GetExplicitValence() + _imph;
  }

  unsigned int OBAtom::ExplicitHydrogenCount(bool ExcludeIsotopes) const
  {
    //If ExcludeIsotopes is true, H atoms with _isotope!=0 are not included.
    //This excludes D, T and H when _isotope exlicitly set to 1 rather than the default 0.
    int numH=0;
    OBAtom *atom;
    OBBondIterator i;
    for (atom = ((OBAtom*)this)->BeginNbrAtom(i);atom;atom = ((OBAtom*)this)->NextNbrAtom(i))
      if (atom->GetAtomicNum() == OBElements::Hydrogen && !(ExcludeIsotopes && atom->GetIsotope()!=0))
        numH++;

    return(numH);
  }

  /**
   *  The returned values count whole lone pairs, so the acid count is the number of
   *  electron pairs desired and the base count is the number of electron pairs
   *  available.
   *
   @verbatim
   Algorithm from:
   Clark, A. M. Accurate Specification of Molecular Structures: The Case for
   Zero-Order Bonds and Explicit Hydrogen Counting. Journal of Chemical Information
   and Modeling, 51, 3149-3157 (2011). http://pubs.acs.org/doi/abs/10.1021/ci200488k
   @endverbatim
  */
  pair<int, int> OBAtom::LewisAcidBaseCounts() const
  {
    // TODO: Is this data stored elsewhere?
    // The number of valence electrons in a free atom
    const int VALENCE[113] = {0,1,2,1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8,9,10,
                              11,12,3,4,5,6,7,8,1,2,3,4,5,6,7,8,9,10,11,12,3,4,5,6,7,8,1,2,
                              4,4,4,4,4,4,4,4,4,4,4,4,4,4,3,4,5,6,7,8,9,10,11,12,3,4,5,6,7,
                              8,1,1,4,4,4,4,4,4,4,4,4,4,4,4,4,4,3,4,5,6,7,8,9,10,11,12};
    // The number of electrons required to make up a full valence shell
    const int SHELL[113]   = {0,2,2,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,18,18,18,18,18,18,
                              18,18,18,18,8,8,8,8,8,8,8,8,18,18,18,18,18,18,18,18,18,18,8,
                              8,8,8,8,8,8,8,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,
                              18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,8,8,18,18,18,18,
                              18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18};

    pair<int, int> counts;
    int N = GetAtomicNum();
    if (N == 0 || N > 112) {
      counts.first = 0;
      counts.second = 0;
    } else {
      int S = SHELL[N];
      int V = VALENCE[N];
      int C = GetFormalCharge();
      int B = GetTotalValence();
      // TODO: Do we actually want to divide by 2 here? (counting pairs instead of single)
      counts.first = (S - V - B + C) / 2;  // Acid: Number of electrons pairs desired
      counts.second = (V - B - C) / 2;     // Base: Number of electrons pairs available
    }
    return counts;
  }

  bool OBAtom::DeleteBond(OBBond *bond)
  {
    OBBondIterator i;
    for (i = _vbond.begin();i != _vbond.end();++i)
      if ((OBBond*)bond == *i)
        {
          _vbond.erase(i);
          return(true);
        }
    return(false);
  }

  bool OBAtom::MatchesSMARTS(const char *pattern)
  {
    OBMol *mol = (OBMol*)((OBAtom*)this)->GetParent();
    vector<vector<int> > mlist;
    vector<vector<int> >::iterator l;

    OBSmartsPattern test;
    test.Init(pattern);
    if (test.Match(*mol))
      {
        mlist = test.GetUMapList();
        for (l = mlist.begin(); l != mlist.end(); ++l)
          if (GetIdx() == mol->GetAtom((*l)[0])->GetIdx())
            return true;
      }
    return false;
  }

  OBBond *OBAtom::BeginBond(OBBondIterator &i)
  {
    i = _vbond.begin();
    return((i == _vbond.end()) ? (OBBond*)NULL : (OBBond*)*i);
  }

  OBBond *OBAtom::NextBond(OBBondIterator &i)
  {
    i++;
    return((i == _vbond.end()) ? (OBBond*)NULL : (OBBond*)*i);
  }

  OBAtom *OBAtom::BeginNbrAtom(OBBondIterator &i)
  {
    i = _vbond.begin();
    return((i != _vbond.end()) ? ((OBBond*) *i)->GetNbrAtom(this):NULL);
  }

  OBAtom *OBAtom::NextNbrAtom(OBBondIterator &i)
  {
    i++;
    return((i != _vbond.end()) ?  ((OBBond*) *i)->GetNbrAtom(this):NULL);
  }

  double OBAtom::GetDistance(OBAtom *b)
  {
    return(( this->GetVector() - b->GetVector() ).length());
  }

  double OBAtom::GetDistance(int b)
  {
    OBMol *mol = (OBMol*)GetParent();
    return(( this->GetVector() - mol->GetAtom(b)->GetVector() ).length());
  }

  double OBAtom::GetDistance(vector3 *v)
  {
    return(( this->GetVector() - *v ).length());
  }

  double OBAtom::GetAngle(OBAtom *b, OBAtom *c)
  {
    vector3 v1,v2;

    v1 = this->GetVector() - b->GetVector();
    v2 = c->GetVector() - b->GetVector();
    if (IsNearZero(v1.length(), 1.0e-3)
      || IsNearZero(v2.length(), 1.0e-3)) {
        return(0.0);
    }

    return(vectorAngle(v1, v2));
  }

  double OBAtom::GetAngle(int b, int c)
  {
    OBMol *mol = (OBMol*)GetParent();
    vector3 v1,v2;

    v1 = this->GetVector() - mol->GetAtom(b)->GetVector();
    v2 = mol->GetAtom(c)->GetVector() - mol->GetAtom(b)->GetVector();

    if (IsNearZero(v1.length(), 1.0e-3)
      || IsNearZero(v2.length(), 1.0e-3)) {
        return(0.0);
    }

    return(vectorAngle(v1, v2));
  }

  bool OBAtom::GetNewBondVector(vector3 &v,double length)
  {
    v = OBBuilder::GetNewBondVector(this, length);
    return(true);
  }

  //! \return a "corrected" bonding radius based on the hybridization.
  //! Scales the covalent radius by 0.95 for sp2 and 0.90 for sp hybrids
  static double CorrectedBondRad(unsigned int elem, unsigned int hyb)
  {
    double rad = OBElements::GetCovalentRad(elem);
    switch (hyb) {
    case 2:
      return rad * 0.95;
    case 1:
      return rad * 0.90;
    default:
      return rad;
    }
  }

  bool OBAtom::HtoMethyl()
  {
    if (GetAtomicNum() != OBElements::Hydrogen)
      return(false);

    obErrorLog.ThrowError(__FUNCTION__,
                          "Ran OpenBabel::HtoMethyl", obAuditMsg);

    OBMol *mol = (OBMol*)GetParent();

    mol->BeginModify();

    SetAtomicNum(6);
    SetType("C3");
    SetHyb(3);

    OBAtom *atom;
    OBBond *bond;
    OBBondIterator i;
    atom  = BeginNbrAtom(i);
    bond = (OBBond *)*i;
    if (!atom)
      {
        mol->EndModify();
        return(false);
      }

    double br1,br2;
    br1 = CorrectedBondRad(6,3);
    br2 = CorrectedBondRad(atom->GetAtomicNum(),atom->GetHyb());
    bond->SetLength(atom,br1+br2);

    OBAtom *hatom;
    br2 = CorrectedBondRad(1,0);
    vector3 v;

    for (int j = 0;j < 3;++j)
      {
        hatom = mol->NewAtom();
        hatom->SetAtomicNum(1);
        hatom->SetType("H");

        GetNewBondVector(v,br1+br2);
        hatom->SetVector(v);
        mol->AddBond(GetIdx(),mol->NumAtoms(),1);
      }

    mol->EndModify();
    return(true);
  }

  static void ApplyRotMatToBond(OBMol &mol,matrix3x3 &m,OBAtom *a1,OBAtom *a2)
  {
    vector<int> children;
    mol.FindChildren(children,a1->GetIdx(),a2->GetIdx());
    children.push_back(a2->GetIdx());

    vector3 v;
    vector<int>::iterator i;
    for (i = children.begin();i != children.end();++i)
      {
        v = mol.GetAtom(*i)->GetVector();
        v -= a1->GetVector();
        v *= m;
        v += a1->GetVector();
        mol.GetAtom(*i)->SetVector(v);
      }

  }

  //! \deprecated This will be removed in future versions of Open Babel
  bool OBAtom::SetHybAndGeom(int hyb)
  {
    obErrorLog.ThrowError(__FUNCTION__,
                          "Ran OpenBabel::SetHybridizationAndGeometry",
                          obAuditMsg);

    //if (hyb == GetHyb()) return(true);
    if (GetAtomicNum() == 1)
      return(false);
    if (hyb == 0 && GetHvyDegree() > 1)
      return(false);
    if (hyb == 1 && GetHvyDegree() > 2)
      return(false);
    if (hyb == 2 && GetHvyDegree() > 3)
      return(false);
    if (hyb == 3 && GetHvyDegree() > 4)
      return(false);

    OBMol *mol = (OBMol*)GetParent();

    OBAtom *nbr;
    vector<OBAtom*> delatm;
    OBBondIterator i;

    for (nbr = BeginNbrAtom(i);nbr;nbr = NextNbrAtom(i))
      if (nbr->GetAtomicNum() == OBElements::Hydrogen)
        delatm.push_back(nbr);

    //delete attached hydrogens
    mol->IncrementMod();
    vector<OBAtom*>::iterator j;
    for (j = delatm.begin();j != delatm.end();++j)
      mol->DeleteAtom(*j);
    mol->DecrementMod();

    double targetAngle;
    if (hyb == 3)
      targetAngle = 109.5;
    else if (hyb == 2)
      targetAngle = 120.0;
    else if (hyb == 1)
      targetAngle = 180.0;
    else
      targetAngle = 0.0;

    if (IsInRing())
      targetAngle = 180.0 - (360.0 / MemberOfRingSize());

    //adjust attached acyclic bond lengths
    double br1,br2, length;
    br1 = CorrectedBondRad(GetAtomicNum(),hyb);
    for (nbr = BeginNbrAtom(i);nbr;nbr = NextNbrAtom(i))
      if (!(*i)->IsInRing())
        {
          br2 = CorrectedBondRad(nbr->GetAtomicNum(),nbr->GetHyb());
          length = br1 + br2;
          if ((*i)->IsAromatic())
            length *= 0.93;
          else if ((*i)->GetBondOrder() == 2)
            length *= 0.91;
          else if ((*i)->GetBondOrder() == 3)
            length *= 0.87;
          ((OBBond*) *i)->SetLength(this, length);
        }

    if (GetExplicitDegree() > 1)
      {
        double angle;
        matrix3x3 m;
        vector3 v1,v2,v3,v4,n,s;
        OBAtom *r1,*r2,*r3,*a1,*a2,*a3,*a4;
        r1 = r2 = r3 = a1 = a2 = a3 = a4 = NULL;

        //find ring atoms first
        for (nbr = BeginNbrAtom(i);nbr;nbr = NextNbrAtom(i))
          if ((*i)->IsInRing()) {
            if (!r1) {
              r1 = nbr;
            }
            else if (!r2) {
              r2 = nbr;
            }
            else if (!r3) {
              r3 = nbr;
            }
          }

        //find non-ring atoms
        for (nbr = BeginNbrAtom(i);nbr;nbr = NextNbrAtom(i))
          if (!(*i)->IsInRing()) {
            if (!a1) {
              a1 = nbr;
            }
            else if (!a2) {
              a2 = nbr;
            }
            else if (!a3) {
              a3 = nbr;
            }
            else if (!a4) {
              a4 = nbr;
            }
          }

        //adjust geometries of heavy atoms according to hybridization
        if (hyb == 1)
          {
            if (a2)
              {
                v1 = a1->GetVector()-GetVector();
                v1.normalize();
                v2 = a2->GetVector()-GetVector();
                v2.normalize();
                n = cross(v1,v2);
                angle = vectorAngle(v1,v2)-targetAngle;
                m.RotAboutAxisByAngle(n,-angle);
                ApplyRotMatToBond(*mol,m,this,a1);
              }
          }
        else if (hyb == 2)
          {
            if (r1 && r2 && a1)
              {
                v1 = r1->GetVector()-GetVector();
                v1.normalize();
                v2 = r2->GetVector()-GetVector();
                v2.normalize();
                v3 = a1->GetVector()-GetVector();
                s = v1+v2;
                s.normalize();
                s *= -1.0;
                n = cross(s,v3);
                angle = vectorAngle(s,v3);
                m.RotAboutAxisByAngle(n,angle);
                ApplyRotMatToBond(*mol,m,this,a1);
              }
            else
              {
                if (a2)
                  {
                    v1 = a1->GetVector()-GetVector();
                    v1.normalize();
                    v2 = a2->GetVector()-GetVector();
                    v2.normalize();
                    n = cross(v1,v2);
                    angle = vectorAngle(v1,v2)-targetAngle;
                    m.RotAboutAxisByAngle(n,-angle);
                    ApplyRotMatToBond(*mol,m,this,a1);
                  }
                if (a3)
                  {
                    v1 = a1->GetVector()-GetVector();
                    v1.normalize();
                    v2 = a2->GetVector()-GetVector();
                    v2.normalize();
                    v3 = a3->GetVector()-GetVector();
                    s = v1+v2;
                    s.normalize();
                    s *= -1.0;
                    n = cross(s,v3);
                    angle = vectorAngle(s,v3);
                    m.RotAboutAxisByAngle(n,angle);
                    ApplyRotMatToBond(*mol,m,this,a3);
                  }
              }
          }
        else if (hyb == 3)
          {
            if (r1 && r2 && r3 && a1)
              {
                v1 = r1->GetVector()-GetVector();
                v1.normalize();
                v2 = r2->GetVector()-GetVector();
                v2.normalize();
                v3 = r3->GetVector()-GetVector();
                v3.normalize();
                v4 = a1->GetVector()-GetVector();
                s = v1 + v2 + v3;
                s *= -1.0;
                s.normalize();
                n = cross(s,v4);
                angle = vectorAngle(s,v4);
                m.RotAboutAxisByAngle(n,angle);
                ApplyRotMatToBond(*mol,m,this,a1);
              }
            else if (r1 && r2 && !r3 && a1)
              {
                v1 = r1->GetVector()-GetVector();
                v1.normalize();
                v2 = r2->GetVector()-GetVector();
                v2.normalize();
                v3 = a1->GetVector()-GetVector();
                s = v1+v2;
                s *= -1.0;
                s.normalize();
                n = cross(v1,v2);
                n.normalize();
#ifndef ONE_OVER_SQRT3
#define ONE_OVER_SQRT3  0.57735026918962576451
#endif //SQRT_TWO_THIRDS
#ifndef SQRT_TWO_THIRDS
#define SQRT_TWO_THIRDS 0.81649658092772603272
#endif //ONE_OVER_SQRT3
                s *= ONE_OVER_SQRT3; //1/sqrt(3)
                n *= SQRT_TWO_THIRDS; //sqrt(2/3)
                s += n;
                s.normalize();
                n = cross(s,v3);
                angle = vectorAngle(s,v3);
                m.RotAboutAxisByAngle(n,angle);
                ApplyRotMatToBond(*mol,m,this,a1);

                if (a2)
                  {
                    v1 = r1->GetVector()-GetVector();
                    v1.normalize();
                    v2 = r2->GetVector()-GetVector();
                    v2.normalize();
                    v3 = a1->GetVector()-GetVector();
                    v3.normalize();
                    v4 = a2->GetVector()-GetVector();
                    s = v1 + v2 + v3;
                    s *= -1.0;
                    s.normalize();
                    n = cross(s,v4);
                    angle = vectorAngle(s,v4);
                    m.RotAboutAxisByAngle(n,angle);
                    ApplyRotMatToBond(*mol,m,this,a2);
                  }
              }
            else
              {
                if (a2)
                  {
                    v1 = a1->GetVector()-GetVector();
                    v1.normalize();
                    v2 = a2->GetVector()-GetVector();
                    v2.normalize();
                    n = cross(v1,v2);
                    angle = vectorAngle(v1,v2)-targetAngle;
                    m.RotAboutAxisByAngle(n,-angle);
                    ApplyRotMatToBond(*mol,m,this,a1);
                  }
                if (a3)
                  {
                    v1 = a1->GetVector()-GetVector();
                    v1.normalize();
                    v2 = a2->GetVector()-GetVector();
                    v2.normalize();
                    v3 = a3->GetVector()-GetVector();
                    s = v1+v2;
                    s *= -1.0;
                    s.normalize();
                    n = cross(v1,v2);
                    n.normalize();
                    s *= ONE_OVER_SQRT3; //1/sqrt(3)
                    n *= SQRT_TWO_THIRDS; //sqrt(2/3)
                    s += n;
                    s.normalize();
                    n = cross(s,v3);
                    angle = vectorAngle(s,v3);
                    m.RotAboutAxisByAngle(n,angle);
                    ApplyRotMatToBond(*mol,m,this,a3);
                  }
              }
          }
      }

    //add hydrogens back to atom
    int impval=1;
    switch (GetAtomicNum())
      {
      case 6:
        if (hyb == 3)
          impval = 4;
        if (hyb == 2)
          impval = 3;
        if (hyb == 1)
          impval = 2;
        break;
      case 7:
        if (hyb == 3)
          impval = 3;
        if (hyb == 2)
          impval = 2;
        if (hyb == 1)
          impval = 1;
        break;
      case 8:
        if (hyb == 3)
          impval = 2;
        if (hyb == 2)
          impval = 2;
        if (hyb == 1)
          impval = 2;
        break;
      case 16:
        if (hyb == 3)
          impval = 2;
        if (hyb == 2)
          impval = 2;
        if (hyb == 1)
          impval = 0;
        break;
      case 15:
        if (hyb == 3)
          impval = 4;
        if (hyb == 2)
          impval = 3;
        if (hyb == 1)
          impval = 2;
        break;
      default:
        impval = 1;
      }

    int hcount = impval-GetHvyDegree();
    if (hcount)
      {
        int k;
        vector3 v;
        OBAtom *atom;
        double brsum = CorrectedBondRad(1,0)+
          CorrectedBondRad(GetAtomicNum(),GetHyb());
        SetHyb(hyb);

        mol->BeginModify();
        for (k = 0;k < hcount;++k)
          {
            GetNewBondVector(v,brsum);
            atom = mol->NewAtom();
            atom->SetAtomicNum(1);
            atom->SetType("H");
            atom->SetVector(v);
            mol->AddBond(atom->GetIdx(),GetIdx(),1);
          }
        mol->EndModify();
      }


    return(true);
  }

  OBBond *OBAtom::GetBond(OBAtom *nbr)
  {
    OBBond *bond;
    vector<OBBond *>::iterator i;
    for (bond = BeginBond(i) ; bond ; bond = NextBond(i))
      if (bond->GetNbrAtom(this) == nbr)
        return bond;
    return NULL;
  }

  /*Now in OBBase
  // OBGenericData methods
  bool OBAtom::HasData(string &s)
  //returns true if the generic attribute/value pair exists
  {
  if (_vdata.empty())
  return(false);

  vector<OBGenericData*>::iterator i;

  for (i = _vdata.begin();i != _vdata.end();++i)
  if ((*i)->GetAttribute() == s)
  return(true);

  return(false);
  }

  bool OBAtom::HasData(const char *s)
  //returns true if the generic attribute/value pair exists
  {
  if (_vdata.empty())
  return(false);

  vector<OBGenericData*>::iterator i;

  for (i = _vdata.begin();i != _vdata.end();++i)
  if ((*i)->GetAttribute() == s)
  return(true);

  return(false);
  }

  bool OBAtom::HasData(unsigned int dt)
  //returns true if the generic attribute/value pair exists
  {
  if (_vdata.empty())
  return(false);

  vector<OBGenericData*>::iterator i;

  for (i = _vdata.begin();i != _vdata.end();++i)
  if ((*i)->GetDataType() == dt)
  return(true);

  return(false);
  }

  OBGenericData *OBAtom::GetData(string &s)
  //returns the value given an attribute
  {
  vector<OBGenericData*>::iterator i;

  for (i = _vdata.begin();i != _vdata.end();++i)
  if ((*i)->GetAttribute() == s)
  return(*i);

  return(NULL);
  }

  OBGenericData *OBAtom::GetData(const char *s)
  //returns the value given an attribute
  {
  vector<OBGenericData*>::iterator i;

  for (i = _vdata.begin();i != _vdata.end();++i)
  if ((*i)->GetAttribute() == s)
  return(*i);

  return(NULL);
  }

  OBGenericData *OBAtom::GetData(unsigned int dt)
  {
  vector<OBGenericData*>::iterator i;
  for (i = _vdata.begin();i != _vdata.end();++i)
  if ((*i)->GetDataType() == dt)
  return(*i);
  return(NULL);
  }

  void OBAtom::DeleteData(unsigned int dt)
  {
  vector<OBGenericData*> vdata;
  vector<OBGenericData*>::iterator i;
  for (i = _vdata.begin();i != _vdata.end();++i)
  if ((*i)->GetDataType() == dt)
  delete *i;
  else
  vdata.push_back(*i);
  _vdata = vdata;
  }

  void OBAtom::DeleteData(vector<OBGenericData*> &vg)
  {
  vector<OBGenericData*> vdata;
  vector<OBGenericData*>::iterator i,j;

  bool del;
  for (i = _vdata.begin();i != _vdata.end();++i)
  {
  del = false;
  for (j = vg.begin();j != vg.end();++j)
  if (*i == *j)
  {
  del = true;
  break;
  }
  if (del)
  delete *i;
  else
  vdata.push_back(*i);
  }
  _vdata = vdata;
  }

  void OBAtom::DeleteData(OBGenericData *gd)
  {
  vector<OBGenericData*>::iterator i;
  for (i = _vdata.begin();i != _vdata.end();++i)
  if (*i == gd)
  {
  delete *i;
  _vdata.erase(i);
  }

  }
  */

  bool OBAtom::IsHbondAcceptorSimple()
  {
    // Changes from Liu Zhiguo
    if (_ele == 8 || _ele == 9)
      return true;
    if (_ele == 7) {
      // N+ ions and sp2 hybrid N with 3 valences should not be Hbond acceptors
      if (!((GetExplicitDegree() == 4 && GetHyb() == 3)
            || (GetExplicitDegree() == 3 && GetHyb() == 2)))
            return true;
    }
    // Changes from Paolo Tosco
    if (_ele == 16 && GetFormalCharge() == -1)
      return true;
    return false;
  }

  // new function, Stefano Forli
  // Incorporate ideas and data from Kubyni and others. 
  // [1] Kubinyi, H. "Changing paradigms in drug discovery.
  //    "SPECIAL PUBLICATION-ROYAL SOCIETY OF CHEMISTRY 304.1 (2006): 219-232.
  //
  // [2] Kingsbury, Charles A. "Why are the Nitro and Sulfone 
  //     Groups Poor Hydrogen Bonders?." (2015).
  //
  // [3] Per Restorp, Orion B. Berryman, Aaron C. Sather, Dariush Ajami 
  //     and Julius Rebek Jr., Chem. Commun., 2009, 5692 DOI: 10.1039/b914171e
  //
  // [4] Dunitz, Taylor. "Organic fluorine hardly ever accepts
  //     hydrogen bonds." Chemistry-A European Journal 3.1 (1997): 83-92.
  //
  // This function has a finer grain than the original
  // implementation, checking also the neighbor atoms. 
  // Accordingly to these rules, the function will return:
  //
  //    aliph-O-aliph ether   -> true   [1]
  //    hydroxy O-sp3         -> true   [1]
  //    aro-O-aliph ether     -> true   [1]
  //    ester O-sp2           -> true   [1]
  //    sulfate O (R-SO3)     -> true   [2]
  //    sulfoxyde O (R-SO-R)  -> true   [2]
  //    organoboron-F (R-BF3) -> true   [3]
  //    ester O-sp3           -> false  [1]
  //    sulfone (R1-SO2-R2 )  -> false  [2]
  //    aro-O-aro             -> false  [1]
  //    aromatic O            -> false  [1]
  //    O-nitro               -> false  [2]
  //    organic F (R-F)       -> false  [4]
  //    
  bool OBAtom::IsHbondAcceptor() {
      if (_ele == 8) {
        // oxygen; this should likely be a separate function
        // something like IsHbondAcceptorOxygen()
        unsigned int aroCount = 0;

        OBBond *bond;
        OBBondIterator i;
        if (IsNitroOxygen()){ // maybe could be a bool option in the function?
          return (false);
          }
        if (IsAromatic()){ // aromatic oxygen (furan) (NO)
          return(false);
          }
        if (IsSulfoneOxygen(this)){ // sulfone (NO)
          return(false);
          }
        FOR_NBORS_OF_ATOM(nbr, this){
          if (nbr->IsAromatic()){ 
            aroCount += 1;
            if (aroCount == 2){ // aromatic ether (aro-O-aro) (NO)
              return(false); 
            }
          }
          else { 
            if (nbr->GetAtomicNum() == OBElements::Hydrogen) { // hydroxyl (YES)
              return(true); 
            }
            else {
              bond = nbr->GetBond(this);
              if ( (bond->IsEster()) && (!(IsCarboxylOxygen()))) 
                 return(false); 
            }
          }
        }
        return(true); // any other oxygen
    } // oxygen END
    // fluorine
    if (_ele == 9 ) {
        OBBondIterator i;
        // organic fluorine (NO)
        for (OBAtom* nbr = BeginNbrAtom(i);nbr;nbr=NextNbrAtom(i))
          if (nbr->GetAtomicNum() == 6)
            return (false);
          else
            return (true);
    };
    if (_ele == 7) {
      // N+ ions and sp2 hybrid N with 3 valences should not be Hbond acceptors
      if (!((GetExplicitDegree() == 4 && GetHyb() == 3)
        || (GetExplicitDegree() == 3 && GetHyb() == 2)))
        return true;
    };
    // Changes from Paolo Tosco
    if (_ele == 16 && GetFormalCharge() == -1){
          return (true); }
    // everything else
    return (false);
  }

  bool OBAtom::IsHbondDonor()
  {
    // return MatchesSMARTS("[$([#8,#7H,#9;!H0])]");
    if (!(_ele == 7 || _ele == 8 || _ele == 9))
      return false;

    OBAtom *nbr;
    OBBondIterator i;
    for (nbr = BeginNbrAtom(i);nbr;nbr = NextNbrAtom(i))
      if (nbr->GetAtomicNum() == OBElements::Hydrogen)
        return true;

    return false;
  }

  bool OBAtom::IsHbondDonorH()
  {
    if (GetAtomicNum() != OBElements::Hydrogen) return(false);

    // Recursive definition -- this atom is an H atom
    // are any neighbors HbondDonors?

    OBAtom *atom;
    OBBond *bond;
    OBBondIterator i;
    for (bond = BeginBond(i);bond;bond = NextBond(i)) {
        atom = bond->GetNbrAtom(this);
      if (atom->IsHbondDonor()) return(true);
      }

    return(false);
  }

  bool OBAtom::IsMetal()
  {
    const unsigned NMETALS = 78;
    const int metals[NMETALS] = {
    3,4,11,12,13,19,20,21,22,23,24,25,26,27,28,29,
    30,31,37,38,39,40,41,42,43,44,45,46,47,48,49,50,55,56,57,58,59,60,61,62,63,
    64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,87,88,89,90,91,
    92,93,94,95,96,97,98,99,100,101,102,103};
    return std::find(metals, metals+78, GetAtomicNum())!=metals+78;
  }

} // end namespace OpenBabel

//! \file atom.cpp
//! \brief Handle OBAtom class

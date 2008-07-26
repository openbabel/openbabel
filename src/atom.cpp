/**********************************************************************
atom.cpp - Handle OBAtom class.
 
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2008 by Geoffrey R. Hutchison
Some portions Copyright (C) 2003 by Michael Banck
 
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
#include <openbabel/babelconfig.h>

#include <openbabel/atom.h>
#include <openbabel/mol.h>
#include <openbabel/molchrg.h>
#include <openbabel/phmodel.h>
#include <openbabel/builder.h>

#include <openbabel/math/matrix3x3.h>

#if !HAVE_STRNCASECMP
extern "C" int strncasecmp(const char *s1, const char *s2, size_t n);
#endif

#ifndef ONE_OVER_SQRT3
#define ONE_OVER_SQRT3  0.57735026918962576451
#endif //SQRT_TWO_THIRDS
#ifndef SQRT_TWO_THIRDS
#define SQRT_TWO_THIRDS 0.81649658092772603272
#endif //ONE_OVER_SQRT3

using namespace std;

namespace OpenBabel
{

  class OBAtomPrivate
  {
    public:
      char                          ele;       //!< atomic number (type char to minimize space, allows for 0..255 elements)
      char                          impval;    //!< implicit valence
      char                          type[6];   //!< atomic type
      short                         fcharge;   //!< formal charge
      unsigned short                isotope;   //!< isotope (0 = most abundant)
      short                         spinmultiplicity;//!< atomic spin, e.g., 2 for radical  1 or 3 for carbene

      OBMol                        *parent;    //!< parent molecule (if any)
      std::vector<OBBond*>          vbond;     //!< bonds to this atom -- assumed to be one of the endpoints

      unsigned short                hyb;       //!< hybridization
      unsigned short                flags;     //!< bitwise flags (e.g. aromaticity)
      double                        pcharge;   //!< partial charge
      OBResidue                    *residue;   //!< parent residue (if applicable)
  };

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

  extern OBAromaticTyper  aromtyper;
  extern OBAtomTyper      atomtyper;
  extern OBPhModel        phmodel;

  //
  // OBAtom protected member functions
  //

  int OBAtom::GetFlag() const
  {
    return d->flags;
  }
      
  void OBAtom::SetFlag(int flag)  
  { 
    d->flags |= flag;   
  }
      
  bool OBAtom::HasFlag(int flag)  
  {  
    return ((d->flags & flag) ? true : false); 
  }

  //
  // OBAtom public member functions
  //

  OBAtom::OBAtom() : d(new OBAtomPrivate)
  {
    m_pos = VZero;
    d->parent = (OBMol*)NULL;
    Clear();
  }

  OBAtom::~OBAtom()
  {
    if (d->residue != NULL)
      {
        d->residue->RemoveAtom(this);
      }
  }

  bool OBAtom::Clear()
  {
    m_cptr = (double**)NULL;
    m_idx = 0;
    m_cidx = 0;
    d->flags=0;
    d->hyb = 0;
    d->ele = (char)0;
    d->isotope = 0;
    d->spinmultiplicity = 0; // CM 18 Sept 2003
    d->impval = 0;
    d->fcharge = 0;
    d->type[0] = '\0';
    d->pcharge = 0.0;
    d->vbond.clear();
    d->vbond.reserve(4);
    d->residue = (OBResidue*)NULL;

    return(OBBase::Clear());
  }

  OBAtom &OBAtom::operator=(OBAtom &src)
  //copy atom information
  //bond info is not copied here as ptrs may be invalid
  {
    m_idx = src.GetIdx();
    Duplicate(&src);
    return(*this);
  }
  
  void OBAtom::Duplicate(OBAtom *src)
  {
    if (!src)
      return;
      
    d->hyb = src->GetHyb();
    d->ele = src->GetAtomicNum();
    d->isotope = src->GetIsotope();
    d->fcharge = src->GetFormalCharge();
    d->spinmultiplicity = src->GetSpinMultiplicity();
    strncpy(d->type,src->GetType(), sizeof(d->type) - 1);
    d->type[sizeof(d->type) - 1] = '\0';
    d->pcharge = src->GetPartialCharge();
    m_pos = src->GetVector();
    d->flags = src->GetFlag();
    d->residue = (OBResidue*)NULL;
    
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
      if (includePandS || (!a1->IsPhosphorus() && !a1->IsSulfur()))
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
    if (!IsHydrogen())
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
    if (!IsHydrogen())
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

  void OBAtom::SetIdx(unsigned int idx)    
  { 
    m_idx = idx; 
    m_cidx = (idx - 1) * 3; 
  }
      
  void OBAtom::SetHyb(int hyb)    
  { 
    d->hyb = hyb; 
  }
      
  void OBAtom::SetAtomicNum(int atomicnum)
  { 
    d->ele = (char)atomicnum; 
  }
      
  void OBAtom::SetImplicitValence(int val)
  { 
    d->impval = (char)val; 
  }
      
  void OBAtom::IncrementImplicitValence()     
  { 
    d->impval++; 
  }
      
  void OBAtom::DecrementImplicitValence()     
  { 
    d->impval--; 
  }
      
  void OBAtom::SetFormalCharge(int fcharge)   
  { 
    d->fcharge = fcharge; 
  }
      
  void OBAtom::SetSpinMultiplicity(short spin)
  { 
    d->spinmultiplicity = spin; 
  }
      
  void OBAtom::SetPartialCharge(double pcharge)
  { 
    d->pcharge = pcharge; 
  }
      
  void OBAtom::SetCoordPtr(double **c)        
  { 
    m_cptr = c; 
    m_cidx = (GetIdx() - 1) * 3; 
  }

  void OBAtom::SetResidue(OBResidue *res)     
  { 
    d->residue=res; 
  }

  void OBAtom::SetParent(OBMol *ptr)
  { 
    d->parent=ptr; 
  }
      
  void OBAtom::UnsetAromatic()                
  { 
    d->flags &= (~(OBAtomFlag::Aromatic)); 
  }
      
  void OBAtom::UnsetStereo()
  {
    d->flags &= ~(OBAtomFlag::StereoC);
    d->flags &= ~(OBAtomFlag::StereoAC);
    d->flags &= ~(OBAtomFlag::ChiralPos);
    d->flags &= ~(OBAtomFlag::ChiralNeg);
    d->flags &= ~(OBAtomFlag::Chiral);
  }
      
  void OBAtom::ClearCoordPtr()     
  { 
    m_cptr = NULL; 
    m_cidx = 0; 
  }
  
  int OBAtom::GetFormalCharge() const 
  { 
    return d->fcharge;
  }
      
  unsigned int OBAtom::GetAtomicNum() const 
  { 
    return (unsigned int) d->ele; 
  }
      
  unsigned short int OBAtom::GetIsotope() const 
  { 
    return d->isotope;    
  }
      
  int OBAtom::GetSpinMultiplicity() const 
  { 
    return d->spinmultiplicity; 
  }
     
  /* 
  unsigned int OBAtom::GetIdx() const 
  { 
    return (int) d->idx;
  }
      
  unsigned int OBAtom::GetCoordinateIdx() const 
  { 
    return (int) d->cidx; 
  }
  */
      
  unsigned int OBAtom::GetValence() const
  {
    return ((d->vbond.empty()) ? 0 : d->vbond.size());
  }

  double OBAtom::x() const 
  {
    if (m_cptr) return((*m_cptr)[m_cidx]);
    else        return m_pos.x();
  }
  double OBAtom::y() const 
  {
    if (m_cptr) return((*m_cptr)[m_cidx+1]);
    else        return m_pos.y();
  }
  double OBAtom::z() const 
  {
    if (m_cptr) return((*m_cptr)[m_cidx+2]);
    else        return m_pos.z();
  }
  
  /*
  double* OBAtom::GetCoordinate()
  {
    if (d->c) return(&(*d->c)[d->cidx]);
    else      return NULL;
  }
  */
  
  OBMol* OBAtom::GetParent()        
  {
    return (OBMol*) d->parent;
  }
      
  OBBondIterator OBAtom::BeginBonds()
  { 
    return d->vbond.begin(); 
  }
  
  OBBondIterator OBAtom::EndBonds()
  { 
    return d->vbond.end();   
  }
  
  void OBAtom::NewResidue()
  {
    if (!d->residue)
      d->residue = new OBResidue;
  }
  
  void OBAtom::DeleteResidue(){
    if (d->residue) {
      delete d->residue;
      d->residue = NULL; // Make sure to clear that a residue existed
    }
  }
  
  void OBAtom::AddBond(OBBond *bond) 
  { 
    d->vbond.push_back(bond); 
  }
  
  void OBAtom::InsertBond(OBBondIterator &i, OBBond *bond)
  {
    d->vbond.insert(i, bond);
  }
  
  void OBAtom::ClearBond() 
  {
    d->vbond.clear();
  }

  bool OBAtom::HasResidue()
  { 
    return (d->residue != NULL);    
  }
  
  void OBAtom::SetVector()
  {
    if (m_cptr)
      m_pos = Eigen::Vector3d( (*m_cptr)[m_cidx], (*m_cptr)[m_cidx+1], (*m_cptr)[m_cidx+2] );
  }

  void OBAtom::SetVector(const Eigen::Vector3d &v)
  {
    if (!m_cptr)
      m_pos = v;
    else
      {
        (*m_cptr)[m_cidx  ] = v.x();
        (*m_cptr)[m_cidx+1] = v.y();
        (*m_cptr)[m_cidx+2] = v.z();
      }
  }

  void OBAtom::SetVector(const double x,const double y,const double z)
  {
    if (!m_cptr)
      m_pos = Eigen::Vector3d(x, y, z);
    else
      {
        (*m_cptr)[m_cidx  ] = x;
        (*m_cptr)[m_cidx+1] = y;
        (*m_cptr)[m_cidx+2] = z;
      }
  }

  void OBAtom::SetType(const char *type)
  {
    strncpy(d->type,type, sizeof(d->type) - 1);
    d->type[sizeof(d->type) - 1] = '\0';
    if (d->ele == 1 && type[0] == 'D')
      d->isotope = 2;
  }

  void OBAtom::SetType(const std::string &type)
  {
    strncpy(d->type,type.c_str(), sizeof(d->type) - 1);
    d->type[sizeof(d->type) - 1] = '\0';
    if (d->ele == 1 && type[0] == 'D')
      d->isotope = 2;
  }

  void OBAtom::SetIsotope(unsigned int iso)
  {
    if (d->ele == 1 && iso == 2)
      SetType("D");
    else if (d->ele == 1 && (iso == 1 || iso == 0))
      SetType("H");

    d->isotope = iso;
  }

  OBResidue *OBAtom::GetResidue()
  {
    return GetResidue(true); // default is to always perceive chains
  }

  OBResidue *OBAtom::GetResidue(bool perception)
  {
    if (d->residue != NULL)
      return d->residue;
    else if (perception && !((OBMol*)GetParent())->HasChainsPerceived())
      {
        ((OBMol*)GetParent())->SetChainsPerceived();
        if ( chainsparser.PerceiveChains(*((OBMol*)GetParent())) )
          return d->residue;
        else
          {
            if (d->residue)
              {
                delete d->residue;
                d->residue = NULL;
              }
            return NULL;
          }
      }
    else
      return NULL;
  }

  double OBAtom::GetAtomicMass() const
  {
    if (d->isotope == 0)
      return etab.GetMass(d->ele);
    else
      return isotab.GetExactMass(d->ele, d->isotope);
  }

  double OBAtom::GetExactMass() const
  {
    return isotab.GetExactMass(d->ele, d->isotope);
  }

  char *OBAtom::GetType()
  {
    OBMol *mol = (OBMol*)GetParent();
    if (mol)
      if (!mol->HasAtomTypesPerceived())
        atomtyper.AssignTypes(*((OBMol*)GetParent()));

    if (strlen(d->type) == 0) // Somehow we still don't have a type!
      {
        char num[6];
        string fromType = ttab.GetFromType(); // save previous types
        string toType = ttab.GetToType();

        ttab.SetFromType("ATN");
        ttab.SetToType("INT");
        snprintf(num, 6, "%d", GetAtomicNum());
        ttab.Translate(d->type, num);
	
        ttab.SetFromType(fromType.c_str());
        ttab.SetToType(toType.c_str());
      }
    if (d->ele == 1 && d->isotope == 2)
      snprintf(d->type, 6, "%s", "D");

    return(d->type);
  }

  unsigned int OBAtom::GetImplicitValence() const
  {
    OBMol *mol = (OBMol*)((OBAtom*)this)->GetParent();
    if (mol && !mol->HasImplicitValencePerceived())
      atomtyper.AssignImplicitValence(*((OBMol*)((OBAtom*)this)->GetParent()));

    return((unsigned int)d->impval);
  }

  unsigned int OBAtom::GetHyb() const
  {
    //hybridization is assigned when atoms are typed
    OBMol *mol = (OBMol*)((OBAtom*)this)->GetParent();
    if (mol && !mol->HasHybridizationPerceived())
      atomtyper.AssignHyb(*mol);

    return(d->hyb);
  }


  unsigned int OBAtom::GetHvyValence() const
  {
    unsigned int count=0;

    OBBond *bond;
    OBBondIterator i;
    for (bond = ((OBAtom*)this)->BeginBond(i); bond; bond = ((OBAtom*)this)->NextBond(i))
      if (!(bond->GetNbrAtom((OBAtom*)this)->IsHydrogen()))
        count++;

    return(count);
  }

  unsigned int OBAtom::GetHeteroValence() const
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
      return(d->pcharge);
    if (!((OBMol*)GetParent())->AutomaticPartialCharge())
      return(d->pcharge);

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

    return(d->pcharge);
  }

  //! Returns true if nitrogen is part of an amide
  bool OBAtom::IsAmideNitrogen()
  {
    if (!IsNitrogen())
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
    if (!IsNitrogen() || !IsAromatic())
      return(false);

    OBAtom *atom;
    OBBondIterator i;

    for (atom = BeginNbrAtom(i);atom;atom = NextNbrAtom(i))
      if (atom->IsOxygen() && !(*i)->IsInRing() && (*i)->GetBondOrder() == 2)
        return(true);

    return(false);
  }

  bool OBAtom::IsCarboxylOxygen()
  {
    if (!IsOxygen())
      return(false);
    if (GetHvyValence() != 1)
      return(false);

    OBAtom *atom;
    OBBond *bond;
    OBBondIterator i;

    atom = NULL;
    for (bond = BeginBond(i);bond;bond = NextBond(i))
      if ((bond->GetNbrAtom(this))->IsCarbon())
        {
          atom = bond->GetNbrAtom(this);
          break;
        }
    if (!atom)
      return(false);
    if (atom->CountFreeOxygens() != 2)
      return(false);

    //atom is connected to a carbon that has a total
    //of 2 attached free oxygens
    return(true);
  }

  bool OBAtom::IsPhosphateOxygen()
  {
    if (!IsOxygen())
      return(false);
    if (GetHvyValence() != 1)
      return(false);

    OBAtom *atom;
    OBBond *bond;
    OBBondIterator i;

    atom = NULL;
    for (bond = BeginBond(i);bond;bond = NextBond(i))
      if ((bond->GetNbrAtom(this))->IsPhosphorus())
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
    if (!IsOxygen())
      return(false);
    if (GetHvyValence() != 1)
      return(false);

    OBAtom *atom;
    OBBond *bond;
    OBBondIterator i;

    atom = NULL;
    for (bond = BeginBond(i);bond;bond = NextBond(i))
      if ((bond->GetNbrAtom(this))->IsSulfur())
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

  bool OBAtom::IsNitroOxygen()
  {
    if (!IsOxygen())
      return(false);
    if (GetHvyValence() != 1)
      return(false);

    OBAtom *atom;
    OBBond *bond;
    OBBondIterator i;

    atom = NULL;
    for (bond = BeginBond(i);bond;bond = NextBond(i))
      if ((bond->GetNbrAtom(this))->IsNitrogen())
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

  bool OBAtom::IsNotCorH()
  {
    switch(GetAtomicNum())
      {
      case 1:
      case 6:
        return(false);
      }
    return(true);
  }

  bool OBAtom::IsAromatic() const
  {
    if (((OBAtom*)this)->HasFlag(OBAtomFlag::Aromatic))
      return(true);

    OBMol *mol = (OBMol*)((OBAtom*)this)->GetParent();

    if (!mol->HasAromaticPerceived())
      {
        aromtyper.AssignAromaticFlags(*mol);
        if (((OBAtom*)this)->HasFlag(OBAtomFlag::Aromatic))
          return(true);
      }

    return(false);
  }

  bool OBAtom::IsInRing() const
  {
    if (((OBAtom*)this)->HasFlag(OBAtomFlag::Ring))
      return(true);

    OBMol *mol = (OBMol*)((OBAtom*)this)->GetParent();
    if (!mol->HasRingAtomsAndBondsPerceived())
      {
        mol->FindRingAtomsAndBonds();
        if (((OBAtom*)this)->HasFlag(OBAtomFlag::Ring))
          return(true);
      }

    return(false);
  }

  bool OBAtom::IsChiral()
  {
    if (HasFlag(OBAtomFlag::Chiral))
      return(true);

    if (!((OBMol*)GetParent())->HasChiralityPerceived())
      {
        ((OBMol*)GetParent())->FindChiralCenters();
        if (HasFlag(OBAtomFlag::Chiral))
          return(true);
      }

    return(false);
  }

  bool OBAtom::IsInRingSize(int size) const
  {
    vector<OBRing*> rlist;
    vector<OBRing*>::iterator i;

    OBMol *mol = (OBMol*)((OBAtom*)this)->GetParent();
    if (!mol->HasSSSRPerceived())
      mol->FindSSSR();

    if (!((OBAtom*)this)->HasFlag(OBAtomFlag::Ring))
      return(false);

    rlist = mol->GetSSSR();
    for (i = rlist.begin();i != rlist.end();++i)
      if ((*i)->IsInRing(GetIdx()) && (*i)->PathSize() == size)
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

  double OBAtom::SmallestBondAngle()
  {
    OBAtom *b, *c;
    Eigen::Vector3d v1, v2;
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
            degrees = VectorAngle(v1, v2);
            if (degrees < minDegrees)
              minDegrees = degrees;
          }
      }
    return minDegrees;
  }

  double OBAtom::AverageBondAngle()
  {
    OBAtom *b, *c;
    Eigen::Vector3d v1, v2;
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
            degrees = VectorAngle(v1, v2);
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
        if (atom->IsOxygen() && atom->GetHvyValence() == 1)
          count++;
      }

    return(count);
  }

  unsigned int OBAtom::BOSum() const
  {
    unsigned int bo;
    unsigned int bosum = 0;
    OBBond *bond;
    OBBondIterator i;

    for (bond = ((OBAtom*)this)->BeginBond(i);bond;bond = ((OBAtom*)this)->NextBond(i))
      {
        bo = bond->GetBondOrder();
        bosum += (bo < 4) ? 2*bo : 3;
      }

    bosum /= 2;
    return(bosum);
  }

  unsigned int OBAtom::ImplicitHydrogenCount() const
  {
    //handles H,C,N,S,O,X
    OBMol *mol = (OBMol*)((OBAtom*)this)->GetParent();
    if (mol && !mol->HasImplicitValencePerceived())
      atomtyper.AssignImplicitValence(*((OBMol*)((OBAtom*)this)->GetParent()));

    // d->impval is assigned by the atomtyper -- same as calling GetImplicitValence()
    int impval = d->impval - GetValence();

    // we need to modify this implicit valence if we're a radical center
    int mult = GetSpinMultiplicity();
    if(mult==2) //radical
      impval-=1;
    else if(mult==1 || mult==3) //carbene
      impval-=2;
    else if(mult>=4) //CH, Catom
      impval -= mult-1;
    return((impval>0)?impval:0);
  }

  unsigned int OBAtom::ExplicitHydrogenCount(bool ExcludeIsotopes) const
  {
    //If ExcludeIsotopes is true, H atoms with d->isotope!=0 are not included.
    //This excludes D, T and H when d->isotope exlicitly set to 1 rather than the default 0.
    int numH=0;
    OBAtom *atom;
    OBBondIterator i;
    for (atom = ((OBAtom*)this)->BeginNbrAtom(i);atom;atom = ((OBAtom*)this)->NextNbrAtom(i))
      if (atom->IsHydrogen() && !(ExcludeIsotopes && atom->GetIsotope()!=0))
        numH++;

    return(numH);
  }

  bool OBAtom::DeleteBond(OBBond *bond)
  {
    OBBondIterator i;
    for (i = d->vbond.begin();i != d->vbond.end();++i)
      if ((OBBond*)bond == *i)
        {
          d->vbond.erase(i);
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
    i = d->vbond.begin();
    return((i == d->vbond.end()) ? (OBBond*)NULL : (OBBond*)*i);
  }

  OBBond *OBAtom::NextBond(OBBondIterator &i)
  {
    i++;
    return((i == d->vbond.end()) ? (OBBond*)NULL : (OBBond*)*i);
  }

  OBAtom *OBAtom::BeginNbrAtom(OBBondIterator &i)
  {
    i = d->vbond.begin();
    return((i != d->vbond.end()) ? ((OBBond*) *i)->GetNbrAtom(this):NULL);
  }

  OBAtom *OBAtom::NextNbrAtom(OBBondIterator &i)
  {
    i++;
    return((i != d->vbond.end()) ?  ((OBBond*) *i)->GetNbrAtom(this):NULL);
  }

  double OBAtom::GetDistance(OBAtom *b)
  {
    return(( this->GetVector() - b->GetVector() ).norm());
  }

  double OBAtom::GetDistance(unsigned int b)
  {
    OBMol *mol = (OBMol*)GetParent();
    return(( this->GetVector() - mol->GetAtom(b)->GetVector() ).norm());
  }

  double OBAtom::GetAngle(OBAtom *b, OBAtom *c)
  {
    Eigen::Vector3d v1,v2;

    v1 = this->GetVector() - b->GetVector();
    v2 = c->GetVector() - b->GetVector();
    if (IsNearZero(v1.norm(), 1.0e-3)
      || IsNearZero(v2.norm(), 1.0e-3)) {
        return(0.0);
    }

    return(VectorAngle(v1, v2));
  }

  double OBAtom::GetAngle(int b, int c)
  {
    OBMol *mol = (OBMol*)GetParent();
    Eigen::Vector3d v1,v2;

    v1 = this->GetVector() - mol->GetAtom(b)->GetVector();
    v2 = mol->GetAtom(c)->GetVector() - mol->GetAtom(b)->GetVector();

    if (IsNearZero(v1.norm(), 1.0e-3)
      || IsNearZero(v2.norm(), 1.0e-3)) {
        return(0.0);
    }

    return(VectorAngle(v1, v2));
  }

  bool OBAtom::HtoMethyl()
  {
    if (!IsHydrogen())
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
    br1 = etab.CorrectedBondRad(6,3);
    br2 = etab.CorrectedBondRad(atom->GetAtomicNum(),atom->GetHyb());
    bond->SetLength(atom,br1+br2);

    OBAtom *hatom;
    br2 = etab.CorrectedBondRad(1,0);
    Eigen::Vector3d v;

    for (int j = 0;j < 3;++j)
      {
        hatom = mol->NewAtom();
        hatom->SetAtomicNum(1);
        hatom->SetType("H");

	v = OBBuilder::GetNewBondVector(hatom, br1 + br2);
        hatom->SetVector(v);
        mol->AddBond(GetIdx(),mol->NumAtoms(),1);
      }

    mol->EndModify();
    return(true);
  }

  static void ApplyRotMatToBond(OBMol &mol, Eigen::Quaternion<double> &m, OBAtom *a1, OBAtom *a2)
  {
    vector<int> children;
    mol.FindChildren(children,a1->GetIdx(),a2->GetIdx());
    children.push_back(a2->GetIdx());
                
    Eigen::Vector3d v;
    vector<int>::iterator i;
    for (i = children.begin();i != children.end();++i)
      {
        v = mol.GetAtom(*i)->GetVector();
        v -= a1->GetVector();
        v = m.toRotationMatrix() * v;
        v += a1->GetVector();
        mol.GetAtom(*i)->SetVector(v);
      }

  }

  bool OBAtom::SetHybAndGeom(int hyb)
  {
    obErrorLog.ThrowError(__FUNCTION__,
                          "Ran OpenBabel::SetHybridizationAndGeometry",
                          obAuditMsg);

    //if (hyb == GetHyb()) return(true);
    if (GetAtomicNum() == 1)
      return(false);
    if (hyb == 0 && GetHvyValence() > 1)
      return(false);
    if (hyb == 1 && GetHvyValence() > 2)
      return(false);
    if (hyb == 2 && GetHvyValence() > 3)
      return(false);
    if (hyb == 3 && GetHvyValence() > 4)
      return(false);

    OBMol *mol = (OBMol*)GetParent();

    OBAtom *nbr;
    vector<OBAtom*> delatm;
    OBBondIterator i;

    for (nbr = BeginNbrAtom(i);nbr;nbr = NextNbrAtom(i))
      if (nbr->IsHydrogen())
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
    br1 = etab.CorrectedBondRad(GetAtomicNum(),hyb);
    for (nbr = BeginNbrAtom(i);nbr;nbr = NextNbrAtom(i))
      if (!(*i)->IsInRing())
        {
          br2 = etab.CorrectedBondRad(nbr->GetAtomicNum(),nbr->GetHyb());
          length = br1 + br2;
          if ((*i)->IsAromatic())
            length *= 0.93;
          else if ((*i)->GetBondOrder() == 2)
            length *= 0.91;
          else if ((*i)->GetBondOrder() == 3)
            length *= 0.87;
          ((OBBond*) *i)->SetLength(this, length);
        }

    if (GetValence() > 1)
      {
        double angle;
        Eigen::Quaternion<double> q;
        Eigen::Vector3d v1,v2,v3,v4,n,s;
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
                n = v1.cross(v2);
                angle = VectorAngle(v1,v2) - targetAngle;
                q = Eigen::AngleAxis<double>(angle * DEG_TO_RAD, n);
                ApplyRotMatToBond(*mol, q, this, a1);
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
                n = s.cross(v3);
                angle = VectorAngle(s,v3);
                q = Eigen::AngleAxis<double>(-angle * DEG_TO_RAD, n);
                ApplyRotMatToBond(*mol, q, this, a1);
              }
            else
              {
                if (a2)
                  {
                    v1 = a1->GetVector()-GetVector();
                    v1.normalize();
                    v2 = a2->GetVector()-GetVector();
                    v2.normalize();
                    n = v1.cross(v2);
                    angle = VectorAngle(v1,v2)-targetAngle;
                    q = Eigen::AngleAxis<double>(angle * DEG_TO_RAD, n);
                    ApplyRotMatToBond(*mol, q, this, a1);
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
                    n = s.cross(v3);
                    angle = VectorAngle(s,v3);
                    q = Eigen::AngleAxis<double>(-angle * DEG_TO_RAD, n);
                    ApplyRotMatToBond(*mol, q, this, a3);
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
                n = s.cross(v4);
                angle = VectorAngle(s,v4);
                q = Eigen::AngleAxis<double>(-angle * DEG_TO_RAD, n);
                ApplyRotMatToBond(*mol, q, this, a1);
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
                n = v1.cross(v2);
                n.normalize();
                s *= ONE_OVER_SQRT3; //1/sqrt(3)
                n *= SQRT_TWO_THIRDS; //sqrt(2/3)
                s += n;
                s.normalize();
                n = s.cross(v3);
                angle = VectorAngle(s,v3);
                q = Eigen::AngleAxis<double>(-angle * DEG_TO_RAD, n);
                ApplyRotMatToBond(*mol, q, this, a1);
 
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
                    n = s.cross(v4);
                    angle = VectorAngle(s,v4);
                    q = Eigen::AngleAxis<double>(-angle * DEG_TO_RAD, n);
                    ApplyRotMatToBond(*mol, q, this, a2);
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
                    n = v1.cross(v2);
                    angle = VectorAngle(v1,v2)-targetAngle;
                    q = Eigen::AngleAxis<double>(angle * DEG_TO_RAD, n);
                    ApplyRotMatToBond(*mol, q, this, a1);
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
                    n = v1.cross(v2);
                    n.normalize();
                    s *= ONE_OVER_SQRT3; //1/sqrt(3)
                    n *= SQRT_TWO_THIRDS; //sqrt(2/3)
                    s += n;
                    s.normalize();
                    n = s.cross(v3);
                    angle = VectorAngle(s,v3);
                    q = Eigen::AngleAxis<double>(-angle * DEG_TO_RAD, n);
                    ApplyRotMatToBond(*mol, q, this, a3);
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

    int hcount = impval-GetHvyValence();
    if (hcount)
      {
        int k;
        Eigen::Vector3d v;
        OBAtom *atom;
        double brsum = etab.CorrectedBondRad(1,0)+
          etab.CorrectedBondRad(GetAtomicNum(),GetHyb());
        SetHyb(hyb);

        mol->BeginModify();
        for (k = 0;k < hcount;++k)
          {
	    v = OBBuilder::GetNewBondVector(this, brsum);
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

  bool OBAtom::IsHbondAcceptor()
  {
    // Changes from Liu Zhiguo
    if (d->ele == 8 || d->ele == 9)
      return true;
    if (d->ele == 7) {
      // N+ ions and sp2 hybrid N with 3 valences should not be Hbond acceptors
      if (!((GetValence() == 4 && GetHyb() == 3) 
            || (GetValence() == 3 && GetHyb() == 2)))
            return true;
    }
    return false;
  }

  bool OBAtom::IsHbondDonor()
  {
    // return MatchesSMARTS("[$([#8,#7H,#9;!H0])]");
    if (!(d->ele == 7 || d->ele == 8 || d->ele == 9))
      return false;

    OBAtom *nbr;
    OBBondIterator i;
    for (nbr = BeginNbrAtom(i);nbr;nbr = NextNbrAtom(i))
      if (nbr->IsHydrogen())
        return true;

    return false;
  }

  bool OBAtom::IsHbondDonorH()
  {
    if (!IsHydrogen()) return(false);

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

  unsigned int OBAtom::KBOSum() const
  {
    OBBond *bond;
    unsigned int bosum;
    OBBondIterator i;

    bosum = GetImplicitValence();

    for (bond = ((OBAtom*)this)->BeginBond(i);bond;bond = ((OBAtom*)this)->NextBond(i))
      {
        if (bond->IsKDouble())
          bosum++;
        else if (bond->IsKTriple())
          bosum += 2;
      }

    return(bosum);
  }

} // end namespace OpenBabel

//! \file atom.cpp
//! \brief Handle OBAtom class

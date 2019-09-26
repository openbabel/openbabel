/********************************************************************** 
mol.cpp - Handle molecules.

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

#include <openbabel/mol.h>
#include <openbabel/bond.h>
#include <openbabel/ring.h>
#include <openbabel/rotamer.h>
#include <openbabel/phmodel.h>
#include <openbabel/bondtyper.h>
#include <openbabel/obiter.h>
#include <openbabel/builder.h>
#include <openbabel/kekulize.h>
#include <openbabel/internalcoord.h>
#include <openbabel/math/matrix3x3.h>
#include <openbabel/obfunctions.h>
#include <openbabel/elements.h>

#include <openbabel/stereo/tetrahedral.h>
#include <openbabel/stereo/cistrans.h>

#include <sstream>
#include <set>

using namespace std;

namespace OpenBabel
{

  extern bool SwabInt;
  extern THREAD_LOCAL OBPhModel  phmodel;
  extern THREAD_LOCAL OBAromaticTyper  aromtyper;
  extern THREAD_LOCAL OBAtomTyper      atomtyper;
  extern THREAD_LOCAL OBBondTyper      bondtyper;
  
  /** \class OBMol mol.h <openbabel/mol.h>
      \brief Molecule Class

      The most important class in Open Babel is OBMol, or the molecule class.
      The OBMol class is designed to store all the basic information
      associated with a molecule, to make manipulations on the connection
      table of a molecule facile, and to provide member functions which
      automatically perceive information about a molecule. A guided tour
      of the OBMol class is a good place to start.

      An OBMol class can be declared:
      \code
      OBMol mol;
      \endcode

      For example:
      \code
      #include <iostream.h>

      #include <openbabel/mol.h>
      #include <openbabel/obconversion.h>
      int main(int argc,char **argv)
      {
      OBConversion conv(&cin,&cout);
      if(conv.SetInAndOutFormats("SDF","MOL2"))
      {
      OBMol mol;
      if(conv.Read(&mol))
      ...manipulate molecule

      conv->Write(&mol);
      }
      return(1);
      }
      \endcode

      will read in a molecule in SD file format from stdin
      (or the C++ equivalent cin) and write a MOL2 format file out
      to standard out. Additionally, The input and output formats can
      be altered using the OBConversion class

      Once a molecule has been read into an OBMol (or created via other methods)
      the atoms and bonds
      can be accessed by the following methods:
      \code
      OBAtom *atom;
      atom = mol.GetAtom(5); //random access of an atom
      \endcode
      or
      \code
      OBBond *bond;
      bond = mol.GetBond(14); //random access of a bond
      \endcode
      or
      \code
      FOR_ATOMS_OF_MOL(atom, mol) // iterator access (see OBMolAtomIter)
      \endcode
      or
      \code
      FOR_BONDS_OF_MOL(bond, mol) // iterator access (see OBMolBondIter)
      \endcode
      It is important to note that atom arrays currently begin at 1 and bond arrays
      begin at 0. Requesting atom 0 (\code
      OBAtom *atom = mol.GetAtom(0); \endcode
      will result in an error, but
      \code
      OBBond *bond = mol.GetBond(0);
      \endcode
      is perfectly valid.
      Note that this is expected to change in the near future to simplify coding
      and improve efficiency.

      The ambiguity of numbering issues and off-by-one errors led to the use
      of iterators in Open Babel. An iterator is essentially just a pointer, but
      when used in conjunction with Standard Template Library (STL) vectors
      it provides an unambiguous way to loop over arrays. OBMols store their
      atom and bond information in STL vectors. Since vectors are template
      based, a vector of any user defined type can be declared. OBMols declare
      vector<OBAtom*> and vector<OBBond*> to store atom and bond information.
      Iterators are then a natural way to loop over the vectors of atoms and bonds.

      A variety of predefined iterators have been created to simplify
      common looping requests (e.g., looping over all atoms in a molecule,
      bonds to a given atom, etc.)

      \code
      #include <openbabel/obiter.h>
      ...
      #define FOR_ATOMS_OF_MOL(a,m)     for( OBMolAtomIter     a(m); a; ++a )
      #define FOR_BONDS_OF_MOL(b,m)     for( OBMolBondIter     b(m); b; ++b )
      #define FOR_NBORS_OF_ATOM(a,p)    for( OBAtomAtomIter    a(p); a; ++a )
      #define FOR_BONDS_OF_ATOM(b,p)    for( OBAtomBondIter    b(p); b; ++b )
      #define FOR_RESIDUES_OF_MOL(r,m)  for( OBResidueIter     r(m); r; ++r )
      #define FOR_ATOMS_OF_RESIDUE(a,r) for( OBResidueAtomIter a(r); a; ++a )
      ...
      \endcode

      These convenience functions can be used like so:
      \code
      #include <openbabel/obiter.h>
      #include <openbabel/mol.h>

      OBMol mol;
      double exactMass = 0.0;
      FOR_ATOMS_OF_MOL(a, mol)
      {
      exactMass +=  a->GetExactMass();
      }
      \endcode

      Note that with these convenience macros, the iterator "a" (or
      whichever name you pick) is declared for you -- you do not need to
      do it beforehand.
  */

  //
  // OBMol member functions
  //
  void  OBMol::SetTitle(const char *title)
  {
    _title = title;
    Trim(_title);
  }

  void  OBMol::SetTitle(std::string &title)
  {
    _title = title;
    Trim(_title);
  }

  const char *OBMol::GetTitle(bool replaceNewlines) const
  {
    if (!replaceNewlines || _title.find('\n')== string::npos )
      return(_title.c_str());

    //Only multiline titles use the following to replace newlines by spaces
    static string title;
    title=_title;
    string::size_type j;
    for ( ; (j = title.find_first_of( "\n\r" )) != string::npos ; ) {
      title.replace( j, 1, " ");
    }

    return(title.c_str());
  }

  bool SortVVInt(const vector<int> &a,const vector<int> &b)
  {
    return(a.size() > b.size());
  }

  bool SortAtomZ(const pair<OBAtom*,double> &a, const pair<OBAtom*,double> &b)
  {
    return (a.second < b.second);
  }

  double OBMol::GetAngle( OBAtom* a, OBAtom* b, OBAtom* c)
  {
    return a->GetAngle( b, c );
  }

  double OBMol::GetTorsion(int a,int b,int c,int d)
  {
    return(CalcTorsionAngle(((OBAtom*)_vatom[a-1])->GetVector(),
                            ((OBAtom*)_vatom[b-1])->GetVector(),
                            ((OBAtom*)_vatom[c-1])->GetVector(),
                            ((OBAtom*)_vatom[d-1])->GetVector()));
  }

  void OBMol::SetTorsion(OBAtom *a,OBAtom *b,OBAtom *c, OBAtom *d, double ang)
  {
    vector<int> tor;
    vector<int> atoms;

    obErrorLog.ThrowError(__FUNCTION__,
                          "Ran OpenBabel::SetTorsion", obAuditMsg);

    tor.push_back(a->GetCoordinateIdx());
    tor.push_back(b->GetCoordinateIdx());
    tor.push_back(c->GetCoordinateIdx());
    tor.push_back(d->GetCoordinateIdx());

    FindChildren(atoms, b->GetIdx(), c->GetIdx());
    int j;
    for (j = 0 ; (unsigned)j < atoms.size() ; j++ )
      atoms[j] = (atoms[j] - 1) * 3;

    double v2x,v2y,v2z;
    double radang,m[9];
    double x,y,z,mag,rotang,sn,cs,t,tx,ty,tz;

    //calculate the torsion angle
    radang = CalcTorsionAngle(a->GetVector(),
                              b->GetVector(),
                              c->GetVector(),
                              d->GetVector()) / RAD_TO_DEG;
    //
    // now we have the torsion angle (radang) - set up the rot matrix
    //

    //find the difference between current and requested
    rotang = ang - radang;

    sn = sin(rotang);
    cs = cos(rotang);
    t = 1 - cs;

    v2x = _c[tor[1]]   - _c[tor[2]];
    v2y = _c[tor[1]+1] - _c[tor[2]+1];
    v2z = _c[tor[1]+2] - _c[tor[2]+2];

   //normalize the rotation vector
    mag = sqrt(SQUARE(v2x)+SQUARE(v2y)+SQUARE(v2z));
    x = v2x/mag;
    y = v2y/mag;
    z = v2z/mag;

    //set up the rotation matrix
    m[0]= t*x*x + cs;
    m[1] = t*x*y + sn*z;
    m[2] = t*x*z - sn*y;
    m[3] = t*x*y - sn*z;
    m[4] = t*y*y + cs;
    m[5] = t*y*z + sn*x;
    m[6] = t*x*z + sn*y;
    m[7] = t*y*z - sn*x;
    m[8] = t*z*z + cs;

    //
    //now the matrix is set - time to rotate the atoms
    //
    tx = _c[tor[1]];
    ty = _c[tor[1]+1];
    tz = _c[tor[1]+2];
    vector<int>::iterator i;
    for (i = atoms.begin(); i != atoms.end(); ++i)
      {
        j = *i;

        _c[j] -= tx;
        _c[j+1] -= ty;
        _c[j+2]-= tz;
        x = _c[j]*m[0] + _c[j+1]*m[1] + _c[j+2]*m[2];
        y = _c[j]*m[3] + _c[j+1]*m[4] + _c[j+2]*m[5];
        z = _c[j]*m[6] + _c[j+1]*m[7] + _c[j+2]*m[8];
        _c[j] = x;
        _c[j+1] = y;
        _c[j+2] = z;
        _c[j] += tx;
        _c[j+1] += ty;
        _c[j+2] += tz;
      }
  }


  double OBMol::GetTorsion(OBAtom *a,OBAtom *b,OBAtom *c,OBAtom *d)
  {
    return(CalcTorsionAngle(a->GetVector(),
                            b->GetVector(),
                            c->GetVector(),
                            d->GetVector()));
  }

  void OBMol::ContigFragList(std::vector<std::vector<int> >&cfl)
  {
    int j;
    OBAtom *atom;
    OBBond *bond;
    vector<OBAtom*>::iterator i;
    vector<OBBond*>::iterator k;
    OBBitVec used,curr,next,frag;
    vector<int> tmp;

    used.Resize(NumAtoms()+1);
    curr.Resize(NumAtoms()+1);
    next.Resize(NumAtoms()+1);
    frag.Resize(NumAtoms()+1);

    while ((unsigned)used.CountBits() < NumAtoms())
      {
        curr.Clear();
        frag.Clear();
        for (atom = BeginAtom(i);atom;atom = NextAtom(i))
          if (!used.BitIsSet(atom->GetIdx()))
            {
              curr.SetBitOn(atom->GetIdx());
              break;
            }

        frag |= curr;
        while (!curr.IsEmpty())
          {
            next.Clear();
            for (j = curr.NextBit(-1);j != curr.EndBit();j = curr.NextBit(j))
              {
                atom = GetAtom(j);
                for (bond = atom->BeginBond(k);bond;bond = atom->NextBond(k))
                  if (!used.BitIsSet(bond->GetNbrAtomIdx(atom)))
                    next.SetBitOn(bond->GetNbrAtomIdx(atom));
              }

            used |= curr;
            used |= next;
            frag |= next;
            curr = next;
          }

        tmp.clear();
        frag.ToVecInt(tmp);
        cfl.push_back(tmp);
      }

    sort(cfl.begin(),cfl.end(),SortVVInt);
  }

  void OBMol::FindAngles()
  {
    //if already has data return
    if(HasData(OBGenericDataType::AngleData))
      return;

    //get new data and attach it to molecule
    OBAngleData *angles = new OBAngleData;
    angles->SetOrigin(perceived);
    SetData(angles);

    OBAngle angle;
    OBAtom *b;
    int unique_angle;

    unique_angle = 0;

    FOR_ATOMS_OF_MOL(atom, this) {
      if(atom->GetAtomicNum() == OBElements::Hydrogen)
        continue;

      b = (OBAtom*) &*atom;

      FOR_NBORS_OF_ATOM(a, b) {
        FOR_NBORS_OF_ATOM(c, b) {
          if(&*a == &*c) {
            unique_angle = 1;
            continue;
          }

          if (unique_angle) {
            angle.SetAtoms((OBAtom*)b, (OBAtom*)&*a, (OBAtom*)&*c);
            angles->SetData(angle);
            angle.Clear();
          }
        }
        unique_angle = 0;
      }
    }

    return;
  }

  void OBMol::FindTorsions()
  {
    //if already has data return
    if(HasData(OBGenericDataType::TorsionData))
      return;

    //get new data and attach it to molecule
    OBTorsionData *torsions = new OBTorsionData;
    torsions->SetOrigin(perceived);
    SetData(torsions);

    OBTorsion torsion;
    vector<OBBond*>::iterator bi1,bi2,bi3;
    OBBond* bond;
    OBAtom *a,*b,*c,*d;

    //loop through all bonds generating torsions
    for(bond = BeginBond(bi1);bond;bond = NextBond(bi1))
      {
        b = bond->GetBeginAtom();
        c = bond->GetEndAtom();
        if(b->GetAtomicNum() == OBElements::Hydrogen || c->GetAtomicNum() == OBElements::Hydrogen)
          continue;

        for(a = b->BeginNbrAtom(bi2);a;a = b->NextNbrAtom(bi2))
          {
            if(a == c)
              continue;

            for(d = c->BeginNbrAtom(bi3);d;d = c->NextNbrAtom(bi3))
              {
                if ((d == b) || (d == a))
                  continue;
                torsion.AddTorsion(a,b,c,d);
              }
          }
        //add torsion to torsionData
        if(torsion.GetSize())
          torsions->SetData(torsion);
        torsion.Clear();
      }

    return;
  }

  void OBMol::FindLargestFragment(OBBitVec &lf)
  {
    int j;
    OBAtom *atom;
    OBBond *bond;
    vector<OBAtom*>::iterator i;
    vector<OBBond*>::iterator k;
    OBBitVec used,curr,next,frag;

    lf.Clear();
    while ((unsigned)used.CountBits() < NumAtoms())
      {
        curr.Clear();
        frag.Clear();
        for (atom = BeginAtom(i);atom;atom = NextAtom(i))
          if (!used.BitIsSet(atom->GetIdx()))
            {
              curr.SetBitOn(atom->GetIdx());
              break;
            }

        frag |= curr;
        while (!curr.IsEmpty())
          {
            next.Clear();
            for (j = curr.NextBit(-1);j != curr.EndBit();j = curr.NextBit(j))
              {
                atom = GetAtom(j);
                for (bond = atom->BeginBond(k);bond;bond = atom->NextBond(k))
                  if (!used.BitIsSet(bond->GetNbrAtomIdx(atom)))
                    next.SetBitOn(bond->GetNbrAtomIdx(atom));
              }

            used |= curr;
            used |= next;
            frag |= next;
            curr = next;
          }

        if (lf.IsEmpty() || lf.CountBits() < frag.CountBits())
          lf = frag;
      }
  }

  //! locates all atoms for which there exists a path to 'end'
  //! without going through 'bgn'
  //! children must not include 'end'
  void OBMol::FindChildren(vector<OBAtom*> &children,OBAtom *bgn,OBAtom *end)
  {
    OBBitVec used,curr,next;

    used |= bgn->GetIdx();
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
            atom = GetAtom(i);
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

  //! locates all atoms for which there exists a path to 'second'
  //! without going through 'first'
  //! children must not include 'second'
  void OBMol::FindChildren(vector<int> &children,int first,int second)
  {
    int i;
    OBBitVec used,curr,next;

    used.SetBitOn(first);
    used.SetBitOn(second);
    curr.SetBitOn(second);

    OBAtom *atom;
    while (!curr.IsEmpty())
      {
        next.Clear();
        for (i = curr.NextBit(-1);i != curr.EndBit();i = curr.NextBit(i))
          {
            atom = GetAtom(i);
            FOR_BONDS_OF_ATOM (bond, atom)
              if (!used.BitIsSet(bond->GetNbrAtomIdx(atom)))
                next.SetBitOn(bond->GetNbrAtomIdx(atom));
          }

        used |= next;
        curr = next;
      }

    used.SetBitOff(first);
    used.SetBitOff(second);
    used.ToVecInt(children);
  }

  /*! \brief Calculates the graph theoretical distance (GTD) of each atom.
   *
   * Creates a vector (indexed from zero) containing, for each atom
   * in the molecule, the number of bonds plus one to the most
   * distant non-H atom.
   *
   * For example, for the molecule H3CC(=O)Cl the GTD value for C1
   * would be 3, as the most distant non-H atom (either Cl or O) is
   * 2 bonds away.
   *
   * Since the GTD measures the distance to non-H atoms, the GTD values
   * for terminal H atoms tend to be larger than for non-H terminal atoms.
   * In the example above, the GTD values for the H atoms are all 4.
   */
  bool OBMol::GetGTDVector(vector<int> &gtd)
  //calculates the graph theoretical distance for every atom
  //and puts it into gtd
  {
    gtd.clear();
    gtd.resize(NumAtoms());

    int gtdcount,natom;
    OBBitVec used,curr,next;
    OBAtom *atom,*atom1;
    OBBond *bond;
    vector<OBAtom*>::iterator i;
    vector<OBBond*>::iterator j;

    next.Clear();

    for (atom = BeginAtom(i);atom;atom = NextAtom(i))
      {
        gtdcount = 0;
        used.Clear();
        curr.Clear();
        used.SetBitOn(atom->GetIdx());
        curr.SetBitOn(atom->GetIdx());

        while (!curr.IsEmpty())
          {
            next.Clear();
            for (natom = curr.NextBit(-1);natom != curr.EndBit();natom = curr.NextBit(natom))
              {
                atom1 = GetAtom(natom);
                for (bond = atom1->BeginBond(j);bond;bond = atom1->NextBond(j))
                  if (!used.BitIsSet(bond->GetNbrAtomIdx(atom1)) && !curr.BitIsSet(bond->GetNbrAtomIdx(atom1)))
                    if (bond->GetNbrAtom(atom1)->GetAtomicNum() != OBElements::Hydrogen)
                      next.SetBitOn(bond->GetNbrAtomIdx(atom1));
              }

            used |= next;
            curr = next;
            gtdcount++;
          }
        gtd[atom->GetIdx()-1] = gtdcount;
      }
    return(true);
  }

  /*!
  **\brief Calculates a set of graph invariant indexes using
  ** the graph theoretical distance, number of connected heavy atoms,
  ** aromatic boolean, ring boolean, atomic number, and
  ** summation of bond orders connected to the atom.
  ** Vector is indexed from zero
  */
  void OBMol::GetGIVector(vector<unsigned int> &vid)
  {
    vid.clear();
    vid.resize(NumAtoms()+1);

    vector<int> v;
    GetGTDVector(v);

    int i;
    OBAtom *atom;
    vector<OBAtom*>::iterator j;
    for (i=0,atom = BeginAtom(j);atom;atom = NextAtom(j),++i)
      {
        vid[i] =  (unsigned int)v[i];
        vid[i] += (unsigned int)(atom->GetHvyDegree()*100);
        vid[i] += (unsigned int)(((atom->IsAromatic()) ? 1 : 0)*1000);
        vid[i] += (unsigned int)(((atom->IsInRing()) ? 1 : 0)*10000);
        vid[i] += (unsigned int)(atom->GetAtomicNum()*100000);
        vid[i] += (unsigned int)(atom->GetImplicitHCount()*10000000);
      }
  }

  static bool OBComparePairSecond(const pair<OBAtom*,unsigned int> &a,const pair<OBAtom*,unsigned int> &b)
  {
    return(a.second < b.second);
  }

  static bool OBComparePairFirst(const pair<OBAtom*,unsigned int> &a,const pair<OBAtom*,unsigned int> &b)
  {
    return(a.first->GetIdx() < b.first->GetIdx());
  }

  //! counts the number of unique symmetry classes in a list
  static void ClassCount(vector<pair<OBAtom*,unsigned int> > &vp,unsigned int &count)
  {
    count = 0;
    vector<pair<OBAtom*,unsigned int> >::iterator k;
    sort(vp.begin(),vp.end(),OBComparePairSecond);
#if 0 // original version

    unsigned int id=0; // [ejk] appease gcc's bogus "might be undef'd" warning
    for (k = vp.begin();k != vp.end();++k)
      {
        if (k == vp.begin())
          {
            id = k->second;
            k->second = count = 0;
          }
        else
          if (k->second != id)
            {
              id = k->second;
              k->second = ++count;
            }
          else
            k->second = count;
      }
    count++;
#else // get rid of warning, moves test out of loop, returns 0 for empty input

    k = vp.begin();
    if (k != vp.end())
      {
        unsigned int id = k->second;
        k->second = 0;
        ++k;
        for (;k != vp.end(); ++k)
          {
            if (k->second != id)
              {
                id = k->second;
                k->second = ++count;
              }
            else
              k->second = count;
          }
        ++count;
      }
    else
      {
        // [ejk] thinks count=0 might be OK for an empty list, but orig code did
        //++count;
      }
#endif
  }

  //! creates a new vector of symmetry classes base on an existing vector
  //! helper routine to GetGIDVector
  static void CreateNewClassVector(vector<pair<OBAtom*,unsigned int> > &vp1,vector<pair<OBAtom*,unsigned int> > &vp2)
  {
    int m,id;
    OBAtom *nbr;
    vector<OBBond*>::iterator j;
    vector<unsigned int>::iterator k;
    vector<pair<OBAtom*,unsigned int> >::iterator i;
    sort(vp1.begin(),vp1.end(),OBComparePairFirst);
    vp2.clear();
    for (i = vp1.begin();i != vp1.end();++i)
      {
        vector<unsigned int> vtmp;
        for (nbr = i->first->BeginNbrAtom(j);nbr;nbr = i->first->NextNbrAtom(j))
          vtmp.push_back(vp1[nbr->GetIdx()-1].second);
        sort(vtmp.begin(),vtmp.end(),OBCompareUnsigned);
        for (id=i->second,m=100,k = vtmp.begin();k != vtmp.end();++k,m*=100)
          id += *k * m;

        vp2.push_back(pair<OBAtom*,unsigned int> (i->first,id));
      }
  }

  /*!
  **\brief Calculates a set of symmetry identifiers for a molecule.
  ** Atoms with the same symmetry ID are symmetrically equivalent.
  ** Vector is indexed from zero
  */
  void OBMol::GetGIDVector(vector<unsigned int> &vgid)
  {
    vector<unsigned int> vgi;
    GetGIVector(vgi);  //get vector of graph invariants

    int i;
    OBAtom *atom;
    vector<OBAtom*>::iterator j;
    vector<pair<OBAtom*,unsigned int> > vp1,vp2;
    for (i=0,atom = BeginAtom(j);atom;atom = NextAtom(j),++i)
      vp1.push_back(pair<OBAtom*,unsigned int> (atom,vgi[i]));

    unsigned int nclass1,nclass2; //number of classes
    ClassCount(vp1,nclass1);

    if (nclass1 < NumAtoms())
      {
        for (i = 0;i < 100;++i) //sanity check - shouldn't ever hit this number
          {
            CreateNewClassVector(vp1,vp2);
            ClassCount(vp2,nclass2);
            vp1 = vp2;
            if (nclass1 == nclass2)
              break;
            nclass1 = nclass2;
          }
      }

    vgid.clear();
    sort(vp1.begin(),vp1.end(),OBComparePairFirst);
    vector<pair<OBAtom*,unsigned int> >::iterator k;
    for (k = vp1.begin();k != vp1.end();++k)
      vgid.push_back(k->second);
  }

  unsigned int OBMol::NumHvyAtoms()
  {
    OBAtom *atom;
    vector<OBAtom*>::iterator(i);
    unsigned int count = 0;

    for(atom = this->BeginAtom(i);atom;atom = this->NextAtom(i))
      {
        if (atom->GetAtomicNum() != OBElements::Hydrogen)
          count++;
      }

    return(count);
  }

  unsigned int OBMol::NumRotors(bool sampleRingBonds)
  {
    OBRotorList rl;
    rl.FindRotors(*this, sampleRingBonds);
    return rl.Size();
  }

  //! Returns a pointer to the atom after a safety check
  //! 0 < idx <= NumAtoms
  OBAtom *OBMol::GetAtom(int idx) const
  {
    if ((unsigned)idx < 1 || (unsigned)idx > NumAtoms())
      {
        obErrorLog.ThrowError(__FUNCTION__, "Requested Atom Out of Range", obDebug);
        return((OBAtom*)NULL);
      }

    return((OBAtom*)_vatom[idx-1]);
  }

  OBAtom *OBMol::GetAtomById(unsigned long id) const
  {
    if (id >= _atomIds.size()) {
      obErrorLog.ThrowError(__FUNCTION__, "Requested atom with invalid id.", obDebug);
      return((OBAtom*)NULL);
    }

    return((OBAtom*)_atomIds[id]);
  }

  OBAtom *OBMol::GetFirstAtom() const
  {
    return((_vatom.empty()) ? (OBAtom*)NULL : (OBAtom*)_vatom[0]);
  }

  //! Returns a pointer to the bond after a safety check
  //! 0 <= idx < NumBonds
  OBBond *OBMol::GetBond(int idx) const
  {
    if (idx < 0 || (unsigned)idx >= NumBonds())
      {
        obErrorLog.ThrowError(__FUNCTION__, "Requested Bond Out of Range", obDebug);
        return((OBBond*)NULL);
      }

    return((OBBond*)_vbond[idx]);
  }

  OBBond *OBMol::GetBondById(unsigned long id) const
  {
    if (id >= _bondIds.size()) {
      obErrorLog.ThrowError(__FUNCTION__, "Requested bond with invalid id.", obDebug);
      return((OBBond*)NULL);
    }

    return((OBBond*)_bondIds[id]);
  }

  OBBond *OBMol::GetBond(int bgn, int end) const
  {
    return(GetBond(GetAtom(bgn),GetAtom(end)));
  }

  OBBond *OBMol::GetBond(OBAtom *bgn,OBAtom *end) const
  {
    OBAtom *nbr;
    vector<OBBond*>::iterator i;

    if (!bgn || !end) return NULL;

    for (nbr = bgn->BeginNbrAtom(i);nbr;nbr = bgn->NextNbrAtom(i))
      if (nbr == end)
        return((OBBond *)*i);

    return(NULL); //just to keep the SGI compiler happy
  }

  OBResidue *OBMol::GetResidue(int idx) const
  {
    if (idx < 0 || (unsigned)idx >= NumResidues())
      {
        obErrorLog.ThrowError(__FUNCTION__, "Requested Residue Out of Range", obDebug);
        return((OBResidue*)NULL);
      }

    return (_residue[idx]);
  }

  std::vector<OBInternalCoord*> OBMol::GetInternalCoord()
  {
    if (_internals.empty())
      {
        _internals.push_back((OBInternalCoord*)NULL);
        for(unsigned int i = 1; i <= NumAtoms(); ++i)
          {
            _internals.push_back(new OBInternalCoord);
          }
        CartesianToInternal(_internals, *this);
      }
    return _internals;
  }

  //! Implements <a href="http://qsar.sourceforge.net/dicts/blue-obelisk/index.xhtml#findSmallestSetOfSmallestRings">blue-obelisk:findSmallestSetOfSmallestRings</a>.
  vector<OBRing*> &OBMol::GetSSSR()
  {
    if (!HasSSSRPerceived())
      FindSSSR();

    OBRingData *rd = 0;
    if (!HasData("SSSR")) {
      rd = new OBRingData();
      rd->SetAttribute("SSSR");
      SetData(rd);
    }

    rd = (OBRingData *) GetData("SSSR");
    rd->SetOrigin(perceived);
    return(rd->GetData());
  }

  vector<OBRing*> &OBMol::GetLSSR()
  {
    if (!HasLSSRPerceived())
      FindLSSR();

    OBRingData *rd = 0;
    if (!HasData("LSSR")) {
      rd = new OBRingData();
      rd->SetAttribute("LSSR");
      SetData(rd);
    }

    rd = (OBRingData *) GetData("LSSR");
    rd->SetOrigin(perceived);
    return(rd->GetData());
  }

  double OBMol::GetMolWt(bool implicitH)
  {
    double molwt=0.0;
    OBAtom *atom;
    vector<OBAtom*>::iterator i;

    double hmass = OBElements::GetMass(1);
    for (atom = BeginAtom(i);atom;atom = NextAtom(i)) {
      molwt += atom->GetAtomicMass();
      if (implicitH)
        molwt += atom->GetImplicitHCount() * hmass;
    }
    return(molwt);
  }

  double OBMol::GetExactMass(bool implicitH)
  {
    double mass=0.0;
    OBAtom *atom;
    vector<OBAtom*>::iterator i;

    double hmass = OBElements::GetExactMass(1, 1);
    for (atom = BeginAtom(i); atom; atom = NextAtom(i)) {
      mass += atom->GetExactMass();
      if (implicitH)
        mass += atom->GetImplicitHCount() * hmass;
    }

    return(mass);
  }

  //! Stochoimetric formula in spaced format e.g. C 4 H 6 O 1
  //! No pair data is stored. Normally use without parameters: GetSpacedFormula()
  //! \since version 2.1
  string OBMol::GetSpacedFormula(int ones, const char* sp, bool implicitH)
  {
    //Default ones=0, sp=" ".
    //Using ones=1 and sp="" will give unspaced formula (and no pair data entry)
    // These are the atomic numbers of the elements in alphabetical order, plus
    // pseudo atomic numbers for D, T isotopes.
    const int NumElements = 118 + 2;
    const int alphabetical[NumElements] = {
      89, 47, 13, 95, 18, 33, 85, 79, 5, 56, 4, 107, 83, 97, 35, 6, 20, 48,
      58, 98, 17, 96, 112, 27, 24, 55, 29, NumElements-1,
      105, 110, 66, 68, 99, 63, 9, 26, 114, 100, 87, 31,
      64, 32, 1, 2, 72, 80, 67, 108, 53, 49, 77, 19, 36, 57, 3, 103, 71, 116, 115, 101,
      12, 25, 42, 109, 7, 11, 41, 60, 10, 113, 28, 102, 93, 8, 118, 76, 15, 91, 82, 46,
      61, 84, 59, 78, 94, 88, 37, 75, 104, 111, 45, 86, 44, 16, 51, 21, 34, 106, 14,
      62, 50, 38, NumElements, 73, 65, 43, 52, 90, 22, 81, 69, 117, 92, 23, 74, 54, 39, 70,
      30, 40 };

    int atomicCount[NumElements];
    stringstream formula;

    for (int i = 0; i < NumElements; ++i)
      atomicCount[i] = 0;

    bool UseImplicitH = (NumBonds()!=0 || NumAtoms()==1);
    // Do not use implicit hydrogens if explicitly required not to
    if (!implicitH) UseImplicitH = false;
    bool HasHvyAtoms = NumHvyAtoms()>0;
    FOR_ATOMS_OF_MOL(a, *this)
      {
        int anum = a->GetAtomicNum();
        if(anum==0)
          continue;
        if(anum > (NumElements-2)) {
          char buffer[BUFF_SIZE];  // error buffer
          snprintf(buffer, BUFF_SIZE, "Skipping unknown element with atomic number %d", anum);
          obErrorLog.ThrowError(__FUNCTION__, buffer, obWarning);
          continue;
        }
        bool IsHiso = anum == 1 && a->GetIsotope()>=2;
        if(UseImplicitH)
          {
            if (anum == 1 && !IsHiso && HasHvyAtoms)
              continue; // skip explicit hydrogens except D,T
            if(anum==1)
              {
                if (IsHiso && HasHvyAtoms)
                  --atomicCount[0]; //one of the implicit hydrogens is now explicit
              }
            else
              atomicCount[0] += a->GetImplicitHCount() + a->ExplicitHydrogenCount();
          }
        if (IsHiso)
          anum = NumElements + a->GetIsotope() - 3; //pseudo AtNo for D, T
        atomicCount[anum - 1]++;
      }

    if (atomicCount[5] != 0) // Carbon (i.e. 6 - 1 = 5)
      {
        if (atomicCount[5] > ones)
          formula << "C" << sp << atomicCount[5] << sp;
        else if (atomicCount[5] == 1)
          formula << "C";

        atomicCount[5] = 0; // So we don't output C twice

        // only output H if there's also carbon -- otherwise do it alphabetical
        if (atomicCount[0] != 0) // Hydrogen (i.e., 1 - 1 = 0)
          {
            if (atomicCount[0] > ones)
              formula << "H" << sp << atomicCount[0] << sp;
            else if (atomicCount[0] == 1)
              formula << "H";

            atomicCount[0] = 0;
          }
      }

    for (int j = 0; j < NumElements; ++j)
      {
        char DT[4] = {'D',0,'T',0};
        const char* symb;
        int alph = alphabetical[j]-1;
        if (atomicCount[ alph ])
          {
            if(alph==NumElements-1)
              symb = DT + 2;//T
            else if (alph==NumElements-2)
              symb = DT; //D
            else
              symb = OBElements::GetSymbol(alphabetical[j]);

            formula << symb << sp;
            if(atomicCount[alph] > ones)
              formula << atomicCount[alph] << sp;
          }
      }

    int chg = GetTotalCharge();
    char ch = chg>0 ? '+' : '-' ;
    chg = abs(chg);
    while(chg--)
      formula << ch << sp;

    string f_str = formula.str();
    return (Trim(f_str));
  }

  //! Stochoimetric formula (e.g., C4H6O).
  //!   This is either set by OBMol::SetFormula() or generated on-the-fly
  //!   using the "Hill order" -- i.e., C first if present, then H if present
  //!   all other elements in alphabetical order.
  string OBMol::GetFormula()
  {
    string attr = "Formula";
    OBPairData *dp = (OBPairData *) GetData(attr);

    if (dp != NULL) // we already set the formula (or it was read from a file)
      return dp->GetValue();

    obErrorLog.ThrowError(__FUNCTION__,
                          "Ran OpenBabel::SetFormula -- Hill order formula",
                          obAuditMsg);

    string sformula = GetSpacedFormula(1, "");

    dp = new OBPairData;
    dp->SetAttribute(attr);
    dp->SetValue( sformula );
    dp->SetOrigin( perceived ); // internal generation
    SetData(dp);
    return sformula;
  }

  void OBMol::SetFormula(string molFormula)
  {
    string attr = "Formula";
    OBPairData *dp = (OBPairData *) GetData(attr);
    if (dp == NULL)
      {
        dp = new OBPairData;
        dp->SetAttribute(attr);
        SetData(dp);
      }
    dp->SetValue(molFormula);
    // typically file input, but this needs to be revisited
    dp->SetOrigin(fileformatInput);
  }

  void OBMol::SetTotalCharge(int charge)
  {
    SetFlag(OB_TCHARGE_MOL);
    _totalCharge = charge;
  }

  //! Returns the total molecular charge -- if it has not previously been set
  //!  it is calculated from the atomic formal charge information.
  //!  (This may or may not be correct!)
  //!  If you set atomic charges with OBAtom::SetFormalCharge()
  //!   you really should set the molecular charge with OBMol::SetTotalCharge()
  int OBMol::GetTotalCharge()
  {
    if(HasFlag(OB_TCHARGE_MOL))
      return(_totalCharge);
    else // calculate from atomic formal charges (seems the best default)
      {
        obErrorLog.ThrowError(__FUNCTION__,
                              "Ran OpenBabel::GetTotalCharge -- calculated from formal charges",
                              obAuditMsg);

        OBAtom *atom;
        vector<OBAtom*>::iterator i;
        int chg = 0;

        for (atom = BeginAtom(i);atom;atom = NextAtom(i))
          chg += atom->GetFormalCharge();
        return (chg);
      }
  }

  void   OBMol::SetTotalSpinMultiplicity(unsigned int spin)
  {
    SetFlag(OB_TSPIN_MOL);
    _totalSpin = spin;
  }

  void OBMol::SetInternalCoord(std::vector<OBInternalCoord*> int_coord) {
    if (int_coord[0] != NULL) {
      std::vector<OBInternalCoord*>::iterator it = int_coord.begin();
      int_coord.insert(it, static_cast<OBInternalCoord*>(NULL));
    }

    if (int_coord.size() != _natoms + 1) {
      string error = "Number of internal coordinates is not the same as";
      error += " the number of atoms in molecule";
      obErrorLog.ThrowError(__FUNCTION__, error, obError);
      return;
    }

    _internals = int_coord;

    return;
  }

  //! Returns the total spin multiplicity -- if it has not previously been set
  //!  It is calculated from the atomic spin multiplicity information
  //!  assuming the high-spin case (i.e. it simply sums the number of unpaired
  //!  electrons assuming no further pairing of spins.
  //!  if it fails (gives singlet for odd number of electronic systems),
  //!  then assign wrt parity of the total electrons.
  unsigned int OBMol::GetTotalSpinMultiplicity()
  {
    if (HasFlag(OB_TSPIN_MOL))
      return(_totalSpin);
    else // calculate from atomic spin information (assuming high-spin case)
      {
        obErrorLog.ThrowError(__FUNCTION__,
                              "Ran OpenBabel::GetTotalSpinMultiplicity -- calculating from atomic spins assuming high spin case",
                              obAuditMsg);

        OBAtom *atom;
        vector<OBAtom*>::iterator i;
        unsigned int unpaired_electrons = 0;
        int chg = GetTotalCharge();
        for (atom = BeginAtom(i);atom;atom = NextAtom(i))
          {
            if (atom->GetSpinMultiplicity() > 1)
              unpaired_electrons += (atom->GetSpinMultiplicity() - 1);
           chg += atom->GetAtomicNum();
          }
        if (chg % 2 != unpaired_electrons %2)
          return ((abs(chg) % 2) + 1);
        else
          return (unpaired_electrons + 1);
      }
  }

  OBMol &OBMol::operator=(const OBMol &source)
  //atom and bond info is copied from src to dest
  //Conformers are now copied also, MM 2/7/01
  //Residue information are copied, MM 4-27-01
  //All OBGenericData incl OBRotameterList is copied, CM 2006
  //Zeros all flags except OB_TCHARGE_MOL, OB_PCHARGE_MOL, OB_HYBRID_MOL
  //OB_TSPIN_MOL, OB_AROMATIC_MOL, OB_CHAINS_MOL and OB_PATTERN_STRUCTURE which are copied
  {
    if (this == &source)
      return *this;

    OBMol &src = (OBMol &)source;
    vector<OBAtom*>::iterator i;
    vector<OBBond*>::iterator j;
    OBAtom *atom;
    OBBond *bond;

    Clear();
    BeginModify();

    _vatom.reserve(src.NumAtoms());
    _atomIds.reserve(src.NumAtoms());
    _vbond.reserve(src.NumBonds());
    _bondIds.reserve(src.NumBonds());

    for (atom = src.BeginAtom(i);atom;atom = src.NextAtom(i))
      AddAtom(*atom);
    for (bond = src.BeginBond(j);bond;bond = src.NextBond(j))
      AddBond(*bond);

    this->_title  = src.GetTitle();
    this->_energy = src.GetEnergy();
    this->_dimension = src.GetDimension();
    this->SetTotalCharge(src.GetTotalCharge()); //also sets a flag
    this->SetTotalSpinMultiplicity(src.GetTotalSpinMultiplicity()); //also sets a flag

    EndModify(); //zeros flags!

    if (src.HasFlag(OB_PATTERN_STRUCTURE))
      this->SetFlag(OB_PATTERN_STRUCTURE);
    if (src.HasFlag(OB_TSPIN_MOL))
      this->SetFlag(OB_TSPIN_MOL);
    if (src.HasFlag(OB_TCHARGE_MOL))
      this->SetFlag(OB_TCHARGE_MOL);
    if (src.HasFlag(OB_PCHARGE_MOL))
      this->SetFlag(OB_PCHARGE_MOL);
    if (src.HasFlag(OB_HYBRID_MOL))
      this->SetFlag(OB_HYBRID_MOL);
    if (src.HasFlag(OB_AROMATIC_MOL))
      this->SetFlag(OB_AROMATIC_MOL);
    if (src.HasFlag(OB_CHAINS_MOL))
      this->SetFlag(OB_CHAINS_MOL);
    //this->_flags = src.GetFlags(); //Copy all flags. Perhaps too drastic a change


    //Copy Residue information
    unsigned int NumRes = src.NumResidues();
    if (NumRes)
      {
        unsigned int k;
        OBResidue *src_res=NULL;
        OBResidue *res=NULL;
        OBAtom *src_atom=NULL;
        OBAtom *atom=NULL;
        vector<OBAtom*>::iterator ii;
        for (k=0 ; k<NumRes ; ++k)
          {
            res = NewResidue();
            src_res = src.GetResidue(k);
            *res = *src_res; //does not copy atoms
            for (src_atom=src_res->BeginAtom(ii) ; src_atom ; src_atom=src_res->NextAtom(ii))
              {
                atom = GetAtom(src_atom->GetIdx());
                res->AddAtom(atom);
                res->SetAtomID(atom,src_res->GetAtomID(src_atom));
                res->SetHetAtom(atom,src_res->IsHetAtom(src_atom));
                res->SetSerialNum(atom,src_res->GetSerialNum(src_atom));
              }
          }
      }

    //Copy conformer information
    if (src.NumConformers() > 1) {
      int k;//,l;
      vector<double*> conf;
      int currConf = -1;
      double* xyz = NULL;
      for (k=0 ; k<src.NumConformers() ; ++k) {
        xyz = new double [3*src.NumAtoms()];
        memcpy( xyz, src.GetConformer(k), sizeof( double )*3*src.NumAtoms() );
        conf.push_back(xyz);

        if( src.GetConformer(k) == src._c ) {
          currConf = k;
        }
      }

      SetConformers(conf);
      if( currConf >= 0 && _vconf.size() ) {
        _c = _vconf[currConf];
      }
    }

    //Copy all the OBGenericData, providing the new molecule, this,
    //for those classes like OBRotameterList which contain Atom pointers
    //OBGenericData classes can choose not to be cloned by returning NULL
    vector<OBGenericData*>::iterator itr;
    for(itr=src.BeginData();itr!=src.EndData();++itr)
      {
        OBGenericData* pCopiedData = (*itr)->Clone(this);
        SetData(pCopiedData);
      }

    if (src.HasChiralityPerceived())
      SetChiralityPerceived();

    return(*this);
  }

  OBMol &OBMol::operator+=(const OBMol &source)
  {
    OBMol &src = (OBMol &)source;
    vector<OBAtom*>::iterator i;
    vector<OBBond*>::iterator j;
    vector<OBResidue*>::iterator k;
    OBAtom *atom;
    OBBond *bond;
    OBResidue *residue;

    BeginModify();

    int prevatms = NumAtoms();

    string extitle(src.GetTitle());
    if(!extitle.empty())
      _title += "_" + extitle;

    // First, handle atoms and bonds
    map<unsigned long int, unsigned long int> correspondingId;
    for (atom = src.BeginAtom(i) ; atom ; atom = src.NextAtom(i)) {
      AddAtom(*atom, true); // forceNewId=true (don't reuse the original Id)
      OBAtom *addedAtom = GetAtom(NumAtoms());
      correspondingId[atom->GetId()] = addedAtom->GetId();
    }
    correspondingId[OBStereo::ImplicitRef] = OBStereo::ImplicitRef;

    for (bond = src.BeginBond(j) ; bond ; bond = src.NextBond(j)) {
      bond->SetId(NoId);//Need to remove ID which relates to source mol rather than this mol
      AddBond(bond->GetBeginAtomIdx() + prevatms,
              bond->GetEndAtomIdx() + prevatms,
              bond->GetBondOrder(), bond->GetFlags());
    }

    // Now update all copied residues too
    for (residue = src.BeginResidue(k); residue; residue = src.NextResidue(k)) {
      AddResidue(*residue);

      FOR_ATOMS_OF_RESIDUE(resAtom, residue)
        {
          // This is the equivalent atom in our combined molecule
          atom = GetAtom(resAtom->GetIdx() + prevatms);
          // So we add this to the last-added residue
          // (i.e., what we just copied)
          (_residue[_residue.size() - 1])->AddAtom(atom);
        }
    }

    // Copy the stereo
    std::vector<OBGenericData*> vdata = src.GetAllData(OBGenericDataType::StereoData);
    for (std::vector<OBGenericData*>::iterator data = vdata.begin(); data != vdata.end(); ++data) {
      OBStereo::Type datatype = ((OBStereoBase*)*data)->GetType();
      if (datatype == OBStereo::CisTrans) {
        OBCisTransStereo *ct = dynamic_cast<OBCisTransStereo*>(*data);
        OBCisTransStereo *nct = new OBCisTransStereo(this);
        OBCisTransStereo::Config ct_cfg = ct->GetConfig();
        ct_cfg.begin = correspondingId[ct_cfg.begin];
        ct_cfg.end = correspondingId[ct_cfg.end];
        for(OBStereo::RefIter ri = ct_cfg.refs.begin(); ri != ct_cfg.refs.end(); ++ri)
          *ri = correspondingId[*ri];
        nct->SetConfig(ct_cfg);
        SetData(nct);
      }
      else if (datatype == OBStereo::Tetrahedral) {
        OBTetrahedralStereo *ts = dynamic_cast<OBTetrahedralStereo*>(*data);
        OBTetrahedralStereo *nts = new OBTetrahedralStereo(this);
        OBTetrahedralStereo::Config ts_cfg = ts->GetConfig();
        ts_cfg.center = correspondingId[ts_cfg.center];
        ts_cfg.from = correspondingId[ts_cfg.from];
        for(OBStereo::RefIter ri = ts_cfg.refs.begin(); ri != ts_cfg.refs.end(); ++ri)
          *ri = correspondingId[*ri];
        nts->SetConfig(ts_cfg);
        SetData(nts);
      }
    }

    // TODO: This is actually a weird situation (e.g., adding a 2D mol to 3D one)
    // We should do something to update the src coordinates if they're not 3D
    if(src.GetDimension()<_dimension)
      _dimension = src.GetDimension();

    EndModify();

    return(*this);
  }

  bool OBMol::Clear()
  {
    if (obErrorLog.GetOutputLevel() >= obAuditMsg)
      obErrorLog.ThrowError(__FUNCTION__,
                            "Ran OpenBabel::Clear Molecule", obAuditMsg);

    vector<OBAtom*>::iterator i;
    vector<OBBond*>::iterator j;
    for (i = _vatom.begin();i != _vatom.end();++i)
      {
        DestroyAtom(*i);
        *i = NULL;
      }
    for (j = _vbond.begin();j != _vbond.end();++j)
      {
        DestroyBond(*j);
        *j = NULL;
      }

    _atomIds.clear();
    _bondIds.clear();
    _natoms = _nbonds = 0;

    //Delete residues
    unsigned int ii;
    for (ii=0 ; ii<_residue.size() ; ++ii)
      {
        DestroyResidue(_residue[ii]);
      }
    _residue.clear();

    //clear out the multiconformer data
    vector<double*>::iterator k;
    for (k = _vconf.begin();k != _vconf.end();++k)
      delete [] *k;
    _vconf.clear();

    //Clear flags except OB_PATTERN_STRUCTURE which is left the same
    _flags &= OB_PATTERN_STRUCTURE;

    _c = (double*) NULL;
    _mod = 0;

    // Clean up generic data via the base class
    return(OBBase::Clear());
  }

  void OBMol::BeginModify()
  {
    //suck coordinates from _c into _v for each atom
    if (!_mod && !Empty())
      {
        OBAtom *atom;
        vector<OBAtom*>::iterator i;
        for (atom = BeginAtom(i);atom;atom = NextAtom(i))
          {
            atom->SetVector();
            atom->ClearCoordPtr();
          }

        vector<double*>::iterator j;
        for (j = _vconf.begin();j != _vconf.end();++j)
          delete [] *j;

        _c = NULL;
        _vconf.clear();

        //Destroy rotamer list if necessary
        if ((OBRotamerList *)GetData(OBGenericDataType::RotamerList))
          {
            delete (OBRotamerList *)GetData(OBGenericDataType::RotamerList);
            DeleteData(OBGenericDataType::RotamerList);
          }
      }

    _mod++;
  }

  void OBMol::EndModify(bool nukePerceivedData)
  {
    if (_mod == 0)
      {
        obErrorLog.ThrowError(__FUNCTION__, "_mod is negative - EndModify() called too many times", obDebug);
        return;
      }

    _mod--;

    if (_mod)
      return;

    // wipe all but whether it has aromaticity perceived or is a reaction
    if (nukePerceivedData)
      _flags = _flags & (OB_AROMATIC_MOL|OB_REACTION_MOL);

    _c = NULL;

    if (Empty())
      return;

    //if atoms present convert coords into array
    double *c = new double [NumAtoms()*3];
    _c = c;

    unsigned int idx;
    OBAtom *atom;
    vector<OBAtom*>::iterator j;
    for (idx=0,atom = BeginAtom(j);atom;atom = NextAtom(j),++idx)
      {
        atom->SetIdx(idx+1);
        (atom->GetVector()).Get(&_c[idx*3]);
        atom->SetCoordPtr(&_c);
      }
    _vconf.push_back(c);

    // Always remove angle and torsion data, since they will interfere with the iterators
    // PR#2812013
    DeleteData(OBGenericDataType::AngleData);
    DeleteData(OBGenericDataType::TorsionData);
  }

  void OBMol::DestroyAtom(OBAtom *atom)
  {
    if (atom)
      {
        delete atom;
        atom = NULL;
      }
  }

  void OBMol::DestroyBond(OBBond *bond)
  {
    if (bond)
      {
        delete bond;
        bond = NULL;
      }
  }

  void OBMol::DestroyResidue(OBResidue *residue)
  {
    if (residue)
      {
        delete residue;
        residue = NULL;
      }
  }

  OBAtom *OBMol::NewAtom()
  {
    return NewAtom(_atomIds.size());
  }

  //! \brief Instantiate a New Atom and add it to the molecule
  //!
  //! Checks bond_queue for any bonds that should be made to the new atom
  //! and updates atom indexes.
  OBAtom *OBMol::NewAtom(unsigned long id)
  {
    //   BeginModify();

    // resize _atomIds if needed
    if (id >= _atomIds.size()) {
      unsigned int size = _atomIds.size();
      _atomIds.resize(id+1);
      for (unsigned long i = size; i < id; ++i)
        _atomIds[i] = (OBAtom*)NULL;
    }

    if (_atomIds.at(id))
      return (OBAtom*)NULL;

    OBAtom *obatom = new OBAtom;
    obatom->SetIdx(_natoms+1);
    obatom->SetParent(this);

    _atomIds[id] = obatom;
    obatom->SetId(id);

#define OBAtomIncrement 100

    if (_natoms+1 >= _vatom.size())
      {
        _vatom.resize(_natoms+OBAtomIncrement);
        vector<OBAtom*>::iterator j;
        for (j = _vatom.begin(),j+=(_natoms+1);j != _vatom.end();++j)
          *j = (OBAtom*)NULL;
      }
#undef OBAtomIncrement


    _vatom[_natoms] = obatom;
    _natoms++;

    if (HasData(OBGenericDataType::VirtualBondData))
      {
        /*add bonds that have been queued*/
        OBVirtualBond *vb;
        vector<OBGenericData*> verase;
        vector<OBGenericData*>::iterator i;
        for (i = BeginData();i != EndData();++i)
          if ((*i)->GetDataType() == OBGenericDataType::VirtualBondData)
            {
              vb = (OBVirtualBond*)*i;
              if (vb->GetBgn() > _natoms || vb->GetEnd() > _natoms)
                continue;
              if (obatom->GetIdx() == static_cast<unsigned int>(vb->GetBgn())
                  || obatom->GetIdx() == static_cast<unsigned int>(vb->GetEnd()))
                {
                  AddBond(vb->GetBgn(),vb->GetEnd(),vb->GetOrder());
                  verase.push_back(*i);
                }
            }

        if (!verase.empty())
          DeleteData(verase);
      }

    // EndModify();

    return(obatom);
  }

  OBResidue *OBMol::NewResidue()
  {
    OBResidue *obresidue = new OBResidue;
    obresidue->SetIdx(_residue.size());
    _residue.push_back(obresidue);
    return(obresidue);
  }

  OBBond *OBMol::NewBond()
  {
    return NewBond(_bondIds.size());
  }

  //! \since version 2.1
  //! \brief Instantiate a New Bond and add it to the molecule
  //!
  //! Sets the proper Bond index and insures this molecule is set as the parent.
  OBBond *OBMol::NewBond(unsigned long id)
  {
    // resize _bondIds if needed
    if (id >= _bondIds.size()) {
      unsigned int size = _bondIds.size();
      _bondIds.resize(id+1);
      for (unsigned long i = size; i < id; ++i)
        _bondIds[i] = (OBBond*)NULL;
    }

    if (_bondIds.at(id))
      return (OBBond*)NULL;

    OBBond *pBond = new OBBond;
    pBond->SetParent(this);
    pBond->SetIdx(_nbonds);

    _bondIds[id] = pBond;
    pBond->SetId(id);

#define OBBondIncrement 100
    if (_nbonds+1 >= _vbond.size())
      {
        _vbond.resize(_nbonds+OBBondIncrement);
        vector<OBBond*>::iterator i;
        for (i = _vbond.begin(),i+=(_nbonds+1);i != _vbond.end();++i)
          *i = (OBBond*)NULL;
      }
#undef  OBBondIncrement

    _vbond[_nbonds] = (OBBond*)pBond;
    _nbonds++;

    return(pBond);
  }

  //! \brief Add an atom to a molecule
  //!
  //! Also checks bond_queue for any bonds that should be made to the new atom
  bool OBMol::AddAtom(OBAtom &atom, bool forceNewId)
  {
    //    BeginModify();

    // Use the existing atom Id unless either it's invalid or forceNewId has been specified
    unsigned long id;
    if (forceNewId)
      id = _atomIds.size();
    else {
      id = atom.GetId();
      if (id == NoId)
        id = _atomIds.size();
    }

    OBAtom *obatom = new OBAtom;
    *obatom = atom;
    obatom->SetIdx(_natoms+1);
    obatom->SetParent(this);

    // resize _atomIds if needed
    if (id >= _atomIds.size()) {
      unsigned int size = _atomIds.size();
      _atomIds.resize(id+1);
      for (unsigned long i = size; i < id; ++i)
        _atomIds[i] = (OBAtom*)NULL;
    }

    obatom->SetId(id);
    _atomIds[id] = obatom;

#define OBAtomIncrement 100

    if (_natoms+1 >= _vatom.size())
      {
        _vatom.resize(_natoms+OBAtomIncrement);
        vector<OBAtom*>::iterator j;
        for (j = _vatom.begin(),j+=(_natoms+1);j != _vatom.end();++j)
          *j = (OBAtom*)NULL;
      }
#undef OBAtomIncrement

    _vatom[_natoms] = (OBAtom*)obatom;
    _natoms++;

    if (HasData(OBGenericDataType::VirtualBondData))
      {
        /*add bonds that have been queued*/
        OBVirtualBond *vb;
        vector<OBGenericData*> verase;
        vector<OBGenericData*>::iterator i;
        for (i = BeginData();i != EndData();++i)
          if ((*i)->GetDataType() == OBGenericDataType::VirtualBondData)
            {
              vb = (OBVirtualBond*)*i;
              if (vb->GetBgn() > _natoms || vb->GetEnd() > _natoms)
                continue;
              if (obatom->GetIdx() == static_cast<unsigned int>(vb->GetBgn())
                  || obatom->GetIdx() == static_cast<unsigned int>(vb->GetEnd()))
                {
                  AddBond(vb->GetBgn(),vb->GetEnd(),vb->GetOrder());
                  verase.push_back(*i);
                }
            }

        if (!verase.empty())
          DeleteData(verase);
      }

    //    EndModify();

    return(true);
  }

  bool OBMol::InsertAtom(OBAtom &atom)
  {
    BeginModify();
    AddAtom(atom);
    EndModify();

    return(true);
  }

  bool OBMol::AddResidue(OBResidue &residue)
  {
    BeginModify();

    OBResidue *obresidue = new OBResidue;
    *obresidue = residue;

    obresidue->SetIdx(_residue.size());

    _residue.push_back(obresidue);

    EndModify();

    return(true);
  }

  bool OBMol::StripSalts(unsigned int threshold)
  {
    vector<vector<int> > cfl;
    vector<vector<int> >::iterator i,max;

    ContigFragList(cfl);
    if (cfl.empty() || cfl.size() == 1)
      {
        return(false);
      }


    obErrorLog.ThrowError(__FUNCTION__, "Ran OpenBabel::StripSalts", obAuditMsg);

    max = cfl.begin();
    for (i = cfl.begin();i != cfl.end();++i)
      {
        if ((*max).size() < (*i).size())
          max = i;
      }

    vector<int>::iterator j;
    vector<OBAtom*> delatoms;
    set<int> atomIndices;
    for (i = cfl.begin(); i != cfl.end(); ++i)
      {
        if (i->size() < threshold || (threshold == 0 && i != max))
          {
            for (j = (*i).begin(); j != (*i).end(); ++j)
              {
                if (atomIndices.find( *j ) == atomIndices.end())
                  {
                    delatoms.push_back(GetAtom(*j));
                    atomIndices.insert(*j);
                  }
              }
          }
      }

    if (!delatoms.empty())
      {
        //      int tmpflags = _flags & (~(OB_SSSR_MOL));
        BeginModify();
        vector<OBAtom*>::iterator k;
        for (k = delatoms.begin(); k != delatoms.end(); ++k)
          DeleteAtom((OBAtom*)*k);
        EndModify();
        //      _flags = tmpflags;  // Gave crash when SmartsPattern::Match()
        // was called susequently
        // Hans De Winter; 03-nov-2010
      }

    return (true);
  }

  // Convenience function used by the DeleteHydrogens methods
  static bool IsSuppressibleHydrogen(OBAtom *atom)
  {
    if (atom->GetIsotope() == 0 && atom->GetHvyDegree() == 1 && atom->GetFormalCharge() == 0
        && !atom->GetData("Atom Class"))
      return true;
    else
      return false;
  }

  bool OBMol::DeletePolarHydrogens()
  {
    OBAtom *atom;
    vector<OBAtom*>::iterator i;
    vector<OBAtom*> delatoms;

    obErrorLog.ThrowError(__FUNCTION__,
                          "Ran OpenBabel::DeleteHydrogens -- polar",
                          obAuditMsg);

    for (atom = BeginAtom(i);atom;atom = NextAtom(i))
      if (atom->IsPolarHydrogen() && IsSuppressibleHydrogen(atom))
        delatoms.push_back(atom);

    if (delatoms.empty())
      return(true);

    IncrementMod();

    for (i = delatoms.begin();i != delatoms.end();++i)
      DeleteAtom((OBAtom *)*i);

    DecrementMod();

    SetSSSRPerceived(false);
    SetLSSRPerceived(false);
    return(true);
  }


  bool OBMol::DeleteNonPolarHydrogens()
  {
    OBAtom *atom;
    vector<OBAtom*>::iterator i;
    vector<OBAtom*> delatoms;

    obErrorLog.ThrowError(__FUNCTION__,
                          "Ran OpenBabel::DeleteHydrogens -- nonpolar",
                          obAuditMsg);


    for (atom = BeginAtom(i);atom;atom = NextAtom(i))
      if (atom->IsNonPolarHydrogen() && IsSuppressibleHydrogen(atom))
        delatoms.push_back(atom);

    if (delatoms.empty())
      return(true);

    /*
      int idx1,idx2;
      vector<double*>::iterator j;
      for (idx1=0,idx2=0,atom = BeginAtom(i);atom;atom = NextAtom(i),++idx1)
      if (atom->GetAtomicNum() != OBElements::Hydrogen)
      {
      for (j = _vconf.begin();j != _vconf.end();++j)
      memcpy((char*)&((*j)[idx2*3]),(char*)&((*j)[idx1*3]),sizeof(double)*3);
      idx2++;
      }
    */

    IncrementMod();

    for (i = delatoms.begin();i != delatoms.end();++i)
      DeleteAtom((OBAtom *)*i);

    DecrementMod();

    SetSSSRPerceived(false);
    SetLSSRPerceived(false);
    return(true);
  }

  bool OBMol::DeleteHydrogens()
  {
    OBAtom *atom;//,*nbr;
    vector<OBAtom*>::iterator i;
    vector<OBAtom*> delatoms,va;

    obErrorLog.ThrowError(__FUNCTION__,
                          "Ran OpenBabel::DeleteHydrogens", obAuditMsg);

    for (atom = BeginAtom(i);atom;atom = NextAtom(i))
      if (atom->GetAtomicNum() == OBElements::Hydrogen && IsSuppressibleHydrogen(atom))
        delatoms.push_back(atom);

    SetHydrogensAdded(false);

    if (delatoms.empty())
      return(true);

    /* decide whether these flags need to be reset
       _flags &= (~(OB_ATOMTYPES_MOL));
       _flags &= (~(OB_HYBRID_MOL));
       _flags &= (~(OB_PCHARGE_MOL));
       _flags &= (~(OB_IMPVAL_MOL));
    */

    IncrementMod();

    // This is slow -- we need methods to delete a set of atoms
    //  and to delete a set of bonds
    // Calling this sequentially does result in correct behavior
    //  (e.g., fixing PR# 1704551)
    OBBondIterator bi;
    for (i = delatoms.begin(); i != delatoms.end(); ++i) {
      OBAtom* nbr = (*i)->BeginNbrAtom(bi);
      if (nbr) // defensive
        nbr->SetImplicitHCount(nbr->GetImplicitHCount() + 1);
      DeleteAtom((OBAtom *)*i);
    }

    DecrementMod();

    SetSSSRPerceived(false);
    SetLSSRPerceived(false);
    return(true);
  }

  bool OBMol::DeleteHydrogens(OBAtom *atom)
  //deletes all hydrogens attached to the atom passed to the function
  {
    OBAtom *nbr;
    vector<OBAtom*>::iterator i;
    vector<OBBond*>::iterator k;
    vector<OBAtom*> delatoms;

    for (nbr = atom->BeginNbrAtom(k);nbr;nbr = atom->NextNbrAtom(k))
      if (nbr->GetAtomicNum() == OBElements::Hydrogen && IsSuppressibleHydrogen(atom))
        delatoms.push_back(nbr);

    if (delatoms.empty())
      return(true);

    IncrementMod();
    for (i = delatoms.begin();i != delatoms.end();++i)
      DeleteHydrogen((OBAtom*)*i);
    DecrementMod();

    SetHydrogensAdded(false);
    SetSSSRPerceived(false);
    SetLSSRPerceived(false);
    return(true);
  }

  bool OBMol::DeleteHydrogen(OBAtom *atom)
  //deletes the hydrogen atom passed to the function
  {
    if (atom->GetAtomicNum() != OBElements::Hydrogen)
      return false;

    unsigned atomidx = atom->GetIdx();

    //find bonds to delete
    OBAtom *nbr;
    vector<OBBond*> vdb;
    vector<OBBond*>::iterator j;
    for (nbr = atom->BeginNbrAtom(j);nbr;nbr = atom->NextNbrAtom(j))
      vdb.push_back(*j);

    IncrementMod();
    for (j = vdb.begin();j != vdb.end();++j)
      DeleteBond((OBBond*)*j); //delete bonds
    DecrementMod();

    int idx;
    if (atomidx != NumAtoms())
      {
        idx = atom->GetCoordinateIdx();
        int size = NumAtoms()-atom->GetIdx();
        vector<double*>::iterator k;
        for (k = _vconf.begin();k != _vconf.end();++k)
          memmove((char*)&(*k)[idx],(char*)&(*k)[idx+3],sizeof(double)*3*size);

      }

    // Deleting hydrogens does not invalidate the stereo objects
    // - however, any explicit refs to the hydrogen atom must be
    //   converted to implicit refs
    OBStereo::Ref id = atom->GetId();
    StereoRefToImplicit(*this, id);

    _atomIds[id] = (OBAtom*)NULL;
    _vatom.erase(_vatom.begin()+(atomidx-1));
    _natoms--;

    //reset all the indices to the atoms
    vector<OBAtom*>::iterator i;
    OBAtom *atomi;
    for (idx=1,atomi = BeginAtom(i);atomi;atomi = NextAtom(i),++idx)
      atomi->SetIdx(idx);

    SetHydrogensAdded(false);

    DestroyAtom(atom);

    SetSSSRPerceived(false);
    SetLSSRPerceived(false);
    return(true);
  }

  /*
  this has become a wrapper for backward compatibility
  */
  bool OBMol::AddHydrogens(bool polaronly, bool correctForPH, double pH)
  {
    return(AddNewHydrogens(polaronly ? PolarHydrogen : AllHydrogen, correctForPH, pH));
  }

  static bool AtomIsNSOP(OBAtom *atom)
  {
    switch (atom->GetAtomicNum()) {
    case OBElements::Nitrogen:
    case OBElements::Sulfur:
    case OBElements::Oxygen:
    case OBElements::Phosphorus:
      return true;
    default:
      return false;
    }
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

  bool OBMol::AddNewHydrogens(HydrogenType whichHydrogen, bool correctForPH, double pH)
  {
    if (!IsCorrectedForPH() && correctForPH)
      CorrectForPH(pH);

    if (HasHydrogensAdded())
      return(true);

    bool hasChiralityPerceived = this->HasChiralityPerceived(); // remember

    /*
    //
    // This was causing bug #1892844 in avogadro. We also want to add hydrogens if the molecule has no bonds.
    //
    if(NumBonds()==0 && NumAtoms()!=1)
    {
    obErrorLog.ThrowError(__FUNCTION__,
    "Did not run OpenBabel::AddHydrogens on molecule with no bonds", obAuditMsg);
    return true;
    }
    */
    if (whichHydrogen == AllHydrogen)
      obErrorLog.ThrowError(__FUNCTION__,
                            "Ran OpenBabel::AddHydrogens", obAuditMsg);
    else if (whichHydrogen == PolarHydrogen)
      obErrorLog.ThrowError(__FUNCTION__,
                            "Ran OpenBabel::AddHydrogens -- polar only", obAuditMsg);
    else
      obErrorLog.ThrowError(__FUNCTION__,
                            "Ran OpenBabel::AddHydrogens -- nonpolar only", obAuditMsg);

    // Make sure we have conformers (PR#1665519)
    if (!_vconf.empty() && !Empty()) {
      OBAtom *atom;
      vector<OBAtom*>::iterator i;
      for (atom = BeginAtom(i);atom;atom = NextAtom(i))
        {
          atom->SetVector();
        }
    }

    SetHydrogensAdded(); // This must come after EndModify() as EndModify() wipes the flags
    // If chirality was already perceived, remember this (to avoid wiping information
    if (hasChiralityPerceived)
      this->SetChiralityPerceived();

    //count up number of hydrogens to add
    OBAtom *atom,*h;
    int hcount,count=0;
    vector<pair<OBAtom*,int> > vhadd;
    vector<OBAtom*>::iterator i;
    for (atom = BeginAtom(i);atom;atom = NextAtom(i))
      {
        if (whichHydrogen == PolarHydrogen && !AtomIsNSOP(atom))
          continue;
        if (whichHydrogen == NonPolarHydrogen && AtomIsNSOP(atom))
          continue;

        hcount = atom->GetImplicitHCount();
        atom->SetImplicitHCount(0);

        if (hcount)
          {
            vhadd.push_back(pair<OBAtom*,int>(atom,hcount));
            count += hcount;
          }
      }

    if (count == 0) {
      // Make sure to clear SSSR and aromatic flags we may have tripped above
      _flags &= (~(OB_SSSR_MOL|OB_AROMATIC_MOL));
      return(true);
    }
    bool hasCoords = HasNonZeroCoords();

    //realloc memory in coordinate arrays for new hydrogens
    double *tmpf;
    vector<double*>::iterator j;
    for (j = _vconf.begin();j != _vconf.end();++j)
      {
        tmpf = new double [(NumAtoms()+count)*3];
        memset(tmpf,'\0',sizeof(double)*(NumAtoms()+count)*3);
        if (hasCoords)
          memcpy(tmpf,(*j),sizeof(double)*NumAtoms()*3);
        delete []*j;
        *j = tmpf;
      }

    IncrementMod();

    int m,n;
    vector3 v;
    vector<pair<OBAtom*,int> >::iterator k;
    double hbrad = CorrectedBondRad(1, 0);

    for (k = vhadd.begin();k != vhadd.end();++k)
      {
        atom = k->first;
        double bondlen = hbrad + CorrectedBondRad(atom->GetAtomicNum(), atom->GetHyb());
        for (m = 0;m < k->second;++m)
          {
            int badh = 0;
            for (n = 0;n < NumConformers();++n)
              {
                SetConformer(n);
                if (hasCoords)
                  {
                    // Ensure that add hydrogens only returns finite coords
                    //atom->GetNewBondVector(v,bondlen);
                    v = OBBuilder::GetNewBondVector(atom,bondlen);
                    if (isfinite(v.x()) || isfinite(v.y()) || isfinite(v.z())) {
                      _c[(NumAtoms())*3]   = v.x();
                      _c[(NumAtoms())*3+1] = v.y();
                      _c[(NumAtoms())*3+2] = v.z();
                    }
                    else {
                      _c[(NumAtoms())*3]   = 0.0;
                      _c[(NumAtoms())*3+1] = 0.0;
                      _c[(NumAtoms())*3+2] = 0.0;
                      obErrorLog.ThrowError(__FUNCTION__,
                                            "Ran OpenBabel::AddHydrogens -- no reasonable bond geometry for desired hydrogen.",
                                            obAuditMsg);
                      badh++;
                    }
                  }
                else
                  memset((char*)&_c[NumAtoms()*3],'\0',sizeof(double)*3);
              }
            if(badh == 0 || badh < NumConformers()) 
              {
                // Add the new H atom to the appropriate residue list
                //but avoid doing perception by checking for existence of residue
                //just in case perception is trigger, make sure GetResidue is called
                //before adding the hydrogen to the molecule
                OBResidue *res = atom->HasResidue() ? atom->GetResidue() : NULL;
                h = NewAtom();
                h->SetType("H");
                h->SetAtomicNum(1);
                string aname = "H";

                if(res) 
                {
                  res->AddAtom(h);
                  res->SetAtomID(h,aname);
                  
                  //hydrogen should inherit hetatm status of heteroatom (default is false)
                  if(res->IsHetAtom(atom)) 
                  {
                      res->SetHetAtom(h, true);
                  }
                }

                int bondFlags = 0;
                AddBond(atom->GetIdx(),h->GetIdx(),1, bondFlags);
                h->SetCoordPtr(&_c);
                OpenBabel::ImplicitRefToStereo(*this, atom->GetId(), h->GetId());
              }
          }
      }

    DecrementMod();

    //reset atom type and partial charge flags
    _flags &= (~(OB_PCHARGE_MOL|OB_ATOMTYPES_MOL|OB_SSSR_MOL|OB_AROMATIC_MOL|OB_HYBRID_MOL));

    return(true);
  }

  bool OBMol::AddPolarHydrogens()
  {
    return(AddNewHydrogens(PolarHydrogen));
  }

  bool OBMol::AddNonPolarHydrogens()
  {
    return(AddNewHydrogens(NonPolarHydrogen));
  }

  bool OBMol::AddHydrogens(OBAtom *atom)
  {
    int hcount = atom->GetImplicitHCount();
    if (hcount == 0)
      return true;

    atom->SetImplicitHCount(0);

    vector<pair<OBAtom*, int> > vhadd;
    vhadd.push_back(pair<OBAtom*,int>(atom, hcount));

    //realloc memory in coordinate arrays for new hydroges
    double *tmpf;
    vector<double*>::iterator j;
    for (j = _vconf.begin();j != _vconf.end();++j)
      {
        tmpf = new double [(NumAtoms()+hcount)*3+10];
        memcpy(tmpf,(*j),sizeof(double)*NumAtoms()*3);
        delete []*j;
        *j = tmpf;
      }

    IncrementMod();

    int m,n;
    vector3 v;
    vector<pair<OBAtom*,int> >::iterator k;
    double hbrad = CorrectedBondRad(1,0);

    OBAtom *h;
    for (k = vhadd.begin();k != vhadd.end();++k)
      {
        atom = k->first;
        double bondlen = hbrad + CorrectedBondRad(atom->GetAtomicNum(),atom->GetHyb());
        for (m = 0;m < k->second;++m)
          {
            for (n = 0;n < NumConformers();++n)
              {
                SetConformer(n);
                //atom->GetNewBondVector(v,bondlen);
                v = OBBuilder::GetNewBondVector(atom,bondlen);
                _c[(NumAtoms())*3]   = v.x();
                _c[(NumAtoms())*3+1] = v.y();
                _c[(NumAtoms())*3+2] = v.z();
              }
            h = NewAtom();
            h->SetType("H");
            h->SetAtomicNum(1);

            int bondFlags = 0;
            AddBond(atom->GetIdx(),h->GetIdx(),1, bondFlags);
            h->SetCoordPtr(&_c);
            OpenBabel::ImplicitRefToStereo(*this, atom->GetId(), h->GetId());
          }
      }

    DecrementMod();
    SetConformer(0);

    //reset atom type and partial charge flags
    //_flags &= (~(OB_PCHARGE_MOL|OB_ATOMTYPES_MOL));

    return(true);
  }

  bool OBMol::CorrectForPH(double pH)
  {
    if (IsCorrectedForPH())
      return(true);
    phmodel.CorrectForPH(*this, pH);

    obErrorLog.ThrowError(__FUNCTION__,
                          "Ran OpenBabel::CorrectForPH", obAuditMsg);

    return(true);
  }

  //! \brief set spin multiplicity for H-deficient atoms
  /**
     If NoImplicitH is true then the molecule has no implicit hydrogens. Individual atoms
     on which ForceNoH() has been called also have no implicit hydrogens.
     If NoImplicitH is false (the default), then if there are any explicit hydrogens
     on an atom then they constitute all the hydrogen on that atom. However, a hydrogen
     atom with its _isotope!=0 is not considered explicit hydrogen for this purpose.
     In addition, an atom which has had ForceImplH()called for it is never considered
     hydrogen deficient, e.g. unbracketed atoms in SMILES.
     Any discrepancy with the expected atom valency is interpreted as the atom being a
     radical of some sort and iits _spinMultiplicity is set to 2 when it is one hydrogen short
     and 3 when it is two hydrogens short and similarly for greater hydrogen deficiency.

     So SMILES C[CH] is interpreted as methyl carbene, CC[H][H] as ethane, and CC[2H] as CH3CH2D.
  **/



  bool OBMol::AssignSpinMultiplicity(bool NoImplicitH)
  {
    // TODO: The following functions simply returns true, as it has been made
    // redundant by changes to the handling of implicit hydrogens, and spin.
    // This needs to be sorted out properly at some point.
    return true;
  }

  // Used by DeleteAtom below. Code based on StereoRefToImplicit
  static void DeleteStereoOnAtom(OBMol& mol, OBStereo::Ref atomId)
  {
    std::vector<OBGenericData*> vdata = mol.GetAllData(OBGenericDataType::StereoData);
    for (std::vector<OBGenericData*>::iterator data = vdata.begin(); data != vdata.end(); ++data) {
      OBStereo::Type datatype = ((OBStereoBase*)*data)->GetType();

      if (datatype != OBStereo::CisTrans && datatype != OBStereo::Tetrahedral) {
        obErrorLog.ThrowError(__FUNCTION__,
            "This function should be updated to handle additional stereo types.\nSome stereochemistry objects may contain explicit refs to hydrogens which have been removed.", obWarning);
        continue;
      }

      if (datatype == OBStereo::CisTrans) {
        OBCisTransStereo *ct = dynamic_cast<OBCisTransStereo*>(*data);
        OBCisTransStereo::Config ct_cfg = ct->GetConfig();
        if (ct_cfg.begin == atomId || ct_cfg.end == atomId ||
            std::find(ct_cfg.refs.begin(), ct_cfg.refs.end(), atomId) != ct_cfg.refs.end())
          mol.DeleteData(ct);
      }
      else if (datatype == OBStereo::Tetrahedral) {
        OBTetrahedralStereo *ts = dynamic_cast<OBTetrahedralStereo*>(*data);
        OBTetrahedralStereo::Config ts_cfg = ts->GetConfig();
        if (ts_cfg.from == atomId ||
            std::find(ts_cfg.refs.begin(), ts_cfg.refs.end(), atomId) != ts_cfg.refs.end())
          mol.DeleteData(ts);
      }
    }
  }

  bool OBMol::DeleteAtom(OBAtom *atom, bool destroyAtom)
  {
    if (atom->GetAtomicNum() == OBElements::Hydrogen)
      return(DeleteHydrogen(atom));

    BeginModify();
    //don't need to do anything with coordinates b/c
    //BeginModify() blows away coordinates

    //find bonds to delete
    OBAtom *nbr;
    vector<OBBond*> vdb;
    vector<OBBond*>::iterator j;
    for (nbr = atom->BeginNbrAtom(j);nbr;nbr = atom->NextNbrAtom(j))
      vdb.push_back(*j);

    for (j = vdb.begin();j != vdb.end();++j)
      DeleteBond((OBBond *)*j); //delete bonds

    _atomIds[atom->GetId()] = (OBAtom*)NULL;
    _vatom.erase(_vatom.begin()+(atom->GetIdx()-1));
    _natoms--;

    //reset all the indices to the atoms
    int idx;
    vector<OBAtom*>::iterator i;
    OBAtom *atomi;
    for (idx=1,atomi = BeginAtom(i);atomi;atomi = NextAtom(i),++idx)
      atomi->SetIdx(idx);

    EndModify();

    // Delete any stereo objects involving this atom
    OBStereo::Ref id = atom->GetId();
    DeleteStereoOnAtom(*this, id);

    if (destroyAtom)
      DestroyAtom(atom);

    SetSSSRPerceived(false);
    SetLSSRPerceived(false);
    return(true);
  }

  bool OBMol::DeleteResidue(OBResidue *residue, bool destroyResidue)
  {
    unsigned short idx = residue->GetIdx();
    _residue.erase(_residue.begin() + idx);

    for ( unsigned short i = idx ; i < _residue.size() ; i++ )
      _residue[i]->SetIdx(i);

    if (destroyResidue)
      DestroyResidue(residue);

    SetSSSRPerceived(false);
    SetLSSRPerceived(false);
    return(true);
  }

  bool OBMol::DeleteBond(OBBond *bond, bool destroyBond)
  {
    BeginModify();

    (bond->GetBeginAtom())->DeleteBond(bond);
    (bond->GetEndAtom())->DeleteBond(bond);
    _bondIds[bond->GetId()] = (OBBond*)NULL;
    _vbond.erase(_vbond.begin() + bond->GetIdx()); // bond index starts at 0!!!
    _nbonds--;

    vector<OBBond*>::iterator i;
    int j;
    OBBond *bondi;
    for (bondi = BeginBond(i),j=0;bondi;bondi = NextBond(i),++j)
      bondi->SetIdx(j);

    EndModify();

    if (destroyBond)
      DestroyBond(bond);

    SetSSSRPerceived(false);
    SetLSSRPerceived(false);
    return(true);
  }

  bool OBMol::AddBond(int first,int second,int order,int flags,int insertpos)
  {
    // Don't add the bond if it already exists
    if (first == second || GetBond(first, second) != NULL)
      return(false);

    //    BeginModify();

    if ((unsigned)first <= NumAtoms() && (unsigned)second <= NumAtoms())
      //atoms exist and bond doesn't
      {
        OBBond *bond = new OBBond;
        if (!bond)
          {
            //EndModify();
            return(false);
          }

        OBAtom *bgn,*end;
        bgn = GetAtom(first);
        end = GetAtom(second);
        if (!bgn || !end)
          {
            obErrorLog.ThrowError(__FUNCTION__, "Unable to add bond - invalid atom index", obDebug);
            return(false);
          }
        bond->Set(_nbonds,bgn,end,order,flags);
        bond->SetParent(this);

        bond->SetId(_bondIds.size());
        _bondIds.push_back(bond);

#define OBBondIncrement 100
        if (_nbonds+1 >= _vbond.size())
          {
            _vbond.resize(_nbonds+OBBondIncrement);
            vector<OBBond*>::iterator i;
            for (i = _vbond.begin(),i+=(_nbonds+1);i != _vbond.end();++i)
              *i = (OBBond*)NULL;
          }
#undef  OBBondIncrement

        _vbond[_nbonds] = (OBBond*)bond;
        _nbonds++;

        if (insertpos == -1)
          {
            bgn->AddBond(bond);
            end->AddBond(bond);
          }
        else
          {
            if (insertpos >= static_cast<int>(bgn->GetExplicitDegree()))
              bgn->AddBond(bond);
            else //need to insert the bond for the connectivity order to be preserved
              {    //otherwise stereochemistry gets screwed up
                vector<OBBond*>::iterator bi;
                bgn->BeginNbrAtom(bi);
                bi += insertpos;
                bgn->InsertBond(bi,bond);
              }
            end->AddBond(bond);
          }
      }
    else //at least one atom doesn't exist yet - add to bond_q
      SetData(new OBVirtualBond(first,second,order,flags));

    //    EndModify();

    return(true);
  }

  bool OBMol::AddBond(OBBond &bond)
  {
    if(!AddBond(bond.GetBeginAtomIdx(),
                   bond.GetEndAtomIdx(),
                   bond.GetBondOrder(),
                   bond.GetFlags()))
      return false;
    //copy the bond's generic data
    OBDataIterator diter;
    for(diter=bond.BeginData(); diter!=bond.EndData();++diter)
      GetBond(NumBonds()-1)->CloneData(*diter);
    return true;
  }

  void OBMol::Align(OBAtom *a1,OBAtom *a2,vector3 &p1,vector3 &p2)
  {
    vector<int> children;

    obErrorLog.ThrowError(__FUNCTION__,
                          "Ran OpenBabel::Align", obAuditMsg);

    //find which atoms to rotate
    FindChildren(children,a1->GetIdx(),a2->GetIdx());
    children.push_back(a2->GetIdx());

    //find the rotation vector and angle
    vector3 v1,v2,v3;
    v1 = p2 - p1;
    v2 = a2->GetVector() - a1->GetVector();
    v3 = cross(v1,v2);
    double angle = vectorAngle(v1,v2);

    //find the rotation matrix
    matrix3x3 m;
    m.RotAboutAxisByAngle(v3,angle);

    //rotate atoms
    vector3 v;
    OBAtom *atom;
    vector<int>::iterator i;
    for (i = children.begin();i != children.end();++i)
      {
        atom = GetAtom(*i);
        v = atom->GetVector();
        v -= a1->GetVector();
        v *= m;   //rotate the point
        v += p1;  //translate the vector
        atom->SetVector(v);
      }
    //set a1 = p1
    a1->SetVector(p1);
  }

  void OBMol::ToInertialFrame()
  {
    double m[9];
    for (int i = 0;i < NumConformers();++i)
      ToInertialFrame(i,m);
  }

  void OBMol::ToInertialFrame(int conf,double *rmat)
  {
    unsigned int i;
    double x,y,z;
    double mi;
    double mass = 0.0;
    double center[3],m[3][3];

    obErrorLog.ThrowError(__FUNCTION__,
                          "Ran OpenBabel::ToInertialFrame", obAuditMsg);

    for (i = 0;i < 3;++i)
      memset(&m[i],'\0',sizeof(double)*3);
    memset(center,'\0',sizeof(double)*3);

    SetConformer(conf);
    OBAtom *atom;
    vector<OBAtom*>::iterator j;
    //find center of mass
    for (atom = BeginAtom(j);atom;atom = NextAtom(j))
      {
        mi = atom->GetAtomicMass();
        center[0] += mi*atom->x();
        center[1] += mi*atom->y();
        center[2] += mi*atom->z();
        mass += mi;
      }

    center[0] /= mass;
    center[1] /= mass;
    center[2] /= mass;

    //calculate inertial tensor
    for (atom = BeginAtom(j);atom;atom = NextAtom(j))
      {
        x = atom->x()-center[0];
        y = atom->y()-center[1];
        z = atom->z()-center[2];
        mi = atom->GetAtomicMass();

        m[0][0] += mi*(y*y+z*z);
        m[0][1] -= mi*x*y;
        m[0][2] -= mi*x*z;
        //        m[1][0] -= mi*x*y;
        m[1][1] += mi*(x*x+z*z);
        m[1][2] -= mi*y*z;
        //        m[2][0] -= mi*x*z;
        //        m[2][1] -= mi*y*z;
        m[2][2] += mi*(x*x+y*y);
      }
    // Fill in the lower triangle using symmetry across the diagonal
    m[1][0] = m[0][1];
    m[2][0] = m[0][2];
    m[2][1] = m[1][2];

    /* find rotation matrix for moment of inertia */
    ob_make_rmat(m,rmat);

    /* rotate all coordinates */
    double *c = GetConformer(conf);
    for(i=0; i < NumAtoms();++i)
      {
        x = c[i*3]-center[0];
        y = c[i*3+1]-center[1];
        z = c[i*3+2]-center[2];
        c[i*3]   = x*rmat[0] + y*rmat[1] + z*rmat[2];
        c[i*3+1] = x*rmat[3] + y*rmat[4] + z*rmat[5];
        c[i*3+2] = x*rmat[6] + y*rmat[7] + z*rmat[8];
      }
  }

  OBMol::OBMol()
  {
    _natoms = _nbonds = 0;
    _mod = 0;
    _totalCharge = 0;
    _dimension = 3;
    _vatom.clear();
    _atomIds.clear();
    _vbond.clear();
    _bondIds.clear();
    _vdata.clear();
    _title = "";
    _c = (double*)NULL;
    _flags = 0;
    _vconf.clear();
    _autoPartialCharge = true;
    _autoFormalCharge = true;
    _energy = 0.0;
  }

  OBMol::OBMol(const OBMol &mol) : OBBase(mol)
  {
    _natoms = _nbonds = 0;
    _mod = 0;
    _totalCharge = 0;
    _dimension = 3;
    _vatom.clear();
    _atomIds.clear();
    _vbond.clear();
    _bondIds.clear();
    _vdata.clear();
    _title = "";
    _c = (double*)NULL;
    _flags = 0;
    _vconf.clear();
    _autoPartialCharge = true;
    _autoFormalCharge = true;
    //NF  _compressed = false;
    _energy = 0.0;
    *this = mol;
  }

  OBMol::~OBMol()
  {
    OBAtom    *atom;
    OBBond    *bond;
    OBResidue *residue;
    vector<OBAtom*>::iterator i;
    vector<OBBond*>::iterator j;
    vector<OBResidue*>::iterator r;
    for (atom = BeginAtom(i);atom;atom = NextAtom(i))
      DestroyAtom(atom);
    for (bond = BeginBond(j);bond;bond = NextBond(j))
      DestroyBond(bond);
    for (residue = BeginResidue(r);residue;residue = NextResidue(r))
      DestroyResidue(residue);

    //clear out the multiconformer data
    vector<double*>::iterator k;
    for (k = _vconf.begin();k != _vconf.end();++k)
      delete [] *k;
    _vconf.clear();
  }

  bool OBMol::HasNonZeroCoords()
  {
    OBAtom *atom;
    vector<OBAtom*>::iterator i;

    for (atom = BeginAtom(i);atom;atom = NextAtom(i))
      if (atom->GetVector() != VZero)
        return(true);

    return(false);
  }

  bool OBMol::Has2D(bool Not3D)
  {
    bool hasX,hasY;
    OBAtom *atom;
    vector<OBAtom*>::iterator i;

    hasX = hasY = false;
    for (atom = BeginAtom(i);atom;atom = NextAtom(i))
      {
        if (!hasX && !IsNearZero(atom->x()))
          hasX = true;
        if (!hasY && !IsNearZero(atom->y()))
          hasY = true;
        if(Not3D && atom->z())
          return false;
      }
    if (hasX || hasY) //was && but this excluded vertically or horizontally aligned linear mols
      return(true);
    return(false);
  }

  bool OBMol::Has3D()
  {
    bool hasX,hasY,hasZ;
    OBAtom *atom;
    vector<OBAtom*>::iterator i;

    hasX = hasY = hasZ = false;
    //    if (this->_c == NULL) **Test removed** Prevented function use during molecule construction
    //      return(false);
    for (atom = BeginAtom(i);atom;atom = NextAtom(i))
      {
        if (!hasX && !IsNearZero(atom->x()))
          hasX = true;
        if (!hasY && !IsNearZero(atom->y()))
          hasY = true;
        if (!hasZ && !IsNearZero(atom->z()))
          hasZ = true;

        if (hasX && hasY && hasZ)
          return(true);
      }
    return(false);
  }

  void OBMol::SetCoordinates(double *newCoords)
  {
    bool noCptr = (_c == NULL); // did we previously have a coordinate ptr
    if (noCptr) {
      _c = new double [NumAtoms()*3];
    }

    // copy from external to internal
    memcpy((char*)_c, (char*)newCoords, sizeof(double)*3*NumAtoms());

    if (noCptr) {
      OBAtom *atom;
      vector<OBAtom*>::iterator i;
      for (atom = BeginAtom(i);atom;atom = NextAtom(i))
        atom->SetCoordPtr(&_c);
      _vconf.push_back(newCoords);
    }
  }

  //! Renumber the atoms according to the order of indexes in the supplied vector
  //! This with assemble an atom vector and call RenumberAtoms(vector<OBAtom*>)
  //! It will return without action if the supplied vector is empty or does not
  //! have the same number of atoms as the molecule.
  //!
  //! \since version 2.3
  void OBMol::RenumberAtoms(vector<int> v)
  {
    if (Empty() || v.size() != NumAtoms())
      return;

    vector <OBAtom*> va;
    va.reserve(NumAtoms());

    vector<int>::iterator i;
    for (i = v.begin(); i != v.end(); ++i)
      va.push_back( GetAtom(*i) );

    this->RenumberAtoms(va);
  }

  //! Renumber the atoms in this molecule according to the order in the supplied
  //! vector. This will return without action if the supplied vector is empty or
  //! does not have the same number of atoms as the molecule.
  void OBMol::RenumberAtoms(vector<OBAtom*> &v)
  {
    if (Empty())
      return;

    obErrorLog.ThrowError(__FUNCTION__,
                          "Ran OpenBabel::RenumberAtoms", obAuditMsg);

    OBAtom *atom;
    vector<OBAtom*> va;
    vector<OBAtom*>::iterator i;

    va = v;

    //make sure all atoms are represented in the vector
    if (va.empty() || va.size() != NumAtoms())
      return;

    OBBitVec bv;
    for (i = va.begin();i != va.end();++i)
      bv |= (*i)->GetIdx();

    for (atom = BeginAtom(i);atom;atom = NextAtom(i))
      if (!bv[atom->GetIdx()])
        va.push_back(atom);

    int j,k;
    double *c;
    double *ctmp = new double [NumAtoms()*3];

    for (j = 0;j < NumConformers();++j)
      {
        c = GetConformer(j);
        for (k=0,i = va.begin();i != va.end(); ++i,++k)
          memcpy((char*)&ctmp[k*3],(char*)&c[((OBAtom*)*i)->GetCoordinateIdx()],sizeof(double)*3);
        memcpy((char*)c,(char*)ctmp,sizeof(double)*3*NumAtoms());
      }

    for (k=1,i = va.begin();i != va.end(); ++i,++k)
      (*i)->SetIdx(k);

    delete [] ctmp;

    _vatom.clear();
    for (i = va.begin();i != va.end();++i)
      _vatom.push_back(*i);

    DeleteData(OBGenericDataType::RingData);
    DeleteData("OpenBabel Symmetry Classes");
    DeleteData("LSSR");
    DeleteData("SSSR");
    UnsetFlag(OB_LSSR_MOL);
    UnsetFlag(OB_SSSR_MOL);
  }

  bool WriteTitles(ostream &ofs, OBMol &mol)
  {
    ofs << mol.GetTitle() << endl;
    return true;
  }

  //check that unreasonable bonds aren't being added
  static bool validAdditionalBond(OBAtom *a, OBAtom *n) 
  {
    if(a->GetExplicitValence() == 5 && a->GetAtomicNum() == 15) 
    {
      //only allow octhedral bonding for F and Cl
      if(n->GetAtomicNum() == 9 || n->GetAtomicNum() == 17)
        return true;
      else
        return false;
    }
    //other things to check?
    return true;
  }

  /*! This method adds single bonds between all atoms
    closer than their combined atomic covalent radii,
    then "cleans up" making sure bonded atoms are not
    closer than 0.4A and the atom does not exceed its valence.
    It implements blue-obelisk:rebondFrom3DCoordinates.

  */
  void OBMol::ConnectTheDots(void)
  {
    if (Empty())
      return;
    if (_dimension != 3) return; // not useful on non-3D structures

    obErrorLog.ThrowError(__FUNCTION__,
                          "Ran OpenBabel::ConnectTheDots", obAuditMsg);

    int j,k,max;
    double maxrad = 0;
    bool unset = false;
    OBAtom *atom,*nbr;
    vector<OBAtom*>::iterator i;
    vector<pair<OBAtom*,double> > zsortedAtoms;
    vector<double> rad;
    vector<int> zsorted;
    vector<int> bondCount; // existing bonds (e.g., from residues in PDB)

    double *c = new double [NumAtoms()*3];
    rad.resize(_natoms);

    for (j = 0, atom = BeginAtom(i) ; atom ; atom = NextAtom(i), ++j)
      {
        bondCount.push_back(atom->GetExplicitDegree());
        //don't consider atoms with a full valance already
        //this is both for correctness (trust existing bonds) and performance
        if(atom->GetExplicitValence() >= OBElements::GetMaxBonds(atom->GetAtomicNum()))
          continue;        
        if(atom->GetAtomicNum() == 7 && atom->GetFormalCharge() == 0 && atom->GetExplicitValence() >= 3)
          continue; 
        (atom->GetVector()).Get(&c[j*3]);
        pair<OBAtom*,double> entry(atom, atom->GetVector().z());
        zsortedAtoms.push_back(entry);
      }
    sort(zsortedAtoms.begin(), zsortedAtoms.end(), SortAtomZ);

    max = zsortedAtoms.size();

    for ( j = 0 ; j < max ; j++ )
      {
        atom   = zsortedAtoms[j].first;
        rad[j] = OBElements::GetCovalentRad(atom->GetAtomicNum());
        maxrad = std::max(rad[j],maxrad);
        zsorted.push_back(atom->GetIdx()-1);
      }

    int idx1, idx2;
    double d2,cutoff,zd;
    for (j = 0 ; j < max ; ++j)
      {
        double maxcutoff = SQUARE(rad[j]+maxrad+0.45);
        idx1 = zsorted[j];
        for (k = j + 1 ; k < max ; k++ )
          {
            idx2 = zsorted[k];

            // bonded if closer than elemental Rcov + tolerance
            cutoff = SQUARE(rad[j] + rad[k] + 0.45);

            zd  = SQUARE(c[idx1*3+2] - c[idx2*3+2]);
            // bigger than max cutoff, which is determined using largest radius,
            // not the radius of k (which might be small, ie H, and cause an early  termination)
            // since we sort by z, anything beyond k will also fail
            if (zd > maxcutoff )
              break;

            d2  = SQUARE(c[idx1*3]   - c[idx2*3]);
            if (d2 > cutoff)
              continue; // x's bigger than cutoff
            d2 += SQUARE(c[idx1*3+1] - c[idx2*3+1]);
            if (d2 > cutoff)
              continue; // x^2 + y^2 bigger than cutoff
            d2 += zd;

            if (d2 > cutoff)
              continue;
            if (d2 < 0.16) // 0.4 * 0.4 = 0.16
              continue;

            atom = GetAtom(idx1+1);
            nbr  = GetAtom(idx2+1);

            if (atom->IsConnected(nbr))
              continue;

            if (!validAdditionalBond(atom,nbr) || !validAdditionalBond(nbr, atom))
              continue;
              
            AddBond(idx1+1,idx2+1,1);
          }
      }

    // If between BeginModify and EndModify, coord pointers are NULL
    // setup molecule to handle current coordinates

    if (_c == NULL)
      {
        _c = c;
        for (atom = BeginAtom(i);atom;atom = NextAtom(i))
          atom->SetCoordPtr(&_c);
        _vconf.push_back(c);
        unset = true;
      }

    // Cleanup -- delete long bonds that exceed max valence
    OBBond *maxbond, *bond;
    double maxlength;
    vector<OBBond*>::iterator l, m;
    int valCount;
    bool changed;
    BeginModify(); //prevent needless re-perception in DeleteBond
    for (atom = BeginAtom(i);atom;atom = NextAtom(i))
      {
        while (atom->GetExplicitValence() > static_cast<unsigned int>(OBElements::GetMaxBonds(atom->GetAtomicNum()))
               || atom->SmallestBondAngle() < 45.0)
          {
            bond = atom->BeginBond(l);
            maxbond = bond;
            // Fix from Liu Zhiguo 2008-01-26
            // loop past any bonds
            // which existed before ConnectTheDots was called
            // (e.g., from PDB resdata.txt)
            valCount = 0;
            while (valCount < bondCount[atom->GetIdx() - 1]) {
              bond = atom->NextBond(l);
              // timvdm: 2008-03-05
              // NextBond only returns NULL if the iterator l == _bonds.end().
              // This was casuing problems as follows:
              // NextBond = 0x????????
              // NextBond = 0x????????
              // NextBond = 0x????????
              // NextBond = 0x????????
              // NextBond = NULL  <-- this NULL was not detected
              // NextBond = 0x????????
              if (!bond) // so we add an additional check
                break;
              maxbond = bond;
              valCount++;
            }
            if (!bond) // no new bonds added for this atom, just skip it
              break;

            // delete bonds between hydrogens when over max valence
            if (atom->GetAtomicNum() == OBElements::Hydrogen)
              {
                m = l;
                changed = false;
                for (;bond;bond = atom->NextBond(m))
                  {
                    if (bond->GetNbrAtom(atom)->GetAtomicNum() == OBElements::Hydrogen)
                      {
                        DeleteBond(bond);
                        changed = true;
                        break;
                      }
                  }
                if (changed)
                  {
                    // bond deleted, reevaluate BOSum
                    continue;
                  }
                else
                  {
                    // reset to first new bond
                    bond = maxbond;
                  }
              }

            maxlength = maxbond->GetLength();
            for (bond = atom->NextBond(l);bond;bond = atom->NextBond(l))
              {
                if (bond->GetLength() > maxlength)
                  {
                    maxbond = bond;
                    maxlength = bond->GetLength();
                  }
              }
            DeleteBond(maxbond); // delete the new bond with the longest length
          }
      }
    EndModify();
    if (unset)
      {
        _c = NULL;
        for (atom = BeginAtom(i);atom;atom = NextAtom(i))
          atom->ClearCoordPtr();
        _vconf.resize(_vconf.size()-1);
      }

    if (_c != NULL)
      delete [] c;
  }

  /*! This method uses bond angles and geometries from current
    connectivity to guess atom types and then filling empty valences
    with multiple bonds. It currently has a pass to detect some
    frequent functional groups. It still needs a pass to detect aromatic
    rings to "clean up."
    AssignSpinMultiplicity(true) is called at the end of the function. The true
    states that there are no implict hydrogens in the molecule.
  */
  void OBMol::PerceiveBondOrders()
  {
    if (Empty())
      return;
    if (_dimension != 3) return; // not useful on non-3D structures

    obErrorLog.ThrowError(__FUNCTION__,
                          "Ran OpenBabel::PerceiveBondOrders", obAuditMsg);

    OBAtom *atom, *b, *c;
    vector3 v1, v2;
    double angle;//, dist1, dist2;
    vector<OBAtom*>::iterator i;
    vector<OBBond*>::iterator j;//,k;

    //  BeginModify();

    // Pass 1: Assign estimated hybridization based on avg. angles
    for (atom = BeginAtom(i);atom;atom = NextAtom(i))
      {
        angle = atom->AverageBondAngle();

        //        cout << atom->GetAtomicNum() << " " << angle << endl;

        if (angle > 155.0)
          atom->SetHyb(1);
        else if (angle <= 155.0 && angle > 115.0)
          atom->SetHyb(2);

        // special case for imines
        if (atom->GetAtomicNum() == OBElements::Nitrogen
            && atom->ExplicitHydrogenCount() == 1
            && atom->GetExplicitDegree() == 2
            && angle > 109.5)
          atom->SetHyb(2);

      } // pass 1

    // Make sure upcoming calls to GetHyb() don't kill these temporary values
    SetHybridizationPerceived();

    // Pass 2: look for 5-member rings with torsions <= 7.5 degrees
    //         and 6-member rings with torsions <= 12 degrees
    //         (set all atoms with at least two bonds to sp2)

    vector<OBRing*> rlist;
    vector<OBRing*>::iterator ringit;
    vector<int> path;
    double torsions = 0.0;

    if (!HasSSSRPerceived())
      FindSSSR();
    rlist = GetSSSR();
    for (ringit = rlist.begin(); ringit != rlist.end(); ++ringit)
      {
        if ((*ringit)->PathSize() == 5)
          {
            path = (*ringit)->_path;
            torsions =
              ( fabs(GetTorsion(path[0], path[1], path[2], path[3])) +
                fabs(GetTorsion(path[1], path[2], path[3], path[4])) +
                fabs(GetTorsion(path[2], path[3], path[4], path[0])) +
                fabs(GetTorsion(path[3], path[4], path[0], path[1])) +
                fabs(GetTorsion(path[4], path[0], path[1], path[2])) ) / 5.0;
            if (torsions <= 7.5)
              {
                for (unsigned int ringAtom = 0; ringAtom != path.size(); ++ringAtom)
                  {
                    b = GetAtom(path[ringAtom]);
                    // if an aromatic ring atom has valence 3, it is already set
                    // to sp2 because the average angles should be 120 anyway
                    // so only look for valence 2
                    if (b->GetExplicitDegree() == 2)
                      b->SetHyb(2);
                  }
              }
          }
        else if ((*ringit)->PathSize() == 6)
          {
            path = (*ringit)->_path;
            torsions =
              ( fabs(GetTorsion(path[0], path[1], path[2], path[3])) +
                fabs(GetTorsion(path[1], path[2], path[3], path[4])) +
                fabs(GetTorsion(path[2], path[3], path[4], path[5])) +
                fabs(GetTorsion(path[3], path[4], path[5], path[0])) +
                fabs(GetTorsion(path[4], path[5], path[0], path[1])) +
                fabs(GetTorsion(path[5], path[0], path[1], path[2])) ) / 6.0;
            if (torsions <= 12.0)
              {
                for (unsigned int ringAtom = 0; ringAtom != path.size(); ++ringAtom)
                  {
                    b = GetAtom(path[ringAtom]);
                    if (b->GetExplicitDegree() == 2 || b->GetExplicitDegree() == 3)
                      b->SetHyb(2);
                  }
              }
          }
      }

    // Pass 3: "Antialiasing" If an atom marked as sp hybrid isn't
    //          bonded to another or an sp2 hybrid isn't bonded
    //          to another (or terminal atoms in both cases)
    //          mark them to a lower hybridization for now
    bool openNbr;
    for (atom = BeginAtom(i);atom;atom = NextAtom(i))
      {
        if (atom->GetHyb() == 2 || atom->GetHyb() == 1)
          {
            openNbr = false;
            for (b = atom->BeginNbrAtom(j); b; b = atom->NextNbrAtom(j))
              {
                if (b->GetHyb() < 3 || b->GetExplicitDegree() == 1)
                  {
                    openNbr = true;
                    break;
                  }
              }
            if (!openNbr && atom->GetHyb() == 2)
              atom->SetHyb(3);
            else if (!openNbr && atom->GetHyb() == 1)
              atom->SetHyb(2);
          }
      } // pass 3

    // Pass 4: Check for known functional group patterns and assign bonds
    //         to the canonical form
    //      Currently we have explicit code to do this, but a "bond typer"
    //      is in progress to make it simpler to test and debug.
    bondtyper.AssignFunctionalGroupBonds(*this);

    // Pass 5: Check for aromatic rings and assign bonds as appropriate
    // This is just a quick and dirty approximation that marks everything
    //  as potentially aromatic

    // This doesn't work perfectly, but it's pretty decent.
    //  Need to have a list of SMARTS patterns for common rings
    //  which would "break ties" on complicated multi-ring systems
    // (Most of the current problems lie in the interface with the
    //   Kekulize code anyway, not in marking everything as potentially aromatic)

    bool needs_kekulization = false; // are there any aromatic bonds?
    bool typed; // has this ring been typed?
    unsigned int loop, loopSize;
    for (ringit = rlist.begin(); ringit != rlist.end(); ++ringit)
      {
        typed = false;
        loopSize = (*ringit)->PathSize();
        if (loopSize == 5 || loopSize == 6 || loopSize == 7)
          {
            path = (*ringit)->_path;
            for(loop = 0; loop < loopSize; ++loop)
              {
                atom = GetAtom(path[loop]);
                if(atom->HasBondOfOrder(2) || atom->HasBondOfOrder(3)
                   || atom->GetHyb() != 2)
                  {
                    typed = true;
                    break;
                  }
              }

            if (!typed)
              for(loop = 0; loop < loopSize; ++loop)
                {
                  //    cout << " set aromatic " << path[loop] << endl;
                  (GetBond(path[loop], path[(loop+1) % loopSize]))->SetAromatic();
                  needs_kekulization = true;
                }
          }
      }

    // Kekulization is neccessary if an aromatic bond is present
    if (needs_kekulization) {
      this->SetAromaticPerceived();
      // First of all, set the atoms at the ends of the aromatic bonds to also
      // be aromatic. This information is required for OBKekulize.
      FOR_BONDS_OF_MOL(bond, this) {
        if (bond->IsAromatic()) {
          bond->GetBeginAtom()->SetAromatic();
          bond->GetEndAtom()->SetAromatic();
        }
      }
      bool ok = OBKekulize(this);
      if (!ok) {
        stringstream errorMsg;
        errorMsg << "Failed to kekulize aromatic bonds in OBMol::PerceiveBondOrders";
        std::string title = this->GetTitle();
        if (!title.empty())
          errorMsg << " (title is " << title << ")";
        errorMsg << endl;
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
        // return false; Should we return false for a kekulization failure?
      }
      this->SetAromaticPerceived(false);
    }

    // Quick pass.. eliminate inter-ring sulfur atom multiple bonds
    // dkoes - I have removed this code - if double bonds are set,
    // we should trust them.  See pdb_ligands_sdf/4iph_1fj.sdf for
    // a case where the charge isn't set, but we break the molecule
    // if we remove the double bond.  Also, the previous code was
    // fragile - relying on the total mol charge being set.  If we
    // are going to do anything, we should "perceive" a formal charge
    // in the case of a ring sulfur with a double bond (thiopyrylium)

    // Pass 6: Assign remaining bond types, ordered by atom electronegativity
    vector<pair<OBAtom*,double> > sortedAtoms;
    vector<double> rad;
    vector<int> sorted;
    int iter, max;
    double maxElNeg, shortestBond, currentElNeg;
    double bondLength, testLength;

    for (atom = BeginAtom(i) ; atom ; atom = NextAtom(i))
      {
        // if atoms have the same electronegativity, make sure those with shorter bonds
        // are handled first (helps with assignment of conjugated single/double bonds)
        shortestBond = 1.0e5;
        for (b = atom->BeginNbrAtom(j); b; b = atom->NextNbrAtom(j))
          {
            if (b->GetAtomicNum()!=1) shortestBond =
                                        std::min(shortestBond,(atom->GetBond(b))->GetLength());
          }
        pair<OBAtom*,double> entry(atom,
                                   OBElements::GetElectroNeg(atom->GetAtomicNum())*1e6+shortestBond);

        sortedAtoms.push_back(entry);
      }
    sort(sortedAtoms.begin(), sortedAtoms.end(), SortAtomZ);

    max = sortedAtoms.size();
    for (iter = 0 ; iter < max ; iter++ )
      {
        atom = sortedAtoms[iter].first;
        // Debugging statement
        //        cout << " atom->Hyb " << atom->GetAtomicNum() << " " << atom->GetIdx() << " " << atom->GetHyb()
        //             << " BO: " << atom->GetExplicitValence() << endl;

        // Possible sp-hybrids
        if ( (atom->GetHyb() == 1 || atom->GetExplicitDegree() == 1)
             && atom->GetExplicitValence() + 2  <= static_cast<unsigned int>(OBElements::GetMaxBonds(atom->GetAtomicNum()))
             )
          {

            // loop through the neighbors looking for a hybrid or terminal atom
            // (and pick the one with highest electronegativity first)
            // *or* pick a neighbor that's a terminal atom
            if (atom->HasNonSingleBond() ||
                (atom->GetAtomicNum() == 7 && atom->GetExplicitValence() + 2 > 3))
              continue;

            maxElNeg = 0.0;
            shortestBond = 5000.0;
            c = NULL;
            for (b = atom->BeginNbrAtom(j); b; b = atom->NextNbrAtom(j))
              {
                currentElNeg = OBElements::GetElectroNeg(b->GetAtomicNum());
                if ( (b->GetHyb() == 1 || b->GetExplicitDegree() == 1)
                     && b->GetExplicitValence() + 2 <= static_cast<unsigned int>(OBElements::GetMaxBonds(b->GetAtomicNum()))
                     && (currentElNeg > maxElNeg ||
                         (IsApprox(currentElNeg,maxElNeg, 1.0e-6)
                          && (atom->GetBond(b))->GetLength() < shortestBond)) )
                  {
                    if (b->HasNonSingleBond() ||
                        (b->GetAtomicNum() == 7 && b->GetExplicitValence() + 2 > 3))
                      continue;

                    // Test terminal bonds against expected triple bond lengths
                    bondLength = (atom->GetBond(b))->GetLength();
                    if (atom->GetExplicitDegree() == 1 || b->GetExplicitDegree() == 1) {
                      testLength = CorrectedBondRad(atom->GetAtomicNum(), atom->GetHyb())
                        + CorrectedBondRad(b->GetAtomicNum(), b->GetHyb());
                      if (bondLength > 0.9 * testLength)
                        continue; // too long, ignore it
                    }

                    shortestBond = bondLength;
                    maxElNeg = OBElements::GetElectroNeg(b->GetAtomicNum());
                    c = b; // save this atom for later use
                  }
              }
            if (c)
              (atom->GetBond(c))->SetBondOrder(3);
          }
        // Possible sp2-hybrid atoms
        else if ( (atom->GetHyb() == 2 || atom->GetExplicitDegree() == 1)
                  && atom->GetExplicitValence() + 1 <= static_cast<unsigned int>(OBElements::GetMaxBonds(atom->GetAtomicNum())) )
          {
            // as above
            if (atom->HasNonSingleBond() ||
                (atom->GetAtomicNum() == 7 && atom->GetExplicitValence() + 1 > 3))
              continue;

            // Don't build multiple bonds to ring sulfurs
            //  except thiopyrylium
            if (atom->IsInRing() && atom->GetAtomicNum() == 16) {
              if (_totalCharge > 1 && atom->GetFormalCharge() == 0)
                atom->SetFormalCharge(+1);
              else
                continue;
            }

            maxElNeg = 0.0;
            shortestBond = 5000.0;
            c = NULL;
            for (b = atom->BeginNbrAtom(j); b; b = atom->NextNbrAtom(j))
              {
                currentElNeg = OBElements::GetElectroNeg(b->GetAtomicNum());
                if ( (b->GetHyb() == 2 || b->GetExplicitDegree() == 1)
                     && b->GetExplicitValence() + 1 <= static_cast<unsigned int>(OBElements::GetMaxBonds(b->GetAtomicNum()))
                     && (GetBond(atom, b))->IsDoubleBondGeometry()
                     && (currentElNeg > maxElNeg || (IsApprox(currentElNeg,maxElNeg, 1.0e-6)) ) )
                  {
                    if (b->HasNonSingleBond() ||
                        (b->GetAtomicNum() == 7 && b->GetExplicitValence() + 1 > 3))
                      continue;

                    if (b->IsInRing() && b->GetAtomicNum() == 16) {
                      if (_totalCharge > 1 && b->GetFormalCharge() == 0)
                        b->SetFormalCharge(+1);
                      else
                        continue;
                    }

                    // Test terminal bonds against expected double bond lengths
                    bondLength = (atom->GetBond(b))->GetLength();
                    if (atom->GetExplicitDegree() == 1 || b->GetExplicitDegree() == 1) {
                      testLength = CorrectedBondRad(atom->GetAtomicNum(), atom->GetHyb())
                        + CorrectedBondRad(b->GetAtomicNum(), b->GetHyb());
                      if (bondLength > 0.93 * testLength)
                        continue; // too long, ignore it
                    }

                    // OK, see if this is better than the previous choice
                    // If it's much shorter, pick it (e.g., fulvene)
                    // If they're close (0.1A) then prefer the bond in the ring
                    double difference = shortestBond - (atom->GetBond(b))->GetLength();
                    if ( (difference > 0.1)
                         || ( (difference > -0.01) &&
                              ( (!atom->IsInRing() || !c || !c->IsInRing() || b->IsInRing())
                                || (atom->IsInRing() && c && !c->IsInRing() && b->IsInRing()) ) ) ) {
                      shortestBond = (atom->GetBond(b))->GetLength();
                      maxElNeg = OBElements::GetElectroNeg(b->GetAtomicNum());
                      c = b; // save this atom for later use
                    } // is this bond better than previous choices
                  }
              } // loop through neighbors
            if (c)
              (atom->GetBond(c))->SetBondOrder(2);
          }
      } // pass 6

    // Now let the atom typer go to work again
    _flags &= (~(OB_HYBRID_MOL));
    _flags &= (~(OB_AROMATIC_MOL));
    _flags &= (~(OB_ATOMTYPES_MOL));
    //  EndModify(true); // "nuke" perceived data

    //Set _spinMultiplicity other than zero for atoms which are hydrogen
    //deficient and which have implicit valency definitions (essentially the
    //organic subset in SMILES). There are assumed to no implicit hydrogens.
    //AssignSpinMultiplicity(true); // TODO: sort out radicals
    }

  void OBMol::Center()
  {
    for (int i = 0;i < NumConformers();++i)
      Center(i);
  }

  vector3 OBMol::Center(int nconf)
  {
    obErrorLog.ThrowError(__FUNCTION__,
                          "Ran OpenBabel::Center", obAuditMsg);

    SetConformer(nconf);

    OBAtom *atom;
    vector<OBAtom*>::iterator i;

    double x=0.0,y=0.0,z=0.0;
    for (atom = BeginAtom(i);atom;atom = NextAtom(i))
      {
        x += atom->x();
        y += atom->y();
        z += atom->z();
      }

    x /= (double)NumAtoms();
    y /= (double)NumAtoms();
    z /= (double)NumAtoms();

    vector3 vtmp;
    vector3 v(x,y,z);

    for (atom = BeginAtom(i);atom;atom = NextAtom(i))
      {
        vtmp = atom->GetVector() - v;
        atom->SetVector(vtmp);
      }

    return(v);
  }


  /*! this method adds the vector v to all atom positions in all conformers */
  void OBMol::Translate(const vector3 &v)
  {
    for (int i = 0;i < NumConformers();++i)
      Translate(v,i);
  }

  /*! this method adds the vector v to all atom positions in the
    conformer nconf. If nconf == OB_CURRENT_CONFORMER, then the atom
    positions in the current conformer are translated. */
  void OBMol::Translate(const vector3 &v, int nconf)
  {
    obErrorLog.ThrowError(__FUNCTION__,
                          "Ran OpenBabel::Translate", obAuditMsg);

    int i,size;
    double x,y,z;
    double *c = (nconf == OB_CURRENT_CONFORMER)? _c : GetConformer(nconf);

    x = v.x();
    y = v.y();
    z = v.z();
    size = NumAtoms();
    for (i = 0;i < size;++i)
      {
        c[i*3  ] += x;
        c[i*3+1] += y;
        c[i*3+2] += z;
      }
  }

  void OBMol::Rotate(const double u[3][3])
  {
    int i,j,k;
    double m[9];
    for (k=0,i = 0;i < 3;++i)
      for (j = 0;j < 3;++j)
        m[k++] = u[i][j];

    for (i = 0;i < NumConformers();++i)
      Rotate(m,i);
  }

  void OBMol::Rotate(const double m[9])
  {
    for (int i = 0;i < NumConformers();++i)
      Rotate(m,i);
  }

  void OBMol::Rotate(const double m[9],int nconf)
  {
    int i,size;
    double x,y,z;
    double *c = (nconf == OB_CURRENT_CONFORMER)? _c : GetConformer(nconf);

    obErrorLog.ThrowError(__FUNCTION__,
                          "Ran OpenBabel::Rotate", obAuditMsg);

    size = NumAtoms();
    for (i = 0;i < size;++i)
      {
        x = c[i*3  ];
        y = c[i*3+1];
        z = c[i*3+2];
        c[i*3  ] = m[0]*x + m[1]*y + m[2]*z;
        c[i*3+1] = m[3]*x + m[4]*y + m[5]*z;
        c[i*3+2] = m[6]*x + m[7]*y + m[8]*z;
      }
  }

  void OBMol::SetEnergies(std::vector<double> &energies)
  {
    if (!HasData(OBGenericDataType::ConformerData))
      SetData(new OBConformerData);
    OBConformerData *cd = (OBConformerData*) GetData(OBGenericDataType::ConformerData);
    cd->SetEnergies(energies);
  }

  vector<double> OBMol::GetEnergies()
  {
    if (!HasData(OBGenericDataType::ConformerData))
      SetData(new OBConformerData);
    OBConformerData *cd = (OBConformerData*) GetData(OBGenericDataType::ConformerData);
    vector<double> energies = cd->GetEnergies();

    return energies;
  }

  double OBMol::GetEnergy(int ci)
  {
    if (!HasData(OBGenericDataType::ConformerData))
      SetData(new OBConformerData);
    OBConformerData *cd = (OBConformerData*) GetData(OBGenericDataType::ConformerData);
    vector<double> energies = cd->GetEnergies();

    if (((unsigned int)ci >= energies.size()) || (ci < 0))
      return 0.0;

    return energies[ci];
  }

  void OBMol::SetConformers(vector<double*> &v)
  {
    vector<double*>::iterator i;
    for (i = _vconf.begin();i != _vconf.end();++i)
      delete [] *i;

    _vconf = v;
    _c = (_vconf.empty()) ? NULL : _vconf[0];

  }

  void OBMol::SetConformer(unsigned int i)
  {
    if (i < _vconf.size())
      _c = _vconf[i];
  }

  void OBMol::CopyConformer(double *c,int idx)
  {
    //    obAssert(!_vconf.empty() && (unsigned)idx < _vconf.size());
    memcpy((char*)c, (char*)_vconf[idx], sizeof(double)*3*NumAtoms());
  }

  // void OBMol::CopyConformer(double *c,int idx)
  // {
  //   obAssert(!_vconf.empty() && (unsigned)idx < _vconf.size());

  //   unsigned int i;
  //   for (i = 0;i < NumAtoms();++i)
  //     {
  //       _vconf[idx][i*3  ] = (double)c[i*3  ];
  //       _vconf[idx][i*3+1] = (double)c[i*3+1];
  //       _vconf[idx][i*3+2] = (double)c[i*3+2];
  //     }
  // }

  void OBMol::DeleteConformer(int idx)
  {
    if (idx < 0 || idx >= (signed)_vconf.size())
      return;

    delete [] _vconf[idx];
    _vconf.erase((_vconf.begin()+idx));
  }

  ///Converts for instance [N+]([O-])=O to N(=O)=O
  bool OBMol::ConvertDativeBonds()
  {
    obErrorLog.ThrowError(__FUNCTION__,
                          "Ran OpenBabel::ConvertDativeBonds", obAuditMsg);

    //Look for + and - charges on adjacent atoms
    OBAtom* patom;
    vector<OBAtom*>::iterator i;
    bool converted = false;
    for (patom = BeginAtom(i);patom;patom = NextAtom(i))
      {
        vector<OBBond*>::iterator itr;
        OBBond *pbond;
        for (pbond = patom->BeginBond(itr);patom->GetFormalCharge() && pbond;pbond = patom->NextBond(itr))
          {
            OBAtom* pNbratom = pbond->GetNbrAtom(patom);
            int chg1 = patom->GetFormalCharge();
            int chg2 = pNbratom->GetFormalCharge();
            if((chg1>0 && chg2<0)|| (chg1<0 && chg2>0))
              {
                //dative bond. Reduce charges and increase bond order
                converted =true;
                if(chg1>0)
                  --chg1;
                else
                  ++chg1;
                patom->SetFormalCharge(chg1);
                if(chg2>0)
                  --chg2;
                else
                  ++chg2;
                pNbratom->SetFormalCharge(chg2);
                pbond->SetBondOrder(pbond->GetBondOrder()+1);
              }
          }
      }
    return converted; //false if no changes made
  }

  static bool IsNotCorH(OBAtom* atom)
  {
    switch (atom->GetAtomicNum())
    {
    case OBElements::Hydrogen:
    case OBElements::Carbon:
      return false;
    }
    return true;
  }

  //This maybe would be better using smirks from a datafile
  bool OBMol::MakeDativeBonds()
  {
    //! Converts 5-valent N to charged form of dative bonds,
    //! e.g. -N(=O)=O converted to -[N+]([O-])=O. Returns true if conversion occurs.
    BeginModify();
    //AddHydrogens();
    bool converted = false;
    OBAtom* patom;
    vector<OBAtom*>::iterator ai;
    for (patom = BeginAtom(ai);patom;patom = NextAtom(ai)) //all atoms
    {
      if(patom->GetAtomicNum() == OBElements::Nitrogen // || patom->GetAtomicNum() == OBElements::Phosphorus) not phosphorus!
        && (patom->GetExplicitValence()==5 || (patom->GetExplicitValence()==4 && patom->GetFormalCharge()==0)))
      {
        // Find the bond to be modified. Prefer a bond to a hetero-atom,
        // and the highest order bond if there is a choice.
        OBBond *bond, *bestbond;
        OBBondIterator bi;
        for (bestbond = bond = patom->BeginBond(bi); bond; bond = patom->NextBond(bi))
        {
          unsigned int bo = bond->GetBondOrder();
          if(bo>=2 && bo<=4)
          {
            bool het = IsNotCorH(bond->GetNbrAtom(patom));
            bool oldhet = IsNotCorH(bestbond->GetNbrAtom(patom));
            bool higherorder = bo > bestbond->GetBondOrder();
            if((het && !oldhet) || (((het && oldhet) || (!het && !oldhet)) && higherorder))
              bestbond = bond;
          }
        }
        //Make the charged form
        bestbond->SetBondOrder(bestbond->GetBondOrder()-1);
        patom->SetFormalCharge(+1);
        OBAtom* at = bestbond->GetNbrAtom(patom);
        at->SetFormalCharge(-1);
        converted=true;
      }
    }
    EndModify();
    return converted;
  }

  /**
   *  This function is useful when writing to legacy formats (such as MDL MOL) that do
   *  not support zero-order bonds. It is worth noting that some compounds cannot be
   *  well represented using just single, double and triple bonds, even with adjustments
   *  to adjacent charges. In these cases, simply converting zero-order bonds to single
   *  bonds is all that can be done.
   *
   @verbatim
   Algorithm from:
   Clark, A. M. Accurate Specification of Molecular Structures: The Case for
   Zero-Order Bonds and Explicit Hydrogen Counting. Journal of Chemical Information
   and Modeling, 51, 3149-3157 (2011). http://pubs.acs.org/doi/abs/10.1021/ci200488k
   @endverbatim
  */
  bool OBMol::ConvertZeroBonds()
  {
    // TODO: Option to just remove zero-order bonds entirely

    // TODO: Is it OK to not wrap this in BeginModify() and EndModify()?
    // If we must, I think we need to manually remember HasImplicitValencePerceived and
    // re-set it after EndModify()

    // Periodic table block for element (1=s, 2=p, 3=d, 4=f)
    const int BLOCKS[113] = {0,1,2,1,1,2,2,2,2,2,2,1,1,2,2,2,2,2,2,1,1,3,3,3,3,3,3,3,3,3,
                             3,2,2,2,2,2,2,1,1,3,3,3,3,3,3,3,3,3,3,2,2,2,2,2,2,1,1,4,4,4,
                             4,4,4,4,4,4,4,4,4,4,4,3,3,3,3,3,3,3,3,3,3,2,2,2,2,2,2,1,1,4,
                             4,4,4,4,4,4,4,4,4,4,4,4,4,3,3,3,3,3,3,3,3,3,3};

    bool converted = false;
    // Get contiguous fragments of molecule
    vector<vector<int> > cfl;
    ContigFragList(cfl);
    // Iterate over contiguous fragments
    for (vector< vector<int> >::iterator i = cfl.begin(); i != cfl.end(); ++i) {
      // Get all zero-order bonds in contiguous fragment
      vector<OBBond*> bonds;
      for(vector<int>::const_iterator j = i->begin(); j != i->end(); ++j) {
        FOR_BONDS_OF_ATOM(b, GetAtom(*j)) {
          if (b->GetBondOrder() == 0 && !(find(bonds.begin(), bonds.end(), &*b) != bonds.end())) {
            bonds.push_back(&*b);
          }
        }
      }
      // Convert zero-order bonds
      while (bonds.size() > 0) {
        // Pick a bond using scoring system
        int bi = 0;
        if (bonds.size() > 1) {
          vector<int> scores(bonds.size());
          for (unsigned int n = 0; n < bonds.size(); n++) {
            OBAtom *bgn = bonds[n]->GetBeginAtom();
            OBAtom *end = bonds[n]->GetEndAtom();
            int score = 0;
            score += bgn->GetAtomicNum() + end->GetAtomicNum();
            score += abs(bgn->GetFormalCharge()) + abs(end->GetFormalCharge());
            pair<int, int> lb = bgn->LewisAcidBaseCounts();
            pair<int, int> le = end->LewisAcidBaseCounts();
            if (lb.first > 0 && lb.second > 0 && le.first > 0 && le.second > 0) {
              score += 100;   // Both atoms are Lewis acids *and* Lewis bases
            } else if ((lb.first > 0 && le.second > 0) && (lb.second > 0 && le.first > 0)) {
              score -= 1000;  // Lewis acid/base direction is mono-directional
            }
            int bcount = bgn->GetImplicitHCount();
            FOR_BONDS_OF_ATOM(b, bgn) { bcount += 1; }
            int ecount = end->GetImplicitHCount();
            FOR_BONDS_OF_ATOM(b, end) { ecount += 1; }
            if (bcount == 1 || ecount == 1) {
              score -= 10; // If the start or end atoms have only 1 neighbour
            }
            scores[n] = score;
          }
          for (unsigned int n = 1; n < scores.size(); n++) {
            if (scores[n] < scores[bi]) {
              bi = n;
            }
          }
        }
        OBBond *bond = bonds[bi];
        bonds.erase(bonds.begin() + bi);
        OBAtom *bgn = bond->GetBeginAtom();
        OBAtom *end = bond->GetEndAtom();
        int blockb = BLOCKS[bgn->GetAtomicNum()];
        int blocke = BLOCKS[end->GetAtomicNum()];;
        pair<int, int> lb = bgn->LewisAcidBaseCounts();
        pair<int, int> le = end->LewisAcidBaseCounts();
        int chg = 0; // Amount to adjust atom charges
        int ord = 1; // New bond order
        if (lb.first > 0 && lb.second > 0 && le.first > 0 && le.second > 0) {
          ord = 2;  // both atoms are amphoteric, so turn it into a double bond
        } else if (lb.first > 0 && blockb == 2 && blocke >= 3) {
          ord = 2;  // p-block lewis acid with d/f-block element: make into double bond
        } else if (le.first > 0 && blocke == 2 && blockb >= 3) {
          ord = 2;  // p-block lewis acid with d/f-block element: make into double bond
        } else if (lb.first > 0 && le.second > 0) {
          chg = -1;  // lewis acid/base goes one way only; charge separate it
        } else if (lb.second > 0 && le.first > 0) {
          chg = 1;  //  no matching capacity; do not charge separate
        }
        // adjust bond order and atom charges accordingly
        bgn->SetFormalCharge(bgn->GetFormalCharge()+chg);
        end->SetFormalCharge(end->GetFormalCharge()-chg);
        bond->SetBondOrder(ord);
        converted = true;
      }
    }
    return converted;
  }

  OBAtom *OBMol::BeginAtom(OBAtomIterator &i)
  {
    i = _vatom.begin();
    return((i == _vatom.end()) ? (OBAtom*)NULL : (OBAtom*)*i);
  }

  OBAtom *OBMol::NextAtom(OBAtomIterator &i)
  {
    ++i;
    return((i == _vatom.end()) ? (OBAtom*)NULL : (OBAtom*)*i);
  }

  OBBond *OBMol::BeginBond(OBBondIterator &i)
  {
    i = _vbond.begin();
    return((i == _vbond.end()) ? (OBBond*)NULL : (OBBond*)*i);
  }

  OBBond *OBMol::NextBond(OBBondIterator &i)
  {
    ++i;
    return((i == _vbond.end()) ? (OBBond*)NULL : (OBBond*)*i);
  }

  //! \since version 2.4
  int OBMol::AreInSameRing(OBAtom *a, OBAtom *b)
  {
    bool a_in, b_in;
    vector<OBRing*> vr;
    vr = GetLSSR();

    vector<OBRing*>::iterator i;
    vector<int>::iterator j;

    for (i = vr.begin();i != vr.end();++i) {
      a_in = false;
      b_in = false;
      // Go through the path of the ring and see if a and/or b match
      // each node in the path
      for(j = (*i)->_path.begin();j != (*i)->_path.end();++j) {
        if ((unsigned)(*j) == a->GetIdx())
          a_in = true;
        if ((unsigned)(*j) == b->GetIdx())
          b_in = true;
      }

      if (a_in && b_in)
        return (*i)->Size();
    }

    return 0;
  }

  vector<OBMol> OBMol::Separate(int StartIndex)
  {
    vector<OBMol> result;
    if( NumAtoms() == 0 )
      return result; // nothing to do, but let's prevent a crash

    OBMolAtomDFSIter iter( this, StartIndex );
    OBMol newMol;
    while( GetNextFragment( iter, newMol ) ) {
      result.push_back( newMol );
      newMol.Clear();
    }

    return result;
  }

  //! \brief Copy part of a molecule to another molecule
  /**
  This function copies a substructure of a molecule to another molecule. The key
  information needed is an OBBitVec indicating which atoms to include and (optionally)
  an OBBitVec indicating which bonds to exclude. By default, only bonds joining
  included atoms are copied.

  When an atom is copied, but not all of its bonds are, by default hydrogen counts are
  adjusted to account for the missing bonds. That is, given the SMILES "CF", if we
  copy the two atoms but exclude the bond, we will end up with "C.F". This behavior
  can be changed by specifiying a value other than 1 for the \p correctvalence parameter.
  A value of 0 will yield "[C].[F]" while 2 will yield "C*.F*" (see \p correctvalence below
  for more information).

  Aromaticity is preserved as present in the original OBMol. If this is not desired,
  the user should call OBMol::SetAromaticPerceived(false) on the new OBMol.

  Stereochemistry is only preserved if the corresponding elements are wholly present in
  the substructure. For example, all four atoms and bonds of a tetrahedral stereocenter
  must be copied.

  Residue information is preserved if the original OBMol is marked as having
  its residues perceived. If this is not desired, either call
  OBMol::SetChainsPerceived(false) in advance on the original OBMol to avoid copying
  the residues (and then reset it afterwards), or else call it on the new OBMol so
  that residue information will be reperceived (when requested).

  Here is an example of using this method to copy ring systems to a new molecule.
  Given the molecule represented by the SMILES string, "FC1CC1c2ccccc2I", we will
  end up with a new molecule represented by the SMILES string, "C1CC1.c2ccccc2".
  \code{.cpp}
  OBBitVec atoms(mol.NumAtoms() + 1); // the maximum size needed
  FOR_ATOMS_OF_MOL(atom, mol) {
    if(atom->IsInRing())
      atoms.SetBitOn(atom->Idx());
  }
  OBBitVec excludebonds(mol.NumBonds()); // the maximum size needed
  FOR_BONDS_OF_MOL(bond, mol) {
    if(!bond->IsInRing())
      excludebonds.SetBitOn(bond->Idx());
  }
  OBMol newmol;
  mol.CopySubstructure(&newmol, &atoms, &excludebonds);
  \endcode

  When used from Python, note that "None" may be used to specify an empty value for
  the \p excludebonds parameter.

  \remark Some alternatives to using this function, which may be preferred in some
          instances due to efficiency or convenience are:
          -# Copying the entire OBMol, and then deleting the unwanted parts
          -# Modifiying the original OBMol, and then restoring it
          -# Using the SMILES writer option -xf to specify fragment atom idxs

  \return A boolean indicating success or failure. Currently failure is only reported
          if one of the specified atoms is not present, or \p atoms is a NULL
          pointer.

  \param newmol   The molecule to which to add the substructure. Note that atoms are
                  appended to this molecule.
  \param atoms    An OBBitVec, indexed by atom Idx, specifying which atoms to copy
  \param excludebonds  An OBBitVec, indexed by bond Idx, specifying a list of bonds
                       to exclude. By default, all bonds between the specified atoms are
                       included - this parameter overrides that.
  \param correctvalence  A value of 0, 1 (default) or 2 that indicates how atoms with missing
                         bonds are handled:
                        0 - Leave the implicit hydrogen count unchanged;
                        1 - Adjust the implicit hydrogen count to correct for
                            the missing bonds;
                        2 - Replace the missing bonds with bonds to dummy atoms
  \param atomorder Record the Idxs of the original atoms. That is, the first element
                   in this vector will be the Idx of the atom in the original OBMol
                   that corresponds to the first atom in the new OBMol. Note that
                   the information is appended to this vector.
  \param bondorder Record the Idxs of the original bonds. See \p atomorder above.

  **/

  bool OBMol::CopySubstructure(OBMol& newmol, OBBitVec *atoms, OBBitVec *excludebonds, unsigned int correctvalence,
                               std::vector<unsigned int> *atomorder, std::vector<unsigned int> *bondorder)
  {
    if (!atoms)
      return false;

    bool record_atomorder = atomorder != (std::vector<unsigned int>*)0;
    bool record_bondorder = bondorder != (std::vector<unsigned int>*)0;
    bool bonds_specified = excludebonds != (OBBitVec*)0;

    newmol.SetDimension(GetDimension());
    // If the parent had aromaticity perceived, then retain that for the fragment
    newmol.SetFlag(_flags & OB_AROMATIC_MOL);
    // The fragment will preserve the "chains perceived" flag of the parent
    newmol.SetFlag(_flags & OB_CHAINS_MOL);
    // We will check for residues only if the parent has chains perceived already
    bool checkresidues = HasChainsPerceived();

    // Now add the atoms
    map<OBAtom*, OBAtom*> AtomMap;//key is from old mol; value from new mol
    for (int bit = atoms->FirstBit(); bit != atoms->EndBit(); bit = atoms->NextBit(bit)) {
      OBAtom* atom = this->GetAtom(bit);
      if (!atom)
        return false;
      newmol.AddAtom(*atom); // each subsequent atom
      if (record_atomorder)
        atomorder->push_back(bit);
      AtomMap[&*atom] = newmol.GetAtom(newmol.NumAtoms());
    }

    //Add the residues
    if (checkresidues) {
      map<OBResidue*, OBResidue*> ResidueMap; // map from old->new
      for (int bit = atoms->FirstBit(); bit != atoms->EndBit(); bit = atoms->NextBit(bit)) {
        OBAtom* atom = this->GetAtom(bit);
        OBResidue* res = atom->GetResidue();
        if (!res) continue;
        map<OBResidue*, OBResidue*>::iterator mit = ResidueMap.find(res);
        OBResidue *newres;
        if (mit == ResidueMap.end()) {
          newres = newmol.NewResidue();
          *newres = *res;
          ResidueMap[res] = newres;
        } else {
          newres = mit->second;
        }
        OBAtom* newatom = AtomMap[&*atom];
        newres->AddAtom(newatom);
        newres->SetAtomID(newatom, res->GetAtomID(atom));
        newres->SetHetAtom(newatom, res->IsHetAtom(atom));
        newres->SetSerialNum(newatom, res->GetSerialNum(atom));
      }
    }

    // Update Stereo
    std::vector<OBGenericData*>::iterator data;
    std::vector<OBGenericData*> stereoData = GetAllData(OBGenericDataType::StereoData);
    for (data = stereoData.begin(); data != stereoData.end(); ++data) {
      if (static_cast<OBStereoBase*>(*data)->GetType() == OBStereo::CisTrans) {
        OBCisTransStereo *ct = dynamic_cast<OBCisTransStereo*>(*data);

        // Check that the entirety of this cistrans cfg occurs in this substructure
        OBCisTransStereo::Config cfg = ct->GetConfig();
        OBAtom* begin = GetAtomById(cfg.begin);
        if (AtomMap.find(begin) == AtomMap.end())
          continue;
        OBAtom* end = GetAtomById(cfg.end);
        if (AtomMap.find(end) == AtomMap.end())
          continue;
        bool skip_cfg = false;
        if (bonds_specified) {
          FOR_BONDS_OF_ATOM(bond, begin) {
            if (excludebonds->BitIsSet(bond->GetIdx())) {
              skip_cfg = true;
              break;
            }
          }
          if (skip_cfg)
            continue;
          FOR_BONDS_OF_ATOM(bond, end) {
            if (excludebonds->BitIsSet(bond->GetIdx())) {
              skip_cfg = true;
              break;
            }
          }
          if (skip_cfg)
            continue;
        }
        for (OBStereo::RefIter ri = cfg.refs.begin(); ri != cfg.refs.end(); ++ri) {
          if (*ri != OBStereo::ImplicitRef && AtomMap.find(GetAtomById(*ri)) == AtomMap.end()) {
            skip_cfg = true;
            break;
          }
        }
        if (skip_cfg)
          continue;

        OBCisTransStereo::Config newcfg;
        newcfg.specified = cfg.specified;
        newcfg.begin = cfg.begin == OBStereo::ImplicitRef ? OBStereo::ImplicitRef : AtomMap[GetAtomById(cfg.begin)]->GetId();
        newcfg.end = cfg.end == OBStereo::ImplicitRef ? OBStereo::ImplicitRef : AtomMap[GetAtomById(cfg.end)]->GetId();
        OBStereo::Refs refs;
        for (OBStereo::RefIter ri = cfg.refs.begin(); ri != cfg.refs.end(); ++ri) {
          OBStereo::Ref ref = *ri == OBStereo::ImplicitRef ? OBStereo::ImplicitRef : AtomMap[GetAtomById(*ri)]->GetId();
          refs.push_back(ref);
        }
        newcfg.refs = refs;

        OBCisTransStereo *newct = new OBCisTransStereo(this);
        newct->SetConfig(newcfg);
        newmol.SetData(newct);
      }
      else if (static_cast<OBStereoBase*>(*data)->GetType() == OBStereo::Tetrahedral) {
        OBTetrahedralStereo *tet = dynamic_cast<OBTetrahedralStereo*>(*data);
        OBTetrahedralStereo::Config cfg = tet->GetConfig();

        // Check that the entirety of this tet cfg occurs in this substructure
        OBAtom *center = GetAtomById(cfg.center);
        std::map<OBAtom*, OBAtom*>::iterator centerit = AtomMap.find(center);
        if (centerit == AtomMap.end())
          continue;
        if (cfg.from != OBStereo::ImplicitRef && AtomMap.find(GetAtomById(cfg.from)) == AtomMap.end())
          continue;
        bool skip_cfg = false;
        if (bonds_specified) {
          FOR_BONDS_OF_ATOM(bond, center) {
            if (excludebonds->BitIsSet(bond->GetIdx())) {
              skip_cfg = true;
              break;
            }
          }
          if (skip_cfg)
            continue;
        }
        for (OBStereo::RefIter ri = cfg.refs.begin(); ri != cfg.refs.end(); ++ri) {
          if (*ri != OBStereo::ImplicitRef && AtomMap.find(GetAtomById(*ri)) == AtomMap.end()) {
            skip_cfg = true;
            break;
          }
        }
        if (skip_cfg)
          continue;

        OBTetrahedralStereo::Config newcfg;
        newcfg.specified = cfg.specified;
        newcfg.center = centerit->second->GetId();
        newcfg.from = cfg.from == OBStereo::ImplicitRef ? OBStereo::ImplicitRef : AtomMap[GetAtomById(cfg.from)]->GetId();
        OBStereo::Refs refs;
        for (OBStereo::RefIter ri = cfg.refs.begin(); ri != cfg.refs.end(); ++ri) {
          OBStereo::Ref ref = *ri == OBStereo::ImplicitRef ? OBStereo::ImplicitRef : AtomMap[GetAtomById(*ri)]->GetId();
          refs.push_back(ref);
        }
        newcfg.refs = refs;

        OBTetrahedralStereo *newtet = new OBTetrahedralStereo(this);
        newtet->SetConfig(newcfg);
        newmol.SetData(newtet);
      }
    }

    // Options:
    // 1. Bonds that do not connect atoms in the subset are ignored
    // 2. As 1. but implicit Hs are added to replace them
    // 3. As 1. but asterisks are added to replace them
    FOR_BONDS_OF_MOL(bond, this) {
      bool skipping_bond = bonds_specified && excludebonds->BitIsSet(bond->GetIdx());
      map<OBAtom*, OBAtom*>::iterator posB = AtomMap.find(bond->GetBeginAtom());
      map<OBAtom*, OBAtom*>::iterator posE = AtomMap.find(bond->GetEndAtom());
      if (posB == AtomMap.end() && posE == AtomMap.end())
        continue;

      if (posB == AtomMap.end() || posE == AtomMap.end() || skipping_bond) {
        switch(correctvalence) {
        case 1:
          if (posB == AtomMap.end() || (skipping_bond && posE != AtomMap.end()))
            posE->second->SetImplicitHCount(posE->second->GetImplicitHCount() + bond->GetBondOrder());
          if (posE == AtomMap.end() || (skipping_bond && posB != AtomMap.end()))
            posB->second->SetImplicitHCount(posB->second->GetImplicitHCount() + bond->GetBondOrder());
          break;
        case 2: {
            OBAtom *atomB, *atomE;
            if (skipping_bond) {
              for(int N=0; N<2; ++N) {
                if (N==0) {
                  if (posB != AtomMap.end()) {
                    atomB = posB->second;
                    atomE = newmol.NewAtom();
                    if (record_atomorder)
                      atomorder->push_back(bond->GetEndAtomIdx());
                  }
                } else if (posE != AtomMap.end()) {
                  atomE = posE->second;
                  atomB = newmol.NewAtom();
                  if (record_atomorder)
                    atomorder->push_back(bond->GetBeginAtomIdx());
                }
                newmol.AddBond(atomB->GetIdx(), atomE->GetIdx(),
                  bond->GetBondOrder(), bond->GetFlags());
                if (record_bondorder)
                  bondorder->push_back(bond->GetIdx());
              }
            }
            else {
              atomB = (posB == AtomMap.end()) ? newmol.NewAtom() : posB->second;
              atomE = (posE == AtomMap.end()) ? newmol.NewAtom() : posE->second;
              if (record_atomorder) {
                if (posB == AtomMap.end())
                  atomorder->push_back(bond->GetBeginAtomIdx());
                else
                  atomorder->push_back(bond->GetEndAtomIdx());
              }
              newmol.AddBond(atomB->GetIdx(), atomE->GetIdx(),
                bond->GetBondOrder(), bond->GetFlags());
              if (record_bondorder)
                bondorder->push_back(bond->GetIdx());
            }
          }
          break;
        default:
          break;
        }
      }
      else {
        newmol.AddBond((posB->second)->GetIdx(), posE->second->GetIdx(),
                       bond->GetBondOrder(), bond->GetFlags());
        if (record_bondorder)
          bondorder->push_back(bond->GetIdx());
      }
    }

    return true;
  }

  bool OBMol::GetNextFragment( OBMolAtomDFSIter& iter, OBMol& newmol ) {
    if( ! iter ) return false;

    // We want to keep the atoms in their original order rather than use
    // the DFS order so just record the information first
    OBBitVec infragment(this->NumAtoms()+1);
    do { //for each atom in fragment
      infragment.SetBitOn(iter->GetIdx());
    } while ((iter++).next());

    bool ok = CopySubstructure(newmol, &infragment);

    return ok;
  }

  // Put the specified molecular charge on a single atom (which is expected for InChIFormat).
  // Assumes all the hydrogen is explicitly included in the molecule,
  // and that SetTotalCharge() has not been called. (This function is an alternative.)
  // Returns false if cannot assign all the charge.
  // Not robust in the general case, but see below for the more common simpler cases.
  bool OBMol::AssignTotalChargeToAtoms(int charge)
  {
    int extraCharge = charge - GetTotalCharge(); //GetTotalCharge() gets charge on atoms

    FOR_ATOMS_OF_MOL (atom, this)
    {
      unsigned int atomicnum = atom->GetAtomicNum();
      if (atomicnum == 1)
        continue;
      int charge = atom->GetFormalCharge();
      unsigned bosum = atom->GetExplicitValence();
      unsigned int totalValence = bosum + atom->GetImplicitHCount();
      unsigned int typicalValence = GetTypicalValence(atomicnum, bosum, charge);
      int diff = typicalValence - totalValence;
      if(diff != 0)
      {
        int c;
        if(extraCharge == 0)
          c = diff > 0 ? -1 : +1; //e.g. CH3C(=O)O, NH4 respectively
        else
          c = extraCharge < 0 ? -1 : 1;
        if (totalValence == GetTypicalValence(atomicnum, bosum, charge + c)) {
          atom->SetFormalCharge(charge + c);
          extraCharge -= c;
        }
      }
    }
    if(extraCharge != 0)
    {
      obErrorLog.ThrowError(__FUNCTION__, "Unable to assign all the charge to atoms", obWarning);
      return false;
    }
    return true;
 }
  /* These cases work ok
   original      charge  result
  [NH4]             +1   [NH4+]
  -C(=O)[O]         -1   -C(=O)[O-]
  -[CH2]            +1   -C[CH2+]
  -[CH2]            -1   -C[CH2-]
  [NH3]CC(=O)[O]     0   [NH3+]CC(=O)[O-]
  S(=O)(=O)([O])[O] -2   S(=O)(=O)([O-])[O-]
  [NH4].[Cl]         0   [NH4+].[Cl-]
  */

} // end namespace OpenBabel

//! \file mol.cpp
//! \brief Handle molecules. Implementation of OBMol.

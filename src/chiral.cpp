/*********************************************************************
chiral.cpp - Detect chiral atoms and molecules.

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
#include <list>

#include <openbabel/mol.h>
#include <openbabel/chiral.h>
#include <openbabel/math/matrix3x3.h>

using namespace std;
namespace OpenBabel
{
  // DEPRECATED
  // Seems to make a vector chirality become filled with array of +/- 1 for chiral atoms.
  void GetChirality(OBMol &mol, std::vector<int> &chirality)
  {
    chirality.resize(mol.NumAtoms()+1);
    fill(chirality.begin(),chirality.end(),0);

    OBAtom *atom;
    vector<OBAtom*>::iterator i;
    for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
      if (atom->IsChiral())
        {
          if (!atom->HasChiralVolume())
            {
              double sv = CalcSignedVolume(mol,atom);
              if (sv < 0.0)
                {
                  chirality[atom->GetIdx()-1] = -1;
                  atom->SetNegativeStereo();
                }
              else if (sv > 0.0)
                {
                  chirality[atom->GetIdx()-1] = 1;
                  atom->SetPositiveStereo();
                }
            }
          else // already calculated signed volume (e.g., imported from somewhere)
            {
              if (atom ->IsPositiveStereo())
                chirality[atom->GetIdx()-1] = 1;
              else
                chirality[atom->GetIdx()-1] = -1;
            }
        }
  }

  // DEPRECATED
  int GetParity4Ref(vector<unsigned int> pref)
  {
    if(pref.size()!=4)return(-1); // should be given a vector of size 4.
    int parity=0;
    for (int i=0;i<3;++i) // do the bubble sort this many times
      {
        for(int j=0;j<3;++j) // iterate across the array 4th element has no
          {                    // right hand neighbour so no need to sort
            if (pref[j+1] < pref[j]) // compare the two neighbors
              {
                unsigned int tmp = pref[j];        // swap a[j] and a[j+1]
                pref[j] = pref[j+1];
                pref[j+1] = tmp;
                parity++; // parity odd will invert stereochem
              }
          }
      } // End Bubble Sort
    return(parity%2);
  }

  // DEPRECATED
  bool CorrectChirality(OBMol &, OBAtom *atm, atomreftype i, atomreftype o)
  {
    if (!atm->HasChiralitySpecified()) // if no chirality defined can't do any more for 0D
      return(false);

    int parityI=0,parityO=0;
    OBChiralData* cd=(OBChiralData*)atm->GetData(OBGenericDataType::ChiralData);
    if ((cd->GetAtom4Refs(input)).size()!=4)return(false); // must have 4 refs
    parityI=GetParity4Ref(cd->GetAtom4Refs(i)); // Gets Atom4Refs used to define the chirality
    parityO=GetParity4Ref(cd->GetAtom4Refs(o));//GetsOutput parity.
    /* switch (CHTYPE)
       {
       case SMILES: // SMILES always uses 1234 atom refs
       parityO=0; // if parityO==parityI then clockwise in input = clockwise in output
       break;
       case MOLV3000: // MOLV3000 uses 1234 unless an H then 123H
       if (atm->GetHvyValence()==3)
       {
       OBAtom *nbr;
       int Hid=1000;// max Atom ID +1 should be used here
       vector<int> nbr_atms;
       vector<OBBond*>::iterator i;
       for (nbr = atm->BeginNbrAtom(i);nbr;nbr = atm->NextNbrAtom(i))
       {
       if (nbr->IsHydrogen()){Hid=nbr->GetIdx();continue;}
       nbr_atms.push_back(nbr->GetIdx());
       }
       sort(nbr_atms.begin(),nbr_atms.end());
       nbr_atms.push_back(Hid);
       int tmp[4];
       for(int i=0;i<4;++i){tmp[i]=nbr_atms[i];}
       parityO=GetParity4Ref(tmp);
       }
       else if (atm->GetHvyValence()==4)
       parityO=0;
       break;
       default:
       parityO=0;
       }*/
    if (parityO==parityI)
      {//cout << "Parity is the same"<<endl;
        return(true);
      }
    else if(parityO!=parityI) // Need to invert the Chirality which has been set
      { //cout << "Parity is Opposite"<<endl;
        if (atm->IsClockwise())
          {atm->UnsetStereo();atm->SetAntiClockwiseStereo();}
        else if (atm->IsAntiClockwise())
          {atm->UnsetStereo();atm->SetClockwiseStereo();}
        else
          return(false);
        return(true);
      }
    return false;
  }

  // DEPRECATED
  //! Calculate the signed volume for an atom.  If the atom has a valence of 3
  //! the coordinates of an attached hydrogen are calculated
  //! Puts attached Hydrogen last at the moment, like mol V3000 format.
  //! If ReZero=false (the default is true) always make pseudo z coords and leave them in mol
  double CalcSignedVolume(OBMol &mol,OBAtom *atm, bool ReZeroZ)
  {
    vector3 tmp_crd;
    vector<unsigned int> nbr_atms;
    vector<vector3> nbr_crds;
    bool use_central_atom = false,is2D=false;
    //   double hbrad = etab.CorrectedBondRad(1,0);

    if (!ReZeroZ || !mol.Has3D()) //give pseudo Z coords if mol is 2D
      {
        vector3 v,vz(0.0,0.0,1.0);
        is2D = true;
        OBAtom *nbr;
        OBBond *bond;
        vector<OBBond*>::iterator i;
        for (bond = atm->BeginBond(i);bond;bond = atm->NextBond(i))
          {
            nbr = bond->GetEndAtom();
            if (nbr != atm)
              {
                v = nbr->GetVector();
                if (bond->IsWedge())
                  v += vz;
                else
                  if (bond->IsHash())
                    v -= vz;

                nbr->SetVector(v);
              }
            else
              {
                nbr = bond->GetBeginAtom();
                v = nbr->GetVector();
                if (bond->IsWedge())
                  v -= vz;
                else
                  if (bond->IsHash())
                    v += vz;

                nbr->SetVector(v);
              }
          }
      }

    if (atm->GetHvyValence() < 3)
      {
        stringstream errorMsg;
        errorMsg << "Cannot calculate a signed volume for an atom with a heavy atom valence of " << atm->GetHvyValence() << endl;
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
        return(0.0);
      }

    // Create a vector with the coordinates of the neighbor atoms
    // Also make a vector with Atom IDs
    OBAtom *nbr;
    vector<OBBond*>::iterator bint;
    for (nbr = atm->BeginNbrAtom(bint);nbr;nbr = atm->NextNbrAtom(bint))
      {
        nbr_atms.push_back(nbr->GetIdx());
      }
    // sort the neighbor atoms to insure a consistent ordering
    sort(nbr_atms.begin(),nbr_atms.end());
    for (unsigned int i = 0; i < nbr_atms.size(); ++i)
      {
        OBAtom *tmp_atm = mol.GetAtom(nbr_atms[i]);
        nbr_crds.push_back(tmp_atm->GetVector());
      }
    /*
    // If we have three heavy atoms we need to calculate the position of the fourth
    if (atm->GetHvyValence() == 3)
    {
    double bondlen = hbrad+etab.CorrectedBondRad(atm->GetAtomicNum(),atm->GetHyb());
    atm->GetNewBondVector(tmp_crd,bondlen);
    nbr_crds.push_back(tmp_crd);
    }
    */
    for(unsigned int j=0;j < nbr_crds.size();++j) // Checks for a neighbour having 0 co-ords (added hydrogen etc)
      {
        // are the coordinates zero to 6 or more significant figures
        if (nbr_crds[j].IsApprox(VZero, 1.0e-6) && use_central_atom==false)
          use_central_atom=true;
        else if (nbr_crds[j].IsApprox(VZero, 1.0e-6))
          {
            obErrorLog.ThrowError(__FUNCTION__, "More than 2 neighbours have 0 co-ords when attempting 3D chiral calculation", obInfo);
          }
      }

    // If we have three heavy atoms we can use the chiral center atom itself for the fourth
    // will always give same sign (for tetrahedron), magnitude will be smaller.
    if(nbr_atms.size()==3 || use_central_atom==true)
      {
        nbr_crds.push_back(atm->GetVector());
        nbr_atms.push_back(mol.NumAtoms()+1); // meed to add largest number on end to work
      }
    OBChiralData* cd=(OBChiralData*)atm->GetData(OBGenericDataType::ChiralData); //Set the output atom4refs to the ones used
    if(cd==NULL)
      {
        cd = new OBChiralData;
        cd->SetOrigin(perceived);
        atm->SetData(cd);
      }
    cd->SetAtom4Refs(nbr_atms,calcvolume);

    //re-zero pseudo-coords
    if (is2D && ReZeroZ)
      {
        vector3 v;
        OBAtom *atom;
        vector<OBAtom*>::iterator k;
        for (atom = mol.BeginAtom(k);atom;atom = mol.NextAtom(k))
          {
            v = atom->GetVector();
            v.SetZ(0.0);
            atom->SetVector(v);
          }
      }

    return(signed_volume(nbr_crds[0],nbr_crds[1],nbr_crds[2],nbr_crds[3]));
  }

  // DEPRECATED
  //! Calculate a signed volume given a set of 4 coordinates
  double signed_volume(const vector3 &a, const vector3 &b, const vector3 &c, const vector3 &d)
  {
    vector3 A,B,C;
    A = b-a;
    B = c-a;
    C = d-a;
    matrix3x3 m(A,B,C);
    return(m.determinant());
  }

  //! \brief Calculate the Graph Potentials of a molecule
  //!
  //! based on
  //! V.E. and Rozenblit, A.B. Golender
  //! <em>Logical and Combinatorial Algorithms for Drug Design</em>. \n
  //! For an example see:
  //! Walters, W. P., Yalkowsky, S. H., \em JCICS, 1996, 36(5), 1015-1017.
  //! <a href="http://dx.doi.org/10.1021/ci950278o">DOI: 10.1021/ci950278o</a>
  void GraphPotentials(OBMol &mol, std::vector<double> &pot)
  {
    double det;

    vector<vector<double> > g,c,h;
    construct_g_matrix(mol,g);
    invert_matrix(g,det);
    construct_c_matrix(mol,c);
    mult_matrix(h,g,c);
    pot.resize(mol.NumAtoms()+1);

    for (unsigned int i = 0; i < mol.NumAtoms();++i)
      pot[i+1] = h[i][0];
  }


  //! Construct the matrix G, which puts each atoms valence+1
  //! on the diagonal and and -1 on the off diagonal if two
  //! atoms are connected.
  void construct_g_matrix(OBMol &mol, std::vector<std::vector<double> > &m)
  {
    unsigned int i,j;

    OBAtom *atm1,*atm2;
    vector<OBAtom*>::iterator aint,bint;

    m.resize(mol.NumAtoms());
    for (i = 0; i < m.size(); ++i)
      m[i].resize(mol.NumAtoms());

    for (atm1 = mol.BeginAtom(aint),i=0;atm1;atm1 = mol.NextAtom(aint),++i)
      for (atm2 = mol.BeginAtom(bint),j=0;atm2;atm2 = mol.NextAtom(bint),++j)
        {
          if (i == j)
            {
              m[i][j] = atm1->GetValence() + 1;
              m[i][j] += (double)atm1->GetAtomicNum()/10.0;
              m[i][j] += (double)atm1->GetHyb()/100.0;
            }
          else
            {
              if (atm1->IsConnected(atm2))
                m[i][j] = -1;
              else
                m[i][j] = 0;
            }
        }
  }

  //! Construct the matrix C, which is simply a column vector
  //! consisting of the valence for each atom
  void construct_c_matrix(OBMol &mol,std::vector<std::vector<double > > &m)
  {
    unsigned int i;
    OBAtom *atm1;
    vector<OBAtom*>::iterator aint;

    m.resize(mol.NumAtoms());
    for (i = 0; i < m.size(); ++i)
      m[i].resize(1);
    for (atm1 = mol.BeginAtom(aint),i=0;atm1;atm1 = mol.NextAtom(aint),++i)
      {
        m[i][0] = atm1->GetValence();
      }
  }

} // namespace OpenBabel

//! \file chiral.cpp
//! \brief Detect chiral atoms and molecules.

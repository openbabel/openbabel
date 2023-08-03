/**********************************************************************
molchrg.cpp - Assign Gasteiger partial charges.

Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2005 by Geoffrey R. Hutchison

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
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/generic.h>
#include <openbabel/oberror.h>
#include <openbabel/molchrg.h>
#include <openbabel/elements.h>

#include <cstring>

using namespace std;
namespace OpenBabel
{

  /*! \class OBGastChrg molchrg.h <openbabel/molchrg.h>
    \brief Assigns Gasteiger partial charges

    The OBGastChrg class is responsible for the assignment of partial
    charges to a molecule according to the Gasteiger charge model (sigma). When
    the partial charge of an atom is requested and partial charges do not
    yet exist for the molecule the OBGastChrg class will automatically be
    called to assign charges for the molecule. If charges have been read
    in for a molecule the following code shows how to force the
    recalculation of partial charges:
    \code
    OBMol mol;

    mol.UnsetPartialChargesPerceived();
    FOR_ATOMS_IN_MOL(atom, mol)
    {
       cout << "atom number = " << atom->GetIdx();
       cout << " charge = " << atom->GetPartialCharge() << endl;
    }
    \endcode
    Formal charges are used as seed values of the initial charge of atoms,
    and the partial charge is propagated to neighboring atoms. For
    example, quaternary amines would have a +1 charge, and the effect of
    the positive charge would be felt by neighboring atoms according to
    the Gasteiger model (sigma).

    For more information, see:
    J. Gasteiger & M. Marsili, "A New Model for Calculating Atomic Charges in Molecules" Tetrahedron Lett., (1978) 3181-3184.
  */

  bool OBGastChrg::AssignPartialCharges(OBMol &mol)
  {
    //InitialPartialCharges(mol);
    obErrorLog.ThrowError(__FUNCTION__,
                          "Ran OpenBabel::AssignPartialCharges", obAuditMsg);

    // Annotate that partial charges come from Gasteiger
    OBPairData *dp = new OBPairData;
    dp->SetAttribute("PartialCharges");
    dp->SetValue("Gasteiger");
    dp->SetOrigin(perceived);
    mol.SetData(dp);

    OBAtom *atom;
    vector<OBAtom*>::iterator i;

    GSVResize(mol.NumAtoms()+1);
    double a,b,c;
    for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
      {
        if (!GasteigerSigmaChi(atom,a,b,c))
          return(false);
        _gsv[atom->GetIdx()]->SetValues(a,b,c,atom->GetPartialCharge());
      }

    double alpha,charge,denom;
    unsigned j;
    int iter;
    OBBond *bond;
    OBAtom *src,*dst;
    vector<OBBond*>::iterator k;
    alpha = 1.0;
    for(iter = 0;iter < OB_GASTEIGER_ITERS;++iter)
      {
        alpha *= OB_GASTEIGER_DAMP;

        for( j=1;j < _gsv.size();++j)
          {
            charge = _gsv[j]->q;
            _gsv[j]->chi = (_gsv[j]->c*charge+_gsv[j]->b)*charge+_gsv[j]->a;
          }

        for (bond = mol.BeginBond(k);bond;bond = mol.NextBond(k))
          {
            src = bond->GetBeginAtom();
            dst = bond->GetEndAtom();

            if (_gsv[src->GetIdx()]->chi >= _gsv[dst->GetIdx()]->chi)
              {
                if (dst->GetAtomicNum() == OBElements::Hydrogen)
                  denom = double(OB_GASTEIGER_DENOM);
                else
                  denom = _gsv[dst->GetIdx()]->denom;
              }
            else
              {
                if (src->GetAtomicNum() == OBElements::Hydrogen)
                  denom = double(OB_GASTEIGER_DENOM);
                else
                  denom = _gsv[src->GetIdx()]->denom;
              }

            charge = (_gsv[src->GetIdx()]->chi - _gsv[dst->GetIdx()]->chi)/denom;
            _gsv[src->GetIdx()]->q -= alpha*charge;
            _gsv[dst->GetIdx()]->q += alpha*charge;
          }
      }

    for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
      atom->SetPartialCharge(_gsv[atom->GetIdx()]->q);

    return(true);
  }

  //! Set initial partial charges in @p mol
  //! Carbonyl O => -0.5
  //! Phosphate O => -0.666
  //! Sulfate O => -0.5
  //! All other atoms are set to have their initial charge from their formal charge
  void OBGastChrg::InitialPartialCharges(OBMol &mol)
  {
    OBAtom *atom;
    vector<OBAtom*>::iterator i;

    for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
      {
        if (atom->IsCarboxylOxygen())
          atom->SetPartialCharge(-0.500);
        else if (atom->IsPhosphateOxygen() &&
                 atom->GetHvyDegree() == 1)
          atom->SetPartialCharge(-0.666);
        else if (atom->IsSulfateOxygen())
          atom->SetPartialCharge(-0.500);
        else
          atom->SetPartialCharge((double)atom->GetFormalCharge());
      }
  }

  bool OBGastChrg::GasteigerSigmaChi(OBAtom *atom,double &a,double &b,double &c )
  {
    int count;
    double val[3] = {0.0,0.0,0.0};

    switch(atom->GetAtomicNum())
      {
      case 1: //H
        val[0] = 0.37;
        val[1] = 7.17;
        val[2] = 12.85;
        break;
      case 6: //C
        if (atom->GetHyb() == 3)
          {
            val[0] = 0.68;
            val[1] = 7.98;
            val[2] = 19.04;
          }
        if (atom->GetHyb() == 2)
          {
            val[0] = 0.98;
            val[1] = 8.79;
            val[2] = 19.62;
          }
        if (atom->GetHyb() == 1)
          {
            val[0] = 1.67;
            val[1] = 10.39;
            val[2] = 20.57;
          }
        break;
      case 7: //N
        if (atom->GetHyb() == 3)
          {
            if (atom->GetExplicitDegree() == 4 || atom->GetFormalCharge())
              {
                val[0] = 0.0;
                val[1] = 0.0;
                val[2] = 23.72;
              }
            else
              {
                val[0] = 2.08;
                val[1] = 11.54;
                val[2] = 23.72;
              }
          }

        if (atom->GetHyb() == 2)
          {
            if (EQ(atom->GetType(),"Npl") || EQ(atom->GetType(),"Nam"))
              {
                val[0] = 2.46;
                val[1] = 12.32;
                val[2] = 24.86;
              }
            else
              {
                val[0] = 2.57;
                val[1] = 12.87;
                val[2] = 24.87;
              }
          }

        if (atom->GetHyb() == 1)
          {
            val[0] = 3.71;
            val[1] = 15.68;
            val[2] = 27.11;
          }
        break;
      case 8: //O
        if (atom->GetHyb() == 3)
          {
            val[0] = 2.65;
            val[1] = 14.18;
            val[2] = 28.49;
          }
        if (atom->GetHyb() == 2)
          {
            val[0] = 3.75;
            val[1] = 17.07;
            val[2] = 31.33;
          }
        break;
      case 9: //F
        val[0] = 3.12;
        val[1] = 14.66;
        val[2] = 30.82;
        break;
      case 15: //P
        val[0] = 1.62;
        val[1] = 8.90;
        val[2] = 18.10;
        break;
      case 16: //S
        count = atom->CountFreeOxygens();
        if (count == 0 || count == 1)
          {
            val[0] = 2.39;
            val[1] = 10.14;
            val[2] = 20.65;
          }
        if (count > 1)
          {
            val[0] = 2.39;
            val[1] = 12.00;
            val[2] = 24.00;
          }
        /*S2? if (count == 0) {val[0] = 2.72;val[1] = 10.88;val[2] = 21.69;}*/
        break;
      case 17: //Cl
        val[0] = 2.66;
        val[1] = 11.00;
        val[2] = 22.04;
        break;
      case 35: //Br
        val[0] = 2.77;
        val[1] = 10.08;
        val[2] = 19.71;
        break;
      case 53: //I
        val[0] = 2.90;
        val[1] = 9.90;
        val[2] = 18.82;
        break;
      case 13: //Al
        val[0] = 1.06;
        val[1] = 5.47;
        val[2] = 11.65;
        break;
      }

    if (!IsNearZero(val[2]))
      {
        a = val[1];
        b = (val[2]-val[0])/2;
        c = (val[2]+val[0])/2 - val[1];
      }
    else
      return(false);

    return(true);
  }

  void OBGastChrg::GSVResize(int size)
  {
    vector <GasteigerState*>::iterator i;
    for (i = _gsv.begin();i != _gsv.end();++i)
      delete *i;
    _gsv.clear();

    for (int j = 0;j < size;++j)
      _gsv.push_back(new GasteigerState);
  }

  OBGastChrg::~OBGastChrg()
  {
    vector <GasteigerState*>::iterator i;
    for (i = _gsv.begin();i != _gsv.end();++i)
      delete *i;
  }

  GasteigerState::GasteigerState()
  {
    a = 0.0;
    b = 0.0;
    c = 0.0;
    denom = 0.0;
    chi = 0.0;
    q = 0.0;
  }

} // end namespace OpenBabel

//! \file molchrg.cpp
//! \brief Assign Gasteiger partial charges.

/**********************************************************************
forcefieldghemical.cpp - Ghemical force field.
 
Copyright (C) 2006 by Tim Vandermeersch <tim.vandermeersch@gmail.com>
 
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
#include <openbabel/mol.h>
#include "forcefieldghemical.h"

using namespace std;

namespace OpenBabel
{
  double OBForceFieldGhemical::E_Bond()
  {
    OBFFParameter *parameter;
    OBAtom *a, *b;
    double e, energy, l, l_ref, force, delta, delta2;
    
    energy = 0.0f;

    FOR_BONDS_OF_MOL(bond, _mol) {
      a = bond->GetBeginAtom();
      b = bond->GetEndAtom();	
      l = bond->GetLength();

      parameter = GetParameter(a->GetType(), b->GetType(), _ffbondparams);
      if (parameter == NULL) {
        parameter = GetParameter("FFFF", a->GetType(), _ffbondparams);

        if (parameter == NULL) {
          parameter = GetParameter("FFFF", b->GetType(), _ffbondparams);

	  if (parameter == NULL) {
            obErrorLog.ThrowError(__FUNCTION__, "Could not find all bond parameters", obError);
            exit(1);
	  }
	}
      }

      l_ref = parameter->dpar1;
      force = parameter->dpar2;

      delta = l - l_ref;
      delta2 = delta * delta;
      e = 0.5f * force * delta2;
      energy += e;

//      std::cout << a->GetType() << "-" << b->GetType() << "  l_ref=" << l_ref << "   l=" << l << "   e=" << e << std::endl;
    }

    return energy;
  }
  
  double OBForceFieldGhemical::E_Angle()
  {
    OBFFParameter *parameter;
    OBAtom *a, *b, *c;
    double e, energy, force, ang, ang_ref, delta, delta2;
    
    energy = 0.0f;

    FOR_ANGLES_OF_MOL(angle, _mol) {
      b = _mol.GetAtom((*angle)[0] + 1);
      a = _mol.GetAtom((*angle)[1] + 1);
      c = _mol.GetAtom((*angle)[2] + 1);
      ang = a->GetAngle(b->GetIdx(), c->GetIdx());

      parameter = GetParameter(a->GetType(), b->GetType(), c->GetType(), _ffangleparams);
      if (parameter == NULL) {
        parameter = GetParameter("FFFF", b->GetType(), c->GetType(), _ffangleparams);
	if (parameter == NULL) {
          parameter = GetParameter("FFFF", b->GetType(), "FFFF", _ffangleparams);
	  if (parameter == NULL) {
            obErrorLog.ThrowError(__FUNCTION__, "Could not find all angle parameters", obError);
            exit(1);
	  }
	}
      }
      
      ang_ref = parameter->dpar1;
      force   = parameter->dpar2;
	
      delta = ang - ang_ref;
      delta2 = delta * delta;
      e = 0.5f * force * delta2;
      energy += e;
    }
	
    return energy;
  }
 
  double OBForceFieldGhemical::E_Torsion() 
  {
    OBFFParameter *parameter;
    OBAtom *a, *b, *c, *d;
    OBBond *bc;
    double e, energy, tor, k, s, n, cosine;
 
    energy = 0.0f;

    FOR_TORSIONS_OF_MOL(t, _mol) {
      a = _mol.GetAtom((*t)[0] + 1);
      b = _mol.GetAtom((*t)[1] + 1);
      c = _mol.GetAtom((*t)[2] + 1);
      d = _mol.GetAtom((*t)[3] + 1);
      tor = _mol.GetTorsion(a->GetIdx(), b->GetIdx(), c->GetIdx(), d->GetIdx());
 
      parameter = GetParameter(a->GetType(), b->GetType(), c->GetType(), d->GetType(), _fftorsionparams);
      if (parameter == NULL) {
        parameter = GetParameter("FFFF", b->GetType(), c->GetType(), d->GetType(), _fftorsionparams);
        if (parameter == NULL) {
          parameter = GetParameter("FFFF", b->GetType(), c->GetType(), "FFFF", _fftorsionparams);
          if (parameter == NULL) {
            obErrorLog.ThrowError(__FUNCTION__, "Could not find all torsion parameters", obError);
            exit(1);
	  }
	}
      }
      
      k = parameter->dpar1;
      s = parameter->dpar2;
      
      bc = _mol.GetBond(b, c);
      n = 4.0f - bc->GetBondOrder();
/*
      if ((b->GetHyb() == 2) && (c->GetHyb() == 2))
        n = 2.0f;
      else if ((b->GetHyb() == 3) || (c->GetHyb() == 3))
        n = 3.0f;
      else 
        n = 1.0f;
*/
      cosine = cos(DEG_TO_RAD * (n * tor));
      e = 0.5f * k * (n + s * cosine);
      //e = 0.5f * k * (1.0f + s * cosine);
      //std::cout << "BO=" << bc->GetBondOrder()<< "   n=" << n << "   tor=" << tor << "  cosine=" << cosine << "   k=" << k << "   s=" << s << "   e=" << e << std::endl;
      energy += e;
    }

    return energy;
  }

  /*
  //  a
  //   \
  //    b---d      plane = a-b-c
  //   /
  //  c
  */
  double OBForceFieldGhemical::E_OOP() 
  {
    OBAtom *a, *b, *c, *d;
    double e, energy, force, angle, angle2;

    energy = 0.0f;
    force = 0.001f;

    FOR_ATOMS_OF_MOL(atom, _mol) {
      b = (OBAtom*) &*atom;

      if ((b->GetHyb() == 2) && (b->GetValence() == 3)) {
          a = NULL;
          c = NULL;
          d = NULL;

	  FOR_NBORS_OF_ATOM(nbr, b) {
	    if (a ==NULL)
	      a = (OBAtom*) &*nbr;
	    else if (c == NULL)
	      c = (OBAtom*) &*nbr;
	    else
	      d = (OBAtom*) &*nbr;
	  }
	  
            angle = PointPlaneAngle(a->GetVector(), b->GetVector(), c->GetVector(), d->GetVector());
	    angle2 = angle * angle;
	    e = 0.5f * force * angle2;
	    energy +=e;

            angle = PointPlaneAngle(d->GetVector(), b->GetVector(), c->GetVector(), a->GetVector());
	    angle2 = angle * angle;
	    e = 0.5f * force * angle2;
	    energy +=e;

            angle = PointPlaneAngle(a->GetVector(), b->GetVector(), d->GetVector(), c->GetVector());
	    angle2 = angle * angle;
	    e = 0.5f * force * angle2;
	    energy +=e;
      }
    }

    return energy;
  }
 
  double OBForceFieldGhemical::E_VDW()
  {
    OBFFParameter *parameter_a, *parameter_b;
    OBAtom *a, *b;
    double e, energy, ra, rb, rab, ka, kb, kab, term, term2, term4, term6, term12;

    energy = 0.0f;

    FOR_PAIRS_OF_MOL(p, _mol) {
      a = _mol.GetAtom((*p)[0]);
      b = _mol.GetAtom((*p)[1]);
      rab = a->GetDistance(b);

      parameter_a = GetParameter(a->GetType(), _ffvdwparams);
      parameter_b = GetParameter(b->GetType(), _ffvdwparams);
      if ((parameter_a == NULL) || (parameter_b == NULL)) {
        obErrorLog.ThrowError(__FUNCTION__, "Could not find all torsion parameters", obError);
        exit(1);
      }

      ra = parameter_a->dpar1;
      ka = parameter_a->dpar2;
      rb = parameter_b->dpar1;
      kb = parameter_b->dpar2;

      kab = sqrt(ka * kb);
      term = rab / (ra + rb);
      term2 = term * term;
      term4 = term2 * term2;
      term6 = term4 * term2;
      term12 = term6 * term6;
      e = kab * ((1.0f / term12) - (2.0f / term6));
      energy += e;
    }

    return energy;
  }


  //***********************************************
  //Make a global instance
  OBForceFieldGhemical theForceFieldGhemical("Ghemical",false);
  //***********************************************

  OBForceFieldGhemical::~OBForceFieldGhemical()
  {
  }

  OBForceFieldGhemical &OBForceFieldGhemical::operator=(OBForceFieldGhemical &src)
  {
    _mol = src._mol;
    return *this;
  }

  bool OBForceFieldGhemical::Setup(OBMol &mol)
  {
    _mol = mol;
    return SetTRPSTypes();
  }
 
  bool OBForceFieldGhemical::ParseParamFile()
  {
    vector<string> vs;
    char buffer[80];
    
    OBFFParameter parameter;

    // open data/ghemical.prm
    string buffer2, subbuffer;
    ifstream ifs1, ifs2, *ifsP;
    buffer2 = BABEL_DATADIR;
    buffer2 += FILE_SEP_CHAR;
    subbuffer = buffer2;
    subbuffer += BABEL_VERSION;
    subbuffer += FILE_SEP_CHAR;
    subbuffer += "ghemical.prm";
    buffer2 += "ghemical.prm";

    ifs1.open(subbuffer.c_str());
    ifsP= &ifs1;
    if (!(*ifsP))
    {
      ifs2.open(buffer2.c_str());
      ifsP = &ifs2;
    }

    while (ifsP->getline(buffer, 80)) {
      tokenize(vs, buffer);

      if (EQn(buffer, "bond", 4)) {
        parameter.clear();
        parameter._a = vs[1];
        parameter._b = vs[2];
        parameter.dpar1 = atof(vs[4].c_str()); // length
        parameter.dpar2 = atof(vs[5].c_str()); // force cte
        _ffbondparams.push_back(parameter);
      }
      if (EQn(buffer, "angle", 5)) {
        parameter.clear();
        parameter._a = vs[1];
        parameter._b = vs[2];
        parameter._c = vs[3];
        parameter.dpar1 = atof(vs[5].c_str()); // angle
        parameter.dpar2 = atof(vs[6].c_str()); // force cte
        _ffangleparams.push_back(parameter);
      }
      if (EQn(buffer, "torsion", 7)) {
        parameter.clear();
        parameter._a = vs[1];
        parameter._b = vs[2];
        parameter._c = vs[3];
        parameter._d = vs[4];
        parameter.dpar1 = atof(vs[6].c_str()); // force cte
        parameter.dpar2 = atof(vs[7].c_str()); // s
        _fftorsionparams.push_back(parameter);
      }
      if (EQn(buffer, "vdw", 3)) {
        parameter.clear();
        parameter._a = vs[1];
        parameter.dpar1 = atof(vs[2].c_str()); // r
        parameter.dpar2 = atof(vs[3].c_str()); // force cte
        _ffvdwparams.push_back(parameter);
      }
    }
	
    if (ifs1)
      ifs1.close();
    if (ifs2)
      ifs2.close();
 
    return 0;
  }
  
  bool OBForceFieldGhemical::SetTRPSTypes()
  {
    std::vector<std::vector<int> > _mlist; //!< match list for atom typing
    std::vector<std::pair<OBSmartsPattern*,std::string> > _vexttyp; //!< external atom type rules
    vector<vector<int> >::iterator j;
    vector<pair<OBSmartsPattern*,string> >::iterator i;
    OBSmartsPattern *sp;
    vector<string> vs;
    char buffer[80];
 
    _mol.SetAtomTypesPerceived();
    
    // open data/ghemical.prm
    string buffer2, subbuffer;
    ifstream ifs1, ifs2, *ifsP;
    buffer2 = BABEL_DATADIR;
    buffer2 += FILE_SEP_CHAR;
    subbuffer = buffer2;
    subbuffer += BABEL_VERSION;
    subbuffer += FILE_SEP_CHAR;
    subbuffer += "ghemical.prm";
    buffer2 += "ghemical.prm";

    ifs1.open(subbuffer.c_str());
    ifsP= &ifs1;
    if (!(*ifsP))
    {
      ifs2.open(buffer2.c_str());
      ifsP = &ifs2;
    }

    cout  << std::endl << "A T O M   T Y P E S" << std::endl << std::endl;
    
    while (ifsP->getline(buffer, 80)) {
      if (EQn(buffer, "atom", 4)) {
      	tokenize(vs, buffer);

        sp = new OBSmartsPattern;
        if (sp->Init(vs[1]))
          _vexttyp.push_back(pair<OBSmartsPattern*,string> (sp,vs[2]));
        else {
          delete sp;
          sp = NULL;
          obErrorLog.ThrowError(__FUNCTION__, " Could not parse EXTTYP line in atom type table from atomtyp.txt", obInfo);
          return false;
        }
        
	for (i = _vexttyp.begin();i != _vexttyp.end();++i) {
          if (i->first->Match(_mol)) {
            _mlist = i->first->GetMapList();
            for (j = _mlist.begin();j != _mlist.end();++j) {
              _mol.GetAtom((*j)[0])->SetType(i->second);
	    }
          }
        }
      }
    }

    cout << "IDX\tTYPE" << std::endl;
    FOR_ATOMS_OF_MOL (a, _mol)
      cout << a->GetIdx() << "\t" << a->GetType() << std::endl;
 
    if (ifs1)
      ifs1.close();
    if (ifs2)
      ifs2.close();
  }
  
  double OBForceFieldGhemical::Energy()
  {
    double energy;
    
    energy = E_Bond();
    energy += E_Angle();
    energy += E_Torsion();
    energy += E_OOP();
    energy += E_VDW();

    return energy;
  }

} // end namespace OpenBabel

//! \file forcefieldghemical.cpp
//! \brief Ghemical force field

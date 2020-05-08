/**********************************************************************
forcefieldmm2.cpp - MM2 force field.

Copyright (C) 2006-2007 by Tim Vandermeersch <tim.vandermeersch@gmail.com>

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
#include <openbabel/locale.h>
#include "forcefieldmm2.h"

using namespace std;

namespace OpenBabel
{
  double OBForceFieldMM2::E_Bond()
  {
    OBFFParameter *parameter;
    OBAtom *a, *b;
    vector3 va, vb, ab, vf;
    double e, energy, l, l_ref, force, delta, delta2, f;

    //sprintf(errbuf, "bondunit=%f  bond_cubic=%f  bond_quartic=%f\n", bondunit, bond_cubic, bond_quartic); // DEBUG

    energy = 0.0;
    if (forces.size() < _mol.NumAtoms() + 1)
      forces.resize(_mol.NumAtoms()+1);

    FOR_BONDS_OF_MOL(bond, _mol) {
      a = bond->GetBeginAtom();
      b = bond->GetEndAtom();
      l = bond->GetLength();

      va = a->GetVector();
      vb = b->GetVector();
      ab = va - vb;

      parameter = GetParameter(atoi(a->GetType()), atoi(b->GetType()), 0, 0, _ffbondparams);
      if (parameter == nullptr) {
        obErrorLog.ThrowError(__FUNCTION__, "Could not find all bond parameters", obError);
        return 0.0;
      }

      l_ref = parameter->_dpar[0];
      force = parameter->_dpar[1];

      delta = l - l_ref;
      delta2 = delta * delta;
      e = bondunit * force * delta2 * (1.0 + bond_cubic * delta + bond_quartic * delta2);
      energy += e;
      f = 2.0 * bondunit * force * delta * (1.0 + 1.5 * bond_cubic * delta + 2.0 * bond_quartic * delta2);
      vf = f * ab.normalize();
      forces[a->GetIdx()] -= vf;
      forces[b->GetIdx()] += vf;

      //sprintf(errbuf, "%sforce=%f  l_ref=%f  l=%f  delta=%f  E=%f (( ab(%f, %f, %f) vf(%f, %f, %f) ))\n", errbuf, force, l_ref, l, delta, e, ab.x(), ab.y(), ab.z(), vf.x(), vf.y(), vf.z()); // DEBUG
      //sprintf(errbuf, "%sab.length=%f\n", errbuf, ab.length()); // DEBUG
    }
    //sprintf(errbuf, "%sBond stretching: E=%f\n", errbuf, energy); // DEBUG
    //obErrorLog.ThrowError(__FUNCTION__, errbuf, obError); // DEBUG

    return energy;
  }

  double OBForceFieldMM2::E_Angle()
  {
    OBFFParameter *parameter;
    OBAtom *a, *b, *c;
    vector3 va, vb, vc, ab, bc, vfa, vfc, normal;
    double e, energy, force, ang, ang_ref, delta, delta2, delta3, delta4, f, dotprod, rab, rab2, rab3, rbc, rbc2, rbc3, root, dotprod_abbc;

    //sprintf(errbuf, "angleunit=%f  angle_sextic=%f\n", angleunit, angle_sextic); // DEBUG

    energy = 0.0;
    if (forces.size() < _mol.NumAtoms() + 1)
      forces.resize(_mol.NumAtoms()+1);

    FOR_ANGLES_OF_MOL(angle, _mol) {
      b = _mol.GetAtom((*angle)[0] + 1);
      a = _mol.GetAtom((*angle)[1] + 1);
      c = _mol.GetAtom((*angle)[2] + 1);
      ang = a->GetAngle(b->GetIdx(), c->GetIdx());

      va = a->GetVector();
      vb = b->GetVector();
      vc = c->GetVector();
      ab = va - vb;
      bc = vb - vc;

      parameter = GetParameter(atoi(a->GetType()), atoi(b->GetType()), atoi(c->GetType()), 0, _ffangleparams);
      if (parameter == nullptr) {
        obErrorLog.ThrowError(__FUNCTION__, "Could not find all angle parameters", obError);
        return 0.0;
      }

      ang_ref = parameter->_dpar[0];
      force   = parameter->_dpar[1];

      delta = ang - ang_ref;
      delta2 = delta * delta;
      delta3 = delta2 * delta;
      delta4 = delta2 * delta2;
      e = angleunit * force * delta2 * (1.0 + angle_sextic * delta4);
      energy += e;
      f = 2.0 * angleunit * force * delta * (1.0 + 3.0 * angle_sextic * delta4);
      dotprod = dot(ab, bc);
      rab = ab.length();
      rab2 = rab * rab;
      rab3 = rab2 * rab;
      rbc = bc.length();
      rbc2 = rbc * rbc;
      rbc3 = rbc2 * rbc;
      dotprod_abbc = dotprod / (rab * rbc);
      root = sqrt(1.0 - dotprod_abbc);
      vfa = rab2 * bc - dotprod * ab;
      vfa = f * vfa / (rab3 * rbc * root);
      vfc = dotprod * bc - rbc2 * ab;
      vfc = f * vfc / (rab * rbc3 * root);
      forces[a->GetIdx()] -= vfa;
      forces[b->GetIdx()] += (vfa + vfc);
      forces[c->GetIdx()] -= vfc;

      //sprintf(errbuf, "%sforce=%f  ang_ref=%f  ang=%f  delta=%f  E=%f   (%d-%d-%d)  (( vfa(%f, %f, %f)  vfc(%f, %f, %f) dotprod=%f))\n",
      //  errbuf, force, ang_ref, ang, delta, e, a->GetIdx(), b->GetIdx(), c->GetIdx(), vfa.x(), vfa.y(), vfa.z(), vfc.x(), vfc.z(), dotprod); // DEBUG
    }
    //sprintf(errbuf, "%sAngle bending: E=%f\n", errbuf, energy); // DEBUG
    //obErrorLog.ThrowError(__FUNCTION__, errbuf, obError); // DEBUG

    return energy;
  }

  double OBForceFieldMM2::E_StrBnd()
  {
    OBFFParameter *parameter;
    OBAtom *a, *b, *c;
    OBBond *b1, *b2;
    double e, energy, force, ang, ang_ref, l1, l2, l_ref1, l_ref2, delta_a, delta_b1, delta_b2;

    //sprintf(errbuf, "stretchbendunit=%f\n", stretchbendunit); // DEBUG

    energy = 0.0;

    FOR_ANGLES_OF_MOL(angle, _mol) {
      b = _mol.GetAtom((*angle)[0] + 1);
      a = _mol.GetAtom((*angle)[1] + 1);
      c = _mol.GetAtom((*angle)[2] + 1);
      ang = a->GetAngle(b->GetIdx(), c->GetIdx());
      b1 = _mol.GetBond(b->GetIdx(), a->GetIdx());
      b2 = _mol.GetBond(b->GetIdx(), c->GetIdx());
      l1 = b1->GetLength();
      l2 = b2->GetLength();

      parameter = GetParameter(atoi(b->GetType()), 0, 0, 0, _ffstretchbendparams);
      if (parameter == nullptr) {
        obErrorLog.ThrowError(__FUNCTION__, "Could not find all stretch-bend parameters", obError);
        return 0.0;
      }
      force = parameter->_dpar[0];

      parameter = GetParameter(atoi(a->GetType()), atoi(b->GetType()), atoi(c->GetType()), 0, _ffangleparams);
      ang_ref = parameter->_dpar[0];

      parameter = GetParameter(atoi(b->GetType()), atoi(a->GetType()), 0, 0, _ffbondparams);
      l_ref1 = parameter->_dpar[0];

      parameter = GetParameter(atoi(b->GetType()), atoi(c->GetType()), 0, 0, _ffbondparams);
      l_ref2 = parameter->_dpar[0];

      delta_a = ang - ang_ref;
      delta_b1 = l1 - l_ref1;
      delta_b2 = l2 - l_ref2;
      e = stretchbendunit * force * delta_a * (delta_b1 + delta_b2);
      energy += e;
      //sprintf(errbuf, "%sforce=%f  ang_ref=%f  ang=%f  l_ref1=%f  l1=%f  l_ref2=%f  l2=%f  E=%f\n", errbuf, force, ang_ref, ang, l_ref1, l1, l_ref2, l2, e); // DEBUG
    }
    //sprintf(errbuf, "%sBend-stretch term: E=%f\n", errbuf, energy); // DEBUG
    //obErrorLog.ThrowError(__FUNCTION__, errbuf, obError); // DEBUG

    return energy;
  }

  double OBForceFieldMM2::E_Torsion()
  {
    OBFFParameter *parameter;
    OBAtom *a, *b, *c, *d;
    double e, energy, tor, v1, v2, v3, cosine, cosine2, cosine3, sine, sine2, sine3, phi1, phi2, phi3, f, rc1, rc12, rc13, rc2, rc22, rc23, dot_c1c2, dot_rc1rc2, root;
    vector3 va, vb, vc, vd, ab, bc, cd, c1, c2, cross_bcc1, cross_bcc2, vfa, vfd;

    //sprintf(errbuf, "torsionunit=%f\n", torsionunit); // DEBUG

    energy = 0.0;

    FOR_TORSIONS_OF_MOL(t, _mol) {
      a = _mol.GetAtom((*t)[0] + 1);
      b = _mol.GetAtom((*t)[1] + 1);
      c = _mol.GetAtom((*t)[2] + 1);
      d = _mol.GetAtom((*t)[3] + 1);
      tor = _mol.GetTorsion(a->GetIdx(), b->GetIdx(), c->GetIdx(), d->GetIdx());

      va = a->GetVector();
      vb = b->GetVector();
      vc = c->GetVector();
      vd = d->GetVector();
      ab = va - vb;
      bc = vb - vc;
      cd = vc - vd;

      parameter = GetParameter(atoi(a->GetType()), atoi(b->GetType()), atoi(c->GetType()), atoi(d->GetType()), _fftorsionparams);
      if (parameter == nullptr) {
        obErrorLog.ThrowError(__FUNCTION__, "Could not find all torsion parameters", obError);
        return 0.0;
      }

      v1 = parameter->_dpar[0];
      v2 = parameter->_dpar[1];
      v3 = parameter->_dpar[2];

      cosine = cos(DEG_TO_RAD * tor);
      sine = sin(DEG_TO_RAD * tor);
      cosine2 = cos(DEG_TO_RAD * 2 * tor);
      sine2 = sin(DEG_TO_RAD * 2 * tor);
      cosine3 = cos(DEG_TO_RAD * 3 * tor);
      sine3 = sin(DEG_TO_RAD * 3 * tor);
      phi1 = 1.0 + cosine;
      phi2 = 1.0 - cosine2;
      phi3 = 1.0 + cosine3;
      e = torsionunit * (v1 * phi1 + v2 * phi2 + v3 * phi3);
      energy += e;
      f = -torsionunit * (- v1 * sine + 2 * v2 * sine2 - 3 * v3 * sine3);
      c1 = cross(ab, bc);
      c2 = cross(bc, cd);
      rc1 = c1.length();
      rc12 = rc1 * rc1;
      rc13 = rc12 * rc1;
      rc2 = c2.length();
      rc22 = rc2 * rc2;
      rc23 = rc22 * rc2;
      cross_bcc1 = cross(bc, c1);
      cross_bcc2 = cross(bc, c2);
      dot_c1c2 = dot(c1, c2);
      dot_rc1rc2 = dot_c1c2 / (rc1 * rc2);
      root = sqrt(1.0 - dot_rc1rc2);
      //        if (dot(b2,c3) > 0.0)
      //        torsion *= -1.0;

      vfa = rc12 * cross_bcc2 - dot_c1c2 * cross_bcc1;
      vfa = f * vfa / (rc13 * rc2 * root);
      vfd = rc22 * cross_bcc1 - dot_c1c2 * cross_bcc2;
      vfd = f * vfd / (rc1 * rc23 * root);

      forces[a->GetIdx()] -= vfa;
      forces[d->GetIdx()] -= vfd;


      //sprintf(errbuf, "%sv1=%f  v2=%f  v3=%f  tor=%f  cosine=%f  cosine2=%f  E=%f\n", errbuf, v1, v2, v3, tor, cosine, cosine2, e); // DEBUG
    }
    //sprintf(errbuf, "%sTorsions: E=%f\n", errbuf, energy); // DEBUG
    //obErrorLog.ThrowError(__FUNCTION__, errbuf, obError); // DEBUG

    return energy;
  }

  //
  //  a
  //   \
  //    b---d      plane = a-b-c
  //   /
  //  c
  //
  double OBForceFieldMM2::E_OOP()
  {
    OBAtom *a, *b, *c, *d;
    double e, energy, force, angle, angle2;
    //sprintf(errbuf, "outplanebendunit=%f\n", outplanebendunit); // DEBUG

    energy = 0.0;

    FOR_ATOMS_OF_MOL(atom, _mol) {
      b = (OBAtom*) &*atom;

      for (int idx=0; idx < _ffoutplanebendparams.size(); idx++) {
        if (atoi(b->GetType()) == _ffoutplanebendparams[idx].a) {
          a = nullptr;
          c = nullptr;
          d = nullptr;

          FOR_NBORS_OF_ATOM(nbr, b) {
            if (a ==nullptr)
              a = (OBAtom*) &*nbr;
            else if (c == nullptr)
              c = (OBAtom*) &*nbr;
            else
              d = (OBAtom*) &*nbr;
          }

          if (atoi(d->GetType()) == _ffoutplanebendparams[idx].b) {
            force = _ffoutplanebendparams[idx]._dpar[0];
            angle = Point2PlaneAngle(d->GetVector(), a->GetVector(), b->GetVector(), c->GetVector());
            //angle = PointPlaneAngle(a->GetVector(), b->GetVector(), c->GetVector(), d->GetVector());
            angle2 = angle * angle;
            e = outplanebendunit * force * angle2;
            energy +=e;

            //sprintf(errbuf, "%sout-of-plane: a-b-c-d  %d-%d-%d-%d  angle=%f  E=%f\n", errbuf, a->GetIdx(), b->GetIdx(), c->GetIdx(), d->GetIdx(), angle, e); // DEBUG
          }
          if (atoi(a->GetType()) == _ffoutplanebendparams[idx].b) {
            // a <-> d
            force = _ffoutplanebendparams[idx]._dpar[0];
            angle = Point2PlaneAngle(a->GetVector(), d->GetVector(), b->GetVector(), c->GetVector());
            //angle = PointPlaneAngle(d->GetVector(), b->GetVector(), c->GetVector(), a->GetVector());
            angle2 = angle * angle;
            e = outplanebendunit * force * angle2;
            energy +=e;

            //sprintf(errbuf, "%sout-of-plane: a-b-c-d  %d-%d-%d-%d  angle=%f  E=%f\n", errbuf, d->GetIdx(), b->GetIdx(), c->GetIdx(), a->GetIdx(), angle, e); // DEBUG
          }
          if (atoi(c->GetType()) == _ffoutplanebendparams[idx].b) {
            // c <-> d
            force = _ffoutplanebendparams[idx]._dpar[0];
            angle = Point2PlaneAngle(c->GetVector(), a->GetVector(), b->GetVector(), d->GetVector());
            //angle = PointPlaneAngle(a->GetVector(), b->GetVector(), d->GetVector(), c->GetVector());
            angle2 = angle * angle;
            e = outplanebendunit * force * angle2;
            energy +=e;

            //sprintf(errbuf, "%sout-of-plane: a-b-c-d  %d-%d-%d-%d  angle=%f  E=%f\n", errbuf, a->GetIdx(), b->GetIdx(), d->GetIdx(), c->GetIdx(), angle, e); // DEBUG
          }
        }
      }
    }
    //sprintf(errbuf, "%sOut-of-plane bending: E=%f\n", errbuf, energy); // DEBUG
    //obErrorLog.ThrowError(__FUNCTION__, errbuf, obError); // DEBUG

    return energy;
  }

  double OBForceFieldMM2::E_VDW()
  {
    OBFFParameter *parameter;
    OBAtom *a, *b;
    vector3 va, vb, ab, vf;
    double e, energy, ra, rb, rab, rr, rrab, rrab2, rrab4, rrab6, rrab7, abrr, eps, epsa, epsb, f;

    //sprintf(errbuf, "a_expterm=%f  b_expterm=%f c_expter=%fm\n", a_expterm, b_expterm, c_expterm); // DEBUG

    energy = 0.0;

    FOR_PAIRS_OF_MOL(p, _mol) {
      a = _mol.GetAtom((*p)[0]);
      b = _mol.GetAtom((*p)[1]);
      rab = a->GetDistance(b);

      va = a->GetVector();
      vb = b->GetVector();
      ab = vb - va;

      parameter = GetParameter(atoi(a->GetType()), atoi(b->GetType()), 0, 0, _ffvdwprparams);
      if (parameter != nullptr) {
        rr = parameter->_dpar[0];
        eps = parameter->_dpar[1];
      } else {
        parameter = GetParameter(atoi(a->GetType()), 0, 0, 0, _ffvdwparams);
        ra = parameter->_dpar[0];
        epsa = parameter->_dpar[1];
        parameter = GetParameter(atoi(b->GetType()), 0, 0, 0, _ffvdwparams);
        rb = parameter->_dpar[0];
        epsb = parameter->_dpar[1];

        //sprintf(errbuf, "%s\n\n     *** v1=%f  v2=%f  v3=%f\n\n", errbuf, v1, v2, v3); // DEBUG
        rr = ra + rb;
        eps = (epsa + epsb) / 2; // is this correct??
      }

      rrab = rr / rab;
      rrab2 = rrab * rrab;
      rrab4 = rrab2 * rrab2;
      rrab6 = rrab4 * rrab2;
      rrab7 = rrab6 * rrab;
      abrr = rab / rr;
      e = eps * (a_expterm * exp(-b_expterm * abrr) - c_expterm * rrab6);
      energy += e;
      f = eps * (a_expterm * b_expterm * exp(-b_expterm * abrr) - 6.0 * c_expterm * rrab7);
      vf = f * ab.normalize();
      //forces[a->GetIdx()] -= vf;
      //forces[b->GetIdx()] += vf;

      //sprintf(errbuf, "%sepsa=%f  epsb=%f epsilon=%f rab=%f  ra=%f  rb=%f  rr=%f  E=%f  f=%f\n", errbuf, epsa, epsb, eps, rab, ra, rb, rr, e, f); // DEBUG
    }
    //sprintf(errbuf, "%sVan der Waals interactions: E=%f\n", errbuf, energy); // DEBUG
    //obErrorLog.ThrowError(__FUNCTION__, errbuf, obError); // DEBUG

    return energy;
  }

  double OBForceFieldMM2::E_Electrostatic()
  {
    OBAtom *a, *b, *c, *d;
    double e, energy, dipole, dipole2, r2, ri2, rk2, rirkr3, dotp, doti, dotk, fik;
    double debye, electric, f;
    vector3 va, vb, vc, vd, ab, cd, q, r;
    int idx;

    //sprintf(errbuf, "dielectric=%f\n", dielectric); // DEBUG

    energy = 0.0;
    debye = 4.8033324;
    electric = 332.05382;
    f = electric / (debye * debye * dielectric);

    FOR_BONDS_OF_MOL(bond, _mol) {
      a = bond->GetBeginAtom();
      b = bond->GetEndAtom();

      idx = GetParameterIdx(atoi(a->GetType()), atoi(b->GetType()), 0, 0, _ffdipoleparams);
      if (idx > 0) {
        dipole = _ffdipoleparams[idx]._dpar[0];

        FOR_BONDS_OF_MOL(bond2, _mol) {
          if (&*bond == &*bond2)
            continue;
          if (bond2->GetIdx() < bond->GetIdx())
            continue;

          c = bond2->GetBeginAtom();
          d = bond2->GetEndAtom();

          idx = GetParameterIdx(atoi(c->GetType()), atoi(d->GetType()), 0, 0, _ffdipoleparams);
          if (idx > 0) {
            dipole2 = _ffdipoleparams[idx]._dpar[0];

            //
            //    dipole1      dipole2
            //       i-----------j           i = center of bond 1
            //        \ a1    a2/            j = center of bond 2
            //         \       /
            //          \     /
            //           \a3 /
            //            \ /
            //             k
            //
            va = a->GetVector();
            vb = b->GetVector();
            vc = c->GetVector();
            vd = d->GetVector();

            ab = va - vb;
            cd = vc - vd;
            ri2 = dot(ab, ab);
            rk2 = dot(cd, cd);
            q = (va + vb) / 2;
            r = (vc + vd) / 2;
            r2 = dot(r, r);
            rirkr3 = sqrt(ri2 * rk2 * r2) * r2;
            dotp = dot(ab, cd);
            doti = dot(ab, r);
            dotk = dot(cd, r);
            fik = f * dipole * dipole2;

            e = fik * (dotp - 3.0 * doti * dotk/r2) / rirkr3;
            energy += e;
            //sprintf(errbuf, "%sdipole=%f  dipole2=%f  r=%f cos_a1=%f  cos_a2=%f  cos_a3=%f  E=%f\n", errbuf, dipole, dipole2, r, cos_a1, cos_a2, cos_a3, e); // DEBUG
          }
        }
      }
    }
    //sprintf(errbuf, "%sDipole-dipole: E=%f\n", errbuf, energy); // DEBUG
    //obErrorLog.ThrowError(__FUNCTION__, errbuf, obError); // DEBUG

    return energy;
  }

  //
  // OBForceFieldMM2 member functions
  //
  //***********************************************
  //Make a global instance
  OBForceFieldMM2 theForceFieldMM2("MM2", false);
  //***********************************************

  OBForceFieldMM2::~OBForceFieldMM2()
  {
  }

  OBForceFieldMM2 &OBForceFieldMM2::operator=(OBForceFieldMM2 &src)
  {
    if (this != &src)
    {
      _mol = src._mol;
      _init = src._init;
      return *this;
    }
  }

  bool OBForceFieldMM2::Setup(OBMol &mol)
  {
    if (!_init) {
      ParseParamFile();
      _init = true;
    }

    _mol = mol;
    SetMM2Types();
    return true;
  }

  bool OBForceFieldMM2::ParseParamFile()
  {
    vector<string> vs;
    char buffer[80];

    OBFFParameter parameter;

    // open data/mm2.prm
    ifstream ifs;
    if (OpenDatafile(ifs, "mm2.prm").length() == 0) {
      obErrorLog.ThrowError(__FUNCTION__, "Cannot open mm2.prm", obError);
      return false;
    }

    // Set the locale for number parsing to avoid locale issues: PR#1785463
    obLocale.SetLocale();

    while (ifs.getline(buffer, 80)) {

      tokenize(vs, buffer);

      if (EQn(buffer, "bondunit", 8)) {
        bondunit = atof(vs[1].c_str());
        continue;
      }
      if (EQn(buffer, "bond-cubic", 10)) {
        bond_cubic = atof(vs[1].c_str());
        continue;
      }
      if (EQn(buffer, "bond-quartic", 12)) {
        bond_quartic = atof(vs[1].c_str());
        continue;
      }
      if (EQn(buffer, "angleunit", 9)) {
        angleunit = atof(vs[1].c_str());
        continue;
      }
      if (EQn(buffer, "angle-sextic", 12)) {
        angle_sextic = atof(vs[1].c_str());
        continue;
      }
      if (EQn(buffer, "strbndunit", 10)) {
        stretchbendunit = atof(vs[1].c_str());
        continue;
      }
      if (EQn(buffer, "torsionunit", 11)) {
        torsionunit = atof(vs[1].c_str());
        continue;
      }
      if (EQn(buffer, "opbendunit", 10)) {
        outplanebendunit = atof(vs[1].c_str());
        continue;
      }
      if (EQn(buffer, "a-expterm", 9)) {
        a_expterm = atof(vs[1].c_str());
        continue;
      }
      if (EQn(buffer, "b-expterm", 9)) {
        b_expterm = atof(vs[1].c_str());
        continue;
      }
      if (EQn(buffer, "c-expterm", 9)) {
        c_expterm = atof(vs[1].c_str());
        continue;
      }
      if (EQn(buffer, "vdwtype", 7))
        continue;
      if (EQn(buffer, "vdw-14-scale", 12))
        continue;
      if (EQn(buffer, "dielectric", 10)) {
        dielectric = atof(vs[1].c_str());
        continue;
      }

      parameter.clear();

      if (EQn(buffer, "bond3", 5))
        continue;
      if (EQn(buffer, "bond4", 5))
        continue;
      if (EQn(buffer, "bond", 4)) {
        parameter.a = atoi(vs[1].c_str());
        parameter.b = atoi(vs[2].c_str());
        parameter._dpar.push_back(atof(vs[4].c_str()));
        parameter._dpar.push_back(atof(vs[3].c_str()));

        _ffbondparams.push_back(parameter);
        continue;
      }
      if (EQn(buffer, "angle3", 6))
        continue;
      if (EQn(buffer, "angle4", 6))
        continue;
      if (EQn(buffer, "angle", 5)) {
        parameter.a = atoi(vs[1].c_str());
        parameter.b = atoi(vs[2].c_str());
        parameter.c = atoi(vs[3].c_str());
        parameter._dpar.push_back(atof(vs[5].c_str()));
        parameter._dpar.push_back(atof(vs[4].c_str()));
        _ffangleparams.push_back(parameter);
        continue;
      }
      if (EQn(buffer, "strbnd", 6)) {
        parameter.a = atoi(vs[1].c_str());
        parameter._dpar.push_back(atof(vs[2].c_str()));
        _ffstretchbendparams.push_back(parameter);
        continue;
      }
      if (EQn(buffer, "torsion4", 8))
        continue;
      if (EQn(buffer, "torsion", 7)) {
        parameter.a = atoi(vs[1].c_str());
        parameter.b = atoi(vs[2].c_str());
        parameter.c = atoi(vs[3].c_str());
        parameter.d = atoi(vs[4].c_str());
        parameter._dpar.push_back(atof(vs[5].c_str()));
        parameter._dpar.push_back(atof(vs[8].c_str()));
        parameter._dpar.push_back(atof(vs[11].c_str()));
        _fftorsionparams.push_back(parameter);
        continue;
      }
      if (EQn(buffer, "opbend", 6)) {
        parameter.a = atoi(vs[1].c_str());
        parameter.b = atoi(vs[2].c_str());
        parameter._dpar.push_back(atof(vs[3].c_str()));
        _ffoutplanebendparams.push_back(parameter);
        continue;
      }
      if (EQn(buffer, "vdwpr", 5)) {
        parameter.a = atoi(vs[1].c_str());
        parameter.b = atoi(vs[2].c_str());
        parameter._dpar.push_back(atof(vs[3].c_str()));
        parameter._dpar.push_back(atof(vs[4].c_str()));
        _ffvdwprparams.push_back(parameter);
        continue;
      }
      if (EQn(buffer, "vdw", 3)) {
        parameter.a = atoi(vs[1].c_str());
        parameter._dpar.push_back(atof(vs[2].c_str()));
        parameter._dpar.push_back(atof(vs[3].c_str()));
        _ffvdwparams.push_back(parameter);
        continue;
      }
      if (EQn(buffer, "dipole", 6)) {
        parameter.a = atoi(vs[1].c_str());
        parameter.b = atoi(vs[2].c_str());
        parameter._dpar.push_back(atof(vs[3].c_str()));
        parameter._dpar.push_back(atof(vs[4].c_str()));
        _ffdipoleparams.push_back(parameter);
        continue;
      }
    }

    if (ifs)
      ifs.close();

    // return the locale to the original one
    obLocale.RestoreLocale();

    return 0;
  }

  bool OBForceFieldMM2::SetMM2Types()
  {
    string atomtype;

    ttab.SetFromType("INT");
    ttab.SetToType("MM2");

    FOR_ATOMS_OF_MOL(atom, _mol) {
      ttab.Translate(atomtype, atom->GetType());
      atom->SetType(atomtype);
    }
    return true;
  }

  double OBForceFieldMM2::Energy()
  {
    double energy;

    energy = E_Bond();
    energy += E_Angle();
    energy += E_StrBnd();
    energy += E_Torsion();
    energy += E_OOP();
    energy += E_VDW();
    energy += E_Electrostatic();

    return energy;
  }

} // end namespace OpenBabel

//! \file forcefieldmm2.cpp
//! \brief MM2 force field

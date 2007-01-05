/**********************************************************************
forcefield.cpp - Handle OBForceField class.
 
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

#include <openbabel/forcefield.h>

#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obiter.h>

using namespace std;

namespace OpenBabel
{
  int OBForceField::GetParameterIdx(int a, int b, int c, int d, std::vector<OBFFParameter> &parameter)
  {
    if (!b)
      for (unsigned int idx=0; idx < parameter.size(); idx++)
        if (a == parameter[idx].a)
	  return idx;

    if (!c)
      for (unsigned int idx=0; idx < parameter.size(); idx++)
        if (((a == parameter[idx].a) && (b == parameter[idx].b)) || ((a == parameter[idx].b) && (b == parameter[idx].a)))
	  return idx;

    if (!d)
      for (unsigned int idx=0; idx < parameter.size(); idx++)
        if (((a == parameter[idx].a) && (b == parameter[idx].b) && (c == parameter[idx].c)) || 
	    ((a == parameter[idx].c) && (b == parameter[idx].b) && (c == parameter[idx].a)))
	  return idx;
 
    for (unsigned int idx=0; idx < parameter.size(); idx++)
      if (((a == parameter[idx].a) && (b == parameter[idx].b) && (c == parameter[idx].c) && (d == parameter[idx].d)) || 
          ((a == parameter[idx].d) && (b == parameter[idx].c) && (c == parameter[idx].b) && (d == parameter[idx].a)))
	return idx;

    return -1;
  }
  
  OBFFParameter* OBForceField::GetParameter(int a, int b, int c, int d, std::vector<OBFFParameter> &parameter)
  {
    OBFFParameter *par;

    if (!b)
      for (unsigned int idx=0; idx < parameter.size(); idx++)
        if (a == parameter[idx].a) {
	  par = &parameter[idx];
	  return par;
	}

    if (!c)
      for (unsigned int idx=0; idx < parameter.size(); idx++)
        if (((a == parameter[idx].a) && (b == parameter[idx].b)) || ((a == parameter[idx].b) && (b == parameter[idx].a))) {
	  par = &parameter[idx];
	  return par;
	}

    if (!d)
      for (unsigned int idx=0; idx < parameter.size(); idx++)
        if (((a == parameter[idx].a) && (b == parameter[idx].b) && (c == parameter[idx].c)) || 
	    ((a == parameter[idx].c) && (b == parameter[idx].b) && (c == parameter[idx].a))) {
	  par = &parameter[idx];
	  return par;
        }

    for (unsigned int idx=0; idx < parameter.size(); idx++)
      if (((a == parameter[idx].a) && (b == parameter[idx].b) && (c == parameter[idx].c) && (d == parameter[idx].d)) || 
          ((a == parameter[idx].d) && (b == parameter[idx].c) && (c == parameter[idx].b) && (d == parameter[idx].a))) {
	par = &parameter[idx];
	return par;
      }

    return NULL;
  }
  
  OBFFParameter* OBForceField::GetParameter(const char* a, const char* b, const char* c, const char* d, std::vector<OBFFParameter> &parameter)
  {
    OBFFParameter *par;
    if (a == NULL)
      return NULL;

    if (b == NULL) {
      std::string _a(a);
      for (unsigned int idx=0; idx < parameter.size(); idx++) 
        if (_a == parameter[idx]._a) {
	  par = &parameter[idx];
	  return par;
	}
      return NULL;
    }
    if (c == NULL) {
      std::string _a(a);
      std::string _b(b);
      for (unsigned int idx=0; idx < parameter.size(); idx++)
        if (((_a == parameter[idx]._a) && (_b == parameter[idx]._b)) || ((_a == parameter[idx]._b) && (_b == parameter[idx]._a))) {
	  par = &parameter[idx];
	  return par;
	}
      return NULL;
    }
    if (d == NULL) {
      std::string _a(a);
      std::string _b(b);
      std::string _c(c);
      for (unsigned int idx=0; idx < parameter.size(); idx++)
        if (((_a == parameter[idx]._a) && (_b == parameter[idx]._b) && (_c == parameter[idx]._c)) || 
	    ((_a == parameter[idx]._c) && (_b == parameter[idx]._b) && (_c == parameter[idx]._a))) {
	  par = &parameter[idx];
	  return par;
        }
      return NULL;
    }
    std::string _a(a);
    std::string _b(b);
    std::string _c(c);
    std::string _d(d);
    for (unsigned int idx=0; idx < parameter.size(); idx++)
      if (((_a == parameter[idx]._a) && (_b == parameter[idx]._b) && (_c == parameter[idx]._c) && (_d == parameter[idx]._d)) || 
          ((_a == parameter[idx]._d) && (_b == parameter[idx]._c) && (_c == parameter[idx]._b) && (_d == parameter[idx]._a))) {
	par = &parameter[idx];
	return par;
      }

    return NULL;
  }


  void OBForceField::SteepestDescent(int steps) 
  {
    double e_n1, e_n2;
    double h;
    //   double fmax;
    vector3 vf;
    std::vector<vector3> old_xyz;

    h = 0.1f;
    old_xyz.resize(_mol.NumAtoms()+1);

    e_n1 = Energy();

    for (int i=0; i<steps; i++) {
      //fmax = GetFmax();
      
      FOR_ATOMS_OF_MOL (a, _mol) {
        //vf = forces[a->GetIdx()];
	vf = NumericalDerivative(a->GetIdx());
        //forces[a->GetIdx()].Set(0.0f, 0.0f, 0.0f);
	//vf = vf / fmax * h;
	vf = vf.normalize();
	vf = vf * h;
	old_xyz[a->GetIdx()] = a->GetVector();
        a->SetVector(a->x() + vf.x(), a->y() + vf.y(), a->z() + vf.z());
        //sprintf(errbuf, "%svf=(%f, %f, %f)\n", errbuf, vf.x(), vf.y(), vf.z()); // DEBUG
      }
      e_n2 = Energy();

      std::cout << "e_n1 e=" << e_n1 << "  e_n2=" << e_n2  << ", h=" << h << std::endl;
      if (e_n2 >= e_n1) {
	h = 0.5f * h;
        FOR_ATOMS_OF_MOL (a, _mol)
          a->SetVector(old_xyz[a->GetIdx()].x(), old_xyz[a->GetIdx()].y(), old_xyz[a->GetIdx()].z());
      }
      if (e_n2 < e_n1) {
        e_n1 = e_n2;
	h = 1.2f * h;
      }
 
    }
    std::cout << std::endl;
  }
  /* 
  double OBForceField::GetFmax(OBMol &mol)
  {
    double fmax;
    std::vector<vector3>::iterator v; 

    for (v = forces.begin(); v != forces.end(); v++) {
      if (v->x() > fmax)
        fmax = v->x();
      if (v->y() > fmax)
        fmax = v->y();
      if (v->z() > fmax)
        fmax = v->z();
    }

    if (fmax < 0)
      fmax = -fmax;
    
    return fmax;
  }
*/
  vector3 OBForceField::NumericalDerivative(int a)
  {
    OBAtom *atom;
    vector3 va, grad;
    double e_orig, e_plus_delta, e_minus_delta, delta, dx, dy, dz;

    delta = 0.01f;

    atom = _mol.GetAtom(a);
    va = atom->GetVector();

    e_orig = Energy();
    atom->SetVector(va.x() + delta, va.y(), va.z());
    e_plus_delta = Energy();
    atom->SetVector(va.x() - delta, va.y(), va.z());
    e_minus_delta = Energy();

    dx = (e_minus_delta + e_plus_delta) / 2;

    if (e_plus_delta > e_minus_delta)
      dx = -dx;

    //e_orig = Energy();
    atom->SetVector(va.x(), va.y() + delta, va.z());
    e_plus_delta = Energy();
    atom->SetVector(va.x(), va.y() - delta, va.z());
    e_minus_delta = Energy();

    dy = (e_minus_delta + e_plus_delta) / 2;
    
    if (e_plus_delta > e_minus_delta)
      dy = -dy;

    //e_orig = Energy();
    atom->SetVector(va.x(), va.y(), va.z() + delta);
    e_plus_delta = Energy();
    atom->SetVector(va.x(), va.y(), va.z() - delta);
    e_minus_delta = Energy();

    dz = (e_minus_delta + e_plus_delta) / 2;

    if (e_plus_delta > e_minus_delta)
      dz = -dz;

    atom->SetVector(va.x(), va.y(), va.z());

    grad.Set(dx, dy, dz);

    //char errbuf[3600]; // DEBUG
    //sprintf(errbuf, "grad=(%f, %f, %f)  e_orig=%f   e+delta=%f  e-delta=%f\n", dx, dy, dz, e_orig, e_plus_delta, e_minus_delta); // DEBUG
    //obErrorLog.ThrowError(__FUNCTION__, errbuf, obError); // DEBUG
    return (grad);
  }

  void OBForceField::UpdateCoordinates(OBMol &mol)
  {
    OBAtom *atom;

    FOR_ATOMS_OF_MOL (a, mol) {
      atom = _mol.GetAtom(a->GetIdx());
      a->SetVector(atom->GetVector());
    }
  }

  double OBForceField::PointPlaneAngle(const vector3 &a, const vector3 &b, const vector3 &c, const vector3 &d)
  {
    vector3 ab, bc, bd, normal;
    double angle;

    ab = a - b;
    bc = b - c;
    bd = b - d;
 
    normal = cross(ab, bc);
    angle = 90.0f - vectorAngle(normal, bd);

    return angle;
  }

} // end namespace OpenBabel


//! \file forcefield.cpp
//! \brief Handle OBForceField class

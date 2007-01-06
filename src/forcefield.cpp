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
      for (unsigned int idx=0; idx < parameter.size(); idx++) {
        if (((_a == parameter[idx]._a) && (_b == parameter[idx]._b)) || ((_a == parameter[idx]._b) && (_b == parameter[idx]._a))) {
	  par = &parameter[idx];
	  return par;
	}
      }
      return NULL;
    }
    if (d == NULL) {
      std::string _a(a);
      std::string _b(b);
      std::string _c(c);
      for (unsigned int idx=0; idx < parameter.size(); idx++) {
        if (((_a == parameter[idx]._a) && (_b == parameter[idx]._b) && (_c == parameter[idx]._c)) || 
	    ((_a == parameter[idx]._c) && (_b == parameter[idx]._b) && (_c == parameter[idx]._a))) {
	  par = &parameter[idx];
	  return par;
        }
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

  int OBForceField::get_nbr (OBAtom* atom, int level) {
    OBAtom *nbr,*nbr2;
    vector<OBEdgeBase*>::iterator i;
  
    if (level == 2)
      if (!get_nbr(atom, 1)) return 0;
    if (level == 3)
      if (!get_nbr(atom, 2)) return 0;
  
    // Find first neighboor
    FOR_NBORS_OF_ATOM(tmp, atom) {
      if (atom->GetIdx() > tmp->GetIdx()) {
        if (level == 1)
          return tmp->GetIdx();
        else {
          nbr = _mol.GetAtom(tmp->GetIdx());
          break;
	}
      }
    }
    if (level == 1) return 0;
  
    // Find second neighboor
    FOR_NBORS_OF_ATOM(tmp, nbr) {
      if (atom->GetIdx() > tmp->GetIdx()) {
        if (level == 2)
          return tmp->GetIdx();
        else {
          nbr2 = _mol.GetAtom(tmp->GetIdx());
          break;
	}
      }
    }
    if (level == 2) return 0;
  
    // Find thirth neighboor
    FOR_NBORS_OF_ATOM(tmp, nbr2) {
      if ((atom->GetIdx() > tmp->GetIdx()) && (nbr->GetIdx() != tmp->GetIdx()))
        return tmp->GetIdx();
    }
    FOR_NBORS_OF_ATOM(tmp, nbr) {
      if ((atom->GetIdx() > tmp->GetIdx()) && (atom->GetIdx() != tmp->GetIdx()) && (nbr2->GetIdx() != tmp->GetIdx()))
        return tmp->GetIdx();
    }
    
    return 0;
  }

  void OBForceField::GenerateCoordinates() 
  {
    //_mol.AddHydrogens(false, true);

    OBAtom *atom, *nbr, *nbr2, *nbr3;
    vector<OBNodeBase*>::iterator i;
    vector<OBEdgeBase*>::iterator j;

    vector<OBInternalCoord*> internals;
    OBInternalCoord *coord;

    coord = new OBInternalCoord();
    internals.push_back(coord);
      
    int torang;
    for (atom = _mol.BeginAtom(i);atom;atom = _mol.NextAtom(i)) {
      coord = new OBInternalCoord();
      nbr = _mol.GetAtom(get_nbr(atom, 1));
      nbr2 = _mol.GetAtom(get_nbr(atom, 2));
      nbr3 = _mol.GetAtom(get_nbr(atom, 3));
        
      if (nbr) {
        coord->_a = _mol.GetAtom(get_nbr(atom, 1));
        OBBond *bond;
        if ( (bond = _mol.GetBond(atom, nbr)) ) {
          coord->_dst = bond->GetEquibLength();
        }
      }

      if (nbr2) {
        coord->_b = _mol.GetAtom(get_nbr(atom, 2));
        if (nbr->GetHyb() == 3)
          coord->_ang = 109;
        if (nbr->GetHyb() == 2)
          coord->_ang = 120;
        if (nbr->GetHyb() == 1)
          coord->_ang = 180;
      }
  
      if (nbr3) {
        // double bestangle, angle, bestscore, score;
        // int nbr_count;
        coord->_c = _mol.GetAtom(get_nbr(atom, 3));
        coord->_tor = torang;
        torang +=60;
      }
            
      internals.push_back(coord);
    }
    
    InternalToCartesian(internals, _mol);
   
    // minimize the created structure
    SteepestDescent(150);
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

    e_n1 = Energy(); // we call Energy instead of GetEnergy 
                     // because coordinates change every step

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

    UnsetEnergyCalculated();
    std::cout << std::endl;
  }

  void OBForceField::ConjugateGradients(int steps)
  {
    double e_n1, e_n2;
    double h, g2g2, g1g1, g2g1;
    bool firststep;
    vector3 grad1, grad2, dir1, dir2;
    std::vector<vector3> old_xyz;

    h = 0.1f;
    firststep = true;
    old_xyz.resize(_mol.NumAtoms()+1);

    e_n1 = Energy();

    for (int i=0; i<steps; i++) {
      if (firststep) {
        FOR_ATOMS_OF_MOL (a, _mol) {
  	  grad1 = NumericalDerivative(a->GetIdx());
	  grad1 = grad1.normalize();
	  grad1 = grad1 * h;
	  old_xyz[a->GetIdx()] = a->GetVector();
          a->SetVector(a->x() + grad1.x(), a->y() + grad1.y(), a->z() + grad1.z());
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
	  firststep = false;
	  dir1 = grad1;
        }
      } else {
        FOR_ATOMS_OF_MOL (a, _mol) {
 	  grad2 = NumericalDerivative(a->GetIdx());
	  grad2 = grad2.normalize();
	  grad2 = grad2 * h;
	  g2g2 = dot(grad2, grad2);
	  g1g1 = dot(grad1, grad1);
	  g2g1 = g2g2 / g1g1;
	  dir2 = grad2 + g2g1 * dir1;
	  old_xyz[a->GetIdx()] = a->GetVector();
          a->SetVector(a->x() + dir2.x(), a->y() + dir2.y(), a->z() + dir2.z());
	  //std::cout << "  dir2=" << dir2 << std::endl;
          //sprintf(errbuf, "%svf=(%f, %f, %f)\n", errbuf, vf.x(), vf.y(), vf.z()); // DEBUG
	  grad1 = grad2;
	  dir1 = dir2;
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
    }
    UnsetEnergyCalculated();
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

    UnsetEnergyCalculated();
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

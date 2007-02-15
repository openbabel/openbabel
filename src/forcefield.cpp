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
#include <openbabel/math/matrix3x3.h>

using namespace std;

namespace OpenBabel
{
  int OBForceField::GetParameterIdx(int a, int b, int c, int d, vector<OBFFParameter> &parameter)
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
  
  OBFFParameter* OBForceField::GetParameter(int a, int b, int c, int d, vector<OBFFParameter> &parameter)
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
  
  OBFFParameter* OBForceField::GetParameter(const char* a, const char* b, const char* c, const char* d, vector<OBFFParameter> &parameter)
  {
    OBFFParameter *par;
    if (a == NULL)
      return NULL;

    if (b == NULL) {
      string _a(a);
      for (unsigned int idx=0; idx < parameter.size(); idx++) 
        if (_a == parameter[idx]._a) {
	  par = &parameter[idx];
	  return par;
	}
      return NULL;
    }
    if (c == NULL) {
      string _a(a);
      string _b(b);
      for (unsigned int idx=0; idx < parameter.size(); idx++) {
        if (((_a == parameter[idx]._a) && (_b == parameter[idx]._b)) || ((_a == parameter[idx]._b) && (_b == parameter[idx]._a))) {
	  par = &parameter[idx];
	  return par;
	}
      }
      return NULL;
    }
    if (d == NULL) {
      string _a(a);
      string _b(b);
      string _c(c);
      for (unsigned int idx=0; idx < parameter.size(); idx++) {
        if (((_a == parameter[idx]._a) && (_b == parameter[idx]._b) && (_c == parameter[idx]._c)) || 
	    ((_a == parameter[idx]._c) && (_b == parameter[idx]._b) && (_c == parameter[idx]._a))) {
	  par = &parameter[idx];
	  return par;
        }
      }
      return NULL;
    }
    string _a(a);
    string _b(b);
    string _c(c);
    string _d(d);
    for (unsigned int idx=0; idx < parameter.size(); idx++)
      if (((_a == parameter[idx]._a) && (_b == parameter[idx]._b) && (_c == parameter[idx]._c) && (_d == parameter[idx]._d)) || 
          ((_a == parameter[idx]._d) && (_b == parameter[idx]._c) && (_c == parameter[idx]._b) && (_d == parameter[idx]._a))) {
	par = &parameter[idx];
	return par;
      }

    return NULL;
  }
  
  bool OBForceField::SetLogFile(ostream* pos)
  {
    if(pos)
       logos = pos;
    else
       logos = &cout;
    
    return true;
  }
  
  bool OBForceField::SetLogLevel(int level)
  {
    loglvl = level; 

    return true;
  }
 
  bool is14(OBAtom *a, OBAtom *b)
  {
    FOR_NBORS_OF_ATOM (nbr, a)
      FOR_NBORS_OF_ATOM (nbr2, &*nbr)
        FOR_NBORS_OF_ATOM (nbr3, &*nbr2)
	  if (b == &*nbr3)
	    return true;
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
    //SteepestDescent(300);
  }

  void OBForceField::DistanceGeometry() 
  {
    int N = _mol.NumAtoms();
    int i = 0;
    int j = 0;
    double matrix[N][N], G[N][N];
    bool is15;
    
    IF_OBFF_LOGLVL_LOW {
      *logos << endl << "D I S T A N C E   G E O M E T R Y" << endl << endl;
    }
 
    // Calculate initial distance matrix
    //
    // - diagonal elements are 0.0
    // - when atoms i and j form a bond, the bond length is used for
    //   upper and lower limit
    // - when atoms i and j are in a 1-3 relationship, the distance
    //   is calculated using the cosine rule: (upper limit = lower limit)
    //      
    //          b_         ab: bond length
    //         /  \_       bc: bond length
    //        /A    \_     ac = sqrt(ab^2 + bc^2 - 2*ab*bc*cos(A))
    //       a--------c
    //
    // - when atoms i anf j are in a 1-4 relationship, the lower distance
    //   limit is calculated using a torsional angle of 0.0. The upper limit 
    //   is calculated using a torsion angle of 180.0.
    //
    //      a       d      ab, bc, cd: bond lengths
    //       \ B C /       
    //        b---c        ad = bc + ab*cos(180-B) + cd*cos(180-C)
    //
    //      a
    //       \ B           delta_x = bc + ab*cos(180-B) + cd*cos(180-C)
    //        b---c        delta_y = ab*sin(180-B) + cd*sin(180-C)
    //           C \
    //              d      ad = sqrt(delta_x^2 + delta_y^2)
    //
    FOR_ATOMS_OF_MOL (a, _mol) {
      i = a->GetIdx() - 1;
      FOR_ATOMS_OF_MOL (b, _mol) {
        j = b->GetIdx() - 1;

	if (&*a == &*b) {
	  matrix[i][j] = 0.0f; // diagonal
	  continue;
	}
        // Find relationship
	is15 = true;
	FOR_NBORS_OF_ATOM (nbr1, _mol.GetAtom(a->GetIdx())) { // 1-2
	  if (&*nbr1 == &*b) {
	    matrix[i][j] = 1.3f;
	    break;
	  }
	  FOR_NBORS_OF_ATOM (nbr2, _mol.GetAtom(nbr1->GetIdx())) { // 1-3
	    if (&*nbr2 == &*b) {
	      matrix[i][j] = sqrt(1.3f*1.3f + 1.3f*1.3f - 2.0f * cos(DEG_TO_RAD*120.0f) * 1.3f*1.3f );
	      is15 = false;
	      break;
	    }
            FOR_NBORS_OF_ATOM (nbr3, &*nbr2) { // 1-4
	      if (&*nbr3 == &*b) {
	        is15 = false;
	        if (i > j) // minimum distance (torsion angle = 0)
	          matrix[i][j] = 1.3f + 1.3f*cos(DEG_TO_RAD*(180.0f-120.0f)) + 1.3f*cos(DEG_TO_RAD*(180.0f-120.0f));
		else {// maximum distance (torsion angle = 180)
		  double delta_x, delta_y;
		  delta_x = 1.3f + 1.3f*cos(DEG_TO_RAD*(180.0f-120.0f)) + 1.3f*cos(DEG_TO_RAD*(180.0f-120.0f));
		  delta_y = 1.3f*sin(DEG_TO_RAD*(180.0f-120.0f)) + 1.3f*sin(DEG_TO_RAD*(180.0f-120.0f));
		  matrix[i][j] = sqrt(delta_x*delta_x + delta_y*delta_y);
		}
	        break;
	      }
              if (i > j && is15) {// minimum distance (sum vdw radii)
	        matrix[i][j] = 1.4f + 1.4f;
              } else if (is15) // maximum distance (torsion angle = 180)
	        matrix[i][j] = 99.0f;
	    }
	  }
	}
      }
    }
   
    // output initial distance matrix
    IF_OBFF_LOGLVL_LOW {
      char logbuf[100];

      *logos << "INITIAL DISTANCE MATRIX" << endl << endl;
      for (i=0; i<N; i++) {
        *logos << "[";
        for (j=0; j<N; j++) {
          sprintf(logbuf, " %8.4f ", matrix[i][j]);
	  *logos << logbuf;
        }
        *logos << "]" << endl;
      }
      *logos << endl;
    }

    // Triangle smoothing
    //FOR_ANGLES_OF_MOL(angle, _mol) {
    int a, b, c;
    bool self_consistent = false;
    while (!self_consistent) {
      self_consistent = true;

      FOR_ATOMS_OF_MOL (_a, _mol) {
        a = _a->GetIdx() - 1;
        FOR_ATOMS_OF_MOL (_b, _mol) {
          if (&*_b == &*_a)
            continue;
          b = _b->GetIdx() - 1;
          FOR_ATOMS_OF_MOL (_c, _mol) {
            if ((&*_c == &*_b) || (&*_c == &*_a))
              continue;
            c = _c->GetIdx() - 1;
      
            double u_ab, u_bc, u_ac; // upper limits
            double l_ab, l_bc, l_ac; // lower limits
      
            // get the upper and lower limits for ab, bc and ac
            if (b > a) {
              u_ab = matrix[a][b];
	      l_ab = matrix[b][a];
            } else {
              u_ab = matrix[b][a];
              l_ab = matrix[a][b];
            }
            if (c > b) {
              u_bc = matrix[b][c];
              l_bc = matrix[c][b];
            } else {
              u_bc = matrix[c][b];
	      l_bc = matrix[b][c];
            }
            if (c > a) {
              u_ac = matrix[a][c];
              l_ac = matrix[c][a];
           } else {
              u_ac = matrix[c][a];
              l_ac = matrix[a][c];
            }

            if (u_ac > (u_ab + u_bc)) { // u_ac <= u_ab + u_bc
              u_ac = u_ab + u_bc;
	      self_consistent = false;
	    }

            if (l_ac < (l_ab - u_bc)) {// l_ac >= l_ab - u_bc
              l_ac = l_ab - u_bc;
      	      self_consistent = false;
	    }

            // store smoothed l_ac and u_ac
            if (c > a) {
              matrix[a][c] = u_ac;
              matrix[c][a] = l_ac;
            } else {
              matrix[c][a] = u_ac;
	      matrix[a][c] = l_ac;
            }
          }
        }
      }
    }

    // output result of triangle smoothing
    IF_OBFF_LOGLVL_LOW {
      char logbuf[100];

      *logos << "TRIANGLE SMOOTHING" << endl << endl;
      for (i=0; i<N; i++) {
        *logos << "[";
        for (j=0; j<N; j++) {
          sprintf(logbuf, " %8.4f ", matrix[i][j]);
	  *logos << logbuf;
        }
        *logos << "]" << endl;
      }
      *logos << endl;
    }
    
    // Generate random distance matrix between lower and upper limits
    FOR_ATOMS_OF_MOL (a, _mol) {
      i = a->GetIdx() - 1;
      FOR_ATOMS_OF_MOL (b, _mol) {
        j = b->GetIdx() - 1;

	if (&*a == &*b) {
	  matrix[i][j] = 0.0f; // diagonal
	  continue;
	}
       
        srand(time(NULL));
	double rand_ab, u_ab, l_ab;
        if (j > i) {
          u_ab = matrix[i][j];
          l_ab = matrix[j][i];
          rand_ab = l_ab + (u_ab - l_ab) * rand()/RAND_MAX;
	  matrix[i][j] = rand_ab;
	  matrix[j][i] = rand_ab;
        } else {
          u_ab = matrix[j][i];
          l_ab = matrix[i][j];
          rand_ab = l_ab + (u_ab - l_ab) * rand()/RAND_MAX;
	  matrix[i][j] = rand_ab;
	  matrix[j][i] = rand_ab;
        }
      }
    }
 
    // output result of triangle smoothing
    IF_OBFF_LOGLVL_LOW {
      char logbuf[100];

      *logos << "RANDOM DISTANCE MATRIX BETWEEN LIMITS" << endl << endl;
      for (i=0; i<N; i++) {
        *logos << "[";
        for (j=0; j<N; j++) {
          sprintf(logbuf, " %8.4f ", matrix[i][j]);
	  *logos << logbuf;
        }
        *logos << "]" << endl;
      }
      *logos << endl;
    }

    // Generate metric matrix
    // (origin = first atom )
    for (i=0; i<N; i++) {
      for (j=0; j<N; j++) {
        G[i][j] = 0.5f * (matrix[0][i]*matrix[0][i] + matrix[0][j]*matrix[0][j] - matrix[i][j]*matrix[i][j]);
      }
    }
    
    // output metric matrix
    IF_OBFF_LOGLVL_LOW {
      char logbuf[100];

      *logos << "METRIC MATRIX" << endl << endl;
      for (i=0; i<N; i++) {
        *logos << "[";
        for (j=0; j<N; j++) {
          sprintf(logbuf, " %8.4f ", G[i][j]);
	  *logos << logbuf;
        }
        *logos << "]" << endl;
      }
      *logos << endl;
    }

    // Calculate eigenvalues and eigenvectors
    double eigenvalues[N];
    double eigenvectors[N][N];
    matrix3x3::jacobi(N, (double *) &G, (double *) &eigenvalues, (double *) &eigenvectors);
    
    // output eigenvalues and eigenvectors
    IF_OBFF_LOGLVL_LOW {
      char logbuf[100];

      *logos << "EIGENVALUES OF METRIC MATRIX" << endl << endl;
      for (i=0; i<N; i++) {
        *logos << eigenvalues[i] << " ";
      }
      *logos << endl << endl;

      *logos << "EIGENVECTORS OF METRIC MATRIX" << endl << endl;
      for (i=0; i<N; i++) {
        *logos << "[";
        for (j=0; j<N; j++) {
          sprintf(logbuf, " %8.4f ", eigenvectors[i][j]);
	  *logos << logbuf;
        }
        *logos << "]" << endl;
      }
      *logos << endl;
    }

    // Assign coordinates
    double xa, ya, za;
    FOR_ATOMS_OF_MOL (a, _mol) {
      i = a->GetIdx() - 1;
      
      if (eigenvectors[i][N-1] > 0)
        xa = sqrt(eigenvalues[N-1] * eigenvectors[i][N-1]);
      else
        xa = -sqrt(eigenvalues[N-1] * -eigenvectors[i][N-1]);

      if (eigenvectors[i][N-2] > 0)
        ya = sqrt(eigenvalues[N-2] * eigenvectors[i][N-2]);
      else
        ya = -sqrt(eigenvalues[N-2] * -eigenvectors[i][N-2]);

      if (eigenvectors[i][N-3] > 0)
        za = sqrt(eigenvalues[N-3] * eigenvectors[i][N-3]);
      else
        za = -sqrt(eigenvalues[N-3] * -eigenvectors[i][N-3]);

      a->SetVector(xa, ya, za);
    }

  }
  
  // LineSearch 
  //
  // atom: coordinates of atom at iteration k (x_k)
  // direction: search direction ( d = -grad(x_0) )
  //
  // ALGORITHM:
  // 
  // step = 1
  // for (i = 1 to 100) {                max steps = 100
  //   e_k = energy(x_k)                 energy of current iteration
  //   x_k = x_k + step * d              update coordinates
  //   e_k+1 = energy(x_k+1)             energy of next iteration
  //   
  //   if (e_k+1 < e_k)
  //     step = step * 1.2               increase step size
  //   if (e_k+1 > e_k) {
  //     x_k = x_k - step * d            reset coordinates to previous iteration
  //     step = step * 0.5               reduce step size
  //   }
  //   if (e_k+1 == e_k)
  //     end                             convergence criteria reached, stop
  // }
  vector3 OBForceField::LineSearch(OBAtom *atom, vector3 &direction)
  {
    double e_n1, e_n2, step;
    vector3 old_xyz, orig_xyz, xyz_k, dir(0.0f, 0.0f, 0.0f);

    step = 0.2f;
    direction.normalize();
    orig_xyz = atom->GetVector();
    
    e_n1 = Energy(); // calculate e_k
    
    for (int i=0; i<100; i++) {
      old_xyz = atom->GetVector();
      
      xyz_k = atom->GetVector() + direction*step;
      atom->SetVector(xyz_k);  // update coordinates
    
      e_n2 = Energy(); // calculate e_k+1
      
      if (e_n2 == e_n1) // convergence criteria
        break;

      if (e_n2 > e_n1) { // decrease stepsize
	step *= 0.5f;
        atom->SetVector(old_xyz);
      }
      if (e_n2 < e_n1) {  // increase stepsize
        e_n1 = e_n2;
	step *= 1.2f;
	if (step > 1.0f)
	  step = 1.0f;
      }
      
    }

    dir = atom->GetVector() - orig_xyz;
    atom->SetVector(orig_xyz);     
    return dir;    
  }
  
  vector3 OBForceField::ValidateLineSearch(OBAtom *atom, vector3 &direction)
  {
    double e_n1, e_n2, step;
    vector3 old_xyz, orig_xyz, xyz_k, dir(0.0f, 0.0f, 0.0f);

    step = 0.2f;
    direction.normalize();
    orig_xyz = atom->GetVector();
    
    e_n1 = atom->x() * atom->x() + 2 * (atom->y() * atom->y()); // e_k
    
    for (int i=0; i<100; i++) {
      old_xyz = atom->GetVector();
      
      xyz_k = atom->GetVector() + direction*step;
      atom->SetVector(xyz_k);  // update coordinates
    
      e_n2 = atom->x() * atom->x() + 2 * (atom->y() * atom->y()); // e_k+1
      
      if (e_n2 == e_n1) // convergence criteria
        break;

      if (e_n2 > e_n1) { // decrease stepsize
	step *= 0.5f;
        atom->SetVector(old_xyz);
      }
      if (e_n2 < e_n1) {  // increase stepsize
        e_n1 = e_n2;
	step *= 1.2f;
	if (step > 1.0f)
	  step = 1.0f;
      }
      
    }

    dir = atom->GetVector() - orig_xyz;
    atom->SetVector(orig_xyz);     
    return dir;    
  }
 
  // used to validate the SteepestDescent implementation using a simple function.
  //
  // f(x,y) = x^2 + 2y^2
  // minimum: (0, 0)
  // df/dx = 2x
  // df/dy = 4y
  //
  void OBForceField::ValidateSteepestDescent(int steps) 
  {
    OBAtom *atom = new OBAtom;
    vector3 grad;
    double e_n1, e_n2;
    char logbuf[100];
    
    atom->SetVector(9.0f, 9.0f, 0.0f);
    e_n1 = atom->x() * atom->x() + 2 * (atom->y() * atom->y());
    
    cout << endl << "V A L I D A T E   S T E E P E S T   D E S C E N T" << endl << endl;
    cout << "STEPS = " << steps << endl << endl;
    cout << "STEP n     E(n)       E(n-1)    " << endl;
    cout << "--------------------------------" << endl;
    
    for (int i=1; i<=steps; i++) {
      grad.Set(-2*atom->x(), -4*atom->y(), 0.0f);
      grad = ValidateLineSearch(&*atom, grad);
      atom->SetVector(atom->x() + grad.x(), atom->y() + grad.y(), 0.0f);
      e_n2 = atom->x() * atom->x() + 2 * (atom->y() * atom->y());
      
      sprintf(logbuf, " %4d    %8.3f    %8.3f", i, e_n2, e_n1);
      cout << logbuf << endl;

      if (fabs(e_n1 - e_n2) < 0.0000001f) {
        cout << "    STEEPEST DESCENT HAS CONVERGED (DELTA E < 0.0000001)" << endl;
        break;
      }

      e_n1 = e_n2;
    }

    UnsetEnergyCalculated();
    IF_OBFF_LOGLVL_LOW
      *logos << endl;
  }

  void OBForceField::ValidateConjugateGradients(int steps)
  {
    OBAtom *atom = new OBAtom;
    vector3 grad1, grad2, dir1, dir2;
    double e_n1, e_n2;
    double g2g2, g1g1, g2g1;
    bool firststep;
    char logbuf[100];

    firststep = true;
    atom->SetVector(9.0f, 9.0f, 0.0f);
    e_n1 = atom->x() * atom->x() + 2 * (atom->y() * atom->y());
 
    cout << endl << "V A L I D A T E   C O N J U G A T E   G R A D I E N T" << endl << endl;
    cout << "STEPS = " << steps << endl << endl;
    cout << "STEP n     E(n)       E(n-1)    " << endl;
    cout << "--------------------------------" << endl;
 
    for (int i=1; i<=steps; i++) {
      if (firststep) {
        grad1.Set(-2*atom->x(), -4*atom->y(), 0.0f);
	dir1 = grad1;
	dir1 = ValidateLineSearch(&*atom, dir1);
        atom->SetVector(atom->x() + dir1.x(), atom->y() + dir1.y(), atom->z() + dir1.z());
        e_n2 = atom->x() * atom->x() + 2 * (atom->y() * atom->y());
      
        sprintf(logbuf, " %4d    %8.3f    %8.3f", i, e_n2, e_n1);
        cout << logbuf << endl;
        
	e_n1 = e_n2;
	dir1 = grad1;
	firststep = false;
      } else {
        grad2.Set(-2*atom->x(), -4*atom->y(), 0.0f);
	g2g2 = dot(grad2, grad2);
	g1g1 = dot(grad1, grad1);
	g2g1 = g2g2 / g1g1;
	dir2 = grad2 + g2g1 * dir1;
	dir2 = ValidateLineSearch(&*atom, dir2);
        atom->SetVector(atom->x() + dir2.x(), atom->y() + dir2.y(), atom->z() + dir2.z());
	grad1 = grad2;
	dir1 = dir2;
        e_n2 = atom->x() * atom->x() + 2 * (atom->y() * atom->y());
	  
        sprintf(logbuf, " %4d    %8.3f    %8.3f", i, e_n2, e_n1);
        cout << logbuf << endl;
        
	if (fabs(e_n1 - e_n2) < 0.0000001f) {
          cout << "    CONJUGATE GRADIENTS HAS CONVERGED (DELTA E < 0.0000001)" << endl;
          break;
        }

	//grad1 = grad2;
	//dir1 = dir2;
	e_n1 = e_n2;
      }
    }
    UnsetEnergyCalculated();
  }

  void OBForceField::SteepestDescent(int steps) 
  {
    double e_n1, e_n2;
    char logbuf[100];
    vector3 grad;

    e_n1 = Energy(); // we call Energy instead of GetEnergy 
                     // because coordinates change every step
    
    IF_OBFF_LOGLVL_LOW {
      *logos << endl << "S T E E P E S T   D E S C E N T" << endl << endl;
      *logos << "STEPS = " << steps << endl << endl;
      *logos << "STEP n     E(n)       E(n-1)    " << endl;
      *logos << "--------------------------------" << endl;
    }
    
    for (int i=1; i<=steps; i++) {
      
      FOR_ATOMS_OF_MOL (a, _mol) {
	//grad = NumericalDerivative(a->GetIdx());
	grad = GetGradient(&*a);
	grad = LineSearch(&*a, grad);
        a->SetVector(a->x() + grad.x(), a->y() + grad.y(), a->z() + grad.z());
      }
      e_n2 = Energy();
      
      IF_OBFF_LOGLVL_LOW {
        sprintf(logbuf, " %4d    %8.3f    %8.3f", i, e_n2, e_n1);
        *logos << logbuf << endl;
      }

      if (fabs(e_n1 - e_n2) < 0.0000001f) {
        *logos << "    STEEPEST DESCENT HAS CONVERGED (DELTA E < 0.0000001)" << endl;
        break;
      }

      e_n1 = e_n2;
    }

    UnsetEnergyCalculated();
    IF_OBFF_LOGLVL_LOW
      *logos << endl;
  }

  void OBForceField::ConjugateGradients(int steps)
  {
    double e_n1, e_n2;
    double g2g2, g1g1, g2g1;
    bool firststep;
    vector<vector3> grad1, dir1;
    vector3 grad2, dir2;
    vector<vector3> old_xyz;
    char logbuf[100];

    firststep = true;

    e_n1 = Energy();
    
    IF_OBFF_LOGLVL_LOW {
      *logos << endl << "C O N J U G A T E   G R A D I E N T" << endl << endl;
      *logos << "STEPS = " << steps << endl << endl;
      *logos << "STEP n     E(n)       E(n-1)    " << endl;
      *logos << "--------------------------------" << endl;
    }
    
    grad1.resize(_mol.NumAtoms() + 1);
    dir1.resize(_mol.NumAtoms() + 1);

    for (int i=0; i<steps; i++) {
      if (firststep) {
        FOR_ATOMS_OF_MOL (a, _mol) {
  	  //grad1 = NumericalDerivative(a->GetIdx());
	  grad2 = GetGradient(&*a);
	  dir2 = grad2;
	  dir2 = LineSearch(&*a, dir2);
          a->SetVector(a->x() + dir2.x(), a->y() + dir2.y(), a->z() + dir2.z());
	  grad1[a->GetIdx()] = grad2;
	  dir1[a->GetIdx()] = grad2;
        }
        e_n2 = Energy();
      
        IF_OBFF_LOGLVL_LOW {
          sprintf(logbuf, " %4d    %8.3f    %8.3f", i, e_n2, e_n1);
          *logos << logbuf << endl;
        }
	
	firststep = false;
	e_n1 = e_n2;
      } else {
        FOR_ATOMS_OF_MOL (a, _mol) {
 	  //grad2 = NumericalDerivative(a->GetIdx());
	  grad2 = GetGradient(&*a);
	  g2g2 = dot(grad2, grad2);
	  g1g1 = dot(grad1[a->GetIdx()], grad1[a->GetIdx()]);
	  g2g1 = g2g2 / g1g1;
	  dir2 = grad2 + g2g1 * dir1[a->GetIdx()];
	  dir2 = LineSearch(&*a, dir2);
          a->SetVector(a->x() + dir2.x(), a->y() + dir2.y(), a->z() + dir2.z());
	  
	  grad1[a->GetIdx()] = grad2;
	  dir1[a->GetIdx()] = dir2;
	  e_n1 = e_n2;
        }
        e_n2 = Energy();
	
	IF_OBFF_LOGLVL_LOW {
          sprintf(logbuf, " %4d    %8.3f    %8.3f", i, e_n2, e_n1);
          *logos << logbuf << endl;
        }
 
	if (fabs(e_n1 - e_n2) < 0.0000001f) {
          cout << "    CONJUGATE GRADIENTS HAS CONVERGED (DELTA E < 0.0000001)" << endl;
          break;
        }

	e_n1 = e_n2;
      }
    }
    UnsetEnergyCalculated();
  }

  vector3 OBForceField::NumericalDerivative(int a)
  {
    OBAtom *atom;
    vector3 va, grad;
    double e_orig, e_plus_delta, e_minus_delta, delta, dx, dy, dz;

    delta = 0.000001f;

    atom = _mol.GetAtom(a);
    va = atom->GetVector();

    e_orig = Energy();
    
    atom->SetVector(va.x() + delta, va.y(), va.z());
    e_plus_delta = Energy();
    atom->SetVector(va.x() - delta, va.y(), va.z());
    e_minus_delta = Energy();

    //dx = (e_plus_delta - e_minus_delta);
    dx = (e_plus_delta - e_orig);
    //cout << "x_" << a << "    " << e_plus_delta << " - " << e_minus_delta << "  (orig=" << e_orig << ")" <<endl;

    atom->SetVector(va.x(), va.y() + delta, va.z());
    e_plus_delta = Energy();
    atom->SetVector(va.x(), va.y() - delta, va.z());
    e_minus_delta = Energy();

    //dy = (e_plus_delta - e_minus_delta);
    dy = (e_plus_delta - e_orig);
    //cout << "y_" << a << "    " << e_plus_delta << " - " << e_minus_delta << "  (orig=" << e_orig << ")" << endl;
    
    atom->SetVector(va.x(), va.y(), va.z() + delta);
    e_plus_delta = Energy();
    atom->SetVector(va.x(), va.y(), va.z() - delta);
    e_minus_delta = Energy();

    //dz = (e_plus_delta - e_minus_delta);
    dz = (e_plus_delta - e_orig);
    //cout << "z_" << a << "    " << e_plus_delta << " - " << e_minus_delta << "  (orig=" << e_orig << ")" << endl;

    atom->SetVector(va.x(), va.y(), va.z());

    dx *= 1.0f / delta;
    dy *= 1.0f / delta;
    dz *= 1.0f / delta;
    grad.Set(-dx, -dy, -dz);

    //cout << "NumericalDerivative(" << a << ")" << endl;
    //cout << "grad: " << grad << endl;
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

} // end namespace OpenBabel


//! \file forcefield.cpp
//! \brief Handle OBForceField class

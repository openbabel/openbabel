/**********************************************************************
forcefield.cpp - Handle OBForceField class.

Copyright (C) 2006-2007 by Tim Vandermeersch <tim.vandermeersch@gmail.com>
 
This file is part of the Open Babel project.
For more information, see <http://openbabel.sourceforge.net/>
 
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

      development status:
      - src/forcefield.cpp
      - LineSearch(): works, but better algoritmh would be better for performance
      - SteepestDescent(): finished
      - ConjugateGradients(): finished
      - GenerateCoordinates(): 
      - SystematicRotorSearch(): not all combinations are tested but works
      - DistanceGeometry(): needs matrix operations

      - src/forcefields/forcefieldghemical.cpp
      - Atom typing: finished
      - Charges: finished
      - Energy terms: finished
      - Analytical gradients: finished
      - Validation: in progress...
        1,4-scaling still needs some work

      - src/forcefields/forcefieldmmff94.cpp
      - Atom typing: needs work
      - Charges: 0%
      - Energy terms:
	    - Bond: finished
	    - Angle: finished
	    - StrBnd: finished
	    - Torsion: finished
	    - OOP: no gradient
	    - VDW: no gradient
	    - Electrostatic: needs charges, charges need correct atom types
      - Validation: in progress... 


***********************************************************************/
#include <openbabel/babelconfig.h>

#include <openbabel/forcefield.h>

#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obiter.h>
#include <openbabel/math/matrix3x3.h>
#include <openbabel/rotamer.h>
#include <openbabel/rotor.h>

using namespace std;

namespace OpenBabel
{
  /** \class OBForceField forcefield.h <openbabel/forcefield.h>
      \brief Base class for molecular mechanics force fields

      The OBForceField class is the base class for molecular mechanics in Open Babel.
      Classes derived from the OBForceField implement specific force fields (Ghemical, MMFF94, ...).
      Other classes such as OBFFParameter, OBFFCalculation and its derived classes are 
      only for internal use. As a user interested in using the available force fields 
      in Open Babel, you don't need these classes. The rest of this short introduction
      is aimed at these users. For information on how to implement additional force
      fields, see the wiki pages or post your questions to the openbabel-devel mailing
      list.

      Before we can start using a force field, we must first select it and set it up.
      This is illustrated in the first example below. The Setup procedure assigns atom
      types, charges and parameters. There are several reasons why this may fail, a log
      message will be written to the logfile before Setup() returns false.

      The force field classes use their own logging functions. You can set the logfile 
      using SetLogFile() and set the log level using SetLogLevel(). If needed you can
      also write to the logfile using OBFFLog(). There are four log levels: 
      BFF_LOGLVL_NONE, OBFF_LOGLVL_LOW, OBFF_LOGLVL_MEDIUM, OBFF_LOGLVL_HIGH.
      See the API documentation to know what kind of output each function writes to
      the logfile for the different log levels.

      Below are two examples which explain the basics. 
      
      This piece of code will output a list of available forcefields to cout:
      \code
      OBPlugin::List("forcefields");
      \endcode

      Calculate the energy for the structure in mol using the Ghemical forcefield.
      \code
      #include <openbabel/forcefield.h>
      #include <openbabel/mol.h>

      // See OBConversion class to fill the mol object.
      OBMol mol;
      // Select the forcefield, this returns a pointer that we 
      // will later use to access the forcefield functions.
      OBForceField* pFF = OBForceField::FindForceField("Ghemical");

      // Make sure we have a valid pointer
      if (!pFF)
      // exit...
     
      // Set the logfile (can also be &cout or &cerr)
      pFF->SetLogFile(&cerr);
      // Set the log level. See indivual functions to know
      // what kind of output each function produces for the
      // different log levels.
      pFF->SetLogLevel(OBFF_LOGLVL_HIGH);

      // We need to setup the forcefield before we can use it. Setup()
      // returns false if it failes to find the atom types, parameters, ...
      if (!pFF->Setup(mol)) {
      cerr << "ERROR: could not setup force field." << endl;
      }
      
      // Calculate the energy. The output will be written to the
      // logfile specified by SetLogFile()
      pFF->Energy();
      \endcode
 
      Minimize the structure in mol using conjugate gradients.
      \code
      #include <openbabel/forcefield.h>
      #include <openbabel/mol.h>

      OBMol mol;
      OBForceField* pFF = OBForceField::FindForceField("Ghemical");

      // Make sure we have a valid pointer
      if (!pFF)
      // exit...
      
      pFF->SetLogFile(&cerr);
      pFF->SetLogLevel(OBFF_LOGLVL_LOW);
      if (!pFF->Setup(mol)) {
      cerr << "ERROR: could not setup force field." << endl;
      }
      
      // Perform the actual minimization, maximum 1000 steps 
      pFF->ConjugateGradients(1000);
      \endcode
  **/

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
  
    if (!nbr) return 0; // ensure nbr is not used uninitialized

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

    if (!nbr2) return 0; // ensure nbr2 is not used uninitialized

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
    OBAtom *atom, *nbr, *nbr2, *nbr3;
    vector<OBNodeBase*>::iterator i;
    vector<OBEdgeBase*>::iterator j;

    vector<OBInternalCoord*> internals;
    OBInternalCoord *coord;

    coord = new OBInternalCoord();
    internals.push_back(coord);
      
    int torang=0;
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
    ConjugateGradients(2500);
  }
  
  void OBForceField::SystematicRotorSearch() 
  {
    OBRotorList rl;
    OBRotamerList rotamers;

    rl.Setup(_mol);
    rotamers.SetBaseCoordinateSets(_mol);
    rotamers.Setup(_mol, rl);
    
    IF_OBFF_LOGLVL_LOW {
      OBFFLog("\nS Y S T E M A T I C   R O T O R   S E A R C H\n\n");
      sprintf(logbuf, "  NUMBER OF ROTATABLE BONDS: %d\n", rl.Size());
      OBFFLog(logbuf);
    }

    if (!rl.Size()) { // only one conformer
      IF_OBFF_LOGLVL_LOW
        OBFFLog("  GENERATED ONLY ONE CONFORMER\n\n");
 
      ConjugateGradients(2500); // final energy minimizatin for best conformation
      
      return;
    }

    std::vector<int> rotorKey(rl.Size() + 1, 0); // indexed from 1

    OBRotorIterator ri;
    OBRotor *rotor = rl.BeginRotor(ri);
    for (int i = 1; i < rl.Size() + 1; ++i, rotor = rl.NextRotor(ri)) { // foreach rotor
      for (unsigned int j = 0; j < rotor->GetResolution().size(); j++) { // foreach torsion
        rotorKey[i] = j;
        rotamers.AddRotamer(rotorKey);
      }
    }

    rotamers.ExpandConformerList(_mol, _mol.GetConformers());
      
    IF_OBFF_LOGLVL_LOW {
      sprintf(logbuf, "  GENERATED %d CONFORMERS\n\n", _mol.NumConformers());
      OBFFLog(logbuf);
      OBFFLog("CONFORMER     ENERGY\n");
      OBFFLog("--------------------\n");
    }
    
    // Calculate energy for all conformers
    char logbuf[100];
    std::vector<double> energies(_mol.NumConformers(), 0.0);
    for (int i = 0; i < _mol.NumConformers(); i++) {
      _mol.SetConformer(i); // select conformer

      // Need a minimization step here, actually

      energies[i] = Energy(); // calculate and store energy
      
      IF_OBFF_LOGLVL_LOW {
        sprintf(logbuf, "   %3d      %8.3f\n", (i + 1), energies[i]);
        OBFFLog(logbuf);
      }
    }

    // Select conformer with lowest energy
    int best_conformer = 0;
    for (int i = 1; i < _mol.NumConformers(); i++) {
      if (energies[i] < energies[best_conformer])
        best_conformer = i;
    }
  
    IF_OBFF_LOGLVL_LOW {
      sprintf(logbuf, "\n  CONFORMER %d HAS THE LOWEST ENERGY\n\n",  best_conformer + 1);
      OBFFLog(logbuf);
    }

    _mol.SetConformer(best_conformer);
    current_conformer = best_conformer;
  }

  void OBForceField::DistanceGeometry() 
  {
    int N = _mol.NumAtoms();
    int i = 0;
    int j = 0;
    //double matrix[N][N], G[N][N];
    vector<vector<double> > matrix(N);
    for(int i=0; i<N; i++)
      matrix[i].reserve(N);
    vector<vector<double> > G(N);
    for(int i=0; i<N; i++)
      G[i].reserve(N);
	
    bool is15;
    
    IF_OBFF_LOGLVL_LOW
      OBFFLog("\nD I S T A N C E   G E O M E T R Y\n\n");
 
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
    //           C \       .
    //              d      ad = sqrt(delta_x^2 + delta_y^2)
    //
    FOR_ATOMS_OF_MOL (a, _mol) {
      i = a->GetIdx() - 1;
      FOR_ATOMS_OF_MOL (b, _mol) {
        j = b->GetIdx() - 1;

        if (&*a == &*b) {
          matrix[i][j] = 0.0; // diagonal
          continue;
        }
        // Find relationship
        is15 = true;
        FOR_NBORS_OF_ATOM (nbr1, _mol.GetAtom(a->GetIdx())) { // 1-2
          if (&*nbr1 == &*b) {
            matrix[i][j] = 1.3;
            break;
          }
          FOR_NBORS_OF_ATOM (nbr2, _mol.GetAtom(nbr1->GetIdx())) { // 1-3
            if (&*nbr2 == &*b) {
              matrix[i][j] = sqrt(1.3*1.3 + 1.3*1.3 - 2.0 * 
                                  cos(DEG_TO_RAD*120.0) * 1.3*1.3 );
              is15 = false;
              break;
            }
            FOR_NBORS_OF_ATOM (nbr3, &*nbr2) { // 1-4
              if (&*nbr3 == &*b) {
                is15 = false;
                if (i > j) // minimum distance (torsion angle = 0)
                  matrix[i][j] = 1.3 + 1.3*cos(DEG_TO_RAD*(180.0 - 120.0)) 
                    + 1.3*cos(DEG_TO_RAD*(180.0 - 120.0));
                else {// maximum distance (torsion angle = 180)
                  double delta_x, delta_y;
                  delta_x = 1.3 + 1.3*cos(DEG_TO_RAD*(180.0-120.0)) 
                    + 1.3*cos(DEG_TO_RAD*(180.0-120.0));
                  delta_y = 1.3*sin(DEG_TO_RAD*(180.0-120.0)) + 1.3*sin(DEG_TO_RAD*(180.0-120.0));
                  matrix[i][j] = sqrt(delta_x*delta_x + delta_y*delta_y);
                }
                break;
              }
              if (i > j && is15) {// minimum distance (sum vdw radii)
                matrix[i][j] = 1.4 + 1.4;
              } else if (is15) // maximum distance (torsion angle = 180)
                matrix[i][j] = 99.0;
            }
          }
        }
      }
    }
   
    // output initial distance matrix
    IF_OBFF_LOGLVL_LOW {

      OBFFLog("INITIAL DISTANCE MATRIX\n\n");
      for (i=0; i<N; i++) {
        OBFFLog("[");
        for (j=0; j<N; j++) {
          sprintf(logbuf, " %8.4f ", matrix[i][j]);
          OBFFLog(logbuf);
        }
        OBFFLog("]\n");
      }
      OBFFLog("\n");
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

      OBFFLog("TRIANGLE SMOOTHING\n\n");
      for (i=0; i<N; i++) {
        OBFFLog("[");
        for (j=0; j<N; j++) {
          sprintf(logbuf, " %8.4f ", matrix[i][j]);
          OBFFLog(logbuf);
        }
        OBFFLog("]\n");
      }
      OBFFLog("\n");
    }
    
    // Generate random distance matrix between lower and upper limits
    FOR_ATOMS_OF_MOL (a, _mol) {
      i = a->GetIdx() - 1;
      FOR_ATOMS_OF_MOL (b, _mol) {
        j = b->GetIdx() - 1;

        if (&*a == &*b) {
          matrix[i][j] = 0.0; // diagonal
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

      OBFFLog("RANDOM DISTANCE MATRIX BETWEEN LIMITS\n\n");
      for (i=0; i<N; i++) {
        OBFFLog("[");
        for (j=0; j<N; j++) {
          sprintf(logbuf, " %8.4f ", matrix[i][j]);
          OBFFLog(logbuf);
        }
        OBFFLog("]\n");
      }
      OBFFLog("\n");
    }

    // Generate metric matrix
    // (origin = first atom )
    for (i=0; i<N; i++) {
      for (j=0; j<N; j++) {
        G[i][j] = 0.5 * (matrix[0][i]*matrix[0][i] + matrix[0][j]*matrix[0][j] - matrix[i][j]*matrix[i][j]);
      }
    }
    
    // output metric matrix
    IF_OBFF_LOGLVL_LOW {
      char logbuf[100];

      OBFFLog("METRIC MATRIX\n\n");
      for (i=0; i<N; i++) {
        OBFFLog("[");
        for (j=0; j<N; j++) {
          sprintf(logbuf, " %8.4f ", G[i][j]);
          OBFFLog(logbuf);
        }
        OBFFLog("]\n");
      }
      OBFFLog("\n");
    }

    // Calculate eigenvalues and eigenvectors
    //double eigenvalues[N];
    vector<double> eigenvalues(N);
    //double eigenvectors[N][N];
    vector<vector<double> > eigenvectors(N);
    for(int i=0; i<N; i++)
      eigenvectors[i].reserve(N);

    matrix3x3::jacobi(N, (double *) &G, (double *) &eigenvalues, (double *) &eigenvectors);
    
    // output eigenvalues and eigenvectors
    IF_OBFF_LOGLVL_LOW {

      OBFFLog("EIGENVALUES OF METRIC MATRIX\n\n");
      for (i=0; i<N; i++) {
        sprintf(logbuf, "%8.4f ", eigenvalues[i]);
        OBFFLog(logbuf);
      }
      OBFFLog("\n");

      OBFFLog("EIGENVECTORS OF METRIC MATRIX\n\n");
      for (i=0; i<N; i++) {
        OBFFLog("[");
        for (j=0; j<N; j++) {
          sprintf(logbuf, " %8.4f ", eigenvectors[i][j]);
          OBFFLog(logbuf);
        }
        OBFFLog("]\n");
      }
      OBFFLog("\n");
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
    vector3 old_xyz, orig_xyz, xyz_k, dir(0.0, 0.0, 0.0);

    // This needs several enhancements
    // 1) Switch to smarter line search method (e.g., Newton's method in 1D)
    //  x(n+1) = x(n) - F(x) / F'(x)   -- can be done numerically
    // 2) Switch to line search for one step of the entire molecule!
    //   This dramatically cuts down on the number of Energy() calls.
    //  (and is more correct anyway)

    step = 0.2;
    direction.normalize();
    orig_xyz = atom->GetVector();
    
    e_n1 = Energy(false); // calculate e_k
    
    unsigned int i;
    for (i=0; i<100; i++) {
      old_xyz = atom->GetVector();
      
      xyz_k = atom->GetVector() + direction*step;
      atom->SetVector(xyz_k);  // update coordinates
    
      e_n2 = Energy(false); // calculate e_k+1
      
      // convergence criteria: A higher precision here 
      // only takes longer with the same result.
      if (IsNear(e_n2, e_n1, 1.0e-3))
        break;

      if (e_n2 > e_n1) { // decrease stepsize
        step *= 0.1;
        atom->SetVector(old_xyz);
      }
      if (e_n2 < e_n1) {  // increase stepsize
        e_n1 = e_n2;
        step *= 2.15;
        if (step > 1.0)
          step = 1.0;
      }
      
    }
    //cout << "LineSearch steps: " << i << endl;

    dir = atom->GetVector() - orig_xyz;
    atom->SetVector(orig_xyz);     

    // cutoff accuracy
    if (dir.length() < 1.0e-8)
      return VZero;

    return dir;
  }
  
  vector3 OBForceField::ValidateLineSearch(OBAtom *atom, vector3 &direction)
  {
    double e_n1, e_n2, step;
    vector3 old_xyz, orig_xyz, xyz_k, dir(0.0, 0.0, 0.0);

    step = 0.2;
    direction.normalize();
    orig_xyz = atom->GetVector();
    
    e_n1 = atom->x() * atom->x() + 2 * (atom->y() * atom->y()); // e_k
    
    for (unsigned int i=0; i<100; i++) {
      old_xyz = atom->GetVector();
      
      xyz_k = atom->GetVector() + direction*step;
      atom->SetVector(xyz_k);  // update coordinates
    
      e_n2 = atom->x() * atom->x() + 2 * (atom->y() * atom->y()); // e_k+1
      
      if (e_n2 == e_n1) // convergence criteria
        break;

      if (e_n2 > e_n1) { // decrease stepsize
        step *= 0.5;
        atom->SetVector(old_xyz);
      }
      if (e_n2 < e_n1) {  // increase stepsize
        e_n1 = e_n2;
        step *= 1.2;
        if (step > 1.0)
          step = 1.0;
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
    
    atom->SetVector(9.0, 9.0, 0.0);
    e_n1 = atom->x() * atom->x() + 2 * (atom->y() * atom->y());
    
    IF_OBFF_LOGLVL_LOW {
      OBFFLog("\nV A L I D A T E   S T E E P E S T   D E S C E N T\n\n");
      sprintf(logbuf, "STEPS = %d\n\n",  steps);
      OBFFLog(logbuf);
      OBFFLog("STEP n     E(n)       E(n-1)    \n");
      OBFFLog("--------------------------------\n");
    }

    for (int i = 1; i <= steps; i++) {
      grad.Set(-2*atom->x(), -4*atom->y(), 0.0);
      grad = ValidateLineSearch(&*atom, grad);
      atom->SetVector(atom->x() + grad.x(), atom->y() + grad.y(), 0.0);
      e_n2 = atom->x() * atom->x() + 2 * (atom->y() * atom->y());
      
      IF_OBFF_LOGLVL_LOW {
        sprintf(logbuf, " %4d    %8.3f    %8.3f\n", i, e_n2, e_n1);
        OBFFLog(logbuf);
      }

      if (IsNear(e_n2, e_n1, 1.0e-7)) {
        IF_OBFF_LOGLVL_LOW
          OBFFLog("    STEEPEST DESCENT HAS CONVERGED (DELTA E < 1.0e-7)\n");
        break;
      }

      e_n1 = e_n2;
    }

    IF_OBFF_LOGLVL_LOW
      OBFFLog("\n");
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
    atom->SetVector(9.0, 9.0, 0.0);
    e_n1 = atom->x() * atom->x() + 2 * (atom->y() * atom->y());
 
    IF_OBFF_LOGLVL_LOW {
      OBFFLog("\nV A L I D A T E   C O N J U G A T E   G R A D I E N T S\n\n");
      sprintf(logbuf, "STEPS = %d\n\n",  steps);
      OBFFLog(logbuf);
      OBFFLog("STEP n     E(n)       E(n-1)    \n");
      OBFFLog("--------------------------------\n");
    }
 
    for (int i = 1; i <= steps; i++) {
      if (firststep) {
        grad1.Set(-2*atom->x(), -4*atom->y(), 0.0);
        dir1 = grad1;
        dir1 = ValidateLineSearch(&*atom, dir1);
        atom->SetVector(atom->x() + dir1.x(), atom->y() + dir1.y(), atom->z() + dir1.z());
        e_n2 = atom->x() * atom->x() + 2 * (atom->y() * atom->y());
      
        IF_OBFF_LOGLVL_LOW {
          sprintf(logbuf, " %4d    %8.3f    %8.3f\n", i, e_n2, e_n1);
          OBFFLog(logbuf);
        }
        
        e_n1 = e_n2;
        dir1 = grad1;
        firststep = false;
      } else {
        grad2.Set(-2*atom->x(), -4*atom->y(), 0.0);
        g2g2 = dot(grad2, grad2);
        g1g1 = dot(grad1, grad1);
        g2g1 = g2g2 / g1g1;
        dir2 = grad2 + g2g1 * dir1;
        dir2 = ValidateLineSearch(&*atom, dir2);
        atom->SetVector(atom->x() + dir2.x(), atom->y() + dir2.y(), atom->z() + dir2.z());
        grad1 = grad2;
        dir1 = dir2;
        e_n2 = atom->x() * atom->x() + 2 * (atom->y() * atom->y());
	  
        IF_OBFF_LOGLVL_LOW {
          sprintf(logbuf, " %4d    %8.3f    %8.3f\n", i, e_n2, e_n1);
          OBFFLog(logbuf);
        }
        
        if (IsNear(e_n2, e_n1, 1.0e-7)) {
          IF_OBFF_LOGLVL_LOW
            OBFFLog("    CONJUGATE GRADIENTS HAS CONVERGED (DELTA E < 1.0e-7)\n");
          break;
        }

        e_n1 = e_n2;
      }
    }
  }
  
  void OBForceField::SteepestDescentInitialize(int steps, double econv, int method) 
  {
    _nsteps = steps;
    _cstep = 0;
    _econv = econv;
    _method = method;

    _e_n1 = Energy();
    
    IF_OBFF_LOGLVL_LOW {
      OBFFLog("\nS T E E P E S T   D E S C E N T\n\n");
      sprintf(logbuf, "STEPS = %d\n\n",  steps);
      OBFFLog(logbuf);
      OBFFLog("STEP n     E(n)       E(n-1)    \n");
      OBFFLog("--------------------------------\n");
    }
  }
 
  bool OBForceField::SteepestDescentTakeNSteps(int n) 
  {
    double e_n2;
    vector3 grad;

    for (int i = 1; i <= n; i++) {
      _cstep++;
      
      FOR_ATOMS_OF_MOL (a, _mol) {
        if (_method & OBFF_ANALYTICAL_GRADIENT)
          grad = GetGradient(&*a);
        else
          grad = NumericalDerivative(&*a);
        grad = LineSearch(&*a, grad);
        a->SetVector(a->x() + grad.x(), a->y() + grad.y(), a->z() + grad.z());
      }
      e_n2 = Energy();
      
      IF_OBFF_LOGLVL_LOW {
        sprintf(logbuf, " %4d    %8.3f    %8.3f\n", i, e_n2, _e_n1);
        OBFFLog(logbuf);
      }

      if (IsNear(e_n2, _e_n1, _econv)) {
        IF_OBFF_LOGLVL_LOW
          OBFFLog("    STEEPEST DESCENT HAS CONVERGED\n");
        return false;
      }
      
      if (_nsteps == _cstep)
        return false;

      _e_n1 = e_n2;
    }

    return true;  // no convergence reached
  }
 
  void OBForceField::SteepestDescent(int steps, double econv, int method) 
  {
    SteepestDescentInitialize(steps, econv, method);
    SteepestDescentTakeNSteps(steps);
  }
  
  void OBForceField::ConjugateGradientsInitialize(int steps, double econv, 
                                                  int method)
  {
    double e_n2;
    vector3 grad2, dir2;

    _cstep = 1;
    _nsteps = steps;
    _econv = econv;
    _method = method;

    ValidateGradients();

    _e_n1 = Energy();
    
    IF_OBFF_LOGLVL_LOW {
      OBFFLog("\nC O N J U G A T E   G R A D I E N T S\n\n");
      sprintf(logbuf, "STEPS = %d\n\n",  steps);
      OBFFLog(logbuf);
      OBFFLog("STEP n     E(n)       E(n-1)    \n");
      OBFFLog("--------------------------------\n");
    }

    _grad1.resize(_mol.NumAtoms() + 1);
    _dir1.resize(_mol.NumAtoms() + 1);

    // Take the first step (same as steepest descent because there is no 
    // gradient from the previous step.
    FOR_ATOMS_OF_MOL (a, _mol) {
      if (_method & OBFF_ANALYTICAL_GRADIENT)
        grad2 = GetGradient(&*a);
      else
        grad2 = NumericalDerivative(&*a);
      dir2 = grad2;
      dir2 = LineSearch(&*a, dir2);
      a->SetVector(a->x() + dir2.x(), a->y() + dir2.y(), a->z() + dir2.z());
      _grad1[a->GetIdx()] = grad2;
      _dir1[a->GetIdx()] = grad2;
    }
    e_n2 = Energy();
      
    IF_OBFF_LOGLVL_LOW {
      sprintf(logbuf, " %4d    %8.3f    %8.3f\n", _cstep, e_n2, _e_n1);
      OBFFLog(logbuf);
    }
 
    _e_n1 = e_n2;
  }
  
  bool OBForceField::ConjugateGradientsTakeNSteps(int n)
  {
    double e_n2;
    double g2g2, g1g1, g2g1;
    vector3 grad2, dir2;
    
    int printStep = min(n, 10); // print results every 10 steps, or fewer

    if (_grad1.size() != (_mol.NumAtoms()+1))
      return false;

    e_n2 = 0.0;
    
    for (int i = 1; i <= n; i++) {
      _cstep++;
      
      FOR_ATOMS_OF_MOL (a, _mol) {
        if (_method & OBFF_ANALYTICAL_GRADIENT)
          grad2 = GetGradient(&*a);
        else
          grad2 = NumericalDerivative(&*a);
        
        // Fletcher-Reeves formula for Beta
        // http://en.wikipedia.org/wiki/Nonlinear_conjugate_gradient_method
        // NOTE: We make sure to reset and use the steepest descent direction
        //   after NumAtoms steps
        if (_cstep % _mol.NumAtoms() != 0) {
        g2g2 = dot(grad2, grad2);
        g1g1 = dot(_grad1[a->GetIdx()], _grad1[a->GetIdx()]);
        g2g1 = g2g2 / g1g1;
        dir2 = grad2 + g2g1 * _dir1[a->GetIdx()];
        } else { // reset conj. direction
          dir2 = grad2;
        }

        dir2 = LineSearch(&*a, dir2);
        a->SetVector(a->x() + dir2.x(), a->y() + dir2.y(), a->z() + dir2.z());
	  
        _grad1[a->GetIdx()] = grad2;
        _dir1[a->GetIdx()] = dir2;
        if (e_n2)
          _e_n1 = e_n2;
      }
      e_n2 = Energy();
	
      IF_OBFF_LOGLVL_LOW {
        if (_cstep % 10 == 0) {
        sprintf(logbuf, " %4d    %8.3f    %8.3f\n", _cstep, e_n2, _e_n1);
        OBFFLog(logbuf);
      }
      }
 
      if (IsNear(e_n2, _e_n1, _econv)) {
        IF_OBFF_LOGLVL_LOW
          OBFFLog("    CONJUGATE GRADIENTS HAS CONVERGED\n");
        return false;
      }

      if (_nsteps == _cstep)
        return false;

      _e_n1 = e_n2;
    }

    return true; // no convergence reached
  }
 
  void OBForceField::ConjugateGradients(int steps, double econv, int method)
  {
    ConjugateGradientsInitialize(steps, econv, method);
    ConjugateGradientsTakeNSteps(steps); // ConjugateGradientsInitialize takes the fisrt step
  }
  
  //  
  //         f(1) - f(0)
  // f'(0) = -----------      f(1) = f(0+h)
  //              h
  //
  vector3 OBForceField::NumericalDerivative(OBAtom *atom, int terms)
  {
    vector3 va, grad;
    double e_orig, e_plus_delta, delta, dx, dy, dz;

    delta = 1.0e-5;

    va = atom->GetVector();

    if (terms & OBFF_ENERGY)
      e_orig = Energy(false);
    else {
      e_orig = 0.0;
      if (terms & OBFF_EBOND)
        e_orig += E_Bond();
      if (terms & OBFF_EANGLE)
        e_orig += E_Angle();
      if (terms & OBFF_ESTRBND)
        e_orig += E_StrBnd();
      if (terms & OBFF_ETORSION)
        e_orig += E_Torsion();
      if (terms & OBFF_EOOP)
        e_orig += E_OOP();
      if (terms & OBFF_EVDW)
        e_orig += E_VDW();
      if (terms & OBFF_EELECTROSTATIC)
        e_orig += E_Electrostatic();
    }
    
    // X direction
    atom->SetVector(va.x() + delta, va.y(), va.z());

    if (terms & OBFF_ENERGY)
      e_plus_delta = Energy(false);
    else {
      e_plus_delta = 0.0;
      if (terms & OBFF_EBOND)
        e_plus_delta += E_Bond();
      if (terms & OBFF_EANGLE)
        e_plus_delta += E_Angle();
      if (terms & OBFF_ESTRBND)
        e_plus_delta += E_StrBnd();
      if (terms & OBFF_ETORSION)
        e_plus_delta += E_Torsion();
      if (terms & OBFF_EOOP)
        e_plus_delta += E_OOP();
      if (terms & OBFF_EVDW)
        e_plus_delta += E_VDW();
      if (terms & OBFF_EELECTROSTATIC)
        e_plus_delta += E_Electrostatic();
    }
    
    dx = (e_plus_delta - e_orig) / delta;
    
    // Y direction
    atom->SetVector(va.x(), va.y() + delta, va.z());
    
    if (terms & OBFF_ENERGY)
      e_plus_delta = Energy(false);
    else {
      e_plus_delta = 0.0;
      if (terms & OBFF_EBOND)
        e_plus_delta += E_Bond();
      if (terms & OBFF_EANGLE)
        e_plus_delta += E_Angle();
      if (terms & OBFF_ESTRBND)
        e_plus_delta += E_StrBnd();
      if (terms & OBFF_ETORSION)
        e_plus_delta += E_Torsion();
      if (terms & OBFF_EOOP)
        e_plus_delta += E_OOP();
      if (terms & OBFF_EVDW)
        e_plus_delta += E_VDW();
      if (terms & OBFF_EELECTROSTATIC)
        e_plus_delta += E_Electrostatic();
    }
 
    dy = (e_plus_delta - e_orig) / delta;
    
    // Z direction
    atom->SetVector(va.x(), va.y(), va.z() + delta);
    
    if (terms & OBFF_ENERGY)
      e_plus_delta = Energy(false);
    else {
      e_plus_delta = 0.0;
      if (terms & OBFF_EBOND)
        e_plus_delta += E_Bond();
      if (terms & OBFF_EANGLE)
        e_plus_delta += E_Angle();
      if (terms & OBFF_ESTRBND)
        e_plus_delta += E_StrBnd();
      if (terms & OBFF_ETORSION)
        e_plus_delta += E_Torsion();
      if (terms & OBFF_EOOP)
        e_plus_delta += E_OOP();
      if (terms & OBFF_EVDW)
        e_plus_delta += E_VDW();
      if (terms & OBFF_EELECTROSTATIC)
        e_plus_delta += E_Electrostatic();
    }
 
    dz = (e_plus_delta - e_orig) / delta;

    // reset coordinates to original
    atom->SetVector(va.x(), va.y(), va.z());

    grad.Set(-dx, -dy, -dz);
    return (grad);
  }
  
  //  
  //         f(2) - 2f(1) + f(0)
  // f'(0) = -------------------      f(1) = f(0+h)
  //                 h^2              f(1) = f(0+2h)
  //
  vector3 OBForceField::NumericalSecondDerivative(OBAtom *atom, int terms)
  {
    vector3 va, grad;
    double e_0, e_1, e_2, delta, dx, dy, dz;

    delta = 1.0e-5;

    va = atom->GetVector();

    // calculate f(0)
    if (terms & OBFF_ENERGY)
      e_0 = Energy(false);
    else {
      e_0 = 0.0;
      if (terms & OBFF_EBOND)
        e_0 += E_Bond(false);
      if (terms & OBFF_EANGLE)
        e_0 += E_Angle(false);
      if (terms & OBFF_ESTRBND)
        e_0 += E_StrBnd(false);
      if (terms & OBFF_ETORSION)
        e_0 += E_Torsion(false);
      if (terms & OBFF_EOOP)
        e_0 += E_OOP(false);
      if (terms & OBFF_EVDW)
        e_0 += E_VDW(false);
      if (terms & OBFF_EELECTROSTATIC)
        e_0 += E_Electrostatic(false);
    }
    
    // 
    // X direction
    //
    
    // calculate f(1)
    atom->SetVector(va.x() + delta, va.y(), va.z());

    if (terms & OBFF_ENERGY)
      e_1 = Energy(false);
    else {
      e_1 = 0.0;
      if (terms & OBFF_EBOND)
        e_1 += E_Bond(false);
      if (terms & OBFF_EANGLE)
        e_1 += E_Angle(false);
      if (terms & OBFF_ESTRBND)
        e_1 += E_StrBnd(false);
      if (terms & OBFF_ETORSION)
        e_1 += E_Torsion(false);
      if (terms & OBFF_EOOP)
        e_1 += E_OOP(false);
      if (terms & OBFF_EVDW)
        e_1 += E_VDW(false);
      if (terms & OBFF_EELECTROSTATIC)
        e_1 += E_Electrostatic(false);
    }
    
    // calculate f(2)
    atom->SetVector(va.x() + 2 * delta, va.y(), va.z());

    if (terms & OBFF_ENERGY)
      e_2 = Energy(false);
    else {
      e_2 = 0.0;
      if (terms & OBFF_EBOND)
        e_2 += E_Bond(false);
      if (terms & OBFF_EANGLE)
        e_2 += E_Angle(false);
      if (terms & OBFF_ESTRBND)
        e_2 += E_StrBnd(false);
      if (terms & OBFF_ETORSION)
        e_2 += E_Torsion(false);
      if (terms & OBFF_EOOP)
        e_2 += E_OOP(false);
      if (terms & OBFF_EVDW)
        e_2 += E_VDW(false);
      if (terms & OBFF_EELECTROSTATIC)
        e_2 += E_Electrostatic(false);
    }
    
    dx = (e_2 - 2 * e_1 + e_0) / (delta * delta);
    
    // 
    // Y direction
    //
    
    // calculate f(1)
    atom->SetVector(va.x(), va.y() + delta, va.z());

    if (terms & OBFF_ENERGY)
      e_1 = Energy(false);
    else {
      e_1 = 0.0;
      if (terms & OBFF_EBOND)
        e_1 += E_Bond(false);
      if (terms & OBFF_EANGLE)
        e_1 += E_Angle(false);
      if (terms & OBFF_ESTRBND)
        e_1 += E_StrBnd(false);
      if (terms & OBFF_ETORSION)
        e_1 += E_Torsion(false);
      if (terms & OBFF_EOOP)
        e_1 += E_OOP(false);
      if (terms & OBFF_EVDW)
        e_1 += E_VDW(false);
      if (terms & OBFF_EELECTROSTATIC)
        e_1 += E_Electrostatic(false);
    }
    
    // calculate f(2)
    atom->SetVector(va.x(), va.y() + 2 * delta, va.z());

    if (terms & OBFF_ENERGY)
      e_2 = Energy(false);
    else {
      e_2 = 0.0;
      if (terms & OBFF_EBOND)
        e_2 += E_Bond(false);
      if (terms & OBFF_EANGLE)
        e_2 += E_Angle(false);
      if (terms & OBFF_ESTRBND)
        e_2 += E_StrBnd(false);
      if (terms & OBFF_ETORSION)
        e_2 += E_Torsion(false);
      if (terms & OBFF_EOOP)
        e_2 += E_OOP(false);
      if (terms & OBFF_EVDW)
        e_2 += E_VDW(false);
      if (terms & OBFF_EELECTROSTATIC)
        e_2 += E_Electrostatic(false);
    }
    
    dy = (e_2 - 2 * e_1 + e_0) / (delta * delta);

    // 
    // Z direction
    //
    
    // calculate f(1)
    atom->SetVector(va.x(), va.y(), va.z() + delta);

    if (terms & OBFF_ENERGY)
      e_1 = Energy(false);
    else {
      e_1 = 0.0;
      if (terms & OBFF_EBOND)
        e_1 += E_Bond(false);
      if (terms & OBFF_EANGLE)
        e_1 += E_Angle(false);
      if (terms & OBFF_ESTRBND)
        e_1 += E_StrBnd(false);
      if (terms & OBFF_ETORSION)
        e_1 += E_Torsion(false);
      if (terms & OBFF_EOOP)
        e_1 += E_OOP(false);
      if (terms & OBFF_EVDW)
        e_1 += E_VDW(false);
      if (terms & OBFF_EELECTROSTATIC)
        e_1 += E_Electrostatic(false);
    }
    
    // calculate f(2)
    atom->SetVector(va.x(), va.y(), va.z() + 2 * delta);

    if (terms & OBFF_ENERGY)
      e_2 = Energy(false);
    else {
      e_2 = 0.0;
      if (terms & OBFF_EBOND)
        e_2 += E_Bond(false);
      if (terms & OBFF_EANGLE)
        e_2 += E_Angle(false);
      if (terms & OBFF_ESTRBND)
        e_2 += E_StrBnd(false);
      if (terms & OBFF_ETORSION)
        e_2 += E_Torsion(false);
      if (terms & OBFF_EOOP)
        e_2 += E_OOP(false);
      if (terms & OBFF_EVDW)
        e_2 += E_VDW(false);
      if (terms & OBFF_EELECTROSTATIC)
        e_2 += E_Electrostatic(false);
    }
    
    dz = (e_2 - 2 * e_1 + e_0) / (delta * delta);

    // reset coordinates to original
    atom->SetVector(va.x(), va.y(), va.z());

    grad.Set(-dx, -dy, -dz);
    return (grad);
  }

  bool OBForceField::UpdateCoordinates(OBMol &mol)
  { 
    OBAtom *atom;
    
    if (_mol.NumAtoms() != mol.NumAtoms())
      return false;
    
    // Copy coordinates for current conformer only
    FOR_ATOMS_OF_MOL (a, _mol) {
      atom = mol.GetAtom(a->GetIdx());
      atom->SetVector(a->GetVector());
    }

    return true;
  }
  
  bool OBForceField::UpdateConformers(OBMol &mol)
  { 
    OBAtom *atom;

    if (_mol.NumAtoms() != mol.NumAtoms())
      return false;
    
    FOR_ATOMS_OF_MOL (a, _mol) {
      atom = mol.GetAtom(a->GetIdx());
      atom->SetVector(a->GetVector());
    }

    //Copy conformer information
    if (_mol.NumConformers() > 1) {
      int k,l;
      vector<double*> conf;
      double* xyz = NULL;
      for (k=0 ; k<_mol.NumConformers() ; ++k) {
        xyz = new double [3*_mol.NumAtoms()];
        for (l=0 ; l<(int) (3*_mol.NumAtoms()) ; ++l)
          xyz[l] = _mol.GetConformer(k)[l];
        conf.push_back(xyz);
      }
      mol.SetConformers(conf);
      mol.SetConformer(current_conformer);
    }
    
    return true;
  }

  vector3 OBForceField::ValidateGradientError(vector3 &numgrad, vector3 &anagrad)
  {
    double errx, erry, errz;

    if (fabs(numgrad.x()) < 1.0)
      errx = numgrad.x() * fabs(numgrad.x() - anagrad.x()) * 100;
    else
      errx = fabs((numgrad.x() - anagrad.x()) / numgrad.x()) * 100;

    if (fabs(numgrad.y()) < 1.0)
      erry = numgrad.y() * fabs(numgrad.y() - anagrad.y()) * 100;
    else
      erry = fabs((numgrad.y() - anagrad.y()) / numgrad.y()) * 100;
    
    if (fabs(numgrad.z()) < 1.0)
      errz = numgrad.z() * fabs(numgrad.z() - anagrad.z()) * 100;
    else
      errz = fabs((numgrad.z() - anagrad.z()) / numgrad.z()) * 100;
    
    errx = fabs(errx);
    erry = fabs(erry);
    errz = fabs(errz);

    return vector3(errx, erry, errz);
  }
  
  double OBForceField::VectorLengthDerivative(vector3 &a, vector3 &b)
  {
    vector3 vab, drab;
    double rab;
    
    vab = a - b;
    rab = vab.length();
    if (rab < 0.1) // atoms are too close to each other
      {
        vab.randomUnitVector();
        vab *= 0.1; // move the atoms a small distance apart
        rab = 0.1;
      }
    drab = vab / rab;

    a = -drab; // -drab/da
    b =  drab; // -drab/db

    return rab;
  }
  
  double OBForceField::VectorAngleDerivative(vector3 &a, vector3 &b, vector3 &c)
  {
    vector3 vab, vcb;
    double theta, rab, rab2, rcb, rcb2, abcb, abcb2;
     
    vab = a - b;
    vcb = c - b;
    rab = vab.length();
    rcb = vcb.length();
    rab2 = rab * rab;
    rcb2 = rcb * rcb;
    
    abcb = dot(vab, vcb) / (rab * rcb);
    abcb2 = 1.0 - abcb * abcb;
    theta = acos(abcb) * RAD_TO_DEG;
    
    if (IsNearZero(abcb2)) {
      a = VZero;
      b = VZero;
      c = VZero;
      return 0.0;
    }

    a = (vcb * rab * rcb - (vab / rab) * dot(vab, vcb) * rcb) / (sqrt(abcb2) * rab2 * rcb2);
    c = -((vcb / rcb) * dot(vab, vcb) * rab - vab * rab * rcb) / (sqrt(abcb2) * rab2 * rcb2);
    b = -a - c;

    a *= (1.0 / DEG_TO_RAD);
    b *= (1.0 / DEG_TO_RAD);
    c *= (1.0 / DEG_TO_RAD);

    return theta;
  }
 
  double OBForceField::VectorTorsionDerivative(vector3 &a, vector3 &b, vector3 &c, vector3 &d)
  {
    vector3 vab, vbc, vcd, vac, vbd, grada, gradb, gradc, gradd;
    vector3 abbc, bccd, bcabbc, bcbccd, cdabbc, cdbccd, acabbc, acbccd, ababbc, abbccd, bdabbc, bdbccd;
    double tor, rabbc, rbccd, rabbc2, rbccd2, rabbc3, rbccd3, abbc_bccd, rabbc_rbccd, abbc_bccd2;

    vab = a - b;
    vbc = b - c;
    vcd = c - d;
    vac = a - c;
    vbd = b - d;
    abbc = cross(vab, vbc);
    bccd = cross(vbc, vcd);
    
    //    cerr << " abbc.length " << abbc.length() << " bccd length " << bccd.length()
    //         << endl;

    double dotAbbcBccd = dot(abbc,bccd);
    tor = RAD_TO_DEG * acos(dotAbbcBccd / (abbc.length() * bccd.length()));
    if (IsNearZero(dotAbbcBccd)) {
      tor = 0.0;
    }
    else if (dotAbbcBccd > 0.0) {
      tor = -tor;
    }
 
    bcabbc = cross(vbc, abbc);
    bcbccd = cross(vbc, bccd);
    cdabbc = cross(vcd, abbc);
    cdbccd = cross(vcd, bccd);
    acabbc = cross(vac, abbc);
    acbccd = cross(vac, bccd);
    ababbc = cross(vab, abbc);
    abbccd = cross(vab, bccd);
    bdabbc = cross(vbd, abbc);
    bdbccd = cross(vbd, bccd);
    rabbc = abbc.length();
    rbccd = bccd.length();
    rabbc2 = rabbc * rabbc;
    rbccd2 = rbccd * rbccd;
    rabbc3 = rabbc2 * rabbc;
    rbccd3 = rbccd2 * rbccd;
    abbc_bccd = dot(abbc, bccd);
    rabbc_rbccd = abbc_bccd / (rabbc * rbccd);
    abbc_bccd2 = sqrt(1.0 - rabbc_rbccd * rabbc_rbccd);

    if (IsNearZero(abbc_bccd2) || IsNan(abbc_bccd2))
      abbc_bccd2 = 1e-7;
    
    /*
      cout << "sqrt(abbc_bccd2) = " << sqrt(abbc_bccd2) << endl;
      cout << "(rabbc*rbccd) = " << (rabbc*rbccd) << endl;
      cout << "dot(abbc,bccd) = " << dot(abbc,bccd) << endl;
      cout << "(rabbc*rbccd3) = " << (rabbc*rbccd3) << endl;
      cout << "(rabbc3*rbccd) = " << (rabbc3*rbccd) << endl;
    */
    
    a = (bcbccd / (rabbc*rbccd) - (bcabbc*abbc_bccd) / (rabbc3*rbccd)) / abbc_bccd2;
    d = (bcabbc / (rabbc*rbccd) - (bcbccd*abbc_bccd) / (rabbc*rbccd3)) / abbc_bccd2;

    b = ( -(cdbccd*abbc_bccd) / (rabbc*rbccd3) + 
          (cdabbc - acbccd) / (rabbc*rbccd) + 
          (acabbc*abbc_bccd) / (rabbc3*rbccd)    ) / abbc_bccd2;
    
    c = (  (bdbccd*abbc_bccd) / (rabbc*rbccd3) + 
           (abbccd - bdabbc) / (rabbc*rbccd) +
           -(ababbc*abbc_bccd) / (rabbc3*rbccd)    ) / abbc_bccd2;

    if (IsNan(a.x()) || IsNan(a.y()) || IsNan(a.z()))
      a = VZero;
    if (IsNan(b.x()) || IsNan(b.y()) || IsNan(b.z()))
      b = VZero;
    if (IsNan(c.x()) || IsNan(c.y()) || IsNan(c.z()))
      c = VZero;
    if (IsNan(d.x()) || IsNan(d.y()) || IsNan(d.z()))
      d = VZero;

    if (dot(abbc, bccd) > 0.0) {
      a = -a;
      b = -b;
      c = -c;
      d = -d;
    }
   
    return tor;  
  }
  
  bool OBForceField::IsInSameRing(OBAtom* a, OBAtom* b)
  {
    bool a_in, b_in;
    vector<OBRing*> vr;
    vr = _mol.GetSSSR();
    
    vector<OBRing*>::iterator i;
    vector<int>::iterator j;
    
    for (i = vr.begin();i != vr.end();i++) {
      a_in = false;
      b_in = false;
      for(j = (*i)->_path.begin();j != (*i)->_path.end();j++) {
        if ((unsigned)(*j) == a->GetIdx())
          a_in = true;
        if ((unsigned)(*j) == b->GetIdx())
          b_in = true;
      }
      
      if (a_in && b_in)
        return true;
    }

    return false;
  }
 
} // end namespace OpenBabel


//! \file forcefield.cpp
//! \brief Handle OBForceField class

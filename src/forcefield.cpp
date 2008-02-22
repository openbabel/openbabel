/**********************************************************************
forcefield.cpp - Handle OBForceField class.

Copyright (C) 2006-2007 by Tim Vandermeersch <tim.vandermeersch@gmail.com>
Some portions Copyright (C) 2007-2008 by Geoffrey Hutchison
 
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
      - LineSearch(): finished
      - SteepestDescent(): finished
      - ConjugateGradients(): finished
      - GenerateCoordinates():  removed, use OBBuilder
      - SystematicRotorSearch(): finished
      - RandomRotorSearch(): finished
      - WeightedRotorSearch: finished
      - DistanceGeometry(): needs matrix operations

      - src/forcefields/forcefieldghemical.cpp
      - Atom typing: finished
      - Charges: finished
      - Energy terms: finished
      - Analytical gradients: finished
      - Validation: finished

      - src/forcefields/forcefieldmmff94.cpp
      - Atom typing: done.
      - Charges: done.
      - Energy terms:
	    - Bond: finished
	    - Angle: finished
	    - StrBnd: small errors
	    - Torsion: finished
	    - OOP: no analytical gradient
	    - VDW: finished
	    - Electrostatic: finished
      - Analytical gradients: finished (except OOP)
      - Validation: http://home.scarlet.be/timvdm/MMFF94_validation_output.gz

      - src/forcefields/forcefielduff.cpp
      - Energy terms: finished
      - OOP: needs validation
      - Gradients: need OOP gradient
      - Validation in progress...


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
      Classes derived from the OBForceField implement specific force fields (Ghemical, 
      MMFF94, UFF, ...).Other classes such as OBFFParameter, OBFFConstraint, 
      OBFFCalculation and its derived classes are only for internal use. As a user 
      interested in using the available force fields in Open Babel, you don't need 
      these classes. The rest of this short introduction is aimed at these users. For 
      information on how to implement additional force fields, see the wiki pages or 
      post your questions to the openbabel-devel mailing list.

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
      
      Minimize the structure in mol using steepest descent and fix the position of atom with index 1.
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
      
      // Set the constraints
      OBFFConstraints constraints;
      constraints.AddAtomConstraint(1);

      // We pass the constraints as argument for Setup()
      if (!pFF->Setup(mol, constraints)) {
        cerr << "ERROR: could not setup force field." << endl;
      }
      
      // Perform the actual minimization, maximum 1000 steps 
      pFF->SteepestDescent(1000);
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
  
  //////////////////////////////////////////////////////////////////////////////////
  //
  // Constraints
  //
  //////////////////////////////////////////////////////////////////////////////////
 
  OBFFConstraints::OBFFConstraints()
  {
    _factor = 50000.0;
  }
 
  void OBFFConstraints::Setup(OBMol& mol)
  {
    vector<OBFFConstraint>::iterator i;
    
    for (i = _constraints.begin(); i != _constraints.end(); ++i) {
      i->a = mol.GetAtom(i->ia);
      i->b = mol.GetAtom(i->ib);
      i->c = mol.GetAtom(i->ic);
      i->d = mol.GetAtom(i->id);
    }
  }
  
  vector3 OBFFConstraints::GetGradient(int a)
  {
    vector<OBFFConstraint>::iterator i;

    vector3 grad(0.0, 0.0, 0.0);
    
    for (i = _constraints.begin(); i != _constraints.end(); ++i)
      if (((*i).ia == a) || ((*i).ib == a)) {
        grad += i->GetGradient(a);
      }

    return grad;
  }
  
  void OBFFConstraints::Clear()
  {
    _constraints.clear();
  }
  
  int OBFFConstraints::Size() const
  {
    return _constraints.size();
  }
  
  double OBFFConstraints::GetConstraintEnergy()
  {
    vector<OBFFConstraint>::iterator i;
    double constraint_energy = 0.0;
        
    for (i = _constraints.begin(); i != _constraints.end(); ++i)
      if ( (i->type == OBFF_CONST_DISTANCE) || (i->type == OBFF_CONST_ANGLE) || (i->type == OBFF_CONST_TORSION) ) {
	vector3 da, db, dc, dd;
        double delta, delta2, rab, theta, phi, sine, dE;
 
	
	switch (i->type) {
          case OBFF_CONST_DISTANCE:
	    if ((i->a == NULL) || ((i->b) == NULL))
	      break;

	    da = (i->a)->GetVector();
            db = (i->b)->GetVector();
	
 	    rab = OBForceField::VectorLengthDerivative(da, db);
            delta = rab - i->constraint_value;
            delta2 = delta * delta;
            constraint_energy += i->factor * delta2;
            dE = 2.0 * i->factor * delta;

	    i->grada = dE * da;
            i->gradb = dE * db; 
            break;
          case OBFF_CONST_ANGLE:
	    if ((i->a == NULL) || ((i->b) == NULL) || ((i->c) == NULL))
	      break;

	    da = (i->a)->GetVector();
            db = (i->b)->GetVector();
            dc = (i->c)->GetVector();
	
 	    theta = OBForceField::VectorAngleDerivative(da, db, dc);
            delta = theta - i->constraint_value;
            delta2 = delta * delta;
            constraint_energy += 0.0002 * i->factor * delta2;
            dE = 0.0004 * i->factor * delta;

	    i->grada = dE * da;
            i->gradb = dE * db;
            i->gradc = dE * dc;
            break;
          case OBFF_CONST_TORSION:
            if ((i->a == NULL) || ((i->b) == NULL) || ((i->c) == NULL) || ((i->d) == NULL))
	      break;

	    /*
	    da = (i->a)->GetVector();
            db = (i->b)->GetVector();
            dc = (i->c)->GetVector();
            dd = (i->d)->GetVector();
	
 	    theta = OBForceField::VectorTorsionDerivative(da, db, dc, dd);
            if (IsNan(theta))
              theta = 1.0e-7;
            
	    if (theta >= 0)
              delta = theta - i->constraint_value;
	    else
              delta = theta + i->constraint_value;
            delta2 = delta * delta;
            constraint_energy += 0.00025 * delta2;
            dE = 0.000125 * delta;

	    if (theta >= 0) {
	      i->grada = dE * da;
              i->gradb = dE * db;
              i->gradc = dE * dc;
              i->gradd = dE * dd;
	    } else {
	      i->grada = -dE * da;
              i->gradb = -dE * db;
              i->gradc = -dE * dc;
              i->gradd = -dE * dd;
	    }
            */
	    break;
        
	  default:
	    break;
        }
      }
    return constraint_energy;
  }
  
  void OBFFConstraints::DeleteConstraint(int index)
  {
    vector<OBFFConstraint> constraints;
    vector<OBFFConstraint>::iterator i;
    int n = 0;

    for (i = _constraints.begin(); i != _constraints.end(); ++n, ++i) {
      if (n == index) {
        _constraints.erase(i);
        break;
      }
    }
  }
  
  void OBFFConstraints::SetFactor(double factor)
  {
    _factor = factor;
  }

  double OBFFConstraints::GetFactor()
  {
    return _factor;
  }

  void OBFFConstraints::AddIgnore(int a)
  {
    OBFFConstraint constraint;
    constraint.type = OBFF_CONST_IGNORE; // constraint type
    constraint.ia   = a; // atom to fix
    _constraints.push_back(constraint);
  }
  
  void OBFFConstraints::AddAtomConstraint(int a)
  {
    OBFFConstraint constraint;
    constraint.type = OBFF_CONST_ATOM; // constraint type
    constraint.ia   = a; // atom to fix
    constraint.factor = _factor;
    _constraints.push_back(constraint);
  }
  
  void OBFFConstraints::AddAtomXConstraint(int a)
  {
    OBFFConstraint constraint;
    constraint.type = OBFF_CONST_ATOM_X; // constraint type
    constraint.ia   = a; // atom to fix
    constraint.factor = _factor;
    _constraints.push_back(constraint);
  }
  
  void OBFFConstraints::AddAtomYConstraint(int a)
  {
    OBFFConstraint constraint;
    constraint.type = OBFF_CONST_ATOM_Y; // constraint type
    constraint.ia   = a; // atom to fix
    constraint.factor = _factor;
    _constraints.push_back(constraint);
  } 
  
  void OBFFConstraints::AddAtomZConstraint(int a)
  {
    OBFFConstraint constraint;
    constraint.type = OBFF_CONST_ATOM_Z; // constraint type
    constraint.ia   = a; // atom to fix
    constraint.factor = _factor;
    _constraints.push_back(constraint);
  }
  
  void OBFFConstraints::AddDistanceConstraint(int a, int b, double length)
  {
    OBFFConstraint constraint;
    constraint.type = OBFF_CONST_DISTANCE; // constraint type
    constraint.ia   = a; // atom a
    constraint.ib   = b; // atom b
    constraint.constraint_value = length; // bond length
    constraint.factor = _factor;
    _constraints.push_back(constraint);
  }
  
  void OBFFConstraints::AddAngleConstraint(int a, int b, int c, double angle)
  {
    OBFFConstraint constraint;
    constraint.type = OBFF_CONST_ANGLE; // constraint type
    constraint.ia   = a; // atom a
    constraint.ib   = b; // atom b
    constraint.ic   = c; // atom c
    constraint.constraint_value = angle; // angle
    constraint.factor = _factor;
    _constraints.push_back(constraint);
  }
  
  void OBFFConstraints::AddTorsionConstraint(int a, int b, int c, int d, double torsion)
  {
    OBFFConstraint constraint;
    constraint.type = OBFF_CONST_TORSION; // constraint type
    constraint.ia   = a; // atom a
    constraint.ib   = b; // atom b
    constraint.ic   = c; // atom c
    constraint.id   = d; // atom d
    constraint.constraint_value = torsion; // torsion
    constraint.factor = _factor;
    _constraints.push_back(constraint);
  }
  
  int OBFFConstraints::GetConstraintType(int index) const
  {
    if (index > _constraints.size())
      return 0;

    return _constraints[index].type; 
  }
  
  double OBFFConstraints::GetConstraintValue(int index) const 
  {
    if (index > _constraints.size())
      return 0;

    return _constraints[index].constraint_value; 
  }
  
  int OBFFConstraints::GetConstraintAtomA(int index) const
  {
    if (index > _constraints.size())
      return 0;

    return _constraints[index].ia; 
  }
  
  int OBFFConstraints::GetConstraintAtomB(int index) const
  {
    if (index > _constraints.size())
      return 0;

    return _constraints[index].ib; 
  }
  
  int OBFFConstraints::GetConstraintAtomC(int index) const
  {
    if (index > _constraints.size())
      return 0;

    return _constraints[index].ic; 
  }
  
  int OBFFConstraints::GetConstraintAtomD(int index) const
  {
    if (index > _constraints.size())
      return 0;

    return _constraints[index].id; 
  }
  
  bool OBFFConstraints::IsIgnored(int index)
  {
    vector<OBFFConstraint>::iterator i;
        
    for (i = _constraints.begin(); i != _constraints.end(); ++i)
      if (i->type == OBFF_CONST_IGNORE)
        if (i->ia == index)
	  return true;
    
    return false;
  }
  
  bool OBFFConstraints::IsFixed(int index)
  {
    vector<OBFFConstraint>::iterator i;
        
    for (i = _constraints.begin(); i != _constraints.end(); ++i)
      if (i->type == OBFF_CONST_ATOM)
        if (i->ia == index)
	  return true;
    
    return false;
  }
  
  bool OBFFConstraints::IsXFixed(int index)
  {
    vector<OBFFConstraint>::iterator i;
        
    for (i = _constraints.begin(); i != _constraints.end(); ++i)
      if (i->type == OBFF_CONST_ATOM_X)
        if ((i->a)->GetIdx() == index)
	  return true;
    
    return false;
  }
  
  bool OBFFConstraints::IsYFixed(int index)
  {
    vector<OBFFConstraint>::iterator i;
        
    for (i = _constraints.begin(); i != _constraints.end(); ++i)
      if (i->type == OBFF_CONST_ATOM_Y)
        if (i->ia == index)
	  return true;
    
    return false;
  }
  
  bool OBFFConstraints::IsZFixed(int index)
  {
     vector<OBFFConstraint>::iterator i;
        
    for (i = _constraints.begin(); i != _constraints.end(); ++i)
      if (i->type == OBFF_CONST_ATOM_Z)
        if (i->ia == index)
	  return true;
    
    return false;
  }
  
  //////////////////////////////////////////////////////////////////////////////////
  //
  // Methods for logging
  //
  //////////////////////////////////////////////////////////////////////////////////
  
  void OBForceField::PrintTypes()
  {
    IF_OBFF_LOGLVL_LOW {
      OBFFLog("\nA T O M   T Y P E S\n\n");
      OBFFLog("IDX\tTYPE\n");
      
      FOR_ATOMS_OF_MOL (a, _mol) {
        sprintf(_logbuf, "%d\t%s\n", a->GetIdx(), a->GetType());
        OBFFLog(_logbuf);
      }
    }
  }
    
  void OBForceField::PrintFormalCharges()
  {
    IF_OBFF_LOGLVL_LOW {
      OBFFLog("\nF O R M A L   C H A R G E S\n\n");
      OBFFLog("IDX\tCHARGE\n");
      
      FOR_ATOMS_OF_MOL (a, _mol) {
        sprintf(_logbuf, "%d\t%f\n", a->GetIdx(), a->GetPartialCharge());
        OBFFLog(_logbuf);
      }
    }
  }

  void OBForceField::PrintPartialCharges()
  {
    IF_OBFF_LOGLVL_LOW {
      OBFFLog("\nP A R T I A L   C H A R G E S\n\n");
      OBFFLog("IDX\tCHARGE\n");
      
      FOR_ATOMS_OF_MOL (a, _mol) {
        sprintf(_logbuf, "%d\t%f\n", a->GetIdx(), a->GetPartialCharge());
        OBFFLog(_logbuf);
      }
    }
  }
  
  void OBForceField::PrintVelocities()
  {
    IF_OBFF_LOGLVL_LOW {
      OBFFLog("\nA T O M   V E L O C I T I E S\n\n");
      OBFFLog("IDX\tVELOCITY\n");
      
      FOR_ATOMS_OF_MOL (a, _mol) {
        sprintf(_logbuf, "%d\t<%8.3f, %8.3f, %8.3f>\n", a->GetIdx(), _velocityPtr[a->GetIdx()], _velocityPtr[a->GetIdx()+1], _velocityPtr[a->GetIdx()+2]);
        OBFFLog(_logbuf);
      }
    }
  }
 
  bool OBForceField::SetLogFile(ostream* pos)
  {
    if(pos)
      _logos = pos;
    else
      _logos = &cout;
    
    return true;
  }
  
  //////////////////////////////////////////////////////////////////////////////////
  //
  // General
  //
  //////////////////////////////////////////////////////////////////////////////////
 
  bool OBForceField::Setup(OBMol &mol)
  {
    if (!_init) {
      ParseParamFile();
      _init = true;
      _velocityPtr = NULL;       
    }    
   
    if (IsSetupNeeded(mol) || !_validSetup) {
      if (_velocityPtr)
        delete [] _velocityPtr;
      _velocityPtr = NULL;       
      
      _mol = mol;
      if (_mol.NumAtoms() && _constraints.Size())
        _constraints.Setup(_mol);

      _mol.UnsetSSSRPerceived();
      
      if (!SetTypes()) {
        _validSetup = false;
        return false;
      }

      SetFormalCharges();
      SetPartialCharges();

      if (!SetupCalculations()) {
        _validSetup = false;
        return false;
      }
    }

    _validSetup = true;
    return true;
  }
  
  bool OBForceField::Setup(OBMol &mol, OBFFConstraints &constraints)
  {
    if (!_init) {
      ParseParamFile();
      _init = true;
      _velocityPtr = NULL;       
    }    
    
    if (IsSetupNeeded(mol)) {
      if (_velocityPtr)
        delete [] _velocityPtr;
      _velocityPtr = NULL;       
 
      _mol = mol;
      _constraints = constraints;
      if (_mol.NumAtoms() && _constraints.Size())
        _constraints.Setup(_mol);

      if (!SetTypes())
        return false;

      SetFormalCharges();
      SetPartialCharges();

      if (!SetupCalculations())
        return false;
    }
   
    return true;
  }
 
  bool OBForceField::SetLogLevel(int level)
  {
    _loglvl = level;

    return true;
  }
  
  // This function might need expanding, bu could really increase performance for
  // avogadro's AutoOpt tool
  bool OBForceField::IsSetupNeeded(OBMol &mol)
  {
    if (_mol.NumAtoms() != mol.NumAtoms())
      return true;
      
    if (_mol.NumBonds() != mol.NumBonds())
      return true;

    FOR_ATOMS_OF_MOL (atom, _mol) {
      if (atom->GetAtomicNum() != (mol.GetAtom(atom->GetIdx()))->GetAtomicNum())
        return true;
      if (atom->GetValence() != (mol.GetAtom(atom->GetIdx()))->GetValence())
        return true;
      if (atom->BOSum() != (mol.GetAtom(atom->GetIdx()))->BOSum())
        return true;
    }

    return false;
  }

  bool OBForceField::GetCoordinates(OBMol &mol)
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
  
  bool OBForceField::GetConformers(OBMol &mol)
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
      mol.SetEnergies(_energies);
      mol.SetConformer(_current_conformer);
    }
    
    return true;
  }
  
  bool OBForceField::SetCoordinates(OBMol &mol)
  { 
    OBAtom *atom;
    
    if (_mol.NumAtoms() != mol.NumAtoms())
      return false;
    
    // Copy coordinates for current conformer only
    FOR_ATOMS_OF_MOL (a, mol) {
      atom = _mol.GetAtom(a->GetIdx());
      atom->SetVector(a->GetVector());
    }

    return true;
  }
 
  bool OBForceField::SetConformers(OBMol &mol)
  { 
    OBAtom *atom;

    if (_mol.NumAtoms() != mol.NumAtoms())
      return false;
    
    FOR_ATOMS_OF_MOL (a, mol) {
      atom = _mol.GetAtom(a->GetIdx());
      atom->SetVector(a->GetVector());
    }

    //Copy conformer information
    if (mol.NumConformers() > 1) {
      int k,l;
      vector<double*> conf;
      double* xyz = NULL;
      for (k=0 ; k<mol.NumConformers() ; ++k) {
        xyz = new double [3*mol.NumAtoms()];
        for (l=0 ; l<(int) (3*mol.NumAtoms()) ; ++l)
          xyz[l] = mol.GetConformer(k)[l];
        conf.push_back(xyz);
      }
      _mol.SetConformers(conf);
      _mol.SetConformer(_current_conformer);
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
 
  //////////////////////////////////////////////////////////////////////////////////
  //
  // Structure generation
  //
  //////////////////////////////////////////////////////////////////////////////////

  int OBForceField::SystematicRotorSearchInitialize(unsigned int geomSteps)
  {
    OBRotorList rl;
    OBRotamerList rotamers;
    _origLogLevel = _loglvl;

    rl.Setup(_mol);
    rotamers.SetBaseCoordinateSets(_mol);
    rotamers.Setup(_mol, rl);
    
    IF_OBFF_LOGLVL_LOW {
      OBFFLog("\nS Y S T E M A T I C   R O T O R   S E A R C H\n\n");
      sprintf(_logbuf, "  NUMBER OF ROTATABLE BONDS: %d\n", rl.Size());
      OBFFLog(_logbuf);
    }

    if (!rl.Size()) { // only one conformer
      IF_OBFF_LOGLVL_LOW
        OBFFLog("  GENERATED ONLY ONE CONFORMER\n\n");
 
      ConjugateGradients(geomSteps); // final energy minimizatin for best conformation
      
      return 1; // there are no more conformers
    }

    OBRotorKeys rotorKeys;
    OBRotorIterator ri;
    OBRotor *rotor = rl.BeginRotor(ri);
    for (int i = 1; i < rl.Size() + 1; ++i, rotor = rl.NextRotor(ri)) // foreach rotor
      rotorKeys.AddRotor(rotor->GetResolution().size());
    
    while (rotorKeys.Next()) 
      rotamers.AddRotamer(rotorKeys.GetKey());

    rotamers.ExpandConformerList(_mol, _mol.GetConformers());
      
    IF_OBFF_LOGLVL_LOW {
      sprintf(_logbuf, "  GENERATED %d CONFORMERS\n\n", _mol.NumConformers());
      OBFFLog(_logbuf);
      OBFFLog("CONFORMER     ENERGY\n");
      OBFFLog("--------------------\n");
    }

    _current_conformer = 0;
    _energies.clear();
      
    return _mol.NumConformers();
  }
  
  bool OBForceField::SystematicRotorSearchNextConformer(unsigned int geomSteps) 
  {
    char _logbuf[100];
    
    if (_current_conformer >=  _mol.NumConformers()) { // done
      // Select conformer with lowest energy
      int best_conformer = 0;
      for (int i = 0; i < _mol.NumConformers(); i++) {
        if (_energies[i] < _energies[best_conformer])
          best_conformer = i;
      }
  
      IF_OBFF_LOGLVL_LOW {
        sprintf(_logbuf, "\n  CONFORMER %d HAS THE LOWEST ENERGY\n\n",  best_conformer + 1);
        OBFFLog(_logbuf);
      }

      _mol.SetConformer(best_conformer);
      _current_conformer = best_conformer;
    
      return false;
    }

    _mol.SetConformer(_current_conformer); // select conformer
    
    _loglvl = OBFF_LOGLVL_NONE;
    ConjugateGradients(geomSteps); // energy minimization for conformer
    _loglvl = _origLogLevel;
    
    _energies.push_back(Energy(false)); // calculate and store energy
      
    IF_OBFF_LOGLVL_LOW {
      sprintf(_logbuf, "   %3d   %20.3f\n", (_current_conformer + 1), _energies[_current_conformer]);
      OBFFLog(_logbuf);
    }
    
    _current_conformer++;
    return true;
  }

  void OBForceField::SystematicRotorSearch(unsigned int geomSteps) 
  {
    if (SystematicRotorSearchInitialize(geomSteps))
      while (SystematicRotorSearchNextConformer(geomSteps)) {}
  }

  void OBForceField::RandomRotorSearchInitialize(unsigned int conformers, unsigned int geomSteps) 
  {
    OBRotorList rl;
    OBRotamerList rotamers;
    OBRotorIterator ri;
    OBRotor *rotor;

    OBRandom generator;
    generator.TimeSeed();
    _origLogLevel = _loglvl;

    if (_mol.GetCoordinates() == NULL)
      return;

    rl.Setup(_mol);
    rotamers.SetBaseCoordinateSets(_mol);
    rotamers.Setup(_mol, rl);
    
    IF_OBFF_LOGLVL_LOW {
      OBFFLog("\nR A N D O M   R O T O R   S E A R C H\n\n");
      sprintf(_logbuf, "  NUMBER OF ROTATABLE BONDS: %d\n", rl.Size());
      OBFFLog(_logbuf);

      int combinations = 1;
      for (rotor = rl.BeginRotor(ri); rotor;
           rotor = rl.NextRotor(ri)) {
        combinations *= rotor->GetResolution().size();
      }
      sprintf(_logbuf, "  NUMBER OF POSSIBLE ROTAMERS: %d\n", combinations);
      OBFFLog(_logbuf);
    }

    if (!rl.Size()) { // only one conformer
      IF_OBFF_LOGLVL_LOW
        OBFFLog("  GENERATED ONLY ONE CONFORMER\n\n");
 
      _loglvl = OBFF_LOGLVL_NONE;
      ConjugateGradients(geomSteps); // energy minimization for conformer
      _loglvl = _origLogLevel;

      return;
    }

    std::vector<int> rotorKey(rl.Size() + 1, 0); // indexed from 1

    for (int c = 0; c < conformers; ++c) {
      rotor = rl.BeginRotor(ri);
      for (int i = 1; i < rl.Size() + 1; ++i, rotor = rl.NextRotor(ri)) {
        // foreach rotor
        rotorKey[i] = generator.NextInt() % rotor->GetResolution().size();
      }
      rotamers.AddRotamer(rotorKey);
    }

    rotamers.ExpandConformerList(_mol, _mol.GetConformers());
      
    IF_OBFF_LOGLVL_LOW {
      sprintf(_logbuf, "  GENERATED %d CONFORMERS\n\n", _mol.NumConformers());
      OBFFLog(_logbuf);
      OBFFLog("CONFORMER     ENERGY\n");
      OBFFLog("--------------------\n");
    }
  
    _current_conformer = 0;
    _energies.clear();
  }

  bool OBForceField::RandomRotorSearchNextConformer(unsigned int geomSteps) 
  {
    char _logbuf[100];
    
    if (_current_conformer >=  _mol.NumConformers()) { // done
      // Select conformer with lowest energy
      int best_conformer = 0;
      for (int i = 0; i < _mol.NumConformers(); i++) {
        if (_energies[i] < _energies[best_conformer])
          best_conformer = i;
      }
  
      IF_OBFF_LOGLVL_LOW {
        sprintf(_logbuf, "\n  CONFORMER %d HAS THE LOWEST ENERGY\n\n",  best_conformer + 1);
        OBFFLog(_logbuf);
      }

      _mol.SetConformer(best_conformer);
      _current_conformer = best_conformer;
      
      return false;
    }
    
    _mol.SetConformer(_current_conformer); // select conformer

    _loglvl = OBFF_LOGLVL_NONE;
    ConjugateGradients(geomSteps); // energy minimization for conformer
    _loglvl = _origLogLevel;
    
    _energies.push_back(Energy(false)); // calculate and store energy
      
    IF_OBFF_LOGLVL_LOW {
      sprintf(_logbuf, "   %3d      %8.3f\n", (_current_conformer + 1), _energies[_current_conformer]);
      OBFFLog(_logbuf);
    }
    
    _current_conformer++;
    return true;
  } 
    
  void OBForceField::RandomRotorSearch(unsigned int conformers, unsigned int geomSteps) 
  {
    RandomRotorSearchInitialize(conformers, geomSteps);
    while (RandomRotorSearchNextConformer(geomSteps)) {}
  }

  void Reweight(std::vector< std::vector <double> > &rotorWeights, std::vector<int> rotorKey, double bonus)
  {
    double fraction, minWeight, maxWeight;
    bool improve = (bonus > 0.0);

    for (int i = 1; i < rotorWeights.size() - 1; ++i) {
      if (improve && rotorWeights[i][rotorKey[i]] > 0.999 - bonus)
        continue;
      if (!improve && rotorWeights[i][rotorKey[i]] < 0.001 - bonus) // bonus < 0
        continue;
      
      // Check to make sure we don't kill some poor weight
      minWeight = maxWeight = rotorWeights[i][0];
      for (int j = 1; j < rotorWeights[i].size(); ++j) {
        if (j == rotorKey[i])
          continue; // we already checked for problems with this entry
        if (rotorWeights[i][j] < minWeight)
          minWeight = rotorWeights[i][j];
        if (rotorWeights[i][j] > maxWeight)
          maxWeight = rotorWeights[i][j];
      }

      fraction = bonus / (rotorWeights[i].size() - 1);
      if (improve && fraction > minWeight) {
        bonus = minWeight / 2.0;
        fraction = bonus / (rotorWeights[i].size() - 1);
      }
      if (!improve && fraction > maxWeight) {
        bonus = (maxWeight - 1.0) / 2.0; // negative "bonus"
        fraction = bonus / (rotorWeights[i].size() - 1);
      }

      for (int j = 0; j < rotorWeights[i].size(); ++j) {
        if (j == rotorKey[i])
          rotorWeights[i][j] += bonus;
        else
          rotorWeights[i][j] -= fraction;
      }
    }
  }
                
  
  void OBForceField::WeightedRotorSearch(unsigned int conformers, unsigned int geomSteps)
  {
    OBRotorList rl;
    OBRotamerList rotamers;
    OBRotorIterator ri;
    OBRotor *rotor;

    OBRandom generator;
    generator.TimeSeed();
    int origLogLevel = _loglvl;

    if (_mol.GetCoordinates() == NULL)
      return;

    rl.Setup(_mol);
    rotamers.SetBaseCoordinateSets(_mol);
    rotamers.Setup(_mol, rl);
    
    IF_OBFF_LOGLVL_LOW {
      OBFFLog("\nW E I G H T E D   R O T O R   S E A R C H\n\n");
      sprintf(_logbuf, "  NUMBER OF ROTATABLE BONDS: %d\n", rl.Size());
      OBFFLog(_logbuf);

      int combinations = 1;
      for (rotor = rl.BeginRotor(ri); rotor;
           rotor = rl.NextRotor(ri)) {
        combinations *= rotor->GetResolution().size();
      }
      sprintf(_logbuf, "  NUMBER OF POSSIBLE ROTAMERS: %d\n", combinations);
      OBFFLog(_logbuf);
    }

    if (!rl.Size()) { // only one conformer
      IF_OBFF_LOGLVL_LOW
        OBFFLog("  GENERATED ONLY ONE CONFORMER\n\n");
 
      _loglvl = OBFF_LOGLVL_NONE;
      ConjugateGradients(geomSteps); // energy minimization for conformer
      _loglvl = origLogLevel;

      return;
    }

    double *initialCoord = new double [_mol.NumAtoms() * 3]; // initial state
    memcpy((char*)initialCoord,(char*)_mol.GetCoordinates(),sizeof(double)*3*_mol.NumAtoms());

    // key for generating particular conformers
    std::vector<int> rotorKey(rl.Size() + 1, 0); // indexed from 1

    // initialize the weights
    std::vector< std::vector <double> > rotorWeights;
    std::vector< double > weightSet;
    rotorWeights.push_back(weightSet); // empty set for unused index 0

    double weight;
    double bestE, worstE, currentE;
    std::vector<double> energies;
    // First off, we test out each rotor position. How good (or bad) is it?
    // This lets us pre-weight the search to a useful level
    // So each rotor is considered in isolation
    rotor = rl.BeginRotor(ri);
    for (int i = 1; i < rl.Size() + 1; ++i, rotor = rl.NextRotor(ri)) {
      rotorKey[i] = 0; // for lack of another default. Unfortunately, this isn't "no rotation"
    }
    
    rotor = rl.BeginRotor(ri);
    int confCount = 0;
    int settings;
    for (int i = 1; i < rl.Size() + 1; ++i, rotor = rl.NextRotor(ri)) {
      // foreach rotor
      energies.clear();
      for (unsigned int j = 0; j < rotor->GetResolution().size(); j++) {
        // foreach rotor position
        _mol.SetCoordinates(initialCoord);
        rotorKey[i] = j;
        rotamers.SetCurrentCoordinates(_mol, rotorKey);
        
        currentE = Energy(false);

        if (j == 0) 
          bestE = worstE = currentE;
        else {
          if (currentE > worstE)
            worstE = currentE;
          else if (currentE < bestE)
            bestE = currentE;
        }
        energies.push_back(currentE);
      }
      rotorKey[i] = 0; // back to the previous setting before we go to another rotor
      
      weightSet.clear();
      // first loop through and calculate the relative populations from Boltzmann
      double totalPop = 0.0;
      for (unsigned int j = 0; j < rotor->GetResolution().size(); j++) {
        currentE = energies[j];
        // add back the Boltzmann population for these relative energies at 300K (assuming kJ/mol)
        energies[j] = exp(-1.0*fabs(currentE - bestE) / 2.5);
        totalPop += energies[j];
      }
      // now set the weights
      for (unsigned int j = 0; j < rotor->GetResolution().size(); j++) {
        if (IsNear(worstE, bestE, 1.0e-3))
          weight = 1 / rotor->GetResolution().size();
        else
          weight = energies[j]/totalPop;
        weightSet.push_back(weight);
      }
      rotorWeights.push_back(weightSet);
    }

    int best_conformer;
    double penalty; // for poor performance
    double randFloat; // generated random number -- used to pick a rotor
    double total; // used to calculate the total probability 
    double *bestCoordPtr = new double [_mol.NumAtoms() * 3]; // coordinates for best conformer
    
    // Start with the current coordinates
    bestE = worstE = Energy(false);
    // We're later going to add this back to the molecule as a new conformer
    memcpy((char*)bestCoordPtr,(char*)_mol.GetCoordinates(),sizeof(double)*3*_mol.NumAtoms());    
    
    // Now we actually test some weightings
    IF_OBFF_LOGLVL_LOW {
      sprintf(_logbuf, "  GENERATED %d CONFORMERS\n\n", conformers);
      OBFFLog(_logbuf);
      OBFFLog("CONFORMER     ENERGY\n");
      OBFFLog("--------------------\n");
    }

    for (int c = 0; c < conformers; ++c) {
      _mol.SetCoordinates(initialCoord);

      // Choose the rotor key based on current weightings
      rotor = rl.BeginRotor(ri);
      for (int i = 1; i < rl.Size() + 1; ++i, rotor = rl.NextRotor(ri)) {
        // foreach rotor
        randFloat = generator.NextFloat();
        total = 0.0;
        rotorKey[i] = 0;
        for (unsigned int j = 0; j < rotor->GetResolution().size(); j++) {
          if (randFloat > total && randFloat < (total+ rotorWeights[i][j])) {
            rotorKey[i] = j;
            break;      
          }
          else
            total += rotorWeights[i][j];
        }
      }
      rotamers.SetCurrentCoordinates(_mol, rotorKey);

      //      _loglvl = OBFF_LOGLVL_NONE;
      //      SteepestDescent(geomSteps); // energy minimization for conformer
      //      _loglvl = origLogLevel;
      currentE = Energy(false);

      IF_OBFF_LOGLVL_LOW {
        sprintf(_logbuf, "   %3d      %8.3f\n", c + 1, currentE);
        OBFFLog(_logbuf);
      }
      
      if (!isfinite(currentE))
        continue;

      if (currentE < bestE) {
        bestE = currentE;
        memcpy((char*)bestCoordPtr,(char*)_mol.GetCoordinates(),sizeof(double)*3*_mol.NumAtoms());
        best_conformer = c;

        // improve this rotorKey
        Reweight(rotorWeights, rotorKey, +0.11);
      } else if (currentE > worstE) { // horrible!
        worstE = currentE;
        
        // penalize this rotorKey
        Reweight(rotorWeights, rotorKey, -0.11);
      } else {
        double slope = -0.2 / (worstE - bestE);
        Reweight(rotorWeights, rotorKey, (currentE - bestE)*slope);
      }
    }

    IF_OBFF_LOGLVL_LOW {
      sprintf(_logbuf, "\n  LOWEST ENERGY: %8.3f\n\n",
              bestE);
      OBFFLog(_logbuf);
    }

    // debugging output to see final weightings of each rotor setting
    IF_OBFF_LOGLVL_HIGH {
      OBFFLog("Final Weights: \n");
      for (int i = 1; i < rotorWeights.size() - 1; ++i) {
        sprintf(_logbuf, " Weight: %d", i);
        OBFFLog(_logbuf);
        for (int j = 0; j < rotorWeights[i].size(); ++j) {
          sprintf(_logbuf, " %8.3f", rotorWeights[i][j]);
          OBFFLog(_logbuf);
        }
        OBFFLog("\n");
     }
   }

    _mol.AddConformer(bestCoordPtr);
    _current_conformer = _mol.NumConformers() - 1;
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
          sprintf(_logbuf, " %8.4f ", matrix[i][j]);
          OBFFLog(_logbuf);
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
          sprintf(_logbuf, " %8.4f ", matrix[i][j]);
          OBFFLog(_logbuf);
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
      char _logbuf[100];

      OBFFLog("RANDOM DISTANCE MATRIX BETWEEN LIMITS\n\n");
      for (i=0; i<N; i++) {
        OBFFLog("[");
        for (j=0; j<N; j++) {
          sprintf(_logbuf, " %8.4f ", matrix[i][j]);
          OBFFLog(_logbuf);
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
      char _logbuf[100];

      OBFFLog("METRIC MATRIX\n\n");
      for (i=0; i<N; i++) {
        OBFFLog("[");
        for (j=0; j<N; j++) {
          sprintf(_logbuf, " %8.4f ", G[i][j]);
          OBFFLog(_logbuf);
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
        sprintf(_logbuf, "%8.4f ", eigenvalues[i]);
        OBFFLog(_logbuf);
      }
      OBFFLog("\n");

      OBFFLog("EIGENVECTORS OF METRIC MATRIX\n\n");
      for (i=0; i<N; i++) {
        OBFFLog("[");
        for (j=0; j<N; j++) {
          sprintf(_logbuf, " %8.4f ", eigenvectors[i][j]);
          OBFFLog(_logbuf);
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
  
  //////////////////////////////////////////////////////////////////////////////////
  //
  // Energy Minimization
  //
  //////////////////////////////////////////////////////////////////////////////////
  
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

    // See the following implementation

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

  double OBForceField::LineSearch(double *currentCoords, double *direction)
  {
    int numCoords = _mol.NumAtoms() * 3;
    double e_n1, e_n2, step, alpha, tempStep;
    double *lastStep = new double [numCoords];

    alpha = 0.0; // Scale factor along direction vector
    step = 0.2;
    double trustRadius = 0.3; // don't move further than 0.3 Angstroms
    
    //e_n1 = Energy(false) + ; // calculate e_k (before any move)
    e_n1 = Energy(false) + _constraints.GetConstraintEnergy();
    // dir = GetGradient(&*a) + _constraints.GetGradient(a->GetIdx());
    
    unsigned int i;
    for (i=0; i < 10; ++i) {
      // Save the current position, before we take a step
      memcpy((char*)lastStep,(char*)currentCoords,sizeof(double)*numCoords);
      
      // Vectorizing this would be a big benefit
      // Need to look up using BLAS or Eigen or whatever
      for (unsigned int c = 0; c < numCoords; ++c) {
        if (isfinite(direction[c])) { 
          // make sure we don't have NaN or infinity
          tempStep = direction[c] * step;

          if (tempStep > trustRadius) // positive big step
            currentCoords[c] += trustRadius;
          else if (tempStep < -1.0 * trustRadius) // negative big step
            currentCoords[c] -= trustRadius;
          else
            currentCoords[c] += direction[c] * step;
        }
      }
    
      //e_n2 = Energy(false); // calculate e_k+1
      e_n2 = Energy(false) + _constraints.GetConstraintEnergy();
      
      // convergence criteria: A higher precision here 
      // only takes longer with the same result.
      if (IsNear(e_n2, e_n1, 1.0e-3))
        break;

      if (e_n2 > e_n1) { // decrease stepsize
        step *= 0.1;
        // move back to the last step
        memcpy((char*)currentCoords,(char*)lastStep,sizeof(double)*numCoords);
      } else if (e_n2 < e_n1) {  // increase stepsize
        e_n1 = e_n2;
        alpha += step; // we've moved some distance
        step *= 2.15;
        if (step > 1.0)
          step = 1.0;
      }
      
    }
    //cout << "LineSearch steps: " << i << endl;

    delete [] lastStep;

    return alpha;
  }
  
  bool OBForceField::DetectExplosion()
  {
    FOR_ATOMS_OF_MOL (atom, _mol) {
      if (!isfinite(atom->GetX()))
        return true;
      if (!isfinite(atom->GetY()))
        return true;
      if (!isfinite(atom->GetZ()))
        return true;
    }
    
    FOR_BONDS_OF_MOL (bond, _mol) {
      if (bond->GetLength() > 30.0)
        return true;
    }
        
    return false;
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
    char _logbuf[100];
    
    atom->SetVector(9.0, 9.0, 0.0);
    e_n1 = atom->x() * atom->x() + 2 * (atom->y() * atom->y());
    
    IF_OBFF_LOGLVL_LOW {
      OBFFLog("\nV A L I D A T E   S T E E P E S T   D E S C E N T\n\n");
      sprintf(_logbuf, "STEPS = %d\n\n",  steps);
      OBFFLog(_logbuf);
      OBFFLog("STEP n     E(n)       E(n-1)    \n");
      OBFFLog("--------------------------------\n");
    }

    for (int i = 1; i <= steps; i++) {
      grad.Set(-2*atom->x(), -4*atom->y(), 0.0);
      grad = ValidateLineSearch(&*atom, grad);
      atom->SetVector(atom->x() + grad.x(), atom->y() + grad.y(), 0.0);
      e_n2 = atom->x() * atom->x() + 2 * (atom->y() * atom->y());
      
      IF_OBFF_LOGLVL_LOW {
        sprintf(_logbuf, " %4d    %8.3f    %8.3f\n", i, e_n2, e_n1);
        OBFFLog(_logbuf);
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
    char _logbuf[100];

    firststep = true;
    atom->SetVector(9.0, 9.0, 0.0);
    e_n1 = atom->x() * atom->x() + 2 * (atom->y() * atom->y());
 
    IF_OBFF_LOGLVL_LOW {
      OBFFLog("\nV A L I D A T E   C O N J U G A T E   G R A D I E N T S\n\n");
      sprintf(_logbuf, "STEPS = %d\n\n",  steps);
      OBFFLog(_logbuf);
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
          sprintf(_logbuf, " %4d    %8.3f    %8.3f\n", i, e_n2, e_n1);
          OBFFLog(_logbuf);
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
          sprintf(_logbuf, " %4d    %8.3f    %8.3f\n", i, e_n2, e_n1);
          OBFFLog(_logbuf);
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

    _e_n1 = Energy() + _constraints.GetConstraintEnergy();
    
    IF_OBFF_LOGLVL_LOW {
      OBFFLog("\nS T E E P E S T   D E S C E N T\n\n");
      sprintf(_logbuf, "STEPS = %d\n\n",  steps);
      OBFFLog(_logbuf);
      OBFFLog("STEP n     E(n)       E(n-1)    \n");
      OBFFLog("--------------------------------\n");
    }
  }
 
  bool OBForceField::SteepestDescentTakeNSteps(int n) 
  {
    int _ncoords = _mol.NumAtoms() * 3;
    int coordIdx;
    double e_n2, alpha;
    double *gradPtr = new double[_ncoords];
    memset(gradPtr, '\0', sizeof(double)*_ncoords);
    vector3 dir;

    for (int i = 1; i <= n; i++) {
      _cstep++;
      
      FOR_ATOMS_OF_MOL (a, _mol) {
        if (!_constraints.IsFixed(a->GetIdx())) {
          coordIdx = (a->GetIdx() - 1) * 3;
          if (_method & OBFF_ANALYTICAL_GRADIENT)
            dir = GetGradient(&*a) + _constraints.GetGradient(a->GetIdx());
          else
            dir = NumericalDerivative(&*a) + _constraints.GetGradient(a->GetIdx());
          
          if (!_constraints.IsXFixed(a->GetIdx()))
	    gradPtr[coordIdx] = dir.x();
          if (!_constraints.IsYFixed(a->GetIdx()))
            gradPtr[coordIdx+1] = dir.y();
          if (!_constraints.IsZFixed(a->GetIdx()))
            gradPtr[coordIdx+2] = dir.z();
        }
      }
      alpha = LineSearch(_mol.GetCoordinates(), gradPtr);
      e_n2 = Energy() + _constraints.GetConstraintEnergy();
      
      IF_OBFF_LOGLVL_LOW {
        if (_cstep % 10 == 0) {
          sprintf(_logbuf, " %4d    %8.3f    %8.3f\n", _cstep, e_n2, _e_n1);
          OBFFLog(_logbuf);
        }
      }

      if (IsNear(e_n2, _e_n1, _econv)) {
        IF_OBFF_LOGLVL_LOW
          OBFFLog("    STEEPEST DESCENT HAS CONVERGED\n");
        delete [] gradPtr;
        return false;
      }
      
      if (_nsteps == _cstep) {
        delete [] gradPtr;
        return false;
      }

      _e_n1 = e_n2;
    }

    delete [] gradPtr;
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
    double e_n2, alpha;
    vector3 grad2, dir2;

    //_cstep = 0;
    _cstep = 0;
    _nsteps = steps;
    _econv = econv;
    _method = method;
    _ncoords = _mol.NumAtoms() * 3;

    _e_n1 = Energy() + _constraints.GetConstraintEnergy();
    
    IF_OBFF_LOGLVL_LOW {
      OBFFLog("\nC O N J U G A T E   G R A D I E N T S\n\n");
      sprintf(_logbuf, "STEPS = %d\n\n",  steps);
      OBFFLog(_logbuf);
      OBFFLog("STEP n     E(n)       E(n-1)    \n");
      OBFFLog("--------------------------------\n");
    }

    if (_grad1 != NULL)
      delete [] _grad1;
    _grad1 = new double[_ncoords];
    memset(_grad1, '\0', sizeof(double)*_ncoords);
    if (_dir1 != NULL)
      delete [] _dir1;
    _dir1 = new double[_ncoords];
    memset(_grad1, '\0', sizeof(double)*_ncoords);

    // Take the first step (same as steepest descent because there is no 
    // gradient from the previous step.
    FOR_ATOMS_OF_MOL (a, _mol) {
      if (!_constraints.IsFixed(a->GetIdx())) 
      {
        if (_method & OBFF_ANALYTICAL_GRADIENT)
          grad2 = GetGradient(&*a) + _constraints.GetGradient(a->GetIdx());
        else
          grad2 = NumericalDerivative(&*a) + _constraints.GetGradient(a->GetIdx());
        dir2 = grad2;

        int coordIdx = (a->GetIdx() - 1) * 3;
        if (!_constraints.IsXFixed(a->GetIdx())) {
          _grad1[coordIdx] = grad2.x();
          _dir1[coordIdx] = grad2.x();
        }
	if (!_constraints.IsYFixed(a->GetIdx())) {
          _grad1[coordIdx + 1] = grad2.y();
          _dir1[coordIdx + 1] = grad2.y();
        }
	if (!_constraints.IsZFixed(a->GetIdx())) {
	  _grad1[coordIdx + 2] = grad2.z();
          _dir1[coordIdx + 2] = grad2.z();
	}
      }
    }
    alpha = LineSearch(_mol.GetCoordinates(), _dir1);
    e_n2 = Energy() + _constraints.GetConstraintEnergy();
      
    IF_OBFF_LOGLVL_LOW {
      sprintf(_logbuf, " %4d    %8.3f    %8.3f\n", _cstep, e_n2, _e_n1);
      OBFFLog(_logbuf);
    }
 
    _e_n1 = e_n2;
  }
  
  bool OBForceField::ConjugateGradientsTakeNSteps(int n)
  {
    double e_n2;
    double g2g2, g1g1, g2g1, alpha;
    vector3 grad2, dir2;
    vector3 grad1, dir1; // temporaries to perform dot product, etc.
    int coordIdx;
    
    int printStep = min(n, 10); // print results every 10 steps, or fewer

    if (_ncoords != _mol.NumAtoms() * 3)
      return false;

    e_n2 = 0.0;
    
    for (int i = 1; i <= n; i++) {
      _cstep++;
      
      FOR_ATOMS_OF_MOL (a, _mol) {
        if (!_constraints.IsFixed(a->GetIdx())) {
          coordIdx = (a->GetIdx() - 1) * 3;

          if (_method & OBFF_ANALYTICAL_GRADIENT)
            grad2 = GetGradient(&*a) + _constraints.GetGradient(a->GetIdx());
          else
            grad2 = NumericalDerivative(&*a) + _constraints.GetGradient(a->GetIdx());
          
          // Fletcher-Reeves formula for Beta
          // http://en.wikipedia.org/wiki/Nonlinear_conjugate_gradient_method
          // NOTE: We make sure to reset and use the steepest descent direction
          //   after NumAtoms steps
          if (_cstep % _mol.NumAtoms() != 0) {
            g2g2 = dot(grad2, grad2);
            grad1 = vector3(_grad1[coordIdx], _grad1[coordIdx+1], _grad1[coordIdx+2]);
            g1g1 = dot(grad1, grad1);
            g2g1 = g2g2 / g1g1;
            dir1 = vector3(_dir1[coordIdx], _dir1[coordIdx+1], _dir1[coordIdx+2]);
            dir2 = grad2 + g2g1 * dir1;
          } else { // reset conj. direction
            dir2 = grad2;
          }
 
          
          if (!_constraints.IsXFixed(a->GetIdx())) {
            _grad1[coordIdx] = grad2.x();
            _dir1[coordIdx] = dir2.x();
          }
          if (!_constraints.IsYFixed(a->GetIdx())) {
   	    _grad1[coordIdx + 1] = grad2.y();
            _dir1[coordIdx + 1] = dir2.y();
          }
          if (!_constraints.IsZFixed(a->GetIdx())) {
            _dir1[coordIdx + 2] = dir2.z();
            _grad1[coordIdx + 2] = grad2.z();
          }
        }
      }
      alpha = LineSearch(_mol.GetCoordinates(), _dir1);
      for (unsigned int j = 0; j < _ncoords; ++j)
        _dir1[j] *= alpha;

      e_n2 = Energy() + _constraints.GetConstraintEnergy();
	
      IF_OBFF_LOGLVL_LOW {
        if (_cstep % 10 == 0) {
          sprintf(_logbuf, " %4d    %8.3f    %8.3f\n", _cstep, e_n2, _e_n1);
          OBFFLog(_logbuf);
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
  
  //////////////////////////////////////////////////////////////////////////////////
  //
  // Molecular Dynamics
  //
  //////////////////////////////////////////////////////////////////////////////////
  // Most OpenBabel MD alogrithms are based on:
  // GROMACS, Groningen Machine for Chemical Simulations, USER MANUAL, version 3.3
  //
  // Quantity:		Symbol:		Unit:
  // length		r		A = 10e-10m
  // mass		m 		amu = 1.6605402 10e-27kg
  // time		t		ps = 10e-12s
  // temperature	T		K
  //
  // force		F		kcal mol^-1 A^-1
  // acceleration	a		kcal mol^-1 A^-1 1000x amu-1 = A ps^-2 (*)
  // velocity		v		A ps^-1
  //
  // The force we calculate comes from the force field. Currently this functions 
  // only works for MMFF94 and so the unit for energy is kcal, but this does not 
  // affect what is explained below. (note: A isn't a SI unit either)
  //
  //     [   kcal   ]          F    [     kcal    ]   [       kcal       ]   [      A      ]   [   A    ]
  // F = [ -------- ]     a = --- = [ ----------- ] = [ ---------------- ] = [ ----------- ] = [ ------ ]
  //     [  mol  A  ]          m    [  mol A amu  ]   [  mol A 10^-27kg  ]   [  10^-24s^2  ]   [  ps^2  ]
  //
  // This means that if we devide the force (comming from the FF) by 1000x its mass in amu, 
  // we get the acceleration in A ps^-2.
  
  // gromacs user manual page 17
  void OBForceField::GenerateVelocities() 
  {
    cout << "OBForceField::GenerateVelocities()" << endl;
    OBRandom generator;
    generator.TimeSeed();
    _ncoords = _mol.NumAtoms() * 3;
    int velocityIdx;
    double velocity, kB;
    kB = 0.00831451 / KCAL_TO_KJ; // kcal/(mol*K)
    
    _velocityPtr = new double[_ncoords];
    memset(_velocityPtr, '\0', sizeof(double)*_ncoords);
    
    FOR_ATOMS_OF_MOL (a, _mol) {
      if (!_constraints.IsFixed(a->GetIdx())) {
        velocityIdx = (a->GetIdx() - 1) * 3;
      
        // add twelve random numbers between 0.0 and 1.0,
        // subtract 6.0 from their sum, multiply with sqrt(kT/m)
	if (!_constraints.IsXFixed(a->GetIdx())) {
          velocity = 0.0;
          for (int i=0; i < 12; ++i) 
            velocity += generator.NextFloat();
          velocity -= 6.0;
          velocity *= sqrt((kB * _T)/ (1000 * a->GetAtomicMass()));
          _velocityPtr[velocityIdx] = velocity; // x10: gromacs uses nm instead of A
        }
        
	if (!_constraints.IsYFixed(a->GetIdx())) {
	  velocity = 0.0;
          for (int i=0; i < 12; ++i)
            velocity += generator.NextFloat();
          velocity -= 6.0;
          velocity *= sqrt((kB * _T)/ (1000 * a->GetAtomicMass()));
          _velocityPtr[velocityIdx+1] = velocity; // idem
        }

	if (!_constraints.IsZFixed(a->GetIdx())) {
	  velocity = 0.0;
          for (int i=0; i < 12; ++i)
	    velocity += generator.NextFloat();
          velocity -= 6.0;
          velocity *= sqrt((kB * _T)/ (1000 * a->GetAtomicMass()));
          _velocityPtr[velocityIdx+2] = velocity; // idem
        }
      }
    }

    CorrectVelocities();
  }

  void OBForceField::CorrectVelocities()
  {
    _ncoords = _mol.NumAtoms() * 3;
    int velocityIdx;
    double velocity, kB, E_kin, E_kin2, factor;
    kB = 0.00831451 / KCAL_TO_KJ; // kcal/(mol*K)
 
    // E_kin = 0.5 * Ndf * kB * T
    E_kin = _ncoords * kB * _T;
    //cout << "E_{kin} = Ndf * kB * T = " << E_kin << endl;
    
    // E_kin = 0.5 * sum( m_i * v_i^2 )
    E_kin2 = 0.0;
    FOR_ATOMS_OF_MOL (a, _mol) {
      velocityIdx = (a->GetIdx() - 1) * 3;
      
      velocity = sqrt( _velocityPtr[velocityIdx] * _velocityPtr[velocityIdx] + 
                       _velocityPtr[velocityIdx+1] * _velocityPtr[velocityIdx+1] + 
                       _velocityPtr[velocityIdx+2] * _velocityPtr[velocityIdx+2] );

      E_kin2 += 1000 * a->GetAtomicMass() * velocity * velocity;
    }
    //cout << "E_{kin} = sum( m_i * v_i^2 ) = " << E_kin2 << endl;

    // correct
    factor = sqrt(E_kin / E_kin2);
    FOR_ATOMS_OF_MOL (a, _mol) {
      velocityIdx = (a->GetIdx() - 1) * 3;
      _velocityPtr[velocityIdx] *= factor; 
      _velocityPtr[velocityIdx+1] *= factor; 
      _velocityPtr[velocityIdx+2] *= factor; 
    }
    
    // E_kin = 0.5 * sum( m_i * v_i^2 )
    E_kin2 = 0.0;
    FOR_ATOMS_OF_MOL (a, _mol) {
      velocityIdx = (a->GetIdx() - 1) * 3;
      
      velocity = sqrt( _velocityPtr[velocityIdx] * _velocityPtr[velocityIdx] + 
                       _velocityPtr[velocityIdx+1] * _velocityPtr[velocityIdx+1] + 
                       _velocityPtr[velocityIdx+2] * _velocityPtr[velocityIdx+2] );

      E_kin2 += 1000 * a->GetAtomicMass() * velocity * velocity;
    }
    //cout << "E_{kin_corr} = sum( m_i * v_i^2 ) = " << E_kin2 << endl;
  }
  
  void OBForceField::MolecularDynamicsTakeNSteps(int n, double T, double timestep, int method)
  {
    int coordIdx;
    double timestep2;
    vector3 force, pos, accel;
    double kB = 0.00831451 / KCAL_TO_KJ; // kcal/(mol*K)
    _method = method;
    _timestep = timestep;
    _T = T;


    timestep2 = 0.5 * _timestep * _timestep;
    
    if (!_velocityPtr)
      GenerateVelocities();
    Energy(true); // compute gradients
 
    for (int i = 1; i <= n; i++) {
      FOR_ATOMS_OF_MOL (a, _mol) {
        if (!_constraints.IsFixed(a->GetIdx())) {
          if (_method & OBFF_ANALYTICAL_GRADIENT)
            force = GetGradient(&*a) + _constraints.GetGradient(a->GetIdx());
          else
            force = NumericalDerivative(&*a) + _constraints.GetGradient(a->GetIdx());

          pos = a->GetVector();
          coordIdx = (a->GetIdx() - 1) * 3;
	  
	  // a(i) = F(i) / m
	  accel = force / (1000 * a->GetAtomicMass());

	  // x(i+1) = x(i) + v(i) /\t + 0.5 a(i) /\t^2
	  pos.SetX(pos.x() + _velocityPtr[coordIdx]   * _timestep + accel.x() * timestep2);
	  pos.SetY(pos.y() + _velocityPtr[coordIdx+1] * _timestep + accel.y() * timestep2);
	  pos.SetZ(pos.z() + _velocityPtr[coordIdx+2] * _timestep + accel.z() * timestep2);
	  a->SetVector(pos);
	  
          // v(i+.5) = v(i) + 0.5 a(i) /\t
	  _velocityPtr[coordIdx]   += 0.5 * accel.x() * _timestep;
	  _velocityPtr[coordIdx+1] += 0.5 * accel.y() * _timestep;
	  _velocityPtr[coordIdx+2] += 0.5 * accel.z() * _timestep;

	  Energy(true); // compute gradients
          
	  if (_method & OBFF_ANALYTICAL_GRADIENT)
            force = GetGradient(&*a) + _constraints.GetGradient(a->GetIdx());
          else
            force = NumericalDerivative(&*a) + _constraints.GetGradient(a->GetIdx());
          
	  // a(i+1) = F(i+1) / m
	  accel = force / (1000 * a->GetAtomicMass());

	  _velocityPtr[coordIdx]   += 0.5 * accel.x() * _timestep;
	  _velocityPtr[coordIdx+1] += 0.5 * accel.y() * _timestep;
	  _velocityPtr[coordIdx+2] += 0.5 * accel.z() * _timestep;

	}
      }
      if (i % 10 == 0) 
        CorrectVelocities();
    }
  }
 
  //////////////////////////////////////////////////////////////////////////////////
  //
  // Vector analyse
  //
  //////////////////////////////////////////////////////////////////////////////////
 
  double OBForceField::VectorLengthDerivative(vector3 &a, vector3 &b)
  {
    vector3 vab, drab;
    double rab;
    
    vab = a - b;
    rab = vab.length();
    if (rab < 0.1) // atoms are too close to each other
      {
        vab.randomUnitVector();
        vab *= 0.1; // move the atoms a small, random distance apart
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
    if (IsNearZero(rab, 1.0e-3))
      rab = 0.01;
    rcb = vcb.length();
    if (IsNearZero(rcb, 1.0e-3))
      rcb = 0.01;
    rab2 = rab * rab;
    rcb2 = rcb * rcb;
    
    abcb = dot(vab, vcb) / (rab * rcb);
    abcb2 = 1.0 - abcb * abcb;
    theta = acos(abcb) * RAD_TO_DEG;
    
    if (IsNearZero(abcb2)) {
      a = VZero;
      b = VZero;
      c = VZero;
      return 180.0;
    }

    a = (vcb * rab * rcb - (vab / rab) * dot(vab, vcb) * rcb) / (sqrt(abcb2) * rab2 * rcb2);
    c = -((vcb / rcb) * dot(vab, vcb) * rab - vab * rab * rcb) / (sqrt(abcb2) * rab2 * rcb2);
    b = -a - c;

    a *= (1.0 / DEG_TO_RAD);
    b *= (1.0 / DEG_TO_RAD);
    c *= (1.0 / DEG_TO_RAD);

    return theta;
  }
 
  double OBForceField::VectorOOPDerivative(vector3 &a, vector3 &b, vector3 &c, vector3 &d)
  {
    // This is adapted from http://scidok.sulb.uni-saarland.de/volltexte/2007/1325/pdf/Dissertation_1544_Moll_Andr_2007.pdf
    // Many thanks to Andreas Moll and the BALLView developers for this

    double angle = Point2PlaneAngle(d, a, b, c) * DEG_TO_RAD; // in the reference, this is x_ijk;l

    vector3 vbc, vbd, vba; // normalized vectors between the atoms
    double rbc, rbd, rba; // distances between the atoms
    double angleABC; // angle between abc (in the reference, this is phi_ijk
    
    // We should add some double-checks for small values here to make sure we don't "explode" the gradient    
    vbc = b - c;
    vbd = b - d;
    vba = b - a;
    rbc = vbc.length();
    if (IsNearZero(rbc, 1.0e-3))
      rbc = 0.01;
    rbd = vbd.length();
    if (IsNearZero(rbd, 1.0e-3))
      rbd = 0.01;
    rba = vba.length();
    if (IsNearZero(rba, 1.0e-3))
      rba = 0.01;
    vbc.normalize();
    vbd.normalize();
    vba.normalize();
    
    angleABC = vectorAngle(vbc, vba) * DEG_TO_RAD;
    
    // A few shortcuts
    double sinChi, sinABC, cosABC;
    sinChi = sin(angle);
    sinABC = sin(angleABC);
    if (IsNearZero(sinABC, 1.0e-3))
      sinABC = 1.0e-3;
    cosABC = cos(angleABC);
    vector3 bcbd = cross(vbc, vbd);
    vector3 bdba = cross(vbd, vba);
    vector3 babd = cross(vba, vbd);
    
    // These terms lack the dE term (or Oijkl in the reference)
    a = (bcbd + (-1.0*vba + (vbc * cosABC*sinChi/sinABC))) / (rba * sinABC);
    c = (bdba + (-1.0*vbc + (vbc * cosABC*sinChi/sinABC))) / (rbc * sinABC);
    d = ((babd / sinABC) - vbd * sinChi) / rbd;
    
    b = -1.0 * (a + c + d);
    
    return angle;
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
    if (IsNearZero(dotAbbcBccd) || !isfinite(tor)) { // stop any NaN or infinity
      tor = 1.0e-3; // rather than NaN
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
      abbc_bccd2 = 1e-5;
    
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

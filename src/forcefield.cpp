/**********************************************************************
forcefield.cpp - Handle OBForceField class.

Copyright (C) 2006-2007 by Tim Vandermeersch <tim.vandermeersch@gmail.com>
Some portions Copyright (C) 2007-2008 by Geoffrey Hutchison

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

#include <set>

#include <openbabel/forcefield.h>

#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/ring.h>
#include <openbabel/obiter.h>
#include <openbabel/math/matrix3x3.h>
#include <openbabel/rotamer.h>
#include <openbabel/rotor.h>
#include <openbabel/grid.h>
#include <openbabel/griddata.h>
#include <openbabel/elements.h>
#include "rand.h"

using namespace std;

namespace OpenBabel
{
#if defined(__CYGWIN__) || defined(__MINGW32__)
  // macro to implement static OBPlugin::PluginMapType& Map()
  PLUGIN_CPP_FILE(OBForceField)
#endif

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
      OBForceField* pFF = OBForceField::FindForceField("MMFF94");

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
      OBForceField* pFF = OBForceField::FindForceField("MMFF94");

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
      OBForceField* pFF = OBForceField::FindForceField("MMFF94");

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
      pFF->ConjugateGradients(1000);
      \endcode

      Minimize a ligand molecule in a binding pocket.
      \code
      #include <openbabel/forcefield.h>
      #include <openbabel/mol.h>

      OBMol mol;

      //
      // Read the pocket + ligand (initial guess for position) into mol...
      //

      OBBitVec pocket; // set the bits with atoms indexes for the pocket to 1...
      OBBitVec ligand; // set the bits with atoms indexes for the ligand to 1...

      OBForceField* pFF = OBForceField::FindForceField("MMFF94");

      // Make sure we have a valid pointer
      if (!pFF)
      // exit...

      pFF->SetLogFile(&cerr);
      pFF->SetLogLevel(OBFF_LOGLVL_LOW);

      // Fix the binding pocket atoms
      OBFFConstraints constraints;
      FOR_ATOMS_OF_MOL (a, mol) {
      if (pocket.BitIsSet(a->GetIdx())
      constraints.AddAtomConstraint(a->GetIdx());
      }

      // Specify the interacting groups. The pocket atoms are fixed, so there
      // is no need to calculate intra- and inter-molecular interactions for
      // the binding pocket.
      pFF->AddIntraGroup(ligand); // bonded interactions in the ligand
      pFF->AddInterGroup(ligand); // non-bonded between ligand-ligand atoms
      pFF->AddInterGroups(ligand, pocket); // non-bonded between ligand and pocket atoms

      // We pass the constraints as argument for Setup()
      if (!pFF->Setup(mol, constraints)) {
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
        if (((a == parameter[idx].a) && (b == parameter[idx].b)) ||
            ((a == parameter[idx].b) && (b == parameter[idx].a)))
          return idx;

    if (!d)
      for (unsigned int idx=0; idx < parameter.size(); idx++)
        if (((a == parameter[idx].a) && (b == parameter[idx].b) && (c == parameter[idx].c)) ||
            ((a == parameter[idx].c) && (b == parameter[idx].b) && (c == parameter[idx].a)))
          return idx;

    for (unsigned int idx=0; idx < parameter.size(); idx++)
      if (((a == parameter[idx].a) && (b == parameter[idx].b) &&
           (c == parameter[idx].c) && (d == parameter[idx].d)) ||
          ((a == parameter[idx].d) && (b == parameter[idx].c) &&
           (c == parameter[idx].b) && (d == parameter[idx].a)))
        return idx;

    return -1;
  }

  OBFFParameter* OBForceField::GetParameter(int a, int b, int c, int d,
                                            vector<OBFFParameter> &parameter)
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
        if (((a == parameter[idx].a) && (b == parameter[idx].b)) ||
            ((a == parameter[idx].b) && (b == parameter[idx].a))) {
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
      if (((a == parameter[idx].a) && (b == parameter[idx].b) &&
           (c == parameter[idx].c) && (d == parameter[idx].d)) ||
          ((a == parameter[idx].d) && (b == parameter[idx].c) &&
           (c == parameter[idx].b) && (d == parameter[idx].a))) {
        par = &parameter[idx];
        return par;
      }

    return NULL;
  }

  OBFFParameter* OBForceField::GetParameter(const char* a, const char* b, const char* c,
                                            const char* d, vector<OBFFParameter> &parameter)
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
        if (((_a == parameter[idx]._a) && (_b == parameter[idx]._b)) ||
            ((_a == parameter[idx]._b) && (_b == parameter[idx]._a))) {
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
      if (((_a == parameter[idx]._a) && (_b == parameter[idx]._b) &&
           (_c == parameter[idx]._c) && (_d == parameter[idx]._d)) ||
          ((_a == parameter[idx]._d) && (_b == parameter[idx]._c) &&
           (_c == parameter[idx]._b) && (_d == parameter[idx]._a))) {
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

  OBFFConstraints OBForceField::_constraints = OBFFConstraints(); // define static data variable
  unsigned int OBForceField::_fixAtom = 0; // define static data variable
  unsigned int OBForceField::_ignoreAtom = 0; // define static data variable

  OBFFConstraints& OBForceField::GetConstraints()
  {
    return _constraints;
  }

  void OBForceField::SetConstraints(OBFFConstraints& constraints)
  {
    if (!(_constraints.GetIgnoredBitVec() == constraints.GetIgnoredBitVec())) {
      _constraints = constraints;
      if (!SetupCalculations()) {
        _validSetup = false;
        return;
      }
    } else {
      _constraints = constraints;
    }

    _constraints.Setup(_mol);
  }

  void OBForceField::SetFixAtom(int index)
  {
    _fixAtom = index;
  }

  void OBForceField::UnsetFixAtom()
  {
    _fixAtom = 0;
  }

  void OBForceField::SetIgnoreAtom(int index)
  {
    _ignoreAtom = index; // remember the index
  }

  void OBForceField::UnsetIgnoreAtom()
  {
    _ignoreAtom = 0;
  }

  bool OBForceField::IgnoreCalculation(int a, int b)
  {
    if (!_ignoreAtom)
      return false;

    if (_ignoreAtom == a)
      return true;
    if (_ignoreAtom == b)
      return true;

    return false;
  }

  bool OBForceField::IgnoreCalculation(int a, int b, int c)
  {
    if (OBForceField::IgnoreCalculation(a, c))
      return true;
    if (_ignoreAtom == b)
      return true;

    return false;
  }

  bool OBForceField::IgnoreCalculation(int a, int b, int c, int d)
  {
    if (OBForceField::IgnoreCalculation(a, b, c))
      return true;
    if (_ignoreAtom == d)
      return true;

    return false;
  }

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
      grad += i->GetGradient(a);

    return grad;
  }

  void OBFFConstraints::Clear()
  {
    _constraints.clear();
    _ignored.Clear();
    _fixed.Clear();
    _Xfixed.Clear();
    _Yfixed.Clear();
    _Zfixed.Clear();
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
      if ( (i->type == OBFF_CONST_DISTANCE) || (i->type == OBFF_CONST_ANGLE) ||
           (i->type == OBFF_CONST_TORSION) ) {
        vector3 da, db, dc, dd;
        double delta, delta2, rab, theta, dE;

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

          da = (i->a)->GetVector();
          db = (i->b)->GetVector();
          dc = (i->c)->GetVector();
          dd = (i->d)->GetVector();

          theta = OBForceField::VectorTorsionDerivative(da, db, dc, dd);
          if (!isfinite(theta))
            theta = 1.0e-7;

          theta = DEG_TO_RAD * (theta + 180.0 - i->constraint_value);
          constraint_energy +=  0.001 * i->factor * (1.0 + cos(theta));
          dE = 0.001 * i->factor * sin(theta);

          i->grada = dE * da;
          i->gradb = dE * db;
          i->gradc = dE * dc;
          i->gradd = dE * dd;
          break;
        default:
          break;
        }
      }
    return constraint_energy;
  }

  void OBFFConstraints::DeleteConstraint(int index)
  {
    vector<OBFFConstraint>::iterator i;
    int n = 0;

    for (i = _constraints.begin(); i != _constraints.end(); ++n, ++i) {
      if (n == index) {
        if (i->type == OBFF_CONST_IGNORE)
          _ignored.SetBitOff(i->ia);
        if (i->type == OBFF_CONST_ATOM)
          _fixed.SetBitOff(i->ia);
        if (i->type == OBFF_CONST_ATOM_X)
          _Xfixed.SetBitOff(i->ia);
        if (i->type == OBFF_CONST_ATOM_Y)
          _Yfixed.SetBitOff(i->ia);
        if (i->type == OBFF_CONST_ATOM_Z)
          _Zfixed.SetBitOff(i->ia);


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
    _ignored.SetBitOn(a);

    OBFFConstraint constraint;
    constraint.type = OBFF_CONST_IGNORE; // constraint type
    constraint.ia   = a; // atom to fix
    _constraints.push_back(constraint);
  }

  void OBFFConstraints::AddAtomConstraint(int a)
  {
    _fixed.SetBitOn(a);

    OBFFConstraint constraint;
    constraint.type = OBFF_CONST_ATOM; // constraint type
    constraint.ia   = a; // atom to fix
    constraint.factor = _factor;
    _constraints.push_back(constraint);
  }

  void OBFFConstraints::AddAtomXConstraint(int a)
  {
    _Xfixed.SetBitOn(a);

    OBFFConstraint constraint;
    constraint.type = OBFF_CONST_ATOM_X; // constraint type
    constraint.ia   = a; // atom to fix
    constraint.factor = _factor;
    _constraints.push_back(constraint);
  }

  void OBFFConstraints::AddAtomYConstraint(int a)
  {
    _Yfixed.SetBitOn(a);

    OBFFConstraint constraint;
    constraint.type = OBFF_CONST_ATOM_Y; // constraint type
    constraint.ia   = a; // atom to fix
    constraint.factor = _factor;
    _constraints.push_back(constraint);
  }

  void OBFFConstraints::AddAtomZConstraint(int a)
  {
    _Zfixed.SetBitOn(a);

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
    if (index >= _constraints.size())
      return 0;

    return _constraints[index].type;
  }

  double OBFFConstraints::GetConstraintValue(int index) const
  {
    if (index >= _constraints.size())
      return 0;

    return _constraints[index].constraint_value;
  }

  int OBFFConstraints::GetConstraintAtomA(int index) const
  {
    if (index >= _constraints.size())
      return 0;

    return _constraints[index].ia;
  }

  int OBFFConstraints::GetConstraintAtomB(int index) const
  {
    if (index >= _constraints.size())
      return 0;

    return _constraints[index].ib;
  }

  int OBFFConstraints::GetConstraintAtomC(int index) const
  {
    if (index >= _constraints.size())
      return 0;

    return _constraints[index].ic;
  }

  int OBFFConstraints::GetConstraintAtomD(int index) const
  {
    if (index >= _constraints.size())
      return 0;

    return _constraints[index].id;
  }

  bool OBFFConstraints::IsIgnored(int index)
  {
    return _ignored.BitIsSet(index);
  }

  bool OBFFConstraints::IsFixed(int index)
  {
    return _fixed.BitIsSet(index);
  }

  bool OBFFConstraints::IsXFixed(int index)
  {
    return _Xfixed.BitIsSet(index);
  }

  bool OBFFConstraints::IsYFixed(int index)
  {
    return _Yfixed.BitIsSet(index);
  }

  bool OBFFConstraints::IsZFixed(int index)
  {
    return _Zfixed.BitIsSet(index);
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
      OBFFLog("IDX\tTYPE\tRING\n");

      FOR_ATOMS_OF_MOL (a, _mol) {
        snprintf(_logbuf, BUFF_SIZE, "%d\t%s\t%s\n", a->GetIdx(), a->GetType(),
          (a->IsInRing() ? (a->IsAromatic() ? "AR" : "AL") : "NO"));
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
        snprintf(_logbuf, BUFF_SIZE, "%d\t%f\n", a->GetIdx(), a->GetPartialCharge());
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
        snprintf(_logbuf, BUFF_SIZE, "%d\t%f\n", a->GetIdx(), a->GetPartialCharge());
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
        snprintf(_logbuf, BUFF_SIZE, "%d\t<%8.3f, %8.3f, %8.3f>\n", a->GetIdx(), _velocityPtr[a->GetIdx()],
                 _velocityPtr[a->GetIdx()+1], _velocityPtr[a->GetIdx()+2]);
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
      _gradientPtr = NULL;
      _grad1 = NULL;
    }

    if (IsSetupNeeded(mol)) {
      _mol = mol;
      _ncoords = _mol.NumAtoms() * 3;

      if (_velocityPtr)
        delete [] _velocityPtr;
      _velocityPtr = NULL;

      if (_gradientPtr)
        delete [] _gradientPtr;
      _gradientPtr = new double[_ncoords];

      if (_mol.NumAtoms() && _constraints.Size())
        _constraints.Setup(_mol);

      _mol.SetSSSRPerceived(false);
      _mol.DeleteData(OBGenericDataType::TorsionData); // bug #1954233

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

    } else {
      if (_validSetup) {
        PrintTypes();
        PrintFormalCharges();
        PrintPartialCharges();
        SetCoordinates(mol);
        return true;
      } else {
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
      _gradientPtr = NULL;
    }

    if (IsSetupNeeded(mol)) {
      _mol = mol;
      _ncoords = _mol.NumAtoms() * 3;

      if (_velocityPtr)
        delete [] _velocityPtr;
      _velocityPtr = NULL;

      if (_gradientPtr)
        delete [] _gradientPtr;
      _gradientPtr = new double[_ncoords];

      _constraints = constraints;
      if (_mol.NumAtoms() && _constraints.Size())
        _constraints.Setup(_mol);

      _mol.SetSSSRPerceived(false);
      _mol.DeleteData(OBGenericDataType::TorsionData); // bug #1954233

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

    } else {
      if (_validSetup) {
        if (!(_constraints.GetIgnoredBitVec() == constraints.GetIgnoredBitVec())) {
          _constraints = constraints;
          if (!SetupCalculations()) {
            _validSetup = false;
            return false;
          }
        } else {
          _constraints = constraints;
        }

        _constraints.Setup(_mol);
        SetCoordinates(mol);
        return true;
      } else {
        return false;
      }
    }

    _validSetup = true;
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
      // if we have Fe or Cu the atom type depends on the formal
      // charge, so a new setup must be forced anyway
      if (atom->GetAtomicNum() == 26 || atom->GetAtomicNum() == 29)
        return true;
      if (atom->GetAtomicNum() != (mol.GetAtom(atom->GetIdx()))->GetAtomicNum())
        return true;
      if (atom->GetExplicitDegree() != (mol.GetAtom(atom->GetIdx()))->GetExplicitDegree())
        return true;
    }
    FOR_BONDS_OF_MOL (bond, _mol) {
      if (bond->GetBondOrder() != (mol.GetBond(bond->GetIdx()))->GetBondOrder())
        return true;
      if (bond->GetBeginAtom()->GetAtomicNum()
          != (mol.GetBond(bond->GetIdx()))->GetBeginAtom()->GetAtomicNum()
          || bond->GetEndAtom()->GetAtomicNum()
          != (mol.GetBond(bond->GetIdx()))->GetEndAtom()->GetAtomicNum())
        return true;
    }

    return false;
  }

  bool OBForceField::GetAtomTypes(OBMol &mol)
  {
    if (_mol.NumAtoms() != mol.NumAtoms())
      return false;

    FOR_ATOMS_OF_MOL (intAtom, _mol) {
      OBAtom *extAtom = mol.GetAtom(intAtom->GetIdx());

      if (extAtom->HasData("FFAtomType")) {
        OBPairData *data = (OBPairData*) extAtom->GetData("FFAtomType");
        data->SetValue(intAtom->GetType());
      } else {
        OBPairData *data = new OBPairData();
       	data->SetAttribute("FFAtomType");
        data->SetValue(intAtom->GetType());
        extAtom->SetData(data);
      }
    }

    return true;
  }

  bool OBForceField::GetPartialCharges(OBMol &mol)
  {
    if (_mol.NumAtoms() != mol.NumAtoms())
      return false;

    ostringstream out;
    FOR_ATOMS_OF_MOL (intAtom, _mol) {
      OBAtom *extAtom = mol.GetAtom(intAtom->GetIdx());

      out.str("");
      out << intAtom->GetPartialCharge();
      if (extAtom->HasData("FFPartialCharge")) {
        OBPairData *data = (OBPairData*) extAtom->GetData("FFPartialCharge");
        data->SetValue(out.str());
      } else {
        OBPairData *data = new OBPairData();
       	data->SetAttribute("FFPartialCharge");
        data->SetValue(out.str());
        extAtom->SetData(data);
      }
    }

    return true;
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

    if (!mol.HasData(OBGenericDataType::ConformerData))
      mol.SetData(new OBConformerData);
    OBConformerData *cd = (OBConformerData*) mol.GetData(OBGenericDataType::ConformerData);
    cd->SetEnergies(_energies);

    vector<vector3> forces;
    vector<vector<vector3> > confForces;
    for (unsigned int i = 0; i < _mol.NumAtoms(); ++i) {
      const int coordIdx = i * 3;
      forces.push_back(vector3(_gradientPtr[coordIdx],
                               _gradientPtr[coordIdx+1], _gradientPtr[coordIdx+2]));
    }
    confForces.push_back(forces);
    cd->SetForces(confForces);

    return true;
  }

  bool OBForceField::GetConformers(OBMol &mol)
  {
    //    OBAtom *atom;

    if (_mol.NumAtoms() != mol.NumAtoms())
      return false;

    /*
      FOR_ATOMS_OF_MOL (a, _mol) {
      atom = mol.GetAtom(a->GetIdx());
      atom->SetVector(a->GetVector());
      }
    */

    //Copy conformer information
    if (_mol.NumConformers() > 0) {
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
      mol.SetConformer(_current_conformer);

      if (!mol.HasData(OBGenericDataType::ConformerData))
        mol.SetData(new OBConformerData);
      OBConformerData *cd = (OBConformerData*) mol.GetData(OBGenericDataType::ConformerData);
      cd->SetEnergies(_energies);

      //mol.SetEnergies(_energies);
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
      SetupPointers(); // update pointers to atom positions in the OBFFCalculation objects
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

  int OBForceField::SystematicRotorSearchInitialize(unsigned int geomSteps, bool sampleRingBonds)
  {
    if (!_validSetup)
      return 0;

    OBRotorList rl;
    OBRotamerList rotamers;
    OBRotorIterator ri;
    OBRotor *rotor;

    _origLogLevel = _loglvl;

    OBBitVec fixed = _constraints.GetFixedBitVec();
    rl.SetFixAtoms(fixed);
    rl.Setup(_mol, sampleRingBonds);
    rotamers.SetBaseCoordinateSets(_mol);
    rotamers.Setup(_mol, rl);

    IF_OBFF_LOGLVL_LOW {
      OBFFLog("\nS Y S T E M A T I C   R O T O R   S E A R C H\n\n");
      snprintf(_logbuf, BUFF_SIZE, "  NUMBER OF ROTATABLE BONDS: %lu\n", (unsigned long)rl.Size());
      OBFFLog(_logbuf);

      unsigned long int combinations = 1;
      for (rotor = rl.BeginRotor(ri); rotor;
           rotor = rl.NextRotor(ri)) {
        combinations *= rotor->GetResolution().size();
      }
      snprintf(_logbuf, BUFF_SIZE, "  NUMBER OF POSSIBLE ROTAMERS: %lu\n", combinations);
      OBFFLog(_logbuf);
    }

    _current_conformer = 0;

    if (!rl.Size()) { // only one conformer
      IF_OBFF_LOGLVL_LOW
        OBFFLog("  GENERATED ONLY ONE CONFORMER\n\n");

      ConjugateGradients(geomSteps); // final energy minimizatin for best conformation

      return 1; // there are no more conformers
    }

    OBRotorKeys rotorKeys;
    rotor = rl.BeginRotor(ri);
    for (unsigned int i = 1; i < rl.Size() + 1; ++i, rotor = rl.NextRotor(ri)) // foreach rotor
      rotorKeys.AddRotor(rotor->GetResolution().size());

    rotamers.AddRotamer(rotorKeys.GetKey());
    while (rotorKeys.Next())
      rotamers.AddRotamer(rotorKeys.GetKey());

    rotamers.ExpandConformerList(_mol, _mol.GetConformers());

    IF_OBFF_LOGLVL_LOW {
      snprintf(_logbuf, BUFF_SIZE, "  GENERATED %d CONFORMERS\n\n", _mol.NumConformers());
      OBFFLog(_logbuf);
      OBFFLog("CONFORMER     ENERGY\n");
      OBFFLog("--------------------\n");
    }

    _energies.clear();

    return _mol.NumConformers();
  }

  bool OBForceField::SystematicRotorSearchNextConformer(unsigned int geomSteps)
  {
    if (!_validSetup)
      return 0;

    if (_current_conformer >=  _mol.NumConformers()) { // done
      // Select conformer with lowest energy
      int best_conformer = 0;
      for (int i = 0; i < _mol.NumConformers(); i++) {
        if (_energies[i] < _energies[best_conformer])
          best_conformer = i;
      }

      IF_OBFF_LOGLVL_LOW {
        snprintf(_logbuf, BUFF_SIZE, "\n  CONFORMER %d HAS THE LOWEST ENERGY\n\n",  best_conformer + 1);
        OBFFLog(_logbuf);
      }

      _mol.SetConformer(best_conformer);
      SetupPointers(); // update pointers to atom positions in the OBFFCalculation objects
      _current_conformer = best_conformer;

      return false;
    }

    _mol.SetConformer(_current_conformer); // select conformer
    SetupPointers(); // update pointers to atom positions in the OBFFCalculation objects

    _loglvl = OBFF_LOGLVL_NONE;
    ConjugateGradients(geomSteps); // energy minimization for conformer
    _loglvl = _origLogLevel;

    _energies.push_back(Energy(false)); // calculate and store energy

    IF_OBFF_LOGLVL_LOW {
      snprintf(_logbuf, BUFF_SIZE, "   %3d   %20.3f\n", (_current_conformer + 1), _energies[_current_conformer]);
      OBFFLog(_logbuf);
    }

    _current_conformer++;
    return true;
  }

  void OBForceField::SystematicRotorSearch(unsigned int geomSteps, bool sampleRingBonds)
  {
    if (SystematicRotorSearchInitialize(geomSteps, sampleRingBonds))
      while (SystematicRotorSearchNextConformer(geomSteps)) {}
  }

  int OBForceField::FastRotorSearch(bool permute)
  {
    if (_mol.NumRotors() == 0)
      return 0;

    int origLogLevel = _loglvl;

    // Remove all conformers (e.g. from previous conformer generators) except for current conformer
    double *initialCoord = new double [_mol.NumAtoms() * 3]; // initial state
    double *store_initial = new double [_mol.NumAtoms() * 3]; // store the initial state
    memcpy((char*)initialCoord,(char*)_mol.GetCoordinates(),sizeof(double)*3*_mol.NumAtoms());
    memcpy((char*)store_initial,(char*)_mol.GetCoordinates(),sizeof(double)*3*_mol.NumAtoms());
    std::vector<double *> newConfs(1, initialCoord);
    _mol.SetConformers(newConfs);

    _energies.clear(); // Wipe any energies from previous conformer generators

    OBRotorList rl;
    OBBitVec fixed = _constraints.GetFixedBitVec();
    rl.SetFixAtoms(fixed);
    rl.SetQuiet();
    rl.Setup(_mol);

    OBRotorIterator ri;
    OBRotamerList rotamerlist;
    rotamerlist.SetBaseCoordinateSets(_mol);
    rotamerlist.Setup(_mol, rl);

    // Start with all of the rotors in their 0 position
    // (perhaps instead I should set them randomly?)
    std::vector<int> init_rotorKey(rl.Size() + 1, 0);
    std::vector<int> rotorKey(init_rotorKey);

    unsigned int j, minj;
    double currentE, minE, best_minE;

    double *verybestconf = new double [_mol.NumAtoms() * 3]; // store the best conformer to date
    double *bestconf = new double [_mol.NumAtoms() * 3]; // store the best conformer to date in the current permutation
    double *minconf = new double [_mol.NumAtoms() * 3];  // store the best conformer for the current rotor
    memcpy((char*)bestconf,(char*)_mol.GetCoordinates(),sizeof(double)*3*_mol.NumAtoms());

    double energy_offset;
    // Can take shortcut later, as 4 components of the energy will be constant
    rotamerlist.SetCurrentCoordinates(_mol, rotorKey);
    SetupPointers();
    energy_offset = E_Bond(false) + E_Angle(false) + E_StrBnd(false) + E_OOP(false);

    // This function relies on the fact that Rotors are ordered from the most
    // central to the most peripheral (due to CompareRotors in rotor.cpp)
    std::vector<OBRotor *> vrotors;
    OBRotor *rotor;
    for (rotor = rl.BeginRotor(ri); rotor; rotor = rl.NextRotor(ri))
      vrotors.push_back(rotor);

    // The permutations are ordered so that the first 2 permutations cover the
    // combinations of 2, and the first 6 permutations cover the combinations of 3
    const char permutations[24*4] = {0,1,2,3, 1,0,2,3, 0,2,1,3, 1,2,0,3, 2,0,1,3, 2,1,0,3,
                                     0,1,3,2, 0,2,3,1, 0,3,1,2, 0,3,2,1, 1,0,3,2, 1,2,3,0,
                                     1,3,0,2, 1,3,2,0, 2,0,3,1, 2,1,3,0, 2,3,0,1, 2,3,1,0,
                                     3,0,1,2, 3,0,2,1, 3,1,0,2, 3,1,2,0, 3,2,0,1, 3,2,1,0};
    const char factorial[5] = {0, 1, 2, 6, 24};

    char num_rotors_to_permute, num_permutations;
    if (permute)
      num_rotors_to_permute = (char)std::min<size_t> (4, vrotors.size());
    else
      num_rotors_to_permute = 1; // i.e. just use the original order
    num_permutations = factorial[num_rotors_to_permute];

    // Initialize reordered_rotors - the order in which to test rotors
    std::vector<unsigned int> reordered_rotors(vrotors.size());
    for (int i=0; i<vrotors.size(); ++i)
      reordered_rotors[i] = i;

    std::set<unsigned int> seen;
    best_minE = DBL_MAX;
    for (int N=0; N<num_permutations; ++N) {
      for (int i=0; i<num_rotors_to_permute; ++i)
        reordered_rotors.at(i) = *(permutations + N*4 + i);

      rotorKey = init_rotorKey;
      _mol.SetCoordinates(store_initial);
      bool quit = false;

      for (int i=0; i<reordered_rotors.size(); ++i) {
        unsigned int idx = reordered_rotors[i];
        rotor = vrotors.at(idx);

        minE = DBL_MAX;

        for (j = 0; j < rotor->GetResolution().size(); j++) { // For each rotor position
          // Note: we could do slightly better by skipping the rotor position we already
          //       tested in the last loop (position 0 at the moment). Note that this
          //       isn't as simple as just changing the loop starting point to j = 1.
          _mol.SetCoordinates(bestconf);
          rotorKey[idx + 1] = j;
          rotamerlist.SetCurrentCoordinates(_mol, rotorKey);
          SetupPointers();

          currentE = E_VDW(false) + E_Torsion(false) + E_Electrostatic(false);

          if (currentE < minE) {
            minE = currentE;
            minj = j;
            memcpy((char*)minconf,(char*)_mol.GetCoordinates(),sizeof(double)*3*_mol.NumAtoms());
          }
        } // Finished testing all positions of this rotor
        rotorKey[idx + 1] = minj;

        if (i==4) { // Check whether this rotorKey has already been chosen
          // Create a hash of the rotorKeys (given that the max value of any rotorKey is 11 from torlib.txt)
          unsigned int hash = rotorKey[1] + rotorKey[2]*12 + rotorKey[3]*12*12 + rotorKey[4]*12*12*12;

          if (seen.find(hash) == seen.end()) // Not seen before
            seen.insert(hash);
          else { // Already seen - no point continuing
            quit = true;
            break;
          }
        }

        memcpy((char*)bestconf,(char*)minconf,sizeof(double)*3*_mol.NumAtoms());
      } // end of this permutation
      if (!quit) {
          if (minE < best_minE) {
            best_minE = minE;
            memcpy((char*)verybestconf,(char*)bestconf,sizeof(double)*3*_mol.NumAtoms());
          }
      }

    } // end of final permutation

    _mol.SetCoordinates(verybestconf);
    SetupPointers();

    delete [] store_initial;
    delete [] bestconf;
    delete [] verybestconf;
    delete [] minconf;

    return true;
  }

  void OBForceField::RandomRotorSearchInitialize(unsigned int conformers, unsigned int geomSteps,
                                                 bool sampleRingBonds)
  {
    if (!_validSetup)
      return;

    OBRotorList rl;
    OBRotamerList rotamers;
    OBRotorIterator ri;
    OBRotor *rotor;

    OBRandom generator;
    generator.TimeSeed();
    _origLogLevel = _loglvl;

    if (_mol.GetCoordinates() == NULL)
      return;

    OBBitVec fixed = _constraints.GetFixedBitVec();
    rl.SetFixAtoms(fixed);
    rl.Setup(_mol, sampleRingBonds);
    rotamers.SetBaseCoordinateSets(_mol);
    rotamers.Setup(_mol, rl);

    IF_OBFF_LOGLVL_LOW {
      OBFFLog("\nR A N D O M   R O T O R   S E A R C H\n\n");
      snprintf(_logbuf, BUFF_SIZE, "  NUMBER OF ROTATABLE BONDS: %lu\n", (unsigned long)rl.Size());
      OBFFLog(_logbuf);

      unsigned long int combinations = 1;
      for (rotor = rl.BeginRotor(ri); rotor;
           rotor = rl.NextRotor(ri)) {
        combinations *= rotor->GetResolution().size();
      }
      snprintf(_logbuf, BUFF_SIZE, "  NUMBER OF POSSIBLE ROTAMERS: %lu\n", combinations);
      OBFFLog(_logbuf);
    }

    _current_conformer = 0;

    if (!rl.Size()) { // only one conformer
      IF_OBFF_LOGLVL_LOW
        OBFFLog("  GENERATED ONLY ONE CONFORMER\n\n");

      _loglvl = OBFF_LOGLVL_NONE;
      ConjugateGradients(geomSteps); // energy minimization for conformer
      _loglvl = _origLogLevel;

      return;
    }

    std::vector<int> rotorKey(rl.Size() + 1, 0); // indexed from 1

    for (unsigned int c = 0; c < conformers; ++c) {
      rotor = rl.BeginRotor(ri);
      for (unsigned int i = 1; i < rl.Size() + 1; ++i, rotor = rl.NextRotor(ri)) {
        // foreach rotor
        rotorKey[i] = generator.NextInt() % rotor->GetResolution().size();
      }
      rotamers.AddRotamer(rotorKey);
    }

    rotamers.ExpandConformerList(_mol, _mol.GetConformers());

    IF_OBFF_LOGLVL_LOW {
      snprintf(_logbuf, BUFF_SIZE, "  GENERATED %d CONFORMERS\n\n", _mol.NumConformers());
      OBFFLog(_logbuf);
      OBFFLog("CONFORMER     ENERGY\n");
      OBFFLog("--------------------\n");
    }

    _energies.clear();
  }

  bool OBForceField::RandomRotorSearchNextConformer(unsigned int geomSteps)
  {
    if (!_validSetup)
      return 0;

    if (_current_conformer >=  _mol.NumConformers()) { // done
      // Select conformer with lowest energy
      int best_conformer = 0;
      for (int i = 0; i < _mol.NumConformers(); i++) {
        if (_energies[i] < _energies[best_conformer])
          best_conformer = i;
      }

      IF_OBFF_LOGLVL_LOW {
        snprintf(_logbuf, BUFF_SIZE, "\n  CONFORMER %d HAS THE LOWEST ENERGY\n\n",  best_conformer + 1);
        OBFFLog(_logbuf);
      }

      _mol.SetConformer(best_conformer);
      SetupPointers(); // update pointers to atom positions in the OBFFCalculation objects
      _current_conformer = best_conformer;

      return false;
    }

    _mol.SetConformer(_current_conformer); // select conformer
    SetupPointers(); // update pointers to atom positions in the OBFFCalculation objects

    _loglvl = OBFF_LOGLVL_NONE;
    ConjugateGradients(geomSteps); // energy minimization for conformer
    _loglvl = _origLogLevel;

    _energies.push_back(Energy(false)); // calculate and store energy

    IF_OBFF_LOGLVL_LOW {
      snprintf(_logbuf, BUFF_SIZE, "   %3d      %8.3f\n", (_current_conformer + 1), _energies[_current_conformer]);
      OBFFLog(_logbuf);
    }

    _current_conformer++;
    return true;
  }

  void OBForceField::RandomRotorSearch(unsigned int conformers, unsigned int geomSteps,
                                       bool sampleRingBonds)
  {
    RandomRotorSearchInitialize(conformers, geomSteps, sampleRingBonds);
    while (RandomRotorSearchNextConformer(geomSteps)) {}
  }

  void Reweight(std::vector< std::vector <double> > &rotorWeights,
                std::vector<int> rotorKey, double bonus)
  {
    double fraction, minWeight, maxWeight;
    bool improve = (bonus > 0.0);

    for (unsigned int i = 1; i < rotorWeights.size() - 1; ++i) {
			if (rotorKey[i] == -1)
				continue; // don't rotate

      if (improve && rotorWeights[i][rotorKey[i]] > 0.999 - bonus)
        continue;
      if (!improve && rotorWeights[i][rotorKey[i]] < 0.001 - bonus) // bonus < 0
        continue;

      // Check to make sure we don't kill some poor weight
      minWeight = maxWeight = rotorWeights[i][0];
      for (unsigned int j = 1; j < rotorWeights[i].size(); ++j) {
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

      for (unsigned int j = 0; j < rotorWeights[i].size(); ++j) {
        if (j == rotorKey[i])
          rotorWeights[i][j] += bonus;
        else
          rotorWeights[i][j] -= fraction;
      }
    }
  }


  void OBForceField::WeightedRotorSearch(unsigned int conformers, unsigned int geomSteps,
                                         bool sampleRingBonds)
  {
    if (!_validSetup)
      return;

    OBRotorList rl;
    OBRotamerList rotamers;
    OBRotorIterator ri;
    OBRotor *rotor;

    OBRandom generator;
    generator.TimeSeed();
    int origLogLevel = _loglvl;

    if (_mol.GetCoordinates() == NULL)
      return;

    double *initialCoord = new double [_mol.NumAtoms() * 3]; // initial state
    memcpy((double*)initialCoord,(double*)_mol.GetCoordinates(),sizeof(double)*3*_mol.NumAtoms());
    vector<double *> newConfs(1, initialCoord);
    _mol.SetConformers(newConfs);
    _current_conformer = 0;
    _mol.SetConformer(_current_conformer);
    SetupPointers();

    OBBitVec fixed = _constraints.GetFixedBitVec();
    rl.SetFixAtoms(fixed);
    rl.Setup(_mol, sampleRingBonds);
    rotamers.SetBaseCoordinateSets(_mol);
    rotamers.Setup(_mol, rl);

    IF_OBFF_LOGLVL_LOW {
      OBFFLog("\nW E I G H T E D   R O T O R   S E A R C H\n\n");
      snprintf(_logbuf, BUFF_SIZE, "  NUMBER OF ROTATABLE BONDS: %lu\n", (unsigned long)rl.Size());
      OBFFLog(_logbuf);

      unsigned long int combinations = 1;
      for (rotor = rl.BeginRotor(ri); rotor;
           rotor = rl.NextRotor(ri)) {
        combinations *= rotor->GetResolution().size();
      }
      snprintf(_logbuf, BUFF_SIZE, "  NUMBER OF POSSIBLE ROTAMERS: %lu\n", combinations);
      OBFFLog(_logbuf);
    }


    _energies.clear(); // Wipe any energies from previous conformer generators

    if (!rl.Size()) { // only one conformer
      IF_OBFF_LOGLVL_LOW
        OBFFLog("  GENERATED ONLY ONE CONFORMER\n\n");

      _loglvl = OBFF_LOGLVL_NONE;
      ConjugateGradients(geomSteps); // energy minimization for conformer
      _loglvl = origLogLevel;
      _energies.push_back(Energy(false));

      return;
    }

    _energies.push_back(Energy(false)); // Store the energy of the original conf

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
    IF_OBFF_LOGLVL_LOW
      OBFFLog("  INITIAL WEIGHTING OF ROTAMERS...\n\n");

    rotor = rl.BeginRotor(ri);
    for (unsigned int i = 1; i < rl.Size() + 1; ++i, rotor = rl.NextRotor(ri)) {
      rotorKey[i] = -1; // no rotation (new in 2.2)
    }

    rotor = rl.BeginRotor(ri);

    for (unsigned int i = 1; i < rl.Size() + 1; ++i, rotor = rl.NextRotor(ri)) {
      // foreach rotor
      energies.clear();
      for (unsigned int j = 0; j < rotor->GetResolution().size(); j++) {
        // foreach rotor position
        _mol.SetCoordinates(initialCoord);
        rotorKey[i] = j;
        rotamers.SetCurrentCoordinates(_mol, rotorKey);
        SetupPointers(); // update pointers to atom positions in the OBFFCalculation objects

        _loglvl = OBFF_LOGLVL_NONE;
        ConjugateGradients(geomSteps); // energy minimization for conformer
        _loglvl = origLogLevel;
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
      rotorKey[i] = -1; // back to the previous setting before we go to another rotor

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

    int best_conformer=-1;
    //    double penalty; // for poor performance
    double randFloat; // generated random number -- used to pick a rotor
    double total; // used to calculate the total probability

    // Start with the current coordinates
    bestE = worstE = Energy(false);

    // Now we actually test some weightings
    IF_OBFF_LOGLVL_LOW {
      snprintf(_logbuf, BUFF_SIZE, "  GENERATED %d CONFORMERS\n\n", conformers);
      OBFFLog(_logbuf);
      OBFFLog("CONFORMER     ENERGY\n");
      OBFFLog("--------------------\n");
    }

    double defaultRotor = 1.0/sqrt((double)rl.Size());
    unsigned c = 0;
    while (c < conformers) {
      _mol.SetCoordinates(initialCoord);

      // Choose the rotor key based on current weightings
      rotor = rl.BeginRotor(ri);
      for (unsigned int i = 1; i < rl.Size() + 1; ++i, rotor = rl.NextRotor(ri)) {
        // foreach rotor
        rotorKey[i] = -1; // default = don't change dihedral
        randFloat = generator.NextFloat();
        if (randFloat < defaultRotor) // should we just leave this rotor with default setting?
          continue;

        randFloat = generator.NextFloat();
        total = 0.0;
        for (unsigned int j = 0; j < rotor->GetResolution().size(); j++) {
          if (randFloat > total && randFloat < (total+ rotorWeights[i][j])) {
            rotorKey[i] = j;
            break;
          }
          else
            total += rotorWeights[i][j];
        }
      }

      //FIXME: for now, allow even invalid ring conformers
      rotamers.SetCurrentCoordinates(_mol, rotorKey);
      ++c;

      SetupPointers(); // update pointers to atom positions in the OBFFCalculation objects

      _loglvl = OBFF_LOGLVL_NONE;
      ConjugateGradients(geomSteps); // energy minimization for conformer
      _loglvl = origLogLevel;
      currentE = Energy(false);
      _energies.push_back(currentE);
      double *confCoord = new double [_mol.NumAtoms() * 3]; // initial state
      memcpy((char*)confCoord,(char*)_mol.GetCoordinates(),sizeof(double)*3*_mol.NumAtoms());
      _mol.AddConformer(confCoord);

      IF_OBFF_LOGLVL_LOW {
        snprintf(_logbuf, BUFF_SIZE, "   %3d      %8.3f\n", c + 1, currentE);
        OBFFLog(_logbuf);
      }

      if (!isfinite(currentE))
        continue;

      if (currentE < bestE) {
        bestE = currentE;
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
      snprintf(_logbuf, BUFF_SIZE, "\n  LOWEST ENERGY: %8.3f\n\n",
               bestE);
      OBFFLog(_logbuf);
    }

    // debugging output to see final weightings of each rotor setting
    IF_OBFF_LOGLVL_HIGH {
      OBFFLog("Final Weights: \n");
      for (unsigned int i = 1; i < rotorWeights.size() - 1; ++i) {
        snprintf(_logbuf, BUFF_SIZE, " Weight: %d", i);
        OBFFLog(_logbuf);
        for (unsigned int j = 0; j < rotorWeights[i].size(); ++j) {
          snprintf(_logbuf, BUFF_SIZE, " %8.3f", rotorWeights[i][j]);
          OBFFLog(_logbuf);
        }
        OBFFLog("\n");
      }
    }

    _current_conformer = best_conformer; // Initial coords are stored in _vconf[0]
    _mol.SetConformer(_current_conformer);
    SetupPointers(); // update pointers to atom positions in the OBFFCalculation objects
  }

  void OBForceField::DistanceGeometry()
  {
    if (!_validSetup)
      return;

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
          snprintf(_logbuf, BUFF_SIZE, " %8.4f ", matrix[i][j]);
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

            if (l_ac < (l_ab - l_bc)) {// l_ac >= l_ab - l_bc
              l_ac = l_ab - l_bc;
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
          snprintf(_logbuf, BUFF_SIZE, " %8.4f ", matrix[i][j]);
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

        srand(static_cast <unsigned int> (time(NULL)));
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
      OBFFLog("RANDOM DISTANCE MATRIX BETWEEN LIMITS\n\n");
      for (i=0; i<N; i++) {
        OBFFLog("[");
        for (j=0; j<N; j++) {
          snprintf(_logbuf, BUFF_SIZE, " %8.4f ", matrix[i][j]);
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
      OBFFLog("METRIC MATRIX\n\n");
      for (i=0; i<N; i++) {
        OBFFLog("[");
        for (j=0; j<N; j++) {
          snprintf(_logbuf, BUFF_SIZE, " %8.4f ", G[i][j]);
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
        snprintf(_logbuf, BUFF_SIZE, "%8.4f ", eigenvalues[i]);
        OBFFLog(_logbuf);
      }
      OBFFLog("\n");

      OBFFLog("EIGENVECTORS OF METRIC MATRIX\n\n");
      for (i=0; i<N; i++) {
        OBFFLog("[");
        for (j=0; j<N; j++) {
          snprintf(_logbuf, BUFF_SIZE, " %8.4f ", eigenvectors[i][j]);
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
  // Cut-off
  //
  //////////////////////////////////////////////////////////////////////////////////

  void OBForceField::UpdatePairsSimple()
  {
    _vdwpairs.Clear();
    _elepairs.Clear();

    const unsigned int numAtoms = _mol.NumAtoms();
    const unsigned int numPairs = numAtoms * (numAtoms - 1) / 2;
    _vdwpairs.Resize(numPairs);
    _elepairs.Resize(numPairs);

    //! \todo set the criteria as squared values
    //  from what the user supplies
    double rvdwSquared = SQUARE(_rvdw);
    double releSquared = SQUARE(_rele);

    unsigned int pairIndex = -1;
    FOR_PAIRS_OF_MOL(p, _mol) {
      ++pairIndex;

      OBAtom *a = _mol.GetAtom((*p)[0]);
      OBAtom *b = _mol.GetAtom((*p)[1]);

      // Check whether or not this interaction is included
      if (HasGroups()) {
        bool isIncludedPair = false;
        for (size_t i=0; i < _interGroup.size(); ++i) {
          if (_interGroup[i].BitIsSet(a->GetIdx()) &&
              _interGroup[i].BitIsSet(b->GetIdx())) {
            isIncludedPair = true;
            break;
          }
        }
        if (!isIncludedPair) {
          for (size_t i=0; i < _interGroups.size(); ++i) {
            if (_interGroups[i].first.BitIsSet(a->GetIdx()) &&
                _interGroups[i].second.BitIsSet(b->GetIdx())) {
              isIncludedPair = true;
              break;
            }
            if (_interGroups[i].first.BitIsSet(b->GetIdx()) &&
                _interGroups[i].second.BitIsSet(a->GetIdx())) {
              isIncludedPair = true;
              break;
            }
          }
        }
        if (!isIncludedPair) {
          continue;
        }
      }

      // Get the distance squared btwn a and b
      // Not currently a method and this can be vectorized
      double ab[3];
      VectorSubtract(a->GetCoordinate(), b->GetCoordinate(), ab);
      double rabSq = 0.0;
      for (int j = 0; j < 3; ++j) {
        rabSq += SQUARE(ab[j]);
      }

      // update vdw pairs
      if (rabSq < rvdwSquared) {
        _vdwpairs.SetBitOn(pairIndex);
      } else {
        _vdwpairs.SetBitOff(pairIndex);
      }
      // update electrostatic pairs
      if (rabSq < releSquared) {
        _elepairs.SetBitOn(pairIndex);
      } else {
        _elepairs.SetBitOff(pairIndex);
      }
    }

    /*
      IF_OBFF_LOGLVL_LOW {
      snprintf(_logbuf, BUFF_SIZE, "UPDATE VDW PAIRS: %d --> %d (VDW), %d (ELE) \n", i+1,
      _vdwpairs.CountBits(), _elepairs.CountBits());
      OBFFLog(_logbuf);
      }
    */
  }

  unsigned int OBForceField::GetNumPairs()
  {
    unsigned int i = 1;
    FOR_PAIRS_OF_MOL(p, _mol)
      i++;

    return i;
  }

  unsigned int OBForceField::GetNumElectrostaticPairs()
  {
    return _elepairs.CountBits();
  }

  unsigned int OBForceField::GetNumVDWPairs()
  {
    return _vdwpairs.CountBits();
  }

  //////////////////////////////////////////////////////////////////////////////////
  //
  // Interaction groups
  //
  //////////////////////////////////////////////////////////////////////////////////


  void OBForceField::AddIntraGroup(OBBitVec &group)
  {
    _intraGroup.push_back(group);
  }

  void OBForceField::AddInterGroup(OBBitVec &group)
  {
    _interGroup.push_back(group);
  }

  void OBForceField::AddInterGroups(OBBitVec &group1, OBBitVec &group2)
  {
    pair<OBBitVec,OBBitVec> groups;
    groups.first = group1;
    groups.second = group2;
    _interGroups.push_back(groups);
  }

  void OBForceField::ClearGroups()
  {
    _intraGroup.clear();
    _interGroup.clear();
    _interGroups.clear();
  }

  bool OBForceField::HasGroups()
  {
    if (!_intraGroup.empty())
      return true;

    if (!_interGroup.empty())
      return true;

    if (!_interGroups.empty())
      return true;

    return false;
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

    dir = atom->GetVector() - orig_xyz;
    atom->SetVector(orig_xyz);

    // cutoff accuracy
    if (dir.length() < 1.0e-8)
      return VZero;

    return dir;
  }

  //
  // Based on the ghemical code (conjgrad.cpp)
  //
  // Implements several enhancements:
  // 1) Switch to smarter line search method (e.g., Newton's method in 1D)
  //  x(n+1) = x(n) - F(x) / F'(x)   -- can be done numerically
  // 2) Switch to line search for one step of the entire molecule!
  //   This dramatically cuts down on the number of Energy() calls.
  //  (and is more correct anyway)
  double OBForceField::Newton2NumLineSearch(double *direction)
  {
    double e_n1, e_n2, e_n3;
    double *origCoords = new double [_ncoords];

    double opt_step = 0.0;
    double opt_e = _e_n1; // get energy calculated by sd or cg
    const double def_step = 0.025; // default step
    const double max_step = 4.5; // don't go too far

    double sum = 0.0;
    for (unsigned int c = 0; c < _ncoords; ++c) {
      if (isfinite(direction[c])) {
        sum += direction[c] * direction[c];
      } else {
        // make sure we don't have NaN or infinity
        direction[c] = 0.0;
      }
    }

    double scale = sqrt(sum);
    if (IsNearZero(scale)) {
      //      cout << "WARNING: too small \"scale\" at Newton2NumLineSearch" << endl;
      scale = 1.0e-70; // try to avoid "division by zero" conditions
    }

    double step = def_step / scale;
    double max_scl = max_step / scale;

    // Save the current position, before we take a step
    memcpy((char*)origCoords,(char*)_mol.GetCoordinates(),sizeof(double)*_ncoords);

    int newton = 0;
    while (true) {
      // Take step X(n) + step
      LineSearchTakeStep(origCoords, direction, step);
      e_n1 = Energy(false) + _constraints.GetConstraintEnergy();

      if (e_n1 < opt_e) {
        opt_step = step;
        opt_e = e_n1;
      }

      if (newton++ > 3)
        break;
      double delta = step * 0.001;

      // Take step X(n) + step + delta
      LineSearchTakeStep(origCoords, direction, step+delta);
      e_n2 = Energy(false) + _constraints.GetConstraintEnergy();

      // Take step X(n) + step + delta * 2.0
      LineSearchTakeStep(origCoords, direction, step+delta*2.0);
      e_n3 = Energy(false) + _constraints.GetConstraintEnergy();

      double denom = e_n3 - 2.0 * e_n2 + e_n1; // f'(x)
      if (denom != 0.0) {
        step = fabs(step - delta * (e_n2 - e_n1) / denom);
        if (step > max_scl) {
          //          cout << "WARNING: damped steplength " << step << " to " << max_scl << endl;
          step = max_scl;
        }
      } else {
        break;
      }
    }

    if (opt_step == 0.0) { // if we still don't have any valid steplength, try a very small step
      step = 0.001 * def_step / scale;

      // Take step X(n) + step
      LineSearchTakeStep(origCoords, direction, step);
      e_n1 = Energy(false) + _constraints.GetConstraintEnergy();

      if (e_n1 < opt_e) {
        opt_step = step;
        opt_e = e_n1;
      }

    }

    // Take optimal step
    LineSearchTakeStep(origCoords, direction, opt_step);

    //    cout << " scale: " << scale << " step: " << opt_step*scale << " maxstep " << max_scl*scale << endl;

    delete [] origCoords;

    return opt_step * scale;
  }

  void OBForceField::LineSearchTakeStep(double* origCoords, double *direction, double step)
  {
    double *currentCoords = _mol.GetCoordinates();

    for (unsigned int c = 0; c < _ncoords; ++c) {
      if (isfinite(direction[c])) {
        currentCoords[c] = origCoords[c] + direction[c] * step;
      }
    }
  }

  double OBForceField::LineSearch(double *currentCoords, double *direction)
  {
    unsigned int numCoords = _mol.NumAtoms() * 3;
    double e_n1, e_n2, step, alpha, tempStep;
    double *lastStep = new double [numCoords];

    alpha = 0.0; // Scale factor along direction vector
    step = 0.2;
    double trustRadius = 0.75; // don't move too far at once

    // The initial energy should be precomputed
    e_n1 = _e_n1; // Energy(false) + _constraints.GetConstraintEnergy();

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
          else if (tempStep < -trustRadius) // negative big step
            currentCoords[c] -= trustRadius;
          else
            currentCoords[c] += tempStep;
        }
      }

      e_n2 = Energy(false) + _constraints.GetConstraintEnergy();

      // convergence criteria: A higher precision here
      // only takes longer with the same result.
      if (IsNear(e_n2, e_n1, 1.0e-3)) {
        alpha += step;
        break;
      }

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

    atom->SetVector(9.0, 9.0, 0.0);
    e_n1 = atom->x() * atom->x() + 2 * (atom->y() * atom->y());

    IF_OBFF_LOGLVL_LOW {
      OBFFLog("\nV A L I D A T E   S T E E P E S T   D E S C E N T\n\n");
      snprintf(_logbuf, BUFF_SIZE, "STEPS = %d\n\n",  steps);
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
        snprintf(_logbuf, BUFF_SIZE, " %4d    %8.3f    %8.3f\n", i, e_n2, e_n1);
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

    delete atom;
  }

  void OBForceField::ValidateConjugateGradients(int steps)
  {
    OBAtom *atom = new OBAtom;
    vector3 grad1, grad2, dir1, dir2;
    double e_n1, e_n2;
    double g2g2, g1g1, g2g1;
    bool firststep;

    firststep = true;
    atom->SetVector(9.0, 9.0, 0.0);
    e_n1 = atom->x() * atom->x() + 2 * (atom->y() * atom->y());

    IF_OBFF_LOGLVL_LOW {
      OBFFLog("\nV A L I D A T E   C O N J U G A T E   G R A D I E N T S\n\n");
      snprintf(_logbuf, BUFF_SIZE, "STEPS = %d\n\n",  steps);
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
          snprintf(_logbuf, BUFF_SIZE, " %4d    %8.3f    %8.3f\n", i, e_n2, e_n1);
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
          snprintf(_logbuf, BUFF_SIZE, " %4d    %8.3f    %8.3f\n", i, e_n2, e_n1);
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

    delete atom;
  }

  void OBForceField::SteepestDescentInitialize(int steps, double econv, int method)
  {
    if (!_validSetup)
      return;

    _nsteps = steps;
    _cstep = 0;
    _econv = econv;
    _gconv = 1.0e-2; // gradient convergence (0.1) squared

    if (_cutoff)
      UpdatePairsSimple(); // Update the non-bonded pairs (Cut-off)

    _e_n1 = Energy() + _constraints.GetConstraintEnergy();

    IF_OBFF_LOGLVL_LOW {
      OBFFLog("\nS T E E P E S T   D E S C E N T\n\n");
      snprintf(_logbuf, BUFF_SIZE, "STEPS = %d\n\n",  steps);
      OBFFLog(_logbuf);
      OBFFLog("STEP n       E(n)         E(n-1)    \n");
      OBFFLog("------------------------------------\n");
      snprintf(_logbuf, BUFF_SIZE, " %4d    %8.3f      ----\n", _cstep, _e_n1);
      OBFFLog(_logbuf);
    }

  }

  bool OBForceField::SteepestDescentTakeNSteps(int n)
  {
    if (!_validSetup)
      return 0;

    _ncoords = _mol.NumAtoms() * 3;
    double e_n2;
    vector3 dir;
    double maxgrad; // for convergence

    for (int i = 1; i <= n; i++) {
      _cstep++;
      maxgrad = 1.0e20;

      FOR_ATOMS_OF_MOL (a, _mol) {
        unsigned int idx = a->GetIdx();
        unsigned int coordIdx = (idx - 1) * 3;

        if (_constraints.IsFixed(idx) || (_fixAtom == idx) || (_ignoreAtom == idx)) {
          _gradientPtr[coordIdx] = 0.0;
          _gradientPtr[coordIdx+1] = 0.0;
          _gradientPtr[coordIdx+2] = 0.0;
        } else {
          if (!HasAnalyticalGradients()) {
            // use numerical gradients
            dir = NumericalDerivative(&*a) + _constraints.GetGradient(a->GetIdx());
          } else {
            // use analytical gradients
            dir = GetGradient(&*a) + _constraints.GetGradient(a->GetIdx());
          }

          // check to see how large the gradients are
          if (dir.length_2() < maxgrad)
            maxgrad = dir.length_2();

          if (!_constraints.IsXFixed(idx))
            _gradientPtr[coordIdx] = dir.x();
          else
            _gradientPtr[coordIdx] = 0.0;

          if (!_constraints.IsYFixed(idx))
            _gradientPtr[coordIdx+1] = dir.y();
          else
            _gradientPtr[coordIdx+1] = 0.0;

          if (!_constraints.IsZFixed(idx))
            _gradientPtr[coordIdx+2] = dir.z();
          else
            _gradientPtr[coordIdx+2] = 0.0;
        }
      }
      // perform a linesearch
      switch (_linesearch) {
      case LineSearchType::Newton2Num:
        Newton2NumLineSearch(_gradientPtr);
        break;
      default:
      case LineSearchType::Simple:
        LineSearch(_mol.GetCoordinates(), _gradientPtr);
        break;
      }
      e_n2 = Energy() + _constraints.GetConstraintEnergy();

      if ((_cstep % _pairfreq == 0) && _cutoff)
        UpdatePairsSimple(); // Update the non-bonded pairs (Cut-off)

      IF_OBFF_LOGLVL_LOW {
        if (_cstep % 10 == 0) {
          snprintf(_logbuf, BUFF_SIZE, " %4d    %8.5f    %8.5f\n", _cstep, e_n2, _e_n1);
          OBFFLog(_logbuf);
        }
      }

      if (IsNear(e_n2, _e_n1, _econv)
          && (maxgrad < _gconv)) { // gradient criteria (0.1) squared
        IF_OBFF_LOGLVL_LOW
          OBFFLog("    STEEPEST DESCENT HAS CONVERGED\n");
        return false;
      }

      if (_nsteps == _cstep) {
        return false;
      }

      _e_n1 = e_n2;
    }

    return true;  // no convergence reached
  }

  void OBForceField::SteepestDescent(int steps, double econv, int method)
  {
    if (steps > 0) {
      SteepestDescentInitialize(steps, econv, method);
      SteepestDescentTakeNSteps(steps);
    }
  }

  void OBForceField::ConjugateGradientsInitialize(int steps, double econv,
                                                  int method)
  {
    if (!_validSetup || steps==0)
      return;

    double e_n2;
    vector3 dir;

    _cstep = 0;
    _nsteps = steps;
    _econv = econv;
    _gconv = 1.0e-2; // gradient convergence (0.1) squared
    _ncoords = _mol.NumAtoms() * 3;

    if (_cutoff)
      UpdatePairsSimple(); // Update the non-bonded pairs (Cut-off)

    _e_n1 = Energy() + _constraints.GetConstraintEnergy();

    IF_OBFF_LOGLVL_LOW {
      OBFFLog("\nC O N J U G A T E   G R A D I E N T S\n\n");
      snprintf(_logbuf, BUFF_SIZE, "STEPS = %d\n\n",  steps);
      OBFFLog(_logbuf);
      OBFFLog("STEP n     E(n)       E(n-1)    \n");
      OBFFLog("--------------------------------\n");
    }

    if (_grad1 != NULL)
      delete [] _grad1;
    _grad1 = new double[_ncoords];
    memset(_grad1, '\0', sizeof(double)*_ncoords);

    // Take the first step (same as steepest descent because there is no
    // gradient from the previous step.
    FOR_ATOMS_OF_MOL (a, _mol) {
      unsigned int idx = a->GetIdx();
      unsigned int coordIdx = (idx - 1) * 3;

      if (_constraints.IsFixed(idx) || (_fixAtom == idx) || (_ignoreAtom == idx)) {
        _gradientPtr[coordIdx] = 0.0;
        _gradientPtr[coordIdx+1] = 0.0;
        _gradientPtr[coordIdx+2] = 0.0;
      } else {
        if (!HasAnalyticalGradients()) {
          // use numerical gradients
          dir = NumericalDerivative(&*a) + _constraints.GetGradient(a->GetIdx());
        } else {
          // use analytical gradients
          dir = GetGradient(&*a) + _constraints.GetGradient(a->GetIdx());
        }

        if (!_constraints.IsXFixed(idx))
          _gradientPtr[coordIdx] = dir.x();
        else
          _gradientPtr[coordIdx] = 0.0;

        if (!_constraints.IsYFixed(idx))
          _gradientPtr[coordIdx+1] = dir.y();
        else
          _gradientPtr[coordIdx+1] = 0.0;

        if (!_constraints.IsZFixed(idx))
          _gradientPtr[coordIdx+2] = dir.z();
        else
          _gradientPtr[coordIdx+2] = 0.0;
      }
    }
    // perform a linesearch
    switch (_linesearch) {
    case LineSearchType::Newton2Num:
      Newton2NumLineSearch(_gradientPtr);
      break;
    default:
    case LineSearchType::Simple:
      LineSearch(_mol.GetCoordinates(), _gradientPtr);
      break;
    }
    e_n2 = Energy() + _constraints.GetConstraintEnergy();

    IF_OBFF_LOGLVL_LOW {
      snprintf(_logbuf, BUFF_SIZE, " %4d    %8.3f    %8.3f\n", 1, e_n2, _e_n1);
      OBFFLog(_logbuf);
    }

    // save the direction and energy
    memcpy(_grad1, _gradientPtr, sizeof(double)*_ncoords);
    _e_n1 = e_n2;
  }

  bool OBForceField::ConjugateGradientsTakeNSteps(int n)
  {
    if (!_validSetup)
      return 0;

    double e_n2;
    double g2g2, g1g1, beta;
    vector3 grad2, dir2;
    vector3 grad1, dir1; // temporaries to perform dot product, etc.
    double maxgrad; // for convergence

    if (_ncoords != _mol.NumAtoms() * 3)
      return false;

    e_n2 = 0.0;

    for (int i = 1; i <= n; i++) {
      _cstep++;
      maxgrad = 1.0e20;

      FOR_ATOMS_OF_MOL (a, _mol) {
        unsigned int idx = a->GetIdx();
        unsigned int coordIdx = (a->GetIdx() - 1) * 3;

        if (_constraints.IsFixed(idx) || (_fixAtom == idx) || (_ignoreAtom == idx)) {
          _grad1[coordIdx] = 0.0;
          _grad1[coordIdx+1] = 0.0;
          _grad1[coordIdx+2] = 0.0;
        } else {
          if (!HasAnalyticalGradients()) {
            // use numerical gradients
            grad2 = NumericalDerivative(&*a) + _constraints.GetGradient(a->GetIdx());
          } else {
            // use analytical gradients
            grad2 = GetGradient(&*a) + _constraints.GetGradient(a->GetIdx());
          }

          // Fletcher-Reeves formula for Beta
          // http://en.wikipedia.org/wiki/Nonlinear_conjugate_gradient_method
          // NOTE: We make sure to reset and use the steepest descent direction
          //   after NumAtoms steps
          if (_cstep % _mol.NumAtoms() != 0) {
            g2g2 = dot(grad2, grad2);
            grad1 = vector3(_grad1[coordIdx], _grad1[coordIdx+1], _grad1[coordIdx+2]);
            g1g1 = dot(grad1, grad1);
            beta = g2g2 / g1g1;
            grad2 += beta * grad1;
          }

          // check to see how large the gradients are
          if (grad2.length_2() < maxgrad)
            maxgrad = grad2.length_2();

          if (!_constraints.IsXFixed(idx))
            _grad1[coordIdx] = grad2.x();
          else
            _grad1[coordIdx] = 0.0;

          if (!_constraints.IsYFixed(idx))
            _grad1[coordIdx+1] = grad2.y();
          else
            _grad1[coordIdx+1] = 0.0;

          if (!_constraints.IsZFixed(idx))
            _grad1[coordIdx+2] = grad2.z();
          else
            _grad1[coordIdx+2] = 0.0;
        }
      }
      // perform a linesearch
      switch (_linesearch) {
      case LineSearchType::Newton2Num:
        Newton2NumLineSearch(_grad1);
        break;
      default:
      case LineSearchType::Simple:
        LineSearch(_mol.GetCoordinates(), _grad1);
        break;
      }
      // save the direction
      memcpy(_gradientPtr, _grad1, sizeof(double)*_ncoords);

      e_n2 = Energy() + _constraints.GetConstraintEnergy();

      if ((_cstep % _pairfreq == 0) && _cutoff)
        UpdatePairsSimple(); // Update the non-bonded pairs (Cut-off)

      if (IsNear(e_n2, _e_n1, _econv)
          && (maxgrad < _gconv)) { // gradient criteria (0.1) squared
        IF_OBFF_LOGLVL_LOW {
          snprintf(_logbuf, BUFF_SIZE, " %4d    %8.3f    %8.3f\n", _cstep, e_n2, _e_n1);
          OBFFLog(_logbuf);
          OBFFLog("    CONJUGATE GRADIENTS HAS CONVERGED\n");
        }
        return false;
      }

      IF_OBFF_LOGLVL_LOW {
        if (_cstep % 10 == 0) {
          snprintf(_logbuf, BUFF_SIZE, " %4d    %8.3f    %8.3f\n", _cstep, e_n2, _e_n1);
          OBFFLog(_logbuf);
        }
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
    if (steps > 1) // ConjugateGradientsInitialize takes the first step
      ConjugateGradientsTakeNSteps(steps);
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
        e_orig += E_Bond(false);
      if (terms & OBFF_EANGLE)
        e_orig += E_Angle(false);
      if (terms & OBFF_ESTRBND)
        e_orig += E_StrBnd(false);
      if (terms & OBFF_ETORSION)
        e_orig += E_Torsion(false);
      if (terms & OBFF_EOOP)
        e_orig += E_OOP(false);
      if (terms & OBFF_EVDW)
        e_orig += E_VDW(false);
      if (terms & OBFF_EELECTROSTATIC)
        e_orig += E_Electrostatic(false);
    }

    // X direction
    atom->SetVector(va.x() + delta, va.y(), va.z());

    if (terms & OBFF_ENERGY)
      e_plus_delta = Energy(false);
    else {
      e_plus_delta = 0.0;
      if (terms & OBFF_EBOND)
        e_plus_delta += E_Bond(false);
      if (terms & OBFF_EANGLE)
        e_plus_delta += E_Angle(false);
      if (terms & OBFF_ESTRBND)
        e_plus_delta += E_StrBnd(false);
      if (terms & OBFF_ETORSION)
        e_plus_delta += E_Torsion(false);
      if (terms & OBFF_EOOP)
        e_plus_delta += E_OOP(false);
      if (terms & OBFF_EVDW)
        e_plus_delta += E_VDW(false);
      if (terms & OBFF_EELECTROSTATIC)
        e_plus_delta += E_Electrostatic(false);
    }

    dx = (e_plus_delta - e_orig) / delta;

    // Y direction
    atom->SetVector(va.x(), va.y() + delta, va.z());

    if (terms & OBFF_ENERGY)
      e_plus_delta = Energy(false);
    else {
      e_plus_delta = 0.0;
      if (terms & OBFF_EBOND)
        e_plus_delta += E_Bond(false);
      if (terms & OBFF_EANGLE)
        e_plus_delta += E_Angle(false);
      if (terms & OBFF_ESTRBND)
        e_plus_delta += E_StrBnd(false);
      if (terms & OBFF_ETORSION)
        e_plus_delta += E_Torsion(false);
      if (terms & OBFF_EOOP)
        e_plus_delta += E_OOP(false);
      if (terms & OBFF_EVDW)
        e_plus_delta += E_VDW(false);
      if (terms & OBFF_EELECTROSTATIC)
        e_plus_delta += E_Electrostatic(false);
    }

    dy = (e_plus_delta - e_orig) / delta;

    // Z direction
    atom->SetVector(va.x(), va.y(), va.z() + delta);

    if (terms & OBFF_ENERGY)
      e_plus_delta = Energy(false);
    else {
      e_plus_delta = 0.0;
      if (terms & OBFF_EBOND)
        e_plus_delta += E_Bond(false);
      if (terms & OBFF_EANGLE)
        e_plus_delta += E_Angle(false);
      if (terms & OBFF_ESTRBND)
        e_plus_delta += E_StrBnd(false);
      if (terms & OBFF_ETORSION)
        e_plus_delta += E_Torsion(false);
      if (terms & OBFF_EOOP)
        e_plus_delta += E_OOP(false);
      if (terms & OBFF_EVDW)
        e_plus_delta += E_VDW(false);
      if (terms & OBFF_EELECTROSTATIC)
        e_plus_delta += E_Electrostatic(false);
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
  // This means that if we divide the force (coming from the FF) by 1000x its mass in amu,
  // we get the acceleration in A ps^-2.

  // gromacs user manual page 17
  void OBForceField::GenerateVelocities()
  {
    cout << "OBForceField::GenerateVelocities()" << endl;
    OBRandom generator;
    generator.TimeSeed();
    _ncoords = _mol.NumAtoms() * 3;
    int velocityIdx;
    double velocity;

    _velocityPtr = new double[_ncoords];
    memset(_velocityPtr, '\0', sizeof(double)*_ncoords);

    FOR_ATOMS_OF_MOL (a, _mol) {
      if (!_constraints.IsFixed(a->GetIdx()) || (_fixAtom == a->GetIdx()) || (_ignoreAtom == a->GetIdx())) {
        velocityIdx = (a->GetIdx() - 1) * 3;

        // add twelve random numbers between 0.0 and 1.0,
        // subtract 6.0 from their sum, multiply with sqrt(kT/m)
        if (!_constraints.IsXFixed(a->GetIdx())) {
          velocity = 0.0;
          for (int i=0; i < 12; ++i)
            velocity += generator.NextFloat();
          velocity -= 6.0;
          velocity *= sqrt((GAS_CONSTANT * _temp)/ (1000 * a->GetAtomicMass()));
          _velocityPtr[velocityIdx] = velocity; // x10: gromacs uses nm instead of A
        }

        if (!_constraints.IsYFixed(a->GetIdx())) {
          velocity = 0.0;
          for (int i=0; i < 12; ++i)
            velocity += generator.NextFloat();
          velocity -= 6.0;
          velocity *= sqrt((GAS_CONSTANT * _temp)/ (1000 * a->GetAtomicMass()));
          _velocityPtr[velocityIdx+1] = velocity; // idem
        }

        if (!_constraints.IsZFixed(a->GetIdx())) {
          velocity = 0.0;
          for (int i=0; i < 12; ++i)
            velocity += generator.NextFloat();
          velocity -= 6.0;
          velocity *= sqrt((GAS_CONSTANT * _temp)/ (1000 * a->GetAtomicMass()));
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
    double velocity, E_kin, E_kin2, factor;

    // E_kin = 0.5 * Ndf * R * T
    E_kin = _ncoords * GAS_CONSTANT * _temp;
    //cout << "E_{kin} = Ndf * R * T = " << E_kin << endl;

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
    if (!_validSetup)
      return;

    int coordIdx;
    double timestep2;
    vector3 force, pos, accel;
    _timestep = timestep;
    _temp = T;


    timestep2 = 0.5 * _timestep * _timestep;

    if (!_velocityPtr)
      GenerateVelocities();
    Energy(true); // compute gradients

    for (int i = 1; i <= n; i++) {
      FOR_ATOMS_OF_MOL (a, _mol) {
        if (!_constraints.IsFixed(a->GetIdx()) || (_fixAtom == a->GetIdx()) || (_ignoreAtom == a->GetIdx())) {
          if (HasAnalyticalGradients())
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

          if (HasAnalyticalGradients())
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

  double OBForceField::VectorBondDerivative(double *pos_i, double *pos_j,
                                            double *force_i, double *force_j)
  {
    double ij[3];
    VectorSubtract(pos_i, pos_j, ij);

    double rij = VectorLength(ij);
    if (rij < 0.1) { // atoms are too close to each other
      vector3 vij;
      vij.randomUnitVector();
      vij *= 0.1; // move the atoms a small, random distance apart
      vij.Get(ij);
      rij = 0.1;
    }
    VectorDivide(ij, rij, force_j);
    VectorMultiply(force_j, -1.0, force_i);

    return rij;
  }


  double OBForceField::VectorDistanceDerivative(const double* const pos_i, const double* const pos_j,
                                                double *force_i, double *force_j)
  {
    VectorSubtract(pos_i, pos_j, force_j);
    const double rij = VectorLength(force_j);
    const double inverse_rij = 1.0 / rij;
    VectorSelfMultiply(force_j, inverse_rij);
    VectorMultiply(force_j, -1.0, force_i);
    return rij;
  }

  double OBForceField::VectorAngleDerivative(vector3 &i, vector3 &j, vector3 &k)
  {
    // This is adapted from http://scidok.sulb.uni-saarland.de/volltexte/2007/1325/pdf/Dissertation_1544_Moll_Andr_2007.pdf
    // Many thanks to Andreas Moll and the BALLView developers for this

    vector3 v1, v2;
    vector3 n1, n2;

    // Calculate the vector between atom1 and atom2,
    // test if the vector has length larger than 0 and normalize it
    v1 = i - j;
    v2 = k - j;

    double length1 = v1.length();
    double length2 = v2.length();

    // test if the vector has length larger than 0 and normalize it
    if (IsNearZero(length1) || IsNearZero(length2)) {
      i = VZero;
      j = VZero;
      k = VZero;
      return 0.0;
    }

    // Calculate the normalized bond vectors
    double inverse_length_v1 = 1.0 / length1;
    double inverse_length_v2 = 1.0 / length2;
    v1 *= inverse_length_v1 ;
    v2 *= inverse_length_v2;

    // Calculate the cross product of v1 and v2, test if it has length unequal 0,
    // and normalize it.
    vector3 c1 = cross(v1, v2);
    double length = c1.length();
    if (IsNearZero(length)) {
      i = VZero;
      j = VZero;
      k = VZero;
      return 0.0;
    }

    c1 /= length;

    // Calculate the cos of theta and then theta
    double costheta = dot(v1, v2);
    double theta;
    if (costheta > 1.0) {
      theta = 0.0;
      costheta = 1.0;
    } else if (costheta < -1.0) {
      theta = 180.0;
      costheta = -1.0;
    } else {
      theta = RAD_TO_DEG * acos(costheta);
    }

    vector3 t1 = cross(v1, c1);
    t1.normalize();
    vector3 t2 = cross(v2, c1);
    t2.normalize();

    i = -t1 * inverse_length_v1;
    k =  t2 * inverse_length_v2;
    j = - (i + k);

    return theta;
  }

  double OBForceField::VectorAngleDerivative(double *pos_i, double *pos_j, double *pos_k,
                                             double *force_i, double *force_j, double *force_k)
  {
    // This is adapted from http://scidok.sulb.uni-saarland.de/volltexte/2007/1325/pdf/Dissertation_1544_Moll_Andr_2007.pdf
    // Many thanks to Andreas Moll and the BALLView developers for this

    // Bond vectors of the three atoms
    double ij[3], jk[3];
    VectorSubtract(pos_i, pos_j, ij);
    VectorSubtract(pos_k, pos_j, jk);

    // length of the two bonds
    double l_ij, l_jk;
    l_ij = VectorLength(ij);
    l_jk = VectorLength(jk);

    if (IsNearZero(l_ij) || IsNearZero(l_jk)) {
      VectorClear(force_i);
      VectorClear(force_j);
      VectorClear(force_k);
      return 0.0;
    }

    // normalize the bond vectors:
    VectorDivide (ij, l_ij, ij);
    VectorDivide (jk, l_jk, jk);

    // Calculate the cross product of v1 and v2, test if it has length unequal 0,
    // and normalize it.
    double c1[3];
    VectorCross(ij, jk, c1);
    double length = VectorLength(c1);
    if (IsNearZero(length)) {
      VectorClear(force_i);
      VectorClear(force_j);
      VectorClear(force_k);
      return 0.0;
    }

    VectorDivide(c1, length, c1);

    // Calculate the cos of theta and then theta
    double cos_ijk = VectorDot(ij, jk);
    double angle_ijk;
    if (cos_ijk > 1.0) {
      angle_ijk = 0.0;
      cos_ijk = 1.0;
    } else if (cos_ijk < -1.0) {
      angle_ijk = 180.0;
      cos_ijk = -1.0;
    } else {
      angle_ijk = RAD_TO_DEG * acos(cos_ijk);
    }

    double t1[3], t2[3];

    VectorCross(ij, c1, t1);
    VectorNormalize(t1);

    VectorCross(jk, c1, t2);
    VectorNormalize(t2);

    VectorDivide(t1, -l_ij, force_i);
    VectorDivide(t2,  l_jk, force_k);

    VectorAdd(force_i, force_k, force_j);
    VectorSelfMultiply(force_j, -1.0);

    return angle_ijk;
  }

  double OBForceField::VectorAngle(double *pos_i, double *pos_j, double *pos_k)
  {
    // This is adapted from http://scidok.sulb.uni-saarland.de/volltexte/2007/1325/pdf/Dissertation_1544_Moll_Andr_2007.pdf
    // Many thanks to Andreas Moll and the BALLView developers for this

    // Bond vectors of the three atoms
    double ij[3], jk[3];
    VectorSubtract(pos_i, pos_j, ij);
    VectorSubtract(pos_k, pos_j, jk);

    // length of the two bonds
    double l_ij, l_jk;
    l_ij = VectorLength(ij);
    l_jk = VectorLength(jk);

    if (IsNearZero(l_ij) || IsNearZero(l_jk)) {
      return 0.0;
    }

    // normalize the bond vectors:
    VectorDivide (ij, l_ij, ij);
    VectorDivide (jk, l_jk, jk);

    // Calculate the cross product of v1 and v2, test if it has length unequal 0,
    // and normalize it.
    double c1[3];
    VectorCross(ij, jk, c1);
    double length = VectorLength(c1);
    if (IsNearZero(length)) {
      return 0.0;
    }

    // Calculate the cos of theta and then theta
    double cos_ijk = VectorDot(ij, jk);
    double angle_ijk;
    if (cos_ijk > 1.0) {
      angle_ijk = 0.0;
      cos_ijk = 1.0;
    } else if (cos_ijk < -1.0) {
      angle_ijk = 180.0;
      cos_ijk = -1.0;
    } else {
      angle_ijk = RAD_TO_DEG * acos(cos_ijk);
    }

    return angle_ijk;
  }

  double OBForceField::VectorOOPDerivative(vector3 &i, vector3 &j, vector3 &k, vector3 &l)
  {
    // This is adapted from http://scidok.sulb.uni-saarland.de/volltexte/2007/1325/pdf/Dissertation_1544_Moll_Andr_2007.pdf
    // Many thanks to Andreas Moll and the BALLView developers for this

    // temp variables:
    double length;
    vector3 delta;

    // normal vectors of the three planes:
    vector3 an, bn, cn;

    // calculate normalized bond vectors from central atom to outer atoms:
    delta = i - j;
    length = delta.length();
    if (IsNearZero(length)) {
      i = VZero;
      j = VZero;
      k = VZero;
      l = VZero;
      return 0.0;
    }
    // normalize the bond vector:
    delta /= length;
    // store the normalized bond vector from central atom to outer atoms:
    const vector3 ji = delta;
    // store length of this bond:
    const double length_ji = length;

    delta = k - j;
    length = delta.length();
    if (IsNearZero(length)) {
      i = VZero;
      j = VZero;
      k = VZero;
      l = VZero;
      return 0.0;
    }
    // normalize the bond vector:
    delta /= length;
    // store the normalized bond vector from central atom to outer atoms:
    const vector3 jk = delta;
    // store length of this bond:
    const double length_jk = length;

    delta = l - j;
    length = delta.length();
    if (IsNearZero(length)) {
      i = VZero;
      j = VZero;
      k = VZero;
      l = VZero;
      return 0.0;
    }
    // normalize the bond vector:
    delta /= length;
    // store the normalized bond vector from central atom to outer atoms:
    const vector3 jl = delta;
    // store length of this bond:
    const double length_jl = length;

    // the normal vectors of the three planes:
    an = cross(ji, jk);
    bn = cross(jk, jl);
    cn = cross(jl, ji);

    // Bond angle ji to jk
    const double cos_theta = dot(ji, jk);
    const double theta = acos(cos_theta);
    // If theta equals 180 degree or 0 degree
    if (IsNearZero(theta) || IsNearZero(fabs(theta - M_PI))) {
      i = VZero;
      j = VZero;
      k = VZero;
      l = VZero;
      return 0.0;
    }

    const double sin_theta = sin(theta);
    const double sin_dl = dot(an, jl) / sin_theta;

    // the wilson angle:
    const double dl = asin(sin_dl);

    // In case: wilson angle equals 0 or 180 degree: do nothing
    if (IsNearZero(dl) || IsNearZero(fabs(dl - M_PI))) {
      i = VZero;
      j = VZero;
      k = VZero;
      l = VZero;
      return RAD_TO_DEG * dl;
    }

    const double cos_dl = cos(dl);

    // if wilson angle equal 90 degree: abort
    if (cos_dl < 0.0001) {
      i = VZero;
      j = VZero;
      k = VZero;
      l = VZero;
      return RAD_TO_DEG * dl;
    }

    l = (an / sin_theta - jl * sin_dl) / length_jl;
    i = ((bn + (((-ji + jk * cos_theta) * sin_dl) / sin_theta)) / length_ji) / sin_theta;
    k = ((cn + (((-jk + ji * cos_theta) * sin_dl) / sin_theta)) / length_jk) / sin_theta;
    j = -(i + k + l);

    return RAD_TO_DEG * dl;
  }

  double OBForceField::VectorOOPDerivative(double *pos_i, double *pos_j, double *pos_k, double *pos_l,
                                           double *force_i, double *force_j, double *force_k, double *force_l)
  {
    // This is adapted from http://scidok.sulb.uni-saarland.de/volltexte/2007/1325/pdf/Dissertation_1544_Moll_Andr_2007.pdf
    // Many thanks to Andreas Moll and the BALLView developers for this

    // vector lengths of the three bonds:
    double ji[3], jk[3], jl[3];
    // normal vectors of the three planes:
    double an[3], bn[3], cn[3];

    // calculate normalized bond vectors from central atom to outer atoms:
    VectorSubtract(pos_i, pos_j, ji);
    // store length of this bond:
    const double length_ji = VectorLength(ji);
    if (IsNearZero(length_ji)) {
      VectorClear(force_i);
      VectorClear(force_j);
      VectorClear(force_k);
      VectorClear(force_l);
      return 0.0;
    }
    // store the normalized bond vector from central atom to outer atoms:
    // normalize the bond vector:
    VectorDivide(ji, length_ji, ji);

    VectorSubtract(pos_k, pos_j, jk);
    const double length_jk = VectorLength(jk);
    if (IsNearZero(length_jk)) {
      VectorClear(force_i);
      VectorClear(force_j);
      VectorClear(force_k);
      VectorClear(force_l);
      return 0.0;
    }
    VectorDivide(jk, length_jk, jk);

    VectorSubtract(pos_l, pos_j, jl);
    const double length_jl = VectorLength(jl);
    if (IsNearZero(length_jl)) {
      VectorClear(force_i);
      VectorClear(force_j);
      VectorClear(force_k);
      VectorClear(force_l);
      return 0.0;
    }
    VectorDivide(jl, length_jl, jl);

    // the normal vectors of the three planes:
    VectorCross(ji, jk, an);
    VectorCross(jk, jl, bn);
    VectorCross(jl, ji, cn);

    // Bond angle ji to jk
    const double cos_theta = VectorDot(ji, jk);
    const double theta = acos(cos_theta);
    // If theta equals 180 degree or 0 degree
    if (IsNearZero(theta) || IsNearZero(fabs(theta - M_PI))) {
      VectorClear(force_i);
      VectorClear(force_j);
      VectorClear(force_k);
      VectorClear(force_l);
      return 0.0;
    }

    const double sin_theta = sin(theta);
    const double sin_dl = VectorDot(an, jl) / sin_theta;

    // the wilson angle:
    const double dl = asin(sin_dl);

    // In case: wilson angle equals 0 or 180 degree: do nothing
    if (IsNearZero(dl) || IsNearZero(fabs(dl - M_PI))) {
      VectorClear(force_i);
      VectorClear(force_j);
      VectorClear(force_k);
      VectorClear(force_l);
      return RAD_TO_DEG * dl;
    }

    const double cos_dl = cos(dl);

    // if wilson angle equal 90 degree: abort
    if (cos_dl < 0.0001) {
      VectorClear(force_i);
      VectorClear(force_j);
      VectorClear(force_k);
      VectorClear(force_l);
      return RAD_TO_DEG * dl;
    }

    double temp[3];
    /* l = (an / sin_theta - jl * sin_dl) / length_jl; */
    VectorDivide(an, sin_theta, an);
    VectorMultiply(jl, sin_dl, temp);
    VectorSubtract(an, temp, force_l);
    VectorDivide(force_l, length_jl, force_l);
    /* i = ((bn + (((-ji + jk * cos_theta) * sin_dl) / sin_theta)) / length_ji) / sin_theta; */
    VectorMultiply(jk, cos_theta, temp);
    VectorSubtract(temp, ji, temp);
    VectorSelfMultiply(temp, sin_dl);
    VectorDivide(temp, sin_theta, temp);
    VectorAdd(bn, temp, force_i);
    VectorSelfMultiply(force_i, (sin_theta/length_ji));
    /* k = ((cn + (((-jk + ji * cos_theta) * sin_dl) / sin_theta)) / length_jk) / sin_theta; */
    VectorMultiply(ji, cos_theta, temp);
    VectorSubtract(temp, jk, temp);
    VectorSelfMultiply(temp, sin_dl);
    VectorDivide(temp, sin_theta, temp);
    VectorAdd(cn, temp, force_k);
    VectorSelfMultiply(force_k, (sin_theta/length_jk));
    /* j = -(i + k + l); */
    VectorAdd(force_i, force_k, temp);
    VectorAdd(force_l, temp, temp);
    VectorMultiply(temp, -1.0, force_j);

    return RAD_TO_DEG * dl;
  }

  double OBForceField::VectorOOP(double *pos_i, double *pos_j, double *pos_k, double *pos_l)
  {
    // This is adapted from http://scidok.sulb.uni-saarland.de/volltexte/2007/1325/pdf/Dissertation_1544_Moll_Andr_2007.pdf
    // Many thanks to Andreas Moll and the BALLView developers for this

    // vector lengths of the three bonds:
    double ji[3], jk[3], jl[3];
    // normal vectors of the three planes:
    double an[3], bn[3], cn[3];

    // calculate normalized bond vectors from central atom to outer atoms:
    VectorSubtract(pos_i, pos_j, ji);
    // store length of this bond:
    const double length_ji = VectorLength(ji);
    if (IsNearZero(length_ji)) {
      return 0.0;
    }
    // store the normalized bond vector from central atom to outer atoms:
    // normalize the bond vector:
    VectorDivide(ji, length_ji, ji);

    VectorSubtract(pos_k, pos_j, jk);
    const double length_jk = VectorLength(jk);
    if (IsNearZero(length_jk)) {
      return 0.0;
    }
    VectorDivide(jk, length_jk, jk);

    VectorSubtract(pos_l, pos_j, jl);
    const double length_jl = VectorLength(jl);
    if (IsNearZero(length_jl)) {
      return 0.0;
    }
    VectorDivide(jl, length_jl, jl);

    // the normal vectors of the three planes:
    VectorCross(ji, jk, an);
    VectorCross(jk, jl, bn);
    VectorCross(jl, ji, cn);

    // Bond angle ji to jk
    const double cos_theta = VectorDot(ji, jk);
    const double theta = acos(cos_theta);
    // If theta equals 180 degree or 0 degree
    if (IsNearZero(theta) || IsNearZero(fabs(theta - M_PI))) {
      return 0.0;
    }

    const double sin_theta = sin(theta);
    const double sin_dl = VectorDot(an, jl) / sin_theta;

    // the wilson angle:
    const double dl = asin(sin_dl);

    return RAD_TO_DEG * dl;
  }

  double OBForceField::VectorTorsionDerivative(vector3 &i, vector3 &j, vector3 &k, vector3 &l)
  {
    // This is adapted from http://scidok.sulb.uni-saarland.de/volltexte/2007/1325/pdf/Dissertation_1544_Moll_Andr_2007.pdf
    // Many thanks to Andreas Moll and the BALLView developers for this

    // Bond vectors of the three atoms
    vector3 ij, jk, kl;
    // length of the three bonds
    double l_ij, l_jk, l_kl;
    // angle between ijk and jkl:
    double angle_ijk, angle_jkl;

    ij = j - i;
    jk = k - j;
    kl = l - k;

    l_ij = ij.length();
    l_jk = jk.length();
    l_kl = kl.length();

    if (IsNearZero(l_ij) || IsNearZero(l_jk) || IsNearZero(l_kl) ) {
      i = VZero;
      j = VZero;
      k = VZero;
      l = VZero;
      return 0.0;
    }

    angle_ijk = DEG_TO_RAD * vectorAngle(ij, jk);
    angle_jkl = DEG_TO_RAD * vectorAngle(jk, kl);

    // normalize the bond vectors:
    ij /= l_ij;
    jk /= l_jk;
    kl /= l_kl;

    double sin_j = sin(angle_ijk);
    double sin_k = sin(angle_jkl);

    double rsj = l_ij * sin_j;
    double rsk = l_kl * sin_k;

    double rs2j = 1. / (rsj * sin_j);
    double rs2k = 1. / (rsk * sin_k);

    double rrj = l_ij / l_jk;
    double rrk = l_kl / l_jk;

    double rrcj = rrj * (-cos(angle_ijk));
    double rrck = rrk * (-cos(angle_jkl));

    vector3 a = cross(ij, jk);
    vector3 b = cross(jk, kl);
    vector3 c = cross(a, b);
    double d1 = dot(c, jk);
    double d2 = dot(a, b);
    double tor = RAD_TO_DEG * atan2(d1, d2);

    i = -a * rs2j;
    l = b * rs2k;

    j = i * (rrcj - 1.) - l * rrck;
    k = -(i + j + l);

    return tor;
  }

  double OBForceField::VectorTorsionDerivative(double *pos_i, double *pos_j, double *pos_k, double *pos_l,
                                               double *force_i, double *force_j, double *force_k, double *force_l)
  {
    // This is adapted from http://scidok.sulb.uni-saarland.de/volltexte/2007/1325/pdf/Dissertation_1544_Moll_Andr_2007.pdf
    // Many thanks to Andreas Moll and the BALLView developers for this

    // Bond vectors of the three atoms
    double ij[3], jk[3], kl[3];
    VectorSubtract(pos_j, pos_i, ij);
    VectorSubtract(pos_k, pos_j, jk);
    VectorSubtract(pos_l, pos_k, kl);

    // length of the three bonds
    double l_ij, l_jk, l_kl;
    l_ij = VectorLength(ij);
    l_jk = VectorLength(jk);
    l_kl = VectorLength(kl);

    if (IsNearZero(l_ij) || IsNearZero(l_jk) || IsNearZero(l_kl) ) {
      VectorClear(force_i);
      VectorClear(force_j);
      VectorClear(force_k);
      VectorClear(force_l);
      return 0.0;
    }

    // normalize the bond vectors:
    VectorDivide (ij, l_ij, ij);
    VectorDivide (jk, l_jk, jk);
    VectorDivide (kl, l_kl, kl);

    // angle between ijk and jkl:
    double angle_ijk, angle_jkl;

    double cos_ijk = VectorDot(ij, jk);
    if (cos_ijk > 1.0) {
      angle_ijk = 0.0;
      cos_ijk = 1.0;
    } else if (cos_ijk < -1.0) {
      angle_ijk = M_PI;
      cos_ijk = -1.0;
    } else {
      angle_ijk = acos(cos_ijk);
    }

    double cos_jkl = VectorDot(jk, kl);
    if (cos_jkl > 1.0) {
      angle_jkl = 0.0;
      cos_jkl = 1.0;
    } else if (cos_jkl < -1.0) {
      angle_jkl = M_PI;
      cos_jkl = -1.0;
    } else {
      angle_jkl = acos(cos_jkl);
    }

    double sin_j = sin(angle_ijk);
    double sin_k = sin(angle_jkl);

    double rsj = l_ij * sin_j;
    double rsk = l_kl * sin_k;

    double rs2j = 1. / (rsj * sin_j);
    double rs2k = 1. / (rsk * sin_k);

    double rrj = l_ij / l_jk;
    double rrk = l_kl / l_jk;

    double rrcj = rrj * (-cos(angle_ijk));
    double rrck = rrk * (-cos(angle_jkl));

    double a[3];
    VectorCross(ij, jk, a);
    double b[3];
    VectorCross(jk, kl, b);
    double c[3];
    VectorCross(a, b, c);
    double d1 = VectorDot(c, jk);
    double d2 = VectorDot(a, b);
    double tor = RAD_TO_DEG * atan2(d1, d2);

    VectorMultiply(a, -rs2j, force_i);
    VectorMultiply(b, rs2k, force_l);

    VectorMultiply(force_i, (rrcj - 1.0), a);
    VectorMultiply(force_l, rrck, b);
    VectorSubtract(a, b, force_j);

    VectorAdd(force_i, force_j, a);
    VectorAdd(a, force_l, b);
    VectorMultiply(b, -1.0, force_k);

    return tor;
  }

  double OBForceField::VectorTorsion(double *pos_i, double *pos_j, double *pos_k, double *pos_l)
  {
    // Bond vectors of the three atoms
    double ij[3], jk[3], kl[3];
    VectorSubtract(pos_j, pos_i, ij);
    VectorSubtract(pos_k, pos_j, jk);
    VectorSubtract(pos_l, pos_k, kl);

    // length of the three bonds
    const double l_ij = VectorLength(ij);
    const double l_jk = VectorLength(jk);
    const double l_kl = VectorLength(kl);

    if (IsNearZero(l_ij) || IsNearZero(l_jk) || IsNearZero(l_kl) ) {
      return 0.0;
    }

    // normalize the bond vectors:
    VectorDivide (ij, l_ij, ij);
    VectorDivide (jk, l_jk, jk);
    VectorDivide (kl, l_kl, kl);

    double a[3];
    VectorCross(ij, jk, a);
    double b[3];
    VectorCross(jk, kl, b);
    double c[3];
    VectorCross(a, b, c);
    double d1 = VectorDot(c, jk);
    double d2 = VectorDot(a, b);
    double tor = RAD_TO_DEG * atan2(d1, d2);

    return tor;
  }

  bool OBForceField::IsInSameRing(OBAtom* a, OBAtom* b)
  {
    bool a_in, b_in;
    vector<OBRing*> vr;
    vr = _mol.GetSSSR();

    vector<OBRing*>::iterator i;
    vector<int>::iterator j;

    for (i = vr.begin();i != vr.end();++i) {
      a_in = false;
      b_in = false;
      for(j = (*i)->_path.begin();j != (*i)->_path.end();++j) {
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

  OBGridData* OBForceField::GetGrid(double step, double padding, const char* type, double pchg)
  {
    cout << "OBForceFieldMMFF94::GetGrid(" << step << ", " << type << ")" << endl;
    OBFloatGrid fgrid;
    fgrid.Init(_mol, step, padding);
    vector3 min;
    unsigned int xDim, yDim, zDim, xyzDim;

    min = fgrid.GetMin();

    xDim = fgrid.GetXdim();
    yDim = fgrid.GetYdim();
    zDim = fgrid.GetZdim();
    xyzDim = xDim * yDim * zDim;

    cout << "xDim = " << xDim << ", yDim = " << yDim << ", zDim = " << zDim << endl;

    // Add the probe atom
    _mol.BeginModify();
    OBAtom *atom = _mol.NewAtom();
    int index = atom->GetIdx();
    _mol.EndModify();
    SetTypes();
    atom->SetType(type);
    atom->SetPartialCharge(pchg);

    SetupCalculations();

    atom = _mol.GetAtom(index);
    double *pos = atom->GetCoordinate();

    vector3 coord;
    double evdw, eele;
    double distance, minDistance;

    OBGridData *grid = new OBGridData;
    vector3 xAxis, yAxis, zAxis;
    xAxis = vector3(step, 0.0, 0.0);
    yAxis = vector3(0.0, step, 0.0);
    zAxis = vector3(0.0, 0.0, step);

    grid->SetNumberOfPoints(xDim, yDim, zDim);
    grid->SetLimits(min, xAxis, yAxis, zAxis);

    // VDW surface
    for (unsigned int i = 0; i < xDim; ++i) {
      coord.SetX(min[0] + i * step);
      for (unsigned int j = 0; j < yDim; ++j) {
        coord.SetY(min[1] + j * step);
        for (unsigned int k = 0; k < zDim; ++k)
          {
            coord.SetZ(min[2] + k * step);
            minDistance = 1.0E+10;
            FOR_ATOMS_OF_MOL (a, _mol) {
              if (a->GetIdx() == atom->GetIdx())
                continue;
              if (a->GetAtomicNum() == OBElements::Hydrogen)
                continue;

              distance = sqrt(coord.distSq(a->GetVector()));

              if (distance < minDistance)
                minDistance = distance;
            } // end checking atoms
            // negative = away from molecule, 0 = vdw surface, positive = inside
            if (minDistance > 1.0) {
              grid->SetValue(i, j, k, 0.0); // outside the molecule
            } else {
              grid->SetValue(i, j, k, 10e99); // inside the molecule
            }
          } // z-axis
      } // y-axis
    } // x-axis


    unsigned int count = 0;
    for (unsigned int i = 0; i < xDim; ++i) {
      coord.SetX(min[0] + i * step);
      for (unsigned int j = 0; j < yDim; ++j) {
        coord.SetY(min[1] + j * step);
        for (unsigned int k = 0; k < zDim; ++k)
          {
            coord.SetZ(min[2] + k * step);

            count++;
            cout << "\r" << count << "/" << xyzDim;

            if (grid->GetValue(i, j, k) == 0.0) {
              pos[0] = coord.x();
              pos[1] = coord.y();
              pos[2] = coord.z();
              evdw = E_VDW(false);
              eele = E_Electrostatic(false);
              grid->SetValue(i, j, k, evdw + eele);
            }
          } // z-axis
      } // y-axis
    } // x-axis

    cout << endl;

    _mol.BeginModify();
    _mol.DeleteAtom(atom);
    _mol.EndModify();

    return grid;
  }

  /**
   * @example obforcefield_energy.cpp
   * Example showing how to compute the enrgy for a molecule.
   */

} // end namespace OpenBabel


//! \file forcefield.cpp
//! \brief Handle OBForceField class

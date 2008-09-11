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
      - SystematicRotorSearch(): finished
      - RandomRotorSearch(): finished
      - WeightedRotorSearch: finished
      - DistanceGeometry(): needs matrix operations (Eigen)

      Constraints:
      - Fix Atom: working 
      - Fix Atom X: working 
      - Fix Atom Y: working 
      - Fix Atom Z: working 
      - Distance: working
      - Angle: working
      - Torsion: working 
      - Chirality: TODO

      src/forcefields/forcefieldghemical.cpp
      - Atom typing: finished
      - Charges: finished
      - Energy terms: finished
      - Analytical gradients: finished
      - Validation: finished

      src/forcefields/forcefieldmmff94.cpp
      - Atom typing: done.
      - Charges: done.
      - Energy terms: finished (small problems with SSSR 
        algorithm not finding all bridged rings)
      - Analytical gradients: finished 
      - Validation: http://home.scarlet.be/timvdm/MMFF94_validation_output.gz

      src/forcefields/forcefielduff.cpp
      - Energy terms: finished
      - OOP: needs validation
      - Gradients: need OOP gradient
      - Validation in progress...


      forcefield.cpp layout:

        1) Constraints                XXXX01
        2) Methods for logging        XXXX02
        3) General                    XXXX03
        4) Structure generation       XXXX04
        5) Cut-off                    XXXX05
        6) Interaction groups         XXXX06
        8) Molecular Dynamics         XXXX08
        9) Misc.                      XXXX09

***********************************************************************/

#include <openbabel/obcutoffgrid.h>

#include <openbabel/obforcefield.h>
#include <openbabel/ffinternal.h>

#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obiter.h>
#include <openbabel/math/vector.h>
#include <openbabel/math/matrix3x3.h>
#include <openbabel/rotamer.h>
#include <openbabel/rotor.h>

#include <openbabel/obminimize.h>

#include <openbabel/babelconfig.h>

using namespace std;

namespace OpenBabel
{
  class OBForceFieldPrivate 
  {
    public:
    // general variables
    OBMol 	mol; //!< Molecule to be evaluated or minimized
    bool 	init; //!< Used to make sure we only parse the parameter file once, when needed
    std::string	parFile; //! < parameter file name
    bool 	validSetup; //!< was the last call to Setup succesfull
    //double	*gradientPtr; //!< pointer to the gradients (used by AddGradient(), minimization functions, ...)
    //std::vector<Eigen::Vector3d> positions;
    //std::vector<Eigen::Vector3d> gradients;
    // logging variables
    std::ostream* logos; //!< Output for logfile
    char 	logbuf[BUFF_SIZE]; //!< Temporary buffer for logfile output
    //int 	loglvl; //!< Log level for output
    int 	origLogLevel;
    // conformer genereation (rotor search) variables
    int 	current_conformer; //!< used to hold i for current conformer (needed by UpdateConformers)
    std::vector<double> energies; //!< used to hold the energies for all conformers
    // minimization variables
    double 	econv, e_n1; //!< Used for conjugate gradients and steepest descent(Initialize and TakeNSteps)
    int 	cstep, nsteps; //!< Used for conjugate gradients and steepest descent(Initialize and TakeNSteps)
    std::vector<Eigen::Vector3d> grad1; //!< Used for conjugate gradients and steepest descent(Initialize and TakeNSteps)
    unsigned int ncoords; //!< Number of coordinates for conjugate gradients
    int         linesearch; //!< LineSearch type
    // molecular dynamics variables
    double 	timestep; //!< Molecular dynamics time step in picoseconds
    double 	temp; //!< Molecular dynamics temperature in Kelvin
    double 	*velocityPtr; //!< pointer to the velocities
    // contraint varibles
    static OBFFConstraints constraints; //!< Constraints
    static unsigned int fixAtom; //!< SetFixAtom()/UnsetFixAtom()
    static unsigned int ignoreAtom; //!< SetIgnoreAtom()/UnsetIgnoreAtom()
    // cut-off variables
    //bool 	cutoff; //!< true = cut-off enabled
    double 	rvdw; //!< VDW cut-off distance
    double 	rele; //!< Electrostatic cut-off distance
    OBBitVec	vdwpairs; //!< VDW pairs that should be calculated
    OBBitVec	elepairs; //!< Electrostatic pairs that should be calculated
    int 	pairfreq; //!< The frequence to update non-bonded pairs
    // group variables
    std::vector<OBBitVec> intraGroup; //!< groups for which intra-molecular interactions should be calculated
    std::vector<OBBitVec> interGroup; //!< groups for which intra-molecular interactions should be calculated
    std::vector<std::pair<OBBitVec, OBBitVec> > interGroups; //!< groups for which intra-molecular 
  };

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
      pFF->SteepestDescent(1000);
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
        if (pocket.BitIsOn(a->GetIdx())
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
  
  void OBForceField::SetParameterFile(const std::string &filename)
  {
    d->parFile = filename;
    d->init = false;
  }
 
  std::string& OBForceField::GetParameterFile()
  {
    return d->parFile;
  }
 
  //////////////////////////////////////////////////////////////////////////////////
  //
  // Constraints  (XXXX01)
  //
  //////////////////////////////////////////////////////////////////////////////////
 
  OBFFConstraints OBForceFieldPrivate::constraints = OBFFConstraints(); // define static data variable
  unsigned int OBForceFieldPrivate::fixAtom = 0; // define static data variable
  unsigned int OBForceFieldPrivate::ignoreAtom = 0; // define static data variable

  OBFFConstraints& OBForceField::GetConstraints() 
  { 
    return d->constraints; 
  }
  
  void OBForceField::SetConstraints(OBFFConstraints& constraints) 
  { 
    if (!(d->constraints.GetIgnoredBitVec() == constraints.GetIgnoredBitVec())) {
      d->constraints = constraints;
      if (!SetupCalculations()) {
        d->validSetup = false;
        return;
      }
    } else {
      d->constraints = constraints;
    }

    d->constraints.Setup(d->mol); 
  }
 
  void OBForceField::SetFixAtom(int index) 
  {
    d->fixAtom = index; 
  }

  void OBForceField::UnsetFixAtom()
  {
    d->fixAtom = 0;
  }

  void OBForceField::SetIgnoreAtom(int index) 
  {
    d->ignoreAtom = index; // remember the index
  }

  void OBForceField::UnsetIgnoreAtom()
  {
    d->ignoreAtom = 0;
  }

  bool OBForceField::IgnoreCalculation(unsigned int a, unsigned int b)
  {
    if (!OBForceFieldPrivate::ignoreAtom)
      return false;

    if (OBForceFieldPrivate::ignoreAtom == a)
      return true;
    if (OBForceFieldPrivate::ignoreAtom == b)
      return true;

    return false;
  }

  bool OBForceField::IgnoreCalculation(unsigned int a, unsigned int b, unsigned int c)
  {
    if (OBForceField::IgnoreCalculation(a, c))
      return true;
    if (OBForceFieldPrivate::ignoreAtom == b)
      return true;
    
    return false;
  }

  bool OBForceField::IgnoreCalculation(unsigned int a, unsigned int b, unsigned int c, unsigned int d)
  {
    if (OBForceField::IgnoreCalculation(a, b, c))
      return true;
    if (OBForceFieldPrivate::ignoreAtom == d)
      return true;
    
    return false;
  }

  OBFFConstraints::OBFFConstraints()
  {
    m_factor = 50000.0;
  }
     
  OBFFConstraints::~OBFFConstraints()
  {
     m_constraints.clear();
     m_ignored.Clear();
     m_fixed.Clear();
     m_Xfixed.Clear();
     m_Yfixed.Clear();
     m_Zfixed.Clear();
  }
    
  OBFFConstraints& OBFFConstraints::operator=(const OBFFConstraints &ai) 
  {
    if (this != &ai) {
      m_constraints = ai.m_constraints;
      m_ignored = ai.m_ignored;
      m_fixed = ai.m_fixed;
      m_Xfixed = ai.m_Xfixed;
      m_Yfixed = ai.m_Yfixed;
      m_Zfixed = ai.m_Zfixed;
    }
    return *this;
  }
 
  void OBFFConstraints::Setup(OBMol& mol)
  {
    vector<OBFFConstraint>::iterator i;
    
    for (i = m_constraints.begin(); i != m_constraints.end(); ++i) {
      i->a = mol.GetAtom(i->ia);
      i->b = mol.GetAtom(i->ib);
      i->c = mol.GetAtom(i->ic);
      i->d = mol.GetAtom(i->id);
    }
  }
  
  /*
  Eigen::Vector3d OBFFConstraints::GetGradient(int a)
  {
    vector<OBFFConstraint>::iterator i;

    Eigen::Vector3d grad(0.0, 0.0, 0.0);
    
    for (i = m_constraints.begin(); i != m_constraints.end(); ++i)
      grad += i->GetGradient(a);

    return grad;
  }
  */
  
  void OBFFConstraints::Clear()
  {
    m_constraints.clear();
    m_ignored.Clear();
    m_fixed.Clear();
    m_Xfixed.Clear();
    m_Yfixed.Clear();
    m_Zfixed.Clear();
  }
  
  int OBFFConstraints::Size() const
  {
    return m_constraints.size();
  }
  
  double OBFFConstraints::GetConstraintEnergy()
  {
    vector<OBFFConstraint>::iterator i;
    double constraint_energy = 0.0;
       
    for (i = m_constraints.begin(); i != m_constraints.end(); ++i)
      if ( (i->type == OBFF_CONST_DISTANCE) || (i->type == OBFF_CONST_ANGLE) || 
           (i->type == OBFF_CONST_TORSION) ) {
        Eigen::Vector3d Fa, Fb, Fc, Fd;
        double delta, delta2, rab, theta, dE;
	
        switch (i->type) {
        case OBFF_CONST_DISTANCE:
          if ((i->a == NULL) || ((i->b) == NULL))
            break;

          rab = OpenBabel::VectorDistanceDerivative((i->a)->GetVector(), (i->b)->GetVector(), i->Fa, i->Fb);
          delta = rab - i->constraint_value;
          delta2 = delta * delta;
          constraint_energy += i->factor * delta2;
          dE = 2.0 * i->factor * delta;

          i->Fa *= dE;
          i->Fb *= dE; 
          break;
        case OBFF_CONST_ANGLE:
          if ((i->a == NULL) || ((i->b) == NULL) || ((i->c) == NULL))
            break;

          theta = OpenBabel::VectorAngleDerivative((i->a)->GetVector(), (i->b)->GetVector(), 
              (i->c)->GetVector(), i->Fa, i->Fb, i->Fc);
          delta = theta - i->constraint_value;
          delta2 = delta * delta;
          constraint_energy += 0.0002 * i->factor * delta2;
          dE = 0.0004 * i->factor * delta;

          i->Fa *= dE;
          i->Fb *= dE;
          i->Fc *= dE;
          break;
        case OBFF_CONST_TORSION:
          if ((i->a == NULL) || ((i->b) == NULL) || ((i->c) == NULL) || ((i->d) == NULL))
            break;

          theta = OpenBabel::VectorTorsionDerivative((i->a)->GetVector(), (i->b)->GetVector(),
              (i->c)->GetVector(), (i->d)->GetVector(), i->Fa, i->Fb, i->Fc, i->Fd);
          if (!isfinite(theta))
            theta = 1.0e-7;
            
          theta = DEG_TO_RAD * (theta + 180.0 - i->constraint_value);
          constraint_energy +=  0.001 * i->factor * (1.0 + cos(theta));
          dE = 0.001 * i->factor * sin(theta);
          
          i->Fa *= dE;
          i->Fb *= dE;
          i->Fc *= dE;
          i->Fd *= dE;
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

    for (i = m_constraints.begin(); i != m_constraints.end(); ++n, ++i) {
      if (n == index) {
        if (i->type == OBFF_CONST_IGNORE)
          m_ignored.SetBitOff(i->ia);
        if (i->type == OBFF_CONST_ATOM)
          m_fixed.SetBitOff(i->ia);
        if (i->type == OBFF_CONST_ATOM_X)
          m_Xfixed.SetBitOff(i->ia);
        if (i->type == OBFF_CONST_ATOM_Y)
          m_Yfixed.SetBitOff(i->ia);
        if (i->type == OBFF_CONST_ATOM_Z)
          m_Zfixed.SetBitOff(i->ia);

 
        m_constraints.erase(i);
        break;
      }
    }
  }
  
  void OBFFConstraints::SetFactor(double factor)
  {
    m_factor = factor;
  }

  double OBFFConstraints::GetFactor()
  {
    return m_factor;
  }

  void OBFFConstraints::AddIgnore(int a)
  {
    m_ignored.SetBitOn(a);

    OBFFConstraint constraint;
    constraint.type = OBFF_CONST_IGNORE; // constraint type
    constraint.ia   = a; // atom to fix
    m_constraints.push_back(constraint);
  }
  
  void OBFFConstraints::AddAtomConstraint(int a)
  {
    m_fixed.SetBitOn(a);

    OBFFConstraint constraint;
    constraint.type = OBFF_CONST_ATOM; // constraint type
    constraint.ia   = a; // atom to fix
    constraint.factor = m_factor;
    m_constraints.push_back(constraint);
  }
  
  void OBFFConstraints::AddAtomXConstraint(int a)
  {
    m_Xfixed.SetBitOn(a);

    OBFFConstraint constraint;
    constraint.type = OBFF_CONST_ATOM_X; // constraint type
    constraint.ia   = a; // atom to fix
    constraint.factor = m_factor;
    m_constraints.push_back(constraint);
  }
  
  void OBFFConstraints::AddAtomYConstraint(int a)
  {
    m_Yfixed.SetBitOn(a);

    OBFFConstraint constraint;
    constraint.type = OBFF_CONST_ATOM_Y; // constraint type
    constraint.ia   = a; // atom to fix
    constraint.factor = m_factor;
    m_constraints.push_back(constraint);
  } 
  
  void OBFFConstraints::AddAtomZConstraint(int a)
  {
    m_Zfixed.SetBitOn(a);

    OBFFConstraint constraint;
    constraint.type = OBFF_CONST_ATOM_Z; // constraint type
    constraint.ia   = a; // atom to fix
    constraint.factor = m_factor;
    m_constraints.push_back(constraint);
  }
  
  void OBFFConstraints::AddDistanceConstraint(int a, int b, double length)
  {
    OBFFConstraint constraint;
    constraint.type = OBFF_CONST_DISTANCE; // constraint type
    constraint.ia   = a; // atom a
    constraint.ib   = b; // atom b
    constraint.constraint_value = length; // bond length
    constraint.factor = m_factor;
    m_constraints.push_back(constraint);
  }
  
  void OBFFConstraints::AddAngleConstraint(int a, int b, int c, double angle)
  {
    OBFFConstraint constraint;
    constraint.type = OBFF_CONST_ANGLE; // constraint type
    constraint.ia   = a; // atom a
    constraint.ib   = b; // atom b
    constraint.ic   = c; // atom c
    constraint.constraint_value = angle; // angle
    constraint.factor = m_factor;
    m_constraints.push_back(constraint);
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
    constraint.factor = m_factor;
    m_constraints.push_back(constraint);
  }
  
  int OBFFConstraints::GetConstraintType(unsigned int index) const
  {
    if (index >= m_constraints.size())
      return 0;

    return m_constraints[index].type; 
  }
  
  double OBFFConstraints::GetConstraintValue(unsigned int index) const 
  {
    if (index >= m_constraints.size())
      return 0;

    return m_constraints[index].constraint_value; 
  }
  
  int OBFFConstraints::GetConstraintAtomA(unsigned int index) const
  {
    if (index >= m_constraints.size())
      return 0;

    return m_constraints[index].ia; 
  }
  
  int OBFFConstraints::GetConstraintAtomB(unsigned int index) const
  {
    if (index >= m_constraints.size())
      return 0;

    return m_constraints[index].ib; 
  }
  
  int OBFFConstraints::GetConstraintAtomC(unsigned int index) const
  {
    if (index >= m_constraints.size())
      return 0;

    return m_constraints[index].ic; 
  }
  
  int OBFFConstraints::GetConstraintAtomD(unsigned int index) const
  {
    if (index >= m_constraints.size())
      return 0;

    return m_constraints[index].id; 
  }
  
  bool OBFFConstraints::IsIgnored(unsigned int index)
  {
    return m_ignored.BitIsSet(index);
  }
  
  bool OBFFConstraints::IsFixed(unsigned int index)
  {
    return m_fixed.BitIsSet(index);
  }
  
  bool OBFFConstraints::IsXFixed(unsigned int index)
  {
    return m_Xfixed.BitIsSet(index);
  }
  
  bool OBFFConstraints::IsYFixed(unsigned int index)
  {
    return m_Yfixed.BitIsSet(index);
  }
  
  bool OBFFConstraints::IsZFixed(unsigned int index)
  {
    return m_Zfixed.BitIsSet(index);
  }
  
  //////////////////////////////////////////////////////////////////////////////////
  //
  // Methods for logging  (XXXX02)
  //
  //////////////////////////////////////////////////////////////////////////////////
  
  void OBForceField::PrintTypes()
  {
    IF_OBFF_LOGLVL_LOW {
      OBFFLog("\nA T O M   T Y P E S\n\n");
      OBFFLog("IDX\tTYPE\n");
      
      FOR_ATOMS_OF_MOL (a, d->mol) {
        snprintf(d->logbuf, BUFF_SIZE, "%d\t%s\n", a->GetIdx(), a->GetType());
        OBFFLog(d->logbuf);
      }
    }
  }
    
  void OBForceField::PrintFormalCharges()
  {
    IF_OBFF_LOGLVL_LOW {
      OBFFLog("\nF O R M A L   C H A R G E S\n\n");
      OBFFLog("IDX\tCHARGE\n");
      
      FOR_ATOMS_OF_MOL (a, d->mol) {
        snprintf(d->logbuf, BUFF_SIZE, "%d\t%f\n", a->GetIdx(), a->GetPartialCharge());
        OBFFLog(d->logbuf);
      }
    }
  }

  void OBForceField::PrintPartialCharges()
  {
    IF_OBFF_LOGLVL_LOW {
      OBFFLog("\nP A R T I A L   C H A R G E S\n\n");
      OBFFLog("IDX\tCHARGE\n");
      
      FOR_ATOMS_OF_MOL (a, d->mol) {
        snprintf(d->logbuf, BUFF_SIZE, "%d\t%f\n", a->GetIdx(), a->GetPartialCharge());
        OBFFLog(d->logbuf);
      }
    }
  }
  
  void OBForceField::PrintVelocities()
  {
    IF_OBFF_LOGLVL_LOW {
      OBFFLog("\nA T O M   V E L O C I T I E S\n\n");
      OBFFLog("IDX\tVELOCITY\n");
      
      FOR_ATOMS_OF_MOL (a, d->mol) {
        snprintf(d->logbuf, BUFF_SIZE, "%d\t<%8.3f, %8.3f, %8.3f>\n", a->GetIdx(), d->velocityPtr[a->GetIdx()], 
            d->velocityPtr[a->GetIdx()+1], d->velocityPtr[a->GetIdx()+2]);
        OBFFLog(d->logbuf);
      }
    }
  }
 
  bool OBForceField::SetLogFile(std::ostream* pos)
  {
    if(pos)
      d->logos = pos;
    else
      d->logos = &cout;
    
    return true;
  }
    
  void OBForceField::OBFFLog(std::string msg)
  {
    if (!d->logos)
      return;
     
    *(d->logos) << msg;
  }

  void OBForceField::OBFFLog(const char *msg)
  {
    if (!d->logos)
      return;
    
    *(d->logos) << msg;
  }
 
 
  //////////////////////////////////////////////////////////////////////////////////
  //
  // General  (XXXX03)
  //
  //////////////////////////////////////////////////////////////////////////////////
 
  bool OBForceField::Setup(OBMol &mol)
  {
    if (!d->init) {
      ParseParamFile();
      d->init = true;
      d->velocityPtr = NULL;       
    }    

    if (IsSetupNeeded(mol)) {
      d->mol = mol;
      d->ncoords = mol.NumAtoms() * 3;
      
      GetPositions().clear();
      GetPositions().reserve(mol.NumAtoms());
      FOR_ATOMS_OF_MOL (atom, mol)
        GetPositions().push_back(atom->GetVector());

      if (d->velocityPtr)
        delete [] d->velocityPtr;
      d->velocityPtr = NULL;       

      GetGradients().resize(d->mol.NumAtoms());
      

      if (d->mol.NumAtoms() && d->constraints.Size())
        d->constraints.Setup(d->mol);

      d->mol.UnsetSSSRPerceived();
      d->mol.DeleteData(OBGenericDataType::TorsionData); // bug #1954233
      
      if (!SetTypes()) {
        d->validSetup = false;
        return false;
      }

      SetFormalCharges();
      SetPartialCharges();

      if (!SetupCalculations()) {
        d->validSetup = false;
        return false;
      }

    } else {
      if (d->validSetup) {
        SetCoordinates(mol);
        return true;
      } else {
        return false;
      }
    }

    d->validSetup = true;
    return true;
  }
  
  bool OBForceField::Setup(OBMol &mol, OBFFConstraints &constraints)
  {
    if (!d->init) {
      ParseParamFile();
      d->init = true;
      d->velocityPtr = NULL;       
    }    

    if (IsSetupNeeded(mol)) {
      d->mol = mol;
      d->ncoords = d->mol.NumAtoms() * 3;

      if (d->velocityPtr)
        delete [] d->velocityPtr;
      d->velocityPtr = NULL;       

      GetPositions().clear();
      GetPositions().reserve(mol.NumAtoms());
      FOR_ATOMS_OF_MOL (atom, mol)
        GetPositions().push_back(atom->GetVector());

      GetGradients().resize(d->mol.NumAtoms());

      d->constraints = constraints;
      if (d->mol.NumAtoms() && d->constraints.Size())
        d->constraints.Setup(d->mol);
      
      d->mol.UnsetSSSRPerceived();
      d->mol.DeleteData(OBGenericDataType::TorsionData); // bug #1954233

      if (!SetTypes()) {
        d->validSetup = false;
        return false;
      }

      SetFormalCharges();
      SetPartialCharges();

      if (!SetupCalculations()) {
        d->validSetup = false;
        return false;
      }

    } else {
      if (d->validSetup) {
        if (!(d->constraints.GetIgnoredBitVec() == constraints.GetIgnoredBitVec())) {
          d->constraints = constraints;
          if (!SetupCalculations()) {
            d->validSetup = false;
            return false;
          }
        } else {
          d->constraints = constraints;
        }

        d->constraints.Setup(d->mol);
        SetCoordinates(mol);
        return true;
      } else {
        return false;
      }
    }

    d->validSetup = true;
    return true;
  }
 
  bool OBForceField::SetLogLevel(int level)
  {
    m_loglvl = level;
    return true;
  }
  
  // This function might need expanding, bu could really increase performance for
  // avogadro's AutoOpt tool
  bool OBForceField::IsSetupNeeded(OBMol &mol)
  {
    if (d->mol.NumAtoms() != mol.NumAtoms())
      return true;
      
    if (d->mol.NumBonds() != mol.NumBonds())
      return true;

    FOR_ATOMS_OF_MOL (atom, d->mol) {
      if (atom->GetAtomicNum() != (mol.GetAtom(atom->GetIdx()))->GetAtomicNum())
        return true;
      if (atom->GetValence() != (mol.GetAtom(atom->GetIdx()))->GetValence())
        return true;
      if (atom->BOSum() != (mol.GetAtom(atom->GetIdx()))->BOSum())
        return true;
    }

    return false;
  }
  
  OBMol* OBForceField::GetMolecule()
  {
    return &(d->mol);
  }
  
  bool OBForceField::GetAtomTypes(OBMol &mol)
  {
    if (d->mol.NumAtoms() != mol.NumAtoms())
      return false;
    
    FOR_ATOMS_OF_MOL (atom, d->mol) {
      if (atom->HasData("FFAtomType")) {
        OBPairData *data = (OBPairData*) atom->GetData("FFAtomType");
	data->SetValue(atom->GetType());
      } else {
        OBPairData *data = new OBPairData();
       	data->SetAttribute("FFAtomType");
	data->SetValue(atom->GetType());
        atom->SetData(data);
      }
    }

    return true; 
  }

  bool OBForceField::GetPartialCharges(OBMol &mol)
  {
    if (d->mol.NumAtoms() != mol.NumAtoms())
      return false;
    
    ostringstream out;
    FOR_ATOMS_OF_MOL (atom, d->mol) {
      out << atom->GetPartialCharge();
      if (atom->HasData("FFPartialCharge")) {
        OBPairData *data = (OBPairData*) atom->GetData("FFPartialCharge");
	data->SetValue(out.str());
      } else {
        OBPairData *data = new OBPairData();
       	data->SetAttribute("FFPartialCharge");
	data->SetValue(out.str());
        atom->SetData(data);
      }
    }

    return true; 
  }


  bool OBForceField::GetCoordinates(OBMol &mol)
  { 
    OBAtom *atom;
    
    if (d->mol.NumAtoms() != mol.NumAtoms())
      return false;
    
    // Copy coordinates for current conformer only
    FOR_ATOMS_OF_MOL (a, d->mol) {
      atom = mol.GetAtom(a->GetIdx());
      atom->SetVector(a->GetVector());
    }

    if (!mol.HasData(OBGenericDataType::ConformerData))
      mol.SetData(new OBConformerData);
    OBConformerData *cd = (OBConformerData*) mol.GetData(OBGenericDataType::ConformerData);
    cd->SetEnergies(d->energies);
    
    vector<vector<Eigen::Vector3d> > confForces;
    /*
    vector<Eigen::Vector3d> forces;
    for (unsigned int i = 0; i < d->mol.NumAtoms(); ++i) {
      const int coordIdx = i * 3;
      forces.push_back(Eigen::Vector3d(d->gradientPtr[coordIdx], 
          d->gradientPtr[coordIdx+1], d->gradientPtr[coordIdx+2]));
    }
    confForces.push_back(forces);
    */
    confForces.push_back(GetGradients());
    cd->SetForces(confForces);

    return true;
  }
  
  bool OBForceField::GetConformers(OBMol &mol)
  { 
    if (d->mol.NumAtoms() != mol.NumAtoms())
      return false;
    
    //Copy conformer information
    if (d->mol.NumConformers() > 1) {
      int k,l;
      vector<double*> conf;
      double* xyz = NULL;
      for (k=0 ; k<d->mol.NumConformers() ; ++k) {
        xyz = new double [3*d->mol.NumAtoms()];
        for (l=0 ; l<(int) (3*d->mol.NumAtoms()) ; ++l)
          xyz[l] = d->mol.GetConformer(k)[l];
        conf.push_back(xyz);
      }
      mol.SetConformers(conf);
      mol.SetConformer(d->current_conformer);
      
      if (!mol.HasData(OBGenericDataType::ConformerData))
        mol.SetData(new OBConformerData);
      OBConformerData *cd = (OBConformerData*) mol.GetData(OBGenericDataType::ConformerData);
      cd->SetEnergies(d->energies);
    }
    
    return true;
  }
  
  bool OBForceField::SetCoordinates(OBMol &mol)
  { 
    OBAtom *atom;
    
    if (d->mol.NumAtoms() != mol.NumAtoms())
      return false;
    
    // Copy coordinates for current conformer only
    FOR_ATOMS_OF_MOL (a, mol) {
      atom = d->mol.GetAtom(a->GetIdx());
      atom->SetVector(a->GetVector());
    }

    return true;
  }
 
  bool OBForceField::SetConformers(OBMol &mol)
  { 
    OBAtom *atom;

    if (d->mol.NumAtoms() != mol.NumAtoms())
      return false;
    
    FOR_ATOMS_OF_MOL (a, mol) {
      atom = d->mol.GetAtom(a->GetIdx());
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
      d->mol.SetConformers(conf);
      d->mol.SetConformer(d->current_conformer);
    }
    
    return true;
  }
  
  Eigen::Vector3d OBForceField::ValidateGradientError(Eigen::Vector3d &numgrad, Eigen::Vector3d &anagrad)
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

    return Eigen::Vector3d(errx, erry, errz);
  }

  void OBForceField::ClearGradients() 
  { 
    // We cannot use memset because IEEE floating point representations
    // are not guaranteed by C/C++ standard, but this loop can be
    // unrolled or vectorized by compilers
    for (unsigned int i = 0; i < GetGradients().size(); ++i)
      GetGradients()[i] = Eigen::Vector3d::Zero();
  }

  //////////////////////////////////////////////////////////////////////////////////
  //
  // Structure generation  (XXXX04)
  //
  //////////////////////////////////////////////////////////////////////////////////

  int OBForceField::SystematicRotorSearchInitialize(unsigned int geomSteps)
  {
    OBRotorList rl;
    OBRotamerList rotamers;
    OBRotorIterator ri;
    OBRotor *rotor;

    d->origLogLevel = m_loglvl;

    rl.Setup(d->mol);
    rotamers.Setup(d->mol, rl);
    
    IF_OBFF_LOGLVL_LOW {
      OBFFLog("\nS Y S T E M A T I C   R O T O R   S E A R C H\n\n");
      snprintf(d->logbuf, BUFF_SIZE, "  NUMBER OF ROTATABLE BONDS: %d\n", rl.Size());
      OBFFLog(d->logbuf);

      unsigned long int combinations = 1;
      for (rotor = rl.BeginRotor(ri); rotor; rotor = rl.NextRotor(ri)) {
        combinations *= rotor->GetTorsionValues().size();
      }
      snprintf(d->logbuf, BUFF_SIZE, "  NUMBER OF POSSIBLE ROTAMERS: %lu\n", combinations);
      OBFFLog(d->logbuf);
    }

    if (!rl.Size()) { // only one conformer
      IF_OBFF_LOGLVL_LOW
        OBFFLog("  GENERATED ONLY ONE CONFORMER\n\n");
 
      OBMinimize mini(this);
      mini.ConjugateGradients(geomSteps); // final energy minimizatin for best conformation
      
      return 1; // there are no more conformers
    }

    OBRotamerKeys keys;
    rotor = rl.BeginRotor(ri);
    for (unsigned int i = 1; i < rl.Size() + 1; ++i, rotor = rl.NextRotor(ri)) // foreach rotor
      keys.AddRotor(rotor->GetTorsionValues().size());
    
    while (keys.Next()) 
      rotamers.AddRotamer(keys.GetKey());

    rotamers.ExpandConformerList(d->mol, d->mol.GetConformers());
      
    IF_OBFF_LOGLVL_LOW {
      snprintf(d->logbuf, BUFF_SIZE, "  GENERATED %d CONFORMERS\n\n", d->mol.NumConformers());
      OBFFLog(d->logbuf);
      OBFFLog("CONFORMER     ENERGY\n");
      OBFFLog("--------------------\n");
    }

    d->current_conformer = 0;
    d->energies.clear();
      
    return d->mol.NumConformers();
  }
  
  bool OBForceField::SystematicRotorSearchNextConformer(unsigned int geomSteps) 
  {
    if (d->current_conformer >=  d->mol.NumConformers()) { // done
      // Select conformer with lowest energy
      int best_conformer = 0;
      for (int i = 0; i < d->mol.NumConformers(); i++) {
        if (d->energies[i] < d->energies[best_conformer])
          best_conformer = i;
      }
  
      IF_OBFF_LOGLVL_LOW {
        snprintf(d->logbuf, BUFF_SIZE, "\n  CONFORMER %d HAS THE LOWEST ENERGY\n\n",  best_conformer + 1);
        OBFFLog(d->logbuf);
      }

      d->mol.SetConformer(best_conformer);
      d->current_conformer = best_conformer;
    
      return false;
    }

    d->mol.SetConformer(d->current_conformer); // select conformer
    
    m_loglvl = OBFF_LOGLVL_NONE;
    OBMinimize mini(this);
    mini.ConjugateGradients(geomSteps); // energy minimization for conformer
    m_loglvl = d->origLogLevel;
    
    d->energies.push_back(Energy(false)); // calculate and store energy
      
    IF_OBFF_LOGLVL_LOW {
      snprintf(d->logbuf, BUFF_SIZE, "   %3d   %20.3f\n", (d->current_conformer + 1), d->energies[d->current_conformer]);
      OBFFLog(d->logbuf);
    }
    
    d->current_conformer++;
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
    d->origLogLevel = m_loglvl;

    if (d->mol.GetCoordinates() == NULL)
      return;

    rl.Setup(d->mol);
    rotamers.Setup(d->mol, rl);
    
    IF_OBFF_LOGLVL_LOW {
      OBFFLog("\nR A N D O M   R O T O R   S E A R C H\n\n");
      snprintf(d->logbuf, BUFF_SIZE, "  NUMBER OF ROTATABLE BONDS: %d\n", rl.Size());
      OBFFLog(d->logbuf);

      unsigned long int combinations = 1;
      for (rotor = rl.BeginRotor(ri); rotor;
           rotor = rl.NextRotor(ri)) {
        combinations *= rotor->GetTorsionValues().size();
      }
      snprintf(d->logbuf, BUFF_SIZE, "  NUMBER OF POSSIBLE ROTAMERS: %lu\n", combinations);
      OBFFLog(d->logbuf);
    }

    if (!rl.Size()) { // only one conformer
      IF_OBFF_LOGLVL_LOW
        OBFFLog("  GENERATED ONLY ONE CONFORMER\n\n");
 
      m_loglvl = OBFF_LOGLVL_NONE;
      OBMinimize mini(this);
      mini.ConjugateGradients(geomSteps); // energy minimization for conformer
      m_loglvl = d->origLogLevel;

      return;
    }

    std::vector<int> rotorKey(rl.Size() + 1, 0); // indexed from 1

    for (unsigned int c = 0; c < conformers; ++c) {
      rotor = rl.BeginRotor(ri);
      for (unsigned int i = 1; i < rl.Size() + 1; ++i, rotor = rl.NextRotor(ri)) {
        // foreach rotor
        rotorKey[i] = generator.NextInt() % rotor->GetTorsionValues().size();
      }
      rotamers.AddRotamer(rotorKey);
    }

    rotamers.ExpandConformerList(d->mol, d->mol.GetConformers());
      
    IF_OBFF_LOGLVL_LOW {
      snprintf(d->logbuf, BUFF_SIZE, "  GENERATED %d CONFORMERS\n\n", d->mol.NumConformers());
      OBFFLog(d->logbuf);
      OBFFLog("CONFORMER     ENERGY\n");
      OBFFLog("--------------------\n");
    }
  
    d->current_conformer = 0;
    d->energies.clear();
  }

  bool OBForceField::RandomRotorSearchNextConformer(unsigned int geomSteps) 
  {
    if (d->current_conformer >=  d->mol.NumConformers()) { // done
      // Select conformer with lowest energy
      int best_conformer = 0;
      for (int i = 0; i < d->mol.NumConformers(); i++) {
        if (d->energies[i] < d->energies[best_conformer])
          best_conformer = i;
      }
  
      IF_OBFF_LOGLVL_LOW {
        snprintf(d->logbuf, BUFF_SIZE, "\n  CONFORMER %d HAS THE LOWEST ENERGY\n\n",  best_conformer + 1);
        OBFFLog(d->logbuf);
      }

      d->mol.SetConformer(best_conformer);
      d->current_conformer = best_conformer;
      
      return false;
    }
    
    d->mol.SetConformer(d->current_conformer); // select conformer

    m_loglvl = OBFF_LOGLVL_NONE;
    OBMinimize mini(this);
    mini.ConjugateGradients(geomSteps); // energy minimization for conformer
    m_loglvl = d->origLogLevel;
    
    d->energies.push_back(Energy(false)); // calculate and store energy
      
    IF_OBFF_LOGLVL_LOW {
      snprintf(d->logbuf, BUFF_SIZE, "   %3d      %8.3f\n", (d->current_conformer + 1), d->energies[d->current_conformer]);
      OBFFLog(d->logbuf);
    }
    
    d->current_conformer++;
    return true;
  } 
    
  void OBForceField::RandomRotorSearch(unsigned int conformers, unsigned int geomSteps) 
  {
    RandomRotorSearchInitialize(conformers, geomSteps);
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
        if ((int)j == rotorKey[i])
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
        if ((int)j == rotorKey[i])
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
    int origLogLevel = m_loglvl;

    if (d->mol.GetCoordinates() == NULL)
      return;

    rl.Setup(d->mol);
    rotamers.Setup(d->mol, rl);
    
    IF_OBFF_LOGLVL_LOW {
      OBFFLog("\nW E I G H T E D   R O T O R   S E A R C H\n\n");
      snprintf(d->logbuf, BUFF_SIZE, "  NUMBER OF ROTATABLE BONDS: %d\n", rl.Size());
      OBFFLog(d->logbuf);

      unsigned long int combinations = 1;
      for (rotor = rl.BeginRotor(ri); rotor;
           rotor = rl.NextRotor(ri)) {
        combinations *= rotor->GetTorsionValues().size();
      }
      snprintf(d->logbuf, BUFF_SIZE, "  NUMBER OF POSSIBLE ROTAMERS: %lu\n", combinations);
      OBFFLog(d->logbuf);
    }

    if (!rl.Size()) { // only one conformer
      IF_OBFF_LOGLVL_LOW
        OBFFLog("  GENERATED ONLY ONE CONFORMER\n\n");
 
      m_loglvl = OBFF_LOGLVL_NONE;
      OBMinimize mini(this);
      mini.ConjugateGradients(geomSteps); // energy minimization for conformer
      m_loglvl = origLogLevel;

      return;
    }

    double *initialCoord = new double [d->mol.NumAtoms() * 3]; // initial state
    memcpy((char*)initialCoord,(char*)d->mol.GetCoordinates(),sizeof(double)*3*d->mol.NumAtoms());

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
    // int settings;
    for (unsigned int i = 1; i < rl.Size() + 1; ++i, rotor = rl.NextRotor(ri)) {
      // foreach rotor
      energies.clear();
      for (unsigned int j = 0; j < rotor->GetTorsionValues().size(); j++) {
        // foreach rotor position
        d->mol.SetCoordinates(initialCoord);
        rotorKey[i] = j;
        rotamers.SetCurrentCoordinates(d->mol, rotorKey);
        
        m_loglvl = OBFF_LOGLVL_NONE;
        OBMinimize mini(this);
        mini.SteepestDescent(geomSteps); // energy minimization for conformer
        m_loglvl = origLogLevel;
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
      for (unsigned int j = 0; j < rotor->GetTorsionValues().size(); j++) {
        currentE = energies[j];
        // add back the Boltzmann population for these relative energies at 300K (assuming kJ/mol)
        energies[j] = exp(-1.0*fabs(currentE - bestE) / 2.5);
        totalPop += energies[j];
      }
      // now set the weights
      for (unsigned int j = 0; j < rotor->GetTorsionValues().size(); j++) {
        if (IsNear(worstE, bestE, 1.0e-3))
          weight = 1 / rotor->GetTorsionValues().size();
        else
          weight = energies[j]/totalPop;
        weightSet.push_back(weight);
      }
      rotorWeights.push_back(weightSet);
    }

    int best_conformer;
//    double penalty; // for poor performance
    double randFloat; // generated random number -- used to pick a rotor
    double total; // used to calculate the total probability 
    double *bestCoordPtr = new double [d->mol.NumAtoms() * 3]; // coordinates for best conformer
    
    // Start with the current coordinates
    bestE = worstE = Energy(false);
    // We're later going to add this back to the molecule as a new conformer
    memcpy((char*)bestCoordPtr,(char*)d->mol.GetCoordinates(),sizeof(double)*3*d->mol.NumAtoms());    
    
    // Now we actually test some weightings
    IF_OBFF_LOGLVL_LOW {
      snprintf(d->logbuf, BUFF_SIZE, "  GENERATED %d CONFORMERS\n\n", conformers);
      OBFFLog(d->logbuf);
      OBFFLog("CONFORMER     ENERGY\n");
      OBFFLog("--------------------\n");
    }

    double defaultRotor = 1.0/sqrt((double)rl.Size());
    for (unsigned int c = 0; c < conformers; ++c) {
      d->mol.SetCoordinates(initialCoord);

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
        for (unsigned int j = 0; j < rotor->GetTorsionValues().size(); j++) {
          if (randFloat > total && randFloat < (total+ rotorWeights[i][j])) {
            rotorKey[i] = j;
            break;      
          }
          else
            total += rotorWeights[i][j];
        }
      }
      rotamers.SetCurrentCoordinates(d->mol, rotorKey);

      m_loglvl = OBFF_LOGLVL_NONE;
      OBMinimize mini(this);
      mini.SteepestDescent(geomSteps); // energy minimization for conformer
      m_loglvl = origLogLevel;
      currentE = Energy(false);

      IF_OBFF_LOGLVL_LOW {
        snprintf(d->logbuf, BUFF_SIZE, "   %3d      %8.3f\n", c + 1, currentE);
        OBFFLog(d->logbuf);
      }
      
      if (!isfinite(currentE))
        continue;

      if (currentE < bestE) {
        bestE = currentE;
        memcpy((char*)bestCoordPtr,(char*)d->mol.GetCoordinates(),sizeof(double)*3*d->mol.NumAtoms());
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
      snprintf(d->logbuf, BUFF_SIZE, "\n  LOWEST ENERGY: %8.3f\n\n",
              bestE);
      OBFFLog(d->logbuf);
    }

    // debugging output to see final weightings of each rotor setting
    IF_OBFF_LOGLVL_HIGH {
      OBFFLog("Final Weights: \n");
      for (unsigned int i = 1; i < rotorWeights.size() - 1; ++i) {
        snprintf(d->logbuf, BUFF_SIZE, " Weight: %d", i);
        OBFFLog(d->logbuf);
        for (unsigned int j = 0; j < rotorWeights[i].size(); ++j) {
          snprintf(d->logbuf, BUFF_SIZE, " %8.3f", rotorWeights[i][j]);
          OBFFLog(d->logbuf);
        }
        OBFFLog("\n");
      }
    }

    d->mol.AddConformer(bestCoordPtr);
    d->current_conformer = d->mol.NumConformers() - 1;
    d->mol.SetConformer(d->current_conformer);
  }

  void OBForceField::DistanceGeometry() 
  {
    int N = d->mol.NumAtoms();
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
    FOR_ATOMS_OF_MOL (a, d->mol) {
      i = a->GetIdx() - 1;
      FOR_ATOMS_OF_MOL (b, d->mol) {
        j = b->GetIdx() - 1;

        if (&*a == &*b) {
          matrix[i][j] = 0.0; // diagonal
          continue;
        }
        // Find relationship
        is15 = true;
        FOR_NBORS_OF_ATOM (nbr1, d->mol.GetAtom(a->GetIdx())) { // 1-2
          if (&*nbr1 == &*b) {
            matrix[i][j] = 1.3;
            break;
          }
          FOR_NBORS_OF_ATOM (nbr2, d->mol.GetAtom(nbr1->GetIdx())) { // 1-3
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
          snprintf(d->logbuf, BUFF_SIZE, " %8.4f ", matrix[i][j]);
          OBFFLog(d->logbuf);
        }
        OBFFLog("]\n");
      }
      OBFFLog("\n");
    }

    // Triangle smoothing
    //FOR_ANGLES_OF_MOL(angle, d->mol) {
    int a, b, c;
    bool self_consistent = false;
    while (!self_consistent) {
      self_consistent = true;

      FOR_ATOMS_OF_MOL (_a, d->mol) {
        a = _a->GetIdx() - 1;
        FOR_ATOMS_OF_MOL (_b, d->mol) {
          if (&*_b == &*_a)
            continue;
          b = _b->GetIdx() - 1;
          FOR_ATOMS_OF_MOL (_c, d->mol) {
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
          snprintf(d->logbuf, BUFF_SIZE, " %8.4f ", matrix[i][j]);
          OBFFLog(d->logbuf);
        }
        OBFFLog("]\n");
      }
      OBFFLog("\n");
    }
    
    // Generate random distance matrix between lower and upper limits
    FOR_ATOMS_OF_MOL (a, d->mol) {
      i = a->GetIdx() - 1;
      FOR_ATOMS_OF_MOL (b, d->mol) {
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
      OBFFLog("RANDOM DISTANCE MATRIX BETWEEN LIMITS\n\n");
      for (i=0; i<N; i++) {
        OBFFLog("[");
        for (j=0; j<N; j++) {
          snprintf(d->logbuf, BUFF_SIZE, " %8.4f ", matrix[i][j]);
          OBFFLog(d->logbuf);
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
          snprintf(d->logbuf, BUFF_SIZE, " %8.4f ", G[i][j]);
          OBFFLog(d->logbuf);
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

    //matrix3x3::jacobi(N, (double *) &G, (double *) &eigenvalues, (double *) &eigenvectors);
    
    // output eigenvalues and eigenvectors
    IF_OBFF_LOGLVL_LOW {

      OBFFLog("EIGENVALUES OF METRIC MATRIX\n\n");
      for (i=0; i<N; i++) {
        snprintf(d->logbuf, BUFF_SIZE, "%8.4f ", eigenvalues[i]);
        OBFFLog(d->logbuf);
      }
      OBFFLog("\n");

      OBFFLog("EIGENVECTORS OF METRIC MATRIX\n\n");
      for (i=0; i<N; i++) {
        OBFFLog("[");
        for (j=0; j<N; j++) {
          snprintf(d->logbuf, BUFF_SIZE, " %8.4f ", eigenvectors[i][j]);
          OBFFLog(d->logbuf);
        }
        OBFFLog("]\n");
      }
      OBFFLog("\n");
    }

    // Assign coordinates
    double xa, ya, za;
    FOR_ATOMS_OF_MOL (a, d->mol) {
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
  // Cut-off  (XXXX05)
  //
  //////////////////////////////////////////////////////////////////////////////////
  
  void OBForceField::EnableCutOff(bool enable)
  {
      m_cutoff = enable;
  }
    
  void OBForceField::SetVDWCutOff(double r)
  {
    d->rvdw = r;
  }
    
  double OBForceField::GetVDWCutOff()
  {
    return d->rvdw;
  }
    
  void OBForceField::SetElectrostaticCutOff(double r)
  {
    d->rele = r;
  }
    
  double OBForceField::GetElectrostaticCutOff()
  {
    return d->rele;
  }
    
  void OBForceField::SetUpdateFrequency(int f)
  {
    d->pairfreq = f;
  }
    
  int OBForceField::GetUpdateFrequency()
  {
    return d->pairfreq;
  } 
 
  void OBForceField::UpdatePairsSimple()
  {
    unsigned int i = 0;

    d->vdwpairs.Clear();
    d->elepairs.Clear();

    FOR_PAIRS_OF_MOL(p, d->mol) {
      OBAtom *a = d->mol.GetAtom((*p)[0]);
      OBAtom *b = d->mol.GetAtom((*p)[1]);
 
      const Eigen::Vector3d ab = a->GetVector() - b->GetVector();
      double rab = ab.norm();
      // update vdw pairs
      if (rab < d->rvdw) {
        d->vdwpairs.SetBitOn(i);
      } else {
        d->vdwpairs.SetBitOff(i);
      }
      // update electrostatic pairs
      if (rab < d->rele) {
        d->elepairs.SetBitOn(i);
      } else {
        d->elepairs.SetBitOff(i);
      }
      
      i++;  
    }
    /* 
    IF_OBFF_LOGLVL_LOW {
      snprintf(d->logbuf, BUFF_SIZE, "UPDATE VDW PAIRS: %d --> %d (VDW), %d (ELE) \n", i+1, 
          d->vdwpairs.CountBits(), d->elepairs.CountBits());
      OBFFLog(d->logbuf);
    }
    */
  }
  
  unsigned int OBForceField::GetNumPairs()
  {
    unsigned int i = 1;
    FOR_PAIRS_OF_MOL(p, d->mol)
      i++;
    
    return i;
  }
  
  OBBitVec& OBForceField::GetVDWPairs()
  {
    return d->vdwpairs;
  }
    
  OBBitVec& OBForceField::GetElePairs()
  {
    return d->elepairs;
  }
 

  //////////////////////////////////////////////////////////////////////////////////
  //
  // Interaction groups  (XXXX06)
  //
  //////////////////////////////////////////////////////////////////////////////////
 

  void OBForceField::AddIntraGroup(OBBitVec &group)
  {
    d->intraGroup.push_back(group);
  }

  void OBForceField::AddInterGroup(OBBitVec &group)
  {
    d->interGroup.push_back(group);
  }

  void OBForceField::AddInterGroups(OBBitVec &group1, OBBitVec &group2)
  {
    pair<OBBitVec,OBBitVec> groups;
    groups.first = group1;
    groups.second = group2;
    d->interGroups.push_back(groups);
  }

  void OBForceField::ClearGroups()
  {
    d->intraGroup.clear();
    d->interGroup.clear();
    d->interGroups.clear();
  }
  
  bool OBForceField::HasGroups()
  {
    if (!d->intraGroup.empty())
      return true;

    if (!d->interGroup.empty())
      return true;

    if (!d->interGroups.empty())
      return true;

    return false;
  }
 
  vector<OBBitVec>& OBForceField::GetIntraGroup()
  {
    return d->intraGroup;
  }

  vector<OBBitVec>& OBForceField::GetInterGroup()
  {
    return d->interGroup;
  }

  vector<pair<OBBitVec, OBBitVec> >& OBForceField::GetInterGroups()
  {
    return d->interGroups;
  }

  bool OBForceField::DetectExplosion()
  {
    FOR_ATOMS_OF_MOL (atom, d->mol) {
      if (!isfinite(atom->GetX()))
        return true;
      if (!isfinite(atom->GetY()))
        return true;
      if (!isfinite(atom->GetZ()))
        return true;
    }
    
    FOR_BONDS_OF_MOL (bond, d->mol) {
      if (bond->GetLength() > 30.0)
        return true;
    }
        
    return false;
  }
  
  //////////////////////////////////////////////////////////////////////////////////
  //
  // Molecular Dynamics  (XXXX08)
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
    d->ncoords = d->mol.NumAtoms() * 3;
    int velocityIdx;
    double velocity, kB;
    kB = 0.00831451 / KCAL_TO_KJ; // kcal/(mol*K)
    
    d->velocityPtr = new double[d->ncoords];
    memset(d->velocityPtr, '\0', sizeof(double) * d->ncoords);
    
    FOR_ATOMS_OF_MOL (a, d->mol) {
      if (!d->constraints.IsFixed(a->GetIdx()) || (d->fixAtom == a->GetIdx()) || (d->ignoreAtom == a->GetIdx())) {
        velocityIdx = (a->GetIdx() - 1) * 3;
      
        // add twelve random numbers between 0.0 and 1.0,
        // subtract 6.0 from their sum, multiply with sqrt(kT/m)
        if (!d->constraints.IsXFixed(a->GetIdx())) {
          velocity = 0.0;
          for (int i=0; i < 12; ++i) 
            velocity += generator.NextFloat();
          velocity -= 6.0;
          velocity *= sqrt((kB * d->temp)/ (1000 * a->GetAtomicMass()));
          d->velocityPtr[velocityIdx] = velocity; // x10: gromacs uses nm instead of A
        }
        
        if (!d->constraints.IsYFixed(a->GetIdx())) {
          velocity = 0.0;
          for (int i=0; i < 12; ++i)
            velocity += generator.NextFloat();
          velocity -= 6.0;
          velocity *= sqrt((kB * d->temp)/ (1000 * a->GetAtomicMass()));
          d->velocityPtr[velocityIdx+1] = velocity; // idem
        }

        if (!d->constraints.IsZFixed(a->GetIdx())) {
          velocity = 0.0;
          for (int i=0; i < 12; ++i)
            velocity += generator.NextFloat();
          velocity -= 6.0;
          velocity *= sqrt((kB * d->temp)/ (1000 * a->GetAtomicMass()));
          d->velocityPtr[velocityIdx+2] = velocity; // idem
        }
      }
    }

    CorrectVelocities();
  }

  void OBForceField::CorrectVelocities()
  {
    d->ncoords = d->mol.NumAtoms() * 3;
    int velocityIdx;
    double velocity, kB, E_kin, E_kin2, factor;
    kB = 0.00831451 / KCAL_TO_KJ; // kcal/(mol*K)
 
    // E_kin = 0.5 * Ndf * kB * T
    E_kin = d->ncoords * kB * d->temp;
    //cout << "E_{kin} = Ndf * kB * T = " << E_kin << endl;
    
    // E_kin = 0.5 * sum( m_i * v_i^2 )
    E_kin2 = 0.0;
    FOR_ATOMS_OF_MOL (a, d->mol) {
      velocityIdx = (a->GetIdx() - 1) * 3;
      
      velocity = sqrt( d->velocityPtr[velocityIdx]   * d->velocityPtr[velocityIdx] + 
                       d->velocityPtr[velocityIdx+1] * d->velocityPtr[velocityIdx+1] + 
                       d->velocityPtr[velocityIdx+2] * d->velocityPtr[velocityIdx+2] );

      E_kin2 += 1000 * a->GetAtomicMass() * velocity * velocity;
    }
    //cout << "E_{kin} = sum( m_i * v_i^2 ) = " << E_kin2 << endl;

    // correct
    factor = sqrt(E_kin / E_kin2);
    FOR_ATOMS_OF_MOL (a, d->mol) {
      velocityIdx = (a->GetIdx() - 1) * 3;
      d->velocityPtr[velocityIdx]   *= factor; 
      d->velocityPtr[velocityIdx+1] *= factor; 
      d->velocityPtr[velocityIdx+2] *= factor; 
    }
    
    // E_kin = 0.5 * sum( m_i * v_i^2 )
    E_kin2 = 0.0;
    FOR_ATOMS_OF_MOL (a, d->mol) {
      velocityIdx = (a->GetIdx() - 1) * 3;
      
      velocity = sqrt( d->velocityPtr[velocityIdx]   * d->velocityPtr[velocityIdx] + 
                       d->velocityPtr[velocityIdx+1] * d->velocityPtr[velocityIdx+1] + 
                       d->velocityPtr[velocityIdx+2] * d->velocityPtr[velocityIdx+2] );

      E_kin2 += 1000 * a->GetAtomicMass() * velocity * velocity;
    }
    //cout << "E_{kin_corr} = sum( m_i * v_i^2 ) = " << E_kin2 << endl;
  }
  
  void OBForceField::MolecularDynamicsTakeNSteps(int n, double T, double timestep)
  {
    int coordIdx;
    double timestep2;
    Eigen::Vector3d force, pos, accel;
    d->timestep = timestep;
    d->temp = T;

    timestep2 = 0.5 * d->timestep * d->timestep;
    
    if (!d->velocityPtr)
      GenerateVelocities();
    Energy(true); // compute gradients
 
    for (int i = 1; i <= n; i++) {
      FOR_ATOMS_OF_MOL (a, d->mol) {
        if (!d->constraints.IsFixed(a->GetIdx()) || (d->fixAtom == a->GetIdx()) || (d->ignoreAtom == a->GetIdx())) {
          if (HasAnalyticalGradients())
            force = GetGradients()[a->GetIdx()-1];
          else
            force = NumericalDerivative(a->GetIdx()-1);

          pos = a->GetVector();
          coordIdx = (a->GetIdx() - 1) * 3;
	  
          // a(i) = F(i) / m
          accel = force / (1000 * a->GetAtomicMass());

          // x(i+1) = x(i) + v(i) /\t + 0.5 a(i) /\t^2
          pos = Eigen::Vector3d(
              pos.x() + d->velocityPtr[coordIdx]   * d->timestep + accel.x() * timestep2,
              pos.y() + d->velocityPtr[coordIdx+1] * d->timestep + accel.y() * timestep2,
              pos.z() + d->velocityPtr[coordIdx+2] * d->timestep + accel.z() * timestep2);
          a->SetVector(pos);
	  
          // v(i+.5) = v(i) + 0.5 a(i) /\t
          d->velocityPtr[coordIdx]   += 0.5 * accel.x() * d->timestep;
          d->velocityPtr[coordIdx+1] += 0.5 * accel.y() * d->timestep;
          d->velocityPtr[coordIdx+2] += 0.5 * accel.z() * d->timestep;

          Energy(true); // compute gradients
          
          if (HasAnalyticalGradients())
            force = GetGradients()[a->GetIdx()-1];
          else
            force = NumericalDerivative(a->GetIdx()-1);
          
          // a(i+1) = F(i+1) / m
          accel = force / (1000 * a->GetAtomicMass());

          d->velocityPtr[coordIdx]   += 0.5 * accel.x() * d->timestep;
          d->velocityPtr[coordIdx+1] += 0.5 * accel.y() * d->timestep;
          d->velocityPtr[coordIdx+2] += 0.5 * accel.z() * d->timestep;

        }
      }
      if (i % 10 == 0) 
        CorrectVelocities();
    }
  }
 
  //////////////////////////////////////////////////////////////////////////////////
  //
  // Misc.  (XXXX09)
  //
  //////////////////////////////////////////////////////////////////////////////////

  bool OBForceField::IsInSameRing(OBAtom* a, OBAtom* b)
  {
    bool a_in, b_in;
    vector<OBRing*> vr;
    vr = d->mol.GetSSSR();
    
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

  OBGridData* OBForceField::GetGrid(double step, double padding, const char* type, double pchg)
  {
    cout << "OBForceFieldMMFF94::GetGrid(" << step << ", " << type << ")" << endl;
    OBFloatGrid fgrid;
    fgrid.Init(d->mol, step, padding);
    Eigen::Vector3d min;
    unsigned int xDim, yDim, zDim, xyzDim;

    min = fgrid.GetMin();

    xDim = fgrid.GetXdim();
    yDim = fgrid.GetYdim();
    zDim = fgrid.GetZdim();
    xyzDim = xDim * yDim * zDim;

    cout << "xDim = " << xDim << ", yDim = " << yDim << ", zDim = " << zDim << endl;

    // Add the probe atom
    d->mol.BeginModify();
    OBAtom *atom = d->mol.NewAtom();
    int index = atom->GetIdx();
    d->mol.EndModify();
    SetTypes();
    atom->SetType(type);
    atom->SetPartialCharge(pchg); 
    
    SetupCalculations();

    atom = d->mol.GetAtom(index);
    double *pos = atom->GetCoordinate();

    Eigen::Vector3d coord;
    double evdw, eele;
    double distance, minDistance;

    OBGridData *grid = new OBGridData;
    Eigen::Vector3d xAxis, yAxis, zAxis;
    xAxis = Eigen::Vector3d(step, 0.0, 0.0);
    yAxis = Eigen::Vector3d(0.0, step, 0.0);
    zAxis = Eigen::Vector3d(0.0, 0.0, step);

    grid->SetNumberOfPoints(xDim, yDim, zDim);
    grid->SetLimits(min, xAxis, yAxis, zAxis);
    
    // VDW surface
    for (unsigned int i = 0; i < xDim; ++i) {
      coord[0] = min[0] + i * step;
      for (unsigned int j = 0; j < yDim; ++j) {
        coord[1] = min[1] + j * step;
        for (unsigned int k = 0; k < zDim; ++k)
        {
          coord[2] = min[2] + k * step;
          minDistance = 1.0E+10;
	  FOR_ATOMS_OF_MOL (a, d->mol) {
            if (a->GetIdx() == atom->GetIdx())
	      continue;
	    if (a->IsHydrogen())
	      continue;

            Eigen::Vector3d v = coord - a->GetVector();
	    distance = v.norm();

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
      coord[0] = min[0] + i * step;
      for (unsigned int j = 0; j < yDim; ++j) {
        coord[1] = min[1] + j * step;
        for (unsigned int k = 0; k < zDim; ++k)
        {
	  coord[2] = min[2] + k * step;
          
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
    
    d->mol.BeginModify();
    d->mol.DeleteAtom(atom);
    d->mol.EndModify();

    return grid;
  }

  //////////////////////////////////////////////////////////////////////////////////
  //
  // Plugin stuff 
  //
  //////////////////////////////////////////////////////////////////////////////////
  
  OBForceField::OBForceField(const char* ID, bool IsDefault) : d(new OBForceFieldPrivate)
  {
    _id=ID;
    
    if (ID && *ID) {
      if (IsDefault || Map().empty()) 
        Default() = this;
        
      Map()[ID]=this;
      PluginMap()[TypeID()] =this;
    }
    
    m_cutoff = false;
    SetAllTermsEnabled(true);
        
    d->init = false;
    d->rvdw = 7.0;
    d->rele = 15.0;
    d->pairfreq = 10;
    d->linesearch = LineSearchType::Simple;
  }
  
  OBForceField* OBForceField::FindType(const char* ID)
  {
    if(!ID || *ID==0) 
      return Default();
      
    return static_cast<OBForceField*>(BaseFindType(Map(),ID));
  }
 
  OBForceField::~OBForceField()
  {
  }


} // end namespace OpenBabel


//! \file forcefield.cpp
//! \brief Handle OBForceField class

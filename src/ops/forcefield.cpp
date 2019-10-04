/**********************************************************************
forcefield.cpp - A OBOp to calculate and minimize the energy using a
                 forcefield (re-wrap of obminimize and obenergy)

Copyright (C) 2006-2007 by Tim Vandermeersch
          (C) 2009 by Chris Morley

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

/******************************************************************************
**** CURRENTLY ONLY SUITABLE FOR USE WITH THE OBABEL COMMANDLINE INTERFACE ****
This allows options to have parameters, e.g. --ff Ghemical
Compile with tools/obabel.cpp rather than tools/babel.cpp

*******************************************************************************/

#include <openbabel/babelconfig.h>
#include <iostream>
#include <cstdlib>
#include<openbabel/op.h>
#include<openbabel/mol.h>
#include<openbabel/forcefield.h>
#include<openbabel/generic.h>

namespace OpenBabel
{
  using namespace std;

  //////////////////////////////////////////////////////////
  //
  //  OpEnergy
  //
  //////////////////////////////////////////////////////////

  class OpEnergy : public OBOp
  {
    public:
      OpEnergy(const char *ID) : OBOp(ID, false) {}

      const char* Description()
      {
        return "ForceField Energy Evaluation (not displayed in GUI)\n"
          "Typical usage: obabel infile.xxx -O outfile.yy --energy --log\n"
          " options:         description\n"
          " --log        output a log of the energies (default = no log)\n"
          " --epsilon #  set the dielectric constant for electrostatics\n"
          " --ff #       select a forcefield (default = MMFF94)\n"
          " The hydrogens are always made explicit before energy evaluation.\n"
          " The energy is put in an OBPairData object \"Energy\" which is\n"
          "   accessible via an SDF or CML property or --append (to title).\n"
          ;
      }

      virtual bool WorksWith(OBBase* pOb) const
      {
        return dynamic_cast<OBMol*>(pOb) != NULL;
      }
      virtual bool Do(OBBase* pOb, const char* OptionText, OpMap* pmap, OBConversion*);
  };

  //////////////////////////////////////////////////////////
  OpEnergy theOpEnergy("energy"); //Global instance

  //////////////////////////////////////////////////////////
  bool OpEnergy::Do(OBBase* pOb, const char* OptionText, OpMap* pmap, OBConversion*)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(!pmol)
      return false;
    pmol->AddHydrogens(false, false);

    bool log = false;

    string ff = "MMFF94";
    double epsilon = 1.0;
    OpMap::const_iterator iter = pmap->find("ff");
    if(iter!=pmap->end())
      ff = iter->second;
    OBForceField* pFF = OBForceField::FindForceField(ff);
    iter = pmap->find("epsilon");
    if (iter!=pmap->end())
      epsilon = atof(iter->second.c_str());

    iter = pmap->find("log");
    if(iter!=pmap->end())
      log=true;

    // set some force field variables
    pFF->SetLogFile(&clog);
    pFF->SetLogLevel(log ? OBFF_LOGLVL_MEDIUM : OBFF_LOGLVL_NONE);

    pFF->SetDielectricConstant(epsilon);
    if (!pFF->Setup(*pmol)) {
      cerr  << "Could not setup force field." << endl;
      return false;
    }

    //Put the energy in a OBPairData object
    OBPairData *dp = new OBPairData;
    dp->SetAttribute("Energy");
    stringstream ss;
    ss << pFF->Energy(false);
    dp->SetValue(ss.str());
    dp->SetOrigin(fileformatInput);
    pmol->SetData(dp);

    return true;
  }

  //////////////////////////////////////////////////////////
  //
  //  OpMinimize
  //
  //////////////////////////////////////////////////////////

  class OpMinimize : public OBOp
  {
    public:
      OpMinimize(const char* ID) : OBOp(ID, false) {}

      const char* Description()
      {
        return "ForceField Energy Minimization (not displayed in GUI)\n"
          "Typical usage: obabel infile.xxx -O outfile.yyy --minimize --steps 1500 --sd\n"
          " options:         description:\n"
          " --log        output a log of the minimization process(default= no log)\n"
          " --crit #     set convergence criteria (default=1e-6)\n"
          " --sd         use steepest descent algorithm (default = conjugate gradient)\n"
          " --newton     use Newton2Num linesearch (default = Simple)\n"
          " --ff #       select a forcefield (default = Ghemical)\n"
          " --steps #    specify the maximum number of steps (default = 2500)\n"
          " --cut        use cut-off (default = don't use cut-off)\n"
          " --epsilon #  relative dielectric constant (default = 1.0)\n"
          " --rvdw #     specify the VDW cut-off distance (default = 6.0)\n"
          " --rele #     specify the Electrostatic cut-off distance (default = 10.0)\n"
          " --freq #     specify the frequency to update the non-bonded pairs (default = 10)\n"
          " The hydrogens are always made explicit before minimization.\n"
          " The energy is put in an OBPairData object \"Energy\" which is\n"
          "   accessible via an SDF or CML property or --append (to title).\n"
          ;
      }

      virtual bool WorksWith(OBBase* pOb) const
      {
        return dynamic_cast<OBMol*>(pOb) != NULL;
      }
      virtual bool Do(OBBase* pOb, const char* OptionText, OpMap* pmap, OBConversion*);
  };

  //////////////////////////////////////////////////////////
  OpMinimize theOpMinimize("minimize"); //Global instance

  //////////////////////////////////////////////////////////
  bool OpMinimize::Do(OBBase* pOb, const char* OptionText, OpMap* pmap, OBConversion*)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(!pmol)
      return false;
    pmol->AddHydrogens(false, false);

    int steps = 2500;
    double crit = 1e-6;
    bool sd = false;
    bool cut = false;
    bool newton = true;
    double epsilon = 1.0;
    double rvdw = 6.0;
    double rele = 10.0;
    int freq = 10;
    bool log = false;

    string ff = "MMFF94";
    OpMap::const_iterator iter = pmap->find("ff");
    if(iter!=pmap->end())
      ff = iter->second;
    OBForceField* pFF = OBForceField::FindForceField(ff);

    iter = pmap->find("sd");
    if(iter!=pmap->end())
      sd=true;

    iter = pmap->find("newton");
    if(iter!=pmap->end())
      newton=true;

    iter = pmap->find("cut");
    if(iter!=pmap->end())
      cut=true;

    iter = pmap->find("crit");
    if(iter!=pmap->end())
      crit = atof(iter->second.c_str());

    iter = pmap->find("steps");
    if(iter!=pmap->end())
      steps = atoi(iter->second.c_str());

    iter = pmap->find("epsilon");
    if(iter!=pmap->end())
      epsilon = atof(iter->second.c_str());

    iter = pmap->find("rvdw");
    if(iter!=pmap->end())
      rvdw = atof(iter->second.c_str());

    iter = pmap->find("rele");
    if(iter!=pmap->end())
      rele = atof(iter->second.c_str());

    iter = pmap->find("pf");
    if(iter!=pmap->end()) {
      freq = atoi(iter->second.c_str());
      if (freq < 1)
        freq = 10; // don't divide by zero
    }

    iter = pmap->find("log");
    if(iter!=pmap->end())
      log=true;

    if (newton)
      pFF->SetLineSearchType(LineSearchType::Newton2Num);

    // set some force field variables
    pFF->SetLogFile(&clog);
    pFF->SetLogLevel(log ? OBFF_LOGLVL_LOW : OBFF_LOGLVL_NONE);
    pFF->SetVDWCutOff(rvdw);
    pFF->SetElectrostaticCutOff(rele);
    pFF->SetUpdateFrequency(freq);
    pFF->SetDielectricConstant(epsilon);
    pFF->EnableCutOff(cut);

    if (!pFF->Setup(*pmol)) {
      cerr  << "Could not setup force field." << endl;
      return false;
    }

    bool done = true;
    if (sd)
      pFF->SteepestDescent(steps, crit);
    else
      pFF->ConjugateGradients(steps, crit);

    pFF->GetCoordinates(*pmol);

    //Put the energy in a OBPairData object
    OBPairData *dp = new OBPairData;
    dp->SetAttribute("Energy");
    stringstream ss;
    ss << pFF->Energy(false);
    dp->SetValue(ss.str());
    dp->SetOrigin(fileformatInput);
    pmol->SetData(dp);

    return true;
  }

}//namespace

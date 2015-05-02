/**********************************************************************
conformer.cpp - A OBOp to calculate and minimize the energy using a
                 forcefield (re-wrap of obminimize and obenergy)

Copyright (C) 2010 by Tim Vandermeersch

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
#include <openbabel/op.h>
#include <openbabel/mol.h>
#include <openbabel/forcefield.h>
#include <openbabel/conformersearch.h>
#include <openbabel/generic.h>

namespace OpenBabel
{
  using namespace std;

  //////////////////////////////////////////////////////////
  //
  //  OpEnergy
  //
  //////////////////////////////////////////////////////////

  class OpConformer : public OBOp
  {
    public:
      OpConformer(const char *ID) : OBOp(ID, false) {}

      const char* Description()
      {
        return "Conformer Searching (not displayed in GUI)\n"
          "Typical usage: obabel infile.xxx -O outfile.yy --conformer --nconf\n"
          " options:             description\n"
          " --log            output a log of the energies (default = no log)\n"
          " --nconf #        number of conformers to generate\n"
          " forcefield based methods for finding stable conformers:\n"
          " --systematic     systematically generate all conformers\n"
          " --fast           fast systematic search from central bonds\n"
          " --random         randomly generate conformers\n"
          " --weighted       weighted rotor search for lowest energy conformer\n"
          " --ff #           select a forcefield (default = MMFF94)\n"
          " genetic algorithm based methods (default):\n"
          " --children #     number of children to generate for each parent (default = 5)\n"
          " --mutability #   mutation frequency (default = 5)\n"
          " --converge #     number of identical generations before convergence is reached\n"
          " --score #        scoring function [rmsd|energy|minrmsd|minenergy] (default = rmsd)\n"
          ;
      }

      virtual bool WorksWith(OBBase* pOb) const
      {
        return dynamic_cast<OBMol*>(pOb) != NULL;
      }
      virtual bool Do(OBBase* pOb, const char* OptionText, OpMap* pmap, OBConversion*);
  };

  //////////////////////////////////////////////////////////
  OpConformer theOpConformer("conformer"); //Global instance

  void getInteger(const std::string &str, int &value)
  {
    std::istringstream iss(str);
    iss >> value;
  }

  //////////////////////////////////////////////////////////
  bool OpConformer::Do(OBBase* pOb, const char* OptionText, OpMap* pmap, OBConversion*)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(!pmol)
      return false;
    pmol->AddHydrogens(false, false);

    OpMap::const_iterator iter;
    bool log = false;
    bool systematic = false;
    bool random = false;
    bool weighted = false;
    bool fast = false;
    int numConformers = 30;

    iter = pmap->find("log");
    if(iter!=pmap->end())
      log=true;

    iter = pmap->find("nconf");
    if(iter!=pmap->end())
      getInteger(iter->second, numConformers);

    iter = pmap->find("systematic");
    if(iter!=pmap->end())
      systematic = true;

    iter = pmap->find("random");
    if(iter!=pmap->end())
      random = true;

    iter = pmap->find("fast");
    if(iter!=pmap->end())
      fast = true;

    iter = pmap->find("weighted");
    if(iter!=pmap->end())
      weighted = true;

    if (systematic || random || fast || weighted) {
      std::string ff = "MMFF94";
      iter = pmap->find("ff");
      if(iter!=pmap->end())
        ff = iter->second;
      OBForceField* pFF = OBForceField::FindForceField(ff);

      // set some force field variables
      pFF->SetLogFile(&clog);
      pFF->SetLogLevel(log ? OBFF_LOGLVL_MEDIUM : OBFF_LOGLVL_NONE);

      // Add cut-offs for faster conformer searching
      // Generally people will perform further optimization on a final conformer
      pFF->EnableCutOff(true);
      pFF->SetVDWCutOff(10.0);
      pFF->SetElectrostaticCutOff(20.0);
      pFF->SetUpdateFrequency(10); // delay updates of non-bonded distances

      if (!pFF->Setup(*pmol)) {
        cerr  << "Could not setup force field." << endl;
        return false;
      }

      // Perform search
      if (systematic) {
        pFF->SystematicRotorSearch(10); // 10 steepest-descent forcfield steps per conformer
      } else if (fast) {
        pFF->FastRotorSearch(true); // permute rotors
      } else if (random) {
        pFF->RandomRotorSearch(numConformers, 10);
      } else if (weighted) {
        pFF->WeightedRotorSearch(numConformers, 10);
      }
      pFF->GetConformers(*pmol);

    } else { // GA-based searching
      int numChildren = 5;
      int mutability = 5;
      int convergence = 5;
      std::string score = "rmsd";

      iter = pmap->find("children");
      if(iter!=pmap->end())
        getInteger(iter->second, numChildren);

      iter = pmap->find("mutability");
      if(iter!=pmap->end())
        getInteger(iter->second, mutability);

      iter = pmap->find("convergence");
      if(iter!=pmap->end())
        getInteger(iter->second, convergence);

      iter = pmap->find("score");
      if(iter!=pmap->end())
        score = iter->second;

      OBConformerSearch cs;
      if (score == "energy")
        cs.SetScore(new OBEnergyConformerScore);
      else if (score == "mine" || score == "minenergy")
        cs.SetScore(new OBMinimizingEnergyConformerScore);
      else if (score == "minr" || score == "minrmsd")
        cs.SetScore(new OBMinimizingRMSDConformerScore);

      if (cs.Setup(*pmol, numConformers, numChildren, mutability, convergence)) {
        cs.Search();
        cs.GetConformers(*pmol);
      }
    }

    return true;
  }


}//namespace

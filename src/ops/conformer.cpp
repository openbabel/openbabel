/**********************************************************************
conformer.cpp - A OBOp to calculate and minimize the energy using a
                 forcefield (re-wrap of obminimize and obenergy)

Copyright (C) 2010 by Tim Vandermeersch
Some portions Copyright (C) 2016 Torsten Sachse

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
#include <vector>
#include <stdio.h>
#include <openbabel/op.h>
#include <openbabel/mol.h>
#include <openbabel/bond.h>
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
          " --rings          sample ring torsions\n"
          " genetic algorithm (GA) based methods (default):\n"
          " --children #     number of children to generate for each parent (default = 5)\n"
          " --mutability #   mutation frequency (default = 5)\n"
          " --convergence #  number of identical generations before convergence is reached\n"
          " --score #        scoring function [rmsd|energy|minrmsd|minenergy] (default = rmsd)\n"
          " customize the filter used to sort out wrong conformers\n"
          " --csfilter #     the filtering algorithm [steric] (default=steric)\n"
          " --cutoff #       absolute distance in Anstroms below which atoms are considered to clash\n"
          " --vdw-factor #   apply this factor to all van-der-Waals radii before detecting clashes\n"
          " --ignore-H       do not detect clashes with hydrogen atoms\n"
          " customize rotors for GA based method (a rotor is given as a-b where a and b are indices)\n"
          " --printrot       print all possible rotors and do not perform a conformer search\n"
          " --onlyrot        declare a comma-separated list of rotors, all but these will be disregarded\n"
          " --norot          declare a comma-separated list of rotors that shall be disregarded\n"
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

  template<typename T>
  bool getValue(const std::string &str, T &value)
  {
    std::istringstream iss(str);
    bool ret = static_cast<bool>(iss >> value);
    return ret;
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
    bool rings = false;
    int numConformers = 30;

    iter = pmap->find("log");
    if(iter!=pmap->end())
      log=true;

    iter = pmap->find("nconf");
    if(iter!=pmap->end())
      getValue<int>(iter->second, numConformers);

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

    iter = pmap->find("ring");
    if(iter!=pmap->end())
        rings = true;

    iter = pmap->find("rings");
    if(iter!=pmap->end())
        rings = true;

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
        pFF->SystematicRotorSearch(10, rings); // 10 steepest-descent forcfield steps per conformer
      } else if (fast) {
        pFF->FastRotorSearch(true); // permute rotors
      } else if (random) {
        pFF->RandomRotorSearch(numConformers, 10, rings);
      } else if (weighted) {
        pFF->WeightedRotorSearch(numConformers, 10, rings);
      }
      pFF->GetConformers(*pmol);

    } else { // GA-based searching
      //default values for score
      int numChildren = 5;
      int mutability = 5;
      int convergence = 5;
      std::string score = "rmsd";
      //default values for filter
      double cutoff = 0.8 ;
      double vdw_factor = 0.5 ;
      bool check_hydrogens = true ;
      std::string filter = "steric";
      bool rotchange = false;
      bool rotorsoff = false;
      OBBitVec extrotors;

      iter = pmap->find("children");
      if(iter!=pmap->end())
        getValue<int>(iter->second, numChildren);

      iter = pmap->find("mutability");
      if(iter!=pmap->end())
        getValue<int>(iter->second, mutability);

      iter = pmap->find("convergence");
      if(iter!=pmap->end())
        getValue<int>(iter->second, convergence);

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

      iter = pmap->find("csfilter");
      if(iter!=pmap->end())
        filter = iter->second;
 
      iter = pmap->find("vdw-factor");
      if(iter!=pmap->end())
        getValue<double>(iter->second, vdw_factor);
 
      iter = pmap->find("cutoff");
      if(iter!=pmap->end())
        getValue<double>(iter->second, cutoff);
      
      iter = pmap->find("ignore-H");
      if(iter!=pmap->end())
        check_hydrogens = false;
 
      if (filter == "steric")
        cs.SetFilter(new OBStericConformerFilter(cutoff, vdw_factor, check_hydrogens));

      iter = pmap->find("printrot");
      if(iter!=pmap->end())
        cs.PrintRotors(true);

      iter = pmap->find("onlyrot");
      if(iter!=pmap->end()){
        rotchange = true;
        rotorsoff = false;
      }
      iter = pmap->find("norot");
      if(iter!=pmap->end()){
        rotchange = true;
        rotorsoff = true;
      }
      if (rotchange){
        if (rotorsoff){
          iter = pmap->find("norot");
        }else{
          iter = pmap->find("onlyrot");
        }
        OBBitVec fixbonds;
        if (rotorsoff){
          fixbonds.SetRangeOff(0,pmol->NumBonds()-1);
        }else{
          fixbonds.SetRangeOn( 0,pmol->NumBonds()-1);
          int end_it = fixbonds.EndBit();
          int it = fixbonds.FirstBit();
          while(it!=end_it){
            OBBond* b = pmol->GetBond(it);
            if (!(b->IsRotor())){
              fixbonds.SetBitOff(it);
            }
            it = fixbonds.NextBit(it);
          }
        }
        std::stringstream input(iter->second);
        std::string segment;
        while (std::getline(input, segment, ',')){
          int p1, p2;
          char additional;
          int converted = sscanf(segment.c_str(), "%d-%d%c", &p1, &p2, &additional);
          if (converted==2 && p1>0 && p2>0){
            if (p1>pmol->NumAtoms() || p2>pmol->NumAtoms()){
              std::cerr << "ERROR at least one of the atom indices " << p1 << " or " << p2 
                        << " is greater than the number of atoms " << pmol->NumAtoms() << "." << std::endl << std::flush;
              return false;
            }
            OBBond* b = pmol->GetBond(p1,p2);
            if (!b){
              std::cerr << "ERROR atoms " << p1 << " and " << p2 << " form no bond." << std::endl << std::flush;
              return false;
            }
            if (!(b->IsRotor())){
                std::cerr << "ERROR bond formed by atoms " << p1 << " and " << p2 << " is not rotable." << std::endl << std::flush;
                return false;
            }
            int idx = b->GetIdx();
            if (rotorsoff){
              fixbonds.SetBitOn(idx);
            }else{
              fixbonds.SetBitOff(idx);
            }
          }else{
              std::cerr << "ERROR parsing '" << segment << "' as rotor (must be i-j where i and j are indices >0)" << std::endl << std::flush;
              return false;
          }
        }
        cs.SetFixedBonds(fixbonds);
      }

      if (cs.Setup(*pmol, numConformers, numChildren, mutability, convergence)) {
        cs.Search();
        cs.GetConformers(*pmol);
      }
    }

    return true;
  }


}//namespace

/**********************************************************************
kekulize.cpp - Alternate algorithm to kekulize a molecule.

Copyright (C) 2004-2006 by Fabien Fontaine
Some portions Copyright (C) 2005-2006 by Geoffrey R. Hutchison
 
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

#include <openbabel/mol.h>

#include <sstream>

#define SINGLE 1
#define DOUBLE 2

using namespace std;

namespace OpenBabel 
{

  ///////////////////////////////////////////////////////////////////////////////
  //! \brief Kekulize aromatic rings without using implicit valence

  //! This new perceive kekule bonds function has been designed to 
  //! handle molecule files without explicit hydrogens such as pdb or xyz.
  //! (It can, of course, easily handle explicit hydrogens too.)
  //! The function does not rely on GetImplicitValence function
  //! The function looks for groups of aromatic cycle 
  //! For each group it tries to guess the number of electrons given by each atom
  //! in order to satisfy the huckel (4n+2) rule
  //! If the huckel rule cannot be satisfied the algorithm try with its best alternative guess
  //! Then it recursively walk on the atoms of the cycle and assign single and double bonds
  void OBMol::NewPerceiveKekuleBonds() 
  {

    if (HasKekulePerceived())  return;
    SetKekulePerceived();

    OBAtom *atom;
    int n, de, minde;
    std::vector<OBAtom*> cycle;
    OBBitVec avisit,cvisit;
    avisit.Resize(NumAtoms()+1);
    cvisit.Resize(NumAtoms()+1);
    OBBond *bond;
    std::vector<OBBond*>::iterator bi;
    std::vector<int> electron;
    int sume, orden;
    int bestatom = 1;
    int bestorden = 99;
    // Init the kekulized bonds
    unsigned i;
    FOR_BONDS_OF_MOL(bond, *this)
      {
        switch (bond->GetBO())
          {
          case 2: bond->SetKDouble(); break;
          case 3: bond->SetKTriple(); break;
          case 1: bond->SetKSingle(); break;
          }
      }

    // Find all the groups of aromatic cycle
    for(i=1; i<= NumAtoms(); i++ ) {
      atom = GetAtom(i);
      if (atom->HasAromaticBond() && !cvisit[i]) { // is new aromatic atom of an aromatic cycle ?

        avisit.Clear();
        electron.clear();
        cycle.clear();
      
        avisit.SetBitOn(i);
        expandcycle(atom, avisit);
        //store the atoms of the cycle(s)
        unsigned int j;
        for(j=1; j<= NumAtoms(); ++j) {
          if ( avisit[j] ) {
            atom = GetAtom(j);
            cycle.push_back(atom);
          }
        }

        // This isn't a real aromatic cycle -- give up
        // Fixes PR#1965566
        // Fixes PR#1784204
        if (cycle.size() < 3)
          continue;

        // At the beginning each atom give one electron to the cycle
        for(j=0; j< cycle.size(); ++j) {
          electron.push_back(1);
        }
      
        // remove one electron if the atom make a double bond out of the cycle
        sume =0;
        for(j=0; j< cycle.size(); ++j) {
          atom = cycle[j];
          for(bond = atom->BeginBond(bi); bond; bond = atom->NextBond(bi)) {
            if ( bond->IsDouble() ) {
              OBAtom *atom2 = bond->GetNbrAtom(atom);
              int fcharge = atom->GetFormalCharge();
              int fcharge2 = atom2->GetFormalCharge();
              if(atom->IsNitrogen() && atom2->IsOxygen() 
                 && fcharge == 0 && fcharge2 == 0) { //n=O to [n+][O-]
                atom->SetFormalCharge(1);
                atom2->SetFormalCharge(-1);
                bond->SetKSingle();
                bond->SetBO(1);
              }
              else {
                electron[j] = 0;
              }
            }
          }
          // count the number of electrons
          sume += electron[j];
        }

        // Save the electron state in case huckel rule is not satisfied
        vector<int> previousElectron = electron;

        // find the ideal number of electrons according to the huckel 4n+2 rule
        minde=99;
        for (i=1; 1; ++i) {
          n = 4 *i +2;
          de = n - sume;
          if (  de < minde )
            minde=de;
          else if ( minde < 0 )
            minde=de;
          else
            break;
        }
      
        stringstream errorMsg;

        //cout << "minde before:" << minde << endl;
        // if huckel rule not satisfied some atoms must give more electrons
        //cout << "minde " << minde << endl;
        while ( minde != 0 ) {
          bestorden=99;
          for(j=0; j< cycle.size(); ++j) {
            if (electron[j] == 1) {
              orden = getorden(cycle[j]);
              if (orden < bestorden) {
                bestorden = orden;
                bestatom = j;
              }
            }
          }
          if (bestorden==99) {  // no electron giving atom found
            errorMsg << "Kekulize: Huckel rule not satisfied for molecule " << GetTitle() << endl;
            obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
            break;             // Huckel rule cannot be satisfied
          }                    // try to kekulize anyway
          else {
            electron[bestatom] += 1;
            minde--;
          }
        }

        if (bestorden == 99) { // Huckel rule not satisfied, just try to get an even number of electron before kekulizing
	
          electron = previousElectron; // restore electon's state
	
          int odd = sume % 2; 
          //cout << "odd:" << odd << endl;
          if(odd) { // odd number of electrons try to add an electron to the best possible atom
            for(j=0; j< cycle.size(); ++j) {
              if (electron[j] == 1) {
                orden = getorden(cycle[j]);
                if (orden < bestorden) {
                  bestorden = orden;
                  bestatom = j;
                }
              }
            }
            if (bestorden==99) {  // no electron giving atom found
              errorMsg << "Kekulize: Cannot get an even number of electron for molecule " << GetTitle() << "\n";
              obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
              break;             // impossible to choose an atom to obtain an even number of electron
            }                    // try to kekulize anyway
            else {
              electron[bestatom] += 1;
            }
          }
        }
      
        //cout << "minde after:" << minde <<endl;
        //for(j=0; j < cycle.size(); ++j) {
        //OBAtom *cycleAtom = cycle[j];
        //cout << "\t" << cycleAtom->GetIdx();
        //}
        //cout << endl;

        //	for(j=0; j < electron.size(); ++j) {
        //	cout << "\t" << electron[j];
        //	}
        //	cout << endl;

        // kekulize the cycle(s)
        start_kekulize(cycle,electron);

        // Set the kekulized cycle(s) as visited
        for(j=1; j<= NumAtoms(); ++j) {
          if (avisit[j])
            cvisit.SetBitOn(j);
        }

      }
    }
    // Double bond have been assigned, set the remaining aromatic bonds to single
    //std::cout << "Set not assigned single bonds\n"; 
    FOR_BONDS_OF_MOL(bond, *this)
      {
        //std::cout << "bond " << bond->GetBeginAtomIdx() << " " << bond->GetEndAtomIdx() << " ";   
        if (bond->GetBO()==5 ) {
          bond->SetKSingle();
          bond->SetBO(1);
          //std::cout << "single\n";
        }
        //else
        //  std::cout << "double\n";
      }

    return;
  }

  ///////////////////////////////////////////////////////////////////////////////////////
  //! \brief Start kekulizing one or a fused set of aromatic ring(s)

  //! The initial electronic state indicates if an atoms must make a double bond or not
  //! Kekulizing is attempted recursively for all the atoms bonded to the first atom
  //! of the cycle. 
  void OBMol::start_kekulize( std::vector <OBAtom*> &cycle, std::vector<int> &electron) {
  
    std::vector<int> initState;
    std::vector<int> currentState;
    std::vector<int> binitState;
    std::vector<int> bcurrentState;
    std::vector<bool> mark;
    unsigned int Idx;
    OBAtom *atom, *atom2;
    OBBond *bond;

    //init the atom arrays
    unsigned i;
    for(i=0;i <NumAtoms()+1; ++i) {
      initState.push_back(-1);
      currentState.push_back(-1);
      mark.push_back(false);
    }
  
    //init the bond arrays with single bonds
    for(i=0;i <NumBonds(); ++i) {
      binitState.push_back(SINGLE);
      bcurrentState.push_back(SINGLE);
    }
  
    //set the electron number
    for(i=0; i< cycle.size(); ++i) {
      atom = cycle[i];
      Idx =  atom->GetIdx();
      if ( electron[i] == 1)
        initState[Idx] = 1; // make 1 double bond
      else
        initState[Idx] = 2; // make 2 single bonds

      currentState[Idx] = initState[Idx];
    }

    std::vector<OBBond*>::iterator b;
    OBAtom *nbr;
  
    // bool second_pass=false;
    // for( i=1; i<= NumAtoms(); ++i) {
    //     if(currentState[i] == 1) { // the atom can make a double bond
    //       atom = GetAtom(i);
    //       //find a neighbour that can make a double bond
    //       // and start kekulize
    //       for (nbr = atom->BeginNbrAtom(b);nbr;nbr = atom->NextNbrAtom(b)) {
    // 	if(currentState[nbr->GetIdx()]==1){
    // 	  if(!expand_kekulize(atom,nbr,currentState,initState, bcurrentState,binitState, mark)) {
    // 	    second_pass=true;
    // 	  }
	    
    // 	}
    //       }
    //     }
    //   }
    bool expand_successful=true;//was previously not initialized. true works but may not be correct
    atom = cycle[0];
    for (nbr = atom->BeginNbrAtom(b);nbr;nbr = atom->NextNbrAtom(b)) {
      if(initState[nbr->GetIdx()] == -1) //neighbor atom not in the cycle, try next one
        continue; 
      //std::cout << "Expand kekulize\n";
      expand_kekulize(atom,nbr,currentState,initState, bcurrentState,binitState, mark) ; 
      //Control that all the electron have been given to the cycle(s)
      expand_successful = true;
      for(unsigned i=0; i< cycle.size(); ++i) {
        atom2 = cycle[i];
        Idx =  atom2->GetIdx();
        //cout << "\t" << currentState[Idx];
        if (currentState[Idx] == 1)
          expand_successful=false;
      }
      //cout << endl;
      if (expand_successful)
        break;
      else {
        unsigned i;
        for(i=0;i <NumAtoms()+1; ++i) {
          currentState[i]=initState[i];
          mark[i]=false;
        }
        for(i=0;i <NumBonds(); ++i) {
          bcurrentState[i]=binitState[i];
        }   
      }
    }
    if (!expand_successful)
      {
        stringstream errorMsg;
        errorMsg << "Kekulize Error for molecule " << GetTitle() << endl;
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
      }

    // Set the double bonds
    // std::cout << "Set double bonds\n";
    for(i=0;i <NumBonds(); ++i) {
      bond = GetBond(i);    
      // std::cout << "bond " << bond->GetBeginAtomIdx() << " " << bond->GetEndAtomIdx() << " ";
      if (bond->GetBO()==5 && bcurrentState[i] == DOUBLE) {
        if ( (bond->GetBeginAtom())->IsSulfur()
             && bond->GetEndAtom()->IsSulfur() ) {
          // no double bonds between aromatic sulfur atoms -- PR#1504089
          continue;
        }

        bond->SetKDouble();
        bond->SetBO(2);
        //std::cout << "double\n";
      }
      //else
      //std::cout << "single\n";
      //else if (bond->IsAromatic() && bond->GetBO() != 2)
      //  bond->SetBO(1);
    }

    return;
  }


  /////////////////////////////////////////////////////////////////////////////////////////
  //! \brief recursively assign single and double bonds according to the electronical state
  //! of the atoms of the current bond 
  int OBMol::expand_kekulize(OBAtom *atom1, OBAtom *atom2, std::vector<int> &currentState, std::vector<int> &initState, std::vector<int> &bcurrentState, std::vector<int> &binitState, std::vector<bool> &mark)
  {
    int done;
    int Idx1=atom1->GetIdx(), Idx2=atom2->GetIdx();
    OBBond *bond;
    std::vector<OBBond*>::iterator i;
    OBAtom *nbr;
    int natom;

    mark[Idx1]= true;
    bond = atom1->GetBond(atom2);
    int bIdx = bond->GetIdx();
  
    //cout << "assign bond state for atoms " << Idx1 << " and " << Idx2 << endl;
    if (currentState[Idx1] == 1 && currentState[Idx2] == 1) {
      currentState[Idx1]=0;
      currentState[Idx2]=0;
      // set bond to double
      //std::cout << "bond " << Idx1 << " " << Idx2 << " double\n";
      bcurrentState[bIdx]=DOUBLE;
    }
    else if ((currentState[Idx1] == 0 && currentState[Idx2] == 1) ||
             (currentState[Idx1] == 2 && currentState[Idx2] == 1) ||
             (currentState[Idx1] == 2 && currentState[Idx2] == 2)) {
      //std::cout << "bond " << Idx1 << " " << Idx2 << " single\n";
      // leave bond to single
    }  
    else if ((currentState[Idx1] == 1 && currentState[Idx2] == 0) ||
             (currentState[Idx1] == 1 && currentState[Idx2] == 2)) {
      mark[Idx1]=false;
      //std::cout << "bond " << Idx1 << " " << Idx2 << " error\n";
      return (0); // error
    }
    else if ((currentState[Idx1] == 0 && currentState[Idx2] == 0)
             || (currentState[Idx1] == 2 && currentState[Idx2] == 0)) { 
      //std::cout << "bond " << Idx1 << " " << Idx2 << " done\n";
      mark[Idx2]=true;
      return (1); //done
    }
    else if (currentState[Idx1] == 0 && currentState[Idx2] == 2) {
      currentState[Idx2]=0;
      //std::cout << "bond " << Idx1 << " " << Idx2 << " leave single\n";
      // leave bond to single
    }
    else {

      stringstream errorMsg;

      errorMsg << "unexpected state:" << "atom " << Idx1 << " " << currentState[Idx1] 
               << " atom " << Idx2 << " " << currentState[Idx2] << endl;
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obDebug);
      return(false);
    }

    //int c;
    //for(c=1; c < currentState.size(); ++c) {
    //cout << c << "\t";
    //}
    //cout << endl;
    //for(c=1; c < currentState.size(); ++c) { 
    //cout << currentState[c] << "\t";
    //}   
    //cout << endl;

    vector<int> previousState = currentState;   // Backup the atom
    vector<int> bpreviousState = bcurrentState; // and the bond states before expanding again

    bool return_false=false;
    // for each neighbor of atom 2 not already kekulized
    for (nbr = atom2->BeginNbrAtom(i);nbr;nbr = atom2->NextNbrAtom(i))
      {
        natom = nbr->GetIdx();
        if(initState[natom] == -1) //neighbor atom not in the cycle, try next one
          continue;            
        if ( !mark[natom] ) {
          done = expand_kekulize(atom2, nbr, currentState, initState, bcurrentState, binitState, mark);
          if ( !done )  // kekulize failed
            return_false =true;
          else 
            return_false =false;
        }
    
      }
    if (return_false) { // no good solution found 
      //cout << "return_false:no good solution\n" << endl;
      //cout << "reset state of " << Idx1 << " and " << Idx2 << " from "  << currentState[Idx1]
      //<< " " << currentState[Idx2] << " to ";
    
      // retrieve the states that might have been changed during kekulize expansion
      currentState = previousState;
    
    
      bcurrentState = bpreviousState;
    
    
      // reset the bond and the atom states
      if (bcurrentState[bIdx] == DOUBLE)
        currentState[Idx1]=initState[Idx1];
    
      currentState[Idx2]=initState[Idx2];
      bcurrentState[bIdx]=binitState[bIdx]; // should be always single
      mark[Idx2]=false;

      //cout << currentState[Idx1] << " " << currentState[Idx2] << endl;
      return (false);
    }
    // atom 2 cannot make any bond, should not have 1 electron 
    if (currentState[Idx2] == 1) {
      // currentState[Idx1]=initState[Idx1];
      //cout << "return true but " << Idx2 << " state = 1\n";
      mark[Idx1]=false;
      return (false);
    }
    else {
      // if we found a good solution, then the state of Idx2 may have shifted from 1 to 0 during the kekulization
      // If it is the case, we should check if there is a remaining unmarked neighbor because it is possible
      // that kekulizing from this neighbor failed just because Idx2 was equal to 1

      if(previousState[Idx2] == 1) {
        // Since now Idx2 is equal to 0 because it kekulized well the kekulizing of the failed neighbor could be successful
        // If there is still an unmarked neighbor try to kekulize it again
        //mark[Idx2]=true;	
        return_false=false;
        //cout << "retry kekulizing from " << Idx2 << endl;
  
        for (nbr = atom2->BeginNbrAtom(i);nbr;nbr = atom2->NextNbrAtom(i))
          {
            natom = nbr->GetIdx();
            if(initState[natom] == -1) //neighbor atom not in the cycle, try next one
              continue;            
            if ( !mark[natom] ) {
              //cout << "atom " << natom << " not marked, expand again" << endl;
              done = expand_kekulize(atom2, nbr, currentState, initState, bcurrentState, binitState, mark);
              if ( !done )  // kekulize failed
                return_false =true;
              else 
                return_false =false;
            }
    
          }
  
        // if we cannot kekulize the remaining neighbor again then we have to return false
        // we do not have to reset the states because the kekulize will fail anyway
        if(return_false) {
          //cout << "rekekulize failed" << endl;
          return(false);
        }
        else {
          //cout << "rekekulized successful" << endl;
          return (true);
        }
	  
      }

    
      //cout << "return_true: good solution" << endl;
      return (true);
    }
  }

  //! Give the priority to give two electrons instead of 1
  int OBMol::getorden( OBAtom *atom) 
  {
    if ( atom->IsSulfur() ) return 1;
    if ( atom->IsOxygen() ) return 2;
    if ( atom->GetAtomicNum() == 34 || atom->GetAtomicNum() == 52 ) return 3;
    if ( atom->IsNitrogen() && atom->GetFormalCharge() == 0 && atom->GetValence() == 3) return 5;
    if ( atom->IsAmideNitrogen() ) return 4;
    if ( atom->IsNitrogen() && atom->GetFormalCharge() == -1) return 6;
    if ( atom->IsNitrogen() && atom->GetFormalCharge() == 0 && atom->IsInRingSize(5) ) return 7;
    if ( atom->IsNitrogen() && atom->GetFormalCharge() == 0 ) return 8;
    if ( atom->IsCarbon() && atom->GetFormalCharge() == -1) return 9;
    //if ( atom->IsCarbon() ) return 9;
  
    return (100); //no atom found
  }

  //! Recursively find the aromatic atoms with an aromatic bond to the current atom
  bool OBMol::expandcycle (OBAtom *atom, OBBitVec &avisit, OBAtom *, int)
  {
    OBAtom *nbr;
    //  OBBond *bond;
    std::vector<OBBond*>::iterator i;
    int natom;
    //for each neighbour atom test if it is in the aromatic ring
    for (nbr = atom->BeginNbrAtom(i);nbr;nbr = atom->NextNbrAtom(i))
      {
        natom = nbr->GetIdx();
        // if (!avisit[natom] && nbr->IsAromatic() && ((OBBond*) *i)->IsAromatic()) {
        if (!avisit[natom] && ((OBBond*) *i)->GetBO()==5 
            && ((OBBond*) *i)->IsInRing()) {
          avisit.SetBitOn(natom);
          expandcycle(nbr, avisit);
        }
      }

    return true;
  }

} // end namespace OpenBabel

//! \file kekulize.cpp
//! \brief Alternate algorithm to kekulize a molecule (OBMol::NewPerceiveKekuleBonds()).

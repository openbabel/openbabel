/* Fabien test for a new algorithm to kekulize molecule */

#include "mol.h"
#include "babelconfig.h"

#define SINGLE 1
#define DOUBLE 2

namespace OpenBabel {

void OBMol::kekulize() 
{
  OBAtom *atom;
  int i,j, n, de, minde;
  std::vector<OBAtom*> cycle;
  OBBitVec avisit,cvisit;
  avisit.Resize(NumAtoms()+1);
  cvisit.Resize(NumAtoms()+1);
  OBBond *bond;
  std::vector<OBEdgeBase*>::iterator bi;
  std::vector<int> electron;
  int BO;
  int sume, orden, bestorden, bestatom;
  // Init the kekulized bonds
  //for( i=0; i< NumBonds(); i++ ) {
  // bond = GetBond(i);
  // BO = bond->GetBO();
  // if (BO == 5 || BO == 1) SetKSingle();
  //  if (BO == 2) SetKDouble();
  //  if (BO == 3) SetKTriple();
  //}

  // Find all the groups of aromatic cycle
  for( i=1; i<= NumAtoms(); i++ ) {
    atom = GetAtom(i);
    if (atom->HasAromaticBond() && !cvisit[i]) { // is new aromatic atom of an aromatic cycle ?
      
      avisit.Clear();
      electron.clear();
      cycle.clear();
      
      avisit.SetBitOn(i);
      expandcycle (atom, avisit);
      //store the atoms of the cycle(s)
      for(j=1; j<= NumAtoms(); j++) {
	if ( avisit[j] ) {
	  atom = GetAtom(j);
	  cycle.push_back(atom);
	}
      }
      // At the begining each atom give one electron to the cycle
      for(j=0; j< cycle.size(); j++) {
	electron.push_back(1);
      }
      
      // remove one electron if the atom make a double bond out of the cycle
      sume =0;
      for(j=0; j< cycle.size(); j++) {
	atom = cycle[j];
	for(bond = atom->BeginBond(bi); bond; bond = atom->NextBond(bi)) {
	  if ( bond->IsDouble() )
	    electron[j] = 0;
	}
	// count the number of electrons
	sume += electron[j];
      }
      
      // find the ideal number of electrons according to the huckel 4n+2 rule
      minde=99;
      for (i=1; 1; i++) {
	n = 4 *i +2;
	de = n - sume;
	if (  de < minde )
	  minde=de;
	else if ( minde < 0 )
	  minde=de;
	else
	  break;
      }
      
      // if huckel rule not satisfied some atoms must give more electrons
      
      while ( minde != 0 ) {
	bestorden=99;
	for(j=0; j< cycle.size(); j++) {
	  if (electron[j] == 1) {
	    orden = getorden(cycle[j]);
	    if (orden < bestorden) {
	      bestorden = orden;
	      bestatom = j;
	    }
	  }
	}
	if (bestorden==99) {  // no electron giving atom found
	  std::cout << "Kekulize Warning: Huckel rule not satisfied for molecule " << GetTitle() << "\n";
	  break;             // Huckel rule cannot be satisfied
	}                    // try to kekulize anyway
	else {
 	  electron[bestatom] += 1;
	  minde--;
	}
      }

      // kekulize the cycle(s)
      start_kekulize(cycle,electron);

      // Set the kekulized cycle(s) as visited
      for(j=1; j<= NumAtoms(); j++) {
	if (avisit[j])
	  cvisit.SetBitOn(j);
      }

    }
  }
  // Double bond have been assigned, set the remaining aromatic bonds to single
  //std::cout << "Set not assigned single bonds\n"; 
  for(i=0;i <NumBonds(); i++) {
    bond = GetBond(i);     
    //std::cout << "bond " << bond->GetBeginAtomIdx() << " " << bond->GetEndAtomIdx() << " ";   
    if (bond->GetBO()==5 ) {
      bond->SetBO(1); 
      //std::cout << "single\n";
    }
    //else
    //  std::cout << "double\n";
  }

  return;
}

void OBMol::start_kekulize( std::vector <OBAtom*> &cycle, std::vector<int> &electron) {
  
  std::vector<int> initState;
  std::vector<int> currentState;
  std::vector<int> binitState;
  std::vector<int> bcurrentState;
  std::vector<bool> mark;
  int i;
  unsigned int Idx;
  OBAtom *atom, *atom2;
  OBBond *bond;

  //init the atom arrays
  for(i=0;i <NumAtoms()+1; i++) {
    initState.push_back(0);
    currentState.push_back(0);
    mark.push_back(false);
  }
  
  //init the bond arrays with single bonds
  for(i=0;i <NumBonds(); i++) {
    binitState.push_back(SINGLE);
    bcurrentState.push_back(SINGLE);
  }
  
  //set the electron number
  for( i=0; i< cycle.size(); i++) {
    atom = cycle[i];
    Idx =  atom->GetIdx();
    initState[Idx] = electron[i];
    currentState[Idx] = initState[Idx];
  }

  std::vector<OBEdgeBase*>::iterator b;
  OBAtom *nbr;
  
  bool second_pass=false;
  // for( i=1; i<= NumAtoms(); i++) {
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
  bool expand_successful;
  atom = cycle[0];
  for (nbr = atom->BeginNbrAtom(b);nbr;nbr = atom->NextNbrAtom(b)) {
    // std::cout << "Expand kekulize\n";
    expand_kekulize(atom,nbr,currentState,initState, bcurrentState,binitState, mark) ; 
    //Control that all the electron have been given to the cycle(s)
    expand_successful = true;
    for( i=0; i< cycle.size(); i++) {
      atom2 = cycle[i];
      Idx =  atom2->GetIdx();
      if (currentState[Idx] == 1)
	expand_successful=false;
    }
    if (expand_successful)
      break;
    else {
      for(i=0;i <NumAtoms()+1; i++) {
	currentState[i]=initState[i];
	mark[i]=false;
      }
      for(i=0;i <NumBonds(); i++) {
	bcurrentState[i]=binitState[i];
      }   
    }
  }
  if (!expand_successful) std::cout << "kekulize error for molecule " << GetTitle() << "\n"; 

  // Set the double bonds
  // std::cout << "Set double bonds\n";
  for(i=0;i <NumBonds(); i++) {
    bond = GetBond(i);    
    // std::cout << "bond " << bond->GetBeginAtomIdx() << " " << bond->GetEndAtomIdx() << " ";
    if (bond->GetBO()==5 && bcurrentState[i] == DOUBLE) {
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

int OBMol::expand_kekulize(OBAtom *atom1, OBAtom *atom2, std::vector<int> &currentState, std::vector<int> &initState, std::vector<int> &bcurrentState, std::vector<int> &binitState, std::vector<bool> &mark)
{
  int done;
  int Idx1=atom1->GetIdx(), Idx2=atom2->GetIdx();
  OBBond *bond;
  std::vector<OBEdgeBase*>::iterator i;
  OBAtom *nbr;
  int natom;

  mark[Idx1]= true;
  bond = atom1->GetBond(atom2);
  int bIdx = bond->GetIdx();
  
  if (currentState[Idx1] == 1 && currentState[Idx2] == 1) {
    currentState[Idx1]=0;
    currentState[Idx2]=0;
    // set bond to double
    //std::cout << "bond " << bond->GetBeginAtomIdx() << " " << bond->GetEndAtomIdx() << "double\n";
    bcurrentState[bIdx]=DOUBLE;
  }
  else if (currentState[Idx1] == 0 && currentState[Idx2] == 1 ||
	   currentState[Idx1] == 2 && currentState[Idx2] == 1 ||
	   currentState[Idx1] == 2 && currentState[Idx2] == 2) {
    //std::cout << "bond " << bond->GetBeginAtomIdx() << " " << bond->GetEndAtomIdx() << "single\n";
    // leave bond to single
  }  
  else if (currentState[Idx1] == 1 && currentState[Idx2] == 0 ||
      currentState[Idx1] == 1 && currentState[Idx2] == 2) {
    mark[Idx1]=false;
    return (0); // error
  }
  else if (currentState[Idx1] == 0 && currentState[Idx2] == 0
      || currentState[Idx1] == 2 && currentState[Idx2] == 0) { 
    return (1); //done
  }
  else if (currentState[Idx1] == 0 && currentState[Idx2] == 2) {
    currentState[Idx2]=0;
    // leave bond to single
  }

  bool return_false=false;
  // for each neighbor of atom 2 not already kekulized
  for (nbr = atom2->BeginNbrAtom(i);nbr;nbr = atom2->NextNbrAtom(i))
   {
     natom = nbr->GetIdx();
     if ( !mark[natom] ) {
       done = expand_kekulize(atom2, nbr, currentState, initState, bcurrentState, binitState, mark);
       if ( !done )  // kekulize failed
	 return_false =true;
       else 
	 return_false =false;
     }
    
   }
  if (return_false) { // no good solution found 
    if (bcurrentState[bIdx] == DOUBLE)
      currentState[Idx1]=initState[Idx1];
    
    currentState[Idx2]=initState[Idx2];
    bcurrentState[bIdx]=binitState[bIdx]; // should be always single}
    mark[Idx2]=false;
    return (false);
  }
  // atom 2 cannot make any bond, should not have 1 electron 
  if (currentState[Idx2] == 1) {
    // currentState[Idx1]=initState[Idx1];
    mark[Idx1]=false;
    return (false);
  }
  else
    return (true);
}

// Give the priority to give two electrons instead of 1
int OBMol::getorden( OBAtom *atom) 
{
  int orden;
  
  if ( atom->IsSulfur() ) return 1;
  if ( atom->IsOxygen() ) return 2;
  if ( atom->GetAtomicNum() == 34 ) return 3;
  if ( atom->IsNitrogen() && atom->GetFormalCharge() == 0 && atom->GetValence() == 3) return 5;
  if ( atom->IsAmideNitrogen() ) return 4;
  if ( atom->IsNitrogen() && atom->GetFormalCharge() == -1) return 6;
  if ( atom->IsNitrogen() && atom->GetFormalCharge() == 0 && atom->IsInRingSize(5) ) return 7;
  if ( atom->IsNitrogen() && atom->GetFormalCharge() == 0 ) return 8;
  if ( atom->IsCarbon() && atom->GetFormalCharge() == -1) return 9;
  //if ( atom->IsCarbon() ) return 9;
  
  return (100); //no atom found
}

// Recursively find the aromatic atoms with an aromatic bond to the current atom
void OBMol::expandcycle (OBAtom *atom, OBBitVec &avisit)
{
  OBAtom *nbr;
  OBBond *bond;
  std::vector<OBEdgeBase*>::iterator i;
  int natom;
  //for each neighbour atom test if it is in the aromatic ring
  for (nbr = atom->BeginNbrAtom(i);nbr;nbr = atom->NextNbrAtom(i))
   {
     natom = nbr->GetIdx();
     // if (!avisit[natom] && nbr->IsAromatic() && ((OBBond*) *i)->IsAromatic()) {
     if (!avisit[natom] && ((OBBond*) *i)->GetBO()==5) {
       avisit.SetBitOn(natom);
       expandcycle(nbr, avisit);
     }
   }
}

} // end namespace OpenBabel

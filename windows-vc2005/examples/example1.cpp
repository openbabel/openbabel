#include <iostream>
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/obiter.h>

using namespace std;
using namespace OpenBabel;

int main(int argc,char **argv){

if(argc<2){
cout << "Usage: ProgrameName InputFileName";
return 1;
}

ifstream ifs(argv[1]);
if(!ifs){
cout << "Cannot open input file";
return 1;
}

OpenBabel::OBConversion conv(&ifs, &cout);
if(conv.SetInAndOutFormats("mol","sdf")){
  OpenBabel::OBMol mol;
  if(conv.Read(&mol)){
    cout << "Mol.wt: "<<mol.GetMolWt()<<endl;//works even with implicit hydrogen
    mol.AddHydrogens(false,false); //ensure that hydrogens are all explicit
    cout << "No. of atoms: " << mol.NumAtoms()<<endl;
    cout << "No. of hvy atoms: " << mol.NumHvyAtoms() << endl;
    cout << "No. of bonds: " << mol.NumBonds() << endl;
    double exactMass = 0.0;
    FOR_ATOMS_OF_MOL(a, mol){
      cout << "iterating..."<<endl;
      exactMass += a->GetExactMass();
    }
    cout << "Exact Mass: " << exactMass << endl;
  }
  else
    cout << "can't conv.read()."<<endl;
}
else
  cout << "Incorrect in/out format." << endl;
return 0;
}

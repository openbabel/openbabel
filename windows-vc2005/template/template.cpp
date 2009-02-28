#include <iostream>
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/obiter.h>

using namespace std;
using namespace OpenBabel;

int main(int argc,char **argv){

if(argc<3){
cout << "Usage: ProgrameName InputFileName OutputFileName";
return 1;
}

ifstream ifs(argv[1]);
if(!ifs){
cout << "Cannot open input file";
return 1;
}

ofstream ofs(argv[2]);
if(!ofs){
cout << "Cannot open output file";
return 1;
}

ofstream fout("out.txt");
OpenBabel::OBConversion conv(&ifs, &ofs);
if(conv.SetInAndOutFormats("mol","sdf")){
OpenBabel::OBMol mol;
OpenBabel::OBBond *bond;
OpenBabel::OBAtom *atom;
if(conv.Read(&mol)){
fout << "Mol.wt: "<<mol.GetMolWt()<<endl;
fout << "No. of atoms: " << mol.NumAtoms()<<endl;
fout << "No. of hvy atoms: " << mol.NumHvyAtoms() << endl;
fout << "No. of bonds: " << mol.NumBonds() << endl;
double exactMass = 0.0;
FOR_ATOMS_OF_MOL(a, mol){
fout << "iterating..."<<endl;
exactMass += a->GetExactMass();
}
fout << "Exact Mass: " << exactMass << endl;
}
else
fout << "can't conv.read()."<<endl;
}
else
fout << "Incorrect in/out format." << endl;
return 0;
}

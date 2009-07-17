#include <iostream>
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/obiter.h>

using namespace std;
using namespace OpenBabel;

int main(int argc,char **argv){
  //OpenBabel::obErrorLog.SetOutputLevel(debug);
  OpenBabel::OBConversion conv;
  if (!conv.SetInFormat("smi")) {
    cout << "Cannot set input format" << endl;
    return 1;
  }
  OpenBabel::OBMol mol;
  string a;
  //conv.ReadFile(&mol, "446.mol");
  conv.ReadString(&mol, "O=c1[nH]c(N)ccn1");
   FOR_BFS_OF_MOL(a, mol)
      {
         // The variable a behaves like OBAtom* when used with -> and * but
         // but needs to be explicitly converted when appearing as a parameter
         // in a function call - use &*a
        cout << a->GetAtomicMass() << " " << a.CurrentDepth() << std::endl;;
  
      }

  if(conv.SetOutFormat("smi")){
    a = conv.WriteString(&mol);
    cout << "Regular smiles:\t" << endl << a << endl;
  }
  /*if (!conv.SetInFormat("mol2")) {
    cout << "Cannot set input format" << endl;
    return 1;
  }
  OpenBabel::OBMol mol;
  conv.ReadFile(&mol, "PR2691618_aa.mol2");
  if(conv.SetOutFormat("mol2")){
    cout << "Mol2 smiles:\t" << endl << conv.WriteString(&mol) << endl;
  }*/
  else
    cout << "Incorrect output format." << endl;
  return 0;
}

//Simple command line to test OBConversion when used in another program

#ifdef WIN32
#pragma warning (disable : 4786)
#endif

#include "babelconfig.h"
#include <conio.h> //for getch
#include "mol.h"
#include "OBConversion.h"
#include "Reaction.h"
#include "ctiter.h"

using namespace OpenBabel;

int main(int argc,char *argv[])
{

/* This example program is one which can use of the full OpenBabel API.

In order to read in a molecule, manipulate it and write it out:

Set up an istream and an ostream, to and from files or elsewhere.
(cin and cout are used in the example)
	
A two stage construction is used to allow error handling
if the format ID is not recognized. This is necessary now that the formats are dynamic
and errors are not caught at compile time.
OBConversion::Read() is a templated function so that any object derived from OBBase
can also be handled, in addition to OBMol, if the format routines are written appropriately.

An executable made from this file together is run with obdll.dll and obconv.dll and 
one or more format dlls, e.g. obformats.obf, in the same directory. (During its construction
it is linked with obdll.lib and obconv.lib to provide entry points in the dlls).

Alternatively, it could be compiled and linked with all the rest of OpenBabel code.
*/
	OBConversion conv(&cin,&cout); //The object that handles input/output

	if(!conv.SetInAndOutFormats("SMI","MOL"))
	{
		cerr << "Invalid format" << endl;
	}
	else
	{
		OBMol mol;

		cout << "Please enter SMILES" <<endl;
		if(!conv.Read(&mol))
			cerr << "Error interpreting this!" << endl;
		else
		{
			//...manipulate the OBMol object

			OBCTIter iter(mol.GetAtom(1));//This iterator accesses the connection table
			//tree, depth first, starting at the specified atom
			for(;iter;++iter)
			{
				cout << iter->GetType();
			}
			cout << endl;

			mol.AddHydrogens();
			cout << "\nThe MW of this molecule is " << mol.GetMolWt() << endl;
			
			// write to cout
			cout <<"Output format is " << conv.GetOutFormat()->Description() << endl;
			conv.Write(&mol);
		}
	}

	//Using a CReaction object
	if(!conv.SetInAndOutFormats("RXN","RXN"))
		cerr << "Invalid format" << endl;
	else
	{
		ifstream ifs;
		ifs.open("FC.rxn");
		if(!ifs)
			cerr << "couldn't open file FC.rxn" << endl;
		OBReaction React;
		if(conv.Read(&React, &ifs))
		{
			
			// you could do something useful with the OBReaction object here

			// write to cout
			cout <<"\nOutput format is " << conv.GetOutFormat()->Description() << endl << endl;
			conv.Write(&React,&cout);
		}
	}

	#ifdef _DEBUG
	 // keep window open
		cout << "Press any key to finish" <<endl;
		getch();
	#endif
	
	return 0;
}
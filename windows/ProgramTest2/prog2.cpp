/*
This example shows how to add automatic format conversion to an existing Windows program
using the DLL form of OpenBabel.

The existing program is assumed in this case to have a native format of SMILES.
The new code allows other formats to be used (recognized by their file extension)
if an appropriate *.obf file is present in the current directory. 
No code changes are necessary in the calling program when new formats added, or
when almost anything else in OpenBabel is changed.
If the format is not recognized the behaviour is as it was before the modification. 
*/

#include <iostream>
#include <sstream>
#include "obconversion.h" //(mol.h not necessary) 
#include <conio.h>
using namespace std;
using namespace OpenBabel;

int main(int argc,char *argv[])
{
	char filename[40];
	cout << "Type in a file name" << endl;
	//*************New code 
	cout << "Originally this had to be SMILES." << endl; 
	cout << endl << "But with the new code any of the following formats are accepted" << endl;
	OBConversion conv;
	Formatpos pos;
	OBFormat* pFormat;
	char* str=NULL;
	while(conv.GetNextFormat(pos,str,pFormat))
		cout << str <<endl;
	//*************
	cin >> filename ;//Original code

	ifstream ifs(filename); 
//	if(!strstr(filename,".smi")) Original code could only handle SMILES
//		cerr << Can't handle this type of file" << endl;
	if(!ifs)
	{
		cerr << "File not opened" << endl;
		return -1;
	}

	//*************New code 
	istream* pIn = &ifs; //just a different name for the input stream
	stringstream newstream;
	OBFormat* inFormat = conv.FormatFromExt(filename);
	OBFormat* outFormat = conv.FindFormat("SMI"); //The supposed native format
	if(inFormat && outFormat)
	{
		conv.SetInAndOutFormats(inFormat,outFormat);
		conv.Convert(pIn,&newstream);
		pIn=&newstream;
		cout << endl << "Here is the input stream converted to SMILES." <<endl ;
	}
	else 
		//error; new features not available; fallback to original functionality 
		cout << "The format conversion was not made. Carry on with the original input stream." <<endl;		
  //*************

	//...Carry on with original code using pIn
	
	cout << pIn->rdbuf() <<endl; //Show the contents of the stream

	#ifdef _DEBUG
	 // keep window open
		cout << "Press any key to finish" <<endl;
		getch();
	#endif

	return 0;
}
/*	This code would be linked with obconv.lib.
	At runtime the following dlls would be in the executable directory (or otherwise accessible):
	obconv.dll, obdll.dll, obformats.obf and possibly other *.obf format files.

*/

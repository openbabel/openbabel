/**********************************************************************
obtanimoto - Open Babel Tanimoto fingerprint calculation

Copyright (C) 2003-2004 Fabien Fontaine
Some portions Copyright (C) 2004-2005 Geoffrey R. Hutchison
 
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

#include "mol.h"
#include "fingerprint.h"
#include "obconversion.h"

using namespace std;
using namespace OpenBabel;

int main (int argc, char **argv)
{

    if (argc != 3)
    {
        cerr << "Usage: obtanimoto <query molecule> <fingerprints>\n";
        exit (1);
    }
    char *program_name = argv[0];
    char *FileIn = argv[1], *FptIn = argv[2];
    ifstream ifs;
    unsigned int i;
    // Find Input filetype
    OBConversion conv;
    OBFormat *format = conv.FormatFromExt(FileIn);
    if(!format || !conv.SetInAndOutFormats(format,format))
    {
        cerr << program_name << ": cannot read input format!" << endl;
        exit (-1);
    }
    //if (extab.CanWriteExtension(FileIn))
    //  outFileType = extab.FilenameToType(FileIn);
    //else
    //  {
    //    cerr << program_name << ": cannot write input format!" << endl;
    //    exit (-1);
    //  }


    OBMol mol;
    // Read the file
    ifs.open(FileIn);
    if (!ifs)
    {
        cerr << program_name << ": cannot read input molecule file!" << endl;
        exit (-1);
    }

    conv.Read(&mol,&ifs);

    ifs.close();
    // hash the molecule fragments
    fingerprint fpt(mol.GetTitle());
    fpt.HashMol(mol);

    ifs.open(FptIn);
    if (!ifs)
    {
        cerr << program_name << ": cannot read input fingerprint file!" << endl;
        exit (-1);
    }
    string line, name;
    OBBitVec fingerprint(1024);
    double t;
    while ( 1 )
    {
        if ( ! getline(ifs,line,'\n') )
            break;
        name= line;
        //erase '>'
        name.erase(0,1);
        // Get the fingerprint line
        if ( !getline(ifs,line,'\n'))
        {
            cerr << "Error: Unable to read second fingerprint line\n.";
            return 1;
        }
        if (line.size() != 1024)
        {
            cerr << "Error: bad fingerprint size: " << line.size();
            return 1;
        }
        // Set the fingerprint bits
        fingerprint.Clear();
        for(i=0; i<line.size(); i++)
        {
            if (line[i] == '1')
                fingerprint.SetBitOn(i);
        }
        t = Tanimoto(fpt.GetFpt(),fingerprint);
        cout << name << " " << t << endl;
    }
    return 0;
}

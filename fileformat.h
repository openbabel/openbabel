/**********************************************************************
Copyright (C) 2000 by Geoffrey Hutchison

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#ifndef FILEFORMAT_H
#define FILEFORMAT_H

#ifdef __sgi
#include <iostream.h>
#include <fstream.h>
#else
#include <iostream>
#include <fstream>
#endif

#include <algorithm>
#include <vector>
#include <string>

using namespace std;

namespace OpenEye {

class OEMol;

class OEFileFormat
{
 public:
  OEFileFormat(void) {}
  ~OEFileFormat(void) {}

  bool ReadMolecule(istream &,OEMol &, char *title="Untitled");
  bool WriteMolecule(ostream &,OEMol &, char *dimension="3D");

 private:
};

// Read Method prototypes
bool ReadAlchemy(istream &, OEMol &, char *defaultTitle="Untitled");
bool ReadAmberPrep(istream &, OEMol &, char *defaultTitle="Untitled");
bool ReadBallAndStick(istream &,OEMol &,char *defaultTitle="Untitled");
bool ReadBinary(istream&,OEMol&);
bool ReadBinary(unsigned char*,OEMol&,int);
bool ReadBinary(istream&,unsigned char **);
bool ReadBiosymCAR(istream &, OEMol &, char *defaultTitle="Untitled");
bool ReadBox(istream &, OEMol &,char *title="Untitled");
bool ReadBox(vector<string>&, OEMol &,char *title="Untitled");
bool ReadCaccrt(istream &, OEMol &, char *defaultTitle="Untitled");
bool ReadCCC(istream &ifs,OEMol &mol,char *title="Untitled");
bool ReadDMol(istream &ifs,OEMol &mol,char *title="Untitled");
bool ReadFeat(istream &, OEMol &, char *defaultTitle="Untitled");
bool ReadGAMESS(istream &, OEMol &, char *defaultTitle="Untitled");
bool ReadGhemical(istream &,OEMol&,char *defaultTitle="Untitled");
bool ReadHIN(istream &, OEMol &, char *defaultTitle="Untitled");
bool ReadJaguar(istream &, OEMol &, char *defaultTitle="Untitled");
bool ReadMacroModel(istream &,OEMol&,char *defaultTitle="Untitled");
bool ReadMol2(istream &,OEMol &,char *title="Untitled");
bool ReadMOPAC(istream &, OEMol &, char *defaultTitle="Untitled");
bool ReadMOPACCartesian(istream &, OEMol &, char *defaultTitle="Untitled");
bool ReadMPQC(istream &, OEMol &, char *defaultTitle="Untitled");
bool ReadNWChem(istream &,OEMol &,char *title="Untitled");
bool ReadPDB(istream &,OEMol &,char *title="Untitled");
bool ReadPDB(vector<string> &,OEMol &,char *title="Untitled");
bool ReadQChem(istream &, OEMol &, char *defaultTitle="Untitled");
bool ReadSDFile(istream &,OEMol &,char *title="Untitled");
bool ReadSmiles(istream &,OEMol &,char *title="Untitled");
bool ReadTerTermPDB(istream &,OEMol &,char *title="Untitled");
bool ReadUnichem(istream &, OEMol &, char *defaultTitle="Untitled");
bool ReadXYZ(istream &,OEMol &,char *defaultTitle="Untitled");
// Add yours here

// Write Method prototypes
bool SmiToMol(OEMol &mol,string &smi, char *title = "");
bool WriteAlchemy(ostream &, OEMol &);
bool WriteBallAndStick(ostream &,OEMol &);
bool WriteBGF(ostream &, OEMol &);
bool WriteBinary(ostream&,OEMol&);
bool WriteBinary(unsigned char*,int&,OEMol&);
bool WriteBox(ostream &ofs,OEMol &mol,float margin);
bool WriteCaccrt(ostream &,OEMol &);
bool WriteCacaoInternal(ostream &,OEMol &);
bool WriteCache(ostream &,OEMol &);
bool WriteChemDraw(ostream &,OEMol &);
bool WriteCSR(ostream &,OEMol &);
bool WriteCSSR(ostream &,OEMol &);
bool WriteDelphiPDB(ostream&,OEMol&);
bool WriteDMol(ostream &,OEMol &);
bool WriteFeat(ostream &,OEMol &);
bool WriteFenskeZmat(ostream &,OEMol &);
bool WriteFixFile(ostream&,OEMol&);
bool WriteGAMESS(ostream &,OEMol &);
bool WriteGhemical(ostream &,OEMol &); 
bool WriteGromos96A(ostream &,OEMol &);
bool WriteGromos96N(ostream &,OEMol &);
bool WriteGaussianCart(ostream &,OEMol &);
bool WriteHIN(ostream &, OEMol &);
bool WriteJaguar(ostream &, OEMol &);
bool WriteMacroModel(ostream &,OEMol &);
bool WriteMol2(ostream &,OEMol &,char *dimension="3D");
bool WriteMOPACCartesian(ostream &,OEMol &);
bool WriteNWChem(ostream &,OEMol&);
bool WritePDB(ostream &, OEMol&);
bool WriteQChem(ostream &,OEMol &);
bool WriteReport(ostream &,OEMol &);
bool WriteSDFile(ostream &,OEMol &,char *dimension="3D");
bool WriteSmiles(ostream &,OEMol &,char *title=NULL);
//test
bool WriteSmiOrderedMol2(ostream &,OEMol &,char *dimension="3D");
//test
bool WriteTheSmiles(OEMol&,char*);
bool WriteTinker(ostream &,OEMol &);
bool WriteTitles(ostream &,OEMol &);
bool WriteUnichem(ostream &,OEMol &);
bool WriteXED(ostream &,OEMol &);
bool WriteXYZ(ostream &,OEMol &);
// Add yours here

}

#endif

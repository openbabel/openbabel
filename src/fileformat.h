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

namespace OpenBabel {

class OBMol;

class OBFileFormat
{
 public:
  OBFileFormat(void) {}
  ~OBFileFormat(void) {}

  bool ReadMolecule(istream &,OBMol &, char *title="Untitled");
  bool WriteMolecule(ostream &,OBMol &, char *dimension="3D");

 private:
};

// Read Method prototypes
bool ReadAlchemy(istream &, OBMol &, char *defaultTitle="Untitled");
bool ReadAmberPrep(istream &, OBMol &, char *defaultTitle="Untitled");
bool ReadBallAndStick(istream &,OBMol &,char *defaultTitle="Untitled");
bool ReadBinary(istream&,OBMol&);
bool ReadBinary(unsigned char*,OBMol&,int);
bool ReadBinary(istream&,unsigned char **);
bool ReadBiosymCAR(istream &, OBMol &, char *defaultTitle="Untitled");
bool ReadBox(istream &, OBMol &,char *title="Untitled");
bool ReadBox(vector<string>&, OBMol &,char *title="Untitled");
bool ReadCaccrt(istream &, OBMol &, char *defaultTitle="Untitled");
bool ReadCCC(istream &ifs,OBMol &mol,char *title="Untitled");
bool ReadChem3d1(istream &ifs,OBMol &mol,char *title="Untitled");
bool ReadChem3d2(istream &ifs,OBMol &mol,char *title="Untitled");
bool ReadDMol(istream &ifs,OBMol &mol,char *title="Untitled");
bool ReadFeat(istream &, OBMol &, char *defaultTitle="Untitled");
bool ReadGAMESS(istream &, OBMol &, char *defaultTitle="Untitled");
bool ReadGhemical(istream &,OBMol&,char *defaultTitle="Untitled");
bool ReadHIN(istream &, OBMol &, char *defaultTitle="Untitled");
bool ReadJaguar(istream &, OBMol &, char *defaultTitle="Untitled");
bool ReadMacroModel(istream &,OBMol&,char *defaultTitle="Untitled");
bool ReadMmads(istream &,OBMol&,char *defaultTitle="Untitled");
bool ReadMol2(istream &,OBMol &,char *title="Untitled");
bool ReadMOPAC(istream &, OBMol &, char *defaultTitle="Untitled");
bool ReadMOPACCartesian(istream &, OBMol &, char *defaultTitle="Untitled");
bool ReadMPQC(istream &, OBMol &, char *defaultTitle="Untitled");
bool ReadNWChem(istream &,OBMol &,char *title="Untitled");
bool ReadPDB(istream &,OBMol &,char *title="Untitled");
bool ReadPDB(vector<string> &,OBMol &,char *title="Untitled");
bool ReadQChem(istream &, OBMol &, char *defaultTitle="Untitled");
bool ReadSDFile(istream &,OBMol &,char *title="Untitled");
bool ReadSmiles(istream &,OBMol &,char *title="Untitled");
bool ReadTerTermPDB(istream &,OBMol &,char *title="Untitled");
bool ReadUnichem(istream &, OBMol &, char *defaultTitle="Untitled");
bool ReadViewMol(istream &, OBMol &, char *defaultTitle="Untitled");
bool ReadXYZ(istream &,OBMol &,char *defaultTitle="Untitled");
// Add yours here

// Write Method prototypes
bool SmiToMol(OBMol &mol,string &smi, char *title = "");
bool WriteAlchemy(ostream &, OBMol &);
bool WriteBallAndStick(ostream &,OBMol &);
bool WriteBGF(ostream &, OBMol &);
bool WriteBinary(ostream&,OBMol&);
bool WriteBinary(unsigned char*,int&,OBMol&);
bool WriteBox(ostream &ofs,OBMol &mol,float margin);
bool WriteCaccrt(ostream &,OBMol &);
bool WriteCacaoInternal(ostream &,OBMol &);
bool WriteCache(ostream &,OBMol &);
bool WriteChem3d1(ostream &,OBMol &);
bool WriteChem3d2(ostream &,OBMol &);
bool WriteChemDraw(ostream &,OBMol &);
bool WriteCSR(ostream &,OBMol &);
bool WriteCSSR(ostream &,OBMol &);
bool WriteDelphiPDB(ostream&,OBMol&);
bool WriteDMol(ostream &,OBMol &);
bool WriteFeat(ostream &,OBMol &);
bool WriteFenskeZmat(ostream &,OBMol &);
bool WriteFixFile(ostream&,OBMol&);
bool WriteGAMESS(ostream &,OBMol &);
bool WriteGhemical(ostream &,OBMol &); 
bool WriteGromos96A(ostream &,OBMol &);
bool WriteGromos96N(ostream &,OBMol &);
bool WriteGaussianCart(ostream &,OBMol &);
bool WriteHIN(ostream &, OBMol &);
bool WriteJaguar(ostream &, OBMol &);
bool WriteMacroModel(ostream &,OBMol &);
bool WriteMmads(ostream &,OBMol &);
bool WriteMol2(ostream &,OBMol &,char *dimension="3D");
bool WriteMOPACCartesian(ostream &,OBMol &);
bool WriteNWChem(ostream &,OBMol&);
bool WritePDB(ostream &, OBMol&);
bool WriteQChem(ostream &,OBMol &);
bool WriteReport(ostream &,OBMol &);
bool WriteSDFile(ostream &,OBMol &,char *dimension="3D");
bool WriteSmiles(ostream &,OBMol &,char *title=NULL);
//test
bool WriteSmiOrderedMol2(ostream &,OBMol &,char *dimension="3D");
//test
bool WriteTheSmiles(OBMol&,char*);
bool WriteTinker(ostream &,OBMol &);
bool WriteTitles(ostream &,OBMol &);
bool WriteUnichem(ostream &,OBMol &);
bool WriteViewMol(ostream &,OBMol &);
bool WriteXED(ostream &,OBMol &);
bool WriteXYZ(ostream &,OBMol &);
// Add yours here

}

#endif

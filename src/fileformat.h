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

  bool ReadMolecule(std::istream &,OBMol &, char *title="Untitled");
  bool WriteMolecule(std::ostream &,OBMol &, char *dimension="3D", char *options="");

 private:
};

// Read Method prototypes
bool ReadAlchemy(std::istream &, OBMol &, char *title="Untitled");
bool ReadAmberPrep(std::istream &, OBMol &, char *title="Untitled");
bool ReadBallAndStick(std::istream &,OBMol &,char *title="Untitled");
bool ReadBinary(std::istream&,OBMol&);
bool ReadBinary(unsigned char*,OBMol&,int);
bool ReadBinary(std::istream&,unsigned char **);
bool ReadBiosymCAR(std::istream &, OBMol &, char *title="Untitled");
bool ReadBox(std::istream &, OBMol &,char *title="Untitled");
bool ReadBox(std::vector<std::string>&, OBMol &,char *title="Untitled");
bool ReadCaccrt(std::istream &, OBMol &, char *title="Untitled");
bool ReadCCC(std::istream &,OBMol &,char *title="Untitled");
bool ReadChem3d1(std::istream &,OBMol &,char *title="Untitled");
bool ReadChem3d2(std::istream &,OBMol &,char *title="Untitled");
bool ReadCML(std::istream &,OBMol &,char *title="Untitled");
bool ReadDMol(std::istream &,OBMol &,char *title="Untitled");
bool ReadFeat(std::istream &, OBMol &, char *title="Untitled");
bool ReadGAMESS(std::istream &, OBMol &, char *title="Untitled");
bool ReadGhemical(std::istream &,OBMol&,char *title="Untitled");
bool ReadHIN(std::istream &, OBMol &, char *title="Untitled");
bool ReadJaguar(std::istream &, OBMol &, char *title="Untitled");
bool ReadMacroModel(std::istream &,OBMol&,char *title="Untitled");
bool ReadMmads(std::istream &,OBMol&,char *title="Untitled");
bool ReadMol2(std::istream &,OBMol &,char *title="Untitled");
bool ReadMOPAC(std::istream &, OBMol &, char *title="Untitled");
bool ReadMOPACCartesian(std::istream &, OBMol &, char *title="Untitled");
bool ReadMPQC(std::istream &, OBMol &, char *title="Untitled");
bool ReadNWChem(std::istream &,OBMol &,char *title="Untitled");
bool ReadPDB(std::istream &,OBMol &,char *title="Untitled");
bool ReadPDB(std::vector<std::string> &,OBMol &,char *title="Untitled");
bool ReadQChem(std::istream &, OBMol &, char *title="Untitled");
bool ReadSDFile(std::istream &,OBMol &,char *title="Untitled");
bool ReadSmiles(std::istream &,OBMol &,char *title="Untitled");
bool ReadTerTermPDB(std::istream &,OBMol &,char *title="Untitled");
bool ReadUnichem(std::istream &, OBMol &, char *title="Untitled");
bool ReadViewMol(std::istream &, OBMol &, char *title="Untitled");
bool ReadXYZ(std::istream &,OBMol &,char *title="Untitled");
// Add yours here

// Write Method prototypes
bool SmiToMol(OBMol &,std::string &smi, char *title = "");
bool WriteAlchemy(std::ostream &, OBMol &);
bool WriteBallAndStick(std::ostream &,OBMol &);
bool WriteBGF(std::ostream &, OBMol &);
bool WriteBinary(std::ostream&,OBMol&);
bool WriteBinary(unsigned char*,int&,OBMol&);
bool WriteBox(std::ostream &,OBMol &,float margin);
bool WriteCaccrt(std::ostream &,OBMol &);
bool WriteCacaoInternal(std::ostream &,OBMol &);
bool WriteCache(std::ostream &,OBMol &);
bool WriteChem3d1(std::ostream &,OBMol &);
bool WriteChem3d2(std::ostream &,OBMol &);
bool WriteChemDraw(std::ostream &,OBMol &);
bool WriteCML(std::ostream &,OBMol &, char *dimension="3D", char *options="");
bool WriteCSR(std::ostream &,OBMol &);
bool WriteCSSR(std::ostream &,OBMol &);
bool WriteDelphiPDB(std::ostream&,OBMol&);
bool WriteDMol(std::ostream &,OBMol &);
bool WriteFeat(std::ostream &,OBMol &);
bool WriteFenskeZmat(std::ostream &,OBMol &);
bool WriteFixFile(std::ostream&,OBMol&);
bool WriteGAMESS(std::ostream &,OBMol &);
bool WriteGhemical(std::ostream &,OBMol &); 
bool WriteGromos96A(std::ostream &,OBMol &);
bool WriteGromos96N(std::ostream &,OBMol &);
bool WriteGaussianCart(std::ostream &,OBMol &);
bool WriteHIN(std::ostream &, OBMol &);
bool WriteJaguar(std::ostream &, OBMol &);
bool WriteMacroModel(std::ostream &,OBMol &);
bool WriteMmads(std::ostream &,OBMol &);
bool WriteMol2(std::ostream &,OBMol &,char *dimension="3D");
bool WriteMOPACCartesian(std::ostream &,OBMol &);
bool WriteNWChem(std::ostream &,OBMol&);
bool WritePDB(std::ostream &, OBMol&);
bool WriteQChem(std::ostream &,OBMol &);
bool WriteReport(std::ostream &,OBMol &);
bool WriteSDFile(std::ostream &,OBMol &,char *dimension="3D");
bool WriteSmiles(std::ostream &,OBMol &,char *title=NULL);
//test
bool WriteSmiOrderedMol2(std::ostream &,OBMol &,char *dimension="3D");
//test
bool WriteTheSmiles(OBMol&,char*);
bool WriteTinker(std::ostream &,OBMol &);
bool WriteTitles(std::ostream &,OBMol &);
bool WriteUnichem(std::ostream &,OBMol &);
bool WriteViewMol(std::ostream &,OBMol &);
bool WriteXED(std::ostream &,OBMol &);
bool WriteXYZ(std::ostream &,OBMol &);
// Add yours here

}

#endif

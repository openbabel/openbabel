/**********************************************************************
fileformat.h - Read and write file formats.

Copyright (C) 2000-2003 by Geoffrey R. Hutchison

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

#ifndef OB_FILEFORMAT_H
#define OB_FILEFORMAT_H

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <string>

namespace OpenBabel {

class OBMol;

//! \brief Central read/write file format class
//! Handles "dispatching" to the various read/write translation methods
//!  for given io_type and options
class OBFileFormat
{
 public:
  //! Read a molecule from the input stream--calls a Read method by the io_type
  static bool ReadMolecule(std::istream&,OBMol&,const char *title="Untitled");
  //! Write a molecule to the output stream--calls a Write method by the io_type
  static bool WriteMolecule(std::ostream&,OBMol&,const char *dimension="3D",
			    const char *options="");

 private:
};

// Read Method prototypes
bool ReadAlchemy(std::istream &, OBMol &, const char *title="Untitled");
bool ReadAmberPrep(std::istream &, OBMol &, const char *title="Untitled");
bool ReadBallAndStick(std::istream &,OBMol &,const char *title="Untitled");
bool ReadBGF(std::istream &,OBMol &,const char *title="Untitled");
bool ReadBinary(std::istream&,OBMol&);
bool ReadBinary(unsigned char *,OBMol&, int);
//bool ReadBinary(istream&,unsigned char **);
bool ReadBiosymCAR(std::istream &, OBMol &, const char *title="Untitled");
bool ReadBox(std::istream &, OBMol &,const char *title="Untitled");
bool ReadBox(std::vector<std::string>&, OBMol &,const char *title="Untitled");
bool ReadCaccrt(std::istream &, OBMol &, const char *title="Untitled");
bool ReadCCC(std::istream &,OBMol &,const char *title="Untitled");
bool ReadChem3d1(std::istream &,OBMol &,const char *title="Untitled");
bool ReadChem3d2(std::istream &,OBMol &,const char *title="Untitled");
bool ReadCML(std::istream &,OBMol &,const char *title="Untitled");
bool ReadCRK2D(std::istream &,OBMol &,const char *title="Untitled");
bool ReadCRK3D(std::istream &,OBMol &,const char *title="Untitled");
bool ReadDMol(std::istream &,OBMol &,const char *title="Untitled");
bool ReadFeat(std::istream &, OBMol &, const char *title="Untitled");
bool ReadGAMESS(std::istream &, OBMol &, const char *title="Untitled");
bool ReadGhemical(std::istream &,OBMol&,const char *title="Untitled");
bool ReadHIN(std::istream &, OBMol &, const char *title="Untitled");
bool ReadJaguar(std::istream &, OBMol &, const char *title="Untitled");
bool ReadMacroModel(std::istream &,OBMol&,const char *title="Untitled");
bool ReadMmads(std::istream &,OBMol&,const char *title="Untitled");
bool ReadMol2(std::istream &,OBMol &,const char *title="Untitled");
bool ReadMOPAC(std::istream &, OBMol &, const char *title="Untitled");
bool ReadMOPACCartesian(std::istream &, OBMol &, const char *title="Untitled");
bool ReadMPQC(std::istream &, OBMol &, const char *title="Untitled");
bool ReadNWChem(std::istream &,OBMol &,const char *title="Untitled");
bool ReadPDB(std::istream &,OBMol &,const char *title="Untitled");
bool ReadPDB(std::vector<std::string> &,OBMol &,const char *title="Untitled");
bool ReadPQS(std::istream &, OBMol &, const char *title="Untitled");
bool ReadQChem(std::istream &, OBMol &, const char *title="Untitled");
bool ReadSDFile(std::istream &,OBMol &,const char *title="Untitled");
bool ReadShelX(std::istream &,OBMol &,const char *title="Untitled");
bool ReadSmiles(std::istream &,OBMol &,const char *title="Untitled");
bool ReadTerTermPDB(std::istream &,OBMol &,const char *title="Untitled");
bool ReadUnichem(std::istream &, OBMol &, const char *title="Untitled");
bool ReadViewMol(std::istream &, OBMol &, const char *title="Untitled");
bool ReadXYZ(std::istream &,OBMol &,const char *title="Untitled");
// Add yours here

// Write Method prototypes
bool SmiToMol(OBMol &,std::string &smi, const char *title = "");
bool WriteAlchemy(std::ostream &, OBMol &);
bool WriteBallAndStick(std::ostream &,OBMol &);
bool WriteBGF(std::ostream &, OBMol &);
bool WriteBinary(std::ostream&,OBMol&);
bool WriteBinary(std::string&,int&,OBMol&); 
bool WriteBox(std::ostream &,OBMol &,double margin);
bool WriteCaccrt(std::ostream &,OBMol &);
bool WriteCacaoInternal(std::ostream &,OBMol &);
bool WriteCache(std::ostream &,OBMol &);
bool WriteChem3d1(std::ostream &,OBMol &);
bool WriteChem3d2(std::ostream &,OBMol &);
bool WriteChemDraw(std::ostream &,OBMol &);
bool WriteCHT(std::ostream &,OBMol &);
bool WriteCML(std::ostream &,OBMol &,
	      const char *dimension="3D", const char *options="");
bool WriteCRK2D(std::ostream &,OBMol &);
bool WriteCRK3D(std::ostream &,OBMol &);
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
bool WriteMol2(std::ostream &,OBMol &,const char *dimension="3D");
bool WriteMOPACCartesian(std::ostream &,OBMol &);
bool WriteNWChem(std::ostream &,OBMol&);
bool WritePDB(std::ostream &, OBMol&);
bool WritePovray(std::ostream &, OBMol&, const char*outputName="Output");
bool WritePQS(std::ostream &, OBMol&);
bool WriteQChem(std::ostream &,OBMol &);
bool WriteReport(std::ostream &,OBMol &);
bool WriteSDFile(std::ostream &,OBMol &,const char *dimension="3D");
bool WriteSmiles(std::ostream &,OBMol &,const char *title=NULL);
bool WriteSmiOrderedMol2(std::ostream &,OBMol &,const char *dimension="3D");
bool WriteTheSmiles(OBMol&,char*);
bool WriteTheSmiles(OBMol&,std::string&);
bool WriteTinker(std::ostream &,OBMol &);
bool WriteTitles(std::ostream &,OBMol &);
bool WriteUnichem(std::ostream &,OBMol &);
bool WriteViewMol(std::ostream &,OBMol &);
bool WriteXED(std::ostream &,OBMol &);
bool WriteXYZ(std::ostream &,OBMol &);
bool WriteZINDO(std::ostream &, OBMol &);
// Add yours here

}

#endif

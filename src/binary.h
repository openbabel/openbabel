/**********************************************************************
Copyright (C) 1998-2000 by OpenEye Scientific Software, Inc.

Some portions Copyright (c) 2001-2002 by Geoffrey R. Hutchison

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

#ifndef OB_BINARY_H
#define OB_BINARY_H

#include <map>

#include "mol.h"
#include "rotor.h"
#include "binary_io.h"
#include "generic.h"

namespace OpenBabel 
{

class OBRotamerList : public OBGenericData
{
  unsigned int _NBaseCoords;
  std::vector<float*> _c;
  std::vector<std::vector<float> > _vres;
  std::vector<unsigned char*> _vrotamer;
  std::vector<std::pair<OBAtom**,std::vector<int> > > _vrotor;

  //Invoking the default copy constructor or
  //assignment operator would cause a serious
  //memory leak.  MM 4/20/01 
  OBRotamerList(const OBRotamerList &cp) {};
  OBRotamerList& operator=(const OBRotamerList& cp) {return *this;}
public:
  OBRotamerList() {_NBaseCoords=0; _type=obRotamerList; _attr="RotamerList";}
  ~OBRotamerList();
  void Setup(OBMol&,OBRotorList&);
  void Setup(OBMol&,unsigned char*,int);
  int NumRotors() {return((_vrotor.empty())?0:_vrotor.size());}
  int NumRotamers() {return((_vrotamer.empty())?0:_vrotamer.size());}
  void AddRotamer(float*);
  void AddRotamer(int *arr);
  void AddRotamer(unsigned char *arr);
  void AddRotamers(unsigned char*,int);
  void GetReferenceArray(unsigned char*);
  void ExpandConformerList(OBMol&,std::vector<float*>&);
  std::vector<unsigned char*>::iterator BeginRotamer() {return(_vrotamer.begin());}
  std::vector<unsigned char*>::iterator EndRotamer() {return(_vrotamer.end());}

  //Support for internal storage of base coordinate sets that rotamers operate on
  //Added by MM 4/19/01
      //I cannot find any reference to _c in any of the member functions of OBRotamerList
      //I'm canibalizing it for internal storage of the base coordinate sets, which I'm
      //assuming was the original intent that was simply never coded.

      //Create a conformer list using the internal base set of coordinates
      std::vector<float*> CreateConformerList(OBMol& mol);

      //Copies the mol's conformers (the coordinates, NOT the pointers) into the object as base coordinates
      void SetBaseCoordinateSets(OBMol& mol) {SetBaseCoordinateSets(mol.GetConformers(),mol.NumAtoms());} 

      //Copies the coordinates in bc, NOT the pointers, into the object
      void SetBaseCoordinateSets(std::vector<float*> bc, unsigned int N); 

      //Returns true if the internal base coordinates are set
      unsigned int NumBaseCoordinateSets() {return _c.size();}

      //Get a pointer to a specific base pointer
      float *GetBaseCoordinateSet(unsigned int i) {return (i<_c.size()) ? _c[i] : NULL;}

      unsigned int NumAtoms() {return _NBaseCoords;}

};

class OBBinaryDBase
{
  std::ifstream _ifs;
  std::vector<std::streampos> _vpos;
 public:
  OBBinaryDBase(char*);
  OBBinaryDBase(std::string&);
  int  Size();
  void GetMolecule(OBMol&,int);
};

int Swab(int);

}

#endif // OB_BINARY_H

/**********************************************************************
obgen.cpp - test program for SMILES 3D coordinate generation
          - using internal coordinates

Copyright (C) 2006 Tim Vandermeersch
Some portions Copyright (C) 2006 Geoffrey R. Hutchison
 
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

// used to set import/export for Cygwin DLLs
#ifdef WIN32
#define USING_OBDLL
#endif

#include <openbabel/babelconfig.h>
#include <openbabel/base.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/forcefield.h>
#include <openbabel/fingerprint.h>
#include <unistd.h>

using namespace std;
using namespace OpenBabel;

// PROTOTYPES /////////////////////////////////////////////////////////////////
//int get_dst_nbr (OBAtom* atom);
//int get_ang_nbr (OBAtom* atom);
int get_nbr (OBAtom* atom, OBMol &mol, int level);


///////////////////////////////////////////////////////////////////////////////
//! \brief  Generate rough 3D coordinates for SMILES (or other 0D files)
//          based on bonding network, rings, atom types, etc.
//
int main(int argc,char **argv)
{
  char *program_name= argv[0];
  int c;
  string basename, filename = "";

  if (argc != 2)
    {
      string err = "Usage: ";
      err += program_name;
      err += " <filename>\n";
      ThrowError(err);
      exit(-1);
    }
  else
    {
      basename = filename = argv[1];
      size_t extPos = filename.rfind('.');

      if (extPos!= string::npos) {
        basename = filename.substr(0, extPos);
      }
    }

  // Find Input filetype
  OBConversion conv;
  OBFormat *format_in = conv.FormatFromExt(filename.c_str());
  OBFormat *format_out = conv.FindFormat("pdb");
    
  if (!format_in || !format_out || !conv.SetInAndOutFormats(format_in, format_out))
    {
      cerr << program_name << ": cannot read input/output format!" << endl;
      exit (-1);
    }

  ifstream ifs;
  ofstream ofs;

  // Read the file
  ifs.open(filename.c_str());
  if (!ifs)
    {
      cerr << program_name << ": cannot read input file!" << endl;
      exit (-1);
    }

  OBMol mol;

  for (c=1;;c++)
    {
      mol.Clear();
      if (!conv.Read(&mol, &ifs))
        break;
      if (mol.Empty())
        break;
      
      mol.AddHydrogens(false, true);

      OBAtom *atom, *nbr, *nbr2, *nbr3, *nbrs;
      vector<OBNodeBase*>::iterator i;
      vector<OBEdgeBase*>::iterator j;

      vector<OBInternalCoord*> internals;
      OBInternalCoord *coord;

      coord = new OBInternalCoord();
      internals.push_back(coord);
      mol.BeginModify();
      
      int torang;
      for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
        {
          coord = new OBInternalCoord();
          nbr = mol.GetAtom(get_nbr(atom, mol, 1));
          nbr2 = mol.GetAtom(get_nbr(atom, mol, 2));
          nbr3 = mol.GetAtom(get_nbr(atom, mol, 3));
          
	  if (nbr) {
            coord->_a = mol.GetAtom(get_nbr(atom, mol, 1));
            OBBond *bond;
            if ( (bond = mol.GetBond(atom, nbr)) ) {
              coord->_dst = bond->GetEquibLength();
            }
          }
          if (nbr2) {
            coord->_b = mol.GetAtom(get_nbr(atom, mol, 2));
            if (nbr->GetHyb() == 3)
              coord->_ang = 109;
            if (nbr->GetHyb() == 2)
	            coord->_ang = 120;
            if (nbr->GetHyb() == 1)
	            coord->_ang = 180;
 
          }
          if (nbr3) {
            double bestangle, angle, bestscore, score;
            int nbr_count;

            coord->_c = mol.GetAtom(get_nbr(atom, mol, 3));
	    coord->_tor = torang;
	    torang +=60;
		    
          }
            
          internals.push_back(coord);
        }

      InternalToCartesian(internals, mol);
      
    //string id;
    //OBForceField* pForceField=NULL;
    //while (OBForceField::GetNextForceField(id, pForceField)) {
    //  cout << id << " -- " << pForceField->Description() << endl;
    //}

    FOR_EACH(OBForceField, iter) {
      cout << iter.ID() << ' ' << iter.Description() << endl;
    }
    FOR_EACH(OBFingerprint, iter) {
      cout << iter.ID() << ' ' << iter.Description() << endl;
    }


    //OBForceField* pFF = OBForceField::FindForceField("MM2");
    OBForceField* pFF = OBForceField::FindForceField("Tripos");

	cout << "MM2:" << endl << endl;
	
	pFF->Setup(mol);
	
	cout << "    bond stretching            = " << pFF->E_Bond() << " kcal/mole" << endl;
	cout << "    angle bending              = " << pFF->E_Angle() << " kcal/mole" << endl;
	cout << "    stretch-bending            = " << pFF->E_StrBnd() << " kcal/mole" << endl;
	cout << "    torsional terms            = " << pFF->E_Torsion() << " kcal/mole" << endl;
	cout << "    out-of-plane bending       = " << pFF->E_OOP() << " kcal/mole" << endl;
	cout << "    Van der Waals interactions = " << pFF->E_VDW() << " kcal/mole" << endl;
	cout << "    Dipole-dipole interactions = " << pFF->E_Electrostatic() << " kcal/mole" << endl;
	
	cout << endl << "Total Energy = " << pFF->Energy() << " kcal/mole" << endl << endl;

	pFF->SteepestDescent(300);
        pFF->UpdateCoordinates(mol);
	
	cout << "    bond stretching            = " << pFF->E_Bond() << " kcal/mole" << endl;
	cout << "    angle bending              = " << pFF->E_Angle() << " kcal/mole" << endl;
	cout << "    stretch-bending            = " << pFF->E_StrBnd() << " kcal/mole" << endl;
	cout << "    torsional terms            = " << pFF->E_Torsion() << " kcal/mole" << endl;
	cout << "    out-of-plane bending       = " << pFF->E_OOP() << " kcal/mole" << endl;
	cout << "    Van der Waals interactions = " << pFF->E_VDW() << " kcal/mole" << endl;
	cout << "    Dipole-dipole interactions = " << pFF->E_Electrostatic() << " kcal/mole" << endl;
	
	cout << endl << "Total Energy = " << pFF->Energy() << " kcal/mole" << endl << endl;
	
      char FileOut[32];
      sprintf(FileOut, "obconfgen.pdb");
      ofs.open(FileOut);
      conv.Write(&mol, &ofs);
      ofs.close();
      conv.Write(&mol, &cout);
    } // end for loop

  return(1);
}

int get_nbr (OBAtom* atom, OBMol &mol, int level) {
  OBAtom *nbr,*nbr2,*nbr3;
  vector<OBEdgeBase*>::iterator i;
  
  if (level == 2)
    if (!get_nbr(atom, mol, 1)) return 0;
  if (level == 3)
    if (!get_nbr(atom, mol, 2)) return 0;
  
  // Find first neighboor
  FOR_NBORS_OF_ATOM(tmp, atom) {
    if (atom->GetIdx() > tmp->GetIdx()) {
      if (level == 1)
        return tmp->GetIdx();
	    else {
        nbr = mol.GetAtom(tmp->GetIdx());
        break;
	    }
    }
  }
  if (level == 1) return 0;
  
  // Find second neighboor
  FOR_NBORS_OF_ATOM(tmp, nbr) {
    if (atom->GetIdx() > tmp->GetIdx()) {
      if (level == 2)
        return tmp->GetIdx();
      else {
        nbr2 = mol.GetAtom(tmp->GetIdx());
        break;
	    }
    }
  }
  if (level == 2) return 0;
  
  // Find thirth neighboor
  FOR_NBORS_OF_ATOM(tmp, nbr2) {
    if ((atom->GetIdx() > tmp->GetIdx()) && (nbr->GetIdx() != tmp->GetIdx()))
      return tmp->GetIdx();
  }
  FOR_NBORS_OF_ATOM(tmp, nbr) {
    if ((atom->GetIdx() > tmp->GetIdx()) && (atom->GetIdx() != tmp->GetIdx()) && (nbr2->GetIdx() != tmp->GetIdx()))
      return tmp->GetIdx();
  }
  return 0;
}



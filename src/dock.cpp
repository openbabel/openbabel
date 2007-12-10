/**********************************************************************
dock.cpp - Docking Class.

Copyright (C) 2007 by Tim Vandermeersch <tim.vandermeersch@gmail.com>
 
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
#include <openbabel/babelconfig.h>

#include <openbabel/dock.h>

#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obiter.h>
#include <openbabel/math/matrix3x3.h>
#include <openbabel/rotamer.h>
#include <openbabel/rotor.h>
#include <openbabel/obconversion.h>
#include <openbabel/builder.h>

using namespace std;

namespace OpenBabel
{
  /** \class OBDock dock.h <openbabel/dock.h>
      \brief Docking Class

      The OBDock class is used for docking.

      Below is and example which explain the basics. 
      
      \code
      OBDock dock;
      dock.SpherePocket(mol, );
      \endcode
  **/


  OBDock::OBDock() 
  {
  }
 
  
  OBDock::~OBDock() 
  {
  }
  
  vector<int> OBDock::CreateLigandFromResidue(OBMol &mol, string &resname)
  {
    vector<int> ligand;

    FOR_RESIDUES_OF_MOL (res, mol) 
      if (resname == res->GetName()) 
	FOR_ATOMS_OF_RESIDUE (atom, &*res)
          ligand.push_back(atom->GetIdx());
    
    return ligand;
  }

  vector<int> OBDock::CreatePocketFromLigand(OBMol &mol, vector<int> &ligand , double radius)
  {
    vector<int> pocket;
    vector<int>::iterator i, j;
    bool is_added, is_in_ligand;

    FOR_ATOMS_OF_MOL (atom, mol) {
      for (i = ligand.begin(); i != ligand.end(); i++) {
        is_in_ligand = false;
	for (j = ligand.begin(); j != ligand.end(); j++)
	  if (atom->GetIdx() == (*j))
	    is_in_ligand = true;
        
	is_added = false;
	for (j = pocket.begin(); j != pocket.end(); j++)
	  if (atom->GetIdx() == (*j))
	    is_added = true;
	
	if (!is_in_ligand && !is_added) {
	  if (atom->GetDistance(mol.GetAtom(*i)) <= radius)
            pocket.push_back(atom->GetIdx());
	}
      }
    }
    
    return pocket;
  }


  vector<int> OBDock::CreatePocketFromResidue(OBMol &mol, string &resname , double radius)
  {
    vector<int> pocket;

    FOR_RESIDUES_OF_MOL (res, mol) 
      if (resname == res->GetName()) 
	FOR_ATOMS_OF_RESIDUE (resatom, &*res)
          FOR_ATOMS_OF_MOL (atom, mol) 
	    if (atom->GetResidue()->GetName() == resname)
	      continue;
            else 
	      if (atom->GetDistance(&*resatom) <= radius)
	        pocket.push_back(atom->GetIdx());
    
    return pocket;
  }
  
  void OBDock::GridDockInitialize(OBMol &mol, vector<int> &ligand, vector<int> &pocket, 
          int optimize_steps, double translate_step, double rotate_step)
  {
    _mol = mol;
    _ligand = ligand;
    _pocket = pocket;
    _optimize_steps = optimize_steps;
    _translate_step = translate_step;
    _rotate_step = rotate_step;
    
    _mol.SetConformer(0);
    _ligand_dimentions = GetDimentions(_ligand);
    _pocket_dimentions = GetDimentions(_pocket);

    cout << "_ligand_dimention = " << _ligand_dimentions << endl;
    cout << "_pocket_dimention = " << _pocket_dimentions << endl;
  }
  
  bool OBDock::GridDockNextPose()
  {
  
  }
  
  vector3 OBDock::GetDimentions(vector<int> &selection)
  {
    vector3 min, max = vector3(0.0, 0.0, 0.0);
    vector<int>::iterator i;
    
    FOR_ATOMS_OF_MOL (atom, _mol) 
      for (i = selection.begin(); i != selection.end(); i++) 
        if (atom->GetIdx() == (*i)) {
	  vector3 pos = atom->GetVector();
	  if (pos.x() < min.x())
	    min.SetX(pos.x());
	  if (pos.y() < min.y())
	    min.SetY(pos.y());
	  if (pos.z() < min.z())
	    min.SetZ(pos.z());
	  if (pos.x() > max.x())
	    max.SetX(pos.x());
	  if (pos.y() > max.y())
	    max.SetY(pos.y());
	  if (pos.z() > max.z())
	    max.SetZ(pos.z());
	} 
    
    return (max - min);
  }


} // end namespace OpenBabel


//! \file dock.cpp
//! \brief Handle OBDock class

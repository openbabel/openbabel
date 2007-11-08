/**********************************************************************
builder.cpp - Class to create structures.

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

#include <openbabel/builder.h>

#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obiter.h>
#include <openbabel/math/matrix3x3.h>
#include <openbabel/rotamer.h>
#include <openbabel/rotor.h>
#include <openbabel/obconversion.h>

/* OBBuilder::GetNewBondVector():
 * - is based on OBAtom::GetNewBondVector()
 * - but: when extending a long chain all the bonds are trans
 *
 * fragments.txt:
 * - fragments should be ordered from complex to simple
 *
 * */

using namespace std;

namespace OpenBabel
{
  OBBuilder::OBBuilder() 
  {
    // open data/fragments.txt
    string buffer2, subbuffer;
    ifstream ifs;
    buffer2 = BABEL_DATADIR;
    buffer2 += FILE_SEP_CHAR;
    subbuffer = buffer2;
    subbuffer += BABEL_VERSION;
    subbuffer += FILE_SEP_CHAR;
    subbuffer += "fragments.txt";
    buffer2 += "fragments.txt";

    ifs.open(subbuffer.c_str());
    if (!ifs) {
      ifs.open(buffer2.c_str());
    }
    if (!ifs) {
      cerr << "cannot find fragments.txt!" << endl;
      return;
    }

    char buffer[80];
    vector<string> vs;
    OBSmartsPattern *sp = NULL;
    vector<vector3> coords;
    bool nextfrag = false;
    int natoms, i = 0;
    while (ifs.getline(buffer, 80)) {
      tokenize(vs, buffer);
      
      if (!nextfrag) {
        nextfrag = true;
        natoms = atoi(vs[0].c_str());
      } else {
        if (i == 0) {
          coords.clear();
   	  sp = new OBSmartsPattern;
          if (!sp->Init(vs[1])) {
            delete sp;
            sp = NULL;
            obErrorLog.ThrowError(__FUNCTION__, " Could not parse SMARTS from contribution data file", obInfo);
          }
	} else if (i <= natoms) {
	  vector3 coord(atof(vs[0].c_str()), atof(vs[1].c_str()), atof(vs[2].c_str()));
	  coords.push_back(coord);
	} else {
	  _fragments.push_back(pair<OBSmartsPattern*, vector<vector3> > (sp, coords));
          nextfrag = false;	
	}
	i++;
      }
      
        
    }
    _fragments.push_back(pair<OBSmartsPattern*, vector<vector3> > (sp, coords));
    
  }
    
  
  OBBuilder::~OBBuilder() 
  {
    _fragments.clear();
  }

  vector3 OBBuilder::GetNewBondVector(OBMol &mol, OBAtom *atom)
  {
    vector3 bond1, bond2, bond3, v1, v2, newbond;
    
    bond1 = VZero;
    bond2 = VZero;
    bond3 = VZero;
    
    //cout << endl << "GetNewBondVector()" << endl;
    
    if (atom == NULL)
      return VZero;

    //cout << "atom : " << atom->GetIdx() << endl;
    //cout << "atom->Valence: " << atom->GetValence() << endl;
    //cout << "atom->Hyb: " << atom->GetHyb() << endl;


    if (atom->GetValence() == 0) {
      newbond = atom->GetVector() + VX * 1.5;
      return newbond;
    }

    if (atom->GetValence() == 1) {
      FOR_NBORS_OF_ATOM (nbr, atom) {
        bond1 = atom->GetVector() - nbr->GetVector();   
        
        FOR_NBORS_OF_ATOM (nbr2, &*nbr)
	  if (&*nbr2 != atom)
            bond2 = nbr->GetVector() - nbr2->GetVector();
      }

      //cout << "bond1: " << bond1 << endl;
      //cout << "bond2: " << bond2 << endl;
      
      bond1 = bond1.normalize();
      if (bond2 == VZero) {
        if (atom->GetHyb() == 1)
          newbond = bond1;
	if (atom->GetHyb() == 2)
	  newbond = bond1 + VY * tan(DEG_TO_RAD*60);
	  //newbond = bond1 + VY*1.7321;
	if (atom->GetHyb() == 3)
	  newbond = bond1 + VY * tan(DEG_TO_RAD*70.5);
	  //newbond = bond1 + VY*2.8239;
	
	newbond = newbond.normalize();
	newbond *= 1.5;
	newbond += atom->GetVector();
	return newbond;
      } else {
        v1 = cross(bond1, bond2);
	v2 = cross(bond1, v1);
        
	if (atom->GetHyb() == 1)
          newbond = bond1;
	if (atom->GetHyb() == 2)
	  newbond = bond1 - v2 * tan(DEG_TO_RAD*60);
	if (atom->GetHyb() == 3)
	  newbond = bond1 - v2 * tan(DEG_TO_RAD*70.5);
	
	newbond = newbond.normalize();
	newbond *= 1.5;
	newbond += atom->GetVector();
	return newbond;
      }
    }
    
    if (atom->GetValence() == 2) {
      FOR_NBORS_OF_ATOM (nbr, atom) {
        if (bond1 == VZero)
	  bond1 = atom->GetVector() - nbr->GetVector();
	else
	  bond2 = atom->GetVector() - nbr->GetVector();
      }

      v1 = bond1 + bond2;
      v1 = v1.normalize();
     
      cerr << "atom idx: " << atom->GetIdx() << ", hyb=" << atom->GetHyb() << endl;

      if (atom->GetHyb() == 2)
        newbond = v1;
      if (atom->GetHyb() == 3) {
        v2 = cross(bond1, bond2);
        v1 = v1.normalize();
	newbond = v2 + v1 * tan(DEG_TO_RAD*35.25);
      }
      
      newbond = newbond.normalize();
      newbond *= 1.5;
      newbond += atom->GetVector();
      return newbond;
    }
    
    if (atom->GetValence() == 3) {
      FOR_NBORS_OF_ATOM (nbr, atom) {
        if (bond1 == VZero)
          bond1 = atom->GetVector() - nbr->GetVector();
	else if (bond2 == VZero)
	  bond2 = atom->GetVector() - nbr->GetVector();
	else
	  bond3 = atom->GetVector() - nbr->GetVector();
      }
	  
      newbond = bond1 + bond2 + bond3;
      newbond = newbond.normalize();
      newbond *= 1.5;
      newbond += atom->GetVector();
      return newbond;
    }
  }

  bool OBBuilder::Build(OBMol &mol)
  {
    vector<int> visit(mol.NumAtoms()+1, 0);
    vector3 molvec, fragvec, rootvec, moldir, fragdir;
    vector<pair<OBSmartsPattern*, vector<vector3 > > >::iterator i;
    vector<vector<int> >::iterator j;
    vector<int>::iterator k, k2, k3;
    vector<vector3>::iterator l;
    vector<vector<int> > _mlist; //!< match list for atom typing
 
    OBMol tmpmol, fragment;

    //cout << "-------------------------------------------------------------" << endl;
    //cout << " OBBuilder::Build()" << endl;
    //cout << "-------------------------------------------------------------" << endl;

    tmpmol = mol;
    
    //FOR_ATOMS_OF_MOL (atoms, tmpmol)
    //  cout << "idx: " << atoms->GetIdx() << ", hyb: " << atoms->GetHyb() << endl;
     
    // delete all bonds in tmpmol
    while (tmpmol.NumBonds())
      tmpmol.DeleteBond(tmpmol.GetBond(0));
    
    tmpmol.SetHybridizationPerceived();
    //FOR_ATOMS_OF_MOL (atoms, tmpmol)
    //  cout << "idx: " << atoms->GetIdx() << ", hyb: " << atoms->GetHyb() << endl;
 
    // iterate over all atoms to place them in 3D space
    FOR_DFS_OF_MOL (a, mol) {
      cerr << "placing atom: " << a->GetIdx() << endl;
      if (visit[a->GetIdx()]) // continue if the atom is already added
        continue;

      // find an atom connected to the current atom that is already visited
      OBAtom *prev = NULL;
      FOR_NBORS_OF_ATOM (nbr, &*a) {
        if (visit[nbr->GetIdx()])
	  prev = &*nbr;
      }
     
      // get the position for the new atom, this is done with GetNewBondVector
      if (prev != NULL) {
        //cout << "prev: " << prev->GetIdx() << endl;

        molvec = GetNewBondVector(tmpmol, tmpmol.GetAtom(prev->GetIdx()));
	moldir = molvec - tmpmol.GetAtom(prev->GetIdx())->GetVector();
      } else {
        molvec = VX;
        moldir = VX;
      }
      
      /*
      if (prev != NULL)
        cout << "  prev: " << prev->GetIdx() << endl;
      else
        cout << "  prev: -" << endl;
      cout << "  molvec: " << molvec << endl;
      cout << "  moldir: " << moldir << endl;
      */
      
      // check if a is in a fragment
      bool foundfrag;
      for (i = _fragments.begin();i != _fragments.end();++i) {
	if (i->first->Match(mol)) { 
          _mlist = i->first->GetMapList();
          
	  for (j = _mlist.begin();j != _mlist.end();++j) {
            foundfrag = false;
            if (visit[(*j)[0]]) // the found match is already added
	      continue;

            int index, index2, counter = 0;
	    for (k = j->begin(); k != j->end(); ++k) { // for all atoms of the fragment
	      index = *k;
	      if (index == a->GetIdx())
	        foundfrag = true;
	    }

	    if (!foundfrag)
	      continue;
	    
	    for (k = j->begin(); k != j->end(); ++k) { // for all atoms of the fragment
	      index = *k;
	      //cout << "k=" << index << endl;

	      if (visit[index])
	        continue;
	      
	      visit[index] = 1; // set visit[] for all atoms of fragment

	      // set coordinates for atoms
	      OBAtom *atom = tmpmol.GetAtom(index);
              atom->SetVector(i->second[counter]);
	      counter++;
            }
	    
	    // add the bonds for the fragment
	    for (k = j->begin(); k != j->end(); ++k) {
	      index = *k;
	      OBAtom *atom1 = mol.GetAtom(index);
	        
	      for (k2 = j->begin(); k2 != j->end(); ++k2) {
	        index2 = *k2;
	        OBAtom *atom2 = mol.GetAtom(index2);
                OBBond *bond = atom1->GetBond(atom2);

                if (bond != NULL) 
                  tmpmol.AddBond(*bond);
	      }
            }
              
	    // 
	    // translate fragment to origin
	    //
	    for (k = j->begin(); k != j->end(); ++k) {
	      index = *k;
	      
	      if (index == a->GetIdx()) {
	        rootvec = tmpmol.GetAtom(a->GetIdx())->GetVector();
	    
	        //cout << "a position: " << rootvec << endl;
		  
	        for (k2 = j->begin(); k2 != j->end(); ++k2) {
	          index = *k2;

                  // translate the fragment
	          vector3 tmpvec = tmpmol.GetAtom(index)->GetVector();
	          tmpvec += molvec - rootvec;
                  tmpmol.GetAtom(index)->SetVector(tmpvec);
		  cerr << "tmpvec origin = " << tmpvec << endl;
     	        }
	      }
	    }
	
            if (prev != NULL) {
	    // 
	    // rotate
	    //  
	    matrix3x3 xymat, xzmat, yzmat;
	    double xyang, yzang, xzang, xyang1, yzang1, xzang1, xyang2, yzang2, xzang2;

  	    cerr << "----------------" << endl;
	    cerr << " R O T A T E" << endl;
	    cerr << "----------------" << endl;
	    cerr << "  moldir = " << moldir << endl;
		  
            fragvec = GetNewBondVector(tmpmol, tmpmol.GetAtom(a->GetIdx()));
	    fragdir = fragvec - tmpmol.GetAtom(a->GetIdx())->GetVector();
	    cerr << "  fragdir = " << fragdir << endl;
            
            cerr << "  --- XY plane ---" << endl;
            xyang = vectorAngle(vector3(moldir.x(), moldir.y(), 0.0), vector3(fragdir.x(), fragdir.y(), 0.0));
	    cerr << "cross: " << cross(vector3(moldir.x(), moldir.y(), 0.0), vector3(fragdir.x(), fragdir.y(), 0.0)) << endl;
	    if (cross(vector3(moldir.x(), moldir.y(), 0.0), vector3(fragdir.x(), fragdir.y(), 0.0)).z() > 0) {
	      cerr << "counterclockwise: " << xyang << endl;
	      xyang = 180 + xyang;
	    } else if (cross(vector3(moldir.x(), moldir.y(), 0.0), vector3(fragdir.x(), fragdir.y(), 0.0)).z() < 0) {
	      cerr << "clockwise: " << xyang << endl;
	      xyang = 180 - xyang;
	    } else
	      xyang = 0.0;
            xymat.SetupRotMat(0.0, 0.0, xyang); 
	    cerr << " xyang = "  << xyang << endl;
            
	    for (k = j->begin(); k != j->end(); ++k) {
	      index = *k;
	      vector3 tmpvec = tmpmol.GetAtom(index)->GetVector();
              tmpvec *= xymat; //apply the rotation
	      tmpmol.GetAtom(index)->SetVector(tmpvec);
	    }
 
            fragvec = GetNewBondVector(tmpmol, tmpmol.GetAtom(a->GetIdx()));
	    fragdir = fragvec - tmpmol.GetAtom(a->GetIdx())->GetVector();
	    cerr << "  fragdir = " << fragdir << endl;
 
            cerr << "  --- XZ plane ---" << endl;
	    xzang = vectorAngle(vector3(moldir.x(), moldir.z(), 0.0), vector3(fragdir.x(), fragdir.z(), 0.0));
	    cerr << "cross: " << cross(vector3(moldir.x(), moldir.z(), 0.0), vector3(fragdir.x(), fragdir.z(), 0.0)) << endl;
	    if (cross(vector3(moldir.x(), moldir.z(), 0.0), vector3(fragdir.x(), fragdir.z(), 0.0)).z() > 0) {
	      cerr << "counterclockwise: " << xzang << endl;
	      xzang = 180 - xzang;
	    } else if (cross(vector3(moldir.x(), moldir.z(), 0.0), vector3(fragdir.x(), fragdir.z(), 0.0)).z() < 0) {
	      cerr << "clockwise: " << xzang << endl;
	      xzang = 180 + xzang;
	    } else
	      xzang = 0.0;
	    xzmat.SetupRotMat(0.0, xzang, 0.0); 
	    cerr << " xzang = "  << xzang << endl;

	    for (k = j->begin(); k != j->end(); ++k) {
	      index = *k;
	      vector3 tmpvec = tmpmol.GetAtom(index)->GetVector();
              tmpvec *= xzmat; //apply the rotation
	      tmpmol.GetAtom(index)->SetVector(tmpvec);
	    }
	    
	    fragvec = GetNewBondVector(tmpmol, tmpmol.GetAtom(a->GetIdx()));
	    fragdir = fragvec - tmpmol.GetAtom(a->GetIdx())->GetVector();
	    cerr << "  fragdir = " << fragdir << endl;

	    cerr << "  --- YZ plane ---" << endl;
	    yzang = vectorAngle(vector3(moldir.y(), moldir.z(), 0.0), vector3(fragdir.y(), fragdir.z(), 0.0));
	    cerr << "cross: " << cross(vector3(moldir.y(), moldir.z(), 0.0), vector3(fragdir.y(), fragdir.z(), 0.0)) << endl;
	    if (cross(vector3(moldir.y(), moldir.z(), 0.0), vector3(fragdir.y(), fragdir.z(), 0.0)).z() > 0) {
	      cerr << "counterclockwise: " << yzang << endl;
	      yzang = 180 + yzang;
	    } else if (cross(vector3(moldir.y(), moldir.z(), 0.0), vector3(fragdir.y(), fragdir.z(), 0.0)).z() < 0) {
	      cerr << "clockwise: " << yzang << endl;
	      yzang = 180 - yzang;
	    } else
	      yzang = 0.0;
	    yzmat.SetupRotMat(yzang, 0.0, 0.0); 
	    cerr << " yzang = "  << yzang << endl;
	
	    for (k = j->begin(); k != j->end(); ++k) {
	      index = *k;
	      vector3 tmpvec = tmpmol.GetAtom(index)->GetVector();
              tmpvec *= yzmat; //apply the rotation
	      tmpmol.GetAtom(index)->SetVector(tmpvec);
	    }


	    fragvec = GetNewBondVector(tmpmol, tmpmol.GetAtom(a->GetIdx()));
	    fragdir = fragvec - tmpmol.GetAtom(a->GetIdx())->GetVector();
	    xyang = vectorAngle(vector3(moldir.x(), moldir.y(), 0.0), vector3(fragdir.x(), fragdir.y(), 0.0));
	    xzang = vectorAngle(vector3(moldir.x(), moldir.z(), 0.0), vector3(fragdir.x(), fragdir.z(), 0.0));
	    yzang = vectorAngle(vector3(moldir.y(), moldir.z(), 0.0), vector3(fragdir.y(), fragdir.z(), 0.0));
	    cerr << "xy angle after rotation: " << xyang << endl;
	    cerr << "xz angle after rotation: " << xzang << endl;
	    cerr << "yz angle after rotation: " << yzang << endl;
    
    
	    // 
	    // translate fragment 
	    //
	    for (k = j->begin(); k != j->end(); ++k) {
	      index = *k;
	      
	      if (index == a->GetIdx()) {
	        rootvec = tmpmol.GetAtom(a->GetIdx())->GetVector();
	    
	        //cout << "a position: " << rootvec << endl;
		  
	        for (k2 = j->begin(); k2 != j->end(); ++k2) {
	          index = *k2;

                  // translate the fragment
	          vector3 tmpvec = tmpmol.GetAtom(index)->GetVector();
	          tmpvec += molvec - rootvec;
                  tmpmol.GetAtom(index)->SetVector(tmpvec);
		  //cout << "tmpvec = " << tmpvec << endl;
     	        }
	      }
	    }
	    
	    // add bond between previous part and added fragment
            OBBond *bond = a->GetBond(prev);
            tmpmol.AddBond(*bond);
            } //if (prev != NULL) {



	  }
        }
      }


      if (foundfrag)
        continue;

      //
      // below is the code to add non-fragment atoms
      //

      visit[a->GetIdx()] = 1;
      
      // add the atom to tmpmol
      OBAtom *atom = tmpmol.GetAtom(a->GetIdx());
      atom->SetVector(molvec);

      // add bond between previous part and added atom
      if (prev != NULL) {
        OBBond *bond = a->GetBond(prev);
        tmpmol.AddBond(*bond);
      }

    }

    mol = tmpmol;

    return true;
  }


} // end namespace OpenBabel


//! \file builder.cpp
//! \brief Handle OBBuilder class

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
    while (ifs.getline(buffer, 80)) {
      tokenize(vs, buffer);
      
      
      if (vs[0] == "INDEX") {
        if (sp != NULL)
	  _fragments.push_back(pair<OBSmartsPattern*, vector<vector3> > (sp, coords));
	
        coords.clear();
	sp = new OBSmartsPattern;
        if (!sp->Init(vs[1])) {
          delete sp;
          sp = NULL;
          obErrorLog.ThrowError(__FUNCTION__, " Could not parse SMARTS from contribution data file", obInfo);
        }
      } else {
        vector3 coord(atof(vs[0].c_str()), atof(vs[1].c_str()), atof(vs[2].c_str()));
	coords.push_back(coord);
      }
 
        
    }
    
  }
    
  
  OBBuilder::~OBBuilder() 
  {
    _fragments.clear();
  }

  bool OBBuilder::AddFragment(OBMol &mol, int ai, OBMol &fragment, int si)
  {
  /* 
    vector<int> transl_idx(mol.NumAtoms()+1, 0);
    
    FOR_ATOMS_OF_MOL (b, fragment) {
      OBAtom *atom = mol.NewAtom();
      atom->SetVector(b->GetVector());

    }
                  
    // iterate over the bonds of fragment and add these bonds to 
    // tmpmol with corrected atoms
    FOR_BONDS_OF_MOL (b, fragment) {
                    b->SetBegin(tmpmol.GetAtom(transl_idx2[b->GetBeginAtomIdx()]));
                    b->SetEnd(tmpmol.GetAtom(transl_idx2[b->GetEndAtomIdx()]));
                    
		    cout << "ADD BOND: " << b->GetBeginAtomIdx() << "-" << b->GetEndAtomIdx() << endl;

                    tmpmol.AddBond(*b);
		  }
	
	          int index3 = 1;
	          for (k2 = j->begin(); k2 != j->end(); ++k2) {
		    index = *k2;
                    transl_idx[index] = transl_idx2[index3];
                    visit[index] = 1;

                    // set the atom type for tmpmol (types from mol). Also set AtomicNum and Hyb.
		    tmpmol.GetAtom(transl_idx[index])->SetType(mol.GetAtom(index)->GetType());
		    tmpmol.GetAtom(transl_idx[index])->SetAtomicNum(mol.GetAtom(index)->GetAtomicNum());
		    tmpmol.GetAtom(transl_idx[index])->SetHyb(mol.GetAtom(index)->GetHyb());
	            index3++;
		  }


	          if (prev != NULL) {
  		  // 
		  // rotate
		  //  
		  matrix3x3 mat;
		  double xyang, yzang, xzang, xyang1, yzang1, xzang1, xyang2, yzang2, xzang2;


		  cout << "----------------" << endl;
		  cout << " R O T A T E" << endl;
		  cout << "----------------" << endl;
		  cout << "  moldir = " << moldir << endl;
		  
		  cout << "  --- XY plane ---" << endl;
                  fragvec = GetNewBondVector(tmpmol, tmpmol.GetAtom(transl_idx[a->GetIdx()]));
		  fragdir = fragvec - tmpmol.GetAtom(transl_idx[a->GetIdx()])->GetVector();
	          xyang = GetRotationAngle(moldir.x(), moldir.y(), fragdir.x(), fragdir.y());
		  mat.SetupRotMat(0.0, 0.0, xyang); 
		  cout << "  fragdir = " << fragdir << endl;
		  cout << " xyang = "  << xyang << endl;

	          for (k2 = j->begin(); k2 != j->end(); ++k2) {
		    index = *k2;
		    vector3 tmpvec = tmpmol.GetAtom(transl_idx[index])->GetVector();
                    tmpvec *= mat; //apply the rotation
		    tmpmol.GetAtom(transl_idx[index])->SetVector(tmpvec);
		  }
	        
		  cout << "  --- XZ plane ---" << endl;
                  fragvec = GetNewBondVector(tmpmol, tmpmol.GetAtom(transl_idx[a->GetIdx()]));
		  fragdir = fragvec - tmpmol.GetAtom(transl_idx[a->GetIdx()])->GetVector();
	          xzang = GetRotationAngle(moldir.x(), moldir.z(), fragdir.x(), fragdir.z());
		  mat.SetupRotMat(0.0, xzang, 0.0); 
		  cout << "  fragdir = " << fragdir << endl;
		  cout << " xzang = "  << xzang << endl;
	          
		  for (k2 = j->begin(); k2 != j->end(); ++k2) {
		    index = *k2;
		    vector3 tmpvec = tmpmol.GetAtom(transl_idx[index])->GetVector();
                    tmpvec *= mat; //apply the rotation
		    tmpmol.GetAtom(transl_idx[index])->SetVector(tmpvec);
		  }

		  cout << "  --- YZ plane ---" << endl;
                  fragvec = GetNewBondVector(tmpmol, tmpmol.GetAtom(transl_idx[a->GetIdx()]));
		  fragdir = fragvec - tmpmol.GetAtom(transl_idx[a->GetIdx()])->GetVector();
		  yzang = GetRotationAngle(moldir.y(), moldir.z(), fragdir.y(), fragdir.z());
		  mat.SetupRotMat(yzang, 0.0, 0.0); 
		  cout << "  fragdir = " << fragdir << endl;
		  cout << " yzang = "  << yzang << endl;
	
	          for (k2 = j->begin(); k2 != j->end(); ++k2) {
		    index = *k2;
		    vector3 tmpvec = tmpmol.GetAtom(transl_idx[index])->GetVector();
                    tmpvec *= mat; //apply the rotation
		    tmpmol.GetAtom(transl_idx[index])->SetVector(tmpvec);
		  }
		  
		  // 
		  // translate
		  //
                  rootvec = tmpmol.GetAtom(transl_idx[a->GetIdx()])->GetVector();
		  
		  for (k2 = j->begin(); k2 != j->end(); ++k2) {
		    index = *k2;

		    // translate the fragment
		    vector3 tmpvec = tmpmol.GetAtom(transl_idx[index])->GetVector();
		    tmpvec += molvec - rootvec;
		    tmpmol.GetAtom(transl_idx[index])->SetVector(tmpvec);
		    cout << "tmpvec = " << tmpvec << endl;
 
		  }
		  } // if (prev != NULL)


		  // add bond between previous part and added fragment
                  if (prev != NULL) {
                    OBBond *bond;
                    bond = a->GetBond(prev);
                    bond->SetBegin(tmpmol.GetAtom(transl_idx[a->GetIdx()]));
                    bond->SetEnd(tmpmol.GetAtom(transl_idx[prev->GetIdx()]));
                    tmpmol.AddBond(*bond);
                  }

*/
    return false; 
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

  int OBBuilder::GetQuadrant(double x, double y) {
    if ((x > 0.0) || (y == 0.0))
      return 1;
    if ((x == 0.0) || (y > 0.0))
      return 2;
    if ((x < 0.0) || (y == 0.0))
      return 3;
    if ((x == 0.0) || (y < 0.0))
      return 4;
    
    if ((x > 0.0) && (y > 0.0))
      return 1;
    if ((x < 0.0) && (y > 0.0))
      return 2;
    if ((x < 0.0) && (y < 0.0))
      return 3;
    if ((x > 0.0) && (y < 0.0))
      return 4;
  }
  
  double OBBuilder::GetRotationAngle(double x, double y, double u, double v)
  {
    double ang, ang1, ang2, diff, q;
    
    ang1 = vectorAngle(vector3(x, y, 0.0), VX);
    ang2 = vectorAngle(vector3(u, v, 0.0), VX);
    //angy1 = vectorAngle(vector3(x, y, 0.0), VY);
    //angy2 = vectorAngle(vector3(u, v, 0.0), VY);
   
    //if (IsNan(ang1))
    //  ang1 = 0.0;

    cout << "           mol  ang1 = " << ang1 << "  quadrant = " << GetQuadrant(x, y) << endl;
    cout << "           frag ang2 = " << ang2 << "  quadrant = " << GetQuadrant(u, v) << endl;
    
    // make corrections to angles so that all angles are between
    // the vector and the x axis in anti-clockwise direction.
    if ((GetQuadrant(x, y) == 3) || GetQuadrant(x, y) == 4)
      ang1 = 180.0 + (180.0 - ang1);
    if ((GetQuadrant(u, v) == 3) || GetQuadrant(u, v) == 4)
      ang2 = 180.0 + (180.0 - ang2);

    if (ang2 > ang1) {
      ang = ang2 - ang1;
      if (ang > 180.0)
        ang -= 180.0;
      else
        ang += 180.0;
    } else {
      ang = ang1 - ang2;
      if (ang > 180.0)
        ang = 360.0 - (ang - 180.0);
      else
        ang = 180.0 - ang;
    }
 
    


    /*
    // set both angles to same quadrant
    q = 90.0 * (GetQuadrant(u, v) - GetQuadrant(x, y));
    
    if ((GetQuadrant(x, y) == GetQuadrant(u, v)) || ((GetQuadrant(u, v) - GetQuadrant(x, y)) == 1)) {
      if (ang1 < ang2)
        ang = 180.0 + ang2 + ang1 + q;
      else
        //ang = 90.0 + ang2 - ang1 + q;
        ang = 90.0 + ang2 + ang1 + q;
    }
    if ((GetQuadrant(u, v) - GetQuadrant(x, y)) == 2) {
      if (ang1 < ang2) 
        ang = -ang2 + ang1 + q;
      else
        ang = ang2 + ang1 + q; 
    }
    */

    return ang;
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
    
    // delete all bonds in tmpmol
    while (tmpmol.NumBonds())
      tmpmol.DeleteBond(tmpmol.GetBond(0));
    
    // iterate over all atoms to place them in 3D space
    FOR_DFS_OF_MOL (a, mol) {
      cout << "placing atom: " << a->GetIdx() << endl;
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
        cout << "prev: " << prev->GetIdx() << endl;

        molvec = GetNewBondVector(tmpmol, tmpmol.GetAtom(prev->GetIdx()));
	moldir = molvec - tmpmol.GetAtom(prev->GetIdx())->GetVector();
      } else {
        molvec = VX;
        moldir = VX;
      }
      
      if (prev != NULL)
        cout << "  prev: " << prev->GetIdx() << endl;
      else
        cout << "  prev: -" << endl;
      cout << "  molvec: " << molvec << endl;
      cout << "  moldir: " << moldir << endl;

      // check if a is in a fragment
      bool foundfrag = false;
      for (i = _fragments.begin();i != _fragments.end();++i) {
	if (i->first->Match(mol)) { 
          _mlist = i->first->GetMapList();
          
	  for (j = _mlist.begin();j != _mlist.end();++j) {
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

	      // add the bonds for the fragment
	      for (k2 = j->begin(); k2 != j->end(); ++k2) {
	        index = *k2;
	        OBAtom *atom1 = mol.GetAtom(index);
	        
		for (k3 = j->begin(); k3 != j->end(); ++k3) {
	          index2 = *k3;
	          OBAtom *atom2 = mol.GetAtom(index2);
                  OBBond *bond = atom1->GetBond(atom2);

                  if (bond != NULL) 
                    tmpmol.AddBond(*bond);
		}
              }

                  // ADDFRAGMENT HERE

	
              //  a->GetNewBondVector(fragvec, 1.5);
	    }
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

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
#include <unistd.h>

using namespace std;
using namespace OpenBabel;

// PROTOTYPES /////////////////////////////////////////////////////////////////
//int get_dst_nbr (OBAtom* atom);
//int get_ang_nbr (OBAtom* atom);
int get_nbr (OBAtom* atom, OBMol &mol, int level);
int get_ring_sign (OBMol &mol, OBAtom* atom_a);
bool is_in_same_ring (OBMol &mol, OBAtom* atom_a, OBAtom* atom_b);


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
  OBFormat *format_out = conv.FindFormat("MOL");
    
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

  for (c=1;;++c)
    {
      mol.Clear();
      if (!conv.Read(&mol, &ifs))
        break;
      if (mol.Empty())
        break;
        
      mol.AddHydrogens(false, true);

      OBAtom *atom, *nbr, *nbr2, *nbr3, *nbrs;
      vector<OBAtom*>::iterator i;
      vector<OBBond*>::iterator j;

      vector<OBInternalCoord*> internals;
      OBInternalCoord *coord;

      coord = new OBInternalCoord();
      internals.push_back(coord);
      mol.BeginModify();
	
      for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
        {
          coord = new OBInternalCoord();
          nbr = mol.GetAtom(get_nbr(atom, mol, 1));
          nbr2 = mol.GetAtom(get_nbr(atom, mol, 2));
          nbr3 = mol.GetAtom(get_nbr(atom, mol, 3));

          cout << "ATOM : " << atom->GetIdx();
          cout << "-" << get_nbr(atom, mol, 1);
          cout << "-" << get_nbr(atom, mol, 2);
          cout << "-" << get_nbr(atom, mol, 3) << endl;

          //
          // CHOOSE DISTANCE
          //
          if (nbr) {
            coord->_a = mol.GetAtom(get_nbr(atom, mol, 1));
            OBBond *bond;
            if ( (bond = mol.GetBond(atom, nbr)) ) {
              coord->_dst = bond->GetEquibLength();
            }
          }
          // 
          // CHOOSE ANGLE
          //
          //if (get_nbr(atom, mol, 2)) {
          if (nbr2) {
            coord->_b = mol.GetAtom(get_nbr(atom, mol, 2));
            if (nbr->GetHyb() == 3)
              coord->_ang = 109;
            if (nbr->GetHyb() == 2)
	            coord->_ang = 120;
            if (nbr->GetHyb() == 1)
	            coord->_ang = 180;
 
            if (nbr->IsInRing()) {
              if (atom->IsInRingSize(6) && nbr->IsInRingSize(6) && nbr2->IsInRingSize(6)) {
                if (nbr->GetHyb() == 3)
                  coord->_ang = 110;	// internal angle of aliphatic 6-ring
                else
                  coord->_ang = 120;	// internal angle of aromatic 6-ring
              }
              if (atom->IsInRingSize(5) && nbr->IsInRingSize(5) && nbr2->IsInRingSize(5))
                coord->_ang = 108;	// internal angle of aromatic 
              if (atom->IsInRingSize(4) && nbr->IsInRingSize(4) && nbr2->IsInRingSize(4))
                coord->_ang = 90;
              if (atom->IsInRingSize(3) && nbr->IsInRingSize(3) && nbr2->IsInRingSize(3))
                coord->_ang = 60;
              if (!is_in_same_ring(mol, atom, nbr2)) {
                if (atom->IsInRingSize(5) || nbr2->IsInRingSize(5))
                  if (!atom->IsHydrogen()) //// ???
                    coord->_ang = 132;
              }
              if (!atom->IsInRing() && nbr->IsInRingSize(5) && nbr2->IsInRingSize(5)) {
                if (nbr->GetHyb() == 3)
                  coord->_ang = 109;
                else
                  coord->_ang = 126;
              }
              /////
              if (!atom->IsInRing() && nbr->IsInRingSize(3))
                coord->_ang = 117;
            } 
          }
          // 
          // CHOOSE TORSION
          //
          if (nbr3) {
            double bestangle, angle, bestscore, score;
            int nbr_count;

            coord->_c = mol.GetAtom(get_nbr(atom, mol, 3));
		    
            if (nbr->GetHyb() == 3) {
              internals.push_back(coord);
              InternalToCartesian(internals, mol);
 
              if (atom->IsInRingSize(6) && nbr->IsInRingSize(6) && nbr2->IsInRingSize(6)) {
                internals[atom->GetIdx()]->_tor = get_ring_sign(mol, atom)*58;
                continue;
	            }
              if (!atom->IsInRing() && nbr->IsInRingSize(6) && nbr2->IsInRingSize(6) && nbr3->IsInRingSize(6)) {
                internals[atom->GetIdx()]->_tor = 177;
                FOR_NBORS_OF_ATOM(tmp, nbr) {
                  if (tmp->GetIdx() >= atom->GetIdx())
                    break;
                  cout << "    NBR=" << tmp->GetIdx() << endl;
                  if (internals[tmp->GetIdx()]->_tor || (tmp->GetIdx() < 4)) {
                    cout << "    YES" << endl;
                    nbr_count++;
                  }
                }
                if (nbr_count == 3)
  		            internals[atom->GetIdx()]->_tor = get_ring_sign(mol, nbr)*65;
                nbr_count = 0;
                continue;
              }

              if (atom->IsInRing() && nbr->IsInRing() && nbr2->IsInRing() && !nbr3->IsInRing()) {
                internals[atom->GetIdx()]->_tor = 300;
                continue;
              }
	                   
              if (atom->IsInRing()) {
                if (!nbr3->IsInRing() || !is_in_same_ring(mol, atom, nbr3)) {
                  internals[atom->GetIdx()]->_tor = 180;
                } else
                  internals[atom->GetIdx()]->_tor = 0;
                continue;
              }
		    
              if (!atom->IsInRing() && nbr->IsInRingSize(3)) {
                internals[atom->GetIdx()]->_tor = 107;
                continue;
              }
		    
              for ( angle=60; angle<360; angle+=120) {
                internals[atom->GetIdx()]->_tor = angle;
                //for (nbrs = nbr->BeginNbrAtom(j);nbrs;nbrs = nbr->NextNbrAtom(j)) {
                FOR_NBORS_OF_ATOM(tmp, nbr) {
                  if (tmp->GetIdx() == nbr3->GetIdx())
                    internals[atom->GetIdx()]->_tor -= 60;
                }
                if (nbr->IsInRing() && !atom->IsInRing())
                  internals[atom->GetIdx()]->_tor -= 60;


                InternalToCartesian(internals, mol);
                score = 0;
                score = atom->GetDistance(nbr3);
                for (nbrs = nbr->BeginNbrAtom(j);nbrs;nbrs = nbr->NextNbrAtom(j)) {
                  //FOR_NBORS_OF_ATOM(tmp, nbr) {
                  if (nbrs->GetX() || nbrs->GetY() || nbrs->GetZ()) {
                    score += atom->GetDistance(nbrs);
                  }
                  if (score > bestscore) {
		                bestscore = score;
                    bestangle = internals[atom->GetIdx()]->_tor;
                  }
                }
              }

              bestscore = 0;
              internals[atom->GetIdx()]->_tor = bestangle;

              continue;
            }
            if (nbr->GetHyb() == 2) {
              internals.push_back(coord);
              InternalToCartesian(internals, mol);

              if (atom->IsInRing()) {
                if (!nbr3->IsInRing() || !is_in_same_ring(mol, atom, nbr3)) {
                  internals[atom->GetIdx()]->_tor = 180;
                } else
                  internals[atom->GetIdx()]->_tor = 0;
                continue;
              }

              if (nbr->IsInRing() && !nbr3->IsInRing()) {
                internals[atom->GetIdx()]->_tor = 0;
                continue;
              }
              if (nbr->IsInRing() && nbr3->IsInRing() && !is_in_same_ring(mol, nbr, nbr3)) {
                internals[atom->GetIdx()]->_tor = 0;
                continue;
              }


              for ( angle=0; angle<360; angle+=180) {
                internals[atom->GetIdx()]->_tor = angle;
                InternalToCartesian(internals, mol);
                score = 0;
                score = atom->GetDistance(get_nbr(atom, mol, 3));
                for (nbrs = nbr->BeginNbrAtom(j);nbrs;nbrs = nbr->NextNbrAtom(j)) {
                  if (nbrs->GetX() || nbrs->GetY() || nbrs->GetZ()) {
                    score += atom->GetDistance(nbrs);
                  }
                  if (score > bestscore) {
		                bestscore = score;
                    bestangle = internals[atom->GetIdx()]->_tor;
                  }
                }
              }

              bestscore = 0;
              internals[atom->GetIdx()]->_tor = bestangle;

              continue;
	
            }
          }
          //cout << atom->GetIdx() << "-" << get_nbr(atom, 1) << "-" << get_nbr(atom, 2) << "-" << get_nbr(atom, 3) << endl;    
            
          internals.push_back(coord);

        }

      // DEBUG
      int k = 0;
      for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i)) {
        k++;
        cout << "ATOM " << k << " - dst: " << internals[k]->_dst << "   ang: " << internals[k]->_ang << "   tor: " << internals[k]->_tor << endl;
      }

      InternalToCartesian(internals, mol);
      mol.AddHydrogens(false, true);

      char FileOut[BUFF_SIZE];
      sprintf(FileOut, "%s-%d.mol", basename.c_str(), c);
      ofs.open(FileOut);
      conv.Write(&mol, &ofs);
      ofs.close();
      //conv.Write(&mol, &cout);
    } // end for loop

  return(1);
}

///////////////////////////////////////////////////////////////////////////////

int get_nbr (OBAtom* atom, OBMol &mol, int level) {
  OBAtom *nbr,*nbr2,*nbr3;
  vector<OBBond*>::iterator i;
  
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

///////////////////////////////////////////////////////////////////////////////

bool is_in_same_ring(OBMol &mol, OBAtom* atom_a, OBAtom* atom_b)
{
  bool a_in, b_in;
  vector<OBRing*> vr;
  vr = mol.GetSSSR();
  
  vector<OBRing*>::iterator i;
  vector<int>::iterator j;
  for (i = vr.begin();i != vr.end();i++) {
    for(j = (*i)->_path.begin();j != (*i)->_path.end();j++) {
	    if (*j == atom_a->GetIdx())
	        a_in = true;
	    if (*j == atom_b->GetIdx())
        b_in = true;
    }
    if (a_in && b_in)
      return true;
    a_in = false;
    b_in = false;
  }
  
  return false;
}

int get_ring_sign (OBMol &mol, OBAtom* atom_a)
{
  int sign;
  vector<OBRing*> vr;
  vr = mol.GetSSSR();
  
  vector<OBRing*>::iterator i;
  vector<int>::iterator j;
  for (i = vr.begin();i != vr.end();i++) {
    sign = 1;
    for(j = (*i)->_path.begin();j != (*i)->_path.end();j++) {
	    if (*j == atom_a->GetIdx()) {
        return sign;
	    }
	    sign = -sign;
    }
  }
  
  return 0;
}

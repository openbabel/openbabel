/**********************************************************************
Copyright (C) 2000 by OpenEye Scientific Software, Inc.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include "mol.h"

using namespace std;

namespace OpenBabel {

  void WriteCharges(ostream &ofs,OBMol &mol)
  {
    unsigned int i;
    OBAtom *atom;
    char buffer[BUFF_SIZE];
    
    for(i = 1;i <= mol.NumAtoms(); i++)
      {
	atom = mol.GetAtom(i);
	sprintf(buffer,"%4s%4d   % 2.10f",
		etab.GetSymbol(atom->GetAtomicNum()),
		i,
		atom->GetPartialCharge());
	
	ofs << buffer << endl;
      }
  }

  void WriteDistanceMatrix(ostream &ofs,OBMol &mol)
  {
    int columns = 7;
    unsigned int max, min = 1;
    unsigned int i,j;
    string type;
    OBAtom *atom, *atom2;
    char buffer[BUFF_SIZE];
    double dst;

    max = columns;
    while (max <= mol.NumAtoms() + columns)
      {
	ofs << endl;
	if (min > mol.NumAtoms()) break;
	atom = mol.GetAtom(min);
	
	sprintf(buffer,"%15s%4d",
		etab.GetSymbol(atom->GetAtomicNum()),
		min);
	ofs << buffer;
	
	for (i = min + 1; ((i < max) && (i <= mol.NumAtoms())); i++)
	  if (i <= mol.NumAtoms())
	    {
	      atom = mol.GetAtom(i);
	      sprintf(buffer,"%7s%4d",
		etab.GetSymbol(atom->GetAtomicNum()),
		i);
	      ofs << buffer;
	    }
	ofs << endl;

	sprintf(buffer,"%14s","");
	ofs << buffer;
	for (i = min; i < max; i++)
	  if (i <= mol.NumAtoms())
	    {
	      sprintf(buffer,"%11s","-----------");
	      ofs << buffer;
	    }
      
	ofs << endl;
	for (i = min; i <= mol.NumAtoms(); i++)
	  {
	    atom = mol.GetAtom(i);
	    sprintf(buffer,"%4s%4d",
		    etab.GetSymbol(atom->GetAtomicNum()),
		    i);
	    ofs << buffer;
	    for (j = min; j < max; j++)
	      if (j <= i)
		{
		  atom2 = mol.GetAtom(j);
		  dst = SQUARE(atom->GetX() - atom2->GetX());
		  dst += SQUARE(atom->GetY() - atom2->GetY());
		  dst += SQUARE(atom->GetZ() - atom2->GetZ());
		  dst = sqrt(dst);
		  sprintf(buffer,"%10.4f ",dst);
		  ofs << buffer;
		}
	    ofs << endl;
	  }
	max += columns - 1;
	min += columns - 1;
      }
    ofs << endl;
}


//  void print_torsions(ums_type *mol,FILE *file1)
//  {
//    int a,b,c,d;
//    int i,j,k;
//    double angle1;
//    int angle_count = 0;
//    torsion_rec *tr;


//    tr = (torsion_rec *)malloc(Atoms * 10 * sizeof(torsion_rec));
//    if (tr == NULL)
//    {
//      printf("Memory Allocation Error\n");
//      exit(0);
//    }
  //    return(CalcTorsionAngle(_atom[a-1]->GetVector(),
  //			       _atom[b-1]->GetVector(),
  //			       _atom[c-1]->GetVector(),
  //			       _atom[d-1]->GetVector()));
//    for (i = 0; i < Bonds; i++)
//    {
//      b = Start(i);
//      c = End(i);
//      for (j = 0; j < Valence(Start(i)); j ++)
//        if (Connection(Start(i),j) != End(i))
//        {
//  	a = Connection(Start(i),j);
//  	for (k = 0; k < Valence(End(i)); k ++)
//  	  if ((Connection(End(i),k) != Start(i)) &&
//  	      (Connection(End(i),k) != a))
//  	  {
//  	    d = Connection(End(i),k);
//  	    tr[angle_count].a = a;
//  	    tr[angle_count].b = b;
//  	    tr[angle_count].c = c;
//  	    tr[angle_count].d = d;
//  	    angle_count ++;
//  	  }
//        }
//    }
//    qsort(tr,angle_count,sizeof(torsion_rec),QSORT_PROTO compare_torsion);
//    for (i = 0; i < angle_count; i++)
//    {
//      angle1 = torsion(Point(tr[i].a),Point(tr[i].b),Point(tr[i].c),Point(tr[i].d));
//      fprintf(file1,"%4d %4d %4d %4d %10.3f\n",tr[i].a,tr[i].b,tr[i].c,tr[i].d,angle1);
//    }
  
//    free(tr); 
//  }	
  

  void WriteAngles(ostream &ofs,OBMol &mol)
  {
    // Alas, we still need to sort these to only list unique entries...
    vector3 v1, v2;
    OBAtom *a, *b, *c, *d;
    OBBond *bond1, *bond2, *bond3;
    vector<OBEdgeBase*>::iterator i, j, k;
    char buffer[BUFF_SIZE];

    for (bond1 = mol.BeginBond(i); bond1; bond1 = mol.NextBond(i))
      {
	b = bond1->GetBeginAtom();
	c = bond1->GetEndAtom();
	ofs << " outer " << endl;

	for (bond2 = b->BeginBond(j); bond2; bond2 = b->NextBond(j))
	  {
	    if (bond2->GetEndAtomIdx() != c->GetIdx() 
		&& bond2->GetEndAtomIdx() != b->GetIdx())
	      {
		a = bond2->GetEndAtom();

		v1 = a->GetVector() - b->GetVector();
		v2 = c->GetVector() - b->GetVector();

		sprintf(buffer,"%4d %4d %4d %4s %4s %4s %10.3f",
			a->GetIdx(),b->GetIdx(),c->GetIdx(),
			a->GetType(),b->GetType(),c->GetType(),
			vectorAngle(v1, v2));
		ofs << buffer << endl;

		for (bond3 = c->BeginBond(k); bond3; bond3 = c->NextBond(k))
		  if (bond3->GetEndAtomIdx() != b->GetIdx()
		      && bond3->GetEndAtomIdx() != c->GetIdx())
		    {
		      d = bond3->GetEndAtom();

		      v1 = b->GetVector() - c->GetVector();
		      v2 = d->GetVector() - c->GetVector();

		      sprintf(buffer,"%4d %4d %4d %4s %4s %4s %10.3f",
			      b->GetIdx(),c->GetIdx(),d->GetIdx(),
			      b->GetType(),c->GetType(),d->GetType(),
			      vectorAngle(v1, v2));
		      ofs << buffer << endl;
		    }
	      }
	  }
      }
  }

  void WriteChiral(ostream &ofs,OBMol &mol)
  {
    OBAtom *atom;
    vector<OBNodeBase*>::iterator i;
    char buffer[BUFF_SIZE];

    for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
      {
	if (atom->IsChiral())
	  {
	    sprintf(buffer,"%4s %5d is chiral: %s",
		    etab.GetSymbol(atom->GetAtomicNum()),
		    atom->GetIdx(),
		    (atom->IsClockwise() ? "clockwise" : "counterclockwise"));
	
	    ofs << buffer << endl;
	  }
      }
  }

bool WriteReport(ostream &ofs,OBMol &mol)
{
  ofs << "FILENAME: " << mol.GetTitle() << endl;
  ofs << "INTERATOMIC DISTANCES" << endl;
  WriteDistanceMatrix(ofs, mol);
  ofs << endl << endl << "ATOMIC CHARGES" << endl;
  WriteCharges(ofs, mol);
  ofs << endl << endl << "BOND ANGLES" << endl;
  WriteAngles(ofs, mol);
  //  ofs << endl << endl << "TORSION ANGLES" << endl;
  //  WriteTorsions(ofs, mol);
  ofs << endl << endl << "CHIRAL ATOMS" << endl;
  WriteChiral(ofs, mol);

  return(true);
}

}

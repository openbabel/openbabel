/**********************************************************************
Copyright (C) 1998-2003 by OpenEye Scientific Software, Inc.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include "mol.h"
#include "math/matrix3x3.h"

using namespace std;

namespace OpenBabel
{

  // ShelX homepage http://shelx.uni-ac.gwdg.de/SHELX/
  //   and 	    http://www.msg.ucsf.edu/local/programs/shelxl/SHELX_97.html

//   while (fgets(the_line,sizeof(the_line), file1) != NULL)
//   {
//     if (count_tokens(the_line,"\n\t ") > 0)
//       if EQ(gettoken(the_line,"\n\t ",1),"CELL")
//       {
// 	found = TRUE;
// 	sscanf(the_line,"%*s%*s%lf%lf%lf%lf%lf%lf",
// 	       &f.A,&f.B,&f.C,&f.Alpha,&f.Beta,&f.Gamma);
// 	fill_orth_matrix(&f,&m);
//       }
//       else
// 	if (found)
// 	{
// 	  if (is_good_shelx_line(the_line))
// 	    Atoms++;
// 	}      
//   }
//   ShowProgress(Atoms,"Reading Atoms");
//   result = initialize_ums(&mol);
//   rewind(file1);
//   i = 0;
//   while (fgets(the_line,sizeof(the_line), file1) != NULL)
//     if (is_good_shelx_line(the_line))
//     {
//       UpdateProgress();
//       i++;
//       sscanf(the_line,"%s %*s %lf %lf %lf",
// 	     Type(i),
// 	     &X(i),
// 	     &Y(i),
// 	     &Z(i));
//       check_shelx_coords(&Point(i)); 
//       clean_atom_type(Type(i));
//       fract_to_cart(&Point(i),&m); 
//     }
//   result = assign_radii(mol);
//   result = assign_bonds(mol);
//   result = assign_types(mol);
//   result = build_connection_table(mol);
//   assign_bond_order(mol);
//   return(TRUE);
// }

// int count_shelx_atoms(char *the_line)
// {
//   int atom_count = 0;
//   int i;
//   int tokens;
//   char the_token[20];  
//   tokens = count_tokens(the_line,"\t\n ");
//   for (i = 1; i <= tokens; i++)
//   {
//     strcpy(the_token,gettoken(the_line,"\t\n ",i));
//     if (isdigit(the_token[0]))
//       atom_count += atoi(the_token);
//   }
//   return(atom_count);
// }

	   
//   if (p->x > 10.0)
//     p->x -= 10.0;
//   if (p->y > 10.0)
//     p->y -= 10.0;
//   if (p->z > 10.0)
//     p->z -= 10.0;


// int is_good_shelx_line(char *the_line)
// {
//   char first[BUFF_SIZE];
//   int possible = FALSE;
//   int has_digit = FALSE;
//   int i;
//   if ((strchr(the_line,'(')) || (strchr(the_line,')')))
//     return(FALSE);
//   if ((count_tokens(the_line,"\n\t ") >= 4) && (isalpha(the_line[0])))
//   {
//     strcpy(first,gettoken(the_line," \n\t",1));
//     if (isdigit(first[0]))
//       return(FALSE);
//     for (i = 0; i < (int) strlen(first); i++)
//     {
//       if (isdigit(first[i]))
//       {
// 	has_digit = TRUE;
// 	break;
//       }
//     }
//     if ((!has_digit) && (strlen(first) > 2))
//       return(FALSE);
//     clean_atom_type(first);
//     if (is_element(first))
//     {
//       return(TRUE);
//     }
//   }
//   return(FALSE);
// }

bool ReadShelX(istream &ifs,OBMol &mol,const char *title)
{
  char buffer[BUFF_SIZE];
  int natoms;
  double A,B,C,Alpha,Beta,Gamma;
  matrix3x3 m;
  
  ifs.getline(buffer,BUFF_SIZE); mol.SetTitle(buffer);
  ifs.getline(buffer,BUFF_SIZE); sscanf(buffer,"%d",&natoms);
  
  while (ifs.getline(buffer,BUFF_SIZE) &&!EQn(buffer,"CELL",4));

  if (!EQn(buffer,"CELL",4)) return(false);
  vector<string> vs;
  tokenize(vs,buffer," \n\t,");
  if (vs.size() != 7) return(false);

  //parse cell values
  A = atof((char*)vs[1].c_str());
  B = atof((char*)vs[2].c_str());
  C = atof((char*)vs[3].c_str());
  Alpha = atof((char*)vs[4].c_str());
  Beta  = atof((char*)vs[5].c_str());
  Gamma = atof((char*)vs[6].c_str());

  OBUnitCell *uc = new OBUnitCell;
  uc->SetData(A, B, C, Alpha, Beta, Gamma);
  mol.SetData(uc);
  m = uc->GetOrthoMatrix();

  int i;
  double x,y,z;
  char type[10];
  OBAtom *atom;
  vector3 v;

  for (i = 1; i <= natoms;i++)
    {
      if (!ifs.getline(buffer,BUFF_SIZE)) return(false);
      tokenize(vs,buffer," \n\t,");
      if (vs.size() < 4) return(false);
      atom = mol.NewAtom();

      x = atof((char*)vs[1].c_str());
      y = atof((char*)vs[2].c_str());
      z = atof((char*)vs[3].c_str());
      v.Set(x,y,z); v *= m;

      strcpy(type,vs[0].c_str());
      atom->SetAtomicNum(etab.GetAtomicNum(type));
      atom->SetVector(v);
    }

  mol.ConnectTheDots();
  mol.PerceiveBondOrders();

  return(true);
}

} // end namespace

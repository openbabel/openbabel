/**********************************************************************
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
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

#ifndef OB_CWRAP_H
#define OB_CWRAP_H

int      ob_dbase_get_number_of_molecules(long int);
int      ob_dbase_get_cmol(int,long int,long int);
int      ob_get_cmol_atom_number(int*,long int);
int      ob_get_cmol_bond_number(int*,long int);
int      ob_get_cmol_coordinates(double**,long int);
int      ob_get_cmol_element(char*,long int);
int      ob_get_cmol_name(char*,long int);
int      ob_get_cmol_conformer_number(int*,long int);
int      ob_get_cmol_conformer(double**,int,long int);
void     ob_delete_cmol(long int);
long int ob_make_dbase(char*);
long int ob_make_cmol();

#endif // OB_CWRAP_H

/**********************************************************************
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#ifndef CWRAP_H
#define CWRAP_H

int      oe_dbase_get_number_of_molecules(long int);
int      oe_dbase_get_cmol(int,long int,long int);
int      oe_get_cmol_atom_number(int*,long int);
int      oe_get_cmol_bond_number(int*,long int);
int      oe_get_cmol_coordinates(float**,long int);
int      oe_get_cmol_element(char*,long int);
int      oe_get_cmol_name(char*,long int);
int      oe_get_cmol_conformer_number(int*,long int);
int      oe_get_cmol_conformer(float**,int,long int);
void     oe_delete_cmol(long int);
long int oe_make_dbase(char*);
long int oe_make_cmol();

#endif //CWRAP_H

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

#include "bitgrid.h"

#ifndef ATOMTYPE_PATH
#define ATOMTYPE_PATH "hbtypes.txt"
#endif

using namespace OpenBabel;
using namespace std;

BitGrid::BitGrid(void)
{
  p.read_rules(ATOMTYPE_PATH);
  fuzzy = false;
}

BitGrid::BitGrid(bool f)
{
  p.read_rules(ATOMTYPE_PATH);
  fuzzy = f;
}

BitGrid::~BitGrid(void)
{
  Clear();
}

void BitGrid::Init(OBMol &box, float spacing)
{
  OBAtom *atom;
  vector<OBNodeBase*>::iterator i;

  for ( atom = box.BeginAtom(i) ; atom ; atom = box.NextAtom(i) )
    {
      if (atom->GetIdx() == 1)
	{
	  xmin = atom->GetX(); xmax = atom->GetX();
	  ymin = atom->GetY(); ymax = atom->GetY();
	  zmin = atom->GetZ(); zmax = atom->GetZ();
	}
      else
	{
	  if (atom->GetX() < xmin) xmin = atom->GetX();
	  if (atom->GetX() > xmax) xmax = atom->GetX();
	  if (atom->GetY() < ymin) ymin = atom->GetY();
	  if (atom->GetY() > ymax) ymax = atom->GetY();
	  if (atom->GetZ() < zmin) zmin = atom->GetZ();
	  if (atom->GetZ() > zmax) zmax = atom->GetZ();
	}
    }

  midx = 0.5 * (xmax+xmin);
  midy = 0.5 * (ymax+ymin);
  midz = 0.5 * (zmax+zmin);
 
  xdim  = 3 + (int) ((xmax-xmin)/spacing);
  ydim  = 3 + (int) ((ymax-ymin)/spacing);
  zdim  = 3 + (int) ((zmax-zmin)/spacing);
  xydim = xdim * ydim;
  size  = xdim * ydim * zdim;

  inv_spa = 1.0 / spacing;
  this->spacing = spacing;

  grid.Clear();
  lipo.Clear();
  don.Clear();
  acc.Clear();

  grid.Resize(size);
  lipo.Resize(size);
  don.Resize(size);
  acc.Resize(size);
}

void BitGrid::Init(float xmi, float ymi, float zmi, float xma, float yma, float zma, float spacing)
{
  xmin = xmi; xmax = xma;
  ymin = ymi; ymax = yma;
  zmin = zmi; zmax = zma;

  midx = 0.5 * (xmax+xmin);
  midy = 0.5 * (ymax+ymin);
  midz = 0.5 * (zmax+zmin);
 
  xdim = 3 + (int) ((xmax-xmin)/spacing);
  ydim = 3 + (int) ((ymax-ymin)/spacing);
  zdim = 3 + (int) ((zmax-zmin)/spacing);
  size = xdim * ydim * zdim;

  inv_spa = 1.0 / spacing;
  this->spacing = spacing;

  grid.Clear();
  lipo.Clear();
  don.Clear();
  acc.Clear();

  grid.Resize(size);
  lipo.Resize(size);
  don.Resize(size);
  acc.Resize(size);
}

void BitGrid::Build(OBMol &mol)
{
  OBAtom *atom;
  vector<OBNodeBase*>::iterator i;

  p.assign_types(mol,types);

  for ( atom = mol.BeginAtom(i) ; atom ; atom = mol.NextAtom(i) )
    SetBits(atom);
}

void BitGrid::Build(OBMol &mol, OBBitVec &bits)
{
  OBAtom *atom;
  vector<OBNodeBase*>::iterator i;

  p.assign_types(mol,types);

  for ( atom = mol.BeginAtom(i) ; atom ; atom = mol.NextAtom(i) )
    if (bits.BitIsOn(atom->GetIdx()))
      SetBits(atom);
}

void BitGrid::Build(OBMol &mol, vector<int> &atoms)
{
  vector<int>::iterator i;

  p.assign_types(mol,types);

  for ( i = atoms.begin() ; i != atoms.end() ; i++ )
    SetBits(mol.GetAtom(*i));
}

void BitGrid::SetBits(OBAtom *atom)
{	
  int i,j,k;
  int xpos   = (int)rint((atom->GetX() - xmin) * inv_spa);
  int ypos   = (int)rint((atom->GetY() - ymin) * inv_spa);
  int zpos   = (int)rint((atom->GetZ() - zmin) * inv_spa);
  int startX = (0 > xpos - 1) ? 0 : xpos - 1;
  int startY = (0 > ypos - 1) ? 0 : ypos - 1;
  int startZ = (0 > zpos - 1) ? 0 : zpos - 1;
  int endX   = (xdim < xpos + 1) ? xdim : xpos + 1;
  int endY   = (ydim < ypos + 1) ? ydim : ypos + 1;
  int endZ   = (zdim < zpos + 1) ? zdim : zpos + 1;

  float rad  = etab.CorrectedVdwRad(atom->GetAtomicNum(),atom->GetHyb());

  if (fuzzy)
    {
      Vector grid_pt, mol_pt = atom->GetVector();
      for ( i = startX ; i <= endX ; i++ )
	{
	  grid_pt.SetX(((float)i * (float)spacing) + xmin);
	  for ( j = startY ; j <= endY ; j++ )
	    {
	      grid_pt.SetY(((float)j * (float)spacing) + ymin);
	      for ( k = startZ ; k <= endZ ; k++ )
		{
		  grid_pt.SetZ(((float)k * (float)spacing) + zmin);
		  if ((grid_pt - mol_pt).length() > rad)
		    continue;
		  
		  int n = i + (ydim*j) + (xydim*k);	 
		  
		  grid.SetBitOn(n);
		  
		  if (types[atom->GetIdx()] == "DON")
		    don.SetBitOn(n);
		  else if ( types[atom->GetIdx()] == "ACC")
		    acc.SetBitOn(n);
		  else
		    lipo.SetBitOn(n);
		}
	    }
	}
    }
  else
    {
      int n = xpos + (ydim*ypos) + (xydim*zpos);
      
      grid.SetBitOn(n);
      
      if (types[atom->GetIdx()] == "DON")
	don.SetBitOn(n);
      else if ( types[atom->GetIdx()] == "ACC")
	acc.SetBitOn(n);
      else
	lipo.SetBitOn(n);
    }
}



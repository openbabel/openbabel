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

#include "mol.h"
#include "grid.h"

using namespace std;

namespace OpenEye {

void OEProxGrid::Setup(OEMol &mol,OEMol &box,float cutoff,float res)
{
  OEAtom *atom;
  vector<OEAtom*>::iterator i;

  for (atom = box.BeginAtom(i);atom;atom = box.NextAtom(i))
      if (atom->GetIdx() == 1)
	{
	  _xmin = atom->GetX(); _xmax = atom->GetX();
	  _ymin = atom->GetY(); _ymax = atom->GetY();
	  _zmin = atom->GetZ(); _zmax = atom->GetZ();
	}
      else
	{
	  if (atom->GetX() < _xmin) _xmin = atom->GetX();
	  if (atom->GetX() > _xmax) _xmax = atom->GetX();
	  if (atom->GetY() < _ymin) _ymin = atom->GetY();
	  if (atom->GetY() > _ymax) _ymax = atom->GetY();
	  if (atom->GetZ() < _zmin) _zmin = atom->GetZ();
	  if (atom->GetZ() > _zmax) _zmax = atom->GetZ();
	}

  _inc = res; // 1/2 angstrom resolution

  _nxinc = (int) floor((_xmax - _xmin)/0.5);
  _nyinc = (int) floor((_ymax - _ymin)/0.5);
  _nzinc = (int) floor((_zmax - _zmin)/0.5);
  _maxinc = _nxinc*_nyinc*_nzinc;

  int j,size = _nxinc*_nyinc*_nzinc;
  cell.resize(size);
  for (unsigned int num = 0; num < cell.size(); num++)
    cell[num].resize(0);
  
  cutoff *= cutoff; //don't do sqrts

  int k,l,m;
  float x,y,z,dx_2,dy_2;
  float *c = mol.GetCoordinates();
  size = mol.NumAtoms()*3;

  for (atom = mol.BeginAtom(i),j=0;atom;atom = mol.NextAtom(i),j+=3)
      if (PointIsInBox(c[j],c[j+1],c[j+2]))
	for (x = _xmin+(_inc/2.0),k=0;k < _nxinc;x+=_inc,k++)
	  if ((dx_2 = SQUARE(c[j]-x)) < cutoff)
	    for (y = _ymin+(_inc/2.0),l=0;l < _nyinc;y+=_inc,l++)
	      if ((dx_2+(dy_2 = SQUARE(c[j+1]-y))) < cutoff)
		for (z = _zmin+(_inc/2.0),m=0;m < _nzinc;z+=_inc,m++)
		  if ((dx_2+dy_2+SQUARE(c[j+2]-z)) < cutoff)
		    cell[(k*_nyinc*_nzinc)+(l*_nzinc)+m].push_back(atom->GetIdx());
  
  _inc = 1/_inc;
}

void OEProxGrid::Setup(OEMol &mol,OEMol &box,float cutoff,vector<bool> &use,
		       float res)
{
  OEAtom *atom;
  vector<OEAtom*>::iterator i;

  for (atom = box.BeginAtom(i);atom;atom = box.NextAtom(i))
      if (atom->GetIdx() == 1)
	{
	  _xmin = atom->GetX(); _xmax = atom->GetX();
	  _ymin = atom->GetY(); _ymax = atom->GetY();
	  _zmin = atom->GetZ(); _zmax = atom->GetZ();
	}
      else
	{
	  if (atom->GetX() < _xmin) _xmin = atom->GetX();
	  if (atom->GetX() > _xmax) _xmax = atom->GetX();
	  if (atom->GetY() < _ymin) _ymin = atom->GetY();
	  if (atom->GetY() > _ymax) _ymax = atom->GetY();
	  if (atom->GetZ() < _zmin) _zmin = atom->GetZ();
	  if (atom->GetZ() > _zmax) _zmax = atom->GetZ();
	}

  _inc = res; // 1/2 angstrom resolution

  _nxinc = (int) floor((_xmax - _xmin)/0.5);
  _nyinc = (int) floor((_ymax - _ymin)/0.5);
  _nzinc = (int) floor((_zmax - _zmin)/0.5);
  _maxinc = _nxinc*_nyinc*_nzinc;

  int j,size = _nxinc*_nyinc*_nzinc;
  cell.resize(size);
  cutoff *= cutoff; //don't do sqrts

  int k,l,m;
  float x,y,z,dx_2,dy_2;
  float *c = mol.GetCoordinates();
  size = mol.NumAtoms()*3;

  for (atom = mol.BeginAtom(i),j=0;atom;atom = mol.NextAtom(i),j+=3)
    if (use[atom->GetIdx()])
      if (PointIsInBox(c[j],c[j+1],c[j+2]))
	for (x = _xmin+(_inc/2.0),k=0;k < _nxinc;x+=_inc,k++)
	  if ((dx_2 = SQUARE(c[j]-x)) < cutoff)
	    for (y = _ymin+(_inc/2.0),l=0;l < _nyinc;y+=_inc,l++)
	      if ((dx_2+(dy_2 = SQUARE(c[j+1]-y))) < cutoff)
		for (z = _zmin+(_inc/2.0),m=0;m < _nzinc;z+=_inc,m++)
		  if ((dx_2+dy_2+SQUARE(c[j+2]-z)) < cutoff)
		    cell[(k*_nyinc*_nzinc)+(l*_nzinc)+m].push_back(atom->GetIdx());
  
  _inc = 1/_inc;
}

vector<int> *OEProxGrid::GetProxVector(float x,float y,float z)
{
  if (x < _xmin || x > _xmax) return(NULL);
  if (y < _ymin || y > _ymax) return(NULL);
  if (z < _zmin || z > _zmax) return(NULL);

  x -= _xmin;y -= _ymin;z -= _zmin;
  int i,j,k,idx;
  i = (int) (x*_inc);j = (int) (y*_inc);k = (int) (z*_inc);
  idx = (i*_nyinc*_nzinc)+(j*_nzinc)+k;
  if (idx >= _maxinc) return(NULL);

  return(&cell[idx]);
}

vector<int> *OEProxGrid::GetProxVector(float *c)
{
  float x,y,z;
  x = c[0];
  y = c[1];
  z = c[2];
  
  if (x < _xmin || x > _xmax) return(NULL);
  if (y < _ymin || y > _ymax) return(NULL);
  if (z < _zmin || z > _zmax) return(NULL);

  x -= _xmin;y -= _ymin;z -= _zmin;
  int i,j,k,idx;
  i = (int) (x*_inc);j = (int) (y*_inc);k = (int) (z*_inc);
  idx = (i*_nyinc*_nzinc)+(j*_nzinc)+k;
  if (idx >= _maxinc) return(NULL);

  return(&cell[idx]);
}

void OEFloatGrid::Init(OEMol &box,float spacing, float pad)
{
  OEAtom *atom;
  vector<OEAtom*>::iterator i;

  for (atom = box.BeginAtom(i);atom;atom = box.NextAtom(i))
      if (atom->GetIdx() == 1)
	{
	  _xmin = atom->GetX(); _xmax = atom->GetX();
	  _ymin = atom->GetY(); _ymax = atom->GetY();
	  _zmin = atom->GetZ(); _zmax = atom->GetZ();
	}
      else
	{
	  if (atom->GetX() < _xmin) _xmin = atom->GetX();
	  if (atom->GetX() > _xmax) _xmax = atom->GetX();
	  if (atom->GetY() < _ymin) _ymin = atom->GetY();
	  if (atom->GetY() > _ymax) _ymax = atom->GetY();
	  if (atom->GetZ() < _zmin) _zmin = atom->GetZ();
	  if (atom->GetZ() > _zmax) _zmax = atom->GetZ();
	}

  _xmin -= pad;   _xmax += pad;
  _ymin -= pad;   _ymax += pad;
  _zmin -= pad;   _zmax += pad;
  float midx,midy,midz;
  int xdim,ydim,zdim;

  midx=0.5*(_xmax+_xmin);
  midy=0.5*(_ymax+_ymin);
  midz=0.5*(_zmax+_zmin);

  xdim=3+(int) ((_xmax-_xmin)/spacing);
  ydim=3+(int) ((_ymax-_ymin)/spacing);
  zdim=3+(int) ((_zmax-_zmin)/spacing);

/* store in Grid */
  _midx=midx;
  _midy=midy;
  _midz=midz;
  _xdim=xdim;
  _ydim=ydim;
  _zdim=zdim;
  _spacing=spacing;
  _halfSpace= _spacing/2.0f;
  _inv_spa=1.0/spacing;
  _val=NULL;
  _ival=NULL;

  int size = _xdim*_ydim*_zdim;
  _val = new float [size];
  memset(_val,'\0',sizeof(float)*size);

  //return(true);
}

void OEFloatGrid::IndexToCoords(int idx, float &x, float &y, float &z)
{
  long int grid_x,grid_y,grid_z;
  
  grid_x = idx % (int)_xdim;
  grid_z = (int)(idx /(_xdim * _ydim));
  grid_y = (int)((idx - (grid_z * _xdim * _ydim))/_xdim);

  x = (grid_x * _spacing + _xmin) + this->_halfSpace;
  y = (grid_y * _spacing + _ymin) + this->_halfSpace;
  z = (grid_z * _spacing + _zmin) + this->_halfSpace;
}

int OEFloatGrid::CoordsToIndex(float &x, float &y, float &z)
{
  int gx=(int)((x-_xmin)*_inv_spa);
  int gy=(int)((y-_ymin)*_inv_spa);
  int gz=(int)((z-_zmin)*_inv_spa);

  return((gz*_ydim*_xdim)+(gy*_xdim)+gx);
}

void OEFloatGrid::CoordsToIndex(int *idx,float *c)
{
  idx[0]=(int)((c[0]-_xmin)*_inv_spa);
  idx[1]=(int)((c[1]-_ymin)*_inv_spa);
  idx[2]=(int)((c[2]-_zmin)*_inv_spa);
}

float OEFloatGrid::Interpolate(float x,float y,float z)
{
  int n,igx,igy,igz;
  float xydim;
  float gx,gy,gz,fgx,fgy,fgz;
  float ax,ay,az,bx,by,bz;
  float AyA,ByA,AyB,ByB,Az,Bz;

  if((x<=_xmin)||(x>=_xmax)) return(0.0);
  if((y<=_ymin)||(y>=_ymax)) return(0.0);
  if((z<=_zmin)||(z>=_zmax)) return(0.0);

  xydim = _xdim*_ydim;

  /* calculate grid voxel and fractional offsets */
  gx=(x-_xmin-_halfSpace)*_inv_spa;   if (gx<0) gx=0; 
  igx=(int)gx;fgx=gx-(float)igx;
  gy=(y-_ymin-_halfSpace)*_inv_spa;   if (gy<0) gy=0;
  igy=(int) gy; fgy= gy - (float) igy;
  gz=(z-_zmin-_halfSpace)*_inv_spa;   if (gz<0) gz=0;
  igz=(int) gz; fgz= gz - (float) igz;
		
  n=(int)(igx+ _xdim*igy + xydim*igz);
  /* calculate linear weightings */
  ax=1.0-fgx; bx=fgx;     ay=1.0-fgy;
  by=fgy;     az=1.0-fgz; bz=fgz;
      
  /* calculate interpolated value */
  AyA=ax*_val[n           ]+bx*_val[(int)(n+1)    ];
  ByA=ax*_val[n+_xdim]+bx*_val[(int)(n+1+_xdim)  ];

  Az=ay*AyA+by*ByA;
		
  AyB=ax*_val[(int)(n     +xydim)]+bx*_val[(int)(n+1     +xydim)];
  ByB=ax*_val[(int)(n+_xdim+xydim)]+bx*_val[(int)(n+1+_xdim+xydim)];
  Bz=ay*AyB+by*ByB;

  return(az*Az+bz*Bz);
}


float OEFloatGrid::InterpolateDerivatives(float x,float y,float z,float *derivatives)
{
  int n,igx,igy,igz;
  float xydim;
  float gx,gy,gz,fgx,fgy,fgz;
  float ax,ay,az,bx,by,bz;
  float AyA,ByA,AyB,ByB,Az,Bz;
	float energy,fx,fy,fz;

  if((x<=_xmin)||(x>=_xmax)) return(0.0);
  if((y<=_ymin)||(y>=_ymax)) return(0.0);
  if((z<=_zmin)||(z>=_zmax)) return(0.0);

  xydim = _xdim*_ydim;

  /* calculate grid voxel and fractional offsets */
  gx=(x-_xmin-_halfSpace)*_inv_spa;   if (gx<0) gx=0; 
  igx=(int)gx;fgx=gx-(float)igx;
  gy=(y-_ymin-_halfSpace)*_inv_spa;   if (gy<0) gy=0;
  igy=(int) gy; fgy= gy - (float) igy;
  gz=(z-_zmin-_halfSpace)*_inv_spa;   if (gz<0) gz=0;
  igz=(int) gz; fgz= gz - (float) igz;
		
  n=(int)(igx+ _xdim*igy + xydim*igz);
  /* calculate linear weightings */
  ax=1.0-fgx; bx=fgx;     ay=1.0-fgy;
  by=fgy;     az=1.0-fgz; bz=fgz;
      
  /* calculate interpolated value */
  AyA=ax*_val[n           ]+bx*_val[(int)(n+1)    ];
  ByA=ax*_val[n+_xdim]+bx*_val[(int)(n+1+_xdim)  ];

  Az=ay*AyA+by*ByA;
		
  AyB=ax*_val[(int)(n      +xydim)]+bx*_val[(int)(n+1      +xydim)];
  ByB=ax*_val[(int)(n+_xdim+xydim)]+bx*_val[(int)(n+1+_xdim+xydim)];
  Bz=ay*AyB+by*ByB;

  energy = az*Az+bz*Bz;

	//calculate derivatives
	
	fz=-Az+Bz;

  Az=-AyA+ByA;
  Bz=-AyB+ByB;

  fy=az*Az+bz*Bz;

  AyA=-_val[n           ]+_val[(int)(n+1)     ];
  ByA=-_val[n+_xdim      ]+_val[(int)(n+1+_xdim)];

  Az=ay*AyA+by*ByA;

  AyB=-_val[(int)(n      +xydim)]+_val[(int)(n+1      +xydim)];
  ByB=-_val[(int)(n+_xdim+xydim)]+_val[(int)(n+1+_xdim+xydim)];

  Bz=ay*AyB+by*ByB;

  fx=az*Az+bz*Bz;

  derivatives[0] += fx;
  derivatives[1] += fy;
  derivatives[2] += fz;

	return(energy);

}


ostream& operator<< ( ostream &os, const  OEFloatGrid& fg) 
{
  os.write((const char*)&fg._xmin,sizeof(float));
  os.write((const char*)&fg._xmax,sizeof(float));
  os.write((const char*)&fg._ymin,sizeof(float));
  os.write((const char*)&fg._ymax,sizeof(float));
  os.write((const char*)&fg._zmin,sizeof(float));
  os.write((const char*)&fg._zmax,sizeof(float));

  os.write((const char*)&fg._midx,sizeof(float));
  os.write((const char*)&fg._midy,sizeof(float));
  os.write((const char*)&fg._midz,sizeof(float));
  os.write((const char*)&fg._inv_spa,sizeof(float));
  os.write((const char*)&fg._spacing,sizeof(float));
  os.write((const char*)&fg._xdim,sizeof(int));
  os.write((const char*)&fg._ydim,sizeof(int));
  os.write((const char*)&fg._zdim,sizeof(int));
  os.write((const char*)&fg._val[0],
	   (sizeof(float)*(fg._xdim*fg._ydim*fg._zdim)));

  return(os);
}

istream& operator>> ( istream &is,OEFloatGrid& fg) 
{
  is.read((char*)&fg._xmin,sizeof(float));
  is.read((char*)&fg._xmax,sizeof(float));
  is.read((char*)&fg._ymin,sizeof(float));
  is.read((char*)&fg._ymax,sizeof(float));
  is.read((char*)&fg._zmin,sizeof(float));
  is.read((char*)&fg._zmax,sizeof(float));

  is.read((char*)&fg._midx,sizeof(float));
  is.read((char*)&fg._midy,sizeof(float));
  is.read((char*)&fg._midz,sizeof(float));
  is.read((char*)&fg._inv_spa,sizeof(float));
  is.read((char*)&fg._spacing,sizeof(float));
  is.read((char*)&fg._xdim,sizeof(int));
  is.read((char*)&fg._ydim,sizeof(int));
  is.read((char*)&fg._zdim,sizeof(int));
  int size = fg._xdim*fg._ydim*fg._zdim;
  fg._val = new float [size];
  size *= sizeof(float);

  is.read((char*)&fg._val[0],size);
  fg._halfSpace= fg._spacing/2.0f;

  return(is);
}

#ifdef FOO
//linear interpolation routine
  float Interpolate(float x,float y,float z)
    {
      int n;
      int igx,igy,igz;
      float scale,gx,gy,gz,fgx,fgy,fgz;
      float qzpypx,qzpynx,qznynx,qznypx,qznyx,qzpyx;

      scale=_inv_spa;

      if((x<=_xmin)||(x>=_xmax)) return(0.0);
      if((y<=_ymin)||(y>=_ymax)) return(0.0);
      if((z<=_zmin)||(z>=_zmax)) return(0.0);

      gx=(x-_xmin-_halfSpace)*scale; igx=(int) gx; fgx= gx - (float) igx;
      gy=(y-_ymin-_halfSpace)*scale; igy=(int) gy; fgy= gy - (float) igy;
      gz=(z-_zmin-_halfSpace)*scale; igz=(int) gz; fgz= gz - (float) igz;

      int xydim=_xdim*_ydim;
      n=igx+ _xdim*igy + xydim*igz;
      qzpypx=fgx*(_val[n+1+_xdim+xydim]-_val[n+_xdim+xydim]) + _val[n+_xdim+xydim];
      qzpynx=fgx*(_val[n+1+xydim] -_val[n+xydim]) + _val[n+xydim];
      qznypx=fgx*(_val[n+1+_xdim] -_val[n+_xdim]) + _val[n+_xdim];
      qznynx=fgx*(_val[n+1] -_val[n]) + _val[n];
      qzpyx =fgy*(-qzpynx+qzpypx)+qzpynx;
      qznyx =fgy*(-qznynx+qznypx)+qznynx;

      return(fgz*(qzpyx-qznyx)+qznyx);
    }
#endif


}

/**********************************************************************
grid.cpp - Handle grids of values.
 
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2005 by Geoffrey R. Hutchison
 
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

#include "mol.h"
#include "grid.h"

using namespace std;

namespace OpenBabel
{

void OBProxGrid::Setup(OBMol &mol,OBMol &box,double cutoff,double res)
{
    OBAtom *atom;
    vector<OBNodeBase*>::iterator i;

    for (atom = box.BeginAtom(i);atom;atom = box.NextAtom(i))
        if (atom->GetIdx() == 1)
        {
            _xmin = atom->GetX();
            _xmax = atom->GetX();
            _ymin = atom->GetY();
            _ymax = atom->GetY();
            _zmin = atom->GetZ();
            _zmax = atom->GetZ();
        }
        else
        {
            if (atom->GetX() < _xmin)
                _xmin = atom->GetX();
            if (atom->GetX() > _xmax)
                _xmax = atom->GetX();
            if (atom->GetY() < _ymin)
                _ymin = atom->GetY();
            if (atom->GetY() > _ymax)
                _ymax = atom->GetY();
            if (atom->GetZ() < _zmin)
                _zmin = atom->GetZ();
            if (atom->GetZ() > _zmax)
                _zmax = atom->GetZ();
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
    double x,y,z,dx_2,dy_2;
    double *c = mol.GetCoordinates();
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

void OBProxGrid::Setup(OBMol &mol,OBMol &box,double cutoff,vector<bool> &use,
                       double res)
{
    OBAtom *atom;
    vector<OBNodeBase*>::iterator i;

    for (atom = box.BeginAtom(i);atom;atom = box.NextAtom(i))
        if (atom->GetIdx() == 1)
        {
            _xmin = atom->GetX();
            _xmax = atom->GetX();
            _ymin = atom->GetY();
            _ymax = atom->GetY();
            _zmin = atom->GetZ();
            _zmax = atom->GetZ();
        }
        else
        {
            if (atom->GetX() < _xmin)
                _xmin = atom->GetX();
            if (atom->GetX() > _xmax)
                _xmax = atom->GetX();
            if (atom->GetY() < _ymin)
                _ymin = atom->GetY();
            if (atom->GetY() > _ymax)
                _ymax = atom->GetY();
            if (atom->GetZ() < _zmin)
                _zmin = atom->GetZ();
            if (atom->GetZ() > _zmax)
                _zmax = atom->GetZ();
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
    double x,y,z,dx_2,dy_2;
    double *c = mol.GetCoordinates();
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

vector<int> *OBProxGrid::GetProxVector(double x,double y,double z)
{
    if (x < _xmin || x > _xmax)
        return(NULL);
    if (y < _ymin || y > _ymax)
        return(NULL);
    if (z < _zmin || z > _zmax)
        return(NULL);

    x -= _xmin;
    y -= _ymin;
    z -= _zmin;
    int i,j,k,idx;
    i = (int) (x*_inc);
    j = (int) (y*_inc);
    k = (int) (z*_inc);
    idx = (i*_nyinc*_nzinc)+(j*_nzinc)+k;
    if (idx >= _maxinc)
        return(NULL);

    return(&cell[idx]);
}

vector<int> *OBProxGrid::GetProxVector(double *c)
{
    double x,y,z;
    x = c[0];
    y = c[1];
    z = c[2];

    return( GetProxVector(x, y, z) );
}

void OBFloatGrid::Init(OBMol &box,double spacing, double pad)
{
    OBAtom *atom;
    vector<OBNodeBase*>::iterator i;

    for (atom = box.BeginAtom(i);atom;atom = box.NextAtom(i))
        if (atom->GetIdx() == 1)
        {
            _xmin = atom->GetX();
            _xmax = atom->GetX();
            _ymin = atom->GetY();
            _ymax = atom->GetY();
            _zmin = atom->GetZ();
            _zmax = atom->GetZ();
        }
        else
        {
            if (atom->GetX() < _xmin)
                _xmin = atom->GetX();
            if (atom->GetX() > _xmax)
                _xmax = atom->GetX();
            if (atom->GetY() < _ymin)
                _ymin = atom->GetY();
            if (atom->GetY() > _ymax)
                _ymax = atom->GetY();
            if (atom->GetZ() < _zmin)
                _zmin = atom->GetZ();
            if (atom->GetZ() > _zmax)
                _zmax = atom->GetZ();
        }

    _xmin -= pad;
    _xmax += pad;
    _ymin -= pad;
    _ymax += pad;
    _zmin -= pad;
    _zmax += pad;
    double midx,midy,midz;
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
    _halfSpace= _spacing/2.0;
    _inv_spa=1.0/spacing;
    _val=NULL;
    _ival=NULL;

    int size = _xdim*_ydim*_zdim;
    _val = new double [size];
    memset(_val,'\0',sizeof(double)*size);

    //return(true);
}

double OBFloatGrid::Inject(double x,double y,double z)
{
  if((x<=_xmin)||(x>=_xmax))
    return(0.0);
  if((y<=_ymin)||(y>=_ymax))
    return(0.0);
  if((z<=_zmin)||(z>=_zmax))
    return(0.0);
  
  int gx=(int)((x-_xmin)*_inv_spa);
  int gy=(int)((y-_ymin)*_inv_spa);
  int gz=(int)((z-_zmin)*_inv_spa);
  
  return(_val[(gz*_ydim*_xdim)+(gy*_xdim)+gx]);
}

void OBFloatGrid::IndexToCoords(int idx, double &x, double &y, double &z)
{
    long int grid_x,grid_y,grid_z;

    grid_x = idx % (int)_xdim;
    grid_z = (int)(idx /(_xdim * _ydim));
    grid_y = (int)((idx - (grid_z * _xdim * _ydim))/_xdim);

    x = ((double)grid_x * _spacing + _xmin) + this->_halfSpace;
    y = ((double)grid_y * _spacing + _ymin) + this->_halfSpace;
    z = ((double)grid_z * _spacing + _zmin) + this->_halfSpace;
}

int OBFloatGrid::CoordsToIndex(double &x, double &y, double &z)
{
    int gx=(int)((x-_xmin)*_inv_spa);
    int gy=(int)((y-_ymin)*_inv_spa);
    int gz=(int)((z-_zmin)*_inv_spa);

    return((gz*_ydim*_xdim)+(gy*_xdim)+gx);
}

void OBFloatGrid::CoordsToIndex(int *idx,double *c)
{
    idx[0]=(int)((c[0]-_xmin)*_inv_spa);
    idx[1]=(int)((c[1]-_ymin)*_inv_spa);
    idx[2]=(int)((c[2]-_zmin)*_inv_spa);
}

double OBFloatGrid::Interpolate(double x,double y,double z)
{
    int n,igx,igy,igz;
    double xydim;
    double gx,gy,gz,fgx,fgy,fgz;
    double ax,ay,az,bx,by,bz;
    double AyA,ByA,AyB,ByB,Az,Bz;

    if((x<=_xmin)||(x>=_xmax))
        return(0.0);
    if((y<=_ymin)||(y>=_ymax))
        return(0.0);
    if((z<=_zmin)||(z>=_zmax))
        return(0.0);

    xydim = _xdim*_ydim;

    /* calculate grid voxel and fractional offsets */
    gx=(x-_xmin-_halfSpace)*_inv_spa;
    if (gx<0)
        gx=0;
    igx=(int)gx;
    fgx=gx-(double)igx;
    gy=(y-_ymin-_halfSpace)*_inv_spa;
    if (gy<0)
        gy=0;
    igy=(int) gy;
    fgy= gy - (double) igy;
    gz=(z-_zmin-_halfSpace)*_inv_spa;
    if (gz<0)
        gz=0;
    igz=(int) gz;
    fgz= gz - (double) igz;

    n=(int)(igx+ _xdim*igy + xydim*igz);
    /* calculate linear weightings */
    ax=1.0-fgx;
    bx=fgx;
    ay=1.0-fgy;
    by=fgy;
    az=1.0-fgz;
    bz=fgz;

    /* calculate interpolated value */
    AyA=ax*_val[n           ]+bx*_val[(int)(n+1)    ];
    ByA=ax*_val[n+_xdim]+bx*_val[(int)(n+1+_xdim)  ];

    Az=ay*AyA+by*ByA;

    AyB=ax*_val[(int)(n     +xydim)]+bx*_val[(int)(n+1     +xydim)];
    ByB=ax*_val[(int)(n+_xdim+xydim)]+bx*_val[(int)(n+1+_xdim+xydim)];
    Bz=ay*AyB+by*ByB;

    return(az*Az+bz*Bz);
}


double OBFloatGrid::InterpolateDerivatives(double x,double y,double z,double *derivatives)
{
    int n,igx,igy,igz;
    double xydim;
    double gx,gy,gz,fgx,fgy,fgz;
    double ax,ay,az,bx,by,bz;
    double AyA,ByA,AyB,ByB,Az,Bz;
    double energy,fx,fy,fz;

    if((x<=_xmin)||(x>=_xmax))
        return(0.0);
    if((y<=_ymin)||(y>=_ymax))
        return(0.0);
    if((z<=_zmin)||(z>=_zmax))
        return(0.0);

    xydim = _xdim*_ydim;

    /* calculate grid voxel and fractional offsets */
    gx=(x-_xmin-_halfSpace)*_inv_spa;
    if (gx<0)
        gx=0;
    igx=(int)gx;
    fgx=gx-(double)igx;
    gy=(y-_ymin-_halfSpace)*_inv_spa;
    if (gy<0)
        gy=0;
    igy=(int) gy;
    fgy= gy - (double) igy;
    gz=(z-_zmin-_halfSpace)*_inv_spa;
    if (gz<0)
        gz=0;
    igz=(int) gz;
    fgz= gz - (double) igz;

    n=(int)(igx+ _xdim*igy + xydim*igz);
    /* calculate linear weightings */
    ax=1.0-fgx;
    bx=fgx;
    ay=1.0-fgy;
    by=fgy;
    az=1.0-fgz;
    bz=fgz;

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


ostream& operator<< ( ostream &os, const  OBFloatGrid& fg)
{
    os.write((const char*)&fg._xmin,sizeof(double));
    os.write((const char*)&fg._xmax,sizeof(double));
    os.write((const char*)&fg._ymin,sizeof(double));
    os.write((const char*)&fg._ymax,sizeof(double));
    os.write((const char*)&fg._zmin,sizeof(double));
    os.write((const char*)&fg._zmax,sizeof(double));

    os.write((const char*)&fg._midx,sizeof(double));
    os.write((const char*)&fg._midy,sizeof(double));
    os.write((const char*)&fg._midz,sizeof(double));
    os.write((const char*)&fg._inv_spa,sizeof(double));
    os.write((const char*)&fg._spacing,sizeof(double));
    os.write((const char*)&fg._xdim,sizeof(int));
    os.write((const char*)&fg._ydim,sizeof(int));
    os.write((const char*)&fg._zdim,sizeof(int));
    os.write((const char*)&fg._val[0],
             (sizeof(double)*(fg._xdim*fg._ydim*fg._zdim)));

    return(os);
}

istream& operator>> ( istream &is,OBFloatGrid& fg)
{
    is.read((char*)&fg._xmin,sizeof(double));
    is.read((char*)&fg._xmax,sizeof(double));
    is.read((char*)&fg._ymin,sizeof(double));
    is.read((char*)&fg._ymax,sizeof(double));
    is.read((char*)&fg._zmin,sizeof(double));
    is.read((char*)&fg._zmax,sizeof(double));

    is.read((char*)&fg._midx,sizeof(double));
    is.read((char*)&fg._midy,sizeof(double));
    is.read((char*)&fg._midz,sizeof(double));
    is.read((char*)&fg._inv_spa,sizeof(double));
    is.read((char*)&fg._spacing,sizeof(double));
    is.read((char*)&fg._xdim,sizeof(int));
    is.read((char*)&fg._ydim,sizeof(int));
    is.read((char*)&fg._zdim,sizeof(int));
    int size = fg._xdim*fg._ydim*fg._zdim;
    fg._val = new double [size];
    size *= (int) sizeof(double);

    is.read((char*)&fg._val[0],size);
    fg._halfSpace= fg._spacing/2.0;

    return(is);
}

#ifdef FOO
//linear interpolation routine
double Interpolate(double x,double y,double z)
{
    int n;
    int igx,igy,igz;
    double scale,gx,gy,gz,fgx,fgy,fgz;
    double qzpypx,qzpynx,qznynx,qznypx,qznyx,qzpyx;

    scale=_inv_spa;

    if((x<=_xmin)||(x>=_xmax))
        return(0.0);
    if((y<=_ymin)||(y>=_ymax))
        return(0.0);
    if((z<=_zmin)||(z>=_zmax))
        return(0.0);

    gx=(x-_xmin-_halfSpace)*scale;
    igx=(int) gx;
    fgx= gx - (double) igx;
    gy=(y-_ymin-_halfSpace)*scale;
    igy=(int) gy;
    fgy= gy - (double) igy;
    gz=(z-_zmin-_halfSpace)*scale;
    igz=(int) gz;
    fgz= gz - (double) igz;

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

} // end namespace OpenBabel

//! \file grid.cpp
//! \brief Handle grids of values.


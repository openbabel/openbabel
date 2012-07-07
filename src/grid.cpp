/**********************************************************************
grid.cpp - Handle grids of values.

Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2006 by Geoffrey R. Hutchison
Some Portions Copyright (C) 2008 by Marcus D. Hanwell

This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/
#include <openbabel/babelconfig.h>

#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/grid.h>

using namespace std;

namespace OpenBabel
{
  void OBGrid::Init(OBMol &box)
  {
    OBAtom *atom;
    vector<OBAtom*>::iterator i;

    for (atom = box.BeginAtom(i); atom; atom = box.NextAtom(i)) {
      if (atom->GetIdx() == 1) {
        _xmin = atom->GetX();
        _xmax = atom->GetX();
        _ymin = atom->GetY();
        _ymax = atom->GetY();
        _zmin = atom->GetZ();
        _zmax = atom->GetZ();
      }
      else {
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
    }
  }

  void OBFloatGrid::Init(OBMol &box, double spacing, double pad)
  {
    OBGrid::Init(box); // handle in the base class

    _xmin -= pad;
    _xmax += pad;
    _ymin -= pad;
    _ymax += pad;
    _zmin -= pad;
    _zmax += pad;

    /* store in Grid */
    _midx=0.5*(_xmax+_xmin);
    _midy=0.5*(_ymax+_ymin);
    _midz=0.5*(_zmax+_zmin);
    _xdim=static_cast<int>((_xmax-_xmin)/spacing) + 1;
    _ydim=static_cast<int>((_ymax-_ymin)/spacing) + 1;
    _zdim=static_cast<int>((_zmax-_zmin)/spacing) + 1;
    _spacing=spacing;
    _halfSpace=_spacing/2.0;
    _inv_spa=1.0/_spacing;
    _ival=NULL;

    // Calculate the size needed, resize the vector and initialise it to 0.0
    int size = _xdim*_ydim*_zdim;
    _values.resize(size, 0.0);
  }

  void OBFloatGrid::SetLimits(const double origin[3], const double v1[3], const double v2[3],
                              const double v3[3])
  {
    // Convert to vectors and call the vector form of the function
    vector3 originv(origin[0], origin[1], origin[2]);
    vector3 x(v1[0], v1[1], v1[2]);
    vector3 y(v2[0], v2[1], v2[2]);
    vector3 z(v3[0], v3[1], v3[2]);
    SetLimits(originv, x, y, z);
  }

  void OBFloatGrid::SetLimits(const vector3& origin, const vector3& x,
                              const vector3& y, const vector3& z)
  {
    // Using vectors instead of arrays of doubles
    _xmin = origin.x();
    _ymin = origin.y();
    _zmin = origin.z();

    // Set the axes
    SetXAxis(x);
    SetYAxis(y);
    SetZAxis(z);

    // Now set the max and mid as well as spacing
    _xmax = _xmin + (_xdim-1) * (x.x() + y.x() + z.x());
    _ymax = _ymin + (_ydim-1) * (x.y() + y.y() + z.y());
    _zmax = _zmin + (_zdim-1) * (x.z() + y.z() + z.z());
    // FIXME If we want support for cubes aligned along arbitrary axes then
    // averaging the x, y and z directions will not work. Just use x for now.
    _spacing = (_xmax - _xmin)/(_xdim-1);
    _halfSpace= _spacing/2.0;
    _inv_spa=1.0/_spacing;
  }

  void OBFloatGrid::SetNumberOfPoints(int nx, int ny, int nz)
  {
    _xdim = nx;
    _ydim = ny;
    _zdim = nz;
    _values.resize(nx*ny*nz);
  }

  void OBFloatGrid::SetXAxis(vector3 v)
  {
    _xAxis = v;
  }

  void OBFloatGrid::SetYAxis(vector3 v)
  {
    _yAxis = v;
  }

  void OBFloatGrid::SetZAxis(vector3 v)
  {
    _zAxis = v;
  }

  double OBFloatGrid::Inject(double x, double y, double z)
  {
    if (_values.empty())
      return 0.0;

    if( x<=_xmin || x>=_xmax
        || y<=_ymin || y>=_ymax
        || z<=_zmin || z>=_zmax ) return 0.0;

    return(_values[CoordsToIndex(x, y, z)]);
  }

  void OBFloatGrid::IndexToCoords(int idx, double &x, double &y, double &z)
  {
    int grid_x, grid_y, grid_z;

    grid_z = idx % _zdim;
    grid_x = static_cast<int>(idx / (_ydim * _zdim));
    grid_y = static_cast<int>((idx - (grid_x * _ydim * _zdim)) / _zdim);

    x = (grid_x * _spacing + _xmin) + _halfSpace;
    y = (grid_y * _spacing + _ymin) + _halfSpace;
    z = (grid_z * _spacing + _zmin) + _halfSpace;
  }

  int OBFloatGrid::CoordsToIndex(double x, double y, double z)
  {
    int gx=static_cast<int>((x-_xmin)*_inv_spa);
    int gy=static_cast<int>((y-_ymin)*_inv_spa);
    int gz=static_cast<int>((z-_zmin)*_inv_spa);

    return((gx*_ydim*_zdim) + (gy*_zdim) + gz);
  }

  void OBFloatGrid::CoordsToIndex(int *idx,double *c)
  {
    idx[0]=static_cast<int>((c[0]-_xmin)*_inv_spa);
    idx[1]=static_cast<int>((c[1]-_ymin)*_inv_spa);
    idx[2]=static_cast<int>((c[2]-_zmin)*_inv_spa);
  }

  double OBFloatGrid::Interpolate(double x, double y, double z)
  {
    if (_values.empty())
      return 0.0;

    int n,igx,igy,igz;
    double yzdim;
    double gx,gy,gz,fgx,fgy,fgz;
    double ax,ay,az,bx,by,bz;
    double AyA,ByA,AyB,ByB,Az,Bz;

    if( x<=_xmin || x>=_xmax
        || y<=_ymin || y>=_ymax
        || z<=_zmin || z>=_zmax ) return 0.0;

    yzdim = _ydim*_zdim;

    /* calculate grid voxel and fractional offsets */
    gx=(x-_xmin-_halfSpace)*_inv_spa;
    if (gx<0)
      gx=0;
    igx=static_cast<int>(gx);
    fgx=gx-static_cast<double>(igx);
    gy=(y-_ymin-_halfSpace)*_inv_spa;
    if (gy<0)
      gy=0;
    igy=static_cast<int>(gy);
    fgy= gy - static_cast<double>(igy);
    gz=(z-_zmin-_halfSpace)*_inv_spa;
    if (gz<0)
      gz=0;
    igz=static_cast<int>(gz);
    fgz= gz - static_cast<double>(igz);

    // Calculate the index of the nearest point
    n=static_cast<int>(igx*yzdim + igy*_zdim + igz);

    if ((n+1+_zdim+yzdim) >= (yzdim*_xdim))
      return 0.0;

    /* calculate linear weightings */
    ax=1.0-fgx;
    bx=fgx;
    ay=1.0-fgy;
    by=fgy;
    az=1.0-fgz;
    bz=fgz;

    /* calculate interpolated value */
    AyA=az*_values[n           ]+bz*_values[(n+1)    ];
    ByA=az*_values[n+_zdim]+bz*_values[(n+1+_zdim)  ];

    Az=ay*AyA+by*ByA;

    AyB=az*_values[static_cast<int>(n     +yzdim)] +
        bz*_values[static_cast<int>(n+1     +yzdim)];
    ByB=az*_values[static_cast<int>(n+_zdim+yzdim)] +
        bz*_values[static_cast<int>(n+1+_zdim+yzdim)];
    Bz=ay*AyB+by*ByB;

    return(ax*Az+bx*Bz);
  }

  std::vector<double> OBFloatGrid::GetDataVector()
  {
    return _values;
  }

  void OBFloatGrid::SetVals(const std::vector<double> & vals)
  {
    _values.clear();
    _values = vals;
  }

  double OBFloatGrid::InterpolateDerivatives(double x,double y,double z,double *derivatives)
  {
    int n,igx,igy,igz;
    double yzdim;
    double gx,gy,gz,fgx,fgy,fgz;
    double ax,ay,az,bx,by,bz;
    double AyA,ByA,AyB,ByB,Az,Bz;
    double energy,fx,fy,fz;

    if( x<=_xmin || x>=_xmax
        || y<=_ymin || y>=_ymax
        || z<=_zmin || z>=_zmax ) return 0.0;

    yzdim = _ydim*_zdim;

    /* calculate grid voxel and fractional offsets */
    gx=(x-_xmin-_halfSpace)*_inv_spa;
    if (gx<0)
      gx=0;
    igx=static_cast<int>(gx);
    fgx=gx-(double)igx;
    gy=(y-_ymin-_halfSpace)*_inv_spa;
    if (gy<0)
      gy=0;
    igy=static_cast<int>(gy);
    fgy= gy - (double) igy;
    gz=(z-_zmin-_halfSpace)*_inv_spa;
    if (gz<0)
      gz=0;
    igz=static_cast<int>(gz);
    fgz= gz - (double) igz;

    n=static_cast<int>(igx*yzdim + igy*_zdim + igz);
    /* calculate linear weightings */
    ax=1.0-fgx;
    bx=fgx;
    ay=1.0-fgy;
    by=fgy;
    az=1.0-fgz;
    bz=fgz;

    /* calculate interpolated value */
    AyA=az*_values[n           ]+bx*_values[n+1    ];
    ByA=az*_values[n+_zdim]+bx*_values[(n+1+_zdim)  ];

    Az=ay*AyA+by*ByA;

    AyB=az*_values[static_cast<int>(n      +yzdim)]+
        bz*_values[static_cast<int>(n+1      +yzdim)];
    ByB=az*_values[static_cast<int>(n+_zdim+yzdim)]+
        bz*_values[static_cast<int>(n+1+_zdim+yzdim)];
    Bz=ay*AyB+by*ByB;

    energy = ax*Az+bx*Bz;

    //calculate derivatives

    fz=-Az+Bz;

    Az=-AyA+ByA;
    Bz=-AyB+ByB;

    fy=az*Az+bz*Bz;

    AyA=-_values[n           ]+_values[(n+1)     ];
    ByA=-_values[n+_zdim      ]+_values[(n+1+_zdim)];

    Az=ay*AyA+by*ByA;

    AyB=-_values[static_cast<int>(n      +yzdim)]+
         _values[static_cast<int>(n+1      +yzdim)];
    ByB=-_values[static_cast<int>(n+_zdim+yzdim)]+
         _values[static_cast<int>(n+1+_zdim+yzdim)];

    Bz=ay*AyB+by*ByB;

    fx=ax*Az+bx*Bz;

    derivatives[0] += fx;
    derivatives[1] += fy;
    derivatives[2] += fz;

    return(energy);
  }

  ostream& operator<< ( ostream &os, const  OBFloatGrid& fg)
  {
    //@todo: this code stores the data in way that depends on
    // the bit representation of floating-point numbers. One can say
    // it's OK because IEEE754 is a widely accepted standard, but then
    // one should at least make it write the data in an endianness-
    // independent way. Of course this comment applies also to
    // the operator>> below.
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
    os.write((const char*)&fg._values[0],
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
    fg._values.resize(size);
    size *= (int) sizeof(double);

    is.read((char*)&fg._values[0],size);
    fg._halfSpace= fg._spacing/2.0;

    return(is);
  }

  void OBProxGrid::Setup(OBMol &mol,OBMol &box,double cutoff,double res)
  {
    this->Init(box); // handle in the base class

    _inc = res; // 1/2 angstrom resolution

    _nxinc = (int) floor((_xmax - _xmin)/0.5);
    _nyinc = (int) floor((_ymax - _ymin)/0.5);
    _nzinc = (int) floor((_zmax - _zmin)/0.5);
    _maxinc = _nxinc*_nyinc*_nzinc;

    int j,size = _nxinc*_nyinc*_nzinc;
    cell.resize(size);
    for (unsigned int num = 0; num < cell.size(); ++num)
      cell[num].resize(0);

    cutoff *= cutoff; //don't do sqrts

    int k,l,m;
    double x,y,z,dx_2,dy_2;
    double *c = mol.GetCoordinates();
    size = mol.NumAtoms()*3;

    OBAtom *atom;
    vector<OBAtom*>::iterator i;
    for (atom = mol.BeginAtom(i),j=0;atom;atom = mol.NextAtom(i),j+=3)
      if (PointIsInBox(c[j],c[j+1],c[j+2]))
        for (x = _xmin+(_inc/2.0),k=0;k < _nxinc;x+=_inc,++k)
          if ((dx_2 = SQUARE(c[j]-x)) < cutoff)
            for (y = _ymin+(_inc/2.0),l=0;l < _nyinc;y+=_inc,++l)
              if ((dx_2+(dy_2 = SQUARE(c[j+1]-y))) < cutoff)
                for (z = _zmin+(_inc/2.0),m=0;m < _nzinc;z+=_inc,++m)
                  if ((dx_2+dy_2+SQUARE(c[j+2]-z)) < cutoff)
                    cell[(k*_nyinc*_nzinc)+(l*_nzinc)+m].push_back(atom->GetIdx());

    _inc = 1/_inc;
  }

  void OBProxGrid::Setup(OBMol &mol,OBMol &box,double cutoff,vector<bool> &use,
                         double res)
  {
    this->Init(box); // handle in the base class

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

    OBAtom *atom;
    vector<OBAtom*>::iterator i;
    for (atom = mol.BeginAtom(i),j=0;atom;atom = mol.NextAtom(i),j+=3)
      if (use[atom->GetIdx()])
        if (PointIsInBox(c[j],c[j+1],c[j+2]))
          for (x = _xmin+(_inc/2.0),k=0;k < _nxinc;x+=_inc,++k)
            if ((dx_2 = SQUARE(c[j]-x)) < cutoff)
              for (y = _ymin+(_inc/2.0),l=0;l < _nyinc;y+=_inc,++l)
                if ((dx_2+(dy_2 = SQUARE(c[j+1]-y))) < cutoff)
                  for (z = _zmin+(_inc/2.0),m=0;m < _nzinc;z+=_inc,++m)
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

} // end namespace OpenBabel

//! \file grid.cpp
//! \brief Handle grids of values.


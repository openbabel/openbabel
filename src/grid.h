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

#ifndef OB_GRID_H
#define OB_GRID_H

#ifdef __sgi
#include <iostream.h>
#else
#include <iostream>
#endif

#include <algorithm>
#include <vector>
#include <string>

#ifndef OBPolarGrid
#define OBPolarGrid 0x01
#endif //OBPolarGrid

#ifndef OBLipoGrid
#define OBLipoGrid 0x02
#endif //OBLipoGrid

namespace OpenBabel {

class OBProxGrid
{
  int _gridtype;
  int _nxinc,_nyinc,_nzinc,_maxinc;
  double _xmin,_xmax,_ymin,_ymax,_zmin,_zmax,_inc;
  std::vector<std::vector<int> > cell;

public:

  OBProxGrid(int gridtype=0){_gridtype=gridtype;}
  ~OBProxGrid(){}
  void Setup(OBMol &,OBMol &,double,double res=0.5);
  void Setup(OBMol &,OBMol &,double,std::vector<bool>&,double res=0.5);
  std::vector<int> *GetProxVector(double,double,double);
  std::vector<int> *GetProxVector(double*);
  // For HasFlag force return type to bool so VC6.0 doesn't complain
  bool LipoGrid() {return((_gridtype&OBLipoGrid) ? true : false);}
  bool PolarGrid() {return(_gridtype&OBPolarGrid);}
  void SetGridType(int gridtype) {_gridtype = gridtype;}
  bool PointIsInBox(double x,double y,double z)
    {
      return (x>=_xmin) && (x<=_xmax) && 
	(y>=_ymin) && (y<=_ymax) &&
	(z>=_zmin) && (z<=_zmax);
    }
};

class OBFloatGrid
{
  double *_val;             //!< doubleing point values
  int   *_ival;            //!< for integer values
  double _midz,_midx,_midy;   //!< center of grid in world coordinates
  int _ydim,_xdim,_zdim;     //!< grid dimensions
  double _spacing,_inv_spa;  //!< spacing between grid points and its invers
  double _xmin,_xmax,_ymin,_ymax,_zmin,_zmax;
  double _halfSpace;         //!< half of the grid spacing */     

public:

    OBFloatGrid() : _halfSpace(0.0) {_val=NULL;_ival=NULL;}
    ~OBFloatGrid() {if (_ival) {delete [] _ival; _ival = NULL;}
                    if (_val) {delete [] _val; _val = NULL;}}
  void Init(OBMol &,double, double pad= 0.0); //!< initialized using boxfile
  bool PointIsInBox(double x,double y,double z)
    {
      return (x>=_xmin) && (x<=_xmax) && 
	(y>=_ymin) && (y<=_ymax) &&
	(z>=_zmin) && (z<=_zmax);
    }
  bool PointIsInBox(double *c)
    {
      return (c[0]>=_xmin) && (c[0]<=_xmax) && 
	(c[1]>=_ymin) && (c[1]<=_ymax) &&
	(c[2]>=_zmin) && (c[2]<=_zmax);
    }
  double GetXmin() const {return(_xmin);}
  double GetYmin() const {return(_ymin);}
  double GetZmin() const {return(_zmin);}
  double GetXmax() const {return(_xmax);}
  double GetYmax() const {return(_ymax);}
  double GetZmax() const {return(_zmax);}
  double GetSpacing() const {return(_spacing);}
  double GetScale() const {return(_inv_spa);}
  double GetHalfSpace() const {return(_halfSpace);}
  int GetXdim() const {return(_xdim);}
  int GetYdim() const {return(_ydim);}
  int GetZdim() const {return(_zdim);}
  void GetMin(double *a) {a[0]=_xmin;a[1]=_ymin;a[2]=_zmin;}
  void GetMax(double *a) {a[0]=_xmax;a[1]=_ymax;a[2]=_zmax;}
  void GetDim(int *a)   {a[0]=_xdim;a[1]=_ydim;a[2]=_zdim;}
  void GetSpacing(double &s) {s=_spacing;}
  vector3 GetMidpointVector() 
    {vector3 v; v.Set(_midx,_midy,_midz); return(v);}
  double *GetVals() {return(_val);}
  void SetVals(double *ptr) {_val = ptr;}
  vector3 Center() { return vector3(_midx,_midy,_midz); } //added by jjc
  friend std::ostream& operator<< ( std::ostream&, const OBFloatGrid& ) ;
  friend std::istream& operator>> ( std::istream&,OBFloatGrid& ) ;

  double Inject(double x,double y,double z)
    {
      if((x<=_xmin)||(x>=_xmax)) return(0.0);
      if((y<=_ymin)||(y>=_ymax)) return(0.0);
      if((z<=_zmin)||(z>=_zmax)) return(0.0);

      int gx=(int)((x-_xmin)*_inv_spa);
      int gy=(int)((y-_ymin)*_inv_spa);
      int gz=(int)((z-_zmin)*_inv_spa);

      return(_val[(gz*_ydim*_xdim)+(gy*_xdim)+gx]);
    }


  void IndexToCoords(int idx, double &x, double &y, double &z);
  void CoordsToIndex(int*,double*);
  int CoordsToIndex(double &x, double &y, double &z);
  double Interpolate(double,double,double);
  double InterpolateDerivatives(double,double,double,double *derivatives);
};

typedef enum { Undefined = -1, PLP, ChemScore } score_t;

class OBScoreGrid
{
protected:

  score_t gridtype;
  bool verbose;

public:

  double score;

  OBScoreGrid(void) { verbose = false; }
  virtual ~OBScoreGrid(void) {}

  void    SetVerbose(bool v)    { verbose = v; }
  void    SetType(score_t type) { gridtype = type; }
  score_t GetType(void)         { return gridtype; }

  virtual void   Clear(void) { }
  virtual double  Eval(double *)    { return -1; }
  virtual double  Eval(OBMol &mol) { return Eval(mol.GetCoordinates()); }
  virtual void   Init(OBMol &, OBMol &, std::string &, double) {}
  virtual void   Setup(OBMol &) {}
  virtual void   Setup(OBMol &, std::vector<int> &) {}
  virtual void   Setup(std::vector<int> &) {}
  virtual void   Config(std::string) {}
  virtual bool   Read(std::string)       { return false; }
  virtual bool   Write(std::string)      { return false; }
  virtual vector3 Center()           { return VZero; }
  virtual vector3 CenterMol(OBMol &) { return VZero; }
};

} // end namespace OpenBabel

#endif // OB_GRID_H

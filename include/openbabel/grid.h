/**********************************************************************
grid.h - Handle grids of values.

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

#ifndef OB_GRID_H
#define OB_GRID_H

#include <openbabel/babelconfig.h>
#include <openbabel/math/vector3.h>
#include <openbabel/base.h>

// Necessary evil for 2.x series -- should use OBMol* below
#include <openbabel/mol.h>

#include <iosfwd>
#include <algorithm>
#include <vector>
#include <string>

namespace OpenBabel
{

  // Forward declaration
  class OBMol;

  //! \class OBGrid grid.h <openbabel/grid.h>
  //! \brief A base grid class
 class OBAPI OBGrid: public OBBase
  {
  protected:
    double _xmin,_xmax,_ymin,_ymax,_zmin,_zmax; //!< the min/max axes in XYZ axes (i.e., the box)

  public:
    OBGrid() {}

    //! \brief Initialize the grid based on a box around the molecule @p box
    //! Subclasses should overload this method -- this only tracks the
    //! dimension of the box itself
    virtual void Init(OBMol &box);

    //! \return the minimum x point of the grid
    double GetXmin() const    { return(_xmin);    }
    //! \return the minimum y point of the grid
    double GetYmin() const    { return(_ymin);    }
    //! \return the minimum z point of the grid
    double GetZmin() const    { return(_zmin);    }
    //! \return the maximum x point of the grid
    double GetXmax() const    { return(_xmax);    }
    //! \return the maximum y point of the grid
    double GetYmax() const    { return(_ymax);    }
    //! \return the maximum z point of the grid
    double GetZmax() const    { return(_zmax);    }

    //! \return whether the supplied XYZ coordinates fall within the box
    bool PointIsInBox(double x,double y,double z)
    {
      return (x>=_xmin) && (x<=_xmax) &&
        (y>=_ymin) && (y<=_ymax) &&
        (z>=_zmin) && (z<=_zmax);
    }
    //! \return true if the point falls within the box
    bool PointIsInBox(double *c)
    {
      return (c[0]>=_xmin) && (c[0]<=_xmax) &&
        (c[1]>=_ymin) && (c[1]<=_ymax) &&
        (c[2]>=_zmin) && (c[2]<=_zmax);
    }

    //! \return true if the point falls within the box
    bool PointIsInBox(vector3 v)
    {
      return (v.x() >= _xmin) && (v.x() <=_xmax) &&
      (v.y()>=_ymin) && (v.y()<=_ymax) &&
      (v.z()>=_zmin) && (v.z()<=_zmax);
    }
  };

  //! \class OBFloatGrid grid.h <openbabel/grid.h>
  //! \brief Handle double precision floating point 3D grids (e.g., charge density around an OBMol)
  //!
  //! Supports input/output and base functionality for simple 3D discrete grids
  //! of some function -- typically around a molecule. Typically you will want
  //! to use OBGridData which uses OBFloatGrid to store its data.
  //! \sa OBGridData
 class OBAPI OBFloatGrid: public OBGrid
  {
  protected:
    std::vector<double> _values;   //!< floating point values
    int   *_ival;             //!< for integer values (deprecated)
    double _midz,_midx,_midy; //!< center of grid in world coordinates
    int _ydim,_xdim,_zdim;    //!< grid dimensions
    double _spacing,_inv_spa; //!< spacing between grid points and its inverse
    double _halfSpace;        //!< half of the grid spacing
    //! Three axes (i.e., translation vectors like a unit cell)
    vector3 _xAxis, _yAxis, _zAxis;

  public:

    OBFloatGrid() : _ival(nullptr), _halfSpace(0.0) {}
    ~OBFloatGrid()
    {
      if (_ival) delete [] _ival;
    }

    //! Initialize the grid using this molecule as a box (plus a padding)
    //! with the supplied spacing between points.
    void Init(OBMol &box,double spacing, double pad=0.0);

    //! \return the minimum point in the grid.
    vector3 GetMin() { return vector3(_xmin, _ymin, _zmin); }

    //! Get the minimum point in the grid.
    //! \deprecated Will be removed.
    //! Use \sa GetMin()
    void GetMin(double *a)
    {
      a[0]=_xmin;
      a[1]=_ymin;
      a[2]=_zmin;
    }

    //! \return the minimum point in the grid.
    vector3 GetMax() { return vector3(_xmax, _ymax, _zmax); }

    //! Get the maximum point in the grid.
    //! \deprecated Will be removed.
    //! \sa GetMax()
    void GetMax(double *a)
    {
      a[0]=_xmax;
      a[1]=_ymax;
      a[2]=_zmax;
    }

    //! \return The grid spacing.
    double GetSpacing() const { return(_spacing); }
    //! Get the grid spacing.
    //! \deprecated Will be removed.
    //! \sa GetSpacing()
    void GetSpacing(double &s)
    {
      s=_spacing;
    }
    //! \return Inverse of the grid spacing.
    double GetScale() const   { return(_inv_spa); }
    //! \return Half of the spacing between grid points.
    double GetHalfSpace() const {return(_halfSpace);}
    //! \return Size of the grid in the x dimension.
    int GetXdim() const       { return(_xdim);    }
    //! \return Size of the grid in the y dimension.
    int GetYdim() const       { return(_ydim);    }
    //! \return Size of the grid in the z dimension.
    int GetZdim() const       { return(_zdim);    }
    //! Get the x, y and z dimensions (must pass an double[3] at least).
    //! \deprecated May be removed in future.
    //! \sa GetXdim() \sa GetYdim() \sa GetZdim()
    void GetDim(int *a)
    {
      a[0]=_xdim;
      a[1]=_ydim;
      a[2]=_zdim;
    }

    //! \return Position of the center of the grid.
    vector3 GetMidpointVector()
    {
      vector3 v;
      v.Set(_midx,_midy,_midz);
      return(v);
    }

    //! \return X axis direction.
    vector3 GetXAxis() const
    {
      return _xAxis;
    }

    //! \return Y axis direction.
    vector3 GetYAxis() const
    {
      return _yAxis;
    }

    //! \return Z axis direction.
    vector3 GetZAxis() const
    {
      return _zAxis;
    }

    //! Sets the number of points in the x, y and z directions.
    void SetNumberOfPoints(int nx, int ny, int nz);

    //! Set the direction of the x axis.
    void SetXAxis(vector3);
    //! Set the direction of the y axis.
    void SetYAxis(vector3);
    //! Set the direction of the z axis.
    void SetZAxis(vector3);

    //! Set the limits (i.e., the origin point and the axes)
    //! NOTE: You must set the number of points first,
    //!       with SetNumberOfPoints
    //!       so the grid spacing can be calculated
    void SetLimits(const vector3& origin, const vector3& x, const vector3& y,
                   const vector3& z);
    //! \deprecated Will be removed.
    //! \sa SetLimits(const vector3& origin, const vector3& x, const vector3& y, const vector3& z)
    void SetLimits(const double origin[3], const double x[3], const double y[3],
                   const double z[3]);

    //! Get a copy of the vector that stores the points in the grid.
    std::vector<double> GetDataVector();
    //! Set the values in the grid to those in the vector passed. Note that the
    //! vector must be of the same dimensions as the grid based on the values
    //! given in SetNumberOfPoints(int nx, int ny, int nz).
    void SetVals(const std::vector<double> & vals);

    //! \return Pointer to the first element of the grid point data stored as a
    //! one dimensional array.
    //! \deprecated Will be removed.
    //! \sa GetDataVector()
    double *GetVals()    {        return(&_values[0]);    }

    //! \return Value at the point in the grid specified by i, j and k.
    double GetValue(int i, int j, int k)
    {
      if (i*_ydim*_zdim + j*_zdim + k > _xdim*_ydim*_zdim)
        return 0.0;
      else
        return _values[i*_ydim*_zdim + j*_zdim + k];
    }

    //! \deprecated Will be removed.
    //! \sa SetVals(const std::vector<double> & vals)
    void SetVals(double *ptr)
    {
     for (int i = 0; i < _xdim*_ydim*_zdim; ++i)
       _values[i] = ptr[i];
    }

    //! Set the value at the grid point specified by i, j and k to val.
    bool SetValue(int i, int j, int k, double val)
    {
      if (i*_ydim*_zdim + j*_zdim + k > _xdim*_ydim*_zdim)
        return false;

      _values[i*_ydim*_zdim + j*_zdim + k] = val;
      return true;
    }

    //! \return Position of the center of the grid.
    vector3 Center()
    {
      return vector3(_midx,_midy,_midz);
    }

    friend std::ostream& operator<< ( std::ostream&, const OBFloatGrid& ) ;
    friend std::istream& operator>> ( std::istream&,OBFloatGrid& ) ;

    //! \return the value at the given point (rounding as needed)
    double Inject(double x,double y,double z);
    void IndexToCoords(int idx, double &x, double &y, double &z);
    void CoordsToIndex(int*,double*);
    int CoordsToIndex(double x, double y, double z);
    //! \return the interpolated value for the given point
    double Interpolate(double,double,double);
    //! \return the interpolated value for the given point and set the dx, dy, dz derivatives
    double InterpolateDerivatives(double,double,double,double *derivatives);
  };

#ifndef OBPolarGrid
#define OBPolarGrid 0x01 /* polar interactions? */
#endif //OBPolarGrid

#ifndef OBLipoGrid
#define OBLipoGrid 0x02 /* lipophilicity? */
#endif //OBLipoGrid

  //! \class OBProxGrid grid.h <openbabel/grid.h>
  //! \brief A grid for determining the proximity of a given point to atoms in an OBMol
  //! \deprecated May be removed in the future, since docking is not a key feature
 class OBAPI OBProxGrid: public OBGrid
  {
  protected:
    int _gridtype;
    int _nxinc,_nyinc,_nzinc,_maxinc;
    double _inc;
    std::vector<std::vector<int> > cell;

  public:

    OBProxGrid(int gridtype=0)
      {
        _gridtype=gridtype;
      }
    ~OBProxGrid()
      {}
    void Setup(OBMol &mol,OBMol &box, double cutoff,double resolution = 0.5);
    void Setup(OBMol &mol,OBMol &box, double cutoff,
               std::vector<bool> &use,double resolution = 0.5);
    std::vector<int> *GetProxVector(double,double,double);
    std::vector<int> *GetProxVector(double*);

    bool LipoGrid()
    {
      return((_gridtype&OBLipoGrid) ? true : false);
    }
    bool PolarGrid()
    {
      return(_gridtype&OBPolarGrid);
    }
    void SetGridType(int gridtype)
    {
      _gridtype = gridtype;
    }
  };

  // scoring function used: PLP = Piecewise Linear Potential or ChemScore algorithm
  typedef enum { Undefined = -1, PLP, ChemScore } score_t;

  //! \class OBScoreGrid grid.h <openbabel/grid.h>
  //! \brief A base class for scoring docking interactions between multiple molecules
  //! \deprecated Will disappear in future versions. Use your own code.
  class OBAPI OBScoreGrid
  {
  protected:
    score_t gridtype;
    bool verbose;

  public:

    double score;

    OBScoreGrid(void)                 {  verbose = false;      }
    virtual ~OBScoreGrid(void) {}

    void    SetVerbose(bool v)        {      verbose = v;      }
    void    SetType(score_t type)     {      gridtype = type;  }
    score_t GetType(void)             {    return gridtype;    }

    virtual void   Clear(void)        { }
    virtual double  Eval(double *)    {       return -1;       }
    virtual double  Eval(OBMol &mol){return Eval(mol.GetCoordinates());}
    virtual void   Init(OBMol &, OBMol &, std::string &, double){}
    virtual void   Setup(OBMol &) {}
    virtual void   Setup(OBMol &, std::vector<int> &){}
    virtual void   Setup(std::vector<int> &) {}
    virtual void   Config(std::string) {}
    virtual bool   Read(std::string)   {      return false;    }
    virtual bool   Write(std::string)  {      return false;    }
    virtual vector3 Center()           {      return VZero;    }
    virtual vector3 CenterMol(OBMol &) {      return VZero;    }
  };

} // end namespace OpenBabel

#endif // OB_GRID_H

//! \file grid.h
//! \brief Handle grids of values.

/**********************************************************************
vector3.h - Handle 3D coordinates.
 
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2006 by Geoffrey R. Hutchison
Some portions Copyright (C) 2006 by Benoit Jacob
 
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

#ifndef OB_VECTOR_H
#define OB_VECTOR_H

#include <ostream>
#include <vector>
#include <math.h>
#include <iostream>

#include <Eigen/Core>
#include <Eigen/Geometry> 
#include <Eigen/Array> // for Vector3d::random()
#include <Eigen/LU> // for Matrix3d::deteminant()

#ifndef RAD_TO_DEG
#define RAD_TO_DEG (180.0/M_PI)
#endif

#ifndef DEG_TO_RAD
#define DEG_TO_RAD (M_PI/180.0)
#endif

namespace OpenBabel
{
  /*! Calculate the angle between vectors
   *  \param ab First vector (ab = a - b)
   *  \param bc Second vector (bc = c - b)
   *  \return The angle abc in degrees
   */
  OBAPI double VectorAngle(const Eigen::Vector3d& ab, const Eigen::Vector3d& bc);
  
  /*! Calculate the torsion angle between vectors
   *  \param a Atom a (coordinates)
   *  \param b Atom b (coordinates)
   *  \param c Atom c (coordinates)
   *  \param d Atom d (coordinates)
   *  \return The torson angle abcd in degrees
   */
  OBAPI double VectorTorsion(const Eigen::Vector3d& a, const Eigen::Vector3d& b, 
      const Eigen::Vector3d& c, const Eigen::Vector3d& d);
 
  /*! Calculate the OOP angle a-b-c-d. b is the central atom, and a-b-c is the plane. 
   *  The OOP angle is given by 90° - arccos(dot(corss(ab,cb),db)/rabbc*rdb).
   *  \param a Atom a (coordinates)
   *  \param b Atom b (coordinates)
   *  \param c Atom c (coordinates)
   *  \param d Atom d (coordinates)
   *  \return The out-of-plane angle abcd in degrees
   */
  OBAPI double VectorOOP(const Eigen::Vector3d &a, const Eigen::Vector3d &b, 
      const Eigen::Vector3d &c, const Eigen::Vector3d &d);
 
  //! Calculate the distance of point a to the plane determined by b,c,d
  OBAPI double Point2Plane(const Eigen::Vector3d& a, const Eigen::Vector3d& b, 
      const Eigen::Vector3d& c, const Eigen::Vector3d& d);
  
  //! Calculate the angle between point a and the plane determined by b,c,d (in degrees)
  OBAPI double Point2PlaneAngle(const Eigen::Vector3d& a, const Eigen::Vector3d& b, 
      const Eigen::Vector3d& c, const Eigen::Vector3d& d);

  /*! Calculate the derivative of a vector length. The vector is given by a - b, 
   *  the length of this vector rab = sqrt(ab.x^2 + ab.y^2 + ab.z^2).
   *  \param a Atom a (coordinates)
   *  \param b Atom b (coordinates)
   *  \param Fa - Return value for the force on atom a
   *  \param Fb - Return value for the force on atom b
   *  \return The distance between a and b (bondlength for bond stretching, separation for vdw, electrostatic)
   */
  OBAPI double VectorBondDerivative(const Eigen::Vector3d &a, const Eigen::Vector3d &b, 
      Eigen::Vector3d &Fa, Eigen::Vector3d &Fb);

  /*! To be used for VDW or Electrostatic interactions. This
   *  is faster than VectorBondDerivative, but does no error checking. 
   */
  OBAPI double VectorDistanceDerivative(const Eigen::Vector3d &a, const Eigen::Vector3d &b, 
      Eigen::Vector3d &Fa, Eigen::Vector3d &Fb);
 
  /*! Calculate the derivative of a angle a-b-c. The angle is given by acos(dot(ab,cb)/rab*rcb).
   *  Used for harmonic (cubic) angle potentials.
   *  \param a Atom a (coordinates)
   *  \param b Atom b (coordinates)
   *  \param c Atom c (coordinates)
   *  \param Fa - Return value for the force on atom a
   *  \param Fb - Return value for the force on atom b
   *  \param Fc - Return value for the force on atom c
   *  \return The angle abc in degrees
   */
  OBAPI double VectorAngleDerivative(const Eigen::Vector3d &a, const Eigen::Vector3d &b, const Eigen::Vector3d &c,
      Eigen::Vector3d &Fa, Eigen::Vector3d &Fb, Eigen::Vector3d &Fc);
 
  /*! Calculate the derivative of a torsion angle a-b-c-d. The torsion angle is given by 
   *  arccos(dot(corss(ab,bc),cross(bc,cd))/rabbc*rbccd).
   *  \param a Atom a (coordinates)
   *  \param b Atom b (coordinates)
   *  \param c Atom c (coordinates)
   *  \param d Atom d (coordinates)
   *  \param Fa - Return value for the force on atom a
   *  \param Fb - Return value for the force on atom b
   *  \param Fc - Return value for the force on atom c
   *  \param Fd - Return value for the force on atom d
   *  \return The tosion angle for atoms a-b-c-d
   */
  OBAPI double VectorTorsionDerivative(const Eigen::Vector3d &a, const Eigen::Vector3d &b, const Eigen::Vector3d &c, 
      const Eigen::Vector3d &d, Eigen::Vector3d &Fa, Eigen::Vector3d &Fb, Eigen::Vector3d &Fc, Eigen::Vector3d &Fd);

  /*! Calculate the derivative of a OOP angle a-b-c-d. b is the central atom, and a-b-c is the plane. 
   *  The OOP angle is given by 90° - arccos(dot(corss(ab,cb),db)/rabbc*rdb).
   *  \param a Atom a (coordinates)
   *  \param b Atom b (coordinates)
   *  \param c Atom c (coordinates)
   *  \param d Atom d (coordinates)
   *  \param Fa - Return value for the force on atom a
   *  \param Fb - Return value for the force on atom b
   *  \param Fc - Return value for the force on atom c
   *  \param Fd - Return value for the force on atom d
   * \return The OOP angle for a-b-c-d
   */
  OBAPI double VectorOOPDerivative(const Eigen::Vector3d &a, const Eigen::Vector3d &b, const Eigen::Vector3d &c, 
      const Eigen::Vector3d &d, Eigen::Vector3d &Fa, Eigen::Vector3d &Fb, Eigen::Vector3d &Fc, Eigen::Vector3d &Fd);
  
  OBAPI void SetTorsion(double *c, int ref[4], double setang, std::vector<int> atoms);
  OBAPI void SetTorsion(double *c, unsigned int ref[4], double setang, std::vector<int> atoms);
 
  //  The global constant Eigen::Vector3d objects
  //! The zero vector: <0.0, 0.0, 0.0>
  extern OBAPI const Eigen::Vector3d VZero;
  //! The x unit vector: <1.0, 0.0, 0.0>
  extern OBAPI const Eigen::Vector3d VX;
  //! The y unit vector: <0.0, 1.0, 0.0>
  extern OBAPI const Eigen::Vector3d VY;
  //! The z unit vector: <0.0, 0.0, 1.0>
  extern OBAPI const Eigen::Vector3d VZ;
}

#endif // OB_VECTOR_H

//! \file
//! \brief Handle 3D coordinates.

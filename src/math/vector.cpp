/**********************************************************************
vector3.cpp - Handle 3D coordinates.
 
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

#include <openbabel/babelconfig.h>

#include <iostream>

#include <openbabel/math/vector.h>
#include <openbabel/obutil.h>

using namespace std;

namespace OpenBabel
{
  const Eigen::Vector3d VZero ( 0.0, 0.0, 0.0 ) ;
  const Eigen::Vector3d VX    ( 1.0, 0.0, 0.0 ) ;
  const Eigen::Vector3d VY    ( 0.0, 1.0, 0.0 ) ;
  const Eigen::Vector3d VZ    ( 0.0, 0.0, 1.0 ) ;

  OBAPI double VectorBondDerivative(const Eigen::Vector3d &a, const Eigen::Vector3d &b, 
      Eigen::Vector3d &Fa, Eigen::Vector3d &Fb)
  {
    Eigen::Vector3d ab = a - b;
    double rab = ab.norm();
 
    if (rab < 0.1) // atoms are too close to each other
      rab = 0.1;

    Fb = ab * (1.0 / rab);
    Fa = -Fb;
    
    return rab;
  }

  OBAPI double VectorDistanceDerivative(const Eigen::Vector3d &a, const Eigen::Vector3d &b,
      Eigen::Vector3d &Fa, Eigen::Vector3d &Fb)
  {
    const Eigen::Vector3d ab = a - b;
    const double rab = ab.norm();
    //const double inverse_rab = 1.0 / rab;
    
    //Fb = ab * inverse_rab;
    Fb = ab * (1.0 / rab);
    Fa = -Fb;
 
    return rab;
  }
  
  OBAPI double VectorAngleDerivative(const Eigen::Vector3d &a, const Eigen::Vector3d &b, const Eigen::Vector3d &c,
      Eigen::Vector3d &Fa, Eigen::Vector3d &Fb, Eigen::Vector3d &Fc)
  {
    // This is adapted from http://scidok.sulb.uni-saarland.de/volltexte/2007/1325/pdf/Dissertation_1544_Moll_Andr_2007.pdf
    // Many thanks to Andreas Moll and the BALLView developers for this

    Eigen::Vector3d v1, v2;

    // Calculate the vector between atom1 and atom2,
    // test if the vector has length larger than 0 and normalize it
    v1 = a - b;
    v2 = c - b;

    const double length1 = v1.norm();
    const double length2 = v2.norm();

    // test if the vector has length larger than 0 and normalize it
    if (IsNearZero(length1) || IsNearZero(length2)) {
      Fa = VZero;
      Fb = VZero;
      Fc = VZero;
      return 0.0;
    }

    // Calculate the normalized bond vectors
    const double inverse_length_v1 = 1.0 / length1;
    const double inverse_length_v2 = 1.0 / length2;
    v1 *= inverse_length_v1 ;
    v2 *= inverse_length_v2;

    // Calculate the cross product of v1 and v2, test if it has length unequal 0,
    // and normalize it.
    Eigen::Vector3d c1 = v1.cross(v2);
    const double length = c1.norm();
    if (IsNearZero(length)) {
      Fa = VZero;
      Fb = VZero;
      Fc = VZero;
      return 0.0;
    }

    c1 /= length;

    // Calculate the cos of theta and then theta
    double costheta = v1.dot(v2);
    double theta;
    if (costheta > 1.0) {
      theta = 0.0;
      costheta = 1.0;
    } else if (costheta < -1.0) {
      theta = 180.0;
      costheta = -1.0;
    } else {
      theta = RAD_TO_DEG * acos(costheta);
    }

    Eigen::Vector3d t1 = v1.cross(c1);
    t1.normalize();
    Eigen::Vector3d t2 = v2.cross(c1);
    t2.normalize();

    Fa = -t1 * inverse_length_v1;
    Fc =  t2 * inverse_length_v2;
    Fb = - (Fa + Fc);

    return theta;
  }
 

  /*! This method calculates the angle between two vectors
     
  \warning If length() of any of the two vectors is == 0.0,
  this method will divide by zero. If the product of the
  length() of the two vectors is very close to 0.0, but not ==
  0.0, this method may behave in unexpected ways and return
  almost random results; details may depend on your particular
  floating point implementation. The use of this method is
  therefore highly discouraged, unless you are certain that the
  length()es are in a reasonable range, away from 0.0 (Stefan
  Kebekus)

  \deprecated This method will probably replaced by a safer
  algorithm in the future.

  \todo Replace this method with a more fool-proof version.

  @returns the angle in degrees (0-360)
  */
  OBAPI double VectorAngle (const Eigen::Vector3d& ab, const Eigen::Vector3d& bc)
  {
    // length of the two bonds
    const double l_ab = ab.norm();
    const double l_bc = bc.norm();
 
    if (IsNearZero(l_ab) || IsNearZero(l_bc)) {
      return 0.0;
    }

    // Calculate the cross product of v1 and v2, test if it has length unequal 0
    const Eigen::Vector3d c1 = ab.cross(bc);
    if (IsNearZero(c1.norm())) {
      return 0.0;
    }

    // Calculate the cos of theta and then theta
    const double dp = ab.dot(bc) / (l_ab * l_bc);
    if (dp > 1.0) {
      return 0.0;
    } else if (dp < -1.0) {
      return 180.0;
    } else {
      return (RAD_TO_DEG * acos(dp));
    }
    
    return 0.0;
  }

  /*!  This function calculates the torsion angle of three vectors, represented
    by four points A--B--C--D, i.e. B and C are vertexes, but none of A--B,
    B--C, and C--D are colinear.  A "torsion angle" is the amount of "twist"
    or torsion needed around the B--C axis to bring A--B into the same plane
    as B--C--D.  The torsion is measured by "looking down" the vector B--C so
    that B is superimposed on C, then noting how far you'd have to rotate
    A--B to superimpose A over D.  Angles are + in theanticlockwise
    direction.  The operation is symmetrical in that if you reverse the image
    (look from C to B and rotate D over A), you get the same answer.
  */
  OBAPI double VectorTorsion(const Eigen::Vector3d& a, const Eigen::Vector3d& b,
      const Eigen::Vector3d& c, const Eigen::Vector3d& d)
  {
    // Bond vectors of the three atoms
    Eigen::Vector3d ab = b - a;
    Eigen::Vector3d bc = c - b;
    Eigen::Vector3d cd = d - c;

    // length of the three bonds
    const double l_ab = ab.norm();
    const double l_bc = bc.norm();
    const double l_cd = cd.norm();
    
    if (IsNearZero(l_ab) || IsNearZero(l_bc) || IsNearZero(l_cd) ) {
      return 0.0;
    }
 
    // normalize the bond vectors:
    ab *= (1.0 / l_ab);
    bc *= (1.0 / l_bc);
    cd *= (1.0 / l_cd);

    const Eigen::Vector3d ca = ab.cross(bc);
    const Eigen::Vector3d cb = bc.cross(cd);
    const Eigen::Vector3d cc = ca.cross(cb);
    const double d1 = cc.dot(bc);
    const double d2 = ca.dot(cb);
    const double tor = RAD_TO_DEG * atan2(d1, d2);
    
    return tor;  
  }


  /* Calculate the distance of point a to the plane determined by b,c,d */
  double Point2Plane(const Eigen::Vector3d& a, const Eigen::Vector3d& b, 
      const Eigen::Vector3d& c, const Eigen::Vector3d& d)
  {
    Eigen::Vector3d ab = a - b;
    Eigen::Vector3d bc = c - b;
    Eigen::Vector3d bd = d - b;
    Eigen::Vector3d normal = bc.cross(bd);
    return fabs( normal.dot( ab ) / normal.norm() );
  }
  
  /* Calculate the angle between point a and the plane determined by b,c,d */
  double Point2PlaneAngle(const Eigen::Vector3d& a, const Eigen::Vector3d& b, 
      const Eigen::Vector3d& c, const Eigen::Vector3d& d)
  {
    Eigen::Vector3d ac, bc, cd, normal;
    double angle;

    ac = a - c;
    bc = b - c;
    cd = c - d;
 
    normal = bc.cross(cd);
    angle = 90.0 - VectorAngle(normal, ac);

    return angle;
  }

  OBAPI double VectorTorsionDerivative(const Eigen::Vector3d &a, const Eigen::Vector3d &b, const Eigen::Vector3d &c,
      const Eigen::Vector3d &d, Eigen::Vector3d &Fa, Eigen::Vector3d &Fb, Eigen::Vector3d &Fc, Eigen::Vector3d &Fd) 
  {
    // This is adapted from http://scidok.sulb.uni-saarland.de/volltexte/2007/1325/pdf/Dissertation_1544_Moll_Andr_2007.pdf
    // Many thanks to Andreas Moll and the BALLView developers for this
    
    // angle between abc and bcd:
    double angle_abc, angle_bcd;
    
    // Bond vectors of the three atoms
    Eigen::Vector3d ab = b - a;
    Eigen::Vector3d bc = c - b;
    Eigen::Vector3d cd = d - c;
    
    // length of the three bonds
    const double l_ab = ab.norm();
    const double l_bc = bc.norm();
    const double l_cd = cd.norm();
    
    if (IsNearZero(l_ab) || IsNearZero(l_bc) || IsNearZero(l_cd) ) {
      Fa = VZero;
      Fb = VZero;
      Fc = VZero;
      Fd = VZero;
      return 0.0;
    }
    
    angle_abc = DEG_TO_RAD * OpenBabel::VectorAngle(ab, bc);
    angle_bcd = DEG_TO_RAD * OpenBabel::VectorAngle(bc, cd);
    
    // normalize the bond vectors:
    ab /= l_ab;
    bc /= l_bc;
    cd /= l_cd;

    const double sin_j = sin(angle_abc);
    const double sin_k = sin(angle_bcd);

    const double rsj = l_ab * sin_j;
    const double rsk = l_cd * sin_k;

    const double rs2j = 1. / (rsj * sin_j);
    const double rs2k = 1. / (rsk * sin_k);
  
    const double rrj = l_ab / l_bc;
    const double rrk = l_cd / l_bc;

    const double rrcj = rrj * (-cos(angle_abc));
    const double rrck = rrk * (-cos(angle_bcd));
    
    const Eigen::Vector3d ca = ab.cross(bc);
    const Eigen::Vector3d cb = bc.cross(cd);
    const Eigen::Vector3d cc = ca.cross(cb);
    double d1 = cc.dot(bc);
    double d2 = ca.dot(cb);
    double tor = RAD_TO_DEG * atan2(d1, d2);

    Fa = -ca * rs2j;
    Fd = cb * rs2k;

    Fb = Fa * (rrcj - 1.) - Fd * rrck;
    Fc = -(Fa + Fb + Fd);
    
    return tor;  
  }
 
  OBAPI double VectorOOP(const Eigen::Vector3d &a, const Eigen::Vector3d &b,
      const Eigen::Vector3d &c,const Eigen::Vector3d &d)
  {
    // This is adapted from http://scidok.sulb.uni-saarland.de/volltexte/2007/1325/pdf/Dissertation_1544_Moll_Andr_2007.pdf
    // Many thanks to Andreas Moll and the BALLView developers for this

    // calculate normalized bond vectors from central atom to outer atoms:
    Eigen::Vector3d ab = a - b;
    // store length of this bond:
    const double length_ab = ab.norm();
    if (IsNearZero(length_ab)) {
      return 0.0;
    }
    // store the normalized bond vector from central atom to outer atoms:
    // normalize the bond vector:
    ab /= length_ab;
		
    Eigen::Vector3d cb = c - b;
    const double length_cb = cb.norm();
    if (IsNearZero(length_cb)) {
      return 0.0;
    }
    cb /= length_cb;
	
    Eigen::Vector3d db = d - b;
    const double length_db = db.norm();
    if (IsNearZero(length_db)) {
      return 0.0;
    }
    db /= length_db;
	
    // the normal vectors of the three planes:
    const Eigen::Vector3d an = ab.cross(cb); 
    const Eigen::Vector3d bn = cb.cross(db); 
    const Eigen::Vector3d cn = db.cross(ab); 

    // Bond angle ji to jk
    const double cos_theta = ab.dot(cb);
    const double theta = acos(cos_theta);
    // If theta equals 180 degree or 0 degree
    if (IsNearZero(theta) || IsNearZero(fabs(theta - M_PI))) {
      return 0.0;
    }
				
    const double sin_theta = sin(theta);
    const double sin_dl = an.dot(db) / sin_theta;

    // the wilson angle:
    const double dl = asin(sin_dl);

    return RAD_TO_DEG * dl;
  }

  OBAPI double VectorOOPDerivative(const Eigen::Vector3d &a, const Eigen::Vector3d &b, const Eigen::Vector3d &c, 
      const Eigen::Vector3d &d, Eigen::Vector3d &Fa, Eigen::Vector3d &Fb, Eigen::Vector3d &Fc, Eigen::Vector3d &Fd)
  {
    // This is adapted from http://scidok.sulb.uni-saarland.de/volltexte/2007/1325/pdf/Dissertation_1544_Moll_Andr_2007.pdf
    // Many thanks to Andreas Moll and the BALLView developers for this

    // temp variables:
    double length;
    Eigen::Vector3d delta;

    // calculate normalized bond vectors from central atom to outer atoms:
    delta = a - b;
    length = delta.norm();
    if (IsNearZero(length)) {
      Fa = VZero;
      Fb = VZero;
      Fc = VZero;
      Fd = VZero;
      return 0.0;
    }
    // normalize the bond vector:
    delta /= length;
    // store the normalized bond vector from central atom to outer atoms:
    const Eigen::Vector3d ab = delta;
    // store length of this bond:
    const double length_ab = length;
		
    delta = c - b;
    length = delta.norm();
    if (IsNearZero(length)) {
      Fa = VZero;
      Fb = VZero;
      Fc = VZero;
      Fd = VZero;
      return 0.0;
    }
    // normalize the bond vector:
    delta /= length;
    // store the normalized bond vector from central atom to outer atoms:
    const Eigen::Vector3d bc = delta;
    // store length of this bond:
    const double length_bc = length;
	
    delta = d - b;
    length = delta.norm();
    if (IsNearZero(length)) {
      Fa = VZero;
      Fb = VZero;
      Fc = VZero;
      Fd = VZero;
      return 0.0;
    }
    // normalize the bond vector:
    delta /= length;
    // store the normalized bond vector from central atom to outer atoms:
    const Eigen::Vector3d bd = delta;
    // store length of this bond:
    const double length_bd = length;
	
    // the normal vectors of the three planes:
    const Eigen::Vector3d an = ab.cross(bc);
    const Eigen::Vector3d bn = bc.cross(bd);
    const Eigen::Vector3d cn = bd.cross(ab);

    // Bond angle ab to bc
    const double cos_theta = ab.dot(bc);
    const double theta = acos(cos_theta);
    // If theta equals 180 degree or 0 degree
    if (IsNearZero(theta) || IsNearZero(fabs(theta - M_PI))) {
      Fa = VZero;
      Fb = VZero;
      Fc = VZero;
      Fd = VZero;
      return 0.0;
    }
				
    const double sin_theta = sin(theta);
    const double sin_dl = an.dot(bd) / sin_theta;

    // the wilson angle:
    const double dl = asin(sin_dl);

    // In case: wilson angle equals 0 or 180 degree: do nothing
    if (IsNearZero(dl) || IsNearZero(fabs(dl - M_PI))) {
      Fa = VZero;
      Fb = VZero;
      Fc = VZero;
      Fd = VZero;
      return RAD_TO_DEG * dl;
    }
				
    const double cos_dl = cos(dl);

    // if wilson angle equal 90 degree: abort
    if (cos_dl < 0.0001) {
      Fa = VZero;
      Fb = VZero;
      Fc = VZero;
      Fd = VZero;
      return RAD_TO_DEG * dl;
    }

    Fd = (an / sin_theta - bd * sin_dl) / length_bd;
    Fa = ((bn + (((-ab + bc * cos_theta) * sin_dl) / sin_theta)) / length_ab) / sin_theta;
    Fc = ((cn + (((-bc + ab * cos_theta) * sin_dl) / sin_theta)) / length_bc) / sin_theta;
    Fb = -(Fa + Fc + Fd);
    
    return RAD_TO_DEG * dl;
  }
  
  OBAPI void SetTorsion(double *c, int ref[4], double setang, std::vector<int> atoms)
  {
    unsigned int cidx[4];
    cidx[0] = (ref[0] - 1) * 3;
    cidx[1] = (ref[1] - 1) * 3;
    cidx[2] = (ref[2] - 1) * 3;
    cidx[3] = (ref[3] - 1) * 3;
    const Eigen::Vector3d va( c[cidx[0]], c[cidx[0]+1], c[cidx[0]+2] );
    const Eigen::Vector3d vb( c[cidx[1]], c[cidx[1]+1], c[cidx[1]+2] );
    const Eigen::Vector3d vc( c[cidx[2]], c[cidx[2]+1], c[cidx[2]+2] );
    const Eigen::Vector3d vd( c[cidx[3]], c[cidx[3]+1], c[cidx[3]+2] );
    // calculate the rotation angle 
    const double ang = setang - DEG_TO_RAD * VectorTorsion(va, vb, vc, vd);
     
    // the angles are the same... 
    if (fabs(ang) < 1e-5)
      return;
    
    // setup the rotation matrix
    const Eigen::Vector3d bc = (vb - vc).normalized();
    Eigen::Quaternion<double> m;
    m = Eigen::AngleAxis<double>(-ang, bc);

    // apply the rotation
    int j;
    for (int i = 0; i < atoms.size(); ++i) {
      j = (atoms[i] - 1) * 3;
        
      Eigen::Vector3d vj( c[j], c[j+1], c[j+2] );
      vj -= vb; // translate so b is at origin
      vj = m.toRotationMatrix() * vj;
      vj += vb; // translate back 

      c[j]   = vj.x();
      c[j+1] = vj.y();
      c[j+2] = vj.z();
    }
  }

  OBAPI void SetTorsion(double *c, unsigned int ref[4], double setang, std::vector<int> atoms)
  {
    unsigned int cidx[4];
    cidx[0] = (ref[0] - 1) * 3;
    cidx[1] = (ref[1] - 1) * 3;
    cidx[2] = (ref[2] - 1) * 3;
    cidx[3] = (ref[3] - 1) * 3;
    const Eigen::Vector3d va( c[cidx[0]], c[cidx[0]+1], c[cidx[0]+2] );
    const Eigen::Vector3d vb( c[cidx[1]], c[cidx[1]+1], c[cidx[1]+2] );
    const Eigen::Vector3d vc( c[cidx[2]], c[cidx[2]+1], c[cidx[2]+2] );
    const Eigen::Vector3d vd( c[cidx[3]], c[cidx[3]+1], c[cidx[3]+2] );
    // calculate the rotation angle 
    const double ang = setang - DEG_TO_RAD * VectorTorsion(va, vb, vc, vd);
     
    // the angles are the same... 
    if (fabs(ang) < 1e-5)
      return;
    
    // setup the rotation matrix
    const Eigen::Vector3d bc = (vb - vc).normalized();
    Eigen::Quaternion<double> m;
    m = Eigen::AngleAxis<double>(-ang, bc);

    // apply the rotation
    int j;
    for (int i = 0; i < atoms.size(); ++i) {
      j = (atoms[i] - 1) * 3;
        
      Eigen::Vector3d vj( c[j], c[j+1], c[j+2] );
      vj -= vb; // translate so b is at origin
      vj = m.toRotationMatrix() * vj;
      vj += vb; // translate back 

      c[j]   = vj.x();
      c[j+1] = vj.y();
      c[j+2] = vj.z();
    }
  }


} // namespace OpenBabel

//! \file vector.cpp
//! \brief Handle 3D coordinates.

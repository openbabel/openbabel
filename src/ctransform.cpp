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

//THIS
#include "ctransform.h"

using namespace std;

namespace OpenBabel {

/*!
**\brief Constructor
*/
OBCoordTrans::OBCoordTrans()
  {
    _trans[0] = 0.0f; _trans[1] = 0.0f; _trans[2] = 0.0f;
    _euler[0] = 0.0f; _euler[1] = 0.0f; _euler[2] = 0.0f;
    _rmat[0] = 1.0f; _rmat[1] = 0.0f; _rmat[2] = 0.0f;
    _rmat[3] = 0.0f; _rmat[4] = 1.0f; _rmat[5] = 0.0f;
    _rmat[6] = 0.0f; _rmat[7] = 0.0f; _rmat[8] = 1.0f;
  }

/*!
**\brief Copy constructor
*/
OBCoordTrans::OBCoordTrans(const OBCoordTrans& cp)
  {
    *this = cp;
  }

/*!
**\brief Copy constructor
*/
OBCoordTrans::~OBCoordTrans()
  {
  }

/*!
**\brief Assignment operator
*/ 
OBCoordTrans& OBCoordTrans::operator=(const OBCoordTrans& cp)
  {
    _trans[0] = cp._trans[0];
    _trans[1] = cp._trans[1];
    _trans[2] = cp._trans[2];
    _euler[0] = cp._euler[0];
    _euler[1] = cp._euler[1];
    _euler[2] = cp._euler[2];
    _rmat[0] = cp._rmat[0];
    _rmat[1] = cp._rmat[1];
    _rmat[2] = cp._rmat[2];
    _rmat[3] = cp._rmat[3];
    _rmat[4] = cp._rmat[4];
    _rmat[5] = cp._rmat[5];
    _rmat[6] = cp._rmat[6];
    _rmat[7] = cp._rmat[7];
    _rmat[8] = cp._rmat[8];
    return *this;
  }

/*!
**\brief Clears the object as if it were just constructed
*/ 
void OBCoordTrans::Clear()
  {
    _trans[0] = 0.0f; _trans[1] = 0.0f; _trans[2] = 0.0f;
    _euler[0] = 0.0f; _euler[1] = 0.0f; _euler[2] = 0.0f;
    _rmat[0] = 1.0f; _rmat[1] = 0.0f; _rmat[2] = 0.0f;
    _rmat[3] = 0.0f; _rmat[4] = 1.0f; _rmat[5] = 0.0f;
    _rmat[6] = 0.0f; _rmat[7] = 0.0f; _rmat[8] = 1.0f;
  }

/*!
**\brief Writes the object to a binary character array
**\param ccc character array to write too (preallocated)
**\return The number of bytes written
*/
unsigned int OBCoordTrans::WriteBinary(char* ccc)
  {
    unsigned int idx=0;
    idx += OB_io_write_binary(&ccc[idx], (char*)&_trans[0], sizeof(float), 3);
    idx += OB_io_write_binary(&ccc[idx], (char*)&_euler[0], sizeof(float), 3);
    return idx;
  }

/*!
**\brief Reads the object from a binary character array
**\param ccc character array to read from (preallocated)
**\return the number of bytes read
*/
unsigned int OBCoordTrans::ReadBinary(char* ccc)
  {
    unsigned int idx=0;
    idx += OB_io_read_binary(&ccc[idx], (char*)&_trans[0], sizeof(float), 3);
    idx += OB_io_read_binary(&ccc[idx], (char*)&_euler[0], sizeof(float), 3);
    EulerToRmatrix(_euler,_rmat);
    return idx;
  }

/*!
**\brief Write the object to an output stream
**\param ostr The output stream
*/
void OBCoordTrans::WriteBinary(ostream& ostr)
  {
    OB_io_write_binary(ostr, (char*) &_trans[0], sizeof(float), 3);
    OB_io_write_binary(ostr, (char*) &_euler[0], sizeof(float), 3);
  }

/*!
**\brief Read the object from an input stream
**\param istr The input stream
*/
void OBCoordTrans::ReadBinary(istream& istr)
  {
    OB_io_read_binary(istr, (char*) &_trans[0], sizeof(float), 3);
    OB_io_read_binary(istr, (char*) &_euler[0], sizeof(float), 3);
    EulerToRmatrix(_euler,_rmat);
  }


/*!
**\brief Changes this transform to the orginal transform followed
**by an additional transform.
**\param ct The additional coordinate transformation.
**\return A transform equivilant to the original transform followed
**by the transform given in ct.  Note that the communative property
**does not apply so this is \b not the same as applying ct then the
**original transform.
*/
OBCoordTrans& OBCoordTrans::operator+=(const OBCoordTrans& ct)
  {
    float xyz[12];
    unsigned int i;
    for (i=0 ; i<12 ; i++) xyz[i] = 0.0f;
    xyz[3*1+0] = 1.0f;
    xyz[3*2+1] = 1.0f;
    xyz[3*3+2] = 1.0f;

    Transform(xyz,4);
    ct.Transform(xyz,4);
    Setup(xyz);

    return *this;
  }

/*!
**\brief Combine two transformations.
**\param ct2 Second transformation
**\return A transformation equivilant to applying (*this) and then ct2.
**\note Obviously the communative property does not apply hence,
**ct1+ct2 is \b not the same as ct2+ct1.
*/
OBCoordTrans OBCoordTrans::operator+(const OBCoordTrans& ct2) const
  {
    float xyz[12];
    unsigned int i;
    for (i=0 ; i<12 ; i++) xyz[i] = 0.0f;
    xyz[3*1+0] = 1.0f;
    xyz[3*2+1] = 1.0f;
    xyz[3*3+2] = 1.0f;

    OBCoordTrans ct;
    Transform(xyz,4);
    ct2.Transform(xyz,4);
    ct.Setup(xyz);
    return ct;
  }


/*!
**\brief Sets up a coordinate transformation from an arbitrary set of coordinates
**in the initial and final reference frame.
**\param init_xyz An array with the coordinates in the initial reference frame.
**\param finial_xyz An array with the coordinate in the final reference frame.
**\param N Number of coordinates.  Note that this procedure \b will deal with
**the case of N = 0,1 or 2.  In the case of N=0 the identity transform is returned.
**In the case of N=1 the appropriate translation, without and translation is
**returned.  In the case of N=2 there are multiple degenerate transformations,
**and a correct, but arbitrary, transformation is returned.  
**\note init_xyz and final_xyz must be identical sets of coordinates except for
**the frame of reference.
*/
void OBCoordTrans::Setup(float *init_xyz, float *final_xyz, unsigned int N)
  {
    Clear();

    float xyz1[12],xyz2[12];
    unsigned int i,j;

    //Get first coordinate
    if (N) {
        for (i=0 ; i<3 ; i++) {
            xyz1[i] = init_xyz[i];
            xyz2[i] = final_xyz[i];
          } 
      }
    else return;

    //Get second coordinate
    if (N>1) {
        float dist,maxdist;
        unsigned int id=1;
        maxdist = 0.0;
        for (i=1 ; i<N ; i++) {
            for (dist=0.0,j=0 ; j<3 ; j++) 
                dist += (xyz1[j]-init_xyz[3*i+j])*(xyz1[j]-init_xyz[3*i+j]);
            dist = sqrt(dist);
            if (dist > maxdist) {maxdist = dist; id = i;} 
          }
        for (i=0 ; i<3 ; i++) {
            xyz1[3+i] = init_xyz[3*id+i];
            xyz2[3+i] = final_xyz[3*id+i];
          }
      }
    else {
        float euler[3],trans[3];
        for (i=0 ; i<3 ; i++) {
            euler[i] = 0.0f;
            trans[i] = xyz2[i] - xyz1[i];
          }
        SetupEulerTranslation(euler,trans);
        return;
      }



    //Get third coordinate
    if (N>2) {
      float mag,maxcross;
      float xx[3],yy[3],cr[3];
      unsigned int ic=1;
      for (j=0 ; j<3 ; j++) xx[j] = xyz1[3+j]-xyz1[j];
      maxcross = 0.0f;
      for (i=1 ; i<N ; i++) {
           for (j=0 ; j<3 ; j++) yy[j] = init_xyz[3*i+j] - xyz1[j];
           cr[0] =  xx[1]*yy[2] - xx[2]*yy[1];
           cr[1] = -xx[0]*yy[2] + xx[2]*yy[0];
           cr[2] =  xx[0]*yy[1] - xx[1]*yy[0];
           mag = sqrt(cr[0]*cr[0] + cr[1]*cr[1] + cr[2]*cr[2]);
           if (mag > maxcross) {maxcross=mag; ic = i;}
         }
        for(i=0 ; i<3 ; i++) {
            xyz1[6+i] = init_xyz[3*ic+i];
            xyz2[6+i] = final_xyz[3*ic+i];
          }
       }
    else {//Deal with case of just two coordinates (just make up an arbitrary non-degenerate third one)
       float xx[3],yy[3];
       xx[0] = xyz1[3+0] - xyz1[0]; xx[1] = xyz1[3+1] - xyz1[1]; xx[2] = xyz1[3+2] - xyz1[2];
       yy[0] = xx[2]; yy[1] = xx[0]; yy[2] = xx[1];
       xyz1[6+0] = yy[0] + xyz1[0]; xyz1[6+1] = yy[1] + xyz1[1]; xyz1[6+2] = yy[2] + xyz1[2];
 
       xx[0] = xyz2[3+0] - xyz2[0]; xx[1] = xyz2[3+1] - xyz2[1]; xx[2] = xyz2[3+2] - xyz2[2];
       yy[0] = xx[2]; yy[1] = xx[0]; yy[2] = xx[1];
       xyz2[6+0] = yy[0] + xyz2[0]; xyz2[6+1] = yy[1] + xyz2[1]; xyz2[6+2] = yy[2] + xyz2[2];
      }
 
    //If we have gotten this far then we have a set of three non-colinear point in two different
    //reference frames (xyz1 and xyz2).  We are now going to convert these coordinates into a set
    //of 4 coordinates such that (c2-c1),(c3-c1) and (c4-c1) are unit vectors of an arbitrary 3rd
    //reference frame.  These coordinates can then used to create transformations from the initial
    //and finial reference frames to this third reference frame.

    //Get the 4 coordinates in the initial frame
        float mag,dot;
        float xx1[3],yy1[3],zz1[3];
        //Normalize x unit vector
        mag = 0.0f;
        for (i=0 ; i<3 ; i++) xx1[i] = xyz1[3+i] - xyz1[i];
        for (i=0 ; i<3 ; i++) mag += xx1[i]*xx1[i];
        mag = sqrt(mag);
        for (i=0 ; i<3 ; i++) {
            xx1[i] /=mag;
            xyz1[3+i] = xx1[i] + xyz1[i];
          }

        //Get the y vector
        dot = 0.0f;
        for (i=0 ; i<3 ; i++) yy1[i] = xyz1[6+i] - xyz1[i];
        for (i=0 ; i<3 ; i++) dot += xx1[i]*yy1[i]; 
        for (i=0 ; i<3 ; i++) yy1[i] -= xx1[i]*dot;

        //Normalize the y vector
        mag = 0.0f;
        for (i=0 ; i<3 ; i++) mag += yy1[i]*yy1[i];
        mag = sqrt(mag);
        for (i=0 ; i<3 ; i++) {
            yy1[i] /=mag;
            xyz1[6+i] = yy1[i] + xyz1[i];
          }

        //Get the z unit vector
        zz1[0] =  xx1[1]*yy1[2] - xx1[2]*yy1[1];
        zz1[1] = -xx1[0]*yy1[2] + xx1[2]*yy1[0];
        zz1[2] =  xx1[0]*yy1[1] - xx1[1]*yy1[0];
        for (i=0 ; i<3 ; i++) xyz1[9+i] = zz1[i] + xyz1[i];  
 
 
    //Get the 4 coordinates in the final reference frame
        float xx2[3],yy2[3],zz2[3];
        //Normalize x unit vector
        mag = 0.0f;
        for (i=0 ; i<3 ; i++) xx2[i] = xyz2[3+i] - xyz2[i];
        for (i=0 ; i<3 ; i++) mag += xx2[i]*xx2[i];
        mag = sqrt(mag);
        for (i=0 ; i<3 ; i++) {
            xx2[i] /=mag;
            xyz2[3+i] = xx2[i] + xyz2[i];
          }
 
        //Get the y vector
        dot = 0.0f;
        for (i=0 ; i<3 ; i++) yy2[i] = xyz2[6+i] - xyz2[i];
        for (i=0 ; i<3 ; i++) dot += xx2[i]*yy2[i];
        for (i=0 ; i<3 ; i++) yy2[i] -= xx2[i]*dot;
 
        //Normalize the y vector
        mag = 0.0f;
        for (i=0 ; i<3 ; i++) mag += yy2[i]*yy2[i];
        mag = sqrt(mag);
        for (i=0 ; i<3 ; i++) {
            yy2[i] /=mag;
            xyz2[6+i] = yy2[i] + xyz2[i];
          }
 
        //Get the z unit vector
        zz2[0] =  xx2[1]*yy2[2] - xx2[2]*yy2[1];
        zz2[1] = -xx2[0]*yy2[2] + xx2[2]*yy2[0];
        zz2[2] =  xx2[0]*yy2[1] - xx2[1]*yy2[0];
        for (i=0 ; i<3 ; i++) xyz2[9+i] = zz2[i] + xyz2[i];
 
    //Get transformations from third reference frame to initial and finial reference frames
    OBCoordTrans cti,ctf;
    cti.Setup(xyz1);
    ctf.Setup(xyz2); 

    //Invert the transformation from the third reference frame to the initial frame
    cti.Invert();

    //Set the transformation between the initial and finial reference frames.  This is
    //done by combining the transformation from the initial frame to the third frame with
    //the transformation from the third frame to the final frame.
    *this = cti + ctf;

    return;
  }

/*!
**\brief Get's an angle given it's sine and cosine.
**\param cs cosine of the angle
**\param sn sine of the angle
**\return Angle in radians
*/
double OBCoordTrans::Angle(double sn, double cs)
  {
    double angle=0.0;
    if (fabs(cs) < fabs(sn)) {
        angle = acos(cs);
        if      (sn < 0.0f && angle < PI) angle = 2.0*PI - angle;
        else if (sn > 0.0f && angle > PI) angle = 2.0*PI - angle;
      }
    else {
        angle = asin(sn);
        if (angle < PI) {
            if      (cs < 0.0f && angle < 0.5*PI) angle = PI - angle;
            else if (cs > 0.0f && angle > 0.5*PI) angle = PI - angle;
          }
        else {
            if (cs < 0.0f && angle > 1.5*PI) angle = 3.0*PI - angle;
            if (cs > 0.0f && angle < 1.5*PI) angle = 3.0*PI - angle;
          }
      }
    return angle;
  }

/*!
**\brief Inverts this objects transformation.  (i.e.,
**instead of converting from reference frame 1 to 2
**it now converts from 2 to 1.)
*/
void OBCoordTrans::Invert()
  {
    float xyz[12];
    unsigned int i;
    for (i=0 ; i<12 ; i++) xyz[i] = 0.0f;
    xyz[3*1+0] = 1.0f;
    xyz[3*2+1] = 1.0f;
    xyz[3*3+2] = 1.0f;

    //Apply reverse translation
    for (i=0 ; i<4 ; i++) {
        xyz[3*i+0] -= _trans[0];
        xyz[3*i+1] -= _trans[1];
        xyz[3*i+2] -= _trans[2];
      }

    //Apply reverse rotation
    ApplyEulerInvert(_euler,xyz,4);

    //Setup new transformation
    Setup(xyz);
  }


/*!
**\brief Returns a \b rotation and \b translation
**(applied in that order) coresponding to this objects
**transformation.
**\param euler A length 3 array that will be 
**returned with euler angles of the rotation.
**The angles are applied in the following order
**around the following axis.  euler[0] rotation
**about the z-axis, euler[1] rotation about the
**x-axis and euler[2] rotation about the z-axis.
**\param trans A length 3 array that will be returned
**with the x,y and z components of the translation.
*/
void OBCoordTrans::GetEulerTranslation(float *euler, float *trans) const
  {
    unsigned int i;
    for (i=0 ; i<3 ; i++) {
        euler[i] = _euler[i];
        trans[i] = _trans[i];
      }
  }

/*!
**\brief Returns a \b translation and \b rotation
**(applied in that order) coresponding to this objects
**transformation.
**\param trans A length 3 array that will be returned
**with the x,y and z components of the translation.
**\param euler A length 3 array that will be
**returned with euler angles of the rotation.
**The angles are applied in the following order
**around the following axis.  euler[0] rotation
**about the z-axis, euler[1] rotation about the
**x-axis and euler[2] rotation about the z-axis.
*/
void OBCoordTrans::GetTranslationEuler(float *trans, float *euler) const
  {
    //Get rotation
    euler[0] = _euler[0];
    euler[1] = _euler[1];
    euler[2] = _euler[2];

    //Get Translation
    trans[0] = 0.0f;
    trans[1] = 0.0f;
    trans[2] = 0.0f;
    Transform(trans,1);
    ApplyEulerInvert(euler,trans,1);
  }

/*!
**\brief Returns a \b rotation and \b translation
**(applied in that order) coresponding to this objects
**transformation.
**\param rmat A length 9 array that will be returned
**with the elements of a rotation matrix.  rmat[3*i+j]
**is the value of the element in the i'th row and j'th column.
**\param trans A length 3 array that will be returned
**with the x,y and z components of the translation.
*/
void OBCoordTrans::GetRmatrixTranslation(float *rmat, float *trans) const
  {
    unsigned int i;
    for (i=0 ; i<3 ; i++) trans[i] = _trans[i];
    for (i=0 ; i<9 ; i++) rmat[i] = _rmat[i];
  }

/*!
**\brief Returns a \b rotation and \b translation
**(applied in that order) coresponding to this objects
**transformation.
**\param rmatrix A Matrix3x3 that will be returned with
**the rotation.
**\param tvec A Vector that will be returned with the translation.
*/
void OBCoordTrans::GetRmatrixTranslation(Matrix3x3& rmatrix, Vector& tvec)
  {
    float rmat[9],trans[3];
    GetRmatrixTranslation(rmat,trans); 
    unsigned int irow,icolumn;
    for (irow=0 ; irow<3 ; irow++) for (icolumn=0 ; icolumn<3 ; icolumn++) {
        rmatrix.Set(irow,icolumn,rmat[3*irow+icolumn]);
      }
    tvec.Set(trans);
  }

/*!
**\brief Returns a \b translation and \b rotation
**(applied in that order) coresponding to this objects
**transformation.
**\param trans A length 3 array that will be returned
**with the x,y and z components of the translation.
**\param rmat A length 9 array that will be returned
**with the elements of a rotation matrix.  rmat[3*i+j]
**is the value of the element in the i'th row and j'th column.
*/
void OBCoordTrans::GetTranslationRmatrix(float *trans, float *rmat) const
  {
    float euler[3];
    GetTranslationEuler(trans,euler);
    EulerToRmatrix(euler,rmat);
  }

/*!
**\brief Returns a \b translation and \b rotation
**(applied in that order) coresponding to this objects
**transformation.
**\param tvec A Vector that will be returned with the translation.
**\param rmatrix A Matrix3x3 that will be returned with
**the rotation.
*/
void OBCoordTrans::GetTranslationRmatrix(Vector& tvec, Matrix3x3& rmatrix)
  {
    float rmat[9],trans[3];
    GetTranslationRmatrix(trans,rmat);
    unsigned int irow,icolumn;
    for (irow=0 ; irow<3 ; irow++) for (icolumn=0 ; icolumn<3 ; icolumn++) {
        rmatrix.Set(irow,icolumn,rmat[3*irow+icolumn]);
      }
    tvec.Set(trans);
  }

/*!
**\brief Sets up a rotation matrix from three euler angles
**\param euler A length 3 array with euler angles.  euler[0]
**is a rotation about the z-axis, euler[1] is a rotation about
**the x-axis and euler[2] is a rotation about the z-axis.  The
**Euler rotations are applied in the order listed.
**\b IMPORTANT: Angles are in \b radians.
**\param rmat A length 9 array representing the rotation matrix.
**rmat[3*i+j] is the value of the element in the i'th row and j'th column.
*/
void OBCoordTrans::EulerToRmatrix(float *euler, float *rmat) const
  {
    float xyz[9];
    unsigned int i,j;
    for (i=0 ; i<9 ; i++) xyz[i] = 0.0f;
    xyz[3*0+0] = 1.0f;
    xyz[3*1+1] = 1.0f;
    xyz[3*2+2] = 1.0f;
    ApplyEuler(euler,xyz,3);
    for (i=0 ; i<3 ; i++) for (j=0 ; j<3 ; j++) rmat[3*i+j] = xyz[3*j+i]; 
  }


/*!
**\brief Sets up this objects transformation based on a specified
**\b rotation and \b translation (applied in that order).
**\param euler A length 3 array with euler angles.  euler[0]
**is a rotation about the z-axis, euler[1] is a rotation about
**the x-axis and euler[2] is a rotation about the z-axis.  The
**Euler rotations are applied in the order listed.
**\b IMPORTANT: Angles are in \b radians.
**\param trans A length 3 array with the x,y and z translations.
**\note This member function is distinct from
**\link SetupTranslationEuler() \endlink.  The transformation setup
**in this function applies the rotation before the translation.
*/
void OBCoordTrans::SetupEulerTranslation(float *euler, float *trans)
  {
    //unsigned int i;
    //for (i=0 ; i<3 ; i++) {_euler[i] = euler[i]; _trans[i] = trans[i];}
    //SetRmatrix();

    float xyz[12];
    unsigned int i;
    for (i=0 ; i<12 ; i++) xyz[i] = 0.0f;
    xyz[3*1+0] = 1.0f;
    xyz[3*2+1] = 1.0f;
    xyz[3*3+2] = 1.0f;
 
    ApplyEuler(euler,xyz,4);
    ApplyTranslation(trans,xyz,4);
    Setup(xyz);
  }

/*!
**\brief Sets up this objects transformation based on a specified
**\b translation and \b rotation (applied in that order).
**\param trans A length 3 array with the x,y and z translations.
**\param euler length 3 array with euler angles.  euler[0]
**is a rotation about the z-axis, euler[1] is a rotation about
**the x-axis and euler[2] is a rotation about the z-axis.  The
**Euler rotations are applied in the order listed.
**\b IMPORTANT: Angles are in \b radians.
**\note This member function is distinct from
**\link SetupEulerTranslation() \endlink.  The transformation setup
**in this function applies the translation before the rotation.
*/
void OBCoordTrans::SetupTranslationEuler(float *trans, float *euler)
  {
    float xyz[12];
    unsigned int i;
    for (i=0 ; i<12 ; i++) xyz[i] = 0.0f;
    xyz[3*1+0] = 1.0f;
    xyz[3*2+1] = 1.0f;
    xyz[3*3+2] = 1.0f;

    ApplyTranslation(trans,xyz,4);
    ApplyEuler(euler,xyz,4);
    Setup(xyz);
  }

/*!
**\brief Sets up this objects transformation based on a specified
**\b rotation and \b translation (applied in that order).
**\param rmat A length 9 array representing the rotation matrix. 
**rmat[3*i+j] is the value of the element in the i'th row and j'th column.
**\param trans A length 3 array with the x,y and z translations.
**\note This member function is distinct from
**\link SetupTranslationRmatrix() \endlink.  This transformation setup
**in this function applies the rotation before the translation.
*/
void OBCoordTrans::SetupRmatrixTranslation(float *rmat, float *trans)
  {
    float xyz[12];
    unsigned int i;
    for (i=0 ; i<12 ; i++) xyz[i] = 0.0f;
    xyz[3*1+0] = 1.0f;
    xyz[3*2+1] = 1.0f;
    xyz[3*3+2] = 1.0f;

    ApplyRmatrix(rmat,xyz,4);
    ApplyTranslation(trans,xyz,4);
    Setup(xyz);
  }

/*!
**\brief Sets up this objects transformation based on a specified
**\b rotation and \b translation (applied in that order).
**\param rmatrix A Matrix3x3 rotation matrix.
**\param tvec A Vector holding the translation
**\note This member function is distinct from
**\link SetupTranslationRmatrix() \endlink.  This transformation setup
**in this function applies the rotation before the translation.
*/
void OBCoordTrans::SetupRmatrixTranslation(Matrix3x3& rmatrix, Vector& tvec)
  {
    float rmat[9],trans[3];
    unsigned int irow,icolumn;
    for (irow=0 ; irow<3 ; irow++) for (icolumn=0 ; icolumn<3 ; icolumn++) {
        rmat[3*irow+icolumn] = rmatrix.Get(irow,icolumn);
      }
    tvec.Get(trans);
    SetupRmatrixTranslation(rmat,trans);
  }

/*!
**\brief Sets up this objects transformation based on a specified
**\b translation and \b rotation (applied in that order).
**\param rmat A length 9 array representing the rotation matrix.
**rmat[3*i+j] is the value of the element in the i'th row and j'th column.
**\param trans A length 3 array with the x,y and z translations.
**\note This member function is distinct from
**\link SetupRmatrixTranslation() \endlink.  This transformation setup
**in this function applies the translation before the rotation.
*/
void OBCoordTrans::SetupTranslationRmatrix(float *trans, float *rmat)
  {
    float xyz[12];
    unsigned int i;
    for (i=0 ; i<12 ; i++) xyz[i] = 0.0f;
    xyz[3*1+0] = 1.0f;
    xyz[3*2+1] = 1.0f;
    xyz[3*3+2] = 1.0f;
 
    ApplyTranslation(trans,xyz,4);
    ApplyRmatrix(rmat,xyz,4);
    Setup(xyz);
  }

/*!
**\brief Sets up this objects transformation based on a specified
**\b translation and \b rotation (applied in that order).
**\param rmatrix A Matrix3x3 rotation matrix.
**\param tvec A Vector holding the translation
**\note This member function is distinct from
**\link SetupRmatrixTranslation() \endlink.  This transformation setup
**in this function applies the translation before the rotation.
*/
void OBCoordTrans::SetupTranslationRmatrix(Vector& tvec, Matrix3x3& rmatrix)
  {
    float rmat[9],trans[3];
    unsigned int irow,icolumn;
    for (irow=0 ; irow<3 ; irow++) for (icolumn=0 ; icolumn<3 ; icolumn++) {
        rmat[3*irow+icolumn] = rmatrix.Get(irow,icolumn);
      }
    tvec.Get(trans);
    SetupTranslationRmatrix(trans,rmat);
  }

/*!
**\brief Applies this objects transformation to a set of coordinates
**\param xyz Array of coordinates to be transformed
**\param N Number of coordinates in xyz array
*/
void OBCoordTrans::Transform(float *xyz, unsigned int N) const
  {
    unsigned int i;
    float x,y,z;
    for (i=0 ; i<N ; i++) {
        x = xyz[3*i+0];
        y = xyz[3*i+1];
        z = xyz[3*i+2];
        xyz[3*i+0] = _rmat[0]*x + _rmat[1]*y + _rmat[2]*z + _trans[0];
        xyz[3*i+1] = _rmat[3]*x + _rmat[4]*y + _rmat[5]*z + _trans[1];
        xyz[3*i+2] = _rmat[6]*x + _rmat[7]*y + _rmat[8]*z + _trans[2];
      } 
  }

/*!
**\brief Applies a translation (no rotation) to a set of coordates.
**param trans A length 3 array holding the x,y and z translation
**\param xyz Array of coordinates to be transformed
**\param N Number of coordinates in xyz array
*/
void OBCoordTrans::ApplyTranslation(float *trans, float *xyz, unsigned int N) const
  {
    unsigned int i;
    for (i=0 ; i<N ; i++) {
        xyz[3*i+0] += trans[0];
        xyz[3*i+1] += trans[1];
        xyz[3*i+2] += trans[2];
      }
  }

/*!
**\brief Applies a rotation (no translation) to a set of coordinates.
**\param rmat A length 9 array representing the rotation matrix. 
**rmat[3*i+j] is the value of the element in the i'th row and j'th column.
**\param xyz Array of coordinates to be transformed
**\param N Number of coordinates in xyz array
*/
void OBCoordTrans::ApplyRmatrix(float *rmat, float *xyz, unsigned int N) const
  {
    unsigned int i;
    float x,y,z;
    for (i=0 ; i<N ; i++) {
        x = rmat[0]*xyz[3*i+0] + rmat[1]*xyz[3*i+1] + rmat[2]*xyz[3*i+2];
        y = rmat[3]*xyz[3*i+0] + rmat[4]*xyz[3*i+1] + rmat[5]*xyz[3*i+2];
        z = rmat[6]*xyz[3*i+0] + rmat[7]*xyz[3*i+1] + rmat[8]*xyz[3*i+2];
        xyz[3*i+0] = x;
        xyz[3*i+1] = y;
        xyz[3*i+2] = z;
      }
  }

/*!
**\brief Applies an euler rotation (no translation) in reverse
**to a set of coordinates.
**\param euler A length 3 array with the euler angles.  euler[0]
**is a rotation about the z-axis, euler[1] is a rotation about
**the x-axis and euler[2] is a rotation about the z-axis.  In
**general the rotations are applied in the order listed, however,
**in this procedure they are applied in reverse order.  
**\b IMPORTANT: Angles are in \b radians.
**\param xyz Array of coordinates to be transformed
**\param N Number of coordinates in xyz array
*/
void OBCoordTrans::ApplyEulerInvert(float *euler, float *xyz, unsigned int N) const
  {
    float cs0 = cos(euler[0]);
    float cs1 = cos(euler[1]);
    float cs2 = cos(euler[2]);
    float sn0 = sin(euler[0]);
    float sn1 = sin(euler[1]);
    float sn2 = sin(euler[2]);

    //Reverse third Euler angle
    unsigned int i;
    float xx,yy;
    for (i=0 ; i<N ; i++) {
        xx = xyz[3*i+0];
        yy = xyz[3*i+1];
        xyz[3*i+0] = cs2*xx - sn2*yy;
        xyz[3*i+1] = cs2*yy + sn2*xx;
      }

    //Reverse second Euler angle
    float zz;
    for (i=0 ; i<N ; i++) {
        yy = xyz[3*i+1];
        zz = xyz[3*i+2];
        xyz[3*i+1] = cs1*yy - sn1*zz;
        xyz[3*i+2] = cs1*zz + sn1*yy;
      }

    //Reverse first Euler angle
    for (i=0 ; i<N ; i++) {
        xx = xyz[3*i+0];
        yy = xyz[3*i+1];
        xyz[3*i+0] = cs0*xx - sn0*yy;
        xyz[3*i+1] = cs0*yy + sn0*xx;
      }
  }

/*!
**\brief Applies an euler rotation (no translation)
**to a set of coordinates.
**\param euler A length 3 array with the euler angles.  euler[0]
**is a rotation about the z-axis, euler[1] is a rotation about
**the x-axis and euler[2] is a rotation about the z-axis.  The
**rotations are applied in the order listed.
**\b IMPORTANT: Angles are in \b radians.
**\param xyz Array of coordinates to be transformed
**\param N Number of coordinates in xyz array
**/
void OBCoordTrans::ApplyEuler(float *euler, float *xyz, unsigned int N) const
  {
    float cs0 = cos(euler[0]);
    float cs1 = cos(euler[1]);
    float cs2 = cos(euler[2]);
    float sn0 = sin(euler[0]);
    float sn1 = sin(euler[1]);
    float sn2 = sin(euler[2]);

    unsigned int i;
    
    //Apply first Euler Angle (rotation about z-axis)
    float xx,yy;
    for (i=0 ; i<N ; i++) {
        xx = xyz[3*i+0];
        yy = xyz[3*i+1];
        xyz[3*i+0] = cs0*xx + sn0*yy;
        xyz[3*i+1] = cs0*yy - sn0*xx;
      }

    //Apply second Euler Angle (rotation about the x-axis)
    float zz;
    for (i=0 ; i<N ; i++) {
        yy = xyz[3*i+1];
        zz = xyz[3*i+2];
        xyz[3*i+1] = cs1*yy + sn1*zz;
        xyz[3*i+2] = cs1*zz - sn1*yy;
      }

    //Apply third Euler Angle (rotation about the z-axis)
    for (i=0 ; i<N ; i++) {
        xx = xyz[3*i+0];
        yy = xyz[3*i+1];
        xyz[3*i+0] = cs2*xx + sn2*yy;
        xyz[3*i+1] = cs2*yy - sn2*xx;
      }

  }

/*!
**\brief Core function to setup the coordinate transformation
**\param in_xyz A length 12 array containing 4 coordinates
**(0,0,0), (1,0,0), (0,1,0) and (0,0,1) from the initial 
**reference frame transformed into the final reference frame.
*/
bool OBCoordTrans::Setup(float *in_xyz)
  {
    //Copy coordinate array
    double xyz[12];
    double *y=&xyz[3*2];
    double *z=&xyz[3*3];
    unsigned int i;
    for (i=0 ; i<12 ; i++) xyz[i] = in_xyz[i];

    //Set translation
    _trans[0] = xyz[0];
    _trans[1] = xyz[1];
    _trans[2] = xyz[2];

    //DEBUG
    //char buffer[1000];
    //cout << "DEBUG : Initial coordinates" << endl;
    //for (i=0 ; i<4 ; i++) {sprintf(buffer,"DEBUG : (%10.6f,%10.6f,%10.6f)",xyz[3*i+0],xyz[3*i+1],xyz[3*i+2]); cout << buffer << endl;} cout << endl;

    //Undo translation
    for (i=0 ; i<4 ; i++) {
        xyz[3*i+0] -= _trans[0]; 
        xyz[3*i+1] -= _trans[1]; 
        xyz[3*i+2] -= _trans[2]; 
      }

    //DEBUG
    //cout << "DEBUG : Undid translation coordinates" << endl;
    //for (i=0 ; i<4 ; i++) {sprintf(buffer,"DEBUG : (%10.6f,%10.6f,%10.6f)",xyz[3*i+0],xyz[3*i+1],xyz[3*i+2]); cout << buffer << endl;} cout << endl;

    //Find the angle of the rotated z unit vector with
    //the y axis IN THE XY PLANE. (i.e., the third Euler angle)
    double sn,cs;
    double mag;
    mag = sqrt(z[0]*z[0] + z[1]*z[1]);
    if (mag > 0.000001) {
        cs = z[1]/mag;
        sn = z[0]/mag;
        _euler[2] = (float) Angle(sn,cs);
      }
    else {
        cs = 1.0f;
        sn = 0.0f;
        _euler[2] = 0.0f;
      }

    //Undo the rotation from the third Euler angle
    double xx,yy;
    for (i=0 ; i<4 ; i++) {
        xx = xyz[3*i+0];
        yy = xyz[3*i+1]; 
        xyz[3*i+0] = cs*xx - sn*yy;
        xyz[3*i+1] = cs*yy + sn*xx;
      }

    //DEBUG
    //cout << "DEBUG : Undid third Euler rotation : " << _euler[2] << endl;
    //for (i=0 ; i<4 ; i++) {sprintf(buffer,"DEBUG : (%10.6f,%10.6f,%10.6f)",xyz[3*i+0],xyz[3*i+1],xyz[3*i+2]); cout << buffer << endl;} cout << endl;

    //Find the angle of the rotated z unit vector with
    //the z axis.  (i.e., the second Euler angle)
    cs = z[2];
    sn = z[1];
    _euler[1] = (float) Angle(sn,cs);

    //Undo the rotation from the second Euler angle
    double zz;
    for (i=0 ; i<4 ; i++) {
        yy = xyz[3*i+1];
        zz = xyz[3*i+2];
        xyz[3*i+1] = cs*yy - sn*zz;
        xyz[3*i+2] = cs*zz + sn*yy;
      } 

    //DEBUG
    //cout << "DEBUG : Undid third Euler rotation : " << _euler[2] << endl;
    //for (i=0 ; i<4 ; i++) {sprintf(buffer,"DEBUG : (%10.6f,%10.6f,%10.6f)",xyz[3*i+0],xyz[3*i+1],xyz[3*i+2]); cout << buffer << endl;} cout << endl;

    //Find the angle of the rotated y unit vector with
    //the y axis (i.e., the first Euler angle)
    cs = y[1];
    sn = y[0];
    _euler[0] = (float) Angle(sn,cs);

    //Find the rotation matrix coresponding to the euler angles
    EulerToRmatrix(_euler,_rmat);

    //Undo the rotation from the first Euler angle
    //This is unnecessary except as a check to make
    //sure we were given valid input.
    for (i=0 ; i<4 ; i++) {
        xx = xyz[3*i+0];
        yy = xyz[3*i+1];
        xyz[3*i+0] = cs*xx - sn*yy;
        xyz[3*i+1] = cs*yy + sn*xx;
      }

    //DEBUG
    //cout << "DEBUG : Undid first Euler rotation : " << _euler[2] << endl;
    //for (i=0 ; i<4 ; i++) {sprintf(buffer,"DEBUG : (%10.6f,%10.6f,%10.6f)",xyz[3*i+0],xyz[3*i+1],xyz[3*i+2]); cout << buffer << endl;} cout << endl;

    bool error=false;
    float tol=0.0001f;
    if (fabs(xyz[0]) > tol) error = true;
    if (fabs(xyz[1]) > tol) error = true;
    if (fabs(xyz[2]) > tol) error = true;
    if (fabs(xyz[3*1+0] - 1.0) > tol) error = true;
    if (fabs(xyz[3*2+1] - 1.0) > tol) error = true;
    if (fabs(xyz[3*3+2] - 1.0) > tol) error = true;
    if (error) {cerr << "WARNING! OBCoordTrans::Setup(float*) probable invalid input" << endl; return false;}

    return true;
  }


}//End OpenEye namespace














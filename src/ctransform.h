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

#ifndef OB_COORDTRANSFORM_INCLUDED
#define OB_COORDTRANSFORM_INCLUDED

#include <math.h>

#include "binary_io.h"
#include "Vector.h"

namespace OpenBabel {

/*!
**\brief An object for storing, manipulating and applying coordinate transformations.
**\par Setup
**The most basic way to setup a transformation is to use the member function Setup(float* xyz),
**which requires you to supply a specific set of 4 coordinates from the initial reference frame
**in the final reference frame.  The function Setup(float *init_xyz, float* final_xyz, unsigned int N)
**will set up a transformation given an arbitrary set of coordinates in the final and initial
**reference frame.  A rotation/translation or a translation/rotation can also be supplied,
**see the various member functions begining with Setup.
**\param Applying the transformation
**The Transform(float *xyz, unsigned int N) member function applies the objects transform
**to a set of coordinates.
**\par Utilities
**The Invert() member function will reverses the current transformation such that it changes
**what was originally the final reference frame into what was originally the initial reference frame.
**A rotation/translation or translation/rotation cooresponding to the objects transformation
**can be extracted, see the various member functions begining with Get.  The operator+(const OBCoordTrans& ct2)
**and operator+=(const OBCoordTrans& ct) are also defined for combining transformation.
**\note \b IMPORTANT : This object can only handle transformations that can be described
**with Euler angles (i.e., No inversions).  Several member function will take rotation
**matricies as input, however, it is assumed that these matricies only apply simple
**rotations that can be described with Euler angles.
*/
class OBCoordTrans {
    protected:
        float _trans[3];
        float _euler[3];
        float _rmat[9];

        double Angle(double sn, double cs);
        void ApplyTranslation(float *trans, float *xyz, unsigned int N) const;
        void ApplyEulerInvert(float *euler, float *xyz, unsigned int N) const;
        void ApplyEuler(float *euler, float *xyz, unsigned int N) const;
        void ApplyRmatrix(float *rmat, float *xyz, unsigned int N) const;
        void EulerToRmatrix(float *euler, float *rmat) const;
    public:
        //Constructor, destructor and copy constructor
        OBCoordTrans();
        OBCoordTrans(const OBCoordTrans& cp);
        ~OBCoordTrans();

        //Operator overloads
        OBCoordTrans& operator=(const OBCoordTrans& cp);
        OBCoordTrans& operator+=(const OBCoordTrans& ct);
        OBCoordTrans operator+(const OBCoordTrans& ct2) const;

        //Clear function
        void Clear();

        //Read/Write to binary array/stream
        unsigned int WriteBinary(char* ccc);
        unsigned int ReadBinary(char* ccc);
        void WriteBinary(std::ostream& ostr);
        void ReadBinary(std::istream& istr);

        //Setup functions
        bool Setup(float *xyz);
        void Setup(float *init_xyz, float *final_xyz, unsigned int N);
        void SetupEulerTranslation(float *euler, float *trans);
        void SetupTranslationEuler(float *trans, float *euler);
        void SetupRmatrixTranslation(float *rmat, float *trans);
        void SetupTranslationRmatrix(float *trans, float *rmat);
        void SetupRmatrixTranslation(Matrix3x3& rmat, Vector& tvec);
        void SetupTranslationRmatrix(Vector& tvec, Matrix3x3& rmat);

        //Retrieve transform
        void GetEulerTranslation(float *euler, float *trans) const;
        void GetTranslationEuler(float *trans, float *euler) const;
        void GetRmatrixTranslation(float *rmat, float *trans) const;
        void GetTranslationRmatrix(float *trans, float *rmat) const;
        void GetRmatrixTranslation(Matrix3x3& rmat, Vector& tvec);
        void GetTranslationRmatrix(Vector& tvec, Matrix3x3& rmat);
        
        //Transformation functions
        void Transform(float *xyz, unsigned int N) const;

        //Invert
        void Invert();
  };

/*!
**\fn OBCoordTrans::_trans()
**\brief Stores the translation part of the transformation.  This
**translation is applied after the rotation.
*/

/*!
**\fn OBCoordTrans::_euler()
**\brief Stores the rotation part of the transformation in euler angles.
**This rotation is applied before the translation.  The elements of
**the array are : euler[0], a rotation about the z-axis; euler[1],
**a rotation about the x-axis; euler[2], a rotation about the z-axis.
**The rotations are applied in the order given.
*/

/*!
**\fn OBCoordTrans::_rmat()
**\brief Stores the rotation part of the transformation as a
**rotation matrix.  This rotation is applied before the translation.
**_rmat[3*i+j] hold the matrix element in the i'th row and j'th column.
*/


} //End OpenEye namespace

#endif


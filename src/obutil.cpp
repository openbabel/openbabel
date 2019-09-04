/**********************************************************************
obutil.cpp - Various utility methods.

Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2006 by Geoffrey R. Hutchison

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
#include <openbabel/math/matrix3x3.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/obiter.h>
#include <openbabel/obutil.h>
#include <openbabel/internalcoord.h>

#include <cstring>

#ifdef HAVE_CONIO_H
#include <conio.h>
#endif

using namespace std;
namespace OpenBabel
{

  /*! \class OBStopwatch obutil.h <openbabel/obutil.h>
    \brief Stopwatch class used for timing length of execution

    The OBStopwatch class makes timing the execution of blocks of
    code to microsecond accuracy very simple. The class effectively
    has two functions, Start() and Elapsed(). The usage of the
    OBStopwatch class is demonstrated by the following code:
    \code
    OBStopwatch sw;
    sw.Start();
    //insert code here
    cout << "Elapsed time = " << sw.Elapsed() << endl;
    \endcode
  */

  //! Deprecated: use the OBMessageHandler class instead
  //! \deprecated Throw an error through the OpenBabel::OBMessageHandler class
  void ThrowError(char *str)
  {
    obErrorLog.ThrowError("", str, obInfo);
  }

  //! Deprecated: use the OBMessageHandler class instead
  //! \deprecated Throw an error through the OpenBabel::OBMessageHandler class
  void ThrowError(std::string &str)
  {
    obErrorLog.ThrowError("", str, obInfo);
  }

  // returns True if a < b, False otherwise.
  bool OBCompareInt(const int &a,const int &b)
  {
    return(a<b);
  }

  // Comparison function (for sorting unsigned ints) returns a < b
  bool OBCompareUnsigned(const unsigned int &a,const unsigned int &b)
  {
    return(a<b);
  }

  //! Comparison for doubles: returns fabs(a - b) < epsilon
  bool IsNear(const double &a, const double &b, const double epsilon)
  {
    return (fabs(a - b) < epsilon);
  }

  //! Comparison for doubles: returns fabs(a) < epsilon
  bool IsNearZero(const double &a, const double epsilon)
  {
    return (fabs(a) < epsilon);
  }

  //! Comparison for nan (not a number)
  bool IsNan(const double &a)
  {
    return ((a) != (a));
  }

  //! Tests whether its argument can be squared without triggering an overflow or
  //! underflow.
  bool CanBeSquared(const double &a)
  {
    if( a == 0 ) return true;
    const double max_squarable_double = 1e150;
    const double min_squarable_double = 1e-150;
    double abs_a = fabs(a);
    return(abs_a < max_squarable_double && abs_a > min_squarable_double);
  }

  //! Utility function: replace the last extension in string &src with new extension char *ext.
  string NewExtension(string &src,char *ext)
  {
    string::size_type pos = (unsigned int)src.find_last_of(".");
    string dst;

    if (pos != string::npos)
      dst = src.substr(0,pos+1);
    else
      {
        dst = src;
        dst += ".";
      }

    dst += ext;
    return(dst);
  }

  //! \return the geometric centroid to an array of coordinates in double* format
  //!  and center the coordinates to the origin. Operates on the first "size"
  //!  coordinates in the array.
  vector3 center_coords(double *c, unsigned int size)
  {
    if (size == 0)
      {
        return(VZero);
      }
		unsigned int i;
    double x=0.0, y=0.0, z=0.0;
    for (i = 0;i < size;++i)
      {
        x += c[i*3];
        y += c[i*3+1];
        z += c[i*3+2];
      }
    x /= (double) size;
    y /= (double) size;
    z /= (double) size;
    for (i = 0;i < size;++i)
      {
        c[i*3]   -= x;
        c[i*3+1] -= y;
        c[i*3+2] -= z;
      }
    vector3 v(x,y,z);
    return(v);
  }

  //! Rotates the coordinate set *c by the transformation matrix m[3][3]
  //!  Operates on the first "size" coordinates in the array.
  void rotate_coords(double *c,double m[3][3],unsigned int size)
  {
    double x,y,z;
    for (unsigned int i = 0;i < size;++i)
      {
        x = c[i*3]*m[0][0] + c[i*3+1]*m[0][1] + c[i*3+2]*m[0][2];
        y = c[i*3]*m[1][0] + c[i*3+1]*m[1][1] + c[i*3+2]*m[1][2];
        z = c[i*3]*m[2][0] + c[i*3+1]*m[2][1] + c[i*3+2]*m[2][2];
        c[i*3] = x;
        c[i*3+1] = y;
        c[i*3+2] = z;
      }
  }

  //! Calculate the RMS deviation between the first N coordinates of *r and *f
  double calc_rms(double *r,double *f, unsigned int N)
  {
    if (N == 0)
      return 0.0; // no RMS deviation between two empty sets

    double d2=0.0;
    for (unsigned int i = 0;i < N;++i)
      {
        d2 += SQUARE(r[i*3] - f[i*3]) +
          SQUARE(r[i*3+1] - f[i*3+1]) +
          SQUARE(r[i*3+2] - f[i*3+2]);
      }

    d2 /= (double) N;
    return(sqrt(d2));
  }

  //! Rotate the coordinates of 'atoms'
  //! such that tor == ang - atoms in 'tor' should be ordered such
  //! that the 3rd atom is the pivot around which atoms rotate
  void SetRotorToAngle(double *c,vector<int> &tor,double ang,vector<int> &atoms)
  {
    double v1x,v1y,v1z,v2x,v2y,v2z,v3x,v3y,v3z;
    double c1x,c1y,c1z,c2x,c2y,c2z,c3x,c3y,c3z;
    double c1mag,c2mag,radang,costheta,m[9];
    double x,y,z,mag,rotang,sn,cs,t,tx,ty,tz;

    //
    //calculate the torsion angle
    //
    v1x = c[tor[0]]   - c[tor[1]];
    v2x = c[tor[1]]   - c[tor[2]];
    v1y = c[tor[0]+1] - c[tor[1]+1];
    v2y = c[tor[1]+1] - c[tor[2]+1];
    v1z = c[tor[0]+2] - c[tor[1]+2];
    v2z = c[tor[1]+2] - c[tor[2]+2];
    v3x = c[tor[2]]   - c[tor[3]];
    v3y = c[tor[2]+1] - c[tor[3]+1];
    v3z = c[tor[2]+2] - c[tor[3]+2];

    c1x = v1y*v2z - v1z*v2y;
    c2x = v2y*v3z - v2z*v3y;
    c1y = -v1x*v2z + v1z*v2x;
    c2y = -v2x*v3z + v2z*v3x;
    c1z = v1x*v2y - v1y*v2x;
    c2z = v2x*v3y - v2y*v3x;
    c3x = c1y*c2z - c1z*c2y;
    c3y = -c1x*c2z + c1z*c2x;
    c3z = c1x*c2y - c1y*c2x;

    c1mag = SQUARE(c1x)+SQUARE(c1y)+SQUARE(c1z);
    c2mag = SQUARE(c2x)+SQUARE(c2y)+SQUARE(c2z);
    if (c1mag*c2mag < 0.01)
      costheta = 1.0; //avoid div by zero error
    else
      costheta = (c1x*c2x + c1y*c2y + c1z*c2z)/(sqrt(c1mag*c2mag));

    if (costheta < -0.999999)
      costheta = -0.999999;
    if (costheta >  0.999999)
      costheta =  0.999999;

    if ((v2x*c3x + v2y*c3y + v2z*c3z) > 0.0)
      radang = -acos(costheta);
    else
      radang = acos(costheta);

    //
    // now we have the torsion angle (radang) - set up the rot matrix
    //

    //find the difference between current and requested
    rotang = ang - radang;

    sn = sin(rotang);
    cs = cos(rotang);
    t = 1 - cs;
    //normalize the rotation vector
    mag = sqrt(SQUARE(v2x)+SQUARE(v2y)+SQUARE(v2z));
    x = v2x/mag;
    y = v2y/mag;
    z = v2z/mag;

    //set up the rotation matrix
    m[0]= t*x*x + cs;
    m[1] = t*x*y + sn*z;
    m[2] = t*x*z - sn*y;
    m[3] = t*x*y - sn*z;
    m[4] = t*y*y + cs;
    m[5] = t*y*z + sn*x;
    m[6] = t*x*z + sn*y;
    m[7] = t*y*z - sn*x;
    m[8] = t*z*z + cs;

    //
    //now the matrix is set - time to rotate the atoms
    //
    tx = c[tor[1]];
    ty = c[tor[1]+1];
    tz = c[tor[1]+2];
    vector<int>::iterator i;
    int j;
    for (i = atoms.begin();i != atoms.end();++i)
      {
        j = *i;
        c[j] -= tx;
        c[j+1] -= ty;
        c[j+2]-= tz;
        x = c[j]*m[0] + c[j+1]*m[1] + c[j+2]*m[2];
        y = c[j]*m[3] + c[j+1]*m[4] + c[j+2]*m[5];
        z = c[j]*m[6] + c[j+1]*m[7] + c[j+2]*m[8];
        c[j] = x;
        c[j+1] = y;
        c[j+2] = z;
        c[j] += tx;
        c[j+1] += ty;
        c[j+2] += tz;
      }
  }

  //! Safely open the supplied filename and return an ifstream, throwing an error
  //! to the default OBMessageHandler error log if it fails.
  bool SafeOpen(std::ifstream &fs, const char *filename)
  {
#ifdef WIN32
    string s(filename);
    if (s.find(".bin") != string::npos)
      fs.open(filename,ios::binary);
    else
#endif

      fs.open(filename);

    if (!fs)
      {
        string error = "Unable to open file \'";
        error += filename;
        error += "\' in read mode";
        obErrorLog.ThrowError(__FUNCTION__, error, obError);
        return(false);
      }

    return(true);
  }


  //! Safely open the supplied filename and return an ofstream, throwing an error
  //! to the default OBMessageHandler error log if it fails.
  bool SafeOpen(std::ofstream &fs, const char *filename)
  {
#ifdef WIN32
    string s(filename);
    if (s.find(".bin") != string::npos)
      fs.open(filename,ios::binary);
    else
#endif

      fs.open(filename);

    if (!fs)
      {
        string error = "Unable to open file \'";
        error += filename;
        error += "\' in write mode";
        obErrorLog.ThrowError(__FUNCTION__, error, obError);
        return(false);
      }

    return(true);
  }

  //! Safely open the supplied filename and return an ifstream, throwing an error
  //! to the default OBMessageHandler error log if it fails.
  bool SafeOpen(std::ifstream &fs, const string &filename)
  {
    return(SafeOpen(fs, filename.c_str()));
  }

  //! Safely open the supplied filename and return an ofstream, throwing an error
  //! to the default OBMessageHandler error log if it fails.
  bool SafeOpen(std::ofstream &fs, const string &filename)
  {
    return(SafeOpen(fs, filename.c_str()));
  }

  //! Shift the supplied string to uppercase
  void ToUpper(std::string &s)
  {
    if (s.empty())
      return;
    unsigned int i;
    for (i = 0;i < s.size();++i)
      if (isalpha(s[i]) && !isdigit(s[i]))
        s[i] = toupper(s[i]);
  }

  //! Shift the supplied char* to uppercase
  void ToUpper(char *cptr)
  {
    char *c;
    for (c = cptr;*c != '\0';++c)
      if (isalpha(*c) && !isdigit(*c))
        *c = toupper(*c);
  }

  //! Shift the supplied string to lowercase
  void ToLower(std::string &s)
  {
    if (s.empty())
      return;
    unsigned int i;
    for (i = 0;i < s.size();++i)
      if (isalpha(s[i]) && !isdigit(s[i]))
        s[i] = tolower(s[i]);
  }

  //! Shift the supplied char* to lowercase
  void ToLower(char *cptr)
  {
    char *c;
    for (c = cptr;*c != '\0';++c)
      if (isalpha(*c) && !isdigit(*c))
        *c = tolower(*c);
  }

  //! Shift the supplied string: lowercase to upper, and upper to lower
  //! \param s - The string to switch case
  //! \param start - The position to start inverting case
  void InvertCase(std::string &s, unsigned int start)
  {
    if (start >= s.size())
      return;
    unsigned int i;
    for (i = start; i < s.size();++i)
      if (isalpha(s[i]) && !isdigit(s[i])) {
        if (isupper(s[i])) s[i] = tolower(s[i]);
        else s[i] = toupper(s[i]);
      }
  }

  //! Shift the supplied char*: lowercase to upper, and upper to lower
  void InvertCase(char *cptr)
  {
    char *c;
    for (c = cptr;*c != '\0';++c)
      if (isalpha(*c) && !isdigit(*c)) {
        if (isupper(*c)) *c = tolower(*c);
        else *c = toupper(*c);
      }
  }

  //! "Clean" the supplied atom type, shifting the first character to uppercase,
  //! the second character (if it's a letter) to lowercase, and terminating with a NULL
  //! to strip off any trailing characters
  void CleanAtomType(char *id)
  {
    id[0] = toupper(id[0]);
    if (isalpha(id[1]) == 0)
      id[1] = '\0';
    else
      {
        id[1] = tolower(id[1]);
        id[2] = '\0';
      }
  }

  //! Transform the supplied vector<OBInternalCoord*> into cartesian and update
  //! the OBMol accordingly. The size of supplied internal coordinate vector
  //! has to be the same as the number of atoms in molecule (+ NULL in the
  //! beginning).
  //! Implements <a href="http://qsar.sourceforge.net/dicts/blue-obelisk/index.xhtml#zmatrixCoordinatesIntoCartesianCoordinates">blue-obelisk:zmatrixCoordinatesIntoCartesianCoordinates</a>
  void InternalToCartesian(std::vector<OBInternalCoord*> &vic,OBMol &mol)
  {
    vector3 n,nn,v1,v2,v3,avec,bvec,cvec;
    double dst = 0.0, ang = 0.0, tor = 0.0;
    OBAtom *atom;
    vector<OBAtom*>::iterator i;
    unsigned int index;

    if (vic.empty())
      return;

    if (vic[0] != NULL) {
      std::vector<OBInternalCoord*>::iterator it = vic.begin();
      vic.insert(it, static_cast<OBInternalCoord*>(NULL));
    }

    if (vic.size() != mol.NumAtoms() + 1) {
      string error = "Number of internal coordinates is not the same as";
      error += " the number of atoms in molecule";
      obErrorLog.ThrowError(__FUNCTION__, error, obError);
      return;
    }

    obErrorLog.ThrowError(__FUNCTION__,
                          "Ran OpenBabel::InternalToCartesian", obAuditMsg);

    for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
      {
        index = atom->GetIdx();

        // make sure we always have valid pointers
        if (index >= vic.size() || !vic[index])
          return;

        if (vic[index]->_a) // make sure we have a valid ptr
          {
            avec = vic[index]->_a->GetVector();
            dst = vic[index]->_dst;
          }
        else
          {
            // atom 1
            atom->SetVector(0.0, 0.0, 0.0);
            continue;
          }

        if (vic[index]->_b)
          {
            bvec = vic[index]->_b->GetVector();
            ang = vic[index]->_ang * DEG_TO_RAD;
          }
        else
          {
            // atom 2
            atom->SetVector(dst, 0.0, 0.0);
            continue;
          }

        if (vic[index]->_c)
          {
            cvec = vic[index]->_c->GetVector();
            tor = vic[index]->_tor * DEG_TO_RAD;
          }
        else
          {
            // atom 3
            cvec = VY;
            tor = 90. * DEG_TO_RAD;
          }

        v1 = avec - bvec;
        v2 = avec - cvec;
        n = cross(v1,v2);
        nn = cross(v1,n);
        n.normalize();
        nn.normalize();

        n  *= -sin(tor);
        nn *= cos(tor);
        v3 = n + nn;
        v3.normalize();
        v3 *= dst * sin(ang);
        v1.normalize();
        v1 *= dst * cos(ang);
        v2 = avec + v3 - v1;

        atom->SetVector(v2);
      }

    // Delete dummy atoms
    vector<OBAtom*> for_deletion;
    FOR_ATOMS_OF_MOL(a, mol)
      if (a->GetAtomicNum() == 0)
        for_deletion.push_back(&(*a));
    for(vector<OBAtom*>::iterator a_it=for_deletion.begin(); a_it!=for_deletion.end(); ++a_it)
      mol.DeleteAtom(*a_it);

  }

  //! Use the supplied OBMol and its Cartesian coordinates to generate
  //! a set of internal (z-matrix) coordinates as supplied in the
  //! vector<OBInternalCoord*> argument.
  //! Implements <a href="http://qsar.sourceforge.net/dicts/blue-obelisk/index.xhtml#cartesianCoordinatesIntoZmatrixCoordinates">blue-obelisk:cartesianCoordinatesIntoZmatrixCoordinates</a>.
  void CartesianToInternal(std::vector<OBInternalCoord*> &vic,OBMol &mol)
  {
    double r,sum;
    OBAtom *atom,*nbr,*ref;
    vector<OBAtom*>::iterator i,j,m;

    obErrorLog.ThrowError(__FUNCTION__,
                          "Ran OpenBabel::CartesianToInternal", obAuditMsg);

    //set reference atoms
    for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
      {
        if      (atom->GetIdx() == 1)
          continue;
        else if (atom->GetIdx() == 2)
          {
            vic[atom->GetIdx()]->_a = mol.GetAtom(1);
            continue;
          }
        else if (atom->GetIdx() == 3)
          {
            if( (atom->GetVector()-mol.GetAtom(2)->GetVector()).length_2()
                <(atom->GetVector()-mol.GetAtom(1)->GetVector()).length_2())
              {
                vic[atom->GetIdx()]->_a = mol.GetAtom(2);
                vic[atom->GetIdx()]->_b = mol.GetAtom(1);
              }
            else
              {
                vic[atom->GetIdx()]->_a = mol.GetAtom(1);
                vic[atom->GetIdx()]->_b = mol.GetAtom(2);
              }
            continue;
          }
        sum=1.0E10;
        ref = mol.GetAtom(1);
        for(nbr = mol.BeginAtom(j);nbr && (i != j);nbr = mol.NextAtom(j))
          {
            r = (atom->GetVector()-nbr->GetVector()).length_2();
            if((r < sum) && (vic[nbr->GetIdx()]->_a != nbr) &&
               (vic[nbr->GetIdx()]->_b != nbr))
              {
                sum = r;
                ref = nbr;
              }
          }

        vic[atom->GetIdx()]->_a = ref;
        if (ref->GetIdx() >= 3)
          {
            vic[atom->GetIdx()]->_b = vic[ref->GetIdx()]->_a;
            vic[atom->GetIdx()]->_c = vic[ref->GetIdx()]->_b;
          }
        else
          {
            if(ref->GetIdx()== 1)
              {
                vic[atom->GetIdx()]->_b = mol.GetAtom(2);
                vic[atom->GetIdx()]->_c = mol.GetAtom(3);
              }
            else
              {//ref->GetIdx()== 2
                vic[atom->GetIdx()]->_b = mol.GetAtom(1);
                vic[atom->GetIdx()]->_c = mol.GetAtom(3);
              }
          }
      }

    //fill in geometries
    unsigned int k;
    vector3 v1,v2;
    OBAtom *a,*b,*c;
    for (k = 2;k <= mol.NumAtoms();++k)
      {
        atom = mol.GetAtom(k);
        a = vic[k]->_a;
        b = vic[k]->_b;
        c = vic[k]->_c;
        v1 = atom->GetVector() - a->GetVector();
        vic[k]->_dst = v1.length();
        if (k == 2)
          continue;

        v2 = b->GetVector()    - a->GetVector();
        vic[k]->_ang = vectorAngle(v1,v2);
        if (k == 3)
          continue;

        vic[k]->_tor = CalcTorsionAngle(atom->GetVector(),
                                        a->GetVector(),
                                        b->GetVector(),
                                        c->GetVector());
      }

    //check for linear geometries and try to correct if possible
    bool done;
    double ang;
    for (k = 2;k <= mol.NumAtoms();++k)
      {
        ang = fabs(vic[k]->_ang);
        if (ang > 5.0 && ang < 175.0)
          continue;
        atom = mol.GetAtom(k);
        done = false;
        for (a = mol.BeginAtom(i);a && a->GetIdx() < k && !done;a = mol.NextAtom(i))
          for (b=mol.BeginAtom(j);b && b->GetIdx()<a->GetIdx() && !done;b = mol.NextAtom(j))
            {
              v1 = atom->GetVector() - a->GetVector();
              v2 = b->GetVector() - a->GetVector();
              ang = fabs(vectorAngle(v1,v2));
              if (ang < 5.0 || ang > 175.0)
                continue;

              // Also check length considerations -- don't bother if the length > 10.0 Angstroms
              if (v1.length_2() > 99.999)
                continue;

              for (c = mol.BeginAtom(m);c && c->GetIdx() < atom->GetIdx();c = mol.NextAtom(m))
                if (c != atom && c != a && c != b)
                  break;
              if (!c)
                continue;

              vic[k]->_a = a;
              vic[k]->_b = b;
              vic[k]->_c = c;
              vic[k]->_dst = v1.length();
              vic[k]->_ang = vectorAngle(v1,v2);
              vic[k]->_tor = CalcTorsionAngle(atom->GetVector(),
                                              a->GetVector(),
                                              b->GetVector(),
                                              c->GetVector());
              if (!isfinite(vic[k]->_tor))
                vic[k]->_tor = 180.0;
              done = true;
            }
      }
  }

  void qtrfit (double *r,double *f,int size, double u[3][3])
  {
    int i;
    double xxyx, xxyy, xxyz;
    double xyyx, xyyy, xyyz;
    double xzyx, xzyy, xzyz;
    double d[4],q[4];
    double c[16],v[16];
    double rx,ry,rz,fx,fy,fz;

    /* generate the upper triangle of the quadratic form matrix */

    xxyx = 0.0;
    xxyy = 0.0;
    xxyz = 0.0;
    xyyx = 0.0;
    xyyy = 0.0;
    xyyz = 0.0;
    xzyx = 0.0;
    xzyy = 0.0;
    xzyz = 0.0;

    for (i = 0; i < size; ++i)
      {
        rx = r[i*3];
        ry = r[i*3+1];
        rz = r[i*3+2];
        fx = f[i*3];
        fy = f[i*3+1];
        fz = f[i*3+2];

        xxyx += fx * rx;
        xxyy += fx * ry;
        xxyz += fx * rz;
        xyyx += fy * rx;
        xyyy += fy * ry;
        xyyz += fy * rz;
        xzyx += fz * rx;
        xzyy += fz * ry;
        xzyz += fz * rz;
      }

    c[4*0+0] = xxyx + xyyy + xzyz;

    c[4*0+1] = xzyy - xyyz;
    c[4*1+1] = xxyx - xyyy - xzyz;

    c[4*0+2] = xxyz - xzyx;
    c[4*1+2] = xxyy + xyyx;
    c[4*2+2] = xyyy - xzyz - xxyx;

    c[4*0+3] = xyyx - xxyy;
    c[4*1+3] = xzyx + xxyz;
    c[4*2+3] = xyyz + xzyy;
    c[4*3+3] = xzyz - xxyx - xyyy;

    /* diagonalize c */

    matrix3x3::jacobi(4, c, d, v);

    /* extract the desired quaternion */

    q[0] = v[4*0+3];
    q[1] = v[4*1+3];
    q[2] = v[4*2+3];
    q[3] = v[4*3+3];

    /* generate the rotation matrix */

    u[0][0] = q[0]*q[0] + q[1]*q[1] - q[2]*q[2] - q[3]*q[3];
    u[1][0] = 2.0 * (q[1] * q[2] - q[0] * q[3]);
    u[2][0] = 2.0 * (q[1] * q[3] + q[0] * q[2]);

    u[0][1] = 2.0 * (q[2] * q[1] + q[0] * q[3]);
    u[1][1] = q[0]*q[0] - q[1]*q[1] + q[2]*q[2] - q[3]*q[3];
    u[2][1] = 2.0 * (q[2] * q[3] - q[0] * q[1]);

    u[0][2] = 2.0 * (q[3] * q[1] - q[0] * q[2]);
    u[1][2] = 2.0 * (q[3] * q[2] + q[0] * q[1]);
    u[2][2] = q[0]*q[0] - q[1]*q[1] - q[2]*q[2] + q[3]*q[3];
  }



  static double Roots[4];

#define ApproxZero 1E-7
#define IsZero(x)  ((double)fabs(x)<ApproxZero)
#ifndef PI
#define PI         3.14159265358979323846226433
#endif
#define OneThird      (1.0/3.0)
#define FourThirdsPI  (4.0*PI/3.0)
#define TwoThirdsPI   (2.0*PI/3.0)

  /*FUNCTION */
  /* receives: the co-efficients for a general
   *           equation of degree one.
   *           Ax + B = 0 !!
   */
  int SolveLinear(double A,double B)
  {
    if( IsZero(A) )
      return( 0 );
    Roots[0] = -B/A;
    return( 1 );
  }

  /*FUNCTION */
  /* receives: the co-efficients for a general
   *           linear equation of degree two.
   *           Ax^2 + Bx + C = 0 !!
   */
  int SolveQuadratic(double A,double B,double C)
  {
    double Descr, Temp, TwoA;

    if( IsZero(A) )
      return( SolveLinear(B,C) );

    TwoA = A+A;
    Temp = TwoA*C;
    Descr = B*B - (Temp+Temp);
    if( Descr<0.0 )
      return( 0 );

    if( Descr>0.0 )
      {
        Descr = sqrt(Descr);
        /* W. Press, B. Flannery, S. Teukolsky and W. Vetterling,
         * "Quadratic and Cubic Equations", Numerical Recipes in C,
         * Chapter 5, pp. 156-157, 1989.
         */
        Temp = (B<0.0)? -0.5*(B-Descr) : -0.5*(B+Descr);
        Roots[0] = Temp/A;
        Roots[1] = C/Temp;

        return( 2 );
      }
    Roots[0] = -B/TwoA;
    return( 1 );
  }

  /*FUNCTION */
  /* task: to return the cube root of the
   *       given value taking into account
   *       that it may be negative.
   */
  double CubeRoot(double X)
  {
    if( X>=0.0 )
      {
        return pow( X, OneThird );
      }
    else
      return -pow( -X, OneThird );
  }

  int SolveCubic(double A,double B,double C,double D)
  {
    double TwoA, ThreeA, BOver3A;
    double Temp, POver3, QOver2;
    double Desc, Rho, Psi;


    if( IsZero(A) )
      {
        return( SolveQuadratic(B,C,D) );
      }

    TwoA = A+A;
    ThreeA = TwoA+A;
    BOver3A = B/ThreeA;
    QOver2 = ((TwoA*BOver3A*BOver3A-C)*BOver3A+D)/TwoA;
    POver3 = (C-B*BOver3A)/ThreeA;


    Rho = POver3*POver3*POver3;
    Desc = QOver2*QOver2 + Rho;

    if( Desc<=0.0 )
      {
        Rho = sqrt( -Rho );
        Psi = OneThird*acos(-QOver2/Rho);
        Temp = CubeRoot( Rho );
        Temp = Temp+Temp;

        Roots[0] = Temp*cos( Psi )-BOver3A;
        Roots[1] = Temp*cos( Psi+TwoThirdsPI )-BOver3A;
        Roots[2] = Temp*cos( Psi+FourThirdsPI )-BOver3A;
        return( 3 );
      }

    if( Desc> 0.0 )
      {
        Temp = CubeRoot( -QOver2 );
        Roots[0] = Temp+Temp-BOver3A;
        Roots[1] = -Temp-BOver3A;
        return( 2 );
      }

    Desc = sqrt( Desc );
    Roots[0] = CubeRoot(Desc-QOver2)-CubeRoot(Desc+QOver2) - BOver3A;

    return( 1 );
  }


#define MAX_SWEEPS 50

  void ob_make_rmat(double a[3][3],double rmat[9])
  {
  	/*
  	onorm, dnorm - hold the sum of diagonals and off diagonals to check Jacobi completion
  	d[3] - holds the diagonals of the input vector, which transofrm to become the Eigenvalues
  	r1, r2 - hold 1st two Eigenvectors
  	v1,v2,v3 - hold orthogonal unit vectors derived from Eigenvectors
  	
  	The junction performs a Jacobi Eigenvalue/vector determination 
  	(https://en.wikipedia.org/wiki/Jacobi_eigenvalue_algorithm) on the supplied
  	Inertial Tensor in a, and returns a unit transform matrix rmat as a row matrix.
  	To work, a must be diagonally symmetric (i.e a[i][j] = a[j][i])
  	v starts out holding the unit matrix (i.e. no transform in co-ordinate frame), 
  	and undergoes the same rotations as applied to the Inertial Tensor in the Jacobi
  	process to arrive at the new co-ordinate frame.
  	Finally, the eigenvalues are sorted in order that the largest principal moment aligns to the 
  	new x-axis
  	*/
    double onorm, dnorm; 
    double b, dma, q, t, c, s,d[3];
    double atemp, vtemp, dtemp,v[3][3];
    double r1[3],r2[3],v1[3],v2[3],v3[3];
    int i, j, k, l;

    memset((char*)d,'\0',sizeof(double)*3);

    for (j = 0; j < 3; ++j)
      {
        for (i = 0; i < 3; ++i)
          v[i][j] = 0.0;

        v[j][j] = 1.0;
        d[j] = a[j][j];
      }

    for (l = 1; l <= MAX_SWEEPS; ++l)
      {
        dnorm = 0.0;
        onorm = 0.0;
        for (j = 0; j < 3; ++j)
          {
            dnorm = dnorm + (double)fabs(d[j]);
            for (i = 0; i <= j - 1; ++i)
              {
                onorm = onorm + (double)fabs(a[i][j]);
              }
          }

        if((onorm/dnorm) <= 1.0e-12)
        /* Completion achieved (i.e. off-diagonals are all 0.0 within error)*/
          goto Exit_now;
        for (j = 1; j < 3; ++j)
          {
            for (i = 0; i <= j - 1; ++i)
              {
                b = a[i][j];
                if(fabs(b) > 0.0)
                  {
                    dma = d[j] - d[i];
                    if((fabs(dma) + fabs(b)) <=  fabs(dma))
                      t = b / dma;
                    else
                      {
                        q = 0.5 * dma / b;
                        t = 1.0/((double)fabs(q) + (double)sqrt(1.0+q*q));
                        if(q < 0.0)
                          t = -t;
                      }
                    c = 1.0/(double)sqrt(t * t + 1.0);
                    s = t * c;
                    a[i][j] = 0.0;
                    /* Perform a Jacobi rotation on the supplied matrix*/
                    for (k = 0; k <= i-1; ++k)
                      {
                        atemp = c * a[k][i] - s * a[k][j];
                        a[k][j] = s * a[k][i] + c * a[k][j];
                        a[k][i] = atemp;
                      }
                    for (k = i+1; k <= j-1; ++k)
                      {
                        atemp = c * a[i][k] - s * a[k][j];
                        a[k][j] = s * a[i][k] + c * a[k][j];
                        a[i][k] = atemp;
                      }
                    for (k = j+1; k < 3; ++k)
                      {
                        atemp = c * a[i][k] - s * a[j][k];
                        a[j][k] = s * a[i][k] + c * a[j][k];
                        a[i][k] = atemp;
                      }
                      /* Rotate the reference frame */
                    for (k = 0; k < 3; ++k)
                      {
                        vtemp = c * v[k][i] - s * v[k][j];
                        v[k][j] = s * v[k][i] + c * v[k][j];
                        v[k][i] = vtemp;
                      }
                    dtemp = c*c*d[i] + s*s*d[j] - 2.0*c*s*b;
                    d[j] = s*s*d[i] + c*c*d[j] +  2.0*c*s*b;
                    d[i] = dtemp;
                  }  /* end if */
              } /* end for i */
          } /* end for j */
      } /* end for l */

  Exit_now:

    /* max_sweeps = l;*/

/* Now sort the eigenvalues and eigenvectors*/
    for (j = 0; j < 3-1; ++j)
      {
        k = j;
        dtemp = d[k];
        for (i = j+1; i < 3; ++i)
          if(d[i] < dtemp)
            {
              k = i;
              dtemp = d[k];
            }

        if(k > j)
          {
            d[k] = d[j];
            d[j] = dtemp;
            for (i = 0; i < 3 ; ++i)
              {
                dtemp = v[i][k];
                v[i][k] = v[i][j];
                v[i][j] = dtemp;
              }
          }
      }

	/* Transfer the 1st two eigenvectors into r1 and r2*/
    r1[0] = v[0][0];
    r1[1] = v[1][0];
    r1[2] = v[2][0];
    r2[0] = v[0][1];
    r2[1] = v[1][1];
    r2[2] = v[2][1];

	/* Generate the 3rd unit vector for the new co-ordinate frame by cross product of r1 and r2*/
    v3[0] =  r1[1]*r2[2] - r1[2]*r2[1];
    v3[1] = -r1[0]*r2[2] + r1[2]*r2[0];
    v3[2] =  r1[0]*r2[1] - r1[1]*r2[0];
	/* Ensure it is normalised |v3|=1 */
    s = (double)sqrt(v3[0]*v3[0] + v3[1]*v3[1] + v3[2]*v3[2]);
    v3[0] /= s;
    v3[1] /= s;
    v3[2] /= s;

	/* Generate the 2nd unit vector for the new co-ordinate frame by cross product of v3 and r1*/
    v2[0] =  v3[1]*r1[2] - v3[2]*r1[1];
    v2[1] = -v3[0]*r1[2] + v3[2]*r1[0];
    v2[2] =  v3[0]*r1[1] - v3[1]*r1[0];
    /* Ensure it is normalised |v2|=1 */
    s = (double)sqrt(v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2]);
    v2[0] /= s;
    v2[1] /= s;
    v2[2] /= s;

	/* Generate the 1st unit vector for the new co-ordinate frame by cross product of v2 and v3*/
    v1[0] =  v2[1]*v3[2] - v2[2]*v3[1];
    v1[1] = -v2[0]*v3[2] + v2[2]*v3[0];
    v1[2] =  v2[0]*v3[1] - v2[1]*v3[0];
     /* Ensure it is normalised |v1|=1 */
    s = (double)sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]);
    v1[0] /= s;
    v1[1] /= s;
    v1[2] /= s;

	/* Transfer to the row matrix form for the result*/
    rmat[0] = v1[0];
    rmat[1] = v1[1];
    rmat[2] = v1[2];
    rmat[3] = v2[0];
    rmat[4] = v2[1];
    rmat[5] = v2[2];
    rmat[6] = v3[0];
    rmat[7] = v3[1];
    rmat[8] = v3[2];
  }

  static int get_roots_3_3(double mat[3][3], double roots[3])
  {
    double rmat[9];

    ob_make_rmat(mat,rmat);

    mat[0][0]=rmat[0];
    mat[0][1]=rmat[3];
    mat[0][2]=rmat[6];
    mat[1][0]=rmat[1];
    mat[1][1]=rmat[4];
    mat[1][2]=rmat[7];
    mat[2][0]=rmat[2];
    mat[2][1]=rmat[5];
    mat[2][2]=rmat[8];

    roots[0]=(double)Roots[0];
    roots[1]=(double)Roots[1];
    roots[2]=(double)Roots[2];

    return 1;
  }

  double superimpose(double *r,double *f,int size)
  {
    int i,j;
    double x,y,z,d2;
    double mat[3][3],rmat[3][3],mat2[3][3],roots[3];

    /* make inertial cross tensor */
    for(i=0;i<3;++i)
      for(j=0;j<3;++j)
        mat[i][j]=0.0;

    for(i=0;i < size;++i)
      {
        mat[0][0]+=r[3*i]  *f[3*i];
        mat[1][0]+=r[3*i+1]*f[3*i];
        mat[2][0]+=r[3*i+2]*f[3*i];
        mat[0][1]+=r[3*i]  *f[3*i+1];
        mat[1][1]+=r[3*i+1]*f[3*i+1];
        mat[2][1]+=r[3*i+2]*f[3*i+1];
        mat[0][2]+=r[3*i]  *f[3*i+2];
        mat[1][2]+=r[3*i+1]*f[3*i+2];
        mat[2][2]+=r[3*i+2]*f[3*i+2];
      }

    d2=mat[0][0]*(mat[1][1]*mat[2][2]-mat[1][2]*mat[2][1])
      -mat[0][1]*(mat[1][0]*mat[2][2]-mat[1][2]*mat[2][0])
      +mat[0][2]*(mat[1][0]*mat[2][1]-mat[1][1]*mat[2][0]);


    /* square matrix= ((mat transpose) * mat) */
    for(i=0;i<3;++i)
      for(j=0;j<3;++j)
        {
          x=mat[0][i]*mat[0][j]+mat[1][i]*mat[1][j]+mat[2][i]*mat[2][j];
          mat2[i][j]=mat[i][j];
          rmat[i][j]=x;
        }
    get_roots_3_3(rmat,roots);

    roots[0]=(roots[0]<0.0001) ? 0.0: (roots[0]);
    roots[1]=(roots[1]<0.0001) ? 0.0: (roots[1]);
    roots[2]=(roots[2]<0.0001) ? 0.0: (roots[2]);

    /* make sqrt of rmat, store in mat*/

    roots[0]=roots[0]<0.0001? 0.0: 1.0/(double)sqrt(roots[0]);
    roots[1]=roots[1]<0.0001? 0.0: 1.0/(double)sqrt(roots[1]);
    roots[2]=roots[2]<0.0001? 0.0: 1.0/(double)sqrt(roots[2]);

    if(d2<0.0)
      {
        if( (roots[0]>=roots[1]) && (roots[0]>=roots[2]) )
          roots[0]*=-1.0;
        if( (roots[1]>roots[0]) && (roots[1]>=roots[2]) )
          roots[1]*=-1.0;
        if( (roots[2]>roots[1]) && (roots[2]>roots[0]) )
          roots[2]*=-1.0;
      }

    for(i=0;i<3;++i)
      for(j=0;j<3;++j)
        mat[i][j]=roots[0]*rmat[i][0]*rmat[j][0]+
          roots[1]*rmat[i][1]*rmat[j][1]+
          roots[2]*rmat[i][2]*rmat[j][2];

    /* and multiply into original inertial cross matrix, mat2 */
    for(i=0;i<3;++i)
      for(j=0;j<3;++j)
        rmat[i][j]=mat[0][j]*mat2[i][0]+
          mat[1][j]*mat2[i][1]+
          mat[2][j]*mat2[i][2];

    /* rotate all coordinates */
    d2 = 0.0;
    for(i=0;i<size;++i)
      {
        x=f[3*i]*rmat[0][0]+f[3*i+1]*rmat[0][1]+f[3*i+2]*rmat[0][2];
        y=f[3*i]*rmat[1][0]+f[3*i+1]*rmat[1][1]+f[3*i+2]*rmat[1][2];
        z=f[3*i]*rmat[2][0]+f[3*i+1]*rmat[2][1]+f[3*i+2]*rmat[2][2];
        f[3*i  ]=x;
        f[3*i+1]=y;
        f[3*i+2]=z;

        x = r[i*3]   - f[i*3];
        y = r[i*3+1] - f[i*3+1];
        z = r[i*3+2] - f[i*3+2];
        d2 += x*x+y*y+z*z;
      }

    d2 /= (double) size;

    return((double)sqrt(d2));
  }

  void get_rmat(double *rvec,double *r,double *f,int size)
  {
    int i,j;
    double x,d2;
    double mat[3][3],rmat[3][3],mat2[3][3],roots[3];

    /* make inertial cross tensor */
    for(i=0;i<3;++i)
      for(j=0;j<3;++j)
        mat[i][j]=0.0;

    for(i=0;i < size;++i)
      {
        mat[0][0]+=r[3*i]  *f[3*i];
        mat[1][0]+=r[3*i+1]*f[3*i];
        mat[2][0]+=r[3*i+2]*f[3*i];
        mat[0][1]+=r[3*i]  *f[3*i+1];
        mat[1][1]+=r[3*i+1]*f[3*i+1];
        mat[2][1]+=r[3*i+2]*f[3*i+1];
        mat[0][2]+=r[3*i]  *f[3*i+2];
        mat[1][2]+=r[3*i+1]*f[3*i+2];
        mat[2][2]+=r[3*i+2]*f[3*i+2];
      }

    d2=mat[0][0]*(mat[1][1]*mat[2][2]-mat[1][2]*mat[2][1])
      -mat[0][1]*(mat[1][0]*mat[2][2]-mat[1][2]*mat[2][0])
      +mat[0][2]*(mat[1][0]*mat[2][1]-mat[1][1]*mat[2][0]);

    /* square matrix= ((mat transpose) * mat) */
    for(i=0;i<3;++i)
      for(j=0;j<3;++j)
        {
          x=mat[0][i]*mat[0][j]+mat[1][i]*mat[1][j]+mat[2][i]*mat[2][j];
          mat2[i][j]=mat[i][j];
          rmat[i][j]=x;
        }
    get_roots_3_3(rmat,roots);

    roots[0]=(roots[0]<0.0001) ? 0.0: (roots[0]);
    roots[1]=(roots[1]<0.0001) ? 0.0: (roots[1]);
    roots[2]=(roots[2]<0.0001) ? 0.0: (roots[2]);

    /* make sqrt of rmat, store in mat*/

    roots[0]=(roots[0]<0.0001) ? 0.0: 1.0/(double)sqrt(roots[0]);
    roots[1]=(roots[1]<0.0001) ? 0.0: 1.0/(double)sqrt(roots[1]);
    roots[2]=(roots[2]<0.0001) ? 0.0: 1.0/(double)sqrt(roots[2]);

    if(d2<0.0)
      {
        if( (roots[0]>=roots[1]) && (roots[0]>=roots[2]) )
          roots[0]*=-1.0;
        if( (roots[1]>roots[0]) && (roots[1]>=roots[2]) )
          roots[1]*=-1.0;
        if( (roots[2]>roots[1]) && (roots[2]>roots[0]) )
          roots[2]*=-1.0;
      }

    for(i=0;i<3;++i)
      for(j=0;j<3;++j)
        mat[i][j]=roots[0]*rmat[i][0]*rmat[j][0]+
          roots[1]*rmat[i][1]*rmat[j][1]+
          roots[2]*rmat[i][2]*rmat[j][2];

    /* and multiply into original inertial cross matrix, mat2 */
    for(i=0;i<3;++i)
      for(j=0;j<3;++j)
        rmat[i][j]=mat[0][j]*mat2[i][0]+
          mat[1][j]*mat2[i][1]+
          mat[2][j]*mat2[i][2];

    rvec[0] = rmat[0][0];
    rvec[1] = rmat[0][1];
    rvec[2] = rmat[0][2];
    rvec[3] = rmat[1][0];
    rvec[4] = rmat[1][1];
    rvec[5] = rmat[1][2];
    rvec[6] = rmat[2][0];
    rvec[7] = rmat[2][1];
    rvec[8] = rmat[2][2];
  }

} // end namespace OpenBabel

//! \file obutil.cpp
//! \brief Various utility methods.

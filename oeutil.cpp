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
#include "Vector.h"
#include "oeutil.h"
#include "oeifstream.h"

#ifdef WIN32
#include <conio.h>
#endif

namespace OpenEye {

void ThrowError(char *str)
{
//#ifdef WIN32
//	_cputs(str);
//#else
//    cerr << str << endl;
    cout << str << endl;
//#endif
}

void ThrowError(string &str)
{
//#ifdef WIN32
//	_cputs(str.c_str());
//	_cputs("\n");
//#else
//    cerr << str << endl;
    cout << str << endl;
//#endif
}

void PauseExit()
{
#ifdef WIN32
	char buf[1];
	_cputs("program complete, hit <enter> to continue");
	_cgets(buf);
	exit(0);
#else
	exit(0);
#endif
}


string NewExtension(string &src,char *ext)
{
  unsigned int pos;
  string dst;
  
  if ((pos = src.find_last_of(".")) != string::npos)
      dst = src.substr(0,pos+1);
  else
    {dst = src; dst += ".";}
  
  dst += ext;
  return(dst);
}

Vector center_coords(float *c,int size)
{
  int i;
  float x=0,y=0,z=0;
  for (i = 0;i < size;i++)
    {
      x += c[i*3];
      y += c[i*3+1];
      z += c[i*3+2];
    }
  x /= (float) size; y /= (float) size; z /= (float) size;
  for (i = 0;i < size;i++)
    {
      c[i*3]   -= x;
      c[i*3+1] -= y;
      c[i*3+2] -= z;
    }
  Vector v(x,y,z);
  return(v);
}

void rotate_coords(float *c,float m[3][3],int size)
{
  int i;
  float x,y,z;
  for (i = 0;i < size;i++)
    {
      x = c[i*3]*m[0][0] + c[i*3+1]*m[0][1] + c[i*3+2]*m[0][2];
      y = c[i*3]*m[1][0] + c[i*3+1]*m[1][1] + c[i*3+2]*m[1][2];
      z = c[i*3]*m[2][0] + c[i*3+1]*m[2][1] + c[i*3+2]*m[2][2];
      c[i*3] = x; c[i*3+1] = y; c[i*3+2] = z;
    }
}

float calc_rms(float *r,float *f,int size)
{
  int i;
  float d2=0.0f;
  for (i = 0;i < size;i++)
    {
      d2 += SQUARE(r[i*3] - f[i*3]) +
	   SQUARE(r[i*3+1] - f[i*3+1]) +
	   SQUARE(r[i*3+2] - f[i*3+2]);
    }

  d2 /= (float) size;
  return(sqrt(d2));
}

void SetRotorToAngle(float *c,vector<int> &tor,float ang,vector<int> &atoms)
     //this function will rotate the coordinates of 'atoms'
     //such that tor == ang - atoms in 'tor' should be ordered such 
     //that the 3rd atom is the pivot around which atoms rotate
{
  float v1x,v1y,v1z,v2x,v2y,v2z,v3x,v3y,v3z;
  float c1x,c1y,c1z,c2x,c2y,c2z,c3x,c3y,c3z;
  float c1mag,c2mag,radang,costheta,m[9];
  float x,y,z,mag,rotang,sn,cs,t,tx,ty,tz;

  //
  //calculate the torsion angle
  //
  v1x = c[tor[0]]   - c[tor[1]];   v2x = c[tor[1]]   - c[tor[2]];
  v1y = c[tor[0]+1] - c[tor[1]+1]; v2y = c[tor[1]+1] - c[tor[2]+1];
  v1z = c[tor[0]+2] - c[tor[1]+2]; v2z = c[tor[1]+2] - c[tor[2]+2];
  v3x = c[tor[2]]   - c[tor[3]];
  v3y = c[tor[2]+1] - c[tor[3]+1];
  v3z = c[tor[2]+2] - c[tor[3]+2];

  c1x = v1y*v2z - v1z*v2y;   c2x = v2y*v3z - v2z*v3y;
  c1y = -v1x*v2z + v1z*v2x;  c2y = -v2x*v3z + v2z*v3x;
  c1z = v1x*v2y - v1y*v2x;   c2z = v2x*v3y - v2y*v3x;
  c3x = c1y*c2z - c1z*c2y;
  c3y = -c1x*c2z + c1z*c2x;
  c3z = c1x*c2y - c1y*c2x; 
  
  c1mag = SQUARE(c1x)+SQUARE(c1y)+SQUARE(c1z);
  c2mag = SQUARE(c2x)+SQUARE(c2y)+SQUARE(c2z);
  if (c1mag*c2mag < 0.01) costheta = 1.0; //avoid div by zero error
  else costheta = (c1x*c2x + c1y*c2y + c1z*c2z)/(sqrt(c1mag*c2mag));

  if (costheta < -0.999999) costheta = -0.999999f;
  if (costheta >  0.999999) costheta =  0.999999f;
			      
  if ((v2x*c3x + v2y*c3y + v2z*c3z) > 0.0) radang = -acos(costheta);
  else                                     radang = acos(costheta);

  //
  // now we have the torsion angle (radang) - set up the rot matrix
  //

  //find the difference between current and requested
  rotang = ang - radang; 

  sn = sin(rotang); cs = cos(rotang);t = 1 - cs;
  //normalize the rotation vector
  mag = sqrt(SQUARE(v2x)+SQUARE(v2y)+SQUARE(v2z));
  x = v2x/mag; y = v2y/mag; z = v2z/mag;
  
  //set up the rotation matrix
  m[0]= t*x*x + cs;     m[1] = t*x*y + sn*z;  m[2] = t*x*z - sn*y;
  m[3] = t*x*y - sn*z;  m[4] = t*y*y + cs;    m[5] = t*y*z + sn*x;
  m[6] = t*x*z + sn*y;  m[7] = t*y*z - sn*x;  m[8] = t*z*z + cs;

  //
  //now the matrix is set - time to rotate the atoms
  //
  tx = c[tor[1]];ty = c[tor[1]+1];tz = c[tor[1]+2];
  vector<int>::iterator i;int j;
  for (i = atoms.begin();i != atoms.end();i++)
    {
      j = *i;
      c[j] -= tx;c[j+1] -= ty;c[j+2]-= tz;
      x = c[j]*m[0] + c[j+1]*m[1] + c[j+2]*m[2];
      y = c[j]*m[3] + c[j+1]*m[4] + c[j+2]*m[5];
      z = c[j]*m[6] + c[j+1]*m[7] + c[j+2]*m[8];
      c[j] = x; c[j+1] = y; c[j+2] = z;
      c[j] += tx;c[j+1] += ty;c[j+2] += tz;
    }
}

bool SafeOpen(ifstream &fs,char *filename)
{
#ifdef WIN32
  string s = filename;
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
      ThrowError(error);
      return(false);
    }

  return(true);
}

bool SafeOpen(oeifstream &fs,char *filename)
{
#ifdef WIN32
	string s = filename;
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
      ThrowError(error);
      return(false);
    }

  return(true);
}

bool SafeOpen(ofstream &fs,char *filename)
{
#ifdef WIN32
	string s = filename;
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
      ThrowError(error);
      return(false);
    }

  return(true);
}

bool SafeOpen(ifstream &fs,string &filename)
{
  return(SafeOpen(fs,(char*)filename.c_str()));
}

bool SafeOpen(oeifstream &fs,string &filename)
{
  return(SafeOpen(fs,(char*)filename.c_str()));
}

bool SafeOpen(ofstream &fs,string &filename)
{
  return(SafeOpen(fs,(char*)filename.c_str()));
}

void ToUpper(string &s)
{
  if (s.empty()) return;
  unsigned int i;
  for (i = 0;i < s.size();i++)
    if (isalpha(s[i]) && !isdigit(s[i]))
      s[i] = toupper(s[i]);
}

void ToUpper(char *cptr)
{
  char *c;
  for (c = cptr;*c != '\0';c++)
    if (isalpha(*c) && !isdigit(*c))
      *c = toupper(*c);
}

void CleanAtomType(char *id)
{
  id[0] = toupper(id[0]);
  id[1] = tolower(id[1]);
  if (isalpha(id[1]) == 0) id[1] = '\0';
  else                     id[2] = '\0';
}   

bool SetInputType(OEMol &mol,string &fname)
{
  io_type format;

  if ((format = extab.FilenameToType((char*)fname.c_str())) == UNDEFINED)
    {
      string err = "Error - Unrecognized input format of file '"; 
      err += fname; err += "'";
      ThrowError(err); 
      return(false);
    }
  mol.SetInputType(format);

  return(true);
}

bool SetOutputType(OEMol &mol,string &fname)
{
  io_type format;

  if ((format = extab.FilenameToType((char*)fname.c_str())) == UNDEFINED)
    {
      string err = "Error - Unrecognized input format of file '"; 
      err += fname; err += "'";
      ThrowError(err); 
      return(false);
    }
  mol.SetOutputType(format);

  return(true);
}


void InternalToCartesian(vector<OEInternalCoord*> &vic,OEMol &mol)
{
  Vector n,nn,v1,v2,v3;
  OEAtom *atom;
  vector<OEAtom*>::iterator i;
  int index;

  for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
  {
    index = atom->GetIdx() - 1;

    if      (index == 0)
      {
	atom->SetVector(0.0f, 0.0f, 0.0f);
	continue;
      }
    else if (index == 1)
      {
	v1.SetX(-vic[index]->_dst);
	atom->SetVector(v1);
	continue;
      }
    else if (index == 2) 
      {
	v1.SetX(-(vic[index]->_dst * cos(vic[index]->_ang)));
	v1.SetZ(-(vic[index]->_dst * sin(vic[index]->_ang)));
	atom->SetVector(v1);
	continue;
      }

    v1 = vic[index]->_a->GetVector() - vic[index]->_b->GetVector();
    v2 = vic[index]->_a->GetVector() - vic[index]->_c->GetVector();
    n = cross(v1,v2);  nn = cross(v1,n);
    n.normalize();     nn.normalize();

    n  *= -sin(vic[index]->_tor);  nn *= cos(vic[index]->_tor); 
    v3 = n + nn; v3.normalize(); v3 *= vic[index]->_dst * sin(vic[index]->_ang);
    v1.normalize();  v1 *= vic[index]->_dst * cos(vic[index]->_ang);
    v2 = vic[index]->_a->GetVector() + v3 - v1;

    atom->SetVector(v2);
  }

  // Delete dummy atoms
  for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
    if (atom->GetAtomicNum() == 0)
      mol.DeleteAtom(atom);
}

void CartesianToInternal(vector<OEInternalCoord*> &vic,OEMol &mol)
{
  float r,sum;
  OEAtom *atom,*nbr,*ref;
  vector<OEAtom*>::iterator i,j,m;

  //set reference atoms
  for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
  {
    if      (atom->GetIdx() == 1) continue;
    else if (atom->GetIdx() == 2)
      {
	vic[atom->GetIdx()]->_a = mol.GetAtom(1);
	continue;
      }
    else if (atom->GetIdx() == 3) 
      {
	vic[atom->GetIdx()]->_a = mol.GetAtom(2);
	vic[atom->GetIdx()]->_b = mol.GetAtom(1);
	continue;
      }

    sum=1.0E10;
    ref = mol.GetAtom(1);
    for(nbr = mol.BeginAtom(j);nbr && i != j;nbr = mol.NextAtom(j))
      {  
	if (nbr->GetIdx() < 3) continue;
	r = (atom->GetVector()-nbr->GetVector()).length_2();
	if((r < sum) && (vic[nbr->GetIdx()]->_a != nbr) && 
	   (vic[nbr->GetIdx()]->_b != nbr))
	  {
	    sum = r;
	    ref = nbr;
	  }
      }

    vic[atom->GetIdx()]->_a = ref;
    vic[atom->GetIdx()]->_b = vic[ref->GetIdx()]->_a;
    vic[atom->GetIdx()]->_c = vic[ref->GetIdx()]->_b;
  }

  //fill in geometries
  int k;
  Vector v1,v2;
  OEAtom *a,*b,*c;
  for (k = 2;k <= (signed)mol.NumAtoms();k++)
    {
      atom = mol.GetAtom(k); 
      a = vic[k]->_a;
      b = vic[k]->_b;
      c = vic[k]->_c;
      if (k == 2)
	{
	  vic[k]->_dst = (atom->GetVector() - a->GetVector()).length();
	  continue;
	}

      v1 = atom->GetVector() - a->GetVector();
      v2 = b->GetVector()    - a->GetVector(); 
      vic[k]->_dst = v1.length();
      vic[k]->_ang = VectorAngle(v1,v2);

      if (k == 3) continue;
      vic[k]->_tor = CalcTorsionAngle(atom->GetVector(),
				      a->GetVector(),
				      b->GetVector(),
				      c->GetVector());
    }

  //check for linear geometries and try to correct if possible
  bool done;
  float ang;
  for (k = 2;k <= (signed)mol.NumAtoms();k++)
    {
      ang = fabs(vic[k]->_ang);
      if (ang > 5.0f && ang < 175.0f) continue;
      atom = mol.GetAtom(k);
      done = false;
      for (a = mol.BeginAtom(i);a && a->GetIdx() < k && !done;a = mol.NextAtom(i))
	for (b=mol.BeginAtom(j);b && b->GetIdx()<a->GetIdx() && !done;b = mol.NextAtom(j))
	  {
	    v1 = atom->GetVector() - a->GetVector();
	    v2 = b->GetVector() - a->GetVector();
	    ang = fabs(VectorAngle(v1,v2));
	    if (ang < 5.0f || ang > 175.0f) continue;
	    
	    for (c = mol.BeginAtom(m);c && c->GetIdx() < atom->GetIdx();c = mol.NextAtom(m))
	      if (c != atom && c != a && c != b)
		break;
	    if (!c) continue;

	    vic[k]->_a = a;
	    vic[k]->_b = b;
	    vic[k]->_c = c;
	    vic[k]->_dst = v1.length();
	    vic[k]->_ang = VectorAngle(v1,v2);
	    vic[k]->_tor = CalcTorsionAngle(atom->GetVector(),
					    a->GetVector(),
					    b->GetVector(),
					    c->GetVector());
	    done = true;
	  }
    }
}

bool OECompareInt(const int &a,const int &b)
{
	return(a<b);
}

bool OECompareUnsigned(const unsigned int &a,const unsigned int &b)
{
	return(a<b);
}


void SmartsLexReplace(string &s,vector<pair<string,string> > &vlex)
{
  size_t j,pos;
  string token,repstr;
  vector<pair<string,string> >::iterator i;

  for (pos = 0,pos = s.find("$",pos);pos < s.size();pos = s.find("$",pos))
    //for (pos = 0,pos = s.find("$",pos);pos != string::npos;pos = s.find("$",pos))
    {
      pos++;
      for (j = pos;j < s.size();j++)
        if (!isalpha(s[j]) && !isdigit(s[j]) && s[j] != '_')
          break;
      if (pos == j) continue;

      token = s.substr(pos,j-pos);
      for (i = vlex.begin();i != vlex.end();i++)
        if (token == i->first)
          {
            repstr = "(" + i->second + ")";
            s.replace(pos,j-pos,repstr);
            j = 0;
			break;
          }
      pos = j;
    }
}

//******************triple template*************************
//based on the STL design of the pair<> template

//comparison
template<class T1, class T2, class T3>
bool operator==(const triple<T1,T2,T3> &a, const triple<T1,T2,T3> &b)
{
	return a.first == b.first && a.second == b.second && a.third == b.third;
}

template<class T1, class T2, class T3>
bool operator< (const triple<T1,T2,T3> &a, const triple<T1,T2,T3> &b)
{
	return a.first < b.first ||
				 ((a.first==b.first) && a.second < b.second) ||
				 ((a.first==b.first && a.second==b.second) && a.third < b.third);
}

template<class T1, class T2, class T3>
bool operator!=(const triple<T1,T2,T3> &a, const triple<T1,T2,T3> &b)
{
	return a.first != b.first || a.second != b.second || a.third != b.third;
}

template<class T1, class T2, class T3>
bool operator<=(const triple<T1,T2,T3> &a, const triple<T1,T2,T3> &b)
{
	return a.first <= b.first ||
				 ((a.first==b.first) && a.second <= b.second) ||
				 ((a.first==b.first && a.second==b.second) && a.third <= b.third);
}

template<class T1, class T2, class T3>
bool operator> (const triple<T1,T2,T3> &a, const triple<T1,T2,T3> &b)
{
	return a.first > b.first ||
				 ((a.first==b.first) && a.second > b.second) ||
				 ((a.first==b.first && a.second==b.second) && a.third > b.third);
}

template<class T1, class T2, class T3>
bool operator>=(const triple<T1,T2,T3> &a, const triple<T1,T2,T3> &b)
{
	return a.first >= b.first ||
				 ((a.first==b.first) && a.second >= b.second) ||
				 ((a.first==b.first && a.second==b.second) && a.third >= b.third);
}

//convenience creation
template<class T1, class T2, class T3>
triple<T1,T2,T3> make_triple (const T1&a, const T2 &b, const T3 &c)
{
	return triple<T1,T2,T3>(a,b,c);
}

//**************quad template********************
//based on the design of the STL pair<> template

//comparison
template <class T1, class T2, class T3, class T4>
bool operator==(const quad<T1,T2,T3,T4> &a, const quad<T1,T2,T3,T4> &b)
{
	return a.first == b.first && a.second == b.second && 
				 a.third == b.third && a.fourth == b.fourth;
}

template <class T1, class T2, class T3, class T4>
bool operator< (const quad<T1,T2,T3,T4> &a, const quad<T1,T2,T3,T4> &b)
{
	return a.first < b.first ||
				 ((a.first==b.first) && a.second < b.second) ||
				 ((a.first==b.first && a.second==b.second) && a.third < b.third) ||
				 ((a.first==b.first && a.second==b.second && a.third==b.third) && a.fourth < b.fourth);
}

template <class T1, class T2, class T3, class T4>
bool operator!=(const quad<T1,T2,T3,T4> &a, const quad<T1,T2,T3,T4> &b)
{
	return a.first != b.first || a.second != b.second ||
				 a.third != b.third || a.fourth != b.fourth;
}

template <class T1, class T2, class T3, class T4>
bool operator<=(const quad<T1,T2,T3,T4> &a, const quad<T1,T2,T3,T4> &b)
{
	return a.first <= b.first ||
				 ((a.first==b.first) && a.second <= b.second) ||
				 ((a.first==b.first && a.second==b.second) && a.third <= b.third) ||
				 ((a.first==b.first && a.second==b.second && a.third==b.third) && a.fourth <= b.fourth);
}

template <class T1, class T2, class T3, class T4>
bool operator> (const quad<T1,T2,T3,T4> &a, const quad<T1,T2,T3,T4> &b)
{
	return a.first > b.first ||
				 ((a.first==b.first) && a.second > b.second) ||
				 ((a.first==b.first && a.second==b.second) && a.third > b.third) ||
				 ((a.first==b.first && a.second==b.second && a.third==b.third) && a.fourth > b.fourth);
}

template <class T1, class T2, class T3, class T4>
bool operator>=(const quad<T1,T2,T3,T4> &a, const quad<T1,T2,T3,T4> &b)
{
	return a.first >= b.first ||
				 ((a.first==b.first) && a.second >= b.second) ||
				 ((a.first==b.first && a.second==b.second) && a.third >= b.third) ||
				 ((a.first==b.first && a.second==b.second && a.third==b.third) && a.fourth >= b.fourth);

}

//convenience creation
template <class T1, class T2, class T3, class T4>
quad<T1,T2,T3,T4> make_quad (const T1&a, const T2 &b, const T3 &c, const T4 &d)
{
	return quad<T1,T2,T3,T4>(a,b,c,d);
}

} // Namespace

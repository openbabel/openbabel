/**********************************************************************
Copyright (C) 1998-2000 by OpenEye Scientific Software, Inc.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h> /* for memset */

#ifdef WIN32
#pragma warning (disable: 4244) // warning: conversion from 'double' to 'float', possible loss of data
#endif

#define MAX_SWEEPS 30

void jacobi(float a[4][4], float *d, float v[4][4]);

void qtrfit (float *r,float *f,int size, float u[3][3])
{
  register int i;
  float xxyx, xxyy, xxyz;
  float xyyx, xyyy, xyyz;
  float xzyx, xzyy, xzyz;
  float d[4],q[4];
  float c[4][4],v[4][4];
  float rx,ry,rz,fx,fy,fz;

/* generate the upper triangle of the quadratic form matrix */

 xxyx = 0.0; xxyy = 0.0; xxyz = 0.0;
 xyyx = 0.0; xyyy = 0.0; xyyz = 0.0;
 xzyx = 0.0; xzyy = 0.0; xzyz = 0.0;
 
 for (i = 0; i < size; i++) 
 {
   rx = r[i*3];   ry = r[i*3+1];   rz = r[i*3+2];
   fx = f[i*3];   fy = f[i*3+1];   fz = f[i*3+2];

   xxyx += fx * rx;xxyy += fx * ry;xxyz += fx * rz;
   xyyx += fy * rx;xyyy += fy * ry;xyyz += fy * rz;
   xzyx += fz * rx;xzyy += fz * ry;xzyz += fz * rz;
 }

 c[0][0] = xxyx + xyyy + xzyz;

 c[0][1] = xzyy - xyyz;
 c[1][1] = xxyx - xyyy - xzyz;

 c[0][2] = xxyz - xzyx;
 c[1][2] = xxyy + xyyx;
 c[2][2] = xyyy - xzyz - xxyx;

 c[0][3] = xyyx - xxyy;
 c[1][3] = xzyx + xxyz;
 c[2][3] = xyyz + xzyy;
 c[3][3] = xzyz - xxyx - xyyy;

/* diagonalize c */

 jacobi(c, d, v);

/* extract the desired quaternion */

 q[0] = v[0][3];
 q[1] = v[1][3];
 q[2] = v[2][3];
 q[3] = v[3][3];

/* generate the rotation matrix */

 u[0][0] = q[0]*q[0] + q[1]*q[1] - q[2]*q[2] - q[3]*q[3];
 u[1][0] = 2.0f * (q[1] * q[2] - q[0] * q[3]);
 u[2][0] = 2.0f * (q[1] * q[3] + q[0] * q[2]);

 u[0][1] = 2.0f * (q[2] * q[1] + q[0] * q[3]);
 u[1][1] = q[0]*q[0] - q[1]*q[1] + q[2]*q[2] - q[3]*q[3];
 u[2][1] = 2.0f * (q[2] * q[3] - q[0] * q[1]);

 u[0][2] = 2.0f *(q[3] * q[1] - q[0] * q[2]);
 u[1][2] = 2.0f * (q[3] * q[2] + q[0] * q[1]);
 u[2][2] = q[0]*q[0] - q[1]*q[1] - q[2]*q[2] + q[3]*q[3];
}

void jacobi(float a[4][4], float *d, float v[4][4])
{
  float onorm, dnorm;
  float b, dma, q, t, c, s;
  float  atemp, vtemp, dtemp;
  register int i, j, k, l;
  int nrot;

  for (j = 0; j <= 3; j++)
  {
    for (i = 0; i <= 3; i++) v[i][j] = 0.0;
    v[j][j] = 1.0;
    d[j] = a[j][j];
  }

  nrot = MAX_SWEEPS;
  for (l = 1; l <= nrot; l++)
  {
    dnorm = 0.0;
    onorm = 0.0;
    for (j = 0; j <= 3; j++)
    {
      dnorm = dnorm + (float)fabs(d[j]);
      for (i = 0; i <= j - 1; i++)
	 onorm = onorm + (float)fabs(a[i][j]);
    }

    if((onorm/dnorm) <= 1.0e-12) goto Exit_now;
    for (j = 1; j <= 3; j++)
    {
      for (i = 0; i <= j - 1; i++)
      {
	b = a[i][j];
	if(fabs(b) > 0.0) {
	  dma = d[j] - d[i];
	  if((fabs(dma) + fabs(b)) <=  fabs(dma))
	     t = b / dma;
	  else
	  {
	    q = 0.5f * dma / b;
	    t = 1.0f/((float)fabs(q) + (float)sqrt(1.0f+q*q));
	    if(q < 0.0f) t = -t;

	  }

	  c = 1.0f/(float)sqrt(t * t + 1.0f);
	  s = t * c;
	  a[i][j] = 0.0f;
	  for (k = 0; k <= i-1; k++)
	  {
	    atemp = c * a[k][i] - s * a[k][j];
	    a[k][j] = s * a[k][i] + c * a[k][j];
	    a[k][i] = atemp;
	  }

	  for (k = i+1; k <= j-1; k++) 
	  {
	    atemp = c * a[i][k] - s * a[k][j];
	    a[k][j] = s * a[i][k] + c * a[k][j];
	    a[i][k] = atemp;
	  }

	  for (k = j+1; k <= 3; k++)
	  {
	    atemp = c * a[i][k] - s * a[j][k];
	    a[j][k] = s * a[i][k] + c * a[j][k];
	    a[i][k] = atemp;
	  }

	  for (k = 0; k <= 3; k++)
	  {
	    vtemp = c * v[k][i] - s * v[k][j];
	    v[k][j] = s * v[k][i] + c * v[k][j];
	    v[k][i] = vtemp;
	  }

	  dtemp = c*c*d[i] + s*s*d[j] - 2.0f*c*s*b;
	  d[j] = s*s*d[i] + c*c*d[j] +  2.0f*c*s*b;
	  d[i] = dtemp;
	}  /* end if */
      } /* end for i */
    } /* end for j */
  } /* end for l */
 
  Exit_now:

  nrot = l;

  for (j = 0; j <= 2; j++)
  {
    k = j;
    dtemp = d[k];
    for (i = j+1; i <= 3; i++)
       if(d[i] < dtemp)
       {
	 k = i;
	 dtemp = d[k];
       }


    if(k > j)
    {
      d[k] = d[j];
      d[j] = dtemp;
      for (i = 0; i <= 3; i++)
      {
	dtemp = v[i][k];
	v[i][k] = v[i][j];
	v[i][j] = dtemp;
      }
    }
  }
}


static double Roots[4];

#define ApproxZero 1E-6
#define IsZero(x)  ((float)fabs(x)<ApproxZero)
#ifndef PI
#define PI         3.1415926536
#endif
#define OneThird      (1.0/3.0)
#define FourThirdsPI  (4.0*PI/3.0)
#define TwoThirdsPI   (2.0*PI/3.0)

#ifdef OLD_RMAT

/*FUNCTION */
/* recieves: the co-efficients for a general
 *           equation of degree one.
 *           Ax + B = 0 !!
 */
static int SolveLinear(double A,double B) 
{
    if( IsZero(A) )
	return( 0 );
    Roots[0] = -B/A;
    return( 1 );
}

/*FUNCTION */
/* recieves: the co-efficients for a general
 *           linear equation of degree two.
 *           Ax^2 + Bx + C = 0 !!
 */
static int SolveQuadratic(double A,double B,double C)
{
    register double Descr, Temp, TwoA;

    if( IsZero(A) )
	return( SolveLinear(B,C) );

    TwoA = A+A;
    Temp = TwoA*C;
    Descr = B*B - (Temp+Temp);
    if( Descr<0.0 )
	return( 0 );

    if( Descr>0.0 )
    {   Descr = sqrt(Descr);
#ifdef ORIG
	Roots[0] = (-B-Descr)/TwoA;
	Roots[1] = (-B+Descr)/TwoA;
#else
	/* W. Press, B. Flannery, S. Teukolsky and W. Vetterling,
	 * "Quadratic and Cubic Equations", Numerical Recipes in C, 
	 * Chapter 5, pp. 156-157, 1989.
	 */
	Temp = (B<0.0)? -0.5*(B-Descr) : -0.5*(B+Descr);
	Roots[0] = Temp/A;
	Roots[1] = C/Temp;
#endif
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
static double CubeRoot(double X)
{
    if( X>=0.0 )
    {   return pow( X, OneThird );
    } else
	return -pow( -X, OneThird );
}

static int SolveCubic(double A,double B,double C,double D)
{
    register double TwoA, ThreeA, BOver3A;
    register double Temp, POver3, QOver2;
    register double Desc, Rho, Psi;


    if( IsZero(A) ){
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
    {   Rho = sqrt( -Rho );
	Psi = OneThird*acos(-QOver2/Rho);
	Temp = CubeRoot( Rho );
	Temp = Temp+Temp;

	Roots[0] = Temp*cos( Psi )-BOver3A;
	Roots[1] = Temp*cos( Psi+TwoThirdsPI )-BOver3A;
	Roots[2] = Temp*cos( Psi+FourThirdsPI )-BOver3A;
	return( 3 );
    }

    if( Desc> 0.0 )
    {   Temp = CubeRoot( -QOver2 );
	Roots[0] = Temp+Temp-BOver3A;
	Roots[1] = -Temp-BOver3A;
	return( 2 );
    }

    Desc = sqrt( Desc );
    Roots[0] = CubeRoot(Desc-QOver2)-CubeRoot(Desc+QOver2) - BOver3A;

    return( 1 );
}
#endif

void oe_make_rmat(float a[3][3],float rmat[9])
{
  float onorm, dnorm;
  float b, dma, q, t, c, s,d[3];
  float atemp, vtemp, dtemp,v[3][3];
  float r1[3],r2[3],v1[3],v2[3],v3[3];
  int i, j, k, l;

  memset((char*)d,'\0',sizeof(float)*3);
  
  for (j = 0; j < 3; j++) 
    {
      for (i = 0; i < 3; i++) v[i][j] = 0.0;

      v[j][j] = 1.0;
      d[j] = a[j][j];
    }
  
  for (l = 1; l <= MAX_SWEEPS; l++) 
    {
      dnorm = 0.0;
      onorm = 0.0;
      for (j = 0; j < 3; j++) 
	{
	  dnorm = dnorm + (float)fabs(d[j]);
	  for (i = 0; i <= j - 1; i++) 
	    {
	      onorm = onorm + (float)fabs(a[i][j]);
	    }
	}
      
      if((onorm/dnorm) <= 1.0e-12) goto Exit_now;
      for (j = 1; j < 3; j++) 
	{
	  for (i = 0; i <= j - 1; i++) 
	    {
	      b = a[i][j];
	      if(fabs(b) > 0.0) 
		{
		  dma = d[j] - d[i];
		  if((fabs(dma) + fabs(b)) <=  fabs(dma)) 
		    t = b / dma;
		  else 
		    {
		      q = 0.5f * dma / b;
		      t = 1.0f/((float)fabs(q) + (float)sqrt(1.0+q*q));
		      if(q < 0.0) t = -t;
		    }
		  c = 1.0f/(float)sqrt(t * t + 1.0f);
		  s = t * c;
		  a[i][j] = 0.0f;
		  for (k = 0; k <= i-1; k++) 
		    {
		      atemp = c * a[k][i] - s * a[k][j];
		      a[k][j] = s * a[k][i] + c * a[k][j];
		      a[k][i] = atemp;
		    }
		  for (k = i+1; k <= j-1; k++) 
		    {
		      atemp = c * a[i][k] - s * a[k][j];
		      a[k][j] = s * a[i][k] + c * a[k][j];
		      a[i][k] = atemp;
		    }
		  for (k = j+1; k < 3; k++) 
		    {
		      atemp = c * a[i][k] - s * a[j][k];
		      a[j][k] = s * a[i][k] + c * a[j][k];
		      a[i][k] = atemp;
		    }
		  for (k = 0; k < 3; k++) 
		    {
		      vtemp = c * v[k][i] - s * v[k][j];
		      v[k][j] = s * v[k][i] + c * v[k][j];
		      v[k][i] = vtemp;
		    }
		  dtemp = c*c*d[i] + s*s*d[j] - 2.0f*c*s*b;
		  d[j] = s*s*d[i] + c*c*d[j] +  2.0f*c*s*b;
		  d[i] = dtemp;
		}  /* end if */
	    } /* end for i */
	} /* end for j */
    } /* end for l */
  
Exit_now:
  
  /* max_sweeps = l;*/
  
  for (j = 0; j < 3-1; j++) 
    {
      k = j;
      dtemp = d[k];
      for (i = j+1; i < 3; i++) 
	if(d[i] < dtemp) 
	  {k = i;dtemp = d[k];}

      if(k > j) 
	{
	  d[k] = d[j];
	  d[j] = dtemp;
	  for (i = 0; i < 3 ; i++) 
	    {
	      dtemp = v[i][k];
	      v[i][k] = v[i][j];
	      v[i][j] = dtemp;
	    }
	}
    }

  r1[0] = v[0][0]; r1[1] = v[1][0]; r1[2] = v[2][0];
  r2[0] = v[0][1]; r2[1] = v[1][1]; r2[2] = v[2][1];

  v3[0] =  r1[1]*r2[2] - r1[2]*r2[1];
  v3[1] = -r1[0]*r2[2] + r1[2]*r2[0];
  v3[2] =  r1[0]*r2[1] - r1[1]*r2[0];
  s = (float)sqrt(v3[0]*v3[0] + v3[1]*v3[1] + v3[2]*v3[2]);
  v3[0] /= s; v3[0] /= s; v3[0] /= s;

  v2[0] =  v3[1]*r1[2] - v3[2]*r1[1];
  v2[1] = -v3[0]*r1[2] + v3[2]*r1[0];
  v2[2] =  v3[0]*r1[1] - v3[1]*r1[0];
  s = (float)sqrt(v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2]);
  v2[0] /= s; v2[0] /= s; v2[0] /= s;

  v1[0] =  v2[1]*v3[2] - v2[2]*v3[1];
  v1[1] = -v2[0]*v3[2] + v2[2]*v3[0];
  v1[2] =  v2[0]*v3[1] - v2[1]*v3[0];
  s = (float)sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]);
  v1[0] /= s; v1[0] /= s; v1[0] /= s;

  rmat[0] = v1[0]; rmat[1] = v1[1]; rmat[2] = v1[2];
  rmat[3] = v2[0]; rmat[4] = v2[1]; rmat[5] = v2[2];
  rmat[6] = v3[0]; rmat[7] = v3[1]; rmat[8] = v3[2];
}

static int get_roots_3_3(float mat[3][3], float roots[3])
{
   float rmat[9];

   oe_make_rmat(mat,rmat);

   mat[0][0]=rmat[0];
   mat[0][1]=rmat[3];
   mat[0][2]=rmat[6];
   mat[1][0]=rmat[1];
   mat[1][1]=rmat[4];
   mat[1][2]=rmat[7];
   mat[2][0]=rmat[2];
   mat[2][1]=rmat[5];
   mat[2][2]=rmat[8];

   roots[0]=(float)Roots[0];
   roots[1]=(float)Roots[1];
   roots[2]=(float)Roots[2];

   return 1;
}

float superimpose(float *r,float *f,int size)
{
  int i,j;
  float x,y,z,d2;
  float mat[3][3],rmat[3][3],mat2[3][3],roots[3];

  /* make inertial cross tensor */
  for(i=0;i<3;i++) 
    for(j=0;j<3;j++)
      mat[i][j]=0.0; 

  for(i=0;i < size;i++)
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
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      {
	x=mat[0][i]*mat[0][j]+mat[1][i]*mat[1][j]+mat[2][i]*mat[2][j];
	mat2[i][j]=mat[i][j];
	rmat[i][j]=x;
      }
  get_roots_3_3(rmat,roots);

  roots[0]=(roots[0]<0.0001) ? 0.0f: (roots[0]);
  roots[1]=(roots[1]<0.0001) ? 0.0f: (roots[1]);
  roots[2]=(roots[2]<0.0001) ? 0.0f: (roots[2]);

  /* make sqrt of rmat, store in mat*/

  roots[0]=roots[0]<0.0001? 0.0f: 1.0f/(float)sqrt(roots[0]);
  roots[1]=roots[1]<0.0001? 0.0f: 1.0f/(float)sqrt(roots[1]);
  roots[2]=roots[2]<0.0001? 0.0f: 1.0f/(float)sqrt(roots[2]);

  if(d2<0.0){
    if( (roots[0]>=roots[1]) && (roots[0]>=roots[2]) ) roots[0]*=-1.0f;
    if( (roots[1]>roots[0]) && (roots[1]>=roots[2]) )  roots[1]*=-1.0f;
    if( (roots[2]>roots[1]) && (roots[2]>roots[0]) )   roots[2]*=-1.0f;
  }

  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      mat[i][j]=roots[0]*rmat[i][0]*rmat[j][0]+
	        roots[1]*rmat[i][1]*rmat[j][1]+
	        roots[2]*rmat[i][2]*rmat[j][2];

  /* and multiply into original inertial cross matrix, mat2 */
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      rmat[i][j]=mat[0][j]*mat2[i][0]+
	         mat[1][j]*mat2[i][1]+
	         mat[2][j]*mat2[i][2];

  /* rotate all coordinates */
  d2 = 0.0;
  for(i=0;i<size;i++)
    {
      x=f[3*i]*rmat[0][0]+f[3*i+1]*rmat[0][1]+f[3*i+2]*rmat[0][2];
      y=f[3*i]*rmat[1][0]+f[3*i+1]*rmat[1][1]+f[3*i+2]*rmat[1][2];
      z=f[3*i]*rmat[2][0]+f[3*i+1]*rmat[2][1]+f[3*i+2]*rmat[2][2];
      f[3*i  ]=x; f[3*i+1]=y; f[3*i+2]=z;

      x = r[i*3]   - f[i*3];
      y = r[i*3+1] - f[i*3+1];
      z = r[i*3+2] - f[i*3+2];
      d2 += x*x+y*y+z*z;
    }

  d2 /= (float) size;

  return((float)sqrt(d2));
}

void get_rmat(float *rvec,float *r,float *f,int size)
{
  int i,j;
  float x,d2;
  float mat[3][3],rmat[3][3],mat2[3][3],roots[3];

  /* make inertial cross tensor */
  for(i=0;i<3;i++) 
    for(j=0;j<3;j++)
      mat[i][j]=0.0; 

  for(i=0;i < size;i++)
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
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      {
	x=mat[0][i]*mat[0][j]+mat[1][i]*mat[1][j]+mat[2][i]*mat[2][j];
	mat2[i][j]=mat[i][j];
	rmat[i][j]=x;
      }
  get_roots_3_3(rmat,roots);

  roots[0]=(roots[0]<0.0001f) ? 0.0f: (roots[0]);
  roots[1]=(roots[1]<0.0001f) ? 0.0f: (roots[1]);
  roots[2]=(roots[2]<0.0001f) ? 0.0f: (roots[2]);

  /* make sqrt of rmat, store in mat*/

  roots[0]=(roots[0]<0.0001f) ? 0.0f: 1.0f/(float)sqrt(roots[0]);
  roots[1]=(roots[1]<0.0001f) ? 0.0f: 1.0f/(float)sqrt(roots[1]);
  roots[2]=(roots[2]<0.0001f) ? 0.0f: 1.0f/(float)sqrt(roots[2]);

  if(d2<0.0){
    if( (roots[0]>=roots[1]) && (roots[0]>=roots[2]) ) roots[0]*=-1.0f;
    if( (roots[1]>roots[0]) && (roots[1]>=roots[2]) )  roots[1]*=-1.0f;
    if( (roots[2]>roots[1]) && (roots[2]>roots[0]) )   roots[2]*=-1.0f;
  }

  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      mat[i][j]=roots[0]*rmat[i][0]*rmat[j][0]+
	        roots[1]*rmat[i][1]*rmat[j][1]+
	        roots[2]*rmat[i][2]*rmat[j][2];

  /* and multiply into original inertial cross matrix, mat2 */
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      rmat[i][j]=mat[0][j]*mat2[i][0]+
	         mat[1][j]*mat2[i][1]+
	         mat[2][j]*mat2[i][2];

  rvec[0] = rmat[0][0]; rvec[1] = rmat[0][1]; rvec[2] = rmat[0][2];
  rvec[3] = rmat[1][0]; rvec[4] = rmat[1][1]; rvec[5] = rmat[1][2];
  rvec[6] = rmat[2][0]; rvec[7] = rmat[2][1]; rvec[8] = rmat[2][2];
}


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
#include "binary_io.h"

using namespace std;

namespace OpenBabel {

//test byte ordering
static int SSINT = 0x00000001;
static unsigned char *SSTPTR = (unsigned char*)&SSINT;
bool SSwabInt = (SSTPTR[0]!=0);
static const int TMPSIZE = sizeof(double);

/******************** reading generic binary data *************************/
unsigned int OB_io_read_binary(char* ccc, char *x, unsigned int size, unsigned int count)
{
    if (x == NULL || size == 0)
    	return 0;

    if (size == 1)
    {
    	memcpy(x, ccc, count);
    	return count;
    }
    else
    {
    	unsigned int sz = size * count;
    	if (!SSwabInt)
    		memcpy(x, ccc, sz);
    	else
    	{
    		unsigned int i,j,k;
    		unsigned char c, tmp[TMPSIZE];

    		for ( i = 0 ; i < sz ; i += size )
    		{
    			memcpy(tmp, &ccc[i], size);
    			for (j = 0, k = size - 1 ; j < k ; j++, k--) {
    				c = tmp[j]; tmp[j] = tmp[k]; tmp[k] = c;
    			}
    			memcpy(&x[i], tmp, size);
    		}
    	}
    	return sz;
    }
}
unsigned int OB_io_read_binary(istream& ifs, char *x, unsigned int size, unsigned int count)
{
    if (x == NULL || size == 0)
    	return 0;

    if (size == 1)
    {
    	ifs.read(x, count);
    	return count;
    }
    else
    {
    	unsigned int sz = size * count;
    	if (!SSwabInt)
    		ifs.read(x, sz);
    	else
    	{
    		unsigned int i,j,k;
    		unsigned char c, tmp[TMPSIZE];

    		for ( i = 0 ; i < count ; i++ )
    		{
    			ifs.read((char*)&tmp[0], size);
    			for (j = 0, k = size - 1 ; j < k ; j++, k--) {
    				c = tmp[j]; tmp[j] = tmp[k]; tmp[k] = c;
    			}
    			memcpy(&x[i*size], tmp, size);
    		}
    	}
    	return sz;
    }
}

/******************** writing generic binary data *************************/
unsigned int OB_io_write_binary(char *ccc, const char *x, unsigned int size, unsigned int count)
{
    if (x == NULL || size == 0)
    	return 0;

    if (size == 1)
    {
    	memcpy(ccc, x, count);
    	return count;
    }
    else
    {
    	unsigned int sz = size * count;
    	if (!SSwabInt)
    		memcpy(ccc, x, sz);
    	else
    	{
    		unsigned i,j,k;
    		unsigned char c, tmp[TMPSIZE];

    		for (i = 0 ; i < sz ; i += size)
    		{
    			memcpy(tmp, &x[i], size);
    			for (j = 0, k = size - 1 ; j < k ; j++, k-- ) {
    				c = tmp[j]; tmp[j] = tmp[k]; tmp[k] = c;
    			}
    			memcpy(&ccc[i], tmp, size);
    		}
    	}
    	return sz;
    }
}
unsigned int OB_io_write_binary(ostream &ofs, const char *x, unsigned int size, unsigned int count)
{
    if (x == NULL || size == 0)
    	return 0;

    if (size == 1)
    {
    	ofs.write(x, count);
    	return count;
    }
    else
    {
    	unsigned int sz = size * count;
    	if (!SSwabInt)
    		ofs.write(x, sz);
    	else
    	{
    		unsigned i,j,k;
    		unsigned char c, tmp[TMPSIZE];

    		for (i = 0 ; i < sz ; i += size)
    		{
    			memcpy(tmp, &x[i], size);
    			for (j = 0, k = size - 1 ; j < k ; j++, k-- ) {
    				c = tmp[j]; tmp[j] = tmp[k]; tmp[k] = c;
    			}
    			ofs.write((char*)&tmp[0], size);
    		}
    	}
    	return sz;
    }
}

/**********************reading STL strings****************************/

unsigned int OB_io_read_binary(char* ccc, string& str)
{
    char *buffer;
    unsigned int i = 0, idx = 0;

    idx    += OB_io_read_binary(ccc, (char*) &i, sizeof(unsigned int), 1);
    buffer  = new char[i+1];
    idx    += OB_io_read_binary(&ccc[idx], buffer, sizeof(char), i);

    buffer[i] = '\0';
    str       = buffer;

	delete [] buffer;

    return idx;
}
unsigned int OB_io_read_binary(istream& ifs, string& str)
{
    char *buffer;
    unsigned int i = 0, idx = 0;

    idx    += OB_io_read_binary(ifs, (char*) &i, sizeof(unsigned int), 1);
    buffer  = new char[i+1];
    idx    += OB_io_read_binary(ifs, buffer, sizeof(char), i);

    buffer[i] = '\0';
    str       = buffer;

	delete [] buffer;

    return idx;
}

/******************** writing STL strings *************************/

unsigned int OB_io_write_binary(char* ccc, const string& str)
{
    unsigned int idx = 0, size = str.size();

    idx += OB_io_write_binary(ccc, (char*) &size, sizeof(unsigned int), 1);
    idx += OB_io_write_binary(&ccc[idx], str.c_str(), sizeof(char), size);

    return idx;
}
unsigned int OB_io_write_binary(ostream& ofs, const string& str)
{
    unsigned int idx = 0, size = str.size();

    idx += OB_io_write_binary(ofs, (char*) &size, sizeof(unsigned int), 1);
    idx += OB_io_write_binary(ofs, str.c_str(), sizeof(char), size);

    return idx;
}

/******************* writing compressed data *************************/
/*!
**\fn OB_io_write_binary_compressed(char*ccc, unsigned int *x, unsigned int NumBits, unsigned int NumInts)
**\brief Writes an unsigned interger array out in a platform independed
**compressed binary format.
**\param ccc A character array that the binary is written to
**\param x An unsigned integer array to be written to ccc.
**\param NumBits The number of bits that will be used to store
**each unsigned integer.  The largest value in x must be less than
**2^NumBits.
**\param NumInts The number of integers in x
**\return The number of bytes written to ccc.  This value should
**always be 1+(NumBits*NumInts)/8 (integer division).
*/
unsigned int OB_io_write_binary_compressed(char*ccc, unsigned int *x, unsigned int NumBits, unsigned int NumInts)
  {
    //If there are no values or x is NULL do nothing
    if (!NumInts || x==NULL) return 0;
 
    //If 32 or more bits are specified simply write uncompressed
    if (NumBits >=32) return OB_io_write_binary(ccc,(const char*) x, sizeof(unsigned int), NumInts);

    bool swap = !SSwabInt; 
    unsigned int bitarraysize = 1+(NumBits*NumInts)/8;
    unsigned char *cc = (unsigned char*) ccc;
    unsigned int i,j;
    for (i=0 ; i<bitarraysize ; i++) cc[i] = 0;
    unsigned char mask1[8],mask2[8];
    mask2[7] = 127;
    mask2[6] = 63;
    mask2[5] = 31;
    mask2[4] = 15;
    mask2[3] = 7;
    mask2[2] = 3;
    mask2[1] = 1;
    mask2[0] = 0;
    for (i=0 ; i<8 ; i++) mask1[i] = ~mask2[i];
    unsigned int bitshift=0;
    unsigned int bytes=0;
    unsigned int value;
    unsigned char *cptr = (unsigned char*)&value;
    for (i=0 ; i<NumInts ; i++) {
        if (swap) {
            cptr[0] = ((unsigned char*)(&x[i]))[3]; 
            cptr[1] = ((unsigned char*)(&x[i]))[2]; 
            cptr[2] = ((unsigned char*)(&x[i]))[1]; 
            cptr[3] = ((unsigned char*)(&x[i]))[0]; 
          }
        else      {value = x[i];}
        for (j=0 ; j<NumBits/8 ; j++) {
            cc[bytes+0+j] |= (cptr[j]<<(bitshift  ))&mask1[bitshift];
            cc[bytes+1+j] |= (cptr[j]>>(8-bitshift))&mask2[bitshift];
          }
        cc[bytes+0+j] |= (cptr[j]<<(bitshift  ))&mask1[bitshift];
        if ((bitshift + NumBits%8) > 7) cc[bytes+1+j] |= (cptr[j]>>(8-bitshift))&mask2[bitshift];
        bitshift += NumBits;
        while (bitshift > 7) {bitshift-=8; bytes++;}
      }
    return bitarraysize;
  }
/*!
**\fn OB_io_write_binary_compressed(ostream& ostr, unsigned int *x, unsigned int NumBits, unsigned int NumInts)
**\brief Writes an unsigned interger array out in a platform independed
**compressed binary format.
**\param ostr An output stream that the binary is written to
**\param x An unsigned integer array to be written to ostr.
**\param NumBits The number of bits that will be used to store
**each unsigned integer.  The largest value in x must be less than
**2^NumBits.
**\param NumInts The number of integers in x
**\return The number of bytes written to ostr.  This value should
**always be 1+(NumBits*NumInts)/8 (integer division).
*/
unsigned int OB_io_write_binary_compressed(ostream& ostr, unsigned int *x, unsigned int NumBits, unsigned int NumInts)
  {
    //If there are no values or x is NULL do nothing
    if (!NumInts || x==NULL) return 0;
 
    //If 32 or more bits are specified simply write uncompressed
    if (NumBits >=32) return OB_io_write_binary(ostr,(const char*) x, sizeof(unsigned int), NumInts);

    unsigned int bitarraysize = 1+(NumBits*NumInts)/8;
    char *cc = new char [bitarraysize];
    unsigned int idx = OB_io_write_binary_compressed(cc,x,NumBits,NumInts);
    ostr.write(cc,bitarraysize);
    delete [] cc;
    return idx;
  }

/*!
**\fn OB_io_write_binary_compressed(char*ccc, float *x, unsigned int NumBits, unsigned int Nfloats)
**\brief Writes a float array out in a platform independed
**compressed binary format.
**\param ccc A character array that the binary is written to
**\param x A float array to be written to ccc.
**\param NumBits The number of bits that will be used to store
**each unsigned integer.  The lost of resolution of the float
**values depends directly on how many bits are used.
**\param Nfloats The number of floats in x.
**\return The number of bytes written to ccc.  This value should
**always be 9+(NumBits*NumInts)/8 (integer division).
**\par Compression Info
**There is a loss of resolution of the floats using this routine.
**The greater the value of NumBits and the smaller the range of values
**in x, the smaller the loss of resolution.  The function 
**OB_io_util_calc_NumBits(float *x, unsigned int N, float res)
**will calculate how many bits are required for a given resolution and
**set of floats.
*/
unsigned int OB_io_write_binary_compressed(char*ccc, float *x, unsigned int NumBits, unsigned int Nfloats)
  {
    unsigned int idx=0;

    //If there are no values of x is NULL do nothing
    if (!Nfloats || x==NULL) return 0;

    //If 32 or more bits are specified simply write uncompressed
    if (NumBits >=32) return OB_io_write_binary(&ccc[idx],(const char*) x, sizeof(unsigned int), Nfloats);



    //Write max and min
    unsigned int i;
    float max = x[0];
    float min = x[0];
    for (i=0 ; i<Nfloats ; i++) {
        if (x[i] > max) max = x[i];
        if (x[i] < min) min = x[i];
      }
    idx += OB_io_write_binary(&ccc[idx],(const char*)&max,sizeof(float),1);
    idx += OB_io_write_binary(&ccc[idx],(const char*)&min,sizeof(float),1);

    //Integerize values and write out 
    unsigned int maxint = (unsigned int) 0; 
    for (i=0 ; i<NumBits ; i++) maxint |= ((unsigned int)1)<<i;
    unsigned int *ix = new unsigned int [Nfloats]; 
    if (max!=min) {
        for (i=0 ; i<Nfloats ; i++) 
            ix[i] = (unsigned int) (((float)maxint*((x[i]-min)/(max-min))) + 0.499999999f);
      }
    else for (i=0 ; i<Nfloats ; i++) ix[i] = (unsigned int) 0;
    idx += OB_io_write_binary_compressed(&ccc[idx],ix,NumBits,Nfloats);
    delete [] ix;  

    return idx;
  }
/*!
**\fn OB_io_write_binary_compressed(ostream& ostr, float *x, unsigned int NumBits, unsigned int Nfloats)
**\brief Writes a float array out in a platform independed
**compressed binary format.
**\param ostr A character array that the binary is written to
**\param x A float array to be written to ostr.
**\param NumBits The number of bits that will be used to store
**each unsigned integer.  The lost of resolution of the float
**values depends directly on how many bits are used.
**\param Nfloats The number of floats in x.
**\return The number of bytes written to ccc.  This value should
**always be 9+(NumBits*NumInts)/8 (integer division).
**\par Compression Info
**There is a loss of resolution of the floats using this routine.
**The greater the value of NumBits and the smaller the range of values
**in x, the smaller the loss of resolution.  The function
**OB_io_util_calc_NumBits(float *x, unsigned int N, float res)
**will calculate how many bits are required for a given resolution and
**set of floats.
*/
unsigned int OB_io_write_binary_compressed(ostream& ostr, float *x, unsigned int NumBits, unsigned int Nfloats)
  {
    unsigned int idx=0;

    //If there are no values of x is NULL do nothing
    if (!Nfloats || x==NULL) return 0;

    //If 32 or more bits are specified simply write uncompressed
    if (NumBits >=32) return OB_io_write_binary(ostr,(const char*) x, sizeof(unsigned int), Nfloats);

    //Write max and min
    unsigned int i;
    float max = x[0];
    float min = x[0];
    for (i=0 ; i<Nfloats ; i++) {
        if (x[i] > max) max = x[i];
        if (x[i] < min) min = x[i];
      }
    idx += OB_io_write_binary(ostr,(const char*)&max,sizeof(float),1);
    idx += OB_io_write_binary(ostr,(const char*)&min,sizeof(float),1);
 
    //Integerize values and write out
    unsigned int maxint = (unsigned int) 0;
    for (i=0 ; i<NumBits ; i++) maxint |= ((unsigned int)1)<<i;
    unsigned int *ix = new unsigned int [Nfloats];
    if (max!=min) {
        for (i=0 ; i<Nfloats ; i++)
            ix[i] = (unsigned int) (((float)maxint*((x[i]-min)/(max-min))) + 0.499999999f);
      }
    else for (i=0 ; i<Nfloats ; i++) ix[i] = (unsigned int) 0;
    idx += OB_io_write_binary_compressed(ostr,ix,NumBits,Nfloats);
    delete [] ix;


    return idx;
  }

/******************* reading compressed data *************************/
/*!
**\fn OB_io_read_binary_compressed(char*ccc, unsigned int *x, unsigned int NumBits, unsigned int NumInts)
**\brief Reads compressed binary data that was written with either
**OB_io_write_binary_compressed(char*ccc, unsigned int *x, unsigned int NumBits, unsigned int NumInts)
**or
**OB_io_write_binary_compressed(ostream& ostr, unsigned int *x, unsigned int NumBits, unsigned int NumInts)
**\param ccc A character array the binary is being read from
**\param x An array of unsigned ints the binary data will be extracted to
**\param NumBits The number of bits used to store each unsigned integer
**\param NumInts The number of integers stored in ccc
**\return The number of bytes read from ccc.  This value should
**always be 1+(NumBits*NumInts)/8 (integer division).
*/
unsigned int OB_io_read_binary_compressed(char*ccc, unsigned int *x, unsigned int NumBits, unsigned int NumInts)
  {
    //If there are no values or x is NULL do nothing
    if (!NumInts || x==NULL) return 0;
 
    //If 32 or more bits are specified simply read uncompressed
    if (NumBits >=32) return OB_io_read_binary(ccc,(char*) x, sizeof(unsigned int), NumInts);

    //Read the values
    bool swap = !SSwabInt;
    unsigned int bitarraysize = 1+(NumBits*NumInts)/8;
    unsigned char *cc = (unsigned char*) ccc;
    unsigned int i,j;
    unsigned int mask = (unsigned int) 0;
    for (i=0 ; i<NumBits ; i++) mask |= ((unsigned int) 1) << i;
    unsigned char mask1[8],mask2[8];
    mask1[0] = 255;
    mask1[1] = 127;
    mask1[2] = 63;
    mask1[3] = 31;
    mask1[4] = 15;
    mask1[5] = 7;
    mask1[6] = 3;
    mask1[7] = 1;
    for (i=0 ; i<8 ; i++) mask2[i] = ~mask1[i];
    unsigned int bitshift=0;
    unsigned int bytes=0;
    unsigned int value;
    unsigned char *cptr = (unsigned char*)&value;
    for (i=0 ; i<NumInts ; i++) {
        value = (unsigned int) 0;
        for (j=0 ; j<NumBits/8 ; j++) {
            cptr[j]  = (cc[bytes+0+j]>>(  bitshift))&mask1[bitshift];
            cptr[j] |= (cc[bytes+1+j]<<(8-bitshift))&mask2[bitshift];
          }
        cptr[j]  = (cc[bytes+0+j]>>(  bitshift))&mask1[bitshift];
        if ((bitshift + NumBits%8) > 7) cptr[j] |= (cc[bytes+1+j]<<(8-bitshift))&mask2[bitshift];
 
        bitshift += NumBits;
        while (bitshift > 7) {bitshift-=8; bytes++;}
        if (swap) {
            ((unsigned char*)(&x[i]))[3] = cptr[0]; 
            ((unsigned char*)(&x[i]))[2] = cptr[1]; 
            ((unsigned char*)(&x[i]))[1] = cptr[2]; 
            ((unsigned char*)(&x[i]))[0] = cptr[3]; 
          }
        else x[i] = value;
        x[i] &= mask;
      }

    return bitarraysize;
  }

/*!
**\fn OB_io_read_binary_compressed(istream& istr, unsigned int *x, unsigned int NumBits, unsigned int NumInts)
**\brief Reads compressed binary data that was written with either
**OB_io_write_binary_compressed(char*ccc, unsigned int *x, unsigned int NumBits, unsigned int NumInts)
**or
**OB_io_write_binary_compressed(ostream& ostr, unsigned int *x, unsigned int NumBits, unsigned int NumInts)
**\param istr An input file stream the binary is being read from
**\param x An array of unsigned ints the binary data will be extracted to
**\param NumBits The number of bits used to store each unsigned integer
**\param NumInts The number of integers stored in istr 
**\return The number of bytes read from istr.  This value should
**always be 1+(NumBits*NumInts)/8 (integer division).
*/
unsigned int OB_io_read_binary_compressed(istream& istr, unsigned int *x, unsigned int NumBits, unsigned int NumInts)
  {
    //If there are no values or x is NULL do nothing
    if (!NumInts || x==NULL) return 0;
 
    //If 32 or more bits are specified simply read uncompressed
    if (NumBits >=32) return OB_io_read_binary(istr,(char*) x, sizeof(unsigned int), NumInts);

    unsigned int bitarraysize = 1+(NumBits*NumInts)/8;
    char *cc = new char [bitarraysize];
    istr.read(cc,bitarraysize);
    unsigned int idx = OB_io_read_binary_compressed(cc,x,NumBits,NumInts);
    delete [] cc;
    return idx;
  }

/*!
**\fn OB_io_read_binary_compressed(char*ccc, float *x, unsigned int NumBits, unsigned int Nfloats)
**\brief Reads compressed binary data that was written with either
**OB_io_write_binary_compressed(char*ccc, float *x, unsigned int NumBits, unsigned int NumInts)
**or
**OB_io_write_binary_compressed(ostream& ostr, float *x, unsigned int NumBits, unsigned int NumInts)
**\param ccc A character array the binary is being read from
**\param x An array of floats the binary data will be extracted to
**\param NumBits The number of bits used to store each float
**\param Nfloats The number of floats to be read 
**\return The number of bytes read from ccc.  This value should
**always be 1+(NumBits*NumInts)/8 (integer division).
*/
unsigned int OB_io_read_binary_compressed(char*ccc, float *x, unsigned int NumBits, unsigned int Nfloats)
  {
    unsigned int idx=0;

    //If there are no values or x is NULL do nothing
    if (!Nfloats || x==NULL) return 0;
 
    //If 32 or more bits are specified simply write uncompressed
    if (NumBits >=32) return OB_io_read_binary(&ccc[idx],(char*) x, sizeof(unsigned int), Nfloats);

    //Read max and min
    float max;
    float min;
    idx += OB_io_read_binary(&ccc[idx],(char*)&max,sizeof(float),1);
    idx += OB_io_read_binary(&ccc[idx],(char*)&min,sizeof(float),1);
 
    //Read in values
    unsigned int i;
    unsigned int maxint = (unsigned int) 0;
    for (i=0 ; i<NumBits ; i++) maxint |= ((unsigned int)1)<<i;
    unsigned int *ix = new unsigned int [Nfloats];
    idx += OB_io_read_binary_compressed(&ccc[idx],ix,NumBits,Nfloats);
    for (i=0 ; i<Nfloats ; i++) x[i] = (max-min)*((float)ix[i]/(float)maxint)+min;
    delete [] ix;

    return idx;
  }

/*!
**\fn OB_io_read_binary_compressed(istream& istr, float *x, unsigned int NumBits, unsigned int Nfloats)
**\brief Reads compressed binary data that was written with either
**OB_io_write_binary_compressed(char*ccc, float *x, unsigned int NumBits, unsigned int NumInts)
**or
**OB_io_write_binary_compressed(ostream& ostr, float *x, unsigned int NumBits, unsigned int NumInts)
**\param istr An input stream  the binary is being read from
**\param x An array of floats the binary data will be extracted to
**\param NumBits The number of bits used to store each float
**\param Nfloats The number of float values to be read
**\return The number of bytes read from istr.  This value should
**always be 1+(NumBits*NumInts)/8 (integer division).
*/
unsigned int OB_io_read_binary_compressed(istream& istr, float *x, unsigned int NumBits, unsigned int Nfloats)
  {
    unsigned int idx=0;

    //If there are no values or x is NULL do nothing
    if (!Nfloats || x==NULL) return 0;
 
    //If 32 or more bits are specified simply write uncompressed
    if (NumBits >=32) return OB_io_read_binary(istr,(char*) x, sizeof(unsigned int), Nfloats);

    //Read max and min
    float max;
    float min;
    idx += OB_io_read_binary(istr,(char*)&max,sizeof(float),1);
    idx += OB_io_read_binary(istr,(char*)&min,sizeof(float),1);
 
    //Read in values
    unsigned int i;
    unsigned int maxint = (unsigned int) 0;
    for (i=0 ; i<NumBits ; i++) maxint |= ((unsigned int)1)<<i;
    unsigned int *ix = new unsigned int [Nfloats];
    idx += OB_io_read_binary_compressed(istr,ix,NumBits,Nfloats);
    for (i=0 ; i<Nfloats ; i++) x[i] = (max-min)*((float)ix[i]/(float)maxint)+min;
    delete [] ix;

    return idx;
  }


/*!
**\fn OB_io_util_calc_NumBits(unsigned int *x, unsigned int N)
**\param x an unsigned integer array
**\paran N The number of values in the unsigned integer array x.
**\return The number of bits per value required to store the
**unsigned integer values in x with the 
**OB_io_write_binary_compressed(char*ccc, unsigned int *x, unsigned int NumBits, unsigned int NumInts)
**and
**OB_io_write_binary_compressed(ostream& ostr, unsigned int *x, unsigned int NumBits, unsigned int NumInts)
** functions.
*/
unsigned int OB_io_util_calc_NumBits(unsigned int *x, unsigned int N)
  {
    if (!N) return 0;
    unsigned int i;
    unsigned int max=0;
    for (i=0 ; i<N ; i++) if (x[i] > max) max = x[i];
    unsigned int mmax=(unsigned int) 0;
    for (i=0 ; i<32 && max > mmax ; i++) mmax |= ((unsigned int)1)<<i;
    return i;
  }

/*!
**\fn OB_io_util_calc_NumBits(float *x, unsigned int N, float res)
**\param x A float array
**\param N The number of values in the float array x.
**\param The resolution desired.
**\returns The number of bits per value required to store float
**the float values in x at a resolution of res (or better) with the
**OB_io_write_binary_compressed(char*ccc, float *x, unsigned int NumBits, unsigned int Nfloats)
**or
**OB_io_write_binary_compressed(istream& istr, float *x, unsigned int NumBits, unsigned int Nfloats)
**functions.
*/
unsigned int OB_io_util_calc_NumBits(float *x, unsigned int N, float res)
  {
    if (!N) return 0;

    //Find max and min
    unsigned int i;
    float max = x[0];
    float min = x[0];
    for (i=0 ; i<N ; i++) {
        if (x[i] > max) max = x[i];
        if (x[i] < min) min = x[i];
      }

    unsigned int imax = (unsigned int) ((max-min)/res+0.5f);
    unsigned int mmax=(unsigned int) 0;
    for (i=0 ; i<32 && imax > mmax ; i++) mmax |= ((unsigned int)1)<<i;
    return i;
  }


/******************* Writing OBBinaryIO objects ************************/
ostream& operator<<(ostream& ostr, const OBBinaryIO& obj)
  {
    unsigned int size = obj.BinarySize();
    OB_io_write_binary(ostr,(const char*)&size,sizeof(unsigned int),1);
    obj.WriteBinary(ostr);
    return ostr;
  }
istream& operator>>(istream& istr, OBBinaryIO& obj)
  {
    unsigned int size;
    OB_io_read_binary(istr,(char*)&size,sizeof(unsigned int),1);
    unsigned int size2;
    size2 = obj.ReadBinary(istr);
    if (size != size2) {
        cerr << "WARNING operator>>(istream& istr, OBBinaryIO& obj) " << endl;
        cerr << "Record size (" << size << ") is not equal to number of bytes read (" << size2 << ")" << endl;
      }
    return istr;
  }


}// End namespace OpenBabel

